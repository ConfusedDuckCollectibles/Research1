import Checker

namespace SubsetCover

/-- Compact literal data from the finite mathematical reduction.
    Rank offsets are measured from the reference low cut. -/
structure Template where
  h : Nat
  width : Nat
  trim : Nat
  tailSize : Nat
  bans : Array Nat

def Template.lowSize (t : Template) : Nat := t.width + 64 - t.trim
def Template.sites (t : Template) : Nat := t.lowSize + t.width + 1

def Template.rank (t : Template) (s : Nat) : Nat :=
  if s < t.lowSize then s + t.trim else t.h + (s - t.lowSize)

def Template.site (t : Template) (r : Nat) : Option Nat :=
  if t.trim ≤ r && r ≤ t.width + 63 then some (r - t.trim)
  else if t.h ≤ r && r ≤ t.h + t.width then some (t.lowSize + r - t.h)
  else none

def sumMask (h mask : Nat) : List Nat → Nat
  | [] => 0
  | d::ds => (if mask % 2 == 1 then h-d else 0) + sumMask h (mask/2) ds

def Template.upper (t : Template) (deficits : List Nat) : List (Nat × Nat) :=
  (List.range (2 ^ t.tailSize)).map (fun m => (m, sumMask t.h m deficits))

/-- tag is offset by tailSize*h so all tags and translates are natural numbers. -/
def Template.incidence (t : Template) (tag : Nat) :
    List (Nat × Nat) → Option Row
  | [] => some []
  | (m,e)::es =>
    let v := tag + e
    let a := t.tailSize * t.h
    if v < a then t.incidence tag es
    else
      let r := v - a
      match t.site r with
      | none => t.incidence tag es
      | some s =>
        if ((t.bans[r]?).getD 0).testBit m then none
        else (t.incidence tag es).map (s :: ·)

def Template.board (t : Template) (deficits : List Nat) : Board :=
  let es := t.upper deficits
  ⟨t.sites, fun s =>
    es.filterMap (fun (_,e) => t.incidence (t.rank s + t.tailSize*t.h - e) es)⟩

structure Entry where
  deficits : List Nat
  certificate : List Nat

def passed (t : Template) (e : Entry) : Bool := check (t.board e.deficits) e.certificate

theorem finite_batch_sound {t : Template} {entries : List Entry}
    (accepted : entries.all (passed t) = true)
    {e : Entry} (member : e ∈ entries) : ¬ Residual (t.board e.deficits) [] := by
  exact check_sound ((List.all_eq_true.mp accepted) e member)

/-- Site buckets are supplied externally. The full reconstructed board
    remains the fallback away from cached sites. -/
structure CachedEntry where
  deficits : List Nat
  certificate : List Nat
  buckets : List (Nat × List Row)

def CachedEntry.board (e : CachedEntry) (t : Template) : Board :=
  let full := t.board e.deficits
  ⟨full.sites, fun s =>
    match e.buckets.find? (fun p => p.1 == s) with
    | some p => p.2
    | none => full.rowsAt s⟩

def CachedEntry.cached (e : CachedEntry) : Bool :=
  e.certificate.all fun s => e.buckets.any fun p => p.1 == s

def CachedEntry.agrees (e : CachedEntry) (t : Template) : Bool :=
  e.buckets.all fun (s, rows) => (t.board e.deficits).rowsAt s == rows

def acceptedCached (t : Template) (e : CachedEntry) : Bool :=
  e.cached && check (e.board t) e.certificate

theorem cached_batch_sound {t : Template} {entries : List CachedEntry}
    (accepted : entries.all (acceptedCached t) = true)
    {e : CachedEntry} (member : e ∈ entries) : ¬ Residual (e.board t) [] := by
  have test := (List.all_eq_true.mp accepted) e member
  exact check_sound (Bool.and_eq_true_iff.mp test).2

/-- This transfer is kernel proved, but its agreement hypotheses are
    checked by a compiled cross-check in the production run, not proved
    by kernel reduction for the entire corpus. -/
theorem CachedEntry.board_eq_of_agrees {e : CachedEntry} {t : Template}
    (ok : e.agrees t = true) : e.board t = t.board e.deficits := by
  cases eq : t.board e.deficits with
  | mk sites rows =>
    simp only [CachedEntry.board, eq]
    congr
    funext s
    cases found : e.buckets.find? (fun p => p.1 == s) with
    | none => rfl
    | some p =>
      have test := List.find?_some found
      have key : p.1 = s := beq_iff_eq.mp test
      have candidate := (List.all_eq_true.mp ok) p (List.mem_of_find?_eq_some found)
      have good : (t.board e.deficits).rowsAt p.1 = p.2 := beq_iff_eq.mp candidate
      rw [eq, key] at good
      exact good.symm

theorem cached_batch_transfer {t : Template} {entries : List CachedEntry}
    (accepted : entries.all (acceptedCached t) = true)
    {e : CachedEntry} (member : e ∈ entries) (agrees : e.agrees t = true) :
    ¬ Residual (t.board e.deficits) [] := by
  rw [← CachedEntry.board_eq_of_agrees agrees]
  exact cached_batch_sound accepted member

#print axioms finite_batch_sound
#print axioms cached_batch_sound
#print axioms cached_batch_transfer

end SubsetCover
