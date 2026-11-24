import Lean

/-!
A small exact-cover refutation checker. Its input board is a site-indexed
incidence relation; the data layer is responsible for the correspondence
between this relation and any outside mathematical problem.
-/
namespace SubsetCover

abbrev Row := List Nat
structure Board where
  sites : Nat
  rowsAt : Nat → List Row

def incompatibleFree (r occupied : Row) : Bool :=
  r.all (fun t => !(occupied.contains t))

def branches (board : Board) (s : Nat) (occupied : Row) : List Row :=
  (board.rowsAt s).filter (fun r => r.contains s && incompatibleFree r occupied)

/-- There is a remaining exact cover using rows registered at their covered
sites. The remaining rows are pairwise disjoint and avoid earlier choices. -/
def Residual (board : Board) (occupied : Row) : Prop :=
  ∃ selected : Row → Prop,
    (∀ r, selected r → ∀ t ∈ r, t ∉ occupied) ∧
    (∀ r r' t, selected r → selected r' → t ∈ r → t ∈ r' → r = r') ∧
    (∀ s < board.sites, s ∉ occupied →
      ∃ r, selected r ∧ r ∈ board.rowsAt s ∧ s ∈ r)

def everyChild (run : Row → List Nat → Option (List Nat)) :
    List Row → List Nat → Option (List Nat)
  | [], cert => some cert
  | r::rs, cert =>
    match run r cert with
    | none => none
    | some pending => everyChild run rs pending

/-- Preorder site certificate. Children are implicit: there must be one and
only one checked subtree for each compatible row through the nominated site. -/
def consume (board : Board) : Nat → Row → List Nat → Option (List Nat)
  | 0, _, _ => none
  | _+1, _, [] => none
  | fuel+1, occupied, s::cert =>
    if s < board.sites ∧ s ∉ occupied then
      everyChild (fun r => consume board fuel (r ++ occupied))
        (branches board s occupied) cert
    else none

def check (board : Board) (cert : List Nat) : Bool :=
  match consume board cert.length [] cert with
  | some [] => true
  | _ => false

theorem everyChild_calls
    {run : Row → List Nat → Option (List Nat)} {rows : List Row}
    {input output : List Nat} (passed : everyChild run rows input = some output)
    {r : Row} (member : r ∈ rows) :
    ∃ a b, run r a = some b := by
  induction rows generalizing input with
  | nil => simp at member
  | cons head tail ind =>
    cases step : run head input with
    | none => simp [everyChild, step] at passed
    | some pending =>
      simp [everyChild, step] at passed
      rcases List.mem_cons.mp member with same | later
      · subst r
        exact ⟨input, pending, step⟩
      · exact ind passed later

theorem remove_row {board : Board} {occupied r : Row}
    {selected : Row → Prop}
    (avoid : ∀ q, selected q → ∀ t ∈ q, t ∉ occupied)
    (pair : ∀ q q' t, selected q → selected q' → t ∈ q → t ∈ q' → q = q')
    (cover : ∀ s < board.sites, s ∉ occupied →
      ∃ q, selected q ∧ q ∈ board.rowsAt s ∧ s ∈ q)
    (picked : selected r) :
    Residual board (r ++ occupied) := by
  refine ⟨fun q => selected q ∧ q ≠ r, ?_, ?_, ?_⟩
  · intro q hq t ht hbad
    rcases List.mem_append.mp hbad with used | old
    · exact hq.2 (pair q r t hq.1 picked ht used)
    · exact avoid q hq.1 t ht old
  · intro q q' t hq hq' hmem hmem'
    exact pair q q' t hq.1 hq'.1 hmem hmem'
  · intro s hs hsfree
    obtain ⟨q, hq, hrow, hsite⟩ :=
      cover s hs (fun old => hsfree (List.mem_append.mpr (Or.inr old)))
    refine ⟨q, ⟨hq, ?_⟩, hrow, hsite⟩
    intro same
    subst q
    exact hsfree (List.mem_append.mpr (Or.inl hsite))

theorem consume_sound {board : Board} {fuel : Nat} {occupied input output : Row}
    (passed : consume board fuel occupied input = some output) :
    ¬ Residual board occupied := by
  induction fuel generalizing occupied input output with
  | zero => simp [consume] at passed
  | succ fuel ind =>
    cases input with
    | nil => simp [consume] at passed
    | cons s cert =>
      by_cases site : s < board.sites ∧ s ∉ occupied
      · simp only [consume, if_pos site] at passed
        intro existsCover
        obtain ⟨selected, avoid, pair, cover⟩ := existsCover
        obtain ⟨r, picked, rowAt, hasSite⟩ := cover s site.1 site.2
        have ok : r ∈ branches board s occupied := by
          have fits : incompatibleFree r occupied = true := by
            simp only [incompatibleFree, List.all_eq_true]
            intro t ht
            cases result : occupied.contains t with
            | false => rfl
            | true => exact False.elim (avoid r picked t ht (List.contains_iff_mem.mp result))
          simp [branches, rowAt, hasSite, fits]
        obtain ⟨a,b,child⟩ := everyChild_calls passed ok
        exact ind child (remove_row avoid pair cover picked)
      · simp [consume, site] at passed

/-- Successful Boolean checking implies no exact cover of the board. -/
theorem check_sound {board : Board} {cert : List Nat}
    (passed : check board cert = true) : ¬ Residual board [] := by
  unfold check at passed
  cases step : consume board cert.length [] cert with
  | none => simp [step] at passed
  | some rest => exact consume_sound step

#print axioms check_sound

end SubsetCover
