#!/usr/bin/env python3
"""Exact CSS contraction/bound controls; abstract GF(2) arithmetic only."""
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import json

BASE=Path(__file__).resolve().parent


def parity(x): return x.bit_count() & 1


def bitlist(x):
    while x:
        t=x & -x
        yield t.bit_length()-1
        x^=t


def sparse(mask):return list(bitlist(mask))


def basis(rows):
    out={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in out: x^=out[p]
            else:
                out[p]=x
                break
    return out


def reduces(x,piv):
    while x:
        p=x.bit_length()-1
        if p not in piv: break
        x^=piv[p]
    return x


def rank(rows): return len(basis(rows))


def all_vectors(rows):
    out=[0]
    for r in rows:
        out += [q^r for q in out]
    return out


def kernel(rows,n):
    echelon=basis(rows)
    vec=[]
    for j in range(n):
        if j in echelon: continue
        x=1<<j
        for p,row in sorted(echelon.items()):
            if parity(row&x): x|=1<<p
        assert all(parity(row&x)==0 for row in rows)
        vec.append(x)
    assert rank(vec)==n-len(echelon)
    return vec


def sector_distance(check,trivial,n,cap=None):
    tt=basis(trivial)
    ker=kernel(check,n)
    if len(ker)<=22:
        tests=[(x.bit_count(),x) for x in all_vectors(ker) if reduces(x,tt)]
        if not tests: return None,None,"exact_kernel"
        d,w=min(tests)
        assert cap is None or d<=cap
        return d,w,"exact_kernel"
    for d in range(1,(cap or n)+1):
        for inds in combinations(range(n),d):
            w=sum(1<<j for j in inds)
            if all(parity(w&z)==0 for z in check) and reduces(w,tt):
                return d,w,"increasing_exact_physical_weights"
    raise AssertionError("cap did not contain a witness")


class Union:
    def __init__(self,n): self.p=list(range(n))
    def r(self,j):
        while self.p[j]!=j:
            self.p[j]=self.p[self.p[j]]
            j=self.p[j]
        return j
    def merge(self,i,j):
        self.p[self.r(j)]=self.r(i)
    def groups(self):
        out={}
        for i in range(len(self.p)):
            out.setdefault(self.r(i),[]).append(i)
        return list(out.values())


def commute(X,Z):
    assert all(parity(a&b)==0 for a in X for b in Z)


def contract(X,Z,n):
    """Intrinsic weight<=2 X subgroup and discard pinned/product blocks."""
    xx=basis(X)
    U=Union(n)
    pins={i for i in range(n) if reduces(1<<i,xx)==0}
    pairs=0
    for i,j in combinations(range(n),2):
        if reduces((1<<i)|(1<<j),xx)==0:
            U.merge(i,j)
            pairs+=1
    remove=[]
    blocks=[]
    for qq in U.groups():
        if pins.intersection(qq):
            assert all(i in pins for i in qq)
            remove.append((qq,"X pinned"))
        else:
            blocks.append(qq)
    masks=[sum(1<<i for i in qq) for qq in blocks]
    A=[sum(parity(row&m)<<j for j,m in enumerate(masks)) for row in X]
    B=[sum(bool(row&m)<<j for j,m in enumerate(masks)) for row in Z]
    for row in Z:
        assert all((row&m) in (0,m) for m in masks)
        assert all(not ((1<<i)&row) for i in pins)
    commute(A,B)
    bb=basis(B)
    fixed={j for j in range(len(blocks)) if reduces(1<<j,bb)==0}
    for j in fixed:
        assert not any(a&(1<<j) for a in A)
        remove.append((blocks[j],"Z fixed"))
    keep=[j for j in range(len(blocks)) if j not in fixed]
    rowmap=lambda rows:[sum(bool(row&(1<<j))<<i for i,j in enumerate(keep)) for row in rows]
    A,B=rowmap(A),rowmap(B)
    blocks=[blocks[j] for j in keep]
    weights=list(map(len,blocks))
    assert n-rank(X)-rank(Z)==len(blocks)-rank(A)-rank(B)
    aa=basis(A)
    assert all(reduces(1<<i,aa) for i in range(len(blocks)))
    assert all(reduces((1<<i)|(1<<j),aa) for i,j in combinations(range(len(blocks)),2))
    return A,B,blocks,weights,remove,pairs


def replay(name,n,X,Z,D,expected,expected_ell=None):
    commute(X,Z)
    rX,rZ=rank(X),rank(Z)
    k=n-rX-rZ
    dX,wx,xmethod=sector_distance(Z,X,n,expected[1])
    dZ,wz,zmethod=sector_distance(X,Z,n,expected[2])
    assert (k,dX,dZ)==expected
    cols=lambda checks:max([0]+[sum((v>>i)&1 for v in checks) for i in range(n)])
    rx_max=max(x.bit_count() for x in X)
    rz_max=max([0]+[z.bit_count() for z in Z])
    A,B,blocks,w,removed,npairs=contract(X,Z,n)
    m=len(w)
    # Independent compressed proof of both distance values (raw physical
    # minimum comes from the search above; no floating point).
    def compressed_d(C,T,weighted):
        tt=basis(T)
        return min(sum(w[i] if weighted else 1 for i in bitlist(v))
                  for v in all_vectors(kernel(C,m)) if reduces(v,tt))
    assert compressed_d(B,A,False)==dX
    assert compressed_d(A,B,True)==dZ
    S=sum((1<<i) for i,z in enumerate(w) if z<D)
    L=((1<<m)-1)^S
    rAS=rank(a&S for a in A)
    ell=S.bit_count()-rAS-(rank(B)-rank(b&L for b in B))
    n_small=sum(w[i] for i in bitlist(S))
    n_large=sum(w[i] for i in bitlist(L))
    assert ell>=0 and ell<=k and ell<=S.bit_count()-rAS
    assert S.bit_count()<=rx_max*rAS
    assert dZ>=D
    if n_small<D: assert ell==0
    # Enumerate only the small contracted sector to independently verify
    # rank-counted logical Z support (not necessarily the full distance).
    az=basis(B)
    z_supports=[v for v in all_vectors(kernel(A,m)) if not(v&L)]
    assert rank(list(B)+z_supports)-rank(B)==ell
    if expected_ell is not None: assert ell==expected_ell
    graph_ok=all((b&L).bit_count()<=2 for b in B)
    mode={"n_after_discard":sum(w),"small_blocks":S.bit_count(),"n_small":n_small,"n_large":n_large,
          "rank_A_on_small":rAS,"ell_Z_on_small":ell,"applies_row_weight":rz_max<3*D,
          "applies_actual_large_check_projection":graph_ok,"bound_logical_X_graph_families":[],
          "fixed_deleted":removed}
    if graph_ok:
        LL=sparse(L)
        ix={j:i for i,j in enumerate(LL)}
        U=Union(len(LL))
        forced=[]
        for b in B:
            inds=sparse(b&L)
            if len(inds)==2: U.merge(ix[inds[0]],ix[inds[1]])
            if len(inds)==1: forced.append(ix[inds[0]])
        AA=basis(A)
        logical=[]
        component_cols=[(cc,sum(1<<LL[i] for i in cc)) for cc in U.groups() if not(set(cc)&set(forced))]
        for _,c in component_cols:
            if reduces(c,AA):logical.append(c)
        assert rank(list(A)+logical)-rank(A)==k-ell
        assert all(v.bit_count()>=dX for v in logical)
        assert n_large>=D*dX*(k-ell)
        # Check independent exact localized X minimum on compressed C_B.
        local=[u.bit_count() for u in all_vectors(kernel(B,m)) if not(u&S) and reduces(u,AA)]
        exact_local=min(local) if local else None
        assert exact_local==(min(v.bit_count() for v in logical) if logical else None)
        mode.update(exact_dX_support_large=exact_local, bound_logical_X_graph_families=[sparse(v) for v in logical],
                    n_bound_after_removals=D*dX*(k-ell)+n_small,
                    n_bound=D*dX*(k-ell)+n_small+sum(len(positions) for positions, _ in removed))
    else:
        assert name=="three_block_parity_boundary" and n_large<D*dX*(k-ell)
    return {"control":name,"n":n,"k":k,"d_X":dX,"d_Z":dZ,"target_D":D,"ranks":[rX,rZ],
            "X_logical_support":sparse(wx),"Z_logical_support":sparse(wz),
            "physical_distance_verification":[xmethod,zmethod],
            "physical_generator_max_row":[rx_max,rz_max],"physical_generator_max_col":[cols(X),cols(Z)],
            "presented_generators":{"X":[sparse(q) for q in X],"Z":[sparse(q) for q in Z]},
            "intrinsic_XX_pairs":npairs,"blocks_after_removals":blocks,"certificate":mode}


def lift(A,B,ws):
    off=0
    blocks=[]
    LL=[]
    for w in ws:
        ids=list(range(off,off+w))
        off+=w
        blocks.append(ids)
        LL.extend((1<<i)|(1<<j) for i,j in zip(ids,ids[1:]))
    masks=[sum(1<<i for i in g) for g in blocks]
    X=LL+[sum(1<<blocks[i][0] for i in bitlist(a)) for a in A]
    Z=[sum(masks[i] for i in bitlist(b)) for b in B]
    return off,X,Z,blocks


def main():
    n1,X1,Z1,_=lift([], [3,6], [6]*3)
    G=[sum(1<<(j-1) for j in range(1,8) if (j>>i)&1) for i in range(3)]
    n2,X2,Z2,bl2=lift(G,G,[2]*7)
    Xs=X1+[x<<n1 for x in X2]; Zs=Z1+[z<<n1 for z in Z2]
    d=n1+n2
    _,l1,_=sector_distance(Z1,X1,n1,3)
    _,l2,_=sector_distance(Z2,X2,n2,3)
    joint=l1|(l2<<n1)
    bX,bZ=Xs+[joint],Zs
    n3,X3,Z3,_=lift([], [7], [3]*3)
    n4,X4,Z4,_=lift([], [], [6])
    # Equivalent row span presenting a light X row as a product with a
    # heavier generator; intrinsic pair contraction must not change.
    Xj=list(Xs)
    Xj[0]^=Xj[-1]
    pinX=X1+[1<<n1,(3<<(n1+1))]
    pinZ=Z1+[(3<<(n1+1))]
    inputs=[
        ("pair_only_asymmetric",n1,X1,Z1,6,(1,3,6),0),
        ("repeated_steane_bias_control",n2,X2,Z2,6,(1,3,6),1),
        ("mixed_check_budget_exact_saturation",d,Xs,Zs,6,(2,3,6),1),
        ("mixed_X_links_clear_small_logicals",d,bX,bZ,6,(1,3,12),0),
        ("basis_change_intrinsic_pair_invariance",d,Xj,Zs,6,(2,3,6),1),
        ("three_block_parity_boundary",n3,X3,Z3,3,(2,2,3),0),
        ("one_sector_protected",n4,X4,Z4,6,(1,1,6),0),
        ("pinned_and_fixed_factors",n1+3,pinX,pinZ,6,(1,3,6),0),
    ]
    arr=[replay(*z) for z in inputs]
    assert arr[2]["blocks_after_removals"]==arr[4]["blocks_after_removals"]
    assert arr[2]["certificate"]["n_bound"]==32
    assert arr[3]["certificate"]["n_bound"]==32
    assert sum(1 for a in arr if a["certificate"]["applies_actual_large_check_projection"])==7
    out={"algorithm":"binary exact; physical sector searches check the compressed rank/weight derivation",
         "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),"controls":arr}
    (BASE/"checks.json").write_text(json.dumps(out,indent=2)+"\n")
    for a in arr:
        c=a["certificate"]
        print(a["control"],a["n"],a["k"],a["d_X"],a["d_Z"],a["ranks"],a["physical_generator_max_row"],a["physical_generator_max_col"],c["ell_Z_on_small"],c.get("n_bound"))


if __name__=="__main__":
    main()
