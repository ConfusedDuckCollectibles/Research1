"""Exact planar and small-profile controls; run with Python 3.10+."""
from itertools import combinations, product
from pathlib import Path

P = 13
INF = P
points = list(product(range(P), repeat=2))
N = len(points)
index = {pt:i for i,pt in enumerate(points)}
base = tuple(index[p] for p in [(0,0), (1,0), (0,1)])
dirs = [[None for _ in points] for _ in points]
vecs = [[None for _ in points] for _ in points]
line = [[0]*(P+1) for _ in points]
for i,(a,b) in enumerate(points):
    for j,(c,d) in enumerate(points):
        if i==j:
            continue
        x,y=(c-a)%P,(d-b)%P
        dirs[i][j]=INF if x==0 else y*pow(x,-1,P)%P
        vecs[i][j]=min((x,y),(-x%P,-y%P))
        line[i][dirs[i][j]]|=1<<j
base_dirs = {dirs[i][j] for i,j in combinations(base,2)}


def find_cap(goal, allowed=range(P+1), monochord=False):
    """Return a cap if one exists up to affine normalization of its triangle."""
    adj = [sum(line[i][s] for s in allowed) for i in range(N)]
    avail = (1<<N)-1
    for i in base:
        avail &= adj[i]
    for i,j in combinations(base,2):
        avail &= ~line[i][dirs[i][j]]
    initial = {dirs[i][j]:vecs[i][j] for i,j in combinations(base,2)}
    def rec(chosen,left,used):
        if len(chosen)==goal:
            return chosen
        while len(chosen)+left.bit_count()>=goal:
            bit=left & -left
            l=bit.bit_length()-1
            left^=bit
            insert={dirs[i][l]:vecs[i][l] for i in chosen}
            assert len(insert)==len(chosen)
            if monochord and any(k in used and used[k]!=v for k,v in insert.items()):
                continue
            nxt=left & adj[l]
            for i in chosen:
                nxt &= ~line[i][dirs[i][l]]
            found=rec((*chosen,l),nxt,used|insert if monochord else used)
            if found:
                return found
        return None
    ans=rec(base,avail,initial)
    return None if ans is None else [points[i] for i in ans]


def masks_at_goal(k,goal):
    return [tuple(sorted(base_dirs|set(extra)))
            for extra in combinations(sorted(set(range(P+1))-base_dirs), k-3)
            if find_cap(goal,base_dirs|set(extra)) is not None]


def projective_canon(slopes):
    best=1<<(P+1)
    for a,b,c,d in product(range(P), repeat=4):
        if (a*d-b*c)%P==0:
            continue
        z=0
        for x in slopes:
            u,v=(a,c) if x==INF else ((a*x+b)%P,(c*x+d)%P)
            y=INF if v==0 else u*pow(v,-1,P)%P
            z|=1<<y
        best=min(best,z)
    return best


def affine_masks(k):
    out=set()
    for bad in combinations(range(P),k):
        out.add(min(sum(1<<((a*x+b)%P) for x in bad)
                    for a in range(1,P) for b in range(P)))
    return out


def low_profile_mass(d):
    minimum=P*P
    for head in range(1,P):
        for tail in product(range(P),repeat=d):
            sm=0
            for x in range(P):
                y=head
                for a in tail:
                    y=(x*y+a)%P
                sm+=y
            minimum=min(minimum,sm)
    return minimum


def check_profiles(path):
    blocks = Path(path).read_text().splitlines()
    cases = {(6,8):(117,3,3), (7,7):(585,12,0), (8,6):(3497,290,84),
             (8,5):(1612,20,0), (9,4):(3562,72,0)}
    seen = set()
    def check_hist(case,hist):
        ex={(6,8):{(1,6):3},(8,6):{(0,6):84}}
        assert hist==ex.get(case,{})
    for line_out in blocks:
        words = line_out.split()
        if words[0]=='L':
            if seen:
                check_hist((k,m),hist)
            k,m,db=map(int,words[1::2])
            assert (k,m) not in seen
            assert db==cases[k,m][0]
            expected_masks=affine_masks(k)
            seen.add((k,m)); got={}; hist={}
        elif words[0]=='mask':
            mask,tuples,survivors=map(int,words[1::2])
            assert mask not in got
            got[mask]=(tuples,survivors)
        elif words[0]=='TOTAL':
            number,tuples,survivors=map(int,words[2::2])
            assert number==len(expected_masks) and set(got)==expected_masks
            assert tuples==sum(x for x,y in got.values())==cases[k,m][1]
            assert survivors==sum(y for x,y in got.values())==cases[k,m][2]
            if (k,m)==(8,6):
                canon={projective_canon(tuple(x for x in range(P+1) if not z>>x&1))
                       for z,(x,y) in got.items() if y}
                assert canon=={95,111}
                print("surviving six-good masks: PGL representatives 95, 111")
        else:
            assert words[0]=='doubles'
            d,peak,n = map(int,words[1::2])
            hist[d,peak]=n
    check_hist((k,m),hist)
    assert seen==set(cases)


def main():
    assert find_cap(8, monochord=True) is None
    assert find_cap(7, monochord=True) is not None
    assert not masks_at_goal(8,8)
    assert not masks_at_goal(5,5)
    cap6=masks_at_goal(6,6)
    assert len(cap6)==10
    assert {projective_canon(s) for s in cap6}=={189}
    assert [low_profile_mass(d) for d in (1,2,3)]==[78,52,39]
    print("planar checks: chord bound 7, no (8,8) or (5,5), six-cap direction orbit 189")
    print("least masses of nonconstant degree 1, 2, 3 profiles: 78, 52, 39")
    import sys
    if len(sys.argv)>1:
        check_profiles(sys.argv[1])
        print("full run matches all independent mask and planar controls")

if __name__=='__main__':
    main()
