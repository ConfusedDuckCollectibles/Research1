#!/usr/bin/env python3
"""Independent, non-kernel check of the finite data's mathematical encoding."""
import itertools
import json
import os
from math import comb
from pathlib import Path

def cube(n):
    return [
        sum(1 << i for i in subset)
        for d in range(n + 1)
        for subset in itertools.combinations(range(n), d)
    ]

def template(n, trim, k):
    h=sum(comb(i,i//2) for i in range(n))
    lex=cube(n)
    ranks={b:r for r,b in enumerate(lex,1)}
    last=1 << (n-1)
    cuts=[ranks[b] for b in lex if not b&last and ranks[b|last]-ranks[b]==h]
    assert len(cuts)==2**(n//2)
    L,U=cuts[0],cuts[-1]
    w=U-L
    bans=[0]*(h+w+1)
    for shift in (0,h):
        for p0,q0 in itertools.pairwise(cuts):
            p,q=p0+shift,q0+shift
            for j in range(k):
                banned_tail=((1 << k)-1)^((1 << j)-1)
                max_prefix=n-k+j
                direct=sum(b < 1 << max_prefix for b in lex[p:q])
                a,b=len(lex)+1-q,len(lex)+1-p
                mirror=sum(b < 1 << max_prefix for b in lex[a:b])
                if direct in (0,q-p):
                    excluded=sum(1<<i for i in range(1<<k)
                                 if (i&banned_tail==0) != (direct==q-p))
                    for r in range(p+1,q+1):
                        bans[r-L]|=excluded
                if mirror in (0,q-p):
                    excluded=sum(1<<i for i in range(1<<k)
                                 if (i&banned_tail==banned_tail) != (mirror==q-p))
                    for r in range(p,q):
                        bans[r-L]|=excluded
    windows=list(range(trim,w+64))+list(range(h,h+w+1))
    return h,w,bans,windows

def bounds(n, h):
    for ds in itertools.combinations(range(2, 77 if n==11 or n==10 else 21), 3 if n!=9 else 2):
        if ds[0]>5 or ds[1]>20:
            continue
        d=tuple(reversed(ds))+(1,0)
        upper=[0]
        for delta in d:
            upper += [b+h-delta for b in upper]
        yield d,upper

def exact_cache(upper, sites, bans, required):
    indices={r:i for i,r in enumerate(sites)}
    result=[]
    for key in sorted(set(required)):
        list_for_site=[]
        pivot=sites[key]
        for start in upper:
            z=pivot-start
            row=[]
            for m,e in enumerate(upper):
                r=z+e
                if r not in indices:
                    continue
                if bans[r] & (1 << m):
                    break
                row.append(indices[r])
            else:
                list_for_site.append(row)
        result.append([key,list_for_site])
    return result

def accepted(cache, cert, ns):
    board={s:[sum(1<<i for i in row) for row in rows] for s,rows in cache}
    def verify(visited, pos):
        if pos==len(cert):
            return None
        s=cert[pos]
        if not 0<=s<ns or visited>>s&1 or s not in board:
            return None
        pos+=1
        for row in board[s]:
            if visited&row or not row>>s&1:
                continue
            pos=verify(visited|row,pos)
            if pos is None:
                return None
        return pos
    return verify(0,0)==len(cert)

def main():
    corpus=json.loads(Path(os.environ.get('SUBSET_CERTIFICATES', Path(__file__).with_name('certificates.json'))).read_text())
    expected={9:(62,4,156,10),10:(3755,455,6329,12),11:(3755,455,7406,26)}
    assert [c['dimension'] for c in corpus]==list(expected)
    for batch in corpus:
        n=batch['dimension']
        trim=3 if n==10 else 0
        k=4 if n==9 else 5
        h,w,bans,sites=template(n,trim,k)
        assert (batch['h'],batch['width'],batch['trim'],batch['tail'],batch['bans'])==(h,w,trim,k,bans)
        actual={tuple(d):(cert,rows) for d,cert,rows in batch['entries']}
        assert len(actual)==len(batch['entries'])
        repeated,dirs=[],set()
        index_count,max_len=0,0
        for d,upper in bounds(n,h):
            if len(upper)!=len(set(upper)):
                repeated.append(d)
                continue
            dirs.add(d)
            cert,cache=actual[d]
            assert cache==exact_cache(upper,sites,bans,cert)
            assert accepted(cache,cert,len(sites))
            assert not any(accepted(cache,c,len(sites))
                           for c in [cert[:-1],cert+[0],[len(sites)]+cert[1:]])
            index_count+=len(cert)
            max_len=max(max_len,len(cert))
        assert dirs==set(actual)
        assert (len(actual),len(repeated),index_count,max_len)==expected[n]
        print(f'Q{n}: {len(actual)} tails, {len(repeated)} repeated floors, '
              f'{index_count} nodes, max {max_len}; exact window, mask bans, '
              f'incidences and {3*len(actual)} negative cases verified.')

if __name__=='__main__':
    main()
