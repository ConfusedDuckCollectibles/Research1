
"use strict";
const p=13,mod=x=>(x%p+p)%p, asc=(a,b)=>a-b;
const B=[0,1,2,3,10,11,12,13,14,16,17,20,28,30,31,56,58,59,73,79,86,88,89,93,94,96,103,109,123,124,126,151,152,154,156,162,165,166,168];
const P=B.map(i=>[Math.floor(i/p),i%p]);
const inv=Array.from({length:p},(_,i)=>i===0?0:Array.from({length:p},(_,j)=>j).find(j=>i*j%p===1));
const pow=(x,k)=>{let y=1;while(k--)y=mod(x*y);return y};
const check=(s,label)=>{if(!s)throw Error(label)};
function fibres(q){
  const f=Array.from({length:p},()=>[]);
  P.forEach(([x,y],i)=>f[q===p?x:mod(y-q*x)].push([q===p?y:x,i]));
  return f;
}
function heights(T){
  // Use first and LAST abscissae to normalize, one canonical representative
  // modulo adding all affine linear functions to the heights.
  const n=T.length, H=Array(n).fill(0), I=Array.from({length:n},(_,i)=>T.map(t=>inv[mod(t-T[i])]));
  const counts=Array(p+1).fill(0), all=Array(p+1).fill(0);
  let valid=0, total=0,witness=null;
  function test(){
    total++;
    let slope=0,cap=true;
    for(let i=0;i<n;i++){
      const used=Array(p).fill(0);
      for(let j=0;j<n;j++)if(i!==j){
        const s=mod((H[j]-H[i])*I[i][j]);
        if(++used[s]>2)cap=false;
        slope|=1<<s;
      }
    }
    let k=0;while(slope){slope &= slope-1;k++}
    all[k]++;
    if(cap){valid++;counts[k]++;if(k===7&&witness===null)witness=H.slice();}
  }
  function go(i){if(i===n-1){test();return;}for(let h=0;h<p;h++){H[i]=h;go(i+1)}}
  go(1);
  return {T,total,valid,directionCountsByValue:counts,minimum:counts.findIndex(Boolean),witness,
    failingSmallDirWithCap:counts.slice(0,7).reduce((a,b)=>a+b,0),
    anySmallDirWithoutCap:all.slice(0,7).reduce((a,b)=>a+b,0)};
}
const sizes=Array.from({length:p+1},(_,q)=>Math.max(...fibres(q).map(b=>b.length)));
const normalSlices=Array.from({length:p+1},(_,q)=>fibres(q)[0].length);
const active=sizes.flatMap((a,q)=>a>=7?[q]:[]);
const remaining=sizes.flatMap((a,q)=>a<7?[q]:[]);
function connected(q,r){
  const K=[...fibres(q),...fibres(r)].flatMap(xs=>xs.length?[[xs[0][1],xs.map(z=>z[1])]]:[]);
  const parent=P.map((_,i)=>i);
  const root=i=>parent[i]===i?i:root(parent[i]);
  for(const [m,all] of K)for(const k of all)parent[root(k)]=root(m);
  return new Set(parent.map((_,i)=>root(i))).size;
}
const C=remaining.flatMap((q,i)=>remaining.slice(i+1).map(r=>({q,r,components:connected(q,r)})));
check(C.length===28&&C.every(x=>x.components===1),"all two-direction unions not connected");
const T=fibres(8)[1].map(z=>z[0]).sort(asc), HGate=heights(T);
check(T.join()==="0,2,4,6,7,11,12"&&HGate.total===371293&&HGate.minimum===7,"q8 height gate");
check(JSON.stringify(active)==="[1,2,4,8,10,13]"&&sizes[2]===8,"wrong heavy slopes");
// Recompute degree-7 mixed-moment gate via row reduction on evaluations of monomials.
function mono(d,parity) {
  const A=[];
  for(let i=0;i<=d;i++)for(let j=0;i+j<=d;j++)if((i+j)%2===parity)A.push([i,j]);
  return A;
}
function rref(A,n){
  const C=A.map(row=>row.map(mod)), piv=[];
  for(let i=0;i<n;i++){
    let r=piv.length;
    let j=C.findIndex((v,k)=>k>=r&&v[i]);
    if(j<0)continue;
    [C[r],C[j]]=[C[j],C[r]];
    C[r]=C[r].map(x=>mod(x*inv[C[r][i]]));
    for(let t=0;t<C.length;t++)if(t!==r&&C[t][i]){
      const c=C[t][i];
      for(let w=0;w<n;w++)C[t][w]=mod(C[t][w]-c*C[r][w]);
    }
    piv.push(i);
  }
  const basis=[];
  for(let j=0;j<n;j++)if(!piv.includes(j)){
    const z=Array(n).fill(0);z[j]=1;
    for(let i=0;i<piv.length;i++)z[piv[i]]=mod(-C[i][j]);
    basis.push(z);
  }
  return {rank:piv.length,basis};
}
const opp=([x,y])=>p*mod(-x)+mod(-y);
check(P.every(q=>B.includes(opp(q))),"not a central base");
const reps=P.filter(([x,y])=>p*x+y<opp([x,y]));
check(reps.length===19,"wrong central decomposition");
const monVec=([a,b])=>reps.map(([x,y])=>mod(pow(x,a)*pow(y,b)));
const U=([x,y])=>reps.map(([a,b])=>+(a===x&&b===y));
const even=rref(mono(6,0).map(monVec),19).basis;
const odd=rref(mono(6,1).map(monVec).concat([U([1,0]),U([0,1])]),19).basis;
check(even.length===3&&odd.length===5,"kernel dims");
let ecount=0,fullOddRank=0,evenObstructed=0;
function evenGate(co){
  const e=reps.map((_,i)=>mod(co.reduce((v,c,j)=>v+c*(even[j][i]),0)));
  const oddCross=mono(5,1).map(monVec).map(m=>odd.map(o=>mod(o.reduce((v,a,i)=>v+a*m[i]*e[i],0))));
  if(rref(oddCross,5).rank===5)fullOddRank++;
  if(mono(5,0).map(monVec).some(m=>mod(e.reduce((v,a,i)=>v+a*a*m[i],0))))evenObstructed++;
  ecount++;
}
for(let a=0;a<p;a++)for(let b=0;b<p;b++)evenGate([1,a,b]);
for(let a=0;a<p;a++)evenGate([0,1,a]);
evenGate([0,0,1]);
check(ecount===183&&fullOddRank===183&&evenObstructed===183,"degree seven counterexample");
const baseMomentFailures=[];
for(let a=0;a<=7;a++)for(let b=0;a+b<=7;b++)if(mod(P.reduce((z,[x,y])=>z+pow(x,a)*pow(y,b),0)))baseMomentFailures.push([a,b]);
check(baseMomentFailures.length===0,"base not moment-flat");
const out={heavyDirections13MeansVertical:active,maximumSliceSizes:sizes,levelZeroSizes:normalSlices,
  remainingDirections:remaining,connectivity28:C,heightGate:HGate,
  moments:{evenDim:even.length,oddDimAfterShear:odd.length,projectiveEven:ecount,fullOddRank,evenObstructed,baseMomentFailures}};
console.log(JSON.stringify(out,null,2));