"use strict";
// Independent F13 mixed-layer replay: variable order from exact cover, full
// stabilizer, and upper-minus-lower projective differences in the z=1 chart.
const p=13,F=Array.from({length:p},(_,x)=>x),mod=x=>(x%p+p)%p, all=8191;
const inv=F.map(x=>F.find(y=>mod(x*y)===1)||0);
const pop=Array.from({length:1<<p},(_,x)=>{let c=0;while(x){x&=x-1;c++;}return c;});
const check=(x,msg)=>{if(!x)throw Error(msg);};
const bit=x=>1<<x;
const rotate=(mask,q)=>{q=mod(q);return ((mask<<q)|(mask>>(13-q)))&all;};
const stringify=f=>f.map(x=>x.toString(16)).join("");
const parse=k=>Array.from(k,c=>parseInt(c,16));
const canon=f=>{
  let best=null;
  for(const a of F){
    let k="";
    for(const x of F)k+=mod(f[mod(x+a)]-f[a]).toString(16);
    if(best===null||k<best)best=k;
  }
  return best;
};
const f=Array(13).fill(0), rem=all^1, whole=new Set, hist={};
const events=Array(14).fill(0), slopes=F.slice(2), lineVal=F.map(x=>F.map(y=>slopes.map(s=>13*s+mod(y-s*x))));
const counters=new Uint8Array(169);
lineVal[0][0].forEach(i=>counters[i]++);
function graphs(U,Y,D,n){
  events[n]++;
  if(U===0){
    const dirs=new Set;
    for(let x=1;x<13;x++)for(let y=0;y<x;y++)dirs.add(mod((f[x]-f[y])*inv[x-y]));
    hist[dirs.size]=(hist[dirs.size]||0)+1;
    whole.add(canon(f));
    return;
  }
  // Choose the current least flexible x from available y and y-x masks.
  let x=-1,Yx=0,best=14;
  for(let z=0;z<13;z++)if(U&bit(z)){
    const cand=Y&rotate(D,z), n=pop[cand];
    if(!n)return;
    if(n<best){x=z;Yx=cand;best=n;if(best===1)break;}
  }
  for(let y=0;y<13;y++)if(Yx&bit(y)){
    const q=lineVal[x][y];
    if(q.some(i=>counters[i]===5))continue;
    q.forEach(i=>counters[i]++);f[x]=y;
    graphs(U^bit(x),Y^bit(y),D^bit(mod(y-x)),n+1);
    q.forEach(i=>counters[i]--);
  }
}
graphs(rem,rem,rem,1);
check(hist[10]===1404&&hist[11]===60346&&whole.size===4750,"independent graph counts");

const triDirs=[[1,0],[0,1],[1,1]], stab=[];
for(let j=0;j<3;j++)for(let k=0;k<3;k++)if(j!==k){
  const r=triDirs[j],s=triDirs[k],t=triDirs[3-j-k];
  for(let a=1;a<13;a++)for(let b=1;b<13;b++){
    const M=[a*r[0],b*s[0],a*r[1],b*s[1]];
    if(mod((M[0]+M[1])*t[1]-(M[2]+M[3])*t[0])===0)stab.push(M);
  }
}
check(stab.length===72,"full linear stabilizer");
const keys=[...whole].sort(), forms=keys.map(parse), roots=[], seen=new Set();
for(let i=0;i<keys.length;i++)if(!seen.has(keys[i])){
  roots.push(i);
  for(const [a,b,c,d] of stab){
    const z=Array(13).fill(-1);
    forms[i].forEach((y,x)=>z[mod(a*x+b*y)]=mod(c*x+d*y));
    check(z.every(x=>x>=0),"stabilizer still graph");
    const k=canon(z);
    check(whole.has(k),"orbit missed graph");seen.add(k);
  }
}
check(roots.length===84&&seen.size===4750,"simultaneous orbits");
let ap=0;
for(let m=0;m<=all;m++)if(pop[m]===6){
  const a=F.filter(x=>m&bit(x));let yes=false;
  for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)if(m&bit(mod(7*(a[i]+a[j]))))yes=true;
  check(yes,"midpoint 3AP");ap++;
}
console.log(JSON.stringify({events,hist,translations:whole.size,stabilizer:stab.length,roots:roots.length,ap}));

function delta(low,high) {
  const h=new Uint16Array(13);
  for(let i=0;i<13;i++)for(let j=0;j<13;j++)
    h[mod(j-i)]|=bit(mod(high[j]-low[i]));
  return h;
}
function num(h){let n=0;for(const x of h)n+=pop[x];return n;}
function transform(raw,d,u=0,v=0) {
  const out=new Uint16Array(13), k=inv[d];
  for(let x=0;x<13;x++){
    const ymask=raw[x];
    const tx=mod(k*(x+u));
    for(let y=0;y<13;y++)if(ymask&bit(y))out[tx]|=bit(mod(k*(y+v)));
  }
  return out;
}
function join(h,g,sx=0,sy=0){
  const r=new Uint16Array(13);
  for(let x=0;x<13;x++)r[x]=h[x]|rotate(g[mod(x-sx)],sy);
  return r;
}
const triples=[];
let shapeTrials=0,allPairSmall=0,partial=0;
const adj=new Map();
for(const i of roots){
  const nb=new Map();adj.set(i,nb);
  for(let j=0;j<forms.length;j++){
    const h=delta(forms[i],forms[j]);
    if(num(h)<=89)nb.set(j,h);
  }
  for(const [j,h01] of nb)for(const [k,h02] of nb){
    shapeTrials++;
    const h12=delta(forms[j],forms[k]);if(num(h12)>89)continue;
    allPairSmall++;
    const q=transform(h02,2);
    for(const u of F)for(const v of F){
      const w=join(h01,q,7*u,7*v);
      if(num(w)>89)continue;partial++;
      const R=join(w,h12,u,v), n=num(R);
      if(n<=89)triples.push({i,j,k,u,v,n,R});
    }
  }
}
check(shapeTrials===47893&&allPairSmall===8371&&partial===253&&triples.length===2,"triple census");
let total=0,pass=0;
for(const T of triples){
  const {i,j,k,u,v,R}=T, counts=[0,0,0,0];
  for(let s=3;s<13;s++)for(const [l,d03] of adj.get(i)){
    const d13=delta(forms[j],forms[l]),d23=delta(forms[k],forms[l]);
    if(num(d13)>89||num(d23)>89)continue;
    // d23 is high minus low, so account for translating its lower endpoint.
    const q03=transform(d03,s), q13=transform(d13,s-1), q23=transform(d23,s-2,-u,-v);
    for(const x of F)for(const y of F){
      total++;counts[0]++;
      const R0=join(R,q03,inv[s]*x,inv[s]*y);
      if(num(R0)>89)continue;counts[1]++;
      const R1=join(R0,q13,inv[s-1]*x,inv[s-1]*y);
      if(num(R1)>89)continue;counts[2]++;
      const R2=join(R1,q23,inv[s-2]*x,inv[s-2]*y);
      if(num(R2)>89)continue;counts[3]++;pass++;
    }
  }
  T.counts=counts;
}
check(total===278850&&pass===0,"fourth-layer census");
const report=triples.map(({i,j,k,u,v,n,counts})=>({g:[i,j,k].map(x=>keys[x]),u,v,n,counts}));
console.log(JSON.stringify({shapeTrials,allPairSmall,partial,report,total,pass}));
