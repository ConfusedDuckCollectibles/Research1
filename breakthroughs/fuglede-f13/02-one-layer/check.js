"use strict";
const p=13,mod=x=>(x%13+13)%13,inv=[0,...Array.from({length:12},(_,i)=>Array.from({length:12},(_,j)=>j+1).find(j=>(i+1)*j%13===1))];
const pc=Array.from({length:8192},(_,i)=>{let c=0;for(;i;i&=i-1)c++;return c;});
const rot=(x,s)=>((x<<s)|(x>>(13-s)))&8191;
const check=(b,s)=>{if(!b)throw Error(s);};
const key=f=>f.map(x=>x.toString(16)).join("");
function canon(f){
 let best=null,ret;
 for(let a=0;a<p;a++){
  const g=Array.from({length:p},(_,i)=>mod(f[(a+i)%p]-f[a])),k=key(g);
  if(best===null||k<best){best=k;ret=g;}
 }
 return [best,ret];
}
const visits=Array(14).fill(0),hist={},forms=new Map(),f=[0];
function dfs(x,u,v,mask){
 visits[x]++;
 if(x===p){
  const j=pc[mask], k=key(f);hist[j]=(hist[j]??0)+1;
  if(!forms.has(k)){const [c,g]=canon(f);forms.set(c,g);}
  return;
 }
 for(let y=1;y<p;y++){
  const d=mod(y-x);if(u>>y&1||v>>d&1)continue;
  let s=mask,ok=true;const q=Array(p).fill(0);
  for(let i=0;i<x;i++){const r=mod((y-f[i])*inv[x-i]);if(++q[r]>=5){ok=false;break;}s|=1<<r;}
  if(ok){f.push(y);dfs(x+1,u|1<<y,v|1<<d,s);f.pop();}
 }
}
dfs(1,1,1,0);
const B=[...forms.values()],ids=new Map(B.map((f,i)=>[key(f),i]));
function orbit(g){
 const seen=new Set([key(g)]),Q=[g];
 for(let j=0;j<Q.length;j++){
  const f=Q[j],gs=[f.map((y,x)=>mod(x-y)),Array.from({length:p},(_,i)=>f.indexOf(i)),Array.from({length:p},(_,i)=>mod(2*f[mod(7*i)]))];
  for(const z of gs){const [k,a]=canon(z);if(!seen.has(k)){seen.add(k);Q.push(a);}}
 }
 return seen;
}
const roots=[],ignore=new Set();
for(const f of B){
 const k=key(f);if(ignore.has(k))continue;roots.push(ids.get(k));
 for(const v of orbit(f)){check(ids.has(v),"shape symmetry");ignore.add(v);}
}
check(visits[13]===61750&&B.length===4750&&roots.length===84&&hist[10]===1404&&hist[11]===60346,"graph census");
function diff(A,B){
 const H=Array(p).fill(0);
 for(let d=0;d<p;d++)for(let i=0;i<p;i++)H[d]|=1<<mod(A[(i+d)%p]-B[i]);
 return H;
}
function scale(H,k){
 const R=Array(p).fill(0);
 for(let x=0;x<p;x++)for(let y=0;y<p;y++)if(H[x]>>y&1)R[mod(x*k)]|=1<<mod(y*k);
 return R;
}
function shift(H,u,v){return Array.from({length:p},(_,x)=>rot(H[mod(x-u)],mod(v)));}
function join(...hs){return Array.from({length:p},(_,x)=>hs.reduce((s,h)=>s|h[x],0));}
const count=H=>H.reduce((s,x)=>s+pc[x],0);
let ap=0;for(let b=0;b<(1<<p);b++)if(pc[b]===6){
 let ok=false;for(let s=0;s<p;s++)if(b>>s&1)for(let d=1;d<p;d++)if((b>>mod(s+d)&1)&&(b>>mod(s+2*d)&1))ok=true;
 check(ok,"six-set contains a 3-AP");ap++;
}
const adjacencies=roots.map(j=>({
 j,adj:B.flatMap((g,i)=>{const H=diff(B[j],g);return count(H)<=89?[[i,H]]:[];})
}));
let trials=0,lastSmall=0,partial=0;const T=[];
for(const {j,adj} of adjacencies)for(const [a,H01] of adj)for(const [b,H02] of adj){
 trials++;const H12=diff(B[a],B[b]);if(count(H12)>89)continue;lastSmall++;
 const q=scale(H02,7);
 for(let u=0;u<p;u++)for(let v=0;v<p;v++){
  const R=join(H01,shift(q,-7*u,-7*v));if(count(R)>89)continue;partial++;
  const H=join(R,shift(H12,-u,-v)),n=count(H);
  if(n<=89)T.push({j,a,b,u,v,H,n,adj});
 }
}
check(T.length===2,"AP-triple controls");
let extensions=0,pass=0;
for(const T0 of T){
 const ins=[0,0,0,0];
 for(let s=3;s<p;s++)for(const [b,H03] of T0.adj){
  const f=B[b],hs=[H03,diff(B[T0.a],f),shift(diff(B[T0.b],f),T0.u,T0.v)];
  if(hs.some(x=>count(x)>89))continue;
  const w=hs.map((H,i)=>scale(H,inv[s-i]));
  for(let x=0;x<p;x++)for(let y=0;y<p;y++){
   extensions++;ins[0]++;
   const q0=join(T0.H,shift(w[0],-inv[s]*x,-inv[s]*y));if(count(q0)>89)continue;ins[1]++;
   const q1=join(q0,shift(w[1],-inv[s-1]*x,-inv[s-1]*y));if(count(q1)>89)continue;ins[2]++;
   const q2=join(q1,shift(w[2],-inv[s-2]*x,-inv[s-2]*y));if(count(q2)>89)continue;ins[3]++;
   pass++;
  }
 }
 T0.ins=ins;
}
check(extensions===278850&&pass===0,"no fourth layer");
console.log(JSON.stringify({visits,hist,translationShapes:B.length,firstShapes:roots.length,sixSets:ap,trials,lastSmall,partial,
 triples:T.map(T0=>({total:T0.n,forms:[T0.j,T0.a,T0.b].map(i=>B[i]),u:T0.u,v:T0.v,adj:T0.adj.length,ext:T0.ins})),extensions,pass}));
