'use strict';
const p=17, md=x=>(x%17+17)%17;
const pw=(x,n)=>{let y=1;while(n--)y=y*x%17;return y;};
const iv=Array.from({length:17},(_,i)=>pw(i,15));
const pop=x=>{let n=0;while(x){x&=x-1;n++;}return n;};
const perm=[x=>x===17?17:md(x+1),x=>x===17?17:x*3%17,x=>x===17?0:x===0?17:iv[x]].map(f=>Array.from({length:18},(_,i)=>f(i)));
const seen=new Uint8Array(1<<18), reps=[], orbitSizes=[];
for(let k=0;k<(1<<18);k++){
 if(seen[k]||pop(k)!==5)continue;
 const q=[k];seen[k]=1;
 for(let i=0;i<q.length;i++)for(const a of perm){
  let m=0;for(let s=0;s<18;s++)if(q[i]&(1<<s))m|=1<<a[s];
  if(!seen[m]){seen[m]=1;q.push(m);}
 }
 // Missing horizontal, vertical, slope one, and two generic missing slopes.
 const r=Math.min(...q.filter(m=>(m&((1<<17)|3))===((1<<17)|3)));
 if(!Number.isFinite(r))throw Error('coverage');
 reps.push(r);orbitSizes.push(q.length);
}
if(orbitSizes.reduce((a,s)=>a+s,0)!==8568)throw Error('orbit census');
const all=new Uint32Array(10);for(let i=0;i<289;i++)all[i>>>5]|=1<<(i&31);
const baseRows=Array.from({length:18},()=>Array.from({length:17},()=>new Uint32Array(10)));
const tt=Array.from({length:289},(_,i)=>Array.from({length:18},(_,s)=>s===17?Math.floor(i/17):md(i%17-s*Math.floor(i/17))));
for(let i=0;i<289;i++)for(let s=0;s<18;s++)baseRows[s][tt[i][s]][i>>>5]|=1<<(i&31);
function solve(mask,cap){
 const bad=Array.from({length:18},(_,i)=>i).filter(i=>mask&(1<<i));
 const good=Array.from({length:18},(_,i)=>i).filter(i=>!(mask&(1<<i)));
 const rows=bad.map(d=>baseRows[d]), pointB=tt.map(a=>bad.map(d=>a[d])), pointG=tt.map(a=>good.map(d=>a[d]));
 const counts=good.map(()=>new Uint8Array(17)), used=new Uint32Array(5);
 const avail=Array.from({length:18},()=>new Uint32Array(10)), ranks=new Uint8Array(18), points=[];
 avail[0].set(all);
 let nodes=0,complete=0,maxDepth=0,first=null;const depthCount=Array(18).fill(0);
 function use(k,v){
  const prev=avail[k],next=avail[k+1],bs=pointB[v],gs=pointG[v];
  next.set(prev);
  for(let b=0;b<5;b++){
   used[b]|=1<<bs[b];const rr=rows[b][bs[b]];
   for(let j=0;j<10;j++)next[j]&=~rr[j];
  }
  for(let g=0;g<good.length;g++)if(++counts[g][gs[g]]===cap){
   const rr=baseRows[good[g]][gs[g]];
   for(let j=0;j<10;j++)next[j]&=~rr[j];
  }
  points.push(v);
  visit(k+1);
  points.pop();
  for(let b=0;b<5;b++)used[b]^=1<<bs[b];
  for(let g=0;g<good.length;g++)counts[g][gs[g]]--;
 }
 function visit(k){
  nodes++;depthCount[k]++;
  if(k>maxDepth)maxDepth=k;
  if(k===17){complete++;if(first===null)first=[...points];return;}
  const V=avail[k];let wh=-1, wi=-1, mn=290;
  for(let b=0;b<5;b++){
   const rem=(1<<17)-1 & ~used[b];
   for(let l=rem;l;l&=l-1){
    const i=31-Math.clz32(l&-l), R=rows[b][i];
    let wt=0;
    for(let j=0;j<10&&wt<mn;j++)wt+=pop(V[j]&R[j]);
    if(!wt)return;
    if(wt<mn){mn=wt;wh=b;wi=i;}
   }
  }
  if(wh<0)throw Error('empty min');
  const R=rows[wh][wi];
  for(let j=0;j<10;j++)for(let x=V[j]&R[j];x;x&=x-1)use(k,j*32+31-Math.clz32(x&-x));
 }
 use(0,0);
 return {mask,bad,cap,nodes,complete,maxDepth,first,depthCount};
}
const p5=reps.map(m=>solve(m,5)); if(p5.some(x=>x.complete))throw Error('counterexample');
const control=solve(reps[0],17); if(control.complete<13)throw Error('positive control');
console.log(JSON.stringify({reps,orbitSizes,p5,control},null,2));
