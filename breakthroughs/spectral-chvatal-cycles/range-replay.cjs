"use strict";
const assert=require("node:assert/strict");
function choose2(n){return n*(n-1)/2;}
function subsets(values,k) {
  if(!k)return [[]];
  return values.flatMap((value,i)=>subsets(values.slice(i+1),k-1).map(a=>[value,...a]));
}
const mask=arr=>arr.reduce((s,i)=>s|(1<<i),0);
function verify(cycles,n) {
  const r=cycles.reduce((s,k)=>s+k,0);
  assert(r>=5&&n>=2*r+1);
  const anchors=Array.from({length:r},(_,i)=>i);
  const tips=anchors.map(i=>i+r);
  const leaves=Array.from({length:n-2*r},(_,i)=>i+2*r);
  const tight=[...anchors,...leaves];
  const old=[...anchors,...tips,...leaves];
  const edges=new Set(),deleted=new Set();
  let offset=0;
  for(const k of cycles) {
    for(let i=0;i<k;i++){
      const a=offset+i,b=offset+(i+1)%k;
      edges.add(mask([a,b]));
      deleted.add(mask([a,b,r+a]));
    }
    offset+=k;
  }
  const U=subsets(old,3).map(mask).filter(a=>!deleted.has(a));
  const h=BigInt(choose2(n-4)),q=BigInt(choose2(n-2)),tipPairs=subsets(tips,2).map(mask);
  const p=BigInt(tipPairs.length),D=h*q*p;
  const has=(a,i)=>Boolean(a&(1<<i));
  const tipMask=mask(tips);
  const weight=BigInt(choose2(n-1));
  const ref=new Map();
  for(const i of tight)ref.set(i,tipPairs.map(t=>U.indexOf(t|(1<<i))));
  assert([...ref.values()].every(v=>v.every(j=>j>=0)));
  const sqrtBridge = U.map(a=>anchors.map(i=>{
    if(has(a,i))return 0n;
    let y=D/q;
    if(tipPairs.includes(a&tipMask)){
      const other=anchors.find(j=>has(a,j));
      if(other!==undefined)y+=BigInt(2-Number(edges.has(mask([i,other]))))*D/(q*p);
    }
    return y;
  }));
  for(const i of anchors)for(const j of tight){
    const v=U.reduce((s,a,u)=>s+(has(a,j)?sqrtBridge[u][i]:0n),0n);
    assert.equal(v,j===i?0n:D);
  }
  const oldEntry=(a,b)=>a===b?D:(a&b)?0n:D/h;
  const demands=U.map((a,u)=>tight.map(i=>D-(i<r?2n*sqrtBridge[u][i]:0n)-
    U.reduce((s,b)=>s+(has(b,i)?oldEntry(a,b):0n),0n)));
  const err=demands.map(row=>row.slice());
  const cols=new Map(tight.map((i,j)=>[i,j]));
  const x=Array.from({length:U.length},()=>new Map());
  function set(u,v,w){
    assert(!(U[u]&U[v]));
    assert(!x[u].has(v));
    if(w){
      x[u].set(v,w);
      x[v].set(u,w);
      for(const [i,j] of cols){
        if(has(U[v],i))err[u][j]-=w;
        if(has(U[u],i))err[v][j]-=w;
      }
    }
  }
  const refRows=new Set([...ref.values()].flat());
  for(let u=0;u<U.length;u++)if(!refRows.has(u))for(const [i,j] of cols)if(err[u][j]){
    const v=ref.get(i).find(k=>!(U[u]&U[k]));
    assert(v!==undefined);
    set(u,v,err[u][j]);
  }
  for(let i=0;i<tight.length;i++)for(let j=i+1;j<tight.length;j++){
    const a=tight[i],b=tight[j],left=ref.get(a),right=ref.get(b),V=[...left,...right];
    const parent=Array(V.length).fill(-1),walk=[0];parent[0]=0;
    for(let k=0;k<walk.length;k++){
      const u=walk[k],lo=u<p?Number(p):0,hi=u<p?V.length:Number(p);
      for(let v=lo;v<hi;v++)if(parent[v]<0&&!(U[V[u]]&U[V[v]])){
        parent[v]=u;walk.push(v);
      }
    }
    assert.equal(walk.length,V.length);
    const c=u=>u<left.length?j:i;
    assert.equal(left.reduce((s,u)=>s+err[u][j],0n),right.reduce((s,u)=>s+err[u][i],0n));
    for(const u of walk.slice(1).reverse())set(V[u],V[parent[u]],err[V[u]][c(u)]);
    assert.equal(err[V[0]][j],0n);
  }
  let nz=0;
  for(let u=0;u<U.length;u++){
    assert(err[u].every(x=>!x));
    for(const [i,j] of cols){
      const v=U.reduce((s,b,k)=>s+(has(b,i)?x[u].get(k)??0n:0n),0n);
      assert.equal(v,demands[u][j]);
    }
    nz+=x[u].size;
  }
  if (cycles.length===1&&r===5&&n===11) {
    const [a,b,c]=[97,386,400].map(v=>U.indexOf(v));
    const entry=(u,v)=>oldEntry(U[u],U[v])+(x[u].get(v)??0n);
    const u=entry(a,b), v=entry(a,c), w=entry(b,c);
    assert.equal(D*D*D-D*(u*u+v*v+w*w)+2n*u*v*w,-19125227520n);
  }
  return {cycles,n,triples:U.length,tight:tight.length,denominator:String(D),editedPairs:nz/2,target:String(weight)};
}
console.log(JSON.stringify([verify([5],11),verify([6],13),verify([3,3],13),verify([3,4],15)]));
