"use strict";
const p=17,F=Array.from({length:p},(_,i)=>i),m=x=>(x%p+p)%p,iv=F.map(i=>F.find(j=>m(i*j)===1)??0);
const Q=[[1,0],[16,0],[0,1],[0,16],[1,1]],code=([x,y])=>17*m(x)+m(y), pts=F.flatMap(x=>F.map(y=>[x,y]));
const exps=F.slice(0,5).flatMap(d=>F.slice(0,d+1).map(i=>[i,d-i])),V=pts.map(([x,y])=>exps.map(([a,b])=>m(x**a*y**b)));
const bs=new Set(Q.flatMap(([x,y])=>Q.filter(([a,b])=>x!==a||y!==b).flatMap(([a,b])=>F.slice(1).map(k=>code([k*(x-a),k*(y-b)])))));
const reps=[...bs].filter(k=>{
 const [x,y]=pts[k];return k===Math.min(k,code([-x,-y]),code([y,x]),code([-y,-x]));
}).sort((a,b)=>a-b);
if(reps.length!==32)throw Error("pair symmetry");
function nullbasis(A,d){
 A=A.map(a=>a.slice());const piv=[];
 for(let i=0;i<d;i++){
  const j=A.findIndex((a,k)=>k>=piv.length&&a[i]);if(j<0)continue;
  const k=piv.length;[A[k],A[j]]=[A[j],A[k]];
  const v=iv[A[k][i]]; for(let t=i;t<d;t++)A[k][t]=m(A[k][t]*v);
  for(let r=0;r<A.length;r++)if(r!==k&&A[r][i]){const a=A[r][i];for(let t=i;t<d;t++) A[r][t]=m(A[r][t]-a*A[k][t]);}
  piv.push(i);
 }
 return F.slice(0,d).filter(i=>!piv.includes(i)).map(j=>F.slice(0,d).map(i=>i===j?1:piv.includes(i)?m(-A[piv.indexOf(i)][j]):0));
}
function idsAt(x,y,t){
 return Q.map(([u,v])=>code([(x+u)*iv[t],(y+v)*iv[t]]));
}
function allforms(B,fun){
 const v=new Int16Array(B[0]?.length??0);
 for(let j=0;j<B.length;j++){
  v.set(B[j]);
  function dfs(k){
   if(k===B.length){
    fun(v);return;
   }
   const q=B[k];
   for(let h=0;h<17;h++){
    dfs(k+1);
    for(let i=0;i<v.length;i++)v[i]=m(v[i]+q[i]);
   }
  }
  dfs(j+1);
 }
}
const key=v=>{const q=iv[v.find(x=>x)];let k=0;for(const t of v)k=17*k+m(q*t);return k;};
function replay(k){
 const [dx,dy]=pts[k],B=nullbasis(idsAt(0,0,1).concat(idsAt(dx,dy,1)).map(i=>V[i]),15);
 const H=V.map(r=>B.map(v=>m(r.reduce((s,x,i)=>s+x*v[i],0))));
 const h=Array(B.length+1).fill(0),A=new Map();
 let add=0;
 for(const t of F.slice(2))for(const [x,y] of pts){
  const J=idsAt(x,y,t).concat([code([x*iv[t-1],y*iv[t-1]]),code([(x-dx)*iv[t-1],(y-dy)*iv[t-1]])]);
  const C=nullbasis(J.map(i=>H[i]),B.length);h[C.length]++;
  allforms(C,v=>{const u=key(v);A.set(u,(A.get(u)??0)+1);add++;});
 }
 const counts={};let max=0;const rem=[];
 for(const [u,c]of A){counts[c]=(counts[c]??0)+1;max=Math.max(max,c);if(c>=22)rem.push([u,c]);}
 return {dx,dy,dim:B.length,h,add,distinct:A.size,max,counts,rem};
}
const chosen=reps.map(replay);

const all=chosen, summary=[];let peak=0,checks=0,found=[];
function decode(k,d){const a=Array(d);for(let i=d-1;i>=0;i--){a[i]=k%17;k=Math.floor(k/17);}if(k)throw Error("long");return a;}
for(const {dx,dy,dim,rem} of all)if(rem.length){
 const B=nullbasis(idsAt(0,0,1).concat(idsAt(dx,dy,1)).map(i=>V[i]),15), base=[[0,0,1],[dx,dy,1]], old=Q.map(([x,y])=>[m(-x),m(-y),0]);
 for(const [u,c]of rem){
  const a=decode(u,dim), f=V.map(v=>m(B.reduce((w,r,i)=>w+a[i]*m(v.reduce((s,x,j)=>s+x*r[j],0)),0))===0);
  const anchors=base.concat(old), pts3=F.slice(1).flatMap(t=>pts.map(([x,y])=>[x,y,t]));
  function compatible(A,B){
   const u=m(A[0]-B[0]),v=m(A[1]-B[1]),w=m(A[2]-B[2]);
   return w?f[code([u*iv[w],v*iv[w]])]:bs.has(code([u,v]));
  }
  function nc(U,V,W){
   const a=U.map((u,i)=>m(u-V[i])),b=U.map((u,i)=>m(u-W[i]));
   return m(a[0]*b[1]-a[1]*b[0])||m(a[0]*b[2]-a[2]*b[0])||m(a[1]*b[2]-a[2]*b[1]);
  }
  const pool=pts3.filter(P=>base.every(B=>!B.every((b,j)=>b===P[j]))&&anchors.every(B=>compatible(B,P))&&base.every((U,i)=>base.slice(i+1).every(V=>nc(U,V,P))));
  if(pool.filter(x=>x[2]!==1).length!==c)throw Error("placement count repeat");
  const adj=pool.map((A,i)=>pool.map((B,j)=>i!==j&&compatible(A,B)&&anchors.every(P=>nc(A,B,P))));
  const F0=pool.map((_,i)=>i),C=[], occ=new Uint8Array(17);occ[1]=2;
  let best=0,nodes=0;
  function color(P){
   let U=P.slice(),out=[],upper=[],i=0;
   while(U.length){
    const L=[],R=[];i++;
    for(const x of U){if(L.some(y=>adj[x][y]))R.push(x);else{L.push(x);out.push(x);upper.push(i);}}
    U=R;
   }
   return {out,upper};
  }
  function dfs(P){
   nodes++;
   best=Math.max(best,C.length+base.length);
   const {out,upper}=color(P);
   for(let j=out.length-1;j>=0;j--){
    if(C.length+base.length+upper[j]<=best)return;
    const x=out[j],t=pool[x][2];
    if(occ[t]===6)continue;
    const next=out.slice(0,j).filter(y=>adj[x][y]&&C.every(z=>nc(pool[z],pool[x],pool[y])));
    C.push(x);occ[t]++;dfs(next);occ[t]--;C.pop();
   }
  }
  dfs(F0);checks+=nodes;peak=Math.max(peak,best);if(best>=28)found.push([dx,dy,u,best]);
  summary.push([dx,dy,u,c,pool.length,best,nodes]);
 }
}
if(all.reduce((s,x)=>s+x.add,0)!==434232||summary.length!==70||peak!==6||checks!==21719||found.length)throw Error('expected no-go statistics');
console.log(JSON.stringify({all,peak,checks,found,summary}));
