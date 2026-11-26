'use strict';
const p=17, mod=x=>(x%17+17)%17, pow=(a,b)=>{let z=1;while(b--)z=z*a%17;return z}, inv=Array.from({length:17},(_,x)=>pow(x,15));
const pc=n=>{let k=0;while(n){n&=n-1;k++;}return k;};
const perms=[s=>s===17?s:(s+1)%17,s=>s===17?s:s*3%17,s=>s===17?0:s===0?17:inv[s]].map(f=>Array.from({length:18},(_,i)=>f(i)));
const seen=new Uint8Array(1<<18), reps=[], sizes=[];
for(let m=0;m<(1<<18);m++){
 if(seen[m]||pc(m)!==9)continue;
 const q=[m];seen[m]=1;
 for(let j=0;j<q.length;j++)for(const pi of perms){
  let n=0;for(let s=0;s<18;s++)if(q[j]&(1<<s))n|=1<<pi[s];
  if(!seen[n]){seen[n]=1;q.push(n);}
 }
 reps.push(Math.min(...q));sizes.push(q.length);
}
if(sizes.reduce((a,b)=>a+b,0)!==48620)throw Error('PGL coverage');
const N=289, bit=Array.from({length:N},(_,i)=>1n<<BigInt(i)), id=new Map(bit.map((b,i)=>[b,i]));
const L=Array.from({length:N},(_,i)=>Array.from({length:18},(_,s)=>s===17?Math.floor(i/17):mod(i%17-s*Math.floor(i/17))));
const lines=Array.from({length:18},()=>Array(17).fill(0n));
for(let i=0;i<N;i++)for(let s=0;s<18;s++)lines[s][L[i][s]]|=bit[i];
const all=(1n<<289n)-1n;
const result=[];
for(const mask of reps){
 const G=Array.from({length:18},(_,i)=>i).filter(i=>!(mask&(1<<i))),B=Array.from({length:18},(_,i)=>i).filter(i=>mask&(1<<i));
 const adj=L.map(c=>B.reduce((a,s)=>a&~lines[s][c[s]],all));
 const count=G.map(()=>new Uint8Array(17));
 const path=[0];for(let g=0;g<G.length;g++)count[g][0]=1;
 let nodes=0, branches=0, maxdepth=1, witness=null;
 function visit(P){
  nodes++;
  const n=path.length;if(n>maxdepth)maxdepth=n;
  if(n===12){witness=path.slice();return true;}
  let Q=P;
  const order=[], colors=[];
  // Greedy proper coloring gives an upper bound on remaining clique size.
  let c=0;
  while(Q){
   c++;let R=Q;
   while(R){
    const x=R&-R,u=id.get(x);
    order.push(u);colors.push(c);Q&=~x;R&=~(x|adj[u]);
   }
  }
  for(let j=order.length-1;j>=0;j--){
   if(colors[j]+n<12)return false;
   const u=order[j];P&=~bit[u];
   let R=P&adj[u];
   for(let g=0;g<G.length;g++){
    const s=G[g],l=L[u][s];
    if(++count[g][l]===5)R&=~lines[s][l];
   }
   path.push(u);branches++;
   const r=visit(R);path.pop();
   for(let g=0;g<G.length;g++)count[g][L[u][G[g]]]--;
   if(r)return true;
  }
  return false;
 }
 visit(adj[0]);result.push({mask,nodes,branches,maxdepth,witness});
}
if(result.some(x=>x.witness))throw Error(JSON.stringify(result));
console.log(JSON.stringify({reps,sizes,result},null,2));
