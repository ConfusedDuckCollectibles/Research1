"use strict";
const p=13,M=(x)=>(x%p+p)%p,inv=[0,1,7,9,10,8,11,2,5,3,4,6,12];
const pp=Array.from({length:p},(_,x)=>Array.from({length:p},(_,y)=>x*p+y)).flat();
const edge=Array.from({length:169},(_,k)=>Array.from({length:169},(_,j)=>{
 const dx=M(Math.floor(k/p)-Math.floor(j/p)),dy=M(k%p-j%p);
 return !dx?13:M(dy*inv[dx]);
}));
const pop=Array.from({length:1<<14},(_,i)=>i===0?0:i.toString(2).replace(/0/g,"").length);
function planar(n,d){
 const Q=new Int16Array(n);Q[0]=0;Q[1]=14;
 // Track actual direction set only, not authored direction masks. To prune
 // early, allow at most three hits on a parallel line through the new point.
 let ans=0,vis=Array(n+1).fill(0);
 function rec(k,start,usedY,mask){
  vis[k]++;
  if(k===n){ans++;return;}
  for(let x=start;x<=p-(n-k);x++)for(let y=2;y<p;y++){
   if((usedY>>y)&1)continue;
   const i=x*p+y;let once=0,twice=0,mm=mask;
   for(let q=0;q<k;q++){
    const j=1<<edge[i][Q[q]];
    if(twice&j){mm=-1;break;} // four collinear including the new point
    twice|=once&j; once|=j; mm|=j;
   }
   if(mm!==-1&&pop[mm]<=d){Q[k]=i;rec(k+1,x+1,usedY|1<<y,mm);}
  }
 }
 rec(2,2,3,2);
 if(ans!==0)throw Error("planar counterexample");
 return {n,d,visits:vis.slice(2),hits:ans};
}
const out=[[7,5],[8,7],[10,9],[12,10]].map(k=>planar(...k));
console.log(JSON.stringify({out},null,2));