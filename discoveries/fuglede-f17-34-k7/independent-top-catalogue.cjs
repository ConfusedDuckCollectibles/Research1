'use strict';
const octic=(()=>{
'use strict';
const p=17, mod=x=>((x%17)+17)%17, pw=(a,b)=>{let r=1;while(b--)r=r*a%17;return r};
const IV=Array.from({length:p},(_,i)=>pw(i,15)), CAP=6;
const xs=Array.from({length:8},(_,i)=>i+1), sq=xs.map(x=>x*x%17), v4=sq.map(x=>pw(x,4));
const lifts=Array.from({length:p},()=>[]);
for(let l=0;l<=CAP;l++)for(let r=0;r<=CAP;r++){lifts[mod(l-r)].push({half:(l+r)*9%17,total:l+r});}
const best=[...lifts].map(L=>Math.min(...L.map(x=>x.total)));
const stats={oddFeasible:0,leaves:0,exact:0}, masks=new Map, catalogue=Array.from({length:17},()=>[]);
function lagrange(S){
 const id=S.join(',');if(masks.has(id))return masks.get(id);
 const nodes=[0,...S.map(i=>sq[i])];
 const out=sq.map(t=>nodes.map((s,j)=>
   nodes.reduce((a,u,k)=>k===j?a:a*mod(t-u)*IV[mod(s-u)]%17,1)));
 masks.set(id,out);return out;
}
for(let a=0;a<17;a++)for(let b=0;b<17;b++)for(let c=0;c<17;c++){
 const od=xs.map((x,i)=>x*(a+b*sq[i]+c*sq[i]*sq[i])%17), need=od.map(x=>2*x%17);
 const mins=need.map(x=>best[x]);
 if(mins.reduce((x,y)=>x+y,0)>34)continue;
 stats.oddFeasible++;
 const I=Array.from({length:8},(_,i)=>i).sort((i,j)=>lifts[need[i]].length-lifts[need[j]].length||j-i).slice(0,3).sort((a,b)=>a-b);
 const J=Array.from({length:8},(_,i)=>i).filter(i=>!I.includes(i)).sort((i,j)=>lifts[need[i]].length-lifts[need[j]].length);
 const bas=34-J.reduce((u,j)=>u+mins[j],0), lam=lagrange(I);
 const r=[0,0,0,0], ev=new Int8Array(8);
 for(let lead=1;lead<17;lead++){
  const d=v4.map(x=>x*lead%17);
  for(let z=0;z<=CAP&&z<=bas;z++){
   r[0]=z;
   for(const v1 of lifts[need[I[0]]]){
    r[1]=mod(v1.half-d[I[0]]);
    for(const v2 of lifts[need[I[1]]]){
     r[2]=mod(v2.half-d[I[1]]);
     for(const v3 of lifts[need[I[2]]]){
      const initial=z+v1.total+v2.total+v3.total;
      if(initial>bas)continue;
      r[3]=mod(v3.half-d[I[2]]);stats.leaves++;
      let cost=initial, i=0;
      for(;i<J.length;i++){
       const j=J[i], W=lam[j];
       ev[j]=(d[j]+W[0]*r[0]+W[1]*r[1]+W[2]*r[2]+W[3]*r[3])%17;
       const pp=(ev[j]+od[j])%17, nn=mod(ev[j]-od[j]);
       if(pp>CAP||nn>CAP)break;
       cost+=pp+nn;if(cost>34)break;
      }
      if(i<J.length||cost!==34)continue;
      ev[I[0]]=v1.half;ev[I[1]]=v2.half;ev[I[2]]=v3.half;
      const vals=Array(17).fill(z);
      xs.forEach((x,j)=>{vals[x]=(ev[j]+od[j])%17;vals[p-x]=mod(ev[j]-od[j]);});
      const cf=vals.map((v,j)=>j===0?z:mod(-vals.reduce((s,u,x)=>s+u*pw(x,16-j),0)));
      if(cf[8]!==lead||cf[7]||cf.slice(9).some(Boolean))throw Error('degree');
      catalogue[lead].push(cf.slice(0,9));stats.exact++;
     }
    }
   }
  }
 }
}
for(const F of catalogue.flat()){ const V=Array.from({length:17},(_,t)=>F.reduceRight((r,x)=>(r*t+x)%17,0)); if(V.some(x=>x>CAP)||V.reduce((s,t)=>s+t,0)!==34)throw Error('profile surface');}
if(stats.exact!==272)throw Error('octic count');
return {stats,perLeading:catalogue.map(row=>row.length)};


})();
const nonic=(()=>{
'use strict';
const p=17, mod=x=>((x%17)+17)%17, pw=(a,b)=>{let r=1;while(b--)r=r*a%17;return r};
const IV=Array.from({length:p},(_,i)=>pw(i,15)), CAP=6;
const xs=Array.from({length:8},(_,i)=>i+1), sq=xs.map(x=>x*x%17), v4=sq.map(x=>pw(x,4));
const lifts=Array.from({length:p},()=>[]);
for(let l=0;l<=CAP;l++)for(let r=0;r<=CAP;r++){lifts[mod(l-r)].push({half:(l+r)*9%17,total:l+r});}
const best=[...lifts].map(L=>Math.min(...L.map(x=>x.total)));
const stats={oddFeasible:0,leaves:0,exact:0}, masks=new Map, catalogue=[];
function lagrange(S){
 const id=S.join(',');if(masks.has(id))return masks.get(id);
 const nodes=[0,...S.map(i=>sq[i])];
 const out=sq.map(t=>nodes.map((s,j)=>
   nodes.reduce((a,u,k)=>k===j?a:a*mod(t-u)*IV[mod(s-u)]%17,1)));
 masks.set(id,out);return out;
}
for(let a=0;a<17;a++)for(let b=0;b<17;b++)for(let c=0;c<17;c++)for(let d=0;d<17;d++){
 const od=xs.map((x,i)=>x*(a+sq[i]*(b+sq[i]*(c+sq[i]*(d+sq[i]))))%17), need=od.map(x=>2*x%17);
 const mins=need.map(x=>best[x]);
 if(mins.reduce((x,y)=>x+y,0)>34)continue;
 stats.oddFeasible++;
 const I=Array.from({length:8},(_,i)=>i).sort((i,j)=>lifts[need[i]].length-lifts[need[j]].length||j-i).slice(0,3).sort((a,b)=>a-b);
 const J=Array.from({length:8},(_,i)=>i).filter(i=>!I.includes(i)).sort((i,j)=>lifts[need[i]].length-lifts[need[j]].length);
 const bas=34-J.reduce((u,j)=>u+mins[j],0), lam=lagrange(I);
 const r=[0,0,0,0], ev=new Int8Array(8);
 {
  const d=Array(8).fill(0);
  for(let z=0;z<=CAP&&z<=bas;z++){
   r[0]=z;
   for(const v1 of lifts[need[I[0]]]){
    r[1]=mod(v1.half-d[I[0]]);
    for(const v2 of lifts[need[I[1]]]){
     r[2]=mod(v2.half-d[I[1]]);
     for(const v3 of lifts[need[I[2]]]){
      const initial=z+v1.total+v2.total+v3.total;
      if(initial>bas)continue;
      r[3]=mod(v3.half-d[I[2]]);stats.leaves++;
      let cost=initial, i=0;
      for(;i<J.length;i++){
       const j=J[i], W=lam[j];
       ev[j]=(d[j]+W[0]*r[0]+W[1]*r[1]+W[2]*r[2]+W[3]*r[3])%17;
       const pp=(ev[j]+od[j])%17, nn=mod(ev[j]-od[j]);
       if(pp>CAP||nn>CAP)break;
       cost+=pp+nn;if(cost>34)break;
      }
      if(i<J.length||cost!==34)continue;
      ev[I[0]]=v1.half;ev[I[1]]=v2.half;ev[I[2]]=v3.half;
      const vals=Array(17).fill(z);
      xs.forEach((x,j)=>{vals[x]=(ev[j]+od[j])%17;vals[p-x]=mod(ev[j]-od[j]);});
      const cf=vals.map((v,j)=>j===0?z:mod(-vals.reduce((s,u,x)=>s+u*pw(x,16-j),0)));
      if(cf[9]!==1||cf[8]||cf.slice(10).some(Boolean))throw Error('degree');
      catalogue.push(cf.slice(0,10));stats.exact++;
     }
    }
   }
  }
 }
}
for(const F of catalogue){ const V=Array.from({length:17},(_,t)=>F.reduceRight((r,x)=>(r*t+x)%17,0)); if(V.some(x=>x>CAP)||V.reduce((s,t)=>s+t,0)!==34)throw Error('profile surface');}
if(stats.exact!==188)throw Error('nonic count');
return {stats,top: Array.from({length:17},(_,i)=>catalogue.filter(P=>P[7]===i).length)};


})();
console.log(JSON.stringify({octic,nonic},null,2));

