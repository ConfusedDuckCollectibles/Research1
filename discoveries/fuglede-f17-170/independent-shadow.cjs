'use strict';
const p=17, mod=x=>(x%17+17)%17, fixed=[1,3], q=(t,a,c)=>mod(a*t*t+c);
const source=Array.from({length:17},()=>[]);
const groups={};
for(let a=1;a<17;a++){
 let out=[];
 for(let c=0;c<17;c++){
  const u=Array.from({length:17},(_,t)=>q(t,a,c));
  const z=u.flatMap((v,i)=>v===0?[i]:[]);
  for(let t=0;t<(1<<z.length);t++){
   const w=u.slice();for(let j=0;j<z.length;j++)if(t&(1<<j))w[z[j]]=17;
   if(w.reduce((x,y)=>x+y,0)!==119)continue;
   const sm=w.reduce((x,y)=>x+(y-7)**2,0);
   if(sm%34)throw Error('normalization');
   out.push({a,c,u:w,score:sm/34});
  }
 }
 source[a]=out;
 if(fixed.includes(a))groups[a]=out.map(w=>[w.c,w.score]);
}
// Enumerate all 240 projective pencil shapes/scales, then all variance-compatible profiles.
let total=0, cpass=0;
for(let r=2;r<17;r++)for(let k=1;k<17;k++){
 const roots=[r*(1-r),r-1,-r,1].map(x=>source[mod(k*x)]);
 for(const a of roots[0])for(const b of roots[1])for(const c of roots[2])for(const z of roots[3]){
   if(a.score+b.score+c.score+z.score!==35)continue;
   total++;
   if(mod(a.c+b.c+c.c+z.c)!==4)continue;
   cpass++;
 }
}
if(total!==96||cpass!==0)throw Error('shadow gate');
console.log(JSON.stringify({groups,variance:total,constant:cpass},null,2));
