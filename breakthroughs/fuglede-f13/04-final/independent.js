"use strict";
const mod=x=>(x%13+13)%13;
const inv=Array.from({length:13},(_,i)=>i===0?0:Array.from({length:13},(_,j)=>j).find(j=>i*j%13===1));
const pop=x=>x?1+pop(x&(x-1)):0;

function audit(mask){
 const s=Array.from({length:13},(_,i)=>i).filter(x=>mask>>x&1);
 if(s.length!==6||s[0]!==0||s[1]!==1)throw Error("gauge coordinates");
 const d=s.map(x=>s.map(y=>inv[mod(x-y)]));
 const byPalette=new Map(), visits=Array(7).fill(0), few={};
 const h=[0,0], models=[];
 function dfs(m){
  const k=h.length;visits[k]++;
  if(k===6){
   if(pop(m)<5){few[pop(m)]=(few[pop(m)]??0)+1;return;}
   models.push({h:h.slice(),D:m});
   // Restore the arbitrary slope h(1)-h(0) without moving the abscissae.
   for(let a=0;a<13;a++){
    const shift=(m<<a|m>>(13-a))&8191;
    const v=h.map((z,i)=>mod(z+a*s[i]));
    const z=byPalette.get(shift)??[];z.push(v);byPalette.set(shift,z);
   }
   return;
  }
  for(let z=0;z<13;z++){
   let m2=m;
   for(let j=0;j<k;j++)m2|=1<<mod((z-h[j])*d[k][j]);
   if(pop(m2)>5)continue;
   h.push(z);dfs(m2);h.pop();
  }
 }
 dfs(1);
 const num={},bad={};
 const ij=s.indexOf(3),kl=s.indexOf(5);
 const pairs=[];
 for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
  if([...byPalette.values()].every(A=>new Set(A.map(a=>mod(a[j]-a[i]))).size===1))
    pairs.push([s[i],s[j]]);
 }
 for(const [m,v] of byPalette){
  num[v.length]=(num[v.length]??0)+1;
  if(v.length>1&&new Set(v.map(a=>mod(a[kl]-a[ij]))).size!==1)bad[m]=v;
  // Fully reconstruct all fifteen pair directions, independently from recursion.
  for(const z of v){
   const mm=new Set;
   for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)mm.add(mod(inv[mod(s[j]-s[i])]*(z[j]-z[i])));
   if(mm.size!==5||[...mm].reduce((a,x)=>a|1<<x,0)!==m)throw Error("recovered palette");
  }
 }
 
 const h0=models[0].h,D0=models[0].D;
 const w=Array.from({length:13},(_,i)=>i).filter(t=>D0>>t&1);
 const affine=[];
 for(let a=1;a<13;a++)for(let b=0;b<13;b++){
  if(w.reduce((m,x)=>m|1<<mod(a*x+b),0)===D0)affine.push([a,b]);
 }
 return {mask,s,visits:visits.slice(2),normalized:models.length,low:few,palettes:byPalette.size,num,invariantPairs:pairs,bad,representative:h0,basePalette:w,stabilizer:affine};

}
for(const m of [111,187,215,311,427]){
 const r=audit(m);
 if(r.low[1]!==1||Object.keys(r.low).length!==1||Object.keys(r.bad).length||r.normalized!==12)
   throw Error("geometry gate");
 console.log(JSON.stringify(r));
}
