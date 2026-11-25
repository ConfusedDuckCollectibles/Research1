"use strict";
const p=13, mod=x=>(x%13+13)%13, inv=[0,...Array.from({length:12},(_,i)=>Array.from({length:12},(_,j)=>j+1).find(j=>(i+1)*j%13===1))];
const pop=n=>n?n%2+pop(n>>1):0;
const maps=[Array.from({length:14},(_,s)=>s===13?13:mod(s+1)),Array.from({length:14},(_,s)=>s===13?13:mod(s*2)),Array.from({length:14},(_,s)=>s===13?0:s===0?13:inv[s])];
const subsets=(n,k)=>Array.from({length:1<<n},(_,i)=>i).filter(i=>pop(i)===k);
function gate(cap){const grids=[];
 const result={q0:0,q1:0,q2:0,admissible:0}, counts={};
 const convert=(a,b,c)=>{
  const r=Array.from({length:13},(_,i)=>mod((a*i+b)*i+c));
  const z=r.map((v,i)=>v===0?i:-1).filter(v=>v!==-1);
  const req=(78-r.reduce((a,b)=>a+b,0))/13;
  if(req<0||!Number.isInteger(req)||req>z.length||(req>0&&cap<13)||r.some(v=>v>cap))return [];
  return subsets(z.length,req).map(m=>r.map((v,i)=>v+(v===0&&(m>>z.indexOf(i)&1)?13:0)));
 };
 const profiles=Array.from({length:13},(_,a)=>Array.from({length:13},(_,d)=>convert(a,0,d)));
 const rep=[], seen=new Set;
 for(const m of subsets(14,4))if(m<8192&&!seen.has(m)){rep.push(m);const queue=[m];seen.add(m);
  for(let j=0;j<queue.length;j++)for(const f of maps){const n=f.reduce((v,k,i)=>v|((queue[j]>>i&1)<<k),0);if(!seen.has(n)){seen.add(n);queue.push(n);}}
 }
 for(const m of rep){
  const s=Array.from({length:13},(_,i)=>i).filter(i=>m>>i&1);
  const w=s.map((t,i)=>inv[s.reduce((z,u,j)=>i===j?z:mod(z*(t-u)),1)]);
  const arg=Array.from({length:169},(_,u)=>s.map(t=>mod(u%13-t*(u/13|0))));
  function evaluate(mat, deg){
   result.admissible++;
   const f=[];
   for(const fs of arg){const n=mat[0][fs[0]]+mat[1][fs[1]]+mat[2][fs[2]]+mat[3][fs[3]]-18;
    if(n!==0&&n!==13&&n!==26)return;
    f.push(n/13);
   }
   const assert=(p)=>{if(!p)throw Error("recovered profile verification");};
   for(let t=0;t<14;t++){
    const r=Array(13).fill(0), i=s.indexOf(t);
    f.forEach((v,z)=>r[t===13?(z/13|0):mod(z%13-t*(z/13|0))]+=v);
    assert(r.every((v,j)=>v===(i<0?6:mat[i][j])));
   }
   const peak=Math.max(...mat.flat()), doubles=f.filter(i=>i===2).length, k=mat.map(r=>r.reduce((s,i)=>s+(i-6)**2,0)/26).sort((a,b)=>a-b);
   const peaks=mat.map(row=>Math.max(...row)).sort((a,b)=>a-b);
   const key=[peak,doubles,...k,...peaks].join(',');
   counts[key]=(counts[key]??0)+1;
   result[deg]++;if(doubles===0)grids.push({s:s.slice(),f});
  }
  // All four quadratic coefficients nonzero, common center translation.
  for(let b=1;b<13;b++){
   const lists=w.map(x=>profiles[mod(b*x)]);
   for(let d=0;d<13;d++)if(lists[0][d].length)for(let e=0;e<13;e++)if(lists[1][e].length)for(let f=0;f<13;f++)if(lists[2][f].length){
    const h=mod(18-d-e-f);if(!lists[3][h].length)continue;
    for(const A of lists[0][d])for(const B of lists[1][e])for(const C of lists[2][f])for(const D of lists[3][h])evaluate([A,B,C,D],"q2");
   }
  }
  if(cap<12)continue;
  for(let a=0;a<13;a++)for(let b=0;b<13;b++)if(a||b){
   const L=w.map((z,i)=>mod(z*(a+b*s[i]))), d=Array(4).fill(0);
   const ids=[0,1,2,3].filter(i=>L[i]);
   const rest=[0,1,2,3].filter(i=>!ids.slice(0,2).includes(i));
   const r0=rest[0],r1=rest[1];
   for(let e=0;e<(L[r0]?13:1);e++){
    d[r0]=e;d[r1]=mod(18-e);
    if(!L[r1]&&d[r1]!==0)continue;
    if(!L[r0]&&d[r0]!==0)continue;
    if(L.some((a,i)=>a===0)&&cap<13)continue;
    const mm=L.map((a,i)=>convert(0,a,d[i]));
    if(mm.some(a=>!a.length))continue;
    for(const A of mm[0])for(const B of mm[1])for(const C of mm[2])for(const D of mm[3])evaluate([A,B,C,D],"q1");
   }
  }
 }
 return {grids,orbits:rep.map(m=>Array.from({length:13},(_,i)=>i).filter(i=>m>>i&1)),result,counts};
}




function canonMask(t){
 let n=8192;
 for(let a=1;a<13;a++)for(let b=0;b<13;b++)n=Math.min(n,t.reduce((z,x)=>z|1<<mod(a*x+b),0));
 return n;
}
function planarShadows(n,budget,cap){
 const reps=new Set, projections=new Map;
 let count=0;
 const a=[[0,0],[1,1]];
 function check(mask,d){
  if(a.length===n){
   count++;const s=a.map(([x,y])=>13*x+y).join(',');if(reps.has(s))return;reps.add(s);
   for(let q=0;q<=13;q++){
    const x=a.map(([x,y])=>q<13?mod(y-q*x):x);
    if(new Set(x).size!==n)continue;
    const c=canonMask(x);if(!projections.has(c))projections.set(c,new Set);
    projections.get(c).add(pop(d));
   }
   return;
  }
  for(let x=a[a.length-1][0]+1;x<=13-(n-a.length);x++)for(let y=2;y<13;y++){
   if(a.some(([u,v])=>v===y))continue;
   const used=Array(13).fill(0);let k=d,ok=true;
   for(const [u,v] of a){const q=mod((y-v)*inv[mod(x-u)]);if(!(mask>>q&1)||++used[q]>=cap){ok=false;break;}k|=1<<q;}
   if(ok){a.push([x,y]);check(mask,k);a.pop();}
  }
 }
 for(let m=2;m<8192;m+=4)if(pop(m)===budget)check(m,2);
 return {count,models:reps.size,proj:[...projections.entries()].map(([k,v])=>[k,[...v].sort((a,b)=>a-b)]).sort(([a],[b])=>a-b)};
}
function profile(grid){
 const rows=Array.from({length:14},()=>Array.from({length:13},()=>[]));
 for(let i=0;i<169;i++)if(grid[i]){
  const x=i/13|0,y=i%13;
  for(let q=0;q<14;q++)rows[q][q===13?x:mod(y-q*x)].push(q===13?y:x);
 }
 return rows;
}
function cubicProfileGate(){
 const weights={};let profiles=0,layer=0,max=0,countMax=0;
 for(let a=0;a<13;a++)for(let b=0;b<13;b++)for(let c=0;c<13;c++)for(let d=0;d<13;d++){
  const row=Array.from({length:13},(_,t)=>mod(((a*t+b)*t+c)*t+d));
  const sum=row.reduce((s,x)=>s+x),z=row.filter(x=>x===0).length,k=(78-sum)/13;
  if(!Number.isInteger(k)||k<0||k>z)continue;
  const multiplicity=subsets(z,k).length;
  if(!(a||b||c||d)){layer+=multiplicity;continue;}
  const n=row.reduce((s,x)=>s+(x-6)**2,0)+13*k;
  if(n%26)throw Error("cubic energy integrality");
  const w=n/26;profiles+=multiplicity;max=Math.max(max,w);
  weights[w]=(weights[w]??0)+multiplicity;
 }
 return {profiles,layer,max,weights};
}
function must(x,msg){if(!x)throw Error(msg);}
function main(){
 const six=planarShadows(6,5,5);
 const zero=[[6,4,5],[10,7,5],[11,9,5],[12,9,5]];
 const negative=zero.map(a=>[a,planarShadows(...a).models]);
 const sixSet=new Set(six.proj.map(r=>r[0]));
 must(negative.every(r=>r[1]===0),"planar negative");
 must(six.count===840&&JSON.stringify([...sixSet])==="[111,187,215,311,427]","six shadow mask");
 console.log(JSON.stringify({negative,six}));
 function labels(f){
  const pp=profile(f);
  return {peaks:pp.map(r=>Math.max(...r.map(t=>t.length))),
  excluded:pp.map(r=>r.filter(t=>t.length===6&&!sixSet.has(canonMask(t))).length)};
 }
 const low=gate(11),high=gate(13);
 must(low.result.q1+low.result.q2===0,"cap11 high pencil");
 must(high.result.q1===0&&high.result.q2===804&&high.grids.length===24,"K10 atlas");
 const keys={};let min=13;
 for(const b of high.grids){
  must(b.f.every(x=>x===0||x===1)&&b.f.reduce((a,b)=>a+b)===78,"Boolean mass");
  const q=labels(b.f);
  const k=q.peaks.slice().sort((a,b)=>a-b).join(",");
  must(k==="6,6,6,6,6,6,6,6,6,6,10,10,10,12","K10 peaks");
  for(let r=0;r<14;r++)if(q.peaks[r]===6){
   must(q.excluded[r]>0,"K10 labelled bound");min=Math.min(min,q.excluded[r]);
  }
  const key=JSON.stringify(q);keys[key]=(keys[key]??0)+1;
 }
 console.log(JSON.stringify({cap11:low.result,q1:high.result.q1,q2:high.result.q2,boolean:high.grids.length,minExcludedSix:min,profiles:Object.entries(keys).map(([s,n])=>[JSON.parse(s),n])}));
 const tri=Array.from({length:169},(_,i)=>+(Math.floor(i/13)+i%13>=13));
 const rr=labels(tri);
 must(rr.peaks.filter(n=>n===12).length===3&&rr.peaks.filter(n=>n===6).length===11,"K11 peaks");
 must(rr.peaks.every((p,i)=>p!==6||rr.excluded[i]>0),"K11 labelled bound");
 console.log(JSON.stringify({triangle:rr,k10_LB:94,k11_LB:96,k10_Umax:89}));
 const c=cubicProfileGate();
 must(c.profiles===23765&&c.layer===1716&&c.max===12,"cubic maximum");
 console.log(JSON.stringify({profiles:c,k11Deficit:129,k11Factors:11,needDivisibility:489%7}));
}
main();
