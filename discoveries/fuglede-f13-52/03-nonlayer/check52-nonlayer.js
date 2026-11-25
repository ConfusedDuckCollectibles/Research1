"use strict";

function plane10Shadow(){
 const profiles=new Map, zero=0;
 let size=0, normalized=0;
 function canon(arr){
  let best=null;
  for(let a=1;a<13;a++)for(let b=0;b<13;b++) {
   let z='';for(let x=0;x<13;x++)z+=arr[(a*x+b)%13];
   if(best===null||z<best)best=z;
  }
  return best;
 }
 const sets=new Set;
 for(let m=2;m<8192;m+=4)if(pop(m)===9) {
  const a=[0,14];
  function dfs(next,used,dmask) {
   if(a.length===10) {
    const sk=a.join(","); if(sets.has(sk))return;
    sets.add(sk);normalized++;
    for(let d=0;d<14;d++){
     const r=Array(13).fill(0);
     for(const z of a)r[d===13?(z/13|0):mod(z%13-d*(z/13|0))]++;
     if(Math.max(...r)>2)continue;
     const c=canon(r),k=profiles.get(c)??{c,n:0,sample:[a.slice(),d]};
     k.n++;profiles.set(c,k);
    }
    return;
   }
   for(let x=next;x<=13-(10-a.length);x++)for(let y=2;y<13;y++)if(!(used>>y&1)){
    const b=13*x+y,c=Array(13).fill(0);let nd=dmask;
    if(a.some(u=>{const k=slope[u][b];nd|=1<<k;return !(m>>k&1)||++c[k]>=4;}))continue;
    a.push(b);dfs(x+1,used|(1<<y),nd);a.pop();
   }
  }
  dfs(2,3,2);
 }
 return {normalized,profiles:[...profiles.values()].sort((a,b)=>a.c.localeCompare(b.c)),canon};
}
function planeShadow(goal,maxD,cap,minDir=0){
 const profiles=new Map, zero=0;
 let size=0, normalized=0;
 function canon(arr){
  let best=null;
  for(let a=1;a<13;a++)for(let b=0;b<13;b++) {
   let z='';for(let x=0;x<13;x++)z+=arr[(a*x+b)%13];
   if(best===null||z<best)best=z;
  }
  return best;
 }
 const sets=new Set;
 for(let m=2;m<8192;m+=4)if(pop(m)===maxD) {
  const a=[0,14];
  function dfs(next,used,dmask) {
   if(a.length===goal) {
    const sk=a.join(","); if(sets.has(sk)||pop(dmask)<minDir)return;
    sets.add(sk);normalized++;
    for(let d=0;d<14;d++){
     const r=Array(13).fill(0);
     for(const z of a)r[d===13?(z/13|0):mod(z%13-d*(z/13|0))]++;
     if(Math.max(...r)>1)continue;
     const c=canon(r),k=profiles.get(c)??{c,n:0,sample:[a.slice(),d]};
     k.n++;profiles.set(c,k);
    }
    return;
   }
   for(let x=next;x<=13-(goal-a.length);x++)for(let y=2;y<13;y++)if(!(used>>y&1)){
    const b=13*x+y,c=Array(13).fill(0);let nd=dmask;
    if(a.some(u=>{const k=slope[u][b];nd|=1<<k;return !(m>>k&1)||++c[k]>=cap;}))continue;
    a.push(b);dfs(x+1,used|(1<<y),nd);a.pop();
   }
  }
  dfs(2,3,2);
 }
 return {normalized,profiles:[...profiles.values()].sort((a,b)=>a.c.localeCompare(b.c)),canon};
}const p=13,mod=n=>(n%13+13)%13,pop=n=>n?n%2+pop(n>>1):0;
const inv=Array.from({length:13},(_,i)=>i?Array.from({length:13},(_,j)=>j).find(j=>i*j%13===1):0);
let maps=[Array.from({length:14},(_,i)=>i===13?13:mod(i+1)),Array.from({length:14},(_,i)=>i===13?13:mod(2*i)),Array.from({length:14},(_,i)=>i===13?0:i===0?13:inv[i])];
const slope=Array.from({length:169},(_,a)=>Array.from({length:169},(_,b)=>((a%13-b%13+13)*inv[(Math.floor(a/13)-Math.floor(b/13)+13)%13])%13));
const planarGate=(goal,directions,cap,one=false)=>{
const masks=[]; for(let m=2;m<8192;m+=4)if(pop(m)===directions)masks.push(m);
const visits=Array(goal+1).fill(0);let model=[];
for(const mask of masks) {
const pts=[0,14];
function search(nextx,yused){
const n=pts.length;visits[n]++;if(n===goal){if(!model.length)model=pts.slice();return one;}
for(let x=nextx;x<=p-(goal-n);x++)for(let y=2;y<p;y++)if(!((yused>>y)&1)){
 const b=x*13+y, u=Array(13).fill(0);
 if(pts.some(a=>{const d=slope[a][b];return !((mask>>d)&1)||++u[d]>=cap;}))continue;
 pts.push(b);const stop=search(x+1,yused|(1<<y));pts.pop();if(stop)return true;
}
}
if(search(2,3))break;
}
return {goal,directions,cap,masks:masks.length,visits:visits.slice(2),model};
};


function orbitIDs(n) {
 const to=Array(1<<14).fill(-1), reps=[];
 for(let z=0;z<1<<14;z++) if(pop(z)===n&&to[z]===-1){
  let c=reps.length; reps.push(z);to[z]=c;
  let a=[z];
  for(let i=0;i<a.length;i++)for(const g of maps){
   let k=g.reduce((v,t,j)=>v|(((a[i]>>j)&1)<<t),0);
   if(to[k]<0){a.push(k);to[k]=c;}
  }
 }
 return {to,reps};
}
function stats(goal,maxD,cap) {
 const table=Array.from({length:15},(_,i)=>i===0?null:orbitIDs(i));
 const peaks={};
 let controls={};
 for(let m=2;m<8192;m+=4)if(pop(m)===maxD) {
  const a=[0,14];
  function dfs(next,used,dmask) {
   if(a.length===goal) {
    let d=pop(dmask),t=table[d].to[dmask],k=d+":"+t;
    peaks[k]=(peaks[k]??0)+1;
    controls[k]??=a.slice();return;
   }
   for(let x=next;x<=13-(goal-a.length);x++)for(let y=2;y<13;y++)if(!(used>>y&1)){
    const b=13*x+y,c=Array(13).fill(0);let nd=dmask;
    if(a.some(u=>{const k=slope[u][b];nd|=1<<k;return !(m>>k&1)||++c[k]>=cap;}))continue;
    a.push(b);dfs(x+1,used|(1<<y),nd);a.pop();
   }
  }
  dfs(2,3,2);
 }
 return {goal,maxD,cap,peaks,controls};
}
function cRow(row){let best=null;for(let a=1;a<13;a++)for(let b=0;b<13;b++){let z='';for(let x=0;x<13;x++)z+=row[(a*x+b)%13];if(best===null||z<best)best=z;}return best;}
function passPeaks(grid,s) {
 let sz=0;
 for(const d of s)for(let j=0;j<13;j++){
 const row=Array.from({length:13},(_,x)=>grid[x*13+mod(j+d*x)]);
 if(row.reduce((a,b)=>a+b,0)!==10)continue;
 sz++;
 const c=cRow(row);if(!allow10.has(c))return {good:false,c};
 }
 return {good:true,sz};
}function run(cap) {
 const pr=Array.from({length:2197},()=>[]);
 for(let a=0;a<13;a++) for(let b=0;b<13;b++) for(let c=0;c<13;c++){
 const q=pr[a*169+b*13+c];if(!a&&!b&&!c)continue;
 for(let d=0;d<13;d++){
 const v=Array.from({length:13},(_,i)=>mod(((a*i+b)*i+c)*i+d));
 const sum=v.reduce((a,b)=>a+b,0);
 if(v.some(x=>x>cap)||sum!==52)continue;
 q.push({v,d,peak:Math.max(...v)});
 }}
 let seen=new Set(),reps=[];
 for(let m=0;m<8192;m++)if(pop(m)===5&&!seen.has(m)){
 const S=[];for(let i=0;i<13;i++)if(m>>i&1)S.push(i);
 const queue=[m];seen.add(m);for(let j=0;j<queue.length;j++)for(const f of maps){
 const k=f.reduce((a,b,i)=>a|(((queue[j]>>i)&1)<<b),0);
 if(!seen.has(k)){queue.push(k);seen.add(k);}
 }
 const w=S.map((s,i)=>inv[S.reduce((a,t,j)=>i===j?a:mod(a*(s-t)),1)]);
 const arg=Array.from({length:169},(_,q)=>S.map(t=>mod(q%13-t*(q/13|0))));
 const stats={S,q3:0,q2:0,mod:0,exact:0,sigs:{},samples:[]};reps.push(stats);
 function test(a,b,c,u,v,z,deg){
 const lists=S.map((s,i)=>pr[169*mod(w[i]*a)+13*mod(w[i]*(b+c*s))+mod(w[i]*(u+v*s+z*s*s))]);if(lists.some(x=>!x.length))return;
 const g=[];
 function dfs(i,total) {
 if(i===5){
  if(mod(total-16))return;stats.mod++;
  const f=[];
  for(const z of arg){
    let v=g.reduce((a,q,j)=>a+q.v[z[j]],0)-16;
    if(v<0||v%13||v>13*Math.floor(5*(cap-4)/13))return;
    f.push(v/13);
  }

  for(let t=0;t<14;t++){
    const r=Array(13).fill(0), j=S.indexOf(t);
    f.forEach((v,q)=>r[t===13?(q/13|0):mod(q%13-t*(q/13|0))]+=v);
    if(!r.every((v,i)=>v===(j<0?4:g[j].v[i])))throw Error("profile recovery");
  }
  stats.exact++;stats[deg]++;
  const k=g.map(x=>x.v.reduce((s,n)=>s+(n-4)**2,0)/26).sort((a,b)=>a-b),counts=[f.filter(x=>x===1).length,f.filter(x=>x===2).length,f.filter(x=>x===3).length];
  const peak=Math.max(...g.map(x=>x.peak)); const ps=passPeaks(f,S);if(ps.good){stats.grids??=[];stats.grids.push(f);}stats.glue??={};const tag=JSON.stringify(ps);stats.glue[tag]=(stats.glue[tag]??0)+1;
  const key=[deg,peak,...counts,...k].join(',');
  stats.sigs[key]=(stats.sigs[key]??0)+1;
  if (stats.samples.length<2)stats.samples.push({coeff:[a,b,c,u,v,z],d:g.map(x=>x.d),key});
  return;
 }
 for(const q of lists[i]){g.push(q);dfs(i+1,total+q.d);g.pop();}
 }
 dfs(0,0);
 }
 for(let a=1;a<13;a++)for(let u=0;u<13;u++)for(let v=0;v<13;v++)for(let z=0;z<13;z++)test(a,0,0,u,v,z,"q3");
 for(let b=0;b<13;b++)for(let c=0;c<13;c++)if(b||c)for(let k=0;k<13;k++)test(0,b,c,c?k:0,0,c?0:k,"q2");
 }
 return reps;
}
function atlas8(cap, only=null) {
 const low=Array.from({length:13},(_,x)=>Array.from({length:5},(_,j)=>x**j%13));
 const app=(a,b,c,d)=>((a*13+b)*13+c)*13+d;
 const prof=Array.from({length:13**4},()=>[]);
 const dom=Array.from({length:13**3},()=>[]);
 for(let a=0;a<13;a++)for(let b=0;b<13;b++)for(let c=0;c<13;c++)for(let d=0;d<13;d++) {
  if(!a&&!b&&!c)continue; // below the least affine mass; excludes constants as well
  const base=low.map(t=>(a*t[4]+b*t[3]+c*t[2]+d*t[1])%13);
  const rows=prof[app(a,b,c,d)];
  for(let e=0;e<13;e++){
   const v=base.map(t=>(t+e)%13); if(v.some(x=>x>cap)||v.reduce((a,b)=>a+b)!==52)continue;
   rows.push({e,v,peak:Math.max(...v)});
  }
  if(rows.length)dom[app(0,a,b,c)].push(d);
 }
 const o=orbitIDs(6), good=orbitIDs(8), result=[];
 for(const m of o.reps) {
  const gid=good.to[16383^m]; if(only!==null&&gid!==only)continue;
  const s=Array.from({length:14},(_,i)=>i).filter(i=>m>>i&1);
  if(s.includes(13))throw Error("orbit normalization requires finite slopes");
  const w=s.map((x,i)=>inv[s.reduce((a,y,j)=>i===j?a:mod(a*(x-y)),1)]);
  const wz=w.map(x=>inv[x]);
  const interp=s.slice(4).map(t=>s.slice(0,4).map((x,i)=>s.slice(0,4).reduce((a,y,j)=>i===j?a:mod(a*(t-y)*inv[mod(x-y)]),1)));
  const arg=Array.from({length:169},(_,q)=>s.map(t=>mod(q%13-t*(q/13|0))));
  const own={goodOrbit:gid,s,loops:0,linearChoices:0,validSix:0,mod:0,exact:0,by:{},samples:[]};result.push(own);
  const cand=[],pick=[];
  function fixed4(a,b,c,d,e,f,tag) {
   own.loops++;
   const b4=s.map((x,i)=>mod(w[i]*a));
   const b3=s.map((x,i)=>mod(w[i]*(b+c*x)));
   const b2=s.map((x,i)=>mod(w[i]*(d+e*x+f*x*x)));
   const choices=b4.map((u,i)=>dom[app(0,u,b3[i],b2[i])]);
   if(choices.some(q=>q.length===0))return;
   for(const q0 of choices[0])for(const q1 of choices[1])for(const q2 of choices[2])for(const q3 of choices[3]){
    own.linearChoices++;
    const vals=[q0*wz[0]%13,q1*wz[1]%13,q2*wz[2]%13,q3*wz[3]%13];
    const q4= w[4]*interp[0].reduce((sum,v,i)=>sum+v*vals[i],0)%13;
    const q5= w[5]*interp[1].reduce((sum,v,i)=>sum+v*vals[i],0)%13;
    if(!prof[app(b4[4],b3[4],b2[4],q4)].length||!prof[app(b4[5],b3[5],b2[5],q5)].length)continue;
    own.validSix++;
    const first=[q0,q1,q2,q3,q4,q5];
    for(let i=0;i<6;i++) cand[i]=prof[app(b4[i],b3[i],b2[i],first[i])];
    searchConstant(0,0,tag,[b4,b3,b2,first]);
   }
  }
  function fixed2(a,b,c,d,e,f,g) {
   own.loops++;
   for(let i=0;i<6;i++){
    const x=s[i];
    cand[i]=prof[app(0,0,mod(w[i]*(a+b*x+c*x*x)),mod(w[i]*(d+e*x+f*x*x+g*x*x*x)))];
    if(!cand[i].length)return;
   }
   own.validSix++;
   searchConstant(0,0,2,[]);
  }
  function searchConstant(i,sum,tag,coeff) {
   if(i<6) {for(const q of cand[i]){pick[i]=q;searchConstant(i+1,sum+q.e,tag,coeff);}return;}
   if(mod(sum-20))return;
   own.mod++;
   const grid=[];
   for(const a of arg) {
    const total=pick[0].v[a[0]]+pick[1].v[a[1]]+pick[2].v[a[2]]+pick[3].v[a[3]]+pick[4].v[a[4]]+pick[5].v[a[5]]-20;
    if(total<0||total%13||total>13*Math.floor(6*(cap-4)/13))return;
    grid.push(total/13);
   }
   for(let t=0;t<14;t++){
    const r=Array(13).fill(0), j=s.indexOf(t);
    grid.forEach((v,q)=>r[t===13?(q/13|0):mod(q%13-t*(q/13|0))]+=v);
    if(!r.every((v,i)=>v===(j<0?4:pick[j].v[i])))throw Error("profile recovery");
   }
   own.exact++; if(passPeaks(grid,s).good){own.grids??=[];own.grids.push(grid);}own.glue??={};const tagRow=JSON.stringify(passPeaks(grid,s));own.glue[tagRow]=(own.glue[tagRow]??0)+1;
   const key=[tag,Math.max(...pick.map(q=>q.peak)),grid.filter(x=>x===2).length].join(":");
   own.by[key]=(own.by[key]??0)+1;
   if(own.samples.length<2)own.samples.push({key,coeff,dc:pick.map(t=>t.e),grid});
  }
  // After translating, eliminate the entire affine Q3 for quartics.
  for(let a=1;a<13;a++)for(let d=0;d<13;d++)for(let e=0;e<13;e++)for(let f=0;f<13;f++)fixed4(a,0,0,d,e,f,4);
  // For cubics, quotient Q2 by Q3 times all affine polynomials.
  for(let b=0;b<13;b++)for(let c=0;c<13;c++)if(b||c)for(let k=0;k<13;k++)fixed4(0,b,c,c?k:0,0,c?0:k,3);
  // For quadratics, quotient Q1 by Q2 times all affine polynomials.
  for(let a=0;a<13;a++)for(let b=0;b<13;b++)for(let c=0;c<13;c++)if(a||b||c)
   for(let j=0;j<13;j++)for(let k=0;k<13;k++)
    if(c)fixed2(a,b,c,j,k,0,0);
    else if(b)fixed2(a,b,c,j,0,0,k);
    else fixed2(a,b,c,0,0,j,k);
 }
 return result;
}const allowed10=plane10Shadow();const allow10=new Set(allowed10.profiles.map(x=>x.c));
function rank(a,n){
 const v=a.map(r=>r.slice());let z=0;
 for(let j=0;j<n;j++) {
   let k=z;while(k<v.length&&!v[k][j])k++;
   if(k===v.length)continue;
   [v[z],v[k]]=[v[k],v[z]];
   for(k=z+1;k<v.length;k++)if(v[k][j]){
    const u=mod(v[k][j]*inv[v[z][j]]);
    for(let i=j;i<n;i++)v[k][i]=mod(v[k][i]-u*v[z][i]);
   }
   z++;
 }
 return z;
}
function quartetRank(e){
 const pts=e.flatMap((v,k)=>v?[k]:[]), goods=[4,6,7,8,9,10,11,12,13];
 const mat=goods.map(d=>Array.from({length:13},(_,t)=>pts.flatMap((p,i)=>(d===13?(p/13|0):mod(p%13-d*(p/13|0)))===t?[i]:[])).flatMap(I=>{
   if(I.length!==4)throw Error("not four");
   const x=i=>d===13?pts[i]%13:(pts[i]/13|0);
   return I.slice(2).map(i=>{
    const r=Array(52).fill(0),[a,b]=I;
    r[a]=mod(x(b)-x(i));r[b]=mod(x(i)-x(a));r[i]=mod(x(a)-x(b));return r;
   });
 }));
 const out={};
 for(let n=2;n<=7;n++){
  const count={};
  for(let m=0;m<512;m++)if(pop(m)===n){
   const t=rank(mat.flatMap((x,i)=>m>>i&1?x:[]),52);
   count[t]=(count[t]??0)+1;
  }
  out[n]=count;
 }
 return out;
}
function twoCollinearGate(e) {
 const pts=e.flatMap((v,k)=>v?[k]:[]);
 const xs=pts.map(p=>p/13|0),ys=pts.map(p=>p%13);
 const good=[4,6,7,8,9,10,11,12,13];
 const mat=good.map(d=>Array.from({length:13},(_,t)=>pts.flatMap((p,i)=>(d===13?xs[i]:mod(ys[i]-d*xs[i]))===t?[i]:[])).flatMap(I=>{
  if(I.length!==4)throw Error("not four");const x=i=>d===13?ys[i]:xs[i];
  const [a,b]=I;
  return I.slice(2).map(i=>{
   const row=Array(52).fill(0);row[a]=mod(x(b)-x(i));row[b]=mod(x(i)-x(a));row[i]=mod(x(a)-x(b));return row;
  });
 }));
 const zero=[0,1,xs.findIndex(x=>x!==xs[0])];
 if(xs[0]!==xs[1]||ys[0]===ys[1]||zero[2]<0)throw Error("gauge points");
 const packed=xs.flatMap((x,i)=>xs.flatMap((y,j)=>i<j?[[i,j,mod(xs[j]-x),mod(ys[j]-ys[i])]]:[]));
 const counts={};let peak=0, low=1000, ge74=0, tests=0, viable=0, lowBad=[];
 function nullFor(rows){
  const m=rows.map(x=>x.slice()),piv=[],free=[];
  for(let j=0;j<52;j++){
   let k=piv.length;while(k<m.length&&!m[k][j])k++;
   if(k===m.length){free.push(j);continue;}
   const r=piv.length;[m[k],m[r]]=[m[r],m[k]];
   const sc=inv[m[r][j]];for(let i=j;i<52;i++)m[r][i]=mod(m[r][i]*sc);
   for(let k=0;k<m.length;k++)if(k!==r&&m[k][j]){
    const t=m[k][j];for(let i=j;i<52;i++)m[k][i]=mod(m[k][i]-t*m[r][i]);
   }
   piv.push(j);
  }
  return free.map(j=>{
   const n=Array(52).fill(0);n[j]=1;
   piv.forEach((v,i)=>n[v]=mod(-m[i][j]));return n;
  });
 }
 const gauge=zero.map(i=>{const r=Array(52).fill(0);r[i]=1;return r;});
 for(let a=0;a<9;a++)for(let b=a+1;b<9;b++){
  const M=mat[a].concat(mat[b],gauge),basis=nullFor(M);
  counts[basis.length]=(counts[basis.length]??0)+1;
  function enumerate(i,h,started){
   if(i<basis.length){
    for(let c=0;c<(started?13:2);c++)enumerate(i+1,c?h.map((v,j)=>mod(v+c*basis[i][j])):h,started||c===1);
    return;
   }
   if(!started)return;
   if(M.some(r=>mod(r.reduce((s,x,j)=>s+x*h[j],0))))throw Error("invalid kernel");
   tests++;
   const directions=new Set;
   for(const [u,v,x,y] of packed){
    const z=mod(h[v]-h[u]);
    if(!x){if(!y)throw Error("duplicate points");directions.add(169+mod(z*inv[y]));}
    else directions.add(mod(y*inv[x])*13+mod(z*inv[x]));
   }
   low=Math.min(low,directions.size);if(directions.size>=74){ge74++;return;}
   let maxPlane=0,maxLine=0;
   for(let m=0;m<183;m++){
    const [a,b,c]=m<169?[1,(m/13|0),m%13]:m<182?[0,1,m-169]:[0,0,1];
    const cc=Array(13).fill(0);for(let i=0;i<52;i++)cc[mod(a*xs[i]+b*ys[i]+c*h[i])]++;
    maxPlane=Math.max(maxPlane,...cc);
   }
   low=Math.min(low,directions.size);lowBad.push([a,b,directions.size,maxPlane]);
   if(maxPlane<=10)viable++;
  }
  enumerate(0,Array(52).fill(0),false);
 }
 return {counts,tests,ge74,low,viable,lowBad};
}
function paletteForce(e) {
 const pts=e.flatMap((v,k)=>v?[k]:[]),X=pts.map(x=>x/13|0),Y=pts.map(x=>x%13),good=[4,6,7,8,9,10,11,12,13];
 const pals=[...new Set(Array.from({length:8192},(_,m)=>m).filter(m=>pop(m)===4).map(m=>{
  const p=Array.from({length:13},(_,v)=>v).filter(v=>m>>v&1);
  let least=8192;
  for(let a=1;a<13;a++)for(let b=0;b<13;b++)least=Math.min(least,p.reduce((s,x)=>s|(1<<mod(a*x+b)),0));
  return least;
 }))].sort((a,b)=>a-b);
 if(pals.length!==7)throw Error("palette reps");
 const lines=good.map(d=>Array.from({length:13},(_,t)=>pts.flatMap((p,i)=>(d===13?X[i]:mod(Y[i]-d*X[i]))===t?[i]:[])));
 const m1=(row,j,k)=>mod(row[k]-row[j]);
 function stencil(I,ab,pm) {
  const hh=Array.from({length:13},(_,i)=>i).filter(x=>pm>>x&1);
  const vectors=[];
  for(const a of hh)for(const b of hh)for(const c of hh){
   const v=[0,ab[1]*a%13,ab[2]*b%13,ab[3]*c%13];
   let good=true;
   for(let i=1;i<4;i++)for(let j=i+1;j<4;j++)if(!(pm>>mod((v[j]-v[i])*inv[m1(ab,i,j)])&1))good=false;
   if(good)vectors.push(v);
  }
  const result=[];
  // A functional on the last two heights after subtracting the interpolating line
  // represents every possible constraint annihilating planar affine functions.
  for(let a=0;a<13;a++)for(let b=0;b<13;b++)if(a===1||a===0&&b===1){
   if(vectors.some(v=>mod(a*(v[2]-ab[2]*v[1]*inv[ab[1]])+b*(v[3]-ab[3]*v[1]*inv[ab[1]]))))continue;
   let r=Array(52).fill(0);
   const d=mod(-(a*ab[2]+b*ab[3])*inv[ab[1]]);
   r[I[1]]=d;r[I[2]]=a;r[I[3]]=b;r[I[0]]=mod(-d-a-b);
   result.push(r);
  }
  return result;
 }
 const local=lines.map((L,j)=>L.map(I=>{
  if(I.length!==4)throw Error("quartet");
  const ab=I.map(i=>good[j]===13?Y[i]:X[i]);const z=ab[0];
  for(let i=0;i<4;i++)ab[i]=mod(ab[i]-z);
  return pals.map(p=>stencil(I,ab,p));
 }));
 const rows=local.map(col=>pals.map((_,p)=>col.flatMap(L=>L[p])));
 return {good,local,rows,ranks:rows.map(ss=>ss.map(r=>rank(r,52)))};
}
function coupling(e){
 const data=paletteForce(e), pts=e.flatMap((v,k)=>v?[k]:[]), good=data.good;
 const X=pts.map(z=>z/13|0),Y=pts.map(z=>z%13);
 function plain(d){
  const pass=[];
  for(let t=0;t<13;t++) {
   const ii=pts.flatMap((_,i)=>(d===13?X[i]:mod(Y[i]-d*X[i]))===t?[i]:[]);
   const x=i=>d===13?Y[i]:X[i];
   const [a,b]=ii;
   for(const i of ii.slice(2)){
    const r=Array(52).fill(0);r[a]=mod(x(b)-x(i));r[b]=mod(x(i)-x(a));r[i]=mod(x(a)-x(b));pass.push(r);
   }
  }
  return pass;
 }
 const el=Array(52).fill(null);
 function add(rows,limit=49) {
  const changed=[];
  for(const source of rows){
   const v=Int8Array.from(source);
   for(let j=0;j<52;j++)if(v[j]){
    if(el[j]){
     const a=v[j],r=el[j];for(let k=j;k<52;k++)v[k]=(v[k]+13*13-a*r[k])%13;
    }else{
     const a=inv[v[j]];for(let k=j;k<52;k++)v[k]=v[k]*a%13;
     changed.push(j);el[j]=v;break;
    }
   }
   if(changed.length===limit)break;
  }
  return changed;
 }
 const undo=arr=>{for(const j of arr)el[j]=null;};
 const red=rows=>{const c=add(rows);const r=c.map(j=>el[j]);undo(c);return r;};
 const dups=data.rows.map(row=>[...new Map(row.map(r=>{const b=red(r);return [b.map(v=>[...v].join(",")).join("/"),b];})).values()]);
 const out={nodes:0,leaves:0,best:0,counts:{},first:[],local:dups.map(x=>x.length)};
 for(let q0=0;q0<9;q0++){
  const c=add(plain(good[q0]));
  const candidates=good.map((_,i)=>i).filter(i=>i!==q0);
  candidates.sort((a,b)=>Math.min(...dups[b].map(x=>x.length))-Math.min(...dups[a].map(x=>x.length)));
  
  const keep=[];
  const gauge=[0,1,X.findIndex(x=>x!==X[0])].map(i=>{const r=Array(52).fill(0);r[i]=1;return r;});
  function dfs(i,n,k){
   out.nodes++;
   if(n===6){
    
    out.leaves++;out.counts[k]=(out.counts[k]??0)+1;
    if(k!==48)throw Error("need general nullspace");
    const changed=add(gauge,3);if(changed.length!==3)throw Error("affine gauge");
    const h=Array(52).fill(0);let f=el.findIndex(x=>!x);
    if(f<0)throw Error("free variable");
    h[f]=1;
    for(let j=51;j>=0;j--)if(el[j])h[j]=mod(-el[j].reduce((s,a,i)=>s+a*h[i],0));
    if(el.some(r=>r&&mod(r.reduce((s,a,i)=>s+a*h[i],0))))throw Error("invalid residual");
    undo(changed);
    const key=h.join(",");
    out.heights??={};
    if(!(key in out.heights)){
     const D=new Set;let pcts=Array(14).fill(0).map(()=>new Set);
     for(let u=0;u<52;u++)for(let v=u+1;v<52;v++){
      const x=mod(X[v]-X[u]),y=mod(Y[v]-Y[u]),z=mod(h[v]-h[u]), q=x?mod(y*inv[x]):13;
      pcts[q].add(x?mod(z*inv[x]):mod(z*inv[y]));
     }
     const sum=pcts.reduce((a,b)=>a+b.size,0);
     let maxPlane=0;
     if(sum<74)for(let m=0;m<183;m++){
      const [a,b,c]=m<169?[1,(m/13|0),m%13]:m<182?[0,1,m-169]:[0,0,1];
      const cc=Array(13).fill(0);for(let j=0;j<52;j++)cc[mod(a*X[j]+b*Y[j]+c*h[j])]++;
      maxPlane=Math.max(maxPlane,...cc);
     }
     out.heights[key]={sum,maxPlane,counts:pcts.map(x=>x.size)};
    }

    if(out.first.length<3)out.first.push({q0,keep:keep.slice(),k});
    return;
   }
   if(8-i<6-n)return;
   const q=candidates[i];
   for(let j=0;j<dups[q].length;j++){
    const c=add(dups[q][j],49-k);
    if(k+c.length!==49){keep.push([q,j]);dfs(i+1,n+1,k+c.length);keep.pop();}
    undo(c);
   }
   if(8-i>6-n)dfs(i+1,n,k);
  }
  dfs(0,0,c.length);
  undo(c);
 }
 return out;
}
function insist(x,why){if(!x) throw Error(why);}
function out(x){console.log(JSON.stringify(x));}
function trim(a,exclusions){return Object.fromEntries(Object.entries(a).filter(([k])=>!exclusions.includes(k)));}
function exactSummary(with8=false) {
 const p9=run(10);
 for(const n of p9)out(trim(n,["grids","samples"]));
 out({plane10:allowed10.normalized, shadows:allowed10.profiles.map(({c,n})=>[c,n])});
 insist(p9.reduce((z,n)=>z+n.exact,0)===40,"K9 projections");
 const G=p9.flatMap(n=>n.grids??[]);
 insist(G.length===4&&G.every(z=>z.reduce((a,b)=>a+b)===52&&z.every(x=>x<=1)),"K9 10-slice gate");
 const E=G[0];
 const scalars=G.map(z=>Array.from({length:12},(_,i)=>i+1).filter(a=>z.every((v,p)=>v===E[(a*(p/13|0)%13)*13+a*(p%13)%13])));
 insist(scalars.every(s=>s.length),"grid homotheties");
 out({scalars,grid:E.flatMap((v,k)=>v?[k]:[])});
 const low7=planeShadow(7,6,4), high9=planeShadow(9,7,4);
 const r7=[];
 const maxima=[];
 for(let d=0;d<14;d++){
  let peak=0;
  for(let t=0;t<13;t++){
   const r=Array.from({length:13},(_,x)=>E[d===13?13*t+x:13*x+mod(t+d*x)]);
   const size=r.reduce((a,b)=>a+b,0);peak=Math.max(peak,size);
   if(size!==7)continue;
   r7.push([d,t,low7.canon(r)]);
  }
  maxima.push(peak);
 }
 insist(low7.normalized===3108&&high9.normalized===0,"finite planar counts");
 insist(r7.length===6&&r7.every(r=>!low7.profiles.some(x=>x.c===r[2])),"medium slice atlas");
 insist(maxima.join(",")==="10,7,10,7,4,9,4,4,4,4,4,4,4,4","slice maxima");
 out({maxima,seven:low7.normalized,sevenShadows:low7.profiles.map(x=>x.c),r7,nineStrong:high9.normalized});
 const q=quartetRank(E);
 insist(q[3][49]===84,"triple rank");out(q);
 const pair=twoCollinearGate(E);
 insist(pair.tests===218&&pair.ge74===218,"two slopes");out(pair);
 const one=coupling(E);
 insist(one.nodes===6779&&one.leaves===2374&&Object.keys(one.heights).length===9,"one slope counts");
 insist(Object.values(one.heights).every(h=>h.sum===63&&h.maxPlane===49),"no cap-10 residual");
 out({...trim(one,["heights","first","best"]),unique:Object.keys(one.heights).length,low:Object.values(one.heights)});
 for(const args of [[9,7,4],[8,7,3]])out(planarGate(...args));
 if(with8){
  const v=atlas8(10,3)[0], c={};
  for(const f of v.grids??[]){
   let max=0;
   for(let t=0;t<14;t++){
    const r=Array(13).fill(0);
    f.forEach((v,z)=>r[t===13?(z/13|0):mod(z%13-t*(z/13|0))]+=v);
    max=Math.max(max,...r);
   }
   const key=[max,Math.max(...f)];c[key]=(c[key]??0)+1;
  }
  out({goodOrbit:v.goodOrbit,s:v.s,exact:v.exact,glue:v.glue,c});
 }
}
exactSummary(typeof process!=="undefined"&&process.argv.includes("--k8"));
