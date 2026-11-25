"use strict";
const p=13, m=x=>(x%p+p)%p,inv=[0,1,7,9,10,8,11,2,5,3,4,6,12];
const all=Array.from({length:169},(_,i)=>[Math.floor(i/p),i%p]);
const check=(b,l)=>{if(!b)throw Error(l);};
const tr=[];for(let a=0;a<p;a++)for(let b=a+1;b<p;b++)for(let c=b+1;c<p;c++)tr.push([a,b,c]);
function step(A){
  const k=A.find(c=>A.filter(x=>x!==c).reduce((s,x)=>s+x,0)%p===(2*c)%p);
  if(k===undefined)return -1;
  const d=m(A.find(a=>a!==k)-k);return Math.min(d,p-d);
}
const sg=tr.map(step);
let tight=0;
for(let i=0;i<tr.length;i++)for(let j=0;j<tr.length;j++){
 const d=new Set(tr[i].flatMap(a=>tr[j].map(b=>m(a-b)))).size;
 check(d>=5,"too small triple differences");
 if(d===5){tight++;check(sg[i]!==-1&&sg[i]===sg[j],"no common AP step");}
}
const quintet=new Set([11,12,0,1,2]);
const stab=[];
for(let a=1;a<p;a++)for(let b=0;b<p;b++)if([...quintet].every(u=>quintet.has(m(a*u+b))))stab.push([a,b]);
check(JSON.stringify(stab)==="[[1,0],[12,0]]","extra affine quintet stabilizer");
// Center nonconstant p-polynomials by translating variable x. Compute integer mass.
const minima=[];
for(let k=1;k<=3;k++){
 let least=1000;const at39=[];
 for(let a=1;a<p;a++)for(let b=0;b<(k===3?p:1);b++)for(let c=0;c<p;c++){
   const val=Array.from({length:p},(_,t)=>m(a*(t**k)+b*t+c));
   const mass=val.reduce((x,y)=>x+y,0);
   least=Math.min(least,mass);
   if(mass===39)at39.push({a,b,c});
 }
 minima.push({degree:k,least,at39});
}
check(minima.map(c=>c.least).join()==="78,52,39","bad residue threshold");
check(minima[2].at39.length===8&&minima[2].at39.every(c=>c.b===0&&[3,9].includes(c.c)),"bad cubic profiles");

function dir(x,y){return x===0?13:m(y*inv[x]);}
const disp=all.map(([x,y])=>all.map(([a,b])=>{
 let dx=m(a-x), dy=m(b-y);
 if(dx===0&&dy===0)return {key:-1,d:-1,sign:0};
 const d=dir(dx,dy);
 let sign=1;if((dx||dy)>6){sign=-1;dx=m(-dx);dy=m(-dy);}
 return {key:13*dx+dy,d,sign};
}));
const dirs=[];for(let j=0;j<p;j++)for(let k=0;k<p;k++)dirs.push([1,j,k]);
for(let j=0;j<p;j++)dirs.push([0,1,j]);dirs.push([0,0,1]);
const idx=new Map(dirs.map((v,i)=>[v.join(),i]));
const per=dirs.map(([x,y,z])=>dirs.flatMap(([a,b,c],j)=>m(a*x+b*y+c*z)===0?[j]:[]));
function proj(x,y,z){
 const f=inv[x||y||z];
 return idx.get([m(x*f),m(y*f),m(z*f)].join());
}
const heightsHisto={}, dims={};let ncaps=0, nheights=0;
function lift(cap){
 const v=Array.from({length:14},()=>[]);
 const rows=[];
 for(let i=0;i<cap.length;i++)for(let j=i+1;j<cap.length;j++){
  const {d,sign}=disp[cap[i]][cap[j]];
  const coeff=Array.from({length:4},(_,k)=>sign*(+(k+3===j)-+(k+3===i)));
  if(v[d].length)rows.push(coeff.map((a,k)=>m(a-v[d][0][k])));
  v[d].push(coeff);
 }
 // Exact row reduction; derive all four-coordinate heights in the nullspace.
 let next=0;const piv=[];
 for(let i=0;i<4;i++){
  const k=rows.findIndex((u,j)=>j>=next&&u[i]);
  if(k<0)continue;
  [rows[k],rows[next]]=[rows[next],rows[k]];
  rows[next]=rows[next].map(a=>m(a*inv[rows[next][i]]));
  for(let j=0;j<rows.length;j++)if(j!==next&&rows[j][i]){
   const c=rows[j][i];rows[j]=rows[j].map((v,k)=>m(v-c*rows[next][k]));
  }
  piv.push(i);next++;
 }
 const non=[0,1,2,3].filter(i=>!piv.includes(i)), t=Array(4).fill(0);
 dims[non.length]=(dims[non.length]||0)+1;
 function finish(k){
  if(k<non.length){for(let x=0;x<p;x++){t[non[k]]=x;finish(k+1)}return;}
  for(let j=0;j<piv.length;j++)t[piv[j]]=m(-non.reduce((s,i)=>s+rows[j][i]*t[i],0));
  const h=[0,0,0,...t], set=new Set([182]); // last PG point is vertical
  for(let a=0;a<7;a++)for(let b=a+1;b<7;b++){
   const [x,y]=all[cap[a]], [u,v]=all[cap[b]];
   for(let u0=-2;u0<=2;u0++)set.add(proj(m(u-x),m(v-y),m(h[b]-h[a]+u0)));
  }
  const count=Array(183).fill(0);
  for(const d of set)for(const k of per[d])count[k]++;
  const largest=Math.max(...count);
  heightsHisto[largest]=(heightsHisto[largest]||0)+1;nheights++;
 }
 finish(0);
}
const cap=[0,13,1], used=Array(14).fill(-1);
for(let i=0;i<3;i++)for(let j=i+1;j<3;j++){
 const {d,key}=disp[cap[i]][cap[j]];used[d]=key;
}
const vis=Array(8).fill(0);
function enumerate(start){
 vis[cap.length]++;
 if(cap.length===7){ncaps++;lift(cap);return;}
 for(let x=start;x<169;x++){
  if(x===cap[0]||x===cap[1]||x===cap[2])continue;
  let seen=0,good=true;const change=[];
  for(const u of cap){
   const {d,key}=disp[u][x], bit=1<<d;
   if(seen&bit || (used[d]!==-1&&used[d]!==key)){good=false;break;}
   seen|=bit;
   if(used[d]===-1)change.push([d,key]);
  }
  if(!good)continue;
  for(const [d,key] of change)used[d]=key;
  cap.push(x);enumerate(x+1);cap.pop();
  for(const [d,key] of change)used[d]=-1;
 }
}
enumerate(0);
check(ncaps===210&&nheights===2730&&heightsHisto[12]===2730,"skeleton finite gate");
console.log(JSON.stringify({tightTriplePairs:tight,quintetAffineStabilizer:stab,minima,skeleton:{caps:ncaps,dims,vis,normalizedHeights:nheights,heightsHisto}},null,2));