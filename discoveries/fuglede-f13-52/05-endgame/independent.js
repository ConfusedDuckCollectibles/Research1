"use strict";
// Independent finite controls for the size-52 low-pencil and AP endgame.
const F=Array.from({length:13},(_,i)=>i);
const mod=x=>(x%13+13)%13;
const inv=x=>F.find(y=>mod(x*y)===1);
function insist(ok,n){ if(!ok) throw Error(n); }
const sample=[0,1,2,3,4], terms=F.slice(5).map(t=>sample.map((x,i)=>{
  let f=1;
  for(let j=0;j<5;j++) if(j!==i) f=mod(f*(t-sample[j])*inv(x-sample[j]));
  return f;
}));
let visits=0,cap=0,mass=0; const wcounts={};const a=[];
function profile() {
  if(a.length<5) {
    for(let y=0;y<=10;y++){a.push(y);profile();a.pop();}
    return;
  }
  visits++;
  const row=a.concat(terms.map(w=>mod(w.reduce((s,x,i)=>s+x*a[i],0))));
  if(row.some(v=>v>10))return;
  cap++;
  if(row.reduce((s,x)=>s+x,0)!==52)return;
  mass++;
  for(let m=0;m<=7;m++)insist(mod(row.reduce((s,x,i)=>s+x*i**m,0))===0,"moment");
  const v=row.reduce((s,x)=>s+x*(x-1)/2,0), k=(v-78)/13;
  insist(Number.isInteger(k) && (k>=2 || (k===0&&row.every(n=>n===4))),"weight");
  wcounts[k]=(wcounts[k]||0)+1;
}
profile();
const table=[];
for(let m=66;m<=73;m++){
  const total=m*(m-1)/2;
  const lax=Math.floor((total-183*6)/22);
  const sharper=Math.floor((total-56*m+1830)/6);
  const eight=Math.min(lax,sharper);
  const forced=2*(110-eight), budget=18*(73-m);
  insist(forced>budget,"budget");
  table.push({m,eight,forced,budget});
}
const fsets=[];
for(let a=0;a<13;a++)for(let b=a+1;b<13;b++)for(let c=b+1;c<13;c++)for(let d=c+1;d<13;d++) fsets.push([a,b,c,d]);
const ap=new Map;
for(const z of F)for(let step=1;step<=6;step++){
  const X=[z,mod(z+step),mod(z+2*step),mod(z+3*step)].sort((x,y)=>x-y);
  insist(!ap.has(X.join(',')),"step unique");
  ap.set(X.join(','),step);
  const cnt=Array(13).fill(0);
  for(const x of X)for(const y of X)if(x!==y)cnt[mod(x-y)]++;
  insist(cnt.every((c,i)=>(c===3)===(i===step||i===13-step)),"difference multiplicity");
  const mean=mod(X.reduce((s,x)=>s+x,0)*inv(4));
  insist(mod(X.reduce((s,x)=>s+(x-mean)**2,0))===mod(5*step**2),"AP variance");
}
const pairCounts={};let tested=0,apEquality=0,apControl=0;
for(let i=0;i<fsets.length;i++)for(let j=0;j<=i;j++){
  const X=fsets[i],Y=fsets[j];let s=0, n=0;
  for(const x of X)for(const y of Y){
    const r=1<<mod(x-y);
    if(!(s&r)){s|=r;n++;}
  }
  tested++; pairCounts[n]=(pairCounts[n]||0)+1;
  const sx=ap.get(X.join(',')), sy=ap.get(Y.join(','));
  const same=sx!==undefined&&sx===sy;
  if(same)apControl++;
  insist(n>=7 && ((n===7)===same),"AP4 equality");
  if(n===7)apEquality++;
}
console.log(JSON.stringify({visits,cap,mass,wcounts,table}));
console.log(JSON.stringify({sets:fsets.length,tested,APsets:ap.size,pairCounts,apEquality,apControl,
  quotientDirections:Math.ceil((13*12/2)/6),minD:1+13*7}));
