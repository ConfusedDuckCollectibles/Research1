"use strict";
const processFlag=typeof process!=="undefined" && process.argv.includes("--full");

function planarMin(R, stopAt=1e8) {
const p=17,I=17,MOD=x=>(x%p+p)%p, inv=Array(p).fill(0);
for(let i=1;i<p;i++)for(let j=1;j<p;j++)if(i*j%p===1)inv[i]=j;
const dt=(a,b)=>a===b?0:a===I?1:b===I?p-1:MOD(a-b);
function masks(k){
 const seen=new Set(), M=[];let scanned=0, covered=0;
 function canon(s){
  scanned++;const base=s.reduce((z,t)=>z|1<<t,0);if(seen.has(base))return;
  const orbit=new Set();
  for(const a of s)for(const b of s)for(const c of s){if(a===b||a===c||b===c)continue; const z=dt(b,c)*inv[dt(b,a)]%p, mask=s.reduce((r,d)=>r|(1<<(d===c?I:dt(d,a)*inv[dt(d,c)]*z%p)),0);orbit.add(mask);}
  for(const r of orbit){if(seen.has(r))throw Error("overlap");seen.add(r);}
  covered+=orbit.size;M.push({s:s.slice(),size:orbit.size});
 }
 const s=[0,1,I];
 function go(j){if(s.length===k){canon(s);return;}for(let i=j;i<=p-(k-s.length);i++){s.push(i);go(i+1);s.pop();}}
 go(2);if(covered!==scanned)throw Error("orbit coverage");
 return {M,scanned,covered};
}
const data=masks(R),out=[];
for(const {s,size} of data.M){
 const sMask=s.reduce((a,b)=>a|1<<b,0);
 const forbid=s.filter(a=>a!==I);
 const allow=Array.from({length:p},(_,i)=>i).filter(v=>!(sMask>>v&1));
 const F=new Int8Array(17).fill(-1);
 const Tr=Array.from({length:17},(_,u)=>Array.from({length:17},(_,x)=>x*u%p));
 const tape=allow.map(()=>new Int8Array(p));
 const changed=Array.from({length:17},()=>new Int32Array(p));
 const fbit=Array.from({length:17},(_,d)=>Array.from({length:17},(_,val)=>forbid.reduce((z,s)=>z|1<<(val+Tr[s][d])%17,0)));
 const dom=new Int32Array(p).fill((1<<p)-1);
 let n=0,max=0,wit=null, timeout=false;
 function dfs(unseen,sz){
  n++;if(n>=stopAt){timeout=true;return true;}
  if(sz>max)max=sz;
  if(unseen===0){wit=Array.from(F);return true;}
  let mi=17,sel=-1,list=0;
  for(let m=unseen;m;m&=m-1){const b=m&-m,x=31-Math.clz32(b),d=dom[x];
    let cnt=0;for(let a=d;a;a&=a-1)cnt++;
    if(cnt<mi){mi=cnt;sel=x;list=d;if(mi<=1)break;}}
  if(!mi)return false;
  unseen&=~(1<<sel);
  for(let i=list;i;i&=i-1){
   const y=31-Math.clz32(i&-i);let good=true;
   for(let j=0;j<allow.length;j++)if(tape[j][MOD(y-Tr[allow[j]][sel])]===5){good=false;break;}
   if(!good)continue;
   F[sel]=y;
   for(let j=0;j<allow.length;j++)tape[j][MOD(y-Tr[allow[j]][sel])]++;
   const ch=changed[sz];
   let bad=false;
   for(let t=unseen;t;t&=t-1){const x=31-Math.clz32(t&-t);ch[x]=dom[x];dom[x]&=~fbit[(x-sel+17)%17][y];if(!dom[x])bad=true;}
   if(!bad&&dfs(unseen,sz+1))return true;
   for(let t=unseen;t;t&=t-1){const x=31-Math.clz32(t&-t);dom[x]=ch[x];}
   for(let j=0;j<allow.length;j++)tape[j][MOD(y-Tr[allow[j]][sel])]--;
   F[sel]=-1;
  }
  return false;
 }
 dom[0]=1;dfs((1<<p)-1,0);
 out.push({s,size,n,max,wit,timeout});
}
return {R,data,out};
}

function atlas(n=12, missingN=6, cap=5, limit=2e8){
const p=17, mod=x=>(x%p+p)%p, U=(1<<p)-1;
function domains(k){
 const seen=new Set(), reps=[], s=[];
 function root(){
  const b=s.reduce((z,x)=>z|1<<x,0);if(seen.has(b))return;
  const a=new Set();
  for(let u=1;u<17;u++)for(let v=0;v<17;v++)a.add(s.reduce((z,x)=>z|1<<((u*x+v)%17),0));
  for(const b of a){if(seen.has(b))throw Error("bad orbit");seen.add(b);}
  reps.push(s.slice());
 }
 function go(i){if(s.length===k){root();return;}for(let j=i;j<=17-(k-s.length);j++){s.push(j);go(j+1);s.pop();}}
 go(0);return {n:seen.size,reps};
}
const X=domains(n), M=domains(missingN);
const nodes=Array(X.reps.length).fill(0),all=[];
const F=new Int8Array(17).fill(-1);
const Tr=Array.from({length:17},(_,u)=>Array.from({length:17},(_,x)=>x*u%p));
let N=0;
for(let j=0;j<M.reps.length;j++){
 const miss=M.reps[j],m=miss.reduce((z,x)=>z|1<<x,0), allow=Array.from({length:17},(_,i)=>i).filter(x=>!(m>>x&1));
 const tape=allow.map(()=>new Int8Array(p));
 const changed=Array.from({length:17},()=>new Int32Array(p));
 const fbit=Array.from({length:17},(_,d)=>Array.from({length:17},(_,val)=>miss.reduce((z,s)=>z|1<<(val+Tr[s][d])%17,0)));
 const dom=new Int32Array(p);
 for(let i=0;i<X.reps.length;i++){
  const x=X.reps[i],xm=x.reduce((z,v)=>z|1<<v,0);dom.fill(U);dom[x[0]]=1;
  let cnt=0, vals=[];
  function dfs(unseen,sz){
   if(++N>limit)throw Error("limit "+N);
   nodes[i]++;
   if(!unseen){cnt++;vals.push(x.map(v=>F[v]));return;}
   let mi=17,sel=-1,list=0;
   for(let m=unseen;m;m&=m-1){const b=m&-m,x=31-Math.clz32(b),d=dom[x];
    let cnt=0;for(let a=d;a;a&=a-1)cnt++;
    if(cnt<mi){mi=cnt;sel=x;list=d;if(mi<=1)break;}}
   if(!mi)return;
   unseen&=~(1<<sel);
   for(let i=list;i;i&=i-1){
    const y=31-Math.clz32(i&-i);let good=true;
    for(let j=0;j<allow.length;j++)if(tape[j][mod(y-Tr[allow[j]][sel])]===cap){good=false;break;}
    if(!good)continue;
    F[sel]=y;
    for(let j=0;j<allow.length;j++)tape[j][mod(y-Tr[allow[j]][sel])]++;
    const ch=changed[sz];
    let bad=false;
    for(let t=unseen;t;t&=t-1){const x=31-Math.clz32(t&-t);ch[x]=dom[x];dom[x]&=~fbit[(x-sel+17)%17][y];if(!dom[x])bad=true;}
    if(!bad)dfs(unseen,sz+1);
    for(let t=unseen;t;t&=t-1){const x=31-Math.clz32(t&-t);dom[x]=ch[x];}
    for(let j=0;j<allow.length;j++)tape[j][mod(y-Tr[allow[j]][sel])]--;
    F[sel]=-1;
   }
  }
  dfs(xm,0);
  if(cnt){
   let fixed=[];
   for(let a=0;a<n;a++)for(let b=a+1;b<n;b++){const v=mod(vals[0][a]-vals[0][b]);if(vals.every(q=>mod(q[a]-q[b])===v))fixed.push([a,b,v]);}
   all.push({i,j,cnt,fixed,vals});
  }
 }
}
return {n,N,nX:X.n,nM:M.n,X:X.reps,M:M.reps,nodes,all};
}

function assert(ok,msg){if(!ok)throw Error(msg||"failed exact certificate");}
function fRank(A,p=17){
 A=A.map(row=>row.slice());const n=A[0]?.length??0;
 let pos=0;
 for(let j=0;j<n && pos<A.length;j++){
  const k=A.findIndex((r,i)=>i>=pos&&r[j]);if(k<0)continue;
  [A[pos],A[k]]=[A[k],A[pos]];
  for(let i=pos+1;i<A.length;i++){
   const a=A[pos][j],b=A[i][j];
   if(b)for(let k=j;k<n;k++) A[i][k]=(a*A[i][k]-b*A[pos][k]+17*17)%17;
  }
  pos++;
 }
 return pos;
}
function checkGraph(x,y,missing,cap){
 const p=17,mod=x=>(x%p+p)%p;
 assert(x.length===y.length&&new Set(x).size===x.length&&y.every(v=>v>=0&&v<p),"graph labels");
 for(let m=0;m<p;m++){
  const z=Array(p).fill(0);
  for(let i=0;i<x.length;i++)z[mod(y[i]-m*x[i])]++;
  assert(Math.max(...z)<=(missing.includes(m)?1:cap),"graph witness invalid");
 }
}
function quadraticGate(){
 let Q=Infinity,W=null,cases=0;
 for(let a=1;a<17;a++)for(let c=0;c<17;c++){
  let s=0;for(let x=0;x<17;x++)s+=(a*x*x+c)%17;cases++;
  if(s<Q){Q=s;W=[a,c];}
 }
 assert(Q===102&&cases===272);
 return {affine:136,quadratic:Q,witness:W,cases,lpSizeCap:135,lpNumerator:289*84,lpDenominator:132};
}
const D=planarMin(5,6000000);
assert(D.data.scanned===105&&D.out.length===4&&D.out.every(d=>!d.timeout&&!d.wit),"full layer fails");
const C=planarMin(4,6000000);
assert(C.out.some(d=>d.wit),"unexpected absence of positive 4-missing control");
for(const u of C.out)if(u.wit)checkGraph(Array.from({length:17},(_,i)=>i),u.wit,u.s,5);
const A=atlas();
const cr={},ranks={},pairs={};
let nShape=0,maxCell=0;
for(const u of A.all){
 assert(new Set(u.vals.map(v=>v.join(","))).size===u.cnt);
 for(const v of u.vals){assert(v[0]===0);checkGraph(A.X[u.i],v,A.M[u.j],5);}
 const span=fRank(u.vals.slice(1).map(v=>v.map((x,i)=>(x-u.vals[0][i]+17)%17)));
 cr[u.cnt]=(cr[u.cnt]??0)+1;ranks[span]=(ranks[span]??0)+1;
 const k=u.cnt+"_"+span;pairs[k]=(pairs[k]??0)+1;
 maxCell=Math.max(maxCell,u.cnt);
 nShape+=u.cnt;
 assert(u.cnt<=4&&span<=2,"palette shape or affine rank");
}
assert(A.nX===6188&&A.X.length===25&&A.nM===12376&&A.M.length===49&&A.N===584348&&A.all.length===207&&nShape===353);
const out={quad:quadraticGate(),full:D.out.map(u=>({s:u.s,n:u.n,max:u.max})),control:C.out.filter(u=>u.wit).map(u=>({s:u.s,wit:u.wit})),atlas:{nX:A.nX,nM:A.nM,rX:A.X.length,rM:A.M.length,nodes:A.N,nonempty:A.all.length,nShape,maxCell,cr,ranks,pairs},...(processFlag?{tables:A.all,X:A.X,M:A.M}:{})};
console.log(JSON.stringify(out));
