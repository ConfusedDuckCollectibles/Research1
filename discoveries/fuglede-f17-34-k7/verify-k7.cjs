"use strict";
function residues(p,dmax){
 let pow=Array.from({length:dmax+1},(_,d)=>Array.from({length:p},(_,x)=>{let s=1;for(let i=0;i<d;i++)s=s*x%p;return s;}));
 let out=[];
 for(let d=1;d<=dmax;d++){
  const hist=new Int32Array(p), C=new Int32Array(d+1), H=new Int32Array(p), v=new Int32Array(p);
  let min=p*p, minw=null,n2=0,n3=0,runs=0, candidates=[];
  if(d===1){
   out.push({d,min:p*(p-1)/2,n2:0,n3:0,runs:p*(p-1),minw:[0,1],candidates:[]});continue;
  }
  function calc(k){
   if(k===d-1){
    hist.fill(0); let z=0;
    for(let x=0;x<p;x++){const r=v[x]%p; hist[r]++;z+=r;}
    runs++;
    for(let c=0;c<p;c++){
     if(z<min){min=z;minw=Array.from(C);minw[0]=c;}
     if(z<=2*p){n2++; if(candidates.length<100000){const a=Array.from(C);a[0]=c;candidates.push(a);}}
     if(z<=3*p)n3++;
     z+=p-p*hist[p-1-c];
    }
    return;
   }
   const X=pow[k]; for(let i=0;i<p;i++){C[k]=i;calc(k+1);for(let x=0;x<p;x++)v[x]+=X[x];}
   for(let x=0;x<p;x++)v[x]-=p*X[x];
  }
  for(let a=1;a<p;a++){ C[d]=a;for(let x=0;x<p;x++)v[x]=a*pow[d][x];calc(1);}
  out.push({d,min,n2,n3,runs,minw,candidates});
 }
 return out;
}

function septic(p,dmax){
 let pow=Array.from({length:dmax+1},(_,d)=>Array.from({length:p},(_,x)=>{let s=1;for(let i=0;i<d;i++)s=s*x%p;return s;}));
 let out=[];
 for(let d=7;d<=dmax;d++){
  const hist=new Int32Array(p), C=new Int32Array(d+1), H=new Int32Array(p), v=new Int32Array(p);
  let min=p*p, minw=null,n2=0,n3=0,runs=0, candidates=[];
  if(d===1){
   out.push({d,min:p*(p-1)/2,n2:0,n3:0,runs:p*(p-1),minw:[0,1],candidates:[]});continue;
  }
  function calc(k){
   if(k===d-1){
    hist.fill(0); let z=0;
    for(let x=0;x<p;x++){const r=v[x]%p; hist[r]++;z+=r;}
    runs++;
    for(let c=0;c<p;c++){
     if(z<min){min=z;minw=Array.from(C);minw[0]=c;}
     if(z<=2*p){n2++; if(candidates.length<100000){const a=Array.from(C);a[0]=c;candidates.push(a);}}
     if(z<=3*p)n3++;
     z+=p-p*hist[p-1-c];
    }
    return;
   }
   const X=pow[k]; for(let i=0;i<p;i++){C[k]=i;calc(k+1);for(let x=0;x<p;x++)v[x]+=X[x];}
   for(let x=0;x<p;x++)v[x]-=p*X[x];
  }
  for(let a=1;a<=1;a++){ C[d]=a;for(let x=0;x<p;x++)v[x]=a*pow[d][x];calc(1);}
  out.push({d,min,n2,n3,runs,minw,candidates});
 }
 return out;
}


function slopeOrbits(L){
 const P=17, mod=x=>(x%P+P)%P, inv=Array(17).fill(0);
 for(let x=1;x<P;x++)for(let y=1;y<P;y++)if(x*y%P===1)inv[x]=y;
 const deterr=(x,a)=>x===17?1:a===17?16:mod(x-a);
 const reps=[],seen=new Set();let normalized=0,sum=0;
 const s=[0,1];
 function root(){
  normalized++;
  const sm=s.reduce((b,x)=>b|(1<<x),0);if(seen.has(sm))return;
  const orbit=new Set();let g=Array.from({length:18},(_,i)=>i).filter(i=>!(sm>>i&1));
  for(const a of s)for(const b of s)if(a!==b)for(const c of g){
   const sc=deterr(b,c)*inv[deterr(b,a)]%17;
   let m=0;for(const x of s){const y=deterr(x,a)*inv[deterr(x,c)]*sc%17;m|=1<<y;}
   if((m&3)!==3)throw Error("normalize");
   orbit.add(m);
  }
  for(const q of orbit)if(seen.has(q))throw Error("orbit collision");else seen.add(q);
  reps.push(s.slice());sum+=orbit.size;
 }
 function dfs(i){if(s.length===L){root();return;}for(let j=i;j<=17-(L-s.length);j++){s.push(j);dfs(j+1);s.pop();}}
 dfs(2);if(sum!==normalized||sum!==seen.size)throw Error("incomplete");
 return {reps,normalized,sum};
}


function cap7gate(reps){
const P=17,R=x=>(x%17+17)%17,iv=new Int8Array(17);
for(let i=1;i<17;i++)for(let j=1;j<17;j++)if(i*j%17===1)iv[i]=j;
const sl=(a,b)=>a[0]===b[0]?17:R((a[1]-b[1])*iv[R(a[0]-b[0])]);
const out=[];
for(const S of reps){
 const good=Array.from({length:18},(_,u)=>!S.includes(u));
 const E=[[0,0],[0,1]],X=[];
 for(let x=1;x<17;x++)for(let y=0;y<17;y++)if(E.every(t=>good[sl([x,y],t)]))X.push([x,y]);
 let n=0,m=2,num=0;const keep=[];
 function go(j){
  n++;m=Math.max(m,E.length);
  if(E.length===7){num++;if(keep.length===0)keep.push(E.slice());return;}
  for(let i=j;i<X.length;i++){
   let z=0,ok=true;
   for(const b of E){const s=sl(X[i],b);if(!good[s]||z>>s&1){ok=false;break;}z|=1<<s;}
   if(!ok)continue;
   E.push(X[i]);go(i+1);E.pop();
  }
 }
 go(0);out.push({S,comp:good.map((u,i)=>u?i:-1).filter(i=>i>=0),n,m,num,witness:keep[0]});
}
return out;
}

function low6(profs,reps){
const p=17,R=x=>(x%17+17)%17,iv=new Int8Array(p),pw=Array.from({length:p},(_,x)=>Array.from({length:11},(_,d)=>{let z=1;for(let i=0;i<d;i++)z=z*x%17;return z;}));
for(let a=1;a<p;a++)for(let b=1;b<p;b++)if(a*b%17===1)iv[a]=b;
const zero=x=>x%17===0;
const bc=Array.from({length:7},(_,d)=>Array.from({length:d+1},(_,j)=>{let r=1;for(let i=1;i<=j;i++)r=r*(d-i+1)*iv[i]%17;return r;}));
const val=(c,x)=>c.reduceRight((v,u)=>(v*x+u)%17,0);
const B=Object.fromEntries(profs.filter(c=>Math.max(...pw.map((_,x)=>val(c,x)))<=7).map(c=>[c[6],c]));
const shift=Object.fromEntries(Object.entries(B).map(([a,c])=>[a,Array.from({length:p},(_,h)=>c.map((_,j)=>c.reduce((v,u,d)=>d<j?v:(v+u*bc[d][j]*pw[h][d-j])%17,0)))]));
const out=[];
for(const s of reps){
 const w=s.map((a,i)=>iv[s.reduce((v,b,j)=>i===j?v:R(v*(a-b)),1)]),win=w.map(x=>iv[x]),power=s.map(x=>pw[x]);
 const cs={s,lead:0,pre:0,sq:0,quartic:0,const:0,formal:0,actual:0};
 const aGroups=[1*win[0]%17,3*win[0]%17];
 const interp=(Z)=>Z.reduce((poly,x,i)=>{
  let t=[1],d=1;for(let j=0;j<Z.length;j++)if(j!==i){const y=Z[j],u=Array(t.length+1).fill(0);for(let k=0;k<t.length;k++){u[k]=R(u[k]-y*t[k]);u[k+1]=R(u[k+1]+t[k]);}t=u;d=R(d*(x-y));}
  poly.push(t.map(x=>x*iv[d]%17));return poly;},[]);
 const L9=interp(s.slice(0,9));
 for(const c0 of aGroups)for(let c1=0;c1<17;c1++)for(let c2=0;c2<17;c2++)for(let c3=0;c3<17;c3++){
  const Q=[c0,c1,c2,c3],qd=c3?3:c2?2:c1?1:0;
  const q=s.map(x=>val(Q,x)),az=q.map((x,i)=>x*w[i]%p);
  if(az.some(u=>!B[u]))continue;
  cs.lead++;
  const A=az.map((u,i)=>B[u][4]*win[i]*q[i]%17);
  let ap=Array(9).fill(0);
  for(let i=0;i<9;i++)for(let j=0;j<9;j++)ap[j]=(ap[j]+A[i]*L9[i][j])%17;
  if(val(ap,s[9])!==A[9]||val(ap,s[10])!==A[10])continue;
  cs.pre++;
  function rem(P){
   let a=P.slice();
   for(let j=a.length-1;j>=qd;j--){const v=a[j]*iv[Q[qd]]%17;for(let h=0;h<qd;h++)a[j-qd+h]=R(a[j-qd+h]-v*Q[h]);}
   return a.slice(0,qd);
  }
  const goal=rem(ap),basis1=Array.from({length:3},(_,k)=>rem([...Array(k+1).fill(0),16,1]));
  const sigma=s.map((x,i)=>modx(x*R(x-1)*iv[6*q[i]%17]));
  function modx(x){return x%17;}
  const g=[];
  for(let h0=0;h0<17;h0++)for(let h1=0;h1<17;h1++)for(let h2=0;h2<17;h2++){
   const a=Array.from({length:qd},(_,j)=>(h0*basis1[0][j]+h1*basis1[1][j]+h2*basis1[2][j])%17);
   const square=Array(Math.max(0,2*qd-1)).fill(0);for(let i=0;i<qd;i++)for(let j=0;j<qd;j++)square[i+j]=(square[i+j]+a[i]*a[j])%17;
   const rt=rem(square);if(rt.some((v,i)=>v!==goal[i]))continue;
   cs.sq++;
   let mass=0;for(let i=0;i<11;i++){const x=s[i],h=(h0+h1*x+h2*power[i][2])*sigma[i]%17;g[i]=shift[az[i]][h];mass+=g[i][0];}
   function check(d){
    for(let j=0;j<=d;j++){let v=0;for(let i=0;i<11;i++)v+=power[i][j]*g[i][d];if(v%17!==0)return false;}
    return true;
   }
   if(!check(4))continue;
   cs.quartic++;
   if(mass%17!==3)continue;cs.const++;
   if(!check(3)||!check(2)||!check(1))continue;cs.formal++;
   let ok=true;
   for(let x=0;x<17&&ok;x++)for(let y=0;y<17;y++){let v=0;for(let i=0;i<11;i++)v+=val(g[i],R(y-s[i]*x));if(![20,37,54].includes(v)){ok=false;break;}}
   if(ok){cs.actual++;if(!cs.first)cs.first=g.slice();}
  }
 }
 out.push(cs);
}
return {eligible:Object.keys(B).length,out};
}

function low7(L6,L7,reps){
const P=17,R=x=>(x%P+P)%P,iv=Array(17).fill(0),power=Array.from({length:17},(_,x)=>Array.from({length:10},(_,d)=>{let z=1;for(let i=0;i<d;i++)z=z*x%17;return z;}));
for(let a=1;a<P;a++)for(let b=1;b<P;b++)if(a*b%17===1)iv[a]=b;
const v=(c,x)=>c.reduceRight((v,u)=>(v*x+u)%17,0);
const cap=c=>Math.max(...power.map((_,i)=>v(c,i)))<=7;
const B6=Object.fromEntries(L6.filter(cap).map(c=>[c[6],c]));
const B7=Array.from({length:17},()=>[]);
for(let u=1;u<17;u++)for(const c of L7.filter(cap)){const f=c.map((k,j)=>k*power[u][j]%17);B7[f[7]].push(f);}
const bc=Array.from({length:8},(_,d)=>Array.from({length:d+1},(_,j)=>{let r=1;for(let i=1;i<=j;i++)r=r*(d-i+1)*iv[i]%17;return r;}));
function shifted(c,h){return c.map((_,j)=>c.reduce((v,u,d)=>d<j?v:(v+u*bc[d][j]*power[h][d-j])%17,0));}
const G6=Array.from({length:17},(_,a)=> B6[a]?power.map((_,h)=>shifted(B6[a],h)):[]);
const G7=B7.map(f=>f.map(c=>power.map((_,h)=>shifted(c,h))));
const out=[];
for(const s of reps){
 const w=s.map((a,i)=>iv[s.reduce((v,b,j)=>i===j?v:R(v*(a-b)),1)]),win=w.map(a=>iv[a]),pow=s.map(t=>power[t]);
 const rec={s,lead:0,base:0,poly:0,square:0,c5:0,c4:0,const:0,form:0,actual:0};
 const cache=new Map();
 function lag(ix){
  const key=ix.join(",");if(cache.has(key))return cache.get(key);
  const a=ix.map(i=>s[i]),co=a.map((x,i)=>{
   let c=[1],d=1;
   for(let j=0;j<a.length;j++)if(j!==i){const y=a[j],b=Array(c.length+1).fill(0);for(let k=0;k<c.length;k++){b[k]=R(b[k]-y*c[k]);b[k+1]=R(b[k+1]+c[k]);}c=b;d=R(d*(x-y));}
   return c.map(x=>x*iv[d]%17);
  });
  const tab=s.map(x=>co.map(c=>v(c,x)));
  const r={co,tab};cache.set(key,r);return r;
 }
 const D=[];
 function Qrun(Q){
  rec.lead++;
  const q=s.map(x=>v(Q,x)),a=q.map((z,i)=>z*w[i]%17);
  const M=Array.from({length:11},(_,i)=>i).filter(i=>q[i]),Z=Array.from({length:11},(_,i)=>i).filter(i=>!q[i]);
  const start=M.slice(0,7),rest=M.slice(7),L=lag(start);
  const values=[],choice=[];
  const list=M.map(i=>B7[a[i]].map(c=>R(c[5]*win[i]*q[i])));
  const mask=list.map(g=>g.reduce((v,x)=>v|1<<x,0));
  const gauge=[s[M[0]],s[M[1]]],A3=[R(gauge[0]*gauge[1]),R(-gauge[0]-gauge[1]),1];
  const need=Array.from({length:11},(_,i)=>(s[i]-gauge[0])*(s[i]-gauge[1]));
  const d=Q[2]?2:Q[1]?1:0;
  function div(A){
   const P=A.slice(),H=Array(Math.max(0,A.length-d)).fill(0);
   for(let i=P.length-1;i>=d;i--){const h=P[i]*iv[Q[d]]%17;H[i-d]=h;for(let j=0;j<=d;j++)P[i-d+j]=R(P[i-d+j]-h*Q[j]);}
   return {H,rem:P.slice(0,d)};
  }
  const sbase0=div(A3).rem,sbase1=div([0,...A3]).rem;
  function sweep(){
   rec.base++;
   for(let j=0;j<rest.length;j++){
    const i=rest[j],Lr=L.tab[i];let z=0;for(let k=0;k<7;k++)z+=values[k]*Lr[k];z%=17;
    if(!(mask[j+7]>>z&1))return;
    D[j]=z;
   }
   rec.poly++;
   const A=Array(7).fill(0);
   for(let i=0;i<7;i++)for(let j=0;j<7;j++)A[j]=(A[j]+values[i]*L.co[i][j])%17;
   const am=div(A).rem;
   for(let u=0;u<17;u++)for(let z=0;z<17;z++){
    const R6=A3.map((x,i)=>R(x*u+(i?A3[i-1]*z:0))).concat(z);
    const coeff=Array.from({length:d},(_,i)=>(u*sbase0[i]+z*sbase1[i])%17);
    const sq=Array(Math.max(0,2*d-1)).fill(0);for(let i=0;i<d;i++)for(let j=0;j<d;j++)sq[i+j]=(sq[i+j]+coeff[i]*coeff[j])%17;
    const rs=div(sq).rem;if(am.some((x,i)=>R(x-2*rs[i])!==0))continue;
    rec.square++;
    const q6=s.map((x,i)=>R(need[i]*(u+z*x)));
    if(Z.some(i=>!q6[i]))continue;
    const S=Array(7).fill(0);for(let i=0;i<4;i++)for(let j=0;j<4;j++)S[i+j]=(S[i+j]+R6[i]*R6[j])%17;
    const q5=div(A.map((x,i)=>R(x-2*S[i]))).H;
    if(q5.slice(5).some(Boolean))continue;
    const g=[],choose=[];
    for(const i of Z){const b=w[i]*q6[i]%17,c=w[i]*v(q5,s[i])%17,h=c*iv[6*b%17]%17;g[i]=G6[b][h];}
    for(let j=0;j<M.length;j++){
     const i=M[j],h=q6[i]*iv[7*q[i]%17]%17;
     if(j<7)g[i]=G7[a[i]][choice[j]][h];
     else choose.push(G7[a[i]].flatMap((hs,k)=>list[j][k]===D[j-7]?[hs[h]]:[]));
    }
    function check(d){
     for(let k=0;k<=d;k++){let x=0;for(let i=0;i<11;i++)x+=g[i][d]*pow[i][k];if(x%17)return false;}
     return true;
    }
    function tail(j){
     if(j<rest.length){for(const c of choose[j]){g[rest[j]]=c;tail(j+1);}return;}
     rec.c5++;
     if(!check(5))throw Error("c5");
     if(!check(4))return;rec.c4++;
     if(g.reduce((n,c)=>n+c[0],0)%17!==3)return;rec.const++;
     if(!check(3)||!check(2)||!check(1))return;
     rec.form++;
     for(let x=0;x<17;x++)for(let y=0;y<17;y++){let z=0;for(let i=0;i<11;i++)z+=v(g[i],R(y-s[i]*x));if(z!==20&&z!==37&&z!==54)return;}
     rec.actual++;
    }
    tail(0);
   }
  }
  function go(j){if(j===7){sweep();return;}for(let k=0;k<list[j].length;k++){choice[j]=k;values[j]=list[j][k];go(j+1);}}
  go(0);
 }
 for(let leadAt=0;leadAt<=2;leadAt++)for(let b=0;b<17;b++)for(let c=0;c<17;c++){
  const a=R(win[leadAt]-b*s[leadAt]-c*s[leadAt]*s[leadAt]);
  if(s.slice(0,leadAt).some(x=>v([a,b,c],x)!==0))continue;
  Qrun([a,b,c]);
 }
 out.push(rec);
}
return {septic:L7.filter(cap).length,sextic:Object.keys(B6).length,out};
}
function eightProfiles(cap=8,total=34){
 const P=17, mod=x=>(x%P+P)%P;
 const inv=new Int32Array(P);for(let i=1;i<P;i++)for(let j=1;j<P;j++)if(i*j%P===1)inv[i]=j;
 const pow=(x,d)=>{let r=1;for(let j=0;j<d;j++)r=r*x%P;return r;};
 const xs=Array.from({length:8},(_,i)=>i+1), ts=xs.map(x=>x*x%P), ts4=ts.map(t=>pow(t,4));
 const L=[1,2,3,4,5,6,7,8]; // reps modulo sign, nth powers = ±1.
 const out=[],hist=[];let odd=0,projected=0,states=0,passesZero=0;
 const lagCache=new Map();
 function lagrange(J){
  const key=J.join(",");
  if(lagCache.has(key))return lagCache.get(key);
  const s=J.map(i=>ts[i]), D=s.map((x,i)=>inv[s.reduce((v,y,j)=>i===j?v:mod(v*mod(x-y)),1)]);
  const a=[0,...ts].map(t=>s.map((x,i)=>s.reduce((v,y,j)=>i===j?v:mod(v*mod(t-y)),D[i])));
  lagCache.set(key,a);return a;
 }
 for(let a=0;a<P;a++)for(let b=0;b<P;b++)for(let c=0;c<P;c++){
  const od=ts.map((t,j)=>mod(xs[j]*(a+b*t+c*t*t)));
  const d=od.map(w=>{const u=2*w%P;return Math.abs(u<=8?u:u-P);});
  const min=d.reduce((u,v)=>u+v,0);if(min>total||d.some(x=>x>cap))continue;
  odd++;
  const J=Array.from({length:8},(_,i)=>i).sort((i,j)=>d[j]-d[i]||i-j).slice(0,4).sort((i,j)=>i-j);
  const K=Array.from({length:8},(_,i)=>i).filter(i=>!J.includes(i));
  projected+=J.reduce((m,i)=>m*(cap-d[i]+1),1)*8;
  const lag=lagrange(J), Y=new Int32Array(4), V=new Int32Array(8);
  for(const l of L){
   function inner(j,left){
    if(j<4){
     const i=J[j], base=d[i]*inv[2]%P;
     for(let k=0;k<=cap-d[i]&&left>=2*k;k++){
      const e=(base+k)%P; V[i]=e;Y[j]=mod(e-l*ts4[i]);inner(j+1,left-2*k);
     }return;
    }
    states++;
    const evalAt=(r)=>{const w=lag[r]; return (w[0]*Y[0]+w[1]*Y[1]+w[2]*Y[2]+w[3]*Y[3])%P;};
    const e0=evalAt(0);if(e0>cap||e0>left)return;
    passesZero++;
    left-=e0;
    for(const i of K){
     const e=(evalAt(i+1)+l*ts4[i])%P;
     const plus=(e+od[i])%P,minus=mod(e-od[i]);
     if(plus>cap||minus>cap)return;
     left-=plus+minus-d[i];if(left<0)return;
     V[i]=e;
    }
    if(left!==0)return;
    const vals=Array(P).fill(e0);for(let i=0;i<8;i++){vals[xs[i]]=(V[i]+od[i])%P;vals[P-xs[i]]=mod(V[i]-od[i]);}
    if(vals.reduce((u,v)=>u+v,0)!==total)throw Error("mass");
    const coef=vals.map((_,j)=>j===0?vals[0]:mod(-vals.reduce((acc,y,x)=>acc+y*pow(x,16-j),0)));
    if(coef[8]!==l || coef[7]!==0 || coef.slice(9).some(Boolean)||vals.some((v,x)=>coef.reduceRight((u,z)=>(u*x+z)%P,0)!==v))throw Error("interp");
    out.push(coef.slice(0,9));
   }
   inner(0,total-min);
  }
 }
 return {odd,projected,states,passesZero,out};
}



function solveOctic(L6,L7,L8,reps){
const P=17,R=x=>(x%P+P)%P,iv=Array(17).fill(0),pw=Array.from({length:17},(_,x)=>Array.from({length:10},(_,d)=>{let z=1;for(let i=0;i<d;i++)z=z*x%17;return z;}));
for(let x=1;x<P;x++)for(let y=1;y<P;y++)if(x*y%17===1)iv[x]=y;
const v=(f,x)=>f.reduceRight((v,u)=>(v*x+u)%17,0), cap=c=>Math.max(...pw.map((_,i)=>v(c,i)))<=6;
const B6=Array.from({length:17},()=>[]), B7=Array.from({length:17},()=>[]), B8=Array.from({length:17},()=>[]);
for(const c of L6.filter(cap))B6[c[6]].push(c);
for(let u=1;u<17;u++)for(const c of L7.filter(cap)){const f=c.map((k,j)=>k*pw[u][j]%17);B7[f[7]].push(f);}
for(const c of L8)for(const u of [1,3]){const f=c.map((k,j)=>k*pw[u][j]%17);B8[f[8]].push(f);}
const bc=Array.from({length:9},(_,d)=>Array.from({length:d+1},(_,j)=>{let r=1;for(let i=1;i<=j;i++)r=r*(d-i+1)*iv[i]%17;return r;}));
function shifted(c,h){return c.map((_,j)=>c.reduce((v,u,d)=>d<j?v:(v+u*bc[d][j]*pw[h][d-j])%17,0));}
const C6=B6.map(a=>a.flatMap(c=>pw.map((_,h)=>shifted(c,h))));
const C7=B7.map(a=>pw.map((_,h)=>a.map(c=>shifted(c,h))));
const C8=B8.map(a=>{
 let ga=Array.from({length:17},()=>Array.from({length:17},()=>[]));
 for(const c of a)for(let h=0;h<17;h++)ga[c[6]][h].push(shifted(c,h));
 return ga;
});
const allowed=B8.map(a=>[...new Set(a.map(c=>c[6]))]);
const out=[];
for(const s of reps){
 const w=s.map((a,i)=>iv[s.reduce((v,b,j)=>i===j?v:R(v*(a-b)),1)]),win=w.map(a=>iv[a]),pp=s.map(x=>pw[x]);
 const C={s,lead:0,base:0,poly:0,phase:0,divisible:0,pool:0,first5:0,c5:0,c4:0,constant:0,formal:0,actual:0};
 const cache=new Map(),I=Array.from({length:11},(_,i)=>i);
 function Lag(J){
  const k=J.join();if(cache.has(k))return cache.get(k);
  const t=J.map(i=>s[i]),co=t.map((x,i)=>{
   let c=[1],d=1;for(let j=0;j<t.length;j++)if(j!==i){
    const y=t[j],b=Array(c.length+1).fill(0);
    for(let k=0;k<c.length;k++){b[k]=R(b[k]-y*c[k]);b[k+1]=R(b[k+1]+c[k]);}
    c=b;d=R(d*(x-y));}
   return c.map(x=>x*iv[d]%17);
  });
  const L=s.map(x=>co.map(c=>v(c,x))),z={co,L};cache.set(k,z);return z;
 }
 for(let r0=0;r0<=1;r0++)for(let lead=1;lead<=8;lead++)for(let q1=0;q1<17;q1++){
  const q0=R(lead*win[r0]-q1*s[r0]);
  if(r0===1&&q0)continue;
  const q=s.map(t=>(q0+q1*t)%17),a=q.map((x,i)=>x*w[i]%17),zero=q.indexOf(0);
  const arm=I.filter(i=>q[i]).sort((i,j)=>allowed[a[i]].length-allowed[a[j]].length||i-j);
  if(arm.some(i=>!allowed[a[i]].length))continue;
  C.lead++;
  const gauge=I.filter(i=>q[i]).slice(0,2).map(i=>s[i]),J=arm.slice(0,5),K=arm.slice(5),L=Lag(J);
  const aa=arm.map(i=>allowed[a[i]].map(z=>z*q[i]*win[i]%17)),bb=aa.map(u=>u.reduce((v,x)=>v|1<<x,0));
  const root=q1?R(-q0*iv[q1]):-1;
  const s0=gauge[0],s1=gauge[1],R0=[R(s0*s1),R(-s0-s1),1],sqR=Array(5).fill(0);
  for(let a=0;a<=2;a++)for(let b=0;b<=2;b++)sqR[a+b]=(sqR[a+b]+R0[a]*R0[b])%17;
  const rw=s.map((t,i)=>v(R0,t));
  const zc=new Int8Array(11),fc=new Int8Array(5);
  function run(){
   C.base++;
   for(let i=0;i<K.length;i++){const j=K[i],r=L.L[j];let zz=0;for(let k=0;k<5;k++)zz+=fc[k]*r[k];zz%=17;if(!(bb[i+5]>>zz&1))return;zc[j]=zz;}
   C.poly++;
   const A=Array(5).fill(0);
   for(let i=0;i<5;i++){zc[J[i]]=fc[i];for(let j=0;j<5;j++)A[j]=(A[j]+fc[i]*L.co[i][j])%17;}
   for(let rho=0;rho<17;rho++){
    C.phase++;
    if(q1&&R(v(A,root)+10*rho*rho*v(sqR,root)))continue;
    const E=A.map((u,i)=>R(u+10*rho*rho*sqR[i]));
    let D=[];
    if(q1){
     for(let i=4;i>=1;i--){D[i-1]=E[i]*iv[q1]%17;E[i-1]=R(E[i-1]-D[i-1]*q0);}
     if(E[0])throw Error("division");
    }else{
     if(E[4])continue;
     D=E.slice(0,4).map(u=>u*iv[q0]%17);
    }
    C.divisible++;
    const pool=new Array(11);
    for(const j of arm){
     const coeff=zc[j]*w[j]*iv[q[j]]%17,h=R(rho*rw[j]*iv[8*q[j]%17]);
     pool[j]=C8[a[j]][coeff][h];
     if(!pool[j].length)throw Error("profile prefix");
    }
    if(zero>=0){
     const lead7=w[zero]*rho*rw[zero]%17,coeff6=w[zero]*v(D,s[zero])%17;
     if(lead7){const h=coeff6*iv[7*lead7%17]%17;pool[zero]=C7[lead7][h];}
     else pool[zero]=C6[coeff6];
     if(!pool[zero].length)continue;
    }
    C.pool++;
    const O=I.slice().sort((a,b)=>pool[a].length-pool[b].length||a-b),J5=O.slice(0,5),K5=O.slice(5),L5=Lag(J5).L;
    const mb=K5.map(i=>pool[i].reduce((a,c)=>a|(1<<(c[5]*win[i]%17)),0)),g=[],f6=[];
    function check(d){
     for(let j=0;j<=d;j++){let n=0;for(let i=0;i<11;i++)n+=g[i][d]*pp[i][j];if(n%17)return false;}
     return true;
    }
    function tail(j){
     if(j<6){const i=K5[j];for(const f of pool[i])if(f[5]*win[i]%17===f6[j]){g[i]=f;tail(j+1);}return;}
     C.c5++;
     if(!check(4))return;C.c4++;
     if(g.reduce((u,v)=>u+v[0],0)%17!==3)return;C.constant++;
     if(!check(3)||!check(2)||!check(1))return;C.formal++;
     for(let x=0;x<17;x++)for(let y=0;y<17;y++){
      let u=0;for(let i=0;i<11;i++)u+=v(g[i],R(y-s[i]*x));
      if(u!==20&&u!==37&&u!==54)return;
     }
     C.actual++;
    }
    function go5(j){
     if(j<5){const i=J5[j];for(const f of pool[i]){g[i]=f;go5(j+1);}return;}
     C.first5++;
     const v=J5.map(i=>g[i][5]*win[i]%17);
     for(let k=0;k<6;k++){const i=K5[k],t=L5[i],val=(v[0]*t[0]+v[1]*t[1]+v[2]*t[2]+v[3]*t[3]+v[4]*t[4])%17;if(!(mb[k]>>val&1))return;f6[k]=val;}
     tail(0);
    }
    go5(0);
   }
  }
  function go(i){if(i===5){run();return;}for(const x of aa[i]){fc[i]=x;go(i+1);}}
  go(0);
 }
 out.push(C);
}
return {octics:L8.length,septics:L7.filter(cap).length,sextics:L6.filter(cap).length,out};
}
function nonics(){
const catalogue=[];
const p=17, R=x=>(x%17+17)%17, pw=(x,k)=>{let y=1;while(k--)y=y*x%17;return y}, iv=Array.from({length:17},(_,i)=>pw(i,15));
const sq=Array.from({length:9},(_,i)=>i*i%17);
const pair=Array.from({length:17},()=>({s:Array(17).fill(-1),e:[],min:99}));
for(let a=0;a<=7;a++)for(let b=0;b<=7;b++){let z=pair[R(a-b)],e=R((a+b)*9);z.s[e]=a+b;z.e.push(e);z.min=Math.min(z.min,a+b);}
const pre=Array.from({length:17**4},(_,i)=>[i%17,Math.floor(i/17)%17,Math.floor(i/289)%17,Math.floor(i/4913)]);
let viable=0,leaves=0,n=0;const hist={},matrix=new Map();
function interp(Q){
 const k=Q.join();
 if(!matrix.has(k))matrix.set(k,sq.map(z=>Q.map(i=>Q.reduce((m,j)=>i===j?m:R(m*(z-sq[j])*iv[R(sq[i]-sq[j])]),1))));
 return matrix.get(k);
}
for(const a of pre){
 const pots=[null];for(let t=1;t<=8;t++){let z=sq[t],o=1;for(let j=3;j>=0;j--)o=R(o*z+a[j]);pots[t]=pair[R(2*t*o)];}
 if(pots.slice(1).reduce((s,v)=>s+v.min,0)>34)continue;
 viable++;
 const Q=Array.from({length:8},(_,i)=>i+1).sort((i,j)=>pots[i].e.length-pots[j].e.length).slice(0,4).sort((a,b)=>a-b);
 const table=interp(Q), ee=Q.map(i=>pots[i].e);
 for(const x0 of ee[0])for(const x1 of ee[1])for(const x2 of ee[2])for(const x3 of ee[3]){
  leaves++;
  const get=(t)=>{const l=table[t];return R(x0*l[0]+x1*l[1]+x2*l[2]+x3*l[3]);};
  const e0=get(0);if(e0>7)continue;
  let s=e0,m=e0,ok=true;
  for(let t=1;t<=8;t++){const e=get(t),y=pots[t].s[e];if(y<0){ok=false;break;}s+=y;if(s>34){ok=false;break;}
    // retrieve the individual pair with this ordinary sum and field midpoint.
    let z=sq[t],o=1;for(let j=3;j>=0;j--)o=R(o*z+a[j]);m=Math.max(m,R(e+t*o),R(e-t*o));
  }
  if(ok&&s===34){n++;hist[m]=(hist[m]??0)+1;
 const ev=Array.from({length:17},(_,x)=>R(get(Math.min(x,17-x))+x*(pw(sq[Math.min(x,17-x)],4)+a.reduce((z,q,j)=>z+q*pw(sq[Math.min(x,17-x)],j),0))));
 if(ev.some(x=>x>7)||ev.reduce((m,x)=>m+x,0)!==34)throw Error("profile");
 const C=ev.map((_,k)=>k===0?ev[0]:R(-ev.reduce((a,u,x)=>a+u*pw(x,16-k),0)));
 if(C[9]!==1||C[8]!==0||C.slice(10).some(Boolean))throw Error("coef");
 catalogue.push(C.slice(0,10));
}
 }
}
let masks=new Set();function gen(i,k,b){if(k===0){masks.add(b);return;}for(let j=i;j<=18-k;j++)gen(j+1,k-1,b|(1<<j));}gen(0,7,0);
const g=[t=>t===17?17:(t+1)%17,t=>t===17?17:3*t%17,t=>t===17?0:t===0?17:iv[t]],os=[];
const act=(n,f)=>{let a=0;for(let i=0;i<18;i++)if(n>>i&1)a|=1<<f(i);return a;};
while(masks.size){const v=masks.values().next().value;const Q=[v];masks.delete(v);for(let j=0;j<Q.length;j++)for(const f of g){const w=act(Q[j],f);if(masks.delete(w))Q.push(w);}os.push(Q.length);}
const result={viable,leaves,n,hist,orbits:os.sort((a,b)=>a-b),masks:os.reduce((a,b)=>a+b,0)};
if(n!==336||viable!==23638||JSON.stringify(hist)!==JSON.stringify({4:7,5:57,6:124,7:148})||os.length!==10||result.masks!==31824)throw Error(JSON.stringify(result));
return {result,catalogue};

}


function solveNonic(L,reps){
const P=17, M=x=>(x%17+17)%17, inverse=Array(17).fill(0);
for(let x=1;x<P;x++)for(let y=1;y<P;y++)if(x*y%P===1)inverse[x]=y;
const XP=Array.from({length:P},(_,x)=>Array.from({length:10},(_,d)=>{let z=1;for(let i=0;i<d;i++)z=z*x%17;return z;}));
const S=Array.from({length:P},()=>[]);
for(let u=1;u<P;u++)for(const c of L){const f=c.map((k,j)=>k*XP[u][j]%17);S[f[9]].push(f);}
const zreturn=[];
for(const s of reps){
 const w=s.map((a,i)=>inverse[s.reduce((v,b,j)=>i===j?v:M(v*(a-b)),1)]),iw=w.map(a=>inverse[a]),D=s.map(x=>XP[x]);
 const f=w.map((a,i)=>S[a].map(c=>c.map(v=>v*iw[i]%P)));
 const list=f.map(c=>{
  const g=Array.from({length:17},()=>Array.from({length:17},()=>[]));
  for(const r of c)g[r[7]][r[6]].push(r);
  return g;
 });
 const opt=list.map(c=>c.map(x=>x.map((pool,v)=>pool.length?v:-1).filter(v=>v>=0)));
 const mask=opt.map(c=>c.map(x=>x.reduce((v,k)=>v|1<<k,0)));
 const C={s,quadratic:0,first4:0,cubic:0,first5:0,quintic:0,quartic:0,constant:0,formal:0,actual:0};
 const vv=new Int8Array(11),g=new Array(11),iota=Array.from({length:11},(_,i)=>i), cache=new Map();
 function Lag(I){
  const sig=I.join();if(cache.has(sig))return cache.get(sig);
  const z=I.map(i=>s[i]),pre=z.map((x,i)=>inverse[z.reduce((v,y,j)=>i===j?v:M(v*(x-y)),1)]);
  const L=s.map(x=>z.map((y,i)=>z.reduce((v,t,j)=>i===j?v:M(v*(x-t)),pre[i])));
  cache.set(sig,L);return L;
 }
 const L0=Lag([0,1,2]);
 const z4=(a,b,c,d,L)=>(a*L[0]+b*L[1]+c*L[2]+d*L[3])%17;
 function check(d){
  for(let j=0;j<=d;j++){let z=0;for(let i=0;i<11;i++)z+=w[i]*g[i][d]*D[i][j];if(z%17)return false;}
  return true;
 }
 for(let a0=0;a0<P;a0++)for(let a1=0;a1<P;a1++)for(let a2=0;a2<P;a2++){
  for(let i=0;i<11;i++){const t=L0[i];vv[i]=(a0*t[0]+a1*t[1]+a2*t[2])%17;}
  C.quadratic++;
  const pos=iota.slice().sort((i,j)=>opt[i][vv[i]].length-opt[j][vv[j]].length||i-j);
  const J=pos.slice(0,4),K=pos.slice(4),l=Lag(J), rows=K.map(i=>l[i]),masks=K.map(i=>mask[i][vv[i]]);
  const A=opt[J[0]][vv[J[0]]],B=opt[J[1]][vv[J[1]]],D1=opt[J[2]][vv[J[2]]],E=opt[J[3]][vv[J[3]]];
  for(const a of A)for(const b of B)for(const c of D1)for(const d of E){
   C.first4++;
   if(!(masks[0]>>z4(a,b,c,d,rows[0])&1)||!(masks[1]>>z4(a,b,c,d,rows[1])&1)||!(masks[2]>>z4(a,b,c,d,rows[2])&1)||!(masks[3]>>z4(a,b,c,d,rows[3])&1)||!(masks[4]>>z4(a,b,c,d,rows[4])&1)||!(masks[5]>>z4(a,b,c,d,rows[5])&1)||!(masks[6]>>z4(a,b,c,d,rows[6])&1))continue;
   C.cubic++;
   const pool=new Array(11);pool[J[0]]=list[J[0]][vv[J[0]]][a];pool[J[1]]=list[J[1]][vv[J[1]]][b];pool[J[2]]=list[J[2]][vv[J[2]]][c];pool[J[3]]=list[J[3]][vv[J[3]]][d];
   for(let j=0;j<7;j++){const i=K[j];pool[i]=list[i][vv[i]][z4(a,b,c,d,rows[j])];}
   const i2=iota.slice().sort((i,j)=>pool[i].length-pool[j].length||i-j),H=i2.slice(0,5),W=i2.slice(5),L=Lag(H);
   const pools=W.map(i=>pool[i].reduce((q,c)=>q|1<<c[5],0));
   const choices=new Int8Array(6);
   const tail=(j)=>{
    if(j<6){const i=W[j],v=choices[j];for(const b of pool[i])if(b[5]===v){g[i]=b;tail(j+1);}return;}
    C.quintic++;
    if(!check(4))return;C.quartic++;
    if(g.reduce((v,c,i)=>v+c[0]*w[i],0)%17!==3)return;C.constant++;
    if(!check(3)||!check(2)||!check(1))return;C.formal++;
    for(let x=0;x<17;x++)for(let y=0;y<17;y++){let v=0;for(let i=0;i<11;i++){
     const t=M(y-s[i]*x);v+=g[i].reduceRight((k,u)=>(k*t+u*w[i])%17,0);}
     if(v!==20&&v!==37&&v!==54)return;
    }
    C.actual++;C.witness=g.map((a,i)=>a.map(v=>v*w[i]%17));
   };
   for(const aa of pool[H[0]])for(const bb of pool[H[1]])for(const cc of pool[H[2]])for(const dd of pool[H[3]])for(const ee of pool[H[4]]){
    C.first5++;
    const A=aa[5],B=bb[5],D1=cc[5],E=dd[5],F=ee[5];
    let good=true;
    for(let j=0;j<6;j++){
     const l=L[W[j]],v=(A*l[0]+B*l[1]+D1*l[2]+E*l[3]+F*l[4])%17;
     if(!(pools[j]>>v&1)){good=false;break;}choices[j]=v;
    }
    if(!good)continue;
    g[H[0]]=aa;g[H[1]]=bb;g[H[2]]=cc;g[H[3]]=dd;g[H[4]]=ee;tail(0);
   }
  }
 }
 zreturn.push(C);
}
return zreturn;
}

const args=typeof process==='undefined'?[]:process.argv.slice(2);
const stages=args.length?args:['g','6','7','8','9'];
const p=slopeOrbits(11);
if(p.reps.length!==10||p.sum!==5005)throw Error("pencil orbit coverage");
const sum=A=>Object.fromEntries(Object.keys(A[0]).filter(k=>typeof A[0][k]==='number').map(k=>[k,A.reduce((v,a)=>v+a[k],0)]));
let l6=null,l7=null;
const catalogues=(once)=>{if(!l6){const r=residues(17,6);if(r.map(x=>x.min).join()!=='136,102,85,51,51,34'||r[5].candidates.length!==16)throw Error('low inventory');l6=r[5].candidates;}
 if(once&&!l7){const z=septic(17,7)[0];if(z.n2!==5)throw Error('septic inventory');l7=z.candidates;}};
for(const k of stages){
  let out=null;
  if(k==='g'){
   out=cap7gate(p.reps);
   if(out.some(o=>o.num)||out.map(o=>o.n).join()!=='125,156,109,111,97,133,101,126,105,117')throw Error('planar cap');
  }else if(k==='6'){
   catalogues(false);out=low6(l6,p.reps);const z=sum(out.out);
   if(z.lead!==53360||z.pre!==170||z.quartic!==152||z.const!==4||z.formal!==0)throw Error('sextic certificate');
  }else if(k==='7'){
   catalogues(true);out=low7(l6,l7,p.reps);const z=sum(out.out);
   if(z.base!==6714090||z.poly!==52833||z.c5!==48537||z.c4!==0)throw Error('septic certificate');
  }else if(k==='8'){
   catalogues(true);const t=eightProfiles(6);
   if(t.out.length!==136||t.odd!==545||t.states!==157560)throw Error('octic inventory');
   out=solveOctic(l6,l7,t.out,p.reps);const z=sum(out.out);
   if(z.c5!==30644||z.c4!==0)throw Error('octic certificate');
  }else if(k==='9'){
   const t=nonics(),q=t.catalogue.filter(c=>Array.from({length:17},(_,x)=>c.reduceRight((v,u)=>(v*x+u)%17,0)).every(y=>y<=6));
   if(t.catalogue.length!==336||q.length!==188)throw Error('nonic inventory');
   out=solveNonic(q,p.reps);const z=sum(out);
   if(z.quintic!==10||z.quartic!==0)throw Error('nonic certificate');
  }else throw Error('use one or more stages g,6,7,8,9');
  console.log(JSON.stringify({stage:k,out}));
}
