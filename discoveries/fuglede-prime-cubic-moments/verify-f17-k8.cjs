"use strict";

function ledger(){
 const oddVal=[]; const out=[];
 for(let m=2;m<=16;m++){
  const quo=17-m, lp=17*quo*quo<256;
  let n=m, fac=[];
  for(let q=2;q<=n;q++)if(n%q===0){let v=0;while(n%q===0){n/=q;v++;}let r=1,j=0;do{r=r*q%17;j++;if(j>16)throw Error("order")}while(r!==16);fac.push([q,v,j]);}
  const norm=m%2===1&&fac.some(([q,v])=>v%2===1);
  if(!lp&&!norm)out.push(m);
  oddVal.push({m,N:17*m,lp,norm,fac});
 }
 if(JSON.stringify(out)!=="[2,4,6,8,9,10,12]")throw Error("ledger");
 return {oddVal,out,sizes:out.map(x=>17*x)};
}

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

function eightSlopes(list){
 const p=17, r=x=>(x%p+p)%p, mul=(a,b)=>a*b%p;
 const inv=Array(17).fill(0);for(let i=1;i<17;i++)for(let j=1;j<17;j++)if(i*j%17===1)inv[i]=j;
 const by=Object.fromEntries(list.map(c=>[c[6],c]));
 let masks=0,center=0,pass0=0,pass2=0,pass4=0,last=[],sqFans=0;
 function check(s){
  masks++;
  const w=s.map((a,i)=>inv[s.reduce((v,b,j)=>i===j?v:r(v*r(a-b)),1)]);
  let any=false;
  for(let a=1;a<17;a++){
   center++; const G=w.map(u=>by[a*u%p]);
   if(G.some(x=>!x))throw Error("incomplete catalogue");
   if(G.reduce((v,c)=>v+c[0],0)%17!==14)continue;
   pass0++; any=true;
   function zero(d){
    for(let k=0;k<=d;k++){
     let x=0;
     for(let i=0;i<8;i++){let v=1;for(let h=0;h<k;h++)v=v*s[i]%p;x=(x+v*G[i][d])%p;}
     if(x!==0)return false;
    }
    return true;
   }
   if(!zero(4))continue; pass4++;
   if(!zero(2))continue; pass2++;
   last.push({s:s.slice(),a,G});
  }
  if(any)sqFans++;
 }
 let s=[];
 function go(i){
  if(s.length===8){check(s);return;}for(let j=i;j<=17-(8-s.length);j++){s.push(j);go(j+1);s.pop();}
 }
 go(0);
 return {masks,center,pass0,pass4,pass2,sqFans,last};
}

function nineSlopes(L6,L7){
 const P=17, mod=x=>(x%P+P)%P, inv=new Int32Array(P), power=Array.from({length:P},(_,x)=>Array.from({length:8},(_,d)=>{let r=1;for(let j=0;j<d;j++)r=r*x%P;return r;}));
 for(let i=1;i<P;i++)for(let j=1;j<P;j++)if(i*j%P===1)inv[i]=j;
 const evaluate=(c,x)=>c.reduceRight((a,v)=>(a*x+v)%P,0);
 const good=arr=>arr.filter(c=>Math.max(...power.map((_,x)=>evaluate(c,x)))<=9);
 L6=good(L6);L7=good(L7);
 const B6=Object.fromEntries(L6.map(c=>[c[6],c]));
 const B7=Array.from({length:P},()=>new Map());
 for(let u=1;u<P;u++){const W=power[u];
  for(const c of L7){const q=c.map((v,d)=>v*W[d]%P);if(!B7[q[7]].has(q[5]))B7[q[7]].set(q[5],[]); B7[q[7]].get(q[5]).push(q);}
 }
 const bc=Array.from({length:8},(_,d)=>Array.from({length:d+1},(_,j)=>{let r=1;for(let u=1;u<=j;u++)r=r*(d-u+1)*inv[u]%P;return r;}));
 const shifted=Array.from({length:17},()=>[]);
 for(let a=1;a<17;a++)if(B6[a])for(let h=0;h<17;h++){
  const co=B6[a],hp=power[h], out=new Int32Array(7);
  for(let j=0;j<=6;j++)for(let d=j;d<=6;d++)out[j]=(out[j]+bc[d][j]*hp[d-j]*co[d])%P;
  shifted[a][h]=out;
 }
 const s=[];
 const C={masks:0,n6free:0,n6const:0,n6deg4:0,n6deg3:0,n6deg2:0,n6deg1:0,n7free:0,n7fits:0,n7const:0,n7deg4:0,n7deg3:0,n7deg2:0,n7deg1:0,surv6:[],surv7:[]};
 function check(){
  C.masks++;
  const w=s.map((a,i)=>inv[s.reduce((v,b,j)=>i===j?v:mod(v*mod(a-b)),1)]),win=w.map(u=>inv[u]),pow=s.map(x=>power[x]);
  function zero(G,d){
   for(let k=0;k<=d;k++){
    let a=0;for(let i=0;i<9;i++)a+=pow[i][k]*G[i][d];
    if(a%P!==0)return false;
   }
   return true;
  }
  const w3=s.slice(0,3).map((a,i)=>inv[s.slice(0,3).reduce((v,b,j)=>i===j?v:mod(v*mod(a-b)),1)]);
  const lag=s.map(x=>s.slice(0,3).map((y,i)=>s.slice(0,3).reduce((v,b,j)=>i===j?v:mod(v*mod(x-b)),w3[i])));
  const g0=B7[w[0]],g1=B7[w[1]],g2=B7[w[2]];
  for(const a of [...g0.values()].flat())for(const b of [...g1.values()].flat())for(const c of [...g2.values()].flat()){
   C.n7free++;
   let G=[a,b,c], choices=[], all=true;
   const A=a[5]*win[0]%P,B=b[5]*win[1]%P,D=c[5]*win[2]%P;
   for(let i=3;i<9;i++){
    const l=lag[i], r=(A*l[0]+B*l[1]+D*l[2])*w[i]%P;
    const q=B7[w[i]].get(r);if(!q){all=false;break;}choices.push(q);
   }
   if(!all)continue;
   function dfs(k,total){
    if(k<6){for(const q of choices[k]){G.push(q);dfs(k+1,total+q[0]);G.pop();}return;}
    C.n7fits++;
    if(total%P!==16)return;
    C.n7const++;
    for(let d=4;d>=1;d--){if(!zero(G,d))break; C["n7deg"+d]++; if(d===1)C.surv7.push({s:s.slice(),G:G.slice()});}
   }
   dfs(0,a[0]+b[0]+c[0]);
  }
  for(const lead0 of [1,3]){
   const A=lead0*win[0]%P;
   for(let B=0;B<P;B++){
    const a=s.map((x,i)=>(A+B*x)*w[i]%P);
    if(a.some(x=>!B6[x]))continue;
    for(let F=0;F<P;F++){
     C.n6free++;
     const G=[],v5=s.map((x,i)=>F*x*mod(x-1)*w[i]%P);
     for(let i=0;i<9;i++){
      const h=v5[i]*inv[6*a[i]%P]%P;
      const out=shifted[a[i]][h];
      if(out[5]!==v5[i])throw Error("translation");
      G.push(out);
     }
     if(G.reduce((v,x)=>v+x[0],0)%P!==16)continue;
     C.n6const++;
     for(let d=4;d>=1;d--){if(!zero(G,d))break; C["n6deg"+d]++; if(d===1)C.surv6.push({s:s.slice(),G:G.map(q=>[...q])});}
    }
   }
  }
 }
 s.push(0,1);
 function go(i){
  if(s.length===9){check();return;}for(let j=i;j<=17-(9-s.length);j++){s.push(j);go(j+1);s.pop();}
 }
 go(2);
 return {L6:L6.length,L7:L7.length,C};
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
function octicGate(profs,reps){
 const P=17, mod=x=>(x%P+P)%P, inv=Array(17).fill(0);
 const pow=Array.from({length:17},(_,x)=>Array.from({length:9},(_,d)=>{let z=1;for(let i=0;i<d;i++)z=z*x%17;return z;}));
 for(let x=1;x<P;x++)for(let y=1;y<P;y++)if(x*y%P===1)inv[x]=y;
 const B=Array.from({length:17},()=>[]);
 const M=Array.from({length:17},()=>Array.from({length:17},()=>Array.from({length:17},()=>[])));
 for(const a of profs)for(const u of [1,3]){const b=a.map((v,d)=>v*pow[u][d]%P);
  B[b[8]].push(b);M[b[8]][b[6]][b[5]].push(b);}
 const cnt={orbits:reps.length,trials:0,first3:0,c6matches:0,first4:0,c5matches:0,full:0,const:0,quartic:0,odd3:0,quad:0,odd1:0,surv:[]};
 for(const s of reps){
  const w=s.map((a,i)=>inv[s.reduce((v,b,j)=>i===j?v:mod(v*mod(a-b)),1)]),wr=w.map(u=>inv[u]);
  const qcache=new Map();
  const Lagrange=J=>{
   const k=J.join(",");if(qcache.has(k))return qcache.get(k);
   const z=J.map(i=>s[i]);const q=s.map(x=>z.map((u,i)=>z.reduce((v,t,j)=>j===i?v:mod(v*mod(x-t)*inv[mod(u-t)]),1)));
   qcache.set(k,q);return q;
  };
  for(let l=1;l<=8;l++){
   cnt.trials++;
   const a=w.map(x=>l*x%P);
   const J=Array.from({length:10},(_,i)=>i).sort((x,y)=>B[a[x]].length-B[a[y]].length||x-y);
   const L3=Lagrange(J.slice(0,3)),L4=Lagrange(J.slice(0,4));
   const G=[];
   for(const a0 of B[a[J[0]]])for(const a1 of B[a[J[1]]])for(const a2 of B[a[J[2]]]){
    cnt.first3++;G[J[0]]=a0;G[J[1]]=a1;G[J[2]]=a2;
    const v=J.slice(0,3).map(i=>G[i][6]*wr[i]%P);
    const rr=L3.map((z,i)=>(z[0]*v[0]+z[1]*v[1]+z[2]*v[2])*w[i]%P);
    const S=J.slice(3).map(i=>B[a[i]].filter(c=>c[6]===rr[i]));
    if(S.some(q=>!q.length))continue;
    cnt.c6matches++;
    for(const a3 of S[0]){
     cnt.first4++;G[J[3]]=a3;
     const v=J.slice(0,4).map(i=>G[i][5]*wr[i]%P);
     const rr5=L4.map((z,i)=>(z[0]*v[0]+z[1]*v[1]+z[2]*v[2]+z[3]*v[3])*w[i]%P);
     const T=J.slice(4).map(i=>M[a[i]][rr[i]][rr5[i]]);
     if(T.some(r=>!r.length))continue;cnt.c5matches++;
     function extend(j,total){
      if(j<6){for(const c of T[j]){G[J[4+j]]=c;extend(j+1,total+c[0]);}return;}
      cnt.full++;if(total%17!==1)return;cnt.const++;
      const names=["odd1","quad","odd3","quartic"];
      for(let d=4;d>=1;d--){
       for(let r=0;r<=d;r++)if(G.reduce((v,g,i)=>v+g[d]*pow[s[i]][r],0)%17!==0)return;
       cnt[names[d-1]]++;
      }
      cnt.surv.push({s,l,G:G.slice()});
     }
     extend(0,J.slice(0,4).reduce((n,i)=>n+G[i][0],0));
    }
   }
  }
 }
 return cnt;
}

function belowOctic(L6,L7,reps){
 const P=17, mod=x=>(x%P+P)%P, inv=Array(17).fill(0);
 const pow=Array.from({length:P},(_,x)=>Array.from({length:8},(_,d)=>{let r=1;for(let i=0;i<d;i++)r=r*x%P;return r;}));
 for(let x=1;x<P;x++)for(let y=1;y<P;y++)if(x*y%P===1)inv[x]=y;
 const value=(c,x)=>c.reduceRight((v,u)=>(v*x+u)%P,0);
 L6=L6.filter(c=>Math.max(...pow.map((_,x)=>value(c,x)))<=8);
 L7=L7.filter(c=>Math.max(...pow.map((_,x)=>value(c,x)))<=8);
 const B6=Object.fromEntries(L6.map(c=>[c[6],c]));
 const B7=Array.from({length:P},()=>[]);
 for(let u=1;u<P;u++)for(const c of L7){const g=c.map((x,i)=>x*pow[u][i]%P);B7[g[7]].push(g);}
 const bc=Array.from({length:8},(_,d)=>Array.from({length:d+1},(_,j)=>{let r=1;for(let i=1;i<=j;i++)r=r*(d-i+1)*inv[i]%P;return r;}));
 function shift(g,h){
  const q=Array(g.length).fill(0);
  for(let j=0;j<g.length;j++)for(let d=j;d<g.length;d++)q[j]=(q[j]+bc[d][j]*g[d]*pow[h][d-j])%P;
  return q;
 }
 const D6=Array.from({length:P},(_,a)=>B6[a]?pow.map((_,h)=>shift(B6[a],h)):[]);
 const D7=Array.from({length:P},(_,a)=>pow.map((_,h)=>B7[a].map(g=>shift(g,h))));
 const C={orbits:reps.length,eligible6:L6.length,eligible7:L7.length,checks6:0,const6:0,quartic6:0,full6:0,first7:0,fits7:0,const7:0,quartic7:0,full7:0};
 for(const s of reps){
  const w=s.map((a,i)=>inv[s.reduce((v,b,j)=>i===j?v:mod(v*mod(a-b)),1)]),wr=w.map(u=>inv[u]);
  function identity(G,j){
   for(let e=0;e<=j;e++){let z=0;for(let i=0;i<10;i++)z+=G[i][j]*pow[s[i]][e];if(z%P!==0)return false;}
   return true;
  }
  for(const lead of [1,3])for(let b=0;b<P;b++)for(let c=0;c<P;c++){
   const a=lead*wr[0]%P;
   const leading=s.map((x,i)=>(a+b*x+c*x*x)*w[i]%P);
   if(leading.some(x=>!B6[x]))continue;
   const F=s.map((x,i)=>6*leading[i]%P===0?0:mod(x*(x-1)*w[i]*inv[6*leading[i]%P]));
   const G=[];
   for(let u=0;u<P;u++)for(let v=0;v<P;v++){
    let sum=0;for(let i=0;i<10;i++){const h=F[i]*(u*s[i]+v)%P;G[i]=D6[leading[i]][h];sum+=G[i][0];}
    C.checks6++;if(sum%P!==1)continue;C.const6++;
    if(!identity(G,4))continue;C.quartic6++;
    if(!identity(G,3)||!identity(G,2)||!identity(G,1))continue;C.full6++;
   }
  }
  const seen=new Set();
  for(let pos=0;pos<2;pos++)for(let b=0;b<P;b++){
   const a=mod(wr[pos]-b*s[pos]);
   if(pos===1&&a!==0)continue;
   const l=s.map((x,i)=>(a+b*x)*w[i]%P);
   if(l[pos]!==1||l.indexOf(l.find(x=>x!==0))!==pos)throw Error("gauge");
   const k0=l.indexOf(0);
   const J=Array.from({length:10},(_,i)=>i).filter(i=>i!==k0);
   const i0=J[0],i1=J[1],z=J.slice(0,4).map(i=>s[i]);
   const den=z.map((x,i)=>inv[z.reduce((v,y,j)=>i===j?v:mod(v*mod(x-y)),1)]);
   const Q=s.map(x=>z.map((y,i)=>z.reduce((v,t,j)=>i===j?v:mod(v*mod(x-t)),den[i])));
   for(let f=0;f<P;f++){
    const g6=s.map((x,i)=>mod(f*(x-s[i0])*(x-s[i1])*w[i]));
    const opts=l.map((t,i)=> t?D7[t][g6[i]*inv[7*t%P]%P]:D6[g6[i]]);
    if(opts.some(x=>!x.length))continue;
    const G=[];
    function head(j){
     if(j<4){const i=J[j];for(const q of opts[i]){G[i]=q;head(j+1);}return;}
     C.first7++;
     const coeff=J.slice(0,4).map(i=>G[i][5]*wr[i]%P), choices=[];
     const remain=Array.from({length:10},(_,i)=>i).filter(i=>!J.slice(0,4).includes(i));
     for(const i of remain){const z=Q[i], c=(z[0]*coeff[0]+z[1]*coeff[1]+z[2]*coeff[2]+z[3]*coeff[3])*w[i]%P;
      const q=opts[i].filter(g=>g[5]===c);if(!q.length)return;choices.push(q);}
     function tail(j,total){
      if(j<remain.length){for(const q of choices[j]){G[remain[j]]=q;tail(j+1,total+q[0]);}return;}
      C.fits7++;if(total%P!==1)return;C.const7++;
      if(!identity(G,4))return;C.quartic7++;
      if(!identity(G,3)||!identity(G,2)||!identity(G,1))return;C.full7++;
     }
     tail(0,J.slice(0,4).reduce((a,i)=>a+G[i][0],0));
    }
    head(0);
   }
  }
 }
 return C;
}
const p=eightProfiles(), s=slopeOrbits(10), a=octicGate(p.out,s.reps);
const b=belowOctic(residues(17,6)[5].candidates,septic(17,7)[0].candidates,s.reps);
if(p.out.length!==284||s.reps.length!==17||s.sum!==6435||a.const!==5652||a.quartic!==0||b.const6!==96335||b.const7!==115||b.quartic6!==0||b.quartic7!==0)throw Error("audit changed");
const {out,...summary}=p;
console.log(JSON.stringify({summary,n:out.length,orbits:s,a,b}));
