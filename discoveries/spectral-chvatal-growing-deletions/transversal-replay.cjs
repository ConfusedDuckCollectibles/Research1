'use strict';
// Independent rational checks of balanced transversal projection and lower layers. No imports.
const ok=(x,label)=>{if(!x)throw Error(label);};
const absB=x=>x<0n?-x:x;
function gcd(a,b){while(b)[a,b]=[b,a%b];return absB(a);}
const lcm=(a,b)=>a/gcd(a,b)*b;
const q=(a,b=1n)=>{a=BigInt(a);b=BigInt(b);ok(b!==0n,'rational denominator');const g=gcd(a,b)*(b<0n?-1n:1n);return [a/g,b/g];};
const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),sub=(a,b)=>add(a,[-b[0],b[1]]);
const mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]);
const sum=a=>a.reduce(add,q(0)),gt=(a,b)=>a[0]*b[1]>b[0]*a[1];
const is=(a,b)=>a[0]===b[0]&&a[1]===b[1];
const show=x=>x.map(String).join('/');
const N=n=>Array.from({length:n},(_,i)=>i);
function solve(A,b){ // exact rational Gauss-Jordan, independent of modular replay
 const n=A.length,h=A.map((row,i)=>[...row.map(v=>Array.isArray(v)?v:q(v)),b[i]]);
 for(let i=0;i<n;i++){
  let p=i;while(p<n&&h[p][i][0]===0n)p++;ok(p<n,'nonsingular incidence fit');[h[i],h[p]]=[h[p],h[i]];
  const z=h[i][i];for(let k=i;k<=n;k++)h[i][k]=div(h[i][k],z);
  for(let j=0;j<n;j++)if(j!==i&&h[j][i][0]!==0n){const y=h[j][i];for(let k=i;k<=n;k++)h[j][k]=sub(h[j][k],mul(y,h[i][k]));}
 }
 return h.map(a=>a[n]);
}

function choose(n,k){if(k<0||k>n)return 0n;let z=1n;for(let j=1;j<=k;j++)z=z*(BigInt(n)-BigInt(k-j))/BigInt(j);return z;}
const bin=(n,k)=>Number(choose(n,k));
function sets(n,k,start=0){if(k===0)return [[]];const z=[];for(let i=start;i+k<=n;i++)for(const s of sets(n,k-1,i+1))z.push([i,...s]);return z;}
function factorPSD(src){
 const H=src.map(x=>x.slice()),n=H.length;let rank=0;
 for(let i=0;i<n;i++){
   const d=H[i][i];ok(d[0]>=0n,'exact real PSD pivot');
   if(d[0]===0n){ok(H[i].slice(i).every(x=>x[0]===0n),'zero pivot row vanishes');continue;}
   rank++;
   for(let j=i+1;j<n;j++)for(let k=j;k<n;k++) H[j][k]=H[k][j]=sub(H[j][k],div(mul(H[j][i],H[k][i]),d));
 }return rank;
}

function cart(s,r){if(!r)return [[]];return N(s).flatMap(j=>cart(s,r-1).map(x=>[j,...x]));}
function top(){
 const s=3,r=4,K=s**(r-1),q0=(s-1)**(r-1),n=r*s,GND=N(n),T=cart(s,r),ix=N(T.length), X=T.map(u=>u.map((a,j)=>j*s+a));
 const fd=4n,g=T.map(([a,b,c,d])=>b===a&&c===a&&d===a?0n:
   b===(a+1)%s&&c===(s-a)%s&&d===(a+2)%s?2n:
   b===(a+2)%s&&c===(a+1)%s&&d===(s-a)%s?3n:4n);
 const dn=g.map(v=>fd*fd-v*v),on=ix.filter(i=>g[i]>0n),off=ix.filter(i=>g[i]===0n);
 const defect=i=>ix.reduce((z,j)=>z+(X[j].includes(i)?dn[j]:0n),0n);
 const d=defect(0), m=BigInt(K)*fd*fd-d;ok(d===35n && GND.every(i=>defect(i)===d),'balanced three matching defects');
 const zn=BigInt(n), z=(i,j)=>BigInt(Math.floor(i/s)===Math.floor(j/s)?r:0)-1n;
 const G=GND.map(i=>GND.map(j=>q(ix.reduce((v,a)=>v+(X[a].includes(i)&&X[a].includes(j)?g[a]*g[a]:0n),0n),fd*fd)));
 const lifted=G.map((R,i)=>R.map((a,j)=>add(a,q(z(i,j),zn))));
 const Gp=GND.map(i=>solve(lifted,GND.map(j=>q(Number(j===i)))).map((a,j)=>sub(a,q(z(i,j),zn))));
 let piTests=0,hTests=0,kTests=0,symTests=0;
 for(const i of GND)for(const j of GND){
   ok(is(sum(GND.map(k=>mul(G[i][k],Gp[k][j]))),sub(q(Number(i===j)),q(z(i,j),zn))),'weighted product pseudo inverse');
   ok(is(Gp[i][j],Gp[j][i]),'selfadjoint GF inverse');symTests++;
   ok(is(sum(GND.map(k=>mul(Gp[i][k],q(z(k,j),zn)))),q(0)),'kernel annihilation');kTests++;
 }
 const gd=Gp.flat().reduce((a,x)=>lcm(a,x[1]),1n);
 const gn=Gp.map(A=>A.map(x=>x[0]*(gd/x[1])));
 const pd=fd*fd*gd;
 const Pi=ix.map(i=>ix.map(j=>g[i]*g[j]*X[i].reduce((a,v)=>a+X[j].reduce((b,w)=>b+gn[v][w],0n),0n)));
 const E=ix.map(i=>ix.map(j=>(i===j?pd:0n)-Pi[i][j]));
 const dense=(A,B)=>ix.map(i=>ix.map(j=>ix.reduce((v,k)=>v+A[i][k]*B[k][j],0n)));
 ok(dense(Pi,Pi).every((row,i)=>row.every((a,j)=>a===Pi[i][j]*pd)),'exact incidence orthogonal projection');
 const disj=(i,j)=>!X[i].some(x=>X[j].includes(x)), zD=BigInt(K*q0);
 const Ln=(i,j)=>BigInt(q0)-(disj(i,j)?BigInt(K):0n);
 const R=ix.map(i=>ix.map(j=>(i===j?fd*fd*zD:0n)-g[i]*g[j]*Ln(i,j)));
 const rD=fd*fd*zD, ER=dense(E,R), ERE=dense(ER,E), uD=pd*pd*rD;
 const RX=ix.map(i=>GND.map(j=>ix.reduce((v,k)=>v+(X[k].includes(j)?R[i][k]*g[k]:0n),0n)));
 const LX=ix.map(i=>GND.map(j=>ix.reduce((v,k)=>v+(X[k].includes(j)?Ln(i,k)*dn[k]:0n),0n)));
 for(const i of ix)for(const j of GND){
   ok(ix.reduce((v,k)=>v+(X[k].includes(j)?Pi[i][k]*g[k]:0n),0n)===(X[i].includes(j)?g[i]*pd:0n),'Pi fixes weighted incidence');piTests++;
   ok(RX[i][j]===g[i]*LX[i][j],'first-order projection identity');hTests++;
 }
 for(const i of GND)for(const j of GND)
   ok(ix.reduce((v,k)=>v+(X[k].includes(i)?g[k]*RX[k][j]:0n),0n)===
   ix.reduce((v,k)=>v+(X[k].includes(i)?dn[k]*((X[k].includes(j)?rD:0n)-LX[k][j]):0n),0n),'second-order projected seed identity');
 function Q0(Mn,Md){return ix.map(i=>ix.map(j=>disj(i,j)?0n:
   (i===j?uD*fd*fd*Mn:0n)-uD*g[i]*g[j]*Md-ERE[i][j]*fd*fd*Mn));}
 const Q=Q0(m,fd*fd), bad=Q0(BigInt(K),1n);
 const qc=(A,i,j)=>ix.reduce((v,k)=>v+(X[k].includes(j)?A[i][k]*g[k]:0n),0n);
 let locks=0,badLock=0;
 const EQE=dense(dense(E,Q),E), AQ=ix.map(i=>ix.map(j=>disj(i,j)?0n:EQE[i][j]));
 for(const i of ix)for(const j of X[i]){
   ok(qc(Q,i,j)===0n && qc(AQ,i,j)===0n,'signed offset and next iterate compatible');locks++;
   if(qc(bad,i,j)!==0n)badLock++;
 }
 ok(badLock===on.length*r,'all positive incidence locks fail with omitted shift');
 for(const i of off)for(const j of ix)ok(Q[i][j]===0n&&AQ[i][j]===0n&&
   R[i][j]===(i===j?rD:0n),'zero rows isolate identity and no corrections');
 // Independently attain actual seed gap via a tensor with three mean-zero factors.
 const y=T.map(a=>a.slice(0,3).reduce((z,j)=>z*[1n,-1n,0n][j],1n));
 const R0=ix.map(i=>ix.map(j=>(i===j?zD:0n)-Ln(i,j)));
 const zY=R0.map(A=>A.reduce((v,a,k)=>v+a*y[k],0n));
 const gap=sub(q(1),q(1,(s-1)**2));
 ok(zY.every((x,i)=>is(q(x,zD),mul(gap,q(y[i])))),'three-factor sharp gap ray');
 return {r,s,K,d:show(q(d,fd*fd)),M:show(q(m,fd*fd)),positiveRows:on.length,zeroRows:off.length,
   pseudoinverseEntries:symTests,kernel: kTests,projectorIdempotence:ix.length**2, incidenceProjectionEntries:piTests,
   firstOrder:hTests,secondOrder:n*n,compatibleInitialAndNext:locks,offsetOmissionFailures:badLock,sharpSeedGap:show(gap)};
}
function lower(r,s,k,useFull){
 const supports=sets(r,k),n=r*s,v=N(n),b=bin(r-1,k-1),c=bin(r-2,k-2),
    eps=q(b,s*(b-c)+k*c),gam=mul(eps,q(b-c)),overlap=(S,T)=>S.filter(x=>T.includes(x)).length;
 const w=(S,T)=>sub(q(1),mul(eps,q(overlap(S,T))));
 const face=supports.flatMap(S=>cart(s,k).map(J=>({S,U:S.map((p,i)=>p*s+J[i])})));
 const mat=(a,b,i,j)=>i===j?q(1):a.U.some(x=>b.U.includes(x))?q(0):
    div(w(a.S,b.S),mul(gam,q((s-1)**overlap(a.S,b.S)*s**(k-overlap(a.S,b.S)))));
 // Independent exact PSD via direct Schur only on modest full cases.
 let psdRank=null;
 if(useFull)psdRank=factorPSD(face.map((a,i)=>face.map((b,j)=>mat(a,b,i,j))));
 // All label sectors, with exact Schur and diagonal star-component nullity.
 let sectorRank=0,blocks=[];
 for(let j=0;j<=k;j++){
   const u=N(j),S=supports.filter(T=>u.every(x=>T.includes(x)));
   const c=div(q((-1)**j),mul(gam,q((s-1)**j)));
   const M=S.map((A,i)=>S.map((B,h)=>add(q(Number(i===h)),mul(c,w(A,B)))));
   const rank=factorPSD(M),mult=bin(r,j)*(s-1)**j;
   sectorRank+=rank*mult;blocks.push({j,rank,mult,dim:S.length});
 }
 if(useFull)ok(sectorRank===psdRank,'direct and sector rational PSD ranks agree');
 const hd=face.reduce((z,a,i)=>face.reduce((u,b,j)=>lcm(u,mat(a,b,i,j)[1]),z),1n);
 let checked=0;
 for(let i=0;i<face.length;i++){const total=v.map(()=>0n);for(let j=0;j<face.length;j++){
   const [a,b]=mat(face[i],face[j],i,j);const m=a*(hd/b);
   if(m)for(const p of face[j].U)total[p]+=m;
 }
 ok(total.every(x=>x===hd),'all lower rank stars map exactly to ones');checked+=n;
 }
 const star=b*s**(k-1),cap=Math.floor(r*s/k),frac=q(face.length,cap);
 ok(is(q(face.length*k,r*s),q(star)),'induced star size');
 if(r*s%k)ok(gt(frac,q(star)),'fractional clique cover strictly worse');
 return {r,s,k,faces:face.length,theta:star,cliqueCover:show(frac),psdRank,sectorRank,starEntries:checked,blocks};
}
function bounds(){
 const large=[];
 for(const r of [3,4,5,9,32,127]){
   const s=BigInt(100*r**3),R=BigInt(r),K=s**(R-1n),Q=(s-1n)**(R-1n);
   const d=q(K,10n*R**3n),mu=sub(q(1),q(1,s-1n));
   const drift=div(q(3,K),sub(q(1),mul(q(3n*R,K),d)));
   const rho=q((2n+4n*R)*4n*R**2n,s),corr=div(mul(q(6n*R**3n,K),d),sub(q(1),rho));
   ok(!gt(q(99,100),q(Q,K))&&!gt(drift,q(4,K))&&!gt(rho,q(1,5))&&!gt(corr,q(3,4))&&gt(mu,corr),'cutoff inequalities');
   large.push({r,s:String(s),rho:show(rho),correction:show(corr),mu:show(mu)});
 }
 // Remove all diagonal-orbit labels through a fixed admissible (r-1)-prefix.
 // Each removed orbit is a perfect matching; all its points have deficit one.
 const r=3,s=2701n,K=s*s,d=s,t=K-d;
 ok(10n*27n*d<=K,'zero-codegree control inside theorem');
 ok((r*Number(s))%2===1,'lower divider obstruction in same control');
 // Four even-parity triples on [2]^3 have degree 2 and no two disjoint.
 const par=cart(2,3).filter(a=>a.reduce((x,y)=>x+y,0)%2===0);
 ok(par.length===4 && N(3).every(j=>N(2).every(i=>par.filter(a=>a[j]===i).length===2))&&
   par.every(a=>par.every(b=>a.some((u,j)=>u===b[j]))),'regularity alone negative control');
 return {large,zeroAdmissibleCodegree:{r,s:String(s),deletedOrbits:String(d),topTheta:String(t),thetaWithLower:String(t+1n+2n*s)},
    regularityAlone:{r:3,s:2,degree:2,theta:4,perfectMatchings:0}};
}
const out={top:top(),lower:[lower(3,3,2,true),lower(5,3,2,true),lower(4,5,3,false),lower(5,3,3,false)],
    analytic:bounds()};
console.log(JSON.stringify(out,null,2));
