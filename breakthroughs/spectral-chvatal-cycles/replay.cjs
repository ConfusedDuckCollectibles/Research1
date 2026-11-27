
"use strict";
const ck=(b,s)=>{if(!b)throw Error(s)};
function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){[a,b]=[b,a%b]}return a;}
function R(p=0,q=1){p=BigInt(p);q=BigInt(q);ck(q!==0n,"division");if(q<0n){p=-p;q=-q}const d=gcd(p,q);return [p/d,q/d]}
const z=R(),o=R(1),f=(x)=>x[0]+"/"+x[1];
const plus=(a,b)=>R(a[0]*b[1]+b[0]*a[1],a[1]*b[1]), neg=a=>R(-a[0],a[1]), sub=(a,b)=>plus(a,neg(b)), mul=(a,b)=>R(a[0]*b[0],a[1]*b[1]),div=(a,b)=>R(a[0]*b[1],a[1]*b[0]),sc=(a,b)=>mul(a,R(b)), sq=a=>mul(a,a);
const abs=a=>a[0]>=0n?a:neg(a),inv=a=>div(o,a),eq=(a,b)=>sub(a,b)[0]===0n,le=(a,b)=>sub(a,b)[0]<=0n;
const arr=(n,g)=>Array.from({length:n},(_,i)=>g(i)),sum=(a)=>a.reduce(plus,z);
const c2=x=>x*(x-1)/2,c3=x=>x*(x-1)*(x-2)/6;
function constants(n,r){
 const K=c2(n-1), a=c2(n-2), b=c2(n-3), q=c2(n-4), L=sub(R(K),sc(r,2));
 const c=plus(R(n-2,K*a),R(12,a*(b-12))),v=R(6*(n-2),b-12);
 const B=sum([sc(plus(R(3,K*K),mul(r,sq(c))),2), div(sc(sq(plus(R(3,K),mul(r,c))),4),L),div(plus(R(3),mul(r,sq(v))),R(K*(K-1)))]);
 const P=plus(sc(plus(div(plus(R(6,K),sc(sub(R(1,q),R(1,K)),9)),R(a)),R(120,a*(a-12)*q)),2),div(sc(plus(o,div(r,R(q))),16),R(a*a)));
 const M=c3(n)*(c3(n)-c3(n-3)),gap=R(q-1,q),kappa=R(a-6-6*(n-4),a-6);
 const lhs=sc(sq(plus(B,P)),M),rhs=sq(mul(gap,kappa));
 return {B,P,M,gap,kappa,lhs,rhs,pass:le(lhs,rhs)&&!eq(lhs,rhs)};
}
function solve(H,b){
 H=H.map((l,i)=>l.concat([b[i]]));const n=b.length;
 for(let j=0;j<n;j++){
  const p=H.findIndex((a,i)=>i>=j&&a[j][0]!==0n);ck(p>=0,"invertible");
  [H[p],H[j]]=[H[j],H[p]];const u=inv(H[j][j]);
  for(let t=j;t<=n;t++)H[j][t]=mul(H[j][t],u);
  for(let i=0;i<n;i++)if(i!==j&&H[i][j][0]!==0n){
   const v=H[i][j];for(let t=j;t<=n;t++)H[i][t]=sub(H[i][t],mul(v,H[j][t]));
  }
 }
 return H.map(a=>a[n]);
}
function triples(n){
 const T=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++)T.push([i,j,k]);return T;
}
function construction(cycles,n){
 const r=cycles.reduce((s,k)=>s+k,0),a=c2(n-2),b=c2(n-3),K=c2(n-1),L=K-2*r;
 const old=arr(n,j=>j),tip=arr(r,i=>i+r),isTip=j=>j>=r&&j<2*r;
 const use=triples(n), removed=new Set();let off=0;
 for(const l of cycles){for(let t=0;t<l;t++){const u=off+t,v=off+(t+1)%l;removed.add([u,v,r+u].sort((a,b)=>a-b).join());}off+=l;}
 const U=use.filter(S=>!removed.has(S.join())), inc=(j,T=U)=>T.map((S,k)=>S.includes(j)?k:-1).filter(k=>k>=0);
 const col=old.map(j=>inc(j)),sumIx=(u,J)=>sum(J.map(i=>u[i])),AT=(u,j)=>sumIx(u,col[j]);
 function adjust(i,base,target){
  const J=old.filter(j=>j!==i),T=U.filter(S=>!S.includes(i)),S=J.map(j=>inc(j,T));
  const G=S.map(L=>J.map(j=>R(L.filter(k=>T[k].includes(j)).length)));
  const demand=J.map(j=>sub(target[j],AT(base,j))),x=solve(G,demand);
  const U1=U.map((S,k)=>S.includes(i)?base[k]:plus(base[k],sum(S.map(j=>x[J.indexOf(j)]))));
  ck(old.every(j=>eq(AT(U1,j),target[j])),"incidence correction");
  return U1;
 }
 const lam=R(K-1,K),C=arr(r,i=>adjust(i,U.map(S=>S.includes(i)?z:R(1,a)),old.map(j=>j===i?z:isTip(j)?lam:o)));
 const V=tip.map(i=>adjust(i,U.map(S=>S.includes(i)?o:z),old.map(j=>j===i?R(K-1):z)));
 let claims=0;
 for(let u=0;u<U.length;u++){
  const Cs=C.map(A=>A[u]),g=sub(o,sc(sum(Cs),2)),Vs=V.map(A=>A[u]),W=Cs.map(x=>sub(x,R(1,K))),h=sum(W);
  ck(U[u].filter(i=>i<r).every(i=>eq(Cs[i],z)),"intersection zeros");
  for(let i=0;i<r;i++){
   ck(le(abs(Cs[i]),plus(R(1,a),R(12,a*(b-12)))),"C entry bound");
   ck(le(abs(W[i]),U[u].includes(i)?R(1,K):plus(R(n-2,K*a),R(12,a*(b-12)))),"W entry bound");
   ck(le(abs(Vs[i]),U[u].includes(i+r)?o:R(6*(n-2),b-12)),"V entry bound");
  }
  for(const j of old){
   const ct=C.map(A=>AT(A,j)),gv=sub(R(col[j].length),sc(sum(ct),2));
   const B0=plus(sc(sum(Cs.map((x,i)=>mul(x,ct[i]))),2),div(mul(g,gv),R(L)));
   const B=plus(B0,isTip(j)?div(V[j-r][u],R(K)):z);
   const want=isTip(j)?plus(lam,div(V[j-r][u],R(K))):j<r?sub(o,sc(Cs[j],2)):o;
   ck(eq(B,want),"B applied to whole incidence");
   if(U[u].includes(j)){ck(eq(B,o),"locked incident compatibility");claims++}
  }
  // Direct cancellation for each diagonal and several cross checks is enough as a finite control.
  for(const v of new Set([u,(u+7)%U.length,(u+33)%U.length])){
   const c0=C.map(A=>A[v]),g0=sub(o,sc(sum(c0),2)),w0=c0.map(x=>sub(x,R(1,K))),h0=sum(w0);
   ck(eq(plus(sc(sum(Cs.map((x,i)=>mul(x,c0[i]))),2),div(mul(g,g0),R(L))),
     sum([R(1,K),sc(sum(W.map((x,i)=>mul(x,w0[i]))),2),div(sc(mul(h,h0),4),R(L))])),"rank cancellation");
  }
 }
 return {cycles,n,triples:U.length,locks:claims,oldStarSums:U.length*n,rankIdentities:3*U.length};
}
function local(n){
 const C=triples(n),L=C.filter(T=>T.some(t=>t<3)),a=c2(n-2),G=arr(n,i=>arr(n,j=>L.filter(T=>T.includes(i)&&T.includes(j)).length));
 for(let i=0;i<n;i++)for(let j=0;j<n;j++)ck(G[i][j]===(i<3&&j<3?(i===j?a:0)+(n-2):i<3||j<3?n-2:(i===j?3*(n-4):0)+3),"locked local incidence Gram");
 const bracket=sub(R(3),R(3*(n-2)**2,a+3*(n-2)));
 ck(eq(bracket,R(3*(7-n),n+3))&&le(bracket,z),"external Schur negative rank-one");
 return {n,local:L.length,bracket:f(bracket)};
}
const B=constants(32,R(31,2));
ck(B.pass,"uniform threshold n32");
const cutoff={n:32,r:f(R(31,2)),...Object.fromEntries(Object.entries(B).map(([k,v])=>[k,Array.isArray(v)?f(v):v]))};
console.log(JSON.stringify({cutoff,controls:[construction([5],11),construction([3,3],13)],local:[11,13,32].map(local)}));

