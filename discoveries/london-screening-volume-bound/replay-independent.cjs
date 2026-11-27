"use strict";
const D=10n**50n;
function exact(n,d=1){
 n=BigInt(n);d=BigInt(d);if(n<0n||d<=0n)throw Error("nonpositive fraction");
 return [n*D/d,(n*D+d-1n)/d];
}
const add=(a,b)=>[a[0]+b[0],a[1]+b[1]],sub=(a,b)=>[a[0]-b[1],a[1]-b[0]];
const mul=(a,b)=>{if(a[0]<0n||b[0]<0n)throw Error("negative multiply");return [a[0]*b[0]/D,(a[1]*b[1]+D-1n)/D];};
const div=(a,b)=>{if(a[0]<0n||b[0]<=0n)throw Error("negative divide");return [a[0]*D/b[1],(a[1]*D+b[0]-1n)/b[0]];};
const mapnum=o=>o.map(x=>Number(x)/Number(D)),one=exact(1),sq=x=>mul(x,x);
function H(V,h,j){
 const z=mul(V,sq(h)),N=9;
 let term=one,sum=one;
 for(let k=1;k<=N;k++){
  term=div(mul(term,z),exact((2*k+j)*(2*k-1+j)));
  sum=add(sum,term);
 }
 const rho=div(z,exact((2*N+2+j)*(2*N+1+j)));
 return [sum[0],add(sum,div(mul(term,rho),sub(one,rho)))[1]];
}
// Solve the constant-coefficient comparison flow p'=V-p² through one cell.
// Wronskian monotonicity (positive u) holds both in initial p and V.
function local(initial,V,h){
 const C=H(V,h,0),X=mul(h,H(V,h,1));
 return div(add(mul(V,X),mul(C,initial)),add(C,mul(X,initial)));
}
function screen(P,l){return div(sub(P,exact(l+1)),add(P,exact(l)));}
function shell(l,N){
 const x=exact(100),s=exact(98),h=exact(2,N),L=exact(l*(l+1));
 let P=div(exact(l+1),s);
 for(let i=0;i<N;i++){
  const a=exact(98*N+2*i,N),b=exact(98*N+2*i+2,N);
  const lo=local([P[0],P[0]],add(one,div(L,sq(b))),h)[0];
  const hi=local([P[1],P[1]],add(one,div(L,sq(a))),h)[1];
  P=[lo,hi];
 }
 return screen(mul(x,P),l);
}
function solid(x,l){
 const A=sq(x),N=170;
 let w=one,M=exact(0),term=one;
 for(let k=1;k<=N;k++){
  term=div(mul(term,A),exact(2*k*(2*k+2*l+1)));
  w=add(w,term); M=add(M,mul(exact(2*k),term));
 }
 const rho=div(A,exact(2*(N+1)*(2*(N+1)+2*l+1))),h=sub(one,rho);
 const wtail=div(mul(term,rho),h);
 const Mtail=mul(exact(2),mul(term,add(div(mul(exact(N),rho),h),div(rho,sq(h)))));
 return screen(add(exact(l+1),div([M[0],M[1]+Mtail[1]],[w[0],w[1]+wtail[1]])),l);
}
const S1=shell(1,1024),S2=shell(2,1024),lo=exact(9654479,100000),hi=exact(9654484,100000);
if(!(solid(lo,1)[1]<S1[0]&&S1[1]<solid(hi,1)[0]))throw Error("solid match bracket");
if(!(exact(9692482007,10n**10n)[1]<S1[0]&&S1[1]<exact(9692482096,10n**10n)[0]))throw Error("explicit shell interval");
if(125000n*(100n**3n-98n**3n)!==7351n*100n**3n)throw Error("fraction check");
const g=[S2[0]-solid(hi,2)[1],S2[1]-solid(lo,2)[0]];
if(!(exact(35916,10n**9n)[1]<g[0]&&g[0]<g[1]&&g[1]<exact(35968,10n**9n)[0]))throw Error("quadrupole bound");
console.log(JSON.stringify({shell1:mapnum(S1),shell2:mapnum(S2),matchedSolidX:["96.54479","96.54484"],gap:mapnum(g),scope:"Exact directed intervals; decimals for display."},null,2));
