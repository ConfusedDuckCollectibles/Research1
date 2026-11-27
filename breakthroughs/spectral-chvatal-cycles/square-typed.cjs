"use strict";
const assert=(b)=>{if(!b)throw Error("assert");};assert.equal=(a,b)=>{assert(a===b);};
const abs = x => x < 0n ? -x : x;
function gcd(a,b) { a=abs(a); b=abs(b); while (b) [a,b]=[b,a%b]; return a; }
function R(a,b=1) {
  a=BigInt(a); b=BigInt(b);
  assert(b);
  if(b<0n){a=-a;b=-b;}
  const d=gcd(a,b);return [a/d,b/d];
}
const add=(a,b)=>R(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const neg=a=>[-a[0],a[1]];
const sub=(a,b)=>add(a,neg(b));
const mul=(a,b)=>R(a[0]*b[0],a[1]*b[1]);
const div=(a,b)=>R(a[0]*b[1],a[1]*b[0]);
const sq=a=>mul(a,a);
const total=(...a)=>a.reduce(add,R(0));
const lt=(a,b)=>a[0]*b[1]<b[0]*a[1];
const maximum=(a,b)=>lt(a,b)?b:a;
function choose(n,k) {
  if(k<0||n<k)return 0;
  let d=1;
  for(let j=1;j<=k;j++)d=d*(n-j+1)/j;
  return d;
}
function ratio(n,r){
  const K=choose(n-1,2),a=choose(n-2,2),b=choose(n-3,2),q=choose(n-4,2),L=K-2*r;
  const c=add(R(n-2,K*a),R(12,a*(b-12))),v=R(6*(n-2),b-12);
  const omega=R(6,n-1),delta=sub(R(1,q),R(1,K));
  const kappa=sub(R(1),R(6*(n-4),a-6)),mu=sub(R(1),R(1,q));
  assert(lt(R(0),kappa));
  const types=[];
  for(let u=0;u<=3;u++)for(let t=0;t<=3-u;t++){
    const l=3-u-t,count=choose(r,u)*choose(r,t)*choose(n-2*r,l);
    if(!count)continue;
    const h=2*u+t, rho=total(mul(R(h),sub(R(1),omega)),mul(R(r),omega));
    const center=sub(R(h),mul(R(r),omega));
    const correction=total(mul([abs(center[0]),center[1]],div(delta,R(a))),mul(R(12,a*(a-12)*q),rho));
    const normError=mul(R(9,a*(a-12)),rho);
    const normMain2=R(1,a*a);
    const squareBound=total(
      mul(mul(R(13,12),normMain2),total(mul(R(2*h),sub(R(1),omega)),mul(R(r),sq(omega)))),
      mul(R(13),sq(normError)));
    const diagonal=total(
      mul(R(2),total(R(u,K*K),mul(R(r-u),sq(c)))),
      mul(R(4,L),sq(total(R(u,K),mul(R(r-u),c)))),
      mul(R(1,K*(K-1)),total(R(t),mul(R(r-t),sq(v)))),
      mul(add(R(1),R(r,q)),squareBound));
    types.push({u,t,l,count,h,correction,diagonal});
  }
  let numerator=R(0);
  let locked=0;
  const m=sub(R(1),div(omega,R(2)));
  for(const A of types)for(const B of types){
    const count=A.count*(B.count-choose(r-A.u,B.u)*choose(r-A.t,B.t)*choose(n-2*r-A.l,B.l));
    locked+=count;
    const signedSum=add(maximum(mul(m,R(A.h)),mul(omega,R(B.h))),maximum(mul(m,R(B.h)),mul(omega,R(A.h))));
    const bound=total(div(add(A.diagonal,B.diagonal),R(2)),A.correction,B.correction,div(signedSum,R(a*q)));
    numerator=add(numerator,mul(R(count),sq(bound)));
  }
  assert.equal(locked,choose(n,3)*(choose(n,3)-choose(n-3,3)));
  return div(numerator,sq(mul(kappa,mu)));
}
function countPartitions(n,min=3) {
  if(!n)return 1;
  let count=0;
  for(let j=min;j<=n;j++)count+=countPartitions(n-j,j);
  return count;
}
const data=Array.from({length:10},(_,i)=>{const n=i+22,Q=ratio(n,4);if(!lt(Q,R(1)))throw Error("r4 threshold");return [n,Q[0].toString(),Q[1].toString()];});console.log(JSON.stringify(data));
