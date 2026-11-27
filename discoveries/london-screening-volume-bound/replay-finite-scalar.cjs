'use strict';
// Independent constant-layer transfer for the printed decimal example.
// This part is floating-point. The second part is an exact witness at b=1/25.
function posRoot(l,q,a){
 let term=1,w=1,d=0;
 for(let j=1;j<=32;j++){
   term*=q*a*a/(2*j*(2*l+2*j+1));
   w+=term; d+=2*j*term/a;
 }
 return [a**(l+1)*w, a**l*((l+1)*w+a*d)];
}
function funcs(l,q,r,left,sgn){
 const k=Math.sqrt(q),x=k*r;
 let term=1,P=1,d=0; 
 for(let j=1;j<=l;j++){
   term*=-sgn*(l+j)*(l-j+1)/(2*j*x);
   P+=term; d-=j*term/x;
 }
 const e=Math.exp(sgn*k*(r-left));
 return [e*P,k*e*(sgn*P+d)];
}
function propagate(l,q,a,b,[u,d]){
 const f=funcs(l,q,a,a,+1),g=funcs(l,q,a,a,-1);
 const det=f[0]*g[1]-f[1]*g[0];
 if(Math.abs(det+2*Math.sqrt(q))>1e-10)throw Error('basis Wronskian');
 const A=(u*g[1]-d*g[0])/det,B=(d*f[0]-u*f[1])/det;
 const F=funcs(l,q,b,a,+1),G=funcs(l,q,b,a,-1);
 return [A*F[0]+B*G[0],A*F[1]+B*G[1]];
}
const r=[0,.5,.625,.75,.875,1];
const base=[0,2,2,2,2], fill=[.01,1.967105046681599,2.0209240832497946,1.9957677210763085,2];
function data(q){
 let p=[1,2,3].map(l=>{
   let y=posRoot(l,q[0],r[1]);
   for(let i=1;i<q.length;i++) y=propagate(l,q[i],r[i],r[i+1],y);
   return y[1]/y[0];
 });
 return {p,S:p.map((p,i)=>(p-i-2)/(p+i+1)),mass:q.reduce((s,c,i)=>s+c*(r[i+1]**3-r[i]**3)/3,0)};
}
const a=data(base),b=data(fill);
const delta={slopes:b.p.map((x,i)=>x-a.p[i]),response:b.S.map((x,i)=>x-a.S[i]),mass:b.mass-a.mass};
if(Math.max(...delta.slopes.slice(0,2).map(Math.abs),Math.abs(delta.mass))>2e-12)throw Error('decimal example discrepancy');
if(!(delta.slopes[2]<-1.1034e-5&&delta.slopes[2]>-1.1035e-5))throw Error('held-out channel');
const rows=[3,5,7].map(k=>[4n,5n,6n].map(i=>(i+1n)**BigInt(k)-i**BigInt(k)));
let det=0n, perm=0n;
for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let k=0;k<3;k++)if(i!==j&&j!==k&&i!==k){
 const v=rows[0][i]*rows[1][j]*rows[2][k];
 perm+=v;det+=((i>j)+(i>k)+(j>k))%2===0?v:-v;
}
// b=1/25. Relative row errors <= b/5 and b/7; each Leibniz product decreases
// by relative at most b*(1/5+1/7)=12/875. All column factors are positive.
const margin=875n*det-12n*perm;
if(margin<=0n)throw Error('not a positive-stiffness Jacobian certificate');
console.log(JSON.stringify({floating:{baseline:a,filled:b,delta},exact:{shellStiffness:'1/25',scaledDetZero:String(det),scaledPermanent:String(perm),errorFactor:'12/875',positiveJacobianMarginTimes875:String(margin)}},null,2));
