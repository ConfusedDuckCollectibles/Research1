// Standard JavaScript exact verifier for the aperture and annular control. No input data.
function G(a,b){if(a<0n)a=-a;if(b<0n)b=-b;while(b){let c=a%b;a=b;b=c;}return a;}
function Q(n,d=1n){n=BigInt(n);d=BigInt(d);if(d<0n){n=-n;d=-d;}if(!d)throw Error("zero denominator");const g=G(n,d);return[n/g,d/g];}
const A=(x,y)=>Q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]), N=x=>[-x[0],x[1]], S=(x,y)=>A(x,N(y)), M=(x,y)=>Q(x[0]*y[0],x[1]*y[1]), D=(x,y)=>Q(x[0]*y[1],x[1]*y[0]), P=(x,k)=>[x[0]**BigInt(k),x[1]**BigInt(k)], lt=(x,y)=>x[0]*y[1]<y[0]*x[1], eq=(x,y)=>x[0]*y[1]===y[0]*x[1], num=x=>Number(x[0])/Number(x[1]);
function ok(v,m){if(!v)throw Error(m);}
function isqrt(a){let lo=0n,hi=1n;while(hi*hi<=a)hi*=2n;while(hi-lo>1n){const m=(lo+hi)>>1n;if(m*m<=a)lo=m;else hi=m;}return lo;}
const one=Q(1), z=Q(0), alpha=Q(1,20), B=10n**18n;
function prob(k){
 const u=BigInt(k)**2n+1n, f=isqrt(4n*B**2n/u**3n);
 return [Q(f,B),Q(f+1n,B)];
}
function testN(k,n){
 const [pLo,pHi]=prob(k);const surviv=[S(one,pHi),S(one,pLo)];
 ok(lt(Q(1,4),P(surviv[0],n-1))&&lt(P(surviv[1],n),Q(1,4)),"strict n-1/n bracket k "+k);
 const pow119=(a)=>A(alpha,S(one,P(a,119)));
 return [pow119(surviv[1]),pow119(surviv[0])];
}
const p5=testN(20,5566),p10=testN(10,703);
ok(lt(Q(79207,1000000),p5[0])&&lt(p5[1],Q(79208,1000000)),"119 power outward bracket");
const x=Q(1,3);let ln2=z;
for(let j=0;j<14;j++)ln2=A(ln2,D(M(Q(2),P(x,2*j+1)),Q(2*j+1)));
const ln2Hi=A(ln2,D(M(Q(2),P(x,29)),M(Q(29),S(one,P(x,2)))));
const Ilo=S(M(Q(4),ln2),Q(11,4)), Ihi=S(M(Q(4),ln2Hi),Q(11,4));
ok(lt(z,Ilo)&&lt(Ilo,Ihi),"nonzero annulus normalization");
const f=[4,-12,13,-6,1].map(n=>Q(n));
const integralPoly=(arr,k)=>arr.reduce((y,c,j)=>A(y,M(c,Q(2**(j+k+1)-1,j+k+1))),z);
ok(eq(integralPoly(f,0),Q(1,30))&&eq(integralPoly(f,1),Q(1,20)),"quartic moments");
const ang=Q(24,25), sqrtLo=Q(97,100);
ok(lt(P(sqrtLo,2),ang),"angular cutoff sqrt lower");
let coshLo=z;let fac=1;
for(let j=0;j<=8;j++){if(j>0)fac*=j;if(j%2===0)coshLo=A(coshLo,Q(3**j,fac));}
ok(lt(Q(10),coshLo),"arcosh(10)<3");
const k=M(D(Q(1,20),Ilo),Q(1,100)), coreUp=M(k,A(D(Q(2),sqrtLo),Q(2)));
const cuspLo=S(one,M(k,A(A(D(Q(2),sqrtLo),Q(6)),D(Q(2),M(ang,sqrtLo)))));
ok(lt(coreUp,Q(9,100))&&lt(Q(77,100),cuspLo),"finite physical separator");
function mulpoly(a,b){const c=Array(a.length+b.length-1).fill(0);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]+=a[i]*b[j];return c;}
const PP=[mulpoly([2],[1,1,1,1]),mulpoly([0,0,7],[1,1,6,6]),mulpoly([0,-2],[1,0,7,0,6])];
let pp=Array(6).fill(0);for(const row of PP)for(let j=0;j<row.length;j++)pp[j]+=row[j];
ok(JSON.stringify(pp)===JSON.stringify([2,0,9,-5,42,30])&&25-4*9*42<0,"Hernquist tracer positivity numerator");
const Aphys=Q(1,4);ok(lt(Q(13,4),D(Q(2),Aphys))&&lt(Q(291,32),D(Q(2),P(Aphys,2))),"cusp remainder physical witness");
const pf=k=>{const [l,u]=prob(k);return [num(l),num(u)];};
console.log(JSON.stringify({checks:"PASS",necessaryN80:{epsilon05:5566,epsilon10:703},epsilon05Prob:pf(20),power119Upper:"<.079208",coreUpper:num(coreUp),cuspLower:num(cuspLo),normalizerBracket:[num(Ilo),num(Ihi)]}));