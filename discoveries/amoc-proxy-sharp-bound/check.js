// Exact Q(sqrt(5)) verifier, standard JavaScript. Rounded outputs are display only.
function gcd(a,b){if(a<0n)a=-a;if(b<0n)b=-b;while(b){const t=a%b;a=b;b=t;}return a;}
function q(a,b=0,d=1){a=BigInt(a);b=BigInt(b);d=BigInt(d);if(d===0n)throw Error("zero");if(d<0n){a=-a;b=-b;d=-d;}const g=gcd(gcd(a,b),d);return[a/g,b/g,d/g];}
const add=(u,v)=>q(u[0]*v[2]+v[0]*u[2],u[1]*v[2]+v[1]*u[2],u[2]*v[2]);
const neg=u=>q(-u[0],-u[1],u[2]), sub=(u,v)=>add(u,neg(v));
const mul=(u,v)=>q(u[0]*v[0]+5n*u[1]*v[1],u[0]*v[1]+u[1]*v[0],u[2]*v[2]);
const div=(u,v)=>mul(u,q(v[0]*v[2],-v[1]*v[2],v[0]**2n-5n*v[1]**2n));
const rat=(n,d=1)=>q(n,0,d), pow=(u,k)=>Array.from({length:k}).reduce(s=>mul(s,u),q(1));
const eq=(u,v)=>{const d=sub(u,v);return d[0]===0n&&d[1]===0n;};
const sgn=a=>a>0?1:a<0?-1:0;
function s(u){const[a,b]=u; const A=sgn(a),B=sgn(b);return A===B||!B?A:!A?B:A*sgn(a*a-5n*b*b);}
const lt=(u,v)=>s(sub(u,v))<0, F=u=>(Number(u[0])+Number(u[1])*Math.sqrt(5))/Number(u[2]), fmt=u=>u.map(String).join(',');
function ok(t,m){if(!t)throw Error(m);}
const a=q(-1,1,4), b=q(1,1,4), phi=q(1,1,2), half=rat(1,2), z=rat(0), one=rat(1), center=half;
ok(eq(sub(b,a),half)&&eq(sub(pow(b,3),pow(a,3)),half),'dual signed moments');
const roots=[q(3,-1,8),q(5,-1,8),q(3,1,8),q(5,1,8)];
const edges=[z,...roots,one];
function moment(k){let m=z;for(let i=0;i<5;i++) m=add(m,mul(rat((-1)**(i+1),k+1),sub(pow(edges[i+1],k+1),pow(edges[i],k+1))));return m;}
for(let k=0;k<=3;k++)ok(eq(moment(k),z),'piecewise exact moment '+k);
function p(x){return roots.reduce((v,t,i)=>lt(t,x)?add(v,mul(q((-1)**i),pow(sub(x,t),2))):v,neg(mul(half,pow(x,2))));}
function pd(x){return roots.reduce((v,t,i)=>lt(t,x)?add(v,mul(q(2*(-1)**i),sub(x,t))):v,neg(x));}
const C=q(-2,1,16);
ok(eq(p(z),z)&&eq(pd(z),z)&&eq(p(one),z)&&eq(pd(one),z)&&eq(pd(center),z),'edge and center derivative');
ok(eq(p(center),C)&&eq(C,mul(rat(1,4),sub(sub(pow(b,2),pow(a,2)),half))),'dual and primal exact value');
const po=p(roots[0]),pi=p(roots[1]);
ok(eq(po,q(-7,3,64))&&eq(pi,q(-11,5,64))&&eq(po,p(roots[3]))&&eq(pi,p(roots[2]))&&lt(po,z)&&lt(z,pi),'staircase pair values');
ok(eq(add(po,mul(phi,pi)),z),'matched salt staircase and heat/affine salt moments');
const rho=q(1025), ff=rat(65,1000000), g=rat(981,100), H=q(1800), M=q(17000000), D=rat(3,10), kap=div(mul(g,D),mul(rho,ff)), MM=q(1000000);
const contrast=div(mul(mul(mul(kap,pow(H,2)),C),one),MM);
const threshold=div(mul(q(2000000),D),mul(contrast,MM));
ok(lt(q(2),contrast)&&lt(contrast,rat(2112,1000)),'2 Sv contrast');
ok(lt(mul(kap,pow(H,2)),mul(q(16),M)),'actual maxima not just fixed depth pair');
const alpha=rat(1,5000),beta=rat(1,1250),tw=div(mul(rho,ff),mul(q(2),g)), piup=rat(22,7);
const deriv0=mul(rho,sub(div(alpha,q(90)),div(beta,q(900))));
const shearDerivativeUpper=mul(tw,div(mul(mul(q(4),pow(piup,3)),M),pow(H,3)));
const edgeDensityJumpUpper=mul(tw,add(div(mul(mul(q(2),pow(piup,2)),M),pow(H,2)),mul(half,kap)));
ok(lt(shearDerivativeUpper,deriv0)&&lt(edgeDensityJumpUpper,half)&&lt(mul(half,D),rat(1,5))&&lt(mul(half,D),mul(phi,rat(1,5))),'strict statically stable throughout');
ok(lt(rat(28414,100000),threshold)&&lt(threshold,rat(28415,100000)),'outward rounded threshold bracket');
function info(r){
 const r1=sub(one,r);
 return div(add(add(mul(pow(r1,2),rat(7,24)),mul(r,rat(7,144))),mul(mul(r,r1),rat(49,288))),sub(one,pow(r,2)));
}
const t=Array.from({length:8},(_,j)=>rat(2*j-7,24));
const sum=arr=>arr.reduce(add,z);
for(const r of [z,rat(3,10),rat(3,5),rat(4,5),rat(99,100)]){
 const direct=div(sub(add(add(pow(t[0],2),pow(t[7],2)),mul(add(one,pow(r,2)),sum(t.slice(1,7).map(y=>pow(y,2))))),mul(mul(q(2),r),sum(t.slice(1).map((y,j)=>mul(y,t[j]))))),sub(one,pow(r,2)));
 ok(eq(direct,info(r)),'tridiagonal AR(1) check');
}
const powers=[z,rat(3,10),rat(3,5),rat(4,5),rat(99,100)].map(r=>[F(r),F(info(r)),.1*Math.sqrt(F(info(r)))/(1.6448536269514722+.8416212335729143)]);
console.log(JSON.stringify({checks:'PASS',constant:fmt(C),contrastSv:F(contrast),twoSvD:F(threshold),smoothLayerDensityDerivative:F(deriv0),maxAdverseDerivativeUpper:F(shearDerivativeUpper),edgeJumpAdverseUpper:F(edgeDensityJumpUpper),iidAndAR1powerThresholds:powers}));