// Exact interval + geometric controls; standard JavaScript, no geodata.
function gcd(a,b){if(a<0n)a=-a;if(b<0n)b=-b;while(b){let m=a%b;a=b;b=m;}return a;}
function q(n,d=1){n=BigInt(n);d=BigInt(d);if(d<0){n=-n;d=-d;}if(!d)throw Error('zero');const a=gcd(n,d);return[n/a,d/a];}
const A=(x,y)=>q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]), N=x=>[-x[0],x[1]], S=(x,y)=>A(x,N(y)), M=(x,y)=>q(x[0]*y[0],x[1]*y[1]), D=(x,y)=>q(x[0]*y[1],x[1]*y[0]), P=(x,n)=>q(x[0]**BigInt(n),x[1]**BigInt(n)), lt=(x,y)=>x[0]*y[1]<y[0]*x[1], eq=(x,y)=>x[0]*y[1]===y[0]*x[1], num=x=>{const u=(x[0]<0n?-x[0]:x[0]).toString(),v=x[1].toString();return (x[0]<0n?-1:1)*Number('0.'+u.slice(0,17))/Number('0.'+v.slice(0,17))*10**(u.length-v.length);}, z=q(0);
function ok(x,s){if(!x)throw Error(s);}
const R=x=>[q(x),q(x)], ap=(u,v)=>[A(u[0],v[0]),A(u[1],v[1])], np=u=>[N(u[1]),N(u[0])], sp=(u,v)=>ap(u,np(v));
const sort=a=>a.sort((x,y)=>lt(x,y)?-1:lt(y,x)?1:0);
function mp(u,v){const a=sort(u.flatMap(x=>v.map(y=>M(x,y))));return[a[0],a[3]];}
function dp(u,v){ok(lt(z,v[0]),'positive denominator');return mp(u,[D(q(1),v[1]),D(q(1),v[0])]);}
function sqr(n){n=BigInt(n);const B=10n**17n;const a=n*B*B;let h=1n;while(h*h<=a)h*=2n;let l=0n;while(h-l>1n){const m=(h+l)>>1n;if(m*m<=a)l=m;else h=m;}return[q(l,B),q(l+1n,B)];}
function arctan(k){let a=z,b;for(let j=0;j<26;j++){b=D(q((-1)**j),M(q(2*j+1),P(q(k),2*j+1)));a=A(a,b);}b=A(a,D(q(1),M(q(53),P(q(k),53))));return[a,b];}
const pi=sp(mp(R(16),arctan(5)),mp(R(4),arctan(239)));
// Coordinates in km. Fn and K below in inverse km.
const d=30,X=20,Sx=15,Sy=15,Y=55,Z=80,L=Z-Y;
ok(4*(Y-Sy)**2>3*(d*d+(X+Sx)**2),'strict y convexity geometry');
function Fn(u,t){const b=d*d+t*t;return dp(R(d*u),mp(mp(R(2*b),pi),sqr(b+u*u)));}
function G(c,t){return sp(sp(Fn(X-c,t),Fn(-X-c,t)),sp(Fn(X-c,t+L),Fn(-X-c,t+L)));}
const lo=dp(mp(R(2),G(Sx,Y)),R(1000)),hi=dp(ap(G(0,Y-Sy),G(0,Y+Sy)),R(1000));
const ratio=dp(lo,hi), circ=mp(lo,R(18000000)), perimeter=R(260000), avgUM=mp(dp(circ,perimeter),R(1e6));
const assert=(v,a,b,den)=>ok(lt(q(a,den),v[0])&&lt(v[1],q(b,den)),[a,b,den].join('/'));
assert(lo,8078389,8078390,1e13);assert(hi,11460748,11460749,1e13);assert(ratio,704874,704875,1e6);
assert(circ,14541101,14541102,1e6);assert(avgUM,5592731,5592732,1e5);
const leak=dp(R(1),sqr(26));assert(leak,196116,196117,1e6);
ok(eq(q(Y-Sy,X+Sx),q(8,7)),'rake wedge');
function lineIntAff(x,y){const ux=C=>[2*C[0]-3*C[1]+5,7*C[0]+11*C[1]-2], a=ux(x),b=ux(y);return q((a[0]+b[0])*(y[0]-x[0])+(a[1]+b[1])*(y[1]-x[1]),2);}
function loopAff(a,b){const v=[[-X,a],[X,a],[X,b],[-X,b]];return v.reduce((c,x,i)=>A(c,lineIntAff(x,v[(i+1)%4])),z);}
ok(eq(loopAff(Y,Z),q((7-(-3))*2*X*L))&&eq(loopAff(Y,Z),loopAff(-Z,-Y)),'affine yaw cancels opposite loops');
for(let j=0;j<=4;j++){
  const at0=j===0?1:0,dp0=j===1?1:0;
  const approx=S(q(1+at0,2),q(j-dp0,12)),exact=q(1,j+1);
  ok(j<=3?eq(approx,exact):eq(S(exact,approx),q(24,720)),'Hermite remainder '+j);
}
const bracket=(a,b,d)=>[a/d,b/d];
const gauss=2.4864748605243865; const stdTargetMm=num(circ[0])/(gauss*Math.sqrt(2)*Math.hypot(40000,25000))*1e3;
console.log(JSON.stringify({checks:'PASS',kMinPerM:bracket(8078389,8078390,1e13),kMaxPerM:bracket(11460748,11460749,1e13),reverseLimit:bracket(704874,704875,1e6),minCirculationM2:bracket(14541101,14541102,1e6),minPerimeterMicrometers:bracket(5592731,5592732,1e5),outOfDiskFraction:bracket(196116,196117,1e6),cornerPrecisionMmApprox:stdTargetMm}));