// Independent exact check of the loop geometry, certified envelope, and sparse countercontrol.
function q(n,d=1n){n=BigInt(n);d=BigInt(d);if(d<=0n)throw Error("bad denominator");
  let a=n<0n?-n:n,b=d;while(b)[a,b]=[b,a%b];return[n/a,d/a];}
const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const neg=a=>[-a[0],a[1]],sub=(a,b)=>add(a,neg(b));
const scale=(a,b)=>q(a[0]*b[0],a[1]*b[1]);
const div=(a,b)=>scale(a,b[0]<0n?q(-b[1],-b[0]):q(b[1],b[0]));
const lt=(a,b)=>a[0]*b[1]<b[0]*a[1];
const eq=(a,b)=>!lt(a,b)&&!lt(b,a);
const exact=(n,d=1)=>[q(n,d),q(n,d)];
const imul=(k,a)=>k>=0?[scale(q(k),a[0]),scale(q(k),a[1])]:[scale(q(k),a[1]),scale(q(k),a[0])];
const plus=(a,b)=>[add(a[0],b[0]),add(a[1],b[1])];
const minus=(a,b)=>plus(a,imul(-1,b));
function idiv(a,b){if(lt(a[0],q(0))||!lt(q(0),b[0]))throw Error("positive interval required");return[div(a[0],b[1]),div(a[1],b[0])];}
function must(x,msg){if(!x)throw Error(msg);}
function atanInv(k,n=32){ // even n, terms 0..n−1 give lower odd truncation
 let r=q(0);const kk=BigInt(k);
 for(let i=0;i<n;i++)r=add(r,q((-1)**i,BigInt(2*i+1)*kk**BigInt(2*i+1)));
 return [r,add(r,q(1,BigInt(2*n+1)*kk**BigInt(2*n+1)))];
}
// tan(2 atan(1/2) - atan(1/7))=1, with angles in the correct quadrant.
const pi=minus(imul(8,atanInv(2,38)),imul(4,atanInv(7,18)));
const tau=1n<<60n;
function root(k){const target=BigInt(k)*tau*tau;
 let x=1n<<BigInt(Math.ceil(target.toString(2).length/2));
 while(true){const u=(x+target/x)>>1n;if(u>=x)break;x=u;}
 must(x*x<=target&&target<(x+1n)*(x+1n),"integer root bracket");
 return [q(x,tau),q(x+1n,tau)];
}
function numer(u,t){
 const v=30*30+t*t;
 return idiv(exact(30*u,2*v),root(v+u*u)); // primitive times pi, no pi yet
}
function band(c,t){ // 0<=c<=20 on the audited benchmark
 const I=(t)=>plus(numer(20-c,t),numer(20+c,t));
 return minus(I(t),I(t+25));
}
function cross(c,t){
 // Integral of P_d(20-c,y)-P_d(20+c,y) across t..t+25, times pi.
 return minus(minus(numer(t+25,20-c),numer(t,20-c)),minus(numer(t+25,20+c),numer(t,20+c)));
}
const Km=idiv(imul(2,band(15,55)),imul(1000,pi));
const Kp=idiv(plus(band(0,40),band(0,70)),imul(1000,pi));
const KxCorner=idiv(plus(band(15,40),band(15,70)),imul(1000,pi));
const KyCorner=idiv(minus(cross(15,40),cross(15,70)),imul(1000,pi));
const wedge=idiv(KyCorner,KxCorner);
const reverse=idiv(Km,Kp),C=imul(18000000,Km),mean=idiv(imul(1000000,C),exact(260000));
const tail=idiv(exact(1),root(26));
const b=(a,l,h,d)=>must(lt(q(l,d),a[0])&&lt(a[1],q(h,d)),"narrow numerical bracket");
b(Km,8078389645n,8078389647n,10n**16n);b(Kp,1146074847n,1146074849n,10n**15n);
b(reverse,70487452,70487453,100000000);b(C,1454110136,1454110137,100000000);
b(mean,5592731293n,5592731294n,100000000);b(tail,196116135,196116136,1000000000);
must(lt(q(0),wedge[0])&&lt(wedge[1],q(7,8)),"explicit y-rake control");
must(4*40**2===6400&&3*(30**2+35**2)===6375,"geometry thresholds");
function edge(f,a,b){const u=f(...a),v=f(...b);return q((u[0]+v[0])*(b[0]-a[0])+(u[1]+v[1])*(b[1]-a[1]),2);}
function loop(f,X,Y,Z){const a=[[-X,Y],[X,Y],[X,Z],[-X,Z]];
 return a.reduce((s,u,j)=>add(s,edge(f,u,a[(j+1)%4])),q(0));}
const affine=(x,y)=>[7*x+3*y+5,11*x-2*y-7];
const grad=(x,y)=>[3*x*x*y*y,2*x*x*x*y]; // gradient of x^3 y^2, curl identically zero
function paired(f,X,Y,Z){return sub(loop(f,X,Y,Z),loop(f,X,-Z,-Y));}
must(eq(paired(affine,3,4,7),q(0))&&eq(loop(affine,3,4,7),q(8*6*3)),"affine exact control");
must(eq(paired(grad,3,4,7),q(-8*3**3*(7**2-4**2))),"non-null vertex estimator for curl-free field");
const hTerm=(n)=>sub(q(n===0?2:1,2),q(n-(n===1?1:0),12));
for(let n=0;n<4;n++)must(eq(hTerm(n),q(1,n+1)),"Hermite exact");
must(eq(sub(q(1,5),hTerm(4)),q(1,30)),"quartic Hermite control");
const printable=x=>{
 const a=x[0].toString(),b=x[1].toString();
 return Number(a.slice(0,16))/Number(b.slice(0,16))*10**(Math.max(0,a.length-16)-Math.max(0,b.length-16));
};
const mid=v=>printable(scale(add(v[0],v[1]),q(1,2)));
console.log(JSON.stringify({checks:"PASS",KminPerM:mid(Km),KmaxPerM:mid(Kp),reverseRatio:mid(reverse),
  circulationM2:mid(C),perimeterMicrometers:mid(mean),tail:mid(tail),cornerCrossRatio:mid(wedge),
  sparseGradientFalseContrast:paired(grad,3,4,7).map(String),gaussianCornerApproxMm:mid(C)/(2.4864748605243865*Math.sqrt(4450000000))*1000}));

