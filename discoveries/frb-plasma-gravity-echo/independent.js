// Independent exact rational certificate of all numerical strict margins in the FRB construction.
// No event data. Decimal displays only.
function q(a,b=1n){a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero denominator");
 if(b<0n){a=-a;b=-b;}let u=a<0n?-a:a,v=b;while(v)[u,v]=[v,u%v];
 return [a/u,b/u];}
const add=(a,b)=>q(a[0]*b[1]+a[1]*b[0],a[1]*b[1]);
const neg=a=>[-a[0],a[1]], sub=(a,b)=>add(a,neg(b));
const mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]), div=(a,b)=>q(a[0]*b[1],a[1]*b[0]);
const square=a=>mul(a,a),pow=(a,k)=>q(a[0]**BigInt(k),a[1]**BigInt(k));
const lt=(a,b)=>a[0]*b[1]<b[0]*a[1], eq=(a,b)=>!lt(a,b)&&!lt(b,a);
const ok=(z,s)=>{if(!z)throw Error(s);};
const zero=q(0),one=q(1),eps=q(1,100),collar=q(1,1000);
function logRatio(x,k){ // log((1+x)/(1-x)), rigorous positive series bound
 let a=zero;
 for(let j=0;j<k;j++)a=add(a,div(mul(q(2),pow(x,2*j+1)),q(2*j+1)));
 return [a,add(a,div(mul(q(2),pow(x,2*k+1)),mul(q(2*k+1),sub(one,square(x)))))];
}
const L4=logRatio(q(3,5),36);
const kLow=div(q(11,50),L4[1]),kHi=div(q(11,50),L4[0]);
ok(lt(q(15869,100000),kLow)&&lt(kHi,q(15870,100000)),"k bracket");
ok(lt(q(4,3),L4[0]),"2 exp(4/3)<8");
const cNeg=q(441,400),cPos=q(1,100),inc=q(11,50),a=q(3,20),b=q(11,10);
ok(eq(square(add(a,one)),add(cNeg,inc))&&eq(square(sub(b,one)),cPos),"exact endpoint tracks");
const delay=div(sub(cNeg,cPos),q(2));ok(eq(delay,q(437,800)),"ray-pair delay");
ok(lt(add(cPos,inc),square(q(48,100))),"R<1.48");
const IL=add(neg(div(square(sub(one,a)),q(8))),div(square(sub(b,one)),q(2)));
const IU=div(sub(square(add(b,one)),square(add(a,one))),q(8));
ok(eq(IL,q(-273,3200))&&eq(IU,q(1235,3200)),"gap barriers exact integrals");
const w=sub(b,a),up=add(IL,add(mul(w,add(mul(q(3,8),eps),q(53,10000))),mul(q(2),collar)));
const down=sub(IU,add(mul(w,q(53,10000)),mul(q(2),collar)));
const target=[mul(q(-3,8),kHi),mul(q(-3,8),kLow)];
ok(lt(up,target[0])&&lt(target[1],down),"exact gap convex interpolation exists");
const plusR=[q(1099,1000),q(1481,1000)],minusR=[q(49,1000),q(151,1000)];
// In the b collar, q_+(r) <= (b-1)/(1-x) via e^x <= 1/(1-x) for x>=0.
const x=div(sub(cPos,square(sub(plusR[0],one))),q(15,100));
ok(lt(div(q(1,10),sub(one,x)),q(102,1000)),"plus-track left collar upper");
function ray(r,s,c,mu){
 const h=r.map(x=>add(x,s)),uLow=sub(one,div(sub(c,square(h[0])),q(15,100)));
 ok(lt(q(98,100),uLow),"annular frequency u>");
 // q target difference A>0 reduced to k*r>2*h^2*g*mu.
 ok(lt(mul(mul(mul(q(2),square(h[1])),q(1,100)),mu[1]),mul(q(15,100),r[0])),"strict cap sign throughout annulus");
 return div(mul(square(r[1]),q(16,100)),mul(mul(mul(mul(q(2),square(h[0])),q(1,100)),mu[0]),q(98,100)));
}
const s5l=q(2236067977499789n,10n**15n),s5u=q(2236067977499790n,10n**15n);
ok(lt(square(s5l),q(5))&&lt(q(5),square(s5u)),"sqrt5 bracket");
const muNeg=[div(sub(div(q(3),s5u),one),q(2)),div(sub(div(q(3),s5l),one),q(2))];
const muPos=muNeg.map(x=>add(one,x));
ok(lt(q(17,100),muNeg[0])&&lt(muNeg[1],q(18,100)),"mu− positive tight bracket");
const capUpper=add(ray(plusR,q(-1),cPos,muPos),ray(minusR,one,cNeg,muNeg));
ok(lt(capUpper,q(2000)),"angular cap upper");
const eTheta=q(1,10000000), loss=mul(mul(q(128),eTheta),q(2000));
ok(lt(loss,minusR[0])&&lt(mul(q(2),mul(eTheta,q(2000))),q(1,100)),"strict no off-axis roots, positive column");
const k_=q(4095),thirdIntensity=div(one,mul(mul(square(k_),q(1,100)),muNeg[0]));
ok(lt(thirdIntensity,q(36,1000000)),"third-to-faint intensity ratio");
const trueV=div(sub(q(3),s5l),q(2)); // actual relative voltage = (3-sqrt(5))/2
const intensityRat=div(muNeg[1],muPos[0]);ok(lt(intensityRat,q(146,1000)),"faint to bright intensity");
const toN=x=>Number(x[0])/Number(x[1]);
const tau=10e-6/toN(delay), c=299792458, au=149597870700, pc=3.0856775814913673e16,re=2.8179403262e-15;
const aAU=Math.sqrt(c*1000*pc*tau)/au,DM=tau*(2*Math.PI*800e6)**2/(2*Math.PI*c*re*pc*1e6);
console.log(JSON.stringify({checks:"PASS",delay:delay.map(String),k:toN(kLow),
 gapEnclosures:{lowUpper:toN(up),targetCoarseBounds:[-0.060,-0.05625],highLower:toN(down)},
 capUpper:toN(capUpper),angularLossBound:toN(loss),pointMassRatioVoltage:toN(trueV),
 thirdMaxRelativeVoltage:Math.sqrt(toN(thirdIntensity)),thirdMaxRelativeIntensity:toN(thirdIntensity),
 illustrative:{tauUs:tau*1e6,sizeAU:aAU,dmPcCm3:DM}}));

