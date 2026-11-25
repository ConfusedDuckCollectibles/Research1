"use strict";
// Synthetic benchmark only. Exact normalized rational numbers, no data or packages.
function gcd(a,b){a=a<0n?-a:a;while(b!==0n)[a,b]=[b,a%b];return a}
function F(n,d=1){n=BigInt(n);d=BigInt(d);if(d<0n){d=-d;n=-n}
  if(!d)throw Error("zero denominator");const g=gcd(n,d);return[n/g,d/g]}
const add=(a,b)=>F(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const mult=(a,b)=>F(a[0]*b[0],a[1]*b[1]);
const negate=a=>[-a[0],a[1]], sub=(a,b)=>add(a,negate(b));
const div=(a,b)=>mult(a,F(b[1],b[0])), square=a=>mult(a,a);
const str=a=>a[0]+"/"+a[1], set=a=>a.map(str), is=(a,n,d=1)=>str(a)===str(F(n,d));
const pos=a=>a[0]>0n;
function ok(x,label){if(!x)throw Error(label)}
function same(a,b,label){ok(JSON.stringify(a)===JSON.stringify(b),label)}
function expression(a,x){return a.reduceRight((y,k)=>add(mult(y,x),k),F(0))}
const P=[F(21,5),F(-8,25),F(1,10)], N=[F(3,10),F(-3,25),F(-1,10)];
const w=F(1,20), base=F(5,4), rcSeries=F(1,20), U0=F(1);
const positions=(x,q)=>[add(x[0],mult(x[2],q)),sub(x[1],mult(x[3],q))];
function voltages(x,q){
  const [p,n]=positions(x,q);
  return [expression(P,p),expression(N,n)];
}
const total=(x,q)=>sub(...voltages(x,q));
// Obtain a quadratic's exact jets by symmetric interpolation, not symbolic composition.
function jets(fn){
  const z=fn(F(0)), lo=fn(negate(w)), hi=fn(w);
  return [z,div(sub(hi,lo),mult(F(2),w)),
    div(sub(add(hi,lo),mult(F(2),z)),square(w))];
}
function state(x){
  const [ep,en]=[0,1].map(k=>jets(q=>voltages(x,q)[k]));
  const v=jets(q=>total(x,q)), caps=[div(U0,x[2]),div(U0,x[3])];
  const L=add(mult(caps[0],x[0]),mult(caps[1],x[1]));
  const window=[negate(w),w].map(q=>positions(x,q));
  ok(window.flat().every(a=>pos(a)&&pos(sub(U0,a))),"infeasible window");
  ok([...caps,L].every(a=>pos(a)&&pos(sub(base,a))),"infeasible pristine losses");
  const derivative=(f,q)=>add(f[1],mult(f[2],q));
  ok([negate(w),w].every(q=>!pos(add(derivative(v,q),F(1,4)))),"voltage not uniformly decreasing");
  const loss=a=>sub(U0,div(a,base));
  const g=q=>[negate(derivative(ep,q)),derivative(en,q)];
  return {v,ep,en,caps,L,window,loss:[loss(L),...caps.map(loss)],g};
}
const A=state([F(19,25),F(1,10),F(5,6),F(1)]);
const B=state([F(9,10),F(6,25),F(1),F(5,6)]);
const q0=F(0), [u,v]=A.g(q0), K=A.v[2], X=[A.ep[2],B.ep[2]], dX=sub(X[1],X[0]);
same(set(A.v),set(B.v),"unequal voltage polynomials");
same(set(A.v),["93189/25000","-7/25","61/180"],"benchmark voltage jets");
same(set(A.g(q0)),["7/50","7/50"],"benchmark split");
same(set(A.g(q0)),set(B.g(q0)),"different central split");
same(set(X),["5/36","1/5"],"curvature roots");
const vertex=(L)=>sub(L[0],div(square(L[1]),mult(F(4),L[2])));
const D=mult(F(2),sub(A.v[0],sub(vertex(P),vertex(N))));
const slopes=X.map(x=>{
  const y=sub(K,x);
  ok(is(add(div(square(u),x),div(square(v),y)),2989,12500),"failed inverse equation");
  return sub(square(div(v,y)),square(div(u,x)));
});
same(set(slopes),["-32879/62500","32879/62500"],"not both locally invertible");
const power=mult(D,K), ab=add(u,v), contrast=sub(u,v);
const disc=mult(sub(power,square(ab)),sub(power,square(contrast)));
ok(str(disc)===str(square(mult(D,dX))),"wrong factored discriminant");
const dAbs=sub(A.ep[0],B.ep[0]), dG=sub(A.g(w)[0],B.g(w)[0]);
const increment=x=>sub(voltages(x,w)[0],voltages(x,q0)[0]);
const dRel=sub(increment([F(9,10),F(6,25),F(1),F(5,6)]),increment([F(19,25),F(1,10),F(5,6),F(1)]));
same(set([dAbs,dG,dRel]),["539/25000","11/3600","11/144000"],"wrong separation");
same(set(A.loss),["119/625","1/25","1/5"],"loss A");
same(set(B.loss),["31/625","1/5","1/25"],"loss B");
// Spherical surface/mean series in z=s*tau: sinh(w)/w divided by
// 3(w*cosh(w)-sinh(w))/w^3. The G kernel has the additional factor 1/s.
const numer=[F(1),F(1,6),F(1,120)], denom=[F(1),F(1,10),F(1,280)], series=[];
for(let j=0;j<3;j++)series.push(sub(numer[j],series.reduce((a,b,k)=>add(a,mult(b,denom[j-k])),F(0))));
same(set(series),["1/1","1/15","-1/525"],"wrong kernel series");
// Uncalibrated control: arbitrary positive-real, q-constant residuals may differ
// between candidate states. Bases C=1/s, Hp=Gp-C, Hn=Gn-C are symbolic.
const U=state([F(22,25),F(2,25),F(1),F(1)]), W=state([F(23,25),F(3,25),F(1),F(1)]);
same(set(U.v),set(W.v),"uncalibrated V mismatch");
same(set(U.g(q0)),["18/125","17/125"],"uncalibrated slopes U");
same(set(W.g(q0)),["17/125","18/125"],"uncalibrated slopes W");
const delta=sub(U.g(q0)[0],W.g(q0)[0]);
function kernelCoeffs(x,q,nuisanceP,nuisanceN){
  const [up,down]=x.g(q);
  return [rcSeries,add(up,down),add(up,nuisanceP),add(down,nuisanceN)];
}
function lineCoeffs(x,a,b){
  const Z0=kernelCoeffs(x,q0,a,b);
  const d=kernelCoeffs(x,U0,a,b).map((z,j)=>sub(z,Z0[j])); // formal coefficient, not a state or measurement at 1
  return [set(Z0),set(d)];
}
const za=lineCoeffs(U,F(0),delta), zb=lineCoeffs(W,delta,F(0));
same(za,zb,"uncalibrated spectra differ at some q");
ok(pos(delta),"non-passive nuisance coefficient");
// Additional pure passive RC oracle at s=i; not an experimental fit.
// Gp-Gn=1/(1+i)-1/(1+2i)=(3/10,-1/10).
const complexGap=[mult(dG,F(3,10)),mult(dG,F(-1,10))];
same(set(complexGap),["11/12000","-11/36000"],"incorrect second-bias complex gap");
const out={calibrated:{jets:set(A.v),centralSplit:set(A.g(q0)),X:set(X),D:str(D),slopes:set(slopes),
  windowA:A.window.map(set),windowB:B.window.map(set),inventory:set([A.L,B.L]),
  lossesA:set(A.loss),lossesB:set(B.loss),absoluteReference:str(dAbs),relativeReference:str(dRel),
  secondSlope:str(dG),atI_RCComplexGap:set(complexGap),hardRadii:set([dAbs,dRel,dG].map(a=>div(a,F(2))))},
  sphericalSurfaceMeanSeries:set(series),
  uncalibrated:{jets:set(U.v),weightsU:set(U.g(q0)),weightsW:set(W.g(q0)),inventory:set([U.L,W.L]),
    windowU:U.window.map(set),windowW:W.window.map(set),positiveBackgroundDifferenceCoefficient:str(delta),
    commonConstantAndLinearCoefficientsInBasisB0C_Hp_Hn:za},
  status:"PASS"};
console.log(JSON.stringify(out,null,2));