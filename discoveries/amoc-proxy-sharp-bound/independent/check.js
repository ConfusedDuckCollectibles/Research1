// Independent AMOC conditional-midpoint replay. Node >=16 or a browser console.
// Only final display approximations use floating point. No outside input or dependency.
function arithmetic(d) {
  d = BigInt(d);
  const gcd = (a,b) => b === 0n ? (a < 0n ? -a : a) : gcd(b,a%b);
  function Q(a,b=0,c=1) {
    a=BigInt(a); b=BigInt(b); c=BigInt(c);
    if(c===0n) throw Error("zero denominator");
    if(c<0n){a=-a;b=-b;c=-c;}
    const g=gcd(gcd(a,b),c); return [a/g,b/g,c/g];
  }
  const plus=(a,b)=>Q(a[0]*b[2]+b[0]*a[2],a[1]*b[2]+b[1]*a[2],a[2]*b[2]);
  const neg=a=>Q(-a[0],-a[1],a[2]);
  const minus=(a,b)=>plus(a,neg(b));
  const times=(a,b)=>Q(a[0]*b[0]+d*a[1]*b[1],a[0]*b[1]+a[1]*b[0],a[2]*b[2]);
  const recip=a=>Q(a[0]*a[2],-a[1]*a[2],a[0]**2n-d*a[1]**2n);
  const div=(a,b)=>times(a,recip(b));
  function pow(a,n) {let p=Q(1);while(n--)p=times(p,a);return p;}
  const signInt=a=>(a>0n)-(a<0n);
  function sign([a,b]) {
    if(a===0n)return signInt(b); if(b===0n||a*b>0n)return signInt(a);
    return signInt(a)*signInt(a*a-d*b*b);
  }
  const eq=(a,b)=>sign(minus(a,b))===0;
  const str=a=>a.map(String).join(",");
  const num=a=>(Number(a[0])+Number(a[1])*Math.sqrt(Number(d)))/Number(a[2]);
  const sum=a=>a.reduce(plus,Q(0));
  return {Q,plus,minus,times,div,pow,sign,eq,str,num,sum};
}
function assert(b,s) {if(!b)throw Error(s);}
function splines(A, knots, curvature) {
  const {Q,plus,minus,times,div,pow,sum}=A;
  let v=Q(0),u=Q(0),mass=Q(0);
  const points=[{v,u,mass}];
  for(let j=0;j<curvature.length;j++){
    const dx=minus(knots[j+1],knots[j]), a=curvature[j];
    mass=sum([mass,times(v,dx),div(times(u,pow(dx,2)),Q(2)),
              div(times(a,pow(dx,3)),Q(6))]);
    v=sum([v,times(u,dx),div(times(a,pow(dx,2)),Q(2))]);
    u=plus(u,times(a,dx)); points.push({v,u,mass});
  }
  return points;
}
const A=arithmetic(5), {Q,plus,minus,times,div,pow,sign,eq,str,num,sum}=A;
const zero=Q(0),one=Q(1),mid=Q(1,0,2),a=Q(-1,1,4),b=Q(1,1,4);
const x=[Q(-1),Q(-1,-1,4),Q(1,-1,4),zero,a,b,one];
const bend=[-1,1,-1,-1,1,-1];
const z=splines(A,x,bend.map(i=>Q(i,0,4))); // t=2r-1
assert([z[6].v,z[6].u,z[6].mass,z[3].u].every(w=>eq(w,zero)), "clamped heat/center");
const dual=div(minus(minus(pow(b,2),pow(a,2)),mid),Q(4));
assert(eq(z[3].v,dual)&&eq(dual,Q(-2,1,16)),"sharp coefficient");
assert(eq(minus(b,a),mid)&&eq(minus(pow(b,3),pow(a,3)),mid),"dual moments");
assert(eq(z[1].v,z[5].v)&&eq(z[2].v,z[4].v)&&sign(z[1].v)<0&&sign(z[2].v)>0,"switch reflection");
const phi=Q(1,1,2);
assert(eq(plus(z[1].v,times(phi,z[2].v)),zero),"common salt jumps");
const rho=Q(1025),f=Q(65,0,1000000),g=Q(981,0,100),d=Q(3,0,10),h=Q(1800),m=Q(17000000),sv=Q(1000000);
const k=div(times(g,d),times(rho,f));
const U=times(k,pow(h,2));
const contrast=div(times(U,dual),sv), threshold=div(times(Q(2),d),contrast);
const peakMargin=minus(times(Q(16),m),U);
assert(sign(peakMargin)>0&&sign(minus(times(m,Q(2)),times(U,dual)))>0,"unique positive full-depth maxima");
const derivative0=times(rho,minus(Q(1,0,5000*90),Q(1,0,1250*900)));
const C=div(times(rho,f),times(Q(2),g)),pi=Q(22,0,7);
const grad=times(C,div(times(times(Q(4),pow(pi,3)),m),pow(h,3)));
const edge=times(C,plus(div(times(times(Q(2),pow(pi,2)),m),pow(h,2)),div(k,Q(2))));
const margins=[minus(derivative0,grad),minus(mid,edge),minus(Q(1,0,5),div(d,Q(2))),minus(div(phi,Q(5)),div(d,Q(2)))];
assert(margins.every(w=>sign(w)>0),"positive smooth and jump margins");
assert(sign(minus(contrast,Q(2111616,0,1000000)))>0&&sign(minus(Q(2111617,0,1000000),contrast))>0,"reported contrast");
assert(sign(minus(threshold,Q(284142,0,1000000)))>0&&sign(minus(Q(284143,0,1000000),threshold))>0,"reported threshold");
const B=arithmetic(2),c=B.Q(2,-1,4),r=B.times(B.Q(2),c);
const off=splines(B,[B.Q(0),c,r,B.Q(1,0,2),B.minus(B.Q(1),c),B.Q(1)],
                    [-1,1,1,-1,1].map(i=>B.Q(i)));
assert([off[5].v,off[5].u,off[5].mass].every(w=>B.eq(w,B.Q(0))),"off-center constraints");
const offAbs=B.times(B.Q(-1),off[2].v);
assert(B.eq(offAbs,B.Q(3,-2,8))&&B.sign(B.minus(offAbs,B.Q(21,0,1000)))>0&&sign(minus(Q(15,0,1000),dual))>0,"off-center larger");
function info(rho){
  const xs=Array.from({length:8},(_,j)=>Q(2*j-7,0,24));
  const v=minus(one,pow(rho,2));
  const innovations=xs.slice(1).map((x,i)=>pow(minus(x,times(rho,xs[i])),2));
  const fromAR=plus(pow(xs[0],2),div(sum(innovations),v));
  const r1=minus(one,rho);
  const closed=div(sum([times(pow(r1,2),Q(7,0,24)),times(rho,Q(7,0,144)),
                         times(times(rho,r1),Q(49,0,288))]),v);
  assert(eq(fromAR,closed),"known-covariance AR(1)");
  return str(fromAR);
}
console.log(JSON.stringify({checks:"PASS",midpointCoefficient:str(dual),offCenterCoefficient:B.str(offAbs),
  saltSwitchValues:[str(z[1].v),str(z[2].v)],contrastSv:num(contrast),twoSvD:num(threshold),
  maximumMargin:str(peakMargin),stableMargins:margins.map(str),ar1:[Q(0),Q(3,0,10),Q(3,0,5),Q(4,0,5),Q(99,0,100)].map(info)}));
