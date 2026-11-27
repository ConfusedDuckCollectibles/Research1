// Node 20+, no dependencies. Exact using BigInt fractions; numeric control last.
function g(a,b){a=a<0n?-a:a; b=b<0n?-b:b; while(b){let t=a%b;a=b;b=t;}return a;}
class F{
  constructor(n,d=1){n=BigInt(n);d=BigInt(d);if(d<0n){n=-n;d=-d;}if(!d)throw Error('zero denominator');const z=g(n,d);this.n=n/z;this.d=d/z;}
  a(b){return new F(this.n*b.d+b.n*this.d,this.d*b.d);}
  s(b){return new F(this.n*b.d-b.n*this.d,this.d*b.d);}
  m(b){return new F(this.n*b.n,this.d*b.d);}
  v(b){return new F(this.n*b.d,this.d*b.n);}
  p(k){return new F(this.n**BigInt(k),this.d**BigInt(k));}
  lt(b){return this.n*b.d<b.n*this.d;}
  eq(b){return this.n===b.n&&this.d===b.d;}
  dcm(){return Number(this.n)/Number(this.d);}
  str(){return this.n+'/'+this.d;}
}
const f=(n,d=1)=>new F(n,d); const one=f(1),three=f(3),x=f(100),s=f(98);
function check(b,msg){if(!b)throw Error(msg);}
const sn=(x,z)=>one.a(three.v(x.p(2))).s(three.v(x).m(one.a(z).v(one.s(z))));
const alpha=s.p(2).s(three.m(s)).a(three).v(s.p(2).a(three.m(s)).a(three));
check(f(0).lt(alpha)&&alpha.lt(one),'reflection factor range');
let exp4lo=one,t=one;
for(let k=1;k<=9;k++){t=t.m(f(4,k));exp4lo=exp4lo.a(t);}
const exp4hi=exp4lo.a(t.m(f(4,10)).v(f(3,5)));
check(f(54).lt(exp4lo)&&exp4hi.lt(f(55)),'e^4 enclosure');
const shellLP=sn(x,alpha.v(f(54))),shellUP=sn(x,alpha.v(f(55)));
check(f(969,1000).lt(shellLP),'shell above .969 of perfectly excluding sphere');
const xlo=f(193,2),xhi=f(483,5),u=f(1,54).p(6);
check(sn(xlo,f(0)).lt(shellLP),'solid below shell at x=96.5');
check(shellUP.lt(sn(xhi,u)),'solid exceeds shell at x=96.6');
const active=one.s(s.v(x).p(3));check(active.eq(f(7351,125000)),'active fraction');
check(f(1,2).lt(f(4,5).p(3)),'half sphere radius less than .8R');
const h=f(1,54).p(10);
const halfGapUpper=f(6,100).m(h).v(one.s(h).p(2));
check(halfGapUpper.lt(f(35,10n**20n)),'half-volume gap');
check(sn(x,x.p(2).s(three.m(x)).a(three).v(x.p(2).a(three.m(x)).a(three))).eq(f(0)),'zero thickness exact');
const numeric=(x,s)=>{const z=(s*s-3*s+3)/(s*s+3*s+3)*Math.exp(-2*(x-s));return 1+3/x**2-3/x*(1+z)/(1-z);};
const ode=(x,s)=>{
  let y=s*s,w=2*s,dt=(x-s)/20000,r=s;
  const accel=(r,y)=>(1+2/r**2)*y;
  for(let j=0;j<20000;j++){
    let a=w,b=accel(r,y),c=w+dt*b/2,d=accel(r+dt/2,y+dt*a/2),
      e=w+dt*d/2,f=accel(r+dt/2,y+dt*c/2),g=w+dt*f,h=accel(r+dt,y+dt*e);
    y+=dt*(a+2*c+2*e+g)/6;
    w+=dt*(b+2*d+2*f+h)/6;
    r=s+(j+1)*dt;
  }
  return (x*w-2*y)/(x*w+y);
};
for(let [x,s] of [[2,1],[17,10],[100,98]])check(Math.abs(numeric(x,s)-ode(x,s))<1e-12,'numerical control based on direct ODE integration');
let a=96.5,b=96.6;
for(let j=0;j<40;j++){const c=(a+b)/2;if(numeric(c,0)<numeric(100,98))a=c;else b=c;}
console.log(JSON.stringify({exact:{active:active.str(),exp4:{lo:exp4lo.str(),hi:exp4hi.str()},shellS:{lo:shellLP.str(),hi:shellUP.str()},matchingSolidX:{lo:xlo.str(),hi:xhi.str()},halfGapUpper:halfGapUpper.str(),halfGapUpperTightRounded:'2.846e-19'},
numericOnly:{S:numeric(100,98),shellDirectODE:ode(100,98),matchingSolidX:(a+b)/2,matchingSolidLambda_over_lambda_min:200/(a+b)}},null,2));
