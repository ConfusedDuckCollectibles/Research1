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

let E=one,t=one;
for(let k=1;k<=25;k++){ t=t.m(f(4,k)); E=E.a(t); }
const EU=E.a(t.m(f(4,26)).v(one.s(f(4,27))));
check(f(54).lt(E)&&EU.lt(f(55)),'exp four enclosure');
const a=s.p(2).s(three.m(s)).a(three).v(s.p(2).a(three.m(s)).a(three));
const xlo=f(96544814,1000000), xhi=f(96544816,1000000), u=f(1,54).p(6);
check(sn(xlo,f(0)).lt(sn(x,a.v(E))),'tight solid x lower');
check(sn(x,a.v(EU)).lt(sn(xhi,u)),'tight solid x upper');
function sol(x){
  let A=three.v(x),B=three.v(x.p(2)),C=f(6).v(x.p(3));
  let kp=one.s(A).a(B),km=one.a(A).a(B);
  return {kp,km,dp:kp.a(B).s(C),dm:f(-1).m(km).s(B).s(C)};
}
const ms=sol(s), b=three.v(s);
const beta=b.m(ms.kp).s(ms.dp).v(ms.dm.s(b.m(ms.km)));
check(f(0).lt(beta)&&beta.lt(one),'shell l2 reflection ratio');
function Sq(x,z){
 const {kp,km,dp,dm}=sol(x);
 let d=kp.a(z.m(km)), n=dp.a(z.m(dm)), P=x.m(n).v(d);
 check(f(0).lt(d)&&f(0).lt(P),'valid positive transfer');
 return P.s(three).v(P.a(f(2)));
}
// S2 decreases in z at fixed x by the basis Wronskian; constant-q solid S2 increases in x by comparison.
const QL=Sq(x,beta.v(E)).s(Sq(xhi,f(-1).m(u)));
const QU=Sq(x,beta.v(EU)).s(Sq(xlo,f(0)));
check(f(35,1000000).lt(QL)&&QL.lt(QU)&&QU.lt(f(37,1000000)),'exact quadrupole interval');
console.log(JSON.stringify({exactDecimalsOnlyForDisplay:{xprimeLo:xlo.dcm(),xprimeHi:xhi.dcm(),
   deltaQ2Lo:QL.dcm(),deltaQ2Hi:QU.dcm()},certified:'0.000035 < shell S2 - dipole-matched sphere S2 < 0.000037'},null,2));
