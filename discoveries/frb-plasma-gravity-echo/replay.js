// Node 20+ or standard JS. Exact rational controls, then separate floating replay.
const check=(t,msg)=>{if(!t)throw Error(msg);};
function gcd(a,b){if(a<0n)a=-a;if(b<0n)b=-b;while(b){const c=a%b;a=b;b=c;}return a;}
const q=(n,d=1n)=>{n=BigInt(n);d=BigInt(d);const g=gcd(n,d);return d>0n?[n/g,d/g]:[-n/g,-d/g];};
const A=(x,y)=>q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]), N=x=>q(-x[0],x[1]), S=(x,y)=>A(x,N(y)), M=(x,y)=>q(x[0]*y[0],x[1]*y[1]), D=(x,y)=>q(x[0]*y[1],x[1]*y[0]), P=(x,e)=>q(x[0]**BigInt(e),x[1]**BigInt(e)), lt=(x,y)=>x[0]*y[1]<y[0]*x[1], eq=(x,y)=>!lt(x,y)&&!lt(y,x);
let log4lo=q(0);
for(let j=0;j<16;j++) log4lo=A(log4lo,q(4n,BigInt(2*j+1)*3n**BigInt(2*j+1)));
const log4hi=A(log4lo,M(q(9,8),q(4n,33n*3n**33n))); // successive summands ratio <1/9
const klo=D(q(11,50),log4hi), khi=D(q(11,50),log4lo);
check(lt(q(15,100),klo)&&lt(khi,q(16,100)),'0.15<k<0.16');
let elo=q(0), ek=q(1);
for(let j=0;j<12;j++){elo=A(elo,ek);ek=D(M(ek,q(4,3)),q(j+1));}
const ehi=A(elo,M(ek,q(39,35)));
check(lt(ehi,q(4)),'universal cutoff derivative <8');
const a=q(3,20), b=q(11,10), B=S(b,a), lo=A(N(D(P(S(q(1),a),2),q(8))),D(P(S(b,q(1)),2),q(2))), hi=D(S(P(A(b,q(1)),2),P(A(a,q(1)),2)),q(8));
check(eq(lo,q(-273,3200))&&eq(hi,q(1235,3200)),'exact envelope integrals');
const lowUpper=A(A(lo,M(B,A(q(3,800),q(53,10000)))),q(2,1000));
const highLower=S(S(hi,M(B,q(53,10000))),q(2,1000));
check(lt(lowUpper,q(-6,100))&&lt(q(-5625,100000),highLower),'strict interpolation bracket with cap costs');
check(lt(q(17,100),D(S(D(q(3),q(223607,100000)),q(1)),q(2)))&&lt(D(S(D(q(3),q(223606,100000)),q(1)),q(2)),q(18,100))&&lt(P(q(223606,100000),2),q(5))&&lt(q(5),P(q(223607,100000),2)),'.17 < mu^- < .18; plus adds one');
const data=[[q(49,1000),q(151,1000),q(1049,1000),q(1151,1000),q(17,100),q(18,100),q(441,400)], [q(1099,1000),q(1481,1000),q(99,1000),q(481,1000),q(117,100),q(118,100),q(1,100)]];
for(const [r0,r1,h0,h1,mu0,mu1,c] of data){
  check(lt(D(S(c,P(h0,2)),q(15,100)),q(2,100)),'u_i>0.98');
  check(lt(M(q(1,100),mu1),D(M(q(15,100),r0),M(q(2),P(h1,2)))),'positive matching A');
  check(lt(M(D(P(r1,2),q(98,100)),D(q(16,100),M(M(q(2,100),P(h0,2)),mu0))),q(2000)),'A < 2000 on disjoint annuli');
}
check(lt(M(q(4*32*2000),q(1,10000000)),q(49,1000)),'no off-axis stationary rays');
check(eq(D(S(q(441,400),q(1,100)),q(2)),q(437,800)),'matched phase/group lag');
check(lt(D(q(1),M(q(4095**2),q(17,10000))),P(q(6,1000),2)),'max ray vs weaker <.006 voltage');
console.log('EXACT rational inequalities PASS');
// Floating implementation of the globally certified profile. Grid samples are not the proof.
const log=Math.log, sqrt=Math.sqrt, exp=Math.exp, k=.22/log(4), ca=1.1025, cb=.01, GG=.01, pm=(3/sqrt(5)+1)/2, mm=pm-1, K=4096, de=.001;
const smf=t=>t<=0||t>=1?0:exp(-1/(t*(1-t)));
function quad(f,a,b,N=32){let z=f(a)+f(b);for(let i=1;i<N;i++)z+=(i&1?4:2)*f(a+(b-a)*i/N);return z*(b-a)/(3*N);}
const zn=quad(smf,0,1,128);
const ss=t=>t<=0?0:t>=1?1:quad(smf,0,t,64)/zn, mix=(x,y,t)=>x+(y-x)*ss(t);
const rP=1+sqrt(.23);
const qp=r=>(r-1)*exp((cb-(r-1)**2)/k), qm=r=>(r+1)*exp((ca-(r+1)**2)/k);
function gap(r,j){
  const up=(r+1)/4, dl=(5*(r-1)+3*sqrt((r-1)**2+1e-4))/8, v=dl+(up-dl)*j;
  if(r<.151)return mix(qm(r),v,(r-.15)/de);
  if(r>1.099)return mix(v,qp(r),(r-1.099)/de);
  return v;
}
const integGap=j=>quad(r=>gap(r,j),.15,.151,32)+quad(r=>gap(r,j),.151,1.099,2048)+quad(r=>gap(r,j),1.099,1.1,32);
const gl=integGap(.01),gh=integGap(.99), w=(-3*k/8-gl)/(gh-gl);
check(0<w&&w<1,'interior gap integral match');
function rder(r){
  if(r<2/K)return K*r;
  if(r<3/K)return mix(K*r,2,K*r-2);
  if(r<.04)return 2;
  if(r<.045)return mix(2,qm(r),(r-.04)/.005);
  if(r<.15)return qm(r);
  if(r<1.1)return gap(r,.01+.98*w);
  if(r<1.49)return qp(r);
  if(r<1.50)return mix(qp(r),0,(r-1.49)/.01);
  return 0;
}
const knots=[0,2/K,3/K,.04,.045,.15,.151,1.099,1.1,1.49,1.50];
const caps=(r,l,a,b,h)=>ss((r-l)/(a-l))*ss((h-r)/(h-b));
function Ap(r){if(r<=1.099||r>=1.481)return 0;const h=r-1,u=exp((h*h-cb)/k);return caps(r,1.099,1.1,rP,1.481)*r*r/u*(k/(2*h*h*GG*pm)-1/r);}
function Am(r){if(r<=.049||r>=.151)return 0;const h=r+1,u=exp((h*h-ca)/k);return caps(r,.049,.05,.15,.151)*r*r/u*(k/(2*h*h*GG*mm)-1/r);}
function n(r){
  let v=1;
  for(let j=0;j<knots.length-1;j++){const l=knots[j],h=Math.min(r,knots[j+1]);if(l>=h)break;v+=quad(rder,l,h,j===6?2048:64);}
  return v*(1-ss(r-2));
}
let dm=0, df=0, dr=0, eps=1e-7;
for(let j=0;j<=32;j++){
  const u=1+3*j/32, p=1+sqrt(cb+k*log(u)), m=sqrt(ca+k*log(u))-1, c=1/(u*K-1);
  check(Math.abs(p-u*rder(p)-1)<1e-11&&Math.abs(m-u*rder(m)+1)<1e-11&&Math.abs(c-u*rder(c)+1)<1e-11,'three specified roots');
  const nm=n(m),np=n(p), d=(m+1)**2/2-(p-1)**2/2;
  dm=Math.max(dm,Math.abs(nm-np));df=Math.max(df,Math.abs(d-u*(nm-np)-437/800),Math.abs(d+u*(nm-np)-437/800));
  const mp=1/((2*(p-1)**2/k)*(p+u*Ap(p))/p**2), me=1/((2*(m+1)**2/k)*(m+u*Am(m))/m**2);
  dr=Math.max(dr,Math.abs(mp-GG*pm),Math.abs(me-GG*mm));
  check(c<m&&m<p&&Ap(m)===0&&Am(p)===0,'ordering/support');
}
check(dm<2e-8&&df<8e-8&&dr<1e-12,'floating quartic-integration/control error');
console.log(JSON.stringify({approximateOnly:true,gapMixWeight:w,maxBranchDmError:dm,maxPhaseOrGroupDiffError:df,maxRayIntensityError:dr,pointMassIntensityRatio:mm/pm,thirdVsWeakerVoltage:1/((K-1)*sqrt(GG*mm)),astrometricSeparationHigh:1.15,astrometricSeparationLow:1.15+sqrt(.23)}));
