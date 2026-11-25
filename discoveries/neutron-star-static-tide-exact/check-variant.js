// Standalone exact rational constants and oracle. Run JS in Node, browser, or V8.
const A=(x,y)=>norm(x[0]*y[1]+y[0]*x[1],x[1]*y[1]), N=x=>[-x[0],x[1]], S=(x,y)=>A(x,N(y)), M=(x,y)=>norm(x[0]*y[0],x[1]*y[1]), D=(x,y)=>norm(x[0]*y[1],x[1]*y[0]), Q=(n,d=1)=>norm(BigInt(n),BigInt(d)), P=(x,n)=>norm(x[0]**BigInt(n),x[1]**BigInt(n)), lt=(x,y)=>x[0]*y[1]<y[0]*x[1], str=x=>x.map(String).join('/'), eq=(x,y)=>x[0]*y[1]===x[1]*y[0], f=x=>Number(x[0])/Number(x[1]);
function norm(n,d){if(d<0n){n=-n;d=-d;} if(d===0n)throw 1; let a=n<0n?-n:n,b=d; while(b){let t=a%b;a=b;b=t;}return[n/a,d/a];}
function ok(c,msg){if(!c)throw Error(msg);}
const one=Q(1), sig=Q(1,10), nu=Q(1,6), a=Q(1,200), eta=Q(1,1000000), PS=Q(4,25), E0=A(one,nu), EU=Q(11,9), FL=Q(1,20), CL=Q(18,25), B=Q(3,4), CR=Q(39,50), h=Q(1,25), H=M(h,P(B,3));
const chi=S(one,M(M(Q(2),sig),EU)), G=D(M(sig,A(EU,M(Q(3),PS))),chi), Gup=Q(23,100), ep=A(EU,PS);
ok(eq(H,Q(27,1600))&&eq(chi,Q(34,45)),'geometry and compactness');
ok(lt(D(H,P(S(CL,eta),3)),FL)&&lt(A(E0,FL),EU),'all profile heights');
ok(lt(G,Gup)&&lt(D(M(ep,Gup),Q(2)),PS),'pressure first-stop');
const sound=D(M(ep,Gup),M(Q(2),nu));ok(lt(sound,one),'sound speed');
const corr=D(A(M(M(Q(2),sig),EU),M(M(Q(3),sig),ep)),chi);
const qcorr=A(D(A(M(M(Q(3),sig),A(M(Q(5),EU),M(Q(9),PS))),M(M(Q(12),sig),EU)),chi), M(Q(4),P(Gup,2)));
const qone=D(M(M(Q(3),sig),S(A(M(Q(5),EU),M(Q(9),PS)),Q(4))),chi);
ok(lt(corr,one)&&lt(qcorr,Q(6))&&lt(qone,Q(3,2))&&lt(D(M(Q(2),EU),Q(3)),S(one,PS)),'Riccati coefficient bounds');
const recip=M(Q(9),A(D(nu,Q(3)),D(M(ep,Gup),Q(2))));
ok(lt(recip,Q(2))&&lt(A(A(Q(9,2),M(Q(3),nu)),M(FL,Q(2))),Q(6)),'finite radius error');
const denomDU=A(M(Q(5),A(D(E0,Q(3)),PS)),M(Q(2,3),D(H,P(CL,3))));
ok(lt(denomDU,Q(3))&&lt(S(P(B,2),P(CL,2)),S(P(CR,2),P(B,2))),'common pressure bracketing');
const gap=M(M(Q(16,5),D(H,Q(3))),S(P(B,2),P(CL,2))), errors=A(M(Q(12),P(a,2)),M(Q(40),eta));
ok(eq(gap,Q(3969,5000000))&&lt(errors,gap),'strict finite IVT');
const X=M(Q(3,2),Gup), expU=A(A(one,X),D(P(X,2),M(Q(2),S(one,D(X,Q(3))))));
ok(lt(P(Q(6,7),2),chi)&&lt(expU,Q(10,7)),'metric ratio');
const stability=M(M(D(Q(5,3),chi),A(Q(4),Gup)),A(D(nu,Q(3)),FL));
ok(lt(stability,Q(197,200))&&eq(D(h,E0),Q(6,175)),'stability and material jump');
const Mstar=M(sig,A(A(one,M(Q(2,5),nu)),M(H,P(a,3))));
const out={checks:'PASS',compactness:str(Mstar),compactnessDecimal:f(Mstar),soundUpper:str(sound),stability:str(stability),stabilityDecimal:f(stability),IVTGap:str(gap),IVTError:str(errors),remainingWitnessGap:str(S(gap,errors)),minJump:str(D(h,E0))};
console.log(JSON.stringify(out));