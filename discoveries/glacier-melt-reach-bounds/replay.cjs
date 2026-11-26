// Reproduce exact controls. No dataset used; published source replay separately.
const abs=x=>x<0n?-x:x;
function gcd(a,b){a=abs(a);b=abs(b);while(b)[a,b]=[b,a%b];return a||1n;}
const q=(n,d=1)=>{n=BigInt(n);d=BigInt(d);if(!d)throw Error("zero denominator");if(d<0n){n=-n;d=-d;}let g=gcd(n,d);return [n/g,d/g];};
const add=(a,b)=>q(a[0]*b[1]+a[1]*b[0],a[1]*b[1]),neg=a=>q(-a[0],a[1]),sub=(a,b)=>add(a,neg(b)),mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]),sq=a=>mul(a,a);
const cmp=(a,b)=>sub(a,b)[0], absq=a=>q(abs(a[0]),a[1]), pos=a=>a[0]>0n?a:q(0);
const assert=(p,m)=>{if(!p)throw Error(m);}, eq=(a,b,m)=>assert(cmp(a,b)===0n,m+":"+a+" != "+b);
const show=a=>a[0]+"/"+a[1], sum=a=>a.reduce(add,q(0));
const hinge=(x,k,h)=>mul(k,pos(sub(x,h))), second=(a,f)=>sum(a.map((x,i)=>mul(q([1,-2,1][i]),f(q(x)))));
const A=[-80,-50,0],B=[-20,10,60],k=q(1,4);
eq(second(A,x=>x),q(20),"A differential tide");
eq(second(B,x=>x),q(20),"B differential tide");
eq(second(A,x=>hinge(x,k,q(20))),q(0),"A fixed nonlinear vertical lift");
eq(second(B,x=>hinge(x,k,q(20))),q(10),"B fixed nonlinear vertical lift");
const atoms=new Map();for (const [arr,sign] of [[A,1],[B,-1]])for(let i=0;i<arr.length;i++)
 atoms.set(arr[i],(atoms.get(arr[i])??0)+sign*[1,-2,1][i]);
const xs=[...atoms.keys()].sort((a,b)=>a-b);
let w=0,C=0, signed=0, fmin=new Map([[xs[0],q(0)]]),fmax=new Map([[xs[0],q(0)]]);
for(let i=0;i<xs.length-1;i++){
 w+=atoms.get(xs[i]);let delta=xs[i+1]-xs[i];
 C+=delta*Math.abs(w);signed+=delta*w;
 fmin.set(xs[i+1],add(fmin.get(xs[i]),w>0?mul(k,q(delta)):q(0)));
 fmax.set(xs[i+1],add(fmax.get(xs[i]),w<0?mul(k,q(delta)):q(0)));
}
assert(C===160&&signed===0,"cumulative L1 and signed area");
const effect=F=>sum(xs.map(x=>mul(q(atoms.get(x)),F.get(x))));
eq(effect(fmin),q(-20),"saturating slope lower");
eq(effect(fmax),q(20),"saturating slope upper");
// On (a,b), Peano kernel is linear; integrate absolute area by splitting at its zero.
function linearAbs(v,w,h){if(v[0]*w[0]>=0n)return mul(h,div(add(absq(v),absq(w)),q(2)));
 return div(mul(h,add(sq(v),sq(w))),mul(q(2),add(absq(v),absq(w))));}
function Pe(s){return sum(xs.map(x=>mul(q(atoms.get(x)),pos(sub(q(x),s)))));}
let pf=sum(xs.slice(1).map((x,i)=>linearAbs(Pe(q(xs[i])),Pe(q(x)),q(x-xs[i]))));
eq(pf,q(3000),"Peano integral for same tides");

const b=q(6),R=q(4),mu=q(3),m2=q(15),peak=q(110);
const g=(l,r=R)=>cmp(l,r)>0n?div(sq(sub(l,r)),mul(q(2),l)):q(0);
const ex=(ls,ps,f)=>sum(ls.map((l,i)=>mul(ps[i],f(l))));
function moments(ls,ps,u,z){eq(sum(ps),q(1),"mass");eq(ex(ls,ps,x=>x),u,"mean");eq(ex(ls,ps,sq),z,"second moment");}
const a=(c,u,z)=>div(sub(mul(c,u),z),sub(c,u));
const tangent=c=>div(mul(c,sub(mul(q(3),R),c)),add(c,R));
const pb=(u,av,c)=>div(sub(u,av),sub(c,av));
const weights=p=>[sub(q(1),p),p];
const extrema=(v,c)=>l=>div(mul(g(c),sq(sub(l,v))),sq(sub(c,v)));
const form=(u,z)=>mul(div(u,q(2)),sq(pos(sub(q(1),div(mul(R,u),z)))));
eq(g(mu),q(0),"deterministic first-moment distal rate");
eq(mul(div(mu,b),g(b)),q(1,6),"first moment upper");
let m=div(m2,mu), pl=div(mu,m);moments([q(0),m],weights(pl),mu,m2);
eq(ex([q(0),m],weights(pl),g),q(3,50),"attained lower");
eq(form(mu,m2),q(3,50),"size biased formula lower");
let aa=a(b,mu,m2), p=pb(mu,aa,b);
eq(aa,q(1),"endpoint lower atom");eq(p,q(2,5),"endpoint upper weight");
assert(cmp(aa,tangent(b))<=0n,"valid endpoint envelope");
moments([aa,b],weights(p),mu,m2);
eq(ex([aa,b],weights(p),g),q(2,15),"attained endpoint upper");
eq(mul(peak,q(3,50)),q(33,5),"dimensional lower");
eq(mul(peak,q(2,15)),q(44,3),"dimensional upper");
moments([q(0),R],weights(div(mu,R)),mu,mul(R,mu));
eq(ex([q(0),R],weights(div(mu,R)),g),q(0),"threshold equality zero");

// Interior branch, exact cubic root and tangent contact; endpoint b would be wrong.
const ci=q(29,5),ai=q(899,245),pi=q(1,2),ui=div(add(ai,ci),q(2)),z=div(add(sq(ai),sq(ci)),q(2));
eq(a(ci,ui,z),ai,"weighted moments imply left support");eq(tangent(ci),ai,"global interior tangency");
eq(sum([mul(sq(ci),ci),neg(mul(mul(q(3),R),sq(ci))),mul(sub(mul(mul(q(4),R),ui),z),ci),neg(mul(R,z))]),q(0),"exact cubic root");
assert(cmp(a(b,ui,z),tangent(b))>0n,"boundary branch fails");
let cert=extrema(ai,ci);
eq(cert(ai),q(0),"left support equality");eq(cert(ci),g(ci),"interior support equality");
moments([ai,ci],weights(pi),ui,z);
let top=ex([ai,ci],weights(pi),g);
eq(top,q(81,580),"interior optimum");
let ab=a(b,ui,z);
assert(cmp(ex([ab,b],weights(pb(ui,ab,b)),g),top)<0n,"nontrivial better than endpoint");
for(let j=0;j<=120;j++){
 let l=q(j,20);assert(cmp(cert(l),g(l))>=0n,"direct rational majorant control");
}
console.log(JSON.stringify({cumulative_abs:C,sharp_tide_contrast_cm:[-20,20],Peano_L1_cm_sq:show(pf),null_hinge_cm:[0,10],first_moment_range: ["0/1","1/6"],two_moment_range:[show(form(mu,m2)),"2/15"],interior_max:show(top)}));
