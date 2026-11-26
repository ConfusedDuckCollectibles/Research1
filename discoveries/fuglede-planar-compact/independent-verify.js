"use strict";
// Independent replay. Seed formulas: Tao Zhang, arXiv:2607.15632v1,
// equations 15, 22–24, 29–30 and the V_j definitions and layers in §4.
const seq=n=>Array.from({length:n},(_,i)=>i), mod=(x,n)=>((x%n)+n)%n;
const prod=(a,b)=>a.flatMap(x=>b.map(y=>[x,y]));
const G=prod(seq(60),seq(12)), H=prod(seq(12),seq(12));
function check(x,m){if(!x)throw new Error(m);}
const key=([x,y])=>mod(x,60)*12+mod(y,12);
const hkey=([x,y])=>mod(x,12)*12+mod(y,12);
const E=[0,5,0,7,9,2,10,4,7,8,10,10].map((x,y)=>[x,y]);
const B=[seq(12).map(x=>[x,0]),seq(12).map(x=>[x,mod(-x,12)]),
  prod(seq(3),prod(seq(2),seq(2))).map(([a,[b,c]])=>[mod(4*a+3*c,12),mod(4*a+6*b,12)])];
const V=[seq(12).map(x=>[0,x]),seq(12).map(x=>[x,x]),
  prod(seq(6),seq(2)).map(([k,j])=>[mod(4*k+j,12),mod(2*k+2*j,12)])];
const CRT=([x,y],t)=>[mod(x+12*(3*(t-x)),60),y];
const layers=(C,I)=>I.flatMap((i,t)=>C[i].map(x=>CRT(x,t)));
const N=layers(B,[0,2,2,1,1]), A=layers(V,[0,0,0,2,1]), comp=E.map(v=>CRT(v,0));
const L=prod(E,seq(5)).map(([[x,y],r])=>[mod(5*x+12*r,60),y]);
const plus=([x,y])=>[mod(2*x+15*y-23,60),mod(x+8*y,12)];
const minusLinear=([x,y])=>[mod(7*x+20*y,60),mod(x+y,12)];
const center=(x,n)=>mod(x+n/2,n)-n/2;
const minus=([x,y])=>[center(7*x+20*y-9,60),center(x+y,12)];
const S=A.map(plus), T=N.map(minus), D=comp.map(minusLinear);
for(const Z of [A,N,L,S,T])check(new Set(Z.map(key)).size===60,"seed distinctness");
for(const f of [plus,minus])check(new Set(G.map(v=>key(f(v)))).size===720,"group map");
const frequencies=new Set(L.map(key));
const mapped=G.filter(([u,v])=>frequencies.has(key([2*u+5*v,3*u+8*v])));
check(mapped.length===60,"contragredient spectrum");
function sumResidues(X,Y,M,N){
 const out=new Set();
 for(const [x,y] of X)for(const [u,v] of Y)out.add(mod(x+u,M)*N+mod(y+v,N));
 return out;
}
check(sumResidues(N,comp,60,12).size===720,"native finite tiling");
check(sumResidues(T,D,60,12).size===720,"mapped finite tiling");
// CRT monomial basis rather than cyclotomic long division.
// zeta60^j=zeta12^(5j)*zeta5^(3j), and both generating factors give degree 16.
const z12=[[1,0,0,0]];
while(z12.length<12){const [a,b,c,d]=z12.at(-1);z12.push([-d,a,b+d,c]);}
const z5=seq(5).map(j=>seq(4).map(i=>j===4?-1:+(i===j)));
const K=seq(60).map(j=>prod(z12[5*j%12],z5[3*j%5]).map(([a,b])=>a*b));
function histogram(X,u,v,q){
 const f=Array(q).fill(0),c=q/12;
 for(const [x,y] of X)f[mod(u*x+c*v*y,q)]++;
 return f;
}
function iszero(f,R){
 const w=Array(R[0].length).fill(0);
 for(let k=0;k<f.length;k++)if(f[k])
  for(let j=0;j<w.length;j++)w[j]+=f[k]*R[k][j];
 return w.every(x=>!x);
}
const Fourier60=X=>G.map(([u,v])=>iszero(histogram(X,u,v,60),K));
const zeroA=Fourier60(A),zeroS=Fourier60(S),zeroN=Fourier60(N),zeroT=Fourier60(T);
let pairs=0;
for(let i=0;i<60;i++)for(let j=i+1;j<60;j++){
 const [a,b]=L[i],[c,d]=L[j], [u,v]=mapped[i],[s,t]=mapped[j];
 check(zeroA[key([a-c,b-d])]&&zeroS[key([u-s,v-t])],"positive finite spectrum");
 pairs++;
}
const common=H.filter(([u,v])=>{
 const a=B.every(X=>iszero(histogram(X,u,v,12),z12));
 check(a===zeroN[key([5*u,v])],"common Fourier zero exact identity");
 return a;
});
function difference(X){
 const D=new Map();
 for(const [a,b] of X)for(const [x,y] of X)D.set((a-x)+","+(b-y),[a-x,b-y]);
 return [...D.values()];
}
const ad=new Set(difference(A).map(key)),sDiff=difference(S),sd=new Set(sDiff.map(key));
check(G.filter(([x,y])=>x%5!==0).every(d=>ad.has(key(d))),"all nonzero CRT layers blocked");
const vd=new Set(V.flatMap(difference).map(hkey));
check(H.every(q=>ad.has(key(CRT(q,0)))===vd.has(hkey(q))),"zero CRT difference layer");
const pgood=new Set(H.filter(q=>!vd.has(hkey(q))).map(hkey));
const ngood=new Set(common.map(hkey));
check(pgood.size===105&&ngood.size===105&&[...pgood].every(x=>ngood.has(x)),"shared Cayley graph");
function rowSolver(good){
 check(seq(12).every(y=>!good.has(hkey([0,y]))),"transversal in x");
 let visits=0;
 function find(rows){ // Fixed x order; independent of original MRV/clique branch.
  const Y=[0];let result=null;
  function recur(i){
   visits++;
   if(i===rows.length){result=rows.map((x,j)=>[x,Y[j]]);return true;}
   for(let y=0;y<12;y++){
    if(!Y.every((b,j)=>good.has(hkey([rows[i]-rows[j],y-b]))))continue;
    Y.push(y);
    if(recur(i+1))return true;
    Y.pop();
   }
   return false;
  }
  recur(1);return result;
 }
 // An 11-subset has exactly one missing x. Translate it to 11, then
 // independently translate y to make the occupied x=0 point (0,0).
 check(find(seq(11))===null,"eleven-center counterexample");
 const eleven=visits;
 let witness=null;
 for(let gap=1;gap<=6&&!witness;gap++)
  witness=find(seq(11).filter(x=>x!==11-gap));
 check(witness?.length===10,"ten-center witness");
 return {eleven,visits,witness};
}
const packing=rowSolver(pgood), dualCommon=rowSolver(ngood);
for(let i=0;i<packing.witness.length;i++)for(let j=0;j<i;j++){
 const [x,y]=packing.witness[i],[u,v]=packing.witness[j];
 check(!ad.has(key(CRT([x-u,y-v],0))),"packing witness");
}
// Directly enumerate the actual difference set in the window by integer
// shifts of S-S. This uses no per-residue witness search.
function enc(x,y){return (x+436)*197+y+98;}
const actual=new Set();
for(const [x,y] of sDiff)for(let i=-7;i<=7;i++)for(let j=-8;j<=8;j++){
 const u=x+60*i,v=y+12*j;
 if(Math.abs(u)<=436&&Math.abs(v)<=98)actual.add(enc(u,v));
}
let conflicts=0, windowDifferences=0;
for(let x=-436;x<=436;x++)for(let y=-98;y<=98;y++){
 const q=sd.has(key([x,y]));
 check(q===actual.has(enc(x,y)),"lattice window conflict equivalence");
 windowDifferences++;conflicts+=q;
}
check(windowDifferences===171981&&conflicts===146755&&437*99>10*4320,"window density gap");
// Interval arithmetic: certified Machin pi bounds; polynomial through x^28
// (the x^29 term is zero); bounded remainder and angle Lipschitz correction.
const abs=x=>x<0n?-x:x;
function gcd(a,b){while(b){[a,b]=[b,a%b];}return abs(a);}
const q=(a,b=1n)=>{const g=gcd(a,b),s=b<0n?-1n:1n;return [s*a/g,s*b/g];};
const add=(a,b)=>q(a[0]*b[1]+a[1]*b[0],a[1]*b[1]);
const neg=([a,b])=>[-a,b];
const mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]);
const lt=(a,b)=>a[0]*b[1]<b[0]*a[1];
function series(n,k){let r=q(0n);
 for(let i=0n;i<k;i++)r=add(r,q((-1n)**i,(2n*i+1n)*n**(2n*i+1n)));return r;
}
const SCALE=10n**11n, PI0=q(314159265358979n,10n**14n),PI1=q(314159265358980n,10n**14n);
const LPI=add(mul(q(16n),series(5n,14n)),mul(q(-4n),series(239n,7n)));
const UPI=add(mul(q(16n),series(5n,15n)),mul(q(-4n),series(239n,6n)));
check(lt(PI0,LPI)&&lt(UPI,PI1)&&lt(PI1,q(22n,7n)),"pi bracket");
const factorial=n=>seq(n).reduce((a,j)=>a*BigInt(j+1),1n);
const rem=q(22n**30n,7n**30n*factorial(30));
const cos=seq(60).map(j=>{
 const phi=mul(PI0,q(BigInt(Math.min(j,60-j)),30n));let s=q(0n),t=q(1n);
 const p=neg(mul(phi,phi));
 for(let i=0n;i<=14n;i++){s=add(s,t);t=mul(t,mul(p,q(1n,(2n*i+1n)*(2n*i+2n))));}
 const lo=add(s,neg(add(rem,add(PI1,neg(PI0)))));
 const [a,b]=mul(lo,q(SCALE));return a/b-(a<0n&&a%b?1n:0n);
});
let count=0, bound=null, minchars=[];
for(const [u,v] of G){
 let f=histogram(T,u,v,60);
 if(zeroT[key([u,v])])continue;
 count++;
 let a=0n;
 for(let i=0;i<60;i++)if(f[i])for(let j=0;j<60;j++)if(f[j])a+=BigInt(f[i]*f[j])*cos[mod(i-j,60)];
 check(625n*a>131n**2n*SCALE,"negative Fourier magnitude threshold");
 if(bound===null||a<bound){bound=a;minchars=[];}
 if(a===bound)minchars.push([u,v]);
}
check(count===191,"number of nonzero transformed coefficients");
const sx=T.reduce((a,[x])=>a+Math.abs(x),0),sy=T.reduce((a,[,y])=>a+Math.abs(y),0);
check(sx===720&&sy===108,"centered seed coordinate moments");
const variation=mul(q(44n,7n),add(q(BigInt(sx),60n*30n),q(BigInt(sy),12n*21n)));
check(variation[0]===1276n&&variation[1]===245n&&lt(variation,q(131n,25n)),"negative variation");
const expanded=T.flatMap(([x,y])=>prod(seq(30),seq(21)).map(([i,j])=>[x+60*i,y+12*j]));
check(expanded.length===37800&&sumResidues(expanded,D,1800,252).size===453600,"expanded direct tiling");
console.log(JSON.stringify({pairs,sharedGraphDegree:pgood.size,packing,dualCommon,windowDifferences,conflicts,windowArea:437*99,nonzero:count,
 lowerSquared:[String(bound),String(SCALE)],minchars,sx,sy,variation:variation.map(String),expandedPeriod:453600},null,2));
