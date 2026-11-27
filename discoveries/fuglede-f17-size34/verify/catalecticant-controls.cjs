'use strict';
// Exact local replay for C-CURVES.md, section C4. Proof of exhaustiveness for
// apolar support types is analytic; these finite examples check conventions.
const p=17,F=[...Array(p).keys()],mod=x=>(x%p+p)%p;
const inv=F.map(x=>F.find(y=>mod(x*y)===1)??0);
const square=new Set(F.slice(1).map(x=>mod(x*x)));
const fac=[1];for(let i=1;i<=6;i++)fac[i]=mod(i*fac[i-1]);
const choose=(n,k)=>mod(fac[n]*inv[fac[k]]*inv[fac[n-k]]);
const canon=x=>x.map(t=>mod(t*inv[x.find(t=>t!==0)]));
const dot=(a,b)=>mod(a.reduce((s,x,i)=>s+x*b[i],0));
function rref(A){
 A=A.map(a=>a.map(mod));const piv=[];let n=A[0].length;
 for(let j=0;j<n;j++){
  let k=A.findIndex((a,i)=>i>=piv.length&&a[j]!==0);if(k<0)continue;
  const i=piv.length;piv.push(j);[A[i],A[k]]=[A[k],A[i]];
  A[i]=A[i].map(x=>mod(x*inv[A[i][j]]));
  for(let k=0;k<A.length;k++)if(k!==i){const c=A[k][j];A[k]=A[k].map((x,l)=>mod(x-c*A[i][l]));}
 }
 const nulls=Array.from({length:n},(_,j)=>j).filter(j=>!piv.includes(j)).map(j=>{
  let a=Array(n).fill(0);a[j]=1;piv.forEach((v,i)=>a[v]=mod(-A[i][j]));return a;
 });
 return {rank:piv.length,nulls};
}
const G=[0,1,2,9,16], hz=[0,16,11,0,6,1,0];
const norm=hz.map((c,i)=>mod(c*inv[choose(6,i)]));
const evalf=(x,y)=>hz.reduce((s,c,i)=>mod(s+c*x**(6-i)*y**i),0);
if(F.filter(t=>evalf(1,t)===0).join()!==G.join()||evalf(0,1))throw Error('six exact roots');
const A2=b=>rref(Array.from({length:5},(_,i)=>b.slice(i,i+3)));
const A3=b=>rref(Array.from({length:4},(_,i)=>b.slice(i,i+4)));
const h=A2(norm), h3=A3(norm);
if(h.rank!==2||h.nulls[0].join()!=='1,1,1'||h3.rank!==2||square.has(mod(-3)))throw Error('hex apolarity');

const pg=[];for(const a of F)for(const b of F)for(const c of F)
 if([a,b,c].find(t=>t!==0)===1)pg.push([a,b,c]);
const G6=(j,k)=>hz.map((c,i)=>{
 let exp=[0,0,0];exp[j]=6-i;exp[k]=i;return {c,exp};
}).filter(x=>x.c);
const Add=(...Ps)=>Ps.flat();
const deriv=(P,i)=>P.filter(a=>a.exp[i]).map(a=>{
 let exp=a.exp.slice(),c=mod(a.c*exp[i]);exp[i]--;return {exp,c};
});
const times=(P,i)=>P.map(a=>{let exp=a.exp.slice();exp[i]++;return {exp,c:a.c};});
const pure=i=>{let exp=[0,0,0];exp[i]=6;return [{exp,c:1}];};
const mons=n=>{let a=[];for(let i=0;i<=n;i++)for(let j=0;j<=n-i;j++)a.push([i,j,n-i-j]);return a};
function tRank(P,n){const coeff=new Map();for(const {exp,c} of P)coeff.set(exp.join(),mod((coeff.get(exp.join())??0)+c));
 return rref(mons(n).map(a=>mons(6-n).map(b=>{
  let e=a.map((x,i)=>x+b[i]);return mod((coeff.get(e.join())??0)*e.reduce((t,i)=>mod(t*fac[i]),1));
 }))).rank;
}
function restrict(P,n){
 let b=pg.filter(x=>dot(x,n)===0).slice(0,2);
 let h=Array(7).fill(0);
 for(const {exp,c} of P){
   let a=[c];for(let i=0;i<3;i++)for(let j=0;j<exp[i];j++){
    const z=Array(a.length+1).fill(0);
    for(let k=0;k<a.length;k++){z[k]=mod(z[k]+a[k]*b[0][i]);z[k+1]=mod(z[k+1]+a[k]*b[1][i]);}
    a=z;
   }
   a.forEach((x,i)=>h[i]=mod(h[i]+x));
 }
 return h.map((c,i)=>mod(c*inv[choose(6,i)]));
}
function check(name,P){
 const z=[],r=[];
 for(const n of pg){const b=restrict(P,n);if(b.every(x=>!x)){z.push(n);continue}
   const s=A2(b);if(s.rank!==2)continue;
   const [a,c,d]=s.nulls[0],D=mod(c*c-4*a*d);
   if(D!==0&&!square.has(D))r.push(n);
 }
 return {name,ranks:[tRank(P,1),tRank(P,2),tRank(P,3)],zeros:z,
  irredBinary:r};
}
const a=G6(0,1),b=G6(0,2),c=G6(1,2);
const pairs=check('two conjugate pairs',Add(a,b));
const jets=check('two double conjugate supports',Add(a,times(deriv(a,0),2)));
const mixed=check('two rational supports and one conjugate pair',Add(pure(0),pure(1),c));
if(pairs.ranks[2]!==4||jets.ranks[2]!==4||mixed.ranks[2]!==4)throw Error('control ranks');
const outside=[...pairs.zeros,...pairs.irredBinary].filter(n=>n[1]!==0&&n[2]!==0);
if(outside.length>2||jets.irredBinary.length>1||mixed.irredBinary.length>4)throw Error('projection controls');
console.log(JSON.stringify({binary:{coefficients:hz,normalized:norm,rank3:h3.rank,ann2:h.nulls[0]},pairs,jets,mixed}));
