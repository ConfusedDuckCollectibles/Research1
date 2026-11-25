"use strict";
// Exact finishing gates for the 39-point F_13^3 spectral exclusion.
// No external data or libraries except the explicitly reproduced base below.
const p = 13;
const B = [0,1,2,3,10,11,12,13,14,16,17,20,28,30,31,56,58,59,73,79,
  86,88,89,93,94,96,103,109,123,124,126,151,152,154,156,162,165,166,168]
  .map(k => [Math.floor(k/p), k%p]);
const mod = x => (x%p+p)%p;
const inv = Array.from({length:p},(_,x) => Array.from({length:p},(_,y)=>y)
    .find(y=>x*y%p===1) || 0);
const pow = (x,k) => { let y=1;while(k--)y=mod(y*x);return y; };
const dot = (a,b) => mod(a.reduce((s,x,i)=>s+x*b[i],0));
const assert = (v,msg) => { if(!v)throw Error(msg); };
const fiber = (q,s) => B.filter(([x,y])=>(q===p?x:mod(y-q*x))===s);
const maxima = Array.from({length:p+1},(_,q)=>
  Math.max(...Array.from({length:p},(_,s)=>fiber(q,s).length)));
assert(maxima.join()==="3,7,8,3,7,3,3,3,7,3,7,3,3,7","base maxima");
const X = fiber(8,1).map(v=>v[0]).sort((a,b)=>a-b);
assert(X.join()==="0,2,4,6,7,11,12","abscissae");

// Any height graph can be sheared by the unique affine function agreeing
// with its first two heights. Do not normalize scale or use a mask quotient.
const y = [0,0], nodes=Array(X.length+1).fill(0);
function popcount(s){let t=0;while(s){s&=s-1;t++;}return t;}
function lift(k, mask){
  nodes[k]++;
  if(k===X.length)return;
  for(let h=0;h<p;h++){
    const tally=Array(p).fill(0);
    let mask2=mask,ok=true;
    for(let i=0;i<k;i++){
      const a=mod((h-y[i])*inv[mod(X[k]-X[i])]);
      mask2|=1<<a;
      if(++tally[a]>=3 || popcount(mask2)>6){ok=false;break;}
    }
    if(ok){y[k]=h;lift(k+1,mask2);}
  }
}
lift(2,1);
assert(nodes[7]===0,"unexpected height graph of ≤6 directions");

// Two-slope graph connectedness: direct flood fill on the 39 encoded points.
const connected=[];
for(const a of [0,3,5])for(const b of [0,3,5])if(a<b){
  const seen=new Set([0]), Q=[0];
  for(let i=0;i<Q.length;i++)for(let j=0;j<B.length;j++)if(!seen.has(j)){
    const [x,y]=B[Q[i]], [u,v]=B[j];
    if([a,b].some(q=>mod(y-v-q*(x-u))===0)){
      seen.add(j);Q.push(j);
    }
  }
  assert(Q.length===39,"disconnected two-slope graph");
  connected.push([a,b,Q.length]);
}

// Independent degree-7 lift finish. Work on one representative of each
// nonzero antipodal pair after translating the base symmetry center to 0.
const R=B.filter(([x,y])=>p*x+y < p*mod(-x)+mod(-y));
assert(R.length===19,"18 or 20 antipodal pairs");
const index=new Set(B.map(([x,y])=>p*x+y));
assert(B.every(([x,y])=>index.has(p*mod(-x)+mod(-y))),"not centrally symmetric");
function monomials(par,d){
  const W=[];
  for(let j=par;j<=d;j+=2)for(let i=0;i<=j;i++)
    W.push(R.map(([x,y])=>mod(pow(x,i)*pow(y,j-i))));
  return W;
}
function kernel(rows,n){
  const A=rows.map(r=>r.slice()), piv=[];
  let z=0;
  for(let j=0;j<n;j++) {
    const u=A.findIndex((r,k)=>k>=z&&r[j]);
    if(u<0)continue;
    [A[u],A[z]]=[A[z],A[u]];
    A[z]=A[z].map(v=>mod(v*inv[A[z][j]]));
    for(let k=0;k<A.length;k++)if(k!==z)
      A[k]=A[k].map((v,i)=>mod(v-A[k][j]*A[z][i]));
    piv.push(j); z++;
  }
  const vecs=[];
  for(let j=0;j<n;j++)if(!piv.includes(j)){
    const v=Array(n).fill(0); v[j]=1;
    for(let i=0;i<piv.length;i++)v[piv[i]]=mod(-A[i][j]);
    vecs.push(v);
  }
  return vecs;
}
const kEven=kernel(monomials(0,6),19);
const constraints=monomials(1,6);
for(const z of [1,13])constraints.push(R.map(([x,y])=>Number(p*x+y===z)));
const kOdd=kernel(constraints,19);
assert(kEven.length===3&&kOdd.length===5,"kernel dimensions");
const bilinear=monomials(1,5), quadratic=monomials(0,5);
let trials=0, oddNonzero=0, quadraticSurvivors=0;
for(let a=0;a<p;a++)for(let b=0;b<p;b++)for(let c=0;c<p;c++){
  if(!(a||b||c))continue;
  trials++;
  const e=Array.from({length:19},(_,i)=>mod(a*kEven[0][i]+b*kEven[1][i]+c*kEven[2][i]));
  const oddConstraints=bilinear.map(w=>kOdd.map(o=>dot(w,e.map((u,i)=>mod(u*o[i])))));
  if(kernel(oddConstraints,5).length)oddNonzero++;
  if(quadratic.every(w=>dot(w,e.map(x=>mod(x*x)))===0))quadraticSurvivors++;
}
assert(trials===2196&&oddNonzero===0&&quadraticSurvivors===0,"asymmetric degree-7 lift");
const output={maxima,q8:X,normalizedLiftNodes:nodes,connected,
evenDimension:kEven.length,oddDimension:kOdd.length,
nonzeroEvenVectors:trials,nontrivialOddKernels:oddNonzero,quadraticSurvivors};

console.log(JSON.stringify(output,null,2));
