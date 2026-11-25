
// Solver-free exact controls for ordinary vertex reconstruction of F-e.
const mod=(x,n)=>(x%n+n)%n;
const same=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
function insist(x,msg){if(!x)throw Error(msg);}
function sorted(x){return [...x].sort((a,b)=>a-b);}
function walkColumns(n,S,kmax){
  let a=Array(n).fill(0n);a[0]=1n;
  const W=[a];
  for(let k=1;k<=kmax;k++){
    a=a.map((_,i)=>S.reduce((s,d)=>s+a[mod(i-d,n)],0n));
    W.push(a);
  }
  return W;
}
function common(n,S,t){
  const set=new Set(S);
  return S.filter(x=>set.has(mod(x-t,n)));
}
function mulS(n,S,u){
  return sorted(S.map(a=>mod(u*a,n)));
}
function fieldControl(p){
  const n=p*p,m=(n-1)/2;
  let connected=0, full=0, partial=0, blowup=0, triFree=0, equalPairs=0;
  for(let mask=1;mask<(1<<m);mask++){
    const S=[];
    for(let x=1;x<=m;x++)if((mask>>(x-1))&1)S.push(x,n-x);
    if(!S.some(x=>x%p!==0))continue;
    S.sort((a,b)=>a-b);connected++;
    const masses=Array(p).fill(0);
    for(const s of S)masses[s%p]++;
    if(masses.includes(p))full++;
    if(masses.slice(1).some(v=>v>0&&v<p))partial++;
    else blowup++;
    const T=S.filter(x=>x<n-x&&common(n,S,x).length===0);
    triFree+=T.length;
    if(!T.length)continue;
    // Enough moments: lambda_j=lambda_-j, so ≤m+1 distinct eigenvalues.
    const W=walkColumns(n,S,m);
    for(let i=0;i<T.length;i++)for(let j=i+1;j<T.length;j++){
      const [s,t]=[T[i],T[j]];
      if(!W.every(a=>a[s]===a[t]))continue;
      equalPairs++;
      if(!S.some(x=>x%p===0)&& !masses.slice(1).some(v=>v>0&&v<p)) {
        // Blow-up: quotient multiplier, acting multiplicatively on the
        // whole original cyclic group, is enough for two unit edges.
      }
      const ok=Array.from({length:n-1},(_,i)=>i+1)
        .some(u=>u%p!==0&&[t,n-t].includes(mod(u*s,n))&&same(mulS(n,S,u),S));
      insist(ok,'unexpected non-multiplier triangle-free pair');
    }
  }
  return {p,connected,full,partial,blowup,triFree,equalPairs};
}
function jointRefinement(n,S,removed){
  const A=removed.map(t=>Array.from({length:n},(_,i)=>
    S.map(a=>mod(i+a,n)).filter(j=>!(i===0&&j===t)&&!(i===t&&j===0))));
  let C=A.map(X=>X.map(N=>N.length));
  for(let step=1;step<=2*n;step++){
    const codes=new Map();
    const D=C.map((c,j)=>c.map((v,i)=>{
      const key=v+':'+sorted(A[j][i].map(x=>c[x])).join(',');
      if(!codes.has(key))codes.set(key,codes.size);
      return codes.get(key);
    }));
    if(!same(sorted(D[0]),sorted(D[1]))){
        const counts=D.map(e=>{const a=Array(codes.size).fill(0);for(const t of e)a[t]++;return a;});
        const i=counts[0].findIndex((w,k)=>w!==counts[1][k]);
        return {step,classes:codes.size,firstMismatch:{color:i,multiplicities:[counts[0][i],counts[1][i]]}};
    }
    if(same(C,D))return null;
    C=D;
  }
  throw Error('refinement did not settle');
}
function mixedControl(){
  const p=17,n=p*p;
  const S=sorted(new Set([1,38,17,68].flatMap(x=>[x,n-x]).concat(
    Array.from({length:n},(_,k)=>k).filter(k=>k%p===6||k%p===p-6))));
  const W=walkColumns(n,S,n-1);
  const mult=38;
  insist(mult*mult%n===n-1,'m^2');
  insist(W.every(row=>row[17]===row[68]),'no exact cospectrality');
  const ref=jointRefinement(n,S,[17,68]);
  insist(ref&&ref.step===3,'not the announced distinction');
  const tFree=S.filter(t=>common(n,S,t).length===0);
  insist(same(tFree,[1,38,251,288]),'triangle flags');
  const distinguishes=W.findIndex(row=>row[1]!==row[38]);
  insist(distinguishes>=0,'unit step distinctions');
  return {p,n,valency:S.length,fullFibres:[6,11],triangularEdges:[17,68],
    commonNeighbours:[17,68].map(t=>common(n,S,t).length),
    triangleFreeSteps:tFree,firstDifferentWalkLengthForUnitSteps:distinguishes,
    unitWalkCounts:[String(W[distinguishes][1]),String(W[distinguishes][38])],
    equalWalkMomentsThrough:n-1,refinement:ref};
}
const report={finiteControls:[fieldControl(3),fieldControl(5)],explicitMixedControl:mixedControl()};
console.log(JSON.stringify(report,null,2));
