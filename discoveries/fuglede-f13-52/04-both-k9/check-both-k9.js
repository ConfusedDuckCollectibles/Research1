"use strict";
const B=[0,1,3,9,15,17,21,22,26,33,35,37,40,45,50,51,52,56,60,63,65,68,73,76,78,79,85,86,91,92,99,101,104,107,108,115,120,122,124,127,130,137,138,140,143,150,152,155,156,163,167,168];

const mod=x=>(x%13+13)%13,pop=n=>n?n%2+pop(n>>1):0;
const inv=Array.from({length:13},(_,i)=>i?Array.from({length:13},(_,j)=>j).find(j=>i*j%13===1):0);
const field=Array.from({length:13},(_,i)=>i);
const monomials=d=>Array.from({length:d+1},(_,a)=>Array.from({length:d+1-a},(_,b)=>[a,b])).flat();
const X=B.map(z=>z/13|0),Y=B.map(z=>z%13);
function rows(s,t){return B.flatMap((_,i)=>mod(Y[i]-s*X[i])===t?[i]:[]);}
function pointRows(d){return monomials(d).map(([a,b])=>B.map((_,i)=>mod(X[i]**a*Y[i]**b)));}
function profileShapes(T,budget=9){
 const N=T.length,v=[0,1],s0=inv[mod(T[1]-T[0])],valid=new Set,shapes=new Map;
 const n=Array(N+1).fill(0);
 for(let mask=2;mask<8192;mask+=2)if(pop(mask)===budget&&(mask>>s0&1)){
  function add(yused,dmask) {
   const i=v.length;n[i]++;
   if(i===N) {
    const k=v.join(',');if(valid.has(k))return;valid.add(k);
    const h=v.map((y,j)=>mod(y-s0*(T[j]-T[0]))),a=h.find(x=>x);
    if(!a)throw Error('affine despite line cap');
    const h2=h.map(x=>mod(x*inv[a])),tag=h2.join(',');
    if(shapes.has(tag)&&shapes.get(tag).d!==pop(dmask))throw Error('direction mismatch');
    shapes.set(tag,{h:h2,d:pop(dmask)});
    return;
   }
   for(let y=2;y<13;y++)if(!(yused>>y&1)){
    const cc=Array(13).fill(0);let m=dmask;
    let ok=true;
    for(let j=0;j<i;j++){
     const d=mod((y-v[j])*inv[mod(T[i]-T[j])]);
     if(!(mask>>d&1)||++cc[d]>=4){ok=false;break;}
     m|=1<<d;
    }
    if(ok){v.push(y);add(yused|(1<<y),m);v.pop();}
   }
  }
  add(3,1<<s0);
 }
 return {T,visits:n.slice(2),arrays:valid.size,shapes:[...shapes.values()].sort((a,b)=>(a.h.join(',')<b.h.join(',')?-1:a.h.join(',')>b.h.join(',')?1:0))};
}
function gauss(A){
 const n=A[0].length,m=A.length;
 const a=A.map((r,i)=>r.concat(Array.from({length:m},(_,j)=>+(i===j))));
 const piv=[],free=[];
 for(let j=0;j<n;j++){
  let k=piv.length;while(k<m&&!a[k][j])k++;
  if(k===m){free.push(j);continue;}
  const z=piv.length;[a[k],a[z]]=[a[z],a[k]];const mu=inv[a[z][j]];
  for(let i=j;i<n+m;i++)a[z][i]=mod(a[z][i]*mu);
  for(let k=0;k<m;k++)if(k!==z&&a[k][j]){
   const v=a[k][j];for(let i=j;i<n+m;i++)a[k][i]=mod(a[k][i]-v*a[z][i]);
  }
  piv.push(j);
 }
 const basis=free.map(j=>{
  const b=Array(n).fill(0);b[j]=1;
  for(let i=0;i<piv.length;i++)b[piv[i]]=mod(-a[i][j]);
  return b;
 });
 const U=a.map(r=>r.slice(n));
 return {piv,free,basis,U,solve(rhs){
  const z=U.map(r=>mod(r.reduce((s,x,i)=>s+x*rhs[i],0)));
  if(z.slice(piv.length).some(x=>x))return null;
  const v=Array(n).fill(0);for(let i=0;i<piv.length;i++)v[piv[i]]=z[i];
  return v;
 }};
}
function setup(){
 const I=rows(0,0),J=rows(2,0),both=[...new Set(I.concat(J))],E=pointRows(6);
 const a=profileShapes(I.map(i=>X[i])), b=profileShapes(J.map(i=>X[i]));
 const M=gauss(E.concat(both.map(i=>B.map((_,j)=>+(i===j)))));
 const anchors=[];
 for(const f of a.shapes)for(const g of b.shapes)for(let t=1;t<13;t++){
  const q={};for(let i=0;i<I.length;i++)q[I[i]]=f.h[i];
  for(let i=0;i<J.length;i++){
   const c=mod(t*g.h[i]);if(J[i] in q&&q[J[i]]!==c)throw Error('intersection matching');q[J[i]]=c;
  }
  const h=M.solve(E.map(()=>0).concat(both.map(i=>q[i])));
  if(h)anchors.push({f,g,t,h});
 }
 return {a,b,I,J,both,rank0:gauss(E).piv.length,M,anchors};
}
function quadraticRelax(M,anchors) {
 const U=pointRows(5),K=M.basis,k=K.length;
 const kv=K.flatMap((u,i)=>K.slice(i).map((v,j)=>u.map((a,z)=>mod(a*v[z]*(j===0?1:2)))));
 function dot(a,b){return mod(a.reduce((s,x,i)=>s+x*b[i],0));}
 const baseQ=U.map(row=>kv.map(v=>dot(row,v)));
 const nQ=gauss(baseQ).piv.length;
 let ranks={},cons={},aug={},survive=[],cert=[];
 for(const a of anchors){
  const hg=a.h,CC=U.map(row=>dot(row,hg.map(x=>mod(x*x))));
  const L=U.map(row=>K.map(u=>mod(2*dot(row,hg.map((v,i)=>mod(v*u[i]))))));
  const mat=baseQ.map((r,i)=>r.concat(L[i]));const H=gauss(mat);
  const rk=H.piv.length; ranks[rk]=(ranks[rk]??0)+1;

  const rhs=CC.map(v=>mod(-v));
  const x=H.solve(rhs);
  if(!x){
   let lambda=H.U.slice(rk).find(r=>dot(r,rhs)!==0);
   if(!lambda)throw Error('no left-annihilator witness');
   const t=inv[dot(lambda,rhs)];lambda=lambda.map(x=>mod(x*t));
   if(dot(lambda,rhs)!==1||mat[0].some((_,i)=>dot(lambda,mat.map(r=>r[i]))!==0))
       throw Error('incorrect inconsistency certificate');
   cert.push(lambda.map(x=>x.toString(16)).join(''));
   aug[rk+1]=(aug[rk+1]??0)+1;
   continue;
  }
  cons[rk]=(cons[rk]??0)+1;
  let r=H.basis.map(v=>v.slice(kv.length));
  survive.push({...a,H,x});
 }
 return {nQ,ranks,cons,aug,survive,cert};
}
function secondReplay(B0){
 const all=Array.from({length:13},(_,i)=>i);
 function sec(s){
  return B0.flatMap((v,i)=>(v%13-s*(v/13|0))%13===0?[i]:[]);
 }
 const I=sec(0),J=sec(2), coords=B0.map(v=>[v/13|0,v%13]);
 function allShapes(indices){
  const T=indices.map(i=>coords[i][0]),a=[0,1],R=new Map,C={seen:0,fin:0};
  const slope=Array.from({length:T.length},(_,i)=>T.map((x,j)=>i!==j?inv[mod(T[i]-x)]:0));
  function dfs(used,D){
   const p=a.length; C.seen++;
   if(p===T.length){
    C.fin++;
    const x0=T[0],m=inv[mod(T[1]-x0)];
    const b=a.map((v,i)=>mod(v-m*(T[i]-x0)));
    if(!b.some(x=>x))throw Error('line');
    const ic=inv[b.find(x=>x)],n=b.map(v=>ic*v%13);
    R.set(n.join(','),n);
    return;
   }
   for(let y=2;y<13;y++)if(!(used>>y&1)){
    let nn=D,ok=true;
    const same=Array(13).fill(0);
    for(let i=0;i<p;i++){
     const q=mod((y-a[i])*slope[p][i]);nn|=1<<q;
     if(++same[q]>=4||pop(nn)>9){ok=false;break;}
    }
    if(ok){a.push(y);dfs(used|1<<y,nn);a.pop();}
   }
  }
  dfs(3,1<<inv[mod(T[1]-T[0])]);
  return {shapes:[...R.values()],...C};
 }
 const u=allShapes(I),v=allShapes(J);
 function momentRows(n){
  const R=[];
  for(let d=0;d<=n;d++) for(let b=0;b<=d;b++)R.push(coords.map(([x,y])=>(x**(d-b)*y**b)%13));
  return R;
 }
 const out=I.concat(J.slice(1)),outSet=new Set(out);
 const free=coords.flatMap((_,i)=>outSet.has(i)?[]:[i]);
 const L=momentRows(6),Q=momentRows(5),A=L.map(r=>free.map(i=>r[i])),B=gauss(A);
 if(B.piv.length!==27||B.basis.length!==6)throw Error('independent linear reduction');
 const K=B.basis.map(a=>coords.map((_,i)=>free.includes(i)?a[free.indexOf(i)]:0));
 const dot=(r,s)=>mod(r.reduce((t,x,j)=>t+x*s[j],0));
 function squareMom(h){return Q.map(r=>dot(r,h.map(x=>x*x%13)));}
 const z=Array(52).fill(0);
 const H=[];
 // Explicit function evaluations, not coefficient products.
 for(let i=0;i<6;i++)for(let j=i;j<6;j++){
  if(i===j)H.push(squareMom(K[i]));
  else{
   const a=squareMom(K[i].map((x,k)=>mod(x+K[j][k]))),b=squareMom(K[i]),c=squareMom(K[j]);
   H.push(a.map((x,k)=>mod(x-b[k]-c[k])));
  }
 }

 function elementaryRank(rows){
  const R=rows.map(x=>x.slice());let z=0;
  for(let j=R[0].length-1;j>=0;j--){
   const i=R.findIndex((r,i)=>i>=z&&r[j]!==0);
   if(i<0)continue;
   [R[i],R[z]]=[R[z],R[i]];
   const s=R[z][j];
   for(let t=z+1;t<R.length;t++)if(R[t][j]){
    const h=R[t][j];
    R[t]=R[t].map((x,i)=>mod(s*x-h*R[z][i]));
   }
   z++;
  }
  return z;
 }
 const count={},fail={};let trials=0,pass=0;
 for(const f of u.shapes)for(const g of v.shapes)for(let t=1;t<13;t++){
  trials++;
  const hp=z.slice();
  I.forEach((j,i)=>hp[j]=f[i]);J.forEach((j,i)=>hp[j]=t*g[i]%13);
  const rest=B.solve(L.map(r=>mod(-dot(r,hp))));
  if(!rest)continue;
  pass++;
  free.forEach((j,i)=>hp[j]=rest[i]);
  if(L.some(r=>dot(r,hp)))throw Error('bad linear solution');
  const C=squareMom(hp);
  const HC=K.map((b,i)=>squareMom(hp.map((x,j)=>(x+b[j])%13)).map((q,j)=>mod(q-C[j]-H[i*6-i*(i-1)/2][j])));
  const mat=C.map((_,i)=>H.map(q=>q[i]).concat(HC.map(q=>q[i])));
  const M=elementaryRank(mat);
  const rankFull=elementaryRank(mat.map((r,i)=>r.concat(C[i])));
  count[M]=(count[M]??0)+1; fail[rankFull]=(fail[rankFull]??0)+1;
  if(rankFull!==M+1)throw Error("non-exclusion candidate");
 }
 return {shapes:[u.shapes.length,v.shapes.length],directArrays:[u.fin,v.fin],nodes:[u.seen,v.seen],unpinnedRank:B.piv.length,trials,pass,count,augmented:fail};
}
const {a,b,M,anchors,rank0}=setup();
const {nQ,ranks,aug,survive,cert}=quadraticRelax(M,anchors);
if(a.shapes.length!==6||b.shapes.length!==6||rank0!==28||M.piv.length!==46||
   nQ!==16||anchors.length!==432||survive.length||cert.length!==432)throw Error("main gate");
const replay=secondReplay(B);
if(replay.trials!==432||replay.pass!==432||
   JSON.stringify(replay.count)!==JSON.stringify(ranks)||
   JSON.stringify(replay.augmented)!==JSON.stringify(aug))throw Error("independent replay");
console.log(JSON.stringify({a,b,linearRank:rank0,pinnedRank:M.piv.length,nullity:M.free.length,anchors:anchors.length}));
console.log(JSON.stringify({pureQuadraticRank:nQ,ranks,augmented:aug,survivors:survive.length,certificates:cert.length,distinct:new Set(cert).size}));
console.log(JSON.stringify(replay));
if(typeof process!=="undefined"&&process.argv.includes("--cert"))console.log(JSON.stringify(cert));
