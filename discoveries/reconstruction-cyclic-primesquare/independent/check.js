// Independent exact checks: number-field residues (not walk-mode computation).
const mod=(x,n)=>(x%n+n)%n;
const asc=(x,y)=>x-y;
function demand(b, label) { if(!b) throw Error(label); }
function vec(exponents,n,p) {
  const d=n/p, c=Array(n).fill(0);
  for(const j of exponents)c[mod(j,n)]++;
  return c.slice(0,n-d).map((a,i)=>a-c[n-d+i%d]).join(",");
}
function algebra(n,p,S) {
  const K=new Map, group=Array(n);
  for(let j=0;j<n;j++){
    const id=vec(S.map(t=>j*t),n,p);
    if(!K.has(id))K.set(id,[]);
    const b=K.get(id); b.push(j); group[j]=b;
  }
  return {group, key:t=>[...K.values()].map(b=>vec(b.map(j=>j*t),n,p)).join(";")};
}
function graph(n,S,t=0){
  return Array.from({length:n},(_,x)=>S.map(s=>(x+s)%n).filter(y=>!(x===0&&y===t)&&!(x===t&&y===0)));
}
function intersection(n,S,t){
  const T=new Set(S);
  return S.filter(s=>T.has(mod(s-t,n)));
}
function perms(n,p,S){
  const T=new Set(S);
  return Array.from({length:n},(_,u)=>u).filter(u=>u%p&&S.every(s=>T.has(u*s%n)));
}
function small(p,a,noFull){
  const n=p**a,d=n/p,m=Math.floor(n/2);
  let scanned=0,edgePairs=0,equalPairs=0,quarterHosts=0,halfHosts=0;
  for(let mask=1;mask<(1<<m);mask++){
    const S=[];
    for(let j=1;j<=m;j++)if((mask>>(j-1))&1){S.push(j);if(j!==n-j)S.push(n-j)}
    if(!S.some(s=>s%p))continue;
    const masses=Array(d).fill(0);
    for(const s of S)masses[s%d]++;
    if(noFull&&masses.includes(p))continue;
    scanned++;
    if(S.includes(n/4))quarterHosts++;
    if(S.includes(n/2))halfHosts++;
    const E=S.filter(t=>t<=n-t&&(noFull||!intersection(n,S,t).length));
    if(E.length<2)continue;
    const {key}=algebra(n,p,S), P=perms(n,p,S), values=E.map(key);
    for(let i=0;i<E.length;i++)for(let j=i+1;j<E.length;j++){
      edgePairs++;
      if(values[i]!==values[j])continue;
      equalPairs++;
      demand(P.some(u=>u*E[i]%n===E[j]),"unexpected alias "+JSON.stringify([n,S,E[i],E[j]]));
    }
  }
  return {n,noFull,scanned,edgePairs,equalPairs,quarterHosts,halfHosts};
}
function bareiss(A) {
  A=A.map(b=>b.slice());
  const n=A.length;
  if(!n)return 1n;
  let last=1n,sign=1n;
  for(let i=0;i<n-1;i++){
    let w=i;
    while(w<n&&!A[w][i])w++;
    if(w===n)return 0n;
    if(w!==i){[A[i],A[w]]=[A[w],A[i]]; sign=-sign;}
    const v=A[i][i];
    for(let j=i+1;j<n;j++)for(let k=i+1;k<n;k++){
      const t=v*A[j][k]-A[j][i]*A[i][k];
      demand(t%last===0n,"nonexact Bareiss quotient");
      A[j][k]=t/last;
    }
    for(let j=i+1;j<n;j++)A[j][i]=0n;
    last=v;
  }
  return sign*A[n-1][n-1];
}
function signs(){
  const n=9,S=[1,2,7,8],z=5n,t=2;
  const mat=G=>G.map((N,i)=>G.map((_,j)=>i===j?z:N.includes(j)?-1n:0n));
  const X=mat(graph(n,S)), base=bareiss(X);
  const M=(I,J)=>X.filter((_,i)=>!I.includes(i)).map(row=>row.filter((_,j)=>!J.includes(j)));
  const adj=bareiss(M([0],[t])), two=bareiss(M([0,t],[0,t])), cut=bareiss(mat(graph(n,S,t)));
  demand(cut===base+2n*adj-two,"wrong determinant sign");
  return {base:String(base), adjugate02:String(adj), endpointsRemoved:String(two), edgeRemoved:String(cut)};
}
function host17(){
  const n=289,p=17,S=[1,38,17,68,288,251,272,221];
  for(let j=0;j<n;j++)if([6,11].includes(j%p))S.push(j);
  const {key,group}=algebra(n,p,S), B=graph(n,S);
  const triFree=S.filter(t=>!intersection(n,S,t).length).sort(asc);
  const third=t=>B[0].reduce((out,u)=>out+B[u].filter(v=>B[v].includes(t)).length,0);
  demand(key(17)===key(68)&&key(1)!==key(38),"wrong exact spectral relation");
  demand(group[1].join()==="1,38,251,288"&&group[p].join()==="17,272","wrong conductor classes");
  demand(triFree.join()==="1,38,251,288"&&third(1)===469&&third(38)===177,"wrong short walks");
  // Isomorphism-invariant: the number of P-profile vertices with exactly
  // two neighbours of Q-profile. A profile lists degree, number of
  // deficient neighbours, and counts of neighbours with 0,1,2 deficient neighbours.
  function count(t) {
    const G=graph(n,S,t), d=G.map(N=>N.length), c=G.map(N=>N.filter(y=>d[y]<S.length).length);
    const k=G.map((N,x)=>[d[x],c[x],...Array.from({length:3},(_,g)=>N.filter(y=>c[y]===g).length)].join(","));
    return G.filter((N,x)=>k[x]==="42,0,42,0,0"&&N.filter(y=>k[y]==="42,0,41,1,0").length===2).length;
  }
  const notIso=[count(17),count(68)];
  demand(notIso.join()==="30,32","missing rooted profile count");
  const G=graph(n,S,17), r=S.length, w=6, card=G.map(N=>N.filter(x=>x!==w));
  const deficient=card.flatMap((N,i)=>i!==w&&N.length===r-2?[i]:[]);
  const attach=card.flatMap((N,i)=>i!==w&&N.length<r?[i]:[]);
  demand(deficient.join()==="0,17"&&attach.join()===[...G[w]].sort(asc).join(),"bad triangle card");
  return {n,degree:r,primitiveClass:group[1],lowerClass:group[p],
    equalSpectralSteps:[17,68],nonisoProfileCounts:notIso,triangleFree:triFree,
    thirdWalks:[third(1),third(38)],triangleCardDeficient:deficient};
}
const out={finite:[small(3,2,false),small(5,2,false),small(2,4,true),small(3,3,true)],rankTwo:signs(),p17:host17()};
demand(out.finite.map(r=>r.scanned).join()==="14,4092,48,4788","unexpected finite host totals");
console.log(JSON.stringify(out,null,2));