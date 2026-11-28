'use strict';
// Independent exact checks of minimum-degree normalization. Node, no imports.
// Modular checks are identity controls, not substitutes for real PSD.
const P=65521, mod=x=>(x%P+P)%P, inv=new Int32Array(P);
inv[1]=1; for(let i=2;i<P;i++)inv[i]=mod(-Math.floor(P/i)*inv[P%i]);
const check=(v,label)=>{if(!v)throw new Error(label);};
const c2=n=>n*(n-1)/2,c3=n=>n*(n-1)*(n-2)/6;
function gcd(a,b){while(b)[a,b]=[b,a%b];return a<0n?-a:a;}
const Q=(a,b=1)=>{a=BigInt(a);b=BigInt(b); const g=gcd(a,b)*(b<0n?-1n:1n);check(b!==0n,'zero rational denominator');return [a/g,b/g];};
const add=(a,b)=>Q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]);
const sub=(a,b)=>add(a,[-b[0],b[1]]);
const mul=(a,b)=>Q(a[0]*b[0],a[1]*b[1]);
const div=(a,b)=>Q(a[0]*b[1],a[1]*b[0]);
const sq=a=>mul(a,a), nsum=(...a)=>a.reduce(add,Q(0));
const g=(a,b)=>a[0]*b[1]>b[0]*a[1];
const asFloat=a=>Number(a[0])/Number(a[1]);

function bound(n,D,m,alpha,d0,e,sigma) {
  const K=Q(c2(n-1)),a=Q(c2(n-2)),b=Q(c2(n-3)),q=Q(c2(n-4)),N=Q(c3(n));
  const M=sub(K,Q(d0)), by6=Q(6*D), a0=sub(a,by6),b0=sub(b,by6),L=sub(M,alpha),tail=sub(M,Q(e));
  const omega=Q(6,n-1), delta=div(Q(3*(n-3)),mul(K,q)),H=add(Q(3*D),mul(m,omega));
  check([L,tail,a0,b0].every(x=>g(x,Q(0)))&&!g(a,M),'strict denominators and M>=a');
  const c=add(div(sub(M,a),mul(M,a)),div(by6,mul(a,b0))),v=div(Q(6*(n-2)),b0);
  const R=nsum(div(Q(3*D),sq(M)),mul(alpha,sq(c)),div(sq(add(div(Q(3*D),M),mul(alpha,c))),L),
    div(add(Q(3*e),mul(sigma,sq(v))),mul(M,tail)));
  const A=add(div(add(div(Q(3*D),K),mul(H,delta)),a),div(mul(by6,H),mul(mul(a,a0),q)));
  const U=add(div(mul(Q(13,12),add(Q(6*D),mul(m,sq(omega)))),sq(a)),
    mul(Q(13*D),sq(div(mul(Q(6),H),mul(a,a0)))));
  const offset=div(Q(d0),mul(K,M));
  const E=nsum(R,mul(Q(2),A),mul(add(Q(1),div(m,q)),U),offset);
  const mu=sub(Q(1),div(Q(1),q)), kappa=sub(Q(1),div(Q(6*(n-4)),sub(a,Q(3*D))));
  const count=mul(N,sub(N,Q(c3(n-3))));
  const ratio=div(mul(sq(E),count),sq(mul(mu,kappa)));
  check(g(mul(mu,kappa),Q(0)) && g(Q(1),ratio),'exact squared sufficient criterion');
  return {squaredRatio:Number(asFloat(ratio).toFixed(10)), offset:offset.map(String)};
}
const uR3=div(Q(3,4),sq(Q(11,25)));
const uR4=nsum(Q(25,2),div(sq(add(div(Q(3,8),Q(11,25)),Q(5,4))),Q(3,10)),div(Q(225,8),sq(Q(11,25))));
const uA3=div(Q(3,8),sq(Q(9,20)));
const uA4=add(div(Q(10),Q(9,20)),div(Q(1,2),mul(sq(Q(9,20)),Q(2,5))));
const uU3=div(mul(mul(Q(11,10),Q(13,12)),Q(3,4)),sq(Q(9,20)));
const uU4=div(mul(mul(Q(11,10),Q(13,12)),Q(49,24)),sq(Q(9,20)));
const uU5=mul(mul(Q(11,10),Q(13,8)),sq(div(Q(4),mul(Q(9,20),Q(2,5)))));
const uOffset=div(Q(1,8),mul(Q(11,25),Q(9,20)));
const coarse=[[uR3,4],[uR4,185],[uA3,2],[uA4,30],[uU3,5],[uU4,13],[uU5,900],[uOffset,1]];
for(const [x,y] of coarse)check(g(Q(y),x),'coefficient round-up');
check(g(Q(76483,90000),Q(1669,2000)),'strict worst-case slack');
function solve(A,b){
  // Plain modular elimination, not an asymptotic approximation.
  const n=b.length,H=A.map((a,i)=>Int32Array.from([...a,b[i]],mod));
  for(let j=0;j<n;j++){
    let p=j;while(p<n&&H[p][j]===0)p++;
    check(p<n,'modular Gram invertibility');
    if(p!==j)[H[p],H[j]]=[H[j],H[p]];
    const v=inv[H[j][j]];
    for(let k=j;k<=n;k++)H[j][k]=mod(H[j][k]*v);
    for(let i=j+1;i<n;i++)if(H[i][j]){
      const w=H[i][j];
      for(let k=j;k<=n;k++)H[i][k]=mod(H[i][k]-w*H[j][k]);
    }
  }
  const x=new Int32Array(n);
  for(let i=n-1;i>=0;i--){let v=H[i][n];for(let j=i+1;j<n;j++)v-=H[i][j]*x[j];x[i]=mod(v);}
  return x;
}
const X=(n)=>[...Array(n).keys()];
const key=(n,t)=>{const [i,j,k]=[...t].sort((a,b)=>a-b);return (i*n+j)*n+k;};
const held=(T,i)=>T.includes(i);
const dot=(a,b)=>mod(a.reduce((s,x,i)=>s+x*b[i],0));
function fixture(n,del,anchors){
  // Build the entire survivor incidence Gram directly from actual triples.
  const deleted=new Set(del.map(t=>key(n,t))),S=[],G=X(n).map(()=>new Int32Array(n));
  check(deleted.size===del.length,'distinct deleted triples');
  for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++){
    const U=[i,j,k];if(deleted.has(key(n,U)))continue;
    S.push(U);for(const a of U)for(const b of U)G[a][b]++;
  }
  const K=c2(n-1), d=G.map((r,i)=>K-r[i]),d0=Math.min(...d),M=K-d0,r=d.map(x=>x-d0);
  const selected=anchors.filter(i=>r[i]>0),aa=r.reduce((s,x,i)=>s+(selected.includes(i)?x:0),0),L=M-aa;
  const lambda=r.map((x,i)=>selected.includes(i)?1:mod(1-x*inv[M]));
  const V=[],C=[], MX=[], MV=[];
  for(const i of X(n)){
    if(r[i]===0)continue;
    const sites=X(n).filter(j=>j!==i);
    // An avoiding Gram computed by subtracting the actual contained-i rows.
    const J=sites.map(j=>sites.map(k=>G[j][k]-(j===k?G[i][j]:1-Number(deleted.has(key(n,[i,j,k]))))));
    const a=c2(n-2), on=selected.includes(i);
    const rhs=sites.map((j,r)=>on?mod(lambda[j]-J[r][r]*inv[a]):mod(-G[i][j]));
    const fit=solve(J,rhs), vals=new Int32Array(n);
    sites.forEach((j,k)=>vals[j]=fit[k]);
    (on?C:V).push({i,vals,base:on?inv[a]:0});
    const col=new Int32Array(n);
    col[i]=on?0:G[i][i];
    for(let r=0;r<sites.length;r++)col[sites[r]]=mod(dot(J[r],fit)+(on?J[r][r]*inv[a]:G[i][sites[r]]));
    for(const j of X(n))check(col[j]===(on?(j===i?0:lambda[j]):(j===i?K-d[i]:0)),'direct aggregate fitted star');
    (on?MX:MV).push(col);
  }
  const cVal=(c,U)=>held(U,c.i)?0:mod(c.base+U.reduce((s,j)=>s+c.vals[j],0));
  const vVal=(v,U)=>held(U,v.i)?1:mod(U.reduce((s,j)=>s+v.vals[j],0));
  const B=(U,W)=>{
    let rank=0, gu=1,gv=1;
    C.forEach(c=>{const u=cVal(c,U),v=cVal(c,W),w=r[c.i];rank=mod(rank+w*u*v);gu=mod(gu-w*u);gv=mod(gv-w*v);});
    V.forEach(v=>{rank=mod(rank+r[v.i]*vVal(v,U)*vVal(v,W)*mod(inv[M]*inv[K-d[v.i]]));});
    return mod(rank+gu*gv*inv[L]);
  };
  const sample=[];
  for(let s=0;s<Math.min(300,S.length);s++)sample.push(S[Math.floor(s*S.length/Math.min(300,S.length))]);
  let starCoordinates=0, spanRows=0;
  for(const U of sample){
    const cu=C.map(c=>cVal(c,U)),vu=V.map(v=>vVal(v,U)),gU=mod(1-cu.reduce((s,x,h)=>s+x*r[C[h].i],0));
    for(const i of X(n)){
      const gx=mod(G[i][i]-MX.reduce((s,x,h)=>s+r[C[h].i]*x[i],0));
      check(gx===mod(lambda[i]*L),'g incidence');
      let v=mod(gU*gx*inv[L]+cu.reduce((s,x,h)=>s+x*r[C[h].i]*MX[h][i],0));
      V.forEach((q,h)=>{v=mod(v+vu[h]*MV[h][i]*r[q.i]*mod(inv[M]*inv[K-d[q.i]]));});
      const ch=C.findIndex(c=>c.i===i),vh=V.findIndex(q=>q.i===i);
      const expect=ch>=0?mod(1-r[i]*cu[ch]):mod(lambda[i]+(vh>=0?r[i]*vu[vh]*inv[M]:0));
      check(v===expect,'coupled fit all old incidence columns');
      if(held(U,i)){check(v===1,'locked incident row-sum compatibility');starCoordinates++;}
      spanRows++;
    }
  }
  // Crossing and all endpoint branches: exact direct stars vs formula.
  let endpointControls=0,paletteControls=0;
  const bValues=[], minSites=X(n).filter(i=>d[i]===d0), maxSite=d.indexOf(Math.max(...d));
  for(const vs of [[],[minSites[0]],[maxSite],anchors,[...new Set([...anchors,minSites[0]])],X(n)]){
    const b=vs.length?Math.min(...vs.map(i=>r[i])):null;
    bValues.push(b);
    for(const [u,v] of [[0,0],[1,0],[0,1],[1,b??0],[3,b??0],[3,4*(b??0)+1]]){
      const S0=n+K, scores=[(vs.length+1)*v,...X(n).map(i=>S0*(u+v)-d[i]*u+(vs.includes(i)?v:0))];
      const expected=S0*(u+v)-d0*u+(vs.length?Math.max(0,v-b*u):0);
      check(Math.max(...scores)===expected,'marker, empty anchors and crossing extrema');
      endpointControls++;
    }
    const Mcolors=n%2?n:n+1, used=[...Array(n+1)].map(()=>new Set());
    function color(F,c){for(const i of F){check(!used[i].has(c),'pair/singleton colors are disjoint classes');used[i].add(c);}paletteControls++;}
    for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)color([i,j],(i+j)%Mcolors);
    for(const i of vs)color([i,n],(2*i)%Mcolors);
    for(let i=0;i<=n;i++){const c=X(n+1).find(c=>!used[i].has(c));check(c!==undefined,'unused singleton color');color([i],c);}
  }
  const D=Math.max(...d), hist={};for(const x of d)hist[x]=(hist[x]??0)+1;
  return {n,S,G,d,d0,M,r,B,stats:{n,deletions:del.length,d0,D,degrees:hist,positiveBridges:selected.length,unbridgedFits:V.length,selectedDegrees:[...new Set(selected.map(i=>d[i]))].sort((a,b)=>a-b),alpha:aa,actualRows:S.length,starCoordinates,spanRows,endpointControls,bValues,paletteControls}};
}

function tinyMask(f){
  const {n,S,G,B,M}=f,N=S.length, K=c2(n-1),q=c2(n-4), basis=X(n), shift=mod(inv[M]-inv[K]);
  check(shift!==0,'positive min-degree needed for omitted-offset negative control');
  const Gi=basis.map(i=>solve(G,basis.map(j=>Number(i===j))));
  const Pi=S.map(U=>Int32Array.from(S,W=>mod(U.reduce((s,i)=>s+W.reduce((t,j)=>t+Gi[i][j],0),0))));
  const F=S.map((U,i)=>Int32Array.from(S,(W,j)=>mod(Number(i===j)-Pi[i][j])));
  const mask=(U,W)=>U.some(i=>held(W,i));
  const Ps=S.map((U,i)=>Int32Array.from(S,(W,j)=>mod(Number(i===j)+(mask(U,W)?0:inv[q])-inv[K])));
  const product=(U,V)=>U.map(u=>Int32Array.from(S,(_,j)=>mod(u.reduce((s,x,i)=>s+x*V[i][j],0))));
  const P0=product(product(F,Ps),F);
  const D0=S.map((U,i)=>Int32Array.from(S,(W,j)=>mask(U,W)?mod(Number(i===j)-B(U,W)-P0[i][j]):0));
  const Wrong=S.map((U,i)=>Int32Array.from(S,(W,j)=>mask(U,W)?mod(D0[i][j]+shift):0));
  const Dhat=product(product(F,D0),F);
  const T=S.map((U,i)=>Int32Array.from(S,(W,j)=>mask(U,W)?Dhat[i][j]:0));
  const incid=(M,i,j)=>mod(M[i].reduce((s,v,k)=>s+(held(S[k],j)?v:0),0));
  let m0=0,m1=0,nonzero=0,nonidentity=0,wrongLocks=0;
  for(let i=0;i<N;i++){
    for(const j of basis){check(incid(P0,i,j)===0,'actual Kneser projection kills all old stars');m0++;}
    for(const j of S[i]){check(incid(D0,i,j)===0&&incid(T,i,j)===0,'compatibility before/after exact mask compression');m1++;
      check(incid(Wrong,i,j)===mod(shift*G[j][j])&&incid(Wrong,i,j)!==0,'omitting K/M offset breaks compatibility');wrongLocks++;}
    for(let j=0;j<N;j++){if(D0[i][j]!==0)nonzero++;if(T[i][j]!==D0[i][j])nonidentity++;}
  }
  check(nonzero>0 && nonidentity>0,'nontrivial mask stress');
  return {n,minDegree:f.d0,rows:N,m0,m1,nonzero,nonidentity,wrongLocks};
}
// Unlike the other replay, these deletions are cyclic translate orbits, not three-partite.
function orbit(n,h){return X(n).map(i=>h.map(j=>(i+j)%n));}
function inputs(n,restore){
 const starts=[[0,1,3],[0,1,4],[0,2,6],[0,2,7]], missing=new Set(restore?starts.slice(0,3).map(t=>key(n,t)):[]);
 return starts.flatMap(h=>orbit(n,h)).filter(t=>!missing.has(key(n,t)));
}
function summarize(f,anchors){
  const z=f.stats, un=f.r.filter((_,i)=>!anchors.includes(i)), e=Math.max(0,...un), su=un.reduce((s,x)=>s+x,0);
  return {...z, unbridgedSurplus:su, rationalCriterion:bound(z.n,z.D,Q(z.deletions),Q(z.alpha),z.d0,e,Q(su))};
}
const anc=X(100).filter(i=>i%4===1||i===4);
const grown=fixture(100,inputs(100,true),anc);
check(grown.stats.d0===9&&grown.stats.D===12&&grown.stats.deletions===397&&grown.stats.selectedDegrees.join(',')==='10,11,12','irregular all-touched boundary fixture');
const regular=fixture(100,inputs(100,false),X(100));
check(regular.stats.d0===12&&regular.stats.D===12&&regular.stats.deletions===400&&regular.stats.positiveBridges===0,'regular no-surplus full anchor fixture');
const tDel=orbit(10,[0,1,3]).filter(t=>key(10,t)!==key(10,[0,1,3])).concat(orbit(10,[0,1,4]));
const tiny=tinyMask(fixture(10,tDel,[0,2,5,8]));
const out={arithmetic:'exact modular and BigInt rational; displayed squared ratios approximate only',coefficients:coarse.map(([x,y])=>[x.map(String),y]),uniformGap:[Q(1669,2000).map(String),Q(76483,90000).map(String)],irregular:summarize(grown,anc),regular:summarize(regular,X(100)),tinyMask:tiny};
console.log(JSON.stringify(out,null,2));
