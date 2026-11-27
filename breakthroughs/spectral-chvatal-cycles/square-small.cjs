"use strict";
const c2=x=>x*(x-1)/2,c3=x=>x*(x-1)*(x-2)/6, A=(n,f)=>Array.from({length:n},(_,i)=>f(i));
function parts(r,m=3){if(!r)return [[]];return A(r-m+1>0?r-m+1:0,j=>m+j).flatMap(k=>parts(r-k,k).map(x=>[k,...x]));}
function tri(n){const S=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++)S.push([i,j,k]);return S;}
const msk=S=>S.reduce((z,t)=>z|1<<t,0);
function structure(n,cycles){
 const r=cycles.reduce((a,b)=>a+b,0),K=c2(n-1),a=c2(n-2),q=c2(n-4),L=K-2*r;
 const D=[];let off=0;
 for(const k of cycles){for(let i=0;i<k;i++)D.push([off+i,off+(i+1)%k,r+off+i]);off+=k;}
 const d=D.map(msk),U=tri(n).filter(t=>!d.includes(msk(t))),u=U.map(msk);
 return {n,r,K,a,q,L,D,d,U,u};
}
function factors(S,O){
 const {n,r,K,a,q,L,D,d,U,u}=S, agg=numeric=>O.make(numeric),J=A(n,i=>i),count=r;
 const old=O.gram(J,D), invAx=old;
 function adjusted(i,fn){
  const b=J.filter(j=>j!==i), G=O.gram(b,D.filter(x=>!x.includes(i))),v=fn(b,G);
  return A(n,j=>j===i?agg(0):v[b.indexOf(j)]);
 }
 const C=A(r,i=>adjusted(i,(b,G)=>O.solve(G,b.map(j=>O.sub(j>=r&&j<2*r?O.rat(K-1,K):agg(1),O.rat(a-D.filter(x=>x.includes(j)&&!x.includes(i)).length,a))))));
 const V=A(r,t=>{const i=t+r;return adjusted(i,(b,G)=>O.solve(G,b.map(j=>agg(-(n-2-D.filter(x=>x.includes(i)&&x.includes(j)).length)))));});
 const E=D.map(Di=>O.solve(invAx,J.map(j=>agg(+Di.includes(j)))));
 const P=d.map((v,i)=>d.map((w,j)=>O.sub(i===j?agg(1):(v&w)?agg(0):O.rat(1,q),O.rat(1,K))));
 return U.map((T,t)=>{
  const sums=(x)=>T.reduce((s,i)=>O.add(s,x[i]),agg(0));
  const w=C.map((x,i)=>O.sub(T.includes(i)?agg(0):O.add(O.rat(1,a),sums(x)),O.rat(1,K)));
  const h=w.reduce(O.add,agg(0));
  const v=V.map((x,i)=>T.includes(i+r)?agg(1):sums(x));
  const e=E.map(sums),b=d.map(x=>O.sub(x&u[t]?agg(0):O.rat(1,q),O.rat(1,K)));
  const eP=P.map(row=>row.reduce((z,x,i)=>O.add(z,O.mul(x,e[i])),agg(0)));
  return [w.concat([h],v,e,b,e),w.map(x=>O.mul(agg(2),x)).concat([O.mul(O.rat(4,L),h)],
          v.map(x=>O.mul(O.rat(1,K*(K-1)),x)),b,e,eP)];
 });
}
function Ogram(J,D,o){
 const m=J.length,a=c2(m-2), b=m-2;
 return J.map(i=>J.map(j=>o((i===j?a:0)+b-D.filter(x=>x.includes(i)&&x.includes(j)).length)));
}
function Gaussian(M,b,OS){const n=b.length,C=M.map((a,i)=>a.concat([b[i]]));
 const {sub,mul,div,zero}=OS;
 for(let j=0;j<n;j++){
  if(zero(C[j][j]))throw Error("positive pivot required");
  const v=C[j][j];
  for(let t=j;t<=n;t++)C[j][t]=div(C[j][t],v);
  for(let i=j+1;i<n;i++){const c=C[i][j];if(!zero(c))for(let t=j;t<=n;t++)C[i][t]=sub(C[i][t],mul(c,C[j][t]));}
 }
 const X=A(n,()=>null);
 for(let j=n-1;j>=0;j--){let a=C[j][n];for(let t=j+1;t<n;t++)a=sub(a,mul(C[j][t],X[t]));X[j]=a;}
 return X;
}
const float={make:x=>x,rat:(x,y)=>x/y,add:(x,y)=>x+y,sub:(x,y)=>x-y,mul:(x,y)=>x*y,div:(x,y)=>x/y,zero:x=>x===0,
 gram:(J,D)=>Ogram(J,D,x=>x),solve:(M,b)=>Gaussian(M,b,float)};
function lockedStruct(S){
 const {n,D,U,u}=S,M=U.length, I0=A(M,()=>[]),J0=A(M,()=>[]),L=[],R=[],W=[];
 for(let i=0;i<M;i++)for(let j=0;j<=i;j++)if(u[i]&u[j]){
  const k=L.length;L.push(i);R.push(j);W.push(i===j?1:2);I0[i].push(k);if(i!==j)J0[j].push(k);
 }
 return {L:Int32Array.from(L),R:Int32Array.from(R),W:Uint8Array.from(W),ld:L.length,I0,J0};
}
function approx(n,cycles,tol=2e-15){
 const S=structure(n,cycles),{U,D}=S,m=U.length,F=factors(S,float),k=F[0][0].length,Lck=lockedStruct(S),{L,R,W,ld}=Lck;
 const G0=float.gram(A(n,j=>j),D),B=A(n,i=>float.solve(G0,A(n,j=>+(i===j)))); // symmetric inverse
 const dot=(u,v)=>{let s=0;for(let k=0;k<ld;k++)s+=W[k]*u[k]*v[k];return s;};
 const rhs=Float64Array.from(L,(i,e)=>{let s=0;for(let j=0;j<k;j++)s+=F[i][0][j]*F[R[e]][1][j];return -s;});
 const trans=Z=>{
  const Q=A(m,()=>new Float64Array(n));
  for(let k=0;k<ld;k++){
   const i=L[k],j=R[k],v=Z[k], a=Q[i],b=Q[j],U0=U[i],U1=U[j];
   a[U1[0]]+=v;a[U1[1]]+=v;a[U1[2]]+=v;
   if(i!==j){b[U0[0]]+=v;b[U0[1]]+=v;b[U0[2]]+=v;}
  }
  const Y=Q.map(q=>Float64Array.from(B,row=>{let s=0;for(let i=0;i<n;i++)s+=q[i]*row[i];return s;}));
  const T=A(n,()=>new Float64Array(n));
  for(let a=0;a<m;a++){const v=U[a];for(const i of v)for(let j=0;j<n;j++)T[i][j]+=Y[a][j];}
  const H=B.map(row=>Float64Array.from(row,(_,j)=>{let s=0;for(let i=0;i<n;i++)s+=row[i]*T[i][j];return s;}));
  return Float64Array.from(Z,(v,k)=>{
   const i=L[k],j=R[k],a=U[i],b=U[j];let s=v;
   for(const q of b)s-=Y[i][q];for(const q of a)s-=Y[j][q];
   for(const x of a)for(const y of b)s+=H[x][y];return s;
  });
 };
 const lin=new Float64Array(ld);
 let d=rhs.slice(),q=rhs.slice(),p=dot(q,q),it=0;
 for(;it<150&&p>tol**2;it++){
  const r=trans(d),z=p/dot(d,r);
  if(!Number.isFinite(z)||z<=0)throw Error("CG "+[n,cycles,it,p,z]);
  for(let i=0;i<ld;i++){lin[i]+=z*d[i];q[i]-=z*r[i];}
  const np=dot(q,q),c=np/p;p=np;
  for(let i=0;i<ld;i++)d[i]=q[i]+c*d[i];
 }
 const actual=trans(lin),r0=Float64Array.from(actual,(v,j)=>v-rhs[j]);
 return {S,Lck,lin,it,ratio:Math.sqrt(dot(lin,lin)),rhs:Math.sqrt(dot(rhs,rhs)),err:Math.sqrt(dot(r0,r0)),rel:Math.max(...r0.map(Math.abs))};
}
function scope(){return A(6,i=>[9+i,[4]]);}

function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b)[a,b]=[b,a%b];return a;}
function lcm(a,b){return a/gcd(a,b)*b;}
const z0=[0n,1n];
function Q(a,b=1n){a=BigInt(a);b=BigInt(b);if(!b)throw Error("division");if(!a)return z0;if(b<0n){a=-a;b=-b;}const g=gcd(a,b);return [a/g,b/g];}
const add=(a,b)=>Q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),sub=(a,b)=>Q(a[0]*b[1]-b[0]*a[1],a[1]*b[1]),
 mul=(a,b)=>Q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>Q(a[0]*b[1],a[1]*b[0]),sq=a=>mul(a,a), zero=a=>a[0]===0n, lt=(a,b)=>a[0]*b[1]<b[0]*a[1],far=a=>a[0]+"/"+a[1];
const abs=a=>a<0n?-a:a,ck=(b,s)=>{if(!b)throw Error(s)}, checkSum=(a)=>a.reduce(add,z0);
const O={make:x=>Q(x),rat:Q,add,sub,mul,div,zero,gram:(J,D)=>Ogram(J,D,x=>Q(x)),solve:(M,b)=>Gaussian(M,b,O)};
function integerMat(H){
 const den=H.flat().reduce((s,v)=>lcm(s,v[1]),1n);
 return {den,H:H.map(a=>a.map(v=>v[0]*(den/v[1])))};
}
function integerFactors(F){
 const k=F[0][0].length;
 const P=A(k,t=>[F.reduce((p,u)=>lcm(p,u[0][t][1]),1n),F.reduce((p,u)=>lcm(p,u[1][t][1]),1n)]);
 const d=P.reduce((s,[a,b])=>lcm(s,a*b),1n);
 return {den:d,l:F.map(row=>row[0].map((x,t)=>x[0]*(P[t][0]/x[1]))),r:F.map(row=>row[1].map((x,t)=>x[0]*(d/P[t][0]/x[1])))};
}

function orbit4(S){
 const {U,n,u}=S, maps=[];
 for(let t=0;t<4;t++)for(const s of [1,3]){
  const g=A(4,i=>(t+s*i)%4), d=A(4,i=>g[(i+1)%4]===(g[i]+1)%4?g[i]+4:g[(i+1)%4]+4);
  const h=g.concat(d);
  maps.push({h,tab:A(256,b=>h.reduce((m,j,i)=>(b>>i&1)?m|1<<j:m,0))});
 }
 // On untouched leaves an orbit is specified only by its common and unique counts.
 const as=U.map(T=>T.filter(i=>i>=8).length), pop=x=>{let r=0;while(x){x&=x-1;r++;}return r;};
 const mapPair=(i,j)=>{
  const a=u[i]&255,b=u[j]&255,k=pop((u[i]&u[j])>>>8),h=as[i]-k,v=as[j]-k;
  let m=10**9;
  for(const {tab}of maps){
   const p=tab[a],q=tab[b];m=Math.min(m,(((256*p+q)*4+k)*4+h)*4+v,(((256*q+p)*4+k)*4+v)*4+h);
  }
  return m;
 };
 const mapRow=(i,j)=>Math.min(...maps.map(({h,tab})=>(tab[u[i]&255]*4+as[i])*9+(j>=8?8:h[j])));
 return {mapPair,mapRow};
}
function orbitRepair(S,Lck,residual,den){
 const {n,r,U,u}=S,{L,R,ld}=Lck,M=u.length,I=A(n,i=>i).filter(i=>i<r||i>=2*r),rhs=A(M,()=>Array(n).fill(0n)),
       z=Array(M).fill(0n), OR=orbit4(S),packed=[],ids=new Map(),varlist=[];
 // Compute exact row residuals, including all tight coordinates.
 for(let k=0;k<ld;k++){
  const a=L[k],b=R[k],w=residual[k];for(const j of U[b])rhs[a][j]-=w;
  if(a!==b)for(const j of U[a])rhs[b][j]-=w;
  z[a]+=abs(w);if(a!==b)z[b]+=abs(w);
 }
 ck(rhs.every((v,i)=>U[i].every(j=>v[j]===0n)),"all-old incident compatibility");
 const b0=new Map(),row=A(M,()=>[]);
 for(let a=0;a<M;a++)for(let b=0;b<a;b++)if(!(u[a]&u[b])){
  const key=OR.mapPair(a,b);if(!ids.has(key)){ids.set(key,ids.size);varlist.push(key);}
  const id=ids.get(key);row[a].push([b,id]);row[b].push([a,id]);packed.push([a,b,id]);
 }
 for(let a=0;a<M;a++)for(const j of I)if(!U[a].includes(j)){
  const key=OR.mapRow(a,j),w=rhs[a][j];
  if(b0.has(key)){ck(b0.get(key).rhs===w,"rhs invariant actual anchor/leaf symmetry");continue;}
  const C=new Map();for(const [b,k]of row[a])if(U[b].includes(j))C.set(k,(C.get(k)??0)+1);
  b0.set(key,{rhs:w,c:C});
 }
 const build=[];let allp=0;
 for(const box of b0.values()){
  let rhs=Q(box.rhs);const C=new Map([...box.c].map(([i,v])=>[i,Q(v)]));
  for(const {key,v,T}of build){const a=C.get(key);if(a&&!zero(a)){
   rhs=sub(rhs,mul(a,v));
   for(const [i,b]of T){const x=sub(C.get(i)??z0,mul(a,b));if(zero(x))C.delete(i);else C.set(i,x);}
  }}
  if(!C.size){ck(zero(rhs),"orbit symmetric range exact");continue;}
  const key=Math.min(...C.keys()),a=C.get(key);for(const [i,q]of C)C.set(i,div(q,a));
  build.push({key,v:div(rhs,a),T:C});allp++;
 }
 const w=A(ids.size,()=>z0);
 for(const {key,v,T}of build.slice().reverse()){
  let a=v;for(const [j,q]of T)if(j!==key)a=sub(a,mul(q,w[j]));
  w[key]=a;
 }
 const ed=w.reduce((a,b)=>lcm(a,b[1]),1n),val=w.map(x=>x[0]*(ed/x[1]));
 let nnz=packed.filter(([a,b,k])=>val[k]!==0n).length;
 const check=A(M,()=>Array(n).fill(0n));
 for(const [a,b,k]of packed)if(val[k]){
  const w=val[k];for(const j of U[b])check[a][j]+=w;for(const j of U[a])check[b][j]+=w;
 }
 ck(check.every((a,i)=>I.every(j=>a[j]===rhs[i][j]*ed)),"full disjoint orbit coupling all old rows");
 const oldAbs=A(M,()=>0n);for(const [a,b,k]of packed)if(val[k]){oldAbs[a]+=abs(val[k]);oldAbs[b]+=abs(val[k]);}
 const norm=z.reduce((a,b,i)=>a>b*ed+oldAbs[i]?a:b*ed+oldAbs[i],0n);
 return {norm:Q(norm,ed*den),nnz,eqs:b0.size,pivots:allp,orbits:ids.size};
}
function certify(n,cycles){
 const {S,Lck,lin,it}=approx(n,cycles);
 const {r,K,L,U,D}=S,M=U.length,cli=Lck,{ld,R,W}=cli,IX=cli.L;
 const factors0=factors(S,O),C0=integerFactors(factors0),X0=A(n,i=>i),G=O.gram(X0,D);
 const Inv=integerMat(A(n,i=>O.solve(G,A(n,j=>Q(+(i===j))))));
 const primitive=Inv.H,gden=Inv.den,Qd=1n<<48n;
 const orbit=orbit4(S),g=new Map();
 const di=Array.from(lin,(t,k)=>{const p=orbit.mapPair(IX[k],R[k]);if(!g.has(p))g.set(p,BigInt(Math.round(t*Number(Qd))));return g.get(p);});
 const mu=Q(S.q-1,S.q),alpha=n===9?Q(1,100):Q(1,4); const dX=A(M,()=>Array(n).fill(0n));
 for(let e=0;e<ld;e++){
  const i=IX[e],j=R[e],v=di[e];
  for(const w of U[j])dX[i][w]+=v;
  if(i!==j)for(const w of U[i])dX[j][w]+=v;
 }
 const Y=dX.map(row=>primitive.map(H=>H.reduce((s,x,j)=>s+x*row[j],0n)));
 const XY=A(n,()=>Array(n).fill(0n));
 for(let i=0;i<M;i++)for(const t of U[i])for(let j=0;j<n;j++)XY[t][j]+=Y[i][j];
 const H=primitive.map(row=>A(n,j=>row.reduce((s,a,k)=>s+a*XY[k][j],0n)));
 const td=Qd*gden*gden,den=lcm(td,C0.den);
 let sumproj=0n;
 const residual=A(ld,k=>{
  const i=IX[k],j=R[k],a=U[i],b=U[j];
  let T=di[k]*gden*gden;
  for(const u of a)T-=Y[j][u]*gden;
  for(const v of b)T-=Y[i][v]*gden;
  for(const u of a)for(const v of b)T+=H[u][v];
  const D0=-C0.l[i].reduce((s,a,t)=>s+a*C0.r[j][t],0n);
  sumproj+=BigInt(W[k])*di[k]*T;
  return D0*(den/C0.den)-T*(den/td);
 });
 const a0=Q(sumproj,Qd*td);
 ck(a0[0]>=0n&&lt(a0,sq(sub(mu,alpha))),"projected Kneser spectral margin");
 const e=orbitRepair(S,Lck,residual,den);
 const ceil=v=>{let s=0n;while(s*s*v[1]<v[0])s++;return s;};
 const gp=ceil(checkSum(factors0.map(J=>sq(sub(Q(L,K),mul(Q(2),J[0][r]))))));
 const vp=ceil(checkSum(factors0.flatMap(J=>J[0].slice(r+1,2*r+1).map(sq))));
 const cp=ceil(Q(3*K)),rp=1n+ceil(Q(r));
 const A0=add(Q(1),mul(Q(cp),add(Q(gp,L),Q(rp*vp,K-1)))),A1=Q(cp,L),A2=Q(cp*rp,K-1);
 const bound=div(Q(1),checkSum([div(sq(A0),alpha),mul(sq(A1),Q(L)),mul(sq(A2),Q(K*(K-1)))]));
 ck(lt(e.norm,bound),"strict exact residual Schur margin");
 return {n,cycles,iterations:it,locked:ld,matrixRows:M,spectralMargin:far(alpha),rankGap:far(bound),
   repairedNorm:far(e.norm),repairNonzeros:e.nnz,denBits:den.toString(2).length,integerProjectedSquared:far(a0),eqs:e.eqs,pivots:e.pivots,orbits:e.orbits,pass:true};
}
console.log(JSON.stringify(scope().map(([n,c])=>certify(n,c))));

