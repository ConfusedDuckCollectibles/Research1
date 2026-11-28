'use strict';
// Independent exact test of heterogeneous triple weights and multiple markers. No imports.
const ok=(x,label)=>{if(!x)throw Error(label);};
const absB=x=>x<0n?-x:x;
function gcd(a,b){while(b)[a,b]=[b,a%b];return absB(a);}
const lcm=(a,b)=>a/gcd(a,b)*b;
const q=(a,b=1n)=>{a=BigInt(a);b=BigInt(b);ok(b!==0n,'rational denominator');const g=gcd(a,b)*(b<0n?-1n:1n);return [a/g,b/g];};
const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),sub=(a,b)=>add(a,[-b[0],b[1]]);
const mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]);
const sq=a=>mul(a,a),sum=a=>a.reduce(add,q(0)),gt=(a,b)=>a[0]*b[1]>b[0]*a[1];
const is=(a,b)=>a[0]===b[0]&&a[1]===b[1],abs=a=>[absB(a[0]),a[1]];
const max=a=>a.reduce((x,y)=>gt(x,y)?x:y,q(0));
const show=x=>x.map(String).join('/'),flo=x=>Number((Number(x[0])/Number(x[1])).toFixed(10));
const N=n=>Array.from({length:n},(_,i)=>i), c2=n=>n*(n-1)/2,c3=n=>n*(n-1)*(n-2)/6;
function solve(A,b){ // exact rational Gauss-Jordan, independent of modular replay
 const n=A.length,h=A.map((row,i)=>[...row.map(v=>Array.isArray(v)?v:q(v)),b[i]]);
 for(let i=0;i<n;i++){
  let p=i;while(p<n&&h[p][i][0]===0n)p++;ok(p<n,'nonsingular incidence fit');[h[i],h[p]]=[h[p],h[i]];
  const z=h[i][i];for(let k=i;k<=n;k++)h[i][k]=div(h[i][k],z);
  for(let j=0;j<n;j++)if(j!==i&&h[j][i][0]!==0n){const y=h[j][i];for(let k=i;k<=n;k++)h[j][k]=sub(h[j][k],mul(y,h[i][k]));}
 }
 return h.map(a=>a[n]);
}
const gram=(S,n)=>{const G=N(n).map(()=>N(n).map(()=>0));for(const U of S)for(const i of U)for(const j of U)G[i][j]++;return G;};
const idxset=U=>U.slice().sort((a,b)=>a-b).join(',');
const dotB=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0n);

function setup(){
 const n=9, K=c2(n-1),q0=c2(n-4), a=c2(n-2),Z=N(n),S=[];
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++)S.push([i,j,k]);
 const tweaks=new Map([['0,1,2',0],['3,4,5',0],['0,3,6',6],['1,4,7',9],['2,5,8',8],['0,5,7',4]]);
 const fd=12n, fps=S.map(t=>BigInt(tweaks.get(idxset(t))??12)), fn=fps.map(t=>q(t,fd)), p=fn.map(sq), rows=N(S.length);
 const D=rows.map(i=>fd*fd-fps[i]*fps[i]), live=rows.filter(i=>fps[i]>0n), zeros=rows.filter(i=>fps[i]===0n);
 const wc=(i,j,avoid=-1)=>sum(rows.filter(z=>S[z].includes(i)&&S[z].includes(j)&&!S[z].includes(avoid)).map(z=>p[z]));
 const G=Z.map(i=>Z.map(j=>wc(i,j))), t=Z.map(i=>G[i][i]), d=t.map(x=>sub(q(K),x)), dmax=max(d);
 const Ginv=Z.map(i=>solve(G,Z.map(j=>q(Number(i===j))))), gd=Ginv.flat().reduce((a,x)=>lcm(a,x[1]),1n);
 const Gnum=Ginv.map(row=>row.map(x=>x[0]*(gd/x[1]))),XGi=S.map(U=>Z.map(j=>U.reduce((s,i)=>s+Gnum[i][j],0n)));
 const mask=(i,j)=>S[i].some(x=>S[j].includes(x)),pd=BigInt(K*q0);
 const Pn=(i,j)=>(i===j?pd:0n)+(mask(i,j)?0n:BigInt(K))-BigInt(q0);
 const Ln=(i,j)=>(i===j?pd:0n)-Pn(i,j);
 const rd=fd*fd*pd,R=rows.map(i=>rows.map(j=>fps[i]*fps[j]*Pn(i,j)+(i===j?D[i]*pd:0n)));
 const pid=fd*fd*gd,Pi=rows.map(i=>rows.map(j=>fps[i]*fps[j]*S[j].reduce((v,k)=>v+XGi[i][k],0n)));
 const dense=(A,B)=>rows.map(i=>rows.map(j=>rows.reduce((v,k)=>v+A[i][k]*B[k][j],0n)));
 const P1=dense(Pi,R), P2=dense(P1,Pi), p0d=pid*pid*rd;
 const P0=rows.map(i=>rows.map(j=>R[i][j]*pid*pid-(P1[i][j]+P1[j][i])*pid+P2[i][j]));
 let first=0,quadratic=0,locks=0,positiveNull=0,zeroComparison=0;
 const RX=rows.map(i=>Z.map(j=>rows.reduce((a,k)=>a+(S[k].includes(j)?R[i][k]*fps[k]:0n),0n)));
 const DX=rows.map(i=>Z.map(j=>rows.reduce((a,k)=>a+(S[k].includes(j)?Ln(i,k)*D[k]:0n),0n)));
 for(const i of rows){
  for(const j of Z){
   ok(RX[i][j]===fps[i]*DX[i][j],'new seed first-order cancellation');first++;
   ok(rows.reduce((v,k)=>v+(S[k].includes(j)?P0[i][k]*fps[k]:0n),0n)===0n,'P0 weighted incidence null');positiveNull++;
  }
  for(const j of rows)if(mask(i,j)){
   ok(R[i][j]===(i===j?rd:0n)-fps[i]*fps[j]*BigInt(q0),'exact seed locked entry');locks++;
  }
 }
 for(const i of Z)for(const j of Z){
  const lhs=rows.reduce((b,k)=>b+(S[k].includes(i)?fps[k]*RX[k][j]:0n),0n);
  const rhs=rows.reduce((b,k)=>b+(S[k].includes(i)?D[k]*((S[k].includes(j)?fd*fd*pd:0n)-DX[k][j]):0n),0n);
  ok(lhs===rhs,'second-order incidence identity, direct sum versus L deficit product');quadratic++;
 }
 for(const i of zeros)for(const j of rows){
  ok(Pi[i][j]===0n&&R[i][j]===(i===j?rd:0n)&&P0[i][j]===(i===j?p0d:0n),'zero rows decouple exactly');zeroComparison++;
 }
 // Verify the exact congruence gap using a signed incidence trade plus a zero-weight basis vector.
 const id=(x,a,b)=>S.findIndex(z=>idxset(z)===idxset([x,a,b]));let trade;
 for(const x of Z)for(const a of Z)for(const b of Z)for(const c of Z)for(const d of Z){
  if(trade||new Set([x,a,b,c,d]).size!==5)continue;
  const v=[id(x,a,b),id(x,c,d),id(x,a,c),id(x,b,d)];
  if(v.every(i=>fps[i]>0n)&&v.some(i=>fps[i]!==fd))trade=v;
 }
 ok(trade,'heterogeneous trade');
 const support=[...trade,zeros[0]], ys=trade.map((i,j)=>div(q(j<2?1:-1),fn[i])).concat([q(1)]);
 for(const j of Z)ok(sum(support.map((i,k)=>S[i].includes(j)?mul(ys[k],fn[i]):q(0)))[0]===0n,'trade in weighted kernel');
 const norm2=sum(ys.map(sq)), ray=sum(support.flatMap((i,k)=>support.map((j,l)=>mul(mul(ys[k],ys[l]),q(R[i][j],rd)))));
 const gap=mul(sub(q(1),q(1,q0)),norm2);
 ok(gt(ray,gap),'new seed gap on mixed positive/zero weighted complement sample');
 return {n,K,a,Z,S,rows,live,zeros,fn,p,wc,G,t,d,dmax,mask,p0d,P0,stats:{n,zeroRows:zeros.length,weights:[...new Set(p.map(show))].sort(),
   first,quadratic,locks,positiveNull,zeroComparison,gapSupport:support.map(i=>S[i]),gap:show(gap),ray:show(ray)}};
}
function multiple(A,scale){
 const {n,K,a,Z,S,rows,live,zeros,fn,p,wc,G,t,d,dmax,mask,p0d,P0}=A;
 const endpoint=[[8,0,q(3,4)],[8,1,q(2,3)],[6,0,q(1,3)],[6,2,q(1,2)],[2,1,q(2,5)],[0,2,q(3,5)],[4,2,q(0)]];
 const edges=endpoint.map(([i,j,z])=>({i,j,w:mul(z,scale)})).filter(e=>e.w[0]>0n), E=N(edges.length);
 const singletons=[q(1),q(3,5),q(2,3),q(2,7),q(0)]; // marker 3 has no bridges; marker 4 omitted
 const w=Z.map(i=>sum(edges.filter(e=>e.i===i).map(e=>e.w))), W=max(w), alpha=sum(w);
 const M=max(Z.map(i=>add(t[i],w[i]))),L=sub(M,alpha),lam=Z.map(i=>div(t[i],sub(M,w[i])));
 const bridged=Z.filter(i=>w[i][0]>0n),slack=Z.filter(i=>gt(q(1),lam[i])),tight=Z.filter(i=>is(lam[i],q(1)));
 ok(tight.length===1&&tight[0]===8&&edges.filter(e=>e.i===8).length===2&&slack.includes(6)&&slack.includes(2),'tight multi-bridge anchor plus slack multi-bridge anchor');
 const C=new Map(),V=new Map();
 function fit(i,bridge){
  const ids=Z.filter(j=>j!==i), A=ids.map(j=>ids.map(k=>wc(j,k,i)));
  const b=ids.map((j,k)=>bridge?sub(lam[j],div(A[k][k],q(a))):mul(q(-1),G[i][j])),y=solve(A,b);
  return S.map(U=>U.includes(i)?q(bridge?0:1):add(bridge?q(1,a):q(0),sum(U.map(j=>y[ids.indexOf(j)]))));
 }
 for(const i of bridged)C.set(i,fit(i,true));
 for(const i of slack)V.set(i,fit(i,false));
 const inc=(col,i)=>sum(live.filter(j=>S[j].includes(i)).map(j=>mul(p[j],col[j])));
 let fitCount=0,rankLocks=0,tightOutputs=0,oldProjLocks=0,wrong=0,dupEquals=0,disjointTwinEntries=0;
 for(const [i,c] of C)for(const j of Z){ok(is(inc(c,j),j===i?q(0):lam[j]),'exact heterogeneous weighted C fit');fitCount++;}
 for(const [i,v] of V)for(const j of Z){ok(is(inc(v,j),j===i?t[i]:q(0)),'exact heterogeneous weighted V fit');fitCount++;}
 const g=rows.map(j=>sub(q(1),sum(bridged.map(i=>mul(w[i],C.get(i)[j]))))),F=new Map(slack.map(i=>[i,div(sub(q(1),lam[i]),t[i])]));
 const Bhat=rows.map(i=>rows.map(j=>mul(mul(fn[i],fn[j]),sum([
  div(mul(g[i],g[j]),L),...bridged.map(k=>mul(w[k],mul(C.get(k)[i],C.get(k)[j]))),
  ...slack.map(k=>mul(F.get(k),mul(V.get(k)[i],V.get(k)[j])))]))));
 const CC=rows.map(i=>E.map(j=>mul(fn[i],C.get(edges[j].i)[i])));
 const oldInc=(H,i,j)=>sum(live.filter(z=>S[z].includes(j)).map(z=>mul(H[i][z],fn[z])));
 for(const i of Z){
  ok(is(inc(g,i),mul(L,lam[i])),'normalized g identity');
  for(const j of rows){
   const old=oldInc(Bhat,j,i),byBridge=sum(E.filter(e=>edges[e].i===i).map(e=>mul(CC[j][e],edges[e].w)));
   const rhs=mul(fn[j],add(mul(lam[i],sub(q(1),C.has(i)?mul(w[i],C.get(i)[j]):q(0))),V.has(i)?mul(sub(q(1),lam[i]),V.get(i)[j]):q(0)));
   ok(is(old,rhs),'fitted rank update');
   if(S[j].includes(i)){ok(is(old,fn[j]),'weighted rank lock');rankLocks++;}
   if(tight.includes(i)&&live.includes(j)){
    ok(is(div(add(old,byBridge),fn[j]),q(1)),'physical triple output at repeated-bridge tight star');tightOutputs++;
   }
  }
  if(tight.includes(i))for(const e of E){
   ok(is(add(sum(live.filter(j=>S[j].includes(i)).map(j=>mul(fn[j],CC[j][e]))),q(Number(edges[e].i===i))),q(1)),'every physical bridge row of tight star');tightOutputs++;
  }
 }
 for(const U of live){
  for(let e=0;e<E.length;e++)for(let f=e+1;f<E.length;f++){
   if(edges[e].i===edges[f].i){ok(is(CC[U][e],CC[U][f]),'duplicated same-anchor column');dupEquals++;}
   else if(edges[e].j!==edges[f].j)disjointTwinEntries++;
  }
  for(const i of S[U]){
   const valid=sum(live.filter(V=>S[V].includes(i)&&mask(U,V)).map(V=>mul(sub(sub(q(Number(U===V)),Bhat[U][V]),q(P0[U][V],p0d)),fn[V])));
   ok(valid[0]===0n,'actual weighted positive-row defect is compatible');oldProjLocks++;
   const counterfactual=mul(fn[U],mul(t[i],sub(div(q(1),M),q(1,K))));
   ok(counterfactual[0]!==0n,'omitting signed heterogeneous center shift breaks every lock');wrong++;
  }
 }
 for(const U of zeros)for(const V of rows)ok(Bhat[U][V][0]===0n&&is(q(P0[U][V],p0d),q(Number(U===V))),'zero triple rows valid prior to restriction');
 // Exact bridge metric s diag(1/w) - J at the full old cap with singleton projections.
 const s=add(q(n),M), QB=E.map(i=>E.map(j=>sub(i===j?div(s,edges[i].w):q(0),q(1))));
 const markerProj=[], qm=[];
 let markerLocks=0,unlockedCross=0,markerNorms=[];
 for(let j=0;j<singletons.length;j++){
  const b=singletons[j];if(b[0]===0n)continue;
  const qa=sum(edges.filter(e=>e.j===j).map(e=>e.w)), coeff=edges.map(e=>e.j===j?mul(q(-1),div(e.w,sub(s,qa))):q(0));
  const product=E.map(i=>sum(E.map(k=>mul(QB[i][k],coeff[k]))));
  for(const i of E){
   if(edges[i].j===j){ok(is(product[i],q(-1)),'marker lock');markerLocks++;}
   else if(product[i][0]!==0n)unlockedCross++;
  }
  const v=sum(E.map(i=>mul(coeff[i],product[i]))), required=sub(div(s,b),q(1));
  ok(is(v,div(qa,sub(s,qa)))&&gt(required,v),'singleton independent orthogonal top-up');
  markerNorms.push(show(v)); markerProj.push(coeff);qm.push(qa);
 }
 const cross=markerProj.flatMap((v,i)=>markerProj.slice(i+1).map((u,j)=>{
   const got=sum(E.flatMap(a=>E.map(b=>mul(mul(v[a],QB[a][b]),u[b]))));
   const p=qm[i],q0=qm[i+1+j];
   ok(is(got,mul(q(-1),div(mul(p,q0),mul(sub(s,p),sub(s,q0))))),'cross-marker disjoint unlocked inner product');
   return show(got);
 }));
 return {scale:show(scale),D:show(dmax),W:show(W),alpha:show(alpha),tight,MMinusK:show(sub(M,q(K))),theta:show(s),
    fitCount,rankLocks,tightOutputs,oldProjLocks,wrong,dupEquals,disjointBridgePairs:disjointTwinEntries/live.length,
    markerLocks,unlockedCross,markerNorms,cross};
}
function constants(){
 const a=q(99,200),u=q(12,25),v=q(39,80),rho=q(99,320);
 const first=div(q(6),mul(a,u)), quad=div(mul(q(9),add(q(1),div(q(1),mul(q(1000),a)))),sq(v));
 const total=div(q(3*126,800),sub(q(1),rho));
 const low=div(mul(q(3*123,2),add(q(1,450),q(1,8000))),sub(q(1),rho));
 const high=div(mul(q(3,2),sum([q(31,99),q(1,99*99),mul(q(35),add(q(1,22500),q(1,400000)))])),sub(q(1),rho));
 ok(gt(q(26),first)&&gt(q(39),quad)&&is(total,q(756,1105))&&gt(q(685,1000),total)&&gt(q(628,1000),low)&&gt(q(684,1000),high),'exact general and older piecewise constants');
 return {first:show(first),quadratic:show(quad),rho:show(rho),correction:show(total),olderLow:show(low),olderHigh:show(high)};
}
const A=setup(), plus=multiple(A,q(1)),minus=multiple(A,q(1,8));
const out={constants:constants(),mixedSeed:A.stats,positiveCenter:plus,negativeCenter:minus,
  padded:{n:1000,sameFinitePattern:true,theta:show(add(q(1000+c2(999)),q(31,36)))}}; 
ok(plus.MMinusK==='31/36'&&minus.MMinusK.startsWith('-'),'opposite exact center shifts');
console.log(JSON.stringify(out,null,2));
