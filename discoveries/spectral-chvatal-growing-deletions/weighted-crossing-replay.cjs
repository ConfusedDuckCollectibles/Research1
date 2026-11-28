'use strict';
// Independent exact check of the budget-free weighted crossing. No imports.
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
// (a,b) encodes a+b sqrt(2); integer-ring arithmetic for large exact matrices.
const Z=[0n,0n], zi=x=>[BigInt(x),0n];
const za=(x,y)=>[x[0]+y[0],x[1]+y[1]],zs=(x,y)=>[x[0]-y[0],x[1]-y[1]];
const zm=(x,y)=>[x[0]*y[0]+2n*x[1]*y[1],x[0]*y[1]+x[1]*y[0]];
const zk=(x,k)=>[x[0]*k,x[1]*k], ze=(x,y)=>x[0]===y[0]&&x[1]===y[1];
const aZ=(x,d=1n)=>[q(x[0],d),q(x[1],d)], AQ=a=>[a,q(0)], QA=i=>AQ(q(i));
const aa=(x,y)=>[add(x[0],y[0]),add(x[1],y[1])], as=(x,y)=>[sub(x[0],y[0]),sub(x[1],y[1])];
const am=(x,y)=>[add(mul(x[0],y[0]),mul(q(2),mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))];
const aeq=(x,y)=>is(x[0],y[0])&&is(x[1],y[1]), sumA=A=>A.reduce(aa,QA(0));
function sign([a,b]){if(a[0]===0n)return Math.sign(Number(b[0]));if(b[0]===0n)return Math.sign(Number(a[0]));
 if(a[0]>0n===b[0]>0n)return a[0]>0n?1:-1;
 const s=gt(sq(a),mul(q(2),sq(b)))?a:b;return s[0]>0n?1:-1;}
const showA=x=>x.map(show);
function exactCrossing(){
 const n=10,K=c2(n-1),q0=c2(n-4),a=c2(n-2), zN=N(n);
 const patch=N(n).filter(i=>i!==0&&i!==9).map(i=>idxset([i,(i+1)%n,(i+3)%n]));
 const S=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++)S.push([i,j,k]);
 const m=S.length, rows=N(m), del=rows.filter(i=>patch.includes(idxset(S[i]))),DS=new Set(del);
 const d=zN.map(i=>del.filter(j=>S[j].includes(i)).length),d0=Math.min(...d);
 // k=1 at three chosen anchors. beta=1/2, f=√2/2 on the eight deleted triples.
 const anchors=[1,4,7], beta=q(1,2), p=rows.map(i=>DS.has(i)?beta:q(1));
 const fn=rows.map(i=>DS.has(i)?[0n,1n]:[2n,0n]), Delta=fn.map(x=>zs([2n,0n],x));
 const f=fn.map(x=>aZ(x,2n)),t=d.map(di=>sub(q(K),mul(q(di),beta))),M=sub(q(K),mul(q(d0),beta));
 const w=zN.map(i=>anchors.includes(i)?beta:q(0)),L=sub(M,sum(w)),lam=zN.map(i=>div(t[i],sub(M,w[i])));
 const slack=zN.filter(i=>gt(q(1),lam[i]));
 ok(d0===1&&Math.min(...anchors.map(i=>d[i]-d0))===1&&slack.includes(4)&&slack.includes(7)&&!slack.includes(1),'exact k=1 tight/slack crossing');
 const wc=(i,j,avoid=-1)=>sum(p.filter((_,r)=>(avoid<0||!S[r].includes(avoid))&&S[r].includes(i)&&S[r].includes(j)));
 const G=zN.map(i=>zN.map(j=>wc(i,j)));
 const C=new Map(), V=new Map();
 function fit(i,bridge){
  const sets=zN.filter(j=>j!==i),GG=sets.map(j=>sets.map(k=>wc(j,k,i)));
  const b=sets.map((j,u)=>bridge?sub(lam[j],div(GG[u][u],q(a))):mul(q(-1),G[i][j]));
  const r=solve(GG,b);
  return S.map(U=>U.includes(i)?q(bridge?0:1):add(bridge?q(1,a):q(0),sum(U.map(v=>r[sets.indexOf(v)]))));
 }
 for(const i of anchors)C.set(i,fit(i,true));
 for(const i of slack)V.set(i,fit(i,false));
 const inc=(v,i)=>sum(rows.filter(r=>S[r].includes(i)).map(r=>mul(v[r],p[r])));
 let incidenceFits=0, normalizedLocks=0;
 for(const [i,c] of C)for(const j of zN){ok(is(inc(c,j),j===i?q(0):lam[j]),'physical weighted C incidence');incidenceFits++;}
 for(const [i,v] of V)for(const j of zN){ok(is(inc(v,j),j===i?t[i]:q(0)),'physical weighted V incidence');incidenceFits++;}
 const g=rows.map(r=>sub(q(1),sum(anchors.map(i=>mul(w[i],C.get(i)[r])))));
 const weight=new Map(slack.map(i=>[i,div(sub(q(1),lam[i]),t[i])]));
 const B=(r,s)=>sum([div(mul(g[r],g[s]),L),...anchors.map(i=>mul(w[i],mul(C.get(i)[r],C.get(i)[s]))),
    ...slack.map(i=>mul(weight.get(i),mul(V.get(i)[r],V.get(i)[s])))]);
 const Bhat=rows.map(r=>rows.map(s=>am(am(f[r],f[s]),AQ(B(r,s)))));
 const tight=zN.filter(i=>is(lam[i],q(1)));
 for(const i of zN){
  ok(is(inc(g,i),mul(L,lam[i])),'weighted g incidence');
  for(const r of rows){
   const v=sumA(rows.filter(s=>S[s].includes(i)).map(s=>am(Bhat[r][s],f[s])));
   const y=am(f[r],addPhysical(r,i));
   ok(aeq(v,y),'direct normalized coupled old incidence');
   if(S[r].includes(i)){ok(aeq(v,f[r]),'weighted normalized lock');normalizedLocks++;}
   if(tight.includes(i)){
    const extra=C.has(i)?am(f[r],AQ(mul(w[i],C.get(i)[r]))):QA(0);
    ok(aeq(aa(v,extra),f[r]),'tight star before physical row rescaling');
   }
  }
  if(tight.includes(i))for(const j of anchors)ok(is(add(inc(C.get(j),i),q(Number(i===j))),q(1)),'tight bridge output');
 }
 function addPhysical(r,i){return AQ(add(mul(lam[i],sub(q(1),C.has(i)?mul(w[i],C.get(i)[r]):q(0))),
     V.has(i)?mul(sub(q(1),lam[i]),V.get(i)[r]):q(0)));}
 const Gi=zN.map(i=>solve(G,zN.map(j=>q(Number(i===j))))),gd=Gi.flat().reduce((d,x)=>lcm(d,x[1]),1n);
 const GiN=Gi.map(v=>v.map(x=>x[0]*(gd/x[1]))), R=S.map(U=>zN.map(i=>U.reduce((s,j)=>s+GiN[j][i],0n)));
 const mask=(i,j)=>S[i].some(x=>S[j].includes(x));
 const pd=BigInt(K*q0), pn=(i,j)=>(i===j?pd:0n)+(mask(i,j)?0n:BigInt(K))-BigInt(q0);
 const pid=4n*gd, Pi=rows.map(i=>rows.map(j=>zk(zm(fn[i],fn[j]),S[j].reduce((x,r)=>x+R[i][r],0n))));
 // Direct dense products for the normalized projection correction, then independent incidence identities.
 const P1=rows.map(i=>rows.map(j=>rows.reduce((x,k)=>za(x,zk(Pi[i][k],pn(k,j))),Z)));
 const P2=rows.map(i=>rows.map(j=>rows.reduce((x,k)=>za(x,zm(P1[i][k],Pi[k][j])),Z)));
 const PXdelta=zN.map(i=>rows.map(j=>del.filter(k=>S[k].includes(i)).reduce((x,k)=>za(x,zk(Delta[k],pn(k,j))),Z)));
 const deltaPXdelta=zN.map(i=>zN.map(j=>del.filter(k=>S[k].includes(j)).reduce((x,k)=>za(x,zm(Delta[k],PXdelta[i][k])),Z)));
 let projectionProducts=0,p0null=0,comp=0,omittedOff=0;
 const p0d=pid*pid*pd, P0=rows.map(i=>rows.map(j=>za(zs(zi(pn(i,j)*pid*pid),zk(za(P1[i][j],P1[j][i]),pid)),P2[i][j])));
 for(const i of rows){
  const p0mom=zN.map(()=>Z),Def=zN.map(()=>QA(0)),Wrong=zN.map(()=>QA(0));
  for(const j of rows){
   const first=zm(fn[i],zN.reduce((x,r)=>zs(x,zk(PXdelta[r][j],R[i][r])),Z));
   const second=zm(zm(fn[i],fn[j]),zN.reduce((x,r)=>zN.reduce((y,s)=>za(y,zk(deltaPXdelta[r][s],R[i][r]*R[j][s])),x),Z));
   ok(ze(first,P1[i][j])&&ze(second,P2[i][j]),'independent direct weighted projection identities');
   projectionProducts++;
   for(const t of S[j])p0mom[t]=za(p0mom[t],zm(P0[i][j],fn[j]));
   if(mask(i,j)){
    const defect=as(as(QA(Number(i===j)),Bhat[i][j]),aZ(P0[i][j],p0d));
    const om=as(am(am(f[i],f[j]),AQ(div(q(1),M))),AQ(q(1,K)));
    for(const t of S[j]){Def[t]=aa(Def[t],am(defect,f[j]));Wrong[t]=aa(Wrong[t],am(aa(defect,om),f[j]));}
   }
  }
  for(const x of p0mom){ok(ze(x,Z),'P0 annihilates every weighted incidence');p0null++;}
  for(const t of S[i]){ok(aeq(Def[t],QA(0)),'actual weighted mask defect is compatible');comp++;if(!aeq(Wrong[t],QA(0)))omittedOff++;}
 }
 ok(omittedOff>0,'weighted offset omission causes genuine incompatibility');
 // A nonvacuous gap test on an exact weighted-complement vector derived from a signed incidence trade.
 const sr=(x,a,b)=>S.findIndex(U=>idxset(U)===idxset([x,a,b]));let trade;
 for(const x of zN)for(const a of zN)for(const b of zN)for(const c of zN)for(const d of zN){
  if(trade||new Set([x,a,b,c,d]).size!==5)continue;
  const v=[sr(x,a,b),sr(x,c,d),sr(x,a,c),sr(x,b,d)];if(v.filter(i=>DS.has(i)).length===1)trade=v;
 }
 ok(trade,'mixed trade');const y=trade.map((i,j)=>am(QA(j<2?1:-1),DS.has(i)?aZ([0n,1n]):QA(1)));
 const norm=sumA(y.map(v=>am(v,v))), Xy=zN.map(k=>sumA(y.filter((_,j)=>S[trade[j]].includes(k))));
 for(const k of zN)ok(aeq(sumA(y.map((v,j)=>S[trade[j]].includes(k)?am(v,f[trade[j]]):QA(0))),QA(0)),'trade is in weighted complement');
 const GF=gram(S,n),GFi=zN.map(i=>solve(GF,zN.map(j=>q(Number(i===j)))));
 const proj=sumA(zN.flatMap(i=>zN.map(j=>am(am(Xy[i],Xy[j]),AQ(GFi[i][j])))));
 const hv=as(QA(1),aZ([0n,1n],2n)), hs=am(hv,hv);
 const major=am(hs,norm), gap=am(AQ(sub(q(1),q(1,q0))),am(as(QA(1),hs),norm));
 const ray=sumA(y.flatMap((v,i)=>y.map((w,j)=>am(am(v,w),AQ(q(pn(trade[i],trade[j]),pd))))));
 ok(sign(proj)>0&&sign(as(major,proj))>=0&&sign(as(ray,gap))>=0,'nontrivial real gap control in Q(sqrt2)');
 return {n,deletions:del.length,degrees:d,anchors,tight,slack,incidenceFits,normalizedLocks,p0null,comp,omittedOff,
    directFirstAndQuadraticProducts:[projectionProducts,projectionProducts],
    gapTrade:trade.map(i=>S[i]),normSquared:showA(norm),unweightedProjectionSquared:showA(proj),gapLowerRay:showA(gap),actualRay:showA(ray)};
}
function constants(){
 const a=q(99,200), l=q(12,25), DN=q(1,400), h=q(8,1000),nn=q(1000);
 const gap=mul(sub(q(1),div(q(1),mul(a,sq(nn)))),sub(q(1),sq(h)));
 const F=mul(div(q(6),l),add(q(1),div(DN,a)));
 const Qc=div(mul(q(72),DN),sq(l)),hh=div(q(2),l),cc=div(DN,mul(a,l));
 const d=sum([mul(q(30),DN),mul(q(2),q(13)),q(1),q(5),q(11,1000)]);
 const rho=q(99,320), corr=div(q(396,1000),sub(q(1),rho)), loose=div(q(396,1000),q(69,100));
 ok(gt(q(13),F)&&is(Qc,q(25,32))&&gt(q(5),hh)&&gt(q(11,1000),cc)&&gt(q(33),d)&&
    gt(gap,q(99,100))&&gt(loose,corr)&&is(loose,q(66,115)),'exact weighted constants');
 return {entryFirst:show(F),entryQuadratic:show(Qc),matrixScaleH:show(hh),matrixScaleGamma:show(cc),
    combinedEntry:show(d),weightedGapLower:show(gap),rho:show(rho),correctionAtExactRho:show(corr),correctionLoose:show(loose)};
}
function large(){
 const s=17,n=60*s,full=[];
 for(let a=0;a<20;a++)for(let i=0;i<3;i++)for(let j=0;j<3;j++)full.push([a,20+(a+i)%20,40+(a+j)%20]);
 const products=full.slice(1),m=products.length*s**3,d=N(60).map(i=>products.filter(x=>x.includes(i)).length*s*s);
 const D=Math.max(...d),d0=Math.min(...d),k=s*s,A=d.filter(x=>x>d0).length*s,K=c2(n-1),S0=n+K;
 const actual=A*k,cap=n*n/8,gamma=q(1,1+k),normalized=mul(q(actual),gamma),effectiveD=mul(q(D),gamma);
 const beta=mul(q(k),gamma), lower=mul(sub(q(1),q(1,c2(n-4))),sub(q(1),sq(gamma)));
 const a0=sub(q(c2(n-2)),mul(q(6),effectiveD)),rho=div(mul(q(18*n),add(q(2),div(q(6*K),a0))),a0);
 const corr=div(mul(q(3*K),div(mul(q(33),gamma),q(n*n))),sub(q(1),rho));
 const Crossing=(1+k)*S0-d0;
 ok(D*400===n*n&&actual>cap&&gt(q(n),normalized)&&is(div(q(actual),q(cap)),q(969,450))&&gt(lower,corr),'analytic above-cap product');
 return {n,m,deletionFraction:show(q(m,c3(n))),D,d0,k,A,actual,oldCap:cap,normalized:show(normalized),effectiveD:show(effectiveD),
    beta:show(beta),crossing:Crossing,markerMass:(A+1)*k,certifiedCorrectionUpper:flo(corr),weightedGapLower:flo(lower)};
}
console.log(JSON.stringify({constants:constants(),weightedSmall:exactCrossing(),aboveBudget:large()},null,2));
