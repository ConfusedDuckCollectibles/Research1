'use strict';
// Independent exact rational check for the fixed-density masked inverse. No imports.
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
 const n=A.length,h=A.map((row,i)=>[...row.map(v=>q(v)),b[i]]);
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
function small(){
 const n=12, K=c2(n-1), a=c2(n-2), q0=c2(n-4), vert=N(n);
 const del=[];
 for(const shift of [0,6])for(const i of [0,1])for(const j of [2,3])for(const k of [4,5]){
  if(i===0&&j===2&&k===4)continue;
  del.push(idxset([shift+i,shift+j,shift+k]));
 }
 const S=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++)
  if(!del.includes(idxset([i,j,k])))S.push([i,j,k]);
 const G=gram(S,n),d=G.map((row,i)=>K-row[i]),d0=Math.min(...d), M=K-d0;
 const w=vert.map(i=>i===1?q(1):i===3?q(1,2):q(0)); // vertex 1 tight; vertex 3 positive but slack
 const allowed=vert.filter(i=>w[i][0]>0n),slacks=vert.filter(i=>gt(q(d[i]-d0),w[i]));
 const lambda=vert.map(i=>div(q(K-d[i]),sub(q(M),w[i]))), L=sub(q(M),sum(w));
 ok(d0===3&&lambda[1][0]===lambda[1][1]&&gt(q(1),lambda[3]),'all-touched tight + positive-slack mixture');
 const C=new Map(),V=new Map();
 function fit(i,bridge){
  const row=vert.filter(j=>j!==i), GG=gram(S.filter(U=>!U.includes(i)),n), mat=row.map(j=>row.map(k=>GG[j][k]));
  const rhs=row.map(j=>bridge?sub(lambda[j],q(GG[j][j],a)):q(-G[i][j]));
  const fit=solve(mat,rhs), u=vert.map(j=>j===i?q(0):fit[row.indexOf(j)]);
  return S.map(U=>U.includes(i)?q(bridge?0:1):add(bridge?q(1,a):q(0),sum(U.map(j=>u[j]))));
 }
 for(const i of allowed)C.set(i,fit(i,true));
 for(const i of slacks)V.set(i,fit(i,false)); // includes positive-weight vertex 3
 const incid=z=>vert.map(i=>sum(z.filter((_,j)=>S[j].includes(i))));
 const CI=new Map(allowed.map(i=>[i,incid(C.get(i))])), VI=new Map(slacks.map(i=>[i,incid(V.get(i))]));
 let directFits=0;
 for(const i of allowed)for(const j of vert){ok(is(CI.get(i)[j],i===j?q(0):lambda[j]),'C fit');directFits++;}
 for(const i of slacks)for(const j of vert){ok(is(VI.get(i)[j],i===j?q(K-d[i]):q(0)),'V fit, including bridge slack');directFits++;}
 const g=S.map((U,i)=>sub(q(1),sum(allowed.map(j=>mul(C.get(j)[i],w[j]))))), gx=incid(g);
 const f=new Map(slacks.map(i=>[i,div(sub(q(1),lambda[i]),q(K-d[i]))]));
 let weightedIncidences=0, locked=0;
 for(let U=0;U<S.length;U++)for(const i of vert){
  const val=sum([div(mul(g[U],gx[i]),L),...allowed.map(j=>mul(w[j],mul(C.get(j)[U],CI.get(j)[i]))),
    ...slacks.map(j=>mul(f.get(j),mul(V.get(j)[U],VI.get(j)[i])))]);
  const target=sum([lambda[i],C.has(i)?mul(q(-1),mul(mul(lambda[i],w[i]),C.get(i)[U])):q(0),
      V.has(i)?mul(sub(q(1),lambda[i]),V.get(i)[U]):q(0)]);
  ok(is(gx[i],mul(L,lambda[i]))&&is(val,target),'independent coupled weighted fit');
  if(S[U].includes(i)){ok(is(val,q(1)),'incident lock');locked++;}
  weightedIncidences++;
 }
 const B=(U,W)=>sum([div(mul(g[U],g[W]),L),...allowed.map(j=>mul(w[j],mul(C.get(j)[U],C.get(j)[W]))),
   ...slacks.map(j=>mul(f.get(j),mul(V.get(j)[U],V.get(j)[W])))]);
 const Ginv=vert.map(i=>solve(G,vert.map(j=>q(Number(i===j)))));
 const den=Ginv.flat().reduce((s,x)=>lcm(s,x[1]),1n);
 const J=Ginv.map(row=>row.map(x=>x[0]*(den/x[1])));
 const R=S.map(U=>vert.map(j=>U.reduce((s,i)=>s+J[i][j],0n)));
 const mask=(i,j)=>S[i].some(v=>S[j].includes(v));
 const Pi=(i,j)=>S[j].reduce((s,v)=>s+R[i][v],0n);
 const Pden=BigInt(K*q0), Pnum=(i,j)=>(i===j?Pden:0n)+(mask(i,j)?0n:BigInt(K))-BigInt(q0);
 const size=S.length, rows=N(size);
 const PSX=rows.map(i=>{const s=vert.map(()=>0n);for(const j of rows){const v=Pnum(i,j);for(const t of S[j])s[t]+=v;}return s;});
 const XPX=vert.map(()=>vert.map(()=>0n));for(const i of rows)for(const t of S[i])for(const j of vert)XPX[t][j]+=PSX[i][j];
 const middle=R.map(r=>XPX.map(x=>dotB(x,r)));
 const allDen=Pden*den*den, p0=(i,j)=>Pnum(i,j)*den*den-den*(dotB(R[i],PSX[j])+dotB(R[j],PSX[i]))+dotB(R[i],middle[j]);
 const offset=sub(q(1,M),q(1,K));
 const batch=[[0,1,2],[0,3,6],[0,4,5],[0,7,8]].map(U=>S.findIndex(T=>idxset(U)===idxset(T)));
 ok(batch.every(i=>i>=0),'compatible signed tile support');
 const aVec=rows.map(i=>BigInt(Number(i===batch[0])-Number(i===batch[2])));
 const bVec=rows.map(i=>BigInt(Number(i===batch[1])-Number(i===batch[3])));
 const ap=rows.map(i=>Pi(i,batch[0])-Pi(i,batch[2])),bp=rows.map(i=>Pi(i,batch[1])-Pi(i,batch[3]));
 const ea=rows.map(i=>den*aVec[i]-ap[i]), eb=rows.map(i=>den*bVec[i]-bp[i]);
 const zeta=sub(q(1),q(0));
 let PiCol=0n,PiQ=0n, PiQPi=0n,contractNum=0n,weightedDefect=q(0),offsetFailures=0,maskRows=0;
 for(const i of rows){
  const ic=vert.map(()=>0n),u=vert.map(()=>0n),v=vert.map(()=>0n),w=vert.map(()=>q(0));let col=0n;
  for(const j of rows){
   const n0=p0(i,j),l=ap[i]*bVec[j]+bp[i]*aVec[j],r=ap[i]*bp[j]+bp[i]*ap[j];
   col+=absB(Pi(i,j)); if(absB(l)>PiQ)PiQ=absB(l);if(absB(r)>PiQPi)PiQPi=absB(r);
   for(const t of S[j])ic[t]+=n0;
   const raw=aVec[i]*bVec[j]+bVec[i]*aVec[j];
   if(!mask(i,j)){ok(raw===0n,'signed tile is locked');continue;}
   const tran=ea[i]*eb[j]+eb[i]*ea[j],d=absB(raw*den*den-tran);
   if(d>contractNum)contractNum=d;
   const err=sub(sub(i===j?zeta:q(0),B(i,j)),q(n0,allDen));
   if(gt(abs(err),weightedDefect))weightedDefect=abs(err);
   for(const t of S[j]){u[t]+=raw;v[t]+=tran;w[t]=add(w[t],err);}
  }
  if(col>PiCol)PiCol=col;
  ok(ic.every(z=>z===0n),'actual projected Kneser block kills incidence');
  for(const t of S[i]){
   ok(u[t]===0n&&v[t]===0n&&w[t][0]===0n,'Q, T(Q), and actual Q0 compatible exactly');
   ok(mul(offset,q(G[t][t]))[0]!==0n,'offset omission breaks each lock');offsetFailures++;maskRows++;
  }
 }
 const mC=q(contractNum,den*den),pc=q(PiCol,den),pl=q(PiQ,den),px=q(PiQPi,den*den);
 ok(gt(q(1),mC)&&gt(mC,q(0)),'nontrivial actual masked correction below one on test tile');
 ok(!gt(px,mul(pc,pl))&&!gt(mC,add(mul(q(2),pl),px)),'direct max and column estimates');
 return {n,m:del.length,minDegree:d0,survivors:size,directFits,weightedIncidences,locked,
   maskRows,offsetFailures,offset:show(offset),weightedDefect:show(weightedDefect),maxPiQ:show(pl),maxPiQPi:show(px),
   maxPiColumnSum:show(pc),normIminusTOnUnitTile:show(mC)};
}
function constants(){
 const a=q(99,200),f=q(12,25),n=q(1000),D=q(1,400),K=q(1,2);
 const c3coef=div(q(1),mul(a,f)),c4coef=add(c3coef,div(q(6),sq(f)));
 const Vconst=div(q(6),f),rho=div(mul(q(18),add(q(2),div(mul(q(6),K),f))),mul(f,n));
 const R=add(div(q(6),sq(f)),add(div(add(q(50),div(q(169),sq(f))),n),mul(q(1024),D)));
 const A=sum([div(q(3),sq(a)),div(q(78),mul(a,n)),div(mul(q(36),D),mul(sq(a),f))]);
 const U=add(div(mul(q(13,12),add(q(6),div(q(49,3),n))),sq(a)),div(mul(q(13*36*36),sq(D)),sq(mul(a,f))));
 const offset=div(q(1),mul(a,f));
 const analytic=div(mul(q(3*123,800),q(1)),sub(q(1),rho));
 for(const [x,y] of [[c3coef,5],[c4coef,32],[Vconst,13],[R,30],[A,14],[U,30],[offset,5]])ok(gt(q(y),x),'coefficient rounded outward');
 ok(is(rho,q(99,320))&&gt(q(123,184),analytic)&&gt(q(99,100),q(123,184)),'strict real spectral bound');
 return {c3coef:show(c3coef),c4coef:show(c4coef),Vconst:show(Vconst),R:show(R),A:show(A),U:show(U),offset:show(offset),rho:show(rho),analytic:show(analytic),reportedLooseCorrection:show(q(123,184))};
}
function productControl(){
 // Dense exact blow-up: at the boundary both the degree and crossing budget saturate.
 const s=17,n=60*s,full=[];
 for(let a=0;a<20;a++)for(let i=0;i<3;i++)for(let j=0;j<3;j++)full.push([a,20+(a+i)%20,40+(a+j)%20]);
 ok(new Set(full.map(idxset)).size===180,'simple block triples');
 const prod=full.slice(1),degrees=N(60).map(i=>prod.filter(t=>t.includes(i)).length*s*s);
 const d0=Math.min(...degrees),D=Math.max(...degrees),m=prod.length*s*s*s;
 const r=s*s,nA=450,alpha=nA*r,eligible=degrees.filter(x=>x>d0).length*s;
 ok(n>=1000&&400*D===n*n&&8*alpha===n*n&&eligible>=nA&&d0===8*r,'saturated large and entirely touched');
 const K=c2(n-1),M=K-d0,a=c2(n-2),S0=n+K,b=r;
 const u=2,v=3*r,crossing=(1+b)*S0-d0;
 const theta=S0*(u+v)-d0*u+Math.max(0,v-b*u);
 const corr=div(mul(q(3*K),q(123*D,n**4)),sub(q(1),q(99,320)));
 ok(M<a&&gt(sub(q(1),q(1,c2(n-4))),corr),'absolute baseline gap and strict positivity at equality in both allowances');
 return {n,D,minDegree:d0,m,totalTripleCount:c3(n),deletedFraction:show(q(m,c3(n))),fractionOfDegreeBudget:show(q(3*m,D*n)),eligibleAnchors:eligible,
   actualAnchors:nA,b,alpha,S0,Ma:a-M,crossing,u,v,theta,certifiedCorrectionUpper:flo(corr)};
}
console.log(JSON.stringify({constants:constants(),small:small(),product:productControl()},null,2));
