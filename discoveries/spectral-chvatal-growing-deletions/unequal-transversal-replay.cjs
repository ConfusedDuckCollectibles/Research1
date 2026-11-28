'use strict';
// Exact rational controls for unbalanced top transversals; portable Node, no imports.
const ok=(x,label)=>{if(!x)throw Error(label);};
const absB=x=>x<0n?-x:x;
function gcd(a,b){while(b)[a,b]=[b,a%b];return absB(a);}
const q=(a,b=1n)=>{a=BigInt(a);b=BigInt(b);ok(b!==0n,'rational denominator');const g=gcd(a,b)*(b<0n?-1n:1n);return [a/g,b/g];};
const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),sub=(a,b)=>add(a,[-b[0],b[1]]);
const mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]);
const sum=a=>a.reduce(add,q(0)),gt=(a,b)=>a[0]*b[1]>b[0]*a[1];
const is=(a,b)=>a[0]===b[0]&&a[1]===b[1];
const show=x=>x.map(String).join('/');
const N=n=>Array.from({length:n},(_,i)=>i);
function solve(A,b){ // exact rational Gauss-Jordan, independent of modular replay
 const n=A.length,h=A.map((row,i)=>[...row.map(v=>Array.isArray(v)?v:q(v)),b[i]]);
 for(let i=0;i<n;i++){
  let p=i;while(p<n&&h[p][i][0]===0n)p++;ok(p<n,'nonsingular incidence fit');[h[i],h[p]]=[h[p],h[i]];
  const z=h[i][i];for(let k=i;k<=n;k++)h[i][k]=div(h[i][k],z);
  for(let j=0;j<n;j++)if(j!==i&&h[j][i][0]!==0n){const y=h[j][i];for(let k=i;k<=n;k++)h[j][k]=sub(h[j][k],mul(y,h[i][k]));}
 }
 return h.map(a=>a[n]);
}

function factorPSD(src){
 const H=src.map(x=>x.slice()),n=H.length;let rank=0;
 for(let i=0;i<n;i++){
   const d=H[i][i];ok(d[0]>=0n,'exact real PSD pivot');
   if(d[0]===0n){ok(H[i].slice(i).every(x=>x[0]===0n),'zero pivot row vanishes');continue;}
   rank++;
   for(let j=i+1;j<n;j++)for(let k=j;k<n;k++) H[j][k]=H[k][j]=sub(H[j][k],div(mul(H[j][i],H[k][i]),d));
 }return rank;
}

function cart(s,r){if(!r)return [[]];return N(s).flatMap(j=>cart(s,r-1).map(x=>[j,...x]));}

const z=q(0), one=q(1),neg=a=>[-a[0],a[1]], abs=a=>a[0]<0n?neg(a):a;
const eq=(A,B)=>A.length===B.length&&A.every((r,i)=>r.length===B[i].length&&r.every((x,j)=>is(x,B[i][j])));
const Gm=(n,k,a)=>N(n).map(i=>N(k).map(j=>a(i,j))),diag=a=>Gm(a.length,a.length,(i,j)=>i===j?a[i]:z);
const at=A=>A[0].map((_,j)=>A.map(R=>R[j])), mplus=(A,B)=>A.map((x,i)=>x.map((a,j)=>add(a,B[i][j]))),
 mminus=(A,B)=>A.map((x,i)=>x.map((a,j)=>sub(a,B[i][j])));
function mm(A,B){const T=at(B);return A.map(a=>T.map(b=>dot(a,b)));}
const dot=(a,b)=>sum(a.map((v,i)=>mul(v,b[i]))),ey=n=>diag(N(n).map(()=>one));
const scale=(A,w)=>A.map(row=>row.map(x=>mul(x,w)));
const flatM=(A,ok_)=>A.every(row=>row.every(ok_)),bound=(A)=>A.reduce((a,b)=>gt(a,b)?a:b,z);
function completion(r,s,p){
 const U=cart(s,r),m=r*s,px=U.map(t=>t.map((v,i)=>i*s+v));
 const deficits=N(m).map(i=>sub(q(s**(r-1)),sum(p.filter((_,k)=>px[k].includes(i)))));
 // exact minimum, not a positive norm
 const mn=deficits.reduce((a,b)=>gt(a,b)?b:a),b=deficits.map(t=>sub(t,mn)),T=sum(b.slice(0,s));
 const pointMass=q(s**(r-1));
 ok(N(r).every(i=>is(sum(b.slice(i*s,(i+1)*s)),T)),"same balancing mass by part");
 const powT=N(r-1).reduce(v=>mul(v,T),one);
 const c=px.map(ik=>T[0]===0n?z:div(ik.reduce((z,i)=>mul(z,b[i]),one),powT));
 for(const i of N(m))ok(is(sum(c.filter((_,j)=>px[j].includes(i))),b[i]),"all completion marginals exact");
 const P=p.map((x,j)=>add(x,c[j])),M=sub(pointMass,mn);
 ok(px[0].length===r && N(m).every(i=>is(sum(P.filter((_,j)=>px[j].includes(i))),M)),"computed constant point marginals");
 const D=bound(deficits),Ds=N(m).map(i=>sum(P.flatMap((v,j)=>px[j].includes(i)?[abs(sub(one,v))]:[]))), Da=bound(Ds);
 ok(!gt(Da,mul(q(2),D)),"absolute signed deficit");
 return {U,px,p,P,c,D,T,mn,M,Ds,Da,deficits};
}
function pseudo(r,s,G){
 const n=r*s,Z=Gm(n,n,(i,j)=>q((Math.floor(i/s)===Math.floor(j/s)?r:0)-1,n)),
  C=mplus(G,Z),Gi=N(n).map(i=>solve(C,N(n).map(j=>q(Number(j===i)))));
 const inv=mminus(Gi,Z);
 ok(eq(mm(G,inv),mminus(ey(n),Z))&&eq(mm(inv,Z),Gm(n,n,()=>z)),"inverse on part-contrast quotient");
 return inv;
}
function geometry(r,s,p,full=false){
 const C=completion(r,s,p),{px,P,D,T,mn,M,Da,c,U}=C,n=r*s,B=U.length,K=s**(r-1),Q=(s-1)**(r-1),
 L=Gm(B,B,(i,j)=>sub(q(1,K),px[i].some(v=>px[j].includes(v))?z:q(1,Q)));
 const X=Gm(B,n,(i,j)=>q(Number(px[i].includes(j)))),Xtr=at(X);
 const H0=mm(mm(X,pseudo(r,s,mm(Xtr,X))),Xtr),
 H=Gm(B,B,(i,j)=>q(s*px[i].filter(v=>px[j].includes(v)).length-(r-1),s*K));
 ok(eq(H,H0)&&eq(mm(H,H),H),"overlap formula equals direct pseudoinverse projector");
 ok(eq(mm(L,X),X)&&eq(mm(H,X),X),"seed and projection fix incidence");
 const BB=mminus(L,H);
 ok(flatM(mm(BB,X),x=>x[0]===0n),"signed overlap kernel kills incidence");
 const schur=BB.map(row=>sum(row.map((x,j)=>mul(abs(x),P[j]))));
 const mSchur=bound(schur);
 const I=N(B).filter(i=>P[i][0]>0n),w=I.map(i=>P[i]),XS=I.map(i=>X[i]),XST=at(XS),W=diag(w);
 const G=mm(mm(XST,W),XS),Gp=pseudo(r,s,G),Nw=mm(mm(XS,Gp),XST),Pi=mm(Nw,W),E=mminus(ey(I.length),Pi);
 ok(eq(mm(Pi,Pi),Pi)&&eq(mm(Pi,XS),XS),"relative projector after singular weighted incidence");
 const Ls=I.map(i=>I.map(j=>L[i][j])),R=mminus(diag(w.map(v=>div(one,v))),Ls),Et=at(E);
 const newSeed=mm(mm(E,R),Et);
 const Delta=diag(P.map(v=>sub(one,v))),DX=mm(Delta,X),LDX=mm(L,DX);
 ok(eq(mm(R,mm(W,XS)),I.map(i=>LDX[i])),"first-order signed seed including omitted zero rows");
 const second=mm(mm(mm(XST,W),R),mm(W,XS));
 ok(eq(second,mminus(mm(Xtr,DX),mm(at(DX),LDX))),"second-order signed seed");
 let negativeRay=null;
 for(let a=0;a<n&&!negativeRay;a++)for(let b=0;b<a&&!negativeRay;b++)for(const s of [-1,1]){
   const v=add(add(second[a][a],second[b][b]),mul(q(2*s),second[a][b]));
   if(v[0]<0n)negativeRay={indexedPoints:[a,b],sign:s,form:show(v)};
 }
 const mask=(i,j)=>px[I[i]].some(u=>px[I[j]].includes(u)),sI=I.length;
 const lock=Gm(sI,sI,(i,j)=>i===j?sub(div(one,w[i]),div(one,M)):neg(div(one,M)));
 const err=Gm(sI,sI,(i,j)=>mask(i,j)?sub(lock[i][j],newSeed[i][j]):z),WXS=mm(W,XS);
 const Trel=Q=>mm(mm(E,Q),Et).map((a,i)=>a.map((x,j)=>mask(i,j)?x:z));
 function compatible(Q){const a=mm(Q,WXS);return I.every((v,i)=>px[v].every(k=>a[i][k][0]===0n));}
 ok(compatible(err)&&compatible(Trel(err)),"relative compatible restriction; first iterate");
 const col=N(sI).map(j=>sum(N(sI).map(i=>abs(Pi[j][i]))));
 const result={r,s,fullTupleCount:B,completedSupport:sI,oldZeroRows:p.filter(v=>v[0]===0n).length,completedZeroRows:B-sI,
  deficit:show(D),minimumDeficit:show(mn),absoluteSignedDeficit:show(Da),addedTotal:show(T),
  marginals:show(M),maxIndividualWeight:show(bound(P)),completionNonzero:c.filter(v=>v[0]!==0n).length,
  overlapProjectorEntries:B*B,annihilationEntries:B*n,weightedOrthogonalProjectorEntries:sI*sI,
  firstOrder:sI*n,secondOrder:n*n,compatibilityChecks:2*sI*r,
  exactSignedSchurNormBound:show(mSchur),relativeProjectorMaxWeightedColumn:show(bound(col)),negativeIncidenceRay:negativeRay
 };
 if(!full)return result;
 const xx=[];
 for(let i=0;i<sI;i++)for(let j=i;j<sI;j++)if(mask(i,j))xx.push([i,j]);
 const aug=[];
 for(const [u,v] of xx){
  aug.push(xx.map(([i,j])=>add(mul(E[u][i],E[v][j]),i!==j?mul(E[u][j],E[v][i]):z)).concat([err[u][v]]));
 }
 for(const u of N(sI))for(const i of px[I[u]]){
  aug.push(xx.map(([a,b])=>add(a===u?WXS[b][i]:z,b===u&&a!==b?WXS[a][i]:z)).concat([z]));
 }
 const sol=linear(aug,xx.length);
 const Qstar=Gm(sI,sI,(i,j)=>mask(i,j)?sol.x[xx.findIndex(([a,b])=>a===Math.min(i,j)&&b===Math.max(i,j))]:z);
 ok(compatible(Qstar)&&eq(Trel(Qstar),err),"complete exact compatible inverse solve");
 const PP=mplus(newSeed,mm(mm(E,Qstar),Et));
 ok(flatM(mm(PP,WXS),x=>x[0]===0n),"corrected matrix kills weighted incidence");
 ok(eq(PP,at(PP))&&PP.every((a,i)=>a.every((x,j)=>!mask(i,j)||is(x,lock[i][j]))),"complete transformed locks");
 const exactRank=factorPSD(PP);
 const oriented=I.filter(i=>p[i][0]>0n),map=oriented.map(i=>I.indexOf(i));
 const orig=Gm(map.length,map.length,(u,v)=>add(PP[map[u]][map[v]],u===v?sub(div(one,p[oriented[u]]),div(one,P[oriented[u]])):z));
 const originalPSD=factorPSD(orig);
 ok(orig.every((a,i)=>a.every((x,j)=>!mask(map[i],map[j])||is(x,i===j?sub(div(one,p[oriented[i]]),div(one,M)):neg(div(one,M))))),"smaller weight certificate by diagonal addition");
 ok(C.deficits.some(d=>is(d,mn)),"actual old star lower bound");
 return {...result,exactCertificate:{maskedDimensions:xx.length,systemRank:sol.rank,positiveCompletedGramRank:exactRank,originalGramRank:originalPSD,correctionNonzero:xx.filter(([i,j])=>Qstar[i][j][0]!==0n).length}};
}
function linear(S,n){
 const A=S.map(r=>r.slice()),idx=[];let r=0;
 for(let c=0;c<n;c++){
  const i=A.findIndex((x,k)=>k>=r&&x[c][0]!==0n);
  if(i<0)continue;
  [A[r],A[i]]=[A[i],A[r]];
  const b=A[r][c];A[r]=A[r].map(x=>div(x,b));
  for(let j=0;j<A.length;j++)if(j!==r&&A[j][c][0]!==0n){
    const b=A[j][c];A[j]=A[j].map((x,k)=>sub(x,mul(b,A[r][k])));
  }
  idx.push(c);r++;
 }
 ok(A.every(x=>x.slice(0,n).some(v=>v[0]!==0n)||x[n][0]===0n),"exact compatible system consistency");
 const x=N(n).map(()=>z);idx.forEach((d,i)=>x[d]=A[i][n]);
 return {x,rank:r};
}
const T3=cart(3,3),bad=new Map([["001",q(1,4)],["010",q(1,2)],["100",z]]);
const hetero=T3.map(U=>U.every(x=>x===U[0])?z:bad.get(U.join(""))??one);
const pure=T3.map(U=>U.join("")==="000"||U.join("")==="111"?z:one);
const tiny=cart(2,3).map((a,i)=>i===3?z:q(4-Math.min(i,7-i),4));
const balanced=cart(2,3).map(()=>q(5,8));
function horizon(r,s,step){
 const R=BigInt(r),S=BigInt(s),ell=S/(BigInt(step)*R**3n),K=S**(R-1n),
  D=ell*(S-1n)**(R-2n),m=ell*(R-1n)*(S-1n)**(R-3n),M=K-m,c=D-m;
 const old=16n*2n**R*D>K;
 ok(D>m&&D*40n*R**3n<=K&&s>=400*r**3,"concentrated construction sits in theorem");
 return {r,s,step,numberResidues:String(ell),D:String(D),m:String(m),M:String(M),largestCompleted:String(1n+c),outsideOldExponentialCap:old};
}
function constants(r){
 const R=BigInt(r),s=400n*R**3n,K=s**(R-1n),d=q(1n,40n*R**3n),da=mul(q(2),d);
 const Q=(s-1n)**(R-1n),rq=q(Q,K),rInv=div(q(3),sub(one,mul(q(3n*R),da)));
 const overlap=add(q(5n*R*(R-1n),2n*s),mul(q(R),d));
 const disjoint=sub(sub(div(one,rq),one),q(R-1n,s));
 const column=mul(q(4n*R),add(q(R-1n,s),d)),rho=mul(q(2n+4n*R),column);
 const norms=div(mul(q(13n*R**3n),d),sub(one,rho));
 ok(!gt(q(99,100),rq)&&!gt(rInv,q(4))&&!gt(disjoint,q(R*(R-1n),s*s))&&
  !gt(rho,q(1,4))&&!gt(norms,q(13,30))&&!gt(overlap,q(1,100)),"cutoff exact inequalities");
 return {r,s:String(s),overlap_bound:show(overlap),compatible_rho:show(rho),correction:show(norms)};
}
const out={heterogeneous:geometry(3,3,hetero),indefiniteControl:geometry(3,3,pure),certificate:geometry(3,2,tiny,true),
 balancedTZero:geometry(3,2,balanced),constants:[3,4,13,14,31,73].map(constants),
 large:[horizon(3,10801,80),horizon(13,900000,45),horizon(14,1097600,80),horizon(21,400*21**3,80)]};
ok(out.large.slice(1).every(x=>x.outsideOldExponentialCap),"new high-rank regime arithmetic");
console.log(JSON.stringify(out,null,2));
