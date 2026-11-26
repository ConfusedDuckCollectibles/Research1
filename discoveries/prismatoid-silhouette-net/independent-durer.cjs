"use strict";
// Independent exact audit: integer point heights and 760-Lipschitz contours;
// full graph angle feasibility by complex fan products, distinct from ear turns.
const P=[[65,0],[63,16],[39,52],[16,63],[-25,60],[-56,33],[-65,0],[-52,-39],[-16,-63],[33,-56]].map(a=>a.map(BigInt));
const n=P.length,M=1000000n,at=(i)=>(i%n+n)%n;
function req(x,m){if(!x)throw Error(m);}
const sub=(a,b)=>a.map((v,i)=>v-b[i]), dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0n);
const det=(a,b)=>a[0]*b[1]-a[1]*b[0], neg=x=>x<0n?-x:x;
const x3=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const pos=(a,b,c)=>det(sub(b,a),sub(c,a)), pr=(a,b)=>a*b>0n;
function root(N){let l=0n,u=1n;while(u*u<=N)u*=2n;while(u-l>1n){const m=(u+l)/2n;if(m*m<=N)l=m;else u=m;}return l;}
function floor(a,b){return a>=0n?a/b:-((-a+b-1n)/b);}
const ear=[];
let crossings=0;
for(let i=0;i<n;i++){
  req(dot(P[i],P[i])===4225n,"radius");
  for(let j=0;j<n;j++)if(j!==i&&j!==at(i+1))req(pos(P[i],P[at(i+1)],P[j])>0n,"strict hull");
  if(i%2===0)for(let j=1;j<n;j+=2){
    const a=P[i],b=P[at(i+2)],c=P[j],d=P[at(j+2)];
    if(pr(-pos(a,b,c),pos(a,b,d))&&pr(-pos(c,d,a),pos(c,d,b)))crossings++;
  }
  const u=sub(P[at(i+1)],P[at(i-1)]), v=sub(P[i],P[at(i-1)]), D=dot(u,u), A=-det(u,v);
  req(D>0n&&A>0n&&dot(u,v)>0n&&dot(u,v)<D,"exterior strict foot");
  ear.push({D,A,out:[u[1],-u[0]]});
}
const obtuse=[];
for(let i=0;i<n;i++)if(dot(sub(P[at(i-1)],P[i]),sub(P[at(i+1)],P[i]))+900n<0n)obtuse.push(i);
req(crossings===10&&obtuse.length,"crossings / obtuse even at h30");
function patch(k,h){
  return P.map((p,i)=>{
    if(i%2===k)return p.map(x=>M*x);
    const {D,A,out}=ear[i],q=root(M*M*(A*A+D*h*h))-M*A;
    // D >= 1, so error in each physical coordinate is < 2/M.
    return p.map((x,a)=>M*x+floor(q*out[a],D));
  });
}
let minDet=null,minAt=null;
for(let h=0n;h<=30n;h++)for(let k=0;k<2;k++){
  const C=patch(k,h);
  for(let j=0;j<8;j++){
    const q=pos(C[8],C[9],C[j]); // approximate determinant scaled M^2
    if(minDet===null||q<minDet){minDet=q;minAt=[Number(h),k,j];}
  }
}
// Coordinate-component errors < 2/M; chord Euclidean-norm errors < 6/M.
// At h <= 30 all developed chord lengths <= 2*65 + 2*30 = 190.
// Thus determinant rounding error < (6/M)*(380+6/M) < 1.
// Each point is 1-Lipschitz in h; |det'| <= 2*190+2*190 = 760.
// Nearest integer h is within 1/2; physical determinant >= minApprox - 381.
req(minDet>381n*M*M,"all-height support");
// Independently recognize facets via all supporting planes at the two test heights.
const expected=[[0,2,4,6,8],[1,3,5,7,9],...P.map((_,i)=>[at(i-1),i,at(i+1)].sort((a,b)=>a-b))].map(f=>f.join()).sort().join(";");
const nb=[-2,-1,1,2], N=P.map((_,i)=>nb.map(d=>at(i+d)));
function configuration(h){
  const X=P.map((p,i)=>[...p,BigInt(i%2)*h]),facets=new Set;
  for(let a=0;a<n;a++)for(let b=a+1;b<n;b++)for(let c=b+1;c<n;c++){
    const w=x3(sub(X[b],X[a]),sub(X[c],X[a]));
    if(w.every(t=>t===0n))continue;
    const d=X.map(x=>dot(w,sub(x,X[a])));
    if(d.every(t=>t>=0n)||d.every(t=>t<=0n)){
      req(w[2]!==0n,"vertical facet");
      facets.add(d.flatMap((t,i)=>t===0n?[i]:[]).join());
    }
  }
  req([...facets].sort().join(";")===expected,"facet census");
  return X;
}
function goodTurns(X){
  // N is in cyclic order around each vertex, with the wrap across its base facet.
  let numGood=0,minCertified=null;
  const keys=(a,b)=>a<b?a+","+b:b+","+a;
  const A=N.map((N0,i)=>{
    const rays=N0.map(j=>sub(X[j],X[i]));
    const f=rays.map((u,j)=>{
      const v=rays[(j+1)%4],s=x3(u,v);
      return [M*dot(u,v),root(M*M*dot(s,s))];
    });
    const fan=(j,k)=>{
      let r=1n,s=0n,ell=0;
      for(;j!==k;j=(j+1)%4){
        const [a,b]=f[j];[r,s]=[r*a-s*b,r*b+s*a];ell++;
      }
      // After division by M the exact factor is cos(theta)+i sin(theta)
      // times positive length product. Each factor has modulus < 34200;
      // |z - z_approx| < 1/M. Telescoping gives the scaled error below.
      const err=BigInt(ell)*(34200n*M)**BigInt(ell-1);
      req(neg(s)>err,"unresolved fan sign");
      const q=(neg(s)-err)/(M**BigInt(ell-1));
      if(minCertified===null||q<minCertified)minCertified=q;
      // Every proper fan summed angle is in (0,2pi), since the full
      // angle at a vertex of a full dimensional convex polytope is < 2pi.
      return s>0n;
    };
    const good=new Set;
    for(let j=0;j<4;j++)for(let k=j+1;k<4;k++){
      const one=fan(j,k), other=fan(k,j);
      if(one&&other)good.add(keys(N0[j],N0[k]));
    }
    numGood+=good.size;
    return good;
  });
  let allCycles=0;const cycles=[];
  function visit(L){
    const last=L.at(-1);
    if(L.length===n){
      if(!N[last].includes(0)||L[1]>last)return;
      allCycles++;
      if(L.every((i,k)=>A[i].has(keys(L[at(k-1)],L[at(k+1)]))))cycles.push(L);
      return;
    }
    for(const j of N[last])if(!L.includes(j))visit([...L,j]);
  }
  visit([0]);
  return {allCycles,numGood,cycles,minCertifiedScaled:minCertified.toString()};
}
const h13=goodTurns(configuration(13n)),h14=goodTurns(configuration(14n));
req(h13.cycles.length===1&&h13.cycles[0].join()==P.map((_,i)=>i).join(),"h=13 only silhouette");
req(h14.cycles.length===0,"h=14 no Hamiltonian quasigeodesic");
console.log(JSON.stringify({crossings,obtuseAt30:obtuse,nIntegerHeightChecks:31,minimumApproxDetScaled:minDet.toString(),detScale:(M*M).toString(),minimumAt:minAt,h13,h14},null,2));
