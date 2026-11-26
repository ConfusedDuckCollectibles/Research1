// Exact example checks for the semiprime cardinality result.
function replaySemiprimeFuglede(log=console.log) {
 const seq=n=>Array.from({length:n},(_,i)=>i);
 const mod=(x,n)=>(x%n+n)%n;
 const gcd=(a,b)=>b?gcd(b,a%b):Math.abs(a);
 const cache=new Map();
 function divide(f,g) {
  const a=f.slice(),d=g.length-1,q=Array(f.length-d).fill(0n);
  for(let i=a.length-1;i>=d;i--){
   const c=q[i-d]=a[i];
   if(c)for(let k=0;k<=d;k++)a[i-d+k]-=c*g[k];
  }
  return [q,a.slice(0,d).every(v=>v===0n)];
 }
 function phi(N) {
  if(cache.has(N))return cache.get(N);
  let f=Array(N+1).fill(0n);f[0]=-1n;f[N]=1n;
  for(let d=1;d<N;d++)if(N%d===0){
   const [q,ok]=divide(f,phi(d));if(!ok)throw Error("Bad factorization");f=q;
  }
  cache.set(N,f);return f;
 }
 function zero(exponents,N){
  const f=Array(N).fill(0n);for(const x of exponents)f[mod(x,N)]++;
  return divide(f,phi(N))[1];
 }
 function orderLevel(x,p,a) {
  x=mod(x,p**a);if(!x)return 0;
  let l=a;while(x%p===0){x/=p;l--;}return l;
 }
 function testCase(c){
  const {p,q,a,b,K,A,B,r,s,label}=c,P=p**a,Q=q**b,N=P*Q*K;
  function zeros(X,u,v,h){
   return zero(X.map(([x,y,z])=>u*x*(N/P)+v*y*(N/Q)+h*(z??0)*(N/K)),N);
  }
  function pure(X,dir){
   return (dir?seq(b):seq(a)).map(i=>i+1).filter(i=>dir?zeros(X,0,q**(b-i),0):zeros(X,p**(a-i),0,0));
  }
  if(A.length!==p*q||B.length!==p*q)throw Error("Cardinality");
  let spectralPairs=0,matchedMixed=0,largerFirst=0,largerSecond=0;
  for(let i=0;i<B.length;i++)for(let j=0;j<i;j++){
   const d=B[i].map((x,k)=>x-(B[j][k]??0));
   if(!zeros(A,d[0],d[1],d[2]??0))throw Error("Nonorthogonal pair "+label);
   spectralPairs++;
   if(mod(d[2]??0,K)===0){
    const u=orderLevel(d[0],p,a),v=orderLevel(d[1],q,b);
    if(u!==r&&v!==s)throw Error("Pair misses required orders");
    if(u===r&&v===s)matchedMixed++;
    if(u>r&&v===s)largerFirst++;
    if(u===r&&v>s)largerSecond++;
   }
  }
  if(pure(A,0).join()!==String(r)||pure(A,1).join()!==String(s)||!zeros(A,p**(a-r),q**(b-s),0))throw Error("Pure or mixed root missing");
  const mirror=[pure(B,0),pure(B,1)];
  if(K===1&&(mirror[0].join()!==String(a-r+1)||mirror[1].join()!==String(b-s+1)))throw Error("Mirror");
  function digits(prime,count,omitted){
   return seq(prime**count).filter(x=>Math.floor(x/prime**omitted)%prime===0);
  }
  const X=digits(p,a,r-1),Y=digits(q,b,s-1),Z=seq(K),covered=new Set();
  for(const [x,y,z=0] of A)for(const u of X)for(const v of Y)for(const w of Z){
   const key=mod(x+u,P)+P*mod(y+v,Q)+P*Q*mod(z+w,K);
   if(covered.has(key))throw Error("Overlap");covered.add(key);
  }
  if(covered.size!==N)throw Error("Gap");
  const counts=seq(K).map(h=>B.filter(c=>mod(c[2]??0,K)===h).length);
  function spectrumCensus(){
   if(K!==1)return null;
   const Lp=Math.max(0,...seq(a).map(i=>i+1).filter(l=>l!==r&&zeros(A,p**(a-l),q**(b-s),0)));
   const Lq=Math.max(0,...seq(b).map(i=>i+1).filter(l=>l!==s&&zeros(A,p**(a-r),q**(b-l),0)));
   const bp=BigInt(p),bq=BigInt(q);
   const fixed=bp**BigInt((p-1)*(r-1)+p*(q-1)*Math.min(Lp,r-1))
     *bq**BigInt((q-1)*(s-1)+q*(p-1)*Math.min(Lq,s-1));
   const count=fixed*(bp**BigInt((q-1)*Math.max(Lp-r,0))+bq**BigInt((p-1)*Math.max(Lq-s,0))-1n);
   const available=seq(P*Q).map(t=>[t%P,(t/P)|0]);
   const root=available.map(([x,y])=>zeros(A,x,y,0));
   function ok(i,j){
    return root[mod(available[i][0]-available[j][0],P)+P*mod(available[i][1]-available[j][1],Q)];
   }
   let exact=0n;
   function visit(depth,C){
    if(depth===p*q){exact++;return;}
    for(let i=0;i<C.length&&C.length-i>=p*q-depth;i++)visit(depth+1,C.slice(i+1).filter(j=>ok(C[i],j)));
   }
   visit(1,seq(P*Q).filter(i=>i&&root[i]));
   if(count!==exact)throw Error("Spectrum count mismatch");
   return {Lp,Lq,normalized:String(count)};
  }
  return {label,N,cardinality:A.length,spectralPairs,r,s,mirror,matchedMixed,largerFirst,largerSecond,counts,complementSize:X.length*Y.length*K,spectra:spectrumCensus()};
 }
 const product=(p,q,fn)=>seq(p).flatMap(i=>seq(q).map(j=>fn(i,j)));
 const plain=(p,q,r,s)=>product(p,q,(i,j)=>[p**(r-1)*i,q**(s-1)*j]);
 const cases=[
  {label:"Cartesian",p:2,q:3,a:2,b:2,K:1,r:2,s:1,A:plain(2,3,2,1),B:plain(2,3,1,2)},
  {label:"Larger q-order only",p:2,q:3,a:2,b:2,K:1,r:1,s:1,A:plain(2,3,1,1),B:product(2,3,(i,j)=>[2*i,i+3*j])},
  {label:"Larger p-order only",p:2,q:3,a:3,b:2,K:1,r:1,s:1,A:plain(2,3,1,1),B:product(2,3,(i,j)=>[j+4*i,3*j])},
  {label:"A gap below matched root",p:2,q:3,a:2,b:2,K:1,r:2,s:2,A:product(2,3,(i,j)=>[2*i,i+3*j]),B:plain(2,3,1,1)},
  {label:"Odd pq with redundant binary factor",p:3,q:5,a:2,b:2,K:2,r:1,s:1,A:product(3,5,(i,j)=>[i,j,i%2]),B:product(3,5,(i,j)=>[3*i,5*j,j%2])},
  {label:"Unused odd 25-coordinate",p:2,q:3,a:2,b:2,K:25,r:2,s:1,A:product(2,3,(i,j)=>[2*i,j,4*i]),B:product(2,3,(i,j)=>[i,3*j,2*j])},
  {label:"Separated primes and binary factor",p:3,q:7,a:2,b:2,K:2,r:1,s:1,A:product(3,7,(i,j)=>[i,j,i%2]),B:product(3,7,(i,j)=>[3*i,7*j,j%2])},
  {label:"Odd pq with two binary digits",p:5,q:7,a:1,b:2,K:4,r:1,s:1,A:product(5,7,(i,j)=>[i,j,i%4]),B:product(5,7,(i,j)=>[i,7*j,j%4])},
 ];
 const tested=cases.map(testCase);
 if(tested[1].matchedMixed||!tested[1].largerSecond||tested[2].matchedMixed||!tested[2].largerFirst)throw Error("Downward closure controls");
 function gridCensus(){
  const p=2,q=3,a=2,b=2,P=4,Q=9,m=6;
  const G=seq(P*Q).map(h=>[h%P,(h/P)|0]),output=[];
  for(let r=1;r<=a;r++)for(let s=1;s<=b;s++){
   const chosen=[0];let count=0,asymmetric=0;
   function ok(i,j){
    return orderLevel(G[i][0]-G[j][0],p,a)===r||orderLevel(G[i][1]-G[j][1],q,b)===s;
   }
   function visit(C){
    if(chosen.length===m){
     count++;
     const B=chosen.map(k=>G[k]),X=new Set(B.map(([x,y])=>x%p**(a-r))),Y=new Set(B.map(([x,y])=>y%q**(b-s)));
     if(X.size>1&&Y.size>1)throw Error("Incompatible varying prefixes");
     if(X.size>1||Y.size>1)asymmetric++;
     if(!zero(B.map(([x])=>x),p**(a-r+1))||!zero(B.map(([,y])=>y),q**(b-s+1)))throw Error("Mirror check");
     return;
    }
    for(let i=0;i<C.length&&C.length-i>=m-chosen.length;i++){
     chosen.push(C[i]);visit(C.slice(i+1).filter(j=>ok(C[i],j)));chosen.pop();
    }
   }
   visit(seq(P*Q).filter(i=>i&&ok(i,0)));
   output.push({r,s,count,asymmetric});
  }
  return output;
 }
 const ambiguous=product(3,5,(i,j)=>[i+3*(j%3),i+5*j,j>=3?1:0]);
 const candidates=[];
 for(let r=1;r<=2;r++)for(let s=1;s<=2;s++){
  const good=ambiguous.every(([x,y,z],i)=>ambiguous.slice(0,i).every(([u,v,w])=>
   z!==w||orderLevel(x-u,3,2)===r||orderLevel(y-v,5,2)===s));
  if(good)candidates.push([r,s]);
 }
 if(JSON.stringify(candidates)!==JSON.stringify([[1,2],[2,1]]))throw Error("Ambiguous cohort example");
 const result={tested,gridCensus:gridCensus(),ambiguousCohortNecessaryOrders:candidates};
 log(result);return result;
}
replaySemiprimeFuglede();

// Independent exhaustive root-mask and anchored-clique census at orders 24 and 36.
function exhaustiveSemiprimeSmall(log=console.log) {
 const seq=n=>Array.from({length:n},(_,i)=>i);
 const gcd=(a,b)=>b?gcd(b,a%b):Math.abs(a);
 const known=new Map();
 function phi(n){
  if(known.has(n))return known.get(n);
  let poly=Array(n+1).fill(0);poly[0]=-1;poly[n]=1;
  for(let d=1;d<n;d++)if(n%d===0){
   const h=phi(d),q=Array(poly.length-h.length+1).fill(0);
   for(let t=poly.length-1;t>=h.length-1;t--){
    const m=q[t-h.length+1]=poly[t];
    for(let j=0;j<h.length;j++)poly[t-h.length+1+j]-=m*h[j];
   }
   if(poly.some(x=>x))throw Error("Cyclotomic polynomial");poly=q;
  }
  known.set(n,poly);return poly;
 }
 function check(p,q,a,b){
  const G=p**a*q**b,m=p*q;
  const orders=seq(G).map(i=>i+1).filter(d=>d!==1&&G%d===0);
  const pos=new Map(orders.map((d,i)=>[d,i]));
  const forms=orders.map(d=>{
   const c=phi(d),r=Array(c.length-1).fill(0);r[0]=1;const R=[];
   for(let i=0;i<d;i++){
    R.push(r.slice());const top=r.at(-1);
    for(let j=r.length-1;j>=0;j--)r[j]=(j===0?0:r[j-1])-top*c[j];
   }
   if(r[0]!==1||r.slice(1).some(x=>x))throw Error("Monomial cycle");
   return seq(G).map(i=>R[i%d]);
  });
  const counts=orders.map((d,k)=>forms[k][0].slice()),byMask=new Map();
  function add(x,mul){
   for(let i=0;i<orders.length;i++){
    const f=forms[i][x],c=counts[i];
    for(let j=0;j<c.length;j++)c[j]+=f[j]*mul;
   }
  }
  let checked=0,spectral=0,certified=0,mirrors=0,totalPairs=0;
  function accepts(mask,d){return d===1?false:(mask>>pos.get(d))&1;}
  const active=[];
  function lookup(mask){
   if(byMask.has(mask))return byMask.get(mask);
   function ok(i,j){return accepts(mask,G/gcd(i-j,G));}
   const rootsP=seq(a).map(i=>i+1).filter(r=>accepts(mask,p**r));
   const rootsQ=seq(b).map(i=>i+1).filter(s=>accepts(mask,q**s));
   let c=0,cert=false;
   if(rootsP.length===1&&rootsQ.length===1)cert=!!accepts(mask,p**rootsP[0]*q**rootsQ[0]);
   function visit(C){
    if(active.length===m-1){
     c++;
     if(!cert)throw Error("Missing tile marks");
     const r=rootsP[0],s=rootsQ[0];
     for(const d of [p**(a-r+1),q**(b-s+1)]){
      const i=pos.get(d),sum=forms[i][0].slice();
      for(const x of active)for(let j=0;j<sum.length;j++)sum[j]+=forms[i][x][j];
      if(sum.some(x=>x))throw Error("Reverse marks");
     }
     mirrors++;
     return;
    }
    for(let i=0;i<C.length&&C.length-i>=m-1-active.length;i++){
     active.push(C[i]);visit(C.slice(i+1).filter(j=>ok(C[i],j)));active.pop();
    }
   }
   visit(seq(G).filter(i=>i&&ok(i,0)));
   if((c>0)!==cert)throw Error("Certificate classification");
   if(cert){
    const r=rootsP[0],s=rootsQ[0];
    const Lp=Math.max(0,...seq(a).map(i=>i+1).filter(l=>l!==r&&accepts(mask,p**l*q**s)));
    const Lq=Math.max(0,...seq(b).map(i=>i+1).filter(t=>t!==s&&accepts(mask,p**r*q**t)));
    const F=p**((p-1)*(r-1)+p*(q-1)*Math.min(Lp,r-1))*q**((q-1)*(s-1)+q*(p-1)*Math.min(Lq,s-1));
    const pred=F*(p**((q-1)*Math.max(Lp-r,0))+q**((p-1)*Math.max(Lq-s,0))-1);
    if(c!==pred)throw Error("Exact normalized count");
   }
   const r={c,cert};byMask.set(mask,r);return r;
  }
  function enumerate(i,remaining){
   if(!remaining){
    checked++;let mask=0;
    for(let j=0;j<counts.length;j++)if(counts[j].every(x=>x===0))mask|=1<<j;
    const r=lookup(mask);
    if(r.c){spectral++;totalPairs+=r.c;}
    if(r.cert)certified++;
    return;
   }
   for(let x=i;x<=G-remaining;x++){
    add(x,1);enumerate(x+1,remaining-1);add(x,-1);
   }
  }
  enumerate(1,m-1);
  return {G,cardinality:m,checked,rootMasks:byMask.size,spectral,certified,totalPairs,uniqueMaskSpectraMirrors:mirrors};
 }
 const result=[check(2,3,3,1),check(2,3,2,2)];
 log(result);return result;
}
exhaustiveSemiprimeSmall();
