// Exact finite checks for Zhang's two seeds and the smaller spectral non-tile lift.
function checkPlanarFuglede(log=console.log) {
 const seq=n=>Array.from({length:n},(_,i)=>i);
 const mod=(x,n)=>(x%n+n)%n;
 function equal(a,b) {if(a!==b)throw Error("Exact check failed: "+a+" versus "+b);}
 const H=seq(144).map(i=>[i%12,Math.floor(i/12)]);
 const G=seq(720).map(i=>[i%60,Math.floor(i/60)]);
 const E=[0,5,0,7,9,2,10,4,7,8,10,10].map((x,y)=>[x,y]);
 const B=[seq(12).map(x=>[x,0]),seq(12).map(x=>[x,mod(-x,12)]),seq(3).flatMap(a=>seq(2).flatMap(b=>seq(2).map(c=>[mod(3*c+4*a,12),mod(4*a+6*b,12)])))];
 const V=[seq(12).map(y=>[0,y]),seq(12).map(x=>[x,x]),seq(2).flatMap(d=>seq(6).map(k=>[mod(4*k+d,12),mod(2*k+2*d,12)]))];
 const crt=(x,y,t)=>[mod(25*x+36*t,60),y];
 const Aminus=[0,2,2,1,1].flatMap((i,t)=>B[i].map(([x,y])=>crt(x,y,t)));
 const Aplus=[0,0,0,2,1].flatMap((i,t)=>V[i].map(([x,y])=>crt(x,y,t)));
 const Bminus=E.map(([x,y])=>crt(x,y,0));
 const Lplus=E.flatMap(([u,v])=>seq(5).map(r=>[mod(5*u+12*r,60),v]));
 const codeG=(x,y)=>mod(x,60)+60*mod(y,12);
 const codeH=(x,y)=>mod(x,12)+12*mod(y,12);
 for(const F of [Aminus,Aplus,Lplus])equal(new Set(F.map(([x,y])=>codeG(x,y))).size,60);
 equal(new Set(Aminus.flatMap(([x,y])=>Bminus.map(([u,v])=>codeG(x+u,y+v)))).size,720);
 function quotient(f,d) {
  const work=f.slice(),n=d.length-1,out=new Array(f.length-n).fill(0);
  for(let j=work.length-1;j>=n;j--){
   const x=work[j];out[j-n]=x;
   if(x)for(let k=0;k<=n;k++)work[j-n+k]-=x*d[k];
  }
  return {quotient:out,remainder:work.slice(0,n)};
 }
 function divisor(n) {
  const phis=[null];
  for(let k=1;k<=n;k++){
   let C=new Array(k+1).fill(0);C[0]=-1;C[k]=1;
   for(let j=1;j<k;j++)if(k%j===0){const z=quotient(C,phis[j]);equal(z.remainder.every(a=>a===0),true);C=z.quotient;}
   phis.push(C);
  }
  return phis[n];
 }
 const phi12=divisor(12),phi60=divisor(60);
 function zero(F,u,v,order) {
  const n=order===12?12:60,weight=order/12;
  let f=new Array(n).fill(0);
  for(const [x,y] of F)f[mod(u*x+weight*v*y,n)]++;
  return quotient(f,n===12?phi12:phi60).remainder.every(a=>a===0);
 }
 const positiveZero=new Uint8Array(720);
 for(const [u,v] of G)positiveZero[codeG(u,v)]=+zero(Aplus,u,v,60);
 let spectralPairs=0;
 for(let i=0;i<Lplus.length;i++)for(let j=0;j<i;j++){
  const [u,v]=Lplus[i],[x,y]=Lplus[j];
  equal(positiveZero[codeG(u-x,v-y)],1);spectralPairs++;
 }
 let FourierChecks=0;
 const commonZero=new Uint8Array(144);
 for(const [u,v] of H){
  const z=B.every(F=>zero(F,u,v,12));
  commonZero[codeH(u,v)]=+z;
  equal(z,zero(Aminus,5*u,v,60));FourierChecks++;
 }
 const AplusD=new Set(Aplus.flatMap(([x,y])=>Aplus.map(([u,v])=>codeG(x-u,y-v))));
 for(const [x,y] of G)if(x%5!==0)equal(AplusD.has(codeG(x,y)),true);
 const forbiddenH=new Set(V.flatMap(F=>F.flatMap(([x,y])=>F.map(([u,v])=>codeH(x-u,y-v)))));
 for(const [x,y] of H)equal(AplusD.has(codeG(...crt(x,y,0))),forbiddenH.has(codeH(x,y)));
 function maxClique(predicate) {
  const adj=H.map(([x,y])=>H.map(([u,v])=>predicate(mod(x-u,12),mod(y-v,12))));
  let best=[],nodes=0;
  // All cliques are translation invariant. One can be shifted to contain index 0.
  function color(P) {
   let todo=P.slice(),out=[],bound=[],color=0;
   while(todo.length){
    color++;
    const rest=[],cl=[];
    for(const v of todo){if(cl.some(u=>adj[u][v]))rest.push(v);else{cl.push(v);out.push(v);bound.push(color);}}
    todo=rest;
   }
   return [out,bound];
  }
  function visit(C,P) {
   nodes++;
   if(C.length>best.length)best=C.slice();
   const [out,bound]=color(P);
   for(let i=out.length-1;i>=0;i--){
    if(C.length+bound[i]<=best.length)return;
    const v=out[i];
    visit(C.concat(v),out.slice(0,i).filter(u=>adj[v][u]));
   }
  }
  visit([0],seq(144).filter(v=>adj[0][v]));
  return {maximum:best.length,witness:best.map(i=>H[i]),nodes};
 }
 const packing=maxClique((x,y)=>!forbiddenH.has(codeH(x,y)));
 equal(packing.maximum,10);
 const commonSpec=maxClique((x,y)=>!!commonZero[codeH(x,y)]);
 if(commonSpec.maximum>=12)throw Error("Common spectrum found");
 // Change representatives by a group automorphism and translate on the torus.
 const move=(x,y)=>[mod(2*x+15*y-23,60),mod(x+8*y,12)];
 equal(new Set(G.map(([x,y])=>codeG(...move(x,y)))).size,720);
 const S=Aplus.map(([x,y])=>move(x,y));
 const deltas=Array.from({length:720},()=>[]);
 for(const [x,y] of S)for(const [u,v] of S){
  const e=x-u,f=y-v;deltas[codeG(e,f)].push([e,f]);
 }
 // Every modular conflict among translation origins in this window is a genuine difference.
 const kx=8,ky=9,Lx=437,Ly=99;let residuesChecked=0,conflictsChecked=0;
 for(let dx=1-Lx;dx<Lx;dx++)for(let dy=1-Ly;dy<Ly;dy++){
  residuesChecked++;
  const array=deltas[codeG(dx,dy)];
  if(!array.length)continue;
  conflictsChecked++;
  if(!array.some(([e,f])=>Math.abs(dx-e)<=60*(kx-1)&&Math.abs(dy-f)<=12*(ky-1)))throw Error("Missing actual difference");
 }
 const squares=60*kx*ky,box=Lx*Ly;
 if(box<=squares*packing.maximum)throw Error("Window density does not contradict tiling");
 const result={positiveFiniteSpectralPairs:spectralPairs,negativeFourierChecks:FourierChecks,packing,commonSpec,residuesChecked,conflictsChecked,squares,kx,ky,Lx,Ly,windowArea:box,threshold:squares*packing.maximum};
 log(result);return result;
}
checkPlanarFuglede();

// Independent search by x-coordinate for an eleven-center packing.
function verifyElevenPacking(log=console.log) {
 const mod=n=>(n%12+12)%12;
 const seq=n=>Array.from({length:n},(_,i)=>i);
 const forbidden=Array.from({length:12},()=>new Set());
 const V=[seq(12).map(y=>[0,y]),seq(12).map(x=>[x,x]),seq(2).flatMap(d=>seq(6).map(k=>[mod(4*k+d),mod(2*k+2*d)]))];
 for(const F of V)for(const [x,y] of F)for(const [u,v] of F)forbidden[mod(x-u)].add(mod(y-v));
 if(forbidden[0].size!==12)throw Error("One point per x");
 const domains=seq(12).map(x=>seq(12).filter(y=>!forbidden[x].has(y)));
 const nodes=[];
 for(let missing=1;missing<12;missing++){
  let count=0;
  function visit(rows) {
   count++;
   if(!rows.length)return true;
   rows.sort((a,b)=>a[1].length-b[1].length);
   const [x,choices]=rows[0];
   for(const y of choices){
    const rest=rows.slice(1).map(([u,ys])=>[u,ys.filter(v=>!forbidden[mod(u-x)].has(mod(v-y)))]);
    if(rest.every(([,ys])=>ys.length)&&visit(rest))return true;
   }
   return false;
  }
  if(visit(seq(12).filter(x=>x!==0&&x!==missing).map(x=>[x,domains[x]])))throw Error("Eleven-center packing found");
  nodes.push(count);
 }
 log(nodes);return nodes;
}
verifyElevenPacking();

// Rational cosine bounds for the quantitative nonspectral tile lift.
function verifyNegativeLift(log=console.log) {
 const abs=n=>n<0n?-n:n;
 function gcd(a,b) { while(b!==0n){const r=a%b;a=b;b=r;}return abs(a); }
 function q(n,d=1n) {const g=gcd(n,d);return [n/g,d/g];}
 const add=([a,b],[c,d])=>q(a*d+b*c,b*d);
 const neg=([a,b])=>[-a,b];
 const mul=([a,b],[c,d])=>q(a*c,b*d);
 const compare=([a,b],[c,d])=>a*d-b*c;
 const lo=q(3141592653n,10n**9n),hi=q(3141592654n,10n**9n);
 function arctan(d,terms) {
  let s=q(0n);
  for(let i=0n;i<terms;i++)s=add(s,q((-1n)**i,(2n*i+1n)*d**(2n*i+1n)));
  return s;
 }
 const piMin=add(mul(q(16n),arctan(5n,12n)),neg(mul(q(4n),arctan(239n,5n))));
 const piMax=add(mul(q(16n),arctan(5n,13n)),neg(mul(q(4n),arctan(239n,6n))));
 if(compare(lo,piMin)>=0n || compare(piMax,hi)>=0n || compare(hi,q(22n,7n))>=0n)throw Error("Pi interval");
 let f24=1n;for(let i=2n;i<=24n;i++)f24*=i;
 if(22n**24n*10n**9n>=7n**24n*f24)throw Error("Taylor remainder");
 const scale=10n**8n,rounds=[];
 for(let j=0;j<60;j++){
  const angle=mul(lo,q(BigInt(Math.min(j,60-j)),30n));
  const minusSq=neg(mul(angle,angle));
  let value=q(1n),term=q(1n);
  for(let i=1n;i<12n;i++){
   term=mul(mul(term,minusSq),q(1n,(2n*i-1n)*2n*i));value=add(value,term);
  }
  const [a,b]=mul(value,q(scale));
  const floor=a/b-(a<0n && a%b!==0n?1n:0n);
  rounds.push(floor-1n);
 }
 const mod=(x,n)=>(x%n+n)%n;
 const seq=n=>Array.from({length:n},(_,i)=>i);
 const B=[seq(12).map(x=>[x,0]),seq(12).map(x=>[x,mod(-x,12)]),seq(3).flatMap(a=>seq(2).flatMap(b=>seq(2).map(c=>[mod(3*c+4*a,12),mod(4*a+6*b,12)])))];
 const A=[0,2,2,1,1].flatMap((i,t)=>B[i].map(([x,y])=>[mod(25*x+36*t,60),y]));
 const phi60=[1,0,1,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,1];
 let nonzero=0,lowest=null;
 for(let u=0;u<60;u++)for(let v=0;v<12;v++){
  const f=Array(60).fill(0);
  for(const [x,y] of A)f[mod(u*x+5*v*y,60)]++;
  const rem=f.slice();
  for(let i=59;i>=16;i--){
   const a=rem[i];
   for(let j=0;j<=16;j++)rem[i-16+j]-=a*phi60[j];
  }
  if(rem.every(a=>a===0))continue;
  nonzero++;
  let lower=0n;
  for(let i=0;i<60;i++)if(f[i])for(let j=0;j<60;j++)if(f[j])lower+=BigInt(f[i]*f[j])*rounds[mod(i-j,60)];
  if(625n*lower<=17161n*scale)throw Error("Nonzero Fourier lower bound");
  if(lowest===null || lower<lowest)lowest=lower;
 }
 const center=(x,n)=>mod(x+n/2,n)-n/2;
 const transform=(x,y)=>[center(7*x+20*y-9,60),center(x+y,12)];
 const all=seq(720).map(i=>transform(i%60,Math.floor(i/60)).join(","));
 if(new Set(all).size!==720)throw Error("Torus map");
 const S=A.map(([x,y])=>transform(x,y));
 const sx=S.reduce((sum,[x])=>sum+Math.abs(x),0);
 const sy=S.reduce((sum,[,y])=>sum+Math.abs(y),0);
 if(sx!==720||sy!==108)throw Error("Centered moment");
 const kx=30,ky=21,nx=60,ny=12;
 const variation=mul(q(44n,7n),add(q(BigInt(sx),BigInt(nx*kx)),q(BigInt(sy),BigInt(ny*ky))));
 if(compare(variation,q(131n,25n))>=0n)throw Error("Fourier variation");
 const period=new Set();
 for(const [x,y] of S)for(let i=0;i<kx;i++)for(let j=0;j<ky;j++)period.add(mod(x+nx*i,nx*kx)+nx*kx*mod(y+ny*j,ny*ky));
 if(period.size!==60*kx*ky)throw Error("Expanded residues");
 const result={nonzero,lowestSquaredBoundNumerator:String(lowest),lowestSquaredBoundDenominator:String(scale),fourierAbsThreshold:"131/25",sx,sy,variationUpperBound:variation.map(String).join("/"),kx,ky,squares:period.size};
 log(result);return result;
}
verifyNegativeLift();
