// Replays finite examples and the exact spectrum counting formula. General claims are proved in the accompanying note.
function verifyPrimePowerResults(log=console.log) {
 const seq=n=>Array.from({length:n},(_,i)=>i),mod=(x,n)=>(x%n+n)%n;
 function gcd(a,b){a=Math.abs(a);b=Math.abs(b);while(b){const r=a%b;a=b;b=r;}return a;}
 function phiDivision(f,g) {
  const work=f.slice(),q=Array(f.length-g.length+1).fill(0n),d=g.length-1;
  for(let i=work.length-1;i>=d;i--){
   const c=q[i-d]=work[i];
   if(c)for(let j=0;j<=d;j++)work[i-d+j]-=c*g[j];
  }
  return [q,work.slice(0,d).every(x=>x===0n)];
 }
 const cache=new Map();
 function phi(n) {
  if(cache.has(n))return cache.get(n);
  let f=Array(n+1).fill(0n);f[0]=-1n;f[n]=1n;
  for(let d=1;d<n;d++)if(n%d===0){
   const [g,z]=phiDivision(f,phi(d));
   if(!z)throw Error("Cyclotomic division");f=g;
  }
  cache.set(n,f);return f;
 }
 function zero(E,n) {
  const f=Array(n).fill(0n);for(const e of E)f[mod(e,n)]++;
  return phiDivision(f,phi(n))[1];
 }
 function v(u,p){let d=0;while(u%p===0){u/=p;d++;}return d;}
 function caseCounts(A,p,extra){
  const dif=A.flatMap((a,i)=>A.slice(0,i).map(b=>Math.abs(a-b)));
  const levels=[...new Set(dif.map(d=>1+v(d,p)))].sort((a,b)=>a-b),a=levels.at(-1),P=p**a;
  const delta=levels.map(r=>dif.filter(t=>t%p**(r-1)===0).map(t=>t/p**(r-1)).reduce(gcd,0));
  const D=delta.at(-1),M=P*D,N=M*extra;
  const Z=seq(N).map(d=>zero(A.map(x=>x*d),N));
  const exact=new Set(),path=[0];
  let nodes=0;
  function visit(C){
   nodes++;
   if(nodes>10000000)throw Error("Exact search guard");
   if(path.length===A.length){exact.add(path.join(","));return;}
   for(let i=0;i<C.length&&C.length-i>=A.length-path.length;i++){
    const c=C[i];path.push(c);visit(C.slice(i+1).filter(u=>Z[mod(u-c,N)]));path.pop();
   }
  }
  visit(seq(N).filter(i=>Z[i]));
  const active=new Set(levels.map(r=>a-r));
  const projected=[];
  function tree(prefixes,depth){
   if(depth===a){projected.push(prefixes.slice().sort((u,v)=>u-v));return;}
   const power=p**depth;
   if(active.has(depth)){tree(prefixes.flatMap(u=>seq(p).map(t=>u+t*power)),depth+1);return;}
   const newPrefix=[];
   function assign(i){
    if(i===prefixes.length){tree(newPrefix,depth+1);return;}
    for(let t=0;t<(prefixes[i]===0?1:p);t++){
     newPrefix.push(prefixes[i]+t*power);assign(i+1);newPrefix.pop();
    }
   }
   assign(0);
  }
  tree([0],0);
  const predicted=new Set(),tags=[0];
  for(const base of projected) {
   function extend(i){
    if(i===base.length){
     predicted.add(base.map((x,j)=>mod((x*D+tags[j]*P)*extra,N)).sort((a,b)=>a-b).join(","));
     return;
    }
    for(let t=0;t<D;t++){
     if(tags.every((s,j)=>{
      const r=a-v(Math.abs(base[i]-base[j]),p);
      return delta[levels.indexOf(r)]%(D/gcd(D,t-s))===0;
     })){
      tags.push(t);extend(i+1);tags.pop();
     }
    }
   }
   extend(1);
  }
  if(exact.size!==predicted.size||[...predicted].some(k=>!exact.has(k)))throw Error("Spectrum family mismatch");
  let power=0,countActive=0;
  for(let d=0;d<a;d++)if(active.has(d))countActive++;else power+=p**countActive-1;
  const formula=p**power*delta.reduce((prod,d,j)=>prod*d**((p-1)*p**(levels.length-j-1)),1);
  if(formula!==exact.size)throw Error("Counting formula");
  const denominators=[...new Set([...exact].map(key=>N/key.split(",").map(Number).reduce(gcd,N)))].sort((x,y)=>x-y);
  if(denominators.join(",")!==seq(D).map(t=>t+1).filter(t=>D%t===0).map(t=>P*t).join(","))throw Error("Least denominators");
  const zerosIn18=seq(18).filter(i=>zero(A.map(x=>x*i),18));
  return {A,p,N,levels,delta,nodes,projected:projected.length,spectra:exact.size,denominators,zerosIn18};
 }
 const crt=(x,c)=>seq(24).find(y=>y%8===x&&y%3===c);
 const sameOne=seq(8).map(x=>crt(x,[0,0,1,1][x%4]));
 const sameTwo=seq(8).map(x=>crt(x,[0,1,2,0][x%4]));
 const A=[
  [[0,6],2,1],[[0,3,6,9],2,1],[[0,3,10,13],2,7],
  [sameOne,2,5],[sameTwo,2,5],
  [seq(9).map(i=>i*2),3,1],[seq(9).map(i=>(i%3)+18*((i/3)|0)),3,5],
  [seq(9).map(i=>(i%3)+12*((i/3)|0)),3,5],
  [seq(5).map(i=>6*i),5,1],[[0,2,8,10],2,3],
 ];
 const examples=A.map(([P,p,m])=>caseCounts(P,p,m));
 function checkTwist(p,a) {
  const m=1<<(a+1),Q=p**a,rows=seq(m).map(i=>seq(a).reduce((x,k)=>x+((i>>k)&1)*p**k,0));
  const cols=seq(m).map(i=>seq(a).reduce((x,k)=>x+((i>>(a-k))&1)*p**k,0));
  const N=2*Q,byStage=[],blocks=[];
  for(let k=0;k<=a;k++){
   const h=p**(a-k),C=cols.map(x=>x%h);
   let checked=0;
   for(let i=0;i<m;i++)for(let j=0;j<i;j++){
    const phases=seq(m).map(c=>{
     let t=(i^j)&c,parity=0;while(t){parity^=1;t&=t-1;}
     return parity*Q+2*(rows[i]-rows[j])*C[c];
    });
    if(!zero(phases,N))throw Error("Prime power twist stage");checked++;
   }
   byStage.push(checked);
  }
  for(let t=1;t<=a;t++){
   const divisor=p**t,rowDivisor=p**(a-t+1),prefixes=[...new Set(cols.map(c=>c%divisor))];
   let checked=0;
   for(let i=0;i<m;i++)for(let j=0;j<i;j++)if((rows[i]-rows[j])%rowDivisor!==0){
    for(const u of prefixes){
     let s=0;for(let c=0;c<m;c++)if(cols[c]%divisor===u){
      let x=(i^j)&c,parity=0;while(x){parity^=1;x&=x-1;}s+=parity?-1:1;
     }
     if(s)throw Error("Complementary prefix-block check");checked++;
    }
   }
   blocks.push(checked);
  }
  return{p,a,dimension:m,byStage,blocks};
 }
 const twists=[[3,1],[3,2],[3,3],[5,2],[5,3],[7,2],[2,2]].map(([p,a])=>checkTwist(p,a));
 const result={examples,twists};log(result);return result;
}
verifyPrimePowerResults();
