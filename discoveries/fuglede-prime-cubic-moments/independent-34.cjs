'use strict';
const p=17, mod=x=>((x%p)+p)%p, pow=(a,d)=>{let r=1;while(d-- >0)r=r*a%p;return r};
const assert=(b,s)=>{if(!b)throw Error(s)};
const g=3, inv=Array.from({length:p},(_,i)=>i?pow(i,p-2):0), xp=Array.from({length:17},(_,i)=>Array.from({length:8},(_,d)=>pow(i,d)));
assert(new Set(Array.from({length:16},(_,i)=>pow(g,i))).size===16,'primitive root');
function leastAndCatalogue(d){
  // First center by t -> t+constant; then use the d-th-power cosets of leading coefficient.
  const reps=Array.from({length:16},(_,i)=>pow(g,i)).filter((a,i)=>i<gcd(d,16));
  const cls=new Map; for(const a of reps)for(let u=1;u<p;u++){const l=a*pow(u,d)%p;if(!cls.has(l))cls.set(l,[a,u]);}
  assert(cls.size===16,'leading coefficient cover');
  let minimum=Infinity, normalized=0;const small=[];
  if(d===1)return {minimum:136,normalized:1,small:[],cls};
  const mx=p**(d-2), cf=new Int8Array(d+1), hist=new Int8Array(p);
  for(const a of reps){
    cf[d]=a;
    for(let code=0;code<mx;code++){
      let r=code;for(let j=1;j<=d-2;j++){cf[j]=r%p;r=Math.floor(r/p)}
      hist.fill(0); let sum=0;
      for(let k=0;k<p;k++){let v=a*xp[k][d];for(let j=1;j<=d-2;j++)v+=cf[j]*xp[k][j];v%=p;hist[v]++;sum+=v;}
      normalized++;
      for(let c=0;c<p;c++){
        if(sum<minimum)minimum=sum;
        if(sum<=34){cf[0]=c;small.push({a,coeff:[...cf],sum});}
        sum+=p*(1-hist[p-1-c]);
      }
    }
  }
  return {minimum,normalized,small,cls};
}
function gcd(x,y){return y?gcd(y,x%y):x;}
const profile=[], cat=Array.from({length:8},()=>[]);
for(let d=1;d<=7;d++){
 const {minimum,normalized,small,cls}=leastAndCatalogue(d);
 for(const n of small)for(let u=1;u<p;u++){
   const coeff=n.coeff.map((v,j)=>v*pow(u,j)%p), a=coeff[d];
   if(!cat[d].some(x=>x.every((w,j)=>w===coeff[j])))cat[d].push(coeff);
 }
 profile.push({d,minimum,normalized,lowNormalized:small.length,lowAllLead:cat[d].length,lowNormalizedBySum:small.reduce((r,x)=>(r[x.sum]=(r[x.sum]||0)+1,r),{})});
}
assert(JSON.stringify(profile.slice(0,6).map(x=>x.minimum))===JSON.stringify([136,102,85,51,51,34]),'minimums');
const sext=new Map(cat[6].map(x=>[x[6],x])); assert(sext.size===16,'one per leading coeff');
assert(cat[6].every(x=>x[1]===0 && x[3]===0 && x[5]===0),'evenness');
const labelStats={nMasks:0,nTuples:0,nConstants:0,nD4:0,firstRejectExponent:[0,0,0,0,0]};
const slopes=[];
function test(){
 labelStats.nMasks++;
 const w=slopes.map(s=>inv[slopes.reduce((v,t)=>s===t?v:v*mod(s-t)%p,1)]);
 for(let lambda=1;lambda<p;lambda++){
  labelStats.nTuples++;
  const cc=w.map(v=>sext.get(v*lambda%p));
  if(cc.reduce((z,x)=>z+x[0],0)%p!==14)continue;
  labelStats.nConstants++;
  let con=true;
  for(let h=0;h<=4;h++){
   if(cc.reduce((z,x,i)=>z+x[4]*xp[slopes[i]][h],0)%p){
    labelStats.firstRejectExponent[h]++;con=false;break;
   }
  }
  if(con)labelStats.nD4++;
 }
}
function dfs(start){if(slopes.length===8)return test(); for(let x=start;x<=p-(8-slopes.length);x++){slopes.push(x);dfs(x+1);slopes.pop()}}
dfs(0);
assert(labelStats.nMasks===24310&&labelStats.nConstants===1360&&labelStats.nD4===0,'label gate');
const primeArith={}, rest=[];
for(let m=2;m<=16;m++){
 if((17*(17-m))**2<16**2*17){primeArith[m]='Malikiosis';continue;}
 let n=m,why=null;
 if(m%2)for(let q=2;q<=n;q++){
  let e=0;while(n%q===0){e++;n/=q;}
  if(e%2&&Array.from({length:16},(_,j)=>pow(q,j)).includes(16)){why='self-conjugate '+q;break;}
 }
 if(why)primeArith[m]=why;else rest.push(m);
}
assert(JSON.stringify(rest)===JSON.stringify([2,4,6,8,9,10,12]),'ledger');

const septicMax=c=>Math.max(...xp.map(x=>c.reduce((s,a,i)=>s+a*x[i],0)%p));
const sepLead=new Map(Array.from({length:16},(_,i)=>[i+1,cat[7].filter(c=>c[7]===i+1&&septicMax(c)<=9)]));
assert([...sepLead.values()].every(x=>x.length===4),'septic cap catalogue');
const eta=Array.from({length:p},(_,s)=>mod(s*(s-1)));
const evalg=(c,x)=>c.reduce((a,b,i)=>a+b*xp[x][i],0)%p;
const stats9={
 masks:0,
 septic:{quadraticAllowances:0,compatibleQuadratics:0,compatibleAssignments:0,constants:0,degree4:0,firstReject:[0,0,0,0,0]},
 sextic:{choices:0,constants:0,degree4:0,firstReject:[0,0,0,0,0]}
};
const S=[0,1];
function degree4(stats,v,g){
 for(let h=0;h<=4;h++){
   if(g.reduce((z,x,i)=>z+x*xp[v[i]][h],0)%p){stats.firstReject[h]++;return}
 }
 stats.degree4++;
}
function test9(){
 stats9.masks++;
 const w=S.map(s=>inv[S.reduce((v,t)=>s===t?v:v*mod(s-t)%p,1)]);
 const sparse=w.map(v=>{
   const a=new Map;
   for(const c of sepLead.get(v)){
     const key=c[5]*inv[v]%p;
     const xs=a.get(key)||[];xs.push(c);a.set(key,xs);
   }
   return a;
 });
 const aStat=stats9.septic, bStat=stats9.sextic;
 const picks=[];
 function trySeptic(i,c0){
   if(i===9){
     aStat.compatibleAssignments++;
     if(c0%p===16){aStat.constants++;degree4(aStat,S,picks);}
     return;
   }
   for(const c of possible[i]){
      picks[i]=c[4];trySeptic(i+1,c0+c[0]);
   }
 }
 const possible=[];
 // Test every possible leading quadratic coefficient and actual first two values.
 for(const a0 of sparse[0].keys())for(const a1 of sparse[1].keys())for(let A=0;A<p;A++){
   aStat.quadraticAllowances++;
   let i=0;
   while(i<9){
      const s=S[i];const key=mod(a0+(a1-a0)*s+A*eta[s]);
      const opts=sparse[i].get(key);if(!opts)break; possible[i++]=opts;
   }
   if(i===9){aStat.compatibleQuadratics++;trySeptic(0,0);}
 }
 // For sextics, choose projective root of c6/w as the uniform vertical label.
 // The remaining homogeneous form makes c6/w constant. A common scale makes it 1 or 3.
 for(const l of [1,3])for(let h=0;h<p;h++){
   bStat.choices++;
   let csum=0;const q=[];
   for(let i=0;i<9;i++){
      const a=l*w[i]%p, c=sext.get(a), x=h*eta[S[i]]%p;
      csum+=evalg(c,x);q[i]=(c[4]+15*a*x*x)%p;
   }
   if(csum%p===16){bStat.constants++;degree4(bStat,S,q);}
 }
}
function choose9(a){if(S.length===9)return test9();for(let s=a;s<=p-(9-S.length);s++){S.push(s);choose9(s+1);S.pop()}}
choose9(2);
assert(stats9.masks===6435&&stats9.septic.compatibleAssignments===18909&&stats9.septic.constants===81&&stats9.septic.degree4===0,'septic K9 gate');
assert(stats9.sextic.choices===218790&&stats9.sextic.constants===12538&&stats9.sextic.degree4===0,'root-at-infinity sextic K9 gate');
console.log(JSON.stringify({profile,sexticCatalogue:[...sext].sort((a,b)=>a[0]-b[0]),leadOneSeptics:cat[7].filter(c=>c[7]===1).map(c=>({c,max:septicMax(c)})),labelStats,stats9,primeArith,rest},null,2));

