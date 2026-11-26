'use strict';
const OCTIC=(()=>{
'use strict';
const p=17, mod=x=>((x%17)+17)%17, pw=(a,b)=>{let r=1;while(b--)r=r*a%17;return r};
const IV=Array.from({length:p},(_,i)=>pw(i,15)), CAP=8;
const xs=Array.from({length:8},(_,i)=>i+1), sq=xs.map(x=>x*x%17), v4=sq.map(x=>pw(x,4));
const lifts=Array.from({length:p},()=>[]);
for(let l=0;l<=CAP;l++)for(let r=0;r<=CAP;r++){lifts[mod(l-r)].push({half:(l+r)*9%17,total:l+r});}
const best=[...lifts].map(L=>Math.min(...L.map(x=>x.total)));
const stats={oddFeasible:0,leaves:0,exact:0}, masks=new Map, catalogue=Array.from({length:17},()=>[]);
function lagrange(S){
 const id=S.join(',');if(masks.has(id))return masks.get(id);
 const nodes=[0,...S.map(i=>sq[i])];
 const out=sq.map(t=>nodes.map((s,j)=>
   nodes.reduce((a,u,k)=>k===j?a:a*mod(t-u)*IV[mod(s-u)]%17,1)));
 masks.set(id,out);return out;
}
for(let a=0;a<17;a++)for(let b=0;b<17;b++)for(let c=0;c<17;c++){
 const od=xs.map((x,i)=>x*(a+b*sq[i]+c*sq[i]*sq[i])%17), need=od.map(x=>2*x%17);
 const mins=need.map(x=>best[x]);
 if(mins.reduce((x,y)=>x+y,0)>34)continue;
 stats.oddFeasible++;
 const I=Array.from({length:8},(_,i)=>i).sort((i,j)=>lifts[need[i]].length-lifts[need[j]].length||j-i).slice(0,3).sort((a,b)=>a-b);
 const J=Array.from({length:8},(_,i)=>i).filter(i=>!I.includes(i)).sort((i,j)=>lifts[need[i]].length-lifts[need[j]].length);
 const bas=34-J.reduce((u,j)=>u+mins[j],0), lam=lagrange(I);
 const r=[0,0,0,0], ev=new Int8Array(8);
 for(let lead=1;lead<17;lead++){
  const d=v4.map(x=>x*lead%17);
  for(let z=0;z<=8&&z<=bas;z++){
   r[0]=z;
   for(const v1 of lifts[need[I[0]]]){
    r[1]=mod(v1.half-d[I[0]]);
    for(const v2 of lifts[need[I[1]]]){
     r[2]=mod(v2.half-d[I[1]]);
     for(const v3 of lifts[need[I[2]]]){
      const initial=z+v1.total+v2.total+v3.total;
      if(initial>bas)continue;
      r[3]=mod(v3.half-d[I[2]]);stats.leaves++;
      let cost=initial, i=0;
      for(;i<J.length;i++){
       const j=J[i], W=lam[j];
       ev[j]=(d[j]+W[0]*r[0]+W[1]*r[1]+W[2]*r[2]+W[3]*r[3])%17;
       const pp=(ev[j]+od[j])%17, nn=mod(ev[j]-od[j]);
       if(pp>CAP||nn>CAP)break;
       cost+=pp+nn;if(cost>34)break;
      }
      if(i<J.length||cost!==34)continue;
      ev[I[0]]=v1.half;ev[I[1]]=v2.half;ev[I[2]]=v3.half;
      const vals=Array(17).fill(z);
      xs.forEach((x,j)=>{vals[x]=(ev[j]+od[j])%17;vals[p-x]=mod(ev[j]-od[j]);});
      const cf=vals.map((v,j)=>j===0?z:mod(-vals.reduce((s,u,x)=>s+u*pw(x,16-j),0)));
      if(cf[8]!==lead||cf[7]||cf.slice(9).some(Boolean))throw Error('degree');
      catalogue[lead].push(cf.slice(0,9));stats.exact++;
     }
    }
   }
  }
 }
}
return {catalogue,stats,perLeading:catalogue.map(a=>a.length),maxPerTopC6:catalogue.map(L=>Math.max(0,...Array.from({length:17},(_,v)=>L.filter(c=>c[6]===v).length)))};

})();
const ORBITS=(()=>{
'use strict';
const p=17, mod=x=>(x%17+17)%17, power=(x,d)=>{let r=1;while(d-- >0)r=r*x%17;return r}, inv=Array.from({length:17},(_,i)=>i?power(i,15):0);
const pop=x=>{let j=0;while(x){x&=x-1;j++;}return j};
const transform=(mask,perm)=>{let r=0;for(let s=0;s<=17;s++)if(mask&(1<<s))r|=1<<perm[s];return r};
const perms=[s=>s===17?17:(s+1)%17,s=>s===17?17:s*3%17,s=>s===17?0:(!s?17:inv[s])].map(f=>Array.from({length:18},(_,i)=>f(i)));
const seen=new Uint8Array(1<<18), models=[], lengths=[];
for(let m=0;m<(1<<18);m++){
 if(seen[m]||pop(m)!==10)continue;
 const fifo=[m];seen[m]=1;
 for(let i=0;i<fifo.length;i++)for(const pi of perms){
  const n=transform(fifo[i],pi);
  if(!seen[n]){seen[n]=1;fifo.push(n);}
 }
 const norm=Math.min(...fifo.filter(x=>x<(1<<17)&&(x&3)===3));
 if(!Number.isFinite(norm))throw Error('missing representative');
 models.push(norm);lengths.push(fifo.length);
}
if(models.length!==17||lengths.reduce((x,y)=>x+y,0)!==43758)throw Error('orbit coverage');
const nm=seen.reduce((a,x,i)=>a+(x&&i<(1<<17)&&(i&3)===3),0);
if(nm!==6435)throw Error('affine census');
return {models,lengths,nm};
})();
const gates=((inputA,inputB)=>{
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

const octic=inputA.catalogue, models=inputB.models;
const comb=(a,b)=>{let v=1;for(let i=1;i<=b;i++)v=v*(a-i+1)/i;return v};
const sf=(c,e,h)=>mod(c.slice(e).reduce((s,a,k)=>s+comb(k+e,e)*a*xp[h][k],0));
const sm=new Map(Array.from({length:16},(_,i)=>[i+1,cat[7].filter(c=>c[7]===i+1&&Math.max(...xp.map(x=>c.reduce((s,a,j)=>s+a*x[j],0)%p))<=8)]));
const stat={octic:{quadraticTrials:0,quadratics:0,cubicTrials:0,cubics:0,assignments:0,constants:0,passes:0,firstReject:[0,0,0,0,0]},
septic:{leadAndShift:0,cubicTrials:0,cubics:0,assignments:0,constants:0,passes:0,firstReject:[0,0,0,0,0]},sextic:{validLd:0,assignments:0,constants:0,passes:0,firstReject:[0,0,0,0,0]}};
function bucket(A,key){const b=new Map;for(const a of A){const k=key(a), z=b.get(k)||[];z.push(a);b.set(k,z)}return b}
const sht=Array.from({length:17},(_,a)=>Array.from({length:17},(_,h)=>a?{
 b0:sf(sext.get(a),0,h),b4:sf(sext.get(a),4,h)
}:null));
for(const m of models){
 const S=Array.from({length:17},(_,i)=>i).filter(s=>m&(1<<s));
 const W=S.map(s=>inv[S.reduce((v,t)=>s===t?v:v*mod(s-t)%p,1)]),Wi=W.map(x=>inv[x]);
 const binomW=S.map(s=>xp[s].slice(0,5)); 
 const fc=new Map;
 function weights(Q){
   const k=Q.join(',');if(fc.has(k))return fc.get(k);
   const ts=Q.map(i=>S[i]);const a=S.map(s=>ts.map((t,i)=>ts.reduce((v,u,j)=>i===j?v:v*mod(s-u)*inv[mod(t-u)]%17,1)));
   fc.set(k,a);return a;
 }
 function testAssignment(O,ss,c0){
  if(ss){
   let b=0;for(const c of O)b+=c.b0;
   c0=b;
  }
  if(c0%p!==1)return; // 2(10-1)=18=1 (mod 17)
  ss.constants++;
  for(let e=0;e<5;e++){
   let y=0;for(let i=0;i<10;i++)y+=O[i].b4*binomW[i][e];
   if(y%p){ss.firstReject[e]++;return;}
  }
  ss.passes++;
 }
 function verifyCubic(A,ss,forceIndex){
  const T=Array.from({length:10},(_,i)=>i).filter(i=>i!==forceIndex).sort((i,j)=>A[i].size-A[j].size||j-i).slice(0,4);
  const Z=weights(T), keys=[];
  const opts=[],chosen=[];
  function assignments(i,b){
   if(i===10){
    ss.assignments++;
    if(b%17===1)testAssignment(chosen,ss,b);
    return;
   }
   for(const c of opts[i]){chosen[i]=c;assignments(i+1,b+c.b0);}
  }
  function rec(j){
   if(j<4){
    for(const k of A[T[j]].keys()){keys[j]=k;rec(j+1);}
    return;
   }
   ss.cubicTrials++;
   for(let i=0;i<10;i++){
    const P=Z[i], q=(P[0]*keys[0]+P[1]*keys[1]+P[2]*keys[2]+P[3]*keys[3])%17;
    const a=A[i].get(q); if(!a)return;
    opts[i]=a;
   }
   ss.cubics++; assignments(0,0);
  }
  rec(0);
 }
 for(let l=1;l<=8;l++){
  const Q=W.map((w,i)=>bucket(octic[w*l%17],c=>c[6]*Wi[i]%17));
  const T=Array.from({length:10},(_,i)=>i).sort((i,j)=>Q[i].size-Q[j].size||j-i).slice(0,3),Z=weights(T),c6=[];
  function match6(j){
   if(j<3){for(const key of Q[T[j]].keys()){c6[j]=key;match6(j+1);}return;}
   stat.octic.quadraticTrials++;
   const P=[];
   for(let i=0;i<10;i++){
    const a=Z[i], v=(a[0]*c6[0]+a[1]*c6[1]+a[2]*c6[2])%17;
    const F=Q[i].get(v);if(!F)return;
    P[i]=bucket(F.map(c=>({b0:c[0],b4:c[4],g5:c[5]*Wi[i]%17})),x=>x.g5);
   }
   stat.octic.quadratics++;verifyCubic(P,stat.octic,-1);
  }
  match6(0);
 }
 // All affine top profile ratios up to common nonzero scalar; seventh powers are all scalars.
 for(const [b,a] of [[1,0],...Array.from({length:17},(_,i)=>[i,1])]){
  const ll=S.map((s,i)=>mod(W[i]*(b+a*s))), iz=ll.indexOf(0);
  const J=Array.from({length:10},(_,i)=>i).filter(i=>i!==iz).slice(0,2);
  const R=S.map(s=>mod((s-S[J[0]])*(s-S[J[1]])));
  for(let F=0;F<17;F++){
   const c6=W.map((w,i)=>w*R[i]*F%17);
   if(iz>=0&&!c6[iz])continue;
   stat.septic.leadAndShift++;
   const A=ll.map((lead,i)=>{
    if(!lead){
     const out=new Map;
     for(let h=0;h<17;h++)out.set((6*c6[i]*h*Wi[i])%17,[sht[c6[i]][h]]);
     return out;
    }
    const h=mod(c6[i]*inv[7*lead%17]);
    const v=sm.get(lead).map(P=>({b0:sf(P,0,h),b4:sf(P,4,h),key:sf(P,5,h)*Wi[i]%17}));
    return bucket(v,P=>P.key);
   });
   verifyCubic(A,stat.septic,iz);
  }
 }
 // For square-class representatives, fix the leading coefficient at the first bad slope.
 const s0=S[0],s1=S[1],R=S.map(s=>mod((s-s0)*(s-s1)));
 for(const l of [1,3])for(let a=0;a<17;a++)for(let b=0;b<17;b++){
  const ld=W.map((w,i)=>mod(w*(l*Wi[0]+a*(S[i]-s0)+b*R[i])));
  if(ld.includes(0))continue;
  stat.sextic.validLd++;
  const H=ld.map((x,i)=>W[i]*inv[6*x%17]*R[i]%17);
  for(let u=0;u<17;u++)for(let t=0;t<17;t++){
   stat.sextic.assignments++;
   const Q=ld.map((x,i)=>sht[x][mod(H[i]*(u*S[i]+t))]);
   testAssignment(Q,stat.sextic,0);
  }
 }
}
if(stat.octic.constants!==5652||stat.septic.constants!==110||stat.sextic.constants!==97361||
    stat.octic.passes||stat.septic.passes||stat.sextic.passes)throw Error('K8 gate '+JSON.stringify(stat));
return stat;

})(OCTIC,ORBITS);
if(OCTIC.stats.exact!==568||OCTIC.stats.oddFeasible!==2625)throw Error('octic census');
const {catalogue,...octics}=OCTIC;
console.log(JSON.stringify({octics,orbits:ORBITS,gates},null,2));

