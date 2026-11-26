'use strict';
const BAS=(()=>{
'use strict';
const p=17, mod=x=>(x%17+17)%17, power=(x,d)=>{let r=1;while(d-- >0)r=r*x%17;return r}, inv=Array.from({length:17},(_,i)=>i?power(i,15):0);
const pop=x=>{let j=0;while(x){x&=x-1;j++;}return j};
const transform=(mask,perm)=>{let r=0;for(let s=0;s<=17;s++)if(mask&(1<<s))r|=1<<perm[s];return r};
const perms=[s=>s===17?17:(s+1)%17,s=>s===17?17:s*3%17,s=>s===17?0:(!s?17:inv[s])].map(f=>Array.from({length:18},(_,i)=>f(i)));
const seen=new Uint8Array(1<<18), models=[], lengths=[];
for(let m=0;m<(1<<18);m++){
 if(seen[m]||pop(m)!==11)continue;
 const fifo=[m];seen[m]=1;
 for(let i=0;i<fifo.length;i++)for(const pi of perms){
  const n=transform(fifo[i],pi);
  if(!seen[n]){seen[n]=1;fifo.push(n);}
 }
 const norm=Math.min(...fifo.filter(x=>x<(1<<17)&&(x&3)===3));
 if(!Number.isFinite(norm))throw Error('missing representative');
 models.push(norm);lengths.push(fifo.length);
}
if(models.length!==10||lengths.reduce((x,y)=>x+y,0)!==31824)throw Error('orbit coverage');
const nm=seen.reduce((a,x,i)=>a+(x&&i<(1<<17)&&(i&3)===3),0);
if(nm!==5005)throw Error('affine census');
return {models,lengths,nm};
})();
const CATALOGUE=(()=>{
'use strict';
const p=17, mod=x=>((x%17)+17)%17, pw=(a,b)=>{let r=1;while(b--)r=r*a%17;return r};
const IV=Array.from({length:p},(_,i)=>pw(i,15)), CAP=6;
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
  for(let z=0;z<=CAP&&z<=bas;z++){
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
for(const F of catalogue.flat()){ const V=Array.from({length:17},(_,t)=>F.reduceRight((r,x)=>(r*t+x)%17,0)); if(V.some(x=>x>CAP)||V.reduce((s,t)=>s+t,0)!==34)throw Error('profile surface');}
if(stats.exact!==272)throw Error('octic count');
return {stats,catalogue};


})();
const stats=((bas,cat8)=>{
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


const masks=bas.models, OCT=cat8.catalogue;
const B7=Array.from({length:17},(_,a)=>cat[7].filter(c=>c[7]===a&&xp.every(x=>c.reduce((s,b,j)=>s+b*x[j],0)%17<=6)));
const stats={mask:0,top:0,validTop:0,Yinterpolants:0,Yfit:0,rootShiftFit:0,fullLower:0,lead5Trials:0,lead5Fits:0,assignments:0,constant:0,passes:0,firstFail:Array(5).fill(0),per:[]};
const choose=(n,k)=>{let z=1;for(let i=1;i<=k;i++)z=z*(n-i+1)/i;return z};
const sc=F=>Array.from({length:17},(_,h)=>Array.from({length:7},(_,j)=>F.slice(j).reduce((a,b,i)=>a+choose(i+j,j)*b*pow(h,i),0)%17));
const oC=OCT.map(v=>v.map(F=>({F,a:sc(F)}))), sC=B7.map(v=>v.map(F=>({F,a:sc(F)}))), xC=Array.from({length:17},(_,i)=>i?{F:sext.get(i),a:sc(sext.get(i))}:null);
function bucket(A,key){const m=new Map;for(const z of A){const q=key(z);const a=m.get(q)||[];a.push(z);m.set(q,a);}return m;}
function dot(U,V){let s=0;for(let i=0;i<U.length;i++)s+=U[i]*V[i];return s%17}
for(const mask of masks){
 stats.mask++;
 const S=Array.from({length:17},(_,i)=>i).filter(i=>mask&(1<<i));
 const W=S.map(s=>inv[S.reduce((q,a)=>s===a?q:q*mod(s-a)%17,1)]), Wi=W.map(a=>inv[a]);
 const sv=S.map(s=>xp[s].slice(0,5)), lagged=new Map;
 function lag(J){
  const key=J.join(',');if(lagged.has(key))return lagged.get(key);
  const T=J.map(i=>S[i]),x=S.map(s=>T.map((t,i)=>T.reduce((a,u,j)=>i===j?a:a*mod(s-u)*inv[mod(t-u)]%17,1)));
  lagged.set(key,x);return x;
 }
 const n0=stats.Yfit, n1=stats.lead5Fits, n2=stats.constant;
 // Choose top R8 up to global sign. Degree one normalized by slope coefficient (constant if slope zero).
 for(let d=1;d<=8;d++)for(const [A,B] of [[d,0],...Array.from({length:17},(_,i)=>[i,d])]){
  stats.top++;
  const R8=S.map(s=>(A+B*s)%17), a8=R8.map((u,i)=>u*W[i]%17), iz=R8.indexOf(0);
  if(a8.some(a=>a&&!oC[a].length))continue;
  stats.validTop++;
  const O=a8.map((a,i)=>a?i:-1).filter(i=>i>=0), C=O.slice(0,4), M=lag(C);
  const ys=a8.map((a,i)=>a?bucket(oC[a],v=>v.F[6]*R8[i]*Wi[i]%17):null);
  const I=O.slice().sort((i,j)=>ys[i].size-ys[j].size||i-j).slice(0,5), L=lag(I), to=[];
  function verifyY(j){
   if(j<5){for(const c of ys[I[j]].keys()){to[j]=c;verifyY(j+1);}return;}
   stats.Yinterpolants++;
   const Y=Array(11);
   for(const i of O){
    const v=dot(L[i],to);
    if(!ys[i].has(v))return;
    Y[i]=v;
   }
   if(iz>=0)Y[iz]=dot(L[iz],to);
   stats.Yfit++;
   const D=S.map(s=>mod((s-S[O[0]])*(s-S[O[1]])));
   for(let h=0;h<17;h++){
    const R7=D.map(a=>a*h%17);
    if(iz>=0&&(Y[iz]+10*R7[iz]*R7[iz])%17)continue;
    stats.rootShiftFit++;
    const Z=C.map(i=>mod((Y[i]+10*R7[i]*R7[i])*inv[R8[i]]));
    const R6=S.map((_,i)=>dot(M[i],Z));
    if(O.some(i=>mod(R8[i]*R6[i]-10*R7[i]*R7[i])!==Y[i]))continue;
    if(iz>=0&&!R7[iz]&&!R6[iz])continue;
    stats.fullLower++;
    const opt=a8.map((a,i)=>{
      if(a){
       const p=R7[i]*inv[8*R8[i]%17]%17;
       return bucket(ys[i].get(Y[i]).map(v=>v.a[p]),v=>v[5]*Wi[i]%17);
      }
      if(R7[i]){
       const p=R6[i]*inv[7*R7[i]%17]%17;
       return bucket(sC[W[i]*R7[i]%17].map(v=>v.a[p]),v=>v[5]*Wi[i]%17);
      }
      return bucket(xC[W[i]*R6[i]%17].a,v=>v[5]*Wi[i]%17);
    });
    const K=Array.from({length:11},(_,i)=>i).sort((i,j)=>opt[i].size-opt[j].size||j-i).slice(0,5), W5=lag(K);
    const ca=[], rows=[], F=[];
    function raw(j,sum){
     if(j<11){for(const a of rows[j]){F[j]=a;raw(j+1,sum+a[0]);}return;}
     stats.assignments++;
     if(sum%17!==3)return;
     stats.constant++;
     for(let c=4;c>=1;c--)for(let k=0;k<=c;k++){
       let z=0;for(let i=0;i<11;i++)z+=F[i][c]*sv[i][k];
       if(z%17){stats.firstFail[c]++;return;}
     }
     stats.passes++;
    }
    function fit5(j){
     if(j<5){for(const a of opt[K[j]].keys()){ca[j]=a;fit5(j+1);}return;}
     stats.lead5Trials++;
     for(let i=0;i<11;i++){const set=opt[i].get(dot(W5[i],ca));if(!set)return;rows[i]=set;}
     stats.lead5Fits++;raw(0,0);
    }
    fit5(0);
   }
  }
  verifyY(0);
 }
 stats.per.push([mask,stats.Yfit-n0,stats.lead5Fits-n1,stats.constant-n2]);
}
assert(stats.validTop===377&&stats.Yfit===2059611&&stats.assignments===30644&&stats.constant===728&&!stats.passes,'highest octic coupling');
return stats;
})(BAS,CATALOGUE);
console.log(JSON.stringify({BAS,octics:CATALOGUE.stats,stats},null,2));

