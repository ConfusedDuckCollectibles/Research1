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
const stats=((bas)=>{
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

const masks=bas.models;
const model7=new Map;
function e(P,h,start=0){return P.slice(start).reduceRight((a,b)=>(a*h+b)%17,0);}
for(let a=1;a<17;a++){
 const list=cat[7].filter(c=>c[7]===a&&xp.every((_,x)=>e(c,x)<=7));
 if(list.length!==3)throw Error('cap septic');
 model7.set(a,new Map(list.map(c=>[c[5],c])));
}
const coeff6=Array.from({length:17},(_,a)=>Array.from({length:17},(_,h)=>{
 if(!a)return null;const c=sext.get(a);return [e(c,h),mod(c[4]+15*a*h*h)];
}));
const coeff7=Array.from({length:17},(_,a)=>{if(!a)return null;return new Map([...model7.get(a)].map(([k,c])=>[k,Array.from({length:17},(_,h)=>[e(c,h),mod(c[4]+5*k*h+a*h*h*h)])]));});
const S6={leading:0,rat:0,quartic:0,constant:0,pass:0,afterC4:0,firstRemaining:[0,0,0,0],degreeTwoWitness:null,per:[]}, S7={leading:0,preChoices:0,pre:0,rootFilters:0,translated:0,polyDivide:0,constant:0,pass:0,per:[]};
const KSEXT=16,KSEPT=15;
for(const mask of masks){
 const X=Array.from({length:17},(_,i)=>i).filter(s=>mask&(1<<s));
 const W=X.map(s=>inv[X.reduce((a,t)=>s===t?a:a*mod(s-t)%p,1)]), WI=W.map(a=>inv[a]);
 const WXP=X.map((s,i)=>xp[s].slice(0,7).map(x=>x*W[i]%p));
 const SX=X.map(s=>xp[s].slice(0,7));
 function dot(U,V){let v=0;for(let i=0;i<U.length;i++)v+=U[i]*V[i];return v%17;}
 function lag(indices){
  const Z=indices.map(i=>X[i]);return X.map(s=>Z.map((t,i)=>Z.reduce((a,u,j)=>i===j?a:a*mod(s-u)*inv[mod(t-u)]%17,1)));
 }
 function c4check(cs){
  for(let j=0;j<=4;j++){let v=0;for(let i=0;i<11;i++)v+=cs[i]*SX[i][j];if(v%17)return false;}
  return true;
 }

function remaining6(cl,sh){
 const Q=cl.map((a,i)=>{
   const P=sext.get(a),h=sh[i]; return [null,
   mod(2*P[2]*h+4*P[4]*h*h*h+6*a*pow(h,5)),
   mod(P[2]+6*P[4]*h*h+15*a*pow(h,4)),
   mod(4*P[4]*h+3*a*h*h*h)];
 });
 for(let j=3;j>=1;j--)for(let k=0;k<=j;k++){
  if(Q.reduce((v,u,i)=>v+u[j]*SX[i][k],0)%17)return j;
 }
 return 0;
}
 const base6=S6.rat, base7=S7.pre;
 for(const A of [1,3])for(let b=0;b<17;b++)for(let c=0;c<17;c++)for(let d=0;d<17;d++){
  const Rat=X.map(s=>mod(A*WI[0]+s*(b+s*(c+s*d))));
  if(Rat.includes(0))continue;
  S6.leading++;
  const al=W.map((w,i)=>w*Rat[i]%17);
  const Z=al.map((a,i)=>Rat[i]*sext.get(a)[4]*WI[i]%17);
  let z0=dot(W,Z),z1=X.reduce((s,x,i)=>s+Z[i]*WXP[i][1],0)%17;
  if(z0||z1)continue;
  S6.rat++;
  const H=Rat.map((a,i)=>inv[6*a%17]*X[i]*(X[i]-1)%17);
  for(let u=0;u<17;u++)for(let v=0;v<17;v++)for(let w=0;w<17;w++){
   S6.quartic++;
   let s=0;const chi=[], hts=[];
   for(let i=0;i<11;i++){
    const h=mod(H[i]*(w+X[i]*(v+X[i]*u))),[c0,c4]=coeff6[al[i]][h];
    s+=c0;chi[i]=c4;hts[i]=h;
   }
   if(s%17!==3)continue;
   S6.constant++;if(c4check(chi)){S6.afterC4++;const j=remaining6(al,hts);S6.firstRemaining[j]++;if(!j){S6.pass++;S6.degreeTwoWitness={X,al,hts};}}
  }
 }
 // Normalize leading polynomial by its highest nonzero coefficient, since seventh-power scaling is bijective.
 for(const [a,b,c] of [[1,0,0],...Array.from({length:17},(_,i)=>[i,1,0]),...Array.from({length:289},(_,i)=>[i%17,Math.floor(i/17),1])]){
  S7.leading++;const R=X.map(s=>(a+s*(b+s*c))%17);
  const al=R.map((r,i)=>r*W[i]%17), sep=al.map((q,i)=>q?i:-1).filter(i=>i>=0), sex=al.map((q,i)=>q?-1:i).filter(i=>i>=0);
  const Y0=Array(11).fill(0), Yset=al.map((q,i)=>q?new Map([...model7.get(q)].map(([c5,P])=>[R[i]*c5*WI[i]%17,c5])):null);
  const T=sep.slice(0,7), YL=lag(T), PL=lag(sep.slice(0,5)), D=X.map(s=>mod((s-X[sep[0]])*(s-X[sep[1]])));
  const out=Array(11);
  const constant=[]
  function fit(j){
   if(j<7){for(const h of Yset[T[j]].keys()){Y0[j]=h;fit(j+1);}return;}
   S7.preChoices++;
   const Ys=Array(11);
   for(const i of sep){
    const q=dot(YL[i],Y0);if(!Yset[i].has(q))return;
    Ys[i]=q;
   }
   S7.pre++;
   for(const i of sex){
    const q=dot(YL[i],Y0);
    Ys[i]=q;if(!q||pow(mod(-q*inv[KSEPT]),8)!==1)return;
   }
   S7.rootFilters++;
   for(let h0=0;h0<17;h0++)for(let h1=0;h1<17;h1++){
    const U=X.map((s,i)=>mod((h0+h1*s)*D[i]));
    if(sex.some(i=>mod(Ys[i]+KSEPT*U[i]*U[i])!==0))continue;
    S7.translated++;
    const P5=sep.map(i=>mod((Ys[i]+KSEPT*U[i]*U[i])*inv[R[i]]));
    if(sep.some((i,j)=>dot(PL[i],P5)!==P5[j]))continue;
    S7.polyDivide++;
    let C=0;const C4=[];
    for(const i of sep){
      const v=coeff7[al[i]].get(Yset[i].get(Ys[i]))[U[i]*inv[7*R[i]%17]%17];
      C+=v[0];C4[i]=v[1];
    }
    for(const i of sex){
      const s=dot(PL[i],P5), h=s*inv[6*U[i]%17]%17;
      const v=coeff6[U[i]*W[i]%17][h];C+=v[0];C4[i]=v[1];
    }
    if(C%17!==3)continue;
    S7.constant++;if(c4check(C4))S7.pass++;
   }
  }
  fit(0);
 }
 S6.per.push([mask,S6.rat-base6]);S7.per.push([mask,S7.pre-base7]);
}

assert(S6.afterC4===6&&S6.firstRemaining[3]===6&&S6.pass===0,'sextic gate');
assert(S7.polyDivide===48537&&S7.constant===2872&&S7.pass===0,'septic gate');
return {S6,S7};

})(BAS);
console.log(JSON.stringify({BAS,...stats},null,2));

