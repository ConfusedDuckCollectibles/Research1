// Nine complete ambient-order controls independent of the Dreamer example fixtures.
// Enumerate integer sets containing 0, exact cyclotomic-zero cliques, and the separate gcd/phase predicate.
// Cases (p,N,k,sets,spectral,spectra): (2,12,2,165,27,57) (2,20,2,969,125,265) (2,36,2,6545,729,2385) (2,18,2,680,0,0) (2,16,3,6435,139,288) (3,18,2,24310,256,520) (3,36,1,595,160,720) (3,40,1,741,0,0) (5,30,1,23751,1296,4991). Not a proof beyond these moduli.
console.log((()=>{
const mod=(a,b)=>(a%b+b)%b,gcd=(a,b)=>b?gcd(b,a%b):a, fac=(n)=>Array.from({length:n},(_,i)=>i+1).filter(i=>n%i===0);
const div=(A,B)=>{let r=A.slice(),q=Array(Math.max(0,r.length-B.length+1)).fill(0);for(let i=q.length-1;i>=0;i--){q[i]=r[i+B.length-1];if(q[i])for(let j=0;j<B.length;j++)r[i+j]-=q[i]*B[j];}return[q,r];};
const cache=new Map;function phi(n){if(cache.has(n))return cache.get(n);const a=Array(n+1).fill(0);a[0]=-1;a[n]=1;let q=a;for(const d of fac(n))if(d!==n){let rr;[q,rr]=div(q,phi(d));if(rr.some(Boolean))throw Error('bad Phi');}cache.set(n,q);return q;}
function eachComb(N,k,cb){let a=[0];function go(s){if(a.length===k){cb(a);return;}for(let t=s;t<=N-(k-a.length);t++){a.push(t);go(t+1);a.pop();}}go(1);}
function audit(p,N,k){
 const m=p**k;let P=1;while(N%(P*p)===0)P*=p;const q=N/P,vp=x=>{let i=0;while(x%p===0){++i;x/=p;}return i;};
 let aCount=0,sCount=0,sets=0,exactCount=0;
 eachComb(N,m,A=>{
  sets++; const z=Array(N).fill(false);const phiN=phi(N);
  for(let d=1;d<N;d++){const r=Array(N).fill(0);for(const a of A)r[d*a%N]++;z[d]=div(r,phiN)[1].every(x=>x===0);}
  let found=0;
  const vals=new Set;for(let i=0;i<m;i++)for(let j=0;j<i;j++)vals.add(vp(A[i]-A[j]));
  const steps=[...vals].sort((x,y)=>x-y).map(v=>v+1);
  const delta=steps.map(v=>{let d=0,t=p**(v-1);for(let i=0;i<m;i++)for(let j=0;j<i;j++)if((A[i]-A[j])%t===0)d=gcd(d,(A[i]-A[j])/t);return d;});
  let candidate=[0];
  function walk(start,pred){
   if(candidate.length===m){if(!pred){
     found++;sCount++;
     const B=candidate.map(x=>x%P),Ap=A.map(x=>x%P);
     if(new Set(B).size!==m||new Set(Ap).size!==m)throw Error('projection collision');
     const phiP=phi(P);
     for(let i=0;i<m;i++)for(let j=0;j<i;j++){
       const dd=candidate[i]-candidate[j],pOrd=P/gcd(P,dd),qOrd=q/gcd(q,dd),idx=steps.findIndex(l=>p**l===pOrd);
       if(idx<0||delta[idx]%qOrd)throw Error('phase classification failure');
       const r=Array(P).fill(0);for(let a of Ap)r[mod((B[i]-B[j])*a,P)]++;
       if(div(r,phiP)[1].some(Boolean))throw Error('projection not spectral');
     }
    }else exactCount++;
    return;
   }
   for(let t=start;t<=N-(m-candidate.length);++t){
    let ok=true;for(const s of candidate){const dd=t-s;if(!pred){if(!z[dd]){ok=false;break;}}
     else {const pOrd=P/gcd(P,dd),qOrd=q/gcd(q,dd),idx=steps.findIndex(l=>p**l===pOrd);if(idx<0||delta[idx]%qOrd){ok=false;break;}}}
    if(ok){candidate.push(t);walk(t+1,pred);candidate.pop();}
   }
  }
  walk(1,false); if(found)aCount++;
  if(steps.length===k&&p**steps[k-1]<=P){
   let prev=exactCount;walk(1,true);if(exactCount-prev!==found)throw Error('reverse phase count '+[p,N,k,A,exactCount-prev,found]);
  } else if(found)throw Error('tree');
 });
 return {p,N,k,sets,spectralSets:aCount,spectra:sCount,phaseModelSpectra:exactCount};
}
return [[2,12,2],[2,20,2],[2,36,2],[2,18,2],[2,16,3],[3,18,2],[3,36,1],[3,40,1],[5,30,1]].map(x=>audit(...x));
})());
