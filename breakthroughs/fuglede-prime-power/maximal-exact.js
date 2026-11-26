// Independent polynomial character-sum and tiling checks for sizes 4 through 60.
const gcd=(a,b)=>b?gcd(b,a%b):a;
function rem(A,B){
 const r=A.slice(); const d=B.length-1;
 for(let i=r.length-1;i>=d;--i){const c=r[i];if(c)for(let j=0;j<=d;j++)r[i-d+j]-=c*B[j];}
 return r.slice(0,d);
}
const exact=new Map();
function phi(n){
 if(exact.has(n))return exact.get(n);
 let ans=Array(n+1).fill(0);ans[n]=1;ans[0]=-1;
 for(let t=1;t<n;t++)if(n%t===0){
   const p=phi(t),d=p.length-1,q=Array(ans.length-d).fill(0),r=ans.slice();
   for(let i=q.length-1;i>=0;--i){q[i]=r[i+d];if(q[i])for(let j=0;j<=d;j++)r[i+j]-=q[i]*p[j];}
   if(r.some(Boolean))throw Error('division');
   ans=q;
 }
 exact.set(n,ans);return ans;
}
function verify(m){
 let p=2;while(m%p)p++;
 let N,A,L;
 if(p===m){
   N=6*p;A=[0,...Array.from({length:p-2},(_,i)=>i+2),p+1];L=[0,p];
 }else{
   let b=2;while(gcd(m,b)!==1)b++;
   N=b*m;A=[];
   for(let j=0;j<p;j++)for(let k=0;k<m/p;k++)A.push(b*j+p*k);
   L=Array.from({length:p},(_,j)=>j*m/p);
 }
 const tiling=Array(N).fill(0);
 for(let a of A)for(let k=0;k<N/m;k++)tiling[(a+k*m)%N]++;
 if(!tiling.every(x=>x===1))throw Error('tiling '+m);
 const zeros=[];
 for(let d=0;d<N;d++){
   const counts=Array(N).fill(0);for(const a of A)counts[a*d%N]++;
   zeros.push(!rem(counts,phi(N)).some(Boolean));
 }
 const ok=(a,b)=>zeros[(a-b+N)%N];
 for(const a of L)for(const b of L)if(a!==b&&!ok(a,b))throw Error('orthogonality '+m);
 for(let t=0;t<N;t++)if(!L.includes(t)&&L.every(a=>ok(t,a)))throw Error('extendible '+m);
 if(L.length>=m)throw Error('not incomplete '+m);
 return {m,N,maximalIncompleteSize:L.length};
}
const checks=Array.from({length:57},(_,i)=>verify(i+4));
console.log(JSON.stringify({checks,allPassed:true}));
