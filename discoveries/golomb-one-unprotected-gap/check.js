
function mod(x,p){return (x%p+p)%p}
function gcd(a,b){while(b){let u=a%b;a=b;b=u}return a}

function singerPower(p,e){
 let q=p**e;const deg=3*e,order=p**deg-1,Q=q*q+q+1;
 let coeff;
 const I=Array(deg).fill(0);I[0]=1;const X=Array(deg).fill(0);X[1]=1;
 function same(a,b){return a.every((x,i)=>x===b[i])}
 function mul(a,b){
   let t=Array(2*deg-1).fill(0);
   for(let i=0;i<deg;i++)for(let j=0;j<deg;j++)t[i+j]+=a[i]*b[j];
   for(let k=2*deg-2;k>=deg;k--){
     let z=mod(t[k],p);
     if(z)for(let l=0;l<deg;l++)t[k-deg+l]-=z*coeff[l];
   }
   return t.slice(0,deg).map(x=>mod(x,p));
 }
 function pow(a,n){let b=I;while(n){if(n%2)b=mul(b,a);a=mul(a,a);n=Math.floor(n/2);}return b;}
 let t=order,fac=[];for(let k=2;k*k<=t;k++)if(t%k===0){fac.push(k);while(t%k===0)t/=k;}if(t>1)fac.push(t);
 for(let k=1;k<=order;k++){
   if(k%p===0)continue;
   coeff=Array.from({length:deg},(_,i)=>Math.floor(k/p**i)%p);
   if(same(pow(X,order),I)&&fac.every(t=>!same(pow(X,order/t),I)))break;
 }
 if(!same(pow(X,order),I)||!fac.every(t=>!same(pow(X,order/t),I)))throw Error("no field");
 const Y=pow(X,q), Z=pow(X,q*q);
 let a=I,b=I,c=I; const D=[];
 for(let r=0;r<Q;r++){
   if(a.every((x,i)=>mod(x+b[i]+c[i],p)===0))D.push(r);
   a=mul(a,X);b=mul(b,Y);c=mul(c,Z);
 }
 if(D.length!==q+1)throw Error();
 const ds=new Set(D.flatMap(x=>D.filter(y=>y!==x).map(y=>mod(x-y,Q))));
 if(ds.size!==Q-1)throw Error();
 return [Q,D,coeff];
}
function forcing(A,k) {
  const n=A.length-1, all=(1n<<BigInt(n))-1n, big=1n<<BigInt(k), small=all^big;
  const span=new Map();for(let i=0;i<n;i++)for(let j=i+1;j<=n;j++)
    span.set(A[j]-A[i],(1n<<BigInt(j))-(1n<<BigInt(i)));
  const U=new Map();
  for(const [d,m] of span) U.set(d,(m&big)?[0n,small^(m&small)]:[m,0n]);
  const arrows=[];
  for(const [d,v] of U) {
    const u=U.get(d-1);
    if(u){
      let pos=(v[0]|u[1]), neg=(v[1]|u[0]);
      // Possible multiplicities of ±2 do not change the sign unless cancelled.
      // Within each transformed interval all coefficients share its sign.
      if(v[0]&&u[0]){pos=v[0]&~u[0];neg=u[0]&~v[0];}
      if(v[1]&&u[1]){pos=u[1]&~v[1];neg=v[1]&~u[1];}
      arrows.push([pos,neg]);
    }
  }
  let zero=0n,round=0;
  while(true){
    let next=zero;
    for(const [pos,neg] of arrows)if((pos&~zero)===0n) next|=neg;
    if(zero===next)break;
    zero=next;round++;
  }
  return [zero===small,round];
}

function shiftOne(A,k) {
  const row=[[0,false,-1,-1]];
  for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++)
    row.push([A[j]-A[i],i<=k&&k<j,i,j]);
  row.sort((a,b)=>a[0]-b[0]);
  let gap=Infinity;
  for(let i=1;i<row.length;i++)if(row[i][1]&&!row[i-1][1])gap=Math.min(gap,row[i][0]-row[i-1][0]-1);
  if(!Number.isFinite(gap))throw Error();
  const opt=A.map((x,i)=>x-(i>k?gap:0));
  if(row.slice(1).some((d,i)=>opt[d[3]]-opt[d[2]]-(i?opt[row[i][3]]-opt[row[i][2]]:0)<1))throw Error("not feasible");
  return {gap,opt};
}
function runExact(P) {
  let all=[0,0,0,0], residual=[],counts={};
  for(const [s,e] of P) {
    const p=s**e;const [m,D]=singerPower(s,e);
    const seen=new Set();const t=[0,0,0,0];
    for(let u=1;u<m;u++)if(gcd(u,m)===1){
      const a=D.map(x=>x*u%m).sort((x,y)=>x-y);
      const g=a.map((x,i)=>mod(a[(i+1)%a.length]-x,m));
      const large=Math.max(...g), cut=Math.max(...g.filter(x=>x!==large));
      const b=a[(g.indexOf(cut)+1)%a.length], A=a.map(x=>mod(x-b,m)).sort((x,y)=>x-y);
      const key=A.join();if(seen.has(key))continue;seen.add(key);
      const k=A.findIndex((v,i)=>A[i+1]-v===large);
      if(k<0)throw Error();
      const {gap,opt}=shiftOne(A,k);
      const isUnique=forcing(opt,k)[0];
      t[0]++;
      if(gap)t[1]++;
      if(isUnique)t[2]++;
      else {t[3]++;residual.push({p,m,A,k,gap,opt});}
    }
    counts[p]=t;all=all.map((v,i)=>v+t[i]);
  }
  return {all,counts,residual};
}

function checkPDS(A,m) {
  const S=new Set(A.flatMap(x=>A.filter(y=>x!==y).map(y=>mod(x-y,m))));
  if(A.length*(A.length-1)+1!==m||S.size!==m-1||S.has(0))throw Error("bad PDS");
}
function exactExamples(){
  let a=[0,4,5];checkPDS(a,7);
  if(JSON.stringify(shiftOne(a,0))!==JSON.stringify({gap:2,opt:[0,2,3]}))throw Error("q2 diag");
  const b=[0,18,23,37,58,62,75,86,177,179,187,199,230,233,259,260,266,275];
  checkPDS(b,307);
  let q17=shiftOne(b,7);
  const A=[0,1,11,19,26,28],M=31, target=[0,1,8,13,17,19], dual=new Map([[1,7],[2,3],[9,2],[10,4],[15,3]]);
  checkPDS(A,M);
  const n=A.length-1, r=[[0,0,0]];
  for(let i=0;i<n;i++)for(let j=i+1;j<=n;j++)r.push([A[j]-A[i],i,j]);
  r.sort((a,b)=>a[0]-b[0]);
  let blocked=Array(n).fill(false), weight=Array(n).fill(0), bound=0, positive=[];
  for(let i=1;i<r.length;i++){
    const [v,a,b]=r[i],[u,c,d]=r[i-1];
    if(target[b]-target[a]-target[d]+target[c]<1)throw Error("bad example");
    const y=dual.get(v)||0;
    const R=weight.map((_,k)=>Number(a<=k&&k<b)-Number(c<=k&&k<d));
    if(v===u+1)for(let k=0;k<n;k++) if(R[k]>0)blocked[k]=true;
    for(let k=0;k<n;k++)weight[k]+=R[k]*y;
    if(y)positive.push(R);
    bound+=y;
  }
  if(!blocked.every(x=>x)||!weight.every(x=>x===1)||bound!==target.at(-1))throw Error("bad certificate");
  function det(M){
    if(!M.length)return 1;
    return M[0].reduce((t,x,i)=>t+(x===0?0:(i%2?-1:1)*x*det(M.slice(1).map(r=>r.filter((_,j)=>j!==i)))),0);
  }
  if(det(positive)===0)throw Error("not unique");
  return {q2shift:2,q17,mod31:{length:bound,det:det(positive),blocked}};
}
const output={control:exactExamples(),sweep:runExact([[2,1],[3,1],[2,2],[5,1],[7,1],[2,3],[3,2],[11,1],[13,1],[2,4],[17,1],[19,1],[23,1],[5,2],[3,3],[29,1],[31,1],[2,5],[37,1],[41,1],[43,1],[47,1],[7,2]])};
if(output.sweep.all.join()!=="3398,30,3398,0")throw Error("not reproduced");
console.log(JSON.stringify(output,null,2));
