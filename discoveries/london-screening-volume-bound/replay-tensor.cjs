"use strict";
const p=1009, m=x=>(x%p+p)%p;
const check=(c,s)=>{if(!c)throw Error(s);};
const E=[[1,0,0],[0,1,0],[0,0,1]], z=x=>x.every(u=>u===0),
  dot=(a,b)=>a.reduce((s,u,i)=>s+u*b[i],0),
  cross=([a,b,c],[x,y,z])=>[b*z-c*y,c*x-a*z,a*y-b*x],
  add=(a,b)=>a.map((x,i)=>x+b[i]), mul=(Q,a)=>Q.map(row=>dot(row,a)),
  outer=(x,y)=>x.map(a=>y.map(b=>a*b));
const pairs=[[0,0],[0,1],[0,2],[1,1],[1,2],[2,2]];
const sym=pairs.map(([a,b])=>E.map((_,i)=>E.map((_,j)=>Number((a===i&&b===j)||(a===j&&b===i)))));
function rank(X){
 const A=X.map(r=>r.map(m));let k=0;
 for(let j=0;j<6;j++){
  const a=A.findIndex((r,i)=>i>=k&&r[j]!==0);
  if(a<0)continue;
  [A[k],A[a]]=[A[a],A[k]];const d=A[k][j];
  for(let v=k+1;v<A.length;v++)if(A[v][j]){
   const c=A[v][j];A[v]=A[v].map((x,i)=>m(x*d-c*A[k][i]));
  }k++;
 }
 return k;
}
function point(x){
 const y=E.map(e=>cross(x,e)).find(v=>!z(v)),t=cross(x,y),r=dot(x,x);
 const S1=outer(y,y).map((a,i)=>a.map((b,j)=>r*b-dot(y,y)*x[i]*x[j]));
 const S2=outer(t,t).map((a,i)=>a.map((b,j)=>r*b-dot(t,t)*x[i]*x[j]));
 const S3=outer(y,t).map((a,i)=>a.map((b,j)=>b+t[i]*y[j]));
 const equations=[];
 for(const S of [S1,S2,S3]){
  const v=mul(S,x);check(z(cross(x,v))&&S.reduce((a,row,i)=>a+row[i],0)===0,"nontrace-free syzygy");
  const cols=E.map(h=>add(cross(h,v),cross(x,mul(S,h))));
  for(const [i,j] of [[0,1],[0,2],[1,2]])equations.push(sym.map(T=>dot(T[i],cols[j])-dot(T[j],cols[i])));
 }
 const nullv=pairs.map(([i,j])=>x[i]*x[j]);
 check(equations.every(row=>dot(row,nullv)===0)&&rank(equations)===5,"local converse rank");
}
let points=0;
for(let x=-3;x<=3;x++)for(let y=-3;y<=3;y++)for(let z0=-3;z0<=3;z0++)if(x||y||z0){point([x,y,z0]);points++;}
function toy(floor){
 const Q=3,W=[1,4,11],top=9*Q,min=9*floor;
 const best=Array(top+1).fill(-1),cnt=Array(top+1).fill(0);
 function rec(i,b,f){
  if(i===9){
   for(let j=b;j<=top;j++)if(f>best[j]){best[j]=f;cnt[j]=1;}else if(f===best[j])cnt[j]++;
   return;
  }
  for(let v=floor;v<=Q;v++)rec(i+1,b+v,f+(i%3===2?0:v*W[Math.floor(i/3)]));
 }
 rec(0,0,0);
 const saturation=3*(2*Q+floor);
 for(let j=min;j<=top;j++){
  let left=j-min,S=W.reduce((x,y)=>x+2*floor*y,0);
  for(let i=2;i>=0;i--){const u=Math.min(left,2*(Q-floor));S+=u*W[i];left-=u;}
  check(S===best[j],"finite sorted budget weight");
  if(j>=saturation){
   let count=0;for(let a=floor;a<=Q;a++)for(let b=floor;b<=Q;b++)for(let c=floor;c<=Q;c++)if(a+b+c+6*Q<=j)count++;
   check(count===cnt[j],"finite radial equality family");
  }
 }
 return {floor,min,saturation,top,saturatedValues:best.slice(saturation),saturatedMaximizers:cnt.slice(saturation)};
}
const shellRational=(a,b)=>{let k=a,l=b;while(l){[k,l]=[l,k%l];}return [a/k,b/k];};
const shell=100**3-98**3;
const h=[shellRational(shell,1e6),shellRational(2*shell,3e6),shellRational(1e6+2*shell,3e6)];
check(JSON.stringify(h)==='[[7351,125000],[7351,187500],[69851,187500]]',"fraction identity");
console.log(JSON.stringify({points,exactLocalSymmetricMatrixRank:5,kernel:"radial span only",finiteAllocation:[toy(0),toy(1)],fractions:h},null,2));
