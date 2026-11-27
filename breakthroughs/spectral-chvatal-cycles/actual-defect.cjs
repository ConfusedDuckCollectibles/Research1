"use strict";
const c2=x=>x*(x-1)/2,c3=x=>x*(x-1)*(x-2)/6, A=(n,f)=>Array.from({length:n},(_,i)=>f(i));
function parts(r,m=3){if(!r)return [[]];return A(r-m+1>0?r-m+1:0,j=>m+j).flatMap(k=>parts(r-k,k).map(x=>[k,...x]));}
function tri(n){const S=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++)S.push([i,j,k]);return S;}
const msk=S=>S.reduce((z,t)=>z|1<<t,0);
function structure(n,cycles){
 const r=cycles.reduce((a,b)=>a+b,0),K=c2(n-1),a=c2(n-2),q=c2(n-4),L=K-2*r;
 const D=[];let off=0;
 for(const k of cycles){for(let i=0;i<k;i++)D.push([off+i,off+(i+1)%k,r+off+i]);off+=k;}
 const d=D.map(msk),U=tri(n).filter(t=>!d.includes(msk(t))),u=U.map(msk);
 return {n,r,K,a,q,L,D,d,U,u};
}
function factors(S,O){
 const {n,r,K,a,q,L,D,d,U,u}=S, agg=numeric=>O.make(numeric),J=A(n,i=>i),count=r;
 const old=O.gram(J,D), invAx=old;
 function adjusted(i,fn){
  const b=J.filter(j=>j!==i), G=O.gram(b,D.filter(x=>!x.includes(i))),v=fn(b,G);
  return A(n,j=>j===i?agg(0):v[b.indexOf(j)]);
 }
 const C=A(r,i=>adjusted(i,(b,G)=>O.solve(G,b.map(j=>O.sub(j>=r&&j<2*r?O.rat(K-1,K):agg(1),O.rat(a-D.filter(x=>x.includes(j)&&!x.includes(i)).length,a))))));
 const V=A(r,t=>{const i=t+r;return adjusted(i,(b,G)=>O.solve(G,b.map(j=>agg(-(n-2-D.filter(x=>x.includes(i)&&x.includes(j)).length)))));});
 const E=D.map(Di=>O.solve(invAx,J.map(j=>agg(+Di.includes(j)))));
 const P=d.map((v,i)=>d.map((w,j)=>O.sub(i===j?agg(1):(v&w)?agg(0):O.rat(1,q),O.rat(1,K))));
 return U.map((T,t)=>{
  const sums=(x)=>T.reduce((s,i)=>O.add(s,x[i]),agg(0));
  const w=C.map((x,i)=>O.sub(T.includes(i)?agg(0):O.add(O.rat(1,a),sums(x)),O.rat(1,K)));
  const h=w.reduce(O.add,agg(0));
  const v=V.map((x,i)=>T.includes(i+r)?agg(1):sums(x));
  const e=E.map(sums),b=d.map(x=>O.sub(x&u[t]?agg(0):O.rat(1,q),O.rat(1,K)));
  const eP=P.map(row=>row.reduce((z,x,i)=>O.add(z,O.mul(x,e[i])),agg(0)));
  return [w.concat([h],v,e,b,e),w.map(x=>O.mul(agg(2),x)).concat([O.mul(O.rat(4,L),h)],
          v.map(x=>O.mul(O.rat(1,K*(K-1)),x)),b,e,eP)];
 });
}
function Ogram(J,D,o){
 const m=J.length,a=c2(m-2), b=m-2;
 return J.map(i=>J.map(j=>o((i===j?a:0)+b-D.filter(x=>x.includes(i)&&x.includes(j)).length)));
}
function Gaussian(M,b,OS){const n=b.length,C=M.map((a,i)=>a.concat([b[i]]));
 const {sub,mul,div,zero}=OS;
 for(let j=0;j<n;j++){
  if(zero(C[j][j]))throw Error("positive pivot required");
  const v=C[j][j];
  for(let t=j;t<=n;t++)C[j][t]=div(C[j][t],v);
  for(let i=j+1;i<n;i++){const c=C[i][j];if(!zero(c))for(let t=j;t<=n;t++)C[i][t]=sub(C[i][t],mul(c,C[j][t]));}
 }
 const X=A(n,()=>null);
 for(let j=n-1;j>=0;j--){let a=C[j][n];for(let t=j+1;t<n;t++)a=sub(a,mul(C[j][t],X[t]));X[j]=a;}
 return X;
}
function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b)[a,b]=[b,a%b];return a;}
function R(a,b=1){a=BigInt(a);b=BigInt(b);if(!b)throw Error("division by zero");if(!a)return Z;if(b<0n){a=-a;b=-b;}const g=gcd(a,b);return [a/g,b/g];}
const Z=[0n,1n],One=[1n,1n],add=(a,b)=>R(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),sub=(a,b)=>R(a[0]*b[1]-b[0]*a[1],a[1]*b[1]),
    mul=(a,b)=>R(a[0]*b[0],a[1]*b[1]),div=(a,b)=>R(a[0]*b[1],a[1]*b[0]),sq=a=>mul(a,a),far=a=>a[0]+"/"+a[1],zero=x=>!x[0];
const O={make:x=>R(x),rat:R,add,sub,mul,div,zero,gram:(J,D)=>Ogram(J,D,x=>R(x)),solve:(M,b)=>Gaussian(M,b,O)};
const F=1n<<26n,FF=F*F,FS=Number(F),SCALE=1n<<22n,DIV=Number(FF/SCALE),MAX=2**52;
function floor(a){const x=a[0]*F,y=a[1];return x>=0n?x/y:(x-y+1n)/y;}
function exact(n,cycles){
 const S=structure(n,cycles),G=factors(S,O),k=G[0][0].length;
 const maxL=A(k,()=>0n),maxR=A(k,()=>0n);
 const Q=G.map(([x,y])=>[x,y].map((p,t)=>Float64Array.from(p,(a,j)=>{
  const z=floor(a),v=z<0n?-z:z,m=t?maxR:maxL;if(v>m[j])m[j]=v;
  if(v>=2n**31n)throw Error("coefficient bound");return Number(z);
 })));
 const B=maxL.reduce((s,x,i)=>s+x*maxR[i],0n);
 if(B>BigInt(MAX))throw Error("dot integer overflow");
 let sum=0n,chunk=0,l=0,peak=0;
 for(let i=0;i<Q.length;i++){
  const x=Q[i][0];
  for(let j=0;j<=i;j++)if(S.u[i]&S.u[j]){
   const y=Q[j][1];let s=0;
   for(let h=0;h<k;h++)s+=x[h]*y[h];
   const t=Math.round(s/DIV),d=(i===j?1:2)*t*t;
   if(d>MAX)throw Error("single square bound");
   if(chunk+d>MAX){sum+=BigInt(chunk);chunk=0;}
   chunk+=d;l+=(i===j?1:2);peak=Math.max(peak,Math.abs(t));
  }
 }
 sum+=BigInt(chunk);
 const e=add(R(maxL.concat(maxR).reduce((s,v)=>s+v,BigInt(k)),FF),R(1n,2n*SCALE));
 const a=R(sum,SCALE*SCALE),b=mul(R(l),sq(e)),u=R(S.q-1,S.q),v=R(S.a-6-6*(n-4),S.a-6),target=sq(mul(u,v)),remain=sub(sub(target,a),b),
 pass=v[0]>0n&&remain[0]>0n&&sub(sq(remain),mul(R(4),mul(a,b)))[0]>0n;
 return {n,cycles,locked:l,peak,sum:sum.toString(),err:far(e),pass,ratio:Number(sq(remain)[0])/Number(sq(remain)[1]),upper:(Math.sqrt(Number(sum))/Number(SCALE)+Math.sqrt(l)*Number(e[0])/Number(e[1]))/(Number(u[0])/Number(u[1])*(Number(v[0])/Number(v[1])))};
}
const done={22:[5],23:[5,6],24:[5,6,7,8],25:[5,6,7,8,9,10]};
function cases(n){const C=[];for(let r=5;2*r+1<=n;r++)if(!(done[n]??[]).includes(r))for(const c of parts(r))C.push(c);return C;}
const low=Number(process.argv[2]??11),high=Number(process.argv[3]??25);
const out=[];
for(let n=low;n<=high;n++)for(const c of cases(n))out.push(exact(n,c));
console.log(JSON.stringify(out));

