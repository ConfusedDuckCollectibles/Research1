"use strict";
// Complete centered capped residue catalogues by one even/odd interpolation method.
const p=17,F=Array.from({length:p},(_,i)=>i),m=x=>(x%p+p)%p,inv=F.map(a=>F.find(b=>m(a*b)===1)||0);
const pow=(x,k)=>k?m(x*pow(x,k-1)):1;
const U=F.slice(0,9).map(t=>m(t*t));
const dot=(a,b)=>m(a.reduce((s,x,i)=>s+x*b[i],0));
const vals=c=>F.map(t=>c.reduceRight((r,x)=>m(r*t+x),0));
const arr=(n,f)=>Array.from({length:n},(_,i)=>f(i));
const cap=4,min=1,max=10;
if(![4,5,6].includes(cap))throw Error("use cap 4,5 or 6");
const D=arr(11,()=>arr(17,()=>[]));
const CC=new Map();
function table(J){
 const w=J.join();if(CC.has(w))return CC.get(w);
 const B=J.map(i=>{
  let a=[1],d=1;
  for(const j of J)if(j!==i){
   const s=Array(a.length+1).fill(0);d=m(d*(U[i]-U[j]));
   a.forEach((q,k)=>{s[k]=m(s[k]-q*U[j]);s[k+1]=m(s[k+1]+q);});a=s;
  }
  return a.map(x=>m(x*inv[d]));
 });
 const T={P:arr(J.length,j=>B.map(a=>a[j])),E:U.map(t=>B.map(a=>a.reduceRight((r,x)=>m(r*t+x),0)))};
 CC.set(w,T);return T;
}
const stats=[];
for(let d=min;d<=max;d++){
 const k=d>>1,odd=d&1,picks=new Set();
 for(let a=1;a<17;a++){
  if(picks.has(a))continue;
  const uni=new Map();for(let u=1;u<17;u++){const b=m(a*pow(u,d));if(!uni.has(b))uni.set(b,u);}
  for(const b of uni.keys())picks.add(b);
  const O=Array(k+(odd?1:-1)).fill(0),p0=U.map(t=>odd?0:m(a*pow(t,k)));
  if(odd)O[k]=a;
  const n=O.length-odd,r={d,a,odd:0,viable:0,leaves:0,yes:0};
  function check(){
   r.odd++;
   const P=U.map((t,i)=>m(i*O.reduceRight((r,x)=>m(r*t+x),0)));
   const A=[arr(cap+1,j=>[j,j])];let weight=0;
   for(let t=1;t<=8;t++){
    let diff=m(2*P[t]);if(diff>8)diff-=17;
    const gap=Math.abs(diff);
    if(gap>cap)return;
    weight+=gap;
    A[t]=arr(cap-gap+1,j=>[m(inv[2]*gap+j-p0[t]),gap+2*j]);
   }
   if(weight>34)return;
   r.viable++;
   const J=F.slice(0,9).sort((a,b)=>A[a].length-A[b].length||a-b).slice(0,k).sort((a,b)=>a-b);
   const K=F.slice(0,9).filter(a=>!J.includes(a)),L=table(J),B=Array(k).fill(0);
   const tail=Array(k+1).fill(K.reduce((w,i)=>w+A[i][0][1],0));
   for(let i=k-1;i>=0;i--)tail[i]=tail[i+1]+A[J[i]][0][1];
   const box=A.map(q=>{const h=Array(17).fill(-1);for(const [v,j] of q)h[v]=j;return h;});
   function dfs(i,h){
    if(h+tail[i]>34)return;
    if(i<k){for(const [b,w] of A[J[i]]){B[i]=b;dfs(i+1,h+w);}return;}
    r.leaves++;
    for(const j of K){const x=box[j][dot(B,L.E[j])];if(x<0||(h+=x)>34)return;}
    if(h!==34)return;
    const c=Array(d+1).fill(0);L.P.forEach((v,j)=>c[2*j]=dot(B,v));O.forEach((x,j)=>c[2*j+1]=x);
    if(!odd)c[d]=a;
    const Z=vals(c);if(c[d]!==a||c[d-1]!==0||Math.max(...Z)>cap||Z.reduce((u,v)=>u+v)!==34)throw Error("catalog control");
    const q=Z.reduce((u,v)=>u+(v-2)**2,0);
    for(const [lead,u] of uni)D[d][lead].push({q,c:c.map((x,j)=>m(x*pow(u,j)))});
    r.yes++;
   }
   dfs(0,0);
  }
  function dfs(i){if(i===n){check();return;}for(const v of F){O[i]=v;dfs(i+1);}}
  dfs(0);stats.push(r);
 }
}
const info={cap,stats,n:D.map(X=>X.reduce((s,R)=>s+R.length,0)),h:D.map(X=>X.reduce((h,R)=>{for(const z of R)h[z.q]=(h[z.q]??0)+1;return h;},{}))};
const evalp=vals; const summary=()=>info;

const good=[0,1,2,9,16,17],S=F.filter(s=>!good.includes(s)),W=S.map(s=>inv[S.reduce((a,t)=>a*(s===t?1:m(s-t))%17,1)]);

const bin=Array.from({length:11},(_,n)=>Array.from({length:n+1},(_,j)=>j===0||j===n?1:0));
for(let n=2;n<=10;n++)for(let j=1;j<n;j++)bin[n][j]=m(bin[n-1][j-1]+bin[n-1][j]);
const caches=Array.from({length:11},()=>Array.from({length:17},()=>Array(17)));
function shifted(d,a,t){
 const M=caches[d][a];if(M[t])return M[t];
 return M[t]=D[d][a].map(({q,c})=>({q,c:c.map((_,j)=>m(c.reduce((s,y,k)=>k<j?s:s+y*bin[k][j]*pow(t,k-j),0)))}));
}
function options(i,d,a,b){
 const A=m(a*W[i]),B=m(b*W[i]);
 if(A)return shifted(d,A,m(B*inv[m(d*A)]));
 if(B)return F.flatMap(t=>shifted(d-1,B,t));
 if(d===9)return F.slice(1).flatMap(a=>F.flatMap(t=>shifted(7,a,t)));
 return [];
}
function tablefor(Dst,ord){
 const raw=ord.slice(0,Dst).map(i=>S[i]),key=raw.join();
 if(tables.has(key))return tables.get(key);
 const T=raw.map(a=>S.map(x=>raw.reduce((s,y)=>y===a?s:m(s*m(x-y)*inv[m(a-y)]),1)));
 tables.set(key,T.map(C=>F.map(v=>C.map(t=>t*v%17))));
 return tables.get(key);
}
const tables=new Map(),pidx=F.slice(0,12);
function grouped(V,j,i){
 const rows=Array(17).fill(null),p=[];
 for(const z of V){
  const u=m((z.c[j]??0)*inv[W[i]]);let G=rows[u];
  if(!G){G=rows[u]={x:u,L:[],low:55,high:0};p.push(G);}
  G.L.push(z);G.low=Math.min(G.low,z.q);G.high=Math.max(G.high,z.q);
 }
 return {rows,p,low:Math.min(...p.map(z=>z.low)),high:Math.max(...p.map(z=>z.high))};
}
function cnt(d){return {d,gauges:0,lead:0,splits:Array(11).fill(0),pass:Array(11).fill(0),mods:0,energy:0,grid:0,models:[]};}
function profiles(d,A,B,H){
 H.gauges++;
 const seats=S.map((_,i)=>options(i,d,A[i],B[i]));
 if(seats.some(Q=>!Q.length))return;
 H.lead++;
 function final(R){
  const chosen=Array(12),rs=R.map(U=>[Math.min(...U.map(x=>x.q)),Math.max(...U.map(x=>x.q))]),tails=Array.from({length:13},(_,j)=>[0,0]);
  for(let j=11;j>=0;j--)for(let t=0;t<2;t++)tails[j][t]=tails[j+1][t]+rs[j][t];
  function dfs(i,c,q){
   if(q+tails[i][0]>510||q+tails[i][1]<510)return;
   if(i===12){
    H.energy++;
    if(c%17!==5)return;
    H.mods++;
    const list=chosen.map(z=>evalp(z.c)), w=[];
    for(const x of F)for(const y of F){
     let v=0;for(let j=0;j<12;j++)v+=list[j][m(y-S[j]*x)];
     if(v!==22&&v!==39)return;
     w.push((v===39)?1:0);
    }
    if(w.reduce((a,b)=>a+b,0)!==34)throw Error("mass");
    for(const s of F.concat(17))for(const t of F){
     let a=0;for(const x of F)a+=s===17?w[17*t+x]:w[17*x+m(s*x+t)];
     const i=S.indexOf(s);if(a!==(i<0?2:list[i][t]))throw Error("reprojection");
    }
    H.grid++;H.models.push({c:chosen.map(z=>z.c),w:w.join("")});return;
   }
   for(const z of R[i]){chosen[i]=z;dfs(i+1,c+z.c[0],q+z.q);}
  }
  dfs(0,0,0);
 }
 function stage(j,X){
  if(j===0){final(X);return;}
  const Q=X.map((a,i)=>grouped(a,j,i)), f=Q.reduce((s,z)=>s+z.low,0), g=Q.reduce((s,z)=>s+z.high,0);
  if(f>510||g<510)return;
  const O=Q.map(z=>z.p.filter(b=>f-z.low+b.low<=510&&g-z.high+b.high>=510));
  if(O.some(x=>!x.length))return;
  const k=11-j,go=pidx.slice().sort((a,b)=>O[a].length-O[b].length), T=tablefor(k,go), y=Array(12), rem=go.slice(k),qL=[0],qH=[0];
  for(let i=11;i>=0;i--){qL.unshift(qL[0]+Q[go[i]].low);qH.unshift(qH[0]+Q[go[i]].high);}
  const R=Array.from({length:k+1},()=>new Int16Array(rem.length));
  function rec(l,a,b){
   if(a+qL[l]>510||b+qH[l]<510)return;
   if(l===k){
    H.splits[j]++;
    let hi=b,lo=a;
    for(let i=0;i<rem.length;i++){
     const u=rem[i], w=Q[u].rows[R[k][i]%17];if(!w)return;
     y[u]=w.L;hi+=w.high;lo+=w.low;
    }
    if(lo>510||hi<510)return;
    H.pass[j]++;stage(j-1,y);return;
   }
   const t=go[l];
   for(const e of O[t]){
    y[t]=e.L;
    const X=R[l],U=R[l+1],Z=T[l][e.x];
    for(let i=0;i<rem.length;i++)U[i]=X[i]+Z[rem[i]];
    rec(l+1,a+e.low,b+e.high);
   }
  }
  rec(0,0,0);
 }
 stage(d-2,seats);
}
function solve10(){
 const H=cnt(10);
 for(const a of [1,3])profiles(10,S.map(()=>a),S.map(()=>0),H);
 return H;
}
function solve9(){
 const H=cnt(9);
 for(const r of F.concat(17)){
  const A=S.map(s=>r===17?1:m(s-r)),[i,j]=pidx.filter(i=>A[i]);
  for(const e of F)profiles(9,A,S.map(s=>m(e*(s-S[i])*(s-S[j]))),H);
 }
 return H;
}
function solve8(){
 const H=cnt(8);
 for(const a of F)for(const b of F)for(const c of F){
  const A=S.map(s=>m(a+b*s+c*s*s)),J=pidx.filter(i=>A[i]);
  if(!J.length)continue;const [i,j]=J;
  const l=m(A[i]*W[i]);if(l>8||J.some(t=>!D[8][m(A[t]*W[t])].length))continue;
  for(const e of F)for(const h of F)profiles(8,A,S.map(s=>m((s-S[i])*(s-S[j])*(e+h*s))),H);
 }
 return H;
}
const mode=(process.argv[2]??"all");
console.log(JSON.stringify({S,W,catalog:summary()}));
for(const d of ["10","9","8","7"])if(mode==="all"||mode===d){
 const q=d==="10"?solve10():d==="9"?solve9():d==="8"?solve8():{d:7,top_bound:408,case:"all profiles must be septic"};if(q.grid)throw Error("counterexample to marked cap-four exclusion");console.log(JSON.stringify(q));
}
