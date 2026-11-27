"use strict";
// Complete centered capped residue catalogues by one even/odd interpolation method.
const p=17,F=Array.from({length:p},(_,i)=>i),m=x=>(x%p+p)%p,inv=F.map(a=>F.find(b=>m(a*b)===1)||0);
const pow=(x,k)=>k?m(x*pow(x,k-1)):1;
const U=F.slice(0,9).map(t=>m(t*t));
const dot=(a,b)=>m(a.reduce((s,x,i)=>s+x*b[i],0));
const vals=c=>F.map(t=>c.reduceRight((r,x)=>m(r*t+x),0));
const arr=(n,f)=>Array.from({length:n},(_,i)=>f(i));
function getCatalog(cap){const min=1,max=10;
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
    for(const [lead,u] of uni)D[d][lead].push({q,peak:Math.max(...Z),c:c.map((x,j)=>m(x*pow(u,j)))});
    r.yes++;
   }
   dfs(0,0);
  }
  function dfs(i){if(i===n){check();return;}for(const v of F){O[i]=v;dfs(i+1);}}
  dfs(0);stats.push(r);
 }
}
const info={cap,stats,n:D.map(X=>X.reduce((s,R)=>s+R.length,0)),h:D.map(X=>X.reduce((h,R)=>{for(const z of R)h[z.q]=(h[z.q]??0)+1;return h;},{}))};
return {info,D};}
const B4=getCatalog(4), B6=getCatalog(6), evalp=vals, summary=()=>({c4:B4.info,c6:B6.info});

const good=[0,1,2,9,16,17],S=F.filter(s=>!good.includes(s)),em=new Set(),W=S.map(s=>inv[S.reduce((a,t)=>a*(s===t?1:m(s-t))%17,1)]);

const bin=Array.from({length:11},(_,n)=>Array.from({length:n+1},(_,j)=>j===0||j===n?1:0));
for(let n=2;n<=10;n++)for(let j=1;j<n;j++)bin[n][j]=m(bin[n-1][j-1]+bin[n-1][j]);
const caches=[0,1].map(()=>Array.from({length:11},()=>Array.from({length:17},()=>Array(17))));
function shifted(d,a,t,ex){
 const M=caches[ex][d][a];if(M[t])return M[t];
 return M[t]=(ex?B6:B4).D[d][a].filter(z=>!ex||z.peak>=5).map(({q,c})=>({q,c:c.map((_,j)=>m(c.reduce((s,y,k)=>k<j?s:s+y*bin[k][j]*pow(t,k-j),0)))}));
}
function options(i,d,a,b){
 const A=m(a*W[i]),B=m(b*W[i]),ex=+(em.has(i));
 if(A)return shifted(d,A,m(B*inv[m(d*A)]),ex);
 if(B)return F.flatMap(t=>shifted(d-1,B,t,ex));
 const all=[];
 for(let j=1;j<=d-2;j++)for(let a=1;a<17;a++)for(const t of F)all.push(...shifted(j,a,t,ex));
 return all;
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
  if(!G){G=rows[u]={x:u,L:[],low:Infinity,high:0};p.push(G);}
  G.L.push(z);G.low=Math.min(G.low,z.q);G.high=Math.max(G.high,z.q);
 }
 return {rows,p,low:Math.min(...p.map(z=>z.low)),high:Math.max(...p.map(z=>z.high))};
}
function cnt(d){return {d,stages:0,gauges:0,lead:0,splits:Array(11).fill(0),pass:Array(11).fill(0),mods:0,energy:0,grid:0,models:[]};}
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
 H.stages++;
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
  const l=m(A[i]*W[i]);if(l>8||J.some(t=>!(em.has(t)?B6:B4).D[8][m(A[t]*W[t])].length))continue;
  for(const e of F)for(const h of F)profiles(8,A,S.map(s=>m((s-S[i])*(s-S[j])*(e+h*s))),H);
 }
 return H;
}

function paletteRep(n) {
 const P=F.concat(17),g0=good.join(), G=[];
 for(const a of F)for(const b of F)for(const c of F)for(const d of F){
  if((a||b||c||d)!==1||m(a*d-b*c)===0)continue;
  const f=P.map(s=>{
   const x=s===17?b:m(a+b*s),y=s===17?d:m(c+d*s);
   return x?m(y*inv[x]):17;
  });
  if(good.map(s=>f[s]).sort((u,v)=>u-v).join()===g0)G.push(f);
 }
 const B=[],seen=new Set(),out=[];
 function dfs(j){
  if(B.length===n){
   const v=B.join();if(seen.has(v))return;
   const Q=[...new Set(G.map(f=>B.map(i=>S.indexOf(f[S[i]])).sort((a,b)=>a-b).join()))];
   Q.forEach(x=>seen.add(x));out.push({idx:B.slice(),size:Q.length});return;
  }
  for(let i=j;i<=S.length-(n-B.length);i++){B.push(i);dfs(i+1);B.pop();}
 }
 dfs(0);return {G:G.length,out,seen:seen.size};
}

const LOW=(()=>{

function p17(x){return (x%17+17)%17;}
const F=[...Array(17).keys()], I=F.map(x=>F.find(y=>p17(x*y)===1)||0), good=[0,1,2,9,16,17],S=F.filter(s=>!good.includes(s));
const W=S.map(s=>I[S.reduce((a,t)=>a*(s===t?1:p17(s-t))%17,1)]), id=F.slice(0,12);
function reps(){const P=[...F,17], G=[];
 for(const a of F)for(const b of F)for(const c of F)for(const d of F){
  if((a||b||c||d)!==1||p17(a*d-b*c)===0)continue;
  const f=P.map(s=>{const x=s===17?b:p17(a+b*s),y=s===17?d:p17(c+d*s);return x?p17(y*I[x]):17;});
  if(good.map(s=>f[s]).sort((u,v)=>u-v).toString()===good.toString())G.push(f);
 }
 const seen=new Set(),out=[];
 for(let i=0;i<12;i++)for(let j=i+1;j<12;j++)for(let k=j+1;k<12;k++){
  const B=[i,j,k];if(seen.has(B.toString()))continue;
  const os=[...new Set(G.map(f=>B.map(i=>S.indexOf(f[S[i]])).sort((a,b)=>a-b).toString()))];
  os.forEach(x=>seen.add(x));out.push({B,size:os.length});
 }
 if(G.length!==12||seen.size!==220)throw Error("orbit coverage");
 return out;
}
const evalp=(P,x)=>P.reduceRight((a,b)=>(a*x+b)%17,0);
const squares=F.filter(c=>c>0&&F.some(x=>x*x%17===c));
function gate(cats) {
 const out=reps().map(({B,size})=>{
  const C=id.filter(i=>!B.includes(i));
  const L=B.map(i=>C.map(j=>C.reduce((a,k)=>j===k?a:p17(a*(S[i]-S[k])*I[p17(S[j]-S[k])]),1)));
  let top7=0,zero7=0,top8=0;const W7=[];
  function test7(a) {
   const v=S.map(x=>evalp(a,x));if(C.some(i=>!v[i]))return;
   top7++;
   const Y=id.map(i=>p17(W[i]**2*v[i]**4));
   const T=L.map(l=>p17(l.reduce((q,t,j)=>q+7*t*Y[C[j]],0)));
   if(B.some((i,j)=>v[i]? T[j]!==12*Y[i]%17&&T[j]!==6*Y[i]%17
                         : !squares.includes(T[j]*I[p17(-15)]%17)))return;
   if(B.some(i=>!v[i]))zero7++;
   W7.push({a,T});
  }
  for(let degree=0;degree<=3;degree++){
   for(let k=0;k<17**degree;k++){
    let h=k,q=[];for(let j=0;j<degree;j++){q.push(h%17);h=Math.floor(h/17);}
    test7([...q,1]);
   }
  }
  for(const a of F)for(const b of F)for(const c of F){
   if(!a&&!b&&!c)continue;
   const v=S.map((s,i)=>W[i]*p17(a+b*s+c*s*s)%17);
   if(C.some(i=>v[i]&&!cats.d8low.includes(v[i]))||
      B.some(i=>v[i]&&!cats.d8high.includes(v[i])))continue;
   top8++;
  }
  return {B:B.map(i=>S[i]),size,top7,pass7:W7.length,zero7,top8,W7};
 });
 return out;
}

function solveLow(record,cats) {
 const B=record.B.map(s=>S.indexOf(s)),C=id.filter(i=>!B.includes(i)),chosen=C.slice(0,6),targets=id.filter(i=>!chosen.includes(i));
 const L=id.map(i=>chosen.map(j=>chosen.reduce((a,k)=>j===k?a:p17(a*(S[i]-S[k])*I[p17(S[j]-S[k])]),1)));
 const aS=S.map(s=>p17((s-S[C[0]])*(s-S[C[1]])));
 const pref=Array(5).fill(0);let count=0, c5hit=0, expT=0, constant=0,grids=0;
 const pow=(x,n)=>x**n%17;
 const bin=Array.from({length:8},(_,n)=>Array.from({length:n+1},(_,j)=>j===0||j===n?1:0));
 for(let j=2;j<=7;j++)for(let r=1;r<j;r++)bin[j][r]=(bin[j-1][r-1]+bin[j-1][r])%17;
 function shift(P,t){return id.slice(0,8).map(j=>p17(P.reduce((s,c,n)=>s+(j<=n?c*bin[n][j]*pow(t,n-j):0),0)));}
 const M=F.map(a=>F.map(b=>p17(a*b)));
 const pg=Array.from({length:8},(_,j)=>S.map(s=>pow(s,j)));
 const adj=F.map(a=>F.map(b=>F.map(c=>S.map((s,i)=>p17(aS[i]*(a+b*s+c*s*s))))));
 const lowCache=F.map((a)=>F.map(b=>(!a?null:shift(cats.c4[1][a][0].c,p17(b*I[p17(7*a)])))));
 const highCache=F.map((a)=>F.map(b=>(!a?null:cats.c6[1][a].filter(z=>z.peak>=5).map(z=>shift(z.c,p17(b*I[p17(7*a)]))))));
 const sixCache=F.map((a)=>F.map(b=>(!a?null:shift(cats.c6[0][a][0].c,p17(b*I[p17(6*a)])))));
 for(const w of record.W7){
   const A=S.map(s=>evalp(w.a,s)),c7=A.map((x,i)=>x*W[i]%17);
   const TL=id.map(i=>A[i]?7*W[i]**2*A[i]**3%17:0);
   const Z=id.filter(i=>A[i]===0);
   for(const a of F)for(const b of F)for(const c of F) {
    const E=adj[a][b][c];count++;
    if(Z.some(i=>p17(-15*E[i]**2)!==w.T[B.indexOf(i)]))continue;
    expT++;
    const r5=C.map(i=>p17(TL[i]+15*I[A[i]]*E[i]**2));
    const y=L.map(K=>p17(K.reduce((s,x,i)=>s+x*r5[i],0)));
    if(C.some((i,j)=>y[i]!==r5[j]))continue;
    const P=id.map(i=>{
      const a=c7[i],b=M[W[i]][E[i]],c=M[W[i]][y[i]];
      if(!B.includes(i))return lowCache[a][b];
      if(!a)return sixCache[b][c];
      return highCache[a][b].find(z=>z[5]===c);
    });
    if(P.some(x=>!x))continue;
    c5hit++;
    let ok=true;
    for(let j=4;j>=1;j--){
     for(let h=0;h<=j;h++)if(p17(P.reduce((s,z,i)=>s+z[j]*pg[h][i],0))){ok=false;break;}
     if(!ok)break;
     pref[j]++;
    }
    if(!ok)continue;
    if(p17(P.reduce((s,z)=>s+z[0],0))!==5)continue;
    constant++;
    const R=P.map(z=>F.map(t=>evalp(z,t)));
    const Q=F.flatMap(x=>F.map(y=>R.reduce((s,z,i)=>s+z[p17(y-S[i]*x)],0)));
    if(Q.every(z=>z===22||z===39)){
     if(Q.filter(z=>z===39).length!==34)throw Error('Boolean size');
     grids++;
    }
   }
 }
 return {B:record.B,count,expT,c5hit,pref,constant,grids};
}


return {gate,solveLow};
})();
function checkCatalog(){
 function gcd(a,b){return b?gcd(b,a%b):a;}
 const low=B4.D.slice(6,8),high=B6.D.slice(6,8);
 function demand(x,msg){if(!x)throw Error(msg);}
 demand(B4.info.n.slice(1,7).every(x=>x===0)&&B6.info.n.slice(1,6).every(x=>x===0),'lowest degrees');
 for(const a of F.slice(1)){
   demand(low[1][a].length===1,'one low centered septic');
   const tr=low[1][a][0];
   demand(tr.q===34&&tr.c[5]===m(7*a**3),'low invariant coefficient');
   const hr=high[1][a].filter(z=>z.peak>=5);
   demand(hr.length===2&&hr.every(v=>v.q===68)&&
     hr.map(z=>m(z.c[5]*inv[m(a**3)])).sort((a,b)=>a-b).toString()==='6,12','two high centered septics');
   demand(high[0][a].length===1&&high[0][a][0].q===68&&high[0][a][0].peak>=5,'one high centered sextic');
 }
 for(const d of [7,9])demand(gcd(d,16)===1,'field gauge');
 return {c4:low,c6:high,d8low:F.filter(a=>B4.D[8][a].length),d8high:F.filter(a=>B6.D[8][a].some(z=>z.peak>=5))};
}
function run(start=0,stop=19,mode='all',out=console.log){
 if(!['all','low','high'].includes(mode)||start<0||stop>19||stop<start)throw Error('usage: node script [start] [stop] [all|low|high]');
 const cats=checkCatalog(), g=paletteRep(3),pre=LOW.gate(cats);
 if(g.G!==12||g.seen!==220||g.out.reduce((s,x)=>s+x.size,0)!==220||g.out.length!==19)throw Error('orbit sizes');
 const ranked=pre.map(({W7,...x})=>x);
 out(JSON.stringify({S, W, cap4:B4.info.n, cap6:B6.info.n,
   rich:B6.D.map(X=>X.reduce((s,C)=>s+C.filter(z=>z.peak>=5).length,0)),
   first:start,last:stop,mode}));
 const got=[];
 for(let i=start;i<stop;i++){
  const r=g.out[i],info=ranked[i],mark=r.idx.map(i=>S[i]);
  if(mark.toString()!==info.B.toString()||info.size!==r.size||info.top8)throw Error('low gates');
  const line={mark,size:r.size};
  if(mode!=='high') {
    const x=LOW.solveLow(pre[i],cats);
    if(x.grids)throw Error('actual low model survived');
    line.pre=info;line.low=x;
  }
  if(mode!=='low'){
    em.clear();r.idx.forEach(i=>em.add(i));
    for(const d of [9,10]){
      const H=d===9?solve9():solve10();
      if(H.grid)throw Error('actual high model survived');
      line[d]=H;
    }
  }
  got.push(line);out(JSON.stringify(line));
 }
 return got;
}
if(typeof module!=='undefined'&&require.main===module)run(Number(process.argv[2]??0),Number(process.argv[3]??19),process.argv[4]??'all');
