'use strict';
const p=17, mod=z=>(z%17+17)%17, all=(1<<17)-1;
const assert=(b,a)=>{if(!b)throw Error(a)};
const pop=n=>{let j=0;while(n){n&=n-1;j++;}return j;};
const bits=n=>{const a=[];while(n){a.push(31-Math.clz32(n&-n));n&=n-1;}return a;};
function affineMasks(size){
 const seen=new Uint8Array(1<<17), rs=[], lens=[];
 for(let m=0;m<=all;m++){
  if(seen[m]||pop(m)!==size)continue;
  const xs=bits(m), o=new Set;
  // Direct canonicalization under all 272 affine maps, no graph-orbit reuse.
  for(let a=1;a<17;a++)for(let b=0;b<17;b++){
   let t=0;for(const x of xs)t|=1<<((a*x+b)%17);
   o.add(t);
  }
  assert(m===Math.min(...o),'canonical minimum');
  for(const t of o)seen[t]=1;
  rs.push(m);lens.push(o.size);
 }
 return {rs,lens};
}
const D=affineMasks(12), P=affineMasks(6);
assert(D.rs.length===25&&D.lens.reduce((a,b)=>a+b,0)===6188,'x masks');
assert(P.rs.length===49&&P.lens.reduce((a,b)=>a+b,0)===12376,'slope masks');
const iv=Array.from({length:17},(_,i)=>{for(let x=1;x<17;x++)if(x*i%17===1)return x;return 0});
function rankAffine(vals){
 if(vals.length<2)return 0;
 const B=[];
 for(let i=1;i<vals.length;i++){
  const r=vals[i].map((x,j)=>mod(x-vals[0][j]));
  for(const {v,k} of B){
    const s=r[k]*iv[v[k]]%17;
    if(s)for(let l=k;l<12;l++)r[l]=mod(r[l]-s*v[l]);
  }
  const k=r.findIndex(x=>x!==0);
  if(k>=0)B.push({v:r,k});
 }
 return B.length;
}
const rows=[],stats={cells:0,nodes:0,nonempty:0,solutions:0,max:0,maxrank:0,sizes:Array(5).fill(0),ranks:Array(3).fill(0)};
for(const xMask of D.rs){
 const X=bits(xMask);
 for(const bMask of P.rs){
  stats.cells++;
  const sales=bits(all^bMask), L=new Uint8Array(17);for(const s of sales)L[s]=1;
  const products=Array.from({length:12},(_,i)=>Array.from({length:12},(_,j)=>{
    let out=0;for(const s of sales)out|=1<<mod(s*(X[j]-X[i]));return out;
  }));
  const y=new Int8Array(12);y.fill(-1);
  const cnt=Array.from({length:17},()=>new Uint8Array(17));
  const buckets=Array.from({length:13},()=>new Uint32Array(12));buckets[0].fill(all);
  let nodes=0;
  const vals=[];
  function assign(k,idx,h){
   const V=buckets[k+1];V.set(buckets[k]);y[idx]=h;V[idx]=0;
   const r=new Int8Array(17);
   for(const a of sales){
    const b=mod(h-a*X[idx]);r[a]=b;
    assert(++cnt[a][b]<=5,'premature cap violation');
   }
   let ok=true;
   for(let j=0;j<12;j++)if(y[j]<0){
    let z=products[idx][j];
    V[j]&=((z<<h)|(z>>>(17-h)))&all;
    for(const a of sales)if(cnt[a][r[a]]===5)V[j]&=~(1<<mod(r[a]+a*X[j]));
    if(!V[j])ok=false;
   }
   if(ok)visit(k+1);
   for(const a of sales)--cnt[a][r[a]];
   y[idx]=-1;
  }
  function visit(k){
   nodes++;
   if(k===12){
    const yy=[...y];
    const line=new Map;
    for(let i=0;i<12;i++)for(let j=0;j<i;j++){
      const s=mod((yy[i]-yy[j])*iv[mod(X[i]-X[j])]);
      assert(L[s],'forbidden pair');
      const b=mod(yy[i]-s*X[i]), key=17*s+b;
      line.set(key,(line.get(key)||0)+1);
      assert(line.get(key)<=10,'six-point line');
    }
    vals.push(yy);
    return;
   }
   const Q=buckets[k];let i=-1,best=18;
   for(let j=0;j<12;j++)if(y[j]<0){
    const n=pop(Q[j]);
    if(n<best){i=j;best=n;}
   }
   assert(i>=0&&best>=1,'invalid DFS');
   for(let H=Q[i];H;H&=H-1)assign(k,i,31-Math.clz32(H&-H));
  }
  assign(0,0,0);
  const af=rankAffine(vals), n=vals.length;
  assert(n<=4&&af<=2,'atlas counterexample: '+JSON.stringify({X,sales,vals,af}));
  stats.max=Math.max(stats.max,n);stats.maxrank=Math.max(stats.maxrank,af);
  stats.solutions+=n;stats.nodes+=nodes;
  stats.sizes[n]++;
  if(n){stats.nonempty++;stats.ranks[af]++;rows.push({xMask,bMask,n,af,heights:vals});}
 }
}
assert(stats.nonempty===207&&stats.solutions===353&&stats.max===4&&stats.maxrank===2,'atlas census');
const rank2=rows.filter(x=>x.af===2).map(x=>[x.xMask,x.bMask,x.n]);
console.log(JSON.stringify({stats,rank2,domains:D.rs,palettes:P.rs},null,2));
