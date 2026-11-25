"use strict";
// Independent fixed-abscissa planar search and projective K10 atlas, F13.
const p=13, F=Array.from({length:p},(_,i)=>i), mod=x=>(x%p+p)%p;
const inv=F.map(x=>F.find(y=>mod(x*y)===1)||0);
const assert=(x,m)=>{if(!x)throw Error(m);};
const vec=(N,f)=>Array.from({length:N},(_,i)=>f(i));
const sum=x=>x.reduce((s,v)=>s+v,0);
const pop=n=>{let c=0;while(n){n&=n-1;c++;}return c;};
const subsets=(n,k)=>vec(1<<n,i=>i).filter(m=>pop(m)===k);
const unpack=m=>F.filter(i=>m>>i&1);
function canon(x){
  let m=1<<13;
  for(let a=1;a<13;a++)for(let b=0;b<13;b++)
    m=Math.min(m,x.reduce((m,t)=>m|(1<<mod(a*t+b)),0));
  return m;
}
function planar(n,bound,cap) {
  const C=new Set(subsets(13,n).map(m=>canon(unpack(m))));
  const good=[], nodes=Array(n+1).fill(0);
  for(const m of [...C].sort((a,b)=>a-b)){
    const x=unpack(m), d=x.map(a=>x.map(b=>inv[mod(a-b)])), h=[0,0];
    let firstLeaf=false;
    function go(scaled,mask){
      const i=h.length;nodes[i]++;
      if(i===n){assert(scaled,"all collinear");firstLeaf=true;return;}
      for(const v of (scaled?F:[0,1])){
        let dirs=mask,ok=true;
        const geometries=Array(13).fill(0);
        for(let j=0;j<i;j++){
          const q=mod((v-h[j])*d[i][j]);dirs|=1<<q;
          if(pop(dirs)>bound||++geometries[q]>=cap){ok=false;break;}
        }
        if(ok){h.push(v);go(scaled||v!==0,dirs);h.pop();}
        if(firstLeaf)return; // per-abscissa feasibility only, not a model count
      }
    }
    go(false,1);
    if(firstLeaf)good.push(m);
  }
  return {n,bound,cap,supportClasses:C.size,good,nodes:nodes.slice(2)};
}

const gates=[[4,3,2],[6,5,2],[10,9,2]].map(q=>planar(...q));
assert(gates.every(x=>x.good.length===0),"cap2 gates");
const b=m=>m===1?0:m%2?m:m-1;
const partitions=[];
for(let i=1;i<=13;i++){
 if(i===13)partitions.push([i]);
 for(let j=i;j<=13-i;j++){
  if(i+j===13)partitions.push([i,j]);
  for(let k=j;k<=13-i-j;k++)if(i+j+k===13)partitions.push([i,j,k]);
 }
}
const bad=partitions.filter(x=>sum(x.map(b))<11);
assert(JSON.stringify(bad)==="[[1,2,10],[1,4,8],[1,6,6]]","worst partitions");
const upgraded=partitions.map(x=>sum(x.map(n=>b(n)+ +([4,6,10].includes(n)))));
assert(Math.min(...upgraded)===11,"upgraded partition directions");
const counts=[1,2,3,4,5].map(c=>{const high=3*c,need=Math.max(0,Math.ceil((213-63*c)/21));return [c,high,need,high+need];});
assert(counts.every((x)=>x[3]>=11),"rich pencils");
console.log(JSON.stringify({gates,partitionCount:partitions.length,bad,upgradedMin:Math.min(...upgraded),counts}));
