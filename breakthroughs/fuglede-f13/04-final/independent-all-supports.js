"use strict";
// Independent direct residual-height gate. No preset slope palettes.
const p=13,F=Array.from({length:p},(_,x)=>x),M=x=>(x%p+p)%p;
const inv=F.map(x=>F.find(y=>M(x*y)===1)??0);
const pc=Array.from({length:8192},(_,x)=>{let n=0;while(x){n++;x&=x-1;}return n;});
const bit=x=>1<<x, decode=m=>F.filter(x=>m&bit(x));
const assert=(q,s)=>{if(!q)throw Error(s);};
const masks=[111,187,215,311,427];
function canonical(x){let b=8192;for(let a=1;a<13;a++)for(let c=0;c<13;c++)b=Math.min(b,x.reduce((m,u)=>m|bit(M(a*u+c)),0));return b;}
const supports=[...new Set(pc.flatMap((n,m)=>n===6?[canonical(decode(m))]:[]))].sort((a,b)=>a-b);
function all(mask) {
  const x=decode(mask), h=[0,0], reduced=[];
  const d=x.map(u=>x.map(v=>inv[M(u-v)]));
  function go(on,palette){
    const n=h.length;
    if(n===6){if(on&&pc[palette]===5)reduced.push({h:h.slice(),palette});return;}
    for(const v of on?F:[0,1]){
      let out=palette;
      for(let j=0;j<n;j++)out|=bit(M((v-h[j])*d[n][j]));
      if(pc[out]<=5){h.push(v);go(on||v!==0,out);h.pop();}
    }
  }
  go(false,1);
  if(!reduced.length)return null;
  const by=new Map;
  for(const {h,palette} of reduced)for(let a=1;a<13;a++)for(let b=0;b<13;b++){
    const ff=h.map((v,i)=>M(a*v+b*x[i]));
    const P=decode(palette).reduce((bmask,v)=>bmask|bit(M(a*v+b)),0);
    if(!by.has(P))by.set(P,new Set);
    by.get(P).add(ff.join(','));
  }
  const stable=[];
  for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
    if([...by.values()].every(fs=>new Set([...fs].map(k=>{const v=k.split(',').map(Number);return M(v[j]-v[i]);})).size===1))
      stable.push([x[i],x[j]]);
  }
  const multiplicities={};
  for(const v of by.values())multiplicities[v.size]=(multiplicities[v.size]||0)+1;
  assert(reduced.length===1,"reduced representative count");
  const mult=[187,427].includes(mask)?2:1;
  assert(by.size===156/mult&&Object.keys(multiplicities).length===1&&multiplicities[mult]===156/mult,"palette count");
  return {mask,x,reduced,multiplicities,stable};
}
const results=supports.map(all).filter(Boolean);
assert(supports.length===14&&results.map(o=>o.mask).join(",")===masks.join(","),"all six-support classes");
for(const o of results){
  const pair=[187,427].includes(o.mask)?[3,5]:o.x.slice(0,2);
  assert(o.stable.some(p=>p.join(',')===pair.join(',')),"fixed stable pair");
}
console.log(JSON.stringify(results));
