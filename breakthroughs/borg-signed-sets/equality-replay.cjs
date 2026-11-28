"use strict";
const ord = n => Array.from({length:n},(_,i)=>i);
function check(r,k,supports,weights){
 const points=ord(r*k), faces=[];
 for(const [idx,S] of supports.entries()){
  const W=weights?weights[idx]:1; if(!W)continue;
  function fill(i,b){
   if(i===r){faces.push({S,b,W});return;}
   if(S&(1<<i))for(let j=0;j<k;j++)fill(i+1,b|(1n<<BigInt(k*i+j)));
   else fill(i+1,b);
  }
  fill(0,0n);
 }
 if(!faces.length)return {opt:0,count:0,nonstar:0};
 const bits=ord(faces.length).map(i=>1n<<BigInt(i)), idx=new Map(bits.map((b,i)=>[b,i]));
 const adj=faces.map((u,i)=>faces.reduce((p,v,j)=>i!==j&&(u.b&v.b)?p|bits[j]:p,0n));
 const star=points.map(x=>faces.reduce((p,v,j)=>v.b&(1n<<BigInt(x))?p+v.W:p,0));
 const opt=Math.max(...star);let best=0,count=0,nonstar=0,enumMax=0;
 function eligible(chosen){
  let common=(1n<<BigInt(r*k))-1n;
  for(const i of chosen)common&=faces[i].b;
  if(common)return false;
  if(k!==2)throw Error("unexpected nonstar");
  const byInc=new Map();
  for(let g=0;g<r;g++){
   if(star[2*g]!==opt)continue;
   const key=supports.filter((S,i)=>(!weights||weights[i]>0)&&(S&(1<<g))).join(",");
   (byInc.get(key)??(byInc.set(key,[]),byInc.get(key))).push(g);
  }
  const selectedMask=chosen.reduce((x,i)=>x|bits[i],0n);
  for(const C of byInc.values()){
   if(C.length<3)continue;
   const cells=faces.filter(f=>f.S&(1<<C[0]));
   if(chosen.some(i=>!(faces[i].S&(1<<C[0]))))continue;
   const crep=f=>C.reduce((z,g,j)=>(f.b&(1n<<BigInt(2*g+1)))?z|(1<<j):z,0);
   const good=new Set(chosen.map(i=>crep(faces[i]))), x=(1<<C.length)-1;
   if(good.size*2!==1<<C.length||[...good].some(n=>good.has(x^n)))continue;
   if(faces.every((f,i)=>Boolean(selectedMask&bits[i])===Boolean((f.S&(1<<C[0]))&&good.has(crep(f)))))return true;
  }
  throw Error("nonstar not classified");
 }
 function search(P,X,q,w){
  if(!P&&!X){enumMax++;if(w>opt)throw Error("star bound");if(w===opt){count++;if(eligible(q))nonstar++;}best=Math.max(best,w);return;}
  let y=P|X,top=-1,pivot=0;while(y){const b=y&-y;y^=b;let z=P&adj[idx.get(b)],n=0;while(z){z&=z-1n;n++;}if(n>top){top=n;pivot=idx.get(b);}}
  let v=P&~adj[pivot];while(v){const b=v&-v;v^=b;const i=idx.get(b);search(P&adj[i],X&adj[i],[...q,i],w+faces[i].W);P^=b;X|=b;}
 }
 search((1n<<BigInt(faces.length))-1n,0n,[],0);
 if(best!==opt||!count)throw Error("no maximum star");
 return {opt,count,nonstar,enumMax};
}
const all=[];for(let mask=1;mask<128;mask++)all.push(check(3,2,ord(7).map(x=>x+1).filter(x=>mask&(1<<(x-1)))));
const summary=key=>all.reduce((s,x)=>s+x[key],0);
console.log(JSON.stringify({binaryThreeAllSupportFamilies:all.length,maximizers:summary("count"),nonstar:summary("nonstar"),maximalCliques:summary("enumMax"),
 fixtures:[
  check(3,3,[3,5,6]),
  check(4,2,[7,15,8,3],[2,3,1,0]),
  check(5,2,[7,15,23,8,16,24],[2,2,1,2,3,1]),
  check(4,3,[3,5,9,14],[2,3,1,1])
 ]},null,2));
