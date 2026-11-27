"use strict";
const F=Array.from({length:17},(_,i)=>i),m=x=>(x%17+17)%17, iv=F.map(a=>F.find(b=>m(a*b)===1)||0);
const I=17, G=[10,12,13,14,16,I], A0=[[0,0],[1,0],[0,1]];
const code=([x,y])=>17*x+y;
const dir=(a,b)=>a[0]===b[0]?I:m((a[1]-b[1])*iv[m(a[0]-b[0])]);
const mask=S=>S.reduce((a,x)=>a|1<<x,0), noncol=(C,p)=>new Set(C.map(q=>dir(p,q))).size===C.length;
const choose=(A,k)=>k?A.flatMap((a,i)=>choose(A.slice(i+1),k-1).map(L=>[a,...L])):[[]];
const orbit=new Set();
for(const a of F)for(const b of F)for(const c of F)for(const d of F){
 if(!m(a*d-b*c)||[a,b,c,d].find(x=>x)!==1)continue;
 orbit.add(mask(G.map(t=>{
   const [x,y]=t===I?[0,1]:[1,t], u=m(a*x+b*y),v=m(c*x+d*y);
   return u===0?I:m(v*iv[u]);
 })));
}
const pts=F.flatMap(x=>F.map(y=>[x,y])).filter(u=>!A0.some(v=>code(u)===code(v)));
const ans=[];
function isHex(C){
 const sx=m(iv[6]*C.reduce((q,p)=>q+p[0],0)),sy=m(iv[6]*C.reduce((q,p)=>q+p[1],0));
 const X=C.map(([x,y])=>[m(x-sx),m(y-sy)]), R=new Set(X.map(code));
 const d=(a,b)=>m(a[0]*b[1]-a[1]*b[0]);
 if(!X.every(p=>R.has(code(p.map(x=>m(-x))))))return false;
 return X.some(u=>X.some(v=>d(u,v)!==0&&R.has(code([m(u[0]+v[0]),m(u[1]+v[1])]))));
}
for(const k of [5,6]){
 let nodes=0,N5=0,N6=0,match=0,hex=0,accepted=0,witness;
 const rest=F.slice(1,16), pal=choose(rest,k-3);
 for(const W of pal){
  const D=mask([0,16,I,...W]), C5=[],C6=[], C=A0.slice();
  const avail=pts.filter(p=>noncol(C,p)&&C.every(q=>D>>dir(p,q)&1));
  function go(j){
   nodes++;
   if(C.length===5) C5.push(C.map(code));
   if(C.length===6){C6.push(C.map(code));if(isHex(C))hex++;return;}
   for(let i=j;i<avail.length;i++){
    const p=avail[i];
    if(!noncol(C,p)||!C.every(q=>D>>dir(p,q)&1))continue;
    C.push(p);go(i+1);C.pop();
   }
  }
  go(0);N5+=C5.length;N6+=C6.length;
  if(C5.length||C6.length){
   if(!orbit.has(D))throw Error("different palette orbit");
   accepted++;
   for(const P of C5){ if(C6.some(Q=>P.every(v=>Q.includes(v))))match++; }
   if(!witness)witness=C6[0];
  }
 }
 ans.push({k,pal:pal.length,nodes,N5,N6,match,hex,accepted,witness});
}
const out={orbit:orbit.size,ans};
console.log(JSON.stringify(out));
