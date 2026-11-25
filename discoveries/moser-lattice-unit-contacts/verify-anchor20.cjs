"use strict";
// Complete exact 5/7-gem attachment gate for the unique 20/54 core.
const check=(x,s)=>{if(!x)throw Error(s);};
const a=(x,y)=>[x[0]+y[0],x[1]+y[1]], s=(x,y)=>[x[0]-y[0],x[1]-y[1]], k=(x,n)=>x.map(z=>z*n);
const m=(x,y)=>[x[0]*y[0]+33n*x[1]*y[1],x[0]*y[1]+x[1]*y[0]];
const sq=x=>m(x,x), zero=x=>x[0]===0n&&x[1]===0n;
function gcd(x,y){if(x<0n)x=-x;if(y<0n)y=-y;while(y)[x,y]=[y,x%y];return x;}
function norm(x,y,d){
 const g=[...x,...y,d].reduce(gcd,0n)*(d<0n?-1n:1n);
 return [x.map(z=>z/g),y.map(z=>z/g),d/g];
}
const add=(u,v)=>norm(a(k(u[0],v[2]),k(v[0],u[2])),a(k(u[1],v[2]),k(v[1],u[2])),u[2]*v[2]);
const sc=(u,n,d=1n)=>norm(k(u[0],n),k(u[1],n),u[2]*d);
const sub=(u,v)=>add(u,sc(v,-1n));
const norm2=u=>a(k(sq(u[0]),3n),sq(u[1]));
const isdist=(u,v,L2)=>{const d=sub(u,v);return zero(s(norm2(d),[432n*L2*d[2]*d[2],0n]));};
const rk=(u,sgn)=>norm(s(u[0],k(u[1],sgn)),a(k(u[0],3n*sgn),u[1]),2n*u[2]);
const key=u=>[...u[0],...u[1],u[2]].join(",");
function basis([a,b,c,d],den){
 return norm([12n*a+6n*b+10n*c+5n*d,-d],[18n*b+15n*d,2n*c+d],den);
}
// Input coordinates scaled as X=12 Re, Y=36 Im / sqrt(3).
function ringRatio(X,Y,D){
 const Db=[D[0],-D[1]],dd=m(D,Db)[0];return norm(m(X,Db),m(Y,Db),dd);
}
function unitCenter(o,U,V){
 const u=sub(U,o),v=sub(V,o),det=s(m(u[0],v[1]),m(v[0],u[1]));
 if(zero(det))return null;
 const nu=norm2(u),nv=norm2(v);
 // nu,u and nv,v have coefficient-vector denominators u[2],v[2].
 const D=k(det,6n*u[2]*v[2]);
 const X=s(k(m(nu,v[1]),v[2]),k(m(nv,u[1]),u[2]));
 const Y=k(s(k(m(u[0],nv),u[2]),k(m(v[0],nu),v[2])),3n);
 if(!zero(s(a(k(sq(X),3n),sq(Y)),k(sq(D),432n))))return null;
 return add(o,ringRatio(X,Y,D));
}
const A=[[1,0],[2,0],[0,1],[1,1],[0,2]],B=A.slice(0,4);
const results=[];
const triUnit=(p,q)=>{const x=p[0]-q[0],y=p[1]-q[1];return x*x+x*y+y*y===1;};
const W=A.flatMap(x=>A.map(y=>basis(x.concat(y).map(BigInt),1n)));
let control=0;
for(let i=0;i<W.length;i++)for(let j=0;j<i;j++)control+=isdist(W[i],W[j],1n);
check(new Set(W.map(key)).size===25&&control===72,"positive 25/72 control");
for(const reflected of [false,true])for(const conjugate of [false,true]){
 const P=A.flatMap(([a,b])=>B.map(([u,v])=>{
   if(reflected) [u,v]=[u,-u-v];
   return basis((conjugate?[3*a+5*u,3*b+5*v,-3*u,-3*v]:[3*a,3*b,3*u,3*v]).map(BigInt),3n);
 }));
 const core=new Set(P.map(key));
 check(core.size===20,"core size");
 const coreEdges=new Set();
 for(let i=0;i<20;i++)for(let j=0;j<i;j++)if(isdist(P[i],P[j],1n))coreEdges.add(i+","+j);
 check(coreEdges.size===54,"core not 20/54");
 const want=new Set(["9,6"]);
 for(let i=0;i<20;i++)for(let j=0;j<i;j++){
   const a=i>>2,b=i&3,c=j>>2,d=j&3;
   if((a===c&&triUnit(B[b],B[d]))||(b===d&&triUnit(A[a],A[c])))want.add(i+","+j);
 }
 check(want.size===54&&[...want].every(x=>coreEdges.has(x)),"core not labelled Cartesian product plus angle edge");
 const centers=new Map();let sourceTriples=0,onCore=0;
 for(let h=0;h<20;h++)for(let i=h+1;i<20;i++)for(let j=i+1;j<20;j++){
   const c=unitCenter(P[h],P[i],P[j]);
   if(c===null)continue;
   check([h,i,j].every(t=>isdist(c,P[t],1n)),"circumcenter verification");
   sourceTriples++;
   if(core.has(key(c))){onCore++;continue;}
   centers.set(key(c),c);
 }
 const contact=c=>P.filter(x=>isdist(c,x,1n)).length;
 const ends=[...centers.values()].filter(c=>contact(c)===3);
 let antipodal=0,nonoverlap=0,correct=0, best=-1;
 for(let i=0;i<ends.length;i++)for(let j=0;j<i;j++)if(isdist(ends[i],ends[j],4n)){
   antipodal++;
   const A=ends[i],E=ends[j],C=sc(add(A,E),1n,2n),R=sub(A,C);
   check(isdist(C,A,1n)&&isdist(C,E,1n),"radial");
   for(const sign of [-1n,1n]){
    const v1=add(C,rk(R,sign)),v2=add(C,rk(rk(R,sign),sign));
    const G=[A,v1,v2,E,C];
    if(new Set(G.map(key)).size!==5||G.some(v=>core.has(key(v))))continue;
    nonoverlap++;
    let e=0;for(let f=0;f<5;f++)for(let g=0;g<f;g++)e+=isdist(G[f],G[g],1n);
    check(e===7,"not a gem");
    const rights=G.map(contact), total=54+e+rights.reduce((u,v)=>u+v,0);
    if(total>best)best=total;
    if(rights.join()==="3,2,2,3,2") correct++;
   }
 }
 results.push({reflected,conjugate,sourceTriples,onCore,allExteriorCenters:centers.size,ends:ends.length,antipodal,nonoverlap,correct,best});
}
check(results.length===4&&results.every(r=>r.sourceTriples===306&&r.onCore===302&&r.ends===4&&r.antipodal===0),
      "unexpected core frontier counts");
console.log(JSON.stringify({control,results}));
