
"use strict";
function checkGridText(raw){
 const p=13, m=x=>(x%p+p)%p, i=[0,1,7,9,10,8,11,2,5,3,4,6,12];
 const rows=raw.split(/\n/).filter(x=>x.startsWith("grid ")).map(x=>x.trim().split(/\s+/).slice(1).map(Number));
 const E="0,1,2,3,10,11,12,13,14,16,17,20,28,30,31,56,58,59,73,79,86,88,89,93,94,96,103,109,123,124,126,151,152,154,156,162,165,166,168";
 const out={rows:rows.length,distinctRows:new Set(rows.map(r=>r.join())).size,uniformCount8:0,peak8:0,canonicalBase:0,centerFailures:0};
 for(const R of rows){
  const P=R.map(z=>[Math.floor(z/p),z%p]), S=new Set(R);
  if(S.size!==39)throw Error("not Boolean 39");
  const blocks=Array.from({length:14},(_,q)=>{
    const A=Array(p).fill(0);
    for(const [x,y] of P)A[q===p?x:m(y-q*x)]++;
    return A;
  });
  if(blocks.filter(A=>A.every(z=>z===3)).length===8)out.uniformCount8++;
  if(Math.max(...blocks.flat())===8)out.peak8++;
  const C=P.filter(([u,v])=>P.every(([x,y])=>S.has(p*m(2*u-x)+m(2*v-y))));
  if(C.length!==1){out.centerFailures++;continue;}
  const [u,v]=C[0], T=P.map(([x,y])=>[m(x-u),m(y-v)]);
  let best=null;
  for(const [a,b] of T)for(const [c,d] of T){
    const z=i[m(a*d-b*c)];
    if(!z)continue;
    const V=T.map(([x,y])=>p*m(z*(d*x-c*y))+m(z*(-b*x+a*y))).sort((l,r)=>l-r);
    if(best===null){best=V;continue;}
    const j=V.findIndex((l,j)=>l!==best[j]);
    if(j!==-1&&V[j]<best[j])best=V;
  }
  if(best.join()===E)out.canonicalBase++;
 }
 if(out.rows!==60||out.uniformCount8!==60||out.peak8!==60||out.canonicalBase!==60||out.centerFailures!==0)throw Error("grid audit failed");
 return out;
}
if(typeof module!=="undefined"&&require.main===module){
  if(process.argv.length!==3)throw Error("usage: node grids.js k8-fixture.txt");
  console.log(JSON.stringify(checkGridText(require("fs").readFileSync(process.argv[2],"utf8")),null,2));
}
