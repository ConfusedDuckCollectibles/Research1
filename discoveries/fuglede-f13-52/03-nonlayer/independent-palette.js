// Independent palette/rank check for the fixed size-52 B shadow.
// Enumerates residual pairs directly, spans them, gauges before branching.
const p=13,m=x=>(x%p+p)%p,inv=Array.from({length:p},(_,x)=>x?Array.from({length:p},(_,y)=>y).find(y=>x*y%p===1):0);
function ok(c,s){if(!c)throw Error(s);}
const points=[0,1,3,9,15,17,21,22,26,33,35,37,40,45,50,51,52,56,60,63,65,68,73,76,78,79,85,86,91,92,99,101,104,107,108,115,120,122,124,127,130,137,138,140,143,150,152,155,156,163,167,168];
const X=points.map(n=>n/13|0),Y=points.map(n=>n%13);
const n=X.length,goods=[4,6,7,8,9,10,11,12,13],indices=Array.from({length:n},(_,i)=>i);
const lines=goods.map(q=>Array.from({length:p},(_,j)=>indices.filter(i=>(q===13?X[i]:m(Y[i]-q*X[i]))===j)).map(J=>{
 const coord=i=>q===13?Y[i]:X[i];
 ok(J.length===4,"quartets");
 return {J,t:J.map(i=>m(coord(i)-coord(J[0])))};
}));
const combos=(a,k)=>{
 let output=[];function r(i,z){
  if(z.length===k){output.push(z.slice());return;}
  for(let j=i;j<=a.length-k+z.length;j++){z.push(a[j]);r(j+1,z);z.pop();}
 }r(0,[]);return output;
};
const v=Array.from({length:p},(_,i)=>i), letters=combos(v,4), mask=a=>a.reduce((z,j)=>z|1<<j,0);
const seen=new Set(),pals=[];
for(const a of letters)if(!seen.has(mask(a))){
 pals.push(a);for(let s=1;s<p;s++)for(let t=0;t<p;t++)seen.add(mask(a.map(x=>m(x*s+t))));
}
ok(pals.length===7&&seen.size===715,"all palette orbits");
function row(L,a,b) {
 const {J,t}=L,u=m(-(a*t[2]+b*t[3])*inv[t[1]]);
 const w=new Int8Array(n);[w[J[0]],w[J[1]],w[J[2]],w[J[3]]]=[m(-a-b-u),u,m(a),m(b)];
 return w;
}
function span(L,palette,globalStats){
 const {t}=L,s=new Set(palette);let found=null;
 for(let u=0;u<p;u++)for(let v=0;v<p;v++){
  if(u===0&&v===0)continue;
  // (u,v) are residual heights after making first two heights zero.
  const H=[0,0,u,v];const d=[];
  for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)d.push(m((H[j]-H[i])*inv[m(t[j]-t[i])]));
  if(!palette.some(shift=>d.every(z=>s.has(m(z+shift)))))continue;
  if(!found)found=[u,v];
  else if(m(found[0]*v-found[1]*u)){globalStats[2]++;return [];}
 }
 if(!found){globalStats[0]++;return [row(L,1,0),row(L,0,1)];}
 globalStats[1]++;return [row(L,found[1],-found[0])];
}
let localStats=[0,0,0];
const sys=lines.map(L=>pals.map(p=>L.flatMap(x=>span(x,p,localStats))));
const plain=lines.map(L=>L.flatMap(x=>[row(x,1,0),row(x,0,1)]));
const factors=Array.from({length:p},(_,a)=>Array.from({length:p},(_,b)=>m(a*b)));
function add(basis,w){
 const r=w.slice();
 for(let j=0;j<n;j++)if(r[j]){
  if(!basis[j]){
   const z=factors[inv[r[j]]];for(let i=j;i<n;i++)r[i]=z[r[i]];
   basis[j]=r;return 1;
  }
  const z=factors[r[j]],b=basis[j];
  for(let i=j;i<n;i++)r[i]=(r[i]+p-z[b[i]])%p;
 }return 0;
}
function extend(basis,dim,rows){
 const z=basis.slice();let k=dim;
 for(const w of rows){k+=add(z,w);if(k===n)break;}
 return [z,k];
}
function canonical(rows){
 let [r,k]=extend(Array(n).fill(null),0,rows);r=r.map(x=>x&&x.slice());
 for(let j=n-1;j>=0;j--)if(r[j]){
  for(let i=0;i<j;i++)if(r[i]?.[j]){
   const z=factors[r[i][j]];
   for(let f=j;f<n;f++)r[i][f]=(r[i][f]+p-z[r[j][f]])%p;
  }
 }
 return [r.filter(x=>x).map(x=>x.join()).join(";"),r.filter(x=>x)];
}
const locals=sys.map(h=>[...new Map(h.map(canonical)).values()]);
// Choose a different noncollinear gauge from the author's first x-adjacent points.
const anc=indices.findLast(i=>X[i]&&Y[i]);
const gauge=[0,2,anc].map(j=>{const w=new Int8Array(n);w[j]=1;return w;});
ok(X[anc]&&X[0]===X[2]&&Y[0]!==Y[2],"independent gauge");
const [base,r0]=extend(Array(n).fill(null),0,gauge);ok(r0===3,"gauge rank");
for(const Q of combos(Array.from({length:9},(_,j)=>j),3))
 ok(extend(base,3,Q.flatMap(j=>plain[j]))[1]===52,"all triple ranks 49 before gauge");
const stats={rounds:0,paths:0,leaves:0,unique:0,K:[],mPlane:[],two:{counts:{},n:0,min:1000}};
const hashes=new Set(), counts=[];
function witness(B,z) {
 const h=z.slice();
 for(let j=n-1;j>=0;j--)if(B[j]){
  let s=0;for(let i=j+1;i<n;i++)s+=B[j][i]*h[i];h[j]=m(-s);
 }
 return h;
}
function measures(h,planes=true) {
 const keys=new Set;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  const a=m(X[j]-X[i]),b=m(Y[j]-Y[i]),c=m(h[j]-h[i]);
  ok(a||b,"boolean projection");
  keys.add(a?p*m(b*inv[a])+m(c*inv[a]):p*p+m(c*inv[b]));
 }
 let mp=0;
 if(planes)for(const a of [[0,0,1],...v.map(x=>[0,1,x]),...v.flatMap(x=>v.map(y=>[1,x,y]))]){
  const F=new Int8Array(p);for(let i=0;i<n;i++)mp=Math.max(mp,++F[m(a[0]*X[i]+a[1]*Y[i]+a[2]*h[i])]);
 }return [keys.size,mp];
}
for(const Q of combos(Array.from({length:9},(_,j)=>j),2)){
 const [B,r]=extend(base,3,Q.flatMap(j=>plain[j]));
 const free=indices.filter(i=>!B[i]),w=Array(n).fill(0);
 stats.two.counts[free.length]=(stats.two.counts[free.length]||0)+1;
 function rec(j,started){
  if(j===free.length){
   if(started){stats.two.n++;stats.two.min=Math.min(stats.two.min,measures(witness(B,w),false)[0]);}return;
  }
  for(let c=0;c<(started?p:2);c++){w[free[j]]=c;rec(j+1,started||c!==0);}
 }rec(0,false);
}
for(let low=0;low<9;low++){
 const [B,k]=extend(base,3,plain[low]);
 for(const six of combos(Array.from({length:9},(_,j)=>j).filter(j=>j!==low),6)){
  stats.rounds++;
  function rec(u,R,k){
   stats.paths++;
   if(u===six.length){
    ok(k===51,"only one nongauge free at every retained leaf");
    stats.leaves++;
    const j=R.findIndex(x=>!x),z=Array(n).fill(0);z[j]=1;
    const h=witness(R,z), key=h.join();
    if(!hashes.has(key)){hashes.add(key);const[a,b]=measures(h);counts.push([a,b]);}
    return;
   }
   for(const P of locals[six[u]]){
    const [K,t]=extend(R,k,P);if(t!==52)rec(u+1,K,t);
   }
  }rec(0,B,k);
 }
}
stats.unique=hashes.size;stats.K=counts.map(t=>t[0]);stats.mPlane=counts.map(t=>t[1]);
ok(stats.rounds===252&&stats.unique===9&&counts.every(t=>t[0]===63&&t[1]===49),"independent palette leaf conclusion");
ok(stats.two.n===218&&stats.two.min===146,"two low independent census");
console.log(JSON.stringify({checks:"PASS",localStats,distinct:locals.map(a=>a.length),anchors:[0,2,anc],stats}));
