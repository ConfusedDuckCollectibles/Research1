// Independent actual-mask planar census for the size-52 nonlayer proof.
const p=13,m=x=>(x%p+p)%p,inv=Array.from({length:p},(_,x)=>x?Array.from({length:p},(_,y)=>y).find(y=>x*y%p===1):0);
const pop=b=>{let a=0;while(b){b&=b-1;a++;}return a;};
function requireThat(x,s){if(!x)throw Error(s);}
function DFS(goal,limit,cap,see) {
 const pt=[[0,0],[1,1]],levels=new Int8Array(p*p);
 levels[p]=2;let visited=0;
 function rec(x0,yused,dm){
  visited++;
  if(pt.length===goal){see(pt,dm);return;}
  for(let x=x0;x<=p-goal+pt.length;x++)next: for(let y=2;y<p;y++){
   if(yused>>y&1)continue;
   let mask=dm;const line=new Map();
   for(const [a,b] of pt){
    const s=m((y-b)*inv[x-a]),h=m(b-s*a),id=p*s+h;
    if((levels[id]||1)>=cap)continue next;
    mask|=1<<s;if(pop(mask)>limit)continue next;
    line.set(id,levels[id]);
   }
   for(const [k,v] of line)levels[k]=(v||1)+1;
   pt.push([x,y]);rec(x+1,yused|1<<y,mask);pt.pop();
   for(const [k,v] of line)levels[k]=v;
  }
 }rec(2,3,2);return visited;
}
function affineCanonical(v) {
 let best=null;
 for(let a=1;a<p;a++)for(let b=0;b<p;b++){
  let t="";for(let x=0;x<p;x++)t+=v[m(a*x+b)];
  if(best===null||t<best)best=t;
 }return best;
}
function histToList(a){return [...a].sort((x,y)=>x[0].localeCompare(y[0]));}
function complete(goal,maxD,cap,maxMultiplicity) {
 let total=0;const profiles=new Map,small=new Map;
 function see(pt,dm) {
  total++;if(pop(dm)<=8)small.set(dm,(small.get(dm)||0)+1);
  for(let d=0;d<=p;d++){
   const v=Array(p).fill(0);pt.forEach(([x,y])=>v[d===p?x:m(y-d*x)]++);
   if(Math.max(...v)>maxMultiplicity)continue;
   const a=affineCanonical(v);profiles.set(a,(profiles.get(a)||0)+1);
  }
 }
 const visits=DFS(goal,maxD,cap,see);
 return {goal,maxD,cap,total,visits,profiles:histToList(profiles),small};
}
const a=complete(10,9,4,2), b=complete(7,6,4,1);
let n9=0,n8=0;
DFS(9,7,4,()=>n9++);DFS(8,7,3,()=>n8++);
function pglCanon(mask) {
 const f=Array.from({length:14},(_,i)=>i).filter(i=>mask>>i&1);let least=1<<14;
 // Act on (coordinate1, coordinate2) so slope is coordinate2/coordinate1.
 for(let a=0;a<p;a++)for(let b=0;b<p;b++)for(let c=0;c<p;c++)for(let d=0;d<p;d++)if(m(a*d-b*c)){
  let mask=0;
  for(const s of f){const x=s===p?b:m(a+b*s),y=s===p?d:m(c+d*s);mask|=1<<(x?m(y*inv[x]):p);}
  least=Math.min(least,mask);
 }return least;
}
const eight=[...a.small].reduce((s,[,count])=>s+count,0);
const orbits=[...new Set([...a.small].map(([mask])=>pglCanon(mask)))];
requireThat(a.total===750&&b.total===3108&&n9===0&&n8===0,"independent planar gates");
requireThat(JSON.stringify(a.profiles)===JSON.stringify([["0000120211012",1200],["0001002201211",1650],["0001011020122",600],["0011011111111",3900]]),"ten-point exact census");
requireThat(!b.profiles.some(x=>x[0]==="0000110101111"),"excluded seven-point shadow");
requireThat(orbits.length===1&&orbits[0]===pglCanon(16383-[0,2,3,4,5,7].reduce((s,x)=>s|(1<<x),0)),"one K8 direction orbit");
console.log(JSON.stringify({checks:"PASS",ten:{total:a.total,visits:a.visits,profiles:a.profiles},seven:{total:b.total,visits:b.visits,profiles:b.profiles},zero:[n9,n8],eight:{count:eight,actualMasks:a.small.size,orbit:orbits[0]}}));
