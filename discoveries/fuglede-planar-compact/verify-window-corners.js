// Exact 2,460-corner check. Seed formulas from Tao Zhang, arXiv:2607.15632v1, §4.
function compactWindowCorners(log = console.log) {
 const range=n=>Array.from({length:n},(_,i)=>i),mod=(x,n)=>(x%n+n)%n;
 const V=[range(12).map(y=>[0,y]),range(12).map(x=>[x,x]),range(6).flatMap(k=>[[mod(4*k,12),mod(2*k,12)],[mod(4*k+1,12),mod(2*k+2,12)]])];
 const S=[0,0,0,2,1].flatMap((i,t)=>V[i].map(([x,y])=>{
  const X=mod(25*x+36*t,60);
  return [mod(2*X+15*y-23,60),mod(X+8*y,12)];
 }));
 const differences=Array.from({length:720},()=>new Set());
 for(const [x,y] of S)for(const [u,v] of S){
  const dx=x-u,dy=y-v,rx=mod(dx,60),ry=mod(dy,12);
  const sx=dx===rx?0:1,sy=dy===ry?0:1;
  differences[rx+60*ry].add(sx+2*sy);
 }
 let residues=0,checks=0;
 for(let rx=0;rx<60;rx++)for(let ry=0;ry<12;ry++){
  const allowed=differences[rx+60*ry];
  if(!allowed.size)continue;
  residues++;
  const x=[Math.ceil((-436-rx)/60),Math.floor((436-rx)/60)];
  const y=[Math.ceil((-98-ry)/12),Math.floor((98-ry)/12)];
  for(const h of x)for(const k of y) {
   checks++;
   if(![...allowed].some(s=>Math.abs(h+(s&1))<=7&&Math.abs(k+(s>>1))<=8))throw Error("Uncovered extreme");
  }
 }
 const witness=[[0,0],[5,11],[2,10],[9,11],[4,11],[11,9],[3,7],[7,10],[10,11],[1,10]];
 for(let i=0;i<witness.length;i++)for(let j=0;j<i;j++){
  const u=mod(witness[i][0]-witness[j][0],12),v=mod(witness[i][1]-witness[j][1],12);
  if(u===0||mod(u-v,12)===0||(mod(u,4)!==2&&mod(v,2)===0&&mod(u+v,3)===0))throw Error("Invalid packing witness");
 }
 const result={residues,checks,packingWitness:witness};
 log(result);return result;
}
compactWindowCorners();
