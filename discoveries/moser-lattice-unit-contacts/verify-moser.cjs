"use strict";
// Pure integer replay. Run node verify-moser.cjs. Nothing downloaded.
const add=([a,b],[c,d])=>[a+c,b+d], sub=([a,b],[c,d])=>[a-c,b-d];
const mul=([a,b],[c,d])=>[a*c+33n*b*d,a*d+b*c], scale=([a,b],n)=>[n*a,n*b];
const sq=x=>mul(x,x), zero=([a,b])=>a===0n&&b===0n;
function gcd(a,b){ a=a<0n?-a:a; b=b<0n?-b:b; while(b) [a,b]=[b,a%b]; return a; }
const require=(yes,why)=>{if(!yes)throw Error(why);};
const pq=([a,b,c,d])=>[6*(a*a+a*b+b*b+c*c+c*d+d*d)+10*(a*c+b*d)+5*(a*d+b*c), b*c-a*d];
const differences=[];
for(let a=-6;a<=6;a++)for(let b=-6;b<=6;b++)for(let c=-6;c<=6;c++)for(let d=-6;d<=6;d++){
  if(!(a||b||c||d))continue;
  const [P,Q]=pq([a,b,c,d]);
  if(P<=24&&(24-P)**2>=33*Q*Q)differences.push([a,b,c,d]);
}
require(differences.length===204,"difference list");
const points=differences.map(v=>{
  const [a,b,c,d]=v.map(BigInt);
  const x=[12n*a+6n*b+10n*c+5n*d,-d], y=[18n*b+15n*d,2n*c+d];
  return {x,y,n:add(scale(sq(x),3n),sq(y))};
});
const centers=new Map();let pairs=0;
for(let i=0;i<points.length;i++)for(let j=i+1;j<points.length;j++){
  const u=points[i],v=points[j],det=sub(mul(u.x,v.y),mul(v.x,u.y));
  if(zero(det))continue;
  const D=scale(det,6n), X=sub(mul(u.n,v.y),mul(v.n,u.y)),Y=scale(sub(mul(u.x,v.n),mul(v.x,u.n)),3n);
  if(!zero(sub(add(scale(sq(X),3n),sq(Y)),scale(sq(D),432n))))continue;
  pairs++;
  const conj=[D[0],-D[1]],den=mul(D,conj)[0],coords=[...mul(X,conj),...mul(Y,conj),den];
  const g=coords.reduce(gcd,0n)*(den<0n?-1n:1n), k=coords.map(x=>x/g),key=k.join(",");
  if(!centers.has(key))centers.set(key,{k,adj:new Set(),pairs:0});
  const c=centers.get(key);
  c.adj.add(i); c.adj.add(j); c.pairs++;
}
const hist=[{},{}],cosets=new Map();let eq=0;
for(const {k,adj,pairs} of centers.values()){
  const [x0,x1,y0,y1,D]=k;
  const d=-x1,c2=y1+x1,b18=y0-15n*d;
  // Coordinates (a,b,c,d) as rational pairs, in that order.
  const coeff=[
    [36n*x0-12n*b18-180n*c2-180n*d,432n*D],
    [b18,18n*D],[c2,2n*D],[d,D],
  ];
  const inM=coeff.every(([p,q])=>p%q===0n);
  const h=hist[+!inM], n=adj.size+1;
  h[n]=(h[n]??0)+1;
  require(pairs===adj.size*(adj.size-1)/2,"complete pair incidence");
  if(!inM){
    require(n===3,"off-lattice high degree");
    const residues=coeff.map(([p,q])=>{
      const r=((3n*p)% (3n*q)+3n*q)%(3n*q);
      require(r%q===0n,"denominator three");
      return Number(r/q);
    });
    cosets.set(residues.join(""),residues);
    const [u,v]=[...adj].map(i=>differences[i]);
    for(const w of [u,v,u.map((x,i)=>x-v[i])])require(pq(w).join(",")==="18,0","non-equilateral triple");
    eq++;
  }
}
const out={differences:differences.length,pairs,centers:centers.size,hist,equilateral:eq,cosets:[...cosets.keys()].sort()};
require(JSON.stringify(out)===JSON.stringify({differences:204,pairs:2460,centers:30,hist:[{"18":18},{"3":12}],equilateral:12,cosets:["0011","0022","1100","2200"]}),"summary");
console.log(JSON.stringify(out));
