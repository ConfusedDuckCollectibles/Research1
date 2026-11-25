"use strict";
// Exact follow-through after verify-moser.cjs. Run node verify-reduction.cjs.
const check=(ok,id)=>{if(!ok)throw Error(id);};
const mod=(x,m)=>((x%m)+m)%m;
const key=w=>w.map(x=>mod(x,3)).join("");
const plus=(u,v)=>u.map((x,i)=>x+v[i]), sub=(u,v)=>u.map((x,i)=>x-v[i]);
const pq=([a,b,c,d])=>[6*(a*a+a*b+b*b+c*c+c*d+d*d)+10*(a*c+b*d)+5*(a*d+b*c),b*c-a*d];
const D=[];
for(let a=-10;a<=10;a++)for(let b=-10;b<=10;b++)for(let c=-10;c<=10;c++)for(let d=-10;d<=10;d++){
  const v=[a,b,c,d], [P,Q]=pq(v);if(P===54&&Q===0)D.push(v);
}
const counts={};for(const v of D)counts[key(v)]=(counts[key(v)]??0)+1;
const C=[[0,0,0,0],[0,0,1,1],[0,0,2,2],[1,1,0,0],[2,2,0,0]];
const table=C.map(a=>C.map(b=>counts[key(sub(a,b))]??0));
check(D.length===42&&JSON.stringify(table)===JSON.stringify([
 [18,3,3,3,3],[3,18,3,0,0],[3,3,18,0,0],[3,0,0,18,3],[3,0,0,3,18]
]),"third-lattice unit cosets");
check(Object.values(counts).filter(x=>x===1).length===12,"no additional multiplicity-two/three coset");
check(D.every(w=>w.slice(0,2).every(x=>x%3===0)||w.slice(2).every(x=>x%3===0)),"unit with two nonzero row residue blocks");
let diamond=[0,0]; // Fourth position in B's off coset, or in M.
for(const b of C.slice(1)){
  const neigh=D.map(w=>plus(b,w)).filter(a=>counts[key(a)]===3);
  check(neigh.length===21,"triple-contact off neighbors");
  for(let i=0;i<neigh.length;i++)for(let j=i+1;j<neigh.length;j++){
    const [a,c]=[neigh[i],neigh[j]];
    if(pq(sub(a,c)).join(",")!=="162,0")continue; // Three known diamond vertices, d(a,c)=sqrt(3).
    const d=sub(plus(a,c),b), k=key(d);
    check(k===key(b)||k==="0000","diamond fourth coset");
    diamond[+(k==="0000")]++;
  }
}
check(JSON.stringify(diamond)==="[72,12]","diamond gate");
const U={20:54,21:57,22:60,23:64,24:68}, small={1:0,2:1,3:3,4:5,5:7};
const out={};
for(let r=1;r<=5;r++){
 const pairs=[];for(let i=0;i<r;i++)for(let j=i+1;j<r;j++)pairs.push([i,j]);
 const hist={};
 for(let g=0;g<1<<pairs.length;g++){
  const edges=pairs.filter((e,i)=>g>>i&1);
  if(edges.length>small[r])continue;
  const A=Array.from({length:r},()=>Array(r).fill(false)),D=Array(r).fill(0);
  for(const [i,j] of edges){A[i][j]=A[j][i]=true;D[i]++;D[j]++;}
  if(D.some(d=>d<2)||pairs.some(([i,j])=>A[i].filter((a,k)=>a&&A[j][k]).length>=3))continue;
  const e=Array.from({length:1<<r},(_,z)=>edges.reduce((n,[i,j])=>n+((z>>i&1)&&(z>>j&1)),0));
  const pop=Array.from({length:1<<r},(_,z)=>z.toString(2).replace(/0/g,"").length);
  if(e.some((q,z)=>q===6&&pop[z]===4))continue;
  const t=Array(r).fill(0);
  function rec(m){
   if(m<r){for(let x=Math.max(0,5-D[m]);x<=3;x++){t[m]=x;rec(m+1);}return;}
   for(let z=1;z<1<<r;z++){
    let sum=0;for(let i=0;i<r;i++)if(z>>i&1)sum+=D[i]+t[i];
    if(73-sum+e[z]>U[25-pop[z]])return;
   }
   // The off-points with three contacts already belong to (1/3)M.
   const known=t.map(x=>x===3);
   let change=true;
   while(change){
    change=false;
    for(let i=0;i<r;i++)if(!known[i]){
      for(let j=0;j<r;j++)if(known[j]&&A[i][j])for(let k=j+1;k<r;k++)if(known[k]&&A[i][k]){
        if(A[j][k]||A[j].some((flag,w)=>flag&&A[k][w]&&known[w]&&w!==i)){known[i]=true;change=true;}
      }
    }
   }
   // A vertex with one contact occurs only as the centre of the gem.
   for(let i=0;i<r;i++)if(t[i]<=1) check(r===5&&edges.length===7&&D[i]===4&&t[i]===1,"unanticipated low contact");
   const unresolved=known.filter(x=>!x).length;
   const K=edges.length+"/"+unresolved;
   hist[K]=(hist[K]??0)+1;
   if(!unresolved) return;
   const center=D.indexOf(4);
   check(r===5&&edges.length===7&&unresolved===3&&center>=0,"non-gem residual");
   const rim=D.map((d,i)=>i).filter(i=>i!==center);
   const ends=rim.filter(i=>D[i]===2), inner=rim.filter(i=>D[i]===3);
   check(ends.length===2&&inner.length===2&&!A[ends[0]][ends[1]]&&A[inner[0]][inner[1]],"not a four-vertex rim path");
   check(ends.every(i=>known[i]&&t[i]===3)&&inner.concat(center).every(i=>!known[i]&&t[i]===2),"gem labels");
  }
  rec(0);
 }
 out[r]=hist;
}
const result={oneThirdUnits:D.length,table,diamond,offPatterns:out};
check(JSON.stringify(out)===JSON.stringify({1:{},2:{},3:{},4:{"4/0":3,"5/0":18},5:{"6/0":210,"5/0":12,"7/0":480,"7/3":60}}),"off-pattern census");
console.log(JSON.stringify(result));
