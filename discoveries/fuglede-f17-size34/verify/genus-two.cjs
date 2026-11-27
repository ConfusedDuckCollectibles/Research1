'use strict';
// Exact F17 smooth genus-2 bound in the range relevant to 34.
// All smooth genus-2 curves are hyperelliptic. If a rational branch value
// exists, send it to infinity and normalize degree-five lead to 1. If none
// and N>=33, a split value exists; send one to infinity and use monic sextic.
const p=17,F=[...Array(p).keys()],M=x=>(x%p+p)%p;
const sq=new Set(F.slice(1).map(x=>M(x*x)));
const inv=F.map(x=>F.find(y=>M(x*y)===1)??0);
const leg=F.map(x=>x===0?0:sq.has(x)?1:-1);
const tr=a=>{while(a.length&&a.at(-1)===0)a.pop();return a;};
const rem=(a,b)=>{
 a=a.slice();while(a.length>=b.length){
   const f=M(a.at(-1)*inv[b.at(-1)]),n=a.length-b.length;
   for(let j=0;j<b.length;j++)a[n+j]=M(a[n+j]-f*b[j]);
   tr(a);
 } return a;
};
const sf=a=>{let b=tr(a.slice(1).map((v,i)=>M(v*(i+1))));while(b.length)[a,b]=[b,rem(a,b)];return a.length===1;};
const rows=[];
for(const d of [5,6]){
 let a=Array(d+1).fill(0);a[d]=1;
 let z={d,total:0,rawHigh:0,squarefreeHigh:0,rawAt32:0,squarefreeAt32:0,control:null};
 function dfs(i){
   if(i<d-1){for(const x of F){a[i]=x;dfs(i+1)}return;}
   z.total++;
   // For degree five, infinity contributes one branch point. For degree six,
   // infinity is split and the current branch intentionally allows no roots.
   let penalty=d===5?1:0;
   for(const x of F){
     let y=1;for(let k=d-1;k>=0;k--)y=M(y*x+a[k]);
     if(d===6&&y===0)return;
     penalty+=1-leg[y];
     if(penalty>4)return;
   }
   if(penalty<=3){
     z.rawHigh++;if(sf(a))z.squarefreeHigh++;
   }else{
     z.rawAt32++;
     if(sf(a)){z.squarefreeAt32++;z.control??=a.slice();}
   }
 }
 dfs(0);if(z.squarefreeHigh)throw Error('genus two 33+ control');rows.push(z);
}
console.log(JSON.stringify(rows));
