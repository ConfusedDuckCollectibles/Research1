'use strict';
const p=17,n=Number(process.argv[2]??7),target=p*n*(p-n), mod=x=>(x%p+p)%p;
const pool=Array.from({length:17},()=>[]);
for(let a=1;a<p;a++)for(let c=0;c<p;c++){
 const vals=Array.from({length:p},(_,x)=>mod(a*x*x+c)),sum=vals.reduce((s,x)=>s+x,0);
 if(sum>p*n)continue;
 const z=vals.flatMap((v,i)=>v===0?[i]:[]);
 const k=(p*n-sum)/p;if(!Number.isInteger(k)||k>z.length)continue;
 for(let m=0;m<1<<z.length;m++){
  let take=z.filter((_,i)=>(m>>i)&1);if(take.length!==k)continue;
  let v=vals.slice();for(let i of take)v[i]=p;
  const ss=v.reduce((s,x)=>s+(x-n)**2,0);
  pool[a].push({v,ss,c,take});
 }
}
let variance=0,constant=0,boolean=0;
for(let r=2;r<p;r++)for(let q=1;q<p;q++){
 const a=mod(r*(1-r)*q), b=mod((r-1)*q), c=mod(-r*q);
 for(let A of pool[a])for(let B of pool[b])for(let C of pool[c])for(let D of pool[q]){
   if(A.ss+B.ss+C.ss+D.ss!==target)continue;
   variance++;
   if(mod(A.c+B.c+C.c+D.c-3*n)!==0)continue;
   constant++;
   let count=0, bad=false;
   for(let x=0;x<p&&!bad;x++)for(let y=0;y<p;y++){
    let u=A.v[x]+B.v[y]+C.v[mod(y-x)]+D.v[mod(y-r*x)]-3*n;
    if(u===p)count++;else if(u!==0){bad=true;break;}
   }
   if(!bad&&count===p*n)boolean++;
 }
}
console.log(JSON.stringify({counts:pool.map(p=>p.length),variance,constant,boolean}));
