"use strict";
function ok(x){if(!x)throw Error("certificate");}
function* projective(p,n){
 for(let j=0;j<n;j++){let z=Array(n).fill(0);z[j]=1;
  function* go(k){if(k===n){yield z.slice();return;}for(let a=0;a<p;a++){z[k]=a;yield* go(k+1);}}
  yield* go(j+1);
 }
}
function sing(p,k){
 const alpha=Array.from({length:4},(_,i)=>+(i<k));let n=0;
 for(const z of projective(p,4)){
  const mult=S=>S.reduce((v,i)=>v*z[i]%p,1);
  const ids=[0,1,2,3];
  const value=alpha.reduce((v,a,i)=>v+a*mult(ids.filter(q=>q!==i)),0)%p;
  const d=ids.map(h=>alpha.reduce((v,a,i)=>h===i?v:v+a*mult(ids.filter(q=>q!==i&&q!==h)),0)%p);
  const s=!value&&d.every(x=>!x);
  const f=z.filter(Boolean).length;
  if(k===4)ok(s===(f===1));
  if(k===3){const q=(z[0]*z[1]+z[0]*z[2]+z[1]*z[2])%p;ok(s===(z[3]===0? q===0 : f===1));}
  if(k<=2){const u=z[2],v=z[3],w=(z[1]+(k===2?z[0]:0))%p;ok(s===([u,v,w].filter(x=>!x).length>=2));}
  n+=+s;
 }
 ok(n===(k<=2?3*p+1:k===3?p+2:4));
 return n;
}
function sharp(p,c){
 const rat=Math.floor(c/2),E=[];
 if(c%2)E.push([0,0,0]);
 for(let j=0;j<3;j++)for(let i=1;i<=rat;i++){const a=[0,0,0];a[j]=i;E.push(a);}
 const P=[...projective(p,3)];let cap=0;
 for(const v of P){const a=Array(p).fill(0);for(const x of E)a[x.reduce((u,t,i)=>u+t*v[i],0)%p]++;cap=Math.max(cap,...a);}
 for(const x of E)for(const y of E)ok(x.some((a,j)=>a===y[j]));
 ok(cap===c&&E.length===Math.floor(3*c/2));
 return [c,E.length];
}
const o={p5:[1,2,3,4].map(k=>sing(5,k)),p17:[1,2,3,4].map(k=>sing(17,k)),sharp:Array.from({length:14},(_,i)=>sharp(17,i+4))};
console.log(JSON.stringify(o));
