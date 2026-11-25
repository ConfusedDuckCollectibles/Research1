// Exact local certificates for the remaining size-52 case; Node.js, no packages.
const p=13, mod=x=>((x%p)+p)%p;
function check(ok,problem){if(!ok)throw new Error(problem);}
const hist={},peak={};
for(let n=0;n<p**5;n++){
  let k=n,co=[];
  for(let j=0;j<5;j++){co.push(k%p);k=Math.floor(k/p);}
  const g=Array.from({length:p},(_,x)=>{let c=0;for(let j=4;j>=0;j--)c=mod(c*x+co[j]);return c;});
  if(g.some(z=>z>10)||g.reduce((a,b)=>a+b,0)!==52) continue;
  const pair=g.reduce((a,b)=>a+b*(b-1)/2,0), w=(pair-78)/p;
  check(Number.isInteger(w)&&w>=0&&w!==1,'profile collision-weight minimum');
  check(w!==0||g.every(x=>x===4),'pointwise uniformity at weight zero');
  hist[w]=(hist[w]??0)+1;
  const m=Math.max(...g);peak[m]=(peak[m]??0)+1;
}
check(hist[0]===1&&hist[2]>0,'profile zero/sharp controls');
const byU=[];
for(let m=66;m<=73;m++){
  const eight=Math.floor((m*(m-1)/2-56*m+1830)/6);
  const available=18*(73-m),needed=2*(110-eight);
  check(needed>available,'direction-incidence deficit at '+m);
  byU.push({m,eight,available,needed});
}
// An independent finite alternative to the application of Cauchy–Davenport/Vosper.
const sets=[], ap=new Map();
const mask=S=>S.reduce((m,x)=>m|(1<<x),0);
for(let x=0;x<p;x++)for(let d=1;d<p;d++){
  const z=mask(Array.from({length:4},(_,i)=>mod(x+i*d))),c=Math.min(d,p-d);
  check(!ap.has(z)||ap.get(z)===c,'arithmetic progression step uniqueness');
  ap.set(z,c);
}
for(let a=0;a<p;a++)for(let b=a+1;b<p;b++)for(let c=b+1;c<p;c++)for(let d=c+1;d<p;d++)sets.push([a,b,c,d]);
function pop(x){let n=0;for(;x;x&=x-1)n++;return n;}
let equal=0;
for(let i=0;i<sets.length;i++)for(let j=i;j<sets.length;j++){
  let differences=0;
  for(const a of sets[i])for(const b of sets[j])differences|=1<<mod(a-b);
  const n=pop(differences);
  check(n>=7,'two four-sets must have at least seven differences');
  if(n===7){
    equal++;
    const d=ap.get(mask(sets[i]));
    check(d!==undefined&&d===ap.get(mask(sets[j])),'exact seven differences force common progression step');
  }
}
check(sets.length===715&&equal>0,'pair-exhaustion control');
console.log(JSON.stringify({profilesByCollisionWeight:hist,profilesByPeak:peak,byU,localSumset:{sets:sets.length,ap4:ap.size,pairs:sets.length*(sets.length+1)/2,exactlySeven:equal}}));
