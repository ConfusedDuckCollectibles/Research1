
const cmp=(a,b)=>a-b;
function rows(a,k) {
  let r=[{d:0,p:0,q:0,c:false}];
  for(let i=0;i<a.length;i++)for(let j=i+1;j<a.length;j++)r.push({d:a[j]-a[i],p:i,q:j,c:i<=k&&k<j});
  return r.sort((x,y)=>x.d-y.d);
}
const gapAt=(r,a)=>a[r.q]-a[r.p];
function preserves(r,a,unit=1) {
  return r.every((u,i)=>i===0||gapAt(u,a)-gapAt(r[i-1],a)>=unit);
}
function diff(v,n) {return (v%n+n)%n}
function pds(a,n) {const d=new Set(a.flatMap(x=>a.filter(y=>x!==y).map(y=>diff(x-y,n))));return d.size===n-1&&a.length*(a.length-1)===n-1}
function floorBridge(a,k){
  const r=rows(a,k);
  let b=Infinity;
  for(let s=0;s<r.length;s++)if(!r[s].c)for(let t=s+1;t<r.length;t++)if(r[t].c)b=Math.min(b,r[t].d-r[s].d-1);
  const witness=a.map((x,i)=>x-(i>k?b:0));
  if(!preserves(r,witness))throw Error("facet violation");
  return [b,witness,r];
}
function certificate(r,y,dual) {
  const n=y.length-1;
  let sum=Array(n).fill(0),rankRows=[],weights=0,blocked=Array(n).fill(false);
  for(let i=1;i<r.length;i++){
    const u=r[i], v=r[i-1], w=dual[u.d]||0;
    let col=sum.map((x,j)=>+(u.p<=j&&j<u.q)-+(v.p<=j&&j<v.q));
    if(u.d-v.d===1)for(let j=0;j<n;j++)blocked[j] ||= col[j]>0;
    if(w) {
      if(gapAt(u,y)-gapAt(v,y)!==1)throw Error("not tight");
      sum=sum.map((x,j)=>x+w*col[j]);
      weights+=w;rankRows.push(col);
    }
  }
  function det(m){return m.length===0?1:m[0].reduce((v,a,j)=>v+(j%2?-1:1)*a*det(m.slice(1).map(r=>r.filter((_,i)=>i!==j))),0)}
  return {sum,weights,rankRows,det:det(rankRows),blocked};
}
function second(a,v) {
  const gaps=a.map((x,i)=>diff(a[(i+1)%a.length]-x,v));
  const g=[...gaps].sort(cmp).at(-2), m=Math.max(...gaps);
  const t=a[(gaps.indexOf(g)+1)%a.length];
  const norm=a.map(x=>diff(x-t,v)).sort(cmp);
  const k=norm.findIndex((p,i)=>i+1<norm.length&&norm[i+1]-p===m);
  const [b,opt,r]=floorBridge(norm,k);
  if (g+b>m)throw Error("beat largest cut");
  return {norm,b,opt,r,k};
}
function enumPDS(v,n){
 let visit=0, map=new Map;
 function rec(a,D){
   if(a.length===n){
     visit++;
     const sub=second(a,v);
     map.set(sub.norm.join(),sub);
     return;
   }
   for(let x=a.at(-1)+1;x<=v-(n-a.length);x++){
     let dd=D,good=true;
     for(const u of a){
       const i=x-u,j=v-i;
       const mask=(1n<<BigInt(i))|(1n<<BigInt(j));
       if(dd&mask){good=false;break}
       dd|=mask;
     }
     if(good)rec([...a,x],dd);
   }
 }
 rec([0],0n);
 return {pdsAnchored:visit,uniqueSecondCuts:map.size,shortened:[...map.values()].filter(s=>s.b>0).length,certifiedUnique:[...map.values()].filter(s=>forcing(s.r,s.opt,s.k)).length};
}

// Same-length perturbations eliminate the exceptional coordinate.
// Explicit signed vectors are independent from the primary bitmask encoding.
function forcing(r,y,k) {
  const m=y.length-1;
  const js=Array.from({length:m},(_,j)=>j).filter(j=>j!==k);
  const coeff=(d,j)=>+(d.p<=j&&j<d.q)-+(d.p<=k&&k<d.q);
  const tight=r.slice(1).flatMap((u,i)=>gapAt(u,y)-gapAt(r[i],y)===1?
    [js.map(j=>coeff(u,j)-coeff(r[i],j))]:[]);
  const bad=new Set;
  for (;;) {
    const size=bad.size;
    for (const row of tight)
      if (row.every((a,j)=>a<=0||bad.has(j)))
        row.forEach((a,j)=>{if(a<0)bad.add(j)});
    if (size===bad.size) return bad.size===js.length;
  }
}
const a0=[0,1,11,19,26,28],a1=[0,1,8,13,17,19],dual={1:7,2:3,9:2,10:4,15:3};
const r=rows(a0,-1);
const x={q2:second([0,4,5],7),q17:second([0,18,23,37,58,62,75,86,177,179,187,199,230,233,259,260,266,275],307),
 multiGap:{pds:pds(a0,31),sameChamber:preserves(r,a1),cert:certificate(r,a1,dual)},
 enumerated:[7,13,21,31].map((v,j)=>enumPDS(v,3+j))};

delete x.q2.r; delete x.q17.r;
x.multiGap.largestCutLength=31-Math.max(...a0.map((a,i)=>diff(a0[(i+1)%a0.length]-a,31)));
const B=[0,9,12,22,23,27,29], Btwice=[0,17,23,42,44,52,56];
const rb=rows(B,-1);
x.BosePrior={preserves:preserves(rb,Btwice,2),blocked:certificate(rb,B,{}).blocked};
if (x.q2.b!==2 || x.q17.b!==1 || !x.multiGap.pds || !x.multiGap.sameChamber ||
    x.multiGap.cert.det!==-1 || x.multiGap.cert.weights!==19 ||
    x.multiGap.largestCutLength!==21 || !x.multiGap.cert.sum.every(a=>a===1) ||
    !x.multiGap.cert.blocked.every(Boolean) || !x.BosePrior.preserves || !x.BosePrior.blocked.every(Boolean) ||
    JSON.stringify(x.enumerated.map(a=>[a.pdsAnchored,a.uniqueSecondCuts,a.shortened,a.certifiedUnique]))!==
    JSON.stringify([[6,2,2,2],[16,4,4,4],[10,2,2,2],[60,10,6,10]]))
  throw Error("unexpected independent audit replay");

console.log(JSON.stringify(x,null,2));