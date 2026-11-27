'use strict';
// Only integer arithmetic on F17. This checks local geometric and binary-form
// ingredients; the global deficiency and blocking arguments are in C-CURVES.md, section C2.
const p=17, r=n=>(n%p+p)%p, range=n=>[...Array(n).keys()], F=range(p);
const expect=(b,s)=>{if(!b)throw Error(s);};
const inv=a=>{expect(r(a)!==0,"inverse of 0");return F.find(b=>r(a*b)===1);};
const dot=(a,b)=>r(a.reduce((s,x,i)=>s+x*b[i],0));
const scale=(a,s)=>a.map(x=>r(x*s));
const canon=a=>scale(a,inv(a.find(x=>x!==0)??0));
const key=a=>a.join(',');
const squares=new Set(F.slice(1).map(x=>r(x*x)));
const type=x=>x===0?'zero':squares.has(x)?'square':'nonsquare';
const pg=[];for(const a of F)for(const b of F)for(const c of F)
  if([a,b,c].find(x=>x!==0)===1)pg.push([a,b,c]);
expect(pg.length===307,'projective points');
const all=[];for(const x of F)for(const y of F)for(const z of F)all.push([x,y,z]);
const q=([x,y,z])=>r(x*y+z*z);
const dual=([a,b,c])=>r(4*a*b+c*c);
const bmap=([x,y,z])=>[r(9*y),r(9*x),z]; // Bx
const Q=c=>all.filter(v=>q(v)===c);
const occ=(pts,n)=>{const a=Array(p).fill(0);for(const v of pts)a[dot(v,n)]++;return a;};
const fingerprint=a=>a.slice().sort((x,y)=>x-y).join(',');
const binary=[];
for(const b of F)for(const c of F){
  const val=[0,1,2,9,16].map(s=>r(1+b*s+c*s*s)).concat(c);
  if(val.every(x=>squares.has(x)))binary.push([b,c]);
}
expect(binary.length===12&&binary.every(([b,c])=>r(b*b-4*c)===0),'G6 sign gate');

const profiles={};
for(const radius of [1,3]){
  const out={zero:0,square:0,nonsquare:0,patterns:{}}, pts=Q(radius);
  for(const n of pg){
    const t=type(dual(n)),a=occ(pts,n),count=z=>a.filter(x=>x===z).length;
    out[t]++; const k=t+":"+fingerprint(a);out.patterns[k]=(out.patterns[k]??0)+1;
    if(radius===3){
      expect(t==='square' ? count(16)===17 : t==='zero' ? a[0]===0&&count(17)===16 :
        count(1)===2&&count(18)===15,'elliptic level types');
    }else{
      expect(t==='square' ? count(33)===2&&count(16)===15 : t==='zero' ?
        a[0]===34&&count(17)===16 : count(18)===17,'split level types');
    }
  }
  profiles[radius]={size:pts.length,classes:[out.zero,out.square,out.nonsquare],
    patterns:Object.values(out.patterns).sort((a,b)=>a-b)};
}

// Pair indices are extracted from the matrix image and kernel, not from a
// ruling occupancy formula. Check the tangent-plane point sets individually.
const dirId=([a,b])=>a===0?17:r(b*inv(a));
const ij=([x,y,z],t)=>{
  const a=r(t+z),b=x,c=y,d=r(t-z);
  const col=a||c?[a,c]:[b,d],row=a||b?[a,b]:[c,d];
  return [dirId(col),dirId([r(-row[1]),row[0]])];
};
const split=Q(1).map(v=>({v,ij:ij(v,1)}));
expect(new Set(split.map(x=>key(x.ij))).size===306&&split.every(x=>x.ij[0]!==x.ij[1]),'grid');
for(const {v,ij:ab} of split){
  expect(key(ij(scale(v,-1),1))===key(ab.slice().reverse()),'antipode');
  const n=bmap(v);
  for(const x of split)expect((dot(n,x.v)===1)===(x.ij[0]===ab[0]||x.ij[1]===ab[1]),'finite tangent');
}
const horizon=pg.filter(v=>q(v)===0),tangents=pg.filter(n=>dual(n)===0).map(n=>[n,[r(2*n[1]),r(2*n[0]),n[2]]]);
expect(horizon.length===18&&tangents.length===18,'horizon');
for(const v of horizon){
  const [i,j]=ij(v,0);expect(i===j,'diagonal');
  for(const x of split)expect((dot(bmap(v),x.v)===0)===(x.ij[0]===i||x.ij[1]===i),'infinite tangent');
}
const cover={zero:{},square:{},nonsquare:{}};
for(const n of pg){
  const t=type(dual(n)),k=tangents.filter(([u,pol])=>dot(n,pol)===0).length;
  expect(k===(t==='zero'?1:t==='square'?2:0),'dual conic tangent cover');
  cover[t][k]=(cover[t][k]??0)+1;
}

const controls=[];
for(let fixture=0;fixture<3;fixture++){
  const allowed=(i,j)=>fixture===0 ? i<17&&j<17&&[1,16].includes(r(j-i)) :
    [1,17].includes((j-i+18)%18)&&!(i===0&&j===1)&&!(fixture===1?i===1&&j===0:i===3&&j===4);
  const pts=split.filter(x=>allowed(...x.ij));expect(pts.length===34,'fixture size');
  const u=range(18).map(i=>2-pts.filter(x=>x.ij[0]===i).length);
  const v=range(18).map(i=>2-pts.filter(x=>x.ij[1]===i).length);
  expect(u.every(x=>x>=0)&&v.every(x=>x>=0),'fixture cap on rulings');
  const T=new Set(), rule=new Set();
  for(const {v:x,ij:[i,j]} of split)if(u[i]+v[j]+Number(allowed(i,j))===2&&
      u[j]+v[i]+Number(allowed(j,i))===2)T.add(key(canon(bmap(x))));
  for(const h of horizon){const [i]=ij(h,0);if(u[i]+v[i]===2)T.add(key(canon(bmap(h))));}
  for(const n of pg)if(type(dual(n))!=='nonsquare'&&occ(pts.map(x=>x.v),n).every(x=>x===2)){
    rule.add(key(n));expect(T.has(key(n)),'actual uniform must meet tangent rule');
  }
  if(fixture===0){
    expect(T.size===17,'one wholly empty row/col potential');
    const [j]=range(18).filter(i=>u[i]+v[i]===4);
    const atInfinity=bmap(horizon.find(x=>ij(x,0)[0]===j));
    const some=bmap(split.find(x=>x.ij[0]===j).v);
    // Any projective line containing two independent points is the unique
    // kernel of a nonzero functional. Find it directly and verify all 17.
    const normal=pg.find(n=>dot(atInfinity,n)===0&&dot(some,n)===0);
    expect([...T].every(s=>dot(s.split(',').map(Number),normal)===0),'one pencil');
  }else expect(T.size<=6,'multiple-defect bound');
  controls.push({deficits:u.map((a,i)=>a+v[i]).filter(Boolean).sort(),
    tangentCandidates:T.size,actuallyUniform:rule.size});
}
// Templates before translation: N=1 (cylinder), N=z (paraboloid), N=z² (cone).
const N=([x,y])=>r(x*x-x*y+y*y);
const cyl=all.filter(v=>N(v)===1), par=all.filter(v=>N(v)===v[2]);
const cone=all.filter(v=>N(v)===r(v[2]*v[2]));
const horiz=pg.filter(x=>x[2]===0);
expect(cyl.length===306&&par.length===289&&cone.length===289,'other quadrics');
for(const n of horiz)expect(occ(cyl,n).filter(Boolean).length<=10,'cylinder horizontal gap');
for(const n of pg.filter(x=>x[2]!==0))expect(occ(par,n).filter(x=>x===1).length===1,'paraboloid singleton');
expect(cone.filter(v=>v[2]===1).length===18,'cone generators');
console.log(JSON.stringify({binary:{tested:289,passing:binary.length,coefficientPairs:binary},profiles,
  checked:{orderedAffineSplitPoints:split.length,finiteTangentChecks:split.length**2,
  infiniteTangentChecks:horizon.length*split.length,dualityChecks:pg.length*tangents.length},
  cover,controls,other:{cylinder:cyl.length,paraboloid:par.length,cone:cone.length}}));
