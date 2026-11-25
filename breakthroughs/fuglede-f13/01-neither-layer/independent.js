"use strict";
// Independent fixed-abscissa planar search and projective K10 atlas, F13.
const p=13, F=Array.from({length:p},(_,i)=>i), mod=x=>(x%p+p)%p;
const inv=F.map(x=>F.find(y=>mod(x*y)===1)||0);
const assert=(x,m)=>{if(!x)throw Error(m);};
const vec=(N,f)=>Array.from({length:N},(_,i)=>f(i));
const sum=x=>x.reduce((s,v)=>s+v,0);
const pop=n=>{let c=0;while(n){n&=n-1;c++;}return c;};
const subsets=(n,k)=>vec(1<<n,i=>i).filter(m=>pop(m)===k);
const unpack=m=>F.filter(i=>m>>i&1);
function canon(x){
  let m=1<<13;
  for(let a=1;a<13;a++)for(let b=0;b<13;b++)
    m=Math.min(m,x.reduce((m,t)=>m|(1<<mod(a*t+b)),0));
  return m;
}
function planar(n,bound,cap) {
  const C=new Set(subsets(13,n).map(m=>canon(unpack(m))));
  const good=[], nodes=Array(n+1).fill(0);
  for(const m of [...C].sort((a,b)=>a-b)){
    const x=unpack(m), d=x.map(a=>x.map(b=>inv[mod(a-b)])), h=[0,0];
    let firstLeaf=false;
    function go(scaled,mask){
      const i=h.length;nodes[i]++;
      if(i===n){assert(scaled,"all collinear");firstLeaf=true;return;}
      for(const v of (scaled?F:[0,1])){
        let dirs=mask,ok=true;
        const geometries=Array(13).fill(0);
        for(let j=0;j<i;j++){
          const q=mod((v-h[j])*d[i][j]);dirs|=1<<q;
          if(pop(dirs)>bound||++geometries[q]>=cap){ok=false;break;}
        }
        if(ok){h.push(v);go(scaled||v!==0,dirs);h.pop();}
        if(firstLeaf)return; // per-abscissa feasibility only, not a model count
      }
    }
    go(false,1);
    if(firstLeaf)good.push(m);
  }
  return {n,bound,cap,supportClasses:C.size,good,nodes:nodes.slice(2)};
}
const six=planar(6,5,5), golden=[111,187,215,311,427];
assert(six.good.join(',')===golden.join(','),"six-shadow gate");
const negatives=[[6,4,5],[10,7,5],[11,9,5],[12,9,5]].map(x=>planar(...x));
assert(negatives.every(q=>q.good.length===0),"planar lower bound");
console.log(JSON.stringify({six,negatives}));

function shade(L,L0) {
  // Each projective bad slope is a linear form on the Boolean grid.
  return L0.map(([x,y])=>mod(L[0]*x+L[1]*y));
}
const allForms=[[1,0],...F.map(t=>[mod(-t),1])];
const pts=F.flatMap(x=>F.map(y=>[x,y]));
const arg=allForms.map(L=>shade(L,pts));
function shapesOf(f){
  const profiles=arg.map(y=>{
    const rows=vec(13,()=>[]);
    for(let j=0;j<pts.length;j++)if(f[j]){
      const t=y[j];
      rows[t].push(y===arg[0]?pts[j][1]:pts[j][0]);
    }
    return {peak:Math.max(...rows.map(r=>r.length)),
      excluded:rows.filter(r=>r.length===6&&!golden.includes(canon(r))).length};
  });
  return profiles;
}
function nullspace(A){
  const R=A.map(r=>r.map(mod)), piv=[];
  for(let j=0;j<R[0].length;j++){
    const z=piv.length,k=R.findIndex((x,i)=>i>=z&&x[j]!==0);
    if(k<0)continue;
    [R[z],R[k]]=[R[k],R[z]];
    R[z]=R[z].map(x=>mod(x*inv[R[z][j]]));
    for(let i=0;i<R.length;i++)if(i!==z){
      const a=R[i][j];R[i]=R[i].map((x,k)=>mod(x-a*R[z][k]));
    }
    piv.push(j);
  }
  return vec(R[0].length,i=>i).flatMap(j=>piv.includes(j)?[]:[
    vec(R[0].length,i=>i===j?1:piv.includes(i)?mod(-R[piv.indexOf(i)][j]):0)
  ]);
}
function rows(a,b,c,cap=13){
  const g=F.map(t=>mod(a*t*t+b*t+c));
  if(g.some(x=>x>cap))return [];
  const k=(78-sum(g))/13, z=g.flatMap((x,i)=>x===0?[i]:[]);
  if(k<0||!Number.isInteger(k)||k>z.length||(k&&cap<13))return [];
  return subsets(z.length,k).map(m=>g.map((x,i)=>x!==0?x:(m>>z.indexOf(i)&1)?13:0));
}
function scale(x,a){return x.map(v=>mod(a*v));}
const lambdas=[];
for(let t=2;t<=12;t++){
  const q=mod(1-t);
  if(t===Math.min(t,q,inv[t],inv[q],mod(t*inv[mod(t-1)]),mod((t-1)*inv[t])))lambdas.push(t);
}
const atlas=[];
for(const cap of [11,13]){
  for(const lam of lambdas){
    const ids=[0,1,2,1+lam], L=ids.map(i=>allForms[i]);
    const qker=nullspace([L.map(([a,b])=>a*a),L.map(([a,b])=>2*a*b),L.map(([a,b])=>b*b)]);
    const lker=nullspace([L.map(([a])=>a),L.map(([a,b])=>b)]);
    assert(qker.length===1&&qker[0].every(x=>x!==0)&&lker.length===2,"atlas kernel");
    const stats={cap,lam,q2:0,q1:0,boolean:0,minSix:13,peakTypes:{}};
    function evaluate(rr,kind){
      const grid=pts.map((_,j)=>sum(rr.map((r,i)=>r[arg[ids[i]][j]]))-18);
      if(grid.some(x=>![0,13,26].includes(x)))return;
      const f=grid.map(x=>x/13);
      for(let l=0;l<14;l++){
        const real=Array(13).fill(0);
        f.forEach((v,j)=>real[arg[l][j]]+=v);
        assert(real.every((x,i)=>x===(ids.includes(l)?rr[ids.indexOf(l)][i]:6)),
               "independent consistency checks");
      }
      stats[kind]++;
      if(f.some(x=>x>1))return;
      stats.boolean++;
      const sh=shapesOf(f), k=sh.map(x=>x.peak).sort((a,b)=>a-b).join(',');
      stats.peakTypes[k]=(stats.peakTypes[k]||0)+1;
      assert(k==="6,6,6,6,6,6,6,6,6,6,10,10,10,12","peak type");
      for(const r of sh.filter(r=>r.peak===6)){assert(r.excluded>0,"missing shape");stats.minSix=Math.min(stats.minSix,r.excluded);}
    }
    function constants(a,b,gauge,kind){
      const c=[0,0,0,0], to=vec(4,i=>i).filter(i=>!gauge.includes(i));
      function fill(j){
        if(j===to.length){
          if(mod(sum(c))!==5)return;
          const r=vec(4,i=>rows(a[i],b[i],c[i],cap));
          for(const A of r[0])for(const B of r[1])for(const C of r[2])for(const D of r[3])evaluate([A,B,C,D],kind);
          return;
        }
        const i=to[j];
        for(const k of (a[i]||b[i]?F:[0])){c[i]=k;fill(j+1);}
      }
      fill(0);
    }
    for(let s=1;s<13;s++)constants(scale(qker[0],s),[0,0,0,0],[],"q2");
    for(let s=0;s<13;s++)for(let t=0;t<13;t++)if(s||t){
      const b=vec(4,i=>mod(s*lker[0][i]+t*lker[1][i]));
      constants([0,0,0,0],b,vec(4,i=>i).filter(i=>b[i]!==0).slice(0,2),"q1");
    }
    atlas.push(stats);
  }
}
assert(atlas.filter(x=>x.cap===11).every(x=>x.q1+x.q2===0),"cap eleven exclusion");
assert(sum(atlas.map(x=>x.q2))===804&&sum(atlas.map(x=>x.q1))===0&&sum(atlas.map(x=>x.boolean))===24,"atlas totals");
const tri=shapesOf(pts.map(([r,s])=>+(r+s>=13)));
assert(tri.filter(q=>q.peak===6&&q.excluded>=1).length===11&&tri.filter(q=>q.peak===12).length===3,"triangle directions");
console.log(JSON.stringify({lambdas,atlas,tri}));

// Second independent census: evaluate every cubic by its four sample values.
// Explicitly enumerate every admissible choice of lifted zeros, rather than
// using only a binomial multiplicity plus a variance shift.
const uni=F.slice(0,4), A=F.slice(4).map(x=>uni.map(y=>{
 let a=1;for(const w of uni)if(w!==y)a=mod(a*(x-w)*inv[mod(y-w)]);return a;
}));
const hist={};let layered=0,profiles=0,visits=0;
for(const a of F)for(const b of F)for(const c of F)for(const d of F){
  visits++;
  const h=[a,b,c,d], g=h.concat(A.map(row=>mod(sum(row.map((r,i)=>r*h[i])))));
  const z=g.flatMap((x,i)=>!x?[i]:[]), excess=(78-sum(g))/13;
  if(!Number.isInteger(excess)||excess<0||excess>z.length)continue;
  for(const m of subsets(z.length,excess)){
    const v=g.slice();z.forEach((i,j)=>{if(m>>j&1)v[i]=13;});
    assert(sum(v)===78&&Math.max(...v)<=13,"profile cap/mass");
    if(v.every(x=>x===13||x===0)){layered++;continue;}
    const k=(sum(v.map(x=>x*(x-1)/2))-195)/13;
    assert(k>=0&&k<=12&&Number.isInteger(k),"weight cap");
    profiles++;hist[k]=(hist[k]||0)+1;
  }
}
assert(visits===13**4&&profiles===23765&&layered===1716&&hist[12]===208,"cubic census");
console.log(JSON.stringify({visits,profiles,layered,hist}));
