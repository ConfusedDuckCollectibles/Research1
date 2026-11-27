'use strict';
// Binary moment and affine secant gates for the quartic-plus-twisted-cubic branch.
const p=17,F=[...Array(p).keys()],M=x=>(x%p+p)%p,S=[...F,17];
const inv=F.map(x=>F.find(y=>M(x*y)===1)??0);
const pow=(x,k)=>{let y=1;while(k--)y=M(x*y);return y};
const norm=a=>{const x=inv[a.find(c=>c!==0)];return a.map(c=>M(c*x));};
const key=a=>a.join(',');
const pd=[];
for(const a of F)for(const b of F)for(const c of F)for(const d of F)
  if([a,b,c,d].find(x=>x!==0)===1)pd.push([a,b,c,d]);
const gens=([a,b,c,d])=>[norm([M(a+b+c+d),M(b+2*c+3*d),M(c+3*d),d]),
  norm([a,M(3*b),M(9*c),M(10*d)]),norm([d,c,b,a])];
const seen=new Set(),out=[];
const mon=(d,i,t)=>t===17?Number(i===d):pow(t,i);
const rank=(A)=>{
 A=A.map(x=>x.map(M));let row=0;
 for(let j=0;j<A[0].length;j++){
  let z=A.findIndex((a,i)=>i>=row&&a[j]!==0);if(z<0)continue;
  [A[row],A[z]]=[A[z],A[row]];const u=inv[A[row][j]];
  for(let k=j;k<A[0].length;k++)A[row][k]=M(A[row][k]*u);
  for(let i=0;i<A.length;i++)if(i!==row&&A[i][j]!==0){
   const u=A[i][j];for(let k=j;k<A[0].length;k++)A[i][k]=M(A[i][k]-u*A[row][k]);}
  row++;
 }return{rank:row,R:A};
};
for(const H of pd){const k=key(H);if(seen.has(k))continue;
 const q=[H];seen.add(k);for(let i=0;i<q.length;i++)for(const h of gens(q[i])){let k=key(h);if(!seen.has(k)){seen.add(k);q.push(h)}}
 const ht=t=>M(H.reduce((s,a,i)=>s+a*mon(3,i,t),0));
 const dom=S.filter(t=>ht(t)!==0);
 const w=dom.map(t=>Array.from({length:16},(_,i)=>M(mon(15,i,t)*pow(inv[ht(t)],5))));
 const diag=Array(16).fill(0),taken=[],ans=[];
 for(const n of [12,13]){
  const st={n,checked:0,ranks:{},nonzero:0,first:null};
  function visit(j,k){
   if(k===0){
    st.checked++; let {rank:r,R}=rank(Array.from({length:11},(_,i)=>diag.slice(i,i+6)));
    st.ranks[r]=(st.ranks[r]??0)+1;
    if(r<6){
      if(r<17-n)throw Error('unexpected large nullity');
      const lead=R.slice(0,r).map(a=>a.findIndex(x=>x!==0)),free=F.slice(0,6).filter(x=>!lead.includes(x));
      const coeffs=free.length===1?[[1]]:F.map(t=>[1,t]).concat([[0,1]]);
      for(const coeff of coeffs){
        const g=Array(6).fill(0);free.forEach((j,i)=>g[j]=coeff[i]);
        lead.forEach((j,i)=>g[j]=M(-free.reduce((s,k,l)=>s+R[i][k]*coeff[l],0)));
        if(taken.every(i=>M(g.reduce((s,a,j)=>s+a*mon(5,j,dom[i]),0))!==0)){
         for(let l=0;l<=10;l++){
          const sum=taken.reduce((s,i)=>s+M(g.reduce((a,b,j)=>a+b*mon(5,j,dom[i]),0))*
             mon(10,l,dom[i])*pow(inv[ht(dom[i])],5),0);
          if(M(sum)!==0)throw Error('direct witness check');
         }
         st.nonzero++;st.first??={g,points:taken.map(i=>dom[i])};
        }
      }
    }
    return;
   }
   for(let i=j;i+k<=dom.length;i++){
    for(let a=0;a<16;a++)diag[a]+=w[i][a];taken.push(i);
    visit(i+1,k-1);taken.pop();for(let a=0;a<16;a++)diag[a]-=w[i][a];
   }
  }
  visit(0,n);ans.push(st);
 }
 out.push({H,orbit:q.length,dom:dom.length,ans});
}
if(seen.size!==5220)throw Error('incomplete forms');
if(out.length!==5||out.slice(1).some(o=>o.ans.some(z=>z.nonzero||Object.keys(z.ranks).join()!=="6")))
  throw Error('non-triple infinity gate');
const pg=[];
for(let a of F)for(let b of F)for(let c of F)if([a,b,c].find(x=>x!==0)===1)pg.push([a,b,c]);
const lines=F.map(u=>F.map(v=>{
 const s=M(u+v),t=M(u*u+u*v+v*v);
 return pg.map(([a,b,c],i)=>M(a+b*s+c*t)===0?i:-1).filter(i=>i>=0);
}));
const hist={};
for(const n of [12,13]){
 const h={},freq=new Uint8Array(307),taken=[];
 function visit(j,k){
   if(!k){const max=Math.max(...freq);h[max]=(h[max]??0)+1;return;}
   for(let i=j;i+k<=p;i++){
     for(const u of taken)for(const a of lines[u][i])freq[a]++;
     taken.push(i);visit(i+1,k-1);taken.pop();
     for(const u of taken)for(const a of lines[u][i])freq[a]--;
   }
 }
 visit(0,n);
 if(Object.keys(h).some(x=>Number(x)<7))throw Error('small projective pencil load');
 hist[n]=h;
}
console.log(JSON.stringify({infinityForms:out,secantPencilMaxima:hist}));
