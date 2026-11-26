
const pop=(x)=>{ let n=0; while(x){x&=x-1;n++;}return n;};
function triples(n){let a=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let k=j+1;k<n;k++)a.push((1<<i)|(1<<j)|(1<<k));return a;}
function mkref(){
 const ans=[7,56,448];
 for(let g=0;g<3;g++) for(let i=3*g;i<3*g+3;i++) for(let j=i+1;j<3*g+3;j++) ans.push((1<<9)|(1<<i)|(1<<j));
 let l=10;
 for(let g=0;g<3;g++)for(let h=g+1;h<3;h++,l++)for(let i=3*g;i<3*g+3;i++)for(let j=3*h;j<3*h+3;j++)ans.push((1<<l)|(1<<i)|(1<<j));
 return ans;
}
function sunflower(a,b,c,d){const x=a&b;return (a&c)===x&&(b&c)===x&&(a&d)===x&&(b&d)===x&&(c&d)===x;}
function valid(a){if(new Set(a).size!==a.length||a.some(e=>pop(e)!==3))return false;
 for(let i=0;i<a.length;i++)for(let j=i+1;j<a.length;j++)for(let k=j+1;k<a.length;k++) for(let l=k+1;l<a.length;l++) if(sunflower(a[i],a[j],a[k],a[l]))return false;
 return true;
}
function oldBad(F){
 const bad=[]; for(let i=0;i<F.length;i++)for(let j=i+1;j<F.length;j++)for(let k=j+1;k<F.length;k++){
  const a=F[i],b=F[j],c=F[k],core=a&b;
  if((a&c)===core&&(b&c)===core)bad.push([i,j,k,core,a|b|c]);
 }
 return bad;
}
function runBudget(k){
 const F=mkref();if(!valid(F))throw Error("bad seed");
 // Only fit at most 3(k+1) fresh labels. Stable through 12 fresh labels.
 const ALL=triples(13+3*(k+1)).filter(t=>!F.includes(t));
 const bad=oldBad(F); const bin=new Map(), full=(1n<<39n)-1n;
 const H=F.map((_,i)=>1n<<BigInt(i));
 let cnt=0, keys=new Map(), branches=0;
 // Invert individually admissible insertions to deletion subsets of exact size k.
 for(const t of ALL){
  const need=bad.filter(q=>(t&q[4])===q[3]);
  const visited=new Set(),dels=new Set();
  function allDel(mask,r){
    if(!r){dels.add(mask);return;}
    for(let i=0;i<39;i++)if(!(mask&H[i]))allDel(mask|H[i],r-1);
  }
  function del(mask,n){
    if(visited.has(mask))return;
    visited.add(mask);
    for(const q of need){
      const B=H[q[0]]|H[q[1]]|H[q[2]];
      if(!(B&mask)){
        if(n<k)for(let i=0;i<3;i++)del(mask|H[q[i]],n+1);
        return;
      }
    }
    allDel(mask,k-n);
  }
  del(0n,0);
  if(dels.size) {
   cnt++; const g=pop(t>>>13);keys.set(g,(keys.get(g)||0)+1);
   for(const d of dels){if(!bin.has(d))bin.set(d,[]);bin.get(d).push(t);}
  }
 }
 let maxPool=0, possible=0;
 for(const [d,A] of bin){maxPool=Math.max(maxPool,A.length);
  if(A.length<k+1)continue;
  possible++;
  const keep=F.filter((x,i)=>!(d&H[i])); const add=[];
  function extend(start){
    if(add.length===k+1)return keep.concat(add);
    for(let p=start;p<=A.length-(k+1-add.length);p++){
      const t=A[p]; let ok=true;
      let C=keep.concat(add);
      // At least one previously added edge in each new check.
      for(let i=keep.length;i<C.length&&ok;i++)
       for(let j=0;j<C.length&&ok;j++)if(j!==i)
        for(let l=j+1;l<C.length&&ok;l++)if(l!==i&&sunflower(t,C[i],C[j],C[l]))ok=false;
      if(!ok)continue;
      branches++;add.push(t);const out=extend(p+1);if(out)return out; add.pop();
    }
  }
  const witness=extend(0); if(witness){if(!valid(witness))throw Error("bad output");return {k,witness};}
 }
 return {k,bad:bad.length,total:ALL.length,individ:cnt,keys:[...keys],bins:bin.size,maxPool,possible,branches};
}
for(let k=0;k<=5;k++)console.log(JSON.stringify(runBudget(k)));
