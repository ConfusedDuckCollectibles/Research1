function assert(ok,label){if(!ok)throw Error(label);}
const gcd=(a,b)=>b?gcd(b,a%b):(a<0n?-a:a);
function f(a,b=1n){
  a=BigInt(a);b=BigInt(b);
  assert(b!==0n,"finite");
  if(b<0n){a=-a;b=-b;}
  const d=gcd(a,b);
  return [a/d,b/d];
}
const add=([a,b],[c,d])=>f(a*d+b*c,b*d);
const neg=([a,b])=>[-a,b];
const sub=(a,b)=>add(a,neg(b));
const mul=([a,b],[c,d])=>f(a*c,b*d);
const div=([a,b],[c,d])=>f(a*d,b*c);
const eq=(a,b)=>sub(a,b)[0]===0n;
const arr=(n,h)=>Array.from({length:n},(_,i)=>h(i));
const mat=(n,m,h)=>arr(n,i=>arr(m,j=>h(i,j)));
const dot=(xs,ys)=>xs.reduce((v,x,i)=>add(v,mul(x,ys[i])),f(0));
function inv(A){
  const n=A.length,M=A.map((row,i)=>[...row,...arr(n,j=>f(+(i===j)))]);
  for(let i=0;i<n;i++){
    const j=M.findIndex((row,j)=>j>=i&&row[i][0]);
    if(j<0)return null;
    [M[i],M[j]]=[M[j],M[i]];
    const p=M[i][i];
    M[i]=M[i].map(x=>div(x,p));
    for(let j=0;j<n;j++)if(i!==j){
      const c=M[j][i];M[j]=M[j].map((x,k)=>sub(x,mul(c,M[i][k])));
    }
  }
  return M.map(row=>row.slice(n));
}
function psd(A){
  const M=A.map(r=>r.slice());
  let rank=0;
  for(let i=0;i<M.length;i++){
    const p=M[i][i];
    assert(p[0]>=0n,"positive pivot");
    if(!p[0]){
      assert(M[i].slice(i).every(x=>!x[0]),"zero pivot row");
      continue;
    }
    rank++;
    for(let j=i+1;j<M.length;j++)for(let k=i+1;k<M.length;k++)
      M[j][k]=sub(M[j][k],div(mul(M[j][i],M[i][k]),p));
  }
  return rank;
}
function run(sizes, supports, full=false){
  const r=sizes.length, number=v=>v.reduce((a,i)=>a*BigInt(sizes[i]),1n);
  const cells=supports.map(v=>({v,w:1n,N:number(v)}));
  const stars=arr(r,i=>cells.reduce((a,c)=>a+(c.v.includes(i)?c.N/BigInt(sizes[i]):0n),0n));
  const M=stars.reduce((a,b)=>a>b?a:b);
  for(let i=0;i<r;i++)if(stars[i]<M){
    const old=cells.find(c=>c.v.length===1&&c.v[0]===i);
    if(old)old.w+=M-stars[i];else cells.push({v:[i],w:M-stars[i],N:BigInt(sizes[i])});
  }
  const gram=mat(r,r,(i,j)=>f(cells.reduce((z,o)=>z+(o.v.includes(i)&&o.v.includes(j)?o.w*o.N:0n),0n)));
  const G=inv(gram);
  assert(G,"nonsingular completed group Gram");
  const parts=cells.map(o=>arr(r,i=>f(+o.v.includes(i))));
  const W=mat(cells.length,cells.length,(i,j)=>sub(f(1,M),dot(parts[i],G.map(row=>dot(row,parts[j])))));
  const negative=W.flat().filter(x=>x[0]<0n).length;
  const negativeIntersecting=W.reduce((n,row,i)=>n+row.filter((x,j)=>x[0]<0n&&cells[i].v.some(k=>cells[j].v.includes(k))).length,0);
  for(const [i,row] of W.entries())for(let j=0;j<r;j++){
    const sum=row.reduce((a,x,k)=>add(a,cells[k].v.includes(j)?mul(x,f(cells[k].N*cells[k].w)):f(0)),f(0));
    assert(eq(sum,f(sizes[j]-Number(cells[i].v.includes(j)))),"all support incidence identities");
  }
  let n=0, starActions=0, rank;
  if(full){
    const starts=arr(r,i=>sizes.slice(0,i).reduce((a,v)=>a+v,0));
    const faces=cells.flatMap((c,i)=>arr(Number(c.N),u=>({
      i, points:c.v.map((part,j)=>starts[part]+Math.floor(u/c.v.slice(0,j).reduce((a,k)=>a*sizes[k],1))%sizes[part])
    })));
    n=faces.length;
    const H=mat(n,n,(i,j)=>{
      const x=faces[i],y=faces[j];
      if(i===j)return f(1,cells[x.i].w);
      if(x.points.some(p=>y.points.includes(p)))return f(0);
      return cells[x.i].v.reduce((value,part)=>cells[y.i].v.includes(part)?mul(value,f(sizes[part],sizes[part]-1)):value,W[x.i][y.i]);
    });
    for(let p=0;p<sizes.reduce((a,b)=>a+b,0);p++){
      const star=faces.map(o=>f(o.points.includes(p)?cells[o.i].w:0n));
      assert(H.every(row=>eq(dot(row,star),f(1))),"full weighted star action");
      starActions+=n;
    }
    rank=psd(H.map(row=>row.map(x=>sub(x,f(1,M)))));
  }
  return {sizes,cap:String(M),supports:cells.length,negative,negativeIntersecting,faces:n,starActions,rank};
}
const small=[[0],[1],[2],[0,1],[1,2]];
const triangles=[[0,1,2],[0,1,3],[0,2,3]];
const pairs=Array.from({length:4},(_,i)=>Array.from({length:4-i-1},(_,j)=>[i,i+1+j])).flat();
for(const i of [0,1,2,3,5,6])pairs.push([i,4]);
const mixed=[...arr(7,i=>[i]),...pairs,...triangles];
const results=[run([2,2,2],small,true),run([2,2,2],[[0],[1],[2],[0,1],[0,2],[1,2],[0,1,2]],true),run([2,3,2],[[0],[1],[2],[0,1],[0,2],[1,2],[0,1,2]],true),run(Array(7).fill(2),mixed,true)];
console.log(JSON.stringify(results,null,2));
