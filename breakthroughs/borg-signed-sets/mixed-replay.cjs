function assert(ok,msg){if(!ok)throw Error(msg);}
function run(sizes, supports){
  const r=sizes.length, str=v=>v.join(",");
  const cells=supports.map(([v,w])=>({v,w}));
  const faces=[];
  for(const [cell,{v,w}] of cells.entries()){
    function add(labels,j){
      if(j===v.length){faces.push({cell,w,labels});return;}
      for(let a=0;a<sizes[v[j]];a++)add([...labels,[v[j],a]],j+1);
    }
    add([],0);
  }
  const stars=sizes.map((_,i)=>cells.reduce((n,o)=>n+(o.v.includes(i)?o.w*o.v.reduce((z,j)=>z*(i===j?1:sizes[j]),1):0),0));
  const M=Math.max(...stars), N=faces.length;
  const bit=Array.from({length:N},(_,i)=>1n<<BigInt(i)), ix=new Map(bit.map((b,i)=>[b,i]));
  const adj=faces.map((x,i)=>faces.reduce((z,y,j)=>i!==j&&x.labels.some(([c,a])=>y.labels.some(([d,b])=>c===d&&a===b))?z|bit[j]:z,0n));
  const classes=new Map();
  for(let i=0;i<r;i++)if(sizes[i]===2&&stars[i]===M){
    const key=str(cells.flatMap((o,j)=>o.v.includes(i)?[j]:[]));
    const c=classes.get(key)||[];
    c.push(i);classes.set(key,c);
  }
  let optimal=0,nonstar=0;
  function check(choice){
    optimal++;
    const selected=faces.filter((_,i)=>(choice&bit[i])!==0n);
    if(selected[0].labels.some(([i,a])=>selected.every(x=>x.labels.some(([j,b])=>i===j&&a===b))))return;
    nonstar++;
    let good=false;
    for(const C of classes.values())if(C.length>=3){
      const all=faces.filter(x=>C.every(i=>x.labels.some(([j])=>i===j)));
      const chosen=all.filter(x=>selected.includes(x));
      if(chosen.length!==selected.length)continue;
      const key=x=>str(C.map(i=>x.labels.find(([j])=>i===j)[1]));
      const flip=x=>x.split(",").map(a=>1-Number(a)).join(",");
      const accept=new Set(chosen.map(key));
      if(all.every(x=>selected.includes(x)===accept.has(key(x)))&&[...new Set(all.map(key))].every(k=>accept.has(k)!==accept.has(flip(k))))good=true;
    }
    assert(good,"mixed nonstar classification");
  }
  function clique(chosen,P,X,value){
    if(!P&&!X){assert(value<=M,"largest star");if(value===M)check(chosen);return;}
    const union=P|X,pivot=ix.get(union&-union);
    let todo=P&~adj[pivot];
    while(todo){
      const b=todo&-todo,j=ix.get(b);
      clique(chosen|b,P&adj[j],X&adj[j],value+faces[j].w);
      P^=b;X|=b;todo^=b;
    }
  }
  clique(0n,(1n<<BigInt(N))-1n,0n,0);
  assert(optimal>0,"some optimum");
  assert(Boolean(nonstar)===[...classes.values()].some(C=>C.length>=3),"binary class iff nonstar");
  return {sizes,N,M,optimal,nonstar};
}
const all=Array.from({length:7},(_,i)=>Array.from({length:3},(_,j)=>j).filter(j=>((i+1)>>j)&1));
let checked=0, optimal=0;
for(const sizes of [[2,2,3],[2,3,4]])for(let mask=1;mask<128;mask++){
  const out=run(sizes,all.flatMap((S,i)=>(mask>>i)&1?[[S,1+((mask+i)%3)]]:[]));
  checked++;
  optimal+=out.optimal;
}
const special=[run([2,2,2,3],[[[0,1,2,3],1]]),run([2,2,2,3],[[[0,1,2,3],2],[[0,1,2],3]])];
console.log(JSON.stringify({checked,optimal,special},null,2));
