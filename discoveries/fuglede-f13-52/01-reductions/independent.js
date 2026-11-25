// Standalone exact audit of size-52 planar gates, local ray shape and four-minors.
const p=13, mod=n=>(n%p+p)%p;
const inv=Array.from({length:p},(_,x)=>x?Array.from({length:p},(_,y)=>y).find(y=>x*y%p===1):0);
const bitCount=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
function requireThat(x,msg){if(!x)throw Error(msg);}
const dir=(a,b)=>a[0]===b[0]?p:mod((b[1]-a[1])*inv[mod(b[0]-a[0])]);
function capAt(a,b){const s=dir(a,b);return [s,s===p?a[0]:mod(a[1]-s*a[0])];}
function combinations(n,k){
  const z=[];
  function walk(a,b){if(b.length===k){z.push(b.slice());return;}
    for(let x=a;x<=n-k+b.length;x++){b.push(x);walk(x+1,b);b.pop();}
  } walk(0,[]);return z;
}
// Enumerate graph points and grow the ACTUAL direction set, not allowed masks.
// Track line occupancies while adding points; lines with first new pair are activated at 2.
function planar(goal,k,cap){
  const visit=[...Array(goal+1)].fill(0);
  const image=[[0,0],[1,1]], levels=new Int8Array(p*p);
  levels[p]=2;
  let witness=null;
  function visitPart(xmin,used,d){
    const n=image.length;++visit[n];
    if(n===goal){witness=image.map(q=>q.slice());return true;}
    for(let x=xmin;x<=p-goal+n;x++){
      next: for(let y=2;y<p;y++){
        if(used>>y&1)continue;
        const b=[x,y], changes=new Map();let actual=d;
        for(const a of image){
          const [s,h]=capAt(a,b); const id=p*s+h;
          if((levels[id]||1)>=cap)continue next;
          actual|=1<<s;
          if(bitCount(actual)>k)continue next;
          if(!changes.has(id))changes.set(id,levels[id]);
        }
        for(const [id,old] of changes)levels[id]=(old||1)+1;
        image.push(b);
        if(visitPart(x+1,used|1<<y,actual))return true;
        image.pop();
        for(const [id,old] of changes)levels[id]=old;
      }
    } return false;
  }
  visitPart(2,3,2);
  if(witness){
    const dd=new Set(),m=new Map();
    for(const [i,j] of combinations(goal,2)){
      const [s,h]=capAt(witness[i],witness[j]);dd.add(s);
      const key=[s,h].join();if(!m.has(key))m.set(key,new Set());
      m.get(key).add(i).add(j);
    }
    requireThat(dd.size<=k&&Math.max(...[...m.values()].map(z=>z.size))<=cap,"positive control");
  }
  return {goal,k,cap,visit:visit.slice(2),witness};
}
// Normalize two selected base points to (0,0),(1,0); one off point to (0,1)
// by a shear and scaling leaving the x-axis fixed. Check exact line intersections.
const mu=[1,3,9];
function characterizeRay(){
  let shapes=0;const bases=[];
  for(const [i,j] of combinations(11,2)){
    const T=[0,1,i+2,j+2], center=mod(T.reduce((a,b)=>a+b)*inv[4]);
    const base=T.map(x=>[x,0]), first=[0,1];
    const slopes=new Set(T.map(t=>dir([t,0],first)));slopes.add(0);
    requireThat(slopes.size===5,"five exact base directions");
    const more=[];
    for(let x=0;x<p;x++)for(let y=1;y<p;y++)if(x!==0||y!==1){
      const b=[x,y];
      if(base.every(a=>slopes.has(dir(a,b)))&&slopes.has(dir(first,b)))more.push(b);
    }
    for(const [i,j] of combinations(more.length,2))if(slopes.has(dir(more[i],more[j]))){
      const tr=[first,more[i],more[j]];
      requireThat(T.includes(center),"centroid selected");
      requireThat(tr.every(([x,y])=>mod(x-center+center*y)===0),"transversal through centroid");
      requireThat(tr.map(([,y])=>y).sort((a,b)=>a-b).join()===mu.join(),"three-element y coset");
      const k=T.find(t=>t!==center);
      requireThat(mu.every(m=>T.includes(mod(center+m*(k-center)))),"three-element base coset");
      shapes++;bases.push(T);
    }
  }
  requireThat(shapes>0,"ray positive control");
  return {shapes,bases};
}
function fourMinors(){
  const tuple=combinations(p,4), perms=[];
  function walk(a){if(a.length===4){
    let v=0;for(let i=0;i<a.length;i++)for(let j=i+1;j<a.length;j++)if(a[i]>a[j])v++;
    perms.push([a.slice(),(-1)**v]);return;
  }for(let j=0;j<4;j++)if(!a.includes(j)){a.push(j);walk(a);a.pop();}}
  walk([]);
  const mask=a=>a.reduce((s,x)=>s|1<<x,0);
  const reps=tuple.filter(a=>Array.from({length:p},(_,j)=>j).every(j=>mask(a)<=mask(a.map(x=>mod(x+j)))));
  let checked=0;
  requireThat(reps.length===55,"cyclic 4-set representatives");
  for(const a of reps)for(const b of reps){
    const cc=new Int8Array(p);
    for(const [ord,sgn] of perms)cc[mod(a.reduce((v,x,i)=>v+x*b[ord[i]],0))]+=sgn;
    requireThat(cc.some(c=>c!==cc[0]),"primitive 13th Fourier 4-minor");
    checked++;
  }
  return {reps:reps.length,checked};
}
const gates=[[8,5,4],[13,10,4],[11,9,4],[10,9,3],[8,6,4],[13,11,4]].map(args=>planar(...args));
requireThat(gates.every((z,i)=>!!z.witness===(i>=4)),"negative and positive planar scope");
console.log(JSON.stringify({checks:"PASS",gates,ray:characterizeRay(),minors:fourMinors()}));
