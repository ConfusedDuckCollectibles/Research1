"use strict";
// Independent direct-height and quotient-space check for the F_13^3 K9 gate.
// Run with Node.js, or a V8 evaluator supplying console.log. No extra packages.
const B = [
  0,1,3,9,15,17,21,22,26,33,35,37,40,45,50,51,52,56,60,63,65,68,73,76,78,79,
  85,86,91,92,99,101,104,107,108,115,120,122,124,127,130,137,138,140,143,150,
  152,155,156,163,167,168
].map(z => [Math.floor(z/13), z%13]);
const field = Array.from({length:13}, (_,i) => i);
const p = x => (x%13+13)%13;
const inverse = [0,...field.slice(1).map(x => field.find(y => x*y%13===1))];
const zero = n => Array(n).fill(0);
const e = (n,j) => zero(n).map((a,i) => +(i===j));
const sum = (x,y,s=1) => x.map((v,i) => p(v+s*y[i]));
const scalar = (x,s) => x.map(v => p(s*v));
const dot = (x,y) => p(x.reduce((s,v,i) => s+v*y[i],0));
const eq = (x,y) => x.length===y.length && x.every((v,i) => v===y[i]);
function invariant(ok,msg) { if(!ok) throw Error(msg); }
const supports = [0,2].map(m => B.flatMap(([x,y],i) => p(y-m*x)===0?[i]:[]));
invariant(B.length===52 && new Set(B.map(String)).size===52,"invalid grid");
invariant(supports.every(a=>a.length===10) &&
          eq(supports[0].filter(i=>supports[1].includes(i)),[0]),"sections");

function census(indices) {
  const x=indices.map(i => B[i][0]), h=[0,0], accepted=[], visits=zero(11);
  const delta=x.map(a => x.map(b=>inverse[p(a-b)]));
  function walk(isScaled,old) {
    const i=h.length;
    visits[i]++;
    if(i===x.length) {
      invariant(isScaled,"nonzero residual");
      accepted.push({h:h.slice(),d:old.size});
      return;
    }
    for(const v of (isScaled?field:[0,1])) {
      const D=new Set(old), newLine = new Map();
      let valid=true;
      for(let j=0;j<i;j++) {
        const slope=p((v-h[j])*delta[i][j]);
        D.add(slope);
        const line=13*slope+p(v-slope*x[i]);
        newLine.set(line,(newLine.get(line)||1)+1);
        if(D.size>9 || newLine.get(line)>4) {valid=false; break;}
      }
      if(valid){h.push(v); walk(isScaled||!!v,D); h.pop();}
    }
  }
  walk(false,new Set([0]));
  return {x,visits:visits.slice(2),shapes:accepted};
}
// Column-space elimination with optional linear origin records. This is not
// row reduction of an augmented moment matrix; it incrementally spans columns.
class Space {
  constructor(n) {this.n=n; this.entries=[];}
  reduce(x, label) {
    let r=x.slice(), a=label?.slice();
    for(const {v,q,origin} of this.entries) if(r[q]) {
      const c=r[q];r=sum(r,v,-c);
      if(a) a=sum(a,origin,-c);
    }
    return [r,a];
  }
  add(x,label) {
    const [r,a]=this.reduce(x,label), q=r.findIndex(Boolean);
    if(q<0) return {dependent:true,relation:a};
    const inv=inverse[r[q]];
    this.entries.push({v:scalar(r,inv),q,origin:a&&scalar(a,inv)});
    return {dependent:false};
  }
  get rank(){return this.entries.length;}
}
const mom = D => Array.from({length:D+1},(_,d)=>
  Array.from({length:d+1},(_,b)=>B.map(([x,y])=>p(x**(d-b)*y**b)))).flat();
const transpose = A => A[0].map((_,i)=>A.map(r=>r[i]));
const P=mom(6), Q=mom(5);
const linearColumns=transpose(P);
const pinned=new Set(supports.flat());
const free=B.flatMap((_,i)=>pinned.has(i)?[]:[i]).reverse();
const allLin=new Space(P.length);
for(const c of linearColumns) allLin.add(c);
const lin=new Space(P.length), kernel=[];
for(const i of free){
  const v=lin.add(linearColumns[i],e(B.length,i));
  if(v.dependent) kernel.push(v.relation);
}
invariant(allLin.rank===28 && lin.rank===27 && kernel.length===6,"linear rank");
invariant(kernel.every(k=>[...pinned].every(i=>!k[i]) && P.every(w=>!dot(w,k))),"kernel");
const moments=(A,x)=>A.map(w=>dot(w,x));
const squares=x=>moments(Q,x.map(v=>p(v*v)));
const Qspace=new Space(Q.length);
const quadratic=[];
for(let i=0;i<kernel.length;i++) for(let j=i;j<kernel.length;j++) {
  const v=moments(Q,kernel[i].map((a,z)=>p(a*kernel[j][z]*(i===j?1:2))));
  quadratic.push(v);Qspace.add(v);
}
invariant(Qspace.rank===16,"pure quadratic rank");
const qreduce=v=>Qspace.reduce(v)[0];
const [c0,c2]=supports.map(census);
invariant(c0.shapes.length===6&&c2.shapes.length===6,"two shape counts");
let n=0,survivors=0;
const hist={}, certs=[];
function validAnchor(h) {
  const rhs=scalar(moments(P,h),-1);
  const [left,negSolution]=lin.reduce(rhs,zero(B.length));
  invariant(!left.some(Boolean),"linear incompatibility");
  const h0=sum(h,negSolution,-1);
  invariant([...pinned].every(i=>h[i]===h0[i]) && P.every(w=>!dot(w,h0)),"bad particular solution");
  const L=kernel.map(k=>moments(Q,k.map((x,i)=>p(2*x*h0[i]))));
  const b=scalar(squares(h0),-1), post=new Space(Q.length);
  for(const v of L)post.add(qreduce(v));
  const [rem]=post.reduce(qreduce(b));
  if(!rem.some(Boolean)) return {rank:Qspace.rank+post.rank, valid:true};
  const k=rem.findIndex(Boolean);
  const leftFunc=Array.from({length:Q.length},(_,j)=>{
    const [r]=post.reduce(qreduce(e(Q.length,j)));
    return p(r[k]*inverse[rem[k]]);
  });
  invariant(dot(leftFunc,b)===1 && quadratic.concat(L).every(c=>!dot(leftFunc,c)),
    "bad explicit certificate");
  return {rank:Qspace.rank+post.rank, valid:false, leftFunc};
}
for(const f of c0.shapes)for(const g of c2.shapes)for(const t of field.slice(1)){
  const h=zero(B.length);
  supports[0].forEach((j,i)=>h[j]=f.h[i]);
  supports[1].forEach((j,i)=>h[j]=p(t*g.h[i]));
  const r=validAnchor(h);
  n++;hist[r.rank]=(hist[r.rank]||0)+1;
  if(r.valid)survivors++;else certs.push(r.leftFunc);
}
// Non-exclusion control: an affine zero height on B really solves all h and h² moments.
// In other words the algorithm does not mark every right-hand side inconsistent.
const control=validAnchor(zero(B.length));
invariant(control.valid&&control.rank===16,"vanishing-height control");
invariant(n===432&&survivors===0&&hist[16]===12&&hist[18]===120&&hist[19]===300,"anchor output");
console.log(JSON.stringify({c0,c2,linear:[allLin.rank,lin.rank,kernel.length]}));
console.log(JSON.stringify({pureQuadraticRank:Qspace.rank,n,hist,survivors,certificates:certs.length,
  certificateKinds:new Set(certs.map(String)).size,control}));
