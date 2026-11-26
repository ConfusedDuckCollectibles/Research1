"use strict";
const b=i=>1<<i, m=a=>a.reduce((s,i)=>s|b(i),0);
const groups=[[0,1,2],[3,4,5],[6,7,8]], pairs=[[0,1,10],[0,2,11],[1,2,12]];
const E=groups.map(m);
for(const a of groups)for(let i=0;i<3;i++)for(let j=i+1;j<3;j++)E.push(m([9,a[i],a[j]]));
for(const [p,q,s] of pairs)for(const a of groups[p])for(const t of groups[q])E.push(m([s,a,t]));
if(E.length!==39||new Set(E).size!==39)throw Error("wrong source witness");
const forb=[];
for(let i=0;i<E.length;i++)for(let j=i+1;j<E.length;j++)for(let k=j+1;k<E.length;k++){
 const c=E[i]&E[j]; if(c===(E[i]&E[k])&&c===(E[j]&E[k]))
  forb.push({union:E[i]|E[j]|E[k],core:c,orig:[i,j,k],mask:(1n<<BigInt(i))|(1n<<BigInt(j))|(1n<<BigInt(k))});
}
const types=[{m:0,n:0}];
for(let i=0;i<13;i++){
 types.push({m:b(i),n:1});
 for(let j=i+1;j<13;j++){
   types.push({m:b(i)|b(j),n:2});
   for(let k=j+1;k<13;k++){
     const tm=b(i)|b(j)|b(k); if(!E.includes(tm))types.push({m:tm,n:3});
   }
 }
}
const pop=x=>{let z=0;while(x){x&=x-1n;z++;}return z;};
const result=[];
for(const t of types){
 const ban=forb.filter(x=>(x.union&t.m)===x.core);
 const sols=new Set,seen=new Set;
 function visit(S,d){
   if(seen.has(S))return;
   seen.add(S);
   let node=null,union=0n,lb=0;
   // Greedy disjoint packing gives a lower bound. A generic exact first-edge branching follows.
   for(const x of ban){
     if((x.mask&S)!==0n)continue;
     node??=x;
     if((union&x.mask)===0n){union|=x.mask;lb++;}
     if(lb+d>5)return;
   }
   if(!node){sols.add(S);return;}
   for(const j of node.orig)visit(S|1n<<BigInt(j),d+1);
 }
 visit(0n,0);
 t.mins=[...sols].sort((a,b)=>pop(a)-pop(b)).filter((s,i,ar)=>!ar.slice(0,i).some(q=>(q&s)===q));
 const min=t.mins.length?pop(t.mins[0]):6;
 result.push({N:t.n,min,sols:t.mins.length,type:t.m,blocks:ban.length});
}
const hist={};for(const a of result){const tag=a.N+":"+a.min; hist[tag]=(hist[tag]||0)+1;}

const indexes=Array.from({length:39},(_,i)=>1n<<BigInt(i));
const twoOrMore=(a,b,c,d)=>{
  const k=a&b;
  return (a&c)===k&&(a&d)===k&&(b&c)===k&&(b&d)===k&&(c&d)===k;
};
// Independent extra check: literal 4-tuples, no sunflower in the reproduced source family.
for(let i=0;i<39;i++)for(let j=i+1;j<39;j++)for(let k=j+1;k<39;k++)for(let l=k+1;l<39;l++)
 if(twoOrMore(E[i],E[j],E[k],E[l]))throw Error("source 4-sunflower");
const Dstats=[];
for(let k=0;k<=5;k++){
 const D=new Map;
 const aim=k+1;
 for(let ti=0;ti<types.length;ti++){
  const t=types[ti], ds=new Set;
  function add(S,d,i0){
    if(!d){ds.add(S);return;}
    for(let i=i0;i<=39-d;i++)if((S&indexes[i])===0n)add(S|indexes[i],d-1,i+1);
  }
  for(const s of t.mins)if(pop(s)<=k)add(s,k-pop(s),0);
  for(const s of ds){
   if(!D.has(s))D.set(s,[]);
   D.get(s).push(ti);
  }
 }
 let large=0, nodes=0, maxAllOld=0, hasFresh=0, maxPartial=0, maxPartialFresh=0;
 let histFreshTypes={}, histLargeOld={};
 for(const [S,Ts] of D){
   const nf=Ts.filter(i=>types[i].n===2).length;
   if(Ts.some(i=>types[i].n<2))throw Error("multiple fresh points possible");
   if(nf)hasFresh++;
   maxAllOld=Math.max(maxAllOld,Ts.length-nf);
   histFreshTypes[nf]=(histFreshTypes[nf]||0)+1;
   if(!nf&&Ts.length<aim)continue;
   large++;
   if(!nf)histLargeOld[Ts.length]=(histLargeOld[Ts.length]||0)+1;
   const A=E.filter((_,i)=>(S&indexes[i])===0n),n=A.length;
   let partial=0;
   function valid(t){
     // Enumerate all previous triples with at least one new set; wholly-old is already handled.
     for(let l=n;l<A.length;l++)for(let i=0;i<l;i++)for(let j=i+1;j<l;j++)
       if(twoOrMore(t,A[i],A[j],A[l]))return false;
     return true;
   }
   function visit(goal,nFresh,nextType,lastLabel){
    nodes++;
    partial=Math.max(partial,goal);
    if(goal===aim)throw Error(JSON.stringify({S:String(S),A: A.slice(n)}));
    for(let j=nextType;j<Ts.length;j++){
     const t=types[Ts[j]];
     if(t.n===3){
       const e=t.m;
       if(valid(e)){A.push(e);visit(goal+1,nFresh,j+1,-1);A.pop();}
     }else{
       // Label by first occurrence independent of support size. Sort repeats within a type.
       const low=(j===nextType)?lastLabel+1:0;
       for(let f=low;f<=nFresh;f++){
         const e=t.m|b(13+f);
         if(valid(e)){A.push(e);visit(goal+1,Math.max(nFresh,f+1),j,f);A.pop();}
       }
     }
    }
   }
   visit(0,0,0,-1);
   maxPartial=Math.max(maxPartial,partial);
   if(nf)maxPartialFresh=Math.max(maxPartialFresh,partial);
 }
 Dstats.push({k,bins:D.size,large,nodes,maxAllOld,hasFresh,histFreshTypes,histLargeOld,maxPartial,maxPartialFresh});
}
console.log(JSON.stringify({Dstats,hist,nForb:forb.length,exactNoForty:true},null,2));
