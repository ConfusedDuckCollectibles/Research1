'use strict';
const p=17, mod=x=>(x%17+17)%17, xyz=i=>[Math.floor(i/17),i%17];
const iv=Array.from({length:17},(_,n)=>Array.from({length:17},(_,i)=>i).find(i=>i*n%17===1)||0);
const pts=Array.from({length:289},(_,i)=>xyz(i));
const lines=Array.from({length:289},(_,i)=>new Uint16Array(289));
for(let i=0;i<289;i++)for(let j=0;j<i;j++){
 const [x,y]=pts[i],[X,Y]=pts[j];
 const a=x===X?17:mod((y-Y)*iv[mod(x-X)]), b=a===17?x:mod(y-a*x);
 lines[i][j]=lines[j][i]=17*a+b;
}
const tri=[0,17,1], base=(1<<0)|(1<<16)|(1<<17);
function solve(cap,target,positive){
 const opt=Array.from({length:15},(_,i)=>i+1);
 const active=new Uint8Array(306), chosen=tri.slice(), inTri=new Set(tri);
 for(let i=0;i<3;i++)for(let j=0;j<i;j++)active[lines[tri[i]][tri[j]]]=1;
 const stats={masks:0,initialCandidates:0,nodes:0,branches:0,complete:0,maxDepth:3,witness:null,witnessDirections:null};
 function check(mask){
  stats.masks++;
  const cand=[];
  for(let v=0;v<289;v++){
   if(inTri.has(v))continue;
   if(tri.every(u=>!active[lines[u][v]]&&(mask&(1<<Math.floor(lines[u][v]/17)))))cand.push(v);
  }
  stats.initialCandidates+=cand.length;
  function dfs(js){
   stats.nodes++;if(chosen.length>stats.maxDepth)stats.maxDepth=chosen.length;
   if(chosen.length===target){
    stats.complete++;
    if(!stats.witness){
     stats.witness=chosen.map(xyz);
     const s=new Set;for(let i=0;i<chosen.length;i++)for(let j=0;j<i;j++)s.add(Math.floor(lines[chosen[i]][chosen[j]]/17));
     stats.witnessDirections=[...s].sort((a,b)=>a-b);
    }
    return positive;
   }
   for(let j=0;j<=js.length-(target-chosen.length);j++){
    const v=js[j], rs=chosen.map(u=>lines[u][v]);
    // js is already compatible with all current points; mark newly completed 2-point lines.
    for(const r of rs){if(active[r])throw Error('bad maintained filter');active[r]=1}
    const rr=[];
    for(let i=j+1;i<js.length;i++){
     const w=js[i], z=lines[v][w];
     if(!(mask&(1<<Math.floor(z/17)))||active[z])continue;
     if(chosen.every(u=>!active[lines[u][w]]))rr.push(w);
    }
    chosen.push(v);stats.branches++;
    const stop=dfs(rr);chosen.pop();
    for(const r of rs)active[r]=0;
    if(stop)return true;
   }
   return false;
  }
  return dfs(cand);
 }
 function masks(j,n,m){
  if(n===0)return check(m);
  for(let i=j;i<=opt.length-n;i++)if(masks(i+1,n-1,m|(1<<opt[i])))return true;
  return false;
 }
 masks(0,cap-3,base);return stats;
}
const seven=solve(7,7,false), control=solve(8,7,true);
if(seven.masks!==1365||seven.complete||!control.complete||control.witnessDirections.length!==8)throw Error('cap or control');
console.log(JSON.stringify({seven,control},null,2));
