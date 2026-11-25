"use strict";
const assert=(cond,msg)=>{if(!cond) throw Error(msg);};
function census(moved) {
  const fixed=40-moved, tau=i=>i^1;
  const edge=(i,j)=>Math.min(i,j)*moved+Math.max(i,j);
  const orbits=new Map();
  for(let i=0;i<moved;i++)for(let j=i+1;j<moved;j++){
    const e=edge(i,j),f=edge(tau(i),tau(j));
    orbits.set(Math.min(e,f),[e,f]);
  }
  const O=[...orbits.values()], byFamily=[];
  const results={fixed,orbitCount:O.length,masks:1<<O.length,triangleFree:0,passesPairNecessities:0,degreeHist:{}};
  for(let mask=0;mask<results.masks;mask++){
    const G=Array.from({length:moved},()=>Array(moved).fill(0));
    for(let j=0;j<O.length;j++)if(mask>>j&1)for(const e of O[j]){
      const u=Math.floor(e/moved),v=e%moved; G[u][v]=G[v][u]=1;
    }
    let triangle=false;
    for(let i=0;i<moved;i++)for(let j=i+1;j<moved;j++)for(let k=j+1;k<moved;k++)
      if(G[i][j]&&G[j][k]&&G[k][i])triangle=true;
    if(triangle)continue;
    results.triangleFree++;
    const ds=G.map(a=>a.reduce((s,x)=>s+x,0));
    const options=[];
    for(let i=0;i<moved;i+=2){
      const internal=G[i][i+1]; let k=0;
      for(let j=0;j<moved;j++)if(j!==i&&j!==i+1&&(G[i][j]||G[i+1][j]))k++;
      const ms=[];
      for(let m=0;m<=9;m++){
        if(m+ds[i]<4||m+ds[i]>9)continue;
        if(internal ? m!==0 : 38-k-m>=28)continue;
        ms.push(m);
      }
      options.push(ms);
    }
    if(options.some(ms=>ms.length===0))continue;
    results.passesPairNecessities++;
    const key=ds.join(",");
    results.degreeHist[key]=(results.degreeHist[key]||0)+1;
    byFamily.push({G,options});
  }
  return {results,byFamily};
}
const maps=[2,4,6].map(census);
assert(maps[0].byFamily.length===0&&maps[1].byFamily.length===0,"involutions with fewer than six moved");
assert(maps[2].byFamily.length===4,"six moved type census");
for(const {G,options} of maps[2].byFamily){
  assert(options.every(ms=>ms.length===1&&ms[0]===7),"fixed-neighborhood size");
  const seen=new Set([0]);
  for(const i of seen)G[i].forEach((v,j)=>{if(v)seen.add(j);});
  assert(seen.size===6&&G.every(a=>a.reduce((s,x)=>s+x,0)===2),"C6");
}
let total=0,max=-1,maxCount=0,witness=null;
for(let a=0;a<=73;a++)for(let b=0;b<=117;b++)for(let c=0;c<=42;c++){
  if(2*a+b>147||b+2*c>117)continue;
  total++;
  const cards=a+2*b+3*c;
  if(cards>max){max=cards;maxCount=0;witness=[a,b,c];}
  if(cards===max)maxCount++;
}
assert(max===249&&3*85>max&&3*84>max,"three-card budget");
const out={movedGraphCensus:maps.map(x=>x.results),budget:{integerTriples:total,max,maximizers:maxCount,firstMaximizer:witness,minimumUsingPublishedInput:255,margin:6}};
console.log(JSON.stringify(out,null,2));
return out;
