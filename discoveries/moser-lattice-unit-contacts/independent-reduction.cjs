"use strict";
function run() {
 const check=(q,m)=>{if(!q)throw Error(m)};
 const rads={20:54,21:57,22:60,23:64,24:68};
 const ranges=n=>[...Array(n).keys()];
 const subsets=n=>ranges(1<<n).map(k=>ranges(n).filter(i=>k>>i&1));
 const pairs=S=>S.flatMap((u,i)=>S.slice(i+1).map(v=>[u,v]));
 const combos=(S,k)=>k===0?[[]]:S.flatMap((u,i)=>combos(S.slice(i+1),k-1).map(A=>[u,...A]));
 const perm=S=>!S.length?[[]]:S.flatMap((u,i)=>perm(S.filter((_,j)=>j!==i)).map(R=>[u,...R]));
 const glyphes=(R,A,t)=>perm(R).map(F=>pairs(F).map(([u,v])=>A[u][v]).join("")+"/"+F.map(u=>t[u]).join("")).sort()[0];
 function graphCensus(n){
  const R=ranges(n), pairsOf=pairs(R), cuts=subsets(n).slice(1), tally={}, smallContactProfiles={}, residualOrbits={};
  let eligible=0, capProfile=0;
  for(let mask=0;mask<1<<pairsOf.length;mask++){
    const A=R.map(()=>R.map(()=>0));let e=0;
    pairsOf.forEach(([u,v],j)=>{if(mask>>j&1){A[u][v]=A[v][u]=1;e++;}});
    const d=A.map(U=>U.reduce((a,b)=>a+b,0));
    if(e>(n===4?5:7)||d.some(x=>x<2))continue;
    if(combos(R,4).some(T=>pairs(T).every(([u,v])=>A[u][v])))continue;
    if(pairsOf.some(([u,v])=>R.filter(w=>A[u][w]&&A[v][w]).length>=3))continue;
    eligible++;
    const cliques=combos(R,3).filter(S=>pairs(S).every(([i,j])=>A[i][j]));
    for(const Q of combos(R,4))if(Q.every(i=>Q.filter(j=>A[i][j]).length>=2))cliques.push(Q);
    // A graph with four vertices, min degree 2, and no K4 is a C4 or a diamond.
    // Both contain a 4-cycle, including when it has a chord.
    const need=cuts.map(T=>[T,73-rads[25-T.length]+pairs(T).reduce((a,[i,j])=>a+A[i][j],0)-T.reduce((a,i)=>a+d[i],0)]);
    function visit(t) {
      if(t.length<n) {for(let k=Math.max(0,5-d[t.length]);k<=3;k++)visit([...t,k]);return;}
      if(need.some(([T,b])=>T.reduce((a,i)=>a+t[i],0)<b))return;
      capProfile++;
      let done=new Set(R.filter(i=>t[i]===3));
      for(let changed=true;changed;){
        changed=false;
        for(const S of cliques){
          const unk=S.filter(i=>!done.has(i));
          if(unk.length===1){done.add(unk[0]);changed=true;}
        }
      }
      const key=e+","+ (n-done.size);
      tally[key]=(tally[key]||0)+1;
      if(done.size<n||t.some(x=>x<2)){
        const c=R.find(i=>d[i]===4);check(n===5&&e===7&&c!==undefined,"gem center");
        const Q=R.filter(i=>i!==c), dd=Q.map(i=>Q.filter(j=>A[i][j]).length), q=Q.filter(i=>Q.filter(j=>A[i][j]).length===1);
        check(dd.slice().sort().join(",")==="1,1,2,2"&&q.every(i=>t[i]===3),"rim is P4 with 3-contact ends");
        if(done.size<n){
          check(t[c]===2&&Q.every(i=>q.includes(i)||t[i]===2),"unpropagated gem t-profile");
          const g=glyphes(R,A,t);residualOrbits[g]=(residualOrbits[g]||0)+1;
        }
        if(t.some(x=>x<2)){
          check(t[c]===1&&Q.every(i=>t[i]>=2),"one contact only at gem hub");
          const g=glyphes(R,A,t);smallContactProfiles[g]=(smallContactProfiles[g]||0)+1;
        }
      }
    }
    visit([]);
  }
  return {eligibleGraphs:eligible,deletionSurvivors:capProfile,tally,smallContactProfiles,residualOrbits};
 }
 const g4=graphCensus(4), g5=graphCensus(5);
 check(JSON.stringify(g4.tally)===JSON.stringify({"4,0":3,"5,0":18}),"four outside");
 const target={"5,0":12,"6,0":210,"7,0":480,"7,3":60};
 check(Object.keys(g5.tally).length===4&&Object.entries(target).every(([k,v])=>g5.tally[k]===v),"five outside");
 const gcd=(a,b)=>b===0?Math.abs(a):gcd(b,a%b);
 const U=new Map();
 for(let s=-2;s<=2;s++)for(let t=-2;t<=2;t++){
   if(!(s>0||(s===0&&t>0))||gcd(s,t)!==1)continue;
   const norm=s*s+s*t+t*t;
   if(norm!==1&&norm!==3)continue;
   for(let u=-7;u<=7;u++)for(let v=-7;v<=7;v++){
     if(norm*(3*u*u+5*u*v+3*v*v)!==27)continue;
     const X=[u*s,u*t,v*s,v*t]; U.set(X.join(","),X);
   }
 }
 const mod3=a=>(a%3+3)%3, key=X=>X.map(mod3).join("");
 const h={};
 for(const a of U.values()){
   const k=key(a); h[k]=(h[k]||0)+1;
   check(k.slice(0,2)==="00"||k.slice(2)==="00","two nonzero blocks forbidden");
 }
 const res=["0000","0011","0022","1100","2200"].map(s=>Array.from(s,Number));
 const unitMatrix=res.map(a=>res.map(b=>h[key(a.map((x,i)=>x-b[i]))]||0));
 check(U.size===42&&Object.values(h).filter(v=>v===1).length===12,"unit residue load");
 check(JSON.stringify(unitMatrix)===JSON.stringify([[18,3,3,3,3],[3,18,3,0,0],[3,3,18,0,0],[3,0,0,18,3],[3,0,0,3,18]]),"unit-residue incidence matrix");
 const A=[[1,0],[2,0],[0,1],[1,1],[0,2]];
 const P=A.flatMap(a=>A.map(b=>a.concat(b))), tally={fiber:0,mixed:0};
 for(const [p,q] of pairs(P)){
   const v=p.map((a,i)=>a-q[i]);
   if(!U.has(v.map(a=>3*a).join(",")))continue;
   tally[(v[0]===0&&v[1]===0)||(v[2]===0&&v[3]===0)?"fiber":"mixed"]++;
 }
 check(tally.fiber===70&&tally.mixed===2,"known 72 witness reverified");
 return {outside:{four:g4,five:g5},normThreeNumerators:U.size,residues:Object.fromEntries(Object.entries(h).sort()),unitMatrix,known25:tally};
}
const out=run();
if(typeof module!=="undefined")console.log(JSON.stringify(out,null,2));
return out;
