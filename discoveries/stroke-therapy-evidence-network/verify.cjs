"use strict";
const direct = "C21/T21";
const base = [
  ["CHANCE","A","C21"], ["POINT","A","C90"], ["THALES","A","T30"],
  ["SOCRATES","A","TM"], ["CHANCE-2","C21","T21"], ["PRINCE","C21","T21"]
];
const key = (a,b) => [a,b].sort().join("/");
function graph(inspires, mergeTicagrelor=false) {
  const adj = new Map(), original = new Map();
  const node = k => mergeTicagrelor && k==="T30" ? "T21" : k;
  for(const [trial, aa, bb] of [...base,["INSPIRES","A",inspires]]) {
    const a=node(aa), b=node(bb), edge=key(a,b);
    if(!adj.has(a)) adj.set(a,new Set());
    if(!adj.has(b)) adj.set(b,new Set());
    if(!original.has(edge)) original.set(edge,[]);
    original.get(edge).push(trial);
    if(edge===direct) continue;
    adj.get(a).add(b); adj.get(b).add(a);
  }
  const seen=new Set(["T21"]), queue=["T21"];
  for(const n of queue) for(const nxt of adj.get(n)) {
    if(!seen.has(nxt)) { seen.add(nxt); queue.push(nxt); }
  }
  return {nodes:adj.size, distinctEdges:original.size, heldOutTrials:original.get(direct),
    reachableFromT21:[...seen].sort(), indirectPath:seen.has("C21")};
}
const published=graph("C90"), reclassified=graph("C21"), coarsePositiveControl=graph("C21",true);
if(published.nodes!==6 || published.distinctEdges!==5 || published.indirectPath ||
    reclassified.nodes!==6 || reclassified.distinctEdges!==5 || reclassified.indirectPath ||
    !coarsePositiveControl.indirectPath) throw Error("Graph check failed");
console.log(JSON.stringify({published,reclassified,coarsePositiveControl},null,2));
