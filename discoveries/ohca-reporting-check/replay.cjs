"use strict";
// Source locators and reading scope: README.md.
// Each rate is [numerator, available denominator]; arm order is standard, intervention.
const eq = (a,b) => JSON.stringify(a) === JSON.stringify(b);
const check = (v,why) => { if (!v) throw Error(why); };
const pct = ([n,d]) => Math.floor((200*n+d)/(2*d));
const fmt = a => a.map(pair => [pair[0],pct(pair)]);
const rows={
  prague:{
    // JAMA Table 1, [arm size, ventricular fibrillation, asystole, PEA].
    standard:[132,84,24,24], intervention:[124,72,31,21],
    day180:[[29,132],[39,124]], shockable180:[[28,84],[35,72]]
  },
  inception:{
    // NEJM p.307 Table 4; Results/Fig.1 supply the eligible denominator of 64.
    day30:[[10,62],[14,70]], month3:[[9,63],[12,68]], month6:[[10,63],[14,70]]
  },
  workshops:{
    pragueFirst:[[29,22],[39,31]],
    // Second row percentages use 1 decimal in Prague column. Store them ×10.
    pragueSecondTenths:[[28,333],[35,486]],
    inceptionFirst:[[9,14],[12,18]], inceptionSecond:[[10,16],[14,20]]
  },
  arrest:{dischargeSurvival:[[1,15],[6,14]], dischargeCPC12ILCOR:[[0,15],[3,14]]}
};
let N=0, shock=0, nons=0;
for(let [n,vf,a,p] of [rows.prague.standard,rows.prague.intervention]){
  check(vf+a+p===n,"exhaustive JAMA arm"); N+=n; shock+=vf; nons+=a+p;
}
const nonanalyzed=264-N;
check(nons===100 && shock===156 && pct([nons,N])===39,"Prague arm sums");
check(nonanalyzed===8 && pct([nons+nonanalyzed,264])===41,"Prague maximal rounding");
check((nons+nonanalyzed)*200<97*264,"Prague upper strictly below 48.5%");
check(eq(fmt(rows.prague.day180),rows.workshops.pragueFirst),"Prague whole sample");
check(eq(rows.prague.shockable180.map(([n,d])=>[n,Math.floor((2000*n+d)/(2*d))]),rows.workshops.pragueSecondTenths),"Prague shockable only");
const labelMatches = printed => Object.entries(rows.inception)
  .filter(([,v])=>eq(fmt(v),printed)).map(([k])=>k);
const inceptionMatches=[
  labelMatches(rows.workshops.inceptionFirst),labelMatches(rows.workshops.inceptionSecond)
];
check(eq(inceptionMatches,[["month3"],["day30","month6"]]),"INCEPTION native row matches");
// Independent aggregate check that integer CPC categories 1..4 cannot all be favorable.
// Find all unordered count vectors whose mean rounds to 2.5 and SD to 0.5.
let cpcHistogram={sample:[],population:[]};
for(let a=0;a<=6;a++)for(let b=0;b<=6-a;b++)for(let c=0;c<=6-a-b;c++){
  let d=6-a-b-c, s=a+2*b+3*c+4*d, q=a+4*b+9*c+16*d;
  if(!(49*6<=20*s && 20*s<51*6))continue;
  for(let [mode,den] of [["sample",30],["population",36]]){
    let V=400*(6*q-s*s);if(81*den<=V&&V<121*den)cpcHistogram[mode].push([a,b,c,d]);
  }
}
for(let v of Object.values(cpcHistogram))check(eq(v,[[0,3,3,0]]),"Aggregate CPC bins");
let reconstructed=[[0,15],[3,14]];
check(eq(rows.arrest.dischargeCPC12ILCOR,reconstructed),"Explicit ILCOR favorable row");
check(eq(fmt(rows.arrest.dischargeSurvival),[[1,7],[6,43]]),"Pooled introduction figures trace to discharge survival");
check(eq(fmt(rows.arrest.dischargeCPC12ILCOR),[[0,0],[3,21]]),"Separate outcome");
console.log(JSON.stringify({prague:{standard:[48,132],intervention:[52,124],total:[nons,N],fullRandomizationCeiling:[nons+nonanalyzed,264]},inception:{inceptionMatches,native:rows.inception},arrest:{...rows.arrest,cpcHistogram}},null,2));
