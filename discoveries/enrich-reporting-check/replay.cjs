'use strict';
// Reporting-only replay; all values are literal published aggregate cells.
// ENRICH Table 1 (p1284), Results (pp1283–1284), Table 2 (p1286),
// Figure 2 (p1285), Table 3 (p1288); guideline PICO 4.1.3 (pp1036–1037).
const assigned=[[102,106],[48,44]], observed=[[100,99],[47,40]], table2=[[99,100],[40,47]], whole=[147,139];
const add=(x,y)=>x.map((u,i)=>u+y[i]), eq=(x,y)=>x.join()===y.join();
const test=(b,s)=>{if(!b)throw Error(s);};
test(eq(add(...observed),whole),'observed partition of whole arm counts');
test(eq(add(...table2),[139,147]),'printed subgroup sums');
test(table2[1][1]>assigned[1][1],'printed subgroup denominator exceeds randomized count');
test(table2.every((x,j)=>eq(x.slice().reverse(),observed[j])),'swap repairs only subgroup observed columns');
const t0=[[14,150],[27,150]], day180figure=[[30,147],[35,139]];
const guidelineProse=[[30,147],[35,139]], mrs30Coverage=[147,146];
test(guidelineProse.every((x,j)=>eq(x,day180figure[j])&&x[0]>t0[j][0]),'guideline counts match 180 days and exceed 30-day totals');
test(eq(day180figure.map(x=>x[1]),whole)&&!eq(whole,mrs30Coverage),'observed day-180 denominator fingerprint');
function gcd(a,b){return b?gcd(b,a%b):a;}
function report([[a,n],[b,m]]){
 const A=a*(m-b),B=(n-a)*b,g=gcd(A,B);
 const se=Math.sqrt(1/a+1/b+1/(n-a)+1/(m-b)), center=Math.log(A/B);
 return {OR:[A/g,B/g],CI95:[-1,1].map(v=>Math.exp(center+v*1.959963984540054*se).toFixed(4))};
}
console.log(JSON.stringify({assigned,observed,table2,whole,mrs30Coverage,t0,day180figure,day30:report(t0),day180Observed:report(day180figure)},null,2));
