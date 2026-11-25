const p=13,mod=n=>(n%13+13)%13,pop=n=>n?n%2+pop(n>>1):0;
const inv=Array.from({length:13},(_,i)=>i?Array.from({length:13},(_,j)=>j).find(j=>i*j%13===1):0);
let maps=[Array.from({length:14},(_,i)=>i===13?13:mod(i+1)),Array.from({length:14},(_,i)=>i===13?13:mod(2*i)),Array.from({length:14},(_,i)=>i===13?0:i===0?13:inv[i])];
function run(cap) {
 const pr=Array.from({length:2197},()=>[]);
 for(let a=0;a<13;a++) for(let b=0;b<13;b++) for(let c=0;c<13;c++){
 const q=pr[a*169+b*13+c];if(!a&&!b&&!c)continue;
 for(let d=0;d<13;d++){
 const v=Array.from({length:13},(_,i)=>mod(((a*i+b)*i+c)*i+d));
 const sum=v.reduce((a,b)=>a+b,0);
 if(v.some(x=>x>cap)||sum!==52)continue;
 q.push({v,d,peak:Math.max(...v)});
 }}
 let seen=new Set(),reps=[];
 for(let m=0;m<8192;m++)if(pop(m)===5&&!seen.has(m)){
 const S=[];for(let i=0;i<13;i++)if(m>>i&1)S.push(i);
 const queue=[m];seen.add(m);for(let j=0;j<queue.length;j++)for(const f of maps){
 const k=f.reduce((a,b,i)=>a|(((queue[j]>>i)&1)<<b),0);
 if(!seen.has(k)){queue.push(k);seen.add(k);}
 }
 const w=S.map((s,i)=>inv[S.reduce((a,t,j)=>i===j?a:mod(a*(s-t)),1)]);
 const arg=Array.from({length:169},(_,q)=>S.map(t=>mod(q%13-t*(q/13|0))));
 const stats={S,q3:0,q2:0,mod:0,exact:0,sigs:{},samples:[]};reps.push(stats);
 function test(a,b,c,u,v,z,deg){
 const lists=S.map((s,i)=>pr[169*mod(w[i]*a)+13*mod(w[i]*(b+c*s))+mod(w[i]*(u+v*s+z*s*s))]);if(lists.some(x=>!x.length))return;
 const g=[];
 function dfs(i,total) {
 if(i===5){
  if(mod(total-16))return;stats.mod++;
  const f=[];
  for(const z of arg){
    let v=g.reduce((a,q,j)=>a+q.v[z[j]],0)-16;
    if(v<0||v%13||v>13*Math.floor(5*(cap-4)/13))return;
    f.push(v/13);
  }

  for(let t=0;t<14;t++){
    const r=Array(13).fill(0), j=S.indexOf(t);
    f.forEach((v,q)=>r[t===13?(q/13|0):mod(q%13-t*(q/13|0))]+=v);
    if(!r.every((v,i)=>v===(j<0?4:g[j].v[i])))throw Error("profile recovery");
  }
  stats.exact++;stats[deg]++;
  const k=g.map(x=>x.v.reduce((s,n)=>s+(n-4)**2,0)/26).sort((a,b)=>a-b),counts=[f.filter(x=>x===1).length,f.filter(x=>x===2).length,f.filter(x=>x===3).length];
  const peak=Math.max(...g.map(x=>x.peak));
  const key=[deg,peak,...counts,...k].join(',');
  stats.sigs[key]=(stats.sigs[key]??0)+1;
  if (stats.samples.length<2)stats.samples.push({coeff:[a,b,c,u,v,z],d:g.map(x=>x.d),key});
  return;
 }
 for(const q of lists[i]){g.push(q);dfs(i+1,total+q.d);g.pop();}
 }
 dfs(0,0);
 }
 for(let a=1;a<13;a++)for(let u=0;u<13;u++)for(let v=0;v<13;v++)for(let z=0;z<13;z++)test(a,0,0,u,v,z,"q3");
 for(let b=0;b<13;b++)for(let c=0;c<13;c++)if(b||c)for(let k=0;k<13;k++)test(0,b,c,c?k:0,0,c?0:k,"q2");
 }
 return reps;
}const slope=Array.from({length:169},(_,a)=>Array.from({length:169},(_,b)=>((a%13-b%13+13)*inv[(Math.floor(a/13)-Math.floor(b/13)+13)%13])%13));
const planarGate=(goal,directions,cap,one=false)=>{
const masks=[]; for(let m=2;m<8192;m+=4)if(pop(m)===directions)masks.push(m);
const visits=Array(goal+1).fill(0);let model=[];
for(const mask of masks) {
const pts=[0,14];
function search(nextx,yused){
const n=pts.length;visits[n]++;if(n===goal){if(!model.length)model=pts.slice();return one;}
for(let x=nextx;x<=p-(goal-n);x++)for(let y=2;y<p;y++)if(!((yused>>y)&1)){
 const b=x*13+y, u=Array(13).fill(0);
 if(pts.some(a=>{const d=slope[a][b];return !((mask>>d)&1)||++u[d]>=cap;}))continue;
 pts.push(b);const stop=search(x+1,yused|(1<<y));pts.pop();if(stop)return true;
}
}
if(search(2,3))break;
}
return {goal,directions,cap,masks:masks.length,visits:visits.slice(2),model};
};

function checkQuadratic(){
 const rows=[];
 for(let a=1;a<13;a++)for(let d=0;d<13;d++){
  const profile=Array.from({length:13},(_,x)=>mod(a*x*x+d));
  const m=profile.reduce((s,v)=>s+v,0), z=profile.filter(x=>!x).length;
  if(m<=52&&m+13*z>=52&&(52-m)%13===0) rows.push([a,d,(52-m)/13,Math.max(...profile)]);
 }
 if(rows.length!==12||!rows.every(([a,d,elev])=>elev===0&&[4,8].includes(d)))throw Error("quadratic profiles");
 return rows;
}
console.log(JSON.stringify({quadratic:checkQuadratic()}));
for(const g of [[8,5,4],[11,9,4],[13,10,4],[10,9,3],[8,6,4],[13,11,4]])console.log(JSON.stringify(planarGate(...g,g[1]===6||g[1]===11)));
for(const cap of [9,10,11,12]) console.log(JSON.stringify({cap,replay:run(cap)}));
