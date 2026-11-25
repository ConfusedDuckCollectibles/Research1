"use strict";
const p=13,mod=x=>(x%13+13)%13,inv=a=>Array.from({length:12},(_,i)=>i+1).find(x=>mod(a*x)===1);
const pop=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
const support=k=>Array.from({length:p},(_,i)=>i).filter(i=>k>>i&1);
const assert=(x,s)=>{if(!x)throw Error(s);};
let total=0;
for(const sm of [111,187,215,311,427]){
 const S=support(sm),d=S.map((a,i)=>S.map((b,j)=>i===j?0:inv(mod(a-b))));
 const ambig=sm===187||sm===427,pair=ambig?[3,5]:[0,1];
 const ia=S.indexOf(pair[0]),ib=S.indexOf(pair[1]),count={},model=[0];
 let palettes=0,maxWidths=0,recursive=0;assert(ia>=0&&ib>=0,'fixed labels exist');
 for(let I=0;I<1<<p;I++)if(pop(I)===5){
  const C=support(I),width=new Set();let hits=0;
  function dfs(k,K){
   recursive++;
   if(k===6){
    if(K===I){hits++;width.add(mod(model[ib]-model[ia]));}return;
   }
   for(const w of C){
    const y=mod(w*S[k]);let M=K,ok=true;
    for(let j=0;j<k;j++){
     const q=1<<mod((y-model[j])*d[k][j]);
     if(!(I&q)){ok=false;break;}
     M|=q;
    }
    if(ok){model.push(y);dfs(k+1,M);model.pop();}
   }
  }
  dfs(1,0);
  count[hits]=(count[hits]??0)+1;if(hits){
   palettes++;total++;maxWidths=Math.max(maxWidths,width.size);
   assert(hits===(ambig?2:1)&&width.size===1,'fixed-pair contrast for support '+sm+' in palette '+I);
  }
 }
 assert(palettes===(ambig?78:156)&&maxWidths===1,'census');
 console.log(JSON.stringify({supportMask:sm,support:S,pair,palettes,count,maxWidths,recursive}));
}
assert(total===624,'nonempty family census');
