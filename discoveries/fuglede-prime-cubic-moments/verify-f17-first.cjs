"use strict";

function ledger(){
 const oddVal=[]; const out=[];
 for(let m=2;m<=16;m++){
  const quo=17-m, lp=17*quo*quo<256;
  let n=m, fac=[];
  for(let q=2;q<=n;q++)if(n%q===0){let v=0;while(n%q===0){n/=q;v++;}let r=1,j=0;do{r=r*q%17;j++;if(j>16)throw Error("order")}while(r!==16);fac.push([q,v,j]);}
  const norm=m%2===1&&fac.some(([q,v])=>v%2===1);
  if(!lp&&!norm)out.push(m);
  oddVal.push({m,N:17*m,lp,norm,fac});
 }
 if(JSON.stringify(out)!=="[2,4,6,8,9,10,12]")throw Error("ledger");
 return {oddVal,out,sizes:out.map(x=>17*x)};
}

function residues(p,dmax){
 let pow=Array.from({length:dmax+1},(_,d)=>Array.from({length:p},(_,x)=>{let s=1;for(let i=0;i<d;i++)s=s*x%p;return s;}));
 let out=[];
 for(let d=1;d<=dmax;d++){
  const hist=new Int32Array(p), C=new Int32Array(d+1), H=new Int32Array(p), v=new Int32Array(p);
  let min=p*p, minw=null,n2=0,n3=0,runs=0, candidates=[];
  if(d===1){
   out.push({d,min:p*(p-1)/2,n2:0,n3:0,runs:p*(p-1),minw:[0,1],candidates:[]});continue;
  }
  function calc(k){
   if(k===d-1){
    hist.fill(0); let z=0;
    for(let x=0;x<p;x++){const r=v[x]%p; hist[r]++;z+=r;}
    runs++;
    for(let c=0;c<p;c++){
     if(z<min){min=z;minw=Array.from(C);minw[0]=c;}
     if(z<=2*p){n2++; if(candidates.length<100000){const a=Array.from(C);a[0]=c;candidates.push(a);}}
     if(z<=3*p)n3++;
     z+=p-p*hist[p-1-c];
    }
    return;
   }
   const X=pow[k]; for(let i=0;i<p;i++){C[k]=i;calc(k+1);for(let x=0;x<p;x++)v[x]+=X[x];}
   for(let x=0;x<p;x++)v[x]-=p*X[x];
  }
  for(let a=1;a<p;a++){ C[d]=a;for(let x=0;x<p;x++)v[x]=a*pow[d][x];calc(1);}
  out.push({d,min,n2,n3,runs,minw,candidates});
 }
 return out;
}

function septic(p,dmax){
 let pow=Array.from({length:dmax+1},(_,d)=>Array.from({length:p},(_,x)=>{let s=1;for(let i=0;i<d;i++)s=s*x%p;return s;}));
 let out=[];
 for(let d=7;d<=dmax;d++){
  const hist=new Int32Array(p), C=new Int32Array(d+1), H=new Int32Array(p), v=new Int32Array(p);
  let min=p*p, minw=null,n2=0,n3=0,runs=0, candidates=[];
  if(d===1){
   out.push({d,min:p*(p-1)/2,n2:0,n3:0,runs:p*(p-1),minw:[0,1],candidates:[]});continue;
  }
  function calc(k){
   if(k===d-1){
    hist.fill(0); let z=0;
    for(let x=0;x<p;x++){const r=v[x]%p; hist[r]++;z+=r;}
    runs++;
    for(let c=0;c<p;c++){
     if(z<min){min=z;minw=Array.from(C);minw[0]=c;}
     if(z<=2*p){n2++; if(candidates.length<100000){const a=Array.from(C);a[0]=c;candidates.push(a);}}
     if(z<=3*p)n3++;
     z+=p-p*hist[p-1-c];
    }
    return;
   }
   const X=pow[k]; for(let i=0;i<p;i++){C[k]=i;calc(k+1);for(let x=0;x<p;x++)v[x]+=X[x];}
   for(let x=0;x<p;x++)v[x]-=p*X[x];
  }
  for(let a=1;a<=1;a++){ C[d]=a;for(let x=0;x<p;x++)v[x]=a*pow[d][x];calc(1);}
  out.push({d,min,n2,n3,runs,minw,candidates});
 }
 return out;
}

function eightSlopes(list){
 const p=17, r=x=>(x%p+p)%p, mul=(a,b)=>a*b%p;
 const inv=Array(17).fill(0);for(let i=1;i<17;i++)for(let j=1;j<17;j++)if(i*j%17===1)inv[i]=j;
 const by=Object.fromEntries(list.map(c=>[c[6],c]));
 let masks=0,center=0,pass0=0,pass2=0,pass4=0,last=[],sqFans=0;
 function check(s){
  masks++;
  const w=s.map((a,i)=>inv[s.reduce((v,b,j)=>i===j?v:r(v*r(a-b)),1)]);
  let any=false;
  for(let a=1;a<17;a++){
   center++; const G=w.map(u=>by[a*u%p]);
   if(G.some(x=>!x))throw Error("incomplete catalogue");
   if(G.reduce((v,c)=>v+c[0],0)%17!==14)continue;
   pass0++; any=true;
   function zero(d){
    for(let k=0;k<=d;k++){
     let x=0;
     for(let i=0;i<8;i++){let v=1;for(let h=0;h<k;h++)v=v*s[i]%p;x=(x+v*G[i][d])%p;}
     if(x!==0)return false;
    }
    return true;
   }
   if(!zero(4))continue; pass4++;
   if(!zero(2))continue; pass2++;
   last.push({s:s.slice(),a,G});
  }
  if(any)sqFans++;
 }
 let s=[];
 function go(i){
  if(s.length===8){check(s);return;}for(let j=i;j<=17-(8-s.length);j++){s.push(j);go(j+1);s.pop();}
 }
 go(0);
 return {masks,center,pass0,pass4,pass2,sqFans,last};
}

function nineSlopes(L6,L7){
 const P=17, mod=x=>(x%P+P)%P, inv=new Int32Array(P), power=Array.from({length:P},(_,x)=>Array.from({length:8},(_,d)=>{let r=1;for(let j=0;j<d;j++)r=r*x%P;return r;}));
 for(let i=1;i<P;i++)for(let j=1;j<P;j++)if(i*j%P===1)inv[i]=j;
 const evaluate=(c,x)=>c.reduceRight((a,v)=>(a*x+v)%P,0);
 const good=arr=>arr.filter(c=>Math.max(...power.map((_,x)=>evaluate(c,x)))<=9);
 L6=good(L6);L7=good(L7);
 const B6=Object.fromEntries(L6.map(c=>[c[6],c]));
 const B7=Array.from({length:P},()=>new Map());
 for(let u=1;u<P;u++){const W=power[u];
  for(const c of L7){const q=c.map((v,d)=>v*W[d]%P);if(!B7[q[7]].has(q[5]))B7[q[7]].set(q[5],[]); B7[q[7]].get(q[5]).push(q);}
 }
 const bc=Array.from({length:8},(_,d)=>Array.from({length:d+1},(_,j)=>{let r=1;for(let u=1;u<=j;u++)r=r*(d-u+1)*inv[u]%P;return r;}));
 const shifted=Array.from({length:17},()=>[]);
 for(let a=1;a<17;a++)if(B6[a])for(let h=0;h<17;h++){
  const co=B6[a],hp=power[h], out=new Int32Array(7);
  for(let j=0;j<=6;j++)for(let d=j;d<=6;d++)out[j]=(out[j]+bc[d][j]*hp[d-j]*co[d])%P;
  shifted[a][h]=out;
 }
 const s=[];
 const C={masks:0,n6free:0,n6const:0,n6deg4:0,n6deg3:0,n6deg2:0,n6deg1:0,n7free:0,n7fits:0,n7const:0,n7deg4:0,n7deg3:0,n7deg2:0,n7deg1:0,surv6:[],surv7:[]};
 function check(){
  C.masks++;
  const w=s.map((a,i)=>inv[s.reduce((v,b,j)=>i===j?v:mod(v*mod(a-b)),1)]),win=w.map(u=>inv[u]),pow=s.map(x=>power[x]);
  function zero(G,d){
   for(let k=0;k<=d;k++){
    let a=0;for(let i=0;i<9;i++)a+=pow[i][k]*G[i][d];
    if(a%P!==0)return false;
   }
   return true;
  }
  const w3=s.slice(0,3).map((a,i)=>inv[s.slice(0,3).reduce((v,b,j)=>i===j?v:mod(v*mod(a-b)),1)]);
  const lag=s.map(x=>s.slice(0,3).map((y,i)=>s.slice(0,3).reduce((v,b,j)=>i===j?v:mod(v*mod(x-b)),w3[i])));
  const g0=B7[w[0]],g1=B7[w[1]],g2=B7[w[2]];
  for(const a of [...g0.values()].flat())for(const b of [...g1.values()].flat())for(const c of [...g2.values()].flat()){
   C.n7free++;
   let G=[a,b,c], choices=[], all=true;
   const A=a[5]*win[0]%P,B=b[5]*win[1]%P,D=c[5]*win[2]%P;
   for(let i=3;i<9;i++){
    const l=lag[i], r=(A*l[0]+B*l[1]+D*l[2])*w[i]%P;
    const q=B7[w[i]].get(r);if(!q){all=false;break;}choices.push(q);
   }
   if(!all)continue;
   function dfs(k,total){
    if(k<6){for(const q of choices[k]){G.push(q);dfs(k+1,total+q[0]);G.pop();}return;}
    C.n7fits++;
    if(total%P!==16)return;
    C.n7const++;
    for(let d=4;d>=1;d--){if(!zero(G,d))break; C["n7deg"+d]++; if(d===1)C.surv7.push({s:s.slice(),G:G.slice()});}
   }
   dfs(0,a[0]+b[0]+c[0]);
  }
  for(const lead0 of [1,3]){
   const A=lead0*win[0]%P;
   for(let B=0;B<P;B++){
    const a=s.map((x,i)=>(A+B*x)*w[i]%P);
    if(a.some(x=>!B6[x]))continue;
    for(let F=0;F<P;F++){
     C.n6free++;
     const G=[],v5=s.map((x,i)=>F*x*mod(x-1)*w[i]%P);
     for(let i=0;i<9;i++){
      const h=v5[i]*inv[6*a[i]%P]%P;
      const out=shifted[a[i]][h];
      if(out[5]!==v5[i])throw Error("translation");
      G.push(out);
     }
     if(G.reduce((v,x)=>v+x[0],0)%P!==16)continue;
     C.n6const++;
     for(let d=4;d>=1;d--){if(!zero(G,d))break; C["n6deg"+d]++; if(d===1)C.surv6.push({s:s.slice(),G:G.map(q=>[...q])});}
    }
   }
  }
 }
 s.push(0,1);
 function go(i){
  if(s.length===9){check();return;}for(let j=i;j<=17-(9-s.length);j++){s.push(j);go(j+1);s.pop();}
 }
 go(2);
 return {L6:L6.length,L7:L7.length,C};
}
const r=residues(17,6), z=septic(17,7)[0];
const a=eightSlopes(r[5].candidates);const b=nineSlopes(r[5].candidates,z.candidates);
if(JSON.stringify(r.map(x=>x.min))!=="[136,102,85,51,51,34]"||r[5].candidates.length!==16||z.candidates.length!==5||a.pass0!==1360||a.pass4!==0||b.C.n6const!==112842||b.C.n7const!==81||b.C.n6deg4!==0||b.C.n7deg4!==0)throw Error("audit changed");
console.log(JSON.stringify(ledger()));
console.log(JSON.stringify(r.map(({candidates,...rest})=>({...rest,candidates:candidates.length}))));
console.log(JSON.stringify(z));
console.log(JSON.stringify(a));
console.log(JSON.stringify(b));
