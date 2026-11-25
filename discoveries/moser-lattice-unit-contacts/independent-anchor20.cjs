"use strict";
function audit(nativeGraph6){
  const check=(q,s)=>{if(!q)throw Error(s)};
  const Z=x=>[BigInt(x),0n], add=([a,b],[c,d])=>[a+c,b+d], sub=([a,b],[c,d])=>[a-c,b-d];
  const mul=([a,b],[c,d])=>[a*c+33n*b*d,a*d+b*c], scale=([a,b],n)=>[a*n,b*n], sq=x=>mul(x,x);
  const eq=([a,b],[c,d])=>a===c&&b===d, zero=x=>eq(x,Z(0)), cj=([a,b])=>[a,-b];
  const gcd=(a,b)=>b===0n?(a<0n?-a:a):gcd(b,a%b);
  const psc=(q,n)=>q.map(x=>scale(x,n)), pk=(q,n)=>q.map(x=>mul(x,n)), padd=(x,y)=>x.map((r,i)=>add(r,y[i])), psub=(x,y)=>x.map((r,i)=>sub(r,y[i]));
  const prod=([a,b],[c,d])=>[sub(mul(a,c),scale(mul(b,d),3n)),add(mul(a,d),mul(b,c))];
  const dot=([a,b],[c,d])=>add(mul(a,c),scale(mul(b,d),3n)), sn=p=>dot(p,p), det=([a,b],[c,d])=>sub(mul(a,d),mul(b,c));
  const key=p=>p.flat().join(",");
  const rationalize=(P,D)=>{
    const q=cj(D), d=mul(D,q)[0], coef=pk(P,q);
    const g=coef.flat().concat([d]).reduce(gcd,0n)*(d<0n?-1n:1n);
    return {p:coef.map(a=>a.map(x=>x/g)),d:d/g};
  };
  const A=[[1,0],[2,0],[0,1],[1,1],[0,2]],B=A.slice(0,4),L=72n,L2=L*L;
  const tri=([a,b])=>[Z(2*a+b),Z(b)];
  const want=(i,j)=>{
    const a=i>>2,b=i&3,c=j>>2,d=j&3;
    const u=(v,w)=>eq(sn(psub(tri(v),tri(w))),Z(4));
    return (a===c&&u(B[b],B[d]))||(b===d&&u(A[a],A[c]))||(i===9&&j===6)||(j===9&&i===6);
  };
  const R=[...Array(20).keys()], wantAdj=R.map(i=>R.map(j=>want(i,j)));
  const result=[];
  for(const reflection of [false,true])for(const sign of [1n,-1n]){
    const w=[Z(15),[0n,sign]], S0=[Z(1),Z(-1)], p=A.flatMap(a=>B.map(b=>{
       const Y=tri(b), RC=reflection?prod(S0,[Y[0],scale(Y[1],-1n)]):psc(Y,2n);
       return padd(psc(tri(a),36n),prod(w,RC));
    }));
    check(new Set(p.map(key)).size===20,"20 distinct core points");
    const adj=p.map(a=>p.map(b=>eq(sn(psub(a,b)),Z(L2))));
    check(adj.every((r,i)=>r.every((x,j)=>x===wantAdj[i][j])),"all labelled product-plus-edge adjacencies");
    let radiusOne=0, internal=0, external=0, equilateralExternal=0;
    const centres=new Map();
    for(let i=0;i<20;i++)for(let j=i+1;j<20;j++)for(let k=j+1;k<20;k++){
      const u=psub(p[j],p[i]),v=psub(p[k],p[i]),D=det(u,v);if(zero(D))continue;
      const lengths=[sn(u),sn(v),sn(psub(u,v))];
      if(!eq(lengths.reduce(mul,Z(1)),scale(sq(D),12n*L2)))continue;
      radiusOne++;
      if(R.some(h=>adj[h][i]&&adj[h][j]&&adj[h][k])){internal++;continue;}
      external++;
      if(lengths.every(d=>eq(d,Z(3n*L2))))equilateralExternal++;
      const [U,V]=lengths,H=dot(u,v),D2=scale(sub(mul(U,V),sq(H)),2n);
      const N=padd(pk(p[i],D2),padd(pk(u,mul(V,sub(U,H))),pk(v,mul(U,sub(V,H)))));
      const c=rationalize(N,D2);
      check([i,j,k].every(h=>eq(sn(psub(c.p,psc(p[h],c.d))),Z(L2*c.d*c.d))),"unit-centre direct check");
      centres.set(key(c.p)+";"+c.d,c);
    }
    const C=[...centres.values()];
    const contacts=C.map(c=>p.filter(x=>eq(sn(psub(c.p,psc(x,c.d))),Z(L2*c.d*c.d))).length);
    let antipodal=0;
    for(let i=0;i<C.length;i++)for(let j=i+1;j<C.length;j++){
      const a=C[i],b=C[j],d=a.d*b.d;
      if(eq(sn(psub(psc(a.p,b.d),psc(b.p,a.d))),Z(4n*L2*d*d)))antipodal++;
    }
    check(radiusOne===306&&internal===302&&external===4&&C.length===4&&contacts.every(x=>x===3)&&antipodal===0,"core census");
    result.push({reflection,rotationSign:Number(sign),radiusOne,internal,external,contacts,equilateralExternal,antipodal});
  }
  const record={coreEdges:wantAdj.flat().filter(Boolean).length/2,result};
  if(nativeGraph6){
    check(nativeGraph6.charCodeAt(0)-63===20,"native row order");
    const bits=[...nativeGraph6.slice(1)].flatMap(q=>[5,4,3,2,1,0].map(k=>(q.charCodeAt(0)-63)>>k&1));
    const G=R.map(()=>R.map(()=>false));let b=0;
    for(let j=1;j<20;j++)for(let i=0;i<j;i++)G[i][j]=G[j][i]=bits[b++]===1;
    check(G.flat().filter(Boolean).length===2*record.coreEdges,"native row edge count");
    const degrees=H=>H.map(a=>a.filter(Boolean).length);
    const d=degrees(wantAdj),e=degrees(G), vis=new Set(), f=R.map(()=>-1);
    const extend=()=>{
      const left=R.filter(i=>f[i]<0);
      if(!left.length)return true;
      const domain=i=>R.filter(j=>!vis.has(j)&&d[i]===e[j]&&R.every(k=>f[k]<0||wantAdj[i][k]===G[j][f[k]]));
      const choice=left.map(i=>[i,domain(i)]).sort((a,b)=>a[1].length-b[1].length)[0];
      const [i,ds]=choice;
      for(const j of ds){f[i]=j;vis.add(j);if(extend())return true;f[i]=-1; vis.delete(j);}
      return false;
    };
    check(extend()&&wantAdj.every((a,i)=>a.every((x,j)=>x===G[f[i]][f[j]])),"native row isomorphic to labelled product");
    record.nativeGraph6Permutation=f;
  }
  return record;
}
const input=typeof process==="undefined"?undefined:process.argv[2];
const output=audit(input);
if(typeof module!=="undefined")console.log(JSON.stringify(output,null,2));
return {audit,output};
