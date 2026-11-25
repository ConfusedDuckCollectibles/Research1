"use strict";
function audit(){
  const assert=(v,s)=>{if(!v)throw Error(s)};
  const I=n=>[BigInt(n),0n], add=([a,b],[c,d])=>[a+c,b+d], neg=([a,b])=>[-a,-b], sub=(x,y)=>add(x,neg(y));
  const mul=([a,b],[c,d])=>[a*c+33n*b*d,a*d+b*c], sc=([a,b],n)=>[a*BigInt(n),b*BigInt(n)];
  const conj=([a,b])=>[a,-b], sq=x=>mul(x,x), equal=([a,b],[c,d])=>a===c&&b===d;
  const z=x=>equal(x,I(0));
  const sgn=x=>x<0n?-1:x>0n?1:0;
  const sgnK=([a,b])=>{
    if(a===0n)return sgn(b);
    if(b===0n||sgn(a)===sgn(b))return sgn(a);
    return sgn(a)*sgn(a*a-33n*b*b);
  };
  const gcd=(a,b)=>b===0n?(a<0n?-a:a):gcd(b,a%b);
  // Point columns correspond to physical (X/36, sqrt(3)*Y/36).
  const cols=[[[36n,0n],[0n,0n]],[[18n,0n],[18n,0n]],[[30n,0n],[0n,2n]],[[15n,-3n],[15n,1n]]];
  const vec=(v)=>cols.reduce(([X,Y],[P,Q],i)=>[add(X,sc(P,v[i])),add(Y,sc(Q,v[i]))],[I(0),I(0)]);
  const norm=([X,Y])=>add(sq(X),sc(sq(Y),3));
  const inProductDisc=([R,S])=>sgnK([R-5184n,S])<=0&&sgnK([R-5184n,-S])<=0;
  const T=[];
  for(let a=-6;a<=6;a++)for(let b=-6;b<=6;b++)for(let c=-6;c<=6;c++)for(let d=-6;d<=6;d++){
    if(a*a+b*b+c*c+d*d>48||(a===0&&b===0&&c===0&&d===0))continue;
    const v=[a,b,c,d],p=vec(v),n=norm(p);
    if(inProductDisc(n))T.push({v,p,n});
  }
  const canon=(N,D)=>{
    let [d,zero]=mul(D,conj(D));assert(d!==0n&&zero===0n,"nonzero rational norm");
    let [a,b]=mul(N,conj(D));
    if(d<0n){a=-a;b=-b;d=-d}
    const g=gcd(gcd(a,b),d);return [a/g,b/g,d/g];
  };
  const scalarCentre=([l,m],[r,s],norm1,norm2)=>{
    const det=sub(mul(l,s),mul(r,m));
    if(z(det))return null;
    const D=sc(det,6), X=sc(sub(mul(norm1,s),mul(norm2,m)),3), Y=sub(mul(l,norm2),mul(r,norm1));
    if(!equal(norm([X,Y]),sc(sq(D),1296)))return null;
    return {X,Y,D,cx:canon(X,D),cy:canon(Y,D)};
  };
  const decode=([X0,X1,Xden],[Y0,Y1,Yden])=>{
    const L=Xden*Yden;x0: {
      const X=X0*Yden,U=X1*Yden,Y=Y0*Xden,V=Y1*Xden;
      if(-U%(3n*L)!==0n)break x0;
      const d=-U/(3n*L);
      if((V-d*L)%(2n*L)!==0n)break x0;
      const c=(V-d*L)/(2n*L);
      if((Y-15n*d*L)%(18n*L)!==0n)break x0;
      const b=(Y-15n*d*L)/(18n*L);
      if((X-(18n*b+30n*c+15n*d)*L)%(36n*L)!==0n)break x0;
      const a=(X-(18n*b+30n*c+15n*d)*L)/(36n*L);
      return [a,b,c,d].map(String);
    }return null;
  };
  const candidates=new Map();let noncollinear=0,successfulPairs=0;
  for(let i=0;i<T.length;i++)for(let j=i+1;j<T.length;j++){
    const a=T[i],b=T[j];
    if(z(sub(mul(a.p[0],b.p[1]),mul(a.p[1],b.p[0]))))continue;
    noncollinear++;
    const s=scalarCentre(a.p,b.p,a.n,b.n);
    if(s===null)continue;
    successfulPairs++;
    candidates.set(s.cx.join(",")+";"+s.cy.join(","),s);
  }
  const on={},off={},controls=[],cosets={};
  for(const s of candidates.values()){
    const verts=[[0,0,0,0]];
    for(const {v,p:[X,Y],n} of T){
      if(equal(sc(add(mul(X,s.X),sc(mul(Y,s.Y),3)),2),mul(n,s.D)))verts.push(v);
    }
    const coeff=decode(s.cx,s.cy), h=coeff===null?off:on;
    h[verts.length]=(h[verts.length]||0)+1;
    if(coeff===null){
      assert(verts.length===3,"off centre has exactly three");
      for(let i=0;i<verts.length;i++)for(let j=i+1;j<verts.length;j++)
        assert(equal(norm(vec(verts[i].map((k,t)=>k-verts[j][t]))),I(3888)),"off triangle is equilateral of side sqrt(3)");
      const triple=decode([3n*s.cx[0],3n*s.cx[1],s.cx[2]],[3n*s.cy[0],3n*s.cy[1],s.cy[2]]);
      assert(triple!==null,"off centre belongs to one-third lattice");
      const residue=triple.map(v=>String((BigInt(v)%3n+3n)%3n)).join("");
      cosets[residue]=(cosets[residue]||0)+1;
      if(controls.length<2)controls.push({X:s.cx.map(String),Y:s.cy.map(String),tripleCentreCoefficients:triple,anchors:verts});
    }
  }
  assert(T.length===204&&on[18]===18&&Object.keys(on).length===1&&Object.keys(off).join(",")==="3","expected census");
  assert(controls[0].anchors.length===3,"sharp positive control");
  assert(Object.keys(cosets).sort().join(",")==="0011,0022,1100,2200"&&Object.values(cosets).every(a=>a===3),"off coset census");
  return {shortDifferences:T.length,noncollinearPairs:noncollinear,successfulPairs,centres:candidates.size,on,off,cosets,firstOffControls:controls};
}
const output=audit();
if(typeof module!=="undefined")console.log(JSON.stringify(output,null,2));
return output;
