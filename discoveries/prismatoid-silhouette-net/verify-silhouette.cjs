"use strict";
// Exact replay of the ten-vertex alternating-prismatoid controls.
// All interval endpoints are integer multiples of 1/S. No floating geometry.
const V=[[65,0],[63,16],[39,52],[16,63],[-25,60],[-56,33],[-65,0],[-52,-39],[-16,-63],[33,-56]].map(a=>a.map(BigInt));
const n=V.length,S=1000000000n;
const check=(x,m)=>{if(!x)throw Error(m);};
const sub=(u,v)=>u.map((x,i)=>x-v[i]),dot=(u,v)=>u.reduce((a,x,i)=>a+x*v[i],0n);
const cross=(u,v)=>u[0]*v[1]-u[1]*v[0],n3=(u,v)=>[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
const abs=x=>x<0n?-x:x, sq=x=>x*x;
const min=a=>a.reduce((x,y)=>x<y?x:y),max=a=>a.reduce((x,y)=>x>y?x:y);
const I=x=>[x*S,x*S],plus=(x,y)=>[x[0]+y[0],x[1]+y[1]],minus=(x,y)=>[x[0]-y[1],x[1]-y[0]];
const floor=(x,y)=>x>=0n?x/y:-((-x+y-1n)/y),ceil=(x,y)=>-floor(-x,y);
function mul(a,b){
 const t=a.flatMap(x=>b.map(y=>x*y));
 return [floor(min(t),S),ceil(max(t),S)];
}
const scale=(t,n)=>n>=0n?[n*t[0],n*t[1]]:[n*t[1],n*t[0]], idiv=(t,d)=>[floor(t[0],d),ceil(t[1],d)];
function sqrt(n){
 if(n<0n)throw Error("negative radicand");
 if(n<2n)return n;
 let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));
 while(true){const y=(x+n/x)/2n;if(y>=x)return x;x=y;}
}
function rsqrt(a,b,den=1n){
 const u=sqrt(a*S*S),v=sqrt(b*S*S);
 return [floor(u,den),ceil(v+(v*v===b*S*S?0n:1n),den)];
}
const at=(i)=>V[(i+n)%n];
const ler=(w,i,j)=>cross(sub(at(j),at(i)),sub(w,at(i)));
for(let i=0;i<n;i++){
 check(dot(V[i],V[i])===4225n,"not on claimed circle");
 for(let j=0;j<n;j++)if(j!==i&&j!==(i+1)%n)check(ler(V[j],i,(i+1)%n)>0n,"not strictly CCW convex");
}
let crossings=0;
for(let i=0;i<n;i+=2)for(let j=1;j<n;j+=2)
 crossings+=Number(ler(at(j),i,i+2)*ler(at(j+2),i,i+2)<0n&&ler(at(i),j,j+2)*ler(at(i+2),j,j+2)<0n);
check(crossings===10,"base crossing count");
let obtuseAfter30=false;
for(let j=0;j<n;j++){
 const l=at(j-1), r=at(j+1), w=at(j), e=sub(r,l),q=sub(w,l),D=dot(e,e),d=dot(q,e);
 check(d>0n&&d<D,"ear foot not in base edge");
 check(cross(e,q)<0n,"ear not exterior");
 if(dot(sub(l,w),sub(r,w))+900n<0n)obtuseAfter30=true;
}
check(obtuseAfter30,"missing obtuse side triangle");
function caps(k,h0,h1,den){
 return V.map((w,j)=>{
   if(j%2===k)return w.map(I);
   const l=at(j-1),r=at(j+1),e=sub(r,l),out=[e[1],-e[0]],q=sub(w,l),D=dot(e,e),foot=dot(q,e),perp=dot(q,out);
   const R=rsqrt(sq(perp*den)+h0*h0*D,sq(perp*den)+h1*h1*D,den);
   return l.map((x,i)=>plus(I(x),idiv(plus(I(foot*e[i]),scale(R,out[i])),D)));
 });
}
const vsub=(u,v)=>u.map((x,i)=>minus(x,v[i]));
const icross=(u,v)=>minus(mul(u[0],v[1]),mul(u[1],v[0]));
function seam(h0,h1,den,i){
 const j=(i+1)%n;let res=null;
 for(const k of [0,1]){
    const C=caps(k,h0,h1,den), e=vsub(C[j],C[i]);
    for(let u=0;u<n;u++)if(u!==i&&u!==j){
      const lower=icross(e,vsub(C[u],C[i]))[0];
      res=res===null?lower:min([lower,res]);
    }
 }
 return res;
}
let margin=null;const den=64n;
for(let h=0n;h<30n*den;h++){
 const q=seam(h,h+1n,den,8);
 check(q>0n,"unverified height interval "+h);
 margin=margin===null?q:min([margin,q]);
}
// For an alternating prismatoid the cap corner on its own base side has
// angles of three sectors: the angle of the base, and the two mixed ear angles.
// These angles are in (0,pi), (0,pi/2), (0,pi/2). Thus their total is in (0,2pi),
// and the interval sign of the imaginary part of the product tests <=pi.
function angleSinSign(i,h){
 const u=sub(at(i-2),at(i)),v=sub(at(i+2),at(i)),w=sub(at(i-1),at(i)),z=sub(at(i+1),at(i));
 const d0=dot(u,v),d1=dot(u,w),d2=dot(v,z),s0=abs(cross(u,v));
 check(d1>0n&&d2>0n,"nonacute base angles");
 const t1=sq(cross(u,w))+sq(h)*dot(u,u),t2=sq(cross(v,z))+sq(h)*dot(v,v);
 const R1=rsqrt(t1,t1),R2=rsqrt(t2,t2);
 return minus(plus(plus(scale(R2,d0*d1),scale(R1,d0*d2)),I(s0*d1*d2)),scale(mul(R1,R2),s0));
}
const signs=Array.from({length:n},(_,i)=>[angleSinSign(i,13n),angleSinSign(i,14n)]);
check(signs.every(x=>x[0][0]>0n),"not all quasigeodesic at h=13");
const nonquasi14=signs.flatMap((x,i)=>x[1][1]<0n?[i]:[]);
check(nonquasi14.join()==="1","strict reflex control at h=14");
const F=[Array.from({length:n/2},(_,i)=>2*i),Array.from({length:n/2},(_,i)=>2*i+1)];
for(let i=0;i<n;i++)F.push([(i+n-1)%n,i,(i+1)%n]);
for(const h of [13n,14n,30n]){
  const P=V.map((w,i)=>[...w,i%2?h:0n]);
  for(const f of F){
    const [a,b,c]=f.map(i=>P[i]),N=n3(sub(b,a),sub(c,a));
    check(N[2]!==0n,"vertical/degenerate facet");
    const s=P.map(p=>dot(N,sub(p,a)));
    check(f.every(i=>s[i]===0n)&&s.flatMap((t,i)=>t===0n?[i]:[]).length===f.length
      &&(s.every(x=>x>=0n)||s.every(x=>x<=0n)),"nonfacetal label");
  }
}
const out={crossings,certifiedHeightIntervals:Number(30n*den),seam:[8,9],
 minSignedCrossLowerNumerator:margin.toString(),scale:S.toString(),nonquasi14,obtuseAfter30};
console.log(JSON.stringify(out));
