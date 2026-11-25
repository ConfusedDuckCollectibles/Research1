// Exact independent checks for dwarf aperture sample and compensated annulus.
// No input files. Only optional output decimals use floating point.
function q(a,b=1n){
  a=BigInt(a);b=BigInt(b);if(b===0n)throw Error("denominator");
  if(b<0n){a=-a;b=-b;}
  let u=a<0n?-a:a,v=b;
  while(v!==0n)[u,v]=[v,u%v];
  return [a/u,b/u];
}
const add=(x,y)=>q(x[0]*y[1]+x[1]*y[0],x[1]*y[1]);
const sub=(x,y)=>add(x,[-y[0],y[1]]);
const mul=(x,y)=>q(x[0]*y[0],x[1]*y[1]);
const div=(x,y)=>q(x[0]*y[1],x[1]*y[0]);
const power=(x,n)=>[x[0]**BigInt(n),x[1]**BigInt(n)];
const lt=(x,y)=>x[0]*y[1]<y[0]*x[1];
const same=(x,y)=>!lt(x,y)&&!lt(y,x);
const sum=xs=>xs.reduce(add,q(0));
const val=x=>Number(x[0])/Number(x[1]);
const must=(a,s)=>{if(!a)throw Error(s);};
const one=q(1),zr=q(0);
function logRatio(x,n){ // log((1+x)/(1-x)), positive x<1
  let m=zr;
  for(let j=0;j<n;j++)m=add(m,mul(q(2,2*j+1),power(x,2*j+1)));
  return [m,add(m,div(mul(q(2,2*n+1),power(x,2*n+1)),sub(one,power(x,2))))];
}
const L4=logRatio(q(3,5),42);
function population(k,N){
  const f=10n**16n, d=(1n+BigInt(k)**2n)**3n, square=4n*f*f;
  let a=0n,b=f;
  while(b-a>1n){const m=(a+b)>>1n;if(m*m*d<=square)a=m;else b=m;}
  const p=[q(a,f),q(b,f)];
  const logs=p.map(x=>logRatio(div(x,sub(q(2),x)),5));
  must(lt(mul(q(N-1),logs[1][1]),L4[0]),"N−1 still below 80% power bound");
  must(lt(L4[1],mul(q(N),logs[0][0])),"N reaches this necessary condition");
  return {p,N};
}
const five=population(20,5566),ten=population(10,703);
function bernoulliLeak(p,n){
  // Inclusion/exclusion truncations of the probability of the union, exact.
  let coeff=1n,lo=zr,hi=null,s=zr;
  for(let j=1;j<=10;j++){
    coeff=coeff*BigInt(n+1-j)/BigInt(j);
    const v=mul(q(coeff),power(p,j));
    s=j%2===1?add(s,v):sub(s,v);if(j===9)hi=s;if(j===10)lo=s;
  }return [lo,hi];
}
const p119=[add(q(1,20),bernoulliLeak(five.p[0],119)[0]),
            add(q(1,20),bernoulliLeak(five.p[1],119)[1])];
must(lt(q(792074,10000000),p119[0])&&lt(p119[1],q(792075,10000000)),"119 upper-bound expression");
const I=L4.map(x=>sub(mul(q(2),x),q(11,4))); // I=2 log 4-11/4
must(lt(q(22588722239781n,1000000000000000n),I[0])&&lt(I[1],q(22588722239782n,1000000000000000n)),"positive normalization");
const poly=[4,-12,13,-6,1];
const integrate=(pow)=>sum(poly.map((c,j)=>q(c*(2**(pow+j+1)-1),pow+j+1)));
must(same(integrate(0),q(1,30))&&same(integrate(1),q(1,20)),"annulus exact moments");
const delta=q(24,25), angleRoot=q(979,1000);
must(lt(power(angleRoot,2),delta),"certified lower sqrt(24/25)");
const partialCosh3=sum([q(1),q(9,2),q(81,24),q(729,720),q(6561,40320)]);
must(lt(q(10),partialCosh3),"arcosh(10)<3");
const factor=div(q(1,2000),I[0]);
const lower=sub(one,mul(factor,sum([div(q(2),angleRoot),q(6),div(q(2),mul(delta,angleRoot))])));
const upper=mul(factor,add(div(q(2),angleRoot),q(2)));
must(lt(upper,q(9,100))&&lt(q(77,100),lower),"finite annular synthetic separation");
// Independent positive-bracket certificate, without the candidate's degree-five expansion.
// For r<=1/4, 2/[r(1+6r^2)]>=64/11>2; for r>=1/4,
// (7r/(1+r^2)-2/(1+r)) has numerator 5r²+7r−2, at least 1/16.
must(lt(q(2),q(64,11))&&same(add(add(q(5,16),q(7,4)),q(-2)),q(1,16)),"Hernquist bracket split");
const A=q(1,4);
const H1=add(mul(q(5),A),q(2));
const H2=add(sum([q(5),mul(q(35),power(A,2)),mul(q(20),A)]),q(6));
must(lt(H1,div(q(2),A))&&lt(div(H2,q(2)),div(q(2),power(A,2))),"Hernquist physical remainder");
console.log(JSON.stringify({checks:"PASS",necessaryN80:[five.N,ten.N],p05:five.p.map(val),bound119:p119.map(val),
  I:"0.022588722239781 < I < 0.022588722239782",coreUpper:val(upper),cuspLower:val(lower),HernquistH1:H1.map(String),HernquistH2:H2.map(String)}));

