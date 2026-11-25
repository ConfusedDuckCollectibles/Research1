// Exact independent rational replay, no ODE sampling or fit.
"use strict";
function gcd(a,b){if(a<0n)a=-a;while(b)[a,b]=[b,a%b];return a}
class Q{
 constructor(n,d=1n){n=BigInt(n);d=BigInt(d);if(!d)throw Error("zero");if(d<0n){n=-n;d=-d}const g=gcd(n,d);this.n=n/g;this.d=d/g}
 add(b){return new Q(this.n*b.d+b.n*this.d,this.d*b.d)}
 sub(b){return new Q(this.n*b.d-b.n*this.d,this.d*b.d)}
 mul(b){return new Q(this.n*b.n,this.d*b.d)}
 div(b){return new Q(this.n*b.d,this.d*b.n)}
 pow(k){k=BigInt(k);return new Q(this.n**k,this.d**k)}
 less(b){return this.n*b.d<b.n*this.d}
 str(){return this.n+"/"+this.d}
}
const R=(a,b=1)=>new Q(a,b), j=x=>x.str(), one=R(1);
const a=R(1,200),eta=R(1,1000000), cm=R(18,25), b=R(3,4), cp=R(39,50);
function enforce(v, label){if(!v)throw Error(label)}
function variant(sigma,nu,h,Ehi,phi,fhi,gu,rq,iu,err,upperDen,vacQ,sharp){
 const E0=one.add(nu), H=h.mul(b.pow(3)), chi=one.sub(R(2).mul(sigma).mul(Ehi)), ep=Ehi.add(phi);
 const gradient=sigma.mul(Ehi.add(R(3).mul(phi))).div(chi);
 const bulk=ep.mul(gu).div(R(2).mul(nu));
 const pre=ep.mul(gu).div(R(2));
 const ebar=Ehi.div(R(3));
 enforce(E0.add(fhi).less(Ehi)&&H.div(cm.sub(eta).pow(3)).less(fhi),"density bounds");
 enforce(gradient.less(gu)&&pre.less(phi)&&bulk.less(one),"TOV and causality");
 // F-1 = (2m/r - 4pi*r^2*(e-p))/D. Signed bounds are sharper
 // than adding its two magnitudes. Q follows from the separated formula.
 const fDist=R(3).mul(sigma).mul(Ehi).div(chi);
 const Qpoly=R(3).mul(sigma).mul(R(5).mul(Ehi).add(R(9).mul(phi)).add(R(12).mul(ebar))).div(chi);
 const Qdist=Qpoly.add(R(4).mul(gu.pow(2)));
 const Qup=R(3).mul(sigma).mul(R(5).mul(Ehi).add(R(9).mul(phi)).sub(R(4))).div(chi);
 enforce(fDist.less(one)&&Qdist.less(rq)&&Qup.less(vacQ),"tidal coefficient envelopes");
 enforce(R(2).mul(ebar).less(one.sub(phi))&&R(6).mul(nu).add(vacQ).less(R(4)),"regular surface barrier");
 const continuation=R(3).mul(fhi.add(nu.mul(a.pow(2)))).add(R(3).add(rq).mul(a.pow(2)).div(R(2)));
 enforce(continuation.less(one),"Riccati continuation");
 const invdiff=R(9).mul(nu.div(R(3)).add(pre));
 const finitePer=R(3).add(rq).div(R(2)).add(R(3).mul(nu)).add(iu.mul(fhi));
 enforce(invdiff.less(iu)&&finitePer.less(err),"finite radius error");
 // The one normalization bound simultaneously controls height and layer mean changes.
 const fDerivative=R(3).mul(fhi).div(cm.sub(eta));
 const smoothErr=R(3).add(R(9).mul(fhi)).add(R(18).div(cm.sub(eta)));
 enforce(fDerivative.less(one)&&smoothErr.less(R(40)),"smooth transition bounds");
 const maxDen=R(5).mul(E0.div(R(3)).add(phi)).add(R(2,3).mul(H.div(cm.pow(3))));
 enforce(maxDen.less(upperDen),"flat step denominator");
 const dL=b.pow(2).sub(cm.pow(2)),dR=cp.pow(2).sub(b.pow(2));
 const gap=R(16,5).mul(H).mul(dL).div(upperDen);
 const error=R(2).mul(err).mul(a.pow(2)).add(R(40).mul(eta));
 enforce(dL.less(dR)&&error.less(gap)&&j(gap)==="3969/5000000","finite IVT signs");
 const metric=sharp?R(5,3):R(25,18), x=R(3,2).mul(gu);
 const expHi=one.add(x).add(x.pow(2).div(R(2).mul(one.sub(x.div(R(3))))));
 if(sharp)enforce(R(6,7).pow(2).less(chi)&&expHi.less(R(10,7)),"sharpened metric");
 else enforce(R(9,10).pow(2).less(chi)&&!R(1,5).less(x),"baseline metric");
 const radial=metric.mul(sharp?one.div(chi):ep.mul(gu).div(sigma))
 .mul(R(4).add(gu)).mul(nu.div(R(3)).add(fhi));
 enforce(radial.less(one),"radial coercivity");
 const M=sigma.mul(one.add(R(2,5).mul(nu)).add(H.mul(a.pow(3))));
 enforce(!M.less(sigma)&&M.less(R(1,2)),"vacuum start");
 // At an interior step TOV gives -[1/p']/[1/g] = 1/Hout - 1/Hin.
 // Check its positive exact value at a valid sample pressure.
 const pout=phi.div(R(2)), outer=E0.sub(nu.mul(a).mul(b).pow(2)), entout=outer.add(pout), entin=entout.add(h);
 const atom=one.div(entout).sub(one.div(entin));
 enforce(R(0).less(atom)&&atom.str()===h.div(entout.mul(entin)).str(),"positive compliance atom");
 function T(c,q){return H.mul(c.pow(2)).div(R(5).mul(E0.div(R(3)).add(q)).add(R(2,3).mul(H.div(c.pow(3)))))}
 for(const q of [R(0),pout,phi]){
  const ys=[cm,b,cp].map(c=>{const z=T(c,q);return R(2).sub(R(5).mul(z).div(one.add(z)));});
  enforce(gap.less(ys[0].sub(ys[1]))&&gap.less(ys[1].sub(ys[2])),"exact flat signs");
 }
 return {compactness:j(M),causalUpper:j(bulk),pressureBootstrapSlack:j(phi.sub(pre)),
   yContinuationUpper:j(continuation),flatGapLower:j(gap),totalErrorUpper:j(error),
   strictWitnessMargin:j(gap.sub(error)),minimumRelativeJump:j(h.div(E0)),
   radialBound:j(radial),radialMargin:j(one.sub(radial))};
}
const base=variant(R(3,40),R(1,9),R(1,30),R(37,32),R(1,10),R(1,25),R(2,15),R(4),R(11,10),R(4),R(5,2),R(3,4),false);
const higher=variant(R(1,10),R(1,6),R(1,25),R(11,9),R(4,25),R(1,20),R(23,100),R(6),R(2),R(6),R(3),R(3,2),true);
enforce(base.radialBound==="54002/54675"&&higher.radialBound==="2679/2720","published numerical constants mismatch");
console.log(JSON.stringify({base,higher,checks:"PASS"},null,2));