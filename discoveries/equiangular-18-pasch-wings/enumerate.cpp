// No floating point. Enumerate all 2^18 sign/inner-product assignments.
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
using I=long long;
using Q=__int128_t;
using V=std::array<I,19>;
static Q dot(const V&a,const V&b){
    Q x=0;
    for(int k=0;k<15;k+=5){
        I sa=0,sb=0;
        for(int j=k;j<k+5;++j){x+=Q(a[j])*b[j];sa+=a[j];sb+=b[j];}
        x+=Q(sa)*sb;
    }
    for(int j=15;j<19;++j)x+=2*Q(a[j])*b[j];
    return x;
}
int main(){
    I d; assert(bool(std::cin>>d));
    std::array<V,19> A;
    std::array<V,21> W;
    for(auto& v:A) for(auto& y:v) assert(bool(std::cin>>y));
    for(auto& v:W) for(auto& y:v) assert(bool(std::cin>>y));
    assert(d>0 && d<(1LL<<45) && (d%6)==0);
    for(auto const& a:A)for(I x:a)assert(-((1LL<<44))<x&&x<(1LL<<44));
    Q D=Q(d)*d;
    V v=A[0];
    for(unsigned i=0;i<(1u<<18);++i){
        if(i){
            int k=__builtin_ctz(i);
            bool set=((i^(i>>1))>>k)&1u;
            for(unsigned j=0;j<19;++j)v[j]+=set ? A[k+1][j] : -A[k+1][j];
        }
        if(dot(v,v)!=3*D)continue;
        V neg=v;
        for(I& y:neg)y=-y;
        neg[18]+=d;
        assert(v[18]==d/2);
        if(v>=neg)continue;
        bool ok=true;
        for(auto const& w:W){Q t=dot(v,w);if(t!=0 && t!=6*Q(d)){ok=false;break;}}
        if(!ok)continue;
        for(I y:v)std::cout<<y<<' ';
        std::cout<<'\n';
    }
}
