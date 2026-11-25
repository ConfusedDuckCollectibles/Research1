#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

using Seq=std::array<int,13>;
using Tri=std::array<int,3>;
constexpr int p=13;
int mod(int x){ x%=p;return x<0?x+p:x; }

Tri transform(Tri t,int a,int b){
    for(int&x:t)x=mod(a*x+b);
    std::sort(t.begin(),t.end());
    return t;
}
bool canonical(Tri t,bool scaling){
    for(int a=1;a<=(scaling?12:1);++a)
        for(int b=0;b<p;++b)
            if(transform(t,a,b)<t)return false;
    return true;
}
bool is_min_rotation(const Seq& s) {
    for(int r=1;r<p;++r) for(int x=0;x<p;++x)
        if(s[(x+r)%p]!=s[x]) {
            if(s[(x+r)%p]<s[x])return false;
            break;
        }
    return true;
}
using Buckets=std::array<std::vector<Seq>,p>;
std::pair<Buckets,int> build(bool scaling,bool rotation) {
    Buckets buckets; int supports=0;
    for(int a=0;a<p;++a)for(int b=a+1;b<p;++b)for(int c=b+1;c<p;++c){
        Tri t{a,b,c}; if(!canonical(t,scaling))continue;
        ++supports;
        const int allowed=(1<<a)|(1<<b)|(1<<c);
        Seq s{};
        auto gen=[&](auto&& self,int x,int s0,int s1,int s2,int used)->void {
            if(x==11) {
                // s0+u+v=0, s1+11u+12v=0.
                int v=mod(11*s0-s1), u=mod(-s0-v);
                if(!(allowed>>u&1)||!(allowed>>v&1))return;
                if((used|(1<<u)|(1<<v))!=allowed)return;
                s[11]=u;s[12]=v;
                if(rotation&&!is_min_rotation(s))return;
                buckets[mod(s2+121*u+144*v)].push_back(s);
                return;
            }
            for(int z:t){s[x]=z;self(self,x+1,s0+z,s1+x*z,s2+x*x*z,used|(1<<z));}
        };
        gen(gen,0,0,0,0,0);
    }
    return {std::move(buckets),supports};
}
int support_for_shift(const Seq& f,const Seq& g,int h) {
    unsigned m=0;
    for(int x=0;x<p;++x)m|=1u<<mod(f[x]-g[(x-h+p)%p]);
    return __builtin_popcount(m);
}
bool line_cap(const Seq& f) {
    for(int m=0;m<p;++m){
        Seq count{};
        for(int x=0;x<p;++x)if(++count[mod(f[x]-m*x)]>4)return false;
    }
    return true;
}
int square_sum(const Seq& f) {
    int n=0;for(int v:f)n+=v*v;return mod(n);
}
void print(const Seq& s) {for(int i:s)std::cout<<' '<<i;}
int main() {
    auto [b0,s0]=build(true,true);
    auto [b1,s1]=build(false,false);
    std::cout<<"supports "<<s0<<' '<<s1<<'\n';
    std::uint64_t pairs=0,second=0,third=0,n11=0,cap=0,valid=0;
    int best=0;
    for(int j=0;j<p;++j) {
        std::cout<<"bucket "<<j<<' '<<b0[j].size()<<' '<<b1[j].size()<<'\n';
        for(const Seq&a:b0[j]) for(const Seq&b:b1[j]){
            ++pairs;
            // For the first two new shifts, stop immediately when a fourth
            // distinct value appears, or at the end if it is not exactly 3.
            unsigned s=0;
            for(int x=0;x<p;++x){
                int y=(x+12)%p;
                s|=1u<<mod(b[x]+b[y]-a[y]);
                if(__builtin_popcount(s)>3)break;
            }
            if(__builtin_popcount(s)!=3)continue;
            ++second;
            s=0;
            for(int x=0;x<p;++x){
                int y=(x+12)%p,z=(x+11)%p;
                s|=1u<<mod(b[x]+b[y]-a[y]+b[z]-a[z]);
                if(__builtin_popcount(s)>3)break;
            }
            if(__builtin_popcount(s)!=3)continue;
            ++third;
            Seq f{},g{};
            for(int x=1;x<p;++x)g[x]=mod(g[x-1]+b[x]-a[x]);
            assert(mod(g[0]-g[12])==mod(b[0]-a[0]));
            for(int x=0;x<p;++x)f[x]=mod(g[x]+a[x]);
            int sf=0,sg=0,xf=0,xg=0;
            for(int x=0;x<p;++x){sf+=f[x];sg+=g[x];xf+=x*f[x];xg+=x*g[x];}
            assert(!mod(sf)&&!mod(sg)&&!mod(xf)&&!mod(xg));
            assert(support_for_shift(f,g,0)==3&&support_for_shift(f,g,1)==3);
            assert(support_for_shift(f,g,2)==3&&support_for_shift(f,g,3)==3);
            int good=0;for(int h=0;h<p;++h)good+=support_for_shift(f,g,h)<=3;
            if(good<11)continue;
            ++n11;
            if(!line_cap(f)||!line_cap(g))continue;
            ++cap;
            if(square_sum(f)||square_sum(g))continue;
            ++valid;
            if(good>best) {
                best=good;
                std::cout<<"witness "<<good<<" f";print(f);std::cout<<" g";print(g);std::cout<<'\n';
            }
        }
    }
    std::cout<<"counts "<<pairs<<' '<<second<<' '<<third<<' '<<n11<<' '<<cap<<' '<<valid<<'\n';
}
