// Independent decision checker: no 32 words in [6]^5 with <=1 agreement.
// No upstream files. Generate each full class row-wise, quotient by adjacent
// group generators, and test the five exact missing-class profiles for 32.
#if defined(NDEBUG)
#error Enable assertions
#endif
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <unordered_set>
#include <vector>
using U=uint64_t;
using W=__uint128_t;
using T=std::array<unsigned,4>;
using Class=std::array<T,6>;

std::array<std::vector<T>,6> tails;
std::vector<T> all;
std::unordered_set<U> pending;
std::array<unsigned,4> used{};
Class growing;

U key(Class const& A) {
    U k=0;
    for(auto const& t:A) for(unsigned j=1;j<4;++j) {
        assert(t[0]<6&&t[j]<6);
        k|=U(t[j]) << (3*(3*t[0]+j-1));
    }
    return k;
}
Class decode(U k) {
    Class A;
    for(unsigned i=0;i<6;++i) {
        A[i][0]=i;
        for(unsigned j=1;j<4;++j) A[i][j]=(k>>(3*(3*i+j-1)))&7;
    }
    return A;
}
void generate(int row) {
    if(row==6) { assert(pending.insert(key(growing)).second); return; }
    for(auto const& t:tails[row]) {
        if((used[1]>>t[1]&1)||(used[2]>>t[2]&1)||(used[3]>>t[3]&1)) continue;
        for(int j=1;j<4;++j) used[j]^=1u<<t[j];
        growing[row]=t; generate(row+1);
        for(int j=1;j<4;++j) used[j]^=1u<<t[j];
    }
}
std::vector<U> quotient() {
    std::vector<U> reps;
    while(!pending.empty()) {
        std::vector<U> orbit{*pending.begin()};
        pending.erase(orbit.front());
        for(size_t k=0;k<orbit.size();++k) {
            Class A=decode(orbit[k]);
            for(unsigned s=0;s<5;++s) {
                Class B=A;
                for(auto& t:B) for(auto& j:t) {
                    if(j==s) j=s+1;
                    else if(j==s+1) j=s;
                }
                U y=key(B);
                if(pending.erase(y)) orbit.push_back(y);
            }
            for(unsigned s=0;s<3;++s) {
                Class B=A;
                for(auto& t:B) std::swap(t[s],t[s+1]);
                U y=key(B);
                if(pending.erase(y)) orbit.push_back(y);
            }
        }
        reps.push_back(*std::min_element(orbit.begin(),orbit.end()));
    }
    std::sort(reps.begin(),reps.end());
    return reps;
}

std::array<W,128> same,cross;
std::array<int,4> sizes;
std::array<U,5> nodes{};
int count(W bits) {
    return std::popcount((U)bits)+std::popcount((U)(bits>>64));
}
int low(W bits) {
    assert(bits);
    U l=(U)bits;
    return l ? std::countr_zero(l) : 64+std::countr_zero((U)(bits>>64));
}
W above(int index) {return (index==127)?W(0):(~W(0)<<(index+1));}
unsigned agree(T const& x,T const& y) {
    unsigned v=0; for(int i=0;i<4;++i) v+=(x[i]==y[i]); return v;
}

bool start_class(int depth,W cross_pool,int prior_min,int profile);
bool fill_class(int depth,int need,W pool,W cross_pool,int this_min,int profile) {
    ++nodes[profile];
    if(need==0) return start_class(depth+1,cross_pool,this_min,profile);
    while(count(pool)>=need) {
        int i=low(pool);pool&=pool-1;
        if(fill_class(depth,need-1,pool&same[i],cross_pool&cross[i],
                      (this_min<0)?i:this_min,profile)) return true;
    }
    return false;
}
bool start_class(int depth,W cross_pool,int prior_min,int profile) {
    if(depth==4) return true;
    if(count(cross_pool)<std::accumulate(sizes.begin()+depth,sizes.end(),0)) return false;
    W pool=cross_pool;
    if(depth>0&&sizes[depth]==sizes[depth-1]) pool&=above(prior_min);
    return fill_class(depth,sizes[depth],pool,cross_pool,-1,profile);
}

int main() {
    for(unsigned a=0;a<6;++a) for(unsigned b=0;b<6;++b)
    for(unsigned c=0;c<6;++c) for(unsigned d=0;d<6;++d) {
        if(std::popcount((1u<<a)|(1u<<b)|(1u<<c)|(1u<<d))!=4) continue;
        T t{a,b,c,d}; tails[a].push_back(t); all.push_back(t);
    }
    assert(all.size()==360);
    generate(0);
    assert(pending.size()==393120);
    std::cout<<"Second full classes "<<pending.size()<<'\n';
    auto reps=quotient();
    assert(reps.size()==56);
    std::cout<<"Symmetry components "<<reps.size()<<'\n';
    std::set<std::array<int,4>> allprof;
    for(int a=0;a<=6;++a) for(int b=0;b<=a;++b)
    for(int c=0;c<=b;++c) for(int d=0;d<=c;++d)
        if(a+b+c+d==20) allprof.insert({a,b,c,d});
    assert(allprof.size()==5);
    unsigned min_pool=360,max_pool=0;
    for(U rep:reps) {
        Class A=decode(rep);
        std::vector<T> pool;
        for(auto const& t:all) {
            if(std::all_of(A.begin(),A.end(),[&](T const& s){return agree(s,t)<=1;}))
                pool.push_back(t);
        }
        unsigned n=pool.size(); assert(n<128);
        min_pool=std::min(min_pool,n);max_pool=std::max(max_pool,n);
        same.fill(0);cross.fill(0);
        for(unsigned i=0;i<n;++i) for(unsigned j=0;j<n;++j) {
            unsigned x=agree(pool[i],pool[j]);
            if(x<=1) cross[i]|=W(1)<<j;
            if(x==0) same[i]|=W(1)<<j;
        }
        int profile=0;
        for(auto const& p:allprof) {
            sizes=p;
            if(start_class(0,(W(1)<<n)-1,-1,profile++)) {
                std::cerr<<"Unexpected 32-word completion at class "<<rep<<'\n';
                return 1;
            }
        }
    }
    assert(min_pool==72&&max_pool==96);
    int i=0;
    for(auto const& p:allprof) {
        std::cout<<"Profiles"; for(int x:p) std::cout<<' '<<x;
        std::cout<<" : "<<nodes[i++]<<" states, infeasible\n";
    }
    std::cout<<"PASS: no 32-word packing; pa(5,6)<=31\n";
}
