// Exact interval-intersection generator; independent checker is verify.py.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <numeric>
#include <vector>
using i64=int64_t;
using i128=__int128_t;
constexpr int n=16;
struct Q { i64 a,b; };
bool less(Q x,Q y){return (i128)x.a*y.b < (i128)y.a*x.b;}
Q reduce(Q x){auto g=std::gcd(x.a,x.b);return {x.a/g,x.b/g};}
Q subtract(Q x,Q y){return {x.a*y.b-y.a*x.b,x.b*y.b};}
std::ostream& operator<<(std::ostream& o,Q q){q=reduce(q);return o<<q.a<<'/'<<q.b;}
using Window = std::array<Q,2>;
using Windows = std::vector<Window>;
// All intervals here are open. Earlier periods end below the given lower endpoint.
Windows install(const Windows& good,int w){
    Windows next;
    for(const auto& J:good){
        int j=(i128(w)*J[0].a)/J[0].b;
        for(;j<w;j++){
            Q left{n*j+1,n*w},right{n*j+n-1,n*w};
            if(!less(left,J[1])) break;
            Q lo=less(left,J[0])?J[0]:left;
            Q hi=less(J[1],right)?J[1]:right;
            if(less(lo,hi)) next.push_back({lo,hi});
        }
    }
    return next;
}
Window longest(const Windows& good){
    assert(!good.empty());
    Window best=good[0];
    for(const auto& J:good)
        if(less(subtract(best[1],best[0]), subtract(J[1],J[0]))) best=J;
    return best;
}
struct Stop {};
struct Audit{
    int K,start;
    std::array<i64,6> visited{};
    std::array<int,5> caps{};
    std::ostream* certificate;
    i64 limit, total=0;
    int roots=0, completed_roots=0, last_max=0;
    std::vector<int> removed;
    Q min_gap{1,1};
    explicit Audit(int k,i64 b,std::ostream *c,int m):K(k),start(m),certificate(c),limit(b){}
    void stop(const char* reason,const std::vector<int>& S,int r) {
        std::cerr<<"STOP "<<reason<<" remaining="<<r<<" removed";
        for(int x:removed) std::cerr<<" "<<x;
        std::cerr<<" installed";
        for(int x:S) std::cerr<<" "<<x;
        std::cerr<<"\n";
        throw Stop{};
    }
    void descend(std::vector<int>& S,const Windows& good,int lo,int r) {
        int depth=K-r;
        if(total>=limit) stop("visit cap",S,r);
        ++total; ++visited[depth];
        if(good.empty()) stop("empty strict-good set",S,r);
        auto J=longest(good);
        Q L=subtract(J[1],J[0]);
        if (!r){
            Q w=reduce({J[0].a*J[1].b+J[1].a*J[0].b,2*J[0].b*J[1].b});
            Q gap{1,1};
            for (int v:S){
                i64 m=(i128(v)*w.a)%w.b;
                Q dist{std::min(m,w.b-m),w.b};
                if(!less({1,n},dist)) stop("midpoint failed",S,r);
                if (less(dist,gap)) gap=dist;
            }
            if (less(gap,min_gap)) min_gap=gap;
            last_max=std::max(last_max,S.back());
            if(certificate) *certificate<<"W "<<w.a<<" "<<w.b<<"\n";
            return;
        }
        if(certificate) *certificate<<"I "<<J[0].a<<" "<<J[0].b<<" "<<J[1].a<<" "<<J[1].b<<"\n";
        i64 num=2*r*(n-2)*L.b,den=n*(n-2*r)*L.a;
        assert(num>0&&den>0);
        i64 cap=num/den;
        // The arithmetic above is exact for speeds up to this defensive limit.
        // Crossing it is UNRESOLVED, never reported as exclusion.
        if(cap>1000000) stop("representation cap",S,r);
        caps[depth]=std::max<i64>(caps[depth],cap);
        assert((i128)cap*den<=num && (i128)(cap+1)*den>num);
        for(int w=lo;w<=cap;w++){
            S.push_back(w);
            descend(S,install(good,w),w+1,r-1);
            S.pop_back();
        }
    }
    bool run(){
        bool pass=false;
        try{
            if(certificate) *certificate<<n<<" "<<K<<"\n";
            for(int mask=start;mask<(1<<(n-1));mask++){
                if (__builtin_popcount((unsigned)mask)!=K) continue;
                std::vector<int> S; removed.clear();
                for(int i=1;i<n;i++)
                    if((mask&1<<(i-1))==0) S.push_back(i);
                    else removed.push_back(i);
                ++roots;
                if(certificate) *certificate<<"R "<<mask<<"\n";
                Windows good{Window{Q{0,1},Q{1,1}}};
                for(int v:S) good=install(good,v);
                descend(S,good,n,K);
                ++completed_roots;
                if(K==5 && completed_roots%128==0)
                    std::cerr<<"PROGRESS completed_roots="<<completed_roots<<" total_visits="<<total<<"\n";
            }
            pass=true;
        }catch(const Stop&){}
        std::cout<<"k="<<K<<" start-mask="<<start<<" status="<<(pass?(start?"COMPLETE SUFFIX":"COMPLETE"):"INCOMPLETE")
                 <<" roots-started="<<roots<<" roots-completed="<<completed_roots<<" visits=[";
        for(int i=0;i<=K;i++) std::cout<<(i?",":"")<<visited[i];
        std::cout<<"] caps=[";
        for(int i=0;i<K;i++) std::cout<<(i?",":"")<<caps[i];
        std::cout<<"] min-witness-gap="<<min_gap<<" largest-leaf-speed="<<last_max<<"\n";
        return pass;
    }
};
int main(int argc,char**argv){
    if(argc<3||argc>5){ std::cerr<<"usage: sweep K MAX-VISITS [certificate.txt [start-mask]]\n"; return 2; }
    const int k=std::stoi(argv[1]); const auto max_visits=std::stoll(argv[2]);
    const int start=argc==5?std::stoi(argv[4]):0;
    if(k<2||k>5||max_visits<=0||start<0||start>=(1<<(n-1))) return 2;
    std::ofstream file;
    if(argc>=4){file.open(argv[3]); if(!file) return 2;}
    bool ok=Audit(k,max_visits,argc>=4?&file:nullptr,start).run();
    if(argc>=4){file.flush(); if(!file) return 2;}
    return ok?0:3;
}
