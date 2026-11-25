// Independent generation: assign the first two symbols by disjoint bit masks;
// recover the two graphs directly for every pair. This does not use the
// two-last-symbol recursion or the difference-array D2/D3 formulas.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;
constexpr int p=13,full=(1<<p)-1;
using Vec=array<unsigned char,p>;
using Tri=array<int,3>;
using Bins=array<vector<Vec>,p>;
int r(int x){x%=p;return x<0?x+p:x;}
int cnt[1<<p],s1[1<<p],s2[1<<p];
Tri transform(Tri q,int l,int d){
    for(int&x:q)x=r(l*x+d);
    sort(q.begin(),q.end());
    return q;
}
bool least(const Vec& v) {
    Vec w=v;
    for(int j=1;j<p;++j){
        rotate(w.begin(),w.begin()+1,w.end());
        if(w<v)return false;
    }
    return true;
}
Bins arrays(bool affine) {
    set<Tri> reps;
    for(int b=1;b<p;++b)for(int c=b+1;c<p;++c){
        Tri q{0,b,c},best=q;
        for(int m=1;m<=(affine?12:1);++m)
            for(int d=0;d<p;++d) best=min(best,transform(q,m,d));
        reps.insert(best);
    }
    assert(reps.size()==(affine?3:22));
    Bins result;
    for(auto [a,b,c]:reps)for(int i=1;i<full;++i){
        int rest=full^i;
        for(int j=(rest-1)&rest;j;j=(j-1)&rest){
            int k=rest^j;
            if((a*cnt[i]+b*cnt[j]+c*cnt[k])%p)continue;
            if((a*s1[i]+b*s1[j]+c*s1[k])%p)continue;
            Vec v;
            for(int x=0;x<p;++x)v[x]=(i>>x&1)?a:(j>>x&1)?b:c;
            if(affine&&!least(v))continue;
            result[(a*s2[i]+b*s2[j]+c*s2[k])%p].push_back(v);
        }
    }
    return result;
}
int image(const Vec&f,const Vec&g,int h){
    unsigned seen=0;
    int y=(p-h)%p;
    for(int x=0;x<p;++x, y=(y+1)%p){
        seen|=1u<<r(f[x]-g[y]);
        if(__builtin_popcount(seen)>3)return 4;
    }
    return __builtin_popcount(seen);
}
int peak(const Vec&f) {
    int high=1;
    for(int k=0;k<p;++k){
        int hist[p]={};
        for(int x=0;x<p;++x)high=max(high,++hist[r(f[x]-k*x)]);
    }
    return high;
}
int main(){
    for(int m=1;m<=full;++m) {
        int x=__builtin_ctz(m), n=m&(m-1);
        cnt[m]=cnt[n]+1;s1[m]=s1[n]+x;s2[m]=s2[n]+x*x;
    }
    Bins a=arrays(true),b=arrays(false);
    uint64_t all=0,second=0,third=0,good=0,cap=0;
    map<array<int,3>,int> types;
    for(int j=0;j<p;++j) {
        cout<<"bucket "<<j<<' '<<a[j].size()<<' '<<b[j].size()<<'\n';
        for(const Vec& x:a[j])for(const Vec& y:b[j]){
            ++all;
            Vec F{},G{};
            F[0]=x[0];
            for(int i=1;i<p;++i){
                G[i]=r(G[i-1]+y[i]-x[i]);F[i]=r(G[i]+x[i]);
            }
            if(image(F,G,2)!=3)continue;
            ++second;
            if(image(F,G,3)!=3)continue;
            ++third;
            int q=0;
            for(int h=0;h<p;++h)q+=image(F,G,h)<=3;
            if(q>=11){
                ++good;int f=peak(F),g=peak(G);cap+=max(f,g)<=4;
                types[{q,min(f,g),max(f,g)}]++;
            }
        }
    }
    cout<<"counts "<<all<<' '<<second<<' '<<third<<' '<<good<<' '<<cap<<'\n';
    for(auto [t,n]:types)cout<<"type "<<t[0]<<' '<<t[1]<<' '<<t[2]<<' '<<n<<'\n';
    assert(all==49895820&&second==70806&&third==11479&&good==5418&&cap==0);
}
