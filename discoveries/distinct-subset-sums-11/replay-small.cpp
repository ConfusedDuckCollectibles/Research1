// Exact numeric local obstruction in dimensions 9 and 10; C++17, no inputs.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
using namespace std;
void need(bool t,const char* m){if(!t){cerr<<m<<endl;exit(1);}}
#ifndef DIM
#define DIM 10
#endif
constexpr int n=DIM, Q=1<<n, H=n==10?274:148, W=n==10?175:49, K=W+63,
    tail=n==10?5:4, TM=(1<<tail)-1, T=K+1+W+1, sumD=n==10?102:26;
static_assert(n==9 || n==10,"compile with DIM=9 or DIM=10");
constexpr int lo=-tail*H,hi=H+K+sumD,domain=hi-lo+1;
struct Template{
    vector<int> canon,ranks,cuts,U;
    array<uint32_t,H+W+1> bans{};
    Template():ranks(Q) {
        for(int k=0;k<=n;++k) gen(k,0,0);
        for(int i=0;i<Q;++i)ranks[canon[i]]=i+1;
        for(int i=0;i<Q;++i){
            int s=canon[i]; if(s&(1<<(n-1)))continue;
            bool good=true;
            for(int b=(n+1)%2;b<n-1;b+=2){int r=(s>>b)&3;if(r!=1&&r!=2)good=false;}
            if(good){ need(ranks[s|(1<<(n-1))]==i+1+H,"shift"); cuts.push_back(i+1); }
        }
        need(cuts.size()==(1u<<(n/2))&&cuts.back()-cuts.front()==W,"span");
        for(int k=0;k<tail-1;++k){
            int best=10000;
            for(int p:cuts)for(int q:cuts)if(p<q)
                for(int r=p;r<q;++r) if(!(canon[r]>>(n-tail+k)))best=min(best,q-p-1);
            need(best==array<int,4>{76,20,5,1}[k+5-tail],"deficits");
        }
        for(int slab=0;slab<=1;++slab)
            for(unsigned i=1;i<cuts.size();++i){
                int a=cuts[i-1]+slab*H,b=cuts[i]+slab*H;
                for(int k=0;k<=tail;++k){
                    int f=TM^((1<<k)-1),direct=0,mirror=0;
                    for(int p=a;p<b;++p){
                        direct+=!((canon[p]>>(n-tail))&f);
                        mirror+=((canon[p-1]>>(n-tail))&f)==f;
                    }
                    for(int m=0;m<=TM;++m){
                        bool sub=!(m&f),sup=(m&f)==f;
                        if((direct==0&&sub)||(direct==b-a&&!sub))
                            for(int r=a+1;r<=b;++r)bans[r-cuts.front()]|=uint32_t(1)<<m;
                        if((mirror==0&&sup)||(mirror==b-a&&!sup))
                            for(int r=a;r<b;++r)bans[r-cuts.front()]|=uint32_t(1)<<m;
                    }
                }
            }
        for(int q=0;q<=K;++q)U.push_back(q);
        for(int q=0;q<=W;++q)U.push_back(H+q);
    }
    void gen(int k,int i,int s){
        if(!k){canon.push_back(s);return;}
        for(;i<=n-k;++i)gen(k-1,i+1,s|(1<<i));
    }
    int idx(int q)const{
        if(q>=0&&q<=K)return q;
        if(q>=H&&q<=H+W)return K+1+q-H;
        return -1;
    }
};
struct Search{
    vector<int> E,X,D;
    vector<vector<int>> R,C;
    array<int,domain> map;
    bool duplicate=false;
    int branches=0;
    Search(const Template& t,const array<int,5>& d):R(T) {
        for(int m=0;m<=TM;++m){
            int e=0;
            for(int i=0;i<tail;++i)if(m&(1<<i))e+=H-d[i];
            E.push_back(e);
        }
        auto sorted=E;sort(sorted.begin(),sorted.end());
        if(adjacent_find(sorted.begin(),sorted.end())!=sorted.end()){duplicate=true;return;}
        for(int a:E)for(int b:E)if(a!=b) D.push_back(a-b);
        sort(D.begin(),D.end());D.erase(unique(D.begin(),D.end()),D.end());
        map.fill(-1);
        for(int r:t.U)for(int e:E){need(r-e>=lo&&r-e<=hi,"domain");map[r-e-lo]=-2;}
        for(int x=lo;x<=hi;++x)if(map[x-lo]==-2){
            vector<int> rows;
            bool bad=false;
            for(int m=0;m<=TM;++m){
                int y=x+E[m];
                if(y>=0&&y<=H+W&&(t.bans[y]&(uint32_t(1)<<m))){bad=true;break;}
                int i=t.idx(y);
                if(i>=0)rows.push_back(i);
            }
            if(bad){map[x-lo]=-1;continue;}
            need(!rows.empty(),"candidate hits a rank");
            int j=int(X.size());
            for(int i:rows)R[i].push_back(j);
            map[x-lo]=j;X.push_back(x);C.push_back(move(rows));
        }
    }
    void choose(int j,vector<unsigned char>& a,vector<unsigned char>& b)const{
        need(a[j],"viable");a[j]=0;
        for(int i:C[j]){need(!b[i],"collision");b[i]=1;}
        for(int d:D){int y=X[j]+d;if(y>=lo&&y<=hi&&map[y-lo]>=0)a[map[y-lo]]=0;}
    }
    bool rec(vector<unsigned char> a,vector<unsigned char> b){
        for(;;){
            int r=-1,q=numeric_limits<int>::max();
            for(int i=0;i<T;++i)if(!b[i]){
                int count=0;
                for(int j:R[i])count+=a[j];
                if(!count)return false;
                if(count<q)r=i,q=count;
            }
            if(r<0)return true;
            if(q==1){for(int j:R[r])if(a[j]){choose(j,a,b);break;}continue;}
            for(int j:R[r])if(a[j]){auto x=a,y=b;++branches;choose(j,x,y);if(rec(move(x),move(y)))return true;}
            return false;
        }
    }
};
int main() {
    Template t;
    int roots=0,repeated=0,decisions=0,unit=0,open=0,peak=0,total=0;
    auto verify=[&](const array<int,5>& deficits) {
        ++roots;
        Search s(t,deficits);
        if(s.duplicate){++repeated;return;}
        if(s.rec(vector<unsigned char>(s.X.size(),1),vector<unsigned char>(T,0)))
            ++open;
        else if(s.branches) ++decisions;
        else ++unit;
        peak=max(peak,s.branches);total+=s.branches;
    };
    for(int c=5;c>=2;--c) for(int b=20;b>c;--b) {
        if(n==9) verify({b,c,1,0,0});
        else for(int a=76;a>b;--a) verify({a,b,c,1,0});
    }
    need(roots==(n==9?66:4210),"unexpected root count");
    cout<<"Dimension "<<n<<", H="<<H<<", low cuts "<<t.cuts.front()<<" to "<<t.cuts.back()
        <<", "<<t.cuts.size()<<" saturation cuts; "<<t.U.size()<<" tested ranks.\n";
    cout<<roots<<" roots: "<<repeated<<" repeated, "<<unit<<" forced, "
        <<decisions<<" branched, "<<open<<" open.\n";
    cout<<total<<" decision children; peak "<<peak<<".\n";
    need(!open,"some tails survived the local relaxation");
    cout<<"PASS\n";
}
