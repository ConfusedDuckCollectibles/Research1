// Standalone symbolic-coordinate replay of the dimension-uniform obstruction.
// C++17, no files or external libraries. Positions are (H coefficient, residual).
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>
using namespace std;

void need(bool t,const char* msg) {
    if(!t) {cerr<<"FAIL: "<<msg<<'\n'; exit(1);}
}
constexpr int n=11, Q=1<<n, H=526, tail_mask=31, W=175, extend=238, tail_d_max=102;
constexpr int width=extend+tail_d_max+1;
constexpr int lowest=-5, highest=1, domains=(highest-lowest+1)*width;

struct Pos {
    int slab, delta;
    bool operator==(const Pos& o) const {return slab==o.slab&&delta==o.delta;}
    bool operator<(const Pos& o) const {return tie(slab,delta)<tie(o.slab,o.delta);}
    Pos operator+(const Pos& o) const {return {slab+o.slab,delta+o.delta};}
    Pos operator-(const Pos& o) const {return {slab-o.slab,delta-o.delta};}
};
int domain_index(Pos x) {
    if (x.slab<lowest||x.slab>highest||x.delta<0||x.delta>=width) return -1;
    return (x.slab-lowest)*width+x.delta;
}
int target_index(Pos x) {
    if (x.delta<0) return -1;
    if (x.slab==0 && x.delta<=extend) return x.delta;
    if (x.slab==1 && x.delta<=W) return extend+1+x.delta;
    return -1;
}
constexpr int target_count=extend+1+W+1;

struct Template {
    array<array<uint32_t,W+1>,2> bans{};
    vector<int> canon, ranks, cuts;
    vector<Pos> targets;
    void generate(int k,int min_index,int chosen) {
        if (k==0) {canon.push_back(chosen);return;}
        for(int i=min_index;i<=n-k;++i) generate(k-1,i+1,chosen|(1<<i));
    }
    Template():ranks(Q) {
        static_assert(H>extend+2*tail_d_max,"disjoint residue slabs");
        for(int k=0;k<=n;++k) generate(k,0,0);
        need(canon.size()==Q,"complete Q13");
        for(int m=1;m<=Q;++m) ranks[canon[m-1]]=m;
        const int max_bit=1<<(n-1);
        for (int m=1;m<=Q;++m) {
            const int S=canon[m-1];
            if (S&max_bit) continue;
            bool saturated=true;
            for(int b=(n+1)%2;b<n-1;b+=2) {
                const int c=(S>>b)&3;
                if (c!=1&&c!=2) saturated=false;
            }
            if(saturated) {
                need(ranks[S|max_bit]-m==H,"exact high-cut shift");
                cuts.push_back(m);
            }
        }
        need(cuts.size()==(1u<<(n/2)) && cuts.back()-cuts.front()==W,"saturation template");
        int anchor=cuts.front();
        // Verify tail-gap witnesses using direct comparisons among the low cuts.
        array<int,4> want={76,20,5,1};
        for(int k=0;k<4;++k) {
            int best=10000;
            for(size_t u=0;u<cuts.size();++u)
                for(size_t v=u+1;v<cuts.size();++v) {
                    int a=cuts[u],b=cuts[v];
                    bool rank_gain=false;
                    // j=n-4+k. Allow through j-1, i.e. lower k tail bits.
                    for(int rank=a+1;rank<=b;++rank)
                        if ((canon[rank-1]>>(n-5+k))==0) rank_gain=true;
                    if(rank_gain) best=min(best,b-a-1);
                }
            need(best==want[k],"exact five-tail deficit bounds");
        }
        for(int stage=0;stage<=1;++stage)
            for(size_t block=1;block<cuts.size();++block) {
                const int a=cuts[block-1]+stage*H, b=cuts[block]+stage*H;
                const int origin=anchor+stage*H;
                for(int k=0;k<=5;++k) {
                    const int beyond=tail_mask^((1<<k)-1);
                    int direct=0, dual=0;
                    for(int r=a+1;r<=b;++r)
                        if (((canon[r-1]>>(n-5))&beyond)==0) ++direct;
                    for(int r=a;r<b;++r)
                        if (((canon[r-1]>>(n-5))&beyond)==beyond) ++dual;
                    for(int t=0;t<=tail_mask;++t) {
                        bool subset=(t&beyond)==0, superset=(t&beyond)==beyond;
                        if ((direct==0&&subset)||(direct==b-a&&!subset))
                            for(int r=a+1;r<=b;++r) bans[stage][r-origin]|=uint32_t(1)<<t;
                        if ((dual==0&&superset)||(dual==b-a&&!superset))
                            for(int r=a;r<b;++r) bans[stage][r-origin]|=uint32_t(1)<<t;
                    }
                }
            }
        for(int r=0;r<=extend;++r) targets.push_back({0,r});
        for(int r=0;r<=W;++r) targets.push_back({1,r});
        need(targets.size()==target_count,"target count");
    }
};

struct Search {
    const Template& cut;
    vector<Pos> E, x, differences;
    array<int,domains> x_at;
    vector<vector<int>> by_target, by_candidate;
    bool duplicate=false;
    long long branches=0;
    Search(const Template& c,const array<int,5>& deficit):cut(c),by_target(target_count) {
        int sum=0; for(int d:deficit)sum+=d;
        need(sum<=tail_d_max,"residue bound");
        for(int m=0;m<=tail_mask;++m) {
            Pos v{0,0};
            for(int i=0;i<5;++i) if (m&(1<<i)) {++v.slab;v.delta-=deficit[i];}
            E.push_back(v);
        }
        auto Q0=E;
        sort(Q0.begin(),Q0.end());
        if(adjacent_find(Q0.begin(),Q0.end())!=Q0.end()) {duplicate=true;return;}
        for(const auto& p:E)for(const auto& q:E)if(!(p==q))differences.push_back(p-q);
        sort(differences.begin(),differences.end());
        differences.erase(unique(differences.begin(),differences.end()),differences.end());
        x_at.fill(-1);
        for(Pos u:cut.targets) for(Pos e:E) {
            Pos v=u-e;
            int i=domain_index(v);
            need(i>=0,"finite coordinate domain");
            x_at[i]=-2; // U - E is the complete search domain before exclusions.
        }
        for(int slab=lowest;slab<=highest;++slab)for(int a=0;a<width;++a) {
            Pos v{slab,a};
            int i=domain_index(v);
            if(x_at[i]!=-2)continue;
            vector<int> hit;
            bool bad=false;
            for(int m=0;m<=tail_mask;++m) {
                Pos z=v+E[m];
                if(z.slab>=0&&z.slab<=1&&z.delta>=0&&z.delta<=W
                    && (cut.bans[z.slab][z.delta]&(uint32_t(1)<<m))) {
                    bad=true; break;
                }
                int k=target_index(z);
                if(k>=0)hit.push_back(k);
            }
            if(bad){x_at[i]=-1;continue;}
            need(!hit.empty(),"eligible candidate covers one of the target ranks");
            int k=int(x.size());
            x_at[i]=k; x.push_back(v); by_candidate.push_back(hit);
            for(int r:hit)by_target[r].push_back(k);
        }
    }
    void choose(int j,vector<unsigned char>& live,vector<unsigned char>& covered) {
        need(live[j],"chosen representative is live");
        live[j]=0;
        for (int r:by_candidate[j]) {
            need(!covered[r],"covered rank collision");
            covered[r]=1;
        }
        for(Pos d:differences) {
            int i=domain_index(x[j]+d);
            if(i>=0 && x_at[i]>=0)live[x_at[i]]=0;
        }
    }
    bool recurse(vector<unsigned char> live,vector<unsigned char> covered) {
        for(;;) {
            int best=-1, num=numeric_limits<int>::max();
            for(int i=0;i<target_count;++i) if(!covered[i]) {
                int cnt=0;
                for(int v:by_target[i]) cnt+=live[v];
                if(cnt==0)return false;
                if(cnt<num){best=i;num=cnt;}
            }
            if(best<0)return true;
            if(num==1) {
                auto it=find_if(by_target[best].begin(),by_target[best].end(),[&](int v){return live[v];});
                need(it!=by_target[best].end(),"unit representative");
                choose(*it,live,covered);
                continue;
            }
            for(int v:by_target[best])if(live[v]) {
                auto ll=live, cc=covered;
                ++branches;
                choose(v,ll,cc);
                if(recurse(move(ll),move(cc)))return true;
            }
            return false;
        }
    }
};
int main() {
    Template table;
    cout<<"Symbolic blocks: "<<table.cuts.front()<<" to "<<table.cuts.back()
        <<"; each shifted by "<<H<<"; "<<table.cuts.size()<<" cuts per block, "<<table.targets.size()<<" local targets."<<endl;
    int duplicate=0, forced=0, branched=0, total=0, peak=0, roots=0;
    for(int c=5;c>=2;--c)for(int b=20;b>c;--b)for(int a=76;a>b;--a) {
        ++roots;
        Search s(table,{a,b,c,1,0});
        if(s.duplicate){++duplicate;continue;}
        need(!s.recurse(vector<unsigned char>(s.x.size(),1),vector<unsigned char>(target_count,0)),
             "symbolic tail survives");
        if(s.branches) ++branched;else ++forced;
        total+=s.branches; peak=max(peak,int(s.branches));
    }
    need(roots==4210,"complete normalized tail count");
    cout<<roots<<" normalized tails: "<<duplicate<<" duplicate, "<<forced<<" unit-only, "
        <<branched<<" branched, no survivors.\n"<<total<<" decision children; maximum "<<peak<<" per root.\nPASS\n";
}
