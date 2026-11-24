// Independent verifier of the 526 equality obstruction. C++17, no input files.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

using namespace std;
constexpr int n=11, Q=1<<n, cap=526;
void need(bool x, const char* msg) {
    if (!x) { cerr<<"FAIL: "<<msg<<'\n'; exit(1); }
}

struct Cube {
    vector<int> order, rank, cuts, profile_cuts;
    int left=0,right=0;
    array<array<int,n+1>,Q+1> counts{};
    Cube() : order(Q), rank(Q) {
        iota(order.begin(),order.end(),0);
        vector<int> w(Q);
        for (int s=0;s<Q;++s) for (int j=0;j<n;++j)
            if (s&(1<<j)) w[s] += 4096-(1<<(n-1-j));
        sort(order.begin(),order.end(),[&](int s,int t){return w[s]<w[t];});
        for (int i=0;i<Q;++i) {
            need(i==0||w[order[i-1]]<w[order[i]],"strict canonical order");
            rank[order[i]]=i+1;
        }
        vector<vector<int>> enters(Q+1);
        for (int s=0;s<Q;++s) {
            int lo=rank[s];
            for (int j=0;j<n;++j) lo=min(lo,rank[s^(1<<j)]);
            enters[lo].push_back(rank[s]);
        }
        int have=0,highest=0,max_boundary=0;
        for (int p=1;p<=Q;++p) {
            for (int v: enters[p]) { ++have; highest=max(highest,v); }
            need(have==highest,"closed canonical 1-neighborhood is a prefix");
            max_boundary=max(max_boundary,have-p);
            if (have-p==cap) cuts.push_back(p);
            for (int j=0;j<=n;++j)
                counts[p][j]=counts[p-1][j]+(order[p-1]<(1<<j));
        }
        need(max_boundary==cap && cuts.size()==32,"sharp canonical boundary");
        left=cuts.front(); right=cuts.back()+cap;
        need(left==674 && right==1375 && left+right==Q+1,"symmetric band");
        profile_cuts=cuts;
        for (int p: cuts) profile_cuts.push_back(p+cap);
        sort(profile_cuts.begin(),profile_cuts.end());
        need(profile_cuts.size()==64 && adjacent_find(profile_cuts.begin(),profile_cuts.end())==profile_cuts.end(),
             "cut count");
        for (size_t i=1;i<cuts.size();++i) need(cuts[i]-cuts[i-1]<=cap,"no band gaps");
    }
    int lower(int j) const {
        int bound=1;
        for (int p: profile_cuts) for (int q: profile_cuts)
            if (counts[p][j]-counts[p][j-1]<counts[q][j-1])
                bound=max(bound,p-q+1);
        return bound;
    }
};

struct Instance {
    const Cube& cube;
    vector<int> at;
    vector<vector<int>> words;
    vector<int> hist;
    vector<int> diffs;
    int low,high,tau;
    vector<array<int,2>> orbit;
    vector<int> orbit_at;
    vector<vector<int>> row_opts,orbit_hits;
    bool high_duplicate=false;
    int extra_nodes=0;
    Instance(const Cube& cc, vector<int> tt):cube(cc),at(move(tt)),words(6),row_opts(cc.right-cc.left+1) {
        words[0]={0};
        for (int i=0;i<5;++i) {
            words[i+1]=words[i];
            for (int u:words[i]) words[i+1].push_back(u+at[i]);
        }
        for (auto& vec:words) sort(vec.begin(),vec.end());
        const auto& full=words.back();
        if (adjacent_find(full.begin(),full.end())!=full.end()) {
            high_duplicate=true;
            return;
        }
        int top_sum=full.back();
        tau=Q+1-top_sum;
        int max_low_sum=6*at.front()-21;
        auto floor_half=[](int x){return x>=0?x/2:(x-1)/2;};
        low=-floor_half(max_low_sum-tau);
        high=floor_half(tau+max_low_sum);
        need(low+high==tau,"reflection domain");
        hist.assign(high-low+2,0);
        const int W=high-low+1;
        orbit_at.assign(W,-1);
        auto mark=[&](int a,int b) {
            a=max(a,low); b=min(b,high);
            if (a<=b) { ++hist[a-low]; --hist[b-low+1]; }
        };
        // Interpret each interval by occupied high translates, independently
        // of the author's bitset implementation.
        for (int j=6;j<=n;++j) {
            vector<pair<int,int>> pts{{low-1,0}};
            for (int p:cube.profile_cuts) pts.emplace_back(p,cube.counts[p][j]);
            pts.emplace_back(high+top_sum,1<<j);
            auto& prefix_high=words[j-6];
            for (size_t k=1;k<pts.size();++k) {
                auto [a,c]=pts[k-1]; auto [b,d]=pts[k];
                need(b>=a && d>=c && d-c<=b-a,"cut count interval");
                for(int e:full) {
                    bool in_prefix=binary_search(prefix_high.begin(),prefix_high.end(),e);
                    if ((in_prefix && d==c) || (!in_prefix && d-c==b-a))
                        mark(a+1-e,b-e);
                }
            }
        }
        for (int u:full) for (int v:full) if(u!=v) diffs.push_back(u-v);
        sort(diffs.begin(),diffs.end());
        diffs.erase(unique(diffs.begin(),diffs.end()),diffs.end());
        partial_sum(hist.begin(),hist.end(),hist.begin());
        for (int a=low;a<tau-a;++a) {
            int b=tau-a;
            if (hist[a-low] || hist[b-low] || binary_search(diffs.begin(),diffs.end(),b-a)) continue;
            int idx=int(orbit.size());
            orbit.push_back({a,b});
            orbit_at[a-low]=orbit_at[b-low]=idx;
            orbit_hits.emplace_back();
            for (int x:{a,b}) for (int e:full) {
                int z=x+e;
                if (z<cube.left || z>cube.right) continue;
                int row=z-cube.left;
                row_opts[row].push_back(idx);
                orbit_hits[idx].push_back(row);
            }
        }
        for (auto& opts:row_opts)
            need(adjacent_find(opts.begin(),opts.end())==opts.end(),"within-orbit coverage is disjoint");
    }
    void select(int k, vector<unsigned char>& viable, vector<unsigned char>& covered) const {
        need(viable[k],"select only a live orbit");
        viable[k]=0;
        for(int row:orbit_hits[k]) {
            need(!covered[row],"selected orbits cannot overlap");
            covered[row]=1;
        }
        for (int x: orbit[k]) for (int d:diffs) {
            int y=x+d;
            if(y<low||y>high) continue;
            int v=orbit_at[y-low];
            if(v>=0) viable[v]=0;
        }
    }
    // A deliberately weaker relaxation: does not impose exactly 32 orbits,
    // or factorability as a six-weight subset-sum set.
    bool can_cover(vector<unsigned char> viable, vector<unsigned char> covered) {
        while (true) {
            int best=-1, qty=numeric_limits<int>::max();
            for (int i=0;i<int(row_opts.size());++i) {
                if (covered[i]) continue;
                int cnt=0;
                for (int k: row_opts[i]) cnt += viable[k];
                if (cnt==0) return false;
                if (cnt<qty) { best=i; qty=cnt; }
            }
            if (best<0) return true;
            if (qty==1) {
                auto it=find_if(row_opts[best].begin(),row_opts[best].end(),[&](int k){return viable[k];});
                need(it!=row_opts[best].end(),"forced coverage orbit");
                select(*it,viable,covered);
                continue;
            }
            for (int k: row_opts[best]) if (viable[k]) {
                auto live=viable, cov=covered;
                select(k,live,cov);
                ++extra_nodes;
                if(can_cover(move(live),move(cov))) return true;
            }
            return false;
        }
    }
};

int main() {
    Cube cube;
    array<int,4> lo={cube.lower(7),cube.lower(8),cube.lower(9),cube.lower(10)};
    need(lo==array<int,4>{450,506,521,525},"verified top-coordinate constraints");
    int roots=0, duplicates=0, forced=0, branched=0, max_nodes=0, total_nodes=0;
    for(int a10=lo[3];a10<cap;++a10)
        for(int a9=lo[2];a9<a10;++a9)
            for(int a8=lo[1];a8<a9;++a8)
                for(int a7=lo[0];a7<a8;++a7) {
                    ++roots;
                    Instance f(cube,{a7,a8,a9,a10,cap});
                    if (f.high_duplicate) { ++duplicates; continue; }
                    need(!f.can_cover(vector<unsigned char>(f.orbit.size(),1),
                                      vector<unsigned char>(f.row_opts.size(),0)),
                         "surviving tail");
                    if (f.extra_nodes) ++branched; else ++forced;
                    max_nodes=max(max_nodes,f.extra_nodes);
                    total_nodes+=f.extra_nodes;
                }
    need(roots==4210,"complete tail count");
    cout<<"Band ["<<cube.left<<','<<cube.right<<"], "<<cube.cuts.size()<<" boundary cuts, "
        <<cube.profile_cuts.size()<<" profile cuts\n";
    cout<<"Lower bounds a7 through a10:";
    for(int v:lo) cout<<' '<<v;
    cout<<"\nRoots "<<roots<<": "<<duplicates<<" duplicate, "<<forced
        <<" forced, "<<branched<<" branched\n"<<"Branch child nodes "
        <<total_nodes<<"; peak per root "<<max_nodes<<"\nPASS\n";
}
