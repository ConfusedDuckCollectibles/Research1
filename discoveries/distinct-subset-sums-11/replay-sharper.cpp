// C++17; integer arithmetic and standard library only; no input files.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

void require(bool test, const char* what) {
    if (!test) { cerr << "FAIL: " << what << '\n'; exit(1); }
}

int binom(int n, int k) {
    int b = 1;
    for (int i = 1; i <= k; ++i) b = b * (n+1-i)/i;
    return b;
}

struct Window {
    int n, k, H=0, L, W, cap, trim, masks;
    vector<int> order, rank, cuts, sites, site_index, beta;
    vector<uint32_t> banned;
    array<int,5> deficit={};

    void lex(int next, int left, int set) {
        if (!left) { order.push_back(set); return; }
        for (int i=next; i<=n-left; ++i) lex(i+1,left-1,set|(1<<i));
    }

    Window(int dimension, int dropped)
      : n(dimension), k(n==9?4:5), cap(n==9?4:5), trim(dropped),
        masks(1<<k), rank(1<<n) {
        for (int i=0; i<n; ++i) H += binom(i,i/2);
        for (int size=0; size<=n; ++size) lex(0,size,0);
        for (unsigned i=0; i<order.size(); ++i) rank[order[i]]=i+1;
        int last=1<<(n-1);
        for (int s : order) {
            if (s&last) continue;
            bool ok=true;
            for (int i=n-3; i>=0; i-=2) {
                int pair=(s>>i)&3;
                if (pair!=1 && pair!=2) ok=false;
            }
            if (ok) {
                cuts.push_back(rank[s]);
                require(rank[s|last]-rank[s]==H,"saturated cuts have the bandwidth shift");
            }
        }
        require(cuts.size()==(1u<<(n/2)), "all saturated cuts");
        L=cuts.front(); W=cuts.back()-L;
        require(W== (n==9?49:175),"expected cut width");
        for (int i=0;i<k-1;++i) {
            int best=1<<20;
            int not_allowed=((1<<n)-1)^((1<<(n-k+i))-1);
            for (int p:cuts) for (int q:cuts) if (q>p && q-p-1<best)
                for (int r=p+1;r<=q;++r)
                    if (!(order[r-1]&not_allowed)) { best=q-p-1; break; }
            deficit[i]=best;
        }
        require(deficit==(n==9?array<int,5>{20,5,1,0,0}:array<int,5>{76,20,5,1,0}),
                "tail deficit witnesses");

        sites.clear();
        site_index.assign(H+W+1,-1);
        for (int r=0;r<=H+W;++r) if ((r>=trim && r<=W+63) || r>=H) {
            site_index[r]=sites.size();
            sites.push_back(r);
        }
        beta.assign(H+W+1,-1);
        beta[0]=0;
        for (int r:sites) for (int s:sites) if (r<s) {
            // A cut p spans two ranks iff s+L-H <= p <= r+L.
            auto it=lower_bound(cuts.begin(),cuts.end(),s+L-H);
            int stretch=it!=cuts.end() && *it<=r+L ? 1 : 2;
            beta[s-r]=max(beta[s-r],stretch);
        }

        // Reconstruct the direct and complement rank restrictions separately.
        banned.assign(H+W+1,0);
        for (int shift : {0,H}) for (unsigned i=1;i<cuts.size();++i) {
            int lo=cuts[i-1]+shift, hi=cuts[i]+shift;
            for (int m=0;m<k;++m) {
                int suffix=(masks-1)^((1<<m)-1);
                for (bool reflected : {false,true}) {
                    int a=lo+!reflected,b=hi+!reflected;
                    int count=0;
                    for (int r=a;r<b;++r) {
                        int t=(order[r-1]>>(n-k))&suffix;
                        count += t==(reflected?suffix:0);
                    }
                    if (count!=0&&count!=b-a) continue;
                    for (int t=0;t<masks;++t) {
                        bool fits=(t&suffix)==(reflected?suffix:0);
                        if (fits!=(count!=0)) for (int r=a;r<b;++r)
                            banned[r-L]|=uint32_t(1)<<t;
                    }
                }
            }
        }
    }

    int Wmask(int a) const {
        int value=0;
        for (int i=0;i<k-1;++i) if ((a>>i)&1) value+= i==k-2?1:2;
        return value;
    }

    void shifts(const vector<int>& e, long long& tested, int& max_coefficient,
                set<int>& bad_distances, int& bad_tails) const {
        bool bad=false;
        for (int p=0;p<masks;++p) for (int q=0;q<masks;++q) {
            int f=(q>>(k-1))-(p>>(k-1));
            int A=Wmask(q&~p)+f, B=Wmask(p&~q)-f;
            for (int d=e[q]-e[p]-1;d<=e[q]-e[p]+1;++d) {
                int u=abs(d);
                if (u>H+W || beta[u]<0) continue;
                ++tested;
                int c=max(A+(d<0?beta[u]:0),B+(d>0?beta[u]:0));
                max_coefficient=max(max_coefficient,c);
                if (c>cap) { bad=true; bad_distances.insert(u); }
            }
        }
        bad_tails+=bad;
    }
};

struct ExactCover {
    vector<vector<int>> hits, rows;
    int branches=0;
    ExactCover(int sites, vector<vector<int>> incidences)
        : hits(move(incidences)), rows(sites) {
        for(unsigned j=0;j<hits.size();++j) for(int i:hits[j]) rows[i].push_back(j);
    }
    void select(int row, vector<uint8_t>& live, vector<uint8_t>& covered) const {
        require(live[row],"select live row");
        for (int col : hits[row]) {
            require(!covered[col],"no overlap with earlier selected rows");
            covered[col]=1;
            for (int j : rows[col]) live[j]=0;  // No off-site collision rules.
        }
    }
    bool search(vector<uint8_t> live, vector<uint8_t> covered) {
        for (;;) {
            int best=-1, score=1<<30, last=-1;
            for (unsigned i=0;i<rows.size();++i) if (!covered[i]) {
                int count=0, last_here=-1;
                for (int r:rows[i]) if (live[r]) { ++count;last_here=r; }
                if (!count) return false;
                if (count<score) { score=count;best=i;last=last_here; }
            }
            if (best<0) return true;
            if (score==1) { select(last,live,covered);continue; }
            for (int r:rows[best]) if (live[r]) {
                auto l=live,c=covered;
                select(r,l,c); ++branches;
                if (search(move(l),move(c))) return true;
            }
            return false;
        }
    }
    bool run() { return search(vector<uint8_t>(hits.size(),1),vector<uint8_t>(rows.size(),0)); }
};

vector<vector<int>> make_rows(const Window& w, const vector<int>& e) {
    int lo=-w.k*w.H, hi=w.H+w.W+102;
    vector<uint8_t> candidates(hi-lo+1);
    for (int r:w.sites) for (int b:e) {
        int idx=r-b-lo;
        require(idx>=0 && idx<=hi-lo,"tag in the enumerated domain");
        candidates[idx]=1;
    }
    vector<vector<int>> rows;
    for (int z=lo;z<=hi;++z) if (candidates[z-lo]) {
        vector<int> incidence;
        bool ok=true;
        for (unsigned m=0;m<e.size();++m) {
            int r=z+e[m];
            if (r<0 || r>w.H+w.W) continue;
            int site=w.site_index[r];
            if (site<0) continue;
            if ((w.banned[r]>>m)&1) {ok=false;break;}
            incidence.push_back(site);
        }
        if (ok) {
            require(!incidence.empty(),"candidate generated from a retained site");
            rows.push_back(move(incidence));
        }
    }
    return rows;
}

void run(int dimension, int trim, bool control=false) {
    Window w(dimension,trim);
    int roots=0, repeated=0, forced=0, branched=0, children=0, peak=0,
        max_coeff=0, bad_tails=0;
    long long shifts=0;
    set<int> bad_distances;
    array<int,5> d={};
    auto process=[&]() {
        ++roots;
        vector<int> e(w.masks);
        for (int m=0; m<w.masks; ++m)
            for (int i=0;i<w.k;++i) if ((m>>i)&1) e[m]+=w.H-d[i];
        w.shifts(e,shifts,max_coeff,bad_distances,bad_tails);
        if (control) return;
        auto sorted=e;
        sort(sorted.begin(),sorted.end());
        if (adjacent_find(sorted.begin(),sorted.end())!=sorted.end()) {++repeated;return;}
        ExactCover cover(w.sites.size(),make_rows(w,e));
        require(!cover.run(),"a distinct-floor tail admits a retained-rank exact cover");
        if (cover.branches) ++branched; else ++forced;
        children+=cover.branches; peak=max(peak,cover.branches);
    };
    d[w.k-2]=1;
    for (int c=2;c<=w.deficit[w.k-3];++c) {
        d[w.k-3]=c;
        for (int b=c+1;b<=w.deficit[w.k-4];++b) {
            d[w.k-4]=b;
            if (w.k==4) process();
            else for (int a=b+1;a<=w.deficit[0];++a) {d[0]=a;process();}
        }
    }
    require(roots==(w.n==9?66:4210),"complete tail enumeration");
    cout << "Q" << w.n << " H=" << w.H << " cuts " << w.L << ".." << w.L+w.W
         << "; " << w.sites.size() << " retained sites (trim " << trim << "). "
         << shifts << " ordered mask/near-shift cases; max coefficient " << max_coeff << ".\n";
    if (control) {
        require(bad_tails==10 && bad_distances==set<int>{447,448,449},
                "untrimmed Q10 reaches 6, as expected");
        cout << "Calibration: " << bad_tails << " tails fail the coefficient-5 bound at absolute shifts";
        for (int i:bad_distances) cout << ' ' << i;
        cout << ".\n";
        return;
    }
    require(!bad_tails,"required near-shift coefficient");
    require(roots==repeated+forced+branched,"all tails excluded");
    cout << roots << " tails: " << repeated << " duplicates, " << forced << " forced, "
         << branched << " branched, 0 open. " << children << " decision children, max " << peak << ".\n";
}

int main() {
    // In this control E={0,2}, two translates {0,2} and {2,4} tile only the
    // retained sites {0,4}. Their collision at 2 must be ignored.
    Window test(9,0);
    test.sites={0,4};
    fill(test.site_index.begin(),test.site_index.end(),-1);
    test.site_index[0]=0; test.site_index[4]=1;
    fill(test.banned.begin(),test.banned.end(),0);
    auto incidence=make_rows(test,{0,2});
    require(incidence==vector<vector<int>>{{0},{0},{1},{1}},"incidence ignores unqueried rank 2");
    ExactCover control(2,move(incidence));
    require(control.run(),"disjointness is required at retained sites only");
    run(9,0);
    run(10,0,true);
    run(10,3);
    run(11,0);
    cout << "PASS\n";
}
