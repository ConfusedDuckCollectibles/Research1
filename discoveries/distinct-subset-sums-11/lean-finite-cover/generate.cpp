// Original C++17 generator adapted from the independent retained-site replay.
// No input files; writes a compact canonical JSON corpus to stdout.
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

struct Certificate {
    const Window& w;
    vector<vector<int>> rows, through;
    vector<int> proof;
    Certificate(const Window& window, const vector<int>& e) : w(window), through(w.sites.size()) {
        int lo=-w.k*w.H, hi=w.H+w.W;
        vector<int> tagToRow(hi-lo+1,-1);
        vector<bool> seen(hi-lo+1);
        for(int u:w.sites)for(int b:e) {
            int tag=u-b,index=tag-lo; require(index>=0&&index<=hi-lo,"tag array index");
            if(seen[index])continue;
            seen[index]=true;
            bool ok=true;
            vector<int> hits;
            for(unsigned m=0;m<e.size();++m){
                int r=tag+e[m];
                if(r<0||r>w.H+w.W||w.site_index[r]<0)continue;
                if(w.banned[r]>>m&1){ok=false;break;}
                hits.push_back(w.site_index[r]);
            }
            if(ok){
                require(!hits.empty(),"retained-target tag");
                tagToRow[index]=rows.size();
                rows.push_back(move(hits));
            }
        }
        // Follow natural mask order, as in the Lean pure row generator.
        for(unsigned s=0;s<w.sites.size();++s) for(int b:e){
            int i=tagToRow[w.sites[s]-b-lo];
            if(i>=0)through[s].push_back(i);
        }
    }
    static bool fits(const vector<int>& row,const vector<bool>& used){
        for(int s:row)if(used[s])return false;
        return true;
    }
    bool emit(vector<bool>& used){
        int site=-1,mincount=int(rows.size())+1;
        for(unsigned s=0;s<used.size();++s)if(!used[s]){
            int count=0;
            for(int r:through[s])count+=fits(rows[r],used);
            if(count<mincount){mincount=count;site=s;}
            if(!count)break;
        }
        if(site<0)return false;
        proof.push_back(site);
        for(int i:through[site]) if(fits(rows[i],used)){
            for(int s:rows[i])used[s]=true;
            if(!emit(used))return false;
            for(int s:rows[i])used[s]=false;
        }
        return true;
    }
    bool read(const vector<int>& cert, unsigned& pos, vector<bool>& used)const{
        if(pos>=cert.size())return false;
        int s=cert[pos++];
        if(s<0||s>=int(used.size())||used[s])return false;
        for(int i:through[s])if(fits(rows[i],used)){
            auto next=used;
            for(int r:rows[i])next[r]=true;
            if(!read(cert,pos,next))return false;
        }
        return true;
    }
    bool passes(const vector<int>& cert)const{
        unsigned pos=0;
        vector<bool> used(w.sites.size());
        return read(cert,pos,used)&&pos==cert.size();
    }
    void run(){
        vector<bool> used(w.sites.size());
        require(emit(used)&&passes(proof),"complete refutation produced and replayed");
        auto bad=proof; bad.pop_back();require(!passes(bad),"truncate");
        bad=proof;bad.push_back(0);require(!passes(bad),"extend");
        bad=proof;bad[0]=w.sites.size();require(!passes(bad),"invalid site");
    }
};
template<class T> void arrayJson(const vector<T>& xs){
    cout<<'[';
    for(unsigned i=0;i<xs.size();++i)cout<<(i?",":"")<<xs[i];
    cout<<']';
}
void generate(int n,int trim){
    Window w(n,trim);
    cout<<"{\"dimension\":"<<n<<",\"h\":"<<w.H<<",\"width\":"<<w.W
        <<",\"trim\":"<<w.trim<<",\"tail\":"<<w.k<<",\"bans\":";
    arrayJson(w.banned);
    cout<<",\"entries\":[";
    int count=0,duplicated=0,indices=0,maxLen=0;
    array<int,5> d={}; d[w.k-2]=1;
    auto entry=[&](){
        vector<int> e(w.masks);
        for(int m=0;m<w.masks;++m)for(int i=0;i<w.k;++i)if(m>>i&1)e[m]+=w.H-d[i];
        auto sorted=e;sort(sorted.begin(),sorted.end());
        if(adjacent_find(sorted.begin(),sorted.end())!=sorted.end()){++duplicated;return;}
        Certificate c(w,e);c.run();
        if(count++)cout<<",\n";
        cout<<'[';arrayJson(vector<int>(d.begin(),d.begin()+w.k));cout<<',';
        arrayJson(c.proof);cout<<",[";
        set<int> cited(c.proof.begin(), c.proof.end());
        bool firstSite=true;
        for(int s:cited) {
            if(!firstSite)cout<<',';
            firstSite=false;cout<<'['<<s<<",[";
            bool firstRow=true;
            for(int r:c.through[s]) {
                if(!firstRow)cout<<',';
                firstRow=false;arrayJson(c.rows[r]);
            }
            cout<<"]]";
        }
        cout<<"]]";
        indices+=c.proof.size();maxLen=max(maxLen,int(c.proof.size()));
    };
    for(int c=2;c<=w.deficit[w.k-3];++c){
        d[w.k-3]=c;
        for(int b=c+1;b<=w.deficit[w.k-4];++b){
            d[w.k-4]=b;
            if(w.k==4)entry();
            else for(int a=b+1;a<=w.deficit[0];++a){d[0]=a;entry();}
        }
    }
    cout<<"]}";
    cerr<<"Q"<<n<<" "<<count<<" certs, "<<duplicated<<" duplicate-floor tails, "
        <<indices<<" indices, max "<<maxLen<<", "<<3*count<<" corruptions rejected.\n";
}
int main(){
    cout<<"[\n";generate(9,0);cout<<",\n";generate(10,3);cout<<",\n";
    generate(11,0);cout<<"\n]\n";
}
