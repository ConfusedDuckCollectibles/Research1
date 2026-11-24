// Independent C++17 falsification search. No code or data imports.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <vector>
using namespace std;
void req(bool p,const char* m){if(!p){cerr<<"FAIL "<<m<<'\n';exit(1);}}
const array<array<int,4>,3> special={{{2,0,0,2},{2,2,0,0},{2,0,0,2}}};
const array<array<int,2>,3> boundary={{{1,2},{1,2},{0,1}}};

struct Machine {
    int n,k,N;
    vector<array<int,2>> edges;
    Machine(int size,int type,const vector<int>& tail)
      :n(size),k(type),N(1<<size),edges(n) {
        for(int i=0;i<4;++i) edges[i/2][i%2]=special[k][i];
        for(int i=0;i<2*(n-2);++i) edges[2+i/2][i%2]=tail[i];
    }
    bool strong() const {
        int a=1,b=1;
        for(int d=0;d<n-1;++d) for(int q=0;q<n;++q)for(int x:edges[q]) {
            if(a>>q&1) a|=1<<x;
            if(b>>x&1) b|=1<<q;
        }
        return a==N-1&&b==N-1;
    }
    int collision() const {
        // Layer propagation uses two bit masks, unlike a product-graph queue.
        array<int,2> r={6,1};
        for(int t=0;t<=2*n;++t) {
            if(r[0]&r[1])return t;
            auto s=r;
            for(int c=0;c<2;++c)for(int v=0;v<n;++v)if(r[c]>>v&1) {
                s[c^1]|=1<<edges[v][0];
                s[c^(k==1)]|=1<<edges[v][1];
            }
            if(s==r)return -1;
            r=s;
        }
        req(false,"finite product graph stabilizes"); return -1;
    }
    int reset() const {
        vector<int> d(N,-1), queue={N-1};d.back()=0;
        for(unsigned i=0;i<queue.size();++i) {
            int x=queue[i];
            if(!(x&(x-1)))return d[x];
            for(int label=0;label<2;++label){
                int y=0;
                for(int z=x;z;z&=z-1)y|=1<<edges[__builtin_ctz(unsigned(z))][label];
                if(d[y]<0){d[y]=d[x]+1;queue.push_back(y);}
            }
        }
        return -1;
    }
    int one_letter_cycles(int l) const {
        vector<int> pass(n,-1);int c=0;
        for(int q=0;q<n;++q) {
            int p=q;
            while(pass[p]<0){pass[p]=q;p=edges[p][l];}
            if(pass[p]==q)++c;
        }
        return c;
    }
    vector<int> indeg() const {
        vector<int> d(n);
        for(auto a:edges)for(int q:a)++d[q];
        return d;
    }
    // Evaluate both orders of x,y and both orders of letters for every
    // witness. Return the counts of excluded and other degree-two witnesses.
    array<int,2> degree_two() const {
        auto d=indeg(); array<int,2> ans={};
        for(int x=0;x<n;++x)for(int y=x+1;y<n;++y)
            for(int s=0;s<n;++s) if(s!=x&&s!=y) {
                int mask=(1<<x)|(1<<y)|(1<<s);
                bool valid=true;
                for(int q=0;q<n;++q)if(!(mask>>q&1)&&d[q]!=2)valid=false;
                for(int q:{x,y})for(int r:edges[q])if(!(mask>>r&1))valid=false;
                if(!valid)continue;
                bool excluded=false;
                for(int rev=0;rev<2;++rev) {
                    int X=rev?y:x,Y=rev?x:y; vector<int> name(n,-1);
                    name[X]=0;name[Y]=1;name[s]=2;
                    array<int,2> border; unsigned take=0;
                    for(int q=0;q<n;++q)if(q!=X&&q!=Y)
                        for(int r:edges[q])if(name[r]>=0) {
                            req(take<2,"two boundary arcs");
                            border[take++]=name[r];
                        }
                    req(take==2,"two boundary arcs");
                    sort(border.begin(),border.end());
                    for(int l=0;l<2;++l) {
                        array<int,4> row={name[edges[X][l]],name[edges[X][1-l]],
                                            name[edges[Y][l]],name[edges[Y][1-l]]};
                        for(int j=0;j<3;++j) if(row==special[j]&&border==boundary[j])excluded=true;
                    }
                }
                ++ans[!excluded];
            }
        return ans;
    }
    int degree_one() const {
        auto d=indeg();int ret=0;
        for(int x=0;x<n;++x)for(int s=0;s<n;++s)if(x!=s){
            bool ok=true;
            for(int rx:edges[x])if(rx!=x&&rx!=s)ok=false;
            for(int j=0;j<n;++j)if(j!=x&&j!=s&&d[j]!=2)ok=false;
            ret+=ok;
        }
        return ret;
    }
    bool uniform_letters_for_all_one_exception_vectors() const {
        // Three distinct double same-label images with disjoint source
        // pairs are sufficient. One target then lies off the exceptional
        // state and neither incoming arc uses it, forcing p=1/2.
        vector<vector<int>> src(n);
        for(int q=0;q<n;++q)for(int l=0;l<2;++l)src[edges[q][l]].push_back(2*q+l);
        vector<int> groups;
        for(int q=0;q<n;++q)if(src[q].size()==2&&(src[q][0]&1)==(src[q][1]&1)) {
            groups.push_back((1<<q)|(1<<(src[q][0]/2))|(1<<(src[q][1]/2)));
        }
        for(int s=0;s<n;++s) {
            bool unaffected=false;
            for(int m:groups)if(!(m>>s&1))unaffected=true;
            if(!unaffected)return false;
        }
        return true;
    }
    array<int,2> letter_ranks() const {
        array<int,2> image{};
        for(auto t:edges)for(int l=0;l<2;++l)image[l]|=1<<t[l];
        return {__builtin_popcount(unsigned(image[0])),__builtin_popcount(unsigned(image[1]))};
    }
};

int even_cost(int n,int k,int j) {
    return k==2?max(j-1,n-j-1):max(j,n-1-j);
}
int count_odd(int n,int k,int h) {
    int result=0;
    int upper=k==2?2*n:2*n-2, from=3;
    for(int m=from;m<upper;m+=2) {
        int lower=h<=2?1:k==2?2*(h/2)+1:2*((h+1)/2)-1;
        if(m>=lower && upper-m>=lower)++result;
    }
    // Also verify the clipped closed expression, including h=1 and h=2.
    int L=k==2?n-1:n-2, f=k==2?n-2*(h/2):n+1-2*((h+1)/2);
    req(result==(h<=2?L:max(0,min(L,f))),"odd mass formula");
    return result;
}
int bound(int n,int k,int tau) {
    int E=0,D=min(n-1,tau+1);
    for(int j=(k==2?2:1);j<=(k==2?n-1:n-2);++j)E+=even_cost(n,k,j);
    int O=0;
    for(int h=1;h<=D;++h)O+=count_odd(n,k,h);
    req(O==(k==2?n*D-D*D/2:n*D-(D*D+1)/2-2),"closed odd-layer sum");
    return E+O;
}

struct Profile {
    vector<vector<uint8_t>> allowed;
    int sum=0;
    Profile(const Machine& A):allowed(A.n) {
        int n=A.n,k=A.k,W=k==2?2*n:2*n-2,start=k==2?3:2;
        allowed[1].resize(W+1);
        for(int j=start;j<W;++j)allowed[1][j]=1;
        sum=W-start;
        int c0=6,c1=1;
        vector<int> weight(A.n,2);weight[0]=(k==2?3:1);weight[1]=1;
        auto mass=[&](int s){int z=0;for(int i=0;i<n;++i)if(s>>i&1)z+=weight[i];return z;};
        for(int h=2;h<n;++h) {
            allowed[h].resize(W+1);
            int free=n-__builtin_popcount(unsigned(c0|c1));
            int w=mass(c0|c1),w0=mass(c0),w1=mass(c1);
            auto fits=[&](int z){return z>=0&&z<=2*free&&!(z&1);};
            for(int j=start;j<W;++j)if(allowed[h-1][j]) {
                bool ok=j&1 ? (!(c0&c1)&&(fits(j-w0)||fits(j-w1)))
                            : (fits(j)||fits(j-w));
                allowed[h][j]=ok;sum+=ok;
            }
            int q0=c0,q1=c1;
            for(int i=0;i<n;++i){
                if(c0>>i&1){q1|=1<<A.edges[i][0]; (k==1?q1:q0)|=1<<A.edges[i][1];}
                if(c1>>i&1){q0|=1<<A.edges[i][0]; (k==1?q0:q1)|=1<<A.edges[i][1];}
            }
            c0=q0;c1=q1;
        }
    }
};

void check_all_sets(const Machine& A,int tau,int rt) {
    const int n=A.n,k=A.k,N=A.N,U=N-1;
    vector<array<int,2>> pre(N);vector<int> w(n,2), mass(N);w[0]=k==2?3:1;w[1]=1;
    vector<int> incoming(n);
    for(int q=0;q<n;++q)for(int p:A.edges[q])incoming[p]+=w[q];
    for(int q=0;q<n;++q)req(incoming[q]==2*w[q],"stationarity");
    for(int s=1;s<N;++s) {
        int q=__builtin_ctz(unsigned(s));mass[s]=mass[s&(s-1)]+w[q];
    }
    for(int s=0;s<N;++s)for(int q=0;q<n;++q)for(int l=0;l<2;++l)
        if(s>>A.edges[q][l]&1)pre[s][l]|=1<<q;
    vector<int> cost(N), effects(N);
    for(int s=1;s<U;++s) {
        vector<int> q={s}, dist(N,-1);dist[s]=0;
        int gain=-1, change=-1, target=-1;
        for(unsigned i=0;i<q.size() && gain<0;++i) {
            int p=q[i];
            if(change<0&&mass[p]!=mass[s])change=dist[p];
            if(mass[p]>mass[s]){gain=dist[p];target=p;break;}
            for(int t:pre[p])if(dist[t]<0){dist[t]=dist[p]+1;q.push_back(t);}
        }
        req(change>=1&&change==gain&&gain<=n-1,"shortest change equals gain by averaging");
        cost[s]=gain;effects[s]=target;
        if(mass[s]&1){
            req(gain<=tau+1,"odd stage bounded by cover collision");
            if(gain>=2)req(((s>>1)&1)==((s>>2)&1),"odd labels of y,s agree");
        } else {
            if(gain>=2){
                int e=s&7;req(e==0||e==7,"even stage aligned");
                int sz=__builtin_popcount(unsigned(s));
                req(gain<=(e?sz-1:n-sz-1),"aligned-stage bound");
            }
        }
    }
    for(int s=1;s<U;++s)req(cost[s]==cost[U^s],"complement shortest change");
    int E=0, start=k==2?3:2;
    Profile p(A);
    vector<set<int>> odd(n);
    for(int s=1;s<U;++s)if(mass[s]>=start) {
        int m=mass[s],t=cost[s];
        for(int h=1;h<=t;++h)req(p.allowed[h][m],"instance profile contains all actual slow masses");
        if(m&1)for(int h=1;h<=t;++h)odd[h].insert(m);
        else req(t<=even_cost(n,k,m/2),"even cost on all qualifying masses");
    }
    for(int h=1;h<n;++h)req(int(odd[h].size())<=count_odd(n,k,h),"number of slow odd masses");
    for(int s=k==2?1:4;s!=U;s=effects[s])E+=cost[s];
    req(rt<=E && E<=p.sum && p.sum<=bound(n,k,tau),"stationary preimage path and profile are within bound");
}

vector<int> free_targets(int n,int k) {
    vector<int> result(boundary[k].begin(),boundary[k].end());
    for(int q=3;q<n;++q) {result.push_back(q);result.push_back(q);}
    return result;
}
void small() {
    for(int n=4;n<=6;++n)for(int k=0;k<3;++k) {
        auto tail=free_targets(n,k);
        unsigned count=0,strong=0,sync=0,infinite=0,nofinite=0,overlap=0;
        int max_tau=0,max_reset=0;
        do {
            ++count;
            Machine A(n,k,tail);
            if(!A.strong())continue;
            ++strong;
            int r=A.reset(),t=A.collision();
            if(r<0) {if(t<0)++infinite;else ++nofinite;continue;}
            ++sync;req(t>=1,"all synchronizing automata have finite collision");
            check_all_sets(A,t,r);
            max_tau=max(t,max_tau);max_reset=max(r,max_reset);
            auto others=A.degree_two();
            req(others[0]>=1,"the stipulated witness is always found");
            if(others[1] || A.degree_one())++overlap;
        } while(next_permutation(tail.begin(),tail.end()));
        cout<<"n="<<n<<", "<<char('A'+k)<<": "<<count<<" maps, "<<strong<<" strong, "<<sync
            <<" sync, "<<infinite<<" nonsync with no collision, "<<nofinite
            <<" nonsync with collision, "<<overlap<<" sync with other Goel witness; max tau,rt "
            <<max_tau<<','<<max_reset<<".\n";
    }
}
void larger() {
    mt19937 g(428498);
    for(int k=0;k<3;++k)for(int n:{13,17}) {
        auto tail=free_targets(n,k);bool found=false;
        for(int trial=0;trial<10000;++trial) {
            shuffle(tail.begin(),tail.end(),g);
            Machine A(n,k,tail);
            if(!A.strong())continue;
            int t=A.collision();
            if(t<1||bound(n,k,t)>(n-1)*(n-1))continue;
            if(A.one_letter_cycles(0)<2||A.one_letter_cycles(1)<2)continue;
            auto d=A.degree_two();
            if(d[1]||A.degree_one()||!A.uniform_letters_for_all_one_exception_vectors())continue;
            auto rank=A.letter_ranks();
            if(max(rank[0],rank[1])>=n-1)continue;
            int rt=A.reset();
            if(rt<0)continue;
            found=true;
            cout<<"Beyond previous finite range: n="<<n<<", "<<char('A'+k)<<", tau="<<t<<", rt="<<rt
                <<", B="<<bound(n,k,t)<<", instance profile="<<Profile(A).sum
                <<", excluded,other degree-two witnesses="<<d[0]<<','<<d[1]
                <<", one-letter cycles="<<A.one_letter_cycles(0)<<','<<A.one_letter_cycles(1)
                <<", one-letter ranks="<<rank[0]<<','<<rank[1]<<":";
            for(int l=0;l<2;++l) {
                cout<<" [";
                for(int q=0;q<n;++q)cout<<(q?",":"")<<A.edges[q][l];
                cout<<"]";
            }
            cout<<'\n';
            break;
        }
        req(found,"find an early-collision witness with no other Goel class or one-cluster letter");
    }
}
void late_control() {
    int n=13; vector<int> a(n),b(n);
    iota(b.begin(),b.end(),0);
    for(int q=2;q<n;++q)a[q]=q+1;
    a[0]=2;a[1]=0;a[n-1]=1;b[1]=2;
    swap(a[n-5],a[n-3]);swap(b[n-5],b[n-3]);
    vector<int> tail;
    for(int q=2;q<n;++q){tail.push_back(a[q]);tail.push_back(b[q]);}
    Machine A(n,0,tail);
    int t=A.collision(),r=A.reset(),f=Profile(A).sum;
    req(A.strong()&&t==9&&r==56&&f==164&&f>(n-1)*(n-1),"late-collision scope control");
    cout<<"Late control: n=13, A, tau="<<t<<", rt="<<r
        <<", B="<<bound(n,0,t)<<", instance profile="<<f<<" exceeds 144.\n";
}
int main(){
    small();larger();late_control();
    array<array<int,3>,5> expected={{{13,3,4},{20,5,6},{30,8,9},{100,29,29},{1000,292,293}}};
    for(auto e:expected) {
        int n=e[0];
        array<int,2> cut={0,0};
        for(int i=0;i<2;++i) {
            int k=2*i;
            for(int t=1;t<n;++t)if(bound(n,k,t)<=(n-1)*(n-1))cut[i]=t;
        }
        req(cut[0]==e[1]&&cut[1]==e[2],"stated finite cutoff");
        cout<<"n="<<n<<", exact cutoff A/B,C "<<cut[0]<<','<<cut[1]<<", uniform "
            <<bound(n,0,n-1)<<','<<bound(n,2,n-1)<<".\n";
    }
    cout<<"PASS\n";
}
