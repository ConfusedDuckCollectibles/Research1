#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <map>
#include <vector>
using namespace std;
constexpr int p=13;
int mod(int x){x%=p;return x<0?x+p:x;}
int inv[p];
struct Profile {
    array<int,9> coefficients;
    array<int,p> values;
    int degree, max_value;
};
vector<Profile> db;
int maxd, bound;
long long histogram[14][14];

// All degree <= d polynomials are uniquely determined by their values at 0,...,d.
// Every value in a tested profile is in [0,bound], with bound < p.
void value_gen(int d) {
 vector<array<int,9>> lag(d+1);
 for(int i=0;i<=d;++i){
    array<int,9> a{};a[0]=1;int den=1;
    for(int j=0;j<=d;++j)if(i!=j){
       den=mod(den*(i-j));
       for(int k=d;k>=0;--k)a[k]=mod((k?a[k-1]:0)-j*a[k]);
    }
    for(int k=0;k<=d;++k)lag[i][k]=mod(a[k]*inv[den]);
 }
 array<int,9> seed{};
 function<void(int,int)> rec=[&](int x,int mass){
    if(mass>26)return;
    if(x<=d){for(int a=0;a<=bound;++a){seed[x]=a;rec(x+1,mass+a);}return;}
    array<int,9>a{};array<int,p> v{};
    int deg=0;
    for(int k=0;k<=d;++k){
       int l=0;for(int i=0;i<=d;++i)l+=lag[i][k]*seed[i];
       a[k]=mod(l);if(a[k])deg=k;
    }
    if(!deg)return;
    int mx=0;
    for(int i=0;i<=d;++i){v[i]=seed[i];mx=max(mx,v[i]);}
    for(int i=d+1;i<p;++i){
       int z=0;for(int k=deg;k>=0;--k)z=mod(z*i+a[k]);
       if(z>bound)return;
       mx=max(mx,z);v[i]=z;mass+=z;if(mass>26)return;
    }
    if(mass==26){
       for(int i=0;i<=d;++i){
           int z=0;for(int k=deg;k>=0;--k)z=mod(z*i+a[k]);
           assert(z==v[i]);
       }
       db.push_back({a,v,deg,mx});
    }
 };
 rec(0,0);
}
int L;
vector<int>s;
vector<int>w;
struct Imat {int val[13][13];};
using Lists=vector<vector<int>>;
long long complete_count, survivors;
int good_mask;
const int SIZE=1<<13;
map<int,Imat> crosses;
void set_slopes(int m){
    s.clear();w.clear();crosses.clear();
    good_mask=((1<<(p+1))-1)^m;
    for(int x=0;x<p;++x)if(m>>x&1)s.push_back(x);
    for(int i=0;i<L;++i){
       int q=1;for(int j=0;j<L;++j)if(i!=j)q=mod(q*(s[i]-s[j]));
       w.push_back(inv[q]);
    }
}
Imat transform(int piv){
    auto it=crosses.find(piv);
    if(it!=crosses.end())return it->second;
    Imat m{};
    for(int j=0;j<L;++j)if(piv>>j&1){
        int den=1;
        for(int k=0;k<L;++k)if(k!=j&&(piv>>k&1))den=mod(den*(s[j]-s[k]));
        for(int i=0;i<L;++i){
            int num=1;
            for(int k=0;k<L;++k)if(k!=j&&(piv>>k&1))num=mod(num*(s[i]-s[k]));
            m.val[i][j]=mod(num*inv[den]*w[i]*inv[w[j]]);
        }
    }
    crosses[piv]=m;
    return m;
}
// The jth coefficients divided by w_i must be evaluations of a
// polynomial of degree <= L-j-2. Use any L-j-1 pivot columns.
void dfs(int j,const Lists& lists){
 if(j==0){
    vector<int> use(L);
    function<void(int,int)> fin=[&](int q,int sm){
       if(q<L){for(auto t:lists[q]){use[q]=t;fin(q+1,sm+db[t].coefficients[0]);}return;}
       ++complete_count;
       if(mod(sm-2*(L-1)))return;
       int doubles=0, mass=0;
       int grid[p][p];
       for(int x=0;x<p;++x)for(int y=0;y<p;++y){
            int val=0;for(int i=0;i<L;++i)val+=db[use[i]].values[mod(y-s[i]*x)];
            int pf=val-2*(L-1);
            if(pf<0||pf>2*p||pf%p)return;
            grid[x][y]=pf/p;
            mass+=grid[x][y];doubles+=(pf==2*p);
       }
       // Independent direct check of each strong survivor, not just
       // interpolation: recalculate all fourteen actual line profiles.
       assert(mass==26);
       for(int m=0;m<=p;++m){
          int g[p]{};
          for(int x=0;x<p;++x)for(int y=0;y<p;++y)
             g[m==p?x:mod(y-m*x)]+=grid[x][y];
          if(good_mask>>m&1) {
             for(int t=0;t<p;++t)assert(g[t]==2);
          } else {
             int i=lower_bound(s.begin(),s.end(),m)-s.begin();
             for(int t=0;t<p;++t)assert(g[t]==db[use[i]].values[t]);
          }
       }
       ++survivors;
       int h=0;for(auto t:use)h=max(h,db[t].max_value);
       assert(doubles<=13);
       ++histogram[doubles][h];
    };
    fin(0,0); return;
 }
 int D=L-j-1;
 vector<array<vector<int>,p>> choices(L);
 vector<pair<int,int>> sz;
 for(int i=0;i<L;++i){
     for(int t:lists[i]) choices[i][db[t].coefficients[j]].push_back(t);
     int a=0;for(const auto&v:choices[i])a+=!v.empty();
     sz.emplace_back(a,i);
 }
 sort(sz.begin(),sz.end());
 int pivot=0;for(int j=0;j<D;++j)pivot|=1<<sz[j].second;
 auto m=transform(pivot);
 vector<int> val(L);
 function<void(int)> walk=[&](int q){
    if(q<D){
       int i=sz[q].second;
       for(int v=0;v<p;++v)if(!choices[i][v].empty()){val[i]=v;walk(q+1);}return;
    }
    Lists v(L);
    for(int i=0;i<L;++i){
       int c=0;for(int j=0;j<D;++j){int z=sz[j].second;c+=m.val[i][z]*val[z];}
       if(choices[i][mod(c)].empty())return;
       v[i]=choices[i][mod(c)];
    }
    dfs(j-1,v);
 };
 walk(0);
}
int orbit_min(int m){
 int ans=m;
 for(int a=1;a<p;++a)for(int b=0;b<p;++b){
    int h=0;
    for(int j=0;j<p;++j)if(m>>j&1)h|=1<<mod(a*j+b);
    ans=min(ans,h);
 }
 return ans;
}
int main(int argc,char**argv){
  assert(argc==3);
  L=atoi(argv[1]);bound=atoi(argv[2]);maxd=L-2;
  assert(L>=6&&L<=9&&bound>=2&&bound<p);
  for(int i=1;i<p;++i)
    for(int j=1;j<p;++j)if(mod(i*j)==1)inv[i]=j;

  value_gen(maxd);
  printf("L %d bound %d profiles %zu\n",L,bound,db.size());
  Lists lists(L);
  // Independent translations of the arguments of the first two profiles
  // set each one's next-to-leading coefficient to zero.
  for(int i=0;i<L;++i)for(int v=0;v<(int)db.size();++v)
    if(i>=2 || !db[v].coefficients[db[v].degree-1])lists[i].push_back(v);
  int num=0;
  long long total_tuples=0,total_survivors=0;
  for(int m=0;m<SIZE;++m)
    if(__builtin_popcount((unsigned)m)==L && orbit_min(m)==m){
      set_slopes(m);
      complete_count=survivors=0;
      dfs(maxd,lists);
      printf("mask %d tuples %lld survivors %lld\n",m,complete_count,survivors);
      total_tuples+=complete_count;
      total_survivors+=survivors;
      ++num;
    }
  printf("TOTAL masks %d tuples %lld survivors %lld\n",num,total_tuples,total_survivors);
  for(int d=0;d<14;++d)for(int b=0;b<14;++b)
    if(histogram[d][b]) printf("doubles %d max %d count %lld\n",d,b,histogram[d][b]);
}
