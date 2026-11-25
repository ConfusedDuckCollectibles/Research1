// Independent exact F13 planar projection gate for spectral size 26.
// Finite differences generate value vectors; direct Gaussian elimination
// on homogeneous binomial moments enforces geometry. No search cutoff.
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <vector>
using namespace std;
constexpr int p=13;
int mod(int a){int x=a%p; return x<0?x+p:x;}
int inv[p];
int powmod(int x,int n){int y=1;while(n--)y=y*x%p;return y;}
int choose(int n,int j){int x=1;for(int i=1;i<=j;i++)x=x*(n-i+1)/i;return x;}
struct Profile {array<int,p> values{},c{};int top=0,peak=0;};
using DB=vector<Profile>;
DB values(int degree,int height,int total,bool at_most=false) {
 DB out;array<int,p> t{};
 array<array<int,p>,p> N{};N[0][0]=1;
 // N[k](x) = binomial(x,k), computed via k! inverse over F13.
 for(int k=1;k<=degree;k++)
  for(int j=0;j<=k;j++)
   N[k][j]=mod(((j?N[k-1][j-1]:0)-(k-1)*N[k-1][j])*inv[k]);
 function<void(int,int)> ext=[&](int n,int sum){
  if(sum>total)return;
  if(n<=degree) {
   for(int z=0;z<=height;++z) {t[n]=z;ext(n+1,sum+z);}
   return;
  }
  int d=degree+1;
  for(int x=d;x<p;x++) {
   int v=0;
   for(int j=0;j<d;j++)v+=((d-j)%2?-1:1)*choose(d,j)*t[x-d+j];
   t[x]=mod(-v);
   if(t[x]>height)return;
   sum+=t[x];
   if(sum>total)return;
  }
  if(!at_most && sum!=total)return;
  Profile F;F.values=t;F.peak=*max_element(t.begin(),t.end());
  auto u=t;
  for(int k=0;k<=degree;k++){
   for(int j=0;j<=k;j++)F.c[j]=mod(F.c[j]+u[0]*N[k][j]);
   for(int j=0;j<degree-k;j++)u[j]=mod(u[j+1]-u[j]);
  }
  for(int j=1;j<=degree;j++)if(F.c[j])F.top=j;
  if(!F.top)return;
  for(int x=0;x<p;x++) {int r=0;for(int j=F.top;j>=0;j--)r=mod(r*x+F.c[j]);assert(r==t[x]);}
  out.push_back(F);
 };
 ext(0,0);
 return out;
}
using Lists=vector<vector<int>>;
struct Gate {
 int n,b,degree; const DB& db;
 vector<int> slopes;
 map<pair<int,int>,array<array<int,9>,9>> table;
 map<pair<int,int>,long long> histogram;
 long long survivors=0, positive_tuples=0;
 Gate(int n,int b,const DB& d):n(n),b(b),degree(n-2),db(d){}
 auto relations(int j,int free_mask)->const array<array<int,9>,9>&{
  auto key=make_pair(j,free_mask);
  auto it=table.find(key);
  if(it!=table.end())return it->second;
  vector<int> free,dep;
  for(int i=0;i<n;i++)(free_mask&(1<<i)?free:dep).push_back(i);
  assert((int)dep.size()==j+1);
  int A[9][9]{};
  for(int h=0;h<=j;h++){
   for(int k=0;k<(int)dep.size();k++)A[h][k]=powmod(slopes[dep[k]],h);
   for(int k=0;k<(int)free.size();k++)A[h][j+1+k]=mod(-powmod(slopes[free[k]],h));
  }
  for(int k=0;k<=j;k++){
   int r=k;while(r<=j&&!A[r][k])r++;
   assert(r<=j);
   for(int l=0;l<n;l++)swap(A[r][l],A[k][l]);
   int z=inv[A[k][k]];
   for(int l=k;l<n;l++)A[k][l]=mod(A[k][l]*z);
   for(int h=0;h<=j;h++)if(h!=k) {
    int t=A[h][k];
    for(int l=k;l<n;l++)A[h][l]=mod(A[h][l]-t*A[k][l]);
   }
  }
  array<array<int,9>,9> K{};
  for(int k=0;k<(int)free.size();k++){
   K[free[k]][free[k]]=1;
   for(int h=0;h<=j;h++)K[dep[h]][free[k]]=A[h][j+1+k];
  }
  return table.emplace(key,K).first->second;
 }
 void consume(const array<int,9>& selected,int constants) {
  ++positive_tuples;
  if(mod(constants-2*(n-1)))return;
  array<int,3> freqs{};int peak=0;
  for(int x=0;x<p;x++)for(int y=0;y<p;y++) {
   int r=-2*(n-1);
   for(int j=0;j<n;j++)r+=db[selected[j]].values[mod(y-slopes[j]*x)];
   if(r<0 || r>2*p || r%p)return;
   ++freqs[r/p];
  }
  if(freqs[1]+2*freqs[2]!=26)return;
  for(int j=0;j<n;j++)peak=max(peak,db[selected[j]].peak);
  ++survivors;
  ++histogram[{freqs[2],peak}];
 }
 void back(int j,const Lists& lists){
  if(j==0){
   array<int,9> sel{};
   auto fini=[&](auto&& fini,int i,int c)->void {
    if(i==n){consume(sel,c);return;}
    for(int k:lists[i]){sel[i]=k;fini(fini,i+1,c+db[k].c[0]);}
   };
   fini(fini,0,0);return;
  }
  vector<array<vector<int>,p>> partition(n);
  vector<pair<int,int>> order;
  for(int i=0;i<n;i++){
   int k=0;for(int q:lists[i]) {
    auto &bucket=partition[i][db[q].c[j]]; k+=bucket.empty(); bucket.push_back(q);
   }
   order.emplace_back(k,i);
  }
  sort(order.begin(),order.end());
  int useful=n-j-1,m=0;
  for(int i=0;i<useful;i++)m|=1<<order[i].second;
  auto& mat=relations(j,m);
  array<int,9> fre{};
  auto fill=[&](auto&& fill,int k)->void {
   if(k<useful){
    int i=order[k].second;
    for(int z=0;z<p;z++)if(!partition[i][z].empty()){fre[i]=z;fill(fill,k+1);}
    return;
   }
   Lists next(n);
   for(int i=0;i<n;i++){
    int v=0;
    for(int k=0;k<useful;k++) {
     int h=order[k].second;v+=mat[i][h]*fre[h];
    }
    if(partition[i][v%p].empty())return;
    next[i]=partition[i][v%p];
   }
   back(j-1,next);
  };
  fill(fill,0);
 }
 void run(int m) {
  slopes.clear();
  table.clear();histogram.clear();survivors=0;positive_tuples=0;
  for(int s=0;s<p;s++)if(m&(1<<s))slopes.push_back(s);
  Lists ls(n);
  for(int i=0;i<n;i++)for(int k=0;k<(int)db.size();k++)
   if(i>=2 || db[k].c[db[k].top-1]==0)ls[i].push_back(k);
  back(degree,ls);
  if(survivors){cout<<"slope mask "<<m<<" positive "<<survivors;
   for(auto [x,n0]:histogram)cout<<" [doubles "<<x.first<<" peak "<<x.second<<" count "<<n0<<']';
   cout<<'\n';
  }
 }
};
int main(int argc,char**argv){
 for(int i=1;i<p;i++)for(int j=1;j<p;j++)if(i*j%p==1)inv[i]=j;
 if(argc==1) {
  auto low=values(3,p-1,26,true);
  assert(low.empty());
  cout<<"nonconstant degree <=3, canonical mass <=26: NONE\n";
  for(int d=1;d<=3;d++){
   int minsum=1000;
   int radix=1;for(int i=0;i<d;i++)radix*=p;
   for(int lead=1;lead<p;lead++)for(int u=0;u<radix;u++){
    array<int,4> c{};c[d]=lead;
    int q=u;for(int i=0;i<d;i++){c[i]=q%p;q/=p;}
    int sum=0;
    for(int x=0;x<p;x++) {int v=0;for(int i=d;i>=0;i--)v=mod(x*v+c[i]);sum+=v;}
    minsum=min(minsum,sum);
   }
   assert(minsum==(d==1?78:d==2?52:39));
   cout<<"degree "<<d<<" canonical minimum "<<minsum<<'\n';
  }
  return 0;
 }
 if(argc!=3)return 2;
 int n=atoi(argv[1]),b=atoi(argv[2]);
 assert(n>=6 && n<=9 && b>=4 && b<=8);
 auto db=values(n-2,b,26);
 Gate g(n,b,db);
 int count=0;long long total=0,p_tuples=0;
 array<bool,1<<p> seen{};
 for(int m=0;m<(1<<p);m++)if(std::popcount((unsigned)m)==n&&!seen[m]) {
  for(int a=1;a<p;a++)for(int b=0;b<p;b++){
   int q=0;for(int s=0;s<p;s++)if(m&(1<<s))q|=1<<((a*s+b)%p);
   assert(q>=m);seen[q]=true;
  }
  count++;g.run(m);total+=g.survivors;p_tuples+=g.positive_tuples;
 }
 cout<<"L "<<n<<" M "<<b<<" DB "<<db.size()<<" AGL representatives "<<count<<" positive tuples "<<p_tuples<<" RESULT "<<total<<"\n";
}
