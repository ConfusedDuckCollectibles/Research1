#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
using namespace std;
using U = uint64_t;
constexpr int MAXC=400;
struct Row {array<int,6> col; array<int,4> marks;};
vector<Row> rows;
vector<U> inc[MAXC], alive[MAXC], blocked;
bool required[MAXC];
int n, L, words, omitted, ncol;
long long nodes=0, cap;

bool capped=false;
int select(int depth,U* allowed, vector<int>& history){
 if(depth==n)return 1;
 if(++nodes>=cap){capped=true;return 0;}
 int best=-1, small=10000000;
 for(int x=1;x<=L;++x)if(required[x]){
   int have=0;const auto& in=inc[x];
   for(int j=0;j<words;++j)have+=__builtin_popcountll(in[j]&allowed[j]);
   if(!have)return 0;
   if(have<small){small=have;best=x;if(small==1)break;}
 }
 assert(best>0);
 for(int j=0;j<words;++j){
   U m=inc[best][j]&allowed[j];
   while(m){
    int r=64*j+__builtin_ctzll(m);m&=m-1;
    const auto& c=rows[r].col;
    for(int x:c)required[x]=false;
    vector<U>& v=alive[depth];
    for(int k=0;k<words;++k){
      U ban=0;
      for(int x:c)ban|=inc[x][k];
      v[k]=allowed[k]&~ban;
    }
    history.push_back(r);
    if(select(depth+1,v.data(),history))return 1;
    history.pop_back();
    for(int x:c)required[x]=true;
    if(capped)return 0;
   }
 }
 return 0;
}
void init(int h, int seed){
 omitted=h;rows.clear();
 for(int a=1;a<=2*n;++a)
 for(int b=a+1;b<=2*n;++b)
 for(int s=2*n+1;s+a+b<=L;++s)
   rows.push_back({{a,b,s,s+a,s+b,s+a+b},{0,a,a+s,a+s+b}});
 if(seed)shuffle(rows.begin(),rows.end(),mt19937(seed));
 words=(rows.size()+63)/64;
 for(int x=1;x<=L;++x){
   required[x]=(x!=h);inc[x].assign(words,0);
 }
 for(int d=0;d<n;++d)alive[d].resize(words);
 for(int r=0;r<int(rows.size());++r)for(int x:rows[r].col)inc[x][r/64]|=1ULL<<(r%64);
}
int main(int argc,char**argv){
 assert(argc>=2);
 n=atoi(argv[1]);assert(2<=n&&n<=50);
 L=6*n;
 cap=(argc>2?atoll(argv[2]):1000000);
 int seed=(argc>3?atoi(argv[3]):0);
 int specified=(argc>4?atoi(argv[4]):-1);
 for(int h=0;h<1;++h){
  if(specified>=0&&h!=specified)continue;
  init(h,seed);
  vector<int> path;vector<U> base(words,~0ULL);
  nodes=0;capped=false;
  int yes=select(0,base.data(),path);
  cout<<"n="<<n<<" omitted="<<h<<" support_rows="<<rows.size()<<" nodes="<<nodes<<" "<<(yes?"FOUND":capped?"CUTOFF":"EXHAUSTED")<<endl;
  if(yes){
   for(int r:path){for(int x:rows[r].marks)cout<<x<<' ';cout<<endl;}
   return 0;
  }
  if(!(n%2))break;
 }
}
