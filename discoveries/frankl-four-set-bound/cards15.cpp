#ifdef NDEBUG
#error "Build with assertions enabled"
#endif
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <vector>
using namespace std;
using W=__uint128_t;
struct Hash {size_t operator()(W x)const{return (uint64_t)x ^ (uint64_t)(x>>64)*0x9e3779b97f4a7c15ULL;}};
int idx[256], c8[512][9], maybe[57][8];
vector<int> tri;
unordered_set<W,Hash> allowed[12];
long long visited,pass,cut8,cutdeg;
vector<int> picked,base;
int del(int x,int omit){int r=0,at=0; for(int i=0;i<9;i++)if(!(omit>>i&1)){r|=(x>>i&1)<<at; ++at;}return r;}
int pc(W x){return popcount((uint64_t)x)+popcount((uint64_t)(x>>64));}
bool good8(W x){int m=pc(x);return m<9 || (m<12&&allowed[m].contains(x));}
void dfs(int first,int need,int mind,const array<int,8>&d,const array<W,8>&eight){
  ++visited;
  for(int i=0;i<8;i++){
    int has=maybe[first][i], nothas=56-first-has;
    if(d[i]+min(need,has)<mind || d[i]+max(0,need-nothas)>15){++cutdeg;return;}
  }
  if(!need){
    ++pass;
    for(int e:base)cout<<e<<' ';for(int e:picked)cout<<e<<' ';cout<<'\n';
    return;
  }
  for(int j=first;j<=56-need;j++){
    int e=tri[j]|256;
    auto x=eight;auto z=d;
    bool ok=true;
    for(int i=0;i<8;i++){
      if(e>>i&1)++z[i];
      else if(!good8(x[i]|=W(1)<<c8[e][i])){++cut8;ok=false;break;}
    }
    if(!ok)continue;
    picked.push_back(e);dfs(j+1,need-1,mind,z,x);picked.pop_back();
  }
}
int main(int argc,char**argv){
  assert(argc==3);
  int k=0;for(int v=0;v<256;v++)if(popcount((unsigned)v)==4)idx[v]=k++;
  assert(k==70);
  for(int e=0;e<512;e++)for(int i=0;i<9;i++)if(!(e>>i&1))c8[e][i]=idx[del(e,1<<i)];
  for(int e=0;e<256;e++)if(popcount((unsigned)e)==3)tri.push_back(e);
  assert(tri.size()==56);
  for(int j=55;j>=0;j--)for(int i=0;i<8;i++)maybe[j][i]=maybe[j+1][i]+(tri[j]>>i&1);
  ifstream f(argv[1]);assert(f);vector<vector<int>> catalog;int m;array<int,12> count{};
  while(f>>m){assert(m>=9&&m<=11);vector<int>a(m);unordered_set<int> seen;for(int&e:a){assert(bool(f>>e));assert(0<=e&&e<256&&popcount((unsigned)e)==4&&seen.insert(e).second);}catalog.push_back(a);++count[m];}
  assert(f.eof()&&count[9]==52&&count[10]==13&&count[11]==3);
  array<int,8> p{0,1,2,3,4,5,6,7};
  do{
    int mp[256]={};
    for(int q=1;q<256;q++){int j=countr_zero((unsigned)q);mp[q]=mp[q&(q-1)]|(1<<p[j]);}
    for(const auto&a:catalog){W s=0;for(int e:a)s|=W(1)<<idx[mp[e]];allowed[a.size()].insert(s);}
  }while(next_permutation(p.begin(),p.end()));
  cerr<<"allowed8 labels: "<<allowed[9].size()<<' '<<allowed[10].size()<<' '<<allowed[11].size()<<'\n';
  // Independently over-approximate the 11 layer from the complete 10
  // layer. The three known negative 11 seeds are not assumed exhaustive.
  auto listed=catalog;
  int ext11=0,extra=0;
  vector<int> mask_at(70);
  for(int e=0;e<256;e++)if(popcount((unsigned)e)==4)mask_at[idx[e]]=e;
  for(const auto& b:listed){
    if(b.size()!=10)continue;
    W family=0;for(int e:b)family|=W(1)<<idx[e];
    for(int j=0;j<70;j++){
      if(family>>j&1)continue;
      W child=family|(W(1)<<j);
      bool ok=true;
      for(int i=0;i<70;i++)if((child>>i&1)&&!allowed[10].contains(child^(W(1)<<i))){ok=false;break;}
      if(!ok)continue;
      ++ext11;
      if(allowed[11].contains(child))continue;
      ++extra;
      vector<int> child_edges=b;child_edges.push_back(mask_at[j]);
      catalog.push_back(child_edges);
      do{
        int mp[256]={};
        for(int q=1;q<256;q++){int i=countr_zero((unsigned)q);mp[q]=mp[q&(q-1)]|(1<<p[i]);}
        W key=0;for(int e:child_edges)key|=W(1)<<idx[mp[e]];
        allowed[11].insert(key);
      }while(next_permutation(p.begin(),p.end()));
    }
  }
  cerr<<"boundary augmentation 10->11 "<<ext11<<" labelled hits; new orbits "<<extra<<"; all 11 labels "<<allowed[11].size()<<'\n';
  assert(ext11==28&&extra==0&&allowed[11].size()==10080);
  ofstream out(argv[2]);assert(out);
  for(const auto& b:catalog){out<<b.size();for(int e:b)out<<' '<<e;out<<'\n';}
  for(const auto& b:catalog){ base=b;
    array<int,8>d{};array<W,8>eight{};
    for(int e:base){
      for(int i=0;i<8;i++)if(e>>i&1)++d[i];else eight[i]|=W(1)<<c8[e][i];
    }
    for(auto x:eight)assert(good8(x));
    int t=15-base.size();dfs(0,t,t,d,eight);
  }
  assert(visited==1079242&&pass==242&&cutdeg==142231&&cut8==15634511);
  cerr<<"visited "<<visited<<" candidate labeled cards "<<pass<<" cuts degree/eight "<<cutdeg<<' '<<cut8<<'\n';
}
