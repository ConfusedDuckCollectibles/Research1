// Second, deliberately simpler completion checker: materialize all c-edge
// graphs on eight vertices and test exact eight-restrictions only at leaves.
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
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
using Degrees=array<uint8_t,8>;
using Config=vector<int>;
struct Graph{Degrees d{};uint32_t bits=0;};
array<vector<Graph>,6> graphs;
array<unordered_set<string>,12> allowed;
vector<int> pairs;
vector<int> omit;
uint8_t restriction[45][1024];
int pop(unsigned x){return popcount(x);}
string key(Config A){
  sort(A.begin(),A.end());
  string result;
  for(int e:A){assert(0<e&&e<256);result.push_back((char)(uint8_t)e);}
  return result;
}
int move_edge(int e,const array<int,8>& p){
  int r=0;
  for(int i=0;i<8;i++)if(e>>i&1)r|=1<<p[i];
  return r;
}
void generate(int at,int left,int target,Graph g){
  if(!left){graphs[target].push_back(g);return;}
  for(int j=at;j+left<=28;j++){
    Graph h=g;h.bits|=1u<<j;
    for(int i=0;i<8;i++)h.d[i]+=pairs[j]>>i&1;
    generate(j+1,left-1,target,h);
  }
}
bool all_restrictions(const Config& A){
  for(int a=0;a<45;a++){
    Config card;
    for(int e:A)if(restriction[a][e])card.push_back(restriction[a][e]);
    int m=card.size();
    if(m>=12||(m>=9&&!allowed[m].contains(key(card))))return false;
  }
  return true;
}
int main(int argc,char**argv){
  assert(argc==3);
  for(int s=0;s<256;s++)if(pop(s)==2)pairs.push_back(s);
  assert(pairs.size()==28);
  for(int i=0;i<10;i++)for(int j=i+1;j<10;j++)omit.push_back((1<<i)|(1<<j));
  assert(omit.size()==45);
  for(int a=0;a<45;a++)for(int e=0;e<1024;e++)if(pop(e)==4&&!(e&omit[a])){
    int r=0,j=0;
    for(int i=0;i<10;i++)if(!(omit[a]>>i&1)){r|=(e>>i&1)<<j;++j;}
    restriction[a][e]=(uint8_t)r;
  }
  for(int c=3;c<=5;c++)generate(0,c,c,Graph{});
  assert(graphs[3].size()==3276&&graphs[4].size()==20475&&graphs[5].size()==98280);
  vector<Config> cat;
  ifstream f(argv[1]);assert(f);
  int m;array<int,12> counts{};
  while(f>>m){
    assert(m>=9&&m<=11);
    Config C(m);set<int> seen;
    for(int& e:C){assert(bool(f>>e));assert(0<e&&e<256&&pop(e)==4&&seen.insert(e).second);}
    cat.push_back(C);++counts[m];
  }
  assert(f.eof()&&counts[9]==52&&counts[10]==13&&counts[11]==3);
  array<int,8> p{0,1,2,3,4,5,6,7};
  do{
    for(const auto& C:cat){
      Config A;
      for(int e:C)A.push_back(move_edge(e,p));
      allowed[C.size()].insert(key(A));
    }
  }while(next_permutation(p.begin(),p.end()));
  assert(allowed[9].size()==389480&&allowed[10].size()==68460&&allowed[11].size()==10080);
  ifstream l(argv[2]);assert(l);
  uint64_t glues=0,cut_degree=0,cut_cards=0,graphs_tested=0,graphs_degree=0;
  int k;
  while(l>>m>>k){
    assert(m>=9&&m<=11&&k>0);
    Config C(m);set<int> seen;
    for(int& e:C){assert(bool(l>>e));assert(0<e&&e<256&&pop(e)==4&&seen.insert(e).second);}
    assert(allowed[m].contains(key(C)));
    vector<Config> L(k);set<Config> unique;
    for(auto& A:L){
      int r;assert(bool(l>>r)&&r==15-m);A.resize(r);
      set<int> z;
      for(int& e:A){assert(bool(l>>e));assert(0<e&&e<256&&pop(e)==3&&z.insert(e).second);}
      sort(A.begin(),A.end());assert(unique.insert(A).second);
    }
    int c=m-6;
    for(int i=0;i<k;i++)for(int j=i;j<k;j++){
      ++glues;
      Config A=C;
      for(int e:L[i])A.push_back(e|256);
      for(int e:L[j])A.push_back(e|512);
      assert(A.size()==24-c);
      Degrees d{};
      for(int e:A)for(int t=0;t<8;t++)d[t]+=e>>t&1;
      int missing=0;bool root_good=true;
      for(int x:d){if(x+c<9)root_good=false;if(x<9)missing+=9-x;}
      if(!root_good||missing>2*c){++cut_degree;continue;}
      if(!all_restrictions(A)){++cut_cards;continue;}
      for(const auto& g:graphs[c]){
        ++graphs_tested;
        bool good=true;
        for(int t=0;t<8;t++)if(d[t]+g.d[t]<9){good=false;break;}
        if(!good)continue;
        ++graphs_degree;
        Config candidate=A;
        for(int t=0;t<28;t++)if(g.bits>>t&1)candidate.push_back(pairs[t]|768);
        assert(candidate.size()==24);
        if(all_restrictions(candidate)){
          cerr<<"SURVIVOR";
          for(int e:candidate)cerr<<' '<<e;
          cerr<<'\n';
          return 1;
        }
      }
    }
  }
  assert(l.eof());
  cout<<"independent low-low: glues "<<glues<<" root degree rejects "<<cut_degree<<" root eight rejects "<<cut_cards
    <<" graph subsets tested "<<graphs_tested<<" completed degree-compatible "<<graphs_degree<<" survivors 0\n";
  assert(glues==12018&&cut_degree==4009&&cut_cards==7509);
  assert(graphs_tested==23542155&&graphs_degree==41091);
}
