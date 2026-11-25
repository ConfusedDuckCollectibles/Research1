// Independent finite-plane assertions for the F13 spectral audit.
// Enumerate normalized caps with direction allowance or with a unique
// displacement per direction. Bitsets preserve all admissible points.
#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>
using namespace std;
constexpr int p=13,n=p*p;
using Bits=bitset<n>;
struct Point { int x,y; };
Point pts[n];
int inv[p],dir[n][n],signed_chord[n][n];
Bits rays[n][p+1],lines[n][n];
constexpr unsigned forced=(1u<<0)|(1u<<(p-1))|(1u<<p);
int mod(int a){int b=a%p;return b<0?b+p:b;}
int id(Point u){return u.x*p+u.y;}
int difference(Point a,Point b){
 int x=mod(a.x-b.x),y=mod(a.y-b.y);
 return min(x*p+y,mod(-x)*p+mod(-y));
}
int slope(Point a,Point b){
 int x=mod(a.x-b.x),y=mod(a.y-b.y);
 return x?mod(y*inv[x]):p;
}
void setup(){
 for(int i=1;i<p;i++)for(int j=1;j<p;j++)if(i*j%p==1)inv[i]=j;
 for(int i=0;i<n;i++)pts[i]={i/p,i%p};
 for(int i=0;i<n;i++)
  for(int j=0;j<n;j++)if(i!=j){
   dir[i][j]=slope(pts[i],pts[j]);
   signed_chord[i][j]=difference(pts[i],pts[j]);
   rays[i][dir[i][j]].set(j);
  }
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j)lines[i][j]=rays[i][dir[i][j]];
}
unsigned transform(unsigned mask,int mode){
 unsigned b=0;
 for(int a=0;a<=p;a++)if(mask&(1u<<a)){
  int x=0;
  if(mode==0)x=a==p?p:(a+1)%p;       // translation
  if(mode==1)x=a==p?p:(2*a)%p;       // primitive scalar
  if(mode==2)x=a==p?0:a==0?p:inv[a]; // inversion
  b|=1u<<x;
 }
 return b;
}
unsigned orbit(unsigned mask){
 set<unsigned> v{mask};
 vector<unsigned> visit{mask};
 for(size_t i=0;i<visit.size();i++)
  for(int g=0;g<3;g++){
   auto w=transform(visit[i],g);
   if(v.insert(w).second)visit.push_back(w);
  }
 return *v.begin();
}
struct Search{
 int goal;
 const vector<Bits>* adj;
 bool mono;
 long long calls=0;
 vector<int> chosen{0,1,p};
 array<int,p+1> slopes{};
 vector<int> certificate;
 Search(int goal,const vector<Bits>* adj,bool mono):goal(goal),adj(adj),mono(mono){slopes.fill(-1);}
 bool back(Bits left, array<int,p+1> chord){
  ++calls;
  if((int)chosen.size()==goal){certificate=chosen;return true;}
  for(int i=0;i<n && left.count()+chosen.size()>=(size_t)goal;i++)
   if(left.test(i)){
    left.reset(i);
    auto c=chord;
    bool ok=true;
    Bits next=left;
    for(auto j:chosen){
     auto s=dir[i][j];
     if(mono && c[s]!=-1 && c[s]!=signed_chord[i][j]){ok=false;break;}
     c[s]=signed_chord[i][j];
     next&=~lines[i][j];
    }
    if(!ok)continue;
    if(adj)next&=(*adj)[i];
    chosen.push_back(i);
    if(back(next,c))return true;
    chosen.pop_back();
   }
  return false;
 }
 bool run(){
  Bits B;B.set();
  for(auto i:chosen){B.reset(i);if(adj)B&=(*adj)[i];}
  for(int i:chosen)for(int j:chosen)if(i<j){
   B&=~lines[i][j];
   slopes[dir[i][j]]=signed_chord[i][j];
  }
  return back(B,slopes);
 }
};
vector<Bits> adjacency(unsigned allowed) {
 vector<Bits>A(n);
 for(int i=0;i<n;i++)for(int j=0;j<=p;j++)if(allowed&(1u<<j))A[i]|=rays[i][j];
 return A;
}
int main(int argc,char**argv) {
 setup();
 for(auto [allowed,sz]: vector<pair<int,int>>{{5,5},{6,6},{8,8}}){
  int hit=0,masks=0;long long calls=0;set<unsigned> reps;
  for(unsigned bits=0;bits<(1u<<(p+1));bits++)
   if((bits&forced)==forced && popcount(bits)==allowed){
    auto a=adjacency(bits);
    Search scan(sz,&a,false);masks++;
    if(scan.run()){
     hit++;
     unsigned realized=0;
     for(int i:scan.certificate)for(int j:scan.certificate)if(i<j)realized|=1u<<dir[i][j];
     assert(popcount(realized)==allowed);
     reps.insert(orbit(realized));
    }
    calls+=scan.calls;
   }
  if(allowed!=6)assert(reps.empty());
  else assert((reps==set<unsigned>{189}));
  cout<<"cap "<<sz<<" max-dir "<<allowed<<" masks "<<masks<<" hits "<<hit<<" calls "<<calls<<" orbits";
  for(auto m:reps)cout<<' '<<m;
  cout<<'\n';
 }
 for(int goal:{7,8}){
  Search scan(goal,nullptr,true);
  bool r=scan.run();
  assert(r==(goal==7));
  cout<<"equal chord cap target "<<goal<<" found "<<r<<" visited "<<scan.calls;
  for(int i:scan.certificate)cout<<" ("<<pts[i].x<<','<<pts[i].y<<')';
  cout<<'\n';
 }
 if(argc>1){
  auto cap=orbit(189);
  assert(cap==189);
  for(int j=1;j<argc;j++){
   unsigned bad=stoul(argv[j]),good=bad^((1u<<(p+1))-1);
   assert(popcount(good)==6);
   cout<<"profile "<<bad<<" good PGL orbit "<<orbit(good)<<'\n';
   assert(orbit(good)!=cap);
  }
 }
}
