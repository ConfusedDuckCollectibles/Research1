#include <array>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
constexpr int p=13;
int m(int x){return (x%p+p)%p;}
using Seq=array<int,p>;
using Trip=array<int,3>;
using Buckets=array<vector<Seq>,p>;
using Labels=array<array<vector<Seq>,p>,p>;
using Hist=array<array<long long,14>,14>;

void halves(int start,int past,const Trip&T,Seq&w,int a,int b,Labels&R){
  if(start==past){R[m(a)][m(b)].push_back(w);return;}
  for(int v:T){w[start]=v;halves(start+1,past,T,w,a+v,b+start*v,R);}
}
int support(const Seq&v){
  unsigned b=0;for(int x:v)b|=(1u<<x);
  return popcount(b);
}
bool minRotate(const Seq&v){
  Seq r=v;for(int h=1;h<p;h++){
    rotate(r.begin(),r.begin()+1,r.end());
    if(r<v)return false;
  }
  return true;
}
pair<Buckets,int> generate(bool affine,bool rotation){
  Buckets ans;set<Trip> visited;int nShapes=0;
  for(int a=0;a<p;a++)for(int b=a+1;b<p;b++)for(int c=b+1;c<p;c++){
    Trip T={a,b,c};if(visited.contains(T))continue;
    ++nShapes;
    for(int sc=1;sc<=(affine?p-1:1);sc++)for(int delta=0;delta<p;delta++){
      Trip V={m(sc*a+delta),m(sc*b+delta),m(sc*c+delta)};
      sort(V.begin(),V.end());visited.insert(V);
    }
    Labels L,R;Seq w{};
    halves(0,6,T,w,0,0,L);w.fill(0);halves(6,p,T,w,0,0,R);
    for(int j=0;j<p;j++)for(int k=0;k<p;k++)for(const Seq&l:L[j][k])
      for(const Seq&r:R[m(-j)][m(-k)]){
        Seq v;
        copy(l.begin(),l.begin()+6,v.begin());
        copy(r.begin()+6,r.end(),v.begin()+6);
        if(support(v)!=3 || (rotation&&!minRotate(v)))continue;
        int q=0;for(int x=0;x<p;x++)q+=x*x*v[x];
        ans[m(q)].push_back(v);
      }
  }
  if(visited.size()!=286)abort();
  return {move(ans),nShapes};
}
int image(const Seq&f,const Seq&g,int h){
  unsigned w=0;
  for(int x=0;x<p;x++){w|=1u<<m(f[x]-g[(x-h+p)%p]);if(popcount(w)>3)return 4;}
  return popcount(w);
}
int maxLine(const Seq&v){
  int best=1;
  for(int x=0;x<p;x++){
    array<int,13> slope{};
    for(int y=x+1;y<p;y++){
      int s=0;while(m(s*(y-x))!=m(v[y]-v[x]))s++;
      best=max(best,1+(++slope[s]));
    }
  }
  return best;
}
int main(){
  auto [A,na]=generate(true,true);auto [B,nb]=generate(false,false);
  cout<<"shapes "<<na<<" "<<nb<<"\n";
  int64_t pairs=0,step2=0,step3=0,good=0,cap=0,constant=0;Hist hist{};
  for(int j=0;j<p;j++){
    cout<<"bucket "<<j<<" "<<A[j].size()<<" "<<B[j].size()<<"\n";
    for(const Seq&a:A[j])for(const Seq&b:B[j]){
      pairs++; Seq f,g;g[0]=0;
      for(int x=1;x<p;x++)g[x]=m(g[x-1]+b[x]-a[x]);
      for(int x=0;x<p;x++)f[x]=m(g[x]+a[x]);
      if(image(f,g,2)!=3)continue;
      ++step2;
      if(image(f,g,3)!=3)continue;
      ++step3;
      if(image(f,g,0)!=3||image(f,g,1)!=3||m(g[0]-g[p-1]-b[0]+a[0]))abort();
      int sf=0,sg=0,mf=0,mg=0;
      for(int x=0;x<p;x++){sf+=f[x];sg+=g[x];mf+=x*f[x];mg+=x*g[x];}
      if(m(sf)||m(sg)||m(mf)||m(mg))abort();
      int n=4;
      for(int s=4;s<p;s++)n+=(image(f,g,s)<=3);
      if(n<11)continue;
      good++;
      for(int h=4;h<p;h++)constant+=(image(f,g,h)==1);
      int k=maxLine(f),l=maxLine(g);if(k>l)swap(k,l);
      hist[k][l]++;
      cap+=(l<=4);
    }
  }
  cout<<"counts "<<pairs<<" "<<step2<<" "<<step3<<" "<<good<<" "<<cap<<"\n";
  cout<<"constant image occurrences "<<constant<<"\n";
  for(int i=0;i<=p;i++)for(int j=i;j<=p;j++)if(hist[i][j])
    cout<<"max-lines "<<i<<" "<<j<<" : "<<hist[i][j]<<"\n";
  if(na!=3||nb!=22||pairs!=49895820||step2!=70806||step3!=11479||good!=5418||cap)abort();
}
