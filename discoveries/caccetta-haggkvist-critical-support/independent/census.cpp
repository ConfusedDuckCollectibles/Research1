// Independent labelled scan. Add vertices, not pairs; check all cyclic triples.
// Exact determinant screening by integer fraction-free elimination.
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
using namespace std;
struct Audit {
  int n;
  array<int,6> outgoing{};
  array<long long,7> all{}, sourcefree{}, singular{};
  long long determinant() {
    long long m[6][6]{};
    for (int i=0;i<n;++i) for(int j=0;j<n;++j)
      if(i!=j) m[i][j]=(outgoing[j]>>i)&1 ? +1 : ((outgoing[i]>>j)&1 ? 0 : -1);
    long long divisor=1, sign=1;
    for(int i=0;i+1<n;++i){
      int p=i;
      while(p<n&&!m[p][i])++p;
      if(p==n)return 0;
      if(p!=i){swap(m[i],m[p]);sign=-sign;}
      const long long pivot=m[i][i];
      for(int j=i+1;j<n;++j)for(int k=i+1;k<n;++k){
        long long x=pivot*m[j][k]-m[j][i]*m[i][k];
        assert(x%divisor==0);
        m[j][k]=x/divisor;
      }
      for(int j=i+1;j<n;++j)m[j][i]=0;
      divisor=pivot;
    }
    return sign*m[n-1][n-1];
  }
  void dfs(int v){
    assert(0<=v && v<=n && v<=6);
    ++all[v];
    if(v==n){
      int columns=0;
      for(int i=0;i<n;++i)columns|=outgoing[i];
      if(columns!=(1<<n)-1)return;
      ++sourcefree[n];
      if(determinant()!=0)return;
      ++singular[n];
      for(int i=0;i<n;++i)cout<<(i?" ":"")<<outgoing[i];
      cout<<'\n';
      return;
    }
    // In masks L and R we choose exactly the predecessors and successors of v.
    // L,R are disjoint; untouched vertices are nonadjacent to v.
    assert(v<6);
    const int previous=(1<<v)-1;
    for(int R=0;R<=previous;++R){
      int reachable=0;
      for(int i=0;i<v;++i)if((R>>i)&1)reachable|=outgoing[i];
      const int free=previous&~(R|reachable);
      int L=free;
      do {
        outgoing[v]=R;
        for(int i=0;i<v;++i)if((L>>i)&1)outgoing[i]|=1<<v;
        dfs(v+1);
        for(int i=0;i<v;++i)outgoing[i]&=~(1<<v);
        outgoing[v]=0;
        if(!L)break;
        L=(L-1)&free;
      } while(true);
    }
  }
};
int main(int argc, char**argv){
  assert(argc==2);
  Audit a{};
  a.n=atoi(argv[1]);
  assert(1<=a.n && a.n<=6);
  a.dfs(0);
  cerr<<"order "<<a.n<<": "<<a.all[a.n]<<" triangle-free, "
      <<a.sourcefree[a.n]<<" source-free, "
      <<a.singular[a.n]<<" singular\n";
}
