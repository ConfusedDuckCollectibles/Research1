// Independent exact replay: orbit representatives + subset path DP.
// Only synthetic complete digraphs, no external data.
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

using U=std::uint32_t;
struct GroupElt { std::array<U,1024> lo{}, hi{}; };
int factorial(int n) { static constexpr int f[]{1,1,2,6,24,120}; return f[n]; }

int main() {
 for(int n=3;n<=5;++n) {
  int m=n*(n-1),mid=m/2;
  std::array<int,25> id;
  id.fill(-1);
  std::vector<std::array<int,2>> arcs;
  // Use destination-major indexing (independent of the source runners).
  for(int v=0;v<n;++v)for(int u=0;u<n;++u)if(u!=v) {
    id[5*u+v]=(int)arcs.size();
    arcs.push_back({u,v});
  }
  std::array<int,5> p{};
  std::iota(p.begin(),p.end(),0);
  std::vector<GroupElt> group;
  do {
   GroupElt g;
   for(U j=1;j<1024;++j) {
    int bit=std::countr_zero(j);
    U old=j&(j-1);
    if(bit<m) {
      auto [a,b]=arcs[bit];
      g.lo[j]=g.lo[old]|(1u<<id[5*p[a]+p[b]]);
    }
    if(bit+10<m) {
      auto [a,b]=arcs[bit+10];
      g.hi[j]=g.hi[old]|(1u<<id[5*p[a]+p[b]]);
    }
   }
   group.push_back(g);
  }while(std::next_permutation(p.begin(),p.begin()+n));
  assert((int)group.size()==factorial(n));
  U all=1u<<m;
  std::vector<bool> seen(all);
  std::vector<int> min_score(m+1,std::numeric_limits<int>::max());
  std::vector<U> min_mask(m+1), census(m+1), n_orbits(m+1);
  for(U mask=0;mask<all;++mask) {
   if(seen[mask])continue;
   const int k=std::popcount(mask);
   ++n_orbits[k];
   for(auto const &g:group) {
     U moved=g.lo[mask&1023]|g.hi[mask>>10];
     assert(moved>=mask && std::popcount(moved)==k);
     if(!seen[moved]) { seen[moved]=true; ++census[k]; }
   }
   U out[5]{};
   for(int e=0;e<m;++e)if(mask&(1u<<e))out[arcs[e][0]]|=1u<<arcs[e][1];
   // dp[visited subset][final vertex] counts ordered simple paths, including
   // length zero. Each fits into exactly (n-|visited|)! full chain orders.
   int dp[32][5]{};
   for(int i=0;i<n;++i)dp[1<<i][i]=1;
   int score=0;
   for(U visited=1;visited<(1u<<n);++visited)
    for(int u=0;u<n;++u)if(dp[visited][u]) {
      score+=dp[visited][u]*factorial(n-std::popcount(visited));
      U available=out[u]&~visited;
      while(available) {
        U bit=available&-available;
        dp[visited|bit][std::countr_zero(bit)] += dp[visited][u];
        available-=bit;
      }
    }
   if(score<min_score[k]) {min_score[k]=score;min_mask[k]=mask;}
  }
  unsigned choose=1;
  for(int k=0;k<=m;++k) {
   assert(census[k]==choose);
   if(k<m) choose=choose*(m-k)/(k+1);
  }
  U n_total=std::accumulate(n_orbits.begin(),n_orbits.end(),0u);
  assert(n_total==(n==3?16:n==4?218:9608));
  // Integral supporting slopes at mid guarantee minima cannot be beaten
  // by a law mixing different mask cardinalities.
  const int slope=(n==3?2:n==4?4:14);
  for(int k=0;k<=m;++k)assert(min_score[k]>=min_score[mid]+slope*(k-mid));
  assert(min_score[mid]==(n==3?10:n==4?40:204));
  std::cout<<"n="<<n<<" raw="<<all<<" permutation orbits="<<n_total
           <<" supporting affine slope="<<slope
           <<" min uniform midpoint="<<min_score[mid]<<'/'<<factorial(n)
           <<"\nMinima by k:";
  for(auto s:min_score)std::cout<<' '<<s;
  std::cout<<"\nExact midpoint mask:";
  for(int e=0;e<m;++e)if(min_mask[mid]&(1u<<e))std::cout<<' '<<arcs[e][0]<<'>'<<arcs[e][1];
  std::cout<<"\n";
 }
}
