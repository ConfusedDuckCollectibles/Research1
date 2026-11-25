// Scoped synthetic exhaustive benchmark. No real kidney data.
#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>
#include <numeric>
#include <limits>
using namespace std;
int fact(int n){int r=1; for(int i=2;i<=n;i++) r*=i; return r;}
void complete(int n){
    int m=n*(n-1);
    unsigned S=1u<<m;
    vector<int> val(S,0);
    vector<int> perm(n); iota(perm.begin(),perm.end(),0);
    vector<vector<int>> idx(n,vector<int>(n,-1));
    int t=0;
    for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(a!=b)idx[a][b]=t++;
    // Every uniform random full permutation is a legal fixed chain.
    // Each successful prefix adds one. The zeta transform counts every
    // prefix mask contained in the actual present-arc mask.
    do{
       unsigned bits=0;
       val[0]++;
       for(int k=1;k<n;k++){
          bits|=1u<<idx[perm[k-1]][perm[k]];
          val[bits]++;
       }
    }while(next_permutation(perm.begin(),perm.end()));
    for(int b=0;b<m;b++)for(unsigned s=0;s<S;s++)if(s&(1u<<b))val[s]+=val[s^(1u<<b)];
    vector<int> mv(m+1,INT32_MAX),counts(m+1,0);
    vector<unsigned> best(m+1,0);
    for(unsigned s=0;s<S;s++){
       int k=__builtin_popcount(s);
       counts[k]++;
       if(val[s]<mv[k]){mv[k]=val[s];best[k]=s;}
    }
    cout<<n<<" "<<m<<" "<<fact(n)<<" "<<S<<"\n";
    for(int k=0;k<=m;k++)cout<<k<<" "<<mv[k]<<" "<<best[k]<<" "<<counts[k]<<"\n";
}
int main(){for(int n=3;n<=5;n++)complete(n);}
