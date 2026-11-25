// Independent enumeration: directly stop each actual chain on its first absent arc.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
using namespace std;

int main(){
    for(int n=3;n<=5;n++){
        vector<array<unsigned,4>> paths;
        vector<int> ord(n); iota(ord.begin(),ord.end(),0);
        auto edge=[n](int u,int v){return u*(n-1)+v-(v>u);};
        do{
            array<unsigned,4> a{};
            for(int k=0;k<n-1;k++)a[k]=1u<<edge(ord[k],ord[k+1]);
            paths.push_back(a);
        }while(next_permutation(ord.begin(),ord.end()));
        int m=n*(n-1);
        unsigned total=1u<<m;
        vector<int> minimum(m+1,100000),counts(m+1,0);
        vector<unsigned> arg(m+1,0);
        for(unsigned mask=0;mask<total;mask++){
            int count=0;
            for(auto const &path: paths){
                count++;
                for(int k=0;k<n-1;k++){
                    if((path[k]&mask)==0) break;
                    count++;
                }
            }
            unsigned h=mask;int k=0;
            while(h){ h&=h-1; k++; }
            counts[k]++;
            if(count<minimum[k]){minimum[k]=count;arg[k]=mask;}
        }
        cout<<n<<" "<<m<<" "<<paths.size()<<" "<<total<<"\n";
        for(int k=0;k<=m;k++)cout<<k<<" "<<minimum[k]<<" "<<arg[k]<<" "<<counts[k]<<"\n";
    }
}
