// Independent normalization: two missing directions become the two axes.
// A non-axis ordered pair becomes (0,0),(1,1); each remaining x is used ≤ once.
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;
const int p=13;
int inv[p];
int mod(int n) {return (n%p+p)%p;}
int slope(int x,int y,int z,int w) {assert(x!=z);return mod((y-w)*inv[mod(x-z)]);}
void gate(int goal, int directions, int cap) {
    array<long long,p+1> visits{};
    int masks=0;
    for (unsigned mask=0;mask<(1u<<p);++mask) {
        if ((mask&3u)!=2u || __builtin_popcount(mask)!=directions) continue;
        ++masks;
        vector<pair<int,int>> a={{0,0},{1,1}};
        auto search = [&](auto&& self,int next_x,unsigned used) -> void {
            ++visits[a.size()];
            if ((int)a.size()==goal) return;
            for (int x=next_x;x<=p-(goal-(int)a.size());++x) {
                for (int y=2;y<p;++y) if (!(used>>y&1u)) {
                    array<int,p> count{};
                    bool ok=true;
                    for (auto [z,w]:a) {
                        int d=slope(x,y,z,w);
                        if (!(mask>>d&1u) || ++count[d]>=cap) {ok=false;break;}
                    }
                    if (ok) {
                        a.emplace_back(x,y);
                        self(self,x+1,used|(1u<<y));
                        a.pop_back();
                    }
                }
            }
        };
        search(search,2,3);
    }
    cout<<goal<<' '<<directions<<' '<<cap<<" masks="<<masks<<" visits=";
    for (int i=2;i<=goal;++i) cout<< (i==2?"":",")<<visits[i];
    cout<<'\n';
    assert(visits[goal]==0);
}
int main() {
    for (int i=1;i<p;++i) for (int j=1;j<p;++j) if (i*j%p==1) inv[i]=j;
    gate(8,5,5);
    gate(8,7,3);
    gate(9,7,4);
}
