#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
#include <algorithm>
#include <bit>
using namespace std;
constexpr int p=13,n=p*p;
using Mask=uint32_t;
array<array<int,n>,n> direction;
int idx(int a,int b){return a*p+b;}
void setup(){
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){
        int x=(p+i/p-j/p)%p,y=(p+i%p-j%p)%p;
        int m=p;
        if(x) for(m=0;m<p && (m*x)%p!=y;m++){}
        assert(!x || m<p);
        direction[i][j]=m;
    }
}
struct Counter{
    int goal,budget,cap;
    vector<int> selected;
    array<uint64_t,11> nodes{}; array<uint64_t,p+2> masks{};
    vector<int> first;
    void walk(int next, Mask used){
        int k=int(selected.size());
        ++nodes[k];
        if(k==goal){
            ++masks[popcount(used)];
            if(first.empty()) first=selected;
            return;
        }
        for(int i=next;i<=n-(goal-k);i++){
            if(find(selected.begin(),selected.end(),i)!=selected.end())continue;
            Mask now=used;
            array<int,p+1> in_line{};
            bool ok=true;
            for(auto j:selected){
                int dir=direction[i][j];
                now|=Mask(1)<<dir;
                if(++in_line[dir]>=cap || popcount(now)>budget){ok=false;break;}
            }
            if(ok){selected.push_back(i); walk(i+1,now); selected.pop_back();}
        }
    }
    uint64_t run(){
        selected={idx(0,0),idx(1,2),idx(3,1)}; // one independent affine normalization
        Mask used=0;
        for(int i=0;i<3;i++)for(int j=0;j<i;j++)used|=Mask(1)<<direction[selected[i]][selected[j]];
        assert(popcount(used)==3);
        walk(0,used);
        cout<<"points "<<goal<<" cap "<<cap<<" budget "<<budget<<": states";
        for(int k=3;k<=goal;k++)cout<<" "<<nodes[k];
        cout<<" ; survivors "<<nodes[goal];
        for(int j=0;j<=p+1;j++) if(masks[j]) cout<<" ["<<j<<" directions: "<<masks[j]<<"]";
        if(!first.empty()){
            cout<<" ; one witness";
            for(auto i:first)cout<<" ("<<i/p<<","<<i%p<<")";
        }
        cout<<endl;
        return nodes[goal];
    }
};
int main(){
    setup();
    assert((Counter{8,5,5}.run()==0));
    assert((Counter{8,7,3}.run()==0));
    assert((Counter{9,7,4}.run()==0));
    // Independently confirm the actual-language positive neighbors, not merely
    // running a test case on an empty recursion.
    assert((Counter{8,6,5}.run()==72));
    assert((Counter{8,7,4}.run()==9));
}
