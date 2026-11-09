// Conservative first-six sign/event search for 28-move 12-flip paths.
// U1,U3,U5,U2/U4 event-runs undercount their actual number of moves;
// the single U6 is relaxed to any jump preserving first-five minor signs.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;
using A3=array<int,3>;using A4=array<int,4>;
struct Rule{uint32_t m[3];int flip[3];};
struct Edge{uint32_t v;uint8_t j,weight;};
vector<A3>P;int pid[7][7][7];vector<Rule>rules;
int orientid(int a,int b,int c,int &sg){int z[3]={a,b,c};sg=(a>b)^(a>c)^(b>c);sort(z,z+3);return pid[z[0]][z[1]][z[2]];}
void init(){
 for(int a=1;a<=6;a++)for(int b=a+1;b<=6;b++)for(int c=b+1;c<=6;c++){pid[a][b][c]=P.size();P.push_back({a,b,c});}
 assert(P.size()==20);
 for(int s=1;s<=6;s++)for(int a=1;a<=6;a++)for(int b=a+1;b<=6;b++)for(int c=b+1;c<=6;c++)for(int d=c+1;d<=6;d++){
  if(s==a||s==b||s==c||s==d)continue;int pairs[3][4]={{a,b,c,d},{a,c,b,d},{a,d,b,c}};Rule q{};
  for(int z=0;z<3;z++){int x,y;int i=orientid(s,pairs[z][0],pairs[z][1],x),j=orientid(s,pairs[z][2],pairs[z][3],y);q.m[z]=(1u<<i)^(1u<<j);q.flip[z]=x^y^(z==1);}
  rules.push_back(q);
 }
 assert(rules.size()==30);
}
bool valid(uint32_t key){for(auto r:rules){int a=__builtin_parity(key&r.m[0])^r.flip[0],b=__builtin_parity(key&r.m[1])^r.flip[1],c=__builtin_parity(key&r.m[2])^r.flip[2];if(a==b&&b==c)return false;}return true;}
int main(int argc,char**argv){
 int cap1=argc>1?atoi(argv[1]):1,cap5=argc>2?atoi(argv[2]):3;assert(cap1<=2&&cap5<=3);init();
 vector<int>index(1<<20,-1);vector<uint32_t>keys;for(uint32_t k=0;k<(1u<<20);k++)if(valid(k)){index[k]=keys.size();keys.push_back(k);}
 int N=keys.size();assert(N==23808);vector<vector<Edge>>edges(N);vector<vector<uint32_t>>groups(1<<10);int ecount=0;
 int p5idx[10],cnt=0;for(int i=0;i<20;i++)if(P[i][2]<=5)p5idx[cnt++]=i;assert(cnt==10);
 auto group=[&](uint32_t k){int z=0;for(int i=0;i<10;i++)z|=((k>>p5idx[i])&1)<<i;return z;};
 for(int u=0;u<N;u++){
  uint32_t k=keys[u];groups[group(k)].push_back(k);
  for(int j=1;j<=5;j++)for(int b=0;b<20;b++){
   A3 I=P[b];bool has=false,blocked=false;for(int i:I){has|=i==j;blocked|=i==j+1;}if(!has||blocked)continue;
   for(int&i:I)if(i==j)i++;sort(I.begin(),I.end());int src=pid[I[0]][I[1]][I[2]];
   if(((k>>b)^(k>>src))&1){uint32_t dest=k^(1u<<b);if(index[dest]>=0){edges[u].push_back({dest,uint8_t(j),uint8_t(b==0)});ecount++;}}
  }
 }
 cout<<"necessary six-column states="<<N<<" event edges="<<ecount<<endl;

 // Independent direct last-label automaton; U2/U3/U4 event runs unbounded.
 int A=cap1+1,B=cap5+1;size_t size=size_t(N)*2*A*B*3;
 vector<int8_t>best(size,-1),ingroup(size,0);vector<int>broadcast(size_t(1<<10)*A*B,-1);deque<int>queue;
 auto encode=[&](int u,int a,int b,int c,int mode){return (((c*A+a)*B+b)*3+mode)*N+u;};
 uint32_t start=0;for(int i=0;i<20;i++)if(((P[i][0]+P[i][1]+P[i][2])&1)==0)start|=1u<<i;
 int root=encode(index[start],0,0,0,0);best[root]=0;ingroup[root]=1;queue.push_back(root);int maximum=0;long long popped=0,bcasts=0;
 auto push=[&](int id,int value){assert(value<=20);if(best[id]>=value)return;best[id]=value;maximum=max(maximum,value);if(!ingroup[id]){ingroup[id]=1;queue.push_back(id);}};
 while(!queue.empty()){
  int x=queue.front();queue.pop_front();ingroup[x]=0;int u=x%N,y=x/N,mode=y%3;y/=3;int b=y%B;y/=B;int a=y%A,c=y/A,val=best[x];
  for(auto e:edges[u]){
   int aa=a+(e.j==1&&mode!=1),bb=b+(e.j==5&&mode!=2),m=e.j==1?1:e.j==5?2:0;
   if(aa<=cap1&&bb<=cap5)push(encode(index[e.v],aa,bb,c,m),val+e.weight);
  }
  if(!c){int g=group(keys[u]),id=(a*B+b)*(1<<10)+g;
   if(broadcast[id]<val){broadcast[id]=val;bcasts++;for(uint32_t v:groups[g])push(encode(index[v],a,b,1,0),val);}
  }
  popped++;
 }
 if(cap1==1&&cap5==3)assert(maximum==11);
 cout<<"PASS independent unbounded interior: U1 runs <= "<<cap1<<", U5 runs <= "<<cap5<<", U6 jumps <= 1; max="<<maximum<<" popped="<<popped<<" broadcasts="<<bcasts<<endl;
}
