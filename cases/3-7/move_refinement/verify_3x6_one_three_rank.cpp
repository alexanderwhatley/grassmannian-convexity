// Direct potential audit for the 3x6, <=1 U1, <=3 U5, <=1 U6 model.
// Independent reconstruction using oriented three-term GP signs. The rank
// byte is only a certificate; every reachable edge inequality is checked.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
struct E{int to,j,w;};
int main(){
 vector<array<int,3>>p;int id[7][7][7]{};for(int a=1;a<=6;a++)for(int b=a+1;b<=6;b++)for(int c=b+1;c<=6;c++){id[a][b][c]=p.size();p.push_back({a,b,c});}
 auto orient=[&](uint32_t s,int a,int b,int c){int parity=(a>b)+(a>c)+(b>c);array<int,3>x={a,b,c};sort(x.begin(),x.end());int z=(s>>id[x[0]][x[1]][x[2]]&1)?1:-1;return parity%2?-z:z;};
 vector<uint32_t>keys;vector<int>rev(1<<20,-1),group;vector<vector<int>>members(1024);vector<int>first;
 for(int a=1;a<=5;a++)for(int b=a+1;b<=5;b++)for(int c=b+1;c<=5;c++)first.push_back(id[a][b][c]);
 for(uint32_t s=0;s<(1u<<20);s++){
  bool good=true;
  for(int z=1;z<=6&&good;z++)for(int a=1;a<=6&&good;a++)for(int b=a+1;b<=6&&good;b++)for(int c=b+1;c<=6&&good;c++)for(int d=c+1;d<=6;d++){
   if(z==a||z==b||z==c||z==d)continue;
   int A=orient(s,z,a,b)*orient(s,z,c,d),B=-orient(s,z,a,c)*orient(s,z,b,d),C=orient(s,z,a,d)*orient(s,z,b,c);
   if(A==B&&B==C){good=false;break;}
  }
  if(!good)continue;int q=0;for(int i=0;i<10;i++)q|=((s>>first[i])&1)<<i;
  rev[s]=keys.size();group.push_back(q);members[q].push_back(keys.size());keys.push_back(s);
 }
 int N=keys.size();assert(N==23808);vector<vector<E>>arcs(N);uint64_t edge_count=0;
 for(int a=0;a<N;a++)for(int j=1;j<=5;j++)for(int t=0;t<20;t++){
  auto A=p[t];if(find(A.begin(),A.end(),j)==A.end()||find(A.begin(),A.end(),j+1)!=A.end())continue;
  for(int &z:A)if(z==j)z++;sort(A.begin(),A.end());int derivative=id[A[0]][A[1]][A[2]];
  if(((keys[a]>>t)&1)==((keys[a]>>derivative)&1))continue;
  int b=rev[keys[a]^(1u<<t)];if(b<0)continue;arcs[a].push_back({b,j,int(t==id[1][2][3])});edge_count++;
 }
 assert(edge_count==115200);const int A=2,B=4,G=2,M=3;size_t SZ=size_t(N)*A*B*G*M;
 vector<int8_t>score(SZ);ifstream f("three_by_six_1u1_3u5_rank.bin",ios::binary);f.read((char*)score.data(),score.size());assert(f.gcount()==int64_t(score.size()));char extra;assert(!f.read(&extra,1));
 auto code=[&](int s,int a,int b,int g,int mode){return (((((g*A+a)*B+b)*M+mode)*N)+s);};
 uint32_t alt=0;for(int i=0;i<20;i++)if((p[i][0]+p[i][1]+p[i][2])%2==0)alt|=1u<<i;int root=rev[alt];assert(root>=0&&score[code(root,0,0,0,0)]==0);
 uint64_t checked=0,reachable=0,jumpchecks=0;int maximum=-1;
 for(int g=0;g<2;g++)for(int a=0;a<A;a++)for(int b=0;b<B;b++)for(int mode=0;mode<M;mode++)for(int s=0;s<N;s++){
  int r=score[code(s,a,b,g,mode)];if(r<0)continue;reachable++;maximum=max(maximum,r);assert(r<=11);
  auto check=[&](int t,int aa,int bb,int gg,int mm,int w){assert(int(score[code(t,aa,bb,gg,mm)])>=r+w);checked++;};
  if(mode)check(s,a,b,g,0,0);
  for(auto e:arcs[s]){
   if(mode==1&&e.j!=1)continue;if(mode==2&&e.j!=5)continue;
   int aa=a,bb=b,mm=mode;if(mode==0&&e.j==1){if(a==1)continue;aa++;mm=1;}
   if(mode==0&&e.j==5){if(b==3)continue;bb++;mm=2;}
   check(e.to,aa,bb,g,mm,e.w);
  }
  if(mode==0&&g==0)for(int t:members[group[s]]){check(t,a,b,1,0,0);jumpchecks++;}
 }
 cout<<"PASS independent 3x6 potential: states="<<N<<" labeled_edges="<<edge_count<<" reachable_certificate_states="<<reachable<<" max="<<maximum<<" checked_inequalities="<<checked<<" U6_jumps="<<jumpchecks<<"\n";
}
