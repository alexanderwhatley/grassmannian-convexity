// Dual-table extension of verify_upper_bound.cpp; shares its graph construction.
// Reconstructs H and checks H <= 9 + m* for the signed-reversed Gale dual.
// Checked-arithmetic data-free verifier of H <= 8+m for Grassmannian (3,7).
// Adds explicit signed-overflow checks to every sparse-polynomial and
// integer-witness arithmetic operation in the original standalone verifier.
// Does not include or execute any original graph program. It reconstructs
// the 50-sign state set from affine-edge angular signs, independently
// derives all 114 feature derivatives by exact polynomial algebra, builds
// the augmented sign constraints directly from their identities, computes
// minimal row variation by positive four-circuits, and checks every edge.
//
// Compile the source and run the executable from any directory:
//   c++ -std=c++17 -O3 -pthread standalone_data_free_checked.cpp -o standalone_data_free_checked
//   ./standalone_data_free_checked
// No input data files are read. All sign states, transitions, and ranks
// are reconstructed and directly verified from scratch.
#ifdef NDEBUG
#error "The verifier requires assertions; compile without -DNDEBUG"
#endif
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>
using namespace std;
using P3=array<int,3>;using F4=array<int,4>;using Edge=pair<int,int>;
vector<P3> plabs;vector<F4> flabs;vector<Edge> edges;int pid[8][8][8],eid[8][8];map<F4,int> fid;map<P3,int> tid;
void labels(){
 for(int a=0;a<8;a++)for(int b=0;b<8;b++){eid[a][b]=-1;for(int c=0;c<8;c++)pid[a][b][c]=-1;}
 for(int a=1;a<=7;a++)for(int b=a+1;b<=7;b++)for(int c=b+1;c<=7;c++){pid[a][b][c]=plabs.size();plabs.push_back({a,b,c});}
 for(int a=1;a<=5;a++)for(int b=a+1;b<=5;b++){eid[a][b]=edges.size();edges.push_back({a,b});}
 for(int a=1;a<=5;a++)for(int b=a+1;b<=5;b++)for(int c=b+1;c<=5;c++){
  tid[{a,b,c}]=tid.size();for(int d=c+1;d<=5;d++){flabs.push_back({a,b,c,d});flabs.push_back({a,c,b,d});flabs.push_back({a,d,b,c});}}
 for(size_t i=0;i<flabs.size();i++)fid[flabs[i]]=i;
 assert(plabs.size()==35&&flabs.size()==15&&edges.size()==10&&tid.size()==10);
}

// Formal affine differences Pj-Pi with P5=0, expressed in e1,...,e4.
array<int,6> bracket(F4 f){
 auto edge=[](int i,int j){array<int,4>x{};for(int a=1;a<=4;a++)x[a-1]=int(j==a)-int(i==a);return x;};
 auto u=edge(f[0],f[1]),v=edge(f[2],f[3]);array<int,6>out{};int k=0;
 for(int a=0;a<4;a++)for(int b=a+1;b<4;b++)out[k++]=u[a]*v[b]-u[b]*v[a];return out;
}
using Rel=vector<pair<int,int>>;
vector<Rel> relations(){
 vector<Rel>r={{{0,1},{3,-1},{6,1}},{{0,1},{11,-1},{14,1}},{{1,1},{4,-1},{9,1}},{{1,1},{8,-1},{13,-1}},{{2,1},{5,-1},{12,-1}},{{2,1},{7,-1},{10,1}},{{3,1},{10,-1},{13,1}},{{4,1},{7,-1},{14,-1}},{{5,1},{8,-1},{11,1}},{{6,1},{9,-1},{12,1}},
 {{0,1},{4,1},{8,-1},{10,-1},{12,-1}},{{0,1},{5,1},{7,-1},{9,1},{13,1}},{{1,1},{3,1},{7,-1},{11,-1},{12,1}},{{1,1},{5,-1},{6,1},{10,-1},{14,-1}},{{2,1},{3,1},{8,-1},{9,-1},{14,1}},{{2,1},{4,-1},{6,1},{11,1},{13,1}}};
 for(auto rel:r)for(int k=0;k<6;k++){int z=0;for(auto [i,s]:rel)z+=s*bracket(flabs[i])[k];assert(z==0);}
 return r;
}
bool valid_f(int f,const vector<Rel>&r){
 for(auto rel:r){int first=0;bool equal=true;for(auto [i,c]:rel){int x=c*((f>>i&1)?1:-1);if(!first)first=x;else equal&=x==first;}if(equal)return false;}return true;
}
// Reconstruct positive rays of all ten edge directions from f/t determinant
// signs. The signs determine exactly one oriented circle up to rotation.
bool direction(int f,int t,array<int,10>&position){
 int s[10][10]={};
 for(int a=0;a<10;a++)for(int b=a+1;b<10;b++){
  auto [i,j]=edges[a];auto[k,l]=edges[b];set<int>u={i,j,k,l};int bit;
  if(u.size()==3){auto it=u.begin();int A=*it++,B=*it++,C=*it++;bit=(t>>tid.at({A,B,C}))&1;}
  else {assert(u.size()==4);bit=(f>>fid.at({i,j,k,l}))&1;}
  s[a][b]=bit?1:-1;s[b][a]=-s[a][b];
 }
 int side[10]={},rank[10]={},order[9];fill(order,order+9,-1);
 for(int a=1;a<10;a++)side[a]=s[0][a];
 for(int a=1;a<10;a++)for(int b=a+1;b<10;b++){
  if(side[a]*side[b]*s[a][b]>0)rank[b]++;else rank[a]++;}
 for(int a=1;a<10;a++){if(rank[a]>=9||order[rank[a]]!=-1)return false;order[rank[a]]=a;}
 position[0]=0;for(int r=0;r<9;r++)position[order[r]]=r+1+(side[order[r]]<0?10:0);
 for(int a=0;a<10;a++)for(int b=a+1;b<10;b++){
  int delta=(position[b]-position[a]+20)%20;if(delta==0||delta==10||(delta<10?1:-1)!=s[a][b])return false;}
 return true;
}
uint64_t sign_flip_mask(int orientation){
 uint64_t out=0;
 for(size_t i=0;i<plabs.size();i++){
  int parity=0;for(int a:plabs[i])if(a<=4)parity^=(orientation>>(a-1))&1;if(parity)out|=1ull<<i;}
 for(size_t i=0;i<flabs.size();i++){
  int parity=0;for(int a:flabs[i])if(a<=4)parity^=(orientation>>(a-1))&1;if(parity)out|=1ull<<(35+i);}
 return out;
}
uint64_t normalized_key(int f,int t,const array<int,10>&pos,int axis6,int axis7){
 auto sdir=[&](int edge,int axis){int delta=(axis-8*pos[edge]+160)%160;assert(delta!=0&&delta!=80);return delta<80;};
 uint64_t key=0;
 for(size_t k=0;k<plabs.size();k++){
  auto p=plabs[k];bool positive;
  if(p[2]<=5)positive=(t>>tid.at(p))&1;
  else if(p[1]<=5)positive=sdir(eid[p[0]][p[1]],p[2]==6?axis6:axis7);
  else positive=true;
  if(positive)key|=1ull<<k;
 }
 for(int k=0;k<15;k++)if(f>>k&1)key|=1ull<<(35+k);
 return key;
}
vector<uint64_t> generate_keys(){
 auto rel=relations();array<uint64_t,16>masks{};for(int i=0;i<16;i++)masks[i]=sign_flip_mask(i);
 vector<uint64_t>keys;keys.reserve(16051200);int nf=0,nd=0;
 for(int x=0;x<(1<<15);x++)if(valid_f(x,rel)){
  nf++;for(int t=0;t<(1<<10);t++){
   array<int,10>pos{};if(!direction(x,t,pos))continue;nd++;
   for(int gap=0;gap<20;gap++)for(int d=0;d<=10;d++){
    int z6=(8*gap+2)%160,z7=(d==0?8*gap+6:d==10?8*(gap+10)+1:8*(gap+d)+4)%160;
    uint64_t key=normalized_key(x,t,pos,z6,z7);for(auto m:masks)keys.push_back(key^m);
   }
  }
 }
 assert(nf==1920&&nd==4560&&keys.size()==16051200);
 sort(keys.begin(),keys.end());assert(adjacent_find(keys.begin(),keys.end())==keys.end());
 cout<<"PASS 16 exterior identities; "<<nf<<" admissible F patterns; "<<nd<<" angular types; "<<keys.size()<<" reconstructed base states"<<endl;
 return keys;
}

// Exact sparse integer polynomials in the 21 entries of the 3-by-7 matrix.
// A monomial is a sorted list of variable indices padded with sentinel 255.
using Mono=array<uint8_t,9>;using Poly=map<Mono,int64_t>;
// These helpers perform the arithmetic in a wider type and assert its range;
// operations already in __int128 use compiler overflow-detection intrinsics.
int64_t narrow64(__int128 x){
 assert(x>=numeric_limits<int64_t>::min()&&x<=numeric_limits<int64_t>::max());return int64_t(x);
}
int64_t add64(int64_t a,int64_t b){return narrow64(__int128(a)+__int128(b));}
int64_t sub64(int64_t a,int64_t b){return narrow64(__int128(a)-__int128(b));}
int64_t mul64(int64_t a,int64_t b){return narrow64(__int128(a)*__int128(b));}
int64_t neg64(int64_t a){return narrow64(-__int128(a));}
__int128 add128(__int128 a,__int128 b){__int128 out;bool bad=__builtin_add_overflow(a,b,&out);assert(!bad);return out;}
__int128 sub128(__int128 a,__int128 b){__int128 out;bool bad=__builtin_sub_overflow(a,b,&out);assert(!bad);return out;}
__int128 mul128(__int128 a,__int128 b){__int128 out;bool bad=__builtin_mul_overflow(a,b,&out);assert(!bad);return out;}
Poly variable(int col,int row){assert(col>=1&&col<=7&&row>=0&&row<3);Mono a; a.fill(255);a[0]=3*(col-1)+row;return {{a,1}};}
void insert_term(Poly&out,Mono m,int64_t c){if(!c)return;auto it=out.find(m);if(it==out.end())out.emplace(m,c);else{it->second=add64(it->second,c);if(it->second==0)out.erase(it);}}
Poly add(Poly a,const Poly&b,int scale=1){for(auto &kv:b)insert_term(a,kv.first,mul64(scale,kv.second));return a;}
Poly neg(Poly a){for(auto&kv:a)kv.second=neg64(kv.second);return a;}
Poly mul(const Poly&a,const Poly&b){
 Poly out;for(auto&x:a)for(auto&y:b){Mono m;m.fill(255);int n=0;for(auto v:x.first)if(v!=255)m[n++]=v;for(auto v:y.first)if(v!=255){assert(n<9);m[n++]=v;}sort(m.begin(),m.begin()+n);insert_term(out,m,mul64(x.second,y.second));}return out;
}
Poly determinant(int a,int b,int c){
 array<int,3>cols={a,b,c},p={0,1,2};Poly out;
 do{int inversions=(p[0]>p[1])+(p[0]>p[2])+(p[1]>p[2]);
  Poly term=mul(mul(variable(cols[0],p[0]),variable(cols[1],p[1])),variable(cols[2],p[2]));out=add(out,term,inversions&1?-1:1);
 }while(next_permutation(p.begin(),p.end()));return out;
}
Poly feature(Edge ab,Edge cd,Edge ef){
 auto[a,b]=ab;auto[c,d]=cd;auto[e,f]=ef;
 return add(mul(determinant(a,b,e),determinant(c,d,f)),mul(determinant(a,b,f),determinant(c,d,e)),-1);
}
Poly derivative(const Poly&p,int j){
 Poly out;for(auto&term:p){int count=0;for(int i=0;i<9;i++)if(term.first[i]!=255&&term.first[i]/3==j-1){
  count++;Mono q=term.first;assert(q[i]<=252);q[i]=uint8_t(q[i]+3);sort(q.begin(),q.end());insert_term(out,q,term.second);
 }assert(count<=1);}
 return out;
}
struct Expr{int a,b,sg;};struct Test{int move,target;Expr e;};
int64_t evaluate(const Poly&p,const array<array<int64_t,3>,7>&v){
 __int128 total=0;for(auto&term:p){__int128 t=term.second;for(auto i:term.first)if(i!=255)t=mul128(t,v[i/3][i%3]);total=add128(total,t);}
 assert(total>=numeric_limits<int64_t>::min()&&total<=numeric_limits<int64_t>::max());return int64_t(total);
}
void generic_separation(const vector<Poly>&v,const vector<Test>&events){
 // A single exact witness suffices to prove each listed resultant is a
 // nonzero polynomial; this is a proof of generic event separability, not
 // an empirical assertion about paths.
 array<array<int64_t,3>,7>M={array<int64_t,3>{28,9,21},{6,17,26},{17,17,-29},{29,4,-4},{-26,-27,1},{17,-29,-12},{21,13,-29}};
 array<int64_t,52>value{};for(int i=0;i<52;i++){value[i]=evaluate(v[i],M);assert(value[i]);}
 bool active[7][52]={};for(auto e:events){assert(!active[e.move][e.target]);active[e.move][e.target]=true;}
 int checked=0;
 for(int j=1;j<=6;j++){
  auto M1=M;for(int k=0;k<3;k++)M1[j-1][k]=add64(M1[j-1][k],M[j][k]);array<int64_t,52>d{};
  for(int i=0;i<52;i++){d[i]=sub64(evaluate(v[i],M1),value[i]);assert((d[i]!=0)==active[j][i]);}
  for(int a=0;a<52;a++)if(active[j][a])for(int b=a+1;b<52;b++)if(active[j][b]){
   assert(sub128(mul128(value[a],d[b]),mul128(value[b],d[a]))!=0);checked++;
  }
 }
 assert(checked==1070);
 cout<<"PASS exact nonzero root-separation witnesses for all "<<checked<<" same-generator feature pairs"<<endl;
}
vector<Test> symbolic_events(){
 vector<Poly>v;for(auto p:plabs)v.push_back(determinant(p[0],p[1],p[2]));
 for(auto f:flabs)v.push_back(feature({f[0],f[1]},{f[2],f[3]},{6,7}));
 v.push_back(feature({1,2},{3,4},{5,6}));v.push_back(feature({1,2},{3,4},{5,7}));assert(v.size()==52);
 for(auto&p:v)assert(!p.empty());
 // Independently establish the two 2-vector dependencies and all cyclic
 // Pluecker identities used by the augmented-sign filter.
 Poly A5=determinant(1,2,5),A6=determinant(1,2,6),A7=determinant(1,2,7);
 Poly B5=determinant(3,4,5),B6=determinant(3,4,6),B7=determinant(3,4,7);
 for(auto a:{array<Poly,3>{A5,A6,A7},array<Poly,3>{B5,B6,B7}}){
  Poly sum=add(add(mul(a[0],v[35]),mul(a[1],v[51]),-1),mul(a[2],v[50]));assert(sum.empty());
 }
 for(int g=0;g<2;g++){
  array<Edge,3>p={Edge{1,2},Edge{3,4},Edge{5,g?7:6}};
  for(int i=0;i<3;i++)assert(feature(p[i],p[(i+1)%3],p[(i+2)%3])==v[50+g]);
 }
 map<Poly,Expr> express;
 auto register_poly=[&](const Poly&p,int a,int b){
  if(!express.count(p))express[p]={a,b,1};Poly q=neg(p);if(!express.count(q))express[q]={a,b,-1};
 };
 for(int i=0;i<52;i++)register_poly(v[i],i,-1);
 for(int i=0;i<35;i++)for(int j=i;j<35;j++)register_poly(mul(v[i],v[j]),i,j);
 vector<Test>out;
 for(int j=1;j<=6;j++)for(int i=0;i<52;i++){
  Poly d=derivative(v[i],j);if(d.empty())continue;auto it=express.find(d);
  if(it==express.end()){cerr<<"unrecognized exact derivative move="<<j<<" feature="<<i<<"\n";abort();}
  out.push_back({j,i,it->second});
 }
 assert(out.size()==114);
 int newevents=0;for(auto e:out)newevents+=e.target>=50;assert(newevents==6);
 cout<<"PASS exact symbolic reconstruction of all "<<out.size()<<" nonzero derivatives and 198 zero derivatives; G identities verified"<<endl;
 generic_separation(v,out);
 return out;
}
Poly one(){Mono a;a.fill(255);return {{a,1}};}
Poly detpoly(const array<Poly,3>&a,const array<Poly,3>&b,const array<Poly,3>&c){
 array<int,3>p={0,1,2};Poly out;do{
  int inversions=(p[0]>p[1])+(p[0]>p[2])+(p[1]>p[2]);out=add(out,mul(mul(a[p[0]],b[p[1]]),c[p[2]]),inversions&1?-1:1);
 }while(next_permutation(p.begin(),p.end()));return out;
}
void check_affine_normalization(){
 array<array<Poly,3>,8>v;
 for(int i=1;i<=5;i++){
  auto x=variable(i,0),y=variable(i,1),z=variable(i,2);v[i]={x,mul(x,y),mul(x,z)};
 }
 v[6]={Poly{},one(),Poly{}};v[7]={Poly{},Poly{},one()};
 for(auto p:plabs){
  int a=p[0],b=p[1],c=p[2];Poly lhs=detpoly(v[a],v[b],v[c]),rhs;
  if(c<=5){
   Poly dyab=add(variable(b,1),variable(a,1),-1),dzab=add(variable(b,2),variable(a,2),-1);
   Poly dyac=add(variable(c,1),variable(a,1),-1),dzac=add(variable(c,2),variable(a,2),-1);
   rhs=add(mul(dyab,dzac),mul(dzab,dyac),-1);
   rhs=mul(mul(mul(rhs,variable(a,0)),variable(b,0)),variable(c,0));
  }else if(b<=5){
   rhs=c==6?add(variable(a,2),variable(b,2),-1):add(variable(b,1),variable(a,1),-1);
   rhs=mul(mul(rhs,variable(a,0)),variable(b,0));
  }else rhs=variable(a,0);
  assert(lhs==rhs);
 }
 for(auto f:flabs){int a=f[0],b=f[1],c=f[2],d=f[3];
  Poly F=add(mul(detpoly(v[a],v[b],v[6]),detpoly(v[c],v[d],v[7])),mul(detpoly(v[a],v[b],v[7]),detpoly(v[c],v[d],v[6])),-1);
  Poly dyab=add(variable(b,1),variable(a,1),-1),dzab=add(variable(b,2),variable(a,2),-1);
  Poly dycd=add(variable(d,1),variable(c,1),-1),dzcd=add(variable(d,2),variable(c,2),-1);
  Poly rhs=add(mul(dyab,dzcd),mul(dzab,dycd),-1);
  for(int i:f)rhs=mul(rhs,variable(i,0));assert(F==rhs);
 }
 cout<<"PASS all 35 minor and 15 feature exact affine-normalization identities"<<endl;
}
int signbit(uint64_t key,int bit){return (key>>bit&1)?1:-1;}
int oriented(uint64_t k,int a,int b,int c){int inv=(a>b)+(a>c)+(b>c);array<int,3>z={a,b,c};sort(z.begin(),z.end());return (inv&1?-1:1)*signbit(k,pid[z[0]][z[1]][z[2]]);}
bool gvalid(uint64_t k,int g){
 int q[9]={signbit(k,pid[1][2][5]),signbit(k,pid[3][4][5]),signbit(k,pid[1][2][6]),signbit(k,pid[3][4][6]),signbit(k,pid[1][2][7]),signbit(k,pid[3][4][7]),(g&1)?1:-1,(g&2)?1:-1,signbit(k,35)};
 auto forced=[](int a,int b,int want){return a==b||a==want;};
 if(!forced(q[0]*q[3],q[2]*q[1],q[6])||!forced(q[0]*q[5],q[4]*q[1],q[7])||!forced(q[2]*q[5],q[4]*q[3],q[8]))return false;
 if(abs(q[0]*q[8]-q[2]*q[7]+q[4]*q[6])==3||abs(q[1]*q[8]-q[3]*q[7]+q[5]*q[6])==3)return false;
 for(int t=0;t<2;t++){array<Edge,3>p={Edge{1,2},Edge{3,4},Edge{5,t?7:6}};int want=(g>>t&1)?1:-1;
  for(int r=0;r<3;r++){auto[A,B]=p[r];auto[C,D]=p[(r+1)%3];auto[E,F]=p[(r+2)%3];int x=oriented(k,A,B,E)*oriented(k,C,D,F),y=oriented(k,A,B,F)*oriented(k,C,D,E);if(x!=y&&x!=want)return false;}}
 return true;
}
struct Circuit {int pa,pb,pc,pd;uint64_t forbidden[16];};
struct Variation {
 uint64_t varies[7]={};vector<Circuit>circuits;vector<array<pair<int,int>,5>>sequences;
 Variation(){
  for(int mask=0;mask<64;mask++){
   int s[7];for(int i=0;i<6;i++)s[i]=(mask>>i&1)?-1:1;s[6]=1;int n=0;for(int i=1;i<7;i++)n+=s[i]!=s[i-1];varies[n]|=1ull<<mask;
  }
  for(int a=1;a<=7;a++)for(int b=a+1;b<=7;b++)for(int c=b+1;c<=7;c++)for(int d=c+1;d<=7;d++){
   Circuit x{pid[b][c][d],pid[a][c][d],pid[a][b][d],pid[a][b][c],{}};
   for(int nib=0;nib<16;nib++)for(int mask=0;mask<64;mask++){
    int s[8];for(int i=1;i<=6;i++)s[i]=(mask>>(i-1)&1)?-1:1;s[7]=1;
    int z[4]={s[a]*((nib&1)?1:-1),-s[b]*((nib&2)?1:-1),s[c]*((nib&4)?1:-1),-s[d]*((nib&8)?1:-1)};
    if(z[0]==z[1]&&z[1]==z[2]&&z[2]==z[3])x.forbidden[nib]|=1ull<<mask;
   }
   for(auto q:x.forbidden)assert(__builtin_popcountll(q)==8);circuits.push_back(x);
  }
  assert(circuits.size()==35);
  for(int a=1;a<=7;a++)for(int b=a+1;b<=7;b++){
   array<pair<int,int>,5>v;int k=0;for(int i=1;i<=7;i++)if(i!=a&&i!=b){array<int,3>q={a,b,i};sort(q.begin(),q.end());v[k++]={pid[q[0]][q[1]][q[2]],i&1};}assert(k==5);sequences.push_back(v);
  }
 }
 int m(uint64_t key)const{
  uint64_t possible=~0ull;for(auto&c:circuits){int nib=((key>>c.pa)&1)|(((key>>c.pb)&1)<<1)|(((key>>c.pc)&1)<<2)|(((key>>c.pd)&1)<<3);possible &=~c.forbidden[nib];}
  for(int i=0;i<=6;i++)if(possible&varies[i])return i;return 100;
 }
 int karp(uint64_t key)const{
  int best=0;for(auto&seq:sequences){int bits=0;for(int i=0;i<5;i++)bits|=(((key>>seq[i].first)&1)^seq[i].second)<<i;best=max(best,__builtin_popcount(unsigned((bits^(bits>>1))&15)));}return 4-best;
 }
};

// Dual is the reversed alternating reorientation of the orthogonal complement.
// Its ordered 4-minor signs agree, up to one global sign, with the original
// 3-minors on the complements of reversed indices.
struct DualVariation {
 vector<array<pair<int,int>,4>> seq;
 DualVariation(){
  for(int a=1;a<=7;a++)for(int b=a+1;b<=7;b++)for(int c=b+1;c<=7;c++){
   array<pair<int,int>,4> q;int z=0;
   for(int i=1;i<=7;i++)if(i!=a&&i!=b&&i!=c){
    bool used[8]={}; for(int x:{a,b,c,i}) used[8-x]=true;
    array<int,3> p{};int j=0;for(int x=1;x<=7;x++)if(!used[x])p[j++]=x;
    assert(j==3);q[z++]={pid[p[0]][p[1]][p[2]],i&1};
   }
   seq.push_back(q);
  }
 }
 int m(uint64_t key)const{
  int best=0;for(auto &q:seq){int bits=0;for(int j=0;j<4;j++)bits|=(((key>>q[j].first)&1)^q[j].second)<<j;best=max(best,__builtin_popcount(unsigned((bits^(bits>>1))&7)));}
  return 3-best;
 }
};

struct Hash {
 size_t cap;vector<uint64_t>key;vector<uint32_t>index;
 explicit Hash(size_t n):cap(n),key(n),index(n){}
 size_t hash(uint64_t x)const{x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return (x^(x>>31))&(cap-1);}
 void insert(uint64_t x,uint32_t v){assert(x);size_t i=hash(x);while(key[i]){assert(key[i]!=x);i=(i+1)&(cap-1);}key[i]=x;index[i]=v;}
 uint32_t find(uint64_t x)const{size_t i=hash(x);while(key[i]){if(key[i]==x)return index[i];i=(i+1)&(cap-1);}return UINT32_MAX;}
};
int main(){
 labels();check_affine_normalization();auto tests=symbolic_events();auto keys=generate_keys();const uint32_t n=keys.size();
 Hash lookup(1ull<<25);for(uint32_t u=0;u<n;u++)lookup.insert(keys[u],u);
 Variation variation;DualVariation dv;vector<uint8_t>minvar(n,255),minstar(n,255),valid(n);vector<thread>threads;uint64_t mcount[8][5]{};
 for(int th=0;th<8;th++)threads.emplace_back([&,th](){
  size_t lo=size_t(n)*th/8,hi=size_t(n)*(th+1)/8;
  for(size_t u=lo;u<hi;u++){
   int m=variation.m(keys[u]),karp=variation.karp(keys[u]);if(m!=karp||m>4){cerr<<"variation discrepancy index="<<u<<" circuit="<<m<<" karp="<<karp<<"\n";abort();}
   mcount[th][m]++;minvar[u]=m; minstar[u]=dv.m(keys[u]);uint8_t gmask=0;for(int g=0;g<4;g++)if(gvalid(keys[u],g))gmask|=1u<<g;valid[u]=gmask;
  }
 });
 for(auto&t:threads)t.join();threads.clear();
 cout<<"PASS four-circuit and Karp m agree on "<<n<<" bases; counts";
 for(int m=0;m<=4;m++){uint64_t total=0;for(int th=0;th<8;th++)total+=mcount[th][m];cout<<' '<<m<<':'<<total;}cout<<endl;

 // Dense indices for all admissible 52-sign states, with no input files.
 vector<uint32_t>dense(size_t(n)*4,UINT32_MAX),code;code.reserve(23246848);
 for(uint32_t u=0;u<n;u++)for(int g=0;g<4;g++)if(valid[u]&(1u<<g)){
  uint32_t c=4*u+g;dense[c]=code.size();code.push_back(c);
 }
 const uint32_t N=code.size();assert(N==23246848);
 // Build the directed graph, retaining one edge per distinct feature flip.
 // Multiple eligible generator labels for that flip do not alter the edge.
 vector<uint8_t>degrees(N,0);vector<vector<uint32_t>>chunks(8);atomic<int>done{0};
 for(int th=0;th<8;th++)threads.emplace_back([&,th](){
  size_t lo=size_t(N)*th/8,hi=size_t(N)*(th+1)/8;auto&out=chunks[th];out.reserve((hi-lo)*7);
  for(size_t id=lo;id<hi;id++){
   uint32_t c=code[id],u=c/4;int g=c%4;uint64_t base=keys[u],key=base|(uint64_t(g)<<50),candidate=0;
   for(auto e:tests){int s=signbit(key,e.target),d=e.e.sg*signbit(key,e.e.a);if(e.e.b>=0)d*=signbit(key,e.e.b);if(s!=d)candidate|=1ull<<e.target;}
   int deg=0;while(candidate){
    int bit=__builtin_ctzll(candidate);candidate&=candidate-1;uint32_t v=u;int h=g;
    if(bit<50){v=lookup.find(base^(1ull<<bit));if(v==UINT32_MAX||!(valid[v]&(1u<<h)))continue;}
    else {h=g^(1<<(bit-50));if(!(valid[v]&(1u<<h)))continue;}
    uint32_t dest=dense[4*v+h];assert(dest<N);out.push_back(dest|(uint32_t(bit==0)<<31));deg++;
   }
   assert(deg<256);degrees[id]=deg;
  }
  cerr<<"graph workers finished "<<++done<<"/8"<<endl;
 });
 for(auto&t:threads)t.join();threads.clear();
 vector<uint32_t>offset(N+1,0);for(uint32_t i=0;i<N;i++)offset[i+1]=offset[i]+degrees[i];uint32_t M=offset[N];
 uint64_t actual=0;for(auto&c:chunks)actual+=c.size();assert(actual==M&&actual<UINT32_MAX);
 vector<uint32_t>arcs;arcs.reserve(M);for(auto&c:chunks){arcs.insert(arcs.end(),c.begin(),c.end());vector<uint32_t>().swap(c);}assert(arcs.size()==M);
 cout<<"PASS constructed "<<N<<" augmented states and "<<M<<" distinct directed feature-flip edges"<<endl;

 // Reverse adjacency and exact iterative Kosaraju SCC decomposition.
 vector<uint32_t>reverse_offset(N+1,0);
 for(uint32_t e:arcs)reverse_offset[(e&0x7fffffffu)+1]++;
 for(uint32_t i=0;i<N;i++)reverse_offset[i+1]+=reverse_offset[i];assert(reverse_offset[N]==M);
 vector<uint32_t>cursor=reverse_offset,reverse_arcs(M);
 for(uint32_t u=0;u<N;u++)for(uint32_t e=offset[u];e<offset[u+1];e++){
  uint32_t v=arcs[e]&0x7fffffffu;reverse_arcs[cursor[v]++]=u;
 }
 vector<uint32_t>().swap(cursor);
 cerr<<"reverse graph constructed"<<endl;
 vector<uint8_t>visited(N,0),next(N,0);vector<uint32_t>order,stack;order.reserve(N);
 for(uint32_t root=0;root<N;root++)if(!visited[root]){
  visited[root]=1;stack.push_back(root);
  while(!stack.empty()){
   uint32_t u=stack.back();if(next[u]<degrees[u]){
    uint32_t v=arcs[offset[u]+next[u]++]&0x7fffffffu;if(!visited[v]){visited[v]=1;stack.push_back(v);}
   }else{stack.pop_back();order.push_back(u);}
  }
 }
 assert(order.size()==N);vector<uint8_t>().swap(visited);vector<uint8_t>().swap(next);
 cerr<<"forward SCC traversal complete"<<endl;
 vector<uint32_t>component(N,UINT32_MAX);uint32_t components=0,largest=0;
 for(uint64_t z=N;z>0;z--){uint32_t root=order[z-1];if(component[root]!=UINT32_MAX)continue;
  uint32_t size=0;component[root]=components;stack.push_back(root);
  while(!stack.empty()){
   uint32_t u=stack.back();stack.pop_back();size++;
   for(uint32_t e=reverse_offset[u];e<reverse_offset[u+1];e++){
    uint32_t v=reverse_arcs[e];if(component[v]==UINT32_MAX){component[v]=components;stack.push_back(v);}
   }
  }
  largest=max(largest,size);components++;
 }
 vector<uint32_t>().swap(order);vector<uint32_t>().swap(stack);vector<uint32_t>().swap(reverse_offset);vector<uint32_t>().swap(reverse_arcs);
 cerr<<"SCC components="<<components<<" largest="<<largest<<endl;
 // Components are topologically ordered; compute the maximal target-weighted
 // path rank in reverse order. A positive internal edge would be a positive
 // cycle and therefore invalid; reject rather than assume it does not occur.
 vector<uint32_t>component_offset(components+1,0);
 for(uint32_t c:component)component_offset[c+1]++;
 for(uint32_t c=0;c<components;c++)component_offset[c+1]+=component_offset[c];assert(component_offset[components]==N);
 vector<uint32_t>fill=component_offset,members(N);
 for(uint32_t u=0;u<N;u++)members[fill[component[u]]++]=u;
 vector<uint32_t>().swap(fill);vector<uint8_t>h(components,0);
 for(uint64_t zz=components;zz>0;zz--){uint32_t c=zz-1;unsigned best=0;
  for(uint32_t i=component_offset[c];i<component_offset[c+1];i++){
   uint32_t u=members[i];for(uint32_t e=offset[u];e<offset[u+1];e++){
    uint32_t x=arcs[e],v=x&0x7fffffffu,w=x>>31,d=component[v];
    if(d==c){if(w){cerr<<"positive target edge in SCC"<<endl;abort();}}
    else {if(d<c){cerr<<"invalid SCC topological ordering"<<endl;abort();}best=max(best,w+unsigned(h[d]));}
   }
  }
  assert(best<=255);h[c]=best;
 }
 int maximum=*max_element(h.begin(),h.end());assert(maximum==12);
 cout<<"PASS computed target rank H from graph: "<<components<<" SCCs; largest="<<largest<<"; maximum H="<<maximum<<endl;
 vector<uint32_t>().swap(component_offset);vector<uint32_t>().swap(members);
 // The following direct exhaustive check is the actual rank certificate; it
 // does not rely on the correctness of the SCC construction or max-path claim.
 uint64_t states[8]{},base_edges[8]{},g_edges[8]{};int maxh[8][5]{}; uint64_t cells[8][5][4]{}; int mh[8][5][4]{};
 for(int th=0;th<8;th++)threads.emplace_back([&,th](){
  size_t lo=size_t(N)*th/8,hi=size_t(N)*(th+1)/8;
  for(size_t id=lo;id<hi;id++){
   uint32_t c=code[id],u=c/4;int g=c%4,m=minvar[u],r=h[component[id]];
   if(r<0||r>8+m){cerr<<"H state bound failure index="<<id<<" m="<<m<<" H="<<r<<"\n";abort();}
   int dvar=minstar[u]; if(r>9+dvar){cerr<<"DUAL BOUND FAILURE "<<u<<" m="<<m<<" mstar="<<dvar<<" H="<<r<<"\n";abort();} cells[th][m][dvar]++; mh[th][m][dvar]=max(mh[th][m][dvar],r); states[th]++;maxh[th][m]=max(maxh[th][m],r);uint64_t base=keys[u],key=base|(uint64_t(g)<<50);
   for(auto e:tests){
    int s=signbit(key,e.target),d=e.e.sg*signbit(key,e.e.a);if(e.e.b>=0)d*=signbit(key,e.e.b);if(s==d)continue;
    uint32_t v=u;int t=g;if(e.target<50){v=lookup.find(base^(1ull<<e.target));if(v==UINT32_MAX||!(valid[v]&(1u<<t)))continue;}
    else {t=g^(1<<(e.target-50));if(!(valid[v]&(1u<<t)))continue;}
    uint32_t dest=dense[4*v+t];assert(dest<N);int w=(e.target==0),r2=h[component[dest]];
    if(r<r2+w){cerr<<"H inequality failure id="<<id<<" dest="<<dest<<" feature="<<e.target<<" move="<<e.move<<" H="<<r<<','<<r2<<"\n";abort();}
    if(e.target<50)base_edges[th]++;else g_edges[th]++;
   }
  }
 });
 for(auto&t:threads)t.join();
 uint64_t ns=0,nb=0,ng=0;for(int th=0;th<8;th++){ns+=states[th];nb+=base_edges[th];ng+=g_edges[th];}
 assert(ns==23246848&&nb==209042944&&ng==21586944);
 cout<<"PASS H in [0,8+m] on all "<<ns<<" admissible states; maxima by m";
 for(int m=0;m<=4;m++){int x=0;for(int th=0;th<8;th++)x=max(x,maxh[th][m]);cout<<' '<<m<<':'<<x;}cout<<endl;
 int occupied=0;for(int a=0;a<=4;a++)for(int b=0;b<=3;b++){uint64_t ncell=0;int rcell=0;for(int th=0;th<8;th++){ncell+=cells[th][a][b];rcell=max(rcell,mh[th][a][b]);}if(ncell){occupied++;assert(rcell==min(8+a,9+b));cout<<"DUAL CELL m="<<a<<" mstar="<<b<<" n="<<ncell<<" maxH="<<rcell<<"\n";}} assert(occupied==13);
 cout<<"PASS H <= 9 + signed-reversed dual variation on all states"<<endl;
 cout<<"PASS direct H rank inequalities on "<<nb<<" base and "<<ng<<" G move-labeled edges"<<endl;
}
