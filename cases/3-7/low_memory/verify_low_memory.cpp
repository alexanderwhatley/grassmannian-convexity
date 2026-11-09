// Independent standalone verifier of Q <= 8+m for Grassmannian (3,7).
// Does not include or execute any original graph program. It reconstructs
// the 50-sign state set from affine-edge angular signs, independently
// derives all 114 feature derivatives by exact polynomial algebra, builds
// the augmented sign constraints directly from their identities, computes
// minimal row variation by positive four-circuits, and checks every edge.
//
// Run from this directory, with ../grassmann_37_continuation present:
//   c++ -std=c++17 -O3 -pthread standalone_q_verifier.cpp -o standalone_q_verifier
//   ./standalone_q_verifier
// Only the sorted rank/K binary arrays are required. directions_filtered.txt
// and keys_by_sorted_key_le64.bin, when present, are extra cross-checks.
#ifdef NDEBUG
#error "Assertions are required for the certificate"
#endif
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <fstream>
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
 string root="../grassmann_37_continuation/";ifstream f(root+"directions_filtered.txt");bool has_direction_reference=bool(f);
 vector<uint8_t>seen(1<<25);int a,b,nread=0;while(f>>a>>b){assert(a>=0&&a<(1<<15)&&b>=0&&b<(1<<10));int code=(a<<10)|b;assert(!seen[code]);seen[code]=1;nread++;}if(has_direction_reference)assert(nread==4560);
 auto rel=relations();array<uint64_t,16>masks{};for(int i=0;i<16;i++)masks[i]=sign_flip_mask(i);
 vector<uint64_t>keys;keys.reserve(16051200);int nf=0,nd=0;
 for(int x=0;x<(1<<15);x++)if(valid_f(x,rel)){
  nf++;for(int t=0;t<(1<<10);t++){
   array<int,10>pos{};if(!direction(x,t,pos))continue;int code=(x<<10)|t;
   if(has_direction_reference&&seen[code]!=1){cerr<<"actual angular direction omitted: "<<x<<' '<<t<<"\n";abort();}seen[code]=2;nd++;
   // Gap coordinates scaled by 8. Axis 7 is inserted CCW between axis 6
   // and its antipode: one choice in each same/opposite gap and 9 between.
   for(int gap=0;gap<20;gap++)for(int d=0;d<=10;d++){
    int z6=(8*gap+2)%160,z7=(d==0?8*gap+6:d==10?8*(gap+10)+1:8*(gap+d)+4)%160;
    uint64_t key=normalized_key(x,t,pos,z6,z7);for(auto m:masks)keys.push_back(key^m);
   }
  }
 }
 if(has_direction_reference)for(auto x:seen)assert(x!=1);assert(nf==1920&&nd==4560&&keys.size()==16051200);
 sort(keys.begin(),keys.end());assert(adjacent_find(keys.begin(),keys.end())==keys.end());
 ifstream check(root+"history_20261108/keys_by_sorted_key_le64.bin",ios::binary);bool compared=bool(check);
 if(check){for(size_t i=0;i<keys.size();i++){uint64_t x=0;for(int b=0;b<8;b++){int c=check.get();assert(c!=EOF);x|=uint64_t(c)<<(8*b);}if(x!=keys[i]){cerr<<"key mismatch at sorted index "<<i<<"\n";abort();}}assert(check.get()==EOF);}
 if(const char *path=getenv("GC37_DUMP_KEYS")){ofstream out(path,ios::binary);assert(out);for(auto x:keys)for(int b=0;b<8;b++)out.put(char(x>>(8*b)));assert(out);}
 cout<<"PASS 16 exterior identities; "<<nf<<" admissible F patterns; "<<nd<<" independent angular types; "<<keys.size()<<" reconstructed sorted keys"<<(compared?" match reference byte-for-byte":" (optional reference not supplied)")<<endl;
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
struct Hash {
 size_t cap;vector<uint64_t>key;vector<uint32_t>index;
 explicit Hash(size_t n):cap(n),key(n),index(n){}
 size_t hash(uint64_t x)const{x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return (x^(x>>31))&(cap-1);}
 void insert(uint64_t x,uint32_t v){assert(x);size_t i=hash(x);while(key[i]){assert(key[i]!=x);i=(i+1)&(cap-1);}key[i]=x;index[i]=v;}
 uint32_t find(uint64_t x)const{size_t i=hash(x);while(key[i]){if(key[i]==x)return index[i];i=(i+1)&(cap-1);}return UINT32_MAX;}
};
int main(){
 labels();check_affine_normalization();auto tests=symbolic_events();auto keys=generate_keys();const uint32_t n=keys.size();
 string dir="../grassmann_37_continuation/history_20261108/";vector<uint8_t>rank(n),old(n/2),minvar(n,255),valid(n);
 auto read=[&](string name,vector<uint8_t>&out){ifstream f(dir+name,ios::binary);assert(f);f.read((char*)out.data(),out.size());assert(f.gcount()==int64_t(out.size()));char ch;assert(!f.read(&ch,1));};
 read("rank_by_sorted_key.bin",rank);read("K_by_sorted_key.bin",old);
 auto inK=[&](uint32_t u,int g){int nib=(old[u/2]>>(4*(u&1)))&15;return (nib>>g)&1;};
 auto P=[&](uint32_t u,int g){return int(rank[u])-2*inK(u,g);};
 auto Qrank=[&](uint32_t u,int g){int p=P(u,g);return p-2*int(p>8+minvar[u]);};
 Hash H(1ull<<25);for(uint32_t u=0;u<n;u++)H.insert(keys[u],u);
 Variation variation;vector<thread>threads;uint64_t mcount[8][5]{};
 for(int th=0;th<8;th++)threads.emplace_back([&,th](){
  size_t lo=size_t(n)*th/8,hi=size_t(n)*(th+1)/8;
  for(size_t u=lo;u<hi;u++){
   int m=variation.m(keys[u]),karp=variation.karp(keys[u]);if(m!=karp||m>4){cerr<<"variation discrepancy index="<<u<<" circuit="<<m<<" karp="<<karp<<"\n";abort();}
   mcount[th][m]++;minvar[u]=m;uint8_t gmask=0;for(int g=0;g<4;g++)if(gvalid(keys[u],g))gmask|=1u<<g;valid[u]=gmask;
  }
 });
 for(auto&t:threads)t.join();threads.clear();
 cout<<"PASS four-circuit and Karp m agree on "<<n<<" bases; counts";
 for(int m=0;m<=4;m++){uint64_t total=0;for(int th=0;th<8;th++)total+=mcount[th][m];cout<<' '<<m<<':'<<total;}cout<<endl;
 uint64_t states[8]{},base_edges[8]{},g_edges[8]{},violations[8][5]{};int maxq[8][5]{};
 for(int th=0;th<8;th++)threads.emplace_back([&,th](){
  size_t lo=size_t(n)*th/8,hi=size_t(n)*(th+1)/8;
  for(size_t u=lo;u<hi;u++){
   uint64_t base=keys[u];int m=minvar[u];
   for(int g=0;g<4;g++)if(valid[u]&(1u<<g)){
    int p=P(u,g),q=Qrank(u,g);if(p<0||p>12||q<0||q>8+m){cerr<<"state rank failure index="<<u<<" g="<<g<<" m="<<m<<" P="<<p<<" Q="<<q<<"\n";abort();}
    states[th]++;violations[th][m]+=p>8+m;maxq[th][m]=max(maxq[th][m],q);
    uint64_t key=base|(uint64_t(g)<<50);
    for(auto e:tests){
     int s=signbit(key,e.target),d=e.e.sg*signbit(key,e.e.a);if(e.e.b>=0)d*=signbit(key,e.e.b);if(s==d)continue;
     if(e.target<50){
      uint32_t v=H.find(base^(1ull<<e.target));if(v==UINT32_MAX||!(valid[v]&(1u<<g)))continue;
      int w=(e.target==0),z=Qrank(v,g);if(q<z+w){cerr<<"base edge failure "<<u<<"->"<<v<<" feature="<<e.target<<" move="<<e.move<<" Q="<<q<<','<<z<<"\n";abort();}base_edges[th]++;
     }else{
      int h=g^(1<<(e.target-50));if(!(valid[u]&(1u<<h)))continue;int z=Qrank(u,h);if(q<z){cerr<<"G edge failure "<<u<<" feature="<<e.target<<"\n";abort();}g_edges[th]++;
     }
    }
   }
  }
 });
 for(auto&t:threads)t.join();
 uint64_t ns=0,nb=0,ng=0,nd=0;for(int th=0;th<8;th++){ns+=states[th];nb+=base_edges[th];ng+=g_edges[th];for(int m=0;m<=4;m++)nd+=violations[th][m];}
 assert(ns==23246848&&nb==209042944&&ng==21586944&&nd==3641);
 cout<<"PASS Q in [0,8+m] on "<<ns<<" admissible augmented states; D="<<nd<<" counts";
 for(int m=0;m<=4;m++){uint64_t c=0;int best=0;for(int th=0;th<8;th++){c+=violations[th][m];best=max(best,maxq[th][m]);}cout<<' '<<m<<':'<<c<<"(maxQ="<<best<<')';}cout<<endl;
 cout<<"PASS Q rank inequalities on "<<nb<<" base and "<<ng<<" G move-labeled edges"<<endl;
}
