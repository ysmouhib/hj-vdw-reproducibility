#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
using namespace std;

static vector<int> primes_factor(int n){vector<int> f;for(int q=2;q*(long long)q<=n;q++)if(n%q==0){f.push_back(q);while(n%q==0)n/=q;}if(n>1)f.push_back(n);return f;}
static int modpow(long long a,long long e,int p){long long z=1;while(e){if(e&1)z=z*a%p;a=a*a%p;e>>=1;}return z;}
static int primitive_root(int p){auto f=primes_factor(p-1);for(int g=2;g<p;g++){bool ok=true;for(int q:f)if(modpow(g,(p-1)/q,p)==1){ok=false;break;}if(ok)return g;}return -1;}

struct Search {
 int p,r,t,g; vector<int> chi,S,bad,weight; vector<unsigned short> cnt; vector<char> inS;
 mt19937 rng; int curBad=0,bestBad=0; long long curCost=0,trials=0,moves=0,restarts=0;
 vector<int> bestS;
 Search(int p_,int r_,int t_,int seed):p(p_),r(r_),t(t_),g(primitive_root(p_)),chi(p_),bad(p_),weight(p_,1),cnt(p_*(r_+1)),inS(p_),rng(seed){
  chi[0]=r;int x=1;for(int i=0;i<p-1;i++){chi[x]=i%r;x=(int)((long long)x*g%p);}
 }
 bool bad_at(int y)const{
  int z=cnt[y*(r+1)+r],need=t-z;
  for(int c=0;c<r;c++)if(cnt[y*(r+1)+c]==need)return true;
  return false;
 }
 void rebuild(){
  fill(cnt.begin(),cnt.end(),0);curBad=0;curCost=0;
  for(int y=0;y<p;y++)for(int s:S)cnt[y*(r+1)+chi[(y+s)%p]]++;
  for(int y=0;y<p;y++){bad[y]=bad_at(y);curBad+=bad[y];if(bad[y])curCost+=weight[y];}
 }
 void random_set(){
  fill(inS.begin(),inS.end(),0);S={0,1};inS[0]=inS[1]=1;fill(weight.begin(),weight.end(),1);
  while((int)S.size()<t){int z=2+rng()%(p-2);if(!inS[z]){inS[z]=1;S.push_back(z);}}
  rebuild();
 }
 int trial(int slot,int z,bool undo=false){ // applies replacement but leaves S/inS unchanged
  int old=S[slot];for(int y=0;y<p;y++){
   int a=chi[(y+(undo?z:old))%p],b=chi[(y+(undo?old:z))%p];if(a==b)continue;
   curBad-=bad[y];if(bad[y])curCost-=weight[y];int base=y*(r+1);cnt[base+a]--;cnt[base+b]++;bad[y]=bad_at(y);curBad+=bad[y];if(bad[y])curCost+=weight[y];
  }trials++;return curCost;
 }
 void accept(int slot,int z){int old=S[slot];inS[old]=0;inS[z]=1;S[slot]=z;moves++;}
 void save(const string&path,double elapsed)const{
  ofstream f(path);f<<"{\"p\":"<<p<<",\"r\":"<<r<<",\"t\":"<<t<<",\"g\":"<<g<<",\"bad\":"<<bestBad<<",\"weights\":[";
  for(int i=0;i<t;i++){if(i)f<<',';f<<bestS[i];}f<<"],\"moves\":"<<moves<<",\"trials\":"<<trials<<",\"restarts\":"<<restarts<<",\"seconds\":"<<elapsed<<"}\n";
 }
 int run(double seconds,const string&path,const vector<int>&initial){
  auto start=chrono::steady_clock::now();auto elapsed=[&](){return chrono::duration<double>(chrono::steady_clock::now()-start).count();};
  if(!initial.empty()){S=initial;fill(inS.begin(),inS.end(),0);for(int x:S)inS[x]=1;rebuild();}else random_set();
  bestBad=curBad;bestS=S;save(path,0);if(bestBad==0){cout<<"FOUND initial\n";return 0;}
  const int CAND=32;long long stale=0;
  while(elapsed()<seconds){
   if(curBad<bestBad){int tracked=curBad;rebuild();if(curBad!=tracked){cerr<<"INCREMENTAL_MISMATCH tracked="<<tracked<<" rebuilt="<<curBad<<"\n";return 3;}bestBad=curBad;bestS=S;stale=0;save(path,elapsed());cerr<<"BEST p="<<p<<" bad="<<bestBad<<" moves="<<moves<<" trials="<<trials<<" sec="<<elapsed()<<"\n";if(bestBad==0){cout<<"FOUND\n";return 0;}}
   if(stale>2500){restarts++;if(restarts%3==0)random_set();else{S=bestS;fill(inS.begin(),inS.end(),0);for(int x:S)inS[x]=1;for(int q=0;q<max(2,t/2);q++){int slot=2+rng()%(t-2),z;do z=2+rng()%(p-2);while(inS[z]);trial(slot,z);accept(slot,z);}rebuild();}stale=0;continue;}
   int y=rng()%p;while(!bad[y])y=rng()%p;
   for(int q=0;q<7;q++){int z=rng()%p;if(bad[z]&&weight[z]>weight[y])y=z;}
   int slot=2+rng()%(t-2);int old=S[slot];
   int chosen=-1;long long chosenScore=1000000000000LL;
   for(int q=0;q<CAND;q++){
    int z;if(q==0)z=(old+1)%p;else if(q==1)z=(old+p-1)%p;else z=2+rng()%(p-2);
    if(z<2||inS[z])continue;
    int score=trial(slot,z);trial(slot,z,true); // undo
    if(score<chosenScore){chosenScore=score;chosen=z;}
   }
   stale++;
   if(chosen>=0&&chosenScore>=curCost&&(rng()%1000)<250){
    for(int z=0;z<p;z++)if(bad[z]){weight[z]++;curCost++;}
    continue;
   }
   if(chosen<0)continue;
   bool take=chosenScore<=curBad || (rng()%1000)<12;
   if(take){trial(slot,chosen);accept(slot,chosen);}else if((rng()%1000)<25){trial(slot,chosen);accept(slot,chosen);}
  }
  save(path,elapsed());cout<<"DONE best="<<bestBad<<" moves="<<moves<<" trials="<<trials<<" restarts="<<restarts<<"\n";return 1;
 }
};
int main(int argc,char**argv){
 if(argc<7){cerr<<"args p r t seconds seed output [comma-separated initial weights]\n";return 2;}
 int p=stoi(argv[1]),r=stoi(argv[2]),t=stoi(argv[3]);double sec=stod(argv[4]);int seed=stoi(argv[5]);string out=argv[6];vector<int>init;
 if(argc>7){string s=argv[7];size_t pos=0;while(pos<s.size()){auto q=s.find(',',pos);if(q==string::npos)q=s.size();init.push_back(stoi(s.substr(pos,q-pos)));pos=q+1;}}
 Search s(p,r,t,seed);return s.run(sec,out,init);
}
