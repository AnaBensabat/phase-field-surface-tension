#include<bits/stdc++.h>
#include <format>
#include <iostream>
#include <string>
#include <string_view>
using namespace std;

vector<vector<double>> dummy(){
 
  return {{1,2},{3,4}}; 
}

int main(){

  auto x = dummy();
  for(auto w:x){
    for(auto z:w) cout<<z;
  }
  
  vector<vector<int>> mat {{1,2,3},{4,5,6}};
  
  ofstream nucleus_info;
  nucleus_info.open ("example.txt");
  for(auto x:mat){
    for(auto y:x) nucleus_info << y <<'\t';
    nucleus_info << '\n';
  }
  nucleus_info.close();

  double eta = 0.1;
  string s = to_string(eta);
  cout<<'\n'<<s<<'\n';
}
