#include <vector>
#include <string.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP mean_shift(SEXP poslist,int step){
  NumericVector pos =as<NumericVector>(poslist);
     List outobj;
     int j=0;

while(j < pos.length()){
if(j==pos.length()-1){break;}
  if(pos[j]+step>=pos[j+1]){
      vector<int> tmp;
      tmp.push_back(j+1);
      tmp.push_back(j+2);    
 int mean = (pos[j]+pos[j+1])/2;
  j++;
if(j==pos.length()-1){break;}
  while(mean+step>=pos[j+1]){
    j++;
    mean = (mean+pos[j])/2;
    tmp.push_back(j+1);
if(j==pos.length()-1){break;}
  }
  outobj.push_back(tmp);
 }
  j++;

}

return wrap(outobj);
}
