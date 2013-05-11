#include "uh.hh"

std::vector<int> uh(const boolMat & M, const boolMat & T, const std::vector<int> & k) {
  if(M.size2() != T.size1() || M.size1() != k.size()) {
    std::cerr << "Error: incompatible arguments to function uh().\n";
    exit(1);
  }
  int g = T.size2();
  std::vector<int> res(T.size2(), 0);
  #pragma omp parallel for schedule(static) shared(M,T,k,res)
  for(int t=0; t < T.size2(); t++) {
    int i=0;
    for(boolMatConstIt1 it1 = M.begin1(); it1 != M.end1(); it1++) {
      bool uniq=true;
      for(boolMatConstIt2 it2 = it1.begin(); it2 != it1.end(); it2++) {
        if(T(it2.index2(),t) == 0) {
          uniq=false;
          break;
        }
      }
      if(uniq) res[t] += k[i];
      i++;
    }
  }
  return(res);
}
