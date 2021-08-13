#include "mathFunc.h"

namespace mathFunc{

  double rootFinder(double (*f) (double), double a, double b, double error){
    //double xNew = (a + b) / 2;
    //double fNew = f(xNew);
    //double fa = f(a);

    //printf("fNew = %.2e\tfOld = %.2e\n", fNew, fa);

    //while(xNew - a > error){
    //  if(fNew * fa > 0){
    //    a = xNew;
    //    fa = fNew;
    //  }
    //  else{b = xNew;}
    //  xNew = (a + b) / 2;
    //  fNew = f(xNew);
    //  printf("xNew = %.2e\n", xNew);
    //}

    // Bisection
    double fa = f(a), fb;
    double x;

    while(fabs(b - a) > error){
      fb = f(b);
      x = b - fb * (b - a) / (fb - fa);
      a = b;
      b = x;
      fa = fb;
    }
    return b;
  }
  double min(double *a, int len){
    return *std::min_element(a, a + len);
  }
  double max(double *a, int len){
    return *std::max_element(a, a + len);
  }
  int minPos(double *a, int len){
    /*
      Find the position of the first minimun element.
    */
    int idx = 0;
    for(int i = 0; i < len; i++){
      if(a[idx] > a[i])
        idx = i;
    }
    return idx;
  }
  int maxPos(double *a, int len){
    /*
      Find the position of the first minimun element.
    */
    int idx = 0;
    for(int i = 0; i < len; i++){
      if(a[idx] < a[i])
        idx = i;
    }
    return idx;
  }
  //double mathFunc::2Dintegral(){};
}
