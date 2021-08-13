#ifndef SPLINE
#define SPLINE

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "def.h"


void genSplineCoefs(REAL *coe1, REAL *coe2, REAL *coe3, REAL *x, int meshNum);

//void spline(REAL *x, REAL var, int meshNum, REAL *Si);
void splineArray(REAL *x, REAL var, int meshNum, REAL *Si, REAL *coe1, REAL *coe2, REAL *coe3);
REAL spline(REAL *x, REAL var, int meshNum, REAL *coe1, REAL *coe2, REAL *coe3, int iRequire);
void spline2(REAL *x, REAL var, int meshNum, REAL *Si);

class splineFunc{
  size_t n, ni; // n: num of given points  ni: num of requiring points
  REAL *xi, *x;
  REAL *yi, *y;
  REAL *coe1;
  REAL *coe2;
  REAL *coe3;

 public:
  REAL (*func)(REAL);

  splineFunc(size_t, size_t);
  ~splineFunc();

  inline void genCoe(void){genSplineCoefs(coe1, coe2, coe3, x, n);}
  inline double si(int i, double p) {return spline(x, p, n, coe1, coe2, coe3, i);}
  inline REAL& cx(size_t i) {if(i > n) printf("ERROR: i out range"); return x[i];}
  inline REAL& cy(size_t i) {if(i > n) printf("ERROR: i out range"); return y[i];}
  inline REAL& cxi(size_t i) {if(i > ni) printf("ERROR: i out range"); return xi[i];}
  inline REAL& cyi(size_t i) {if(i > ni) printf("ERROR: i out range"); return yi[i];}

  void printArr(REAL*array, size_t len);
  void genYFun(void);
  void genYIFun(void);
  void plotTable(std::ofstream&);
};

#endif /* SPLINE */
