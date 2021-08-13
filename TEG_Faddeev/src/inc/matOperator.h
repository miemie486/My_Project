#ifndef MATOPT
#define MATOPT

#include <string>
#include <iostream>
#include <tgmath.h>
#include "def.h"
#include "mkl.h"

namespace MatOpt{
  void linspace(double, double, int, double*);
  void printMat(double*, int, int, std::string);
  double findEigenVal();

  template <class CR>
  CR findDet(CR*, MKL_INT);

  // double findDet(double*, int);
  int findOneEigenFunc(int n, double* a, double wIn, double* vr);
  int findEigenValue(int n, double* a, double* wr, double* wi);
  int findEigenV(int n, double* a, double *wr, double* wi, double *vr);
  void printEigenvectors( std::string desc, int n, double* wi, double* v);
  void homoEqSolver(int n, double* a, double *wr, double* wi, double *sol, double errorLevel);
}

#endif /* MATOPT */
