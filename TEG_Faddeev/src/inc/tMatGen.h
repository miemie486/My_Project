#ifndef TMATGEN
#define TMATGEN

#include <stdio.h>
#include "def.h"
#include "gaussIntegral.h"

extern "C"{
  void get_FOST_NCH_CPP(int* L, int* S, int* J, int* NCH);
}

REAL infIntegral(REAL ( *func ) ( REAL x, REAL *), size_t meshNum, REAL mambda, REAL *para = NULL);

REAL inILambda(REAL l, REAL energy, REAL lambda);

void newTMat(REAL intL, REAL intS, REAL intJ, size_t NQ, REAL *eArray, REAL mambda, size_t NP, REAL *pn, REAL *wpn, REAL *tMat, REAL lambda = 1000);

#endif /* TMATGEN */
