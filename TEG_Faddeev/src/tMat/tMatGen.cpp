/* TMatrix for pionless theory. This is useless */
#include "tMatGen.h"

REAL fLambda(REAL x, REAL lambda){
  return exp(-x * x / lambda / lambda);
}

REAL inILambda(REAL l, REAL energy, REAL lambda){
  return massN * l * l / ( massN * energy - l * l) * exp( - 2 * l * l / ( lambda * lambda )) / 2.0 / M_PI / M_PI;
}

REAL funcinILambda(REAL x, REAL *para){
  return inILambda(x, para[0], para[1]);
}

REAL infIntegral(REAL (*func) ( REAL, REAL *), size_t meshNum, REAL mambda, REAL *para){
  REAL *xn = new REAL[meshNum], *wn = new REAL[meshNum];
  REAL sum {0};

  // createGauss(xn, wn, meshNum, mambda, mambda);
  createGaussInf(xn, wn, meshNum, mambda);
  for(size_t i = 0; i < meshNum; i++) sum += func(xn[i], para) * wn[i];

  return sum;
}

void newTMat(REAL intL, REAL intS, REAL intJ, size_t NQ, REAL *eArray, REAL mambda, size_t NP, REAL *pn, REAL *wpn, REAL *tMat, REAL lambda){
  REAL (*func) (REAL, REAL*) = &funcinILambda;
  REAL ILambda {0}, a {0};
  REAL *para = new REAL[2], *regulator = new REAL[NP];
  REAL C {0};

  if( intS == 0 && intL == 0 && intJ == 0)  // S, L, J
    a = -23.714;
  else if( intS == 1 && intL == 0 && intJ == 1 )
    a = 5.42;
  else {printf("Channel does exist!\n");}

  C = 4 * M_PI / massN * a * _HBARC_;

  para[1] = lambda;

  // createGauss(pn, wpn, NP, mambda, mambda);
  createGaussInf(pn, wpn, NP, mambda);
  for(size_t idx = 0; idx < NP; idx ++){
    regulator[idx] = fLambda(pn[idx], lambda);
    //printf("p[%lu] = %.2e\n", idx, pn[idx]);
    //printf("regulator[%lu] = %.2e\n", idx, regulator[idx]);
  }
  //printf("\n");

  for(size_t idx = 0; idx < NQ; idx ++){  // index of energyArray
    para[0] = eArray[idx];
    ILambda = infIntegral(func, NP, mambda, para);
    for(size_t idy = 0; idy < NP; idy ++) // index of pn1
      for(size_t idz = 0; idz < NP; idz ++){ // index of pn2
        tMat[idx * NP * NP + idy * NP + idz] = regulator[idy] * regulator[idz] / (1.0 / C + ILambda);
      }
  }
}
