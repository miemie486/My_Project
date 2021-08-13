#ifndef DELTA
#define DELTA

#include "iostream"
#include "def.h"
#include "3nj.h"
#include "mkl.h"
using namespace std;

void get_Upara_list(REAL opJ_, REAL parity_);

int SerialofIdx_U_list(REAL opJ_, REAL lambda_p, REAL I_p, REAL lambda, REAL I);

//equ.(40) in codenote
COMPLEX U_lambda_sigma(REAL q0, REAL opJ, REAL lambda_prime, REAL sigma_prime, REAL lambda, REAL sigma, COMPLEX * U_array);

COMPLEX S_lambda_sigma(REAL q0, REAL opJ, REAL lambda_prime, REAL sigma_prime, REAL lambda, REAL sigma, COMPLEX * U_array);

//generate the matrix of S_lambda_sigma
void SMatGen(COMPLEX * Smat, int N_Smat, REAL opJ, REAL _parity_, REAL q0, COMPLEX * U_array);

void Gen_deltaArray(COMPLEX * deltaArray_, int N_delta, COMPLEX * Mixangle, int N_mixangle, COMPLEX * Smat_, int N_Smat);

//Reorganize U matrix and the matrix of eigenvalue to satisfy the convention given by Hubber.
void Reorganize_U(int N_line, COMPLEX * U, COMPLEX * Ev);

#endif

