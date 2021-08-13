#include "matOperator.h"
#include "mkl.h"  // mkl.h must appear after inclusion of def.h
#include <stdio.h>
#include <iostream>

void MatOpt::linspace(double x0, double x1, int N, double* mat){
  double dx = (x1 - x0) / N;
  mat[0] = x0;
  for(int iter = 1; iter < N; iter++)
    mat[iter] = mat[iter - 1] + dx;
}

void MatOpt::printMat(double* mat, int NR, int NC, std::string str){
  std::cout << str << std::endl;
  for(int iterRow = 0; iterRow < NR; iterRow++){
    for(int iterCol = 0; iterCol < NC; iterCol++)
      printf("%12.2e", mat[iterRow * NR + iterCol]);
    printf("\n");
  }
}

//double MatOpt::findEigenVal(){
//}

template<>
double MatOpt::findDet<double>(double* kMat, MKL_INT size){

  size_t ul_size = (size_t) size;
  MKL_INT  lda {size}, sizeMKL {size};
  MKL_INT *ipiv = new MKL_INT[size];
  double det {1.0};

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, sizeMKL, sizeMKL, kMat, lda, ipiv);

  for (size_t idx = 0; idx < ul_size; idx++) {
    if((size_t) ipiv[idx] != idx + 1)
      det *= - kMat[idx * ul_size + idx];
    else
      det *= kMat[idx * ul_size + idx];
  }

  delete[] ipiv;

  return det;
}

template<>
MKL_Complex16 MatOpt::findDet<MKL_Complex16>(MKL_Complex16* kMat, MKL_INT size){

  size_t ul_size = (size_t) size;
  MKL_INT  lda {size}, sizeMKL {size};
  MKL_INT *ipiv = new MKL_INT[size];
  MKL_Complex16 det {1.0};

  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, sizeMKL, sizeMKL, kMat, lda, ipiv);

  for (size_t idx = 0; idx < ul_size; idx++) {
    if((size_t) ipiv[idx] != idx + 1)
      det *= - kMat[idx * ul_size + idx];
    else
      det *= kMat[idx * ul_size + idx];
  }

  delete[] ipiv;

  return det;
}

int MatOpt::findOneEigenFunc(int n, double* a, double wIn, double* vr){
  /*
    Find the right eigenvector of the particular real eigenvalue. Note that the input matrix a will be modified.
    Input:
    n: the size of target matrix a (n * n)
    a(n*n): target matrix
    wIn: Eigenvalue
    vr(n): Eigenvector.

    Ref:
    https://software.intel.com/en-us/mkl-developer-reference-c-nonsymmetric-eigenvalue-problems-lapack-computational-routines
    https://software.intel.com/en-us/mkl-developer-reference-c-hsein#D8D6221D-7621-40C5-B494-3AA71854D6E7
   */
  // All the following variable names consist with the ref.
  MKL_INT ilo {1}, ihi {n};
  MKL_INT m = 1, ldh = n, mm = m, ifaill[mm], ifailr[mm];
  MKL_INT ldvl = 1, ldvr = mm;
  double* tau = new double[n];
  lapack_logical select[n] {0};
  double vl[ldvl * n] {0};
  double wr[n] {0}, wi[n] {0};
  MKL_INT nMKL = (MKL_INT) n;

  wr[0] = wIn;

  select[0] = 1;

  LAPACKE_dgehrd(LAPACK_ROW_MAJOR, n, ilo, ihi, a, n, tau);
  MKL_INT info = LAPACKE_dhsein(LAPACK_ROW_MAJOR, 'R', 'N', 'N', select, nMKL,
    a, ldh, wr, wi, vl, ldvl, vr, ldvr, mm,  &m, ifaill, ifailr);

  MKL_INT lda = n, ldc = 1;
  LAPACKE_dormhr (LAPACK_ROW_MAJOR, 'L', 'T', nMKL, 1, ilo, ihi, a, lda, tau, vr, ldc);
  //double ifailrD[mm];
  //for(int i = 0; i < mm; i++){
  //  ifailrD[i] = ifailr[i];
  //}

  //MatOpt::printMat(ifailrD, 1, mm, "ifailr");
  delete[] tau;
  return info;
}

int MatOpt::findEigenValue(int n, double* a, double *wr, double* wi){
  /*
    Find all eigenvalues. Note that the input matrix a will be modified.
    Input:
    n: the size of target matrix a (n * n)
    a(n*n): target matrix
    wr, wi (n): real part and imaginary part of the eigenvalues

    Ref:
    https://software.intel.com/en-us/mkl-developer-reference-c-nonsymmetric-eigenvalue-problems-lapack-computational-routines
    https://software.intel.com/en-us/mkl-developer-reference-c-hseqr
   */
  // All the following variable names consist with the ref.
  MKL_INT ilo {1}, ihi {n}, ldz {n}, ldh {n};
  double* tau = new double[n];
  double *z = new double[(size_t) n * (size_t) ldz];
  MKL_INT nMKL = (MKL_INT) n;

  LAPACKE_dgehrd(LAPACK_ROW_MAJOR, nMKL, ilo, ihi, a, n, tau);

  // a is the hessenberg matrix and also be named as h in the second reference.
  MKL_INT info = LAPACKE_dhseqr(LAPACK_ROW_MAJOR, 'E', 'N', nMKL, 1, nMKL, a, ldh, wr, wi, z, ldz);

  delete[] tau;
  delete[] z;
  return info;
}

int MatOpt::findEigenV(int n, double* a, double *wr, double* wi, double *vr){
  /*
    Computes the eigenvalues and left and right eigenvectors of a general matrix.
    Input:
    n: the size of target matrix a (n * n)
    a(n*n): target matrix
    wr, wi (n): real part and imaginary part of the eigenvalues
    vr (n*n): eigenvectors. Please see the ref for the detail.

    Ref:
    https://software.intel.com/en-us/onemkl-developer-reference-c-geev
  */
  // All the following variable names consist with the ref.
  MKL_INT nMKL = (MKL_INT) n;
  lapack_int ldvl = nMKL;
  double* vl = new double[(size_t) ldvl * (size_t) n];

  lapack_int info =  LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', nMKL, a, nMKL,
    wr, wi, vl, ldvl, vr, nMKL);

  delete[] vl;
  return info;
}

void MatOpt::printEigenvectors( std::string desc, int n, double* wi, double* v){
  /*
    Print out the eigenvectors of a matrix calculated by findEigenV
    desc: description
    n: size of the matrix
    wi (n): imaginary part of the eigenvalues calculated by findEigenV
    v (n*n): eigenvectors calculated by findEigenV
  */
  size_t i, j, ul_n = (size_t) n;
  size_t ldv = ul_n;
  printf( "\n %s\n", desc.c_str() );
  for( i = 0; i < ul_n; i++ ) {
    j = 0;
    while( j < ul_n ) {
      if( wi[j] == (double) 0.0 ) {
        printf( " %10.6f", v[i*ldv+j] );
        j++;
      } else {
        printf( " (%10.6f,%10.6f)", v[i*ldv+j], v[i*ldv+(j+1)] );
        printf( " (%10.6f,%10.6f)", v[i*ldv+j], -v[i*ldv+(j+1)] );
        j += 2;
      }
    }
    printf( "\n" );
  }
}

void MatOpt::homoEqSolver(int n, double* a, double *wr, double* wi, double *sol, double errorLevel){
  /*
    Find the solution of a homogeneous equations which is represented by matrix a.
    !! Use QR decomposition would be the better choice

    Input:
    n: the size of target matrix a (n * n)
    a(n*n): target matrix
    errorLevel: The tolerent error for eigenvalue. For example, errorLevel = 7, then every eigenvalue that smaller than 10^(-7) is condidered 0.

    Output:
    wr, wi (n): Real part and imaginary part of the eigenvalues
    sol (n): Eigenvectors.
   */

  size_t firstReal = 0, ul_n = (size_t) n;

  double *vr = new double[ul_n * ul_n];
  MatOpt::findEigenV(n, a, wr, wi, vr);
  double mechineError = pow(0.1, errorLevel);
  // Find the position of the minimun real eigenvalue
  for(size_t i = 0; i < ul_n; i++){
    if(wi[i] < mechineError){
      firstReal = i;
      break;
    }
  } // Find the first real eigenvalue
  for(size_t i = firstReal; i < ul_n; i++){
    if(wi[i] < mechineError && fabs(wr[firstReal]) > fabs(wr[i]))
      firstReal = i;
  } // Find the minimum
  for(size_t j = 0; j < ul_n; j++){
    sol[j] = vr[firstReal + j * n];
  }
  delete[] vr;

  // Output the solution
  //printf("\neigenValue: %10.2e + i%.2e", wr[firstReal], wi[firstReal]);
}
