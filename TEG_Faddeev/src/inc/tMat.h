#ifndef TMAT
#define TMAT

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <sstream>

#include "def.h"
#include "gaussIntegral.h"
#include "msg.h"

#define NUM_PARA_TWOBODY 10  // Fixed number of callback two-body parameters

typedef struct {REAL dr; REAL di;} Cmplx;

extern "C"{
  // Initialize global variables used by nnscat
  void init_nnscat_chengdu_C(void);

  // void get_FOST_tmtrx_DLSPR(int* L, int* S, int* J, int* N_mE, double* mE_lst, double* Mambda, int* Nmesh, double* msh, double*  wght, double* tm_lst, int* NCH);
  void get_FOST_NCH_CPP(int* L, int* S, int* J, int* NCH);
  void get_tmtrx_chengdu_C(int* L, int* S, int* J, int* N_mE, double* mE_lst,
    int* Nmesh, double* msh, double* wght, char* pot_name, int* lenPotName,
    REAL* extrapara, double* tm_lst, int* NCH1);

  /*extrapara[0] = uptoQn
    extrapara[i] = x_i and so on
    => V(x1, x2, ...) = V^(0) + \sum_i x_i * V^(i)
  */
  void get_tmtrx_cmplx_chengdu_C(int* L, int* S, int* J, int* N_mE, Cmplx* mE_lst,
    int* Nmesh, Cmplx* _facCntr, double* msh, double* wght,
    char* pot_name, int* lenPotName, REAL* extrapara, Cmplx* tm_lst, int* NCH1);

  /* get_gnr_tmtrx_cmplx_chengdu_C
  Input:
    pp_lst[k], k = 0..(N_pp - 1)
    p_lst[i],  i = 0..(N_p  - 1)
    mE_lst[r], r = 0..(N_mE - 1)
    NCH1 = 1 (uncoupled two-body channel) or 2 (coupled two-body channel)
    The rest is similar to get_tmtrx_cmplx_chengdu_C
  Output:
    tm_lst[SerialIdx_tMStrg(k, i, r, a, b)] = t_{ab}(p_k, p_i, mE_r)
  */
  void get_gnr_tmtrx_cmplx_chengdu_C(int* L, int* S, int* J,
    int* N_mE, Cmplx* mE_lst, int* N_pp, Cmplx* pp_lst, int* N_p, Cmplx* p_lst,
    int* Nmesh, Cmplx* _facCntr, double* msh, double* wght,
    char* pot_name, int* lenPotName, REAL* extrapara, Cmplx* tm_lst, int* NCH1);

  /* get_generic_h2wfvals_cmplx_chengdu_C
    return analytically continued deuteron w.f. values at p_lst[i]
  Input:
    Mesh must be Cmplx!
    p_lst[i], i = 0..(N_p - 1)
    p_lst can, but doesn't have to, be identical to Cmsh
  Output:
    wfval_lst[i] = S-wave w.f. value at p_lst[i]
    wfval_lst[i+N_p] = S-wave w.f. value at p_lst[i]
  */
  void get_generic_h2wfvals_cmplx_chengdu_C(int* Nmesh, Cmplx* Cmsh,
    Cmplx* Cwght, int* N_p, Cmplx* p_lst, char* pot_name, int* lenPotName,
    REAL* extrapara, Cmplx* H2BE, Cmplx* wfval_lst);
}

void cnvrtCmplxToComplex(Cmplx& _strct, COMPLEX& _obj);
void cnvrtComplexToCmplx(COMPLEX& _obj, Cmplx& _strct);

template <class CR>
class tMatrix{
  /*
  Generate and store two-body t-matrices and the deuteron wave function
  */
protected:
  int intL, intS, intJ;
  size_t NmE, NPK, NPI, Nmesh, Nwfp;
  size_t dim_tMStorage;
  // t_{ab}(p_k, p_i; mE_r) with k = [1..NPK] & i = [1..NPI] & r = [1..NmE]
  // stored in 1D array tMatrix.tMStorage
  CR *mEArray = NULL, *tMStorage = NULL, *wfvArray = NULL;
  bool wfv_set = false;
  CR H2BE;
  int NCH;
  std::string potName;
  CR facCntr {1.0}; // contour variable: q -> facCntr * q

public:

  int numPara {1};
  double potPara[NUM_PARA_TWOBODY] {0.0}; //para[0] = uptoQn, 0 by default

  // constructors for generating t-matrices
  tMatrix(int _intL, int _intS, int _intJ, size_t _NmE, size_t _NP);
  tMatrix(int _intL, int _intS, int _intJ, size_t _NmE, size_t _NPK, size_t _NPI);

  // constructor for generating deuteron w.f.
  tMatrix(size_t _Nwfp);
  ~tMatrix();
  tMatrix(const tMatrix &other) {printf("ERROR: COPY IS FORBBIDDEN"); abort();};

  void loadPotPara(int _numPara, REAL* _para) {
    numPara = (NUM_PARA_TWOBODY > _numPara) ? _numPara : NUM_PARA_TWOBODY;
    for (int ii = 0; ii < _numPara; ii++)
      potPara[ii] = _para[ii];
  };

  void setFacContour(CR _facCntr) {
    facCntr = _facCntr;
  };

  CR FacContour() {
    return facCntr;
  };

  void copyMEArray(CR* _mNElst);
  size_t SerialIdx_tMStrg(size_t k, size_t i, size_t r, size_t a, size_t b) const;
  CR readTMat_Mat(int lpin, int lpinp, size_t k, size_t i,
    size_t r, int l_alpha, int lcpld_alpha) const;
  CR readH2wfVal(int l, size_t idx_p) const;
  CR readH2BE() const {
    if (wfv_set) {
      return H2BE;
    } else {
      std::cout << "readH2BE: H2BE not been calculated." << std::endl;
      return 0.0;
    }
  };
  int readNCH() const {return NCH;};
  void setPotName(std::string _potName){potName = _potName;};
  void chengduT(REAL *pn, REAL *wpn);
  void chengduT_generic(size_t Nmesh, REAL *pn, REAL *wpn, CR *pk, CR *pi);
  void getChengduH2wf(size_t Nmesh, REAL *pn, REAL *wpn, CR *_pArray);

};

template <class CR>
tMatrix<CR>::tMatrix(int _intL, int _intS, int _intJ, size_t _NmE, size_t _NP){

  intL = _intL;
  intS = _intS;
  intJ = _intJ;
  NmE = _NmE;
  NPK = _NP;
  NPI = _NP;
  get_FOST_NCH_CPP(&intL, &intS, &intJ, &NCH);
  mEArray = new CR[(size_t) NmE];
  dim_tMStorage = (size_t) NCH * (size_t) NPK * (size_t) NCH * (size_t) NPI * (size_t) NmE;
  tMStorage = new CR[dim_tMStorage];
}

template <class CR>
tMatrix<CR>::tMatrix(int _intL, int _intS, int _intJ, size_t _NmE, size_t _NPK,
  size_t _NPI){

  intL = _intL;
  intS = _intS;
  intJ = _intJ;
  NmE = _NmE;
  NPK = _NPK;
  NPI = _NPI;
  get_FOST_NCH_CPP(&intL, &intS, &intJ, &NCH);
  mEArray = new CR[(size_t) NmE];
  dim_tMStorage = (size_t) NCH * (size_t) NPK * (size_t) NCH * (size_t) NPI * (size_t) NmE;
  tMStorage = new CR[dim_tMStorage];
}

template <class CR>
tMatrix<CR>::tMatrix(size_t _Nwfp) {
  Nwfp = _Nwfp;
  wfvArray = new CR[2*(size_t) Nwfp];
}

template <class CR>
tMatrix<CR>::~tMatrix(){
  delete[] mEArray;
  delete[] tMStorage;
  delete[] wfvArray;
}

template <class CR>
void tMatrix<CR>::copyMEArray(CR* _mNElst) {
  for (size_t i = 0; i < NmE; i++) {
    mEArray[i] = _mNElst[i];
  }
}

/* Location of t_{ab}(p_k, p_i; mE_r) in 1D array tMStorage[], i.e.,
   tMStorage[SerialIdx_tMStrg(k, i, r, a, b)] = t_{ab}(p_k, p_i; mE_r)
   a = 1, 2 and b = 1, 2
*/
template <class CR>
size_t tMatrix<CR>::SerialIdx_tMStrg(size_t k, size_t i, size_t r, size_t a, size_t b) const {

  if ((a > (size_t) NCH) || (b > (size_t)NCH)) {
    std::cout << "tMatrix<CR>::SerialIdx_tMStrg: a > NCH or b > NCH, wrong result returned"
      << std::endl;
  }
  size_t loc = r*NCH*NPI*NCH*NPK + ((b - 1)*NPI + i)*NCH*NPK + (a - 1)*NPK + k;
  if (loc < dim_tMStorage) {
    return loc;
  } else {
    std::cout << "ERROR! tMatrix::SerialIdx_tMStrg: out of range " << loc << " >= " << dim_tMStorage;
    return 0;
  }
}

template <class CR>
CR tMatrix<CR>::readTMat_Mat(int l, int lp, size_t k, size_t i, size_t r,
  int l_alpha, int lcpld_alpha) const {
  /* Read tMat_{l, lp}(p_k, p_i, mEArray[r])

    l: quantum number l in Faddeev eq
    lp: quantum number l p in Faddeev eq
    l_alpha: quantum number l in the current channel
    lcpld_alpha: coupling quantum number l in the current channel
    p_k, p_i, E-3*q_l^2/4/m: see (3.67), Elster_Notes
    NP: number of p P points
    NCH: uncoupled--1, coupled--2
  */
  CR out {0.0};

  if (NCH == 1) {
    out = tMStorage[SerialIdx_tMStrg(k, i, r, 1, 1)];
  } else {
    if (NCH == 2) {
      // Let l_lower be the min{l_alpha, lcpld_alpha}
      int l_lower = (l_alpha < lcpld_alpha) ? l_alpha : lcpld_alpha;
      if (l_lower < 0) {
        std::cout << "Error=> readTMat_Mat: lcpld < 0 for NCH = 2" << std::endl;
      }
      size_t a = (l == l_lower) ? 1 : 2;
      size_t b = (lp == l_lower) ? 1 : 2;
      out = tMStorage[SerialIdx_tMStrg(k, i, r, a, b)];

      // debug
      if (! (SerialIdx_tMStrg(k, i, r, a, b) < dim_tMStorage)) {
        std::cout << "Error => readTMat_Mat: out of range "
          << SerialIdx_tMStrg(k, i, r, a, b) << " >= " << dim_tMStorage
          << std::endl;
      }
    } else {
      printf("readTMat ERROR: NCH = 1 or 2");
    }
  }

  return out;
}

template <class CR>
CR tMatrix<CR>::readH2wfVal(int l, size_t idx_p) const {

  if (! wfv_set) {
    std::cout << "readH2wfVal: w.f. not been built" << std::endl;
    return 0.0;
  }
  if (l == 0) {
    return wfvArray[idx_p];
  } else {
    if (l == 2) {
      return wfvArray[idx_p + Nwfp];
    } else {
      std::cout << "l = " << l << " is not part of the deuteron" << std::endl;
      return 0.0;
    }
  }
}

#endif /* TMAT */
