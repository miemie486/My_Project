#include "tMat.h"
#include <iomanip>

using namespace std;

void cnvrtCmplxToComplex(Cmplx& _strct, COMPLEX& _obj) {
  _obj = COMPLEX(_strct.dr, _strct.di);
}

void cnvrtComplexToCmplx(COMPLEX& _obj, Cmplx& _strct) {
  _strct.dr = _obj.real();
  _strct.di = _obj.imag();
}

/* Sample code for using tMatrix::chengduT
...
...
// Obtain t_{l, lp}(f*pn[k], f*pn[i], mNElst[r]) where pn[] is REAL array
tMatrix<COMPLEX> tMat(l, s, j, NmE, NP); // NmE = length of mNElst; NP = length of pn[]
tMat.setFacContour(facCntr); // input: facCntr = f
tMat.copyMEArray(mNElst); // input: mNElst[r] = two-body CM energy
tMat.setPotName(potName);
tMat.loadPotPara(numPara2B, para_2B);
tMat.chengduT(pn, wpn); // input: pn[] and wpn[] are REAL mesh points and weights on [0, \infty)
...
...
// alpha[] has all the quantum numbers for the current channel
// To be consistent, l must equal alpha[_CHNIDX_l_]
tval = tMat.readTMat_Mat(l, lp, k, i, r, alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]);

*/

template <>
void tMatrix<REAL>::chengduT(REAL *pn, REAL *wpn){ // Read potential from nnscat
  // Convert string to char array.
  int NmEint = (int) NmE;
  int NPInt = (int) NPK;
  int lenPotName = potName.length();
  char potNameArr[lenPotName];
  strcpy(potNameArr, potName.c_str());

  // For REAL, ignore facCntr
  get_tmtrx_chengdu_C(&intL, &intS, &intJ, &NmEint, mEArray,
    &NPInt, pn, wpn, potNameArr, &lenPotName, potPara, tMStorage, &NCH);

}

template <>
void tMatrix<COMPLEX>::chengduT(REAL *pn, REAL *wpn){
  int NmEint = (int) NmE;
  int NPInt = (int) NPK;
  int lenPotName = potName.length();
  char potNameArr[lenPotName];
  strcpy(potNameArr, potName.c_str());

  size_t sizeT = (size_t) NCH * (size_t) NPInt * (size_t) NCH
    * (size_t) NPInt * (size_t) NmEint;
  Cmplx *tMatCmplx = new Cmplx[sizeT];
  Cmplx *mEArrayCmplx = new Cmplx[(size_t) NmEint];
  Cmplx _facCntrCmplx;

  cnvrtComplexToCmplx(facCntr, _facCntrCmplx);
  for (int i = 0; i < NmEint; i++)
    cnvrtComplexToCmplx(mEArray[i], mEArrayCmplx[i]);
  get_tmtrx_cmplx_chengdu_C(&intL, &intS, &intJ, &NmEint, mEArrayCmplx,
    &NPInt, &_facCntrCmplx, pn, wpn, potNameArr, &lenPotName, potPara, tMatCmplx, &NCH);
  for (size_t i = 0; i < sizeT; i++)
    cnvrtCmplxToComplex(tMatCmplx[i], tMStorage[i]);

  delete[] tMatCmplx;
  delete[] mEArrayCmplx;

}

/* Sample code for using tMatrix<CR>::chengduT_generic()
...
...
// Obtain t_{l, lp}(pk[k], pi[i], mNElst[r])
// where pk[] and pi[] are arbitrary complex array and their lengths are *not*
// necessarily identical
tMatrix<COMPLEX> tMat(l, s, j, NmE, NPK, NPI);
  // NmE = length of mNElst, NPK = length of pk[], NPI = length of pi[]
tMat.setFacContour(facCntr); // input: facCntr = f
tMat.copyMEArray(mNElst); // input: mNElst[r] = two-body CM energy
tMat.setPotName(potName);
tMat.loadPotPara(numPara2B, para_2B);
tMat.chengduT_generic(Nmesh, pn, wpn, pk, pi);
  // input:
  // pn[] and wpn[] are REAL mesh points and weights on [0, \infty)
  // Nmesh is the number of mesh points
  // pn[] and wpn[] are auxiliary data structures used to evaluate
  // t_{l, lp}(pk[k], pi[i], mNElst[r]). pn[] is independent of pk[] or pi[]
  // input: pk[] = array of p', pi[] = array of p
...
...
// alpha[] has all the quantum numbers for the current channel
// To be consistent, l must equal alpha[_CHNIDX_l_]
tval = tMat.readTMat_Mat(l, lp, k, i, r, alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]);

*/

template <>
void tMatrix<COMPLEX>::chengduT_generic(size_t Nmesh, REAL *pn, REAL *wpn, COMPLEX *pk, COMPLEX *pi){
  /*
  Input:
  Nmesh, pn[], wpn[]: mesh pts and weights are to be used by get_gnr_tmtrx_cmplx_chengdu_C
  pk[], pi[]: final and initial complex momenta of t-matrix. Warning
  Output:
  tMStorage[SerialIdx_tMStrg(k, i, r, a, b)] = t_{ab}(pk[k], pi[i], mEarray[r])
  */

  int NmEint = (int) NmE;
  int NPKint = (int) NPK;
  int NPIint = (int) NPI;
  int Nmesh_int = (int) Nmesh;
  int lenPotName = potName.length();
  char potNameArr[lenPotName];
  strcpy(potNameArr, potName.c_str());

  size_t sizeT = (size_t) NCH * (size_t) NPKint * (size_t) NCH
    * (size_t) NPIint * (size_t) NmEint;
  Cmplx *tMatCmplx = new Cmplx[sizeT];
  Cmplx *mEArrayCmplx = new Cmplx[(size_t) NmEint];
  Cmplx *pkArrayCmplx = new Cmplx[(size_t) NPKint];
  Cmplx *piArrayCmplx = new Cmplx[(size_t) NPIint];
  Cmplx _facCntrCmplx;

  cnvrtComplexToCmplx(facCntr, _facCntrCmplx);
  for (int i = 0; i < NmEint; i++)
    cnvrtComplexToCmplx(mEArray[i], mEArrayCmplx[i]);
  for (int i = 0; i < NPKint; i++)
    cnvrtComplexToCmplx(pk[i], pkArrayCmplx[i]);
  for (int i = 0; i < NPIint; i++)
    cnvrtComplexToCmplx(pi[i], piArrayCmplx[i]);
  get_gnr_tmtrx_cmplx_chengdu_C(&intL, &intS, &intJ,
    &NmEint, mEArrayCmplx, &NPKint, pkArrayCmplx, &NPIint, piArrayCmplx,
    &Nmesh_int, &_facCntrCmplx, pn, wpn, potNameArr, &lenPotName, potPara, tMatCmplx, &NCH);

  for (size_t i = 0; i < sizeT; i++)
    cnvrtCmplxToComplex(tMatCmplx[i], tMStorage[i]);

  delete[] tMatCmplx;
  delete[] mEArrayCmplx;

}

/* Sample code for using tMatrix<COMPLEX>::getChengduH2wf
...
...
// Obtain wfval_{l}(plst[i])
// where plst[] is an arbitrary complex array
tMatrix<COMPLEX> h2WaveFunc(Nwfp);
  // Nwfp = length of p[]
h2WaveFunc.setFacContour(facCntr); // input: facCntr = f
h2WaveFunc.setPotName(potName);
h2WaveFunc.loadPotPara(numPara2B, para_2B);
h2WaveFunc.getChengduH2wf(Nmesh, pn, wpn, plst);
  // input:
  // pn[] and wpn[] are REAL mesh points and weights on [0, \infty)
  // Nmesh is the number of mesh points
  // pn[] and wpn[] are auxiliary data structures used to evaluate wave function
  // S- and D-wave wf values at plst[] (complex) will be evaluated and stored
  // internally
...
...
BE = h2WaveFunc.readH2BE();
varphi = h2WaveFunc.readH2wfVal(l, idx_p);
  // return w.f. value for S-wave (l = 0) or D-wave (l = 2) at plst[idx_p]
*/

template <>
void tMatrix<COMPLEX>::getChengduH2wf(size_t Nmesh, REAL *pn, REAL *wpn, COMPLEX *_pArray){

  int lenPotName = potName.length();
  char potNameArr[lenPotName];
  strcpy(potNameArr, potName.c_str());

  Cmplx *Cmsh  = new Cmplx[(size_t) Nmesh];
  Cmplx *Cwght = new Cmplx[(size_t) Nmesh];
  Cmplx *plst  = new Cmplx[(size_t) Nwfp];
  Cmplx *wfvlst  = new Cmplx[2*(size_t) Nwfp];
  Cmplx _H2BE;
  COMPLEX tmpc;

  for (size_t i = 0; i < Nmesh; i++) {
    tmpc = facCntr*pn[i];
    cnvrtComplexToCmplx(tmpc, Cmsh[i]);
    tmpc = facCntr*wpn[i];
    cnvrtComplexToCmplx(tmpc, Cwght[i]);
  }
  for (size_t i = 0; i < Nwfp; i++) {
    cnvrtComplexToCmplx(_pArray[i], plst[i]);
  }

  int Nmesh_int = (int) Nmesh;
  int Nwfp_int  = (int) Nwfp;
  get_generic_h2wfvals_cmplx_chengdu_C(&Nmesh_int, Cmsh, Cwght, &Nwfp_int, plst,
    potNameArr, &lenPotName, potPara, &_H2BE, wfvlst);
  for (size_t i = 0; i < 2*Nwfp; i++) {
    cnvrtCmplxToComplex(wfvlst[i], wfvArray[i]);
  }
  cnvrtCmplxToComplex(_H2BE, H2BE);

  stringstream ss;
  ss << "-->tMatrix::getChengduH2wf: H2BE = " << setprecision(8) << H2BE;
  show_message_fdv(ss.str(), _VBSLVL_HIGH_);

  wfv_set = true;
  delete[] Cmsh;
  delete[] Cwght;
  delete[] plst;
  delete[] wfvlst;

}
