#ifndef GVVP
#define GVVP

/* Calculate g_vvp, Ref. Stadler, 1991 */
#include <iostream>

#include "3nj.h"
#include "def.h"
#include "channels.h"

/* Rosetta stone of notations for QMN
  Stadler  Elster   codenote      Meaning
   L        l         l          Orbit momentem of the pairs
   T        T         T          Isospin of the pairs
   opT      opT       opT        Total isospin
   l        lambda    lambda     Orbit momentem of the third paticle
   S        s         s          Spin of the pair
   I        j         j          Total momentem of the pair
   j        J         I          Total momentem of the third paticle
                      -1         Coupled channels' pair's orbit angular momentem
                      opT        Total isospin
 */
// notation of Stadler is followed here

extern unsigned int _KMAX_Gqqx_;

COMPLEX PiMmntm(COMPLEX q1, COMPLEX q2, REAL x);
COMPLEX PiMmntm(COMPLEX q1, REAL q2, REAL x);
COMPLEX PiMmntm(REAL q1, COMPLEX q2, REAL x);
REAL PiMmntm(REAL q1, REAL q2, REAL x);

void mapGArray(size_t& Index_of_GValue, int idxQ1, int idxQ2, int idxX,
  int NQ, int NX);

void breakup_serial_idxQ1Q2X(int serial_idx, int NQ, int NX,
  int& idxQ1, int& idxQ2, int& idxX);

REAL gvvp(REAL L, REAL Lp, REAL L1, REAL L2, REAL L1p, REAL L2p, REAL T, REAL Tp, REAL opT, REAL l, REAL S, REAL I, REAL j, REAL lp, REAL Sp, REAL Ip, REAL jp, REAL k, REAL opJ);

template <class CR>
CR Gqqx(REAL L, REAL Lp, REAL T, REAL Tp, REAL opT, REAL l, REAL S, REAL I,
  REAL j, REAL lp, REAL Sp, REAL Ip, REAL jp, CR q, CR qp, REAL x, REAL opJ){
  /* Please see the meaning of all these quantum numbers in channels.h. The
  naming scheme of symbols follows that of Stadler 1991.
  */

  REAL L1 {0}, L2 {L}, L1p {0}, L2p {Lp};
  CR subsum {0}, sum {0}, sumTemp {0};
  REAL pl[2] {0}, err {1.0}, temp, gValue {0.0};
  int k;

  pl[0] = 1;
  pl[1] = x;
  CR pow_q {0}, pow_qp {0};
  for(k = 0; k <= 1; k++){
    L2 = L;
    subsum = 0;
    for(L1 = 0; L1 <= L; L1++){
      L2p = Lp;
      for(L1p = 0; L1p <= Lp; L1p++){
        gValue = gvvp(L, Lp, L1, L2, L1p, L2p, T, Tp, opT, l, S, I, j, lp, Sp,
          Ip, jp, k, opJ);
        if(abs(L1 + L1p) < 1.0e-10) {
          pow_q = pow(q, 0);
        } else {
          pow_q = pow(q, L1 + L1p);
        }
        if(abs(L2 + L2p) < 1.0e-10) {
          pow_qp = pow(qp, 0);
        } else {
          pow_qp = pow(qp, L2+L2p);
        }
        subsum += pow_q*pow_qp*gValue;
        L2p -= 1;
      }
      L2 -= 1;
    }
    sum += subsum*pl[k];
  }

  sumTemp = subsum;
  k = 1;

  /* k normally stops at 6 for _ERROR_Gqqx_ < 1.0e-5 */
  while(err > _ERROR_Gqqx_){
    L2 = L;
    subsum = 0;
    k += 1;

    temp = pl[1];
    pl[1] = ((2*k-1)*x*pl[1]-(k-1)*pl[0])/k;
    pl[0] = temp;

    for(L1 = 0; L1 <= L; L1++){
      L2p = Lp;
      for(L1p = 0; L1p <= Lp; L1p++){
        gValue = gvvp(L, Lp, L1, L2, L1p, L2p, T, Tp, opT, l, S, I, j, lp, Sp,
          Ip, jp, k, opJ);
        if(abs(L1 + L1p) < 1.0e-10) {
          pow_q = pow(q, 0);
        } else {
          pow_q = pow(q, L1 + L1p);
        }
        if(abs(L2 + L2p) < 1.0e-10) {
          pow_qp = pow(qp, 0);
        } else {
          pow_qp = pow(qp, L2 + L2p);
        }
        subsum += pow_q*pow_qp*gValue;
        L2p -= 1;
      }
      L2 -= 1;
    }
    err = abs(sumTemp - subsum)/abs(sumTemp);
    sumTemp = subsum;
    sum += subsum*pl[1];
  }

  return sum;
}

template <class CR>
CR Gqqx(int kmax, REAL L, REAL Lp, REAL T, REAL Tp, REAL opT, REAL l, REAL S,
  REAL I, REAL j, REAL lp, REAL Sp, REAL Ip, REAL jp, CR q, CR qp, REAL x,
  REAL opJ){

  REAL L1 {0}, L2 {L}, L1p {0}, L2p {Lp};
  CR subsum {0}, sum {0};
  REAL pl[2] {0}, temp, gValue {0.0};
  int k;

  pl[0]=1;
  pl[1]=x;

  for(k=0; k<=1; k++){
    L2=L;
    subsum=0;
    for(L1=0; L1<=L; L1++){
      L2p=Lp;
      for(L1p=0; L1p<=Lp; L1p++){
        gValue = gvvp(L, Lp, L1, L2, L1p, L2p, T, Tp, opT, l, S, I, j, lp, Sp,
          Ip, jp, k, opJ);
        subsum+=pow(q, L1+L1p)*pow(qp, L2+L2p)*gValue;
        L2p-=1;
      }
      L2-=1;
    }
    sum+=subsum*pl[k];
  }

  for (k = 2; k<= kmax; k++){

    L2=L;
    subsum=0;

    temp=pl[1];
    pl[1]=((2*k-1)*x*pl[1]-(k-1)*pl[0])/k;
    pl[0]=temp;

    for(L1=0; L1<=L; L1++){
      L2p=Lp;
      for(L1p=0; L1p<=Lp; L1p++){
        gValue = gvvp(L, Lp, L1, L2, L1p, L2p, T, Tp, opT, l, S, I, j, lp, Sp,
          Ip, jp, k, opJ);
        subsum+=pow(q, L1+L1p)*pow(qp, L2+L2p)*gValue;
        L2p-=1;
      }
      L2-=1;
    }
    sum+=subsum*pl[1];
  }

  return sum;
}

// Wrapper for Gqqx
// In: alpha[] and beta[] are 3N quantum number arrays (see channels.h)
// Out: G_{alpha, beta}(q, qp, x) [Eq.A1 of Stadler '91]
template <class CR>
CR Gqqx_3Nchn(REAL *alpha, REAL *beta, CR q, CR qp, REAL x) {

  if (abs(alpha[_CHNIDX_opJ_] - beta[_CHNIDX_opJ_]) > 1.0e-4
    || abs(alpha[_CHNIDX_opT_] - beta[_CHNIDX_opT_]) > 1.0e-4) {
    return 0.0;
  }

  // return Gqqx<CR>(alpha[_CHNIDX_l_], beta[_CHNIDX_l_],
  return Gqqx<CR>(_KMAX_Gqqx_, alpha[_CHNIDX_l_], beta[_CHNIDX_l_],
    alpha[_CHNIDX_T_], beta[_CHNIDX_T_], alpha[_CHNIDX_opT_],
    alpha[_CHNIDX_lmbd_], alpha[_CHNIDX_s_], alpha[_CHNIDX_j_],
    alpha[_CHNIDX_I_], beta[_CHNIDX_lmbd_], beta[_CHNIDX_s_],
    beta[_CHNIDX_j_], beta[_CHNIDX_I_],
    q, qp, x, alpha[_CHNIDX_opJ_]);
}

// G_{alpha, alphapp}(qn[i], qn[j], xn[h]) -> GArray
template <class CR>
void GArrayGen_tmplt(CR facCntr, int NQ, REAL* qn, int NX, REAL *xn, REAL *alpha,
  REAL *alphapp, CR *GArray){

  REAL alphabar[_N_QMN_CHN_];
  copy_3Nchn(alpha, alphabar);

  #pragma omp parallel for
  for (size_t idx_serial = 0; idx_serial < (size_t) NQ*NQ*NX; idx_serial++) {
    int idxQ1, idxQ2, idxX;
    breakup_serial_idxQ1Q2X(idx_serial, NQ, NX, idxQ1, idxQ2, idxX);
    GArray[idx_serial]
      = Gqqx_3Nchn(alphabar, alphapp, facCntr*qn[idxQ1], facCntr*qn[idxQ2], xn[idxX]);
  }

}

void GArrayGen(REAL facCntr, int NQ, REAL* qn, int NX, REAL *xn, REAL *alpha,
  REAL *alphapp, REAL *GArray);

void GArrayGen(COMPLEX facCntr, int NQ, REAL* qn, int NX, REAL *xn, REAL *alpha,
  REAL *alphapp, COMPLEX *GArray);

#endif /* GVVP */
