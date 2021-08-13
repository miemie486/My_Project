// g Factor.
#include "gAndG.h"
#include "iostream"
#include <fstream>

/* Calculate g_vvp, Ref. Stadler, 1991 Eqs.(A.1) and (A.2) */

const REAL spin_sng  = _SPIN_N_; // spin of the single particle in the bra (s)
const REAL s_prm_sng = _SPIN_N_; // spin of the single particle in the ket (s')
const REAL t_sngl    = _ISOSPIN_N_; // isospin of the single particle (t)
unsigned int _KMAX_Gqqx_ {10};

using namespace std;

COMPLEX PiMmntm(COMPLEX q1, COMPLEX q2, REAL x) {
  return sqrt(q1*q1 + q2*q2 + 2.0*x*q1*q2);
}

COMPLEX PiMmntm(COMPLEX q1, REAL q2, REAL x) {
  return sqrt(q1*q1 + q2*q2 + 2.0*x*q1*q2);
}

COMPLEX PiMmntm(REAL q1, COMPLEX q2, REAL x) {
  return sqrt(q1*q1 + q2*q2 + 2.0*x*q1*q2);
}

REAL PiMmntm(REAL q1, REAL q2, REAL x) {
  return sqrt(q1*q1 + q2*q2 + 2.0*x*q1*q2);
}

inline void breakup_serial_idxQ1Q2X(int serial_idx, int NQ, int NX,
  int& idxQ1, int& idxQ2, int& idxX){

  int remainder = serial_idx;
  idxX = remainder % NX;
  remainder = (remainder - idxX)/NX;
  idxQ2 = remainder % NQ;
  idxQ1 = remainder / NQ;
}

void mapGArray(size_t& Index_of_GValue, int idxQ1, int idxQ2, int idxX,
  int NQ, int NX) {

  Index_of_GValue = idxQ1 * NQ * NX + idxQ2 * NX + idxX;
}


// 5th row, (A.2)
inline REAL gvvp5(REAL L1, REAL L2, REAL L, REAL opL, REAL l, REAL f, REAL L2p, REAL L1p, REAL Lp, REAL lp, REAL fp, REAL k ){
  return sixJ(L1, L2, L, opL, l, f)*sixJ(L2p, L1p, Lp, opL, lp, fp)*sixJ(f, L2, opL, fp, L1p, k);
}

inline REAL gvvp4(int L1, int l, int L2p, int lp, int k, int L2, int L1p, int opL, REAL L, REAL Lp){
  REAL sum {0};
  int lim[3] {0}, minf {0}, maxf {0}, minfp {0}, maxfp {0};
  int f, fp;

  lim[0]=abs(L1-l);
  lim[1]=abs(k-L1p);
  lim[2]=abs(opL-L2);
  minf= *max_element(begin(lim),end(lim));

  lim[0]=abs(L2p-lp);
  lim[1]=abs(k-L2);
  lim[2]=abs(opL-L1p);
  minfp= *max_element(begin(lim),end(lim));

  lim[0]=L1+l;
  lim[1]=k+L1p;
  lim[2]=opL+L2;
  maxf= *min_element(begin(lim),end(lim));

  lim[0]=L2p+lp;
  lim[1]=k+L2;
  lim[2]=opL+L1p;
  maxfp= *min_element(begin(lim),end(lim));

  for(f=minf; f<=maxf; f++){
    for(fp=minfp; fp<=maxfp; fp++){
      sum+=CGm0(L1, l, f)*CGm0(L2p, lp, fp)*CGm0(k, L2, fp)*CGm0(k, L1p, f)*gvvp5(L1, L2, L, opL, l, f, L2p, L1p, Lp, lp, fp, k);
    }
  }
  return sum;
}

// Third row, (A2)
inline REAL gvvp3(REAL L, REAL l, REAL S, REAL I, REAL j, REAL Lp, REAL lp, REAL Sp, REAL Ip, REAL jp, REAL L1, REAL L1p, REAL L2, REAL L2p, REAL k, REAL opJ){
  REAL sum {0};
  REAL min[3] {0}, max[3] {0};
  REAL minsb, maxsb, minlb, maxlb;
  REAL opS, opL;

  //All possible value of ornamental penmanship S(opS). These limitations can be found by decomposing 9j symbol into several 6j symbols. Please see for the restriction rules of 6j symbol, 9j symbol and CG. Restriction which involves both opL and opS is included in the for loop below.
  min[0]=abs(spin_sng-Sp);
  min[1]=abs(spin_sng-S);
  min[2]=abs(s_prm_sng-Sp);
  minsb= *max_element(begin(min),end(min));
  //cout<< minsb<< endl;

  max[0]=spin_sng+Sp;
  max[1]=spin_sng+S;
  max[2]=s_prm_sng+Sp;
  maxsb= *min_element(begin(max), end(max));
  //cout<< maxsb<< endl;

  //All possible value of ornamental penmanship L
  min[0]=abs(l-L);
  min[1]=abs(lp-Lp);
  minlb = (min[0] > min[1]) ? min[0] : min[1];

  max[0]=l+L;
  max[1]=lp+Lp;
  maxlb = (max[0] < max[1]) ? max[0] : max[1];

  for(opS=minsb; opS<=maxsb; opS++){
    for(opL=minlb; opL<=maxlb; opL++){
      if(opL+opS>= opJ && abs(opL-opS)<= opJ){ // Restriction both involves opS and opL
        sum+=(2*opL+1)*(2*opS+1)*sixJ(spin_sng, spin_sng, S, spin_sng, opS, Sp)*nineJ(Lp, lp, opL, Sp, s_prm_sng, opS, Ip, jp, opJ)*gvvp4(L1, l, L2p, lp, k, L2, L1p, opL, L, Lp)*nineJ(L, l, opL, S, spin_sng, opS, I, j, opJ);
      }
    }
  }

  return sum;
}

inline REAL gvvp2(REAL L, REAL Lp, REAL L1, REAL L2, REAL L1p, REAL L2p, REAL T, REAL Tp, REAL opT, REAL l, REAL S, REAL I, REAL j, REAL lp, REAL Sp, REAL Ip, REAL jp, REAL k, REAL opJ){

  REAL sum {0};

  sum=sqrt(1.*facArray[int(2*L+1)]*facArray[int(2*Lp+1)]/facArray[int(2*L1)]/facArray[int(2*L2)]/facArray[int(2*L1p)]/facArray[int(2*L2p)])*sixJ(t_sngl, t_sngl, T, t_sngl, opT, Tp)*gvvp3(L, l, S, I, j, Lp, lp, Sp, Ip, jp, L1, L1p, L2, L2p, k, opJ);

  return sum;
}

inline REAL hat(REAL x){return sqrt(2*x+1);}

REAL gvvp(REAL L, REAL Lp, REAL L1, REAL L2, REAL L1p, REAL L2p, REAL T, REAL Tp, REAL opT, REAL l, REAL S, REAL I, REAL j, REAL lp, REAL Sp, REAL Ip, REAL jp, REAL k, REAL opJ){

  REAL sum {1};

  /*
  Compared with (A.19) of Glockle's QM Few-Body (83), the difference in phase
  factor arises from opposite sign in def. of \vec{q}
  */
  // sum= int(l+lp+1)% 2 ==0? 1:-1;  // Stadler '91
  sum = -1; // Glockle's Few-body '83
  sum *= hat(I)*hat(Ip)*hat(j)*hat(jp)*hat(L)*hat(Lp)*hat(S)*hat(Sp)*hat(l)*hat(lp)*hat(T)*hat(Tp)*hat(k)*hat(k)*pow(.5, L1+L2p)*gvvp2(L, Lp, L1, L2, L1p, L2p, T, Tp, opT, l, S, I, j, lp, Sp, Ip, jp, k, opJ);

  return sum;
}

void GArrayGen(REAL facCntr, int NQ, REAL* qn, int NX, REAL *xn, REAL *alpha,
  REAL *alphapp, REAL *GArray){

  GArrayGen_tmplt<REAL>(facCntr, NQ, qn, NX, xn, alpha, alphapp, GArray);
}

void GArrayGen(COMPLEX facCntr, int NQ, REAL* qn, int NX, REAL *xn, REAL *alpha,
  REAL *alphapp, COMPLEX *GArray){

  GArrayGen_tmplt<COMPLEX>(facCntr, NQ, qn, NX, xn, alpha, alphapp, GArray);
}

