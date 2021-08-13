/* Operations for identity - kMat, include finding the eigenvalue (binding energy),  eigenvector (wave function), etc.(3.63) in elster's note will become a homogeneous equation, i.e. (identity - kMat) * Phi = 0 */
#include "kMatGen.h"
// #include "mkl.h"  // mkl.h always appears at last among the included
#include "mkl_lapacke.h"


using namespace std;

// Equation numbers refer to codenote [January 21, 2021]

template <>
REAL kMatrix<REAL>::findBE(REAL energy0, REAL energy1, REAL error, std::string str){
  // findout the binding energy given the initial intavel energy0 to energy1.
  // This task in fact treating the determinent of (identity - kMat) as a function of three body energy, and find out its root. We know that the correct three body energy (binding energy) will have a non-trivial wave function, which means equation 0 = (identity - kMat) phi will have non-trivial solution. And this implies the determinent of (identity - kMat) with binding energy is 0. Error is the difference of the last iterate energy and the new one
  double a, b;

  if(energy0 < energy1){
    a = energy0;
    b = energy1;
  }
  else{
    a = energy1;
    b = energy0;
  }

  IterFindZero = 0;
  time_t now;

  // Newton root finder
  if (str == "Newton") {
    double xNew = (a + b) / 2;

    // verbose
    // cout << "-->Try H3BE = " << xNew << endl << flush;
    printf("-->Try H3BE = %-22.15e\n", xNew);
    fflush(stdout);

    double fNew = findKHmgnsDetFunc(xNew);

    cout << "-->  detK = " << fNew << endl;
    now = time(0);
    cout << _TIMER_COL_ << ctime(&now) << flush;

    // cout << "-->Try H3BE = " << a << endl << flush;
    printf("-->Try H3BE = %-22.15e\n", a);
    fflush(stdout);

    double fa = findKHmgnsDetFunc(a);

    cout << "-->  detK = " << fa << endl;
    now = time(0);
    cout << _TIMER_COL_ << ctime(&now) << flush;

    IterFindZero = IterFindZero + 2;

    while((xNew - a > error) && IterFindZero < _MAX_ITER_){
      if(fNew * fa > 0) {
        a = xNew;
        fa = fNew;
      } else {
        b = xNew;
      }
      xNew = (a + b) / 2;

      // verbose
      // cout << "-->Try H3BE = " << xNew << endl << flush;
      printf("-->Try H3BE = %-22.15e\n", xNew);
      fflush(stdout);

      fNew = findKHmgnsDetFunc(xNew);

      cout << "-->  detK = " << fNew << endl;
      now = time(0);
      cout << _TIMER_COL_ << ctime(&now) << flush;
    }

    if (IterFindZero >= _MAX_ITER_) {
      cout << "kMatrix<REAL>::findBE : iteration = " << IterFindZero << endl;
    }
    return xNew;

  } else { // Secant root finder

    double fa, fb, x;

    // verbose
    // cout << "-->Try H3BE = " << a << endl << flush;
    printf("-->Try H3BE = %-22.15e\n", a);
    fflush(stdout);

    fa = findKHmgnsDetFunc(a);

    // cout << "-->  detK = " << fa << endl;
    printf("-->  detK = %-22.15e\n", fa);
    now = time(0);
    cout << _TIMER_COL_ << ctime(&now) << flush;
    IterFindZero++;

    while( (fabs(b - a) > error) && IterFindZero < _MAX_ITER_){

      // verbose
      // cout << "-->Try H3BE = " << b << endl << flush;
      printf("-->Try H3BE = %-22.15e\n", b);
      fflush(stdout);

      fb = findKHmgnsDetFunc(b);

      // cout << "-->  detK = " << fb << endl;
      printf("-->  detK = %-22.15e\n", fb);
      now = time(0);
      cout << _TIMER_COL_ << ctime(&now) << flush;
      IterFindZero++;

      x = b - fb * (b - a) / (fb - fa);
      if(x >= 0) x = -1;
      a = b;
      b = x;
      fa = fb;
    }

    if (IterFindZero >= _MAX_ITER_) {
      cout << "kMatrix<REAL>::findBE : iteration = " << IterFindZero << endl;
    }
    return a;
  }
}

template <>
REAL kMatrix<COMPLEX>::findBE(REAL energy0, REAL energy1, REAL error, std::string str){

  double a, b;

  if(energy0 < energy1){
    a = energy0;
    b = energy1;
  }
  else{
    a = energy1;
    b = energy0;
  }

  // Newton root finder
  if(str == "Newton"){
    double xNew = (a + b) / 2;
    COMPLEX fNew_cmplx = findKHmgnsDetFunc(xNew);
    double fNew = fNew_cmplx.real();
    double fa = findKHmgnsDetFunc(a).real();

    while(xNew - a > error){
      if(fNew * fa > 0){
        a = xNew;
        fa = fNew;
      } else {
        b = xNew;
      }
      xNew = (a + b) / 2;
      fNew_cmplx = findKHmgnsDetFunc(xNew);
      fNew = fNew_cmplx.real();
      // verbose
      cout << "-->Try H3BE = " << xNew << ";    detK = " << fNew_cmplx << endl;
    }
    return xNew;
  }
  // Secant root finder
  else{
    COMPLEX fb_cmplx;
    double fa = findKHmgnsDetFunc(a).real();
    double fb, x;

    while(abs(b - a) > error){
      fb_cmplx = findKHmgnsDetFunc(b);
      // verbose
      cout << "-->Try H3BE = " << b << ";    detK = " << fb_cmplx << endl;
      fb = fb_cmplx.real();
      x = b - fb * (b - a) / (fb - fa);
      if(x >= 0) x = -1;
      a = b;
      b = x;
      fa = fb;
    }
    return a;
  }
}

// Find zero of Re(detK(z)), given imaginary part of z, z = (x, ImE) and return x.
template <>
REAL kMatrix<COMPLEX>::findZeroReDetK(REAL ReE0, REAL ReE1, REAL ImE, REAL error,
  COMPLEX &det) {

  double real_a, real_b;

  if(ReE0 < ReE1){
    real_a = ReE0;
    real_b = ReE1;
  } else{
    real_a = ReE1;
    real_b = ReE0;
  }

  COMPLEX fa_cmplx, fb_cmplx, a, b;

  a = complex<double>(real_a, ImE);
  b = complex<double>(real_b, ImE);
  double fa, fb, x;

  // verbose
  printf("-->Try H3BE = (%-17.10e, %-17.10e)\n", a.real(), a.imag());
  fa_cmplx = findKHmgnsDetFunc(a);
  IterFindZero++;
  fa = fa_cmplx.real();
  printf("--------> detK = (%-17.10e, %-17.10e)\n", fa_cmplx.real(), fa_cmplx.imag());

  // verbose
  time_t now = time(0);
  cout << _TIMER_COL_ << ctime(&now) << flush;
  fflush(stdout);

  while(abs(b.real() - a.real()) > error && IterFindZero < _MAX_ITER_) {

    // verbose
    printf("-->Try H3BE = (%-17.10e, %-17.10e)\n", b.real(), b.imag());
    fb_cmplx = findKHmgnsDetFunc(b);
    IterFindZero++;
    printf("--------> detK = (%-17.10e, %-17.10e)\n", fb_cmplx.real(), fb_cmplx.imag());

    // verbose
    now = time(0);
    cout << _TIMER_COL_ << ctime(&now) << flush;
    fflush(stdout);

    fb = fb_cmplx.real();
    x = b.real() - fb * (b.real() - a.real()) / (fb - fa);

    // debug
    // if(x > 0.0 || x < -1.0e3) x = -2.3;

    a = b;
    b = complex<double>(x, ImE);
    fa = fb;
  }

  det = fb_cmplx;
  return a.real();

}

template <>
COMPLEX kMatrix<COMPLEX>::findResPole(REAL ReE0, REAL ReE1,
  REAL ImE0, REAL ImE1, REAL errorRe, REAL errorIm, COMPLEX &finalDetK) {

  COMPLEX detK;
  double f0, f1, xzero, prevx, y, y0, y1, deltaE = 1.0e-3;

  IterFindZero = 0;

  y0 = ImE0;
  xzero = findZeroReDetK(ReE0, ReE1, y0, errorRe, detK);
  f0 = detK.imag();
  y1 = ImE1;
  prevx = xzero;

  while(abs(y0 - y1) > errorIm && IterFindZero < _MAX_ITER_){

    xzero = findZeroReDetK(prevx, prevx + deltaE, y1, errorRe, detK);
    f1 = detK.imag();
    y = y1 - f1 * (y1 - y0) / (f1 - f0);

    prevx = xzero;
    y0 = y1;
    y1 = y;
    f0 = f1;
  }

  if (IterFindZero >= _MAX_ITER_) {
    cout << "-->kMatrix::findResPole : iteration = " << IterFindZero << endl;
  }

  finalDetK = detK;
  return complex<double>(prevx, y0);
}

// Solve for tilde_T in Eq.(18)
template <>
void kMatrix<COMPLEX>::gen_tilde_T(COMPLEX energyLoc, REAL q0, REAL  input_lamda, REAL input_I){

  show_message_timestamp_fdv("-->kMatrix::gen_tilde_T: begins", _VBSLVL_HIGH_);

  delete [] tilde_T;
  tilde_T = new COMPLEX[dimKMat];

  show_message_timestamp_fdv(
    "---->kMatrix::gen_tilde_T: starts to call linear equation solver", _VBSLVL_HIGH_);

  MKL_INT lda = dimKMat, nrhs = 1, ldb = 1, info;
  MKL_INT ipiv[dimKMat];
  info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, dimKMat, nrhs, kMat, lda, ipiv, tNd_tildeT, ldb);

  show_message_timestamp_fdv(
    "---->kMatrix::gen_tilde_T: finished calling linear equation solver", _VBSLVL_HIGH_);

  if (info > 0) {
    cout << "kMatrix::gen_tilde_T: LAPACK_zgesv reports failure" << endl;
    exit(-1);
  } else {
    for (size_t i = 0; i < dimKMat; i++)
      tilde_T[i] = tNd_tildeT[i];
  }

  show_message_timestamp_fdv("-->kMatrix::gen_tilde_T: ends", _VBSLVL_HIGH_);
}

// Calculate TAMP_omega [Eq.(32)]
template <>
void kMatrix<COMPLEX>::genTAMP_omega(COMPLEX E3CM, REAL mN,
  REAL q0_p, REAL q0, REAL  input_lamda, REAL  input_I) {

  if (! flag_mNE_array_set) {
    show_message_fdv("ERROR (kMatrix::genTAMP_omega): mNE_qr_array not ready",
      _VBSLVL_LOW_);
    return;
  }
  if (! flag_omegas_array_set) {
    show_message_fdv("ERROR (kMatrix::genTAMP_omega): omega_2_array not ready",
      _VBSLVL_LOW_);
    return;
  }
  show_message_timestamp_fdv("-->kMatrix::genTAMP_omega: begins", _VBSLVL_HIGH_);

  delete [] TAMP_omega;
  TAMP_omega = new COMPLEX[nChannels*NQ*NX];

  // show_message_timestamp_fdv(
  //   "---->kMatrix::genTAMP_omega: starts to build two_body tMat", _VBSLVL_HIGH_);
  // tMatrix<COMPLEX>* ptr_tMat[_MAX_3NCHN_] {NULL};
  // for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha ++) {
  //   ptr_tMat[idx_alpha] = tmat_2body(channels[idx_alpha], NQ, mNE_qr_array,
  //     NQ*NX, omega_2_array, NP, PF_pn);
  // }
  // show_message_timestamp_fdv(
  //   "---->kMatrix::genTAMP_omega: finished building two_body tMat", _VBSLVL_HIGH_);

  // generate int_omega
  COMPLEX *int_omega = new COMPLEX[NQ*NX*nChannels] {0.0};

  for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha++) {

    show_message_timestamp_fdv(
      "---->kMatrix::genTAMP_omega: starts to build two_body tMat", _VBSLVL_HIGH_);
    tMatrix<COMPLEX>* ptr_tMat = tmat_2body(channels[idx_alpha], NQ, mNE_qr_array,
      NQ*NX, omega_2_array, NP, PF_pn);
    show_message_timestamp_fdv(
      "---->kMatrix::genTAMP_omega: finished building two_body tMat", _VBSLVL_HIGH_);

    REAL alpha[_N_QMN_CHN_];
    copy_3Nchn(channels[idx_alpha], alpha);
    REAL lprime[2] = {alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]};
    size_t nch_2body = ptr_tMat->readNCH();

    #pragma omp parallel for
    for (size_t idx_rh = 0; idx_rh < NQ*NX; idx_rh++) {
      size_t r, h;
      mapSerialIdxToQX(idx_rh, r, h);
      COMPLEX subsum {0.0};
      for (size_t idx_alphapp = 0; idx_alphapp < nChannels; idx_alphapp++) {
        for (size_t idx_lprime = 0; idx_lprime < nch_2body; idx_lprime++) { // Loop of l prime
          for (size_t idx_ni = 0; idx_ni < NQ*NP; idx_ni++) {
            size_t n, i;
            mapSerialIdxToPQ(idx_ni, i, n);
            COMPLEX tMatFactor =
              ptr_tMat->readTMat_Mat(
                (int) alpha[_CHNIDX_l_], (int) lprime[idx_lprime],
                idx_rh, i, r, alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]);
            /* Skip if the tMat element is smaller than _tMat_ZERO_ */
            // if ( abs(tMatFactor) < _tMat_ZERO_ ) continue;
            COMPLEX mpart {0.0};
            for(size_t m = 0; m < NP; m++) {
              size_t idx_alphabar = (idx_lprime == 0 ? idx_alpha : alphabarArray[idx_alpha]);
              mpart += intzTable_Inhmgns[serialIdxKMat(idx_alphabar, i, r, idx_alphapp, m, n)]
                *tilde_T[serialIdx_by_alphaPQ(idx_alphapp, m, n)];
            }
            subsum += pow(facCntr, 3)*wqn[n]*qn[n]*qn[n]*tMatFactor*mpart;
          }
        }
      }
      int_omega[serialIdx_by_alphaQX(idx_alpha, r, h)] = subsum;
    }
    delete ptr_tMat;
  }

  // for (size_t idx1d_chn = 0; idx1d_chn < nChannels*nChannels; idx1d_chn++) {
  //   size_t idx_alpha = idx1d_chn / nChannels;
  //   size_t idx_alphapp = idx1d_chn % nChannels;

  //   REAL alpha[_N_QMN_CHN_];
  //   copy_3Nchn(channels[idx_alpha], alpha);
  //   REAL lprime[2] = {alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]};
  //   size_t nch_2body = ptr_tMat[idx_alpha]->readNCH();

  //   for (size_t idx_rh = 0; idx_rh < NQ*NX; idx_rh++) {
  //     size_t r, h;
  //     mapSerialIdxToQX(idx_rh, r, h);
  //     for (size_t idx_lprime = 0; idx_lprime < nch_2body; idx_lprime++){ // Loop of l prime
  //       for (size_t idx_ni = 0; idx_ni < NQ*NP; idx_ni++){
  //         size_t n, i;
  //         mapSerialIdxToPQ(idx_ni, i, n);
  //         COMPLEX tMatFactor =
  //           ptr_tMat[idx_alpha]->readTMat_Mat(
  //             (int) alpha[_CHNIDX_l_], (int) lprime[idx_lprime],
  //             idx_rh, i, r, alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]);
  //         /* Skip if the tMat element is smaller than _tMat_ZERO_ */
  //         // if ( abs(tMatFactor) < _tMat_ZERO_ ) continue;
  //         COMPLEX mpart {0.0};
  //         for(size_t m = 0; m < NP; m++) {
  //           size_t idx_alphabar = (idx_lprime == 0 ? idx_alpha : alphabarArray[idx_alpha]);
  //           mpart += intzTable_Inhmgns[serialIdxKMat(idx_alphabar, i, r, idx_alphapp, m, n)]
  //             *tilde_T[serialIdx_by_alphaPQ(idx_alphapp, m, n)];
  //         }
  //         int_omega[serialIdx_by_alphaQX(idx_alpha, r, h)]
  //           += pow(facCntr, 3)*wqn[n]*qn[n]*qn[n]*tMatFactor*mpart;
  //       }
  //     }
  //   }
  // }

  for(size_t idx = 0; idx < nChannels*NQ*NX; idx++)
    TAMP_omega[idx] = tNd_omega[idx] + int_omega[idx];

  // for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha ++) {
  //   delete ptr_tMat[idx_alpha];
  // }

  delete[] int_omega;

  show_message_timestamp_fdv("-->kMatrix::genTAMP_omega: ends", _VBSLVL_HIGH_);
}

