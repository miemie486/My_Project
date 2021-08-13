#ifndef KMATINHTPP
#define KMATINHTPP

#include <iomanip>

// Equation numbers refer to codenote [January 21, 2021]

// Pi_1 and Pi_2 arrays
template <class CR>
void kMatrix<CR>::genPiArrays(REAL q0) {

  show_message_fdv("-->kMatrix::genPiArrays: begins", _VBSLVL_HIGH_);

  for (size_t i = 0; i < NQ*NX; i++) {
    size_t r, h;
    mapSerialIdxToQX(i, r, h);
    Pi_1_array[i] = PiMmntm(0.5*facCntr*qn[r], q0, xn[h]);
    Pi_2_array[i] = PiMmntm(facCntr*qn[r], 0.5*q0, xn[h]);
  }
  flag_Pis_array_set = true;

  show_message_fdv("-->kMatrix::genPiArrays: ends", _VBSLVL_HIGH_);
}

template <class CR>
void kMatrix<CR>::gen_omega_arrays(REAL q0_p) {

  show_message_fdv("-->kMatrix::gen_omega_arrays: begins", _VBSLVL_HIGH_);

  for (size_t i = 0; i < NQ*NX; i++) {
    size_t r, h;
    mapSerialIdxToQX(i, r, h);
    omega_1_array[i] = PiMmntm(0.5*q0_p, facCntr*qn[r], xn[h]);
    omega_2_array[i] = PiMmntm(q0_p, 0.5*facCntr*qn[r], xn[h]);
  }
  flag_omegas_array_set = true;
  show_message_fdv("-->kMatrix::gen_omega_arrays: ends", _VBSLVL_HIGH_);
}

// OKed by Bingwei 07/21/20
template <class CR>
void kMatrix<CR>::genIntzTable_Inhmgns(CR energy) {

  if (! GTableAllocFlag) {
    cout << "-->kMatrix::genIntzTable_Inhmgns : Error GqqxTable not yet constructed " << endl;
    return;
  }
  show_message_timestamp_fdv("-->kMatrix::genIntzTable_Inhmgns: begins",
    _VBSLVL_LOW_);

  show_message_fdv("-->kMatrix::genIntzTable_Inhmgns -> about to alloc MEM (GB) "
    + to_string(double(dimKMat*dimKMat*sizeof(CR))/(_SIZE_GB_)),
    _VBSLVL_HIGH_);
  delete[] intzTable_Inhmgns;
  intzTable_Inhmgns = new CR[dimKMat*dimKMat];

  #pragma omp parallel for
  for (size_t idx_ir = 0; idx_ir < NP*NQ; idx_ir++) {
    for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha++) {
      for (size_t idx_beta = 0; idx_beta < nChannels; idx_beta++) {
        size_t i, r;
        mapSerialIdxToPQ(idx_ir, i, r);
        for (size_t m = 0; m < NP; m++) {
          for (size_t n = 0; n < NQ; n++) {
            intzTable_Inhmgns[serialIdxKMat(idx_alpha, i, r, idx_beta, m, n)]
              = intz_Inhmgns(energy, avg_mN,
              channels[idx_alpha][_CHNIDX_l_], channels[idx_alpha],
              channels[idx_beta], GqqxTable[idx_alpha][idx_beta], i, m, r, n);
          }
        }
      }
    }
  }

  show_message_timestamp_fdv("-->kMatrix::genIntzTable_Inhmgns: ends",
    _VBSLVL_LOW_);
}

template <class CR>
CR kMatrix<CR>::intz_Inhmgns(CR E3CM, REAL MASS, REAL lprime,
  REAL alpha[], REAL alphapp[], CR *GArray, size_t& i, size_t& m, size_t& r, size_t& n) {

  CR z {0.0}, sum {0.0}, gFactor {0.0},fctr_rn {0.0};
  size_t index_GValue;
  REAL si {0.0}, sm {0.0};
  REAL pi1, pi2;

  for(size_t iX=0; iX< NX; iX++){
    mapGArray(index_GValue, r, n, iX, NQ, NX);
    gFactor = GArray[index_GValue];
    if( abs(gFactor) > _INTZ_EPS_ ) {
      fctr_rn = 1.0/(E3CM - facCntr*facCntr*(qn[r]*qn[r] + qn[n]*qn[n]
        + xn[iX]*qn[r]*qn[n]) / MASS);
      pi1 = PiMmntm(0.5*qn[r], qn[n], xn[iX]);
      pi2 = PiMmntm(qn[r], 0.5*qn[n], xn[iX]);
      si = spline(pn, pi1, NP, spln_coe1, spln_coe2, spln_coe3, i);
      sm = spline(pn, pi2, NP, spln_coe1, spln_coe2, spln_coe3, m);
      z = fctr_rn*gFactor*sm*si/(pow(pi1*facCntr, lprime)
        * pow(pi2*facCntr, alphapp[_CHNIDX_l_]));
      sum += wxn[iX]*z;
    }
  }

  return sum;
}

template <class CR>
CR kMatrix<CR>::findKInhmgnsDetFunc(REAL energyLoc)
{
  // A function that pack find Det
  genKMat_Inhmgns(energyLoc, avg_mN);
  return MatOpt::findDet<CR>(kMat, dimKMat);
}

// generate tNd_tildeT for Inhomogeneous equation [Eq.(24)]
template <class CR>
void kMatrix<CR>::gen_tNd_tildeT(CR E3, REAL q0, REAL lambda, REAL I) {

  if (! flag_Pis_array_set) {
    show_message_fdv("ERROR (kMatrix::gen_tNd_tildeT): Pi_1/2_array not ready",
      _VBSLVL_LOW_);
    return;
  }
  if (! flag_mNE_array_set) {
    show_message_fdv("ERROR (kMatrix::gen_tNd_tildeT): mNE_qr_array not ready",
      _VBSLVL_LOW_);
    return;
  }
  if (! flag_h2wf_tNd_set) {
    show_message_fdv("ERROR (kMatrix::gen_tNd_tildeT): h2wf_tNd not ready",
      _VBSLVL_LOW_);
    return;
  }
  show_message_timestamp_fdv("-->kMatrix::gen_tNd_tildeT: begins", _VBSLVL_HIGH_);

  size_t N_tNd = dimKMat;

  delete[] tNd_tildeT;
  tNd_tildeT = new CR[N_tNd] {0.0};

  for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha ++) {

    REAL alpha[_N_QMN_CHN_];
    copy_3Nchn(channels[idx_alpha], alpha);

    show_message_timestamp_fdv("-->kMatrix::gen_tNd_tildeT: starts to build two_body tMat",
      _VBSLVL_HIGH_);
    tMatrix<CR>* ptr_tMat = tmat_2body(channels[idx_alpha], NQ, mNE_qr_array,
      NP, PF_pn, NQ*NX, Pi_1_array);
    show_message_timestamp_fdv("-->kMatrix::gen_tNd_tildeT: finished building two_body tMat",
      _VBSLVL_HIGH_);

    size_t nch_2body = ptr_tMat->readNCH();
    REAL lprime[2] = {alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]};

    #pragma omp parallel for
    for (size_t idx_kr = 0; idx_kr < NP*NQ; idx_kr++) {// discrete label other than alpha
      size_t k, r;
      mapSerialIdxToPQ(idx_kr, k, r);
      REAL alpha_dtrn[_N_QMN_CHN_];
      gen_alphaNd_3Nchn(alpha_dtrn, lambda, I, opJ);

      for(size_t idx_lprime = 0; idx_lprime < nch_2body; idx_lprime++){ // Loop of l prime
        size_t idx_alphabar = (idx_lprime == 0 ? idx_alpha : alphabarArray[idx_alpha]);
        for(size_t h = 0; h < NX; h++) {
          size_t idx_rh = serialIdx_by_QX(r, h);
          CR Pi_1 = Pi_1_array[idx_rh];
          CR Pi_2 = Pi_2_array[idx_rh];
          CR z1 = 1.0/pow(Pi_1, lprime[idx_lprime]);

          CR tMatFactor =
            ptr_tMat->readTMat_Mat(alpha[_CHNIDX_l_], lprime[idx_lprime],
            k, idx_rh, r, alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]);

          CR ldpart {0.0};
          for (unsigned ld = 0; ld <= 2; ld += 2) {
            alpha_dtrn[_CHNIDX_l_] = (REAL) ld;
            alpha_dtrn[_CHNIDX_lcpld_] = (ld == 0 ? 2.0 : 0.0);
            ldpart += Gqqx_3Nchn<CR>(channels[idx_alphabar], alpha_dtrn,
              facCntr*qn[r], q0, xn[h]) * h2wf_tNd[ld/2*NX*NQ + idx_rh]
              / (ld == 0 ? 1.0 : pow(Pi_2, ld));
          }

          tNd_tildeT[serialIdx_by_alphaPQ(idx_alpha, k, r)]
            += wxn[h] * tMatFactor * z1 * ldpart;
        }
      }
    }
    delete ptr_tMat;
  }

  show_message_timestamp_fdv("-->kMatrix::gen_tNd_tildeT: ends", _VBSLVL_HIGH_);
}

// generate tNd_omega for Inhomogeneous equation [Eq.(30)]
template <class CR>
void kMatrix<CR>::gen_tNd_omega(CR E3, REAL q0_p, REAL q0, REAL lambda, REAL I) {

  if (! flag_Pis_array_set) {
    show_message_fdv("ERROR (kMatrix::gen_tNd_omega): Pi_1/2_array not ready",
      _VBSLVL_LOW_);
    return;
  }
  if (! flag_mNE_array_set) {
    show_message_fdv("ERROR (kMatrix::gen_tNd_omega): mNE_qr_array not ready",
      _VBSLVL_LOW_);
    return;
  }
  if (! flag_omegas_array_set) {
    show_message_fdv("ERROR (kMatrix::gen_tNd_omega): omega_2_array not ready",
      _VBSLVL_LOW_);
    return;
  }
  if (! flag_h2wf_tNd_set) {
    show_message_fdv("ERROR (kMatrix::gen_tNd_omega): h2wf_tNd not ready",
      _VBSLVL_LOW_);
    return;
  }
  show_message_timestamp_fdv("-->kMatrix::gen_tNd_omega: begins", _VBSLVL_HIGH_);

  size_t N_tNd = NQ*NX*nChannels;
  delete[] tNd_omega;
  tNd_omega = new CR[N_tNd] {0.0};

  for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha ++) {

    REAL alpha[_N_QMN_CHN_];
    copy_3Nchn(channels[idx_alpha], alpha);

    show_message_timestamp_fdv("-->kMatrix::gen_tNd_omega: starts to build two_body tMat",
      _VBSLVL_HIGH_);
    tMatrix<CR>* ptr_tMat = tmat_2body(channels[idx_alpha], NQ, mNE_qr_array,
      NQ*NX, omega_2_array, NQ*NX, Pi_1_array);
    show_message_timestamp_fdv("-->kMatrix::gen_tNd_omega: finished building two_body tMat",
      _VBSLVL_HIGH_);

    size_t nch_2body = ptr_tMat->readNCH();
    REAL lprime[2] = {alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]};

    #pragma omp parallel for
    for (size_t idx_rl = 0; idx_rl < NQ*NX; idx_rl++) {// discrete label other than alpha
      size_t r, l;
      mapSerialIdxToQX(idx_rl, r, l);
      REAL alpha_dtrn[_N_QMN_CHN_];
      gen_alphaNd_3Nchn(alpha_dtrn, lambda, I, opJ);

      for(size_t idx_lprime = 0; idx_lprime < nch_2body; idx_lprime++){ // Loop of l prime
        size_t idx_alphabar = (idx_lprime == 0 ? idx_alpha : alphabarArray[idx_alpha]);
        for(size_t h = 0; h < NX; h++) {
          size_t idx_rh = serialIdx_by_QX(r, h);
          CR Pi_1 = Pi_1_array[idx_rh];
          CR Pi_2 = Pi_2_array[idx_rh];
          CR z1 = 1.0/pow(Pi_1, lprime[idx_lprime]);

          CR tMatFactor =
            ptr_tMat->readTMat_Mat(alpha[_CHNIDX_l_], lprime[idx_lprime],
            idx_rl, idx_rh, r, alpha[_CHNIDX_l_], alpha[_CHNIDX_lcpld_]);

          CR ldpart {0.0};
          for (unsigned ld = 0; ld <= 2; ld += 2) {
            alpha_dtrn[_CHNIDX_l_] = (REAL) ld;
            alpha_dtrn[_CHNIDX_lcpld_] = (ld == 0 ? 2.0 : 0.0);
            ldpart += Gqqx_3Nchn<CR>(channels[idx_alphabar], alpha_dtrn,
              facCntr*qn[r], q0, xn[h]) * h2wf_tNd[ld/2*NX*NQ + idx_rh]
              / (ld == 0 ? 1.0 : pow(Pi_2, ld));
          }

          tNd_omega[serialIdx_by_alphaQX(idx_alpha, r, l)]
            += wxn[h] * tMatFactor * z1 * ldpart;
        }
      }
    }
    delete ptr_tMat;
  }

  show_message_timestamp_fdv("-->kMatrix::gen_tNd_omega: ends", _VBSLVL_HIGH_);

}


// deuteron wf. on phi_ld(Pi_2)
// Prerequisite: genPiArrays()
// return H2 binding energy
template <class CR>
CR kMatrix<CR>::gen_h2wf_tNd()
{

  if (! flag_Pis_array_set) {
    show_message_fdv("ERROR (kMatrix::gen_h2wf_tNd): Pi_2_array not ready",
      _VBSLVL_LOW_);
    return 100.0;
  }
  show_message_fdv("-->kMatrix::gen_h2wf_tNd: begins", _VBSLVL_HIGH_);

  size_t NPi = NQ*NX;
  CR BE2H = gen_deutwf_generic(NPi, Pi_2_array, h2wf_tNd);
  flag_h2wf_tNd_set = true;
  return BE2H;
  show_message_fdv("-->kMatrix::gen_h2wf_tNd: ends", _VBSLVL_HIGH_);
}

// Return H2 binding energy
template <class CR>
CR kMatrix<CR>::gen_deutwf_generic(size_t NPi, CR *Pi, CR *wfArray)
{

  if (! AuxMeshSetFlag) {
    cout << "kMatrix::gen_deutwf_generic: AuxMesh not ready! nothing to be done";
    return 10.0;
  }

  tMatrix<COMPLEX> h2WaveFunc(NPi);
  h2WaveFunc.setFacContour(facCntr);
  h2WaveFunc.setPotName(potName);
  h2WaveFunc.loadPotPara(numPara2B, para_2B);
  h2WaveFunc.getChengduH2wf(NAuxMesh, pn_aux, wpn_aux, Pi);

  for(size_t i = 0; i < NPi; i++){
    wfArray[i] = h2WaveFunc.readH2wfVal(0, i);
    wfArray[NPi + i] = h2WaveFunc.readH2wfVal(2, i);
  }
  return h2WaveFunc.readH2BE();
}

template  <class CR>
tMatrix<CR>* kMatrix<CR>::tmat_2body(REAL *alpha, size_t n_mNE, CR* mNE_lst,
  size_t n_pp, CR* pp_lst, size_t n_p,  CR* p_lst) {

  if (! AuxMeshSetFlag) {
    cout << "kMatrix::tmat_2body: AuxMesh not ready! null pointer returned";
    return NULL;
  }

  int intL, intS, intJ;
  intL = (int) alpha[_CHNIDX_l_]; // L
  intS = (int) alpha[_CHNIDX_s_]; // S
  intJ = (int) alpha[_CHNIDX_j_]; // J (I)

  tMatrix<CR>* ptr_tMat = new tMatrix<CR>(intL, intS, intJ, n_mNE, n_pp, n_p);
  ptr_tMat->setFacContour(facCntr);
  ptr_tMat->copyMEArray(mNE_lst);
  ptr_tMat->setPotName(potName);
  ptr_tMat->loadPotPara(numPara2B, para_2B);
  ptr_tMat->chengduT_generic(NAuxMesh, pn_aux, wpn_aux, pp_lst, p_lst);

  return ptr_tMat;
}

// Deuteron wf values on Pi_2[x_h] [Eqs.(35), (36)]
// Pi_2 = sqrt(q0_p^2 + q0^2/4 + q0_p*q0*x)
template <class CR>
void kMatrix<CR>::gen_PInvG0_phi_initial(REAL q0, REAL q0_p, CR phi_ld []) {

  CR *Pi = new CR[NX];
  for(size_t iX=0; iX<NX; iX++) {
    Pi[iX] = PiMmntm(q0_p, 0.5*q0, xn[iX]);
  }
  gen_deutwf_generic(NX, Pi, phi_ld);
  delete [] Pi;
}

// Deuteron wf values on Pi_1[x_h] [Eqs.(35), (36)]
// Pi_1 = sqrt(q0_p^2/4 + q0^2 + q0_p*q0*x)
template <class CR>
void kMatrix<CR>::gen_PInvG0_phi_final(REAL q0, REAL q0_p, CR phi_ld []) {

  gen_PInvG0_phi_initial(2.0*q0, 0.5*q0_p, phi_ld);
}

// [Eq.(36)]
template <class CR>
CR kMatrix<CR>::get_PInvG0(REAL energyLoc, REAL q0_p, REAL input_lamdap, REAL input_Ip,
  REAL q0, REAL input_lamda, REAL input_I) {

  show_message_timestamp_fdv("-->kMatrix::get_PInvG0: begins", _VBSLVL_HIGH_);

  CR *phi_ld_Pi1 = new CR[2*NX], *phi_ld_Pi2 = new CR[2*NX];
  if (! xmeshFlag) {
    cout << "-->kMatrix::get_PInvG0: xmesh not allocated" << endl;
    return 0.0;
  }

  gen_PInvG0_phi_initial(q0, q0_p, phi_ld_Pi2);
  gen_PInvG0_phi_final(q0, q0_p, phi_ld_Pi1);

  REAL alpha_d[_N_QMN_CHN_], alpha_d_prime[_N_QMN_CHN_];
  gen_alphaNd_3Nchn(alpha_d, input_lamda, input_I, opJ);
  gen_alphaNd_3Nchn(alpha_d_prime, input_lamdap, input_Ip, opJ);

  CR PG {0.0};
  for (int ld_p = 0; ld_p <= 2; ld_p += 2) {
    alpha_d_prime[_CHNIDX_l_] = (REAL) ld_p;
    alpha_d_prime[_CHNIDX_lcpld_] = (ld_p == 0 ? 2.0 : 0.0);

    for (size_t h = 0; h < NX; h++) {
      CR Pi_1 = PiMmntm(0.5*q0_p, q0, xn[h]);
      CR Pi_2 = PiMmntm(q0_p, 0.5*q0, xn[h]);
      CR ldpart {0.0};
      for (int ld = 0; ld <= 2; ld += 2) {
        alpha_d[_CHNIDX_l_] = (REAL) ld;
        alpha_d[_CHNIDX_lcpld_] = (ld == 0 ? 2.0 : 0.0);
        ldpart += phi_ld_Pi2[ld/2*NX + h]/pow(Pi_2, ld)
          *Gqqx_3Nchn<CR>(alpha_d_prime, alpha_d, q0_p, q0, xn[h]);
      }
      PG += wxn[h]*(energyLoc - Pi_1*Pi_1/avg_mN - 0.75*q0_p*q0_p/avg_mN)*
        phi_ld_Pi1[ld_p/2*NX + h]/pow(Pi_1, ld_p)*ldpart;
    }
  }
  delete [] phi_ld_Pi1;
  delete [] phi_ld_Pi2;

  show_message_timestamp_fdv("-->kMatrix::get_PInvG0: ends", _VBSLVL_HIGH_);
  return PG;
}

// [Eq.(37)]
template <class CR>
CR kMatrix<CR>::get_PT(REAL energyLoc, REAL q0_p, REAL input_lamdap, REAL input_Ip, REAL q0, REAL input_lamda, REAL input_I){

  show_message_timestamp_fdv("-->kMatrix::get_PT: begins", _VBSLVL_HIGH_);

  REAL alpha[_N_QMN_CHN_], alpha_d_prime[_N_QMN_CHN_];
  gen_alphaNd_3Nchn(alpha_d_prime, input_lamdap, input_Ip, opJ);
  CR PT {0.0};

  // phi_ld(omega_1)
  CR *phi_ld = new CR[2*NQ*NX];
  gen_deutwf_generic(NQ*NX, omega_1_array, phi_ld);

  for (int ld_p = 0; ld_p <= 2; ld_p += 2) {
    alpha_d_prime[_CHNIDX_l_] = (REAL) ld_p;
    alpha_d_prime[_CHNIDX_lcpld_] = (ld_p == 0 ? 2.0 : 0.0);

     CR rpart {0.0};
      for(size_t r=0; r<NQ; r++) {
        CR hpart {0.0};
        for(size_t h=0; h<NX; h++) {
          CR apart {0.0};
          for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha++) {
            copy_3Nchn(channels[idx_alpha], alpha);
            apart += TAMP_omega[serialIdx_by_alphaQX(idx_alpha, r, h)]
              /pow(PiMmntm(q0_p, 0.5*facCntr*qn[r], xn[h]), alpha[_CHNIDX_l_])
              *Gqqx_3Nchn<CR>(alpha_d_prime, alpha, q0_p, facCntr*qn[r], xn[h]);
          }
          hpart += wxn[h]*phi_ld[ld_p/2*NQ*NX + r*NX + h]
            /pow(PiMmntm(0.5*q0_p, facCntr*qn[r], xn[h]), ld_p)*apart;
        }
        rpart += wqn[r]*qn[r]*qn[r]*hpart;
      }
      PT += pow(facCntr, 3)*rpart;
  }
  delete [] phi_ld;

  show_message_timestamp_fdv("-->kMatrix::get_PT: ends", _VBSLVL_HIGH_);
  return PT;
}

// template <class CR>
// CR kMatrix<CR>::get_U(REAL energyLoc, REAL q0_p, REAL input_lamdap, REAL input_Ip, REAL q0, REAL input_lamda, REAL input_I){

//   facLabel();

//   genPiArrays(q0);
//   gen_omega_arrays(q0_p);
//   gen_h2wf_tNd();

//   genIntzTable_Inhmgns(energyLoc);

//   CR PG = get_PInvG0(energyLoc, q0_p, input_lamdap, input_Ip, q0, input_lamda, input_I);

//   // debug
//   cout << "--->PG = " << PG << endl;

//   genKMat_Inhmgns(energyLoc, avg_mN);
//   gent_Nd(energyLoc, avg_mN, q0, input_lamda, input_I);
//   gen_tilde_T(energyLoc, q0, input_lamda, input_I);

//   gent_Nd(energyLoc, avg_mN, q0, input_lamda, input_I, true, q0_p);

//   genTAMP_omega(energyLoc, avg_mN, q0_p, q0, input_lamda, input_I);

//   CR PT = get_PT(energyLoc, q0_p, input_lamdap, input_Ip, q0, input_lamda, input_I);

//   // debug
//   cout << "--->PT = " << PT << endl;

//   return PG + PT;
// }


// Calculate U_{lambda' I', lambda I} [on the l.h.s. of Eq.(33)]
template <class CR>
CR kMatrix<CR>::get_U_elastic(REAL q0, REAL lambda_prime, REAL I_prime, REAL lambda, REAL I){

  stringstream ss;
  ss << "-->kMatrix::get_U_elastic: (lambda_p, I_p; lambda, I) = ("
    << lambda_prime << ", " << I_prime << "; " << lambda << ", " << I << ")";
  show_message_fdv(ss.str(), _VBSLVL_HIGH_);

  genPiArrays(q0);
  gen_omega_arrays(q0);
  CR BE2H = gen_h2wf_tNd();

  REAL E3 = BE2H.real() + 0.75*q0*q0 / avg_mN;
  gen_mNE_qr_array(E3);
  genIntzTable_Inhmgns(E3);

  CR PG = get_PInvG0(E3, q0, lambda_prime, I_prime, q0, lambda, I);
  ss.str("");
  ss << "-->kMatrix::get_U_elastic: PinvG0 = " << setprecision(8) << PG;
  show_message_fdv(ss.str(), _VBSLVL_NORMAL_);

  genKMat_Inhmgns(E3, avg_mN);
  gen_tNd_tildeT(E3, q0, lambda, I);
  gen_tilde_T(E3, q0, lambda, I);
  // release kMat to save memory **
  delete[] kMat;
  kMat = NULL;
  delete[] tNd_tildeT;
  tNd_tildeT = NULL;

  gen_tNd_omega(E3, q0, q0, lambda, I);
  genTAMP_omega(E3, avg_mN, q0, q0, lambda, I);

  CR PT = get_PT(E3, q0, lambda_prime, I_prime, q0, lambda, I);
  ss.str("");
  ss << "-->kMatrix::get_U_elastic: PT = " << setprecision(8) << PT;
  show_message_fdv(ss.str(), _VBSLVL_NORMAL_);

  return PG + PT;
}



#endif
