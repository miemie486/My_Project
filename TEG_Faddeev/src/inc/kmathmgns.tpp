#ifndef KMATHMGNSTPP
#define KMATHMGNSTPP

// Equation numbers refer to codenote [January 21, 2021]

/* Fill out big table intzTable_Hmgns(...) for all values of alphabar, alpha,
qr, qn, pi and pm [Eq.(9)]
Prerequisite : genGTable()
*/
template <class CR>
void kMatrix<CR>::genIntzTable_Hmgns() {

  if (! GTableAllocFlag) {
    cout << "-->kMatrix::genIntzTable_Hmgns : Error GqqxTable not yet constructed " << endl;
    return;
  }
  show_message_timestamp_fdv("-->kMatrix::genIntzTable_Hmgns: begins ",
    _VBSLVL_NORMAL_);

  show_message_fdv("-->kMatrix::genIntzTable_Hmgns -> about to alloc MEM (GB) "
    + to_string(double(dimKMat*dimKMat*sizeof(CR))/(_SIZE_GB_)),
    _VBSLVL_HIGH_);
  delete[] intzTable_Hmgns;
  intzTable_Hmgns = new CR[dimKMat*dimKMat];

  #pragma omp parallel for
  for (size_t idx_ir = 0; idx_ir < NP*NQ; idx_ir++) {
    for (size_t idx_alpha = 0; idx_alpha < nChannels; idx_alpha++) {
      for (size_t idx_beta = 0; idx_beta < nChannels; idx_beta++) {
        size_t i, r;
        mapSerialIdxToPQ(idx_ir, i, r);
        for (size_t m = 0; m < NP; m++) {
          for (size_t n = 0; n < NQ; n++) {
            intzTable_Hmgns[serialIdxKMat(idx_alpha, i, r, idx_beta, m, n)]
              = intz_Hmgns(channels[idx_alpha][_CHNIDX_l_], channels[idx_alpha],
              channels[idx_beta], GqqxTable[idx_alpha][idx_beta], i, m, r, n);
          }
        }
      }
    }
  }

  show_message_timestamp_fdv("-->kMatrix::genIntzTable_Hmgns: ends ",
    _VBSLVL_NORMAL_);
}

/*
Given 3-body CM energy E3CM, fill out kMat_Hmgns and return its determinant
[Eq.(11)]
*/
template <class CR>
CR kMatrix<CR>::findKHmgnsDetFunc(CR E3CM){

  gen_mNE_qr_array(E3CM);
  genKMat_Hmgns(E3CM, avg_mN);
  return MatOpt::findDet<CR>(kMat, dimKMat);
}

/*
x integration needed for homogeneous-equation kMat elements [Eq.(9)]
*/
template <class CR>
CR kMatrix<CR>::intz_Hmgns(REAL lprime, REAL alpha[], REAL alphapp[], CR *GArray,
  size_t& i, size_t& m, size_t& r, size_t& n) {

  CR z {0.0}, sum {0.0}, gFactor {0.0};
  size_t index_GValue;
  REAL si {0.0}, sm {0.0};
  REAL pi1, pi2;

  for(size_t iX=0; iX< NX; iX++){
    mapGArray(index_GValue, r, n, iX, NQ, NX);
    gFactor = GArray[index_GValue];
    if( abs(gFactor) > _INTZ_EPS_ ) {
      // pi1 = sqrt(0.25*qn[r]*qn[r]+qn[n]*qn[n]+qn[r]*qn[n]*xn[iX]);
      // pi2 = sqrt(qn[r]*qn[r]+0.25*qn[n]*qn[n]+qn[r]*qn[n]*xn[iX]);
      pi1 = PiMmntm(0.5*qn[r], qn[n], xn[iX]);
      pi2 = PiMmntm(qn[r], 0.5*qn[n], xn[iX]);
      si = spline(pn, pi1, NP, spln_coe1, spln_coe2, spln_coe3, i);
      sm = spline(pn, pi2, NP, spln_coe1, spln_coe2, spln_coe3, m);
      z = gFactor*sm*si/(pow(pi1*facCntr, lprime) * pow(pi2*facCntr, alphapp[_CHNIDX_l_]));
      sum += wxn[iX]*z;
    }
  }

  return sum;
}

#endif
