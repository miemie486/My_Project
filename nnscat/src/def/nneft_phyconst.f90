! Bingwei 02/03/2011
! - added PC_hbarc_cube
! Bingwei 12/16/2010
! Physical constants
! * The default units are MeV and MeV^(-1), unless noted otherwise.

module nneft_phyconst

  use nneft_type
  implicit none

  ! PC = Physical Constant
  ! PC_default_* : the real-world values of various constant
  ! Essential constants
  integer, parameter  :: PC_default_GAUSSIAN_POWER = 4
  integer, parameter  :: PC_default_GAUSS_CMPLX = 4
  integer             :: GAUSSIAN_POWER = PC_default_GAUSSIAN_POWER
  integer             :: PC_GAUSS_CMPLX = PC_default_GAUSS_CMPLX

  ! #####  Begin MeV or MeV^(-1)  #####
  real(NER), parameter  ::                  &
  & PC_default_hbarc = 197.3286_NER,        &
!   & PC_default_gA = 1.26_NER,              &
!   & PC_default_gA = 1.289_NER,              &
   & PC_default_gA = 1.29_NER,              &
!   & PC_default_mN = 938.27_NER ,            &
  & PC_default_mN = 939.0_NER ,             &
!   & PC_default_fpi = 92.2_NER ,             &
  & PC_default_fpi = 92.4_NER ,             &
  & PC_default_mpi = 138.0_NER ,            &
  ! & PC_default_mpi = 139.57_NER ,           &
!   & PC_delta_hA = 0.800_NER,                &
!   & PC_delta_hA = 1.000_NER,                &
  & PC_delta_hA = 1.400_NER,                &
!   & PC_delta  =  293.0_NER,                 &
  & PC_delta  =  293.73_NER,                 &
!   & PC_delta  =  320.0_NER,                 &
  & PC_default_deltacutoff = 1200.0_NER

  real(NER)  :: PC_delta_cutoff = PC_default_deltacutoff


  ! Values of ci's
  real(NER), parameter  ::                  &
  ! Epelbaum et al, 1999
  & PC_Epel99_c1 = -0.81E-3_NER,            &
  & PC_Epel99_c2 = 3.28E-3_NER,             &
  & PC_Epel99_c3 = -4.7E-3_NER,             &
  & PC_Epel99_c4 = 3.4E-3_NER,              &
  ! Krebs' Q^2 delta-less fit
  & PC_Kreb_c1 = -0.57E-3_NER,              &
  ! & PC_Kreb_c2 = 2.84E-3_NER,               &
  & PC_Kreb_c3 = -3.87E-3_NER,              &
  & PC_Kreb_c4 = 2.89E-3_NER,               &
  ! Valderrama
  & PC_Pavon_c1 = -0.81E-3_NER,             &
  ! & PC_Pavon_c2 = 3.28E-3_NER,              &
  & PC_Pavon_c3 = -3.4E-3_NER,              &
  & PC_Pavon_c4 = 3.4E-3_NER,               &
  ! Roy Steiner equation, Table II of PRC 96, 024004
  & PC_RSE2_c1 = -0.74E-3_NER,              &
  ! & PC_RSE2_c2 = 0.0_NER,                   &
  & PC_RSE2_c3 = -3.61E-3_NER,              &
  & PC_RSE2_c4 = 2.44E-3_NER,               &
  & PC_RSE3_c1 = -1.07E-3_NER,              &
  ! & PC_RSE3_c2 = 3.20E-3_NER,               &
  & PC_RSE3_c3 = -5.32E-3_NER,              &
  & PC_RSE3_c4 = 3.56E-3_NER,               &
  & PC_RSE4_c1 = -1.10E-3_NER,              &
  ! & PC_RSE4_c2 = 3.57E-3_NER,               &
  & PC_RSE4_c3 = -5.54E-3_NER,              &
  & PC_RSE4_c4 = 4.17E-3_NER,               &
  & PC_default_c1 = PC_RSE2_c1,             &
  ! & PC_default_c2 = PC_RSE2_c2,             &
  & PC_default_c3 = PC_RSE2_c3,             &
  & PC_default_c4 = PC_RSE2_c4,             &
  ! add those parameters followed at oct 15 2018 by pray to calculate deltafull potential
  ! which come from "reconciling threshold and subthreshold expansions for pion-nucleon scattering"
  ! arXiv:1610.08978v2 [nucl-th] 27 Apr 2017
  ! 1-3 for the NLO sigma_2  4-6 for N2LO sigma_3  7-9 for N3LOsigma_4
  & PC_DELTAF1_c1 = -0.74E-3_NER,              &
  & PC_DELTAF1_c2 = -0.49E-3_NER,              &
  & PC_DELTAF1_c3 = -0.65E-3_NER,              &
  & PC_DELTAF1_c4 =  0.96E-3_NER,              &
  & PC_DELTAF2_c1 = -0.69E-3_NER,              &
  & PC_DELTAF2_c2 =  0.81E-3_NER,              &
  & PC_DELTAF2_c3 = -0.44E-3_NER,              &
  & PC_DELTAF2_c4 =  0.64E-3_NER,              &
  & PC_DELTAF3_c1 = -0.69E-3_NER,              &
  & PC_DELTAF3_c2 =  0.40E-3_NER,              &
  & PC_DELTAF3_c3 = -0.49E-3_NER,              &
  & PC_DELTAF3_c4 =  0.64E-3_NER,              &
  & PC_DELTAF4_c1 = -1.25E-3_NER,              &
  & PC_DELTAF4_c2 =  1.37E-3_NER,              &
  & PC_DELTAF4_c3 = -2.41E-3_NER,              &
  & PC_DELTAF4_c4 =  1.66E-3_NER,              &
  & PC_DELTAF5_c1 = -1.24E-3_NER,              &
  & PC_DELTAF5_c2 =  0.79E-3_NER,              &
  & PC_DELTAF5_c3 = -2.49E-3_NER,              &
  & PC_DELTAF5_c4 =  1.67E-3_NER,              &
  & PC_DELTAF6_c1 = -1.12E-3_NER,              &
  & PC_DELTAF6_c2 =  1.02E-3_NER,              &
  & PC_DELTAF6_c3 = -2.27E-3_NER,              &
  & PC_DELTAF6_c4 =  1.21E-3_NER,              &
  & PC_DELTAF7_c1 = -1.11E-3_NER,              &
  & PC_DELTAF7_c2 =  1.52E-3_NER,              &
  & PC_DELTAF7_c3 = -1.99E-3_NER,              &
  & PC_DELTAF7_c4 =  1.88E-3_NER,              &
  & PC_DELTAF8_c1 = -1.11E-3_NER,              &
  & PC_DELTAF8_c2 =  1.29E-3_NER,              &
  & PC_DELTAF8_c3 = -2.15E-3_NER,              &
  & PC_DELTAF8_c4 =  1.94E-3_NER,              &
  & PC_DELTAF9_c1 = -1.10E-3_NER,              &
  & PC_DELTAF9_c2 =  1.20E-3_NER,              &
  & PC_DELTAF9_c3 = -2.19E-3_NER,              &
  & PC_DELTAF9_c4 =  1.77E-3_NER,              &
  & PC_DELTAF10_c1 = -0.74E-3_NER,              &
  & PC_DELTAF10_c2 =  1.81E-3_NER,              &
  & PC_DELTAF10_c3 = -3.61E-3_NER,              &
  & PC_DELTAF10_c4 =  2.17E-3_NER,              &
  & PC_DELTAF11_c1 = -1.08E-3_NER,              &
  & PC_DELTAF11_c2 =  3.26E-3_NER,              &
  & PC_DELTAF11_c3 = -5.39E-3_NER,              &
  & PC_DELTAF11_c4 =  3.62E-3_NER,              &
  & PC_DELTAF12_c1 = -1.11E-3_NER,              &
  & PC_DELTAF12_c2 =  3.17E-3_NER,              &
  & PC_DELTAF12_c3 = -5.67E-3_NER,              &
  & PC_DELTAF12_c4 =  4.35E-3_NER

  ! Delta-isobar related parameters
  real(NER), parameter  ::                  &
  & PC_default_c1_del = -0.00057_NER ,      &
  & PC_default_c2_del = -0.00192_NER ,      &
  & PC_default_c3_del = -0.00088_NER ,      &
  & PC_default_c4_del =  0.0005_NER ,       &
  & PC_default_hA = 1.75_NER,               &
  & PC_default_b3 = 0.0_NER,                &
  & PC_default_b8 = 0.0_NER,                &
  & PC_default_cr = 1.0_NER,                &
  & PC_default_delta = 293.0_NER,           &
  & PC_default_delta_Ha = 1.400_NER
  ! Misc
  real(NER), parameter  ::                  &
  & PC_default_Lambda_inner = 1000.0_NER
  ! #####  End MeV or MeV^(-1)  #####

  ! ! ######  Begin GeV or GeV^(-1)  ######
  ! real(NER), parameter  ::                  &
  ! & PC_default_gA = 1.29_NER,               &
  ! & PC_default_mN = 939.0E-03_NER ,         &
  ! & PC_default_fpi = 92.4E-03_NER ,         &
  ! & PC_default_mpi = 138.0E-03_NER ,        &
  ! & PC_default_hbarc = 197.3286E-03_NER

  ! ! Values of ci's
  ! real(NER), parameter  ::                  &
  ! & PC_Pavon_c1 = -0.81_NER,             &
  ! & PC_Pavon_c3 = -3.4_NER,              &
  ! & PC_Pavon_c4 = 3.4_NER,               &
  ! ! Roy Steiner equation, Table II of PRC 96, 024004
  ! & PC_RSE2_c1 = -0.74_NER,                 &
  ! & PC_RSE2_c3 = -3.61_NER,                 &
  ! & PC_RSE2_c4 = 2.44_NER,                  &
  ! & PC_RSE3_c1 = -1.07_NER,                 &
  ! & PC_RSE3_c3 = -5.32_NER,                 &
  ! & PC_RSE3_c4 = 3.56_NER,                  &
  ! & PC_RSE4_c1 = -1.10_NER,                 &
  ! & PC_RSE4_c3 = -5.54_NER,                 &
  ! & PC_RSE4_c4 = 4.17_NER,                  &
  ! & PC_default_c1 = PC_RSE2_c1,             &
  ! & PC_default_c3 = PC_RSE2_c3,             &
  ! & PC_default_c4 = PC_RSE2_c4

  ! ! Misc
  ! real(NER), parameter  ::                  &
  ! & PC_default_Lambda_inner = 1.0_NER
  ! ! ######  End GeV or GeV^(-1)  ######

  real(NER), parameter   ::                               &
  & PC_default_gA_sqr = PC_default_gA * PC_default_gA,             &
  & PC_default_mpi_sqr = PC_default_mpi * PC_default_mpi,          &
  & PC_default_fpi_sqr = PC_default_fpi * PC_default_fpi,          &
  & PC_default_gA_4th = PC_default_gA_sqr * PC_default_gA_sqr,     &
  & PC_default_mpi_4th = PC_default_mpi_sqr * PC_default_mpi_sqr,  &
  & PC_default_fpi_4th = PC_default_fpi_sqr * PC_default_fpi_sqr,  &
  & PC_default_hbarc_cube = PC_default_hbarc * PC_default_hbarc * PC_default_hbarc

  ! Variables for physical constants
  real(NER) ::                        &
  & PC_gA = PC_default_gA,            &
  & PC_mN = PC_default_mN,            &
  & PC_fpi = PC_default_fpi,          &
  & PC_mpi = PC_default_mpi,          &
  & PC_hbarc = PC_default_hbarc,      &
  ! & PC_c1_del = PC_default_c1_del,    &
  ! & PC_c2_del = PC_default_c2_del,    &
  ! & PC_c3_del = PC_default_c3_del,    &
  ! & PC_c4_del = PC_default_c4_del,    &
  & PC_c1 = PC_default_c1,            &
  ! & PC_c2 = PC_default_c2,            &
  & PC_c3 = PC_default_c3,            &
  & PC_c4 = PC_default_c4,            &
  & PC_hA = PC_default_hA,            &
  & PC_hA1 = PC_delta_hA,             &
  ! & PC_b3 = PC_default_b3,            &
  ! & PC_b8 = PC_default_b8,            &
  ! & PC_cr = PC_default_cr,            &
  ! & PC_delta = PC_default_delta,      &
  & PC_Lambda_inner = PC_default_Lambda_inner

  ! real(NER) ::                    &
  ! & PC_gA_sqr = PC_default_gA_sqr,   &
  ! & PC_mpi_sqr = PC_default_mpi_sqr, &
  ! & PC_fpi_sqr = PC_default_fpi_sqr, &
  ! & PC_gA_4th = PC_default_gA_4th,   &
  ! & PC_mpi_4th = PC_default_mpi_4th, &
  ! & PC_fpi_4th = PC_default_fpi_4th, &
  ! & PC_hbarc_cube = PC_default_hbarc_cube

  type :: struct_pc

    real(NER) ::                    &
    & gA = PC_default_gA,           &
    & mN = PC_default_mN,           &
    & fpi = PC_default_fpi,         &
    & mpi = PC_default_mpi,         &
    & hbarc = PC_default_hbarc,     &
    ! & c1_del = PC_default_c1_del,   &
    ! & c2_del = PC_default_c2_del,   &
    ! & c3_del = PC_default_c3_del,   &
    ! & c4_del = PC_default_c4_del,   &
    & c1 = PC_default_c1,           &
    ! & c2 = PC_default_c2,           &
    & c3 = PC_default_c3,           &
    & c4 = PC_default_c4,           &
    ! & hA = PC_default_hA,           &
    ! & b3 = PC_default_b3,           &
    ! & b8 = PC_default_b8,           &
    ! & cr = PC_default_cr,           &
    ! & delta = PC_default_delta,     &
    & Lambda_inner = PC_default_Lambda_inner

  end type struct_pc

  type  :: strct_cis

    real(NER) ::            &
    & c1 = PC_default_c1,   &
    & c3 = PC_default_c3,   &
    & c4 = PC_default_c4

  end type strct_cis

  type  :: Delta_paras

    real(NER)  ::           &
    & c1 = PC_DELTAF1_c1,   &
    & c2 = PC_DELTAF1_c2,   &
    & c3 = PC_DELTAF1_c3,   &
    & c4 = PC_DELTAF1_c4,   &
    & ha = PC_delta_hA,     &
    & del= PC_delta

  end type Delta_paras

end module nneft_phyconst
