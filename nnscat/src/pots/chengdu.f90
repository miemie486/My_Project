! Long 01/01/2021 (V0.41)
! - Added TPE-less version
! Rui 23/11/2020 (V0.40)
! add cutoff 900/1300/1400/6400 for MMWLY
! notice:: parameters of 3s1 and 1s0 at 6400 fitted with Nmesh=200
! Rui 16/11/2020 (V0.39)
! add cutoff 500/1800/2200/4000/4800 for MMWLY
! and modified 1200/1600/2000/2400 N2LO parameters in 1s0
! notice:: parameters of 3s1 and 1s0 at 4000/4800 fitted with Nmesh=200
! Rui 28/10/2020 (V0.38)
! modified 1s0 and 3s1 parameters which fitted by scattering length and binding energy
! Rui 09/10/2020
! add cutoff 1000/1300/1500/2500/3000 in 3s1 and 1s0 at LO
! Bingwei 09/08/2020 (V0.37)
! - Rui added cutoff 3600/4000/4400/4800 in 3s1 and 1s0 at LO
! Rui 04/09/2020
! add cutoff 1200/2000/2400/2800 for MMWLY
! Bingwei Long 08/09/2020 (V0.36)
! - Added diagnostic pot
! Bingwei Long 08/09/2020 (V0.35)
! - Fixed a bug in get_chengdu_MMWLY that caused to fail to standardize 3D1 to
!   3S1 for MMWLY_N2LO
! Rui Peng 07/30/2020 (V0.34)
! - N2LO MMWLY paras fitted by fit_mmwly_smpl
! Rui Peng 07/05/2020 (V0.33)
! - Added MWPLL_400, tested with Faddeev code
! Bingwei Long 06/20/2020 (V0.32)
! - For coupled channels, making sure L = J+1 outputs the same 2X2 matrix
!   as L = J-1
! Bingwei Long  05/18/2020 (V0.31)
! - Energy-dependent complex potentials
! Bingwei Long  03/26/2020 (V0.30)
! - Minor tweaks to variable naming scheme and comments
! Bingwei Long  01/12/2020 (V0.291)
! - Changed pionless_1000 paras so that 3S1 scat length = 5.42fm
! Bingwei Long  01/12/2020 (V0.29)
! - Added chengdu_hotpot_mixed, tested
! Bingwei Long  01/04/2020 (V0.28)
! - Added complex version of MMWLY routines
! Bingwei Long  01/04/2020 (V0.27)
! - Added MMWLY_350
! Rui Peng  01/04/2020 (V0.26)
! - Updated MMWLY NLO 1S0 short-range paras
! Rui Peng  01/02/2020 (V0.25)
! - MMWLY 1S0 fitted to scattering length (LO) and kcm = 53.075 Mev and
!   phs = 62.9169
! Rui Peng  12/24/2019 (V0.24)
! - MMWLY 1S0 fitted to scattering length
! Rui Peng  12/24/2019 (V0.232)
! - Factorized DLSPR and MMWLY routines
! Bingwei Long 12/03/2019 (V0.231)
! - Fixed a bug in MMWLY
! Bingwei Long 11/30/2019
! - Added MMWLY LO
! Bingwei Long 10/27/2019
! - Added generic interface: chengdu_dispatch, chengdu_hotpot, chengdu_hotpot_Edep
! Bingwei Long 10/22/2019
! - Added chengdu_LBWDIB_600
! Bingwei Long 10/21/2019
! - Fitted pionless_1000 3S1 to a = 5.41 fm (added as comments)
! Bingwei Long 10/20/2019
! - Fitted pionless_1000 3S1 to deuteron BE
! - Modified chengdu.README, show_mtrx.f90
! Bingwei Long 10/19/2019
! - Fixed bugs in chengdu.Makefile, chengdu_DLSPR_350 and chengdu_DLSPR_400
! Rui Peng 10/17/2019
! - Added Pionless LO for 1000 MeV
! Rui Peng 10/07/2019
! - Added 350 and 400 MeV (3S1 fitted to B_deuteron)
! Bingwei Long  10/04/2019
!   - Added init_PCs4chengdu
! Bingwei Long  07/02/2019
!   - Added 800, 1600, 3200 MeV LO values
! Bingwei Long  03/18/2019
!   - Improved comments
!   - Modified to accommodate changes made to interface of *_epwrap()
! Bingwei Long  01/28/2019

! Feautures: Deltaless
!
! All units are in MeV or (MeV)^{-1}
! %%%%%%%%%%%%
! Inputs
! %%%%%%%%%%%%
! uptoQn: number of order
! L, S, J: quantum # of partial-wave, as in ^{2S+1}L_J
!   Non-vanishing channel at LO: 1S0(LSJ = 000), 3S1 - 3D1(LSJ = 011)
!     , 3P0(LSJ = 110)
!   For coupled channels, L = lower one of two coupled angular momenta
! po, pi: magnitude of outgoing and incoming c.m. momenta
!
! %%%%%%%%%%%%
! Outputs
! %%%%%%%%%%%%
! potval is the value of partial-wave potentials
! uptoQn = 0 ==> VLO
! uptoQn = 1 ==> VNLO (incremental value)
! Storage of potval:
!   2-by-2 matrix
!   For uncoupled channels, e.g., 1S0(LSJ = 000) or 3P0(LSJ = 110)
!     potval(1, 1) = <L, po | V | L, pi > and (1, 2), (2, 1), and (2, 2) are
!     zero
!   For coupled channels, e.g., 3S1 - 3D1 (LSJ = 011)
!     potval(1, 1) = < L, po   | V | L, pi >
!     potval(1, 2) = < L, po   | V | L+2, pi >
!     potval(2, 1) = < L+2, po | V | L, pi >
!     potval(2, 2) = < L+2, po | V | L+2, pi >
!
! Remarks:
! * For coupled channels, 2X2 matrix of V is the same regardless of L (J-1 or J+1)
! * potval(2, 1) = < L+2, po | V | L, pi > = < L, pi   | V | L+2, po >
! * Normalization of V can be explained as follows.
! If V were so weak as to validate the Born approximation, the relation of
! V in uncoupled channels and the respective scattering length 'a' would have
! been
! < L, p | V | L, p> = 2/pi * a + ... for p -> 0

module chengdu

  use potwrap
  use cmplx_epmod
  use nopion_epmod
  use testwrap
  use vbar_wfs
  use vbar_3p0_wfs
  implicit none

  private
  public  :: gnrc_chengdu_pot_real, gnrc_chengdu_pot_cmplx,                    &
    & chengdu_hotpot, chengdu_hotpot_basic, chengdu_hotpot_mixed,              &
    & chengdu_hotpot_is_Edep, chengdu_hotpot_Edep, chengdu_dispatch,           &
    & chengdu_cmplx_hotpot_mixed, chengdu_cmplx_hotpot, chengdu_cmplx_hotpot_Edep, &
    & chengdu_pionless_1000, read_chengdu_version,                             &
    & chengdu_DLSPR_350, chengdu_DLSPR_400, chengdu_DLSPR,                     &
    & chengdu_DLSPR_800, chengdu_DLSPR_1600, chengdu_DLSPR_3200,               &
    & chengdu_LBWDIB_600, chengdu_MMWLY_350, chengdu_MMWLY_400,                &
    & chengdu_MMWLY_600, chengdu_MMWLY_800, chengdu_MMWLY_1600,                &
    & chengdu_MMWLY_3200,                                                      &
    & chengdu_MWPLL_400,                                                       &
    & LO_chengdu_epwrap, NLO_chengdu_epwrap, N2LO_chengdu_epwrap

  character(len = 30)  :: verstr = "Version 0.41"

  integer, parameter    :: CHENGDU_REGTYPE = REGTYPE_GAUSSIAN, NCH = 2
  character(len = 256)  :: hotpot_name = "chengdu_MMWLY_600"

  interface chengdu_hotpot
    module procedure chengdu_hotpot_basic, chengdu_hotpot_mixed
  end interface chengdu_hotpot

  abstract interface
    subroutine gnrc_chengdu_pot_real(L, S, J, uptoQn, po, pi, potval)
      import                  :: NER, NCH
      implicit none
      integer, intent(in)     :: L, S, J, uptoQn
      real(NER), intent(in)   :: po, pi
      real(NER), intent(out)  :: potval(1:NCH, 1:NCH)
    end subroutine gnrc_chengdu_pot_real
    subroutine gnrc_chengdu_pot_cmplx(L, S, J, uptoQn, po, pi, potval)
      import                  :: NEC, NCH
      implicit none
      integer, intent(in)     :: L, S, J, uptoQn
      complex(NEC), intent(in)   :: po, pi
      complex(NEC), intent(out)  :: potval(1:NCH, 1:NCH)
    end subroutine gnrc_chengdu_pot_cmplx
  end interface

  real(NER)             :: dlspr1s0(1:14), mmwly1s0(1:6), Cs3s1(1:7),   &
    & Cs3p0(1:5), Cs1p1(1:4), Cs3p1(1:4), Cs3p2(1:5), Cs1d2(1:4),       &
    & Cs3d2(1:4), CS3d3(1:4)

  real(NER), parameter ::  &
    & spr1s0_350_set2_para(1:14) = (/                                       &
    &   0.42647716568097287E+002_NER,   -0.15817891988985733E+000_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_350_set2_para(1:7) =  (/                                        &
    &  -0.54973804635129386E-002_NER,    -0.17006900114863951E-002_NER,     &
    &  -0.29022475006695189E-007_NER,    0.15748583965631266E-008_NER,      &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_350_set2_para(1:5) =  (/                                        &
    & 0.28737080761182509E-007_NER,      0.22221189688459573E-007_NER,      &
    & -0.10425692732865342E-013_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_350_set2_para(1:5) =  (/                                        &
    & -0.16800861064631123E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_350_set2_para(1:4) =  (/                                        &
    & 0.74536008661367232E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_350_set2_para(1:4) =  (/                                        &
    & 0.90918701761135266E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_350_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_350_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_350_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & spr1s0_400_set2_para(1:14) = (/                                       &
    &   0.39275245500005816E+002_NER,   -0.13862263192892474E+000_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_400_set2_para(1:7) =  (/                                        &
    &  -0.45148506612549544E-002_NER,    -0.21129650647474007E-002_NER,     &
    &  -0.30723565621962921E-007_NER,    0.74587285877961200E-009_NER,      &
    &   0.0_NER, 0.0_NER, 0.0_NER/),                                        &
    &                                                                       &
    & Cs3p0_400_set2_para(1:5) =  (/                                        &
    & 0.30164174135677295E-007_NER,       0.22401580364683590E-007_NER,     &
    & -0.26593399259140350E-013_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_400_set2_para(1:5) =  (/                                        &
    & -0.10344282401234271E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_400_set2_para(1:4) =  (/                                        &
    & 0.22521280536709626E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_400_set2_para(1:4) =  (/                                        &
    & 0.65967556790677961E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_400_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_400_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_400_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter  ::  &
    & spr1s0_600_set2_para(1:14) = (/                                       &
    &   0.28476583881316405E+002_NER,   -0.94720748744753225E-001_NER,      &
    &   0.20650483002051136E+002_NER,   -0.97698777628877326E-001_NER,      &
    &   0.95630015528180693E-003_NER,   -0.11694692770857738E+002_NER,      &
    &   0.86955645956220531E-001_NER,   -0.27572312790623797E-002_NER,      &
    &   0.51272513824358814E-007_NER,   -0.15365873931283828E+002_NER,      &
    &   0.24216603991719265E+000_NER,   -0.25138848487176349E-002_NER,      &
    &   0.55453661564745701E-007_NER,    0.61631676811145938E-012_NER/),    &
    &                                                                       &
    & Cs3s1_600_set2_para(1:7) =  (/                                        &
    &  -0.14442429421521084E-002_NER,     -0.37723180431405998E-002_NER,     &
    &  -0.37905375782350260E-007_NER,    -0.36289549203893510E-008_NER,     &
    &   0.0_NER, 0.0_NER, 0.0_NER/),                                        &
    &                                                                       &
    & Cs3p0_600_set2_para(1:5) =  (/                                        &
    & 0.47638609114073750E-007_NER, -0.67617335664002759E-007_NER,          &
    &-0.18955801253529542E-012_NER, 0.0_NER, 0.0_NER/),                     &
    &                                                                       &
    & Cs3p2_600_set2_para(1:5) =  (/                                        &
    & -0.50060906511174551E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_600_set2_para(1:4) =  (/                                        &
    & 0.33105259048813496E-011_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_600_set2_para(1:4) =  (/                                        &
    & 0.42675011907583266E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_600_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_600_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_600_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & spr1s0_800_set2_para(1:14) = (/                                       &
    &   0.26900213333936705E+002_NER,   -0.85152284790577140E-001_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_800_set2_para(1:7) =  (/                                        &
    &  0.29003872118011503E-002_NER,      0.51387312371913661E-002_NER,     &
    & -0.20820565625816905E-007_NER,     -0.10716199065261968E-007_NER,     &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & Cs3p0_800_set2_para(1:5) =  (/                                        &
    & -0.39975686839150306E-007_NER, 0.39434636364644249E-006_NER,          &
    & -0.15159307417491144E-012_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_800_set2_para(1:5) =  (/                                        &
    & -0.38652287513101609E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_800_set2_para(1:4) =  (/                                        &
    & 0.10324794006950635E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_800_set2_para(1:4) =  (/                                        &
    & 0.39761956440587437E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_800_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_800_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_800_set2_para(1:4) =  (/                                        &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & spr1s0_1600_set2_para(1:14) = (/                                      &
    &   0.25741081111240778E+002_NER,   -0.75410791126733698E-001_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_1600_set2_para(1:7) =  (/                                       &
    & -0.46750075253685118E-002_NER,     -0.26429843062470332E+000_NER,      &
    & -0.17488247733838489E-007_NER,     0.76178244824950753E-007_NER,      &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_1600_set2_para(1:5) =  (/                                       &
    & 0.31968616606084914E-009_NER, 0.34410486157143277E-007_NER,           &
    & 0.43012680194265140E-014_NER, 0.0_NER, 0.0_NER/),                     &
    &                                                                       &
    & Cs3p2_1600_set2_para(1:5) =  (/                                       &
    & -0.22686718892723393E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_1600_set2_para(1:4) =  (/                                       &
    & 0.32015497600836161E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_1600_set2_para(1:4) =  (/                                       &
    & 0.38760693975285904E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_1600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_1600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_1600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & spr1s0_3200_set2_para(1:14) = (/                                      &
    &   0.25530308969449205E+002_NER,   -0.72062514706775244E-001_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_3200_set2_para(1:7) =  (/                                       &
    & -0.12838596032815119E-001_NER,     -0.41514797159834611E+001_NER,      &
    & -0.28469173742695489E-006_NER,    -0.38151492832344918E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_3200_set2_para(1:5) =  (/                                       &
    & 0.11323845800451927E-007_NER, -0.10981246749302965E-003_NER,          &
    & -0.30553581583820761E-012_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_3200_set2_para(1:5) =  (/                                       &
    & -0.13118398457593397E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_3200_set2_para(1:4) =  (/                                       &
    & 0.43167705688310954E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_3200_set2_para(1:4) =  (/                                       &
    & 0.38715739925675247E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_3200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_3200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_3200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

! fitting LO by scatlength
! ! empirical scatlength = -24.740(20)
! fitting N2LO by phase shifts and scatlength
  real(NER), parameter  ::  &
    & MMWLY1s0_350_para(1:6) = (/                                            &
    & -0.49527443017614902E-002_NER,       -0.11880991849022070E-002_NER,    &
    &  0.23605908822236607E-007_NER,       -0.27697377278486411E-003_NER,    &
    &  0.49356919946805012E-007_NER,       -0.38218511945315221E-012_NER/)

  real(NER), parameter  ::  &
    & MMWLY1s0_400_para(1:6) = (/                                            &
    & -0.46345112093629202E-002_NER,       -0.13812599401879487E-002_NER,    &
    &  0.21086806723513865E-007_NER,       -0.42546287376399490E-003_NER,    &
    &  0.50805943780941848E-007_NER,       -0.24703699288896164E-012_NER/)

  real(NER), parameter  ::  &
    & MMWLY1s0_450_para(1:6) = (/                                            &
    & -0.43891487903421457E-002_NER,       -0.15417924905996663E-002_NER,    &
    &  0.18652151984588734E-007_NER,       -0.59527537315732279E-003_NER,    &
    &  0.52329511860140106E-007_NER,       -0.17297146771694753E-012_NER/)

  real(NER), parameter  ::  &
    & MMWLY1s0_500_para(1:6) = (/                                            &
    & -0.41940501683613270E-002_NER,       -0.16799644919596995E-002_NER,    &
    &  0.16499403193895097E-007_NER,       -0.42783779918269781E-003_NER,    &
    &  0.48512066078090494E-007_NER,       -0.56991625015055438E-013_NER/)

  real(NER), parameter  ::  &
    & MMWLY1s0_600_para(1:6) = (/                                            &
    & -0.39028006210685457E-002_NER,       -0.19113791264257359E-002_NER,    &
    &  0.13071651933734054E-007_NER,       -0.10699886131465934E-002_NER,    &
    &  0.54788515175561186E-007_NER,       -0.56812585457311099E-013_NER/)

  real(NER), parameter  ::  &
    & MMWLY1s0_800_para(1:6) = (/                                            &
    & -0.35387941376915466E-002_NER,       -0.22691865795908422E-002_NER,    &
    &  0.87320913973441971E-008_NER,       -0.13380262130199116E-002_NER,    &
    &  0.54503555720890117E-007_NER,        0.98670327253654256E-014_NER/)

  real(NER), parameter  ::  &
    & MMWLY1s0_1600_para(1:6) = (/                                           &
    & -0.29771946913657479E-002_NER,       -0.31504605152330055E-002_NER,    &
    &  0.29719305488873427E-008_NER,        0.82615385721146739E-002_NER,    &
    &  0.42209907382510694E-007_NER,        0.49717691251764593E-013_NER/)

  real(NER), parameter  ::  &
    & MMWLY1s0_3200_para(1:6) = (/                                           &
    & -0.26722872208465821E-002_NER,       -0.40356982089648340E-002_NER,    &
    &  0.91813647153278955E-009_NER,        0.12707269864860637E+000_NER,    &
    &  0.15885530867285035E-007_NER,        0.30293767664959401E-013_NER/)

  real(NER), parameter  ::  &
    & MWPLL1s0_400_para(1:9) = (/                                           &
    & -0.53511021078353396E-002_NER,       0.12803625440773314E-007_NER,    &
    & -0.82297860339877651E-003_NER,       0.19279741957230324E-007_NER,    &
    & -0.26036529272854275E-012_NER,      -0.59859675954648383E-003_NER,    &
    &  0.20097824336626228E-007_NER,      -0.78630209320857908E-012_NER,    &
    &  0.21196299898834966E-017_NER/)


!   add 1200/2000/2400/2800 mmwly for each channel
  real(NER), parameter ::  &
    & MMWLY1s0_1200_para(1:6) = (/                                          &
    & -0.31686486633502316E-002_NER,       -0.27797190987699032E-002_NER,    &
    &  0.47143394436464893E-008_NER,       -0.12703874994129903E-004_NER,    &
    &  0.51570029181387921E-007_NER,       0.39876817667954619E-013_NER/), &
    &                                                                       &
    & Cs3s1_1200_set2_para(1:7) =  (/                                       &
    & -0.12668825908618342E-001_NER,     -0.11638980324086622E+000_NER,      &
    & -0.16041302940490850E-006_NER,    -0.21455137390411132E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_1200_set2_para(1:5) =  (/                                       &
    & -0.15607073425249616E-008_NER, 0.34000981167047608E-007_NER,          &
    & 0.46861907507408023E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_1200_set2_para(1:5) =  (/                                       &
    & -0.28347593612560600E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_1200_set2_para(1:4) =  (/                                       &
    & 0.24528518288591582E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_1200_set2_para(1:4) =  (/                                       &
    & 0.38881833048003478E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_1200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_1200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_1200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_2000_para(1:6) = (/                                          &
    & -0.28585563548057182E-002_NER,       -0.34398508249289133E-002_NER,    &
    &  0.20545654031144380E-008_NER,        0.22969786306947099E-001_NER,    &
    &  0.35322089529800873E-007_NER,        0.44914594573730232E-013_NER/), &
    &                                                                       &
    & Cs3s1_2000_set2_para(1:7) =  (/                                       &
    & -0.21944206584707093E-002_NER,     -0.22132248752054690E+000_NER,      &
    & -0.57952346722395748E-007_NER,    0.16911932370415890E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_2000_set2_para(1:5) =  (/                                       &
    & 0.66001753049767985E-009_NER, 0.34982399750526136E-007_NER,          &
    & 0.39960393941440373E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_2000_set2_para(1:5) =  (/                                       &
    & -0.19018684915151714E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_2000_set2_para(1:4) =  (/                                       &
    & 0.36501645137517515E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_2000_set2_para(1:4) =  (/                                       &
    & 0.38731283232071887E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_2000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_2000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_2000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_2400_para(1:6) = (/                                          &
    & -0.27772405429084906E-002_NER,       -0.36746014859360466E-002_NER,    &
    &  0.15096483095713912E-008_NER,        0.48118369074706902E-001_NER,    &
    &  0.27430148434906396E-007_NER,        0.40352022876378480E-013_NER/), &
    &                                                                       &
    & Cs3s1_2400_set2_para(1:7) =  (/                                       &
    & 0.62485155018945579E-003_NER,     -0.57704915109189869E+000_NER,      &
    & -0.31820898927098933E-007_NER,    0.28353391897097949E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_2400_set2_para(1:5) =  (/                                       &
    & 0.78982690154086563E-009_NER, 0.36063565542175693E-007_NER,          &
    & 0.49520900530735847E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_2400_set2_para(1:5) =  (/                                       &
    & -0.16453613477427865E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_2400_set2_para(1:4) =  (/                                       &
    & 0.39475303126497257E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_2400_set2_para(1:4) =  (/                                       &
    & 0.38721524985168739E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_2400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_2400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_2400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_2800_para(1:6) = (/                                          &
    & -0.27177795837334149E-002_NER,       -0.38699626349011809E-002_NER,   &
    &  0.11583722015198334E-008_NER,        0.81780108982680622E-001_NER,   &
    &  0.21595014997734121E-007_NER,        0.34781718707060716E-013_NER/), &
    &                                                                       &
    & Cs3s1_2800_set2_para(1:7) =  (/                                       &
    & 0.13118803810359370E-001_NER,     -0.49181247813847295E+003_NER,       &
    & -0.61618035689136642E-005_NER,     0.75291419082984288E-005_NER,       &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & Cs3p0_2800_set2_para(1:5) =  (/                                       &
    & 0.10708384950704921E-008_NER, 0.11086086197938928E-006_NER,          &
    & 0.13406769002348014E-013_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_2800_set2_para(1:5) =  (/                                       &
    & -0.14564568023302546E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_2800_set2_para(1:4) =  (/                                       &
    & 0.41588735858637634E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_2800_set2_para(1:4) =  (/                                       &
    & 0.38717575165736458E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_2800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_2800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_2800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_3600_para(1:6) = (/                                          &
    & -0.26362959663116291E-002_NER,       0.0_NER,    &
    &  0.0_NER,        0.0_NER,    &
    &  0.0_NER,        0.0_NER/), &
    &                                                                       &
    & Cs3s1_3600_set2_para(1:7) =  (/                                       &
    &  -0.60711E-002_NER,     0.0_NER,      &
    & 0.0_NER,    0.0_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_3600_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER,          &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_3600_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_3600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_3600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_3600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_3600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_3600_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_4000_para(1:6) = (/                                          &
    & -0.26070755680917880E-002_NER,       -0.43029463919601009E-002_NER,    &
    &  0.61902751186230983E-009_NER,        0.25615052136714667E+000_NER,    &
    &  0.53424399259917738E-008_NER,        0.23544131783273070E-013_NER/), &
    &                                                                       &
    & Cs3s1_4000_set2_para(1:7) =  (/                                       &
    & -0.40755330091299634E-002_NER,     0.42453412215091859E+001_NER,      &
    & -0.26322264625611011E-006_NER,    -0.29241243825314195E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_4000_set2_para(1:5) =  (/                                       &
    & -0.68226727804105308E-010_NER, 0.12251631727563841E-007_NER,          &
    & 0.14899010347483228E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_4000_set2_para(1:5) =  (/                                       &
    & -0.11054821533732302E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_4000_set2_para(1:4) =  (/                                       &
    & 0.45369621102558078E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_4000_set2_para(1:4) =  (/                                       &
    & 0.38714271667061442E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_4000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_4000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_4000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_4400_para(1:6) = (/                                          &
    & -0.25828581699516551E-002_NER,       0.0_NER,    &
    &  0.0_NER,        0.0_NER,    &
    &  0.0_NER,        0.0_NER/), &
    &                                                                       &
    & Cs3s1_4400_set2_para(1:7) =  (/                                       &
    & -0.28236E-002_NER,     0.0_NER,      &
    & 0.0_NER,    0.0_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_4400_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER,          &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_4400_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_4400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_4400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_4400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_4400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_4400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_4800_para(1:6) = (/                                          &
    & -0.25624472971879611E-002_NER,       -0.45099223506560072E-002_NER,    &
    &  0.44631438051037783E-009_NER,        0.48007349893574780E+000_NER,    &
    &  -0.11524141107766023E-007_NER,       0.20302751004355259E-013_NER/), &
    &                                                                       &
    & Cs3s1_4800_set2_para(1:7) =  (/                                       &
    & -0.16823198672076066E-002_NER,     0.47063841243881953E+002_NER,      &
    & -0.78233682395786781E-006_NER,    -0.14587045841015296E-005_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_4800_set2_para(1:5) =  (/                                       &
    & 0.62810050992021014E-010_NER, 0.42364241684694717E-007_NER,          &
    & 0.73024210093575117E-015_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_4800_set2_para(1:5) =  (/                                       &
    & -0.96563762246778555E-009_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_4800_set2_para(1:4) =  (/                                       &
    & 0.46832298547720555E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_4800_set2_para(1:4) =  (/                                       &
    & 0.38713769075135848E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_4800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_4800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_4800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_1000_para(1:6) = (/                                          &
    & -0.33181124691830394E-002_NER,       -0.25483306792893298E-002_NER,    &
    &  0.62547444546995307E-008_NER,        -0.12698064376194659E-003_NER,    &
    &  0.49136102168858503E-007_NER,        0.50524822618569185E-013_NER/), &
    &                                                                       &
    & Cs3s1_1000_set2_para(1:7) =  (/                                       &
    & 4.8765294272879282E-002_NER,     0.0_NER,      &
    & 0.0_NER,    0.0_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_1000_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER,          &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_1000_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_1000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_1000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_1000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_1000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_1000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_900_para(1:6) = (/                                          &
    & -0.34165812068935051E-002_NER,       -0.24160479691425134E-002_NER,    &
    &  0.73356910679664988E-008_NER,        -0.40869272925127828E-003_NER,    &
    &  0.49136102168858503E-007_NER,        0.50524822618569185E-013_NER/), &
    &                                                                       &
    & Cs3s1_900_set2_para(1:7) =  (/                                       &
    & 0.87113709647920132E-002_NER,     0.20411114836416330E+000_NER,      &
    & 0.41512542500445431E-007_NER,    -0.56839605771366384E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_900_set2_para(1:5) =  (/                                       &
    & -0.13014271883685450E-007_NER, 0.65928898230536317E-007_NER,          &
    & -0.14956121924349247E-013_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_900_set2_para(1:5) =  (/                                       &
    & -0.35270309883036896E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_900_set2_para(1:4) =  (/                                       &
    & 0.14834499214720965E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_900_set2_para(1:4) =  (/                                       &
    & 0.39326731971430270E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_900_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_900_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_900_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_1300_para(1:6) = (/                                          &
    & -0.31103986808504844E-002_NER,       -0.28822911677480923E-002_NER,    &
    &  0.41537183148985765E-008_NER,        0.22302549819547893E-002_NER,    &
    &  0.46974998073529568E-007_NER,        0.50754331837276121E-013_NER/), &
    &                                                                       &
    & Cs3s1_1300_set2_para(1:7) =  (/                                       &
    & -8.7385234733682038E-003_NER,     0.14229623777691733E+000_NER,      &
    & -0.18763314724310322E-006_NER,    -0.21292199994538430E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_1300_set2_para(1:5) =  (/                                       &
    & -0.71415922789162087E-009_NER, 0.33784552503496068E-007_NER,          &
    & 0.47750851067555052E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_1300_set2_para(1:5) =  (/                                       &
    & -0.26660889759993569E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_1300_set2_para(1:4) =  (/                                       &
    & 0.26826485606564429E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_1300_set2_para(1:4) =  (/                                       &
    & 0.38831443509701096E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_1300_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_1300_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_1300_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_1500_para(1:6) = (/                                          &
    & -0.30160255620977276E-002_NER,       0.0_NER,    &
    &  0.0_NER,        0.0_NER,    &
    &  0.0_NER,        0.0_NER/), &
    &                                                                       &
    & Cs3s1_1500_set2_para(1:7) =  (/                                       &
    & -5.5699456605899055E-003_NER,     0.0_NER,      &
    & 0.0_NER,    0.0_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_1500_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER,          &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_1500_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_1500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_1500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_1500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_1500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_1500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_2500_para(1:6) = (/                                          &
    & -0.27607180465463346E-002_NER,       0.0_NER,    &
    &  0.0_NER,        0.0_NER,    &
    &  0.0_NER,        0.0_NER/), &
    &                                                                       &
    & Cs3s1_2500_set2_para(1:7) =  (/                                       &
    & 1.8042039073015606E-003_NER,     0.0_NER,      &
    & 0.0_NER,    0.0_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_2500_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER,          &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_2500_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_2500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_2500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_2500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_2500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_2500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_3000_para(1:6) = (/                                          &
    & -0.26936204822600824E-002_NER,       0.0_NER,    &
    &  0.0_NER,        0.0_NER,    &
    &  0.0_NER,        0.0_NER/), &
    &                                                                       &
    & Cs3s1_3000_set2_para(1:7) =  (/                                       &
    & -5.6168006619612187E-002_NER,     0.0_NER,      &
    & 0.0_NER,    0.0_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_3000_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER,          &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_3000_set2_para(1:5) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_3000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_3000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_3000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_3000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_3000_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
!     & MMWLY1s0_500_para(1:6) = (/                                          &
!     & -0.26936204822600824E-002_NER,       0.0_NER,    &
!     &  0.0_NER,        0.0_NER,    &
!     &  0.0_NER,        0.0_NER/), &
!     &                                                                       &
    & Cs3s1_500_set2_para(1:7) =  (/                                       &
    & -0.29222864817793378E-002_NER,     -0.30443673269880955E-002_NER,      &
    & -0.35422271311161108E-007_NER,     -0.33049294424165896E-009_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_500_set2_para(1:5) =  (/                                       &
    & 0.34434390126552547E-007_NER, 0.15866399060768903E-007_NER,          &
    & -0.66122887358886656E-013_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_500_set2_para(1:5) =  (/                                       &
    & -0.63403863705893584E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_500_set2_para(1:4) =  (/                                       &
    & -0.11946227863389642E-009_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_500_set2_para(1:4) =  (/                                       &
    & 0.48060688801085891E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_500_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_1800_para(1:6) = (/                                          &
    & -0.29117205777111277E-002_NER,       -0.33032767476037948E-002_NER,    &
    &  0.24485872054651786E-008_NER,        0.13302883775371241E-001_NER,    &
    &  0.40600979163933854E-007_NER,        0.44947563848805843E-013_NER/), &
    &                                                                       &
    & Cs3s1_1800_set2_para(1:7) =  (/                                       &
    & -0.33316408182589127E-002_NER,     -0.18748706802369042E+000_NER,      &
    & -0.52448300820408648E-007_NER,      0.25958922853975199E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_1800_set2_para(1:5) =  (/                                       &
    & 0.54847055164191530E-009_NER, 0.34823781823562049E-007_NER,          &
    & 0.40477075571456276E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_1800_set2_para(1:5) =  (/                                       &
    & -0.20675537750543899E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_1800_set2_para(1:4) =  (/                                       &
    & 0.34510856663214677E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_1800_set2_para(1:4) =  (/                                       &
    & 0.38741671955770587E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_1800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_1800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_1800_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_2200_para(1:6) = (/                                          &
    & -0.28144571255542625E-002_NER,       -0.35629348428037094E-002_NER,    &
    &  0.17500851589453367E-008_NER,        0.31598893287444585E-001_NER,    &
    &  0.33843595304747309E-007_NER,        0.40299720600767652E-013_NER/), &
    &                                                                       &
    & Cs3s1_2200_set2_para(1:7) =  (/                                       &
    & -0.98735381605867772E-003_NER,     -0.32219040606618982E+000_NER,      &
    & -0.53455148717113709E-007_NER,      0.16917641780342937E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_2200_set2_para(1:5) =  (/                                       &
    & 0.72738283590019237E-009_NER, 0.35049778409732377E-007_NER,          &
    & 0.42253040515554155E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_2200_set2_para(1:5) =  (/                                       &
    & -0.17631306747142849E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_2200_set2_para(1:4) =  (/                                       &
    & 0.38125683109519305E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_2200_set2_para(1:4) =  (/                                       &
    & 0.38725234381961927E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_2200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_2200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_2200_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_1400_para(1:6) = (/                                          &
    & -0.30600398721280560E-002_NER,       -0.29777019711699979E-002_NER,    &
    &  0.36896581100845415E-008_NER,        0.37926068416106878E-002_NER,    &
    &  0.45536569558327257E-007_NER,        0.50717306387187488E-013_NER/), &
    &                                                                       &
    & Cs3s1_1400_set2_para(1:7) =  (/                                       &
    & -0.67972155786817472E-002_NER,     0.10530600046579308E+001_NER,      &
    & -0.50890645284624999E-006_NER,      -0.69237464679038240E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_1400_set2_para(1:5) =  (/                                       &
    & -0.20567466947773754E-009_NER, 0.33919778081774643E-007_NER,          &
    & 0.46493792980760116E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_1400_set2_para(1:5) =  (/                                       &
    & -0.25177751473873542E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_1400_set2_para(1:4) =  (/                                       &
    & 0.28802178581073601E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_1400_set2_para(1:4) =  (/                                       &
    & 0.38798463693538108E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_1400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_1400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_1400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & MMWLY1s0_6400_para(1:6) = (/                                          &
    & -0.25051121960678767E-002_NER,       -0.48122147564501087E-002_NER,    &
    &  0.26417854442064662E-009_NER,        0.99340809418944942E+000_NER,    &
    &  -0.18589030758697550E-007_NER,       0.12737078480312943E-013_NER/), &
    &                                                                       &
    & Cs3s1_6400_set2_para(1:7) =  (/                                       &
    & -0.20213503934403580E-001_NER,     -0.18390101892012171E+003_NER,      &
    & -0.13954980761466109E-004_NER,     -0.14466176405133300E-004_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_6400_set2_para(1:5) =  (/                                       &
    & 0.11177324067944568E-009_NER, 0.46126167431461877E-007_NER,          &
    & 0.13486092587141029E-014_NER, 0.0_NER, 0.0_NER/),                    &
    &                                                                       &
    & Cs3p2_6400_set2_para(1:5) =  (/                                       &
    & -0.78864193038394094E-009_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),  &
    &                                                                       &
    & Cs3p1_6400_set2_para(1:4) =  (/                                       &
    & 0.48655392525006389E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1p1_6400_set2_para(1:4) =  (/                                       &
    & 0.38713462469008114E-008_NER, 0.0_NER, 0.0_NER, 0.0_NER/),            &
    &                                                                       &
    & Cs1d2_6400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d2_6400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/),                                 &
    &                                                                       &
    & Cs3d3_6400_set2_para(1:4) =  (/                                       &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & WOTPE3s1_800_para(1:7) = (/                                           &
    &  0.29003872118011503E-002_NER,    -0.13370417535128715E-001_NER,      &
    &  0.13883318499795604E-007_NER,     0.16814122475763160E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_1000_para(1:7) = (/                                          &
    &  0.48765294272804342E-001_NER,     0.11087319415123418E+001_NER,      &
    &  0.10236236937585180E-006_NER,    -0.12903651202445730E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_1200_para(1:7) = (/                                          &
    & -0.12668825908618342E-001_NER,    -0.29763547219582014E-001_NER,      &
    & -0.44281731090811596E-007_NER,    -0.97068151000320441E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_1300_para(1:7) = (/                                          &
    & -0.87385234733682038E-002_NER,     0.66560316692963584E-001_NER,      &
    & -0.51515889611697102E-007_NER,    -0.96264819217019688E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_1400_para(1:7) = (/                                          &
    & -0.67972155786817472E-002_NER,     0.41592405689932110E+000_NER,      &
    & -0.17078616738660048E-006_NER,    -0.27550310098404955E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_1600_para(1:7) = (/                                          &
    & -0.46750075253685118E-002_NER,    -0.12784306743582069E+000_NER,      &
    &  0.31824237258764145E-007_NER,     0.40824378957278200E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_2000_para(1:7) = (/                                          &
    & -0.21944206584707093E-002_NER,    -0.13550053621688371E+000_NER,      &
    &  0.18513675745254923E-007_NER,     0.21428304945017605E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_2200_para(1:7) = (/                                          &
    & -0.98735381605867772E-003_NER,    -0.23378222713067850E+000_NER,      &
    &  0.23736873929914491E-007_NER,     0.28362917420481648E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_2400_para(1:7) = (/                                          &
    &  0.62485155018945579E-003_NER,    -0.40999120464946254E+000_NER,      &
    &  0.32747550207470601E-007_NER,     0.37039643324965198E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_2800_para(1:7) = (/                                          &
    &  0.13118803810359370E-001_NER,    -0.12112320686487544E+003_NER,      &
    & -0.13817734329847372E-005_NER,     0.19179560038465822E-005_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_3200_para(1:7) = (/                                          &
    & -0.12838596032815119E-001_NER,    -0.16816166753210682E+001_NER,      &
    & -0.10089798772883788E-006_NER,    -0.18715221493839820E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_4000_para(1:7) = (/                                          &
    & -0.40755330091299634E-002_NER,     0.28450091151755035E+001_NER,      &
    & -0.10172612563303106E-006_NER,    -0.18184227297155841E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_4800_para(1:7) = (/                                          &
    & -0.16823198672076066E-002_NER,     0.29208616069861893E+002_NER,      &
    & -0.41378864807727690E-006_NER,    -0.88878109730960812E-006_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_400_para(1:7) = (/                                          &
    & -0.45148506612549544E-002_NER,    -0.21327102441374320E-002_NER,      &
    &  0.15055828254460183E-007_NER,     0.14188883072635077E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/),                                         &
    &                                                                       &
    & WOTPE3s1_600_para(1:7) = (/                                          &
    & -0.14442429421521084E-002_NER,    -0.46146175474309570E-002_NER,      &
    &  0.10549357099250846E-007_NER,     0.12725937868281851E-007_NER,      &
    &  0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter ::  &
    & WOTPE1s0_800_para(1:6) = (/                                            &
    & -0.35387941376915466E-002_NER,       -0.22691807277402604E-002_NER,    &
    &  0.87320693514134414E-008_NER,       -0.23246842289800073E-002_NER,    &
    &  0.19044544660693673E-007_NER,       -0.70811127342129122E-013_NER/),  &
    & WOTPE1s0_1000_para(1:6) = (/                                           &
    & -0.33181124691830394E-002_NER,       -0.25483306792893298E-002_NER,    &
    &  0.62547444546995307E-008_NER,       -0.48295900919165221E-002_NER,    &
    &  0.23859768813464021E-007_NER,       -0.63946308594804029E-013_NER/),  &
    & WOTPE1s0_1200_para(1:6) = (/                                           &
    & -0.31686486633502316E-002_NER,       -0.27797190987699032E-002_NER,    &
    &  0.47143394436464893E-008_NER,       -0.58115869591442966E-002_NER,    &
    &  0.19975482435472168E-007_NER,       -0.32299179762461904E-013_NER/),  &
    & WOTPE1s0_1400_para(1:6) = (/                                           &
    & -0.30600398721280560E-002_NER,       -0.29777019711699979E-002_NER,    &
    &  0.36896581100845415E-008_NER,       -0.81096123559724047E-002_NER,    &
    &  0.20128793885675478E-007_NER,       -0.24326583771108012E-013_NER/),  &
    & WOTPE1s0_1600_para(1:6) = (/                                           &
    & -0.29771946913657479E-002_NER,       -0.31504605152330055E-002_NER,    &
    &  0.29719305488873427E-008_NER,       -0.12050803588911889E-001_NER,    &
    &  0.22496855881263524E-007_NER,       -0.23244657000454164E-013_NER/),  &
    & WOTPE1s0_2000_para(1:6) = (/                                           &
    & -0.28585563548057182E-002_NER,       -0.34398508249289133E-002_NER,    &
    &  0.20545654031144380E-008_NER,       -0.64456883897513750E-002_NER,    &
    &  0.77612380165128946E-008_NER,        0.15615982725719953E-014_NER/),  &
    & WOTPE1s0_2200_para(1:6) = (/                                           &
    & -0.28144571255542625E-002_NER,       -0.35629348428037094E-002_NER,    &
    &  0.17500851589453367E-008_NER,       -0.87627478477502411E-002_NER,    &
    &  0.85147773244301861E-008_NER,        0.32436838688756727E-015_NER/),  &
    & WOTPE1s0_2400_para(1:6) = (/                                           &
    & -0.27772405429084906E-002_NER,       -0.36746014859360466E-002_NER,    &
    &  0.15096483095713912E-008_NER,       -0.80939869842212679E-002_NER,    &
    &  0.65021029412225862E-008_NER,        0.16768986935922078E-014_NER/),  &
    & WOTPE1s0_2800_para(1:6) = (/                                           &
    & -0.27177795837334149E-002_NER,       -0.32414183866159040E-002_NER,    &
    &  0.97026648371686457E-009_NER,       -0.86831931048748888E-002_NER,    &
    &  0.48918110129926344E-008_NER,        0.50905835017193730E-015_NER/),  &
    & WOTPE1s0_3200_para(1:6) = (/                                           &
    & -0.26722872208465821E-002_NER,       -0.40482615077028763E-002_NER,    &
    &  0.92099418329107802E-009_NER,       -0.18958389825670198E-001_NER,    &
    &  0.82121595533336700E-008_NER,       -0.26004845101515215E-015_NER/),  &
    & WOTPE1s0_4000_para(1:6) = (/                                           &
    & -0.26070755680917880E-002_NER,       -0.43029463919601009E-002_NER,    &
    &  0.61902751186230983E-009_NER,       -0.19344208319205660E-001_NER,    &
    &  0.50820222563057666E-008_NER,        0.47657359204835328E-015_NER/),  &
    & WOTPE1s0_4800_para(1:6) = (/                                           &
    & -0.25624472971879611E-002_NER,       -0.45099223506560072E-002_NER,    &
    &  0.44631438051037783E-009_NER,        0.62964915856283049E-002_NER,    &
    & -0.17391373689254603E-008_NER,        0.15368046320172716E-014_NER/),  &
    & WOTPE1s0_400_para(1:6) = (/                                            &
    & -0.46345112093629202E-002_NER,       -0.13812599401879487E-002_NER,    &
    &  0.21086806723513865E-007_NER,       -0.31201176173854529E-003_NER,    &
    &  0.13337776683270184E-007_NER,       -0.28333596126007128E-012_NER/),  &
    & WOTPE1s0_600_para(1:6) = (/                                            &
    & -0.39028006210685457E-002_NER,       -0.19113791264257359E-002_NER,    &
    &  0.13071651933734054E-007_NER,       -0.19211577654154400E-002_NER,    &
    &  0.26107387193442305E-007_NER,       -0.21303402517078763E-012_NER/)

  real(NER), parameter ::  &
    & WOTPE3p0_800_set2_para(1:5) =  (/                                      &
    & -0.39975686839150306E-007_NER, 0.51596589380781465E-006_NER,           &
    & -0.21109408774405163E-012_NER, 0.0_NER, 0.0_NER/),                     &
    & WOTPE3p0_1000_set2_para(1:5) =  (/                                     &
    & -0.60120587086032020E-008_NER, 0.11414186906350170E-007_NER,           &
    & -0.67519961016966236E-014_NER, 0.0_NER, 0.0_NER/),                     &
    & WOTPE3p0_1200_set2_para(1:5) =  (/                                     &
    & -0.15607073425249616E-008_NER, 0.89059783305349159E-009_NER,           &
    & -0.47467961807755260E-015_NER, 0.0_NER, 0.0_NER/),                     &
    & WOTPE3p0_1400_set2_para(1:5) =  (/                                     &
    & -0.20567466947773754E-009_NER, -0.44205036410298323E-009_NER,          &
    & 0.20355725724188493E-015_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_1600_set2_para(1:5) =  (/                                     &
    & 0.31968616606084914E-009_NER, -0.67064563544192988E-009_NER,           &
    & 0.27179874369641636E-015_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_2000_set2_para(1:5) =  (/                                     &
    & 0.66001753049767985E-009_NER, -0.71288376202444310E-009_NER,           &
    & 0.26127428086141483E-015_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_2200_set2_para(1:5) =  (/                                     &
    & 0.72738283590019237E-009_NER, -0.70599504671998807E-009_NER,           &
    & 0.28978035995060151E-015_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_2400_set2_para(1:5) =  (/                                     &
    & 0.78982690154086563E-009_NER, -0.56227008995922139E-009_NER,           &
    & 0.37468374176881846E-015_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_2800_set2_para(1:5) =  (/                                     &
    & 0.10708384950704921E-008_NER,  0.86667963204304351E-008_NER,           &
    & 0.13677810633530822E-014_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_3200_set2_para(1:5) =  (/                                     &
    & 0.11323845800451927E-007_NER, -0.14074802862188733E-004_NER,           &
    & -0.41039823197544460E-013_NER, 0.0_NER, 0.0_NER/),                     &
    & WOTPE3p0_4000_set2_para(1:5) =  (/                                     &
    & -0.68226727804105308E-010_NER,  -0.71019194513511619E-008_NER,         &
    & 0.19984156396339779E-015_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_4800_set2_para(1:5) =  (/                                     &
    & 0.62810050992021014E-010_NER, -0.24943224090949278E-008_NER,           &
    & 0.74113456328065310E-016_NER, 0.0_NER, 0.0_NER/),                      &
    & WOTPE3p0_400_set2_para(1:5) =  (/                                      &
    & 0.30164174135677295E-007_NER,  0.29078598732007712E-008_NER,           &
    & -0.68972755152071681E-013_NER, 0.0_NER, 0.0_NER/),                     &
    & WOTPE3p0_600_set2_para(1:5) =  (/                                      &
    & 0.47638609114073750E-007_NER, -0.10727377466375298E-006_NER,           &
    & -0.28593803662354769E-012_NER, 0.0_NER, 0.0_NER/)


  contains

  subroutine read_chengdu_version(verinfo)

    character(len = *), intent(inout) :: verinfo

    verinfo = verstr

  end subroutine

  subroutine init_PCs4chengdu()

    PC_gA = 1.29_NER
    PC_mN = 939.0_NER
    PC_mpi = 138.0_NER
    PC_fpi = 92.4_NER

  end subroutine init_PCs4chengdu

  subroutine chengdu_dispatch(pot_name)

    character(len = *), intent(in)  :: pot_name

    hotpot_name = pot_name

  end subroutine chengdu_dispatch

  function chengdu_hotpot_is_Edep(L, S, J)

    logical :: chengdu_hotpot_is_Edep
    integer, intent(in)     :: L, S, J

    chengdu_hotpot_is_Edep = .false.
    if (L == 0 .and. S == 0 .and. J == 0) then
      if (hotpot_name == "chengdu_LBWDIB_600")          &
        & chengdu_hotpot_is_Edep = .true.
      if (hotpot_name == "chengdu_LBWDIB_cmplx_600") then
        chengdu_hotpot_is_Edep = .true.
      end if
    end if

  end function chengdu_hotpot_is_Edep

  subroutine chengdu_hotpot_basic(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    select case(hotpot_name)
      case ("chengdu_DLSPR")
        call chengdu_DLSPR(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_DLSPR_350")
        call chengdu_DLSPR_350(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_DLSPR_400")
        call chengdu_DLSPR_400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_DLSPR_800")
        call chengdu_DLSPR_800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_DLSPR_1600")
        call chengdu_DLSPR_1600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_DLSPR_3200")
        call chengdu_DLSPR_3200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_pionless_1000")
        call chengdu_pionless_1000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_350")
        call chengdu_MMWLY_350(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_400")
        call chengdu_MMWLY_400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_600")
        call chengdu_MMWLY_600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_800")
        call chengdu_MMWLY_800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_1600")
        call chengdu_MMWLY_1600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_3200")
        call chengdu_MMWLY_3200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_1200")
        call chengdu_MMWLY_1200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_2000")
        call chengdu_MMWLY_2000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_2400")
        call chengdu_MMWLY_2400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_2800")
        call chengdu_MMWLY_2800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_3600")
        call chengdu_MMWLY_3600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_4000")
        call chengdu_MMWLY_4000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_4400")
        call chengdu_MMWLY_4400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_4800")
        call chengdu_MMWLY_4800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_1000")
        call chengdu_MMWLY_1000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_1300")
        call chengdu_MMWLY_1300(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_1500")
        call chengdu_MMWLY_1500(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_2500")
        call chengdu_MMWLY_2500(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_3000")
        call chengdu_MMWLY_3000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_500")
        call chengdu_MMWLY_500(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_1800")
        call chengdu_MMWLY_1800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_2200")
        call chengdu_MMWLY_2200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_900")
        call chengdu_MMWLY_900(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_1400")
        call chengdu_MMWLY_1400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_6400")
        call chengdu_MMWLY_6400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_VLbar_1200")
        call chengdu_VLbar_1200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_VLbar_1300")
        call chengdu_VLbar_1300(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_VLbar_1400")
        call chengdu_VLbar_1400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_VLbar_1600")
        call chengdu_VLbar_1600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_400")
        call chengdu_WOTPE_400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_600")
        call chengdu_WOTPE_600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_800")
        call chengdu_WOTPE_800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_1000")
        call chengdu_WOTPE_1000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_1200")
        call chengdu_WOTPE_1200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_1300")
        call chengdu_WOTPE_1300(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_1400")
        call chengdu_WOTPE_1400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_1600")
        call chengdu_WOTPE_1600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_2000")
        call chengdu_WOTPE_2000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_2200")
        call chengdu_WOTPE_2200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_2400")
        call chengdu_WOTPE_2400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_2800")
        call chengdu_WOTPE_2800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_3200")
        call chengdu_WOTPE_3200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_4000")
        call chengdu_WOTPE_4000(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_WOTPE_4800")
        call chengdu_WOTPE_4800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MWPLL_400")
        call chengdu_MWPLL_400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_600_diagnose")
        call chengdu_MMWLY_600_diagnose(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_600_only_1S0_HO")
        call chengdu_MMWLY_600_only_1S0_HO(L, S, J, uptoQn, po, pi, potval)
      case default
        write (standard_error_unit, '(a)')  &
          & trim(hotpot_name)//" not recognized"
        potval = 0.0_NER
    end select

  end subroutine chengdu_hotpot_basic

  subroutine chengdu_hotpot_mixed(L, S, J, uptoQn, para, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: para(:), po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    integer :: ii
    real(NER) :: tmpval(1:2, 1:2)

    potval = 0.0_NER
    call chengdu_hotpot_basic(L, S, J, 0, po, pi, potval)
    do ii = 1, uptoQn
      call chengdu_hotpot_basic(L, S, J, ii, po, pi, tmpval)
      potval = potval + para(ii)*tmpval
    end do

  end subroutine chengdu_hotpot_mixed

  subroutine chengdu_cmplx_hotpot(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)       :: L, S, J, uptoQn
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    select case(hotpot_name)
      case ("chengdu_MMWLY_cmplx_350")
        call chengdu_MMWLY_cmplx_350(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_cmplx_400")
        call chengdu_MMWLY_cmplx_400(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_cmplx_600")
        call chengdu_MMWLY_cmplx_600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_cmplx_800")
        call chengdu_MMWLY_cmplx_800(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_cmplx_1600")
        call chengdu_MMWLY_cmplx_1600(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_MMWLY_cmplx_3200")
        call chengdu_MMWLY_cmplx_3200(L, S, J, uptoQn, po, pi, potval)
      case ("chengdu_LBWDIB_cmplx_600")
        if (L == 0 .and. S == 0 .and. J == 0) then
          write (standard_error_unit, '(a)')  &
          & trim(hotpot_name)//" is energy dependent in 1S0"
        else
          call chengdu_MMWLY_cmplx_600(L, S, J, uptoQn, po, pi, potval)
        end if
      case default
        write (standard_error_unit, '(a)')  &
          & trim(hotpot_name)//" not recognized"
        potval = 0.0_NER
    end select

  end subroutine chengdu_cmplx_hotpot

  subroutine chengdu_cmplx_hotpot_mixed(L, S, J, uptoQn, para, po, pi, potval)

    integer, intent(in)       :: L, S, J, uptoQn
    real(NER), intent(in)     :: para(:)
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    integer       :: ii
    complex(NEC)  :: tmpval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call chengdu_cmplx_hotpot(L, S, J, 0, po, pi, potval)
    do ii = 1, uptoQn
      call chengdu_cmplx_hotpot(L, S, J, ii, po, pi, tmpval)
      potval = potval + para(ii)*tmpval
    end do

  end subroutine chengdu_cmplx_hotpot_mixed

  subroutine chengdu_hotpot_Edep(L, S, J, uptoQn, Ecm, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi, Ecm
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    select case(hotpot_name)
      case ("chengdu_LBWDIB_600")
        call chengdu_LBWDIB_600(L, S, J, uptoQn, Ecm, po, pi, potval)
      case default
        write (standard_error_unit, '(a)')  &
          & trim(hotpot_name)//" not recognized"
        potval = 0.0_NER
    end select

  end subroutine chengdu_hotpot_Edep

  subroutine chengdu_cmplx_hotpot_Edep(L, S, J, uptoQn, Ecm, po, pi, potval)

    integer, intent(in)         :: L, S, J, uptoQn
    complex(NEC), intent(in)    :: Ecm, po, pi
    complex(NEC), intent(out)   :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    select case(hotpot_name)
      case ("chengdu_LBWDIB_cmplx_600")
        call chengdu_LBWDIB_cmplx_600(L, S, J, uptoQn, Ecm, po, pi, potval)
      case default
        write (standard_error_unit, '(a)')  &
          & trim(hotpot_name)//" not recognized"
        potval = 0.0_NER
    end select

  end subroutine chengdu_cmplx_hotpot_Edep

  subroutine chengdu_cmplx_hotpot_mixed_Edep(L, S, J, uptoQn, para, Ecm, po, pi, potval)

    integer, intent(in)         :: L, S, J, uptoQn
    real(NER), intent(in)       :: para(:)
    complex(NEC), intent(in)    :: Ecm, po, pi
    complex(NEC), intent(out)   :: potval(1:NCH, 1:NCH)

    integer :: ii
    complex(NEC) :: tmpval(1:2, 1:2)

    potval = 0.0_NER
    do ii = 0, uptoQn
      call chengdu_cmplx_hotpot_Edep(L, S, J, uptoQn, Ecm, po, pi, tmpval)
      if (ii == 0) then
        potval = tmpval
      else
        potval = potval + para(ii)*tmpval
      end if
    end do

  end subroutine chengdu_cmplx_hotpot_mixed_Edep

  subroutine init_dlspr_para(para1s0, para3s1, para3p0, para1p1, &
    & para3p1, para3p2, para1d2, para3d2, para3d3)

    real(NER),intent(in) :: para1s0(1:14), para3s1(1:7), para3p0(1:5), &
    & para1p1(1:4), para3p1(1:4), para3p2(1:5), para1d2(1:4),          &
    & para3d2(1:4), para3d3(1:4)

    dlspr1s0 = para1s0
    Cs3s1    = para3s1
    Cs3p0    = para3p0
    Cs3p2    = para3p2
    Cs1p1    = para1p1
    Cs3p1    = para3p1
    Cs1d2    = para1d2
    Cs3d2    = para3d2
    Cs3d3    = para3d3

  end subroutine init_dlspr_para

  subroutine get_chengdu_DLSPR(L, S, J, uptoQn, po, pi, Mambda, potval)

    integer,intent(in)    :: L, S, J, uptoQn
    real(NER),intent(in)  :: po, pi, Mambda
    real(NER),intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)
    real(NER) :: para0(1:1)

    para0 = 0.0_NER
    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_DLSPR : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
      lsj%L = lsj%J - 1
      call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_sprbl_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & dlspr1s0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
              & CHNID_3D2_PP, CHNID_1D2_PP)
            call OPE_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
              & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine get_chengdu_DLSPR

  subroutine init_mmwly_para(mmwlypara, para3s1, para3p0, para1p1, &
    & para3p1, para3p2, para1d2, para3d2, para3d3)

    real(NER), intent(in) :: mmwlypara(1:6), para3s1(1:7), para3p0(1:5), &
    & para1p1(1:4), para3p1(1:4), para3p2(1:5), para1d2(1:4),            &
    & para3d2(1:4), para3d3(1:4)

    mmwly1s0 = mmwlypara
    Cs3s1    = para3s1
    Cs3p0    = para3p0
    Cs3p2    = para3p2
    Cs1p1    = para1p1
    Cs3p1    = para3p1
    Cs1d2    = para1d2
    Cs3d2    = para3d2
    Cs3d3    = para3d3

  end subroutine init_mmwly_para

  subroutine get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi, Mambda
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)
    real(NER) :: para0(1:1)

    para0 = 0.0_NER
    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_MMWLY : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
      lsj%L = lsj%J - 1
      call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VNLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
              & CHNID_3D2_PP, CHNID_1D2_PP)
            call OPE_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
              & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(2)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VN2LO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P1_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p1, po, pi, tmppotval)
          case (CHNID_1P1_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs1p1, po, pi, tmppotval)
          case (CHNID_3P2_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p2, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine get_chengdu_MMWLY

! use Vbar constructed by VLO_withC0 and spuriours bound state in 3s1

  subroutine get_chengdu_VLbar(L, S, J, uptoQn, po, pi, Mambda, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi, Mambda
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)
    real(NER) :: para0(1:1)

    para0 = 0.0_NER
    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_MMWLY : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
      lsj%L = lsj%J - 1
      call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VbarLO_withc0_chengdu_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VbarLO_3p0_withc0_chengdu_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VNLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
              & CHNID_3D2_PP, CHNID_1D2_PP)
            call OPE_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
              & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(2)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VN2LO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P1_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p1, po, pi, tmppotval)
          case (CHNID_1P1_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs1p1, po, pi, tmppotval)
          case (CHNID_3P2_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p2, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine get_chengdu_VLbar

! use Vbar constructed by VLO_withC0 and spuriours bound state in 3s1

  subroutine get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi, Mambda
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)
    real(NER) :: para0(1:1)

    para0 = 0.0_NER
    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_MMWLY : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
      lsj%L = lsj%J - 1
      call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VNLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
              & CHNID_3D2_PP, CHNID_1D2_PP)
            call OPE_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
              & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(2)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VN2LO_MMWLY_SHRG_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VN2LO_withc0_SHRG_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VN2LO_withc0_SHRG_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P1_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p1, po, pi, tmppotval)
          case (CHNID_1P1_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs1p1, po, pi, tmppotval)
          case (CHNID_3P2_PP)
            call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p2, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine get_chengdu_WOTPE

  subroutine get_chengdu_MMWLY_cmplx(L, S, J, uptoQn, po, pi, Mambda, potval)

    integer, intent(in)       :: L, S, J, uptoQn
    real(NER), intent(in)     :: Mambda
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol)  :: lsj
    integer       :: chnid
    complex(NEC)  :: tmppotval(1:2, 1:2)
    real(NER)     :: para0(1:1)

    para0 = 0.0_NER
    call init_PCs4chengdu()
    ! debug
    PC_GAUSS_CMPLX = 4
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_MMWLY_cmplx : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
        lsj%L = lsj%J - 1
        call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_MMWLY_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VLO_withc0_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VLO_withc0_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VNLO_MMWLY_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP,     &
              & CHNID_3D2_PP, CHNID_1D2_PP)
            call OPE_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
              & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine get_chengdu_MMWLY_cmplx

  subroutine chengdu_DLSPR_350(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 350.0_NER
    call init_dlspr_para(spr1s0_350_set2_para, Cs3s1_350_set2_para, Cs3p0_350_set2_para,     &
      & Cs1p1_350_set2_para, Cs3p1_350_set2_para, Cs3p2_350_set2_para, Cs1d2_350_set2_para   &
      & , Cs3d2_350_set2_para, Cs3d3_350_set2_para)
    call get_chengdu_DLSPR(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_DLSPR_350

  subroutine chengdu_DLSPR_400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 400.0_NER
    call init_dlspr_para(spr1s0_400_set2_para, Cs3s1_400_set2_para, Cs3p0_400_set2_para,     &
      & Cs1p1_400_set2_para, Cs3p1_400_set2_para, Cs3p2_400_set2_para, Cs1d2_400_set2_para   &
      & , Cs3d2_400_set2_para, Cs3d3_400_set2_para)
    call get_chengdu_DLSPR(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_DLSPR_400

  subroutine chengdu_DLSPR(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 600.0_NER
    call init_dlspr_para(spr1s0_600_set2_para, Cs3s1_600_set2_para, Cs3p0_600_set2_para,   &
      & Cs1p1_600_set2_para, Cs3p1_600_set2_para, Cs3p2_600_set2_para, Cs1d2_600_set2_para &
      & , Cs3d2_600_set2_para, Cs3d3_600_set2_para)
    call get_chengdu_DLSPR(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_DLSPR

  subroutine chengdu_DLSPR_800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 800.0_NER
    call init_dlspr_para(spr1s0_800_set2_para, Cs3s1_800_set2_para, Cs3p0_800_set2_para,   &
      & Cs1p1_800_set2_para, Cs3p1_800_set2_para, Cs3p2_800_set2_para, Cs1d2_800_set2_para &
      & , Cs3d2_800_set2_para, Cs3d3_800_set2_para)
    call get_chengdu_DLSPR(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_DLSPR_800

  subroutine chengdu_DLSPR_1600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1600.0_NER
    call init_dlspr_para(spr1s0_1600_set2_para, Cs3s1_1600_set2_para,          &
      & Cs3p0_1600_set2_para, Cs1p1_1600_set2_para, Cs3p1_1600_set2_para,      &
      & Cs3p2_1600_set2_para, Cs1d2_1600_set2_para, Cs3d2_1600_set2_para,      &
      & Cs3d3_1600_set2_para)
    call get_chengdu_DLSPR(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_DLSPR_1600

  subroutine chengdu_DLSPR_3200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 3200.0_NER
    call init_dlspr_para(spr1s0_3200_set2_para, Cs3s1_3200_set2_para, Cs3p0_3200_set2_para,     &
      & Cs1p1_3200_set2_para, Cs3p1_3200_set2_para, Cs3p2_3200_set2_para, Cs1d2_3200_set2_para  &
      & , Cs3d2_3200_set2_para, Cs3d3_3200_set2_para)
    call get_chengdu_DLSPR(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_DLSPR_3200

  subroutine chengdu_pionless_1000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda
    real(NER) :: pionlesspara1s0(1:1), pionlesspara3s1(1:1)

    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_pionless_1000 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
        lsj%L = lsj%J - 1
        call convert_lsj_to_chnid(lsj, chnid)
    end if

    ! pionlesspara1s0(1) =  -1.56315E-003_NER   ! Fitted to a = -23.7 fm
    ! ! pionlesspara3s1(1) =  -1.8275E-003_NER    ! Fitted to a = +5.42 fm
    ! pionlesspara3s1(1) =  -1.7862E-003_NER    ! Fitted to Bd = 2.22 MeV

    pionlesspara1s0(1) =  -1.56315E-003_NER   ! Fitted to a = -23.7 fm
    pionlesspara3s1(1) =  -1.7560E-003_NER    ! Fitted to a = +5.42 fm
    ! pionlesspara3s1(1) =  -1.7862E-003_NER    ! Fitted to Bd = 2.22 MeV

    Mambda = 1000.0_NER
    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_withc0_pionless_NP_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
          & pionlesspara1s0, po, pi, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_pionless_NP_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
          & pionlesspara3s1, po, pi, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_pionless_1000

  subroutine chengdu_LBWDIB_600(L, S, J, uptoQn, Ecm, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi, Ecm
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda, para(1:20)

    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_LBWDIB_600 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
        lsj%L = lsj%J - 1
        call convert_lsj_to_chnid(lsj, chnid)
    end if

    para(1) =  0.52106222046242991E+002_NER
    para(2) = -0.82992540492249800E-001_NER

    Mambda = 600.0_NER
    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_LBWDIB_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          & para, Ecm, po, pi, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          & Cs3s1_600_set2_para, po, pi, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          & Cs3p0_600_set2_para, po, pi, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_LBWDIB_600

  subroutine chengdu_LBWDIB_cmplx_600(L, S, J, uptoQn, Ecm, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    complex(NEC), intent(in)   :: po, pi, Ecm
    complex(NEC), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    complex(NEC)  :: tmppotval(1:2, 1:2)
    real(NER)     :: Mambda, para(1:20)

    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_LBWDIB_cmplx_600 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
        lsj%L = lsj%J - 1
        call convert_lsj_to_chnid(lsj, chnid)
    end if

    para(1) =  0.52106222046242991E+002_NER
    para(2) = -0.82992540492249800E-001_NER

    Mambda = 600.0_NER
    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_LBWDIB_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          & para, Ecm, po, pi, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          & Cs3s1_600_set2_para, po, pi, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_cmplx_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          & Cs3p0_600_set2_para, po, pi, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_LBWDIB_cmplx_600

  subroutine chengdu_MMWLY_350(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 350.0_NER
    call init_mmwly_para(MMWLY1s0_350_para, Cs3s1_350_set2_para,               &
      & Cs3p0_350_set2_para, Cs1p1_350_set2_para, Cs3p1_350_set2_para,         &
      & Cs3p2_350_set2_para, Cs1d2_350_set2_para, Cs3d2_350_set2_para,         &
      & Cs3d3_350_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_350

  subroutine chengdu_MMWLY_400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 400.0_NER
    call init_mmwly_para(MMWLY1s0_400_para, Cs3s1_400_set2_para, Cs3p0_400_set2_para,       &
      & Cs1p1_400_set2_para, Cs3p1_400_set2_para, Cs3p2_400_set2_para, Cs1d2_400_set2_para  &
      & , Cs3d2_400_set2_para, Cs3d3_400_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_400

  subroutine chengdu_MMWLY_600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 600.0_NER
    call init_mmwly_para(MMWLY1s0_600_para, Cs3s1_600_set2_para, Cs3p0_600_set2_para,       &
      & Cs1p1_600_set2_para, Cs3p1_600_set2_para, Cs3p2_600_set2_para, Cs1d2_600_set2_para  &
      & , Cs3d2_600_set2_para, Cs3d3_600_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_600

  subroutine chengdu_MMWLY_800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 800.0_NER
    call init_mmwly_para(MMWLY1s0_800_para, Cs3s1_800_set2_para, Cs3p0_800_set2_para,       &
      & Cs1p1_800_set2_para, Cs3p1_800_set2_para, Cs3p2_800_set2_para, Cs1d2_800_set2_para  &
      & , Cs3d2_800_set2_para, Cs3d3_800_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_800

  subroutine chengdu_MMWLY_1600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1600.0_NER
    call init_mmwly_para(MMWLY1s0_1600_para, Cs3s1_1600_set2_para, Cs3p0_1600_set2_para,        &
      & Cs1p1_1600_set2_para, Cs3p1_1600_set2_para, Cs3p2_1600_set2_para, Cs1d2_1600_set2_para  &
      & , Cs3d2_1600_set2_para, Cs3d3_1600_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_1600

  subroutine chengdu_MMWLY_3200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 3200.0_NER
    call init_mmwly_para(MMWLY1s0_3200_para, Cs3s1_3200_set2_para, Cs3p0_3200_set2_para,        &
      & Cs1p1_3200_set2_para, Cs3p1_3200_set2_para, Cs3p2_3200_set2_para, Cs1d2_3200_set2_para  &
      & , Cs3d2_3200_set2_para, Cs3d3_3200_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_3200
!add mmwly 1200/2000/2400/2800

  subroutine chengdu_MMWLY_1200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1200.0_NER
    call init_mmwly_para(MMWLY1s0_1200_para, Cs3s1_1200_set2_para, Cs3p0_1200_set2_para,        &
      & Cs1p1_1200_set2_para, Cs3p1_1200_set2_para, Cs3p2_1200_set2_para, Cs1d2_1200_set2_para  &
      & , Cs3d2_1200_set2_para, Cs3d3_1200_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_1200

  subroutine chengdu_MMWLY_2000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2000.0_NER
    call init_mmwly_para(MMWLY1s0_2000_para, Cs3s1_2000_set2_para, Cs3p0_2000_set2_para,        &
      & Cs1p1_2000_set2_para, Cs3p1_2000_set2_para, Cs3p2_2000_set2_para, Cs1d2_2000_set2_para  &
      & , Cs3d2_2000_set2_para, Cs3d3_2000_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_2000

  subroutine chengdu_MMWLY_2400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2400.0_NER
    call init_mmwly_para(MMWLY1s0_2400_para, Cs3s1_2400_set2_para, Cs3p0_2400_set2_para,        &
      & Cs1p1_2400_set2_para, Cs3p1_2400_set2_para, Cs3p2_2400_set2_para, Cs1d2_2400_set2_para  &
      & , Cs3d2_2400_set2_para, Cs3d3_2400_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_2400

  subroutine chengdu_MMWLY_2800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2800.0_NER
    call init_mmwly_para(MMWLY1s0_2800_para, Cs3s1_2800_set2_para, Cs3p0_2800_set2_para,        &
      & Cs1p1_2800_set2_para, Cs3p1_2800_set2_para, Cs3p2_2800_set2_para, Cs1d2_2800_set2_para  &
      & , Cs3d2_2800_set2_para, Cs3d3_2800_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_2800

  subroutine chengdu_MMWLY_3600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 3600.0_NER
    call init_mmwly_para(MMWLY1s0_3600_para, Cs3s1_3600_set2_para, Cs3p0_3600_set2_para,        &
      & Cs1p1_3600_set2_para, Cs3p1_3600_set2_para, Cs3p2_3600_set2_para, Cs1d2_3600_set2_para  &
      & , Cs3d2_3600_set2_para, Cs3d3_3600_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_3600

  subroutine chengdu_MMWLY_4000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 4000.0_NER
    call init_mmwly_para(MMWLY1s0_4000_para, Cs3s1_4000_set2_para, Cs3p0_4000_set2_para,        &
      & Cs1p1_4000_set2_para, Cs3p1_4000_set2_para, Cs3p2_4000_set2_para, Cs1d2_4000_set2_para  &
      & , Cs3d2_4000_set2_para, Cs3d3_4000_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_4000

  subroutine chengdu_MMWLY_4400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 4400.0_NER
    call init_mmwly_para(MMWLY1s0_4400_para, Cs3s1_4400_set2_para, Cs3p0_4400_set2_para,        &
      & Cs1p1_4400_set2_para, Cs3p1_4400_set2_para, Cs3p2_4400_set2_para, Cs1d2_4400_set2_para  &
      & , Cs3d2_4400_set2_para, Cs3d3_4400_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_4400

  subroutine chengdu_MMWLY_4800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 4800.0_NER
    call init_mmwly_para(MMWLY1s0_4800_para, Cs3s1_4800_set2_para, Cs3p0_4800_set2_para,        &
      & Cs1p1_4800_set2_para, Cs3p1_4800_set2_para, Cs3p2_4800_set2_para, Cs1d2_4800_set2_para  &
      & , Cs3d2_4800_set2_para, Cs3d3_4800_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_4800

  subroutine chengdu_MMWLY_1000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1000.0_NER
    call init_mmwly_para(MMWLY1s0_1000_para, Cs3s1_1000_set2_para, Cs3p0_1000_set2_para,        &
      & Cs1p1_1000_set2_para, Cs3p1_1000_set2_para, Cs3p2_1000_set2_para, Cs1d2_1000_set2_para  &
      & , Cs3d2_1000_set2_para, Cs3d3_1000_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_1000

  subroutine chengdu_MMWLY_1300(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1300.0_NER
    call init_mmwly_para(MMWLY1s0_1300_para, Cs3s1_1300_set2_para, Cs3p0_1300_set2_para,        &
      & Cs1p1_1300_set2_para, Cs3p1_1300_set2_para, Cs3p2_1300_set2_para, Cs1d2_1300_set2_para  &
      & , Cs3d2_1300_set2_para, Cs3d3_1300_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_1300

  subroutine chengdu_MMWLY_1500(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1500.0_NER
    call init_mmwly_para(MMWLY1s0_1500_para, Cs3s1_1500_set2_para, Cs3p0_1500_set2_para,        &
      & Cs1p1_1500_set2_para, Cs3p1_1500_set2_para, Cs3p2_1500_set2_para, Cs1d2_1500_set2_para  &
      & , Cs3d2_1500_set2_para, Cs3d3_1500_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_1500

  subroutine chengdu_MMWLY_2500(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2500.0_NER
    call init_mmwly_para(MMWLY1s0_2500_para, Cs3s1_2500_set2_para, Cs3p0_2500_set2_para,        &
      & Cs1p1_2500_set2_para, Cs3p1_2500_set2_para, Cs3p2_2500_set2_para, Cs1d2_2500_set2_para  &
      & , Cs3d2_2500_set2_para, Cs3d3_2500_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_2500

  subroutine chengdu_MMWLY_3000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 3000.0_NER
    call init_mmwly_para(MMWLY1s0_3000_para, Cs3s1_3000_set2_para, Cs3p0_3000_set2_para,        &
      & Cs1p1_3000_set2_para, Cs3p1_3000_set2_para, Cs3p2_3000_set2_para, Cs1d2_3000_set2_para  &
      & , Cs3d2_3000_set2_para, Cs3d3_3000_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_3000

  subroutine chengdu_MMWLY_500(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 500.0_NER
    call init_mmwly_para(MMWLY1s0_500_para, Cs3s1_500_set2_para, Cs3p0_500_set2_para,        &
      & Cs1p1_500_set2_para, Cs3p1_500_set2_para, Cs3p2_500_set2_para, Cs1d2_500_set2_para  &
      & , Cs3d2_500_set2_para, Cs3d3_500_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_500

  subroutine chengdu_MMWLY_1800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1800.0_NER
    call init_mmwly_para(MMWLY1s0_1800_para, Cs3s1_1800_set2_para, Cs3p0_1800_set2_para,        &
      & Cs1p1_1800_set2_para, Cs3p1_1800_set2_para, Cs3p2_1800_set2_para, Cs1d2_1800_set2_para  &
      & , Cs3d2_1800_set2_para, Cs3d3_1800_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_1800

  subroutine chengdu_MMWLY_2200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2200.0_NER
    call init_mmwly_para(MMWLY1s0_2200_para, Cs3s1_2200_set2_para, Cs3p0_2200_set2_para,        &
      & Cs1p1_2200_set2_para, Cs3p1_2200_set2_para, Cs3p2_2200_set2_para, Cs1d2_2200_set2_para  &
      & , Cs3d2_2200_set2_para, Cs3d3_2200_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_2200

  subroutine chengdu_MMWLY_900(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 900.0_NER
    call init_mmwly_para(MMWLY1s0_900_para, Cs3s1_900_set2_para, Cs3p0_900_set2_para,        &
      & Cs1p1_900_set2_para, Cs3p1_900_set2_para, Cs3p2_900_set2_para, Cs1d2_900_set2_para  &
      & , Cs3d2_900_set2_para, Cs3d3_900_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_900

  subroutine chengdu_MMWLY_1400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1400.0_NER
    call init_mmwly_para(MMWLY1s0_1400_para, Cs3s1_1400_set2_para, Cs3p0_1400_set2_para,        &
      & Cs1p1_1400_set2_para, Cs3p1_1400_set2_para, Cs3p2_1400_set2_para, Cs1d2_1400_set2_para  &
      & , Cs3d2_1400_set2_para, Cs3d3_1400_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_1400

  subroutine chengdu_MMWLY_6400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 6400.0_NER
    call init_mmwly_para(MMWLY1s0_6400_para, Cs3s1_6400_set2_para, Cs3p0_6400_set2_para,        &
      & Cs1p1_6400_set2_para, Cs3p1_6400_set2_para, Cs3p2_6400_set2_para, Cs1d2_6400_set2_para  &
      & , Cs3d2_6400_set2_para, Cs3d3_6400_set2_para)
    call get_chengdu_MMWLY(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_6400

! ============================== !
! mmwly with Vbar routines
! ============================== !

  subroutine chengdu_VLbar_1200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1200.0_NER
    call init_mmwly_para(MMWLY1s0_1200_para, Cs3s1_1200_set2_para, Cs3p0_1200_set2_para,        &
      & Cs1p1_1200_set2_para, Cs3p1_1200_set2_para, Cs3p2_1200_set2_para, Cs1d2_1200_set2_para  &
      & , Cs3d2_1200_set2_para, Cs3d3_1200_set2_para)
    call get_chengdu_VLbar(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_VLbar_1200

  subroutine chengdu_VLbar_1300(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1300.0_NER
    call init_mmwly_para(MMWLY1s0_1300_para, Cs3s1_1300_set2_para, Cs3p0_1300_set2_para,        &
      & Cs1p1_1300_set2_para, Cs3p1_1300_set2_para, Cs3p2_1300_set2_para, Cs1d2_1300_set2_para  &
      & , Cs3d2_1300_set2_para, Cs3d3_1300_set2_para)
    call get_chengdu_VLbar(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_VLbar_1300

  subroutine chengdu_VLbar_1400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1400.0_NER
    call init_mmwly_para(MMWLY1s0_1400_para, Cs3s1_1400_set2_para, Cs3p0_1400_set2_para,        &
      & Cs1p1_1400_set2_para, Cs3p1_1400_set2_para, Cs3p2_1400_set2_para, Cs1d2_1400_set2_para  &
      & , Cs3d2_1400_set2_para, Cs3d3_1400_set2_para)
    call get_chengdu_VLbar(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_VLbar_1400

  subroutine chengdu_VLbar_1600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1600.0_NER
    call init_mmwly_para(MMWLY1s0_1600_para, Cs3s1_1600_set2_para, Cs3p0_1600_set2_para,        &
      & Cs1p1_1600_set2_para, Cs3p1_1600_set2_para, Cs3p2_1600_set2_para, Cs1d2_1600_set2_para  &
      & , Cs3d2_1600_set2_para, Cs3d3_1600_set2_para)
    call get_chengdu_VLbar(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_VLbar_1600

  subroutine chengdu_WOTPE_400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 400.0_NER
    call init_mmwly_para(WOTPE1s0_400_para, WOTPE3s1_400_para, WOTPE3p0_400_set2_para,        &
      & Cs1p1_400_set2_para, Cs3p1_400_set2_para, Cs3p2_400_set2_para, Cs1d2_400_set2_para  &
      & , Cs3d2_400_set2_para, Cs3d3_400_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_400

  subroutine chengdu_WOTPE_600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 600.0_NER
    call init_mmwly_para(WOTPE1s0_600_para, WOTPE3s1_600_para, WOTPE3p0_600_set2_para,        &
      & Cs1p1_600_set2_para, Cs3p1_600_set2_para, Cs3p2_600_set2_para, Cs1d2_600_set2_para  &
      & , Cs3d2_600_set2_para, Cs3d3_600_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_600

  subroutine chengdu_WOTPE_800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 800.0_NER
    call init_mmwly_para(WOTPE1s0_800_para, WOTPE3s1_800_para, WOTPE3p0_800_set2_para,        &
      & Cs1p1_800_set2_para, Cs3p1_800_set2_para, Cs3p2_800_set2_para, Cs1d2_800_set2_para  &
      & , Cs3d2_800_set2_para, Cs3d3_800_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_800

  subroutine chengdu_WOTPE_1000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1000.0_NER
    call init_mmwly_para(WOTPE1s0_1000_para, WOTPE3s1_1000_para, WOTPE3p0_1000_set2_para,        &
      & Cs1p1_1000_set2_para, Cs3p1_1000_set2_para, Cs3p2_1000_set2_para, Cs1d2_1000_set2_para  &
      & , Cs3d2_1000_set2_para, Cs3d3_1000_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_1000

  subroutine chengdu_WOTPE_1200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1200.0_NER
    call init_mmwly_para(WOTPE1s0_1200_para, WOTPE3s1_1200_para, WOTPE3p0_1200_set2_para,        &
      & Cs1p1_1200_set2_para, Cs3p1_1200_set2_para, Cs3p2_1200_set2_para, Cs1d2_1200_set2_para  &
      & , Cs3d2_1200_set2_para, Cs3d3_1200_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_1200

  subroutine chengdu_WOTPE_1300(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1300.0_NER
    call init_mmwly_para(MMWLY1s0_1300_para, WOTPE3s1_1300_para, Cs3p0_1300_set2_para,        &
      & Cs1p1_1300_set2_para, Cs3p1_1300_set2_para, Cs3p2_1300_set2_para, Cs1d2_1300_set2_para  &
      & , Cs3d2_1300_set2_para, Cs3d3_1300_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_1300

  subroutine chengdu_WOTPE_1400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1400.0_NER
    call init_mmwly_para(WOTPE1s0_1400_para, WOTPE3s1_1400_para, WOTPE3p0_1400_set2_para,        &
      & Cs1p1_1400_set2_para, Cs3p1_1400_set2_para, Cs3p2_1400_set2_para, Cs1d2_1400_set2_para  &
      & , Cs3d2_1400_set2_para, Cs3d3_1400_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_1400

  subroutine chengdu_WOTPE_1600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 1600.0_NER
    call init_mmwly_para(WOTPE1s0_1600_para, WOTPE3s1_1600_para, WOTPE3p0_1600_set2_para,        &
      & Cs1p1_1600_set2_para, Cs3p1_1600_set2_para, Cs3p2_1600_set2_para, Cs1d2_1600_set2_para  &
      & , Cs3d2_1600_set2_para, Cs3d3_1600_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_1600

  subroutine chengdu_WOTPE_2000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2000.0_NER
    call init_mmwly_para(WOTPE1s0_2000_para, WOTPE3s1_2000_para, WOTPE3p0_2000_set2_para,        &
      & Cs1p1_2000_set2_para, Cs3p1_2000_set2_para, Cs3p2_2000_set2_para, Cs1d2_2000_set2_para  &
      & , Cs3d2_2000_set2_para, Cs3d3_2000_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_2000

  subroutine chengdu_WOTPE_2200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2200.0_NER
    call init_mmwly_para(WOTPE1s0_2200_para, WOTPE3s1_2200_para, WOTPE3p0_2200_set2_para,        &
      & Cs1p1_2200_set2_para, Cs3p1_2200_set2_para, Cs3p2_2200_set2_para, Cs1d2_2200_set2_para  &
      & , Cs3d2_2200_set2_para, Cs3d3_2200_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_2200

  subroutine chengdu_WOTPE_2400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2400.0_NER
    call init_mmwly_para(WOTPE1s0_2400_para, WOTPE3s1_2400_para, WOTPE3p0_2400_set2_para,        &
      & Cs1p1_2400_set2_para, Cs3p1_2400_set2_para, Cs3p2_2400_set2_para, Cs1d2_2400_set2_para  &
      & , Cs3d2_2400_set2_para, Cs3d3_2400_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_2400

  subroutine chengdu_WOTPE_2800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 2800.0_NER
    call init_mmwly_para(WOTPE1s0_2800_para, WOTPE3s1_2800_para, WOTPE3p0_2800_set2_para,        &
      & Cs1p1_2800_set2_para, Cs3p1_2800_set2_para, Cs3p2_2800_set2_para, Cs1d2_2800_set2_para  &
      & , Cs3d2_2800_set2_para, Cs3d3_2800_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_2800

  subroutine chengdu_WOTPE_3200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 3200.0_NER
    call init_mmwly_para(WOTPE1s0_3200_para, WOTPE3s1_3200_para, WOTPE3p0_3200_set2_para,        &
      & Cs1p1_3200_set2_para, Cs3p1_3200_set2_para, Cs3p2_3200_set2_para, Cs1d2_3200_set2_para  &
      & , Cs3d2_3200_set2_para, Cs3d3_3200_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_3200

  subroutine chengdu_WOTPE_4000(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 4000.0_NER
    call init_mmwly_para(WOTPE1s0_4000_para, WOTPE3s1_4000_para, WOTPE3p0_4000_set2_para,        &
      & Cs1p1_4000_set2_para, Cs3p1_4000_set2_para, Cs3p2_4000_set2_para, Cs1d2_4000_set2_para  &
      & , Cs3d2_4000_set2_para, Cs3d3_4000_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_4000

  subroutine chengdu_WOTPE_4800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 4800.0_NER
    call init_mmwly_para(WOTPE1s0_4800_para, WOTPE3s1_4800_para, WOTPE3p0_4800_set2_para,        &
      & Cs1p1_4800_set2_para, Cs3p1_4800_set2_para, Cs3p2_4800_set2_para, Cs1d2_4800_set2_para  &
      & , Cs3d2_4800_set2_para, Cs3d3_4800_set2_para)
    call get_chengdu_WOTPE(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_WOTPE_4800

!===========================!
! MWPLL routines
!===========================!

  subroutine chengdu_MWPLL_400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)       :: L, S, J, uptoQn
    real(NER)                 :: Mambda
    real(NEC), intent(in)  :: po, pi
    real(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol)  :: lsj
    integer       :: chnid
    real(NEC)  :: tmppotval(1:2, 1:2)
    real(NER) :: para0(1:1)

    para0 = 0.0_NER

    Mambda = 400.0_NER
    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_MWPLL : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
        lsj%L = lsj%J - 1
        call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_sepble_smpl_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
                & MWPLL1s0_400_para, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1_400_set2_para, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0_400_set2_para, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VNLO_sepble_smpl_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & MWPLL1s0_400_para, po, pi, tmppotval)
          case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
              & CHNID_3D2_PP, CHNID_1D2_PP)
            call OPE_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
              & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      ! case(2)
      !   select case(chnid)
      !     case (CHNID_1S0_PP)
      !       call VN2LO_sepble_smpl_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
      !          & MWPLL1s0_400_para, po, pi, tmppotval)
!           case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
!               & CHNID_3D2_PP, CHNID_1D2_PP)
        ! end select
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_MWPLL_400

!===========================!
! Complex MMWLY routines
!===========================!

  subroutine chengdu_MMWLY_cmplx_350(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)       :: L, S, J, uptoQn
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER)     :: Mambda
    complex(NEC)  :: tmppotval(1:2, 1:2)

    Mambda = 350.0_NER
    call init_mmwly_para(MMWLY1s0_350_para, Cs3s1_350_set2_para,               &
      & Cs3p0_350_set2_para, Cs1p1_350_set2_para, Cs3p1_350_set2_para,         &
      & Cs3p2_350_set2_para, Cs1d2_350_set2_para, Cs3d2_350_set2_para,         &
      & Cs3d3_350_set2_para)
    call get_chengdu_MMWLY_cmplx(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_cmplx_350

  subroutine chengdu_MMWLY_cmplx_400(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER)     :: Mambda
    complex(NEC)  :: tmppotval(1:2, 1:2)

    Mambda = 400.0_NER
    call init_mmwly_para(MMWLY1s0_400_para, Cs3s1_400_set2_para, Cs3p0_400_set2_para,       &
      & Cs1p1_400_set2_para, Cs3p1_400_set2_para, Cs3p2_400_set2_para, Cs1d2_400_set2_para  &
      & , Cs3d2_400_set2_para, Cs3d3_400_set2_para)
    call get_chengdu_MMWLY_cmplx(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_cmplx_400

  subroutine chengdu_MMWLY_cmplx_600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER)     :: Mambda
    complex(NEC)  :: tmppotval(1:2, 1:2)

    Mambda = 600.0_NER
    call init_mmwly_para(MMWLY1s0_600_para, Cs3s1_600_set2_para, Cs3p0_600_set2_para,       &
      & Cs1p1_600_set2_para, Cs3p1_600_set2_para, Cs3p2_600_set2_para, Cs1d2_600_set2_para  &
      & , Cs3d2_600_set2_para, Cs3d3_600_set2_para)
    call get_chengdu_MMWLY_cmplx(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_cmplx_600

  subroutine chengdu_MMWLY_cmplx_800(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER)     :: Mambda
    complex(NEC)  :: tmppotval(1:2, 1:2)

    Mambda = 800.0_NER
    call init_mmwly_para(MMWLY1s0_800_para, Cs3s1_800_set2_para, Cs3p0_800_set2_para,       &
      & Cs1p1_800_set2_para, Cs3p1_800_set2_para, Cs3p2_800_set2_para, Cs1d2_800_set2_para  &
      & , Cs3d2_800_set2_para, Cs3d3_800_set2_para)
    call get_chengdu_MMWLY_cmplx(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_cmplx_800

  subroutine chengdu_MMWLY_cmplx_1600(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER)     :: Mambda
    complex(NEC)  :: tmppotval(1:2, 1:2)

    Mambda = 1600.0_NER
    call init_mmwly_para(MMWLY1s0_1600_para, Cs3s1_1600_set2_para, Cs3p0_1600_set2_para,        &
      & Cs1p1_1600_set2_para, Cs3p1_1600_set2_para, Cs3p2_1600_set2_para, Cs1d2_1600_set2_para  &
      & , Cs3d2_1600_set2_para, Cs3d3_1600_set2_para)
    call get_chengdu_MMWLY_cmplx(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_cmplx_1600

  subroutine chengdu_MMWLY_cmplx_3200(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    complex(NEC), intent(in)  :: po, pi
    complex(NEC), intent(out) :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER)     :: Mambda
    complex(NEC)  :: tmppotval(1:2, 1:2)

    Mambda = 3200.0_NER
    call init_mmwly_para(MMWLY1s0_3200_para, Cs3s1_3200_set2_para, Cs3p0_3200_set2_para,        &
      & Cs1p1_3200_set2_para, Cs3p1_3200_set2_para, Cs3p2_3200_set2_para, Cs1d2_3200_set2_para  &
      & , Cs3d2_3200_set2_para, Cs3d3_3200_set2_para)
    call get_chengdu_MMWLY_cmplx(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_cmplx_3200

!   pray these subs for getting phase shifts with the chengdu potential

  subroutine LO_chengdu_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call chengdu_hotpot(L, S, J, 0, p1, p2, potval)
!     debug
!   print *, 'into LO_chengdu_epwrap'
!   print *, 'LSJ=', L, S, J, regtype, Mambda, Cs
  end subroutine LO_chengdu_epwrap

  subroutine NLO_chengdu_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call chengdu_hotpot(L, S, J, 1, p1, p2, potval)
  end subroutine NLO_chengdu_epwrap

  subroutine N2LO_chengdu_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)

    integer, intent(in)     :: L, S, j, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, Cs(:)
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    potval = 0.0_NER
    call chengdu_hotpot(L, S, J, 2, p1, p2, potval)
  end subroutine N2LO_chengdu_epwrap

!=======================================
! Pots for diagnosing power counting
!=======================================

  subroutine get_chengdu_MMWLY_diagnose(L, S, J, uptoQn, po, pi, Mambda, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi, Mambda
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)
    real(NER) :: para0(1:1)

    para0 = 0.0_NER
    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_MMWLY_diagnose : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
      lsj%L = lsj%J - 1
      call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VNLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          ! case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
          !     & CHNID_3D2_PP, CHNID_1D2_PP)
          !   call OPE_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
          !     & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(2)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VN2LO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          ! case (CHNID_3P0_PP)
          !   call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs3p0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          ! case (CHNID_3P1_PP)
          !   call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs3p1, po, pi, tmppotval)
          ! case (CHNID_1P1_PP)
          !   call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs1p1, po, pi, tmppotval)
          ! case (CHNID_3P2_PP)
          !   call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs3p2, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine get_chengdu_MMWLY_diagnose

  subroutine chengdu_MMWLY_600_diagnose(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 600.0_NER
    call init_mmwly_para(MMWLY1s0_600_para, Cs3s1_600_set2_para, Cs3p0_600_set2_para,       &
      & Cs1p1_600_set2_para, Cs3p1_600_set2_para, Cs3p2_600_set2_para, Cs1d2_600_set2_para  &
      & , Cs3d2_600_set2_para, Cs3d3_600_set2_para)
    call get_chengdu_MMWLY_diagnose(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_600_diagnose

  subroutine get_chengdu_MMWLY_only_1S0_HO(L, S, J, uptoQn, po, pi, Mambda, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi, Mambda
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)
    real(NER) :: para0(1:1)

    para0 = 0.0_NER
    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'get_chengdu_MMWLY_diagnose : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_CPLD .and. L > J) then
      lsj%L = lsj%J - 1
      call convert_lsj_to_chnid(lsj, chnid)
    end if

    select case(uptoQn)
      case(0)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          case (CHNID_3S1_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3s1, po, pi, tmppotval)
          case (CHNID_3P0_PP)
            call VLO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
              & Cs3p0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(1)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VNLO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          ! case (CHNID_3P1_PP, CHNID_1P1_PP, CHNID_3D3_PP, CHNID_3P2_PP, &
          !     & CHNID_3D2_PP, CHNID_1D2_PP)
          !   call OPE_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,         &
          !     & para0, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case(2)
        select case(chnid)
          case (CHNID_1S0_PP)
            call VN2LO_MMWLY_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,   &
              & mmwly1s0, po, pi, tmppotval)
          ! case (CHNID_3P0_PP)
          !   call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs3p0, po, pi, tmppotval)
          ! case (CHNID_3S1_PP)
            ! call VN2LO_withc0_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
            !   & Cs3s1, po, pi, tmppotval)
          ! case (CHNID_3P1_PP)
          !   call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs3p1, po, pi, tmppotval)
          ! case (CHNID_1P1_PP)
          !   call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs1p1, po, pi, tmppotval)
          ! case (CHNID_3P2_PP)
          !   call VN2LO_smplper_epwrap(lsj%L, lsj%S, lsj%J, CHENGDU_REGTYPE, Mambda,  &
          !     & Cs3p2, po, pi, tmppotval)
          case default
            tmppotval = 0.0_NER
        end select
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine get_chengdu_MMWLY_only_1S0_HO

  subroutine chengdu_MMWLY_600_only_1S0_HO(L, S, J, uptoQn, po, pi, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: po, pi
    real(NER), intent(out)  :: potval(1:NCH, 1:NCH)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2), Mambda

    Mambda = 600.0_NER
    call init_mmwly_para(MMWLY1s0_600_para, Cs3s1_600_set2_para, Cs3p0_600_set2_para,       &
      & Cs1p1_600_set2_para, Cs3p1_600_set2_para, Cs3p2_600_set2_para, Cs1d2_600_set2_para  &
      & , Cs3d2_600_set2_para, Cs3d3_600_set2_para)
    call get_chengdu_MMWLY_only_1S0_HO(L, S, J, uptoQn, po, pi, Mambda, potval)

  end subroutine chengdu_MMWLY_600_only_1S0_HO

end module chengdu

