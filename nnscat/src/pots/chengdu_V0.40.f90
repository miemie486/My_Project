! Bingwei Long 08/03/2019
!   - Ver 0.4
!   - Added NLO for chengdu_DLSPR*
! Bingwei Long 08/03/2019
!   - Ver 0.32
!   - Added show_mtrx.f90 as example
!   - Added chengdu.README
! Bingwei Long 07/25/2019
!   - Ver 0.30
!   - Used Shaowei's contacts for pert. P waves
! Bingwei Long 07/16/2019
!   - Added chengdu_DLSPR_generic in refactoring
! Bingwei Long  07/02/2019
!   - Added 800, 1600, 3200 MeV LO values
! Bingwei Long  03/18/2019
!   - Improved comments
!   - Modified to accommodate changes made to interface of *_epwrap()
! Bingwei Long  01/28/2019

module chengdu

  use potwrap
  implicit none

  private

  ! Routine names for potential matrix elements. Find explanation in comments
  ! above their code.
  public  :: chengdu_DLSPR_generic, chengdu_DLSPR, chengdu_DLSPR_800,          &
    & chengdu_DLSPR_1600, chengdu_DLSPR_3200

  integer, parameter ::  CHENGDU_REGTYPE = REGTYPE_GAUSSIAN
  integer, parameter ::  MAX_NPARA = 30, MAX_NFMLY = 10

  integer, parameter ::           &
    & NPARA_N3LO_spr1S0 = 14,     &
    & NPARA_N3LO_3S1 = 7,         &
    & NPARA_N3LO_3P0 = 5,         &
    & NPARA_N3LO_PWAVES = 2

  integer, parameter    ::  NFMLY_dltless = 4

  real(NER), parameter  ::  Mmbd_dltless(1:NFMLY_dltless) = (/              &
    & 600.0_NER,                                                            &
    & 800.0_NER,                                                            &
    & 1600.0_NER,                                                           &
    & 3200.0_NER /)

  integer, parameter  :: cis_setID_dltless(1:NFMLY_dltless) =  (/           &
    & 2,                                                                    &
    & 2,                                                                    &
    & 2,                                                                    &
    & 2 /)

  real(NER), parameter :: phony_paras(1:30) = 0.0_NER

  real(NER), dimension(1:NFMLY_dltless, 1:NPARA_N3LO_spr1S0), parameter     &
    & :: spr1S0_dltless = &
    & reshape( (/         &
    & 0.28476583881316405E+002_NER,   -0.94720748744753225E-001_NER,      &
    & 0.20650483002051136E+002_NER,   -0.97698777628877326E-001_NER,      &
    & 0.95630015528180693E-003_NER,   -0.11694692770857738E+002_NER,      &
    & 0.86955645956220531E-001_NER,   -0.27572312790623797E-002_NER,      &
    & 0.51272513824358814E-007_NER,   -0.15365873931283828E+002_NER,      &
    & 0.24216603991719265E+000_NER,   -0.25138848487176349E-002_NER,      &
    & 0.55453661564745701E-007_NER,    0.61631676811145938E-012_NER,      &
    ! 600.0 ends here
    & 0.26900213333936705E+002_NER,   -0.85152284790577140E-001_NER,      &
    & 0.20815050385250149E+002_NER,   -0.92906542271354697E-001_NER,      &
    & 0.85676399438280849E-003_NER,    0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    ! 800.0 ends here
    & 0.25741081111240778E+002_NER,   -0.75410791126733698E-001_NER,      &
    & 0.19005783674206164E+002_NER,   -0.76594455043054649E-001_NER,      &
    & 0.58079932223686253E-003_NER,    0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    ! 1600 ends here
    & 0.25530308969449205E+002_NER,   -0.72062514706775244E-001_NER,      &
    & 0.21227509375846985E+002_NER,   -0.93179562856372411E-001_NER,      &
    & 0.55183479526630817E-003_NER,    0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER,                           &
    & 0.0_NER,                         0.0_NER                            &
    ! 3200 ends here
    & /), shape(spr1S0_dltless))

  real(NER), dimension(1:NFMLY_dltless, 1:NPARA_N3LO_3S1), parameter      &
    & :: Cs3S1_dltless = reshape( (/                                      &
    & -0.12458041537666564E-002_NER, 0.0_NER, 0.0_NER, 0.0_NER,           &
    & 0.0_NER, 0.0_NER, 0.0_NER,                                          &
    ! 600 ends here
    & 0.31185721559687796E-002_NER, 0.0_NER, 0.0_NER, 0.0_NER,            &
    & 0.0_NER, 0.0_NER, 0.0_NER,                                          &
    ! 800 ends here
    & -0.45770398757953171E-002_NER, 0.0_NER, 0.0_NER, 0.0_NER,           &
    & 0.0_NER, 0.0_NER, 0.0_NER,                                          &
    ! 1600 ends here
    & -0.12038131670927033E-001_NER, 0.0_NER, 0.0_NER, 0.0_NER,           &
    & 0.0_NER, 0.0_NER, 0.0_NER                                           &
    ! 3200 ends here
    & /), shape(Cs3S1_dltless))

  real(NER), dimension(1:NFMLY_dltless, 1:NPARA_N3LO_3P0), parameter      &
    & :: Cs3P0_dltless =    &
    & reshape( (/           &
    & 0.47499661698324306E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER,  &
    ! 600 ends here
    & -0.39852401392998397E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, &
    ! 800 ends here
    & 0.32109431145360149E-009_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER,  &
    ! 1600 ends here
    & 0.11866368639056420E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER   &
    ! 3200 ends here
    & /), shape(Cs3P0_dltless))

  real(NER), dimension(1:NFMLY_dltless, 1:NPARA_N3LO_PWAVES), parameter        &
    & :: Cs1P1_dltless =    &
    & reshape( (/           &
    & 4.2339204028700879E-009, 8.0038133082506403E-009,                        &
    ! 600 ends here
    & 3.9887140152981041E-009, 8.3462807917658315E-009,                        &
    ! 800 ends here
    & 3.9076834859914674E-009, 8.9475674238727265E-009,                        &
    ! 1600 ends here
    & 3.9043472269067605E-009, 9.2636353487866599E-009                         &
    ! 3200 ends here
    & /), shape(Cs1P1_dltless))

  real(NER), dimension(1:NFMLY_dltless, 1:NPARA_N3LO_PWAVES), parameter        &
    & :: Cs3P1_dltless =    &
    & reshape( (/           &
    & -1.7477045629225740E-010, -1.6746565542979603E-008,                      &
    ! 600 ends here
    & 1.0990360062986895E-009, -1.6462595788485315E-008,                       &
    ! 800 ends here
    & 3.2868414934285669E-009, -1.0859556344724893E-008,                       &
    ! 1600 ends here
    & 4.3912359225933715E-009, 5.2651340671094915E-009                         &
    ! 3200 ends here
    & /), shape(Cs3P1_dltless))

  real(NER), dimension(1:NFMLY_dltless, 1:NPARA_N3LO_PWAVES), parameter        &
    & :: Cs3P2_dltless =    &
    & reshape( (/           &
    & -4.0603072057458670E-009, 8.2678556817969446E-010,                      &
    ! 600 ends here
    & -2.8326731914154320E-009, 1.1228037831553700E-009,                      &
    ! 800 ends here
    & -1.2767286153515018E-009, 1.7739395934620311E-009,                      &
    ! 1600 ends here
    & -3.3841229410501060E-010, 2.1816438465765831E-009                        &
    ! 3200 ends here
    & /), shape(Cs3P2_dltless))

contains

  ! %%%%%%
  ! In
  ! %%%%%%
  ! uptoQn: number of order
  ! uptoQn = 0 -> LO
  ! uptoQn = 1 -> NLO, smaller than LO by O(Q/Mhi), where Mhi is the breakdown
  ! scale of chiral EFT.
  ! And so on...
  !
  ! L, S, J: quantum # of partial-wave, as in ^{2S+1}L_J
  !   Non-vanishing channel at LO: 1S0(LSJ = 000), 3S1 - 3D1(LSJ = 011)
  !     , 3P0(LSJ = 110)
  !   For coupled channels, L = lower one of two coupled angular momenta
  ! p1, p2: magnitude of outgoing and incoming c.m. momenta
  !
  ! %%%%%%
  ! Out
  ! %%%%%%
  ! potval is the value of partial-wave potentials
  ! Storage of potval:
  !   2-by-2 matrix
  !   For uncoupled channels, e.g., 1S0(LSJ = 000) or 3P0(LSJ = 110)
  !     potval(1, 1) = <L, p1 | V | L, p2 > and (1, 2), (2, 1), and (2, 2) are
  !     zero
  !   For coupled channels, e.g., 3S1 - 3D1 (LSJ = 011)
  !     potval(1, 1) = < L, p1   | V | L, p2 >
  !     potval(1, 2) = < L, p1   | V | L+2, p2 >
  !     potval(2, 1) = < L+2, p1 | V | L, p2 >
  !     potval(2, 2) = < L+2, p1 | V | L+2, p2 >
  !
  ! Remarks:
  ! * All units are MeV or (MeV)^{-1}
  !
  ! * Normalization
  !   If V were so weak as to validate the Born approximation, the relation of
  !   V_1S0 and the scattering length 'a' would have been
  !   < L, p | V | L, p> = 2/pi * a + ... for p -> 0
  !   This may be larger than V in other conventions by a factor of mN.
  !
  ! * Sign
  ! For attractive V, sign of V is minus

  ! Force family: Deltaless, separable 1S0
  ! forceID = 1 - 4 : individual index of force from this family
  ! * 1: cutoff = 600 MeV
  ! * 2: cutoff = 800 MeV
  ! * 3: cutoff = 1600 MeV
  ! * 4: cutoff = 3200 MeV
  subroutine chengdu_DLSPR_generic(forceID, L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: forceID, L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)

    lsj%L = L
    lsj%S = S
    lsj%J = J
    tmppotval = 0.0_NER

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_DLSPR : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    call Select_cis2PE_vtpe1(cis_setID_dltless(forceID))

    if (uptoQn == 0) then
      select case(chnid)
        case (CHNID_1S0_PP)
          call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE,                      &
          & Mmbd_dltless(forceID),                                             &
          & spr1S0_dltless(forceID, 1:NPARA_N3LO_spr1S0), p1, p2, tmppotval)
        case (CHNID_3S1_PP)
          call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE,                     &
          & Mmbd_dltless(forceID),                                             &
          & Cs3S1_dltless(forceID, 1:NPARA_N3LO_3s1), p1, p2, tmppotval)
        case (CHNID_3P0_PP)
          call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE,                     &
          & Mmbd_dltless(forceID),                                             &
          & Cs3P0_dltless(forceID, 1:NPARA_N3LO_3p0), p1, p2, tmppotval)
        case default
          tmppotval = 0.0_NER
      end select
    end if

    if (uptoQn == 1) then
      select case(chnid)
        case (CHNID_1S0_PP)
          call VNLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE,                     &
          & Mmbd_dltless(forceID),                                             &
          & spr1S0_dltless(forceID, 1:NPARA_N3LO_spr1S0), p1, p2, tmppotval)
        case (CHNID_3S1_PP, CHNID_3P0_PP)
          tmppotval = 0.0_NER
        case default
          call OPE_epwrap(L, S, J, CHENGDU_REGTYPE,                            &
          & Mmbd_dltless(forceID), phony_paras, p1, p2, tmppotval)
      end select
    end if

    if (uptoQn == 2) then
      select case(chnid)
        case (CHNID_1S0_PP)
          call VN2LO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE,                    &
          & Mmbd_dltless(forceID),                                             &
          & spr1S0_dltless(forceID, 1:NPARA_N3LO_spr1S0), p1, p2, tmppotval)
        case (CHNID_3S1_PP)
          call VN2LO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE,                   &
          & Mmbd_dltless(forceID),                                             &
          & Cs3S1_dltless(forceID, 1:NPARA_N3LO_3s1), p1, p2, tmppotval)
        case (CHNID_3P0_PP)
          call VN2LO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE,                   &
          & Mmbd_dltless(forceID),                                             &
          & Cs3P0_dltless(forceID, 1:NPARA_N3LO_3p0), p1, p2, tmppotval)
        case (CHNID_1P1_PP)
          call VN2LO_smplper_epwrap(L, S, J, CHENGDU_REGTYPE,                  &
          & Mmbd_dltless(forceID),                                             &
          & Cs1P1_dltless(forceID, 1:NPARA_N3LO_PWAVES), p1, p2, tmppotval)
        case (CHNID_3P1_PP)
          call VN2LO_smplper_epwrap(L, S, J, CHENGDU_REGTYPE,                  &
          & Mmbd_dltless(forceID),                                             &
          & Cs3P1_dltless(forceID, 1:NPARA_N3LO_PWAVES), p1, p2, tmppotval)
        case (CHNID_3P2_PP)
          call VN2LO_smplper_epwrap(L, S, J, CHENGDU_REGTYPE,                  &
          & Mmbd_dltless(forceID),                                             &
          & Cs3P2_dltless(forceID, 1:NPARA_N3LO_PWAVES), p1, p2, tmppotval)
        case default
          tmppotval = 0.0_NER
      end select
    end if

    if (uptoQn == 3) then
      select case(chnid)
        case (CHNID_1S0_PP)
          call VN3LO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE,                    &
          & Mmbd_dltless(forceID),                                             &
          & spr1S0_dltless(forceID, 1:NPARA_N3LO_spr1S0), p1, p2, tmppotval)
        case (CHNID_3S1_PP)
          call VN3LO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE,                   &
          & Mmbd_dltless(forceID),                                             &
          & Cs3S1_dltless(forceID, 1:NPARA_N3LO_3s1), p1, p2, tmppotval)
        case (CHNID_3P0_PP)
          call VN3LO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE,                   &
          & Mmbd_dltless(forceID),                                             &
          & Cs3P0_dltless(forceID, 1:NPARA_N3LO_3p0), p1, p2, tmppotval)
        case (CHNID_1P1_PP)
          call VN3LO_pope_smpl_epwrap(L, S, J, CHENGDU_REGTYPE,                &
          & Mmbd_dltless(forceID),                                             &
          & Cs1P1_dltless(forceID, 1:NPARA_N3LO_PWAVES), p1, p2, tmppotval)
        case (CHNID_3P1_PP)
          call VN3LO_pope_smpl_epwrap(L, S, J, CHENGDU_REGTYPE,                &
          & Mmbd_dltless(forceID),                                             &
          & Cs3P1_dltless(forceID, 1:NPARA_N3LO_PWAVES), p1, p2, tmppotval)
        case (CHNID_3P2_PP)
          call VN3LO_pope_smpl_epwrap(L, S, J, CHENGDU_REGTYPE,                &
          & Mmbd_dltless(forceID),                                             &
          & Cs3P2_dltless(forceID, 1:NPARA_N3LO_PWAVES), p1, p2, tmppotval)
        case default
          call TPE0_epwrap(L, S, J, CHENGDU_REGTYPE,                           &
          & Mmbd_dltless(forceID), phony_paras, p1, p2, tmppotval)
          tmppotval = 0.0_NER
      end select
    end if

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_DLSPR_generic


! Legacy subroutines

  ! forceID = 1, cutoff = 600 MeV, 1S0 separable, deltaless
  subroutine chengdu_DLSPR(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    call chengdu_DLSPR_generic(1, L, S, J, uptoQn, p1, p2, potval)

  end subroutine chengdu_DLSPR

  ! forceID = 2, cutoff = 800 MeV, 1S0 separable, deltaless
  subroutine chengdu_DLSPR_800(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    call chengdu_DLSPR_generic(2, L, S, J, uptoQn, p1, p2, potval)

  end subroutine chengdu_DLSPR_800

  ! forceID = 3, cutoff = 1600 MeV, 1S0 separable, deltaless
  subroutine chengdu_DLSPR_1600(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    call chengdu_DLSPR_generic(3, L, S, J, uptoQn, p1, p2, potval)

  end subroutine chengdu_DLSPR_1600

  ! forceID = 4, cutoff = 3200 MeV, 1S0 separable, deltaless
  subroutine chengdu_DLSPR_3200(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    call chengdu_DLSPR_generic(4, L, S, J, uptoQn, p1, p2, potval)

  end subroutine chengdu_DLSPR_3200

end module chengdu
