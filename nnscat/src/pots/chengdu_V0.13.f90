! Bingwei Long  10/04/2019
!   - Added init_PCs4chengdu
! Bingwei Long  07/02/2019
!   - Added 800, 1600, 3200 MeV LO values
! Bingwei Long  03/18/2019
!   - Improved comments
!   - Modified to accommodate changes made to interface of *_epwrap()
! Bingwei Long  01/28/2019

module chengdu

  ! Ver. 0.13

  use potwrap
  implicit none

  private
  public  :: chengdu_DLSPR, chengdu_DLSPR_800, chengdu_DLSPR_1600,    &
    & chengdu_DLSPR_3200

  integer, parameter, private ::  CHENGDU_REGTYPE = REGTYPE_GAUSSIAN

  real(NER), parameter, private  ::  &
    & Mambda_1 = 600.0_NER,                                                 &
    &                                                                       &
    & spr1s0_Mmbd1_set2_para(1:14) = (/                                     &
    &   0.28476583881316405E+002_NER,   -0.94720748744753225E-001_NER,      &
    &   0.20650483002051136E+002_NER,   -0.97698777628877326E-001_NER,      &
    &   0.95630015528180693E-003_NER,   -0.11694692770857738E+002_NER,      &
    &   0.86955645956220531E-001_NER,   -0.27572312790623797E-002_NER,      &
    &   0.51272513824358814E-007_NER,   -0.15365873931283828E+002_NER,      &
    &   0.24216603991719265E+000_NER,   -0.25138848487176349E-002_NER,      &
    &   0.55453661564745701E-007_NER,    0.61631676811145938E-012_NER/),    &
    &                                                                       &
    & Cs3s1_Mmbd1_set2_para(1:7) =  (/                                      &
    & -0.12458041537666564E-002_NER, 0.0_NER, 0.0_NER, 0.0_NER,             &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd1_set2_para(1:5) =  (/                                      &
    & 0.47499661698324306E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter, private  ::  &
    & Mambda_2 = 800.0_NER,                                                 &
    &                                                                       &
    & spr1s0_Mmbd2_set2_para(1:14) = (/                                     &
    &   0.26900213333936705E+002_NER,   -0.85152284790577140E-001_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_Mmbd2_set2_para(1:7) =  (/                                      &
    & 0.31185721559687796E-002_NER, 0.0_NER, 0.0_NER, 0.0_NER,             &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd2_set2_para(1:5) =  (/                                      &
    & -0.39852401392998397E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter, private  ::  &
    & Mambda_3 = 1600.0_NER,                                                 &
    &                                                                       &
    & spr1s0_Mmbd3_set2_para(1:14) = (/                                     &
    &   0.25741081111240778E+002_NER,   -0.75410791126733698E-001_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_Mmbd3_set2_para(1:7) =  (/                                      &
    & -0.45770398757953171E-002_NER, 0.0_NER, 0.0_NER, 0.0_NER,             &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd3_set2_para(1:5) =  (/                                      &
    & 0.32109431145360149E-009_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter, private  ::  &
    & Mambda_4 = 3200.0_NER,                                                 &
    &                                                                       &
    & spr1s0_Mmbd4_set2_para(1:14) = (/                                     &
    &   0.25530308969449205E+002_NER,   -0.72062514706775244E-001_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_Mmbd4_set2_para(1:7) =  (/                                      &
    & -0.12038131670927033E-001_NER, 0.0_NER, 0.0_NER, 0.0_NER,             &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd4_set2_para(1:5) =  (/                                      &
    & 0.11866368639056420E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

contains

  subroutine init_PCs4chengdu()

    PC_gA = 1.29_NER
    PC_mN = 939.0_NER
    PC_mpi = 138.0_NER
    PC_fpi = 92.4_NER

  end subroutine init_PCs4chengdu

  ! Feautures: Deltaless; separable 1s0
  !
  ! All units are in MeV or (MeV)^{-1}
  ! %%%%%%
  ! In
  ! %%%%%%
  ! uptoQn: number of order
  !   As of 01/29/2019, uptoQn = 0 (LO)
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
  ! * Normalization of V can be explained as follows.
  ! If V were so weak as to validate the Born approximation, the relation of
  ! V in uncoupled channels and the respective scattering length 'a' would have
  ! been
  ! < L, p | V | L, p> = 2/pi * a + ... for p -> 0

  subroutine chengdu_DLSPR(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)

    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_DLSPR : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_1,  &
          & spr1s0_Mmbd1_set2_para, p1, p2, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_1,  &
          & Cs3s1_Mmbd1_set2_para, p1, p2, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_1,  &
          & Cs3p0_Mmbd1_set2_para, p1, p2, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_DLSPR

  subroutine chengdu_DLSPR_800(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)

    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_DLSPR_800 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_2,  &
          & spr1s0_Mmbd2_set2_para, p1, p2, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_2,  &
          & Cs3s1_Mmbd2_set2_para, p1, p2, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_2,  &
          & Cs3p0_Mmbd2_set2_para, p1, p2, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_DLSPR_800

  subroutine chengdu_DLSPR_1600(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)

    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_DLSPR_1600 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_3,  &
          & spr1s0_Mmbd3_set2_para, p1, p2, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_3,  &
          & Cs3s1_Mmbd3_set2_para, p1, p2, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_3,  &
          & Cs3p0_Mmbd3_set2_para, p1, p2, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_DLSPR_1600

  subroutine chengdu_DLSPR_3200(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)

    call init_PCs4chengdu()
    lsj%L = L
    lsj%S = S
    lsj%J = J

    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')  &
        & 'chengdu_DLSPR_3200 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_4,  &
          & spr1s0_Mmbd4_set2_para, p1, p2, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_4,  &
          & Cs3s1_Mmbd4_set2_para, p1, p2, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_4,  &
          & Cs3p0_Mmbd4_set2_para, p1, p2, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_DLSPR_3200

end module chengdu
