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

module chengdu

  ! Version 0.17

  use potwrap
  use nopion_epmod
  implicit none

  private
  public  :: chengdu_DLSPR_350, chengdu_DLSPR_400, chengdu_DLSPR,   &
    & chengdu_DLSPR_800, chengdu_DLSPR_1600, chengdu_DLSPR_3200,    &
    & chengdu_pionless_1000

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
    & -1.4410457830873011E-003_NER, 0.0_NER, 0.0_NER, 0.0_NER,              &
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
    & 2.9151500814276572E-003_NER, 0.0_NER, 0.0_NER, 0.0_NER,               &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd2_set2_para(1:5) =  (/                                      &
    & -0.39852401392998397E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter, private  ::  &
    & Mambda_3 = 1600.0_NER,                                                &
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
    & -4.6744212081189567E-003_NER, 0.0_NER, 0.0_NER, 0.0_NER,              &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd3_set2_para(1:5) =  (/                                      &
    & 0.32109431145360149E-009_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  real(NER), parameter, private  ::  &
    & Mambda_4 = 3200.0_NER,                                                &
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
    & -1.2736819542172845E-002_NER, 0.0_NER, 0.0_NER, 0.0_NER,              &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd4_set2_para(1:5) =  (/                                      &
    & 0.11866368639056420E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

    real(NER), parameter, private  ::  &
    & Mambda_5 = 350.0_NER,                                                 &
    &                                                                       &
    & spr1s0_Mmbd5_set2_para(1:14) = (/                                     &
    &   0.42647716568097287E+002_NER,   -0.15817891988985733E+000_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_Mmbd5_set2_para(1:7) =  (/                                      &
    & -5.4945264707960580E-003_NER, 0.0_NER, 0.0_NER, 0.0_NER,              &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd5_set2_para(1:5) =  (/                                      &
    & 0.28737080761182509E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

    real(NER), parameter, private  ::  &
    & Mambda_6 = 400.0_NER,                                                 &
    &                                                                       &
    & spr1s0_Mmbd6_set2_para(1:14) = (/                                     &
    &   0.39275245500005816E+002_NER,   -0.13862263192892474E+000_NER,      &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER,                           &
    &   0.0_NER,                         0.0_NER/),                         &
    &                                                                       &
    & Cs3s1_Mmbd6_set2_para(1:7) =  (/                                      &
    & -4.5126679738286763E-003_NER, 0.0_NER, 0.0_NER, 0.0_NER,              &
    & 0.0_NER, 0.0_NER, 0.0_NER/),                                          &
    &                                                                       &
    & Cs3p0_Mmbd6_set2_para(1:5) =  (/                                      &
    & 0.30164174135677295E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

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
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_1,   &
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
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_2,   &
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
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_3,   &
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
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_4,   &
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

  subroutine chengdu_DLSPR_350(L, S, J, uptoQn, p1, p2, potval)

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
        & 'chengdu_DLSPR_350 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_5,   &
          & spr1s0_Mmbd5_set2_para, p1, p2, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_5,  &
          & Cs3s1_Mmbd5_set2_para, p1, p2, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_5,  &
          & Cs3p0_Mmbd5_set2_para, p1, p2, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_DLSPR_350

  subroutine chengdu_DLSPR_400(L, S, J, uptoQn, p1, p2, potval)

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
        & 'chengdu_DLSPR_400 : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)

    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_sprbl_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_6,   &
          & spr1s0_Mmbd6_set2_para, p1, p2, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_6,  &
          & Cs3s1_Mmbd6_set2_para, p1, p2, tmppotval)
      case (CHNID_3P0_PP)
        call VLO_withc0_epwrap(L, S, J, CHENGDU_REGTYPE, Mambda_6,  &
          & Cs3p0_Mmbd6_set2_para, p1, p2, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_DLSPR_400

  subroutine chengdu_pionless_1000(L, S, J, uptoQn, p1, p2, potval)

    integer, intent(in)     :: L, S, J, uptoQn
    real(NER), intent(in)   :: p1, p2
    real(NER), intent(out)  :: potval(:, :)

    type(lsj_symbol) :: lsj
    integer :: chnid
    real(NER) :: tmppotval(1:2, 1:2)
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

    pionlesspara1s0(1) =  -1.56315E-003_NER   ! Fitted to a = -23.7 fm
    pionlesspara3s1(1) =  -1.7862E-003_NER    ! Fitted to Bd = 2.22 MeV

    select case(chnid)
      case (CHNID_1S0_PP)
        call VLO_withc0_pionless_NP_epwrap(L, S, J, CHENGDU_REGTYPE, 1000.0_NER,  &
          & pionlesspara1s0, p1, p2, tmppotval)
      case (CHNID_3S1_PP)
        call VLO_withc0_pionless_NP_epwrap(L, S, J, CHENGDU_REGTYPE, 1000.0_NER,  &
          & pionlesspara3s1, p1, p2, tmppotval)
      case default
        tmppotval = 0.0_NER
    end select

    potval(1:2, 1:2) = tmppotval(1:2, 1:2)

  end subroutine chengdu_pionless_1000

end module chengdu
