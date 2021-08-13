module eft_tmtrx

  use util_gauleg,    only: feb_tanabsc
  use potwrap
  use chengdu
  use nneft_lsesv
  use nneft_loops
  use eft_h2wfs
  ! use omp_lib
  implicit none

  integer, parameter  ::  NUM_PARA_TWOBODY = 10

contains

  ! Full Off-shell Sub-Threshold t-matrix
  ! For a list of mE < 0, calculate t(p1, p2; mE) where mE = mN * E
  ! In:
  !   L, S, J as in (2S+1) L_J
  !   Non-vanishing channel: 1S0(LSJ = 000), 3S1 - 3D1(LSJ = 011)
  !     , 3P0(LSJ = 110)
  !   For other channel, t = 0.0
  !   N_mE = number of mE's
  !   mE_lst = array of mE
  !   Mambda = Cutoff value; use 600.0_NER if you don't know anything about it
  !   Nmesh = number of mesh points
  ! Out:
  !   msh = list of mesh points in p
  !   wght = list of weights to be used with msh
  !   tm_lst = list of T-matrix elements
  ! Normalization of t-matrix matches Elster's notes Eq.(1.94)
  !   S = 1 - i\pi k mN t(k, k, mE) and S = e^(i2\delta) for uncoupled
  !   channels, when mE > 0 and k^2 = mE
  ! Storage of tm_lst:
  ! For uncoupled channels, e.g., 1S0(LSJ = 000) or 3P0(LSJ = 110)
  !   T(msh(i), msh(j); mE_lst(k)) = tm_lst(i, j, k)
  ! For coupled channels, e.g., 3S1 - 3D1 (LSJ = 011)
  !   T_{l l'}(p1, p2; mE) is a 2 X 2 matrix where
  !   the value of l(l') can be 0 or 2, hence a 2 X 2 matrix. Or equivalently, T
  !   can be labeled by T_{ab}(p1, p2; mE) where
  !   a(b) = 1 -> l(l') = 0 and a(b) = 2 -> l(l') = 2. Now, back to tm_lst
  !   T_{ab}(msh(i), msh(j); mE_lst(k)) =
  !         tm_lst((a-1)*Nmesh + i, (b-1)*Nmesh + j, k)
  subroutine get_FOST_tmtrx_DLSPR(L, S, J, N_mE, mE_lst, Mambda, Nmesh, msh,   &
    & wght, tm_lst)
    ! Calculate T_matrix. Fixed--size tm_lst is used as input instead of tm_lst(:, :, :) because it's more convenience for passing to C++. Ref: https://stackoverflow.com/questions/29132606/the-run-time-aborting-when-calling-c-sub-from-fortran

    integer, intent(in)     :: L, S, J, N_mE, Nmesh
    real(NER), intent(in)   :: Mambda, mE_lst(:)
    real(NER), intent(out)  :: msh(:), wght(:), tm_lst(:, :, :)

    integer :: ii, jj, ll, mm, chnid, NCH, lsesv_flag
    real(NER) :: qsqr(1:Nmesh), potval(1:2, 1:2)
    real(NER), allocatable :: tmpTM(:, :), tmpVM(:, :)
    type(lsj_symbol) :: lsj

    lsj%L = L
    lsj%S = S
    lsj%J = J
    if (.not. is_lsj_valid(lsj)) then
       write (standard_error_unit, '(a)')   &
        & 'get_FOST_tmtrx_DLSPR : L, S, J invalid'
       return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
       NCH = 1
    else
       NCH = 2
    end if

    allocate(tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh), tmpVM(1:NCH*Nmesh, 1:NCH*Nmesh))
    call feb_tanabsc(Nmesh,                                           &
      & Mambda*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), msh, wght)
    do ii = 1, Nmesh
       qsqr(ii) = msh(ii)*msh(ii)
    end do

    ! For energy-independent pots
    do ii = 1, Nmesh
       do jj = 1, Nmesh
          call chengdu_DLSPR(L, S, J, 0, msh(ii), msh(jj), potval)
          ! call chengdu_DLSPR_350(L, S, J, 0, msh(ii), msh(jj), potval)
          ! call chengdu_DLSPR_400(L, S, J, 0, msh(ii), msh(jj), potval)
          ! call chengdu_DLSPR_800(L, S, J, 0, msh(ii), msh(jj), potval)
          ! call chengdu_DLSPR_1600(L, S, J, 0, msh(ii), msh(jj), potval)
          ! call chengdu_DLSPR_3200(L, S, J, 0, msh(ii), msh(jj), potval)
          call SetLargeVM(tmpVM, NCH, Nmesh, ii, jj, potval)
       end do
    end do

    do ii = 1, N_mE
      call lsesv_nch_full_ST_cpld(NCH, mE_lst(ii), Nmesh, tmpVM, msh, wght,   &
        & qsqr, tmpTM, lsesv_flag)
      tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii) = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
    end do

    ! ! For energy-dependent pots
    ! do ii = 1, N_mE
    !   do ll = 1, Nmesh
    !      do mm = 1, Nmesh
    !         call chengdu_LBWDIB_600(L, S, J, 0, msh(ll), msh(mm), mE_lst(ii)/PC_mN, potval)
    !         call SetLargeVM(tmpVM, NCH, Nmesh, ll, mm, potval)
    !      end do
    !   end do
    !   call lsesv_nch_full_ST_cpld(NCH, mE_lst(ii), Nmesh, tmpVM, msh, wght,   &
    !     & qsqr, tmpTM, lsesv_flag)
    !   tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii) = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
    ! end do

    deallocate(tmpTM, tmpVM)

  end subroutine get_FOST_tmtrx_DLSPR

  ! Tangent mesh to be used
  ! pii = TanMeshC*tan(\pi/4*(xi+1)) where ti are GAUSS-LEGENDRE abscissas
  ! num_pole : number of poles in mE
  ! poltlst : list of poles in mE
  subroutine get_tmtrx_polelst_chengdu(L, S, J, N_mE, mE_lst, TanMeshC, Nmesh, msh,  &
    & wght, pot_name, tm_lst, num_pole, polelst)

    integer, intent(in)     :: L, S, J, N_mE, Nmesh
    integer, intent(out)    :: num_pole
    real(NER), intent(in)   :: TanMeshC, mE_lst(:)
    character(len = *), intent(in)  :: pot_name
    real(NER), intent(out)  :: msh(:), wght(:), tm_lst(:, :, :), polelst(:)

    integer :: ii, jj, ll, mm, chnid, NCH, lsesv_flag
    real(NER) :: qsqr(1:Nmesh), potval(1:2, 1:2)
    real(NER), allocatable :: tmpTM(:, :), tmpVM(:, :)
    type(lsj_symbol) :: lsj

    ! For now, poles are not reported
    num_pole = 0

    lsj%L = L
    lsj%S = S
    lsj%J = J
    if (.not. is_lsj_valid(lsj)) then
       write (standard_error_unit, '(a)')   &
        & 'get_tmtrx_chengdu : L, S, J invalid'
       return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
       NCH = 1
    else
       NCH = 2
    end if

    allocate(tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh), tmpVM(1:NCH*Nmesh, 1:NCH*Nmesh))
    call feb_tanabsc(Nmesh, TanMeshC, msh, wght)
    do ii = 1, Nmesh
      qsqr(ii) = msh(ii)*msh(ii)
    end do

    call chengdu_dispatch(pot_name)
    if (chengdu_hotpot_is_Edep(L, S, J)) then
      ! For energy-dependent pots
      do ii = 1, N_mE
        do ll = 1, Nmesh
           do mm = 1, Nmesh
              call chengdu_hotpot_Edep(L, S, J, 0, msh(ll), msh(mm),           &
                & mE_lst(ii)/PC_mN, potval)
              call SetLargeVM(tmpVM, NCH, Nmesh, ll, mm, potval)
           end do
        end do
        call lsesv_nch_full_ST_cpld(NCH, mE_lst(ii), Nmesh, tmpVM, msh, wght,  &
          & qsqr, tmpTM, lsesv_flag)
        tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii)      &
          & = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
      end do
    else
      ! For energy-independent pots
      do ii = 1, Nmesh
         do jj = 1, Nmesh
            call chengdu_hotpot_basic(L, S, J, 0, msh(ii), msh(jj), potval)
            call SetLargeVM(tmpVM, NCH, Nmesh, ii, jj, potval)
         end do
      end do
      do ii = 1, N_mE
        call lsesv_nch_full_ST_cpld(NCH, mE_lst(ii), Nmesh, tmpVM, msh, wght,  &
          & qsqr, tmpTM, lsesv_flag)
        tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii)      &
          & = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
      end do
    end if

    deallocate(tmpTM, tmpVM)

  end subroutine get_tmtrx_polelst_chengdu

  ! extrapara(1) = uptoQn
  ! extrapara(2) = x value for NLO
  ! extrapara(3) = x value for NNLO
  ! and so on
  ! 'potname' are predefined in chengdu.f90
  subroutine get_tmtrx_chengdu(L, S, J, N_mE, mE_lst, Nmesh,    &
    & msh, wght, pot_name, extrapara, tm_lst)

    integer, intent(in)     :: L, S, J, N_mE, Nmesh
    real(NER), intent(in)   :: mE_lst(:), extrapara(1:NUM_PARA_TWOBODY)
    character(len = *), intent(in)  :: pot_name
    real(NER), intent(in)   :: msh(:), wght(:)
    real(NER), intent(out)  :: tm_lst(:, :, :)

    integer :: ii, jj, ll, mm, chnid, NCH, lsesv_flag
    real(NER) :: qsqr(1:Nmesh), potval(1:2, 1:2)
    real(NER), allocatable :: tmpTM(:, :), tmpVM(:, :)
    type(lsj_symbol) :: lsj

    lsj%L = L
    lsj%S = S
    lsj%J = J
    if (.not. is_lsj_valid(lsj)) then
      write (standard_error_unit, '(a)')   &
        & 'get_tmtrx_para_chengdu : L, S, J invalid'
      return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
      NCH = 1
    else
      NCH = 2
    end if

    allocate(tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh), tmpVM(1:NCH*Nmesh, 1:NCH*Nmesh))
    do ii = 1, Nmesh
      qsqr(ii) = msh(ii)*msh(ii)
    end do

    call chengdu_dispatch(pot_name)
    if (chengdu_hotpot_is_Edep(L, S, J)) then
      ! For energy-dependent pots
      do ii = 1, N_mE
        do ll = 1, Nmesh
           do mm = 1, Nmesh
              call chengdu_hotpot_Edep(L, S, J, 0, msh(ll), msh(mm),           &
                & mE_lst(ii)/PC_mN, potval)
              call SetLargeVM(tmpVM, NCH, Nmesh, ll, mm, potval)
           end do
        end do
        call lsesv_nch_full_ST_cpld(NCH, mE_lst(ii), Nmesh, tmpVM, msh, wght,  &
          & qsqr, tmpTM, lsesv_flag)
        tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii)      &
          & = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
      end do
    else
      ! For energy-independent pots
      do ii = 1, Nmesh
         do jj = 1, Nmesh
            call chengdu_hotpot_mixed(L, S, J, int(extrapara(1)),       &
              & extrapara(2:NUM_PARA_TWOBODY), msh(ii), msh(jj), potval)
            call SetLargeVM(tmpVM, NCH, Nmesh, ii, jj, potval)
         end do
      end do
      do ii = 1, N_mE
        call lsesv_nch_full_ST_cpld(NCH, mE_lst(ii), Nmesh, tmpVM, msh, wght,  &
          & qsqr, tmpTM, lsesv_flag)
        tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii)      &
          & = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
      end do

    end if

    deallocate(tmpTM, tmpVM)

  end subroutine get_tmtrx_chengdu

  ! Return t(facCntr*msh(i), facCntr*msh(j); mE_lst(k))
  subroutine get_tmtrx_cmplx_chengdu(L, S, J, N_mE, mE_lst, Nmesh,   &
    & facCntr, msh, wght, pot_name, extrapara, tm_lst)

    integer, intent(in)             :: L, S, J, N_mE, Nmesh
    complex(NEC), intent(in)        :: facCntr, mE_lst(:)
    character(len = *), intent(in)  :: pot_name
    real(NER), intent(in)           :: msh(:), wght(:), extrapara(1:NUM_PARA_TWOBODY)
    complex(NEC), intent(out)       :: tm_lst(:, :, :)

    integer       :: ii, jj, ll, mm, chnid, NCH, lsesv_flag
    complex(NEC)  :: Cmesh(1:Nmesh), Cwght(1:Nmesh)
    complex(NEC)  :: qsqr(1:Nmesh), potval(1:2, 1:2)
    complex(NEC), allocatable :: tmpTM(:, :), tmpVM(:, :)
    type(lsj_symbol) :: lsj

    lsj%L = L
    lsj%S = S
    lsj%J = J
    if (.not. is_lsj_valid(lsj)) then
       write (standard_error_unit, '(a)')   &
        & 'get_tmtrx_cmplx_chengdu : L, S, J invalid'
       return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
       NCH = 1
    else
       NCH = 2
    end if

    allocate(tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh), tmpVM(1:NCH*Nmesh, 1:NCH*Nmesh))
    do ii = 1, Nmesh
      Cmesh(ii) = facCntr*msh(ii)
      Cwght(ii) = facCntr*wght(ii)
      qsqr(ii) = Cmesh(ii)*Cmesh(ii)
    end do

    call chengdu_dispatch(pot_name)
    if (chengdu_hotpot_is_Edep(L, S, J)) then
      ! For energy-dependent pots
      do ii = 1, N_mE
        do ll = 1, Nmesh
           do mm = 1, Nmesh
              call chengdu_cmplx_hotpot_Edep(L, S, J, 0, mE_lst(ii)/PC_mN,     &
                & Cmesh(ll), Cmesh(mm), potval)
              call SetLargeVM_cmplx(tmpVM, NCH, Nmesh, ll, mm, potval)
           end do
        end do
        call lsesv_nch_full_cmplx_cpld(NCH, mE_lst(ii), Nmesh, tmpVM,          &
          & Cmesh, Cwght, qsqr, tmpTM, lsesv_flag)
        tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii)                                   &
          & = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
      end do
    else
      ! For energy-independent pots
      do ii = 1, Nmesh
         do jj = 1, Nmesh
            call chengdu_cmplx_hotpot_mixed(L, S, J, int(extrapara(1)),        &
              & extrapara(2:NUM_PARA_TWOBODY), Cmesh(ii), Cmesh(jj), potval)
            call SetLargeVM_cmplx(tmpVM, NCH, Nmesh, ii, jj, potval)
         end do
      end do
      do ii = 1, N_mE
        call lsesv_nch_full_cmplx_cpld(NCH, mE_lst(ii), Nmesh, tmpVM, Cmesh,   &
          & Cwght, qsqr, tmpTM, lsesv_flag)
        tm_lst(1:NCH*Nmesh, 1:NCH*Nmesh, ii)      &
          & = tmpTM(1:NCH*Nmesh, 1:NCH*Nmesh)/PC_mN
      end do
    end if

    deallocate(tmpTM, tmpVM)

  end subroutine get_tmtrx_cmplx_chengdu


  ! Evaluate tmatrix(pp_lst(i), p_lst(j); mE_lst(k))
  ! In:
  !  pp_lst(1:Npp)
  !  p_lst(1:Np)
  !  facCntr*msh(1:Nmesh), facCntr*wght(1:Nmesh) will be used as complex contour
  !  to evaluate tmatrix
  subroutine get_generic_tmtrx_cmplx_chengdu(L, S, J, N_mE, mE_lst,            &
    & Npp, pp_lst, Np, p_lst, Nmesh, facCntr, msh, wght,                       &
    & pot_name, extrapara, tm_lst)

    integer, intent(in)             :: L, S, J, N_mE, Npp, Np, Nmesh
    complex(NEC), intent(in)        :: facCntr, mE_lst(:), pp_lst(:), p_lst(:)
    character(len = *), intent(in)  :: pot_name
    real(NER), intent(in)           :: msh(:), wght(:), extrapara(1:NUM_PARA_TWOBODY)
    complex(NEC), intent(out)       :: tm_lst(:, :, :)

    integer       :: ii, jj, kk, chnid, NCH, lsesv_flag
    complex(NEC)  :: Cmesh(1:Nmesh), Cwght(1:Nmesh)
    complex(NEC)  :: Cqsqr(1:Nmesh), potval(1:2, 1:2)
    complex(NEC), allocatable :: tmpTMq_p(:, :), tmpVMq_q(:, :),               &
      & tmpVMpp_q(:, :), tmpVMpp_p(:, :), tmpTMpp_p(:, :), tmpVMq_p(:, :),     &
      & loop_pp_p(:, :)
    type(lsj_symbol) :: lsj
    ! logical       :: tmpVMs_not_set

    lsj%L = L
    lsj%S = S
    lsj%J = J
    if (.not. is_lsj_valid(lsj)) then
       write (standard_error_unit, '(a)')   &
        & 'get_generic_tmtrx_cmplx_chengdu : L, S, J invalid'
       return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
      NCH = 1
    else
      NCH = 2
    end if

    do ii = 1, Nmesh
      Cmesh(ii) = facCntr*msh(ii)
      Cwght(ii) = facCntr*wght(ii)
      Cqsqr(ii) = Cmesh(ii)*Cmesh(ii)
    end do
    ! tmpVMs_not_set = .true.
    call chengdu_dispatch(pot_name)

    if (chengdu_hotpot_is_Edep(L, S, J)) then
      call set_Edep_tmatrx()
    else
      call set_EIndep_tmatrx()
    end if

  contains

    subroutine set_Edep_tmatrx()

      !$omp parallel private(ii, jj, kk, tmpTMq_p, tmpTMpp_p, tmpVMq_q, tmpVMpp_p, tmpVMpp_q, tmpVMq_p, loop_pp_p)
      allocate(tmpTMq_p(1:NCH*Nmesh, 1:NCH*Np), tmpTMpp_p(1:NCH*Npp, 1:NCH*Np),  &
        & tmpVMq_q(1:NCH*Nmesh, 1:NCH*Nmesh), tmpVMpp_q(1:NCH*Npp, 1:NCH*Nmesh), &
        & tmpVMpp_p(1:NCH*Npp, 1:NCH*Np), tmpVMq_p(1:NCH*Nmesh, 1:NCH*Np),       &
        & loop_pp_p(1:NCH*Npp, 1:NCH*Np))
      !$omp do
      do kk = 1, N_mE

        do ii = 1, Nmesh
          do jj = 1, Nmesh
            call chengdu_cmplx_hotpot_Edep(L, S, J, 0, mE_lst(kk)/PC_mN,       &
              & Cmesh(ii), Cmesh(jj), potval)
              call SetLargeVM_cmplx(tmpVMq_q, NCH, Nmesh, ii, jj, potval)
          end do
          do jj = 1, Npp
            call chengdu_cmplx_hotpot_Edep(L, S, J, 0, mE_lst(kk)/PC_mN,       &
              & pp_lst(jj), Cmesh(ii), potval)
            call SetLargeVM_generic_cmplx(tmpVMpp_q, NCH, Npp, Nmesh, jj, ii, potval)
          end do
          do jj = 1, Np
            call chengdu_cmplx_hotpot_Edep(L, S, J, 0, mE_lst(kk)/PC_mN,       &
              & Cmesh(ii), p_lst(jj), potval)
            call SetLargeVM_generic_cmplx(tmpVMq_p, NCH, Nmesh, Np, ii, jj, potval)
          end do
        end do
        do ii = 1, Npp
          do jj = 1, Np
            call chengdu_cmplx_hotpot_Edep(L, S, J, 0, mE_lst(kk)/PC_mN,       &
              & pp_lst(ii), p_lst(jj), potval)
            call SetLargeVM_generic_cmplx(tmpVMpp_p, NCH, Npp, Np, ii, jj, potval)
          end do
        end do

        call lsesv_nch_generic_cmplx(NCH, mE_lst(kk), Nmesh, Np, tmpVMq_q,     &
        & tmpVMq_p, Cmesh, Cwght, Cqsqr, tmpTMq_p, lsesv_flag)

        ! debug
        if (lsesv_flag == 1) then
          print *, "get_generic_tmtrx_cmplx_chengdu: lsesv_nch_generic_cmplx reports error."
        end if

        call get_loop_VpG0Tk_cmplx_nch_nrow_ncol(mE_lst(kk), NCH, Npp, Np,       &
          & Nmesh, Cwght, Cqsqr, tmpVMpp_q, tmpTMq_p, loop_pp_p)

        tm_lst(1:NCH*Npp, 1:NCH*Np, kk) = (tmpVMpp_p(1:NCH*Npp, 1:NCH*Np)        &
          & + loop_pp_p(1:NCH*Npp, 1:NCH*Np))/PC_mN

      end do
      !$omp end do

      deallocate(tmpTMq_p, tmpTMpp_p, tmpVMq_q, tmpVMpp_p, tmpVMpp_q, tmpVMq_p,  &
        & loop_pp_p)
      !$omp end parallel

    end subroutine set_Edep_tmatrx

    subroutine set_EIndep_tmatrx()

      allocate(tmpVMq_q(1:NCH*Nmesh, 1:NCH*Nmesh),    &
        & tmpVMpp_q(1:NCH*Npp, 1:NCH*Nmesh),          &
        & tmpVMpp_p(1:NCH*Npp, 1:NCH*Np),             &
        & tmpVMq_p(1:NCH*Nmesh, 1:NCH*Np))

      do ii = 1, Nmesh
        do jj = 1, Nmesh
          call chengdu_cmplx_hotpot(L, S, J, 0, Cmesh(ii), Cmesh(jj), potval)
          call SetLargeVM_cmplx(tmpVMq_q, NCH, Nmesh, ii, jj, potval)
        end do
        do jj = 1, Npp
          call chengdu_cmplx_hotpot(L, S, J, 0, pp_lst(jj), Cmesh(ii), potval)
          call SetLargeVM_generic_cmplx(tmpVMpp_q, NCH, Npp, Nmesh, jj, ii, potval)
        end do
        do jj = 1, Np
          call chengdu_cmplx_hotpot(L, S, J, 0, Cmesh(ii), p_lst(jj), potval)
          call SetLargeVM_generic_cmplx(tmpVMq_p, NCH, Nmesh, Np, ii, jj, potval)
        end do
      end do
      do ii = 1, Npp
        do jj = 1, Np
          call chengdu_cmplx_hotpot(L, S, J, 0, pp_lst(ii), p_lst(jj), potval)
          call SetLargeVM_generic_cmplx(tmpVMpp_p, NCH, Npp, Np, ii, jj, potval)
        end do
      end do

      !$omp parallel private(tmpTMq_p, tmpTMpp_p, loop_pp_p, lsesv_flag)
      allocate(tmpTMq_p(1:NCH*Nmesh, 1:NCH*Np),   &
          & tmpTMpp_p(1:NCH*Npp, 1:NCH*Np),       &
          & loop_pp_p(1:NCH*Npp, 1:NCH*Np))

      ! debug
      ! print *, "get_generic_tmtrx_cmplx_chengdu : num_threads = ", omp_get_thread_num()

      !$omp do
      do kk = 1, N_mE

        call lsesv_nch_generic_cmplx(NCH, mE_lst(kk), Nmesh, Np, tmpVMq_q,     &
        & tmpVMq_p, Cmesh, Cwght, Cqsqr, tmpTMq_p, lsesv_flag)

        ! debug
        if (lsesv_flag == 1) then
          print *, "get_generic_tmtrx_cmplx_chengdu: lsesv_nch_generic_cmplx reports error."
        end if

        call get_loop_VpG0Tk_cmplx_nch_nrow_ncol(mE_lst(kk), NCH, Npp, Np,       &
          & Nmesh, Cwght, Cqsqr, tmpVMpp_q, tmpTMq_p, loop_pp_p)

        tm_lst(1:NCH*Npp, 1:NCH*Np, kk) = (tmpVMpp_p(1:NCH*Npp, 1:NCH*Np)        &
          & + loop_pp_p(1:NCH*Npp, 1:NCH*Np))/PC_mN

      end do
      !$omp end do
      deallocate(tmpTMq_p, tmpTMpp_p, loop_pp_p)
      !$omp end parallel

      deallocate(tmpVMq_q, tmpVMpp_p, tmpVMpp_q, tmpVMq_p)

    end subroutine set_EIndep_tmatrx


    ! subroutine get_tmatrix_local()

    !   call lsesv_nch_generic_cmplx(NCH, mE_lst(kk), Nmesh, Np, tmpVMq_q,       &
    !     & tmpVMq_p, Cmesh, Cwght, Cqsqr, tmpTMq_p, lsesv_flag)

    !   ! debug
    !   if (lsesv_flag == 1) then
    !     print *, "get_generic_tmtrx_cmplx_chengdu: lsesv_nch_generic_cmplx reports error."
    !   end if

    !   call get_loop_VpG0Tk_cmplx_nch_nrow_ncol(mE_lst(kk), NCH, Npp, Np,       &
    !     & Nmesh, Cwght, Cqsqr, tmpVMpp_q, tmpTMq_p, loop_pp_p)

    !   tm_lst(1:NCH*Npp, 1:NCH*Np, kk) = (tmpVMpp_p(1:NCH*Npp, 1:NCH*Np)        &
    !     & + loop_pp_p(1:NCH*Npp, 1:NCH*Np))/PC_mN

    ! end subroutine get_tmatrix_local

  end subroutine get_generic_tmtrx_cmplx_chengdu

  subroutine get_h2wfs_cmplx_chengdu(Nmesh, Cmsh, Cwght, pot_name, extrapara,  &
    & h2BE, Cwf)

    integer, intent(in)       :: Nmesh
    complex(NEC), intent(in)  :: Cmsh(:), Cwght(:)
    real(NER), intent(in)     :: extrapara(1:NUM_PARA_TWOBODY)
    complex(NEC), intent(out) :: h2BE, Cwf(:)
    character(len = *), intent(in)  :: pot_name

    call chengdu_dispatch(pot_name)
    call get_h2wfs_cmplx(Nmesh, Cmsh, Cwght, PC_mN, chengdu_cmplx_hotpot,      &
      & h2BE, Cwf)

  end subroutine get_h2wfs_cmplx_chengdu

  subroutine get_generic_h2wfvals_cmplx_chengdu(Nmesh, Cmsh, Cwght, Np, p_lst, &
    & pot_name, extrapara, h2BE, wfval_lst)

    integer, intent(in)       :: Nmesh, Np
    complex(NEC), intent(in)  :: Cmsh(:), Cwght(:), p_lst(:)
    real(NER), intent(in)     :: extrapara(1:NUM_PARA_TWOBODY)
    complex(NEC), intent(out) :: h2BE, wfval_lst(:)
    character(len = *), intent(in)  :: pot_name

    complex(NEC)  :: Cwf(2*Nmesh), tmpBE

    call chengdu_dispatch(pot_name)
    call get_h2wfs_cmplx(Nmesh, Cmsh, Cwght, PC_mN, chengdu_cmplx_hotpot,      &
      & tmpBE, Cwf)

    call chengdu_dispatch(pot_name)
    call get_h2wfvals_cmplx(Nmesh, Cmsh, Cwght, PC_mN, chengdu_cmplx_hotpot,   &
      & tmpBE, Cwf, Np, p_lst, wfval_lst)

    h2BE = tmpBE

  end subroutine get_generic_h2wfvals_cmplx_chengdu


  subroutine get_h2wfs_real_potwrap(Mambda, para, regtype, Nmesh, Cmsh, Cwght, pot_name, extrapara,  &
    & h2BE, Cwf)

    integer, intent(in)       :: Nmesh, regtype
    real(NER), intent(in)     :: Cmsh(:), Cwght(:), Mambda, para(:)
    real(NER), intent(in)     :: extrapara(1:NUM_PARA_TWOBODY)
    real(NER), intent(out)    :: h2BE, Cwf(:)
    procedure(gnrc_eft_pot)   :: pot_name

    call get_h2wfs_real(Mambda, para, regtype, Nmesh, Cmsh, Cwght, PC_mN, pot_name, &
      & h2BE, Cwf)

  end subroutine get_h2wfs_real_potwrap


  subroutine get_h2wfs_vbd_real_potwrap(Mambda, para, regtype, Nmesh, Cmsh, Cwght, pot_name, extrapara,  &
    & h2BE, Cwf)

    integer, intent(in)       :: Nmesh, regtype
    real(NER), intent(in)     :: Cmsh(:), Cwght(:), Mambda, para(:)
    real(NER), intent(in)     :: extrapara(1:NUM_PARA_TWOBODY)
    real(NER), intent(out)    :: h2BE, Cwf(:)
    procedure(gnrc_eft_pot)   :: pot_name

    call get_h2wfs_vbd_real(Mambda, para, regtype, Nmesh, Cmsh, Cwght, PC_mN, pot_name, &
      & h2BE, Cwf)

  end subroutine get_h2wfs_vbd_real_potwrap

  subroutine get_generic_h2wfvals_real_potwrap(Mambda, para, regtype, Nmesh, Cmsh, Cwght, Np, p_lst, &
    & pot_name, extrapara, h2BE, wfval_lst)

    integer, intent(in)       :: Nmesh, Np, regtype
    real(NER), intent(in)     :: Cmsh(:), Cwght(:), p_lst(:), Mambda, para(:)
    real(NER), intent(in)     :: extrapara(1:NUM_PARA_TWOBODY)
    real(NER), intent(out)    :: h2BE, wfval_lst(:)
    procedure(gnrc_eft_pot)   :: pot_name

    real(NER)    :: Cwf(2*Nmesh), tmpBE

    call get_h2wfs_real(Mambda, para, regtype, Nmesh, Cmsh, Cwght, PC_mN, pot_name,      &
      & tmpBE, Cwf)
    call get_h2wfvals_real(Mambda, para, regtype, Nmesh, Cmsh, Cwght, PC_mN, pot_name,   &
      & tmpBE, Cwf, Np, p_lst, wfval_lst)

    h2BE = tmpBE

  end subroutine get_generic_h2wfvals_real_potwrap

end module eft_tmtrx
