
module eft_tmtrx_C

  use eft_tmtrx
  use iso_c_binding
  implicit none

contains

    ! Initialize global functions/variables in nnscat that will facilitate
    ! eft_tmtrx routines
    subroutine init_nnscat_chengdu_C()  &
      & bind(c, name='init_nnscat_chengdu_C')

      call initmshlg_pwo()

    end subroutine init_nnscat_chengdu_C

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
  subroutine get_FOST_tmtrx_DLSPR_C(L, S, J, N_mE, mE_lst, Mambda, Nmesh, msh,   &
    & wght, tm_lst, NCH1) bind(c, name='get_FOST_tmtrx_DLSPR')
    ! Calculate T_matrix. Fixed--size tm_lst is used as input instead of tm_lst(:, :, :) because it's more convenience for passing to C++. Ref: https://stackoverflow.com/questions/29132606/the-run-time-aborting-when-calling-c-sub-from-fortran

    integer, intent(in)     :: L, S, J, N_mE, Nmesh, NCH1
    real(NER), intent(in)   :: Mambda, mE_lst(N_mE)
    real(NER), intent(out)  :: msh(Nmesh), wght(Nmesh)
    real(NER), intent(out)  :: tm_lst(NCH1 * Nmesh, NCH1 * Nmesh, N_mE)

    call get_FOST_tmtrx_DLSPR(L, S, J, N_mE, mE_lst, Mambda, Nmesh, msh,   &
         & wght, tm_lst)

  end subroutine get_FOST_tmtrx_DLSPR_C

  subroutine get_FOST_NCH(L, S, J, NCH) bind(c, name='get_FOST_NCH_CPP')
    ! Calculate the size of tm_lst which C++ requires.

    integer, intent(in)     :: L, S, J
    integer, intent(out)     :: NCH

    integer :: chnid
    type(lsj_symbol) :: lsj

    lsj%L = L
    lsj%S = S
    lsj%J = J
    if (.not. is_lsj_valid(lsj)) then
       write (standard_error_unit, '(a)')  &
            & 'get_FOST_NCH : L, S, J invalid'
       return
    end if
    call convert_lsj_to_chnid(lsj, chnid)
    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
       NCH = 1
    else
       NCH = 2
    end if

  end subroutine get_FOST_NCH

! Under construction (08/16/2020)
! Tangent mesh to be used
! pii = TanMeshC*tan(\pi/4*(xi+1)) where ti are GAUSS-LEGENDRE abscissas
! num_pole : number of poles in mE
! poltlst : list of poles in mE. Maximum 20 poles.
! pot_name : name of the potential, with maximum length 20 chars.
! lenPotName : number of char(1 byte) in pot_name
  ! subroutine get_tmtrx_polelst_chengdu_C(L, S, J, N_mE, mE_lst, TanMeshC,      &
  !   & Nmesh, msh, wght, pot_name, lenPotName, tm_lst, num_pole, polelst, NCH1) &
  !   & bind(c, name='get_tmtrx_polelst_chengdu_C')


  !   integer, intent(in)     :: L, S, J, N_mE, Nmesh, NCH1, lenPotName
  !   integer, intent(out)    :: num_pole
  !   real(NER), intent(in)   :: TanMeshC, mE_lst(N_mE)
  !   character(kind=c_char), intent(in)  :: pot_name(lenPotName)
  !   real(NER), intent(out)  :: msh(Nmesh), wght(Nmesh), polelst(20)
  !   real(NER), intent(out)  :: tm_lst(NCH1 * Nmesh, NCH1 * Nmesh, N_mE)

  !   character(kind=c_char), target :: pot_nameT(lenPotName)
  !   character(len = lenPotName), dimension(:), pointer:: pot_nameP

  !   pot_nameT = pot_name
  !   ! Convert char array to string
  !   call c_f_pointer (c_loc(pot_nameT), pot_nameP, [1])

  !   call get_tmtrx_polelst_chengdu(L, S, J, N_mE, mE_lst,       &
  !     & TanMeshC, Nmesh, msh, wght,  &
  !     & pot_nameP(1), tm_lst, num_pole, polelst)

  ! end subroutine get_tmtrx_polelst_chengdu_C

  subroutine get_tmtrx_chengdu_C(L, S, J, N_mE, mE_lst, Nmesh,  &
    & msh, wght, pot_name, lenPotName, extrapara, tm_lst, NCH1)                &
    & bind(c, name='get_tmtrx_chengdu_C')

    integer, intent(in)     :: L, S, J, N_mE, Nmesh, NCH1, lenPotName
    real(NER), intent(in)   :: mE_lst(N_mE),                                   &
      & extrapara(NUM_PARA_TWOBODY)
    character(kind=c_char), intent(in)  :: pot_name(lenPotName)
    real(NER), intent(in)   :: msh(Nmesh), wght(Nmesh)
    real(NER), intent(out)  :: tm_lst(NCH1 * Nmesh, NCH1 * Nmesh, N_mE)

    character(kind=c_char), target :: pot_nameT(lenPotName)
    character(len = lenPotName), dimension(:), pointer:: pot_nameP

    pot_nameT = pot_name
    ! Convert char array to string
    call c_f_pointer (c_loc(pot_nameT), pot_nameP, [1])

    call get_tmtrx_chengdu(L, S, J, N_mE, mE_lst, Nmesh, msh,   &
      & wght, pot_nameP(1), extrapara, tm_lst)

  end subroutine get_tmtrx_chengdu_C

  subroutine get_tmtrx_cmplx_chengdu_C(L, S, J, N_mE, mE_lst, Nmesh,   &
    & facCntr, msh, wght, pot_name, lenPotName, extrapara, tm_lst, NCH1) &
    & bind(c, name='get_tmtrx_cmplx_chengdu_C')
    ! My C wrapper use complex double only, if you wanna use float, please write a short interface  in tMatrix::chengduT() in tMat.cpp which can switch double to float.
    integer, intent(in)       :: L, S, J, N_mE, Nmesh, NCH1, lenPotName
    ! real(NER), intent(in)     :: TanMeshC
    complex(NEC), intent(in)  :: facCntr, mE_lst(N_mE)
    character(kind=c_char), intent(in)  :: pot_name(lenPotName)
    real(NER), intent(in)     :: msh(Nmesh), wght(Nmesh), extrapara(NUM_PARA_TWOBODY)
    complex(NEC), intent(out)  :: tm_lst(NCH1 * Nmesh, NCH1 * Nmesh, N_mE)

    character(kind=c_char), target :: pot_nameT(lenPotName)
    character(len = lenPotName), dimension(:), pointer:: pot_nameP

    pot_nameT = pot_name
    ! Convert char array to string
    call c_f_pointer (c_loc(pot_nameT), pot_nameP, [1])

    call get_tmtrx_cmplx_chengdu(L, S, J, N_mE, mE_lst, Nmesh, facCntr,&
      & msh, wght, pot_nameP(1), extrapara, tm_lst)

  end subroutine get_tmtrx_cmplx_chengdu_C

  subroutine get_gnr_tmtrx_cmplx_chengdu_C(L, S, J, N_mE, mE_lst,           &
    & Npp, pp_lst, Np, p_lst, Nmesh, facCntr, msh, wght,                    &
    & pot_name, lenPotName, extrapara, tm_lst, NCH1)                        &
    & bind(c, name='get_gnr_tmtrx_cmplx_chengdu_C')
    ! My C wrapper use complex double only, if you wanna use float, please write a short interface in tMatrix::chengduT() in tMat.cpp which can switch double to float.
    integer, intent(in)       :: L, S, J, N_mE, Npp, Np, Nmesh, NCH1, lenPotName
    ! real(NER), intent(in)     :: TanMeshC
    complex(NEC), intent(in)  :: facCntr, mE_lst(N_mE), pp_lst(Npp), p_lst(Np)
    character(kind=c_char), intent(in)  :: pot_name(lenPotName)
    real(NER), intent(in)     :: msh(Nmesh), wght(Nmesh), extrapara(NUM_PARA_TWOBODY)
    complex(NEC), intent(out)  :: tm_lst(NCH1 * Npp, NCH1 * Np, N_mE)

    character(kind=c_char), target :: pot_nameT(lenPotName)
    character(len = lenPotName), dimension(:), pointer:: pot_nameP

    pot_nameT = pot_name
    ! Convert char array to string
    call c_f_pointer (c_loc(pot_nameT), pot_nameP, [1])

    call get_generic_tmtrx_cmplx_chengdu(L, S, J, N_mE, mE_lst,            &
      & Npp, pp_lst, Np, p_lst, Nmesh, facCntr, msh, wght,                 &
      & pot_nameP(1), extrapara, tm_lst)

  end subroutine get_gnr_tmtrx_cmplx_chengdu_C

  subroutine get_generic_h2wfvals_cmplx_chengdu_C(Nmesh, Cmsh, Cwght,          &
    & Np, p_lst, pot_name, lenPotName, extrapara, H2BE, wfval_lst)             &
    & bind(c, name='get_generic_h2wfvals_cmplx_chengdu_C')

    integer, intent(in)       :: Nmesh, Np, lenPotName
    complex(NEC), intent(in)  :: Cmsh(Nmesh), Cwght(Nmesh), p_lst(Np)
    real(NER), intent(in)     :: extrapara(1:NUM_PARA_TWOBODY)
    complex(NEC), intent(out) :: H2BE, wfval_lst(2*Np)
    character(kind=c_char), intent(in)  :: pot_name(lenPotName)

    character(kind=c_char), target :: pot_nameT(lenPotName)
    character(len = lenPotName), dimension(:), pointer:: pot_nameP
    integer :: ii

    pot_nameT = pot_name
    ! Convert char array to string
    call c_f_pointer (c_loc(pot_nameT), pot_nameP, [1])

    call get_generic_h2wfvals_cmplx_chengdu(Nmesh, Cmsh, Cwght, Np, p_lst,     &
      & pot_nameP(1), extrapara, H2BE, wfval_lst)

    ! debug
    ! print *, "Np = ", Np
    ! do ii = Np - 9, Np
    !   print *, "ii = ", ii, "  ", wfval_lst(ii+Np)
    ! end do

  end subroutine get_generic_h2wfvals_cmplx_chengdu_C

end module eft_tmtrx_C
