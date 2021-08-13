! Rui Peng 12/20/2020

module vbar_3p0_wfs

  use nneft_type
  use util_rootfinder
  use nneft_lsesv
  use nneft_eigens
  use potwrap
  use util_gauleg
  implicit none

  private
  real(NER), parameter :: FINDH2BE_TOL = 1.0E-9_NER
  integer, parameter   :: NCH = 1, L = 1, S = 1, J = 0
  logical              :: Rmsh_inited  = .false., Cwf_inited = .false., spuBED_inited = .false.
  real(NER)            :: spuBDE
  integer              :: Nmesh_spu, int_Mambda = 0
  real(NER), allocatable          :: Rmsh_spu(:), Rwght_spu(:), Cwf_spu(:)

  public  :: VbarLO_3p0_withc0_chengdu_epwrap, get_h2wfs_real_wrap, get_eta_3p0_real

contains

! the real type with potwrap potential

  subroutine get_h2wfvals_real_wrap(Mambda, para, regtype, Nmesh, Rmsh, Rwght, AvgMN, pot, &
    h2BE, Rwf, Np, plst, wfval_lst)

    procedure(gnrc_eft_pot)   :: pot
    integer, intent(in)       :: Nmesh, Np, regtype
    real(NER), intent(in)     :: AvgMN
    real(NER), intent(in)     :: Rmsh(:), Rwght(:), h2BE, Rwf(:), plst(:)
    real(NER), intent(in)     :: Mambda, para(:)
    real(NER), intent(out)    :: wfval_lst(:)

    integer :: ii, jj
    real(NER)  :: potval(2, 2), p, tmpwf, tmpprdct

    do ii = 1, Np
      tmpprdct = 0.0_NER
      p = plst(ii)
      do jj = 1, Nmesh
        call pot(L, S, J, regtype, Mambda, para, p, Rmsh(jj), potval)
        tmpwf = Rwf(jj)
        tmpprdct = tmpprdct + Rwght(jj)*(Rmsh(jj)*Rmsh(jj))*potval(1,1)*tmpwf
      end do
      wfval_lst(ii) = tmpprdct/((h2BE - p*p/AvgMN)*AvgMN)
    end do

  end subroutine get_h2wfvals_real_wrap

! For a given potwrap potential with real type

  subroutine get_h2wfs_real_wrap(Mambda, para, regtype, Nmesh, Rmsh, Rwght, AvgMN, pot, h2BE, Cwf)

    procedure(gnrc_eft_pot)         :: pot
    integer, intent(in)             :: Nmesh, regtype
    real(NER), intent(in)           :: AvgMN
    real(NER), intent(in)           :: Rmsh(:), Rwght(:)
    real(NER), intent(in)           :: Mambda, para(:)
    real(NER), intent(out)          :: h2BE, Cwf(:)

    real(NER)     :: AMat(1:2*Nmesh,1:2*Nmesh)
    real(NER)     :: eta, normfac, Gamma(1:2*Nmesh)
    integer       :: dd, flag, uptoQn, ii
    real(NER)     :: x1, x2, Reh2BE

    x1 = -1500.0_NER
    x2 = -100.0_NER
    if(int(Mambda) .lt. 1000) x2 = -500.0_NER
    uptoQn = 0

    call findroot_h2wfs(eta_no, x1, x2, FINDH2BE_TOL, Reh2BE)
    normfac = 0.0_NER
    do ii = 1, Nmesh
      normfac = normfac + Gamma(ii)*Gamma(ii)*Rwght(ii)*Rmsh(ii)*Rmsh(ii)
    end do
    do ii = 1, NCH*Nmesh
      Cwf(ii) = Gamma(ii)/sqrt(normfac)
    end do
    h2BE = Reh2BE

  contains

    function eta_no(E)

      real(NER), intent(in) :: E

      real(NER)             :: eta_no
      integer               :: ss, tt
      real(NER)    :: potval(1:2, 1:2), tmp(1:2, 1:2)
      real(NER)    :: p,k

      do ss=1, Nmesh
        p=Rmsh(ss)
        do tt=1, Nmesh
          k=Rmsh(tt)
          call pot(L, S, J, regtype, Mambda, para, p, k, potval)
          tmp = Rwght(tt)*(k*k)/(E-p*p/AvgMN)*potval/AvgMN
          call SetLargeVM(AMat, NCH, Nmesh, ss, tt, tmp)
        end do
      end do

      call get_eigens_matrix_real_modeigens(AMat, NCH*Nmesh, eta, Gamma, flag)
      eta_no = eta - 1.0_NER

    end function eta_no

  end subroutine get_h2wfs_real_wrap

  subroutine findroot_h2wfs(func, x1, x2, tol, rslt)

    REAL(NER),intent(in)  :: x1, x2, tol
    REAL(NER),intent(out) :: rslt
    INTERFACE
      FUNCTION func(x)
        USE nneft_type
        IMPLICIT NONE
        REAL(NER), INTENT(IN) :: x
        REAL(NER) :: func
      END FUNCTION func
    END INTERFACE
    rslt = zbrent(func, x1, x2, tol)

  end subroutine findroot_h2wfs

  subroutine get_eta_3p0_real(Mambda, para, regtype, Nmesh, Rmsh, Rwght, AvgMN, pot, h2BE, Cwf, tmp_e, eta1)

    procedure(gnrc_eft_pot)         :: pot
    integer, intent(in)             :: Nmesh, regtype
    real(NER),    intent(in)        :: AvgMN
    real(NER), intent(in)           :: Rmsh(:), Rwght(:)
    real(NER), intent(in)           :: Mambda, para(:)
    real(NER), intent(out)          :: h2BE, Cwf(:), tmp_e(100), eta1(100)

    real(NER)     :: AMat(1:2*Nmesh,1:2*Nmesh)
    real(NER)     :: eta, normfac, Gamma(1:2*Nmesh)
    integer       :: dd, flag, uptoQn, ii
    real(NER)     ::  Reh2BE

    uptoQn = 0

    do ii = 1, 100
      tmp_e(ii) = -4000.0_NER*ii/100
!       tmp_e(ii) = 445.94_NER-0.05_NER*ii/100
      eta1(ii) = eta_no(tmp_e(ii))
    end do
  contains

    function eta_no(E)

      real(NER), intent(in) :: E

      real(NER)             :: eta_no
      integer               :: ss, tt
      real(NER)    :: potval(1:2, 1:2), tmp(1:2, 1:2)
      real(NER)    :: p,k

      do ss=1, Nmesh
        p=Rmsh(ss)
        do tt=1, Nmesh
          k=Rmsh(tt)
          call pot(L, S, J, regtype, Mambda, para, p, k, potval)
          tmp = Rwght(tt)*(k*k)/(E-p*p/AvgMN)*potval/AvgMN
          call SetLargeVM(AMat, NCH, Nmesh, ss, tt, tmp)
        end do
      end do

      call get_eigens_matrix_real_modeigens(AMat, NCH*Nmesh, eta, Gamma, flag)
      eta_no = eta - 1.0_NER

    end function eta_no

  end subroutine get_eta_3p0_real

!   subroutine get_h2wfs_vbd_real_wrap(Mambda, para, regtype, Nmesh, Rmsh, Rwght, AvgMN, pot, h2BE, Cwf)

!     procedure(gnrc_eft_pot)         :: pot
!     integer, intent(in)             :: Nmesh, regtype
!     real(NER), intent(in)           :: AvgMN
!     real(NER), intent(in)           :: Rmsh(:), Rwght(:)
!     real(NER), intent(in)           :: Mambda, para(:)
!     real(NER), intent(out)          :: h2BE, Cwf(:)

!     real(NER)     :: AMat(1:NCH*Nmesh,1:NCH*Nmesh)
!     real(NER)     :: eta, normfac, Gamma(1:NCH*Nmesh)
!     integer       :: dd, flag, uptoQn, ii
!     real(NER)     :: x1, x2, Reh2BE

!     x1 = -3000.0_NER
!     x2 = -1400.0_NER
!     uptoQn = 0

!     call findroot_h2wfs(eta_no, x1, x2, FINDH2BE_TOL, Reh2BE)
!     normfac = 0.0_NER
!     do ii = 1, Nmesh
!       normfac = normfac + (Gamma(ii)*Gamma(ii)          &
!         & + Gamma(ii + Nmesh)*Gamma(ii + Nmesh))*Rwght(ii)*Rmsh(ii)*Rmsh(ii)
!     end do
!     do ii = 1, NCH*Nmesh
!       Cwf(ii) = Gamma(ii)/sqrt(normfac)
!     end do
!     h2BE = Reh2BE

!   contains

!     function eta_no(E)

!       real(NER), intent(in) :: E

!       real(NER)             :: eta_no
!       integer               :: ss, tt
!       real(NER)    :: potval(1:2, 1:2), tmp(1:2, 1:2)
!       real(NER)    :: p,k

!       do ss=1, Nmesh
!         p=Rmsh(ss)
!         do tt=1, Nmesh
!           k=Rmsh(tt)
!           call pot(L, S, J, regtype, Mambda, para, p, k, potval)
!           tmp = Rwght(tt)*(k*k)/(E-p*p/AvgMN)*potval/AvgMN
!           call SetLargeVM(AMat, NCH, Nmesh, ss, tt, tmp)
!         end do
!       end do

!       call get_eigens_matrix_real_modeigens(AMat, NCH*Nmesh, eta, Gamma, flag)
!       eta_no = eta - 1.0_NER

!     end function eta_no

!   end subroutine get_h2wfs_vbd_real_wrap

  subroutine get_spuBDE_real(Mambda, para, regtype, Nmesh, Rmsh, Rwght, AvgMN, pot)

    procedure(gnrc_eft_pot)         :: pot
    integer, intent(in)             :: Nmesh, regtype
    real(NER), intent(in)           :: AvgMN
    real(NER), intent(in)           :: Rmsh(:), Rwght(:)
    real(NER), intent(in)           :: Mambda, para(:)

    real(NER)           :: h2BE, Cwf(NCH*Nmesh)

    allocate(Cwf_spu(NCH*Nmesh))

    call get_h2wfs_real_wrap(Mambda, para, regtype, Nmesh, Rmsh, Rwght, AvgMN, pot, h2BE, Cwf)
    spuBDE = h2BE
    Cwf_spu = Cwf

    Cwf_inited = .true.
    spuBED_inited = .true.
  end subroutine get_spuBDE_real

  subroutine init_Rmesh_wrap(N, Rmsh, Rwght)

    integer, intent(in)       ::N
    real(NER), intent(in)     :: Rmsh(1:N), Rwght(1:N)

    Nmesh_spu = N
    allocate(Rmsh_spu(N), Rwght_spu(N))
    Rmsh_spu = Rmsh
    Rwght_spu = Rwght
    Rmsh_inited = .true.
  end subroutine init_Rmesh_wrap

  subroutine relase_spu_flag_wrap()

    if (allocated(Rmsh_spu)) deallocate(Rmsh_spu)
    if (allocated(Rwght_spu)) deallocate(Rwght_spu)
    if (allocated(Cwf_spu)) deallocate(Cwf_spu)
    Cwf_inited = .false.
    spuBED_inited = .false.
    Rmsh_inited = .false.
  end subroutine relase_spu_flag_wrap

  subroutine VbarLO_3p0_withc0_chengdu_epwrap(L, S, J, regtype, Mambda, para, p1, p2, potval)

    integer, intent(in)     :: L, S, J, regtype
    real(NER), intent(in)   :: Mambda, p1, p2, para(:)
    real(NER), intent(out)  :: potval(1:2, 1:2)

    real(NER)    :: v_withc0(1:2, 1:2) , v_spu(1:2, 1:2)
    real(NER)    :: wfval1(2), wfval2(2), tmpplst1(1), tmpplst2(1)
    integer      :: N
    real(NER), allocatable    :: Rmsh(:), Rwght(:)

    if (int_Mambda /= int(Mambda)) then
      call relase_spu_flag_wrap()
      int_Mambda = int(Mambda)
      N = 100
      allocate(Rmsh(N), Rwght(N))
      call feb_tanabsc(N, 0.5_NER*Mambda, Rmsh, Rwght)
      call init_Rmesh_wrap(N, Rmsh, Rwght)
      deallocate(Rmsh, Rwght)
    end if

    if(.not. spuBED_inited) then
      call get_spuBDE_real(Mambda, para, regtype, Nmesh_spu, Rmsh_spu, Rwght_spu, PC_mN, VLO_withc0_epwrap)
    end if
    tmpplst1(1) = p1
    tmpplst2(1) = p2
    call get_h2wfvals_real_wrap(Mambda, para, regtype, Nmesh_spu, Rmsh_spu, Rwght_spu, PC_mN, VLO_withc0_epwrap, &
      & spuBDE, Cwf_spu, 1, tmpplst1, wfval1)
    call get_h2wfvals_real_wrap(Mambda, para, regtype, Nmesh_spu, Rmsh_spu, Rwght_spu, PC_mN, VLO_withc0_epwrap, &
      & spuBDE, Cwf_spu, 1, tmpplst2, wfval2)
    v_spu = 0.0_NER
    v_spu(1,1) = wfval1(1)*wfval2(1)
!     v_spu = -200*spuBDE*v_spu
     v_spu = -2*PC_mN*spuBDE*v_spu

    call VLO_withc0_epwrap(L, S, J, regtype, Mambda, para, p1, p2, v_withc0)
    potval = v_spu + v_withc0

  end subroutine VbarLO_3p0_withc0_chengdu_epwrap

end module vbar_3p0_wfs