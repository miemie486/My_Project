

program fit_Tmtrx_smpl

  use drv_fitting
  implicit none


  integer               :: npar_Qn(1:10)
  real(NER)             :: initial_fraction, tolerance
  character(len = 256)  :: fitphs_fname, inpar_fname, out_fname, errmsg
  logical               :: succeeded
  character(len = 78), dimension(1:4), parameter  :: parname_Qn = (/                    &
    & '# Mambda  C0                                                                  ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0                                                          ', &
    & '# Mambda  C0  C2  D0  C3  D1                                                  '  &
    & /)
  integer :: regtype
  class(obj_smplchn), allocatable   :: chn
  real(NER)                   :: t0, t1

  call cpu_time(t0)

  fitphs_fname = 'fit1s0_phs.in'
  inpar_fname = 'fit1s0_inpar.in'
  out_fname = 'fit1s0_result.out'
  npar_Qn(1:3) = NPAR_LY1S0_SMPL
  initial_fraction = 1.0E-3_NER
  tolerance = 1.0E-7_NER
  errmsg = 'fit_Tmtrx_1s0: Something went wrong'


    PC_gA = 1.29_NER
    PC_mN = 939.0_NER
    PC_mpi = 138.0_NER
    PC_fpi = 92.4_NER

!   call init_smpl()
!   call PrintOutFitSngl_fitting(GetNPWithC0Sngl_del_fitpwa,                         &
!     & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
!     & initial_fraction, tolerance, succeeded)
  call PrintOutFitSnglTmtrx_fitting(Get1s0Tmtrx_smpl_fitpwa,                         &
    & fitphs_fname, inpar_fname, npar_Qn, out_fname, parname_Qn,               &
    & initial_fraction, tolerance, succeeded)

  if (.not. succeeded) then
    write(standard_error_unit, '(a)') errmsg
  else
    call cpu_time(t1)
    print '(a, f9.2, a)', 'running time = ', t1 - t0, ' seconds.'
  end if
contains

!     subroutine init_smpl()


!     allocate(obj_smplchn::chn)

!     ! NLO pot does *not* vanishes
!     chn%vanishing(1) = .false.
!     chn%npar(1:3) = NPAR_LY1S0_SMPL

!     chn%pVfunc_V0 => VLO_withc0_epwrap

!     chn%pVfunc_V1 => VNLO_MMWLY_epwrap
!     chn%V1_cheb_base_n = 4
!     chn%V1_cheb_bnd = 1.0E-3_NER

!     ! chn%pVfunc_V2 => VN2LO_sprbl_epwrap
!     ! chn%V2_cheb_base_n = V2_BASEN_sprbl_mvfs
!     ! chn%V2_cheb_bnd = V2_BND_sprbl_mvfs

!     ! chn%pVfunc_V3 => VN3LO_sprbl_epwrap
!     ! chn%V3_cheb_base_n = V3_BASEN_sprbl_mvfs
!     ! chn%V3_cheb_bnd = V3_BND_sprbl_mvfs

!   end subroutine init_smpl

  subroutine Get1s0Tmtrx_smpl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,      &
    & num_paras, klist, numk, amplitude)

    real(NER), intent(in)       :: Mambda, paralst(:), klist(:)
    integer, intent(in)         :: Nmesh, chnid, uptoQn, num_paras, numk
    real(NER), intent(out)      :: amplitude(:, :)

    integer :: regtype, ii
    class(obj_smplchn), pointer :: ptrchn

    regtype = REGTYPE_GAUSSIAN
    allocate(obj_smplchn::ptrchn)
    call create_smplchn(ptrchn, chnid, uptoQn)
    call SetUpVfunc_mmwly_smpl(ptrchn)


!     call ptrchn%load_inputs(num_paras, paralst)
!     call ptrchn%init(Nmesh, regtype, Mambda)

!     do ii = 1, numk
!       call ptrchn%update_k(klist(ii))
!       call ptrchn%get_allTQn()
!       amplitude(ii, 1: ptrchn%uptoQn+1) = abs(ptrchn%TQn(1, 1, ptrchn%uptoQn+1))
! !   print *, "ii=",ii,"k=", klist(ii), 'Tmtrx=', amplitude(ii,1)
!     end do

    call get_phase_as_value_for_klist_smpl(ptrchn, Nmesh, regtype, Mambda,     &
      & klist, numk, paralst, num_paras, amplitude)

    call ptrchn%erase()
    deallocate(ptrchn)

  end subroutine Get1s0Tmtrx_smpl_fitpwa

end program fit_Tmtrx_smpl