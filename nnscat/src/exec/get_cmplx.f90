
program get_cmplx

  ! use nneft_type
  use nneft_lsesv
  use util_gauleg
  use nneft_loops
  use cmplx_epmod
  use potwrap
  use eft_tmtrx
  use eft_phaseconv
  implicit none

  integer             ::  ii, jj, L, S, J, flag, GAUN = 2
  integer             :: N = 200, CH = 1, regtype = REGTYPE_GAUSSIAN, nk
  complex(NEC)        :: potval(2, 2), intgrl_c, f, loop_sngl, mE_lst(1:10),   &
    & klst(1:10), tmlst(2, 2, 1)
  real(NER), allocatable     :: msh(:), wght(:), qsqr(:)
  complex(NEC), allocatable  :: Cqsqr(:), Cmsh(:), Cwght(:), tm_msh(:, :, :)

  complex(NEC), allocatable  :: Vknch(:, :), Vk1nch(:, :), Tknch(:, :),   &
    & VMnch(:, :), loopnch(:, :), Vknrow(:, :)

  complex(NEC), allocatable  :: Vk(:), Vk1(:), Tk(:), VM(:, :)
  complex(NEC)  :: loop

  real(NER)     :: k, mE, Mambda, LMRatio = 3.0_NER,    &
    & CMRatio = 0.7, phi = 3.0, realV(2, 2), intgrl,    &
    & extrapara(1:NUM_PARA_TWOBODY), phs(1:10)

  character(len = 1024)   :: cmdbuff

  real(NER), parameter  ::  &
    & MMWLY1s0_600_para(1:6) = (/                                            &
    & -0.39028006210685457E-002_NER,       -0.19095743809213869E-002_NER,    &
    &  0.13059371072198318E-007_NER,        0.0_NER,                         &
    &  0.0_NER,                             0.0_NER/)

  real(NER), parameter  ::  &
    & Cs3s1_600_set2_para(1:7) = (/ -1.4410457830873011E-003_NER, 0.0_NER,   &
    & 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  ! real(NER), parameter  ::  &
  !   & Cs3p0_600_set2_para(1:5) =  (/                                        &
  !   & 0.47499661698324306E-007_NER, 0.0_NER, 0.0_NER, 0.0_NER, 0.0_NER/)

  CH = 2
  k = 10.0_NER
  mE = k*k
  Mambda = 600.0_NER
  ! CMRatio = 0.6_NER
  LMRatio = 1.5_NER
  ! phi = 1.0_NER

  call get_command_argument(1, cmdbuff)
  read(cmdbuff, *) N
  call get_command_argument(2, cmdbuff)
  read(cmdbuff, *) CMRatio
  call get_command_argument(3, cmdbuff)
  read(cmdbuff, *) phi
  call get_command_argument(4, cmdbuff)
  read(cmdbuff, *) GAUN
  call get_command_argument(5, cmdbuff)
  read(cmdbuff, *) mE
  if (mE > 0.0_NER) then
    k = sqrt(mE)
  end if

  f = exp(-IMGUNT_NE*phi/180.0*PI_NE)
  ! L = 0
  ! S = 0
  ! J = 0
  L = 0
  S = 1
  J = 1

  allocate(Cqsqr(N), Cmsh(N), Cwght(N), tm_msh(CH*N, CH*N, 10), msh(N), wght(N), qsqr(N))
  allocate(Vk(N), Vk1(N), Tk(N), VM(N*CH, N*CH))
  allocate(Vknch(N*CH, CH), Vknrow(CH, N*CH), Vk1nch(N*CH, CH), Tknch(N*CH, CH),      &
    & VMnch(N*CH, N*CH), loopnch(CH, CH))

  print *, "phi = ", phi
  ! print *, "LMRatio = ", LMRatio
  print *, "N = ", N
  ! print *, "GAUSS_N = ", GAUN
  print *, "mE, k = ", mE, k

  call init_mesh()

  ! klst(1) = (65.001475478876003_NER, -5.6868922168599827_NER)
  ! klst(1) = (10.825270090785095_NER, -0.94708841254994580_NER)

  call generic_tmtrx()

  deallocate(Cqsqr, Cmsh, Cwght, msh, wght, qsqr, tm_msh)
  deallocate(Vk, Vk1, Tk, VM)
  deallocate(Vknch, Vk1nch, Vknrow, Tknch, VMnch, loopnch)

contains

  subroutine generic_tmtrx()

    klst(1) = cmplx(k, 0.0_NER, NEC)
    mE_lst(1) = cmplx(mE, 0.0_NER, NEC)
    call get_generic_tmtrx_cmplx_chengdu(L, S, J, 1, mE_lst(1:1), 1, klst(1:1),   &
      & 1, klst(1:1), N, f, msh, wght, "chengdu_MMWLY_cmplx_600", extrapara, tmlst)
    ! print *, "p = ", klst(1)
    ! print *, "tm = ", tmlst(1:2, 1:2, 1)*PC_mN
    call t2d1d2e_np_series(k, tmlst(1:2, 1:2, 1)*PC_mN, phs, flag)
    if (flag == UNITARITY_BROKEN) then
      print *, "UNITARITY_BROKEN"
    else
      print *, phs(1), phs(2), phs(3)
    end if

  end subroutine generic_tmtrx

  subroutine show_tmtrx_full_offshell(mNE)

    ! integer, intent(in)       :: NmNE
    complex(NEC), intent(in)  :: mNE

    ! integer :: ii
    complex(NEC)  :: local_mE_lst(1:10)

    local_mE_lst(1) = mNE;
    call get_tmtrx_cmplx_chengdu(L, S, J, 1, local_mE_lst(1:1), N,   &
      & f, msh, wght, "chengdu_MMWLY_cmplx_600", extrapara, tm_msh)
    print *, "f*msh(nk) = ", f*msh(nk)
    print *, "tm_msh = ", tm_msh(nk, nk, 1)*PC_mN, tm_msh(nk, nk+N, 1)*PC_mN,     &
    & tm_msh(nk+N, nk, 1)*PC_mN, tm_msh(nk+N, nk+N, 1)*PC_mN

  end subroutine show_tmtrx_full_offshell

  real(NER) function arctanh(x)
    real(NER), intent(in)  :: x
    arctanh = 0.5*(log((1.0_NER + x)/(x - 1.0_NER)))
  end function arctanh

  real(NER) function intgrnd(x)

    real(NER), intent(in)  :: x

    intgrnd = (exp(-x**GAUN/(Mambda**GAUN))-exp(-k**GAUN/(Mambda**GAUN)))/(mE - x*x)

  end function intgrnd

  real(NER) function intgrnd_bare(x)

    real(NER), intent(in)  :: x

    intgrnd_bare = exp(-x**GAUN/(Mambda**GAUN))/(mE - x*x)

  end function intgrnd_bare

  complex(NEC) function intgrnd_cmplx(x)

    complex(NEC), intent(in)  :: x

    intgrnd_cmplx = exp(-x**GAUN/(Mambda**GAUN))/(mE - x*x)

  end function intgrnd_cmplx

  subroutine get_intgr_cmplx(func, s)

    interface
      function func(x)
        import NEC
        complex(NEC)  :: func
        complex(NEC), intent(in)  :: x
      end function func
    end interface
    complex(NEC)  :: s

    integer :: ii

    s = 0.0_NER

    do ii = 1, N
      s = s + Cwght(ii)*func(Cmsh(ii))
    end do

  end subroutine get_intgr_cmplx

  subroutine get_intgr(func, s)

    interface
      function func(x)
        import NER
        real(NER)  :: func
        real(NER), intent(in)  :: x
      end function func
    end interface
    real(NER)  :: s

    integer :: ii

    s = 0.0_NER

    do ii = 1, N
      s = s + wght(ii)*func(msh(ii))
    end do

  end subroutine get_intgr

  subroutine init_mesh()

    call feb_tanabsc(N, Mambda*CMRatio, msh, wght)
    print *, "TanMesh in use, ", "CMRatio = ", CMRatio

    ! call feb_glabsc_2pt(N, msh, wght, 0.0_NER, Mambda*LMRatio)
    ! print *, "LinMesh in use"

    do ii = 1, N
      Cmsh(ii) = msh(ii)*f
      Cwght(ii) = wght(ii)*f

      ! debug
      ! print *, "ii = ", ii, ", Cmsh(ii) = ", Cmsh(ii)

      Cqsqr(ii) = Cmsh(ii)*Cmsh(ii)
    end do
  end subroutine init_mesh

  subroutine init_pots()

    complex(NEC)  :: pii, pjj

    PC_gA = 1.29_NER
    PC_mN = 939.0_NER
    PC_mpi = 138.0_NER
    PC_fpi = 92.4_NER

    do ii = 1, N
      ! call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda,                &
        ! & MMWLY1s0_600_para, Cmsh(ii), cmplx(k, 0.0_NER, NEC), potval)
      call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda,                &
        & Cs3s1_600_set2_para, Cmsh(ii), cmplx(k, 0.0_NER, NEC), potval)

      Vk(ii) = potval(1, 1)
      call SetTkFromTQ(Vknch, CH, N, ii, potval)

      ! debug
      ! print *, "ii, Vk(ii) = ", ii, Vk(ii, 1)

      ! call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda,                &
      !   &  MMWLY1s0_600_para, cmplx(k, 0.0_NER, NEC), Cmsh(ii), potval)
      call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda,                &
        &  Cs3s1_600_set2_para, cmplx(k, 0.0_NER, NEC), Cmsh(ii), potval)

      Vk1(ii) = potval(1, 1)
      call SetTkFromTQ(Vk1nch, CH, N, ii, potval)
      call SetLargeVM_generic_cmplx(Vknrow, CH, 1, N, 1, ii, potval)

      ! debug
      ! print *, "ii, Vk1(ii) = ", ii, Vk1(ii, 1)

      do jj = 1, N
        ! call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda,              &
        !   &  MMWLY1s0_600_para, Cmsh(ii), Cmsh(jj), potval)
        call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda,              &
          &  Cs3s1_600_set2_para, Cmsh(ii), Cmsh(jj), potval)

        VM(ii, jj) = potval(1, 1)*Cmsh(jj)**2/(mE - Cmsh(jj)**2)
        call SetLargeVM_cmplx(VMnch, CH, N, ii, jj, potval)

        ! debug
        ! print *, "ii = ", ii, ", jj = ", jj, ", VM = ", VM(ii, jj)
        ! call VLO_withc0_epwrap(L, S, J, regtype, Mambda,                    &
        !   &  MMWLY1s0_600_para, msh(ii), msh(jj), realV)
        ! print *, "realV = ", realV(1, 1)

      end do
    end do

  end subroutine init_pots


  subroutine get_loop_VkG0Tk_cmplx_sngl(mE, N, meshpts, wghts, lsqr, Vk, Tk,   &
    & loop_sum)

    integer, intent(in)       :: N
    complex(NEC), intent(in)  :: mE
    complex(NEC), intent(in)  :: meshpts(:), wghts(:), lsqr(:), Vk(:)
    complex(NEC), intent(in)  :: Tk(:)
    ! complex(NEC), intent(out) :: loop_sum(2, 2)
    complex(NEC), intent(out) :: loop_sum

    integer     :: ii

    loop_sum = 0.0_NER
    do ii = 1, N
      loop_sum = loop_sum + wghts(ii)/(mE - lsqr(ii))*lsqr(ii)*Vk(ii)*Tk(ii)
    end do

  end subroutine get_loop_VkG0Tk_cmplx_sngl

  ! Reserved for get_cmplx.f90
  ! \int_0^\infty dl l^2 Vk(l)*Tk(l)/(k^2 - l^2 + i\epsilon)
  subroutine get_loop_VkG0Tk_cmplx_nch(mE, NCH, N, meshpts, wghts, lsqr, Vk, Tk, loop_sum)

    integer, intent(in)       :: N, NCH
    complex(NEC), intent(in)  :: mE
    complex(NEC), intent(in)  :: meshpts(:), wghts(:), lsqr(:), Vk(:, :)
    complex(NEC), intent(in)  :: Tk(:, :)
    ! complex(NEC), intent(out) :: loop_sum(2, 2)
    complex(NEC), intent(out) :: loop_sum(:, :)

    integer     :: ii
    complex(NEC)  :: v(NCH, NCH), tmp(NCH, NCH)

    loop_sum(1:NCH, 1:NCH) = 0.0_NER
    do ii = 1, N
      call GetTQFromTk(Vk, NCH, N, ii, v)
      call GetTQFromTk(Tk, NCH, N, ii, tmp)

      ! debug
      ! print *, "ii", ii
      ! print *, "v = ", v
      ! print *, "Tk = ", tmp

      loop_sum(1:NCH, 1:NCH) = loop_sum(1:NCH, 1:NCH) + wghts(ii)/(mE - lsqr(ii))*lsqr(ii)*matmul(v, tmp)

    end do

  end subroutine get_loop_VkG0Tk_cmplx_nch

  ! call get_intgr_cmplx(intgrnd_cmplx, intgrl_c)
  ! call get_intgr(intgrnd, intgrl)
  ! print *, "CmplxInt: ", intgrl_c
  ! print *, "LinMesh : ", intgrl + exp(-k**GAUN/(Mambda**GAUN))/k * ( arctanh(Mambda/k*LMRatio)   &
  !   & - 0.5*PI_NE*IMGUNT_NE )
  ! print *, "TanMesh : ", intgrl + exp(-k**GAUN/(Mambda**GAUN))/k * (- 0.5*PI_NE*IMGUNT_NE)
  ! call get_intgr(intgrnd_bare, intgrl)
  ! print *, "Bare    : ", intgrl

  ! call init_pots()

  ! call lsesv_nonsing_half_cmplx_sngl(N, Vk, VM, Cwght, Tk, flag)
  ! call get_loop_VkG0Tk_cmplx_sngl(cmplx(mE, 0.0_NER, NEC), N, Cmsh, Cwght,    &
  !   Cqsqr, Vk1, Tk, loop_sngl)
  ! call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda, MMWLY1s0_600_para,   &
  !   & cmplx(k, 0.0_NER, NEC), cmplx(k, 0.0_NER, NEC), potval)
  ! call VLO_withc0_epwrap(L, S, J, regtype, Mambda, MMWLY1s0_600_para,   &
  !   & k, k, realV)
  ! print *, "flag = ", flag
  ! print *, "potval = ", potval(1, 1)
  ! print *, "realV = ", realV(1, 1)
  ! print *, "loop = ", loop_sngl
  ! print *, "potval + loop = ", potval(1, 1) + loop_sngl

  ! call lsesv_nch_half_cmplx_cpld(CH, cmplx(mE, 0.0_NER, NEC), N, VMnch, Vknch,  Cmsh, Cwght, Cqsqr, Tknch, flag)
  ! call get_loop_VkG0Tk_cmplx_nch(cmplx(mE, 0.0_NER, NEC), CH, N, Cmsh, Cwght, &
  !   & Cqsqr, Vk1nch, Tknch, loopnch)
  ! call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda, MMWLY1s0_600_para,    &
  !   & cmplx(k, 0.0_NER, NEC), cmplx(k, 0.0_NER, NEC), potval)

  ! call lsesv_nch_generic_cmplx(CH, cmplx(mE, 0.0_NER, NEC), N, 1, VMnch, Vknch, Cmsh, Cwght, Cqsqr, Tknch, flag)
  ! call get_loop_VpG0Tk_cmplx_nch_nrow_ncol(cmplx(mE, 0.0_NER, NEC), CH, 1, 1, N, Cwght, Cqsqr, Vknrow, Tknch, loopnch)
  ! call VLO_withc0_cmplx_epwrap(L, S, J, regtype, Mambda, Cs3s1_600_set2_para,  &
  !   & cmplx(k, 0.0_NER, NEC), cmplx(k, 0.0_NER, NEC), potval)


  ! print *, "VM = ", VM(1:3, 1:3)
  ! print *, "Vk = ", Vk(10, 1)
  ! print *, "flag = ", flag
  ! print *, "Tk = ", Tk(1:10, 1)
  ! print *, "potval = ", potval
  ! print *, "loop = ", loopnch
  ! print *, "potval + loop = ", potval + loopnch

  ! klst(1) = cmplx(k, 0.0_NER, NEC)

  ! nk = 10
  ! klst(1) = cmplx(msh(nk), 0.0_NER, NEC)
  ! klst(1) = f*msh(nk)


end program get_cmplx
