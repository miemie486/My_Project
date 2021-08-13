
program get_h2wfs

  use eft_tmtrx
  use util_gauleg
  implicit none

  integer, parameter   :: Nmax = 500
  integer   :: N, ii, jj, NPi
  real(NER) :: Mambda, CMRatio, phi
  real(NER) :: msh(Nmax), wght(Nmax), extrapara(NUM_PARA_TWOBODY)
  complex(NEC)  :: f, h2BE_R, h2BE_C, norm
  complex(NEC)  :: Cmsh(Nmax), Cwght(Nmax), Cwf(2*Nmax), wfvals(2*Nmax), plst(Nmax)
  complex(NEC)  :: Rmsh(Nmax), Rwght(Nmax)
  character(len=256) :: tmpstr

  Mambda = 400.0_NER
  CMRatio = 0.5_NER
  phi = 10.0_NER
  ! phi = 0.0_NER
  f = exp(-IMGUNT_NE*phi/180.0*PI_NE)

  call get_command_argument(1, tmpstr)
  read (tmpstr, *) N

  call feb_tanabsc(N, Mambda*CMRatio, msh, wght)
  do ii = 1, N
    Cmsh(ii) = f*msh(ii)
    Cwght(ii) = f*wght(ii)
    Rmsh(ii) = msh(ii)
    Rwght(ii) = wght(ii)
    ! plst(ii) = cmplx(msh(ii), 0.0_NER, NEC)
  end do

  ! call get_h2wfs_cmplx_chengdu(N, Cmsh, Cwght, "chengdu_MMWLY_cmplx_400",      &
  !   & extrapara, h2BE, Cwf)

  NPi = N
  plst(1:N) = Cmsh(1:N)

  ! NPi = 5
  ! plst(1:NPi) = (/0.1, 0.2, 0.5, 1.0, 2.0/)
  ! plst(1:NPi) = (/10.0, 20.0, 50.0, 100.0, 200.0/)

  ! get_generic_h2wfvals_cmplx_chengdu(Nmesh, Cmsh, Cwght, Np, p_lst, &
  !   & pot_name, extrapara, h2BE, wfval_lst)
  call get_generic_h2wfvals_cmplx_chengdu(N, Cmsh, Cwght, NPi, plst,    &
    & "chengdu_MMWLY_cmplx_400", extrapara, h2BE_C, wfvals)

  call get_generic_h2wfvals_cmplx_chengdu(N, Rmsh, Rwght, N, Cmsh,    &
    & "chengdu_MMWLY_cmplx_400", extrapara, h2BE_R, Cwf)

  norm = 0.0_NER
  do ii = 1, N
    norm = norm + (Cwf(ii)**2 + Cwf(ii+N)**2)*Cmsh(ii)*Cmsh(ii)*Cwght(ii)
    ! write (*, '(i4, 1x, e11.3e3, 1x, e11.3e3, 1x, e11.3e3, 1x, e11.3e3, 1x, e11.3e3, 1x, e11.3e3)') ii, real(Cmsh(ii)), aimag(Cmsh(ii)), real(Cwf(ii)), aimag(Cwf(ii)), real(Cwf(ii+N)), aimag(Cwf(ii+N))
  end do

  write (*, '(a)') "real(Cwf), real(wfvals), aimag(Cwf), aimag(wfvals)"
  do ii = 1, N
    write (*, '(a4, i4, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3)')  &
      & '3S1', ii, real(Cwf(ii)), real(wfvals(ii)), aimag(Cwf(ii)), aimag(wfvals(ii))
  end do
  print *
  write (*, '(a)') "real(Cwf), real(wfvals), aimag(Cwf), aimag(wfvals)"
  do ii = 1, N
    write (*, '(a4, i4, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3)')  &
      & '3D1', ii, real(Cwf(ii+N)), real(wfvals(ii+N)), aimag(Cwf(ii+N)), aimag(wfvals(ii+N))
  end do

  ! print *, "plst_Re, plst_Im, wfval_Re, wfval_Im"
  ! do ii = 1, NPi
  !   write (*, '(a4, i4, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3)')  &
  !     & '3S1', ii, real(plst(ii)), aimag(plst(ii)), real(wfvals(ii)), aimag(wfvals(ii))
  ! end do
  ! do ii = 1, NPi
  !   write (*, '(a4, i4, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3, 1x, e13.5e3)')  &
  !     & '3D1', ii,  real(plst(ii)), aimag(plst(ii)),                   &
  !     & real(wfvals(ii + N)), aimag(wfvals(ii + N))
  ! end do

    ! , real(Cwf(ii+N)), aimag(Cwf(ii+N))

  print *, "phi = ", phi
  print *, "Nmesh = ", N
  print *, "h2BE_R = ", h2BE_R
  print *, "h2BE_C = ", h2BE_C
  print *, "norm = ", norm


end program get_h2wfs
