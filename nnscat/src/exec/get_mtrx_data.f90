! pray 06/05/2020
! getting potential datas

program get_mtrx_data

  use nneft_type
  use chengdu
  use util_gauleg, only : feb_glabsc_2pt
  implicit none

  integer, parameter        :: N_mambda = 34
  real(NER), allocatable    :: pout(:), pin(:), weig(:), potval(:,:,:,:)
  integer                   :: uptoQn, L, S, J, N, cmdcnt, kk
  real(NER)                 :: K_MAX(N_mambda)
  real(NER)                 :: v(1:2,1:2)
  integer                   :: ii, jj, ios, ios1
  character(len = 1024)     :: cmdbuff, tmp
  character(len = 20)       :: pot(N_mambda), cLL, cSS, cJJ, cuptoQn, cK_MAX, Pversion, cNmesh
  character(len = 64)      :: fname
!   character(len = 1024)

!   read (*, IOSTAT=ios) N, L, S, J, uptoQn, K_MAX

  cmdbuff = ''
  cmdcnt = command_argument_count()
  do ii = 1, cmdcnt
    call get_command_argument(ii, tmp)
    cmdbuff = trim(cmdbuff)//' '//trim(tmp)
  end do
  read (cmdbuff, *, IOSTAT=ios) N, L, S, J, uptoQn
!   read (cmdbuff, *, IOSTAT=ios) pot, N, L, S, J, uptoQn, K_MAX
!   K_MAX(1:10)=(/700.0_NER,800.0_NER,1200.0_NER,1600.0_NER,3200.0_NER,6400.0_NER,2400.0_NER,4000.0_NER,4800.0_NER,5600.0_NER/)
  K_MAX(1:N_mambda)=(/350.0_NER,400.0_NER,600.0_NER,800.0_NER,1600.0_NER,3200.0_NER,1200.0_NER,2000.0_NER,2400.0_NER,2800.0_NER, &
    & 500.0_NER,1800.0_NER,2200.0_NER,4000.0_NER,4800.0_NER, 900.0_NER, 1300.0_NER, 1400.0_NER, 6400.0_NER,     &
    & 1000.0_NER, 1200.0_NER, 1300.0_NER, 1400.0_NER, 1600.0_NER, 2000.0_NER, 2200.0_NER, 2400.0_NER, 800.0_NER, 2800.0_NER, 3200.0_NER, &
    & 4000.0_NER, 4800.0_NER, 400.0_NER, 600.0_NER/)
  pot(1:N_mambda)=(/'chengdu_MMWLY_350 ','chengdu_MMWLY_400 ','chengdu_MMWLY_600 ','chengdu_MMWLY_800 ','chengdu_MMWLY_1600','chengdu_MMWLY_3200'&
    & ,'chengdu_MMWLY_1200','chengdu_MMWLY_2000','chengdu_MMWLY_2400','chengdu_MMWLY_2800','chengdu_MMWLY_500 ','chengdu_MMWLY_1800', &
    & 'chengdu_MMWLY_2200','chengdu_MMWLY_4000','chengdu_MMWLY_4800','chengdu_MMWLY_900 ','chengdu_MMWLY_1300','chengdu_MMWLY_1400', &
    & 'chengdu_MMWLY_6400','chengdu_WOTPE_1000','chengdu_WOTPE_1200','chengdu_WOTPE_1300','chengdu_WOTPE_1400','chengdu_WOTPE_1600','chengdu_WOTPE_2000', &
    & 'chengdu_WOTPE_2200','chengdu_WOTPE_2400','chengdu_WOTPE_800 ','chengdu_WOTPE_2800','chengdu_WOTPE_3200','chengdu_WOTPE_4000','chengdu_WOTPE_4800', &
    & 'chengdu_WOTPE_400 ','chengdu_WOTPE_600 '/)

  if  (ios /= 0) then
    write (standard_error_unit, '(a)')  &
      & 'get_mtrx_data: Something wrong with the arguments'
    write (standard_error_unit, '(a)')  &
      & 'Usage: get_mtrx_data N L S J uptoQn'
    stop
  end if

  allocate(pout(1:N), pin(1:N), weig(1:N), potval(1:N,1:N,1:2,1:2))

  write (cLL, '(i1.1)') int(L)
  write (cSS, '(i1.1)') int(S)
  write (cJJ, '(i1.1)') int(J)
  write (cuptoQn, '(i1.1)') uptoQn
  call read_chengdu_version(Pversion)

  do kk = 1,N_mambda
    fname = trim(pot(kk))//'_'//trim(cLL)//trim(cSS)//trim(cJJ)//'_'//trim(cuptoQn)//'.mtx'
    ios1 = int(57)+ kk
    open(unit=ios1,file=fname)
    pot(kk) = trim(pot(kk))
    K_MAX(kk) = 3.0_NER*K_MAX(kk)
    call feb_glabsc_2pt(N, pin, weig, 0.0_NER, K_MAX(kk))
    pout(1:N) = pin(1:N)
    write (cNmesh, '(i3.3)') int(N)
    if(K_MAX(kk) <= 9999.0) then
      write (cK_MAX, '(i4.4)') int(K_MAX(kk))
    else
      write (cK_MAX, '(i5.5)') int(K_MAX(kk))
    end if
    do ii = 1,N
      do jj = 1,N
        ! Or, firstly set the name of the potential,
        call chengdu_dispatch(pot(kk))
        ! secondly, call the generic routine
        call chengdu_hotpot(L, S, J, uptoQn, pout(jj), pin(ii), v)
!       call chengdu_DLSPR_350(L, S, J, uptoQn, pout(ii), pin(jj), v)
        potval(ii,jj,1:2,1:2) = v(1:2,1:2)
      end do
    end do
    if (L /= j .and. j /= 0) then  ! is coupled channel?
      write(ios1,'(a)')"# "//trim(Pversion)
      write(ios1,'(a)')"# Potential name = "//pot(kk)
      write(ios1,'(a)')"# Nmesh = "//trim(cNmesh)
      write(ios1,'(a)')"# L = "//trim(cll)
      write(ios1,'(a)')"# S = "//trim(cSS)
      write(ios1,'(a)')"# J = "//trim(cJJ)
      write(ios1,'(a)')"# Order = "//cuptoQn
      write(ios1,'(a)')"# K_MAX = "//trim(cK_MAX)//" MeV"
      write(ios1,'(a)')"# For coupled channels, e.g., 3S1 - 3D1 (LSJ = 011)"
      write(ios1,'(a)')"# storage = D(1:2×N，1:2×N) (row, column)"
      write(ios1,'(a)')"# D(i,j)     =  <pout(j),L|V|L,pin(i)>"
      write(ios1,'(a)')"# D(i,j+N)   =  <pout(j),L|V|L+2,pin(i)>"
      write(ios1,'(a)')"# D(i+N,j)   =  <pout(j),L+2|V|L,pin(i)>"
      write(ios1,'(a)')"# D(i+N,j+N) =  <pout(j),L+2|V|L+2,pin(i)>"
      write(ios1,'(a)')"# Normalization"
      write(ios1,'(a)')"# If V were so weak as to validate the Born approximation, the relation of V"
      write(ios1,'(a)')"# in 1S0 and the respective scattering length 'a' would have been"
      write(ios1,'(a)')"# <p|V|p> = 2/pi * a + ... for p->0, in unit of MeV^(-1)"
      write(ios1,'(a)')"# Momenta are in unit of MeV"
      write(ios1,'(a)')"# The value of pin and pout, where pin(1:N) = pout(1:N)"
      write(ios1,*)pin
      write(ios1,'(a)')"# The weights of pin and pout"
      write(ios1,*)weig
      write(ios1,'(a)')"# Matrix elements start here "
      do ii = 1, N
        write(ios1,*)potval(ii,1:N,1,1),potval(ii,1:N,1,2)
      end do
      do ii = 1, N
        write(ios1,*)potval(ii,1:N,2,1),potval(ii,1:N,2,2)
      end do
    else
      write(ios1,'(a)')"# "//trim(Pversion)
      write(ios1,'(a)')"# Potential name = "//pot(kk)
      write(ios1,'(a)')"# Nmesh = "//trim(cNmesh)
      write(ios1,'(a)')"# L = "//trim(cll)
      write(ios1,'(a)')"# S = "//trim(cSS)
      write(ios1,'(a)')"# J = "//trim(cJJ)
      write(ios1,'(a)')"# Order = "//cuptoQn
      write(ios1,'(a)')"# K_MAX = "//trim(cK_MAX)//" MeV"
      write(ios1,'(a)')"# For uncoupled channels, e.g., 3P0 (LSJ = 110)"
      write(ios1,'(a)')"# storage = D(1:N，1:N) (row, column)"
      write(ios1,'(a)')"# D(i,j)  = <pout(j),L|V|L,pin(i)>"
      write(ios1,'(a)')"# Normalization"
      write(ios1,'(a)')"# If V were so weak as to validate the Born approximation, the relation of V"
      write(ios1,'(a)')"# in 1S0 and the respective scattering length 'a' would have been"
      write(ios1,'(a)')"# <p|V|p> = 2/pi * a + ... for p->0, in unit of MeV^(-1)"
      write(ios1,'(a)')"# Momenta are in unit of MeV"
      write(ios1,'(a)')"# The value of pin and pout, where pin(1:N) = pout(1:N)"
      write(ios1,*)pin
      write(ios1,'(a)')"# The weights of pin and pout"
      write(ios1,*)weig
      write(ios1,'(a)')"# Matrix elements start here "
      do ii = 1, N
        write(ios1,*)potval(ii,1:N,1,1)
      end do
    end if
      close(ios1)
! !       debug
    print *, pot(kk)
    print *, fname
    print *, K_MAX(kk)
  end do!   debug
!   ii = 5
!   jj = ii +1
! !   print *, "pin=", pin(ii), "pout=", pout(ii+1), "potential=", potval(ii,ii+1,1:2,1:2)
!   write (cLL, '(i1.1)') int(L)
!   write (cSS, '(i1.1)') int(S)
!   write (cJJ, '(i1.1)') int(J)
!   write (cuptoQn, '(i1.1)') uptoQn
!   write (cK_MAX, '(i4.4)') int(K_MAX)


!       ios1 = int(57)+ uptoQn
!       open(unit=ios1,file=fname)
! !       open(unit=ios1,file="test.out")
! !       debug
! print *, fname

!     if (L /= j .and. j /= 0) then  ! is coupled channel?
!       write(ios1,'(a)')"# "//trim(Pversion)
!       write(ios1,'(a)')"# Potential name = "//pot
!       write(ios1,'(a)')"# L = "//trim(cll)
!       write(ios1,'(a)')"# S = "//trim(cSS)
!       write(ios1,'(a)')"# J = "//trim(cJJ)
!       write(ios1,'(a)')"# Order = "//cuptoQn
!       write(ios1,'(a)')"# K_MAX = "//trim(cK_MAX)//" MeV"
!       write(ios1,'(a)')"# For coupled channels, e.g., 3S1 - 3D1 (LSJ = 011)"
!       write(ios1,'(a)')"# storage = D(1:2×N，1:2×N) (row, column)"
!       write(ios1,'(a)')"# D(i,j)     =  <pout(j),L|V|L,pin(i)>"
!       write(ios1,'(a)')"# D(i,j+N)   =  <pout(j),L|V|L+2,pin(i)>"
!       write(ios1,'(a)')"# D(i+N,j)   =  <pout(j),L+2|V|L,pin(i)>"
!       write(ios1,'(a)')"# D(i+N,j+N) =  <pout(j),L+2|V|L+2,pin(i)>"
!       write(ios1,'(a)')"# Normalization"
!       write(ios1,'(a)')"# If V were so weak as to validate the Born approximation, the relation of V"
!       write(ios1,'(a)')"# in 1S0 and the respective scattering length 'a' would have been"
!       write(ios1,'(a)')"# <p|V|p> = 2/pi * a + ... for p->0, in unit of MeV^(-1)"
!       write(ios1,'(a)')"# Momenta are in unit of MeV"
!       write(ios1,'(a)')"# The value of pin and pout, where pin(1:N) = pout(1:N)"
!       write(ios1,*)pin
!       write(ios1,'(a)')"# Matrix elements start here "
!       do ii = 1, N
!         write(ios1,*)potval(ii,1:N,1,1),potval(ii,1:N,1,2)
!       end do
!       do ii = 1, N
!         write(ios1,*)potval(ii,1:N,2,1),potval(ii,1:N,2,2)
!       end do
!     else
!       write(ios1,'(a)')"# "//trim(Pversion)
!       write(ios1,'(a)')"# Potential name = "//pot
!       write(ios1,'(a)')"# L = "//trim(cll)
!       write(ios1,'(a)')"# S = "//trim(cSS)
!       write(ios1,'(a)')"# J = "//trim(cJJ)
!       write(ios1,'(a)')"# Order = "//cuptoQn
!       write(ios1,'(a)')"# K_MAX = "//trim(cK_MAX)//" MeV"
!       write(ios1,'(a)')"# For uncoupled channels, e.g., 3P0 (LSJ = 110)"
!       write(ios1,'(a)')"# storage = D(1:N，1:N) (row, column)"
!       write(ios1,'(a)')"# D(i,j)  = <pout(j),L|V|L,pin(i)>"
!       write(ios1,'(a)')"# Normalization"
!       write(ios1,'(a)')"# If V were so weak as to validate the Born approximation, the relation of V"
!       write(ios1,'(a)')"# in 1S0 and the respective scattering length 'a' would have been"
!       write(ios1,'(a)')"# <p|V|p> = 2/pi * a + ... for p->0, in unit of MeV^(-1)"
!       write(ios1,'(a)')"# Momenta are in unit of MeV"
!       write(ios1,'(a)')"# The value of pin and pout, where pin(1:N) = pout(1:N)"
!       write(ios1,*)pin
!       write(ios1,'(a)')"# Matrix elements start here "
!       do ii = 1, N
!         write(ios1,*)potval(ii,1:N,1,1)
!       end do
!     end if
!       close(ios1)
    !debug
!     ii = 5
!     jj = ii +1
!     print *, "potential(ii,jj)=", potval(ii,jj, 1,1), potval(ii,jj, 1,2), potval(ii,jj, 2,1), potval(ii,jj, 2,2)
!     print *, "potential(ii,jj)=", potval(ii,jj, 1:2,1:2)

end program get_mtrx_data
