!get swave potential  pray 2019/12/07

program get_swave_potential
	use chengdu
  use potwrap
	implicit none
	integer, parameter  :: N = 100
  real, parameter     :: pr_hbarc = 197.3286_NER
	integer         :: L, S, J, ii, jj
	real(NER)       :: mambda, p1(1:N), p2(1:N), para(1:2), max_p


  real(NER), allocatable    :: pout(:), pin(:), weig(:), potval(:,:,:,:)
  integer                   :: uptoQn, cmdcnt, kk
  real(NER)                 :: K_MAX(6)
  real(NER)                 :: v(1:2,1:2)
  integer                   :: ios, ios1
  character(len = 1024)     :: cmdbuff, tmp
  character(len = 20)       :: pot(6), cLL, cSS, cJJ, cuptoQn, cK_MAX, Pversion, chnstr
  character(len = 64)      :: fname
  type(lsj_symbol)         :: lsj
  logical                  :: succ_flag



!   max_p = 500.0_NER
  cmdbuff = ''
  cmdcnt = command_argument_count()
  do ii = 1, cmdcnt
    call get_command_argument(ii, tmp)
    cmdbuff = trim(cmdbuff)//' '//trim(tmp)
  end do
  read (cmdbuff, *, IOSTAT=ios) chnstr, uptoQn
  call convert_text_to_lsj(trim(chnstr), lsj, succ_flag)
!   read (cmdbuff, *, IOSTAT=ios) pot, N, L, S, J, uptoQn, K_MAX
  K_MAX(1:6)=(/350.0_NER,400.0_NER,600.0_NER,800.0_NER,1600.0_NER,3200.0_NER/)
  pot(1:6)=(/'chengdu_MMWLY_350 ','chengdu_MMWLY_400 ','chengdu_MMWLY_600 ','chengdu_MMWLY_800 ','chengdu_MMWLY_1600','chengdu_MMWLY_3200'/)

  if  (ios /= 0) then
    write (standard_error_unit, '(a)')  &
      & 'get_swave_data: Something wrong with the arguments'
    write (standard_error_unit, '(a)')  &
      & 'Usage: get_swave_data chnstr uptoQn'
    stop
  end if

  allocate(pout(1:N), pin(1:N), weig(1:N), potval(1:N,1:N,1:2,1:2))

  write (cLL, '(i1.1)') int(lsj%L)
  write (cSS, '(i1.1)') int(lsj%S)
  write (cJJ, '(i1.1)') int(lsj%J)
  write (cuptoQn, '(i1.1)') uptoQn
  call read_chengdu_version(Pversion)

  do kk = 1,6
    fname = trim(pot(kk))//'_'//trim(chnstr)//'_'//trim(cuptoQn)//'.mtx'
    ios1 = int(57)+ kk
    open(unit=ios1,file=fname)
    pot(kk) = trim(pot(kk))
!     call feb_glabsc_2pt(N, pin, weig, 0.0_NER, K_MAX(kk))
!     pout(1:N) = pin(1:N)
    write (cK_MAX, '(i4.4)') int( K_MAX(kk))
    do ii = 1, N
      max_p = K_MAX(kk)
      pin(ii) = ii * max_p/N
      do jj = 1, N
        pout(jj) = jj * max_p/N
        call chengdu_dispatch(pot(kk))
        ! secondly, call the generic routine
        call chengdu_hotpot(lsj%L, lsj%S, lsj%J, uptoQn, pout(jj), pin(ii), v)
!       call chengdu_DLSPR_350(L, S, J, uptoQn, pout(ii), pin(jj), v)
        potval(ii,jj,1:2,1:2) = v(1:2,1:2)*pr_hbarc
      end do
    end do
    if (lsj%L /= lsj%j .and. lsj%j /= 0) then  ! is coupled channel?
      write(ios1,'(a)')"# coupled"
      write(ios1,'(a)')"# "//trim(Pversion)
      write(ios1,'(a)')"# Potential name = "//pot(kk)
      write(ios1,'(a)')"# L = "//trim(cll)
      write(ios1,'(a)')"# S = "//trim(cSS)
      write(ios1,'(a)')"# J = "//trim(cJJ)
      write(ios1,'(a)')"# Order = "//cuptoQn
      write(ios1,'(a)')"# K_MAX = "//trim(cK_MAX)//" MeV"
      write(ios1,'(a)')"# For coupled channels, e.g., 3S1 - 3D1 (LSJ = 011)"
      write(ios1,'(a)')"# Normalization"
      write(ios1,'(a)')"# If V were so weak as to validate the Born approximation, the relation of V"
      write(ios1,'(a)')"# in 1S0 and the respective scattering length 'a' would have been"
      write(ios1,'(a)')"# <p|V|p> = 2/pi * a + ... for p->0, in unit of MeV^(-1)"
      write(ios1,'(a)')"# Momenta are in unit of MeV"
      write(ios1,'(a)')"# pin  pout  v(1, 1)  v(1, 2)  v(2, 2)"
!       write(ios1,'(a)')"# The value of pin and pout, where pin(1:N) = pout(1:N)"
!       write(ios1,*)pin
!       write(ios1,'(a)')"# Matrix elements start here "
      do ii = 1, N
        do jj =1, N
          write(ios1,*)pin(ii), pout(jj),potval(ii,jj,1,1),potval(ii,jj,1,2),potval(ii,jj,2,2)
        end do
      end do
    else
      write(ios1,'(a)')"# uncoupled"
      write(ios1,'(a)')"# "//trim(Pversion)
      write(ios1,'(a)')"# Potential name = "//pot(kk)
      write(ios1,'(a)')"# L = "//trim(cll)
      write(ios1,'(a)')"# S = "//trim(cSS)
      write(ios1,'(a)')"# J = "//trim(cJJ)
      write(ios1,'(a)')"# Order = "//cuptoQn
      write(ios1,'(a)')"# K_MAX = "//trim(cK_MAX)//" MeV"
      write(ios1,'(a)')"# For uncoupled channels, e.g., 3P0 (LSJ = 110)"
      write(ios1,'(a)')"# Normalization"
      write(ios1,'(a)')"# If V were so weak as to validate the Born approximation, the relation of V"
      write(ios1,'(a)')"# in 1S0 and the respective scattering length 'a' would have been"
      write(ios1,'(a)')"# <p|V|p> = 2/pi * a + ... for p->0, in unit of MeV^(-1)"
      write(ios1,'(a)')"# Momenta are in unit of MeV"
      write(ios1,'(a)')"# pin  pout  v(1, 1)"
!       write(ios1,'(a)')"# The value of pin and pout, where pin(1:N) = pout(1:N)"
!       write(ios1,*)pin
!       write(ios1,'(a)')"# Matrix elements start here "
      do ii = 1, N
        do jj = 1, N
          write(ios1,*) pin(ii), pout(jj), potval(ii,jj,1,1)
        end do
      end do
    end if
      close(ios1)
! !       debug
    print *, pot(kk)
    print *, fname
!     print *, K_MAX(kk)
  end do
!   character(len = 64)      :: fname
! ! 3s1
!   L = 0
!   S = 1
!   J = 1
! !   mambda = 350.0_NER
! !   para(1) = -5.4945264707960580E-003_NER
!    max_p = 800.0_NER


!   open(unit=257,file = "3s1_potential")
!   open(unit=157,file = "1s0_potential")
!   open(unit=107,file = "3p0_potential")
!   open(unit=57, file = "3s1_onshell")
!   open(unit=27, file = "1s0_onshell")
!   open(unit=7, file = "3p0_onshell")

!   do ii = 1, N
!   	p1(ii) = ii * max_p/N
!   	do jj = 1, N
!       p2(jj) = jj * max_p/N
!       call chengdu_dlspr_800(L, S, J, 0, p1(ii), p2(jj), potval)
!       s31(ii, jj) = potval(1, 1)
!       s311(ii, jj) = potval(1, 2)
!       s312(ii, jj) = potval(2, 1)
!       s313(ii, jj) = potval(2, 2)

!   		write(257,*) p1(ii), p2(jj), s31(ii, jj)*197.0_NER
!   		if (ii == jj) then
!   			write(57,*) p1(ii), s31(ii,jj)*197.0_NER, s311(ii, jj)*197.0_NER, s312(ii, jj)*197.0_NER, s313(ii, jj)*197.0_NER
!   		end if
!     end do
!   end do

! ! 1s0
!   L = 0
!   S = 0
!   J = 0
! !   mambda = 350.0_NER
! !   para(1) = 0.42647716568097287E+002_NER
! !   para(2) = -0.15817891988985733E+000_NER
! !   max_p = 350.0_NER

!   do ii = 1, N
!     p1(ii) = ii * max_p/N
!     do jj = 1, N
!       p2(jj) = jj * max_p/N
!       call chengdu_dlspr_800(L, S, J, 0, p1(ii), p2(jj), potval)
!       s10(ii, jj) = potval(1, 1)
!       write(157,*) p1(ii), p2(jj), s10(ii, jj)*197.0_NER
!       if (ii == jj) then
!         write(27,*) p1(ii), s10(ii,jj)*197.0_NER
!       end if
!     end do
!   end do

! ! 3p0
!   L = 1
!   S = 1
!   J = 0
! !   mambda = 350.0_NER
! !   para(1) = 0.42647716568097287E+002_NER
! !   para(2) = -0.15817891988985733E+000_NER
! !   max_p = 350.0_NER

!   do ii = 1, N
!     p1(ii) = ii * max_p/N
!     do jj = 1, N
!       p2(jj) = jj * max_p/N
!       call chengdu_dlspr_800(L, S, J, 0, p1(ii), p2(jj), potval)
!       p30(ii, jj) = potval(1, 1)
!       write(107,*) p1(ii), p2(jj), p30(ii, jj)*197.0_NER
!       if (ii == jj) then
!         write(7,*) p1(ii), p30(ii,jj)*197.0_NER
!       end if
!     end do
!   end do

!   close(257)
!   close(157)
!   close(107)
!   close(57)
!   close(27)
!   close(7)
end program