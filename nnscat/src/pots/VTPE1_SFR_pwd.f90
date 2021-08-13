! pray August/11/2019
! NLO TPE potentail with spectral-function regularization(SFR)
! the potential come from Epelbuam_2004(Eur. Phys. J. A 19, 125-137 (2004) DOI 10.1140/epja/i200310096-0),
!  there are no polynomial and mN-corresponding terms in it, which represents the only difference to Kaiser.
! We used (2.4) and (2.8), replace L(q) with LAq(q) which in (2.22)

module VTPE1_SFR_pwd

	use feb_type,         only: standard_error_unit
  use nneft_type,       only: NER, NEC, PI_NE
  use nneft_phyconst,   only: PC_gA, PC_fpi, PC_mpi, PC_delta_cutoff ,strct_cis, &
  &                           PC_Pavon_c1, PC_Pavon_c3, PC_Pavon_c4,             &
  &                           PC_RSE2_c1, PC_RSE2_c3, PC_RSE2_c4,                &
  &                           PC_RSE3_c1, PC_RSE3_c3, PC_RSE3_c4,                &
  &                           PC_RSE4_c1, PC_RSE4_c3, PC_RSE4_c4,                &
  &                           PC_default_c1, PC_default_c3, PC_default_c4,       &
  &                           PC_DELTAF1_c1, PC_DELTAF1_c3, PC_DELTAF1_c4,       &
  &                           PC_DELTAF2_c1, PC_DELTAF2_c3, PC_DELTAF2_c4,       &
  &                           PC_DELTAF3_c1, PC_DELTAF3_c3, PC_DELTAF3_c4,       &
  &                           PC_DELTAF4_c1, PC_DELTAF4_c3, PC_DELTAF4_c4,       &
  &                           PC_DELTAF5_c1, PC_DELTAF5_c3, PC_DELTAF5_c4,       &
  &                           PC_DELTAF6_c1, PC_DELTAF6_c3, PC_DELTAF6_c4,       &
  &                           PC_DELTAF7_c1, PC_DELTAF7_c3, PC_DELTAF7_c4,       &
  &                           PC_DELTAF8_c1, PC_DELTAF8_c3, PC_DELTAF8_c4,       &
  &                           PC_DELTAF9_c1, PC_DELTAF9_c3, PC_DELTAF9_c4,       &
  &                           PC_DELTAF10_c1, PC_DELTAF10_c3, PC_DELTAF10_c4,    &
  &                           PC_DELTAF11_c1, PC_DELTAF11_c3, PC_DELTAF11_c4,    &
  &                           PC_DELTAF12_c1, PC_DELTAF12_c3, PC_DELTAF12_c4
!   use OPE_pwd,          only: ONEOVEREIGHTPI_pwo, initmshlg_pwo, Nlg_pwo,        &
!   & mshlg_inited_pwo, PLtbl => lgndr_table_pwo, mshz => mshlg_pwo,               &
!   & wghtz => wghtslg_pwo
  use OPE_pwd

  implicit none

!   real(NER),private    :: cut = PC_delta_cutoff
  type(strct_cis)      :: cisSFR

contains

  subroutine set_cisSFR_vtpe1(c1, c3, c4)

  	real(NER), intent(in) :: c1, c3, c4

    cisSFR%c1 = c1
    cisSFR%c3 = c3
    cisSFR%c4 = c4

  end subroutine set_cisSFR_vtpe1

  subroutine Select_cisSFR_vtpe1(set_n)

    integer, intent(in) :: set_n

    select case (set_n)
      case (1)
        call set_cisSFR_vtpe1(PC_Pavon_c1, PC_Pavon_c3, PC_Pavon_c4)
      case (2)
        call set_cisSFR_vtpe1(PC_RSE2_c1, PC_RSE2_c3, PC_RSE2_c4)
      case (3)
        call set_cisSFR_vtpe1(PC_RSE3_c1, PC_RSE3_c3, PC_RSE3_c4)
      case (4)
        call set_cisSFR_vtpe1(PC_RSE4_c1, PC_RSE4_c3, PC_RSE4_c4)
      case (5)
        call set_cisSFR_vtpe1(PC_DELTAF1_c1, PC_DELTAF1_c3, PC_DELTAF1_c4)
      case (6)
        call set_cisSFR_vtpe1(PC_DELTAF2_c1, PC_DELTAF2_c3, PC_DELTAF2_c4)
      case (7)
        call set_cisSFR_vtpe1(PC_DELTAF3_c1, PC_DELTAF3_c3, PC_DELTAF3_c4)
      case (8)
        call set_cisSFR_vtpe1(PC_DELTAF4_c1, PC_DELTAF4_c3, PC_DELTAF4_c4)
      case (9)
        call set_cisSFR_vtpe1(PC_DELTAF5_c1, PC_DELTAF5_c3, PC_DELTAF5_c4)
      case (10)
        call set_cisSFR_vtpe1(PC_DELTAF6_c1, PC_DELTAF6_c3, PC_DELTAF6_c4)
      case (11)
        call set_cisSFR_vtpe1(PC_DELTAF7_c1, PC_DELTAF7_c3, PC_DELTAF7_c4)
      case (12)
        call set_cisSFR_vtpe1(PC_DELTAF8_c1, PC_DELTAF8_c3, PC_DELTAF8_c4)
      case (13)
        call set_cisSFR_vtpe1(PC_DELTAF9_c1, PC_DELTAF9_c3, PC_DELTAF9_c4)
      case (14)
        call set_cisSFR_vtpe1(PC_DELTAF10_c1, PC_DELTAF10_c3, PC_DELTAF10_c4)
      case (15)
        call set_cisSFR_vtpe1(PC_DELTAF11_c1, PC_DELTAF11_c3, PC_DELTAF11_c4)
      case (16)
        call set_cisSFR_vtpe1(PC_DELTAF12_c1, PC_DELTAF12_c3, PC_DELTAF12_c4)
      case default
        write (standard_error_unit, '(a)')  &
        & 'Select_cisSFR_vtpe1: unknown id of ci set. Default set pre-defined will be used'
    end select

  end subroutine Select_cisSFR_vtpe1

  subroutine get_cisSFR_vtpe1(cisout)

    type (strct_cis), intent(out) :: cisout

    cisout%c1 = cisSFR%c1
    cisout%c3 = cisSFR%c3
    cisout%c4 = cisSFR%c4

  end subroutine get_cisSFR_vtpe1

  subroutine get_W_V_SFR(q, Vc, Wt, Ws)

  	real(NER), intent(in)  :: q
  	real(NER), intent(out) :: Vc, Ws, Wt
  	real(NER)              :: Acq, cut, qsqr, ga_4th, mpsqr, fpi_4th, c1, c3, c4

  	cut = PC_delta_cutoff
  	qsqr = q*q
  	mpsqr = PC_mpi*PC_mpi
  	ga_4th = PC_gA*PC_gA*PC_gA*PC_gA
  	fpi_4th = PC_fpi*PC_fpi*PC_fpi*PC_fpi
  	c1 = cisSFR%c1
  	c3 = cisSFR%c3
  	c4 = cisSFR%c4
!debug
!  print *, "sfrTPE=",c1, c3,c4
  	Vc = 0.0_NER
  	Wt = 0.0_NER
  	Ws = 0.0_NER
  	Acq =  1.0_NER/(2.0_NER*q)*atan(q*(cut-2.0_NER*PC_mpi)/(qsqr+2.0_NER*cut*PC_mpi))

  	Vc = -3.0_NER*PC_gA*PC_gA/16.0_NER/PI_NE/fpi_4th*(2.0_NER*mpsqr*(2.0_NER*c1-c3)-c3*qsqr) &
  	&    *(2.0_NER*mpsqr+qsqr)*Acq
  	Wt = -PC_gA*PC_gA*c4/32.0_NER/PI_NE/fpi_4th*(4.0_NER*mpsqr+qsqr)*Acq
  	Ws = wt * (-qsqr)

  end subroutine get_W_V_SFR

  subroutine get_U1_SFR(j, q, Ut, Us, Uc)

  	real(NER), intent(out)    :: Uc, Us, Ut
  	real(NER), intent(in)     :: q
  	integer, intent(in)       :: j
  	real(NER)                 :: Vc, Wt, Ws

    call get_W_V_SFR(q, Vc, Wt, Ws)
    if (mod(j, 2) .eq. 1) then
    	Uc = Vc
    	Us = -3.0*Ws
    	Ut = -3.0*wt
    else
    	Uc = Vc
    	Us = Ws
    	Ut = Wt
    end if
  end subroutine get_U1_SFR

  function TPESFR1_j0j(j,p,k)

    real(NER)               :: TPESFR1_j0j,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, q
    if (j .lt. 0) then
      TPESFR1_j0j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if
    regsum = 0.0_NER
    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
        q      = sqrt(p*p+k*K-2.0_NER*p*k*z)
      call get_U1_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)*(Uc-3.0_NER*Us-(p*p+k*k- &
      &      2.0_NER*z*p*k)*Ut)
    end do

    TPESFR1_j0j=regsum

  end function TPESFR1_j0j




   function TPESFR1_j1j(j,p,k)
    real(NER)               :: TPESFR1_j1j,Ut,Us,Uc, Vc, Wt, Ws
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, q

    if (j .lt. 1) then
      TPESFR1_j1j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr=p*p

    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
      q      = sqrt(p*p+k*k-2.0_NER*p*k*z)
    call get_W_V_SFR(q, Vc, Wt, Ws)
    if (mod(j, 2) .eq. 0) then
    	Uc = Vc
    	Us = -3.0*Ws
    	Ut = -3.0*wt
    else
    	Uc = Vc
    	Us = Ws
    	Ut = Wt
    end if
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*((lgndr_table_pwo(j+1, ii)+lgndr_table_pwo(j-1, ii))*(-2.0_NER)*p*k*Ut+ &
      &      (Uc+Us+Ut*(psqr+k*k+2.0_NER*p*k*z))*lgndr_table_pwo(j, ii))*wghtslg_pwo(ii)
    end do

    TPESFR1_j1j=regsum

  end function TPESFR1_j1j




  function TPESFR1_jpp(j,p,k)
    real(NER)               :: TPESFR1_jpp,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 0) then
      TPESFR1_jpp = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr=p*p
    ksqr=k*k

    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
        q      = sqrt(p*p+k*K-2.0_NER*p*k*z)
      call get_U1_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*2.0_NER/(2.0_NER*j+1.0_NER)*Ut*lgndr_table_pwo(j, ii) &
      &      +(Uc+Us+1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j+1, ii))
    end do

    TPESFR1_jpp=regsum

  end function TPESFR1_jpp


  function TPESFR1_jmm(j,p,k)
    real(NER)               :: TPESFR1_jmm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 1) then
      TPESFR1_jmm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr=p*p
    ksqr=k*k

    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
        q      = sqrt(p*p+k*K-2.0_NER*p*k*z)
      call get_U1_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*(-2.0_NER)/(2.0_NER*j+1.0_NER)*Ut*lgndr_table_pwo(j, ii) &
      &      +(Uc+Us-1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j-1, ii))
    end do

    TPESFR1_jmm=regsum

  end function TPESFR1_jmm


  function TPESFR1_jpm(j,p,k)
    real(NER)               :: TPESFR1_jpm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr , q

    if (j .lt. 1) then
      TPESFR1_jpm = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    regsum = 0.0_NER
    psqr=p*p
    ksqr=k*k

    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
        q      = sqrt(p*p+k*K-2.0_NER*p*k*z)
      call get_U1_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*sqrt(j*(j+1.0_NER))/(2.0_NER*j+1.0_NER)*((-4.0_NER)*Ut*p*k* &
      &      lgndr_table_pwo(j, ii)+ksqr*2.0_NER*Ut*lgndr_table_pwo(j-1, ii)+2.0_NER *psqr*Ut*lgndr_table_pwo(j+1, ii))
      TPESFR1_jpm=regsum
    end do

   end function TPESFR1_jpm





end module VTPE1_SFR_pwd
