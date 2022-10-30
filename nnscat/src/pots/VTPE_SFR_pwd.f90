! pray  Aug/07/2019
! deltaless TPE potential with SFR in Epelbuam_04 Eur. Phys. J. A 19, 125-137 (2004)
!  DOI: 10.1140/epja/i2003-10096-0
! there are no polynomial and mN-corresponding terms in it, which represents the only difference to Kaiser.
! We used (2.4) and (2.8), replace L(q) with LAq(q) which in (2.22)

module VTPE_SFR_pwd
	use nneft_type,       only: NER, NEC, PI_NE
  use nneft_phyconst,   only: PC_gA, PC_fpi, PC_mpi, PC_delta_cutoff
  use OPE_pwd
  use util_gauleg,   only:feb_glabsc_2pt

  implicit none

  logical, private              :: inited_spc = .false.
  integer, parameter, private   :: spc_N = 1000
  real(NER), private            :: u_ii(1:spc_N), weig_ii(1:spc_N)
!   real(NER),parameter,private   :: func_factor = PC_default_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)

contains

  subroutine init_spc()

    real(NER)           :: M, cut

      M = 2.0_NER * PC_mpi
      cut = PC_delta_cutoff
      if (inited_spc) then
      	return
      end if
      call feb_glabsc_2pt(spc_N, u_ii, weig_ii, M, cut)
      inited_spc = .true.

  end subroutine init_spc

! integral

!     subroutine get_V_W_SFR(q, V_T, V_S, W_C)

!     	real(NER),   intent(in)   :: q
!       real(NER),   intent(out)  :: V_T, V_S, W_C
!       real(NER)                 :: regsum1, regsum2, regsum3
!       integer                   :: ii

!       if (.not. inited_spc) then
!       	call init_spc()
!       end if

!       regsum1 = 0.0_NER
!       regsum2 = 0.0_NER
!       regsum3 = 0.0_NER

!       do ii = 1, spc_N
!       	regsum1 = weig_ii(ii)*(2.0_NER/PI_NE*u_ii(ii)/(u_ii(ii)*u_ii(ii)+q*q)*3.0_NER*PC_gA*PC_gA*PC_gA*PC_gA/(128.0_NER &
!       		       & * PI_NE*PC_fpi*PC_fpi*PC_fpi*PC_fpi)*sqrt(u_ii(ii)*u_ii(ii)-4.0_NER*PC_mpi*PC_mpi)/u_ii(ii)) &
!                  & + regsum1
!       	regsum2 = weig_ii(ii)*(2.0_NER/PI_NE*u_ii(ii)/(u_ii(ii)*u_ii(ii)+q*q)*3.0_NER*PC_gA*PC_gA*PC_gA*PC_gA/(128.0_NER &
!       		       & * PI_NE*PC_fpi*PC_fpi*PC_fpi*PC_fpi)*sqrt(u_ii(ii)*u_ii(ii)-4.0_NER*PC_mpi*PC_mpi)/u_ii(ii) &
!                  & * 1.0_NER/(u_ii(ii)*u_ii(ii))) + regsum2
!         regsum3 = weig_ii(ii)*(2.0_NER/PI_NE*u_ii(ii)/(u_ii(ii)*u_ii(ii)+q*q)/768.0_NER/PI_NE/(PC_fpi*PC_fpi*PC_fpi*PC_fpi) &
!                  & * (4.0_NER*PC_mpi*PC_mpi*(5.0_NER*PC_gA*PC_gA*PC_gA*PC_gA-4.0_NER*PC_gA*PC_gA-1.0_NER) &
!                  & - u_ii(ii)*u_ii(ii)*(23.0_NER*PC_gA*PC_gA*PC_gA*PC_gA-10.0_NER*PC_gA*PC_gA-1.0_NER)+48.0_NER* &
!                  & PC_gA*PC_gA*PC_gA*PC_gA*PC_mpi*PC_mpi*PC_mpi*PC_mpi/(4.0_NER*PC_mpi*PC_mpi-u_ii(ii)*u_ii(ii))) &
!                  & *sqrt(u_ii(ii)*u_ii(ii)-4.0_NER*PC_mpi*PC_mpi)/u_ii(ii)) + regsum3
!       end do
!       V_S = regsum1
!       V_T = regsum2
!       W_C = regsum3

!     end subroutine get_V_W_SFR


!Lcq replace Lq in TPE0

  subroutine get_V_W_SFR(q, V_T, V_S, W_C)

    	real(NER),   intent(in)   :: q
      real(NER),   intent(out)  :: V_T, V_S, W_C
      real(NER)                 :: Lcq, w, s,cut

      V_S = 0.0_NER
      V_T = 0.0_NER
      W_C = 0.0_NER
      cut = PC_delta_cutoff
      s = sqrt(cut*cut - 4.0_NER * PC_mpi*PC_mpi)
      w = sqrt(q * q + 4.0_NER *PC_mpi*PC_mpi)
      Lcq = w / (2.0_NER * q)*log((cut*cut*w*w+q*q*s*s+2.0_NER*cut*q*w*s)/(4.0_NER*PC_mpi*PC_mpi*(cut*cut+q*q)))

      W_C = -1.0_NER/(384.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi*PC_fpi*PC_fpi)*Lcq*(4.0_NER*PC_mpi*PC_mpi*(5.0_NER     &
      	  & *PC_gA*PC_gA*PC_gA*PC_gA-4.0_NER*PC_gA*PC_gA-1.0_NER)+q*q*(23.0_NER*PC_gA*PC_gA*PC_gA*PC_gA-10.0_NER &
      	  & *PC_gA*PC_gA-1.0_NER)+48.0_NER*PC_gA*PC_gA*PC_gA*PC_gA*PC_mpi*PC_mpi*PC_mpi*PC_mpi/(4.0_NER*PC_mpi*PC_mpi+q*q))
      V_T = -3.0_NER*PC_gA*PC_gA*PC_gA*PC_gA/(64.0_NER*PI_NE*PI_NE*PC_fpi*PC_fpi*PC_fpi*PC_fpi)*Lcq
      V_S = V_T * (-q*q)

  end subroutine get_V_W_SFR

  subroutine get_U_SFR(j,q,Ut,Us,Uc)
    real(NER)            :: Uc,Us,Ut,V_T, V_S, W_C
    integer , intent(in) :: j
    real(NER),intent(in) :: q



    call get_V_W_SFR(q, V_T, V_S, W_C)
    if (mod(j, 2) .eq. 1) then
       Uc=-3.0_NER*W_C
       Us=V_S
       Ut=V_T
    else
       Uc=W_C
       Us=V_S
       Ut=V_T
    end if
  end subroutine get_U_SFR


  function TPESFR_j0j(j,p,k)

    real(NER)               :: TPESFR_j0j,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, q
    if (j .lt. 0) then
      TPESFR_j0j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if
    regsum = 0.0_NER
    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
        q      = sqrt(p*p+k*K-2.0_NER*p*k*z)
      call get_U_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)*(Uc-3.0_NER*Us-(p*p+k*k- &
      &      2.0_NER*z*p*k)*Ut)
    end do

    TPESFR_j0j=regsum

  end function TPESFR_j0j




   function TPESFR_j1j(j,p,k)
    real(NER)               :: TPESFR_j1j,Ut,Us,Uc,W_C,V_T,V_S
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, q

    if (j .lt. 1) then
      TPESFR_j1j = 0.0_NER
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
    call get_V_W_SFR(q, V_T, V_S, W_C)
    if (mod(j, 2) .eq. 0) then
       Uc=-3.0_NER*W_C
       Us=V_S
       Ut=V_T
    else
       Uc=W_C
       Us=V_S
       Ut=V_T
    end if
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*((lgndr_table_pwo(j+1, ii)+lgndr_table_pwo(j-1, ii))*(-2.0_NER)*p*k*Ut+ &
      &      (Uc+Us+Ut*(psqr+k*k+2.0_NER*p*k*z))*lgndr_table_pwo(j, ii))*wghtslg_pwo(ii)
    end do

    TPESFR_j1j=regsum

  end function TPESFR_j1j




  function TPESFR_jpp(j,p,k)
    real(NER)               :: TPESFR_jpp,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 0) then
      TPESFR_jpp = 0.0_NER
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
      call get_U_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*2.0_NER/(2.0_NER*j+1.0_NER)*Ut*lgndr_table_pwo(j, ii) &
      &      +(Uc+Us+1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j+1, ii))
    end do

    TPESFR_jpp=regsum

  end function TPESFR_jpp


  function TPESFR_jmm(j,p,k)
    real(NER)               :: TPESFR_jmm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 1) then
      TPESFR_jmm = 0.0_NER
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
      call get_U_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*(-2.0_NER)/(2.0_NER*j+1.0_NER)*Ut*lgndr_table_pwo(j, ii) &
      &      +(Uc+Us-1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j-1, ii))
    end do

    TPESFR_jmm=regsum

  end function TPESFR_jmm


  function TPESFR_jpm(j,p,k)
    real(NER)               :: TPESFR_jpm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr , q

    if (j .lt. 1) then
      TPESFR_jpm = 0.0_NER
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
      call get_U_SFR(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*sqrt(j*(j+1.0_NER))/(2.0_NER*j+1.0_NER)*((-4.0_NER)*Ut*p*k* &
      &      lgndr_table_pwo(j, ii)+ksqr*2.0_NER*Ut*lgndr_table_pwo(j-1, ii)+2.0_NER *psqr*Ut*lgndr_table_pwo(j+1, ii))
      TPESFR_jpm=regsum
    end do

   end function TPESFR_jpm

end module VTPE_SFR_pwd

