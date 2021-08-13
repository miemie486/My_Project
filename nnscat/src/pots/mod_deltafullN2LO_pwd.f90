! Bingwei Long 07/16/2018

! Subleading Deltaful TPE
module mod_deltafullN2LO_pwd

  ! Put behind "only:" the name of variables, functions and subs your code is
  ! going to refer to in this module
  use nneft_type,       only: NER, NEC, PI_NE
  ! Check nneft_phyconst for name of common physics parameters.

  use nneft_phyconst,   only: PC_mN, PC_gA, PC_fpi, PC_mpi, Delta_paras,PC_delta,PC_hA1,      &
  &                         PC_DELTAF1_c1, PC_DELTAF1_c2, PC_DELTAF1_c3, PC_DELTAF1_c4,       &
  &                         PC_DELTAF2_c1, PC_DELTAF2_c2, PC_DELTAF2_c3, PC_DELTAF2_c4,       &
  &                         PC_DELTAF3_c1, PC_DELTAF3_c2, PC_DELTAF3_c3, PC_DELTAF3_c4,       &
  &                         PC_DELTAF4_c1, PC_DELTAF4_c2, PC_DELTAF4_c3, PC_DELTAF4_c4,       &
  &                         PC_DELTAF5_c1, PC_DELTAF5_c2, PC_DELTAF5_c3, PC_DELTAF5_c4,       &
  &                         PC_DELTAF6_c1, PC_DELTAF6_c2, PC_DELTAF6_c3, PC_DELTAF6_c4,       &
  &                         PC_DELTAF7_c1, PC_DELTAF7_c2, PC_DELTAF7_c3, PC_DELTAF7_c4,       &
  &                         PC_DELTAF8_c1, PC_DELTAF8_c2, PC_DELTAF8_c3, PC_DELTAF8_c4,       &
  &                         PC_DELTAF9_c1, PC_DELTAF9_c2, PC_DELTAF9_c3, PC_DELTAF9_c4
  use ope_pwd
  use util_gauleg,   only:feb_glabsc_2pt
  use mod_deltafullNLO_pwd,   only:Dc_q , get_LAH, dps




  implicit none

  ! real(NER) ,parameter, private ::  del=  292.0_NER    ,  PC_delta_cutoff=  1200.0_NER
!    real(NER), parameter, private ::  PC_delta_cutoff= PC_delta_PC_delta_cutoffoff

   ! integer              :: i=5

   ! real(NER), dimension(1:6)   :: c_1, c_2 ,c_3, c_4

  ! data c_1/1.,2.,3.,4.,5.,6./


   ! these parameters from "reconciling threshold and subthreshold expansions for pion-nucleon scattering"
   ! arXiv:1610.08978v2 [nucl-th] 27 Apr 2017
  ! 1 2 3 for NLO at N-N PI-N and covariant respectively,  4 5 6 for N2LO at N-N PI-N and covariant respectively

    ! data c_1 /-0.74_NER,-0.69_NER,-0.69_NER,-1.25_NER,-1.24_NER,-1.12_NER/
    ! data c_2 /-0.49_NER,0.81_NER,0.40_NER,1.37_NER,0.79_NER,1.02_NER/
    ! data C_3 /-0.65_NER,-0.44_NER,-0.49_NER,-2.41_NER,-2.49_NER,-2.27_NER/
    ! data C_4 /0.96_NER,0.64_NER,0.64_NER,1.66_NER,1.67_NER,1.21_NER/

!      type (Delta_paras)  :: dps



  ! logical             :: mshlg_inited_pwo = .false.
  ! integer, parameter  :: Nlg_pwo = 100
  ! real(NER) :: mshlg_pwo(1:Nlg_pwo), wghtslg_pwo(1:Nlg_pwo)
  ! integer, parameter  :: lgrank_pwo = 20
  ! real(NER) :: lgndr_table_pwo(0:lgrank_pwo, 1:Nlg_pwo)



  ! Define your own physics parameters here if you do not find them in
  ! nneft_phyconst

contains

  subroutine get_deltaparas_N2LO(para)

    type(Delta_paras), intent(out) :: para

    para%c1 = dps%c1
    para%c2 = dps%c2
    para%c3 = dps%c3
    para%c4 = dps%c4
    para%ha = dps%ha
    para%del= dps%del

  end subroutine get_deltaparas_N2LO

  ! Add below your function that generate partial-wave matrix elements
  ! For example:


    subroutine set_cis2PE_vdel(c1, c2, c3, c4 , ha, del)

    real(NER), intent(in) :: c1, c2, c3, c4 , ha ,del

    dps%c1 = c1
    dps%c2 = c2
    dps%c3 = c3
    dps%c4 = c4
    dps%ha = ha
    dps%del= del

    end subroutine set_cis2PE_vdel

    subroutine Select_cis2PE_vdel(set_n)

    integer, intent(in) :: set_n

    select case (set_n)
      case (5)
        call set_cis2PE_vdel(PC_DELTAF1_c1, PC_DELTAF1_c2, PC_DELTAF1_c3, PC_DELTAF1_c4,PC_hA1,PC_delta)
      case (6)
        call set_cis2PE_vdel(PC_DELTAF2_c1, PC_DELTAF2_c2, PC_DELTAF2_c3, PC_DELTAF2_c4,PC_hA1,PC_delta)
      case (7)
        call set_cis2PE_vdel(PC_DELTAF3_c1, PC_DELTAF3_c2, PC_DELTAF3_c3, PC_DELTAF3_c4,PC_hA1,PC_delta)
      case (8)
        call set_cis2PE_vdel(PC_DELTAF4_c1, PC_DELTAF4_c2, PC_DELTAF4_c3, PC_DELTAF4_c4,PC_hA1,PC_delta)
      case (9)
        call set_cis2PE_vdel(PC_DELTAF5_c1, PC_DELTAF5_c2, PC_DELTAF5_c3, PC_DELTAF5_c4,PC_hA1,PC_delta)
      case (10)
        call set_cis2PE_vdel(PC_DELTAF6_c1, PC_DELTAF6_c2, PC_DELTAF6_c3, PC_DELTAF6_c4,PC_hA1,PC_delta)
      case (11)
        call set_cis2PE_vdel(PC_DELTAF7_c1, PC_DELTAF7_c2, PC_DELTAF7_c3, PC_DELTAF7_c4,PC_hA1,PC_delta)
      case (12)
        call set_cis2PE_vdel(PC_DELTAF8_c1, PC_DELTAF8_c2, PC_DELTAF8_c3, PC_DELTAF8_c4,PC_hA1,PC_delta)
      case (13)
        call set_cis2PE_vdel(PC_DELTAF9_c1, PC_DELTAF9_c2, PC_DELTAF9_c3, PC_DELTAF9_c4,PC_hA1,PC_delta)
      ! case default
      !   write (standard_error_unit, '(a)')  &
      !   & 'Select_cis2PE_vtpe1: unknown id of ci set. Default set pre-defined will be used'
        ! call set_cis2PE_vtpe1(PC_default_c1, PC_default_c3, PC_default_c4)
    end select

  end subroutine Select_cis2PE_vdel



   subroutine get_v_w(q,wt,ws,vc)
      real(NER), intent(in)   :: q

      real(NER)              ::sigma,LCq,Acq,Hcq,Dcq

      real(NER), intent(out)  ::wt,ws,vc
      real(NER)   ::W, hasqr, PC_fpi4, pisqr,  wsqr ,gasqr ,wv_fac ,PC_mpisqr, qsqr
      real(NER)        :: c1,c2,c3,c4,del,ha
      c1 = dps%c1
      c2 = dps%c2
      c3 = dps%c3
      c4 = dps%c4
      del= dps%del
      ha = dps%ha
!debug
!   print *, "vdel=", c1,c2,c3,c4
      qsqr=q*q
      hasqr=ha*ha
      PC_fpi4=PC_fpi*PC_fpi*PC_fpi*PC_fpi
      pisqr=PI_NE *PI_NE
      w=sqrt(q*q+4.0_NER*PC_mpi*PC_mpi)
      wsqr=w*w
      gasqr=PC_gA*PC_gA
      PC_mpisqr=PC_mpi*PC_mpi
      wv_fac=-gasqr*hasqr/(pisqr*PC_fpi4)

      call Dc_q(q,Dcq)
      call get_LAH(q,sigma,LCq,Acq,Hcq)

      vc=wv_fac/gasqr*del/18.0_NER*(6*sigma*(4.0_NER*c1*PC_mpisqr-2.0_NER*c2*del*del-c3*(2.0_NER &
      &  *del*del+sigma))*Dcq+(-24.0_NER*c1*PC_mpisqr+c2*(wsqr-6.0_NER*sigma)+6.0_NER*c3*(2*del* &
      &  del+sigma))*LCq)
      wt=wv_fac/gasqr*del/72.0_NER*c4*((wsqr-4.0_NER*del*del)*Dcq-2.0_NER*LCq)
      ws=-qsqr*wt

  end subroutine get_v_w

  subroutine get_U_from_v_w(j,q,Ut,Us,Uc)
    real(NER)            :: Uc,Us,Ut,wt,ws,vc
    integer , intent(in) :: j
    real(NER),intent(in) :: q

    call get_v_w(q,wt,ws,vc)
    if (mod(j,2).eq.1) then
       Uc=vc
       Us=-3.0_NER*ws
       Ut=-3.0_NER*wt
    else
       Uc=vc
       Us=ws
       Ut=wt
    end if
  end subroutine get_U_from_v_w


  function TPEdel_sub_j0j(j,p,k)

    real(NER)               :: TPEdel_sub_j0j,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, q
    if (j .lt. 0) then
      TPEdel_sub_j0j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if
    regsum = 0.0_NER
    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
        q      = sqrt(p*p+k*K-2.0_NER*p*k*z)
      call get_U_from_v_w(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)*(Uc-3.0_NER &
      &  *Us-(p*p+k*k-2.0_NER*z*p*k)*Ut)
    end do

    TPEdel_sub_j0j=regsum

  end function TPEdel_sub_j0j


  function TPEdel_sub_j1j(j,p,k)
    real(NER)               :: TPEdel_sub_j1j,Ut,Us,Uc,wt,ws,vc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, q

    if (j .lt. 1) then
      TPEdel_sub_j1j = 0.0_NER
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
          call get_v_w(q,wt,ws,vc)
            if (mod(j,2).eq.0) then
              Uc=vc
              Us=-3.0_NER*ws
              Ut=-3.0_NER*wt
           else
              Uc=vc
              Us=ws
              Ut=wt
           end if
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*((lgndr_table_pwo(j+1, ii)+lgndr_table_pwo(j-1, ii)) &
      &  *(-2.0_NER)*p*k*Ut+(Uc+Us+Ut*(psqr+k*k+2.0_NER*p*k*z))*lgndr_table_pwo(j, ii))*wghtslg_pwo(ii)
    end do

    TPEdel_sub_j1j=regsum

  end function TPEdel_sub_j1j


  function TPEdel_sub_jpp(j,p,k)
    real(NER)               :: TPEdel_sub_jpp,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 0) then
      TPEdel_sub_jpp = 0.0_NER
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
      call get_U_from_v_w(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*2.0_NER/(2.0_NER*j+1.0_NER) &
      &  *Ut*lgndr_table_pwo(j, ii)+(Uc+Us+1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j+1, ii))
    end do

    TPEdel_sub_jpp=regsum

  end function TPEdel_sub_jpp


  function TPEdel_sub_jmm(j,p,k)
    real(NER)               :: TPEdel_sub_jmm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 1) then
      TPEdel_sub_jmm = 0.0_NER
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
      call get_U_from_v_w(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*(-2.0_NER)/(2.0_NER*j+1.0_NER) &
      &  *Ut*lgndr_table_pwo(j, ii)+(Uc+Us-1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j-1, ii))
    end do

    TPEdel_sub_jmm=regsum

  end function TPEdel_sub_jmm


  function TPEdel_sub_jpm(j,p,k)
    real(NER)               :: TPEdel_sub_jpm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr , q

    if (j .lt. 1) then
      TPEdel_sub_jpm = 0.0_NER
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
      call get_U_from_v_w(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*sqrt(j*(j+1.0_NER))/(2.0_NER*j+ &
      &  1.0_NER)*((-4.0_NER)*Ut*p*k*lgndr_table_pwo(j, ii)+ksqr*2.0_NER*Ut*lgndr_table_pwo(j-1, ii)&
      &  +2.0_NER *psqr*Ut*lgndr_table_pwo(j+1, ii))
      TPEdel_sub_jpm=regsum
    end do

   end function TPEdel_sub_jpm


end module mod_deltafullN2LO_pwd
