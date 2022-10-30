! Rui Peng circa 08/2018

! Leading Deltaful TPE
module mod_deltafullNLO_pwd
  ! Put behind "only:" the name of variables, functions and subs your code is
  ! going to refer to in this module
  use nneft_type,       only: NER, NEC, PI_NE
  ! Check nneft_phyconst for name of common physics parameters.
  use nneft_phyconst,   only: PC_mN, PC_gA, PC_fpi, PC_mpi,  Delta_paras, PC_delta_cutoff
  use ope_pwd
  use util_gauleg,   only:feb_glabsc_2pt





  implicit none

  logical, private              :: inited_Dcq = .false.
  integer, parameter, private   :: dcq_N = 100
  real(NER), private            :: u_i(1:dcq_N), weig_i(1:dcq_N), Ldelq

  type (Delta_paras)            :: dps

  ! logical             :: mshlg_inited_pwo = .false.
  ! integer, parameter  :: Nlg_pwo = 100
  ! real(NER) :: mshlg_pwo(1:Nlg_pwo), wghtslg_pwo(1:Nlg_pwo)
  ! integer, parameter  :: lgrank_pwo = 20
  ! real(NER) :: lgndr_table_pwo(0:lgrank_pwo, 1:Nlg_pwo)



  ! Define your own physics parameters here if you do not find them in
  ! nneft_phyconst

contains

  ! Add below your function that generate partial-wave matrix elements
  ! For example:

  ! subroutine set_delta0_paras(input_delta,cutoff)

  !   real(NER), intent(in) :: input_delta

  !   del = input_delta

  ! end subroutine set_delta0_paras

  subroutine set_2PEparas_vdel(ha, del)
    real(NER), intent(in)  :: ha, del

    dps%ha = ha
    dps%del = del

  end subroutine set_2PEparas_vdel

  subroutine init_Dcq()

    real(NER)                           :: M, cut, q, w, s, del, cutsqr, qsqr

    del = dps%del
    M = 2.0_NER * PC_mpi
    cut = PC_delta_cutoff
    cutsqr = cut * cut
    if (inited_Dcq) then
      return
    end if
    call feb_glabsc_2pt(dcq_N,u_i,weig_i,M,cut)
    q = 2.0_NER * sqrt(del*del - PC_mpi*PC_mpi)
    qsqr = q*q
    w = sqrt(q*q + 4.0_NER*PC_mpi*PC_mpi)
    s = sqrt(cut*cut - 4.0_NER*PC_mpi*PC_mpi)
    Ldelq = w/(2.0_NER*q)*log((cutsqr*w*w+qsqr*s*s+2.0_NER*cut*q*w*s)/(4.0_NER*PC_mpi*PC_mpi*(cutsqr+qsqr)))

    inited_Dcq = .true.

  end subroutine init_Dcq

  subroutine get_LAH(q,sigma,LCq,Acq,Hcq)

      real(NER), intent(in)   :: q



      real(NER), intent(out)::sigma,LCq,Acq,Hcq
      real(NER) ::W, S, qsqr, PC_mpisqr, delsqr, cutsqr,seita
      real(NER) ::del, ha,cut

      del = dps%del
      ha = dps%ha
      cut = PC_delta_cutoff
      qsqr=q*q
      PC_mpisqr=PC_mpi*PC_mpi
      delsqr=del*del
      cutsqr=cut*cut
      seita=cut-2.0_NER*PC_mpi
      w=sqrt(qsqr+4.0_NER*PC_mpisqr)
      s=sqrt(cutsqr-4.0_NER*PC_mpisqr)

      sigma=2.0_NER*PC_mpisqr+qsqr-2.0_NER*delsqr

      if (seita  .gt. 0.0_NER)  then
      LCq= w/(2.0_NER*q)*log((cutsqr*w*w+qsqr*s*s+2.0_NER*cut*q*w*s)/(4.0_NER*PC_mpisqr*(cutsqr+qsqr)))
      Acq= 1.0_NER/(2.0_NER*q)*atan(q*(cut-2.0_NER*PC_mpi)/(qsqr+2.0_NER*cut*PC_mpi))
      Hcq= 2.0_NER*sigma/(w*w-4.0_NER*delsqr)*(LCq-Ldelq)
      else
      Lcq=0
      Acq=0
      Hcq=0
      end if
  end subroutine get_LAH


   subroutine Dc_q(q,Dcq)
        real(NER), intent(in)    :: q

        real(NER), intent(out)   :: Dcq
        real(NER)                :: W,S,u,qsqr,regsum,PC_mpisqr,del
        integer     :: ii

      if (.not. inited_Dcq) then
        call init_Dcq()
      end if

        del = dps%del
        regsum=0.0_NER
        qsqr=q*q
        PC_mpisqr=PC_mpi*PC_mpi

        do ii = 1, dcq_N
          u= u_i(ii)
          regsum=regsum+1.0_NER/(del*(u*u+qsqr))*atan(sqrt(u*u-4.0_NER*PC_mpisqr)/2.0_NER/del)*weig_i(ii)
        end do

        Dcq=regsum

    end subroutine Dc_q

   subroutine get_w_v(q,wc,wt,ws,vc,vt,vs)

      real(NER), intent(in)   :: q

      real(NER)              ::w_c1,w_c2,w_c3,v_c2,v_c3,v_t2,sigma,LCq,Acq,Hcq,Dcq
      real(NER)              ::v_t3,w_t2,w_t3,w_s2 ,w_s3,v_s2,v_s3
      real(NER), intent(out)  :: wc,wt,ws,vc,vt,vs
      real(NER)   ::W, hasqr, PC_fpi4, pisqr,  wsqr ,gasqr ,wv_fac ,PC_mpisqr, qsqr
      real(NER)   :: ha, del

      ha = dps%ha
      del = dps%del

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

      w_c1=-hasqr/(216.0_NER*pisqr*PC_fpi4)*((6.0_NER*sigma - wsqr )*LCq+12.0_NER*del*del*sigma*Dcq)
      w_c2=wv_fac/216.0_NER*((12.0_NER*del*del-20.0_NER*PC_mpisqr-11.0_NER*qsqr)*LCq+6.0_NER*sigma *sigma *Dcq)
      w_c3=wv_fac*hasqr/486.0_NER/gasqr*((12.0_NER*sigma -wsqr)*Lcq +3.0_NER*sigma*(Hcq+(8.0_NER*del*del-sigma)*Dcq))
      w_t2=wv_fac/144.0_NER/del*PI_NE *wsqr*Acq
      w_t3=wv_fac*hasqr/1296.0_NER/gasqr*(2.0_NER*Lcq+(4.0_NER*del*del+wsqr)*Dcq)
      w_s2=-qsqr*w_t2
      w_s3=-qsqr*w_t3
      v_c2=wv_fac*PI_NE /12.0_NER/del*(2.0_NER*PC_mpisqr+qsqr)*(2.0_NER*PC_mpisqr+qsqr)*Acq
      v_c3=wv_fac*hasqr/27.0_NER/gasqr*(-4.0_NER*del*del*Lcq+sigma*(Hcq+(sigma+8.0_NER*del*del)*Dcq))
      v_t2=wv_fac/48.0_NER*(-2.0_NER*Lcq+(wsqr-4.0_NER*del*del)*Dcq)
      v_t3=wv_fac*hasqr/gasqr/216.0_NER*(6.0_NER*Lcq+(12.0_NER*del*del-wsqr)*Dcq)
      v_s2=-qsqr*v_t2
      v_s3=-qsqr*v_t3

      wc=W_c1+w_c2+w_c3
      wt=w_t2+w_t3
      ws=w_s2+w_s3
      vc=v_c2+v_c3
      vt=v_t2+v_t3
      vs=v_s2+v_s3


    end subroutine get_w_v

  subroutine get_U_from_w_v(j,q,Ut,Us,Uc)
    real(NER)            :: Uc,Us,Ut,wc,wt,ws,vc,vt,vs
    integer , intent(in) :: j
    real(NER),intent(in) :: q



    call get_w_v(q,wc,wt,ws,vc,vt,vs)
    if (mod(j, 2) .eq. 1) then
       Uc=vc-3.0_NER*wc
       Us=vs-3.0_NER*ws
       Ut=vt-3.0_NER*wt
    else
       Uc=vc+wc
       Us=vs+ws
       Ut=vt+wt
    end if
  end subroutine get_U_from_w_v


  function TPEdel_j0j(j,p,k)

    real(NER)               :: TPEdel_j0j,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, q
    if (j .lt. 0) then
      TPEdel_j0j = 0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if
    regsum = 0.0_NER
    do ii = 1, Nlg_pwo
      z      = mshlg_pwo(ii)
        q      = sqrt(p*p+k*K-2.0_NER*p*k*z)
      call get_U_from_w_v(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)*(Uc-3.0_NER*Us-(p*p+k*k- &
      &      2.0_NER*z*p*k)*Ut)
    end do

    TPEdel_j0j=regsum

  end function TPEdel_j0j




   function TPEdel_j1j(j,p,k)
    real(NER)               :: TPEdel_j1j,Ut,Us,Uc,wc,wt,ws,vc,vt,vs
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, q

    if (j .lt. 1) then
      TPEdel_j1j = 0.0_NER
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
         call get_w_v(q,wc,wt,ws,vc,vt,vs)
           if (mod(j, 2) .eq. 0) then
             Uc=vc-3.0_NER*wc
             Us=vs-3.0_NER*ws
             Ut=vt-3.0_NER*wt
           else
             Uc=vc+wc
             Us=vs+ws
             Ut=vt+wt
           end if
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*((lgndr_table_pwo(j+1, ii)+lgndr_table_pwo(j-1, ii))*(-2.0_NER)*p*k*Ut+ &
      &      (Uc+Us+Ut*(psqr+k*k+2.0_NER*p*k*z))*lgndr_table_pwo(j, ii))*wghtslg_pwo(ii)
    end do

    TPEdel_j1j=regsum

  end function TPEdel_j1j




  function TPEdel_jpp(j,p,k)
    real(NER)               :: TPEdel_jpp,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 0) then
      TPEdel_jpp = 0.0_NER
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
      call get_U_from_w_v(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*2.0_NER/(2.0_NER*j+1.0_NER)*Ut*lgndr_table_pwo(j, ii) &
      &      +(Uc+Us+1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j+1, ii))
    end do

    TPEdel_jpp=regsum

  end function TPEdel_jpp


  function TPEdel_jmm(j,p,k)
    real(NER)               :: TPEdel_jmm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr, q

    if (j .lt. 1) then
      TPEdel_jmm = 0.0_NER
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
      call get_U_from_w_v(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*(p*k*(-2.0_NER)/(2.0_NER*j+1.0_NER)*Ut*lgndr_table_pwo(j, ii) &
      &      +(Uc+Us-1.0_NER/(2.0_NER*j+1.0_NER)*(psqr+ksqr)*(-Ut))*lgndr_table_pwo(j-1, ii))
    end do

    TPEdel_jmm=regsum

  end function TPEdel_jmm


  function TPEdel_jpm(j,p,k)
    real(NER)               :: TPEdel_jpm,Ut,Us,Uc
    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: j

    integer     :: ii
    real(NER)   :: regsum,  z, psqr, ksqr , q

    if (j .lt. 1) then
      TPEdel_jpm = 0.0_NER
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
      call get_U_from_w_v(j,q,Ut,Us,Uc)
      regsum=regsum+PC_mN/(4.0_NER*PI_NE*PI_NE)*wghtslg_pwo(ii)*sqrt(j*(j+1.0_NER))/(2.0_NER*j+1.0_NER)*((-4.0_NER)*Ut*p*k* &
      &      lgndr_table_pwo(j, ii)+ksqr*2.0_NER*Ut*lgndr_table_pwo(j-1, ii)+2.0_NER *psqr*Ut*lgndr_table_pwo(j+1, ii))
      TPEdel_jpm=regsum
    end do

   end function TPEdel_jpm



end module mod_deltafullNLO_pwd
