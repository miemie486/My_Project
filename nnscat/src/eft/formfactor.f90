!wenchao shi    2019/10/26

module formfactor
  use util_cbcspn
  use wave_func_deute
  use util_cheb
  implicit none

  real(NER),parameter :: xdn_2=-1.0
  real(NER),parameter :: xup_2=1.0
  integer,parameter   :: N_2=200

  real(NER) :: q,lam

  real(NER) :: Mambda_common=800.0
  real(NER) :: C_3s1_LO_common(1:4)=(/2.9148117222315998E-003 ,  0.0,  0.0,  0.0/)
  real(NER) :: C_3s1_N2LO_common(1:4)=(/0.30631172064254268E-002  ,  0.13610613710601530E-001  , &
  &-0.37147265609656745E-007  , -0.22990080481946310E-007/)
  real(NER) :: lambda_common

  contains
 
  function f(p,q,z)

    real(NER) :: p,q,z,f
    f=sqrt(p*p+p*q*z+0.25*q*q)
    
  end function f


  function g(p,q,z)

    real(NER) :: p,q,z,g
    g=(p*p+0.5*p*q*z)/(p*sqrt(p*p+p*q*z+0.25*q*q))

  end function g


  function G_C(q)
  
    real(NER)    :: q,G_C
    real(NER)    :: sum
    integer      :: jj,kk,ll
    real(NER)    :: p,z
    real(NER)    :: x(N),weig(N),x2(N_2),weig2(N_2)

    real(NER)    :: x_input,yu,yw
    real(NER)    :: y2a_1(N),y2a_2(N),yau(N),yaw(N)
    complex(NEC) :: yau_com(N),yaw_com(N)
    REAL(NER) :: yp1,ypn
    real(NER)    :: norm_factor

    open(unit=11,file="deuteron/eigenvector.out")
    do ll=1,N
      read(11,*) yau_com(ll)

      read(11,*) yaw_com(ll)
    end do
    close(11)
   do ll=1,N
    yau(ll)=real(yau_com(ll))
    yaw(ll)=real(yaw_com(ll))
   end do

   call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)
   call feb_glabsc_2pt(N_2, x2, weig2, xdn_2, xup_2)

   norm_factor=0.0_NER
   do ll=1,N
      norm_factor=norm_factor+weig(ll)*x(ll)*x(ll)*( yau(ll)*yau(ll)+yaw(ll)*yaw(ll) )
   end do

   do ll=1,N
      yau(ll)=x(ll)*yau(ll)/sqrt(norm_factor)
      yaw(ll)=x(ll)*yaw(ll)/sqrt(norm_factor)
   end do

   yp1=0.0_NER
   ypn=0.0_NER
   sum=0.0_NER
   do jj=1,N
     do kk=1,N_2
        p=x(jj)
        z=x2(kk)
        x_input=f(p,q,z)

        call spline(x,yau,N,yp1,ypn,y2a_1)
        call spline(x,yaw,N,yp1,ypn,y2a_2)
        call splint(x,yau,y2a_1,N,x_input,yu)
        call splint(x,yaw,y2a_2,N,x_input,yw)
        sum=sum+p*p*(  3.0_NER*yu*yau(jj)/(f(p,q,z)*p)+((9.0_NER/2.0_NER)*g(p,q,z)*g(p,q,z)-(3.0_NER/2.0_NER))&
        &*yw*yaw(jj)/(f(p,q,z)*p)  )*weig(jj)*weig2(kk)
      end do
   end do

     G_C=(1.0_NER/6.0_NER)*sum*G_ES(q)


  end function G_C


 function G_ES(q)

   integer :: ll,dd
   integer,parameter :: Np=25,Nn=22
   real(NER) :: q_Ep(Np),G_Ep(Np),q_En(Nn),G_En(Nn)
   real(NER)  :: yp1,ypn,y2a_1(Np),y2a_2(Nn)
   real(NER)   :: q,G_ES
   real(NER)  :: G_Eplx,G_Enlx

   open(unit=12,file="deuteron/G_Ep.in")
   open(unit=13,file="deuteron/G_En.in")
   do ll=1,Np
     read(12,*) q_Ep(ll),G_Ep(ll)
   end do
   do dd=1,Nn
     read(13,*) q_En(dd),G_En(dd)
   end do
   close(12)
   close(13)

   yp1=0.0_NER
   ypn=0.0_NER

   call spline(q_Ep,G_Ep,Np,yp1,ypn,y2a_1)
   call spline(q_En,G_En,Nn,yp1,ypn,y2a_2)
   call splint(q_Ep,G_Ep,y2a_1,Np,q,G_Eplx)
   call splint(q_En,G_En,y2a_2,Nn,q,G_Enlx)
   
   G_ES=G_Eplx/( (1.0_NER+1.408450704E-6_NER*q*q)*(1.0_NER+1.408450704E-6_NER*q*q) )+G_Enlx

  end function G_ES


  subroutine getG_C_lamb(lamb,G_C_lamb)
!     integer,parameter :: q_num=10
!     real(NER)    :: G_C_lamb(q_num),lamb,q,G_C
    real(NER)    :: G_C_lamb,lamb,G_C
    real(NER)    :: sum
    integer      :: jj,kk,ll,bb,ff,hh
    real(NER)    :: p,z
    real(NER)    :: x(N),weig(N),x2(N_2),weig2(N_2)

    real(NER)    :: x_input,yu,yw
    real(NER)    :: y2a_1(N),y2a_2(N),yau(N),yaw(N)
    complex(NEC) :: yau_com(N),yaw_com(N)
    REAL(NER) :: yp1,ypn
    real(NER)    :: norm_factor
    real(NER)    :: G_C_lam_zanc(21),lam_input(21)
    integer,parameter :: lam_num=21


   open(unit=19,file="deuteron/lam.in")
   do hh=1,21 
   read(19,*) lam_input(hh)
   end do
   close(19)
   call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)
   call feb_glabsc_2pt(N_2, x2, weig2, xdn_2, xup_2)
   yp1=0.0_NER
   ypn=0.0_NER


!   open(unit=15,file="deuteron/qlist.in") 
!   do bb=1,q_num
!    read(15,*) q

   open(unit=16,file="deuteron/lam_eigenvector.out")
   do ff=1,21
   
    do ll=1,N
      read(16,*) yau_com(ll)

      read(16,*) yaw_com(ll)
    end do

   do ll=1,N
    yau(ll)=real(yau_com(ll))
    yaw(ll)=real(yaw_com(ll))
   end do

   norm_factor=0.0_NER
   do ll=1,N
      norm_factor=norm_factor+weig(ll)*x(ll)*x(ll)*( yau(ll)*yau(ll)+yaw(ll)*yaw(ll) )
   end do

   do ll=1,N
      yau(ll)=x(ll)*yau(ll)/sqrt(norm_factor)
      yaw(ll)=x(ll)*yaw(ll)/sqrt(norm_factor)
   end do

  
   sum=0.0_NER
   do jj=1,N
     do kk=1,N_2
        p=x(jj)
        z=x2(kk)
        x_input=f(p,q,z)

        call spline(x,yau,N,yp1,ypn,y2a_1)
        call spline(x,yaw,N,yp1,ypn,y2a_2)
        call splint(x,yau,y2a_1,N,x_input,yu)
        call splint(x,yaw,y2a_2,N,x_input,yw)
        sum=sum+p*p*(  3.0_NER*yu*yau(jj)/(f(p,q,z)*p)+((9.0_NER/2.0_NER)*g(p,q,z)*g(p,q,z)-(3.0_NER/2.0_NER))&
        &*yw*yaw(jj)/(f(p,q,z)*p)  )*weig(jj)*weig2(kk)
      end do
   end do
     G_C_lam_zanc(ff)=(1.0_NER/6.0_NER)*sum*G_ES(q)

    end do
   
      call spline(lam_input, G_C_lam_zanc,lam_num,yp1,ypn,y2a_1)
      call splint(lam_input,G_C_lam_zanc,y2a_1,lam_num,lamb,G_C)
        
!      G_C_lamb(bb)=G_C
     G_C_lamb=G_C
  close(16)

  end subroutine getG_C_lamb




  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! method three
  subroutine get_A_matrix_M_C_l_common(E,A_matrix)
     
    real(NER), intent(in) :: E
!     real(NER)             :: eta_no_lam
    complex(NEC)          :: A_matrix(1:2*N,1:2*N)
!     complex(NEC)          :: eta, Gamma(1:2*N)
!     integer               :: flag



    integer                                  :: L, S, j, regtype
!     real(NER)                                :: CsLO(1:4),CsN2LO(1:4)
    real(NER),dimension(1:2,1:2)             :: potvalLO,potvalN2LO,potval
!     real(NER),intent(in)                     :: p1, p2get_A_matrix_lamE(E,lambd,A_matrix)
     

!     real(NER),dimension(1:2,1:2),intent(out) :: V_3s1_3d1
!     real(NER)                                :: MambdaN2LO



!     real(NER)    :: E
!     complex(NEC) :: A_matrix(1:2*N,1:2*N)
    integer   :: ss,tt
    real(NER) :: A_mm(1:N,1:N),A_mp(1:N,1:N),A_pm(1:N,1:N),A_pp(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N) 



!   MambdaN2LO=800
!   data CsLO /2.9148117222315998E-003 ,  0.0,  0.0,  0.0/
!   data CsN2LO /0.30631172064254268E-002  ,  0.13610613710601530E-001  , &
!   &-0.37147265609656745E-007  , -0.22990080481946310E-007/
  
! data   Mambda,Cs(1)    /    450 , -3.6804098078605466E-003     /  
   
    L=0
    S=1
    j=1
    regtype=REGTYPE_GAUSSIAN

    call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)

   do ss=1,N
      do tt=1,N
      p=x(ss)
      k=x(tt)

 call VLO_withc0_epwrap(L, S, j, regtype, Mambda_common, C_3s1_LO_common, p,k, potvalLO)
 call VN2LO_withc0_epwrap(L, S, j, regtype, Mambda_common, C_3s1_N2LO_common, p, k , potvalN2LO)
    potval(1,1)=potvalLO(1,1)+potvalN2LO(1,1)*lambda_common
    potval(1,2)=potvalLO(1,2)+potvalN2LO(1,2)*lambda_common
    potval(2,1)=potvalLO(2,1)+potvalN2LO(2,1)*lambda_common
    potval(2,2)=potvalLO(2,2)+potvalN2LO(2,2)*lambda_common


      A_mm(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*potval(1,1)*weig(tt)
      A_mp(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*potval(1,2)*weig(tt)
      A_pm(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*potval(2,1)*weig(tt)
      A_pp(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*potval(2,2)*weig(tt)
      end do 
    end do 

   do ss=1,N
      do tt=1,N
      A_matrix(2*ss-1,2*tt-1)=dcmplx(A_mm(ss,tt),0.0_NER)
      A_matrix(2*ss-1,2*tt)=dcmplx(A_mp(ss,tt),0.0_NER)
      end do 
   end do

   do ss=1,N
     do tt=1,N
     A_matrix(2*ss,2*tt-1)=cmplx(A_pm(ss,tt),0.0_NER)
     A_matrix(2*ss,2*tt)=cmplx(A_pp(ss,tt),0.0_NER)
     end do 
   end do
 end subroutine get_A_matrix_M_C_l_common


  subroutine get_Bd_3s1(x1,x2,bindenergy)

  REAL(NER) :: x1,x2
  REAL(NER) :: bindenergy
 
  complex(NEC) :: A_matrix(1:2*N,1:2*N)
  complex(NEC)  :: eta, Gamma(1:2*N)
  integer       :: flag
  integer :: dd

  call get_root(eta_M_C_l_common,x1,x2,bindenergy)

 end subroutine get_Bd_3s1

function eta_M_C_l_common(E)
          
    real(NER), intent(in) :: E
    real(NER)             :: eta_M_C_l_common
    complex(NEC)          :: A_matrix(1:2*N,1:2*N)
    complex(NEC)          :: eta, Gamma(1:2*N)
    integer               :: flag
    
    call get_A_matrix_M_C_l_common(E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag)   
    eta_M_C_l_common=real(eta)-1.0_NER

  end function eta_M_C_l_common


  subroutine get_G_C_common(lambda, q, G_C_M_C_l_common)
  
    real(NER),intent(in)    :: q,lambda
    real(NER)    :: G_C_M_C_l_common
    real(NER)    :: sum
    integer      :: jj,kk,ll
    real(NER)    :: p,z
    real(NER)    :: x(N),weig(N),x2(N_2),weig2(N_2)

    real(NER)    :: x_input,yu,yw
    real(NER)    :: y2a_1(N),y2a_2(N),yau(N),yaw(N)
    complex(NEC) :: yau_com(N),yaw_com(N)
    REAL(NER) :: yp1,ypn
    real(NER)    :: norm_factor
    
    real(NER)             :: eta_M_C_l_common
    complex(NEC)          :: A_matrix(1:2*N,1:2*N)
    complex(NEC)          :: eta, Gamma(1:2*N)
    integer               :: flag
    REAL(NER) :: x1,y1
    REAL(NER) :: bindenergy
    integer :: dd

  lambda_common=lambda
  x1=-3.0_NER
  y1=-1.0_NER
  call get_Bd_3s1(x1,y1,bindenergy)
  call get_A_matrix_M_C_l_common(bindenergy,A_matrix) 
  call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag) 
!   open(unit=99,file="deuteron/eigenvector_Mam_C_lam_common.out") 
!   do dd=1,2*N
!   write(99,*) Gamma(dd)
!   end do
!   close(99)



!     open(unit=11,file="deuteron/eigenvector_Mam_C_lam_common.out") 
    do ll=1,N
       yau_com(ll)=Gamma(2*ll-1)
       yaw_com(ll)=Gamma(2*ll)
    end do
!     close(11)
   do ll=1,N
    yau(ll)=real(yau_com(ll))
    yaw(ll)=real(yaw_com(ll))
   end do

   call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)
   call feb_glabsc_2pt(N_2, x2, weig2, xdn_2, xup_2)

   norm_factor=0.0_NER
   do ll=1,N
      norm_factor=norm_factor+weig(ll)*x(ll)*x(ll)*( yau(ll)*yau(ll)+yaw(ll)*yaw(ll) )
   end do

   do ll=1,N
      yau(ll)=x(ll)*yau(ll)/sqrt(norm_factor)
      yaw(ll)=x(ll)*yaw(ll)/sqrt(norm_factor)
   end do

   yp1=0.0_NER
   ypn=0.0_NER
   sum=0.0_NER
   do jj=1,N
     do kk=1,N_2
        p=x(jj)
        z=x2(kk)
        x_input=f(p,q,z)

        call spline(x,yau,N,yp1,ypn,y2a_1)
        call spline(x,yaw,N,yp1,ypn,y2a_2)
        call splint(x,yau,y2a_1,N,x_input,yu)
        call splint(x,yaw,y2a_2,N,x_input,yw)
        sum=sum+p*p*(  3.0_NER*yu*yau(jj)/(f(p,q,z)*p)+((9.0_NER/2.0_NER)*g(p,q,z)*g(p,q,z)-(3.0_NER/2.0_NER))&
        &*yw*yaw(jj)/(f(p,q,z)*p)  )*weig(jj)*weig2(kk)
      end do
   end do

     G_C_M_C_l_common=(1.0_NER/6.0_NER)*sum*G_ES(q)


  end subroutine get_G_C_common

  subroutine G_C_taylor_lambda(q_in)
    real(NER) :: q_in
    real(NER)                  :: a, b
    integer                    :: nl=4
    real(NER)    :: gl(4)
    
     open(unit=73,file="deuteron/G_C_lambda_taylor.out")
      b=0.001
      a=-b
      call cheb_taylor_real(a, b, nl, gl, G_C_lambda)
      write(73,*) q_in,gl

  contains
    function G_C_lambda(lambda)

      real(NER),intent(in) :: lambda
      real(NER) :: G_C_M_C_l_common, G_C_lambda
      call get_G_C_common(lambda, q_in, G_C_M_C_l_common)
      G_C_lambda=G_C_M_C_l_common

    end function G_C_lambda
 end subroutine G_C_taylor_lambda
! method three
 ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module formfactor 












