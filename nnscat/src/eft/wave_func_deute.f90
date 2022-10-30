!wenchao shi    2019/10/24

module wave_func_deute

	use util_gauleg
	use nneft_phyconst
  use rel_ope_pwd
  use vtpe0_pwd
  use vtpe1_pwd
  use potwrap
  use nneft_eigens
  use chengdu
  use util_rootfinder
	implicit none

  real(NER),parameter    :: xdn=0.0
  real(NER),parameter    :: xup=10000.0
  integer,parameter      :: N=200
!   real(NER),parameter    :: tan_M=350
 abstract interface
    subroutine potwrap_potval(L, S, j, regtype, Mambda, para, p1, p2, potval)
      import                  :: NER
      implicit none
      integer, intent(in)     :: L, S, j, regtype
      real(NER), intent(in)   :: Mambda, p1, p2, para(:)
      real(NER), intent(out)  :: potval(1:2, 1:2)
    end subroutine potwrap_potval
  end interface

 contains

  subroutine get_A_matrix(E,A_matrix)

    real(NER)    :: E
    complex(NEC) :: A_matrix(1:2*N,1:2*N)
    integer   :: ss,tt
    real(NER) :: A_mm(1:N,1:N),A_mp(1:N,1:N),A_pm(1:N,1:N),A_pp(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N)
    real(NER),dimension(1:2,1:2) :: V_3s1_3d1

    call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)

   do ss=1,N
      do tt=1,N
      p=x(ss)
      k=x(tt)
      call get_V_array(p,k,V_3s1_3d1)
      A_mm(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(1,1)*weig(tt)
      A_mp(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(1,2)*weig(tt)
      A_pm(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(2,1)*weig(tt)
      A_pp(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(2,2)*weig(tt)
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
     A_matrix(2*ss,2*tt-1)=dcmplx(A_pm(ss,tt),0.0_NER)
     A_matrix(2*ss,2*tt)=dcmplx(A_pp(ss,tt),0.0_NER)
     end do
   end do

  end subroutine get_A_matrix


  subroutine get_V_array(p1,p2,V_3s1_3d1)

   real(NER), intent(in)   :: p1, p2
   real(NER),dimension(1:2,1:2),intent(out) :: V_3s1_3d1
   integer     :: L, S, J, uptoQn
   real(NER),dimension(1:2,1:2) :: potval

   uptoQn=0
    L=0
    S=1
    j=1

    call chengdu_DLSPR_800(L, S, J, uptoQn, p1, p2, potval)
    V_3s1_3d1=potval

  end subroutine get_V_array

  subroutine get_root(func,x1,x2,tt)

    REAL(NER),intent(in)  :: x1,x2
    REAL(NER),intent(out) :: tt
    REAL(NER),PARAMETER   :: tol=1.0E-10
    INTERFACE
            FUNCTION func(x)
            USE nneft_type
            IMPLICIT NONE
            REAL(NER), INTENT(IN) :: x
            REAL(NER) :: func
            END FUNCTION func
        END INTERFACE
    tt=zbrent(func,x1,x2,tol)

  end subroutine get_root


  function eta_(E)

    real(NER), intent(in) :: E
    real(NER)             :: eta_
    complex(NEC)          :: A_matrix(1:2*N,1:2*N)
    complex(NEC)          :: eta, Gamma(1:2*N)
    integer               :: flag

    call get_A_matrix(E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag)
    eta_=real(eta)-1.0_NER

  end function eta_


  subroutine get_A_matrix_nowrap(E,A_matrix)

    real(NER)    :: E
    complex(NEC) :: A_matrix(1:2*N,1:2*N)
    integer   :: ss,tt
    real(NER) :: A_mm(1:N,1:N),A_mp(1:N,1:N),A_pm(1:N,1:N),A_pp(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N)
    real(NER),dimension(1:2,1:2) :: V_3s1_3d1
    call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)

   do ss=1,N
      do tt=1,N
      p=x(ss)
      k=x(tt)
      call get_V_array_nowrap(p,k,V_3s1_3d1)
      A_mm(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(1,1)*weig(tt)
      A_mp(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(1,2)*weig(tt)
      A_pm(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(2,1)*weig(tt)
      A_pp(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_3s1_3d1(2,2)*weig(tt)
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
     A_matrix(2*ss,2*tt-1)=dcmplx(A_pm(ss,tt),0.0_NER)
     A_matrix(2*ss,2*tt)=dcmplx(A_pp(ss,tt),0.0_NER)
     end do
   end do

  end subroutine get_A_matrix_nowrap


  subroutine get_V_array_nowrap(p1,p2,V_3s1_3d1)

    integer                                  :: L, S, j, regtype
    real(NER)                                :: Cs(1:4)
    real(NER),dimension(1:2,1:2)             :: potval
    real(NER),intent(in)                     :: p1, p2
    real(NER),dimension(1:2,1:2),intent(out) :: V_3s1_3d1
    real(NER)                                :: Mambda


  !     Cs(1)=2.9151499918865895E-003_NER
    Cs(1)=-1.4410458018551903E-003_NER
    Cs(2)=0.0_NER
    Cs(3)=0.0_NER
    Cs(4)=0.0_NER

    L=0
    S=1
    j=1
    regtype=REGTYPE_GAUSSIAN
!     Mambda=800.0_NER
    Mambda=600.0_NER

    call VLO_withc0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)
    V_3s1_3d1=potval

  end subroutine get_V_array_nowrap

  function eta_nowrap(E)

    real(NER), intent(in) :: E
    real(NER)             :: eta_nowrap
    complex(NEC)          :: A_matrix(1:2*N,1:2*N)
    complex(NEC)          :: eta, Gamma(1:2*N)
    integer               :: flag

    call get_A_matrix_nowrap(E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag)
    eta_nowrap=real(eta)-1.0_NER

  end function eta_nowrap


  subroutine get_A_matrix_sngl(E,A_matrix)

    real(NER)    :: E
    complex(NEC) :: A_matrix(1:N,1:N)
    integer   :: ss,tt
    real(NER) :: A(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N)
    real(NER) :: V_sngl
    call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)

   do ss=1,N
      do tt=1,N
      p=x(ss)
      k=x(tt)
      call get_V_array_sngl(p,k,V_sngl)
      A(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_sngl*weig(tt)
      end do
    end do

   do ss=1,N
      do tt=1,N
      A_matrix(ss,tt)=dcmplx(A(ss,tt),0.0_NER)
      end do
   end do

  end subroutine get_A_matrix_sngl


  subroutine get_V_array_sngl(p1,p2,V_sngl)

   real(NER), intent(in)   :: p1, p2
   real(NER),intent(out) :: V_sngl
   integer     :: L, S, J, uptoQn
   real(NER),dimension(1:2,1:2) :: potval

   uptoQn=0
    L=0
    S=0
    j=0

    call chengdu_DLSPR_1600(L, S, J, uptoQn, p1, p2, potval)

    V_sngl=potval(1,1)

  end subroutine get_V_array_sngl

  function eta_sngl(E)

    real(NER), intent(in) :: E
    real(NER)             :: eta_sngl
    complex(NEC)          :: A_matrix(1:N,1:N)
    complex(NEC)          :: eta, Gamma(1:N)
    integer               :: flag

    call get_A_matrix_sngl(E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, N, eta, Gamma, flag)
    eta_sngl=real(eta)-1.0_NER

  end function eta_sngl

  subroutine get_A_matrix_sngl_nowrap(E,A_matrix)

    real(NER)    :: E
    complex(NEC) :: A_matrix(1:N,1:N)
    integer   :: ss,tt
    real(NER) :: A(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N)
    real(NER) :: V_sngl
    call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)

   do ss=1,N
      do tt=1,N
      p=x(ss)
      k=x(tt)
      call get_V_array_sngl_nowrap(E,p,k,V_sngl)
      A(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_sngl*weig(tt)
      end do
    end do

   do ss=1,N
      do tt=1,N
      A_matrix(ss,tt)=dcmplx(A(ss,tt),0.0_NER)
      end do
   end do

  end subroutine get_A_matrix_sngl_nowrap


  subroutine get_V_array_sngl_nowrap(E,p1,p2,V_sngl)

   real(NER), intent(in)   :: E,p1, p2
   real(NER),intent(out) :: V_sngl

    integer    :: L, S, j, regtype
    real(NER)  :: Mambda, para(1:2)
    real(NER)  :: potval(1:2, 1:2)

    real(NER) :: vope(1:2, 1:2)
    data   Mambda,para(1),para(2)   &
    &/      800.00 ,   0.59100655067130575E+002,   -0.72482410781848861E-001      /
    L=0
    S=0
    j=0
    regtype=REGTYPE_GAUSSIAN


    call VLO_LBWDIB_epwrap(L, S, j, regtype, Mambda, para, E, p1, p2, potval)

    V_sngl=potval(1,1)

  end subroutine get_V_array_sngl_nowrap

 function eta_sngl_nowrap(E)

    real(NER), intent(in) :: E
    real(NER)             :: eta_sngl_nowrap
    complex(NEC)          :: A_matrix(1:N,1:N)
    complex(NEC)          :: eta, Gamma(1:N)
    integer               :: flag

    call get_A_matrix_sngl_nowrap(E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, N, eta, Gamma, flag)
    eta_sngl_nowrap=real(eta)-1.0_NER

  end function eta_sngl_nowrap


  subroutine get_A_matrix_sngl_nowrap_Mambda_C(Mambda,E, para, A_matrix)

    real(NER), intent(in)    :: E
    complex(NEC) :: A_matrix(1:N,1:N)
    integer   :: ss,tt
    real(NER) :: A(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N)
    real(NER) :: V_sngl
    real(NER)  :: Mambda, para(:)
    call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N, tan_M*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)

   do ss=1,N
      do tt=1,N
      p=x(ss)
      k=x(tt)
      call get_V_array_sngl_nowrap_Mambda_C(Mambda,para, E, p,k,V_sngl)
      A(ss,tt)=(k*k)/(E-p*p/PC_mN)*(1.0_NER/PC_mN)*V_sngl*weig(tt)
      end do
    end do

   do ss=1,N
      do tt=1,N
      A_matrix(ss,tt)=dcmplx(A(ss,tt),0.0_NER)
      end do
   end do

  end subroutine get_A_matrix_sngl_nowrap_Mambda_C

  subroutine get_V_array_sngl_nowrap_Mambda_C(Mambda,para, E, p1,p2,V_sngl)

   real(NER), intent(in)   :: E, p1, p2
   real(NER),intent(out) :: V_sngl

    integer    :: L, S, j, regtype
    real(NER)  :: Mambda, para(:)
    real(NER)  :: potval(1:2, 1:2)

    real(NER) :: vope(1:2, 1:2)

    L=0
    S=0
    j=0
    regtype=REGTYPE_GAUSSIAN


    call VLO_LBWDIB_epwrap(L, S, j, regtype, Mambda, para, E, p1, p2, potval)

    V_sngl=potval(1,1)

  end subroutine get_V_array_sngl_nowrap_Mambda_C











  subroutine get_A_matrix_fit(Mambda,C_3s1,E,A_matrix)

    real(NER),intent(in)   :: Mambda,C_3s1,E
    complex(NEC) :: A_matrix(1:2*N,1:2*N)
    integer   :: ss,tt
    real(NER) :: A_mm(1:N,1:N),A_mp(1:N,1:N),A_pm(1:N,1:N),A_pp(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N)
    real(NER),dimension(1:2,1:2) :: V_3s1_3d1

    call feb_glabsc_2pt(N, x, weig, xdn, xup)
!     call feb_tanabsc(N,Mambda*abs(log(0.5_NER))**(1.0_NER/real(GAUSSIAN_POWER)), x, weig)

   do ss=1,N
      do tt=1,N
      p=x(ss)
      k=x(tt)
      call get_V_array_fit(Mambda,C_3s1,p,k,V_3s1_3d1)
      A_mm(ss,tt)=(k*k)/(E-p*p/PC_mN)*V_3s1_3d1(1,1)*weig(tt)
      A_mp(ss,tt)=(k*k)/(E-p*p/PC_mN)*V_3s1_3d1(1,2)*weig(tt)
      A_pm(ss,tt)=(k*k)/(E-p*p/PC_mN)*V_3s1_3d1(2,1)*weig(tt)
      A_pp(ss,tt)=(k*k)/(E-p*p/PC_mN)*V_3s1_3d1(2,2)*weig(tt)
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
     A_matrix(2*ss,2*tt-1)=dcmplx(A_pm(ss,tt),0.0_NER)
     A_matrix(2*ss,2*tt)=dcmplx(A_pp(ss,tt),0.0_NER)
     end do
   end do

  end subroutine get_A_matrix_fit

  subroutine get_V_array_fit(Mambda,C_3s1,p1,p2,V_3s1_3d1)

    integer                                  :: L, S, j, regtype
    real(NER)                                :: Cs(1:4)
    real(NER),dimension(1:2,1:2)             :: potval
    real(NER),intent(in)                     :: p1, p2, C_3s1
    real(NER),dimension(1:2,1:2),intent(out) :: V_3s1_3d1
    real(NER),intent(in)         :: Mambda

    Cs(1)=C_3s1
    Cs(2)=0.0_NER
    Cs(3)=0.0_NER
    Cs(4)=0.0_NER

    L=0
    S=1
    j=1
    regtype=REGTYPE_GAUSSIAN


    call VLO_withc0_epwrap(L, S, j, regtype, Mambda, Cs, p1, p2, potval)
    V_3s1_3d1=(1.0_NER/PC_mN)*potval

  end subroutine get_V_array_fit

  function eta_C_8(C_3s1)

    real(NER),parameter   :: E=-2.225
    real(NER),intent(in)  :: C_3s1
    real(NER)             :: eta_C_8
    complex(NEC)          :: A_matrix(1:2*N,1:2*N)
    complex(NEC)          :: eta, Gamma(1:2*N)
    integer               :: flag
    real(NER)             :: Mambda=500.0

    call get_A_matrix_fit(Mambda,C_3s1,E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag)
    eta_C_8=real(eta)-1.0_NER

  end function eta_C_8
end module wave_func_deute










