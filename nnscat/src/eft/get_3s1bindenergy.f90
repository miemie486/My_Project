!wenchao shi    2019/12/6

module  get_3s1bindenergy

  use wave_func_deute
  ! use wave_func_deute, only : get_root, potwrap_potval
  ! use nneft_phyconst
  ! use util_gauleg

  implicit none

contains
subroutine get_Bd_3s1(Mambda, para,pot,x1,x2,bindenergy)
  ! procedure(potwrap_potval) :: pot
  procedure(gnrc_eft_pot) :: pot

  REAL(NER) :: x1,x2
  REAL(NER) :: bindenergy
  real(NER)   :: Mambda, para(:)

  complex(NEC) :: A_matrix(1:2*N,1:2*N)
  complex(NEC)  :: eta, Gamma(1:2*N)
  integer       :: flag
  integer :: dd

!   x1=-3.0_NER
!   x2=-1.0_NER
  call get_root(eta_no,x1,x2,bindenergy)

contains
  function eta_no(E)

    real(NER), intent(in) :: E
    real(NER)             :: eta_no
    complex(NEC)          :: A_matrix(1:2*N,1:2*N)
    complex(NEC)          :: eta, Gamma(1:2*N)
    integer               :: flag
    integer                                  :: L, S, j, regtype
    real(NER),dimension(1:2,1:2)             :: potval
    integer   :: ss,tt
    real(NER) :: A_mm(1:N,1:N),A_mp(1:N,1:N),A_pm(1:N,1:N),A_pp(1:N,1:N)
    real(NER) :: p,k
    real(NER) :: x(N),weig(N)

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

call pot(L, S, j, regtype, Mambda, para, p, k, potval)

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
     A_matrix(2*ss,2*tt-1)=dcmplx(A_pm(ss,tt),0.0_NER)
     A_matrix(2*ss,2*tt)=dcmplx(A_pp(ss,tt),0.0_NER)
     end do
   end do
    call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag)
    eta_no=real(eta)-1.0_NER

  end function eta_no

 end subroutine get_Bd_3s1

end module  get_3s1bindenergy

