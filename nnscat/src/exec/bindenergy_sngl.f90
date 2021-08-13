!wenchao shi    2019/10/24

program bindenergy_sngl

  use wave_func_deute
  use nneft_eigens
  use util_rootfinder
  implicit none

!   REAL(NER) :: x1,x2
!   REAL(NER) :: bindenergy
!   complex(NEC) :: A_matrix(1:N,1:N)
!   complex(NEC)  :: eta, Gamma(1:N)
!   integer       :: flag
!   integer :: dd
! !   x1=-100.0_NER
! !   x2=-40.0_NER
!   x1=-40.0_NER
!   x2=-1.0_NER
!   call get_root(eta_sngl,x1,x2,bindenergy)
!   write(*,*) bindenergy


    real(NER) :: E
    complex(NEC)          :: A_matrix(1:N,1:N)
    complex(NEC)          :: eta, Gamma(1:N)
    integer               :: flag,ll
 
    do ll=1,30
    read(*,*) E
    call get_A_matrix_sngl(E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, N, eta, Gamma, flag)  
    write(*,*) E,real(eta)
    enddo

end program bindenergy_sngl
