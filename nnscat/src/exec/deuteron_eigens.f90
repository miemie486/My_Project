!wenchao shi    2019/10/24

program deuteron_eigens

  use wave_func_deute
  use nneft_eigens
  use util_rootfinder
  implicit none


!   integer,parameter   ::    unit_eigenvectorout=102
!   REAL(NER) :: x1,x2
!   REAL(NER) :: bindenergy
 
!   complex(NEC) :: A_matrix(1:2*N,1:2*N)
!   complex(NEC)  :: eta, Gamma(1:2*N)
!   integer       :: flag
!   integer :: dd

!   x1=-3.0_NER
!   x2=-1.0_NER
!   call get_root(eta_,x1,x2,bindenergy)
!   write(*,*) bindenergy

!   call get_A_matrix(bindenergy,A_matrix)
!   call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag)   
!   open(unit=unit_eigenvectorout,file="deuteron/eigenvector.out") 
!   do dd=1,2*N
!   write(unit_eigenvectorout,*) Gamma(dd)
!   end do


!   integer,parameter   ::    unit_eigenvectorout=102
!   REAL(NER) :: x1,x2
!   REAL(NER) :: bindenergy
 
!   complex(NEC) :: A_matrix(1:2*N,1:2*N)
!   complex(NEC)  :: eta, Gamma(1:2*N)
!   integer       :: flag
!   integer :: dd

!   x1=-3.0_NER
!   x2=-1.0_NER
!   call get_root(eta_nowrap,x1,x2,bindenergy)
!   write(*,*) bindenergy

!   call get_A_matrix_nowrap(bindenergy,A_matrix)
!   call get_eigens_matrix_cmplx_modeigens(A_matrix, 2*N, eta, Gamma, flag)  
!   open(unit=unit_eigenvectorout,file="deuteron/eigenvector.out") 
!   do dd=1,2*N
!   write(unit_eigenvectorout,*) Gamma(dd)
!   end do








!!!!!!!   subroutine get_A_matrix_fit(Mambda,C_3s1,E,A_matrix)
!   integer,parameter   ::    unit_C_3s1out=103
  REAL(NER) :: x1,x2
  REAL(NER) :: C_3s1
!     PC_gA = 1.29_NER
!     PC_mN = 939.0_NER
!     PC_mpi = 138.0_NER
!     PC_fpi = 92.4_NER
    
!!!!!! cutoff=500
  x1=-0.52458041537666564E-002
  x2=-0.00458041537666564E-002
  call get_root(eta_C_8,x1,x2,C_3s1)
!   open(unit=unit_C_3s1out,file="deuteron/C_3s1.out")
  write(*,*) "500",C_3s1
end program deuteron_eigens



