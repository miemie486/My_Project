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
!   call get_root(eta_sngl_nowrap,x1,x2,bindenergy)
!   write(*,*) bindenergy


    real(NER) :: E
    complex(NEC)          :: A_matrix(1:N,1:N)
    complex(NEC)          :: eta, Gamma(1:N)
    integer               :: flag,ll

    do ll=1,30
    read(*,*) E
    call get_A_matrix_sngl_nowrap(E,A_matrix)
    call get_eigens_matrix_cmplx_modeigens(A_matrix, N, eta, Gamma, flag)
    write(*,*) E,real(eta)
    enddo


!   integer,parameter   ::    unit_Elistin=333,unit_Mambda_Cin=666,unit_E_etaout=999
!     real(NER)  :: Mambda(12), para(0:2)
!     complex(NEC)          :: A_matrix(1:N,1:N)
!     complex(NEC)          :: eta, Gamma(1:N)
!     integer               :: flag,ll,yy
!     open(unit=unit_E_etaout,file="sngl/E_eta.out")
!     open(unit=unit_Mambda_Cin,file="sngl/Mambda_C.in")

!     do ll=1,12
!     read(unit_Mambda_Cin,*) Mambda(ll),para(1),para(2)
!     write(unit_E_etaout,*) "# Mambda=", Mambda(ll)
!     open(unit=unit_Elistin,file="sngl/Elist.in")
!     do yy=1,30
!     read(unit_Elistin,*) para(0)
!     call get_A_matrix_sngl_nowrap_Mambda_C(Mambda(ll), E, para, A_matrix)
!     call get_eigens_matrix_cmplx_modeigens(A_matrix, N, eta, Gamma, flag)
!     write(unit_E_etaout,*) para(0),real(eta)
!     end do
!     close(unit_Elistin)
!     end do

end program bindenergy_sngl
