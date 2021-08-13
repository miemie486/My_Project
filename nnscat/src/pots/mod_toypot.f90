! Bingwei Long
! 01/09/2014
! Toy-model potentials

module mod_toypot

  use feb_type,       only: standard_error_unit
  use ope_pwd
  use util_gauleg,    only: feb_glabsc_2pt
  implicit none


contains


  ! Partial-wave decomposition of a local potential in spatial coordinates
  ! Units in MeV or MeV^(-1)
  subroutine pwd_coord(p, k, l, Vr, potval)

    real(NER), intent(in)   :: p, k
    integer, intent(in)     :: l
    real(NER), intent(out)  :: potval
    interface
      function Vr(r)
        import NER
        real(NER), intent(in)   :: r
        real(NER)               :: Vr
      end function
    end interface

    integer, parameter      :: ngp = 100
    real(NER)   :: rdn, rup
    integer                 :: ii
    real(NER)               :: ri(ngp), wghti(ngp), reg

    rdn = 0.0_NER
    rup = 50.0_NER/PC_mpi
    call feb_glabsc_2pt(ngp, ri, wghti, rdn, rup)
    reg = 0.0_NER
    do ii = 1, ngp
      select case (l)
        case (0)
          reg = reg + wghti(ii)*sin(ri(ii)*p)*sin(ri(ii)*k)*Vr(ri(ii))
!                case (1)
!                    reg = reg + wghti(ii)*sin(ri(ii)*p)*sin(ri(ii)*k)*Vr(ri(ii))/(p*k)
        case default
          write (standard_error_unit, '(a)') 'pwd_coord: l is too high.'
      end select
    end do
    potval = 4.0_NER*PI_NE*reg/(p*k)

  end subroutine pwd_coord


  ! Singular attractive potential in coordinate space
  function singattrpot_toy(r)

    real(NER), intent(in)   :: r
    real(NER)               :: singattrpot_toy

    real(NER), parameter    :: g = -1.0_NER, M_NN = 50.0_NER, mm = M_NN
    real(NER)               :: Rc

    Rc = 1.0E-2_NER/PC_mpi
    singattrpot_toy = g/(M_NN*PC_mN*(r+Rc)**2*r)*exp(-mm*r)

  end function singattrpot_toy


end module mod_toypot
