! Bingwei Long 12/14/2010
! Complex quadrature

module mod_cmplx_quad

  implicit none

  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))

  integer, parameter  :: NER = DP
  integer, parameter  :: NEC = DPC
  ! integer, parameter  :: NER = SP
  ! integer, parameter  :: NEC = SPC

  integer, parameter :: standard_error_unit = 0
  integer, parameter :: standard_input_unit = 5
  integer, parameter :: standard_output_unit = 6

  !EPS is the relative precision
  real(DP), parameter, private :: EPS_D = 3.0E-15_DP
  real(SP), parameter, private :: EPS_S = 3.0E-6_SP
  real(NER), parameter :: PI_NE      = 3.141592653589793238462643383279502884197_NER
  complex(NEC), parameter :: IMGUNT_NE = (0.0_NER, 1.0_NER)


contains


!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!* For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************

  subroutine feb_glabsc_std(ngp, xabsc, weig)

    integer, intent(in)                         :: ngp  ! # of Gauss Points
    ! real(NER), dimension(1:ngp), intent(out)    :: xabsc, weig
    real(NER), intent(out)    :: xabsc(1:ngp), weig(1:ngp)

    integer     :: i, j, m
    real(NER)    :: p1, p2, p3, pp, z, z1, EPS

    if (ngp < 1) then
      write (standard_error_unit, '(a)')  &
        'feb_glabsc_std: Number of Gaussian points must be at least 1. Nothing will be done.'
      return
    end if
    if (NER == DP) then
      EPS = EPS_D
    else
      EPS = EPS_S
    end if
    m = int((ngp + 1) / 2)
    ! Roots are symmetric in the interval - so only need to find half of them
    ! Loop over the desired roots
    do i = 1, m
      z = cos( PI_NE * (i-0.25) / (ngp+0.5) )
      ! Starting with the above approximation to the ith root,
      ! we enter the main loop of refinement by NEWTON'S method
      do
        p1 = 1.0
        p2 = 0.0
        ! Loop up the recurrence relation to get the Legendre
        ! polynomial evaluated at z
        do j = 1, ngp
          p3 = p2
          p2 = p1
          p1 = ((2.0*j-1.0) * z * p2 - (j-1.0)*p3) / j
        end do

        ! p1 is now the desired Legendre polynomial. We next compute pp,
        ! its derivative, by a standard relation involving also p2, the
        ! polynomial of one lower order.
        pp = ngp*(z*p1-p2)/(z*z-1.0)
        z1 = z
        z  = z1 - p1/pp                               ! Newton's Method

        if (abs(z-z1) < EPS) exit
      end do

      ! Roots will be bewteen -1.0 & 1.0 and symmetric about the origin
      xabsc(i)       =  - z
      xabsc(ngp+1-i) =  + z
      ! Compute the weight and its symmetric counterpar
      weig(i)       = 2.0/((1.0-z*z)*pp*pp)
      weig(ngp+1-i) = weig(i)
    end do                                               ! i loop

  end subroutine feb_glabsc_std


! GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
! with given lower and upper bounds xdn, xup
! [xdn, xup] is *linearly* mapped onto [-1, 1]

  subroutine feb_glabsc_2pt(ngp, xabsc, weig, xdn, xup)

    ! # of Gauss Points
    integer, intent(in)                         :: ngp
    ! bounds
    real(NER), intent(in)                       :: xdn, xup
    real(NER), dimension(1:ngp), intent(out)    :: xabsc, weig

    integer     :: i, m
    real(NER)   :: z, xm, xl

    m = int((ngp + 1) / 2)
    xm = 0.5 * (xup + xdn)
    xl = 0.5 * (xup - xdn)
    call feb_glabsc_std(ngp, xabsc, weig)

    do i = 1, m
      ! Roots will be bewteen xdn & xup and symmetric about the origin
      z = - xabsc(i)
      xabsc(i) = xm - xl * z
      xabsc(ngp+1-i) = xm + xl * z
      ! Compute the weight and its symmetric counterpar
      weig(i) = weig(i) * xl
      weig(ngp+1-i) = weig(i)
    end do                                               ! i loop

  end subroutine feb_glabsc_2pt

  subroutine build_zigzag_mesh(numseg, segN, lpts, mshc, wghtc)

    integer, intent(in) :: numseg, segN(:)
    complex(NEC), intent(in) :: lpts(:)
    complex(NEC), intent(out) :: mshc(:), wghtc(:)

    integer :: ii, is, tN, pos
    real(NER) :: x, w
    real(NER), allocatable  :: mshr(:), wghtr(:)
    complex(NEC)    :: l

    pos = 1
    do is = 1, numseg
      if (allocated(mshr)) deallocate(mshr)
      if (allocated(wghtr)) deallocate(wghtr)
      allocate(mshr(1:segN(is)), wghtr(1:segN(is)))
      call feb_glabsc_2pt(segN(is), mshr, wghtr, 0.0_NER, 1.0_NER)
      do ii = 1, segN(is)
        x = mshr(ii)
        w = wghtr(ii)
        l = sn(lpts(is), lpts(is+1), x)
        mshc(pos) = l
        wghtc(pos) = w*drvtv_sn(lpts(is), lpts(is+1), x)
        pos = pos + 1
      end do
      deallocate(mshr, wghtr)
    end do

  end subroutine build_zigzag_mesh

  complex(NEC) function sn(l1, l2, x)

    complex(NEC), intent(in):: l1, l2
    real(NER), intent(in)   :: x

    sn = x*l2 + (1.0_NER - x)*l1   ! linear in l

  end function sn


  complex(NEC) function drvtv_sn(l1, l2, x)

    complex(NEC), intent(in):: l1, l2
    real(NER), intent(in)   :: x

    drvtv_sn = l2 - l1

  end function drvtv_sn


end module mod_cmplx_quad
