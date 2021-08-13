! Bingwei Long 12/14/2010
! procedures of implementing Gauss-Legendre algorithsm
! references: Numerical recipe


module util_gauleg

  use nneft_type
  implicit none

  !EPS is the relative precision
  real(DP), parameter, private :: EPS_D = 3.0E-15_DP
  real(SP), parameter, private :: EPS_S = 3.0E-6_SP

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
    real(NER), dimension(1:ngp), intent(out)    :: xabsc, weig

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


! GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
! with given lower and upper bounds x1, x3
! [x1, x2] is mapped onto [-1, 0] and [x2, x3] onto [0, 1]
! using (A.5) and (A.6) in Jerry's thesis

  subroutine feb_glabsc_3pt(ngp, xabsc, weig, x1, x2, x3)

    integer, intent(in)                         :: ngp  ! # of Gauss Points
    real(NER), intent(in)                       :: x1, x2, x3   ! bounds, middle
    real(NER), dimension(1:ngp), intent(out)    :: xabsc, weig

    integer     :: i, m
    real(NER)   :: xi, reg1, reg2, reg3 ! temporary vars

    call feb_glabsc_std(ngp, xabsc, weig)
    do i = 1, ngp
      xi = xabsc(i)
      reg1 = x1 - x3 + xi * (x1 + x3 - 2.0E0 * x2)
      reg2 = x1 * x3
      reg3 = x2 * (x1 + x3)
      xabsc(i) = (x2 * (x1 - x3) + xi * (2.0E0 * reg2 - reg3)) / reg1
      weig(i) = 2.0E0 * weig(i) * (x1 - x3) * (reg2 + x2 * x2 - reg3) / (reg1 * reg1)
    end do

  end subroutine feb_glabsc_3pt

  ! ! t = exp(-x) where t are GAUSS-LEGENDRE abscissas
  ! subroutine feb_expabsc_2pt(ngp, xabsc, weig, xdn, xup)

  !     ! # of Points
  !     integer, intent(in)                         :: ngp
  !     ! bounds
  !     real(NER), intent(in)                       :: xdn, xup
  !     real(NER), intent(out)    :: xabsc(:), weig(:)

  !     integer :: i, m
  !     real(NER)   :: xi, ti, x1, x2, wxi, tmpxabsc(1:ngp), tmpweig(1:ngp)

  !     if (x1 < 0 .or. x2 < 0) then
  !         write (standard_error_unit, '(a)') 'feb_expabsc_2pt: both upper and lower limits must be positive. Will take their absolute values.'
  !     end if
  !     x1 = abs(xdn)
  !     x2 = abs(xup)
  !     if (x1 > x2) then
  !         write (standard_error_unit, '(a)') 'feb_expabsc_2pt: xdn must be smaller than xup. Will swap them.'
  !     end if
  !     if (x1 > x2) then
  !         xi = x1
  !         x1 = x2
  !         x2 = xi
  !     end if

  !     call feb_glabsc_2pt(ngp, tmpxabsc, tmpweig, exp(-x2), exp(-x1))
  !     do i = 1, ngp
  !         ti = tmpxabsc(i)
  !         wxi = tmpweig(i)/ti
  !         xi = -log(ti)
  !         xabsc(i) = xi
  !         weig(i) = wxi
  !     end do

  ! end subroutine feb_expabsc_2pt

  ! t = log(x) where t are GAUSS-LEGENDRE abscissas
  subroutine feb_logabsc_2pt(ngp, xabsc, weig, xdn, xup)

    integer, intent(in)     :: ngp
    real(NER), intent(in)   :: xdn, xup
    real(NER), intent(out)  :: xabsc(:), weig(:)

    integer :: i, m
    real(NER)   :: xi, ti, x1, x2, wxi, tmpxabsc(1:ngp), tmpweig(1:ngp)

    if (xdn < 0.0_NER .or. xup < 0.0_NER) then
      write (standard_error_unit, '(a)') 'feb_logabsc_2pt: both upper and lower limits must be positive. Will take their absolute values.'
    end if
    x1 = abs(xdn)
    x2 = abs(xup)
    if (x1 > x2) then
      write (standard_error_unit, '(a)') 'feb_logabsc_2pt: xdn must be smaller than xup. Will swap them.'
    end if
    if (x1 > x2) then
      xi = x1
      x1 = x2
      x2 = xi
    end if

    call feb_glabsc_2pt(ngp, tmpxabsc, tmpweig, log(x1), log(x2))
    do i = 1, ngp
      xi = exp(tmpxabsc(i))
      weig(i) = tmpweig(i)*xi
      xabsc(i) = xi
    end do

  end subroutine feb_logabsc_2pt

  ! xi = M*tan(\pi/4*(ti+1)) where ti are GAUSS-LEGENDRE abscissas
  ! wxi = M*\pi/4 * wti/cos(\pi/4*(ti+1))^2
  ! 0 < xi < \infty
  subroutine feb_tanabsc(ngp, M, xabsc, weig)

    integer, intent(in)     :: ngp
    real(NER), intent(in)   :: M
    real(NER), intent(out)  :: xabsc(:), weig(:)

    integer     :: ii
    real(NER)   :: xi, wxi, ti(1:ngp), wti(1:ngp)

    call feb_glabsc_std(ngp, ti, wti)
    forall (ii = 1:ngp)
      xabsc(ii) = M*tan(0.25_NER*PI_NE*(ti(ii)+1.0_NER))
      weig(ii) = M*0.25_NER*PI_NE * wti(ii)              &
        & / (cos(0.25_NER*PI_NE*(ti(ii)+1.0_NER)))**2
    end forall

  end subroutine feb_tanabsc

  subroutine feb_glabsc_cmplx_line(initpos, ngp, zabsc, weig, z1, z2)

    integer, intent(in)                        :: initpos, ngp
    complex(NEC), intent(in)                   :: z1, z2
    complex(NEC), intent(out)    :: zabsc(:), weig(:)

    integer         :: ii
    real(NER)       :: dist, x(1:ngp), w(1:ngp)
    complex(NEC)    :: unitn

    dist = abs(z2 - z1)
    unitn = (z2 - z1)/dist
    call feb_glabsc_2pt(ngp, x, w, 0.0_NER, dist)
    do ii = 1, ngp
      weig(initpos + ii - 1) = unitn*w(ii)
      zabsc(initpos + ii - 1) = z1 + unitn*x(ii)
    end do

  end subroutine feb_glabsc_cmplx_line

  ! subroutine feb_expabsc_cmplx_line(ngp, zabsc, weig, z1, z2, a)

  !     integer, intent(in)         :: ngp
  !     complex(NEC), intent(in)    :: z1, z2
  !     real(NER), intent(in)       :: a
  !     complex(NEC), intent(out)   :: zabsc(:), weig(:)

  !     integer         :: ii
  !     real(NER)       :: dist, x(1:ngp), w(1:ngp)
  !     complex(NEC)    :: unitn

  !     dist = abs(z2 - z1)
  !     unitn = (z2 - z1)/dist
  !     call feb_expabsc_2pt(ngp, x, w, a, dist+a)
  !     do ii = 1, ngp
  !         weig(ii) = unitn*w(ii)
  !         zabsc(ii) = z1 + unitn*(x(ii) - a)
  !     end do

  ! end subroutine feb_expabsc_cmplx_line

  ! z = z1 + unitn*(exp(t) - a)
  subroutine feb_logabsc_cmplx_line(initpos, ngp, zabsc, weig, z1, z2, a)

    integer, intent(in)         :: initpos, ngp
    complex(NEC), intent(in)    :: z1, z2
    real(NER), intent(in)       :: a
    complex(NEC), intent(out)   :: zabsc(:), weig(:)

    integer         :: ii
    real(NER)       :: dist, x(1:ngp), w(1:ngp)
    complex(NEC)    :: unitn

    dist = abs(z2 - z1)
    unitn = (z2 - z1)/dist
    call feb_logabsc_2pt(ngp, x, w, a, dist+a)
    do ii = 1, ngp
      weig(initpos + ii - 1) = unitn*w(ii)
      zabsc(initpos + ii - 1) = z1 + unitn*(x(ii) - a)
    end do

  end subroutine feb_logabsc_cmplx_line


end module util_gauleg
