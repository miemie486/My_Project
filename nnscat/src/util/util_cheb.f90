! Bingwei Long 10/24/2013
! Utility subroutines
! Chebyshev coefficients and polynomial approximations

module util_cheb

  use nneft_type
  use util_nr
  implicit none

  abstract interface
      function func_real(x)
        use nneft_type
        implicit none
        real(NER), intent(in)   :: x
        real(NER)               :: func_real
      end function func_real
      function func_cmplx(x)
        use nneft_type
        implicit none
        real(NER), intent(in)   :: x
        complex(NEC)            :: func_cmplx
      end function func_cmplx
      ! multi-valued function
      subroutine func_multi(x, rslts)
        use nneft_type
        implicit none
        real(NER), intent(in)                   :: x
        real(NER), dimension(:), intent(out)    :: rslts
      end subroutine func_multi
      ! multi-valued and 2-variable function
      subroutine func_multi_xy(x, y, rslts)
        use nneft_type
        implicit none
        real(NER), intent(in)                   :: x, y
        real(NER), dimension(:), intent(out)    :: rslts
      end subroutine func_multi_xy
  end interface

  interface chebft
    module procedure chebft_real, chebft_multi_real
  end interface chebft
  interface cheb_taylor
    module procedure cheb_taylor_real, cheb_taylor_cmplx, cheb_taylor_multi, cheb_taylor_multi_xy
  end interface cheb_taylor

contains


  function dfridr(func, x, h, err)
    !~ use nrtype; use nrutil, only : assert, geop, iminloc

    real(NER), intent(in)   :: x, h
    real(NER), intent(out)  :: err
    real(NER)               :: dfridr

    interface
      function func(x)
        import :: NER
        implicit none
        real(NER), intent(in) :: x
        real(NER) :: func
      end function func
    end interface

    integer, parameter :: ntab=10
    real(NER), parameter :: con=1.4_NER, con2=con*con, big=huge(x), safe=2.0_NER
    integer :: ierrmin, i, j
    real(NER) :: hh
    real(NER), dimension(ntab-1) :: errt, fac
    real(NER), dimension(ntab, ntab) :: a
    if (h == 0.0_NER) then
      write (standard_error_unit, '(a)') 'dfridr: h is too small'
      dfridr = 0.0_NER
      return
    end if
    hh=h
    a(1, 1)=(func(x+hh)-func(x-hh))/(2.0_NER*hh)
    err=big
    fac(1:ntab-1)=geop(con2, con2, ntab-1)
    do i=2, ntab
      hh=hh/con
      a(1, i)=(func(x+hh)-func(x-hh))/(2.0_NER*hh)
      do j=2, i
        a(j, i)=(a(j-1, i)*fac(j-1)-a(j-1, i-1))/(fac(j-1)-1.0_NER)
      end do
      errt(1:i-1)=max(abs(a(2:i, i)-a(1:i-1, i)), abs(a(2:i, i)-a(1:i-1, i-1)))
      ierrmin=iminloc(errt(1:i-1))
      if (errt(ierrmin) <= err) then
        err=errt(ierrmin)
        dfridr=a(1+ierrmin, i)
      end if
      if (abs(a(i, i)-a(i-1, i-1)) >= safe*err) return
    end do

  end function dfridr


  ! Taylor expansion a real function around x = 0, return the coefficients
  ! g(i) such that f(x) = \sum_{i=1}^n g(i)*x^{i-1}. [a, b] is the interval
  ! in which the expansion is expected to approximate the original function.

  subroutine cheb_taylor_real(a, b, n, g, func)

    real(NER), intent(in)                   :: a, b
    integer, intent(in)                     :: n
    real(NER), dimension(:), intent(out)    :: g
    procedure(func_real)                    :: func

    real(NER), dimension(n) :: c, d

    call chebft_real(a, b, n, c, func)
    call chebpc(c, d, n)
    call pcshft(a, b, d, g, n)

  end subroutine cheb_taylor_real


  ! Taylor expansion of a multi-valued function f_k around x = 0, return the
  ! coefficients g(k, i) such that f_k(x) = \sum_{i=1}^n g(k, i)*x^{i-1}

  subroutine cheb_taylor_multi(a, b, m, n, g, func)

    real(NER), intent(in) :: a,b
    integer, intent(in) :: m, n
    real(NER), dimension(:, :), intent(out) :: g
    procedure(func_multi)   :: func

    integer :: ii
    real(NER), dimension(1:n) :: d
    real(NER), dimension(1:m, 1:n) :: c

    call chebft_multi_real(a, b, m, n, c, func)
    do ii = 1, m
      call chebpc(c(ii, 1:n), d, n)
      call pcshft(a, b, d, g(ii, 1:n), n)
    end do

  end subroutine cheb_taylor_multi


  ! Taylor expansion of a complex valued function around x = 0, return the
  ! coefficients g(i) such that f(x) = \sum_{i=1}^n g(i)*x^{i-1}

  subroutine cheb_taylor_cmplx(a, b, n, g, func)

    real(NER), intent(in)                   :: a,b
    integer, intent(in)                     :: n
    complex(NEC), dimension(:), intent(out) :: g
    procedure(func_cmplx)                   :: func

    real(NER), dimension(1:2, 1:n)  :: gri

    call cheb_taylor_multi(a, b, 2, n, gri, func_re_im)
    g(1:n) = cmplx(gri(1, 1:n), gri(2, 1:n), NEC)

  contains

    subroutine func_re_im(x, rslts)

      real(NER), intent(in)               :: x
      real(NER), dimension(:), intent(out):: rslts

      complex(NEC)    :: zz

      zz = func(x)
      rslts(1) = real(zz)
      rslts(2) = aimag(zz)

    end subroutine func_re_im

  end subroutine cheb_taylor_cmplx


  ! Taylor expansion of multi-valued, 2-variable function f_k(x, y) around
  ! x, y = 0. Return g(k, i, j) such that
  ! f_k(x, y) = \sum_{i=1}^{n_1} \sum_{j=1}^{n_2} g(k, i, j) x^{i-1} y^{j-1}

  subroutine cheb_taylor_multi_xy(ax, bx, n1, ay, by, n2, m, g, func)

    real(NER), intent(in)                       :: ax, bx, ay, by
    integer, intent(in)                         :: m, n1, n2
    real(NER), dimension(:, :, :), intent(out)  :: g
    procedure(func_multi_xy)                    :: func

    integer                                 :: ii, jj
    real(NER)                               :: x_public
    real(NER), dimension(1:n1)              :: d
    real(NER), dimension(1:m*n2)            :: xi
    real(NER), dimension(1:m*n2, 1:n1)      :: c

    call chebft_multi_real(ax, bx, m*n2, n1, c, func_xi)
    do ii = 1, m
      do jj = 1, n2
        call chebpc(c((ii-1)*n2+jj, 1:n1), d, n1)
        call pcshft(ax, bx, d, g(ii, 1:n1, jj), n1)
      end do
    end do

  contains

    subroutine func_xi(x, fnval)

      real(NER), intent(in) :: x
      real(NER), dimension(:), intent(out) :: fnval

      integer :: fx_ii, fx_jj
      real(NER), dimension(1:m, 1:n2) :: xi_x

      x_public = x
      call cheb_taylor_multi(ay, by, m, n2, xi_x, func_fixed_x)
      do fx_ii = 1, m
        do fx_jj = 1, n2
          fnval((fx_ii-1)*n2+fx_jj) = xi_x(fx_ii, fx_jj)
        end do
      end do

    end subroutine func_xi

    subroutine func_fixed_x(y, ff_fnval)

      real(NER), intent(in) :: y
      real(NER), dimension(:), intent(out) :: ff_fnval

      call func(x_public, y, ff_fnval)

    end subroutine func_fixed_x

  end subroutine cheb_taylor_multi_xy


  ! chebyshev fit: given a function func, lower and upper limits of the
  ! interval [a,b], and a maximum degree n, this routine computes the n
  ! coefficients c_k. This routine is to be used with moderately large n
  ! (e.g., 30 or 50), the array of c’s subsequently to be truncated at the
  ! smaller value m such that cm+1 and subsequent elements are negligible.

  subroutine chebft_real(a, b, n, c, func)

    real(NER), intent(in) :: a,b
    integer, intent(in) :: n
    real(NER), dimension(:), intent(out) :: c
    procedure(func_real)   :: func

    real(NER) :: bma,bpa
    real(NER), dimension(n) :: theta, fn
    integer :: ii

    bma=0.5_NER*(b-a)
    bpa=0.5_NER*(b+a)
    theta(:)=PI_NE*arth(0.5_NER, 1.0_NER,n)/n
    do ii = 1, n
      fn(ii) = func(real(cos(theta(ii))*bma+bpa, NER))
    end do
    c(:)=matmul(cos(outerprod(arth(0.0_NER,1.0_NER,n), n, theta, n)), &
      fn(:))*2.0_NER/n

  end subroutine chebft_real


  ! chebft for multiple valued functions
  ! n : degrees; m : number of functions

  subroutine chebft_multi_real(a, b, m, n, c, func)

    real(NER), intent(in)   :: a,b
    integer, intent(in)     :: m, n
    real(NER), dimension(:, :), intent(out) :: c
    procedure(func_multi)   :: func

    real(NER) :: bma,bpa
    real(NER), dimension(n) :: theta
    real(NER), dimension(1:m, 1:n)  :: fnval
    integer :: ii, jj

    bma=0.5_NER*(b-a)
    bpa=0.5_NER*(b+a)
    theta(:)=PI_NE*arth(0.5_NER, 1.0_NER,n)/n
    do ii = 1, n
      call func(real(cos(theta(ii))*bma+bpa, NER), fnval(1:m, ii))
    end do
    do ii = 1, m
      c(ii, 1:n)=matmul(cos(outerprod(arth(0.0_NER,1.0_NER,n), n, theta, n)), &
        fnval(ii, 1:n))*2.0_NER/n
    end do

  end subroutine chebft_multi_real


  ! chebyshev polynomial coefficients. given a coefficient array c(1:n) of
  ! length n, this routine geNERates a coefficient array d(1:n)

  subroutine chebpc(c, d, n)

    real(NER), dimension(:), intent(in) :: c
    real(NER), dimension(:), intent(out) :: d
    integer, intent(in)                 :: n

    integer :: j
    real(NER), dimension(n) :: dd,sv

    d=0.0_NER
    dd=0.0_NER
    d(1)=c(n)
    do j=n-1,2,-1
      sv(2:n-j+1)=d(2:n-j+1)
      d(2:n-j+1)=2.0_NER*d(1:n-j)-dd(2:n-j+1)
      dd(2:n-j+1)=sv(2:n-j+1)
      sv(1)=d(1)
      d(1)=-dd(1)+c(j)
      dd(1)=sv(1)
    end do
    d(2:n)=d(1:n-1)-dd(2:n)
    d(1)=-dd(1)+0.5_NER*c(1)

  end subroutine chebpc


  ! polynomial coefficient shift. given a coefficient array d(1:n), this
  ! routine geNERates a coefficient array g(1:n), where g(k) is the (k-1)-th
  ! derivative.

  subroutine pcshft(a, b, d, g, n)

    real(NER), intent(in) :: a, b
    real(NER), dimension(:), intent(in) :: d
    real(NER), dimension(:), intent(out) :: g
    integer, intent(in)     :: n

    integer :: j
    real(NER), dimension(n) :: dd
    real(NER) :: x

    dd=d*geop(1.0_NER,2.0_NER/(b-a),n)
    x=-0.5_NER*(a+b)
    g(1)=dd(n)
    g(2:n)=0.0_NER
    do j=n-1,1,-1
      g(2:n+1-j)=g(2:n+1-j)*x+g(1:n-j)
      g(1)=g(1)*x+dd(j)
    end do

  end subroutine pcshft


end module util_cheb
