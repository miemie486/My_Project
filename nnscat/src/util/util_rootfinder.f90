! Bingwei Long 06/15/2017
! Bingwei Long 06/18/2011
! NR subroutines to find roots

module util_rootfinder

  use nneft_type
  implicit none

  integer, parameter, private :: NPAR_ARTH=16, NPAR2_ARTH=8
  integer, parameter  :: TWOVARS_MAXITER_REACHED = 1, TWOVARS_MAXITER_SUCCESSFUL = 0, &
  TWOVARS_MAXITER_OTHER = 2

  contains

  subroutine findroot_2vars(f, g, x0, y0, eps, maxit, x, y, iter, ierror)

    use nneft_type
    implicit none

    integer, intent(in)     :: maxit
    integer, intent(out)    :: iter, ierror
    real(NER), intent(in)   :: x0, y0, eps
    real(NER), intent(out)  :: x, y
    interface
      function f(x1, x2)
        use nneft_type
        implicit none
        real(NER), intent(in)   :: x1, x2
        real(NER)   :: f
      end function f
      function g(x1, x2)
        use nneft_type
        implicit none
        real(NER), intent(in)   :: x1, x2
        real(NER)   :: g
      end function g
    end interface

    real(NER)   ::  a, b, c, d, t, xm, xn, p, q, h
    real(NER), parameter    :: local_EPS = 1.0E-6_NER

    h = 0.01_NER
    ierror=TWOVARS_MAXITER_REACHED; iter=0
    x=x0; y=y0
  100 iter=iter+1
    if (iter > maxit) return
    a=f(x+h,y)
    b=g(x+h,y)
    a=(a-f(x-h,y))*0.5_NER/h
    b=(b-g(x-h,y))*0.5_NER/h
    c=f(x,y+h)
    d=g(x,y+h)
    c=(c-f(x,y-h))*0.5_NER/h
    d=(d-g(x,y-h))*0.5_NER/h
    t=a*d-b*c
    if (abs(t)<local_EPS) then
      ierror=TWOVARS_MAXITER_OTHER
      return
    end if
    xm=f(x,y)
    xn=g(x,y)
    p=(xm*d-xn*c)/t
    q=(xn*a-xm*b)/t
    x=x-p; y=y-q
    if (abs(p)+abs(q) > eps) goto 100
    ierror=TWOVARS_MAXITER_SUCCESSFUL
    return
  end subroutine findroot_2vars

  subroutine findroot_cmplx_var(f, z0, eps, maxit, z, iter, ierror)

    use nneft_type
    implicit none

    integer, intent(in)     :: maxit
    integer, intent(out)    :: iter, ierror
    ! real(NER), intent(in)   :: x0, y0, eps
    ! real(NER), intent(out)  :: x, y
    real(NER), intent(in)       :: eps
    complex(NEC), intent(in)    :: z0
    complex(NEC), intent(out)   :: z
    interface
      function f(z)
        use nneft_type
        implicit none
        complex(NEC), intent(in)   :: z
        complex(NEC)   :: f
      end function f
      ! function g(x1, x2)
      !     use nneft_type
      !     implicit none
      !     real(NER), intent(in)   :: x1, x2
      !     real(NER)   :: g
      ! end function g
    end interface

    real(NER)   ::  a, b, c, d, t, xm, xn, p, q, h
    real(NER), parameter    :: local_EPS = 1.0E-9_NER
    complex(NEC)    :: w, wxm, wxp, wym, wyp

    h = 0.01_NER
    ierror=TWOVARS_MAXITER_REACHED; iter=0
    ! x=x0; y=y0
    z = z0
  100 iter=iter+1
    if (iter > maxit) return
    ! a=f(x+h,y)
    ! b=g(x+h,y)
    wxp = f(z+h)
    ! a=(a-f(x-h,y))*0.5_NER/h
    ! b=(b-g(x-h,y))*0.5_NER/h
    wxm = f(z-h)
    a = (real(wxp - wxm))*0.5_NER/h
    b = (aimag(wxp - wxm))*0.5_NER/h
    ! c=f(x,y+h)
    ! d=g(x,y+h)
    wyp = f(z+IMGUNT_NE*h)
    ! c=(c-f(x,y-h))*0.5_NER/h
    ! d=(d-g(x,y-h))*0.5_NER/h
    wym = f(z-IMGUNT_NE*h)
    c = (real(wyp - wym))*0.5_NER/h
    d = (aimag(wyp - wym))*0.5_NER/h
    t=a*d-b*c
    if (abs(t)<local_EPS) then
      ierror=TWOVARS_MAXITER_OTHER
      return
    end if
    ! xm=f(x,y)
    ! xn=g(x,y)
    w = f(z)
    xm = real(w)
    xn = aimag(w)
    p=(xm*d-xn*c)/t
    q=(xn*a-xm*b)/t
    ! x=x-p; y=y-q
    z = z - cmplx(p, q, NEC)
    if (abs(p)+abs(q) > eps) goto 100
    ierror=TWOVARS_MAXITER_SUCCESSFUL
    return
  end subroutine findroot_cmplx_var

! Given a function func and an initial guessed range x1 to x2, the routine expands the range
! geometrically until a root is bracketed by the returned values x1 and x2 (in which case
! succes returns as .true.) or until the range becomes unacceptably large (in which case
! succes returns as .false.).

  subroutine zbrac(func,x1,x2,succes)

    USE nneft_type
    IMPLICIT NONE

    REAL(NER), INTENT(INOUT)    :: x1,x2
    logical, INTENT(OUT)        :: succes
    INTERFACE
      FUNCTION func(x)
        USE nneft_type
        IMPLICIT NONE
        REAL(NER), INTENT(IN)   :: x
        REAL(NER)               :: func
      END FUNCTION func
    END INTERFACE
    integer, PARAMETER      :: NTRY=500
    REAL(NER), PARAMETER    :: FACTOR=1.6
    integer :: j
    REAL(NER) :: f1,f2

    f1=func(x1)
    f2=func(x2)
    succes=.true.
    do j=1,NTRY
      if ((f1 > 0.0 .and. f2 < 0.0) .or. (f1 < 0.0 .and. f2 > 0.0)) RETURN
      if (abs(f1) < abs(f2)) then
        x1=x1+FACTOR*(x1-x2)
        f1=func(x1)
      else
        x2=x2+FACTOR*(x2-x1)
        f2=func(x2)
      end if
    end do
    succes=.false.

  END SUBROUTINE zbrac


! Given a function func defined on the interval from x1-x2 subdivide the interval into n
! equally spaced segments, and search for zero crossings of the function. nb is returned as
! the number of bracketing pairs xb1(1:nb), xb2(1:nb) that are found. xb1 and xb2 are
! pointers to arrays of length nb that are dynamically allocated by the routine.
  SUBROUTINE zbrak(func,x1,x2,n,xb1,xb2,nb)
    USE nneft_type
    IMPLICIT NONE
    integer, INTENT(IN) :: n
    integer, INTENT(OUT) :: nb
    REAL(NER), INTENT(IN) :: x1,x2
    REAL(NER), allocatable:: xb1(:), xb2(:)
    INTERFACE
      FUNCTION func(x)
        USE nneft_type
        IMPLICIT NONE
        REAL(NER), INTENT(IN) :: x
        REAL(NER) :: func
      END FUNCTION func
    END INTERFACE
    integer :: i
    REAL(NER) :: dx
    REAL(NER), DIMENSION(0:n) :: f,x
    logical, DIMENSION(1:n) :: mask
    ! logical, SAVE :: init=.true.

    ! if (init) then
    !     init=.false.
    !     nullify(xb1,xb2)
    ! end if
    if (allocated(xb1)) deallocate(xb1)
    if (allocated(xb2)) deallocate(xb2)
    dx=(x2-x1)/n     ! Determine the spacing appropriate to the mesh.
    x=x1+dx*arth_i(0,1,n+1)
    do i=0,n         ! Evaluate the function at the mesh points.
      f(i)=func(x(i))
    end do
    mask=f(1:n)*f(0:n-1) <= 0.0   ! Record where the sign changes occur.
    nb=count(mask)                ! Number of sign changes.
    allocate(xb1(nb),xb2(nb))
    xb1(1:nb)=pack(x(0:n-1),mask) ! Store the bounds of each bracket.
    xb2(1:nb)=pack(x(1:n),mask)
  END SUBROUTINE zbrak

  SUBROUTINE zbrak_pointer(func,x1,x2,n,xb1,xb2,nb)
    USE nneft_type
    IMPLICIT NONE
    integer, INTENT(IN) :: n
    integer, INTENT(OUT) :: nb
    REAL(NER), INTENT(IN) :: x1,x2
    REAL(NER), DIMENSION(:), POINTER :: xb1,xb2
    INTERFACE
      FUNCTION func(x)
        USE nneft_type
        IMPLICIT NONE
        REAL(NER), INTENT(IN) :: x
        REAL(NER) :: func
      END FUNCTION func
    END INTERFACE
    integer :: i
    REAL(NER) :: dx
    REAL(NER), DIMENSION(0:n) :: f,x
    logical, DIMENSION(1:n) :: mask
    logical, SAVE :: init=.true.

    if (init) then
      init=.false.
      nullify(xb1,xb2)
    end if
    if (associated(xb1)) deallocate(xb1)
    if (associated(xb2)) deallocate(xb2)
    dx=(x2-x1)/n     ! Determine the spacing appropriate to the mesh.
    x=x1+dx*arth_i(0,1,n+1)
    do i=0,n         ! Evaluate the function at the mesh points.
      f(i)=func(x(i))
    end do
    mask=f(1:n)*f(0:n-1) <= 0.0   ! Record where the sign changes occur.
    nb=count(mask)                ! Number of sign changes.
    allocate(xb1(nb),xb2(nb))
    xb1(1:nb)=pack(x(0:n-1),mask) ! Store the bounds of each bracket.
    xb2(1:nb)=pack(x(1:n),mask)
  END SUBROUTINE zbrak_pointer

! Using Brent’s method, find the root of a function func known to lie between x1 and x2.
! The root, returned as zbrent, will be refined until its accuracy is tol.
! Parameters: Maximum allowed number of iterations, and machine floating-point precision.
  FUNCTION zbrent(func,x1,x2,tol)

    USE nneft_type
    IMPLICIT NONE
    REAL(NER), INTENT(IN) :: x1,x2,tol
    REAL(NER) :: zbrent
    INTERFACE
      FUNCTION func(x)
      USE nneft_type
      IMPLICIT NONE
      REAL(NER), INTENT(IN) :: x
      REAL(NER) :: func
      END FUNCTION func
    END INTERFACE
    integer, PARAMETER :: ITMAX=100
    REAL(NER), PARAMETER :: EPS=epsilon(x1)
    integer :: iter
    REAL(NER) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

    a=x1
    b=x2
    fa=func(a)
    fb=func(b)
    if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
      zbrent = 0.0
      return
    end if
    c=b
    fc=fb
    do iter=1,ITMAX
      if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
        c=a                              ! Rename a, b, c and adjust bounding interval d.
        fc=fa
        d=b-a
        e=d
      end if
      if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
      end if
      tol1=2.0*EPS*abs(b)+0.5*tol     ! Convergence check.
      xm=0.5*(c-b)
      if (abs(xm) <= tol1 .or. fb == 0.0) then
        zbrent=b
        RETURN
      end if
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s=fb/fa                            ! Attempt inverse quadratic interpolation.
        if (a == c) then
          p=2.0*xm*s
          q=1.0-s
        else
          q=fa/fc
          r=fb/fc
          p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
          q=(q-1.0)*(r-1.0)*(s-1.0)
        end if
        if (p > 0.0) q=-q                 ! Check whether in bounds.
        p=abs(p)
        if (2.0*p < min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
          e=d                       ! Accept interpolation.
          d=p/q
        else
          d=xm                      ! Interpolation failed; use bisection.
          e=d
        end if
      else                              ! Bounds decreasing too slowly; use bisection.
        d=xm
        e=d
      end if
        a=b                               ! Move last best guess to a.
        fa=fb
        b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )   ! Evaluate new trial root.
        fb=func(b)
    end do
    zbrent=b

  END FUNCTION zbrent

! Using the secant method, find the root of a function func thought to lie between x1 and
! x2 . The root, returned as rtsec , is refined until its accuracy is ± xacc .

  FUNCTION rtsec(func,x1,x2,xacc)

    USE nneft_type
    IMPLICIT NONE
    REAL(NER), INTENT(IN) :: x1,x2,xacc
    REAL(NER) :: rtsec
    INTERFACE
      FUNCTION func(x)
      USE nneft_type
      IMPLICIT NONE
      REAL(NER), INTENT(IN) :: x
      REAL(NER) :: func
      END FUNCTION func
    END INTERFACE
    INTEGER, PARAMETER :: MAXIT=100  ! Parameter: MAXIT is the maximum allowed number of iterations.
    INTEGER :: j
    REAL(NER) :: dx,f,fl,xl,ftmp

    fl=func(x1)
    f=func(x2)
    if (abs(fl) < abs(f)) then
! Pick the bound with the smaller function value as the most recent guess.
      rtsec=x1
      xl=x2
      ftmp=fl
      fl=f
      f=ftmp
    else
      xl=x1
      rtsec=x2
    end if
    do j=1,MAXIT
      dx=(xl-rtsec)*f/(f-fl)
      xl=rtsec
      fl=f
      rtsec=rtsec+dx
      f=func(rtsec)
      if (abs(dx) < xacc .or. f == 0.0) RETURN
    end do
    write(*,*) 'rtsec: exceed maximum iterations, STOP'
    STOP
  END FUNCTION rtsec

  subroutine smart_brak(func, x1, x2, n, tol, xroots, numroots)
    INTERFACE
      FUNCTION func(x)
        import :: NER
        IMPLICIT NONE
        REAL(NER), INTENT(IN) :: x
        REAL(NER) :: func
      END FUNCTION func
    END INTERFACE
    real(NER), intent(in)   :: x1, x2, tol
    integer, intent(in)     :: n
    real(NER), intent(out)  :: xroots(:)
    integer, intent(out)    :: numroots

    real(NER), allocatable    :: xb1(:), xb2(:)
    integer :: nb, ii

    call zbrak(func, x1, x2, n, xb1, xb2, nb)
    if (nb == 0) then
      numroots = 0
      if (allocated(xb1)) deallocate(xb1)
      if (allocated(xb2)) deallocate(xb2)
      return
    end if
    do ii = 1, nb
      xroots(ii) = zbrent(func, xb1(ii), xb2(ii), tol)
    end do
    numroots = nb
    if (allocated(xb1)) deallocate(xb1)
    if (allocated(xb2)) deallocate(xb2)

  end subroutine smart_brak


  FUNCTION arth_d(first,increment,n)
    REAL(NER), INTENT(IN) :: first,increment
    integer, INTENT(IN) :: n
    REAL(NER), DIMENSION(n) :: arth_d
    integer :: k,k2
    REAL(NER) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
      do k=2,n
        arth_d(k)=arth_d(k-1)+increment
      end do
    else
      do k=2,NPAR2_ARTH
        arth_d(k)=arth_d(k-1)+increment
      end do
      temp=increment*NPAR2_ARTH
      k=NPAR2_ARTH
      do
        if (k >= n) exit
        k2=k+k
        arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
        temp=temp+temp
        k=k2
      end do
    end if
  END FUNCTION arth_d


  FUNCTION arth_i(first,increment,n)
    integer, INTENT(IN) :: first,increment,n
    integer, DIMENSION(n) :: arth_i
    integer :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
      do k=2,n
        arth_i(k)=arth_i(k-1)+increment
      end do
    else
      do k=2,NPAR2_ARTH
        arth_i(k)=arth_i(k-1)+increment
      end do
      temp=increment*NPAR2_ARTH
      k=NPAR2_ARTH
      do
        if (k >= n) exit
        k2=k+k
        arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
        temp=temp+temp
        k=k2
      end do
    end if
  END FUNCTION arth_i


end module util_rootfinder
