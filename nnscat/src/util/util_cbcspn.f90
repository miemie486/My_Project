!Wencho Shi    2019/11/3


module util_cbcspn
  use nneft_type

contains

SUBROUTINE spline(x,y,n,yp1,ypn,y2)

  INTEGER n,NMAX
  REAL(NER) yp1,ypn,x(n),y(n),y2(n)
  PARAMETER (NMAX=500)
     ! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
     ! x1 < x2 < ::: < xN, and given values yp1 and ypn for the rst derivative of the interpolating
     ! function at points 1 and n, respectively, this routine returns an array y2(1:n) of
     ! length n which contains the second derivatives of the interpolating function at the tabulated
     ! points xi. If yp1 and/or ypn are equal to 1 1030 or larger, the routine is signaled to set
     ! the corresponding boundary condition for a natural spline, with zero second derivative on
     ! that boundary.
     ! Parameter: NMAX is the largest anticipated value of n.
    INTEGER i,k
    REAL(NER) p,qn,sig,un,u(NMAX)

   if (yp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
   else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.
    y2(i)=(sig-1.)/p
    u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (ypn.gt..99e30) then
      qn=0.
      un=0.
    else
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do  k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
    enddo
    return
 END

 SUBROUTINE splint(xa,ya,y2a,n,x,y)
    INTEGER  n
    REAL(NER) x,y,xa(n),y2a(n),ya(n)
   ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
   ! xai 's in order), and given the array y2a(1:n), which is the output from spline above,
   ! and given a value of x, this routine returns a cubic-spline interpolated value y.
    INTEGER   k,khi,klo
    REAL(NER) a,b,h
    klo=1
    khi=n

    1 if (khi-klo.gt.1) then
     k=(khi+klo)/2
    if(xa(k).gt.x)then
     khi=k
    else
     klo=k
    endif
    goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) write(*,*) 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    return
  END


end module util_cbcspn
