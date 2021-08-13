! Bingwei Long
! Subroutines for interpolation and extrapolation

module util_inter_extra

    use nneft_type
    implicit none


    INTERFACE outerdiff
        MODULE PROCEDURE outerdiff_r, outerdiff_i
    END INTERFACE outerdiff


contains


    ! Given same-size arrays x and y containing a tabulated function yi = f(xi),
    ! this routine returns a same-size array of coefficients c_j, such that
    ! yi = \sum_j c_j x^(j−1)

    subroutine polcoe(n, x, y, c)

        integer, intent(in)                     :: n
        real(NER), dimension(1:n), intent(in)   :: x, y
        real(NER), dimension(1:n), intent(out)  :: c


        integer :: i,k
        real(NER), dimension(n) :: s
        real(NER), dimension(n, n) :: a

        s    = 0.0
        s(n) = -x(1)
        do i=2,n
            s(n+1-i:n-1)=s(n+1-i:n-1)-x(i)*s(n+2-i:n)
            s(n)=s(n)-x(i)
        end do
        a=outerdiff(x, x)
        c = product(a,dim=2,mask=a /= 0.0)
        a(:,1)=-s(1)/x(:)
        do k=2,n
            a(:,k)=-(s(k)-a(:,k-1))/x(:)
        end do
        s = y/c
        c = matmul(s, a)

    END subroutine polcoe


    FUNCTION outerdiff_r(a,b)
    real(NER), dimension(:), intent(in) :: a,b
    real(NER), dimension(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
        spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerdiff_r


    FUNCTION outerdiff_i(a,b)
    INTEGER, dimension(:), intent(in) :: a,b
    INTEGER, dimension(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
        spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerdiff_i


end module util_inter_extra
