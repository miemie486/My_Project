! Bingwei Long 07/03/2011

module util_spcfnc
! Special functions

    use nneft_type
    implicit none


    contains

    ! Legendre functions
    ! lgndr_n(n, x) = P_n(x)
    function lgndr_n(n, x)

        real(NER)               :: lgndr_n
        integer, intent(in)     :: n
        real(NER), intent(in)   :: x

        integer     :: ii
        real(NER)   :: p0, p1, pf

        p0 = 1.0_NER
        p1 = x
        pf = x
        do ii = 2, n
            pf = ((2.0_NER*ii-1.0_NER)*x*p1 - (ii-1.0_NER)*p0)/ii
            p0 = p1
            p1 = pf
        end do
        lgndr_n = pf

    end function lgndr_n

    ! Legendre functions evaluated at x. Return an array:
    ! pnarray = {P_0(x), P_1(x), P_2(x), ... , P_n(x)}
    subroutine lgndr_array(n, x, pnarray)

        integer, intent(in)         :: n
        real(NER), intent(in)       :: x
        real(NER), dimension(0:)    :: pnarray

        integer     :: ii
        real(NER)   :: p0, p1, pf

        pnarray(0) = 1.0_NER
        pnarray(1) = x
        p0 = 1.0_NER
        p1 = x
        do ii = 2, n
            pf = 2.0_NER*x*p1 - p0 - (x*p1 - p0)/ii
            pnarray(ii) = pf
            p0 = p1
            p1 = pf
        end do

    end subroutine lgndr_array


end module util_spcfnc
