! Bingwei Long 10/24/2013
! Utility subroutines from Numerical Recipies

module util_nr

    use nneft_type
    implicit none

    INTEGER, PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8
    INTEGER, PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2

contains

    FUNCTION outerprod(a, n, b, m)

        REAL(NER), DIMENSION(:), INTENT(IN) :: a, b
        integer, intent(in)                 :: n, m
        REAL(NER), DIMENSION(1:n, 1:m) :: outerprod

        outerprod = spread(a,dim=2,ncopies=m) * &
            spread(b,dim=1,ncopies=n)

    END FUNCTION outerprod


    FUNCTION arth(first,increment,n)

        REAL(NER), INTENT(IN) :: first,increment
        INTEGER, INTENT(IN) :: n
        REAL(NER), DIMENSION(n) :: arth

        INTEGER :: k,k2
        REAL(NER) :: temp

        if (n > 0) arth(1)=first
        if (n <= NPAR_ARTH) then
            do k=2,n
                arth(k)=arth(k-1)+increment
            end do
        else
            do k=2,NPAR2_ARTH
                arth(k)=arth(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                if (k >= n) exit
                k2=k+k
                arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
                temp=temp+temp
                k=k2
            end do
        end if

    END FUNCTION arth


    FUNCTION geop(first,factor,n)

        REAL(NER), INTENT(IN) :: first,factor
        INTEGER, INTENT(IN) :: n
        REAL(NER), DIMENSION(n) :: geop

        INTEGER :: k,k2
        REAL(NER) :: temp
        if (n > 0) geop(1)=first
        if (n <= NPAR_GEOP) then
            do k=2,n
                geop(k)=geop(k-1)*factor
            end do
        else
            do k=2,NPAR2_GEOP
                geop(k)=geop(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                if (k >= n) exit
                k2=k+k
                geop(k+1:min(k2,n))=temp*geop(1:min(k,n-k))
                temp=temp*temp
                k=k2
            end do
        end if

    END FUNCTION geop


    function iminloc(arr)

        real(NER), dimension(:), intent(in) :: arr

        integer, dimension(1) :: imin
        integer :: iminloc

        imin=minloc(arr(:))
        iminloc=imin(1)

    end function iminloc


end module util_nr
