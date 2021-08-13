! Bingwei Long 12/05/2010
! subroutines to a number of small things

module util_gadgets

    use nneft_type
    use nneft_phyconst
    implicit none


    contains


    function function_zero()

        real(NER) :: function_zero

        function_zero = 0.0

    end function function_zero


! ! convert threshold T(p, 0; 0)/0 to \alpha in fm^3
! ! both T and p are in MeV

!     function t2alpha_fm3(T, k)

!         real(NER)                   :: t2alpha_fm3
!         complex(NEC), intent(in)    :: T
!         real(NER), intent(in)       :: k

!         t2alpha_fm3 = real(T) / (k*k) * PC_mN * PC_hbarc ** 3

!     end function t2alpha_fm3

!     ! Compute TLab from kcm (relativisitic)
!     function k2Tlab(k)

!         real(NER)               :: k2Tlab
!         real(NER), intent(in)   :: k

!         k2Tlab = (2.0 * k * k + PC_mN * PC_mN) / PC_mN - PC_mN

!     end function k2Tlab


!     function k2Ecm(k)

!         real(NER)               :: k2Ecm
!         real(NER), intent(in)   :: k
!         k2Ecm = 2.0 * (sqrt(k*k + PC_mN * PC_mN) - PC_mN)

!     end function k2Ecm


!     function k2Ecm_nonrel(k)

!         real(NER)               :: k2Ecm_nonrel
!         real(NER), intent(in)   :: k

!         k2Ecm_nonrel = k*k/PC_mN

!     end function k2Ecm_nonrel

!     ! Compute TLab from kcm (nonrelativistic)
!     function Tlab2k(Tl)

!         real(NER)               :: Tlab2k
!         real(NER), intent(in)   :: Tl

!         Tlab2k = sqrt(PC_mN * Tl * 0.5)

!     end function Tlab2k


end module util_gadgets
