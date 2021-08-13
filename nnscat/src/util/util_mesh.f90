! Bingwei Long 04/09/2012
! utility procedures for mesh

! Obsolete: 08/10/2018
module util_mesh

	use nneft_type
	use util_gauleg, only : feb_glabsc_2pt, feb_tanabsc
	implicit none

    ! in the case of Gaussian regulator, the ratio of the mesh upper limit to Mambda
    real(NER), parameter :: PC_ratio_Lambda_Mambda = 1.5
    real(NER), parameter, private :: large_ratio_max = 1000.0

	contains


    ! set up external msh and wghts *regardless* of k

    subroutine set_mesh_gnr(N, start, Lambda, msh, wghts)

        integer, intent(in) 		:: N
        real(NER), intent(in)       :: start, Lambda
        real(NER), intent(inout)    :: msh(:), wghts(:)

        call  feb_glabsc_2pt(N, msh, wghts, start, Lambda)

    end subroutine set_mesh_gnr


    ! fill in Vk(i) with ppi^nu

    subroutine set_Vk_power_sngl(N, msh, k, nu, Vk)

        integer, intent(in)     :: N, nu
        ! real(NER), dimension(1:N), intent(in)       :: msh
        ! real(NER), dimension(1:N+1), intent(out)    :: Vk
        real(NER), intent(in)   :: k, msh(:)
        real(NER), intent(out)  :: Vk(:)

        integer     :: ii
        real(NER)   :: ppi

        forall (ii = 1:N) Vk(ii) = msh(ii)**nu
        Vk(N+1) = k**nu

    end subroutine set_Vk_power_sngl


    ! fill in outside Vk(i) with Vfnc(p, k)

    subroutine set_Vk_gnr_sngl(N, msh, k, Vfnc, Vk)

        integer, intent(in)                         :: N
        ! real(NER), intent(in)                       :: k
        ! real(NER), dimension(1:N), intent(in)       :: msh
        ! real(NER), dimension(1:N+1), intent(out)    :: Vk
        real(NER), intent(in)   :: k, msh(:)
        real(NER), intent(out)  :: Vk(:)

        interface
            function Vfnc(p1, p2)
                use nneft_type
                implicit none
                real(NER)              :: Vfnc
                real(NER), intent(in)  :: p1, p2
            end function Vfnc
        end interface

        integer     :: ii
        real(NER)   :: ppi

        do ii = 1, N
            ppi = msh(ii)
            Vk(ii) = Vfnc(ppi, k)
        end do
        Vk(N+1) = Vfnc(k, k)

    end subroutine set_Vk_gnr_sngl


    subroutine set_VM_gnr_sngl(N, msh, Vfnc, VM)

        integer, intent(in)         :: N
        real(NER), intent(in)       :: msh(:)
        real(NER), intent(inout)    :: VM(:, :)

        interface
            function Vfnc(p1, p2)
                use nneft_type
                implicit none
                real(NER)                :: Vfnc
                real(NER), intent(in)    :: p1, p2
            end function Vfnc
        end interface

        integer     :: ii, jj
        real(NER)  :: ppi, ppj

        do ii = 1, N
            ppi = msh(ii)
            do jj = 1, ii - 1
                ppj = msh(jj)
                VM(ii, jj) = Vfnc(ppi, ppj)
                VM(jj, ii) = VM(ii, jj)
            end do
            VM(ii, ii) = Vfnc(ppi, ppi)                     ! set up the diagonal elements
        end do

    end subroutine set_VM_gnr_sngl

! build the verge of VM, VM(N_AS+1, j) = VM(j, N_AS+1) = Vfnc(ppj, k) = Vfnc(k, ppj)
! Vfnc is assumed to be symmtric

    subroutine set_VM_edge_gnr_sngl(N, msh, k, Vfnc, VM)

        integer, intent(in)         :: N
        real(NER), intent(in)       :: k
        ! real(NER), dimension(1:N), intent(in)                   :: msh
        real(NER), intent(in)       :: msh(:)
        real(NER), intent(inout)    :: VM(:, :)

        interface
            function Vfnc(p1, p2)
                use nneft_type
                implicit none
                real(NER)                :: Vfnc
                real(NER), intent(in)    :: p1, p2
            end function Vfnc
        end interface

        integer     :: jj
        real(NER)   :: ppj

        do jj = 1, N
            ppj = msh(jj)
            VM(N+1, jj) = Vfnc(ppj, k)
            VM(jj, N+1) = VM(N+1, jj)
        end do
        VM(N+1, N+1) = Vfnc(k, k)

    end subroutine set_VM_edge_gnr_sngl


    subroutine set_VM_gnr_cpld(N, msh, Vfnc, VM)

        integer, intent(in)         :: N
        ! real(NER), dimension(1:N), intent(in)                   :: msh
        ! real(NER), dimension(1:N+1, 1:N+1, 2, 2), intent(inout) :: VM
        real(NER), intent(in)       :: msh(:)
        real(NER), intent(inout)    :: VM(:, :, :, :)


        interface
            subroutine Vfnc(p1, p2, vab)
                use nneft_type
                implicit none
                real(NER), intent(in)                   :: p1, p2
                real(NER), dimension(2, 2), intent(out) :: vab
            end subroutine Vfnc
        end interface

        integer                     :: ii, jj
        real(NER)                   :: ppi, ppj
        real(NER), dimension(2, 2)  :: vab

        do ii = 1, N
            ppi = msh(ii)
            do jj = 1, ii - 1
                ppj = msh(jj)
                call Vfnc(ppi, ppj, vab)
                VM(1:2, 1:2, ii, jj) = vab
                VM(1:2, 1:2, jj, ii) = transpose(vab)
            end do
            ! set up the diagonal elements
            call Vfnc(ppi, ppi, vab)
            VM(1:2, 1:2, ii, ii) = vab
        end do

    end subroutine set_VM_gnr_cpld


    ! build the verge of VM.

    subroutine set_VM_edge_gnr_cpld(N, msh, k, Vfnc, VM)

        integer, intent(in)                                     :: N
        real(NER), intent(in)                                   :: k
        ! real(NER), dimension(1:N), intent(in)                   :: msh
        real(NER), intent(in)       :: msh(:)
        ! real(NER), dimension(1:N+1, 1:N+1, 2, 2), intent(inout) :: VM
        real(NER), intent(inout)    :: VM(:, :, :, :)

        interface
            subroutine Vfnc(p1, p2, vab)
                use nneft_type
                implicit none
                real(NER), intent(in)                   :: p1, p2
                real(NER), dimension(2, 2), intent(out) :: vab
            end subroutine Vfnc
        end interface

        integer                     :: jj
        real(NER)                   :: ppj
        real(NER), dimension(2, 2)  :: vab

        do jj = 1, N
            ppj = msh(jj)
            call Vfnc(k, ppj, vab)
            VM(1:2, 1:2, N+1, jj) = vab
            VM(1:2, 1:2, jj, N+1) = transpose(vab)
        end do
        call Vfnc(k, k, vab)
        VM(1:2, 1:2, N+1, N+1) = vab

    end subroutine set_VM_edge_gnr_cpld


    subroutine set_Vk_gnr_cpld(N, msh, k, Vfnc, Vk)

        integer, intent(in)                             :: N
        real(NER), intent(in)                           :: k
        ! real(NER), dimension(1:N), intent(in)           :: msh
        real(NER), intent(in)       :: msh(:)
        ! real(NER), dimension(1:N+1, 2, 2), intent(out)  :: Vk
        real(NER), intent(out)  :: Vk(:, :, :)

        interface
            subroutine Vfnc(p1, p2, vab)
                use nneft_type
                implicit none
                real(NER), intent(in)                   :: p1, p2
                real(NER), dimension(2, 2), intent(out) :: vab
            end subroutine Vfnc
        end interface

        integer                     :: ii
        real(NER)                   :: ppi
        real(NER), dimension(2, 2)  :: vab

        do ii = 1, N
            ppi = msh(ii)
            call Vfnc(ppi, k, vab)
            Vk(1:2, 1:2, ii) = vab
        end do
        call Vfnc(k, k, vab)
        Vk(1:2, 1:2, N+1) = vab

    end subroutine set_Vk_gnr_cpld


    subroutine check_mesh_crash(k, Lambda, N, msh, wght, crashed)

        real(NER), intent(in)                           :: k, Lambda
        integer, intent(in)                             :: N
        ! real(NER), dimension(1:N), intent(in)           :: msh, wght
        real(NER), intent(in)       :: msh(:), wght(:)
        logical, intent(out)                            :: crashed

        integer     :: ii
        real(NER)   :: ratio1, ratio2
        logical     :: bisected

        bisected = .false.
        ii = 1
        if (k .lt. msh(1) .and. k .gt. 0.0) then
            ratio1 = 0.0
            ratio2 = abs(wght(1)/(k - msh(1)))
            bisected = .true.
        end if
        do ii = 2, N
            if (k .lt. msh(ii) .and. k .gt. msh(ii-1)) then
                ratio1 = abs(wght(ii-1)/(k - msh(ii-1)))
                ratio2 = abs(wght(ii)/(k - msh(ii)))
                bisected = .true.
                exit
            end if
        end do
        if (ii .eq. N+1) then
            if (k .lt. Lambda .and. k .gt. msh(N)) then
                ratio1 = abs(wght(N)/(k - msh(N)))
                ratio2 = 0.0
                bisected = .true.
            end if
        end if
        if (bisected) then
            if (ratio1 .gt. large_ratio_max .or. ratio2 .gt. large_ratio_max) then
                crashed = .true.
                write (standard_error_unit, '(a)') '(check_mesh_crash): k crashes onto mesh points'
                if (ii .eq. 1) then
                    write (standard_error_unit, '(f9.2, 2x, f9.2, 2x, f9.2)') 0.0, k, msh(1)
                else
                    if (ii .eq. N+1) then
                        write (standard_error_unit, '(f9.2, 2x, f9.2, 2x, f9.2)') msh(N), k, Lambda
                    else
                        write (standard_error_unit, '(f9.2, 2x, f9.2, 2x, f9.2)') msh(ii-1), k, msh(ii)
                    end if
                end if
                write (standard_error_unit, '(f8.2, 2x, f8.2)') ratio1, ratio2
            else
                crashed = .false.
            end if
        else
            write (standard_error_unit, '(a)') '(check_mesh_crash): cannot bisect'
            crashed = .false.
        end if

    end subroutine check_mesh_crash


end module util_mesh
