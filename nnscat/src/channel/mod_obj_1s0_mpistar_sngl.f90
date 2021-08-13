! Bingwei Long 05/26/2014
! Class: 1S0 with mpistar counting, ie, expansion around the unitarity limit


module mod_obj_1s0_mpistar_sngl

    use mod_obj_1s0LY
    implicit none

    type, extends(obj_1s0LY)  :: obj_1s0_mpistar_sngl

        type(TM_sngl)   :: Tpq_yukawa
        real(NER)       :: mpisqr

    contains

        procedure   :: init             => init_1s0_mpistar_sngl
        procedure   :: release          => release_1s0_mpistar_sngl
        procedure   :: load_inputs      => load_inputs_1s0_mpistar_sngl
        procedure   :: get_num_paras    => get_num_paras_1s0_mpistar_sngl
        procedure   :: get_C0star       => get_C0star_1s0_mpistar_sngl

    end type obj_1s0_mpistar_sngl

contains


    subroutine Vfunc_1s0_mpi_sngl(self, p1, p2, potval)

        class(obj_sngl), intent(inout)    :: self
        real(NER), intent(in)          :: p1, p2
        real(NER), intent(out)         :: potval

        select type (self)
            class is (obj_1s0_mpistar_sngl)
                potval = yukawa_1s0_mpisqr(self%mpisqr, p1, p2)
            class default
                potval = yukawa_1s0_mpisqr(PC_mpi*PC_mpi, p1, p2)
        end select
        if (self%regtype .ne. REGTYPE_NONE .and. self%regtype .ne. REGTYPE_SHARP)   &
            potval = potval * self%regulator(p1, p2)

    end subroutine Vfunc_1s0_mpi_sngl


    subroutine init_1s0_mpistar_sngl(self, N, regtype, Mambda)

        class(obj_1s0_mpistar_sngl), intent(inout)  :: self
        integer, intent(in)                         :: N, regtype
        real(NER), intent(in)                       :: Mambda

        call init_sngl(self, N, regtype, Mambda)
        call self%create_Tk(self%Tk)
        call self%create_vstruct(self%vs_vl0)
        call self%init_symmtr_vstruct(self%vs_vl0, Vfunc_1s0_mpi_sngl)
        call self%create_VM(self%VMVLO)
        call self%create_VM_sngl(self%VM_poly0)
        call self%init_VMpoly0(self%VM_poly0)
        call self%create_TM(self%Tpq_yukawa)
        call self%get_C0star(self%Cs(1))

    end subroutine init_1s0_mpistar_sngl


    subroutine release_1s0_mpistar_sngl(self)

        class(obj_1s0_mpistar_sngl), intent(inout)  :: self

        if (self%Tpq_yukawa%created) call self%release_TM(self%Tpq_yukawa)
        ! call self%obj_1s0LY%release()
        call release_1s0LY(self)

    end subroutine release_1s0_mpistar_sngl


    subroutine load_inputs_1s0_mpistar_sngl(self, num_inputs, inputs)

        class(obj_1s0_mpistar_sngl), intent(inout)      :: self
        integer, intent(in)                             :: num_inputs
        real(NER), dimension(:), intent(in)             :: inputs

        self%inparas(1:num_inputs+1) = inputs(1:num_inputs+1)
        self%mpisqr = self%inparas(1)**2

    end subroutine load_inputs_1s0_mpistar_sngl


    subroutine get_num_paras_1s0_mpistar_sngl(self, num_cnttrm)

        class(obj_1s0_mpistar_sngl), intent(inout)  :: self
        integer, intent(out)                :: num_cnttrm

        select case (self%uptoQn)
            case (0)
                ! mpistar
                num_cnttrm = 1
!            case (1)
!                ! C0, C1, D0
!                num_cnttrm = 3
!            case (2)
!                ! C0, C1, D0, C2, D1, E0
!                num_cnttrm = 6
            case default
                write (standard_error_unit, '(a, i2)')  &
                    'get_num_paras_1s0LY: I don''t know how to handle order ', self%uptoQn
                num_cnttrm = 0
        end select

    end subroutine get_num_paras_1s0_mpistar_sngl


    ! C0star is the value that generates a pole at k = 0. Assuming that the pure
    ! Yukawa doesn't generate a pole
    subroutine get_C0star_1s0_mpistar_sngl(self, C0star)

        class(obj_1s0_mpistar_sngl), intent(inout)  :: self
        real(NER), intent(out)                      :: C0star

        complex(NER)            :: Ik
        real(NER), parameter    :: local_EPS = 1.0E-3_NER

        call self%update_k(0.001_NER)
        call self%get_Tpq_from_VM(self%vs_vl0%VM, self%Tpq_yukawa)
        call self%get_loop_Ik(1, 1, self%Tpq_yukawa, Ik)
        if (abs(aimag(Ik)/real(Ik)) > local_EPS)  &
            write (standard_error_unit, '(a)') 'get_C0star_1s0_mpistar_sngl: Warning! Ik not real'
        C0star = 1.0_NER/real(Ik)
        print *, C0star

    end subroutine get_C0star_1s0_mpistar_sngl


end module mod_obj_1s0_mpistar_sngl
