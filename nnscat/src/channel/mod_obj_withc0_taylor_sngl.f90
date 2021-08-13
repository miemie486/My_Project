! Bingwei Long 10/23/2013
! Class obj_withc0_taylor_sngl (Single channels where OPE is nonperturbative and
! C0 is at LO.)
! Still under construction 01/14/2014

! Abandoned, 07/14/2018
module mod_obj_withc0_taylor_sngl

    use mod_obj_taylor_sngl
    implicit none

    type, extends(obj_woc0_taylor_sngl) :: obj_withc0_taylor_sngl

        type(vstruct_sngl)  :: vs_ope, vs_2pi0, vs_dummy

    contains

        procedure   :: init                 => init_withc0_taylor_sngl
        procedure   :: release              => release_withc0_taylor_sngl
        procedure   :: get_paras_filename  => get_paras_filename_withc0_taylor_sngl
        procedure   :: update_k             => update_k_withc0_taylor_sngl
!         procedure   :: get_num_paras       => get_num_paras_withc0_taylor_sngl
        procedure   :: get_allTQn           => get_allTQn_withc0_taylor_sngl

    end type obj_withc0_taylor_sngl


contains


    subroutine init_withc0_taylor_sngl(self, N, regtype, Mambda, Lambda)

        class(obj_withc0_taylor_sngl), intent(inout)  :: self
        integer, intent(in)                         :: N, regtype
        real(NER), intent(in)                       :: Mambda, Lambda

        ! call self%obj_taylor_sngl%init(N, regtype, Mambda, Lambda)
        call init_taylor_sngl(self, N, regtype, Mambda, Lambda)
        ! debug
        self%uptoQn = 2
        call self%create_vstruct(self%vs_ope)
        call self%create_vstruct(self%vs_2pi0)
        call self%create_vstruct(self%vs_dummy)
        call self%init_symmtr_vstruct(self%vs_ope, Vfunc_OPE_sngl)
        ! debug
        call self%init_symmtr_vstruct(self%vs_2pi0, Vfunc_OPE_sngl)

    end subroutine init_withc0_taylor_sngl


    subroutine release_withc0_taylor_sngl(self)

        class(obj_withc0_taylor_sngl), intent(inout) :: self

        call self%release_vstruct(self%vs_ope)
        call self%release_vstruct(self%vs_2pi0)
        call self%release_vstruct(self%vs_dummy)
        ! call self%obj_taylor_sngl%release()
        call release_taylor_sngl(self)

    end subroutine release_withc0_taylor_sngl


    subroutine get_paras_filename_withc0_taylor_sngl(self, fname)

        class(obj_withc0_taylor_sngl), intent(inout)   :: self
        character(len = *), intent(out)     :: fname

        fname = 'run_'//self%chnstr//'_Lambda.in'

    end subroutine get_paras_filename_withc0_taylor_sngl


    subroutine update_k_withc0_taylor_sngl(self, k)

        class(obj_withc0_taylor_sngl), intent(inout) :: self
        real(NER), intent(in)               :: k

        ! call self%obj_taylor_sngl%update_k(k)
        call update_k_sngl(self, k)
        call self%update_symmtr_edge_vstruct(self%vs_ope)
        call self%update_symmtr_edge_vstruct(self%vs_2pi0)
        self%Tk%updated = .false.

    end subroutine update_k_withc0_taylor_sngl


    subroutine get_allTQn_withc0_taylor_sngl(self)

        class(obj_withc0_taylor_sngl), intent(inout) :: self

        complex(NEC), dimension(1:10, 1:10)    :: tmp_amp

        ! call self%obj_sngl%get_allTQn()
        call get_allTQn_sngl(self)
        ! debug
        call self%get_n1_ins_vs1_n2_ins_vs2_default(self%vs_dummy%VM, 1, self%vs_ope%VM, 1, self%vs_2pi0%VM, tmp_amp)
        self%TQn(1) = tmp_amp(1, 1)
        self%TQn(2) = tmp_amp(1, 2)
        ! debug. the following is true only when vs_ope is identical to vs_2pi0
        self%TQn(3) = tmp_amp(2, 2)*0.5_NER
        self%TQn_updated = .true.

    end subroutine get_allTQn_withc0_taylor_sngl


end module mod_obj_withc0_taylor_sngl
