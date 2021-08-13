! Bingwei Long
! Class: dedicated to study behaviors of residual counterterms

module mod_obj_withc0_residual_sngl

    use mod_obj_withc0_sngl
    implicit none

    type, extends(obj_withc0_sngl)    :: obj_withc0_residual_sngl

    contains

        procedure   :: init                     => init_withc0_residual_sngl
        procedure   :: get_num_paras            => get_num_paras_withc0_residual_sngl
        procedure   :: get_paras_filename       => get_paras_filename_withc0_residual_sngl
        procedure   :: get_allTQn               => get_allTQn_withc0_residual_sngl
        procedure   :: get_C1D0_by_phase        => get_C1D0_by_phase_withc0_residual_sngl
        procedure   :: output_phase_to_string   => output_phase_to_string_withc0_residual_sngl
        procedure   :: get_output_formatted_text_for_klist  => get_output_formatted_text_for_klist_withc0_residual_sngl

    end type obj_withc0_residual_sngl


contains

    subroutine init_withc0_residual_sngl(self, N, regtype, Mambda)

        class(obj_withc0_residual_sngl), intent(inout)  :: self
        integer, intent(in)                             :: N, regtype
        real(NER), intent(in)                           :: Mambda

        self%uptoQn = 0
        call init_withc0_sngl(self, N, regtype, Mambda)
        call self%create_Tk(self%Tk)
        call self%create_vstruct(self%vs_vl0)
        call self%init_symmtr_vstruct(self%vs_vl0, Vfunc_OPE_sngl)
        call self%create_VM(self%VMVLO)
        call self%create_VM_sngl(self%VM_poly0)
        call self%init_VMpoly0(self%VM_poly0)
        if (.not. self%inpar_set) &
            write (standard_error_unit, '(a)') 'init_withc0_residual_sngl: Warning! C0 has yet been defined.'
        call self%get_C1D0_by_phase(PHASE_INCREMENTAL)

    end subroutine init_withc0_residual_sngl


    subroutine get_num_paras_withc0_residual_sngl(self, num_cnttrm)

        class(obj_withc0_residual_sngl), intent(inout)  :: self
        integer, intent(out)                            :: num_cnttrm

        num_cnttrm = 1

    end subroutine get_num_paras_withc0_residual_sngl


    subroutine get_paras_filename_withc0_residual_sngl(self, fname)

        class(obj_withc0_residual_sngl), intent(inout)  :: self
        character(len = *), intent(out)                 :: fname

        fname = 'run_'//self%chnstr//'_Cs_Q0.in'

    end subroutine get_paras_filename_withc0_residual_sngl

    ! inparas(2) = kcm_1
    ! inparas(3) = phase_1
    ! inparas(4) = kcm_2
    ! inparas(5) = phase_2

    subroutine get_C1D0_by_phase_withc0_residual_sngl(self, phase_option)

        class(obj_withc0_residual_sngl), intent(inout)  :: self
        integer, intent(in)                             :: phase_option

        integer         :: tmp_uptoQn
        ! , uflag(1:MAX_NUMBER_ORDER)
        real(NER)       :: a(1:2, 1:2), b(1:2)
        complex(NEC)    :: CT, DT

        call self%update_k(self%inparas(2))
        tmp_uptoQn = self%uptoQn
        self%uptoQn = 0
        call self%get_allTQn()
        if (phase_option == PHASE_INCREMENTAL) then
            b(1) = self%inparas(3)
        else
            call self%convert_TQn_to_phase_quiet()
            b(1) = real(delta2t_nlo(self%inparas(3)-self%phsn(1), self%phsn(1), self%k))
        end if
        call self%get_CDT(self%Tk, CT, DT)
        a(1, 1) = real(CT)
        a(1, 2) = real(DT)

        call self%update_k(self%inparas(4))
        call self%get_allTQn()
        if (phase_option == PHASE_INCREMENTAL) then
            b(2) = self%inparas(5)
        else
            call self%convert_TQn_to_phase_quiet()
            b(1) = real(delta2t_nlo(self%inparas(3)-self%phsn(1), self%phsn(1), self%k))
        end if

        call self%convert_TQn_to_phase_quiet()
        b(2) = real(delta2t_nlo(self%inparas(5)-self%phsn(1), self%phsn(1), self%k))
        call self%get_CDT(self%Tk, CT, DT)
        a(2, 1) = real(CT)
        a(2, 2) = real(DT)

        ! C1 = Cs(2), D0 = Cs(3)
        self%Cs(3) = (b(2)*a(1, 1) - b(1)*a(2, 1))/(a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))
        self%Cs(2) = (b(1)-a(1, 2)*self%Cs(3))/a(1, 1)

        self%uptoQn = tmp_uptoQn

    end subroutine get_C1D0_by_phase_withc0_residual_sngl


    subroutine get_allTQn_withc0_residual_sngl(self)

        class(obj_withc0_residual_sngl), intent(inout)   :: self

        integer         :: N
        ! , tmp_uptoQn
        complex(NEC)    :: CT, DT

        ! call self%obj_sngl%get_allTQn()
        call get_allTQn_channel(self)
        N = self%N
        if (.not. self%VMVLO%inited) call self%init_VMVLO(self%Cs(1))
        call self%get_Tk_from_VM(self%VMVLO, self%Tk)
        self%TQn(1) = self%Tk%ptr(N+1)
        if (self%uptoQn == 0) then
            self%TQn_updated = .true.
            return
        end if

        call self%get_CDT(self%Tk, CT, DT)
        ! C1 = Cs(2), D0 = Cs(3)
        self%TQn(2) = self%Cs(2)*CT + self%Cs(3)*DT
        ! self%uptoQn = tmp_uptoQn
        self%TQn_updated = .true.

    end subroutine get_allTQn_withc0_residual_sngl


    ! Output the phases in the incremental form

    subroutine output_phase_to_string_withc0_residual_sngl(self, outstr)

        class(obj_withc0_residual_sngl), intent(inout)  :: self
        character(len = *), intent(out)                 :: outstr

        ! call self%obj_channel%output_phase_to_string(outstr)
        call output_phase_to_string_channel(self, outstr)
        call write_phsn_to_text(self%uptoQn+1, self%phsn, outstr)

    end subroutine output_phase_to_string_withc0_residual_sngl


    subroutine get_output_formatted_text_for_klist_withc0_residual_sngl(self, N, regtype, Mambda, kk, numk, inputs, num_inputs, output_option, output_text)

        class(obj_withc0_residual_sngl), intent(inout)  :: self
        integer, intent(in)                             :: regtype, N, numk, num_inputs
        real(NER), intent(in)                           :: Mambda
        real(NER), dimension(1:), intent(in)            :: kk, inputs
        integer, intent(in)                             :: output_option
        character(len=*), dimension(1:), intent(inout)  :: output_text

        integer     :: ii
        real(NER)   :: B0, err_B0, B2, err_B2, kcot

        if (output_option == OUTOPT_PHASE) then
            ! call self%obj_withc0_sngl%get_output_formatted_text_for_klist(N, regtype, Mambda, Lambda, kk, numk, inputs, num_inputs, output_option, output_text)
            call get_output_formatted_text_for_klist_channel(self, N, regtype, Mambda, kk, numk, inputs, num_inputs, output_option, output_text)
            return
        end if
        ! Frozen
        ! if (output_option == OUTOPT_B0B2) then
        !     call self%load_inputs(num_inputs, inputs)
        !     call self%init(N, regtype, Mambda, Lambda)
        !     do ii = 1, numk
        !         call self%get_B0(kk(ii), B0, err_B0)
        !         call self%get_B2(kk(ii), B2, err_B2)
        !         ! Frozen
        !         ! call self%get_kcot(kk(ii), kcot)
        !         write (output_text(ii), '(f9.3, 2x, f9.3, 2x, f9.3, 2x, f9.3, 2x, f9.3, 2x, f9.3)') kcot, atan(kk(ii)/kcot)*180.0_NER/PI_NE, B0, err_B0, B2, err_B2
        !     end do
        !     call self%release()
        ! end if

    end subroutine get_output_formatted_text_for_klist_withc0_residual_sngl


end module mod_obj_withc0_residual_sngl
