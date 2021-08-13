! Bingwei Long

module mod_obj_dibcpld

    use mod_obj_cpld
    implicit none

    ! Coupled channels where a dibaryon couples to the lower partial wave
    type, extends(obj_cpld) :: obj_dibcpld

    contains

        ! procedure   :: get_num_paras        => get_num_paras_dibcpld
        procedure   :: set_EdepCs           => set_EdepCs_dibcpld
        procedure   :: update_k             => update_k_dibcpld
        procedure   :: get_allTQn           => get_allTQn_dibcpld

    end type obj_dibcpld

contains

    ! Frozen
    ! subroutine get_num_paras_dibcpld(self, num_cnttrm)

    !     class(obj_dibcpld), intent(inout)   :: self
    !     integer, intent(out)                :: num_cnttrm

    !     select case (self%uptoQn)
    !         case (0)
    !             num_cnttrm = 2
    !         ! case (1)
    !         !     num_cnttrm = 5
    !         ! case (2)
    !         !     num_cnttrm = 9
    !         case default
    !             write (standard_error_unit, '(a, i2)')  &
    !                 'get_num_paras_dibcpld: I don''t know how to handle order ', self%uptoQn
    !             num_cnttrm = 0
    !     end select

    ! end subroutine get_num_paras_dibcpld

    subroutine update_k_dibcpld(self, k)

        class(obj_dibcpld), intent(inout)   :: self
        real(NER), intent(in)               :: k

        call self%set_EdepCs(k)
        self%VMVLO%inited = .false.
        call update_k_cpld(self, k)

    end subroutine update_k_dibcpld

    subroutine set_EdepCs_dibcpld(self, k)

        class(obj_dibcpld), intent(inout)   :: self
        real(NER), intent(in)               :: k

        real(NER)   :: Ecm

        ! Delta = inparas(1), y = inparas(2)
        ! Cs(1) = y/(Ecm + Delta)
        Ecm = k2Ecm(k)
        self%Cs(1) = self%inparas(2)/(Ecm + self%inparas(1))
        if (self%uptoQn == 0) return

        self%Cs(2) = self%inparas(3)/(Ecm + self%inparas(1))**2 + self%inparas(4)/(Ecm + self%inparas(1)) + self%inparas(5)

    end subroutine set_EdepCs_dibcpld

    subroutine get_allTQn_dibcpld(self)

        class(obj_dibcpld), intent(inout)   :: self

        integer         :: N
        real(NER)       :: VL(2, 2)
        complex(NEC)    :: L1(2, 2), L2(2, 2), CT(2, 2), DT(2, 2), ET(2, 2)

        call get_allTQn_channel(self)
        if (.not. self%VMVLO%inited) call self%init_VMVLO(self%Cs(1))
        call self%get_Tk_from_VM(self%VMVLO, self%Tk)
        N = self%N
        self%TQn(1:2, 1:2, 1) = self%Tk%ptr(1:2, 1:2, N+1)
        if (self%uptoQn == 0) then
            self%TQn_updated = .true.
            return
        end if

        call self%get_CDET(self%Tk, CT, DT, ET)
        self%TQn(1:2, 1:2, 3) = self%Cs(2)*CT
        if (self%uptoQn == 1) then
            self%TQn_updated = .true.
            return
        end if

        ! VL(1:2, 1:2) = self%vs_vl2%VM%ptr(1:2, 1:2, N+1, N+1)
        ! call self%get_1loop_os_VMTk(self%vs_vl2%VM, self%Tk, L1)
        ! call self%get_2loops_TkVMTk(self%vs_vl2%VM, self%Tk, L2)
        ! call self%get_CDET(self%Tk, CT, DT, ET)
        ! self%TQn(1:2, 1:2, 3) = VL + L1 + transpose(L1) + L2 + self%Cs(2)*CT + self%Cs(3)*DT  &
        !     + self%Cs(4)*ET
        ! if (self%uptoQn == 2) then
        !     self%TQn_updated = .true.
        !     return
        ! end if

    end subroutine get_allTQn_dibcpld

end module mod_obj_dibcpld
