! Bingwei Long 07/30/2015
! Non-typical

program testrun

    use drv_pwphase
    implicit none

    integer, parameter          :: N = 100, regtype = REGTYPE_GAUSSIAN
    integer                     :: chnid, uptoQn = 0
    real(NER)                   :: ratio_Lambda_to_Mambda = 2.0
    class(obj_channel), pointer :: ptr_chn => null()
    character(len = 3)          :: chnstr
    character(len = 64)         :: cmdbuff
    logical                     :: fxk_is_needed = .false., succeeded
    type(lsj_symbol)            :: lsj

    call get_command_argument(1, cmdbuff)
    read (cmdbuff, '(a3)') chnstr
    call get_command_argument(2, cmdbuff)
    read (cmdbuff, '(i1)') uptoQn
    call get_command_argument(3, cmdbuff)
    select case (cmdbuff)
        case ('fxk')
            fxk_is_needed = .true.
        case default
    end select

    ! call allocate_create_channel_class(ptr_chn, command_str, uptoQn)
    !!! debug
    call convert_text_to_lsj(chnstr, lsj, succeeded)
    call convert_lsj_to_chnid(lsj, chnid)
    select case (chnid)
        case (CHNID_3P0_PP)
            allocate(obj_withc0_sngl::ptr_chn)
        case (CHNID_3P1_PP)
            allocate(obj_withc0_sngl::ptr_chn)
        case (CHNID_1S0_PP)
            allocate(obj_1s0LY::ptr_chn)
        case (CHNID_3S1_PP)
            allocate(obj_cpld::ptr_chn)
        case default
            write (standard_error_unit, '(a, a3)')  &
                'testrun: Don''t know how to handle', chnstr
            stop
    end select
    call create_channel(ptr_chn, chnid, uptoQn)

    call get_output_for_Lambda_list(ptr_chn, N, regtype, ratio_Lambda_to_Mambda, OUTOPT_PHASE, fxk_is_needed)
    call ptr_chn%erase()
    deallocate(ptr_chn)

end program testrun
