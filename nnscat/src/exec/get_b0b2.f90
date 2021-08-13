! Bingwei Long 04/17/2014
! Calculate B0 and B2 in channels where OPE is nonperturbative

program get_b0b2

    use drv_pwphase
    implicit none

    integer, parameter          :: N = 100, regtype = REGTYPE_GAUSSIAN
    integer                     :: chnid, uptoQn = 0
    real(NER)                   :: ratio_Lambda_to_Mambda = 1.5
    class(obj_channel), pointer :: ptr_chn => null()
    character(len = 3)          :: chnstr
    character(len = 64)         :: cmdbuff
    logical                     :: succeeded

    call get_command_argument(1, cmdbuff)
    read (cmdbuff, '(a3)') chnstr
    uptoQn = 0
    ! call get_command_argument(2, cmdbuff)
    ! read (cmdbuff, '(i1)') uptoQn

    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (.not. succeeded) then
        write (standard_error_unit, '(a, a)') 'Error in converting ', chnstr
        stop
    end if

    allocate(obj_withc0_residual_sngl::ptr_chn)
    call create_channel(ptr_chn, chnid, uptoQn)
    call get_output_for_Lambda_list(ptr_chn, N, regtype, ratio_Lambda_to_Mambda, OUTOPT_B0B2, .false.)
    call ptr_chn%erase()
    deallocate(ptr_chn)

end program get_b0b2
