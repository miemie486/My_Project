! Bingwei Long 11/05/2013
! Calculate the phase shifts of peripherial waves where OPE is perturbative.

program get_ope_ladder

    use drv_pwphase
    implicit none

    integer, parameter          :: N = 100, regtype = REGTYPE_GAUSSIAN
    integer                     :: chnid, uptoQn = 0
    real(NER)                   :: ratio_Lambda_to_Mambda = 1.5, time_start, time_end
    class(obj_channel), pointer :: ptr_chn => null()
    character(len = 3)          :: chnstr
    character(len = 64)         :: cmdbuff
    logical                     :: succeeded

    call get_command_argument(1, cmdbuff)
    read (cmdbuff, '(a3)') chnstr
    call get_command_argument(2, cmdbuff)
    read (cmdbuff, '(i2)') uptoQn

    call cpu_time(time_start)

    call convert_text_to_chnid(chnstr, chnid, succeeded)
    if (type_of_channel_PP(chnid) == CHTYPE_SNGL) then
        ! debug
!         allocate(obj_woc0_taylor_sngl::ptr_chn)
        allocate(obj_pope_sngl::ptr_chn)
    else
        allocate(obj_pope_cpld::ptr_chn)
    end if
    call create_channel(ptr_chn, chnid, uptoQn)
    call get_output_for_Lambda_list(ptr_chn, N, regtype, ratio_Lambda_to_Mambda, OUTOPT_PHASE, .true.)
    call ptr_chn%erase()
    deallocate(ptr_chn)

    call cpu_time(time_end)
    print '(a, f9.2, a)', 'running time = ', time_end - time_start, ' seconds.'

end program get_ope_ladder
