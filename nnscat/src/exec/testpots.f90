! Bingwei 03/08/2011
! test born term of a potential

program testpots

    use eft_potspwd
    use ope_pwd
    use vtpe0_pwd
    use vtpe1_pwd
    use util_mesh
    implicit none

    real(NER), parameter :: nfac = PC_mN, Lambda = 1000.0, localEPS = 1.E-2
    integer, parameter :: N = 100
    real(NER) :: v, k, k0 = 0.001, dk  = 0.001, ppi, ppj, eps, t1, t2, v1, v2, v3
    real(NER), dimension(1:N) :: msh, wghts
    real(NER), dimension(1:2, 1:2) :: vab
    integer :: ii, jj, numk, argnum
    character (len=10)  :: cmdbuff
    type(flags_struct)  :: flagstrc
    logical             :: opt_valid

    call set_default_flags(flagstrc)
    dfltN = flagstrc%mesh_N
    pots_var_PP = POTS_WPC

    argnum = command_argument_count()
    do ii = 1, argnum
        call get_command_argument(ii, cmdbuff)
        call word_to_flag(trim(cmdbuff), flagstrc, opt_valid)
        if (.not. opt_valid) then
            select case (cmdbuff)
                case ('--help', '-h')
                    call output_usage()
                    stop
                case default
                    if (cmdbuff .ne. '') print '(a, a)', trim(cmdbuff), ' is an unknown parameter.'
            end select
        end if
    end do
    if (flagstrc%chnid .eq. CHNID_UNDEFINED) then
        print '(a)', 'un-identifiable channel'
        stop
    end if
    call convert_chnid_to_lsj(flagstrc%chnid, lsj, lgc_succ)
    call convert_lsj_to_text(lsj, chnstr)
    if (chnstr .eq. '') then
        print '(a)', 'no known channel is specified.'
        stop
    end if
    pots_var_PP = flagstrc%lrpots
    if (flagstrc%mesh_N .lt. 5) then
        print ('(a)'), 'N must be at least 5'
        stop
    end if
    misc%is_dibaryon = flagstrc%is_dibaryon

    numk = 30
    call set_mesh_gnr(N, Lambda, msh, wghts)
!~     call cpu_time(t1)

    do ii = 1, N, 50
        ppi = msh(ii)
        do jj = 1, N, 50
            ppj = msh(jj)
            v1 = OPE_1S0_lbw(ppi, ppj)
            print *, 'OPE', ppi, ppj, v1
            v1 = VLNLO_1S0_lbw(ppi, ppj)
            print *, 'VLNLO', ppi, ppj, v1
            v1 = VLNNLO_1S0_lbw(ppi, ppj)
            print *, 'VLNNLO', ppi, ppj, v1
       end do
   end do


!~    print *, t2-t1, '  seconds'

end program testpots
