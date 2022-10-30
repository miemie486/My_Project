! Bingwei Long 06/27/2012
! Bingwei Long 05/23/2012
! solve for the values of C's

program slvcs

    use nneft_type
    use eft_potspwd
    use obj_pwamp
    use util_mesh
    use util_optstruct
    use util_io
    use util_rootfinder
    use util_gadgets
    use eft_phaseconv
    implicit none

    ! external subroutines from MINUIT
    external MNINIT, MNEXCM, MNPARM, MNSETI, MNPOUT

    integer,    parameter   :: funt1_glb = 23, numsegs_default = 1000
    real(NER),  parameter   :: goal_default = 1.0E-5, rel_tol_default = 1.0E-7

    integer             :: eflag, ii, numLambdas, numnij, argnum, &
                           numsegs = numsegs_default, N, N_default, solver_selector
    real(NER)           :: Lmbd, time_start, time_end, goal = goal_default, rel_tol = rel_tol_default, meshupper
    character (len=128) :: cmdbuff, fname1
    character (len=3)   :: chnstr
    logical             :: opt_valid, lgc_succ, to_fit = .false.

    real(NER), dimension(1:200)         :: Lambda
    real(NER), dimension(1:50, 1:2)     :: nijphs
    real(NER), dimension(1:50, 1:4)     :: nijmxphs
    type(neobj_pwamp)                   :: pwamp
    type(opts_pwamp)                    :: misc_opts
    type(flags_struct)                  :: flagstrc
    type(lsj_symbol)                    :: lsj

    call cpu_time(time_start)

    call set_default_flags(flagstrc)
    N_default = 100

    argnum = command_argument_count()
    do ii = 1, argnum
        call get_command_argument(ii, cmdbuff)
        call word_to_flag(trim(cmdbuff), flagstrc, opt_valid)
        if (opt_valid) cycle
        select case (cmdbuff)
            case ('--fit', '-f')
                to_fit = .true.
            case default
                if (cmdbuff(1:6) .eq. '--help' .or.  cmdbuff(1:2) .eq. '-h') then
                    call output_usage()
                    stop
                end if
                if (cmdbuff(1:5) .eq. '-seg=') then
                    read(cmdbuff(6:), '(i8)') numsegs
                    if (numsegs .le. 0) numsegs = numsegs_default
                    cycle
                end if
                if (cmdbuff(1:6) .eq. '-goal=') then
                    read(cmdbuff(7:), '(e10.2e3)') goal
                    goal = abs(goal)
                    cycle
                end if
                if (cmdbuff(1:5) .eq. '-tol=') then
                    read(cmdbuff(6:), '(e10.2e3)') rel_tol
                    rel_tol = abs(rel_tol)
                    cycle
                end if
                if (cmdbuff .ne. '') then
                    print '(a, a)', trim(cmdbuff), ' is an unknown parameter.'
                    stop
                end if
        end select
    end do
    if (flagstrc%chnid .eq. CHNID_UNDEFINED) then
        print '(a)', 'un-identifiable channel'
        call output_usage()
        stop
    end if
    call convert_chnid_to_lsj(flagstrc%chnid, lsj, lgc_succ)
    call convert_lsj_to_text(lsj, chnstr)
    if (chnstr .eq. '') then
        print '(a)', 'no known channel is specified.'
        stop
    end if
    N = flagstrc%mesh_N
    if (N .lt. 5) then
        print ('(a)'), 'N must be at least 5'
        stop
    end if
    misc_opts%is_dibaryon = flagstrc%is_dibaryon
    misc_opts%is_pert     = flagstrc%is_pert

    call goto_solver()
    call cpu_time(time_end)
    print '(a, f9.2, a)', 'running time = ', time_end - time_start, ' seconds.'


    contains


    subroutine output_usage()

        print '(a)', 'slvcs CHNID [-gaussian|-g] [-sharp|-s] [-tol=tolerance] [-goal=precision] &
        &[-seg=num_of_segments] [-N=num_of_mesh_points] [--dibaryon|-d] [--nondibaryon|-n]'
        print '(a, e10.2e3, a, e10.2e3, a, i3)', 'default options: -i -g -n -seg=1000 -tol=', rel_tol_default, ' -goal=',  goal_default, ' -N=', N_default

    end subroutine output_usage


    subroutine output_regtype(funt, regtype)

        integer, intent(in) :: funt, regtype

        select case (regtype)
            case (REGTYPE_GAUSSIAN)
                write (funt, '(a)') '# Gaussian regulator'
            case (REGTYPE_SHARP)
                write (funt, '(a)') '# Sharp regulator'
        end select

    end subroutine output_regtype


    subroutine output_dontknow()
        print '(a)', 'don''t know how to handle this'
    end subroutine output_dontknow


    subroutine output_pwa_pts(funt, num, Tlab_pts, phs_pts)

        integer, intent(in)                 :: funt, num
        real(NER), dimension(:), intent(in) :: Tlab_pts, phs_pts

        integer :: ii

        do ii = 1, num
            write (funt, '(a, f6.2, a, f9.4)')  '# Tlab = ', Tlab_pts(ii), ',  phs = ', phs_pts(ii)
        end do

    end subroutine output_pwa_pts


    subroutine goto_solver()

        if (flagstrc%order .eq. 0) then
            if (to_fit) then
                call fit_c0()
            else
                call solve_c0()
            end if
        else
            if (flagstrc%is_dibaryon) then
                select case(flagstrc%order)
                    case (1)
                        call solve_Dy_Q1()
                    case (2)
                        call solve_Dy_Q1Q2()
                    case default
                        call output_dontknow()
                end select
            else
                select case (flagstrc%chnid)
                    case (CHNID_1S0_PP)
                        if (flagstrc%order .eq. 2) then
                            call solve_1s0_Q1Q2()
                            return
                        end if
                    case (CHNID_1P1_PP, CHNID_3P1_PP)
                        if (flagstrc%order .eq. 3) then
                            call solve_c0c1_Q2Q3()
                            return
                        end if
                    case (CHNID_3S1_PP)
                        if (flagstrc%order .eq. 2) then
                            call solve_cde_cpld_Q2()
                            return
                        end if
                    case (CHNID_3P2_PP)
                        if (flagstrc%is_pert .and. flagstrc%order .eq. 3) then
                            call solve_ksw_cpld_Q3()
                            return
                        else
                            if (flagstrc%order .eq. 2) then
                                call solve_cde_cpld_Q2()
                                return
                            end if
                        end if
                end select
                call output_dontknow()
            end if
        end if

    end subroutine goto_solver


    function meshul(regtype, Mambda)

        real(NER)               :: meshul
        integer, intent(in)     :: regtype
        real(NER), intent(in)   :: Mambda

        select case(regtype)
            case (REGTYPE_GAUSSIAN)
                meshul = Lmbd*PC_ratio_Lambda_Mambda
            case default
                meshul = Lmbd
        end select

    end function meshul


    subroutine read_nijphs(funt, fname, numrec)

        integer, intent(in) :: funt
        character(len = *), intent(in)  :: fname
        integer, intent(out):: numrec

        integer :: flag

        call read2cols_io(funt, fname, nijphs(:, 1), nijphs(:, 2), numrec, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

    end subroutine read_nijphs


    subroutine read_nijmxphs(funt, fname, numrec)

        integer, intent(in) :: funt
        character(len = *), intent(in)  :: fname
        integer, intent(out):: numrec

        integer :: flag

        call read_ncols_io(funt, fname, 4, nijmxphs, numrec, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

    end subroutine read_nijmxphs


    subroutine solve_c0()

        integer             :: ii, jj, bn, cnt
        real(NER)           :: k, C0, C01, C02
        character(len = 30) :: known_chns, fname
        logical             :: bracketed
        real(NER), dimension(1:200)         :: C0i
        real(NER), dimension(:), pointer    :: xb1, xb2

        known_chns = '1s0_3s1_3p0_3p2'
        if (index(known_chns, chnstr) .eq. 0) then
            print '(a, a)', 'do not know how to solve for C0 of ', chnstr
            return
        end if

        ! if there are more than one PWA points, only the first will be used
        ! for coupled channels, delta_1, the phase of the lower partial wave will be used
        fname = '1pt_'//chnstr//'_phs.in'
        call read_nijphs(funt1_glb, fname, numnij)
        k = Tlab2k(nijphs(1, 1))
        print *, 'k (MeV) = ', k, 'Phase (deg) = ', nijphs(1, 2)
        misc_opts%is_dibaryon = .false.

        fname = 'slvcs_'//chnstr//'_Lambda.in'
        ! C0i is the suggested initial value
        call read2cols_io(funt1_glb, fname, Lambda, C0i, numLambdas, eflag)
        if (eflag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if
        fname = 'slv_'//chnstr//'_C0.out'
        open (funt1_glb, file = fname, status = 'replace')
        write (funt1_glb, '(a, f6.2)') '# fitted point : Tlab = ', nijphs(1, 1)
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda     C0'

        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp_ope_VM(flagstrc%chnid, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
            call update_k_pwamp(k, pwamp)
            C01 = C0i(ii)
            C02 = -C01
            call zbrak_pointer(diff_pwa_lo, C01, C02, numsegs, xb1, xb2, bn)
            if (bn .eq. 0) then
                print *, 'couldn''t bracket at Lambda = ', Lmbd
            end if
            cnt = 0
            do jj = 1, bn
                C0 = zbrent(diff_pwa_lo, xb1(jj), xb2(jj), rel_tol*C0i(ii))
                if (abs(diff_pwa_lo(C0)/nijphs(1, 2)) .lt. goal) then
                    print *, Lmbd, C0, diff_pwa_lo(C0)
                    write (funt1_glb, *) Lmbd, C0, diff_pwa_lo(C0)
                    cnt = cnt + 1
                end if
            end do
            deallocate(xb1, xb2)
            if (cnt .gt. 1) then
                print *, 'there are more than 1 solution for Lambda = ', Lmbd
            end if
            if (cnt .eq. 0) then
                print ('(a, x, f8.2)'), 'couldn''t find any solution for Lambda =', Lmbd
            end if
            call destroy_pwamp(pwamp)
        end do

        close (funt1_glb)

    end subroutine solve_c0


    function diff_pwa_lo(C)

        real(NER)               :: diff_pwa_lo
        real(NER), intent(in)   :: C

        integer, dimension(1:1)     :: uflag
        real(NER), dimension(1:1)   :: Cn

        Cn(1) = C
        call set_Cs_pwamp(0, Cn, pwamp)
        call update_VM0_from_ope_pwamp(pwamp)
        call get_allTQn_pwamp(0, pwamp)
        call convert_TQn_to_phsn_quiet_pwamp(0, uflag, pwamp)
        diff_pwa_lo = pwamp%phsn(1) - nijphs(1, 2)

    end function diff_pwa_lo


    subroutine fit_c0()

        integer             :: ii, bn, cnt, errflag, ivarbl
        real(NER)           :: C0, C0upper, C0lower, C0err
        character(len = 30) :: known_chns, fname, prname
        logical             :: bracketed
        real(NER), dimension(1:5)   :: argmn
        real(NER), dimension(1:200) :: C0i

        known_chns = '1s0_3s1_3p0_3p2'
        if (index(known_chns, chnstr) .eq. 0) then
            print '(a, a)', 'do not know how to solve for C0 of ', chnstr
            return
        end if
        fname = 'npts_'//chnstr//'_phs.in'
        call read_nijphs(funt1_glb, fname, numnij)
        misc_opts%is_dibaryon = .false.

        fname = 'fitc0_'//chnstr//'_Lambda.in'
        ! C0i is the suggested initial value
        call read2cols_io(funt1_glb, fname, Lambda, C0i, numLambdas, eflag)
        if (eflag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if
        fname = 'fit_'//chnstr//'_C0.out'
        open (funt1_glb, file = fname, status = 'replace')
        call output_pwa_pts(funt1_glb, numnij, nijphs(1:numnij, 1), nijphs(1:numnij, 2))
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda     C0'

        call MNINIT (5, 6, 7)
        call MNSETI ('Fit C0 of '//chnstr)
        if (errflag .ne. 0) then
            print *, 'Unable to define C0 with NPARM'
            return
        end if

        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp_ope_VM(flagstrc%chnid, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
            call MNPARM (1, 'C0_'//chnstr, C0i(ii), abs(C0i(ii)*1.0E-3), 0.0, 0.0, errflag)
            argmn(1) = 2
            call MNEXCM (chisqr_c0, 'SET STRATEGY', argmn, 1, errflag, function_zero)
            ! max number of calls
            argmn(1) = 200
            ! tolerence
            argmn(2) = C0i(ii)*1.0E-8
            call MNEXCM (chisqr_c0, 'MIGRAD', argmn, 2, errflag, function_zero)
            call MNPOUT (1, prname, C0, C0err, C0upper, C0lower, ivarbl)
            if (errflag .eq. 0) then
                write (funt1_glb, *) Lmbd, C0
            end if
            call destroy_pwamp(pwamp)
        end do

        close (funt1_glb)

    end subroutine fit_c0


    subroutine chisqr_c0(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

        integer, intent(in) :: NPAR, IFLAG
        real(NER), dimension(1:NPAR), intent(in) :: XVAL
        real(NER), intent(out) :: FVAL
        real(NER), dimension(1:NPAR), intent(out) :: GRAD

        external FUTIL

        integer     :: ii
        real(NER)   :: regsum
        integer, dimension(1:1)     :: uflag
        real(NER), dimension(1:1)   :: Cn

        regsum = 0.0
        Cn(1)  = XVAL(1)
        do ii = 1, numnij
            call update_k_pwamp(Tlab2k(nijphs(ii, 1)), pwamp)
            call set_Cs_pwamp(0, Cn, pwamp)
            call update_VM0_from_ope_pwamp(pwamp)
            call get_allTQn_pwamp(0, pwamp)
            call convert_TQn_to_phsn_quiet_pwamp(0, uflag, pwamp)
            regsum = regsum + (nijphs(ii, 2) - pwamp%phsn(1))**2
        end do
        FVAL =  regsum

    end subroutine chisqr_c0


! T_nlo = CT*(D1/(Ecm+Delta)^2 + y1/(Ecm+Delta) + C0)

    subroutine solve_Dy_Q1()

        external cgesv, zgesv

        real(NER)           :: k, Ecm
        character(len = 99) :: fname
        integer             :: ii, jj, flag
        complex(NEC)        :: CT
        integer, dimension(3)           :: ipiv, uflag
        real(NER), dimension(2)         :: cnttrm_1s0
        real(NER), dimension(50)        :: Delta, y, Lambda
        complex(NEC), dimension(3)      :: BB
        complex(NEC), dimension(3, 3)   :: AA

        if (flagstrc%chnid .ne. CHNID_1S0_PP) return
        misc_opts%is_dibaryon = .true.
        fname = '3pt_1s0_phs.in'
        call read_nijphs(funt1_glb, fname, numnij)
        fname = 'slv_1s0_Dy.in'
        call read3cols_io(funt1_glb, fname, Lambda, Delta, y, numLambdas, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

        fname = 'slv_1s0_Dy_Q1.out'
        open (funt1_glb, file = fname, status = 'replace')
        write (funt1_glb, '(a, 3(f6.2, 2x))') '# fitted point : Tlab = ', nijphs(1:3, 1)
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda  Delta y  D1  y1  C0'

        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp_ope_VM(CHNID_1S0_PP, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
            cnttrm_1s0(1) = Delta(ii)
            cnttrm_1s0(2) = y(ii)
            do jj = 1, 3
                k = Tlab2k(nijphs(jj, 1))
                Ecm = k2Ecm_nonrel(k)
                call update_k_pwamp(k, pwamp)
                call set_Cs_pwamp(0, Ecm, cnttrm_1s0, pwamp)
                call update_VM0_from_ope_pwamp(pwamp)
                call get_allTQn_pwamp(0, pwamp)
                call convert_TQn_to_phsn_quiet_pwamp(0, uflag, pwamp)
                BB(jj) = delta2t_nlo(nijphs(jj, 2) - pwamp%phsn(1), pwamp%phsn(1), k)
                call get_CT_tvm_sngl(CT, pwamp%tvms)
                AA(jj, 1) = CT/(Ecm+Delta(ii))**2
                AA(jj, 2) = CT/(Ecm+Delta(ii))
                AA(jj, 3) = CT
            end do

            if (NER .eq. SP) then
                call cgesv(3, 1, AA, 3, ipiv, BB, 3, flag)
            else
                call zgesv(3, 1, AA, 3, ipiv, BB, 3, flag)
            end if
            write(funt1_glb, '(f7.2, 8(2x, e26.17e3))'), Lmbd, Delta(ii), y(ii), real(BB(1:3))
            call destroy_pwamp(pwamp)
        end do

        close(funt1_glb)

    end subroutine solve_Dy_Q1


! T_nnlo = VLR_nlo(k, k) + 2*L1 + L2 + (xi0 + VS2)*CT + D0*DT
! where VS2 = Delta1**2/y/(E+Delta)**3 + Delta2/(E+Delta)**2 + y2/(E+Delta) + C1

    subroutine solve_Dy_Q1Q2()

        external cgesv, zgesv

        real(NER)           :: k, Ecm
        character(len = 99) :: fname
        integer             :: ii, jj, flag
        complex(NEC)        :: CT, DT, xi0, L1, L2
        integer, dimension(4)           :: ipiv, uflag
        real(NER), dimension(5)         :: cnttrm_1s0
        real(NER), dimension(50)        :: Delta, y, Delta1, y1, C0, Lambda
        complex(NEC), dimension(4)      :: BB
        complex(NEC), dimension(4, 4)   :: AA

        if (flagstrc%chnid .ne. CHNID_1S0_PP) return
        misc_opts%is_dibaryon = .true.

        call solve_Dy_Q1()

        fname = '4pt_1s0_phs.in'
        call read_nijphs(funt1_glb, fname, numnij)
        fname = 'slv_1s0_Dy_Q1.out'
        call read6cols_io(funt1_glb, fname, Lambda, Delta, y, Delta1, y1, C0, numLambdas, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

        fname = 'slv_1s0_Dy_Q2.out'
        open (funt1_glb, file = fname, status = 'replace')
        write (funt1_glb, '(a, 4(f6.2, 1x))') '# fitted point : Tlab = ', nijphs(1:4, 1)
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda  Delta  y  Delta1  y1  C0  Delta2  y2  C1  D0'

        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp_ope_VM(CHNID_1S0_PP, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
            cnttrm_1s0(1) = Delta(ii)
            cnttrm_1s0(2) = y(ii)
            cnttrm_1s0(3) = Delta1(ii)
            cnttrm_1s0(4) = y1(ii)
            cnttrm_1s0(5) = C0(ii)
            call initVM12_pwamp(2, pwamp)
            do jj = 1, 4
                k = Tlab2k(nijphs(jj, 1))
                Ecm = k2Ecm_nonrel(k)
                call update_k_pwamp(k, pwamp)
                call set_Cs_pwamp(1, Ecm, cnttrm_1s0, pwamp)
!                 call set_VS_dibaryon_PP(CHNID_1S0_PP, 1, cnttrm_1s0, k)
                call update_VM0_from_ope_pwamp(pwamp)
                call get_Tpq_pwamp(pwamp)
                call get_DCT_tvm_sngl(DT, CT, pwamp%tvms)
                pwamp%TQn(2) = pwamp%tvms%cms%C1 * CT
                call convert_TQn_to_phsn_quiet_pwamp(1, uflag, pwamp)
                BB(jj) = delta2t_nnlo(nijphs(jj, 2) - pwamp%phsn(1) - pwamp%phsn(2), pwamp%phsn(2), pwamp%phsn(1), k)
                call get_C_xi0_sngl(pwamp%tvms%cms%C1, xi0, pwamp%tvms)
                BB(jj) = BB(jj) - (xi0 + Delta1(ii)**2/y(ii)/(Ecm+Delta(ii))**3)*CT

                call set_symmetric_nthVM_edge_tvm_sngl(VMID_VLNLO_tvm, VLNLO_tvm_sngl, pwamp%tvms)
                call get_L1_nthVM_tvm_sngl(VMID_VLNLO_tvm, L1, pwamp%tvms)
                call get_L2_nthVM_tvm_sngl(VMID_VLNLO_tvm, L2, pwamp%tvms)
                BB(jj) = BB(jj) - VLNLO_tvm_sngl(k, k, pwamp%tvms) - 2.0*L1 - L2

                AA(jj, 1) = CT/(Ecm+Delta(ii))**2
                AA(jj, 2) = CT/(Ecm+Delta(ii))
                AA(jj, 3) = CT
                AA(jj, 4) = DT
            end do

            if (NER .eq. SP) then
                call cgesv(4, 1, AA, 4, ipiv, BB, 4, flag)
            else
                call zgesv(4, 1, AA, 4, ipiv, BB, 4, flag)
            end if
            write(funt1_glb, '(f7.2, 12(2x, e26.17e3))'), Lmbd, Delta(ii), y(ii), Delta1(ii), y1(ii), C0(ii), real(BB(1:4))
            call destroy_pwamp(pwamp)
        end do

        close(funt1_glb)

    end subroutine solve_Dy_Q1Q2


    subroutine solve_1s0_Q1Q2()

        external cgesv, zgesv

        real(NER)           :: k
        character(len = 99) :: fname
        integer             :: ii, jj, flag
        complex(NEC)        :: L1_Q2, L2_Q2, CT, DT, ET, xi0, xi2, xi4 ! I0, I2, I4, J00, J02, J22,

        integer, dimension(3)           :: ipiv, uflag
        real(NER), dimension(2)         :: cnttrm_1s0
        real(NER), dimension(50)        :: C0, Lambda
        complex(NEC), dimension(2)      :: BB_Q1
        complex(NEC), dimension(3)      :: BB_Q2
        complex(NEC), dimension(2, 2)   :: AA_Q1
        complex(NEC), dimension(3, 3)   :: AA_Q2

        if (flagstrc%chnid .ne. CHNID_1S0_PP) return
        misc_opts%is_dibaryon = .false.
        fname = '3pt_1s0_phs.in'
        call read_nijphs(funt1_glb, fname, numnij)
        fname = 'slv_1s0_C0.in'
        call read2cols_io(funt1_glb, fname, Lambda, C0, numLambdas, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

        fname = 'slv_1s0_Q2.out'
        open (funt1_glb, file = fname, status = 'replace')
        write (funt1_glb, '(a, 3(f6.2, 2x))') '# fitted point : Tlab = ', nijphs(1:3, 1)
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda  C0  C1  D0  C2  D1  E0'

        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp(CHNID_1S0_PP, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
            cnttrm_1s0(1) = C0(ii)
            call initVM12_pwamp(2, pwamp)
            call set_Cs_pwamp(0, cnttrm_1s0, pwamp)
!             call set_VS_normal_PP(CHNID_1S0_PP, 0, cnttrm_1s0)

            do jj = 1, 2
                k = Tlab2k(nijphs(jj, 1))
                call update_k_pwamp(k, pwamp)
                call get_Tpq_pwamp(pwamp)
                call convert_TQn_to_phsn_quiet_pwamp(0, uflag, pwamp)
                BB_Q1(jj) = delta2t_nlo(nijphs(jj, 2) - pwamp%phsn(1), pwamp%phsn(1), k)
                call get_DCT_tvm_sngl(DT, CT, pwamp%tvms)
                AA_Q1(jj, 1) = CT; AA_Q1(jj, 2) = DT
            end do
            if (NER .eq. SP) then
                call cgesv(2, 1, AA_Q1, 2, ipiv, BB_Q1, 2, flag)
            else
                call zgesv(2, 1, AA_Q1, 2, ipiv, BB_Q1, 2, flag)
            end if

            do jj = 1, 3
                k = Tlab2k(nijphs(jj, 1))
                call update_k_pwamp(k, pwamp)
                call get_Tpq_pwamp(pwamp)
                call get_EDCT_tvm_sngl(ET, DT, CT, pwamp%tvms)
                AA_Q2(jj, 1) = CT; AA_Q2(jj, 2) = DT; AA_Q2(jj, 3) = ET
                pwamp%TQn(2) = CT*real(BB_Q1(1)) + DT*real(BB_Q1(2))
                call convert_TQn_to_phsn_quiet_pwamp(1, uflag, pwamp)
                BB_Q2(jj) = delta2t_nnlo(nijphs(jj, 2) - pwamp%phsn(1) - pwamp%phsn(2), pwamp%phsn(2), pwamp%phsn(1), k)

                call get_CD_xis_sngl(real(BB_Q1(1)), real(BB_Q1(2)), xi0, xi2, xi4, pwamp%tvms)
                BB_Q2(jj) = BB_Q2(jj) - (xi0*CT + xi2*DT + xi4*ET)

                call set_symmetric_nthVM_edge_tvm_sngl(VMID_VLNLO_tvm, VLNLO_tvm_sngl, pwamp%tvms)
                call get_L1_nthVM_tvm_sngl(VMID_VLNLO_tvm, L1_Q2, pwamp%tvms)
                call get_L2_nthVM_tvm_sngl(VMID_VLNLO_tvm, L2_Q2, pwamp%tvms)
                BB_Q2(jj) = BB_Q2(jj) - VLNLO_tvm_sngl(k, k, pwamp%tvms) - 2.0*L1_Q2 - L2_Q2
            end do

            if (NER .eq. SP) then
                call cgesv(3, 1, AA_Q2, 3, ipiv, BB_Q2, 3, flag)
            else
                call zgesv(3, 1, AA_Q2, 3, ipiv, BB_Q2, 3, flag)
            end if
            write(funt1_glb, '(f7.2, 8(2x, e26.17e3))'), Lmbd, C0(ii), real(BB_Q1(1:2)), real(BB_Q2(1:3))
            call destroy_pwamp(pwamp)
        end do

        close(funt1_glb)

    end subroutine solve_1s0_Q1Q2


! T_nlo = VLR_nlo(k, k) + 2*L1 + L2 + CT*C0

    subroutine solve_c0c1_Q2Q3()

        real(NER)           :: k
        integer             :: ii, flag
        complex(NEC)        :: L1_Q2, L2_Q2, L1_Q3, L2_Q3, CT, C0, C1, BB0, BB1
        character(len = 30) :: known_chns
        character(len = 99) :: fname
        integer, dimension(3)           :: uflag
        real(NER), dimension(50)        :: Lambda

        known_chns = '_1p1_3p1'
        if (index(known_chns, chnstr) .eq. 0) then
            print '(a, a)', 'do not know how to solve for C0 & C1 of ', chnstr
            return
        end if

        fname = '1pt_'//chnstr//'_phs.in'
        call read_nijphs(funt1_glb, fname, numnij)
        fname = 'slv_'//chnstr//'_Lambda.in'
        call read1col_io(funt1_glb, fname, Lambda, numLambdas, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

        fname = 'slv_'//chnstr//'_Q3.out'
        open (funt1_glb, file = fname, status = 'replace')
        write (funt1_glb, '(a, f6.2)') '# fitted point : Tlab = ', nijphs(1, 1)
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda  C0  C1'

        k = Tlab2k(nijphs(1, 1))
        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp(flagstrc%chnid, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
            call initVM012_pwamp(3, pwamp)
            call update_k_pwamp(k, pwamp)
            call get_allTQn_pwamp(0, pwamp)
            call convert_TQn_to_phsn_quiet_pwamp(0, uflag, pwamp)
            BB0 = delta2t_nlo(nijphs(1, 2) - pwamp%phsn(1), pwamp%phsn(1), k)
            BB1 = 0.0
            call set_symmetric_nthVM_edge_tvm_sngl(VMID_VLNLO_tvm, VLNLO_tvm_sngl, pwamp%tvms)
            call set_symmetric_nthVM_edge_tvm_sngl(VMID_VLNNLO_tvm, VLNNLO_tvm_sngl, pwamp%tvms)
            call get_L1_nthVM_tvm_sngl(VMID_VLNLO_tvm, L1_Q2, pwamp%tvms)
            call get_L2_nthVM_tvm_sngl(VMID_VLNLO_tvm, L2_Q2, pwamp%tvms)
            call get_L1_nthVM_tvm_sngl(VMID_VLNNLO_tvm, L1_Q3, pwamp%tvms)
            call get_L2_nthVM_tvm_sngl(VMID_VLNNLO_tvm, L2_Q3, pwamp%tvms)
            call get_CT_tvm_sngl(CT, pwamp%tvms)

            BB0 = BB0 - VLNLO_tvm_sngl(k, k, pwamp%tvms) - 2.0*L1_Q2 - L2_Q2
            BB1 = BB1 - VLNNLO_tvm_sngl(k, k, pwamp%tvms) - 2.0*L1_Q3 - L2_Q3
            C0 = BB0/CT
            C1 = BB1/CT
            write(funt1_glb, '(f7.2, 2(2x, e26.17e3))'), Lmbd, real(C0), real(C1)
            call destroy_pwamp(pwamp)
        end do

        close(funt1_glb)

    end subroutine solve_c0c1_Q2Q3


    subroutine solve_cde_cpld_Q2()

        integer             :: ii, flag
        real(NER)           :: k1, k2, k3
        complex(NEC)        :: expfac, expfac_d1, expfac_d2
        character(len = 30) :: known_chns
        character(len = 99) :: fname
        integer, dimension(3)           :: ipiv, uflag
        real(NER), dimension(10)        :: Cn
        real(NER), dimension(50)        :: Lambda, C0
        real(NER), dimension(2, 2)      :: vab
        real(NER), dimension(3, 3)      :: AA
        real(NER), dimension(3)         :: XX
        complex(NEC), dimension(2, 2)   :: L1_Q2, L2_Q2, BB, tab
        complex(NEC), dimension(3, 2, 2):: CDET

        known_chns = '_3s1_3p2'
        if (index(known_chns, chnstr) .eq. 0) then
            print '(a, a)', 'do not know how to solve for C, D & E of ', chnstr
            return
        end if

        fname = '2pt_'//chnstr//'_phs.in'
        call read_nijmxphs(funt1_glb, fname, numnij)
        fname = 'slv_'//chnstr//'_C0.in'
        call read2cols_io(funt1_glb, fname, Lambda, C0, numLambdas, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

        fname = 'slv_'//chnstr//'_Q2.out'
        open (funt1_glb, file = fname, status = 'replace')
        call output_pwa_pts(funt1_glb, numnij, nijmxphs(:, 1), nijmxphs(:, 2))
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda  C0  C2  D0  E0'
        Cn = 0.0

        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp(flagstrc%chnid, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
            Cn(1) = C0(ii)
            call set_Cs_pwamp(0, Cn, pwamp)
!             call set_VS_normal_PP(flagstrc%chnid, 0, Cn)
            call initVM012_pwamp(2, pwamp)

! the real and imaginary parts of T_nlo(1, 1) are independent, and must be used
! as separate inputs
! fit to d1 and epsilon at the first pwa point and d1 at the second
            k1 = Tlab2k(nijmxphs(1, 1))
            call update_k_pwamp(k1, pwamp)
            call get_allTQn_pwamp(0, pwamp)
            call convert_TQn_to_phsn_quiet_pwamp(0, uflag, pwamp)
            expfac_d1 = exp(-IMGUNT_NE*pwamp%phsn(1)*PIONINETY_NE)
            expfac_d2 = exp(-IMGUNT_NE*pwamp%phsn(2)*PIONINETY_NE)
            call get_L1L2_nthVM_tvm_cpld(VMID_VLNLO_tvm, VLNLO_tvm_cpld, L1_Q2, L2_Q2, pwamp%tvmc)
            call VLNLO_tvm_cpld(k1, k1, vab, pwamp%tvmc)
            BB = - vab - L1_Q2 - transpose(L1_Q2) - L2_Q2
            call get_EDCT_tvm_cpld(CDET(3, 1:2, 1:2), CDET(2, 1:2, 1:2), CDET(1, 1:2, 1:2), pwamp%tvmc)
            AA(1, 1:3) = aimag(-IMGUNT_NE*PI_NE*k1*expfac_d1*CDET(1:3, 1, 1))
            AA(2, 1:3) = real(-IMGUNT_NE*PI_NE*k1*expfac_d1*CDET(1:3, 1, 1))
            XX(1) = (nijmxphs(1, 2) - pwamp%phsn(1))*PIONINETY_NE*cos(pwamp%phsn(3)*PIONINETY_NE) + aimag(-IMGUNT_NE*PI_NE*k1*expfac_d1*BB(1, 1))
            XX(2) = -(nijmxphs(1, 4) - pwamp%phsn(3))*PIONINETY_NE*sin(pwamp%phsn(3)*PIONINETY_NE) + real(-IMGUNT_NE*PI_NE*k1*expfac_d1*BB(1, 1))

            k2 = Tlab2k(nijmxphs(2, 1))
            call update_k_pwamp(k2, pwamp)
            call get_allTQn_pwamp(0, pwamp)
            call convert_TQn_to_phsn_quiet_pwamp(0, uflag, pwamp)
            expfac_d1 = exp(-IMGUNT_NE*pwamp%phsn(1)*PIONINETY_NE)
            call get_L1L2_nthVM_tvm_cpld(VMID_VLNLO_tvm, VLNLO_tvm_cpld, L1_Q2, L2_Q2, pwamp%tvmc)
            call VLNLO_tvm_cpld(k2, k2, vab, pwamp%tvmc)
            BB = - vab - L1_Q2 - transpose(L1_Q2) - L2_Q2
            call get_EDCT_tvm_cpld(CDET(3, 1:2, 1:2), CDET(2, 1:2, 1:2), CDET(1, 1:2, 1:2), pwamp%tvmc)
            AA(3, 1:3) = aimag(-IMGUNT_NE*PI_NE*k2*expfac_d1*CDET(1:3, 1, 1))
            XX(3) = (nijmxphs(2, 2) - pwamp%phsn(1))*PIONINETY_NE*cos(pwamp%phsn(3)*PIONINETY_NE) + aimag(-IMGUNT_NE*PI_NE*k2*expfac_d1*BB(1, 1))

            if (NER .eq. SP) then
                call sgesv(3, 1, AA, 3, ipiv, XX, 3, flag)
            else
                call dgesv(3, 1, AA, 3, ipiv, XX, 3, flag)
            end if

            write(funt1_glb, '(f7.2, 4(2x, e26.17e3))'), Lmbd, C0(ii), XX(1:3)
            call destroy_pwamp(pwamp)
        end do

        close(funt1_glb)

    end subroutine solve_cde_cpld_Q2


    subroutine solve_ksw_cpld_Q3()

        integer             :: ii, flag
        real(NER)           :: k, C0
        character(len = 30) :: known_chns
        character(len = 99) :: fname
        integer, dimension(3)           :: uflag
        real(NER), dimension(1)         :: Cn
        real(NER), dimension(50)        :: Lambda, delta
        real(NER), dimension(2, 2)      :: vab
        complex(NEC), dimension(2, 2)      :: s3_fit, s3_bb

        known_chns = '_3p2'
        if (index(known_chns, chnstr) .eq. 0) then
            print '(a, a)', 'do not know how to solve for C, D & E of ', chnstr
            return
        end if

        fname = '1pt_'//chnstr//'_phs.in'
        call read_nijmxphs(funt1_glb, fname, numnij)
        fname = 'slv_'//chnstr//'_Lambda.in'
        call read1col_io(funt1_glb, fname, Lambda, numLambdas, flag)
        if (flag .ne. 0) then
            print '(a, a, a)', 'couldn''t open "', trim(fname), '"'
            return
        end if

        fname = 'slv_'//chnstr//'_pert_Q3.out'
        open (funt1_glb, file = fname, status = 'replace')
        call output_pwa_pts(funt1_glb, numnij, nijmxphs(:, 1), nijmxphs(:, 2))
        call output_regtype(funt1_glb, flagstrc%regtype)
        write (funt1_glb, '(a)') '# Lambda  C0'
        Cn = 0.0
        k  = Tlab2k(nijmxphs(1, 1))

        do ii = 1, numLambdas
            Lmbd = Lambda(ii)
            meshupper = meshul(flagstrc%regtype, Lmbd)
            call create_pwamp(flagstrc%chnid, flagstrc%regtype, Lmbd, N, meshupper, misc_opts, pwamp)
!~             Cn(1) = 0.0
!~             call set_VS_normal_PP(flagstrc%chnid, 3, Cn)
            call initVM012_pert_pwamp(3, pwamp)
            pwamp%tvmc%cms%C0 = 0.0
            call update_k_pwamp(k, pwamp)
            call get_allTQn_pwamp(3, pwamp)
            call convert_TQn_to_phsn_quiet_pwamp(2, uflag, pwamp)
            delta(1:9) = pwamp%phsn(1:9)
!~             delta(10)  = nijmxphs(1, 2) - sum(pwamp%phsn(1:9:3))
!~             delta(11)  = nijmxphs(1, 3) - sum(pwamp%phsn(2:9:3))
!~             delta(12)  = nijmxphs(1, 4) - sum(pwamp%phsn(3:9:3))
            delta(10)  = nijmxphs(1, 2) - (delta(1)+delta(4)+delta(7))
            delta(11)  = nijmxphs(1, 3) - (delta(2)+delta(5)+delta(8))
            delta(12)  = nijmxphs(1, 4) - (delta(3)+delta(6)+delta(9))
!~             !!! debug
!~             print *, delta(10:12)
            call get_nth_s_cpld(3, delta, s3_fit)
!~             s3_bb = -IMGUNT_NE*k*PI_NE*reshape(pwamp%TQn(13:16), (/2, 2/))
            C0 = (aimag(s3_fit(1, 1)) + k*PI_NE*pwamp%TQn(13))/(-regulator_tvm_aux(k, k, pwamp%tvmc%cms)*k*k*k*PI_NE)

            write(funt1_glb, '(f7.2, 2x, e26.17e3)'), Lmbd, C0
            call destroy_pwamp(pwamp)
        end do

        close(funt1_glb)

    end subroutine solve_ksw_cpld_Q3


end program slvcs


