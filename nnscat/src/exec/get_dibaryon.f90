! Bingwei Long 06/13/2013
! Compute properties of the 1S0 dibaryon

program get_dibaryon

    use drivers_pwphase
    implicit none

    integer, parameter :: funt = 23, numlmbds_max = 80, numcnts_max = 20, N = 100
    integer :: num_cnts, numks, chnid, uptoQn, nLmbds, flag, ii, jj, Nsteps, zz, aa, num_rsteps
    real(NER)           :: B0, B1, step_B, B, tmpsum, Delta, sigma_ysqr
    complex(NEC)        :: Sigma, Sigma1, Sinv, Sinv1
    type(opts_pwamp)    :: misc_opts
    type(neobj_pwamp)   :: pwamp
    logical             :: lgc_succ
    character(len = 3)     :: chnstr
    character(len = 64)    :: fname, cmdbuff

    real(NER), dimension(1:N)        :: Gamma
    real(NER), dimension(1:numlmbds_max)        :: Lambda
    real(NER), dimension(1:numlmbds_max, 1:numcnts_max+1)   :: Lambda_cnttrms = 0.0


    chnid = CHNID_1S0_PP
    chnstr = '1s0'
    uptoQn = 0
    misc_opts%is_dibaryon = .true.
    misc_opts%is_yukawa   = .true.
    call read_1s0_cnttrms('DIBARYON', uptoQn, funt, num_cnts, nLmbds, Lambda_cnttrms, lgc_succ)
    if (.not. lgc_succ) stop
    Lambda(1:nLmbds) = Lambda_cnttrms(1:nLmbds, 1)

    B0 = 24.0
    step_B = 1.0
    Nsteps = 20

    do ii = 1, nLmbds
        print '(a, f8.1)', 'Lambda = ', Lambda(ii)
        call setup_pwamp(chnid, uptoQn, REGTYPE_GAUSSIAN, Lambda(ii), N, Lambda(ii)*PC_ratio_Lambda_Mambda,     &
            misc_opts, Lambda_cnttrms(ii, 2:num_cnts+1), num_cnts, pwamp)
        Delta = Lambda_cnttrms(ii, 2)
        sigma_ysqr = Lambda_cnttrms(ii, 3) 
        do jj = 1, Nsteps
            B = B0 + step_B*(jj-1)
            B1 = B - 0.01
            call get_GreenFunction_tvm_sngl(complex(-PC_mN*B, 0.0), Sigma, pwamp%tvms)
            call get_GreenFunction_tvm_sngl(complex(-PC_mN*B1, 0.0), Sigma1, pwamp%tvms)
            Sinv  = (-B+Delta)/sigma_ysqr - Sigma
            Sinv1 = (-B1+Delta)/sigma_ysqr - Sigma1
            write (*, *) B, real((Sinv1 - Sinv)/(-B1 + B)), real(Sinv)
!             write (*, *) B, abs((-B+Delta)/sigma_ysqr - Sigma)
        end do
        call destroy_pwamp(pwamp)
    end do

!     contains

!     function v1s0_Q1(p1, p2)

!         real(NER)   :: v1s0_Q1
!         real(NER), intent(in)   :: p1, p2

!         v1s0_Q1 = C0_1S0_PP + C1_1S0_PP
        
!     end function v1s0_Q1


end program get_dibaryon
