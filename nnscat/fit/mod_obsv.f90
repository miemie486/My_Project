module mod_obsv

! //--------------------------------------------------------------------------//
! Given phase shifts as inputs, calculate various observables
! //--------------------------------------------------------------------------//

    use nneft_type
    use util_mathrle
    implicit none

    real(NER), parameter :: MeV_to_barn = 389379.65_NER

! -----\\ References \\-----
!
! WGQF: Walter Gl\"ockle, The Quantum Mechanical Few_Body Problem
!
! BLW78: J. Bystricky, F. Lehar, and P. Winternitz, Formalism of nucleon-nucleon
! elastic scattering experiments, Journal de Physique, 1978, 39 (1), pp.1-32
!
! NPARA13: R. Navarro Perez, J. E. Amaro, and E. Ruiz Arriola, Coarse-grained
! potential analysis of neutron-proton and proton-proton scattering below the
! pion production threshold, Phys. Rev. C 88, 064002 (2013)

! --------------------------

    ! Structure defined by MeiYin, now obsolete
    type :: phaseDat
        real(NER)                   ::  energy
        real(NER), allocatable      ::  unphld0(:), unphld1(:), cpld(:,:)
    end type phaseDat

    ! Structure that encapsulates all the physical constants, such as the
    ! average nucleon mass, proton mass, neutron mass, pion mass, etc.
    type struct_phys_const
        ! mN: average nucleon mass
        real(NER)   :: mN

    end type struct_phys_const

    type :: struct_Expt_Data
        character(len=6)    :: Ref
        character(len=4)    :: obsName
        real(NER)           :: acm, elab, val, error
        integer				:: ii    ! ii is the the sequence number of elab list.
    end type struct_Expt_Data

    integer, parameter :: MAXVALJ_OBS = 30

    ! struct_phase_allwaves includes phase shifts at Elab for channels below
    ! certain cutoff J's. Jmax0 is the cutoff J for S=0 (uncoupled); Jmax1 is
    ! the cutoff J for S=1 (uncoupled); Jmaxcp is the cutoff J for coupled
    ! channels. Note that the arrays for storing phase shifts are *not*
    ! dynamically allocated, for the purpose of simplicity.
    ! For S = 0:
    ! struct_phase_allwaves%unps0(J) = \theta{^1(J)_J}.
    ! For example, when J = 0
    ! struct_phase_allwaves%unps0(1) = \theta(^1S_0)
    ! For S = 1:
    ! struct_phase_allwaves%unps1(J) = \theta ^3(L)_J. When J > 0, L = J; when
    ! J = 0, L = 1.
    ! For example, when J = 0
    ! struct_phase_allwaves%unps1(0) = \theta ^3(P)_0.
    ! When J = 1,
    ! struct_phase_allwaves%unps1(1) = \theta ^3(P)_1.
    ! When J >= 1,
    ! cpldps(J, 1) = \theta ^3(J-1)_J;
    ! cpldps(J, 2) = \theta ^3(J+1)_J;
    ! cpldps(J, 3) = mixing angle
    type :: struct_phase_allwaves

        real(NER)      ::  Elab
        integer        ::  Jmax0, Jmax1, Jmaxcp
        real(NER)      ::  unps0(0:MAXVALJ_OBS), unps1(0:MAXVALJ_OBS), cpldps(1:MAXVALJ_OBS, 1:3)

    end type struct_phase_allwaves

    ! Structure defined by MeiYin, now obsolete
    ! Wolfenstein parameters at Elab
    ! Reference: WGQF Eq. (2.53).
    ! Parity conservation and isospin invariance assumed
    ! x: CM angle
    type :: struct_Wlfp

        real(NER)   :: Elab, x
        complex(NEC) :: a, c, m, g, h

    end type struct_Wlfp

    type :: struct_Saclay

        real(NER)  :: Elab, x
        complex(NEC) :: a, b, c, d, e

    end type struct_Saclay

contains

    ! Given lab kinetic energy Tl, calculating the CM momentum k, using the
    ! nonrelativistic energy-momentum relation
    real(NER) function Tlab2kcm_nr_isoinv(Tl, mN)

        real(NER), intent(in)   :: Tl, mN

        Tlab2kcm_nr_isoinv = sqrt(mN * Tl * 0.5_NER)

    end function Tlab2kcm_nr_isoinv

    real(NER) function Tlab2kcm(Tl, mN)

        real(NER), intent(in)   :: Tl, mN

        Tlab2kcm = Tlab2kcm_nr_isoinv(Tl, mN)

    end function Tlab2kcm

    ! Calculating S matrix elements SmatrixEle_obs from phase shifts
    ! Reference: WGQF Eqs. (2.139) and (2.140)
    ! SmatrixEle_obs = S_l'sls^{J}
    ! uncpld_phase0: uncoupled phase shift for S=0;
    ! uncpld_phase1: uncoupled phase shift for s=1;
    ! cpld_phase: coupled phase shift
    ! The data layout of uncpld_phase0, uncpld_phase1 and cpld_phase are the
    ! same as that of similarly named members of struct_phase_allwaves
    complex(NEC) function SmatrixEle_obs(J, l, lp, s, uncpld_phase0, uncpld_phase1, cpld_phase)

        integer,intent(in)          :: J, lp, l, s
        real(NER), intent(in)       :: uncpld_phase0(:), uncpld_phase1(:), cpld_phase(:, :)

        complex(NEC)                :: n

        n = (0.0_NER, 1.0_NER)
        if (s == 0 .and. (l == J  .and. lp == J) .and. J >= 0) then
            SmatrixEle_obs = exp(2.0_NER * n * deg2rad(uncpld_phase0(J+1)))
                    write(*,*) "s",SmatrixEle_obs
        else if(((J >= 1 .and. l == J  .and. lp == J) .or. (J == 0 .and. l == 1 .and. lp == 1)).and. s == 1) then
            SmatrixEle_obs = exp(2.0_NER * n * deg2rad(uncpld_phase1(J+1)))
        else if(J >= 1 .and. lp == J-1 .and. l == J-1 .and. s == 1) then
            SmatrixEle_obs = cos(2.0_NER*deg2rad(cpld_phase(J,3))) * exp(2.0_NER*n*deg2rad(cpld_phase(J,1)))
        else if(J >= 1 .and. lp == J+1 .and. l == J-1 .and. s == 1) then
            SmatrixEle_obs = n * sin(2.0_NER*deg2rad(cpld_phase(J,3))) * exp(n*(deg2rad(cpld_phase(J,1))+deg2rad(cpld_phase(J,2))))
        else if(J >= 1 .and. lp == J-1 .and. l == J+1 .and. s == 1) then
            SmatrixEle_obs = n * sin(2.0_NER*deg2rad(cpld_phase(J,3))) * exp(n*(deg2rad(cpld_phase(J,1))+deg2rad(cpld_phase(J,2))))
        else if(J >= 1 .and. lp == J+1 .and. l == J+1 .and. s == 1) then
            SmatrixEle_obs = cos(2.0_NER*deg2rad(cpld_phase(J,3))) * exp(2.0_NER*n*deg2rad(cpld_phase(J,2)))
        else
            SmatrixEle_obs = 0.0_NER
            write(*,*) "smatrixele",j,l,lp,s,"=0"
        end if

    end function SmatrixEle_obs

    ! Calculating M-matrix elements Mm_s_obs from phase shifts, by using
    ! SmatrixEle_obs.
    ! Reference: WGQF Eq. (2.115)
    ! Mm_s_obs = M_{{m_s}'{m_s}}^S
    ! s: total spin
    ! x: CM angle
    ! mN: mass of average nucleon
    ! paw: phase shifts used as the inputs

    ! There are two function Mm_s_obs below, they both work well. The first Mm_s_obs do not call function SmatrixEle_obs

    complex(NEC) function Mm_s_obs(s, m_sp, m_s, x, mN, paw)  ! m_sp is m_s'

        integer, intent(in)                     :: s, m_sp, m_s
        real(NER), intent(in)                   :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        integer                             :: ls, l, J
        complex(NEC)                        :: n, Mm
        real(NER)                           :: k

        complex(NEC)                        :: Smatrix(2,2),L1(1,2),L2(2,1)
        complex(NEC)                        :: Mcp(1,1), Mun

        k = Tlab2kcm(paw%Elab, mN)
        n = (0.0_NER,1.0_NER)
        Mm_s_obs = (0.0_NER, 0.0_NER)
        Mm = (0.0_NER, 0.0_NER)

        if (s == 0) then


            do J = 0 , paw%Jmax0

                Mm = (exp(2.0_NER * n * deg2rad(paw%unps0(J)))-1.0_NER) * spherhms(J, 0, Cos(deg2rad(x))) *sqrt(2.0_NER*J+1.0_NER)
                Mm_s_obs = Mm_s_obs + Mm

            end do

            Mm_s_obs = Mm_s_obs * sqrt(PI_NE) / (n*k)


        else if(s == 1) then

            Mm_s_obs = (exp(2.0_NER * n * deg2rad(paw%unps1(0))) - 1.0_NER) * spherhms(1, m_s-m_sp, cos(deg2rad(x))) * &
                       & CG_cff(1,1,0, m_s-m_sp, m_sp) * CG_cff(1,1,0, 0, m_s) * sqrt(3.0_NER)

            do J = 1, max(paw%Jmax0, paw%Jmax1, paw%Jmaxcp)

                L1(1,1) = spherhms(J-1, m_s-m_sp, cos(deg2rad(x))) * CG_cff(J-1,1,J,m_s-m_sp,m_sp)
                L1(1,2) = spherhms(J+1, m_s-m_sp, cos(deg2rad(x))) * CG_cff(J+1,1,J,m_s-m_sp,m_sp)

                Smatrix(1,1) = cos(2.0_NER*deg2rad(paw%cpldps(J,3))) * exp(2.0_NER*n*deg2rad(paw%cpldps(J,1))) -1.0_NER
                Smatrix(1,2) = -n*sin(2.0_NER*deg2rad(paw%cpldps(J,3)))*exp(n*(deg2rad(paw%cpldps(J,1))+deg2rad(paw%cpldps(J,2))))
                Smatrix(2,1) = Smatrix(1,2)
                Smatrix(2,2) = cos(2.0_NER*deg2rad(paw%cpldps(J,3))) * exp(2.0_NER*n*deg2rad(paw%cpldps(J,2))) -1.0_NER

                L2(1,1) = sqrt(2.0_NER*(J-1)+1.0_NER) * CG_cff(J-1, 1, J, 0, m_s)
                L2(2,1) = sqrt(2.0_NER*(J+1)+1.0_NER) * CG_cff(J+1, 1, J, 0, m_s)

                Mcp = matmul(L1,matmul(Smatrix,L2))
                Mun = (exp( 2.0_NER * n * deg2rad(paw%unps1(J))) - 1.0_NER)* &
                      & spherhms(J, m_s-m_sp, cos(deg2rad(x))) * CG_cff(J,1,J,m_s-m_sp,m_sp) * &
                      & sqrt(2.0_NER*J+1.0_NER) * CG_cff(J, 1, J, 0, m_s)
                Mm_s_obs = Mm_s_obs + Mcp(1,1) + Mun

            end do

            Mm_s_obs = Mm_s_obs * sqrt(PI_NE) / (n*k)

        end if

        end function Mm_s_obs


!     complex(NEC) function Mm_s_obs(s, m_sp, m_s, x, mN, paw)

!         integer, intent(in)                     :: s, m_sp, m_s
!         real(NER), intent(in)                   :: x, mN
!         type(struct_phase_allwaves), intent(in) :: paw

!         integer                             :: ls, l, J
!         complex(NEC)                        :: n, Mm
!         real(NER)                           :: k

!         k = Tlab2kcm(paw%Elab, mN)
!         n = (0.0_NER,1.0_NER)
!         Mm_s_obs = (0.0_NER, 0.0_NER)
!         Mm = (0.0_NER, 0.0_NER)

!         if (s == 0) then
!             do J = 0, paw%Jmax0
!                 Mm = CG_cff(J, 0, J, 0, 0) * spherhms(J, 0, cos(deg2rad(x))) * &
!                 &(SmatrixEle_obs(J, J, J, 0, paw%unps0, paw%unps1, paw%cpldps) &
!                     & - 1.0_NER) * CG_cff(J, 0, J, 0, 0) * sqrt(PI_NE*(2.0_NER*J+1.0_NER))/(n*k)
!                 Mm_s_obs = Mm_s_obs + Mm
!             end do
!         else if(s == 1) then
!                 do J = 0, max(paw%Jmax0, paw%Jmax1, paw%Jmaxcp)
!                     do l = abs(J-1), J+1
!                         do ls = abs(J-1), J+1
!                             if(valid_lsj(l,ls,J)) then

!                             Mm = CG_cff(ls, 1, J, m_s-m_sp, m_sp) * &
!                             & spherhms(ls, m_s-m_sp, cos(deg2rad(x))) * &
!                             & (SmatrixEle_obs(J, ls,l, 1, paw%unps0, paw%unps1, paw%cpldps)  - dirac_func(ls, l)) * &
!                             & CG_cff(l, 1, J, 0, m_s) * &
!                             & n**(l-ls) * &
!                             & sqrt(PI_NE*(2.0_NER * l + 1.0_NER)) /(n*k)

!                                 if (ABS(Mm ) .LT. 0.000000000000001) then
!                                     write (40,*) "mm=0", J, l, ls, m_s, m_sp, Mm
!                                 else
!                                     write(40,*)"mm", J, l, ls, m_s, m_sp, Mm
!                                 end if
!                             else
!                             Mm = 0.0_NER
!                             end if

!                             Mm_s_obs = Mm_s_obs+Mm

!                         end do
!                     end do
!                 end do
!         end if
!         write(40,*) "MMM", s, m_sp, m_s , Mm_s_obs ,paw%Jmax0
!     end function Mm_s_obs

!-----------
! Calculating Wolfeinstein parameters from phase shifts
! Reference: WGQF Eq.(2.96)
! According to Mei Yin, the fifth line of Eq. (2.96) has a wrong sign
! h = 1/8 cos\theta (M^1_{11} + M^1_{1 -1} + ...
! --> h = 1/8 cos\theta (M^1_{11} - M^1_{1 -1} + ...
! x: CM angle


    function valid_lsj(l,lp,j)

        logical :: valid_lsj
        integer, intent(in) :: l,lp,j

        if (l==lp .and. (l==J .or. l==J+1 .or. l==J-1) ) then
        valid_lsj = .true.
        write (40,*) "T", l,lp,j
        else if ( (l==J+1 .and. lp==J-1 ) .or. (l==J-1 .and. lp==J+1) ) then
        valid_lsj = .true.
        write (40,*) "T", l,lp,j
        else
        valid_lsj = .false.
        write (40,*) "F", l,lp,j
        end if

    end function valid_lsj


    subroutine Wlfp_para_obs(Wlfp, x, mN, paw)

        type(struct_Wlfp), intent(out)  :: Wlfp
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Wlfp%Elab = paw%Elab
        Wlfp%x = x
        call Wlfp_pa_obs(Wlfp%a, x, mN, paw)
        call Wlfp_pm_obs(Wlfp%m, x, mN, paw)
        call Wlfp_pc_obs(Wlfp%c, x, mN, paw)
        call Wlfp_pg_obs(Wlfp%g, x, mN, paw)
        call Wlfp_ph_obs(Wlfp%h, x, mN, paw)

    end subroutine Wlfp_para_obs

    subroutine Wlfp_pa_obs(Wlfp_a, x, mN, paw)

        complex(NEC),intent(out)        :: Wlfp_a
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Wlfp_a = (Mm_s_obs(1, 1, 1, x, mN, paw) + Mm_s_obs(1, -1, -1, x, mN, paw) + Mm_s_obs(0, 0, 0, x, mN, paw) &
            & + Mm_s_obs(1, 0, 0, x, mN, paw))/4.0_NER

    end subroutine Wlfp_pa_obs

    subroutine Wlfp_pm_obs(Wlfp_m, x, mN, paw)

        complex(NEC),intent(out)        :: Wlfp_m
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Wlfp_m = (-Mm_s_obs(1, 1, -1, x, mN, paw) - Mm_s_obs(1, -1, 1, x, mN, paw) + Mm_s_obs(1, 0, 0, x, mN, paw) &
            & - Mm_s_obs(0, 0, 0, x, mN, paw))/4.0_NER

    end subroutine Wlfp_pm_obs

    subroutine Wlfp_pc_obs(Wlfp_c, x, mN, paw)

        complex(NEC),intent(out)        :: Wlfp_c
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        complex(NEC)        :: n

        n = (0.0_NER, 1.0_NER)
        Wlfp_c = n * sqrt(2.0_NER) * (Mm_s_obs(1, 1, 0, x, mN, paw)- Mm_s_obs(1, 0, 1, x, mN, paw))/4.0_NER

    end subroutine Wlfp_pc_obs

    subroutine Wlfp_pg_obs(Wlfp_g, x, mN, paw)

        complex(NEC),intent(out)        :: Wlfp_g
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Wlfp_g = (Mm_s_obs(1, 1, 1, x, mN, paw) + Mm_s_obs(1, 1, -1, x, mN, paw) -  Mm_s_obs(0, 0, 0, x, mN, paw))/4.0_NER

    end subroutine Wlfp_pg_obs

    subroutine Wlfp_ph_obs(Wlfp_h, x, mN, paw)

        complex(NEC),intent(out)        :: Wlfp_h
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Wlfp_h = (Mm_s_obs(1, 1, 1,x, mN, paw) - Mm_s_obs(1, 1, -1, x, mN, paw) &
                &- Mm_s_obs(1, 0, 0, x, mN, paw)) * cos(deg2rad(x))/4.0_NER + &
                &(Mm_s_obs(1, 1, 0, x , mN, paw) + Mm_s_obs(1, 0, 1, x, mN, paw)) * &
                &sin(deg2rad(x))/(2.0_NER * sqrt(2.0_NER))

    end subroutine Wlfp_ph_obs

!calculation of differential cross section; NOTE

!In the code we use natrual units so that unit of M-matrix element is MeV^(-)1 and unit of cross section is MeV^(-2).
!In experiments the unit is mbarn(1mbarn = 10^-3 barn = 10^-31 m^2).
!c^2 * hbar^2 = (299792458 m/s)^2 * (6.58211951e-22 MeV*s) = 3.8937965e-26 m^2 / MeV^(-2) = 389379.65 mbarn / MeV(-2).
!So if we want to change unit from MeV(-2) to mbarn, we should times 389379.65. We define MeV_to_barn=389379.65

    function diff_crosec_obs(mN, Wlfp)

        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp
        real(NER)                       :: diff_crosec_obs

        diff_crosec_obs = MeV_to_barn * (abs(Wlfp%a)**2 + 2.0_NER*abs(Wlfp%c)**2 + &
                            &abs(Wlfp%m)**2_NER + 2.0_NER*abs(Wlfp%g)**2 + 2.0_NER*abs(Wlfp%h)**2)

    end function diff_crosec_obs


!calculation of polarization(P); NOTE

    subroutine plrz_obs(polariz, mN, Wlfp)

        real(NER), intent(out)          :: polariz
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        polariz = MeV_to_barn * 2.0_NER * real(conjg(Wlfp%c) * (Wlfp%a+Wlfp%m))/ diff_crosec_obs(mN, Wlfp)

    end subroutine plrz_obs

!calculation of depolarization(D); NOTE

    subroutine deplrz_obs(depolariz, mN, Wlfp)

        real(NER), intent(out)          :: depolariz
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        depolariz = MeV_to_barn * (abs(Wlfp%a)**2 + abs(Wlfp%m)**2 + 2.0_NER*abs(Wlfp%c)**2&
                    & - 2.0_NER*abs(Wlfp%g)**2 - 2.0_NER*abs(Wlfp%h)**2)&
            &/diff_crosec_obs(mN, Wlfp)

    end subroutine deplrz_obs

    subroutine deplrz_tr_obs(depolarizTr, mN, Wlfp)

        real(NER), intent(out)          :: depolarizTr
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        depolarizTr = 1.0_NER - MeV_to_barn * (abs(Wlfp%a-Wlfp%m)**2 + 4*abs(Wlfp%h)**2)/diff_crosec_obs(mN, Wlfp)

    end subroutine deplrz_tr_obs


!calculation of spin parameter R'; NOTE

    subroutine spinp_Rp_obs(SpRp, mN, Wlfp)

        real(NER), intent(out)          :: SpRp
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        SpRp = MeV_to_barn * (sin(deg2rad(Wlfp%x)/2.0_NER) * (abs(Wlfp%a)**2 - abs(Wlfp%m)**2&
                & - abs(Wlfp%g-Wlfp%h)**2 + abs(Wlfp%g+Wlfp%h)&
                &**2) + cos(deg2rad(Wlfp%x)/2.0_NER) * aimag(2.0_NER*Wlfp%c * &
                &(conjg(Wlfp%a) - conjg(Wlfp%m))))/diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_Rp_obs

!calculation of spin parameter R; NOTE

    subroutine spinp_R_obs(SpR, mN, Wlfp)

        real(NER), intent(out)          :: SpR
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        SpR = MeV_to_barn * (cos(deg2rad(Wlfp%x)/2.0_NER) * (abs(Wlfp%a)**2 - abs(Wlfp%m)**2&
            & + abs(Wlfp%g-Wlfp%h)**2 - abs(Wlfp%g+Wlfp%h)&
            &**2) - sin(deg2rad(Wlfp%x)/2.0_NER) * aimag(2.0_NER * Wlfp%c * &
            &(conjg(Wlfp%a) - conjg(Wlfp%m))))/diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_R_obs

!calculation of spin parameter A'; NOTE

    subroutine spinp_Ap_obs(SpAp, mN, Wlfp)

        real(NER), intent(out)          :: SpAp
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        SpAp = MeV_to_barn * (cos(deg2rad(Wlfp%x)/2.0_NER)*(abs(Wlfp%a)**2 - &
                &abs(Wlfp%m)**2 - abs(Wlfp%g-Wlfp%h)**2 + abs(Wlfp%g+Wlfp%h)&
            &**2) - sin(deg2rad(Wlfp%x)/2.0_NER) * aimag(2.0_NER * Wlfp%c * &
            &(conjg(Wlfp%a) - conjg(Wlfp%m))))/diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_Ap_obs

!calculation of spin parameter A; NOTE

    subroutine spinp_A_obs(SpA, mN, Wlfp)

        real(NER), intent(out)          :: SpA
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        SpA = MeV_to_barn * (-sin(deg2rad(Wlfp%x)/2.0_NER)*(abs(Wlfp%a)**2 - &
            &abs(Wlfp%m)**2 + abs(Wlfp%g-Wlfp%h)**2 - abs(Wlfp%g+Wlfp%h)&
            &**2) - cos(deg2rad(Wlfp%x)/2.0_NER) * aimag(2.0_NER * Wlfp%c * (conjg(Wlfp%a)&
            & - conjg(Wlfp%m))))/diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_A_obs

!calculation of spin parameter for target:At, Rt, Atp, Rtp; NOTE

    subroutine spinp_At_obs(SpAt, mN, Wlfp)

        real(NER), intent(out)          :: SpAt
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        complex(NEC)                    :: n

        n = (0.0_NER, 1.0_NER)

        SpAt = MeV_to_barn * (4.0_NER * real(n * Wlfp%c * conjg(Wlfp%g)) * sin(deg2rad(Wlfp%x)/2.0_NER)&
            & + 2.0_NER * real((Wlfp%a + Wlfp%m) * conjg(Wlfp%g) + (Wlfp%a&
            &- Wlfp%m) * conjg(Wlfp%h)) * cos(deg2rad(Wlfp%x)/2.0_NER)) /diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_At_obs

    subroutine spinp_Rt_obs(SpRt, mN, Wlfp)

        real(NER), intent(out)          :: SpRt
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        complex(NEC)                    :: n

        n = (0.0_NER, 1.0_NER)
        SpRt = MeV_to_barn * (-4.0_NER * real(n * Wlfp%c * conjg(Wlfp%g)) * cos(deg2rad(Wlfp%x)/2.0_NER)&
            & + 2.0_NER * real((Wlfp%a + Wlfp%m) * conjg(Wlfp%g)&
            & + (Wlfp%a - Wlfp%m) * conjg(Wlfp%h)) * sin(deg2rad(Wlfp%x)/2.0_NER)) /diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_Rt_obs

    subroutine spinp_Apt_obs(SpApt, mN, Wlfp)

        real(NER), intent(out)          :: SpApt
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        complex(NEC)                    :: n

        n = (0.0_NER, 1.0_NER)
        SpApt = MeV_to_barn * (4.0_NER * real(n * Wlfp%c * conjg(Wlfp%g)) * cos(deg2rad(Wlfp%x)/2.0_NER)&
            & - 2.0_NER * real((Wlfp%a + Wlfp%m) * conjg(Wlfp%g)&
            & - (Wlfp%a - Wlfp%m) * conjg(Wlfp%h)) * sin(deg2rad(Wlfp%x)/2.0_NER)) /diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_Apt_obs

    !Acorrding to Hoshizaki, the first term of At,Rt,Apt is Re(icg*), Meiying's code is Re(icg). Now it is modified.

    subroutine spinp_Rpt_obs(SpRpt, mN, Wlfp)

        real(NER), intent(out)          :: SpRpt
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        complex(NEC)                    :: n

        n = (0.0_NER, 1.0_NER)
        SpRpt = MeV_to_barn * (4.0_NER * real(n * Wlfp%c * Wlfp%g) * sin(deg2rad(Wlfp%x)&
            &/2.0_NER) + 2.0_NER * real((Wlfp%a + Wlfp%m) * conjg(Wlfp%g)&
            & - (Wlfp%a- Wlfp%m) * conjg(Wlfp%h)) * cos(deg2rad(Wlfp%x)/2.0_NER)) /diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_Rpt_obs


!calculation of spin parameter Axx; NOTE;

    subroutine spinp_Axx_obs(SpAxx, mN, Wlfp)

        real(NER), intent(out)          :: SpAxx
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        complex(NEC)                    :: n

        n = (0.0_NER, 1.0_NER)
        SpAxx = MeV_to_barn * 2.0_NER * real((Wlfp%a - Wlfp%m) * conjg(Wlfp%g) - (Wlfp%a + Wlfp%m)&
            & * conjg(Wlfp%h) * cos(deg2rad(Wlfp%x)) - 2.0_NER*&
            &n * Wlfp%c * conjg(Wlfp%h) * sin(deg2rad(Wlfp%x))) / diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_Axx_obs

!calculation of spin parameter Azz; NOTE;

    subroutine spinp_Azz_obs(SpAzz, mN, Wlfp)

        real(NER), intent(out)          :: SpAzz
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        complex(NEC)                    :: n

        n = (0.0_NER, 1.0_NER)
        SpAzz = MeV_to_barn * 2.0_NER * real((Wlfp%a - Wlfp%m) * conjg(Wlfp%g) + &
            &(Wlfp%a + Wlfp%m) * conjg(Wlfp%h) * cos(deg2rad(Wlfp%x)) + 2.0_NER* &
            &n * Wlfp%c * conjg(Wlfp%h) * sin(deg2rad(Wlfp%x))) / diff_crosec_obs(mN, Wlfp)

    end subroutine spinp_Azz_obs

!calculation pf spin parameter Axz; NOTE;

    subroutine spinp_Axz_obs(SpAxz, mN, Wlfp)

        real(NER), intent(out)          :: SpAxz
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        complex(NEC)                    :: n

        n = (0.0_NER, 1.0_NER)
        SpAxz = -(MeV_to_barn * 2.0_NER * real((Wlfp%a + Wlfp%m) * conjg(Wlfp%h) * sin(deg2rad(Wlfp%x))&
            & - 2.0_NER * n * Wlfp%c * conjg(Wlfp%h) * &
            &cos(deg2rad(Wlfp%x)))/diff_crosec_obs(mN, Wlfp))

    end subroutine spinp_Axz_obs

    subroutine spinp_Ayy_obs(SpAyy, mN, Wlfp)

        real(NER), intent(out)          :: SpAyy
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        SpAyy = 1 - MeV_to_barn * ((abs(Wlfp%a-Wlfp%m)**2 + 4.0_NER*abs(Wlfp%g)**2)/diff_crosec_obs(mN, Wlfp))

    end subroutine spinp_Ayy_obs

    subroutine SGT_obs(SGT, mN, Wlfp)

        real(NER),intent(out)           :: SGT
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        real(NER)                       :: k

        k = Tlab2kcm(Wlfp%Elab, mN)
        SGT = MeV_to_barn * 4.0_NER * PI_NE * (aimag(Wlfp%a)) / k

    end subroutine SGT_obs

    subroutine SGTL_obs(SGTL, mN, Wlfp)

        real(NER),intent(out)           :: SGTL
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        real(NER)                       :: k

        k = Tlab2kcm(Wlfp%Elab, mN)
        SGTL = -MeV_to_barn * 8.0_NER * PI_NE * (aimag(Wlfp%g-Wlfp%h)+2.0_NER*aimag(Wlfp%h)) / k

    end subroutine SGTL_obs

    subroutine SGTT_obs(SGTT, mN, Wlfp)

        real(NER),intent(out)           :: SGTT
        real(NER),intent(in)            :: mN
        type(struct_Wlfp), intent(in)   :: Wlfp

        real(NER)                       :: k

        k = Tlab2kcm(Wlfp%Elab, mN)
        SGTT = -MeV_to_barn * 8.0_NER * PI_NE * (aimag(Wlfp%g-Wlfp%h)) / k

    end subroutine SGTT_obs


! Below is Saclay parameter

    subroutine Saclay_para_obs(Saclay, x, mN, paw)

        type(struct_Saclay), intent(out)  :: Saclay
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        saclay%elab = paw%Elab
        Saclay%x = x

        call Saclay_pa_obs(saclay%a, x, mN, paw)
        call Saclay_pb_obs(saclay%b, x, mN, paw)
        call Saclay_pc_obs(saclay%c, x, mN, paw)
        call Saclay_pd_obs(saclay%d, x, mN, paw)
        call Saclay_pe_obs(saclay%e, x, mN, paw)

    end subroutine Saclay_para_obs

    subroutine Saclay_pa_obs(Saclay_a, x, mN, paw)

        complex(NEC),intent(out)        :: Saclay_a
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Saclay_a = (Mm_s_obs(1, 1, 1, x, mN, paw) + Mm_s_obs(1, 0, 0, x, mN, paw) - Mm_s_obs(1, 1, -1, x, mN, paw))/2.0_NER

    end subroutine Saclay_pa_obs

    subroutine Saclay_pb_obs(Saclay_b, x, mN, paw)

        complex(NEC),intent(out)        :: Saclay_b
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Saclay_b = (Mm_s_obs(1, 1, 1, x, mN, paw) + Mm_s_obs(0, 0, 0, x, mN, paw) + Mm_s_obs(1, 1, -1, x, mN, paw))/2.0_NER

    end subroutine Saclay_pb_obs

    subroutine Saclay_pc_obs(Saclay_c, x, mN, paw)

        complex(NEC),intent(out)        :: Saclay_c
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Saclay_c = (Mm_s_obs(1, 1, 1, x, mN, paw) - Mm_s_obs(0, 0, 0, x, mN, paw) + Mm_s_obs(1, 1, -1, x, mN, paw))/2.0_NER

    end subroutine Saclay_pc_obs

    subroutine Saclay_pd_obs(Saclay_d, x, mN, paw)

        complex(NEC),intent(out)        :: Saclay_d
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw

        Saclay_d = -(Mm_s_obs(1, 1, 0, x, mN, paw) + Mm_s_obs(1, 0, 1, x, mN, paw))/(sqrt(2.0_NER)*sin(deg2rad(x)))

    end subroutine Saclay_pd_obs

    subroutine Saclay_pe_obs(Saclay_e, x, mN, paw)

        complex(NEC),intent(out)        :: Saclay_e
        real(NER),intent(in)            :: x, mN
        type(struct_phase_allwaves), intent(in) :: paw
        complex(NEC)        :: n = (0.0_NER,1.0_NER)

        Saclay_e = n *(Mm_s_obs(1, 1, 0, x, mN, paw) - Mm_s_obs(1, 0, 1, x, mN, paw)) /sqrt(2.0_NER)

    end subroutine Saclay_pe_obs

    function DSG_sac(mN, Saclay)

        real(NER),intent(in)            :: mN
        type(struct_Saclay), intent(in)   :: Saclay
        real(NER)                       :: DSG_sac

        DSG_sac = MeV_to_barn * (abs(saclay%a)**2 + abs(saclay%b)**2 + &
                            &abs(saclay%c)**2 + abs(saclay%d)**2 + abs(saclay%e)**2) / 2.0_NER

    end function DSG_sac



    function P_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: P_sac

    	P_sac = real(conjg(saclay%a) * saclay%e)/ I0000(mN, Saclay)

    end function P_sac

    function D_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: D_sac

        D_sac = (abs(saclay%a)**2 + abs(saclay%b)**2 - abs(saclay%c)**2 &
                 - abs(saclay%d)**2 + abs(saclay%e)**2)/2.0_NER/I0000(mN, Saclay)

    end function D_sac

    function AYY_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: AYY_sac

    	AYY_sac = (abs(saclay%a)**2 - abs(saclay%b)**2 - abs(saclay%c)**2 &
                & + abs(saclay%d)**2 + abs(saclay%e)**2)/2.0_NER/I0000(mN, Saclay)

    end function AYY_sac

    function DT_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: DT_sac

    	DT_sac = (abs(saclay%a)**2 - abs(saclay%b)**2 + abs(saclay%c)**2 &
                 - abs(saclay%d)**2 + abs(saclay%e)**2)/2.0_NER/I0000(mN, Saclay)

    end function DT_sac

    function AP_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: AP_sac

    	AP_sac = (real( conjg(saclay%a) * saclay%b - conjg(saclay%c) * saclay%d ) * cos(deg2rad(Saclay%x)/2.0_NER) &
    			 - aimag( conjg(saclay%b) * saclay%e) * sin(deg2rad(Saclay%x)/2.0_NER))/I0000(mN, Saclay)

    end function AP_sac

    function APT_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: APT_sac

    	APT_sac = (real( conjg(saclay%a) * saclay%c - conjg(saclay%b) * saclay%d ) * cos(deg2rad(Saclay%x)/2.0_NER) &
    			 - aimag( conjg(saclay%c) * saclay%e) * sin(deg2rad(Saclay%x)/2.0_NER))/I0000(mN, Saclay)

    end function APT_sac

    function AT_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: AT_sac

    	AT_sac = (real( - conjg(saclay%a) * saclay%c - conjg(saclay%b) * saclay%d ) * sin(deg2rad(Saclay%x)/2.0_NER) &
    			 - aimag( conjg(saclay%c) * saclay%e) * cos(deg2rad(Saclay%x)/2.0_NER))/I0000(mN, Saclay)

    end function AT_sac

    function AXX_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: AXX_sac

    	AXX_sac = (real( conjg(saclay%a) * saclay%d) * cos(deg2rad(Saclay%x)) + real( conjg(saclay%b) * saclay%c) &
    			  -aimag( conjg(saclay%d) * saclay%e) * sin(deg2rad(Saclay%x)))/I0000(mN, Saclay)

    end function AXX_sac

    function AZX_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: AZX_sac

    	AZX_sac = (real( conjg(saclay%a) * saclay%d) * sin(deg2rad(Saclay%x)) &
    			  +aimag( conjg(saclay%b) * saclay%e) * cos(deg2rad(Saclay%x)))/I0000(mN, Saclay)

    end function AZX_sac

    function AZZ_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: AZZ_sac

    	AZZ_sac = (-real( conjg(saclay%a) * saclay%d) * cos(deg2rad(Saclay%x)) + real( conjg(saclay%b) * saclay%c) &
    			  +aimag( conjg(saclay%d) * saclay%e) * sin(deg2rad(Saclay%x)))/I0000(mN, Saclay)

    end function AZZ_sac

    function D0SK_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: D0SK_sac

    	D0SK_sac = (real( conjg(saclay%a) * saclay%b + conjg(saclay%c) * saclay%d) * sin(deg2rad(Saclay%x)/2.0_NER) &
    				+ aimag( conjg(saclay%b) * saclay%e) * cos(deg2rad(Saclay%x)/2.0_NER))/I0000(mN, Saclay)

    end function D0SK_sac

    function NNKK_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: NNKK_sac

    	NNKK_sac = (- real( conjg(saclay%d) * saclay%e) * cos(deg2rad(Saclay%x)) &
                    - aimag( conjg(Saclay%a) * saclay%d) * sin(deg2rad(Saclay%x)))/I0000(mN, Saclay)

    end function NNKK_sac

    function NNSK_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: NNSK_sac

    	NNSK_sac = (- real( conjg(saclay%d) * saclay%e) * sin(deg2rad(Saclay%x)) +aimag( conjg(saclay%a) * saclay%d) &
    				* cos(deg2rad(Saclay%x)) + aimag( conjg(saclay%b) * saclay%c))/I0000(mN, Saclay)

    end function NNSK_sac

    function NSKN_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: NSKN_sac

    	NSKN_sac = real( conjg(saclay%c) * saclay%e) * sin(deg2rad(Saclay%x)/2.0_NER) &
    				+aimag( conjg(saclay%a) * saclay%c - conjg(saclay%b) * saclay%d) * cos(deg2rad(Saclay%x)/2.0_NER)

    end function NSKN_sac

    function NSSN_sac(mN, Saclay)

        real(NER),intent(in)            :: mN
        type(struct_Saclay),intent(in)  :: Saclay
        real(NER)                       :: NSSN_sac

        NSSN_sac = -aimag( conjg(saclay%a) * saclay%c) * sin(deg2rad(Saclay%x)/2.0_NER) &
                    - aimag( conjg(saclay%b) * saclay%d) * sin(-deg2rad(Saclay%x)/2.0_NER)

    end function NSSN_sac


    function R_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: R_sac

    	R_sac = real( conjg(saclay%a) * saclay%b + conjg(saclay%c) * saclay%d) * cos(deg2rad(Saclay%x)/2.0_NER) &
    			-aimag( conjg(saclay%b) * saclay%e) * sin(deg2rad(Saclay%x)/2.0_NER)

    end function R_sac

    function RP_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: RP_sac

    	RP_sac = real( conjg(saclay%a) * saclay%b - conjg(saclay%c) * saclay%d) * sin(deg2rad(Saclay%x)/2.0_NER) &
    			+aimag( conjg(saclay%b) * saclay%e) * cos(deg2rad(Saclay%x)/2.0_NER)

    end function RP_sac

    function RPT_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: RPT_sac

    	RPT_sac = real( -conjg(saclay%a) * saclay%c + conjg(saclay%b) * saclay%d) * sin(deg2rad(Saclay%x)/2.0_NER) &
    			-aimag( conjg(saclay%c) * saclay%e) * cos(deg2rad(Saclay%x)/2.0_NER)

    end function RPT_sac

    function RT_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: RT_sac

    	RT_sac = real( -conjg(saclay%a) * saclay%c - conjg(saclay%b) * saclay%d) * cos(deg2rad(Saclay%x)/2.0_NER) &
    			-aimag( conjg(saclay%c) * saclay%e) * sin(deg2rad(Saclay%x)/2.0_NER)

    end function RT_sac

    function SGT_sac(mN, Saclay)

    	real(NER),intent(in)			:: mN
    	type(struct_Saclay),intent(in)	:: Saclay
    	real(NER)						:: SGT_sac
        real(NER)                       :: k

        k = Tlab2kcm(saclay%Elab, mN)
        SGT_sac = MeV_to_barn * 2.0_NER * PI_NE * aimag(saclay%a+saclay%b) / k

    end function SGT_sac

    function I0000(mN, Saclay)

        real(NER),intent(in)            :: mN
        type(struct_Saclay), intent(in)   :: Saclay
        real(NER)                       :: I0000

        I0000 = (abs(saclay%a)**2 + abs(saclay%b)**2 + abs(saclay%c)**2 + abs(saclay%d)**2 + abs(saclay%e)**2) / 2.0_NER

    end function I0000

    function Cnn00(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cnn00

        Cnn00 = (abs(saclay%a)**2-abs(saclay%b)**2-abs(saclay%c)**2+abs(saclay%d)**2+abs(saclay%e)**2)/2.0_NER/I0000(mN,Saclay)

    end function Cnn00

    function Dn0n0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Dn0n0

        Dn0n0 = (abs(saclay%a)**2+abs(saclay%b)**2-abs(saclay%c)**2-abs(saclay%d)**2+abs(saclay%e)**2)/2.0_NER/I0000(mN,Saclay)

    end function Dn0n0

    function K0nn0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: K0nn0

        K0nn0 = (abs(saclay%a)**2-abs(saclay%b)**2+abs(saclay%c)**2-abs(saclay%d)**2+abs(saclay%e)**2)/2.0_NER/I0000(mN,Saclay)

    end function K0nn0

    function Cllll(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cllll

        Cllll = (abs(saclay%a)**2-abs(saclay%b)**2+abs(saclay%c)**2-abs(saclay%d)**2+abs(saclay%e)**2)/2.0_NER/I0000(mN,Saclay)

    end function Cllll

    function Pn000(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Pn000

        Pn000 = real(conjg(saclay%a) * saclay%e)/I0000(mN,Saclay)

    end function Pn000

    function Clllm(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Clllm

        Clllm = aimag(conjg(saclay%a)*saclay%e)/I0000(mN,Saclay)

    end function Clllm

    function Clnl0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Clnl0

        Clnl0 = real(conjg(saclay%b)*saclay%e)/I0000(mN,Saclay)

    end function Clnl0

    function Dl0m0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Dl0m0

        Dl0m0 = aimag(conjg(saclay%b)*saclay%e)/I0000(mN,Saclay)

    end function Dl0m0

    function Cnll0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cnll0

        Cnll0 = real(conjg(saclay%c)*saclay%e)/I0000(mN,Saclay)

    end function Cnll0

    function K0lm0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: K0lm0

        K0lm0 = aimag(conjg(saclay%c)*saclay%e)/I0000(mN,Saclay)

    end function K0lm0

    function Clln0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Clln0

        Clln0 = - real(conjg(saclay%d)*saclay%e)/I0000(mN,Saclay)

    end function Clln0

    function Clm00(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Clm00

        Clm00 = aimag(conjg(saclay%d)*saclay%e)/I0000(mN,Saclay)

    end function Clm00

    function Dm0m0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Dm0m0

        Dm0m0= real(conjg(saclay%a)*saclay%b + conjg(saclay%c)*saclay%d)/I0000(mN,Saclay)

    end function Dm0m0

    function Cmnl0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cmnl0

        Cmnl0= aimag(conjg(saclay%a)*saclay%b + conjg(saclay%c)*saclay%d)/I0000(mN,Saclay)

    end function Cmnl0

    function Dl0l0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Dl0l0

        Dl0l0= real(conjg(saclay%a)*saclay%b - conjg(saclay%c)*saclay%d)/I0000(mN,Saclay)

    end function Dl0l0

    function Clnm0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Clnm0

        Clnm0= - aimag(conjg(saclay%a)*saclay%b - conjg(saclay%c)*saclay%d)/I0000(mN,Saclay)

    end function Clnm0

    function K0mm0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: K0mm0

        K0mm0= real(conjg(saclay%a)*saclay%c + conjg(saclay%b)*saclay%d)/I0000(mN,Saclay)

    end function K0mm0

    function Cnml0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cnml0

        Cnml0= - aimag(conjg(saclay%a)*saclay%c + conjg(saclay%b)*saclay%d)/I0000(mN,Saclay)

    end function Cnml0

    function K0ll0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: K0ll0

        K0ll0= real(conjg(saclay%a)*saclay%c - conjg(saclay%b)*saclay%d)/I0000(mN,Saclay)

    end function K0ll0

    function Cnlm0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cnlm0

        Cnlm0= - aimag(conjg(saclay%a)*saclay%c - conjg(saclay%b)*saclay%d)/I0000(mN,Saclay)

    end function Cnlm0

    function Cmm00(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cmm00

        Cmm00= real(conjg(saclay%a)*saclay%d + conjg(saclay%b)*saclay%c)/I0000(mN,Saclay)

    end function Cmm00

    function Clmn0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Clmn0

        Clmn0= - aimag(conjg(saclay%a)*saclay%c + conjg(saclay%b)*saclay%d)/I0000(mN,Saclay)

    end function Clmn0

    function Cll00(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cll00

        Cll00= - real(conjg(saclay%a)*saclay%d - conjg(saclay%b)*saclay%c)/I0000(mN,Saclay)

    end function Cll00

    function Cmln0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Cmln0

        Cmln0= - aimag(conjg(saclay%a)*saclay%d - conjg(saclay%b)*saclay%c)/I0000(mN,Saclay)

    end function Cmln0




    function Dsp0s0(mN, Saclay)

        real(NER),intent(in)    :: mN
        type(struct_Saclay), intent(in) :: Saclay
        real(NER)                       :: Dsp0s0

        Dsp0s0 = (cos(deg2rad(Saclay%x)/2.0_NER) * (Real(conjg(saclay%a)*saclay%b) + real(conjg(saclay%c)*saclay%d)) - &
                & sin(deg2rad(Saclay%x)/2.0_NER) * aimag(conjg(saclay%b)*saclay%e))/I0000(mN,Saclay)

    end function Dsp0s0


end module mod_obsv