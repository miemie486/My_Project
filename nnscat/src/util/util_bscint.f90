! Bingwei Long 12/05/2010
! basic integrals

! Obsolete: 08/10/2018
module util_bscint

    use nneft_type
    implicit none

    contains

    ! If TanMesh is used, beta_k_bscint must always return 0.0
    function beta_k_bscint(k, Lambda)

        complex(NEC)            :: beta_k_bscint
        real(NER), intent(in)   :: k, Lambda

        ! beta_k_bscint = 0.5_NER * (log((Lambda+k) / (Lambda-k)) - IMGUNT_NE * PI_NE)
        beta_k_bscint = - 0.5_NER * IMGUNT_NE * PI_NE

    end function beta_k_bscint

    ! Calculate principle integral using sharp cutoff
    ! \mathcal P \int_0^\Lambda dx f(x)/(y^2 - x^2), where Lambda > y > 0
    ! f is given as an (1:N+1) array, fcn(ii), evaluated already over the mesh points, ie,
    ! fcn(ii) = f(meshpts(ii)) and fcn(N+1) = f(y)

    function prncpl_shrp(y, fcn, Lambda, N, meshpts, wghts)

        complex(NEC)                :: prncpl_shrp
        integer, intent(in)         :: N
        real(NER), intent(in)       :: y, Lambda, meshpts(:), wghts(:)
        ! real(NER), dimension(1:N), intent(in)       :: meshpts, wghts
        ! complex(NEC), dimension(1:N+1), intent(in)  :: fcn
        complex(NEC), intent(in)    :: fcn(:)

        integer         :: ii
        real(NER)       :: ysqr, xx
        complex(NEC)    :: regsum, fy

        regsum  = (0.0_NER, 0.0_NER)
        fy      = fcn(N+1)

        do ii = 1, N
            xx      = meshpts(ii)
            regsum  = regsum + (fcn(ii) - fy) / (y*y - xx*xx) * wghts(ii)
        end do

        prncpl_shrp = regsum + 0.5_NER * fy / y * log((Lambda+y) / (Lambda-y))

    end function prncpl_shrp

    ! Calculate singular integral with +i\epsilon prescription
    ! \int_0^\Lambda dx f(x)/(y^2 - x^2 + i\epsilon), where Lambda > y > 0

    function plusepsln_shrp(y, fcn, Lambda, N, meshpts, wghts)

        complex(NEC)                :: plusepsln_shrp
        integer, intent(in)         :: N
        real(NER), intent(in)       :: y, Lambda, meshpts(:), wghts(:)
        ! real(NER), dimension(1:N), intent(in)       :: meshpts, wghts
        complex(NEC), intent(in)    :: fcn(:)
        ! complex(NEC), dimension(1:N+1), intent(in)  :: fcn

        plusepsln_shrp = prncpl_shrp(y, fcn, Lambda, N, meshpts, wghts) - 0.5_NER * IMGUNT_NE * PI_NE / y * fcn(N+1)

    end function plusepsln_shrp


! same as plusepsln_shrp, but the denominator (y^2 - x^2) is input by "den". the result is put out by "qdr".

    subroutine get_IntPlusEPS(y, fcn, den, Lambda, N, meshpts, wghts, qdr)

        integer, intent(in)         :: N
        real(NER), intent(in)       :: y, Lambda, den(:), meshpts(:), wghts(:)
        ! real(NER), dimension(1:N), intent(in)       :: den, meshpts, wghts
        complex(NEC), intent(in)    :: fcn(:)
        ! complex(NEC), dimension(1:N+1), intent(in)  :: fcn

        complex(NEC), intent(out) :: qdr

        integer         :: ii
        complex(NEC)    :: regsum, fy

        regsum  = (0.0_NER, 0.0_NER)
        fy      = fcn(N+1)

        do ii = 1, N
            regsum = regsum + (fcn(ii) - fy) / den(ii) * wghts(ii)
        end do

        qdr = regsum + 0.5_NER*fy/y*(log((Lambda+y)/(Lambda-y)) - IMGUNT_NE*PI_NE)

    end subroutine get_IntPlusEPS


end module util_bscint
