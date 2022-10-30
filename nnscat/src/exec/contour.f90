program contour

    use mod_cmplx_quad
    implicit none

    integer, parameter  :: maxpoints = 100000, maxsegs = 10
    integer :: ii, numseg, segN(1:maxsegs), total_N
    complex(NEC)    :: tmpsum, l(1:maxpoints), w(1:maxpoints), cntr_pts(1:maxsegs+1)
    real(NER)       :: delta

    ! numseg = 1
    ! segN(1) = 20000
    ! cntr_pts(1) = (0.0_NER, 0.0_NER)
    ! cntr_pts(2) = (2.0_NER, 0.0_NER)

    numseg = 3
    cntr_pts(1) = (0.0_NER, 0.0_NER)
    segN(1) = 500
    delta = 1.0E-6_NER
    cntr_pts(2) = cmplx(1.0_NER - delta, 0.0_NER, NEC)
    segN(2) = 10
    cntr_pts(3) = cmplx(1.0_NER + delta, 0.0_NER, NEC)
    segN(3) = 500
    cntr_pts(4) = (2.0_NER, 0.0_NER)

    call build_zigzag_mesh(numseg, segN, cntr_pts, l, w)

    ! total_N = sum(segN(1:numseg))
    total_N = 0
    do ii = 1, numseg
        total_N = total_N + segN(ii)
    end do
    tmpsum = 0.0_NER
    do ii = 1, total_N
        tmpsum = tmpsum + f(l(ii))*w(ii)
    end do
    print *, tmpsum

contains

    complex(NEC) function f(t)
        complex(NEC), intent(in)    :: t

        f = 1.0_NER/(sqrt(t - (1.0_NER, 0.0_NER)))
        ! f = 1.0/(t - (1.0_NER, 0.0_NER))

    end function f

end program contour
