! Bingwei Long 05/23/2012
! solve for Delta & y for 1S0
! C0(E) = y/(E+Delta)

program slvDy

    use util_io
    use util_gadgets
    use nneft_type
    use nneft_phyconst
    use eft_potspwd
    implicit none

    integer, parameter  :: funt1 = 23, funt2 = 24
    integer             :: ii, numLmbd, flag
    real(NER)           :: ecm_1, ecm_2, d, y, PC_CgA, tlab1 = 25.0, tlab2 = 50.0, k

    real(NER), dimension(50) :: C0_1, C0_2, Lmbd_1, Lmbd_2
    character(len=80) :: fname1, fname2

    fname1 = 'slvc0_1s0_C0_5MeV.out'
    fname2 = 'slvc0_1s0_C0_15MeV.out'
    tlab1  = 5.0
    tlab2  = 15.0
    k      = tlab2k(tlab1)
    ecm_1  = k*k/PC_mN
    k      = tlab2k(tlab2)
    ecm_2  = k*k/PC_mN
    PC_CgA = PC_gA*PC_gA*PC_mN/(8.0*PI_NE*PI_NE*PC_fpi*PC_fpi)
    print *, PC_CgA
!~     CgA = PC_gA*PC_gA*PC_mN/(8.0*PI_NE*PI_NE*PC_fpi*PC_fpi)
    call read2cols_io(funt1, fname1, C0_1, Lmbd_1, numLmbd, flag)
    call read2cols_io(funt2, fname2, C0_2, Lmbd_2, numLmbd, flag)
    C0_1 = C0_1 + PC_CgA
    C0_2 = C0_2 + PC_CgA
    do ii = 1, numLmbd
        if (abs((Lmbd_1(ii)-Lmbd_2(ii))/Lmbd_1(ii)) .gt. 0.01) print *, Lmbd_1(ii), Lmbd_2(ii)
        d = -(C0_1(ii)*ecm_1 - C0_2(ii)*ecm_2)/(C0_1(ii) - C0_2(ii))
        y = C0_1(ii)*(ecm_1 + d)
        print *, d, y, Lmbd_1(ii)
    end do

end program slvDy
