
program test_lgndrtbl

    use ope_pwd
    implicit none

    integer     :: ii, N
    real(NER)   :: z


    N = Nlg_pwo
    call initmshlg_pwo()
    do ii = 1, N
        z = mshlg_pwo(ii)
!~         print *, 'z = ', z
!~         print *, lgndr_table_pwo(0, ii), '1.0'
!~         print *, lgndr_table_pwo(1, ii), z
!~         print *, lgndr_table_pwo(2, ii), 0.5*(3.0*z*z - 1.0)
        if (abs(lgndr_table_pwo(1, ii) - z)/z .gt. 1.0E-3) print *, z, lgndr_table_pwo(1, ii), z
        if (abs(lgndr_table_pwo(2, ii) - 0.5*(3.0*z*z - 1.0))/lgndr_table_pwo(2, ii) .gt. 1.0E-3) &
            & print *, lgndr_table_pwo(2, ii), 0.5*(3.0*z*z - 1.0)
    end do

end program test_lgndrtbl
