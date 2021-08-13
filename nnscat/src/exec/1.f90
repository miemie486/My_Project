!main program to call fit subroutine

program main 

    use mod_fit  

    implicit none 

    integer,parameter         :: Nmesh=50, uptoQn=0 
    real(NER)       :: Mambda=800.0_NER

    real(NER)       :: minchisqr

    real			:: t1,t2 ! calculate run time.


    real(NER) :: unknpara(4)=(/0.63645639631157650E+002_NER , -0.77966138474560367E-001_NER, &
            	-4.218958142353386E-008_NER, 3.155774938761064E-003_NER /)



    call cpu_time(t1)  

    call fitdata(Nmesh,uptoQn,Mambda,unknpara,minchisqr)

    call cpu_time(t2)
    write(*,*) t2-t1
     
end program main
