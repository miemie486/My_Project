!wenchao shi    2019/10/24

program deuteron_form
	use formfactor
	use util_cheb
	implicit none
! 	real(NER) :: qq
! 	integer :: ll
! 	integer,parameter  ::  unit_qlistin=104
! 	integer,parameter  ::  unit_G_Cout=105

	
! 	open(unit=unit_qlistin,file="deuteron/qlist.in")  
! 	open(unit=unit_G_Cout,file="deuteron/G_C.out")
!   do ll=1,40
! 	 read(unit_qlistin,*) qq
! 	 write(unit_G_Cout,*) qq,G_C(q)
!   end do



    real(NER)                  :: a=-1.0, b=1.0,lamb=0.1
    integer                    :: nl=3
    real(NER)    :: gl(3)
    integer :: ww

    open(unit=17,file="deuteron/qlist.in")
    open(unit=18,file="deuteron/G_C_LONLO.out")

do ww=1,43
   read(17,*) q
   call cheb_taylor_real(a, b, nl, gl, G_C_)
   write(18,*) q,gl

end do
 close(17)
 close(18)

 contains
 function G_C_(lamb)

    real(NER),intent(in)    :: lamb
    real(NER)    :: G_C_,G_C_lamb

    call getG_C_lamb(lamb,G_C_lamb)
    G_C_=G_C_lamb

  end function G_C_


end program deuteron_form
