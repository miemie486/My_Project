

module rel_ope_pwd


  use nneft_type,       only: NER, NEC, PI_NE

  use nneft_phyconst,   only: PC_mN, PC_gA, PC_fpi, PC_mpi
  use ope_pwd

  implicit none



contains



 function WT_rel(q0sqr,qsqr)
    real(NER) ::  WT_rel
    real(NER),intent(in) :: q0sqr,qsqr
    real(NER) :: factor,PC_mpi_sqr

    factor=(PC_gA/PC_fpi)*(PC_gA/PC_fpi)*0.25_NER
    PC_mpi_sqr = PC_mpi * PC_mpi

    WT_rel=-factor*q0sqr/((qsqr +PC_mpi_sqr)*(qsqr +PC_mpi_sqr))
 end function WT_rel


 function WR_rel(q0,qsqr)
    real(NER) ::  WR_rel,z
    real(NER),intent(in) :: q0,qsqr
    real(NER) :: factor,PC_mpi_sqr,m_N
    factor=(PC_gA/PC_fpi)*(PC_gA/PC_fpi)*0.25_NER
    PC_mpi_sqr = PC_mpi * PC_mpi
    m_N=PC_mN

    WR_rel=factor*q0/m_N/(qsqr+PC_mpi_sqr )
 end function WR_rel



 function rel_OPE_j0j(j,k,p)
    real(NER) :: rel_OPE_j0j
    real(NER),intent(in) :: k,p
    integer,intent(in) :: j
    real(NER) :: relation_factor
    integer :: ii
    real(NER) :: regsum,z,q0sqr,qsqr
    relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
    q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))


    if(j<0) then
      rel_OPE_j0j=0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if


      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         z=mshlg_pwo(ii)
         qsqr=k*k+p*p-2.0_NER*k*p*z
          regsum=regsum+2.0_NER*PI_NE*WT_rel(q0sqr,qsqr)*(-(k*k+p*p)+2.0_NER*k*p*z)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)
       end do

   if(mod(j,2)==1) then

      regsum=-3.0_NER*regsum

    end if

    rel_OPE_j0j=relation_factor*regsum
  end function rel_OPE_j0j



 function rel_OPE_j1j(j,k,p)
    real(NER) :: rel_OPE_j1j
    real(NER),intent(in) :: k,p
    integer,intent(in) :: j
    real(NER) :: relation_factor
    integer :: ii
    real(NER) :: regsum,z,q0sqr,qsqr
    relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
    q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))


    if(j<1) then
      rel_OPE_j1j=0.0_NER
      return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if


     regsum=0.0_NER
     do ii=1,Nlg_pwo,1
          z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
        regsum=regsum+2.0_NER*PI_NE*WT_rel(q0sqr,qsqr)*((k*k+p*p)*lgndr_table_pwo(j, ii)&
        &-2.0_NER/(2.0_NER*j+1.0_NER)*j*k*p*lgndr_table_pwo(j+1, ii)&
       &-2.0_NER/(2.0_NER*j+1.0_NER)*(j+1.0_NER)*k*p*lgndr_table_pwo(j-1, ii))*wghtslg_pwo(ii)
      end do

    if(mod(j,2)==0) then

      regsum=-3.0_NER*regsum
    end if

     rel_OPE_j1j=relation_factor*regsum
 end function rel_OPE_j1j



 function rel_OPE_jpp(j,k,p)
   real(NER) :: rel_OPE_jpp
     real(NER),intent(in) :: k,p
      integer,intent(in) :: j
   real(NER) :: relation_factor
        integer :: ii
     real(NER) :: regsum,z,q0sqr,qsqr
      relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
   q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))


   if(j<0) then
        rel_OPE_jpp=0.0_NER
       return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if


       regsum=0.0_NER
         do ii=1,Nlg_pwo,1
       z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
      regsum=regsum+2.0_NER*PI_NE*WT_rel(q0sqr,qsqr)*&
       &1.0_NER/(2.0_NER*j+1.0_NER)*(-(k*k+p*p)*lgndr_table_pwo(j+1, ii)&
     &+2.0_NER*k*p*lgndr_table_pwo(j, ii) )*wghtslg_pwo(ii)
      end do

    if(mod(j,2)==1) then

      regsum=-3.0_NER*regsum

    end if

    rel_OPE_jpp=relation_factor*regsum
 end function rel_OPE_jpp



 function rel_OPE_jmm(j,k,p)
    real(NER) :: rel_OPE_jmm
    real(NER),intent(in) :: k,p
       integer,intent(in) :: j
     real(NER) :: relation_factor
   integer :: ii
       real(NER) :: regsum,z,q0sqr,qsqr
     relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
     q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))


    if(j<1) then
       rel_OPE_jmm=0.0_NER
      return
   end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if


      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
      z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
         regsum=regsum+2.0_NER*PI_NE*WT_rel(q0sqr,qsqr)*&
      &1.0_NER/(2.0_NER*j+1.0_NER)*((k*k+p*p)*lgndr_table_pwo(j-1, ii)&
     &-2.0_NER*k*p*lgndr_table_pwo(j, ii) )*wghtslg_pwo(ii)
      end do

    if(mod(j,2)==1) then

        regsum=-3.0_NER*regsum

    end if

      rel_OPE_jmm=relation_factor*regsum
 end function rel_OPE_jmm



 function rel_OPE_jpm(j,k,p)
    real(NER) :: rel_OPE_jpm
   real(NER),intent(in) :: k,p
       integer,intent(in) :: j
        real(NER) :: relation_factor
       integer :: ii
       real(NER) :: regsum,z,q0sqr,qsqr,q0
       relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
       q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))
        q0=(p*p-k*k)/(2.0_NER*PC_mN)


    if(j<1) then
       rel_OPE_jpm=0.0_NER
         return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if


       regsum=0.0_NER
        do ii=1,Nlg_pwo,1
       z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
         regsum=regsum+(-2.0_NER*PI_NE)*&
       &2.0_NER/(2.0_NER*j+1.0_NER)*sqrt(j*(j+1.0_NER))*(WT_rel(q0sqr,qsqr)*(p*p*lgndr_table_pwo(j+1, ii)&
       &+k*k*lgndr_table_pwo(j-1, ii)-2.0_NER*k*p*lgndr_table_pwo(j, ii) )&
       &- WR_rel(q0,qsqr)*k*p*(lgndr_table_pwo(j+1, ii)-lgndr_table_pwo(j-1, ii)))*wghtslg_pwo(ii)
       end do


    if(mod(j,2)==1) then

        regsum=-3.0_NER*regsum

    end if

       rel_OPE_jpm=relation_factor*regsum
  end function rel_OPE_jpm



  function rel_OPE_jmp(j,k,p)
   real(NER) :: rel_OPE_jmp
       real(NER),intent(in) :: k,p
       integer,intent(in) :: j
     real(NER) :: relation_factor
      integer :: ii
         real(NER) :: regsum,z,q0sqr,qsqr,q0
       relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
       q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))
    q0=(p*p-k*k)/(2.0_NER*PC_mN)

    
    if(j<1) then
       rel_OPE_jmp=0.0_NER
         return
    end if
    if (.not. mshlg_inited_pwo) then
      call initmshlg_pwo()
    end if

    
    regsum=0.0_NER
      do ii=1,Nlg_pwo,1
        z=mshlg_pwo(ii)
         qsqr=k*k+p*p-2.0_NER*k*p*z
        regsum=regsum+(-2.0_NER*PI_NE)*&
      &2.0_NER/(2.0_NER*j+1.0_NER)*sqrt(j*(j+1.0_NER))*(WT_rel(q0sqr,qsqr)*(k*k*lgndr_table_pwo(j+1, ii)&
     &+p*p*lgndr_table_pwo(j-1, ii)-2.0_NER*k*p*lgndr_table_pwo(j, ii) )&
      &+ WR_rel(q0,qsqr)*k*p*(lgndr_table_pwo(j+1, ii)-lgndr_table_pwo(j-1, ii)))*wghtslg_pwo(ii)
      end do

   if(mod(j,2)==1) then

       regsum=-3.0_NER*regsum


    end if

       rel_OPE_jmp=relation_factor*regsum
  end function rel_OPE_jmp



end module rel_ope_pwd
