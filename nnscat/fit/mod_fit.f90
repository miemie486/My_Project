module mod_fit

use mod_obsv
use drv_fitpwa


implicit none

real(NER)       :: mN = 938.0_NER
external futil

type :: elab_chain
    real(NER) :: elab
    type(elab_chain),pointer :: next
end type elab_chain

type :: data_chain
    character(len=6)    :: Ref
    character(len=4)    :: obsName
    real(NER)           :: acm, elab, val, error
    type(data_chain),pointer    :: next
end type data_chain



!there are 640 different elab in 0-300mev data.
real(NER),allocatable       :: elist(:)




type(struct_phase_allwaves)     :: phase_sac
type(struct_Saclay)             :: saclay







contains


subroutine fitdata(Nmesh,uptoQn,Mambda,unknpara,minchisqr)

    real(NER) :: argmn(1:5), parup, parlow, parerr
    integer, intent(in)         :: Nmesh, uptoQn
    real(NER), intent(in)       :: Mambda
    real(NER), intent(inout)    :: unknpara(:)
    real(NER), intent(out)      :: minchisqr
    integer,parameter           :: NPAR=4
    integer                     :: numk,numdata
    real(NER),allocatable       :: klist(:)
    character(len=6)            :: paraname(NPAR)
    integer                     :: nn
    real(NER)                   :: para(4)
    integer                     :: errflag
    real(NER)                   :: limit = 0.0_NER
    real(NER)                   :: limit1,limit2
    type(struct_expt_data),allocatable      :: expt_data(:)
    type(struct_phase_allwaves),allocatable :: phase(:)

    real(NER)                   :: FVAL
    real(NER)                   :: GRAD(1:NPAR)
    integer                     :: IFLAG
    real(NER)                   :: XVAL(NPAR)

    character(len=20)           ::data_filename="all_data.dat"
!     external FUTIL
    paraname(1) = "delta"
    paraname(2) = "y"
    paraname(3) = "c0_3p0"
    paraname(4) = "c0_3s1"

    limit1 = 0.     ! limit1 is lower limit of data to fit, limit2 is upper
    limit2 = 300.


    para=(/ 0.59100655067130575E+002_NER , -0.72482410781848861E-001_NER , &
            -4.218958142353386E-008_NER, 3.155774938761064E-003_NER /)

    xval = unknpara

    call get_numdata(numdata,limit1,limit2)
    allocate(expt_data(numdata),stat=errflag)
    call get_data(expt_data,numdata,limit1,limit2)
    call get_numk(numk,expt_data,numdata)
    allocate(elist(numk),stat=errflag)
    call get_elist(elist,numk,expt_data,numdata)
    allocate(phase(numk))
    allocate(klist(numk))

    klist=sqrt(elist*mN*0.5_NER)

    write(*,*) numk,numdata
    write(*,*) elist(1),elist(2),klist(1),klist(2)


     call chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)



!     call MNINIT (5, 6, 7)
!     call MNSETI ('Fit')
!     argmn(1) = 1.0_NER
!     call MNEXCM (chisqr, 'SET PRINTOUT', argmn, 1, errflag, FUTIL)

!     do nn= 1,NPAR
!         CALL MNPARM(nn,paraname(nn),para(nn),abs(para(nn))*0.001,limit,limit,errflag)

!     end do



!     argmn(1) = 2.
!     call MNEXCM (chisqr, 'SET STRATEGY', argmn, 1, errflag, FUTIL)

!     argmn(1) = 10000.
!     argmn(2) = 0.001
!     call MNEXCM(chisqr, 'migrad',argmn,2,errflag,futil)












contains
    subroutine chisqr(NPAR, GRAD, FVAL, XVAL, IFLAG, FUTIL)

        integer, intent(in)     :: NPAR, IFLAG
        real(NER), intent(in)   :: XVAL(1:NPAR)
        real(NER), intent(out)  :: FVAL
        real(NER), intent(out)  :: GRAD(1:NPAR)

        integer                 :: chnid
        real(NER)               :: paralst(NPAR)
        real(NER)               :: phaselstun(numk,uptoQn+1)
        type(triplet_phs_cpld)  :: phaselstcp(numk,uptoQn+1)
        type(struct_phase_allwaves) :: phase(1:numk)
        real(NER)               :: acm,edata,xdata,errdata
        character(len=4)        :: obsname
!         integer,parameter       :: numdata=35
        real(NER)               :: chisquare(numdata)
        integer                 :: ii,kk
        integer                 :: istat

        real(NER)               :: para1s0(2),para3s1(1),para3p0(1)
        real(NER)               :: phaselst3s1(numk:3)
        real :: t1,t2
        external FUTIL


        phase%Jmax0 = 0
        phase%Jmax1 = 0
        phase%Jmaxcp = 0

!1s0
        para1s0=xval(1:2)
        chnid = 100

        call cpu_time(t1)
        call Get1s0dib_fitpwa(Nmesh, Mambda, chnid, uptoQn, para1s0,           &
         & 2, klist, numk, phaselstun)

        do ii=1,numk
            phase(ii)%unps0(0)=phaselstun(ii,uptoQN+1)
        end do

!3p0
        para3p0=xval(3)
        chnid = 310

        call GetNPWOC0Sngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, para3p0,       &
        & 1, klist, numk, phaselstun)

        do ii=1,numk
            phase(ii)%unps1(0)=phaselstun(ii,uptoQN+1)
        end do



!3s1
        para3s1=xval(4)
        chnid = 301



        call GetNPCpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, para3s1,           &
        & 1, klist, numk, phaselstcp)

        do ii=1,numk
            phase(ii)%cpldps(1,1)=phaselstcp(ii,uptoQn+1)%d1
            phase(ii)%cpldps(1,2)=phaselstcp(ii,uptoQn+1)%d2
            phase(ii)%cpldps(1,3)=phaselstcp(ii,uptoQn+1)%e
        end do

        call cpu_time(t2)
        write(*,*) "calculate phases time ", t2-t1









! !       1S0
!         chnid = 100
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps0(0)=phaselstun(ii,uptoQN+1)

!         end do

! !       1P1
!         chnid = 111
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps0(1)=sum(phaselstun(ii,1:uptoQN+1))
!         end do

! !       1D2
!         chnid = 122
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps0(2)=sum(phaselstun(ii,1:uptoQN+1))
!         end do

! !       1F3
!         chnid = 133
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps0(3)=sum(phaselstun(ii,1:uptoQN+1))
!         end do

! !       3P0
!         chnid = 310
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps1(0)=sum(phaselstun(ii,1:uptoQN+1))
!         end do

! !       3P1
!         chnid = 111
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps1(1)=sum(phaselstun(ii,1:uptoQN+1))
!         end do

! !       3D2
!         chnid = 322
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps1(2)=sum(phaselstun(ii,1:uptoQN+1))
!         end do

! !       3F3
!         chnid = 333
!         paralst()=XVAL()
!         num_paras=

!         call GetSngl_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstun)

!         do ii=1,numk
!             phase(ii)%unps1(3)=sum(phaselstun(ii,1:uptoQN+1))
!         end do

! !       3S1-3D1
!         chnid = 301
!         paralst()=XVAL()
!         num_paras=

!         call GetCpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstcp)

!         do ii=1,numk
!             phase(ii)%cpldps(1,1)=sum(phaselstcp(ii,1:uptoQn+1)%d1)
!             phase(ii)%cpldps(1,2)=sum(phaselstcp(ii,1:uptoQn+1)%d2)
!             phase(ii)%cpldps(1,3)=sum(phaselstcp(ii,1:uptoQn+1)%e)
!         end do

! !       3P2-3F2
!         chnid = 312
!         paralst()=XVAL()
!         num_paras=

!         call GetCpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstcp)

!         do ii=1,numk
!             phase(ii)%cpldps(2,1)=sum(phaselstcp(ii,1:uptoQn+1)%d1)
!             phase(ii)%cpldps(2,2)=sum(phaselstcp(ii,1:uptoQn+1)%d2)
!             phase(ii)%cpldps(2,3)=sum(phaselstcp(ii,1:uptoQn+1)%e)
!         end do

! !       3D3-3G3
!         chnid = 323
!         paralst()=XVAL()
!         num_paras=

!         call GetCpld_fitpwa(Nmesh, Mambda, chnid, uptoQn, paralst,           &
!         & num_paras, klist, numk, phaselstcp)

!         do ii=1,numk
!             phase(ii)%cpldps(3,1)=sum(phaselstcp(ii,1:uptoQn+1)%d1)
!             phase(ii)%cpldps(3,2)=sum(phaselstcp(ii,1:uptoQn+1)%d2)
!             phase(ii)%cpldps(3,3)=sum(phaselstcp(ii,1:uptoQn+1)%e)
!         end do

!        open(unit=101,file="300mev.txt",status='old')

            phase_sac%jmax0=0
            phase_sac%jmax1=0
            phase_sac%jmaxcp=1

        call cpu_time(t1)
        do ii=1,numdata

!         read(330,*) expt_data(nn)%obsname,expt_data(nn)%energy,expt_data(nn)%angle,&
!                     expt_data(nn)%val,expt_data(nn)%error,expt_data(nn)%ref

            kk=BINARY_SEARCH(elist,numk,expt_data(ii)%elab)

!             write(*,*) "kk= ", kk


            phase_sac%unps0(0) = phase(kk)%unps0(0)
            phase_sac%unps1(0) = phase(kk)%unps1(0)
            phase_sac%cpldps(1,1)=phase(kk)%cpldps(1,1)
            phase_sac%cpldps(1,2)=phase(kk)%cpldps(1,2)
            phase_sac%cpldps(1,3)=phase(kk)%cpldps(1,3)

            phase_sac%Elab=expt_data(ii)%elab

            !write phaseshift on screen
!             write(6,"(f6.2,5f7.2)")phase_sac%elab,phase_sac%unps0(0),phase_sac%unps1(0), &
!             phase_sac%cpldps(1,1),phase_sac%cpldps(1,2),phase_sac%cpldps(1,3)


            call Saclay_para_obs(Saclay, acm, mN, phase_sac)




            select case (trim(expt_data(ii)%obsname))

            case("DSG")
                chisquare(ii)=(DSG_sac(mN,Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
!                 write(*,*) expt_data(ii)%val,dsg_sac(mn,saclay),(DSG_sac(mN,Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2

            case("P","A")
                chisquare(ii)=(P_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("AYY")
                chisquare(ii)=(AYY_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("D")
                chisquare(ii)=(D_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("DT")
                chisquare(ii)=(DT_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("AT")
                chisquare(ii)=(AT_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("AZZ")
                chisquare(ii)=(AZZ_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("D0SK")
                chisquare(ii)=(D0SK_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("NNKK")
                chisquare(ii)=(NNKK_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("NSKN")
                chisquare(ii)=(NSKN_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("NSSN")
                chisquare(ii)=(NSSN_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("R")
                chisquare(ii)=(R_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("RPT")
                chisquare(ii)=(RPT_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("RT")
                chisquare(ii)=(RT_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case("SGT")
                chisquare(ii)=(SGT_sac(mN, Saclay)-expt_data(ii)%val)**2/expt_data(ii)%error**2
            case default
                write(6,*) "wrong obsname. data number is ",ii ,expt_data(ii)%obsname

            end select
        end do
        call cpu_time(t2)
        write(*,*) "fit time ",t2-t1
!         close(101)


        FVAL=sum(chisquare(1:numdata))/numdata
        write(*,*) "FVAL=",FVAL


        write(*,"(a7,f10.7,2x,f10.7,2x,f10.7,2x,f10.7)") "XVAL=",xval(1) ,xval(2),xval(3),xval(4)

    end subroutine chisqr

    function BINARY_SEARCH(A,N,KEY)

        integer   :: N,BINARY_SEARCH
        real(NER) :: KEY,A(N)
        integer   ::L,R,M

        L=1
        R=N
        M=(L+R)/2

        if ( (KEY < A(L)) .OR. (KEY > A(R)) ) then
            BINARY_SEARCH = 0
            return
        end if

        do while( L <= R )
            if ( ABS(KEY - A(M))<0.00001 ) then
                BINARY_SEARCH = M
                return
            else if ( KEY < A(M) ) then
                R=M-1
                M=(L+R)/2
            else if ( KEY > A(M) ) then
                L=M+1
                M=(L+R)/2
            end if
        end do

        BINARY_SEARCH = 0
        return
    end function



subroutine FUTIL()
    write(6,*) "HELLO WORLD"  ! If no "6", there will be an unknow error.
end subroutine FUTIL


subroutine get_numdata(numdata,limit1,limit2)

    character(len=4)    :: obsname
    real(NER)           :: elab,acm,val,error
    character(len=6)    :: ref
    integer             :: numdata
    real(NER)           :: limit1,limit2
    integer             :: ii , istat

    numdata=0

    open(unit=330,file=data_filename,iostat=istat)

    do
        read(330,*,iostat=istat) obsname,elab

        if (istat /= 0) exit
        if (elab >= limit1 .and. elab <= limit2) numdata = numdata + 1
    end do

    close(330)
end subroutine get_numdata

subroutine get_data(expt_data,numdata,limit1,limit2)

    character(len=4)    :: obsname
    real(NER)           :: elab,acm,val,error
    character(len=6)    :: ref
    type(struct_Expt_Data)  :: expt_data(:)
    integer             :: numdata
    real(NER)           :: limit1,limit2
    integer             :: ii , istat

    open(unit=330,file=data_filename)

    ii=0
    do
        read(330,*,iostat=istat) obsname,elab,acm,val,error,ref
        if (istat /= 0 ) exit
        if(elab >= limit1 .and. elab <= limit2) then
            ii = ii + 1
            expt_data(ii)%obsname = obsname
            expt_data(ii)%elab = elab
            expt_data(ii)%acm = acm
            expt_data(ii)%val = val
            expt_data(ii)%error = error
            expt_data(ii)%ref = ref

        end if
    end do
    close(330)


end subroutine get_data

subroutine get_numk(numk,expt_data,numdata)

    type(data_chain),pointer    :: head,tail,ptr,ptr1,ptr2
    real(NER)   :: elab
    character(len=4)    :: obsname
    integer             :: numk
    type(struct_Expt_Data)  :: expt_data(:)
    integer                 :: numdata
    integer                 :: ii,istat

    nullify(head,tail,ptr,ptr1,ptr2)
    numk = 0

    do ii=1,numdata

        elab = expt_data(ii)%elab

        if (.not. associated(head)) then
            allocate(ptr,stat=istat)
            numk=numk+1
            ptr%elab = elab
            head => ptr
            tail => head
            nullify(ptr%next)

        else if (elab < head%elab-0.0001) then
            allocate(ptr,stat=istat)
            numk=numk+1
            ptr%elab = elab
            ptr%next => head
            head => ptr

        else if (elab > tail%elab+0.0001) then
            allocate(ptr,stat=istat)
            numk=numk+1
            ptr%elab = elab
            tail%next => ptr
            tail => ptr
            nullify(ptr%next)

        else
            ptr1 => head
            if ( .not. (elab-ptr1%elab)<=0.0001) then
                ptr2 => ptr1%next
                do
                    if (abs(elab-ptr2%elab)<=0.0001) exit
                    if ( (elab > ptr1%elab +0.0001) .and. (elab < ptr2%elab - 0.0001)) then
                        allocate(ptr,stat=istat)
                        numk=numk+1
                        ptr%elab = elab
                        ptr%next  => ptr2
                        ptr1%next => ptr
                        exit
                    end if
                    ptr1 => ptr2
                    ptr2 => ptr1%next
                end do
            end if
        end if
    end do
end subroutine get_numk


subroutine get_elist(elist,numk,expt_data,numdata)
    integer     :: numk
    real(NER)   :: elist(:)
    real(NER)   :: elab
    type(data_chain),pointer :: head,tail,ptr,ptr1,ptr2
    integer     :: ii , istat
    type(struct_Expt_Data)  :: expt_data(:)
    integer     :: numdata

    nullify(head,tail,ptr,ptr1,ptr2)

    do ii=1,numdata

        elab = expt_data(ii)%elab

            if (.not. associated(head)) then

                allocate(ptr,stat=istat)

                ptr%elab = elab
                head => ptr
                tail => head
                nullify(ptr%next)

            else if (elab < head%elab-0.0001) then
                allocate(ptr,stat=istat)
                ptr%elab = elab
                ptr%next => head
                head => ptr

            else if (elab > tail%elab+0.0001) then
                allocate(ptr,stat=istat)
                ptr%elab = elab
                tail%next => ptr
                tail => ptr
                nullify(ptr%next)

            else
                ptr1 => head
                if ( .not. (elab-ptr1%elab)<=0.0001) then
                    ptr2 => ptr1%next
                    do
                        if (abs(elab-ptr2%elab)<=0.0001) exit
                        if ( (elab > ptr1%elab +0.0001) .and. (elab < ptr2%elab - 0.0001)) then
                            allocate(ptr,stat=istat)
                            ptr%elab = elab
                            ptr%next  => ptr2
                            ptr1%next => ptr
                            exit
                        end if
                        ptr1 => ptr2
                        ptr2 => ptr1%next
                    end do
                end if
            end if
    end do

    ptr => head
    do ii=1,numk
        elist(ii) = head%elab
        head => head%next
    end do

end subroutine get_elist


end subroutine fitdata




end module

