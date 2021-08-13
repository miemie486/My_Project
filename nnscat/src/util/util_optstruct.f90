! Bingwei Long 06/11/2012

! Obsolete on 08/10/2018
module util_optstruct

    use nneft_type
    use eft_potspwd

    implicit none

    type :: flags_struct

        logical :: is_dibaryon = .false., &
                &  two_mn_o_pi = .false., &
                &  is_pert     = .false.
        integer :: order = 0, chnid = CHNID_UNDEFINED, regtype = REGTYPE_GAUSSIAN, &
                &  pots_type = POTS_WPC, mesh_N

! is there a factor of 2m_N/pi in LSE?

    end type flags_struct


    contains


    subroutine set_default_flags(flags)

        type(flags_struct), intent(out)   :: flags

        flags%is_dibaryon = .false.
        flags%is_pert     = .false.
        flags%two_mn_o_pi = .false.
        flags%order       = 0
        flags%chnid       = CHNID_UNDEFINED
        flags%regtype     = REGTYPE_GAUSSIAN
        flags%pots_type   = POTS_WPC
        flags%mesh_N      = 100

    end subroutine set_default_flags


    subroutine word_to_flag(word, flags, hitted)

        character(len = *), intent(in)      :: word
        type(flags_struct), intent(inout)   :: flags
        logical, intent(out)                :: hitted

        type(lsj_symbol)    :: lsj
        integer             :: tmpchnid, ios, tmpN
        logical             :: succ_flag

        tmpchnid = CHNID_UNDEFINED
        call convert_text_to_lsj(word, lsj, succ_flag)
        if (succ_flag) then
            call convert_lsj_to_chnid(lsj, tmpchnid)
            flags%chnid = tmpchnid
            hitted = .true.
            return
        end if
        if ((len(word) .eq. 1) .and. (ichar(word(1:1)) .ge. ichar('0')) .and. (ichar(word(1:1)) .le. ichar('9'))) then
            read(word, '(i1)') flags%order
            hitted = .true.
            return
        end if
        hitted = .true.
        select case (word)
            case ('--gaussian', '-g')
                flags%regtype = REGTYPE_GAUSSIAN
                return
            case ('--sharp', '-s')
                flags%regtype = REGTYPE_SHARP
                return
            case ('--tmop', '-t')
                flags%two_mn_o_pi = .true.
                return
            case ('--dibaryon', '-d')
                flags%is_dibaryon = .true.
                return
            case ('--nondibaryon', '-n')
                flags%is_dibaryon = .false.
                return
!~             case ('--dwave_pert', '-dt')
!~                 flags%pots_type = POTS_DWAVE_PERT
!~                 return
            case ('--pert', '-p')
                flags%is_pert = .true.
                return
!~             case ('--pots_wpc', '-w')
!~                 flags%pots_type = POTS_WPC
!~                 return
            case default
                hitted = .false.
        end select
        if (word(1:3) .eq. '-N=') then
            read(word(4:), '(i5)', iostat = ios) tmpN
            if (ios .eq. 0) then          ! succeed reading N
                hitted = .true.
                flags%mesh_N = tmpN
                return
            end if
        end if
        hitted = .false.

    end subroutine word_to_flag


! un-processed options will be left in out_str

    subroutine sentence_to_flags(src_str, out_str, flags, hitted_once)

        character(len = *), intent(in)      :: src_str
        character(len = *), intent(out)     :: out_str
        type(flags_struct), intent(inout)   :: flags
        logical, intent(out)                :: hitted_once

        integer :: i1, pend, cnt
        character(len = 100) :: word, tmp_str
        logical :: word_hitted

        tmp_str = trim(adjustl(src_str))
        out_str = ''
        cnt = 0
        hitted_once = .false.
        do
            cnt = cnt + 1
            pend = len(tmp_str)
            i1 = scan(tmp_str, ' ')
            if (i1 .eq. 0) then
                word = trim(tmp_str)
            else
                word = tmp_str(1:i1-1)
            end if
            call word_to_flag(word, flags, word_hitted)
            if (.not. word_hitted) then
                out_str = trim(out_str)//''//trim(word)
            else
                if (.not. hitted_once) hitted_once = .true.
            end if
            if (i1 .eq. 0) exit
            if (cnt .gt. 1000) then
                print '(a)', 'sentence_to_flags: too many cycles'
                exit
            end if
            tmp_str = trim(adjustl(tmp_str(i1:pend)))
        end do

    end subroutine sentence_to_flags


end module util_optstruct
