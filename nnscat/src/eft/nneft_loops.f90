! Bingwei Long 10/06/2018
!   - Removed Lambda
! major rev 04/10/2012
! Bingwei Long 12/05/2010

module nneft_loops

  use nneft_type
  use nneft_lsesv, only : idx1d_by_ch_p
  use eft_potspwd, only : regltr_half_PP
  implicit none

  integer, parameter :: VM_IS_SYMMETRIC_LOOP = 10,                &
    &                   VM_IS_NOT_SYMMETRIC_LOOP = 20
  complex, parameter, private :: BETAK = - 0.5_NER*IMGUNT_NE*PI_NE

  contains

  ! Vp(NCH*Np, NCH*Nmesh) G0 T(NCH*Nmesh, NCH*Nk) = loop(NCH*Np, NCH*Nk)
  subroutine get_loop_VpG0Tk_cmplx_nch_nrow_ncol(mE, NCH, Np, Nk, Nmesh, Cwghts, Clsqr, Vp, Tk, loop)

    integer, intent(in)       :: NCH, Np, Nk, Nmesh
    complex(NEC), intent(in)  :: mE
    complex(NEC), intent(in)  :: Cwghts(:), Clsqr(:), Vp(:, :)
    complex(NEC), intent(in)  :: Tk(:, :)
    complex(NEC), intent(out) :: loop(:, :)

    integer     :: ii, jj, kk, aa, bb, cc, pos1, pos2, pos_cq
    complex(NEC)  :: tmp


    loop(1:NCH*Np, 1:NCH*Nk) = 0.0_NER
    do ii = 1, Np
      do aa = 1, NCH
        pos1 = idx1d_by_ch_p(aa, ii, Np)
        do kk = 1, Nk
          do bb = 1, NCH
          pos2 = idx1d_by_ch_p(bb, kk, Nk)
            tmp = 0.0_NER
            do jj = 1, Nmesh
              do cc = 1, NCH
                pos_cq = idx1d_by_ch_p(cc, jj, Nmesh)
                tmp = tmp + Vp(pos1, pos_cq) * Clsqr(jj)           &
                  &  /(mE - Clsqr(jj)) * Cwghts(jj) * Tk(pos_cq, pos2)
              end do
            end do
            loop(pos1, pos2) = tmp
          end do
        end do
      end do
    end do

  end subroutine get_loop_VpG0Tk_cmplx_nch_nrow_ncol

  ! If the upper limit is \infty, the integrand is assumed to be regularized in
  ! the UV region.

  ! \int_0^\infty dl l^2 Vk1(l)*Vk2(l)/(k^2 - l^2 + i\epsilon)
  function loop_1loop_VkVk_sngl(k, N, meshpts, wghts, lsqr, Vk1, Vk2)

    complex(NEC)              :: loop_1loop_VkVk_sngl
    integer, intent(in)       :: N
    real(NER), intent(in)     :: k
    real(NER), intent(in)     :: meshpts(:), wghts(:), lsqr(:), Vk1(:), Vk2(:)

    integer     :: ii
    real(NER)   :: ksqr
    complex(NEC):: k2VV
    real(NER)   :: intgrnd(1:N)

    ksqr   = k*k
    k2VV   = ksqr*Vk1(N+1)*Vk2(N+1)
    forall (ii = 1:N)
      intgrnd(ii) = wghts(ii)*(lsqr(ii)*Vk1(ii)*Vk2(ii) - k2VV)/(ksqr - lsqr(ii))
    end forall

    loop_1loop_VkVk_sngl = sum(intgrnd) + k*BETAK*Vk1(N+1)*Vk2(N+1)

  end function loop_1loop_VkVk_sngl

  ! \int_0^\infty dl l^2 Vk1(l)*Vk2(l)/(mE - l^2)
  function loop_1loop_VkVk_cmplxE_sngl(mE, N, meshpts, wghts, lsqr, Vk1, Vk2)

    complex(NEC)                :: loop_1loop_VkVk_cmplxE_sngl
    complex(NEC), intent(in)    :: mE
    integer, intent(in)         :: N
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), Vk1(:), Vk2(:)

    integer                         :: ii
    complex(NEC), dimension(1:N)    :: intgrnd

    forall (ii = 1:N)
      intgrnd(ii) = wghts(ii)*lsqr(ii)*Vk1(ii)*Vk2(ii)/(mE - lsqr(ii))
    end forall

    loop_1loop_VkVk_cmplxE_sngl = sum(intgrnd)

  end function loop_1loop_VkVk_cmplxE_sngl

  ! \int_0^\infty dl l^2 Vk(l)*Tk(l)/(k^2 - l^2 + i\epsilon)
  function loop_1loop_VkTk_sngl(k, N, meshpts, wghts, lsqr, Vk, Tk)

    complex(NEC)                                :: loop_1loop_VkTk_sngl
    integer, intent(in)                         :: N
    real(NER), intent(in)                       :: k
    real(NER), intent(in)     :: meshpts(:), wghts(:), lsqr(:), Vk(:)
    complex(NEC), intent(in)  :: Tk(:)

    integer     :: ii
    real(NER)   :: ksqr
    complex(NEC):: k2VT
    complex(NEC), dimension(1:N) :: intgrnd

    ksqr   = k*k
    k2VT   = ksqr*Vk(N+1)*Tk(N+1)
    forall (ii = 1:N)
      intgrnd(ii) = wghts(ii)*(lsqr(ii)*Vk(ii)*Tk(ii) - k2VT)/(ksqr - lsqr(ii))
    end forall

    loop_1loop_VkTk_sngl = sum(intgrnd) + k*BETAK*Vk(N+1)*Tk(N+1)

  end function loop_1loop_VkTk_sngl

  ! \int_0^\infty dl l^2 Tk(l)*Tk(l)/(k^2 - l^2 + i\epsilon)
  subroutine get_1loop_Tk1G0Tk2_slp(k, N, meshpts, wghts, lsqr, Tk1,   &
    & Tk2, val)

    integer, intent(in)       :: N
    real(NER), intent(in)     :: k
    real(NER), intent(in)     :: meshpts(:), wghts(:), lsqr(:)
    complex(NEC), intent(in)  :: Tk1(:), Tk2(:)
    complex(NEC), intent(out) :: val

    integer       :: ii
    real(NER)     :: ksqr
    complex(NEC)  :: k2TT

    ksqr   = k*k
    k2TT   = ksqr*Tk1(N+1)*Tk2(N+1)
    val = k*BETAK*Tk1(N+1)*Tk2(N+1)
    do ii = 1, N
      val = val + wghts(ii)*(lsqr(ii)*Tk1(ii)*Tk2(ii) - k2TT)/(ksqr - lsqr(ii))
    end do

  end subroutine get_1loop_Tk1G0Tk2_slp

  ! subroutine get_1loop_Tk1G0Tk2_slp(k, N, meshpts, wghts, lsqr, Tk1,   &
  !   & Tk2, val)

  !   integer, intent(in)       :: N
  !   real(NER), intent(in)     :: k
  !   real(NER), intent(in)     :: meshpts(:), wghts(:), lsqr(:)
  !   complex(NEC), intent(in)  :: Tk1(:), Tk2(:)
  !   complex(NEC), intent(out) :: val

  !   integer       :: ii
  !   real(NER)     :: ksqr
  !   complex(NEC)  :: k2TT
  !   complex(NEC)  :: intgrnd(1:N)

  !   ksqr   = k*k
  !   k2TT   = ksqr*Tk1(N+1)*Tk2(N+1)
  !   forall (ii = 1:N)
  !     intgrnd(ii) = wghts(ii)*(lsqr(ii)*Tk1(ii)*Tk2(ii) - k2TT)/(ksqr - lsqr(ii))
  !   end forall

  !   val = sum(intgrnd) + k*BETAK*Tk1(N+1)*Tk2(N+1)

  ! end subroutine get_1loop_Tk1G0Tk2_slp

  ! \int_0^\infty dl l^2 VM(k, l)*Tk(l)/(k^2 - l^2 + i\epsilon)
  function loop_1loop_os_VMTk_sngl(k, N, meshpts, wghts, lsqr, VM, Tk)

    complex(NEC)              :: loop_1loop_os_VMTk_sngl
    integer, intent(in)       :: N
    real(NER), intent(in)     :: k, meshpts(:), wghts(:), lsqr(:), VM(:, :)
    complex(NEC), intent(in)  :: Tk(:)

    integer       :: ii
    real(NER)     :: ksqr
    complex(NEC)  :: k2VT, intgrnd(1:N)

    ksqr   = k*k
    k2VT   = ksqr*VM(N+1, N+1)*Tk(N+1)
    forall (ii = 1:N)
      intgrnd(ii) = wghts(ii)*(lsqr(ii)*VM(N+1, ii)*Tk(ii) - k2VT)/(ksqr - lsqr(ii))
    end forall

    loop_1loop_os_VMTk_sngl = sum(intgrnd) + k*BETAK*VM(N+1, N+1)*Tk(N+1)

  end function loop_1loop_os_VMTk_sngl

  ! Tk_new(p') = \int_0^\infty dl l^2 VM(p', l)*Tk(l)/(k^2 - l^2 + i\epsilon)
  subroutine loop_1loop_Tk_VMTk_sngl(k, N, meshpts, wghts, lsqr, VM, Tk)

    integer, intent(in)                             :: N
    real(NER), intent(in)                           :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), VM(:, :)
    complex(NEC), intent(inout)   :: Tk(:)

    integer         :: ii, jj
    real(NER)       :: ksqr
    complex(NEC)    :: betaok
    complex(NEC)    :: k2ViT(1+N), grnd(N+1, N)

    ksqr = k*k
    betaok = BETAK/k
    forall (ii = 1:N+1)
      k2ViT(ii) = ksqr*VM(ii, N+1)*Tk(N+1)
    end forall
    forall (ii = 1:N+1, jj = 1:N)
      grnd(ii, jj) = wghts(jj)*(lsqr(jj)*VM(ii, jj)*Tk(jj) - k2ViT(ii))/(ksqr - lsqr(jj))
    end forall
    forall (ii = 1:N+1)
      Tk(ii) = sum(grnd(ii, 1:N)) + betaok*k2ViT(ii)
    end forall

  end subroutine loop_1loop_Tk_VMTk_sngl


  ! Tk(p') = \int dl l^2/(k^2-l^2) VM(p', l) VM(l, k)
  ! Tested on  08/03/2018  12:03am
  subroutine loop_1loop_Tk_VMVM_sngl(k, N, meshpts, wghts, lsqr, VM, Tk)

    integer, intent(in)         :: N
    real(NER), intent(in)       :: k, meshpts(:), wghts(:), lsqr(:), VM(:, :)
    complex(NEC), intent(out)   :: Tk(:)

    integer         :: ii, jj
    real(NER)       :: ksqr
    complex(NEC)    :: betaok
    complex(NEC)    :: k2ViV(1+N), grnd(N+1, N)

    ksqr = k*k
    betaok = BETAK/k
    forall (ii = 1:N+1)
      k2ViV(ii) = ksqr*VM(ii, N+1)*VM(N+1, N+1)
    end forall
    forall (ii = 1:N+1, jj = 1:N)
      grnd(ii, jj) = wghts(jj)*(lsqr(jj)*VM(ii, jj)*VM(jj, N+1) - k2ViV(ii))/(ksqr - lsqr(jj))
    end forall
    forall (ii = 1:N+1)
      Tk(ii) = sum(grnd(ii, 1:N)) + betaok*k2ViV(ii)
    end forall

  end subroutine loop_1loop_Tk_VMVM_sngl

  ! TM(p', p; k) = \int dl l^2 T1(p', l) (k^2-l^2+i0)^(-1) T2(l, p)
  ! Tested on  08/03/2018  12:03am
  subroutine get_1loop_TM1G0TM2_slp(k, N, meshpts, wghts, lsqr, T1, T2, TM)

    integer, intent(in)         :: N
    real(NER), intent(in)       :: k, meshpts(:), wghts(:), lsqr(:)
    complex(NEC), intent(in)    :: T1(:, :), T2(:, :)
    complex(NEC), intent(out)   :: TM(:, :)

    integer         :: ii, jj, ll
    real(NER)       :: ksqr
    complex(NEC)    :: betaok, tmp,k2T1iT2j
    ! complex(NEC), dimension(1:N+1, 1:N+1) :: k2T1iT2j
    ! complex(NEC), dimension(1:N+1, 1:N+1, 1:N) :: grnd

    ksqr = k*k
    betaok = BETAK/k
    do ii = 1, N+1
      do jj = 1, N+1
        k2T1iT2j = ksqr*T1(ii, N+1)*T2(N+1, jj)
        tmp = betaok * k2T1iT2j
        do ll = 1, N
          tmp = tmp + wghts(ll)*(lsqr(ll)*T1(ii, ll)*T2(ll, jj) - k2T1iT2j)/(ksqr - lsqr(ll))
        end do
        TM(ii, jj) = tmp
      end do
    end do

  end subroutine get_1loop_TM1G0TM2_slp

  ! subroutine get_1loop_TM1G0TM2_slp(k, N, meshpts, wghts, lsqr, T1, T2, TM)

  !   integer, intent(in)         :: N
  !   real(NER), intent(in)       :: k, meshpts(:), wghts(:), lsqr(:)
  !   complex(NEC), intent(in)    :: T1(:, :), T2(:, :)
  !   complex(NEC), intent(out)   :: TM(:, :)

  !   integer         :: ii, jj, ll
  !   real(NER)       :: ksqr
  !   complex(NEC)    :: betaok
  !   complex(NEC), dimension(1:N+1, 1:N+1) :: k2T1iT2j
  !   complex(NEC), dimension(1:N+1, 1:N+1, 1:N) :: grnd

  !   ksqr = k*k
  !   betaok = BETAK/k
  !   forall (ii = 1:N+1, jj = 1:N+1)
  !     k2T1iT2j(ii, jj) = ksqr*T1(ii, N+1)*T2(N+1, jj)
  !   end forall
  !   forall (ii = 1:N+1, jj = 1:N+1, ll = 1:N)
  !     grnd(ii, jj, ll) = wghts(ll)*(lsqr(ll)*T1(ii, ll)*T2(ll, jj) - k2T1iT2j(ii, jj))/(ksqr - lsqr(ll))
  !   end forall
  !   forall (ii = 1:N+1, jj = 1:N+1)
  !     TM(ii, jj) = sum(grnd(ii, jj, 1:N)) + betaok*k2T1iT2j(ii, jj)
  !   end forall

  ! end subroutine get_1loop_TM1G0TM2_slp

  ! \int_0^\infty dl' l'^2 dl l^2 Tk(l')/(k^2 - l^2 + i0)*VM(l', l)*Tk(l)/(k^2 - l^2 + i0)
  ! VM is assumed to be symmetric
  function loop_2loops_TkVMTk_sym_sngl(k, N, meshpts, wghts, lsqr, VM, Tk)

    complex(NEC)              :: loop_2loops_TkVMTk_sym_sngl
    integer, intent(in)       :: N
    real(NER), intent(in)     :: k, meshpts(:), wghts(:), lsqr(:), VM(:, :)
    complex(NEC), intent(in)  :: Tk(:)

    real(NER)     :: ksqr, wokml(1:N)
    integer       :: ii, jj
    complex(NEC)  :: beta, k2T, k2TV, k4TTV, l2TiVi(1:N), l2Ti(1:N),           &
      &              intgrnd_NX1(1:N), intgrnd_NXN(1:N, 1:N)

    beta   = BETAK
    ksqr   = k*k
    k2T    = ksqr*Tk(N+1)
    k2TV   = k2T*VM(N+1, N+1)
    k4TTV  = k2T*k2TV
    intgrnd_NXN = 0.0_NER
    forall (ii = 1:N)
      l2Ti(ii)   = lsqr(ii)*Tk(ii)
      l2TiVi(ii) = l2Ti(ii)*VM(ii, N+1)
      wokml(ii)  = wghts(ii)/(ksqr - lsqr(ii))
    end forall
    forall (ii = 1:N)
      forall (jj = 1:ii-1)
      intgrnd_NXN(ii, jj) = wokml(ii)*wokml(jj)*(l2Ti(ii)*VM(ii, jj)*l2Ti(jj) - k2T*(l2TiVi(ii) + l2TiVi(jj)) + k4TTV)
      end forall
    end forall
    forall (ii = 1:N)
      intgrnd_NXN(ii, ii) = wokml(ii)*wokml(ii)*(0.5_NER*(l2Ti(ii)*l2Ti(ii)*VM(ii, ii) + k4TTV) - k2T*l2TiVi(ii))
      intgrnd_NX1(ii) = wokml(ii)*(l2TiVi(ii) - k2TV)
    end forall

    loop_2loops_TkVMTk_sym_sngl = 2.0_NER*sum(intgrnd_NXN) + 2.0_NER*k*beta*Tk(N+1)*sum(intgrnd_NX1) + beta*beta*k2TV*Tk(N+1)

  end function loop_2loops_TkVMTk_sym_sngl

  ! Tpq is assumed to be symmetric
  function loop_2loops_TkTpqTk_sym_sngl(k, N, meshpts, wghts, lsqr, Tpq, Tk)

    complex(NEC)              :: loop_2loops_TkTpqTk_sym_sngl
    integer, intent(in)       :: N
    real(NER), intent(in)     :: k, meshpts(:), wghts(:), lsqr(:), Tpq(:, :)
    complex(NEC), intent(in)  :: Tk(:)

    real(NER)                       :: ksqr
    integer                         :: ii, jj
    complex(NEC)                    :: beta, k2T, k2TV, k4TTV
    real(NER), dimension(1:N)       :: wokml
    complex(NEC), dimension(1:N)    :: l2TiVi, l2Ti, intgrnd_NX1
    complex(NEC), dimension(1:N, 1:N) :: intgrnd_NXN

    beta   = BETAK
    ksqr   = k*k
    k2T    = ksqr*Tk(N+1)
    k2TV   = k2T*Tpq(N+1, N+1)
    k4TTV  = k2T*k2TV
    intgrnd_NXN = 0.0_NER
    forall (ii = 1:N)
      l2Ti(ii)   = lsqr(ii)*Tk(ii)
      l2TiVi(ii) = l2Ti(ii)*Tpq(ii, N+1)
      wokml(ii)  = wghts(ii)/(ksqr - lsqr(ii))
    end forall
    forall (ii = 1:N)
      forall (jj = 1:ii-1)
      intgrnd_NXN(ii, jj) = wokml(ii)*wokml(jj)*(l2Ti(ii)*Tpq(ii, jj)*l2Ti(jj) - k2T*(l2TiVi(ii) + l2TiVi(jj)) + k4TTV)
      end forall
    end forall
    forall (ii = 1:N)
      intgrnd_NXN(ii, ii) = wokml(ii)*wokml(ii)*(0.5_NER*(l2Ti(ii)*l2Ti(ii)*Tpq(ii, ii) + k4TTV) - k2T*l2TiVi(ii))
      intgrnd_NX1(ii) = wokml(ii)*(l2TiVi(ii) - k2TV)
    end forall

    loop_2loops_TkTpqTk_sym_sngl = 2.0_NER*sum(intgrnd_NXN) + 2.0_NER*k*beta*Tk(N+1)*sum(intgrnd_NX1) + beta*beta*k2TV*Tk(N+1)

  end function loop_2loops_TkTpqTk_sym_sngl


  ! Tpq is assumed to be symmetric
  function loop_2loops_VkTpqVk_sym_sngl(k, N, meshpts, wghts, lsqr, Vk, Tpq)

    complex(NEC)                                    :: loop_2loops_VkTpqVk_sym_sngl
    integer, intent(in)                             :: N
    real(NER), intent(in)                           :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), Vk(:)
    complex(NEC), intent(in)  :: Tpq(:, :)

    real(NER)                       :: ksqr
    integer                         :: ii, jj
    complex(NEC)                    :: beta, k2V, k2VT, k4VVT
    real(NER), dimension(1:N)       :: wokml
    complex(NEC), dimension(1:N)    :: l2ViTi, l2Vi, intgrnd_NX1
    complex(NEC), dimension(1:N, 1:N) :: intgrnd_NXN

    beta   = BETAK
    ksqr   = k*k
    k2V    = ksqr*Vk(N+1)
    k2VT   = k2V*Tpq(N+1, N+1)
    k4VVT  = k2V*k2VT
    intgrnd_NXN = 0.0_NER
    forall (ii = 1:N)
      l2Vi(ii)   = lsqr(ii)*Vk(ii)
      l2ViTi(ii) = l2Vi(ii)*Tpq(ii, N+1)
      wokml(ii)  = wghts(ii)/(ksqr - lsqr(ii))
    end forall
    forall (ii = 1:N)
      forall (jj = 1:ii-1)
        intgrnd_NXN(ii, jj) = wokml(ii)*wokml(jj)*(l2Vi(ii)*Tpq(ii, jj)*l2Vi(jj) - k2V*(l2ViTi(ii) + l2ViTi(jj)) + k4VVT)
      end forall
    end forall
    forall (ii = 1:N)
      intgrnd_NXN(ii, ii) = wokml(ii)*wokml(ii)*(0.5_NER*(l2Vi(ii)*l2Vi(ii)*Tpq(ii, ii) + k4VVT) - k2V*l2ViTi(ii))
      intgrnd_NX1(ii) = wokml(ii)*(l2ViTi(ii) - k2VT)
    end forall

    loop_2loops_VkTpqVk_sym_sngl = 2.0_NER*sum(intgrnd_NXN) + 2.0_NER*k*beta*Vk(N+1)*sum(intgrnd_NX1) + beta*beta*k2VT*Vk(N+1)

  end function loop_2loops_VkTpqVk_sym_sngl


  ! Tpq is assumed to be symmetric
  ! \int_0^\infty dl' l'^2 dl l^2 Vk(l')/(mE - l^2)*Tpq(l', l)*Vk(l)/(mE - l^2)
  ! mE is on the physical sheet and off the positive half of the real axis
  function loop_2loops_VkTpqVk_sym_cmplxE_sngl(mE, N, meshpts, wghts, lsqr, Vk, Tpq)

    complex(NEC)                                        :: loop_2loops_VkTpqVk_sym_cmplxE_sngl
    complex(NEC), intent(in)                            :: mE
    integer, intent(in)                                 :: N
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), Vk(:)
    complex(NEC), intent(in)  :: Tpq(:, :)

    integer                         :: ii, jj
    complex(NEC), dimension(1:N)    :: wokml
    complex(NEC), dimension(1:N)    :: l2Vi
    complex(NEC), dimension(1:N, 1:N) :: intgrnd_NXN

    intgrnd_NXN = 0.0_NER
    forall (ii = 1:N)
      l2Vi(ii)   = lsqr(ii)*Vk(ii)
      wokml(ii)  = wghts(ii)/(mE - lsqr(ii))
    end forall
    forall (ii = 1:N)
      forall (jj = 1:ii-1)
        intgrnd_NXN(ii, jj) = wokml(ii)*wokml(jj)*l2Vi(ii)*Tpq(ii, jj)*l2Vi(jj)
      end forall
    end forall
    forall (ii = 1:N)
      intgrnd_NXN(ii, ii) = wokml(ii)*wokml(ii)*0.5_NER*l2Vi(ii)*l2Vi(ii)*Tpq(ii, ii)
    end forall

    loop_2loops_VkTpqVk_sym_cmplxE_sngl = 2.0_NER*sum(intgrnd_NXN)

  end function loop_2loops_VkTpqVk_sym_cmplxE_sngl

  function loop_2loops_Vk1TpqVk2_asym_sngl(k, N, meshpts, wghts, lsqr, Vk1, Vk2, Tpq)

    complex(NEC)                                    :: loop_2loops_Vk1TpqVk2_asym_sngl
    integer, intent(in)                             :: N
    real(NER), intent(in)                           :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), Vk1(:), Vk2(:)
    complex(NEC), intent(in)  :: Tpq(:, :)

    real(NER)                       :: ksqr
    integer                         :: ii, jj
    complex(NEC)                    :: beta, k2V1, k2V2, k2V1T, k2TV2, k4V1TV2
    real(NER), dimension(1:N)       :: wokml
    complex(NEC), dimension(1:N)    :: l2V1i, l2V1iTi, l2V2i, l2TiV2i, intgrnd_NX1
    complex(NEC), dimension(1:N, 1:N) :: intgrnd_NXN

    beta   = BETAK
    ksqr   = k*k
    k2V1   = ksqr*Vk1(N+1)
    k2V2   = ksqr*Vk2(N+1)
    k2V1T  = k2V1*Tpq(N+1, N+1)
    k2TV2  = Tpq(N+1, N+1)*k2V2
    k4V1TV2 = k2V1*k2TV2
    forall (ii = 1:N)
      l2V1i(ii)  = lsqr(ii)*Vk1(ii)
      l2V1iTi(ii)= l2V1i(ii)*Tpq(ii, N+1)
      l2V2i(ii)  = lsqr(ii)*Vk2(ii)
      l2TiV2i(ii)= Tpq(N+1, ii)*l2V2i(ii)
      wokml(ii)  = wghts(ii)/(ksqr - lsqr(ii))
    end forall
    forall (ii = 1:N, jj = 1:N)
      intgrnd_NXN(ii, jj) = wokml(ii)*wokml(jj)*(l2V1i(ii)*Tpq(ii, jj)*l2V2i(jj) &
        & - l2V1iTi(ii)*k2V2 - k2V1*l2TiV2i(jj) + k4V1TV2)
    end forall
    forall (ii = 1:N)
      intgrnd_NX1(ii) = wokml(ii)*(k2V1*l2TiV2i(ii) + l2V1iTi(ii)*k2V2 - 2.0_NER*k4V1TV2)
    end forall

    loop_2loops_Vk1TpqVk2_asym_sngl = sum(intgrnd_NXN) + beta/k*sum(intgrnd_NX1) + beta*beta*k4V1TV2/ksqr

  end function loop_2loops_Vk1TpqVk2_asym_sngl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                  !
!                      coupled-channel loops                       !
!                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! \int_0^\infty dl l^2 VM(k, l)*Tk(l)/(k^2 - l^2 + i\epsilon)
  ! VM(k, l) and Tk(l) are 2X2 matrices
  subroutine loop_1loop_os_VMTk_cpld(k, N, meshpts, wghts, lsqr, VM, Tk, L1)

    integer, intent(in)                                     :: N
    real(NER), intent(in)                                   :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), VM(:, :, :, :)
    complex(NEC), intent(in)  :: Tk(:, :, :)
    complex(NEC), intent(out)  :: L1(1:2, 1:2)

    integer                             :: ii, jj
    real(NER)                           :: ksqr
    complex(NEC), dimension(1:2, 1:2)   :: k2VT
    complex(NEC), dimension(1:2, 1:2, N):: grnd

    ksqr   = k*k
    k2VT   = ksqr*matmul(VM(1:2, 1:2, N+1, N+1), Tk(1:2, 1:2, N+1))
    forall (ii = 1:N)
      grnd(1:2, 1:2, ii) = wghts(ii)/(ksqr - lsqr(ii))*(lsqr(ii)*matmul(VM(1:2, 1:2, N+1, ii), Tk(1:2, 1:2, ii)) - k2VT)
    end forall
    forall (ii = 1:2, jj = 1:2)
      L1(ii, jj) = sum(grnd(ii, jj, 1:N))
    end forall

    L1 = L1 + BETAK/k*k2VT

  end subroutine loop_1loop_os_VMTk_cpld

  ! Frozen
  ! ! \int_0^\infty dl l^2 VMa(k, l)*VMb(l, k)/(k^2 - l^2 + i\epsilon)
  ! ! VMa(k, l) and VMb(l, k) are 2X2 matrices
  ! subroutine Get1LoopVMVMOSCpld_loop(k, N, meshpts, wghts, lsqr, VMa, VMb, L1)

  !   integer, intent(in)    :: N
  !   real(NER), intent(in)  :: k, meshpts(:), wghts(:), lsqr(:),        &
  !   & VMa(:, :, :, :), VMb(:, :, :, :)
  !   complex(NEC), intent(out)  :: L1(1:2, 1:2)

  !   integer       :: ii, jj
  !   real(NER)     :: ksqr
  !   complex(NEC)  :: k2VT(1:2, 1:2), grnd(1:2, 1:2, N)

  !   ksqr   = k*k
  !   k2VT   = ksqr*matmul(VMa(1:2, 1:2, N+1, N+1), VMb(1:2, 1:2, N+1, N+1))
  !   forall (ii = 1:N)
  !     grnd(1:2, 1:2, ii) = wghts(ii)/(ksqr - lsqr(ii))*(lsqr(ii)               &
  !     &   *matmul(VMa(1:2, 1:2, N+1, ii), VMb(1:2, 1:2, ii, N+1)) - k2VT)
  !   end forall
  !   forall (ii = 1:2, jj = 1:2)
  !     L1(ii, jj) = sum(grnd(ii, jj, 1:N))
  !   end forall

  !   L1 = L1 + BETAK/k*k2VT

  ! end subroutine Get1LoopVMVMOSCpld_loop


  ! \int_0^\infty dl' l'^2 dl l^2 Tk(l')/(k^2 - l^2 + i0)*VM(l', l)*Tk(l)/(k^2 - l^2 + i0)
  ! VM is assumed to be symmetric
  subroutine loop_2loops_TkVMTk_sym_cpld(k, N, meshpts, wghts, lsqr, VM, Tk, L2)

    integer, intent(in)                                     :: N
    real(NER), intent(in)                                   :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), VM(:, :, :, :)
    complex(NEC), intent(in)  :: Tk(:, :, :)
    complex(NEC), intent(out)  :: L2(1:2, 1:2)

    real(NER)                           :: ksqr
    integer                             :: ii, jj, aa, bb
    complex(NEC)                        :: betaok
    complex(NEC), dimension(1:2, 1:2)   :: k2T, k4TVT
    real(NER), dimension(1:N)           :: wokml
    complex(NEC), dimension(1:2, 1:2, 1:N)  :: k2l2TVNiTi, l2Ti, grnd_1d
    complex(NEC), dimension(1:2, 1:2, 1:N, 1:N) :: grnd_2d

    betaok = BETAK/k
    ksqr   = k*k
    k2T    = ksqr*Tk(1:2, 1:2, N+1)
    k4TVT  = matmul(matmul(k2T, VM(1:2, 1:2, N+1, N+1)), k2T)
    forall (ii = 1:N)
      l2Ti(1:2, 1:2, ii)       = lsqr(ii)*Tk(1:2, 1:2, ii)
      k2l2TVNiTi(1:2, 1:2, ii) = matmul(k2T, matmul(VM(1:2, 1:2, N+1, ii), l2Ti(1:2, 1:2, ii)))
      wokml(ii) = wghts(ii)/(ksqr - lsqr(ii))
    end forall
    forall (ii = 1:N)
      forall (jj = 1:ii-1)
        grnd_2d(1:2, 1:2, ii, jj) = wokml(ii)*wokml(jj)*(matmul(transpose(l2Ti(1:2, 1:2, ii)), matmul(VM(1:2, 1:2, ii, jj), l2Ti(1:2, 1:2, jj))) - (k2l2TVNiTi(1:2, 1:2, jj) + transpose(k2l2TVNiTi(1:2, 1:2, ii))) + k4TVT)
        grnd_2d(1:2, 1:2, jj, ii) = transpose(grnd_2d(1:2, 1:2, ii, jj))
      end forall
    end forall
    forall (ii = 1:N)
      grnd_2d(1:2, 1:2, ii, ii) = wokml(ii)*wokml(ii)*(matmul(transpose(l2Ti(1:2, 1:2, ii)), matmul(VM(1:2, 1:2, ii, ii), l2Ti(1:2, 1:2, ii))) - (k2l2TVNiTi(1:2, 1:2, ii) + transpose(k2l2TVNiTi(1:2, 1:2, ii))) + k4TVT)
      grnd_1d(1:2, 1:2, ii) = wokml(ii)*(transpose(k2l2TVNiTi(1:2, 1:2, ii)) + k2l2TVNiTi(1:2, 1:2, ii) - 2.0_NER*k4TVT)
    end forall
    forall (ii = 1:2, jj = 1:2)
      L2(ii, jj) = sum(grnd_2d(ii, jj, 1:N, 1:N)) + betaok*sum(grnd_1d(ii, jj, 1:N)) + betaok*betaok*k4TVT(ii, jj)
    end forall

  end subroutine loop_2loops_TkVMTk_sym_cpld


  ! \int dl l^2/(k^2 - l^2 + i0) Vk(l)*T_{a b}(l, k; k)
  function loop_1loop_VkTk_ab_cpld(k, N, meshpts, wghts, lsqr, Vk, Tk, a, b)

    complex(NEC)                                :: loop_1loop_VkTk_ab_cpld
    integer, intent(in)                         :: N, a, b
    real(NER), intent(in)                       :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), Vk(:)
    complex(NEC), intent(in)  :: Tk(:, :, :)

    integer     :: ii
    real(NER)   :: ksqr
    complex(NEC):: regsum, k2VT

    if (a .lt. 1 .or. a .gt. 2 .or. b .lt. 1 .or. b .gt. 2) then
      loop_1loop_VkTk_ab_cpld = 0.0_NER
      print ('(a)'), 'loop_1loop_VkTk_ab_cpld: Tk is not 2X2'
      return
    end if
    regsum = 0.0_NER
    ksqr   = k*k
    k2VT   = ksqr*Vk(N+1)*Tk(a, b, N+1)
    do ii = 1, N
      regsum = regsum + wghts(ii)/(ksqr - lsqr(ii))*(lsqr(ii)*Vk(ii)*Tk(a, b, ii) - k2VT)
    end do

    loop_1loop_VkTk_ab_cpld = regsum + BETAK/k*k2VT

  end function loop_1loop_VkTk_ab_cpld


  ! T(p', k; k) = \int_0^\infty dl l^2 VM(p', l)*VM(l, k)/(k^2 - l^2 + i\epsilon)
  subroutine loop_1loop_Tk_VMVM_cpld(k, N, meshpts, wghts, lsqr, VM, Tk)

    integer,                                 intent(in)     :: N
    real(NER),                               intent(in)     :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:), VM(:, :, :, :)
    complex(NEC), intent(out)  :: Tk(:, :, :)

    integer         :: ii, jj, aa, bb
    real(NER)       :: ksqr
    complex(NEC)    :: betaok
    complex(NEC), dimension(1:2, 1:2, 1:N+1)        :: k2ViV
    complex(NEC), dimension(1:2, 1:2, 1:N+1, 1:N)   :: grnd

    ksqr   = k*k
    betaok = BETAK/k
    forall (ii = 1:N+1)
      k2ViV(1:2, 1:2, ii)  = ksqr*matmul(VM(1:2, 1:2, ii, N+1), VM(1:2, 1:2, N+1, N+1))
    end forall
    forall (ii = 1:N+1, jj = 1:N)
      grnd(1:2, 1:2, ii, jj) = wghts(jj)/(ksqr - lsqr(jj))*(lsqr(jj)*matmul(VM(1:2, 1:2, ii, jj), VM(1:2, 1:2, jj, N+1)) - k2ViV(1:2, 1:2, ii))
    end forall
    forall (ii = 1:N+1, aa = 1:2, bb = 1:2)
      Tk(aa, bb, ii) = sum(grnd(aa, bb, ii, 1:N)) + betaok*k2ViV(aa, bb, ii)
    end forall

  end subroutine loop_1loop_Tk_VMVM_cpld


  !   \int_0^\infty dl l^2 T(k, l; k)*T(l, k; k)/(k^2 - l^2 + i\epsilon)
  ! = \int_0^\infty dl l^2 transpose(T(l, k; k))*T(l, k; k)/(k^2 - l^2 + i\epsilon)
  subroutine loop_1loop_TkTk_cpld(k, N, meshpts, wghts, lsqr, Tk, rslt)

    integer, intent(in)                             :: N
    real(NER), intent(in)                           :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), lsqr(:)
    complex(NEC), intent(in)    :: Tk(:, :, :)
    complex(NEC), intent(out)   :: rslt(1:2, 1:2)

    integer                         :: ii, jj
    real(NER)                       :: ksqr
    complex(NEC), dimension(1:2, 1:2)   :: k2TT
    complex(NEC), dimension(1:2, 1:2, N):: grnd

    ksqr   = k*k
    ! NOTE: transpose(Tk(N+1)) == Tk(N+1)
    k2TT   = ksqr*matmul(Tk(1:2, 1:2, N+1), Tk(1:2, 1:2, N+1))
    forall (ii = 1:N)
      grnd(1:2, 1:2, ii) = wghts(ii)/(ksqr - lsqr(ii))*(lsqr(ii)*matmul(transpose(Tk(1:2, 1:2, ii)), Tk(1:2, 1:2, ii)) - k2TT)
    end forall
    forall (ii = 1:2, jj = 1:2)
      rslt(ii, jj) = sum(grnd(ii, jj, 1:N))
    end forall

    rslt = rslt + BETAK/k*k2TT

  end subroutine loop_1loop_TkTk_cpld


end module nneft_loops


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!                       recycle
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




