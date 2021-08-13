! Bingwei Long 12/14/2012

module nneft_lsesv

  use nneft_type
  implicit none

  integer, parameter :: lsesv_SUCCEEDED = 0, lsesv_FAILED = -1

  private :: lsesv_EPS
  real(NER) :: lsesv_EPS = 1.0E-15_NER

  contains

  ! aa : index of channel
  ! ii : index of momentum mesh point
  ! NP  : number of p-mesh points
  ! map (aa, ii) to a 1d index
  pure function idx1d_by_ch_p(aa, ii, NP)

    integer, intent(in) :: aa, ii, NP
    integer             :: idx1d_by_ch_p

    idx1d_by_ch_p = (aa-1)*NP + ii

  end function idx1d_by_ch_p

  ! Extract potval(1:NCH, 1:NCH) from LargeVM that stores pot values for each
  ! p and ch
  ! LargeVM((a-1)*NP + i, (b-1)*NP + j) = SmallVM(msh(i), msh(j); E)
  ! If potM has the (NP+1)-th row and (NP+1)-th column, use
  ! call ExtractSmallVM(LargeVM, NCH, NP+1, ii, jj, SmallVM)
  subroutine ExtractSmallVM(LargeVM, NCH, NP, ii, jj, SmallVM)

    real(NER), intent(in)     :: LargeVM(:, :)
    integer, intent(in)       :: NCH, NP, ii, jj
    real(NER), intent(inout)  :: SmallVM(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        ! SmallVM(aa, bb) = LargeVM((aa-1)*NP + ii, (bb-1)*NP + jj)
        SmallVM(aa, bb) = LargeVM(idx1d_by_ch_p(aa, ii, NP), idx1d_by_ch_p(bb, jj, NP))
      end do
    end do

  end subroutine ExtractSmallVM

  ! If potM has the (NP+1) row and column, use
  ! call SetLargeVM(LargeVM, NCH, NP+1, ii, jj, SmallVM)
  subroutine SetLargeVM(LargeVM, NCH, NP, ii, jj, SmallVM)

    real(NER), intent(inout) :: LargeVM(:, :)
    integer, intent(in)      :: NCH, NP, ii, jj
    real(NER), intent(inout) :: SmallVM(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        LargeVM(idx1d_by_ch_p(aa, ii, NP), idx1d_by_ch_p(bb, jj, NP)) = SmallVM(aa, bb)
      end do
    end do

  end subroutine SetLargeVM

  subroutine SetLargeVM_cmplx(LargeVM, NCH, NP, ii, jj, SmallVM)

    complex(NEC), intent(inout) :: LargeVM(:, :)
    integer, intent(in)         :: NCH, NP, ii, jj
    complex(NEC), intent(inout) :: SmallVM(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        LargeVM(idx1d_by_ch_p(aa, ii, NP), idx1d_by_ch_p(bb, jj, NP)) = SmallVM(aa, bb)
      end do
    end do

  end subroutine SetLargeVM_cmplx

  ! Array rank:
  !   LargeVM(NCH*Nrow, NCH*Ncol)
  !   SmallVM(1:NCH, 1:NCH)
  subroutine SetLargeVM_generic_cmplx(LargeVM, NCH, Nrow, Ncol, ii, jj, SmallVM)

    complex(NEC), intent(inout) :: LargeVM(:, :)
    integer, intent(in)         :: NCH, Nrow, Ncol, ii, jj
    complex(NEC), intent(inout) :: SmallVM(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        LargeVM(idx1d_by_ch_p(aa, ii, Nrow), idx1d_by_ch_p(bb, jj, Ncol)) = SmallVM(aa, bb)
      end do
    end do

  end subroutine SetLargeVM_generic_cmplx

  subroutine ExtractSmallTM(LargeTM, NCH, NP, ii, jj, SmallTM)

    complex(NEC), intent(in)     :: LargeTM(:, :)
    integer, intent(in)          :: NCH, NP, ii, jj
    complex(NEC), intent(inout)  :: SmallTM(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        SmallTM(aa, bb) = LargeTM(idx1d_by_ch_p(aa, ii, NP), idx1d_by_ch_p(bb, jj, NP))
      end do
    end do

  end subroutine ExtractSmallTM

  subroutine SetLargeTM(LargeTM, NCH, NP, ii, jj, SmallTM)

    complex(NEC), intent(inout) :: LargeTM(:, :)
    integer, intent(in)         :: NCH, NP, ii, jj
    complex(NEC), intent(inout) :: SmallTM(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        LargeTM(idx1d_by_ch_p(aa, ii, NP), idx1d_by_ch_p(bb, jj, NP)) = SmallTM(aa, bb)
      end do
    end do

  end subroutine SetLargeTM

  ! Extract TQ from Tk
  ! TQ(a, b)
  subroutine GetTQFromTk(Tk, NCH, NP, ii, TQ)

    complex(NEC), intent(in)     :: Tk(:, :)
    integer, intent(in)          :: NCH, NP, ii
    complex(NEC), intent(inout)  :: TQ(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        TQ(aa, bb) = Tk(idx1d_by_ch_p(aa, ii, NP), bb)
      end do
    end do

  end subroutine GetTQFromTk

  subroutine SetTkFromTQ(Tk, NCH, NP, ii, TQ)

    complex(NEC), intent(out) :: Tk(:, :)
    integer, intent(in)       :: NCH, NP, ii
    complex(NEC), intent(in)  :: TQ(:, :)

    integer :: aa, bb

    do aa = 1, NCH
      do bb = 1, NCH
        Tk(idx1d_by_ch_p(aa, ii, NP), bb) = TQ(aa, bb)
      end do
    end do

  end subroutine SetTkFromTQ

  ! Lippmann-Schwinger equation of "T-matrix", T(p', p; k) with +i\epsilon prescription.
  !
  ! T(p', p; k) = v(p', p; E) + \int^\infty_0 dq*q^2*v(p', q; E)*T(q, p; k)/(k^2 - q^2 + i\epsilon)
  ! with k=sqrt(2\mu E) where \mu is the reduced mass
  !
  ! T is related to the phase shift by
  !        T = -2*e^(i\delta) sin\delta / (\pi k)
  ! and to the S-matrix by
  !        S = 1 - i\pi k T
  !
  ! Inputs:  k,
  ! N: number of mesh pts, excluding T(k, p; k)
  ! meshpts: array with p' mesh points, excluding p' = k
  ! wghts: array with weights for numerical integration
  ! potV: potV(i) = v(meshpts(i), p; E) for i =< N and potV(N+1) = v(k, p; E)
  ! potM: potM(i, j) = v(meshpts(i), meshpts(j); E) for i,j =< N
  !   potM(N+1, j) = v(k, meshpts(j); E) for j =< N
  !   potM(i, N+1) = v(meshpts(i), k; E) for i =< N
  !   potM(N+1, N+1) = v(k, k; E)
  ! Outputs:  tmcolp with tmcolp(i) = T(meshpts(i), p; k) ( i = 1:N)
  !           and tmcolp(N+1) = T(k, p; k)
  !           flag : 0, success; -1, zgesv reports error
  ! To get the on-shell T, solve first the half on-shell T by setting p = k and
  ! then T(k, k; k) = tmcolp
  !
  ! note:
  ! * instead of passing the potential v, the calling function
  !   is asked to fill out potV and potM.
  ! * v is real and/or non-symmetric.
  ! * v is assumed to be regularized properly so that the integral converges.
  ! * The mesh is assumed to properly map onto [0, \infty] (e.g., TanMesh)
  ! * LAPACK routine zgesv is used.
  !
  ! Upon discretization, the integral equation can be written as
  ! T(ii) = potV(ii) + \sum_{jj=1}^N w(jj) * (
  !   q(jj)^2/(k^2-q(jj)^2)*potM(ii, jj)*T(ii) - k^2/(k^2-q(ii)^2)*T(N+1) )
  !   + k^2*potM(ii, N+1)*T(N+1) * \int^\infty_0 dqq^2 1/(k^2-q^2+i\epsilon)
  ! Evaluating the integral in the third line and sorting out the coefficients
  ! of T(ii) for ii = 1 to N+1
  ! T(ii) = potV(ii) + \sum_{jj=1}^{N+1} coefD(jj)*potM(ii, jj)*T(jj)
  ! where
  ! coefD(ii)     =  wghts(ii) * q(ii)^2/(k^2 - q(ii)^2) for ii < N+1
  ! coefD(N+1)    = -k^2 * \sum_{ii=1}^N * wghts(ii)/(k^2 - q(ii)^2) - i\pi/2*k
  !
  ! Continue to rewrite the equation as a linear equation
  ! \sum_{ii=1}^{N+1} coefF(ii, jj)*T(jj) = potV(ii)   for ii = 1 to N+1
  ! where
  ! coefF(ii, jj) = \delta_{ii, jj} - coefD(jj)*potM(ii, jj)

  subroutine lsesv_half_sngl(k, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    integer, intent(in)         :: N
    integer, intent(out)        :: flag
    real(NER), intent(in)       :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:)

    ! loop indexes
    integer                                 :: ii, jj
    ! zgesv error status
    integer                                 :: info
    ! pivot indecies, needed by zgesv
    integer, dimension(1:N+1)               :: ipiv
    ! temporary vars
    real(NER)                               :: ksqr, wokml, regsum
    complex(NEC), dimension(1:N+1)          :: coefD            ! explaind above
    complex(NEC), dimension(1:N+1, 1:N+1)   :: coefF            ! explaind above

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *)       &
        & "lsesv_half_sngl: N = ", N, " is not positive. nothing done"
      flag = lsesv_FAILED
      return
    end if

    ! coefD(ii)     =  wghts(ii) * q(ii)^2/(k^2 - q(ii)^2) for ii <= N
    ! coefD(N+1)    = -k^2 * \sum_{ii=1}^N * wghts(ii)/(k^2 - q(ii)^2) - i\pi/2*k
    ksqr = k * k
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (ksqr - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do
    ! LinMesh
    ! coefD(N+1) = - ksqr*regsum - k * (IMGUNT_NE*PIO2_NE)
    ! TanMesh
    coefD(N+1) = - ksqr*regsum - k * IMGUNT_NE*PIO2_NE

    forall (ii = 1:N+1, jj = 1:N+1)
      coefF(ii, jj) = - coefD(jj) * potM(ii, jj)
    end forall
    ! This acts like the Kronecker delta
    do ii = 1, N+1
      coefF(ii, ii) = coefF(ii, ii) + 1.0_NER
    end do

  ! call zgesv of LAPACK to solve A*X = B. B contains the result, X.
  ! be careful about the meaning of LDA and LDB
  ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    tmcolp(1:N+1) = potM(1:N+1, N+1)
    if (NEC .eq. DPC) then
      call zgesv(N+1, 1, coefF, N+1, ipiv, tmcolp, N+1, info)
    else
      call cgesv(N+1, 1, coefF, N+1, ipiv, tmcolp, N+1, info)
    end if
    if (info .ne. 0) then               ! zgesv reports error
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_half_sngl

  ! solve for full off-shell T-matrix
  subroutine lsesv_full_sngl(k, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    integer, intent(in)         :: N
    integer, intent(out)        :: flag
    real(NER), intent(in)       :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:, :)


    ! loop indexes
    integer                                 :: ii, jj
    ! zgesv error status
    integer                                 :: info
    ! pivot indecies, needed by zgesv
    integer, dimension(1:N+1)               :: ipiv
    ! temporary vars
    real(NER)                               :: ksqr, wokml, regsum
    complex(NEC), dimension(1:N+1)          :: coefD            ! explaind above
    complex(NEC), dimension(1:N+1, 1:N+1)   :: coefF            ! explaind above

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *)       &
        & "lsesv_full_sngl: N = ", N, " is not positive. nothing done"
      flag = lsesv_FAILED
      return
    end if
    if (k <= lsesv_EPS) then
      write(standard_error_unit, *)       &
        & "lsesv_full_sngl: k = ", k, " is not positive. nothing done"
      flag = lsesv_FAILED
      return
    end if

  ! coefD(ii)     =  wghts(ii) * q(ii)^2/(k^2 - q(ii)^2) for ii <= N
  ! coefD(N+1)    = -k^2 * \sum_{ii=1}^N * wghts(ii)/(k^2 - q(ii)^2) - i\pi/2*k
    ksqr = k * k
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (ksqr - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do
    ! LinMesh
    ! coefD(N+1) = - ksqr*regsum - k * (IMGUNT_NE*PIO2_NE)
    ! TanMesh
    coefD(N+1) = - ksqr*regsum - k * IMGUNT_NE*PIO2_NE

    forall (ii = 1:N+1, jj = 1:N+1)
      coefF(ii, jj) = - coefD(jj) * potM(ii, jj)
    end forall
    ! This acts like the Kronecker delta
    do ii = 1, N+1
      coefF(ii, ii) = coefF(ii, ii) + 1.0_NER
    end do

  ! call zgesv of LAPACK to solve A*X = B. B contains the result, X.
  ! be careful about the meaning of LDA and LDB
  ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    tmcolp(1:N+1, 1:N+1) = potM(1:N+1, 1:N+1)
    if (NEC .eq. DPC) then
      call zgesv(N+1, N+1, coefF, N+1, ipiv, tmcolp, N+1, info)
    else
      call cgesv(N+1, N+1, coefF, N+1, ipiv, tmcolp, N+1, info)
    end if
    if (info .ne. 0) then               ! zgesv reports error
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_full_sngl

  ! LS eqn solvers for coupled-channel problem (NCH channels)
  ! Lippmann-Schwinger equation of "T-matrix", T(p', p; k) with +i\epsilon
  ! prescription.
  !
  ! T_n'n(p', p; k) = v_n'n(p', p; E) + \sum_n" \int^\infty_0 dq*q^2*v_n'n"(p', q; E)
  !       *T_n"n(q, p; k)/(k^2 - q^2 + i\epsilon) with k=sqrt(2\mu E) where \mu is the reduced mass
  ! where n, n' and n" = 1, 2 ... NCH are the channel index.
  ! Inputs:
  !       k > 0,
  !       N: number of mesh pts, excluding T(k, p; k)
  !       NCH: number of channels that are coupled
  !       meshpts: array with p' mesh points, excluding p' = k
  !       wghts: array with weights for numerical integration
  ! potM:
  !       potM((a-1)*(N+1) + i, (b-1)*(N+1) + j) =
  !         v_ab(meshpts(i), meshpts(j); E) for i,j =< N
  !       potM(a*(N+1), (b-1)*(N+1)+j) = v_ab(k, meshpts(j); E) for j =< N
  !       potM((a-1)*(N+1)+i, b*(N+1)) = v_ab(meshpts(i), k; E) for i =< N
  !       potM(a*(N+1), b*(N+1)) = v_ab(k, k; E)
  ! Outputs:
  !       tmcolp((a-1)*(N+1) + i, n) = T_an(meshpts(i), p; k) (i = 1:N)
  !         and tmcolp(a*(N+1))  = T_an(k, p; k)
  !       flag : 0, success; 1, zgesv reports wrong
  ! To obtain the on-shell T, solve first the half on-shell T by setting p = k
  ! and then T_an(k, k; k) = tmcolp(a*(N+1), n)
  !
  ! note:
  ! * instead of passing the potential v, the calling function
  !   is asked to fill out potM.
  ! * v is real and/or non-symmetric.
  ! * LAPACK routine zgesv is used.
  subroutine lsesv_nch_half_cpld(NCH, k, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    external zgesv
    external cgesv

    integer, intent(in)         :: NCH, N
    integer, intent(out)        :: flag
    real(NER), intent(in)       :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:, :)

    integer :: ii, jj, aa, bb, pos, posp
    integer :: info
    integer :: ipiv(1:NCH*(N+1))
    real(NER) :: ksqr, wokml, regsum
    complex(NEC)  :: coefD(1:N+1), coefF(1:NCH*(N+1), 1:NCH*(N+1))

    flag = lsesv_SUCCEEDED
    if (N <= 0 .or. NCH <= 0) then
      write(standard_error_unit, *)       &
        & "lsesv_nch_half_cpld: N or NCH <= 0"
      flag = lsesv_FAILED
      return
    end if
    if (k <= lsesv_EPS) then
      write(standard_error_unit, *)       &
        & "lsesv_nch_half_cpld: k = ", k, " is not positive. nothing done"
      flag = lsesv_FAILED
      return
    end if

    ! assemble coefD according to Eq (A.10) in Jerry's thsis
    ksqr = k * k
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (ksqr - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do
    coefD(N+1) = -ksqr*regsum - k*IMGUNT_NE*PIO2_NE

    ! assemble coefF according to Eq (A.14) in Jerry's thsis
    do ii = 1, N+1
      do jj = 1, N+1
        do aa = 1, NCH
          do bb = 1, NCH
            pos = (aa-1)*(N+1) + ii
            posp = (bb-1)*(N+1) + jj
            coefF(pos, posp) = -coefD(jj) * potM(pos, posp)
          end do
        end do
      end do
    end do

    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, N+1
      do aa = 1, NCH
        pos = (aa-1)*(N+1) + ii
        coefF(pos, pos) = coefF(pos, pos) + 1.0_NER
      end do
    end do

    do bb = 1, NCH
      posp = bb*(N+1)
      do ii = 1, N+1
        do aa = 1, NCH
          pos = (aa-1)*(N+1) + ii
          tmcolp(pos, bb) = potM(pos, posp)
        end do
      end do
    end do

    ! call zgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NEC .eq. DPC) then
      call zgesv(NCH*(N+1), NCH, coefF, NCH*(N+1), ipiv, tmcolp, NCH*(N+1), info)
    else
      call cgesv(NCH*(N+1), NCH, coefF, NCH*(N+1), ipiv, tmcolp, NCH*(N+1), info)
    end if
    if (info /= 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nch_half_cpld

  subroutine lsesv_nch_half_cmplx_cpld(NCH, mE, N, potM, Vk, meshpts, wghts,    &
    & qsqr, tmcolp, flag)

    external zgesv
    external cgesv

    integer, intent(in)         :: NCH, N
    integer, intent(out)        :: flag
    complex(NEC), intent(in)    :: mE
    complex(NEC), intent(in)    :: meshpts(:), wghts(:), qsqr(:), potM(:, :), Vk(:, :)
    complex(NEC), intent(out)   :: tmcolp(:, :)

    integer :: ii, jj, aa, bb, pos, posp
    integer :: info
    integer :: ipiv(1:NCH*N)
    complex(NEC)  :: wokml, coefD(1:N), coefF(1:NCH*N, 1:NCH*N)

    flag = lsesv_SUCCEEDED
    if (N <= 0 .or. NCH <= 0) then
      write(standard_error_unit, *) "lsesv_nch_half_cmplx_cpld: N or NCH <= 0"
      flag = lsesv_FAILED
      return
    end if

    ! assemble coefD according to Eq (A.10) in Jerry's thsis
    do ii = 1, N
      wokml = wghts(ii) / (mE - qsqr(ii))
      coefD(ii) =  qsqr(ii) * wokml
    end do

    ! assemble coefF according to Eq (A.14) in Jerry's thsis
    do ii = 1, N
      do aa = 1, NCH
        pos = idx1d_by_ch_p(aa, ii, N)
        do jj = 1, N
          do bb = 1, NCH
            posp = idx1d_by_ch_p(bb, jj, N)
            coefF(pos, posp) = -coefD(jj) * potM(pos, posp)
          end do
        end do
      end do
    end do

    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, N
      do aa = 1, NCH
        pos = idx1d_by_ch_p(aa, ii, N)
        coefF(pos, pos) = coefF(pos, pos) + 1.0_NER
      end do
    end do

    tmcolp(1:N*NCH, 1:NCH) = Vk(1:N*NCH, 1:NCH)

    ! call zgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NEC .eq. DPC) then
      call zgesv(NCH*N, NCH, coefF, NCH*N, ipiv, tmcolp, NCH*N, info)
    else
      call cgesv(NCH*N, NCH, coefF, NCH*N, ipiv, tmcolp, NCH*N, info)
    end if
    if (info /= 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nch_half_cmplx_cpld

  ! T(p_i, k_j) = V(p_i, k_j) + \sum_n V(p_i, q_n) G_0(q_n) T(q_n, k_j)
  ! i,n = 1..Np; j = 1..Nk
  subroutine lsesv_nch_generic_cmplx(NCH, mE, Np, Nk, potM, Vk, meshpts, wghts, &
    & qsqr, tmcolp, flag)

    external zgesv
    external cgesv

    integer, intent(in)         :: NCH, Np, Nk
    integer, intent(out)        :: flag
    complex(NEC), intent(in)    :: mE
    complex(NEC), intent(in)    :: meshpts(:), wghts(:), qsqr(:), potM(:, :), Vk(:, :)
    complex(NEC), intent(out)   :: tmcolp(:, :)

    integer :: ii, jj, aa, bb, pos, posp
    integer :: info
    integer :: ipiv(1:NCH*Np)
    complex(NEC)  :: wokml, coefD(1:Np), coefF(1:NCH*Np, 1:NCH*Np)

    flag = lsesv_SUCCEEDED
    if (Np <= 0 .or. Nk <= 0 .or. NCH <= 0) then
      write(standard_error_unit, *) "lsesv_nch_generic_cmplx: Np, Nk, or NCH <= 0"
      flag = lsesv_FAILED
      return
    end if

    ! assemble coefD according to Eq (A.10) in Jerry's thsis
    do ii = 1, Np
      wokml = wghts(ii) / (mE - qsqr(ii))
      coefD(ii) =  qsqr(ii) * wokml
    end do

    ! assemble coefF according to Eq (A.14) in Jerry's thsis
    do ii = 1, Np
      do aa = 1, NCH
        pos = idx1d_by_ch_p(aa, ii, Np)
        do jj = 1, Np
          do bb = 1, NCH
            posp = idx1d_by_ch_p(bb, jj, Np)
            coefF(pos, posp) = -coefD(jj) * potM(pos, posp)
          end do
        end do
      end do
    end do

    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, Np
      do aa = 1, NCH
        pos = idx1d_by_ch_p(aa, ii, Np)
        coefF(pos, pos) = coefF(pos, pos) + 1.0_NER
      end do
    end do

    tmcolp(1:Np*NCH, 1:Nk*NCH) = Vk(1:Np*NCH, 1:Nk*NCH)

    ! call zgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NEC .eq. DPC) then
      call zgesv(NCH*Np, NCH*Nk, coefF, NCH*Np, ipiv, tmcolp, NCH*Np, info)
    else
      call cgesv(NCH*Np, NCH*Nk, coefF, NCH*Np, ipiv, tmcolp, NCH*Np, info)
    end if
    if (info /= 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nch_generic_cmplx

  ! k > 0
  subroutine lsesv_nch_full_cpld(NCH, k, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    external zgesv
    external cgesv

    integer, intent(in)         :: NCH, N
    integer, intent(out)        :: flag
    real(NER), intent(in)       :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:, :)

    integer :: ii, jj, aa, bb, pos, posp
    integer :: info
    integer :: ipiv(1:NCH*(N+1))
    real(NER) :: ksqr, wokml, regsum
    complex(NEC)  :: coefD(1:N+1), coefF(1:NCH*(N+1), 1:NCH*(N+1))

    flag = lsesv_SUCCEEDED
    if (N <= 0 .or. NCH <= 0) then
      write(standard_error_unit, *) "lsesv_nch_full_cpld: N or NCH <= 0"
      flag = lsesv_FAILED
      return
    end if
    if (k <= lsesv_EPS) then
      write(standard_error_unit, *)       &
        & "lsesv_nch_full_cpld: k = ", k, " is not positive. nothing done"
      flag = lsesv_FAILED
      return
    end if

    ! assemble coefD according to Eq (A.10) in Jerry's thsis
    ksqr = k * k
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (ksqr - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do
    coefD(N+1) = -ksqr*regsum - k*IMGUNT_NE*PIO2_NE

    ! assemble coefF according to Eq (A.14) in Jerry's thsis
    do ii = 1, N+1
      do jj = 1, N+1
        do aa = 1, NCH
          do bb = 1, NCH
            pos = (aa-1)*(N+1) + ii
            posp = (bb-1)*(N+1) + jj
            coefF(pos, posp) = -coefD(jj) * potM(pos, posp)
          end do
        end do
      end do
    end do

    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, N+1
      do aa = 1, NCH
        pos = (aa-1)*(N+1) + ii
        coefF(pos, pos) = coefF(pos, pos) + 1.0_NER
      end do
    end do

    tmcolp(1:NCH*(N+1), 1:NCH*(N+1)) = potM(1:NCH*(N+1), 1:NCH*(N+1))
    ! call zgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NEC .eq. DPC) then
      call zgesv(NCH*(N+1), NCH*(N+1), coefF, NCH*(N+1), ipiv, tmcolp, NCH*(N+1), info)
    else
      call cgesv(NCH*(N+1), NCH*(N+1), coefF, NCH*(N+1), ipiv, tmcolp, NCH*(N+1), info)
    end if
    if (info /= 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nch_full_cpld

  ! Sub-threshold: mE < 0
  ! Note: potM, tmcolp -> (1:N, 1:N)
  subroutine lsesv_nch_full_ST_cpld(NCH, mE, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    external dgesv
    external sgesv

    integer, intent(in)      :: NCH, N
    integer, intent(out)     :: flag
    real(NER), intent(in)    :: mE
    real(NER), intent(in)    :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    real(NER), intent(out)   :: tmcolp(:, :)

    integer :: ii, jj, aa, bb, pos, posp
    integer :: info
    integer :: ipiv(1:NCH*N)
    real(NER) :: wokml, regsum
    real(NER) :: coefD(1:N), coefF(1:NCH*N, 1:NCH*N)

    flag = lsesv_SUCCEEDED
    if (N <= 0 .or. NCH <= 0) then
      write(standard_error_unit, *) "lsesv_nch_full_ST_cpld: N or NCH <= 0"
      flag = lsesv_FAILED
      return
    end if

    if (mE > 0.0_NER) then
      flag = lsesv_FAILED
      write (standard_error_unit, '(a)')  &
        & 'lsesv_nch_full_ST_cpld: mE appears to be positive.'
        flag = lsesv_FAILED
      return
    end if
    ! assemble coefD according to Eq (A.10) in Jerry's thsis
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (mE - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do

    ! assemble coefF according to Eq (A.14) in Jerry's thsis
    do ii = 1, N
      do jj = 1, N
        do aa = 1, NCH
          do bb = 1, NCH
            pos = (aa-1)*N + ii
            posp = (bb-1)*N + jj
            coefF(pos, posp) = -coefD(jj) * potM(pos, posp)
          end do
        end do
      end do
    end do

    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, N
      do aa = 1, NCH
        pos = (aa-1)*N + ii
        coefF(pos, pos) = coefF(pos, pos) + 1.0_NER
      end do
    end do

    tmcolp(1:NCH*N, 1:NCH*N) = potM(1:NCH*N, 1:NCH*N)
    ! call sgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NER .eq. DP) then
      call dgesv(NCH*N, NCH*N, coefF, NCH*N, ipiv, tmcolp, NCH*N, info)
    else
      call sgesv(NCH*N, NCH*N, coefF, NCH*N, ipiv, tmcolp, NCH*N, info)
    end if
    if (info /= 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nch_full_ST_cpld

  subroutine lsesv_nch_full_cmplx_cpld(NCH, mNtimesE, N, potM, meshpts, wghts,  &
    & qsqr, tmcolp, flag)

    external zgesv
    external cgesv

    integer, intent(in)         :: NCH, N
    integer, intent(out)        :: flag
    complex(NEC), intent(in)    :: mNtimesE
    complex(NEC), intent(in)    :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:, :)

    integer :: ii, jj, aa, bb, pos, posp
    integer :: info
    integer :: ipiv(1:NCH*N)
    complex(NEC) :: wokml
    ! , regsum
    complex(NEC)  :: coefD(1:N), coefF(1:NCH*N, 1:NCH*N)

    flag = lsesv_SUCCEEDED
    if (N <= 0 .or. NCH <= 0) then
      write(standard_error_unit, *) "lsesv_nch_full_cmplx_cpld: N or NCH <= 0"
      flag = lsesv_FAILED
      return
    end if
    ! assemble coefD according to Eq (A.10) in Jerry's thsis
    ! regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (mNtimesE - qsqr(ii))
      ! regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do

    ! assemble coefF according to Eq (A.14) in Jerry's thsis
    do ii = 1, N
      do jj = 1, N
        do aa = 1, NCH
          do bb = 1, NCH
            pos = (aa-1)*N + ii
            posp = (bb-1)*N + jj
            coefF(pos, posp) = -coefD(jj) * potM(pos, posp)
          end do
        end do
      end do
    end do

    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, N
      do aa = 1, NCH
        pos = (aa-1)*N + ii
        coefF(pos, pos) = coefF(pos, pos) + 1.0_NER
      end do
    end do

    tmcolp(1:NCH*N, 1:NCH*N) = potM(1:NCH*N, 1:NCH*N)
    ! call zgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NEC .eq. DPC) then
      call zgesv(NCH*N, NCH*N, coefF, NCH*N, ipiv, tmcolp, NCH*N, info)
    else
      call cgesv(NCH*N, NCH*N, coefF, NCH*N, ipiv, tmcolp, NCH*N, info)
    end if
    if (info /= 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nch_full_cmplx_cpld

  ! /-------------------------------------------------------------/
  ! Below are less commonly used subs, perhaps not stable
  ! /-------------------------------------------------------------/

  ! 2-channel LSE solver
  ! assuming the regulator is a sharp cutoff, \theta(\Lambda - qii),
  ! where \theta is the step function
  ! tmcolp(ii, n, n') = T_{n n'}(pii, k)
  subroutine lsesv_2ch_half_cpld(k, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    external zgesv
    external cgesv

    integer, intent(in)         :: N
    integer, intent(out)        :: flag
    real(NER), intent(in)       :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :, :, :)
    complex(NEC), intent(out)   :: tmcolp(:, :, :)

    ! loop indecies
    integer :: ii, jj
    ! zgesv error status
    integer :: info
    ! pivot indecies, needed by zgesv
    integer, dimension(1:2*(N+1)) :: ipiv
    ! temporary vars
    real(NER) :: ksqr, wokml, regsum
    ! coefD: Coefficient vector D defined in Appendix A of Jerry's thesis
    complex(NEC), dimension(1:N+1) :: coefD
    ! coefF: Coefficient matrix F defined in Appendix A of Jerry's thesis
    complex(NEC), dimension(1:2, 1:2, 1:N+1, 1:N+1) :: coefF
    complex(NEC), dimension(1:2*(N+1), 1:2*(N+1)) :: reshaped_coefF
    complex(NEC), dimension(1:2*(N+1), 1:2) :: reshaped_tmcolp

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *) "lsesv_2ch_half_cpld: N <= 0"
      flag = lsesv_FAILED
      return
    end if
    if (k <= lsesv_EPS) then
      write(standard_error_unit, *)       &
        & "lsesv_2ch_half_cpld: k = ", k, " is not positive. nothing done"
      flag = lsesv_FAILED
      return
    end if

    ! assemble coefD according to Eq (A.10) in Jerry's thsis
    ksqr = k * k
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (ksqr - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do
    ! coefD(N+1) = -ksqr*regsum - k*(IMGUNT_NE*PIO2_NE)
    coefD(N+1) = -ksqr*regsum - k*IMGUNT_NE*PIO2_NE

    ! assemble coefF according to Eq (A.14) in Jerry's thsis
    forall (ii = 1:N+1, jj = 1:N+1)
      coefF(1:2, 1:2, ii, jj) = - coefD(jj) * potM(1:2, 1:2, ii, jj)
    end forall
    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, N+1
      coefF(1, 1, ii, ii) = coefF(1, 1, ii, ii) + 1.0_NER
      coefF(2, 2, ii, ii) = coefF(2, 2, ii, ii) + 1.0_NER
    end do

    reshaped_coefF(1:N+1, 1:N+1) = coefF(1, 1, 1:N+1, 1:N+1)
    reshaped_coefF(1:N+1, N+2:2*N+2) = coefF(1, 2, 1:N+1, 1:N+1)
    reshaped_coefF(N+2:2*N+2, 1:N+1) = coefF(2, 1, 1:N+1, 1:N+1)
    reshaped_coefF(N+2:2*N+2, N+2:2*N+2) = coefF(2, 2, 1:N+1, 1:N+1)

    forall (ii = 1:N+1, jj = 1:2)
      reshaped_tmcolp(ii, jj)     = potM(1, jj, ii, N+1)
      reshaped_tmcolp(N+1+ii, jj) = potM(2, jj, ii, N+1)
    end forall

    ! call zgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NEC .eq. DPC) then
      call zgesv(2*N+2, 2, reshaped_coefF, 2*N+2, ipiv, reshaped_tmcolp, 2*N+2, info)
    else
      call cgesv(2*N+2, 2, reshaped_coefF, 2*N+2, ipiv, reshaped_tmcolp, 2*N+2, info)
    end if
    forall (ii = 1:N+1, jj = 1:2)
      tmcolp(1, jj, ii) = reshaped_tmcolp(ii, jj)
      tmcolp(2, jj, ii) = reshaped_tmcolp(N+1+ii, jj)
    end forall
    if (info .ne. 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_2ch_half_cpld

  ! tmcolp(ii, jj, n, n') = T_{n n'}(pii, pjj)
  subroutine lsesv_2ch_full_cpld(k, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    external zgesv
    external cgesv

    integer, intent(in)         :: N
    integer, intent(out)        :: flag
    real(NER), intent(in)       :: k
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :, :, :)
    complex(NEC), intent(out)   :: tmcolp(:, :, :, :)

    ! loop indecies
    integer :: ii, jj
    ! zgesv error status
    integer :: info
    ! pivot indecies, needed by zgesv
    integer, dimension(1:2*(N+1)) :: ipiv
    ! temporary vars
    real(NER) :: ksqr, wokml, regsum
    ! coefD: Coefficient vector D defined in Appendix A of Jerry's thesis
    complex(NEC), dimension(1:N+1) :: coefD
    ! coefF: Coefficient matrix F defined in Appendix A of Jerry's thesis
    complex(NEC), dimension(1:2, 1:2, 1:N+1, 1:N+1) :: coefF
    complex(NEC), dimension(1:2*(N+1), 1:2*(N+1)) :: reshaped_coefF
    complex(NEC), dimension(1:2*(N+1), 1:2*(N+1)) :: reshaped_tmcolp

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *) "lsesv_2ch_full_cpld: N <= 0"
      flag = lsesv_FAILED
      return
    end if
    if (k <= lsesv_EPS) then
      write(standard_error_unit, *)       &
        & "lsesv_2ch_full_cpld: k = ", k, " is not positive. nothing done"
      flag = lsesv_FAILED
      return
    end if

    ! assmble coefD according to Eq (A.10) in Jerry's thsis
    ksqr = k * k
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (ksqr - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do
    ! coefD(N+1) = -ksqr*regsum - k*(IMGUNT_NE*PIO2_NE)
    coefD(N+1) = -ksqr*regsum - k*IMGUNT_NE*PIO2_NE

    ! assmble coefF according to Eq (A.14) in Jerry's thsis
    forall (ii = 1:N+1, jj = 1:N+1)
      coefF(1:2, 1:2, ii, jj) = - coefD(jj) * potM(1:2, 1:2, ii, jj)
    end forall
    ! This acts like the Kronecker delta in (A.14)
    do ii = 1, N+1
      coefF(1, 1, ii, ii) = coefF(1, 1, ii, ii) + 1.0_NER
      coefF(2, 2, ii, ii) = coefF(2, 2, ii, ii) + 1.0_NER
    end do

    reshaped_coefF(1:N+1, 1:N+1)         = coefF(1, 1, 1:N+1, 1:N+1)
    reshaped_coefF(1:N+1, N+2:2*N+2)     = coefF(1, 2, 1:N+1, 1:N+1)
    reshaped_coefF(N+2:2*N+2, 1:N+1)     = coefF(2, 1, 1:N+1, 1:N+1)
    reshaped_coefF(N+2:2*N+2, N+2:2*N+2) = coefF(2, 2, 1:N+1, 1:N+1)

    reshaped_tmcolp(1:N+1, 1:N+1)         = potM(1, 1, 1:N+1, 1:N+1)
    reshaped_tmcolp(1:N+1, N+2:2*N+2)     = potM(1, 2, 1:N+1, 1:N+1)
    reshaped_tmcolp(N+2:2*N+2, 1:N+1)     = potM(2, 1, 1:N+1, 1:N+1)
    reshaped_tmcolp(N+2:2*N+2, N+2:2*N+2) = potM(2, 2, 1:N+1, 1:N+1)

    ! call zgesv to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    if (NEC .eq. DPC) then
      call zgesv(2*N+2, 2*N+2, reshaped_coefF, 2*N+2, ipiv, reshaped_tmcolp, 2*N+2, info)
    else
      call cgesv(2*N+2, 2*N+2, reshaped_coefF, 2*N+2, ipiv, reshaped_tmcolp, 2*N+2, info)
    end if
    tmcolp(1, 1, 1:N+1, 1:N+1) = reshaped_tmcolp(1:N+1, 1:N+1)
    tmcolp(1, 2, 1:N+1, 1:N+1) = reshaped_tmcolp(1:N+1, N+2:2*N+2)
    tmcolp(2, 1, 1:N+1, 1:N+1) = reshaped_tmcolp(N+2:2*N+2, 1:N+1)
    tmcolp(2, 2, 1:N+1, 1:N+1) = reshaped_tmcolp(N+2:2*N+2, N+2:2*N+2)
    if (info .ne. 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_2ch_full_cpld

  !   T(p', p; E) where E is complex and is on the physical plane, while p' and p are still real and > 0
  !   When mE is real, it must be negative.

  subroutine lsesv_full_cmplxE_sngl(mE, Lambda, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    complex(NEC), intent(in)    :: mE
    real(NER), intent(in)       :: Lambda
    integer, intent(in)         :: N
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:, :)
    integer, intent(out)        :: flag

    ! loop indexes
    integer                                 :: ii, jj
    ! zgesv error status
    integer                                 :: info
    ! pivot indecies, needed by zgesv
    integer, dimension(1:N)                 :: ipiv
    ! temporary vars
    complex(NEC)                            :: wokml
    ! , regsum
    complex(NEC), dimension(1:N)            :: coefD            ! explaind above
    complex(NEC), dimension(1:N, 1:N)       :: coefF            ! explaind above

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *) "lsesv_full_cmplxE_sngl: N <= 0"
      flag = lsesv_FAILED
      return
    end if
    ! coefD(ii)     =  wghts(ii) * q(ii)^2/(mE - q(ii)^2) for ii <= N
    ! mE must be imaginary or off the positive half of the real axis
    if (aimag(mE/Lambda) .lt. 1.0E-10_NER .and. real(mE) .gt. 0.0_NER) then
      flag = lsesv_FAILED
      return
    end if
    ! regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (mE - qsqr(ii))
      ! regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do

    forall (ii = 1:N, jj = 1:N)
      coefF(ii, jj) = - coefD(jj) * potM(ii, jj)
    end forall
    ! This acts like the Kronecker delta
    do ii = 1, N+1
      coefF(ii, ii) = coefF(ii, ii) + 1.0_NER
    end do

    ! call zgesv of LAPACK to solve A*X = B. B contains the result, X.
    ! be careful about the meaning of LDA and LDB
    ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    tmcolp(1:N+1, 1:N+1) = potM(1:N+1, 1:N+1)
    if (NEC .eq. DPC) then
      call zgesv(N, N, coefF, N, ipiv, tmcolp, N+1, info)
    else
      call cgesv(N, N, coefF, N, ipiv, tmcolp, N+1, info)
    end if
    if (info .ne. 0) then               ! zgesv reports error
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_full_cmplxE_sngl


  ! Lippmann-Schwinger equation of "T-matrix", T(p', p; k) with +i\epsilon
  ! prescription and SHARP cutoff
  !
  ! T(p', p; k) = v(p', p; E) + \int^\Lambda_0 dq*q^2*v(p', q; E)*T(q, p; k)/(k^2 - q^2 + i\epsilon)
  ! with k=sqrt(2\mu E) where \mu is the reduced mass
  subroutine lsesv_half_sngl_sharp(k, Lambda, N, potM, meshpts, wghts, qsqr, tmcolp, flag)

    integer, intent(in)         :: N
    integer, intent(out)        :: flag
    real(NER), intent(in)       :: k, Lambda
    real(NER), intent(in)       :: meshpts(:), wghts(:), qsqr(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:)

    ! loop indexes
    integer                                 :: ii, jj
    ! zgesv error status
    integer                                 :: info
    ! pivot indecies, needed by zgesv
    integer, dimension(1:N+1)               :: ipiv
    ! temporary vars
    real(NER)                               :: ksqr, wokml, regsum
    complex(NEC), dimension(1:N+1)          :: coefD            ! explaind above
    complex(NEC), dimension(1:N+1, 1:N+1)   :: coefF            ! explaind above

    flag = lsesv_SUCCEEDED

    ! coefD(ii)     =  wghts(ii) * q(ii)^2/(k^2 - q(ii)^2) for ii <= N
    ! coefD(N+1)    = -k^2 * \sum_{ii=1}^N * wghts(ii)/(k^2 - q(ii)^2) - i\pi/2*k + k/2*ln((Lambda+k)/(Lambda-k))
    ksqr = k * k
    regsum = 0.0_NER
    do ii = 1, N
      wokml = wghts(ii) / (ksqr - qsqr(ii))
      regsum = regsum + wokml
      coefD(ii) =  qsqr(ii) * wokml
    end do
    coefD(N+1) = - ksqr*regsum - k * (IMGUNT_NE*PIO2_NE - 0.5_NER*log((Lambda+k)/(Lambda-k)))

    forall (ii = 1:N+1, jj = 1:N+1)
      coefF(ii, jj) = - coefD(jj) * potM(ii, jj)
    end forall
    ! This acts like the Kronecker delta
    do ii = 1, N+1
      coefF(ii, ii) = coefF(ii, ii) + 1.0_NER
    end do

  ! call zgesv of LAPACK to solve A*X = B. B contains the result, X.
  ! be careful about the meaning of LDA and LDB
  ! ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    tmcolp(1:N+1) = potM(1:N+1, N+1)
    if (NEC .eq. DPC) then
      call zgesv(N+1, 1, coefF, N+1, ipiv, tmcolp, N+1, info)
    else
      call cgesv(N+1, 1, coefF, N+1, ipiv, tmcolp, N+1, info)
    end if
    if (info .ne. 0) then               ! zgesv reports error
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_half_sngl_sharp


  ! T(q) = V(q) + \int^Lambda_0 dl M(q, l) T(l), where q > 0 and M(q, l) does
  ! not have singularity for l > 0
  subroutine lsesv_nonsing_half_sngl(N, V, potM, wghts, tmcolp, flag)

    integer, intent(in)     :: N
    integer, intent(out)    :: flag
    real(NER), intent(in)   :: wghts(:), V(:), potM(:, :)
    real(NER), intent(out)  :: tmcolp(:)

    integer     :: ii, jj, info, ipiv(1:N)
    real(NER)   :: coefF(1:N, 1:N)

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *) "lsesv_nonsing_half_sngl: N <= 0"
      flag = lsesv_FAILED
      return
    end if

    forall (ii = 1:N, jj = 1:N)
      coefF(ii, jj) = - wghts(jj) * potM(ii, jj)
    end forall
    do ii = 1, N
      coefF(ii, ii) = coefF(ii, ii) + 1.0_NER
    end do
    tmcolp(1:N) = V(1:N)
    if (NER .eq. DP) then
      call dgesv(N, 1, coefF, N, ipiv, tmcolp, N, info)
    else
      call sgesv(N, 1, coefF, N, ipiv, tmcolp, N, info)
    end if
    if (info .ne. 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nonsing_half_sngl


  ! T(q) = V(q) + \int_{Cl} dl M(q, l) T(l), where M(q, l) does
  ! not have singularity on the contour Cl
  subroutine lsesv_nonsing_half_cmplx_sngl(N, V, potM, wghts, tmcolp, flag)

    integer, intent(in)         :: N
    integer, intent(out)        :: flag
    complex(NEC), intent(in)    :: wghts(:), V(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:)

    integer         :: ii, jj, info, ipiv(1:N)
    complex(NEC)    :: coefF(1:N, 1:N)

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *) "lsesv_nonsing_half_cmplx_sngl: N <= 0"
      flag = lsesv_FAILED
      return
    end if

    forall (ii = 1:N, jj = 1:N)
      coefF(ii, jj) = - wghts(jj) * potM(ii, jj)
    end forall
    do ii = 1, N
      coefF(ii, ii) = coefF(ii, ii) + 1.0_NER
    end do
    tmcolp(1:N) = V(1:N)
    if (NEC .eq. DPC) then
      call zgesv(N, 1, coefF, N, ipiv, tmcolp, N, info)
    else
      call cgesv(N, 1, coefF, N, ipiv, tmcolp, N, info)
    end if
    if (info .ne. 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_nonsing_half_cmplx_sngl


  ! T_i = V_i + \sum_j M_ij T_j
  subroutine lsesv_barebone_half_cmplx_sngl(N, V, potM, tmcolp, flag)

    integer, intent(in)         :: N
    integer, intent(out)        :: flag
    complex(NEC), intent(in)    :: V(:), potM(:, :)
    complex(NEC), intent(out)   :: tmcolp(:)

    integer         :: ii, jj, info, ipiv(1:N)
    complex(NEC)    :: coefF(1:N, 1:N)

    flag = lsesv_SUCCEEDED
    if (N <= 0) then
      write(standard_error_unit, *) "lsesv_barebone_half_cmplx_sngl: N <= 0"
      flag = lsesv_FAILED
      return
    end if

    forall (ii = 1:N, jj = 1:N)
      coefF(ii, jj) = - potM(ii, jj)
    end forall
    do ii = 1, N
      coefF(ii, ii) = coefF(ii, ii) + 1.0_NER
    end do
    tmcolp(1:N) = V(1:N)
    if (NEC .eq. DPC) then
      call zgesv(N, 1, coefF, N, ipiv, tmcolp, N, info)
    else
      call cgesv(N, 1, coefF, N, ipiv, tmcolp, N, info)
    end if
    if (info .ne. 0) then
      flag = lsesv_FAILED
    end if

  end subroutine lsesv_barebone_half_cmplx_sngl

end module nneft_lsesv
