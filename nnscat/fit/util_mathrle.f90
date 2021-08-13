! Ying Mei before 06/18/2017
! Bingwei Long 08/15/2017

module util_mathrle

	use nneft_type
	implicit none

contains

	function plgndr_s(l,m,x)

		integer, intent(in)  		:: l, m
		real(NER), intent(in) 	  	:: x

		real(NER) 		:: plgndr_s
		integer		 	:: ll
		real(NER) 		:: pll,pmm,pmmp1,somx2

		pmm=1.0_NER
		if (m >= 0 .and. l >=m .and. abs(x) <= 1.0_NER) then
			if (m > 0) then
				somx2 = sqrt((1.0_NER-x)*(1.0_NER+x))
				pmm = product(arth(1.0_NER,2.0_NER,m))*somx2**m
				if (mod(m,2) == 1) pmm =- pmm
			end if

			if (l == m) then
				plgndr_s = pmm
			else
				pmmp1 = x*(2.0_NER*m+1.0_NER)*pmm
				if (l == m+1) then
					plgndr_s = pmmp1
				else
					do ll = m+2,l
						pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
						pmm = pmmp1
						pmmp1 = pll
					end do
					plgndr_s = pll
				end if
			end if
		else if(l <= abs(m) .and. abs(x)>1.0_NER) then
			plgndr_s = 0
		end if
		return
	end function plgndr_s

	function spherhms(l, m, x)
		integer, intent(in)  		:: l, m
		real(NER), intent(in) 	  	:: x

		integer 				:: i
		real(NER)				:: factorial(0:99), spherhms


		factorial(0) = 1.0

		do I = 1, 99
			factorial(I) = I*factorial(I-1)
		end do


		if (m >= 0 .and. m <= l) then
			spherhms = sqrt((2*l+1)*factorial(l-m)/((4*PI_NE)*factorial(l+m))) * plgndr_s(l, m, x)
		else if (m < 0 .and. abs(m) <= l) then
			spherhms = (-1)**ABS(m) * sqrt((2*l+1)*factorial(l-abs(m))/((4*PI_NE)*factorial(l+abs(m)))) * plgndr_s(l, abs(m), x)
		else
			spherhms = 0
		end if
	end function spherhms

	function dirac_func(l,ls)
		integer, intent(in) 		:: l, ls
		integer    					:: dirac_func

		if (l == ls) then
			dirac_func = 1
		else
			dirac_func = 0
		end if

		return
	end function dirac_func


	function CG_cff(J1,J2,J3,M1,M2)

		integer,intent(in)		:: J1, J2, J3, M1, M2
		real(NER)				:: CG_cff, sumk, term, FACT(0:99)

		integer					:: I, K, M3

		FACT(0) = 1.0

		do I = 1, 99
			FACT(I) = I*FACT(I-1)
		end do

		M3 = M1 + M2

		if	((J3 .LT. ABS(J1-J2)) .OR. &
            &(J3 .GT. (J1+J2))    .OR.  &
            &(ABS(M1) .GT. J1)    .OR.  &
            &(ABS(M2) .GT. J2)    .OR.  &
            &(ABS(M3) .GT. J3)) then
         	CG_cff = 0.0
      	ELSE

			CG_cff = SQRT((J3+J3+1)/FACT(J1+J2+J3+1))
			CG_cff = CG_cff * SQRT(FACT(J1+J2-J3)*FACT(J2+J3-J1)*FACT(J3+J1-J2))
			CG_cff = CG_cff * SQRT(FACT(J1+M1)*FACT(J1-M1)*FACT(J2+M2)*FACT(J2-M2)*FACT(J3+M3)*FACT(J3-M3))
			sumk = 0.0
			do K = 0, 99
				if (J1+J2-J3-K .LT. 0.0D0) CYCLE
        	    if (J3-J1-M2+K .LT. 0.0D0) CYCLE
        	    if (J3-J2+M1+K .LT. 0.0D0) CYCLE
        	    if (J1-M1-K    .LT. 0.0D0) CYCLE
        	    if (J2+M2-K    .LT. 0.0D0) CYCLE
        	    term = FACT(J1+J2-J3-K)*FACT(J3-J1-M2+K)*FACT(J3-J2+M1+K)*FACT(J1-M1-K)*  &
        	       FACT(J2+M2-K)*FACT(K)
        	    if (MOD(K,2) .EQ. 1) term = -term
        	    sumk = sumk + 1.0D0/term
        	 end do
        	 CG_cff = CG_cff * sumk
      end if
  end function CG_cff

	! function deg_rad(x)

	! 	real(NER),intent(in)	:: x
	! 	real(NER)				:: deg_rad

	! 	deg_rad = x/180 * PI_NE

	! 	return
	! end function deg_rad

	real(NER) function deg2rad(x)

		real(NER),intent(in)	:: x
		real(NER)				:: deg_rad

		deg2rad = x/180 * PI_NE

	end function deg2rad

	FUNCTION arth(first,increment,n)
	REAL(DP), INTENT(IN) :: first,increment
	INTEGER(I4B), INTENT(IN) :: n
	REAL(DP), DIMENSION(n) :: arth
	INTEGER(I4B) :: k,k2
	REAL(DP) :: temp
	if (n > 0) arth(1)=first
	if (n <= 16) then
		do k=2,n
			arth(k)=arth(k-1)+increment
		end do
	else
		do k=2,8
			arth(k)=arth(k-1)+increment
		end do
		temp=increment*8
		k=16
		do
			if (k >= n) exit
			k2=k+k
			arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
	END FUNCTION arth


end module util_mathrle
