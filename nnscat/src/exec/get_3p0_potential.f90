!get 3p0 potential  pray

program get_3p0_potential
	use potwrap
	implicit none
	integer, parameter  :: N = 400, regtype = REGTYPE_GAUSSIAN
	integer         :: L, S, J, ii, jj
	real(NER)       :: mambda, p1(1:N), p2(1:N), para(1:2), max_p
	real(NER)       :: potval(1:2,1:2), potential1(1:N,1:N), potential2(1:N,1:N), potential3(1:N,1:N)
  real(NER)       :: cco(1:N,1:N), cpoco(1:N,1:N)

  L = 1
  S = 1
  J = 0
  mambda = 800.0_NER
  para(1) = 0.295E-007_NER
  max_p = 400.0_NER


  open(unit=57, file = "ope_potential")
  open(unit=257,file = "all_potential")
  open(unit=157,file = "dig_potential")

  do ii = 1, N
  	p1(ii) = ii * max_p/N
  	do jj = 1, N
      p2(jj) = jj * max_p/N
      call OPE_epwrap(L, S, J, regtype, mambda, para, p1(ii), p2(jj), potval)
      potential1(ii, jj) = potval(1, 1)
      call VLO_withc0_epwrap(L, S, J, regtype, mambda, para, p1(ii), p2(jj), potval)
      potential2(ii, jj) = potval(1, 1)
      call VPNQ0_epwrap(L, S, j, regtype, Mambda, para, p1(ii), p2(jj), potval)
      potential3(ii, jj) = para(1) * potval(1, 1)
      write(57,*) p1(ii), p2(jj), potential1(ii, jj)
  		write(257,*) p1(ii), p2(jj), potential2(ii, jj)
  		if (ii == jj) then
  		  cco(ii, jj) = abs( potential3(ii, jj) / potential1(ii, jj) * 100.0_NER )
  		  cpoco(ii, jj)  = abs( potential2(ii, jj) / potential1(ii, jj) * 100.0_NER )
  			write(157,*) p1(ii), potential1(ii,jj), potential2(ii, jj), potential3(ii, jj), cco(ii, jj), cpoco(ii, jj)
  		end if
    end do
  end do
  close(57)
  close(257)
  close(157)
end program