	SUBROUTINE avevar(datain,n,ave,var)
	IMPLICIT NONE
	integer n,j
	double precision, INTENT(in) :: datain(n)
	double precision, INTENT(OUT) :: ave,var
!	Given array datain, returns its mean as ave and its variance as var.
!	Amir Shahmoradi, June 5, 2009, 10:26 PM, MTU
	double precision s, ep
	ave = 0
	do j = 1,n
		ave = ave + datain(j)
	end do
!	ave=sum(datain(:))/n
	ave=ave/dble(n)
	var = 0
	ep = 0
	do j = 1,n
		s = datain(j) - ave
		ep = ep + s
		var = var + s*s
	end do
!	s(:)=datain(:)-ave
!	var=dot_product(s,s)
	var = (var - (ep**2)/(dble(n)))/(dble(n)-1.0d0) 	!   Corrected two-pass formula (14.1.8).
	END SUBROUTINE avevar