! Given the one_letter amino acid name, this function spits out the corresponding maximum RSA value taken from Tien et al. paper.
! AMIR SHAHMORADI, SATURDAY 9:56 PM, OCT 26, 2013, ICMB, UT AUSTIN

function MRSA_finder(aa)
use aa_dictionary
implicit none
integer :: i
character(len=1), intent(in) :: aa
real*8 :: MRSA_finder
do i=1,aa_num
	if (aa == aa_dict(i)) then
		MRSA_finder = max_rsa(i)
		return
	end if
end do
write(*,*) 'The Amino Acid letter does not exist in the Amino Acid dictionary!'
write(*,*) aa
stop
end function MRSA_finder
