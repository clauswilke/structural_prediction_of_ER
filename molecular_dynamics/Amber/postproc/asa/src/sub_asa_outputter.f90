! THIS SUBROUTINE TAKES IN THE NAME OF A DSSP FILE (dssp_out) AND WRITES THE CALCULATED RSAs FOR THE CORRESPONDING PROTEIN AMINO ACIDS INTO A FILE WHICH ITS UNIT IS GIVEN BY rsa_file_unit.
! AMIR SHAHMORADI, SATURDAY 4:00 PM, OCT 26, 2013, ICMB, UT AUSTIN

subroutine asa_outputter(dssp_out,asa_file_unit,nres,frame)
implicit none
integer, intent(in) :: asa_file_unit,nres,frame	! WHICH frame OF THE MD TRAJECTORY IS INPUT.
character(len=300), intent(in) :: dssp_out
integer :: i,ios,asa_int
integer, dimension(nres) :: res_num
real*8, dimension(nres) :: asa
character(len=2), dimension(nres) :: chain_id,res_name
character(len=20) :: dummy_char
character(len=100) :: asa_output_format,res_output_format,resnum_output_format
!real*8 :: MRSA_finder		! function

! NOW PREPARE THE ASA OUTPUT FILE.

!write(*,*) 'asa_file_unit: ',asa_file_unit	
write(dummy_char,*) nres
asa_output_format = trim(adjustl('(1I15,' // trim(adjustl(dummy_char)) // 'E15.5)'))
resnum_output_format = trim(adjustl('(1A15,' // trim(adjustl(dummy_char)) // 'I15)'))
res_output_format = trim(adjustl('(1A15,' // trim(adjustl(dummy_char)) // 'A15)'))
!write(*,*) asa_output_format,res_output_format

! NOW READ THE ASA VALUES FROM DSSP OUTPUT FILE.

open (unit=999,file=trim(adjustl(dssp_out)),status='old')

do i = 1,25; read(999,*); end do

read(999,'(2I5,2A2,1A20,1I4)') res_num(1),res_num(1),chain_id(1),res_name(1),dummy_char,asa_int;  asa(1) = real(asa_int)
if (res_num(1)/=1) then; write(*,*) 'Sum Tin Wong! : in asa_outputter.f90, res_num(1) /= 1'; stop; end if

! NORMALIZE THE ASA VALUES.
!asa(1) = asa(1)/MRSA_finder(trim(adjustl(res_name(1))))

!if (rsa(i)>1.5d0) then
!	write(*,*) 'RSA for residue number',res_num(i),'is greater than 1: ',rsa(i)
!	!write(*,*) 'Press Enter to continue...'
!	!read(*,*)
!end if

do i=2,nres
	do
		read(999,'(2I5,2A2,1A20,1I4)') res_num(i),res_num(i),chain_id(i),res_name(i),dummy_char,asa_int;  asa(i) = real(asa_int)
		if (trim(adjustl(res_name(i)))=='!') then
			!write(*,*) 'TER'
			cycle
		end if
		exit
	end do
	!if (trim(adjustl(res_name(i)))/='!') rsa(i) = rsa(i)/MRSA_finder(trim(adjustl(res_name(i))))
	!if (rsa(i)>1.5d0) then
	!	write(*,*) 'RSA for residue number',res_num(i),'is greater than 1: ',rsa(i)
	!	!write(*,*) 'Press Enter to continue...'
	!	!read(*,*)
	!end if
end do
write(*,*) res_name(nres)

! NOW WRITE OUT THE ASA FILE.

if (frame==0) then
	!if (mod(frame,100)==0) write(*,'(1A6,1I6)') 'frame:',frame
	write(asa_file_unit,res_output_format) 'res_name',(res_name(i),i=1,nres)
	write(asa_file_unit,resnum_output_format) 'res_num',(res_num(i),i=1,nres)
	write(asa_file_unit,res_output_format) 'chain_id',(chain_id(i),i=1,nres)
	write(asa_file_unit,asa_output_format) frame,(asa(i),i=1,nres)
else
	!if (mod(frame,100)==0) write(*,'(1A6,1I6)') 'frame:',frame
	write(asa_file_unit,asa_output_format) frame,(asa(i),i=1,nres)
end if

close(999)

end subroutine asa_outputter
