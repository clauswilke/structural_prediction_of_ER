! Amir Shahmoradi, Friday 7:20 PM, October 25, 2013, WilkeLab, UT Austin.
! This version of the code only recognizes Na+ (or SOD) and Cl- (or CLA) as the ions in the pdb file. The pdb file 

!module aa_dictionary
!implicit none
!integer, parameter :: aa_num = 20
!character(len=1), dimension(aa_num) :: aa_dict = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
!real*8, dimension(aa_num) :: max_rsa = [129.,274.,195.,193.,167.,225.,223.,104.,224.,197.,201.,236.,224.,240.,159.,155.,172.,285.,263.,174.]
!end module aa_dictionary

PROGRAM dynamic_asa
IMPLICIT NONE
INTEGER :: i,j,ios=0,frame=0
INTEGER :: natoms=0,nions=0,nres	    ! total number of protein atoms and ions and residues in the Amber pdb file
CHARACTER(LEN=300) :: pdb_in			! pdb file must be free of any non amino acid atoms, except Na+ and Cl- ions that are automatically handled by the code.
CHARACTER(LEN=300) :: crd_in			! The correspodning Amber mdcrd file for the given pdb file .
CHARACTER(LEN=300) :: pdb_out			! pdb file must be free of any non amino acid atoms, except Na+ and Cl- ions that are automatically handled by the code.
CHARACTER(LEN=300) :: asa_out			! The output file that contains the asa information from MD trajectory for the given pdb structure.
CHARACTER(LEN=300) :: command			! The command to be executed in shell.
CHARACTER(LEN=6) :: record
CHARACTER(LEN=4) :: atom_name
CHARACTER(LEN=3) :: res_name
CHARACTER(LEN=1) :: chain_id,alt_loc_ind,res_code
INTEGER :: atom_num,res_num
DOUBLE PRECISION, DIMENSION(3) :: crd
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mdcrd
DOUBLE PRECISION :: occupancy,bfactor,dummy
DOUBLE PRECISION, dimension(3) :: box
LOGICAL :: protein_stat=.true.
!CHARACTER(LEN=8) :: temp_pdb = 'temp.pdb', dssp_out = 'dssp.out'
CHARACTER(LEN=15) :: temp_pdb
CHARACTER(LEN=16) :: dssp_out
INTEGER :: asa_file_unit=22
logical :: dssp_file_stat=.true., pdb_file_stat=.true.
character(len=6) :: dummy_string


call random_seed()

if (command_argument_count()/=4) then
	write(*,*)
	write(*,*) "Incorrect number of input arguments on the command line."
	write(*,*) "Correct use:"
	write(*,*) "./a.out <input pdb file> <input mdcrd file> <output pdb file> <output ASA file>"
	write(*,*)
	stop
end if
call get_command_argument(1,pdb_in)
call get_command_argument(2,crd_in)
call get_command_argument(3,pdb_out)
call get_command_argument(4,asa_out)


open(unit=11,file=trim(adjustl(pdb_in)),status='old')
open(unit=12,file=trim(adjustl(crd_in)),status='old')
open(unit=21,file=trim(adjustl(pdb_out)),status='replace')
open(unit=asa_file_unit,file=trim(adjustl(asa_out)),status='replace')

! FIRST DETERMINE THEW NUMBER OF ATOMS IN THE PDB FILE (EXCLUDING THE IONS: ONLY Na+ & Cl- ARE CONSIDERED AT THE MOMENT).

do
	read(11,'(1A4)') record
	if (trim(adjustl(record))=='ATOM') exit
	!write(*,'(1A4)') record
	cycle
end do

backspace(unit=11)

do
	read(11,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',iostat=ios) record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,(crd(j),j=1,3),occupancy,bfactor
	if (ios<0) then
		exit
	elseif (ios>0) then
		write(*,*) 'Sum Tin Wong! : PDB file broken.'; stop
	else
		if (trim(adjustl(atom_name))=='Na+' .or. trim(adjustl(atom_name))=='Cl-' .or. trim(adjustl(atom_name))=='SOD' .or. trim(adjustl(atom_name))=='CLA') then
			nions = nions + 1
			protein_stat = .false.
		elseif (trim(adjustl(record))=='TER' .and. protein_stat) then
			write(21,'(A3)') 'TER'	! write(*,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',iostat=ios) record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,(crd(j),j=1,3),occupancy,bfactor
		elseif (trim(adjustl(record))=='END') then
			write(21,'(A3)') 'END'
			exit
		elseif (protein_stat) then
			if (chain_id == " ") chain_id = 'X'
			write(21,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)') record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,(crd(j),j=1,3),occupancy,bfactor
			nres = res_num
			natoms = natoms + 1
		end if
		cycle
	end if
end do

close(11); close(21)
write(*,*); write(*,*) "natoms: ",natoms; write(*,*) "nions: ",nions; write(*,*) "nres: ",nres

! NOW CREATE THE TEMPORARY FILENAMES FOR DSSP AND PDB FILES.

do while (pdb_file_stat)
	call random_number(dummy)
	write(dummy_string,'(1I6)') 100000+nint(100000.*dummy)
	temp_pdb = 'temp_' // trim(adjustl(dummy_string)) // '.pdb'
	inquire(file=temp_pdb,exist=pdb_file_stat)
	!file_stat = .true.
	!print *, ran_file
	!print *, file_stat
end do
do while (dssp_file_stat)
	call random_number(dummy)
	write(dummy_string,'(1I6)') 100000+nint(100000.*dummy)
	dssp_out = 'temp_' // trim(adjustl(dummy_string)) // '.dssp'
	inquire(file=dssp_out,exist=dssp_file_stat)
	!file_stat = .true.
	!print *, ran_file
	!print *, file_stat
end do

! NOW CALCULATE THE ASA FOR THE ORIGINAL PDB FILE (WITH THE CRYSTAL STRUCTURE COORDINATES).

command = 'cp ' // trim(adjustl(pdb_out)) // ' ' // temp_pdb
write(*,*); write(*,*) 'CREATE A TEMPORARY PDB FILE FOR DSSP: '; write(*,*) trim(adjustl(command))
call system(trim(adjustl(command)))
command = 'dssp -i ' // temp_pdb // ' -o ' // dssp_out
write(*,*); write(*,*) 'CALCULATE ASA USING DSSP: '; write(*,*) trim(adjustl(command))
call system(trim(adjustl(command)))
call asa_outputter(dssp_out,asa_file_unit,nres,frame)
!call sleep(300)

! NOW READ AMBER MD TRAJECTORIES.

allocate(mdcrd(natoms+nions,3))
read(12,*)

do
	frame = frame + 1
	read(12,'(10F8.3)',iostat=ios) ((mdcrd(i,j),j=1,3), i = 1,natoms+nions)
	if (ios < 0) then
		exit
	elseif (ios > 0) then
		write(*,*) 'Sum Tin Wong! : Amber MDCRD file broken.'; stop
	else	! if (frame==1) then	!mod(frame,1) == 0) then
		i=0
		open(unit=21,file=trim(adjustl(pdb_out)),status='old')
		open(unit=23,file=temp_pdb,status='replace')		! Temporary pdb file containg the latest MD trajectory for DSSP ASA calculations.
		do	! WRITE OUT THE NEW PDB FILE IN UNIT=23
			read(21,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',iostat=ios) record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,(crd(j),j=1,3),occupancy,bfactor
			if (ios<0) then
				exit
			elseif (ios>0) then
				write(*,*) 'Sum Tin Wong! : PDB file broken.'; stop
			else
				if (trim(adjustl(record))=='TER') then
					write(23,'(A3)') 'TER'	! write(*,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)',iostat=ios) record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,(crd(j),j=1,3),occupancy,bfactor
				elseif (trim(adjustl(record))=='END') then
					write(23,'(A3)') 'END'
					exit
				else
					i = i + 1	! i counts the atoms, written so far, it is used as the index of mdcrd array.
					if (trim(adjustl(record))/='ATOM') then; write(*,*) 'Something is horribly wrong!: record = ',record; stop; endif
					write(23,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)') record,atom_num,atom_name,alt_loc_ind,res_name,chain_id,res_num,res_code,(mdcrd(i,j),j=1,3),occupancy,bfactor
				end if
			end if
			cycle
		end do
		close(21); close(23)
		call system(trim(adjustl(command)))
		call asa_outputter(dssp_out,asa_file_unit,nres,frame)
	end if
	read(12,'(10F8.3)') (box(j),j=1,3)
	if (mod(frame,10)==0) then
		write(*,'(1A6,1I6)') 'frame:',frame
		write(*,'(1A9,10F8.3)') 'box size:',(box(j),j=1,3)
	end if
	cycle
end do

! NOW REMOVE TEMPORARY FILES.

!command = 'rm ' // 'temp*'
!write(*,*); write(*,*) 'REMOVE ALL TEMPORARY FILES: '; write(*,*) trim(adjustl(command))
!call system(trim(adjustl(command)))

END PROGRAM dynamic_asa

include 'sub_asa_outputter.f90'
!include 'MRSA_finder.f90'
!include 'random_filename.f90'

