!	Amir Shahmoradi, Wednesday 6:48 PM, January 22, 2014, WilkeLab, ICMB, UT Austin.
!	ATTN: I wrote this code more than 6 months ago and I decided to update this code today, however, it seems too complex and time-consuming at the moment.
!		  So, I have decided to leave it in the working format that it is at the moment. Perhaps some time in future I will update it, if needed.
!	This program reads in a given pdb file, the corresponding Amber mdcrd file, and a pdb specification file (dical.in).
!	Then calculates the phi & psi dihedral angles for each of the sites in the protein.
!	Then outputs the average and standard deviation of the two angles for each site in two separate files.
!	It also outputs a newly written pdb file that contains a dummy chain number.
!	This file will be later used by DSSP software for the calculation of RSA.
!	Amir Shahmoradi, Thursday 10:39 AM, August 1, 2013, WilkeLab, UT Austin.
!	Amir Shahmoradi, Wednesday 3:27 AM, August 7, 2013, WilkeLab, ICMB, UT Austin.
!	I have made some modifications to the code to reduce the amount of hard-coding in the code.
!   Input data for ntot values in the previous code were all wrong, leading to wrong values for the calculated dihedral angles.
!	Amir Shahmoradi, Wednesday 3:27 AM, August 7, 2013, WilkeLab, ICMB, UT Austin.
PROGRAM dical
IMPLICIT NONE
INTEGER :: i=1,j,ii,frame,Ncounter=0,Ccounter=0,CAcounter=0,CBcounter=0
INTEGER :: chi1_4th_atom=0,ALA_GLY_count=0,non_ALA_GLY_count=0,convention=360
INTEGER :: natoms,ntot,nres,nframes
! natom: total number of protein atoms in the Amber pdb file
! ntot: total number of protein atoms in the Amber pdb file + any extra ions or water atoms or ... .
DOUBLE PRECISION, PARAMETER :: pi = 3.14159265359d0
CHARACTER(LEN=256) :: pdb_in,mdcrd_filename,dihedral_out,box_size_file
!CHARACTER(LEN=256), PARAMETER :: pdb_in = 'C:\Users\Amir\Documents\GitHub\structural_prediction_of_ER\molecular_dynamics\Amber\postproc\pdb\1RD8_X.pdb'
!CHARACTER(LEN=256), PARAMETER :: pdb_out = '1RD8_X.pdb'
!CHARACTER(LEN=256), PARAMETER :: mdcrd_filename = '1RD8_image.mdcrd'
!CHARACTER(LEN=256), PARAMETER :: dihedral_out = '1RD8_dihedral.txt'
!CHARACTER(LEN=100) :: header_format,row_format
!CHARACTER(LEN=1), PARAMETER :: chain_id = 'X', Alt_loc_ind
CHARACTER(LEN=6), DIMENSION(:), allocatable :: record
CHARACTER(LEN=4), DIMENSION(:), allocatable :: atom_name
CHARACTER(LEN=3), DIMENSION(:), allocatable :: res_name
CHARACTER(LEN=1), DIMENSION(:), allocatable :: chain_id,alt_loc_ind,res_code
INTEGER, DIMENSION(:), allocatable :: atom_num,res_num
INTEGER, DIMENSION(:), allocatable :: Npos,Cpos,CApos,CBpos,chi1_4th_atom_pos
!DOUBLE PRECISION, DIMENSION(nres) :: psi,phi
!DOUBLE PRECISION :: mean_sin_psi,mean_cos_psi,mean_psi,std_psi
!DOUBLE PRECISION :: mean_sin_phi,mean_cos_phi,mean_phi,std_phi
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: crd
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: mdcrd
DOUBLE PRECISION, DIMENSION(:),   allocatable :: occupancy,bfactor
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: phi	! , cos_psi,sin_psi,cos_phi,sin_phi,
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: psi	! , cos_psi,sin_psi,cos_phi,sin_phi,
DOUBLE PRECISION, DIMENSION(:,:), allocatable :: chi1
DOUBLE PRECISION, DIMENSION(3) :: BOX	! simulation box size.
DOUBLE PRECISION :: dihedral	! returns dihedral angle in degrees, according to the requested convention: [-180,180] or [0,360].
DOUBLE PRECISION :: ave_phi,var_phi,ave_psi,var_psi,ave_chi1,var_chi1

INTERFACE
	SUBROUTINE circstat(angles,ave,var)
		DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: angles
		DOUBLE PRECISION, INTENT(OUT) :: ave,var
	END SUBROUTINE circstat
END INTERFACE

namelist /pdb_1RD8/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_2JLY/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_2Z83/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_3GOL/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_3GSZ/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_3I5K/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_3LYF/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_4IRY/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,
namelist /pdb_4AQF_B/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file	!natoms,ntot,nres,nframes,
namelist /pdb_4GHA_A/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_1AJ8_A/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_1AOR_A/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_1CTS_A/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_1MP9_A/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_2JLY_temp_50/   	natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_2JLY_temp_100/  	natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_2JLY_temp_200/  	natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_2JLY_temp_450/  	natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_1AOR_A_temp_373/	natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file  	!natoms,ntot,nres,nframes,
namelist /pdb_2FP7_B/ 			natoms,ntot,nres,nframes,pdb_in,mdcrd_filename,dihedral_out,box_size_file   !natoms,ntot,nres,nframes,

!write(*,*) pdb_in
!call get_command_argument(1,pdb_in)
!call get_command_argument(2,mdcrd_filename)
!call get_command_argument(3,pdb_out)
!call get_command_argument(4,dihedral_out)
open(unit=10,file='dical.in',status='old')
!write(header_format,'(1A1,1I3,1A4)') '(',nres+1,'A10)'
!write(row_format,'(1A6,1I3,1A7)') '(1I10,',nres,',F10.3)'
!write(*,*) header_format
!write(*,*) row_format

!read(10,NML=pdb_1RD8)
!read(10,NML=pdb_2JLY)
!read(10,NML=pdb_2Z83)
!read(10,NML=pdb_3GOL)
!read(10,NML=pdb_3GSZ)
!read(10,NML=pdb_3I5K)
!read(10,NML=pdb_3LYF)
!read(10,NML=pdb_4IRY)
!read(10,NML=pdb_4AQF_B)
!read(10,NML=pdb_4GHA_A)
!read(10,NML=pdb_1AJ8_A)
!read(10,NML=pdb_1AOR_A)
!read(10,NML=pdb_1CTS_A)
!read(10,NML=pdb_1MP9_A)
!read(10,NML=pdb_2JLY_temp_50)
!read(10,NML=pdb_2JLY_temp_100)
!read(10,NML=pdb_2JLY_temp_200)
!read(10,NML=pdb_2JLY_temp_450)
!read(10,NML=pdb_1AOR_A_temp_373)
read(10,NML=pdb_2FP7_B)

allocate(record(natoms),atom_name(natoms),res_name(natoms),chain_id(natoms),alt_loc_ind(natoms), &
         res_code(natoms),atom_num(natoms),res_num(natoms),occupancy(natoms),bfactor(natoms))
allocate(Npos(nres),Cpos(nres),CApos(nres),CBpos(nres),chi1_4th_atom_pos(nres))
allocate(crd(natoms,3),mdcrd(ntot,3),phi(2:nres,nframes),psi(1:nres-1,nframes),chi1(1:nres,nframes))
CBpos=-9999; chi1_4th_atom_pos=-9999; chi1=-9999.

open(unit=11,file=trim(adjustl(pdb_in)),status='old')
open(unit=12,file=trim(adjustl(mdcrd_filename)),status='old')
open(unit=22,file=trim(adjustl(dihedral_out)),status='replace')
open(unit=24,file=box_size_file,status='replace')

! READ IN PROTEIN PDB FILE:
!read(11,*)

do while (i<=natoms)
	read(11,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)') &
	record(i),atom_num(i),atom_name(i),alt_loc_ind(i),res_name(i),chain_id(i), &
	res_num(i),res_code(i),(crd(i,j),j=1,3),occupancy(i),bfactor(i)
	chain_id(i) = 'X'
	!alt_loc_ind(i) = ' '
	if (trim(adjustl(record(i)))=='TER') then
		write(*,'(A3)') 'TER'
		!write(21,'(A3)') 'TER'
		!write(*,*) 'yes'
	else
		!write(21,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2)') &
		!record(i),atom_num(i),atom_name(i),alt_loc_ind(i),res_name(i),chain_id(i), &
		!res_num(i),res_code(i),(crd(i,j),j=1,3),occupancy(i),bfactor(i)
		if (trim(adjustl(atom_name(i)))=='N') then
			Ncounter = Ncounter + 1
			Npos(Ncounter) = atom_num(i)
		elseif (trim(adjustl(atom_name(i)))=='C') then
			Ccounter = Ccounter + 1
			Cpos(Ccounter) = atom_num(i)
		elseif (trim(adjustl(atom_name(i)))=='CA') then
			CAcounter = CAcounter + 1
			CApos(CAcounter) = atom_num(i)
		! for chi1 dihedral angle, find the fourth atom:
		elseif (trim(adjustl(atom_name(i)))=='CB' .and. res_name(i)/='ALA') then
			CBcounter = CBcounter + 1
			CBpos(CBcounter) = atom_num(i)
		elseif (trim(adjustl(atom_name(i)))=='CG' .or. trim(adjustl(atom_name(i)))=='OG' .or. &
		trim(adjustl(atom_name(i)))=='CG1' .or. trim(adjustl(atom_name(i)))=='OG1' .or. &
		trim(adjustl(atom_name(i)))=='SG') then
			non_ALA_GLY_count = non_ALA_GLY_count + 1
			chi1_4th_atom = chi1_4th_atom + 1
			chi1_4th_atom_pos(chi1_4th_atom) = atom_num(i)
		elseif (trim(adjustl(atom_name(i)))=='HA2' .or. trim(adjustl(atom_name(i)))=='HB1') then
			if (res_name(i) /= 'ALA' .and. res_name(i) /= 'GLY') then
				write(*,*) 'Residue recognized as ALA/GLY according to atom name: ', atom_name(i)
				write(*,*) 'However the residue name does not correspond to ALA/GLY: ', res_name(i)
				write(*,*) 'Press enter to end the program'
				read(*,*)
				STOP
			end if
			CBcounter = CBcounter + 1
			chi1_4th_atom = chi1_4th_atom + 1
			ALA_GLY_count = ALA_GLY_count + 1
		end if
		if (res_name(i)=='SEC') then
			write(*,*) 'Selenocysteine (SEC) was found among the protein amino acids.'
			write(*,*) 'Press enter to end the program'
			read(*,*)
			STOP
		end if
		i = i + 1
	end if
end do

!write(21,'(A3)') 'TER'
!write(21,'(A3)') 'END'
write(*,*) 'Ncounter:  ', Ncounter
write(*,*) 'Ccounter:  ', Ccounter
write(*,*) 'CAcounter: ', CAcounter

if (non_ALA_GLY_count+ALA_GLY_count/=nres .or. CBcounter/=chi1_4th_atom) then
	write(*,*) 'non_ALA_GLY_count+ALA_GLY_count/=nres OR CBcounter/=non_ALA_GLY_count'
	write(*,*) 'non_ALA_GLY_count: ', non_ALA_GLY_count
	write(*,*) 'ALA_GLY_count: ', ALA_GLY_count
	write(*,*) 'CBcounter: ', CBcounter
	write(*,*) 'Press enter to end the program'
	read(*,*)
	STOP
end if



read(12,*)
write(22,'(8A20)') 'res_Num','res_name','mean_phi','var_phi','mean_psi','var_psi', &
				   'mean_chi1','var_chi1'
!write(23,'(4A20)') 'res_Num','res_name','mean_phi','var_phi'
! START CALCULATING DIHEDRAL ANGLES:
do frame = 1,nframes
	read(12,'(10F8.3)') ((mdcrd(i,j),j=1,3), i = 1,ntot)
	do i = 1,nres-1
		phi(i+1,frame) = dihedral(mdcrd(Cpos(i),1:3),mdcrd(Npos(i+1),1:3),&
								  mdcrd(CApos(i+1),1:3),mdcrd(Cpos(i+1),1:3))
		psi(i,frame) = dihedral(mdcrd(Npos(i),1:3),mdcrd(CApos(i),1:3),&
								mdcrd(Cpos(i),1:3),mdcrd(Npos(i+1),1:3))	! *pi/180.
		if (chi1_4th_atom_pos(i) > 0) then
			if (CBpos(i)<0) then
				write(*,*) 'CBpos(i)<0, while chi1_4th_atom_pos(i) > 0 for residue number: ', i
				write(*,*) 'CBpos(i): ', CBpos(i)
				write(*,*) 'chi1_4th_atom_pos(i): ', chi1_4th_atom_pos(i)
				write(*,*) 'Press enter to end the program'
				read(*,*)
				stop
			end if
			chi1(i,frame) = dihedral(mdcrd(Npos(i),1:3),mdcrd(CApos(i),1:3),&
									 mdcrd(CBpos(i),1:3),mdcrd(chi1_4th_atom_pos(i),1:3))	! *pi/180.
		end if
		!cos_psi(i,frame) = dcos(psi(i,frame))
		!sin_psi(i,frame) = dsin(psi(i,frame))
		!cos_phi(i,frame) = dcos(phi(i,frame))
		!sin_phi(i,frame) = dsin(phi(i,frame))
	end do
	read(12,'(10F8.3)')  BOX(1), BOX(2), BOX(3)	! Simulation box size.
	!write(*,'(10F8.3)')  BOX(1), BOX(2), BOX(3)	! Simulation box size.
	write(24,'(10F8.3)')  BOX(1), BOX(2), BOX(3)	! Simulation box size.
	if (mod(frame,100)==0) write(*,*) 'Reading mdcrd frame number ',frame
end do

if (any(phi<0.) .or. any(psi<0.)) then
	write(*,*) 'Something is fishy here. psi and phi angles should not be negative according to the assumed range convention [0,360].'
	write(*,*) 'Press enter to continue ...'
	read(*,*)
end if

do ii=1,nres
	do j=1,nframes
		if (chi1(ii,j)<0. .and. chi1(ii,j) >-9998.) then
			write(*,*) 'Chi1 is calculated incorrectly for residue & frame numbers: ', ii,j
			write(*,*) 'Chi1: ', chi1(ii,j)
		end if
	end do
end do

!write(23,'(1I20,3A20)') res_num(Npos(1)),res_name(Npos(1)),'NA','NA'	! phi angle does not exist for the first residue
! First calculate the psi angle for the first residue. phi angle does not exist for the first residue.
call circstat(psi(1,1:nframes),ave_psi,var_psi)

if (any(chi1(1,1:nframes)<0.)) then
	write(22,'(1I20,3A20,2F20.4,2A20)') res_num(Npos(1)),res_name(Npos(1)),&
	'NA','NA',ave_psi,var_psi,'NA','NA'		! phi angle does not exist for the first residue
else
	call circstat(chi1(1,1:nframes),ave_chi1,var_chi1)
	write(22,'(1I20,3A20,4F20.4)') res_num(Npos(1)),res_name(Npos(1)),&
	'NA','NA',ave_psi,var_psi,ave_chi1,var_chi1		! phi angle does not exist for the first residue
end if

do ii = 2,nres-1
	call circstat(phi(ii,1:nframes),ave_phi,var_phi)
	call circstat(psi(ii,1:nframes),ave_psi,var_psi)
	if (any(chi1(ii,1:nframes)<0.)) then
		write(22,'(1I20,1A20,4F20.4,2A20)') res_num(Npos(ii)),res_name(Npos(ii)),&
		ave_phi,var_phi,ave_psi,var_psi,'NA','NA'
	else
		call circstat(chi1(ii,1:nframes),ave_chi1,var_chi1)
		write(22,'(1I20,1A20,6F20.4)') res_num(Npos(ii)),res_name(Npos(ii)),&
		ave_phi,var_phi,ave_psi,var_psi,ave_chi1,var_chi1
	end if
	!write(22,'(1I20,1A20,4F20.4)') res_num(Npos(i)),res_name(Npos(i)),ave_phi,var_phi,ave_psi,var_psi
	!mean_sin_psi = sum(sin_psi(i,1:nframes))/dble(nframes)
	!mean_cos_psi = sum(cos_psi(i,1:nframes))/dble(nframes)
	!mean_psi = datan2(mean_sin_psi,mean_cos_psi)*180./pi
	!write(23,'(1I20,1A20,2F20.4)') res_num(Npos(i+1)),res_name(Npos(i+1)),ave,var
	!mean_sin_phi = sum(sin_phi(i,1:nframes))/dble(nframes)
	!mean_cos_phi = sum(cos_phi(i,1:nframes))/dble(nframes)
	!mean_phi = datan2(mean_sin_phi,mean_cos_phi)
end do

call circstat(phi(nres,1:nframes),ave_phi,var_phi)		! psi angle does not exist for the last residue.

if (any(chi1(nres,1:nframes)<0.)) then
	write(22,'(1I20,1A20,2F20.4,4A20)') res_num(Npos(ii)),res_name(Npos(ii)),&
	ave_phi,var_phi,'NA','NA','NA','NA'
else
	call circstat(chi1(nres,1:nframes),ave_chi1,var_chi1)
	write(22,'(1I20,1A20,2F20.4,2A20,2F20.4)') res_num(Npos(ii)),res_name(Npos(ii)),&
	ave_phi,var_phi,'NA','NA',ave_chi1,var_chi1
end if
!write(22,'(1I20,1A20,2F20.4,2A20)') res_num(Npos(nres)),res_name(Npos(nres)),ave_phi,var_phi,'NA','NA'	! phi angle does not exist for the first residue

END PROGRAM dical


!   *************************************************************************************************************
!   *************************************************************************************************************
!   *************************************************************************************************************
!   *************************************************************************************************************
!   *************************************************************************************************************
	!include 'dihedral.f90'
	!include 'circstat.f90'
!   *************************************************************************************************************
!   *************************************************************************************************************
!   *************************************************************************************************************
!   *************************************************************************************************************
!   *************************************************************************************************************
!	This function takes 3D coordinates of four points as the input and calculates the dihedral in radians.
!	V.1 Amir Shahmoradi, Wednesday 8:33 PM, July 24, 2013, WilkeLab, UT Austin.
!	V.2 Amir Shahmoradi, Saturday 3:21 AM, August 3, 2013, WilkeLab, UT Austin.
	FUNCTION dihedral(crd1,crd2,crd3,crd4,convention)
	IMPLICIT NONE
	INTEGER, INTENT(IN), OPTIONAL :: convention
	DOUBLE PRECISION, PARAMETER :: pi = 3.14159265359d0
	DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: crd1,crd2,crd3,crd4
	DOUBLE PRECISION, DIMENSION(3) :: b1,b2,b3,cross_b1b2,cross_b2b3
	DOUBLE PRECISION :: dihedral
	DOUBLE PRECISION :: len_b2	! len_b1,len_b2,len_b3
	INTEGER, SAVE :: counter=0
	counter = counter + 1
	b1 = crd2 - crd1
	b2 = crd3 - crd2
	b3 = crd4 - crd3
	!len_b1 = dsqrt((crd1(1)-crd2(1))**2 + (crd1(2)-crd2(2))**2 + (crd1(3)-crd2(3))**2)
	len_b2 = dsqrt((crd2(1)-crd3(1))**2 + (crd2(2)-crd3(2))**2 + (crd2(3)-crd3(3))**2)
	!len_b3 = dsqrt((crd3(1)-crd4(1))**2 + (crd3(2)-crd4(2))**2 + (crd3(3)-crd4(3))**2)
	call cross(b1,b2,cross_b1b2)
	call cross(b2,b3,cross_b2b3)
	!print *, present(convention)
	!read(*,*)
	if (present(convention)) then
		if (convention==180) then
			if (counter==1) then
				write(*,*) 'Requested dihedral angle convention is 180.'
				write(*,*) 'Setting the angle range to [-180,180] ...'
				write(*,*) 
			end if
			dihedral = datan2(dot_product(len_b2*b1,cross_b2b3),dot_product(cross_b1b2,cross_b2b3))*180./pi	! in degrees
		elseif (convention==360) then
			if (counter==1) then
				write(*,*) 'Requested dihedral angle convention is 360.'
				write(*,*) 'Setting the angle range to [0,360] ...'
				write(*,*) 
			end if
			dihedral = datan2(dot_product(len_b2*b1,cross_b2b3),dot_product(cross_b1b2,cross_b2b3))*180./pi	! in degrees
			if (dihedral < 0.) dihedral = 360 + dihedral
		else
			if (counter==1) then
				write(*,*) 'Inappropriate input dihedral convention: ', convention
				write(*,*) 'Input dihedral convention must be either 180 or 360.'
				write(*,*) 'Dihedral angle range will be reset to default: 360 --> [0,360]'
				write(*,*) 
			end if
			dihedral = datan2(dot_product(len_b2*b1,cross_b2b3),dot_product(cross_b1b2,cross_b2b3))*180./pi	! in degrees
			if (dihedral < 0.) dihedral = 360 + dihedral
		end if
	else
		if (counter==1) then
			write(*,*) 'No dihedral angle convention available in the input.'
			write(*,*) 'Dihedral angle range will be set to default: [0,360].'
		end if
		dihedral = datan2(dot_product(len_b2*b1,cross_b2b3),dot_product(cross_b1b2,cross_b2b3))*180./pi	! in degrees
		if (dihedral < 0.) dihedral = 360 + dihedral
	end if
	!if (mod(counter,20000)==0) write(*,*) 'dihedral angle calculation in progress...', counter, 'dihedral measures done.'
	END FUNCTION dihedral
	
	SUBROUTINE cross(a,b,c)
	! CALCULATES THE CROSS PRODUCT OF TWO VECTORS.
	DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: a, b
	DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: c
	c(1) = a(2) * b(3) - a(3) * b(2)
	c(2) = a(3) * b(1) - a(1) * b(3)
	c(3) = a(1) * b(2) - a(2) * b(1)
	END SUBROUTINE cross



	
!	ATTN: this code is set to output the average of the angles in the range [-180,180].
!	It can be changed to [0,360] by uncommenting line 26 below.
!	Amir Shahmoradi, Saturday 1:42 AM, August 4, 2013, Wilke Lab., ICMB, Universty of Texas at Austin.
!	Given an array of angles with the range [0,360], this subroutine spits out the circular average of the angles (ave) in degrees, also the circular variance (var) which should be in the range [0,1].
!	NOTE: If the range of input angles is [-180,180], then the angles will be tranformed to the range [0,360] in the output.
!	Amir Shahmoradi, Saturday 6:12 PM, August 3, 2013, Wilke Lab., Universty of Texas at Austin.
	SUBROUTINE circstat(angles,ave,var)
	IMPLICIT NONE
	INTEGER :: i,n
	DOUBLE PRECISION, PARAMETER :: pi = 3.14159265359d0
	DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: angles
	DOUBLE PRECISION, DIMENSION(size(angles)) :: sines,cosines
	DOUBLE PRECISION :: mean_sines,mean_cosines
	DOUBLE PRECISION, INTENT(OUT) :: ave,var
	n = size(angles)
	! First make sure all angles are defined in the range [0,360].
	if (any(angles<0.)) then
		write(*,*) 'The range of the input array of angles is not [0,360].'
		write(*,*) 'The input array of angles will be tranformed to the range [0,360].'
	end if
	forall (i=1:n, angles(i) < 0.) angles(i) = angles(i) + 360.
	sines = dsin(angles*pi/180.)
	cosines = dcos(angles*pi/180.)
	mean_sines = sum(sines(1:n))/dble(n)
	mean_cosines = sum(cosines(1:n))/dble(n)
	ave = datan2(mean_sines,mean_cosines)*180./pi
	!if (ave<0.d0) ave = ave + 360
	var = 1.d0 - dsqrt(mean_sines**2+mean_cosines**2)
	if (var > 1.d0 .or. var < 0.d0) then
		write(*,*) 'Calculated circular variance: ', var
		write(*,*) 'Something is fishy here. Circular variance cannot be out of [0,1] range.'
		write(*,*) 'Press enter to continue ...'
		read(*,*)
	end if
	END SUBROUTINE circstat