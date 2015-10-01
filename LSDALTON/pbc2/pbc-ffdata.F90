!> ???
MODULE pbcffdata
use files
	USE precision
	USE typedef
	IMPLICIT NONE
	PUBLIC

	!> ???
	INTEGER, PARAMETER :: min_qlmax = 3, max_qlmax = 21
	!> ???
	INTEGER, PARAMETER :: min_TWlmax = 5, max_TWlmax = 50
	
	!> ???
	TYPE ffdata_t
		!> ???
		LOGICAL :: reset_T_diverg 
		!> ???
		LOGICAL :: square_intermediates 
		!> Lmax for multipole moments
		INTEGER :: q_lmax     
		!> Lmax for T and W tensors
		INTEGER :: TW_lmax    
		!> ???
		REAL(realk), POINTER :: cell_locmom(:) 
		!> ???
		REAL(realk), POINTER :: Tlattice(:,:) 
	END TYPE ffdata_t

	TYPE(ffdata_t), SAVE :: ffdata

CONTAINS 
	!> \author JR
	!> \date 2013
	!> \brief ???
	!> \param tlattice 			???
	!> \param lmax					???
	!> \param mattxt 				???
	!> \param lupri 				Print unit for output datafile.
	SUBROUTINE read_pbc_tlattice(tlattice,lmax,mattxt,lupri)
		IMPLICIT NONE
		INTEGER,INTENT(IN) :: lmax,lupri
		REAL(realk),INTENT(INOUT) ::tlattice((lmax+1)**2,(lmax+1)**2)
		CHARACTER(LEN=18),INTENT(IN) :: mattxt
		! local
		INTEGER :: i,j,dimin,dimtest
		INTEGER :: iunit

		dimin=(lmax+1)**2
		iunit=-1

		write(*,*) mattxt
		call lsOPEN(IUNIT,mattxt,'OLD','UNFORMATTED')
		read(iunit) dimtest!,nil
		if(dimtest .ne. dimin) then
			write(*,*) 'ERROR: dims do not match: input dim=', &
				& dimin,'read dim=',dimtest
			call lsquit('In read_pbc_tlattice: dims do not match',lupri)
		endif

		do j=1,dimin
			read(iunit) (tlattice(j,i),i=1,dimin)
		enddo

		CALL lsCLOSE(IUNIT,'KEEP')

	END SUBROUTINE read_pbc_tlattice
END MODULE pbcffdata
