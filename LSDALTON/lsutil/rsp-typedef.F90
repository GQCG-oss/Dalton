MODULE RSP_type_module
 use matrix_module
!FIXME: does the use of definition 'real' only for freq and Output
!       have consequencies on the accuracy of the result????
Type Equation_item
 real(realk)      :: freq
 CHARACTER(len=8) :: Glabel
 integer          :: rsp_location
end Type Equation_item
type QR_ITEM
 integer             :: number,g_index,f_index
 real(realk)         :: Output,LROutput(2)
 type(Equation_item) :: Equations(3)
end type QR_ITEM
!Derived type(s) used in PDBS subroutines: see below for explanation
Type LOOP_item
 CHARACTER(len=8)    :: Alabel,Blabel
end type LOOP_item
type PDBS_ITEM
 CHARACTER(len=6)    :: PropLabel      
 integer             :: number,Xindex,PDBSvecs
 real(realk)         :: freq
 !type(LOOP_item)     :: LOOP(9)
 type(LOOP_item)     :: LOOP(6)    !6 should suffice
 integer             :: iatom
end type PDBS_ITEM
type PDBSvec
 CHARACTER(len=8)    :: Label
 real(realk)         :: freq
 Logical             :: A1
end type PDBSvec
! Perturbation dependent basis set calculation:
! In PDBS-PROPERTY-Interface.f90 we will use a PDBS_STACK(PDBSitems), which is
! a vector of PDBS_item's. Each PDBS_item collects information to run either:
!  An MCD calculation for a given excitation 
!  An Herzberg-Teller calculation for a given excitation 
!  An Excited state gradient for a given excitation 
!  An Raman calculation for a given frequency 
!  A  Verdet calculation for a given frequency 
! The dimension PDBSitems is determined in config.f90 from data in LSDALTON.INP
!
!  So if you want to do a Raman calculation for 2 freq. you have 
!  2 PDBS_item's in the PDBS_STACK
!
!  The PDBS_ITEM is a derived type which contains:
!
!  Proplabel:
!    being either 'RAMAN ','VERDET','MCD   ','IR int' or 'ExGrad'
!    and is at the moment only used to define how to do the Build_AVERAGES 
!    The Proplabel can be accessed as 
!        PDBS_STACK(PDBS_number)%Proplabel='RAMAN '
!    where the 'PDBS_number' is the index of the vector you are interested in
!    (we basically loop PDBS_number=1,PDBSitems)
!
!  number:
!   Indicates how many Xbar one needs to calculate, hence 
!   how many times one should loop over. It relates to the number
!   of different pairs of <<A,B>> required (as its RHS depends on A and B)
!
!  Xindex:
!   Used for MCD, Excited state gradient, Herzberg-Teller calculations 
!   to indicate which excitationvector/exfreq one is working with
! 
!  freq:
!   The value of frequency for which to calculate the (rsp) A and/or B vector
!   - For Raman and Verdet the value assigned to freq is the 
!     laser frequency specified by the user in the input.
!   - NB: For Excited state gradient, MCD and Herzberg-Teller  the excitation 
!     frequency is stored on freq when loading the excitation vector 
!     from file.
!
!  PDBSvecs:
!   The number of linear response vectors needed:
!   - MCD/Herzberg-Teller: 1 response vector  (+ 1 excitation vector)
!   - Exc. state gradient: 0 response vectors (only exc vectors required)
!   - Raman/Verdet       : 2 response vectors (N^A,N^B)
!  The value of PDBSvecs is used to decide whether to add certain contributions

!  LOOP:
!    Is a vector of a derived type called LOOP_item
!    And regulates which things should be looped over (outer loop J)
!
!    E.g. VERDET/RAMAN this means that it contains the 2 pairs of labels for 
!    each of the 'number'=3/3*iatom loops
!
!  'XDIPLEN ','YDIPLEN '
!  'YDIPLEN ','ZDIPLEN '
!  'ZDIPLEN ','XDIPLEN '
!
!  iatom:
!    indicates whether or not it is a magnetic or 
!    geometrical derivative. It is set equal to
!    iatom = 0                 => Magnetic derivatives
!    iatom = Number of atoms   => geometrical derivatives
!-----------------------------------------------------------------------


! the type PDBSvec keeps track on which vectors has been written on disk 
! such that no vectors are calculated twice 

contains

!Added to avoid "has no symbols" linking warning
subroutine RSP_type_module_void()
end subroutine RSP_type_module_void

end Module RSP_type_module

