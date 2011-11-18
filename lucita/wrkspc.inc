!
!     LUCITA work space dummy array: WORK
!
!     note by s. knecht - mai 2011
!
!     the WORK array declared below is used only inside LUCITA and independent from 
!     the Dalton WORK array. 
!     historically, LUCITA was part of a standalone program LUCIA which had its own 
!     memory manager. 
!     in the startup of each LUCITA run we thus calculate
!     the offset from LUCITA-WORK to DALTON-WORK and take then piece-by-piece
!     from the allocated memory, alway referring via pointers 
!     to the pre-calculated offset.
!
!     as detected with gfortran-4.6 (x86-64, mac system), 
!     it appears that the original declaration of the LUCITA-WORK as WORK(1) 
!     may lead to the assumption of a scalar variable and the memory manager setup
!
!     address x: start marker (0.123456789d0)
!     address y: end marker   (0.987654321d0)
!
!     did not work anymore.
!     The solution was thus to declare WORK as an array of size of (at least) 2:
 
      COMMON/BIGGY/WORK(2)
