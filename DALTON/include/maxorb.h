!
!     File: maxorb.h
!
!     MXSHEL = maximum number of shells (insert shell definition here).
!              (if modified: also  change MXSHEL for __CVERSION__ in infpar.h)
!     MXPRIM = maximum number of primitives.
!     MXCORB = maximum number of orbitals (possibly contracted).
!     MAXOCC = maximum number of occupied orbitals
!
!     IF you change any of these parameters you should rebuild with "make".
!
      INTEGER    MXSHEL, MXPRIM, MXCORB, MXORBT, MAXOCC
      PARAMETER (MXSHEL = 1500, MXPRIM = 15000 )
      PARAMETER (MXCORB = 5000, MXORBT = MXCORB*(MXCORB + 1)/2 )
      PARAMETER (MAXOCC = 1500 )

!     MXCORB_CC = max number of orbitals in CC modules
!     (normally less than MXCORB because of a lot of static allocations
!      in the CC module for address pointers)

      INTEGER    MXCORB_CC
      PARAMETER (MXCORB_CC = 600 )

! -- end of maxorb.h --
