!     File : ccom.h -- permanent info about basis functions
!         (l = J + 1, i.e. e.g. p-orbital is J = 2)
!         NHTYP - maximum angular quantum number + 1 for ALL orbitals
!         KHK(J) - number of spherical (cartesian) components for given J
!         KCK(J) - number of Cartesian components for given J
!         NHKOFF(J) - offset for components in list of l-functions
!         GTOTYP(index(K,J)) - label for K component of J-type orbital
!         DOCART - Cartesian basis fu., if false:  spherical or own def.
!         SPH(J) - true if cartesian to spherical/own basis needed
!         SPHNRM - true if all basis functions are normalized
!                  (only true for spherical, not for cartesian or own)
!
!     NHTYP  is set in BASINP(herrdn.F),
!     KHK, KCK, SPH are set in SPHINP(herrdn.F), based on user input
!     NHKOFF is set in BASPAR(herrdn.F)
!     DOCART are set in line 4 from .mol READI1(herrdn.F)
!     SPHNRM true if spherical basis funcitons (herrdn.F)
!     GTOTYP is set in CARLAB, SPHLAB or input (herrdn.F)


      REAL*8  THRS
      INTEGER NHTYP,  KHK(MXQN), KCK(MXQN), NHKOFF(MXQN), CCOMlast
      LOGICAL DOCART, SPH(MXQN), SPHNRM
      COMMON /CCOM/ THRS,                                               &
     &              NHTYP, KHK, KCK, NHKOFF,                            &
     &              DOCART, SPH, SPHNRM,                                &
     &   CCOMlast !  Very important:
      !  Always keep CCOMlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

      CHARACTER*4 GTOTYP
      INTEGER     CCOMClast
      COMMON /CCOMC/ GTOTYP(MXQN*(MXQN+1)*(MXQN+2)/6),                  &
     &   CCOMClast !  Very important:
      !  Always keep CCOMClast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
!     --- end of ccom.h ---
