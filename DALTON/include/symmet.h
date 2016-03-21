! --- FILE: symmet.h ---
!
! Symmetry information for basis functions,
! mostly generated in abacus/herrdn.F and abacus/hergp.F.
!
! Below: 0 .le. i .and. i .le. MAXREP (irrep index); i1 = i + 1 = irrep no.
!        1 .le. j .and. j .le. KMAX
!
! FMULT(i) =
! PT(i)   =
! MAXREP  = 2**NSYMOP - 1; NSYM = MAXREP + 1 is number of irreps
! MAXOPR  = maximum number of operations necessary to loop over [0..7]
! MULT(i) =
! ISYMAX(q,1) = irrep of q-axis
! ISYMAX(q,2) = irrep of rotation around q-axis
! ISYMAO(,) =
! NPARSU(8)   : offset pointer for symmetry dependent AOs for given irrep
! NAOS(8)     : Number of AO'S for given irrep
! NPARNU(8,8) : offset pointer from non-symmetric operators for given irrep
! ...
! IPTCNT(:,:,1:2) - symmetry adapted nuclear coordinates/magnetic moments
! NCRREP(:,1:2) - number of symmetry adapted coordinates/magnetic moments
! IPTCOR - points from symmetry coordinate to generating Cartesian coordinat
! ...
! NCOS(8,-2)  : Number of AOs for density fitting class 2 for given irrep
! NCOS(8,-1)  : Number of AOs for density fitting class 1 for given irrep
! NCOS(8, 0)  : Number of AOs for Huckel for given irrep
! NCOS(8, 1)  : Number of AOs for Large(DIRAC) / Main(LMULBS) for given irrep
! NCOS(8, 2)  : Number of AOs for Small(DIRAC) / Auxiliary(LMULBS) for given irrep
! ...
! ICLASS(j) =
! ICNTAO(j) =
! PCM_IGEN(4) : contains the number of generators and the generators of the group 
!		to be passed to PCMSolver		 
!
#include "mxbsets.h"
!
!     length of SYMMTI (for MPI) :
!     call getbytespan(I1_SYMMTI, SYMMTILAST, SizeInBytes)
!     This returns the number of bytes to the variable SizeInBytes. It includes I1_symmti but
!     excludes symmtilast (which is a variable that is never used for anything
!     but to keep track of common block lengths).
!     Assuming I1_symmti is a dummy variable, this means you can also get the size in bytes if you call:
!     call getbytesize(maxrep, symmtilast, SizeInBytes)

      INTEGER I1_SYMMTI, MAXREP, MAXOPR, MULT(0:7), ISYMAX(3, 2),       &
     &        ISYMAO(MXQN, MXAQN), NPARSU(8), NAOS(8), NPARNU(8, 8),    &
     &        IPTSYM(MXCORB, 0:7), IPTCNT(3*MXCENT, 0:7, 2),            &
     &        NCRREP(0:7, 2), IPTCOR(3*MXCENT, 2), NAXREP(0:7, 2),      &
     &        IPTAX(3, 2), IPTXYZ(3, 0:7, 2), IPTNUC(MXCENT, 0:7),      &
     &        NROTS, NINVC, NREFL, IXVAL(0:7, 0:7),                     &
     &        ICLASS(MXCORB), ICNTAO(MXCORB), IPAR(0:7), JSOP(0:7),     &
     &        NCOS(8,-mxbsets:mxbsets), ICOS(8,-mxbsets:mxbsets),       &
     &        I2COSX(8,8,-mxbsets:mxbsets,-mxbsets:mxbsets),            &
     &        SYMMTIlast,  SYMMTRlast,                                  &
     &        PCM_IGEN(4), I2_SYMMTI 

      COMMON /SYMMTI/ I1_SYMMTI, MAXREP, MAXOPR, MULT, ISYMAX, ISYMAO,  &
     &                NPARSU, NAOS,  NPARNU, IPTSYM, IPTCNT, NCRREP,    &
     &                IPTCOR, NAXREP, IPTAX, IPTXYZ, IPTNUC,            &
     &                NROTS,  NINVC,  NREFL, IXVAL,  ICLASS, ICNTAO,    &
     &                IPAR,   JSOP,   NCOS,  ICOS,   I2COSX, PCM_IGEN,  &
     &                I2_SYMMTI,                                        & !     OBS! ALWAYS add new variables BEFORE the end tag: I2_SYMMTI
     &   SYMMTIlast !  Very important:
      !  Always keep SYMMTIlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.


      INTEGER LEN_SYMMTR
      PARAMETER (LEN_SYMMTR = 16) ! length of SYMMTR, used in MPI calls 
      ! Beyer says: Note that this distance is hardcoded. If for whatever 
      ! reason the length of the block changes you can still use 
      ! getbytesize(fmult, symmtrlast, SizeInBytes) to get the actual 
      ! size of the common block.

      REAL*8 FMULT(0:7), PT(0:7)
      COMMON /SYMMTR/ FMULT, PT,                                        &
     &   SYMMTRlast !  Very important:
      !  Always keep SYMMTRlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.

! --- end of symmet.h ---
