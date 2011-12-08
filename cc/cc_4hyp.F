!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
c /* deck CC_4HYP */
*=====================================================================*
       SUBROUTINE CC_4HYP(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of frequency-dependent fourth
*             hyperpolarizabilities for the Couple Cluster models
*
*                        CCS, CC2, CCSD
*
*     Written by Christof Haettig, maj 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccr1rsp.h>
#include <cc5rinf.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_4HYP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )

      CHARACTER*10 MODEL
      LOGICAL LADD
      INTEGER LWORK
      INTEGER KHYPPOL, NBHYPPOL
      INTEGER KEND0, LEND0, MXTRAN3, MXTRAN4, MXVEC2, MXVEC3
      INTEGER K0HTRAN, K0HDOTS, K0HCONS, N0HTRAN, MX0HTRAN, MX0HDOTS
      INTEGER K1HTRAN, K1HDOTS, K1HCONS, N1HTRAN, MX1HTRAN, MX1HDOTS
      INTEGER K0GTRAN, K0GDOTS, K0GCONS, N0GTRAN, MX0GTRAN, MX0GDOTS
      INTEGER K1GTRAN, K1GDOTS, K1GCONS, N1GTRAN, MX1GTRAN, MX1GDOTS
      INTEGER K1FATRAN,K1FADOTS,K1FACONS,N1FATRAN,MX1FATRAN,MX1FADOTS
      INTEGER KDTRAN,  KDDOTS,  KDCONS,  NDTRAN,  MXDTRAN,  MXDDOTS
      INTEGER KCTRAN,  KCDOTS,  KCCONS,  NCTRAN,  MXCTRAN,  MXCDOTS
      INTEGER KBTRAN,  KBDOTS,  KBCONS,  NBTRAN,  MXBTRAN,  MXBDOTS
      INTEGER KBATRAN, KBADOTS, KBACONS, NBATRAN, MXBATRAN, MXBADOTS
      INTEGER K0FTRAN, K0FDOTS, K0FCONS, N0FTRAN, MX0FTRAN, MX0FDOTS
      INTEGER KCHI3,   KCHIDOTS,KCHICONS,NXTRAN,  MXCHI3T,  MXCHI3D
      INTEGER IOPT

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION WORK(LWORK)
#endif

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<<  OUTPUT FROM COUPLED CLUST',
     &                                'ER PENTIC RESPONSE  >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<     CALCULATION OF CC HYPER',
     &                                'POLARIZABILITIES    >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD) ) THEN
         CALL QUIT('CC_4HYP called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPR5HYP .GT. 10) WRITE(LUPRI,*) 'CC_4HYP Workspace:',LWORK
  
*---------------------------------------------------------------------*
* allocate & initialize work space for hyperpolarizabilities
*---------------------------------------------------------------------*
      NBHYPPOL = N5ROPER * N5RFREQ

      KHYPPOL = 1
      KEND0   = KHYPPOL + 2*NBHYPPOL


      MXVEC2  = (NLRTLBL+1) * NLRTLBL / 2
      MXVEC3  = (NLRTLBL+2) * MXVEC2  / 3
      MXTRAN3 = NLRTLBL * NLRTLBL * NLRTLBL
      MXTRAN4 = NLRTLBL * NLRTLBL * NLRTLBL * NLRTLBL
      MXVEC2  = MIN(MXVEC2,2*180*NBHYPPOL)  ! 180 = # perm. for B^AB{O}
      MXVEC3  = MIN(MXVEC3,2*180*NBHYPPOL)  ! 180 = # perm. for B^AB{O}
      MXTRAN3 = MIN(MXTRAN3,2*180*NBHYPPOL) ! 180 = # perm. for B^AB{O}
      MXTRAN4 = MIN(MXTRAN4,2*180*NBHYPPOL) ! 180 = # perm. for B^AB{O}

      MX0HTRAN  = 5 * MXTRAN4
      MX1HTRAN  = 5 * MXTRAN4
      MX0GTRAN  = 4 * MXTRAN4
      MX1GTRAN  = 4 * MXTRAN4
      MX1FATRAN = MXDIM_FATRAN * MXTRAN4
      MXDTRAN   = 5 * MXTRAN4
      MXCTRAN   = 4 * MXTRAN4
      MXBTRAN   = 3 * MXTRAN4
      MXBATRAN  = MXDIM_BATRAN * MXTRAN4
      MX0FTRAN  = 3 * MXTRAN3
      MXCHI3T   = 1 * MXTRAN3

      MX0HDOTS  = MXTRAN4 * MXVEC2
      MX1HDOTS  = MXTRAN4 * MXVEC2
      MX0GDOTS  = MXTRAN4 * MXVEC2
      MX1GDOTS  = MXTRAN4 * MXVEC2
      MX1FADOTS = MXTRAN4 * MXVEC2
      MXDDOTS   = MXTRAN4 * MXVEC2
      MXCDOTS   = MXTRAN4 * MXVEC2
      MXBDOTS   = MXTRAN4 * MXVEC2
      MXBADOTS  = MXTRAN4 * MXVEC2
      MX0FDOTS  = MXTRAN3 * MXVEC3
      MXCHI3D   = MXTRAN3 * MXVEC3


      K0HTRAN  = KEND0
      K0HDOTS  = K0HTRAN  + MX0HTRAN
      K0HCONS  = K0HDOTS  + MX0HDOTS
      KEND0    = K0HCONS  + MX0HDOTS

      K1HTRAN  = KEND0
      K1HDOTS  = K1HTRAN  + MX1HTRAN
      K1HCONS  = K1HDOTS  + MX1HDOTS
      KEND0    = K1HCONS  + MX1HDOTS

      K0GTRAN  = KEND0
      K0GDOTS  = K0GTRAN  + MX0GTRAN
      K0GCONS  = K0GDOTS  + MX0GDOTS
      KEND0    = K0GCONS  + MX0GDOTS

      K1GTRAN  = KEND0
      K1GDOTS  = K1GTRAN  + MX1GTRAN
      K1GCONS  = K1GDOTS  + MX1GDOTS
      KEND0    = K1GCONS  + MX1GDOTS

      K1FATRAN = KEND0
      K1FADOTS = K1FATRAN + MX1FATRAN
      K1FACONS = K1FADOTS + MX1FADOTS
      KEND0    = K1FACONS + MX1FADOTS

      KDTRAN   = KEND0
      KDDOTS   = KDTRAN   + MXDTRAN
      KDCONS   = KDDOTS   + MXDDOTS
      KEND0    = KDCONS   + MXDDOTS

      KCTRAN   = KEND0
      KCDOTS   = KCTRAN   + MXCTRAN
      KCCONS   = KCDOTS   + MXCDOTS
      KEND0    = KCCONS   + MXCDOTS

      KBTRAN   = KEND0
      KBDOTS   = KBTRAN   + MXBTRAN
      KBCONS   = KBDOTS   + MXBDOTS
      KEND0    = KBCONS   + MXBDOTS

      KBATRAN  = KEND0
      KBADOTS  = KBATRAN  + MXBATRAN
      KBACONS  = KBADOTS  + MXBADOTS
      KEND0    = KBACONS  + MXBADOTS

      K0FTRAN  = KEND0
      K0FDOTS  = K0FTRAN  + MX0FTRAN
      K0FCONS  = K0FDOTS  + MX0FDOTS
      KEND0    = K0FCONS  + MX0FDOTS

      KCHI3    = KEND0
      KCHIDOTS = KCHI3    + MXCHI3T
      KCHICONS = KCHIDOTS + MXCHI3D
      KEND0    = KCHICONS + MXCHI3D

      LEND0    = LWORK - KEND0

      IF (LEND0.LT.0) THEN
        WRITE (LUPRI,*) 'KEND0,LEND0:',KEND0,LEND0
        CALL QUIT('Insufficient memory in CC_4HYP.')
      END IF

* initialize hyperpolarizabilities and the lists of contributions:
      CALL DZERO(WORK(1),KEND0-1)
     
      IF (LOCDBG) THEN
        CALL DCOPY(MX0HDOTS,8.8888d88,0,WORK(K0HCONS),1)
        CALL DCOPY(MX1HDOTS,8.8888d88,0,WORK(K1HCONS),1)
        CALL DCOPY(MX0GDOTS,8.8888d88,0,WORK(K0GCONS),1)
        CALL DCOPY(MX1GDOTS,8.8888d88,0,WORK(K1GCONS),1)
        CALL DCOPY(MX1FADOTS,8.8888d88,0,WORK(K1FACONS),1)
        CALL DCOPY(MXDDOTS,8.8888d88,0,WORK(KDCONS),1)
        CALL DCOPY(MXCDOTS,8.8888d88,0,WORK(KCCONS),1)
        CALL DCOPY(MXBDOTS,8.8888d88,0,WORK(KBCONS),1)
        CALL DCOPY(MXBADOTS,8.8888d88,0,WORK(KBACONS),1)
        CALL DCOPY(MX0FDOTS,8.8888d88,0,WORK(K0FCONS),1)
        CALL DCOPY(MXCHI3D,8.8888d88,0,WORK(KCHICONS),1)
      END IF

*---------------------------------------------------------------------*
* set up lists for H, G and F transformations and ETA vectors:
*---------------------------------------------------------------------*

      LADD = .FALSE.

      CALL CC5R_SETUP(MXTRAN3, MXVEC3, MXTRAN4, MXVEC2,
     &     WORK(K0HTRAN), WORK(K0HDOTS), WORK(K0HCONS), N0HTRAN,
     &     WORK(K1HTRAN), WORK(K1HDOTS), WORK(K1HCONS), N1HTRAN,
     &     WORK(K0GTRAN), WORK(K0GDOTS), WORK(K0GCONS), N0GTRAN,
     &     WORK(K1GTRAN), WORK(K1GDOTS), WORK(K1GCONS), N1GTRAN,
     &     WORK(K1FATRAN),WORK(K1FADOTS),WORK(K1FACONS),N1FATRAN,
     &     WORK(KDTRAN),  WORK(KDDOTS),  WORK(KDCONS),  NDTRAN,
     &     WORK(KCTRAN),  WORK(KCDOTS),  WORK(KCCONS),  NCTRAN,
     &     WORK(KBTRAN),  WORK(KBDOTS),  WORK(KBCONS),  NBTRAN,
     &     WORK(KBATRAN), WORK(KBADOTS), WORK(KBACONS), NBATRAN,
     &     WORK(K0FTRAN), WORK(K0FDOTS), WORK(K0FCONS), N0FTRAN,
     &     WORK(KCHI3),   WORK(KCHIDOTS),WORK(KCHICONS),NXTRAN,
     &     WORK(KHYPPOL), NBHYPPOL, LADD )

*---------------------------------------------------------------------*
* calculate H matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_HMAT('L0 ','R2 ','R1 ','R1 ','R2 ',N0HTRAN, MXVEC2,
     &              WORK(K0HTRAN),WORK(K0HDOTS),WORK(K0HCONS),
     &              WORK(KEND0), LEND0, IOPT )

      IOPT = 5
      CALL CC_HMAT('L1 ','R1 ','R1 ','R1 ','R2 ',N1HTRAN, MXVEC2,
     &              WORK(K1HTRAN),WORK(K1HDOTS),WORK(K1HCONS),
     &              WORK(KEND0), LEND0, IOPT )

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_GMATRIX('L0 ','R2 ','R2 ','R2 ',N0GTRAN, MXVEC2,
     &              WORK(K0GTRAN),WORK(K0GDOTS),WORK(K0GCONS),
     &              WORK(KEND0), LEND0, IOPT )

      IOPT = 5
      CALL CC_GMATRIX('L1 ','R2 ','R1 ','R2 ',N1GTRAN, MXVEC2,
     &              WORK(K1GTRAN),WORK(K1GDOTS),WORK(K1GCONS),
     &              WORK(KEND0), LEND0, IOPT )

*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      CALL CCQR_FADRV('L1 ','o1 ','R2 ','R2 ',N1FATRAN, MXVEC2,
     &                 WORK(K1FATRAN),WORK(K1FADOTS),
     &                 WORK(K1FACONS),
     &                 WORK(KEND0), LEND0, 'DOTP' )

*---------------------------------------------------------------------*
* calculate D matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_DMAT(WORK(KDTRAN),NDTRAN,
     &             'R1 ','R1 ','R1 ','R1 ',IOPT,'L2 ',
     &             WORK(KDDOTS),WORK(KDCONS),MXVEC2,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate C matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_CMAT(WORK(KCTRAN),NCTRAN,'R2 ','R1 ','R1 ',IOPT,'L2 ',
     &             WORK(KCDOTS),WORK(KCCONS),MXVEC2,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate B matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_BMAT(WORK(KBTRAN),NBTRAN,'R2 ','R2 ',IOPT,'L2 ',
     &             WORK(KBDOTS),WORK(KBCONS),MXVEC2,
     &             .FALSE.,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate B{O} matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_BAMAT(WORK(KBATRAN),NBATRAN,'o1 ','R2 ','R1 ',IOPT,'L2 ',
     &            WORK(KBADOTS),WORK(KBACONS),MXVEC2,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_FMATRIX(WORK(K0FTRAN),N0FTRAN,'L0 ','R3 ',IOPT,'R3 ',
     &                WORK(K0FDOTS),WORK(K0FCONS),MXVEC3,
     &                WORK(KEND0), LEND0)

*---------------------------------------------------------------------*
* calculate chi3 vector contributions:
*---------------------------------------------------------------------*
      CALL CC_DOTDRV('X3 ','R3 ',NXTRAN,MXVEC3,
     &               WORK(KCHI3), WORK(KCHIDOTS), WORK(KCHICONS),
     &               WORK(KEND0), LEND0 )

*---------------------------------------------------------------------*
* collect contributions and add them to hyperpolarizabilities
*---------------------------------------------------------------------*

      LADD = .TRUE.

      CALL CC5R_SETUP(MXTRAN3, MXVEC3, MXTRAN4, MXVEC2,
     &     WORK(K0HTRAN), WORK(K0HDOTS), WORK(K0HCONS), N0HTRAN,
     &     WORK(K1HTRAN), WORK(K1HDOTS), WORK(K1HCONS), N1HTRAN,
     &     WORK(K0GTRAN), WORK(K0GDOTS), WORK(K0GCONS), N0GTRAN,
     &     WORK(K1GTRAN), WORK(K1GDOTS), WORK(K1GCONS), N1GTRAN,
     &     WORK(K1FATRAN),WORK(K1FADOTS),WORK(K1FACONS),N1FATRAN,
     &     WORK(KDTRAN),  WORK(KDDOTS),  WORK(KDCONS),  NDTRAN,
     &     WORK(KCTRAN),  WORK(KCDOTS),  WORK(KCCONS),  NCTRAN,
     &     WORK(KBTRAN),  WORK(KBDOTS),  WORK(KBCONS),  NBTRAN,
     &     WORK(KBATRAN), WORK(KBADOTS), WORK(KBACONS), NBATRAN,
     &     WORK(K0FTRAN), WORK(K0FDOTS), WORK(K0FCONS), N0FTRAN,
     &     WORK(KCHI3),   WORK(KCHIDOTS),WORK(KCHICONS),NXTRAN,
     &     WORK(KHYPPOL), NBHYPPOL, LADD )

*---------------------------------------------------------------------*
* print output & return:
*---------------------------------------------------------------------*

      CALL CC5HYPPRT(WORK(KHYPPOL))

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_4HYP                              *
*=====================================================================*
c /* deck CC5HYPPRT */
*=====================================================================*
       SUBROUTINE CC5HYPPRT(HYPERPOL)
*---------------------------------------------------------------------*
*
*    Purpose: print fourth hyperpolarizabilities 
*
*    Written by Christof Haettig in maj 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <cc5rinf.h>
#include <ccroper.h>


      CHARACTER*10 BLANKS
      CHARACTER*80 STRING
      INTEGER IFREQ, IOPER, ISYOLD, ISYTOT, IDX


#if defined (SYS_CRAY)
      REAL HYPERPOL(N5RFREQ,N5ROPER,2)
      REAL HALF
#else
      DOUBLE PRECISION HYPERPOL(N5RFREQ,N5ROPER,2)
      DOUBLE PRECISION HALF
#endif
      PARAMETER (HALF = 0.5D0)

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      BLANKS = '          '
      STRING = ' RESULTS FOR THE FOURTH HYPERPOLARIZABILITIES '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CC5HYPPRT called for an unknown Coupled '//
     &             'Cluster model.')
      END IF

      IF (IPR5HYP.GT.5) THEN
       WRITE(LUPRI,'(/1X,6(1X,A," operator",5X),4X,A,8X,A,/,121("-"))')
     &     'A','B','C','D','E','F','epsilon','(asy. Resp.)'
      ELSE
       WRITE(LUPRI,'(/1X,6(1X,A," operator",5X),4X,A,/,110("-"))')
     &     'A','B','C','D','E','F','epsilon'
      END IF

      ISYOLD = 1

      DO IOPER = 1, N5ROPER
         ISYTOT = 1

         DO IDX = 1, 6
           ISYTOT = MULD2H(ISYTOT,ISYOPR(I5ROP(IOPER,IDX)))
         END DO

         IFREQ = 1
         IF (ISYTOT.EQ.1) THEN
          IF (IPR5HYP.GT.5) THEN
            WRITE(LUPRI,'(/1X,6(A7,F8.4,1X),G16.8," (",G10.3,")")')
     &        (LBLOPR(I5ROP(IOPER,IDX))(1:7),FREQ5(IFREQ,IDX),IDX=1,6),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2)),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)-HYPERPOL(IFREQ,IOPER,2))
          ELSE
            WRITE(LUPRI,'(/1X,6(A7,F8.4,1X),G16.8)')
     &        (LBLOPR(I5ROP(IOPER,IDX))(1:7),FREQ5(IFREQ,IDX),IDX=1,6),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
          END IF
          ISYOLD = 1
         ELSE
          IF (IPR5HYP.GT.5) THEN
            IF (ISYOLD.EQ.1) WRITE(LUPRI,*)
            WRITE(LUPRI,'(1X,6(A7,A8,1X),6X,A,7X," (",5X,A,6X,")")')
     &        (LBLOPR(I5ROP(IOPER,IDX))(1:7),'   -.-  ',IDX=1,6),
     &        '---',
     &        '---'
          ELSE IF (IPR5HYP.GT.0) THEN
            IF (ISYOLD.EQ.1) WRITE(LUPRI,*)
            WRITE(LUPRI,'(1X,6(A7,A8,1X),6X,A,7X)')
     &        (LBLOPR(I5ROP(IOPER,IDX))(1:7),'   -.-  ',IDX=1,6),
     &        '---'
          END IF
          ISYOLD = 0
         END IF

         DO IFREQ = 2, N5RFREQ
          IF (ISYTOT.EQ.1) THEN
           IF (IPR5HYP.GT.5) THEN
            WRITE(LUPRI,'(1X,6(7X,F8.4,1X),G16.8," (",G10.3,")")')
     &        (FREQ5(IFREQ,IDX),IDX=1,6),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2)),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)-HYPERPOL(IFREQ,IOPER,2))
           ELSE
            WRITE(LUPRI,'(1X,6(7X,F8.4,1X),G16.8)')
     &        (FREQ5(IFREQ,IDX),IDX=1,6),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
           END IF
          END IF
         END DO

      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC5HYPPRT                            *
*---------------------------------------------------------------------*
c /* deck CC5R_SETUP */
*=====================================================================*
      SUBROUTINE CC5R_SETUP(MXTRAN3, MXVEC3, MXTRAN4, MXVEC2,
     &                      I0HTRAN,  I0HDOTS,  W0H,  N0HTRAN,
     &                      I1HTRAN,  I1HDOTS,  W1H,  N1HTRAN,
     &                      I0GTRAN,  I0GDOTS,  W0G,  N0GTRAN,
     &                      I1GTRAN,  I1GDOTS,  W1G,  N1GTRAN,
     &                      I1FATRAN, I1FADOTS, W1FA, N1FATRAN,
     &                      IDTRAN,   IDDOTS,   WD,   NDTRAN,
     &                      ICTRAN,   ICDOTS,   WC,   NCTRAN,
     &                      IBTRAN,   IBDOTS,   WB,   NBTRAN,
     &                      IBATRAN,  IBADOTS,  WBA,  NBATRAN,
     &                      I0FTRAN,  I0FDOTS,  W0F,  N0FTRAN,
     &                      IXTRAN,   IXDOTS, WCHI, NXTRAN,
     &                      HYPPOL,   MXHYPPOL, LADD     )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC5R section
*                - list of H^0  matrix transformations 
*                - list of H^A  matrix transformations 
*                - list of G^0  matrix transformations 
*                - list of G^A  matrix transformations 
*                - list of F^A{O}  matrix transformations 
*                - list of D matrix transformations 
*                - list of C matrix transformations 
*                - list of B matrix transformations 
*                - list of B{O} matrix transformations 
*                - list of F^0  matrix transformations 
*                - list of chi3 dot products
*
*     LADD = .FALSE.  --> set lists
*     LADD = .TRUE.   --> add contributions to hyperpolarizabilities
*
*     Written by Christof Haettig, maj 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cc5rinf.h>
#include <cc5perm.h>
#include <ccroper.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC5R_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ORD1, ORD2, ORD3
      PARAMETER ( ORD1 = 6 ) ! order (nb. of operators)
      PARAMETER ( ORD2 = (ORD1-1)*ORD1/2 )
      PARAMETER ( ORD3 = (ORD1-2)*ORD2/3 )

      LOGICAL LADD

      INTEGER IOP(ORD1), IR1(ORD1), IL1(ORD1), ISY(ORD1)
      INTEGER IR2(ORD2), IL2(ORD2)
      INTEGER IR3(ORD3), IX3(ORD3)

      INTEGER MXVEC3, MXTRAN3, MXVEC2, MXTRAN4, MXHYPPOL

      INTEGER I0HTRAN(5,MXTRAN4)
      INTEGER I0HDOTS(MXVEC2,MXTRAN4)
      INTEGER I1HTRAN(5,MXTRAN4)
      INTEGER I1HDOTS(MXVEC2,MXTRAN4)

      INTEGER I0GTRAN(4,MXTRAN4)
      INTEGER I0GDOTS(MXVEC2,MXTRAN4)
      INTEGER I1GTRAN(4,MXTRAN4)
      INTEGER I1GDOTS(MXVEC2,MXTRAN4)

      INTEGER I1FATRAN(MXDIM_FATRAN,MXTRAN4)
      INTEGER I1FADOTS(MXVEC2,MXTRAN4)

      INTEGER IDTRAN(5,MXTRAN4)
      INTEGER IDDOTS(MXVEC2,MXTRAN4)
      INTEGER ICTRAN(4,MXTRAN4)
      INTEGER ICDOTS(MXVEC2,MXTRAN4)
      INTEGER IBTRAN(3,MXTRAN4)
      INTEGER IBDOTS(MXVEC2,MXTRAN4)
      INTEGER IBATRAN(MXDIM_BATRAN,MXTRAN4)
      INTEGER IBADOTS(MXVEC2,MXTRAN4)

      INTEGER I0FTRAN(3,MXTRAN3)
      INTEGER I0FDOTS(MXVEC3,MXTRAN3)

      INTEGER IXTRAN(MXTRAN3)
      INTEGER IXDOTS(MXVEC3,MXTRAN3)

      INTEGER N0HTRAN, N0GTRAN, N0FTRAN
      INTEGER N1HTRAN, N1GTRAN, N1FATRAN
      INTEGER NDTRAN,  NCTRAN,  NBTRAN,  NBATRAN, NXTRAN
      INTEGER N5RRESF
      INTEGER MXVH0, MXVH1, MXVG0, MXVG1, MXVF0, MXVFA1
      INTEGER MXVD,  MXVC,  MXVB,  MXVBA, MXCHI

      INTEGER IVEC, ITRAN, I, IDX, IDXAB, IDXA, IDXB, IDXC, IDXABC
      INTEGER ISYML, ISYM1, ISYM2, ISYM3, ISYTOT
      INTEGER IFREQ, IOPER
      INTEGER P, ISIGN

#if defined (SYS_CRAY)
      REAL SIGN, SUM, SSUM
      REAL HYPPOL(MXHYPPOL)
      REAL H0CON(45), H1CON(60)
      REAL G0CON(15), G1CON(90)
      REAL F0CON(10), F1ACON(90), CHICON(20)
      REAL DCON(15), CCON(90), BCON(45), BACON(180)
      REAL W0H(MXVEC2,MXTRAN4), W1H(MXVEC2,MXTRAN4)
      REAL W0G(MXVEC2,MXTRAN4), W1G(MXVEC2,MXTRAN4)
      REAL W1FA(MXVEC2,MXTRAN4), WD(MXVEC2,MXTRAN4)
      REAL WC(MXVEC2,MXTRAN4), WB(MXVEC2,MXTRAN4)
      REAL WBA(MXVEC2,MXTRAN4), W0F(MXVEC3,MXTRAN3)
      REAL WCHI(MXVEC3,MXTRAN3)
#else
      DOUBLE PRECISION SIGN, SUM, SSUM
      DOUBLE PRECISION HYPPOL(MXHYPPOL)
      DOUBLE PRECISION H0CON(45), H1CON(60)
      DOUBLE PRECISION G0CON(15), G1CON(90)
      DOUBLE PRECISION F0CON(10), F1ACON(90), CHICON(20)
      DOUBLE PRECISION DCON(15), CCON(90), BCON(45), BACON(180)
      DOUBLE PRECISION W0H(MXVEC2,MXTRAN4), W1H(MXVEC2,MXTRAN4)
      DOUBLE PRECISION W0G(MXVEC2,MXTRAN4), W1G(MXVEC2,MXTRAN4)
      DOUBLE PRECISION W1FA(MXVEC2,MXTRAN4), WD(MXVEC2,MXTRAN4)
      DOUBLE PRECISION WC(MXVEC2,MXTRAN4), WB(MXVEC2,MXTRAN4)
      DOUBLE PRECISION WBA(MXVEC2,MXTRAN4), W0F(MXVEC3,MXTRAN3)
      DOUBLE PRECISION WCHI(MXVEC3,MXTRAN3)
#endif


* external functions:
      INTEGER IR3TAMP
      INTEGER IR2TAMP
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER ICHI3


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN
        N0HTRAN  = 0
        N1HTRAN  = 0
        N0GTRAN  = 0
        N1GTRAN  = 0
        N1FATRAN = 0
        NDTRAN   = 0
        NCTRAN   = 0
        NBTRAN   = 0
        NBATRAN  = 0
        N0FTRAN  = 0
        NXTRAN    = 0
      END IF

      MXVH0  = 0
      MXVH1  = 0
      MXVG0  = 0
      MXVG1  = 0
      MXVFA1 = 0
      MXVD   = 0
      MXVC   = 0
      MXVB   = 0
      MXVBA  = 0
      MXVF0  = 0
      MXCHI  = 0

      N5RRESF  = 0

      CALL DZERO(HYPPOL,MXHYPPOL)
 
*---------------------------------------------------------------------*
* start loop over all requested third hyperpolarizabilities
*---------------------------------------------------------------------*
 
      DO IOPER = 1, N5ROPER
        ISYTOT = 1
        DO IDX = 1, 6
          IOP(IDX) = I5ROP(IOPER,IDX)
          ISY(IDX) = ISYOPR(IOP(IDX))
          ISYTOT = MULD2H(ISYTOT,ISY(IDX))
        END DO

        IF (LOCDBG) THEN
          CALL DCOPY (45,999.99d0,0,H0CON,1)
          CALL DCOPY (60,999.99d0,0,H1CON,1)
          CALL DCOPY (15,999.99d0,0,G0CON,1)
          CALL DCOPY (90,999.99d0,0,G1CON,1)
          CALL DCOPY (10,999.99d0,0,F0CON,1)
          CALL DCOPY (90,999.99d0,0,F1ACON,1)
          CALL DCOPY (20,999.99d0,0,CHICON,1)
          CALL DCOPY (60,999.99d0,0,DCON,1)
          CALL DCOPY (90,999.99d0,0,CCON,1)
          CALL DCOPY (45,999.99d0,0,BCON,1)
          CALL DCOPY (180,999.99d0,0,BACON,1)
        END IF
 
      IF (ISYTOT .EQ. 1) THEN

      DO IFREQ = 1, N5RFREQ

        N5RRESF = N5RRESF + 1
      
      DO ISIGN = 1, -1, -2
        SIGN = DBLE(ISIGN)
  
        DO IDX = 1, 6
         IL1(IDX) = IL1ZETA(LBLOPR(IOP(IDX)),.FALSE.,
     &                       SIGN*FREQ5(IFREQ,IDX),ISYM1)
         IR1(IDX) = IR1TAMP(LBLOPR(IOP(IDX)),.FALSE.,
     &                       SIGN*FREQ5(IFREQ,IDX),ISYM1)
        END DO

       IDXAB = 0
       DO IDXB = 2, 6
       DO IDXA = 1, IDXB-1
        IDXAB = IDXAB + 1
        IR2(IDXAB) = 
     &   IR2TAMP(LBLOPR(IOP(IDXA)),.FALSE.,SIGN*FREQ5(IFREQ,IDXA),ISYM1,
     &           LBLOPR(IOP(IDXB)),.FALSE.,SIGN*FREQ5(IFREQ,IDXB),ISYM2)

        IL2(IDXAB) = 
     &   IL2ZETA(LBLOPR(IOP(IDXA)),        SIGN*FREQ5(IFREQ,IDXA),ISYM1,
     &           LBLOPR(IOP(IDXB)),        SIGN*FREQ5(IFREQ,IDXB),ISYM2)
       END DO
       END DO

        IDXABC = 0
        DO IDXC = 3, 6
        DO IDXB = 2, IDXC-1
        DO IDXA = 1, IDXB-1
          IDXABC = IDXABC + 1
          IR3(IDXABC) = 
     &         IR3TAMP( LBLOPR(IOP(IDXA)),SIGN*FREQ5(IFREQ,IDXA),ISYM1,
     &                  LBLOPR(IOP(IDXB)),SIGN*FREQ5(IFREQ,IDXB),ISYM2,
     &                  LBLOPR(IOP(IDXC)),SIGN*FREQ5(IFREQ,IDXC),ISYM3)

          IX3(IDXABC) = 
     &           ICHI3( LBLOPR(IOP(IDXA)),SIGN*FREQ5(IFREQ,IDXA),ISYM1,
     &                  LBLOPR(IOP(IDXB)),SIGN*FREQ5(IFREQ,IDXB),ISYM2,
     &                  LBLOPR(IOP(IDXC)),SIGN*FREQ5(IFREQ,IDXC),ISYM3)
        END DO
        END DO
        END DO

*---------------------------------------------------------------------*
* set up list of H^0 matrix transformations, 45 permutations
*---------------------------------------------------------------------*
      DO P = 1, 45
        CALL CC_SETH2112(I0HTRAN,I0HDOTS,MXTRAN4,MXVEC2,
     &    0,IR2(KP2(P)),IR1(K1(P)),IR1(K2(P)),IR2(KP3(P)),ITRAN,IVEC)
        N0HTRAN  = MAX(N0HTRAN,ITRAN)
        MXVH0    = MAX(MXVH0,IVEC)
        H0CON(P) = W0H(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of H^A matrix transformations, 60 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CC_SETH1112(I1HTRAN,I1HDOTS,MXTRAN4,MXVEC2,
     &    IL1(J1(P)),IR1(J2(P)),IR1(J3(P)),IR1(J4(P)),IR2(JP(P)),
     &    ITRAN,IVEC)
        N1HTRAN  = MAX(N1HTRAN,ITRAN)
        MXVH1    = MAX(MXVH1,IVEC)
        H1CON(P) = W1H(IVEC,ITRAN)

        CALL CC_SETH1112(I1HTRAN,I1HDOTS,MXTRAN4,MXVEC2,
     &    IL1(J2(P)),IR1(J1(P)),IR1(J3(P)),IR1(J4(P)),IR2(JP(P)),
     &    ITRAN,IVEC)
        N1HTRAN     = MAX(N1HTRAN,ITRAN)
        MXVH1       = MAX(MXVH1,IVEC)
        H1CON(P+15) = W1H(IVEC,ITRAN)

        CALL CC_SETH1112(I1HTRAN,I1HDOTS,MXTRAN4,MXVEC2,
     &    IL1(J3(P)),IR1(J1(P)),IR1(J2(P)),IR1(J4(P)),IR2(JP(P)),
     &    ITRAN,IVEC)
        N1HTRAN     = MAX(N1HTRAN,ITRAN)
        MXVH1       = MAX(MXVH1,IVEC)
        H1CON(P+30) = W1H(IVEC,ITRAN)

        CALL CC_SETH1112(I1HTRAN,I1HDOTS,MXTRAN4,MXVEC2,
     &    IL1(J4(P)),IR1(J1(P)),IR1(J2(P)),IR1(J3(P)),IR2(JP(P)),
     &    ITRAN,IVEC)
        N1HTRAN     = MAX(N1HTRAN,ITRAN)
        MXVH1       = MAX(MXVH1,IVEC)
        H1CON(P+45) = W1H(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of G^0 matrix transformations, 15 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CCQR_SETG(I0GTRAN,I0GDOTS,MXTRAN4,MXVEC2,   
     &              0,IR2(IP1(P)),IR2(IP2(P)),IR2(IP3(P)),ITRAN,IVEC)
        N0GTRAN  = MAX(N0GTRAN,ITRAN)
        MXVG0    = MAX(MXVG0,IVEC)
        G0CON(P) = W0G(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of G^A matrix transformations, 90 permutations
*---------------------------------------------------------------------*
      DO P = 1, 45
        CALL CC_SETG212(I1GTRAN,I1GDOTS,MXTRAN4,MXVEC2,   
     &     IL1(K1(P)),IR2(KP2(P)),IR1(K2(P)),IR2(KP3(P)),ITRAN,IVEC)
        N1GTRAN  = MAX(N1GTRAN,ITRAN)
        MXVG1    = MAX(MXVG1,IVEC)
        G1CON(P) = W1G(IVEC,ITRAN)

        CALL CC_SETG212(I1GTRAN,I1GDOTS,MXTRAN4,MXVEC2,   
     &     IL1(K2(P)),IR2(KP2(P)),IR1(K1(P)),IR2(KP3(P)),ITRAN,IVEC)
        N1GTRAN     = MAX(N1GTRAN,ITRAN)
        MXVG1       = MAX(MXVG1,IVEC)
        G1CON(P+45) = W1G(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of F^A{O} matrix transformations, 90 permutations
*---------------------------------------------------------------------*
      DO P = 1, 45
        CALL CCQR_SETFA(I1FATRAN,I1FADOTS,MXTRAN4,MXVEC2,
     &     IL1(K1(P)),IOP(K2(P)),IR2(KP2(P)),IR2(KP3(P)),ITRAN,IVEC)
        N1FATRAN  = MAX(N1FATRAN,ITRAN)
        MXVFA1    = MAX(MXVFA1,IVEC)
        F1ACON(P) = W1FA(IVEC,ITRAN)

        CALL CCQR_SETFA(I1FATRAN,I1FADOTS,MXTRAN4,MXVEC2,
     &     IL1(K2(P)),IOP(K1(P)),IR2(KP2(P)),IR2(KP3(P)),ITRAN,IVEC)
        N1FATRAN     = MAX(N1FATRAN,ITRAN)
        MXVFA1       = MAX(MXVFA1,IVEC)
        F1ACON(P+45) = W1FA(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of D matrix transformations, 15 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CC_SETD1111(IDTRAN,IDDOTS,MXTRAN4,MXVEC2, 
     &         IL2(JP(P)),IR1(J1(P)),IR1(J2(P)),IR1(J3(P)),IR1(J4(P)),
     &         ITRAN,IVEC)
        NDTRAN  = MAX(NDTRAN,ITRAN)
        MXVD    = MAX(MXVD,IVEC)
        DCON(P) = WD(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of C matrix transformations, 90 permutations
*---------------------------------------------------------------------*
      DO P = 1, 45
        CALL CC_SETC211(ICTRAN,ICDOTS,MXTRAN4,MXVEC2, 
     &       IL2(KP2(P)),IR2(KP3(P)),IR1(K1(P)),IR1(K2(P)),ITRAN,IVEC)
        NCTRAN  = MAX(NCTRAN,ITRAN)
        MXVC    = MAX(MXVC,IVEC)
        CCON(P) = WC(IVEC,ITRAN)

        CALL CC_SETC211(ICTRAN,ICDOTS,MXTRAN4,MXVEC2, 
     &       IL2(KP3(P)),IR2(KP2(P)),IR1(K1(P)),IR1(K2(P)),ITRAN,IVEC)
        NCTRAN     = MAX(NCTRAN,ITRAN)
        MXVC       = MAX(MXVC,IVEC)
        CCON(P+45) = WC(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of B matrix transformations, 45 permutations
*---------------------------------------------------------------------*
      DO P = 1, 45
        CALL CC_SETB11(IBTRAN,IBDOTS,MXTRAN4,MXVEC2, 
     &       IL2(KP1(P)),IR2(KP2(P)),IR2(KP3(P)),ITRAN,IVEC)
        NBTRAN  = MAX(NBTRAN,ITRAN)
        MXVB    = MAX(MXVB,IVEC)
        BCON(P) = WB(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of B{O} matrix transformations, 180 permutations
*---------------------------------------------------------------------*
      DO P = 1, 45
        CALL CC_SETBA12(IBATRAN,IBADOTS,MXTRAN4,MXVEC2, 
     &       IL2(KP2(P)),IOP(K1(P)),IR2(KP3(P)),IR1(K2(P)),ITRAN,IVEC)
        NBATRAN  = MAX(NBATRAN,ITRAN)
        MXVBA    = MAX(MXVBA,IVEC)
        BACON(P) = WBA(IVEC,ITRAN)

        CALL CC_SETBA12(IBATRAN,IBADOTS,MXTRAN4,MXVEC2, 
     &       IL2(KP3(P)),IOP(K1(P)),IR2(KP2(P)),IR1(K2(P)),ITRAN,IVEC)
        NBATRAN     = MAX(NBATRAN,ITRAN)
        MXVBA       = MAX(MXVBA,IVEC)
        BACON(P+45) = WBA(IVEC,ITRAN)

        CALL CC_SETBA12(IBATRAN,IBADOTS,MXTRAN4,MXVEC2, 
     &       IL2(KP2(P)),IOP(K2(P)),IR2(KP3(P)),IR1(K1(P)),ITRAN,IVEC)
        NBATRAN     = MAX(NBATRAN,ITRAN)
        MXVBA       = MAX(MXVBA,IVEC)
        BACON(P+90) = WBA(IVEC,ITRAN)

        CALL CC_SETBA12(IBATRAN,IBADOTS,MXTRAN4,MXVEC2, 
     &       IL2(KP3(P)),IOP(K2(P)),IR2(KP2(P)),IR1(K1(P)),ITRAN,IVEC)
        NBATRAN      = MAX(NBATRAN,ITRAN)
        MXVBA        = MAX(MXVBA,IVEC)
        BACON(P+135) = WBA(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of F^0 matrix transformations, 10 permutations
*---------------------------------------------------------------------*
      DO P = 1, 10
        CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN3,MXVEC3, 
     &                 0,IR3(IT1(P)),IR3(IT2(P)),ITRAN,IVEC)
        N0FTRAN  = MAX(N0FTRAN,ITRAN)
        MXVF0    = MAX(MXVF0,IVEC)
        F0CON(P) = W0F(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of chi3 vector dot products, 20 permutations
*---------------------------------------------------------------------*
      DO P = 1, 10
        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN3,MXVEC3,
     &                 IX3(IT1(P)),IR3(IT2(P)),ITRAN,IVEC)
        NXTRAN    = MAX(NXTRAN,ITRAN)
        MXCHI     = MAX(MXCHI,IVEC)
        CHICON(P) = WCHI(IVEC,ITRAN)

        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN3,MXVEC3,
     &                 IX3(IT2(P)),IR3(IT1(P)),ITRAN,IVEC)
        NXTRAN       = MAX(NXTRAN,ITRAN)
        MXCHI        = MAX(MXCHI,IVEC)
        CHICON(P+10) = WCHI(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* add all 660 contributions to the hyperpolarizabilities:
*---------------------------------------------------------------------*
      IF (LADD) THEN
        IDX = N5RFREQ*(IOPER-1) + IFREQ
        IF (ISIGN.EQ.-1) IDX = IDX + N5RFREQ*N5ROPER

        DO P = 1, 10
          HYPPOL(IDX) = HYPPOL(IDX) + F0CON(P) + CHICON(P)+CHICON(P+10)
        END DO

        DO P = 1, 15
          HYPPOL(IDX) = HYPPOL(IDX) + G0CON(P) + DCON(P) +
     &              H1CON(P) + H1CON(P+15) + H1CON(P+30) + H1CON(P+45)
        END DO

        DO P = 1, 45
         HYPPOL(IDX) = HYPPOL(IDX) + 
     &    H0CON(P)  + G1CON(P)     + G1CON(P+45) +
     &    F1ACON(P) + F1ACON(P+45) + CCON(P) + CCON(P+45) + BCON(P) +
     &    BACON(P)  + BACON(P+45)  + BACON(P+90) + BACON(P+135)
        END DO

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'IOP:',IOP
          WRITE (LUPRI,*) 'IDX',IDX
          WRITE (LUPRI,*) 'HYPPOL:',HYPPOL(IDX)

          SSUM= 0.0d0

          SUM = 0.0d0
          DO I = 1, 10
            SUM = SUM + F0CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'F0CON:',F0CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 20
            SUM = SUM + CHICON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'CHICON:',CHICON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 15
            SUM = SUM + G0CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'G0CON:',G0CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 15
            SUM = SUM + DCON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'DCON:',DCON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 60
            SUM = SUM + H1CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'H1CON:',H1CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 45
            SUM = SUM + H0CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'H0CON:',H0CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 90
            SUM = SUM + G1CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'G1CON:',G1CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 90
            SUM = SUM + F1ACON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'F1ACON:',F1ACON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 90
            SUM = SUM + CCON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'CCON:',CCON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 45
            SUM = SUM + BCON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'BCON:',BCON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 180
            SUM = SUM + BACON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'BACON:',BACON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          WRITE (LUPRI,*) 'HYPPOL:',HYPPOL(IDX)
        END IF

      END IF

*---------------------------------------------------------------------*
* end loop over all requested hyperpolarizabilities
*---------------------------------------------------------------------*
      END DO
      END DO
      END IF
      END DO

*---------------------------------------------------------------------*
* print statistics and lists: 
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN

* general statistics:
      WRITE(LUPRI,'(/,/3X,A,I4,A)') 'For the requested',N5RRESF,
     &      ' pentic response functions '
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',N0HTRAN,  ' H^0 matrix transformations ',
     &   ' - ',N1HTRAN,  ' H^A matrix transformations ',
     &   ' - ',N0GTRAN,  ' G^0 matrix transformations ',
     &   ' - ',N1GTRAN,  ' G^A matrix transformations ',
     &   ' - ',N1FATRAN, ' F^A{O} matrix transformations ',
     &   ' - ',NDTRAN,   ' D matrix transformations ',
     &   ' - ',NCTRAN,   ' C matrix transformations ',
     &   ' - ',NBTRAN,   ' B matrix transformations ',
     &   ' - ',NBATRAN,  ' B{O} matrix transformations ',
     &   ' - ',N0FTRAN,  ' F^0 matrix transformations ',
     &   ' - ',NXTRAN,    ' dot products with chi3 vectors' 
      WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'


* H^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of H^0 matrix transformations:'
        DO ITRAN = 1, N0HTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I0HTRAN(I,ITRAN),I=1,4),(I0HDOTS(I,ITRAN),I=1,MXVH0)
        END DO
        WRITE (LUPRI,*)
      END IF

* H^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of H^A matrix transformations:'
        DO ITRAN = 1, N1HTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I1HTRAN(I,ITRAN),I=1,4),(I1HDOTS(I,ITRAN),I=1,MXVH1)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^0 matrix transformations:'
        DO ITRAN = 1, N0GTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I0GTRAN(I,ITRAN),I=1,3),(I0GDOTS(I,ITRAN),I=1,MXVG0)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^A matrix transformations:'
        DO ITRAN = 1, N1GTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I1GTRAN(I,ITRAN),I=1,3),(I1GDOTS(I,ITRAN),I=1,MXVG1)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^0 matrix transformations:'
        DO ITRAN = 1, N0FTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (I0FTRAN(I,ITRAN),I=1,2),(I0FDOTS(I,ITRAN),I=1,MXVF0)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^A{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^A{O} matrix transformations:'
        DO ITRAN = 1, N1FATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (I1FATRAN(I,ITRAN),I=1,5),(I1FADOTS(I,ITRAN),I=1,MXVFA1)
        END DO
        WRITE (LUPRI,*)
      END IF

* D matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of D matrix transformations:'
        DO ITRAN = 1, NDTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (IDTRAN(I,ITRAN),I=1,4),(IDDOTS(I,ITRAN),I=1,MXVD)
        END DO
        WRITE (LUPRI,*)
      END IF

* C matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of C matrix transformations:'
        DO ITRAN = 1, NCTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (ICTRAN(I,ITRAN),I=1,3),(ICDOTS(I,ITRAN),I=1,MXVC)
        END DO
        WRITE (LUPRI,*)
      END IF

* B matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of B matrix transformations:'
        DO ITRAN = 1, NBTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IBTRAN(I,ITRAN),I=1,2),(IBDOTS(I,ITRAN),I=1,MXVB)
        END DO
        WRITE (LUPRI,*)
      END IF

* B{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of B{A} matrix transformations:'
        DO ITRAN = 1, NBATRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IBATRAN(I,ITRAN),I=1,3),(IBADOTS(I,ITRAN),I=1,MXVBA)
        END DO
        WRITE (LUPRI,*)
      END IF

* chi vector dot products:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of chi vector dot products:'
        DO ITRAN = 1, NXTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &     IXTRAN(ITRAN),(IXDOTS(I,ITRAN),I=1,MXCHI)
        END DO
        WRITE (LUPRI,*)
      END IF

      END IF ! (.NOT. LADD)
*---------------------------------------------------------------------*
* return
*---------------------------------------------------------------------*

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC5R_SETUP                           *
*---------------------------------------------------------------------*
c /* deck CC_SETBA12 */
*=====================================================================*
       SUBROUTINE CC_SETBA12(IBATRAN,IBADOTS,MXTRAN,MXVEC,
     &                       IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintains a list of dot products of B{O} matrix 
*             transformations with two right amplitude vectors:
*                        L * (B{O} * T^A * T^B)
*             assumes that T^A and T^B belong to different lists
*
*             IBATRAN - list of F matrix transformations
*             IBADOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IBADOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             IOPER  - index of property operator 
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*
*             ITRAN - index in IBATRAN list
*             IVEC  - second index in IBADOTS list
*
*    Written by Christof Haettig, june 1997.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IBATRAN(MXDIM_BATRAN,MXTRAN)
      INTEGER IBADOTS(MXVEC,MXTRAN)

      LOGICAL LFNDB
      INTEGER IZETAV, IOPER, ITAMPA, ITAMPB
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LBATST, LBAEND
      INTEGER IB, IA, IO
      LBATST(ITRAN,IO,IA,IB) = IBATRAN(1,ITRAN).EQ.IO 
     &       .AND. IBATRAN(2,ITRAN).EQ.IA .AND. IBATRAN(3,ITRAN).EQ.IB 
      LBAEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &      (IBATRAN(1,ITRAN)+IBATRAN(2,ITRAN)+IBATRAN(3,ITRAN)).LE.0 


*---------------------------------------------------------------------*
* set up list of B{A} matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDB  = LBATST(ITRAN,IOPER,ITAMPA,ITAMPB)

      DO WHILE ( .NOT. (LFNDB.OR.LBAEND(ITRAN)))
       ITRAN = ITRAN + 1
       LFNDB  = LBATST(ITRAN,IOPER,ITAMPA,ITAMPB)
      END DO

      IF (.NOT.LFNDB) THEN
        IBATRAN(1,ITRAN) = IOPER
        IBATRAN(2,ITRAN) = ITAMPA
        IBATRAN(3,ITRAN) = ITAMPB
        IBATRAN(4,ITRAN) = 0
        IBATRAN(5,ITRAN) = 0
        IBATRAN(6,ITRAN) = 0
        IBATRAN(7,ITRAN) = 0
        IBATRAN(8,ITRAN) = 0
        ITAMP = IZETAV
      ELSE 
        IF (LFNDB) ITAMP = IZETAV
      END IF

      IVEC = 1
      DO WHILE (IBADOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &            IBADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IBADOTS(IVEC,ITRAN) = ITAMP

C     WRITE (LUPRI,*) 'CC_SETBA12>',IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC
*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC:',
     &              IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC
        IDX = 1
        DO WHILE ( .NOT. LBAEND(IDX) )
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') 'CC_SETBA12>',
     &       (IBATRAN(I,IDX),I=1,4),(IBADOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETBA12.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETBA12                           *
*---------------------------------------------------------------------*
c /* deck CC_SETB11 */
*=====================================================================*
       SUBROUTINE CC_SETB11(IBTRAN,IBDOTS,MXTRAN,MXVEC,
     &                      IZETAV,ITAMPA,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of B matrix transformations
*             assumes that T^A and T^B members of the same list
*
*             IBTRAN - list of B matrix transformations
*             IBDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IBDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*
*             ITRAN - index in IBTRAN list
*             IVEC  - second index in IBDOTS list
*
*    Written by Christof Haettig, june 1997.
*
*=====================================================================*
      IMPLICIT NONE  

#include <priunit.h>
      INTEGER MXVEC, MXTRAN
      INTEGER IBTRAN(3,MXTRAN)
      INTEGER IBDOTS(MXVEC,MXTRAN)

      LOGICAL LFND
      INTEGER IZETAV, ITAMPA, ITAMPB
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LBTST, LBEND
      INTEGER IB, IA
      LBTST(ITRAN,IA,IB) = 
     &      ( IBTRAN(1,ITRAN).EQ.IA .AND. IBTRAN(2,ITRAN).EQ.IB ) .OR.
     &      ( IBTRAN(1,ITRAN).EQ.IB .AND. IBTRAN(2,ITRAN).EQ.IA )
      LBEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (IBTRAN(1,ITRAN)+IBTRAN(2,ITRAN)).LE.0

*---------------------------------------------------------------------*
* set up list of B matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFND = LBTST(ITRAN,ITAMPA,ITAMPB)

      DO WHILE ( .NOT. (LFND.OR.LBEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFND = LBTST(ITRAN,ITAMPA,ITAMPB)
      END DO

      IF (.NOT.LFND) THEN
        IBTRAN(1,ITRAN) = ITAMPA
        IBTRAN(2,ITRAN) = ITAMPB
        IBTRAN(3,ITRAN) = 0
        ITAMP = IZETAV
      ELSE 
        IF (LFND) ITAMP = IZETAV
      END IF

      IVEC = 1
      DO WHILE (IBDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           IBDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IBDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB:',ITAMPA,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LBEND(IDX) )
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') 'CC_SETB11>',
     &       (IBTRAN(I,IDX),I=1,3),(IBDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETB11.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETB11                            *
*---------------------------------------------------------------------*
c /* deck CC_SETD1111 */
*=====================================================================*
       SUBROUTINE CC_SETD1111(IDTRAN,IDDOTS,MXTRAN,MXVEC,IZETAV,
     &                        ITAMPA,ITAMPB,ITAMPC,ITAMPD,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of D matrix transformations
*             assumes that T^A, T^B, T^C, T^D members of the same list
*
*             IDTRAN - list of D matrix transformations
*             IDDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IDDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*             ITAMPC - index of amplitude vector C
*             ITAMPD - index of amplitude vector D
*
*             ITRAN - index in IDTRAN list
*             IVEC  - second index in IDDOTS list
*
*    Written by Christof Haettig, june 1997.
*
*=====================================================================*
      IMPLICIT NONE  

#include <priunit.h>
      INTEGER MXVEC, MXTRAN
      INTEGER IDTRAN(5,MXTRAN)
      INTEGER IDDOTS(MXVEC,MXTRAN)

      LOGICAL LFND
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC, ITAMPD
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LDTST, LDXTST, LDEND
      INTEGER IA, IB, IC, ID, IA1, IB1, IC1, ID1

      LDXTST(ITRAN,IA1,IB1,IC1,ID1) = 
     &        IDTRAN(1,ITRAN).EQ.IA1 .AND. IDTRAN(2,ITRAN).EQ.IB1 .AND.
     &        IDTRAN(3,ITRAN).EQ.IC1 .AND. IDTRAN(4,ITRAN).EQ.ID1  

      LDTST(ITRAN,IA,IB,IC,ID) = 
     &   LDXTST(ITRAN,IA,IB,IC,ID) .OR. LDXTST(ITRAN,IA,IB,ID,IC) .OR.
     &   LDXTST(ITRAN,IB,IA,IC,ID) .OR. LDXTST(ITRAN,IB,IA,ID,IC) .OR.
     &   LDXTST(ITRAN,IA,IC,IB,ID) .OR. LDXTST(ITRAN,IA,IC,ID,IB) .OR.
     &   LDXTST(ITRAN,IC,IA,IB,ID) .OR. LDXTST(ITRAN,IC,IA,ID,IB) .OR.
     &   LDXTST(ITRAN,IA,ID,IC,IB) .OR. LDXTST(ITRAN,IA,ID,IB,IC) .OR.
     &   LDXTST(ITRAN,ID,IA,IC,IB) .OR. LDXTST(ITRAN,ID,IA,IB,IC) .OR.
     &   LDXTST(ITRAN,ID,IB,IC,IA) .OR. LDXTST(ITRAN,ID,IB,IA,IC) .OR.
     &   LDXTST(ITRAN,IB,ID,IC,IA) .OR. LDXTST(ITRAN,IB,ID,IA,IC) .OR.
     &   LDXTST(ITRAN,ID,IC,IB,IA) .OR. LDXTST(ITRAN,ID,IC,IA,IB) .OR.
     &   LDXTST(ITRAN,IC,ID,IB,IA) .OR. LDXTST(ITRAN,IC,ID,IA,IB) .OR.
     &   LDXTST(ITRAN,ID,IA,IC,IB) .OR. LDXTST(ITRAN,ID,IA,IB,IC) .OR.
     &   LDXTST(ITRAN,IA,ID,IC,IB) .OR. LDXTST(ITRAN,IA,ID,IB,IC) 

      LDEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (IDTRAN(1,ITRAN)+IDTRAN(2,ITRAN)+
     &    IDTRAN(3,ITRAN)+IDTRAN(4,ITRAN)).LE.0

*---------------------------------------------------------------------*
* set up list of D matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFND = LDTST(ITRAN,ITAMPA,ITAMPB,ITAMPC,ITAMPD)

      DO WHILE ( .NOT. (LFND.OR.LDEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFND = LDTST(ITRAN,ITAMPA,ITAMPB,ITAMPC,ITAMPD)
      END DO

      IF (.NOT.LFND) THEN
        IDTRAN(1,ITRAN) = ITAMPA
        IDTRAN(2,ITRAN) = ITAMPB
        IDTRAN(3,ITRAN) = ITAMPC
        IDTRAN(4,ITRAN) = ITAMPD
        IDTRAN(5,ITRAN) = 0
        ITAMP = IZETAV
      ELSE 
        ITAMP = IZETAV
      END IF

      IVEC = 1
      DO WHILE (IDDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           IDDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IDDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC,ITAMPD:',
     &              ITAMPA,ITAMPB,ITAMPC,ITAMPD
        IDX = 1
        DO WHILE ( .NOT. LDEND(IDX) )
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') 'CC_SETD1111>',
     &       (IDTRAN(I,IDX),I=1,5),(IDDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETD1111.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*             END OF SUBROUTINE CC_SETD1111                           *
*---------------------------------------------------------------------*
c /* deck CC_SETH2112 */
*=====================================================================*
      SUBROUTINE CC_SETH2112(IHTRAN,IHDOTS,MXTRAN,MXVEC,IZETAV,
     &                       ITAMPA,ITAMPB,ITAMPC,ITAMPD,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of H matrix 
*             transformations with right amplitude vectors:  
*                       (Z*H*T^A*T^B*T^C) * T^D
*    N.B.: assumes that T^A and T^D vectors are members of one list
*          and that T^B and T^C vectors members of a second list
*
*             IHTRAN - list of H matrix transformations
*             IHDOTS - list of vectors it should be dottet on
*             
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimesion for IHDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*             ITAMPC - index of amplitude vector C
*             ITAMPD - index of amplitude vector D
*
*             ITRAN - index in IHTRAN list
*             IVEC  - second index in IHDOTS list
*
*    Written by Christof Haettig, june 1997.
*
*=====================================================================*
      IMPLICIT NONE  

#include <priunit.h>
      INTEGER MXVEC, MXTRAN
      INTEGER IHTRAN(5,MXTRAN)
      INTEGER IHDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDA, LFNDD
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC, ITAMPD
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LHTST, LHEND
      INTEGER IL, IA, IB,IC
      LHTST(ITRAN,IL,IA,IB,IC) = IHTRAN(1,ITRAN).EQ.IL .AND. (
     &     (IHTRAN(2,ITRAN).EQ.IA .AND. IHTRAN(3,ITRAN).EQ.IB
     &                            .AND. IHTRAN(4,ITRAN).EQ.IC) .OR. 
     &     (IHTRAN(2,ITRAN).EQ.IA .AND. IHTRAN(3,ITRAN).EQ.IC
     &                            .AND. IHTRAN(4,ITRAN).EQ.IB) )
      LHEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &               (IHTRAN(1,ITRAN)+IHTRAN(2,ITRAN)+
     &                IHTRAN(3,ITRAN)+IHTRAN(4,ITRAN) ).LE.0

*---------------------------------------------------------------------*
* set up list of H matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDA = LHTST(ITRAN,IZETAV,ITAMPD,ITAMPB,ITAMPC)
      LFNDD = LHTST(ITRAN,IZETAV,ITAMPA,ITAMPB,ITAMPC)

      DO WHILE ( .NOT. (LFNDA.OR.LFNDD.OR.LHEND(ITRAN)) )
        ITRAN = ITRAN + 1
        LFNDA = LHTST(ITRAN,IZETAV,ITAMPD,ITAMPB,ITAMPC)
        LFNDD = LHTST(ITRAN,IZETAV,ITAMPA,ITAMPB,ITAMPC)
      END DO

      IF (.NOT.(LFNDA.OR.LFNDD)) THEN
        IHTRAN(1,ITRAN) = IZETAV
        IHTRAN(2,ITRAN) = ITAMPA
        IHTRAN(3,ITRAN) = ITAMPB
        IHTRAN(4,ITRAN) = ITAMPC
        IHTRAN(5,ITRAN) = 0
        ITAMP = ITAMPD
      ELSE 
        IF (LFNDA) ITAMP = ITAMPA
        IF (LFNDD) ITAMP = ITAMPD
      END IF

      IVEC = 1
      DO WHILE (IHDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           IHDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IHDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'IZETAV,ITAMPA,ITAMPB,ITAMPC,ITAMPD:',
     &              IZETAV,ITAMPA,ITAMPB,ITAMPC,ITAMPD
        IDX = 1
        DO WHILE (.NOT. LHEND(IDX))
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') 'CC_SETH2112>',
     &       (IHTRAN(I,IDX),I=1,4),(IHDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETH2112.')
      END IF
      
      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETH2112                          *
*---------------------------------------------------------------------*
c /* deck CC_SETC211 */
*=====================================================================*
       SUBROUTINE CC_SETC211(ICTRAN,ICDOTS,MXTRAN,MXVEC,
     &                       IZETAV,ITAMPA,ITAMPB,ITAMPC,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of C matrix transformations
*             assumes that T^A vectors on one list and 
*             T^B and T^C vectors members of a second list
*
*             ICTRAN - list of C matrix transformations
*             ICDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for ICDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*             ITAMPC - index of amplitude vector C
*
*             ITRAN - index in ICTRAN list
*             IVEC  - second index in ICDOTS list
*
*    Written by Christof Haettig, june 1997.
*
*=====================================================================*
      IMPLICIT NONE  

#include <priunit.h>
      INTEGER MXVEC, MXTRAN
      INTEGER ICTRAN(4,MXTRAN)
      INTEGER ICDOTS(MXVEC,MXTRAN)

      LOGICAL LFND
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LCTST, LCEND
      INTEGER IB, IA, IC
      LCTST(ITRAN,IA,IB,IC) =  ICTRAN(1,ITRAN).EQ.IA .AND. (
     &      ( ICTRAN(2,ITRAN).EQ.IB .AND. ICTRAN(3,ITRAN).EQ.IC ) .OR.
     &      ( ICTRAN(2,ITRAN).EQ.IC .AND. ICTRAN(3,ITRAN).EQ.IB ) )
      LCEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (ICTRAN(1,ITRAN)+ICTRAN(2,ITRAN)+ICTRAN(3,ITRAN)).LE.0

*---------------------------------------------------------------------*
* set up list of B matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFND = LCTST(ITRAN,ITAMPA,ITAMPB,ITAMPC)

      DO WHILE ( .NOT. (LFND.OR.LCEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFND = LCTST(ITRAN,ITAMPA,ITAMPB,ITAMPC)
      END DO

      IF (.NOT.LFND) THEN
        ICTRAN(1,ITRAN) = ITAMPA
        ICTRAN(2,ITRAN) = ITAMPB
        ICTRAN(3,ITRAN) = ITAMPC
        ICTRAN(4,ITRAN) = 0
        ITAMP = IZETAV
      ELSE 
        ITAMP = IZETAV
      END IF

      IVEC = 1
      DO WHILE (ICDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           ICDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO
     
      IF (.FALSE. .AND. ICDOTS(IVEC,ITRAN).NE.ITAMP) THEN
        WRITE (LUPRI,'(A,5I5)') 'CC_SETC211> add ',
     &        ITAMPA,ITAMPB,ITAMPC,IZETAV
      END IF

      ICDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC:',ITAMPA,ITAMPB,ITAMPC
        IDX = 1
        DO WHILE ( .NOT. LCEND(IDX) )
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') 'CC_SETC211>',
     &       (ICTRAN(I,IDX),I=1,4),(ICDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETC211.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*             END OF SUBROUTINE CC_SETC211                            *
*---------------------------------------------------------------------*
