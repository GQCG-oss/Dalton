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
c /* deck CC_XOPA */
*=====================================================================*
       SUBROUTINE CC_XOPA(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of first-order transition properties
*             (transition moments and oscillator strengths)
*             for transitions between two excited states with the 
*             Coupled Cluster models
*
*                        CCS, CC2, CCSD, CC3
*
*             and partially with SCF and CIS
*
*     Written by Christof Haettig winter 2002/2003.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>
#include <ccxopainf.h>
#include <ccsdinp.h>
#include <dummy.h>
#include <second.h>
#include <ccexcinf.h>
#include <ccorb.h>

* local parameters:
      CHARACTER*(16) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_XOPA> ')

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0, TIM1, TIMF, TIMXE1, TIMXE2
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0, TIM1, TIMF, TIMXE1, TIMXE2
#endif

      LOGICAL LADD
      INTEGER NBOPA, MXFTRAN, MXATRAN, MXXTRAN, MXFVEC, MXAVEC, MXXVEC,
     &        NFTRAN, NXE1TRAN, NXE2TRAN, NSTATES,
     &        KRESULT, KFTRAN, KFDOTS, KFCONS, KEND0, LEND0,
     &        KE1TRAN, KE1DOTS, KE1CONS,
     &        KX2TRAN, KX2DOTS, KX2CONS,
     &        IOPT, IORDER, ISYM

* external functions: none

*---------------------------------------------------------------------*
* print header for second-order property section:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '******************************',
     & '*                                   ',
     &                               '                             *',
     & '*<<<<<<    OUTPUT FROM COUPLED CLUST',
     &                               'ER LINEAR RESPONSE    >>>>>>>*',
     & '*<<<<<<  CALCULATION OF ONE-PHOTON A',
     &                               'ABSORPTION STRENGHTS  >>>>>>>*',
     & '*<<<<<<     FOR EXCITED TO EXCITED S',
     &                               'TATE TRANSITIONS      >>>>>>>*',
     & '*                                   ',
     &                               '                             *',
     & '************************************',
     &                               '******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CC_XOPA called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_XOPA Workspace:',LWORK
  
      TIM0  = SECOND()

*---------------------------------------------------------------------*
* allocate & initialize work space for property contributions:
*---------------------------------------------------------------------*
      ! maximum number of transition moments to compute
      NBOPA   = 2 * NQR2OP * NXQR2ST

      ! number of excited states
      NSTATES = 0
      DO ISYM = 1, NSYM
        NSTATES = NSTATES + NCCEXCI(ISYM,1)
      END DO

      ! maximum number of transformations or vector calculations
      ! NSTATES * NQR2OP   LE x Eta{X} transformations
      ! NQR2OP             Xi{X} vectors
      ! 2*NXQR2ST          LE x B x RE transformations   
      MXATRAN = NSTATES * NQR2OP
      MXXTRAN = NQR2OP
      MXFTRAN = 2*NXQR2ST

      ! maximum number of vectors to dot on 
      ! NSTATES    RE vectors dotted on a LE x Eta{X} transformation
      ! 2*NXQR2ST  N2 vectors dotted on a Xi{X} vector
      ! NQR2OP     R1 vectors dotted on a LE x B x RE transformation
      MXAVEC   = NSTATES
      MXXVEC   = 2*NXQR2ST
      MXFVEC   = NQR2OP


      KRESULT  = 1
      KEND0    = KRESULT  + NBOPA
               
      KFTRAN   = KEND0
      KFDOTS   = KFTRAN   + MXFTRAN * MXDIM_FTRAN
      KFCONS   = KFDOTS   + MXFVEC  * MXFTRAN
      KEND0    = KFCONS   + MXFVEC  * MXFTRAN

      KE1TRAN  = KEND0
      KE1DOTS  = KE1TRAN  + MXATRAN * MXDIM_XEVEC
      KE1CONS  = KE1DOTS  + MXAVEC  * MXATRAN
      KEND0    = KE1CONS  + MXAVEC  * MXATRAN

      KX2TRAN  = KEND0
      KX2DOTS  = KX2TRAN  + MXXTRAN * MXDIM_XEVEC
      KX2CONS  = KX2DOTS  + MXXVEC  * MXXTRAN
      KEND0    = KX2CONS  + MXXVEC  * MXXTRAN

      LEND0 = LWORK - KEND0
      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_XOPA. (1)')
      END IF

      CALL DZERO(WORK(KRESULT),NBOPA)

*---------------------------------------------------------------------*
* set up lists for F transformations, ETA{O} and Xi{O} vectors:
*---------------------------------------------------------------------*
      LADD = .FALSE.

      CALL CCXOPA_SETUP(WORK(KFTRAN),WORK(KFDOTS),WORK(KFCONS),
     &                  NFTRAN,MXFTRAN,MXFVEC,
     &                  WORK(KE1TRAN),WORK(KE1DOTS),WORK(KE1CONS),
     &                  NXE1TRAN,MXATRAN,MXAVEC,
     &                  WORK(KX2TRAN),WORK(KX2DOTS),WORK(KX2CONS),
     &                  NXE2TRAN,MXXTRAN,MXXVEC,
     &                  WORK(KRESULT),NBOPA,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KFCONS),MXFVEC*NFTRAN)

      IOPT = 5
      CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'LE ','RE ',IOPT,'R1 ',
     &                WORK(KFDOTS),WORK(KFCONS),MXFVEC,
     &                WORK(KEND0), LEND0)

      TIMF = SECOND() - TIM1

      IF (NFTRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NFTRAN,' F matrix transformations:',TIMF
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate LE x A{O} x RE contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KE1CONS),MXAVEC*NXE1TRAN)

      IOPT   = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KE1TRAN), NXE1TRAN, IOPT, IORDER, 'LE ',
     &               '---',DUMMY,DUMMY,
     &               'RE ',WORK(KE1DOTS),WORK(KE1CONS),
     &               .FALSE.,MXAVEC, WORK(KEND0), LEND0 )

      TIMXE1 = SECOND() - TIM1
      IF (NXE1TRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NXE1TRAN,' A{X} matrix transformations:',
     & TIMXE1
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate N2 x Xksi{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KX2CONS),MXXVEC*NXE2TRAN)

      IOPT   = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KX2TRAN), NXE2TRAN, IOPT, IORDER, '---',
     &               'N2 ',WORK(KX2DOTS),WORK(KX2CONS),
     &               '---',IDUMMY,DUMMY,
     &               .FALSE.,MXXVEC, WORK(KEND0), LEND0 )

      TIMXE2 = SECOND() - TIM1
      IF (NXE2TRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NXE2TRAN,' O1/X1 vector calculation:',TIMXE2
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* collect contributions and sum them up to the final results:
*---------------------------------------------------------------------*
      LADD = .TRUE.

      CALL CCXOPA_SETUP(WORK(KFTRAN),WORK(KFDOTS),WORK(KFCONS),
     &                  NFTRAN,MXFTRAN,MXFVEC,
     &                  WORK(KE1TRAN),WORK(KE1DOTS),WORK(KE1CONS),
     &                  NXE1TRAN,MXATRAN,MXAVEC,
     &                  WORK(KX2TRAN),WORK(KX2DOTS),WORK(KX2CONS),
     &                  NXE2TRAN,MXXTRAN,MXXVEC,
     &                  WORK(KRESULT),NBOPA,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBOPA,' quadratic response func.:', SECOND() - TIM0

*---------------------------------------------------------------------*
* print one-photon absorption properties and return:
*---------------------------------------------------------------------*
      CALL  CCOPAPRT(WORK(KRESULT),.TRUE.,NQR2OP,NXQR2ST)

      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_XOPA                              *
*=====================================================================*
c /* deck ccxopa_setup */
*=====================================================================*
      SUBROUTINE CCXOPA_SETUP(IFTRAN,  IFDOTS,  FCONS,  
     &                        NFTRAN,  MXFTRAN, MXFVEC,
     &                        IEATRAN, IEADOTS, EACONS, 
     &                        NXE1TRAN,MXATRAN, MXAVEC,
     &                        IXE2TRAN,IX2DOTS, X2CONS, 
     &                        NXE2TRAN,MXXTRAN, MXXVEC,
     &                        RESULT,  MXOPA,   LADD,   WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC first-order transition moments
*         - list of B matrix transformations with eigenvectors
*         - list of A{X} matrix transformations with eigenvectors
*         - list of XKSI vector contractions with Nbar multipliers
*
*     Written by Christof Haettig, Oct 2003 
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>
#include <ccxopainf.h>
#include <ccroper.h>
#include <ccexci.h>
#include <ccsdinp.h>
#include <ccorb.h>

* local parameters:
      CHARACTER*(22) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCXOPA_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD
      INTEGER MXOPA,MXFTRAN,MXFVEC,MXATRAN,MXAVEC,MXXTRAN,MXXVEC

      INTEGER IFTRAN(MXDIM_FTRAN,MXFTRAN)
      INTEGER IFDOTS(MXFVEC,MXFTRAN)
      INTEGER IEATRAN(MXDIM_XEVEC,MXATRAN)
      INTEGER IEADOTS(MXAVEC,MXATRAN)
      INTEGER IXE2TRAN(MXDIM_XEVEC,MXXTRAN)
      INTEGER IX2DOTS(MXXVEC,MXXTRAN)

      INTEGER NFTRAN, NXE1TRAN, NXE2TRAN, LWORK

#if defined (SYS_CRAY)
      REAL RESULT(MXOPA)
      REAL FCONS(MXFVEC,MXFTRAN)
      REAL EACONS(MXAVEC,MXATRAN)
      REAL X2CONS(MXXVEC,MXXTRAN)
      REAL WORK(LWORK)
      REAL ZERO, SIGN, EIGVI, EIGVF
      REAL WIAF, WXINIF, WIBF
#else
      DOUBLE PRECISION RESULT(MXOPA)
      DOUBLE PRECISION FCONS(MXFVEC,MXFTRAN)
      DOUBLE PRECISION EACONS(MXAVEC,MXATRAN)
      DOUBLE PRECISION X2CONS(MXXVEC,MXXTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, SIGN, EIGVI, EIGVF
      DOUBLE PRECISION WIAF, WXINIF, WIBF
#endif
      PARAMETER (ZERO = 0.0D0)

      CHARACTER LABEL*(8)
      LOGICAL LORX, LPDBS
      INTEGER ITRAN, I, IRSD, IRSDX, ISTATEI, ISTATEF, ISYMI, ISYMF,
     &        ISTISY, ISTFSY, IOP, IOPER, ISYMO, ISYME, ITURN,
     &        IKAP, MXEAVEC, MXE2VEC, IN2VEC, IR1VEC, MFVEC, 
     &        ITMIF, IVEC, NBOPA, IDUM

* external functions:
      INTEGER IR1TAMP
      INTEGER IN2AMP

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      DO ITRAN = 1, MXATRAN
       IEATRAN(1,ITRAN)  = 0
       IEATRAN(2,ITRAN)  = 0
       IEATRAN(3,ITRAN)  = -1
       IEATRAN(4,ITRAN)  = -1
       IEATRAN(5,ITRAN)  = 0
       DO IVEC = 1, MXAVEC
        IEADOTS(IVEC,ITRAN) = 0
       END DO
      END DO

      DO ITRAN = 1, MXXTRAN
       IXE2TRAN(1,ITRAN)  = 0
       IXE2TRAN(2,ITRAN)  = 0
       IXE2TRAN(3,ITRAN)  = -1
       IXE2TRAN(4,ITRAN)  = -1
       IXE2TRAN(5,ITRAN)  = 0
       DO IVEC = 1, MXXVEC
        IX2DOTS(IVEC,ITRAN) = 0
       END DO
      END DO

      DO ITRAN = 1, MXFTRAN
       DO I = 1, 3
        IFTRAN(I,ITRAN)  = 0
       END DO
       DO IVEC = 1, MXFVEC
        IFDOTS(IVEC,ITRAN)  = 0
       END DO
      END DO

      NFTRAN   = 0
      NXE1TRAN = 0
      NXE2TRAN = 0

      NBOPA   = 0
      MFVEC   = 0
      MXE2VEC = 0
      MXEAVEC = 0
 

*---------------------------------------------------------------------*
* start loop over all requested transition moments:
*---------------------------------------------------------------------*
      DO IRSDX  = 1, 2*NXQR2ST
       ITURN = 1 + (IRSDX-1)/NXQR2ST
       IRSD  = IRSDX - (ITURN-1)*NXQR2ST

       IF (ITURN.EQ.1) THEN
         ISTATEI = IQR2ST(IRSD,1)
         ISTATEF = IQR2ST(IRSD,2)
       ELSE IF (ITURN.EQ.2) THEN
         ! switch state indeces (and thereby also the sign of the freqs)
         ! to get the conjugated transition moments
         ISTATEI = IQR2ST(IRSD,2)
         ISTATEF = IQR2ST(IRSD,1)
       ELSE
         CALL QUIT('Error in CCXOPA_SETUP')
       END IF

       ISYMI   = ISYEXC(ISTATEI)
       ISYMF   = ISYEXC(ISTATEF)
       ISYME   = MULD2H(ISYMI,ISYMF)
       ISTISY  = ISTATEI - ISYOFE(ISYMI)
       ISTFSY  = ISTATEF - ISYOFE(ISYMF)
       EIGVI   = EIGVAL(ISTATEI)
       EIGVF   = EIGVAL(ISTATEF)

       IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'CCXOPA_SETUP:'
         WRITE(LUPRI,*) 'ITURN,IRSD:',ITURN,IRSD
         WRITE(LUPRI,*) 'ISTATEI,ISTATEF:',ISTATEI,ISTATEF
         WRITE(LUPRI,*) 'ISYMI,ISYMF:',ISYMI,ISYMF
         WRITE(LUPRI,*) 'ISTISY,ISTFSY:',ISTISY,ISTFSY
         WRITE(LUPRI,*) 'EIGVI,EIGVF:',EIGVI,EIGVF
       END IF

       DO IOP = 1, NQR2OP
        IOPER = IQR2OP(IOP)
        LORX  = .FALSE.
        ISYMO = ISYOPR(IOPER)
        LABEL = LBLOPR(IOPER)
        LPDBS = LPDBSOP(IOPER)
        IKAP  = 0

        IF (LPDBS) CALL QUIT('perturbation-dependent basis sets not '//
     &              'implemented in CCXOPA_SETUP.')


        IF (ISYMO.EQ.ISYME) THEN 

          NBOPA = NBOPA + 1

          IF (NBOPA.GT.MXOPA) THEN
             CALL QUIT('NBOPA out of range in CCXOPA_SETUP.')
          END IF

*---------------------------------------------------------------------*
*         in all cases we need LE x A{X} x RE 
*---------------------------------------------------------------------*
          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXATRAN,MXAVEC,
     &                  ISTATEI,IOPER,IKAP,0,0,0,ISTATEF,ITRAN,IVEC)
          NXE1TRAN = MAX(NXE1TRAN,ITRAN)
          MXEAVEC  = MAX(MXEAVEC, IVEC)
          WIAF     = EACONS(IVEC,ITRAN)

*---------------------------------------------------------------------*
*         add N2 * Xksi{X} or LE * B * RE * R1, depending on QR22N1
*---------------------------------------------------------------------*
          WXINIF = ZERO
          WIBF   = ZERO

          IF (.NOT.CIS) THEN
            IF (QR22N1) THEN
              IN2VEC=IN2AMP(ISTATEI,-EIGVI,ISYMI,ISTATEF,+EIGVF,ISYMF)
              CALL CC_SETXE('Xi ',IXE2TRAN,IX2DOTS,MXXTRAN,MXXVEC,
     &                      0,IOPER,IKAP,0,0,0,IN2VEC,ITRAN,IVEC)
              NXE2TRAN = MAX(NXE2TRAN,ITRAN)
              MXE2VEC  = MAX(MXE2VEC, IVEC)
              WXINIF   = X2CONS(IVEC,ITRAN)
            ELSE
              IR1VEC = IR1TAMP(LABEL,LORX,EIGVI-EIGVF,IDUM)
              CALL CC_SETF12(IFTRAN,IFDOTS,MXFTRAN,MXFVEC,
     &                       ISTATEI,ISTATEF,IR1VEC,ITRAN,IVEC)
              NFTRAN = MAX(NFTRAN,ITRAN)
              MFVEC  = MAX(MFVEC, IVEC)
              WIBF   = FCONS(IVEC,ITRAN)
            END IF
          END IF

*---------------------------------------------------------------------*
*          add contributions together:
*---------------------------------------------------------------------*
           IF (LADD) THEN

              ITMIF = (NQR2OP*(IRSD-1) + IOP-1)*2 + ITURN

              RESULT(ITMIF) = WIAF + WXINIF + WIBF

              IF (LOCDBG) THEN
                 WRITE (LUPRI,*) 'ISTATEI, EIGVI:',ISTATEI,EIGVI
                 WRITE (LUPRI,*) 'ISTATEF, EIGVF:',ISTATEF,EIGVF
                 WRITE (LUPRI,*) 'OPERATOR:',LABEL
                 WRITE (LUPRI,*) 'IDX = ',ITMIF
                 WRITE (LUPRI,*) 'L^i A{X} x R^f :',WIAF
                 WRITE (LUPRI,*) 'N^if x Xksi{X}:',WXINIF
                 WRITE (LUPRI,*) 'L^i x B x R^f x R^X:',WIBF
                 WRITE (LUPRI,*) 'Total result:',RESULT(ITMIF)
              END IF

           END IF

*---------------------------------------------------------------------*
*       end loop over transition moments
*---------------------------------------------------------------------*

        END IF
       END DO
      END DO

      IF      (MFVEC.GT.MXFVEC) THEN
         CALL QUIT('MFVEC has been out of bounds in CCXOPA_SETUP.')
      ELSE IF (MXEAVEC.GT.MXAVEC) THEN
         CALL QUIT('MXEAVEC has been out of bounds in CCXOPA_SETUP.')
      ELSE IF (MXE2VEC.GT.MXXVEC) THEN
         CALL QUIT('MXE2VEC has been out of bounds in CCXOPA_SETUP.')
      ELSE IF (NFTRAN.GT.MXFTRAN) THEN
         CALL QUIT('NFTRAN has been out of bounds in CCXOPA_SETUP.')
      ELSE IF (NXE1TRAN.GT.MXATRAN) THEN
         CALL QUIT('NXE1TRAN has been out of bounds in CCXOPA_SETUP.')
      ELSE IF (NXE2TRAN.GT.MXXTRAN) THEN
         CALL QUIT('NXE2TRAN has been out of bounds in CCXOPA_SETUP.')
      END IF

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      IF ((.NOT.LADD) .OR. LOCDBG) THEN
       WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NBOPA,
     &      ' transition moments'
       WRITE(LUPRI,'((8X,A,I3,A))') 
     & ' - ',NFTRAN,  ' F matrix transformations with RE vectors',
     & ' - ',NXE1TRAN,' A{X} matrix transformations with LE vectors',
     & ' - ',NXE2TRAN,' extra XKSI vector calculations '
       WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'
      END IF

      IF (LOCDBG) THEN

         ! F matrix transformations:
         WRITE(LUPRI,*)'List of F matrix transformations:'
         DO ITRAN = 1, NFTRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IFTRAN(I,ITRAN),I=1,2),(IFDOTS(I,ITRAN),I=1,MFVEC)
         END DO
         WRITE(LUPRI,*)

         ! LE x A{X} vector calculations:
         WRITE(LUPRI,*) 'List of A{O} matrix transformations:'
         DO ITRAN = 1, NXE1TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IEATRAN(I,ITRAN),I=1,5),(IEADOTS(I,ITRAN),I=1,MXEAVEC)
         END DO
         WRITE(LUPRI,*)

         ! extra Xi{O} vector calculations:
         WRITE(LUPRI,*) 'List of extra Xi{O} vector calculations:'
         DO ITRAN = 1, NXE2TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IXE2TRAN(I,ITRAN),I=1,5),(IX2DOTS(I,ITRAN),I=1,MXE2VEC)
         END DO
         WRITE(LUPRI,*)

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCXOPA_SETUP                         *
*---------------------------------------------------------------------*
