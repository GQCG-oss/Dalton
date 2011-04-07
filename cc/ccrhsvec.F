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
c /* deck CCRHSVEC */
*=====================================================================*
      SUBROUTINE CCRHSVEC(TYPE,LABEL,ISYMS,ISTAT,EIGV,ISYMO,
     &                    FREQS,LORX,ICAU,NVEC,MAXVEC,IOFFV,
     &                    WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: calculate right hand side vectors for higher-order
*             coupled cluster amplitude response equations,
*             left and right excited state response equations
*
*             if called for ORDER=n, the solutions for ORDER=n-1
*             must be available on file
*
*     implemented:  T:   ORDER = 1, 2, 3, 4 
*                   RE:  ORDER = 1, 2
*                   CR:  ORDER = 1, 2
*
*    Written by Christof Haettig maj 1997, extension to RE july '97
*                                          extension to CR march '98
*                                          extension to O1 jan '99
*    orb.-relax. or derivatives by Christof Haettig, Aug '99.        
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cclists.h>
#include <dummy.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky

* local parameters:
      CHARACTER*(18) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCRHSVEC> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )



      CHARACTER TYPE*(*), LISTR*3

      INTEGER NVEC, MAXVEC, IOFFV, LWORK
      INTEGER ISYMS(MAXVEC,*), ISYMO(MAXVEC,*)
      INTEGER ISTAT(MAXVEC,*), ICAU(MAXVEC,*)
      LOGICAL LORX(MAXVEC,*)

      CHARACTER*8 LABEL(MAXVEC,*)

#if defined (SYS_CRAY)
      REAL FREQS(MAXVEC,*), EIGV(MAXVEC,*)
      REAL WORK(LWORK)
      REAL ZERO, RDUM
      REAL XNORM, DDOT
#else
      DOUBLE PRECISION FREQS(MAXVEC,*), EIGV(MAXVEC,*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, RDUM
      DOUBLE PRECISION XNORM, DDOT
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER MODEL*(10), MODELW*(10)
      CHARACTER APROXR12*(3)
      LOGICAL NEW_RHS
      INTEGER IOPT, ISYM, IVEC, MPERM, NSTAT, ORDER, IDUM, IOPTE
      INTEGER MXTRAN,MXDTRAN,MXCTRAN,MXBTRAN,MXBATRAN,MXAATRAN,MXXETRAN
      INTEGER MXCATRAN
      INTEGER KDTRAN,KCTRAN,KB1TRAN,KB2TRAN,KBA1TRAN,KAA1TRAN,KXETRAN
      INTEGER NDTRAN,NCTRAN,NB1TRAN,NB2TRAN,NBA1TRAN,NAA1TRAN,NXETRAN
      INTEGER KCATRAN, KBA2TRAN, KAA2TRAN
      INTEGER NCATRAN, NBA2TRAN, NAA2TRAN
      INTEGER KEND0, LEND0, LMAX1, LMAX2, KRHS1, KRHS2, KEND1, LEND1
      INTEGER KLHS1, KLHS2, KEND2, LEND2, IDXVEC 
      INTEGER KRHSR12, LMAXR12, IOPTR12, MODLEN

* external functions:
      INTEGER ILSTSYM, ILRCAMP


*---------------------------------------------------------------------*
* check number of required rhs vectors, if zero return immediatly:
*---------------------------------------------------------------------*
      IF (NVEC.EQ.0) RETURN

*---------------------------------------------------------------------*
* print header for rhs vector section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '------------------------------------',
     &                               '-------------------------------',
     & '|          OUTPUT FROM AMPLITUDE RHS',    
     &                               ' VECTOR SECTION               |',
     & '------------------------------------',
     &                               '-------------------------------' 
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CCRHSVEC called for unknown Coupled Cluster model.')
      END IF

      NEW_RHS = .FALSE.
  
      IF (TYPE(1:3).EQ.'O1 ') THEN
        ORDER = 1
        NSTAT = 0
        MPERM = 1
      ELSE IF (TYPE(1:2).EQ.'O2') THEN
        ORDER = 2
        NSTAT = 0
        MPERM = 2
        ! compute complete O2 vector in B matrix module
        NEW_RHS = .TRUE.
      ELSE IF (TYPE(1:2).EQ.'O3') THEN
        ORDER = 3
        NSTAT = 0
        MPERM = 3
      ELSE IF (TYPE(1:2).EQ.'O4') THEN
        ORDER = 4
        NSTAT = 0
        MPERM = 12
        WRITE (LUPRI,*) 'warning: rhs vectors ',TYPE(1:2),
     &       ' not tested!!!.'
      ELSE IF (TYPE(1:3).EQ.'EO1') THEN
        ORDER = 1
        NSTAT = 1
        MPERM = 1
        ! compute complete EO1 vector in B matrix module
        NEW_RHS = .TRUE.
      ELSE IF (TYPE(1:3).EQ.'EO2') THEN
        ORDER = 2
        NSTAT = 1
        MPERM = 2
        WRITE (LUPRI,*) 'warning: rhs vectors ',TYPE(1:3),
     &       ' not tested!!!.'
      ELSE IF (TYPE(1:3).EQ.'CO1') THEN
        ORDER = 1
        NSTAT = 0
        MPERM = 1
      ELSE IF (TYPE(1:3).EQ.'CO2') THEN
        ORDER = 2
        NSTAT = 0
        MPERM = 2
      ELSE
        WRITE (LUPRI,*) 'rhs vectors ',TYPE(1:2),' not implemented.'
        CALL QUIT('required rhs vectors not implemented.')
      END IF

* Cholesky check: only CC2 O1 has been implemented
      IF (CHOINT .AND. CC2) THEN
         IF (TYPE(1:2).NE.'O1') THEN
            WRITE (LUPRI,*)
     &      'rhs vectors ',TYPE(1:2),' not implemented for Cholesky.'
            CALL QUIT('required rhs vectors not implemented.')
         ENDIF
      ENDIF

* print some debug/info output
      IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
        WRITE(LUPRI,*) 'CCRHSVEC Workspace:',LWORK
      END IF
  
*---------------------------------------------------------------------*
* allocate & initialize work space for lists
*---------------------------------------------------------------------*

      MXTRAN   = MPERM * NVEC

      MXDTRAN  = MXDIM_DTRAN  * MXTRAN
      MXCTRAN  = MXDIM_CTRAN  * MXTRAN
      MXBTRAN  = MXDIM_BTRAN  * MXTRAN
      MXBATRAN = MXDIM_BATRAN * MXTRAN
      MXCATRAN = MXDIM_CATRAN * MXTRAN
      MXAATRAN = MXDIM_AATRAN * MXTRAN
      MXXETRAN = MXDIM_XEVEC  * MXTRAN

      KDTRAN   = 1
      KCTRAN   = KDTRAN   + MXDTRAN
      KB1TRAN  = KCTRAN   + MXCTRAN
      KB2TRAN  = KB1TRAN  + MXBTRAN
      KCATRAN  = KB2TRAN  + MXBTRAN
      KBA1TRAN = KCATRAN  + MXCATRAN
      KBA2TRAN = KBA1TRAN + MXBATRAN
      KAA1TRAN = KBA2TRAN + MXBATRAN
      KAA2TRAN = KAA1TRAN + MXAATRAN
      KXETRAN  = KAA2TRAN + MXAATRAN
      KEND0    = KXETRAN  + MXXETRAN
      LEND0    = LWORK    - KEND0


      IF (LEND0 .LT. 0 ) THEN
        WRITE (LUPRI,*) 'Insufficient work space in CCRHSVEC.'
        WRITE (LUPRI,*) 'KEND0, LEND0, LWORK:',KEND0,LEND0,LWORK
        WRITE (LUPRI,*) 'MXTRAN:',MXTRAN
        CALL QUIT('Insufficient work space in CCRHSVEC.')
      END IF

*---------------------------------------------------------------------*
* set up lists for D, C, B, B{O} and A{O} transformations:
*---------------------------------------------------------------------*
      CALL CC_RHS_SETUP(TYPE,NSTAT,ORDER,LABEL,ISTAT,EIGV,ISYMO,FREQS,
     &                  LORX, ICAU, NVEC, MAXVEC, IOFFV, MXTRAN, 
     &                  NEW_RHS,
     &                  WORK(KDTRAN),  NDTRAN,
     &                  WORK(KCTRAN),  NCTRAN,
     &                  WORK(KB1TRAN), NB1TRAN,
     &                  WORK(KB2TRAN), NB2TRAN,
     &                  WORK(KCATRAN), NCATRAN,
     &                  WORK(KBA1TRAN),NBA1TRAN,
     &                  WORK(KBA2TRAN),NBA2TRAN,
     &                  WORK(KAA1TRAN),NAA1TRAN,
     &                  WORK(KAA2TRAN),NAA2TRAN,
     &                  WORK(KXETRAN), NXETRAN )

*---------------------------------------------------------------------*
* initialize rhs vector files:
*---------------------------------------------------------------------*
      IF( TYPE(1:3).NE.'O1 ' .AND. TYPE(1:3).NE.'X1 ' .and.
     &    TYPE(1:3).NE.'CO1'                                ) THEN

        LMAX1 = 0
        LMAX2 = 0
        LMAXR12 = 0
        DO ISYM = 1, NSYM
          LMAX1 = MAX(LMAX1,NT1AM(ISYM))
          LMAX2 = MAX(LMAX2,NT2AM(ISYM))
          IF (CCR12) LMAXR12 = MAX(LMAXR12,NTR12AM(ISYM))
        END DO

        KRHS1 = KEND0
        KRHS2 = KRHS1 + LMAX1
        KRHSR12 = KRHS2 + LMAX2
        KEND1 = KRHSR12 + LMAXR12
        LEND1 = LWORK - KEND1

        IF (LEND1 .LT. 0 ) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CCRHSVEC. (2)'
          WRITE (LUPRI,*) 'KEND1, LEND1, LWORK:',KEND1,LEND1,LWORK
          CALL QUIT('Insufficient work space in CCRHSVEC. (2)')
        END IF

        CALL DZERO(WORK(KRHS1),LMAX1)
        IF (.NOT.CCS) CALL DZERO(WORK(KRHS2),LMAX2)
        IF (CCR12) CALL DZERO(WORK(KRHSR12),LMAXR12)

        IF (CCS) THEN
           MODEL = 'CCS       '
           IOPT  = 1
        ELSE IF (CC2) THEN
           MODEL = 'CC2       '
           IOPT  = 3
        ELSE IF (CCSD) THEN
           MODEL = 'CCSD      '
           IOPT  = 3
        ELSE IF (CC3) THEN
           MODEL = 'CC3       '
           ! intialize usual and effective rhs vector
           IOPT  = 3
           IOPTE = 24
        ELSE
           CALL QUIT('Unknown coupled cluster model in CCRHSVEC.')
        END IF
        IF (CCR12) THEN
          APROXR12 = '   '
          IOPTR12 = 32
        END IF
        CALL CCSD_MODEL(MODELW,MODLEN,10,MODEL,10,APROXR12)

        DO IVEC = IOFFV+1, IOFFV+NVEC
         ISYM = ILSTSYM(TYPE,IVEC)
         CALL CC_WRRSP(TYPE,IVEC,ISYM,IOPT,MODELW,IDUMMY,
     &                 WORK(KRHS1),WORK(KRHS2),WORK(KEND1),LEND1)
         IF (CCR12) THEN
           CALL CC_WRRSP(TYPE,IVEC,ISYM,IOPTR12,MODELW,IDUMMY,
     &                   IDUMMY,WORK(KRHSR12),WORK(KEND1),LEND1)
         END IF
         IF (CCSDT) THEN
           CALL CC_WRRSP(TYPE,IVEC,ISYM,IOPTE,MODELW,IDUMMY,
     &                   WORK(KRHS1),WORK(KRHS2),WORK(KEND1),LEND1)
           CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPTE,MODEL,
     &                   WORK(KRHS1),WORK(KRHS2))
         END IF
        END DO
      
      END IF

*---------------------------------------------------------------------*
* calculate D matrix contributions:
*---------------------------------------------------------------------*
      IF (TYPE(1:2).EQ.'O4') THEN
        IOPT = 4
        CALL CC_DMAT(WORK(KDTRAN),NDTRAN,
     &               'R1 ','R1 ','R1 ','R1 ',IOPT,TYPE,
     &               IDUM, RDUM, 0, WORK(KEND0), LEND0 )
      END IF

      IF (LOCDBG .AND. TYPE(1:3).NE.'O1 '.AND. TYPE(1:3).NE.'X1 '
     &           .AND. TYPE(1:3).NE.'CO1'                        ) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'NORM^2 of RHS vectors after D matrix terms:'
        DO IVEC = IOFFV+1, IOFFV+NVEC
          IOPT = 3
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          XNORM = DDOT(NT1AM(ISYM),WORK(KRHS1),1,WORK(KRHS1),1)
          IF (.NOT. CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KRHS2),1,WORK(KRHS2),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF


*---------------------------------------------------------------------*
* calculate C matrix contributions:
*---------------------------------------------------------------------*
      IF      (TYPE(1:2).EQ.'O4') THEN
        IOPT = 4
        CALL CC_CMAT(WORK(KCTRAN),NCTRAN,'R2 ','R1 ','R1 ',IOPT,TYPE,
     &               IDUM, RDUM, 0, WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
        IOPT = 4
        CALL CC_CMAT(WORK(KCTRAN),NCTRAN,'R1 ','R1 ','R1 ',IOPT,TYPE,
     &               IDUM, RDUM, 0, WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:3).EQ.'EO2' ) THEN
        IOPT = 4
        CALL CC_CMAT(WORK(KCTRAN),NCTRAN,'R1 ','R1 ','RE ',IOPT,TYPE,
     &               IDUM, RDUM, 0, WORK(KEND0), LEND0 )
      END IF

      IF (LOCDBG .AND. TYPE(1:3).NE.'O1 '.AND. TYPE(1:3).NE.'X1 '
     &           .AND. TYPE(1:3).NE.'CO1'                        ) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'NORM^2 of RHS vectors after C matrix terms:'
        DO IVEC = IOFFV+1, IOFFV+NVEC
          IOPT = 3
          IF (CC3) IOPT = 24
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          XNORM = DDOT(NT1AM(ISYM),WORK(KRHS1),1,WORK(KRHS1),1)
          IF (.NOT. CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KRHS2),1,WORK(KRHS2),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate B matrix contributions:
*---------------------------------------------------------------------*
      IF      ( TYPE(1:2).EQ.'O4' ) THEN
        IOPT = 4
        CALL CC_BMAT(WORK(KB1TRAN), NB1TRAN,'R3 ','R1 ',IOPT,TYPE,
     &               IDUM, RDUM, 0, .FALSE.,WORK(KEND0), LEND0  )
        IOPT = 4
        CALL CC_BMAT(WORK(KB2TRAN), NB2TRAN,'R2 ','R2 ',IOPT,TYPE,
     &               IDUM, RDUM, 0, .FALSE.,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
        IOPT = 4
        CALL CC_BMAT(WORK(KB1TRAN), NB1TRAN,'R2 ','R1 ',IOPT,TYPE,
     &               IDUM, RDUM, 0, .FALSE.,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:2).EQ.'O2' ) THEN
        IOPT = 4
        CALL CC_BMAT(WORK(KB1TRAN), NB1TRAN,'R1 ','R1 ',IOPT,TYPE,
     &               IDUM, RDUM, 0, NEW_RHS,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:3).EQ.'EO2') THEN
        IOPT = 4
        CALL CC_BMAT(WORK(KB1TRAN), NB1TRAN,'R2 ','RE ',IOPT,TYPE,
     &               IDUM, RDUM, 0, .FALSE.,WORK(KEND0), LEND0 ) 
        IOPT = 4
        CALL CC_BMAT(WORK(KB2TRAN), NB2TRAN,'R1 ','ER1',IOPT,TYPE,
     &               IDUM, RDUM, 0, .FALSE.,WORK(KEND0), LEND0 ) 
      ELSE IF ( TYPE(1:3).EQ.'EO1') THEN
        IOPT = 4
        CALL CC_BMAT(WORK(KB1TRAN), NB1TRAN,'R1 ','RE ',IOPT,TYPE,
     &               IDUM, RDUM, 0, NEW_RHS,WORK(KEND0), LEND0 ) 
      ELSE IF ( TYPE(1:3).EQ.'CO2' ) THEN
        IOPT = 4
        CALL CC_BMATRIX(WORK(KB1TRAN), NB1TRAN,'RC ','RC ',IOPT,TYPE,
     &                  IDUM, RDUM, 0, .FALSE.,WORK(KEND0), LEND0 ) 
      END IF

      IF (LOCDBG .AND. TYPE(1:3).NE.'O1 '.AND. TYPE(1:3).NE.'X1 '
     &           .AND. TYPE(1:3).NE.'CO1'                        ) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'NORM^2 of RHS vectors after B matrix terms:'
        DO IVEC = IOFFV+1, IOFFV+NVEC
          IOPT = 3
          IF (CC3) IOPT = 24
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          XNORM = DDOT(NT1AM(ISYM),WORK(KRHS1),1,WORK(KRHS1),1)
          IF (.NOT. CCS) 
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KRHS2),1,WORK(KRHS2),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate C{O} matrix contributions:
*---------------------------------------------------------------------*
      IF      ( TYPE(1:2).EQ.'O4' .AND. NCATRAN.NE.0) THEN
        IOPT = 4
c       CALL CC_CAMAT(WORK(KCATRAN),NCATRAN,'o1 ','R1 ','R1 ','R1 ',
c    &                IOPT, TYPE, IDUM, RDUM, 0, WORK(KEND0), LEND0 )
        CALL QUIT('cc_camat routine not yet implememted.')
      END IF

*---------------------------------------------------------------------*
* calculate B{O} matrix contributions:
*---------------------------------------------------------------------*
      IF      ( TYPE(1:2).EQ.'O4' ) THEN
        IOPT = 4
        CALL CC_BAMAT(WORK(KBA1TRAN),NBA1TRAN,'o1 ','R2 ','R1 ',IOPT,
     &                TYPE, IDUM, RDUM, 0,WORK(KEND0), LEND0 )
        IOPT = 4
        CALL CC_BAMAT(WORK(KBA2TRAN),NBA2TRAN,'o2 ','R1 ','R1 ',IOPT,
     &                TYPE, IDUM, RDUM, 0,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
        IOPT = 4
        CALL CC_BAMAT(WORK(KBA1TRAN),NBA1TRAN,'o1 ','R1 ','R1 ',IOPT,
     &                TYPE, IDUM, RDUM, 0,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:3).EQ.'EO2' ) THEN
        IOPT = 4
        CALL CC_BAMAT(WORK(KBA1TRAN),NBA1TRAN,'o1 ','R1 ','RE ',IOPT,
     &                TYPE, IDUM, RDUM, 0,WORK(KEND0), LEND0 )
      END IF

      IF (LOCDBG .AND. TYPE(1:3).NE.'O1 '.AND. TYPE(1:3).NE.'X1 '
     &           .AND. TYPE(1:3).NE.'CO1'                        ) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'NORM^2 of RHS vectors after B{O} matrix terms:'
        DO IVEC = IOFFV+1, IOFFV+NVEC
          IOPT = 3
          IF (CC3) IOPT = 24
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          XNORM = DDOT(NT1AM(ISYM),WORK(KRHS1),1,WORK(KRHS1),1)
          IF (.NOT. CCS)
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KRHS2),1,WORK(KRHS2),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate A{O} matrix contributions:
*---------------------------------------------------------------------*
      IF      ( TYPE(1:2).EQ.'O4' ) THEN
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA1TRAN),NAA1TRAN,'o1 ','R3 ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA2TRAN),NAA2TRAN,'o2 ','R2 ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA1TRAN),NAA1TRAN,'o1 ','R2 ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA2TRAN),NAA2TRAN,'o2 ','R1 ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:2).EQ.'O2' .AND. (.NOT.NEW_RHS)) THEN
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA1TRAN),NAA1TRAN,'o1 ','R1 ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:3).EQ.'EO2' ) THEN
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA1TRAN),NAA1TRAN,'o1 ','ER1',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA2TRAN),NAA2TRAN,'o2 ','RE ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:3).EQ.'EO1' .AND. (.NOT.NEW_RHS)) THEN
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA1TRAN),NAA1TRAN,'o1 ','RE ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
      ELSE IF ( TYPE(1:3).EQ.'CO2' ) THEN
        IOPT = 4
        CALL CC_AAMAT(WORK(KAA1TRAN),NAA1TRAN,'o1 ','RC ',IOPT,TYPE,
     &                IDUMMY,DUMMY,1,WORK(KEND0), LEND0 )
      END IF

      IF (LOCDBG .AND. TYPE(1:3).NE.'O1 '.AND. TYPE(1:3).NE.'X1 '
     &           .AND. TYPE(1:3).NE.'CO1'                        ) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'NORM^2 of RHS vectors after A{O} matrix terms:'
        DO IVEC = IOFFV+1, IOFFV+NVEC
          IOPT = 3
          IF (CC3) IOPT = 24
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          XNORM = DDOT(NT1AM(ISYM),WORK(KRHS1),1,WORK(KRHS1),1)
          IF (.NOT. CCS)
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KRHS2),1,WORK(KRHS2),1)
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate Xi{O} vector contributions:
*---------------------------------------------------------------------*

Cholesky
*
*   the Cholesky routine is *far* less general, hence most of the
*   input is implicit: all that's needed are the operator indices
*   in KXETRAN.

      IF (CHOINT .AND. TYPE(1:2).EQ.'O1') THEN
         CALL CC_CHOXI0(WORK(KXETRAN),NXETRAN,WORK(KEND0),LEND0)
         GOTO 1234
      END IF
Cholesky

      IF ( TYPE(1:3).EQ.'O1 '.OR. TYPE(1:3).EQ.'X1 ') THEN
        IOPT = 3
        CALL CC_XIETA(WORK(KXETRAN),NXETRAN,IOPT, ORDER, 'L0 ',
     &                'O1 ', IDUM, RDUM, 'X1 ', IDUM, RDUM,
     &                .FALSE.,0, WORK(KEND0),LEND0)
      ELSE IF ( TYPE(1:2).EQ.'O2' .OR. TYPE(1:2).EQ.'X2' ) THEN
        IOPT = 3
        CALL CC_XIETA(WORK(KXETRAN),NXETRAN,IOPT, ORDER, 'L0 ',
     &                'O2 ', IDUM, RDUM, 'X2 ', IDUM, RDUM,
     &                .FALSE.,0, WORK(KEND0),LEND0)
      ELSE IF ( TYPE(1:3).EQ.'CO1' ) THEN
        IOPT = 3
        CALL CC_XIETA(WORK(KXETRAN),NXETRAN,IOPT, ORDER, 'L0 ',
     &                'RC ', IDUM, RDUM, '---', IDUM, RDUM,
     &                .TRUE.,0, WORK(KEND0),LEND0)
      END IF

 1234 CONTINUE     ! From Cholesky

      IF (LOCDBG) THEN
        LMAX1 = 0
        LMAX2 = 0
        LMAXR12 = 0
        DO ISYM = 1, NSYM
          LMAX1 = MAX(LMAX1,NT1AM(ISYM))
          LMAX2 = MAX(LMAX2,NT2AM(ISYM))
          IF (CCR12) THEN
            LMAXR12 = MAX(LMAXR12,NTR12AM(ISYM))
          END IF
        END DO

        KRHS1 = KEND0
        KRHS2 = KRHS1 + LMAX1
        KRHSR12 = KRHS2 + LMAX2
        KEND1 = KRHSR12 + LMAXR12
        LEND1 = LWORK - KEND1

        IF (LEND1 .LT. 0 ) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CCRHSVEC. (3)'
          WRITE (LUPRI,*) 'KEND1, LEND1, LWORK:',KEND1,LEND1,LWORK
          CALL QUIT('Insufficient work space in CCRHSVEC. (3)')
        END IF

        WRITE (LUPRI,*) MSGDBG,
     &       'NORM^2 of RHS vectors after Xi{O} matrix terms:'
        LISTR = TYPE(1:3)
        IF (TYPE.EQ.'CO1') LISTR ='RC '
        DO IVEC = IOFFV+1, IOFFV+NVEC
          IOPT = 3
          IF (CC3) IOPT = 24
          ISYM = ILSTSYM(LISTR,IVEC)
          IF (TYPE.EQ.'CO1') THEN
           WRITE(LUPRI,*) 'Cauchy order:',ICAU(IVEC,1)
           IDXVEC=ILRCAMP(LABEL(IVEC,1),ICAU(IVEC,1)-1,ISYM)
          ELSE
           IDXVEC= IVEC
          END IF
          CALL CC_RDRSP(LISTR,IDXVEC,ISYM,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          XNORM = DDOT(NT1AM(ISYM),WORK(KRHS1),1,WORK(KRHS1),1)
Chol      IF (.NOT. CCS)
          IF ((.NOT. CCS) .AND. (.NOT. (CHOINT.AND.CC2)))
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KRHS2),1,WORK(KRHS2),1)
          IF (CCR12) THEN
            CALL CC_RDRSP(LISTR,IDXVEC,ISYM,IOPTR12,MODEL,
     &                  DUMMY,WORK(KRHSR12))
            XNORM = XNORM + DDOT(NTR12AM(ISYM),WORK(KRHSR12),1,
     &                           WORK(KRHSR12),1)
          END IF
          WRITE (LUPRI,*) MSGDBG, IVEC,XNORM
        END DO
      END IF

*---------------------------------------------------------------------*
* test (static) EO1 vectors by calculating the excited state FOP's
*---------------------------------------------------------------------*
      IF (LOCDBG .AND. TYPE(1:3).EQ.'EO1') THEN
        LMAX1 = 0
        LMAX2 = 0
        DO ISYM = 1, NSYM
          LMAX1 = MAX(LMAX1,NT1AM(ISYM))
          LMAX2 = MAX(LMAX2,NT2AM(ISYM))
        END DO

        KLHS1 = KEND1
        KLHS2 = KLHS1 + LMAX1
        KEND2 = KLHS2 + LMAX2
        LEND2 = LWORK - KEND2

        IF (LEND2 .LT. 0 ) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CCRHSVEC. (4)'
          WRITE (LUPRI,*) 'KEND2, LEND2, LWORK:',KEND2,LEND2,LWORK
          CALL QUIT('Insufficient work space in CCRHSVEC. (4)')
        END IF

        WRITE (LUPRI,*) MSGDBG, 'excited state first order properties:'
        DO IVEC = IOFFV+1, IOFFV+NVEC
        IF (ISYMO(IVEC,1).EQ.1 .AND. FREQS(IVEC,1).EQ.ZERO) THEN
          IOPT = 3
          ISYM = ILSTSYM(TYPE,IVEC)
          CALL CC_RDRSP(TYPE,IVEC,ISYM,IOPT,MODEL,
     &                  WORK(KRHS1),WORK(KRHS2))
          CALL CCLR_DIASCL(WORK(KRHS2),0.5d0,ISYM)
          CALL CC_RDRSP('LE',ISTAT(IVEC,1),ISYMS(IVEC,1),IOPT,MODEL,
     &                  WORK(KLHS1),WORK(KLHS2))
          XNORM = DDOT(NT1AM(ISYM),WORK(KLHS1),1,WORK(KRHS1),1)
          IF (.NOT. CCS)
     &     XNORM = XNORM+DDOT(NT2AM(ISYM),WORK(KLHS2),1,WORK(KRHS2),1)
          WRITE (LUPRI,'(A,I3,2X,F12.8,2X,A,2X,F12.8)') MSGDBG,
     &              ISTAT(IVEC,1),EIGV(IVEC,1),LABEL(IVEC,1),XNORM
        ELSE
          WRITE (LUPRI,'(A,I3,2X,F12.8,2X,A,2X,F12.8)') MSGDBG,
     &              ISTAT(IVEC,1),EIGV(IVEC,1),LABEL(IVEC,1),ZERO
        END IF
        END DO
      END IF

*---------------------------------------------------------------------*
* that's it:
*---------------------------------------------------------------------*

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CCRHSVEC                             *
*=====================================================================*
c /* deck CC_RHS_SETUP */
*=====================================================================*
      SUBROUTINE CC_RHS_SETUP(TYPE,NSTAT,ORDER,LAB,ISTAT,
     &                        EIGV,ISYMO,FREQ,LORX,ICAU,
     &                        NVEC,MAXVEC,IOFFV,MXTRAN, 
     &                        NEW_RHS,
     &                        IDTRAN,  NDTRAN,
     &                        ICTRAN,  NCTRAN,
     &                        IB1TRAN, NB1TRAN,
     &                        IB2TRAN, NB2TRAN,
     &                        ICATRAN, NCATRAN,
     &                        IBA1TRAN,NBA1TRAN,
     &                        IBA2TRAN,NBA2TRAN,
     &                        IAA1TRAN,NAA1TRAN,
     &                        IAA2TRAN,NAA2TRAN,
     &                        IXETRAN, NXETRAN )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CCRHSVEC section
*                - list of D matrix transformations 
*                - list of C matrix transformations 
*                - list of B matrix transformations 
*                - list of B{O} matrix transformations 
*                - list of A{O} matrix transformations 
*                - list of Xi{O} vector calculations 
*
*     Written by Christof Haettig, maj 1997.
*     O1, O2, O3, O4, EO1, EO2 with one perturbation including
*     orb.-relax. or derivatives by Christof Haettig, Aug '99.        
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccroper.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(22) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_RHS_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER MXORD, MXORD2, MXORD3
      PARAMETER (MXORD  = 4)
      PARAMETER (MXORD2 = MXORD *(MXORD-1)/2 )
      PARAMETER (MXORD3 = MXORD2*(MXORD-2)/3 )

      INTEGER MXTRAN, NSTAT, ORDER, MAXVEC, NVEC, IOFFV
 
      CHARACTER*(*) TYPE
    
      CHARACTER*(8) LAB(MAXVEC,*)
      INTEGER ISTAT(MAXVEC,*), ICAU(MAXVEC,*), ISYMO(MAXVEC,*)
      LOGICAL LORX(MAXVEC,*), NEW_RHS

#if defined (SYS_CRAY)
      REAL FREQ(MAXVEC,*), EIGV(MAXVEC)
#else
      DOUBLE PRECISION FREQ(MAXVEC,*), EIGV(MAXVEC)
#endif

      INTEGER IDTRAN( MXDIM_DTRAN, MXTRAN)
      INTEGER ICTRAN( MXDIM_CTRAN, MXTRAN)
      INTEGER IB1TRAN(MXDIM_BTRAN, MXTRAN)
      INTEGER IB2TRAN(MXDIM_BTRAN, MXTRAN)
      INTEGER ICATRAN(MXDIM_CATRAN,MXTRAN)
      INTEGER IBA1TRAN(MXDIM_BATRAN,MXTRAN)
      INTEGER IBA2TRAN(MXDIM_BATRAN,MXTRAN)
      INTEGER IAA1TRAN(MXDIM_AATRAN,MXTRAN)
      INTEGER IAA2TRAN(MXDIM_AATRAN,MXTRAN)
      INTEGER IXETRAN(MXDIM_XEVEC,MXTRAN)

      INTEGER NDTRAN,NCTRAN,NB1TRAN,NB2TRAN,NBA1TRAN,NAA1TRAN,NXETRAN
      INTEGER NBA2TRAN,NAA2TRAN,NCATRAN

      INTEGER IOP(MXORD), IOP2(MXORD2)
      INTEGER IR1(MXORD), IR2(MXORD2), IR3(MXORD3)
      INTEGER IEX, IE1(MXORD), IET1(MXORD), IET2(MXORD2)
      INTEGER ISYMS, ISYM(MXORD), IRELAX(MXORD)

      INTEGER A, B, C, D
      PARAMETER (A = 1, B = 2, C = 3, D = 4)
      INTEGER AB, AC, AD, BC, BD, CD
      PARAMETER (AB = 1, AC = 2, BC = 3, AD = 4, BD = 5, CD = 6)
      INTEGER ABC, ABD, ACD, BCD 
      PARAMETER (ABC = 1, ABD = 2, ACD = 3, BCD = 4) 
     

      INTEGER NS2A, NS3A, NS4A, NP3AB, NP4AB, NT4ABC
      PARAMETER (NS2A = 2, NS3A = 3, NS4A = 4)
      PARAMETER (NP3AB = 3, NP4AB = 6, NT4ABC = 4)

      INTEGER ISA(NS4A), ISB(NS4A), ISC(NS4A), ISD(NS4A)
      INTEGER IPAB(NP4AB), IPCD(NP4AB)
      INTEGER IPA(NP4AB), IPB(NP4AB), IPC(NP4AB), IPD(NP4AB)
      INTEGER ITABC(NT4ABC), ITD(NT4ABC)

      DATA ISA  / A, B, C, D/
      DATA ISB  / B, A, A, A/
      DATA ISC  / C, C, B, B/
      DATA ISD  / D, D, D, C/

      DATA IPAB / AB, AC, BC, AD, BD, CD /
      DATA IPA  / A,  A,  B,  A,  B,  C  /
      DATA IPB  / B,  C,  C,  D,  D,  D  /     
      DATA IPCD / CD, BD, AD, BC, AC, AB /
      DATA IPC  / C,  B,  A,  B,  A,  A  /
      DATA IPD  / D,  D,  D,  C,  C,  B  /

      DATA ITABC / ABC, ABD, ACD, BCD /
      DATA ITD   / D,   C,   B,   A   /

      CHARACTER*8 LABSOP
      INTEGER IDXA, IDXB, IDXC, ITRAN, IDX, IDXAB, IDXABC, IVEC
      INTEGER NRELAX, ISGNSOP, ISYSOP, INUM  

* external functions:
      INTEGER IROPER
      INTEGER IROPER2 
      INTEGER IETA1
      INTEGER ICHI2
      INTEGER IR1KAPPA
      INTEGER IR1TAMP
      INTEGER IR2TAMP
      INTEGER IR3TAMP
      INTEGER IER1AMP
      INTEGER IER2AMP
      INTEGER ILRCAMP
      INTEGER ICR2AMP


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      NDTRAN   = 0
      NCTRAN   = 0
      NB1TRAN  = 0
      NB2TRAN  = 0
      NCATRAN  = 0
      NBA1TRAN = 0
      NBA2TRAN = 0
      NAA1TRAN = 0
      NAA2TRAN = 0
      NXETRAN  = 0

*---------------------------------------------------------------------*
* start loop over all requested rhs-vectors:
*---------------------------------------------------------------------*
 
      DO IVEC = IOFFV+1, IOFFV+NVEC

* eigenvectors that contribute:
      IF (NSTAT.EQ.1) THEN
        IEX = ISTAT(IVEC,1)
      END IF

* first-order operators:
      DO IDXA = 1, ORDER
        IOP(IDXA)=IROPER(LAB(IVEC,IDXA),ISYM(IDXA))
      END DO

* relaxation flags:
      IF (TYPE(1:1).EQ.'O'  .OR. TYPE(1:1).EQ.'X' .OR.
     &    TYPE(1:2).EQ.'EO' .OR. TYPE(1:2).EQ.'EX'     ) THEN
        NRELAX = 0
        DO IDXA = 1, ORDER
         IF ( LORX(IVEC,IDXA) ) THEN
           IRELAX(IDXA) = IR1KAPPA(LAB(IVEC,IDXA),
     &                            FREQ(IVEC,IDXA),ISYM(IDXA))
           NRELAX = NRELAX + 1
         ELSE
           IRELAX(IDXA) = 0
         END IF
        END DO
      ELSE
        NRELAX = 0
        DO IDXA = 1, ORDER
         IRELAX(IDXA) = 0
        END DO
      END IF

      IF (NRELAX.GT.1) THEN
         CALL QUIT('NRELAX TOO LARGE IN CC_RHS_SETUP.')
      END IF
                                                                
* second-order operators that contribute:
      IF (     (TYPE(1:1).EQ.'O'  .AND. ORDER.GE.2)
     &    .OR. (TYPE(1:2).EQ.'EO' .AND. ORDER.GE.2) ) THEN
        IDXAB  = 0
        DO IDXB = 2, ORDER
        DO IDXA = 1, IDXB-1
         IDXAB = IDXAB + 1
         IF (IRELAX(IDXA).GT.1 .OR. LPDBSOP(IOP(IDXA)) .OR.
     &       IRELAX(IDXB).GT.1 .OR. LPDBSOP(IOP(IDXB))      ) THEN
           INUM        = IROPER2(LAB(IVEC,IDXA),LAB(IVEC,IDXB),
     &                           LABSOP,ISGNSOP,ISYSOP)
           IOP2(IDXAB) = IROPER(LABSOP,ISYSOP)
         ELSE
           IOP2(IDXAB) = -1
         END IF
        END DO
        END DO
      END IF           
 
* first-order vectors that contribute:
      IF (     (TYPE(1:1).EQ.'O'  .AND. ORDER.GT.1)
     &    .OR. (TYPE(1:2).EQ.'EO' .AND. ORDER.GE.1) ) THEN
        DO IDXA = 1, ORDER
         IR1(IDXA)=IR1TAMP(LAB(IVEC,IDXA),LORX(IVEC,IDXA),
     &                     FREQ(IVEC,IDXA),ISYM(IDXA))
        END DO
      END IF
      IF (TYPE(1:3).EQ.'O1 ') THEN
        DO IDXA = 1, ORDER
         IET1(IDXA) = IETA1(LAB(IVEC,IDXA),LORX(IVEC,IDXA),
     &                      FREQ(IVEC,IDXA),ISYM(IDXA))
        END DO
      END IF
      IF (TYPE(1:2).EQ.'EO' .AND. ORDER.GT.1) THEN
        call quit('Sonia: please define LPROJ in IER1AMP call')
        DO IDXA = 1, ORDER
         IE1(IDXA)=IER1AMP(ISTAT(IVEC,1),EIGV(IVEC),ISYMS,
     &                     LAB(IVEC,IDXA),FREQ(IVEC,IDXA),ISYM(IDXA))
        END DO
      END IF
      IF (TYPE(1:2).EQ.'CO' .AND. ORDER.GT.1) THEN
        DO IDXA = 1, ORDER
         IR1(IDXA)=ILRCAMP(LAB(IVEC,IDXA),ICAU(IVEC,IDXA),ISYM(IDXA))
        END DO
      END IF

* second-order vectors that contribute:
      IF (     (TYPE(1:1).EQ.'O'  .AND. ORDER.GT.2)
     &    .OR. (TYPE(1:2).EQ.'EO' .AND. ORDER.GE.2) ) THEN
        IDXAB  = 0
        DO IDXB = 2, ORDER
        DO IDXA = 1, IDXB-1
         IDXAB = IDXAB + 1
         IR2(IDXAB) =
     &        IR2TAMP(LAB(IVEC,IDXA),LORX(IVEC,IDXA),
     &                   FREQ(IVEC,IDXA),ISYM(IDXA),
     &                LAB(IVEC,IDXB),LORX(IVEC,IDXB),
     &                   FREQ(IVEC,IDXB),ISYM(IDXB))
        END DO
        END DO
      END IF
      IF (TYPE(1:2).EQ.'O2' .AND. IOP2(AB).GT.0) THEN
         IET2(IDXAB) = ICHI2(LAB(IVEC,1),LORX(IVEC,1),
     &                              FREQ(IVEC,1),ISYM(1),
     &                       LAB(IVEC,2),LORX(IVEC,2),
     &                              FREQ(IVEC,2),ISYM(2))
      END IF                                              

* third-order vectors that contribute:
      IF (ORDER .GT. 3) THEN
       IDXABC = 0
       DO IDXC = 3, ORDER
       DO IDXB = 2, IDXC-1
       DO IDXA = 1, IDXB-1
        IDXABC = IDXABC + 1
        IR3(IDXABC) =
     &            IR3TAMP(LAB(IVEC,IDXA),FREQ(IVEC,IDXA),ISYM(IDXA),
     &                    LAB(IVEC,IDXB),FREQ(IVEC,IDXB),ISYM(IDXB),
     &                    LAB(IVEC,IDXC),FREQ(IVEC,IDXC),ISYM(IDXC))
       END DO
       END DO
       END DO
      END IF

*---------------------------------------------------------------------*
* set up list of D matrix transformations:
*---------------------------------------------------------------------*
        IF ( TYPE(1:2).EQ.'O4' ) THEN
          NDTRAN = NDTRAN + 1
          IDTRAN(1,NDTRAN) = IR1(A)
          IDTRAN(2,NDTRAN) = IR1(B)
          IDTRAN(3,NDTRAN) = IR1(C)
          IDTRAN(4,NDTRAN) = IR1(D)
          IDTRAN(5,NDTRAN) = IVEC
        END IF

*---------------------------------------------------------------------*
* set up list of C matrix transformations:
*---------------------------------------------------------------------*
        IF      ( TYPE(1:2).EQ.'O4' ) THEN
          DO IDX = 1, NP4AB
            NCTRAN = NCTRAN + 1
            ICTRAN(1,NCTRAN) = IR2(IPAB(IDX))
            ICTRAN(2,NCTRAN) = IR1(IPC(IDX))
            ICTRAN(3,NCTRAN) = IR1(IPD(IDX))
            ICTRAN(4,NCTRAN) = IVEC
          END DO
        ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
          NCTRAN = NCTRAN + 1
          ICTRAN(1,NCTRAN) = IR1(1)
          ICTRAN(2,NCTRAN) = IR1(2)
          ICTRAN(3,NCTRAN) = IR1(3)
          ICTRAN(4,NCTRAN) = IVEC
        ELSE IF ( TYPE(1:3).EQ.'EO2' ) THEN
          NCTRAN = NCTRAN + 1
          ICTRAN(1,NCTRAN) = IR1(1)
          ICTRAN(2,NCTRAN) = IR1(2)
          ICTRAN(3,NCTRAN) = IEX
          ICTRAN(4,NCTRAN) = IVEC
        END IF

*---------------------------------------------------------------------*
* set up list of B matrix transformations
*---------------------------------------------------------------------*
        IF      ( TYPE(1:2).EQ.'O4' ) THEN
          DO IDX = 1, NT4ABC
            NB1TRAN = NB1TRAN + 1
            IB1TRAN(1,NB1TRAN) = IR3(ITABC(IDX))
            IB1TRAN(2,NB1TRAN) = IR1(ITD(IDX))
            IB1TRAN(3,NB1TRAN) = IVEC
          END DO

          DO IDX = 1, NP4AB
            NB2TRAN = NB2TRAN + 1
            IB2TRAN(1,NB2TRAN) = IR2(IPAB(IDX))
            IB2TRAN(2,NB2TRAN) = IR2(IPCD(IDX))
            IB2TRAN(3,NB2TRAN) = IVEC
          END DO
        ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
          DO IDX = 1, NP3AB
            NB1TRAN = NB1TRAN + 1
            IB1TRAN(1,NB1TRAN) = IR2(IPAB(IDX))
            IB1TRAN(2,NB1TRAN) = IR1(IPC(IDX))
            IB1TRAN(3,NB1TRAN) = IVEC
          END DO
        ELSE IF ( TYPE(1:2).EQ.'O2' ) THEN
          NB1TRAN = NB1TRAN + 1
          IB1TRAN(1,NB1TRAN) = IR1(1)
          IB1TRAN(2,NB1TRAN) = IR1(2)
          IB1TRAN(3,NB1TRAN) = IVEC
        ELSE IF ( TYPE(1:3).EQ.'EO2' ) THEN
          NB1TRAN = NB1TRAN + 1
          IB1TRAN(1,NB1TRAN) = IR2(1)
          IB1TRAN(2,NB1TRAN) = IEX
          IB1TRAN(3,NB1TRAN) = IVEC

          DO IDX = 1, NS2A
            NB2TRAN = NB2TRAN + 1
            IB2TRAN(1,NB2TRAN) = IR1(ISA(IDX))
            IB2TRAN(2,NB2TRAN) = IE1(ISB(IDX))
            IB2TRAN(3,NB2TRAN) = IVEC
          END DO
        ELSE IF ( TYPE(1:3).EQ.'EO1' ) THEN
          NB1TRAN = NB1TRAN + 1
          IB1TRAN(1,NB1TRAN) = IR1(1)
          IB1TRAN(2,NB1TRAN) = IEX
          IB1TRAN(3,NB1TRAN) = IVEC
        ELSE IF ( TYPE(1:3).EQ.'CO2' ) THEN
          NB1TRAN = NB1TRAN + 1
          IB1TRAN(1,NB1TRAN) = IR1(1)
          IB1TRAN(2,NB1TRAN) = IR1(2)
          IB1TRAN(3,NB1TRAN) = IVEC
        END IF

*---------------------------------------------------------------------*
* set up list of C{O} matrix transformations:
*---------------------------------------------------------------------*
        IF      ( TYPE(1:2).EQ.'O4' ) THEN
          DO IDX = 1, NS4A
            IF (IRELAX(ISA(IDX)).GT.0) THEN
              NCATRAN = NCATRAN + 1
              ICATRAN(1,NCATRAN) = IOP(ISA(IDX))
              ICATRAN(2,NCATRAN) = IR1(ISB(IDX))
              ICATRAN(3,NCATRAN) = IR1(ISC(IDX))
              ICATRAN(4,NCATRAN) = IR1(ISD(IDX))
              ICATRAN(5,NCATRAN) = IVEC
              ICATRAN(6,NCATRAN) = IRELAX(ISA(IDX))
              ICATRAN(7,NCATRAN) = 0
              ICATRAN(8,NCATRAN) = 0
              ICATRAN(9,NCATRAN) = 0
            END IF
          END DO
        END IF                

*---------------------------------------------------------------------*
* set up list of B{O} matrix transformations:
*---------------------------------------------------------------------*
        IF      ( TYPE(1:2).EQ.'O4' ) THEN
          DO IDX = 1, NP4AB
            NBA1TRAN = NBA1TRAN + 1
            IBA1TRAN(1,NBA1TRAN) = IOP(IPC(IDX))
            IBA1TRAN(2,NBA1TRAN) = IR2(IPAB(IDX))
            IBA1TRAN(3,NBA1TRAN) = IR1(IPD(IDX))
            IBA1TRAN(4,NBA1TRAN) = IVEC
            IBA1TRAN(5,NBA1TRAN) = IRELAX(IPC(IDX))
            IBA1TRAN(6,NBA1TRAN) = 0
            IBA1TRAN(7,NBA1TRAN) = 0
            IBA1TRAN(8,NBA1TRAN) = 0

            NBA1TRAN = NBA1TRAN + 1
            IBA1TRAN(1,NBA1TRAN) = IOP(IPD(IDX))
            IBA1TRAN(2,NBA1TRAN) = IR2(IPAB(IDX))
            IBA1TRAN(3,NBA1TRAN) = IR1(IPC(IDX))
            IBA1TRAN(4,NBA1TRAN) = IVEC
            IBA1TRAN(5,NBA1TRAN) = IRELAX(IPD(IDX))
            IBA1TRAN(6,NBA1TRAN) = 0
            IBA1TRAN(7,NBA1TRAN) = 0
            IBA1TRAN(8,NBA1TRAN) = 0

            IF (IOP2(IPAB(IDX)).GT.0) THEN
              NBA2TRAN = NBA2TRAN + 1
              IBA2TRAN(1,NBA2TRAN) = IOP2(IPAB(IDX))
              IBA2TRAN(2,NBA2TRAN) = IR1(IPC(IDX))
              IBA2TRAN(3,NBA2TRAN) = IR1(IPD(IDX))
              IBA2TRAN(4,NBA2TRAN) = IVEC
              IBA2TRAN(5,NBA2TRAN) = IRELAX(IPA(IDX))
              IBA2TRAN(6,NBA2TRAN) = IRELAX(IPB(IDX))
              IBA2TRAN(7,NBA2TRAN) = 0
              IBA2TRAN(8,NBA2TRAN) = 0
            END IF                        
          END DO
        ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
          DO IDX = 1, NS3A
            NBA1TRAN = NBA1TRAN + 1
            IBA1TRAN(1,NBA1TRAN) = IOP(ISA(IDX))
            IBA1TRAN(2,NBA1TRAN) = IR1(ISB(IDX))
            IBA1TRAN(3,NBA1TRAN) = IR1(ISC(IDX))
            IBA1TRAN(4,NBA1TRAN) = IVEC
            IBA1TRAN(5,NBA1TRAN) = IRELAX(ISA(IDX))
            IBA1TRAN(6,NBA1TRAN) = 0
            IBA1TRAN(7,NBA1TRAN) = 0
            IBA1TRAN(8,NBA1TRAN) = 0
          END DO
        ELSE IF ( TYPE(1:3).EQ.'EO2' ) THEN
          DO IDX = 1, NS2A
            NBA1TRAN = NBA1TRAN + 1
            IBA1TRAN(1,NBA1TRAN) = IOP(ISA(IDX))
            IBA1TRAN(2,NBA1TRAN) = IR1(ISB(IDX))
            IBA1TRAN(3,NBA1TRAN) = IEX
            IBA1TRAN(4,NBA1TRAN) = IVEC
            IBA1TRAN(5,NBA1TRAN) = IRELAX(ISA(IDX))
            IBA1TRAN(6,NBA1TRAN) = 0
            IBA1TRAN(7,NBA1TRAN) = 0
            IBA1TRAN(8,NBA1TRAN) = 0
          END DO
        END IF

*---------------------------------------------------------------------*
* set up list of A{O} vector calculations:
*---------------------------------------------------------------------*
        IF      ( TYPE(1:2).EQ.'O4' ) THEN
          DO IDX = 1, NT4ABC
            NAA1TRAN = NAA1TRAN + 1
            IAA1TRAN(1,NAA1TRAN) = IOP(ITD(IDX))
            IAA1TRAN(2,NAA1TRAN) = IR3(ITABC(IDX))
            IAA1TRAN(3,NAA1TRAN) = IVEC
            IAA1TRAN(4,NAA1TRAN) = IRELAX(ITD(IDX))
            IAA1TRAN(5,NAA1TRAN) = 0
            IAA1TRAN(6,NAA1TRAN) = 0
            IAA1TRAN(7,NAA1TRAN) = 0
          END DO
          DO IDX = 1, NP4AB
            IF (IOP2(IPAB(IDX)).GT.0) THEN
              NAA2TRAN = NAA2TRAN + 1
              IAA2TRAN(1,NAA2TRAN) = IOP2(IPAB(IDX))
              IAA2TRAN(2,NAA2TRAN) = IR2(IPCD(IDX))
              IAA2TRAN(3,NAA2TRAN) = IVEC
              IAA2TRAN(4,NAA2TRAN) = IRELAX(IPA(IDX))
              IAA2TRAN(5,NAA2TRAN) = IRELAX(IPB(IDX))
              IAA2TRAN(6,NAA2TRAN) = 0
              IAA2TRAN(7,NAA2TRAN) = 0
            END IF
          END DO                       
        ELSE IF ( TYPE(1:2).EQ.'O3' ) THEN
          DO IDX = 1, NP3AB
            NAA1TRAN = NAA1TRAN + 1
            IAA1TRAN(1,NAA1TRAN) = IOP(IPC(IDX))
            IAA1TRAN(2,NAA1TRAN) = IR2(IPAB(IDX))
            IAA1TRAN(3,NAA1TRAN) = IVEC
            IAA1TRAN(4,NAA1TRAN) = IRELAX(IPC(IDX))
            IAA1TRAN(5,NAA1TRAN) = 0
            IAA1TRAN(6,NAA1TRAN) = 0
            IAA1TRAN(7,NAA1TRAN) = 0

            IF (IOP2(IPAB(IDX)).GT.0) THEN
              NAA2TRAN = NAA2TRAN + 1
              IAA2TRAN(1,NAA2TRAN) = IOP2(IPAB(IDX))
              IAA2TRAN(2,NAA2TRAN) = IR1(IPC(IDX))
              IAA2TRAN(3,NAA2TRAN) = IVEC
              IAA2TRAN(4,NAA2TRAN) = IRELAX(IPA(IDX))
              IAA2TRAN(5,NAA2TRAN) = IRELAX(IPB(IDX))
              IAA2TRAN(6,NAA2TRAN) = 0
              IAA2TRAN(7,NAA2TRAN) = 0
            END IF                            
          END DO
        ELSE IF ( TYPE(1:2).EQ.'O2' ) THEN
          DO IDX = 1, NS2A
            NAA1TRAN = NAA1TRAN + 1
            IAA1TRAN(1,NAA1TRAN) = IOP(ISB(IDX))
            IAA1TRAN(2,NAA1TRAN) = IR1(ISA(IDX))
            IAA1TRAN(3,NAA1TRAN) = IVEC
            IAA1TRAN(4,NAA1TRAN) = IRELAX(ISB(IDX))
            IAA1TRAN(5,NAA1TRAN) = 0
            IAA1TRAN(6,NAA1TRAN) = 0
            IAA1TRAN(7,NAA1TRAN) = 0
          END DO
        ELSE IF ( TYPE(1:3).EQ.'EO2' ) THEN
          DO IDX = 1, NS2A
            NAA1TRAN = NAA1TRAN + 1
            IAA1TRAN(1,NAA1TRAN) = IOP(ISB(IDX))
            IAA1TRAN(2,NAA1TRAN) = IE1(ISA(IDX))
            IAA1TRAN(3,NAA1TRAN) = IVEC
            IAA1TRAN(4,NAA1TRAN) = IRELAX(ISB(IDX))
            IAA1TRAN(5,NAA1TRAN) = 0
            IAA1TRAN(6,NAA1TRAN) = 0
            IAA1TRAN(7,NAA1TRAN) = 0
          END DO
          IF (IOP2(AB).GT.0) THEN
            NAA2TRAN = NAA2TRAN + 1
            IAA2TRAN(1,NAA2TRAN) = IOP2(AB)
            IAA2TRAN(2,NAA2TRAN) = IEX
            IAA2TRAN(3,NAA2TRAN) = IVEC
            IAA2TRAN(4,NAA2TRAN) = IRELAX(A)
            IAA2TRAN(5,NAA2TRAN) = IRELAX(B)
            IAA2TRAN(6,NAA2TRAN) = 0
            IAA2TRAN(7,NAA2TRAN) = 0
          END IF                               
        ELSE IF ( TYPE(1:3).EQ.'EO1' ) THEN
          NAA1TRAN = NAA1TRAN + 1
          IAA1TRAN(1,NAA1TRAN) = IOP(1)
          IAA1TRAN(2,NAA1TRAN) = IEX
          IAA1TRAN(3,NAA1TRAN) = IVEC
          IAA1TRAN(4,NAA1TRAN) = IRELAX(1)
          IAA1TRAN(5,NAA1TRAN) = 0
          IAA1TRAN(6,NAA1TRAN) = 0
          IAA1TRAN(7,NAA1TRAN) = 0
        ELSE IF ( TYPE(1:3).EQ.'CO2' ) THEN
          DO IDX = 1, NS2A
            IF (ICAU(IVEC,ISB(IDX)).EQ.0) THEN
              NAA1TRAN = NAA1TRAN + 1
              IAA1TRAN(1,NAA1TRAN) = IOP(ISB(IDX))
              IAA1TRAN(2,NAA1TRAN) = IR1(ISA(IDX))
              IAA1TRAN(3,NAA1TRAN) = IVEC
              IAA1TRAN(4,NAA1TRAN) = IRELAX(ISB(IDX))
              IAA1TRAN(5,NAA1TRAN) = 0
              IAA1TRAN(6,NAA1TRAN) = 0
              IAA1TRAN(7,NAA1TRAN) = 0
            END IF
          END DO
        END IF

*---------------------------------------------------------------------*
* set up list of Xi{O} vector calculations:
* Note, that we set up here a list for the simultaneous calculation
* of the first-order xi "O1" and the first-order eta "X1" vectors.
* Xi and eta vectors are only precalculated for orbital relaxed
* "operators" or for field-dependent basis sets. For simple unrelaxed
* one-electron perturbations they are calculated on the fly when needed
*---------------------------------------------------------------------*
        IF ( TYPE(1:3).EQ.'O1 ') THEN
C         IF ( IRELAX(A).EQ.1 .OR. LPDBSOP(IOP(A)) ) THEN
            NXETRAN = NXETRAN + 1
            IXETRAN(1,NXETRAN) = IOP(A)
            IXETRAN(2,NXETRAN) = 0      ! L0 for first-order ETA vec.
            IXETRAN(3,NXETRAN) = IVEC
            IXETRAN(4,NXETRAN) = IET1(A)
            IXETRAN(5,NXETRAN) = IRELAX(A)
            IXETRAN(6,NXETRAN) = 0
            IXETRAN(7,NXETRAN) = 0
            IXETRAN(8,NXETRAN) = 0
C         END IF
        ELSE IF ( TYPE(1:2).EQ.'O2' ) THEN
          IF ( IOP2(AB).GT.0 ) THEN
            NXETRAN = NXETRAN + 1
            IXETRAN(1,NXETRAN) = IOP2(AB)
            IXETRAN(2,NXETRAN) = 0      ! L0 for second-order ETA vec.
            IXETRAN(3,NXETRAN) = IVEC
            IXETRAN(4,NXETRAN) = IET2(AB)
            IXETRAN(5,NXETRAN) = IRELAX(A)
            IXETRAN(6,NXETRAN) = IRELAX(B)
            IXETRAN(7,NXETRAN) = 0
            IXETRAN(8,NXETRAN) = 0
          END IF                    
        ELSE IF ( TYPE(1:3).EQ.'CO1') THEN
            NXETRAN = NXETRAN + 1
            IXETRAN(1,NXETRAN) = IOP(A)
            IXETRAN(2,NXETRAN) = 0      ! L0 for first-order ETA vec.
            IXETRAN(3,NXETRAN) = IVEC
            IXETRAN(4,NXETRAN) = -1
            IXETRAN(5,NXETRAN) = IRELAX(A)
            IXETRAN(6,NXETRAN) = 0
            IXETRAN(7,NXETRAN) = 0
            IXETRAN(8,NXETRAN) = 0
        END IF

*---------------------------------------------------------------------*
* end loop over all requested rhs vectors
*---------------------------------------------------------------------*
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      WRITE(LUPRI,'(/,/3X,A,I3,I2,3A)') 'For the requested',NVEC,ORDER,
     &      'th.-order amplitude rhs vectors "',TYPE,'".'
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',NDTRAN,            ' D matrix transformations ',
     &   ' - ',NCTRAN,            ' C matrix transformations ',
     &   ' - ',NB1TRAN+NB2TRAN,   ' B matrix transformations ',
     &   ' - ',NCATRAN,           ' C{O} matrix transformations ',
     &   ' - ',NBA1TRAN+NBA2TRAN, ' B{O} matrix transformations ',
     &   ' - ',NAA1TRAN+NAA2TRAN, ' A{O} matrix transformations ',
     &   ' - ',NXETRAN,           'Xi{O} vector calculations    ' 
      IF (NEW_RHS) WRITE(LUPRI,'(14X,A)') 
     &   '(A{O} matrix included in B matrix)'
      WRITE(LUPRI,'(3X,A/,/)') 'will be performed.'


* D matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of D matrix transformations:'
        DO ITRAN = 1, NDTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (IDTRAN(IDX,ITRAN),IDX=1,5)
        END DO
        WRITE (LUPRI,*)
      END IF

* C matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of C matrix transformations:'
        DO ITRAN = 1, NCTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (ICTRAN(IDX,ITRAN),IDX=1,4)
        END DO
        WRITE (LUPRI,*)
      END IF

* B matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of B matrix transformations (type1):'
        DO ITRAN = 1, NB1TRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IB1TRAN(IDX,ITRAN),IDX=1,3)
        END DO
        WRITE (LUPRI,*) 'List of B matrix transformations (type2):'
        DO ITRAN = 1, NB2TRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IB2TRAN(IDX,ITRAN),IDX=1,3)
        END DO
        WRITE (LUPRI,*)
      END IF

* C{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of C{O} matrix transformations:'
        DO ITRAN = 1, NCATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (ICATRAN(IDX,ITRAN),IDX=1,5)
        END DO
        WRITE (LUPRI,*)
      END IF

* B{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of B{O} matrix transformations (type1):'
        DO ITRAN = 1, NBA1TRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (IBA1TRAN(IDX,ITRAN),IDX=1,4)
        END DO
        WRITE (LUPRI,*) 'List of B{O} matrix transformations (type 2):'
        DO ITRAN = 1, NBA2TRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (IBA2TRAN(IDX,ITRAN),IDX=1,4)
        END DO                          
        WRITE (LUPRI,*)
      END IF

* A{O} matrix calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of A{O} matrix transformations (type1):'
        DO ITRAN = 1, NAA1TRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IAA1TRAN(IDX,ITRAN),IDX=1,3)
        END DO
        WRITE (LUPRI,*) 'List of A{O} matrix transformations (type 2):'
        DO ITRAN = 1, NAA2TRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IAA2TRAN(IDX,ITRAN),IDX=1,3)
        END DO                
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

* Xi{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of Xi{O} vector calculations:'
        DO ITRAN = 1, NXETRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IXETRAN(IDX,ITRAN),IDX=1,2)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF


      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_RHS_SETUP                         *
*---------------------------------------------------------------------*
