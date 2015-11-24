!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
c===========================================================================
c/* Comdeck comments */
c  950505-hjaaj: changed print such that no print if IPRINT .le. 2
c  *** PRT_COMMENTS
c  910207  PSORG/PSORT error for cases in which NBLOCK is very large
c          but NUMP is relatively small (many centres with small
c          basis sets).  This is a failure to reset the slice and pass
c          counts (if necessary) if the number of blocks in one
c          slice exceeds MAXLOC.  This involves extending the tests
c          in the loop nest in PSORG and changing pointer assignments
c          in PSORT
c  910105  Corrected error in PSORG (PSORT was called with MEMS
c          as available memory instead of MEMMAX) and two errors
c          in PSORT.  First was test on computed ITAB(NTAB), which
c          was reset if .LT. 0 --- should have been .LE. 0.
c          Second was that if LASTAD() was -1, the whole write of
c          the P block was skipped --- must write zero block.  We
c          need to look at compression for this.  Also, HJJ:s comment
c          of 890529 should be looked at closely.
c  *** hj_comments
c  900508: Corrected error in PSORT. IPQRS set to 0 at 1000
c  890529: something is very fishy about PSORT for multiple passes,
c          (I think), from a quick glance it seems as if IPQRS is
c          only compared to a lower limit:
c                  IF (IPQRS .LT. IPQRST) GOTO 40
c          which must give problems if multiple passes.
c          Determination of NPASS in PSIZE was wrong, has been corrected.
c===========================================================================
C  /* Deck ptran */
      SUBROUTINE PTRAN(NODPTR,WORK,LWORK,JPRINT,ANTI,PANTI,
     &                 DIA2SO,ZFS2EL)
C
C            ******      M O L E C U L E - P T R A N      ******
C
C        Transformation of MCSCF and CI 2-matrix from an MO basis
C        to a symmetry orbital basis.
C
C     ANTI : P(ij,kl) = -P(ji,kl) = -P(ij,lk)
C     PANTI: P(ij,kl) = -P(kl,ij)
C
C     DIA2SO, ZFS2EL: get corresponding special 2-matrix with MAKPVM
C
C        ******      VERSION 1.0      RELEASE 820423      ******
C        (Peter Taylor)
C        ******      VERSION 2.0      RELEASE 140731      ******
C        (Hans Joergen Aa. Jensen)
C
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
      LOGICAL   NODPTR, ANTI, PANTI, DIA2SO, ZFS2EL
      REAL*8    WORK(LWORK)
C
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NBAS(8), NASHT, NBAST,
     1                IASH(8), IBAS(8), IASYM(MXCORB), IOSYM(MXCORB),
     2                LSMAX, MOFF(MXCORB), NCMOT, MXBAS1, MXBAS2, NNASHX
      COMMON /PTRINF/ THR1, THR2, DAWEIG, SQWEIG, IPRINT
      COMMON /PTRSYM/ NSYM, MUL(8,8), NENT(8,8), ISYOFF(8,8)
C
C
      CALL QENTER('PTRAN')
C
C     Initialize /PTRFIL/
C
      LUSCR  = -9999
      LUDA   = -9999
      LUPSO  = -9999
C
      LDAMAX = 8*1024
      MX1BUF = LDAMAX/2 - 2
      MX2BUF = MX1BUF
C
C     Initialize /PTRINF/
C
      THR1   = 1.D-10
      THR2   = 1.D-10
      DAWEIG = 1.D0
      SQWEIG = 1.D0
      TIMSTR = SECOND()
      IPRINT = JPRINT
      IF (IPRINT .GT. 02) WRITE(LUPRI,2000)
      CALL PTRMOS
      IF (IPRINT .GT. 02) CALL PTPRNT
C
C     *************************************
C     ***** FIRST HALF-TRANSFORMATION *****
C     *************************************
C
      TIMOLD = SECOND()
C
C     ALLOCATE WORKING STORAGE
C
C     IT, P2, and P4V
C
      ISCM1 = 0
      ISCM2 = ISCM1 + NBAST
      ISCM3 = ISCM2 + MXBAS2
      ISCMT = ISCM3 + MXBAS1
C
C     CMO  and PV
C
      NTU   = (NASHT*(NASHT+1))/2
      NTUVX = NTU*NTU
      ILCM1 = ISCMT
      ILCM2 = ILCM1 + NCMOT
      ILCM3 = ILCM2 + NTUVX
C
C     IBUF and BUF
C
      ILCM4 = ILCM3 + MX1BUF
      ILCMT = ILCM4 + MX1BUF + 1
      IF (ILCMT.GT.LWORK) CALL STOPIT('PTRAN','PT1TRA',ILCMT,LWORK)
C
      KWRK = ILCMT
      LWRK = LWORK  - ILCMT
      CALL PT1TRA (WORK(ILCM1+1), WORK(ILCM2+1), WORK(ISCM2+1),
     &             WORK(ISCM3+1), WORK(ILCM3+1), WORK(ILCM4+1),
     &             WORK(ISCM1+1), WORK(KWRK),LWRK,ANTI,PANTI,
     &             DIA2SO,ZFS2EL)
C     CALL PT1TRA (CMO,PV,P2,P4V,IBUF,BUF,IT,ANTI,PANTI,DIA2SO,ZFS2EL)
C
      IF (IPRINT .GT. 02) THEN
         TIMNEW = SECOND()
         TIME = TIMNEW - TIMOLD
         WRITE(LUPRI,3000) TIME
         CALL FLSHFO(LUPRI)
         TIMOLD = TIMNEW
      END IF
 3000 FORMAT(/,'  Time in first half-transformation is ',
     *       F7.2,' seconds')
C
C     ************************************************
C     ***** SORTING OF HALF-TRANSFORMED ELEMENTS *****
C     ************************************************
C
C
C     REALLOCATE CORE FOR SORT OF HALF-TRANSFORMED 2-MATRIX ELEMENTS
C
C     Work space for IT (SCM1) is reused.
C
C     IBUF and BUF
C
      ISCM3 = ISCM2 + MX1BUF
      ISCMT = ISCM3 + MX1BUF + 1
      IF (ISCMT .GT. LWORK) CALL STOPIT('PTRAN','PT1SOR',ISCMT,LWORK)
C
      ILCM1 = ISCMT
      MEMS  = LWORK - ILCM1
      MEMT  = MEMS  - NCMOT
C
      CALL PT1SOR(WORK(ILCM1+1), WORK(ILCM1+1), WORK(ISCM2+1),
     *            WORK(ISCM3+1), WORK(ISCM1+1))
C     CALL PT1SOR(ISORT,SORT,IBUF,BUF,IT)
C
      IF (IPRINT .GT. 02) THEN
         TIMNEW = SECOND()
         TIME   = TIMNEW - TIMOLD
         WRITE(LUPRI,3100) TIME
         CALL FLSHFO(LUPRI)
         TIMOLD = TIMNEW
      END IF
 3100 FORMAT(/,'  Time for sorting half-transformed matrix is ',
     *       F7.2,' seconds')
C
C     **************************************
C     ***** SECOND HALF-TRANSFORMATION *****
C     **************************************
C
C
C     REALLOCATE STORAGE
C
C     Work space for IT (SCM1) is reused.
C
C     IBIN AND BIN
C
      ISCM3 = ISCM2 + MXABUF
      ISCM4 = ISCM3 + MXABUF + 2
C
C     IBUF and BUF
C
      ISCM5 = ISCM4 + MX2BUF
      ISCM6 = ISCM5 + MX2BUF + 1
C
C     P4T and P2
C
      ISCM7 = ISCM6 + MXBAS1
      ISCMT = ISCM7 + MXBAS2
      IF (ISCMT .GT. LWORK) CALL STOPIT('PTRAN','PT2TRA',ISCMT,LWORK)
C
C     CMO and PV (which uses rest of work space).
C
      ILCM1 = ISCMT
      ILCM2 = ILCM1 + NCMOT
      MEMS  = LWORK - ISCMT
      MEMT  = MEMS  - NCMOT
C
      CALL PT2TRA (WORK(ILCM1+1), WORK(ISCM6+1), WORK(ISCM7+1),
     *             WORK(ILCM2+1), WORK(ISCM2+1), WORK(ISCM3+1),
     *             WORK(ISCM4+1), WORK(ISCM5+1), WORK(ISCM1+1),
     *             ANTI,PANTI)
C     CALL PT2TRA (CMO,P4T,P2,PV,IBIN,BIN,IBUF,BUF,IT,ANTI,PANTI)
C
      TIMNEW = SECOND()
      TIME = TIMNEW - TIMOLD
      IF (IPRINT .GT. 02) WRITE(LUPRI,3200) TIME
 3200 FORMAT(/,'  Time in second half-transformation is ',
     *       F7.2,' seconds')
C
C     Test densities
C
      IF (NODPTR) CALL PTRTST(WORK,LWORK,ANTI,PANTI,DIA2SO,ZFS2EL,
     &                        MX2BUF,IPRINT)
C
C     ***********************************
C     ***** Transformation complete *****
C     ***********************************
C
      CALL GPCLOSE(LUDA,'DELETE')
      IF (IPRINT .GT. 1) THEN
         TIMEND = SECOND()
         TIME = TIMEND - TIMSTR
         WRITE(LUPRI,3300) TIME
      END IF
 3300 FORMAT(/,'  CPU time in PTRAN is',F8.2,' seconds')
      CALL QEXIT('PTRAN')
      RETURN
C
 2000 FORMAT(//,25X,'M O L E C U L E - P T R A N',//,
     1       10X,' Transformation of MCSCF and CI 2-matrices from a',
     2           ' molecular orbital',/,
     3       10X,' to a symmetry orbital basis',//,
     4       30X,'Peter R. Taylor',/,
     5       30X,'CSIRO Division of Chemical Physics',/,
     6       30X,'P.O. Box 160',/,
     7       30X,'Clayton, Victoria, 3168',/,
     8       30X,'Australia',/,
     9       50X,'Feb - Apr 1982',///)
      END
C  /* Deck ptrmos */
      SUBROUTINE PTRMOS
C
C....    READ ORBITAL SYMMETRY DATA FROM INTERFACE
C
C
C....    ******      VERSION 1.0      RELEASE 820423      ******
C
#include "implicit.h"
C
#include "priunit.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NBAS(8), NASHT, NBAST,
     1                IASH(8), IBAS(8), IASYM(MXCORB), IOSYM(MXCORB),
     2                LSMAX, MOFF(MXCORB), NCMOT, MXBAS1, MXBAS2, NNASHX
      COMMON /PTRINF/ THR1, THR2, DAWEIG, SQWEIG, IPRINT
      COMMON /PTRSYM/ NSYM, MUL(8,8), NENT(8,8), ISYOFF(8,8)
      INTEGER    NDUM(8)
#include "abainf.h"
#include "inftap.h"
#include "ccom.h"
      CALL QENTER('PTRMOS')
      THR1 = THRS
      THR2 = THRS
C
C     Modified to read full symmetry data     881006  PRT
C     Used from LUSIFC: NCMOT_I, NNASHX, NSYM, MUL, NISH, NASH, NORB, NBAS
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) NISHI,NASHI,NOCCI,NORBI,NBASI,NCONF,NWOPT,NWOPH,
     *              NCDETS,NCMOT_I,NNASHX,NNASHY,NNORBT,N2ORBT,
     *              NSYM, MUL, NDUM, NDUM,
     *              NISH, NASH, NORB, NBAS
C
C....    COMPUTE ORBITAL DATA
C
      NASHT = 0
      NBAST = 0
      NBASX = 0
      NCMOT = 0
      MXBAS1 = 0
      MXBAS2 = 0
      NAX = 0
      DO 10 ISYM = 1,NSYM
         IASH(ISYM) = NASHT
         IBAS(ISYM) = NBAST
         NBASI = NBAS(ISYM)
      IF (NBASI .EQ. 0) GOTO 10
         NASHI = NASH(ISYM)
         NORBI = NORB(ISYM)
         NASHT = NASHT + NASHI
         NBAST = NBAST + NBASI
         NCMOT = NCMOT + NBASI*NORBI
         MXBAS1 = MAX(NBASI,MXBAS1)
         DO 20 I = 1,NBASI
            NBASX = NBASX + 1
            IOSYM(NBASX) = ISYM
 20      CONTINUE
         DO 30 JSYM = 1,ISYM
            NBASJ = NBAS(JSYM)
         IF (NBASJ .EQ. 0) GOTO 30
            NBASIJ = NBASI*NBASJ
            IF (ISYM .EQ. JSYM) NBASIJ = (NBASI*(NBASI+1))/2
            MXBAS2 = MAX(NBASIJ,MXBAS2)
 30      CONTINUE
         DO 40 I = 1,NASHI
            NAX = NAX + 1
            IASYM(NAX) = ISYM
 40      CONTINUE
 10   CONTINUE
      IF (NCMOT .NE. NCMOT_I) THEN
         WRITE (LUPRI,'(//A/,(A,I7))')
     *      ' PTRMOS ERROR, NCMOT .ne. NCMOT_I',
     *      ' NCMOT =',NCMOT, 'NCMOT_I =',NCMOT_I
         CALL QUIT('PTRMOS ERROR, NCMOT(calc.) .ne. NCMOT(from file).')
      END IF
      CALL QEXIT('PTRMOS')
      RETURN
      END
C  /* Deck ptprnt */
      SUBROUTINE PTPRNT
C
C....    PRINT DATA READ FROM INPUT AND FROM INTERFACE, TOGETHER
C....    WITH DATA INITIALIZED IN COMMON
C
C
C....    ******      VERSION 1.0      RELEASE 820423      ******
C
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
C
#include "inftap.h"
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NBAS(8), NASHT, NBAST,
     1                IASH(8), IBAS(8), IASYM(MXCORB), IOSYM(MXCORB),
     2                LSMAX, MOFF(MXCORB), NCMOT, MXBAS1, MXBAS2, NNASHX
      COMMON /PTRINF/ THR1, THR2, DAWEIG, SQWEIG, IPRINT
      COMMON /PTRSYM/ NSYM, MUL(8,8), NENT(8,8), ISYOFF(8,8)
      CALL QENTER('PTPRNT')
C
      MXAB = (LDAMAX - 2)/2
      WRITE(LUPRI,2050) 'MCSCF'
      ! write(lupri,*) 'ldamax',ldamax
      WRITE(LUPRI,2200) MX1BUF, MXAB, MX2BUF
      WRITE(LUPRI,2300) MAXCHN
      WRITE(LUPRI,2500) SQWEIG, DAWEIG
      WRITE(LUPRI,2600) THR1, THR2
      WRITE(LUPRI,2700) IPRINT
      WRITE(LUPRI,2800) NSYM
      DO, I = 1,NSYM
        WRITE(LUPRI,2850) (MUL(I,J),J = 1,NSYM)
      END DO
      WRITE(LUPRI,2900) (I,I = 1,NSYM)
      WRITE(LUPRI,2910) (NISH(I),I = 1,NSYM)
      WRITE(LUPRI,2920) (NASH(I),I = 1,NSYM)
      WRITE(LUPRI,2930) (NBAS(I),I = 1,NSYM)
      IF (IPRINT .GE. 25) THEN
         WRITE(LUPRI,3000) (I,IASYM(I),I = 1,NASHT)
         WRITE(LUPRI,3100) (I,IOSYM(I),I = 1,NBAST)
      END IF
      IF (IPRINT .GE. 5) THEN
         WRITE(LUPRI,3200) (IASH(I),I = 1,NSYM)
         WRITE(LUPRI,3300) (IBAS(I),I = 1,NSYM)
      END IF
      CALL FLSHFO(LUPRI)
      CALL QEXIT('PTPRNT')
      RETURN
 2050 FORMAT(//,'    Transformation of ',A,' 2-matrix')
 2100 FORMAT(//,'    Unit numbers:',
     1       //,'    Function',23X,'Unit',
     2       //,4X,'input',26X,I4,/,
     3       4X,'output',25X,I4,/,
     4       4X,'interface',22X,I4,/,
     5       4X,'half-transformed 2-matrix',6X,I4,/,
     6       4X,'direct access sort',13X,I4,/,
     7       4X,'transformed 2-matrix',11X,I4)
 2200 FORMAT(//,'    Number of 2-matrix elements per buffer:',//,
     2        I16,/,
     3        I16,' (default maximum only)',/,
     4        I16)
 2300 FORMAT(//,'    Maximum number of bin sort chains',I5)
 2500 FORMAT(//,'    I/O weight factors ',
     1          '(expressed relative to DA write):',//,
     2        6X,'sequential read',2X,F10.4,/,
     3        6X,'DA read',10X,F10.4)
 2600 FORMAT(//,'    Thresholds for 2-matrix elements:',//,
     1        6X,'during first  half-transformation',1P,D18.8,/,
     2        6X,'during second half-transformation',1P,D18.8)
 2700 FORMAT(//,'    Print level is',I5)
 2800 FORMAT(//,'    Symmetry data:',//,
     1        6X,'Order of group',I4,//,
     2        10X,'Multiplication Table',/)
 2850 FORMAT(15X,8I4)
 2900 FORMAT(//,'    Orbital data:',//,4X,'Irrep',15X,8I5)
 2910 FORMAT(4X,'Inactive MOs',8X,8I5)
 2920 FORMAT(4X,'Active MOs',10X,8I5)
 2930 FORMAT(4X,'Symmetry orbitals',3X,8I5)
 3000 FORMAT(//,'    Active orbital symmetries:',//,
     1        10X,'Orbital',4x,'Symmetry',//,
     2        (10X,I5,6X,I5))
 3100 FORMAT(//,'    Symmetry orbitals:',//,
     1        10X,'Orbital',4X,'Symmetry',//,
     2        (10X,I5,6X,I5))
 3200 FORMAT(//,'    Active orbital symmetry offset vector:',/,4X,8I5)
 3300 FORMAT(//,'    All orbitals   symmetry offset vector:',/,4X,8I5)
      END
C  /* Deck pt1tra */
      SUBROUTINE PT1TRA(CMO,PV,P2,P4V,IBUF,BUF,IT,WORK,LWORK,ANTI,
     &                  PANTI,DIA2SO,ZFS2EL)
C
C....    FIRST HALF-TRANSFORMATION OF AN MCSCF 2-MATRIX.
C....    THE LAST TWO 2-MATRIX INDICES, V AND X, ARE TRANSFORMED TO
C....    SYMMETRY ORBITAL INDICES LAMBDA AND SIGMA. THE 2-MATRIX PV
C....    IS STORED AS IN THE CASSCF PROGRAM, AS A LOWER LOWER TRIANGLE
C....    SUPERMATRIX WITH ZEROES IN PLACE. NO SYMMETRY BLOCKING IS
C....    USED.
C
C....    ******      VERSION 1.0      RELEASE 820423      ******
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
      INTEGER*8  L_FIELD, FAC_1, FAC_2, FAC_3, I8_2
      PARAMETER (L_FIELD = 16, I8_2 = 2,
     &           FAC_1 = I8_2**(3*L_FIELD),
     &           FAC_2 = I8_2**(2*L_FIELD),
     &           FAC_3 = I8_2**   L_FIELD     )
      PARAMETER (DHALF=0.5D0, D0=0.0D0, D1=1.0D0)
      PARAMETER (D2 = 2.D0)
C
      LOGICAL   ANTI, PANTI, DIA2SO, ZFS2EL, FOUND
      INTEGER   T, U, V, X, TU, VX, TUVX, SIGMA, SIGMB, ALPHA, ALPHAP,
     *          IT(*)
      INTEGER*8 IBUF(*), IND
      REAL*8    PV(*), CMO(*), BUF(*), P2(*), P4V(*)
C
#include "priunit.h"
#include "inftap.h"
#include "abainf.h"
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NBAS(8), NASHT, NBAST,
     1                IASH(8), IBAS(8), IASYM(MXCORB), IOSYM(MXCORB),
     2                LSMAX, MOFF(MXCORB), NCMOT, MXBAS1, MXBAS2, NNASHX
      COMMON /PTRINF/ THR1, THR2, DAWEIG, SQWEIG, IPRINT
      COMMON /PTRSYM/ NSYM, MUL(8,8), NENT(8,8), ISYOFF(8,8)
C
      CALL QENTER('PT1TRA')
C
C     Initialize counters and rewind half-transformed 2-matrix
C     output file.
C
      DO, I = 1,NBAST
        IT(I) = I*(I-1)/2
      END DO
      L1BUF = 2*MX1BUF + 1
      NBUF  = 0
      NPTOT = 0
      CALL GPOPEN(LUSCR,'PTRSCR','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      REWIND LUSCR
C
C     Read PMAT and MOs.  This form will handle symmetry as is.
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL RD_SIRIFC('CMO',FOUND,CMO,WORK,LWORK)
      IF (.NOT.FOUND)
     &   CALL QUIT('PT1TRA error: CMO not found on SIRIFC')
      IF (ANTI.OR.DIA2SO.OR.ZFS2EL) THEN
         CALL MAKPVM(PV,WORK,LWORK,ANTI,PANTI,DIA2SO,ZFS2EL,IPRINT)
      ELSE
         CALL RD_SIRIFC('PV',FOUND,PV,WORK,LWORK)
         IF (.NOT.FOUND)
     &        CALL QUIT('PT1TRA error: PV not found on SIRIFC')
      END IF
      IF (IPRINT .LT. 10) GOTO 800
      WRITE(LUPRI,2100)
      TU = 0
      DO 810 T = 1,NASHT
        DO 820 U = 1,T
          TU = TU + 1
          VX = 0
          DO 830 V = 1,T
            LV = V
            IF (V .EQ. T) LV = U
            DO 840 X = 1,LV
              VX = VX + 1
              TUVX = (VX - 1)*NNASHX + TU
              WRITE(LUPRI,2200) T, U, V, X, PV(TUVX)
  840       CONTINUE
  830     CONTINUE
  820   CONTINUE
  810 CONTINUE
  800 CONTINUE
      IEND = 0
      IAC = 0
      IF (IPRINT .GE. 10) WRITE(LUPRI,2600)
      DO 710 I = 1,NSYM
        NORBI = NORB(I)
        IF (NORBI .EQ. 0) GOTO 710
        NASHI = NASH(I)
        NISHI = NISH(I)
        NOCCI = NISHI + NASHI
        NBASI = NBAS(I)
        DO 720 J = 1,NORBI
          IBEG = IEND + 1
          IEND = IEND + NBASI
          IF (J .GT. NOCCI .OR. J .LE. NISHI) GOTO 720
          IAC = IAC + 1
          MOFF(IAC) = IBEG - 1
          IF (IPRINT .GE. 10)
     1        WRITE(LUPRI,2500) I,J,IAC,MOFF(IAC),
     2                           (CMO(K),K = IBEG,IEND)
  720   CONTINUE
  710 CONTINUE
C
C....    OUTER DRIVER LOOPS OVER FIRST TWO INDICES (ACTIVE): THESE
C....    ARE FULL RANGES WITH SYMMETRY DETERMINATION INSIDE THE LOOP
C
      DO 10 T = 1,NASHT
        ISYMT = IASYM(T)
        ITT = IT(T)
        DO 20 U = 1,T
          ISYMU = IASYM(U)
          TU = ITT + U
          ISYMTU = MUL(ISYMT,ISYMU)
          DO 40 ALPHA = 1,NSYM
            ALPHAP = MUL(ALPHA,ISYMTU)
            IF (ALPHAP .LT. ALPHA) GOTO 40
            NAV = NASH(ALPHA)
            IF (NAV .EQ. 0) GOTO 40
            NAVB = IASH(ALPHA)
            NAVE = NAVB + NAV
            NAVB = NAVB + 1
            NBASL = NBAS(ALPHA)
            IBASL = IBAS(ALPHA)
            NAX = NASH(ALPHAP)
            IF (NAX .EQ. 0) GOTO 40
            NAXB = IASH(ALPHAP)
            NAXE = NAXB + NAX
            NAXB = NAXB + 1
            NBASP = NBAS(ALPHAP)
            IBASS = IBAS(ALPHAP)
            NCLR = NBASL*NBASP
            IF (ALPHA .EQ. ALPHAP) NCLR = IT(NBASL) + NBASL
            DO 45 LS = 1,NCLR
              P2(LS) = D0
   45       CONTINUE
C
C....    THIRD INDEX DRIVER LOOP
C
            DO 50 V = NAVB,NAVE
C
C....    CLEAR STORAGE FOR QUARTER-TRANSFORMED INTEGRALS
C
              DO 55 SIGMA = 1,NBASP
                P4V(SIGMA) = D0
   55         CONTINUE
C
C....    MO ORIGIN OFFSET
C
              MOFFV = MOFF(V)
              IVT = IT(V)
C
C....    FOURTH INDEX DRIVER LOOP
C
              DO 60 X = NAXB, NAXE
                MOFFX = MOFF(X)
C
C....    PICK UP DENSITY MATRIX ELEMENT
C
                VX = IVT + X
                IF (X .GT. V) VX = IT(X) + V
                TUVX = (VX - 1)*NNASHX + TU
                PVAL = PV(TUVX)
                IF (ANTI .AND. X .GT. V) PVAL = - PVAL
                IF (ABS(PVAL) .LT. THR1) GOTO 60
C
C....    INNER LOOP FOR FIRST QUARTER-TRANSFORMATION
C
                DO 70 SIGMA = 1,NBASP
                  P4V(SIGMA) = P4V(SIGMA) + CMO(MOFFX+SIGMA)*PVAL
   70           CONTINUE
   60         CONTINUE
C
C....    SECOND QUARTER-TRANSFORMATION
C
              LS = 0
              DO 80 SIGMA = 1,NBASP
                LSTART = 1
                P4VS = P4V(SIGMA)
                IF (ALPHA .EQ. ALPHAP)LSTART = SIGMA
                IF (ABS(P4VS) .GE. THR1) GOTO 85
                LS = LS + NBASL - LSTART + 1
                GOTO 80
   85           CONTINUE
                DO 90 LAMBDA = LSTART, NBASL
                  LS = LS + 1
                  P2(LS) = P2(LS) + P4VS*CMO(MOFFV+LAMBDA)
   90           CONTINUE
   80         CONTINUE
   50       CONTINUE
C
C....    P2(L,S) CONTAINS PV(T,U,LA,SA') FOR FIXED T, U, L IN A AND
C....    S IN A'. WRITE THESE ELEMENTS OUT
C
            LS = 0
            IF (IPRINT .GE. 40) WRITE(LUPRI,2400)
            DO 100 SIGMA = 1,NBASP
              LSTART = 1
              IF (ALPHA .EQ. ALPHAP) LSTART = SIGMA
              DO 110 LAMBDA = LSTART,NBASL
                LS = LS + 1
                LAMBDB = LAMBDA + IBASL
                SIGMB  = SIGMA  + IBASS
                PVAL = P2(LS)
                IF (ABS(PVAL) .LT. THR1) GOTO 110
                IF (LAMBDB .GE. SIGMB) GOTO 105
                ILY = LAMBDB
                LAMBDB = SIGMB
                SIGMB = ILY
                IF (ANTI) PVAL = - PVAL
  105           CONTINUE
                IND = T*FAC_1 + U*FAC_2 + LAMBDB*FAC_3 + SIGMB
                NBUF = NBUF + 1
                IBUF(NBUF) = IND
                BUF(NBUF) = PVAL
                IF (IPRINT .GE. 40)
     1             WRITE (LUPRI,2300) T, U, ALPHA, LAMBDA, ALPHAP,
     2                                 SIGMA, PVAL
                IF (NBUF .GE. MX1BUF) THEN
                   IBUF(L1BUF) = NBUF
                   CALL WRITT(LUSCR,L1BUF,IBUF)
                   NPTOT = NPTOT + NBUF
                   NBUF = 0
                END IF
  110         CONTINUE
  100       CONTINUE
C
C....    END SYMMETRY DRIVER LOOP
C
   40     CONTINUE
C
C....    END T, U LOOPS
C
   20   CONTINUE
   10 CONTINUE
C
C....    EMPTY BUFFER IF NECESSARY AND WRITE END OF FILE
C
      IF (NBUF .GT. 0) THEN
         IBUF(L1BUF) = NBUF
         CALL WRITT(LUSCR,L1BUF,IBUF)
         NPTOT = NPTOT + NBUF
      END IF
      IF (IPRINT .GT. 2) WRITE(LUPRI,2000) NPTOT
      CALL QEXIT('PT1TRA')
      RETURN
 2000 FORMAT(//,'    Number of half-transformed matrix elements is',I12)
 2100 FORMAT(//,'    MCSCF 2-matrix',/)
 2200 FORMAT(4X,4I5,F12.8)
 2300 FORMAT(4X,I5,4X,I5,4X,I2,4X,I5,4X,I2,4X,I5,4X,F12.8)
 2400 FORMAT(8X,'T',8X,'U',5X,'A',8X,'L',4X,'A`',8X,'S',7X,'element')
 2500 FORMAT(/,4X,'symmetry',I2,'   active orbital',I5,
     1       '   active index',I5,'  offset',I5,//,
     2       (4X,5F12.8))
 2600 FORMAT(//,'    Active molecular orbitals:')
      END
C  /* Deck pt1sor */
      SUBROUTINE PT1SOR(ISORT,SORT,IBUF,BUF,IT)
C
C....    SINGLE PASS BIN SORT OF PTULS ON LS VALUE. BLOCKING IS
C....    DETERMINED FOR SYMMETRY-BLOCKED STORAGE OF TU FOR EACH LS.
C....    THE ACTUAL BLOCK LENGTHS ARE THUS DETERMINED BY THE NASH VALUES
C....    THIS IS DONE TO FACILITATE LATER EXPANSION TO CI 2-MATRICES.
C
C
C....    ******      VERSION 1.0      RELEASE 820423      ******
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
      LOGICAL OLDDX
C
      INTEGER   IT(*), ALPHA, BETA, GAMMA, SIGMA
      INTEGER*8 ISORT(*), IBUF(*)
      REAL*8     SORT(*),  BUF(*)
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NBAS(8), NASHT, NBAST,
     1                IASH(8), IBAS(8), IASYM(MXCORB), IOSYM(MXCORB),
     2                LSMAX, MOFF(MXCORB), NCMOT, MXBAS1, MXBAS2, NNASHX
      COMMON /PTRINF/ THR1, THR2, DAWEIG, SQWEIG, IPRINT
      COMMON /PTRSYM/ NSYM, MUL(8,8), NENT(8,8), ISYOFF(8,8)

      INTEGER*8  L_FIELD, IBIT_FIELD
      PARAMETER (L_FIELD = 16, IBIT_FIELD = (2**L_FIELD-1))
      INTEGER*8 I, J

      INTEGER*8 LABV, LABEL
      INTEGER   LOCV
      GET_FIELD(LABV,LOCV)
     &   = IAND(ISHFT(LABV,-L_FIELD*LOCV),IBIT_FIELD)

      CALL QENTER('PT1SOR')
C
C....    COMPUTE ARRAY (EXPANDED TO A SQUARE) OF BLOCK LENGTHS
C
      NENMAX = 0
      DO 10 ALPHA = 1,NSYM
        DO 20 BETA = 1,ALPHA
          MULAB = MUL(ALPHA,BETA)
          NTOT = 0
          DO 30 GAMMA = 1,NSYM
          NAG = NASH(GAMMA)
          IF (NAG .EQ. 0) GOTO 30
          MULABG = MUL(MULAB,GAMMA)
C
C....    THE FOURTH SYMMETRY MUST BE MULABG, OR THE BLOCK VANISHES BY
C....    SYMMETRY
C
          ISYOFF(MULAB,GAMMA) = NTOT
          NAD = NASH(MULABG)
          IF (NAD .EQ. 0) GOTO 30
          NBLK = NAG*NAD
          IF (MULAB .EQ. 1) NBLK = IT(NAG) + NAG
          NTOT = NTOT + NBLK
   30   CONTINUE
        NENT(ALPHA,BETA) = NTOT
        NENT(BETA,ALPHA) = NTOT
        IF (NTOT .GT. NENMAX) NENMAX = NTOT
   20 CONTINUE
   10 CONTINUE
C
C....    CHECK THAT THE LARGEST ENTRY IS NOT TOO LARGE.
C
      IF (NENMAX .GT. MEMT)
     1    CALL QUIT('INSUFFICIENT TRANSFORMATION CORE - PT1SOR')
C
C....    DETERMINE SORT PARAMETERS
C
      NBLOCK = MEMT/NENMAX
      IF (NBLOCK .GT. MAXADR) NBLOCK = MAXADR
C
C....    NBLOCK IS NUMBER OF BLOCKS TO HOLD SIMULTANEOUSLY
C
      NTOTP = (NBAST*(NBAST+1))/2
      IF (NBLOCK .GT. NTOTP) NBLOCK = NTOTP
C
C....    NTOTP IS TOTAL NUMBER OF LS VALUES
C
      NCHAIN = (NTOTP + NBLOCK - 1)/NBLOCK
C
C....    NCHAIN IS THE NUMBER OF CHAINS
C
      IF (NCHAIN .GT. MAXCHN) CALL QUIT('Too many chains - PT1SOR')
      MXABUF = (MEMS/NCHAIN - 2)/2
      MXABUF = MIN(MXABUF, (LDAMAX - 2)/2)
      LABUF  = 2*MXABUF + 2
      LABUF1 = LABUF - 1
      CALL GPOPEN(LUDA,'ABAPTRAN1.DA','NEW','DIRECT',' ',
     &   IRAT*LABUF,OLDDX)
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,2100) MEMS, MEMT, NBLOCK, NCHAIN, LABUF
         WRITE (LUPRI,'(/4X,A)') 'ISYOFF array:'
         DO, I = 1,NSYM
           WRITE(LUPRI,2850) (ISYOFF(I,J),J = 1,NSYM)
         END DO
         WRITE (LUPRI,'(/4X,A)') 'NENT array:'
         DO, I = 1,NSYM
           WRITE(LUPRI,2850) (NENT(I,J),J = 1,NSYM)
         END DO
         CALL FLSHFO(LUPRI)
      END IF
 2850 FORMAT(15X,8I4)
C
C....    NOW SET UP BUFFER ADDRESSES, AND PRESET DISK ADDRESS AND
C....    BUFFER LENGTH LOCATIONS IN THE BUFFERS
C
      NOFF = 0
      DO 100 ICH = 1,NCHAIN
        IADR(ICH) = NOFF
        NOFF = NOFF + LABUF
        ISORT(NOFF) = -1
        ISORT(NOFF-1) = 0
  100 CONTINUE
C
C....    LOOP OVER PTULS
C
      REWIND LUSCR
      IDISK = 1
C     IDISK = 0
 1000 CONTINUE
      CALL READI8(LUSCR,L1BUF,IBUF)
      LENGTH = IBUF(L1BUF)
      IF (LENGTH .EQ. 0) GOTO 1000
      IF (LENGTH .LT. 0) GOTO 2000
      DO 200 I = 1,LENGTH
        LABEL = IBUF(I)
        LAMBDA = GET_FIELD(LABEL,1)
        SIGMA  = GET_FIELD(LABEL,0)
        LS = IT(LAMBDA) + SIGMA
C
C....    FIND APPROPRIATE BIN
C
        LLS = (LS+NBLOCK-1)/NBLOCK
        IOFF = IADR(LLS)
        NBUF = ISORT(IOFF+LABUF1) + 1
        ISORT(IOFF+NBUF) = LABEL
        SORT((IOFF+MXABUF) + NBUF) = BUF(I)
C
        ISORT(IOFF+LABUF1) = NBUF
        IF (NBUF .LT. MXABUF) GOTO 200
C
C....    THIS SORT BIN IS FULL
C
        IDISKO = IDISK
        CALL WRITDX(LUDA,IDISK,IRAT*LABUF,ISORT(IOFF+1))
        IDISK = IDISK + 1
        IF (IPRINT .GE. 90) WRITE(LUPRI,2200) IDISKO, NBUF,
     1     (ISORT(IOFF+K),SORT((IOFF+MXABUF)+K), K = 1,NBUF)
        ISORT(IOFF+LABUF1) = 0
        ISORT(IOFF+LABUF) = IDISKO
  200 CONTINUE
      GOTO 1000
 2000 CONTINUE
C
C....    NOW EMPTY BUFFERS, AND WRITE LAST ADDRESSES
C
      DO 230 I = 1,NCHAIN
        IOFF   = IADR(I)
        IDISKO = IDISK
        CALL WRITDX(LUDA,IDISK,IRAT*LABUF,ISORT(IOFF+1))
        IDISK = IDISK + 1
        IF (IPRINT .GE. 90) WRITE(LUPRI,2200) IDISKO, NBUF,
     2     (ISORT(IOFF+K),SORT((IOFF+MXABUF)+K), K = 1,NBUF)
        LASTAD(I) = IDISKO
  230 CONTINUE
      CALL QEXIT('PT1SOR')
      RETURN
 2100 FORMAT(//,'    Sort Parameters:',//,
     *       6X,'Sorting memory            ',I8,/,
     *       6X,'Transformation memory     ',I8,/,
     *       6X,'Number of blocks per chain',I8,/,
     *       6X,'Number of chains          ',I8,/,
     *       6X,'Length of DA buffer       ',I8)
 2200 FORMAT(/'    Bin contents at disk address ',I12,'   length',I6,
     *       /(2(1X,I15,F16.8,5X)))
      END
C  /* Deck pt2tra */
      SUBROUTINE PT2TRA(CMO,P4T,P2,PHALF,IBIN,BIN,IBUF,BUF,IT,
     &                  ANTI,PANTI)
C
C....    Second half-transformation of 2-matrix.  Elements PHALF TULS
C....    are transformed to PHALF MNLS.
C
C
C....    ******      VERSION 1.0      RELEASE 820423      ******
C
#include "implicit.h"
#include "iratdef.h"
#include "dummy.h"
C
      PARAMETER (D0 = 0.0D0)
#include "maxorb.h"
#include "mxcent.h"
C
      INTEGER   SIGMA, T, U, TA, UA, TM, UM, ALPHA, BETA, GAMMA,
     *          ALPHAP, SIGMB, TUA, IT(*)
      INTEGER*8 IBUF(*), IBIN(*)
      LOGICAL   LTRIAN, LALBE, ANTI, PANTI, FOUND
      REAL*8    CMO(*), PHALF(*), P4T(*), P2(*), BUF(*), BIN(*)
C
#include "priunit.h"
#include "abainf.h"
#include "inftap.h"
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NBAS(8), NASHT, NBAST,
     1                IASH(8), IBAS(8), IASYM(MXCORB), IOSYM(MXCORB),
     2                LSMAX, MOFF(MXCORB), NCMOT, MXBAS1, MXBAS2, NNASHX
      COMMON /PTRINF/ THR1, THR2, DAWEIG, SQWEIG, IPRINT
      COMMON /PTRSYM/ NSYM, MUL(8,8), NENT(8,8), ISYOFF(8,8)
      COMMON /PTRTOT/ NPTOT

      INTEGER*8  L_FIELD, IBIT_FIELD, FAC_1, FAC_2, FAC_3, I8_2
      PARAMETER (L_FIELD = 16, IBIT_FIELD = (2**L_FIELD-1), I8_2 = 2)
      PARAMETER (FAC_1 = I8_2**(3*L_FIELD),
     &           FAC_2 = I8_2**(2*L_FIELD),
     &           FAC_3 = I8_2**   L_FIELD     )
      INTEGER*8 I, J

C
      INTEGER*8 LABV, LABEL, MUX, NUX, LSPACK
      INTEGER   LOCV
      GET_FIELD(LABV,LOCV)
     &   = IAND(ISHFT(LABV,-L_FIELD*LOCV),IBIT_FIELD)

      UNPAKK(IJ) = SQRT(2.0D0*IJ + 0.25D0) + 0.4999D0
C
      CALL QENTER('PT2TRA')
C
C....     Retrieve MO coefficients
C
      CALL RD_SIRIFC('CMO',FOUND,CMO,BUF,MX2BUF)
      IF (.NOT.FOUND) CALL QUIT('PT2TRA error: CMO not found on SIRIFC')
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,2600)
         IAC = 0
         IEND = 0
         DO 710 I = 1,NSYM
           NORBI = NORB(I)
           IF (NORBI .EQ. 0) GOTO 710
           NASHI = NASH(I)
           NISHI = NISH(I)
           NOCCI = NISHI + NASHI
           NBASI = NBAS(I)
           DO 720 J = 1,NORBI
             IBEG = IEND + 1
             IEND = IEND + NBASI
             IF (J .GT. NOCCI .OR. J .LE. NISHI) GOTO 720
             IAC = IAC + 1
             WRITE(LUPRI,2500) I,J,IAC,MOFF(IAC),
     &                          (CMO(K),K = IBEG,IEND)
  720      CONTINUE
  710    CONTINUE
      END IF
C
      L2BUF = 2*MX2BUF + 1
      IBUF(L2BUF) = 0
      NPTOT = 0
C
C....    DRIVER LOOP BASED ON THE NUMBER OF BIN SORT CHAINS
C
      CALL GPOPEN(LUPSO,'PTRPSO','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      REWIND LUPSO
      NTOTP = (NBAST*(NBAST+1))/2
      LSMAX = 0
      DO 100 ICH = 1,NCHAIN
C
C     Kludge to avoid excessive real memory use on virtual memory
C     systems when running small cases.
C
      NCLTOP = MIN(MEMT,NTOTP**2)
      DO 105 I = 1,NCLTOP
        PHALF(I) = D0
  105 CONTINUE
C
C....    LS INDEX RANGE
C
      LSMIN = LSMAX + 1
      LSMAX = LSMAX + NBLOCK
      IF (LSMAX .GT. NTOTP) LSMAX = NTOTP
C
C....    COMPUTE ARRAY STARTING ADDRESSES FOR TU ELEMENTS
C
      NBEG = 0
      DO 110 LS = LSMIN,LSMAX
         IND    = LS - LSMIN + 1
         LAMBDA = UNPAKK(LS)
         SIGMA  = LS - IT(LAMBDA)
         ALPHA  = IOSYM(LAMBDA)
         BETA   = IOSYM(SIGMA)
         NTOT   = NENT(ALPHA,BETA)
         IADR(IND) = -1
         IF (NASH(ALPHA) .EQ. 0) GOTO 110
         IF (NASH(BETA)  .EQ. 0) GOTO 110
         IADR(IND) = NBEG
         NBEG   = NTOT + NBEG
  110 CONTINUE
      LSMIN1 = LSMIN - 1
      IF (IPRINT .GE. 30) THEN
         WRITE(LUPRI,2700)
         DO 116 LS = LSMIN,LSMAX
            LSP = LS - LSMIN1
            WRITE(LUPRI,2800) LS, IADR(LSP)
  116    CONTINUE
      END IF
C
C....    NOW READ BACK THE BINS IN THIS CHAIN
C
C....    SET DISK ADDRESS
C
      IDISK  = LASTAD(ICH)
  120 CONTINUE
      CALL READDX(LUDA,IDISK,IRAT*LABUF,IBIN)
      IDISK  = IBIN(LABUF)
      LENGTH = IBIN(LABUF-1)
      IF (LENGTH .EQ. 0) THEN
         IF (IDISK .EQ. -1) GO TO 131
         GO TO 120
      END IF
C
C....    LOOP OVER ELEMENTS IN THIS BIN AND DISTRIBUTE THEM TO THE
C....    DESIRED LOCATIONS
C
      DO 130 I = 1,LENGTH
        PVAL   = BIN(I)
        LABEL  = IBIN(I)
        SIGMA  = GET_FIELD(LABEL,0)
        LAMBDA = GET_FIELD(LABEL,1)
        U      = GET_FIELD(LABEL,2)
        T      = GET_FIELD(LABEL,3)
        LS     = IT(LAMBDA) + SIGMA
        LSP    = LS - LSMIN1
        IOFF   = IADR(LSP)
C
C....    NOW DETERMINE SUB-ADDRESS FOR TU BASED ON TU AND LS SYMMETRY
C
        ALPHA  = IOSYM(LAMBDA)
        BETA   = IOSYM(SIGMA)
        ISYMLS = MUL(ALPHA,BETA)
        ISYMT  = IASYM(T)
        ISYMU  = IASYM(U)
        TA     = T - IASH(ISYMT)
        UA     = U - IASH(ISYMU)
        IF (ISYMT .EQ. ISYMU) THEN
           TUA = IT(TA) + UA
        ELSE
           TUA = NASH(ISYMU)*(TA-1) + UA
        END IF
        TUA = TUA + ISYOFF(ISYMLS,ISYMT)
        PHALF(IOFF+TUA) = PVAL
  130 CONTINUE
      IF (IDISK .NE. -1) GOTO 120
  131 CONTINUE
C
C....    PHALF BLOCKS NOW IN CORE. BEGIN TRANSFORMATION
C
      DO 160 LS = LSMIN,LSMAX
      LSP = LS - LSMIN1
      IADLSP = IADR(LSP)
      IF (IADLSP .LT. 0) GOTO 160
      LAMBDA = UNPAKK(LS)
      SIGMA = LS - IT(LAMBDA)
      LSPACK = LAMBDA*FAC_3 + SIGMA
      BETA = IOSYM(LAMBDA)
      GAMMA = IOSYM(SIGMA)
      LAMBDB = LAMBDA - IBAS(BETA)
      SIGMB = SIGMA - IBAS(GAMMA)
      ISYMLS = MUL(BETA,GAMMA)
      LTRIAN = ISYMLS .EQ. 1
        DO 170 ALPHA = 1,BETA
          LALBE = ALPHA .EQ. BETA
          ALPHAP = MUL(ISYMLS,ALPHA)
          IF (ALPHAP .GT. ALPHA) GOTO 170
          ISOFF = ISYOFF(ISYMLS,ALPHA) + IADLSP
          NAT = NASH(ALPHA)
          IF (NAT .EQ. 0) GOTO 170
          NAU = NASH(ALPHAP)
          IF (NAU .EQ. 0) GOTO 170
          NATB  = IASH(ALPHA)
          NATE  = NATB + NAT
          NATB1 = NATB
          NATB  = NATB + 1
          NBASM = NBAS(ALPHA)
          IBASM = IBAS(ALPHA)
          NAUB  = IASH(ALPHAP)
          NAUE  = NAUB + NAU
          NAUB1 = NAUB
          NAUB  = NAUB + 1
          NBASN = NBAS(ALPHAP)
          IBASN = IBAS(ALPHAP)
C
C....    CLEAR TRANSFORMED VECTOR
C
          NCLR = NBASM*NBASN
          IF (LTRIAN) NCLR = IT(NBASM) + NBASM
          DO 175 MN = 1,NCLR
            P2(MN) = D0
  175     CONTINUE
C
C....    FIRST ACTIVE INDEX DRIVER LOOP
C
          DO 180 T = NATB,NATE
            TM = T - NATB1
C
C....    CLEAR QUARTER-TRANSFORMED VECTOR
C
            DO 185 NU = 1,NBASN
              P4T(NU) = D0
  185       CONTINUE
C
C....    PHALF BLOCK AND MO OFFSETS
C
            MOFFT = MOFF(T)
            IOMT = NAU*(TM-1)
            IF (LTRIAN) IOMT = IT(TM)
            IOMT = IOMT + ISOFF
C
C....    SECOND INDEX DRIVER LOOP
C
            DO 190 U = NAUB,NAUE
              UM = U - NAUB1
              MOFFU = MOFF(U)
              IF (LTRIAN) GOTO 200
              INDP = IOMT + UM
              PVAL = PHALF(INDP)
              GOTO 210
  200         CONTINUE
              INDP = IOMT + UM
              IF (UM .GT. TM) INDP = ISOFF + IT(UM) + TM
              PVAL = PHALF(INDP)
              IF (ANTI .AND. UM .GT. TM) PVAL = - PVAL
  210         CONTINUE
              IF (ABS(PVAL) .LT. THR2) GOTO 190
C
C....    TRANSFORM SECOND INDEX
C
              DO 220 NU = 1,NBASN
                P4T(NU) = P4T(NU) + CMO(MOFFU+NU)*PVAL
  220         CONTINUE
  190       CONTINUE
C
C....    NOW COMPLETE PROCESS BY FINAL QUARTER-TRANSFORMATION (FIRST
C....    INDEX)
C
            MN = 0
            DO 230 NU = 1,NBASN
              MSTART = 1
              IF (LALBE) GOTO 231
              IF (.NOT. LTRIAN) GOTO 232
              MSTART = NU
              GOTO 232
  231         CONTINUE
              IF (LTRIAN) GOTO 233
              MSTART = LAMBDB
              IF (NU .LT. SIGMB) MSTART = LAMBDB + 1
              IF (MSTART .GT. NBASM) GOTO 230
              GOTO 232
  233         CONTINUE
              MSTART = MAX0(NU,LAMBDB)
              IF (NU .LT. SIGMB) MSTART = LAMBDB + 1
              IF (MSTART .GT. NBASM) GOTO 230
  232         CONTINUE
              P4TN = P4T(NU)
              IF (ABS(P4TN) .GE. THR2) GOTO 235
              MN = MN + NBASM - MSTART + 1
              GOTO 230
  235         CONTINUE
              DO 240 MU = MSTART,NBASM
                MN = MN + 1
                P2(MN) = P2(MN) + P4TN*CMO(MOFFT+MU)
  240         CONTINUE
  230       CONTINUE
  180     CONTINUE
C
C....    THIS GIVES PHALF M A  N A'  L S  FOR FIXED A A' L S. WRITE
C....    THESE OUT
C
          MN = 0
          IF (IPRINT .GE. 70) WRITE(LUPRI,2400)
          DO 250 NU = 1,NBASN
              MSTART = 1
              IF (LALBE) GOTO 251
              IF (.NOT. LTRIAN) GOTO 252
              MSTART = NU
              GOTO 252
  251         CONTINUE
              IF (LTRIAN) GOTO 253
              MSTART = LAMBDB
              IF (NU .LT. SIGMB) MSTART = LAMBDB + 1
              IF (MSTART .GT. NBASM) GOTO 250
              GOTO 252
  253         CONTINUE
              MSTART = MAX0(NU,LAMBDB)
              IF (NU .LT. SIGMB) MSTART = LAMBDB + 1
              IF (MSTART .GT. NBASM) GOTO 250
  252         CONTINUE
            DO 260 MU = MSTART,NBASM
              MN = MN + 1
              PVAL = P2(MN)
            IF (ABS(PVAL) .LT. THR2) GOTO 260
C
C....    CONSTRUCT LABEL FROM FULL INDICES
C
              MUX = MU + IBASM
              NUX = NU + IBASN
              IF (MUX .LT. NUX) THEN
                 MUX = NU + IBASN
                 NUX = MU + IBASM
                 IF (ANTI) PVAL = - PVAL
              END IF
              IF ((MUX.LT.LAMBDA .OR. MUX.EQ.LAMBDA.AND.NUX.LT.SIGMA)
     &             .AND. PANTI) PVAL=-PVAL
              LABEL = MUX*FAC_1 + NUX*FAC_2 + LSPACK
              NOBUF = IBUF(L2BUF) + 1
              BUF( NOBUF) = PVAL
              IBUF(NOBUF) = LABEL
              IBUF(L2BUF) = NOBUF
              IF (IPRINT .GE. 70) WRITE(LUPRI,2300)
     *           LAMBDA, SIGMA, ALPHA, MU, ALPHAP, NU, PVAL
              IF (NOBUF .GE. MX2BUF) THEN
                 CALL WRITT(LUPSO,L2BUF,IBUF)
                 IBUF(L2BUF) = 0
                 NPTOT = NPTOT + NOBUF
              END IF
  260       CONTINUE
  250     CONTINUE
C
C        END LOOP OVER SYMMETRIES FOR THIS LS PAIR
C
  170   CONTINUE
C
C....    END LS LOOP FOR THIS CHAIN
C
  160 CONTINUE
C
C....    END LOOP OVER CHAINS
C
  100 CONTINUE
C
C....    EMPTY BUFFER IF REQUIRED
C
      NOBUF = IBUF(L2BUF)
      IF (NOBUF .NE. 0 .OR. NPTOT.EQ.0) CALL WRITT(LUPSO,L2BUF,IBUF)
      NPTOT = NPTOT + NOBUF
      REWIND LUPSO
      IF (IPRINT .GT. 2) WRITE(LUPRI,2000) NPTOT
      CALL QEXIT('PT2TRA')
      RETURN
 2000 FORMAT(//,' Number of two-electron density',
     *          ' elements transformed:',I8)
 2300 FORMAT(4X,I5,4X,I5,4X,I2,4X,I5,4X,I2,4X,I5,4X,F12.8)
 2400 FORMAT(8X,'L',8X,'S',5X,'A',8X,'M',4X,'A`',8X,'N',7X,'element')
 2500 FORMAT(/,4X,'symmetry',I2,'   active orbital',I5,
     1       '   active index',I5,'  offset',I5,//,
     2       (4X,5F12.8))
 2600 FORMAT(//,'    Active molecular orbitals:')
 2700 FORMAT(//,'    Index vector for pair offsets')
 2800 FORMAT(6X,I6,6X,I8)
      END
C  /* Deck psorg */
      SUBROUTINE PSORG(WORK,IWORK,LWORK,NCONTS,IPRINT,ANTI,PANTI)
C
C     Organize sorting of symmetry-adapted second-order reduced
C     density matrix to match calculation of distinct two-electron
C     integral derivatives.
C
C     ANTI : P(ij,kl) = -P(ji,kl) = -P(ij,lk)
C     PANTI: P(ij,kl) = -P(kl,ij)
C
C     ******      VERSION 1.0   RELEASE 880925      ******
C     ******      VERSION 2.0   RELEASE 140731      ******
C
#include "implicit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "aovec.h"
C
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
C
C     If MAXBUF is changed be sure to make the corresponding change
C     of MX2BUF in PTRAN (used in PT2TRA) and in DISINT.
C
      PARAMETER ( MINMEM = 50000,
     *            MAXLOC = 4000,
     *            MAXPAS = 12,    IBUFMX = 600,
     *            IBUFMN = 200,   MXMAX = 2000 )
C
      LOGICAL   ANTI, PANTI, OLDDX
      INTEGER   P, Q, R, S, PQST1, PQST, PQFIN, PQ, RS,
     *          RSMIN1, RSMIN, RSMAX
      INTEGER*8 IWORK(LWORK)
      REAL*8     WORK(LWORK)
      INTEGER   NCONTS(MXSHEL,MXAOVC,2)
#include "ccom.h"
#include "blocks.h"
#include "symmet.h"
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      COMMON /PTRORB/ NISH(8), NASH(8), NORB(8), NBAS(8), NASHT, NBAST,
     1                IASH(8), IBAS(8), IASYM(MXCORB), IOSYM(MXCORB),
     2                LSMAX, MOFF(MXCORB), NCMOT, MXBAS1, MXBAS2, NNASHX
      COMMON /PTRTOT/ NPTOTR
C

      CALL QENTER('PSORG')
      MAXBUF = MX2BUF ! fetch from ptrbuf.h
      IF (IPRINT .GT. 2) WRITE (LUPRI, 1000)
CHJ-860802:  find MEMS = available space for SORT(*), ISORT(*)
      MEMS  = LWORK / 2
      IF (MEMS .LT. MINMEM)
     *   CALL STOPIT('PSORG',' ',2*MINMEM,2*MEMS)
      LBUF  = 2*MAXBUF +  1
C
C           *************************************
C           **  COMPUTE BLOCK LENGTHS FOR SORT **
C           *************************************
C
      MEMMAX = (MEMS-LBUF-2*MAXADR-3*MAXLOC)
      NBLOCK = 0
      NPS = 0
      NPSMAX = 0
      MEMTOT = 0
      NSLICE = 1
      NUMP = 0
      DO 10 P = 1,MAXSHL
         KHKTP = KHKTSH(P)
         MULP = MULT(ISTBSH(P))
         NSETP = NSETSH(P,1)
         NRCP = 0
         DO 15 IP = 1,NSETP
            NRCPI = NCONTS(P,IP,1)
            IF (NRCPI .LT. 0) NRCPI = 1
            NRCP = NRCP + NRCPI
15       CONTINUE
         NRCSH(P) = NRCP
         DO 20 Q = 1,P
            KHKTQ = KHKTSH(Q)
            MULQ = MULT(ISTBSH(Q))
            NSETQ = NSETSH(Q,1)
            NRCQ = 0
            DO 25 IQ = 1,NSETQ
               NRCQI = NCONTS(Q,IQ,1)
               IF (NRCQI .LT. 0) NRCQI = 1
               NRCQ = NRCQ + NRCQI
25          CONTINUE
            DO 30 R = 1,P
               KHKTR = KHKTSH(R)
               MULR = MULT(ISTBSH(R))
               NSETR = NSETSH(R,1)
               NRCR = 0
               DO 35 IR = 1,NSETR
                  NRCRI = NCONTS(R,IR,1)
                  IF (NRCRI .LT. 0) NRCRI = 1
                  NRCR = NRCR + NRCRI
35             CONTINUE
               MXS = R
               IF (P .EQ. R) MXS = Q
               DO 40 S = 1,MXS
                  KHKTS = KHKTSH(S)
                  MULS = MULT(ISTBSH(S))
                  NSETS = NSETSH(S,1)
                  NRCS = 0
                  DO 45 IS = 1,NSETS
                     NRCSI = NCONTS(S,IS,1)
                     IF (NRCSI .LT. 0) NRCSI = 1
                     NRCS = NRCS + NRCSI
45                CONTINUE
                  NBLOCK = NBLOCK + 1
                  NPS = NPS + 1
                  LENGTH = NRCP*NRCQ*NRCR*NRCS*KHKTP*KHKTQ*KHKTR*KHKTS
     1                    *MULP*MULQ*MULR*MULS
                  NUMP = NUMP + LENGTH
                  IF (MEMTOT + LENGTH .GT. MEMMAX .OR.
     1                NPS .GE. (MAXLOC-1)) THEN
                     NSLICE = NSLICE + 1
                     MEMTOT = LENGTH
                     NPSMAX = MAX(NPS,NPSMAX)
                     NPS = 0
                  ELSE
                     MEMTOT = MEMTOT + LENGTH
                  ENDIF
40             CONTINUE
30          CONTINUE
20       CONTINUE
10    CONTINUE
      CALL PSIZE(MEMMAX,NBLOCK,NSLICE,NUMP,NPASS,MXBUF,MXMAX,
     1           LBUF,IPRINT)
C
C     **************************************
C     **   DYNAMIC ALLOCATION FOR PSORT   **
C     **************************************
C
      NBINS  = (NSLICE+NPASS-1)/NPASS
      NBINS  = MIN(NBINS,MAXADR)
      NPS    = MAX(NPSMAX,NPS)
      NPS1   = NPS + 1
      NPOINT = 2*NPS + 2
      LDABUG = MEMS/NBINS
      LDABUG = MIN(LDABUG,LDAMAX)
      LDABUF = (LDABUG-2)/2
      LDABUF = (LDABUF/2)*2
      LDABUG = LDABUF*2 + 2
      IF (IPRINT .GT. 2)
     &   WRITE(LUPRI,2000) MAXSHL, NSLICE, NPASS, NBINS,
     &   LDABUF, NPS
      CALL GPOPEN(LUDA,'ABAPTRAN2.DA','NEW','DIRECT',' ',
     &   IRAT*LDABUG,OLDDX)
      JINDGEN= 1
      JBUF   = JINDGEN+ 4*NBAST
      JRBUF  = JBUF   + MAXBUF
      JITAB  = JBUF   + LBUF
      JNTOP  = JITAB  + NPS
      JPOINT = JNTOP  + MAXADR
      JBIN   = JPOINT + NPOINT + 1
      JSORT  = JBIN   + LDABUG + 1
      LSORT  = LWORK  - (JSORT - 1)
      CALL PSORT_SET(IWORK(JINDGEN),IPRINT)
      CALL PSORT(IWORK(JBUF),IWORK(JRBUF),LBUF,
     1           IWORK(JBIN),IWORK(JBIN),LDABUF,
     2           IWORK(JPOINT),IWORK(JITAB),
     3           IWORK(JNTOP),IWORK(JSORT),NPS,NPS1,IWORK(JINDGEN),
     4           MEMMAX,NPTOT,IPRINT,LSORT,ANTI,PANTI)
      CALL GPCLOSE(LUDA ,'DELETE')
      CALL GPCLOSE(LUSCR,'DELETE')
      CALL GPCLOSE(LUPSO,'DELETE')
      IF (IPRINT .GE. 02) WRITE(LUPRI,2100) NPTOT
      CALL QEXIT('PSORG')
      RETURN
 1000 FORMAT (12X,'----------  P S O R T ----------',//)
 2000 FORMAT('    Statistics for bin sort of P-elements:',//,
     1       '    Number of shells                    ',I6,/,
     2       '    Number of bin chains                ',I6,/,
     3       '    Number of passes                    ',I6,/,
     4       '    Number of buffers per pass          ',I6,/,
     5       '    Number of elements per buffer       ',I6,/,
     6       '    Maximum number of P blocks in a pass',I6,//)
 2100 FORMAT('  Number of two-electron density',
     *       ' elements sorted:',I8,/)
      END
C  /* Deck psize */
      SUBROUTINE PSIZE(MEMMAX,NBLOCK,NSLICE,NUMP,NPASS,MXBUF,MXMAX,
     1                 LSQBUF,IPRINT)
#include "implicit.h"
C
C     ************************************************
C     **  COMPUTE BUFFER SIZE AND NUMBER OF PASSES  **
C     ************************************************
C
      INTEGER   IXBUF(10), IOOPS(10)
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
#include "priunit.h"
#include "ioweight.h"
      CALL QENTER('PSIZE')
      SQBUF = LSQBUF
      DO 10 I = 1,10
         IXBUF(I)= MIN(MXMAX,(MEMMAX*I - 2*NBLOCK)/NSLICE)
         XBUF  = IXBUF(I)
         DAOPS = 4*NUMP*DAWEIG
         DAOPS = DAOPS/XBUF
         SQOPS = 2*NUMP*I*SQWEIG
         SQOPS = SQOPS/SQBUF
         IOOPS(I) = DAOPS + SQOPS
10    CONTINUE
      IF (IPRINT .GT. 5) THEN
         WRITE(LUPRI,'(A/A/,(3I15))')
     1    '    I/O Requests for different numbers of passes ',
     2    '           Pass   Buffer Size        I/O Ops',
     3    (I, IXBUF(I), IOOPS(I), I = 1,10)
      ENDIF
      IMAX = 0
      IOMAX = 10**9
      DO 20 I = 1,10
         IF (IOOPS(I) .LT. IOMAX) THEN
            IMAX = I
            IOMAX = IOOPS(I)
         ENDIF
20    CONTINUE
      NPASS = IMAX
      MXBUF = IXBUF(IMAX)
      CALL QEXIT('PSIZE')
      RETURN
      END
C  /* Deck PSORT */
      SUBROUTINE PSORT(IBUF,BUF,LBUF,IBIN,BIN,LDABUF,IPOINT,
     1                ITAB,NTOP,SORT,NPOINT,NPS1, INDGEN,
     2                MEMMAX,NPTOT,IPRINT,LSORT,ANTI,PANTI)
#include "implicit.h"
#include "iratdef.h"
#include "aovec.h"
C
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "blocks.h"
#include "symmet.h"
#include "ptrbuf.h"
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
      INTEGER*8 IBUF(*), IBIN(*)
      REAL*8     BUF(*),  BIN(*), SORT(LSORT)
      INTEGER*8 ITAB(*), NTOP(*), IPOINT(NPS1,2)
      INTEGER*8 INDGEN(4,*)
      INTEGER   P, Q, R, S, PQ, RS, PQRS
      LOGICAL   PEQQ, REQS, PQEQRS, ANTI, PANTI
      INTEGER   LAB1, LAB2, LAB3, LAB4
      INTEGER*8 ILOC, INUM, IRAYO
#include "inftap.h"
#include "chrnos.h"
C
C     *****************************************************************
C     **  STATEMENT FUNCTIONS FOR PACKING, UNPACKING AND ADDRESSING  **
C     *****************************************************************
C
      integer*1 i8byte(8)
      integer*8 i8
      equivalence (i8byte, i8)
      INTEGER*8 L_FIELD, IBIT_FIELD
      PARAMETER (L_FIELD = 16, IBIT_FIELD = (2**L_FIELD-1))
      INTEGER   I4, J4, K4, L4
      INTEGER*8 I, J, K, L


      INTEGER*8 IPART
      IPART(I,J4) = IAND(ISHFT(I,-L_FIELD*(J4-1)),IBIT_FIELD)

      INTEGER   ISHLL,ICONTR,IANG,IRREP
      ISHLL(I4)   = INDGEN(4,I4)
      ICONTR(I4)  = INDGEN(3,I4)
      IANG(I4)    = INDGEN(2,I4)
      IRREP(I4)   = INDGEN(1,I4)

      ICR(I4,J4,K4,L4)  = NRCP*(NRCQ*(NRCR*(L4-1)+K4-1)+J4-1)+I4
      IAR(I4,J4,K4,L4)  = KHKTP*(KHKTQ*(KHKTR*(L4-1)+K4-1)+J4-1)+I4
      IRR(I4,J4,K4,L4)  = MULP*(MULQ*(MULR*L4+K4)+J4)+I4+1
      IADDRS(I4,J4,K4)  = NRCTOT*(MULTOT*(K4-1)+J4-1)+I4
      LBIN(I4)   = NDIMI*I4
      LROFF(I4)  = NDIMR*(I4-1)
      LIOFF(I4)  = NDIMI*(I4-1) + NOFFR
      LISOFF(I4) = NDIMI*(I4-1)
C
      CALL QENTER('PSORT')
C
      LDABUG = 2*LDABUF + 2
C
C     CONSTANT OFFSETS AND ADDRESSING
C
      NDIMI = 2*LDABUF + 2
      NDIMR = 2*LDABUF + 2
      NOFFR = LDABUF
C
C     ***************************************
C     **   BIN SORT OF P-MATRIX ELEMENTS   **
C     ***************************************
C
      ITOP = (MAXSHL*(MAXSHL+1))/2
      ITOP = (ITOP*(ITOP+1))/2
      IF (ANTI .OR. PANTI) THEN
         LUPMAT = LUPAS
      ELSE
         LUPMAT = LUPAO
      END IF
      CALL REWSPL(LUPMAT)
C
C     *****************************************
C     **   BEGIN PASSES OVER INTEGRAL FILE   **
C     *****************************************
C
      NPTOT  = 0
      IPQRST = 1
1000  CONTINUE
      IRAYO  = 0
      IPQRS  = 0
      IBINO  = 1
      NTOP(1) = ITOP
      DO 10 P = 1,MAXSHL
         KHKTP = KHKTSH(P)
         MULP  = MULT(ISTBSH(P))
         NRCP  = NRCSH(P)
         DO 20 Q = 1,P
            KHKTQ = KHKTSH(Q)
            MULQ  = MULT(ISTBSH(Q))
            NRCQ  = NRCSH(Q)
            DO 30 R = 1,P
               KHKTR = KHKTSH(R)
               MULR  = MULT(ISTBSH(R))
               NRCR  = NRCSH(R)
               MXS   = R
               IF (P .EQ. R) MXS = Q
               DO 40 S = 1,MXS
                  KHKTS = KHKTSH(S)
                  MULS  = MULT(ISTBSH(S))
                  NRCS  = NRCSH(S)
                  IPQRS = IPQRS + 1
                  IF (IPQRS .LT. IPQRST) GOTO 40
                  NUM = NRCP*NRCQ*NRCR*NRCS*KHKTP*KHKTQ*KHKTR*KHKTS
     1                    *MULP*MULQ*MULR*MULS
                  LOC = IPQRS - IPQRST + 1
                  IF (IRAYO + NUM .LT. MEMMAX) THEN
                     IPOINT(LOC,1) = IBINO
                     IPOINT(LOC,2) = IRAYO
                     IRAYO = IRAYO + NUM
                     IPOINT(LOC+1,2) = IRAYO
                     IF (LOC .GE. NPOINT) GOTO 50
                  ELSE IF (IBINO+1 .LT. MAXADR) THEN
                     NTOP(IBINO) = IRAYO
                     IBINO = IBINO + 1
                     IRAYO = 0
                     IPOINT(LOC,1) = IBINO
                     IPOINT(LOC,2) = IRAYO
                     IRAYO = IRAYO + NUM
                     IPOINT(LOC+1,2) = IRAYO
                     IF (LOC .GE. NPOINT) GOTO 50
                  ELSE
                     IPQRS = IPQRS - 1
                     NTOP(IBINO) = IRAYO
                     GOTO 50
                  ENDIF
40             CONTINUE
30          CONTINUE
20       CONTINUE
10    CONTINUE
50    CONTINUE
      NPOINU = MIN(LOC,NPOINT)
C
C     ***********************************************************
C     **   SORT RANGE IS NOW DEFINED, BEGIN READING P-MATRIX   **
C     ***********************************************************
C
      IPQRSF = IPQRS
      NCHAIN = IBINO
      IDISK  = 1
      DO, ICHAIN = 1,NCHAIN
         IBIN(LBIN(ICHAIN))   = 0
         IBIN(LBIN(ICHAIN)-1) = -1
      END DO

      CALL REWSPL(LUPSO)
      JREC = 0
100   CONTINUE
      JREC = JREC + 1
      CALL READI8(LUPSO,LBUF,IBUF)
      NUT = IBUF(LBUF)
      IF (NUT .EQ. 0) GOTO 100
      IF (NUT .LT. 0) GOTO 200
      DO 70 I = 1,NUT
         PVAL = BUF(I)
         LAB4 = IPART(IBUF(I),1)
         LAB3 = IPART(IBUF(I),2)
         LAB2 = IPART(IBUF(I),3)
         LAB1 = IPART(IBUF(I),4)
         P = ISHLL(LAB1)
         Q = ISHLL(LAB2)
         R = ISHLL(LAB3)
         S = ISHLL(LAB4)
         IF (P .LT. Q) THEN
            IH = LAB1
            LAB1 = LAB2
            LAB2 = IH
            IH = P
            P  = Q
            Q  = IH
            IF (ANTI) PVAL = - PVAL
         ENDIF
         IF (R .LT. S) THEN
            IH   = LAB3
            LAB3 = LAB4
            LAB4 = IH
            IH   = R
            R    = S
            S    = IH
            IF (ANTI) PVAL = - PVAL
         ENDIF
         IF (P .LT. R .OR. (P .EQ. R .AND. Q .LT. S)) THEN
            IH   = LAB1
            LAB1 = LAB3
            LAB3 = IH
            IH   = LAB2
            LAB2 = LAB4
            LAB4 = IH
            IH   = P
            P    = R
            R    = IH
            IH   = Q
            Q    = S
            S    = IH
            IF (PANTI) PVAL = -PVAL
         ENDIF
         PQ = (P*(P-1))/2 + Q
         RS = (R*(R-1))/2 + S
         PQRS = (PQ*(PQ-1))/2 + RS
C
C        CHECK SPECIAL CASES
C
         PEQQ   = P .EQ. Q
         REQS   = R .EQ. S
         PQEQRS = (P .EQ. R) .AND. (Q .EQ. S)
         ICSPQ  = 0
         ICSRS  = 0
         ICSPR  = 0
         IF (PQRS .LT. IPQRST .OR. PQRS .GT. IPQRSF) GOTO 70
C
C        PQRS IS IN RANGE FOR THIS PASS
C
         IND   = PQRS - IPQRST + 1
         IBINO = IPOINT(IND,1)
         IRAYO = IPOINT(IND,2)
C
C        EXTRACT SHELL COMPONENT INFORMATION FROM LABELS AND
C        COMPUTE SORT ADDRESS AND OFFSET
C
         MULP = MULT(ISTBSH(P))
         MULQ = MULT(ISTBSH(Q))
         MULR = MULT(ISTBSH(R))
         MULS = MULT(ISTBSH(S))
         NRCP = NRCSH(P)
         NRCQ = NRCSH(Q)
         NRCR = NRCSH(R)
         NRCS = NRCSH(S)
         KHKTP = KHKTSH(P)
         KHKTQ = KHKTSH(Q)
         KHKTR = KHKTSH(R)
         KHKTS = KHKTSH(S)
         IP = ICONTR(LAB1)
         IQ = ICONTR(LAB2)
         IR = ICONTR(LAB3)
         IS = ICONTR(LAB4)
         KP = IANG(LAB1)
         KQ = IANG(LAB2)
         KR = IANG(LAB3)
         KS = IANG(LAB4)
         NP = IRREP(LAB1)
         NQ = IRREP(LAB2)
         NR = IRREP(LAB3)
         NS = IRREP(LAB4)
         NRCTOT = ICR(NRCP,NRCQ,NRCR,NRCS)
         MULTOT = IRR(MULP-1,MULQ-1,MULR-1,MULS-1)
700      CONTINUE
         ILOC = IRAYO + IADDRS(ICR(IP,IQ,IR,IS),
     1                         IRR(NP,NQ,NR,NS),
     2                         IAR(KP,KQ,KR,KS))
         INUM = IBIN(LBIN(IBINO)) + 1
         BIN( LROFF(IBINO)+INUM) = PVAL
         IBIN(LIOFF(IBINO)+INUM) = ILOC
         IBIN(LBIN(IBINO))       = INUM
         IF (INUM .EQ. LDABUF) THEN
            CALL WRITDX(LUDA,IDISK,IRAT*LDABUG,IBIN(LISOFF(IBINO)+1))
            IBIN(LBIN(IBINO))   = 0
            IBIN(LBIN(IBINO)-1) = IDISK
            IDISK = IDISK + 1
         ENDIF
C
C        CHECK FOR COINCIDENCES AND BRANCH BACK
C
         IF (PEQQ .AND. ICSPQ .EQ. 0) THEN
            IF (ANTI) PVAL = - PVAL
            IH = IP
            IP = IQ
            IQ = IH
            IH = KP
            KP = KQ
            KQ = IH
            IH = NP
            NP = NQ
            NQ = IH
            ICSPQ = 1
            GOTO 700
         ENDIF
         IF (REQS .AND. ICSRS .EQ. 0) THEN
            IF (ANTI) PVAL = - PVAL
            IH = IR
            IR = IS
            IS = IH
            IH = KR
            KR = KS
            KS = IH
            IH = NR
            NR = NS
            NS = IH
            ICSRS = 1
C
C           Need to reset ICSPQ here
C
            ICSPQ = 0
            GOTO 700
         ENDIF
         IF (PQEQRS .AND. ICSPR .EQ. 0) THEN
C
C           At this point we must already have interchanged
C           both the PQ shell information and the RS.  We
C           therefore need only to interchange the pairs here.
C
            IF (PANTI) PVAL = -PVAL
            IH = IP
            IP = IR
            IR = IH
            IH = IQ
            IQ = IS
            IS = IH
            IH = KP
            KP = KR
            KR = IH
            IH = KQ
            KQ = KS
            KS = IH
            IH = NP
            NP = NR
            NR = IH
            IH = NQ
            NQ = NS
            NS = IH
            ICSPQ = 0
            ICSRS = 0
            ICSPR = 1
            GOTO 700
         ENDIF
70    CONTINUE
      GOTO 100
200   CONTINUE
      DO 75 ICHAIN = 1,NCHAIN
         INUM = IBIN(LBIN(ICHAIN))
         IF (INUM .GT. 0) THEN
            CALL WRITDX(LUDA,IDISK,IRAT*LDABUG,IBIN(LISOFF(ICHAIN)+1))
            IBIN(LBIN(ICHAIN)-1) = IDISK
            IDISK = IDISK + 1
         ENDIF
         LASTAD(ICHAIN) = IBIN(LBIN(ICHAIN)-1)
75    CONTINUE
      REWIND LUSCR
C
C     PREPARE OFFSET AND BLOCK LENGTH TABLES
C
      DO 400 ICHAIN = 1,NCHAIN
         NTAB = 0
         DO 410 IOFF = 1,NPOINU
            IF (IPOINT(IOFF,1) .EQ. ICHAIN) THEN
               NTAB = NTAB + 1
               ITAB(NTAB) = IPOINT(IOFF+1,2) - IPOINT(IOFF,2)
               IF (ITAB(NTAB) .LE. 0)
     1             ITAB(NTAB) = NTOP(ICHAIN) - IPOINT(IOFF,2)
            ENDIF
410      CONTINUE
         WRITE(LUSCR) NTAB, (ITAB(I),I = 1,NTAB)
400   CONTINUE
C
C     ****************************************************
C     **   SORT FILE FOR THIS PASS IS COMPLETE          **
C     **   READ BACK CHAINS AND WRITE SORTED P-MATRIX   **
C     ****************************************************
C
      REWIND LUSCR
      JPQRS = IPQRST - 1
      DO 80 ICHAIN = 1,NCHAIN
         READ(LUSCR) NTAB, (ITAB(I), I = 1,NTAB)
C
C        INITIALIZE SORT SPACE TO ZERO
C
         IFIN = 0
         DO 83 I = 1,NTAB
            IFIN = IFIN + ITAB(I)
83       CONTINUE
         IF (IFIN .GT. LSORT) CALL STOPIT('PSORT',' ',IFIN,LSORT)
         CALL DZERO(SORT,IFIN)
         IF (LASTAD(ICHAIN) .EQ. -1) GOTO 95
C
         IDISK  = LASTAD(ICHAIN)
   85    CONTINUE
         CALL READDX(LUDA,IDISK,IRAT*LDABUG,IBIN)
         LENGTH = IBIN(LBIN(1))
         DO 90 I = 1,LENGTH
            ILOC = IBIN(LIOFF(1)+I)
            SORT(ILOC) = BIN(I)
   90    CONTINUE
         IDISK = IBIN(LBIN(1)-1)
         IF (IDISK .NE. -1) GOTO 85
C
C        ARRAYS ARE SET UP - WRITE OUT TO SEQUENTIAL FILE
C
95       CONTINUE
         IFIN = 0
         DO 120 I = 1,NTAB
            NBL   = ITAB(I)
            IST   = IFIN  + 1
            IEND  = IFIN  + MAX(4,NBL)
            IFIN  = IFIN  + NBL
            WRITE(LUPMAT) NBL
            WRITE(LUPMAT) (SORT(J), J = IST,IEND)
            IF (IPRINT .GT. 30) THEN
               JPQRS = JPQRS + 1
               WRITE(LUPRI,'(A,I10,/,(4X,1P,5D12.4))')
     &          ' SO 2-matrix block for shell quadruplet', JPQRS,
     &          (SORT(J),J = IST,IFIN)
            END IF
  120    CONTINUE
         NPTOT = NPTOT + IFIN
   80 CONTINUE
      IPQRST = IPQRSF + 1
      IF (IPQRSF .LT. ITOP) GOTO 1000
      CALL QEXIT('PSORT')
      RETURN
      END
C  /* Deck PSORT_SET */
      SUBROUTINE PSORT_SET(INDGEN,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      INTEGER*8 INDGEN(4,*)
      INTEGER   IRREP(MXCORB)
#include "shells.h"
#include "symmet.h"
#include "ptrbuf.h"

C
      CALL QENTER('PSORT_SET')
      JORB = 0
      CALL IZERO(IRREP,MXCORB)
      DO 100 IREP = 0, MAXREP
         IK = 0
         DO 200 I = 1, KMAX
            MULA   = ISTBAO(I)
            NUMCFA = NUMCF(I)
            DO 300 K = 1, KHKT(I)
               IK = IK + 1
               IF(IAND(MULA,IEOR(IREP,ISYMAO(NHKT(I),K))).EQ.0)THEN
                  JORB = JORB + 1
                  INDGEN(1,JORB) = IRREP(IK)
                  INDGEN(2,JORB) = K
                  INDGEN(3,JORB) = NUMCFT(I)
                  INDGEN(4,JORB) = IPTSHL(I)
                  IRREP(IK)      = IRREP(IK) + 1
               END IF
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      CALL QEXIT('PSORT_SET')
      RETURN
      END
C  /* Deck ptrtst */
      SUBROUTINE PTRTST(WORK,LWORK,ANTI,PANTI,DIA2SO,ZFS2EL,
     &                  MAXBUF,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL   ANTI, PANTI, DIA2SO, ZFS2EL
      REAL*8    WORK(LWORK)
#include "inforb.h"
C
      CALL QENTER('PTRTST')
      LBUF  = 2*MAXBUF + 1
C
      KPVAO  = 1
      KPVTRA = KPVAO  + N2BASX*N2BASX
      KIBUF  = KPVTRA + N2BASX*N2BASX
      KRBUF  = KIBUF  + MAXBUF
      KLAST  = KIBUF  + LBUF
      LWRK   = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('PTRTST',' ',KLAST,LWORK)
C
C     Noddy transformation of densities to AO basis
C
      CALL PTRNOD(WORK(KPVAO),WORK(KLAST),LWRK,ANTI,PANTI,DIA2SO,
     &            ZFS2EL,IPRINT)
C
C     Comparison with transformed integrals
C
      CALL PVCMP(WORK(KPVAO),WORK(KPVTRA),WORK(KIBUF),WORK(KRBUF),LBUF,
     &           ANTI,PANTI,IPRINT)
C
      CALL QEXIT('PTRTST')
      RETURN
      END
C  /* Deck ptrnod */
      SUBROUTINE PTRNOD(PVAO,WORK,LWORK,ANTI,PANTI,DIA2SO,ZFS2EL,IPRINT)
#include "implicit.h"
#include "priunit.h"
      LOGICAL   ANTI, PANTI, DIA2SO, ZFS2EL, FOUND
      REAL*8    PVAO(N2BASX,N2BASX), WORK(LWORK)
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
#include "inftap.h"
#include "inforb.h"
C
      CALL QENTER('PTRNOD')
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (IDUM, I = 1,9), NCMOT,NNASHX
C
      KCMO  = 1
      KTRA  = KCMO  + NCMOT
      KPV   = KTRA  + NASHT*NBAST
      KPVMO = KPV   + NNASHX*NNASHX
      KLAST = KPVMO + N2ASHX*N2ASHX
      LWRK  = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('PTRNOD',' ',KLAST,LWORK)
C
      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KLAST),LWRK)
      IF (.NOT.FOUND) CALL QUIT('PTRNOD error: CMO not found on SIRIFC')
      IF (ANTI .OR. DIA2SO .OR. ZFS2EL) THEN
         CALL MAKPVM(WORK(KPV),WORK(KLAST),LWRK,ANTI,PANTI,DIA2SO,
     &      ZFS2EL, IPRINT)
      ELSE
         CALL RD_SIRIFC('PV',FOUND,WORK(KPV),WORK(KLAST),LWRK)
         IF (.NOT.FOUND)
     &        CALL QUIT('PTRNOD error: PV not found on SIRIFC')
      END IF
      CALL PTRNO1(WORK(KPV),WORK(KPVMO),PVAO,WORK(KCMO),
     &            WORK(KTRA),ANTI,PANTI,DIA2SO,ZFS2EL,IPRINT)
      CALL QEXIT('PTRNOD')
      RETURN
      END
C  /* Deck ptrno1 */
      SUBROUTINE PTRNO1(PV,PVMO,PVAO,CMO,TRA,ANTI,PANTI,DIA2SO,ZFS2EL,
     &                  IPRINT)
#include "implicit.h"
#include "priunit.h"
      LOGICAL   ANTI, PANTI, DIA2SO, ZFS2EL
      INTEGER   T, U, V, X, TU, VX
      REAL*8    PV(NNASHX,NNASHX), CMO(NCMOT),
     &          PVMO(NASHT,NASHT,NASHT,NASHT),
     &          PVAO(NBAST,NBAST,NBAST,NBAST),
     &          TRA(NBAST,NASHT)
#include "inforb.h"
      CALL QENTER('PTRNO1')
C
C     Matrix in MO basis
C
      TU = 0
      DO 100 T = 1, NASHT
      DO 100 U = 1, T
         TU = TU + 1
         VX = 0
         DO 200 V = 1, NASHT
         DO 200 X = 1, V
            VX = VX + 1
	       IF (DIA2SO) THEN
		  PVMO(T,U,V,X) = PV(TU,VX)
		  PVMO(U,T,V,X) = PV(TU,VX)
		  PVMO(T,U,X,V) = PV(TU,VX)
		  PVMO(U,T,X,V) = PV(TU,VX)
	       ELSE
		  PVMO(T,U,V,X) = PV(TU,VX)
		  PVMO(U,T,X,V) = PV(TU,VX)
		  IF (ANTI) THEN 
		     PVMO(U,T,V,X) = - PV(TU,VX)
		     PVMO(T,U,X,V) = - PV(TU,VX)
		  ELSE
		     PVMO(U,T,V,X) = PV(TU,VX)
		     PVMO(T,U,X,V) = PV(TU,VX)
		  END IF
	       END IF
  200    CONTINUE
  100 CONTINUE
      CALL HEADER('PV matrix in MO basis in PTRNO1',-1)
      CALL OUTPUT(PVMO,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
C
C     Matrix in AO basis
C
C
      CALL DZERO(TRA,NASHT*NBAST)
      DO 500 ISYM = 1, NSYM
         CALL HEADER('CMO in PTRNO1',-1)
         CALL OUTPUT(CMO(ICMO(ISYM)+1),1,NBAS(ISYM),1,NORB(ISYM),
     &               NBAS(ISYM),NORB(ISYM),1,LUPRI)
         IT = ICMO(ISYM) + NBAS(ISYM)*NISH(ISYM)
         DO 600 T = 1, NASH(ISYM)
            DO 700 I = 1, NBAS(ISYM)
               IT = IT + 1
               TRA(IBAS(ISYM) + I,IASH(ISYM) + T) = CMO(IT)
  700       CONTINUE
  600    CONTINUE
  500 CONTINUE
C
      CALL HEADER('TRA in PTRNO1',-1)
      CALL OUTPUT(TRA,1,NBAST,1,NASHT,NBAST,NASHT,1,LUPRI)
C
      CALL DZERO(PVAO,N2BASX*N2BASX)
      DO 300 I = 1, NBAST
      DO 300 J = 1, NBAST
      DO 300 K = 1, NBAST
      DO 300 L = 1, NBAST
         DO 400 T = 1, NASHT
         DO 400 U = 1, NASHT
         DO 400 V = 1, NASHT
         DO 400 X = 1, NASHT
            PVAO(I,J,K,L) = PVAO(I,J,K,L) + PVMO(T,U,V,X)
     &                                *TRA(I,T)*TRA(J,U)
     &                                *TRA(K,V)*TRA(L,X)
  400    CONTINUE
  300 CONTINUE
      CALL HEADER('PV matrix in AO basis in PTRNO1',-1)
      CALL OUTPUT(PVAO,1,N2BASX,1,N2BASX,N2BASX,N2BASX,1,LUPRI)
      CALL QEXIT('PTRNO1')
      RETURN
      END
C  /* Deck pvcmp */
      SUBROUTINE PVCMP(PVAO,PVTRA,IBUF,BUF,LBUF,ANTI,PANTI,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D1 = 1.0D0)
      LOGICAL   ANTI, PANTI
      INTEGER*8 IBUF(*), P,Q,R,S, NUT
      REAL*8    PVAO(NBAST,NBAST,NBAST,NBAST), BUF(*),
     &          PVTRA(NBAST,NBAST,NBAST,NBAST)
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
#include "ccom.h"
#include "inforb.h"

      INTEGER*8  L_FIELD, IBIT_FIELD
      PARAMETER (L_FIELD = 16, IBIT_FIELD = (2**L_FIELD-1))
      INTEGER*8 I, J

C
      CALL QENTER('PVCMP')
C
      SGN = D1
      PSGN = D1
      IF (ANTI) SGN = - D1
      IF (PANTI) PSGN = - D1
      CALL DZERO(PVTRA,NBAST**4)
C
      REWIND LUPSO
  100 CONTINUE
      CALL READI8(LUPSO,LBUF,IBUF)
      NUT = IBUF(LBUF)
      IF (NUT .EQ. 0) GOTO 100
      IF (NUT .LT. 0) GOTO 200
      DO 300 I = 1,NUT
         P = IAND(ISHFT(IBUF(I),-3*L_FIELD),IBIT_FIELD)
         Q = IAND(ISHFT(IBUF(I),-2*L_FIELD),IBIT_FIELD)
         R = IAND(ISHFT(IBUF(I),  -L_FIELD),IBIT_FIELD)
         S = IAND(       IBUF(I)           ,IBIT_FIELD)
         PVTRA(P,Q,R,S) =      BUF(I)
         PVTRA(P,Q,S,R) = SGN*BUF(I)
         PVTRA(Q,P,R,S) = SGN*BUF(I)
         PVTRA(Q,P,S,R) =      BUF(I)
         PVTRA(R,S,P,Q) =      BUF(I)*PSGN
         PVTRA(R,S,Q,P) = SGN*BUF(I)*PSGN
         PVTRA(S,R,P,Q) = SGN*BUF(I)*PSGN
         PVTRA(S,R,Q,P) =      BUF(I)*PSGN
  300 CONTINUE
      GO TO 100
  200 CONTINUE
C
      DIFMAX = 0.0D0
      DO 400 P = 1, NBAST
      DO 400 Q = 1, NBAST
      DO 400 R = 1, NBAST
      DO 400 S = 1, NBAST
         DIFFER = PVTRA(P,Q,R,S) - PVAO(P,Q,R,S)
         DIFMAX = MAX(DIFMAX,ABS(DIFFER))
         IF (ABS(DIFFER) .GT. THRS) THEN
            WRITE (LUPRI,'(2X,3E12.6,10X,4I5)')
     &         PVTRA(P,Q,R,S), PVAO(P,Q,R,S),DIFFER,P,Q,R,S
         END IF
  400 CONTINUE
      CALL HEADER('PV matrix in AO basis from PT2TRA',-1)
      CALL OUTPUT(PVTRA,1,N2BASX,1,N2BASX,N2BASX,N2BASX,1,LUPRI)
      CALL DAXPY(N2BASX**2,-D1,PVAO,1,PVTRA,1)
      CALL HEADER('PV difference matrix',-1)
      CALL OUTPUT(PVTRA,1,N2BASX,1,N2BASX,N2BASX,N2BASX,1,LUPRI)
      WRITE (LUPRI,'(2X,A,E20.12)') ' Largest difference found:',DIFMAX
      CALL QEXIT('PVCMP')
      RETURN
      END
C  /* Deck makpvm */
      SUBROUTINE MAKPVM(PV,WORK,LWORK,ANTI,PANTI,DIA2SO,ZFS2EL,IPRINT)
#include "implicit.h"
#include "priunit.h"
      LOGICAL   ANTI, PANTI, DIA2SO, ZFS2EL
      REAL*8    PV(*), WORK(LWORK)
#include "inflin.h"
#include "infdim.h"
#include "inforb.h"
      CALL QENTER('MAKPVM')
      KCINDX = 1
      KCREF  = KCINDX + LCINDX
      KUDV   = KCREF  + NCONRF
      KPVSQ  = KUDV   + N2ASHX
      KWRK   = KPVSQ  + N2ASHX*N2ASHX
      LWRK   = LWORK  - KWRK + 1
      IF (KWRK.GT.LWORK) CALL STOPIT('MAKPVM',' ',LWRK,LWORK)
      IF (DIA2SO) THEN
         CALL MAKPV2(PV,WORK(KCINDX),WORK(KCREF),WORK(KUDV),WORK(KPVSQ),
     &               WORK(KWRK),LWRK,PANTI,IPRINT)
      ELSE IF (ZFS2EL) THEN
         CALL MAKPVQ(PV,WORK(KCINDX),WORK(KCREF),WORK(KUDV),WORK(KPVSQ),
     &            WORK(KWRK),LWRK,IPRINT)
      ELSE
         CALL MAKPV1(PV,WORK(KCINDX),WORK(KCREF),WORK(KUDV),WORK(KPVSQ),
     &               WORK(KWRK),LWRK,ANTI,IPRINT)
      END IF
      CALL QEXIT('MAKPVM')
      RETURN
      END
C  /* Deck makpv1 */
      SUBROUTINE MAKPV1(PV,XINDX,CREF,DV,PVSQ,WORK,LWORK,ANTI,IPRINT)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      LOGICAL   ANTI
      INTEGER   T,U,V,X, TU, VX
      REAL*8    XINDX(*)
      REAL*8    PV(NNASHX,NNASHX), PVSQ(NASHT,NASHT,NASHT,NASHT),
     &          CREF(NCONRF), DV(N2ASHX), WORK(LWORK)
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
#include "inflin.h"
#include "infdim.h"
#include "inforb.h"
#include "inftap.h"
C
      CALL QENTER('MAKPV1')
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI8(LUSIFC,NCONRF,CREF)
C
      CALL GETCIX(XINDX,LSYMRF,LSYMRF,WORK,LWORK,0)
C
      ISPIN1 = 0
      ISPIN2 = 0
      NVAR   = 1
      CALL RSPDM(LSYMRF,LSYMRF,NCONRF,NCONRF,CREF,CREF,DV,PVSQ,
     &           ISPIN1,ISPIN2,.FALSE.,.FALSE.,XINDX,WORK,NVAR,
     &           LWORK)
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Complete PV matrix MAKPVM',-1)
         CALL OUTPUT(PVSQ,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
      END IF
C
      SGN = D1
      IF (ANTI) SGN = - 1
      TU = 0
      DO 100 T = 1, NASHT
      DO 100 U = 1, T
         TU = TU + 1
         VX = 0
         DO 200 V = 1, NASHT
         DO 200 X = 1, V
            VX = VX + 1
            PV(TU,VX) = DP5*(PVSQ(T,U,V,X) + SGN*PVSQ(T,U,X,V))
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('PV matrix MAKPVM',-1)
         CALL OUTPUT(PV,1,NNASHX,1,NNASHX,NNASHX,NNASHX,-1,LUPRI)
      END IF
      CALL QEXIT('MAKPV1')
      RETURN
      END
C  /* Deck makpv2 */
      SUBROUTINE MAKPV2(PV,XINDX,CREF,DV,PVSQ,WORK,LWORK,PANTI,IPRINT)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP5 = 0.5D0, D1=1.0D0, D2 = 2.0D0, DP25 = 0.25D0)
      LOGICAL   PANTI
      INTEGER   T,U,V,X, TU, VX
      REAL*8    XINDX(*)
      REAL*8    PV(NNASHX,NNASHX), PVSQ(NASHT,NASHT,NASHT,NASHT),
     &          CREF(NCONRF), DV(N2ASHX), WORK(LWORK)
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
#include "inflin.h"
#include "infdim.h"
#include "inforb.h"
#include "inftap.h"
C
      CALL QENTER('MAKPV2')
C
      KFREE  = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KPV,N2ASHX*N2ASHX,WORK,KFREE,LFREE)
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI8(LUSIFC,NCONRF,CREF)
C
      CALL GETCIX(XINDX,LSYMRF,LSYMRF,WORK(KFREE),LFREE,0)
C
      ISPIN1 = 1
      ISPIN2 = 0
      CALL RSPDM(LSYMRF,LSYMRF,NCONRF,NCONRF,CREF,CREF,DV,PVSQ,
     &           ISPIN1,ISPIN2,.FALSE.,.FALSE.,XINDX,WORK,KFREE,
     &           LFREE)
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Square PV matrix for electron 1 in MAKPV2',-1)
         CALL OUTPUT(PVSQ,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
      END IF
C
C F90:  PV = PV + 2*RESHAPE(PV,(/NASHT,NASHT,NASHT,NASHT/),ORDER=(/3,4,1,2/))
C
      CALL MTRSP(N2ASHX,N2ASHX,PVSQ,N2ASHX,WORK(KPV),N2ASHX)
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Square PV matrix for electron 2 in MAKPV2',-1)
         CALL OUTPUT(WORK(KPV),1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
      END IF
      CALL DAXPY(N2ASHX*N2ASHX,D2,WORK(KPV),1,PVSQ,1)
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Complete square PV matrix in MAKPV2',-1)
         CALL OUTPUT(PVSQ,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
      END IF
      CALL MEMREL('MAKPV2',WORK,1,1,KFREE,LFREE)
C
      SGN = D1
      IF (PANTI) SGN = -D1
      TU = 0
      DO 100 T = 1, NASHT
      DO 100 U = 1, T
         TU = TU + 1
         VX = 0
         DO 200 V = 1, NASHT
         DO 200 X = 1, V
            VX = VX + 1
            PV(TU,VX) = DP25*(( PVSQ(T,U,V,X) + PVSQ(T,U,X,V))
     &                  + SGN*( PVSQ(V,X,U,T) + PVSQ(X,V,U,T)))
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('PV matrix MAKPVM',-1)
         CALL OUTPUT(PV,1,NNASHX,1,NNASHX,NNASHX,NNASHX,1,LUPRI)
      END IF
      CALL QEXIT('MAKPV2')
      RETURN
      END
C  /* Deck mem_psort */
      SUBROUTINE MEM_PSORT(MEMMIN,NCONTS)
C
C     Determine minimum amount of memory needed for the PSORT sorting program
C     K.Ruud, July-00
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
      INTEGER   NCONTS(MXSHEL,MXAOVC,2)
      INTEGER   P, Q, R, S
#include "symmet.h"
#include "blocks.h"
C
      MEMMIN = 0
C
      DO 10 P = 1,MAXSHL
         KHKTP = KHKTSH(P)
         MULP  = MULT(ISTBSH(P))
         NSETP = NSETSH(P,1)
         NRCP = 0
         DO IP = 1, NSETP
            NRCPI = NCONTS(P,IP,1)
            IF (NRCPI .LT. 0) NRCPI = 1
            NRCP = NRCP + NRCPI
         END DO
         DO 20 Q = 1,P
            KHKTQ = KHKTSH(Q)
            MULQ  = MULT(ISTBSH(Q))
            NSETQ = NSETSH(Q,1)
            NRCQ  = 0
            DO IQ = 1, NSETQ
               NRCQI = NCONTS(Q,IQ,1)
               IF (NRCQI .LT. 0) NRCQI = 1
               NRCQ = NRCQ + NRCQI
            END DO
            DO 30 R = 1,P
               KHKTR = KHKTSH(R)
               MULR  = MULT(ISTBSH(R))
               NSETR = NSETSH(R,1)
               NRCR  = 0
               DO IR = 1, NSETR
                  NRCRI = NCONTS(R,IR,1)
                  IF (NRCRI .LT. 0) NRCRI = 1
                  NRCR = NRCR + NRCRI
               END DO
               MXS   = R
               IF (P .EQ. R) MXS = Q
               DO 40 S = 1,MXS
                  KHKTS = KHKTSH(S)
                  MULS  = MULT(ISTBSH(S))
                  NSETS = NSETSH(S,1)
                  NRCS  = 0
                  DO IS = 1, NSETS
                     NRCSI = NCONTS(S,IS,1)
                     IF (NRCSI .LT. 0) NRCSI = 1
                     NRCS = NRCS + NRCSI
                  END DO
                  NUM = NRCP*NRCQ*NRCR*NRCS*KHKTP*KHKTQ*KHKTR*KHKTS
     &                    *MULP*MULQ*MULR*MULS
                  MEMMIN = MAX(NUM,MEMMIN)
 40            CONTINUE
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
      RETURN
      END
C  /* Deck makpvq */
      SUBROUTINE MAKPVQ(PV,XINDX,CREF,DV,PVSQ,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
      PARAMETER (DM1=-1.0D0, DP25=0.25D0)
      INTEGER   T,U,V,X, TU, VX
      REAL*8    XINDX(*)
      REAL*8    PV(NNASHX,NNASHX), PVSQ(NASHT,NASHT,NASHT,NASHT),
     &          CREF(NCONRF), DV(N2ASHX), WORK(LWORK)
      COMMON /PTRFIL/ LUSCR, LUDA, LUPSO
#include "inflin.h"
#include "infdim.h"
#include "inforb.h"
#include "inftap.h"
#include "infinp.h"
C
      CALL QENTER('MAKPVQ')
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI8(LUSIFC,NCONRF,CREF)
C
      IF (.NOT. HSROHF) CALL GETCIX(XINDX,LSYMRF,LSYMRF,WORK,LWORK,0)
C
C Create 2-eletron quintet spin density 2sz1*sz2 - sx1*sx2 - sy1sy2
C (zero component of tensor product of rank 2)
C spherical descripiton: 2s(0)*s(0) + s(1)s(-1) + s(-1)s(1)
C                                   - s(+)s(-)/2 - s(-)s(+)/2
C 
C The first term is a standard call to RSPDM
C
C     2sz1*sz2 -> 0.5d--(p,q,r,s)
C
C In evaluating the last two terms with RSPDM we use
C
C     pa qb rb sa + pb qa ra sb - E(p,s)d(r,q)
C   = - pa sa rb qb - pb sb ra qa
C and
C     aabb + bbaa = [ (aa+bb)(aa+bb) - (aa-bb)(aa-bb) ] / 2
C 
C     d(p,q,r,s) = 0.5* d--(p,q,r,s) 
C
      KFREE=1
      LFREE=LWORK
      N4=N2ASHX*N2ASHX
      CALL MEMGET('REAL',KPVTMP,N4,WORK,KFREE,LFREE)
      ISPIN1 = 0
      ISPIN2 = 0
      NVAR   = 1
      CALL RSPDM(LSYMRF,LSYMRF,NCONRF,NCONRF,CREF,CREF,DV,PVSQ,
     &           ISPIN1,ISPIN2,.FALSE.,.FALSE.,XINDX,WORK,KFREE,
     &           LFREE)

      ISPIN1 = 1
      ISPIN2 = 1
      CALL RSPDM(LSYMRF,LSYMRF,NCONRF,NCONRF,CREF,CREF,DV,WORK(KPVTMP),
     &           ISPIN1,ISPIN2,.FALSE.,.FALSE.,XINDX,WORK,KFREE,
     &           LFREE)
      IF (IPRINT.GT.20) THEN
         CALL HEADER('MAKPVQ: singlet density',-1)
         CALL OUTPUT(PVSQ,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
         CALL HEADER('MAKPVQ: triplet density',-1)
         CALL OUTPUT(WORK(KPVTMP),1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,
     &      LUPRI)
      END IF

      CALL DAXPY(N4,DM1,WORK(KPVTMP),1,PVSQ,1)
      CALL DSCAL(N4,DP25,PVSQ,1)
      PVSQ=RESHAPE(PVSQ,(/NASHT,NASHT,NASHT,NASHT/),ORDER=(/1,4,3,2/))
      CALL DAXPY(N4,DP5,WORK(KPVTMP),1,PVSQ,1)
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Complete PV matrix MAKPVM',-1)
         CALL OUTPUT(PVSQ,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
      END IF
      CALL MEMREL('MAKPVQ',WORK,1,1,KFREE,LFREE)
C
      TU = 0
      DO 100 T = 1, NASHT
      DO 100 U = 1, T
         TU = TU + 1
         VX = 0
         DO 200 V = 1, NASHT
         DO 200 X = 1, V
            VX = VX + 1
            PV(TU,VX) = DP5*(PVSQ(T,U,V,X) + PVSQ(T,U,X,V))
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('PV matrix MAKPVM',-1)
         CALL OUTPUT(PV,1,NNASHX,1,NNASHX,NNASHX,NNASHX,1,LUPRI)
      END IF
      CALL QEXIT('MAKPVQ')
      RETURN
      END
      SUBROUTINE PERM24(N,P,Q)
#include "implicit.h"
      REAL*8    P(N,N,N,N), Q(N,N,N,N)
      DO J=1,N
      DO K=1,N
      DO L=1,N
      DO I=1,N
         Q(I,L,K,J) = P(I,J,K,L)
      END DO
      END DO
      END DO
      END DO
      N4=N*N*N*N
      CALL DCOPY(N4,Q,1,P,1)
      END
