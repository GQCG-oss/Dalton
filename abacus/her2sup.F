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
#ifdef REV_LOG
c===========================================================================
c031021-hjaaj: removed OLDSUP
c950605-vebjornb:
c      Added test to check if LUSUPM already open to avoid error message.
c940818-hjaaj:
c- new label: IRIX for Silicon Graphics
c- Removed check for duplicate integrals (reactivate with *DEFINE DUPLCHCK)
c- move symmetry check in SUPSRT in front of canonical index ordering
c- passed NOSYM to SUPSRT, skip symmetry check if NOSYM
c- changed def. of ZERO,ONE,ON2,...,TH8 from DATA to PARAMETER
c- SUPSRT: use that hermit already has put the 4 indices in canonical
c  order (when SUPSRT's def. of np,nq,nr,ns was reversed).
c- SUP1,SUP,SUP1D,SUPD: reverse order of np,nq,nr,ns in label
c- removed N8,N16,N24 from /CSPINT/ (not used)
c- maximize NST1 for LWORK
cnov 90-hjaaj:
cSUPDRV: LINT = 1706 instead of 540
cSUPNWR: fixed errors for NOSYM false
c900319-hjaaj: corrected error in NPLAST for SUP1D and SUPD
c- rewritten SUPNWR to write up to NP=NST-1 as type 1 if advantageous,
c  used new routine NDXGTA to determine this.
c- added WRSUP2 to correcte error in SUPNWR
c900112-hjaaj: new routines SUPOWR and SUPNWR for writing to LUSUPM
cwhen OLDSUP and .NOT.OLDSUP, resp.
c- use new variable NBAS1 = true NBAS(1) for NST in SUPSRT
c900111-hjaaj: added NPLAST to CSPINT
c- Format of LUSUPM if .not.oldsup:
crec 1 : '*******','date    ','time    ','PSUPMAT '
crec 2 : nosym,nsym,nbas(1:8)
cfor i = 3 until finished
crec ia: ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS
crec ib: if (ityp.eq.1) P(ibef+1:ibef+npqrs)
c        if (ityp.eq.2) P(1:npqrs),INDP(1:npqrs) and np1=npl, nq1=nql
crec n : -2
c900109-hjaaj: s/SRTINT/SUPSRT/;s/ONEL/SUPONE/;s/RDINFO/SUPRDI/;
c   s/FRMDRV/SUPDRV/; implemented NOSYM option;
c   s/SUPA/SUP1D/; s/SUPB/SUPD/; implemented OLDSUP option;
c   implemented NEWSUP = .NOT. OLDSUP format (label "PSUPMAT ");
c   removed saving IJ,KL in INDP in SUP1,SUP,SUP1D,SUPD; just
c   allocate IT(NP+1) for INDP when not OLDSUP
c===========================================================================
#endif
C  /* Deck frmsup */
      SUBROUTINE FORMSUP(WORK,LWORK,NOSYM,THRESH,IPRINT)
C
C --- PROGRAM FORMSUP (form super matrices for Fock matrices)
C
C Written by Olav Kvalheim, University of Bergen, Norway.
C
C Revised
C        5-Jul-1984 hjaaj (changed all var. to max 6 char.)
C       17-Nov-1984 hjaaj SUPSRT error corrected: must use IRAT
C                         instead of 2 (because of e.g. FPS)
C       23-Sep-1986 tuh   adaption to HERMIT
C       12-Jun-1987 hjaaj increased NBMAX from 150 to 255 and
C                         made appropriate changes in commons etc.
C       13-Jun-1987 hjaaj changed packing to accept NBMAX .gt. 127
C        (VAX gives integer overflow on e.g. LABEL = NS*2**24, when NS .gt. 127)
C        9-Jan-1990 tuh+hjaaj: NOSYM parameter;
C          Dec 1991 - work space in call, tuh
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      LOGICAL NOSYM
      DIMENSION WORK(LWORK)
C
      CALL QENTER('FORMSUP')
C
      IF (IPRINT .GT. 1) THEN
         CALL TITLER('Output from FORMSUP','*',125)
         WRITE (LUPRI,'(A/A,1P,D10.2/A,L10)')
     &   '    Precalculated two-electron integrals'
     &   //' are transformed to P-supermatrix elements.',
     &   '    Threshold for discarding integrals            :',THRESH,
     &   '    Include integrals for general density matrices:',NOSYM
      ELSE
         WRITE (LUPRI,'(A/A,1P,D10.2,A)')
     &   '     (Precalculated two-electron integrals'
     &   //' are transformed to P-supermatrix elements.',
     &   '      Threshold for discarding integrals :',THRESH,' )'
      END IF
C
C*** Allocate memory for IT vector in IWORK
C
      NBMAX = 255
      LWRK2 = (NBMAX*(NBMAX+1))/2
#if defined (VAR_DUPLCHCK)
C and allocate memory for IBIT vector in IWORK
      LWRK2 = LWRK2 + LWORK / 32
#endif
      LWRK2 = LWRK2 / IRAT + 1
      LWRK1 = LWORK - LWRK2
      KWRK1 = 1
      KWRK2 = KWRK1 + LWRK1
      LIWRK2 = IRAT*LWRK2
      CALL SUPDRV(NOSYM,WORK(KWRK1),LWRK1,WORK(KWRK2),LIWRK2,
     &            THRESH,IPRINT)
      CALL QEXIT('FORMSUP')
      RETURN
      END
C  /* Deck supdrv */
      SUBROUTINE SUPDRV(NOSYM,WORK,LWORK,IWORK,LIWORK,THRESH,
     &                  IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION WORK(LWORK), IWORK(LIWORK)
      LOGICAL   NOSYM, OPTST
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
      COMMON/CSPREA/THRQ,BUF(600),DUM4(4)
#include <inftap.h>
C
C*** DEFAULTS
C
      NWORK = LIWORK
      DST1  = SECOND()
      NBMAX = 255
      NUT   = 0
      THRQ = THRESH
      KSUP = 0
C
      LINT   = 1706
      LDA2   = (IRAT+1)*LINT
      IDABUF = IRAT*LINT
      LDA21  = LDA2+1
      LDA22  = LDA21+1
      LRDA22 = (LDA22-1)/IRAT + 1
      CALL SUPRDI(NOSYM)
C
      NUT  = 0
C
C***  SIZE OF FIRST BLOCK on LUDASP
C     NST = starting NP for second block
C     NST1 = last NP in first block
C     NST1 max. equal to nbas1, ie. to nbas(1) before NOSYM reset to nbast
C
      NST1  = NBAS1
      LAVAIL = LWORK - LRDA22
C     ... max. memory for INDP
      IJ = (NST1*(NST1+1))/2
      IJ = (IJ-1)/IRAT + 1
      LAVAIL = LAVAIL - IJ
C     IJ = nst1*(nst1+1)/2 = # (p q) distributions in first block
C     NPMAX = ij*(ij+1)/2 = # unique integrals in first block
C     To have all integrals of first block in memory NPMAX .le. LAVAIL
C     thus:
C          IJ = int(sqrt(2*lavail+1)-.5) is max IJ which fit in memory
C          NST1 = int(sqrt(2*IJ)+1)-.5) is max NST1 which fit in memory
      RTMP = 2*LAVAIL+1
      RTMP = SQRT(RTMP) - DP5
      IJ   = RTMP
      RTMP = 2*IJ+1
      RTMP = SQRT(RTMP) - DP5
      NST1 = RTMP
      NST1 = MAX(15,NST1)
      NST1 = MIN(NBAS1,NST1)
      NST  = NST1 + 1
C
C*** SET UP DYNAMICAL STORAGE FOR SORT ROUTINE
C
C allocate space for buffer with 600 real + 602 integer elements
      LENBUF  = ((600*IRAT+602)-1)/IRAT + 1
      LW1   = 1
      LW2   = LW1+LENBUF
      LSORT = LWORK-LENBUF
      DST   = SECOND()
      CALL SUP_SORT(NOSYM,WORK(LW1),WORK(LW1),WORK(LW2),WORK(LW2),LSORT,
     &            IWORK(1),IPRINT)
      DFIN  = SECOND()
      DTOT  = DFIN-DST
      IF (IPRINT .GT. 1) WRITE(LUPRI,340) DTOT
 340  FORMAT(/,'    TIME IN SUP_SORT :',F9.2,' SEC.')
C
C*** START TRANSFERRING INTEGRALS TO SUPERMATRIX FORM.
C*** BEGIN WITH FIRST BLOCK CONTAINING ALL INDICES BETWEEN
C*** 1 AND NST-1. ALL INDICES IN FIRST IRREP.
C
C     first NBW elements of IWORK contains IT vector
C     (initialized by SUP_SORT)
C
      NBW1 = NBW+1
C
C*** FIRST FIND NUMBER OF P-SUPER MATRIX ELEMENTS
C
      IBEF   = 0
      NPLAST = 0
      NST1   = NST-1
      IJ     = IWORK(NST1+1)
      NPMAX  = IWORK(IJ+1)
      LW1    = 1
      LW2    = LRDA22+1
      IJ     = (IJ-1)/IRAT + 1
      NPM    = LRDA22 + NPMAX + IJ
      LSUP   = LWORK-LRDA22-IJ
      ICHAIN = 1
      DST    = SECOND()
      JBUF = 256 + 144*IRAT
      IF (NPM.LE.LWORK) THEN
         LW3 = LW2+NPMAX
         CALL SUP1(ICHAIN,NPMAX,WORK(LW1),WORK(LW1),WORK(LW2),
     &             WORK(LW3),IBEF,IWORK(1),IWORK(NBW1),JBUF)
         DFIN = SECOND()
         DTOT = DFIN-DST
         IF (IPRINT .GT. 1) WRITE(LUPRI,350) DTOT
 350     FORMAT(/,'    TIME IN SUP1  :',F9.2,' SEC.')
      ELSE
         LW3 = LW2+LSUP
         CALL SUP1D(ICHAIN,NPMAX,WORK(LW1),WORK(LW1),WORK(LW2),
     &             WORK(LW3),IBEF,IWORK(1),IWORK(NBW1),JBUF)
         DFIN = SECOND()
         DTOT = DFIN-DST
         IF (IPRINT .GT. 1) WRITE(LUPRI,351) DTOT
 351     FORMAT(/,'    TIME IN SUP1D :',F9.2,' SEC.')
      END IF
      IBEF = IBEF+NPMAX
C
C*** CONTINUE WITH REST OF FIRST SYMMETRY
C
      NPLAST = NST1
      IFIRST = NBAS(1)-NPLAST
      IF (IFIRST.EQ.0) GO TO 100
      DST = DFIN
      DTOTD = D0
      DO 60 I = 1,IFIRST
         ICHAIN = ICHAIN+1
         NPMX   = NPMAX
         NP     = NST1+I
         IJ     = IWORK(NP+1)
         NPMAX  = IWORK(IJ+1)
         NPMM   = NPMAX-NPMX
         IJ     = (IJ-1)/IRAT + 1
         NPM    = LRDA22 + NPMM + IJ
         LSUP   = LWORK - LRDA22 - IJ
         IF (NPM.LE.LWORK) THEN
            LW3 = LW2+NPMM
            CALL SUP1(ICHAIN,NPMM,WORK(LW1),WORK(LW1),WORK(LW2),
     &                WORK(LW3),IBEF,IWORK(1),IWORK(NBW1),JBUF)
         ELSE
            DTOTD = DTOTD - SECOND()
            LW3 = LW2+LSUP
            CALL SUP1D(ICHAIN,NPMM,WORK(LW1),WORK(LW1),WORK(LW2),
     &                WORK(LW3),IBEF,IWORK(1),IWORK(NBW1),JBUF)
            DTOTD = DTOTD + SECOND()
         END IF
         IBEF=IBEF+NPMM
         NPLAST = NP
  60  CONTINUE
      DFIN = SECOND()
      DTOT = DFIN-DST - DTOTD
      IF (IPRINT .GT. 1) THEN
         WRITE(LUPRI,350) DTOT
         IF (DTOTD .NE. D0) WRITE(LUPRI,351) DTOTD
      END IF
C
C*** NOW TRANSFER REST OF THE INTEGRALS
C
 100  CONTINUE
      IF(NSYM.EQ.1) GO TO 250
      DST = DFIN
      DTOTD = D0
      NSTA = IBAS(2)+1
      DO 200 NP = NSTA,NBAST
         ICHAIN = ICHAIN+1
C
C*** NUMBER OF SUPER MATRIX ELEMENTS
C
         KK = 0
         ISP = IS(NP)
         NPR = NP-IBAS(ISP)
         ISR = ISP-1
         DO 190 I = 1,ISR
            NBI = NBAS(I)+1
 190        KK = KK+IWORK(NBI)
         I1 = IWORK(NPR+1)
         I2 = IWORK(NPR)
         IJ = KK + I1
         KOMB = I1-I2
         I1 = IWORK(I1+1)
         I2 = IWORK(I2+1)
         NPMAX = KK*KOMB+I1-I2
         IJ = (IJ-1)/IRAT + 1
         NPM = LRDA22+NPMAX+IJ
         LSUP = LWORK - LRDA22 - IJ
         IF (NPM.LE.LWORK) THEN
            LW3 = LW2+NPMAX
            CALL SUP(ICHAIN,NPMAX,WORK(LW1),WORK(LW1),WORK(LW2),
     &               WORK(LW3),IWORK(1),IWORK(NBW1),IBEF,NP,IPRINT,
     &               JBUF)
         ELSE
            DTOTD = DTOTD - SECOND()
            LW3 = LW2+LSUP
            CALL SUPD(ICHAIN,NPMAX,WORK(LW1),WORK(LW1),WORK(LW2),
     &                WORK(LW3),IWORK(1),IWORK(NBW1),IBEF,NP,JBUF)
            DTOTD = DTOTD + SECOND()
         END IF
         IBEF = IBEF+NPMAX
         NPLAST = NP
 200  CONTINUE
      DFIN = SECOND()
      DTOT = DFIN-DST - DTOTD
      IF (IPRINT .GT. 1) THEN
         WRITE(LUPRI,370) DTOT
         IF (DTOTD .NE. D0) WRITE(LUPRI,371) DTOTD
      END IF
 370  FORMAT(/,'    TIME IN SUP   :',F9.2,' SEC.')
 371  FORMAT(/,'    TIME IN SUPD  :',F9.2,' SEC.')
 250  CONTINUE
C
C*** FINALLY BUF MUST BE EMPTIED IF IT STILL CONTAINS INFORMATION
C
      WRITE(LUSUPM) -2,0,0,0,0,0,0
      CALL GPCLOSE(LUDASP,' ')
      DFIN1 = SECOND()
      DTOT1 = DFIN1-DST1
      IF (IPRINT .GT. 1) WRITE(LUPRI,380) DTOT1
 380  FORMAT(/,'    TOTAL TIME    :',F9.2,' SEC.')
      END
C  /* Deck suprdi */
      SUBROUTINE SUPRDI(NOSYM)
#include <implicit.h>
#include <priunit.h>
      LOGICAL NOSYM
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
#include <inftap.h>
      CHARACTER*8 LAB123(3), OLDLBL, NEWLBL
      DATA LAB123 /'********','        ','        '/
      DATA OLDLBL, NEWLBL /'TWOELSUP','PSUPMAT '/
C
      IF (LUONEL .LT. 0)
     &CALL GPOPEN(LUONEL,'AOONEINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUONEL
      READ(LUONEL)
      READ(LUONEL) NSYM,(NBAS(I),I = 1,NSYM),POTNUC
      CALL GPCLOSE(LUONEL,'KEEP')
      NBAST = 0
      DO 10 I = 1,NSYM
  10     NBAST = NBAST+NBAS(I)
      CALL GETDAT(LAB123(2),LAB123(3))
      CALL REWSPL(LUSUPM)
      WRITE (LUSUPM) LAB123, NEWLBL
      WRITE (LUSUPM) NOSYM, NSYM, NBAS
C
C***  If (NOSYM) then reset symmetry index arrays
C
      NBAS1 = NBAS(1)
      IF (NOSYM) THEN
         CALL IZERO(NBAS,8)
         NSYM = 1
         NBAS(1) = NBAST
      END IF
      RETURN
      END
C  /* Deck supsrt */
      SUBROUTINE SUP_SORT(NOSYM,RINT,INT,SCR,ISCR,LSORT,IT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
      COMMON/CSPREA/THRQ,BUF(600),DUM4(4)
      DIMENSION RINT(*),INT(*),SCR(*),ISCR(*),IT(*)
      LOGICAL NOSYM, OLDDX 
      DIMENSION NSOINT(255), NAOS(8)
C
#include <iratdef.h>
#include <inftap.h>
#include <ibtfun.h>
C
      CALL QENTER('SUP_SORT')
C
      IF (LUINTA .LT. 0) CALL GPOPEN(LUINTA,'AOTWOINT',
     &   'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
C
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) MSYM,(NAOS(I),I=1,8),LBUF,NIBUF,NBITS
C
C     Consistency check
C
      NERR = 0
      IF (NSYM .NE. MSYM) NERR = NERR + 1
      DO ISYM = 1,NSYM
         IF (NAOS(ISYM) .NE. NBAS(ISYM)) NERR = NERR + 1
      END DO
      IF (NERR .GT. 0) THEN
         WRITE(LUPRI,'(//A//A,2I5//A/)')
     &      'SUP_SORT ERROR: Info on AOTWOINT inconsistent with input.',
     &      'NSYM from input and from AOTWOINT:',NSYM,MSYM,
     &      'NBAS(*) from input and from AOTWOINT:'
         WRITE(LUPRI,1050)(I,I = 1,NSYM)
         WRITE(LUPRI,1051)(NBAS(I),I = 1,NSYM)
         WRITE(LUPRI,1051)(NAOS(I),I = 1,MSYM)
         CALL QUIT('Info on AOTWOINT inconsistent with input')
      END IF
      IF (LBUF .NE. 600) THEN
        CALL QUIT('Sorry, SUP_SORT only works for LBUF=600 on AOTWOINT')
      END IF
      IF (NIBUF .NE. 1) THEN
        CALL QUIT('Sorry, SUP_SORT only works for NIBUF=1 on AOTWOINT')
      END IF
C
C***  TEST NSYM
C
      IF (NSYM.LT.1 .OR. NSYM.GT.8) THEN
         WRITE(LUPRI,1100) NSYM
 1100    FORMAT(/' SUP_SORT: NO SORT -- WRONG NUMBER OF SYMMETRIES',I3/)
         CALL QUIT('SUP_SORT ERROR: wrong number of symmetries')
      END IF
      IF (IPRINT .GE. 3) THEN
         WRITE(LUPRI,1050)(I,I = 1,NSYM)
         WRITE(LUPRI,1051)(NBAS(I),I = 1,NSYM)
      END IF
 1050 FORMAT(/'  Symmetry                    ',8I10)
 1051 FORMAT( '  Basis functions per symmetry',8I10/)
C
      NN = 0
      DO 1000 I = 1,NSYM
 1000    NN = NN + NBAS(I)
      IF (NN .GT. NBMAX) THEN
         WRITE(LUPRI,1002) NN,NBMAX
 1002    FORMAT(/,' NUMBER OF BASIS FUNCTIONS IS ',I5,' WHICH EXCEEDS',
     &          /,' THE LIMITATION OF ',I5)
         CALL QUIT('SUP_SORT ERROR: max. # basis funcions exceeded')
      END IF
C
C***  COMPUTE VECTOR IT(I)=(I*(I-1))/2
C
      NBW = (NBAST*NBAST+NBAST)/2 + 1
      IF (NBW.GT.NWORK) GO TO 950
      NWW = NWORK-NBW-1
      K = 0
      DO 3 I = 1,NBW
         IT(I) = K
         K     = K+I
   3  CONTINUE
C
C***  SYMMETRY LABELS AND NUMBER OF PREVIOUS BASIS FUNCTIONS
C
      K  = 0
      DO 5 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         DO 4 I = 1,NBASI
            IS(K+I) = ISYM
   4     CONTINUE
         K = K + NBASI
   5  CONTINUE
C
      IBAS(1) = 0
      IIBAS(1)  = 0
      DO 6 ISYM = 2,NSYM
         NBI = NBAS(ISYM-1)
         IBAS(ISYM) = IBAS(ISYM-1) + NBI
         IIBAS(ISYM)  = IIBAS(ISYM-1) + IT(NBI+1)
   6  CONTINUE
C
C***  OPEN UNIT LUDASP FOR SORTED INTEGRALS
C
      LUDASP = -1
      CALL GPOPEN(LUDASP,' ','UNKNOWN','DIRECT',' ',LDA22,OLDDX)
C
C***  NUMBER OF BATCHES
C
      NST1 = NST-1
      NBATCH = NBAST-NST1+1
      DO 10 I = 1,NST1
  10     IBATCH(I) = 1
      K = 1
      DO 11 I = NST,NBAST
         K = K+1
  11     IBATCH(I) = K
C
C***  PRESET LASTAD
C
      DO 12 I = 1,NBMAX
         LASTAD(I) = -1
         NREC(I) = 0
         NSOINT(I) = 0
  12  CONTINUE
      IINT = IRAT*LINT
C
C***  NUMBER OF BATCHES PER PASS
C
      NBP = LSORT/(LDA22/IRAT)
      IF(NBP.GT.NBATCH) NBP = NBATCH
C
C***  NUMBER OF STEPS
C
      NSTEP = (NBATCH+NBP-1)/NBP
      IF (IPRINT .GE. 3) WRITE(LUPRI,1200) NBATCH,LDA22,LSORT,NSTEP,NBP
 1200 FORMAT(/,' NBATCH:',I5,' LDA22:',I5,' LSORT:',I10,' NSTEP:',I5,
     &       /,' BATCHES PER PASS:',I5)
C
C***  BEGIN LOOP OVER INTEGRAL FILE
C
      NINTAO = 0
      NINTSO = 0
      NCHAIN = 0
      IDISK  = 1
C
      NBP2 = 0
      DO 50 ISTEP = 1,NSTEP
      NBP1 = NBP2+1
      NBP2 = NBP2+NBP
      IF(NBP2.GT.NBATCH) NBP2 = NBATCH
C
C***  INITIALIZE THE BATCHES
C
      KBATCH = 0
      DO 15 I = 1,NBP
         ISCR(KBATCH+LDA22-1) = -1
         ISCR(KBATCH+LDA22)   =  0
         KBATCH = KBATCH+LDA22
  15  CONTINUE
Cha
      IAOINT = IRAT*600
      LAOINT = IAOINT+600+1
C
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
      JBUF = 608
C
  20  CONTINUE
#if defined (VAR_SPLITFILES)
      CALL READSI(LUINTA,LAOINT,INT,JBUF)
#else
      CALL READI(LUINTA,LAOINT,INT)
#endif
      NINT = INT(LAOINT)
      IF(NINT.EQ. 0) GO TO 20
      IF(NINT.EQ.-1) GO TO 35
C
C***  LOOP OVER INTEGRALS IN THIS BUFFER
C
      IBITMASK  = 2**NBITS - 1
C     IBITSHR_4 =       0 ! shift for index 4 in (12|34) -- not needed!
      IBITSHR_3 =   NBITS ! shift for index 3 in (12|34)
      IBITSHR_2 = 2*NBITS ! shift for index 2 in (12|34)
      IBITSHR_1 = 3*NBITS ! shift for index 1 in (12|34)
      DO 30 I = 1,NINT
         LABEL  = INT(IAOINT+I)
         VALUE  = RINT(I)
         NINTAO = NINTAO+1
      IF(ABS(VALUE).LE.THRQ) GO TO 30
C
         NP = IBTAND(IBTSHR(LABEL,IBITSHR_1),IBITMASK)
C        HERMIT: We know NP is largest index, thus we don't need NQ,NR,NS
C                unless we check symmetry below
C
C***  ONLY INTEGRALS WITH AT MOST TWO DIFFERENT SYMMETRY INDICES
C***  ARE TO BE KEPT
C
         IF (NOSYM) GO TO 25
         ISP = IS(NP)
         NQ = IBTAND(IBTSHR(LABEL,IBITSHR_2),IBITMASK)
         IF (IS(NQ).EQ.ISP) GO TO 25
         NR = IBTAND(IBTSHR(LABEL,IBITSHR_3),IBITMASK)
         IF (IS(NR).EQ.ISP) GO TO 25
         NS = IBTAND(       LABEL,    IBITMASK)
         IF (IS(NS).EQ.ISP) GO TO 25
      GO TO 30
  25  CONTINUE
C
C        INTRODUCE CANONICAL ORDERING OF INTEGRAL INDICES
C
C        940819-hjaaj: HERMIT already has canonical ordering.
C
         KBATCH = IBATCH(NP)
      IF (KBATCH.LT.NBP1.OR.KBATCH.GT.NBP2) GO TO 30
         KBP = KBATCH-NBP1
         NINTSO = NINTSO+1
C
C***  ALLOCATE INTEGRAL AND LABEL TO BATCH
C
         NSOINT(KBATCH) = NSOINT(KBATCH)+1
         IPOS   = LDA22*KBP
         LENGTH = ISCR(IPOS+LDA22)
         LENGTH = LENGTH+1
         IRPOS  = IPOS/IRAT
         SCR(IRPOS+LENGTH) = VALUE
         ISCR(IPOS+IINT+LENGTH) = LABEL
         ISCR(IPOS+LDA22) = LENGTH
C
         IF (LENGTH.GE.LINT) THEN
C           ... THIS BATCH IS NOW FULL AND MUST BE EMPTIED
            NCHAIN = NCHAIN+1
            IDO = IDISK
            NREC(KBATCH) = NREC(KBATCH)+1
            CALL WRITDX(LUDASP,IDISK,LDA22,ISCR(IPOS+1))
            ISCR(IPOS+LDA22-1) = IDO
            ISCR(IPOS+LDA22)   = 0
            IDISK = IDISK+1
         END IF
C
C     THIS COMPLETES THE LOOP OVER THIS BUFFER OF AO INTEGRALS
C
   30 CONTINUE
C
C     GO BACK FOR THE NEXT BUFFER
C
      GO TO 20
C
C     THE AO FILE IS READ
C
   35 CONTINUE
C
C     FINALLY BATCHES STILL CONTAINING INFORMATION MUST BE EMPTIED
C
      IPOS = -LDA22
      DO 45 I = NBP1,NBP2
         IPOS = IPOS+LDA22
         IDO  = IDISK
      IF (ISCR(IPOS+LDA22).EQ.0) GO TO 40
         NCHAIN = NCHAIN+1
         NREC(I) = NREC(I)+1
         CALL WRITDX(LUDASP,IDISK,LDA22,ISCR(IPOS+1))
         LASTAD(I) = IDO
         IDISK = IDISK+1
      GO TO 45
   40 CONTINUE
         IDO = ISCR(IPOS+LDA22-1)
         LASTAD(I) = IDO
   45 CONTINUE
C
C     END OF LOOP OVER INTEGRAL FILE
C
   50 CONTINUE
C
      IF (IPRINT .GE. 2) WRITE(LUPRI,1250) NINTAO,NINTSO,NCHAIN
 1250 FORMAT(//10X,'NUMBER OF INTEGRALS READ',I12,
     *        /10X,'NUMBER OF INTEGRALS SORTED',I10,
     *        /10X,'NUMBER OF BUFFERS ON LUDASP',I9)
      IF (IPRINT .GE. 3) THEN
         WRITE(LUPRI,1260) (NREC(I),I = 1,NBATCH)
         WRITE(LUPRI,1261) (NSOINT(I),I = 1,NBATCH)
      END IF
 1260 FORMAT(//10X,'NUMBER OF RECORDS IN EACH CHAIN'/(5X,8I9))
 1261 FORMAT(/10X,'NUMBER OF INTEGRALS IN EACH CHAIN'/(5X,8I9))
      CALL GPCLOSE(LUINTA,'KEEP')
      CALL QEXIT('SUP_SORT')
      RETURN
C
C     ERROR MESSAGE
C
 950  WRITE(LUPRI,9500) NBW
 9500 FORMAT(' SUP_SORT: DIMENSION OF IWORK MUST BE INCREASED TO',I10,
     &      /' PLUS ALLOWANCE FOR BITVECTOR IN SUP AND SUP1')
      CALL QUIT('SUP_SORT ERROR: insufficient memory')
      END
C  /* Deck sup1 */
      SUBROUTINE SUP1(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IBEF,IT,IBIT,
     &                JBUF)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
C***  TRANSFORMATION OF FIRST SYMMETRY
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
      COMMON/CSPREA/THRQ,BUF(600),DUM4(4)
#include <inftap.h>
      DIMENSION P(*),INDP(*),IBUF(*),BBUF(*),IBIT(*),IT(*)
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
C
#include <ibtfun.h>
C
#if defined (VAR_DUPLCHCK)
      JO(L) = IBTAND(IBTSHR(IBIT((L+31)/32),((L+31)/32*32-L)),1)
#endif
C
C
C***  PREPARE TEMPORARY STORAGE OF P-SUPERMATRIX
C
      DO 10 I = 1,NPMAX
  10     P(I) = ZERO
#if defined (VAR_DUPLCHCK)
C
C***  PREPARE BIT VECTOR FOR CHECK OF UNIQUE INTEGRALS
C
      KBMAX = NPMAX/32+1
      IF (KBMAX.GT.NWW) GO TO 950
      DO 20 I = 1,KBMAX
  20     IBIT(I) = 0
#endif
C
      IADR = LASTAD(ICHAIN)
      IF (IADR.EQ.-1) GO TO 900
  50  CALL READDX(LUDASP,IADR,LDA22,IBUF)
      IADR = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF (LENGTH.EQ.0) GO TO 500
      DO 400 I = 1,LENGTH
         LDAI  = IDABUF+I
         IBL   = IBUF(LDAI)
         NP    = IBTAND(IBTSHR(IBL,IBITSHR_1),IBITMASK)
         NQ    = IBTAND(IBTSHR(IBL,IBITSHR_2),IBITMASK)
         NR    = IBTAND(IBTSHR(IBL,IBITSHR_3),IBITMASK)
         NS    = IBTAND(       IBL           ,IBITMASK)
         IJ    = IT(NP)+NQ
         KL    = IT(NR)+NS
         IJKL  = IT(IJ)+KL-IBEF
C
C***  ALL INDICES ARE IN IRREP 1 AND IN CANONICAL ORDER
#if defined (VAR_DUPLCHCK)
C***  FIRST CHECK INTEGRAL HASN'T OCCURRED BEFORE
C
      IF (JO(IJKL).EQ.1) GO TO 400
         IK32  = (IJKL+31)/32
         IK322 = IK32*32-IJKL
         IBIT(IK32) = IBTOR(IBIT(IK32),IBTSHL(1,IK322))
#endif
C
C***  CONTRIBUTIONS TO P(NP,NQ,NR,NS),P(NP,NR,NQ,NS) AND
C***  P(NP,NS,NQ,NR). CHECK ORDERING OF NQ-NS AND NQ-NR.
C***  FIRST FIND NP,NQ,NR,NS INDEX
C
      IF (NP.EQ.NQ) GO TO 180
      IF (NR.EQ.NP) GO TO 140
      IF (NR.NE.NQ) GO TO 80
      IF (NS.EQ.NR) GO TO 60
C
C***  PQ,QS
C
      P   (IJKL) = P(IJKL) + TH4*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
      IJ = IT(NP)+NS
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  CONTINUE
      P   (IJKL) = ON2*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      P(IJKL)    = P(IJKL)+BBUF(I)
      IJ = IT(NP)+NR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      P   (IJKL) = P(IJKL)+BBUF(I)
      IJ = IT(NP)+NR
      IF (NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NS)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON4*BBUF(I)
 110  IJ = IT(NP)+NS
      IF(NQ.GT.NR) GO TO 115
C
C***  CONTRIBUTION TO PS,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  CONTINUE
      P   (IJKL) = P(IJKL)+TH4*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
      IJ = IT(NP)+NR
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      P   (IJKL) = P(IJKL)+TH4*BBUF(I)
      IJ = IT(NP)+NP
C
C***  CONTRIBUTION TO PP,QS   CANONICAL ORDERING ASSUMED
C
      KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  CONTINUE
      P   (IJKL) = P(IJKL)+TH8*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
      IJ = IT(NP)+NP
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      P   (IJKL) = P(IJKL)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
      IJ = IT(NP)+NR
      KL = IT(NP)+NS
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  CONTINUE
      P   (IJKL) = P(IJKL)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
      IJ = IT(NP)+NR
      KL = IJ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-ON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      P   (IJKL) = ON2*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  CONTINUE
      P   (IJKL) = ON4*BBUF(I)
 400  CONTINUE
 500  IF(IADR.NE.-1) GO TO 50
C
C***  Transfer P to disk
C
      CALL SUPNWR(P,INDP,IT,IBEF,NPMAX,JBUF)
C
 900  CONTINUE
      RETURN
#if defined (VAR_DUPLCHCK)
 950  NWT = NBW+KBMAX
      WRITE(LUPRI,960) NWT
 960  FORMAT(' SUP1: DIMENSION OF IWORK MUST BE INCREASED TO',I10)
      CALL QUIT('Insfspce')
#endif
      END
C  /* Deck sup */
      SUBROUTINE SUP(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IT,IBIT,
     &               IBEF,NP,IPRINT,JBUF)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
#include <inftap.h>
      DIMENSION P(*),INDP(*),IBUF(*),BBUF(*),IT(*),IBIT(*)
C
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
      COMMON/CSPREA/THRQ,BUF(600),DUM4(4)
C
#include <ibtfun.h>
C
#if defined (VAR_DUPLCHCK)
      JO(L) = IBTAND(IBTSHR(IBIT((L+31)/32),((L+31)/32*32-L)),1)
#endif
C
C***  PREPARE TEMPORARY STORAGE OF P-SUPERMATRIX
C
      DO 10 I = 1,NPMAX
  10     P(I) = ZERO
#if defined (VAR_DUPLCHCK)
C
C***  PREPARE BITVECTOR
C
      II    = IT(NP+1)
      ITP   = IT(NP)
      IMAX  = IT(II+1)
      JMAX  = IT(ITP+1)
      KMAX  = IMAX-JMAX
      KBMAX = KMAX/32+1
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI,'(/A)') ' *** test output from SUP'
         WRITE (LUPRI,'(/A,I15)') ' >>> NP    ',NP
         WRITE (LUPRI,'( A,I15)') ' >>> NPMAX ',NPMAX
         WRITE (LUPRI,'( A,I15)') ' >>> II    ',II
         WRITE (LUPRI,'( A,I15)') ' >>> ITP   ',ITP
         WRITE (LUPRI,'( A,I15)') ' >>> IMAX  ',IMAX
         WRITE (LUPRI,'( A,I15)') ' >>> JMAX  ',JMAX
         WRITE (LUPRI,'( A,I15)') ' >>> KMAX  ',KMAX
         WRITE (LUPRI,'( A,I15)') ' >>> KBMAX ',KBMAX
         WRITE (LUPRI,'( A,I15)') ' >>> NWW   ',NWW
      END IF
      IF(KBMAX.GT.NWW) GO TO 950
      DO 15 I = 1,KBMAX
  15     IBIT(I) = 0
#endif
C
C*** NP-INDEX HANDLING OUTSIDE OF LOOP
C
      ISP = IS(NP)
      NPR = NP-IBAS(ISP)
      MPP = IT(NPR)+IIBAS(ISP)
C
      IADR = LASTAD(ICHAIN)
      IF(IADR.EQ.-1) GO TO 900
  50  CALL READDX(LUDASP,IADR,LDA22,IBUF)
      IADR   = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF(LENGTH.EQ.0) GO TO 410
      DO 400 I = 1,LENGTH
      LDAI = IDABUF+I
      IBL  = IBUF(LDAI)
C
      NQ = IBTAND(IBTSHR(IBL,IBITSHR_2),IBITMASK)
      NR = IBTAND(IBTSHR(IBL,IBITSHR_3),IBITMASK)
      NS = IBTAND(       IBL,           IBITMASK)
#if defined (VAR_DUPLCHCK)
C
C***  FIRST CHECK IF INTEGRAL HAS OCCURRED BEFORE
C
      IJ   = ITP+NQ
      KL   = IT(NR)+NS
      IJKL = IT(IJ)+KL-JMAX
      IF(JO(IJKL).EQ.1) GO TO 400
      IK31 = (IJKL+31)/32
      IK32 = IK31*32-IJKL
      IBIT(IK31) = IBTOR(IBIT(IK31),IBTSHL(1,IK32))
#endif
C
      ISQ = IS(NQ)
      ISR = IS(NR)
      ISS = IS(NS)
C
      NQR = NQ-IBAS(ISQ)
      NRR = NR-IBAS(ISR)
      NSR = NS-IBAS(ISS)
      IF(ISP.NE.ISQ) GO TO 250
      IF(ISP.NE.ISR) GO TO 350
C
C***  ISP = ISQ = ISR = ISS
C
      IJ = MPP+NQR
      KL = IT(NRR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(NP.EQ.NQ) GO TO 180
      IF(NR.EQ.NP) GO TO 140
      IF(NR.NE.NQ) GO TO 80
      IF(NS.EQ.NR) GO TO 60
C
C***  PQ,QS
C
      P   (ID) = P(ID)+TH4*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
      IJ = MPP+NSR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  CONTINUE
      P(ID)    = ON2*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      P(ID)    = P(ID)+BBUF(I)
      IJ       = MPP+NRR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NRR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      P(ID)    = P(ID)+BBUF(I)
      IJ = MPP+NRR
      IF(NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NSR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON4*BBUF(I)
 110  CONTINUE
      IJ = MPP+NSR
      IF(NQ.GT.NR) GO TO 115
C
C***  CONTRIBUTION TO PS,RQ
C
      KL = IT(NRR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  CONTINUE
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  CONTINUE
      P(ID)    = P(ID)+TH4*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
      IJ = MPP+NRR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      P(ID)    = P(ID)+TH4*BBUF(I)
      IJ = MPP+NPR
C
C***  CONTRIBUTION TO PP,QS  CANONICAL ORDERING ASSUMED
C
      KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  CONTINUE
      P(ID)    = P(ID)+TH8*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
      IJ = MPP+NPR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      P(ID)    = P(ID)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
      IJ = MPP+NRR
      KL = IT(NPR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  CONTINUE
      P(ID)    = P(ID)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
      IJ = MPP+NRR
      KL = IJ
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-ON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      P(ID)    = ON2*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  CONTINUE
      P(ID)    = ON4*BBUF(I)
      GO TO 400
C
C***  ISP = ISR, ISP>ISQ, ISQ = ISS
C
C***  P(NP,NR,NQ,NS)
C
 250  IF(NQ.GE.NS) GO TO 252
      NX  = NQR
      NQR = NSR
      NSR = NX
 252  IJ  = MPP+NRR
      KL  = IT(NQR)+NSR+IIBAS(ISQ)
      ID  = IT(IJ)+KL-IBEF
      IF(NQ.EQ.NS) GO TO 264
      IF(NP.EQ.NR) GO TO 264
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
 264  P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  ISP = ISQ, ISP>ISR, ISR = ISS
C
 350  IJ = MPP+NQR
      KL = IT(NRR)+NSR+IIBAS(ISR)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)+BBUF(I)
 400  CONTINUE
C
 410  IF(IADR.NE.-1) GO TO 50
C
C***  Transfer P to disk
C
      CALL SUPNWR(P,INDP,IT,IBEF,NPMAX,JBUF)
C
 900  CONTINUE
      RETURN
#if defined (VAR_DUPLCHCK)
 950  KK = NBW+KBMAX
      WRITE(LUPRI,960) KK
 960  FORMAT(' SUP: DIMENSION OF IWORK MUST BE INCREASED TO',I10)
      CALL QUIT('Insfspce')
#endif
      END
C  /* Deck sup1d */
      SUBROUTINE SUP1D(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IBEF,
     &                 IT,IBIT,JBUF)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
C***  TRANSFORMATION OF FIRST SYMMETRY
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
      COMMON/CSPREA/THRQ,BUF(600),DUM4(4)
#include <inftap.h>
      DIMENSION P(*),INDP(*),IBUF(*),BBUF(*),IBIT(*),IT(*)
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
C
#include <ibtfun.h>
C
#if defined (VAR_DUPLCHCK)
      JO(L) = IBTAND(IBTSHR(IBIT((L+31)/32),((L+31)/32*32-L)),1)
#endif
C
C*** FIND NUMBER OF STEPS
C
      NSTEP = NPMAX/LSUP+1
      NB2 = 0
      IRED = 0
      DO 800 II = 1,NSTEP
      NB1 = NB2+1
      NB2 = NB1+LSUP-1
      IF(NB2.GT.NPMAX) NB2 = NPMAX
      NBB = NB2-NB1+1
C
C***  PREPARE TEMPORARY STORAGE OF P-SUPERMATRIX
C
      DO 10 I = 1,NBB
  10  P(I) = ZERO
#if defined (VAR_DUPLCHCK)
C
C***  PREPARE BIT VECTOR FOR CHECK OF UNIQUE INTEGRALS
C
      KBMAX = NPMAX/32+1
      IF(KBMAX.GT.NWW) GO TO 950
      DO 20 I = 1,KBMAX
  20  IBIT(I) = 0
#endif
C
      IADR = LASTAD(ICHAIN)
      IF(IADR.EQ.-1) GO TO 900
  50  CALL READDX(LUDASP,IADR,LDA22,IBUF)
      IADR = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF(LENGTH.EQ.0) GO TO 500
      DO 400 I = 1,LENGTH
      LDAI = IDABUF+I
      IBL = IBUF(LDAI)
      NP = IBTAND(IBTSHR(IBL,IBITSHR_1),IBITMASK)
      NQ = IBTAND(IBTSHR(IBL,IBITSHR_2),IBITMASK)
      NR = IBTAND(IBTSHR(IBL,IBITSHR_3),IBITMASK)
      NS = IBTAND(       IBL,           IBITMASK)
      IJ    = IT(NP)+NQ
      KL    = IT(NR)+NS
      IJKL  = IT(IJ)+KL-IBEF
      ID    = IJKL-IRED
C
C***  ALL INDICES ARE IN IRREP 1 AND IN CANONICAL ORDER
#if defined (VAR_DUPLCHCK)
C***  FIRST CHECK INTEGRAL HASN'T OCCURRED BEFORE
C
      IF(JO(IJKL).EQ.1) GO TO 400
      IK32  = (IJKL+31)/32
      IK322 = IK32*32-IJKL
      IBIT(IK32) = IBTOR(IBIT(IK32),IBTSHL(1,IK322))
#endif
C
C***  CONTRIBUTIONS TO P(NP,NQ,NR,NS),P(NP,NR,NQ,NS) AND
C***  P(NP,NS,NQ,NR). CHECK ORDERING OF NQ-NS AND NQ-NR.
C***  FIRST FIND NP,NQ,NR,NS INDEX
C
      IF(NP.EQ.NQ) GO TO 180
      IF(NR.EQ.NP) GO TO 140
      IF(NR.NE.NQ) GO TO 80
      IF(NS.EQ.NR) GO TO 60
C
C***  PQ,QS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 55
      P(ID) = P(ID)+TH4*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
  55  IJ = IT(NP)+NS
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      P(ID) = ON2*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 82
      P(ID) = P(ID)+BBUF(I)
  82  IJ = IT(NP)+NR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 101
      P(ID) = P(ID)+BBUF(I)
 101  IJ = IT(NP)+NR
      IF(NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NS)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 110
      ID = IJKL-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 110
      ID = IJKL-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
 110  IJ = IT(NP)+NS
      IF(NQ.GT.NR) GO TO 115
C
C***  CONTRIBUTION TO PS,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 125
      P(ID) = P(ID)+TH4*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
 125  IJ = IT(NP)+NR
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 142
      P(ID) = P(ID)+TH4*BBUF(I)
 142  IJ = IT(NP)+NP
C
C***  CONTRIBUTION TO PP,QS   CANONICAL ORDERING ASSUMED
C
      KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 162
      P(ID) = P(ID)+TH8*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
 162  IJ = IT(NP)+NP
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 182
      P(ID) = P(ID)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
 182  IJ = IT(NP)+NR
      KL = IT(NP)+NS
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 202
      P(ID) = P(ID)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
 202  IJ = IT(NP)+NR
      KL = IJ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-ON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      P(ID) = ON2*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      P(ID) = ON4*BBUF(I)
 400  CONTINUE
 500  IF(IADR.NE.-1) GO TO 50
C
C
C
C***  Transfer P to disk
C
      IOFF = IBEF + IRED
C     Note: A PQ distribution may be divided into two records
      CALL SUPNWR(P,INDP,IT,IOFF,NBB,JBUF)
C
C
 800  IRED = IRED+NBB
 900  CONTINUE
      RETURN
#if defined (VAR_DUPLCHCK)
 950  NWT = NBW+KBMAX
      WRITE(LUPRI,960) NWT
 960  FORMAT(' SUP1D: DIMENSION OF IWORK MUST BE INCREASED TO',I10)
      CALL QUIT('Insfspce')
#endif
      END
C  /* Deck supd */
      SUBROUTINE SUPD(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IT,IBIT,
     &               IBEF,NP,JBUF)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
      COMMON/CSPREA/THRQ,BUF(600),DUM4(4)
#include <inftap.h>
      DIMENSION P(*),INDP(*),IBUF(*),BBUF(*),IT(*),IBIT(*)
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
C
#include <ibtfun.h>
C
#if defined (VAR_DUPLCHCK)
      JO(L) = IBTAND(IBTSHR(IBIT((L+31)/32),((L+31)/32*32-L)),1)
#endif
C
C*** FIND NUMBER OF STEPS
C
      NSTEP = NPMAX/LSUP+1
      NB2 = 0
      IRED = 0
      DO 800 III = 1,NSTEP
      NB1 = NB2+1
      NB2 = NB1+LSUP-1
      IF(NB2.GT.NPMAX) NB2 = NPMAX
      NBB = NB2-NB1+1
C
C***  PREPARE TEMPORARY STORAGE OF P-SUPERMATRIX
C
      DO 10 I = 1,NBB
  10  P(I) = ZERO
#if defined (VAR_DUPLCHCK)
C
C***  PREPARE BITVECTOR
C
      II = IT(NP+1)
      ITP = IT(NP)
      IMAX = IT(II+1)
      JMAX = IT(ITP+1)
      KMAX = IMAX-JMAX
      KBMAX = KMAX/32+1
      IF(KBMAX.GT.NWW) GO TO 950
      DO 15 I = 1,KBMAX
  15     IBIT(I) = 0
#endif
C
C*** NP-INDEX HANDLING OUTSIDE OF LOOP
C
      ISP = IS(NP)
      NPR = NP-IBAS(ISP)
      MPP = IT(NPR)+IIBAS(ISP)
C
      IADR = LASTAD(ICHAIN)
      IF(IADR.EQ.-1) GO TO 900
  50  CALL READDX(LUDASP,IADR,LDA22,IBUF)
      IADR = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF(LENGTH.EQ.0) GO TO 410
      DO 400 I = 1,LENGTH
      LDAI = IDABUF+I
      IBL = IBUF(LDAI)
C
      NQ = IBTAND(IBTSHR(IBL,IBITSHR_2),IBITMASK)
      NR = IBTAND(IBTSHR(IBL,IBITSHR_3),IBITMASK)
      NS = IBTAND(       IBL,           IBITMASK)
#if defined (VAR_DUPLCHCK)
C
C***  FIRST CHECK IF INTEGRAL HAS OCCURRED BEFORE
C
      IJ = ITP+NQ
      KL = IT(NR)+NS
      IJKL = IT(IJ)+KL-JMAX
      IF(JO(IJKL).EQ.1) GO TO 400
      IK31 = (IJKL+31)/32
      IK32 = IK31*32-IJKL
      IBIT(IK31) = IBTOR(IBIT(IK31),IBTSHL(1,IK32))
#endif
C
      ISQ = IS(NQ)
      ISR = IS(NR)
      ISS = IS(NS)
C
      NQR = NQ-IBAS(ISQ)
      NRR = NR-IBAS(ISR)
      NSR = NS-IBAS(ISS)
      IF(ISP.NE.ISQ) GO TO 250
      IF(ISP.NE.ISR) GO TO 350
C
C***  ISP = ISQ = ISR = ISS
C
      IJ = MPP+NQR
      KL = IT(NRR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(NP.EQ.NQ) GO TO 180
      IF(NR.EQ.NP) GO TO 140
      IF(NR.NE.NQ) GO TO 80
      IF(NS.EQ.NR) GO TO 60
C
C***  PQ,QS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 55
      ID = ID-IRED
      P(ID) = P(ID)+TH4*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
  55  IJ = MPP+NSR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = ON2*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 82
      ID = ID-IRED
      P(ID) = P(ID)+BBUF(I)
  82  IJ = MPP+NRR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NRR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 101
      ID = ID-IRED
      P(ID) = P(ID)+BBUF(I)
 101  IJ = MPP+NRR
      IF(NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NSR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 110
      ID = ID-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 110
      ID = ID-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
 110  CONTINUE
      IJ = MPP+NSR
      IF(NQ.GT.NR) GO TO 115
C
C***  CONTRIBUTION TO PS,RQ
C
      KL = IT(NRR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  CONTINUE
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 122
      ID = ID-IRED
      P(ID) = P(ID)+TH4*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
 122  IJ = MPP+NRR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 142
      ID = ID-IRED
      P(ID) = P(ID)+TH4*BBUF(I)
 142  IJ = MPP+NPR
C
C***  CONTRIBUTION TO PP,QS  CANONICAL ORDERING ASSUMED
C
      KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 162
      ID = ID-IRED
      P(ID) = P(ID)+TH8*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
 162  IJ = MPP+NPR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 182
      ID = ID-IRED
      P(ID) = P(ID)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
 182  IJ = MPP+NRR
      KL = IT(NPR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 202
      ID = ID-IRED
      P(ID) = P(ID)+BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
 202  IJ = MPP+NRR
      KL = IJ
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-ON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = ON2*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = ON4*BBUF(I)
      GO TO 400
C
C***  ISP = ISR, ISP>ISQ, ISQ = ISS
C
C***  P(NP,NR,NQ,NS)
C
 250  IF(NQ.GE.NS) GO TO 252
      NX = NQR
      NQR = NSR
      NSR = NX
 252  IJ = MPP+NRR
      KL = IT(NQR)+NSR+IIBAS(ISQ)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      IF(NQ.EQ.NS) GO TO 264
      IF(NP.EQ.NR) GO TO 264
      P(ID) = P(ID)-ON4*BBUF(I)
      GO TO 400
 264  P(ID) = P(ID)-ON2*BBUF(I)
      GO TO 400
C
C***  ISP = ISQ, ISP>ISR, ISR = ISS
C
 350  IJ = MPP+NQR
      KL = IT(NRR)+NSR+IIBAS(ISR)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)+BBUF(I)
 400  CONTINUE
C
 410  IF(IADR.NE.-1) GO TO 50
C
C
C***  Transfer P to disk
C
      IOFF = IBEF + IRED
C     Note: A PQ distribution may be divided into two records
      CALL SUPNWR(P,INDP,IT,IOFF,NBB,JBUF)
C
C
 800  IRED = IRED+NBB
 900  CONTINUE
      RETURN
#if defined (VAR_DUPLCHCK)
 950  KK = NBW+KBMAX
      WRITE(LUPRI,960) KK
 960  FORMAT('SUPD: DIMENSION OF IWORK MUST BE INCREASED TO',I10)
      CALL QUIT('Insfspce')
#endif
      END
C  /* Deck supnwr */
      SUBROUTINE SUPNWR(P,INDP,IT,IBEF,NPMAX,JBUF)
C
C     Jan 90 Hans Joergen Aa. Jensen
C     Write to LUSUPM
C
C     On input P contains P(IBEF+1:IBEF+NPMAX)
C
C     WRITE (LUSUPM) ITYP,NPF,NQF,NPL,NQL,IBEF,NPQRS
C     Two format types:
C     ITYP=1: WRITE (LUSUPM) P(IBEF+1:IBEF+NPQRS)
C     ITYP=2: WRITE (LUSUPM) P(1:NPQRS),INDP(1:NPQRS)
C
C     note: IBEF=IT(NPQF),       NPQF = IT(NPF) + NQF
C           NPQRS=IT(NPQL)-IOFF, NPQL = IT(NPL) + NQL
C
#include <implicit.h>
      DIMENSION P(1), INDP(1), IT(1)
C
#include <iratdef.h>
      PARAMETER ( ITYP1 = 1, ITYP2 = 2)
C
#include <inftap.h>
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,NBMAX,NWORK,IIBAS(8),
     &      LBUF,NIBUF,NBITS,IBITSHR_1,IBITSHR_2,IBITSHR_3, IBITMASK,
     A      NBAS1,NPLAST,NST,IND(600),NUT,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LASTAD(255),
     C      NREC(255),IBATCH(255),IS(255),NBW,NWW,LSUP
      COMMON/CSPREA/THRQ,BUF(600),DUM4(4)
C
      IOFF  = IBEF
      IF (NPLAST .EQ. 0) THEN
         NPL   = NPLAST
         NPLAST= MIN(NST-1,15)
C        if NPLAST = 15 then NPQL = 120 and NPQRS = 7260
C
  100    NPF   = NPL+1
         NQF   = 1
         NPL   = NPLAST
         NQL   = NPL
         NPQF  = IT(NPF)+NQF
         NPQL  = IT(NPL)+NQL
         NPQRS = IT(NPQL+1) - IOFF
         IF ( (IOFF + NPQRS) .GT. (IBEF+NPMAX) ) THEN
C           ... all elements not available, switch to type 2
            NPLAST = NPF
            GO TO 200
         END IF
         IF ( NPF .GT. 1 .AND.
     &        NDXGTA(NPQRS,THRQ,P(IOFF+1),1) .LE. NPQRS/2 ) THEN
C           ... less than half the elements .ge. THRQ,
C               switch to type 2
            NPLAST = NPF
            GO TO 200
         END IF
         CALL WRSUP1(LUSUPM,P(IOFF+1),NPQRS,ITYP1,NPF,NQF,NPL,
     &               NQL,IOFF,NPQL,JBUF)
         IOFF  = IOFF + NPQRS
         IF (NPLAST + 1 .LE. NST-1) THEN
            NPLAST = NPLAST + 1
            GO TO 100
         END IF
      END IF
  200 IF (IOFF .EQ. (IBEF+NPMAX)) GO TO 900
C
C
      NPF = NPLAST
      DO 740 NP = NPF,NBAST
C     ... we begin with NPLAST, because NPLAST may not be finished
C         if SUP1D or SUPD was calling routine.
      ISP = IS(NP)
      NPR = NP - IBAS(ISP)
      DO 730 NQR = 1,NPR
         NPQR   = IIBAS(ISP) + IT(NPR) + NQR
         NPQROW = IT(NPQR)
       IF (NPQROW+NPQR .LE. IOFF) GO TO 730
       IF (NPQROW .GE. (IBEF+NPMAX)) GO TO 900
C      (  Finished with this block of integrals )
         NRSOUT = 0
         NPQRS  = NPQROW
C        ... to avoid compiler messages
         DO 720 NR = 1,NP
         ISR = IS(NR)
         NRR = NR - IBAS(ISR)
         INDPR = IT(NR) + IBAS(ISR)
         NSEND  = NRR
         IF (NR .EQ. NP) NSEND = NQR
         DO 710 NSR = 1,NSEND
            NRSR = IIBAS(ISR) + IT(NRR) + NSR
            NPQRS = NPQROW + NRSR
          IF (NPQRS .LE. IOFF) GO TO 710
          IF (NPQRS .GT. (IBEF+NPMAX)) GO TO 721
          IF(ABS(P(NPQRS-IBEF)).GT.THRQ) THEN
            NRSOUT = NRSOUT + 1
            P(NRSOUT) = P(NPQRS-IBEF)
            INDP(NRSOUT) = INDPR + NSR
          END IF
  710    CONTINUE
  720    CONTINUE
  721    CONTINUE
         IF (NRSOUT .GT. 0) THEN
            NQ   = IBAS(ISP) + NQR
            NPQ  = IT(NP) + NQ
            IOFPQ= IT(NPQ)
            CALL WRSUP2(LUSUPM,P,INDP,ITYP2,NP,NQ,
     &                  IOFPQ,NRSOUT,NPQ,JBUF)
         END IF
         IF (NPQRS .GT. (IBEF+NPMAX)) GO TO 900
  730 CONTINUE
      NPLAST = NP
  740 CONTINUE
C
  900 CONTINUE
      RETURN
      END
C  /* Deck wrsup1 */
      SUBROUTINE WRSUP1(LUSUPM,P,NPQRS,ITYP1,NPF,NQF,NPL,NQL,IOFF,
     &                  NPQL,JBUF)
#include <implicit.h>
#include <priunit.h>
      DIMENSION P(NPQRS)
#if defined (VAR_SPLITFILES)
#include <dummy.h>
#include <iratdef.h>
#include <2gbdef.h>
      CHARACTER*80 FNNAME, FNNM2
#include <chrnos.h>
C
      IPASS = 1
 5    CONTINUE
      IF ((((JBUF + NPQRS*IRAT) .GT. I2GB) .AND. (IPASS .EQ. 2)) .OR.
     &     (((JBUF + 8*IRAT) .GT. I2GB) .AND. (IPASS .EQ. 1))) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
         INQUIRE(UNIT=LUSUPM,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         CALL GPCLOSE(LUSUPM,'KEEP')
         I = LN - 1
         IF (FNNAME(I:I) .NE. '-') THEN
            FNNM2 = FNNAME(1:LN)//'-0'
            LN = LN + 2
         ELSE
            READ(FNNAME(LN:),'(I1)') INUM
            INUM = INUM + 1
            IF (INUM .GT. 9) THEN
               WRITE (LUPRI,'(/A)') ' DALTON needs to split '//
     &              ' a file more than 11 times.',
     &              ' This is currently not supported'
               CALL QUIT('Too many splittings of a file')
            END IF
            FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
         END IF
         CALL GPOPEN(LUSUPM,FNNM2(1:LN),'UNKNOWN',' ',' ',
     &        IDUMMY,.FALSE.)
         JBUF = 0
      END IF
      IF (IPASS .EQ. 1) THEN
         WRITE (LUSUPM) ITYP1,NPF,NQF,NPL,NQL,IOFF,NPQRS,NPQL
         JBUF = JBUF + 8*IRAT
         IPASS = 2
         GO TO 5
      END IF
      WRITE (LUSUPM) P
      JBUF = JBUF + NPQRS*IRAT
#else
      WRITE (LUSUPM) ITYP1,NPF,NQF,NPL,NQL,IOFF,NPQRS,NPQL
      WRITE (LUSUPM) P
#endif
      RETURN
      END
C  /* Deck wrsup2 */
      SUBROUTINE WRSUP2(LUSUPM,P,INDP,ITYP2,NP,NQ,IOFPQ,NRSOUT,
     &                  NPQ,JBUF)
C
C 900319-hjaaj
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION P(NRSOUT), INDP(NRSOUT)
#if defined (VAR_SPLITFILES)
#include <dummy.h>
#include <iratdef.h>
#include <2gbdef.h>
      CHARACTER*80 FNNAME, FNNM2
#include <chrnos.h>
C
      IPASS = 1
 5    CONTINUE
      IF ((((JBUF + NRSOUT*(1 + IRAT) .GT. I2GB) .AND. (IPASS .EQ. 2)))
     &     .OR. (((JBUF + 8*IRAT) .GT. I2GB) .AND. (IPASS .EQ. 1))) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
         INQUIRE(UNIT=LUSUPM,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         CALL GPCLOSE(LUSUPM,'KEEP')
         I = LN - 1
         IF (FNNAME(I:I) .NE. '-') THEN
            FNNM2 = FNNAME(1:LN)//'-0'
            LN = LN + 2
         ELSE
            READ(FNNAME(LN:),'(I1)') INUM
            INUM = INUM + 1
            IF (INUM .GT. 9) THEN
               WRITE (LUPRI,'(/A)') ' DALTON needs to split '//
     &              ' a file more than 11 times.',
     &              ' This is currently not supported'
               CALL QUIT('Too many splittings of a file')
            END IF
            FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
         END IF
         CALL GPOPEN(LUSUPM,FNNM2(1:LN),'UNKNOWN',' ',' ',
     &        IDUMMY,.FALSE.)
         JBUF = 0
      END IF
      IF (IPASS .EQ. 1) THEN
         WRITE (LUSUPM) ITYP2,NP,NQ,NP,NQ,IOFPQ,NRSOUT,NPQ
         JBUF = JBUF + 8*IRAT
         IPASS = 2
         GO TO 5
      END IF
      WRITE (LUSUPM) P,INDP
      JBUF = JBUF + NRSOUT*(1 + IRAT)
#else 
      WRITE (LUSUPM) ITYP2,NP,NQ,NP,NQ,IOFPQ,NRSOUT,NPQ
      WRITE (LUSUPM) P,INDP
#endif
      RETURN
      END
