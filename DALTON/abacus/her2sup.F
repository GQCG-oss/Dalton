#define HER2SUP_DEBUG -1
!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
c===========================================================================
C old REV_LOG
c160207-hjaaj: removed unused, obsolete IND(600),BUF(600) from common
c031021-hjaaj: removed OLDSUP
c940818-hjaaj:
c- Removed check for duplicate integrals (reactivate with *DEFINE DUPLCHCK)
c- move symmetry check in SUPSRT in front of canonical index ordering
c- passed NOSYM to SUPSRT, skip symmetry check if NOSYM
c- SUPSRT: use that hermit already has put the 4 indices in canonical
c  order (when SUPSRTs def. of np,nq,nr,ns was reversed).
c- SUP1,SUP,SUP1D,SUPD: reverse order of np,nq,nr,ns in label
cnov 90-hjaaj:
cSUPDRV: LINT = 1706 instead of 540
c900319-hjaaj: corrected error in NPLAST for SUP1D and SUPD
c- rewritten SUPNWR to write up to NP=NST-1 as type 1 if advantageous,
c  used new routine NDXGTA to determine this.
c- added WRSUP2 to correcte error in SUPNWR
c900112-hjaaj: new routines SUPOWR and SUPNWR for writing to LUSUPM
cwhen OLDSUP and .NOT.OLDSUP, resp.
c900111-hjaaj: added NPLAST to CSPINT
c- Format of LUSUPM if .not.oldsup:
crec 1 : '*******','date    ','time    ','PXSUPMAT'
crec 2 : xfac,nosym,nsym,nbas(1:8)
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
c===========================================================================
C  /* Deck formsup */
      SUBROUTINE FORMSUP(WORK,LWORK,NOSYM,HFXFAC,THRESH,IPRINT)
C
C --- PROGRAM FORMSUP (form super matrices for Fock matrices)
C
C Written by Olav Kvalheim, University of Bergen, Norway.
C
C Revised
C        5-Jul-1984 hjaaj (changed all var. to max 6 char.)
C       23-Sep-1986 tuh   adaption to HERMIT
C        9-Jan-1990 tuh+hjaaj: NOSYM parameter;
C          Dec 1991 - work space in call, tuh
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      LOGICAL NOSYM
      REAL*8  WORK(LWORK)
C
      CALL QENTER('FORMSUP')
C
#if HER2SUP_DEBUG > 0
      iprint = max(iprint,HER2SUP_DEBUG) ! DEBUG
      write (lupri,*) 'HER2SUP DEBUG, print level set to',IPRINT
#endif

      IF (IPRINT .GT. 1) THEN
         CALL TITLER('Output from FORMSUP','*',125)
         WRITE (LUPRI,'(A,2(/A,1P,D10.2)/A,L10)')
     &   '    Precalculated AO two-electron integrals'
     &   //' are transformed to P-supermatrix elements.',
     &   '    Factor on exchange                            :',HFXFAC,
     &   '    Threshold for discarding integrals            :',THRESH,
     &   '    Include integrals for general density matrices:',NOSYM
      ELSE IF (HFXFAC .NE. 0.0D0) THEN
         WRITE (LUPRI,'(A/A,1P,D10.2,A)')
     &   '     (Precalculated AO two-electron integrals'
     &   //' are transformed to P-supermatrix elements.',
     &   '      Threshold for discarding integrals :',THRESH,' )'
      ELSE
         WRITE (LUPRI,'(A/A,1P,D10.2,A)')
     &   '     (Precalculated AO two-electron integrals'
     &   //' are resorted to new file.',
     &   '      Threshold for discarding integrals :',THRESH,' )'
      END IF
      LWRK2 = 0
#if defined (VAR_DUPLCHCK)
C Allocate memory for integer*4 IBIT vector in WORK
C (allocation is big enough for all of WORK used for P vector,
C  so no tests for out of bounds needed.)
      LWRK2 = LWRK2 + LWORK / 32
      LWRK2 = LWRK2 / 2 + 1 ! "/2" because integer*4
#endif
      LWRK1 = LWORK - LWRK2
      KWRK1 = 1
      KWRK2 = KWRK1 + LWRK1
      CALL SUPDRV(NOSYM,WORK(KWRK1),LWRK1,WORK(KWRK2),
     &            HFXFAC,THRESH,IPRINT)
      CALL QEXIT('FORMSUP')
      RETURN
      END
C  /* Deck supdrv */
      SUBROUTINE SUPDRV(NOSYM,WORK,LWORK,IBIT_VEC,
     &                  HFXFAC,THRESH,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION WORK(LWORK), IBIT_VEC(*)
      LOGICAL   NOSYM
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ
#include "inftap.h"

      INTEGER, allocatable :: IS(:), LASTAD(:)

!     statement function
      IT( I ) = (I*(I-1))/2

      CALL GETTIM(DST1,WST1)
C
C*** DEFAULTS
C
      XFAC  = HFXFAC
      THRQ  = THRESH
      KSUP  = 0

      CALL SUP_RDINFO(NOSYM,IPRINT)
C
      LINT   = 1706 ! NOTE: LINT *must* be an even number!!!
      ! calculate off-sets for integer*4 ISCR
      IDABUF = 2*LINT ! 2 for integer*4
      LDA2   = IDABUF + NIBUF*LINT
      LDA21  = LDA2+1
      LDA22  = LDA21+1
      LRDA22 = LDA22/2 ! 2 because integer*4
      IF (2*LRDA22 .NE. LDA22)
     &   CALL QUIT('Program error FORMSUP: 2*LRDA22 .ne. LDA22')
C
C***  SYMMETRY LABELS
C
      allocate( IS(NBAST), LASTAD(NBAST) )
      DO ISYM = 1,NSYM
         DO I = IBAS(ISYM)+1,IBAS(ISYM)+NBAS(ISYM)
            IS(I) = ISYM
         END DO
      END DO

C
C***  SIZE OF FIRST BLOCK on LUDASP
C     NST = starting NP for second block
C     NST1 = last NP in first block
C
      NST1  = NBAS(1)
      LAVAIL = LWORK - LRDA22
C     ... IJ = max. memory for integer*4 INDP
      IJ = (NST1*(NST1+1))/2
      IJ = (IJ-1)/2 + 1
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
      NST1 = MIN(NBAS(1),NST1)
      NST  = NST1 + 1
C
C*** SET UP DYNAMICAL STORAGE FOR SORT ROUTINE
C
C allocate space for buffer with LBUF real*8 + NIBUF*LBUF + 2 integer*4 elements
      LENBUF = LBUF*(2+NIBUF)+2
      LENBUF = (LENBUF-1)/2 + 1 ! divide by 2 because integer*4 -> real*8
      KFREE  = 1
      LFREE  = LWORK
      CALL MEMGET2('REAL','BUFFER',KW1,LENBUF,WORK,KFREE,LFREE)
      KW2   = KFREE
      LSORT = LFREE
      CALL GETTIM(DST,WST)
      CALL SUP_SORT(NOSYM,WORK(KW1),WORK(KW1),WORK(KW2),WORK(KW2),LSORT,
     &            IS,LASTAD,IPRINT)
      CALL MEMCHK('after SUP_SORT',WORK,1)
      CALL GETTIM(DFIN,WFIN)
      DTOT  = DFIN-DST
      WTOT  = WFIN-WST
      IF (IPRINT .GT. 1) WRITE(LUPRI,340) DTOT,WTOT
 340  FORMAT(/'    CPU & WALL TIMES IN SUP_SORT :',2(F9.2,' SEC.')/)
#if HER2SUP_DEBUG > 10
      write(lupri,'(/A)') 'LASTAD array:'
      write(lupri,'(5I8,5X,5I8)') LASTAD(1:NBAST)
#endif
C
C*** START TRANSFERRING INTEGRALS TO SUPERMATRIX FORM.
C*** BEGIN WITH FIRST BLOCK CONTAINING ALL INDICES BETWEEN
C*** 1 AND NST-1. ALL INDICES IN FIRST IRREP (i.e. irrep(NP) = 1).
C
C
C*** FIRST FIND NUMBER OF P-SUPER MATRIX ELEMENTS
C
      CALL GETTIM(DST,WST)
      DTOTD = D0
      WTOTD = D0
      IBEF   = 0
      NPLAST = 0
      NST1   = NST-1
      IJ     = IT(NST)
      NPMAX  = IT(IJ+1)
      KW1    = 1
      KW2    = KW1 + LRDA22
      IJ     = (IJ-1)/2 + 1 ! = memory for integer*4 INDP
      NPM    = LRDA22 + NPMAX + IJ
      LSUP   = LWORK-LRDA22-IJ
      ICHAIN = 1
      IF (NPM.LE.LWORK) THEN
         KW3 = KW2+NPMAX
         CALL SUP1(ICHAIN,NPMAX,WORK(KW1),WORK(KW1),WORK(KW2),
     &             WORK(KW3),IBEF,IBIT_VEC(1),IS,LASTAD)
      ELSE
         CALL GETTIM(DST1D,WST1D)
         KW3 = KW2+LSUP
         CALL SUP1D(ICHAIN,NPMAX,WORK(KW1),WORK(KW1),WORK(KW2),
     &             WORK(KW3),IBEF,IBIT_VEC(1),IS,LASTAD)
         CALL GETTIM(DFIN1D,WFIN1D)
         DTOTD = DTOTD + DFIN1D - DST1D
         WTOTD = WTOTD + WFIN1D - WST1D
      END IF
      IBEF = IBEF+NPMAX
C
C*** CONTINUE WITH REST OF FIRST SYMMETRY
C
      NPLAST = NST1
      IFIRST = NBAS(1)-NPLAST
      IF (IFIRST.EQ.0) GO TO 100
      DO 60 I = 1,IFIRST
         ICHAIN = ICHAIN+1
         NPMX   = NPMAX
         NP     = NST1+I
         IJ     = IT(NP+1)
         NPMAX  = IT(IJ+1)
         NPMM   = NPMAX-NPMX
         IJ     = (IJ-1)/2 + 1
         NPM    = LRDA22 + NPMM + IJ
         LSUP   = LWORK - LRDA22 - IJ
         IF (NPM.LE.LWORK) THEN
            KW3 = KW2+NPMM
            CALL SUP1(ICHAIN,NPMM,WORK(KW1),WORK(KW1),WORK(KW2),
     &                WORK(KW3),IBEF,IBIT_VEC(1),IS,LASTAD)
         ELSE
            CALL GETTIM(DST1D,WST1D)
            KW3 = KW2+LSUP
            CALL SUP1D(ICHAIN,NPMM,WORK(KW1),WORK(KW1),WORK(KW2),
     &                WORK(KW3),IBEF,IBIT_VEC(1),IS,LASTAD)
            CALL GETTIM(DFIN1D,WFIN1D)
            DTOTD = DTOTD + DFIN1D - DST1D
            WTOTD = WTOTD + WFIN1D - WST1D
         END IF
         IBEF=IBEF+NPMM
         NPLAST = NP
  60  CONTINUE
      CALL GETTIM(DFIN,WFIN)
      DTOT = DFIN - DST - DTOTD
      WTOT = WFIN - WST - WTOTD
      IF (IPRINT .GT. 1) THEN
         WRITE(LUPRI,350) DTOT,WTOT
         IF (DTOTD .NE. D0) WRITE(LUPRI,351) DTOTD,WTOTD
      END IF
 350  FORMAT('    TIME IN SUP1     :',2(F9.2,' SEC.'))
 351  FORMAT('    TIME IN SUP1D    :',2(F9.2,' SEC.'))
C
C*** NOW TRANSFER REST OF THE INTEGRALS
C    (irrep(NP) > 1)
C
 100  CONTINUE
      IF(NSYM.EQ.1) GO TO 250
      DST = DFIN
      WST = WFIN
      DTOTD = D0
      WTOTD = D0
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
 190        KK = KK+IT(NBI)
         I1 = IT(NPR+1)
         I2 = IT(NPR)
         IJ = KK + I1
         KOMB = I1-I2
         I1 = IT(I1+1)
         I2 = IT(I2+1)
         NPMAX = KK*KOMB+I1-I2
         IJ = (IJ-1)/2 + 1
         NPM = LRDA22+NPMAX+IJ
         LSUP = LWORK - LRDA22 - IJ
         IF (NPM.LE.LWORK) THEN
            KW3 = KW2+NPMAX
            CALL SUP(ICHAIN,NPMAX,WORK(KW1),WORK(KW1),WORK(KW2),
     &               WORK(KW3),IBIT_VEC(1),IBEF,NP,IS,LASTAD,IPRINT)
         ELSE
            CALL GETTIM(DST1D,WST1D)
            KW3 = KW2+LSUP
            CALL SUPD(ICHAIN,NPMAX,WORK(KW1),WORK(KW1),WORK(KW2),
     &                WORK(KW3),IBIT_VEC(1),IBEF,NP,IS,LASTAD)
            CALL GETTIM(DFIN1D,WFIN1D)
            DTOTD = DTOTD + DFIN1D - DST1D
            WTOTD = WTOTD + WFIN1D - WST1D
         END IF
         IBEF = IBEF+NPMAX
         NPLAST = NP
 200  CONTINUE
      CALL GETTIM(DFIN,WFIN)
      DTOT = DFIN-DST - DTOTD
      WTOT = WFIN-WST - WTOTD
      IF (IPRINT .GT. 1) THEN
         WRITE(LUPRI,370) DTOT,WTOT
         IF (DTOTD .NE. D0) WRITE(LUPRI,371) DTOTD,WTOTD
      END IF
 370  FORMAT('    TIME IN SUP      :',2(F9.2,' SEC.'))
 371  FORMAT('    TIME IN SUPD     :',2(F9.2,' SEC.'))
 250  CONTINUE
C
C*** FINALLY signal end of integrals
C
      WRITE(LUSUPM) -2,0,0,0,0,0,0
      CALL GPCLOSE(LUDASP,' ')
      CALL GETTIM(DFIN,WFIN)
      DTOT = DFIN-DST1 - DTOTD
      WTOT = WFIN-WST1 - WTOTD
      IF (IPRINT .GT. 1) WRITE(LUPRI,380) DTOT,WTOT
 380  FORMAT(/'    TOTAL CPU & WALL TIMES      :',2(F9.2,' SEC.'))
      deallocate(IS,LASTAD)
      RETURN
      END
C  /* Deck sup_rdinfo */
      SUBROUTINE SUP_RDINFO(NOSYM,IPRINT)
#include "implicit.h"
#include "priunit.h"
      LOGICAL NOSYM
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ
#include "inftap.h"

      INTEGER NAOS(8)

      CHARACTER*8 LAB123(3), NEWLBL
      DATA LAB123 /'********','        ','        '/
      DATA NEWLBL /'PXSUPMAT'/

!     statement function
      IT( I ) = (I*(I-1))/2

! 1. Read symmetry and basis set size from LUONEL

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

      IF (NBAS(1) .EQ. 0) NOSYM = .TRUE.
      ! very special case, no basis functions in symmetry 1:
      ! activate NOSYM because FORMSUP does not work if no
      ! basis functions in symmetry 1.

      CALL GETDAT(LAB123(2),LAB123(3))
      REWIND LUSUPM
      WRITE (LUSUPM) LAB123, NEWLBL
      WRITE (LUSUPM) XFAC, POTNUC, NOSYM, NSYM, NBAS

! 2. Read AO buffer info etc. from LUNITA

      IF (LUINTA .LT. 0) CALL GPOPEN(LUINTA,'AOTWOINT',
     &   'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)

      REWIND LUINTA
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) MSYM,NAOS,LBUF,NIBUF,NBITS,LENINT4
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

      LAOINT4 = 2*LBUF+NIBUF*LBUF+1
      IF (LAOINT4 .NE. LENINT4) THEN
         WRITE(LUPRI,*)
     &   'Fatal error LAOINT4 .ne. LENINT4',LAOINT4,LENINT4
         CALL QUIT('LAOINT4 .ne. LENINT4')
      END IF
C
      IF (IPRINT .GE. 3) THEN
         WRITE(LUPRI,1050)(I,I = 1,NSYM)
         WRITE(LUPRI,1051)(NBAS(I),I = 1,NSYM)
      END IF
 1050 FORMAT(/'  Symmetry                    ',8I10)
 1051 FORMAT( '  Basis functions per symmetry',8I10/)

C***  If (NOSYM) then reset symmetry index arrays

      IF (NOSYM) THEN
         NBAS(1:8) = 0
         NSYM = 1
         NBAS(1) = NBAST
      END IF
C
C***  NUMBER OF PREVIOUS BASIS FUNCTIONS
C
      IBAS(1) = 0
      IIBAS(1)  = 0
      DO 6 ISYM = 2,NSYM
         NBI = NBAS(ISYM-1)
         IBAS(ISYM) = IBAS(ISYM-1) + NBI
         IIBAS(ISYM)  = IIBAS(ISYM-1) + IT(NBI+1)
   6  CONTINUE

      RETURN
      END
C  /* Deck sup_sort */
      SUBROUTINE SUP_SORT(NOSYM,RINT,INT4,SCR,ISCR,LSORT,
     &   IS,LASTAD,IPRINT)
#include "implicit.h"
#include "priunit.h"
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ
      DIMENSION RINT(*),SCR(*)
      INTEGER*4 INT4(*),ISCR(*), LABEL1, LABEL2, NUMINT
      LOGICAL   NOSYM, OLDDX 
      INTEGER   IS(*), LASTAD(*)
      INTEGER, allocatable :: NREC(:),IBATCH(:),NSOINT(:)
C
#include "inftap.h"
#include "ibtdef.h"

      CALL QENTER('SUP_SORT')
      allocate ( NREC(NBAST),IBATCH(NBAST),NSOINT(NBAST) )
C
C***  OPEN UNIT LUDASP FOR SORTED INTEGRALS
C
      LUDASP = -1
#ifdef VAR_INT64
      LDA22_DX = LDA22 / 2 ! LDA22 is in integer*4
#else
      LDA22_DX = LDA22
#endif
      CALL GPOPEN(LUDASP,'AO2SORTED.DA','UNKNOWN','DIRECT',
     &   ' ',LDA22_DX,OLDDX)
C
C***  BATCH NUMBER FOR EACH AO
C
      NST1 = NST-1
      NBATCH = NBAST-NST1+1
      IF (NST1 .EQ. 0) THEN
         K = 0
      ELSE
         IBATCH(1:NST1) = 1
         K = 1
      END IF
      DO I = NST,NBAST
         K = K+1
         IBATCH(I) = K
      END DO
C
C***  PRESET LASTAD
C
      DO 12 I = 1,NBAST
         LASTAD(I) = -1
         NREC(I) = 0
         NSOINT(I) = 0
  12  CONTINUE
      IINT = 2*LINT ! 2 because always integer*4
C
C***  NUMBER OF BATCHES PER PASS
C
      NBP = LSORT/(LDA22/2)
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
C
      IAOINT  = 2*LBUF
      LAOINT4 = IAOINT+NIBUF*LBUF+1
C
      REWIND LUINTA
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
C
  20  CONTINUE
      CALL READI4(LUINTA,LAOINT4,INT4)
      NUMINT = INT4(LAOINT4)
      IF(NUMINT.EQ. 0) GO TO 20
      IF(NUMINT.EQ.-1) GO TO 35
C
C***  LOOP OVER INTEGRALS IN THIS BUFFER
C
      NINTAO = NINTAO + NUMINT
      DO 30 I = 1,NUMINT
         VALUE  = RINT(I)
      IF(ABS(VALUE).LE.THRQ) GO TO 30
         IF (NIBUF .EQ. 1) THEN
            LABEL1 = INT4(IAOINT+I)
            NP = IAND(ISHFT(LABEL1,-24),IBT08)
         ELSE
            LABEL1 = INT4(IAOINT+2*I-1)
            LABEL2 = INT4(IAOINT+2*I)
            NP = IAND(ISHFT(LABEL1,-16),IBT16)
         END IF
         KBATCH = IBATCH(NP)
      IF (KBATCH.LT.NBP1.OR.KBATCH.GT.NBP2) GO TO 30
C
C        HERMIT: We know NP is largest index, thus we don't need NQ,NR,NS
C                unless we check symmetry below
C
C***  ONLY INTEGRALS WITH AT MOST TWO DIFFERENT SYMMETRY INDICES
C***  ARE TO BE KEPT
C
         IF (NOSYM) GO TO 25
         IF (NIBUF .EQ. 1) THEN
            NQ = IAND(ISHFT(LABEL1,-16),IBT08)
            NR = IAND(ISHFT(LABEL1, -8),IBT08)
            NS = IAND(      LABEL1     ,IBT08)
         ELSE
            NQ = IAND(      LABEL1     ,IBT16)
            NR = IAND(ISHFT(LABEL2,-16),IBT16)
            NS = IAND(      LABEL2     ,IBT16)
         END IF
         ISP = IS(NP)
         IF (IS(NQ).EQ.ISP) GO TO 25
         IF (IS(NR).EQ.ISP) GO TO 25
         IF (IS(NS).EQ.ISP) GO TO 25
      GO TO 30
  25  CONTINUE
C
         KBP = KBATCH-NBP1
         NINTSO = NINTSO+1
C
C***  ALLOCATE INTEGRAL AND LABEL TO BATCH
C
         NSOINT(KBATCH) = NSOINT(KBATCH)+1
         IPOS   = LDA22*KBP
         LENGTH = ISCR(IPOS+LDA22)
         LENGTH = LENGTH+1
         IRPOS  = IPOS/2 ! "/2" because ISCR is always integer*4
         SCR(IRPOS+LENGTH) = VALUE
         IF (NIBUF .EQ. 1) THEN
            ISCR(IPOS+IINT+LENGTH) = LABEL1
         ELSE
            ISCR(IPOS+IINT+2*LENGTH-1) = LABEL1
            ISCR(IPOS+IINT+2*LENGTH  ) = LABEL2
         END IF
         ISCR(IPOS+LDA22) = LENGTH
C
         IF (LENGTH.GE.LINT) THEN
C           ... THIS BATCH IS NOW FULL AND MUST BE EMPTIED
            NCHAIN = NCHAIN+1
            IDO = IDISK
            NREC(KBATCH) = NREC(KBATCH)+1
            CALL WRITDX(LUDASP,IDISK,LDA22_DX,ISCR(IPOS+1))
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
         CALL WRITDX(LUDASP,IDISK,LDA22_DX,ISCR(IPOS+1))
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

      deallocate ( NREC,IBATCH,NSOINT )
      CALL QEXIT('SUP_SORT')
      RETURN
      END
C  /* Deck sup1 */
      SUBROUTINE SUP1(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IBEF,IBIT,IS,LASTAD)
C***  TRANSFORMATION OF FIRST SYMMETRY
#include "implicit.h"
#include "priunit.h"
#include "ibtdef.h"
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ
#include "inftap.h"
      DIMENSION P(*),BBUF(*)
      INTEGER*4 INDP(*), IBUF(*), IBL1, IBL2, IBIT(*)
      INTEGER   IS(*), LASTAD(*)
C

!     statement function
      IT( I ) = (I*(I-1))/2

C
#if defined (VAR_DUPLCHCK)
      JO(L) = IAND(ISHFT(IBIT((L+31)/32),-((L+31)/32*32-L)),1)
#endif
C
#if HER2SUP_DEBUG > 1
      write (lupri,*) 'SUP1: XFAC =',XFAC
#endif
C
C***  PREPARE TEMPORARY STORAGE OF P-SUPERMATRIX
C
      P(1:NPMAX) = 0.0D0
      CFAC  = 1.0D0 ! Coulomb factor; XFAC is exchange factor
      FPQQS = CFAC - XFAC*0.25D0
      FPQQQ = CFAC - XFAC*0.5D0
      FPQRQ = CFAC - XFAC*0.25D0
      FPQPS = CFAC - XFAC*0.25D0
      FPQPQ = CFAC*0.5D0 - XFAC*0.125
      FPPPS = CFAC - XFAC*0.50D0
      FPPPP = CFAC*0.5D0 - XFAC*0.25D0
      XON2 = XFAC*ON2
      XON4 = XFAC*ON4
      XON8 = XFAC*ON8
#if defined (VAR_DUPLCHCK)
C
C***  PREPARE BIT VECTOR FOR CHECK OF UNIQUE INTEGRALS
C
      KBMAX = NPMAX/32+1
      DO 20 I = 1,KBMAX
  20     IBIT(I) = 0
#endif
C
      IADR = LASTAD(ICHAIN)
#if HER2SUP_DEBUG > 1
      write (lupri,*) 'SUP1: ICHAIN, IADR =',ICHAIN, IADR
#endif
      IF (IADR.EQ.-1) GO TO 900
  50  CALL READDX(LUDASP,IADR,LDA22_DX,IBUF)
      IADR = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF (LENGTH.EQ.0) GO TO 500
      DO 400 I = 1,LENGTH
         IF (NIBUF .EQ. 1) THEN
            IBL1  = IBUF(IDABUF+I)
            NP    = IAND(ISHFT(IBL1,-24),IBT08)
            NQ    = IAND(ISHFT(IBL1,-16),IBT08)
            NR    = IAND(ISHFT(IBL1, -8),IBT08)
            NS    = IAND(      IBL1     ,IBT08)
         ELSE
            IBL1  = IBUF(IDABUF+2*I-1)
            IBL2  = IBUF(IDABUF+2*I)
            NP    = IAND(ISHFT(IBL1,-16),IBT08)
            NQ    = IAND(      IBL1     ,IBT08)
            NR    = IAND(ISHFT(IBL2,-16),IBT08)
            NS    = IAND(      IBL2     ,IBT08)
         END IF
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
         IBIT(IK32) = IOR(IBIT(IK32),ISHFT(1,IK322))
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
      P   (IJKL) = P(IJKL) + FPQQS*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
      IJ = IT(NP)+NS
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL) - XON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  CONTINUE
      P   (IJKL) = FPQQQ*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      P(IJKL)    = P(IJKL)+CFAC*BBUF(I)
      IJ = IT(NP)+NR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      P   (IJKL) = P(IJKL)+CFAC*BBUF(I)
      IJ = IT(NP)+NR
      IF (NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NS)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON4*BBUF(I)
 110  IJ = IT(NP)+NS
      IF(NQ.GT.NR) GO TO 115
C
C***  CONTRIBUTION TO PS,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  CONTINUE
      P   (IJKL) = P(IJKL)+FPQRQ*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
      IJ = IT(NP)+NR
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      P   (IJKL) = P(IJKL)+FPQPS*BBUF(I)
      IJ = IT(NP)+NP
C
C***  CONTRIBUTION TO PP,QS   CANONICAL ORDERING ASSUMED
C
      KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  CONTINUE
      P   (IJKL) = P(IJKL)+FPQPQ*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
      IJ = IT(NP)+NP
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      P   (IJKL) = P(IJKL)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
      IJ = IT(NP)+NR
      KL = IT(NP)+NS
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  CONTINUE
      P   (IJKL) = P(IJKL)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
      IJ = IT(NP)+NR
      KL = IJ
      IJKL = IT(IJ)+KL-IBEF
      P   (IJKL) = P(IJKL)-XON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      P   (IJKL) = FPPPS*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  CONTINUE
      P   (IJKL) = FPPPP*BBUF(I)
 400  CONTINUE
 500  IF(IADR.NE.-1) GO TO 50
C
C***  Transfer P to disk
C
      CALL SUPNWR(P,INDP,IBEF,NPMAX,IS)
C
 900  CONTINUE
      RETURN
      END
C  /* Deck sup */
      SUBROUTINE SUP(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IBIT,
     &               IBEF,NP,IS,LASTAD,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "ibtdef.h"
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
#include "inftap.h"
      DIMENSION P(*),BBUF(*)
      INTEGER*4 INDP(*), IBUF(*), IBL1, IBL2, IBIT(*)
      INTEGER   IS(*), LASTAD(*)
C
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ
C

!     statement function
      IT( I ) = (I*(I-1))/2

C
#if defined (VAR_DUPLCHCK)
      JO(L) = IAND(ISHFT(IBIT((L+31)/32),-((L+31)/32*32-L)),1)
#endif
C
C***  PREPARE TEMPORARY STORAGE OF P-SUPERMATRIX
C
      P(1:NPMAX) = 0.0D0
      CFAC  = 1.0D0 ! Coulomb factor; XFAC is exchange factor
      FPQQS = CFAC - XFAC*0.25D0
      FPQQQ = CFAC - XFAC*0.5D0
      FPQRQ = CFAC - XFAC*0.25D0
      FPQPS = CFAC - XFAC*0.25D0
      FPQPQ = CFAC*0.5D0 - XFAC*0.125
      FPPPS = CFAC - XFAC*0.50D0
      FPPPP = CFAC*0.5D0 - XFAC*0.25D0
      XON2 = XFAC*ON2
      XON4 = XFAC*ON4
      XON8 = XFAC*ON8
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
         WRITE (LUPRI,'(/A,I15)') ' NP    ',NP
         WRITE (LUPRI,'( A,I15)') ' NPMAX ',NPMAX
         WRITE (LUPRI,'( A,I15)') ' II    ',II
         WRITE (LUPRI,'( A,I15)') ' ITP   ',ITP
         WRITE (LUPRI,'( A,I15)') ' IMAX  ',IMAX
         WRITE (LUPRI,'( A,I15)') ' JMAX  ',JMAX
         WRITE (LUPRI,'( A,I15)') ' KMAX  ',KMAX
         WRITE (LUPRI,'( A,I15)') ' KBMAX ',KBMAX
      END IF
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
  50  CALL READDX(LUDASP,IADR,LDA22_DX,IBUF)
      IADR   = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF(LENGTH.EQ.0) GO TO 410
      DO 400 I = 1,LENGTH
C
         IF (NIBUF .EQ. 1) THEN
            IBL1  = IBUF(IDABUF+I)
            NQ    = IAND(ISHFT(IBL1,-16),IBT08)
            NR    = IAND(ISHFT(IBL1, -8),IBT08)
            NS    = IAND(      IBL1     ,IBT08)
         ELSE
            IBL1  = IBUF(IDABUF+2*I-1)
            IBL2  = IBUF(IDABUF+2*I)
            NQ    = IAND(      IBL1     ,IBT08)
            NR    = IAND(ISHFT(IBL2,-16),IBT08)
            NS    = IAND(      IBL2     ,IBT08)
         END IF
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
      IBIT(IK31) = IOR(IBIT(IK31),ISHFT(1,IK32))
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
      P   (ID) = P(ID)+FPQQS*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
      IJ = MPP+NSR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  CONTINUE
      P(ID)    = FPQQQ*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      P(ID)    = P(ID)+CFAC*BBUF(I)
      IJ       = MPP+NRR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NRR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      P(ID)    = P(ID)+CFAC*BBUF(I)
      IJ = MPP+NRR
      IF(NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NSR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON4*BBUF(I)
 110  CONTINUE
      IJ = MPP+NSR
      IF(NQ.GT.NR) GO TO 115
C
C***  CONTRIBUTION TO PS,RQ
C
      KL = IT(NRR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  CONTINUE
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  CONTINUE
      P(ID)    = P(ID)+FPQRQ*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
      IJ = MPP+NRR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      P(ID)    = P(ID)+FPQPS*BBUF(I)
      IJ = MPP+NPR
C
C***  CONTRIBUTION TO PP,QS  CANONICAL ORDERING ASSUMED
C
      KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  CONTINUE
      P(ID)    = P(ID)+FPQPQ*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
      IJ = MPP+NPR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      P(ID)    = P(ID)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
      IJ = MPP+NRR
      KL = IT(NPR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  CONTINUE
      P(ID)    = P(ID)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
      IJ = MPP+NRR
      KL = IJ
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)-XON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      P(ID)    = FPPPS*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  CONTINUE
      P(ID)    = FPPPP*BBUF(I)
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
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
 264  P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  ISP = ISQ, ISP>ISR, ISR = ISS
C
 350  IJ = MPP+NQR
      KL = IT(NRR)+NSR+IIBAS(ISR)
      ID = IT(IJ)+KL-IBEF
      P(ID)    = P(ID)+CFAC*BBUF(I)
 400  CONTINUE
C
 410  IF(IADR.NE.-1) GO TO 50
C
C***  Transfer P to disk
C
      CALL SUPNWR(P,INDP,IBEF,NPMAX,IS)
C
 900  CONTINUE
      RETURN
      END
C  /* Deck sup1d */
      SUBROUTINE SUP1D(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IBEF,IBIT,
     &   IS,LASTAD)
C***  TRANSFORMATION OF FIRST SYMMETRY
#include "implicit.h"
#include "priunit.h"
#include "ibtdef.h"
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ
#include "inftap.h"
      DIMENSION P(*),BBUF(*)
      INTEGER*4 INDP(*), IBUF(*), IBL1, IBL2, IBIT(*)
      INTEGER   IS(*), LASTAD(*)
C

!     statement function
      IT( I ) = (I*(I-1))/2

C
#if defined (VAR_DUPLCHCK)
      JO(L) = IAND(ISHFT(IBIT((L+31)/32),-((L+31)/32*32-L)),1)
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
      P(1:NBB) = 0.0D0
      CFAC  = 1.0D0 ! Coulomb factor; XFAC is exchange factor
      FPQQS = CFAC - XFAC*0.25D0
      FPQQQ = CFAC - XFAC*0.5D0
      FPQRQ = CFAC - XFAC*0.25D0
      FPQPS = CFAC - XFAC*0.25D0
      FPQPQ = CFAC*0.5D0 - XFAC*0.125
      FPPPS = CFAC - XFAC*0.50D0
      FPPPP = CFAC*0.5D0 - XFAC*0.25D0
      XON2 = XFAC*ON2
      XON4 = XFAC*ON4
      XON8 = XFAC*ON8
#if defined (VAR_DUPLCHCK)
C
C***  PREPARE BIT VECTOR FOR CHECK OF UNIQUE INTEGRALS
C
      KBMAX = NPMAX/32+1
      DO 20 I = 1,KBMAX
  20  IBIT(I) = 0
#endif
C
      IADR = LASTAD(ICHAIN)
      IF(IADR.EQ.-1) GO TO 900
  50  CALL READDX(LUDASP,IADR,LDA22_DX,IBUF)
      IADR = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF(LENGTH.EQ.0) GO TO 500
      DO 400 I = 1,LENGTH
         IF (NIBUF .EQ. 1) THEN
            IBL1  = IBUF(IDABUF+I)
            NP    = IAND(ISHFT(IBL1,-24),IBT08)
            NQ    = IAND(ISHFT(IBL1,-16),IBT08)
            NR    = IAND(ISHFT(IBL1, -8),IBT08)
            NS    = IAND(      IBL1     ,IBT08)
         ELSE
            IBL1  = IBUF(IDABUF+2*I-1)
            IBL2  = IBUF(IDABUF+2*I)
            NP    = IAND(ISHFT(IBL1,-16),IBT08)
            NQ    = IAND(      IBL1     ,IBT08)
            NR    = IAND(ISHFT(IBL2,-16),IBT08)
            NS    = IAND(      IBL2     ,IBT08)
         END IF
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
      IBIT(IK32) = IOR(IBIT(IK32),ISHFT(1,IK322))
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
      P(ID) = P(ID)+FPQQS*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
  55  IJ = IT(NP)+NS
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      P(ID) = FPQQQ*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 82
      P(ID) = P(ID)+CFAC*BBUF(I)
  82  IJ = IT(NP)+NR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 101
      P(ID) = P(ID)+CFAC*BBUF(I)
 101  IJ = IT(NP)+NR
      IF(NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NS)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 110
      ID = IJKL-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 110
      ID = IJKL-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
 110  IJ = IT(NP)+NS
      IF(NQ.GT.NR) GO TO 115
C
C***  CONTRIBUTION TO PS,RQ
C
      KL = IT(NR)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  KL = IT(NQ)+NR
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 125
      P(ID) = P(ID)+FPQRQ*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
 125  IJ = IT(NP)+NR
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 142
      P(ID) = P(ID)+FPQPS*BBUF(I)
 142  IJ = IT(NP)+NP
C
C***  CONTRIBUTION TO PP,QS   CANONICAL ORDERING ASSUMED
C
      KL = IT(NQ)+NS
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 162
      P(ID) = P(ID)+FPQPQ*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
 162  IJ = IT(NP)+NP
      KL = IT(NQ)+NQ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 182
      P(ID) = P(ID)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
 182  IJ = IT(NP)+NR
      KL = IT(NP)+NS
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 202
      P(ID) = P(ID)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
 202  IJ = IT(NP)+NR
      KL = IJ
      IJKL = IT(IJ)+KL-IBEF
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      ID = IJKL-IRED
      P(ID) = P(ID)-XON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      P(ID) = FPPPS*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  IF(IJKL.LT.NB1.OR.IJKL.GT.NB2) GO TO 400
      P(ID) = FPPPP*BBUF(I)
 400  CONTINUE
 500  IF(IADR.NE.-1) GO TO 50
C
C
C
C***  Transfer P to disk
C
      IOFF = IBEF + IRED
C     Note: A PQ distribution may be divided into two records
      CALL SUPNWR(P,INDP,IOFF,NBB,IS)
C
C
 800  IRED = IRED+NBB
 900  CONTINUE
      RETURN
      END
C  /* Deck supd */
      SUBROUTINE SUPD(ICHAIN,NPMAX,BBUF,IBUF,P,INDP,IBIT,
     &               IBEF,NP,IS,LASTAD)
#include "implicit.h"
#include "priunit.h"
#include "ibtdef.h"
      PARAMETER (ZERO=0.D0, ONE=1.D0, ON2=0.5D0, ON4=0.25D0)
      PARAMETER (ON8=0.125D0, TH4=0.75D0, TH8=0.375D0)
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ
#include "inftap.h"
      DIMENSION P(*),BBUF(*)
      INTEGER*4 INDP(*), IBUF(*), IBL1, IBL2, IBIT(*)
      INTEGER   IS(*), LASTAD(*)
C

!     statement function
      IT( I ) = (I*(I-1))/2

C
#if defined (VAR_DUPLCHCK)
      JO(L) = IAND(ISHFT(IBIT((L+31)/32),-((L+31)/32*32-L)),1)
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
      P(1:NBB) = 0.0D0
      CFAC  = 1.0D0 ! Coulomb factor; XFAC is exchange factor
      FPQQS = CFAC - XFAC*0.25D0
      FPQQQ = CFAC - XFAC*0.5D0
      FPQRQ = CFAC - XFAC*0.25D0
      FPQPS = CFAC - XFAC*0.25D0
      FPQPQ = CFAC*0.5D0 - XFAC*0.125
      FPPPS = CFAC - XFAC*0.50D0
      FPPPP = CFAC*0.5D0 - XFAC*0.25D0
      XON2 = XFAC*ON2
      XON4 = XFAC*ON4
      XON8 = XFAC*ON8
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
  50  CALL READDX(LUDASP,IADR,LDA22_DX,IBUF)
      IADR = IBUF(LDA21)
      LENGTH = IBUF(LDA22)
      IF(LENGTH.EQ.0) GO TO 410
      DO 400 I = 1,LENGTH
C
         IF (NIBUF .EQ. 1) THEN
            IBL1  = IBUF(IDABUF+I)
            NQ    = IAND(ISHFT(IBL1,-16),IBT08)
            NR    = IAND(ISHFT(IBL1, -8),IBT08)
            NS    = IAND(      IBL1     ,IBT08)
         ELSE
            IBL1  = IBUF(IDABUF+2*I-1)
            IBL2  = IBUF(IDABUF+2*I)
            NQ    = IAND(      IBL1     ,IBT08)
            NR    = IAND(ISHFT(IBL2,-16),IBT08)
            NS    = IAND(      IBL2     ,IBT08)
         END IF
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
      IBIT(IK31) = IOR(IBIT(IK31),ISHFT(1,IK32))
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
      P(ID) = P(ID)+FPQQS*BBUF(I)
C
C***  CONTRIBUTION TO PS,QQ
C
  55  IJ = MPP+NSR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  PQ,QQ
C
  60  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = FPQQQ*BBUF(I)
      GO TO 400
C
  80  IF(NS.NE.NR) GO TO 100
C
C***  PQ,RR
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 82
      ID = ID-IRED
      P(ID) = P(ID)+CFAC*BBUF(I)
  82  IJ = MPP+NRR
      IF(NR.GT.NQ) GO TO 85
C
C***  CONTRIBUTION TO PR,QR
C
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
  85  CONTINUE
C
C***  CONTRIBUTION TO PR,RQ
C
      KL = IT(NRR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
 100  IF(NS.EQ.NQ) GO TO 120
C
C***  PQ,RS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 101
      ID = ID-IRED
      P(ID) = P(ID)+CFAC*BBUF(I)
 101  IJ = MPP+NRR
      IF(NQ.GT.NS) GO TO 105
C
C***  CONTRIBUTION TO PR,SQ
C
      KL = IT(NSR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 110
      ID = ID-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 110
C
C***  CONTRIBUTION TO PR,QS
C
 105  KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 110
      ID = ID-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
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
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  CONTRIBUTION TO PS,QR
C
 115  CONTINUE
      KL = IT(NQR)+NRR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  PQ,RQ
C
 120  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 122
      ID = ID-IRED
      P(ID) = P(ID)+FPQRQ*BBUF(I)
C
C***  CONTRIBUTION TO PR,QQ
C
 122  IJ = MPP+NRR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
 140  IF(NS.EQ.NQ) GO TO 160
C
C***  PQ,PS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 142
      ID = ID-IRED
      P(ID) = P(ID)+FPQPS*BBUF(I)
 142  IJ = MPP+NPR
C
C***  CONTRIBUTION TO PP,QS  CANONICAL ORDERING ASSUMED
C
      KL = IT(NQR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  PQ,PQ
C
 160  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 162
      ID = ID-IRED
      P(ID) = P(ID)+FPQPQ*BBUF(I)
C
C***  CONTRIBUTION TO PP,QQ
C
 162  IJ = MPP+NPR
      KL = IT(NQR)+NQR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
 180  IF(NR.EQ.NQ) GO TO 220
      IF(NS.EQ.NR) GO TO 200
C
C***  PP,RS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 182
      ID = ID-IRED
      P(ID) = P(ID)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PS
C
 182  IJ = MPP+NRR
      KL = IT(NPR)+NSR+IIBAS(ISP)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
C
C***  PP,RR
C
 200  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 202
      ID = ID-IRED
      P(ID) = P(ID)+CFAC*BBUF(I)
C
C***  CONTRIBUTION TO PR,PR
C
 202  IJ = MPP+NRR
      KL = IJ
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)-XON8*BBUF(I)
      GO TO 400
C
 220  IF(NS.EQ.NR) GO TO 240
C
C***  PP,PS
C
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = FPPPS*BBUF(I)
      GO TO 400
C
C***  PP,PP
C
 240  IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = FPPPP*BBUF(I)
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
      P(ID) = P(ID)-XON4*BBUF(I)
      GO TO 400
 264  P(ID) = P(ID)-XON2*BBUF(I)
      GO TO 400
C
C***  ISP = ISQ, ISP>ISR, ISR = ISS
C
 350  IJ = MPP+NQR
      KL = IT(NRR)+NSR+IIBAS(ISR)
      ID = IT(IJ)+KL-IBEF
      IF(ID.LT.NB1.OR.ID.GT.NB2) GO TO 400
      ID = ID-IRED
      P(ID) = P(ID)+CFAC*BBUF(I)
 400  CONTINUE
C
 410  IF(IADR.NE.-1) GO TO 50
C
C
C***  Transfer P to disk
C
      IOFF = IBEF + IRED
C     Note: A PQ distribution may be divided into two records
      CALL SUPNWR(P,INDP,IOFF,NBB,IS)
C
C
 800  IRED = IRED+NBB
 900  CONTINUE
      RETURN
      END
C  /* Deck supnwr */
      SUBROUTINE SUPNWR(P,INDP,IBEF,NPMAX,IS)
C
C     Jan 90 Hans Joergen Aa. Jensen
C     Write to LUSUPM
C
C     On input P(1:NPMAX) contains Psupmat(IBEF+1:IBEF+NPMAX)
C
C     WRITE (LUSUPM) ITYP,NPF,NQF,NPL,NQL,IBEF,NPQRS
C     Two format types:
C     ITYP=1: WRITE (LUSUPM) Psupmat(IBEF+1:IBEF+NPQRS)
C     ITYP=2: WRITE (LUSUPM) P(1:NPQRS),INDP(1:NPQRS)
C
C     note: IBEF=IT(NPQF),       NPQF = IT(NPF) + NQF
C           NPQRS=IT(NPQL)-IOFF, NPQL = IT(NPL) + NQL
C
#include "implicit.h"
      REAL*8    P(*)
      INTEGER*4 INDP(*)
      INTEGER   IS(*)
C
      PARAMETER ( ITYP1 = 1, ITYP2 = 2)
C
#include "inftap.h"
      COMMON/CSPINT/NBAS(8),IBAS(8),NBAST,NSYM,IIBAS(8),
     &      LBUF,NIBUF,NBITS, NPLAST,NST,
     B      LINT,LDA2,IDABUF,LDA21,LDA22,LDA22_DX, LSUP
      COMMON/CSPREA/XFAC,THRQ

!     statement function
      IT( I ) = (I*(I-1))/2

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
     &               NQL,IOFF,NPQL)
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
     &                  IOFPQ,NRSOUT,NPQ)
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
      SUBROUTINE WRSUP1(LUSUPM,P,NPQRS,ITYP1,NPF,NQF,NPL,NQL,IOFF,NPQL)
#include "implicit.h"
#include "priunit.h"
      DIMENSION P(NPQRS)
      WRITE (LUSUPM) ITYP1,NPF,NQF,NPL,NQL,IOFF,NPQRS,NPQL
      WRITE (LUSUPM) P
#if HER2SUP_DEBUG > 5
      write(lupri,*) 'type 1; NPF,NQF,NPL,NQL,IOFF,NPQRS,NPQL',
     &                        NPF,NQF,NPL,NQL,IOFF,NPQRS,NPQL
#if HER2SUP_DEBUG > 25
      write(lupri,*) 'P(1:NPQRS) array type 1'
      write(lupri,'(1P,10D12.5)') P
#endif
#endif
      RETURN
      END
C  /* Deck wrsup2 */
      SUBROUTINE WRSUP2(LUSUPM,P,INDP,ITYP2,NP,NQ,IOFPQ,NRSOUT,NPQ)
C
C 900319-hjaaj
C
#include "implicit.h"
#include "priunit.h"
      REAL*8    P(NRSOUT)
      INTEGER*4 INDP(NRSOUT)
      WRITE (LUSUPM) ITYP2,NP,NQ,NP,NQ,IOFPQ,NRSOUT,NPQ
      WRITE (LUSUPM) P,INDP
#if HER2SUP_DEBUG > 5
      write(lupri,*) 'type 2; NP,NQ,NP,NQ,IOFPQ,NRSOUT,NPQ',
     &                        NP,NQ,NP,NQ,IOFPQ,NRSOUT,NPQ
#endif
      RETURN
      END
!  -- end of abacus/her2sup.F --
