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
C     FILE : intsort.F
C
#define VAR_AOCANORDER
#ifdef OLD_REV_LOG
!===========================================================================
!961213-hjaaj
!    removed check for LSORT .ge. 8*1024**2 (see MAERKE)
!    : This check is not necessary, but instead we outght to
!    : check that buffer length is greater than a minimum !!
!961112-hjaaj
!    delete SORINT after use
!951025-ov
!    Optional subroutine called from DALTON
!    Renamed subroutine ORDER (conflict with DALTON subroutine) to N_ORDER
!    Renamed common block SORT (conflict with DALTON subroutine) to N_SORT
!951013-ov
!    Replaced update with cpp
!--- version n05 ---
!9412230hjaaj:
!  - SORTAO: initialize IPRINT,IPRFIO,INTSYM w/o DATA because in common
!940727-hjaaj
!  - SORTA: new code for finding canonical order
!  - SORTA: new namelist parameter DELAO, default true
!           if (DELAO) then delete AOTWOINT after sorting.
!--- version n04 ---
!940407-hjaaj
!  - corrected error in WRITE(LUPRI,9000) in ORDER
!  - removed FPS code
!  - moved LSORT parameter to cdk LSORT
!  - modified some formats to print higher integer numbers
!901221-hjaaj
!  - corrected error in writing to LUORD
!  - added VCLR and READI
!901023-hjaaj-l04
!  - use XI4LAB(1:2) instead of XI2LAB(3:4) for LUSORT
!890826-hjaaj-t10.1
!  - Absolute addressing of IADUT in subroutine ORDER
!    giving absolute addressing of each symmetry quadruplet on file LUORD
!    Initialize IADUT to -2 (error code if accessed in TR2CTL)
!    (old relative addressing may be obtained with *DEFINE RELADR)
!890619-hjaaj-t9.1
!  - Check for sufficient RINT allocation in ORDER
!  (and increased LRINT to LRECAO in cdk file)
!  - Transfer LUPRI to ORDER via common /SORT/
!890606-hjaaj
!  removed DARECL on LUORD unless *DEFINE DARECL
!890510-hjaaj
!  corrected code in DECK ORDER for INTSYM .ne. 1
!  -- note that INTSYM .ne. 1 has never been tested yet
!  stop if NBAST .gt. MXCORB or 2*NBUF .gt. LRINT
!890501-hjaaj
!  idea: dynamic LUORD buffer length?
!890104-hjaaj
!  IPRINT,INTSYM,THRQ parameters + more
!  save INTSYM as IDATA(2) on LUORD
!  set KEEP(I)=1 if NBAS(I)=0  
!===========================================================================
#endif
C ======================================================================
      SUBROUTINE SORTAO(SCR,LSORT)
C ======================================================================
C     AUTHOR ORIGINAL VERSION: B. ROOS IN ENSILLRE AUG 1981.
C
C     THIS PROGRAM PERFORMS A SORT OF THE AO INTEGRAL FILE (LUINTA)
C     PRIOR TO THE TRANSFORMATION OF THE AO TWO-ELECTRON INTEGRALS
C     THE SORT HAS TWO PURPOSES: IT PRODUCES A LIST OF INTEGRALS ORDERED
C     AFTER SYMMETRY BLOCKS AND AT THE SAME TIME SKIPS INTEGRALS WHICH A
C     NOT TO BE USED IN THE SUBSEQUENT CI OR CASSCF RUN FOR SYMMETRY
C     REASONS. INTEGRAL LABELS ARE SORTED CANONICALLY AND WRITTEN AS
C     RELATIVE PAIR INDICES FOR EACH SYMMETRY BLOCK.
C     MODIFIED IN APRIL 1984 TO ALLOW CANONICAL ORDERING OF THE
C     INTEGRALS IN EACH SYMMETRY BLOCK. THE ORDERED LIST OF INTEGRALS
C     IS ALSO SQUARED ( PQRS AND RSPQ IF PQ NOT EQUAL TO RS) FOR
C     EFFECTIVE USE IN CASSCF TRANSFORMATION.
C     IBM-3090 VERSION IN NOVEMBER 1986 (B. ROOS)
C
C     ARRAYS IN USE:
C     RINT   INTEGRAL BUFFER FOR AO FILE and lusort file
C     SCR    SORT WORK AREA OF LENGTH LSORT
C     IBATCH BATCH NUMBER FOR EACH CANONICALLY ORDERED SYMMETRY BLOCK
C     LASTAD GIVES LAST ADDRESS FOR EACH CHAIN (ONE CHAIN FOR EACH
C            SYMMETRY BLOCK).
C     IS     SYMMETRY LABEL FOR ALL AO'S  MAX NUMBER OF AO'S IS 200.
C
C     OUTPUT IS ON FILE LUSORT IN CHAINS. EACH CHAIN CORRESPONDS TO A
C     SYMMETRY BLOCK OF INTEGRALS. EACH BATCH IN SCR HAS LENGTH LBATCH
C     (MAX = LBMAX) WHERE LBATCH IS OBTAINED AS LSORT/NBP.
C     THE FOLLOWING NUMBER OF BATCHES ARE USED:
C     NSYM=1     NBATCH=1
C     NSYM=2     NBATCH=4
C     NSYM=4     NBATCH=19
C     NSYM=8     NBATCH=106
C     NBATCH IS HOWEVER REDUCED IF ALL SYMMETRY BLOCKS ARE NOT USED.
C     NBP IS THE OPTIMUM NUMBER OF BATCHES USED IN EACH PASS OVER
C     THE AO INTEGRAL FILE. THE NUMBER OF PASSES IS OPTIMIZED WITH
C     RESPECT TO THE TOTAL NUMBER OF I/O ACCESSES (COUNTING ALSO
C     ONE READ OF THE SORTED FILE IN THE TRANSFORMATION PROGRAM).
C     The number of passes is always one in this version.
C     THE MAXIMUM NUMBER OF SYMMETRY BLOCKS IS 666 (D2H SYMMETRY)
C     AT MOST 106 OF THEM CORRESPOND HOWEVER TO NON-ZERO INTEGRALS.
C
C     EXTERNAL REFERENCES: DAFILE - DIRECT ACCESS FILE I/O ROUTINE
C
C          ********** IBM-3090 RELEASE 86 11 28 **********
C
C                     general version for IBM, Alliant, ...
C                     Release May 89 -- hjaaj
C
C
#include <implicit.h>
#include <dummy.h>
#include <n_sort.h>
#include <inftra.h>
C
      DIMENSION SCR(LSORT)
C
      PARAMETER (LRECL = 8*LRECAO)
C     ... record length on LUINTA and LUORD
      PARAMETER (THRQD  = 1.D-15)
C     *********** THRESHOLD FOR ZERO INTEGRAL
      PARAMETER (D0 = 0.0D0)
C
C
      INTEGER*8 NPQRS, LPQRS, NINTAO, NINTSO
C
      REAL*8    XLABEL
      INTEGER*4 XI4LAB(2)
      INTEGER*2 XI2LAB(4)
      EQUIVALENCE (XLABEL, XI4LAB, XI2LAB)
      COMMON /UNITAR/ IUNTAB(1:99)
C
#include <cbisor.h>
#include <priunit.h>
#include <inftap.h>
#include <dacodes.h>
C
#include <ibtdef.h>
C
      WRITE(LUPRI,1000)
 1000 FORMAT(10X,'AO INTEGRAL SORT PROGRAM ACTIVE - RELEASE 86 12 01',
     &  /10X,'Author: Bjoern Roos',
     &  /10X,'Revised by hjaaj for Sirius release Jan89/Apr94/Sep08')
C
      IF (IRAT .NE. 2 ) THEN
         NWARN = NWARN + 1
         WRITE(LUPRI,'(/,3(/A))') 'WARNING-WARNING-WARNING'
     *  ,' This version of INTSORT only tested on computers with'
     *  ,' REAL*8, INTEGER*4, and INTEGER*2.'
     &  ,' --- This may not work!!!! ---'
CHJ  *  ,' You must change the packing scheme.'
CHJ      CALL QUIT('Error in INTSORT, packing scheme will not work.')
      END IF
C
C     Initialize
C
      IPRINT = ISPRINT
      IPRFIO = ISPRFIO
      INTSYM = ISNTSYM
C
C     Read input
C
      CALL GETTIM(TCPU0,TWALL0)
C
C     MINUS ONE AS A R*8 PACKED NUMBER
      XLABEL = D0
      XI4LAB(1) = -1
      XMIN1  = XLABEL
C
C     SET RECORD LENGTHS AND NAMES ON LOGICAL UNITS
C
C
C...  OPEN MOLECULE AO INTEGRAL FILES
C
      CALL GPOPEN(LUINTA,'AOTWOINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) NSYM, (NBAS(I),I=1,8), LBUF, NIBUF, NBITS
C
C     LBUF  = number of AO integrals per buffer
C             (= 600 in old MOLECULE format)
C
      NBUF2 = IRAT*LBUF + NIBUF*LBUF + 1
      NBUF  = (NBUF2-1)/IRAT + 1
C     ********** NBUF  = AO INTEGRAL BUFFER LENGTH
C                NBUF2 = buffer length, counted in integers
      IF (2*NBUF .GT. LRINT) THEN
         WRITE (LUPRI,'(/A/2(A,I8))')
     *      ' Insufficient length of RINT in common N_SORT',
     *      ' Length is',LRINT,', need',2*NBUF
         CALL QUIT('ERROR: common block N_SORT too short.')
      END IF
C
C
      LBMAX = LRECL/8
C     ********** MAX BATCH AREA IN USE
C                double precision: 8 bytes per floating point word.
C
C     SET START ADDRESS ON LUSORT
C
      IDISK = 0
C
C
C     INPUT SECTION
C
C
C   locate **INTSORT in DALTON input stream
C  
      DO I = 1,8
         KEEP(I) = ISKEEP(I)
      END DO
C
Chj-890104 exclude symmetries with no orbitals
      DO 140 I = 1,NSYM
         IF (NBAS(I) .EQ. 0) KEEP(I) = 1
  140 CONTINUE
C
C     NSYM= NUMBER OF IRREDUCIBLE REPRESENTATIONS OF THE MOLECULAR POINT
C           GROUP IN USE.
C     NBAS(I) = NUMBER OF AO BASIS FUNCTIONS OF SYMMETRY SPECIES I.
C
C     KEEP(I) = 1 IF INTEGRALS CONTAINING ORBITALS OF SYMMETRY SPECIES
C                 I ARE NOT USED
C     KEEP(I) = 0 OTHERWISE.
C
C     THE FOLLOWING THREE CLASSES MAY OCCUR WHERE KEEP(I) = 1 FOR SOME I
C     1. FIRST ORDER CASSCF: IF SYMMETRY SPECIES I DOES NOT CONTAIN ANY
C        ACTIVE ORBITALS.
C     2. SECOND ORDER CASSCF: IF SYMMETRY SPECIES I DOES NOT INCLUDE ANY
C        OCCUPIED (INACTIVE OR ACTIVE) ORBITALS.
C     3. CI RUNS: IF THERE ARE NO ORBITALS IN USE OF SYMMETRY SPECIES I.
C
C
C  changed by Robert Berger on 01.06.2001
C  If we use 'I' here as a counting index, - for reasons I don't really
C  understand - g77 with O3 optimization will stumble later when we try to
C  compute the VECTOR IT(I) = I*(I-1)/2. Thus, I changed 'I' to 'J'.
C
C     WRITE(LUPRI,1150) NSYM,(I,I = 1,NSYM)
C
      WRITE(LUPRI,1150) NSYM,(J,J = 1,NSYM)
C
      WRITE(LUPRI,1151) (NBAS(I),I = 1,NSYM)
      WRITE(LUPRI,1152) (KEEP(I),I = 1,NSYM)
 1150 FORMAT(//' NUMBER OF SYMMETRIES',I4,
     *        /' SYMMETRY SPECIES',20X,8I4)
 1151 FORMAT(  ' NUMBER OF BASIS FUNCTIONS',11X,8I4)
 1152 FORMAT(  ' SYMMETRIES DELETED',18X,8I4)
C
      WRITE(LUPRI,'(/A,1P,D10.2)')
     *   ' Threshold for omitting integrals            :',THRQ2
      IF (ABS(THRQ2-THRQD).GT.1.0D-3*THRQD) WRITE(LUPRI,'(1P,A,D10.2)')
     *   ' Threshold has been reset, default value was :',THRQD
      IF (INTSYM .NE. 1) WRITE(LUPRI,'(/A,I5)')
     *   ' Non-symmetric integrals of symmetry',INTSYM
      IF (IPRINT .GT. 0) WRITE(LUPRI,'(/A,I5)')
     *   ' Print level            ',IPRINT
      IF (IPRFIO .GT. 0) WRITE(LUPRI,'(/A,I5)')
     *   ' Print level for fastio ',IPRFIO
C
      IF (IPRFIO .GE. 20) THEN
         CALL FASTIO('TRACE=ON')
         CALL FASTIO('ERRFIX=FIRM')
      END IF
C
C     TEST ON LIMITS OF NSYM
C
      IF (NSYM.LT.1 .OR. NSYM.GT.8) THEN
         WRITE(LUPRI,1100) NSYM
         CALL QUIT('INTSORT ERROR: invalid number of symmetries')
      END IF
 1100 FORMAT(//' NUMBER OF SYMMETRY SPECIES',I5,' IS NOT WITHIN'
     *      ,' ACCEPTED LIMITS'/' NO SORT WILL BE PERFORMED')
C
C     COMPUTE VECTOR IT(I) = I*(I-1)/2
C
      K = 0
      DO 3 I = 1,MXCORB
        IT(I) = K
        K = K+I
    3 CONTINUE
C
C     COMPUTE SYMMETRY LABELS FOR ALL ORBITALS
C     AND NUMBER OF BASIS FUNCTIONS OF EARLIER SYMMETRIES
C
      K = 0
      IBAS(1) = 0
      NBAST   = 0
      DO 5 ISYM = 1,NSYM
        NBASI = NBAS(ISYM)
        NBAST = NBAST + NBASI
        IF(ISYM.NE.1) IBAS(ISYM) = IBAS(ISYM-1)+NBAS(ISYM-1)
        DO 4 I = 1,NBASI
          K = K+1
          IS(K) = ISYM
    4   CONTINUE
    5 CONTINUE
C
      IF (NBAST .GT. MXCORB) THEN
         WRITE (LUPRI,'(/A/2(A,I8))')
     *      ' Insufficient length of IT(*) and IS(*) in common N_SORT',
     *      ' Length is',MXCORB,', need',NBAST
         CALL QUIT('IT(*) and IS(*) too short in common N_SORT')
      END IF
C
C     COMPUTE BATCH NUMBER FOR EACH SYMMETRY BLOCK
      LSPQRS = 0
      NBATCH = 0
      LINTM  = 0
      DO 8 I = 1,666
         IBATCH(I) = 0
    8 CONTINUE
C
      LPQRS = 0
      nbpq_max = 0
      DO 10 ISP = 1,NSYM
         NBP = NBAS(ISP)
         KEEPP = KEEP(ISP)
      DO 10 ISQ = 1,ISP
         NBQ = NBAS(ISQ)
         NBPQ = NBP*NBQ
         IF(ISP.EQ.ISQ) NBPQ = IT(NBP)+NBP
         nbpq_max = max(nbpq_max,nbpq)
         ISPQ = IEOR(ISP-1,ISQ-1)
         KEEPQ = KEEP(ISQ)
      DO 10 ISR = 1,ISP
         NBR = NBAS(ISR)
         ISPQR = IEOR(ISPQ,ISR-1)
         KEEPR = KEEP(ISR)
         ISSM = ISR
         IF(ISR.EQ.ISP) ISSM = ISQ
      DO 10 ISS = 1,ISSM
         NBS = NBAS(ISS)
         NBRS = NBR*NBS
         IF(ISR.EQ.ISS) NBRS = IT(NBR)+NBR
         IF(ISP.EQ.ISR) THEN
           NPQRS = (NBPQ**2+NBPQ)/2
         ELSE
           NPQRS = NBPQ*NBRS
         END IF
         LSPQRS = LSPQRS+1
C
C     Check if this block is totally symmetric
C
         ISPQRS = IEOR(ISPQR,ISS-1) + 1
      IF (ISPQRS .NE. INTSYM) GO TO 10
         KEEPS = KEEP(ISS)
         KEEPT = KEEPP+KEEPQ+KEEPR+KEEPS
      IF(KEEPT.NE.0) GO TO 10
         LPQRS = MAX(LPQRS,NPQRS)
         NBATCH = NBATCH+1
         IBATCH(LSPQRS) = NBATCH
   10 CONTINUE
C
Chj-890102; stop if no integrals
      IF (NBATCH .EQ. 0) THEN
         WRITE (LUPRI,*) 'Number of batches in INTSORT =',NBATCH
         WRITE (LUPRI,*) 'No integrals kept! Probably input error '//
     &                   'in KEEP.'
         CALL QUIT(
     &   'INTSORT: No integrals kept! Probably input error in KEEP.')
      END IF
C
C     PRESET LASTAD (LAST ADDRESS IN EACH CHAIN) 
C        AND NREC   (NUMBER OF RECORDS IN EACH CHAIN)
C
      DO 9 I=1,106
         LASTAD(I)=-1
         NREC(I)=0
         NSOINT(I)=0
         NSOBAT(I)=0
    9 CONTINUE
C
C     COMPUTE OPTIMUM NUMBER OF PASSES OVER THE AO FILE.
C
      LSORTX=MIN(LSORT,NBATCH*LBMAX)
C     P=NBATCH*(LBUF-1)
C     Q=LSORT
C     NSTEP=2*SQRT(P/Q)
C     IF(NSTEP.EQ.0) NSTEP=1
      NSTEP=1
C     ******* ALWAYS ONE STEP IN THIS VERSION
C     NUMBER OF BLOCKS PER PASS
      NBP=(NBATCH+NSTEP-1)/NSTEP
C     BATCH SIZE
      LBATCH=MIN( LBMAX, LSORT/NBP )
      LINT = (LBATCH-2)/2
      IF (LPQRS .LT. LINT) LINT = LPQRS
C     ... done this way instead of MIN(LINT,LPQRS)
C         because LPQRS int*8, while LINT may be int*4.
      LBATCH=2*LINT+2
      LBATC2=IRAT*LBATCH
      WRITE(LUPRI,1200) NBATCH,LBATCH,LSORTX,NSTEP
 1200 FORMAT(/' NUMBER OF BATCHES',I12/
     *        ' BATCH LENGTH     ',I12/
     *        ' TOTAL SORT AREA  ',I12/
     *        ' NUMBER OF STEPS  ',I12)
C
C     NBP IS THE NUMBER OF BATCHES  LBATCH IS THE LENGTH OF EACH BATCH
C
C     SET RECORD LENGTH FOR LUSORT AND NAME IT
C
      LREC1=((8*LBATCH+1023)/1024)*1024
      IF (LREC1.GT.LRECAO) LREC1=LRECAO
      CALL DAOPEN(LUSORT,'SORINT')
#if defined (VAR_DARECL)
      CALL DARECL(LUSORT,LREC1)
#endif
      IF (IPRFIO .GE. 1) CALL FASTIO('STATUS')
C
C
C      BEGIN LOOP OVER INTEGRAL FILE
C
      NINTAO=0
      NINTSO=0
      NCHAIN=0
C
      NBP2=0
      DO 50 ISTEP=1,NSTEP
      NBP1=NBP2+1
      NBP2=NBP2+NBP
      IF(NBP2.GT.NBATCH) NBP2=NBATCH
      WRITE(LUPRI,1130) ISTEP,NBP1,NBP2
 1130 FORMAT(/' STEP',I3,'  NBP RANGE',2I5)
C
C     INITIALIZE THE BATCHES
C
      DO 15 I=1,NBP
         SCR(I*LBATCH-1)=XMIN1
         SCR(I*LBATCH)  =D0
   15 CONTINUE
C
C     ASYNCHRONOUS READ OF AO INTEGRAL FILE
C
      KBUF1=1
      KBUF2=NBUF + 1
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
#if defined (VAR_SPLITFILES)
      JBUF = 608
      CALL READSI(LUINTA,NBUF2,RINT(KBUF1),JBUF) 
#else
      CALL READI(LUINTA,NBUF2,RINT(KBUF1)) 
#endif
   20 CONTINUE
      IIND1 = IRAT*(KBUF1-1 + LBUF)
      IIND2 = IIND1 + LBUF
C     ..... IIND2 for when NIBUF .eq. 2
      NINT  = IINT(IIND1 + NIBUF*LBUF + 1)
#if defined (VAR_SPLITFILES)
      IF(NINT.GE.0) CALL READSI(LUINTA,NBUF2,RINT(KBUF2),JBUF) 
#else
      IF(NINT.GE.0) CALL READI(LUINTA,NBUF2,RINT(KBUF2)) 
#endif
      IF(NINT.EQ.-1) GO TO 35
C     WRITE(LUPRI,'(Z16)') RINT(NBUF+KBUF1-1)
C
C     LOOP OVER INTEGRALS IN THIS BUFFER
C
      NINTAO = NINTAO + NINT
      DO 30 I=1,NINT
      VALUE=RINT(I+KBUF1-1)
C     XLABEL=RINT(I+IBUF+KBUF1-1)
C     WRITE(LUPRI,'(Z8,2Z16)') I,VALUE,XLABEL
      IF (ABS(VALUE).LT.THRQ2) GO TO 30
C
      IF (NIBUF .EQ. 1) THEN
         ILABEL = IINT(IIND1 + I)
#if defined (VAR_INT64)
         NP =IBTAND(IBTSHR(ILABEL,48),IBT16)
         NQ =IBTAND(IBTSHR(ILABEL,32),IBT16)
         NR =IBTAND(IBTSHR(ILABEL,16),IBT16)
         NS =IBTAND(       ILABEL,    IBT16) 
      ELSE
         CALL QUIT('Programming error: NIBUF .ne. 1 for VAR_INT64')
#else
         NP =IBTAND(IBTSHR(ILABEL,24),IBT08)
         NQ =IBTAND(IBTSHR(ILABEL,16),IBT08)
         NR =IBTAND(IBTSHR(ILABEL, 8),IBT08)
         NS =IBTAND(       ILABEL,    IBT08) 
      ELSE
         ILABEL = IINT(IIND1 + I)
         NP =IBTAND(IBTSHR(ILABEL,16),IBT16)
         NQ =IBTAND(       ILABEL,    IBT16) 
         ILABEL = IINT(IIND2 + I)
         NR =IBTAND(IBTSHR(ILABEL,16),IBT16)
         NS =IBTAND(       ILABEL,    IBT16) 
#endif
      END IF
#if defined VAR_AOCANORDER
C     --- indices are assumed to be in canonical order
C         from integral program
#else
C
C     INTRODUCE CANONICAL ORDERING OF INTEGRAL INDICES
C     ... rewritten 940727-hjaaj
      LP=MAX(NP,NQ)
      LQ=MIN(NP,NQ)
      LR=MAX(NR,NS)
      LS=MIN(NR,NS)
      IF (LP .LT. LR .OR. (LP .EQ. LR .AND. LQ .LT. LS)) THEN
         NP=LR
         NQ=LS
         NR=LP
         NS=LQ
      ELSE
         NP=LP
         NQ=LQ
         NR=LR
         NS=LS
      END IF
#endif
C
      ISP=IS(NP)
      ISQ=IS(NQ)
      ISR=IS(NR)
      ISS=IS(NS)
      KEEPT=KEEP(ISP)+KEEP(ISQ)+KEEP(ISR)+KEEP(ISS)
      IF(KEEPT.NE.0) GO TO 30
      LSPQ=IT(ISP)+ISQ
      LSRS=IT(ISR)+ISS
      LSPQRS=IT(LSPQ)+LSRS
C
C     BATCH NUMBER FOR THIS INTEGRAL
C
      KBATCH=IBATCH(LSPQRS)
C     WRITE(LUPRI,6666) NP,NQ,NR,NS,ISP,ISQ,ISR,ISS,VALUE,KBATCH
C6666 FORMAT(1X,4I3,5X,4I2,5X,F10.6,2X,I5)
      IF(KBATCH.EQ.0) GO TO 900
      IF(KBATCH.LT.NBP1.OR.KBATCH.GT.NBP2) GO TO 30
      KBP=KBATCH-NBP1
      NINTSO=NINTSO+1
C
C     RESET LABEL FOR THIS INTEGRAL AS RELATIVE PAIR INDEX LABEL
C     FOR EACH SYMMETRY BLOCK
C
      NPR=NP-IBAS(ISP)
      NQR=NQ-IBAS(ISQ)
      NRR=NR-IBAS(ISR)
      NSR=NS-IBAS(ISS)
      NPQR=(NPR-1)*NBAS(ISQ)+NQR
      IF(ISP.EQ.ISQ) NPQR=IT(NPR)+NQR
      NRSR=(NRR-1)*NBAS(ISS)+NSR
      IF(ISR.EQ.ISS) NRSR=IT(NRR)+NSR
      XLABEL = D0
      XI4LAB(2) = NPQR
      XI4LAB(1) = NRSR
C
C     ALLOCATE INTEGRALS AND LABEL TO BATCH
C
      NSOINT(KBATCH)=NSOINT(KBATCH)+1
      NSOBAT(KBATCH)=NSOBAT(KBATCH)+1
      LENGTH=NSOBAT(KBATCH)
      IPOS=LBATCH*KBP
      SCR(IPOS+LENGTH)=VALUE
      SCR(IPOS+LINT+LENGTH)=XLABEL
      IF (LENGTH.EQ.LINT) THEN
C
C        THIS BATCH IS NOW FULL AND MUST BE EMPTIED
C
         XLABEL = D0
         XI4LAB(1) = LENGTH
         SCR(IPOS+LBATCH) = XLABEL
         NCHAIN=NCHAIN+1
C        WRITE(LUPRI,1120) NCHAIN
C1120    FORMAT(20X,'NCHAIN',I5)
         IDO=IDISK
         NREC(KBATCH)=NREC(KBATCH)+1
         CALL DAWRITE(LUSORT,SCR(IPOS+1),LBATC2,IDISK)
         XLABEL = D0
         XI4LAB(1) = IDO
         SCR(IPOS+LBATCH-1) = XLABEL
         SCR(IPOS+LBATCH)=D0
         NSOBAT(KBATCH) = 0
      END IF
C
C     THIS COMPLETES THE LOOP OVER THIS BUFFER OF AO INTEGRALS
C
   30 CONTINUE
      KBUF1=NBUF+2-KBUF1
      KBUF2=NBUF+2-KBUF2
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
      IPOS=-LBATCH
      DO 45 I=NBP1,NBP2
      IPOS=IPOS+LBATCH
      LENGTH=NSOBAT(I)
      IF (LENGTH.EQ.0) THEN
         XLABEL = SCR(IPOS+LBATCH-1)
         IDO = XI4LAB(1)
         LASTAD(I)=IDO
      ELSE
         XLABEL = D0
         XI4LAB(1) = LENGTH
         SCR(IPOS+LBATCH) = XLABEL
         NCHAIN=NCHAIN+1
         NREC(I)=NREC(I)+1
C        WRITE(LUPRI,1121) NCHAIN
C1121    FORMAT(25X,'NCHAIN',I5)
         IDO=IDISK
         CALL DAWRITE(LUSORT,SCR(IPOS+1),LBATC2,IDISK)
         LASTAD(I)=IDO
      END IF
   45 CONTINUE
C
C     END OF LOOP OVER INTEGRAL FILE
C
   50 CONTINUE
C
      WRITE(LUPRI,1250) NINTAO,NINTSO,NCHAIN
      WRITE(LUPRI,1260) (NREC(I),I=1,NBATCH)
      WRITE(LUPRI,1261) (NSOINT(I),I=1,NBATCH)
 1250 FORMAT(/' NUMBER OF INTEGRALS READ  ',I15,
     *       /' NUMBER OF INTEGRALS SORTED',I15,
     *       /' NUMBER OF BUFFERS ON SORINT',I14)
 1260 FORMAT(/' NUMBER OF RECORDS IN EACH BATCH'/,(10I8))
 1261 FORMAT(/' NUMBER OF INTEGRALS IN EACH BATCH'/,(10I8))
      IF (DELAO) THEN
         CALL GPCLOSE(LUINTA,'DELETE')
         WRITE (LUPRI,*) 'AOTWOINT deleted.'
      ELSE
         CALL GPCLOSE(LUINTA,'KEEP')
      END IF
C
      IF (IPRFIO .GE. 1) CALL FASTIO('STATUS')
C
C     ORDER INTEGRALS IN EACH BLOCK CANONICALLY. ALSO SQUARE
C
      CALL N_ORDER(SCR(1),LUSORT,LSORT)
      IF (IPRFIO .GE. 1) CALL FASTIO('STATUS')
      CALL DARMOV(LUSORT)
C
      CALL GETTIM(TCPU,TWALL)
      TCPU  = (TCPU -TCPU0 )/60.0D0
      TWALL = (TWALL-TWALL0)/60.0D0
      WRITE(LUPRI,1300) TCPU,TWALL
 1300 FORMAT(/' CPU  TIME USED (MINUTES)',F9.3,
     *       /' WALL TIME USED (MINUTES)',F9.3)
      RETURN
C
C     ERROR MESSAGE
C
  900 WRITE(LUPRI,9000) ISP,ISQ,ISR,ISS,NP,NQ,NR,NS
 9000 FORMAT(//10X,'SYMMETRY LABELS',4I4,' NOT CONSISTENT FOR'
     *        ,' INTEGRAL',4I4,/10X,'PROGRAM STOPS HERE')
      END
C ======================================================================
      SUBROUTINE N_ORDER(X,LUSORT,LWORK)
C ======================================================================
C
C     PURPOSE: STARTING FROM A SYMMETRY SORTED LIST OF TWO-ELECTRON
C              INTEGRALS AN ORDERED LIST IS CREATED WHICH IS QUADRATIC
C              IN THE PAIR INDICES: BOTH PQ AND RS RUN OVER ALL ALLOWED
C              SYMMETRY PAIRS INDEPENDENTLY.
C
C     IBM-3090 RELEASE 86 11 28
C
#include <implicit.h>
      DIMENSION X(LWORK)
#include <n_sort.h>
      PARAMETER (LBUFM=LRECAO, D0 = 0.0D0)
C     *********** MAXIMUM SIZE OF OUTPUT BUFFER CORRESPONDS TO 28K BYTES
      REAL*8    XLABEL
      INTEGER*4 XI4LAB(2)
      INTEGER*2 XI2LAB(4)
      EQUIVALENCE (XLABEL, XI4LAB, XI2LAB)
      COMMON /UNITAR/ IUNTAB(1:99)
#include <priunit.h>
#if defined (VAR_RELADR3)
      PARAMETER ( LIADUT = 176)
#else
      PARAMETER ( LIADUT = 1296)
#endif
      DIMENSION IADUT(2*LIADUT),NRECUT(LIADUT)
      EQUIVALENCE (IADUT(LIADUT+1),NRECUT(1))
C
      IF (IPRINT.GE.5) WRITE(LUPRI,1000)
 1000 FORMAT(
     & ' AO INTEGRAL ORDERING PROGRAM ACTIVE - RELEASE 86 12 01',
     &/' SYMMETRY   BASIS     CHAIN LADX NRECX LPQ NPQ NBUF LBUF',
     &   '   LX   NOINTX NOINT IOUTCH NRECUT IADUT',
     &/'   BLOCK   FUNCTIONS   IN')
      LBATC2=IRAT*LBATCH
      IF (LBATCH .GT. LRINT) THEN
         WRITE (LUPRI,'(/A/2(A,I8))')
     *      ' Insufficient length of RINT in common N_SORT',
     *      ' Length is',LRINT,', need',LBATCH
         CALL QUIT('Insufficient allocation for RINT in N_ORDER')
      END IF
C
C     Initialize IADUT to error indicator (code -2)
C
      DO 10 I = 1,2*LIADUT
         IADUT(I) = -2
   10 CONTINUE
C
C     PREWRITE TWO RECORDS ON OUTPUT UNIT LUORD
C
      LUORD = -1
      CALL DAOPEN(LUORD,'AOORDINT')
#if defined (VAR_DARECL)
      LRECL = 8*LRECAO
      CALL DARECL(LUORD,LRECL)
#endif
      IADR  =0
      CALL DAWRITE(LUORD,IDATA,18,IADR)
      CALL DAWRITE(LUORD,IADUT,2*LIADUT,IADR)
      NBORD =0
C
      INCH  =0
#if defined (VAR_RELADR)
      IOUTCH=0
#endif
C
C     START FOURFOLD LOOP OVER SYMMETRIES
C
      MSPQ=0
      DO 104 NSP=1,NSYM
       KEEPP=KEEP(NSP)
       NBP  =NBAS(NSP)
       DO 103 NSQ=1,NSP
        KEEPQ=KEEP(NSQ)
        NBQ  =NBAS(NSQ)
        NSPQ =IEOR(NSP-1,NSQ-1)
        MSPQ =MSPQ+1
        MSRS =0
        DO 102 NSR=1,NSYM
         KEEPR=KEEP(NSR)
         NSPQR=IEOR(NSPQ,NSR-1)
         NBR  =NBAS(NSR)
         DO 101 NSS=1,NSR
          MSRS   = MSRS+1
          NSPQRS = IEOR(NSPQR,NSS-1) + 1
         IF (NSPQRS .NE. INTSYM) GO TO 101
          KEEPS  = KEEP(NSS)
          KEEPT  = KEEPP+KEEPQ+KEEPR+KEEPS
         IF (KEEPT.NE.0) GO TO 101
          NBS    = NBAS(NSS)
#if defined (VAR_RELADR)
          IOUTCH = IOUTCH+1
#else
          IOUTCH = ( (NSP**2 - NSP) / 2 + NSQ - 1 ) * 36
     &           +   (NSR**2 - NSR) / 2 + NSS
#endif
          ISW    = MSPQ-MSRS
          INCH   = IT(MAX(MSPQ,MSRS))+MIN(MSPQ,MSRS)
          LADX   = LASTAD(IBATCH(INCH))
          NRECX  = NREC(  IBATCH(INCH))
          NOINTX = NSOINT(IBATCH(INCH))
          IF (IPRINT .GE. 5) 
     *    WRITE(LUPRI,*) '% INCH,IBATCH(INCH),LADX',INCH,IBATCH(INCH),
     &                    LADX
          IF(LADX.EQ.-1) GO TO 101
          NOINT = 0
C
C         START SORT OF SYMMETRY BLOCK NSP,NSQ,NSR,NSS
C
          IF(NSP.EQ.NSQ) THEN
           NPQM = (NBP+NBP**2)/2
          ELSE
           NPQM = NBP*NBQ
          ENDIF
          IF(NSR.EQ.NSS) THEN
           NRS  = (NBR+NBR**2)/2
          ELSE
           NRS  = NBR*NBS
          ENDIF
          IF (IPRINT .GE. 5)
     &       WRITE(LUPRI,*) '% ISW,NPQM,NRS',ISW,NPQM,NRS
C
C         NUMBER OF PQ VALUES PER BUFFER: LPQ
C
          LPQ = MIN(NPQM,LBUFM/NRS)
          IF (LPQ.LT.1) THEN
             WRITE(LUPRI,9000) LBUFM,NRS
 9000 FORMAT(' Buffersize too small in N_ORDER in INTSORT'
     &      /' LBUFM=',I15,' but NRS=',I15,
     &      /' Buffersize reset to NRS')
c            CALL QUIT('INTSORT ERROR: buffersize too small in N_ORDER')
             LPQ = 1
          END IF
          NRECUT(IOUTCH) = LPQ
C         LENGTH OF EACH BUFFER: LBUF
          LBUF  = LPQ*NRS
          LBUF2 = IRAT*LBUF
C         NUMBER OF PQ-VALUES IN CORE: NPQ
          NPQ   = MIN(NPQM,(LWORK/LBUF)*LPQ)
          IF (NPQ .LE. 0) CALL ERRWRK('INTSORT.N_ORDER',-LBUF,LWORK)
C         NUMBER OF BUFFERS IN CORE: NBUF
          NBUF  = (NPQ-1)/LPQ+1
          LX    = NBUF*LBUF
C         LX IS THE CORE SIZE ACTUALLY USED
C
C         START SORT IN STEPS OF NPQ PQ VALUES
C
          IADUT(IOUTCH) = IADR
C
          IPQE  = 0
   20     IPQS  = IPQE+1
          IPQE  = MIN(NPQM,IPQE+NPQ)
          NBUFX = (IPQE-IPQS)/LPQ+1
C         NBUFX IS NUMBER OF BUFFERS ACTUALLY USED IN THIS STEP
          LX    = NBUFX*LBUF
          CALL DZERO(X,LX)
          IDA   = LADX
          IF (IPRINT .GE. 5) THEN
             WRITE(LUPRI,*) '% IPQS,IPQE,NBUFX,LX',IPQS,IPQE,NBUFX,LX
             WRITE(LUPRI,*) '% NOINTX,IDA',NOINTX,IDA
          END IF
C
   21     CALL DAREAD(LUSORT,RINT,LBATC2,IDA)
          XLABEL = RINT(LBATCH)
          NINT   = XI4LAB(1)
          XLABEL = RINT(LBATCH-1)
          IDA    = XI4LAB(1)
          IF (IPRINT .GE. 5) WRITE(LUPRI,*) '% NINT,  IDA',NINT,IDA
          IF (NINT.EQ.0) GO TO 30
C
C         LOOP OVER INTEGRALS IN THIS BUFFER
C
          DO 23 I=1,NINT
           XLABEL=RINT(I+LINT)
           IPQ = XI4LAB(2)
           IRS = XI4LAB(1)
           IF (ISW.GE.0) THEN
            IF(IPQ.GE.IPQS.AND.IPQ.LE.IPQE) THEN
             II=NRS*(IPQ-IPQS)+IRS
             X(II)=RINT(I)
             NOINT=NOINT+1
            ENDIF
           ENDIF
           IF (ISW.LE.0) THEN
            IF(IRS.GE.IPQS.AND.IRS.LE.IPQE) THEN
             II=NRS*(IRS-IPQS)+IPQ
             X(II)=RINT(I)
             NOINT=NOINT+1
            ENDIF
           ENDIF
   23     CONTINUE
   30     IF(IDA.NE.-1) GO TO 21
C
C         WRITE OUT THESE INTEGRALS
C
          IST=1
          DO 40 I=1,NBUFX
           CALL DAWRITE(LUORD,X(IST),LBUF2,IADR)
           IST=IST+LBUF
   40     CONTINUE
          NBORD = NBORD + NBUFX
          IF(IPQE.NE.NPQM) GO TO 20
C
          IF (IPRINT.GE.5)
     *    WRITE(LUPRI,1100) NSP,NSQ,NSR,NSS,NBP,NBQ,NBR,NBS,INCH,LADX,
     *                  NRECX,LPQ,NPQ,NBUF,LBUF,LX,NOINTX,NOINT,
     *                  IOUTCH,NRECUT(IOUTCH),IADUT(IOUTCH)
 1100     FORMAT(1X,4I2,4I3,I4,I7,2I5,I4,2I5,3I7,I4,2I7)
  101    CONTINUE
  102   CONTINUE
  103  CONTINUE
  104 CONTINUE
C
C     WRITE AUXILIARY DATA IN FIRST BLOCK ON OUTPUT FILE
C
      IDATA(1)=NSYM
      IDATA(2)=INTSYM
      DO 110 ISYM=1,NSYM
       IDATA(ISYM+2) = NBAS(ISYM)
       IDATA(ISYM+10)= KEEP(ISYM)
  110 CONTINUE
      IADR=0
      CALL DAWRITE(LUORD,IDATA,18,IADR)
      CALL DAWRITE(LUORD,IADUT,2*LIADUT,IADR)
      WRITE (LUPRI,'(/A,I8)') ' Number of buffers on AOORDINT:',NBORD
C
C     CALL PRPQRS
      CALL DACLOS(LUORD)
      RETURN
C
C
      END
