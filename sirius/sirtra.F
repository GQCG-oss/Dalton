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
#ifdef OLD_REV_OLG
!===========================================================================
!
!---
!971201-hjaaj: TRACTL: make check for transformation already there
!              also work when AOSRTINT.DA has been deleted !!!
!951219-hjaaj: TRACTL: defined DM1=-1.0D0 (was undefined); never
!   modify input KTRLVL (was ITRLVL, which is abs(KTRLVL))
!950830-hjaaj: LRECL in *4 words for direct access under IRIX
!941115-hjaaj: SIR_INTOPEN: included FLAG(34) in def. of OLDTRA
!940927-hjaaj: SIR_INTOPEN: s/FLAG(8)/DOMP2/ and s/FLAG(9)/DORSP/
!940909-hjaaj: added PARAGON define label (same code as AIX)
!940517-hjaaj: Ide: gem integraler paa LUINTM uden labels (afhaengigt af ITRLVL)
!   Brug ITRLVL til at bestemme LBINTM (ikke NMORBT men f.eks. NMOCCT + ...)
!   Beregn ogsaa LBINTD v.hj.af LBMXSQ ??
!MAERKE 930305: usikkert,indfoer USEDRC el. LVLDRC i parameterliste til TRACTL
!MAERKE 900111: reconsider in SIR_INTOPEN (allow flag(14) also when flag(34?)
!-------- Revision n05 ends here ---------
!940517-hjaaj
!SORTA,SORTB: changed NTRI to NMBASX (was = NMBASX); removed old TAYLSPACE
!SIR_INTOPEN: deleted write DUMMY for CRAY; attempt open old LUINTM file if FLAG(14)
!        true, irrespective of FLAG(34); fixed ICASE print error on Cray
!        by extending IF (OPINTM) THEN to after ICASE print
!TRAAB : removed unused DUMMY parameter statement
!-------- Revision n04 ends here ---------
!931201-hjaaj
!TRACTL: corrected error for TRAAB; LRBUF for 600 but LBINTM could be up to 1024
! new parameter LBMXSQ specifying max buffer length for seq. files, used for
! calculating LBINTM and as buffer length for LUHALF
!930720-hjaaj
!TRACTL: new parameter PROCC false in CALL PRORB
!921203-hjaaj revision n04 (based on revision m02)
!- added define label DEC (same selection as AIX); removed define label FPS
!SORTA : LDAMIN = LDAMAX / 4 (was LDAMAX / 2)
!TRACTL: print KWORK,MWORK when IPRTRA .gt. 5 (was .ge. 20);
!        SAVE LRBUF; LRBUF defined for 600 integrals/buffer in AO file
!910319-hjaaj
!SORTA,SORTB: reduced NBLOCK to max NTRI, thus only calling DZERO with
!  what is used later; at the same time this gives the most even
!  distribution of distributions per chain.
!TRACD,TRAAB: insert symmetry check in write-to-file loop, gprof showed
!  for Neon that a lot of time was used in skipping zero integrals.
!901028-hjaaj
! TRACTL: new MFRO(1:NSYM) = 0 because wrong integral indices are
!    written to LUINTM if MFRO(isym) .ne. 0
!901025-hjaaj
! TRACTL: JTER=1 only if IPRTRA.gt.1, print TIME IN SORTA if iprtra.gt.0
!900801-hjaaj
! Moved writing of pre-MOLTWOEL information on LUINTM from TRAAB to TRACTL;
!   added NSYM,NORB,NBAS and CMO to information.
! TRACTL: check ITRLVL,CMO against LUINTM information, abandon transformation
!   if the needed mo integrals already available.
! SIR_INTOPEN: revised ITRLVL check against LUINTM for levels 2 and 3;
!   rationalized file openings for different computer types;
!   write warning if flag(14) and integrals not found.
!900307-hjaaj
! TRACTL: empty LUINTM in beginning of transformation to give more disk
! space for temporary files during the transformation.
!900305-hjaaj
! TRACTL: only call drcctl if ITRLVL.gt.0 (not ci) and usedrc
! MAERKE: usikkert,indfoer USEDRC el. LVLDRC i parameterliste til TRACTL
!900302-hjaaj
! TRACTL: reset /CBGETD/ to tell any H2AC on disk now obsolete.
!900226-hjaaj
! TRACTL: if (usedrc) then call DRCCTL
!900115-hjaaj
! NXTH2M,NXTH2D: needtp(ityp) .lt. 0 now means at least one distribution
! of that type needed, but not all (in that case needtp(ityp) .gt. 0)
!900111-hjaaj
! SIR_INTOPEN: local logical variable OPINTM, added FLAG(8:9) [RESPONSE and
!   MP2] to cases when LUINTM is opened.
!900108-hjaaj
! NXTH2D: read LBINTD,LVLDRC after label DRCINFO on LUINTD
!891208-hjaaj
! SORTA,TRACD,TRAAB: print threshold for discarding integrals
! TRACTL: if (thrp.le.0) thrp=1.d-12 (was 1.d-9)
!         if (thrq.le.0) thrq=1.d-12 (was 1.d-8)
!891116-hjaaj
! NXTH2M,NXTH2D: more tests that buffers are allocated/have not been
!   corrupted. Release buffer allocation when no more buffers.
!890710-hjaaj
! Write label 'END INTM' at end of LUINTM in TRAAB
! Write total time in tracd,sortb,traab in tractl if iprtra .gt.0
!    and only write TRANSFORMATION level each time if iprtra .gt. 0
!    and corrected 890615 tractl output error (wrt. ITRSAV) in TRACTL
!890220-hjaaj
! traab,tracd -- NOSYNC directive for read direct access to TRI(*)
!===========================================================================
#endif
C  /* Deck trauth */
      SUBROUTINE TRAUTH(IWUNIT)
C
C     891117-hjaaj
C     Author of two-electron integral transformation routines.
C
      WRITE(IWUNIT,100)
  100 FORMAT(/T2,'Two-electron integral transformation:'/
     *   T5,'Peter Taylor,            University of Lund,    Sweden'/,
     *   T5,'Hans Joergen Aa. Jensen, University of Aarhus,  Denmark')
      RETURN
      END
C  /* Deck sir_intopen */
      SUBROUTINE SIR_INTOPEN
C
C  3-Nov-1989 Hans Joergen Aa. Jensen
C
C  Open mo integral units required by this Sirius calculation,
C  as specified in command input.
C
#include <implicit.h>
#include <dummy.h>
C
#include <iratdef.h>
#include <lbmxsq.h>
C
C Used from common blocks:
C   INFINP : FLAG(*),ITRLVL,ITRFIN,ICI0
C   INFORB : NNORBT
C   INFDIM : LBUFMA
C   INFTRA : LSRTAO
C   INFTAP : LU*
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infdim.h>
#include <inftra.h>
#include <inftap.h>
C
      LOGICAL OPINTM, FNDLAB, OLDTRA, FEXIST
C
      IF (NEWTRA) THEN
         CALL SIR_INTOPEN_NEWTRA
         RETURN
      END IF
C
C *** (1) Specify LSRTAO = 0 to tell TRACTL to sort ao integrals
C     on first call; if old LUINTM LSRTAO is reset to -1,
C     signalling TRACTL to check if LUDA exist.
C
C     Find if we should open LUINTM:
C     need LUINTM in MC/CI optimization or for WRGEOM
C     or for MP2 or for RESPONSE
C
      LSRTAO = 0
      OPINTM = .NOT.FLAG(34) .OR. (FLAG(25) .AND. ABS(ITRFIN).LE.10)
     &     .OR. DOMP2 .OR. DORSP
C
C
      IF (OPINTM) THEN
         OLDTRA = FLAG(14) .AND. .NOT.FLAG(34)
C        ... FLAG(14) means .OLD TRANSFORMATION has been set in input
C            or, if also FLAG(34), that transformation is not needed for
C            next module (which will then be RHF).
         NBINTM = (NNORBT-1)/LBMXSQ + 1
         LBINTM = (NNORBT-1)/NBINTM + 1
C
C Always attempt to open LUINTM to see if AOSRTINT.DA information
C is available
C
         ICASE = 1
         CALL GPINQ(FNINTM,'EXIST',FEXIST)
         IF (FEXIST) THEN
            IF (LUINTM .LE. 0) CALL GPOPEN(LUINTM,FNINTM,'OLD',' ',
     &                                     'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND (LUINTM)
            IF (.NOT.FNDLAB('MOLTWOEL',LUINTM)) THEN
               ICASE    = 2
               FLAG(14) = .FALSE.
            ELSE
               ICASE    = 0
               LSRTAO   = -1
C     ... use of old AOSRTINT.DA may be possible
            END IF
            IF (FLAG(14)) THEN
               KTRLVL = ABS(ITRLVL)
               CALL REWSPL(LUINTM)
               READ (LUINTM)
               READ (LUINTM) LBINTM, JTRLVL
               IF (JTRLVL .EQ. 3 .AND. KTRLVL .EQ. 2) KTRLVL =  7
               IF (JTRLVL .EQ. 2 .AND. KTRLVL .EQ. 3) KTRLVL =  2
               IF (JTRLVL .LT. KTRLVL) THEN
                  FLAG(14) = .FALSE.
                  ICASE = 3
               END IF
            END IF
         ELSE
            FLAG(14) = .FALSE.
C           no old luintm opened, then old integrals not available
            IF (LUINTM .LE. 0) CALL GPOPEN(LUINTM,FNINTM,'UNKNOWN',' ',
     &                                     'UNFORMATTED',IDUMMY,.FALSE.)
         END IF
         IF (OLDTRA .AND. ICASE .GT. 0) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A,I2,3(/A))')
     *         ' WARNING SIR_INTOPEN, old integral transformation not'
     *       //' found as expected, ICASE =',ICASE,
     *         ' - ICASE = 1: MO integral file MOTWOINT does not exist',
     *         ' - ICASE = 2: no MO integrals on MOTWOINT',
     *         ' - ICASE = 3: transformation level on MOTWOINT'
     *       //' not sufficient'
            IF (ICASE .EQ. 3) WRITE(LUPRI,'(A,2I5)')
     *         ' ITRLVL requested & ITRLVL on MOTWOINT :',ITRLVL,JTRLVL
         END IF
         LBUFMA = MAX(LBUFMA,NNORBT,LBINTM)
      END IF
C
C     END OF SIR_INTOPEN (Open MO integral files)
C
      RETURN
      END
C  /* Deck sorta */
      SUBROUTINE SORTA(LUDA,RINT,IINT,SCR,ISCR,MEMSX,MEMTX,IT)
C
C....   THIS SUBROUTINE PERFORMS A SORT OF THE AO INTEGRAL FILE PRIOR
C....    TO THE FIRST HALF OF THE TWO-ELECTRON TRANSFORMATION. THE
C....    ALGORITHM INVOLVES MULTIPLE PASSES OVER THE AO INTEGRAL FILE
C....    AND THE NUMBER OF PASSES AND THE BUFFER SIZE ON THE BACK-
C....    CHAINED DIRECT ACCESS SORT FILE IS DETERMINED BY BOTH THE
C....    AVAILABLE CORE DURING THE TRANSFORMATION AND BY THE REQUIRE-
C....    MENT THAT THE TOTAL NUMBER OF I/O REQUESTS SHOULD BE A
C....    MINIMUM. THE ROUTINE ASSUMES THAT THE ENTIRE SORT FILE
C....    CAN BE HELD ON SYSTEM DISK.
C
C....                            PETER R. TAYLOR   MAR *79
C
C....    FORMAL PARAMETERS
C
C....       RINT --  START ADDRESS OF INTEGRAL BUFFER FOR AO FILE
C...        IINT --  (AS RINT - ALLOWS INTEGER ADDRESSING)
C....       SCR  --  START OF SCRATCH WORK AREA
C....       ISCR --  (AS SCR - ALLOWS INTEGER ADDRESSING)
C....       MEMS --  LENGTH OF SCRATCH SORT AREA
C....       MEMT --    **   *    **    TRANSFORMATION AREA
C
#include <implicit.h>
      DIMENSION IINT(*),ISCR(*),SCR(*),RINT(*)
      INTEGER   IT(*)
#include <iratdef.h>
#include <dummy.h>
C
      PARAMETER (IBUFL = 500)
      INTEGER   IBUF(IBUFL), NAOS(8)
C
      LOGICAL FEXIST, OLDDX
C
#include <priunit.h>
#include "inttra.h"
#include <inftra.h>
#include <inftap.h>
#include <ibtdef.h>
C
      CALL QENTER('SORTA ')
C
C....    PARAMETERS USED HERE
C
C....    THRQ: THRESHOLD TO DISCARD INTEGRALS DURING SORT
C              -- FROM /INFTRA/ , DEFINED IN INPUT --
C
C....    LDAMAX: MAXIMUM DA BUFFER LENGTH (-- now set in TRACTL --)
C
C
C....    MEMS is the memory available for the sort in this routine.
C        MEMT is the memory available for the integral transformation
C             to follow after the sort.
C
C        Determine number of chains from MEMT:
C          NBLOCK is max number of simult. distrib. in transf. step
C          NCHAIN is number of chains for this blocking in the
C                 transformation step
C
C
      MEMS   = MEMSX
      MEMT   = MEMTX
      NBLOCK = MEMT/NMBASX
      NCHAIN = 1 + (NMBASX-1)/NBLOCK
      NBLOCK = 1 + (NMBASX-1)/NCHAIN
C
C....    FIND NSTEP, THE NUMBER OF PASSES OVER SEQUENTIAL AO FILE
C
C
C     here LDABUF = DA buffer length for NSTEP = 1
C
      LDABUF = (IRAT*MEMS) / ((1+IRAT)*NCHAIN) - 1
      IF (LDABUF .GE. LDAMAX) THEN
         NSTEP  = 1
      ELSE
         LDAMIN = LDAMAX / 4
         NSTEP  = LDAMIN / LDABUF + 1
      END IF
C
C....    THIS DEFINES THE NUMBER OF PASSES AS NSTEP.
C....    NOW DETERMINE THE REMAINING SORT STATISTICS
C
C....    NBATCH: THE NUMBER OF BATCH BUFFERS USED IN EACH PASS
C                HERE IN SORTA
C
      NBATCH = 1 + (NCHAIN-1)/NSTEP
C
C....    CHECK STATISTICS FOR OVERFLOW OF MAXIMUM INSTALLED VALUES
C
      IF (NCHAIN .GT. MAXCHA) THEN
         WRITE(LUPRI,200) NSTEP,NBLOCK,NCHAIN,NBATCH
         WRITE(LUPRI,153) MAXCHA
         CALL QUIT('SORTA: TOO MANY CHAINS')
      END IF
  153 FORMAT(/' TRACTL.SORTA: TOO MANY CHAINS, MORE THAN',I5)
C
#if defined (VAR_OLDCODE)
      IF (NBATCH .GT. IBUFL) THEN
         WRITE(LUPRI,200) NSTEP,NBLOCK,NCHAIN,NBATCH
         WRITE(LUPRI,151) IBUFL
         CALL QUIT('SORTA: TOO MANY BUFFERS')
      END IF
  151 FORMAT(/' TRACTL.SORTA: MORE THAN',I5,' BUFFERS PER PASS',
     1       /' Increase IBUFL in SORTA (and probably also in SORTB).')
#else
      IF (NBATCH .GT. IBUFL) THEN
         NSTEP  = (NCHAIN-1)/IBUFL + 1
         NBATCH = 1 + (NCHAIN-1)/NSTEP
      END IF
#endif
C
      IF (JTER.EQ.1) WRITE(LUPRI,200) NSTEP,NBLOCK,NCHAIN,NBATCH
C
C....    NUMBER OF PAIR INDICES PER CHAIN (EXCEPT THE LAST)
C
      NPQ  = NBLOCK
      NDIV = NBATCH*NBLOCK
C
C....    BUFFER LENGTH
C
      LDABUF = (MEMS/(NBATCH+NBATCH/IRAT+1)) - 1
      LDABUF = MIN(LDABUF,LDAMAX)
      LDABUF = (LDABUF/IRAT) * IRAT
      IF (JTER.EQ.1) WRITE (LUPRI,201) LDABUF
C
      IDABUF = IRAT*LDABUF
      LBUF   = IDABUF + LDABUF + 2
      NDAREC = (2 * LPQRS - NMBASX)/ LDABUF + 1
C
C...  OPEN AO INTEGRAL FILE AO2INTFILE_LABEL, typically = "AOTWOINT"
C
      CALL GPOPEN(LUINTA,AO2INTFILE_LABEL,'OLD',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) NSYMAO, (NAOS(I),I=1,8),LBUFAO,NIBUFAO,NBITSAO
C
      IF (NIBUFAO .NE. 1) THEN
         CALL QUIT('Only NIBUF.eq.1 on AOTWOINT has been implemented')
      END IF
C
      IBITMASK  = 2**NBITSAO - 1
      IBITSHR_1 = 3*NBITSAO
      IBITSHR_2 = 2*NBITSAO
      IBITSHR_3 =   NBITSAO
C
C...  OPEN DIRECT ACCESS FILE FOR SORTED AO INTEGRALS "AOSRTINT.DA"
C
      CALL GPINQ('AOSRTINT.DA','EXIST',FEXIST)
      IF (FEXIST) THEN
         IF (LUDA .LE. 0) CALL GPOPEN(LUDA,'AOSRTINT.DA','OLD','DIRECT',
     &                                ' ',LBUF,OLDDX)
         CALL GPCLOSE(LUDA,'DELETE')
      END IF
      CALL GPOPEN(LUDA,'AOSRTINT.DA','NEW','DIRECT',' ',LBUF,OLDDX)
C
C....    PREPARE LOOK-UP TABLE OF BUFFER STARTING ADDRESSES
C
      KAP = 0
      DO 25 K = 1,NBATCH
         IBUF(K) = KAP
         KAP = KAP + LBUF
25    CONTINUE
C
C....   AO INTEGRAL FILE SPEC; MOLECULE FORMAT
C
      LRINT  = IRAT*LBUFAO
      LENINT = LRINT + LBUFAO + 1
C
C....    BEGIN THE LOOP OVER PASSES THROUGH AO INTEGRAL FILE
C
      ISUM  = 0
      IDISK = 1
      NFIN  = 0
      DO 100 ISQ = 1,NSTEP
         NST  = NFIN + 1
         NFIN = NFIN + NDIV
         NOB  = NBATCH
C
C....    RESET THE NUMBER OF BUFFERS(IF NECESSARY) FOR THE LAST PASS
C
         IF (NFIN.GT.NMBASX) THEN
            NOB = (NMBASX-NST+NBLOCK)/NBLOCK
            NFIN = NMBASX
         END IF
C
C....    NOB NOW GIVES THE NUMBER OF BUFFERS IN USE ON THIS PASS
C
C....    INITIALIZE BUFFERS
C
         KAP = 0
         DO 30 K = 1,NOB
            KAP = KAP + LBUF
            ISCR(KAP-1) = -1
            ISCR(KAP)   = 0
30       CONTINUE
C
C....    BEGIN LOOP OVER THE INTEGRAL FILE
C
C....    REWIND IT FIRST
C
         CALL REWSPL(LUINTA)
         CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
         JBUF = 608
1        CONTINUE
#if defined (VAR_SPLITFILES)
         CALL READSI(LUINTA,LENINT,IINT,JBUF)
#else
         CALL READI(LUINTA,LENINT,IINT)
#endif
         NINT = IINT(LENINT)
         IF(NINT.EQ.0)GOTO 1
         IF(NINT.EQ.-1)GOTO 2
C
C....    LOOP OVER INTEGRALS
C
         DO 3 M = 1,NINT
C
C....    PICK UP INTEGRAL, LABEL AND UNPACK THE LATTER
C....    NOTE THAT STANDARD ORDER IS NOT CERTAIN IN SYMINT, SO
C....    P.LT.Q ETC  IS TESTED FOR.
C
            LABEL = IINT(M+LRINT)
            VALUE = RINT(M)
C
C....    CHECK INTEGRAL AGAINST THRESHOLD
C
         IF (ABS(VALUE).LT.THRQ) GOTO 3
            ISS=IBTAND(IBTSHR(LABEL,IBITSHR_1),IBITMASK)
            IR =IBTAND(IBTSHR(LABEL,IBITSHR_2),IBITMASK)
            IQ =IBTAND(IBTSHR(LABEL,IBITSHR_3),IBITMASK)
            IP =IBTAND(       LABEL           ,IBITMASK)
            IF (IP.GE.IQ) THEN
               IPQ = IT(IP) + IQ
            ELSE
               IPQ = IT(IQ) + IP
            END IF
            IF (IR.GE.ISS) THEN
               IRS = IT(IR) + ISS
            ELSE
               IRS = IT(ISS) + IR
            END IF
C
C....    CHECK PQ FIRST
C
            IF (IPQ.GE.NST .AND. IPQ.LE.NFIN) THEN
C
C....    PQ IS IN RANGE. DETERMINE APPROPRIATE BUFFER
C
               IPQS   = IPQ - NST + 1
               IPQB   = (IPQS+NPQ-1)/NPQ
C
C....    ALLOCATE INTEGRAL AND LABEL TO BUFFER
C
               IBAD   = IBUF(IPQB)
               IRBAD  = IBAD/IRAT
               LENGTH = ISCR(IBAD+LBUF) + 1
               SCR(IRBAD+LENGTH) = VALUE
C              ISCR(IBAD+IDABUF+LENGTH) = IPQ*(2**16) + IRS
               ISCR(IBAD+IDABUF+LENGTH) =
     *            IBTOR( IBTSHL( IPQ , 16 ) , IRS )
               ISCR(IBAD+LBUF) = LENGTH
               IF (LENGTH .EQ. LDABUF) THEN
C
C....    THIS BUFFER IS NOW FULL AND MUST BE EMPTIED
C
                  ISUM = ISUM + LDABUF
                  IDO  = IDISK
                  CALL WRITDX(LUDA,IDISK,LBUF,ISCR(IBAD+1))
                  ISCR(IBAD+LBUF-1) = IDO
                  ISCR(IBAD+LBUF)   = 0
                  IDISK = IDISK + 1
               END IF
            END IF
C
C....    NOW CHECK RS
C
         IF (IPQ.EQ.IRS) GOTO 3
            IF (IRS.GE.NST .AND. IRS.LE.NFIN) THEN
               IRSS = IRS - NST + 1
               IRSB = (IRSS+NPQ-1)/NPQ
C
C....    ALLOCATE TO APPROPRIATE BUFFER
C
               IBAD   = IBUF(IRSB)
               IRBAD  = IBAD/IRAT
               LENGTH = ISCR(IBAD+LBUF) + 1
               SCR(IRBAD+LENGTH) = VALUE
C              ISCR(IBAD+IDABUF+LENGTH) = IRS*(2**16) + IPQ
               ISCR(IBAD+IDABUF+LENGTH) =
     *            IBTOR( IBTSHL( IRS , 16 ) , IPQ )
               ISCR(IBAD+LBUF) = LENGTH
               IF (LENGTH .EQ. LDABUF) THEN
C
C....    BUFFER IS FULL AND MUST BE EMPTIED
C
                  ISUM = ISUM + LDABUF
                  IDO = IDISK
                  CALL WRITDX(LUDA,IDISK,LBUF,ISCR(IBAD+1))
                  ISCR(IBAD+LBUF-1) = IDO
                  ISCR(IBAD+LBUF)   = 0
                  IDISK = IDISK + 1
               END IF
            END IF
    3    CONTINUE
C
C....    THIS COMPLETES THE LOOP OVER THIS BUFFER OF AO INTEGRALS.
C....    LOOP BACK TO READ ANOTHER
C
         GOTO 1
    2    CONTINUE
C
C....    AT THIS POINT END-OF-FILE HAS BEEN ENCOUNTERED ON THE AO
C...     INTEGRAL FILE, THUS SIGNIFYING THE END OF THIS PASS.
C....    IT IS NECESSARY TO WRITE LAST DA BACK-CHAIN ADDRESSES TO
C....    THE LAST ADDRESS TABLE, EMPTYING THE BUFFERS IF THEY STILL
C....    CONTAIN INFORMATION.
C
         KAP = -LBUF
         IOFF = (ISQ-1)*NBATCH
         DO 6 K = 1,NOB
            KAP = KAP + LBUF
            IDO = IDISK
            IF (ISCR(KAP+LBUF) .NE. 0) THEN
C
C....    THIS BUFFER CONTAINS INFORMATION SO IT IS EMPTIED
C
               ISUM = ISUM + ISCR(KAP+LBUF)
               CALL WRITDX(LUDA,IDISK,LBUF,ISCR(KAP+1))
               LASTAD(IOFF+K) = IDO
               IDISK = IDISK + 1
            ELSE
               IDO = ISCR(KAP+LBUF-1)
               LASTAD(IOFF+K) = IDO
            END IF
    6    CONTINUE
C
C....    END OF THIS PASS. CHECK TIMING
C
100   CONTINUE
      CALL GPCLOSE(LUINTA,'KEEP')
      IF (JTER.EQ.1) THEN
         PNZ = (LPQRS*2-NMBASX) * 100
         PNZ = ISUM / PNZ
         WRITE (LUPRI,202) NDAREC,(IDISK-1),PNZ
         WRITE (LUPRI,203) THRQ
      END IF
      CALL QEXIT('SORTA ')
      RETURN
C
C....    FORMAT STATEMENTS
C
  200 FORMAT(/' FIRST HALF-SORT STATISTICS'
     1       /'  NUMBER OF PASSES OVER AO FILE',T33,I4
     2       /'  NUMBER OF BLOCKS PER CHAIN',T30,I7
     3       /'  TOTAL NUMBER OF CHAINS',T30,I7
     4       /'  NUMBER OF BUFFERS PER PASS',T30,I7/)
  201 FORMAT('  NUMBER OF INTEGRALS PER DIRECT ACCESS BUFFER ',I4)
  202 FORMAT(/'  NUMBER OF DA RECORDS FOR ALL INTEGRALS',I6
     1       /'  NUMBER OF DA RECORDS ACTUALLY USED    ',I6
     2       /'  PERCENT NON-ZERO AO INTEGRALS',F15.2)
  203 FORMAT(/'  THRESHOLD FOR DISCARDING INTEGRALS',1P,D10.2)
      END
C  /* Deck sortb */
      SUBROUTINE SORTB(ITRLVL,LUHALF,LUDA2,RINT,IINT,SCR,ISCR,MEMSX,
     &                 MEMTX,IT)
C
C....    THIS ROUTINE SORTS HALF-TRANSFORMED INTEGRALS FOR THE SECOND
C....    HALF-TRANSFORMATION. AS IN *SORTA*, THE SORT IS OF MULTI-PASS
C....    TYPE, AND THE OPTIMIZING ALGORITHM IS AGAIN DESIGNED TO
C....    MINIMIZE THE NUMBER OF DISK ACCESSES.
C
C....    NOTE THAT THIS VERSION IS DESIGNED TO COPE ONLY WITH MO
C....    INDICES OF THE FORM IJ, OR IA, WHERE I,J ARE OCCUPIED AND
C....    A IS A VIRTUAL INDEX.
C
C....    FORMAL PARAMETERS AND EXTERNALS REQUIRED ARE THE SAME AS
C....    FOR *SORTA* (Q.V.)
C
C....                        PETER R. TAYLOR   LUND   MAR *79
C
#include <implicit.h>
      PARAMETER (IBUFL = 500)
#include <iratdef.h>
      DIMENSION IINT(*),ISCR(*),SCR(*),RINT(*),IBUF(IBUFL)
      INTEGER   IT(*)
      LOGICAL FEXIST
#include <lbmxsq.h>
C
#include <priunit.h>
#include "inttra.h"
#include <inftra.h>
#include <inftap.h>
#include <ibtdef.h>
C
      CALL QENTER('SORTB ')
C
C....    SCRATCH CORE AVAILABILITY
C
      MEMS = MEMSX
      MEMT = MEMTX
      JBUF = 0
C
C....    NUMBER OF MO PAIRS
C        MAXD = number of orbitals for full integral transformation
C               else number of occupied orbitals.
C
C
      IF (ITRLVL.EQ.0) THEN
         NCDA   = NMASHX
         MAXD   = MOCCT
      ELSE IF (ITRLVL .EQ. 10) THEN
         NCDA   = NMORBX
         MAXD   = MORBT
      ELSE
         NCDA   = NMOCCX + MOCCT*MSSHT
         MAXD   = MOCCT
      END IF
C
C     NBLOCK is max number of simult. distrib. in transf. step
C     NCHAIN is number of chains for this blocking in the
C            transformation step
C
      NBLOCK = MEMT/NMBASX
      NCHAIN = 1 + (NCDA-1)/NBLOCK
      NBLOCK = 1 + (NCDA-1)/NCHAIN
C
C....    DETERMINE NUMBER OF PASSES
C
C
C     here LDABUG = DA buffer length for NSTEP = 1
C
      LDABUG = (IRAT*MEMS) / ((1+IRAT)*NCHAIN) - 1
      IF (LDABUG .GE. LDAMAX) THEN
         NSTEP = 1
      ELSE
         LDAMIN = LDAMAX / 4
         NSTEP  = LDAMIN / LDABUG + 1
      END IF
C
C....    NOW DETERMINE THE REMAINING DATA
C
      NBATCH = (NCHAIN+NSTEP-1)/NSTEP
C
C....    CHECK FOR OVERFLOW
C
      IF (NCHAIN .GT. MAXCHA) THEN
         WRITE(LUPRI,200) NSTEP,NBLOCK,NCHAIN,NBATCH
         WRITE(LUPRI,153) MAXCHA
         CALL QUIT('TRACTL.SORTB: TOO MANY CHAINS')
      END IF
  153 FORMAT(/' TRACTL.SORTB: TOO MANY CHAINS; MORE THAN',I5)
C
#if defined (VAR_OLDCODE)
      IF (NBATCH .GT. IBUFL) THEN
         WRITE(LUPRI,200) NSTEP,NBLOCK,NCHAIN,NBATCH
         WRITE(LUPRI,151) IBUFL
         CALL QUIT('TRACTL.SORTB: TOO MANY BUFFERS')
      END IF
  151 FORMAT(/' TRACTL.SORTB: MORE THAN',I5,' BUFFERS PER PASS',
     1       /' Increase IBUFL in SORTB.')
#else
      IF (NBATCH .GT. IBUFL) THEN
         NSTEP  = (NCHAIN-1)/IBUFL + 1
         NBATCH = 1 + (NCHAIN-1)/NSTEP
      END IF
#endif
C
      IF(JTER.EQ.1)WRITE(LUPRI,200)NSTEP,NBLOCK,NCHAIN,NBATCH
C
C....    NO. OF INDICES PER CHAIN
C
      NCD    = NBLOCK
      NDIV   = NBATCH*NBLOCK
C
C....    BUFFER LENGTH
C
      LDABUG = (MEMS/(1+NBATCH+NBATCH/IRAT)) - 1
      LDABUG = MIN(LDABUG,LDAMAX)
      LDABUG = (LDABUG/IRAT) * IRAT
      IF (JTER.EQ.1) WRITE (LUPRI,201) LDABUG
C
      IDABUG = IRAT*LDABUG
      LBUF   = IDABUG + LDABUG + 2
      NDAREC = LPQCD / LDABUG + 1
C
C     delete any existing LUDA2
C
      CALL GPINQ('AOMOTRA.DA','EXIST',FEXIST)
      IF (FEXIST) THEN
         IF (LUDA2 .LE. 0) CALL GPOPEN(LUDA2,'AOMOTRA.DA','OLD',
     &                                'DIRECT','UNFORMATTED',LBUF,OLDDX)
         CALL GPCLOSE(LUDA2,'DELETE')
      END IF
      CALL GPOPEN(LUDA2,'AOMOTRA.DA','NEW','DIRECT','UNFORMATTED',
     &            LBUF,OLDDX)
C
C....    PREPARE LOOK-UP TABLE OF BUFFER ADDRESSES
C
      KAP = 0
      DO 25 K = 1,NBATCH
         IBUF(K) = KAP
         KAP = KAP + LBUF
25    CONTINUE
C
C....    HALF-TRANSFORMED INTEGRAL FILE SPEC
C
      LRINT  = IRAT*LBMXSQ
      LENINT = LRINT + LBMXSQ + 2
C
C....    BEGIN PASSES
C
      ISUM  = 0
      IDISK = 1
      NFIN  = 0
      DO 100 ISQ = 1,NSTEP
         NST  = NFIN + 1
         NFIN = NFIN + NDIV
         NOB  = NBATCH
C
C....    CHECK TO SEE WHETHER IT IS NECESSARY TO RESET THESE VALUES
C....    ON THE LAST PASS
C
         IF (NFIN.GT.NCDA) THEN
            NFIN = NCDA
            NOB  = (NFIN-NST+NBLOCK)/NBLOCK
         END IF
C
C....    NOB GIVES THE ACTUAL NUMBER OF BUFFERS USED IN THIS PASS
C
C....    INITIALIZE BUFFERS
C
         KAP = 0
         DO 30 K = 1,NOB
            KAP = KAP + LBUF
            ISCR(KAP-1) = -1
            ISCR(KAP)   = 0
   30    CONTINUE
C
C....    BEGIN READING THE INTEGRAL FILE
         CALL REWSPL(LUHALF)
    1    CONTINUE
#if defined (VAR_SPLITFILES)
         CALL READSI(LUHALF,LENINT,IINT,JBUF)
#else
         CALL READI(LUHALF,LENINT,IINT)
#endif
         NINT = IINT(LENINT-1)
         IF(NINT.EQ.0)GOTO 1
         IF(NINT.EQ.-1)GOTO 2
C
C....    LOOP OVER THE INTEGRALS IN THIS BUFFER
C
         DO 3 M = 1,NINT
C
C....    PICK UP INTEGRAL AND LABEL, UNPACK IPQ, IC, AND ID
C
            LABEL = IINT(M+LRINT)
            VALUE = RINT(M)
            IPQ=IBTAND(IBTSHR(LABEL,16),IBT16)
            IC =IBTAND(IBTSHR(LABEL, 8),IBT08)
            ID =IBTAND(       LABEL,    IBT08)
C
C....    PACK IC AND ID TO THE INDEX FORM USED IN BUFFER ALLOCATION
C
            IF (IC.GT.MAXD) THEN
               ICDX = MOCCT*IC - NMOCCX + ID
            ELSE
               ICDX = IT(IC) + ID
            END IF
         IF (ICDX.LT.NST.OR.ICDX.GT.NFIN) GOTO 3
C
C....    IN RANGE. PICK UP BUFFER ADDRESS
C
            ICDB   = (ICDX - NST)/NCD + 1
C
C....    ALLOCATE INTEGRAL AND PACK LABEL
C
            IBAD   = IBUF(ICDB)
            IRBAD  = IBAD/IRAT
            LENGTH = ISCR(IBAD+LBUF) + 1
            SCR(IRBAD+LENGTH) = VALUE
C           ISCR(IBAD+IDABUG+LENGTH) = IPQ*2**16 + ICDX
            ISCR(IBAD+IDABUG+LENGTH) =
     *         IBTOR( IBTSHL( IPQ , 16 ) , ICDX )
            ISCR(IBAD+LBUF) = LENGTH
            IF (LENGTH.GE.LDABUG) THEN
C
C....    THIS BUFFER IS FULL AND MUST BE EMPTIED
C
               ISUM = ISUM + LDABUG
               IDO  = IDISK
               CALL WRITDX(LUDA2,IDISK,LBUF,ISCR(IBAD+1))
               ISCR(IBAD+LBUF-1) = IDO
               ISCR(IBAD+LBUF)   = 0
               IDISK= IDISK + 1
            END IF
    3    CONTINUE
C
C....    LOOP OVER THIS BUFFER COMPLETE. READ ANOTHER
C
         GOTO 1
    2    CONTINUE
C
C....    END OF HALF-TRANSFORMED INTEGRAL FILE. EMPTY BUFFERS OF
C....    LAST RECORDS
C
         KAP  = -LBUF
         IOFF = (ISQ-1)*NBATCH
         DO 6 K = 1,NOB
            KAP = KAP + LBUF
            IF (ISCR(KAP+LBUF).NE.0) THEN
C
C....    NEEDS EMPTYING
C
               ISUM = ISUM + ISCR(KAP+LBUF)
               IDO  = IDISK
               CALL WRITDX(LUDA2,IDISK,LBUF,ISCR(KAP+1))
               IDISK = IDISK + 1
            ELSE
               IDO = ISCR(KAP+LBUF-1)
            END IF
            LASTAE(IOFF+K) = IDO
    6    CONTINUE
C
C....    CHECK TIMING AND THEN END THIS PASS
C
  100 CONTINUE
      IF (JTER.EQ.1) THEN
         PNZ = LPQCD * 100
         PNZ = ISUM / PNZ
         WRITE (LUPRI,202) NDAREC,(IDISK-1),PNZ
      END IF
      CALL QEXIT('SORTB ')
      RETURN
C
C....    FORMAT STATEMENTS
C
  200 FORMAT(/,' SECOND HALF-SORT STATISTICS'/
     1  '  NUMBER OF PASSES OVER FILE',T33,I4/
     2  '  NUMBER OF BLOCKS PER CHAIN',T30,I7/
     3  '  TOTAL NUMBER OF CHAINS',T30,I7/
     4  '  NUMBER OF BUFFERS PER PASS',T30,I7/)
  201 FORMAT(' NUMBER OF INTEGRALS PER DIRECT ACCESS BUFFER ',I4)
  202 FORMAT(/'  NUMBER OF DA RECORDS FOR ALL INTEGRALS',I6
     1       /'  NUMBER OF DA RECORDS ACTUALLY USED    ',I6
     2       /'  PERCENT NON-ZERO INTEGRALS',F18.2)
      END
C  /* Deck tracd */
#include <vdir.h>
      SUBROUTINE TRACD(ITRLVL,LUHALF,LUDA,CMOT,BUF,IBUF,OUTB,IOUTB,PQCD,
     &                 PQRD,TRI,ITSMO,ITSAO,ITSX,IT)
C
C Revisions
C   26-Jul-1984 Hans Joergen Aa. Jensen
C   25-Nov-1984 hjaaj (corrected code for ITRLVL.eq.2)
C   30-Jul-1986 hjaaj (use CMOT instead of CMO)
C   10-Feb-1989 hjaaj (ITRLVL=10)
C
C     CASSCF:TRANSFORMATION SECTION.TRANSFORMATION OF LAST TWO INDICES
C
C     THIS SUBROUTINE PERFORMS THE TRANSFORMATION OF THE SECOND INDEX
C     PAIR IN THE TWO-ELECTRON INTEGRALS (PQ/RS). THE RESULT IS A LIST
C     OF HALF TRANSFORMED INTEGRALS (PQ/CD),STORED ON UNIT LUHALF.
C
C     CALLED FROM TRACTL
C
C     CMOT = CMO(transposed)
C
C          ********** RELEASE 79 03 28 **********
C
#include <implicit.h>
      DIMENSION CMOT(*),BUF(*),IBUF(*),OUTB(*),IOUTB(*),
     *          PQCD(MORBT,*),PQRD(*),TRI(*)
      INTEGER   ITSMO(*), ITSAO(*), ITSX(*), IT(*)
#include <iratdef.h>
#include <lbmxsq.h>
C
#include <priunit.h>
#include "inttra.h"
#include <inftra.h>
#include <inftap.h>
C
      LOGICAL SKIPR
C
#include <ibtdef.h>
C
      CALL QENTER('TRACD ')
CHJ-START-840726
      IF (ITRLVL.GE.2 .AND. ITRLVL.LT.4) THEN
         NXISHT = MORBT-MSSHT-MASHT
      ELSE
         NXISHT = 0
      END IF
CHJ-END-840726
      IF (IPRTRA .GE. 100) THEN
         WRITE (LUPRI,*) 'Test output from TRACD'
         WRITE (LUPRI,*) 'ITRLVL =',ITRLVL
         WRITE (LUPRI,*) 'NXISHT =',NXISHT
      END IF
C     ***** AUXILIARY PARAMETERS *****
C     ***** LOUT22: LENGTH OF BUFFER FOR PQCD
C     ***** LDA22 : LENGTH OF DIRECT ACCESS BUFFER
C     ***** LPQCD : NUMBER OF INTEGRALS (PQ/CD)
      IOUT   = 0
      LOUT   = LBMXSQ
      LOUTI  = IRAT*LOUT
      LOUT2  = LOUTI + LOUT
      LOUT21 = LOUT2 + 1
      LOUT22 = LOUT2 + 2
      IDABUF = IRAT*LDABUF
      LDA2   = IDABUF + LDABUF
      LDA21  = LDA2   + 1
      LDA22  = LDA21  + 1
      LPQCD  = 0
      JBUF   = 0
C
C     ***** REWIND FILE FOR HALF TRANSFORMED INTEGRALS LUHALF *****
C
      CALL REWSPL(LUHALF)
C
C     ***** START TRANSFORMATION   *****
C     ***** LOOP OVER ALL PAIRS PQ *****
C
      NPQTOT = (MBAST*MBAST+MBAST)/2
      LPQ    = MIN(NPQ,NPQTOT)
      LTRIPQ = LPQ*NMBASX
      IF (LTRIPQ .GT. LTRI) CALL ERRWRK('TRACD for TRI(pq)',LTRIPQ,LTRI)
C     This will only happen if TRACTL is called with (much) less work
C     space than was available when SORTA was called (first call).
      ICHAIN = 0
      IPQ    = 0
      LPQ    = NPQ
      ISROLD = 0
      SKIPR  = .FALSE.
      DO 100 NP = 1,MBAST
         ISP = ITSAO(NP)
         DO 100 NQ = 1,NP
            ISQ  = ITSAO(NQ)
            ISPQ = MUL(ISP,ISQ)
            IPQ  = IPQ + 1
            IPQL16 = IBTSHL(IPQ,16)
            LPQ  = LPQ + 1
            IF (LPQ .GT. NPQ) THEN
               ICHAIN = ICHAIN + 1
               LPQ = MIN(NPQ,NPQTOT+1-IPQ)
               LTRIPQ = LPQ*NMBASX
               CALL DZERO(TRI,LTRIPQ)
               LPQ = 1
C     ***** FIND LAST ADDRESS OF NEW CHAIN AND START READ *****
               IADR = LASTAD(ICHAIN)
         IF (IADR.LE.0) GO TO 100
               IPQRS0 = NMBASX*(-1 - NPQ*(ICHAIN-1))
  103          CALL READDX(LUDA,IADR,LDA22,IBUF)
                  IADR   = IBUF(LDA21)
                  LENGTH = IBUF(LDA22)
#include <nosync.h>
                  DO 104 I = 1,LENGTH
                     LDAI = IDABUF + I
                     IBL  = IBUF(LDAI)
                     MPQ  = IBTAND(IBTSHR(IBL,16),IBT16)
                     MRS  = IBTAND(       IBL,    IBT16)
C                    ***** FIND ADDRESS AND DISTRIBUTE INTEGRAL *****
                     IPQRS = IPQRS0 + NMBASX*MPQ + MRS
                     TRI(IPQRS) = BUF(I)
  104             CONTINUE
               IF (IADR.NE.-1) GO TO 103
            END IF
C
C     ***** A NEW CHAIN IS NOW TRANSFERRED INTO CORE. START *****
C     ***** TRANSFORMATION OF LAST INDICES FOR THIS PQ      *****
C
C     ***** SET THE ARRAY PQCD TO ZERO *****
C     ***** SET START ADDRESS IN TRI FOR INTEGRALS WITH *****
C     ***** FIRST PAIR INDEX EQUAL TO IPQ.              *****
C     ***** LOOP OVER THIRD INDEX R    *****
C
            CALL DZERO(PQCD,MORBT*MORBT)
            NTPQ = NMBASX*(LPQ-1)
            DO 150 NR = 1,MBAST
               ISR = ITSAO(NR)
               IF (ISR .NE. ISROLD) THEN
                  ISROLD = ISR
                  SKIPR  = .FALSE.
C                 ***** DETERMINE SYMMETRY FOR FOURTH INDEX *****
                  ISS = MUL(ISPQ,ISR)
                  IF (ISS .GT. ISR) THEN
                     IF (ITRLVL .EQ. 0 .OR. ITRLVL .EQ. 10) GO TO 149
                  END IF
C
C                 ***** NUMBER OF BASIS FUNCTIONS AND ORBITALS *****
                  IF (ITRLVL .EQ. 10) THEN
                     NENDD = MORB(ISS)
                  ELSE
                     NENDD = MOCC(ISS)
                  END IF
                  IF (ITRLVL.EQ.0) THEN
                     MASHC = MASH(ISR)
                     MASHD = MASH(ISS)
                     IF (MASHC.EQ.0 .OR. MASHD.EQ.0) GO TO 149
C      V------------------------------------------------
                     NENDC=MOCC(ISR)
                  ELSE
                     NENDC=MORB(ISR)
                     IF (NENDC.EQ.0 .OR. NENDD.EQ.0) GO TO 149
C      V------------------------------------------------
                  END IF
                  MORBR=MORB(ISR)
                  MORBS=MORB(ISS)
                  MBASS=MBAS(ISS)
C                 ***** NUMBER OF ORBITALS IN PRECEDING SYMMETRIES ***
C                 ***** NUMBER OF BASIS FUNCTIONS OF EARLIER SYMMETRIES
                  JORBC=JORB(ISR)
                  JORBD=JORB(ISS)
                  JBASR=JBAS(ISR)
                  JBASS=JBAS(ISS)
C                 ***** SET LOOP LIMITS FOR ORBITAL C *****
                  NCRMAX=NENDC
                  IF (ITRLVL.EQ.0) THEN
                     NCRMIN = MISH(ISR)+1
                     NDRMIN = MISH(ISS)+1
                  ELSE IF (ITRLVL.EQ.10) THEN
                     NCRMIN = 1
                     NDRMIN = 1
                  ELSE IF (ITRLVL.GE.2 .AND. ITRLVL.LT.4) THEN
                     IF (ISS.GT.ISR) THEN
                        NCRMIN = MOCC(ISR)+1
                        NDRMIN = MISH(ISS)+1
C                      (we only need (ai/ integrals
C                       where ITSMO(a) = ITSMO(i) )
                     ELSE
                        NCRMIN = 1
                        NDRMIN = 1
C                      (we only need (ii/ integrals when
C                       ISP=ISQ=ISR=ISS, but we always
C                       need (iu/ and (ui/, so no gain here )
                     END IF
                  ELSE
C                 -- ITRLVL.EQ.1 or 4 or 5:
                     IF (ISS.GT.ISR) THEN
                        NCRMIN = MOCC(ISR)+1
                     ELSE
                        NCRMIN = 1
                     END IF
                     NDRMIN = 1
                  END IF
                  NNDR = 1 + NENDD - NDRMIN
                  IF (NCRMIN.GT.NCRMAX) GO TO 149
                  IF (NNDR .LE. 0) GO TO 149
            ELSE IF (SKIPR) THEN
               GO TO 150
            END IF
C
C     ***** START ADDRESSES FOR MO'S   *****
C     ***** SET THE ARRAY PQRD TO ZERO *****
            IMOR = JCMO(ISR) + NCRMIN + (NR-JBASR-1)*MORBR
            IMOS = JCMO(ISS) + NDRMIN
            CALL DZERO(PQRD,NENDD)
C           ***** LOOP OVER RELATIVE INDEX S IN SYMMETRY ISS *****
            IPQRS1 = NTPQ + NR
            IPQRS2 = NTPQ + IT(NR)
            DO 130 NS = (JBASS+1),(JBASS+MBASS)
C              ***** COMPUTE POSITION OF APPROPRIATE INTEGRAL *****
               IF (NS.GT.NR) THEN
                  IPQRS = IPQRS1 + IT(NS)
               ELSE
                  IPQRS = IPQRS2 + NS
               END IF
               PQRS = TRI(IPQRS)
C              ***** LOOP OVER MO'S OF SYMMETRY ISS *****
               IF (ABS(PQRS).GE.THRP)
     *         CALL DAXPY(NNDR,PQRS,CMOT(IMOS),1,PQRD(NDRMIN),1)
               IMOS = IMOS + MORBS
  130       CONTINUE
C
C     ***** PQRD NOW CONTAINS INTEGRALS (PQ/RD) FOR A GIVEN *****
C     ***** TRIPLE PQR AND ALL SYMMETRY ALLOWED VALUES OF D *****
C     ***** USE THESE TO COMPUTE CONTRIBUTIONS TO (PQ/CD)   *****
C
C
C           ***** LOOP OVER ORBITALS C *****
            DO 140 NCR = NCRMIN,NCRMAX
               CMOR = CMOT(IMOR)
               IMOR = IMOR + 1
            IF (ABS(CMOR) .LE. THRP) GO TO 140
               NC   = NCR + JORBC
C              ***** SET LOOP LIMITS FOR ORBITALS D *****
               IF (ISS.EQ.ISR) THEN
                  NDRMAX = MIN(NCR,NENDD)
               ELSE
                  NDRMAX = NENDD
               END IF
               DO 135 NDR = NDRMIN,NDRMAX
                  PQCD(JORBD+NDR,NC) = PQCD(JORBD+NDR,NC)
     *                               + CMOR*PQRD(NDR)
  135          CONTINUE
  140       CONTINUE
C           END DO NCR
            GO TO 150
C
  149       CONTINUE
C           ***** NOTHING NEEDED FROM THIS SYMMETRY OF R (C)
            SKIPR = .TRUE.
C
  150    CONTINUE
C        END DO NR
C
C        ***** ALL INTEGRALS (PQ/CD) FOR A GIVEN PAIR PQ HAVE *****
C        ***** BEEN CREATED. WRITE THEM ON UNIT LUHALF WITH   *****
C        ***** INDICES IPQ AND ICD (REORDERED PAIR INDEX)     *****
C
         IF (ITRLVL.EQ.0) THEN
            ICMAX=MASHT
         ELSE
            ICMAX=MORBT
         END IF
         DO 170 IC = 1,ICMAX
            NC = ITSX( IC )
            ISC = ITSMO( NC )
            ISD = MUL(ISC,ISPQ)
            IF (ITRLVL.EQ.10) THEN
C           ... full integral transformation
               IDMIN = 1
               IDMAX = IC
            ELSE IF (IC.GT.MOCCT) THEN
               IF (ITRLVL.EQ.2 .AND. ISC.EQ.ISP .AND. ISP.EQ.ISQ) THEN
                  IDMIN = 1
C
C               ( (ai/ distributions only needed for (ai/ai) integrals,
C                ITRLVL.EQ.2, if ITRLVL.eq.3 then (ai/ai) are neglected
C                             if ITRLVL.eq.4 all (ai/bj) are included)
C
               ELSE
                  IDMIN = NXISHT + 1
               END IF
               IDMAX = MIN(IC,MOCCT)
            ELSE
C           -- IC .LE. MOCCT
               IF (IC.LE.NXISHT) THEN
                  IF (ITRLVL.EQ.3) GO TO 170
C                ( if ITRLVL.eq.3 neglect (ii/aa) integrals )
                  IDMIN = IC
C                ( only (ii/ inac-inac distributions needed )
               ELSE
                  IDMIN = 1
               END IF
               IDMAX = MIN(IC,MOCCT)
            END IF
            DO 160 ID = IDMIN,IDMAX
               ND  = ITSX( ID )
            IF (ITSMO( ND ) .NE. ISD) GO TO 160
               IF (NC .EQ. ND) THEN
                  P   = PQCD(NC,NC)
               ELSE
                  P   = PQCD(NC,ND) + PQCD(ND,NC)
               END IF
            IF (ABS(P).LT.THRP) GO TO 160
               LPQCD = LPQCD + 1
               IOUT  = IOUT  + 1
               IF (IOUT.GT.LOUT) THEN
C                 ***** WRITE THIS BUFFER *****
                  IOUT = 1
                  IOUTB(LOUT21) = LOUT
#if defined (VAR_SPLITFILES)
                  CALL WRITSI(LUHALF,LOUT22,IOUTB,JBUF)
#else
                  CALL WRITI(LUHALF,LOUT22,IOUTB)
#endif
                  IF (IPRTRA .GE. 100) THEN
                     WRITE(LUPRI,*)
     &               'Writing buffer for NP =',NP,' NQ = ', NQ
                     WRITE(LUPRI,*) 'IC =',IC,' ID = ', ID
                     WRITE(LUPRI,*) 'First 10 integrals'
                     WRITE(LUPRI,'(5F14.7)') (OUTB(I),I=1,10)
                  END IF
               END IF
               OUTB(IOUT) = P
C              IOUTB(LOUTI+IOUT) = IPQ*2**16 + IC*2**8 + ID
               ICD               = IBTOR(IBTSHL(IC,  8),ID)
               IOUTB(LOUTI+IOUT) = IBTOR(IPQL16,ICD)
C              ... IPQL16 = IBTSHL(IPQ,16)
  160       CONTINUE
  170    CONTINUE
  100 CONTINUE
C        END DO NQ
C     END DO NP
C
C     ***** WRITE LAST BUFFER *****
C
      IOUTB(LOUT21)=IOUT
#if defined (VAR_SPLITFILES)
      CALL WRITSI(LUHALF,LOUT22,IOUTB,JBUF)
#else
      CALL WRITI(LUHALF,LOUT22,IOUTB)
#endif
      IF (IPRTRA .GE. 100) THEN
         WRITE(LUPRI,*) 'Writing last buffer in TRACD'
         WRITE(LUPRI,*) 'First 10 integrals, length =',IOUT
         LPRI = MIN(10,IOUT)
         WRITE(LUPRI,'(5F14.7)') (OUTB(I),I=1,LPRI)
      END IF
C
C     ***** WRITE END BUFFER *****
C
      IOUTB(LOUT21)=-1
#if defined (VAR_SPLITFILES)
      CALL WRITSI(LUHALF,LOUT22,IOUTB,JBUF)
#else
      CALL WRITI(LUHALF,LOUT22,IOUTB)
#endif
      CALL REWSPL(LUHALF)
      IF (JTER.EQ.1) THEN
         WRITE(LUPRI,1000) LUHALF,LPQCD
         WRITE(LUPRI,1100) THRP
      END IF
 1000 FORMAT(/,' NUMBER OF INTEGRALS (PQ/CD) WRITTEN ON UNIT',
     1       I3,' IS',I10)
 1100 FORMAT(/'  THRESHOLD FOR DISCARDING INTEGRALS',1P,D10.2)
      CALL QEXIT('TRACD ')
      RETURN
      END
C  /* Deck traab */
#include <vdir.h>
      SUBROUTINE TRAAB(ITRLVL,LUDA2,CMOT,BUF,IBUF,OUTB,IOUTB,ABCD,PBCD,
     &                 TRI,ITSMO,ITSAO,ITSW,ITSX,IT)
C
C Revised for SIRIUS fall 1983, see below.
C Revisions:
C   25-Nov-1984 hjaaj (corrected code for ITRLVL.eq.2)
C   30-Jul-1986 hjaaj (use CMOT instead of CMO)
C                      CMOT = CMO(transposed)
C   28-Jun-1987 hjaaj (ITRLVL=0, removed AB .ge. CD restriction)
C
C     TRAAB: TRANSFORMATION SECTION-TRANSFORM TWO FIRST INDICES
C
C     THIS SUBROUTINE PERFORMS THE TRANSFORMATION OF THE FIRST
C     INDEX PAIR IN THE TWO-ELECTRON INTEGRALS (PQ/CD).THE
C     RESULT IS A LIST OF TRANSFORMED INTEGRALS (AB/CD),STORED
C     ON UNIT LUINTM IN BUFFERS OF LENGTH 1344
C     (671 INTEGRALS FOLLOWED BY 671 INDICES PLUS A POSITION
C     GIVING THE NUMBER OF INTEGRALS IN THE BUFFER.LAST POSITION
C     IS EMPTY).
C     INDEX PACKING: NA*2**24+NB*2**16+NC*2**8+ND
C
C     NOTE:CANONICAL ORDERING IS NOT ASSURED IN
C          PARTIAL TRANSFORMATION
C
C          ********** RELEASE 79 04 03 **********
C
C Revision 14-Oct-1983 by Hans Jorgen Aa. Jensen, Uppsala
C    This version writes also (AB/CD) where (AB) .lt. (CD),
C    which makes the one-index transformation in NEO easier
C    (could of course also be achieved with a SORTC similar to SORTA)
C    Also start new record when (CD) increments, this means only the
C    (CD) index of the first integral of each record needs to be
C    checked and maybe decoded.
C
C Last revisions  5-Sep-1984 hjaaj /
C                27-Jul-1984 hjaaj /
C                23-Mar-1984 hjaaj /
C                11-Dec-1983 ha/hjaaj
C
#include <implicit.h>
      DIMENSION CMOT(*),BUF(*),IBUF(*),OUTB(*),IOUTB(*),
     *          ABCD(MORBT,*),PBCD(*),TRI(*)
      INTEGER   ITSMO(*), ITSAO(*), ITSW(*), ITSX(*), IT(*)
#include <iratdef.h>
C
#include <priunit.h>
#include <inforb.h>
#include "inttra.h"
#include <inftra.h>
#include <inftap.h>
C
      CHARACTER*8 TABLE(2),LAB123(3)
C
#include <ibtdef.h>
C
      DATA TABLE /'MOLTWOEL','END INTM'/
      DATA LAB123/'********','********','********'/
C
      CALL QENTER('TRAAB ')
      DST=SECOND()
      IF (IPRTRA .GE. 100) THEN
         WRITE (LUPRI,*) 'Test output from TRAAB'
         WRITE (LUPRI,*) 'ITRLVL =',ITRLVL
      END IF
C
      CALL GETDAT(LAB123(2),LAB123(3))
C     ... place date in LAB123(2) and time in LAB123(3)
C
      IF (IPRTRA .GE. 200) THEN
      WRITE (LUPRI,1008)
 1008 FORMAT(/,' *** TRAAB-INFORMATION, two-electron molecular ',
     *   'integrals'/2X,'Int. no.',4X,'C',4X,'D',4X,'A',4X,'B',
     *   10X,'Value')
      END IF
CHJ-START-840726
      IF (ITRLVL.GE.2 .AND. ITRLVL.LT.4) THEN
        NXISHT = MORBT-MASHT-MSSHT
      ELSE
        NXISHT = 0
      END IF
CHJ-END-840726
C ************** LENGTH OF BUFFER FOR ABCD on LUINTM
      LOUT  =LBINTM
      LOUTI =IRAT*LOUT
      LOUT2 =LOUTI+LOUT
      LOUT21=LOUT2+1
      LOUT22=LOUT2+2
C ************* LENGTH OF DIRECT ACCESS BUFFER on LUDA2
      IDABUG=IRAT*LDABUG
      LDA2  =IDABUG+LDABUG
      LDA21 =LDA2+1
      LDA22 =LDA21+1
C ***** initialize counter for NUMBER OF INTEGRALS (AB/CD)
      LABCD =0
      IOUT  =0
      JBUF = 2521*IRAT + 256 + NCMOT
C
C     LUINTM has been positioned before call of TRAAB (900801-hjaaj)
C
      WRITE(LUINTM)LAB123,TABLE(1)
C
      NRINTM = 0
      IF (LPQCD .EQ. 0) THEN
         LABCD = 0
         GO TO 8000
      END IF
C
C     ***** START TRANSFORMATION *****
C     ***** LOOP OVER ALL PAIRS CD *****
      ICD   =0
      LTRICD=NCD*NMBASX
      IF (LTRICD .GT. LTRI) CALL ERRWRK('TRAAB for TRI(cd)',LTRICD,LTRI)
C     should never happen.
      ICHAIN=0
      LCD   =NCD
      IF (ITRLVL.EQ.0) THEN
         ICMAX=MASHT
      ELSE
         ICMAX=MORBT
      END IF
      IF (ITRLVL.LE.5) THEN
         MAXID = MOCCT
      ELSE
         MAXID = MORBT
C        ... full integral transformation
      END IF
      DO 100 IC=1,ICMAX
         NC    = ITSX(IC)
         ISC   = ITSMO(NC)
         ITISC = IT(ISC)
         NCR=NC-JORB(ISC)
CHJ      IDMAX=IC
CHJ      IF (IC.GT.MOCCT) IDMAX=MOCCT
         IF (IC.LE.NXISHT) THEN
            IDMIN = IC
            IDMAX = IC
         ELSE
            IDMIN = 1
            IDMAX = MIN(MAXID,IC)
         END IF
         DO 100 ID=1,IDMAX
            ICD=ICD+1
            LCD=LCD+1
            IF (LCD.LE.NCD) GO TO 110
C           ***** ONE CHAIN COMPLETED. READ AND DISTRIBUTE A NEW CHAIN
            LCD=1
            ICHAIN=ICHAIN+1
C           ***** FIND LAST ADDRESS OF NEW CHAIN AND START READ *****
            IADR=LASTAE(ICHAIN)
            IF(IADR.LE.0) GO TO 100
            CALL DZERO(TRI,LTRICD)
            MCDOFF = 1 + NCD*ICHAIN - NCD
  103       CALL READDX(LUDA2,IADR,LDA22,IBUF)
               IADR=IBUF(LDA21)
               LENGTH=IBUF(LDA22)
            IF(LENGTH.EQ.0) GO TO 106
#include <nosync.h>
               DO 104 I=1,LENGTH
                  IBL=IBUF(IDABUG+I)
                  MPQ=IBTAND(IBTSHR(IBL,16),IBT16)
                  MCD=IBTAND(       IBL,    IBT16)
C                 ***** FIND ADDRESS AND DISTRIBUTE INTEGRAL *****
                  IPQCD=NMBASX*(MCD-MCDOFF)+MPQ
                  TRI(IPQCD)=BUF(I)
  104          CONTINUE
  106       IF(IADR.NE.-1) GO TO 103
CHJ 1 ! IADR = 1 means non-zero integrals in this chain
            IADR = 1
C           ***** A NEW CHAIN IS NOW TRANSFERRED INTO CORE. START *****
C           ***** TRANSFORMATION OF FIRST INDICES FOR THIS CD     *****
  110       CONTINUE
CHJ 2
         IF (IADR.LE.0) GO TO 100
CHJ-S-841018
C This patch (instead of DO 100 ID = IDMIN,IDMAX) was necessary
C to get right CHAIN in core (keep ICD and LCD in agreement
C with the counting in SORTB).
         IF (ID.LT.IDMIN) GO TO 100
CHJ-E-841018
            ND =ITSX(ID)
            ISD=ITSMO(ND)
CHJ-S-841125
            IF (IC.GT.MOCCT .AND. ID.LE.NXISHT) THEN
               IF (ISC.NE.ISD) GO TO 100
C              ( we only need (ai/ distributions for (ai/ai),
C                sym(a)=sym(i) )
            END IF
CHJ-E-841125
            NDR = ND - JORB(ISD)
            ISCD=MUL(ISC,ISD)
            IF (ISD.LE.ISC) THEN
               NSCD=ITISC+ISD
            ELSE
               NSCD=IT(ISD)+ISC
            END IF
            IF (NC.GE.ND) THEN
C              INDCD = NC*2**24 + ND*2**16
               INDCD = NC*2**8  + ND
            ELSE
C              INDCD = ND*2**24 + NC*2**16
               INDCD = ND*2**8  + NC
            END IF
            INDCD = IBTSHL( INDCD , 16 )
C           ***** SET THE ARRAY ABCD TO ZERO *****
            CALL DZERO(ABCD,MORBT*MORBT)
C           ***** SET START ADDRESS IN TRI FOR INTEGRALS WITH *****
C           ***** FIRST PAIR INDEX EQUAL TO ICD.              *****
            NTCD=NMBASX*(LCD-1)
C           ***** LOOP OVER FIRST INDEX P *****
            DO 150 NP=1,MBAST
               ISPA=ITSAO(NP)
C              ***** DETERMINE SYMMETRY FOR SECOND INDEX *****
               ISQB=MUL(ISPA,ISCD)
C              ***** NUMBER OF ORBITALS AND BASIS FUNCTIONS *****
               IF(ITRLVL.EQ.0.AND.MASH(ISQB).EQ.0)GO TO 150
               MORBA=MORB(ISPA)
               IF(MORBA.EQ.0) GO TO 150
               MORBB=MORB(ISQB)
               IF(MORBB.EQ.0) GO TO 150
               MOCCA=MOCC(ISPA)
               MOCCB=MOCC(ISQB)
               MORBP=MORB(ISPA)
               MORBQ=MORB(ISQB)
               MBASQ=MBAS(ISQB)
C              ***** DETERMINE LOOP LIMITS FOR FIRST INDEX A *****
               IF (ITRLVL .GE. 5) THEN
C                 ... general index on A and B.
                  IF (ISPA.LT.ISQB) GO TO 150
                  NARMIN = 1
                  NARMAX = MORBA
               ELSE IF (IC.LE.MOCCT) THEN
C              ***** OCCUPIED-OCCUPIED PAIR CD *****
                  IF(ISPA.LT.ISQB.AND.ITRLVL.NE.0) GO TO 150
C                 ... A .ge. B condition,
C                     for ITRLVL=0  A may be .lt. B when A secondary.
C                     (870628-hjaaj)
                  IF (IC.LE.NXISHT) THEN
C                 (  we now know IC .eq. ID,
C                    this distribution is for (ii/aa) integrals )
                     IF (ITRLVL.EQ.3 .OR. ISPA.NE.ISC) GO TO 150
C                  ( we only need (ii/aa) with sym(a) .eq. sym(i) )
C                  ( if ITRLVL.eq.3, neglect (ii/aa) integrals )
                     NARMIN = MOCCA+1
                     NARMAX = MORBA
                     NBRMIN = NARMIN
                     NBRMAX = NARMAX
                     GO TO 145
                  ELSE
                     NARMIN=1
                     NARMAX=MORBA
                  END IF
               ELSE
C              ***** VIRTUAL-OCCUPIED PAIR CD *****
CHJ-S-841125
                  IF (ID.LE.NXISHT) THEN
C                 ( symmetry tested above, this is (ai/,
C                   used for (ai/ai) )
C                 ( if ITRLVL.eq.3, neglect (ai/ai) integrals )
                     IF (ITRLVL.EQ.3 .OR. ISPA.NE.ISC) GO TO 150
                     NARMIN = NCR
                     NARMAX = NCR
                     NBRMIN = NDR
                     NBRMAX = NDR
                     GO TO 145
                  ELSE IF (ISPA.LT.ISQB) THEN
C                 ( (au/ case )
                     NARMIN = MOCCA+1
                  ELSE
                     NARMIN = 1
                  END IF
CHJ-E-841125
                  NARMAX=MORBA
#if defined (VAR_OLDCAS)
                  IF (ITRLVL.LT.2) THEN
                     IF (ISPA.LT.ISC) GO TO 150
                     IF (ISPA.EQ.ISC) NARMIN=NCR
                  END IF
#endif
                  IF(NARMIN.GT.NARMAX) GO TO 150
               END IF
C
C     *****SET ABSOLUTE LOOP LIMITS FOR B *****
               IF (ITRLVL.GE.5) THEN
                  NBRMIN = 1
                  NBRMAX = MORBB
               ELSE IF (ITRLVL.EQ.0) THEN
                  NBRMIN = MISH(ISQB) + 1
                  NBRMAX = MOCCB
               ELSE
                  NBRMIN = 1
                  NBRMAX = MORBB
                  IF (IC.GT.MOCCT) NBRMAX = MOCCB
               END IF
               IF (NBRMAX.LT.NBRMIN) GO TO 150
CHJ-S-841125
  145          CONTINUE
CHJ-E-841125
               NNBR = 1 + NBRMAX - NBRMIN
C     ***** NUMBER OF ORBITALS OF PRECEDING SYMMETRIES *****
               JORBA=JORB(ISPA)
               JORBB=JORB(ISQB)
               JBASP=JBAS(ISPA)
               JBASQ=JBAS(ISQB)
               NSAB = IT(MAX(ISPA,ISQB)) + MIN(ISPA,ISQB)
C     ***** START ADDRESSES FOR MO*S *****
C     ***** SET THE ARRAY PBCD TO ZERO *****
               IMOP = JCMO(ISPA) + NARMIN-1 + (NP-JBASP-1)*MORBP
               IMOQ = JCMO(ISQB) + NBRMIN
               CALL DZERO(PBCD(NBRMIN),NNBR)
               NCDITP = NTCD + IT(NP)
               NCDP   = NTCD + NP
C     ***** LOOP OVER RELATIVE INDICES Q IN SYMMETRY ISQB *****
      DO 130 NQ = (JBASQ+1),(JBASQ+MBASQ)
C        ***** FIND POSITION OF INTEGRAL (PQ/CD) *****
         IF (NP.GE.NQ) THEN
            IPQCD = NCDITP + NQ
         ELSE
            IPQCD = NCDP   + IT(NQ)
         END IF
         PQCD=TRI(IPQCD)
C        ***** LOOP OVER MO'S OF SYMMETRY ISQB *****
C        ***** AND TRANSFORM SECOND INDEX      *****
         IF (ABS(PQCD) .GE. THRP)
     *      CALL DAXPY(NNBR,PQCD,CMOT(IMOQ),1,PBCD(NBRMIN),1)
         IMOQ = IMOQ + MORBQ
  130 CONTINUE
C     ***** PBCD NOW CONTAINS INTEGRALS (PB/CD) FOR A *****
C     ***** GIVEN TRIPLE PCD AND ALL APPROPRIATE B    *****
C     ***** LOOP OVER ORBITALS A AND ADD CONTRIBUTIONS*****
C     ***** TO INTEGRALS (AB/CD) FOR THIS TRIPLE      *****
      DO 140 NAR=NARMIN,NARMAX
         IMOP = IMOP + 1
         CMOP = CMOT(IMOP)
      IF (ABS(CMOP) .LE. THRP) GO TO 140
C        ***** SET ABSOLUTE INDEX *****
         NA   = NAR  + JORBA
C        ***** SET LOOP LIMITS FOR ORBITAL B *****
         NBRS = NBRMIN
         NBRE = NBRMAX
         IF (ITRLVL.GE.2) THEN
            IF (ISPA.EQ.ISQB) NBRE = MIN(NBRE,NAR)
         ELSE IF (ITRLVL.EQ.0) THEN
            IF (ITSW(NA).LE.MASHT) THEN
C              ... A is active
               IF(ISPA.LT.ISQB) GO TO 140
#if defined (VAR_OLDCAS)
               IF(NSAB.LT.NSCD) GO TO 140
               IF(ISPA.EQ.ISC.AND.NAR.LT.NCR)GO TO 140
               IF(NA.EQ.NC)NBRS=NDR
#endif
               IF(ISPA.EQ.ISQB)NBRE=NAR
            END IF
         ELSE
            IF(NSAB.LT.NSCD.AND.NAR.LE.MOCCA) NBRS=MOCCB+1
            IF(NA.EQ.NC) NBRS=NDR
            IF(ISPA.EQ.ISC.AND.NAR.LT.NCR) NBRS=MOCCB+1
            IF(ISPA.EQ.ISQB.AND.IC.LE.MOCCT) NBRE=NAR
C          (HJ: why IC.le.MOCCT? I have removed it for ITRLVL.eq.2 or 3;
C             IA.le.MOCCT would make more sense)
         END IF
CHJ-S-841125
         IF (IC.LE.NXISHT) THEN
C        ( (ii/aa) case, symmetry etc. has already been checked)
            NBRS = NAR
            NBRE = NAR
         END IF
CHJ-E-841125
C        ***** LOOP OVER ORBITALS B AND ADD CONTRIBUTION *****
C        ***** TO INTEGRALS (AB/CD)                      *****
         DO 135 NBR = NBRS,NBRE
!           that could be done with daxpy, could not it?
            ABCD(JORBB+NBR,NA) = ABCD(JORBB+NBR,NA) + CMOP*PBCD(NBR)
  135    CONTINUE
  140 CONTINUE
  150 CONTINUE
C     ***** ALL INTEGRALS (AB/CD) FOR A GIVEN PAIR CD HAVE *****
C     ***** BEEN CREATED. WRITE THEM ON UNIT LUINTM WITH   *****
C     ***** INDICES A,B,C AND D.                           *****
CHJ: new index order: C,D; A,B
      DO 170 NA=1,MORBT
         INDCDA = IBTOR(INDCD , IBTSHL( NA, 8 ) )
         ISA = ITSMO( NA )
         ISB = MUL( ISA, ISCD )
         NBST  = JORB(ISB) + 1
         NBEND = JORB(ISB) + MORB(ISB)
C        ... not NBEND = MIN(NA,NBEND)
C            because e.g. (ph/ph) integrals do not follow NB .le. NA
         DO 160 NB = NBST,NBEND
            P = ABCD(NB,NA)
         IF (ABS(P).LT.THRP) GO TO 160
            LABCD=LABCD+1
            IOUT =IOUT+1
            IF (IPRTRA .GE. 200) THEN
               WRITE (LUPRI,1010)LABCD,NC,ND,NA,NB,P
 1010          FORMAT(I10,4I5,F20.15)
            END IF
            IF (IOUT.GT.LOUT) THEN
               IOUT=1
               IOUTB(LOUT21)=LOUT
C              ***** WRITE THIS BUFFER *****
#if defined (VAR_SPLITFILES)
               CALL WRITSI(LUINTM,LOUTI,IOUTB,JBUF)
               CALL WRITSI(LUINTM,LOUT+2,IOUTB(LOUTI + 1),JBUF)
#else
               CALL WRITI(LUINTM,LOUT22,IOUTB)
#endif
               NRINTM = NRINTM + 1
               IF (IPRTRA .GE. 100 .AND. IPRTRA .LT. 200) THEN
                  WRITE(LUPRI,*)
     &            ' Writing buffer for NC = ',NC,' ND = ',ND
                  WRITE(LUPRI,*) ' NA = ',NA,' NB = ',NB
                  WRITE(LUPRI,*) 'Buffer no.',NRINTM
                  WRITE(LUPRI,*) ' First 10 integrals'
                  WRITE(LUPRI,'(5F14.7)') (OUTB(I),I=1,10)
               END IF
            END IF
            OUTB(IOUT) = P
            IOUTB(LOUTI+IOUT) = IBTOR( INDCDA , NB )
  160    CONTINUE
  170 CONTINUE
CHJ 3 *** Going to next (CD), empty this buffer ***
      IOUTB(LOUT21) = IOUT
#if defined (VAR_SPLITFILES)
      CALL WRITSI(LUINTM,LOUTI,IOUTB,JBUF)
      CALL WRITSI(LUINTM,LOUT+2,IOUTB(LOUTI + 1),JBUF)
#else
      CALL WRITI(LUINTM,LOUT22,IOUTB)
#endif
      NRINTM = NRINTM + 1
      IF (IPRTRA .GE. 100 .AND. IPRTRA .LT. 200) THEN
         WRITE(LUPRI,*)
     &   ' Writing last buffer for NC = ',NC,' ND = ',ND
         WRITE(LUPRI,*) 'Buffer no.',NRINTM
         WRITE(LUPRI,*) ' First 10 integrals, length =',IOUT
         LPRI = MIN(10,IOUT)
         WRITE(LUPRI,'(5F14.7)') (OUTB(I),I=1,LPRI)
      END IF
      IOUT = 0
  100 CONTINUE
C     ***** WRITE LAST BUFFER *****
 8000 CONTINUE
      IOUTB(LOUT21)=-1
#if defined (VAR_SPLITFILES)
      CALL WRITSI(LUINTM,LOUTI,IOUTB,JBUF)
      CALL WRITSI(LUINTM,LOUT+2,IOUTB(LOUTI + 1),JBUF)
#else
      CALL WRITI(LUINTM,LOUT22,IOUTB)
#endif
      WRITE(LUINTM)LAB123,TABLE(2)
      CALL REWSPL(LUINTM)
      DFIN=SECOND()
      DTIM=DFIN-DST
      IF (JTER.EQ.1) THEN
         WRITE(LUPRI,49)DTIM
         WRITE(LUPRI,1000) LUINTM,LABCD,NRINTM,LBINTM
         WRITE(LUPRI,1100) THRP
      END IF
   49 FORMAT(' TIME IN TRAAB',F10.2)
 1000 FORMAT(/,' NUMBER OF TRANSFORMED INTEGRALS WRITTEN ON UNIT'
     1       ,I3,' IS',I10,
     2       /,' USING',I6,' + 1 RECORDS WITH BUFFER LENGTH',I6)
 1100 FORMAT(/'  THRESHOLD FOR DISCARDING INTEGRALS',1P,D10.2)
      CALL QEXIT('TRAAB ')
      RETURN
      END
C  /* Deck tractl */
      SUBROUTINE TRACTL(KTRLVL,CMO,WORK,KWORK)
C
C Revisions
C   050125-hjaaj: new TRACTL_1 for dynamic allocation
C
C Input:
C  KTRLVL, abs(KTRLVL) is transformation level, DELDA = (KTRLVL.lt.0)
C  CMO,   mo coefficients
C  WORK,  work array
C  KWORK, actual length of work array
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*), WORK(KWORK)
#include <iratdef.h>
#include <thrzer.h>
#include <lbmxsq.h>
C
C Used from common blocks:
C
#include <priunit.h>
#include <inforb.h>
#include <exeinf.h>
#include "inttra.h"
#include <inftra.h>
#include <inftap.h>
#include <gnrinf.h>
C
      DIMENSION MFRO(8)
      LOGICAL VFIRST, FOPEN, DELDA, FEXIST
      SAVE VFIRST, IPRTRA_SAVE, LDBUF, LRBUF, ITRSAV
      DATA VFIRST /.TRUE./
C
C     VFIRST is is used to control setting of print level, see below.
C
C     MAKE SURE BLOCK DATA MODULE FOR INFTRA IS INCLUDED :
C
      EXTERNAL BDTRA
C
      IF (VFIRST) THEN
C        if SIRIUS has been called, IPRTRA_SAVE will have value
C        from SIRIUS in all TRACTL calls, if SIRIUS has
C        not been called before ABACUS is called, then
C        the value specified in RHSINP for DERTRA is
C        used. This makes sure that the SIRIUS value is
C        always used if it has been defined, and the a
C        value is defined otherwise.
C
         IPRTRA_SAVE = IPRTRA
         VFIRST = .FALSE.
      ELSE
         IPRTRA = IPRTRA_SAVE
      END IF
C
C
      IF (NEWTRA) THEN
         CALL N_TRACTL(KTRLVL,CMO,WORK,KWORK)
      ELSE
         KFRSAV = 1
         KFREE  = 1
         LFREE  = KWORK
         CALL MEMGET('INTE',KITSMO,NBAST,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KITSAO,NBAST,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KITSW ,NBAST,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KITSX ,NBAST,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KIT   ,NBAST+1,WORK,KFREE,LFREE)
         CALL TRACTL_1(KTRLVL,CMO,WORK(KITSMO),WORK(KITSAO),
     &        WORK(KITSW),WORK(KITSX),WORK(KIT),WORK(KFREE),LFREE)
         CALL MEMREL('TRACTL_1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      END IF
      RETURN
      END
C  /* Deck tractl_1 */
      SUBROUTINE TRACTL_1(KTRLVL,CMO,ITSMO,ITSAO,ITSW,ITSX,IT,
     &                    WORK,KWORK)
C
C Revisions
C   18-Apr-1984 hjaaj / 4-Dec-1983 hjaaj
C   13-Dec-1984 hjaaj (ITRLVL=3 option)
C   28-Aug-1986 hjaaj (test on LTRI, ITRLVL=4 option)
C   18-May-1987 hjaaj (ITRLVL=5 option)
C   19-Dec-1988 hjaaj (ITRLVL,CMO in parameter list)
C   10-Feb-1989 hjaaj (ITRLVL=10, full integral transf.)
C   890615-hjaaj disabled if (first) for orb.spec., so it can be used
C                in combination runs.
C   900801-hjaaj: added CMO to LUINTM
C   050125-hjaaj: TRACTL divided in TRACTL and TRACTL_1 for dynamic allocation
C                 of IT, ITSMO, ITSAO, ITSW, ITSX
C
C Input:
C  KTRLVL, abs(KTRLVL) is transformation level, DELDA = (KTRLVL.lt.0)
C  CMO,   mo coefficients
C  WORK,  work array
C  KWORK, actual length of work array
C  MWORK, desired length of work array
C         (dependent on size of physical memory)
C
C  FTRCTL in exeinf.h forces integral sort and transformation,
C      FTRCTL is true if and only if Hermit has calculated
C      a new AOTWOINT file since last call of TRACTL.
C
C
C  SUBROUTINE CALLS:
C          SORTA (BIN SORT OF INTEGRALS (PQ/RS))
C          TRACD (TRANSFORMATION OF (PQ/RS) INTO (PQ/CD))
C          SORTB (BIN SORT OF INTEGRALS (PQ/CD))
C          TRAAB (TRANSFORMATION OF (PQ/CD) INTO (AB/CD))
C          MOLLAB
C
C  Based on:
C     CASSCF:TRANSFORMATION SECTION.CONTROL ROUTINE
C      ********** RELEASE 79 04 20 **********
C
#include <implicit.h>
#include <dummy.h>
      INTEGER   ITSMO(*), ITSAO(*), ITSW(*), ITSX(*), IT(*)
      DIMENSION CMO(*), WORK(KWORK), MISH_TEST(8), MASH_TEST(8)
#include <iratdef.h>
#include <thrzer.h>
#include <lbmxsq.h>
      PARAMETER (DM1 = -1.0D0)
C
C Used from common blocks:
C CBGETDIS : IAD*
C INFDIM   : MWORK
C
#include <priunit.h>
#include <inforb.h>
#include <exeinf.h>
#include "infdim.h"
#include "inttra.h"
#include <inftra.h>
#include <inftap.h>
#include <gnrinf.h>
#include <cbgetdis.h>
C
      DIMENSION MFRO(8)
      LOGICAL FOPEN, DELDA, FEXIST
      SAVE LDBUF, LRBUF, ITRSAV
C
      CALL QENTER('TRACTL_1')
      ITRLVL = ABS(KTRLVL)
      DELDA  = (KTRLVL .LT. 0)
      IF (IPRTRA .GT. 5) THEN
         WRITE(LUPRI,*) 'TRACTL_1(KTRLVL,CMO,WORK,KWORK,MWORK)'
         WRITE(LUPRI,*) '====================================='
         WRITE(LUPRI,*) 'KTRLVL = ', KTRLVL
         WRITE(LUPRI,*) 'KWORK  = ', KWORK
         WRITE(LUPRI,*) 'MWORK  = ', MWORK
      END IF
      IF (IPRTRA .GE. 20) THEN
         WRITE(LUPRI,*) ' CMO array '
         CALL PRORB(CMO,.FALSE.,LUPRI)
      END IF
C
      IF (FTRCTL) THEN
         ITRSAV = 999
         LSRTAO = 0
C        ... force new integral sort and transformation
         FTRCTL = .FALSE.
      END IF
      IF (ITRLVL .NE. ITRSAV) THEN
         ITRSAV = ITRLVL
         IF (IPRTRA .GT. 1) THEN
            JTER   = 1
         ELSE
            JTER   = 2
         END IF
      END IF
      LWORK = MIN(KWORK,MWORK)
      IF (JTER.EQ.1) WRITE(LUPRI,903) LWORK
  903 FORMAT(/' WORKING AREA IN TRANSFORMATION SECTION ',I10)
C     ***** SET THRESHOLDS *****
C     *********** THRP: THRESHOLD USED IN TRACD AND TRAAB
C     *********** THRQ: THRESHOLD USED IN SORTA FOR INTEGRALS (PQ/RS)
      IF (THRP.LE.0.D0) THEN
         WRITE(LUPRI,'(/A,1P,D10.2,A,D10.2)')
     &      ' INFO: Changed THRP threshold used in transformation from',
     &      THRP,' to',1.D-15
         THRP=1.D-15
      END IF
      IF (THRQ.LE.0.D0) THEN
         WRITE(LUPRI,'(/A,1P,D10.2,A,D10.2)')
     &      ' INFO: Changed THRQ threshold used in AO sort from',THRP,
     &      ' to',1.D-15
         THRQ=1.D-15
      END IF
C
C     ***** Reset /CBGETD/ because any H2AC on disk will be
C           obsolete after transformation.
C
      IADINT = -1
      IADH2  = -1
      IADH2X = -1
      IADH2D = -1
      LUHALF = -1
      LUDA   = -1
      LUDA2  = -1
C
C     ***** TYPE OF TRANSFORMATION                             *****
C     ***** ITRLVL=0:FIRST INDEX,ALL ORBITALS                  *****
C     *****         OTHER INDICES,ACTIVE ORBITALS ONLY         *****
C                   (CAS-SCI, CI, MC gradient transformation)
C     ***** ITRLVL=1:TWO INDICES,ALL ORBITALS                  *****
C     *****         TWO INDICES,PRIMARY ORBITALS ONLY          *****
C                   (CAS-NR transformation)
C
C           ITRLVL=2:SIRIUS transformation
C           ITRLVL=3:SIRIUS transformation, neglect (ii/aa) and (ia/ia)
C           ITRLVL=4:SIRIUS transformation, include (ij/ab) and (ia/jb)
C
C           ITRLVL=5:Three indices, all orbitals
C                    One   index,   occupied orbitals
C
C           ITRLVL=10: Full integral transformation
C
      IF (IPRTRA.GT.0) WRITE (LUPRI,1000) ITRLVL
 1000 FORMAT(/,' TRANSFORMATION SECTION: PARAMETER =',I4)
      IF (ITRLVL.GT.5 .AND. ITRLVL .NE. 10) THEN
         WRITE (LUPRI,'(//A)')
     1      ' TRACTL_1, parameter outside valid range'
         CALL QTRACE(LUPRI)
         CALL QUIT('ERROR, INVALID TRANSFORMATION PARAMETER (TRACTL_1)')
      END IF
C
C     ***** Assign INTTRA from INFORB *****
C     Note: index on LUITM are wrong if MFRO(ISYM) .ne. 0 !!
      MSYM = NSYM
      DO 5 ISYM = 1,NSYM
         MFRO(ISYM) = 0
         MISH(ISYM) = NISH(ISYM) - MFRO(ISYM)
         MASH(ISYM) = NASH(ISYM)
         MSSH(ISYM) = NSSH(ISYM)
         MBAS(ISYM) = NBAS(ISYM)
    5 CONTINUE
      IF (IPRTRA .GE. 5) THEN
         WRITE(LUPRI,*) 'MSYM',MSYM
         WRITE(LUPRI,*) 'MFRO',MFRO(1:8)
         WRITE(LUPRI,*) 'MISH',MISH(1:8)
         WRITE(LUPRI,*) 'MASH',MASH(1:8)
         WRITE(LUPRI,*) 'MSSH',MSSH(1:8)
         WRITE(LUPRI,*) 'MBAS',MBAS(1:8)
      END IF
C
C     If first call of TRACTL_1 set up orbital information in /INTTRA/
C        ... 890615 - hjaaj: in new combination runs
C            (e.g. RHF-MP2-CI) number of active orbitals will change
C            therefore we now must recalculate this information
C            each time.
C
C     ***** SET UP ORBITAL DATA *****
         MOCCT=0
         MSSHT=0
         MASHT=0
         MBAST=0
         NMORBT=0
         NMBAST=0
         DO 10 ISYM=1,NSYM
            MOCC(ISYM)=MISH(ISYM)+MASH(ISYM)
            MORB(ISYM)=MOCC(ISYM)+MSSH(ISYM)
            MOCCT=MOCCT+MOCC(ISYM)
            MSSHT=MSSHT+MSSH(ISYM)
            MASHT=MASHT+MASH(ISYM)
            MBAST=MBAST+MBAS(ISYM)
            NMORBT=NMORBT+(MORB(ISYM)**2+MORB(ISYM))/2
            NMBAST=NMBAST+(MBAS(ISYM)**2+MBAS(ISYM))/2
   10    CONTINUE
C
         MORBT=MOCCT+MSSHT
         IF (MORBT .GT. MBAST) CALL QUIT('MORBT .gt. MBAST :(')
CHJ-S-841208
         IF (MORBT.GT.255 .OR. MBAST.GT.255) THEN
            WRITE (LUPRI,*)
     *         'TRACTL_1-ERROR; TRACTL_1 is limited to 255 orbitals'
            WRITE (LUPRI,*)
     *         'NORBT =',MORBT,', NBAST =',MBAST
            CALL QTRACE(LUPRI)
            CALL QUIT('TRACTL_1: TOO MANY ORBITALS')
         END IF
CHJ-E-841208
         NMBASX=MBAST*(MBAST+1)/2
         NMORBX=MORBT*(MORBT+1)/2
         NMASHX=MASHT*(MASHT+1)/2
         M2ORBX=MORBT*MORBT
C        ***** NUMBER OF ORBITALS AND BASIS FUNCTIONS OF EARLIER *****
C        ***** SYMMETRIES. START ADDRESSES FOR MO*S OF GIVEN     *****
C        ***** SYMMETRY.                                         *****
         JBAS(1) = 0
         JORB(1) = 0
         JCMO(1) = 0
         DO 20 ISYM = 2,NSYM
            ISYM1 = ISYM - 1
            JBAS(ISYM) = JBAS(ISYM1) + MBAS(ISYM1)
            JORB(ISYM) = JORB(ISYM1) + MORB(ISYM1)
            JCMO(ISYM) = JCMO(ISYM1) + MORB(ISYM1)*MBAS(ISYM1)
   20    CONTINUE
C        ***** COMPUTE SYMMETRY INDEX ITSAO and ITSMO *****
         II = 0
         JJ = 0
         DO 30 ISYM = 1,NSYM
            MBASI = MBAS(ISYM)
            MORBI = MORB(ISYM)
            DO 25 I = 1,MBASI
               II = II + 1
               ITSAO(II) = ISYM
            IF (I.GT.MORBI) GO TO 25
               JJ = JJ + 1
               ITSMO(JJ) = ISYM
   25       CONTINUE
   30    CONTINUE
C
C        MCMOT : SIZE OF THE MO-SPACE
C        NNOCX : NUMBER OF OCCUPIED PAIRS
C        MAXCHA: MAXIMUM NUMBER OF DA CHAINS
C        LDAMAX: MAXIMUM SIZE OF DA RECORDS
C        LRBUF : Max size of AO/MO integral records
C
         MCMOT = 0
         DO 50 ISYM = 1,NSYM
   50       MCMOT = MCMOT + MORB(ISYM)*MBAS(ISYM)
         NMOCCX = MOCCT*(MOCCT+1)/2
C
         LDAMAX = MAX(1706,NMBAST)
         LDAMAX = ((LDAMAX+IRAT-1)/IRAT)*IRAT
C        ... 1706 gives buffer length of 20kB when IRAT=2
         LDBUF  = (LDAMAX/IRAT)*(IRAT+1) + (1+1/IRAT)
         LRBUF  = (LBMXSQ/IRAT)*(IRAT+1) + (1+1/IRAT)
C
C        All orbital information now found and written in /INTTRA/.
C        ... 890615 - hjaaj: in new combination runs
C            (e.g. RHF-MP2-CI) number of active orbitals will change
C            therefore we now must recalculate this information
C            each time.
C
C     ***** REORDERING INDICES FOR MOLECULAR ORBITALS   *****
C     ***** ITSW REORDERS TO PRIMARY-SECONDARY ORDERING *****
C     ***** ITSX REORDERS BACK                          *****
C
      IF (ITRLVL.EQ.0) GO TO 42
      IF (ITRLVL.EQ.10) GO TO 43
      IF (ITRLVL.EQ.1 .OR. ITRLVL.GE.4) GO TO 199
C
C     --- ITRLVL.EQ.2 .OR. ITRLVL.EQ.3:
C
      I    = 0
      IO   = 0
      IACT = MOCCT - MASHT
      ISEC = MOCCT
      DO 180 ISYM = 1,NSYM
        MISHI = MISH(ISYM)
        DO 120 K = 1,MISHI
          I  = I  + 1
          IO = IO + 1
          ITSW(I)  = IO
          ITSX(IO) = I
  120   CONTINUE
        MASHI = MASH(ISYM)
        DO 140 K = 1,MASHI
          I    = I    + 1
          IACT = IACT + 1
          ITSW(I)    = IACT
          ITSX(IACT) = I
  140   CONTINUE
        MSSHI = MSSH(ISYM)
        DO 160 K = 1,MSSHI
          I    = I    + 1
          ISEC = ISEC + 1
          ITSW(I)    = ISEC
          ITSX(ISEC) = I
  160   CONTINUE
  180 CONTINUE
      GO TO 44
C
C     --- ITRLVL.EQ.1 .OR. ITRLVL.EQ.4 .OR. ITRLVL.EQ.5
C
  199 CONTINUE
      I    = 0
      ISEC = MOCCT
      IO   = 0
      DO 40 ISYM = 1,NSYM
         MOCCI = MOCC(ISYM)
         DO 32 KOCC = 1,MOCCI
            I  = I  + 1
            IO = IO + 1
            ITSW(I)  = IO
            ITSX(IO) = I
   32    CONTINUE
         MSSHI = MSSH(ISYM)
         DO 37 KSSH = 1,MSSHI
            I    = I    + 1
            ISEC = ISEC + 1
            ITSW(I)    = ISEC
            ITSX(ISEC) = I
   37    CONTINUE
   40 CONTINUE
      GO TO 44
C
C     --- ITRLVL.EQ.0:
C
   42 I    = 0
      ISEC = MASHT
      IO   = 0
      DO 400 ISYM = 1,NSYM
         MISHI = MISH(ISYM)
         DO 320 KISH = 1,MISHI
            ISEC = ISEC + 1
            I    = I    + 1
            ITSW(I)    = ISEC
            ITSX(ISEC) = I
  320    CONTINUE
         MASHI = MASH(ISYM)
         DO 340 KASH = 1,MASHI
            IO=IO+1
            I=I+1
            ITSW(I)=IO
            ITSX(IO)=I
  340    CONTINUE
         MSSHI=MSSH(ISYM)
         DO 360 KSSH=1,MSSHI
            ISEC=ISEC+1
            I=I+1
            ITSW(I)=ISEC
            ITSX(ISEC)=I
  360    CONTINUE
  400 CONTINUE
      GO TO 44
C
   43 CONTINUE
C     ... ITRLVL .EQ. 10
      DO 430 I = 1,MORBT
         ITSW(I) = I
         ITSX(I) = I
  430 CONTINUE
C
   44 CONTINUE
C
C     ***** SET UP DYNAMIC STORAGE FOR TRACD *****
C
      LW1   = 1
      LW2   = LW1 + MCMOT
      LW3   = LW2 + LDBUF
      LW4   = LW3 + LRBUF
      LW5   = LW4 + M2ORBX
      LW6   = LW5 + MORBT
      LTRI  = LWORK - LW6 + 1
      NPQMIN = 1 + (NMBASX-1)/MAXCHA
      LTRIMIN = NPQMIN*NMBASX
      IF (LTRI.LT.LTRIMIN) THEN
         LWORK  = LW6 - 1 + LTRIMIN
         IF (LWORK.LE.KWORK) THEN
            LTRI  = LTRIMIN
            MWORK = MAX(LWORK,MWORK)
            WRITE(LUPRI,9030) LWORK
         ELSE
            WRITE(LUPRI,9000) LTRI
            WRITE(LUPRI,'(/A,6I9)')
     *      ' LW1, LW2, ..., LW6 :',LW1,LW2,LW3,LW4,LW5,LW6
            CALL QTRACE(LUPRI)
            CALL QUIT('TRACTL_1: INSUFFICIENT SPACE FOR TRACD')
         END IF
 9030    FORMAT(/' TRACTL_1, not enough work space for TRACD,',
     *          /'         work space increased to',I10/)
 9000    FORMAT(/' WORK SPACE IN TRACD,',I10,', IS MUCH TOO SMALL.')
      END IF
C
C     Allocate sorting area for SORTA/SORTB.
C     hjaaj nov 2001: Now allow SORTA/SORTB to use all memory
C     (i.e. KWORK) for sorting (LWORK is the memory available
C     for TRACD, based on MWORK))
C
      KWORK1 = 1 + LRBUF
      LWORK1 = KWORK - KWORK1
C
C calculate LPQRS,  max. number of ao integrals
C
      LPQRS = 0
      DO 11 I = 1,NSYM
         DO 12 J = 1,I
            NSIJ = MUL(I,J)
            NBIJ = MBAS(I)*MBAS(J)
            IF (I.EQ.J) NBIJ = MBAS(I)*(MBAS(I)+1)/2
            DO 13 K = 1,I
               LMAX = K
               IF (K.EQ.I) LMAX = J
               DO 14 L = 1,LMAX
                  NSKL = MUL(K,L)
               IF(NSKL.NE.NSIJ)GO TO 14
                  NBKL = MBAS(K)*MBAS(L)
                  IF(K.EQ.L)NBKL = MBAS(K)*(MBAS(K)+1)/2
                  ITERM = NBIJ*NBKL
                  IF (I.EQ.K.AND.J.EQ.L) ITERM = NBIJ*(NBIJ+1)/2
                  LPQRS = LPQRS + ITERM
   14          CONTINUE
   13       CONTINUE
   12    CONTINUE
   11 CONTINUE
      IF (JTER.EQ.1) WRITE (LUPRI,902) LPQRS
  902 FORMAT(' MAXIMUM NUMBER OF AO INTEGRALS = ',I10)
C
      NBINTM = (NMORBT-1)/LBMXSQ + 1
      LBINTM = (NMORBT-1)/NBINTM + 1
      LBINTM = MAX(4,LBINTM)
C     ... max distributions for any symmetry =
C         number of distributions of symmetry 1 = NMORBT (890215-hjaaj)
C
C     ***** CALLING SEQUENCE FOR FIRST SORT *****
C
C....    LOOP TO SET UP TRIANGULAR INDEXING ARRAY
C
      ITLIM = MBAST + 1
      II = 0
      DO 45 I = 1,ITLIM
         IT(I) = II
         II = II + I
   45 CONTINUE
C
      DST = SECOND()
C
C     900801-hjaaj: check if transformation already done
C     971201-hjaaj: code moved to here (was too late)
C        LSRTAO .ne. 0 should cover all situations where
C        it is possible that MO integrals already exist
C        on LUINTM. NOTE: this requires that SIR_INTOPEN has been
C        called before first call to TRACTL_1 !!!!!!
C        SIR_INTOPEN will set LSRTAO = -1 if LUINTM exists with
C        a 'MOLTWOEL' label. Whenever TRACTL_1 has performed
C        a transformation LSRTAO will be = +1.
C     
C   
      IF (LSRTAO .NE. 0) THEN
         CALL REWSPL(LUINTM)
         READ  (LUINTM) NPQ,LDABUF,(LASTAD(I),I=1,2500)
         READ  (LUINTM) LBINT1,JTRLVL
         IF (JTRLVL .EQ. 3 .AND. ITRLVL .EQ. 2) JTRLVL = -1
         IF (JTRLVL .EQ. 2 .AND. ITRLVL .EQ. 3) JTRLVL =  3
         IF (JTRLVL .GE. ITRLVL) THEN
C           read CMO matrix from LUINTM
C           and subtract from input CMO matrix
            READ  (LUINTM) (WORK(I),I=1,NCMOT)
            CALL DAXPY(NCMOT,DM1,CMO,1,WORK,1)
            I = IDAMAX(NCMOT,WORK,1)
            READ (LUINTM) MISH_TEST(1:8), MASH_TEST(1:8)
            MISH_TEST(1:8) = MISH_TEST(1:8) - MISH(1:8)
            MASH_TEST(1:8) = MASH_TEST(1:8) - MASH(1:8)
            J = 0
            DO ISYM = 1,8
               J = J + ABS(MISH_TEST(ISYM)) + ABS(MASH_TEST(ISYM))
            END DO
            IF (ABS(WORK(I)) .LE. THRZER .AND. J .EQ. 0) THEN
               IF (IPRTRA.GE.0) WRITE (LUPRI,'(/A/A)')
     *            ' TRACTL_1: Integral transformation abandoned,',
     *            ' the required MO integrals are already available.'
               LBINTM = LBINT1
               GO TO 9999
            END IF
         END IF
      END IF
      IF (LSRTAO .LT. 0) THEN
C     ... check if old LUDA exist, if yes then check if file
C         is consistent with LUINTM information /890302-hjaaj
         CALL REWSPL(LUINTM)
         READ  (LUINTM) NPQ,LDABUF,(LASTAD(I),I=1,2500)
         CALL REWSPL(LUINTM)
         LBUF = (IRAT + 1)*LDABUF + 2
         CALL GPINQ('AOSRTINT.DA','EXIST',FEXIST)
         IF (FEXIST) THEN
            IF (LUDA .LE. 0) CALL GPOPEN(LUDA,'AOSRTINT.DA','OLD',
     &                                   'DIRECT',' ',LBUF,OLDDX)
C           LUDA exists, try to read last record acc. to LUINTM ...
            IPQ  = (NMBASX-1)/NPQ + 1
            IADR = LASTAD(1)
            DO 110 I = 2,IPQ
               IADR = MAX(IADR,LASTAD(IPQ))
  110       CONTINUE
            LBUF = LBUF/IRAT
            READ (LUDA,REC=IADR,IOSTAT=IODA) (WORK(I),I=1,LBUF)
            IF (IODA .NE. 0) THEN
C        ... LUDA not consistent with LUINTM information
               CALL GPCLOSE(LUDA,'DELETE')
               LSRTAO = 0
            ELSE
               LSRTAO = 1
               WRITE(LUPRI,'(/A)')
     *              ' Old direct access AO file found, SORTA skipped.'
            END IF
         ELSE
C        ... LUDA does not exist
            LSRTAO = 0
         END IF
      END IF
      IF (LSRTAO.EQ.0) THEN
C        hjaaj 11-jan-2007: make sure AOTWOINT exists
         CALL MAKE_AOTWOINT(WORK,KWORK)

C     HJ-840710 was: IF (JTER.EQ.1) THEN
         CALL GPINQ(FNINTM,'EXIST',FEXIST)
         if(FEXIST) CALL GPCLOSE(LUINTM,'DELETE')
         CALL GPOPEN(LUINTM,FNINTM,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
C        ... empty old mo integral file to save disk space during
C            the transformation
C        Now done by deleting the file, just in case it has been split
C        K.Ruud, April 19 2000
C
         CALL SORTA(LUDA,WORK,WORK,WORK(KWORK1),WORK(KWORK1),
     *              LWORK1,LTRI,IT)
C        CALL SORTA(RINT,IINT,SCR,ISCR,MEMSX,MEMTX,IT)
         DFIN = SECOND()
         DTIM = DFIN - DST
         DST  = DFIN
         IF (IPRTRA.GT.0) WRITE (LUPRI,46) DTIM
   46    FORMAT(/' TIME IN SORTA',F10.2)
         CALL REWSPL(LUINTM)
         WRITE (LUINTM) NPQ,LDABUF,(LASTAD(I),I=1,2500)
         CALL REWSPL(LUINTM)
         LSRTAO = 1
      ELSE
         CALL REWSPL(LUINTM)
         READ  (LUINTM) NPQ,LDABUF,(LASTAD(I),I=1,2500)
         CALL REWSPL(LUINTM)
C
         LTRIPQ = NPQ*NMBASX
         LBUF = (IRAT + 1)*LDABUF + 2
         IF (LTRI .LT. LTRIPQ) THEN
C           This will only happen if TRACTL_1 is called with (much) less
C           work space than was available when SORTA was called.
            LTRI = KWORK - LW6 + 1
            IF (LTRI .LT. LTRIPQ)
     *      CALL ERRWRK('TRACTL_1.TRACD, too large chains',LTRIPQ,LTRI)
         END IF
         IF (LUDA .LE. 0) CALL GPOPEN(LUDA,'AOSRTINT.DA','OLD',
     &                                'DIRECT','UNFORMATTED',LBUF,OLDDX)
      END IF
C
      DSTTRA = DST
C
      CALL GPOPEN(LUHALF,' ',' ',' ','UNFORMATTED',IDUMMY,.FALSE.)
C
C
C     ***** first transpose CMO
C
      DO 720 ISYM=1,NSYM
         MORBI=MORB(ISYM)
      IF(MORBI.EQ.0) GO TO 720
         MFROI=MFRO(ISYM)
         MBASI=NBAS(ISYM)
         I1   = 1   + ICMO(ISYM) + MFROI*MBASI
         I2   = LW1 + JCMO(ISYM)
         CALL MTRSP(MBASI,MORBI,CMO(I1),MBASI,WORK(I2),MORBI)
  720 CONTINUE
C
C     ***** CALLING SEQUENCE FOR FIRST TRANSFORMATION STEP *****
C
      CALL TRACD(ITRLVL,LUHALF,LUDA,WORK(LW1),WORK(LW2),WORK(LW2),
     &           WORK(LW3),WORK(LW3),WORK(LW4),WORK(LW5),WORK(LW6),
     &           ITSMO,ITSAO,ITSX,IT)
      DFIN = SECOND()
      DTIM = DFIN - DST
      DST  = DFIN
      IF (JTER.EQ.1) WRITE (LUPRI,47) DTIM
   47 FORMAT(' TIME IN TRACD',F10.2)
C
      IF (DELDA) THEN
         CALL GPCLOSE(LUDA,'DELETE')
         LSRTAO = 0
      ELSE
         CALL GPCLOSE(LUDA,'KEEP')
      END IF
C
      IF (LPQCD.EQ.0) THEN
         WRITE (LUPRI,'(/A)') ' --- NO NON-ZERO MO INTEGRALS ---'
         GO TO 7000
C        7000: CALL TRAAB TO OPEN LUINTM AND WRITE INFORMATION
C              ABOUT NO INTEGRALS.
      END IF
C
C     ***** SET UP DYNAMIC STORAGE FOR TRAAB *****
C           870522: is now the same as for TRACD /hjaaj
C
C
C     ***** CALLING SEQUENCE FOR SECOND SORT *****
C
      CALL SORTB(ITRLVL,LUHALF,LUDA2,WORK,WORK,WORK(KWORK1),
     1           WORK(KWORK1),LWORK1,LTRI,IT)
      CALL GPCLOSE(LUHALF,'DELETE')
      DFIN = SECOND()
      DTIM = DFIN - DST
      IF (JTER.EQ.1) WRITE (LUPRI,48) DTIM
   48 FORMAT(' TIME IN SORTB',F10.2)
C
C     ***** CALL SEQUENCE FOR SECOND TRANSFORMATION STEP *****
C
 7000 CONTINUE
C
C     ***** first transpose CMO
C
      DO 820 ISYM=1,NSYM
         MORBI=MORB(ISYM)
      IF(MORBI.EQ.0) GO TO 820
         MFROI=MFRO(ISYM)
         MBASI=NBAS(ISYM)
         I1   = 1   + ICMO(ISYM) + MFROI*MBASI
         I2   = LW1 + JCMO(ISYM)
         CALL MTRSP(MBASI,MORBI,CMO(I1),MBASI,WORK(I2),MORBI)
  820 CONTINUE
C
C ***** save information about AO DA file (LUDA) on LUINTM
C
      CALL REWSPL(LUINTM)
      READ  (LUINTM)
C 831014-hjaaj: maybe here also save LUINTM buffer length
C 860731-hjaaj: done
      WRITE (LUINTM) LBINTM,ITRLVL,NSYM,NORB,NBAS
      WRITE (LUINTM) (CMO(I),I=1,NCMOT)
      WRITE (LUINTM) MISH(1:8), MASH(1:8)
C
      CALL TRAAB(ITRLVL,LUDA2,WORK(LW1),WORK(LW2),WORK(LW2),
     &           WORK(LW3),WORK(LW3),WORK(LW4),WORK(LW5),WORK(LW6),
     &           ITSMO,ITSAO,ITSW,ITSX,IT)
      IF (LPQCD.GT.0) CALL GPCLOSE(LUDA2,'DELETE')
C
      DFIN = SECOND()
      DTIM = DFIN - DSTTRA
      IF (IPRTRA .GT. 0) WRITE (LUPRI,49) DTIM
   49 FORMAT(' TOTAL TIME IN TRACD,SORTB,TRAAB',F10.2)
C
C set JTER to 2 to suppress printing on next call of TRACTL_1/hj-840710
      IF (IPRTRA .LT. 10) JTER = 2
C
C
 9999 CALL FLSHFO(LUPRI)
C
C Sort mo integrals to Dirac form, if requested
C 971201-hjaaj: after 9999 so we are sure Dirac form is present
C    if needed and LUINTM is OK
C
      IF (USEDRC .AND. ITRLVL .GT. 0) THEN
         IF (ITRLVL .GE. 4) THEN
            LVLDRC = 1
C           ... all occ-occ distributions
         ELSE
            LVLDRC = 0
C           ... all active-active distributions
         END IF
         CALL DRCCTL(LVLDRC,CMO,WORK,LWORK)
         CALL FLSHFO(LUPRI)
      END IF
      CALL QEXIT('TRACTL_1')
      RETURN
C
C
C     End of TRACTL_1
C
      END
C  /* Deck bdtra */
      BLOCK DATA BDTRA
C
C     Define MUL in /INTTRA/
C
#include "inttra.h"
C
C
C     MULTIPLICATION TABLE FOR SYMMETRIES
C
      DATA MUL/1,2,3,4,5,6,7,8,
     *         2,1,4,3,6,5,8,7,
     *         3,4,1,2,7,8,5,6,
     *         4,3,2,1,8,7,6,5,
     *         5,6,7,8,1,2,3,4,
     *         6,5,8,7,2,1,4,3,
     *         7,8,5,6,3,4,1,2,
     *         8,7,6,5,4,3,2,1/
C
      END
C  /* Deck nxth2m */
      SUBROUTINE NXTH2M(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C
C  Written by Hans Joergen Aa. Jensen October 1989
C  This version is interface routine for old integral transformation.
C
C NOTE: The space allocated in WRK must not be touched outside
C       until all desired distributions have been read.
C
C Purpose:
C    Read next Mulliken two-electron integral distribution (**|cd)
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C    (if needtp(itypcd) .gt. 0 all distributions of that type needed;
C     if needtp(itypcd) .lt. 0 at least one distribution needed;
C     if needtp(itypcd) .eq. 0 no distributions of that type needed).
C
C Usage:
C    Set IDIST = 0 before first call of NXTH2M.
C    DO NOT CHANGE IDIST or WRK(KFREE1:KFREE2-1) in calling routine
C    until last distribution has been read (signalled by IDIST .eq. -1)
C    Prototype code:
C     IDIST = 0
C     define NEEDTP(1:6)
C 100 CALL NXTH2M(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C     IF (IDIST .GT. 0) THEN
C        KW1 = KFREE
C        LW1 = LFREE
C        use (**|cd) distribution in H2CD as desired
C        WRK(KW1:KW1-1+LW1) may be used
C        GO TO 100
C     END IF
C
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
      DIMENSION H2CD(*),NEEDTP(6),WRK(LFREE)
C
C Used from common blocks:
C   INFTAP : LUINTM,LBINTM
C
#include <inftap.h>
#include <inftra.h>
C   
C
      SAVE KBUF, KIBUF, KNEXT
      DATA KNEXT /-1/
C
      IF (NEWTRA) THEN
         CALL N_NXTH2M(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
         RETURN
      END IF
      CALL QENTER('NXTH2M')
C
C     If first read (IDIST .eq. 0) then
C        get buffer length and which distributions have been saved
C        (lvlmol=0 for CI; lvlmol.ge.4 for inact-inact and inact-act;
C         lvlmol=10 for full integral transformation)
C     and allocate space for buffers.
C
      IF (IDIST .EQ. 0) THEN
         CALL REWSPL(LUINTM)
         READ (LUINTM,ERR=998,END=999)
         READ (LUINTM,ERR=998,END=999) LBINTM, LVLMOL
         LVLMIN = 0
         IF (NEEDTP(1).NE.0 .OR. NEEDTP(2).NE.0. .OR.
     &       NEEDTP(4).NE.0 .OR. NEEDTP(5).NE.0) LVLMIN = 2
         IF (NEEDTP(1).GT.0 .OR. NEEDTP(4).GT.0) LVLMIN = 4
C        level 4 is needed if all inac-inac or all inac-sec
C        distributions must be available (needtp .gt. 0).
         IF (NEEDTP(6).NE.0) LVLMIN = 10
         IF (LVLMOL .LT. LVLMIN) THEN
            WRITE (LUPRI,*)
     &         'NXTH2M error, LVLMOL on LUINTM too small'
            WRITE (LUPRI,*) 'LVLMOL =',LVLMOL
            WRITE (LUPRI,*) 'need   :',LVLMIN
            WRITE (LUPRI,*) 'NEEDTP :',(NEEDTP(J),J=1,6)
            CALL QTRACE(LUPRI)
            CALL QUIT('NXTH2M error, LVLMOL too small')
         END IF
C
C        We need to keep BUF and IBUF consecutively in memory. 
C        IBUF(LBINTM + 1) = LENGTH of record
C        K.Ruud, April 2000
C
         LBINTT = LBINTM*(IRAT + 1) + 2
         CALL MEMGET('INTE',KBUF ,LBINTT,WRK,KFREE,LFREE)
         KIBUF = KBUF + LBINTM
         KNEXT = KFREE
      ELSE
C        ... check that BUF,IBUF has not been destroyed by calling
C            routine.
         IF (KNEXT.EQ. -1  ) THEN
            WRITE (LUPRI,*)
     &         'NXTH2M error, IDIST must be zero in first call'
            WRITE (LUPRI,*) 'IDIST =',IDIST
            CALL QTRACE(LUPRI)
            CALL QUIT('NXTH2M error, IDIST must be zero in first call')
         END IF
         IF (KFREE.LT.KNEXT) THEN
            WRITE (LUPRI,*)
     &         'NXTH2M error, KFREE lower than buffer allocation'
            WRITE (LUPRI,*) 'KFREE ',KFREE
            WRITE (LUPRI,*) 'KBUF  ',KBUF
            WRITE (LUPRI,*) 'KIBUF ',KIBUF
            WRITE (LUPRI,*) 'KNEXT ',KNEXT,
     &         ' ( next avail. address after buffers)'
         END IF
         CALL MEMCHK('IDIST .ne. 0 MEMCHK in NXTH2M',WRK,KBUF)
         IF (KFREE.LT.KNEXT) THEN
            CALL QTRACE(LUPRI)
            CALL QUIT('NXTH2M error:KFREE lower than buffer allocation')
         END IF
      END IF
C
      CALL NX2H2M(IC,ID,H2CD,NEEDTP,WRK(KBUF),WRK(KIBUF),IDIST)
C
C     If no more buffers (IDIST .lt. 0) then release buffer space
C
      IF (IDIST .LT. 0) THEN
         CALL MEMREL('Releasing buffer space in NXTH2M',WRK,KBUF,
     &               KBUF,KFREE,LFREE)
         KNEXT = -1
      END IF
      CALL QEXIT('NXTH2M')
      RETURN
  998 CONTINUE
         WRITE (LUPRI,*) 'ERROR reading MOTWOINT in NXTH2M'
         CALL QUIT('ERROR reading MOTWOINT in NXTH2M')
  999 CONTINUE
         WRITE (LUPRI,*) 'END OF FILE reading MOTWOINT in NXTH2M'
         CALL QUIT('END OF FILE reading MOTWOINT in NXTH2M')
      END
C  /* Deck nx2h2m */
      SUBROUTINE NX2H2M(IC,ID,H2CD,NEEDTP,BUF,IBUF,IDIST)
C
C  Written by Hans Joergen Aa. Jensen October 1989
C
C Purpose:
C
C    Read next Mulliken two-electron integral distribution (**|cd)
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C
C Input:
C       NEEDTP(6); non-zero for needed (cd) distribution types
C                  if .gt. 0 all distributions of that type needed
C       IDIST; .eq. 0 first read
C              .gt. 1 intermediate read
C              .lt. 0 end-of-file has been reached previously
C Output:
C       H2CD(NORBT,NORBT); H2CD(a,b) = (ab|cd)
C       IC,ID; value of c and d
C       IDIST; .gt. 0 when next distribution IC,ID available in H2CD
C              = -1 when no more distributions
C Scratch:
C       BUF, IBUF
C
C ****************************************************************
C
#include <implicit.h>
#include <iratdef.h>
#include <inftap.h>
      DIMENSION H2CD(NORBT,NORBT),BUF(LBINTM), IBUF(LBINTM)
#include <priunit.h>
      DIMENSION NEEDTP(6)
C
C
C Used from common blocks:
C   INFORB : NORBT
C   INFIND : JTACT,JTSEC,ISW(*),...
C   INFTAP : LUINTM,LBINTM
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
#include <orbtypdef.h>
C
      SAVE INDCDN, LENGTH, JBUF
C
#include <ibtdef.h>
C
C ****************************************************************
C
C     If first read (IDIST .EQ. 0)
C     then setup for reading MO integrals ...
C
      IF (IDIST .EQ. 0) THEN
         JBUF = 2521*IRAT + 256 + NCMOT
         CALL REWSPL(LUINTM)
         CALL MOLLAB('MOLTWOEL',LUINTM,LUPRI)
 100     CONTINUE
#if defined (VAR_SPLITFILES)
         CALL READSI(LUINTM,LBINTM*IRAT,BUF,JBUF)
         CALL READSI(LUINTM,LBINTM + 2,IBUF,JBUF)
         LENGTH = IBUF(LBINTM + 1)
#else
         READ (LUINTM) BUF, IBUF, LENGTH
#endif
         IF (LENGTH.EQ.0) GOTO 100
         IF (LENGTH.EQ.-1) THEN
            WRITE (LUPRI,'(/A)')
     &         ' **** NXTH2M-ERROR *** no MO integrals on LUINTM'
            CALL QTRACE(LUPRI)
            CALL QUIT('*** ERROR ***(NXTH2M) no mo integrals on LUINTM')
         END IF
         INDCDN = IBTAND(IBTSHR(IBUF(1),16),IBT16)
      END IF
C
C *** Initialize H2CD
C
      CALL DZERO(H2CD,N2ORBX)
C
C *** Read next distribution which is needed according to NEEDTP(6)
C     into H2CD
C
C  ITYPCD values:  1=i*i :  2=t*i : 3=t*t : 4=a*i : 5=a*t : 6=a*a
C                  0 for not wanted type.
C
C  The CD distributions are stored by the present transformation
C  program with IC.ge.ID
C
  200 CONTINUE
      IF (INDCDN .EQ. -1) THEN
C        Last distribution has been read
         IDIST = -1
         GO TO 9999
      END IF
      INDCD = INDCDN
      IDIST = IDIST + 1
      IC = IBTAND(IBTSHR(INDCD,8),IBT08)
      ID = IBTAND(       INDCD,   IBT08)
      ITYPC  = IOBTYP(IC)
      ITYPD  = IOBTYP(ID)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF ( NEEDTP(ITYPCD) .EQ. 0 ) THEN
         ITYPCD = 0
      ELSE
         GO TO 400
C        ... first buffer is already in BUF,IBUF
      END IF
C
 300  CONTINUE
C
#if defined (VAR_SPLITFILES)
      CALL READSI(LUINTM,LBINTM*IRAT,BUF,JBUF)
      CALL READSI(LUINTM,LBINTM+2,IBUF,JBUF)
      LENGTH = IBUF(LBINTM + 1)
#else 
      READ (LUINTM) BUF, IBUF, LENGTH
#endif
      IF (LENGTH.EQ.0) GO TO 300
      IF (LENGTH.EQ.-1) THEN
         INDCDN = -1
      ELSE
         INDCDN = IBTAND(IBTSHR(IBUF(1),16),IBT16)
      END IF
      IF (INDCDN.NE.INDCD) THEN
C        ... finished reading the INDCD distribution
C            if not needed (ITYPCD .eq. 0) go find type of
C            new INDCDN distribution
         IF (ITYPCD.EQ.0) THEN
            GO TO 200
         ELSE
            GO TO 9999
         END IF
      END IF
      IF (ITYPCD.EQ.0) GO TO 300
C
  400 DO 450 I = 1,LENGTH
         IA = IBTAND(IBTSHR(IBUF(I),8),IBT08)
         IB = IBTAND(       IBUF(I),   IBT08)
         H2CD(IA,IB) = BUF(I)
         H2CD(IB,IA) = BUF(I)
  450 CONTINUE
      GO TO 300
C
C*******************************************************************
C
C End of subroutine NXTH2M
C
 9999 CONTINUE
      RETURN
      END
C  /* Deck nxth2d */
      SUBROUTINE NXTH2D(IC,ID,H2CD,NEEDTP,LUINTD,WRK,KFREE,LFREE,IDIST)
C
C  Written by Hans Joergen Aa. Jensen October 1989
C  This version is interface routine for old integral transformation.
C
C NOTE: The space allocated in WRK must not be touched outside
C       until all desired distributions have been read.
C
C Purpose:
C    Read next Dirac two-electron integral distribution <**|cd>
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C    (if needtp(itypcd) .gt. 0 all distributions of that type needed;
C     if needtp(itypcd) .lt. 0 at least one distribution needed;
C     if needtp(itypcd) .eq. 0 no distributions of that type needed).
C
C Usage:
C    Set IDIST = 0 before first call of NXTH2D.
C    DO NOT CHANGE IDIST or WRK(KFREE1:KFREE2-1) in calling routine
C    until last distribution has been read (signalled by IDIST .eq. -1)
C    Prototype code:
C     IDIST = 0
C     define NEEDTP(1:6)
C 100 CALL NXTH2D(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C     IF (IDIST .GT. 0) THEN
C        KW1 = KFREE
C        LW1 = LFREE
C        use (**|cd) distribution in H2CD as desired
C        WRK(KW1:KW1-1+LW1) may be used
C        GO TO 100
C     END IF
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      DIMENSION H2CD(*),NEEDTP(6),WRK(LFREE)
C
C Used from common blocks:
C   INFTAP : LUINTD,LBINTD
C
#include <inftap.h>
#include <inftra.h>
C
      SAVE KBUF, KIBUF, KNEXT
      DATA KNEXT /-1/
C
      IF (NEWTRA) THEN
         CALL N_NXTH2D(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
         RETURN
      END IF
      CALL QENTER('NXTH2D')
C
C     If first read (IDIST .eq. 0) then
C        get buffer length and which distributions have been saved
C        (lvldrc=0: act-act; =1: occ-occ; else orb-orb)
C     and allocate space for buffers.
C
      IF (IDIST .EQ. 0) THEN
         REWIND LUINTD
         CALL MOLLAB('DRCINFO ',LUINTD,LUPRI)
         READ (LUINTD) LBINTD, LVLDRC
         LVLMIN = 0
         IF (NEEDTP(1).NE.0 .OR. NEEDTP(2).NE.0) LVLMIN = 1
         IF (NEEDTP(4).NE.0 .OR. NEEDTP(5).NE.0 .OR.
     &       NEEDTP(6).NE.0) LVLMIN = 2
         IF (LVLDRC .LT. LVLMIN) THEN
            WRITE (LUPRI,*)
     &         'NXTH2D error, LVLDRC on LUINTD too small'
            WRITE (LUPRI,*) 'LVLDRC =',LVLDRC
            WRITE (LUPRI,*) 'need   :',LVLMIN
            CALL QTRACE(LUPRI)
            CALL QUIT('NXTH2D error, LVLDRC too small')
         END IF
         CALL MEMGET('REAL',KBUF ,LBINTD,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KIBUF,LBINTD,WRK,KFREE,LFREE)
         KNEXT = KFREE
      ELSE
C        ... check that BUF,IBUF has not been destroyed by calling
C            routine.
         IF (KNEXT.EQ. -1  ) THEN
            WRITE (LUPRI,*)
     &         'NXTH2D error, IDIST must be zero in first call'
            WRITE (LUPRI,*) 'IDIST =',IDIST
            CALL QTRACE(LUPRI)
            CALL QUIT('NXTH2D error, IDIST must be zero in first call')
         END IF
         IF (KFREE.LT.KNEXT) THEN
            WRITE (LUPRI,*)
     &         'NXTH2D error, KFREE lower than buffer allocation'
            WRITE (LUPRI,*) 'KFREE ',KFREE
            WRITE (LUPRI,*) 'KBUF  ',KBUF
            WRITE (LUPRI,*) 'KIBUF ',KIBUF
            WRITE (LUPRI,*) 'KNEXT ',KNEXT,
     &         ' ( next avail. address after buffers)'
         END IF
         CALL MEMCHK('IDIST .ne. 0 MEMCHK in NXTH2D',WRK,KBUF)
         IF (KFREE.LT.KNEXT) THEN
            CALL QTRACE(LUPRI)
            CALL QUIT('NXTH2D error:KFREE lower than buffer allocation')
         END IF
      END IF
      CALL NX2H2D(LUINTD,IC,ID,H2CD,NEEDTP,WRK(KBUF),WRK(KIBUF),IDIST)
C
C     If no more buffers (IDIST .lt. 0) then release buffer space
C
      IF (IDIST .LT. 0) THEN
         CALL MEMREL('Releasing buffer space in NXTH2D',WRK,KBUF,
     &               KBUF,KFREE,LFREE)
         KNEXT = -1
      END IF
      CALL QEXIT('NXTH2D')
      RETURN
      END
C  /* Deck nx2h2d */
      SUBROUTINE NX2H2D(LUINTD,IC,ID,H2CD,NEEDTP,BUF,IBUF,IDIST)
C
C  Written by Hans Joergen Aa. Jensen October 1989
C
C Purpose:
C
C    Read next Dirac two-electron integral distribution <**|cd>
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C
C Input:
C       NEEDTP(6); non-zero for needed (cd) distribution types
C                  if .gt. 0 all distributions of that type needed
C       IDIST; .eq. 0 first read
C              .gt. 1 intermediate read
C              .lt. 0 end-of-file has been reached previously
C Output:
C       H2CD(NORBT,NORBT); H2CD(a,b) = <ab|cd>
C       IC,ID; value of c and d
C       IDIST; .gt. 0 when next distribution IC,ID available in H2CD
C              = -1 when no more distributions
C Scratch:
C       BUF, IBUF
C
C ****************************************************************
C
#include <implicit.h>
#include <inftap.h>
      DIMENSION H2CD(NORBT,NORBT),BUF(LBINTD),IBUF(LBINTD)
      DIMENSION NEEDTP(6)
C
C
C Used from common blocks:
C   INFORB : NORBT
C   INFIND : JTACT,JTSEC,ISW(*),...
C   INFTAP : LUINTD,LBINTD
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
#include <orbtypdef.h>
C
      SAVE INDCDN, LENGTH
C
#include <ibtdef.h>
C
C ****************************************************************
C
C     If first read (IDIST .EQ. 0)
C     then setup for reading MO integrals ...
C
      IF (IDIST .EQ. 0) THEN
         REWIND LUINTD
         CALL MOLLAB('DRCTWOEL',LUINTD,LUPRI)
  100    READ (LUINTD) BUF,IBUF,LENGTH
         IF (LENGTH.EQ.0) GOTO 100
         IF (LENGTH.EQ.-1) THEN
            WRITE (LUPRI,'(/A)')
     &         ' **** NXTH2D-ERROR *** no MO integrals on LUINTD'
            CALL QTRACE(LUPRI)
            CALL QUIT('*** ERROR ***(NXTH2D) no mo integrals on LUINTD')
         END IF
         INDCDN = IBTAND(IBTSHR(IBUF(1),16),IBT16)
      END IF
C
C *** Initialize H2CD
C
      CALL DZERO(H2CD,N2ORBX)
C
C *** Read next distribution which is needed according to NEEDTP(6)
C     into H2CD
C
C  ITYPCD values:  1=i*i :  2=t*i : 3=t*t : 4=a*i : 5=a*t : 6=a*a
C                  0 for not wanted type.
C
C  The CD distributions are stored by the present transformation
C  program with IC.ge.ID
C
  200 CONTINUE
      IF (INDCDN .EQ. -1) THEN
C        Last distribution has been read
         IDIST = -1
         GO TO 9999
      END IF
      INDCD = INDCDN
      IDIST = IDIST + 1
      IC = IBTAND(IBTSHR(INDCD,8),IBT08)
      ID = IBTAND(       INDCD,   IBT08)
      ITYPC  = IOBTYP(IC)
      ITYPD  = IOBTYP(ID)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF ( NEEDTP(ITYPCD) .EQ. 0 ) THEN
         ITYPCD = 0
      ELSE
         GO TO 400
C        ... first buffer is already in BUF,IBUF
      END IF
C
  300 READ (LUINTD) BUF,IBUF,LENGTH
      IF (LENGTH.EQ.0) GO TO 300
      IF (LENGTH.EQ.-1) THEN
         INDCDN = -1
      ELSE
         INDCDN = IBTAND(IBTSHR(IBUF(1),16),IBT16)
      END IF
      IF (INDCDN.NE.INDCD) THEN
C        ... finished reading the INDCD distribution
C            if not needed (ITYPCD .eq. 0) go find type of
C            new INDCDN distribution
         IF (ITYPCD.EQ.0) THEN
            GO TO 200
         ELSE
            GO TO 9999
         END IF
      END IF
      IF (ITYPCD.EQ.0) GO TO 300
C
  400 DO 450 I = 1,LENGTH
         IA = IBTAND(IBTSHR(IBUF(I),8),IBT08)
         IB = IBTAND(       IBUF(I),   IBT08)
         H2CD(IA,IB) = BUF(I)
  450 CONTINUE
      GO TO 300
C
C*******************************************************************
C
C End of subroutine NXTH2D
C
 9999 CONTINUE
      RETURN
      END
