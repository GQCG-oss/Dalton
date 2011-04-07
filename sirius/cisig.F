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
===========================================================================
Revision 1.2  2000/05/05 11:19:33  hjj
bug fix: now always calculate ci diag. in CIDAG, previously it was only
done if LUIT2 .ge. 0 (must be a leftover, the code LUIT2 .lt. 0 was not
used anywhere, rather UPDGRD expected diagonal to be calculated with
LUIT2 .lt. 0 !!!).
===========================================================================
#endif
C  /* Deck tpblm2 */
      SUBROUTINE TPBLM2(A,AT,LBLR,NBLR,LBLC,NBLC,IFLAG,IORD)
C
C A BLOCKED MATRIX A IS GIVEN .
C
C THE BLOCK STRUCTURE CAN BE OF THE FOLLOWING TYPES
C IORD = 1 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER COLUMNS OF THIS BLOCK
C             LOOP OVER ROWS OF THIS BLOCK
C
C IORD = 2 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER ROWS OF THIS BLOCK
C             LOOP OVER COLUMNS OF THIS BLOCK
C
C     FOR IORD = 2 ARE THE INDIVIDUAL BLOCKS THUS TRANSPOSED
C
C THE COMBINATION OF TWO BLOCKS IABL AND IBBL ARE ALLOWED
C IF IFLAG(IABL,IBBL) = 1
C
C TRANSPOSE THE INDIVIDUAL BLOCKS OF THIS MATRIX TO GIVE AT
C THE ORDER OF THE BLOCKS ARE NOT CHANGED
C
C JEPPE OLSEN , NOVEMBER 1988
C
#include <implicit.h>
      DIMENSION A(*),AT(*)
      DIMENSION LBLR(NBLR),LBLC(NBLC)
      DIMENSION IFLAG(NBLR,NBLC)
C
      IOFF = 1
      DO 200 IBLR = 1, NBLR
        DO 100 IBLC = 1, NBLC
          IF(IFLAG(IBLR,IBLC) .EQ. 1 ) THEN
            LR = LBLR(IBLR)
            LC = LBLC(IBLC)
            IF( IORD .EQ. 1 ) THEN
              CALL TRPMAT(A(IOFF),LR,LC,AT(IOFF))
            ELSE IF( IORD .EQ. 2 ) THEN
              CALL TRPMAT(A(IOFF),LC,LR,AT(IOFF))
            END IF
            IOFF = IOFF + LR * LC
          END IF
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck prblm2 */
      SUBROUTINE PRBLM2(A,LBLR,NBLR,LBLC,NBLC,IFLAG,IORD)
C
C A BLOCKED MATRIX A IS GIVEN .
C
C THE BLOCK STRUCTURE CAN BE OF THE FOLLOWING TYPES
C IORD = 1 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER COLUMNS OF THIS BLOCK
C             LOOP OVER ROWS OF THIS BLOCK
C
C IORD = 2 :
C     LOOP OVER BLOCK OF ROWS
C       LOOP OVER BLOCK OF COLUMNS ALLOWED FOR GIVEN ROW BLOCK
C           LOOP OVER ROWS OF THIS BLOCK
C             LOOP OVER COLUMNS OF THIS BLOCK
C
C     FOR IORD = 2 ARE THE INDIVIDUAL BLOCKS THUS TRANSPOSED
C
C THE COMBINATION OF TWO BLOCKS IABL AND IBBL ARE ALLOWED
C IF IFLAG(IABL,IBBL) = 1
C
C PRINT THIS MATRIX !
C
C JEPPE OLSEN , NOVEMBER 1988
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION A(*)
      DIMENSION LBLR(NBLR),LBLC(NBLC)
      DIMENSION IFLAG(NBLR,NBLC)
C
      IOFF = 1
      DO 200 IBLR = 1, NBLR
        DO 100 IBLC = 1, NBLC
          IF(IFLAG(IBLR,IBLC) .EQ. 1 ) THEN
            WRITE(LUPRI,*) ' BLOCK INDICES ',IBLR,IBLC
            LR = LBLR(IBLR)
            LC = LBLC(IBLC)
            IF(IORD.EQ.1) THEN
              CALL WRTMAT(A(IOFF),LR,LC,LR,LC,1)
            ELSE IF( IORD .EQ. 2 ) THEN
              CALL WRTMAT(A(IOFF),LC,LR,LC,LR,1)
            END IF
            IOFF = IOFF + LR * LC
          END IF
  100   CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck traci */
      SUBROUTINE TRACI(NCVEC,CVEC,LCVEC,UMO,XNDXCI,WRK,LFREE,IPRTCI)
C
C Driver for TRACI using determinant routines
C Version with CSF business
C
C  JULY 15 1988 Jeppe Olsen
C
C
C MOTECC-90: The algorithms used in this module, TRACI, are
C            described in Chapter 8 Section D.5 of MOTECC-90
C            "Counter Rotations of CI coefficients"
C
#include <implicit.h>
      DIMENSION CVEC(LCVEC,NCVEC), UMO(*), XNDXCI(*), WRK(LFREE)
C
C Used from common blocks:
C  PRIUNIT : IPRSTAT
C  INFINP : LSYM
C  INFORB : MULD2H(8,8), NASHT, N2ASHX
C  INFVAR : NCONF
C  DETBAS : KIASTR,KIBSTR, ...
C  STRNUM : NAEL,NBEL, ...
C  CIINFO : NDTASM
C  CSFBAS : Pointers for CSF information
C  CBREOR : SLREOR
C
#include <maxorb.h>
#include <mxpdim.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <detbas.h>
#include <strnum.h>
#include <csfbas.h>
#include <ciinfo.h>
#include <cbreor.h>
      LOGICAL CSFEXP
C
C
      CSFEXP = .NOT.FLAG(27)
      IF (CSFEXP) THEN
        NCDET = NDTASM(LSYM)
        IF (NCONF .NE. NCSASM(LSYM)) THEN
           WRITE (LUPRI,*) 'TRACI ERROR, NCONF .ne. NCSASM(LSYM):',
     *        NCONF, NCSASM(LSYM)
           CALL QUIT('TRACI ERROR, NCONF .ne. NCSASM(LSYM)')
        END IF
      ELSE
        NCDET = NCONF
        IF (NCONF .NE. NDTASM(LSYM)) THEN
           WRITE (LUPRI,*) 'TRACI ERROR, NCONF .ne. NDTASM(LSYM):',
     *        NCONF, NDTASM(LSYM)
           CALL QUIT('TRACI ERROR, NCONF .ne. NDTASM(LSYM)')
        END IF
      END IF
c Local memory
      KLSOT  = 1
      KLC2   = KLSOT + N2ASHX
      KLFREE = KLC2  + MAX(N2ASHX, NCDET )
      IF(CSFEXP) THEN
        KLDET1 = KLFREE
        KLDET2 = KLDET1 + NCDET
        KLFREE = KLDET2 + NCDET
      END IF
C..   note: C2 only needs max of nia*nib (see TRACI1)
      IF ( KLFREE-1 .GT. LFREE ) CALL ERRWRK('TRACI',KLFREE-1,LFREE)
c orbital tranformation matrix
      IF (.NOT. SLREOR) THEN
         CALL DCOPY(N2ASHX,UMO,1,WRK(KLSOT),1)
      ELSE
C        reorder from Sirius order to Lunar order
         CALL REORMT(UMO,WRK(KLSOT),NASHT,NASHT,XNDXCI(KLTSOB),
     &               XNDXCI(KLTSOB) )
      END IF
      CALL SOTMAT(NASHT,WRK(KLSOT),IFAIL)
      IF ( IFAIL .NE. 0 ) THEN
         WRITE(LUPRI,'(///,A,I5,/)')
     &        ' TRACI : ERROR IN SOTMAT, IFAIL = ',IFAIL
         CALL QUIT(' TRACI : ERROR IN SOTMAT ')
      END IF
      IF (IPRTCI .GT. 10 ) THEN
         WRITE(LUPRI,'(/A)') ' SINGLE ORBITAL TRANSFORMATION MATRIX '
         CALL WRTMAT(WRK(KLSOT),NASHT,NASHT,NASHT,NASHT,0)
      END IF
C..TRANSPOSE SOT MATRIX IN ACCORDANCE WITH LUNAR DESCRIPTION
      CALL DGETRN(WRK(KLSOT),NASHT,NASHT)
c Store CI offsets in C arrays
      IPRXXX = MAX(0,IPRTCI - 10)
      CALL CIOFF(LSYM,1,XNDXCI,IPRXXX)
c     CALL CIOFF(IREFSM,ICORHC,XNDXCI,NTEST)
c
      DO 450 IVEC = 1, NCVEC
        IF(CSFEXP) THEN
c
c CSF to DET transformation, rotate in DET basis, DET to CSF transformation
c
          CALL COPVEC(CVEC(1,IVEC),WRK(KLDET1),NCONF)
          CALL CSDTVC(WRK(KLDET1),WRK(KLDET2),1,XNDXCI(KDTOC),
     &                XNDXCI(KICTS(1)),LSYM,0,IPRXXX)
C         CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,IREFSM,ICOPY,NTEST)
          CALL TRACI2(WRK(KLSOT),NASHT,LSYM,NAEL,XNDXCI(KIASTR),
     &                XNDXCI(KTAIJ),XNDXCI(KTATO),XNDXCI(KTASYM),
     &                XNDXCI(KSTASA),XNDXCI(KSTBAA),NBEL,
     &                XNDXCI(KIBSTR),XNDXCI(KTBIJ),
     &                XNDXCI(KTBTO),XNDXCI(KTBSYM),XNDXCI(KSTASB),
     &                XNDXCI(KSTBAB),NASTR,NBSTR,XNDXCI(KIANNI),
     &                WRK(KLDET2),XNDXCI(KCOFF),WRK(KLC2),
     &                NOCTPA,XNDXCI(KNSSOA),XNDXCI(KISSOA),
     &                XNDXCI(KTPFSA),NOCTPB,XNDXCI(KNSSOB),
     &                XNDXCI(KISSOB),
     &                XNDXCI(KTPFSB),XNDXCI(KIOCOC),XNDXCI(KICOOS),
     &                XNDXCI(KCDTAS),MULD2H,MAXSYM,IPRTCI)
          CALL CSDTVC(WRK(KLDET1),WRK(KLDET2),2,XNDXCI(KDTOC),
     &                XNDXCI(KICTS(1)),LSYM,0,IPRXXX)
C         CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,IREFSM,ICOPY,NTEST)
          CALL COPVEC(WRK(KLDET1),CVEC(1,IVEC),NCONF)
        ELSE
c
c Rotate in DET expansion
c
          CALL TRACI2(WRK(KLSOT),NASHT,LSYM,NAEL,XNDXCI(KIASTR),
     &                XNDXCI(KTAIJ),XNDXCI(KTATO),XNDXCI(KTASYM),
     &                XNDXCI(KSTASA),XNDXCI(KSTBAA),NBEL,
     &                XNDXCI(KIBSTR),XNDXCI(KTBIJ),
     &                XNDXCI(KTBTO),XNDXCI(KTBSYM),XNDXCI(KSTASB),
     &                XNDXCI(KSTBAB),NASTR,NBSTR,XNDXCI(KIANNI),
     &                CVEC(1,IVEC),XNDXCI(KCOFF),WRK(KLC2),
     &                NOCTPA,XNDXCI(KNSSOA),XNDXCI(KISSOA),
     &                XNDXCI(KTPFSA),NOCTPB,XNDXCI(KNSSOB),
     &                XNDXCI(KISSOB),
     &                XNDXCI(KTPFSB),XNDXCI(KIOCOC),XNDXCI(KICOOS),
     &                XNDXCI(KCDTAS),MULD2H,MAXSYM,IPRTCI)
        END IF
  450 CONTINUE
C
      RETURN
      END
C  /* Deck traci2 */
      SUBROUTINE TRACI2(T,NORB,ICSM,
     &                  NAEL,IASTR,TAIJ,TATO,TASYM,NSTASA,ISTBAA,
     &                  NBEL,IBSTR,TBIJ,TBTO,TBSYM,NSTASB,ISTBAB,
     &                  NASTR,NBSTR, IANNI,C,ICOFF,C2,
     &                  NOCTPA,NSTAOS,ISTAOS,ITPAST,
     &                  NOCTPB,NSTBOS,ISTBOS,ITPBST,
     &                  IOCOC,ICOOS,NCDTAS,SYMPRO,MAXSYM,NTEST )
C
C  COUNTER ROTATE CI - COEFFICIENTS CORRESPONDING
C  TO ORBITAL ROTATIONS DEFINED IN T
C  INPUT CI VECTOR IS C AND OUTPUT CI VECTOR OVERWRITERS C
C  C2 IS SCRATCH ABLE TO HOLD LARGEST SYMMETRY BLOCK
C  SYMMETRY OF CIN AND COUT ARE ASSUMED TO BE IDENTICAL AND EQUAL TO
C  ICSM
C
C d2h RAS Version , Jeppe Olsen November 1988
C ( Debugged ! Error corrected 890809-hjaaj)
c
c New format of single excitations, October 1990
c
#include <implicit.h>
#include <priunit.h>
      INTEGER SYMPRO(8,8)
      INTEGER TAIJ(*),TATO(*)
      INTEGER TBIJ(*),TBTO(*)
      INTEGER TASYM(MAXSYM+1,*),TBSYM(MAXSYM+1,*)
      INTEGER NSTASA(MAXSYM),NSTASB(MAXSYM)
      INTEGER ISTBAA(MAXSYM),ISTBAB(MAXSYM)
      INTEGER ICOFF(MAXSYM),IANNI(NORB ** 2 )
      INTEGER IASTR(NAEL,NASTR),IBSTR(NBEL,NBSTR)
      INTEGER NSTAOS(NOCTPA,*),ISTAOS(NOCTPA,*)
      INTEGER NSTBOS(NOCTPB,*),ISTBOS(NOCTPB,*)
      INTEGER IOCOC(NOCTPA,NOCTPB),ICOOS(NOCTPB,NOCTPA,MAXSYM)
      INTEGER NCDTAS(*), ITPAST(*),ITPBST(*)
C
      DIMENSION T(*), C(*),C2(*)
C LENGTH OF C2 MUST BE AT LEAST
C
      IF( NTEST .GE. 20 ) WRITE(LUPRI,*) ' >>>> INFORMATION FROM TRACI2'
C
C
      DO 2000 IASM = 1, MAXSYM
       IF( NTEST .GE. 25 ) WRITE(LUPRI,*) ' INFO FOR SYMMETRY... ',IASM
       KDET   = NCDTAS(IASM)
       IBSM   = SYMPRO(IASM,ICSM)
       IOFF11 = ICOOS(1,1,IASM)
C.. TRANSFORM THE BETA STRINGS
       IF(NTEST.GE.25) WRITE(LUPRI,*) ' TRANSFORMATION OF BETA STRINGS'
C..    TRANSFORM ORBITAL IORB
       DO 900 IORB = 1, NORB
        CALL SETVEC(C2,0.0D0,KDET)
        IF(NTEST .GE. 25)WRITE(LUPRI,*)'  *** ORBITAL TO BE ROTATED ',
     &                                 IORB
        DO  890 IATP =  1, NOCTPA
        DO  880 IBTP =  1, NOCTPB
         IF( IOCOC(IATP,IBTP) .NE. 1 ) GOTO  880
         IOFF =   ICOOS (IBTP,IATP,IASM)
C
         NIA = NSTAOS(IATP,IASM)
         NIB = NSTBOS(IBTP,IBSM)
         IF ( NTEST .GE. 25 ) THEN
           WRITE(LUPRI,*) ' INITIAL CI BLOCK '
           CALL WRTMAT(C(IOFF),NIA,NIB,NIA,NIB,0)
         END IF
C
         IBSTRT = ISTBOS(IBTP,IBSM)
         IASTRT = ISTAOS(IATP,IASM)
           DO 800 IB = IBSTRT, IBSTRT+NIB-1
             IOFFI = IOFF + (IB-IBSTRT)*NIA
C IS ORBITAL IORB OCCUPIED IN IB ?
             IOCC = 0
             DO 750 IEL = 1, NBEL
               IF(IBSTR(IEL,IB) .EQ. IORB ) IOCC = 1
  750        CONTINUE
C
             IF ( IOCC .EQ. 0 ) THEN
               CALL VECSUM(C2(IOFFI-IOFF11+1),C2(IOFFI-IOFF11+1),
     &                     C(IOFFI),1.0D0,1.0D0,NIA)
             ELSE
               DO 700 IEX = TBSYM(1,IB),TBSYM(2,IB)-1
                 IJEX = TBIJ(IEX)
                 IF(IANNI(IJEX) .EQ. IORB ) THEN
                   JB = TBTO(IEX)
                   IF( JB .GT. NBSTR ) THEN
                      JB = JB - NBSTR
                      FACTOR = -T(IJEX)
                   ELSE
                      FACTOR =  T(IJEX)
                   END IF
C
                   JBTP = ITPBST(JB)
                   IF( IOCOC(IATP,JBTP) .NE. 1 ) GOTO 700
                   JBEFF = JB - ISTBOS(JBTP,IBSM)+1
                   IOFFJ = ICOOS(JBTP,IATP,IASM) - IOFF11
     &                   + (JBEFF-1)*NIA + 1
                   CALL VECSUM(C2(IOFFJ),C2(IOFFJ),
     &                         C(IOFFI),1.0D0,FACTOR,NIA)
                 END IF
  700          CONTINUE
             END IF
  800      CONTINUE
C
  880    CONTINUE
  890    CONTINUE
         CALL COPVEC(C2,C(IOFF11),KDET)
         IF( NTEST .GE. 50 ) THEN
           WRITE(LUPRI,*) ' ROTATED CI BLOCK '
           CALL PRBLM2(C(IOFF11),NSTAOS(1,IASM),NOCTPA,
     &          NSTBOS(1,IBSM),NOCTPB,IOCOC,1)
         END IF
  900    CONTINUE
         IF( NTEST .GE. 25 .AND. NTEST .LT. 50) THEN
           WRITE(LUPRI,*) ' ROTATED CI BLOCK '
           CALL PRBLM2(C(IOFF11),NSTAOS(1,IASM),NOCTPA,
     &          NSTBOS(1,IBSM),NOCTPB,IOCOC,1)
         END IF
C
C      TPBLM2(A,AT,LBLR,NBLR,LBLC,NBLC,IFLAG,IORD)
C      PRBLM2(A,LBLR,NBLR,LBLC,NBLC,IFLAG,IORD)
C.. BETA STRINGS HAVE NOW BEEN TRANSFORMED , TRANSFORM ALPHA QQ
C   STRINGS
C.. TRANSPOSE TO EASE INDEXING
        IF(NTEST .GE. 25)WRITE(LUPRI,*)' TRANSFORMATION OF ALPHA '//
     &                                 'STRINGS'
C        CALL TRPMAT(C2,NIA,NIB,C(IOFF))
         CALL TPBLM2(C2,C(IOFF11),NSTAOS(1,IASM),NOCTPA,
     &        NSTBOS(1,IBSM),NOCTPB,IOCOC,1)
         DO 1900 IORB = 1, NORB
          IF( NTEST .GE. 25 )WRITE(LUPRI,*)' ORBITAL TO BE ROTATED ',
     &                                       IORB
          CALL SETVEC(C2,0.0D0,KDET)
          DO 1890 IATP = 1, NOCTPA
          DO 1880 IBTP = 1, NOCTPB
           IF( IOCOC(IATP,IBTP) .NE. 1 ) GOTO 1880
           NIA = NSTAOS(IATP,IASM)
           NIB = NSTBOS(IBTP,IBSM)
           IASTRT = ISTAOS(IATP,IASM)
           IBSTRT = ISTBOS(IBTP,IBSM)
           IOFF = ICOOS(IBTP,IATP,IASM)
           IOFFR = IOFF - IOFF11 + 1
           DO 1800 IA = IASTRT, IASTRT+NIA-1
             IOFFI = IOFFR + (IA-IASTRT)*NIB
C IS ORBITAL IORB OCCUPIED IN IA ?
             IOCC = 0
             DO 1750 IEL = 1, NAEL
               IF(IASTR(IEL,IA) .EQ. IORB ) IOCC = 1
 1750        CONTINUE
C
             IF ( IOCC .EQ. 0 ) THEN
               CALL VECSUM(C2(IOFFI),C2(IOFFI),
     &                     C(IOFFI+IOFF11-1),1.0D0,1.0D0,NIB)
             ELSE
               DO 1700 IEX = TASYM(1,IA),TASYM(2,IA)-1
                 IJEX = TAIJ(IEX)
                 IF(IANNI(IJEX) .EQ. IORB ) THEN
                   JA = TATO(IEX)
                   IF( JA .GT. NASTR ) THEN
                      JA = JA - NASTR
                      FACTOR = -T(IJEX)
                   ELSE
                      FACTOR =  T(IJEX)
                   END IF
                   JATP = ITPAST(JA)
C
                   IF(IOCOC(JATP,IBTP) .NE. 1 ) GOTO 1700
C
                   JASTRT = ISTAOS(JATP,IASM)
                   IOFFJ = ICOOS(IBTP,JATP,IASM)-IOFF11
     &                   + (JA-JASTRT)*NIB + 1
                   CALL VECSUM(C2(IOFFJ),C2(IOFFJ),C(IOFFI+IOFF11-1),
     &                         1.0D0,FACTOR,NIB)
C
                 END IF
 1700          CONTINUE
             END IF
 1800      CONTINUE
 1880      CONTINUE
 1890      CONTINUE
           CALL COPVEC(C2,C(IOFF11),KDET)
 1900    CONTINUE
C        CALL TRPMAT(C2,NIB,NIA,C(IOFF) )
         CALL TPBLM2(C2,C(IOFF11),NSTAOS(1,IASM),NOCTPA,
     &        NSTBOS(1,IBSM),NOCTPB,IOCOC,2)
         IF( NTEST .GE. 20 ) THEN
           WRITE(LUPRI,*) ' ROTATED CI BLOCK FOR IASM ...', IASM
           CALL PRBLM2(C(IOFF11),NSTBOS(1,IBSM),NOCTPB,
     &          NSTAOS(1,IBSM),NOCTPB,IOCOC,1)
           END IF
 2000 CONTINUE
C
      RETURN
      END
C  /* Deck typstr */
      SUBROUTINE TYPSTR(STRING,NEL,NSTRIN,ITYPE,IAB,NTEST)
C
C  OCCUPATION TYPE OF STRINGS
C
#include <implicit.h>
#include <priunit.h>
      INTEGER STRING(NEL,NSTRIN),ITYPE(NSTRIN)
C
      DO 100 ISTRIN = 1, NSTRIN
        ITYPE(ISTRIN) = IOCTYP(STRING(1,ISTRIN),IAB,1)
  100 CONTINUE
C
      IF ( NTEST .GT. 8 ) THEN
         WRITE(LUPRI,*) ' TYPSTR: OCCUPATION TYPES OF STRINGS '
         CALL IWRTMA(ITYPE,1,NSTRIN,1,NSTRIN)
      END IF
C
      RETURN
      END
C  /* Deck cidiag */
      SUBROUTINE CIDIAG(ICSYM,NOH2,FCAC,H2AC,XNDXCI,DIAGC,WRK,LFREE)
C
C Written 18-jan-1988 J.O.
C ( after cidiag of hjaaj )
C NOH2 parameter 900717-hjaaj
C
C Part of routines for determinant based ci
C
C Purpose:
C
C  Calculate the diagonal of the CI matrix and,
C  if LUIT2 .gt. 0, save it on LUIT2.
C
C
#include <implicit.h>
      DIMENSION FCAC(*),H2AC(*),DIAGC(*)
      DIMENSION XNDXCI(*), WRK(LFREE)
      LOGICAL   NOH2, FNDLAB
#include <iratdef.h>
C
C Used from common blocks:
C   INFINP : FLAG(*)
C   INFORB : N2ASHX
C   INFVAR : NCONF
C   INFTAP : LUIT2
C   INFPRI : IPRDIA
C   CIINFO : NDTASM(*)
C   SPINFO : ?
C   CSFBAS : CSF core allocation for XNDXCI
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
#include <ciinfo.h>
#include <spinfo.h>
#include <csfbas.h>
C
C *** local variables:
      LOGICAL     CSFEXP
      CHARACTER*8 TABLE(4)
      DATA TABLE/'********','CIDIAG1 ','CIDIAG2 ',
     &           'EODATA  '/
C
      CSFEXP = (.NOT.FLAG(27) .AND. ICSYM .EQ. LSYM )
C     ... 980820-hjaaj: only CSF for ref.sym. LSYM,
C         determinants for other symmetries
C
C *** Construct diagonal of CI matrix
      KFIJ  = 1
      KFJI  = KFIJ  + N2ASHX
      KW1   = KFJI  + N2ASHX
      LW1   = LFREE - KW1
      IF (NOH2) THEN
         CALL DZERO(WRK(KFIJ),2*N2ASHX)
      ELSE
         CALL GETFIJ(WRK(KFIJ),WRK(KFJI),H2AC)
      END IF
      IF (.NOT.CSFEXP) THEN
         CALL DZERO (DIAGC,NCONF)
         CALL CIDIAD(DIAGC,FCAC,WRK(KFIJ),WRK(KFJI),XNDXCI,ICSYM,
     *               WRK(KW1),LW1)
      ELSE
C
C.. Change to CSF format
C
         KDETDG = KW1
         KW1    = KDETDG + NDTASM(ICSYM)
         LW1    = LFREE  - KW1
         CALL CIDIAD(WRK(KDETDG),FCAC,WRK(KFIJ),WRK(KFJI),XNDXCI,
     *               ICSYM,WRK(KW1),LW1)
         CALL CSDIAG(DIAGC,WRK(KDETDG),NCNFTP(1,ICSYM),NTYP,
     &               XNDXCI(KICTS(1)),NDTFTP,NCSFTP,0,IPRDIA)
C        CSDIAG(CSFDIA,DETDIA,NCNFTP,NTYP,
C    &          ICTSDT,NDTFTP,NCSFTP,IFLAG,NTEST)
      END IF
C
      IF (LUIT2 .GT. 0) THEN
C
C *** Write diagonal of CI matrix on LUIT2 (label: 'CIDIAG2 ' )
C
         NC4 = MAX(4,NCONF)
         REWIND LUIT2
         IF( FNDLAB(TABLE(4),LUIT2) ) THEN
            BACKSPACE LUIT2
         ELSE
            REWIND LUIT2
         END IF
         WRITE (LUIT2) TABLE(1),TABLE(1),TABLE(1),TABLE(3)
         CALL WRITT(LUIT2,NC4,DIAGC)
         WRITE (LUIT2) TABLE(1),TABLE(1),TABLE(1),TABLE(4)
      END IF
C
C *** End of subroutine CIDIAG
C
      RETURN
      END
C  /* Deck cidiad */
      SUBROUTINE CIDIAD(DIAG,FCAC,FIJ,FJI,XNDXCI,ICSYM,WRK,LFREE)
C
C ***  CONSTRUCT DIAGONAL PART OF CI MATRIX
C
C      SOME INPUT
C                 FCAC : inactive Fock matrix
C                        (modified one body integrals)
C                 FIJ  : (II/JJ) integrals
C                 FJI  : (IJ/IJ) integrals
C
#include <implicit.h>
      DIMENSION DIAG(*), FCAC(*), FIJ(*), FJI(*), WRK(LFREE)
      DIMENSION XNDXCI(*)
#include <priunit.h>
C
C
C Used from common blocks:
C   INFINP : 
C   INFORB : MULD2H(8,8),NASHT,N2ASHX
C   INFIND : NSM(NASHT)
C   INFPRI : IPRDIA
C
C   DETBAS : core allocation
C   STRNUM : NAEL,NBEL,NASTR,NBSTR,?
C   CIINFO : NDTASM(*),ICOMBI
C
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
#include <mxpdim.h>
#include <detbas.h>
#include <strnum.h>
#include <ciinfo.h>
C
      CALL GETTIM(T0,W0)
      ECORE  = 0.0D0
c Store CI offsets in C arrays
      CALL CIOFF(ICSYM,1,XNDXCI,IPRDIA)
C     CALL CIOFF(IREFSM,ICORHC,XNDXCI,NTEST)
C
C
      KFREE = 1
      CALL MEMADD(KXA,   NASHT,  KFREE,2)
      CALL MEMADD(KXB,   NASHT,  KFREE,2)
      CALL MEMADD(KSCR,  2*NASHT,KFREE,2)
      CALL MEMADD(KH1DIA,NASHT,  KFREE,2)
      CALL MEMADD(KRJ,   N2ASHX, KFREE,2)
      CALL MEMADD(KRK,   N2ASHX, KFREE,2)
      CALL MEMADD(KLSCR, N2ASHX, KFREE,2)
      IF (IPRDIA .GT. 15) THEN
         WRITE (LUPRI,'(/2A/8I8)') ' CIDIAD : local pointers',
     *      ' KXA KXB KSCR KH1DIA KRJ KRK KFREE LFREE',
     *        KXA,KXB,KSCR,KH1DIA,KRJ,KRK,KFREE,LFREE
         WRITE (LUPRI,*) '  ICSYM IN CIDIAD ', ICSYM
      END IF
      IF (KFREE .GT. LFREE) CALL ERRWRK('CIDIAD',KFREE,LFREE)
C
C     Pack information as wanted by CIDIA4
C
      DO 100 I = 1,NASHT
         WRK(KXA-1+I) = FCAC( I*(I+1)/2 )
  100 CONTINUE
      CALL REORMT(WRK(KXA),WRK(KH1DIA),NASHT,1,XNDXCI(KLTSOB),1)
C     ... diagonal of modified one-body integrals
C
      CALL REORMT(FIJ,WRK(KRJ),NASHT,NASHT,XNDXCI(KLTSOB),
     &            XNDXCI(KLTSOB) )
      CALL REORMT(FJI,WRK(KRK),NASHT,NASHT,XNDXCI(KLTSOB),
     &            XNDXCI(KLTSOB) )
C     ... Coulomb and exchange integrals, FIJ and FJI
C
C
C
      IF ( IPRDIA .GE. 15 ) THEN
         WRITE(LUPRI,'(/A)')
     *      ' CIDIAD : Diagonal of FCAC, the modified one body matrix'
         WRITE(LUPRI,'(5F15.8)') (WRK(KH1DIA-1+I),I=1,NASHT)
         WRITE(LUPRI,'(/A)')
     *      ' CIDIAD : Coulomb  integrals RJ(u,v) = (uu/vv)'
         CALL OUTPUT(WRK(KRJ),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A)')
     *      ' CIDIAD : Exchange integrals RK(u,v) = (uv/uv)'
         CALL OUTPUT(WRK(KRK),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
C
C
      NCDET = NDTASM(ICSYM)
      CALL GETTIM(T1,W1)
         CALL CIDIA4(NAEL,NASTR,XNDXCI(KIASTR),
     &               NBEL,NBSTR,XNDXCI(KIBSTR),
     &               NASHT,DIAG,NCDET,XNDXCI(KSTASA),XNDXCI(KSTASB),
     &               MAXSYM,WRK(KH1DIA),
     &               XNDXCI(KSTBAA),XNDXCI(KSTBAB),
     &               ICSYM,WRK(KXA),WRK(KXB),WRK(KSCR),WRK(KRJ),
     &               WRK(KRK),MULD2H,XNDXCI(KNSSOA),XNDXCI(KNSSOB),
     &               XNDXCI(KIOCOC),NOCTPA,NOCTPB,XNDXCI(KISSOA),
     &               XNDXCI(KISSOB),ECORE,XNDXCI(KICOOS),IPRDIA )
C        CALL CIDIA4(NAEL,NASTR,IASTR,NBEL,NBSTR,IBSTR,
C    &               NORB,DIAG,NDET,NSTASA,NSTASB,
C    &               MAXSYM,H,ISTBAA,ISTBAB,
C    &               ICSYM,XA,XB,SCR,RJ,RK,
C    &               SYMPRO,NSSOA,NSSOB,IOCOC,NOCTPA,NOCTPB,
C    &               ISSOA,ISSOB,ECORE,ICOOS,NTEST)
         CALL GETTIM(T2,W2)
         IF (IPRDIA .GE. 1) WRITE (LUPRI,'(/A,2I10/A,2F10.3)')
     *   ' CIDIAD.cidia4 : KFREE, LFREE =           ',KFREE,LFREE,
     *   '                 CPU and wall time (sec) :',T2-T1,W2-W1
C
C... end of cidiad so :
      RETURN
      END
C  /* Deck cisigc */
      SUBROUTINE CISIGC(NSIM,BCVECS,SCVECS,LSCVEC,
     *                  FCAC,H2AC,XNDXCI,WRK,LFREE)
C
C Written 18-JAN-1988 by J.O.
C   ( from hjaaj's cisigc )
C
C Purpose:
C   Compute CI sigma vector(s) for direct CI,
C   the inactive energy EMY is not included.
C   calls lunar determinant routines
C
C Output:
C   SCVECS; CI sigma vector(s)
C
C Scratch:
C   WRK(LFREE)
C
C
#include <implicit.h>
      DIMENSION BCVECS(NCONST,*),FCAC(*),H2AC(*),SCVECS(LSCVEC,*)
      DIMENSION XNDXCI(*), WRK(LFREE)
#include <iratdef.h>
      LOGICAL   NOH2, IH8SM
      PARAMETER ( NOH2 = .FALSE. , IH8SM = .TRUE. , D0 = 0.0D0 )
C     NOH2  ... 2-electron part to be included
C     IH8SM ... integrals have 8-fold symmetry
C
C Used from common blocks:
C   INFORB : NASHT, N2ASHX
C   INFLIN : LSYMST,NCONST
C   INFTIM : NCALLS,TIMCPU,TIMWAL    ! IDTIM is index for these
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <inflin.h>
C ---
      PARAMETER (IDTIM = 2)
#include <inftim.h>
C
C
C *** SCVECS(I,K) = SUM(J) OF L_CI(I,J)*BCVECS(J,K)
C
C
      CALL GETTIM(T0,W0)
C
C  ** One- and two-electron contributions
C
      CALL GETTIM(T1,W1)
      CALL GETTIM(T2,W2)
C
C     Unpack FCAC
C
      KUFCAC = 1
      KW1    = KUFCAC + N2ASHX
      LW1    = LFREE  - KW1
      CALL DSPTSI(NASHT,FCAC,WRK(KUFCAC))
C     ... unpack FCAC using CALL DSPTSI(N,ASP,ASI)
      ISPIN1 = 0
      ISPIN2 = 0
C     ... singlet-singlet coupling of 2-electron integrals
      DO 100 ISIM = 1, NSIM
         CALL DZERO(SCVECS(1,ISIM),NCONST)
         CALL CISIGD(LSYMST,LSYMST,NCONST,NCONST,
     &               BCVECS(1,ISIM),SCVECS(1,ISIM), WRK(KUFCAC),H2AC,
     &               NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK(KW1),LW1)
C        CALL CISIGD(ICSYM,IHCSYM,NCDET,NHCDET, C,HC, FCAC,H2AC,
C    &               NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK,LFREE)
  100 CONTINUE
C
      CALL GETTIM(T3,W3)
C
C
C     Acuumulate timing for TIMOPT
C
      NCALLS(  IDTIM) = NCALLS(  IDTIM) + 1
      NVECS (  IDTIM) = NVECS (  IDTIM) + NSIM
      TIMCPU(1,IDTIM) = TIMCPU(1,IDTIM) + T3 - T0
      TIMCPU(2,IDTIM) = TIMCPU(2,IDTIM) + T1 - T0
      TIMCPU(3,IDTIM) = TIMCPU(3,IDTIM) + T2 - T1
      TIMCPU(4,IDTIM) = TIMCPU(4,IDTIM) + T3 - T2
      TIMWAL(1,IDTIM) = TIMWAL(1,IDTIM) + W3 - W0
      TIMWAL(2,IDTIM) = TIMWAL(2,IDTIM) + W1 - W0
      TIMWAL(3,IDTIM) = TIMWAL(3,IDTIM) + W2 - W1
      TIMWAL(4,IDTIM) = TIMWAL(4,IDTIM) + W3 - W2
C
C
C
C *** End of subroutine CISIGC
C
      RETURN
      END
C  /* Deck ciprp */
      SUBROUTINE CIPRP(NSIM,CREF,SCVECS,LSCVEC,PRPAC,XNDXCI,WRK,LFREE)
C
C Written Dec 1989 hjaaj
C
C Purpose:
C   Compute CI sigma vector(s) for NSIM one-electron operators
C   SCVECS(I) = SUM(J) <I|PRP|J> * CREF(J)
C
C   calls lunar determinant routines
C
C Output:
C   SCVECS; CI sigma vector(s)
C
C Scratch:
C   WRK(LFREE)
C
C
#include <implicit.h>
      DIMENSION CREF(NCONRF),PRPAC(N2ASHX,NSIM),SCVECS(LSCVEC,NSIM)
      DIMENSION XNDXCI(*), WRK(LFREE)
#include <iratdef.h>
      LOGICAL   NOH2, IH8SM
      PARAMETER ( NOH2 = .TRUE. , IH8SM = .FALSE. , D2 = 2.0D0 )
#include <thrzer.h>
C     NOH2  ... No 2-electron part
C     IH8SM ... integrals do not have 8-fold symmetry
C
C Used from common blocks:
C   INFORB : N2ASHX
C   INFLIN : LSYMRF,NCONRF,LSYMST,NCONST
C
#include <maxorb.h>
#include <inforb.h>
#include <inflin.h>
C
C
C *** SCVECS(I) = SUM(J) <I|PRP|J> * CREF(J)
C
C
      KW1    = 1
      LW1    = LFREE  - KW1
C
C     Set spin couplings
C
      ISPIN1 = 0
      ISPIN2 = 0
C     ... singlet-singlet coupling of 2-electron integrals
      DO 100 ISIM = 1, NSIM
         CALL DZERO(SCVECS(1,ISIM),NCONST)
         CALL CISIGD(LSYMRF,LSYMST,NCONRF,NCONST, CREF,SCVECS(1,ISIM),
     &               PRPAC(1,ISIM),DUMMY,NOH2,IH8SM,
     &               XNDXCI,ISPIN1,ISPIN2,WRK(KW1),LW1)
C        CALL CISIGD(ICSYM,IHCSYM,NCDET,NHCDET, C,HC, FCAC,H2AC,
C    &               NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK,LFREE)
  100 CONTINUE
C
C
C *** End of subroutine CIPRP
C
      RETURN
      END
C  /* Deck cisigd */
      SUBROUTINE CISIGD(ICSYM,IHCSYM,NCVEC,NHCVEC,C,HC,
     *                  FCAC,H2AC,ONLYH1,IH8SM,
     *                  XNDXCI,JSPIN1,JSPIN2,WRK,LFREE)
C
C  MASTER ROUTINE FOR DIRECT CI CALCULATION
C
C MOTECC-90: The algorithms used in this module, CISIGD,
C            are described in Chapter 8 Section D.3 of
C            MOTECC-90 "Direct CI for RAS Expansions"
C
C
C  SOME INPUT :
C        C  : CI vector, of length NCVEC
C     FCAC  : ONE-ELECTRON HAMILTONIAN MODIFIED FOR CORE ELECTRONS
C             *** FULL MATRIX ***
C
C     XNDXCI: ARRAY CONTAINING STRING INFORMATION
C             AS OBTAINED IN DETINF
C OUTPUT :
C        HC : H TIMES C, OF SYMMETRY IHCSYM AND LENGTH NHCVEC
C
C
#include <implicit.h>
      DIMENSION C(*), HC(*)
      DIMENSION FCAC(*), H2AC(*), WRK(LFREE), XNDXCI(*)
      LOGICAL   ONLYH1,IH8SM
C
C   THRSML : THRESHOLD FOR INTEGRALS AND CI COEFFICIENTS TO BE ZAPPED
C
      PARAMETER ( D1 = 1.0D0, THRSML = 1.0D-14 )
C
#include <priunit.h>
C
C Used from common blocks:
C   INFINP : FLAG(27), LSYM,ISPIN
C   INFPRI : IPRSIG
C
C   CIINFO : NDTASM(*)
C   CSFBAS : CSF core allocation for XNDXCI
C   CBESPN : ISPIN1, ISPIN2
C
#include <maxorb.h>
#include <infinp.h>
#include <infpri.h>
C
#include <ciinfo.h>
#include <csfbas.h>
#include <cbespn.h>
C
C
      LOGICAL CSFEXP
C
C     Transfer spin-couplings to CBESPN
C
      ISPIN1 = JSPIN1
      ISPIN2 = JSPIN2
C
C
C     980820-hjaaj: now use CSF for LSYM and det.s for other sym.s
C
      CSFEXP = .NOT.FLAG(27) .AND.
     &         (ICSYM .EQ. LSYM .OR. IHCSYM .EQ. LSYM)
C
      IF (CSFEXP) THEN
         IERR = 0
         ICCSF = 0
         IF ( ICSYM.EQ.LSYM ) THEN
C           ABACUS and RESPONSE use CSF for singlet and det for triplet,
C           we find what is the case for this call by comparing
C           to NCSASM and NDTASM:
            IF ( NCVEC .EQ. NCSASM( ICSYM) ) THEN
               ICCSF = 1
            ELSE IF ( NCVEC .NE. NDTASM( ICSYM) ) THEN
               IERR=IERR+1
            END IF
         ELSE
            IF ( NCVEC .NE. NDTASM( ICSYM) ) IERR=IERR+1
         END IF
         IHCCSF = 0
         IF (IHCSYM.EQ.LSYM ) THEN
C           ABACUS and RESPONSE use CSF for singlet and det for triplet,
C           we find what is the case for this call by comparing
C           to NCSASM and NDTASM:
            IF (NHCVEC .EQ. NCSASM(IHCSYM) ) THEN
               IHCCSF = 1
            ELSE IF (NHCVEC .NE. NDTASM(IHCSYM) ) THEN
               IERR=IERR+1
            END IF
         ELSE
            IF (NHCVEC .NE. NDTASM(IHCSYM) ) IERR=IERR+1
         END IF
         IF (IERR .GT. 0) THEN
            WRITE (LUPRI,*) 'CSF ERROR in CISIGD, LSYM =',LSYM
            WRITE (LUPRI,*)' ISPIN, ISPIN1, ISPIN2:',
     &           ISPIN,ISPIN1,ISPIN2
            WRITE (LUPRI,*)
     *         ' ICSYM,  NCVEC, NCSASM( ICSYM), NDTASM( ICSYM):',
     *           ICSYM,  NCVEC, NCSASM( ICSYM), NDTASM( ICSYM)
            WRITE (LUPRI,*)
     *         'IHCSYM, NHCVEC, NCSASM(IHCSYM), NDTASM(IHCSYM):',
     *          IHCSYM, NHCVEC, NCSASM(IHCSYM), NDTASM(IHCSYM)
            CALL QUIT('NCVEC/NHCVEC disagree with NCSASM(:) in cisigd')
         END IF
         IF (ICCSF .EQ. 0 .AND. IHCCSF .EQ. 0) CSFEXP = .FALSE.
      END IF
C
      IF (CSFEXP) THEN
C        .. change from CSf basis to determinant basis
         NCDET  = NDTASM(ICSYM)
         NHCDET = NDTASM(IHCSYM)
         NDET   = MAX(NCDET,NHCDET)
         KFREE  = 1
         CALL MEMADD(KCDET,NDET,KFREE,2)
         CALL MEMADD(KHCDET,NDET,KFREE,2)
         IF (ICCSF .EQ. 1) THEN
            CALL COPVEC(C,WRK(KHCDET),NCVEC)
            CALL CSDTVC(WRK(KHCDET),WRK(KCDET),1,XNDXCI(KDTOC),
     &                  XNDXCI(KICTS(1)),ICSYM,0,IPRSIG)
         ELSE
            CALL COPVEC(C,WRK(KCDET),NCVEC)
         END IF
         CALL SETVEC(WRK(KHCDET),0.0D0,NHCDET)
         LW1   = LFREE - KFREE
         CALL CISGD2(ICSYM,IHCSYM,NCDET,NHCDET,WRK(KCDET),WRK(KHCDET),
     *               FCAC,H2AC,ONLYH1,IH8SM,XNDXCI,WRK(KFREE),LW1)
         IF (IHCCSF .EQ. 1) THEN
            CALL CSDTVC(WRK(KCDET),WRK(KHCDET),2,XNDXCI(KDTOC),
     &                  XNDXCI(KICTS(1)),IHCSYM,0,IPRSIG)
            CALL DAXPY(NHCVEC,D1,WRK(KCDET),1,HC,1)
         ELSE
            CALL DAXPY(NHCVEC,D1,WRK(KHCDET),1,HC,1)
         END IF
      ELSE
         IF (NCVEC  .NE. NDTASM( ICSYM) .OR.
     *       NHCVEC .NE. NDTASM(IHCSYM)) THEN
            WRITE (LUPRI,*) 'ERROR in CISIGD'
            WRITE (LUPRI,*) ' ICSYM,  NCVEC, NDTASM( ICSYM):',
     *                        ICSYM,  NCVEC, NDTASM( ICSYM)
            WRITE (LUPRI,*) 'IHCSYM, NHCVEC, NDTASM(IHCSYM):',
     *                       IHCSYM, NHCVEC, NDTASM(IHCSYM)
            CALL QUIT('NCVEC/NHCVEC disagree with NDTASM(:) in cisigd')
         END IF
         CALL CISGD2(ICSYM,IHCSYM,NCVEC,NHCVEC,C,HC,
     *               FCAC,H2AC,ONLYH1,IH8SM,XNDXCI,WRK,LFREE)
      END IF
      RETURN
      END
C  /* Deck cisgd2 */
      SUBROUTINE CISGD2(ICSYM,IHCSYM,NCDET,NHCDET,C,HC,
     *                  FCAC,H2AC,ONLYH1,IH8SM,XNDXCI,WRK,LFREE)
C
C  MASTER ROUTINE FOR DIRECT CI CALCULATION
C
C  SOME INPUT :
C        C  : VECTOR TO BE MULTIPLIED WITH H,
C             SYMMETRY OF C IS ICSYM, LENGTH IS NCDET
C     FCAC  : ONE-ELECTRON HAMILTONIAN MODIFIED FOR CORE ELECTRONS
C             *** FULL MATRIX ***
C
C     XNDXCI: ARRAY CONTAINING STRING INFORMATION
C             AS OBTAINED IN DETINF
C OUTPUT :
C       HC  : H TIMES C, OF SYMMETRY IHCSYM AND LENGTH NHCDET
C
C
#include <implicit.h>
      DIMENSION C(*), HC(*), FCAC(*)
      DIMENSION H2AC(*), WRK(LFREE), XNDXCI(*)
      LOGICAL   ONLYH1,IH8SM,PERMSM
C
C   THRSML : THRESHOLD FOR INTEGRALS AND CI COEFFICIENTS TO BE ZAPPED
C
      PARAMETER ( THRSML = 1.0D-14 )
C
#include <priunit.h>
#include <mxsmob.h>
C
C Used from common blocks:
C   INFINP : FLAG(66)
C   INFORB : MULD2H(8,8), NASHT, N2ASHX,NNASHX
C   INFPRI : IPRSIG
C
C   STRNUM : EQUAL,NAEL,NASTR,NAEXCI, NBEL,...
C   DETBAS : core allocation for XNDXCI
C   MXBLK  : MXSASM,MXVBLK
C   CBESPN : ISPIN1,ISPIN2
C   SPINFO : MS2,?
C CBGETDIS : DISTYP, IADINT
C   INFSPI : ISGN1,ISGN2
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
#include <mxpdim.h>
#include <strnum.h>
#include <detbas.h>
#include <mxblk.h>
#include <cbespn.h>
#include <spinfo.h>
#include <cbgetdis.h>
#include <infspi.h>
C
C
      CALL GETTIM(T0,W0)
C
C    *** Set some control variables ***
C
C  IDOH2 : <> 0  ' NORMAL ' DIRECT CI WITH ONE- AND TWO- BODY TERMS
C        :  = 0             DIRECT CI WITH ONE-          BODY TERMS ONLY
C
C.Spin permutation simplifications
C 980826-hjaaj: added ICSYM.eq.IHCSYM test, PERMSM only works
C  for totally symmetric operators.
C 990427-hjaaj: tests show PERMSM does not always work in this version,
C  we therefore never use PERMSM.
Chj   IF(MS2.EQ.0.AND..NOT.FLAG(66) .AND. ICSYM.EQ.IHCSYM
Chj  &  .AND.ISPIN1.EQ.0.AND.ISPIN2.EQ.0) THEN
Chj     PERMSM = .TRUE.
Chj     IF(MOD((MULTS+1)/2,2).EQ.1) THEN
Chj       PSIGN = 1.0D0
Chj     ELSE
Chj       PSIGN = -1.0D0
Chj     END IF
Chj   ELSE
        PERMSM = .FALSE.
        PSIGN = 999999999.D0
Chj   END IF
C
      IF (ONLYH1) THEN
         IDOH2 = 0
      ELSE
         IDOH2 = 1
      END IF
C
C     INTFRM = 1 : use GETINT
C            = 3 : get integrals directly from H2AC
C
      IF (DISTYP .EQ. 1 .AND. IADINT .LT. 0) THEN
         INTFRM = 3
      ELSE
         INTFRM = 1
      END IF
C
C     Make sure ISGN1/ISGN2 is defined for the DISTYP's where
C     they are needed in GETIN2. They are now only defined in
C     RSP2GR /950524-hjaaj
C
      ISGN1 = (-1)**ISPIN1
      ISGN2 = (-1)**ISPIN2
C
C core energy neglected
      ECORE  = 0.0D0
c. Store determinant CI off-sets in C arrays and HC arrays
      CALL CIOFF(ICSYM ,1,XNDXCI,IPRSIG)
      CALL CIOFF(IHCSYM,2,XNDXCI,IPRSIG)
C     CALL CIOFF(IREFSM,ICORHC,XNDXCI,NTEST)
C
C** 2 : DO H TIMES C
c
      KFREE = 1
      CALL MEMADD(KUFCAC,N2ASHX, KFREE,2)
      CALL MEMADD(KRIJKL,N2ASHX*MXSMOB, KFREE,2)
      CALL MEMADD(KVEC1,MXSASM,KFREE,2)
      CALL MEMADD(KVEC2,MXSASM*MXSMOB,KFREE,2)
      CALL MEMADD(KVEC3,MXSASM*MXSMOB,KFREE,2)
      CALL MEMADD(KINDX1,MXSASM,KFREE,1)
      CALL MEMADD(KINDX2,MXSASM,KFREE,1)
      CALL MEMADD(KL,MXSASM*MXSMOB,KFREE,1)
      CALL MEMADD(KR,MXSASM*MXSMOB,KFREE,1)
      CALL MEMADD(KC2,MXVBLK,KFREE,2)
      CALL MEMADD(KWIJKL,N2ASHX,KFREE,2)
      CALL MEMADD(KINDE2,MXPST*MXPTP,KFREE,1)
      CALL MEMADD(KF3,MXPST*MXPTP,KFREE,2)
      IF (KFREE .GT. LFREE) CALL ERRWRK('CISIGD',KFREE,LFREE)
c
c     Reorder FCAC from Sirius order to Lunar order.
c
      CALL REORMT(FCAC,WRK(KUFCAC),NASHT,NASHT,XNDXCI(KLTSOB),
     &            XNDXCI(KLTSOB))
C
      CALL GETTIM(T1,W1)
C
      CALL CISIG9(XNDXCI(KTAIJ),XNDXCI(KTATO),NAEL,NASTR,NAEXCI,
     &            XNDXCI(KTBIJ),XNDXCI(KTBTO),NBEL,NBSTR,NBEXCI,
     &     NASHT,C,HC,NCDET,NHCDET,XNDXCI(KSTASA),XNDXCI(KSTASB),
     &     XNDXCI(KCOFF),XNDXCI(KHCOFF),MAXSYM,XNDXCI(KISSYM),
     &     XNDXCI(KSTBAA),XNDXCI(KSTBAB),XNDXCI(KTASYM),XNDXCI(KTBSYM),
     &     WRK(KUFCAC),WRK(KRIJKL),ICSYM,IHCSYM,PERMSM,
     &     XNDXCI(KIOCOC),NOCTPA,NOCTPB,
     &     XNDXCI(KICSO), XNDXCI(KIHCSO),XNDXCI(KICOOS),XNDXCI(KIHOOS),
     &     XNDXCI(KNSSOA),XNDXCI(KISSOA),XNDXCI(KNSSOB),XNDXCI(KISSOB),
     &     XNDXCI(KKLTP),
     &     XNDXCI(KICREA),XNDXCI(KIANNI),WRK(KVEC1),WRK(KVEC2),
     &     WRK(KC2),WRK(KINDX1),WRK(KINDX2),WRK(KL),
     &     WRK(KR),PSIGN,WRK(KVEC3),XNDXCI(KTPFSA),XNDXCI(KTPFSB),
     &     XNDXCI(KCDTAS),XNDXCI(KHDTAS),XNDXCI(KKLCAN),XNDXCI(KTPFOB),
     &     HC,IDOH2,ECORE,H2AC,WRK(KWIJKL),IPRSIG,XNDXCI(KIASTR),
     &     XNDXCI(KIBSTR),MXSASM,XNDXCI(KLTSOB),XNDXCI(KSTLOB),
     &     ISPIN1,ISPIN2,IH8SM,INTFRM,WRK(KINDE2),WRK(KF3))
      CALL GETTIM(T2,W2)
      IF (IPRSIG .GE. 1) WRITE (LUPRI,'(/A,2I10/A,2F10.3)')
     *   ' CISIGD : KFREE, LFREE =           ',KFREE,LFREE,
     *   '          CPU and wall time (sec) :',T2-T1,W2-W1
C
      RETURN
      END
C  /* Deck cisigo */
      SUBROUTINE CISIGO(NOSIM,SOVECS,CREF,EMYX,FXCAC,H2XAC,
     *                  XNDXCI,WRK,LFREE)
C
C  Written 18-Jan-1988 J.O.
C  (After cisigo for casguga of hjaaj )
C
C Parameter list:
C   NOSIM  number of simultaneous orbital expansion vectors
C   CREF   reference CI-vector
C   SOVECS NOSIM sigma vectors of orbital trial vectors
C   EMYX   the inactive "energy" from the one-index transformed
C          "Hamiltonian"
C   FXCAC  1-ind. transf. inactive Fock matrix
C   H2XAC  1-ind. transf. active 2-el. integrals
C   XNDXCI CI information
C
C Output:
C   SOVECS(I) = SUM(sr) of K(I,sr)*BOVECS(sr)
C
C Scratch:
C   WRK work area containing :
C
C
#include <implicit.h>
C
      DIMENSION CREF(NCONRF), SOVECS(NVARPT,NOSIM),EMYX(NOSIM),XNDXCI(*)
      DIMENSION FXCAC(NNASHX,*),H2XAC(NNASHX,NNASHX,*)
      DIMENSION WRK(LFREE)
#include <iratdef.h>
      PARAMETER (IDTIM = 10)
      PARAMETER ( D2 = 2.0D0 )
      LOGICAL   NOH2, IH8SM
      PARAMETER ( NOH2 = .FALSE. , IH8SM = .TRUE. )
C     NOH2  ... 2-electron part to be included
C     IH8SM ... integrals have 8-fold symmetry
C
C
C Used from common blocks:
C
C   INFORB : NASHT,NNASHX,N2ASHX
C   INFLIN : LSYMRF,LSYMPT,LSYMST,NCONRF,NCONST
C   INFTIM : NCALLS,TIMCPU,TIMWAL    ! IDTIM is index for these
C
#include <maxorb.h>
#include <inforb.h>
#include <inflin.h>
#include <inftim.h>
C
C
      KUFCAC = 1
      KW1    = KUFCAC + N2ASHX
      LW1    = LFREE  - KW1
C
      ISPIN1 = 0
      ISPIN2 = 0
C     ... singlet-singlet coupling of 2-electron integrals
C
C
C *** Calculate sigma vectors with modified integrals
C
      CALL GETTIM(T0,W0)
      CALL GETTIM(T1,W1)
      CALL GETTIM(T2,W2)
      CALL GETTIM(T3,W3)
C
C
      DO 300 IOSIM = 1, NOSIM
C
         CALL DSPTSI(NASHT,FXCAC(1,IOSIM),WRK(KUFCAC))
C        ... unpack FXCAC using CALL DSPTSI(N,ASP,ASI)
         CALL DZERO(SOVECS(1,IOSIM),NCONST)
         CALL CISIGD(LSYMRF,LSYMST,NCONRF,NCONST,
     &               CREF,SOVECS(1,IOSIM),WRK(KUFCAC),H2XAC(1,1,IOSIM),
     &               NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK(KW1),LW1)
C        CALL CISIGD(ICSYM,IHCSYM,NCDET,NHCDET, C,HC, FCAC,H2AC,
C    &               NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK,LFREE)
C
C
C --- Add inactive energy contributions
C --- multiply orbital sigma vectors by 2 to get final SOVECS
C
      IF (LSYMPT .EQ. 1) THEN
         DO 200 I = 1,NCONST
            SOVECS(I,IOSIM) = D2*(SOVECS(I,IOSIM) + CREF(I)*EMYX(IOSIM))
  200    CONTINUE
      ELSE
         CALL DSCAL(NCONST,D2,SOVECS(1,IOSIM),1)
      END IF
C
  300 CONTINUE
C
C
      CALL GETTIM(T4,W4)
C
      NCALLS(  IDTIM) = NCALLS(  IDTIM) + 1
C!!!      NVECS (  IDTIM) = NVECS (  IDTIM) + NOSIM
      TIMCPU(1,IDTIM) = TIMCPU(1,IDTIM) + T4 - T0
      TIMCPU(2,IDTIM) = TIMCPU(2,IDTIM) + T1 - T0
      TIMCPU(3,IDTIM) = TIMCPU(3,IDTIM) + T2 - T1
      TIMCPU(4,IDTIM) = TIMCPU(4,IDTIM) + T3 - T2
      TIMCPU(5,IDTIM) = TIMCPU(5,IDTIM) + T4 - T3
      TIMWAL(1,IDTIM) = TIMWAL(1,IDTIM) + W4 - W0
      TIMWAL(2,IDTIM) = TIMWAL(2,IDTIM) + W1 - W0
      TIMWAL(3,IDTIM) = TIMWAL(3,IDTIM) + W2 - W1
      TIMWAL(4,IDTIM) = TIMWAL(4,IDTIM) + W3 - W2
      TIMWAL(5,IDTIM) = TIMWAL(5,IDTIM) + W4 - W3
C
C
C *** End of subroutine CISIGO
C
      RETURN
      END
