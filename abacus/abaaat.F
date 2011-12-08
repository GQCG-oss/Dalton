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
C  /* Deck aatinp */
C***********************************************************************
      SUBROUTINE AATINP(WORD)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 10)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <cbiaat.h>
C
#include <abainf.h>
C
      SAVE SET
      DATA TABLE /'.SKIP  ','.PRINT ','.STOP  ','.INTPRI',
     &            '.NONUC ','.NOELC ','.NODBDR','.XXXXXX',
     &            '.NOSEC ','.NODDY '/
C
      NEWDEF = (WORD .EQ. '*AAT   ')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 ITABLE = 1, NTABLE
                  IF (TABLE(ITABLE) .EQ. WORD) THEN
                     GOTO (1,2,3,4,5,6,7,8,9,10), ITABLE
                  END IF
 200           CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                  CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                       LUPRI)
                  GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &              '" not recognized in AATINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                    LUPRI)
               CALL QUIT('Illegal keyword in AATINP.')
 1             CONTINUE
                  SKIP = .TRUE.
                  GOTO 100
 2             CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
                  GOTO 100
 3             CONTINUE
                  CUT = .TRUE.
                  GOTO 100
 4             CONTINUE
                  READ (LUCMD,*) INTPRI
                  IF (INTPRI .EQ. IPRDEF) ICHANG = ICHANG -1
                  GOTO 100
 5             CONTINUE
                  NONUC = .TRUE.
                  GOTO 100
 6             CONTINUE
                  NOELC = .TRUE.
                  GOTO 100
 7             CONTINUE
                  NODBDR = .TRUE.
                  GOTO 100
 8             CONTINUE
                  GOTO 100
 9             CONTINUE
                  NOSEC  = .TRUE.
                  GOTO 100
10             CONTINUE
                  NODDY  = .TRUE.
                  GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "', WORD,
     &              '" not recognized in AATINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',
     &                    LUPRI)
               CALL QUIT('Illegal keyword in AATINP.')
            END IF
      END IF
 300  CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for AAT:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' AAT skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &           ' Print level in AAT:', IPRINT
            IF (INTPRI .NE. IPRDEF) WRITE (LUPRI,'(A,I5)')
     &           ' Print level of integrals in AAT:', INTPRI
            IF (NONUC) WRITE (LUPRI,'(A)')
     &           ' Nuclear contributions to AAT are NOT calculated'
            IF (NOELC) WRITE (LUPRI,'(A)')
     &           ' Electronic contributions to AAT are NOT calculated'
            IF (NODBDR) WRITE (LUPRI,'(A)')
     &           ' DBDR contribution to AAT is NOT calculated'
            IF (NOSEC) WRITE (LUPRI,'(A)')
     &           ' Orbital second order  contribution to AAT is NOT'
     &          //' calculated'
            IF (NODDY) WRITE (LUPRI,'(A)')
     &           ' AATs are calculated using noddy code for comparison.'
         END IF
         IF (CUT) THEN
            WRITE (LUPRI,'(/,A,/)')
     &           ' Program is stopped after AAT'
         END IF
      END IF
      RETURN
      END
C  /* Deck aatini */
      SUBROUTINE AATINI
C
C     Initialize CBIAAT
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbiaat.h>
C
      IPRINT = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      INTPRI = IPRDEF
      NONUC  = .FALSE.
      NOELC  = .FALSE.
      NODBDR = .FALSE.
      NOSEC  = .FALSE.
      NODDY  = .FALSE.
      RETURN
      END
C  /* Deck aatdrv */
C***********************************************************************
      SUBROUTINE AATDRV(WORK,LWORK,PASS)
C***********************************************************************
C
C     Calculation of Atomic Axial Tensor (AAT) in cartesian coordinates
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      LOGICAL PASS
      DIMENSION WORK(LWORK)
#include <nuclei.h>
#include <cbiaat.h>
#include <abainf.h>
C
      IF (SKIP .OR. .NOT.VCD) RETURN
      CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPRINT .GT. 2) CALL TITLER('Output from AATDRV','*',103)
C
C     Calculate Atomic Axial Tensor (AAT) in cartesian coordinates
C
      KCSTRA = 1
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('AATDRV',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST + 1
      CALL AATDR1(WORK(KCSTRA),WORK(KSCTRA),WORK(KLAST),LWRK)
C
      CALL TIMER('AATDRV',TIMSTR,TIMEND)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/A/)')
     &        ' Program stopped after AATDRV as required.'
         WRITE (LUPRI,'(/A/)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in AATDRV) *****')
      END IF
      RETURN
      END
C  /* Deck aatdr1 */
C***********************************************************************
      SUBROUTINE AATDR1(CSTRA,SCTRA,WORK,LWORK)
C***********************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION WORK(LWORK), CSTRA(*), SCTRA(*)
#include <nuclei.h>
#include <inforb.h>
#include <aatens.h>
#include <cbiaat.h>
C
      CALL DZERO(AATNUC,9*NUCDEP)
      CALL DZERO(AATORB,9*NUCDEP)
      CALL DZERO(AATCI ,9*NUCDEP)
      CALL DZERO(AAT2ND,9*NUCDEP)
C
C     Calculate and add nuclear part Atomic Axial Tensor (AAT)
C
      IF (.NOT. NONUC) CALL NUCAAT(CSTRA,SCTRA,IPRINT)
C
C     Calculate and add electronic part of AAT
C
      IF (.NOT. NOELC) CALL AATELC(CSTRA,SCTRA,WORK,LWORK)
C
      CALL DZERO(AATTOT,9*NUCDEP)
      CALL ADDAAT(AATNUC)
      CALL ADDAAT(AATORB)
      CALL ADDAAT(AATCI )
      CALL ADDAAT(AAT2ND)
      IF (IPRINT .GT. 0) THEN
         CALL HEADER('Total AAT in AATDR1',-1)
         CALL FCPRI(AATTOT,'AAT',CSTRA,SCTRA)
      END IF
      RETURN
      END
C  /* Deck nucaat */
C***********************************************************************
      SUBROUTINE NUCAAT(CSTRA,SCTRA,IPRINT)
C***********************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (DP25 = 0.25D0)
      INTEGER X, Y, Z
      DIMENSION CSTRA(*), SCTRA(*)
#include <symmet.h>
#include <nuclei.h>
#include <aatens.h>
      NEXT(I) = MOD(I,3) + 1
C
      DO 100 IATOM = 1, NUCIND
         FAC = DP25*FMULT(ISTBNU(IATOM))*CHARGE(IATOM)
         DO 200 X = 1, 3
            Y  = NEXT(X)
            Z  = NEXT(Y)
            IY = IPTCNT(3*(IATOM - 1) + Y,ISYMAX(X,2),1)
            IZ = IPTCNT(3*(IATOM - 1) + Z,ISYMAX(X,2),1)
            IF (IY .GT. 0) AATNUC(IPTAX(X,2),IY) = - FAC*CORD(Z,IATOM)
            IF (IZ .GT. 0) AATNUC(IPTAX(X,2),IZ) =   FAC*CORD(Y,IATOM)
  200    CONTINUE
  100 CONTINUE
C
      IF (IPRINT .GT. 0) THEN
         CALL HEADER('Nuclear contribution to AAT',-1)
         CALL FCPRI(AATNUC,'AAT',CSTRA,SCTRA)
      END IF
      RETURN
      END
C  /* Deck aatelc */
C***********************************************************************
      SUBROUTINE AATELC(CSTRA,SCTRA,WORK,LWORK)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION CSTRA(*), SCTRA(*), WORK(LWORK)
#include <inforb.h>
C
      KCMO  = 1
      KDTOT = KCMO  + NCMOT
      KLAST = KDTOT + N2ORBX
      IF (KLAST .GT. LWORK) CALL STOPIT('AATELC',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
      CALL AATEL1(WORK(KCMO),WORK(KDTOT),CSTRA,SCTRA,WORK(KLAST),LWRK)
      RETURN
      END
C  /* Deck aatel1 */
C***********************************************************************
      SUBROUTINE AATEL1(CMO,DTOT,CSTRA,SCTRA,WORK,LWORK)
C***********************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
      LOGICAL OLDDX, FOUND
      CHARACTER*8 LABINT(3*MXCOOR)
      DIMENSION CMO(*), DTOT(*), CSTRA(*), SCTRA(*), WORK(LWORK)
#include <nuclei.h>
#include <inftap.h>
#include <inforb.h>
#include <cbiaat.h>
#include <abainf.h>
C
C     Open direct access file
C
      LUAAT = -1
      CALL GPOPEN(LUAAT,'ABAAAT','UNKNOWN','DIRECT',' ',IRAT*N2ORBX,
     &            OLDDX)
C
C     Read MO coefficients
C
      CALL RD_SIRIFC('CMO',FOUND,CMO,WORK,LWORK)
      IF (.NOT.FOUND) CALL QUIT('AATEL1 error: CMO not found on SIRIFC')
C
      IF (NODIFC) THEN
C
C        Construct DB, DR
C
         IF (.NOT. NODBDR) CALL DBDR(CMO,WORK,LWORK,LUAAT,LABINT,
     &                               INTPRI,IPRINT)
      END IF
C
C     Construct one-electron density matrix
C
      CALL MKDENS(DTOT,WORK,LWORK,IPRINT)
C
C     Construct first orbital contribution and CI contribution
C
      CALL AATOPC(DTOT,LUAAT,CSTRA,SCTRA,WORK,LWORK,IPRINT)
C
C     Second-order contribution
C
      IF (.NOT. NOSEC) CALL AATSEC(DTOT,CMO,LABINT,CSTRA,SCTRA,WORK,
     &                             LWORK,INTPRI,IPRINT)
C
C     Noddy code
C
      IF (NODDY) CALL AATNOD(DTOT,CMO,LABINT,CSTRA,SCTRA,WORK,LWORK,
     &                       INTPRI,IPRINT)
C
      CALL GPCLOSE(LUAAT,'DELETE')
      RETURN
      END
C  /* Deck dbdr */
C***********************************************************************
      SUBROUTINE DBDR(CMO,WORK,LWORK,LUAAT,LABINT,INTPRI,IPRINT)
C***********************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
      CHARACTER LABINT(*)*8
      DIMENSION CMO(*), WORK(LWORK)
#include <nuclei.h>
#include <inforb.h>
C
      KDBR   = 1
      KINTRP = KDBR + N2ORBX
      KINTDR = KINTRP + (3*NUCDEP + 1)/IRAT
      KWRK   = KINTDR + (3*NUCDEP + 1)/IRAT
      IF (KWRK .GT. LWORK) CALL STOPIT('DBDR',' ',LWORK,KWRK)
      LWRK = LWORK - KWRK + 1
C
      CALL DBDR1(WORK(KDBR),CMO,LABINT,WORK(KINTRP),WORK(KINTDR),
     &           LUAAT,WORK(KWRK),LWRK,INTPRI,IPRINT)
C
      RETURN
      END
C  /* Deck dbdr1 */
C***********************************************************************
      SUBROUTINE DBDR1(DBR,CMO,LABINT,INTREP,INTADR,LUAAT,WORK,
     &                 LWORK,INTPRI,IPRINT)
C***********************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
#include <inforb.h>
#include <nuclei.h>
#include <symmet.h>
#include <abainf.h>
C
      PARAMETER (D05 = 0.5D0)
      INTEGER X
      CHARACTER*8 LABINT(*)
      DIMENSION DBR(NORBT,NORBT), CMO(*), INTREP(*), INTADR(*),
     &          WORK(LWORK)
C
      CALL QENTER('DBDR1 ')
C
C     Construct symmetrized half B-derivative overlap
C     matrix (DB) and write to file
C
      NCOMP  = 0
      NPATOM = 0
      CALL GET1IN(DUMMY,'HBDO   ',NCOMP,WORK,LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &            .FALSE.,DUMMY,INTPRI)
C
      DO 300 X = 1, 3
         ISCOOR = 3*NUCDEP + IPTAX(X,2)
         IF (NOLOND) THEN
            CALL DZERO(DBR,NORBT*NORBT)
         ELSE
            KSYMP = -1
            CALL PRPGET(LABINT(X),CMO,DBR,KSYMP,ASM,WORK,LWORK,IPRINT)
            IF (KSYMP.NE.INTREP(ISCOOR) + 1)
     &      CALL QUIT('ERROR: unexpected symmetry of property matrix')
            CALL DSCAL(N2ORBX,D05,DBR,1)
         ENDIF
         CALL WRITDX(LUAAT,ISCOOR,IRAT*N2ORBX,DBR)
         IF (IPRINT .GT. 10) THEN
            CALL HEADER('MO matrix '//LABINT(X)//' in DBDR1',-1)
            CALL OUTPUT(DBR,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
  300 CONTINUE
C
C     Construct antisymmetrized half R-derivative overlap
C     matrix (DR) and write to file
C
      NCOMP  = 0
      NPATOM = 0
      CALL GET1IN(DUMMY,'HDO    ',NCOMP,WORK,LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,
     &            DUMMY,.FALSE.,DUMMY,INTPRI)
C
      DO 200 ISCOOR = 1, 3*NUCDEP
         KSYMP = -1
         CALL PRPGET(LABINT(ISCOOR),CMO,DBR,KSYMP,ASM,WORK,LWORK,IPRINT)
         IF (KSYMP.NE.INTREP(ISCOOR) + 1)
     &      CALL QUIT('ERROR: unexpected symmetry of property matrix')
         CALL WRITDX(LUAAT,ISCOOR,IRAT*N2ORBX,DBR)
         IF (IPRINT .GT. 10) THEN
            CALL HEADER('MO matrix '//LABINT(ISCOOR)//' in DBDR1',-1)
            CALL OUTPUT(DBR,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
  200 CONTINUE
      CALL QEXIT('DBDR1 ')
      RETURN
      END
C  /* Deck mkdens */
C***********************************************************************
      SUBROUTINE MKDENS(DTOT,WORK,LWORK,IPRINT)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
      DIMENSION DTOT(NORBT,NORBT), WORK(LWORK)
C
      KDV   = 1
      KLAST = KDV + NNASHX
      IF (KLAST .GT. LWORK) CALL STOPIT('MKDENS',' ',LWORK,KLAST)
      CALL MKDEN1(DTOT,WORK(KDV),IPRINT)
      RETURN
      END
C  /* Deck mkden1 */
C***********************************************************************
      SUBROUTINE MKDEN1(DTOT,DV,IPRINT)
C***********************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
C
#include <inftap.h>
#include <inforb.h>
#include <infind.h>
#include <iratdef.h>
C
      PARAMETER (D2 = 2.D0)
      INTEGER   U, U1, U2, V, V1, V2, UV
      DIMENSION DV(NNASHX), DTOT(NORBT,NORBT)
      LOGICAL   FOUND
C
      IF (NASHT .GT. 0) THEN
         CALL RD_SIRIFC('DV',FOUND,DV,DTOT,N2ORBX)
         IF (.NOT.FOUND)
     &      CALL QUIT('MKDEN1 error: DV not found on SIRIFC')
         IF (IPRINT .GT. 10) THEN
            CALL AROUND('DV matrix in MKDEN1')
            CALL OUTPAK(DV,NASHT,1,LUPRI)
         END IF
      END IF
C
C     Construct total density matrix
C
      CALL DZERO(DTOT,N2ORBX)
C
      DO 100 ISYMU = 1, NSYM
         NASHU = NASH(ISYMU)
         NORBQ = NORB(ISYMU)
         DO 110 I = 1, NISH(ISYMU)
            INACT = IORB(ISYMU) + I
            DTOT(INACT,INACT) = D2
  110    CONTINUE
         IF (NASHU.GT.0 .AND. NORBQ.GT.0) THEN
            IOFFU1 = IORB(ISYMU) + NISH(ISYMU)
            IOFFU2 = IASH(ISYMU)
            DO 120 U = 1, NASHU
               U1 = IOFFU1 + U
               U2 = IOFFU2 + U
               DO 130 V = 1, NASHU
                  V1 = IOFFU1 + V
                  V2 = IOFFU2 + V
                  IF (U2 .GT. V2) THEN
                     UV = IROW(U2) + V2
                  ELSE
                     UV = IROW(V2) + U2
                  END IF
                  DTOT(U1,V1) = DTOT(U1,V1) + DV(UV)
 130           CONTINUE
 120        CONTINUE
         END IF
 100  CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL AROUND('DTOT in MKDEN1')
         CALL OUTPUT(DTOT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck aatsec */
C***********************************************************************
      SUBROUTINE AATSEC(DTOT,CMO,LABINT,CSTRA,SCTRA,WORK,LWORK,
     &                  INTPRI,IPRINT)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
      CHARACTER*8 LABINT(*)
      DIMENSION DTOT(*), CMO(*), CSTRA(*), SCTRA(*), WORK(LWORK)
#include <nuclei.h>
#include <inforb.h>
      KDBR   = 1
      KDS1   = KDBR   + N2ORBX
      KINTRP = KDS1   + 3*N2ORBX
      KINTDR = KINTRP + (9*NUCDEP + 1)/IRAT
      KWRK   = KINTDR + (9*NUCDEP + 1)/IRAT
      IF (KWRK .GT. LWORK) CALL STOPIT('AATSEC',' ',LWORK,KWRK)
      LWRK = LWORK - KWRK + 1
      CALL AATSE1(DTOT,WORK(KDBR),WORK(KDS1),CMO,LABINT,
     &            WORK(KINTRP),WORK(KINTDR),CSTRA,SCTRA,WORK(KWRK),
     &            LWRK,INTPRI,IPRINT)
      RETURN
      END
C  /* Deck aatse1 */
C***********************************************************************
      SUBROUTINE AATSE1(DTOT,DBR,DS1,CMO,LABINT,INTREP,INTADR,CSTRA,
     &                  SCTRA,WORK,LWORK,INTPRI,IPRINT)
C***********************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
C
#include <inforb.h>
#include <nuclei.h>
#include <aatens.h>
#include <symmet.h>
#include <abainf.h>
C
      INTEGER X, XS
      CHARACTER*8 LABINT(*)
      DIMENSION DTOT(NORBT,NORBT), DBR(NORBT,NORBT), CMO(*),
     &          INTREP(*), INTADR(*), WORK(LWORK), DS1(N2ORBX,3),
     &          CSTRA(*), SCTRA(*)
C
      CALL QENTER('AATSE1')
      KLAST = 1
C
C     Construct half B-derivative overlap matrix and write to file
C
      CALL DZERO(AAT2ND,3*MXCOOR)
      IF (.NOT.NOLOND) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'HDOBR  ',NCOMP,WORK(KLAST),LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &            .FALSE.,DUMMY,INTPRI)
C
         DO 100 IREP = 0, MAXREP
         DO 110 IATOM = 1, NUCIND
C           IF (DOATOM(IATOM)) THEN
               DO 120 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     DO 130 X = 1, 3
                     IF (ISYMAX(X,2) .EQ. IREP) THEN
                        LSCOOR = 3*(ISCOOR - 1) + X
                        KSYMP = -1
                        CALL PRPGET(LABINT(LSCOOR),CMO,DBR,KSYMP,ASM,
     &                              WORK(KLAST),LWORK,IPRINT)
                        IF (KSYMP.NE.1) CALL QUIT(
     &                  'ERROR: unexpected symmetry of property matrix')
                        XS = IPTAX(X,2)
                        AAT2ND(XS,ISCOOR) = - DDOT(N2ORBX,DTOT,1,DBR,1)
                     END IF
 130                 CONTINUE
                  END IF
 120           CONTINUE
C           END IF
 110        CONTINUE
 100     CONTINUE
C
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('HDOBR contribution to AAT2ND in AATSE1',-1)
            CALL FCPRI(AAT2ND,'AAT',CSTRA,SCTRA)
         END IF
C
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'S1MAGR ',NCOMP,WORK(KLAST),LWORK,LABINT,
     &               INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &               .FALSE.,DUMMY,INTPRI)
         DO 200 X = 1, 3
            KSYMP = -1
            CALL PRPGET(LABINT(X),CMO,DBR,KSYMP,ASM,
     &                  WORK(KLAST),LWORK,IPRINT)
            IF (KSYMP.NE.INTREP(X) + 1)
     &      CALL QUIT('ERROR: unexpected symmetry of property matrix')
            CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
     &                 DTOT,NORBT,
     &                 DBR,NORBT,0.D0,
     &                 DS1(1,X),NORBT)
 200     CONTINUE
C
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'SQHDOL ',NCOMP,WORK(KLAST),LWORK,LABINT,
     &               INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &               .FALSE.,DUMMY,INTPRI)
         DO 300 IREP = 0, MAXREP
            DO 310 IATOM = 1, NUCIND
C           IF (DOATOM(IATOM)) THEN
               DO 320 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                  IF (ISCOOR .GT. 0) THEN
                     KSYMP = -1
                     CALL PRPGET(LABINT(ISCOOR),CMO,DBR,KSYMP,ASM,
     &                           WORK(KLAST),LWORK,IPRINT)
                     IF (KSYMP.NE.IREP + 1) CALL QUIT(
     &                  'ERROR: unexpected symmetry of property matrix')
                     DO 330 X = 1, 3
                     IF (ISYMAX(X,2) .EQ. IREP) THEN
                        XS = IPTAX(X,2)
                        AAT2ND(XS,ISCOOR) = AAT2ND(XS,ISCOOR)
     &                            + DDOT(N2ORBX,DS1(1,X),1,DBR,1)
                     END IF
 330                 CONTINUE
                  END IF
 320           CONTINUE
C           END IF
 310        CONTINUE
 300     CONTINUE
      END IF
C
      IF (IPRINT .GT. 0) THEN
         CALL HEADER('Total AAT2ND in AATSE1',-1)
         CALL FCPRI(AAT2ND,'AAT',CSTRA,SCTRA)
      END IF
      CALL QEXIT('AATSE1')
      RETURN
      END
C  /* Deck upwop */
C***********************************************************************
      SUBROUTINE UPWOP(PACKED,UNPACK,ANTI)
C***********************************************************************
C
C     Based on UPKWOP
C
#include <implicit.h>
      LOGICAL ANTI
#include <inforb.h>
#include <maxorb.h>
#include <infvar.h>
#include <inflin.h>
      DIMENSION UNPACK(NORBT,NORBT),PACKED(NWOPT)
C
      CALL DZERO(UNPACK,N2ORBX)
      IF (ANTI) THEN
         DO 100 IG = 1,NWOPPT
            I = JWOP(1,IG)
            J = JWOP(2,IG)
            UNPACK(I,J) =  PACKED(IG)
            UNPACK(J,I) = -PACKED(IG)
  100    CONTINUE
      ELSE
         DO 200 IG = 1,NWOPPT
            I = JWOP(1,IG)
            J = JWOP(2,IG)
            UNPACK(I,J) = PACKED(IG)
            UNPACK(J,I) = PACKED(IG)
  200    CONTINUE
      END IF
      RETURN
      END
C  /* Deck addaat */
C***********************************************************************
      SUBROUTINE ADDAAT(AMAT)
C***********************************************************************
#include <implicit.h>
#include <mxcent.h>
      DIMENSION AMAT(3,MXCOOR)
#include <nuclei.h>
#include <aatens.h>
      DO 100 I = 1, 3
         DO 200 J = 1, 3*NUCDEP
            AATTOT(I,J) = AATTOT(I,J) + AMAT(I,J)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck aatnod */
C***********************************************************************
      SUBROUTINE AATNOD(DTOT,CMO,LABINT,CSTRA,SCTRA,WORK,LWORK,
     &                  INTPRI,IPRINT)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <iratdef.h>

#include <inforb.h>
#include <inflin.h>
#include <nuclei.h>

      DIMENSION DTOT(NORBT,NORBT),CMO(*),CSTRA(*),SCTRA(*),WORK(LWORK)
      CHARACTER*8 LABINT(*)
C
      NCART  = 3*NUCDEP
C
      KINTRP = 1      + 3*NCART
      KINTDR = KINTRP + (3*NCART + 1)/IRAT
      KRSPVC = KINTDR + (3*NCART + 1)/IRAT
      KHDOBR = KRSPVC + NVARPT
      KDB    = KHDOBR + NORBT*NORBT*3*NCART
      KDR    = KDB    + NORBT*NORBT*3
      KS1MAG = KDR    + NORBT*NORBT*NCART
      KSQHDO = KS1MAG + NORBT*NORBT*3
      KKAPPR = KSQHDO + NORBT*NORBT*NCART
      KKAPPB = KKAPPR + NORBT*NORBT*NCART
      KAATOK = KKAPPB + NORBT*NORBT*3
      KAATOD = KAATOK + 3*NCART
      KAATS1 = KAATOD + 3*NCART
      KAATS2 = KAATS1 + 3*NCART
      KAATOT = KAATS2 + 3*NCART
      KAATAL = KAATOT + 3*NCART
      KAATNC = KAATAL + 3*NCART
      KLAST  = KAATNC + 3*NCART
      IF (KLAST .GT. LWORK) CALL STOPIT('AATNOD',' ',LWORK,KLAST)
      LWRK   = LWORK - KLAST + 1
C
      CALL AATNO1(DTOT,CMO,LABINT,WORK(KINTRP),WORK(KINTDR),
     &            WORK(KRSPVC),WORK(KLAST),LWRK,WORK(KHDOBR),WORK(KDB),
     &            WORK(KDR),WORK(KS1MAG),WORK(KSQHDO),WORK(KKAPPR),
     &            WORK(KKAPPB),WORK(KAATOK),WORK(KAATOD),WORK(KAATS1),
     &            WORK(KAATS2),WORK(KAATOT),WORK(KAATAL),CSTRA,SCTRA,
     &            NCART,INTPRI,IPRINT)
C
      RETURN
      END
C  /* Deck aatno1 */
C***********************************************************************
      SUBROUTINE AATNO1(DTOT,CMO,LABINT,INTREP,INTADR,RSPVEC,WORK,LWORK,
     &                  HDOBR,DB,DR,S1MAGL,SQHDOL,RKAPPA,BKAPPA,AATOK,
     &                  AATOD,AATS1,AATS2,AATOT,AATALL,CSTRA,SCTRA,
     &                  NCART,INTPRI,IPRINT)
C***********************************************************************
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <iratdef.h>

#include <inftap.h>
#include <inforb.h>
#include <inflin.h>
#include <abainf.h>
#include <aatens.h>

      PARAMETER (D05 = 0.5D0)
      INTEGER B, R
      CHARACTER*8 LABINT(*)
      DIMENSION DTOT(NORBT,NORBT), HDOBR(NORBT,NORBT,3,NCART),
     &          CMO(*), INTREP(*), DB(NORBT,NORBT,3),
     &          DR(NORBT,NORBT,NCART), S1MAGL(NORBT,NORBT,3),
     &          SQHDOL(NORBT,NORBT,NCART), RSPVEC(NVARPT),
     &          RKAPPA(NORBT,NORBT,NCART), BKAPPA(NORBT,NORBT,3),
     &          AATOK(3,NCART), AATOD(3,NCART), AATS1(3,NCART),
     &          AATS2(3,NCART), AATOT(3,NCART), AATALL(3,NCART),
     &          CSTRA(*), SCTRA(*), WORK(LWORK)
C
C     DB integrals
C     ============
C
      NCOMP  = 0
      NPATOM = 0
      CALL GET1IN(DUMMY,'HBDO   ',NCOMP,WORK,LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &            .FALSE.,DUMMY,INTPRI)
C
      DO 100 B = 1, 3
         KSYMP = -1
         CALL PRPGET(LABINT(B),CMO,DB(1,1,B),KSYMP,ASM,
     &               WORK,LWORK,IPRINT)
         CALL DSCAL(N2ORBX,D05,DB(1,1,B),1)
  100 CONTINUE
      CALL AATNPR(DB,LABINT,3,IPRINT)
C
C     DR integrals
C     ============
C
      NCOMP  = 0
      NPATOM = 0
      CALL GET1IN(DUMMY,'HDO    ',NCOMP,WORK,LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &            .FALSE.,DUMMY,INTPRI)
      DO 200 R = 1, NCART
        KSYMP = -1
        CALL PRPGET(LABINT(R),CMO,DR(1,1,R),KSYMP,ASM,WORK,LWORK,IPRINT)
  200 CONTINUE
      CALL AATNPR(DR,LABINT,NCART,IPRINT)
C
C     HDOBR integrals
C     ===============
C
      NCOMP  = 0
      NPATOM = 0
      CALL GET1IN(DUMMY,'HDOBR  ',NCOMP,WORK,LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &            .FALSE.,DUMMY,INTPRI)
      DO 300 R = 1, NCART
      DO 300 B = 1, 3
         ICOOR = 3*(R - 1) + B
         KSYMP = -1
         CALL PRPGET(LABINT(ICOOR),CMO,HDOBR(1,1,B,R),KSYMP,ASM,
     &               WORK,LWORK,IPRINT)
  300 CONTINUE
      CALL AATNPR(HDOBR,LABINT,3*NCART,IPRINT)
C
C     S1MAGR
C     ======
C
      NCOMP  = 0
      NPATOM = 0
      CALL GET1IN(DUMMY,'S1MAGL ',NCOMP,WORK,LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &            .FALSE.,DUMMY,INTPRI)
      DO 400 B = 1, 3
         KSYMP = -1
         CALL PRPGET(LABINT(B),CMO,S1MAGL(1,1,B),KSYMP,ASM,
     &               WORK,LWORK,IPRINT)
 400  CONTINUE
      CALL AATNPR(S1MAGL,LABINT,3,IPRINT)
C
C     SQHDOL
C     ======
C
      NCOMP  = 0
      NPATOM = 0
      CALL GET1IN(DUMMY,'SQHDOL ',NCOMP,WORK,LWORK,LABINT,
     &            INTREP,INTADR,IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,
     &            .FALSE.,DUMMY,INTPRI)
      DO 500 R = 1, NCART
         KSYMP = -1
         CALL PRPGET(LABINT(R),CMO,SQHDOL(1,1,R),KSYMP,ASM,
     &               WORK,LWORK,IPRINT)
 500  CONTINUE
      CALL AATNPR(SQHDOL,LABINT,NCART,IPRINT)
C
C     Responses
C     =========
C
      DO 600 B = 1, 3
         CALL READDX(LURDI,2*(NCART+B) - 1,IRAT*NVARPT,RSPVEC)
         CALL UPWOP(RSPVEC(NCONST+1),BKAPPA(1,1,B),.FALSE.)
         LABINT(B) = 'BKAPPA  '
  600 CONTINUE
      CALL AATNPR(BKAPPA,LABINT,3,IPRINT)
      DO 610 R = 1, NCART
         CALL READDX(LURDR,2*R-1,IRAT*NVARPT,RSPVEC)
         CALL UPWOP(RSPVEC(NCONST+1),RKAPPA(1,1,R),.TRUE.)
         LABINT(R) = 'RKAPPA  '
  610 CONTINUE
      CALL AATNPR(RKAPPA,LABINT,NCART,IPRINT)
C
      DO 700 B = 1, 3
      DO 700 R = 1, NCART
         TERMO1 = 0.0D0
         TERMO2 = 0.0D0
         TERMS1 = 0.0D0
         TERMS2 = 0.0D0
         DO 800 I = 1, NORBT
         DO 800 J = 1, NORBT
            TERMS1 = TERMS1 - DTOT(I,J)*HDOBR(I,J,B,R)
            DO 900 K = 1, NORBT
               TERMO1 = TERMO1 - DTOT(I,J)*(RKAPPA(J,K,R) - DR(J,K,R))
     &                                    *BKAPPA(K,I,B)
               TERMO2 = TERMO2 - DTOT(I,J)*(RKAPPA(J,K,R) - DR(J,K,R))
     &                                    *DB(K,I,B)
               TERMS2 = TERMS2 - DTOT(I,J)*SQHDOL(I,K,R)*S1MAGL(J,K,B)
  900       CONTINUE
  800    CONTINUE
         IF (NOLOND) THEN
            TERMO2 = 0.0D0
            TERMS1 = 0.0D0
            TERMS2 = 0.0D0
         END IF
         AATOK(B,R)  = TERMO1
         AATOD(B,R)  = TERMO2
         AATS1(B,R)  = TERMS1
         AATS2(B,R)  = TERMS2
         AATOT(B,R)  = TERMO1 + TERMO2 + TERMS1 + TERMS2
         AATALL(B,R) = TERMO1 + TERMO2 + TERMS1 + TERMS2 + AATNUC(B,R)
  700 CONTINUE
C
      CALL HEADER('AATOK in AATNOD',-1)
      CALL FCPRI(AATOK,'AAT',CSTRA,SCTRA)
      CALL HEADER('AATOD in AATNOD',-1)
      CALL FCPRI(AATOD,'AAT',CSTRA,SCTRA)
      CALL HEADER('AATS1 in AATNOD',-1)
      CALL FCPRI(AATS1,'AAT',CSTRA,SCTRA)
      CALL HEADER('AATS2 in AATNOD',-1)
      CALL FCPRI(AATS2,'AAT',CSTRA,SCTRA)
      CALL HEADER('Total electronic AAT in AATNOD',-1)
      CALL FCPRI(AATOT,'AAT',CSTRA,SCTRA)
      CALL HEADER('Total AAT in AATNOD',-1)
      CALL FCPRI(AATALL,'AAT',CSTRA,SCTRA)
C
      RETURN
      END
C  /* Deck aatnpr */
C***********************************************************************
      SUBROUTINE AATNPR(AVEC,LABINT,NTYPE,IPRINT)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
      CHARACTER*8 LABINT(*)
      DIMENSION AVEC(N2ORBX,NTYPE)
#include <inforb.h>
C
      IF (IPRINT .GT. 10) THEN
         DO 100 I = 1, NTYPE
            CALL HEADER(LABINT(I)//' from AATNOD',-1)
            CALL OUTPUT(AVEC(1,I),1,NORBT,1,NORBT,
     &                  NORBT,NORBT,1,LUPRI)
  100    CONTINUE
      END IF
      RETURN
      END
C  /* Deck aatopc */
C***********************************************************************
      SUBROUTINE AATOPC(DTOT,LUAAT,CSTRA,SCTRA,WORK,LWORK,IPRINT)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION WORK(LWORK), DTOT(*)
#include <inforb.h>
#include <infdim.h>
#include <inflin.h>
C
      KDB    = 1
      KDR    = KDB    + N2ORBX
      KRSPVB = KDR    + N2ORBX
      KRSPVR = KRSPVB + NVARMA
      KDENB  = KRSPVR + NVARMA
      KDENR  = KDENB  + N2ASHX
      KREFVC = KDENR  + N2ASHX
      KDK    = KREFVC + NCONRF
      KDTOTK = KDK    + N2ORBX
      KCINDX = KDTOTK + N2ORBX
      KLAST  = KCINDX + LCINDX
      LWORK1 = LWORK  - KLAST
      IF (KLAST .GT. LWORK) CALL STOPIT('ORBAAT',' ',KLAST,LWORK)
      CALL GETCIX(WORK(KCINDX),LSYMRF,LSYMRF,WORK(KLAST),LWORK1,0)
      CALL AATOP1(DTOT,WORK(KDB),WORK(KDR),WORK(KRSPVB),WORK(KRSPVR),
     &            WORK(KDENB),WORK(KDENR),WORK(KREFVC),WORK(KDK),
     &            WORK(KDTOTK),WORK(KCINDX),WORK(KLAST),LWORK1,
     &            LUAAT,CSTRA,SCTRA,IPRINT)
      RETURN
      END
C  /* Deck aatop1 */
C***********************************************************************
      SUBROUTINE AATOP1(DTOT,DB,DR,RSPVCB,RSPVCR,DENB,DENR,REFVEC,DK,
     &                  DTOTDK,XNDXCI,WORK,LWORK,LUAAT,CSTRA,
     &                  SCTRA,IPRINT)
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
#include <dummy.h>
      PARAMETER (DM1 = -1.0D0, D0 = 0.0D0, D1 = 1.0D0, DP1 = 0.1D0,
     &           DP25 = 0.25D0)
      LOGICAL ANTI, USEPSM, NORHO2, TDM, CSFEXP
      INTEGER X, XS
      DIMENSION DTOT(*),DB(N2ORBX),DR(N2ORBX),RSPVCB(*),RSPVCR(*)
      DIMENSION DENB(N2ASHX),DENR(N2ASHX),REFVEC(NCONRF),DK(N2ORBX)
      DIMENSION DTOTDK(N2ORBX),XNDXCI(LCINDX),CSTRA(*),SCTRA(*),
     &          WORK(LWORK)
#include <nuclei.h>
#include <inftap.h>
#include <inforb.h>
#include <inflin.h>
#include <infdim.h>
#include <aatens.h>
#include <infinp.h>
#include <abainf.h>
C
C     Open response files
C     ===================
C
      CALL GPOPEN(LURDR,ABARDR,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
      CALL GPOPEN(LURDI,ABARDI,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
C
      DO 100 ISYM = 1, NSYM
         DO 200 X = 1, 3
         IF (ISYMAX(X,2) .EQ. ISYM - 1) THEN
            XS = IPTAX(X,2)
C
C           Set variables needed for MC calculations
C           ========================================
C
            CALL ABAVAR(ISYM,.FALSE.,IPRINT,WORK,LWORK)
         IF (NVARPT .EQ. 0) GO TO 200
C
C           Read B-field response vector and DB matrix
C           ==========================================
C
            IREC = 2*(3*NUCDEP + XS) - 1
            CALL READDX(LURDI,IREC,IRAT*NVARPT,RSPVCB)
            IF (NODIFC) CALL READDX(LUAAT,IREC,IRAT*N2ORBX,DB)
            IF (IPRINT. GT. 10) THEN
               CALL HEADER('B-field response vector read in AATOP1',-1)
               CALL PRKAP(NWOPPT,RSPVCB(1+NCONST),DP1,LUPRI)
               IF (IPRINT .GT. 25) THEN
                  WRITE (LUPRI,*) 'CSF part of B-field response vector'
                  CALL OUTPUT(RSPVCB,1,NCONST,1,1,NCONST,1,1,LUPRI)
               END IF
               IF (NODIFC) THEN
                  CALL HEADER('DB matrix in AATOP1',-1)
                  CALL OUTPUT(DB,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
            END IF
C
C           Set up for construction of orbital part of AAT
C           ==============================================
C
            IF (NWOPPT.GT.0) THEN
               IF (NODIFC) THEN
                  CALL UPWOP(RSPVCB(NCONST+1),DK,.FALSE.)
               ELSE
                  CALL UPWOP(RSPVCB(NCONST+1),DB,.FALSE.)
               ENDIF
               IF (IPRINT .GT. 10) THEN
                  CALL HEADER('B-field response vector written as '
     &                         //'KappaB matrix in AATOP1',-1)
                  IF (NODIFC) THEN
                     CALL OUTPUT(DK,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                  ELSE
                     CALL OUTPUT(DB,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                  ENDIF
               END IF
            ELSE
               IF (NODIFC) THEN
                  CALL DZERO(DK,N2ORBX)
               ELSE
                  CALL DZERO(DB,N2ORBX)
               END IF
            END IF
            IF (NODIFC) CALL DAXPY(N2ORBX,D1,DK,1,DB,1)
            IF (IPRINT .GT. 15) THEN
               CALL HEADER('KappaB plus DB in AATOP1',-1)
               CALL OUTPUT(DB,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
            IF (IPRINT .GT. 25) THEN
               CALL AROUND('DTOT in AATOP1')
               CALL OUTPUT(DTOT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
            CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
     &                 DTOT,NORBT,
     &                 DB,NORBT,0.D0,
     &                 DTOTDK,NORBT)
            IF (IPRINT .GT. 15) THEN
               CALL HEADER('DTOTDK in AATOP1',-1)
               CALL OUTPUT(DTOTDK,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
C
C           Get transition density matrix
C           =============================
C
            IF (NASHT .GE. 1) CALL GETREF(REFVEC,NCONRF)
            IF (NASHT .EQ. 1) THEN
               DENB(1) = D0
C              as RSPVCB(1) = D0 for NCONF .eq. 1
C
            ELSE IF (NASHT.GT.1) THEN
C
C              Set RHOTTYP and USEPSM
C              RHOTYP = for non-symmetrized density matrix
C              (PV(ij,kl) usablly .ne. PV(ij,lk))
C              USEPSM tells densid if it is allowed to use permutation
C              symmetry for Ms = 0.
C              =======================================================
C
               RHOTYP = 2
               IF (FLAG(66)) THEN
                  USEPSM = .FALSE.
               ELSE
                  USEPSM = .TRUE.
               END IF
               ISPIN1 = 0
               ISPIN2 = 0
               TDM    = .TRUE.
               NORHO2 = .TRUE.
               KFREE  = 1
               LFREE  = LWORK
               CSFEXP = .NOT.FLAG(27)
               CALL DZERO(DENB,N2ASHX)
               CALL DENSID(LSYMRF,LSYMST,NCONRF,NCONST,REFVEC,RSPVCB,
     &                      DENB,DUMMY,RHOTYP,CSFEXP,USEPSM,NORHO2,
     &                      ISPIN1,ISPIN2,TDM,XNDXCI,WORK,KFREE,LFREE)
C              CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
C                          RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,
C                          TDM,XNDXCI,WORK,KFREE,LFREE)
               IF (IPRINT .GT. 10) THEN
                  CALL HEADER('DENB matrix in AATOP1',-1)
                  CALL OUTPUT(DENB,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
               END IF
            ENDIF
C
            DO 300 ICOOR = 1, 3*NUCIND
            ISCOOR = IPTCNT(ICOOR,ISYM-1,1)
            IF (ISCOOR .GT. 0) THEN
C
C              Read geometric response vector and DR matrix
C              ============================================
C
               IREC = 2*ISCOOR - 1
               CALL READDX(LURDR,IREC,IRAT*NVARPT,RSPVCR)
               IF(NODIFC) CALL READDX(LUAAT,ISCOOR,IRAT*N2ORBX,DR)
               IF (IPRINT. GT. 15) THEN
                  CALL HEADER('Geometric response vector read in'
     &                        //' AATOP1',-1)
                  CALL PRKAP(NWOPPT,RSPVCR(1+NCONST),DP1,LUPRI)
                  IF (IPRINT .GT. 25) THEN
                     WRITE (LUPRI,'(/A)')
     &                  ' CSF part of geometric response vector'
                     CALL OUTPUT(RSPVCR,1,NCONST,1,1,NCONST,1,1,LUPRI)
                  END IF
                  IF (NODIFC) THEN
                     CALL HEADER('DR matrix in AATOP1',-1)
                     CALL OUTPUT(DR,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                  END IF
               END IF
C
C              Construction of first orbital part of AAT
C              =========================================
C
               IF (NWOPPT.GT.0) THEN
                  IF (NODIFC) THEN
                     CALL UPWOP(RSPVCR(NCONST+1),DK,.TRUE.)
                  ELSE
                     CALL UPWOP(RSPVCR(NCONST+1),DR,.TRUE.)
                  END IF
                  IF (IPRINT .GT. 15) THEN
                     CALL HEADER('Geometrical response vector written '
     &                            //'as matrix in AATOP1',-1)
                     IF (NODIFC) THEN
                        CALL OUTPUT(DK,1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                              LUPRI)
                     ELSE
                        CALL OUTPUT(DR,1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                              LUPRI)
                     ENDIF
                  END IF
               ELSE
                  IF (NODIFC) THEN
                     CALL DZERO(DK,N2ORBX)
                  ELSE
                     CALL DZERO(DR,N2ORBX)
                  END IF
               END IF
               IF (NODIFC) THEN
                  CALL DAXPY(N2ORBX,DM1,DK,1,DR,1)
                  CALL DSCAL(N2ORBX,DM1,DR,1)
               END IF
               IF (IPRINT .GT. 15) THEN
                  CALL HEADER('KappaR minus DR in AATOP1',-1)
                  CALL OUTPUT(DR,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
               AATORB(XS,ISCOOR) = - DDOT(N2ORBX,DTOTDK,1,DR,1)
               IF (IPRINT .GT. 15) THEN
                  WRITE(LUPRI,'(A,I5,I5,F12.6)')
     &               ' XS, ISCOOR, AATORB(XS,ISCOOR) :',
     &                 XS, ISCOOR, AATORB(XS,ISCOOR)
               END IF
C
C              Construction of CI part of AAT
C              ==============================
C
C
C                 Get transition density matrix
C                 =============================
C
               IF (NASHT .EQ. 1) THEN
                  DENR(1) = D0
C                 as RSPVCB(1) = D0 for NCONF .eq. 1
               ELSE IF (NASHT .GT.1) THEN
                  CALL DZERO(DENR,N2ASHX)
                  CALL DENSID(LSYMRF,LSYMST,NCONRF,NCONST,REFVEC,RSPVCR,
     &                        DENR,DUMMY,RHOTYP,CSFEXP,USEPSM,NORHO2,
     &                        ISPIN1,ISPIN2,TDM,XNDXCI,WORK,KFREE,LFREE)
C                 CALL DENSID(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR,RHO1,RHO2,
C                             RHOTYP,CSFEXP,USEPSM,NORHO2,ISPIN1,ISPIN2,
C                             TDM,XNDXCI,WORK,KFREE,LFREE)
               IF (IPRINT .GT. 10) THEN
                  CALL HEADER('DENR matrix in AATOP1',-1)
                  CALL OUTPUT(DENR,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
               END IF
               END IF
               IF (NASHT.GT.1) THEN
                  CALL PRPACT(DR,DENB,ACTDR,ISYM,IPRINT)
                  CALL PRPACT(DB,DENR,ACTDB,ISYM,IPRINT)
C                 CALL PRPACT(PRPMO,DVT,ONEACT,LSYMOP,IPRINT)
                  AATCI(XS,ISCOOR) = AATCI(XS,ISCOOR) - ACTDR - ACTDB
     &                              - DDOT(NCONST,RSPVCB,1,RSPVCR,1)
                  IF (IPRINT .GE. 5) THEN
                     WRITE(LUPRI,*)
     &               'Individual CI contributions to AATCI(R,B)'
                     WRITE(LUPRI,*) 'R, B     ',ISCOOR,XS
                     WRITE(LUPRI,*) 'CI1-DR   ',-ACTDR
                     WRITE(LUPRI,*) 'CI1-DB   ',-ACTDB
                     WRITE(LUPRI,*) 'CI1 total',-ACTDR-ACTDB
                     WRITE(LUPRI,*) 'CI2      ',
     &                  -DDOT(NCONST,RSPVCB,1,RSPVCR,1)
                     WRITE(LUPRI,*) '+CI1+CI2 ',AATCI(XS,ISCOOR)
                  END IF
               END IF
C
            END IF
  300       CONTINUE
         ENDIF
  200    CONTINUE
  100 CONTINUE
C
C     Print first orbital part and CI part of AAT
C     ===========================================
C
      IF (IPRINT .GT. 0) THEN
         CALL HEADER('Total AATORB in AATOP1',-1)
         CALL FCPRI(AATORB,'AAT',CSTRA,SCTRA)
         CALL HEADER('Total AATCI in AATOP1',-1)
         CALL FCPRI(AATCI,'AAT',CSTRA,SCTRA)
      END IF
      IF (DOWALK) THEN
         CALL GPCLOSE(LURDR,'KEEP')
      ELSE
         CALL GPCLOSE(LURDR,'DELETE')
      END IF
      CALL GPCLOSE(LURDI,'DELETE')
C
      RETURN
      END
C  /* Deck prpact */
C***********************************************************************
      SUBROUTINE PRPACT(PRPMO,DVT,ONEACT,LSYMOP,IPRINT)
C***********************************************************************
C
C CALCULATE  TRANSITION VALUE OF ACTIVE PART OF ONE ELECTRON OPERATOR
C
C INPUT  : PRPMO ( ONE ELECTRON MATRIX )
C          DVT   ( ACTIVE PART OF TRANSION DENSITY )
C
C OUTPUT : ONEACT ( TRANSITION VALUE )
C
#include <implicit.h>
C
#include <priunit.h>
#include <inforb.h>
#include <infdim.h>
C
      DIMENSION PRPMO(NORBT,2),DVT(NASHDI,2)
C
      PARAMETER (  D0 = 0.0D0 )
      PARAMETER ( BIGLIM = 100000.0D0, SMLLIM = 0.01D0 )
C
      ONEACT = D0
      DO 50 ISYM = 1,NSYM
         JSYM = MULD2H(LSYMOP,ISYM)
         NASHI = NASH(ISYM)
         NISHI = NISH(ISYM)
         IORBI = IORB(ISYM)
         IASHI = IASH(ISYM)
         NASHJ = NASH(JSYM)
         NISHJ = NISH(JSYM)
         IORBJ = IORB(JSYM)
         IASHJ = IASH(JSYM)
         DO 60 IACT = 1,NASHI
            DO 70 JACT = 1,NASHJ
               ONEACT = ONEACT + DVT(IASHI+IACT,IASHJ+JACT) *
     *               PRPMO(IORBI+NISHI+IACT,IORBJ+NISHJ+JACT)
 70         CONTINUE
 60      CONTINUE
 50   CONTINUE
C
      IF (IPRINT.GT.10) THEN
         IF (ABS(ONEACT) .GT. SMLLIM .AND. ABS(ONEACT) .LT. BIGLIM)
     *                                                       THEN
            WRITE(LUPRI,'(/5X,A,F15.8)')
     *      ' ACTIVE PART OF ONE ELECTRON PROPERTY :',ONEACT
         ELSE
            WRITE(LUPRI,'(/5X,A,1P,D15.8)')
     *      ' ACTIVE PART OF ONE ELECTRON PROPERTY :',ONEACT
         END IF
      END IF
C
C END OF ACTPRP
C
      RETURN
      END
