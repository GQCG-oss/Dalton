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
C  /* Deck opgctl */
      SUBROUTINE OPGCTL(ORTOPG,RHSOPG,CMO,DV,SOINT,PRPACT,PRPICT,WORK,
     &                  LWORK,IREP,NODC,NODV,IPRINT)
C
C     tuh Dec 1989
C
C     Control routine for calculation of reorthonormalization terms
C     (ORTOPG) and right-hand sides of response equation (RHSOPG) for
C     the molecular gradient of an arbitrary one-electron property.
C     Input are SO integrals SOINT of symmetry IREP and molecular
C     orbitals (CMO) and one-electron active density (DV).
C
#include <implicit.h>
#include <mxcent.h>
#include <nuclei.h>
      DIMENSION ORTOPG(3*NUCDEP), RHSOPG(NVARPT), CMO(NCMOT),
     *          DV(NASHT,NASHT), SOINT(NBAST,NBAST), WORK(LWORK)
      LOGICAL NODC, NODV
#include <inforb.h>
#include <inflin.h>
C
      CALL QENTER('OPGCTL')
      KMOINT = 1
      KFCKMO = KMOINT + NORBT*NORBT
      KWRK   = KFCKMO + NORBT*NORBT
      LWRK   = LWORK - KWRK + 1
      IF (KWRK .GT. LWORK) CALL STOPIT('OPGCTL',' ',KWRK,LWORK)
      CALL OPGCT1(ORTOPG,RHSOPG,CMO,DV,SOINT,WORK(KMOINT),WORK(KFCKMO),
     &            PRPACT,PRPICT,WORK(KWRK),LWRK,IREP,NODC,NODV,
     &            IPRINT)
      CALL QEXIT('OPGCTL')
      RETURN
      END
C  /* Deck opgct1 */
      SUBROUTINE OPGCT1(ORTOPG,RHSOPG,CMO,DV,SOINT,DMOINT,FOCKMO,PRPACT,
     &                  PRPICT,WORK,LWORK,IREP,NODC,NODV,IPRINT)
C
C     December 1989, tuh
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
C
      LOGICAL NODC, NODV
      DIMENSION CMO(NCMOT), DV(NASHT,NASHT),
     *          SOINT(NBAST,NBAST), DMOINT(NORBT,NORBT),
     *          ORTOPG(3*NUCDEP), RHSOPG(NVARPT),
     *          FOCKMO(NORBT,NORBT), WORK(LWORK)
C
#include <abainf.h>
#include <inforb.h>
#include <inflin.h>
C
C     ***********************************
C     ***** Fock matrix in MO basis *****
C     ***********************************
C
      CALL MO1TRA(CMO,DMOINT,SOINT,WORK,LWORK,IREP,IPRINT)
      CALL OPGFCK(DMOINT,FOCKMO,DV,IREP,NODC,NODV,IPRINT)
C
C     *********************************************
C     ***** Reorthonormalization contribution *****
C     *********************************************
C
      IF (DIPDER .OR. QPGRAD) THEN
         CALL OPGORT(CMO,FOCKMO,ORTOPG,WORK,LWORK,NODC,NODV,
     &            IREP,IPRINT)
      END IF
C
C     ***************************
C     ***** Right-hand side *****
C     ***************************
C
      CALL OPGRHS(RHSOPG,FOCKMO,DMOINT,PRPACT,PRPICT,WORK,LWORK,
     &            IPRINT)
C
      RETURN
      END
C  /* Deck opgfck */
      SUBROUTINE OPGFCK(DMOINT,FOCK,DV,IREP,NODC,NODV,IPRINT)
C
C     Purpose: Construction of property Fock matrix in MO basis
C
C     tuh 131289
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER(D2=2.0D0)
C
      INTEGER Q, U, U1, U2
      LOGICAL NODC, NODV
      DIMENSION DMOINT(NORBT,NORBT), FOCK(NORBT,NORBT), DV(NASHT,NASHT)
C
#include <inforb.h>
C
      ISYMPT = IREP + 1
      IF (IPRINT .GE. 5) THEN
         CALL TITLER('Output from OPGFCK','*',103)
         WRITE (LUPRI,'(/A,I5)') ' ISYMPT ', ISYMPT
         IF (IPRINT .GE. 10) THEN
            CALL AROUND('Integrals in OPGFCK')
            CALL OUTPUT(DMOINT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            CALL AROUND('DV matrix in OPGFCK')
            CALL OUTPUT(DV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
      END IF
C
      CALL DZERO(FOCK,N2ORBX)
C
C     Inactive part
C
      IF (.NOT.NODC) THEN
         DO 100 ISYMI = 1, NSYM
            ISYMQ = MULD2H(ISYMI,ISYMPT)
            NISHI = NISH(ISYMI)
            NORBQ = NORB(ISYMQ)
            IF (NISHI.GT.0 .AND. NORBQ.GT.0) THEN
               ISTRI = IORB(ISYMI) + 1
               IOFFQ = IORB(ISYMQ)
               DO 110 Q = IOFFQ + 1, IOFFQ + NORBQ
                  CALL DCOPY(NISHI,DMOINT(ISTRI,Q),1,FOCK(ISTRI,Q),1)
                  CALL DSCAL(NISHI,D2,FOCK(ISTRI,Q),1)
 110           CONTINUE
            END IF
 100     CONTINUE
      END IF
C
C     Active part
C
      IF (.NOT.NODV) THEN
         DO 200 ISYMU = 1, NSYM
            ISYMQ = MULD2H(ISYMU,ISYMPT)
            NASHU = NASH(ISYMU)
            NORBQ = NORB(ISYMQ)
            IF (NASHU.GT.0 .AND. NORBQ.GT.0) THEN
               IOFFU1 = IORB(ISYMU) + NISH(ISYMU)
               IOFFU2 = IASH(ISYMU)
               IOFFQ  = IORB(ISYMQ)
               DO 210 U = 1, NASHU
                  U1 = IOFFU1 + U
                  U2 = IOFFU2 + U
                  DO 220 Q = IOFFQ + 1, IOFFQ + NORBQ
                     FOCK(U1,Q) = DDOT(NASHU,DV (IOFFU2 + 1, U2),1,
     *                                 DMOINT(IOFFU1 + 1, Q ),1)
 220              CONTINUE
 210           CONTINUE
            END IF
 200     CONTINUE
      END IF
C
      IF (IPRINT .GE. 7) THEN
         CALL AROUND('Property Fock matrix in OPGFCK')
         CALL OUTPUT(FOCK,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck opgort */
      SUBROUTINE OPGORT(CMO,FOCKMO,ORTOPG,WORK,LWORK,
     &                  NODC,NODV,IREP,IPRINT)
C
C     Dec 89, tuh
C
C     This subroutine calculates the reorthonormalization contribution
C     to the one-electron property gradient by
C
C      (1)  transforming the property Fock matrices to (contravariant)
C           SO basis
C      (2)  contracting these matrices with the differentiated overlap
C           matrices (SO basis)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
      DIMENSION CMO(NCMOT), FOCKMO(NORBT,NORBT), WORK(LWORK),
     &          ORTOPG(3*NUCDEP)
#include <inforb.h>
C
      CALL QENTER('OPGORT')
      IF (IPRINT .GT. 5) CALL TITLER('Output from OPGORT','*',103)
      KSDSP  = 1
      KSDSQ  = KSDSP  + 3*N2BASX
      KSDSQT = KSDSQ  + N2ORBX
      KWRK   = KSDSQT + N2ORBX
      LWRK   = LWORK  - KWRK + 1
      IF (KWRK .GE. LWORK) CALL STOPIT('OPGORT',' ',KWRK,LWORK)
      CALL OPGOR1(CMO,FOCKMO,ORTOPG,WORK(KSDSP),
     &            WORK(KSDSQ),WORK(KWRK),LWRK,WORK(KSDSQT),
     &            IREP,IPRINT)
      CALL QEXIT('OPGORT')
      RETURN
      END
C  /* Deck opgor1 */
      SUBROUTINE OPGOR1(CMO,FOCKMO,ORTOPG,SDSP,SDSQ,WORK,LWORK,
     &                  SDSQT,IREP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
      PARAMETER (D2 = 2.0D0)
      LOGICAL DERIV(3), DOTD
      DIMENSION CMO(NCMOT), FOCKMO(NORBT,NORBT), ORTOPG(3*NUCDEP),
     *          SDSP(N2BASX,3),
     *          SDSQ(NORBT,NORBT), SDSQT(NORBT,NORBT), WORK(LWORK)
      CHARACTER*4 KEY
#include <abainf.h>
#include <symmet.h>
#include <inforb.h>
C
      DATA DERIV /3*.TRUE./
C
      DOTD = .FALSE.
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('MO Fock matrix in OPGOR1',-1)
         CALL OUTPUT(FOCKMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C     ***** Calculate reorthonormalization terms *****
C
      CALL DZERO(ORTOPG,3*NUCDEP)
      DO 100 IATOM = 1, NUCIND
         IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I5)') ' IATOM ', IATOM
         KEY = 'DMAT'
         IF (NODIFC) KEY = 'OMAT'
         CALL GETSD(DERIV,CMO,SDSP,WORK,LWORK,IATOM,.FALSE.,
     *              KEY,DOTD,IPRINT,IPRINT)
         DO 200 ICOOR = 1, 3
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I5)') ' ICOOR ', ICOOR
            ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
            IF (ISCOOR .GT. 0) THEN
               IF (IPRINT.GT.5) WRITE (LUPRI,'(/A,I5)') ' ISCOOR',ISCOOR
               CALL EXFDTD(SDSP(1,ICOOR),SDSQ,IREP+1,IPRINT)
               CALL TRPMAT(SDSQ,NORBT,NORBT,SDSQT)
               ORTOPG(ISCOOR) = - DDOT(N2ORBX,SDSQT,1,FOCKMO,1)
               IF (.NOT. NODIFC) ORTOPG(ISCOOR) = D2*ORTOPG(ISCOOR)
            END IF
  200    CONTINUE
  100 CONTINUE
C
C     ***** Print *****
C
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('ORTOPG in OPGOR1',-1)
         CALL OUTPUT(ORTOPG,1,1,1,3*NUCDEP,1,3*NUCDEP,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck opgrhs */
      SUBROUTINE OPGRHS(RHSOPG,FOCKMO,DMOINT,PRPACT,PRPICT,WORK,LWORK,
     &                  IPRINT)
C
C     Jan 1990 tuh
C
C     Calculate right-hand side for response equations for
C     one-electron properties
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
      DIMENSION RHSOPG(NVARPT), FOCKMO(NORBT,NORBT),
     &          DMOINT(NORBT,NORBT), WORK(LWORK)
C
#include <inftap.h>
#include <infinp.h>
#include <inforb.h>
#include <infdim.h>
#include <inflin.h>
C
C
      CALL QENTER('OPGRHS')
      CALL TIMER('START ',TIME1,TIME2)
      IF (IPRINT .GT. 4) CALL TITLER('Output from OPGRHS','*',103)
C
C     *******************************************
C     ***** Construction of right-hand side *****
C     *******************************************
C
      CALL DZERO(RHSOPG,NVARPT)
C
C     (A) Construct orbital part of right-hand side
C     ---------------------------------------------
C
      CALL OPGORB(FOCKMO,RHSOPG(NCONST + 1))
C
C     Average rotation gradients if SUPSYM and sym 1:
C
      IF (LSYMPT.EQ.1) CALL AVERAG(RHSOPG(NCONST + 1),NWOPPT,1)
C
C     (B) Construct configuration part of right-hand side
C     ---------------------------------------------------
C
      PRPACT = D0
      IF (NCONST .GT. 1) THEN
C
C        Work space allocation:
C
         KCREF  = 1
         KACAC  = KCREF  + NCONRF
         KCINDX = KACAC  + NASHT*NASHT
         KWRK   = KCINDX + LCINDX
         LWRK   = LWORK  - KWRK + 1
         IF (KWRK .GE. LWORK) CALL STOPIT('OPGRHS','GETCIX',KWRK,LWORK)
C
C        CI index vector:
C
         CALL GETCIX(WORK(KCINDX),LSYMRF,LSYMST,WORK(KWRK),LWRK,0)
C
C        Active part of property integrals:
C
         CALL GETAC1(DMOINT,WORK(KACAC))
C
         IF (NCONRF .GT. 0) THEN
            REWIND LUSIFC
            CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
            READ (LUSIFC)
            READ (LUSIFC)
            READ (LUSIFC)
            CALL READI(LUSIFC,IRAT*NCONRF,WORK(KCREF))
         END IF
C
         CALL CIPRP(1,WORK(KCREF),RHSOPG,NVARPT,WORK(KACAC),
     *              WORK(KCINDX),WORK(KWRK),LWRK)
C        CALL CIPRP(NSIM,CREF,SCVECS,LSCVEC,PRPAC,CINDEX,WORK,LFREE)
C
         IF (LSYMPT .EQ. 1) THEN
            PRPACT = DDOT(NCONST,WORK(KCREF),1,RHSOPG,1)
            CALL DAXPY(NCONST,(-PRPACT),WORK(KCREF),1,RHSOPG,1)
         ELSE
            PRPACT = D0
         END IF
         CALL DSCAL(NCONST,D2,RHSOPG,1)
      END IF
C
C     ***** Print right-hand side *****
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Orbital part of right-hand side',-1)
         IF (IPRINT .GT. 20) THEN
            PRFAC = 0.0D0
         ELSE
            PRFAC = 0.1D0
         END IF
         CALL PRKAP(NWOPPT,RHSOPG(NCONST + 1),PRFAC,LUPRI)
         IF (IPRINT .GT. 20) THEN
            CALL HEADER('Configuration part of right-hand side',-1)
            CALL OUTPUT(RHSOPG,1,1,1,NCONST,1,NCONST,1,LUPRI)
         END IF
      END IF
C
C     One-electron property calculated from CI part of RHS
C
      IF (LSYMPT .EQ. 1) THEN
         PRPICT = D0
         DO 100 ISYM = 1, NSYM
            IORB1 = IORB(ISYM) + 1
            NISHI = NISH(ISYM)
            IF (NISHI .GT. 0) THEN
               PRPICT = PRPICT + DSUM(NISHI,DMOINT(IORB1,IORB1),NORBT+1)
            END IF
  100    CONTINUE
         PRPICT = D2*PRPICT
         PRPCI  = PRPACT + PRPICT
         IF (IPRINT .GT. 5) THEN
            WRITE (LUPRI,'(/,A,D24.12)') ' PRPICT ', PRPICT
            WRITE (LUPRI,'(  A,D24.12)') ' PRPACT ', PRPACT
            WRITE (LUPRI,'(  A,D24.12)') ' PRPCI  ', PRPCI
         END IF
      ELSE
         PRPICT = D0
         PRPACT = D0
      END IF
C
      IF (IPRINT .GT. 1) CALL TIMER('OPGRHS',TIME1,TIME2)
      CALL QEXIT('OPGRHS')
      RETURN
C
C end of OPGRHS
C
      END
C  /* Deck opgorb */
      SUBROUTINE OPGORB(FOCKMO,ORBGRD)
C
C Dec-1989 hjaaj & tuh
C
C Purpose:
C   To add the orbital property gradients in ORBGRD
C
C Input:
C   The total property Fock matrix FOCKMO
C
C Output:
C   Property orbital gradients in ORBGRD
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <infvar.h>
      PARAMETER (D2 = 2.0D0)
      DIMENSION FOCKMO(NORBT,NORBT), ORBGRD(NWOPT)
C
C Used from common blocks:
C   INFVAR : NWOPT,JWOP(2,*)
C   INFORB : NORBT
C
#include <inforb.h>
C
      DO 100 IG = 1, NWOPT
         K = JWOP(1,IG)
         L = JWOP(2,IG)
         ORBGRD(IG) = ORBGRD(IG) + D2*(FOCKMO(K,L) - FOCKMO(L,K))
  100 CONTINUE
      RETURN
      END
C  /* Deck mo1tra */
      SUBROUTINE MO1TRA(CMO,SQMO,SQSO,WORK,LWORK,IREPO,IPRINT)
C
C     This routine transforms derivative SO matrices to MO basis.
C
C     tuh 010190
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <inforb.h>
      DIMENSION CMO(NCMOT), WORK(LWORK), SQMO(NORBT,NORBT), SQSO(N2BASX)
      PARAMETER (D1 =1.0D0)
C
C     ***** Print Section *****
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from MO1TRA',-1)
         WRITE (LUPRI,'(A,I5)') ' IREPO  ', IREPO
         WRITE (LUPRI,'(A,I5)') ' NORBT  ', NORBT
         WRITE (LUPRI,'(A,I5)') ' NBAST  ', NBAST
         WRITE (LUPRI,'(A,I5)') ' NCMOT  ', NCMOT
         WRITE (LUPRI,'(A,I10)') ' LWORK  ', LWORK
         IF (IPRINT .GT. 15) THEN
            CALL HEADER('Square SO matrix',-1)
            CALL OUTPUT(SQSO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
C
C     ***** Transform matrix from SO to MO basis *****
C
      CALL DZERO(SQMO,N2ORBX)
C
C     Loop over irreps for orbitals
C
      ISYMO = IREPO + 1
      DO 100 ISYMA = 1, NSYM
      IF (NORB(ISYMA) .GT. 0) THEN
         ISYMB = MULD2H(ISYMA,ISYMO)
         IF ((ISYMA.GE.ISYMB) .AND. (NORB(ISYMB).GT.0)) THEN
C
C           Print symmetries and orbitals
C
            IF (IPRINT.GT.15) THEN
               WRITE (LUPRI,'(A,3I3)') ' ISYMA/B/O',ISYMA,ISYMB,ISYMO
               WRITE (LUPRI,'(/A,I5,A)')
     *            ' MO coefficients for symmetry', ISYMA
               CALL OUTPUT(CMO(ICMO(ISYMA)+1),1,NBAS(ISYMA),1,
     *            NORB(ISYMA),NBAS(ISYMA),NORB(ISYMA),1,LUPRI)
               WRITE (LUPRI,'(/A,I5,A)')
     *            ' MO coefficients for symmetry', ISYMB
               CALL OUTPUT(CMO(ICMO(ISYMB)+1),1,NBAS(ISYMB),1,
     *            NORB(ISYMB),NBAS(ISYMB),NORB(ISYMB),1,LUPRI)
            END IF
C
C           Transform matrix block(s)
C
            CALL UTHV(CMO(ICMO(ISYMA)+1),SQSO,CMO(ICMO(ISYMB)+1),
     *                ISYMA,ISYMB,NBAS(ISYMA),NBAS(ISYMB),SQMO,WORK)
C
C           Print transformed matrix thus far
C
            IF (IPRINT.GT.25) THEN
               WRITE(LUPRI,'(/4A)')' Unfinished matrix in MO basis'
               CALL OUTPUT(SQMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
 100  CONTINUE
      IF (ISYMO .GT. 1) CALL TRANSX(SQMO,SQMO,NORBT,NORBT,D1,IPRINT)
      IF (IPRINT.GT.15) THEN
         CALL HEADER('Matrix in MO basis',-1)
         CALL OUTPUT(SQMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
