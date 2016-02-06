c -*-mode:fortran; fortran-continuation-string: "&" -*-
!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C  /* Deck getfd */
      SUBROUTINE GETFD(IATOM,NOORBI,TWOH1,NOH1,NOH2,NODC,NODV,IPRINT,
     &                 DERIV,EMYD,CMO,DV,FCD,FVD,FDFTD,WORK,LWORK,AOMAT,
     &                 DOREAL,INTPRI)
C
C  3-Nov-1988 hjaaj -- interface routine for GETFD
C
#include "implicit.h"
      LOGICAL   TWOH1, NOH1, NOH2, NODC, NODV, DERIV(3),
     *          NOORBI, AOMAT, DOREAL
      DIMENSION CMO(*), DV(*), FCD(N2BASX,3), FVD(N2ORBX,3),
     *          EMYD(3), FDFTD(N2BASX,3), WORK(LWORK)
#include "priunit.h"
#include "iratdef.h"
C
C Used from common blocks:
C   INFORB : NNBASX,...
C
#include "inforb.h"
C
      CALL QENTER('GETFD')
      KFCAOX = 1
      KFCAOY = KFCAOX + NNBASX
      KFCAOZ = KFCAOY + NNBASX
      IF (NODV) THEN
         KFVAOX = KFCAOZ + NNBASX
         KFVAOY = KFVAOX
         KFVAOZ = KFVAOX
         KDVAO  = KFVAOX
         KDCAO  = KFVAOX
      ELSE
         KFVAOX = KFCAOZ + NNBASX
         KFVAOY = KFVAOX + NNBASX
         KFVAOZ = KFVAOY + NNBASX
         KDVAO  = KFVAOZ + NNBASX
         KDCAO  = KDVAO  + NNBASX
      END IF
      KDCFO  = KDCAO  + NNBASX
      KIINDX4= KDCFO  + NNBASX
      KGETFD = KIINDX4+ 4*600/IRAT
      LGETFD = LWORK  - KGETFD + 1
      IF (KGETFD.GT.LWORK) CALL STOPIT('GETFD',' ',KGETFD,LWORK)
C
      CALL GETFD1(IATOM,NOORBI,TWOH1,NOH1,NOH2,NODC,NODV,IPRINT,
     &            DERIV,EMYD,CMO,DV,FCD,FVD,FDFTD,
     &            WORK(KFCAOX),WORK(KFCAOY),WORK(KFCAOZ),
     &            WORK(KFVAOX),WORK(KFVAOY),WORK(KFVAOZ),
     &            WORK(KDCAO), WORK(KDVAO), WORK(KDCFO),
     &            WORK(KIINDX4),WORK(KGETFD),LGETFD,AOMAT,DOREAL,
     &            INTPRI)
      CALL QEXIT('GETFD')
      RETURN
      END
C  /* Deck getfd1 */
      SUBROUTINE GETFD1(IATOM,NOORBI,TWOH1,NOH1,NOH2,NODC,NODV,IPRINT,
     &                  DERIV,EMYD,CMO,DV,FCD,FVD,FDFTD,FCSOX,FCSOY,
     &                  FCSOZ,FVSOX,FVSOY,FVSOZ,DCSO,DVSO,DCFOLD,
     &                  IINDX4,WORK,LWORK,AOMAT,DOREAL,INTPRI)
C
C     This subroutine calculates the AO-differentiated Fock (inactive
C     and active) in MO basis.
C
C     tuh Sep 1988
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "iratdef.h"
      PARAMETER (D0 = 0.00 D00, D1 = 1.0D0, DP5 = 0.50D00)
      LOGICAL TWOH1, NOH1, NOH2, NODC, NODV, DERIV(3), DERX, DERY, DERZ,
     *        NOORBI, AOMAT, DOFV, DOREAL
#include "inforb.h"
      DIMENSION IINDX4(4,600), CMO(*), DV(*),
     *          FCD(N2BASX,3), FVD(N2ORBX,3), FDFTD(N2BASX,3),
     *          FCSOX(NNBASX), FCSOY(NNBASX), FCSOZ(NNBASX),
     *          FVSOX(NNBASX), FVSOY(NNBASX), FVSOZ(NNBASX),
     *          DCSO(NNBASX),  DVSO(NNBASX),
     *          EMYD(3), DCFOLD(*), WORK(LWORK)
#include "ccom.h"
#include "abainf.h"
#include "dftcom.h"
#include "nuclei.h"
#include "inftap.h"
#include "eribuf.h"

      DOFV = (NISHT .GT. 0) .AND. (NASHT .GT. 0) .AND. (.NOT.NOORBI)
      DERX = DERIV(1)
      DERY = DERIV(2)
      DERZ = DERIV(3)
      IF (IPRINT .GT. 05 .OR. AOMAT) THEN
         CALL TITLER('Output from GETFD','*',103)
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, '(A,3L5)')' DERIV ', DERX, DERY, DERZ
         WRITE (LUPRI, '(A,L5)') ' NOH1  ', NOH1
         WRITE (LUPRI, '(A,L5)') ' NOH2  ', NOH2
         WRITE (LUPRI, '(A,L5)') ' NODC  ', NODC
         WRITE (LUPRI, '(A,L5)') ' NODV  ', NODV
         WRITE (LUPRI, '(A,L5)') ' AOMAT ', AOMAT
      END IF
      CALL ERIBUF_INI
      LBUF = 600
C
      LENINT4 = 2*LBUF + NIBUF*LBUF + 1  ! length in integer*4 units
C
C
C     ***************************************************
C     ***** Set up active and inactive one-electron *****
C     ***** density matrices in SO basis            *****
C     ***************************************************
C
      CALL DCDVSO(CMO,DV,DCSO,DVSO,DCFOLD,IPRINT,NODC,NODV)
C
C     *******************************************************
C     ***** One-electron contributions to Fock matrices *****
C     *******************************************************
C
C     Inactive FOCK matrix
C
      CALL DZERO(FCD,3*N2BASX)
      IF (NOH1) THEN
         EMY1X = D0
         EMY1Y = D0
         EMY1Z = D0
      ELSE
C
C     Due to changes in GETSD, UFCSOX/Y/Z are now square, K.Ruud, Nov.-93
C     Note that in order to save memory, we use FCD for UFCDSOX/Y/Z
C
         CALL ONEDRL('HMAT',FCD(1,1),FCD(1,2),FCD(1,3),IATOM,DERX,
     &               DERY,DERZ,WORK,LWORK,IPRINT,INTPRI)
C
C        DFT contribution
C
         IF (DFTADD) THEN
            CALL DZERO(FDFTD,3*N2BASX)
            IF (IATOM.GT.0) THEN
               CALL DFTHED(IATOM,FCD,FDFTD,WORK,LWORK,IPRINT)
            ELSE IF (IATOM .EQ. 0 .AND. .NOT.NOLOND) THEN
               CALL DFTBRHS(FCD,WORK,LWORK,IPRINT)
            END IF
         END IF
C
c        IF (DFTADD .AND. IATOM .EQ. 0 .AND. .NOT.NOLOND) THEN
c            CALL DFTBRHS(FCD,WORK,LWORK,IPRINT)
c        END IF
C
         IF (IATOM .GT. 0) THEN ! nuclear derivative - symmetric matrix
            CALL DSITSP(NBAST,FCD(1,1),FCSOX)
            CALL DSITSP(NBAST,FCD(1,2),FCSOY)
            CALL DSITSP(NBAST,FCD(1,3),FCSOZ)
            IF (DERX) EMY1X = DDOT(NNBASX,DCFOLD,1,FCSOX,1)
            IF (DERY) EMY1Y = DDOT(NNBASX,DCFOLD,1,FCSOY,1)
            IF (DERZ) EMY1Z = DDOT(NNBASX,DCFOLD,1,FCSOZ,1)
         ELSE ! magnetic field - antisymmetric matrix
            CALL DGETAP(NBAST,FCD(1,1),FCSOX)
            CALL DGETAP(NBAST,FCD(1,2),FCSOY)
            CALL DGETAP(NBAST,FCD(1,3),FCSOZ)
         END IF
      END IF
C
C     Active Fock matrix
C
      IF (DOFV) THEN
         CALL DZERO(FVSOX,NNBASX)
         CALL DZERO(FVSOY,NNBASX)
         CALL DZERO(FVSOZ,NNBASX)
      END IF
C
C     *******************************************************
C     ***** Two-electron contributions to Fock matrices *****
C     *******************************************************
C
      IF ((NISHT.GT.0) .AND. (.NOT.(NOH2.OR.NOORBI))
     &                 .AND. (.NOT.(IATOM.EQ.0 .AND. NOLOND))) THEN
C
C        Field Derivatives
C        =================
C
         IF (IATOM.EQ.0) THEN
C
            IF (FCKDDR) THEN
               KREAD = 1
               LLAST = KREAD + NNBASX
               IF (LLAST.GT.LWORK) CALL STOPIT('GETFD1','1',LLAST,LWORK)
               IREC = 0
               IF (DOREAL) THEN
                  IREC = 3
                  IF (EXPFCK) IREC = 3*NUCIND
               END IF
               IF (DOFV) IREC = 2*IREC
               LNGTH = IRAT*NNBASX
               CALL READDX(LUDFCK,IREC+1,LNGTH,WORK(KREAD))
               CALL DAXPY(NNBASX,D1,WORK(KREAD),1,FCSOX,1)
               CALL READDX(LUDFCK,IREC+2,LNGTH,WORK(KREAD))
               CALL DAXPY(NNBASX,D1,WORK(KREAD),1,FCSOY,1)
               CALL READDX(LUDFCK,IREC+3,LNGTH,WORK(KREAD))
               CALL DAXPY(NNBASX,D1,WORK(KREAD),1,FCSOZ,1)
               IF (DOFV) THEN
                  CALL READDX(LUDFCK,IREC+4,LNGTH,FVSOX)
                  CALL READDX(LUDFCK,IREC+5,LNGTH,FVSOY)
                  CALL READDX(LUDFCK,IREC+6,LNGTH,FVSOZ)
               END IF
            ELSE
               KFCSQ = 1
               KDCSQ = KFCSQ + 3*NBAST*NBAST
               LLAST = KDCSQ +   NBAST*NBAST
               IF (NASHT .GT. 0) THEN
                  KFVSQ = LLAST
                  KDVSQ = KFVSQ + 3*NBAST*NBAST
                  LLAST = KDVSQ +   NBAST*NBAST
               END IF
               LVLAST = LLAST + 600*(IRAT + 1) + 1
C
               IF (LVLAST.GT.LWORK) 
     &              CALL STOPIT('GETFD1',' ',LVLAST,LWORK)
C
               CALL FCKMG2(FCSOX,FVSOX,DCSO,DVSO,WORK(KFCSQ),
     &                     WORK(KFVSQ),WORK(KDCSQ),WORK(KDVSQ),IINDX4,
     &                     IPRINT,.TRUE.,WORK(LLAST))
            END IF
C
C        Geometrical Derivatives
C        =======================
C
         ELSE
            IF (FCKDDR) THEN
               IREC = 0
               IF (EXPFCK) IREC = 3*(IATOM - 1)
               IF (DOFV)   IREC = 2*IREC
               KREAD = 1
               LLAST = KREAD + NNBASX
               IF (LLAST.GT.LWORK) CALL STOPIT('GETFD1','2',LLAST,LWORK)
               LNGTH = IRAT*NNBASX
               CALL READDX(LUDFCK,IREC+1,LNGTH,WORK(KREAD))
               CALL DAXPY(NNBASX,D1,WORK(KREAD),1,FCSOX,1)
               CALL READDX(LUDFCK,IREC+2,LNGTH,WORK(KREAD))
               CALL DAXPY(NNBASX,D1,WORK(KREAD),1,FCSOY,1)
               CALL READDX(LUDFCK,IREC+3,LNGTH,WORK(KREAD))
               CALL DAXPY(NNBASX,D1,WORK(KREAD),1,FCSOZ,1)
               IF (DOFV) THEN
                  CALL READDX(LUDFCK,IREC+4,LNGTH,FVSOX)
                  CALL READDX(LUDFCK,IREC+5,LNGTH,FVSOY)
                  CALL READDX(LUDFCK,IREC+6,LNGTH,FVSOZ)
               END IF
            ELSE
               CALL REWSPL(LU2DER)
               KINT  = 1
               KIINT = KINT + LBUF
 150           CONTINUE
               CALL READI4(LU2DER,LENINT4,WORK(KINT))
               CALL AOLAB4(WORK(KIINT),LBUF,NIBUF,NBITS,IINDX4,LENGTH)
               IF (IPRINT .GE. 05) WRITE (LUPRI, '(/,A,I5)')
     &               ' Two-electron integral buffer read, LENGTH =',
     &               LENGTH
                  IF (LENGTH .GT. 0) THEN
                     CALL FCKTWO(FCSOX,FVSOX,DCSO,DVSO,WORK(KINT),
     &                           IINDX4,LENGTH,IPRINT)
                  ELSE IF (LENGTH .LT. 0 ) THEN
                     GO TO 300
                  END IF
               GO TO 150
C
  300          CONTINUE
           END IF
        END IF
      END IF
C
C     ********************************************************
C     ***** Transform matrices to MO basis and calculate *****
C     ***** inactive contributions to molecular gradient *****
C     ********************************************************
C
C     x direction
C     ===========
C
      KFSOSQ = 1
      KLAST  = KFSOSQ + N2BASX
      LFREE  = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('GETFD1',' ',KLAST,LWORK)
      IF (DERX) THEN
         IF (IPRINT .GE. 5 .OR. AOMAT) CALL TITLER
     *      ('Differentiation in x direction',' ',200)
         ISIM = 1
         IF (IATOM .GT. 0) THEN
            IF (AOMAT) THEN
               CALL AROUND('Inactive Fock matrix diff. in x direction')
               CALL TRSOAO(FCSOX,WORK,LWORK,NBAST,IATOM,IPRINT)
               IF (DOFV) THEN
                  CALL AROUND('Active Fock matrix diff. in x direction')
                  CALL TRSOAO(FVSOX,WORK,LWORK,NBAST,IATOM,IPRINT)
               END IF
            END IF
            EMYD(ISIM) = DP5*(EMY1X + DDOT(NNBASX,DCFOLD,1,FCSOX,1))
         END IF
         IF (IATOM .EQ. 0) THEN
            CALL DAPTGE(NBAST,FCSOX,WORK(KFSOSQ))
         ELSE
            CALL DSPTSI(NBAST,FCSOX,WORK(KFSOSQ))
         END IF
         CALL ONETRA(CMO,FCD(1,ISIM),WORK(KFSOSQ),WORK(KLAST),LFREE,
     &               IATOM,1,IPRINT)
         IF (DOFV) THEN
            IF (IATOM .GT. 0) THEN
               CALL DSPTSI(NBAST,FVSOX,WORK(KFSOSQ))
            ELSE
               CALL DAPTGE(NBAST,FVSOX,WORK(KFSOSQ))
            END IF
            CALL ONETRA(CMO,FVD(1,ISIM),WORK(KFSOSQ),WORK(KLAST),LFREE,
     &                  IATOM,1,IPRINT)
         END IF
         IF (IPRINT .GE. 5) THEN
            WRITE (LUPRI,'(A,F24.10,/)')
     *        ' Inactive contribution to molecular gradient:',EMYD(ISIM)
            CALL HEADER('Inactive SO Fock matrix',-1)
            CALL OUTPAK(FCSOX,NBAST,1,LUPRI)
            CALL HEADER('Inactive MO Fock matrix',-1)
            CALL OUTPUT(FCD(1,ISIM),1,NORBT,1,NORBT,NBAST,NBAST,1,LUPRI)
            IF (DOFV) THEN
               CALL HEADER('Active SO Fock matrix',-1)
               CALL OUTPAK(FVSOX,NBAST,1,LUPRI)
               CALL HEADER('Active MO Fock matrix',-1)
               CALL OUTPUT(FVD(1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                     LUPRI)
            END IF
         END IF
      END IF
C
C     y direction
C     ===========
C
      IF (DERY) THEN
         IF (IPRINT .GE. 5 .OR. AOMAT) CALL TITLER
     *      ('Differentiation in y direction',' ',200)
         ISIM = 2
         IF (IATOM .GT. 0) THEN
            IF (AOMAT) THEN
               CALL AROUND('Inactive Fock matrix diff. in y direction')
               CALL TRSOAO(FCSOY,WORK,LWORK,NBAST,IATOM,IPRINT)
               IF (DOFV) THEN
                  CALL AROUND('Active Fock matrix diff. in y direction')
                  CALL TRSOAO(FVSOY,WORK,LWORK,NBAST,IATOM,IPRINT)
               END IF
            END IF
            EMYD(ISIM) = DP5*(EMY1Y + DDOT(NNBASX,DCFOLD,1,FCSOY,1))
         END IF
         IF (IATOM .GT. 0) THEN
            CALL DSPTSI(NBAST,FCSOY,WORK(KFSOSQ))
         ELSE
            CALL DAPTGE(NBAST,FCSOY,WORK(KFSOSQ))
         END IF
         CALL ONETRA(CMO,FCD(1,ISIM),WORK(KFSOSQ),WORK(KLAST),LFREE,
     &               IATOM,2,IPRINT)
         IF (DOFV) THEN
            IF (IATOM .GT. 0) THEN
               CALL DSPTSI(NBAST,FVSOY,WORK(KFSOSQ))
            ELSE
               CALL DAPTGE(NBAST,FVSOY,WORK(KFSOSQ))
            END IF
            CALL ONETRA(CMO,FVD(1,ISIM),WORK(KFSOSQ),WORK(KLAST),LFREE,
     &                  IATOM,2,IPRINT)
         END IF
         IF (IPRINT .GE. 5) THEN
            WRITE (LUPRI,'(A,F24.10,/)')
     *        ' Inactive contribution to molecular gradient:',EMYD(ISIM)
            CALL HEADER('Inactive SO Fock matrix',-1)
            CALL OUTPAK(FCSOY,NBAST,1,LUPRI)
            CALL HEADER('Inactive MO Fock matrix',-1)
            CALL OUTPUT(FCD(1,ISIM),1,NORBT,1,NORBT,NBAST,NBAST,1,LUPRI)
            IF (DOFV) THEN
               CALL HEADER('Active SO Fock matrix',-1)
               CALL OUTPAK(FVSOY,NBAST,1,LUPRI)
               CALL HEADER('Active MO Fock matrix',-1)
               CALL OUTPUT(FVD(1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                     LUPRI)
            END IF
         END IF
      END IF
C
C     z direction
C     ===========
C
      IF (DERZ) THEN
         IF (IPRINT .GE. 5 .OR. AOMAT) CALL TITLER
     *      ('Differentiation in z direction',' ',200)
         ISIM = 3
         IF (IATOM .GT. 0) THEN
            IF (AOMAT) THEN
               CALL AROUND('Inactive Fock matrix diff. in z direction')
               CALL TRSOAO(FCSOZ,WORK,LWORK,NBAST,IATOM,IPRINT)
               IF (DOFV) THEN
                  CALL AROUND('Active Fock matrix diff. in z direction')
                  CALL TRSOAO(FVSOZ,WORK,LWORK,NBAST,IATOM,IPRINT)
               END IF
            END IF
            EMYD(ISIM) = DP5*(EMY1Z + DDOT(NNBASX,DCFOLD,1,FCSOZ,1))
         END IF
         IF (IATOM .GT. 0) THEN
            CALL DSPTSI(NBAST,FCSOZ,WORK(KFSOSQ))
         ELSE
            CALL DAPTGE(NBAST,FCSOZ,WORK(KFSOSQ))
         END IF
         CALL ONETRA(CMO,FCD(1,ISIM),WORK(KFSOSQ),WORK(KLAST),LFREE,
     &               IATOM,3,IPRINT)
         IF (DOFV) THEN
            IF (IATOM .EQ. 0) THEN
               CALL DAPTGE(NBAST,FVSOZ,WORK(KFSOSQ))
            ELSE
               CALL DSPTSI(NBAST,FVSOZ,WORK(KFSOSQ))
            END IF
            CALL ONETRA(CMO,FVD(1,ISIM),WORK(KFSOSQ),WORK(KLAST),LFREE,
     &                  IATOM,3,IPRINT)
         END IF
         IF (IPRINT .GE. 5) THEN
            WRITE (LUPRI,'(A,F24.10,/)')
     *        ' Inactive contribution to molecular gradient:',EMYD(ISIM)
            CALL HEADER('Inactive SO Fock matrix',-1)
            CALL OUTPAK(FCSOZ,NBAST,1,LUPRI)
            CALL HEADER('Inactive MO Fock matrix',-1)
            CALL OUTPUT(FCD(1,ISIM),1,NORBT,1,NORBT,NBAST,NBAST,1,LUPRI)
            IF (DOFV) THEN
               CALL HEADER('Active SO Fock matrix',-1)
               CALL OUTPAK(FVSOZ,NBAST,1,LUPRI)
               CALL HEADER('Active MO Fock matrix',-1)
               CALL OUTPUT(FVD(1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                     LUPRI)
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck dcdvso */
      SUBROUTINE DCDVSO(CMO,DV,DCSO,DVSO,DCFOLD,IPRINT,NODC,NODV)
C
C     This subroutine calculates the inactive and active one-electron
C     density matrices in SO basis (contravariant).
C
C     tuh 310888
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0)
      LOGICAL NODC, NODV
      INTEGER R, S, U, V, UV, UR, US, VR, VS, RS
      DIMENSION CMO(*), DV(*), DCSO(*), DVSO(*), DCFOLD(*)
#include "inforb.h"
C
C     ***** Print Section *****
C
      IF (IPRINT .GT. 03) THEN
         WRITE (LUPRI, '(//,A)') ' ----- SUBROUTINE DCDVSO  ----- '
         WRITE (LUPRI, '(A,L5)') ' NODC ', NODC
         WRITE (LUPRI, '(A,L5)') ' NODV ', NODV
         WRITE (LUPRI, '(A,8I5)') ' NISH ', (NISH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NASH ', (NASH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NOCC ', (NOCC(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NORB ', (NORB(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NBAS ', (NBAS(I),I = 1,NSYM)
         IF (IPRINT .GE. 05) THEN
            CALL HEADER('Occupied molecular orbitals',0)
            IEND = 0
            DO 1000 ISYM = 1,NSYM
               IF (NBAS(ISYM) .EQ. 0) GOTO 1000
               IF (NOCC(ISYM) .EQ. 0) GOTO 1100
               WRITE (LUPRI, '(//,A,I5,/)') ' Symmetry ', ISYM
               IENDI = 0
               DO 1200 I = 1, NOCC(ISYM)
                  WRITE (LUPRI, '(/,A,I5,/)')
     *                     ' Molecular orbital ', I
                  WRITE (LUPRI, '(6F12.6)')
     *               (CMO(IEND+IENDI+J), J = 1, NBAS(ISYM))
                  IENDI = IENDI + NBAS(ISYM)
1200           CONTINUE
1100           CONTINUE
               IEND = IEND + NORB(ISYM)*NBAS(ISYM)
1000        CONTINUE
            CALL HEADER('Active density matrix (MO basis)',-1)
            CALL OUTPAK(DV,NASHT,1,LUPRI)
         END IF
      END IF
C
C     ***** Construct contravariant SO matrices *****
C
      CALL DZERO(DCSO,NNBASX)
      CALL DZERO(DCFOLD,NNBASX)
      CALL DZERO(DVSO,NNBASX)
      ISEND  = 0
      ICEND  = 0
      DO 110 ISYM = 1,NSYM
         IOFF = IBAS(ISYM)
         DO 100 R = 1, NBAS(ISYM)
            DO 200 S = 1,R
C
C              (I) Inactive contribution DCRS
C
               DCRS = D0
               ICENDI = 0
               DO 300 I = 1, NISH(ISYM)
                  DCRS = DCRS + CMO(ICEND+ICENDI+R)*CMO(ICEND+ICENDI+S)
                  ICENDI = ICENDI + NBAS(ISYM)
  300          CONTINUE
               DCRS = DCRS + DCRS
               IF (NODC) DCRS = D0
C
C              (II) Active contribution DVRS
C
               DVRS = D0
               IF (.NOT. NODV) THEN
                  IASHI = IASH(ISYM)
                  UV = ((IASHI + 1)*(IASHI + 2))/2
                  IDVEND = ICEND + NBAS(ISYM)*NISH(ISYM)
                  ICENDU = IDVEND
                  DO 400 U = 1,NASH(ISYM)
                     ICENDV = IDVEND
                     DO 410 V = 1, U
                        DUV = DV(UV)
                        IF (ABS(DUV) .GT. D0) THEN
                           TMP = CMO(ICENDU+R)*CMO(ICENDV+S)
                           IF (U .NE. V) TMP = TMP
     *                         + CMO(ICENDU+S)*CMO(ICENDV+R)
                           DVRS = DVRS + DUV*TMP
                        END IF
                        UV = UV + 1
                        ICENDV = ICENDV + NBAS(ISYM)
  410                CONTINUE
                     UV = UV + IASHI
                     ICENDU = ICENDU + NBAS(ISYM)
  400             CONTINUE
               END IF
               RS = (IOFF + R)*(IOFF + R - 1)/2 + IOFF + S
               DCSO(RS) = DCRS
               DVSO(RS) = DVRS
               IF (R .EQ. S) THEN
                  DCFOLD(RS) = DCRS
               ELSE
                  DCFOLD(RS) = DCRS + DCRS
               END IF
  200       CONTINUE
  100    CONTINUE
         ISEND  = ISEND  + (NBAS(ISYM)*(NBAS(ISYM) + 1))/2
         ICEND  = ICEND  + NORB(ISYM)*NBAS(ISYM)
110   CONTINUE
C
C     ***** Print Section *****
C
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('Inactive SO density matrix',-1)
         CALL OUTPAK(DCSO,NBAST,1,LUPRI)
         CALL HEADER('Inactive SO density matrix (folded)',-1)
         CALL OUTPAK(DCFOLD,NBAST,1,LUPRI)
         IF (.NOT. NODV) THEN
            CALL HEADER('Active SO density matrix',-1)
            CALL OUTPAK(DVSO,NBAST,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck fcktwo */
      SUBROUTINE FCKTWO(FC,FV,DC,DV,BUF,IINDX4,LENGTH,IPRINT)
C
C     This subroutine adds derivative two-electron integrals to the
C     active and inactive Fock matrices.
C
C     The 14 classes of integral labels have been taken from Duke,
C     Chem. Phys. Lett. 13 (1972) 76
C
C     Note: Canonical ordering of indices assumed
C
C     tuh 120988
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (TWO = 2.00 D00, THRHLF = 1.50 D00,
     *           HALF = 0.50D00, IBIT08 = 255)
      INTEGER P, Q, R, S
      LOGICAL FIRST, ACTIVE
#include "inforb.h"
      DIMENSION BUF(600), IINDX4(4,600), ITRI(MXCORB),
     *          FC(3*NNBASX), FV(3*NNBASX),
     *          DC(NNBASX),   DV(NNBASX)
      SAVE FIRST, ITRI, ACTIVE, IOFF
      DATA FIRST /.TRUE./

C
      IF (FIRST) THEN
         DO 100 I = 1, NBAST
            ITRI(I) = I*(I - 1)/2
  100    CONTINUE
         ACTIVE = NASHT .GT. 0
         FIRST = .FALSE.
      END IF
C
      DO 200 I = 1, LENGTH
         S    = IINDX4(4,I)
         IF (S .EQ. 0) THEN
            IOFF = (IINDX4(3,I) - 1)*NNBASX
            IF (IPRINT .GE. 25) THEN
               ICOOR = IINDX4(3,I)
               IREPE = IINDX4(2,I)
               WRITE (LUPRI,'(A,2I5)') ' ICOOR, IREPE', ICOOR,IREPE
            END IF
         ELSE
            DINT = BUF(I)
            P = IINDX4(1,I)
            Q = IINDX4(2,I)
            R = IINDX4(3,I)
            IF (IPRINT .GE. 25) WRITE (LUPRI,'(1X,F12.6,5X,4I5)')
     *                                  DINT,P,Q,R,S
            IF (Q .LT. R) THEN
               IF (Q .LT. S) THEN
                  IF (R .EQ. S) THEN
C
C                    p > r = s > q
C
                     II = ITRI(R) + R
                     IK = ITRI(P) + R
                     IL = ITRI(R) + Q
                     KL = ITRI(P) + Q
                     FC(II+IOFF) = FC(II+IOFF) +  TWO*DC(KL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) - HALF*DC(IK)*DINT
                     FC(IK+IOFF) = FC(IK+IOFF) - HALF*DC(IL)*DINT
                     FC(KL+IOFF) = FC(KL+IOFF) +      DC(II)*DINT
                     IF (ACTIVE) THEN
                        FV(II+IOFF) = FV(II+IOFF) +  TWO*DV(KL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) - HALF*DV(IK)*DINT
                        FV(IK+IOFF) = FV(IK+IOFF) - HALF*DV(IL)*DINT
                        FV(KL+IOFF) = FV(KL+IOFF) +      DV(II)*DINT
                     END IF
                  ELSE
C
C                    p > r > s > q
C
                     IJ = ITRI(P) + Q
                     IK = ITRI(P) + R
                     IL = ITRI(P) + S
                     JK = ITRI(R) + Q
                     JL = ITRI(S) + Q
                     KL = ITRI(R) + S
                     FC(IJ+IOFF) = FC(IJ+IOFF) +  TWO*DC(KL)*DINT
                     FC(IK+IOFF) = FC(IK+IOFF) - HALF*DC(JL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) - HALF*DC(JK)*DINT
                     FC(JK+IOFF) = FC(JK+IOFF) - HALF*DC(IL)*DINT
                     FC(JL+IOFF) = FC(JL+IOFF) - HALF*DC(IK)*DINT
                     FC(KL+IOFF) = FC(KL+IOFF) +  TWO*DC(IJ)*DINT
                     IF (ACTIVE) THEN
                        FV(IJ+IOFF) = FV(IJ+IOFF) +  TWO*DV(KL)*DINT
                        FV(IK+IOFF) = FV(IK+IOFF) - HALF*DV(JL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) - HALF*DV(JK)*DINT
                        FV(JK+IOFF) = FV(JK+IOFF) - HALF*DV(IL)*DINT
                        FV(JL+IOFF) = FV(JL+IOFF) - HALF*DV(IK)*DINT
                        FV(KL+IOFF) = FV(KL+IOFF) +  TWO*DV(IJ)*DINT
                     END IF
                  END IF
               ELSE IF (Q .EQ. S) THEN
                  IF (P .EQ. R) THEN
C
C                    p = r > q = s
C
                     II = ITRI(P) + P
                     IJ = ITRI(P) + Q
                     JJ = ITRI(Q) + Q
                     FC(II+IOFF) = FC(II+IOFF) -   HALF*DC(JJ)*DINT
                     FC(IJ+IOFF) = FC(IJ+IOFF) + THRHLF*DC(IJ)*DINT
                     FC(JJ+IOFF) = FC(JJ+IOFF) -   HALF*DC(II)*DINT
                     IF (ACTIVE) THEN
                        FV(II+IOFF) = FV(II+IOFF) -   HALF*DV(JJ)*DINT
                        FV(IJ+IOFF) = FV(IJ+IOFF) + THRHLF*DV(IJ)*DINT
                        FV(JJ+IOFF) = FV(JJ+IOFF) -   HALF*DV(II)*DINT
                     END IF
                  ELSE
C
C                    p > r > q = s
C
                     IJ = ITRI(P) + Q
                     IL = ITRI(P) + R
                     JJ = ITRI(Q) + Q
                     JL = ITRI(R) + Q
                     FC(IJ+IOFF) = FC(IJ+IOFF) + THRHLF*DC(JL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) -   HALF*DC(JJ)*DINT
                     FC(JJ+IOFF) = FC(JJ+IOFF) -        DC(IL)*DINT
                     FC(JL+IOFF) = FC(JL+IOFF) + THRHLF*DC(IJ)*DINT
                     IF (ACTIVE) THEN
                        FV(IJ+IOFF) = FV(IJ+IOFF) + THRHLF*DV(JL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) -   HALF*DV(JJ)*DINT
                        FV(JJ+IOFF) = FV(JJ+IOFF) -        DV(IL)*DINT
                        FV(JL+IOFF) = FV(JL+IOFF) + THRHLF*DV(IJ)*DINT
                     END IF
                  END IF
               ELSE
                  IF (P .EQ. R) THEN
C
C                    p = r > q > s
C
                     IJ = ITRI(P) + Q
                     IL = ITRI(Q) + S
                     JJ = ITRI(P) + P
                     JL = ITRI(P) + S
                     FC(IJ+IOFF) = FC(IJ+IOFF) + THRHLF*DC(JL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) -   HALF*DC(JJ)*DINT
                     FC(JJ+IOFF) = FC(JJ+IOFF) -        DC(IL)*DINT
                     FC(JL+IOFF) = FC(JL+IOFF) + THRHLF*DC(IJ)*DINT
                     IF (ACTIVE) THEN
                        FV(IJ+IOFF) = FV(IJ+IOFF) + THRHLF*DV(JL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) -   HALF*DV(JJ)*DINT
                        FV(JJ+IOFF) = FV(JJ+IOFF) -        DV(IL)*DINT
                        FV(JL+IOFF) = FV(JL+IOFF) + THRHLF*DV(IJ)*DINT
                     END IF
                  ELSE
C
C                    p > r > q > s
C
                     IJ = ITRI(P) + Q
                     IK = ITRI(P) + R
                     IL = ITRI(P) + S
                     JK = ITRI(R) + Q
                     JL = ITRI(Q) + S
                     KL = ITRI(R) + S
                     FC(IJ+IOFF) = FC(IJ+IOFF) +  TWO*DC(KL)*DINT
                     FC(IK+IOFF) = FC(IK+IOFF) - HALF*DC(JL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) - HALF*DC(JK)*DINT
                     FC(JK+IOFF) = FC(JK+IOFF) - HALF*DC(IL)*DINT
                     FC(JL+IOFF) = FC(JL+IOFF) - HALF*DC(IK)*DINT
                     FC(KL+IOFF) = FC(KL+IOFF) +  TWO*DC(IJ)*DINT
                     IF (ACTIVE) THEN
                        FV(IJ+IOFF) = FV(IJ+IOFF) +  TWO*DV(KL)*DINT
                        FV(IK+IOFF) = FV(IK+IOFF) - HALF*DV(JL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) - HALF*DV(JK)*DINT
                        FV(JK+IOFF) = FV(JK+IOFF) - HALF*DV(IL)*DINT
                        FV(JL+IOFF) = FV(JL+IOFF) - HALF*DV(IK)*DINT
                        FV(KL+IOFF) = FV(KL+IOFF) +  TWO*DV(IJ)*DINT
                     END IF
                  END IF
               END IF
            ELSE IF (Q .EQ. R) THEN
               IF (P .EQ. Q) THEN
                  IF (R .EQ. S) THEN
C
C                    p = q = r = s
C
                     II = ITRI(P) + P
                     FC(II+IOFF) = FC(II+IOFF) + HALF*DC(II)*DINT
                     FV(II+IOFF) = FV(II+IOFF) + HALF*DV(II)*DINT
                  ELSE
C
C                    p = q = r > s
C
                     II = ITRI(P) + P
                     IL = ITRI(P) + S
                     FC(II+IOFF) = FC(II+IOFF) +      DC(IL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) + HALF*DC(II)*DINT
                     IF (ACTIVE) THEN
                        FV(II+IOFF) = FV(II+IOFF) +      DV(IL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) + HALF*DV(II)*DINT
                     END IF
                  END IF
               ELSE
                  IF (R .EQ. S) THEN
C
C                    p > q = r = s
C
                     II = ITRI(Q) + Q
                     IL = ITRI(P) + Q
                     FC(II+IOFF) = FC(II+IOFF) +      DC(IL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) + HALF*DC(II)*DINT
                     IF (ACTIVE) THEN
                        FV(II+IOFF) = FV(II+IOFF) +      DV(IL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) + HALF*DV(II)*DINT
                     END IF
                  ELSE
C
C                    p > q = r > s
C
                     IJ = ITRI(P) + Q
                     IL = ITRI(P) + S
                     JJ = ITRI(Q) + Q
                     JL = ITRI(Q) + S
                     FC(IJ+IOFF) = FC(IJ+IOFF) + THRHLF*DC(JL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) -   HALF*DC(JJ)*DINT
                     FC(JJ+IOFF) = FC(JJ+IOFF) -        DC(IL)*DINT
                     FC(JL+IOFF) = FC(JL+IOFF) + THRHLF*DC(IJ)*DINT
                     IF (ACTIVE) THEN
                        FV(IJ+IOFF) = FV(IJ+IOFF) + THRHLF*DV(JL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) -   HALF*DV(JJ)*DINT
                        FV(JJ+IOFF) = FV(JJ+IOFF) -        DV(IL)*DINT
                        FV(JL+IOFF) = FV(JL+IOFF) + THRHLF*DV(IJ)*DINT
                     END IF
                  END IF
               END IF
            ELSE
               IF (P .EQ. Q) THEN
                  IF (R .EQ. S) THEN
C
C                    p = q > r = s
C
                     II = ITRI(P) + P
                     IK = ITRI(P) + R
                     KK = ITRI(R) + R
                     FC(II+IOFF) = FC(II+IOFF) +      DC(KK)*DINT
                     FC(IK+IOFF) = FC(IK+IOFF) - HALF*DC(IK)*DINT
                     FC(KK+IOFF) = FC(KK+IOFF) +      DC(II)*DINT
                     IF (ACTIVE) THEN
                        FV(II+IOFF) = FV(II+IOFF) +      DV(KK)*DINT
                        FV(IK+IOFF) = FV(IK+IOFF) - HALF*DV(IK)*DINT
                        FV(KK+IOFF) = FV(KK+IOFF) +      DV(II)*DINT
                     END IF
                  ELSE
C
C                    p = q > r > s
C
                     II = ITRI(P) + P
                     IK = ITRI(P) + R
                     IL = ITRI(P) + S
                     KL = ITRI(R) + S
                     FC(II+IOFF) = FC(II+IOFF) +  TWO*DC(KL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) - HALF*DC(IK)*DINT
                     FC(IK+IOFF) = FC(IK+IOFF) - HALF*DC(IL)*DINT
                     FC(KL+IOFF) = FC(KL+IOFF) +      DC(II)*DINT
                     IF (ACTIVE) THEN
                        FV(II+IOFF) = FV(II+IOFF) +  TWO*DV(KL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) - HALF*DV(IK)*DINT
                        FV(IK+IOFF) = FV(IK+IOFF) - HALF*DV(IL)*DINT
                        FV(KL+IOFF) = FV(KL+IOFF) +      DV(II)*DINT
                     END IF
                  END IF
               ELSE
                  IF (R .EQ. S) THEN
C
C                    p > q > r = s
C
                     II = ITRI(R) + R
                     IK = ITRI(P) + R
                     IL = ITRI(Q) + R
                     KL = ITRI(P) + Q
                     FC(II+IOFF) = FC(II+IOFF) +  TWO*DC(KL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) - HALF*DC(IK)*DINT
                     FC(IK+IOFF) = FC(IK+IOFF) - HALF*DC(IL)*DINT
                     FC(KL+IOFF) = FC(KL+IOFF) +      DC(II)*DINT
                     IF (ACTIVE) THEN
                        FV(II+IOFF) = FV(II+IOFF) +  TWO*DV(KL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) - HALF*DV(IK)*DINT
                        FV(IK+IOFF) = FV(IK+IOFF) - HALF*DV(IL)*DINT
                        FV(KL+IOFF) = FV(KL+IOFF) +      DV(II)*DINT
                     END IF
                  ELSE
C
C                    p > q > r > s
C
                     IJ = ITRI(P) + Q
                     IK = ITRI(P) + R
                     IL = ITRI(P) + S
                     JK = ITRI(Q) + R
                     JL = ITRI(Q) + S
                     KL = ITRI(R) + S
                     FC(IJ+IOFF) = FC(IJ+IOFF) +  TWO*DC(KL)*DINT
                     FC(IK+IOFF) = FC(IK+IOFF) - HALF*DC(JL)*DINT
                     FC(IL+IOFF) = FC(IL+IOFF) - HALF*DC(JK)*DINT
                     FC(JK+IOFF) = FC(JK+IOFF) - HALF*DC(IL)*DINT
                     FC(JL+IOFF) = FC(JL+IOFF) - HALF*DC(IK)*DINT
                     FC(KL+IOFF) = FC(KL+IOFF) +  TWO*DC(IJ)*DINT
                     IF (ACTIVE) THEN
                        FV(IJ+IOFF) = FV(IJ+IOFF) +  TWO*DV(KL)*DINT
                        FV(IK+IOFF) = FV(IK+IOFF) - HALF*DV(JL)*DINT
                        FV(IL+IOFF) = FV(IL+IOFF) - HALF*DV(JK)*DINT
                        FV(JK+IOFF) = FV(JK+IOFF) - HALF*DV(IL)*DINT
                        FV(JL+IOFF) = FV(JL+IOFF) - HALF*DV(IK)*DINT
                        FV(KL+IOFF) = FV(KL+IOFF) +  TWO*DV(IJ)*DINT
                     END IF
                  END IF
               END IF
            END IF
         END IF
  200 CONTINUE
      RETURN
      END
C  /* Deck fckmg2 */
      SUBROUTINE FCKMG2(FC,FV,DC,DV,FCSQ,FVSQ,DCSQ,DVSQ,IINDX4,
     &                 IPRINT,SYMMET,WORK)
C
C     ================================================================
C     | This subroutine adds magnetic-field derivative two-electron  |
C     | integrals to the active and inactive Fock matrices.          |
C     ================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (DP5 = 0.50D00, IBIT08 = 255)
      LOGICAL SYMMET
      INTEGER P, Q, R, S, X
      DIMENSION FC(3*NNBASX), FV(3*NNBASX), DC(NNBASX), DV(NNBASX),
     &          FCSQ(NBAST,NBAST,3), FVSQ(NBAST,NBAST,3),
     &          DCSQ(NBAST,NBAST), DVSQ(NBAST,NBAST),
     &          IINDX4(4,600), WORK(*)
#include "dftcom.h"
#include "nuclei.h"
#include "inftap.h"
#include "eribuf.h"
#include "inforb.h"

C
      LBUF = 600
      CALL ERIBUF_INI
C
      LENINT4 = 2*LBUF + NIBUF*LBUF + 1  ! length in integer*4 units
      KINT   = 1
      KIINT  = KINT + LBUF
C
      CALL DZERO(FCSQ,3*N2BASX)
      IF (NASHT .GE. 1) CALL DZERO(FVSQ,3*N2BASX)
C
C     Form full (squared) active and inactive density matrices
C     ========================================================
C
      IF (SYMMET) THEN
         CALL DSPTSI(NBAST,DC,DCSQ)
         IF (NASHT .GT. 0) THEN
            CALL DSPTSI(NBAST,DV,DVSQ)
         END IF
      END IF
C
C     Read integrals into buffers and process the integrals
C     =====================================================
C
      CALL REWSPL(LU2DER)
      CALL MOLLAB('AO2MGINT',LU2DER,LUPRI)
C
 300  CONTINUE
         CALL READI4(LU2DER,LENINT4,WORK(KINT))
         CALL AOLAB4(WORK(KIINT),LBUF,NIBUF,NBITS,IINDX4,LENGTH)
         IF (IPRINT .GE. 05) WRITE (LUPRI, '(/,A,I5)')
     *      ' Two-electron integral buffer read, LENGTH =', LENGTH
         IF (LENGTH .GT. 0) THEN
C
            DO 400 I = 1, LENGTH
               S     = IINDX4(4,I)
               IF (S .EQ. 0) THEN
                  X = IINDX4(3,I)
                  IF (IPRINT .GE. 25) THEN
                     IREPE = IINDX4(2,I)
                     WRITE (LUPRI,'(A,2I5)') ' ICOOR, IREPE', X,IREPE
                  END IF
               ELSE
                  CINT = WORK(KINT + I - 1)
                  P = IINDX4(1,I)
                  Q = IINDX4(2,I)
                  R = IINDX4(3,I)
C                 IF (IPRINT .GE. 25) WRITE (LUPRI,'(1X,F12.6,5X,4I5)')
C    &                                        CINT,P,Q,R,S
C
                  IF (P.EQ.R .AND. S.EQ.Q) CINT = DP5*CINT
                  EINT = DP5*HFXFAC*CINT
C
C     Add Coulomb contributions to inactive fock matrix
C     =================================================
C
                  FCSQ(P,Q,X) = FCSQ(P,Q,X) - CINT*DCSQ(R,S)
                  FCSQ(Q,P,X) = FCSQ(Q,P,X) + CINT*DCSQ(S,R)
                  FCSQ(R,S,X) = FCSQ(R,S,X) - CINT*DCSQ(P,Q)
                  FCSQ(S,R,X) = FCSQ(S,R,X) + CINT*DCSQ(Q,P)
C
C     Add Exchange contributions to inactive fock matrix
C     ==================================================
C
                  FCSQ(P,S,X) = FCSQ(P,S,X) + EINT*DCSQ(Q,R)
                  FCSQ(S,P,X) = FCSQ(S,P,X) - EINT*DCSQ(R,Q)
                  FCSQ(R,Q,X) = FCSQ(R,Q,X) + EINT*DCSQ(S,P)
                  FCSQ(Q,R,X) = FCSQ(Q,R,X) - EINT*DCSQ(P,S)
C
                  IF (NASHT.GE.1) THEN
C
C     Add Coulomb contributions to active fock matrix
C     ===============================================
C
                     FVSQ(P,Q,X) = FVSQ(P,Q,X) - CINT*DVSQ(R,S)
                     FVSQ(Q,P,X) = FVSQ(Q,P,X) + CINT*DVSQ(S,R)
                     FVSQ(R,S,X) = FVSQ(R,S,X) - CINT*DVSQ(P,Q)
                     FVSQ(S,R,X) = FVSQ(S,R,X) + CINT*DVSQ(Q,P)
C
C     Add Exchange contributions to active fock matrix
C     ================================================
C
                     FVSQ(P,S,X) = FVSQ(P,S,X) + EINT*DVSQ(Q,R)
                     FVSQ(S,P,X) = FVSQ(S,P,X) - EINT*DVSQ(R,Q)
                     FVSQ(R,Q,X) = FVSQ(R,Q,X) + EINT*DVSQ(S,P)
                     FVSQ(Q,R,X) = FVSQ(Q,R,X) - EINT*DVSQ(P,S)
C
                  END IF
C
               END IF
  400       CONTINUE
         ELSE IF (LENGTH .LT. 0 ) THEN
            GO TO 500
         END IF
      GO TO 300
C
C     Pack active and inactive fock matrices in lower triangular form
C     ===============================================================
C
  500 CONTINUE
C
      IF (SYMMET) THEN
         IPCK = 0
         IF (NASHT .GT. 0) THEN
            DO 800 X = 1,3
               DO 700 J = 1, NBAST
                  DO 600 I = 1, J
                     IPCK     = IPCK + 1
                     FC(IPCK) = FC(IPCK) + FCSQ(I,J,X)
                     FV(IPCK) = FV(IPCK) + FVSQ(I,J,X)
 600              CONTINUE
 700           CONTINUE
 800        CONTINUE
         ELSE
            DO 810 X = 1, 3
               DO 710 J = 1, NBAST
                  DO 610 I = 1, J
                     IPCK     = IPCK + 1
                     FC(IPCK) = FC(IPCK) + FCSQ(I,J,X)
 610              CONTINUE
 710           CONTINUE
 810        CONTINUE
         END IF
      END IF
C
      RETURN
      END
C  /* Deck trsoao */
      SUBROUTINE TRSOAO(FCSO,WORK,LWORK,IDIM,IATOM,IPRINT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION FCSO(*), WORK(LWORK)
C
      IF (IPRINT .GT. 2) CALL TITLER('Subroutine TRSOAO','*',103)
C
C     Allocation
C
      KFCAO = 1
      KTRMT = KFCAO + IDIM*(IDIM + 1)/2
      KWRK  = KTRMT + IDIM*IDIM
      LMAX  = KWRK  + IDIM
      IF (LMAX.GT.LWORK) CALL STOPIT('TRSOAO',' ',LMAX,LWORK)
C
C     Set up AO to SO transformation matrix
C
      CALL TRAORB(WORK(KTRMT),IDIM,IATOM,IPRINT)
C     
C     Transform one-electron matrix
C
      CALL DZERO(WORK(KFCAO),IDIM*(IDIM + 1)/2)
      CALL UTHU(FCSO,WORK(KFCAO),WORK(KTRMT),
     *          WORK(KWRK),IDIM,IDIM)
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Matrix transformed from SO to AO basis',-1)
         CALL OUTPAK(WORK(KFCAO),IDIM,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck traoso */
      SUBROUTINE TRAOSO(FCAO,WORK,LWORK,IDIM,IATOM,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION FCAO(*), WORK(LWORK)
#include "aosotr.h"
C
      IF (IPRINT .GT. 2) CALL TITLER('Subroutine TRAOSO','*',103)
C
C     Allocation
C
      KFCSO = 1
      KTRMT = KFCSO + IDIM*(IDIM + 1)/2
      KWRK  = KTRMT + IDIM*IDIM
      LMAX  = KWRK  + IDIM
      IF (LMAX.GT.LWORK) CALL STOPIT('TRAOSO',' ',LMAX,LWORK)
C
C     Construct transformation matrix
C
      CALL TRAORB(WORK(KTRMT),IDIM,IATOM,IPRINT)
      CALL DGETRN(WORK(KTRMT),IDIM,IDIM)
C
C     Transform one-electron matrix
C
      CALL DZERO(WORK(KFCSO),IDIM*(IDIM + 1)/2)
      CALL UTHU(FCAO,WORK(KFCSO),WORK(KTRMT),
     *          WORK(KWRK),IDIM,IDIM)
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Matrix transformed from AO to SO basis',-1)
         CALL OUTPAK(WORK(KFCSO),IDIM,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck traorb */
      SUBROUTINE TRAORB(ORBTRA,IDIM,IATOM,IPRINT)
C
C     Sets up transformation matrices between SO and AO orbitals
C     tuh Sep 07 1988
C
C     IATOM >  0 : SO -> AO transformation, modified for atom IATOM
C     IATOM =  0 : SO -> AO transformation
C     IATOM = -1 : AO -> SO transformation
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0)
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
      DIMENSION ORBTRA(IDIM,IDIM)

C
      IF (IATOM .GT. 0) THEN
         FACTOR = SQRT(FMULT(ISTBNU(IATOM)))
      ELSE IF (IATOM .EQ. 0 .OR. IATOM .EQ. -1) THEN
         FACTOR = 1.0D0
      ELSE
         CALL QENTER('TRAORB')
         CALL QUIT('IATOM .lt. -1')
      END IF
      CALL DZERO(ORBTRA,IDIM*IDIM)
      DO 100 I = 1, IDIM
         DO 200 J = 1, MAXREP + 1
            IF (IATOM .LE. 0) THEN
               K = ITRAN(I,J)
            ELSE
               K = IAOAO(ITRAN(I,J))
            END IF
            IF (K .GT. 0) THEN
               ORBTRA(I,K) = CTRAN(I,J)
            END IF
  200    CONTINUE
  100 CONTINUE
      IF (IATOM .GE. 0) THEN
         DO 300 I = 1, IDIM
            SUM = D0
            DO 400 J = 1, IDIM
               SUM = SUM + ABS(ORBTRA(I,J))
 400        CONTINUE
            SUM = 1.0D0/(FACTOR*SUM)
            DO 500 J = 1, IDIM
               ORBTRA(I,J) = ORBTRA(I,J)*SUM
 500        CONTINUE
 300     CONTINUE
      END IF
C
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('SO to AO transformation matrix',-1)
         WRITE(LUPRI,*) '   for IATOM = ',IATOM
         CALL OUTPUT(ORBTRA,1,IDIM,1,IDIM,IDIM,IDIM,-1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck testfd */
      SUBROUTINE TESTFD(IATOM,ICOOR,IPRINT,NOH1,NOH2,NODC,NODV,NOPV,
     *                  FTD,FCD,FVD,QD,DV)
C
C     This subroutine calculates the total AO differentiated Fock
C     matrix in MO basis and tests by taking the trace of the Fock
C     matrix and comparing it with the appropriate integral
C     contributions to the molecular gradient.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.D00, TWO = 2.D00)
      INTEGER P, Q, T, U, TU, QU
      LOGICAL NOH1, NOH2, NODC, NODV, NOPV
      DIMENSION FTD(NOCCT,NORBT), FCD(NNORBT), FVD(NNORBT),
     *          QD(NASHDI,NORBT), DV(*)
C
C Used from common blocks:
C  INFDIM: NASHDI
C
#include "inforb.h"
#include "infdim.h"
#include "energy.h"
      ITRI(I) = I*(I - 1)/2
C
      CALL DZERO(FTD(1,1),NOCCT*NORBT)
C
C     Inactive block
C
      IF (.NOT. NODC) THEN
         DO 200 I = 1, NISHT
            DO 210 Q = 1, NORBT
               IF (I .GE. Q) THEN
                  IQ = ITRI(I) + Q
               ELSE
                  IQ = ITRI(Q) + I
               END IF
               FTD(I,Q) = TWO*(FCD(IQ) + FVD(IQ))
  210       CONTINUE
  200    CONTINUE
      END IF
C
C     Active block
C
      DO 300 T = 1, NASHT
         DO 310 Q = 1, NORBT
            FTDTQ = D0
            IF (.NOT. NODV) THEN
               DO 400 U = 1, NASHT
                  IU = NISHT + U
                  IF (Q .GE. IU) THEN
                     QU = ITRI(Q) + IU
                  ELSE
                     QU = ITRI(IU) + Q
                  END IF
                  IF (T .GE. U) THEN
                     TU = ITRI(T) + U
                  ELSE
                     TU = ITRI(U) + T
                  END IF
                  FTDTQ = FTDTQ + DV(TU)*FCD(QU)
  400          CONTINUE
            END IF
            IF (.NOT. (NOH2 .OR. NOPV)) THEN
               FTDTQ = FTDTQ + QD(T,Q)
            END IF
            FTD(NISHT + T,Q) = FTDTQ
  310    CONTINUE
  300 CONTINUE
      IF (IPRINT .GE. 10) THEN
         CALL HEADER('AO derivative of total Fock matrix',-1)
         CALL OUTPUT(FTD(1,1),1,NOCCT,1,NORBT,NOCCT,NORBT,1,LUPRI)
      END IF
      TRACE = D0
      DO 500 P = 1, NOCCT
         TRACE = TRACE + FTD(P,P)
  500 CONTINUE
      IAC = 3*(IATOM - 1) + ICOOR
      GRDINT = D0
      IF (.NOT. NOH1) THEN
         GRDINT = GRDINT + GRADKE(IAC) + GRADNA(IAC)
      END IF
      IF (.NOT. NOH2) THEN
         GRDINT = GRDINT + TWO*GRADEE(IAC)
      END IF
      CALL HEADER('FOCK MATRIX TEST:',0)
      WRITE (LUPRI,'(8X,A,/)')
     *      ' GRD(KE+NA+2*EE) TRACE(FTD)      DIFFERENCE'
      WRITE (LUPRI,'(5X,3F16.10)') GRDINT, TRACE, GRDINT - TRACE
      RETURN
      END
C  /* Deck testtr */
      SUBROUTINE TESTTR(IATOM,ICOOR,IOPSYM,PV,H2D,QD,IPRINT)
#include "implicit.h"
#include "priunit.h"
C
C     *****************************************************
C     ***  Take trace of Q matrices and PV.H2D product  ***
C     ***  as a check of transformation.  This version  ***
C     ***  handles symmetry.                            ***
C     *****************************************************
C
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "iratdef.h"
      PARAMETER (D0 = 0.D00, D1 = 1.D00, TWO = 2.D00)
      PARAMETER (HALF = 0.5D00)
      INTEGER T, U, V, X, TU, VX
      DIMENSION PV(NNASHX,NNASHX), H2D(NNASHX,NNASHX)
      DIMENSION QD(NORBT,NASHT)
#include "symmet.h"
#include "inforb.h"
#include "energy.h"
      ITRI(I) = I*(I - 1)/2
C
      WRITE (LUPRI,'(/,A,/)')
     *       ' ---------- Output from TESTTR ---------- '
      IAC = IPTCNT(3*(IATOM-1)+ICOOR,IOPSYM-1,1)
      GRDINT = GRADEE(IAC)
C
C     ***** H2D elements *****
C
      IF (IPRINT .GT. 60) THEN
         CALL HEADER('Two-electron density elements',0)
         CALL OUTPUT(PV,1,NNASHX,1,NNASHX,NNASHX,NNASHX,1,LUPRI)
         CALL HEADER('Derivative MO two-electron integrals',0)
         CALL OUTPUT(H2D(1,1),1,NNASHX,1,NNASHX,
     *               NNASHX,NNASHX,1,LUPRI)
      END IF
      TRACE = D0
C
C     Extra factor of half for new definition of 2-matrix
C
      DO 100 T = 1, NASHT
         DO 110 U = 1, T
            TU = ITRI(T) + U
            FACTU = HALF
            IF (T .NE. U) FACTU = TWO*FACTU
            DO 120 V = 1, NASHT
               DO 130 X = 1, V
                  VX = ITRI(V) + X
                  FACVX = D1
                  IF (V .NE. X) FACVX = TWO
                  TRACE = TRACE + FACTU*FACVX*PV(TU,VX)*H2D(TU,VX)
130            CONTINUE
120         CONTINUE
110      CONTINUE
100   CONTINUE
      CALL HEADER('Test on H2D matrix:',0)
      WRITE (LUPRI,'(9X,A,/)')
     *       ' GRADEE(PV)     PV*H2D         Difference '
      WRITE (LUPRI,'(5X,3F15.10)') GRDINT, TRACE, GRDINT - TRACE
      WRITE (LUPRI,'(/,A,/)')
     *       ' Note: Only PV must be included in TWOEXP!'
C
C     ***** QD elements *****
C
      IF (IPRINT .GT. 60) THEN
         CALL HEADER('QD matrix',-1)
         CALL OUTPUT(QD(1,1),1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
      END IF
C
C     Loop over symmetries here and accumulate only active-active terms
C
      TRACE = D0
      IOFF = 0
      JOFF = 0
      DO 300 ISYM = 1,NSYM
         NASHI = NASH(ISYM)
         NISHI = NISH(ISYM)
         NORBI = NORB(ISYM)
         IF (NASHI .NE. 0) THEN
            DO 310 I = 1, NASHI
               TRACE = TRACE + QD(IOFF+NISHI+I,JOFF+I)
310         CONTINUE
         ENDIF
         IOFF = IOFF + NORBI
         JOFF = JOFF + NASHI
300   CONTINUE
      GRDINT = TWO*GRDINT
      CALL HEADER('Test on QD matrix:',0)
      WRITE (LUPRI,'(9X,A,/)')
     *       ' 2*GRADEE(PV)   Trace(QD)      Difference '
      WRITE (LUPRI,'(5X,3F15.10)') GRDINT, TRACE, GRDINT - TRACE
      WRITE (LUPRI,'(/,A,/)')
     *       ' Note: Only PV must be included in TWOEXP!'
      RETURN
      END
C  /* Deck testsd */
      SUBROUTINE TESTSD(TD,FT,IATOM,DERIV,IPRINT)
C
C     This routine performs a test on differentiated overlap matrices by
C     calculating the gradient reorthonormalization in MO basis and
C     comparing it with the results obtained in AO basis.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "iratdef.h"
      PARAMETER (D0 = 0.D00, TWO = 2.D00)
      LOGICAL DERIV(3)
      DIMENSION TD(NNORBX,3), FT(*), TRACE(3)
#include "symmet.h"
#include "inforb.h"
#include "energy.h"
C
      CALL HEADER('Reorthonormalization of gradient from'//
     *            ' overlap matrices in MO basis',0)
      IF (IPRINT .GT. 5) THEN
         DO 50 ISYM = 1, NSYM
            WRITE (LUPRI,'(A,I5)') ' Irrep', ISYM
            NORBI = NORB(ISYM)
            ISTRI = I2ORB(ISYM) + 1
            CALL HEADER('Total Fock matrix',-1)
            CALL OUTPUT(FT(ISTRI),1,NORBI,1,NORBI,NORBI,NORBI,1,LUPRI)
   50    CONTINUE
      END IF
      DO 100 ICOOR = 1, 3
         IF (DERIV(ICOOR)) THEN
            IF (IPRINT .GT. 5) THEN
               CALL HEADER('TD matrix',-1)
               CALL OUTPAK(TD(1,ICOOR),NORBT,1,LUPRI)
            END IF
            TRAC = D0
            DO 150 ISYM = 1, NSYM
               ISTR = IORB(ISYM) + 1
               IEND = IORB(ISYM) + NORB(ISYM)
               IJ   = I2ORB(ISYM)
               DO 200 I = ISTR, IEND
                  DO 300 J = ISTR, IEND
                     IJ = IJ + 1
                     IJTD = (MAX(I,J)*(MAX(I,J) - 3))/2 + I + J
                     TRAC = TRAC + FT(IJ)*TD(IJTD,ICOOR)
  300             CONTINUE
  200          CONTINUE
  150       CONTINUE
            TRACE(ICOOR) = TWO*TRAC
         END IF
  100 CONTINUE
      IOFF = 3*(IATOM - 1)
      CALL HEADER
     *   ('Coordinate     AO basis       MO basis       Difference',10)
      DO 500 ICOOR = 1, 3
         IF (DERIV(ICOOR)) THEN
            JCOOR = IPTCNT(IOFF + ICOOR,0,1)
            IF (JCOOR .GT. 0) THEN
               FS   = GRADFS(JCOOR)
               TRAC = TRACE(ICOOR)
               WRITE (LUPRI,'(5X,I10,5X,3F15.10)')
     *               JCOOR, FS, TRAC, FS - TRAC
            END IF
         END IF
  500 CONTINUE
      RETURN
      END
C  /* Deck onetra */
      SUBROUTINE ONETRA(CMO,SQMO,SQSO,WORK,LWORK,IATOM,ICOOR,IPRINT)
C
C     This routine transforms derivative SO matrices to MO basis.
C
C     tuh 010988
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "symmet.h"
#include "inforb.h"
#include "infdim.h"
      DIMENSION CMO(NCMOT), WORK(LWORK),
     *          SQMO(NORBT,NORBT), SQSO(N2BASX)
      LOGICAL   NONSYM
      PARAMETER (D1 =1.0D0, D2 = 2.0D0, DP5 = 0.5D0)
C
C     ***** Print Section *****
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from ONETRA',-1)
         WRITE (LUPRI,'(A,I5)') ' IATOM  ', IATOM
         WRITE (LUPRI,'(A,I5)') ' ICOOR  ', ICOOR
         WRITE (LUPRI,'(A,I5)') ' NORBT  ', NORBT
         WRITE (LUPRI,'(A,I5)') ' NBAST  ', NBAST
         WRITE (LUPRI,'(A,2I5)')' NNBASX, N2BASX ', NNBASX, N2BASX
         WRITE (LUPRI,'(A,2I5)')' NNORBX, N2ORBX ', NNORBX, N2ORBX
         WRITE (LUPRI,'(A,I5)') ' NCMOT  ', NCMOT
         WRITE (LUPRI,'(A,I10)')' LWORK  ', LWORK
      END IF
C
C     ***** Make square matrix *****
C
      IF (IPRINT .GT. 15) THEN
         CALL HEADER('Square SO matrix',-1)
         CALL OUTPUT(SQSO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     ***** Transform matrix from SO to MO basis *****
C
      CALL DZERO(SQMO,N2ORBX)
C
C     Loop over irreps for orbitals
C
      NONSYM = .FALSE.
      JCOOR = MAX(1,3*(IATOM - 1) + ICOOR)
C     ... hjaaj: to avoid referring to IPTCNT(-2,...) below when IATOM.eq.0
      DO 100 IREPA = 1,NSYM
      IF (NORB(IREPA) .GT. 0) THEN
         DO 200 IREPB = 1, IREPA
         IF (NORB(IREPB) .GT. 0) THEN
C
C           Check if operator belongs to this symmetry
C
            IREPO = MULD2H(IREPA,IREPB)
            IF ((IATOM.EQ.0 .AND. (IREPO-1).EQ.ISYMAX(ICOOR,2)) .OR.
     &          (IATOM.GT.0 .AND. IPTCNT(JCOOR,IREPO-1,1).GT.0)) THEN
               NONSYM = NONSYM .OR. IREPO .GT. 1
C
C              Print symmetries and orbitals
C
               IF (IPRINT.GT.15) THEN
                  WRITE (LUPRI,'(A,3I3)') ' IREPA/B/O',IREPA,IREPB,IREPO
                  WRITE (LUPRI,'(/A,I5,A)')
     *               ' MO coefficients for symmetry', IREPA
                  CALL OUTPUT(CMO(ICMO(IREPA)+1),1,NBAS(IREPA),1,
     *               NORB(IREPA),NBAS(IREPA),NORB(IREPA),1,LUPRI)
                  WRITE (LUPRI,'(/A,I5,A)')
     *               ' MO coefficients for symmetry', IREPB
                  CALL OUTPUT(CMO(ICMO(IREPB)+1),1,NBAS(IREPB),1,
     *               NORB(IREPB),NBAS(IREPB),NORB(IREPB),1,LUPRI)
               END IF
C
C              Transform matrix block(s)
C
               CALL UTHV(CMO(ICMO(IREPA)+1),SQSO,CMO(ICMO(IREPB)+1),
     *                   IREPA,IREPB,NBAS(IREPA),NBAS(IREPB),SQMO,WORK)
               IF (NONSYM) THEN
                  CALL UTHV(CMO(ICMO(IREPB)+1),SQSO,CMO(ICMO(IREPA)+1),
     &                 IREPB,IREPA,NBAS(IREPB),NBAS(IREPA),SQMO,WORK)
               END IF
C
C              Print transformed matrix thus far
C
               IF (IPRINT.GT.25) THEN
                  WRITE(LUPRI,'(/4A)')' Unfinished matrix in MO basis'
                  CALL OUTPUT(SQMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
            END IF
         END IF
 200     CONTINUE
      END IF
 100  CONTINUE
      IF (IPRINT.GT.15) THEN
         CALL HEADER('Matrix in MO basis',-1)
         CALL OUTPUT(SQMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck transx */
       SUBROUTINE TRANSX(A,B,NA,NB,ANTSYM,IPRINT)
C
C Modified from TRANSA, tuh 2-Sep-1988
C TRANSA written 1-Feb-1984 PJ
C
C Put the lower triangle of A into the upper triangle of B
C (if A = B, the matrix will become (anti)symmetric
C  based on the lower triangle)
C
#include "implicit.h"
#include "priunit.h"
C
C
       DIMENSION A(NA,*), B(NB,*)
       IF (IPRINT.GT.15) THEN
          WRITE(LUPRI,150) NA,NB,ANTSYM
 150      FORMAT(/' Dimension of matrix A is ',I5,I5,
     &           /' ANTSYM =',F10.2)
          WRITE(LUPRI,'(/A)')' Matrix to be (anti)symmetrized '
          CALL OUTPUT(A(1,1),1,NA,1,NB,NA,NB,1,LUPRI)
       END IF
       DO 200 I=1,NB
         DO 100 J=1,I
           B(J,I)=A(I,J)*ANTSYM
 100     CONTINUE
 200   CONTINUE
       IF (IPRINT.GT.15) THEN
          WRITE(LUPRI,'(/A)')' Matrix after (anti)symmetrization '
          CALL OUTPUT(B(1,1),1,NB,1,NA,NB,NA,1,LUPRI)
       END IF
       RETURN
      END
C  /* Deck getsd */
      SUBROUTINE GETSD(DERIV,CMO,TD,WORK,LWORK,IATOM,DIAGTD,
     &                 KEY,DOTD,IPRINT,INTPRI)
C
C     This subroutine calculates the AO-differentiated overlap matrices
C     in MO basis.
C
C
C     Due to the need of square matrices in constructing right-hand
C     sides when using differentiated creation operators, all matrices
C     in this routine is square, all packing is to be done outside this
C     routine, K.Ruud, Nov.-93
C     
C     tuh 051088
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      LOGICAL DERIV(3), DIAGTD, DIVIDE, DOTD
      PARAMETER (D1 = 1.0D0)
      CHARACTER*4 KEY
#include "abainf.h"
#include "inforb.h"
      DIMENSION CMO(*), TD(N2BASX,3), WORK(LWORK)
      CALL QENTER('GETSD')
      IF (IPRINT .GT. 05) THEN
         CALL TITLER('Output from GETSD','*',103)
      END IF
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, '(A,3L5)') ' DERIV ', (DERIV(I),I=1,3)
         WRITE (LUPRI, '(A,I10)') ' LWORK ', LWORK
         WRITE (LUPRI, '(A,I5)')  ' IATOM ', IATOM
         WRITE (LUPRI, '(A,L5)')  ' DIAGTD', DIAGTD
      END IF
      CALL DZERO(TD,3*N2BASX)
      CALL ONEDRL(KEY,TD(1,1),TD(1,2),TD(1,3),IATOM,
     &            DERIV(1),DERIV(2),DERIV(3),WORK,LWORK,IPRINT,INTPRI)
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('SO overlap matrix - x direction',-1)
         CALL OUTPUT(TD(1,1),1,NBAST,1,NBAST,NBAST,NBAST,
     &               1,LUPRI)
         CALL HEADER('SO overlap matrix - y direction',-1)
         CALL OUTPUT(TD(1,2),1,NBAST,1,NBAST,NBAST,NBAST,
     &               1,LUPRI)
         CALL HEADER('SO overlap matrix - z direction',-1)
         CALL OUTPUT(TD(1,3),1,NBAST,1,NBAST,NBAST,NBAST,
     &               1,LUPRI)
      END IF 
      KSMOSQ = 1
      KGETSD = KSMOSQ +   N2BASX
      LGETSD = LWORK  -   KGETSD + 1
      IF (KGETSD.GT.LWORK) CALL STOPIT('GETSD',' ',KGETSD,LWORK)
      DO 100 ICOOR = 1, 3
         IF (DERIV(ICOOR)) THEN
            IF (IPRINT .GE. 10) THEN
               IF (ICOOR. EQ. 1) THEN
                  CALL TITLER('Differentiation in x direction',' ',200)
               ELSE IF (ICOOR. EQ. 2) THEN
                  CALL TITLER('Differentiation in y direction',' ',200)
               ELSE
                  CALL TITLER('Differentiation in z direction',' ',200)
               END IF
            END IF
            CALL ONETRA(CMO,WORK(KSMOSQ),TD(1,ICOOR),WORK(KGETSD),
     &                  LGETSD,IATOM,ICOOR,IPRINT)
            CALL DCOPY(N2ORBX,WORK(KSMOSQ),1,TD(1,ICOOR),1)
            IF (DOTD) THEN
               IF (KEY .EQ. 'SMAT' .OR. KEY .EQ. 'OMAT') THEN
                  DIVIDE = .TRUE.
               ELSE
                  DIVIDE = .FALSE.
               END IF
               CALL GETTD(TD(1,ICOOR),NORBT,DIAGTD,DIVIDE,IPRINT)
               IF (IPRINT .GE. 10) THEN
                  CALL HEADER('TD matrix',-1)
                  CALL OUTPUT(TD(1,ICOOR),1,NORBT,1,NORBT,NORBT,NORBT,
     &                        1,LUPRI)
               END IF
            ELSE
               IF (IPRINT .GE. 10) THEN
                  CALL HEADER('SO overlap matrix',-1)
                  CALL OUTPUT(TD(1,ICOOR),1,NORBT,1,NORBT,NORBT,NORBT,
     &                        1,LUPRI)
               END IF
            END IF
         END IF
  100 CONTINUE
      CALL QEXIT('GETSD')
      RETURN
      END
C  /* Deck onedrl */
      SUBROUTINE ONEDRL(KEY,OMATX,OMATY,OMATZ,IATOM,DERX,DERY,DERZ,WORK,
     &                  LWORK,IPRINT,INTPRI)
C
         use pelib_interface, only: use_pelib, pelib_ifc_london
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0)
      INTEGER TMAT
      LOGICAL DERX, DERY, DERZ
      CHARACTER*4 KEY
      CHARACTER*7 KEY2
#include "inforb.h"
#include "gnrinf.h"
#include "abainf.h"
#include "cbisol.h"
#include "infinp.h"
#include "pcmlog.h"
#include "qmmm.h"
#include "efield.h"
#include "infpar.h"
      DIMENSION OMATX(N2BASX), OMATY(N2BASX), OMATZ(N2BASX),
     &          WORK(LWORK)
      IF (IPRINT .GT. 5) CALL TITLER('Output from ONEDRL','*',103)
C
C     Field derivatives
C     =================
C
      IF (IATOM .EQ. 0) THEN
         KMATS = 1
         KWRK  = KMATS + 3*N2BASX
         LWRK  = LWORK - KWRK + 1
         IF (KWRK .GT. LWORK) CALL STOPIT('ONEDRL','GET1PR',KWRK,LWORK)
         KEY2 = 'S1MAG  '
         IF (KEY .EQ. 'BMAT') KEY2 = 'S1MAGR '
         IF ((KEY .EQ. 'HMAT') .AND. (.NOT. NOLOND)) KEY2 = 'MAGMOM '
         IF ((KEY .EQ. 'HMAT') .AND. NOLOND)         KEY2 = 'ANGMOM '
         NCOMP = 3
         IF (KEY2 .EQ. 'S1MAGR ') THEN
            CALL GET1PR(WORK(KMATS),KEY2,NCOMP,'SQUARE',.FALSE.,
     &                  WORK(KWRK),LWRK,IDUMMY,INTPRI)
         ELSE
            CALL GET1PR(WORK(KMATS),KEY2,NCOMP,'SQUARE',.TRUE.,
     &                  WORK(KWRK),LWRK,IDUMMY,INTPRI)
         END IF
         IF (KEY2.EQ. 'ANGMOM') CALL DSCAL(3*N2BASX,DP5,WORK(KMATS),1)
C
C     Add contribution from electric field
C
         IF (NFIELD .GT. 0 .AND. (KEY .EQ. 'HMAT')
     &                     .AND. (.NOT. NOLOND)) THEN
            KTMAT = KWRK
            MWRK  = KTMAT + 3*N2BASX
            IF (MWRK .GT. LWRK) CALL STOPIT('ONEDRL','CM1   ',MWRK,
     &                                       LWRK)
            LLEFT = LWRK - MWRK + 1
            DO 556 IFIELD = 1, NFIELD
               IF (EFIELD(IFIELD) .NE. D0) THEN
                  IF (LFIELD(IFIELD)(2:7) .EQ. 'DIPLEN') THEN
                     IF (LFIELD(IFIELD) .EQ. 'XDIPLEN ') THEN
                        FIELD1 = 'X-FIELD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'YDIPLEN ') THEN
                        FIELD1 = 'Y-FIELD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'ZDIPLEN ') THEN
                        FIELD1 = 'Z-FIELD'
                     ELSE
                        WRITE (LUPRI,'(/,3A,/)') 'Field type ',
     &                       LFIELD(IFIELD),
     &                       ' not implemented for magnetic properties'
                        CALL QUIT('Illegal field type for '//
     &                            'magnetic properties')
                     END IF
                     NCOMP = 3
                     CALL GET1PR(WORK(KTMAT),'CM1    ',NCOMP,'SQUARE',
     &                           .TRUE.,WORK(MWRK),LLEFT,IDUMMY,INTPRI)
                  ELSE IF (LFIELD(IFIELD)(3:7) .EQ. 'THETA' .OR.
     &                     LFIELD(IFIELD)(3:8) .EQ. 'SECMOM') THEN
                     IF (LFIELD(IFIELD) .EQ. 'XXTHETA ') THEN
                        FIELD3 = 'XX-FGRD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'YYTHETA ') THEN
                        FIELD3 = 'YY-FGRD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'ZZTHETA ') THEN
                        FIELD3 = 'ZZ-FGRD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'XYTHETA ') THEN
                        FIELD3 = 'XY-FGRD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'XZTHETA ') THEN
                        FIELD3 = 'XZ-FGRD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'YZTHETA ') THEN
                        FIELD3 = 'YZ-FGRD'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'XXSECMOM') THEN
                        FIELD3 = 'XX-2MOM'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'YYSECMOM') THEN
                        FIELD3 = 'YY-2MOM'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'ZZSECMOM') THEN
                        FIELD3 = 'ZZ-2MOM'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'XYSECMOM') THEN
                        FIELD3 = 'XY-2MOM'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'XZSECMOM') THEN
                        FIELD3 = 'XZ-2MOM'
                     ELSE IF (LFIELD(IFIELD) .EQ. 'YZSECMOM') THEN
                        FIELD3 = 'YZ-2MOM'
                     ELSE
                        WRITE (LUPRI,'(/,3A,/)') 'Field type ',
     &                       LFIELD(IFIELD),
     &                       ' not implemented for magnetic properties'
                        CALL QUIT('Illegal field type for '//
     &                            'magnetic properties')
                     END IF
                     NCOMP = 3
                     CALL GET1PR(WORK(KTMAT),'QDBINT ',NCOMP,'SQUARE',
     &                           .TRUE.,WORK(MWRK),LLEFT,IDUMMY,INTPRI)
                  END IF
                  CALL DAXPY(N2BASX*3,EFIELD(IFIELD),WORK(KTMAT),1,
     &                    WORK(KMATS),1)
               END IF
 556        CONTINUE
         END IF
C
         IF ((SOLVNT .OR. PCM) .AND. (KEY .EQ. 'HMAT')
     &                         .AND. (.NOT. NOLOND)) THEN
            KTMAT = KWRK
            MWRK  = KTMAT + 3*N2BASX
            IF (MWRK .GT. LWRK) CALL STOPIT('ONEDRL','MAGSOL',MWRK,
     &                                      LWRK)
            LLEFT = LWRK - MWRK + 1
            IF (SOLVNT) THEN
               CALL MAGSOL(WORK(KTMAT),3,IPRINT,INTPRI,
     &                     WORK(MWRK),LLEFT)
            ELSE
               CALL MAGPCMFIRST(WORK(KTMAT),3,IPRINT,INTPRI,WORK(MWRK),
     &                          LLEFT)
            END IF
            CALL DAXPY(N2BASX*3,D1,WORK(KTMAT),1,WORK(KMATS),1)
         END IF
C
         IF (QM3 .AND. (KEY .EQ. 'HMAT') .AND. .NOT. NOLOND) THEN
            KTMAT = KWRK
            MWRK  = KTMAT + 3*N2BASX
            IF (MWRK .GT. LWRK) CALL STOPIT('ONEDRL','QM3-1',MWRK,LWRK)
            LLEFT = LWRK - MWRK + 1
C     Calls QM3FIRST_P if parallel and with #nodes>1, Arnfinn nov. -08
            IF (NODTOT.GE.1) THEN
               CALL QM3FIRST_P(WORK(KTMAT),3,IPRINT,INTPRI,WORK(MWRK),
     &                       LLEFT)
            ELSE
               CALL QM3FIRST(WORK(KTMAT),3,IPRINT,INTPRI,WORK(MWRK),
     &                       LLEFT)
            ENDIF
            CALL DAXPY(N2BASX*3,D1,WORK(KTMAT),1,WORK(KMATS),1)
         END IF

         IF (QMMM .AND. (KEY .EQ. 'HMAT') .AND. .NOT. NOLOND) THEN
            IF (IPQMMM .GT. 1) CALL TIMER('START ',TIMSTR,TIMEND)
            KTMAT = KWRK
            MWRK  = KTMAT + 3*N2BASX
            IF (MWRK .GT. LWRK) CALL STOPIT('ONEDRL','QMMM-1',MWRK,LWRK)
            LLEFT = LWRK - MWRK + 1
            CALL QMMMFIRST(WORK(KTMAT),IPRINT,INTPRI,WORK(MWRK),LLEFT)
            CALL DAXPY(N2BASX*3,D1,WORK(KTMAT),1,WORK(KMATS),1)
            IF (IPQMMM .GT. 1) CALL TIMER('QMMMFIRST',TIMSTR,TIMEND)
         END IF

         IF (USE_PELIB() .AND. (KEY.EQ.'HMAT') .AND. .NOT. NOLOND) THEN
            KTMAT = KWRK
            MWRK = KTMAT + 3 * N2BASX
            IF (MWRK > LWRK) CALL STOPIT('ONEDRL', 'PEQM', MWRK, LWRK)
            CALL PELIB_IFC_LONDON(WORK(KTMAT))
            CALL DAXPY(3*N2BASX, D1, WORK(KTMAT), 1, WORK(KMATS), 1)
         END IF

         IF (DERX) CALL DCOPY(N2BASX,WORK(KMATS         ),1,OMATX,1)
         IF (DERY) CALL DCOPY(N2BASX,WORK(KMATS+  N2BASX),1,OMATY,1)
         IF (DERZ) CALL DCOPY(N2BASX,WORK(KMATS+2*N2BASX),1,OMATZ,1)
C
C     Geometrical derivatives
C     =======================
C
      ELSE
         IF (LWORK.LT.(NNBASX + N2BASX))
     &      CALL STOPIT('ONEDRL','ONEDER',LWORK,N2BASX)
         IF (DERX) THEN
            CALL ONEDER(KEY,OMATX,WORK(1),WORK(1+NNBASX),
     &                  IATOM,1,IPRINT)
         END IF
         IF (DERY) THEN
            CALL ONEDER(KEY,OMATY,WORK(1),WORK(1+NNBASX),
     &                  IATOM,2,IPRINT)
         END IF
         IF (DERZ) THEN
            CALL ONEDER(KEY,OMATZ,WORK(1),WORK(1+NNBASX),
     &                  IATOM,3,IPRINT)
         END IF
      END IF
C
C     Print
C     =====
C
      IF (IPRINT .GE. 15) THEN
         WRITE (LUPRI,'(2X,2A,I3)') 'Integral type and atom:',KEY,IATOM
         IF (DERX) THEN
            CALL HEADER('X component in ONEDRL',-1)
            CALL OUTPUT(OMATX,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DERY) THEN
            CALL HEADER('Y component in ONEDRL',-1)
            CALL OUTPUT(OMATY,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
         IF (DERZ) THEN
            CALL HEADER('Z component in ONEDRL',-1)
            CALL OUTPUT(OMATZ,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck oneder */
      SUBROUTINE ONEDER(KEY,ONEMAT,TRONEM,WORK,ICENT,ICOOR,IPRINT)
C
C     130988 tuh
C
C     This subroutine retrieves one-electron differentiated matrices
C     (one-electron Hamiltonian or overlap)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
      LOGICAL DNULL, TOTSYM, FNDLAB
      CHARACTER*4 KEY
      CHARACTER*8 LABEL
#include "nuclei.h"
#include "symmet.h"
#include "inforb.h"
#include "inftap.h"
#include "infinp.h"
#include "chrnos.h"
      DIMENSION WORK(N2BASX), ONEMAT(N2BASX), TRONEM(NNBASX)

C
      IF (IPRINT .GE. 10) CALL TITLER('Output from ONEDER','*',103)
      CALL DZERO(TRONEM,NNBASX)
      CALL DZERO(WORK,N2BASX)
      DO 100 IREPO = 0, MAXREP
         TOTSYM = IREPO .EQ. 0
         ISCOOR = IPTCNT(3*(ICENT-1)+ICOOR,IREPO,1)
         IF (ISCOOR .GT. 0) THEN
            IF (KEY .EQ. 'OMAT' .OR. KEY .EQ. 'HMAT') THEN
               DO 200 IREPA = 0, MAXREP
                  IREPB = IEOR(IREPO,IREPA)
                  IF (IREPA .GE. IREPB) THEN
                     NDIMA = NBAS(IREPA + 1)
                     NDIMB = NBAS(IREPB + 1)
                     NDIM  = NDIMA*NDIMB
                     IF (NDIM .GT. 0) THEN
                        IF (TOTSYM) NDIM = NDIMA*(NDIMA + 1)/2
C
C                       ***** Get elements *****
C
                        CALL GETONE(KEY,ISCOOR,IREPA,NDIM,WORK,DNULL)
C
C                       ***** Print section *****
C
                        IF (IPRINT .GE. 15) THEN
                           WRITE (LUPRI,'(//A,2X,A,I5)')
     *                          ' Coordinate and symmetry: ',
     *                          NAMEX(3*ICENT - 3 + ICOOR), IREPO
                           WRITE (LUPRI,'(/A,2I5)')
     *                          ' Symmetry of orbitals:', IREPA, IREPB
                           WRITE (LUPRI,'(A,2I5)')
     *                          ' Number of orbitals:  ', NDIMA, NDIMB
                           IF (.NOT.DNULL) THEN
                              IF (TOTSYM) THEN
                                 CALL OUTPAK(WORK,NDIMA,1,LUPRI)
                              ELSE
                                 CALL OUTPUT(WORK,1,NDIMB,1,NDIMA,
     *                                NDIMB,NDIMA,1,LUPRI)
                              END IF
                           END IF
                        END IF
C
C                       ***** Add elements to matrix *****
C
                        IF (.NOT.DNULL) THEN
                           IJ = 0
                           IOFFA = IBAS(IREPA + 1)
                           IOFFB = IBAS(IREPB + 1)
                           DO 300 I = IOFFA + 1, IOFFA + NDIMA
                              DO 400 J = IOFFB + 1, MIN(IOFFB + NDIMB,I)
                                 IJ = IJ + 1
                                 TRONEM((I*(I - 1))/2 + J) = WORK(IJ)
 400                          CONTINUE
 300                       CONTINUE
                        END IF
                     END IF
                  END IF
 200           CONTINUE
               IF (NFIELD .GT. 0 .AND. KEY .EQ.'HMAT') THEN
                  CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
                  DO 556 IFIELD = 1, NFIELD
                     IF (EFIELD(IFIELD) .NE. 0.D0) THEN
                        IF (LFIELD(IFIELD)(2:8) .NE. 'DIPLEN ') THEN
                           WRITE (LUPRI,'(/,3A,/)') 'Field type ',
     &                          LFIELD(IFIELD),
     &                  ' not implemented for geometric derivatives'
                           CALL QUIT('Illegal field type for '//
     &                          'geometric derivative')
                        END IF
C
C     Determine label to search for
C
                        LABEL = CHRNOS(ISCOOR/100)
     &                       //CHRNOS(MOD(ISCOOR,100)/10)
     &                       //CHRNOS(MOD((MOD(ISCOOR,100)),10))
     &                       //'DPG '// LFIELD(IFIELD)(1:1)
                        REWIND LUPROP
                        IF (.NOT. FNDLAB(LABEL,LUPROP)) THEN
                           WRITE(LUPRI,'(1X,3A)')'Label ',
     &                          LABEL, ' not found on file AOPROPER'
                           CALL QUIT('Error in ONEDER')
                        END IF
                        CALL READT(LUPROP,NNBASX,WORK)
                        CALL DAXPY(NNBASX,-EFIELD(IFIELD),
     &                       WORK,1,TRONEM,1)
                     END IF
 556              CONTINUE
                  CALL GPCLOSE(LUPROP,'KEEP')
               END IF
            ELSE
               CALL GETDMA(ISCOOR,ONEMAT,WORK,DNULL)
            END IF
         END IF
 100  CONTINUE
      IF (KEY .EQ. 'OMAT' .OR. KEY .EQ. 'HMAT') THEN
         CALL DSPTSI(NBAST,TRONEM,ONEMAT)
      END IF
      IF (IPRINT .GE. 10) THEN
         IF (KEY .EQ. 'HMAT') THEN
            CALL HEADER('One-electron differentiated Hamiltonian',-1)
         ELSE
            CALL HEADER('Differentiated overlap matrix',-1)
         END IF
         CALL OUTPUT(ONEMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck gettd */
      SUBROUTINE GETTD(TD,NORBT,DIAGTD,DIVIDE,IPRINT)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D1 = 1.0D0, DP5 = 0.5D00, DP25 = 0.25D00)
      LOGICAL DIAGTD, DIVIDE
      DIMENSION TD(NORBT,NORBT)
C
      IF (DIAGTD) THEN
         CALL DZERO(TD,NORBT*NORBT)
         DO 100 I = 1, NORBT
            TD(I,I) = DP25
  100    CONTINUE
      ELSE
         FACTOR = - D1
         IF (DIVIDE) FACTOR = - DP5
         CALL DSCAL(NORBT*NORBT,FACTOR,TD,1)
      END IF
      RETURN
      END
C  /* Deck getdma */
      SUBROUTINE GETDMA(ISCOOR,ONEMAT,WORK,DNULL)
C
C     Collect the half differentiated overlap matrice needed when using
C     new connections. K.Ruud June-94
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (D1 = 1.0D0)
      CHARACTER*8 LABEL
      DIMENSION WORK(N2BASX), ONEMAT(N2BASX)
      LOGICAL DNULL, FNDLAB, OPTST
#include "inforb.h"
#include "inftap.h"
#include "chrnos.h"
C
      CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
C
C     Determine label to search for
C
      LABEL = 'SQHDR'//CHRNOS(ISCOOR/100)
     &     //CHRNOS(MOD(ISCOOR,100)/10)
     &     //CHRNOS(MOD((MOD(ISCOOR,100)),10))
      REWIND LUPROP
      IF (.NOT. FNDLAB(LABEL,LUPROP)) THEN
         WRITE(LUPRI,'(1X,3A)')
     &        'Label ', LABEL, ' not found on file AOPROPER'
         CALL QUIT('Error in GETSD')
      END IF
      CALL READT(LUPROP,N2BASX,WORK)
      CALL DAXPY(N2BASX,-D1,WORK,1,ONEMAT,1)
      IF (DNRM2(N2BASX,WORK,1) .GT. 0.0D0) THEN
         DNULL = .FALSE.
      ELSE
         DNULL = .TRUE.
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
      RETURN
      END
