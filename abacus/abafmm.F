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
C  /* Deck fmminp */
      SUBROUTINE FMMINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 4)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <cbifmm.h>
#include <abainf.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.STOP  ','.LMAX  '/
C
      NEWDEF = (WORD .EQ. '*FMM   ')
C
      LMAXDF = 2
      LMAX   = LMAXDF
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FMMINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FMMINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  CUT  = .TRUE.
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) LMAX
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in FMMINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in FMMINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for FMMDRV:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' FMMDRV skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in FMMDRV:',IPRINT
            END IF
            IF (LMAX .NE. LMAXDF) THEN
               WRITE (LUPRI,'(A,I5)')
     &            ' Order of multipole expansion: ',LMAX
            END IF
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after FMMDRV.'
            END IF
         END IF
      END IF
      RETURN
      END
C
C  /* Deck fmmini */
      SUBROUTINE FMMINI
C     Initialize /CBIFMM/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbifmm.h>
C     CHARACTER OUTCEN*60                          !maw
C     PARAMETER (OUTCEN =                         !maw
C    &     "/people/mark/desktop/code/data/FMM_g-data.out")
C
C     Connecting output for Gaussian product co-ordinate data
C      to be written in MOMINT subroutine
C
C     OPEN (UNIT=18, FILE=OUTCEN)                    !maw
C
      IPRINT = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      LMAXDF = 2
      LMAX   = 2
      RETURN
      END
C
C  /* Deck fmmdrv */
      SUBROUTINE FMMDRV(WORK,LWORK)
C
C     M. Watson & T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <maxmom.h>
#include <iratdef.h>
#include <inforb.h>
#include <cbifmm.h>
      DIMENSION WORK(LWORK)
C
      IF (SKIP) RETURN
      IF (NSYM .GT. 1) THEN
         WRITE (LUPRI,'(2X,2A)')
     &      ' Symmetry not implemented in FMMDRV.',
     &      ' Program stopped.'
         CALL QUIT('Symmetry not implemented in FMMDRV.')
      END IF
      IPRINT = 0
C
      NCOMP = (LMAX + 1)**2
C
C     Allocations
C
      KINTS = 1
      KDMAT = KINTS + NCOMP*NNBASX
      KAOCR = KDMAT + NNBASX
      KAOEX = KAOCR + 3*NBAST
      KLAST = KAOEX + NBAST
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMDRV',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
      CALL FMMDR1(WORK(KINTS),WORK(KDMAT),WORK(KAOCR),WORK(KAOEX),LMAX,
     &            NCOMP,WORK(KLAST),LWRK,IPRINT)
C
      RETURN
      END
C
C  /* Deck fmmdr1 */
      SUBROUTINE FMMDR1(SPHINT,DENMAT,AOCOOR,AOEXP,LMAX,NCOMP,
     &                  WORK,LWORK,IPRINT)
C
C     M. Watson & T. Helgaker
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxmom.h>
#include <maxorb.h>
#include <aovec.h>
C
      DIMENSION SPHINT(NNBASX,NCOMP), DENMAT(NNBASX),
     &          AOCOOR(NBAST,3), AOEXP(NBAST),
     &          WORK(LWORK)
C
#include <inforb.h>
#include <inftap.h>
#include <infinp.h>
#include <nuclei.h>
#include <cbiher.h>
#include <cbieri.h>
#include <clsfmm.h>
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMDR1','*',103)
C
C     Calculate nonclassical energy
C     =============================
C
      CALL FMMNCL(WORK,LWORK,IPRINT)
C
C     AO density matrix
C     =================
C
      CALL FMMORB(AOCOOR,AOEXP,IPRINT)
C
C     AO density matrix
C     =================
C
      CALL FMMDNS(DENMAT,WORK,LWORK,IPRINT)
C
C     Multipole moment integrals
C     ==========================
C
      CALL FMMINT(SPHINT,LMAX,WORK,LWORK,IPRINT)
C
C     Expectation values
C     ==================
C
      KDISM = 1
      KLAST = KDISM + 8*NNBASX
      IF (KLAST .GT. LWORK) CALL STOPIT('SPHAVE',' ',KLAST,LWORK)
      CALL SPHAVE(SPHINT,DENMAT,LMAX,WORK(KDISM),IPRINT)
C
      RETURN
      END
C  /* Deck fmmncl */
      SUBROUTINE FMMNCL(WORK,LWORK,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxmom.h>
#include <maxorb.h>
#include <aovec.h>
C
      PARAMETER (D0 = 0.0D0, D100 = 1.0D2)
      DIMENSION WORK(LWORK)
C
#include <inforb.h>
#include <inftap.h>
#include <infinp.h>
#include <nuclei.h>
#include <cbiher.h>
#include <cbieri.h>
#include <clsfmm.h>
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMNCL','*',103)
C
      KDMAT  = 1
      KHMAT  = KDMAT  + N2BASX
      KCCFBT = KHMAT  + NNBAST
      KINDXB = KCCFBT + MXCONT*MXPRIM
      KLAST  = KINDXB + 8*MXSHEL*MXCONT
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMDR1','2EL',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
C     Density matrix
C     ==============
C
      CALL GETDMT(WORK(KDMAT),1,WORK(KLAST),LWRK,.FALSE.,.TRUE.,.TRUE.,
     &            IPRINT)
C
C     One-electron interactions
C     =========================
C
C     total
C
      LUONEL = -1
      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
      CALL READIN(WORK(KLAST),LWRK,.TRUE.)
      CALL ONEDRV(WORK(KLAST),LWRK,IPRINT,.FALSE.,0,.FALSE.,.TRUE.,
     &            .TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)
c     call onedrv(work,lwork,iprint,propty,maxdif,difint,nodc,
C    &            nodv,difdip,difqdp,hfonly,nclone,pcm)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KHMAT),NNBAST)
      CALL GET1EL(ERGONE,WORK(KDMAT),WORK(KHMAT),IPRINT)
      CALL RDONEL('KINETINT',.TRUE.,WORK(KHMAT),NNBAST)
      CALL GET1EL(ERGKIN,WORK(KDMAT),WORK(KHMAT),IPRINT)
C
C     nonclassical
C
      LUONEL = -1
      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
      CALL READIN(WORK(KLAST),LWRK,.TRUE.)
      CALL ONEDRV(WORK(KLAST),LWRK,IPRINT,.FALSE.,0,.FALSE.,.TRUE.,
     &            .TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.)
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KHMAT),NNBAST)
      CALL GET1EL(ERGNCL,WORK(KDMAT),WORK(KHMAT),IPRINT)
C
C     contributions
C
      EL1TOT = ERGONE - ERGKIN
      EL1NON = ERGNCL - ERGKIN
      EL1CLS = EL1TOT - EL1NON
C
C     Two-electron interactions
C     =========================
C
C     total
C
      LUINTA = -1
      CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL ER2INT(WORK(KCCFBT),WORK(KINDXB),WORK(KLAST),LWRK)
      CALL GETERG(ALLTOT,COUTOT,WORK(KDMAT),WORK(KLAST),LWRK)
C
C     nonclassical
C
      LUINTA = -1
      CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      NCLERI = .TRUE.
      CALL ER2INT(WORK(KCCFBT),WORK(KINDXB),WORK(KLAST),LWRK)
      NCLERI = .FALSE.
      CALL GETERG(ALLNCL,COUNCL,WORK(KDMAT),WORK(KLAST),LWRK)
C
C     contributions
C
      ERGEXC = ALLTOT - COUTOT
      EL2TOT = COUTOT
      EL2NON = COUNCL
      EL2CLS = EL2TOT - EL2NON
C
C     Print
C     =====
C
      ELNNON = EL1NON + EL2NON
      ELNCLS = EL1CLS + EL2CLS + POTNUC
      ELNTOT = ELNNON + ELNCLS
C
      ETOT = ELNTOT + ERGKIN + ERGEXC
C
      CALL HEADER('Decomposition of Coulomb energy',-1)
      WRITE (LUPRI,'(2X,17X,2A,/)') '        nonclassical',
     &             '         classical             total'
      WRITE (LUPRI,'(2X,A,3F18.8)')
     &       ' nuclear attraction',EL1NON,EL1CLS,EL1TOT
      WRITE (LUPRI,'(2X,A,3F18.8)')
     &       ' electron repulsion',EL2NON,EL2CLS,EL2TOT
      WRITE (LUPRI,'(2X,A,3F18.8)')
     &       ' nuclear repulsion ',D0,POTNUC,POTNUC
      WRITE (LUPRI,'(2X,A,3F18.8)')
     &       ' total Coulomb     ',ELNNON,ELNCLS,ELNTOT
      CALL HEADER('Decomposition of total electronic energy',-1)
      WRITE (LUPRI,'(2X,A,F18.8,9X,F7.2,A)')
     &       ' kinetic energy    ',ERGKIN, D100*ERGKIN/ETOT,' %'
      WRITE (LUPRI,'(2X,A,F18.8,9X,F7.2,A)')
     &       ' Coulomb energy    ',ELNTOT, D100*ELNTOT/ETOT,' %'
      WRITE (LUPRI,'(2X,A,F18.8,9X,F7.2,A)')
     &       ' exchange energy   ',ERGEXC, D100*ERGEXC/ETOT,' %'
      WRITE (LUPRI,'(2X,A,F18.8,9X,F7.2,A)')
     &       ' total energy      ',ETOT, D100,' %'
C
C
      RETURN
      END
C  /* Deck fmmdns */
      SUBROUTINE FMMDNS(DENMAT,WORK,LWORK,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
      LOGICAL NODC, NODV
      DIMENSION DENMAT(NNBASX), WORK(LWORK)
#include <inforb.h>
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMDNS','*',103)
C
      KFOCK = 1
      KLAST = KFOCK + NNBASX
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMDNS',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
      NODC = .FALSE.
      NODV = .TRUE.
      CALL DZERO(DENMAT,NNBASX)
      CALL DZERO(WORK(KFOCK),NNBASX)
      CALL DSOFSO(DENMAT,WORK(KFOCK),WORK(KLAST),LWRK,IPRINT,NODC,NODV)
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Density matrix in FMMDNS',-1)
         CALL OUTPAK(DENMAT,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck fmmint */
      SUBROUTINE FMMINT(SPHINT,LMAX,WORK,LWORK,IPRINT)
C
C     T. Helgaker
C
#include <implicit.h>
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxmom.h>
      LOGICAL      TOFILE, TRIMAT, DOINT(2,2)
      DIMENSION    SPHINT(NNBASX,*), WORK(LWORK)
      CHARACTER*10 LABINT(2*LMAX+1)
#include <inforb.h>
#include <nuclei.h>
#include <cbiher.h>
#include <orgcom.h>
C
      IF (IPRINT .GT. 2) CALL TITLER('Output from FMMINT','*',103)
C
      NCMPS = (LMAX + 1)**2
C
      KINTR = 1
      KINTA = KINTR + (NCMPS + IRAT - 1)/IRAT
      KLAST = KINTA + (NCMPS + IRAT - 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('FMMINT',' ',KLAST,LWORK)
      LWRK  = LWORK - KLAST + 1
C
      NCOMP = 1
      DOINT(1,1) = .TRUE.
      DOINT(2,1) = .FALSE.
      DOINT(1,2) = .FALSE.
      DOINT(2,2) = .FALSE.
      TRIANG = .TRUE.
      TOFILE = .FALSE.
      PROPRI = IPRINT .GT. 5
C
      KFREE = 1
      LFREE = LWRK
      IOFF = 0
      DO IORDER = 0, LMAX
         NCOMPS = 2*IORDER + 1
         CALL PR1IN1(WORK(KLAST),KFREE,LFREE,WORK(KINTR),WORK(KINTA),
     &               LABINT,'SPHMOML',IORDER,MPQUAD,TRIANG,PROPRI,
     &               IPRINT,SPHINT(1,IOFF+1),NCOMP,TOFILE,'TRIANG',
     &               DOINT,DUMMY,.FALSE.,DUMMY)
         IF (IPRINT .GT. 5)
     &      CALL FMMSPR(SPHINT(1,IOFF+1),LABINT,WORK(KINTR),NCOMPS)
         IOFF = IOFF + NCOMPS
      END DO
      RETURN
      END
C
C  /* Deck fmmspr */
      SUBROUTINE FMMSPR(SPHINT,LABINT,INTREP,NCOMPS)
C
C     T. Helgaker
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
      CHARACTER LABINT(*)*8
      DIMENSION SPHINT(NNBASX,*), INTREP(*)
#include <inforb.h>
C
      DO ICOMP = 1, NCOMPS
         CALL AROUND('Integrals of operator: '//LABINT(ICOMP))
         WRITE (LUPRI,'(A,I2)') ' Symmetry of operator:',
     &          INTREP(ICOMP) + 1
         CALL OUTPAK(SPHINT(1,ICOMP),NBAST,1,LUPRI)
      END DO
      RETURN
      END
C  /* Deck sphave */
      SUBROUTINE SPHAVE(SPHINT,DENMAT,LMAX,DIS,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxmom.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION SPHINT(NNBASX,*), DENMAT(NNBASX), DIS(*)
#include <ccom.h>
#include <primit.h>
#include <nuclei.h>
#include <shells.h>
#include <symmet.h>
#include <inforb.h>
#include <clsfmm.h>
C
      CALL HEADER('Output from SPHAVE',-1)
C
      LUINTM = -1
      CALL GPOPEN(LUINTM,'ABAFMM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUINTM
C
C     Nuclear information
C     ===================
C
      WRITE (LUINTM) NUCIND
      DO I = 1, NUCIND
         WRITE (LUINTM) CHARGE(I), CORD(1,I), CORD(2,I),  CORD(3,I)
      END DO
C
C     Electronic information
C     ======================
C
      ERFAC = ERFCIV(THRCLS)
      WRITE (LUINTM) NBAST, LMAX, THRCLS, ERFAC
C
      ICOMP = 0
      DO L = 0, LMAX
      MM = 0
      DO M = -L, L
         ICOMP = ICOMP + 1
         EXPVAL = D0
         IORBA = 0
         IJ = 0
         DO ISHELA = 1,KMAX
         DO ICOMPA = 1, KHKT(ISHELA)
            IORBA = IORBA + 1
            IORBB = 0
            DO ISHELB = 1, KMAX
            DO ICOMPB = 1, KHKT(ISHELB)
               IORBB = IORBB + 1
               IF (IORBB .LE. IORBA) THEN
                  IJ = IJ + 1
                  EXPA = PRIEXP(ISHELA)
                  EXPB = PRIEXP(ISHELB)
                  EXPP = EXPA + EXPB
                  PINV = D1/EXPP
                  PX = PINV*(EXPA*CENT(ISHELA,1,1)
     &                     + EXPB*CENT(ISHELB,1,1))
                  PY = PINV*(EXPA*CENT(ISHELA,2,1)
     &                     + EXPB*CENT(ISHELB,2,1))
                  PZ = PINV*(EXPA*CENT(ISHELA,3,1)
     &                     + EXPB*CENT(ISHELB,3,1))
                  DSPH = DENMAT(IJ)*SPHINT(IJ,ICOMP)
                  EXPVAL = EXPVAL + DSPH
                  WRITE (LUINTM) L, MM, ICOMP, IJ, EXPP, PX,PY,PZ, DSPH
               END IF
            END DO
            END DO
         END DO
         END DO
C
         WRITE (LUPRI,'(2X,A,2I3,2X,F12.6)')
     &      ' Expectation value of spherical harmonic ', L,M,EXPVAL
         IF (MM .GT. 0) THEN
            MM = -MM
         ELSE
            MM = -MM + 1
         END IF
      END DO
      END DO
C
      WRITE (LUINTM) -1
      CALL GPCLOSE(LUINTM,'KEEP')
C
      CALL FMMZER(DIS,NNBASX,IPRINT)
      RETURN
      END
C  /* Deck fmmorb */
      SUBROUTINE FMMORB(AOCOOR,AOEXP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
C
      DIMENSION AOCOOR(NBAST,3), AOEXP(NBAST)
C
#include <ccom.h>
#include <primit.h>
#include <nuclei.h>
#include <shells.h>
#include <symmet.h>
#include <inforb.h>
C
#include <ibtfun.h>
C
      IORBA = 0
      DO ISHELA = 1,KMAX
         KHKTA = KHKT(ISHELA)
         DO ICOMP = 1, KHKTA
            IORBA = IORBA + 1
            AOCOOR(IORBA,1) = CENT(ISHELA,1,1)
            AOCOOR(IORBA,2) = CENT(ISHELA,2,1)
            AOCOOR(IORBA,3) = CENT(ISHELA,3,1)
            AOEXP(IORBA) = PRIEXP(JSTRT(ISHELA) + NUMCF(ISHELA))
         END DO
      END DO
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from FMMORB',-1)
         DO I = 1, NBAST
            WRITE (LUPRI,'(2X,I5,2X,3F12.6,5X,F12.6)')
     &            I, (AOCOOR(I,J),J=1,3), AOEXP(I)
Chj: see FMMINI
C           WRITE (18,'(2X,I5,2X,3F24.16,2X,F20.12)')
C    &            I, (AOCOOR(I,J),J=1,3), AOEXP(I)
         END DO
      END IF
      RETURN
      END
C  /* Deck geterg */
      SUBROUTINE GETERG(ERGTOT,ERGCOU,DMAT,WORK,LWORK)
C
C     T. Helgaker
C
#include <implicit.h>
#include <iratdef.h>
      PARAMETER (D0 = 0.0D0, D4 = 4.0D0)
      DIMENSION DMAT(N2BASX), WORK(LWORK), MBAS(8)
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
C
#include <memint.h>
C
      ERGTOT = D0
      ERGCOU = D0
C
      LUINTA = -1
      CALL GPOPEN(LUINTA,'AOTWOINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
      CALL REWSPL(LUINTA)
      CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
      READ (LUINTA) MSYM, MBAS, LBUF, NIBUF, NBITS
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
C
      KINT  = 1
      KIINT = KINT + LBUF
      KIIN4 = KINT + LENINT
      KLAST = KIIN4 + 4*LBUF
      IF (KLAST .GT. LWORK) CALL STOPIT('GETERG','LUINTA',KLAST,LWORK)
C
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
    1 CONTINUE
         CALL READI(LUINTA,LENINT,WORK(KINT))
         CALL AOLAB4(WORK(KIINT),LBUF,NIBUF,NBITS,WORK(KIIN4),NINT)
         IF (NINT .EQ. 0)  GOTO 1
         IF (NINT .EQ. -1) GOTO 2
         CALL ADDERG(ERGTOT,ERGCOU,DMAT,WORK(KINT),WORK(KIIN4),NINT)
      GOTO 1
C
    2 CONTINUE
      CALL GPCLOSE(LUINTA,'DELETE')
      ERGTOT = D4*ERGTOT
      ERGCOU = D4*ERGCOU
C     WRITE (LUPRI,'(1X,A,2F12.6)')
C    &    ' Two-electron total and Coulomb energies: ',ERGTOT,ERGCOU
      RETURN
      END
C  /* Deck adderg */
      SUBROUTINE ADDERG(ERGTOT,ERGCOU,DMAT,BUF,IBUF,LENGTH)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.50D00, DP25 = 0.25D0)
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION DMAT(NBAST,NBAST), BUF(LENGTH), IBUF(4,LENGTH)
#include <ibtfun.h>
C
      DO 100 INT = 1, LENGTH
         DINT = BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DP5*DINT
         IF (R.EQ.S)              DINT = DP5*DINT
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DP5*DINT
         EINT = DP25*DINT
         ERGTOT = ERGTOT + DMAT(P,Q)*DMAT(R,S)*DINT
     &                   - DMAT(P,R)*DMAT(Q,S)*EINT
     &                   - DMAT(P,S)*DMAT(R,Q)*EINT
         ERGCOU = ERGCOU + DMAT(P,Q)*DMAT(R,S)*DINT
  100 CONTINUE
      RETURN
      END
C  /* Deck get1el */
      SUBROUTINE GET1EL(ERGONE,DMAT,HMAT,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.D0, D2 = 2.0D0)
      DIMENSION DMAT(NBAST,NBAST), HMAT(NNBAST)
#include <inforb.h>
      ERGONE = D0
      IJ = 0
      DO I = 1, NBAST
         DO J = 1, I - 1
            IJ = IJ + 1
            ERGONE = ERGONE + D2*DMAT(I,J)*HMAT(IJ)
         END DO
         IJ = IJ + 1
         ERGONE = ERGONE + DMAT(I,I)*HMAT(IJ)
      END DO
      RETURN
      END
C  /* Deck fmmzer */
      SUBROUTINE FMMZER(DIS,NDIM,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
      DIMENSION DIS(NDIM,8)
C
      LUINTM = -1
      CALL GPOPEN(LUINTM,'ABAFMM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUINTM
C
C     Nuclear information
C     ===================
C
      CALL HEADER('Nuclear information',-1)
      READ (LUINTM) NUCIND
      DO I = 1, NUCIND
         READ (LUINTM) CHRG, CRDX, CRDY, CRDZ
         WRITE (LUPRI,'(2X,4F18.6)') CHRG, CRDX, CRDY, CRDZ
      END DO
C
C     Electronic information
C     ======================
C
      CALL HEADER('Electronic information',-1)
      READ (LUINTM) NBAST, LMAX, THRCLS, ERFAC
C
      WRITE (LUPRI,'(2X,A,I12)')      ' NBAST  ', NBAST
      WRITE (LUPRI,'(2X,A,I12)')      ' LMAX   ', LMAX
      WRITE (LUPRI,'(2X,A,1P,E13.6)') ' THRCLS', THRCLS
      WRITE (LUPRI,'(2X,A,1P,F13.6)') ' ERFAC ', ERFAC
C
      DO L =  0, LMAX
      DO M = -L, L
         DO IORBA = 1, NBAST
         DO IORBB = 1, IORBA
             READ (LUINTM) LL, MM, ICOMP, IJ, EXPP, PX,PY,PZ, DSPH
             IF (L.EQ.0) THEN
               DIS(IJ,1) = PX
               DIS(IJ,2) = PY
               DIS(IJ,3) = PZ
               DIS(IJ,4) = ERFAC/SQRT(EXPP)
               DIS(IJ,5) = DSPH
             ELSE IF (L.EQ.1) THEN
               IF (M.EQ.-1) DIS(IJ,6) = DSPH
               IF (M.EQ. 0) DIS(IJ,7) = DSPH
               IF (M.EQ. 1) DIS(IJ,8) = DSPH
             END IF
         END DO
         END DO
      END DO
      END DO
C
      ERG0 = D0
      ERGP = D0
      ERGQ = D0
      DO I = 1, NBAST*(NBAST+1)/2
      DO J = 1, NBAST*(NBAST+1)/2
         RX = DIS(I,1)-DIS(J,1)
         RY = DIS(I,2)-DIS(J,2)
         RZ = DIS(I,3)-DIS(J,3)
         R  = SQRT(RX**2 + RY**2 + RZ**2)
         R3 = R**3
         IF (R .GT. DIS(I,4) + DIS(J,4)) THEN
            P00 =   DIS(I,5)
            P1M =   DIS(I,8)
            P10 =   DIS(I,6)
            P1P =   DIS(I,7)
            Q00 =   DIS(J,5)
            Q1M =   DIS(J,8)
            Q10 =   DIS(J,6)
            Q1P =   DIS(J,7)
            C0000 = P00*Q00/R
            C1M00 = RY*P1M*Q00/R3
            C1000 = RZ*P10*Q00/R3
            C1P00 = RX*P1P*Q00/R3
            C001M = RY*P00*Q1M/R3
            C0010 = RZ*P00*Q10/R3
            C001P = RX*P00*Q1P/R3
            ERG0  = ERG0 + C0000
            ERGP  = ERGP - C1M00 - C1000 - C1P00
            ERGQ  = ERGQ + C001M + C0010 + C001P
         END IF
      END DO
      END DO
      ERG0 = DP5*ERG0
      ERG1 = DP5*ERGP + DP5*ERGQ
      WRITE (LUPRI,'(/2X,A,3F15.8)')
     &      ' First-order classical interaction:',ERG0,ERG1,ERG0+ERG1
      RETURN
      END
C  /* Deck erfciv */
      FUNCTION ERFCIV(THR)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, THRMAX = 1.0D-15)
      DIMENSION ERFCI(0:15)
      DATA (ERFCI(I),I=0,15)/0.00000D0, 1.16309D0, 1.82139D0, 2.32675D0,
     &                       2.75106D0, 3.12341D0, 3.45891D0, 3.76656D0,
     &                       4.05224D0, 4.32001D0, 4.57282D0, 4.81292D0,
     &                       5.04203D0, 5.26151D0, 5.47248D0, 5.67585D0/
      ARG = MAX(MAX(ABS(THR),THRMAX),D0)
      ERFCIV = ERFCI(-NINT(DLOG10(ARG)))
      RETURN
      END
