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
#include "single.h"
C
C  /* Deck onesop */
      SUBROUTINE ONESOP(STHMAT,DENMAT,FOCMAT,FACINT,COORC,WORK,LWORK,
     &                  IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,
     &                  CORBZ0,DIFDIP,DIFQDP,SECDER,NATOMC,TOLOG,TOLS,
     &                  JSYMC,JCENTC,NCENTC,SIGNC,GEXP,NNBASX,FCM,TLMD,
     &                  HESSKE,HESSNA,HESFS2,HSOLT2,HSOLNN,NCLONE,PCM)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
C
#include "mxcent.h"
#include "pcmdef.h"
#include "pcm.h"
Clf#include "pcmlog.h"
C
      LOGICAL   SECDER, DIFDIP, PROPTY, NCLONE, DIFQDP, PCM
      DIMENSION WORK(LWORK)
      DIMENSION STHMAT(NNBASX,3), DENMAT(NNBASX), FOCMAT(NNBASX),
     &          FACINT(*), COORC(3,*), JSYMC(*),
     &          JCENTC(*), NCENTC(*), SIGNC(3,*),
     &          GEXP(*),   FCM(*), TLMD(*)
#ifndef PRG_DIRAC
      DIMENSION HESSKE(*), HESSNA(*), HESFS2(*), HSOLT2(*), HSOLNN(*)
#endif
#include "cbisol.h"
#include "onecom.h"
#include "expopt.h"
C
      CALL QENTER('ONESOP')
C
C     Memory allocation
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
      IF (PROPTY .AND. .NOT. EXPGRA) THEN
         NADER = 28*NATOMC*KCKTAB
      ELSE
         NADER = KCKTAB
      END IF
      CALL MEMGET('REAL',KSDER0, 2*KCKTAB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSDER1, 6*KCKTAB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSDER2,12*KCKTAB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KADER ,    NADER,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSINT0,   KCKTAB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDINT1, 9*KCKTAB,WORK,KFREE,LFREE)
      IF (DIFQDP) THEN
         CALL MEMGET('REAL',KSINT1, 3*KCKTAB,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KQINT1,18*KCKTAB,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KSINT1,        0,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KQINT1,        0,WORK,KFREE,LFREE)
      END IF
      CALL MEMGET('REAL',KDSHEL,   KHKTAB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFSHEL,   KHKTAB,WORK,KFREE,LFREE)
      IF (EXPGRA) THEN
         LXPDER = 4*KCKTAB*NUCA*NUCB
         CALL MEMGET('REAL',KXPDER,   LXPDER,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KXPDER,        0,WORK,KFREE,LFREE)
      END IF
      IF (SOLVNT) THEN
         IF (MAXDIF .EQ. 0) THEN
            CALL MEMGET('REAL',KRLMINT,  LMNTOT*KCKTAB,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KRLMTAB,              0,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KLMNO  ,              0,WORK,KFREE,LFREE)
         ELSE IF (MAXDIF .EQ. 1) THEN
            CALL MEMGET('REAL',KRLMINT,7*LMNTOT*KCKTAB,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KRLMTAB,              0,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KLMNO  ,              0,WORK,KFREE,LFREE)
         ELSE
            CALL MEMGET('REAL',KRLMINT,7*LMNTOT*KCKTAB,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KRLMTAB,      21*KCKTAB,WORK,KFREE,LFREE)
            CALL MEMGET('INTE',KLMNO  ,   3*LMNTOT + 1,WORK,KFREE,LFREE)
         END IF
      ELSE IF (PCM) THEN
         CALL MEMGET('REAL',KRLMINT,7*NTS*KCKTAB,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KRLMTAB,           0,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KLMNO  ,           0,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KRLMINT, 0,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KRLMTAB, 0,WORK,KFREE,LFREE)
         CALL MEMGET('INTE',KLMNO  , 0,WORK,KFREE,LFREE)
      END IF
C
      CALL ONESO1(STHMAT,WORK(KSDER0),WORK(KSDER1),WORK(KSDER2),
     &            WORK(KADER),WORK(KSINT0),WORK(KSINT1),WORK(KDINT1),
     &            WORK(KQINT1),DENMAT,FOCMAT,FACINT,COORC,
     &            WORK(KFREE),LFREE,IPRINT,PROPTY,MAXDIF,IDENB0,
     &            CORBX0,CORBY0,CORBZ0,DIFDIP,DIFQDP,
     &            SECDER,NATOMC,TOLOG,TOLS,JSYMC,JCENTC,NCENTC,SIGNC,
     &            GEXP,NNBASX,WORK(KDSHEL),WORK(KFSHEL),WORK(KRLMINT),
     &            WORK(KRLMTAB),FCM,TLMD,WORK(KLMNO),
     &            HESSKE,HESSNA,HESFS2,HSOLT2,HSOLNN,NCLONE,
     &            WORK(KXPDER),PCM)
      CALL MEMREL('ONESOP',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('ONESOP')
      RETURN
      END
C  /* Deck oneso1 */
      SUBROUTINE ONESO1(STHMAT,STDER0,STDER1,STDER2,ADER,SINT0,SINT1,
     &                  DINT1,QINT1,DENMAT,FOCMAT,FACINT,COORC,
     &                  WORK, LWORK,IPRINT,PROPTY,MAXDIF,IDENB0,
     &                  CORBX0,CORBY0,CORBZ0,DIFDIP,DIFQDP,SECDER,
     &                  NATOMC,TOLOG,TOLS,JSYMC,JCENTC,NCENTC,
     &                  SIGNC,GEXP,NNBASX,DSHELL,FSHELL,
     &                  RLMINT,RLMTAB,FCM,TLMD,LMNO,
     &                  HESSKE,HESSNA,HESFS2,HSOLT2,HSOLNN,
     &                  NCLONE,EXPDER,PCM)
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "iratdef.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (D1 = 1.0D0, DP5 = 0.5D0)
C
      LOGICAL SECDER, DIFDIP, PROPTY, FULMAT, ANTI, NCLONE, DIFQDP, PCM
      DIMENSION WORK(LWORK)
#ifdef PRG_DIRAC
      PARAMETER(NDER1=6,NDER2=6)
C     ... note that STDER2 is not implemented yet, and will require more than 6
#else
      PARAMETER(NDER1=3,NDER2=6)
#endif
      DIMENSION STHMAT(NNBASX,3), STDER0(KCKTAB,2),
     &          STDER1(KCKTAB,NDER1,2), STDER2(KCKTAB,NDER2,2), ADER(*),
     &          EXPDER(*), DSHELL(KHKTAB), FSHELL(KHKTAB),
     &          DENMAT(NNBASX),FOCMAT(NNBASX),
     &          SINT0(KCKTAB), SINT1(KCKTAB,3), DINT1(KCKTAB,3,3),
     &          QINT1(KCKTAB,3,6),
     &          FACINT(*), COORC(3,*), JSYMC(*),
     &          JCENTC(*), NCENTC(*), SIGNC(3,*), GEXP(*),
     &          RLMINT(*), RLMTAB(*), FCM(*), TLMD(*), LMNO(*)
      DIMENSION HESSKE(*), HESSNA(*), HESFS2(*), HSOLT2(*), HSOLNN(*)
#include "ccom.h"
#include "onecom.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbisol.h"
#include "pcmdef.h"
#include "pcm.h"
#include "expopt.h"
C
#include "energy.h"
#include "ibtfun.h"
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
      IF (IPRINT .GT. 4) CALL TITLER('Output from ONESO1','*',103)
C
C     *****************************************
C     ***** Loop over symmetry operations *****
C     *****************************************
C
      IDENB = IDENB0 - KHKTB
      DO 100 ISYMOP = 0, MAXOPR
      IF(IBTAND(ISYMOP,MAB) .EQ. 0) THEN
         IDENB  = IDENB + KHKTB
         ICENTB = NUCNUM(NCENTB,ISYMOP+1)
         ONECEN = ICENTA .EQ. ICENTB
         SIGNBX = PT(IBTAND(ISYMAX(1,1),ISYMOP))
         SIGNBY = PT(IBTAND(ISYMAX(2,1),ISYMOP))
         SIGNBZ = PT(IBTAND(ISYMAX(3,1),ISYMOP))
         CORBX  = SIGNBX*CORBX0
         CORBY  = SIGNBY*CORBY0
         CORBZ  = SIGNBZ*CORBZ0
         IF (IPRINT .GE. 05) WRITE (LUPRI, 1010) ISYMOP
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI,'(A,F12.6)') ' CORBX ', CORBX
            WRITE (LUPRI,'(A,F12.6)') ' CORBY ', CORBY
            WRITE (LUPRI,'(A,F12.6)') ' CORBZ ', CORBZ
         END IF
C
C        **********************************************
C        ***** Calculation of Cartesian integrals *****
C        **********************************************
C
         CALL ONEPRM(STDER0,STDER1,STDER2,ADER,SINT0,SINT1,DINT1,QINT1,
     &               RLMINT,RLMTAB,FCM,WORK,LWORK,IPRINT,PROPTY,MAXDIF,
     &               NATOMC,TOLOG,TOLS,SECDER,DIFDIP,DIFQDP,
     &               FACINT,COORC,GEXP,JCENTC,NCENTC,NCLONE,EXPDER,PCM)
C
C        *************************************************
C        ***** Transform to spherical harmonic basis *****
C        *************************************************
C
         IF (SPHRAB) THEN
            IF (EXPGRA) THEN
               CALL SPHRM1(EXPDER,EXPDER,4*NUCA*NUCB,WORK,LWORK,
     &                     PROPTY,IPRINT)
            ELSE IF (PROPTY) THEN
               IF (ONECEN) THEN
                  CALL SPHRM1(ADER,ADER,10*NATOMC,WORK,LWORK,PROPTY,
     &                        IPRINT)
               ELSE
                  CALL SPHRM1(STDER0,STDER0,2,WORK,LWORK,PROPTY,IPRINT)
                  CALL SPHRM1(STDER1,STDER1,2*NDER1,
     &                        WORK,LWORK,PROPTY,IPRINT)
                  IF (SECDER) THEN
                     CALL SPHRM1(STDER2,STDER2,2*NDER2,
     &                           WORK,LWORK,PROPTY,IPRINT)
                  END IF
                  CALL SPHRM1(ADER,ADER,28*NATOMC,WORK,LWORK,PROPTY,
     &                        IPRINT)
                  IF (DIFDIP) THEN
                     CALL SPHRM1(DINT1,DINT1,9,WORK,LWORK,PROPTY,IPRINT)
                  END IF
                  IF (DIFQDP) THEN
                    CALL SPHRM1(QINT1,QINT1,18,WORK,LWORK,PROPTY,IPRINT)
                  END IF
               END IF
               CALL SPHRM1(SINT0,SINT0,1,WORK,LWORK,PROPTY,IPRINT)
               IF (DIFQDP) THEN
                  CALL SPHRM1(SINT1,SINT1,3,WORK,LWORK,PROPTY,IPRINT)
               END IF
               IF (SOLVNT) THEN
                  CALL SPHRM1(RLMINT,RLMINT,7*LMNTOT,WORK,LWORK,PROPTY,
     &                        IPRINT)
                  IF (SECDER) THEN
                     CALL SPHRM1(RLMTAB,RLMTAB,21,WORK,LWORK,PROPTY,
     &                           IPRINT)
                  END IF
               ELSE IF (PCM) THEN
                  CALL SPHRM1(RLMINT,RLMINT,7*NTS,WORK,LWORK,PROPTY,
     &                        IPRINT)
               END IF
            ELSE
               CALL SPHRM1(STDER0,STDER0,2,WORK,LWORK,PROPTY,IPRINT)
               CALL SPHRM1(ADER,ADER,1,WORK,LWORK,PROPTY,IPRINT)
               IF (SOLVNT) THEN
                  CALL SPHRM1(RLMINT,RLMINT,LMNTOT,WORK,LWORK,PROPTY,
     &                        IPRINT)
               END IF
            END IF
         END IF
C
C        ******************************
C        ***** Expectation values *****
C        ******************************
C
         IF (PROPTY) THEN
C
C           One-electron Hamiltonian integrals and reorthonormalization
C           ===========================================================
C
C           Collect density and Fock elements
C
            MAXCMP = 0
            DO 200 IORBA = IDENA + 1, IDENA + KHKTA
               DO 250 IORBB = IDENB + 1, IDENB + KHKTB
                  FAC = D1
                  IF (ONECEN .AND. LDIAG .AND. IORBB.NE.IORBA) FAC = DP5
                  MAXCMP = MAXCMP + 1
                  IORBAB = ITRI(IORBA,IORBB)
                  DSHELL(MAXCMP) = FAC*DENMAT(IORBAB)
                  FSHELL(MAXCMP) = FAC*FOCMAT(IORBAB)
  250          CONTINUE
  200       CONTINUE
C
            IF (EXPGRA) THEN
C
C              Orbital-exponent gradient 
C
               CALL AVEEXP(EXPDER,DSHELL,FSHELL)
            ELSE IF (ONECEN) THEN
C
C              Nuclear attraction
C
               CALL AVENA1(ADER,NATOMC,SECDER,NCENTC,JCENTC,MAXCMP,
     &                     JSYMC,SIGNC,DSHELL,HESSNA)
            ELSE
C
C              Kinetic energy and reorthonormalization
C
               CALL AVEKFS(STDER0,STDER1,STDER2,ISYMOP,MAXCMP,SECDER,
     &                     DSHELL,FSHELL,HESSKE,HESFS2)
C
C              Nuclear attraction
C
               CALL AVENA2(ADER,NATOMC,ISYMOP,SECDER,NCENTC,MAXCMP,
     &                     JCENTC,JSYMC,SIGNC,DSHELL,HESSNA)
            END IF
C
C           Dipole gradient
C           ===============
C
            IF (DIFDIP) THEN
               CALL AVEDIP(SINT0,DINT1,ISYMOP,DSHELL,MAXCMP)
            END IF
C
#ifndef PRG_DIRAC
C           Second moment gradient
C           ======================
C
            IF (DIFQDP) THEN
               CALL AVEQDP(SINT1,QINT1,ISYMOP,DSHELL,MAXCMP)
            END IF
C
C           Solvent contributions
C           =====================
            IF (SOLVNT) THEN
                CALL AVESOL(RLMINT,RLMTAB,FCM,MAXDIF,DIFDIP,ISYMOP,
     &                      MAXCMP,DSHELL,TLMD,LMNO,IPRINT,HSOLT2)
             ELSE IF (PCM) THEN
                CALL AVEPCM(RLMINT,MAXCMP,DSHELL,IPRINT)
            END IF
#endif
         END IF
C
C        *******************************************
C        ***** Transform integrals to SO basis *****
C        *******************************************
C
         FULMAT = .TRUE.
         ANTI   = .FALSE.
C
C        Overlap integrals
C        =================
C
         CALL SYM1S(STDER0(1,1),STHMAT(1,1),ISYMOP,MULA,MULB,NHKTA,
     &              NHKTB,KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
     &              IPRINT)
C
C        Kinetic energy integrals
C        ========================
C
         CALL SYM1S(STDER0(1,2),STHMAT(1,3),ISYMOP,MULA,MULB,NHKTA,
     &              NHKTB,KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
     &              IPRINT)
C
C        Nuclear attraction integrals
C        ============================
C
         CALL SYM1S(ADER,STHMAT(1,2),ISYMOP,MULA,MULB,NHKTA,NHKTB,
     &              KHKTA,KHKTB,HKAB,LDIAG,FULMAT,DUMMY,IDUMMY,
     &              IPRINT)
C
C        **************************************************
C        ***** Write differentiated integrals on file *****
C        **************************************************
C
         IF (PROPTY .AND. (SECDER .OR. DIFDIP .OR. DIFQDP)) THEN
            CALL WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
     &                  ISYMOP,NATOMC,IPRINT,PCM)
         END IF
      END IF
  100 CONTINUE
      RETURN
 1010 FORMAT (//,2X,'***************************************',
     *         /,2X,'******** Symmetry operation ',I2,' ********',
     *         /,2X,'***************************************',/)
      END
C  /* Deck wd1sym */
      SUBROUTINE WD1SYM(STDER1,ADER,RLMINT,FCM,WORK,LWORK,JSYMC,JCENTC,
     &                  ISYMOP,NATOMC,IPRINT,PCM)
#include "implicit.h"
#include "priunit.h"
      DIMENSION STDER1(KCKTAB,3,2), ADER(*), JSYMC(*), JCENTC(*),
     &          RLMINT(*), FCM(*), WORK(LWORK)
      LOGICAL PCM
#include "onecom.h"
#include "mxcent.h"
#include "pcmdef.h"
#include "pcm.h"
clf#include "pcmlog.h"
C
#ifdef PRG_DIRAC
C     HJAaj Nov 2002: this routine is not programmed for Dirac (yet)
      CALL QUIT('WD1SYM called in Dirac, but not programmed for Dirac')
#endif
      KOMAT = 1
      IF (PCM) THEN
         KQTOT = KOMAT + 6*KCKTAB
         KLAST = KQTOT + NTSIRR
      ELSE
         KLAST  = KOMAT + 6*KCKTAB
         KQTOT  = KLAST
      ENDIF
      IF (KLAST .GT. LWORK) CALL STOPIT('WD1SYM',' ',KLAST,LWORK)
      CALL WD1SY1(STDER1,ADER,RLMINT,FCM,WORK(KOMAT),JSYMC,JCENTC,
     &            ISYMOP,NATOMC,IPRINT,WORK(KQTOT),PCM)
      RETURN
      END
C  /* Deck wd1sy1 */
      SUBROUTINE WD1SY1(STDER1,ADER,RLMINT,FCM,OMAT,JSYMC,JCENTC,ISYMOP,
     &                  NATOMC,IPRINT,QTOT,PCM)
#include "implicit.h"
#include "priunit.h"
      DIMENSION STDER1(KCKTAB,3,2), ADER(KCKTAB,NATOMC,*), JSYMC(*),
     &          JCENTC(*), OMAT(KCKTAB,3,2), RLMINT(KCKTAB,LMNTOT,7),
     &          FCM(*),QTOT(*)
      LOGICAL PCM
#include "cbisol.h"
#include "onecom.h"
#include "ader.h"
#include "mxcent.h"
#include "pcmdef.h"
#include "pcm.h"
      IF (IPRINT .GT. 5) THEN
         NADER = 28*NATOMC
         IF (ONECEN) NADER = 10*NATOMC
         CALL HEADER('First derivative overlap matrix',-1)
         CALL OUTPUT(STDER1(1,1,1),1,KHKTAB,1,3,KCKTAB,3,1,LUPRI)
         CALL HEADER('First derivative kinetic energy matrix',-1)
         CALL OUTPUT(STDER1(1,1,2),1,KHKTAB,1,3,KCKTAB,3,1,LUPRI)
         CALL HEADER('ADER',-1)
         CALL OUTPUT(ADER,1,KHKTAB,1,NADER,KCKTAB,NADER,1,LUPRI)
      END IF
C
C     Write differentiated integrals on file
C
C     Overlap and kinetic energy integrals
C     ====================================
C
      IF (.NOT.ONECEN) THEN
         CALL DCOPY(6*KCKTAB,STDER1,1,OMAT,1)
         CALL DRSYM1(OMAT(1,1,1),OMAT(1,1,2),NCENTA,NCENTB,ISYMOP,MULA,
     &               MULB,NHKTA,NHKTB,KHKTA,KHKTB,KHKTAB,KCKTAB,HKAB,
     &               LDIAG,IPRINT)
      END IF
C
C     Nuclear attraction integrals
C     ============================
C
      DO 200 IATOMC = 1, NATOMC
         DO 300 ICMPAB = 1, KHKTAB
            OMAT(ICMPAB,1,1) = ADER(ICMPAB,IATOMC,IA000X)
            OMAT(ICMPAB,2,1) = ADER(ICMPAB,IATOMC,IA000Y)
            OMAT(ICMPAB,3,1) = ADER(ICMPAB,IATOMC,IA000Z)
            IF (.NOT.ONECEN) THEN
               OMAT(ICMPAB,1,2) = ADER(ICMPAB,IATOMC,IA0X00)
               OMAT(ICMPAB,2,2) = ADER(ICMPAB,IATOMC,IA0Y00)
               OMAT(ICMPAB,3,2) = ADER(ICMPAB,IATOMC,IA0Z00)
            END IF
  300    CONTINUE
         ISYMC  = JSYMC(IATOMC)
         ICENTC = JCENTC(IATOMC)
         CALL DASYM1(OMAT(1,1,1),OMAT(1,1,2),ONECEN,NCENTA,NCENTB,
     &               ICENTC,ISYMOP,ISYMC,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &               KHKTB,KHKTAB,KCKTAB,HKAB,LDIAG,IPRINT)
  200 CONTINUE
C
C     Solvent contribution
C     ====================
C
      IF (SOLVNT) THEN
        DO 400 ICMPAB = 1, KHKTAB
          OMAT(ICMPAB,1,1)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,5),KCKTAB)
          OMAT(ICMPAB,2,1)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,6),KCKTAB)
          OMAT(ICMPAB,3,1)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,7),KCKTAB)
          OMAT(ICMPAB,1,2)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,2),KCKTAB)
          OMAT(ICMPAB,2,2)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,3),KCKTAB)
          OMAT(ICMPAB,3,2)=DDOT(LMNTOT,FCM,1,RLMINT(ICMPAB,1,4),KCKTAB)
  400   CONTINUE
        CALL DASYM1(OMAT(1,1,1),OMAT(1,1,2),ONECEN,NCENTA,NCENTB,
     &              NCNTCV,ISYMOP,0,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &              KHKTB,KHKTAB,KCKTAB,HKAB,LDIAG,IPRINT)
      ELSE IF (PCM) THEN
#ifdef PRG_DIRAC
         call quit('WD1SY1: PCM not implemented in DIRAC')
#else
C
CLF Do we have to multiply this by the multiplicity?
C
         DO ITS = 1, NTSIRR
            QTOT(ITS) = QSE(ITS) + QSN(ITS)
         END DO
         DO 410 ICMPAB = 1, KHKTAB
          OMAT(ICMPAB,1,1)=DDOT(NTSIRR,QTOT,1,RLMINT(ICMPAB,1,5),KCKTAB)
          OMAT(ICMPAB,2,1)=DDOT(NTSIRR,QTOT,1,RLMINT(ICMPAB,1,6),KCKTAB)
          OMAT(ICMPAB,3,1)=DDOT(NTSIRR,QTOT,1,RLMINT(ICMPAB,1,7),KCKTAB)
          OMAT(ICMPAB,1,2)=DDOT(NTSIRR,QTOT,1,RLMINT(ICMPAB,1,2),KCKTAB)
          OMAT(ICMPAB,2,2)=DDOT(NTSIRR,QTOT,1,RLMINT(ICMPAB,1,3),KCKTAB)
          OMAT(ICMPAB,3,2)=DDOT(NTSIRR,QTOT,1,RLMINT(ICMPAB,1,4),KCKTAB)
 410     CONTINUE
         CALL DASYM1(OMAT(1,1,1),OMAT(1,1,2),ONECEN,NCENTA,NCENTB,
     &        NCNTCV,ISYMOP,0,MULA,MULB,NHKTA,NHKTB,KHKTA,
     &        KHKTB,KHKTAB,KCKTAB,HKAB,LDIAG,IPRINT)
#endif
      END IF
      RETURN
      END
C  /* Deck dsym1 */
      SUBROUTINE DSYM1(DENMAT,FOCMAT,DSO,FSO,DOFOCK,NBAST,IPRINT)
C
C     Take packed density matrix in symmetry orbital basis (NNBAST)
C     and generate density matrix (NNBASX) over distinct pairs of AOs
C
C                                          880418  PRT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      DIMENSION DSO(*), FSO(*), DENMAT(*), FOCMAT(*)
      LOGICAL DOFOCK
#include "shells.h"
#include "pincom.h"
#include "symmet.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) CALL HEADER('Subroutine DSYM1',-1)
C
C     Loop over all irreps in molecule
C
      NNBASX = NBAST*(NBAST + 1)/2
      CALL DZERO(DENMAT,NNBASX)
      IF (DOFOCK) CALL DZERO(FOCMAT,NNBASX)

      ISOFF = 0
      ISTR  = 1
      DO 100 IREP = 0, MAXREP
         NAOSI = NAOS(IREP+1)
         IF (NAOSI .EQ. 0) GOTO 110
         DO 200 I = ISTR,ISTR + NAOSI - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            IOFF = KSTRT(IA)
            MULA = ISTBAO(IA)
            INDA = IOFF + NA
            DO 300 J = ISTR,I
               IB     = IBTAND(IBTSHR(IPIND(J),16),65535)
               NB     = IBTAND(IBTSHR(IPIND(J), 8),  255)
               JOFF   = KSTRT(IB)
               NHKTB  = NHKT(IB)
               KHKTB  = KHKT(IB)
               MULB   = ISTBAO(IB)
               MAB    = IBTOR(MULA,MULB)
               KAB    = IBTAND(MULA,MULB)
               HKAB   = FMULT(KAB)
               ISOFF  = ISOFF + 1
               DSYMIJ = DSO(ISOFF)
               IF (DOFOCK) FSYMIJ = FSO(ISOFF)
               INDB   = JOFF + NB - KHKTB
               DO 400 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MAB) .NE. 0) GOTO 400
                  INDB = INDB + KHKTB
C
C                 Weight and parity factor
C
                  FAC = HKAB*
     *                  PT(IBTAND(ISYMOP,IBTXOR(IREP,ISYMAO(NHKTB,NB))))
                  INDM = MAX(INDA,INDB)
                  IND  = (INDM*(INDM - 3))/2 + INDA + INDB
                  DENMAT(IND) = DENMAT(IND) + FAC*DSYMIJ
                  IF (DOFOCK) FOCMAT(IND) = FOCMAT(IND) + FAC*FSYMIJ
400            CONTINUE
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NAOSI
100   CONTINUE
      IF (ISTR - 1 .NE. NBAST) THEN
         CALL QENTER('DSYM1  ')
         WRITE (LUPRI,*)
     &      'ERROR in DSYM1: NBAST =',NBAST,', sum NAOS =',ISTR-1
         CALL QUIT('parameter NBAST does not agree with NAOS array')
      END IF
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Total density matrix in DSYM1 '//
     &     '(sym. distinct AO basis)',-1)
         CALL OUTPAK(DENMAT,NBAST,1,LUPRI)
         IF (DOFOCK) THEN
            CALL HEADER('Total Fock matrix in DSYM1 '//
     &           '(sym. distinct AO basis)',-1)
            CALL OUTPAK(FOCMAT,NBAST,1,LUPRI)
         END IF
      END IF
      RETURN
      END
C  /* Deck drsym1 */
      SUBROUTINE DRSYM1(SAODER,TAODER,ICENTA,ICENTB,ISYMOP,MULA,MULB,
     &                  NHKTA,NHKTB,KHKTA,KHKTB,KHKTAB,KCKTAB,HKAB,
     &                  LDIAG,IPRINT)
C
C     Arrange calculation of symmetry-adapted integral derivatives
C     from distinct AO integral derivatives
C                                                PRT & TUH  880428
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "csym1.h"
#include "ccom.h"
      DIMENSION SAODER(KCKTAB,*), TAODER(KCKTAB,*)
      LOGICAL LDIAG, FULMAT, ANTI
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DRSYM1',-1)
         WRITE (LUPRI,'(A,2I5)') ' ICENTA/B ', ICENTA, ICENTB
         WRITE (LUPRI,'(A, I5)') ' ISYMOP   ', ISYMOP
         WRITE (LUPRI,'(A,2I5)') ' NHKTA/B  ', NHKTA, NHKTB
         WRITE (LUPRI,'(A,2I5)') ' KHKTA/B  ', KHKTA, KHKTB
         WRITE (LUPRI,'(A,2I5)') ' MULA/B   ', MULA, MULB
         WRITE (LUPRI,'(A,F12.6)') ' HKAB   ', HKAB
         WRITE (LUPRI,'(A, L5)') ' LDIAG    ', LDIAG
      END IF
      FULMAT = .FALSE.
      ANTI   = .FALSE.
      NMATS = 3*NUCDEP*(MAXREP+1)
      DO 10 ICL = 1,2
C
C        Determine factors to account for use of transl. invariance
C
         IF (ICL .EQ. 1) THEN
            JCENT = ICENTA
         ELSE
            JCENT = ICENTB
         ENDIF
         MULJ = ISTBNU(JCENT)
C
C     Run over Cartesian directions
C
      DO 20 JCDIR = 1,3
         ISYTYJ = ISYMAX(JCDIR,1)
C
C     Run over irrep's of the differentiation operator
C
      DO 30 IREPD = 0,MAXREP
         IF (IBTAND(MULJ,IBTXOR(IREPD,ISYTYJ)) .EQ. 0) THEN
            IF (ICL .EQ. 1) THEN
               FAC = HKAB
            ELSE
               FAC = - HKAB*PT(IBTAND(ISYTYJ,ISYMOP))
     *                     *PT(IBTAND(IREPD, ISYMOP))
            END IF
            IF (IPRINT .GT. 20) THEN
               WRITE (LUPRI, '(/A,3I5)') ' ICL, JCDIR, IREPD ',
     *                                     ICL, JCDIR, IREPD
               WRITE (LUPRI, '(A,2I5)') ' MULJ, ISYTYJ ', MULJ, ISYTYJ
               WRITE (LUPRI, '(A,F12.6)') ' FAC ', FAC
            END IF
            IMAT0 = (IPTCNT(3*(JCENT-1)+JCDIR,IREPD,1)-1)*(MAXREP+1)+1
            IF (IMAT0 .GE. 2**10) THEN
               WRITE (LUPRI,*) 'Not enough space in label to store '//
     &              'derivative direction of one-electron integrals'
               WRITE (LUPRI,*) 'Please contact dalton-admin@kjemi.'//
     &              'uio.no for assistance'
               CALL QUIT('Not enough space in label for diff.direction')
            END IF
            IF (IREPD .EQ. 0) THEN
C
C              Overlap matrix - totally symmetric perturbation
C
               CALL SYM1S(SAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                    NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                    THRS,IMAT0,IPRINT)
C
C              Kinetic energy  - totally symmetric perturbation
C
               IMAT0 = NMATS + IMAT0
               CALL SYM1S(TAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                    NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                    THRS,IMAT0,IPRINT)
            ELSE
C
C              Overlap matrix - non-symmetric perturbation
C
               CALL SYM1N(SAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                    MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                    FULMAT,ANTI,THRS,IMAT0,IPRINT)
C
C              Kinetic energy  - non-symmetric perturbation
C
               IMAT0 = NMATS + IMAT0
               CALL SYM1N(TAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                    MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                    FULMAT,ANTI,THRS,IMAT0,IPRINT)
            END IF
         END IF
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
C  /* Deck dasym1 */
      SUBROUTINE DASYM1(CAODER,AAODER,ONECEN,ICENTA,ICENTB,ICENTC,
     &                  ISYMOP,JSYMOP,MULA,MULB,NHKTA,NHKTB,
     &                  KHKTA,KHKTB,KHKTAB,KCKTAB,HKAB,LDIAG,
     &                  IPRINT)
C
C     Arrange calculation of symmetry-adapted integral derivatives
C     from distinct AO nuclear attraction integral derivatives
C                                                PRT & TUH  880502
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      LOGICAL LDIAG, ONECEN
      DIMENSION CAODER(KCKTAB,*), AAODER(KCKTAB,*)
      LOGICAL FULMAT, ANTI
#include "ccom.h"
#include "csym1.h"
#include "nuclei.h"
#include "symmet.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine DASYM1',-1)
         WRITE (LUPRI,'(A,2I5)') ' ICENTA/B ', ICENTA, ICENTB
         WRITE (LUPRI,'(A, I5)') ' ISYMOP   ', ISYMOP
         WRITE (LUPRI,'(A, I5)') ' ICENTC   ', ICENTC
         WRITE (LUPRI,'(A, I5)') ' JSYMOP   ', JSYMOP
         WRITE (LUPRI,'(A,2I5)') ' NHKTA/B  ', NHKTA, NHKTB
         WRITE (LUPRI,'(A,2I5)') ' KHKTA/B  ', KHKTA, KHKTB
         WRITE (LUPRI,'(A,2I5)') ' MULA/B   ', MULA, MULB
         WRITE (LUPRI,'(A,F12.6)') ' HKAB   ', HKAB
         WRITE (LUPRI,'(A, L5)') ' LDIAG    ', LDIAG
         WRITE (LUPRI,'(A, L5)') ' ONECEN   ', ONECEN
      END IF
      FULMAT = .FALSE.
      ANTI   = .FALSE.
      NMATS = 3*NUCDEP*(MAXREP+1)
      IF (ONECEN) THEN
         ICLMX = 2
      ELSE
         ICLMX = 3
      ENDIF
      DO 10 ICL = 1,ICLMX
C
C        Determine factors to account for use of transl. invariance
C
         IF (ICL .EQ. 1) THEN
            JCENT = ICENTC
         ELSE IF (ICL .EQ. 2) THEN
            JCENT = ICENTA
         ELSE
            JCENT = ICENTB
         ENDIF
         MULJ = ISTBNU(JCENT)
C
C     Run over Cartesian directions
C
      DO 20 JCDIR = 1,3
         ISYTYJ = ISYMAX(JCDIR,1)
         IF (ICL .EQ. 3) THEN
            DO 50 I = 1, KHKTAB
               CAODER(I,JCDIR) = CAODER(I,JCDIR) + AAODER(I,JCDIR)
50          CONTINUE
         ENDIF
C
C     Run over irreps of the differentiation operator
C
      DO 30 IREPD = 0,MAXREP
         IF (IBTAND(MULJ,IBTXOR(IREPD,ISYTYJ)) .EQ. 0) THEN
            IF (ONECEN) THEN
               IF (ICL .EQ. 1) THEN
                  FAC = HKAB*PT(IBTAND(ISYTYJ,JSYMOP))
     *                      *PT(IBTAND(IREPD ,JSYMOP))
               ELSE
                  FAC = - HKAB
               ENDIF
            ELSE
               IF (ICL .EQ. 1) THEN
                  FAC = HKAB*PT(IBTAND(ISYTYJ,JSYMOP))
     *                      *PT(IBTAND(IREPD ,JSYMOP))
               ELSE IF (ICL .EQ. 2) THEN
                  FAC = HKAB
               ELSE
                  FAC = - HKAB*PT(IBTAND(ISYTYJ,ISYMOP))
     *                        *PT(IBTAND(IREPD, ISYMOP))
               END IF
            END IF
            IF (IPRINT .GT. 20) THEN
               WRITE (LUPRI, '(/A,3I5)') ' ICL, JCDIR, IREPD ',
     *                                     ICL, JCDIR, IREPD
               WRITE (LUPRI, '(A,2I5)') ' MULJ, ISYTYJ ', MULJ, ISYTYJ
               WRITE (LUPRI, '(A,F12.6)') ' FAC ', FAC
            END IF
            IMAT0 = NMATS
     *            + (IPTCNT(3*(JCENT-1)+JCDIR,IREPD,1)-1)*(MAXREP+1)+1
            IF (IREPD .EQ. 0) THEN
               IF (ICL .EQ. 2 .AND. .NOT.ONECEN) THEN
                  CALL SYM1S(AAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                       NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                       THRS,IMAT0,IPRINT)
               ELSE
                  CALL SYM1S(CAODER(1,JCDIR),DUMMY,ISYMOP,MULA,MULB,
     *                       NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,FULMAT,
     *                       THRS,IMAT0,IPRINT)
               END IF
            ELSE
               IF (ICL .EQ. 2 .AND. .NOT.ONECEN) THEN
                  CALL SYM1N(AAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                       MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                       FULMAT,ANTI,THRS,IMAT0,IPRINT)
               ELSE
                  CALL SYM1N(CAODER(1,JCDIR),DUMMY,IREPD,ISYMOP,MULA,
     *                       MULB,NHKTA,NHKTB,KHKTA,KHKTB,FAC,LDIAG,
     *                       FULMAT,ANTI,THRS,IMAT0,IPRINT)
               END IF
            END IF
         END IF
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
C  /* Deck sym1s */
      SUBROUTINE SYM1S(AO,SO,KB,MULA,MULB,NHKTA,NHKTB,KHKTA,KHKTB,
     *                 HKAB,LDIAG,FULMAT,THRESH,IMAT0,IPRINT)
C
C     Take block of distinct AO integral (derivatives) and
C     generate symmetrized contributions to SO integral
C     (derivatives) for the totally symmetric case
C                                          880407  PRT
C     Modified tuh 880819
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (IREPO = 0)
      LOGICAL LDIAG, FULMAT
      DIMENSION AO(*), SO(*)
#ifdef PRG_DIRAC
      PARAMETER (LUITMP = 48)
#else
#include "inftap.h"
#endif
#include "csym1.h"
#include "symmet.h"
#include "symind.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine SYM1S',-1)
         WRITE (LUPRI, '(A, I5)') ' KB ', KB
         WRITE (LUPRI, '(A,2I5)') ' NHKTA/B ', NHKTA, NHKTB
      END IF
C
C     Loop over all irreps in molecule
C
      DO 100 IREP = 0, MAXREP
         IF (FULMAT) THEN
            INDOFF = NPARSU(IREP + 1)
         ELSE
            IMAT = IMAT0 + IREP
         END IF
C
C        Loop over AOs which are of symmetry IREP in stabilizer MULA
C
         DO 200 NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
            NAT = KHKTB*(NA - 1)
            IF (LDIAG) THEN
               KHKTBB = NA
            ELSE
               KHKTBB = KHKTB
            END IF
            DO 300 NB = 1,KHKTBB
            IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
C
C              Weight and parity factor
C
               FAC = HKAB*PT(IBTAND(KB,IBTXOR(IREP,ISYMAO(NHKTB,NB))))
C
C              Locate SO integrals to which AO's contribute
C
               INDA = INDFA(IREP + 1,NA)
               INDB = INDFB(IREP + 1,NB)
               INDM = MAX(INDA,INDB)
               RINT = FAC*AO(NAT+NB)
               IF (FULMAT) THEN
                  IND  = INDOFF + (INDM*(INDM - 3))/2 + INDA + INDB
                  SO(IND) = SO(IND) + RINT
               ELSE
                  IF (ABS(RINT) .GT. THRESH) THEN
                     IND  = (INDM*(INDM - 3))/2 + INDA + INDB
                     INDMAX = MAX(IND,INDMAX)
                     LABEL = IND*2**10 + IMAT
                     IF (IPRINT .GT. 20) THEN
                       WRITE (LUPRI,'(A,F12.6,2I3,2I2,I4,I2,I5)')
     *                     'SYM1S - NA/B,IREPA/B,IND,IREPO,IMAT',
     *                     RINT, NA, NB, IREP, IREP, IND, IREPO,
     *                     IMAT
                     END IF
                     IF (LENGTH .EQ. 600) THEN
                        WRITE (LUITMP) BUF, IBUF, LENGTH
                        IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *                     ' Buffer of length',LENGTH,
     *                     ' has been written in SYM1S.'
                        LENGTH = 0
                     ENDIF
                     LENGTH = LENGTH + 1
                     BUF (LENGTH) = RINT
                     IBUF(LENGTH) = LABEL
                  ENDIF
               ENDIF
            END IF
300         CONTINUE  ! DO 300 NB = 1,KHKTBB
         END IF
200      CONTINUE  ! DO 200 NA = 1, KHKTA
100   CONTINUE  ! DO 100 IREP = 0, MAXREP
      RETURN
      END
C /* Deck sym1ev */
      SUBROUTINE SYM1EV(AO,DMAT,EXPVAL,FACSYM,MULTA,MULTB,IPRINT)
C
C     Expectation values, K.Ruud, Pisa June-05
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      DIMENSION AO(*), DMAT(*)
#include "onecom.h"
#include "symmet.h"
#include "ibtfun.h"
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 1)/2 + MIN(I,J)
C
      DO JREP = 0, MAXREP
         DO NA = 1, KHKTA
            IF (IBTAND(MULA,IBTXOR(JREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
               NAT = KHKTB*(NA - 1)
               IF (LDIAG) THEN
                  KHKTBB = NA
               ELSE
                  KHKTBB = KHKTB
               END IF
               DO NB = 1, KHKTBB
                  IF (IBTAND(MULB,IBTXOR(JREP,ISYMAO(NHKTB,NB))).EQ.0) 
     &                 THEN
                     FAC = FACSYM*PT(IBTAND(0,
     &                     IBTXOR(JREP,ISYMAO(NHKTB,NB))))/(MULTA*MULTB)
                     IORBAB = ITRI(IDENA+NA,IDENB+NB)
                     EXPVAL = EXPVAL + FAC*AO(NAT+NB)*DMAT(IORBAB)
                  END IF
               END DO
            END IF
         END DO
      END DO
      RETURN
      END
C  /* Deck sym1n */
      SUBROUTINE SYM1N(AO,SO,IREPO,KB,MULA,MULB,NHKTA,NHKTB,
     *                 KHKTA,KHKTB,HKAB,LDIAG,FULMAT,ANTI,THRESH,
     *                 IMAT0,IPRINT)
C
C     Take block of distinct AO integral (derivatives) and
C     generate symmetrized contributions to SO integral
C     (derivatives) over non-symmetric operators
C                                          880408  PRT
C     Modified tuh 880819
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      LOGICAL LDIAG, FULMAT, ANTI
      DIMENSION AO(*), SO(*)
#ifdef PRG_DIRAC
      PARAMETER (LUITMP = 48)
#else
#include "inftap.h"
#endif
#include "csym1.h"
#include "symmet.h"
#include "symind.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine SYM1N',-1)
         WRITE (LUPRI, '(A,2I5)') ' IREPO, KB ', IREPO, KB
         WRITE (LUPRI, '(A,2I5)') ' NHKTA/B ', NHKTA, NHKTB
      END IF
C
C     Loop over irreps for first basis function - those for second
C     are obtained from operator symmetry IREPO
C
      DO 100 IREPA = 0, MAXREP
         IREPB = IBTXOR(IREPO,IREPA)
         IF (FULMAT) INDOFF = NPARNU(IREPO+1,MAX(IREPA,IREPB)+1)
         IF (ANTI .AND. (IREPA .LT. IREPB)) THEN
            FAB = - HKAB
         ELSE
            FAB = HKAB
         END IF
C
C        Loop over AOs which are of symmetry IREPA in stabilizer MULA
C
         DO 200 NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))).EQ.0) THEN
            NAT = KHKTB*(NA - 1)
            IF (LDIAG) THEN
               KHKTBB = NA
            ELSE
               KHKTBB = KHKTB
            ENDIF
            DO 300 NB = 1,KHKTBB
            IF (NA.EQ.NB .AND. LDIAG .AND. IREPA .LT. IREPB) GOTO 300
            IF (IBTAND(MULB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))).EQ.0) THEN
C
C              Weight and parity factor
C
               FAC = FAB*PT(IBTAND(KB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
C
C              Locate SO integrals to which AOs contribute
C
               INDA = INDFA(IREPA + 1,NA)
               INDB = INDFB(IREPB + 1,NB)
               RINT = FAC*AO(NAT+NB)
               IF (FULMAT) THEN
                  IF (IREPA .GE. IREPB) THEN
                     IND  = INDOFF + NAOS(IREPB+1)*(INDA-1) + INDB
                  ELSE
                     IND  = INDOFF + NAOS(IREPA+1)*(INDB-1) + INDA
                  ENDIF
                  SO(IND) = SO(IND) + RINT
               ELSE
                  IF (ABS(RINT) .GT. THRESH) THEN
                     IF (IREPA .GE. IREPB) THEN
                        IND  = NAOS(IREPB+1)*(INDA-1) + INDB
                     ELSE
                        IND  = NAOS(IREPA+1)*(INDB-1) + INDA
                     ENDIF
                     INDMAX = MAX(IND,INDMAX)
                     LABEL  = IND*2**10 + IMAT0 + MAX(IREPA,IREPB)
                     IF (IPRINT .GT. 20) THEN
                       WRITE (LUPRI,'(A,F12.6,2I3,2I2,I4,I2,I5)')
     *                     'SYM1N - NA/B,IREPA/B,IND,IREPO,IMAT',
     *                     RINT, NA, NB, IREPA, IREPB, IND, IREPO,
     *                     IMAT0 + MAX(IREPA,IREPB)
                     END IF
                     IF (LENGTH .EQ. 600) THEN
                        WRITE (LUITMP) BUF, IBUF, LENGTH
                        IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *                     ' Buffer of length',LENGTH,
     *                     ' has been written in SYM1N.'
                        LENGTH = 0
                     ENDIF
                     LENGTH = LENGTH + 1
                     BUF(LENGTH)  = RINT
                     IBUF(LENGTH) = LABEL
                  ENDIF
               ENDIF
            END IF
300         CONTINUE
         END IF
200      CONTINUE
100   CONTINUE
      RETURN
      END
C  /* Deck sym1nc */
      SUBROUTINE SYM1NC(AO,DMAT,EXPVAL,IREPO,KB,MULA,MULB,NHKTA,
     *                  NHKTB,KHKTA,KHKTB,HKAB,LDIAG,ANTI,IPRINT)
C
C     Transform AO integrals to SO basis and contract with non-totally
C     symmetric density matrix provided in SO basis. Based on SYM1N.
C     K.Ruud, August 2005
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      LOGICAL LDIAG, FULMAT, ANTI
      DIMENSION AO(*), DMAT(*)
#include "inftap.h"
#include "csym1.h"
#include "symmet.h"
#include "symind.h"
#include "ibtfun.h"
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Subroutine SYM1NC',-1)
         WRITE (LUPRI, '(A,2I5)') ' IREPO, KB ', IREPO, KB
         WRITE (LUPRI, '(A,2I5)') ' NHKTA/B ', NHKTA, NHKTB
      END IF
C
C     Loop over irreps for first basis function - those for second
C     are obtained from operator symmetry IREPO
C
      DO 100 IREPA = 0, MAXREP
         IREPB = IBTXOR(IREPO,IREPA)
         JOFFA = 0
         JOFFB = 0
         IF (IREPA .GT. 0) JOFFA = JOFFA + ISUM(IREPA,NAOS,1)
         IF (IREPB .GT. 0) JOFFB = JOFFB + ISUM(IREPB,NAOS,1)
         INDOFF = NPARNU(IREPO+1,MAX(IREPA,IREPB)+1)
         IF (ANTI .AND. (IREPA .LT. IREPB)) THEN
            FAB = - HKAB
         ELSE
            FAB = HKAB
         END IF
C
C        Loop over AOs which are of symmetry IREPA in stabilizer MULA
C
         DO 200 NA = 1, KHKTA
         IF (IBTAND(MULA,IBTXOR(IREPA,ISYMAO(NHKTA,NA))).EQ.0) THEN
            NAT = KHKTB*(NA - 1)
            IF (LDIAG) THEN
               KHKTBB = NA
            ELSE
               KHKTBB = KHKTB
            ENDIF
            DO 300 NB = 1,KHKTBB
            IF (NA.EQ.NB .AND. LDIAG .AND. IREPA .LT. IREPB) GOTO 300
            IF (IBTAND(MULB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))).EQ.0) THEN
C
C              Weight and parity factor
C
               FAC = FAB*PT(IBTAND(KB,IBTXOR(IREPB,ISYMAO(NHKTB,NB))))
C
C              Locate SO integrals to which AOs contribute
C
               INDA = INDFA(IREPA + 1,NA)
               INDB = INDFB(IREPB + 1,NB)
               JNDA = JOFFA + INDA
               JNDB = JOFFB + INDB
               IF (JNDA .GT. JNDB) THEN
                  INDD = JNDA*(JNDA - 1)/2 + JNDB
               ELSE
                  INDD = JNDB*(JNDB - 1)/2 + JNDA
               END IF
               RINT = FAC*AO(NAT+NB)
               IF (IREPA .GE. IREPB) THEN
                  IND  = INDOFF + NAOS(IREPB+1)*(INDA-1) + INDB
               ELSE
                  IND  = INDOFF + NAOS(IREPA+1)*(INDB-1) + INDA
               ENDIF
               EXPVAL = EXPVAL - RINT*DMAT(INDD)
            END IF
300         CONTINUE
         END IF
200      CONTINUE
100   CONTINUE
      RETURN
      END
C  /* Deck sphrm1 */
      SUBROUTINE SPHRM1(CI,SPI,NTYPE,WORK,LWORK,PROPTY,IPRINT)
#include "implicit.h"
      LOGICAL PROPTY
      DIMENSION CI(KCKTAB,*), SPI(KCKTAB,*), WORK(LWORK)
#include "onecom.h"
      KTMP  = 1
      KTMQ  = KTMP + KCKTA*KCKTB
      KLAST = KTMQ + KCKTA*KHKTB
      IF (KLAST .GT. LWORK) CALL STOPIT('SPHRM1',' ',KLAST,LWORK)
      DO 100 I = 1, NTYPE
         CALL SPHRMX(CI(1,I),SPI(1,I),WORK(KTMP),WORK(KTMQ),PROPTY,
     &               IPRINT)
  100 CONTINUE
      RETURN
      END
C  /* Deck sphrmx */
      SUBROUTINE SPHRMX(CI,SPI,TMPINT,HALF,PROPTY,IPRINT)
C
C     Transform a block of Cartesian integrals to spherical harmonics
C
C                                          920511  PRT
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION CI(*), SPI(*), HALF(*), TMPINT(*)
      LOGICAL PROPTY
#include "onecom.h"
#include "sphtrm.h"
#include "symmet.h"
#include "ibtfun.h"
C
      IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI,'(A,2I6)') ' KCKT? ', KCKTA, KCKTB
         WRITE(LUPRI,'(A,2L6)') ' SPHR? ', SPHRA, SPHRB
         WRITE(LUPRI,'(A,4I6)') ' KHKT? ', KHKTA, KHKTB
         CALL HEADER('Cartesian integrals in SPHRMX',-1)
         IJ = 0
         DO 100 I = 1, KCKTA
            DO 110 J = 1,KCKTB
               IJ = IJ + 1
               WRITE(LUPRI,'(2I4,D13.6)') I,J, CI(IJ)
 110        CONTINUE
 100     CONTINUE
      END IF
C
      IF (SPHRA .OR. SPHRB) THEN
C
C        Collect integrals to be transformed
C        ===================================
C
         ICOFF = 0
         DO 200 ICOMPA = 1,KCKTA
            DO 220 ICOMPB = 1,KCKTB
               ICOFF = ICOFF + 1
               TMPINT((ICOMPA-1)*KCKTB + ICOMPB) = CI(ICOFF)
  220       CONTINUE
  200    CONTINUE
C
C        Transform second index (B) if required
C        ======================================
C
         IF (SPHRB) THEN
            CALL MXMA(CSP(ISPADR(NHKTB)),1,KHKTB,TMPINT,1,KCKTB,HALF,1,
     &                KHKTB,KHKTB,KCKTB,KCKTA)
ckr            CALL DGEMM('N','N',KHKTB,KCKTA,KCKTB,D1,CSP(ISPADR(NHKTB)),
ckr     &           KHKTB,TMPINT,KCKTB,D0,HALF,KHKTB)
         ELSE
            CALL DCOPY(KCKTA*KCKTB,TMPINT,1,HALF,1)
         END IF
C
C        Transpose half transformed integrals
C        ====================================
C
         DO 300 ICOMPA = 1,KCKTA
            DO 310 ICOMPB = 1,KHKTB
               TMPINT((ICOMPB-1)*KCKTA + ICOMPA)
     &               = HALF((ICOMPA-1)*KHKTB + ICOMPB)
  310       CONTINUE
  300    CONTINUE
C
C        Transform first index (A) if required
C        =====================================
C
         IF (SPHRA) THEN
            CALL MXMA(CSP(ISPADR(NHKTA)),1,KHKTA,TMPINT,1,KCKTA,HALF,1,
     &                KHKTA,KHKTA,KCKTA,KHKTB)
ckr            CALL DGEMM('N','N',KHKTA,KHKTB,KCKTA,D1,CSP(ISPADR(NHKTA)),
ckr     &           KHKTA,TMPINT,KCKTA,D0,HALF,KHKTA)
         ELSE
            CALL DCOPY(KCKTA*KHKTB,TMPINT,1,HALF,1)
         END IF
C
C        Collect transformed integrals
C        =============================
C
         ISOFF = 0
         DO 400 ICOMPA = 1,KHKTA
            DO 410 ICOMPB = 1,KHKTB
               ISOFF = ISOFF + 1
               SPI(ISOFF) = HALF((ICOMPB-1)*KHKTA + ICOMPA)
  410        CONTINUE
  400     CONTINUE
      END IF
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Spherical integrals in SPHRMX',-1)
         IJ = 0
         DO 500 I = 1, KHKTA
            DO 510 J = 1,KHKTB
               IJ = IJ + 1
               WRITE(LUPRI,'(A,(2I4,D13.6))') '++',I,J,SPI(IJ)
 510        CONTINUE
 500     CONTINUE
      END IF
      RETURN
      END
C  -- end of her1sym.F --
