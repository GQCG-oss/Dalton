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
#ifdef OLD_REV_LOG
!===========================================================================
! look for HJMAERKE; (1) symort of UMAT should be done
!
!Revision 1.2  2001/01/17 10:41:32  vebjornb
!Calls to *MPA*B* in arhpack.F have been replaced with DGEMM calls
!
!Revision 1.4  2000/05/18 13:46:43  hjj
!Bugfix in AVERSS, it is strange if this bug hasn't caused errors before!
! KLWOP was referenced out of bounds for inactive orbitals (a check was
! missing), this could lead to any value and in a test today it lead to
! bus error.
!
!Revision 1.3  2000/05/01 13:13:37  hjj
!(1) removed AVEOLD, AVECMO, AVECHO connected with old KAVER option
!    (now removed)
!(2) dynamically allocate KLWOP in AVERSS
!
!revision 1.2 dec 1999 hjj
!improved error messages for super symmetry errors
!
!961213-hjaaj:
!AVESET: disable SUPSYM if no supersymmetry found
!        HJMAERKE: this means that debug SUPSYM for abelian group
!        cannot be performed without commenting new code in AVESET!
!941215-hjaaj:
!renamed ICKDIA to NOFDIA; moved NOFDIA to gp
!940909-hjaaj: CALL SCALSUB2 after SUBROUTINE for CRAY (see 940427 comment)
! changed CALL SCALSUB to CALL SCALSUB1 before SUBROUTINE
!940823-hjaaj:
!AVESE1: print if IPHCHA .gt. 0 and IPRAVE .ge. 4
!940818-hjaaj:
!AVERSS,AVEPHA,AVEDEG: more print before certain STOPs to aid user
!940815-hjaaj:
!AVESE1,AVEDEG,AVEGTU,AVEPHA: use THRSSY from INFINP (new in INFINP)
!   instead of THREQL parameter; include INFINP
!AVESE1: Check and reset THRSSY if needed; more print if disentangling
!940427-hjaaj: SCALAR directives for Cray (as SCALSUB does not work
!  for Cray because directive has to be after SUBROUTINE statement)
!931028-hjaaj
!AVEPRT: rettet print
!930715-hjaaj
!AVESE0: define MXDGSS and NINFSS(*,3)
!AVEDEG: define MXDGSS
!AVESET: removed call to AVEPRT (moved to PRTINP)
!930625-hjaaj
!  inserted some more print before two stop statements.
!===========================================================================
#endif
C  /* Deck averag */
      SUBROUTINE AVERAG (VECS,NVDIM,NVSIM)
C
C 7-Jul-1992 hjaaj+hh
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION VECS(NVDIM,NVSIM)
C
C Used from common blocks:
C  INFINP : SUPSYM
C  INFVAR : NWOPT,JWOPSY
C
#include <maxorb.h>
#include <infinp.h>
#include <infvar.h>
#include <infpri.h>
C
      CALL QENTER('AVERAG')
      IF (.NOT. SUPSYM) GO TO 9999
      IF (NWOPT  .EQ. 0) GO TO 9999
      IF (JWOPSY .NE. 1) GO TO 910
      IF (NVDIM  .LT. NWOPT) GO TO 920
      IF (SUPSYM) THEN
         CALL AVERSS(VECS,NVDIM,NVSIM)
      END IF
 9999 CALL QEXIT('AVERAG')
      RETURN
C
C *** Error sections
C
  910 CONTINUE
      WRITE (LUERR,9010) JWOPSY
      CALL QTRACE(LUERR)
      CALL QUIT('AVERAG ERROR, operator symmetry JWOPSY .ne. 1')
 9010 FORMAT(/' ERROR-AVERAG, operator symmetry .ne. 1; JWOPSY =',I8)
C
  920 CONTINUE
      WRITE (LUERR,9020) NVDIM,NWOPT
      CALL QTRACE(LUERR)
      CALL QUIT('AVERAG ERROR, vector length lt NWOPT')
 9020 FORMAT(/' ERROR-AVERAG, vector length',I8,' is less than NWOPT',
     *   I8)
C
      END
C  /* Deck averss */
      SUBROUTINE AVERSS (VECS,NVDIM,NVSIM)
C
C 7-Jul-1992 Hans Joergen Aa. Jensen + Hinne Hettema
C
C Purpose: Average degenerate super symmetries
C          (which are different components of a degenerate irrep).
C          This corresponds to doing a microcanonical average.
C
#include <implicit.h>
      DIMENSION VECS(NVDIM,NVSIM)
C
C Used from common blocks:
C   INFORB : NSSYM,NORBSS(),IORBSS(),NINFSS(*,3), NOCCT,NORBT
C   INFIND : ISW(*),ISSORD(*)
C   INFVAR : NWOPT
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infind.h>
#include <inforb.h>
#include <infvar.h>
#include <infpri.h>
C
      PARAMETER (MAXDEG = 20, D1 = 1.0D0)
      DIMENSION ISSDEG(MAXDEG),ISSOFF(MAXDEG),IDEGV(MAXDEG)
C
Chjaaj-may2000:
C     Dynamic allocation of KLWOP with adjustable run-time dimensions
      DIMENSION KLWOP(NOCCT,NORBT)
C
      CALL QENTER('AVERSS')
      CALL MAKE_KLWOP(KLWOP)
      IF (IPRAVE .GT. 100) THEN
         WRITE(LUPRI,'(/A)') ' AVERSS, Incoming orbital vector(s):'
         CALL OUTPUT(VECS,1,NWOPT,1,NVSIM,NVDIM,NVSIM,-1,LUPRI)
      END IF
C
      DO 800 ISSYMR = 1,NSSYM
C        Skip this supsym if not degenerate or not "root supsym"
         IF (NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF (NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
         NDEG = 0
         DO 100 ISSYM = ISSYMR,NSSYM
            IF (NINFSS(ISSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISSDEG(NDEG) = ISSYM
               ISSOFF(NDEG) = IORBSS(ISSYM) - IORBSS(ISSYMR)
            END IF
  100    CONTINUE
         IF (IPRAVE .GT. 5 .OR. NDEG .NE. NINFSS(ISSYMR,1)) THEN
            WRITE(LUPRI,'(/A,I3,A/A,(T35,8I5))')
     &         ' AVERSS averaging ',NDEG,' degenerate super symmetries',
     &         ' namely sup.sym.s no.',(ISSDEG(IDEG),IDEG=1,NDEG)
         END IF
         IF (NDEG .NE. NINFSS(ISSYMR,1)) THEN
            WRITE(LUPRI,'(/A)')
     &      ' AVERSS-software error: inconsistent degeneracy in NINFSS'
            CALL AVEDMP(LUPRI)
            CALL QUIT('AVERSS: inconsistent degeneracy in NINFSS')
         END IF
C
         AVFAC = NDEG
         AVFAC = D1 / AVFAC
C
         IOFFSR = IORBSS(ISSYMR)
         DO 480 IMOSS = IOFFSR+1,IOFFSR+NORBSS(ISSYMR)
            IMO = ISSORD(IMOSS)
            IMOW = ISW(IMO)
            IF (IMOW .GT. NOCCT) GO TO 480
            DO 470 JMOSS = IMOSS+1,IOFFSR+NORBSS(ISSYMR)
               JMO = ISSORD(JMOSS)
               JMOW = ISW(JMO) - NISHT
               IF (JMOW .LE. 0) GO TO 470
Chjaaj         ... this is an inactive orbital
               IF (KLWOP(IMOW,JMOW) .EQ. 0) GO TO 470
               IDEGV(1) = KLWOP(IMOW,JMOW)
               NERR = 0
               DO 410 IDEG = 2,NDEG
                  KMO = ISSORD(IMOSS + ISSOFF(IDEG))
                  LMO = ISSORD(JMOSS + ISSOFF(IDEG))
                  LMOW = ISW(LMO) - NISHT
                  IF (LMOW .LE. 0) THEN
                     IDEGV(IDEG) = 0
                     NERR = NERR + 1
                     WRITE(LUPRI,'(/A)')
     &       'AVERSS error: degenerate orbitals not same orbital class.'
                  ELSE
                     IDEGV(IDEG) = KLWOP( ISW(KMO) , LMOW )
                     IF (IDEGV(IDEG) .EQ. 0) NERR = NERR + 1
                  END IF
  410          CONTINUE
               IF (IPRAVE .GT. 10 .OR. NERR .GT. 0) THEN
                  IF (NERR .GT. 0) THEN
                     WRITE(LUPRI,'(/A/A/)')
     &         ' AVERSS error: not all degenerate rotations included',
     &         ' Probable error: some but not all'//
     &            ' of the degenerate orbitals have been frozen.'
                  END IF
                  WRITE(LUPRI,'(A,(T35,4I10))')
     &               ' AVERSS averaging orb.rot. no.s',
     &               (IDEGV(IDEG),IDEG=1,NDEG)
                  IF (NERR .GT. 0) CALL AVEDMP(LUPRI)
               END IF
               IF (NERR .GT. 0) THEN
                  CALL QTRACE(LUPRI)
                  CALL QUIT(
     &            'AVERSS error: not all deg. rotations included')
               END IF
               DO 460 IVSIM = 1,NVSIM
                  XKAPIJ = VECS(IDEGV(1),IVSIM)
                  DO 420 IDEG = 2,NDEG
                     XKAPIJ = XKAPIJ + VECS(IDEGV(IDEG),IVSIM)
  420             CONTINUE
                  XKAPIJ = AVFAC*XKAPIJ
                  DO 430 IDEG = 1,NDEG
                     VECS(IDEGV(IDEG),IVSIM) = XKAPIJ
  430             CONTINUE
  460          CONTINUE
  470       CONTINUE
  480    CONTINUE
  800 CONTINUE
C
      IF (IPRAVE .GT. 100) THEN
         WRITE(LUPRI,'(/A)') ' AVERSS, Averaged orbital vector(s):'
         CALL OUTPUT(VECS,1,NWOPT,1,NVSIM,NVDIM,NVSIM,-1,LUPRI)
      END IF
C
      CALL QEXIT('AVERSS')
      RETURN
      END
C  /* Deck aveset */
      SUBROUTINE AVESET(CMO,WRK,KFREE,LFREE)
C
C 2-Jul-1992 hjaaj : core allocation for AVESE1
C
#include <implicit.h>
      DIMENSION CMO(*), WRK(*)
C
C Used from common blocks:
C  INFINP: SUPSYM
C  INFORB: N2ORBX, MXSSYM
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
      CALL QENTER('AVESET')
      IF (SUPSYM) THEN
         IF (IPRAVE .GE. 3) THEN
            WRITE(LUPRI,'(//A)') ' TEST OUTPUT FROM AVESET'
         END IF
         KFRSAV = KFREE
         LFRSAV = LFREE
         CALL MEMGET('REAL',KTKMO,N2ORBX,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KNDEGS,MXSSYM,WRK,KFREE,LFREE)
         CALL AVESE1(CMO,WRK(KTKMO),WRK(KNDEGS),
     *               WRK,KFREE,LFREE)
         CALL MEMREL('AVESET',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C        Check if supersymmetry is identical to "D2h" symmetry
C
         IF (NSSYM .LE. NSYM .AND. MXDGSS .EQ. 1) THEN
            DO 200 ISYM = 1,NSYM
               JSYM = 0
               DO 100 ISSYM = 1,NSSYM
                  IF (NINFSS(ISSYM,3) .EQ. ISYM) JSYM = JSYM + 1
  100          CONTINUE
               IF (JSYM .GT. 1) GO TO 300
  200       CONTINUE
            IF (IPRAVE .GE. 1) THEN
               WRITE(LUPRI,'(/A)')' AVESET: No "super symmetry" found.'
            END IF
            SUPSYM = .FALSE.
            GO TO 300
         END IF
  300    CONTINUE
      END IF
      IF (.NOT.SUPSYM) THEN
         CALL AVESE0
      END IF
      CALL QEXIT('AVESET')
      RETURN
      END
C  /* Deck avese0 */
      SUBROUTINE AVESE0
C
C 7-Jul-1992 hjaaj
C Set super symmetry information as symmetry information
C (when super symmetry is not to be used()
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB : NSSYM, NORBSS(), IORBSS(); NSYM, NORB(), IORB(), NORBT
C  INFIND : ISSMO(), ISSORD(); ISMO()
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('AVESE0')
      NSSYM = NSYM
      DO 100 ISYM = 1,NSYM
         NORBSS(ISYM) = NORB(ISYM)
         IORBSS(ISYM) = IORB(ISYM)
         NINFSS(ISYM,1) = 1
         NINFSS(ISYM,2) = ISYM
         NINFSS(ISYM,3) = ISYM
  100 CONTINUE
      DO 200 IMO = 1,NORBT
         ISSMO(IMO)  = ISMO(IMO)
         ISSORD(IMO) = IMO
  200 CONTINUE
      MXDGSS = 1
      CALL QEXIT('AVESE0')
      RETURN
      END
C  /* Deck avese1 */
      SUBROUTINE AVESE1(CMO,TKMO,NDEGSS,WRK,KFREE,LFREE)
C
C Written 29-Jun-1992 Hans Joergen Aa. Jensen and Hinne Hettema
C
C Purpose:
C  Find super symmetry and degeneracy of molecular orbitals in CMO
C  Disentangle mixed degenerate orbitals
C  Coordinate relative phases of degenerate orbitals
C
C Input:
C  CMO; molecular orbitals
C
C Output:
C
#include <implicit.h>
      DIMENSION CMO(NCMOT), TKMO(NORBT,NORBT), WRK(*)
      DIMENSION NDEGSS(*)
C     DIMENSION NDEGSS(MXSSYM)
C
#include <threql.h>
      PARAMETER (TMXSSY = 1.0D-4)
C
C Used from common blocks:
C   INFINP : THRSSY
C   INFORB : NSYM,NSSYM,NORBT,NORB(),?
C   INFIND : ISSMO(),ISSORD(),?
C   INFPRI : IPRAVE,NINFO,?
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
C
      CALL QENTER('AVESE1')
      IF (THRSSY .GT. TMXSSY) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A/2(A,1P,D10.2))') ' AVESE1-WARNING: SUPSYM'//
     &      ' threshold for degeneracy and "zero" elements too big;',
     &      '                 value reset from',THRSSY,' to',TMXSSY
         THRSSY = TMXSSY
      END IF
      IF (THRSSY .LT. THREQL) THEN
C        ... THREQL must be .gt. accuracy of diagonalization routine
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/A/2(A,1P,D10.2))') ' AVESE1-INFO:'//
     &      ' SUPSYM threshold for degeneracy and "zero" elements',
     &      '              is reset from',THRSSY,' to',THREQL
         THRSSY = THREQL
      END IF
      KFRSAV = KFREE
      ITURN = 0
  100 ITURN = ITURN + 1
      KFREE1 = KFREE
C
C Step 0: Obtain TKMO matrix for finding super symmetries
C         TKMO must be the matrix of a totally symmetric
C         operator which is very unlikely to have
C         accidental zeroes (we use the kinetic energy).
C
      CALL AVETK(TKMO,CMO,WRK,KFREE,LFREE)
C     Note: this can give problems for FINITE FIELD (NFIELD.gt.0)
C     because the TKMO will reflect the symmetry without the finite field(s)!
C
C Step 1: Find super symmetry of each MO
C         Define ISSMO(i)  = super sym. of m.o. no. i
C                ISSORD(i) = sup. sym. reordering array
C                NSSYM     = no. of super symmetries
C                NORBSS(*) = no. of mo's in each sup.sym.
C                IORBSS(*) = pointer array to ISSORD
C
      CALL IZERO(ISSMO,NORBT)
      CALL IZERO(ISSORD,NORBT)
      NSSYMX = 0
      DO 180 ISYM = 1,NSYM
         IF (NORB(ISYM) .EQ. 0) GO TO 180
         IOFFMO = IORB(ISYM)
         DO 170 IMO = IOFFMO+1,IOFFMO+NORB(ISYM)
            IF (ISSMO(IMO) .EQ. 0) THEN
               NSSYMX = NSSYMX + 1
               IF (NSSYMX .GT. MXSSYM) GO TO 901
C     v-------- error exit
               ISSMO(IMO) = NSSYMX
               NORBSS(NSSYMX) = 1
            END IF
            ISSYM = ISSMO(IMO)
            DO 160 JMO = IMO+1,IOFFMO+NORB(ISYM)
            IF (ABS(TKMO(JMO,IMO)) .GT. THRSSY) THEN
            IF (ISSMO(JMO) .EQ. 0) THEN
C -- We have a new orbital in the supersymmetry rep.
               ISSMO(JMO) = ISSYM
               NORBSS(ISSYM) = NORBSS(ISSYM) + 1
            ELSE IF (ISSMO(JMO) .NE. ISSYM) THEN
C -- We have a conflict and must collapse the two supersym.s
               JSSYMH = MAX(ISSMO(JMO),ISSYM)
               JSSYML = MIN(ISSMO(JMO),ISSYM)
               IF (IPRAVE .GE. 4) THEN
                  WRITE(LUPRI,'(/A,2I4)')
     *             ' Collapsing super symmetries :',JSSYML,JSSYMH
               END IF
               DO 140 KMO = IOFFMO+1,IOFFMO+NORB(ISYM)
                  IF (ISSMO(KMO) .EQ. JSSYMH) ISSMO(KMO) = JSSYML
  140          CONTINUE
               ISSYM  = JSSYML
               NORBSS(JSSYML) = NORBSS(JSSYML) + NORBSS(JSSYMH)
               NORBSS(JSSYMH) = 0
            END IF
            END IF
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
C
C  in case of debugging, get out some print
C
      IF ( IPRAVE .GE. 10) THEN
         NSSYM = NSSYMX
         CALL IZERO(IORBSS,NSSYM)
         WRITE(LUPRI,'(/A/A/)')
     *      ' After initializing in AVESE1 (step 1)',
     *      ' ====================================='
         CALL AVEDMP(LUPRI)
      END IF
C
C Step 2: Find NSSYM, IORBSS(*) and ISSORD(*)
C
      NSSYM = 0
      IMOSS = 0
      IORBSS(1) = 0
      DO 260 ISSYMX = 1,NSSYMX
      IF (NORBSS(ISSYMX) .GT. 0) THEN
         NSSYM = NSSYM + 1
         IF (NSSYM .NE. ISSYMX) THEN
            NORBSS(NSSYM) = NORBSS(ISSYMX)
            DO 220 IMO = 1,NORBT
               IF (ISSMO(IMO) .EQ. ISSYMX) ISSMO(IMO) = NSSYM
  220       CONTINUE
         END IF
         IORBSS(NSSYM) = IMOSS
         DO 240 IMO = 1,NORBT
            IF (ISSMO(IMO) .EQ. NSSYM) THEN
               IMOSS = IMOSS + 1
               ISSORD(IMOSS) = IMO
            END IF
  240    CONTINUE
         IF ( (IMOSS - IORBSS(NSSYM)) .NE. NORBSS(NSSYM) ) THEN
            GO TO 902
         END IF
      END IF
  260 CONTINUE
C
C  in case of debugging, get out some print
C
      IF ( IPRAVE .GE. 3) THEN
         WRITE(LUPRI,'(/A/A/)')
     *      ' After clearing up in AVESE1 (step 2)',
     *      ' ===================================='
         CALL AVEDMP(LUPRI)
      END IF
C
C Step 3: Find degeneracies
C
      CALL MEMGET('INTE',KINDX1,NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDX2,NORBT,WRK,KFREE,LFREE)
      CALL AVEDEG(TKMO,WRK(KINDX1),WRK(KINDX2),NDEGSS)
C
C Step 4: Disentangle degenerate super symmetries
C
      NDIS = 0
      DO 480 ISSYM = 1,NSSYM
         NDEGI = NDEGSS(ISSYM)
         IF (NDEGI .GT. 1) THEN
            NDIS = NDIS + 1
         IF (ITURN .EQ. 1) THEN
            ISYM  = NINFSS(ISSYM,3)
            NORBI = NORB(ISYM)
            NBASI = NBAS(ISYM)
            JCMOI = ICMO(ISYM) + 1
            NSSMOI = NORBSS(ISSYM)
            NSSMON = NSSMOI / NDEGI
            IF (NDEGI*NSSMON .NE. NSSMOI) THEN
               CALL QUIT('AVESE1 error: NDEGI*NSSMON .ne. NSSMOI')
            END IF
            KFREE4 = KFREE
            CALL MEMGET('REAL',KUMAT,NDEGI*NDEGI,WRK,KFREE,LFREE)
            CALL MEMGET('INTE',KIMODE,NSSMOI,WRK,KFREE,LFREE)
            CALL MEMGET('REAL',KWRK,NBASI*NDEGI,WRK,KFREE,LFREE)
            CALL AVEDIS(ISSYM,ISYM,NDEGI,NORBI,NBASI,NSSMOI,NSSMON,
     *                  CMO(JCMOI),TKMO,WRK(KINDX1),
     *                  WRK(KUMAT),WRK(KIMODE),WRK(KWRK))
            CALL MEMREL('AVESE1.AVEDIS',WRK,1,KFREE4,KFREE,LFREE)
         END IF
         END IF
  480 CONTINUE
C
      CALL MEMREL('AVESE1.step 4',WRK,1,KFREE1,KFREE,LFREE)
C
      IF (NDIS .GT. 0) THEN
         IF (ITURN .EQ. 1) THEN
            IF (IPRAVE .GE. 3) WRITE(LUPRI,'(/A)')
     &         ' AVESE1: new super symmetry analysis after'//
     &         ' disentangling degenerate super symmetries.'
            GO TO 100
         END IF
         WRITE(LUPRI,'(/A)')
     &      ' AVESE1 error: disentangling failed, sorry!'
         CALL QUIT('AVESE1 error: disentangling failed, sorry!')
      END IF
C
C Step 5: Check (and change) relative phases of MO's
C
      IF (IPRAVE .GE. 4) THEN
         WRITE(LUPRI,'(/A)')
     &      ' AVESE1 : checking relative phases of degenerate MOs'
      END IF
      IPHCHA = 0
      DO 580 ISSYM = 1, NSSYM
         IDEG   = NINFSS(ISSYM,1)
         ISSYMR = NINFSS(ISSYM,2)
         IF (IDEG .GT. 1 .AND. ISSYM .EQ. ISSYMR) THEN
C        --- we have a degeneracy to other supersym.(s)
            CALL AVEPHA(ISSYM,CMO,TKMO,IPHCHA)
         END IF
  580 CONTINUE
      IF (IPRAVE .GE. 4 .AND. IPHCHA .GT. 0) THEN
         WRITE(LUPRI,'(A,I5)')
     &      ' AVESE1 : # of phase shifts of degenerate MOs :',IPHCHA
      END IF
C
C *** end of subroutine AVESE1
C
      CALL QEXIT('AVESE1')
      RETURN
C
  901 CONTINUE
      CALL QUIT('AVESE1: NSSYM .gt. MXSSYM in AVESET code 901')
  902 CONTINUE
      CALL QUIT('AVESE1: Inconsistent NORBSS(NSSYM) code 902')
      END
C  /* Deck aveord */
      SUBROUTINE AVEORD()
C
C Aug. 2004 hjaaj - create ISSORD() from ISSMO()
C
C The order of the super symmetry orbitals may be changed
C in ORDRSS and ORD2SS calls, thus we need to remake ISSORD()
C to be sure it is correct.
C
#include <implicit.h>
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
#include <infind.h>
#include <inforb.h>
#include <infpri.h>
C
      CALL QENTER('AVEORD')
      IMOSS = 0
      IORBSS(1) = 0
      DO ISSYM = 1,NSSYM
         IF (IORBSS(ISSYM) .NE. IMOSS) THEN
            CALL AVEDMP(LUPRI)
            CALL QUIT('AVEORD: Inconsistent IORBSS(ISSYM)')
         END IF
         DO IMO = 1,NORBT
            IF (ISSMO(IMO) .EQ. ISSYM) THEN
               IMOSS = IMOSS + 1
               ISSORD(IMOSS) = IMO
            END IF
         END DO
         IF ( (IMOSS - IORBSS(ISSYM)) .NE. NORBSS(ISSYM) ) THEN
            CALL AVEDMP(LUPRI)
            CALL QUIT('AVEORD: Inconsistent NORBSS(ISSYM)')
         END IF
      END DO
C
C  In case of debugging, get out some print
C
      IF ( IPRAVE .GE. 3) THEN
         WRITE(LUPRI,'(/A/A/)')
     *      ' After remake of ISSORD() in AVEORD  ',
     *      ' ===================================='
         CALL AVEDMP(LUPRI)
      END IF
      CALL QEXIT('AVEORD')
      RETURN
      END
C  /* Deck avetk */
      SUBROUTINE AVETK(TKMO,CMO,WRK,KFREE,LFREE)
C
C 920629-hjaaj:
C Purpose: obtain TKMO(NORBT,NORBT)
C
#include <implicit.h>
      DIMENSION CMO(*), TKMO(NORBT,NORBT), WRK(*)
C
      CHARACTER*8 TKLABEL
      PARAMETER (TKLABEL = 'KINETINT')
C
C Used from common blocks:
C   INFORB : NORBT,NNBAST,...
C   INFPRI : IPRAVE
C
#include <priunit.h>
#include <inforb.h>
#include <infpri.h>
C
      LOGICAL FOUND
C
      CALL QENTER('AVETK ')
      KFRSAV = KFREE
      CALL MEMGET('REAL',KTKAO,NNBAST,WRK,KFREE,LFREE)
      FOUND = .FALSE.
      CALL RDONEL(TKLABEL,FOUND,WRK(KTKAO),NNBAST)
      IF (.NOT. FOUND) THEN
         CALL QUIT('AVETK error: TK integrals not found')
      END IF
C
C
      CALL MEMGET('REAL',KTKPK,NNORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCR1,NNORBX,WRK,KFREE,LFREE)
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         NBASI = NBAS(ISYM)
         JTKAO = KTKAO + IIBAS(ISYM)
         JTKPK = KTKPK + IIORB(ISYM)
         JCMO  = 1     + ICMO(ISYM)
C
         CALL UTHU(WRK(JTKAO),WRK(JTKPK),CMO(JCMO),WRK(KSCR1),
     *             NBASI,NORBI)
C        CALL UTHU(H,HT,U,WRK,NAO,NMO)
  200 CONTINUE
C
      CALL PKSYM1(WRK(KSCR1),WRK(KTKPK),NORB,NSYM,-1)
C     CALL PKSYM1(ASP,APK,NDIM,NBLK,IWAY)
C
      CALL DSPTSI(NORBT,WRK(KSCR1),TKMO)
C     CALL DSPTSI(N,ASP,ASI)
C
      CALL MEMREL('AVETK',WRK,1,KFRSAV,KFREE,LFREE)
C
      IF ( IPRAVE .GE. 15 ) THEN
         WRITE(LUPRI,'(3(/A))')
     &      ' Matrix used for determining supersymmetries',
     &      ' (the kinetic energy matrix)',
     &      ' ==========================================='
         CALL OUTPUT(TKMO,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL QEXIT('AVETK ')
      RETURN
      END
C  /* Deck avedeg */
      SUBROUTINE AVEDEG(TKMO,INDX1,INDX2,NDEGSS)
C
C   920630 hjaaj
C   Purpose : find degeneracies, and construct
C             INDX1(*) = labeling degenerate orbitals
C             INDX2(*) = degeneracies of labels in INDX1
C
#include <implicit.h>
      DIMENSION TKMO(NORBT,NORBT),INDX1(NORBT),INDX2(NORBT)
      DIMENSION NDEGSS(*)
C     DIMENSION NDEGSS(MXSSYM)
C
C Used from common blocks:
C  INFINP : THRSSY
C  INFORB : NORBT, NSSYM, NINFSS(MXSSYM,3),MXDGSS
C  INFIND : ISSMO()
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      CALL QENTER('AVEDEG')
      CALL IZERO(INDX1,NORBT)
      CALL IZERO(INDX2,NORBT)
      IMOX = 0
      NDEGMX = 0
      DO 300 IMO = 1,NORBT
      IF (INDX1(IMO) .EQ. 0) THEN
         IMOX = IMOX + 1
         INDX1(IMO) = IMOX
         INDX2(IMOX) = 1
         TKMOI = TKMO(IMO,IMO)
         DO 200 JMO = IMO+1,NORBT
            IF (ABS(TKMO(JMO,JMO)-TKMOI) .LT. THRSSY) THEN
               INDX1(JMO)  = IMOX
               INDX2(IMOX) = INDX2(IMOX) + 1
            END IF
  200    CONTINUE
         NDEGMX = MAX(NDEGMX,INDX2(IMOX))
      END IF
  300 CONTINUE
      MXDGSS = NDEGMX
      NORXT = IMOX
C
C     Define NINFSS(i,1) = degeneracy of MO's in supsym_i
C            NINFSS(i,2) = "root" supsym of this degeneracy
C            NINFSS(i,3) = Abelian symmetry of supsym_i
C            NDEGSS(i  ) = no. of deg. MO's in this supsym
C
      CALL IZERO(NINFSS,MXSSYM*3)
      CALL IZERO(NDEGSS,MXSSYM)
      IDGERR = 0
      DO 500 IMO = 1,NORBT
         ISSYM = ISSMO(IMO)
         IMOX  = INDX1(IMO)
         NDEGI = INDX2(IMOX)
         IF (NINFSS(ISSYM,1) .EQ. 0) THEN
            NINFSS(ISSYM,1) = NDEGI
            NINFSS(ISSYM,2) = ISSYM
            NINFSS(ISSYM,3) = ISMO(IMO)
            NDEGSS(ISSYM)   = 1
            DO 400 JMO = IMO+1,NORBT
            IF (INDX1(JMO) .EQ. IMOX) THEN
               JSSYM = ISSMO(JMO)
               NINFSS(JSSYM,1) = NDEGI
               NINFSS(JSSYM,2) = ISSYM
               NINFSS(JSSYM,3) = ISMO(JMO)
               NDEGSS(JSSYM)   = NDEGSS(JSSYM) + 1
            END IF
  400       CONTINUE
         ELSE IF (NINFSS(ISSYM,1) .NE. NDEGI) THEN
            WRITE(LUPRI,'(A,4I4)')
     &           ' Error for IMO,ISSYM, NINFSS(ISSYM,1),NDEGI  = ',
     &           IMO,ISSYM, NINFSS(ISSYM,1),NDEGI
            IDGERR = IDGERR + 1
         END IF
  500 CONTINUE
C
C Print NINFSS for debugging
C
      IF ( IPRAVE .GE. 3 .OR. IDGERR .GT. 0) THEN
         IF (IDGERR .GT. 0) THEN
            WRITE(LUPRI,'(7(/A))')
     &      ' AVEDEG: Degeneracy error, dump of NINFSS follows below.',
     &      ' Probable cause:'//
     &         ' nuclear geometry has only approximate super symmetry.',
     &      ' Some possible actions:',
     &' (1) include more digits in nuclear coordinates',
     &' (2) attempt to force super symmetry by making .THRSSY larger',
     &' (3) attempt to avoid degeneracy error by making .THRSSY smaller'
     &,' (4) disable supersymmetry by removing .SUPSYM'
         END IF
         WRITE(LUPRI,'(4(/A))')
     &   ' NINFSS(i,1) = degeneracy of MOs in supsym_i',
     &   ' NINFSS(i,2) = root supsym of this degeneracy',
     &   ' NINFSS(i,3) = Abelian symmetry of supsym_i',
     &   ' NDEGSS(i)   = no. of deg. supsyms currently in this supsym'
         WRITE(LUPRI,'(/A/A)')
     &      '  ISS         NINFSS 1 to 3            NDEGSS',
     &      ' ============================================'
         DO 710 ISSYM = 1, NSSYM
            WRITE(LUPRI,'(1X,I4,4I10)')
     &         ISSYM,(NINFSS(ISSYM,J),J=1,3),NDEGSS(ISSYM)
  710    CONTINUE
         WRITE(LUPRI,'(A/)')
     &      ' ============================================'
      END IF
C
C     print for debuggers
C
      IF ( IPRAVE .GE. 10 .OR. IDGERR .GT. 0 ) THEN
         WRITE(LUPRI,'(A,F15.12/)')
     &      ' Degeneracy pattern found with THRSSY =',THRSSY
         WRITE(LUPRI,'("    I  ISSMO INDX1 INDX2       TKMO(I,I)")')
       WRITE(LUPRI,'(" ============================================")')
         DO 700 IMO = 1, NORBT
            WRITE(LUPRI,'(I5,3I6,F25.15)')
     &         IMO, ISSMO(IMO), INDX1(IMO), INDX2(IMO), TKMO(IMO,IMO)
  700    CONTINUE
       WRITE(LUPRI,'(" ============================================")')
      END IF
      IF (IDGERR .GT. 0) THEN
         CALL QUIT('AVEDEG: Degeneracy error ...')
      END IF
C
      CALL QEXIT('AVEDEG')
      RETURN
      END
C  /* Deck avedis */
      SUBROUTINE AVEDIS(ISSYM,ISYM,NDEGI,NORBI,NBASI,NSSMOI,NSSMON,
     *                  CMOI,TKMO,INDX1,UMAT,IMODEG,WRK)
C
#include <implicit.h>
C
      DIMENSION CMOI(NBASI,NORBI), TKMO(NORBT,NORBT)
      DIMENSION INDX1(NORBT),UMAT(NDEGI,NDEGI)
      DIMENSION IMODEG(NDEGI,NSSMON), WRK(NBASI,NDEGI)
C
C Used from common blocks
C  INFORB: NSSYM, NORBT, NORBSS(), IORBSS()
C  INFIND: ISMO(),ISSORD()
C  INFPRI: IPRAVE
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      CALL QENTER('AVEDIS')
      IF (IPRAVE .GE. 3) THEN
         WRITE(LUPRI,'(/A,I5/A,I5)')
     *      ' AVEDIS disentangling mixed degenerate components in'//
     *      ' supersymmetry:',ISSYM,' Abelian symmetry:',ISYM
      END IF
      IMOXP = 0
      KMOSS = 0
      IMOSSE = IORBSS(ISSYM)+NORBSS(ISSYM)
      DO 180 IMOSS = IORBSS(ISSYM)+1,IMOSSE
         IMO = ISSORD(IMOSS)
         IMOX = INDX1(IMO)
         IF (IMOX .GT. IMOXP) THEN
            IMOXP = IMOX
            KMOSS = KMOSS + 1
            IDEG = 1
            IMODEG(1,KMOSS) = IMO
            DO 120 JMOSS = IMOSS+1,IMOSSE
               JMO = ISSORD(JMOSS)
               IF (INDX1(JMO) .EQ. IMOX) THEN
                  IDEG = IDEG + 1
                  IMODEG(IDEG,KMOSS) = JMO
               END IF
  120       CONTINUE
         END IF
  180 CONTINUE
      IF (KMOSS .NE. NSSMON) THEN
         CALL QUIT('AVEDIS error: KMOSS .ne. NSSMON')
      END IF
C
      IOFFMO = IORB(ISYM)
C     check that first block is diagonal
      CALL AVEGTU(NDEGI,NORBT,1,1,
     *            IMODEG(1,1),IMODEG(1,1),UMAT,TKMO,WRK)
      DO 280 KMOSS = 2,NSSMON
C        a) check that KMOSS block is diagonal
         CALL AVEGTU(NDEGI,NORBT,KMOSS,KMOSS,
     *               IMODEG(1,KMOSS),IMODEG(1,KMOSS),UMAT,TKMO,WRK)
C        b) find UMAT to transform the KMOSS degenerate orbital set to
C           the same irrep rows as the first deg. orbital set
         CALL AVEGTU(NDEGI,NORBT,1,KMOSS,
     *               IMODEG(1,1),IMODEG(1,KMOSS),UMAT,TKMO,WRK)
         CALL DZERO(WRK,NDEGI*NBASI)
C        c) do the transformation
         DO 240 IDEG = 1,NDEGI
            IMO = IMODEG(IDEG,KMOSS) - IOFFMO
            DO 220 JDEG = 1,NDEGI
               CALL DAXPY(NBASI,UMAT(JDEG,IDEG),
     *                    CMOI(1,IMO),1,WRK(1,JDEG),1)
  220       CONTINUE
  240    CONTINUE
         DO 260 IDEG = 1,NDEGI
            IMO = IMODEG(IDEG,KMOSS) - IOFFMO
            CALL DCOPY(NBASI,WRK(1,IDEG),1,CMOI(1,IMO),1)
  260    CONTINUE
  280 CONTINUE
C
      CALL QEXIT('AVEDIS')
      RETURN
      END
C  /* Deck avegtu */
      SUBROUTINE AVEGTU(NDEGI,NORBT,IMOSS,KMOSS,IMOI,IMOK,
     *                  UMAT,TKMO,UTUMAT)
C
#include <implicit.h>
      DIMENSION IMOI(NDEGI), IMOK(NDEGI)
      DIMENSION UMAT(NDEGI,NDEGI), TKMO(NORBT,NORBT)
      DIMENSION UTUMAT(NDEGI,NDEGI)
      PARAMETER (D1 = 1.0D0)
C
C Used from common blocks:
C   INFINP : THRSSY
C   INFPRI : IPRAVE
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infpri.h>
C
      CALL QENTER('AVEGTU')
      DO 120 K = 1,NDEGI
         DO 110 I = 1,NDEGI
            UMAT(I,K) = TKMO( IMOI(I) , IMOK(K) )
  110    CONTINUE
  120 CONTINUE
C
C  -- Do some printout for debuggers
C
      IF (IPRAVE .GE. 10) THEN
         WRITE(LUPRI,'(/A/A/A,3I5)')
     *      ' Test output from AVEGTU',
     *      ' =======================',
     *      ' IMOSS, KMOSS, NDEGI =',IMOSS,KMOSS,NDEGI
         WRITE(LUPRI,'(A,20I4)') ' IMOI :',(IMOI(I),I=1,NDEGI)
         WRITE(LUPRI,'(A,20I4)') ' IMOK :',(IMOK(I),I=1,NDEGI)
         WRITE(LUPRI,'(/A/A/)') ' UMAT found in AVEGTU',
     *                         ' ===================='
         CALL OUTPUT(UMAT,1,NDEGI,1,NDEGI,NDEGI,NDEGI,1,LUPRI)
      END IF
C
      IF (IMOSS .EQ. KMOSS) THEN
C     ... check that UMAT is diagonal
         NERR = NOFDIA(NDEGI,NDEGI,UMAT,THRSSY)
         IF (NERR .GT. 0) THEN
            WRITE(LUPRI,'(2A,I5)') ' AVEGTU error, ',
     *         'no. of non-zero off-diagonal elements in U =',NERR
            WRITE(LUPRI,'(A,I4)') ' -- IMOSS = ', IMOSS
            WRITE(LUPRI,'(A,I4)') ' -- KMOSS = ', KMOSS
            WRITE(LUPRI,'(/A)') ' The U matrix :'
            CALL OUTPUT(UMAT,1,NDEGI,1,NDEGI,NDEGI,NDEGI,1,LUPRI)
            CALL QUIT('AVEGTU error for IMOSS .eq. KMOSS')
         END IF
      ELSE
C     ... IMOSS .ne. KMOSS
C         a) check that matrix is x*orthogonal matrix
         ITURN = 0
  300    ITURN = ITURN + 1
         CALL DGEMM('T','N',NDEGI,NDEGI,NDEGI,1.D0,
     &              UMAT,NDEGI,
     &              UMAT,NDEGI,0.D0,
     &              UTUMAT,NDEGI)
         NERR = NOFDIA(NDEGI,NDEGI,UTUMAT,THRSSY)
         IF (NERR .GT. 0 .OR. IPRAVE .GE. 10) THEN
            IF (NERR .GT. 0) WRITE(LUPRI,'(2A,I5)') ' AVEGTU error, ',
     *         'no. of non-zero off-diagonal elements in U^(T)U =',NERR
            WRITE(LUPRI,'(A,I4)') ' -- IMOSS = ', IMOSS
            WRITE(LUPRI,'(A,I4)') ' -- KMOSS = ', KMOSS
            WRITE(LUPRI,'(/A)') ' The U matrix :'
            CALL OUTPUT(UMAT,1,NDEGI,1,NDEGI,NDEGI,NDEGI,1,LUPRI)
            WRITE(LUPRI,'(/A)') ' The UTU matrix :'
            CALL OUTPUT(UTUMAT,1,NDEGI,1,NDEGI,NDEGI,NDEGI,1,LUPRI)
            IF (NERR .GT. 0) 
     &           CALL QUIT('AVEGTU error for IMOSS .ne. KMOSS')
         END IF
C
C HJMAERKE 940818-hjaaj: we ought to force UMAT orthogonal with
C   symmetric orthonormalization; with the current code we may
C   get deviations from orthonormality up to THRSSY
         X2 = UTUMAT(1,1)
         DO 310 I = 2,NDEGI
            X2 = X2 * UTUMAT(I,I)
  310    CONTINUE
C        det(U) = sqrt(X2) = X2 ** ( .5 )
         XEXP = 2*NDEGI
         XEXP = - D1 / XEXP
         XINV = X2 ** ( XEXP )
         CALL DSCAL(NDEGI*NDEGI,XINV,UMAT,1)
         IF (ABS(D1-XINV) .GT. THRSSY) THEN
            IF (IPRAVE .GE. 10 .OR. ITURN .GE. 3) THEN
               WRITE(LUPRI,'(A,I4)') ' -- ITURN = ', ITURN
               WRITE(LUPRI,'(A,1D20.10)') ' -- XINV  = ', XINV
               WRITE(LUPRI,'(A,I4)') ' -- IMOSS = ', IMOSS
               WRITE(LUPRI,'(A,I4)') ' -- KMOSS = ', KMOSS
               WRITE(LUPRI,'(/A)') ' The U matrix :'
               CALL OUTPUT(UMAT,1,NDEGI,1,NDEGI,NDEGI,NDEGI,1,LUPRI)
            END IF
            IF (ITURN .LT. 3) GO TO 300
            WRITE(LUPRI,'(A)')
     *      ' AVEGTU error : failed to normalize UMAT'
            CALL QUIT('AVEGTU error: could not normalize UMAT')
         END IF
      END IF
      CALL QEXIT('AVEGTU')
      RETURN
      END
C  /* Deck avepha */
      SUBROUTINE AVEPHA(ISSYMR,CMO,TKMO,IPHCHA)
C
C  2-Jul-1992 hjaaj+hh
C  Purpose: make relative phases the same in each block
C           of a set of degenerate mo.s.
C           A block corresponds to orbitals following
C           a row of an irrep.
C
#include <implicit.h>
      DIMENSION CMO(NCMOT), TKMO(NORBT,NORBT)
C
C Used from common blocks:
C   INFINP : THRSSY
C   INFORB : NCMOT,NORBT,NORBSS(), IORBSS(),...,NINFSS(*,3),...
C   INFIND : ISSORD()
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      CALL QENTER('AVEPHA')
      IF (NINFSS(ISSYMR,2) .NE. ISSYMR) THEN
         WRITE(LUPRI,'(/A)')
     *      ' AVEPHA error: ISSYMR is not a "root" symmetry'
         WRITE(LUPRI,'(A,I4)') ' ISSYMR           = ', ISSYMR
         WRITE(LUPRI,'(A,I4)') ' NINFSS(ISSYMR,2) = ', NINFSS(ISSYMR,2)
         CALL QUIT('AVEPHA error: ISSYMR is not a "root" symmetry')
      END IF
      IOSSI = IORBSS(ISSYMR)
      NOSSI = NORBSS(ISSYMR)
      DO 800 JSSYM = 1,NSSYM
         JSSYMR = NINFSS(JSSYM,2)
         IF (JSSYMR .EQ. ISSYMR) THEN
            IF (JSSYM .LT. ISSYMR) THEN
               CALL QUIT('AVEPHA error: JSSYM .lt ISSYMR')
            END IF
            IF (NORBSS(JSSYM) .NE. NOSSI) THEN
               WRITE(LUPRI,'(/A/A,2I5/A)')
     *            ' AVEHPA error: Diff. no. of orb.s in deg. supsym.s',
     *            ' JSSYMR, JSSYM =',JSSYMR,JSSYM,
     *            ' Dump of super symmetry information :'
               CALL AVEDMP(LUPRI)
               CALL QUIT(
     &         'AVEHPA error: Diff. no. of orb.s in deg. supsym.s')
            END IF
            ISYMJ = NINFSS(JSSYM,3)
            ICMOJ = ICMO(ISYMJ)
            NBASJ = NBAS(ISYMJ)
            IORBJ = IORB(ISYMJ)
            IOSSJ = IORBSS(JSSYM)
            IMO1  = ISSORD(IOSSI + 1)
            JMO1  = ISSORD(IOSSJ + 1)
            DO 700 K = 1,NOSSI
               IMO = ISSORD(IOSSI + K)
               JMO = ISSORD(IOSSJ + K)
               IF (ABS(TKMO(IMO,IMO) - TKMO(JMO,JMO)) .GT. THRSSY) THEN
                  WRITE (LUPRI,'(/A,2I5)')
     &               ' AVEPHA error: different ordering'//
     &               ' of degenerate orbitals in supsyms:',ISSYMR,JSSYM
                  CALL QUIT(
     &            'AVEPHA error: different ordering of deg. orb.s')
               END IF
               IF      (ABS(TKMO(IMO,IMO1)) .LT. THRSSY) THEN
                  WRITE(LUPRI,'(/A/A,2I5/A,D10.2)')
     &            ' AVEPHA error: TKMO(IMO,IMO1) is zero',
     &            ' IMO, IMO1      =',IMO,IMO1,
     &            ' TKMO(IMO,IMO1) =',TKMO(IMO,IMO1)
                  CALL QTRACE(LUERR)
                  CALL QUIT('AVEPHA error: TKMO(IMO,IMO1) is zero')
               ELSE IF (ABS(TKMO(IMO,IMO1)+TKMO(JMO,JMO1)) .LT. THRSSY)
     *         THEN
                  IPHCHA = IPHCHA + 1
                  IF (IPRAVE .GE. 4) THEN
                     WRITE(LUPRI,'(A,I4)')
     *                  ' AVEPHA changing phase of mo no.',JMO
                  END IF
                  JOFF = ICMOJ + NBASJ*(JMO-IORBJ-1)
                  DO 600 J = 1,NBASJ
                     CMO(JOFF+J) = -CMO(JOFF+J)
  600             CONTINUE
               ELSE IF (ABS(TKMO(IMO,IMO1)-TKMO(JMO,JMO1)) .GT. THRSSY)
     *         THEN
                  WRITE(LUPRI,'(/A,2(/A,2I5/A,D20.10))')
     &            ' AVEPHA error: TKMO(IMO,IMO1) not TKMO(JMO,JMO1)',
     &            ' IMO, IMO1      =',IMO,IMO1,
     &            ' TKMO(IMO,IMO1) =',TKMO(IMO,IMO1),
     &            ' JMO, JMO1      =',JMO,JMO1,
     &            ' TKMO(JMO,JMO1) =',TKMO(JMO,JMO1)
                  CALL QTRACE(LUERR)
                  CALL QUIT('AVEPHA error: error in deg. orb.s code 2')
               END IF
  700       CONTINUE
         END IF
  800 CONTINUE
      CALL QEXIT('AVEPHA')
      RETURN
      END
C  /* Deck aveudv */
#if defined (VAR_NOTUSED  )
#include <scalsub1.h>
      SUBROUTINE AVEUDV(UDV)
C
C     09-Jul-1992 Hinne Hettema
C
C     Purpose : perform an averaging of the one-electron density
C     matrix. First we identify an IMO which is degenerate with
C     other MO's, then we average the matrix rows and colums.
C
#include <implicit.h>
#include <scalsub2.h>
      DIMENSION UDV(NASHT,*)
C
C     Used from common blocks
C     INFORB    : NSYM ,NSSYM, NORBSS(*), IORBSS(*), NINFSS(*,3)
C     INFIND    : ICH(*), ISSMO(*), ISSORD(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      PARAMETER ( MAXDEG = 20, D0 = 0.0D0, D1 = 1.0D0 )
      DIMENSION ISSOFF(MAXDEG),IDEGV(MAXDEG,2)
C
      CALL QENTER('AVEUDV')
C
      DO 800 ISSYMR = 1, NSSYM
C     -- skip if supersym not degenerate or not 'root' supersym
         IF ( NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF ( NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
C
         NDEG  = 0
         DO 100 JSSYM = ISSYMR, NSSYM
C        -- skip if not same 'root' supersymmetry
C           (because then not degenerate at all)
            IF ( NINFSS(JSSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISSOFF(NDEG) = IORBSS(JSSYM) - IORBSS(ISSYMR)
            END IF
  100    CONTINUE
         IF ( NDEG .NE. NINFSS(ISSYMR,1) ) THEN
            WRITE(LUERR,'(A)') ' AVEUDV error : inconsistent'//
     &            ' degeneracies'
            WRITE(LUERR,'(A,I4)') ' - ISSYMR           = ', ISSYMR
            WRITE(LUERR,'(A,I4)') ' - NDEG             = ', NDEG
            WRITE(LUERR,'(A,I4)') ' - NINFSS(ISSYMR,1) = ',
     &                 NINFSS(ISSYMR,1)
            CALL QUIT(' AVEUDV : inconsistent degeneracy in NINFSS')
         END IF
C
C     -- Compute averaging factor
         AVFAC = NDEG
         AVFAC = D1 / AVFAC
C
C     -- in the following, use supersymmetry info to get IMO
         IRSYM  = NINFSS(ISSYMR,3)
         NASHI  = NASH(IRSYM)
         IOFFSR = IORBSS(ISSYMR)
         IOFFDV = IASH(IRSYM) + 1
         DO 700 IMOSS = IOFFSR+1, IOFFSR+NORBSS(ISSYMR)
            IMO  = ISSORD(IMOSS)
            IMOW = ICH(IMO)
            IF (IMOW .LE. 0) GO TO 700
            DO 600 JMOSS = IMOSS, IOFFSR+NORBSS(ISSYMR)
               JMO = ISSORD(JMOSS)
               JMOW = ICH(JMO)
               IF ( JMOW .LE. 0 ) GO TO 600
               DO 300 IDEG = 2, NDEG
                  KMO = ISSORD(IMOSS + ISSOFF(IDEG))
                  KMOW = ICH(KMO)
                  IF (KMOW .LE. 0 ) THEN
                     CALL QUIT('AVEUDV code 300.1: ICH(KMO).le.0')
                  END IF
                  LMO = ISSORD(JMOSS + ISSOFF(IDEG))
                  LMOW = ICH(LMO)
                  IF (LMOW .LE. 0 ) THEN
                     CALL QUIT('AVEUDV code 300.2: ICH(LMO).le.0')
                  END IF
                  IDEGV(IDEG,1) = KMOW
                  IDEGV(IDEG,2) = LMOW
  300          CONTINUE
C           -- now collect all elements which should be made equal
               XUDVIJ = UDV(IMOW,JMOW)
               DO 400 IDEG = 2, NDEG
                  KMOW = IDEGV(IDEG,1)
                  LMOW = IDEGV(IDEG,2)
                  XUDVIJ = XUDVIJ + UDV(KMOW,LMOW)
  400          CONTINUE
               XUDVIJ = XUDVIJ * AVFAC
C           -- now distribute all elements which should be made equal
               UDV(IMOW,JMOW) = XUDVIJ
               IF (IMOW .NE. JMOW) UDV(JMOW,IMOW) = XUDVIJ
               DO 500 IDEG = 2, NDEG
                  KMOW = IDEGV(IDEG,1)
                  LMOW = IDEGV(IDEG,2)
                  UDV(KMOW,LMOW) = XUDVIJ
                  IF(KMOW .NE. LMOW) UDV(LMOW,KMOW) = XUDVIJ
  500          CONTINUE
  600       CONTINUE
  700    CONTINUE
C
  800 CONTINUE
C
C     Zero elements which are not totally symmetric
C     with respect to the supersymmetry.
C
      DO 980 ISYM = 1,NSYM
         IASHI = IASH(ISYM)
         NASHI = NASH(ISYM)
         DO 960 IMOW = IASHI+1,IASHI+NASHI
            IMO = ISX(NISHT+IMOW)
            DO 950 JMOW = IMOW+1,IASHI+NASHI
               JMO = ISX(NISHT+JMOW)
               IF (ISSMO(JMO) .NE. ISSMO(IMO)) THEN
                  UDV(IMOW,JMOW) = D0
                  UDV(JMOW,IMOW) = D0
               END IF
  950       CONTINUE
  960    CONTINUE
  980 CONTINUE
C
      CALL QEXIT('AVEUDV')
C
      RETURN
      END
#endif
C  /* Deck aveh1 */
#include <scalsub1.h>
      SUBROUTINE AVEH1(H1)
C
C     12-Jul-1992 Hinne Hettema
C
C     Purpose : perform an averaging of a one-electron
C     matrix. First we identify an IMO which is degenerate with
C     other MO's, then we average the matrix rows and colums.
C
#include <implicit.h>
#include <scalsub2.h>
      DIMENSION H1(NORBT,*)
C
C     Used from common blocks
C     INFORB    : NSYM ,NSSYM, NORBSS(*), IORBSS(*), NINFSS(*,3)
C     INFIND    : ICH(*), ISSMO(*), ISSORD(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      PARAMETER ( MAXDEG = 20, D0 = 0.0D0, D1 = 1.0D0 )
      DIMENSION ISSOFF(MAXDEG),IDEGV(MAXDEG,2)
C
      CALL QENTER('AVEH1 ')
C
      IF ( IPRAVE .GE. 10 ) THEN
         WRITE(LUPRI,'(//A)') ' H1 matrix before averaging'
         WRITE(LUPRI,'(A)')   ' =========================='
         CALL OUTPUT(H1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      DO 800 ISSYMR = 1, NSSYM
C     -- skip if supersym not degenerate or not 'root' supersym
         IF ( NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF ( NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
C
         NDEG  = 0
         DO 100 JSSYM = ISSYMR, NSSYM
C        -- skip if not same 'root' supersymmetry
C           (because then not degenerate at all)
            IF ( NINFSS(JSSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISSOFF(NDEG) = IORBSS(JSSYM) - IORBSS(ISSYMR)
            END IF
  100    CONTINUE
         IF ( NDEG .NE. NINFSS(ISSYMR,1) ) THEN
            WRITE(LUERR,'(A)') ' AVEH1 error : inconsistent'//
     &            ' degeneracies'
            WRITE(LUERR,'(A,I4)') ' - ISSYMR           = ', ISSYMR
            WRITE(LUERR,'(A,I4)') ' - NDEG             = ', NDEG
            WRITE(LUERR,'(A,I4)') ' - NINFSS(ISSYMR,1) = ',
     &                 NINFSS(ISSYMR,1)
            CALL QUIT(' AVEH1 : inconsistent degeneracy in NINFSS')
         END IF
C
C     -- Compute averaging factor
         AVFAC = NDEG
         AVFAC = D1 / AVFAC
C
C     -- in the following, use supersymmetry info to get IMO
         IRSYM  = NINFSS(ISSYMR,3)
         NORBI  = NORB(IRSYM)
         IOFFSR = IORBSS(ISSYMR)
         DO 700 IMOSS = IOFFSR+1, IOFFSR+NORBSS(ISSYMR)
            IMO  = ISSORD(IMOSS)
            DO 600 JMOSS = IMOSS, IOFFSR+NORBSS(ISSYMR)
               JMO = ISSORD(JMOSS)
               IDEGV(1,1) = IMO
               IDEGV(1,2) = JMO
               DO 300 IDEG = 2, NDEG
                  IDEGV(IDEG,1) = ISSORD(IMOSS + ISSOFF(IDEG))
                  IDEGV(IDEG,2) = ISSORD(JMOSS + ISSOFF(IDEG))
  300          CONTINUE
C           -- now collect all elements which should be made equal
               XH1IJ = H1(IMO,JMO)
               DO 400 IDEG = 2, NDEG
                  KMO = IDEGV(IDEG,1)
                  LMO = IDEGV(IDEG,2)
                  XH1IJ = XH1IJ + H1(KMO,LMO)
  400          CONTINUE
               XH1IJ = XH1IJ * AVFAC
C           -- now distribute all elements which should be made equal
               H1(IMO,JMO) = XH1IJ
               IF (IMO .NE. JMO) H1(JMO,IMO) = XH1IJ
               DO 500 IDEG = 2, NDEG
                  KMO = IDEGV(IDEG,1)
                  LMO = IDEGV(IDEG,2)
                  H1(KMO,LMO) = XH1IJ
                  IF(KMO .NE. LMO) H1(LMO,KMO) = XH1IJ
  500          CONTINUE
  600       CONTINUE
  700    CONTINUE
C
  800 CONTINUE
C
C     Zero elements which are not totally symmetric
C     with respect to the supersymmetry.
C
      DO 980 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         DO 960 IMO = IORBI+1,IORBI+NORBI
            DO 950 JMO = IMO+1,IORBI+NORBI
               IF (ISSMO(JMO) .NE. ISSMO(IMO)) THEN
                  H1(IMO,JMO) = D0
                  H1(JMO,IMO) = D0
               END IF
  950       CONTINUE
  960    CONTINUE
  980 CONTINUE
C
      IF ( IPRAVE .GE. 10 ) THEN
         WRITE(LUPRI,'(//A)') ' H1 matrix after averaging'
         WRITE(LUPRI,'(A)')   ' ========================='
         CALL OUTPUT(H1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      CALL QEXIT('AVEH1 ')
C
      RETURN
      END
C  /* Deck avedv */
#include <scalsub1.h>
      SUBROUTINE AVEDV(DV)
C
C     09-Jul-1992 Hinne Hettema
C
C     Purpose : perform an averaging of the one-electron density
C     matrix. First we identify an IMO which is degenerate with
C     other MO's, then we average the matrix rows and colums.
C
#include <implicit.h>
#include <scalsub2.h>
      DIMENSION DV(NNASHX)
C
C     Used from common blocks
C     INFORB    : NSYM ,NSSYM, NORBSS(*), IORBSS(*), NINFSS(*,3)
C     INFIND    : IROW(*), ICH(*), ISSMO(*), ISSORD(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      PARAMETER ( MAXDEG = 20, D0 = 0.0D0, D1 = 1.0D0 )
      DIMENSION ISSOFF(MAXDEG),IDEGV(MAXDEG)
C
      CALL QENTER('AVEDV ')
C
      IF (IPRAVE .GE. 5) THEN
         WRITE(LUPRI,'(/A)')
     *      ' AVEDV test output: DV matrix before averaging'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF
C
      DO 800 ISSYMR = 1, NSSYM
C     -- skip if supersym not degenerate or not 'root' supersym
         IF ( NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF ( NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
C
         NDEG  = 0
         DO 100 JSSYM = ISSYMR, NSSYM
C        -- skip if not same 'root' supersymmetry
C           (because then not degenerate at all)
            IF ( NINFSS(JSSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISSOFF(NDEG) = IORBSS(JSSYM) - IORBSS(ISSYMR)
            END IF
  100    CONTINUE
         IF ( NDEG .NE. NINFSS(ISSYMR,1) ) THEN
            WRITE(LUERR,'(A)') ' AVEDV error : inconsistent'//
     &            ' degeneracies'
            WRITE(LUERR,'(A,I4)') ' - ISSYMR           = ', ISSYMR
            WRITE(LUERR,'(A,I4)') ' - NDEG             = ', NDEG
            WRITE(LUERR,'(A,I4)') ' - NINFSS(ISSYMR,1) = ',
     &                 NINFSS(ISSYMR,1)
            CALL QUIT(' AVEDV : inconsistent degeneracy in NINFSS')
         END IF
C
C     -- Compute averaging factor
         AVFAC = NDEG
         AVFAC = D1 / AVFAC
C
C     -- in the following, use supersymmetry info to get IMO
         IRSYM  = NINFSS(ISSYMR,3)
         NASHI  = NASH(IRSYM)
         IOFFSR = IORBSS(ISSYMR)
         IOFFDV = IASH(IRSYM) + 1
         DO 700 IMOSS = IOFFSR+1, IOFFSR+NORBSS(ISSYMR)
            IMO  = ISSORD(IMOSS)
            IMOW = ICH(IMO)
            IF (IMOW .LE. 0) GO TO 700
            DO 600 JMOSS = IOFFSR+1, IMOSS
               JMO = ISSORD(JMOSS)
               JMOW = ICH(JMO)
               IF ( JMOW .LE. 0 ) GO TO 600
               IDEGV(1) = IROW(IMOW) + JMOW
               DO 300 IDEG = 2, NDEG
                  KMO = ISSORD(IMOSS + ISSOFF(IDEG))
                  KMOW = ICH(KMO)
                  IF (KMOW .LE. 0 ) THEN
                     CALL QUIT('AVEDV code 300.1')
                  END IF
                  LMO = ISSORD(JMOSS + ISSOFF(IDEG))
                  LMOW = ICH(LMO)
                  IF (LMOW .LE. 0 ) THEN
                     CALL QUIT('AVEDV code 300.2')
                  END IF
                  IDEGV(IDEG) = IROW(KMOW) + LMOW
  300          CONTINUE
C           -- now collect all elements which should be made equal
               XDVIJ = DV(IDEGV(1))
               DO 400 IDEG = 2, NDEG
                  XDVIJ = XDVIJ + DV(IDEGV(IDEG))
  400          CONTINUE
               XDVIJ = XDVIJ * AVFAC
C           -- now distribute all elements which should be made equal
               DO 500 IDEG = 1, NDEG
                  DV(IDEGV(IDEG)) = XDVIJ
  500          CONTINUE
  600       CONTINUE
  700    CONTINUE
C
  800 CONTINUE
C
C     Zero elements which are not totally symmetric
C     with respect to the supersymmetry.
C
      DO 980 ISYM = 1,NSYM
         IASHI = IASH(ISYM)
         NASHI = NASH(ISYM)
         DO 960 IMOW = IASHI+2,IASHI+NASHI
            IMO = ISX(NISHT+IMOW)
            DO 950 JMOW = IASHI+1,IMOW-1
               JMO = ISX(NISHT+JMOW)
               IF (ISSMO(JMO) .NE. ISSMO(IMO)) THEN
                  IDV = IROW(IMOW) + JMOW
                  DV(IDV) = D0
               END IF
  950       CONTINUE
  960    CONTINUE
  980 CONTINUE
C
      IF (IPRAVE .GE. 5) THEN
         WRITE(LUPRI,'(/A)')
     *      ' AVEDV test output: DV matrix after averaging'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF
      CALL QEXIT('AVEDV ')
C
      RETURN
      END
C  /* Deck aveh1p */
#if defined (VAR_NOTUSED  )
#include <scalsub1.h>
      SUBROUTINE AVEH1P(H1P)
C
C     12-Jul-1992 Hinne Hettema
C
C     Purpose : perform an averaging of a one-electron
C     matrix. First we identify an IMO which is degenerate with
C     other MO's, then we average the matrix rows and colums.
C
#include <implicit.h>
#include <scalsub2.h>
      DIMENSION H1P(NNORBX)
C
C     Used from common blocks
C     INFORB    : NSYM ,NSSYM, NORBSS(*), IORBSS(*), NINFSS(*,3)
C     INFIND    : ISSMO(*), ISSORD(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      PARAMETER ( MAXDEG = 20, D0 = 0.0D0, D1 = 1.0D0 )
      DIMENSION ISSOFF(MAXDEG),IDEGV(MAXDEG)
C
      CALL QENTER('AVEH1P')
C
      IF (IPRAVE .GE. 5) THEN
         WRITE(LUPRI,'(/A)')
     *      ' AVEH1P test output: H1P matrix before averaging'
         CALL OUTPAK(H1P,NORBT,LUPRI)
      END IF
C
      DO 800 ISSYMR = 1, NSSYM
C     -- skip if supersym not degenerate or not 'root' supersym
         IF ( NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF ( NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
C
         NDEG  = 0
         DO 100 JSSYM = ISSYMR, NSSYM
C        -- skip if not same 'root' supersymmetry
C           (because then not degenerate at all)
            IF ( NINFSS(JSSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISSOFF(NDEG) = IORBSS(JSSYM) - IORBSS(ISSYMR)
            END IF
  100    CONTINUE
         IF ( NDEG .NE. NINFSS(ISSYMR,1) ) THEN
            WRITE(LUERR,'(A)') ' AVEH1P error : inconsistent'//
     &            ' degeneracies'
            WRITE(LUERR,'(A,I4)') ' - ISSYMR           = ', ISSYMR
            WRITE(LUERR,'(A,I4)') ' - NDEG             = ', NDEG
            WRITE(LUERR,'(A,I4)') ' - NINFSS(ISSYMR,1) = ',
     &                 NINFSS(ISSYMR,1)
            CALL QUIT(' AVEH1P : inconsistent degeneracy in NINFSS')
         END IF
C
C     -- Compute averaging factor
         AVFAC = NDEG
         AVFAC = D1 / AVFAC
C
C     -- in the following, use supersymmetry info to get IMO
         IRSYM  = NINFSS(ISSYMR,3)
         NORBI  = NORB(IRSYM)
         IOFFSR = IORBSS(ISSYMR)
         DO 700 IMOSS = IOFFSR+1, IOFFSR+NORBSS(ISSYMR)
            IMO  = ISSORD(IMOSS)
            DO 600 JMOSS = IOFFSR+1, IMOSS
               JMO = ISSORD(JMOSS)
               IDEGV(1) = IROW(IMO) + JMO
               DO 300 IDEG = 2, NDEG
C                 920714-hjaaj: KMO var. is used to be safe
C                 xlf version 2.2 has made errors in a similar statem.
C                 to the version with IROW(IMO+ISSOFF(IDEG))
                  KMO = ISSORD(IMOSS + ISSOFF(IDEG))
                  IDEGV(IDEG) = IROW(KMO) + ISSORD(JMOSS + ISSOFF(IDEG))
  300          CONTINUE
C           -- now collect all elements which should be made equal
               XH1PIJ = H1P(IDEGV(1))
               DO 400 IDEG = 2, NDEG
                  XH1PIJ = XH1PIJ + H1P(IDEGV(IDEG))
  400          CONTINUE
               XH1PIJ = XH1PIJ * AVFAC
C           -- now distribute all elements which should be made equal
               DO 500 IDEG = 1, NDEG
                  H1P(IDEGV(IDEG)) = XH1PIJ
  500          CONTINUE
  600       CONTINUE
  700    CONTINUE
C
  800 CONTINUE
C
C     Zero elements which are not totally symmetric
C     with respect to the supersymmetry.
C
      DO 980 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         DO 960 IMO = IORBI+2,IORBI+NORBI
            DO 950 JMO = IORBI+1,IMO-1
               IF (ISSMO(JMO) .NE. ISSMO(IMO)) THEN
                  IJH1P = IROW(IMO) + JMO
                  H1P(IJH1P) = D0
               END IF
  950       CONTINUE
  960    CONTINUE
  980 CONTINUE
C
      IF (IPRAVE .GE. 5) THEN
         WRITE(LUPRI,'(/A)')
     *      ' AVEH1P test output: H1P matrix after averaging'
         CALL OUTPAK(H1P,NORBT,1,LUPRI)
      END IF
C
      CALL QEXIT('AVEH1P')
C
      RETURN
      END
#endif
C  /* Deck avefck */
#include <scalsub1.h>
      SUBROUTINE AVEFCK(FCK)
C
C     12-Jul-1992 Hinne Hettema
C
C     Purpose : perform an averaging of a one-electron
C     matrix. First we identify an IMO which is degenerate with
C     other MO's, then we average the matrix rows and colums.
C
#include <implicit.h>
#include <scalsub2.h>
      DIMENSION FCK(NNORBT)
C
C     Used from common blocks
C     INFORB    : NSYM ,NSSYM, NORBSS(*), IORBSS(*), NINFSS(*,3)
C     INFIND    : ISSMO(*), ISSORD(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      PARAMETER ( MAXDEG = 20, D0 = 0.0D0, D1 = 1.0D0 )
      DIMENSION ISMDEG(MAXDEG),ISSOFF(MAXDEG),IDEGV(MAXDEG)
C
      CALL QENTER('AVEFCK')
C
      IF (IPRAVE .GE. 7) THEN
         WRITE(LUPRI,'(/A)')
     *      ' AVEFCK test output: FCK matrix before averaging'
         CALL OUTPKB(FCK,NORB,NSYM,1,LUPRI)
      END IF
C
      DO 800 ISSYMR = 1, NSSYM
C     -- skip if supersym not degenerate or not 'root' supersym
         IF ( NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF ( NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
C
         NDEG  = 0
         DO 100 JSSYM = ISSYMR, NSSYM
C        -- skip if not same 'root' supersymmetry
C           (because then not degenerate at all)
            IF ( NINFSS(JSSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISMDEG(NDEG) = NINFSS(JSSYM,3)
               ISSOFF(NDEG) = IORBSS(JSSYM) - IORBSS(ISSYMR)
            END IF
  100    CONTINUE
         IF ( NDEG .NE. NINFSS(ISSYMR,1) ) THEN
            WRITE(LUERR,'(A)') ' AVEFCK error : inconsistent'//
     &            ' degeneracies'
            WRITE(LUERR,'(A,I4)') ' - ISSYMR           = ', ISSYMR
            WRITE(LUERR,'(A,I4)') ' - NDEG             = ', NDEG
            WRITE(LUERR,'(A,I4)') ' - NINFSS(ISSYMR,1) = ',
     &                 NINFSS(ISSYMR,1)
            CALL QUIT(' AVEFCK : inconsistent degeneracy in NINFSS')
         END IF
C
C     -- Compute averaging factor
         AVFAC = NDEG
         AVFAC = D1 / AVFAC
C
C     -- in the following, use supersymmetry info to get IMO
         IRSYM  = ISMDEG(1)
         IORBI  = IORB(IRSYM)
         NORBI  = NORB(IRSYM)
         IIORBI = IIORB(IRSYM)
         IOFFSR = IORBSS(ISSYMR)
         DO 700 IMOSS = IOFFSR+1, IOFFSR+NORBSS(ISSYMR)
            IMO  = ISSORD(IMOSS) - IORBI
            DO 600 JMOSS = IOFFSR+1, IMOSS
               JMO = ISSORD(JMOSS) - IORBI
               IDEGV(1) = IIORBI + IROW(IMO) + JMO
               DO 300 IDEG = 2, NDEG
                  JRSYM  = ISMDEG(IDEG)
                  IORBJ  = IORB(JRSYM)
                  IIORBJ = IIORB(JRSYM)
C                 920714-hjaaj: KMO var. is used to be safe
C                 xlf version 2.2 has made errors in a similar statem.
C                 to the version with IROW(IMO+ISSOFF(IDEG))
                  KMO = ISSORD(IMOSS + ISSOFF(IDEG)) - IORBJ
                  LMO = ISSORD(JMOSS + ISSOFF(IDEG)) - IORBJ
                  IDEGV(IDEG) = IIORBJ + IROW(KMO) + LMO
  300          CONTINUE
C           -- now collect all elements which should be made equal
               XFCKIJ = FCK(IDEGV(1))
               DO 400 IDEG = 2, NDEG
                  XFCKIJ = XFCKIJ + FCK(IDEGV(IDEG))
  400          CONTINUE
               XFCKIJ = XFCKIJ * AVFAC
C           -- now distribute all elements which should be made equal
               DO 500 IDEG = 1, NDEG
                  FCK(IDEGV(IDEG)) = XFCKIJ
  500          CONTINUE
  600       CONTINUE
  700    CONTINUE
C
  800 CONTINUE
C
C     Zero elements which are not totally symmetric
C     with respect to the supersymmetry.
C
      DO 980 ISYM = 1,NSYM
         IORBI  = IORB(ISYM)
         NORBI  = NORB(ISYM)
         IIORBI = IIORB(ISYM)
         DO 960 IMO = 2,NORBI
            DO 950 JMO = 1,IMO-1
               IF (ISSMO(IORBI+JMO) .NE. ISSMO(IORBI+IMO)) THEN
                  IJFCK = IIORBI + IROW(IMO) + JMO
                  FCK(IJFCK) = D0
               END IF
  950       CONTINUE
  960    CONTINUE
  980 CONTINUE
C
      IF (IPRAVE .GE. 7) THEN
         WRITE(LUPRI,'(/A)')
     *      ' AVEFCK test output: FCK matrix after averaging'
         CALL OUTPKB(FCK,NORB,NSYM,1,LUPRI)
      END IF
C
      CALL QEXIT('AVEFCK')
C
      RETURN
      END
C  /* Deck aveprt */
      SUBROUTINE AVEPRT(LUPRI)
C
C     7-Jul-92 Hinne Hettema and Hans Joergen Aa. Jensen.
C     Revised 15-Jul-1993 hjaaj
C     Purpose : write supersymmetry information
C
#include <implicit.h>
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      PARAMETER ( ISTEP = 10 )
C
      CALL QENTER('AVEPRT')
C
C     Check if supersymmetry is identical to "D2h" symmetry
C
      IF (NSSYM .LE. NSYM .AND. MXDGSS .EQ. 1) THEN
         DO 200 ISYM = 1,NSYM
            JSYM = 0
            DO 100 ISSYM = 1,NSSYM
               IF (NINFSS(ISSYM,3) .EQ. ISYM) JSYM = JSYM + 1
  100       CONTINUE
            IF (JSYM .GT. 1) GO TO 300
  200    CONTINUE
         WRITE(LUPRI,'(/5X,A)') 'No supersymmetry found.'
         GO TO 9999
      END IF
  300 CONTINUE
C
      WRITE(LUPRI,'(/5X,A)') 'Supersymmetry specification'
      WRITE(LUPRI,'( 5X,A)') '==========================='
      IF (MXDGSS .EQ. 1)
     &WRITE(LUPRI,'(  5X,A)') '(only one dimensional irreps)'
      ICOUNT = 0
      DO 500 ISSYM = 1, NSSYM, ISTEP
         ICOUNT = ICOUNT + 1
         ISTRT = ISSYM
         IEND  = MIN(NSSYM, ICOUNT * ISTEP)
         ILEN = IEND - ISTRT + 1
C
         WRITE(LUPRI,1100)
     &            'Supersymmetry number',
     &            (I, I=ISTRT, IEND)
         IF (MXDGSS .GT. 1) WRITE(LUPRI,1100)
     &            'degenerate with sup.sym.',
     &            (NINFSS(I,2), I = ISTRT, IEND)
         WRITE(LUPRI,1100)
     &            'and located in abelian symmetry ',
     &            (NINFSS(I,3), I=ISTRT, IEND)
         WRITE(LUPRI,1110)
     &            ('  --', I = 1, ILEN)
         WRITE(LUPRI,1100)
     &            'Total number of orbitals',
     &            (NORBSS(I), I = ISTRT, IEND)
         WRITE(LUPRI,'()')
  500 CONTINUE
 9999 CALL QEXIT('AVEPRT')
      RETURN
 1100 FORMAT(5X,A,T36,10I4)
 1110 FORMAT(5X,  T36,10A4)
      END
C  /* Deck avecph */
      SUBROUTINE AVECPH(IPHCHA,CMO,WRK,KFRSAV,LFRSAV)
C
C 920715-hjaaj+hh
C
C Purpose: change phases such that degenerate MO's have
C          the same relative phase.
C
#include <implicit.h>
      DIMENSION CMO(NCMOT), WRK(*)
C
C Used from common blocks:
C   INFORB : NCMOT,N2ORBX,NINFSS(*,4), NSSYM
C
#include <inforb.h>
C
      CALL QENTER('AVECPH')
      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET('REAL',KTKMO,N2ORBX,WRK,KFREE,LFREE)
      CALL AVETK(WRK(KTKMO),CMO,WRK,KFREE,LFREE)
      IPHCHA = 0
      DO 100 ISSYM = 1, NSSYM
         IDEG   = NINFSS(ISSYM,1)
         ISSYMR = NINFSS(ISSYM,2)
         IF (IDEG .GT. 1 .AND. ISSYM .EQ. ISSYMR) THEN
C        --- we have a degeneracy to other supersym.(s)
            CALL AVEPHA(ISSYMR,CMO,WRK(KTKMO),IPHCHA)
         END IF
  100 CONTINUE
      CALL MEMREL('AVECPH',WRK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('AVECPH')
      RETURN
      END
C  /* Deck avechk */
      SUBROUTINE AVECHK(WRK,LWRK,INPERR)
C
C 14-Jul-1992 Hinne Hettema and Hans Joergen Aa. Jensen.
C
C Purpose: Check supersymmetry information which has been built up
C          up till now.
C
#include <implicit.h>
C
C Used from common blocks:
C   INFINP : SUPSYM
C   INFORB : NSSYM,NORBSS(),IORBSS(),NINFSS(*,3)
C   INFIND : ISW(*),ISSORD(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infind.h>
#include <inforb.h>
#include <infpri.h>
#include <orbtypdef.h>
C
      PARAMETER (MAXDEG = 20, D1 = 1.0D0)
      DIMENSION ISSDEG(MAXDEG),ISSOFF(MAXDEG), WRK(LWRK)
C
      CALL QENTER('AVECHK')
      IF (.NOT.SUPSYM) GO TO 9999
C
      DO 800 ISSYMR = 1,NSSYM
C        Skip this supsym if not degenerate or not "root supsym"
         IF (NINFSS(ISSYMR,1) .EQ. 1) GO TO 800
         IF (NINFSS(ISSYMR,2) .NE. ISSYMR) GO TO 800
         NDEG = 0
         DO 100 ISSYM = ISSYMR,NSSYM
            IF (NINFSS(ISSYM,2) .EQ. ISSYMR) THEN
               NDEG = NDEG + 1
               ISSDEG(NDEG) = ISSYM
               ISSOFF(NDEG) = IORBSS(ISSYM) - IORBSS(ISSYMR)
            END IF
  100    CONTINUE
         IF (NDEG .NE. NINFSS(ISSYMR,1)) THEN
            INPERR = INPERR + 1
            WRITE(LUPRI,'(/A/A,I4/A,I4,A,I4)')
     &      ' Super symmetry error: inconsistent degeneracy in NINFSS',
     &      '    for "root" super symmetry',ISSYMR,
     &      '    Found:',NDEG,';   Expected:',NINFSS(ISSYMR,1)
         END IF
C
         IOFFSR = IORBSS(ISSYMR)
         ITYPA = 0
         DO 480 IMOSS = IOFFSR+1,IOFFSR+NORBSS(ISSYMR)
            IMO = ISSORD(IMOSS)
            ITYP = IOBTYP(IMO)
            IF (ITYP .EQ. JTACT) THEN
               ITYPA = IACTYP( ICH(IMO) )
            END IF
            DO 410 IDEG = 2,NDEG
               KMO = ISSORD(IMOSS + ISSOFF(IDEG))
               KTYP = IOBTYP(KMO)
               IF ( ITYP .NE. KTYP ) THEN
                  INPERR = INPERR + 1
                  WRITE (LUPRI,'(/A,2I5/A,2X,A,2X,A/A,I4)')
     &         ' Super symmetry error: degenerate orbitals',IMO,KMO,
     &         '    are not of same type',COBTYP(ITYP),COBTYP(KTYP),
     &         '    for "root" super symmetry',ISSYMR
               ELSE IF (ITYP .EQ. JTACT) THEN
                  KTYPA = IACTYP( ICH(KMO) )
                  IF (ITYPA .NE. KTYPA) THEN
                     INPERR = INPERR + 1
                     WRITE (LUPRI,'(/A,2I5/A,2I5/A,I4)')
     &      ' Super symmetry error: degenerate active orbitals',IMO,KMO,
     &      '    are not in same RAS shell type',ITYPA,KTYPA,
     &      '    for "root" super symmetry',ISSYMR
                  END IF
               END IF
  410       CONTINUE
  480    CONTINUE
  800 CONTINUE
C
 9999 CONTINUE
      CALL QEXIT('AVECHK')
      RETURN
      END
C  /* Deck avedmp */
      SUBROUTINE AVEDMP(LUWDMP)
C
C     2-Jul-92 Hinne Hettema
C     Purpose : dump supersymmetry information to output for debugging
C
#include <implicit.h>
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('AVEDMP')
      WRITE(LUWDMP,'(/A/A)')
     *   '    Dump of super symmetry information',
     *   ' ========================================'
      WRITE(LUWDMP,'(/A,I4)')
     *   ' Number of supersymmetries (NSSYM) is ',NSSYM
      WRITE(LUWDMP,'(/A/A)')
     *   '  ISS    NORBSS    IORBSS',
     *   ' =============================='
      DO 100 I = 1, NSSYM
         WRITE(LUWDMP,'(1X,I4,2I10)') I, NORBSS(I), IORBSS(I)
  100 CONTINUE
      WRITE(LUWDMP,'(A)') ' =============================='
      WRITE(LUWDMP,'(/A/A)')
     *   '    I     ISSMO    ISSORD      ISMO',
     *   ' ========================================'
      DO 200 I = 1, NORBT
         WRITE(LUWDMP,'(1X,I4,3I10)') I, ISSMO(I), ISSORD(I), ISMO(I)
  200 CONTINUE
      WRITE(LUWDMP,'(A/)') ' ========================================'
      CALL QEXIT('AVEDMP')
      RETURN
      END
