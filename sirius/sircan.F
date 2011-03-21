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
C  /* Deck sircan */
      SUBROUTINE SIRCAN (IPRCAN,CREF,CMO,OCC,WRK,KFRSAV,LFRSAV)
C
C Written  25-Nov-1984 by Hans Agren, finished  2-Jun-1985 hjaaj
C Revision 2 finished 19-Jul-1986 hjaaj (FD=FC+FV, TRACI)
C
C Purpose:
C   SIRCAN directs the calculation of canonical or natural orbitals and
C   natural occupation numbers from the MCSCF wave function
C
C Output:
C   CREF : CI-vector in natural orbitals
C   CMO  : canonical/natural orbitals
C   OCC  : occupation numbers for natural orbitals
C
#include <implicit.h>
      DIMENSION CREF(*), CMO(*), OCC(*), WRK(*)
C
C -- local variables and constants
C
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0 )
#include <dummy.h>
C
C Used from common blocks:
C   INFINP : ISTATE,NROOTS,THRGRD,FLAG(*),LSYM,?
C   INFORB : NNASHX,NASHT,?
C   INFVAR : NCONF,?
C   INFDIM : NBASMA,LCINDX,?
C   INFTAP : LUIT1,?
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C
      CHARACTER STAR8*8, KEYCNO*6, RTNLBL(2)*8, LAB123(3)*8
      LOGICAL REDOCI, ONLYNO, ONLYFD
      DATA LAB123/'********','********','********'/
      DATA STAR8 /'********'/
C
      CALL QENTER('SIRCAN')
      IF (IPRCAN .GT. 0)
     *   WRITE (LUW4,'(//A/)') ' <<<<< OUTPUT FROM SIRCAN >>>>>'
      REDOCI = FLAG(45)
      ONLYNO = FLAG(46)
      ONLYFD = FLAG(48)
      IF (ONLYNO) THEN
         KEYCNO = 'ONLYNO'
      ELSE IF (ONLYFD) THEN
         KEYCNO = 'ONLYFD'
      ELSE
         KEYCNO = 'FD+NO '
      END IF
C
C ******************************************************************
C Step 0: WRK space allocation:
C
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C .. 1.1 .. (SIRCNO)
      CALL MEMGET('REAL',KCINDX,LCINDX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KAOCC,NASHT,WRK,KFREE,LFREE)
C
C Step 1: Get CI index information (CINDX)
C         Read molecular orbitals (CMO) from disk
C         Read CI-vector (CREF) from disk
C
      IF (NASHT .GT. 1 .AND. .NOT. HSROHF)
     *   CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KFREE),LFREE,0)
      CALL READMO (CMO,9,DUMMY,1)
      CALL RDCREF (CREF)
C
C
C
      KFREE1 = KFREE
      CALL SIRCNO (KEYCNO,1,1,CREF,CMO,WRK(KAOCC),WRK(KCINDX),
     *             WRK,KFREE,LFREE)
C     CALL SIRCNO (KEYCNO,NCVECS,ICREF,CVECS,CMO,AOCC,INDXCI,
C    *             WRK,KFREE,LFREE)
C
      CALL MEMREL('SIRCAN.SIRCNO',WRK,KFRSAV,KFREE1,KFREE,LFREE)
C
C
C
      IF (NCONF .GT. 1 .AND. .NOT.ONLYFD .AND. IPRCAN .GT. 0)
     *   WRITE (LUW4,'(//A/A)')
     *   ' Occupancies of natural active orbitals',
     *   ' --------------------------------------'
      CALL DZERO(OCC,NORBT)
      IF (NASHT .EQ. 1) WRK(KAOCC) = NACTEL
      IAOCC = KAOCC-1
      DO 300 ISYM = 1,NSYM
C
         JOCC1  = IORB(ISYM)
         NISHI  = NISH(ISYM)
         DO 271 J = 1,NISHI
            OCC(JOCC1+J) = D2
  271    CONTINUE
C
         NASHI  = NASH(ISYM)
      IF (NASHI.EQ.0) GO TO 300
         JOCC1  = JOCC1 + NISHI
         OCCSUM = D0
         DO 275 I=1,NASHI
            IAOCC = IAOCC + 1
            OCC(JOCC1 + I) = WRK(IAOCC)
            OCCSUM = OCCSUM + OCC(JOCC1 + I)
  275    CONTINUE
C
         IF (NCONF.GT.1 .AND. .NOT.ONLYFD .AND. IPRCAN .GT. 0) THEN
            WRITE (LUW4,3040) ISYM, OCCSUM, (OCC(JOCC1+I),I=1,NASHI)
         END IF
  300 CONTINUE
C
 3040 FORMAT(/' Symmetry',I3,'  Total occupancy this symmetry',F14.9,
     &      //,(5F14.9))
C
C     Print out natural orbitals
C
      IF (P6FLAG(59)) THEN
         IF (ONLYFD .OR. (NASHT .LE. 1 .AND. .NOT.ONLYNO)) THEN
            WRITE(LUPRI,'(//A/)') ' Canonical orbitals in SIRCAN :'
         ELSE
            WRITE(LUPRI,'(//A/)') ' Natural orbitals in SIRCAN :'
         END IF
         CALL PRORB(CMO,.TRUE.,LUPRI)
      END IF
C
C  Save transformed CI vectors and canonical/natural orbitals
C
C
      IF (IPRCAN .GT. 0) WRITE (LUPRI,'(//A/A)')
     *   ' SIRCAN: canonical/natural orbitals and corresponding',
     *   '         CI coefficents saved on the Sirius restart file.'
      CALL SIRSAV ('CANSAVE',CMO,DUMMY,DUMMY,DUMMY,DUMMY,
     *             WRK(KCINDX),CREF,NCONF)
C     CALL SIRSAV (KEYWRD,CMO,IBNDX,REDL,EVEC,XKAP,INDXCI,WRK,LFREE)
      IF (NASHT .EQ. 1 .OR. (.NOT.ONLYFD .AND. NASHT.GT.1)) THEN
         CALL MOLLAB('NEWORB  ',LUIT1,lupri)
         READ (LUIT1)
         CALL GETDAT(LAB123(2),LAB123(3))
         WRITE (LUIT1) LAB123(1),LAB123(2),'(SIRCAN)','NATOCCAC'
         NAT4 = MAX(NASHT,4)
         CALL WRITT(LUIT1,NAT4,WRK(KAOCC))
         WRITE (LUIT1) LAB123,'EODATA  '
         REWIND LUIT1
      END IF
C
C
C Step 5: Redo CI calculation for natural orbitals.
C         >>> This is a test option <<<
C
      IF (REDOCI .AND. NCONF .GT. 1) THEN
         IF (FLAG(47)) THEN
C            860719-hjaaj : 47 - transform CREF the old way (for test)
            ICIST  = 1
         ELSE
            WRITE (LUPRI,'(//A/)') ' >> calling CICTL for test of TRACI'
            ICIST  = 2
         END IF
C
         IF (.NOT.FLAG(34)) FLAG(14) = .FALSE.
C        if (int.transf. needed in optimization)
         NCROOT = NROOTS
         MAXITC = MAXMAC
         THRCIX = DP5*THRGRD
         CALL MEMGET('REAL',KECI,NCROOT,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KICROO,NCROOT,WRK,KFREE,LFREE)
         CALL MEMGET('WORK',KWRK,LWRK,WRK,KFREE,LFREE)
         ISTACS = ISTACI
         ISTACI = 0
         CALL CICTL(ICIST,NCROOT,MAXITC,THRCIX,CMO,WRK(KCINDX),
     *              WRK(KECI),ICONV,WRK(KICROO),WRK(KWRK),LWRK)
         ISTACI = ISTACS
         CALL MEMREL('SIRCAN.CICTL',WRK,KFRSAV,KWRK,KFREE,LFREE)
C
C        get the CREF for NO's calculated by CICTL
C
         IF (FLAG(47)) THEN
            REWIND LUIT1
            CALL MOLLAB('STARTVEC',LUIT1,lupri)
            DO 980 I = 1,(ISTATE-1)
               READ (LUIT1)
  980       CONTINUE
            CALL READT (LUIT1,NCONF,CREF)
         END IF
      END IF
C
      CALL MEMREL('SIRCAN',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SIRCAN')
      RETURN
C
C end of SIRCAN
C
      END
C  /* Deck sircno */
      SUBROUTINE SIRCNO (KEYCNO,NCVECS,ICREF,CVECS,CMO,AOCC,
     *                   INDXCI,WRK,KFRSAV,LFRSAV)
C
C Written  26-Jul-1986 hjaaj.0
C l.r. 19-Jul-1992 hjaaj (call ORDRSS)
C
C Purpose:
C   SIRCNO calculates Fock canonical and/or natural orbitals
C   and transforms the CI vector accordingly.
C
C MOTECC-90: The purpose of this module, SIRCNO, and the algorithms used
C            are described in Chapter 8 Section E.4 of MOTECC-90
C            "Transformation to Natural and Fock Type Orbitals"
C
C Input:
C   CVECS : NCVECS CI-vectors, no. ICREF determines transformation
C   CMO   : orbitals corresponding to the CI expansion in CVECS
C
C Output:
C   CVECS : CI-vectors in canonical/natural orbitals for vector ICREF
C   CMO   : canonical/natural orbitals
C   AOCC  : occupation of natural active orbitals
C
C
#include <implicit.h>
#include <infvar.h>
      CHARACTER*6 KEYCNO
      DIMENSION CVECS(NCONF,*),  CMO(*), AOCC(*), INDXCI(*), WRK(*)
C
C -- local variables and constants
C
      LOGICAL   DFT_SPINDNS_SAVE
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0 )
      PARAMETER ( THRCAN = 1.0D-10 )
C
C Used from common blocks:
C   GNRINF : WRINDX
C   INFINP : FLAG(*),NACTEL,SUPSYM,LNOROT,?
C   INFVAR : NCONF,?
C   INFORB : NASHT,?
C   INFIND : IROW(*),ISSMO(*),...
C   INFDIM : NBASMA,?
C   INFTAP : 
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <gnrinf.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
#include <dftcom.h>
C
      DIMENSION DVZRAS(3)
      LOGICAL TRTONO, GETOCC, TRTOFC, TRTOCO, TRCVEC, TSTFLG
C
      CALL QENTER('SIRCNO')
      DFT_SPINDNS_SAVE = DFT_SPINDNS
      DFT_SPINDNS = .FALSE.  ! no spin densities needed for canonical orbitals
      IF (IPRCNO .GE. 2) CALL GETTIM (T0,W0)
      IF (IPRCNO .GE. 1) THEN
         WRITE (LUPRI,'(/2A/)')
     &      '   <<< OUTPUT FROM SIRCNO >>>    Keyword = ',KEYCNO
      END IF
      IF (KEYCNO .EQ. 'ONLYNO' .OR. KEYCNO .EQ. 'NATORB') THEN
         IF (IPRCNO .GE. 5) WRITE (LUPRI,'(A/)')
     *      ' Transformation of active orbitals to natural orbitals.'
         TRTONO = .TRUE.
         TRTOFC = .FALSE.
         TRTOCO = .FALSE.
      ELSE IF (KEYCNO .EQ. 'ONLYFD' .OR. KEYCNO .EQ. 'ONLYFC' .OR.
     &         KEYCNO .EQ. 'ONLYCO' .OR. KEYCNO .EQ. 'CANORB') THEN
         IF (IPRCNO .GE. 5) WRITE (LUPRI,'(A/)')
     *      ' Transformation of all orbitals to Fock type orbitals.'
         TRTONO = .FALSE.
         TRTOFC = .TRUE.
         TRTOCO = .TRUE.
      ELSE
         IF (IPRCNO .GE. 5) WRITE (LUPRI,'(A/A/)')
     *      ' Transformation of active orbitals to natural orbitals',
     *      ' and inactive/secondary orbitals to Fock type orbitals.'
         TRTONO = .TRUE.
         TRTOFC = .TRUE.
         TRTOCO = .FALSE.
      END IF
      GETOCC = TRTONO
      TRCVEC = (TRTONO .OR. TRTOCO)
      IF ( TRCVEC ) THEN
         TRCVEC = .FALSE.
         DO 100 ISYM = 1,NSYM
            IF (NASH(ISYM) .GT. 1) TRCVEC = .TRUE.
  100    CONTINUE
         IF (IPRCNO .GE. 5 .AND. .NOT.TRCVEC) WRITE (LUPRI,'(A/A/)')
     *      ' Transformation of active orbitals',
     *      ' abandoned, because nash(isym) .le. 1, all symmetries.'
      END IF
      IF ( TRCVEC .AND. NCONF.EQ.1 ) THEN
         TRCVEC = .FALSE.
         IF ( IPRCNO .GE. 5 ) WRITE (LUPRI,'(A/A/)')
     *      ' Transformation of active orbitals to natural orbitals',
     *      ' abandoned, because only one configuration.'
      END IF
      IF (.NOT. TRCVEC) THEN
         TRTONO = .FALSE.
         TRTOCO = .FALSE.
      END IF
      IF ( TRCVEC .AND. (ICREF .LE. 0 .OR. ICREF .GT. NCVECS) ) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//A,I3,A,I3/A/)')
     *   ' *** WARNING FROM SIRCNO: ICREF illegal, ICREF =',ICREF,
     *   ', NCVECS =',NCVECS,
     *   '     canonical/natural orbital transformation abandoned.'
         GO TO 9999
      END IF
C
C     Exit if LNOROT .true., that is, if orbitals are frozen with
C     .NO ROT (old input) or .FREEZE (new input).
C     SIRCNO will in this version mix the such frozen orbitals with the
C     rest, orbitals frozen with .FROZEN should be handled correctly.
C     891212/901029(handling of NFRO(nsym)/hjaaj
C     901029/hjaaj: also test if frozen orbitals would not be rotated
C     anyway, this makes it possible to use TRACI with, e.g. .CORHOL)
      IF ( LNOROT ) THEN
         TSTFLG = .FALSE.
         DO 110 I = 1,NORBT
            IF (NOROT(I) .EQ. 1) THEN
               ISYMI = ISMO(I)
               ITYPI = IOBTYP(I)
               IF (ITYPI .EQ. JTINAC) THEN
                  IF ((NISH(ISYMI)-NFRO(ISYMI)) .GT. 1) TSTFLG = .TRUE.
               ELSE IF (ITYPI .EQ. JTACT) THEN
                  IACT = ICH(I)
                  IACT = IACTYP(IACT)
                  IF (IACT .EQ. 1) THEN
                     IF (NAS1(ISYMI) .GT. 1) TSTFLG = .TRUE.
                  ELSE IF (IACT .EQ. 2) THEN
                     IF (NAS2(ISYMI) .GT. 1) TSTFLG = .TRUE.
                  ELSE
                     IF (NAS3(ISYMI) .GT. 1) TSTFLG = .TRUE.
                  END IF
               ELSE IF (ITYPI .EQ. JTSEC) THEN
                  IF (NSSH(ISYMI) .GT. 1) TSTFLG = .TRUE.
               END IF
C              (no problem if type is frozen (JTFRO))
            END IF
  110    CONTINUE
         IF (TSTFLG) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//A/A/)')
     *      ' *** WARNING FROM SIRCNO: orbitals frozen with .FREEZE'//
     *      ' are not compatible with',
     *      '     CNO transformation, '//
     *      ' canonical/natural orbital transformation abandoned.'
            GO TO 9999
         END IF
      END IF
C
C ******************************************************************
C Step 0: WRK space allocation:
C
      LWRK   = LFRSAV + KFRSAV - 1
C .. 1.1 .. (MAKDV, TRACI)
      KDV    = KFRSAV
      KWRK1  = KDV    + NNASHX
      LWRK1  = LWRK   - KWRK1
C .. 1.2 .. (FCKMAT)
      KFC    = KWRK1
      KFV    = KFC    + NNORBT
      KWRK2  = KFV    + NNORBT
      IF (SRHYBR) KWRK2 = KWRK2 + NNORBT
C     ... space for extra matrix for SRHYBR
C         but we do not use this matrix in SIRCNO /hjaaj oct 04
      LWRK2  = LWRK   - KWRK2
C .. 1.3  (2*JACO)
      KBIG   = KWRK2
      KJBIG  = KBIG   + MAX(NASHT,NBASMA)
      KCMO2  = KJBIG  + MAX(NASHT,NBASMA)
      KWRK3  = KCMO2  + NCMOT
C Test if WRK is big enough:
      KWRK = MAX(KWRK1,KWRK2,KWRK3)
      IF (KWRK.GT.LWRK) CALL ERRWRK('SIRCNO',-KWRK,LWRK)
C
C     Calculate active density matrix DV (in MAKDV).
C     Average DV if SUPSYM
C
      IF (NASHT .GT. 1 .AND. .NOT. HSROHF) THEN
         CALL MAKDV(CVECS(1,ICREF),WRK(KDV),INDXCI,WRK(KWRK1),LWRK1)
C        CALL MAKDV(CREF,DV,INDXCI,WRK,LFREE)
         IF (SUPSYM) CALL AVEDV(WRK(KDV))
      ELSE IF (HSROHF) THEN
         DO I = 1, NASHT
            II = KDV + I*(I+1)/2 - 1
            WRK(II) = D1
         END DO 
      ELSE IF (NASHT .EQ. 1) THEN
         WRK(KDV) = NACTEL
         IF (GETOCC) AOCC(1) = WRK(KDV)
      END IF
      IF (TRTONO) THEN
         TRTONO = .FALSE.
         DO 40 J = 2, NASHT
            JDV = KDV + IROW(J)
            DO 30 I = 1, (J-1)
               IF (ABS(WRK(JDV+I)) .GT. THRCAN) TRTONO = .TRUE.
   30       CONTINUE
   40    CONTINUE
         IF (IPRCNO .GE. 5 .AND. .NOT.TRTONO)
     *      WRITE (LUPRI,'(A/A,1P,D10.1/)')
     *      ' Transformation of active orbitals to natural orbitals'//
     *      ' abandoned:',' max off-diagonal element in'//
     *      ' density matrix is less than',THRCAN
      END IF
      IF (GETOCC .AND. .NOT.TRTONO) THEN
C     ... orbitals are already NO's
         DO 77 I = 1,NASHT
            AOCC(I) = WRK(KDV-1+IROW(I+1))
   77    CONTINUE
      END IF
C
C
C Step 2: Call FCKMAT to get Fock matrices
C         If active orbitals then add FC and FV.
C
C
      IF ( .NOT. TRTOFC ) GO TO 4000
      CALL FCKMAT((NASHT.EQ.0),WRK(KDV),CMO,EMY,WRK(KFC),WRK(KFV),
     *            WRK(KWRK2),LWRK2)
C     CALL FCKMAT(ONLYFC,DV,CMO,EMY,FC,FV,WRK,LFREE)
      IF (NASHT .GT. 0) CALL DAXPY(NNORBT,D1,WRK(KFV),1,WRK(KFC),1)
      IF (SUPSYM) CALL AVEFCK(WRK(KFC))
C
C     If (TRTOCO) then extract FDAC=FCAC+FVAC for TRACI to canon. orb.
C
      IF (TRTOCO) CALL GETAC(WRK(KFC),WRK(KDV))
C
C Step 3: Diagonalize Fock matrix (FD=FC+FV) to get canonical
C         inactive orbitals and secondary orbitals
C
C
      IF (IPRCNO .GE. 10) WRITE (LUPRI,3010)
 3010 FORMAT(//' Eigenvalues of modified FC + FV Fock matrix',
     *        /' (coupling betw. inactive, active, and secondary',
     *         ' blocks are zeroed)')
      CALL DCOPY(NCMOT,CMO,1,WRK(KCMO2),1)
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         IORBI = IORB(ISYM)
         NBASI = NBAS(ISYM)
         JFC   = KFC   + IIORB(ISYM)
         JCMO2 = KCMO2 + ICMO(ISYM)
C
C        Zero couplings between inactive-active-secondary blocks
C        in Fock matrix and between frozen orbitals and all other.
C
         NFROI = NFRO(ISYM)
         NISHI = NISH(ISYM)
         NOCCI = NOCC(ISYM)
         NASHI = NOCCI - NISHI
         IF (NFROI .GT. 0) THEN
C           zero off-diagonal frozen-frozen elements
            DO I = 2,NFROI
               JSTA = JFC + IROW(I)
               JEND = JSTA + I - 2
               DO J = JSTA,JEND
                  WRK(J) = D0
               END DO
            END DO
         END IF
         IF (NISHI .GT. 0) THEN
C           zero inactive-(active+secondary) block
C           Note that any frozen orbitals are also "inactive".
            CALL DSPZERO(NORBI,WRK(JFC),1,NISHI,NISHI+1,NORBI)
C           call dspzero(n,asp, nrsta,nrend, ncsta,ncend)
         END IF
         IF (NASHI .GT. 0) THEN
C           zero active-secondary block
            CALL DSPZERO(NORBI,WRK(JFC),NISHI+1,NOCCI,NOCCI+1,NORBI)
         END IF
C
C
         CALL JACO(WRK(JFC),WRK(JCMO2),NORBI,NORBI,NBASI,
     *             WRK(KBIG),WRK(KJBIG))
C        CALL JACO (F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
         DO 175 I=1,NORBI
            WRK(KBIG-1 + I) = WRK(JFC-1 + IROW(I+1))
  175    CONTINUE
C
C        skip frozen orbitals
C
         JCMO2 = JCMO2 + NFROI*NBASI
         JCMO  = 1 + ICMO(ISYM) + NFROI*NBASI
         NISHI_NF = NISHI - NFROI
         IF (NISHI_NF .GT. 1) THEN
C           order inactive orbitals which are not frozen
            CALL ORDRSS(WRK(JCMO2),WRK(KBIG+NFROI),
     &                  ISSMO(IORBI+NFROI+1),NISHI_NF,NBASI)
            CALL DCOPY(NISHI_NF*NBASI,WRK(JCMO2),1,CMO(JCMO),1)
         END IF
         JCMO2 = JCMO2 + NISHI_NF*NBASI
         JCMO  = JCMO  + NISHI_NF*NBASI
         IF (.NOT. GETOCC .AND. .NOT. TRTOCO .AND.
     &       NASHI .GT. 1 .AND. .NOT. SUPSYM) THEN
C           order active orbitals
            CALL ORDER(WRK(JCMO2),WRK(KBIG+NISHI),NASHI,NBASI)
            CALL DCOPY(NASHI*NBASI,WRK(JCMO2),1,CMO(JCMO),1)
         END IF
         NSSHI = NORBI - NOCCI
         IF (NSSHI .GT. 1) THEN
            JCMO2 = JCMO2 + NASHI*NBASI
            JCMO  = JCMO  + NASHI*NBASI
C           order secondary orbitals
            CALL ORDRSS(WRK(JCMO2),WRK(KBIG+NOCCI),
     &                  ISSMO(IORBI+NOCCI+1),NSSHI,NBASI)
            CALL DCOPY(NSSHI*NBASI,WRK(JCMO2),1,CMO(JCMO),1)
         END IF
C
         IF (IPRCNO .GE. 10)
     *   WRITE (LUPRI,'(/A,I3,//,(5F14.9))')
     *   ' Symmetry',ISYM, (WRK(KBIG-1 + I),I=1,NORBI)
  200 CONTINUE
      IF (SUPSYM) CALL AVEORD()
C     ... remake ISSORD() as ISSMO() may have changed in ORDRSS
C
C
C Step 4: Diagonalize 1-el density matrix to get active natural
C         orbitals
C
 4000 CONTINUE
      IF (TRTONO .OR. TRTOCO) THEN
C
C     Transform mo coefficients to natural or "canonical" orbitals
C     and counter-transform ci coefficients in TRACI.
C
C     If (TRTONO) then WRK(KDV) contains DV        matrix
C                 else WRK(KDV) contains FCAC+FVAC matrix
C
         IF (IPRCNO .GE. 8) THEN
            IF (TRTONO) THEN
               WRITE (LUPRI,'(/A)') ' DV matrix before TRACI :'
            ELSE
               WRITE (LUPRI,'(/A)')
     &            ' FDAC = FCAC+FVAC matrix before TRACI :'
            END IF
            CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
         END IF
C
C        if RAS, zero blocks between RAS1, RAS2, and RAS3 spaces
C        in RASNO
C
         IF (MCTYPE .EQ. 2) THEN
            CALL RASNO(WRK(KDV),DVZRAS)
            IF (IPRCNO .GE. 8) THEN
               IF (TRTONO) THEN
                  WRITE (LUPRI,'(//A)')
     *            ' RAS modified DV matrix before TRACI :'
               ELSE
                  WRITE (LUPRI,'(//A)')
     *            ' RAS modified FDAC = FCAC+FVAC matrix before TRACI :'
               END IF
               WRITE (LUPRI,'(1P,3(/A,D15.6))')
     &           ' Norm of zeroed RAS1-RAS2 block:',DVZRAS(1),
     &           ' Norm of zeroed RAS1-RAS3 block:',DVZRAS(2),
     &           ' Norm of zeroed RAS2-RAS3 block:',DVZRAS(3)
               CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
            END IF
         END IF
C
C        do NO transformation w/o symmetry blocking (for TRACI)
         KUNO   = KWRK1
         KBIG   = KUNO   + N2ASHX
         KJBIG  = KBIG   + NASHT
         CALL DUNIT(WRK(KUNO),NASHT)
         CALL JACO(WRK(KDV),WRK(KUNO),NASHT,NASHT,NASHT,
     *             WRK(KBIG),WRK(KJBIG))
C
C        select positive phase of transformation vectors
C        (920715-hjaaj: we assume that degenerate orbitals
C         will get equivalent phase for SUPSYM)
C
         JUNO   = KUNO
         DO 410 I = 1,NASHT
            IF (TRTONO) AOCC(I) = WRK(KDV-1 + IROW(I+1))
            JUNOMX = JUNO-1 + IDAMAX(NASHT,WRK(JUNO),1)
            IF (WRK(JUNOMX) .LT. D0) CALL DSCAL(NASHT,DM1,WRK(JUNO),1)
            JUNO   = JUNO + NASHT
  410    CONTINUE
CHJ-S-860827
C        JOCC = KDV
C        JUNO = KUNO
C        DO 420 ISYM = 1,NSYM
C           NASHI = NASH(ISYM)
C        IF (NASHI .EQ. 0) GO TO 420
C           CALL ORDER2(WRK(JUNO),WRK(JOCC),NASHI,NASHT)
C           860827-hj: Do not reorder active orbitals;
C                      1) reduce may have been specified
C                      2) we may get small diagonal elements in
C                         UNO(NASHT,NASHT) below; TRACI only works
C                         for non-zero diagonal elements.
C           JOCC = JOCC + NASHI
C           JUNO = JUNO + NASHI*NASHT
C 420    CONTINUE
CHJ-E
         IF (TRTONO .AND. IPRCNO .GE. 1) THEN
            IF (MCTYPE .EQ. 2) THEN
               WRITE (LUPRI,'(/A)')
     &            ' "Occupations" of RAS pseudo-natural orbitals:'
               WRITE (LUPRI,'(1P,3(/A,D15.6))')
     &           ' Norm of zeroed RAS1-RAS2 block:',DVZRAS(1),
     &           ' Norm of zeroed RAS1-RAS3 block:',DVZRAS(2),
     &           ' Norm of zeroed RAS2-RAS3 block:',DVZRAS(3)
            ELSE
               WRITE (LUPRI,'(/A)')
     &            ' Occupations of CAS natural orbitals:'
            END IF
            DO 440 ISYM = 1,NSYM
               IF (NASH(ISYM) .GT. 0) THEN
                  WRITE (LUPRI,'(/A,I3/)') ' Symmetry',ISYM
                  IST  = IASH(ISYM) + 1
                  IEND = IASH(ISYM) + NASH(ISYM)
                  WRITE (LUPRI,'(5F14.9)') (AOCC(I), I = IST,IEND)
               END IF
  440       CONTINUE
         END IF
         IF (IPRCNO .GE. 7) THEN
            IF (TRTONO) THEN
               WRITE (LUPRI,'(/A)')
     *         ' Transformation to NO (active orb.)'
            ELSE
               WRITE (LUPRI,'(/A)')
     *         ' Transformation to Fock-type orbitals (active orb.)'
            END IF
            CALL OUTPUT(WRK(KUNO),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
C
C900801-hjaaj: If transformation is unit matrix, skip TRACI.
C
         CALL DUNIT(WRK(KBIG),NASHT)
         CALL DAXPY(N2ASHX,DM1,WRK(KUNO),1,WRK(KBIG),1)
         I = IDAMAX(N2ASHX,WRK(KBIG),1)
         IF (ABS(WRK(KBIG-1+I)) .LE. THRCAN) THEN
            IF (IPRCNO .GE. 5) WRITE (LUPRI,'(/A/A,1P,D10.1/)')
     *         ' Transformation of active orbitals'//
     *         ' and CI vector abandoned',
     *         ' because max deviation from unit transformation'//
     *         ' is ',ABS(WRK(KBIG-1+I))
            GO TO 9999
         END IF
C
C
         DO 460 ISYM = 1,NSYM
            NASHI  = NASH(ISYM)
         IF (NASHI.EQ.0) GO TO 460
            NBASI  = NBAS(ISYM)
            ICMOA  = 1    + ICMO(ISYM) + NISH(ISYM)*NBASI
            JUNO   = KUNO + IASH(ISYM)*NASHT + IASH(ISYM)
C
            CALL DCOPY(NASHI*NBASI,CMO(ICMOA),1,WRK(KBIG),1)
            CALL DGEMM('N','N',NBASI,NASHI,NASHI,1.D0,
     &                 WRK(KBIG),NBASI,
     &                 WRK(JUNO),NASHT,0.D0,
     &                 CMO(ICMOA),NBASI)
C
  460    CONTINUE
C
C
C Step 5: Redo CI calculation for natural orbitals.
C
         IF (.NOT. FLAG(47)) THEN
C         FLAG(47)860719-hjaaj: transform CVECS the old way (for test)
C
            CALL TRACI(NCVECS,CVECS,NCONF,WRK(KUNO),INDXCI,
     *                 WRK(KBIG),(LWRK-KBIG),IPRCNO)
C           CALL TRACI(NCVEC,CVEC,LCVEC,UMO,INDXCI,WRK,LFREE,IPRCNO)
            IF (IPRCNO .GT. 8 .AND. TRTONO) THEN
               WRITE (LUPRI,'(//A,I4/)')
     *         '  >> SIRCNO TRACI test, calling MAKDV'//
     *         ' with CI vector = no.',ICREF
               CALL MAKDV(CVECS(1,ICREF),WRK(KDV),
     *                    INDXCI,WRK(KWRK1),LWRK1)
               WRITE (LUPRI,'(//A)') ' DV matrix after TRACI :'
               CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
            END IF
         END IF
C
      END IF
C     ... end if (trtono .or. trtofc)
C
 9999 CONTINUE
      IF (SUPSYM) THEN
         KFREE = 1
         LFREE = LWRK
         WRINDX = .TRUE.
         CALL SIRSET(WRK(KFRSAV),LWRK,.FALSE.)
C        ... redo JWOP in case ISSMO should have changed in ORDRSS /hjaaj aug 04
         CALL AVECPH(IPHCHA,CMO,WRK(KFRSAV),KFREE,LFREE)
         IF (IPHCHA .GT. 0) THEN
            WRITE(LUPRI,'(/A,I4,A/A)')
     &      ' SIRCNO ERROR:',IPHCHA,' degenerate m.o.(s) have',
     &      '    no longer the same relative phase'
            CALL QUIT('SIRCNO: sup.sym. MO phase error')
         END IF
      END IF
      IF (IPRCNO .GE. 2) THEN
         CALL GETTIM (TEND,WEND)
         WRITE (LUPRI,'(/A,F9.2)') ' CPU time used in SIRCNO :',TEND-T0
      END IF
      DFT_SPINDNS = DFT_SPINDNS_SAVE
      CALL QEXIT('SIRCNO')
      RETURN
C
C     End of SIRCNO.
C
      END
C  /* Deck rasno */
      SUBROUTINE RASNO(DV,DVZRAS)
C
C 30-Jan-1989 hjaaj
C
C Zero blocks in one-electron density matrix DV
C connecting different RAS spaces (RAS1, RAS2, RAS3).
C
C
#include <implicit.h>
      DIMENSION DV(NNASHX),DVZRAS(3)
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C   INFORB : NASHT, NNASHX
C   INFIND : IROW(*),IACTYP(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      DIMENSION DVZTMP(3,3)
C
      CALL DZERO(DVZTMP,9)
      DO 200 J = 1,NASHT
         IACTJ = IACTYP(J)
         JR = IROW(J)
         DO 100 I = 1,J-1
            IACTI = IACTYP(I)
            IF (IACTI .NE. IACTJ) THEN
               DVZTMP(IACTI,IACTJ) = DVZTMP(IACTI,IACTJ)
     &            + DV(JR + I) ** 2
               DV(JR + I) = D0
            END IF
  100    CONTINUE
  200 CONTINUE
      DVZRAS(1) = SQRT(DVZTMP(2,1) + DVZTMP(1,2))
      DVZRAS(2) = SQRT(DVZTMP(3,1) + DVZTMP(1,3))
      DVZRAS(3) = SQRT(DVZTMP(3,2) + DVZTMP(2,3))
      RETURN
      END
C  /* Deck getno */
      SUBROUTINE GETNO(CVEC,ICSYM,OCCNO,UNO,CMO,ORDER,NATORB,NOAVER,
     *                 MOISNO,INDXCI,WRK,KFRSAV,LFRSAV,LUPRI,IPRINT)
C
C 29-Apr-1989 Hans Joergen Aa. Jensen
C
C Purpose:
C  Calculate natural orbital occupation numbers and, if NATORB is true,
C  transform orbitals to natural orbitals.
C
#include <implicit.h>
      DIMENSION CVEC(*), OCCNO(*), UNO(NASHT,*), CMO(*), INDXCI(*),
     *          WRK(*)
      LOGICAL   ORDER, NATORB, MOISNO, NOAVER
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, D2 = 2.0D0)
#include <thrzer.h>
C
C Used from common blocks:
C  GNRINF: WRINDX
C  INFINP: LSYM
C  INFORB: ICMO(*),...
C  INFIND: ISSMO()
C
#include <maxash.h>
#include <maxorb.h>
#include <gnrinf.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('GETNO ')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      IF (IPRINT .GE. 2) THEN
         WRITE (LUPRI,'(//A//A,I2/)')
     *      ' <<< output from getno >>>',
     *      ' Symmetry of state vector :',ICSYM
      END IF
C
C     Calculate natural orbital transformation and
C     natural orbital occupation numbers
C
      CALL MEMGET('REAL',KDV  ,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KBIG ,NASHT ,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KJBIG,NASHT ,WRK,KFREE,LFREE)
      CALL MEMGET('WORK',KWRK1,LWRK1 ,WRK,KFREE,LFREE)
      IF (NASHT .EQ. 1) THEN
         WRK(KDV) = D1
      ELSE IF (HSROHF) THEN
         CALL DZERO(WRK(KDV),NNASHX)
         DO I = 1, NASHT
            II = KDV + I*(I+1)/2 - 1
            WRK(II) = D1
         END DO 
      ELSE IF (NASHT .GT. 1) THEN
         LSYMSV = LSYM
         LSYM   = ICSYM
         CALL MAKDV(CVEC,WRK(KDV),INDXCI,WRK(KWRK1),LWRK1)
         LSYM   = LSYMSV
      END IF
C
      IF (NASHT .GT. 0) THEN
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI,'(//A)') ' GETNO: 1-electron density matrix'
            CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
         END IF
C
         IF ( SUPSYM .AND. .NOT.NOAVER ) THEN
            CALL AVEDV(WRK(KDV))
            IF (IPRINT .GE. 10) THEN
               WRITE(LUPRI,'(//A/A)')
     *            ' Density matrix adapted to supersymmetry',
     *            ' ======================================='
               CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
            END IF
         END IF
         DVMAX = D0
         DO 220 I = 2,NASHT
            II = (I*I-I)/2
            DO 210 J = 1,I-1
               DVMAX = MAX(DVMAX,ABS(WRK(KDV-1+II+J)))
  210       CONTINUE
  220    CONTINUE
         IF (DVMAX .LE. THRZER) THEN
            MOISNO = .TRUE.
            IF (IPRINT .GE. 2) THEN
               WRITE(LUPRI,'(/A)')
     *           ' GETNO: The orbitals are already natural orbitals.'
            ENDIF
         ELSE
            MOISNO = .FALSE.
         END IF
         CALL DUNIT(UNO,NASHT)
      ELSE
         MOISNO = .TRUE.
      END IF
C
      IF (.NOT. MOISNO) THEN
         CALL JACO(WRK(KDV),UNO,NASHT,NASHT,NASHT,WRK(KBIG),WRK(KJBIG))
         DO 240 I = 1,NASHT
            JUNOMX = IDAMAX(NASHT,UNO(1,I),1)
            IF (UNO(JUNOMX,I) .LT. D0) CALL DSCAL(NASHT,DM1,UNO(1,I),1)
  240    CONTINUE
      END IF
C
C     extract natural occupation numbers
C     select positive phase of transformation vectors
C
      JOCC = 0
      JACT = 0
      DO 360 ISYM = 1,NSYM
         DO 320 I = 1,NISH(ISYM)
            JOCC = JOCC + 1
            OCCNO(JOCC) = D2
  320    CONTINUE
         DO 340 I = 1,NASH(ISYM)
            JOCC = JOCC + 1
            JACT = JACT + 1
            OCCNO(JOCC) = WRK(KDV-1 + (JACT*JACT+JACT)/2)
  340    CONTINUE
  360 CONTINUE
      IF (ORDER) THEN
C     ... reasons why reordering of active orbitals may not be desired:
C         1) reduce may have been specified and TRACI will be called.
C         2) we may get small diagonal elements in
C            UNO(NASHT,NASHT) below; TRACI only works
C            for non-zero diagonal elements.
         DO 420 ISYM = 1,NSYM
            NASHI = NASH(ISYM)
         IF (NASHI .EQ. 0) GO TO 420
            JACT = IASH(ISYM) + 1
            JORB = IORB(ISYM) + 1
            JOCCA= IOCC(ISYM) + NISH(ISYM) + 1
            CALL ORD2SS(UNO(1,JACT),OCCNO(JOCCA),
     &                  ISSMO(JORB),NASHI,NASHT)
  420    CONTINUE
         IF (SUPSYM) CALL AVEORD()
C        ... regenerate ISSORD() as ISSMO() may have been changed in ORD2SS
         IF (MOISNO) THEN
            DO 430 I = 1,NASHT
  430          IF (UNO(I,I) .NE. D1) MOISNO = .FALSE.
            IF (IPRINT .GE. 2 .AND. .NOT. MOISNO) THEN
               WRITE(LUPRI,'(/A)')
     *           ' GETNO: The natural orbitals are reordered.'
            ENDIF
        END IF
      END IF
      IF (IPRINT .GE. 1) THEN
         WRITE (LUPRI,'(//A/A)')
     *      ' Occupancies of natural orbitals',
     *      ' -------------------------------'
         IF (SUPSYM .AND. .NOT. NOAVER) WRITE (LUPRI,'(/A)')
     *      ' NOTE: degenerate orbitals have been averaged.'
         DO 520 ISYM = 1,NSYM
            NOCCI = NOCC(ISYM)
            IF (NOCCI .EQ. 0) THEN
               WRITE (LUPRI,3030) ISYM
            ELSE
               OCCSUM = DSUM(NOCCI,OCCNO(IOCC(ISYM)+1),1)
               IOCCI = IOCC(ISYM)
               IORBI = IORB(ISYM)
               IF (SUPSYM) THEN
                  WRITE (LUPRI,3045) ISYM,OCCSUM,
     &               (OCCNO(IOCCI+I),ISSMO(IORBI+I),I=1,NOCCI)
               ELSE
                  WRITE (LUPRI,3040) ISYM,OCCSUM,
     &               (OCCNO(IOCCI+I),I=1,NOCCI)
               END IF
            END IF
  520    CONTINUE
C  HJAAJ MAERKE TODO 921106: if (supsym) print also sorted after supsym
C  HJAAJ MAERKE TODO 071108: if (supsym) also OCCSUM for each supersym
      END IF
 3030 FORMAT(/' Symmetry',I3,' -- No occupied orbitals')
 3040 FORMAT(/' Symmetry',I3,' -- Total occupation in this symmetry is',
     &          F14.9,//,(5F14.9))
 3045 FORMAT(/' Symmetry',I3,' (supersymmetry in parenthesis)',
     &       /' -- Total occupation in this symmetry is',F14.9,
     *       //,4(F14.9,' (',I2,')'))
      IF (IPRINT .GT. 5 .AND. NASHT .GT. 0 .AND. .NOT.MOISNO) THEN
         DO 620 ISYM = 1,NSYM
         IF (NASH(ISYM) .GT. 0) THEN
            WRITE (LUPRI,'(/A,I3)')
     *         ' Transformation to NO (active orb.) for symmetry',ISYM
            IST  = IASH(ISYM) + 1
            IEND = IASH(ISYM) + NASH(ISYM)
            CALL OUTPUT(UNO,IST,IEND,IST,IEND,NASHT,NASHT,1,LUPRI)
         END IF
  620    CONTINUE
      END IF
C
C     If NATORB, transform CMO to natural orbitals.
C
      IF (.NOT.MOISNO) THEN
      IF (.NOT.NATORB) THEN
         IF (IPRINT .GE. 5) WRITE (LUPRI,'(/A)')
     *      ' GETNO: the orbital coefficients are NOT transformed to'//
     *      ' natural orbitals.'
      ELSE
         IF (IPRINT .GE. 5) WRITE (LUPRI,'(/A)')
     *      ' GETNO: the orbital coefficients are transformed to'//
     *      ' natural orbitals.'
         IF (LWRK1 .LT. NCMOT) CALL ERRWRK('GETNO.NATORB',NCMOT,LWRK1)
         DO 800 ISYM = 1,NSYM
            NASHI  = NASH(ISYM)
         IF (NASHI.EQ.0) GO TO 800
            NBASI  = NBAS(ISYM)
            ICMOA  = 1 + ICMO(ISYM) + NISH(ISYM)*NBASI
            JUNO   = 1 + IASH(ISYM)
C
            CALL DCOPY(NASHI*NBASI,CMO(ICMOA),1,WRK(KWRK1),1)
            CALL DGEMM('N','N',NBASI,NASHI,NASHI,1.D0,
     &                 WRK(KWRK1),NBASI,
     &                 UNO(JUNO,JUNO),NASHT,0.D0,
     &                 CMO(ICMOA),NBASI)
C
  800    CONTINUE
         MOISNO = .TRUE.
      END IF
      END IF
C
C
      CALL MEMREL('GETNO',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      IF (ORDER .AND. SUPSYM) THEN
         WRINDX = .TRUE.
         CALL SIRSET(WRK(KFREE),LFREE,.FALSE.)
C        ... redo JWOP in case ISSMO should have changed in ORDRSS /hjaaj aug 04
      END IF
      CALL QEXIT('GETNO ')
      RETURN
      END
