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
C File:    sirgp.F
C Purpose: general purpose routines for the Sirius module
C
C  /* Deck errwrk */
      SUBROUTINE ERRWRK (STRING,LNEED,LAVAIL)
C
C Version 6-Mar-1985 by hjaaj
C
      CHARACTER*(*) STRING
C
#include <priunit.h>
C
      CALL QENTER('ERRWRK')
      IF (LNEED .GE. 0) THEN
         WRITE (LUPRI,1010) STRING,LNEED,LAVAIL
      ELSE
         WRITE (LUPRI,1020) STRING,-LNEED,LAVAIL
      END IF
      CALL QTRACE(LUPRI)
      CALL QUIT('ERROR, insufficient work space in memory')
C
 1010 FORMAT(/'  FATAL ERROR, insufficient core for ',A,
     *       /T16,'( Need:',I10,', available:',I10,' )')
 1020 FORMAT(/'  FATAL ERROR, insufficient core for ',A,
     *      //T16,'Need      :',I10,' + uncalculated amount',
     *       /T16,'Available :',I10)
      END
C  /* Deck symort */
      SUBROUTINE SYMORT(CMO,SAO,NAO,NMO,WRK,LFREE,IRETUR)
C
C 28-Apr-1985 Hans Jorgen Aa. Jensen
C Revised  9-Aug-1989 hjaaj (exit if diverging)
C          1-Dec-1995 hjaaj (new exit for numerical problems,
C                            see comments)
C
C Purpose:
C  Symmetrical orthonormalization of orbitals.
C
C Input:
C  CMO(*), non-orthogonal orbitals
C  SAO(*), ao overlap matrix
C
C Output:
C  CMO(*), normalized and symmetrically orthogonalized orbitals.
C
C Scratch:
C  WRK(LFREE)
C
#include <implicit.h>
      DIMENSION CMO(NAO,*), SAO(*), WRK(*)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D3 = 3.0D0, DMP5 = -0.5D0 )
      PARAMETER ( ITSMAX = 20, TMIN = 1.D-6, CNVFAC = 0.1D0)
      PARAMETER ( DIFMAX=1.0D-12, DIFMIN=1.0D-4)
C
#include <priunit.h>
#include <infpri.h>
C
      CALL QENTER('SYMORT')
C
C Symmetric Orthogonalization with N-R iterative algorithm:
C
C  C(p+1) = 1/2( 3C(p) - C(p) [Ct(p) Sao C(p)])
C         = 3/2 C(p) - 1/2 C(p) Smo(p)
C         = C(p) [-1/2] [Smo(p) - 3 I]
C
      NCMOI = NAO*NMO
      NNMO  = NMO*(NMO+1)/2
      KSMO  = 1
      KW    = KSMO + NNMO
      LNEED = KW + NAO*NMO - 1
      IF (LNEED .GT. LFREE) CALL ERRWRK('SYMORT',LNEED,LFREE)
C
C     Normalize initial vectors
C     (norm must be less than sqrt(3) for proper convergence,
C      if norm .gt. sqrt(5) the algorithm will diverge)
C     951201: this loop is moved inside iteration loop.
C
C     951201/hjaaj:
C     Numerical round-off limit may be reached long before
C     DELMAX .lt. DIFMAX, for example,
C     DIFMAX used to be 1.0D-8, but that caused problems
C     (divergence from e.g. 1.4D-8 to 2.9D-8) for basis sets
C     with eigenvalues of overlap matrix down around 1.0D-6).
C     New exit without error print if DELMAX decreases
C     too slowly when DELMAX .lt. DIFMIN (=1.0D-4).
C     NOTE: because of these round-off problems you should
C     always follow SYMORT with Gram-Schmidt orthonormalizations
C     to come as close as possible to full orthonormalization.
C     By the way, test cases show that iterations (may) converge
C     even though DELMAX increases in initial iterations, but this
C     implementation will exit if DELMAX increases.
C     Test program also indicates that cubic convergence is
C     obtained if CMO vectors are renormalized in each iteration,
C     therefore 10-loop is moved inside iteration loop.
C
C     Symmetric orthonormalization, using Newton-Raphson iterations
C
      PRVMAX = D0
C     ... to avoid compiler messages
      ITS = 0
  100 CONTINUE
         ITS = ITS + 1
         DO 10 I = 1,NMO
            CALL MPAPV(NAO,SAO,CMO(1,I),WRK)
            T = DDOT(NAO,WRK,1,CMO(1,I),1)
            IF (T .LT. TMIN) THEN
               IMO = I
               GO TO 8820
            END IF
            T = D1 / SQRT(T)
            CALL DSCAL(NAO,T,CMO(1,I),1)
   10    CONTINUE
C
         CALL UTHU(SAO,WRK(KSMO),CMO,WRK(KW),NAO,NMO)
C        CALL UTHU(H,HT,U,WRK,NAO,NMO)
         II = KSMO - 1
         DO 200 I = 1,NMO
            II = II + I
            WRK(II) = WRK(II) - D3
  200    CONTINUE
         CALL DSCAL(NNMO,DMP5,WRK(KSMO),1)
C
C
         CALL DZERO(WRK(KW),NCMOI)
         IJ  = KSMO - 1
         IC1 = KW
         DO 400 I = 1,NMO
            JC1 = KW
            DO 300 J = 1,I-1
               IJ = IJ + 1
               FAC = WRK(IJ)
               CALL DAXPY(NAO,FAC,CMO(1,I),1,WRK(JC1),1)
               CALL DAXPY(NAO,FAC,CMO(1,J),1,WRK(IC1),1)
               JC1 = JC1 + NAO
  300       CONTINUE
            IJ = IJ + 1
            FAC = WRK(IJ)
            CALL DAXPY(NAO,FAC,CMO(1,I),1,WRK(IC1),1)
            IC1 = IC1 + NAO
  400    CONTINUE
C
         DELMAX = D0
         IJ = KW - 1
         DO 550 J = 1,NMO
            DO 540 I = 1,NAO
               IJ = IJ + 1
               DEL  = ABS(CMO(I,J)-WRK(IJ))
               DELMAX = MAX(DELMAX,DEL)
  540       CONTINUE
  550    CONTINUE
         CALL DCOPY (NCMOI,WRK(KW),1,CMO,1)
      IF (DELMAX.GT.DIFMAX) THEN
         IF (ITS .GE. ITSMAX) GO TO 8810
         IF (ITS .GT. 1) THEN
            IF (DELMAX .LT. DIFMIN) THEN
               IF (DELMAX .GT. CNVFAC*PRVMAX) GO TO 8840
C              ... convergence cannot be quadratic in this case,
C                  exit because we probably have numeric lin.dep.
            END IF
            IF (DELMAX.GT.PRVMAX) GO TO 8830
         END IF
         PRVMAX = DELMAX
         GO TO 100
      END IF
C
      IF (P6FLAG(12)) WRITE (LUPRI,2211) ITS,DELMAX
 2211 FORMAT(/' (SYMORT) Symmetric orthogonalization in',
     *        I3,' iterations.',/T11,'Max. error: ',F15.10)
C
      IRETUR = 0
      GO TO 9999
C
 8810 IRETUR = ITS
      WRITE (LUPRI,2211) ITS,DELMAX
      GO TO 9999
C
 8820 IRETUR = -IMO
      WRITE (LUPRI,2311) IMO, T
 2311 FORMAT(/' (SYMORT) ***ERROR*** norm of input vector no.',I4,
     *        ' is',1P,D10.2)
      GO TO 9999
C
 8830 WRITE (LUPRI,2411) ITS,DELMAX,PRVMAX
      IRETUR = 9999
 2411 FORMAT(/' (SYMORT) Symmetric orthogonalization diverging',
     *       /' Iteration no.',I5,'  Max. error :',1P,T40,D16.8,
     *       /' Max. error previous iteration :',T40,D16.8)
      GO TO 9999
C
 8840 IF (P6FLAG(12)) WRITE (LUPRI,2511) ITS,DELMAX,PRVMAX
      IRETUR = 8888
 2511 FORMAT(/' (SYMORT) Symmetric orthogonalization diverging',
     &       /' because of num. lin.dep. causing round-off errors',
     *       /' Iteration no.',I5,'  Max. error :',1P,T40,D16.8,
     *       /' Max. error previous iteration :',T40,D16.8)
C
C End of "SYMORT"
C
 9999 CALL QEXIT('SYMORT')
      RETURN
      END
C  /* Deck vsymtr */
      SUBROUTINE VSYMTR(LEN,VEC,THREQL)
C
C  5-Jan-1989  Hans Joergen Aa. Jensen
C  Last revision 26-Nov-1992 hjaaj: do not compare zero elements
C
C  Purpose: symmetrize vector so that numerical
C           inaccuracy won't break symmetry (e.g. in a CI vector)
C
C  Quite slow because n ** 2 process
C
#include <implicit.h>
      DIMENSION VEC(LEN)
      PARAMETER (D0 = 0.0D0)
C
      CALL QENTER('VSYMTR')
      IF (THREQL .LE. D0) GO TO 9999
      DO 200 ICONF = 1,LEN-1
         CTEST = VEC(ICONF)
         IF (ABS(CTEST) .LE. THREQL) THEN
            VEC(ICONF) = D0
            GO TO 200
         END IF
         DO 100 JCONF = ICONF+1, LEN
         IF      (ABS(CTEST-VEC(JCONF)) .LE. THREQL) THEN
            VEC(JCONF) = CTEST
         ELSE IF (ABS(CTEST+VEC(JCONF)) .LE. THREQL) THEN
            VEC(JCONF) = -CTEST
         END IF
  100    CONTINUE
  200 CONTINUE
C
 9999 CALL QEXIT('VSYMTR')
      RETURN
      END
C  /* Deck upkwop */
      SUBROUTINE UPKWOP(NWOPT,JWOP,WOP,UPWOP)
C
C Written 1-Nov-1983 by Hans Jorgen Aa. Jensen in Uppsala
C Revisions
C  31-Jan-1984 hjaaj
C  31-Oct-1989 hjaaj : N2ORBX instead of N2ORBT (thus symemtry
C     independent); zero UPWOP
C
C This subroutine unpacks ("UPK") an orbital rotation type vector
C ("WOP") to matrix form.
C This could for instance be the orbital gradient, or
C the kappa vector for use in one-index transformation.
C
#include <implicit.h>
      DIMENSION UPWOP(NORBT,NORBT),WOP(NWOPT)
      INTEGER JWOP(2,NWOPT)
C
C Used from common blocks:
C   INFORB : NORBT,N2ORBX
C
#include <inforb.h>
C
      CALL DZERO(UPWOP,N2ORBX)
      DO 200 IG = 1,NWOPT
         I = JWOP(1,IG)
         J = JWOP(2,IG)
         UPWOP(I,J) =  WOP(IG)
         UPWOP(J,I) = -WOP(IG)
  200 CONTINUE
C
      RETURN
      END
C  /* Deck pkwop */
      SUBROUTINE PKWOP(NWOPT,JWOP,WOP,UPWOP)
C
C Written 30-Jul-1986 by Hans Jorgen Aa. Jensen
C Revised 31-Oct-1989 hjaaj (removed symmetry blocking of UPWOP)
C
C This subroutine packs ("PK") an orbital rotation type vector
C ("WOP") from matrix form (the reverse of UPKWOP)
C
#include <implicit.h>
      DIMENSION UPWOP(NORBT,NORBT),WOP(NWOPT)
      INTEGER JWOP(2,NWOPT)
C
C Used from common blocks:
C   INFORB : NORBT
C
#include <inforb.h>
C
      DO 200 IG = 1,NWOPT
         I = JWOP(1,IG)
         J = JWOP(2,IG)
         WOP(IG) = UPWOP(I,J)
  200 CONTINUE
C
      RETURN
      END
C  /* Deck getwop */
      LOGICAL FUNCTION GETWOP(XLABEL)
C
C Purpose: read new orbital rotation specification into INFVAR.
C
C 9-Nov-1985 hjaaj
C 1-May-2000 hjaaj: Removed reading of KLWOP (KLWOP not in infvar.h any more)
C
#include <implicit.h>
      CHARACTER*8 XLABEL
C
C Used from common blocks:
C   INFORB: MULD2H()
C   INFVAR: MAXOCC,NWOPT,NVAR,NWOPH,NVARH,JWOP(2,*),JWOPSY
C   INFLIN: LSYMRF,LSYMPT,LSYMST,NCONST,NWOPPT,NVARPT
C   INFDIM: NWOPDI,NWOPMA,NCONMA,NVARMA
C   INFPRI: NWARN
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infvar.h>
#include <inflin.h>
#include <infdim.h>
#include <infpri.h>
C
      LOGICAL FNDLAB
C
      LUNIT = -1
      CALL GPOPEN(LUNIT,'LUINDF','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUNIT
      IF ( FNDLAB(XLABEL,LUNIT) ) THEN
         GETWOP = .TRUE.
         READ (LUNIT) NWOPT, NWOPH, JWOPSY
         CALL READI  (LUNIT,(2*NWOPH),JWOP)
         NWOPDI = MAX(1,NWOPT)
         NWOPMA = MAX(NWOPMA,NWOPT,NWOPH)
         NVAR   = NCONF + NWOPT
         NVARH  = NCONF + NWOPH
         NVARMA = NCONMA+ NWOPMA
         IF (LSYMST .NE. MULD2H(LSYMRF,JWOPSY)) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//A/A/)')
     &      ' GETWOP WARNING: NVARPT not defined because NCONST',
     &      '                 corresponds to wrong symmetry.'
            WRITE (LUPRI,*) 'JWOPSY,LSYMPT',JWOPSY,LSYMPT
            WRITE (LUPRI,*) 'LSYMRF,LSYMPT,LSYMST',LSYMRF,LSYMPT,LSYMST
            WRITE (LUPRI,*) 'NCONRF,NCONST,NCONF ',NCONRF,NCONST,NCONF
            WRITE (LUPRI,*) 'NWOPPT,NWOPT ,NVARPT',NWOPPT,NWOPT ,NVARPT
            CALL QTRACE(LUPRI)
         ELSE
            LSYMPT = JWOPSY
            NWOPPT = NWOPT
            NVARPT = NCONST + NWOPPT
         END IF
      ELSE
         GETWOP = .FALSE.
      END IF
      CALL GPCLOSE(LUNIT,'KEEP')
C
      RETURN
      END
C  /* Deck make_klwop */
      SUBROUTINE MAKE_KLWOP(KLWOP)
C
C 1-May-2000 hjaaj
C
C Purpose: Generate KLWOP array from JWOP
C
#include <implicit.h>
      DIMENSION   KLWOP(NOCCT,NORBT)
      CHARACTER*8 XLABEL
C
C Used from common blocks:
C   INFORB: NOCCT,NORBT
C   INFIND: ISW(*)
C   INFVAR: NWOPT,JWOP(2,*)
C
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
C
      NKLWOP = NOCCT*NORBT
      CALL IZERO(KLWOP,NKLWOP)
C
      DO IG = 1,NWOPT
         KW = ISW(JWOP(1,IG))
         LW = ISW(JWOP(2,IG))-NISHT
         IF (KLWOP(KW,LW) .EQ. 0) THEN
            KLWOP(KW,LW)  = IG
         ELSE
            WRITE (LUPRI,'(/A,I10,A,I10//A)')
     &      'Fatal error, duplicate entry in JWOP found for elements',
     &      KLWOP(KW,LW),' and',IG,'Dump of JWOP:'
            WRITE (LUPRI,'(3I10)') (I,JWOP(1,I),JWOP(2,I),I=1,NWOPT)
            CALL QUIT
     &      ('MAKE_KLWOP: Fatal error, duplicate entry in JWOP.')
         END IF
      END DO
C
      RETURN
      END
C  /* Deck prkap */
      SUBROUTINE PRKAP (NKAPT,XKAP,PRFAC,LUPRIN)
C
C Jan 90 hjaaj (descendant of PRWOP)
C
C Purpose:
C  Write orbital operator array XKAP to unit LUPRIN.
C
C  if PRFAC  .gt. 0, write only elements with absolute value
C                    .ge. PRFAC * (max. absolute value)
C  if PRFAC  .eq. 0, write all elements;
C  if PRFAC  .lt. 0, write only elements with
C                    value .le. PRFACX * min. pos value
C                    where PRFACX = MAX(-PRFAC, -1.0/PRFAC)
C
#include <implicit.h>
      DIMENSION XKAP(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFVAR : NWOPT,NWOPH,JWOPSY,JWOP(2,*)
C
#include <maxorb.h>
#include <infvar.h>
C
C Check if NKAPT legal
C
      IF (NKAPT .NE. NWOPT .AND. NKAPT .NE. NWOPH) THEN
         WRITE (LUPRIN,'(//A/A,3I10)')
     &      ' ERROR in PRKAP, NKAPT .ne. NWOPT,NWOPH',
     &      ' NKAPT, NWOPT, NWOPH :',NKAPT,NWOPT,NWOPH
         CALL QUIT(' ERROR in PRKAP, NKAPT .ne. NWOPT,NWOPH')
      END IF
C
C Any elements?
C
      IF (NKAPT.LE.0) THEN
         WRITE (LUPRIN,8000) JWOPSY
         GO TO 9999
      END IF
 8000 FORMAT(/5X,'Orbital operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C Print XKAP as requested:
C
      IF (PRFAC .GT. D0) THEN
         KMAX   = IDAMAX(NKAPT,XKAP,1)
         THPKAP = PRFAC*ABS(XKAP(KMAX))
         WRITE (LUPRIN,1010) JWOPSY,PRFAC
         INPR = 0
         DO 100 IG = 1,NKAPT
            IF ( ABS(XKAP(IG)) .GT. THPKAP ) THEN
               K = JWOP(1,IG)
               L = JWOP(2,IG)
               WRITE (LUPRIN,1100) IG,K,L,XKAP(IG)
            ELSE
               INPR = INPR + 1
            END IF
  100    CONTINUE
         IF (INPR.GT.0) WRITE (LUPRIN,1210) INPR,THPKAP
      ELSE IF (PRFAC .EQ. D0) THEN
         WRITE (LUPRIN,1020) JWOPSY
         INPR = 0
         DO 200 IG = 1,NKAPT
            IF ( XKAP(IG) .NE. D0 ) THEN
               K = JWOP(1,IG)
               L = JWOP(2,IG)
               WRITE (LUPRIN,1100) IG,K,L,XKAP(IG)
            ELSE
               INPR = INPR + 1
            END IF
  200    CONTINUE
         IF (INPR.GT.0) WRITE (LUPRIN,1220) INPR
      ELSE
C     (PRFAC .LT. D0)
         KMIN   = IDAMIN(NKAPT,XKAP,1)
         PRFACX = MAX( -PRFAC, -D1/PRFAC)
         THPKAP = PRFACX*ABS(XKAP(KMIN))
         WRITE (LUPRIN,1030) JWOPSY,THPKAP
         INPR = 0
         DO 300 IG = 1,NKAPT
            IF ( XKAP(IG) .LE. THPKAP ) THEN
               K = JWOP(1,IG)
               L = JWOP(2,IG)
               WRITE (LUPRIN,1100) IG,K,L,XKAP(IG)
            ELSE
               INPR = INPR + 1
            END IF
  300    CONTINUE
         IF (INPR.GT.0) WRITE (LUPRIN,1230) INPR,THPKAP
      END IF
C
C
 9999 CONTINUE
      RETURN
C
 1010 FORMAT(/5X,'Orbital operator symmetry =',I2,/,
     *  5X,'(only elements abs greater than',1P,D8.1,
     *     ' times max abs value)'//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1020 FORMAT(/5X,'Orbital operator symmetry =',I2,//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1030 FORMAT(/5X,'Orbital operator symmetry =',I2,/,
     *  5X,'(only elements less than',1P,D10.2,')'//
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1100 FORMAT(I15,I7,I3,F20.10)
 1210 FORMAT(/1X,I8,' elements with absolute value less than',
     *       1P,D9.2,' not printed.')
 1220 FORMAT(/1X,I8,' zero elements not printed.')
 1230 FORMAT(/1X,I8,' elements with value greater than',
     *       1P,D9.2,' not printed.')
C
C -- end of subroutine PRKAP
C
      END
C  /* Deck prwop */
      SUBROUTINE PRWOP (WOP,IUNITX)
C
C Written 11-Apr-1984 -- hjaaj
C Last revision 17-Oct-1984 hjaaj (IUNITX.lt.0 option)
C                9-May-1984 hjaaj
C
C Purpose:
C  Write orbital operator array WOP to unit IUNIT.
C
C  IUNIT = IABS(IUNITX)
C
C  If IUNITX .lt. 0, write all elements;
C  if IUNITX .gt. 0, write only elements with absolute value
C                    .ge. PRFAC * (max. absolute value)
C
#include <implicit.h>
      DIMENSION WOP(*)
      PARAMETER (PRFAC = 0.1D0, D0 = 0.0D0)
C
C Used from common blocks:
C   INFVAR : NWOPT,JWOPSY,JWOP(2,*)
C
#include <maxorb.h>
#include <infvar.h>
C
      IUNIT = IABS(IUNITX)
C
C Any elements?
C
      IF (NWOPT.LE.0) GO TO 300
C
C 1. find maximum absolute value in WOP
C
      IF (IUNITX.GT.0) THEN
         MAX = IDAMAX(NWOPT,WOP,1)
         THPWOP = PRFAC*ABS(WOP(MAX))
         WRITE (IUNIT,1000) JWOPSY,PRFAC
      ELSE
         WRITE (IUNIT,1010) JWOPSY
         THPWOP = D0
      END IF
      INPR = 0
      DO 200 IG = 1,NWOPT
         IF (ABS(WOP(IG)).GE.THPWOP) THEN
            K = JWOP(1,IG)
            L = JWOP(2,IG)
            WRITE (IUNIT,1100) IG,K,L,WOP(IG)
         ELSE
            INPR = INPR + 1
         END IF
  200 CONTINUE
      IF (INPR.GT.0) WRITE (IUNIT,1200) INPR,THPWOP
C
      RETURN
C
  300 CONTINUE
      WRITE (IUNIT,3000) JWOPSY
      RETURN
C
 1000 FORMAT(/5X,'Orbital operator symmetry =',I2,/,
     *  5X,'(only elements abs greater than',1P,D8.1,
     *     ' times max abs value)'//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1010 FORMAT(/5X,'Orbital operator symmetry =',I2,//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1100 FORMAT(I15,I7,I3,F20.10)
 1200 FORMAT(/1X,I8,' elements with absolute value less than',
     *       1P,D9.2,' not printed.')
 3000 FORMAT(/5X,'Orbital operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C -- end of subroutine PRWOP
C
      END
C  /* Deck readmo */
      SUBROUTINE READMO(CMO,JRDMO,WRK,LFREE)
C
C Last revision 7-Jan-1985 hjaaj -- ver. 2.6
C  8-Jul-1992 hjaaj (removed H1VIRT and FCMO calls)
C  6-May-1994 hjaaj (removed IRDMO=8 call H1MO)
C
C Purpose:Read MO coefficients from IUNIT or generate MO coefficients.
C         IRDMO = abs(JRDMO) determines how they are read.
C
C negative JRDMO means NO STOP but return in JRDMO number of errors
C
C IRDMO    format
C -----    ------
C   3      normal formatted input: (4F18.14)
C   4      formatted input: (6F12.8)
C   7      from LUIT1 with label "OLDORB"
C   9      from LUIT1 with label "NEWORB"
C  10      from SIRIFC (LUSIFC)
C
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*),WRK(*)
      LOGICAL   FOUND
C
#include <priunit.h>
#include <inforb.h>
#include <inftap.h>
C
      LOGICAL      NOSTOP
      CHARACTER*4  LABEL(104), IFLGC, ISTAR, IEVEC
      CHARACTER*8  MOLABL
      CHARACTER*10 PFMT
      DATA ISTAR/'****'/, IEVEC/'EVEC'/
C
      CALL QENTER('READMO')
C
      IRDMO = ABS(JRDMO)
      IF (JRDMO .LT. 0) THEN
         JRDMO = 0
         NOSTOP = .TRUE.
      ELSE
         NOSTOP = .FALSE.
      END IF
      IF (IRDMO.LT.1 .OR. IRDMO.GT.10) THEN
         WRITE (LUPRI,10000) IRDMO
         WRITE (LUERR ,10000) IRDMO
10000    FORMAT(/' *** READMO-ERROR *** IRDMO =',I4,' out of range')
         CALL QTRACE(LUPRI)
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('*** ERROR *** illegal option in READMO')
         END IF
      END IF
      GO TO (9000,9000,300,300,8000,8000,700,8000,700,1000),IRDMO
C
C-- IRDMO=3 or 4:
C
  300 IF (IRDMO .EQ. 3) THEN
         PFMT = '(4F18.14)'
      ELSE
         PFMT = '(6F12.8)'
      END IF
C
      JLUCMD = LUCMD
      IF (JLUCMD .LT. 0) CALL GPOPEN(LUCMD,'DALTON.INP',
     &   'OLD',' ','FORMATTED',IDUMMY,.FALSE.)
      REWIND (LUCMD)
  310 READ (LUCMD,'(A)',END=400,ERR=400) MOLABL
      IF (MOLABL .NE. '**MOLORB' .AND. MOLABL .NE. '**NATORB') GO TO 310
      IEND=0
      DO 340 ISYM=1,NSYM
         NBASI=NBAS(ISYM)
      IF(NBASI.EQ.0) GO TO 340
         NORBI=NORB(ISYM)
         DO 335 NI=1,NBASI
         IF (NI.LE.NORBI) THEN
            IST=IEND+1
            IEND=IEND+NBASI
            READ(LUCMD,PFMT) (CMO(I),I=IST,IEND)
         ELSE
            READ(LUCMD,PFMT) (DUM,I=1,NBASI)
         END IF
  335    CONTINUE
  340 CONTINUE
      IF (JLUCMD .LT. 0) CALL GPCLOSE(LUCMD,'KEEP')
      GO TO 9000
  400 CONTINUE
         WRITE (LUPRI,'(//A)')
     *   ' >>> FATAL ERROR, label "**MOLORB" not found on SIRIUS input.'
         CALL QTRACE(LUPRI)
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT(
     &         'ERROR (READMO) label "**MOLORB" not found in input.')
         END IF
C
C-- IRDMO = 7 or 9:
C
  700 CONTINUE
      JLUIT1 = LUIT1
      IF (JLUIT1 .LE. 0) CALL GPOPEN(LUIT1,'SIRIUS.RST','OLD',' ',
     &   'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUIT1
      IF (IRDMO.EQ.7) THEN
         MOLABL = 'OLDORB  '
      ELSE
         MOLABL = 'NEWORB  '
      END IF
      LBLERR = -1
      CALL MOLLAB(MOLABL,LUIT1,LBLERR)
      IF (LBLERR .GE. 0) THEN
         CALL READT(LUIT1,NCMOT,CMO)
      ELSE IF (LBLERR .EQ. -1) THEN
         WRITE (LUPRI,'(//3A)')
     & ' >>> FATAL ERROR, MO label "',MOLABL,'" not found on SIRIUS.RST'
         CALL QTRACE(LUPRI)
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('ERROR (READMO) MO coefficients not found.')
         END IF
      ELSE
         WRITE (LUPRI,'(//A/2A)')
     &   ' >>> FATAL ERROR, could not read SIRIUS.RST',
     &   '                  when searching for mo label :',MOLABL
         CALL QTRACE(LUPRI)
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('ERROR (READMO) COULD NOT READ MO FILE.')
         END IF
      END IF
      IF (JLUIT1 .LE. 0) CALL GPCLOSE(LUIT1,'KEEP')
      GO TO 9000
C
C-- IRDMO = 10:
C
 1000 CONTINUE
      CALL RD_SIRIFC('CMO',FOUND,CMO,WRK,LFREE)
      IF (.NOT.FOUND) THEN
         IF (NOSTOP) THEN
            JRDMO = JRDMO + 1
         ELSE
            CALL QUIT('ERROR (READMO) MO''s not found on SIRIFC')
         END IF
      END IF
      GO TO 9000
C
C--   IRDMO 5,6,8 not defined any more.
C
 8000 CONTINUE
C
      WRITE (LUPRI,'(//A/A)')
     *   ' >>> FATAL INTERNAL ERROR, IRDMO options 5, 6, 8',
     *   '     is not defined in READMO since May-95'
      IF (NOSTOP) THEN
         JRDMO = JRDMO + 1
      ELSE
         CALL QUIT('ERROR (READMO) options 5,6,8 is not defined')
      END IF
C
C--
C
 9000 CONTINUE
C
      CALL QEXIT('READMO')
      RETURN
      END
C  /* Deck neworb */
      SUBROUTINE NEWORB(LABEL3,CMO,REWIT1)
C
C  19-Oct-1989; 3-Aug-1990 (REWIT1; hjaaj)
C
C  If (REWIT1) then
C     Write molecular orbitals with label 'NEWORB  ' at
C     beginning of LUIT1
C  else
C     Overwrite 'NEWORB  ' on LUIT1
C  end if
C
#include <implicit.h>
      CHARACTER*8 LABEL3
      DIMENSION CMO(NCMOT)
      LOGICAL   REWIT1
C
#include <dummy.h>
C
C Used from common blocks:
C  INFORB : NCMOT
C  INFTAP : LUIT1
C
#include <inforb.h>
#include <inftap.h>
C
      LOGICAL FNDLAB
      CHARACTER*8 MOLBL(3),LABTBL(3)
      DATA MOLBL /'********','        ','        '/
      DATA LABTBL/'BASINFO ','NEWORB  ','EODATA  '/
C
      CALL QENTER('NEWORB')
C
C     Write orbitals to LUIT1
C
      CALL GETDAT(MOLBL(2),MOLBL(3))
      NCMOT4 = MAX(4,NCMOT)
      IF (REWIT1) THEN
         CALL NEWIT1
      ELSE
         REWIND LUIT1
         IF (.NOT.FNDLAB(LABTBL(1),LUIT1)) THEN
C           no "BASINFO ",
C           must be a new SIRIUS.RST file, save BASINFO on it:
            CALL NEWIT1
            GO TO 100
         END IF
         REWIND LUIT1
         IF (FNDLAB(LABTBL(2),LUIT1)) THEN
C           found "NEWORB  ",
C           backspace before writing new "NEWORB  " label
            BACKSPACE LUIT1
            GO TO 100
         END IF
         REWIND LUIT1
         IF (FNDLAB(LABTBL(3),LUIT1)) THEN
C           no "NEWORB  " but found "EODATA  ",
C           backspace before writing "NEWORB  " label
            BACKSPACE LUIT1
            GO TO 100
         END IF
C
C        no "NEWORB  " and no "EODATA  "
C        is positioned at end-of-file on LUIT1 after last FNDLAB
C        /hjaaj June 09
      END IF

  100 WRITE (LUIT1) MOLBL(1),MOLBL(2),LABEL3,LABTBL(2)
C     ...label 'NEWORB  '
      CALL WRITT(LUIT1,NCMOT4,CMO)
      WRITE (LUIT1) MOLBL,LABTBL(3)
C     ...label 'EODATA  '
      REWIND LUIT1
C
      CALL QEXIT('NEWORB')
      RETURN
      END
C  /* Deck newit1 */
      SUBROUTINE NEWIT1
C
C 2-Dec-1992 Hans Joergen Aa. Jensen
C
C Write basis set information at beginning of LUIT1
C 950825-hjaaj: also save NRHF(*) and IOPRHF for use with AUTOCC
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB: NSYM, NORB, NBAS, NRHF
C  SCBRHF: IOPRHF
C  INFTAP: LUIT1
C
#include <inforb.h>
#include <scbrhf.h>
#include <inftap.h>
C
C
      INTEGER*4   NSYM_i4, NBAS_i4(8),NORB_i4(8), NRHF_i4(8), IOPRHF_i4
      ! hjaaj Apr 2011: such that Dalton with 8 byte integers can read a
      ! SIRIUS.RST from a Dalton with 4 byte integers, and vice versa.
C
      CHARACTER*8 BASINFO(4)
      DATA BASINFO(1)/'********'/, BASINFO(4)/'BASINFO '/
C
      CALL GETDAT(BASINFO(2),BASINFO(3))
      REWIND LUIT1
      WRITE (LUIT1) BASINFO
      NSYM_i4    = NSYM
      NBAS_i4(:) = NBAS(:)
      NORB_i4(:) = NORB(:)
      NRHF_i4(:) = NRHF(:)
      IOPRHF_i4  = IOPRHF
      WRITE (LUIT1) NSYM_i4,NBAS_i4,NORB_i4,NRHF_i4,IOPRHF_i4
      RETURN
      END
C  /* Deck getac */
      SUBROUTINE GETAC(FC,FCAC)
C
C  7-Feb-1987 Hans Joergen Aa. Jensen
C
C Module : SIRGP
C
C Purpose: Extract block with active-active orbital indices out
C          of symmetry packed matrix (e.g. the Fock matrix).
C
#include <implicit.h>
      DIMENSION FC(*), FCAC(*)
C
C INFORB : NASHT, NNASHX, NSYM, NASH(*), IIORB(*),
C          NISH(*), NOCC(*)
C INFIND : IROW(*), ICH(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('GETAC ')
      IF (NASHT .EQ. 0) GO TO 9999
C
      CALL DZERO(FCAC,NNASHX)
      DO 300 ISYM=1,NSYM
      IF (NASH(ISYM).EQ.0) GO TO 300
         NISHI = NISH(ISYM)
         NOCI  = NOCC(ISYM)
         IJFC  = IIORB(ISYM) + IROW(NISHI+1)
         IORBI = IORB(ISYM)
         DO 290 I = (NISHI+1),NOCI
            NI     = ICH(IORBI+I)
            IROWNI = IROW(NI)
            IJFC   = IJFC + NISHI
            DO 280 J = (NISHI+1),I
               IJFC = IJFC + 1
               NJ   = ICH(IORBI+J)
               IF (NI.GE.NJ) THEN
                  NIJ = IROWNI   + NJ
               ELSE
                  NIJ = IROW(NJ) + NI
               END IF
               FCAC(NIJ) = FC(IJFC)
  280       CONTINUE
  290    CONTINUE
  300 CONTINUE
C
 9999 CALL QEXIT('GETAC ')
      RETURN
C
C     End of GETAC.
C
      END
C  /* Deck getac1 */
      SUBROUTINE GETAC1(HH,HHAC)
C
C  2-Nov-1989 Hans Joergen Aa. Jensen
C
C Module : SIRGP
C
C Purpose: Extract block with active-active orbital indices out
C          of full matrix.
C
#include <implicit.h>
      DIMENSION HH(NORBT,NORBT), HHAC(NASHT,NASHT)
C
C INFORB : NORBT, NASHT, N2ASHX
C INFIND : IOBTYP(*), ICH(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('GETAC1')
      IF (NASHT .EQ. 0) GO TO 9999
C
      DO 200 J = 1,NORBT
      IF (IOBTYP(J) .EQ. JTACT) THEN
         NJ = ICH(J)
         DO 100 I = 1,NORBT
         IF (IOBTYP(I) .EQ. JTACT) THEN
            NI = ICH(I)
            HHAC(NI,NJ) = HH(I,J)
         END IF
  100    CONTINUE
      END IF
  200 CONTINUE
C
 9999 CALL QEXIT('GETAC1')
      RETURN
C
C     End of GETAC1.
C
      END
C  /* Deck getac2 */
      SUBROUTINE GETAC2(HH,HHAC)
C
C  6-May-1990 Hans Joergen Aa. Jensen
C
C Module : SIRGP
C
C Purpose: Extract NNASHX block with active-active orbital indices out
C          of NNORBX matrix.
C
#include <implicit.h>
      DIMENSION HH(NNORBX), HHAC(NNASHX)
C
C INFORB : NORBT, NASHT, NNORBX,NNASHX
C INFIND : IOBTYP(*), ICH(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('GETAC2')
      IF (NASHT .EQ. 0) GO TO 9999
C
      DO 200 J = 1,NORBT
      IF (IOBTYP(J) .EQ. JTACT) THEN
         NJ  = ICH(J)
         IRNJ= IROW(NJ)
         IRJ = IROW(J)
         DO 100 I = 1,J
         IF (IOBTYP(I) .EQ. JTACT) THEN
            NI = ICH(I)
            HHAC(IRNJ+NI) = HH(IRJ+I)
         END IF
  100    CONTINUE
      END IF
  200 CONTINUE
C
 9999 CALL QEXIT('GETAC2')
      RETURN
C
C     End of GETAC1.
C
      END
C  /* Deck pksym1 */
      SUBROUTINE PKSYM1(ASP,APK,NDIM,NBLK,IWAY)
C
C Jan 90/May 96 Hans Joergen Aa. Jensen
C
C If IWAY .ge. 0 then
C    (in this case ASP and APK may overlap as long as 
C     loc(APK) .le. loc(ASP))
C    If IWAY .eq. 2 then
C       Pack NNDIMX matrix ASP of symmetry 1 to ASP in NNDIMT format
C       i.e. call pksym1(asp,asp,ndim,nblk,+2)
C    else
C       Pack NNDIMX matrix ASP of symmetry 1 to APK in NNDIMT format
C    end if
C else
C    Unpack from APK to ASP
C
#include <implicit.h>
      DIMENSION ASP(*), APK(*), NDIM(NBLK)
      IF (IWAY .LT. 0) THEN
         NDIMT  = ISUM(NBLK,NDIM,1)
         NNDIMX = (NDIMT*NDIMT+NDIMT)/2
         CALL DZERO(ASP,NNDIMX)
      END IF
      IF (IWAY .EQ. 2) THEN
C         block 1 is not moved if equivalence (asp,apk)
         IBLK1  = 2
         IDIMI  = NDIM(1)
         ISPOFF = IDIMI*(IDIMI+1)/2
         IPKOFF = ISPOFF
      ELSE
         IBLK1  = 1
         IDIMI  = 0
         ISPOFF = 0
         IPKOFF = 0
      END IF
      DO 800 IBLK = IBLK1,NBLK
         NDIMI = NDIM(IBLK)
         DO 600 J = 1,NDIMI
            ISPOFF = ISPOFF + IDIMI
            IF (IWAY .GE. 0) THEN
C              no synchronization problems as long as
C              loc(apk) .le. loc(asp) (or no overlap)
               DO I = 1,J
                  APK(IPKOFF+I) = ASP(ISPOFF+I)
               END DO
            ELSE
               DO I = 1,J
                  ASP(ISPOFF+I) = APK(IPKOFF+I)
               END DO
            END IF
            ISPOFF = ISPOFF + J
            IPKOFF = IPKOFF + J
  600    CONTINUE
         IDIMI = IDIMI + NDIMI
  800 CONTINUE
      RETURN
      END
C  /* Deck rdonel */
      SUBROUTINE RDONEL(KEY,FOUND,A,NA)
C
C 31-Jul-1987 hjaaj
C
C     Purpose : Read LUONEL
C
C     Note: if FOUND true, it must not be updated.
C
C
#include "implicit.h"
#include "priunit.h"
C     Global array
      DIMENSION A(NA,*)
      CHARACTER KEY*(*)
      LOGICAL   FOUND, RDBUF
C
C Used from common blocks:
C   INFINP : TITMOL(2), POTNUC, CENT(*), TYPE(*)
C   INFORB : NBAS(8), NSYM, NBAST, NNBAST
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "inftap.h"
#include "r12int.h"
C     Local arrays
      PARAMETER ( LBUF = 600 )
      DIMENSION BUF(LBUF), DBUF(LBUF,3), QBUF(LBUF,6), IBUF(LBUF)
      DIMENSION ITMPMO(2*MXCORB)
      EQUIVALENCE (BUF,DBUF,QBUF,ITMPMO)
      LOGICAL   OPEND, FEXIST, ERRSTP
C
      CALL QENTER('RDONEL')
      ERRSTP = FOUND
      IF (.NOT.FOUND) FOUND  = .TRUE.
C     Note: if FOUND true, it must not be updated.
      RDBUF  = .FALSE.
C
      CALL GPINQ(FNONEL,'OPENED',OPEND)
      IF (.NOT.OPEND) THEN
         CALL GPINQ(FNONEL,'EXIST',FEXIST)
       IF (FEXIST) THEN
         CALL GPOPEN(LUONEL,FNONEL,'OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         OPEND = .TRUE.
         IF (LBONEL .LT. 0) LBONEL = LBUF
       ELSE
         WRITE (LUPRI,'(/3A)')
     *      ' RDONEL ERROR: LUONEL file ',FNONEL,' does not exist'
         WRITE (LUERR ,'(/3A)')
     *      ' RDONEL ERROR: LUONEL file ',FNONEL,' does not exist'
         IF (ERRSTP) THEN
            CALL QUIT('RDONEL ERROR: AOONEINT file DOES NOT EXIST')
         ELSE
            FOUND = .FALSE.
            GO TO 9999
         END IF
       END IF
      END IF
C
      IF     (KEY.EQ.'OPEN    ') THEN
C        ... file is opened automatically first time RDONEL is called.
      ELSE IF(KEY.EQ.'CLOSE   ') THEN
         CALL GPCLOSE(LUONEL,'KEEP')
      ELSE IF(KEY.EQ.'MLCLINFO') THEN
         REWIND LUONEL
         READ (LUONEL) TITMOL
         READ (LUONEL) NSYM,(NBAS(ISYM),ISYM=1,NSYM),POTNUC
         NBAST  = 0
         NNBAST = 0
         DO 110 I = 1,NSYM
            NBAST  = NBAST  + NBAS(I)
            NNBAST = NNBAST + NBAS(I)*(NBAS(I)+1)/2
  110    CONTINUE
         DO 120 I = NSYM+1,8
            NBAS(I) = 0
  120    CONTINUE
         CALL MOLLAB('SYMINPUT',LUONEL,LUPRI)
         READ (LUONEL) NBT,(ITMPMO(I),I=1,2*NBT)
C        READ (LUONEL) NBT,(ICENT(I),I=1,NBT),(ITYPE(I),I=1,NBT)
         IF (NBT .NE. NBAST) THEN
            WRITE (LUPRI,'(//A,2(/A,I5))')
     *      ' >>> RDONEL ERROR',
     *      '     NBAST from first record on AOONEINT :',NBAST,
     *      '     NBAST from "SYMINPUT"   on AOONEINT :',NBT
            CALL QUIT('RDONEL ERROR, two different NBAST on AOONEINT')
         END IF
         REWIND LUONEL
         DO 130 I=1,NBAST
            WRITE(CENT(I),'(A4)') ITMPMO(I)
            WRITE(TYPE(I),'(A4)') ITMPMO(NBAST+I)
  130    CONTINUE
      ELSE IF(KEY.EQ.'OVERLAP ') THEN
C        Read overlap, in symmetry blocked triangular form
         REWIND LUONEL
         CALL MOLLAB('OVERLAP ',LUONEL,LUPRI)
         RDBUF = .TRUE.
      ELSE IF(KEY.EQ.'ONEHAMIL') THEN
C        Read one electron part, in symmetry blocked triangular form
         REWIND LUONEL
         CALL MOLLAB('ONEHAMIL',LUONEL,LUPRI)
         RDBUF = .TRUE.
      ELSE IF(KEY.EQ.'KINETINT') THEN
         REWIND LUONEL
         CALL MOLLAB('KINETINT',LUONEL,LUPRI)
         RDBUF = .TRUE.
      ELSE IF(KEY.EQ.'DIPOLMOM') THEN
         REWIND LUONEL
         CALL MOLLAB('DIPOLMOM',LUONEL,LUPRI)
         IF (LBONEL .NE. LBUF) THEN
            WRITE(LUPRI,'(/A,2(/A,I8))')
     *         ' ERROR RDONEL, LBONEL .ne. local LBUF for DIPOLMOM',
     *         ' LBONEL =',LBONEL,' LBUF   =',LBUF
            CALL QUIT('RDONEL ERROR, LBONEL .ne. LBUF for key DIPOLMOM')
         END IF
  501    READ (LUONEL) DBUF,IBUF,NUT
         IF (NUT.EQ.0)   GO TO 501
         IF (NUT.LT.0)   GO TO 502
         IF (NOAUXB) THEN
            DO I=1,NUT
               CALL INOAUX(IBUF(I),IR12)
               IF (IR12 .GT. 0) THEN
                  A(IR12,1) = DBUF(I,1)
                  A(IR12,2) = DBUF(I,2)
                  A(IR12,3) = DBUF(I,3)
               END IF
            END DO 
         ELSE
            DO I=1,NUT
               J = IBUF(I)
               A(J,1) = DBUF(I,1)
               A(J,2) = DBUF(I,2)
               A(J,3) = DBUF(I,3)
            END DO 
         END IF
         GO TO 501
  502    CONTINUE
      ELSE IF(KEY.EQ.'QUADRP  ') THEN
         REWIND LUONEL
         CALL MOLLAB('QUADRP  ',LUONEL,LUPRI)
         IF (LBONEL .NE. LBUF) THEN
            WRITE(LUPRI,'(/A,2(/A,I8))')
     *         ' ERROR RDONEL, LBONEL .ne. local LBUF for QUADRP  ',
     *         ' LBONEL =',LBONEL,' LBUF   =',LBUF
            CALL QUIT('RDONEL ERROR, LBONEL .ne. LBUF for key QUADRP')
         END IF
  601    READ (LUONEL) QBUF,IBUF,NUT
         IF (NUT.EQ.0)   GO TO 601
         IF (NUT.LT.0)   GO TO 602
         IF (NOAUXB) THEN
            DO I=1,NUT
               CALL INOAUX(IBUF(I),IR12)
               IF (IR12 .GT. 0) THEN
                  DO ICOMP = 1,6
                     A(IR12,ICOMP) = QBUF(I,ICOMP)
                  END DO
               END IF
            END DO
         ELSE
            DO I=1,NUT
               J = IBUF(I)
               DO ICOMP = 1,6
                  A(J,ICOMP) = QBUF(I,ICOMP)
               END DO
            END DO
         END IF
         GO TO 601
  602    CONTINUE
      ELSE IF(KEY.EQ.'ISORDK  ') THEN
         WRITE(LUPRI,'(2A)') KEY, ' IS NOT IMPLEMENTED YET IN RDONEL'
         IF (ERRSTP) THEN
            CALL QUIT('RDONEL ERROR, KEY ISORDK NOT IMPLEMENTED YET.')
         ELSE
            FOUND = .FALSE.
         END IF
      ELSE IF(KEY.EQ.'SCFINP  ') THEN
         WRITE(LUPRI,'(2A)') KEY, ' IS NOT IMPLEMENTED YET IN RDONEL'
         IF (ERRSTP) THEN
            CALL QUIT('RDONEL ERROR, KEY SCFINP NOT IMPLEMENTED YET.')
         ELSE
            FOUND = .FALSE.
         END IF
      ELSE
         WRITE(LUPRI,'(2A)') KEY, ' IS AN INVALID KEY TO RDONEL'
         IF (ERRSTP) THEN
            CALL QUIT('ERROR RDONEL, INVALID KEY')
         ELSE
            FOUND = .FALSE.
         END IF
      END IF
      IF (RDBUF) THEN
         IF (LBONEL .GT. LBUF) THEN
            WRITE(LUPRI,'(/A,2(/A,I8))')
     *         ' ERROR RDONEL, LBONEL .GT. local LBUF',
     *         ' LBONEL =',LBONEL,' LBUF   =',LBUF
            CALL QUIT('ERROR RDONEL, LBONEL .gt. LBUF')
         END IF
         CALL DZERO(A,NNBAST)
 2100    READ (LUONEL) (BUF(I),I=1,LBONEL),(IBUF(I),I=1,LBONEL),LENGTH
         IF (NOAUXB) THEN
            DO I = 1, LENGTH
               CALL INOAUX(IBUF(I),IR12)
               IF (IR12 .GT. 0) A( IR12 , 1 ) = BUF(I)
            END DO
         ELSE
            DO I = 1, LENGTH
               A( IBUF(I) , 1 ) = BUF(I)
            END DO
         END IF
         IF (LENGTH .GE. 0) GO TO 2100
         REWIND LUONEL
      END IF
      CALL GPCLOSE(LUONEL,'KEEP')
C
 9999 CALL QEXIT('RDONEL')
      RETURN
      END
C  /* Deck ijkaux */
      SUBROUTINE IJKAUX(II,JJ,KK,LL)
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
C ... only called from cc/ routines,
C     thus NSYM from ccorb.h and not from inforb.h /hjaaj checked aug 04
#include <maxorb.h>
#include <r12int.h>
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DIMENSION IJ(MXCORB)
      SAVE IJ
      IF (FIRST) THEN
         KL = 0
         MN = 0
         DO ISYM = 1, NSYM
            DO I = 1, MBAS1(ISYM) + MBAS2(ISYM)
               KL = KL + 1
               IF (I .LE. MBAS1(ISYM)) THEN
                  MN = MN + 1
                  IJ(KL) = MN
               ELSE
                  IJ(KL) = -1
               END IF
            END DO
         END DO
         FIRST = .FALSE.
      END IF
      II = IJ(II)
      JJ = IJ(JJ)
      KK = IJ(KK)
      LL = IJ(LL)
      RETURN
      END
C  /* Deck inoaux */
      SUBROUTINE INOAUX(II,JJ)
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
C ... only used in RDONEL when called from cc/ routines,
C     thus NSYM from ccorb.h and not from inforb.h /hjaaj checked aug 04
#include <maxorb.h>
#include <r12int.h>
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DIMENSION IJ(MXORBT)
      SAVE IJ
      IF (FIRST) THEN
         KL = 0
         MN = 0
         DO ISYM = 1, NSYM
            DO I = 1, MBAS1(ISYM) + MBAS2(ISYM)
               DO J = 1, I
                  KL = KL + 1
                  IF (I. LE. MBAS1(ISYM) .AND. J .LE. MBAS1(ISYM)) THEN
                     MN = MN + 1
                     IJ(KL) = MN
                  ELSE
                     IJ(KL) = -1
                  END IF
               END DO
            END DO
         END DO
         FIRST = .FALSE.
      END IF
      JJ = IJ(II)
      RETURN
      END
C  /* Deck rdsupm */
      SUBROUTINE RDSUPM(NSIM,FMAT,DMAT,WORK,LFREE)
C
C     Jan 90, Hans Joergen Aa. Jensen
C     l.r. 2-Apr-1997 hjaaj
C     (FMAT and DMAT now (NBAST,NBAST))
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FMAT(*), DMAT(*), WORK(*)
C
C Used from common blocks:
C  INFORB : NBAST,N2BASX
C  INFTAP : LUSUPM
C  PRIUNIT : IPRSTAT
C  CBIHRS : parameters for FORMSUP
C
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <cbihrs.h>
C
C     Local variables:
C
      LOGICAL FNDLAB

      CALL QENTER('RDSUPM')
C
      IF (LUSUPM .EQ. -1) THEN
C     ... special code for AOSUPMAT not to be used
         WRITE (LUPRI,'(/A/A)')
     &      'PROGRAMMING ERROR: RDSUPM called with LUSUPM = -1',
     &      'Please report to dalton-admin.'
         CALL QTRACE(LUPRI)
         CALL QUIT('Programming ERROR: RDSUPM called with LUSUPM.eq.-1')
      END IF
C
      I = 0
 100     IF (LUSUPM .LE. 0) CALL GPOPEN(LUSUPM,FNSUPM,
     &      'UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
         I = I + 1
         CALL REWSPL(LUSUPM)
         IF (.NOT. FNDLAB('PSUPMAT ',LUSUPM)) THEN
            IF (I.GT.1) CALL QUIT('CALL FORMSUP failed')
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL FORMSUP(WORK,LFREE,NOSSUP,THRSUP,IPRSUP)
            CALL TIMER('FORMSUP',TIMSTR,TIMEND)
            CALL FLSHFO(LUPRI)
            GO TO 100
         END IF
C
C     Fold density matrices
C
      IF (LFREE .LT. N2BASX) CALL ERRWRK('RDSUPM',-N2BASX,LFREE)
      DO I = 1,NSIM
         JDG = 1 + (I-1)*N2BASX
         CALL DCOPY(N2BASX,DMAT(JDG),1,WORK,1)
         CALL DGEFSP(NBAST,WORK,DMAT(JDG))
      END DO
C
C     note equivalence(WORK, IWORK)
C
      CALL RDSUPP(NSIM,FMAT,DMAT,WORK,WORK,LFREE)
C
C     square Fock matrices, restore DCAO
C
      DO I = 1,NSIM
         JDG = 1 + (I-1)*N2BASX
         CALL DCOPY(NNBASX,FMAT(JDG),1,WORK,1)
         CALL DSPTGE(NBAST,WORK,FMAT(JDG))
         CALL DCOPY(NNBASX,DMAT(JDG),1,WORK,1)
         CALL DUNFLD(NBAST,WORK,DMAT(JDG))
      END DO
C
      CALL GPCLOSE(LUSUPM,'KEEP')
      CALL QEXIT('RDSUPM')
      RETURN
      END
C  /* Deck rdsupp */
#include <vdir.h>
      SUBROUTINE RDSUPP(NSIM,FMAT,DMAT,WORK,IWORK,LFREE)
C
C Jan 90, Hans Joergen Aa. Jensen
C
C NOTE: WORK and IWORK must be equivalenced.
C
#include <implicit.h>
#include <iratdef.h>
#include <thrzer.h>
      PARAMETER (D0 = 0.0D0)
C
#if defined (VAR_SPLITFILES)
#include <dummy.h>
#include <2gbdef.h>
      CHARACTER*80 FNNAME, FNNM2
      DIMENSION ITMP(7)
#include <chrnos.h>
#endif
      DIMENSION FMAT(N2BASX,NSIM), DMAT(N2BASX,NSIM)
      DIMENSION WORK(LFREE), IWORK(IRAT*LFREE)
C
C Used from common blocks:
C
C   INFORB : N2BASX,MULD2H(8,8),...
C   INFIND : IROW(:),ISAO(:)
C   INFVAR : JWOPSY
C   INFTAP : LUSUPM
C   PRIUNIT : LUERR, IPRSTAT
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL NOSYM
      DIMENSION MBAS(8)
C
      CALL QENTER('RDSUPP')
      CALL REWSPL(LUSUPM)
      CALL MOLLAB('PSUPMAT ',LUSUPM,LUPRI)
C
C *** Test if LUSUPM file is OK
C
      READ (LUSUPM) NOSYM, MSYM, MBAS
      NERR = 0
      IF (JWOPSY .NE. 1 .AND. .NOT.NOSYM) NERR = NERR + 1
      IF (MSYM .NE. NSYM) THEN
         NERR = NERR + 10
      ELSE
         DO 10 ISYM = 1,NSYM
   10    IF (MBAS(ISYM) .NE. NBAS(ISYM)) NERR = NERR + 10
      END IF
      IF (NERR .GT. 0) THEN
         WRITE (LUPRI,'(//5X,A)')
     &      'RDSUPP ERROR: AOSUPINT inconsistent with Sirius input'
         IF (MOD(NERR,10) .NE. 0) WRITE (LUPRI,'(/5X,2A)')
     &      'AOSUPINT must be sorted with NOSYM when ',
     &      'perturbation symmetry .ne. 1'
         IF (NERR .GE. 10) WRITE (LUPRI,'(/5X,A,2I3/,(5X,A,8I5))')
     &      'AOSUPINT vs. AOONEINT, NSYM =',MSYM,NSYM,
     &      'NBAS from AOSUPINT :',(MBAS(I),I=1,8),
     &      'NBAS from AOONEINT :',(NBAS(I),I=1,8)
         CALL QTRACE(LUPRI)
         CALL QUIT('RDSUPP ERROR: AOSUPINT inconsistent with AOONEINT')
      END IF
      IF (IPRSTAT .GT. 5) THEN
         WRITE (LUERR,*) 'Test output from RDSUPP, NOSYM =',NOSYM
      END IF
C
C
C *** Start loop over psupmat file
C
      IREC  = 0
      IREC1 = 0
      IREC2 = 0
      JBUF = 256 + 144*IRAT
      MBUF = 0
  100 CONTINUE
         IREC = IREC + 1
#if defined (VAR_SPLITFILES)
C
         IF ((JBUF + 8*IRAT) .GT. I2GB) THEN
C
C     Ooops, file will be overfull, need to open a new one......
C
            INQUIRE(UNIT=LUSUPM,NAME=FNNAME)
            LN = 1
 14         CONTINUE
            IF (FNNAME(LN:LN) .NE. ' ') THEN
               LN = LN + 1
               GOTO 14
            END IF
            LN = LN - 1
            CALL GPCLOSE(LUSUPM,'KEEP')
            I = LN - 1
            IF (FNNAME(I:I) .NE. '-') THEN
               FNNM2 = FNNAME(1:LN)//'-0'
               LN = LN + 2
            ELSE
               READ(FNNAME(LN:),'(I1)') INUM
               INUM = INUM + 1
               IF (INUM .GT. 9) THEN
                  WRITE (LUPRI,'(/A)') ' DALTON needs to split '//
     &                 ' a file more than 11 times.',
     &                 ' This is currently not supported'
                  CALL QUIT('Too many splittings of a file')
               END IF
               FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
            END IF
            CALL GPOPEN(LUSUPM,FNNM2(1:LN),'UNKNOWN',' ',' ',
     &           IDUMMY,.FALSE.)
            JBUF = 0
         END IF
         CALL READSI(LUSUPM,7,ITMP,MBUF)
         ITYP = ITMP(1)
         NP1  = ITMP(2)
         NQ1  = ITMP(3)
         NPL  = ITMP(4)
         NQL  = ITMP(5)
         IBEF = ITMP(6)
         NPQRS= ITMP(7)
         JBUF = JBUF + 8*IRAT
#else
         READ (LUSUPM) ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS
#endif
         IF (IPRSTAT .GT. 20) THEN
            WRITE (LUERR,*) 'ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS'
            WRITE (LUERR,'(7I10)') ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS
         END IF
C        ITYP = -2 finished
C                1 P(IBEF+1:IBEF+NPQRS)
C                2 P(1:NPQRS), INDP(1:NPQRS) with (RS) addresses
         IF (ITYP .EQ. -2) GO TO 900
C
         ISP1 = ISAO(NP1)
         ISPL = ISAO(NPL)
         ISQ1 = ISAO(NQ1)
         ISQL = ISAO(NQL)
         IF (ISP1.EQ.ISPL .AND. ISQ1.EQ.ISQL) THEN
            ISPQ = MULD2H(ISP1,ISQ1)
            IF (ISPQ .NE. JWOPSY) THEN
               READ (LUSUPM)
#if defined (VAR_SPLITFILES)
               IF (ITYP .EQ. 1) THEN
                  JBUF = JBUF + NPQRS
               ELSE
                  JBUF = JBUF + NPQRS*(1 + IRAT)
               END IF
#endif
               GO TO 100
            END IF
         END IF
C
         IF (ITYP .EQ. 1) THEN
            IREC1 = IREC1 + 1
            LENP  = IRAT*NPQRS
            IF (LENP .GT. LFREE) GO TO 810
#if defined (VAR_SPLITFILES)
            CALL READSI(LUSUPM,LENP,WORK,JBUF)
#else
            CALL READI(LUSUPM,LENP,WORK)
#endif
            DO 280 ISIM = 1,NSIM
               IOFF = 0
               DO 260 IP = NP1,NPL
                  ISP = ISAO(IP)
                  IF (IP .EQ. NPL) THEN
                     NQEND = NQL
                  ELSE
                     NQEND = IP
                  END IF
                  IF (IP .EQ. NP1) THEN
                     NQBEG = NQ1
                  ELSE
                     NQBEG = 1
                  END IF
                  DO 240 IQ = NQBEG, NQEND
                     ISQ = ISAO(IQ)
                     ISPQ = MULD2H(ISQ,ISP)
                     IPQ = IROW(IP) + IQ
                   IF (ISPQ .EQ. JWOPSY) THEN
                     DPQ = DMAT(IPQ,ISIM)
C                    note NRS = IPQ
                     IF (ABS(DPQ) .GT. THRZER)
     &               CALL DAXPY(IPQ,DPQ,WORK(IOFF+1),1,FMAT(1,ISIM),1)
                     FMAT( IPQ , ISIM ) = FMAT( IPQ , ISIM ) +
     &                  DDOT(IPQ,WORK(IOFF+1),1,DMAT(1,ISIM),1)
                   END IF
                     IOFF = IOFF + IPQ
  240             CONTINUE
                  NQ1 = 1
  260          CONTINUE
              IF (IOFF .NE. NPQRS .OR. NQL .GT. NPL) THEN
               WRITE (LUPRI,'(//5X,A,3(/5X,A,I12)/5X,A/5X,5I5,2I12)')
     &         'RDSUPP ERROR, inconsistent type 1 specifications',
     &         '   in record pair no.',IREC,
     &         '   Number of integrals from NP1,NQ1,NPL,NQL:',IOFF,
     &         '   Number of integrals from NPQRS          :',NPQRS,
     &         '   dump of ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS:',
     &         ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS
               CALL QTRACE(LUPRI)
               CALL QUIT(
     &            'RDSUPP ERROR, inconsistent type 1 spec. on LUSUPM')
              END IF
  280       CONTINUE
         ELSE IF (ITYP .EQ. 2) THEN
            IREC2 = IREC2 + 1
            LENP  = (IRAT+1)*NPQRS
            IF (LENP .GT. LFREE) GO TO 820
#if defined (VAR_SPLITFILES)
            CALL READSI(LUSUPM,LENP,WORK,JBUF)
#else
            CALL READI(LUSUPM,LENP,WORK)
#endif
            INDP  = IRAT*NPQRS
            IPQ   = IROW(NP1) + NQ1
            DO 380 ISIM = 1,NSIM
               DPQ = DMAT(IPQ,ISIM)
               IF (ABS(DPQ) .GT. THRZER) THEN
#include <ivdep.h>
                  DO 320 IRS = 1,NPQRS
                     FMAT( IWORK(INDP+IRS) , ISIM ) =
     &               FMAT( IWORK(INDP+IRS) , ISIM ) + DPQ * WORK(IRS)
  320             CONTINUE
               END IF
               FPQ = D0
               DO 340 IRS = 1,NPQRS
                  FPQ = FPQ + DMAT( IWORK(INDP+IRS) ,ISIM ) * WORK(IRS)
  340          CONTINUE
C Use first line below if P(pq,pq) has been divided by two
C otherwise use second line.
               FMAT( IPQ , ISIM ) = FMAT( IPQ , ISIM ) + FPQ
C              FMAT( IPQ , ISIM ) = FPQ
  380       CONTINUE
         ELSE
            WRITE (LUPRI,'(//5X,A,I5/5X,A,I12/5X,A/5X,5I5,2I12)')
     &         'RDSUPP ERROR, unrecognized type on LUSUPM:',ITYP,
     &         '   record pair no.',IREC,
     &         '   dump of ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS:',
     &         ITYP,NP1,NQ1,NPL,NQL,IBEF,NPQRS
            CALL QTRACE(LUPRI)
            CALL QUIT('RDSUPP ERROR, unrecognized type on AOSUPINT')
         END IF
      GO TO 100
C
  810 CONTINUE
         CALL ERRWRK('RDSUPP.TYP1',LENP,LFREE)
  820 CONTINUE
         CALL ERRWRK('RDSUPP.TYP2',LENP,LFREE)
  900 CONTINUE
      IF (IPRSTAT .GT. 5) THEN
         WRITE (LUERR,*) 'Total number of records :',IREC
         WRITE (LUERR,*) 'No. of records of type 1:',IREC1
         WRITE (LUERR,*) 'No. of records of type 2:',IREC2
      END IF
C
      CALL QEXIT('RDSUPP')
      RETURN
      END
C  /* Deck sirh1 */
      SUBROUTINE SIRH1(H1AO,WRK,LWRK)
C
C Copyright 29_nov-1990 Hans Joergen Aa. Jensen
C
#include <implicit.h>
      DIMENSION H1AO(*), WRK(LWRK)
C
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C  INFINP : NFIELD, EFIELD, LFIELD
C  INFORB : NNBAST,...
C
#include <maxorb.h>
#include <priunit.h>
#include <thrzer.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
      LOGICAL ANTSYM
C
      CALL QENTER('SIRH1')
      IF (P6FLAG(18)) WRITE (LUPRI,'(/A)')
     &   ' >>>>> One-electron Hamiltonian matrix from SIRH1 :'
      CALL RDONEL('ONEHAMIL',.TRUE.,H1AO,NNBAST)
C
C     Add field-dependent terms
C
      IF (NFIELD .LE. 0) GO TO 800
      KPRPAO = 1
      KH1AO  = KPRPAO + NNBASX
      LNEED  = KH1AO  + NNBAST
      IF (LNEED .GT. LWRK) CALL ERRWRK('SIRH1.RDPROP',LNEED,LWRK)
      DO 100 IFIELD = 1,NFIELD
         IF (P6FLAG(18)) WRITE (LUPRI,'(A,1P,G15.6,A,A8)')
     &      ' Field strength :',EFIELD(IFIELD),
     &      ' Field type : ',LFIELD(IFIELD)
      IF (EFIELD(IFIELD) .NE. D0) THEN
         CALL RDPROP(LFIELD(IFIELD),WRK,ANTSYM)
         IF (.NOT. ANTSYM) THEN
            CALL PKSYM1(WRK(KPRPAO),WRK(KH1AO),NBAS,NSYM,1)
            DPX = DASUM(NNBASX,WRK(KPRPAO),1)
            DPT = DASUM(NNBAST,WRK( KH1AO),1)
            IF (ABS(DPX - DPT) .GT. THRZER*DPX) THEN
               WRITE (LUPRI,'(/3A,/2A,2(/A,1P,D25.15))')
     &            ' FATAL ERROR in SIRH1: ',LFIELD(IFIELD),
     &            ' is not totally symmetric,',
     &            ' and finite field is only allowed with',
     &            ' totally symmetric operators.',
     &            ' Absolute sum of lower triangle of property matrix:',
     &            DPX,
     &            ' Absolute sum of the totally symmetric elements   :',
     &            DPT
               CALL QUIT(
     &         'ERROR: finite field with non tot.sym. operator')
            END IF
            CALL DAXPY(NNBAST,EFIELD(IFIELD),WRK(KH1AO),1,H1AO,1)
         ELSE
            WRITE (LUPRI,'(/3A/A)') ' FATAL ERROR in SIRH1: ',
     &   LFIELD(IFIELD),' is antisymmetric (i.e. imaginary)',
     &   ' and finite field with imaginary operators is not allowed.'
            CALL QUIT('ERROR: finite field with imaginary operator')
         END IF
      END IF
  100 CONTINUE
  800 IF (P6FLAG(18)) CALL OUTPKB(H1AO,NBAS,NSYM,1,LUPRI)
      CALL QEXIT('SIRH1')
      RETURN
      END
C  /* Deck rdprop */
      SUBROUTINE RDPROP(WORD,PRPAO,ANTSYM)
C
C Written 28-Nov-1990 HJAaJ
C
C Read AO property integrals with label WORD from LUPROP
C Output:
C  PRPAO  : AO property integrals
C  ANTSYM : true if integral matrix is antisymmetric
C           false if matrix is symmetric
C
#include <implicit.h>
#include <dummy.h>
      CHARACTER*8 WORD
      DIMENSION PRPAO(*)
      LOGICAL   ANTSYM
C
C -- common blocks:
C  INFORB : NNBASX, N2BASX
C
#include <priunit.h>
#include <inftap.h>
#include <inforb.h>
C
      CHARACTER*8 RTNLBL(2)
      LOGICAL OPEND, FNDLB2
C
      CALL QENTER('RDPROP')
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL GPOPEN(LUPROP,'AOPROPER',
     &            'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUPROP
      IF ( FNDLB2(WORD,RTNLBL,LUPROP)) THEN
         IF (RTNLBL(2).EQ.'SQUARE  ') THEN
            CALL READT(LUPROP,N2BASX,PRPAO)
            ANTSYM = .FALSE.
         ELSE IF (RTNLBL(2).EQ.'SYMMETRI') THEN
            ANTSYM = .FALSE.
            CALL READT(LUPROP,NNBASX,PRPAO)
         ELSE IF (RTNLBL(2).EQ.'ANTISYMM') THEN
            ANTSYM = .TRUE.
            CALL READT(LUPROP,NNBASX,PRPAO)
         ELSE
            WRITE (LUPRI,'(//3A/3(A,1X))')
     &         ' >>> RDPROP error: property "',WORD,
     &         '" on AOPROPER has no valid antisymmetry label.',
     &         '     Labels read on file:',RTNLBL(1),RTNLBL(2)
            CALL QTRACE(LUPRI)
            CALL QUIT('RDPROP error: No antisymmetry label on AOPROPER')
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> RDPROP error: property "',WORD,
     *      '" not found on AOPROPER.'
         CALL QTRACE(LUPRI)
         CALL QUIT('RDPROP error: property not found on AOPROPER.')
      END IF
      CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('RDPROP')
      RETURN
      END
C  /* Deck rd_sirifc */
      SUBROUTINE RD_SIRIFC(KEY,FOUND,AMAT,WRK,LWRK)
C
C Read specific information from Sirius interface file
C
C 10-Dec-1992 HJAaJ
C Revision 2000/05/24 HJAaJ: new options (read CMO and others)
C Revision 2011/05/02 HJAaJ: new keys FC, FCDIAG
C
C The following records have been written by WR_SIRIFC:
C
C    0) label LBSIFC ("SIR IPH ")
C    1) POTNUC,EMY,EACTIV,EMCSCF
C    2) NISHT,NASHT,...
C    3) CMO
C    4) CREF
C    5) DV
C    6) FOCK
C    7) PV
C    8) FC
C    9) FV
C   10) FCAC
C   11) H2AC
C   12) GORB
C   If (GUGA) then
C      13) label "CIDIAG1 "
C      14) CIDIAG1
C   end if
C   15) label "CIDIAG2 "
C   16) CIDIAG2
C   17) label "ORBDIAG "
C   18) ORBDIAG
C
C    *) label "SIRFLAGS"
C    *) FLAG(1:NFLAG)
C    *) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C   If (fields) then
C    *) label "EXTFIELD"
C    *) NFIELD
C    *) (EFIELD(I),I=1,NFIEL4)
C    *) (LFIELD(I),I=1,NFIEL4)
C   end if
C   If (solvent) then
C    *) label "SOLVINFO"
C    *) EPSOL,EPSTAT,RSOL(1:3),LSOLMX
C    *) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C    *) ERLM(LM,1), LM = 1,NLMSOL)
C    *) (TRLM(LM), LM = 1,NLMSOL)    where TRLM(i) = ERLM(i,2)
C    *) NSYM, NBAS
C    *) label "SOLVTMAT"
C    *) TMAT(1:NNORBX)
C   end if
C   If (CSF's) then
C    *) label "CREFDETS"
C    *) CREF in dets
C   end if
C    *) label 'TRCCINT"
C    *) NSYM, NORBT, ...
C    *) FCDIA, ISMO
C    *) CMO
C
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER*(*) KEY
      LOGICAL   FOUND
      DIMENSION AMAT(*), WRK(LWRK)
C
C Used from common blocks:
C   INFINP : NLMSOL
C   INFORB : NSYM,NCMOT,NNORBT,NNASHX
C   INFOPT : POTNUC, EMY
C   INFTAP : LUSIFC,FNSIFC,LBSIFC
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
C
      LOGICAL CLSIFC, FNDLAB
C
      CALL QENTER('RD_SIRIFC')
      IF (LUSIFC .GT. 0) THEN
         CLSIFC = .FALSE.
      ELSE
         CLSIFC = .TRUE.
         CALL GPOPEN(LUSIFC,FNSIFC,
     &      'UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      END IF
      REWIND LUSIFC
      IF (KEY .EQ. 'CMO') THEN
         IF (.NOT. FNDLAB(LBSIFC,LUSIFC)) THEN
            REWIND LUSIFC
            IF (.NOT. FNDLAB('CIRESPON',LUSIFC)) GO TO 8888
         ENDIF
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT(LUSIFC,NCMOT,AMAT)
      ELSE IF (KEY .EQ. 'DV') THEN
         IF (.NOT. FNDLAB(LBSIFC,LUSIFC)) GO TO 8888
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT(LUSIFC,NNASHX,AMAT)
      ELSE IF (KEY .EQ. 'PV') THEN
         IF (.NOT. FNDLAB(LBSIFC,LUSIFC)) GO TO 8888
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT(LUSIFC,NNASHX*NNASHX,AMAT)
      ELSE IF (KEY .EQ. 'FOCK') THEN
         IF (.NOT. FNDLAB(LBSIFC,LUSIFC)) GO TO 8888
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT(LUSIFC,N2ORBT,AMAT)
      ELSE IF (KEY .EQ. 'FC for MP2') THEN ! 3-May-2011 hjaaj, specifically for use in MP2FCK
         IF (.NOT. FNDLAB(LBSIFC,LUSIFC)) GO TO 8888
         ! rec no. 1
!        WRITE (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,
!    &               ISTATE,ISPIN,NACTEL,LSYM
         READ (LUSIFC) POTNUC_IFC, EMY ! we need EMY for SIRMP2
         IF (POTNUC .NE. POTNUC_IFC)
     &      CALL QUIT('POTNUC on SIRIFC and input are different')
         ! rec no. 2
!        WRITE (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,
!    &               NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,
!    &               NSYM,MULD2H, NRHF,NFRO,
!    &               NISH,NASH,NORB,NBAS,
!    &               NELMN1, NELMX1, NELMN3, NELMX3, MCTYPE,
!    &               NAS1, NAS2, NAS3
         READ (LUSIFC) NISHT_IFC,NASHT_IFC,idum,NORBT_IFC,NBAST_IFC,
     &               (idum,i=1,9),NSYM_IFC
         IF (NISHT .ne. NISHT_IFC .or. NASHT .ne. NASHT_IFC .or.
     &       NORBT .ne. NORBT_IFC .or. NBAST .ne. NBAST_IFC .or.
     &       NSYM  .ne. NSYM_IFC ) CALL
     &      QUIT('SIRIFC info not consistent with this calculation.')

         DO IREC = 3,7 ! skip records 3:7
            READ (LUSIFC)
         END DO
         CALL READT(LUSIFC,NNORBT,AMAT) ! FC matrix is in rec. no. 8
      ELSE IF (KEY .EQ. 'SOLTLM') THEN
         IF (.NOT. FNDLAB('SOLVINFO',LUSIFC)) GO TO 8888
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT(LUSIFC,NLMSOL,AMAT)
      ELSE IF (KEY .EQ. 'SOLVTMAT') THEN
         IF (.NOT. FNDLAB('SOLVTMAT',LUSIFC)) GO TO 8888
         CALL READT(LUSIFC,NNORBT,AMAT)
      ELSE IF (KEY .EQ. 'FCDIAG') THEN
         IF (.NOT. FNDLAB('TRCCINT ',LUSIFC)) GO TO 8888
         READ (LUSIFC) NSYM_SIRIFC
         IF (NSYM .NE. NSYM_SIRIFC) THEN
            CALL QUIT('RD_SIRIFC - key FCDIAG : NSYM inconsistent')
         END IF
         CALL READT(LUSIFC,NORBT,AMAT)
      ELSE
         WRITE(LUPRI,*) 'RD_SIRIFC error, invalid key word ',KEY
         CALL QUIT('RD_SIRIFC error, invalid key word ')
      END IF
      FOUND = .TRUE.
      GO TO 9999
C
 8888 FOUND = .FALSE.
 9999 IF (CLSIFC) CALL GPCLOSE (LUSIFC,'KEEP')
      CALL QEXIT('RD_SIRIFC')
      RETURN
      END
C  /* Deck rdcref */
      SUBROUTINE RDCREF(CREF)
C
C 26-Sep-1990 Hans Joergen Aa. Jensen
C Read reference CI vector from LUIT1
C
#include <implicit.h>
      DIMENSION CREF(*)
      PARAMETER ( D1 = 1.0D0 )
C
C Used from common blocks:
C   INFINP : ISTATE
C   INFVAR : NCONF
C   INFTAP : LUIT1
C   INFPRI : NWARN
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
C
      CHARACTER*8 RTNLBL(2), STAR8
      DATA STAR8/'********'/
C
      IF ( NCONF .GT. 1 ) THEN
         REWIND LUIT1
         CALL MOLLB2 ('STARTVEC',RTNLBL,LUIT1,LUPRI)
         IF (RTNLBL(1) .NE. STAR8) THEN
            READ (RTNLBL(1),'(2I4)') NRS,I
            IF (ABS(I) .NE. ISTATE) THEN
            IF (I .EQ. 0 .AND. RTNLBL(2).EQ.'CISAVE  ') THEN
               IF (NRS .LT. ISTATE) THEN
                  WRITE (LUPRI,'(//2(A,I4/))')
     &    ' ERROR from RDCREF:       ISTATE specified in input:',ISTATE,
     &    ' is greater than # of CI vectors on LUIT1 (CISAVE) :',NRS
                  CALL QUIT('RDCREF: ISTATE CI vector not available')
               END IF
               I = ISTATE
            ELSE
               NWARN = NWARN + 1
               WRITE (LUPRI,'(//2(A,I4/))')
     *      ' WARNING from RDCREF: ISTATE specified in input:',ISTATE,
     *      '                      ISTATE read from LUIT1   :',ABS(I)
            END IF
            END IF
         ELSE
C           this is an old LUIT1 file (before 5-Aug-1986)
            I = ISTATE
         END IF
         IF (I .GT. 0) THEN
            DO 110 I = 1,(ISTATE-1)
               READ (LUIT1)
  110       CONTINUE
C        else only reference vector saved
         END IF
         CALL READT (LUIT1,NCONF,CREF)
      ELSE
         CREF(1) = D1
      END IF
C
      RETURN
      END
C  /* Deck sirphp */
      SUBROUTINE SIRPHP(DIAGL,EACTVN,XNDXCI,FCAC,H2AC,WRK,LWRK)
C
C Oct 1990 hjaaj
C
C Purpose: set up diagonal and PHP matrix for SIRIUS calculations.
C          (Called by SIRNEO and SIRNR)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DIAGL(*),XNDXCI(*),FCAC(*),H2AC(*),WRK(*)
C
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFINP : LSYM, FLAG()
C   INFVAR : NCONF, NWOPT
C   INFDIM : LPHPMX, MAXPHP
C   INFTAP : LUIT2
C   INFPRI : IPRI6
C
#include <maxorb.h>
#include <mxcent.h>
Clf#include <pcmdef.h>
Clf#include <pcm.h>
#include <pcmlog.h>
#include <infinp.h>
#include <infvar.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C
      CHARACTER*8 LABIT2(3)
      DATA LABIT2/'CIDIAG2 ','ORBDIAG ','SOLVDIAG'/
C
      CALL QENTER('SIRPHP')
C
C     Get diagonal of Hessian / L-matrix
C     (divided by two for PHP routines)
C
      IF (NCONF .GT. 1) THEN
         REWIND LUIT2
C        Find label CIDIAG2, for diagonal of CI Ham. matrix.
C        Hessian diagonal is 2 * (CIDIAG2 - EACTVN).
         CALL MOLLAB(LABIT2(1),LUIT2,LUPRI)
         CALL READT(LUIT2,NCONF,DIAGL)
         IF (FLAG(16) .OR. PCM) THEN
C        IF (SOLVNT) THEN
            MAXPHP = 0
CPHPMAERKE  ... PHP not implemented for SOLVNT in this version.
            IF (LWRK .LT. NCONF) CALL ERRWRK('SIRPHP',NCONF,LWRK)
            CALL MOLLAB(LABIT2(3),LUIT2,LUPRI)
            CALL READT (LUIT2,NCONF,WRK)
            CALL DAXPY (NCONF,DP5,WRK,1,DIAGL,1)
         END IF
         DO 100 I = 1,NCONF
  100       DIAGL(I) = DIAGL(I) - EACTVN
      ELSE
         DIAGL(1) = D0
      END IF
C
      IF (NWOPT .GT. 0) THEN
         REWIND LUIT2
C        Find label ORBDIAG, for diagonal of orbital Hessian.
         CALL MOLLAB(LABIT2(2),LUIT2,LUPRI)
         CALL READT(LUIT2,NWOPT,DIAGL(1+NCONF))
      END IF
      IF (MAXPHP .GT. 0 .AND. NCONF .GT. 1) THEN
         ESHIFT = -EACTVN
         IPWAY  = 1
         IPRPHP = MAX(0,IPRI6-4)
         CALL PHPGET(LSYM,NCONF,XNDXCI,FCAC,H2AC,DIAGL,ESHIFT,
     *               IPWAY,IPRPHP,WRK,LWRK)
      END IF
      CALL SRWPHP(DIAGL,LPHPMX,.TRUE.)
C     CALL SRWPHP(PHPVEC,LPHPMX,PHPWRT)
      CALL QEXIT('SIRPHP')
      RETURN
      END
C  /* Deck srwphp */
      SUBROUTINE SRWPHP(PHPVEC,LPHPMX,PHPWRT)
C
C 27-Oct-1990 Hans Joergen Aa. Jensen
C
C Sirius routine for Reading or Writing PHP information.
C
C Hessian DIAGONAL + PHP INFORMATION AS LAST RECORD ON LUIT2
C
C  PHPWRT = .true.  write information
C         = .false. read information
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION PHPVEC(*)
      LOGICAL PHPWRT
C
C Information used from common blocks;
C
C   INFTAP : LUIT2
C
#include <inftap.h>
C
      LOGICAL FNDLAB
      CHARACTER*8 LAB123(3),LBPHP,LBEOD
      DATA LAB123(1),LBPHP,LBEOD/'********','PHPINFO ','EODATA  '/
C
      REWIND LUIT2
C
C    read or write diagonal php information
C
      IF (PHPWRT) THEN
         CALL GETDAT(LAB123(2),LAB123(3))
         IF( FNDLAB(LBPHP,LUIT2) ) THEN
            BACKSPACE LUIT2
         ELSE
            REWIND LUIT2
            IF( FNDLAB(LBEOD,LUIT2) ) BACKSPACE LUIT2
         END IF
         WRITE (LUIT2) LAB123,LBPHP
         LPHPM4 = MAX(4,LPHPMX)
         CALL WRITT(LUIT2,LPHPM4,PHPVEC)
         WRITE (LUIT2) LAB123,LBEOD
      ELSE
         IF (.NOT.FNDLAB(LBPHP,LUIT2)) THEN
            CALL QTRACE(LUPRI)
            CALL QUIT('ERROR SRWPHP: label PHPINFO not found on LUIT2')
         END IF
         CALL READT(LUIT2,LPHPMX,PHPVEC)
      END IF
C
C *** END OF SRWPHP
C
      RETURN
      END
C  /* Deck upkcmo */
      SUBROUTINE UPKCMO(CMO,UCMO)
C
C  Copyright 6-May-1990 Hans Joergen Aa. Jensen
C
C  Unpack CMO from symmetry packed form to  UCMO
C
#include <implicit.h>
      DIMENSION CMO(*), UCMO(NBAST,NORBT)
C
C Used from common blocks:
C  INFORB : NSYM,NBAST,NORBT,...
C
#include <inforb.h>
C
      CALL DZERO(UCMO,NBAST*NORBT)
      DO 400 ISYM = 1,NSYM
         IF (NORB(ISYM) .GT. 0)
     &   CALL MCOPY(NBAS(ISYM),NORB(ISYM),CMO(ICMO(ISYM)+1),NBAS(ISYM),
     &              UCMO(IBAS(ISYM)+1,IORB(ISYM)+1),NBAST)
C        CALL MCOPY(NROWA,NCOLA,A,NRDIMA,B,NRDIMB)
  400 CONTINUE
      RETURN
      END
