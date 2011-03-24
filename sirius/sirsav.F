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
C  Principal author of this file: Hans Jørgen Aa. Jensen
C
C  /* Deck sirsav */
      SUBROUTINE SIRSAV (KEYWRD,CMO,IBNDX,REDL,EVEC,XKAP,INDXCI,
     *                   WRK,LFREE)
C
C  Written by Hans Agren 4-Apr-84
C  Revisions
C   14-May-1984 hjaaj
C   19-Nov-1984 hjaaj (treat CI (NWOPT=0) as special case)
C    6-Dec-1984 hjaaj (cleaned upd writing of CI-vec to LUIT1)
C    7-Jan-1985 hjaaj (write RESTART label on LUIT1 and info)
C    7-Mar-1985 hjaaj (included IBNDX)
C   28-Mar-1985 hjaaj (ICTL parameter)
C   19-Jun-1985 ha    (ICTL=5, merge of CISAVE)
C   14-jan-1986 hjaaj (ICTL=6 N-R save, SIRSAV restructured)
C    3-aug-1987 ha    (save orbital trial vectors for ICTL=1, NEOSAVE)
C
C Purpose: To write present orbitals, the updated CI-vectors and
C          the reduced L-matrix on unit LUIT1 (To be used for
C          a backstep or for another start)
C
C MOTECC-90: The purpose of this module, SIRSAV, and the algorithms used
C            are described in Chapter 8 Section C.2 of MOTECC-90
C            "Step Calculation and Step comntrol"
C
C Input:
C  KEYWRD; transformed to:
C  ICTL; =1 save MC optimization information, when NEO used.
C        =2 save geometry walk information for MC
C           (called from ABACUS, no reduced L)
C        =3 save MC optimization information, when update used.
C        =4 save MC information, after SIRCAN
C        =5 save CI information, after CICTL
C        =6 save MC optimization information, when N-R used.
C
#include <implicit.h>
C
      CHARACTER*(*) KEYWRD
      DIMENSION CMO(*),IBNDX(*),EVEC(*),REDL(*),XKAP(*),INDXCI(*)
      DIMENSION WRK(LFREE)
#include <ibndxdef.h>
C
      PARAMETER (THRTT = 1.0D-4, D0 = 0.0D0, D1 = 1.0D0)
#include <dummy.h>
C
C Used from common blocks:
C   INFINP : ISTATE,FLAG(*),?
C   INFVAR : NCONF,NWOPT,NVAR
C   INFORB : NCMOT,N2BAST,NASHT,?
C   INFOPT : EMCOLD,DEPRED,BETA,GAMMA,RTRUST,ITMAC,NREDL,?
C   INFDIM : NORBMA,NBASMA,IVEX4,INTOT, MAXRL
C   INFTAP : LUIT1,?
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infopt.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL DOCNO, DOXKAP, FNDLAB
      CHARACTER*6 CNOLBL
      CHARACTER*8 TABLE(9), STAKEY, RSTKEY, BLANK8, LAB123(3)
      DATA TABLE /'EODATA  ', 'OLDORB  ','STARTVEC', 'LREDUCED',
     *            'ROTKAPPA', 'NEWORB  ','RESTART ', 'GEOWALK ',
     &            'NATOCCAC'/
      DATA BLANK8/'        '/
      DATA LAB123/'********','********','********'/
C
      CALL QENTER('SIRSAV')
      CALL GETDAT(LAB123(2),LAB123(3))
C     ... place date in LAB123(2) and time in LAB123(3)
C
C
      IF (ISTATE .LE. 0 .AND. KEYWRD .NE. 'CISAVE') THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//2A/A,I4)')
     *      ' WARNING from SIRSAV with keyword ',KEYWRD,
     *      ' ------- ISTATE set to 1, was',ISTATE
         WRITE (LUERR,'(//2A/A,I4)')
     *      ' WARNING from SIRSAV with keyword ',KEYWRD,
     *      ' ------- ISTATE set to 1, one',ISTATE
         ISTATE = 1
      END IF
      DOCNO = FLAG(15)
      IF (KEYWRD .EQ. 'NEOSAVE') THEN
         ICTL = 1
         NRS  = MIN (NROOTS,NREDL,NCONF)
         JSTA = ISTATE
      ELSE IF (KEYWRD .EQ. 'GEOSAVE') THEN
         ICTL = 2
         NRS  = 1
         JSTA = -ISTATE
      ELSE IF (KEYWRD .EQ. 'UPDSAVE') THEN
         ICTL = 3
         NRS  = 1
         JSTA = -ISTATE
      ELSE IF (KEYWRD .EQ. 'CANSAVE') THEN
         ICTL = 4
         NRS  = 1
         JSTA = -ISTATE
      ELSE IF (KEYWRD .EQ. 'CISAVE') THEN
         ICTL = 5
         NRS  = MIN (NROOTS,NREDL,NCONF)
         JSTA = ISTATE
         IF (JSTA .EQ. 0 .AND. NROOTS .GT. 1 .AND. DOCNO) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(3(/A))')
     &      ' WARNING from SIRSAV with keyword CISAVE',
     &      '    Canonical/natural orbital transformation abandoned',
     &      '    because # CI ROOTS .gt. 1 and no ".STATE" selected'
            DOCNO = .FALSE.
         END IF
      ELSE IF (KEYWRD .EQ. 'NRSAVE') THEN
         ICTL = 6
         NRS  = 1
         JSTA = -ISTATE
      ELSE
         WRITE (LUPRI,'(//2A)') ' >>> ERROR (SIRSAV), KEYWRD = ',KEYWRD
         WRITE (LUERR,'(//2A)') ' >>> ERROR (SIRSAV), KEYWRD = ',KEYWRD
         CALL QTRACE(LUERR)
         CALL QUIT('*** ERROR (SIRSAV) KEYWRD illegal.')
      END IF
C
      RSTKEY = BLANK8
      LENKEY = MIN(8,LEN(KEYWRD))
      RSTKEY(1:LENKEY) = KEYWRD(1:LENKEY)
      WRITE (STAKEY,'(2I4)') NRS,JSTA
      NC4 = MAX(4,NCONF)
      NW4 = MAX(4,NWOPT)
      NV4 = MAX(4,NVAR)
      NCMOT4 = MAX(4,NCMOT)
C
C Section 1: Save old molecular orbitals
C
      CALL NEWIT1
      WRITE (LUIT1) LAB123(1),LAB123(2),RSTKEY, TABLE(2)
      CALL WRITT(LUIT1,NCMOT4,CMO)
C
C Section 2: (A) Save restart information on LUIT1.
C                Label: RESTART
C            (B) Write reduced L-matrix on LUIT1.
C                Label: LREDUCED
      IF (ICTL .EQ. 1 .OR. ICTL .EQ. 6) THEN
C     ... if (neosave or nrsave) then
         WRITE (LUIT1) LAB123(1),LAB123(2),RSTKEY,TABLE(7)
         WRITE (LUIT1) EMCOLD,DEPRED,BETA,GAMMA,RTRUST,ITMAC,
     *                 NREDL,(IBNDX(I),I=1,NREDL)
         WRITE (LUIT1) EMY,EACTIV,EMCOLD,EMCSCF,
     *                 FLAG(51),FLAG(52),FLAG(53),FLAG(39)
         WRITE (LUIT1) LAB123,TABLE(4)
         NNRED4 = MAX( 4 , (NREDL*(NREDL+1)/2) )
         IF (ICTL .EQ. 1) THEN
            CALL WRITT(LUIT1,NNRED4,REDL(1))
         ELSE
Chj-aug99:  NR; we must therefore create REDL matrix in WRK(KREDL)
C           so backup to NEO and use in UPDATE is possible
C           REDG = REDL(1:MAXRL)
C           REDH = REDL(MAXRL+1:)
C
            KREDL = 1
            KEND  = KREDL + NREDL*(NREDL+1)/2
            IF (KEND .GT. LFREE) 
     &         CALL ERRWRK('SIRSAV for REDL',KEND,LFREE)
            WRK(KREDL) = D0
Chj:        redl(1,1) = 0
            NREDH = NREDL - 1
            IREDH = MAXRL
            IREDL = KREDL
            DO I = 1, NREDH
               WRK(IREDL + 1) = REDL(I)
Chj:           redl(1+i,1) = redg(i)
               DO J = 1, I
                  WRK(IREDL + 1 + J) = REDL(IREDH + J)
Chj:              redl(1+i,j) = redh(i,j)
               END DO
               IREDH = IREDH + I
               IREDL = IREDL + 1 + I
            END DO
            CALL WRITT(LUIT1,NNRED4,WRK(KREDL))
         END IF
C
         IF (NWOPT.GT.0) THEN
            WRITE (LUIT1) LAB123,TABLE(5)
            CALL WRITT(LUIT1,NW4,XKAP)
         END IF
      ELSE IF (ICTL.EQ.2) THEN
         WRITE (LUIT1) LAB123,TABLE(8)
         WRITE (LUIT1) EMCOLD,DEPRED,REJWMI,REJWMA
      ELSE IF (ICTL.EQ.3) THEN
         WRITE (LUIT1) LAB123(1),LAB123(2),RSTKEY,TABLE(7)
         WRITE (LUIT1) EMCOLD,DUMMY,DUMMY,DUMMY
         IF (NWOPT.GT.0) THEN
            WRITE (LUIT1) LAB123,TABLE(5)
            CALL WRITT(LUIT1,NW4,XKAP)
         END IF
      END IF
C
C Section 3: Save trial vectors and new (rotated) orbitals
C
      GO TO (100,200,300,400,500,600), ICTL
  100 CONTINUE
C
C       "NEOSAVE"
C
C       If save for NEO MC optimization, then:
C
C       Calculate new starting vectors from EVEC and BVECS.
C       (calculate two more than we use so we easily can
C       increment NSIM by one or two, if desired)
C       Calculate XKAP matrix (returned to calling program)
C
C     OUTPUT:
C
      JEVEC  = (ISTATE-1)*NREDL
      FAC    = GAMMA / ABS( BETA*EVEC(JEVEC+1) )
      STPLEN = FAC * SQRT(D1 - EVEC(JEVEC+1) * EVEC(JEVEC+1))
      STPC = D0
      STPO = D0
      DO 1010 K = 2,NREDL
         IF (IBNDX(K).EQ.JBCNDX) THEN
            STPC = STPC + EVEC(JEVEC+K)*EVEC(JEVEC+K)
         ELSE
            STPO = STPO + EVEC(JEVEC+K)*EVEC(JEVEC+K)
         END IF
 1010 CONTINUE
      STPC = FAC * SQRT(STPC)
      STPO = FAC * SQRT(STPO)
C
      IF (P4FLAG(10)) THEN
         WRITE (LUW4,'(//3(/A,F15.10),/2(/A,2F15.10),/)')
     *      ' *** Step length  :',STPLEN,
     *      '  -  CI step      :',STPC,
     *      '  -  Orbital step :',STPO,
     *      ' *** Gamma, level shift :',GAMMA,SHFLVL,
     *      '  -  Beta,  "a1"        :',BETA,EVEC(JEVEC+1)
      END IF
C
C     Prepare for 2:
C
C     NRS  = number of CI vectors to be written to LUIT1
      NSIM = (LFREE/NVAR)-1
      NSIM = MIN(NSIM,NRS)
      IF (NSIM .LE. 0) CALL ERRWRK('SIRSAV (for NSIM = 1)',2*NVAR,LFREE)
      IF (NWOPT.GT.0) THEN
         CALL DZERO(XKAP,NWOPT)
         FAC0 = (GAMMA / BETA) / EVEC(JEVEC+1)
      ELSE
         FAC0 = D0
C        ... to avoid compiler messages
      END IF
      DOXKAP = .TRUE.
      KBVEC  = 1
      KY     = KBVEC + MAX(NCONF,NWOPT)
C
      WRITE (LUIT1) LAB123(1), STAKEY, RSTKEY, TABLE(3)
C
C     Now get them:
C
      JOFF = 0
 1100 CONTINUE
      REWIND LUIT3
      CALL DZERO(WRK(KY),NSIM*NVAR)
      DO 1200 K = 1,NREDL
         IF (IBNDX(K).EQ.JBCNDX) THEN
            CALL READT(LUIT3,NCONF,WRK(KBVEC))
            JYVEC = KY
            LYVEC = NCONF
         ELSE IF (IBNDX(K).EQ.JBONDX) THEN
            CALL READT(LUIT3,NWOPT,WRK(KBVEC))
            JYVEC = KY + NCONF
            LYVEC = NWOPT
CHJ- husk at fjerne dette "DOXKAP" senere ...
            IF (DOXKAP) THEN
               JEVEC = (ISTATE-1)*NREDL + K
               FAC = FAC0 * EVEC(JEVEC)
               CALL DAXPY(NWOPT,FAC,WRK(KBVEC),1,XKAP,1)
            END IF
         ELSE
            WRITE (LUERR,'(/A,2I5)')
     *         ' SIRSAV: illegal IBNDX(K); K, IBNDX(K) =',K,IBNDX(K)
            WRITE (LUERR,'(/A,/,(5I5,5X,5I5))')
     *         ' List of all IBNDX values:',(IBNDX(I),I=1,NREDL)
            CALL QTRACE(LUERR)
            CALL QUIT('SIRSAV: illegal IBNDX value')
         END IF
         DO 1300 J = 1,NSIM
            JEVEC = (JOFF+J-1)*NREDL + K
            FAC = EVEC(JEVEC)
            IF (K.EQ.1) THEN
C
C              special treatment for EVEC(1,ISTATE), i.e. a0,
C              instead of multiplying EVEC(2:NREDL,ISTATE) by
C              gamma / (beta*a0) we multiply EVEC(1,ISTATE) by
C              (beta*a0) / gamma; this is equivalent because
C              the vector is normalized below.
C
               IF (NWOPT.GT.0 .AND. (JOFF+J).EQ.ISTATE) THEN
                  FAC = FAC * (BETA/GAMMA)
               END IF
            END IF
            CALL DAXPY(LYVEC,FAC,WRK(KBVEC),1,WRK(JYVEC),1)
            JYVEC = JYVEC + NVAR
 1300    CONTINUE
 1200 CONTINUE
      DOXKAP = .FALSE.
C
C     We now have eigenvectors no. JOFF+1 to JOFF+NSIM.
C
      DO 1400 J = 1,NSIM
         KYJ = KY+(J-1)*NVAR
C
C        Make linear correction for change in orbital part
C        because of orbital transformation by XKAP.
C        ( /istate(n+1)> = exp(-xkap) /istate(n)'>,
C          where the prime denotes changed CI coefficients,
C          /j(n+1)> = exp(-x(j)) /j(n)'>, where j .ne. istate. Combining
C          we get /j(n+1)> = exp(-x(j)) exp(+xkap) /j(n+1)">,
C          to first order = exp(- (x(j)-xkap)) /j(n)>, where
C          /j(n+1)"> = exp(-xkap) /j(n)>. The modified vector
C          x(j) - xkap is thus a guess for the orbital part of
C          NEO eigenvector no. j in next iteration).
C
         DN = DNORM2(NCONF,WRK(KYJ),1)
         IF (DN.NE.D1 .AND. DN.NE.D0) THEN
            IF (DN .LE. THRTT) THEN
               DN = D1 / DN
               CALL DSCAL(NCONF,DN,WRK(KYJ),1)
               DN = DNORM2(NCONF,WRK(KYJ),1)
            END IF
            DN = D1 / DN
            CALL DSCAL(NCONF,DN,WRK(KYJ),1)
         END IF
         IF (.NOT. DOCNO) THEN
C        ... if we transform to Fock type/canonical orbitals (FLAG(15)),
C            then this transformation would further change the
C            x(j) - xkap vector in an unknown way. Thus we skip the
C            orbital part of the eigenvectors if flag(15).
C
            DO 1405 I = 1,NWOPT
               WRK(KYJ + NCONF + I) = WRK(KYJ + NCONF + I) - XKAP(I)
 1405       CONTINUE
            CALL WRITT(LUIT1,NV4,WRK(KYJ))
         ELSE
            CALL WRITT(LUIT1,NC4,WRK(KYJ))
         END IF
C
 1400 CONTINUE
C
      JOFF = JOFF + NSIM
      NSIM = MIN(NSIM,NRS-JOFF)
      IF (NSIM.GT.0) GO TO 1100
C     ^------------------------
C
      GO TO 3000
C ------------------------------------------------------
C     "GEOSAVE" (200)
C     "UPDSAVE" (300)
C     "CANSAVE" (400)
C
  200 CONTINUE
  300 CONTINUE
  400 CONTINUE
C
C        save new CREF vector and go to next section:
C        (new CREF vector is in WRK)
C
         WRITE (LUIT1) LAB123(1), STAKEY, RSTKEY, TABLE(3)
         CALL WRITT(LUIT1,NC4,WRK)
      IF (KEYWRD .NE. 'CANSAVE') THEN
         GO TO 3000
      ELSE
C
C        "CANSAVE" only
C
C        Save restart information on LUIT1.
C        Label: RESTART
C
         WRITE (LUIT1) LAB123(1),LAB123(2),RSTKEY,TABLE(7)
         WRITE (LUIT1) EMCSCF,DUMMY,DUMMY,DUMMY
C
C        save orbitals with label NEWORB
C
         LAB123(3) = 'CANORB  '
         WRITE (LUIT1) LAB123,TABLE(6)
         CALL WRITT(LUIT1,NCMOT4,CMO)
         GO TO 9000
      END IF
C        ----------------------------------------------------
  500 CONTINUE
C
C        ******************************
C        ***** <<<<< CISAVE >>>>> *****
C        ******************************
C
C
C        Calculate CI vectors from EVEC And BVECS.
C        save them with label STARTVEC
C
C        (calculate the NROOTS specified in input)
C        NRS   = number of CI vectors to be written to LUIT1
C
         KBVEC = 1
         KY    = KBVEC + NCONF
         NSIM  = (LFREE/NCONF)-1
         NSIM  = MIN(NSIM,NRS)
         JOFF  = 0
C
         WRITE (LUIT1) LAB123(1), STAKEY, RSTKEY, TABLE(3)
 5100    CONTINUE
         REWIND LUIT3
         CALL DZERO(WRK(KY),NSIM*NCONF)
         DO 5200 K=1,NREDL
           CALL READT(LUIT3,NCONF,WRK(KBVEC))
           DO 5300 J=1,NSIM
             JEVEC = (JOFF+J-1)*NREDL + K
             CALL DAXPY(NCONF,EVEC(JEVEC),WRK(KBVEC),1,
     *                  WRK(KY+(J-1)*NCONF),1)
 5300      CONTINUE
 5200    CONTINUE
C
C        We now have eigenvectors no. JOFF+1 to JOFF+NSIM.
         DO 5400 J = 1,NSIM
            DN = DNRM2(NCONF,WRK(KY+(J-1)*NCONF),1)
            IF (DN.NE.D1) THEN
               DN = D1 / DN
               CALL DSCAL(NCONF,DN,WRK(KY+(J-1)*NCONF),1)
            END IF
            CALL WRITT(LUIT1,NC4,WRK(KY+(J-1)*NCONF))
 5400    CONTINUE
         JOFF = JOFF + NSIM
         NSIM = MIN(NSIM,NRS-JOFF)
         IF (NSIM.GT.0) GO TO 5100
C
C
C        save orbitals with label NEWORB
C
         LAB123(3) = ' CIORB  '
         WRITE (LUIT1) LAB123,TABLE(6)
         CALL WRITT(LUIT1,NCMOT4,CMO)
C
C        Go to 4000 to transform to CNO orbitals, if requested.
C
         GO TO 4000
C        ----------------------------------------------------
  600 CONTINUE
C
C        ******************************
C        ***** <<<<< NRSAVE >>>>> *****
C        ******************************
C
C
C
C      If save for NR MC optimization, then:
C
C       Section 6a: Calculate new reference vector from EVEC and BVECS.
C       Section 6b: Calculate XKAP matrix (returned to calling program)
C
      WRITE (LUIT1) LAB123(1), STAKEY, RSTKEY, TABLE(3)
C
      NREDH  = NREDL - 1
      STPLEN = DNORM2(NREDH,EVEC,1)
      STPC = D0
      STPO = D0
      DO 6010 K = 1,NREDH
         IF (IBNDX(1+K).EQ.JBCNDX) THEN
            STPC = STPC + EVEC(K)*EVEC(K)
         ELSE
            STPO = STPO + EVEC(K)*EVEC(K)
         END IF
 6010 CONTINUE
      STPC = SQRT(STPC)
      STPO = SQRT(STPO)
C
      IF (P4FLAG(10)) THEN
         WRITE (LUW4,'(//3(/A,F15.10))')
     *      ' *** Step length  :',STPLEN,
     *      '  -  CI step      :',STPC,
     *      '  -  Orbital step :',STPO
      END IF
C
C     Prepare for 6b:
C
      KBVEC = 1
      KY    = KBVEC + MAX(NCONF,NWOPT)
C
C     6a: Now get the new reference vector
C
      REWIND LUIT3
      IF (NCONF .GT. 1) THEN
         CALL READT(LUIT3,NCONF,WRK(KY))
      ELSE
         READ (LUIT3)
         WRK(KY) = D1
      END IF
      IF (NWOPT.GT.0) THEN
         CALL DZERO(XKAP,NWOPT)
      END IF
      DO 6200 K = 1,NREDH
         IF (IBNDX(1+K).EQ.JBCNDX) THEN
            CALL READT(LUIT3,NCONF,WRK(KBVEC))
            FAC = EVEC(K)
            CALL DAXPY(NCONF,FAC,WRK(KBVEC),1,WRK(KY),1)
         ELSE IF (IBNDX(1+K).EQ.JBONDX) THEN
            CALL READT(LUIT3,NWOPT,WRK(KBVEC))
            FAC = EVEC(K)
            CALL DAXPY(NWOPT,FAC,WRK(KBVEC),1,XKAP,1)
         ELSE
            WRITE (LUERR,'(/A,2I5)')
     *         ' NRSAVE: illegal IBNDX(K); K, IBNDX(K) =',K,IBNDX(K)
            WRITE (LUERR,'(/A,/,(5I5,5X,5I5))')
     *         ' List of all IBNDX values:',(IBNDX(I),I=1,NREDL)
            CALL QTRACE(LUERR)
            CALL QUIT('NRSAVE: illegal IBNDX value')
         END IF
 6200 CONTINUE
C
      IF (NCONF .GT. 1) THEN
         DN = DNORM2(NCONF,WRK(KY),1)
         IF (DN.NE.D1) THEN
            IF (DN .LE. THRTT) THEN
               DN = D1 / DN
               CALL DSCAL(NCONF,DN,WRK(KY),1)
               DN = DNORM2(NCONF,WRK(KY),1)
            END IF
            DN = D1 / DN
            CALL DSCAL(NCONF,DN,WRK(KY),1)
         END IF
      END IF
      CALL WRITT(LUIT1,NC4,WRK(KY))
C
      GO TO 3000
C ------------------------------------------------------
C
C Section 3:     Rotate orbitals and write the new rotated orbitals
C                on LUIT1 with label NEWORB.
C
 3000 CONTINUE
      CALL MO_ROTATE(XKAP,CMO,WRK,LFREE)
      IF (NCONF .GT. 1) THEN
         LAB123(3) = '(MCOPT) '
      ELSE
         LAB123(3) = '(HFOPT) '
      END IF
      WRITE (LUIT1) LAB123,TABLE(6)
      CALL WRITT(LUIT1,NCMOT,CMO)
      REWIND LUIT1
C
C
C Section 4:     Transform to canonical/natural orbitals, if requested
C
 4000 IF (DOCNO) THEN
         KCVECS = 1
         KAOCC  = KCVECS + NRS*NCONF
         KW31   = KAOCC  + NASHT
         LW31   = LFREE  - KW31
         IF (LW31 .LE. 0) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//A/A/A,I8,A,I8/)')
     *      ' WARNING SIRSAV, insufficient space for transformation',
     *      ' ==============  to natural orbitals, nothing done.',
     *      '                 available:',LFREE,' need more than:',LW31
            GO TO 9000
         END IF
         IF (NCONF .GT. 1) THEN
            REWIND (LUIT1)
C           search for "STARTVEC"
            CALL MOLLAB(TABLE(3),LUIT1,lupri)
            DO 9010 I = 1,NRS
               JCVECS = KCVECS + (I-1)*NCONF
               CALL READT(LUIT1,NCONF,WRK(JCVECS))
 9010       CONTINUE
         ELSE
C        ... case: one configuration
            WRK(KCVECS) = D1
         END IF
         IF (FLAG(48)) THEN
            CNOLBL = 'ONLYFD'
         ELSE IF (FLAG(46)) THEN
            CNOLBL = 'ONLYNO'
         ELSE
            CNOLBL = 'FD+NO '
         END IF
         IF (JSTA .LE. 0) THEN
            ICREF = 1
         ELSE
            ICREF = ISTATE
         END IF
         CALL SIRCNO(CNOLBL,NRS,ICREF,WRK(KCVECS),CMO,
     *               WRK(KAOCC),INDXCI,WRK(KW31),1,LW31)
C        CALL SIRCNO(KEYCNO,NCVECS,ICREF,CVECS,CMO,AOCC,
C    *               INDXCI,WRK,KFREE,LFREE)
C
         REWIND (LUIT1)
         CALL MOLLAB(TABLE(3),LUIT1,lupri)
C        search for "STARTVEC"
         RSTKEY = 'CNOSAVE '
         BACKSPACE LUIT1
C        write startvec with new vectors
         WRITE (LUIT1) LAB123(1), STAKEY, RSTKEY, TABLE(3)
         DO 9022 I = 1,NRS
            JCVECS = KCVECS + (I-1)*NCONF
            CALL WRITT(LUIT1,NC4,WRK(JCVECS))
 9022    CONTINUE
C        write the natural orbitals to LUIT1 label NEWORB
         CALL GETDAT(LAB123(2),LAB123(3))
         LAB123(3) = '(CNOORB)'
         WRITE (LUIT1) LAB123,TABLE(6)
         CALL WRITT(LUIT1,NCMOT4,CMO)
         IF (NASHT.EQ.1 .OR.
     &      (CNOLBL.NE.'ONLYFD' .AND. NASHT.GT.0)) THEN
            IF (NASHT.EQ.1) WRK(KAOCC) = NACTEL
            WRITE (LUIT1) LAB123,TABLE(9)
            NAT4 = MAX(NASHT,4)
            CALL WRITT(LUIT1,NAT4,WRK(KAOCC))
         END IF
         IF (.NOT.FLAG(34)) FLAG(14) = .FALSE.
C        if (int.transf. needed in optimization) then new transf. needed
      END IF
C
C *** End of subroutine SIRSAV
C     Empty buffers to LUIT1 by rewinding
C     (so we do not loose the restart information if for instance
C      the system crashes)
C
 9000 CONTINUE
      REWIND LUIT1
      IF ( .NOT.FNDLAB(TABLE(1),LUIT1) ) THEN
         CALL GETDAT(LAB123(2),LAB123(3))
         WRITE (LUIT1) LAB123,TABLE(1)
      END IF
      REWIND LUIT1
C
      CALL QEXIT('SIRSAV')
      RETURN
C     end of SIRSAV
      END
C  /* Deck MO_rotate */
      SUBROUTINE MO_ROTATE(XKAP,CMO,WRK,LFRSAV)
C
C Written by Hans Agren 6-Apr-1984
C Revisions:
C   18-May-1984 hjaaj
C    9-Nov-1984 hjaaj (if kaver.gt.0, average XKAP)
C    7-Jan-1985 hjaaj (write NEWORB to LUIT1)
C   31-Oct-1989 hjaaj (use matrix multiply, UXKAP not symmetry blocked)
C   12-Sep-1990 hjaaj (do not write NEWORB to LUIT1; print XKAP section)
C
C Purpose: To construct new set of orbitals from the Kappa matrix
C
C Input : XKAP the packed Kappa matrix
C         CMO old orbitals
C
C Output: CMO new orbitals
C
C Scratch:WRK
C
#include <implicit.h>
      DIMENSION XKAP(*),CMO(*), WRK(*)
C
C local parameters:
C
      PARAMETER ( D0=0.D0, DP5=0.5D0,  D1=1.D0, D1P5 = 1.5D0 )
      PARAMETER ( DIVERG = 2.0D0 , THREQL = 1.D-12 , THNORM = 1.D-04 )
C
C Used from common blocks:
C   INFINP : ??
C   INFVAR : NWOPT,JWOPSY,JWOP(2,*)
C   INFORB : NSYM,N2ORBT,NCMOT,...
C   INFDIM : NORBMA
C   INFPRI : P6FLAG()
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
C
C
      CALL QENTER('MO_ROTATE')
      IF (NWOPT.EQ.0) GO TO 9999
      IF (JWOPSY .NE. 1) THEN
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR, MO_ROTATE called with JWOPSY.ne.1')
      END IF
C
C Average rotation matrix.
C (AVERAG returns immediately if no averageing requested)
C
      CALL AVERAG(XKAP,NWOPT,1)
C
C Print XKAP, if so requested.
C
      IF (P6FLAG(29) .AND. NWOPT .GT. 0) THEN
         WRITE (LUPRI,'(/A/A)') ' MO_ROTATE: Orbital rotation vector',
     *                         ' -------------------------------'
         IF (IPRI6.GT.100) THEN
            PRFAC = 0.0D0
         ELSE
            PRFAC = 0.1D0
         END IF
         CALL PRKAP (NWOPT,XKAP,PRFAC,LUPRI)
      END IF
C
      KUXKAP = 1
      KWRK0  = KUXKAP + N2ORBX
      KNEED  = KWRK0  + 3*NORBMA*NBASMA
      IF (KNEED .GT. LFRSAV) CALL ERRWRK('MO_ROTATE_1',KNEED,LFRSAV)
C
C Unpack kappa vector and add unit matrix (giving 1+X)
C
      CALL UPKWOP(NWOPT,JWOP,XKAP,WRK(KUXKAP))
      DO 50 I = 1,NORBT
         WRK(KUXKAP-1 + (I-1)*NORBT + I) = D1
   50 CONTINUE
C
C Section 1 **********************
C Start symmetry loop
C
      DO 500 ISYM = 1,NSYM
         NORBI  = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 500
         IORBI  = IORB(ISYM)
         N2ORBI = N2ORB(ISYM)
         NBASI  = NBAS(ISYM)
         ICMOI  = ICMO(ISYM)
         NCMOI  = NORBI*NBASI
C
C        Symmetric Orthogonalization with N-R iterative algorithm:
C
C            C(p+1) = 1/2( 3C(p) - C(p) S Ct(p) C(p))
C
C                   = 3/2 C(p) + ( -1/2 C(p) Ct(p) ) C(p)
C
C            (here S = unit matrix)
C
         KCPNEW = KWRK0
         KCPOLD = KCPNEW + N2ORBI
         KWRK1  = KCPOLD + N2ORBI
         CALL MCOPY(NORBI,NORBI,WRK(KUXKAP+IORBI*NORBT+IORBI),NORBT,
     *              WRK(KCPNEW),NORBI)
C        CALL MCOPY(NROWA,NCOLA,A,NRDIMA,B,NRDIMB)
C
         DMX = D0
         ITS = 0
  100    CONTINUE
            ITS = ITS + 1
            CALL DCOPY(N2ORBI,WRK(KCPNEW),1,WRK(KCPOLD),1)
            CALL DSCAL(N2ORBI,D1P5,WRK(KCPNEW),1)
            CALL DGEMM('N','T',NORBI,NORBI,NORBI,1.D0,
     &                 WRK(KCPOLD),NORBI,
     &                 WRK(KCPOLD),NORBI,0.D0,
     &                 WRK(KWRK1),NORBI)
            CALL DSCAL(N2ORBI,-DP5,WRK(KWRK1),1)
            CALL DGEMM('N','N',NORBI,NORBI,NORBI,1.D0,
     &                 WRK(KWRK1),NORBI,
     &                 WRK(KCPOLD),NORBI,1.D0,
     &                 WRK(KCPNEW),NORBI)
C
            DMXM1 = DMX
            DMX   = D0
            DO 160 I = 0,N2ORBI-1
               DEL = ABS(WRK(KCPOLD+I)-WRK(KCPNEW+I))
               DMX = MAX(DMX,DEL)
  160       CONTINUE
            IF (DMX .GT. DIVERG*DMXM1 .AND. DMXM1 .GT. D0) THEN
               WRITE (LUPRI,2212) ISYM,ITS,DMX,DMXM1
               WRITE (LUPRI,3300)
               CALL OUTPUT(WRK(KCPNEW),1,NORBI,1,NORBI,
     *                     NORBI,NORBI,1,LUPRI)
               CALL QUIT(
     &         'MO_ROTATE, divergence in symmetric orthonormalization')
            END IF
            IF (DMX .GT. THREQL) GO TO 100
C        ^----------------------------------
C
C     Printflag here to check rotation matrix
C     (e.g. Norot specified elements should contain zeroes)
C
         IF (P6FLAG(12)) THEN
            WRITE (LUPRI,2211) ISYM,ITS,DMX,DMXM1
            WRITE (LUPRI,3300)
            CALL OUTPUT(WRK(KCPNEW),1,NORBI,1,NORBI,NORBI,NORBI,1,LUPRI)
         END IF
 2211 FORMAT(/' (MO_ROTATE) Symmetric orthogonalization of symmetry',I2,
     *        ' in',I3,' iterations.',
     *        /T11,'Max. error this and previous iteration',1P,2D15.2)
 2212 FORMAT(/' (MO_ROTATE) Symmetric orthogonalization of symmetry',I2,
     *        ' diverging, iteration no.',I3,
     *        /T11,'Max. error this and previous iteration',1P,2D15.2)
 3300 FORMAT(/' MO_ROTATE: UXKAP, the symmetrically orthogonalized',
     *        ' rotation matrix')
C
C        Rotate orbitals
C        CMOnew(q,r) = CMOold(q,s) Cpnew(r,s)
C
         KCMONW = KCPOLD
         CALL DGEMM('N','T',NBASI,NORBI,NORBI,1.D0,
     &              CMO(ICMOI+1),NBASI,
     &              WRK(KCPNEW),NORBI,0.D0,
     &              WRK(KCMONW),NBASI)
         CALL DCOPY(NCMOI,WRK(KCMONW),1,CMO(ICMOI+1),1)
C
  500 CONTINUE
C
C
C Section 2 **********************
C Read overlap matrix S from LUONEL
C
      KSOVLP = 1
      KWRK2  = KSOVLP + NNBAST
      KNEED  = KWRK2  + 2*NBASMA
      IF (KNEED .GT. LFRSAV) CALL ERRWRK('MO_ROTATE_2',KNEED,LFRSAV)
      CALL RDONEL('OVERLAP ',.TRUE.,WRK(KSOVLP),NNBAST)
C
C Section 3 **********************
C Gram-Schmidt orthogonalize orbitals symmetry by symmetry
C
      DO 400 ISYM=1,NSYM
         NORBI= NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 400
         NBASI= NBAS(ISYM)
         ISTS = KSOVLP + IIBAS(ISYM)
         ISTC = 1 + ICMO(ISYM)
         CALL NORM(WRK(ISTS),CMO(ISTC),NBASI,NORBI,WRK(KWRK2),
     *             THNORM,IRETUR)
         IF(IRETUR.GE.1)GO TO 4000
  400 CONTINUE
C
 9999 CALL QEXIT('MO_ROTATE')
      RETURN
C
 4000 WRITE(LUERR,4010)IRETUR
 4010 FORMAT(/' MO_ROTATE FATAL ERROR, linear dependency in norm',I5)
      CALL QTRACE(LUERR)
      CALL QUIT('ERROR, linear dependency (MO_ROTATE)')
C
C ** End of subroutine MO_ROTATE
C
      END
