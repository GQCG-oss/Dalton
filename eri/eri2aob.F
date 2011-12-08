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
C  /* Deck setaob */
      SUBROUTINE SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
      DIMENSION CCFBT(MXPRIM*MXCONT), INDXBT(MXSHEL*MXCONT,0:7),
     &          WORK(LWORK)
#include <shells.h>
#include <symmet.h>
#include <aobtch.h>
#include <veclen.h> 
#include <ibtfun.h>
C
C     Number of basis functions
C
      NBASE  = 0
      DO K = 0, MAXREP
      DO I = 1, KMAXTOT
         IF (IBTAND(K,ISTBAO(I)).EQ.0) NBASE = NBASE + KHKT(I)
      END DO
      END DO
C
      NODD = 1 - MOD(NBASE,2)
C
      NGAB = 0
      DO I = 1, KMAXTOT
         NGAB = NGAB + MULT(ISTBAO(I))*NUCO(I)
      END DO
C
      IORBRP(0) = 1
      DO IREPA = 1, MAXREP
         IORBRP(IREPA) = IORBRP(IREPA-1) + NAOS(IREPA)
      END DO
C
C     AO batches
C
      NB = -1
      NC = -1
      ND = -1
      DO ISHELL = 1, KMAXTOT
         IF (NB.NE.NBCH(ISHELL) .OR. NC.NE.NCENT(ISHELL) .OR.
     &       ND.NE.LCLASS(ISHELL) ) THEN
            NB = NBCH(ISHELL)
            NC = NCENT(ISHELL)
            ND = LCLASS(ISHELL)
            CALL AOBCH(NB,ISHELL,CCFBT,INDXBT,WORK,LWORK,IPRINT)
         END IF
         DO IREPA = 0, MAXREP
            IORBRP(IREPA) = IORBRP(IREPA)
     &        + MLTCMP(NHKT(ISHELL),KHKT(ISHELL),ISTBAO(ISHELL),IREPA)
         END DO
      END DO
C
C     Identity contraction matrices
C
      CALL CNTTYP(CCFBT,.TRUE.,IPRINT)
C
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Common block AOBTCH in SETAOB',-1)
         WRITE (LUPRI,'(1X,A,I3)') ' Number of batches NAOBCH:',NAOBCH
         CALL HEADER(
     &       '   #  NORB NHKT KCKT KHKT NPRF NCTF KCMT'//
     &       ' NCNT KCLS ISTB MULT',1)
         DO I = 1, NAOBCH
            WRITE (LUPRI,'(1X,12I5)')
     &          I,
     &          NORBBT(I),
     &          NHKTBT(I), KCKTBT(I), KHKTBT(I),
     &          NPRFBT(I), NCTFBT(I), KCMTBT(I),
     &          NCNTBT(I), KCLSBT(I),
     &          ISTBBT(I), MULTBT(I)
         END DO
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Exponents and contraction coefficients',-1)
         DO I = 1, NAOBCH
            WRITE (LUPRI,'(/,1X,A,I2,A,2I3,/)')
     &         ' Exp. and cont. coef. for AO batch',I,
     &         ' of dimensions',NPRFBT(I),NCTFBT(I)
            KEXP = KEXPBT(I) - 1
            KPRI = KCCFBT(I) - 1
            DO J = 1, NPRFBT(I)
               WRITE (LUPRI,'(2X,E8.3,7(2X,E8.2))') EXPBT(KEXP + J),
     &             (CCFBT(KPRI + J + NPRFBT(I)*(K-1)),K=1,NCTFBT(I))
            END DO
         END DO
         CALL HEADER('Coordinates for AO batches',-1)
         DO I = 1, NAOBCH
            WRITE (LUPRI,'(2X,I5,5X,3F12.6)')
     &          I, CORXBT(I), CORYBT(I), CORZBT(I)
         END DO
         CALL HEADER('Orbital indices for AO batches',-1)
         DO IREPA = 0, MAXREP
            WRITE (LUPRI,'(2X,A,I5)') 'Symmetry ',IREPA
            DO I = 1, NAOBCH
               WRITE (LUPRI,'(2X,I5,3X,(22I3))')
     &             I, (INDXBT(KNDXBT(I) - 1 + J,IREPA),J=1,NCTFBT(I))
            END DO
         END DO
      END IF
C
C     Sort AO batches
C
      CALL AOBSRT(IPRINT)
C
C     MAXQN, KQNBT, NQNBT
C
      MAXQN = 0
      NHKOLD = 0
      DO I = 1, NAOBCH
         NHKTA = NHKTBT(KAOSRT(I))
         IF (NHKTA.NE.NHKOLD) THEN
            KQNBT(NHKTA) = I
            NHKOLD = NHKTA
         END IF
         MAXQN = MAX(MAXQN,NHKTA)
      END DO
      DO I = 1, MAXQN - 1
         NQNBT(I) = KQNBT(I+1) - KQNBT(I)
      END DO
      NQNBT(MAXQN) = NAOBCH - KQNBT(MAXQN) + 1
C
      IF (IPRINT .GT. 3) THEN
         WRITE (LUPRI,'(/,1X,A,I3)')
     &      ' Highest angular momentum:',MAXQN - 1
         WRITE (LUPRI,'(1X,A,10I3)')
     &      ' Start addresses for ang. mom.:     ',(KQNBT(I),I=1,MAXQN)
         WRITE (LUPRI,'(1X,A,10I3)')
     &      ' Number of AO batches for ang. mom.:',(NQNBT(I),I=1,MAXQN)
      END IF
C
      RETURN
      END
C  /* Deck aobch */
      SUBROUTINE AOBCH(NB,ISHELL,CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      DIMENSION CCFBT(*), INDXBT(*), WORK(LWORK)
#include <shells.h>
C
      NUCA  = NUCO(ISHELL)
      NRCA  = NRCO(ISHELL)
      JSTRA = JSTRT(ISHELL)
C
      KMAT   = 1
      KDONE  = KMAT   + NUCA*NRCA
      KDIMP  = KDONE  + NUCA
      KDIMC  = KDIMP  + NUCA
      KPNTC  = KDIMC  + NUCA
      KPNTP  = KPNTC  + NUCA*NRCA
      KLAST  = KPNTP  + NUCA*NUCA
      IF (KLAST .GT. LWORK) CALL STOPIT('AOBCH',' ',KLAST,LWORK)
C
      CALL AOBCH1(NB,NMAT,WORK(KMAT),NUCA,NRCA,JSTRA,WORK(KDONE),
     &            WORK(KDIMP),WORK(KDIMC),WORK(KPNTC),WORK(KPNTP),
     &            ISHELL,CCFBT,INDXBT,IPRINT)
C
      RETURN
      END
C  /* Deck aobch1 */
      SUBROUTINE AOBCH1(NB,NMAT,MAT,NUCA,NRCA,JSTRA,IDONE,NDIMP,NDIMC,
     &                  IPNTC,IPNTP,ISHELL,CCFBT,INDXBT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
      PARAMETER (THR = 1.0D-17)
C
      DIMENSION MAT(NUCA,NRCA), IDONE(NUCA), NDIMP(NUCA), NDIMC(NUCA),
     &          IPNTC(NUCA,NRCA),IPNTP(NUCA,NUCA)
      DIMENSION CCFBT(MXPRIM*MXCONT), INDXBT(MXSHEL*MXCONT,0:7)
C
#include <symmet.h>
#include <aobtch.h>
#include <shells.h>
#include <primit.h>
#include <erisel.h>
C
      SAVE IADRP, IADRC, IADRO, NBTCH
C
C     Construct matrix MAT containing 1 and 0 only, then fill
C     out all elements within each contraction matrix
C
      CALL IZERO(MAT,NUCA*NRCA)
C
      DO I = 1, NUCA
      DO J = 1, NRCA
         IF (ABS(PRICCF(JSTRA+I,J)).GT.THR) MAT(I,J) = 1
      END DO
      END DO
C
      DO I = 1, NRCA
      DO J = 1, NRCA
         IF (INPRD(MAT(1,I),MAT(1,J),NUCA) .GT. 0) THEN
            DO K = 1, NUCA
               MAT(K,I) = MAX(MAT(K,I),MAT(K,J))
               MAT(K,J) = MAX(MAT(K,I),MAT(K,J))
            END DO
         END IF
      END DO
      END DO
C
C     NMAT, NDIMP(), NDIMC(), IPNTP(), IPNTC()
C
      NMAT = 0
      CALL IZERO(IDONE,NUCA)
      DO I = 1, NRCA
         IF (IDONE(I).EQ.0) THEN
C
            NMAT = NMAT + 1
            DO J = I, NRCA
               IDONE(J) =  INPRD(MAT(1,I),MAT(1,J),NUCA)
            END DO
C
C           First non-zero primitive
C
            KCONT = I
            DO J = 1, NUCA
               IF (MAT(J,KCONT) .EQ. 1) THEN
                  KPRIM = J
                  GO TO 100
               END IF
            END DO
  100       CONTINUE
C
C           IPNTP and NDIMP
C
            IPRI = 0
            DO J = 1, NUCA
            IF (MAT(J,KCONT) .EQ. 1) THEN
               IPRI = IPRI + 1
               IPNTP(NMAT,IPRI) = JSTRA + J
            END IF
            END DO
            NDIMP(NMAT) = IPRI
C
C           IPNTC and NDIMC
C
            ICNT = 0
            DO J = 1, NRCA
            IF (MAT(KPRIM,J) .EQ. 1) THEN
               ICNT = ICNT + 1
               IPNTC(NMAT,ICNT) = J
            END IF
            END DO
            NDIMC(NMAT) = ICNT
C
         END IF
      END DO
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from AOBCH1',-1)
         WRITE (LUPRI,'(1X,A,I3,A)')
     &      ' Input contraction matrix for block', NB,':'
         CALL OUTPUT(PRICCF(JSTRA+1,1),1,NUCA,1,NRCA,
     &               MXPRIM,NRCA,1,LUPRI)
         WRITE (LUPRI,'(/,1X,A,I3)') ' Number of AO batches:',NMAT
         DO I = 1, NMAT
            WRITE (LUPRI,'(1X,A,I3,A,2I3)')
     &      ' Dimensions (prim. - cont.) for batch',I,':',
     &        NDIMP(I),NDIMC(I)
            WRITE (LUPRI,'(1X,A,17I3/,(29X,17I3))')
     &      ' Contributing primitives:  ',(IPNTP(I,J),J = 1,NDIMP(I))
            WRITE (LUPRI,'(1X,A,17I3/,(29X,17I3))')
     &      ' Contributing contracted:  ',(IPNTC(I,J),J = 1,NDIMC(I))
         END DO
      END IF
C
C     Assign to COMMON /AOBTCH/
C
      IF (ISHELL .EQ. 1) THEN
         IADRP = 1
         IADRC = 1
         IADRO = 1
         NBTCH = 0
      END IF
C
      DO IMAT = 1, NMAT
         IBATCH = NBTCH + IMAT
         CORXBT(IBATCH) = CENT(ISHELL,1,1)
         CORYBT(IBATCH) = CENT(ISHELL,2,1)
         CORZBT(IBATCH) = CENT(ISHELL,3,1)
         NHKTBT(IBATCH) = NHKT(ISHELL)
         KCKTBT(IBATCH) = KCKT(ISHELL)
         KHKTBT(IBATCH) = KHKT(ISHELL)
         NPRFBT(IBATCH) = NDIMP(IMAT)
         NCTFBT(IBATCH) = NDIMC(IMAT)
         ISTBBT(IBATCH) = ISTBAO(ISHELL)
         MULTBT(IBATCH) = MULT(ISTBAO(ISHELL))
         NCNTBT(IBATCH) = NCENT(ISHELL)
         KCLSBT(IBATCH) = LCLASS(ISHELL)
         NORBBT(IBATCH) = NDIMC(IMAT)*KHKT(ISHELL)*MULT(ISTBAO(ISHELL))
C        Basis-set identifier (WK/UniKA/04-11-2002).
         MBIDBT(IBATCH) = MBSID(ISHELL)
C
C        exponents
C
         KEXPBT(IBATCH) = IADRP
         DO IPRI = 1, NDIMP(IMAT)
            EXPBT(IADRP) = PRIEXP(IPNTP(IMAT,IPRI))
            IADRP = IADRP + 1
         END DO
C
C        contraction coefficients
C
         KCCFBT(IBATCH) = IADRC
         DO ICNT = 1, NDIMC(IMAT)
         DO IPRI = 1, NDIMP(IMAT)
            CCFBT(IADRC) = PRICCF(IPNTP(IMAT,IPRI),IPNTC(IMAT,ICNT))
            IADRC = IADRC + 1
         END DO
         END DO
C
C        orbital indices
C
         KNDXBT(IBATCH) = IADRO
         DO ICNT = 1, NDIMC(IMAT)
            DO IREPA = 0, MAXREP
               INDXBT(IADRO,IREPA) = IORBRP(IREPA)+(IPNTC(IMAT,ICNT)-1)
     &           *MLTCMP(NHKT(ISHELL),KHKT(ISHELL),ISTBAO(ISHELL),IREPA)
            END DO
            IADRO = IADRO + 1
         END DO
C
      END DO
C
      NBTCH = NBTCH + NMAT
      NAOBCH = NBTCH
C
      RETURN
      END
C  /* Deck inprd */
      FUNCTION INPRD(IVEC1,IVEC2,NDIM)
#include <implicit.h>
      DIMENSION IVEC1(NDIM), IVEC2(NDIM)
      INNER = 0
      DO 100 I = 1, NDIM
         INNER = INNER + IVEC1(I)*IVEC2(I)
  100 CONTINUE
      INPRD = INNER
      RETURN
      END
C  /* Deck aobsrt */
      SUBROUTINE AOBSRT(IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxcent.h>
      LOGICAL AOBTGT
#include <aobtch.h>
C
      DO 100 I = 1, NAOBCH
         KAOSRT(I) = I
  100 CONTINUE
C
      DO 200 I = 2, NAOBCH
         IAO = KAOSRT(I)
         DO 300 J = I - 1, 1, -1
            IF (AOBTGT(IAO,KAOSRT(J))) GO TO 400
            KAOSRT(J+1) = KAOSRT(J)
  300    CONTINUE
         J = 0
  400    KAOSRT(J+1) = IAO
  200 CONTINUE
C
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Sorted AO batch list in AOBSRT',-1)
         WRITE(LUPRI,'(1X,A,/)') '   #       old      NHKT      NPRF'//
     &                         ' NCTF      ISTB MULT      NCNT'
         DO 500 I = 1, NAOBCH
            K = KAOSRT(I)
            WRITE (LUPRI,'(1X,I4,5X,I5,5X,I5,5X,2I5,5X,2I5,5X,I5)')
     &         I, K, NHKTBT(K), NPRFBT(K), NCTFBT(K),
     &         ISTBBT(K), MULTBT(K), NCNTBT(K)
  500    CONTINUE
      END IF
      RETURN
      END
C  /* Deck aobtgt */
      LOGICAL FUNCTION AOBTGT(IAO,JAO)
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aobtch.h>
C
C     Sorts by...
C
C     ... class
C
C     AOBTGT = KCLSBT(IAO) .GT. KCLSBT(JAO)
C     IF (KCLSBT(IAO) .EQ. KCLSBT(JAO)) THEN
C
C        ... angular momentum
C
         AOBTGT = NHKTBT(IAO) .GT. NHKTBT(JAO)
         IF (NHKTBT(IAO) .EQ. NHKTBT(JAO)) THEN
C
C           ... number of primitives
C
            AOBTGT = NPRFBT(IAO) .LT. NPRFBT(JAO)
            IF (NPRFBT(IAO) .EQ. NPRFBT(JAO)) THEN
C
C              ... number of contracted
C
               AOBTGT = NCTFBT(IAO) .GT. NCTFBT(IAO)
               IF (NCTFBT(IAO) .EQ. NCTFBT(JAO)) THEN
C
C                 ... multiplicity
C
                  AOBTGT = ISTBBT(IAO) .GT. ISTBBT(JAO)
                  IF (ISTBBT(IAO) .EQ. ISTBBT(JAO)) THEN
                      AOBTGT = NCNTBT(IAO) .GT. NCNTBT(JAO)
                  END IF
               END IF
            END IF
         END IF
C     END IF
C
      RETURN
      END
C  /* Deck mltcmp */
      FUNCTION MLTCMP(NHKTA,KHKTA,MULA,IREPA)
#include <implicit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <ibtfun.h>
C
      MLTCMP = 0
      DO ICMPA = 1, KHKTA
         IPARA = IBTXOR(IREPA,ISYMAO(NHKTA,ICMPA))
         IF (IBTAND(MULA,IPARA) .EQ. 0) MLTCMP = MLTCMP + 1
      END DO
C
      RETURN
      END
C  /* Deck cnttyp */
      SUBROUTINE CNTTYP(CCFBT,SAMCOF,IPRINT)
C
C     T. Helgaker Jan 2001
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0)
      PARAMETER (THR = 1.0D-17)
      LOGICAL SAMCOF
      DIMENSION CCFBT(*)
#include <aobtch.h>
C
      ITYPE = 0
      CALL IZERO(KCMTBT,NAOBCH)
      DO 100 I = 1, NAOBCH
      IF (KCMTBT(I).EQ.0) THEN
         ITYPE = ITYPE + 1
         KCMTBT(I) = ITYPE
         NPRIMI = NPRFBT(I)
         NCONTI = NCTFBT(I)
         KPRI = KCCFBT(I) - 1
         DO 200 J = I + 1, NAOBCH
         IF (KCMTBT(J).EQ.0) THEN
            IF (NPRIMI.EQ.NPRFBT(J) .AND. NCONTI.EQ.NCTFBT(J)) THEN
               IF (SAMCOF) THEN
                  DIFMAX = D0
                  KPRJ = KCCFBT(J) - 1
                  DO K = 1, NPRIMI*NCONTI
                     DIFMAX=MAX(DIFMAX,ABS(CCFBT(KPRI+K)-CCFBT(KPRJ+K)))
                  END DO
                  IF (DIFMAX .LT. THR) KCMTBT(J) = ITYPE
               ELSE
                  KCMTBT(J) = ITYPE
               END IF
            END IF
         END IF
  200    CONTINUE
      END IF
  100 CONTINUE
C
      RETURN
      END
