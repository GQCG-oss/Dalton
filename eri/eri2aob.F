C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id: eri2aob.F,v 1.2 2001-02-12 18:17:53 vebjornb Exp $
C
#include <single.h>
C  /* Deck setaob */
      SUBROUTINE SETAOB(CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
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
      DO 100 K = 0, MAXREP
      DO 100 I = 1, KMAX
         IF (IBTAND(K,ISTBAO(I)).EQ.0) NBASE = NBASE + KHKT(I)
  100 CONTINUE
CNECgh980808
      NODD = 1 - MOD(NBASE,2) 
C
      NGAB = 0
      DO 110 I = 1, KMAX
         NGAB = NGAB + MULT(ISTBAO(I))*NUCO(I)
  110 CONTINUE
C
C     AO batches
C
      NB = 0
      NC = 0
      IORB = 1
      IORBRP(0) = 1
      DO 120 IREPA = 1, MAXREP
         IORBRP(IREPA) = IORBRP(IREPA-1) + NAOS(IREPA)
  120 CONTINUE
      DO 200 ISHELL = 1, KMAX
         IF (NB.NE.NBCH(ISHELL) .OR. NC.NE.NCENT(ISHELL)) THEN
            NB = NBCH(ISHELL)
            NC = NCENT(ISHELL)
            CALL AOBCH(NB,ISHELL,IORB,CCFBT,INDXBT,WORK,LWORK,IPRINT)
         END IF
C        IORB = IORB + MULT(ISTBAO(ISHELL))*KHKT(ISHELL)
         DO 210 IREPA = 0, MAXREP
            IORBRP(IREPA) = IORBRP(IREPA)
     &        + MLTCMP(NHKT(ISHELL),KHKT(ISHELL),ISTBAO(ISHELL),IREPA)
  210    CONTINUE
  200 CONTINUE
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Common block AOBTCH in SETAOB',-1)
         WRITE (LUPRI,'(1X,A,I3)') ' Number of batches NAOBCH:',NAOBCH
         CALL HEADER(
     &       '   #       NORBBT NHKT KCKT KHKT      NPRF NCTF      '//
     &       'NCNT      ISTB MULT',5)
         DO 900 I = 1, NAOBCH
            WRITE (LUPRI,'(5X,I5,I5,5X,3I5,5X,2I5,5X,I5,5X,2I5)')
     &          I,
     &          NORBBT(I),
     &          NHKTBT(I), KCKTBT(I), KHKTBT(I),
     &          NPRFBT(I), NCTFBT(I),
     &          NCNTBT(I),
     &          ISTBBT(I), MULTBT(I)
  900    CONTINUE
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Exponents and contraction coefficients',-1)
         DO 910 I = 1, NAOBCH
            WRITE (LUPRI,'(/,1X,A,I2,A,2I3,/)')
     &         ' Exp. and cont. coef. for AO batch',I,
     &         ' of dimensions',NPRFBT(I),NCTFBT(I)
            KEXP = KEXPBT(I) - 1
            KPRI = KCCFBT(I) - 1
            DO 920 J = 1, NPRFBT(I)
               WRITE (LUPRI,'(2X,E8.3,7(2X,E8.2))') EXPBT(KEXP + J),
     &             (CCFBT(KPRI + J + NPRFBT(I)*(K-1)),K=1,NCTFBT(I))
  920       CONTINUE
  910    CONTINUE
         CALL HEADER('Coordinates for AO batches',-1)
         DO 930 I = 1, NAOBCH
            WRITE (LUPRI,'(2X,I5,5X,3F12.6)')
     &          I, CORXBT(I), CORYBT(I), CORZBT(I)
  930    CONTINUE
         CALL HEADER('Orbital indices for AO batches',-1)
         DO 940 IREPA = 0, MAXREP
            WRITE (LUPRI,'(2X,A,I5)') 'Symmetry ',IREPA
            DO 950 I = 1, NAOBCH
               WRITE (LUPRI,'(2X,I5,3X,(22I3))')
     &             I, (INDXBT(KNDXBT(I) - 1 + J,IREPA),J=1,NCTFBT(I))
  950       CONTINUE
  940    CONTINUE
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
      DO 300 I = 1, NAOBCH
         NHKTA = NHKTBT(KAOSRT(I))
         IF (NHKTA.NE.NHKOLD) THEN
            KQNBT(NHKTA) = I
            NHKOLD = NHKTA
         END IF
         MAXQN = MAX(MAXQN,NHKTA)
  300 CONTINUE
      DO 400 I = 1, MAXQN - 1
         NQNBT(I) = KQNBT(I+1) - KQNBT(I)
  400 CONTINUE
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
      SUBROUTINE AOBCH(NB,ISHELL,IORB,CCFBT,INDXBT,WORK,LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
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
     &            ISHELL,IORB,CCFBT,INDXBT,IPRINT)
C
      RETURN
      END
C  /* Deck aobch1 */
      SUBROUTINE AOBCH1(NB,NMAT,MAT,NUCA,NRCA,JSTRA,IDONE,NDIMP,NDIMC,
     &                  IPNTC,IPNTP,ISHELL,IORB,CCFBT,INDXBT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxorb.h>
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
      DO 100 I = 1, NUCA
      DO 100 J = 1, NRCA
         IF (ABS(PRICCF(JSTRA+I,J)).GT.THR) MAT(I,J) = 1
  100 CONTINUE
C
      DO 200 I = 1, NRCA
      DO 200 J = 1, NRCA
       IF (INPRD(MAT(1,I),MAT(1,J),NUCA) .GT. 0) THEN
            DO 210 K = 1, NUCA
               MAT(K,I) = MAX(MAT(K,I),MAT(K,J))
               MAT(K,J) = MAX(MAT(K,I),MAT(K,J))
  210       CONTINUE
         END IF
  200 CONTINUE
C
C     NMAT, NDIMP(), NDIMC(), IPNTP(), IPNTC()
C
      NMAT = 0
      CALL IZERO(IDONE,NUCA)
      DO 300 I = 1, NRCA
         IF (IDONE(I).EQ.0) THEN
C
            NMAT = NMAT + 1
            DO 400 J = I, NRCA
               IDONE(J) =  INPRD(MAT(1,I),MAT(1,J),NUCA)
  400       CONTINUE
C
C           First non-zero primitive
C
            KCONT = I
            DO 410 J = 1, NUCA
            IF (MAT(J,KCONT) .EQ. 1) THEN
               KPRIM = J
               GO TO 410
            END IF
  410       CONTINUE
C
C           IPNTP and NDIMP
C
            IPRI = 0
            DO 420 J = 1, NUCA
            IF (MAT(J,KCONT) .EQ. 1) THEN
               IPRI = IPRI + 1
               IPNTP(NMAT,IPRI) = JSTRA + J
            END IF
  420       CONTINUE
            NDIMP(NMAT) = IPRI
C
C           IPNTC and NDIMC
C
            ICNT = 0
            DO 430 J = 1, NRCA
            IF (MAT(KPRIM,J) .EQ. 1) THEN
               ICNT = ICNT + 1
               IPNTC(NMAT,ICNT) = J
            END IF
  430       CONTINUE
            NDIMC(NMAT) = ICNT
C
         END IF
  300 CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from AOBCH1',-1)
         WRITE (LUPRI,'(1X,A,I3,A)')
     &      ' Input contraction matrix for block', NB,':'
         CALL OUTPUT(PRICCF(JSTRA+1,1),1,NUCA,1,NRCA,
     &               MXPRIM,NRCA,1,LUPRI)
         WRITE (LUPRI,'(/,1X,A,I3)') ' Number of AO batches:',NMAT
         DO 900 I = 1, NMAT
            WRITE (LUPRI,'(1X,A,I3,A,2I3)')
     &      ' Dimensions (prim. - cont.) for batch',I,':',
     &        NDIMP(I),NDIMC(I)
            WRITE (LUPRI,'(1X,A,17I3/,(29X,17I3))')
     &      ' Contributing primitives:  ',(IPNTP(I,J),J = 1,NDIMP(I))
            WRITE (LUPRI,'(1X,A,17I3/,(29X,17I3))')
     &      ' Contributing contracted:  ',(IPNTC(I,J),J = 1,NDIMC(I))
  900    CONTINUE
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
      DO 500 IMAT = 1, NMAT
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
C
         NORBBT(IBATCH)   = KHKT(ISHELL)*MULT(ISTBAO(ISHELL))
C
         KEXPBT(IBATCH) = IADRP
         DO 600 IPRI = 1, NDIMP(IMAT)
            EXPBT(IADRP) = PRIEXP(IPNTP(IMAT,IPRI))
            IADRP = IADRP + 1
  600    CONTINUE
C
         KCCFBT(IBATCH) = IADRC
         DO 610 ICNT = 1, NDIMC(IMAT)
         DO 610 IPRI = 1, NDIMP(IMAT)
            CCFBT(IADRC) = PRICCF(IPNTP(IMAT,IPRI),IPNTC(IMAT,ICNT))
            IADRC = IADRC + 1
  610    CONTINUE
C
         KNDXBT(IBATCH) = IADRO
C        DO 620 ICNT = 1, NDIMC(IMAT)
C           INDXBT(IADRO) = IORB
C    &          + (IPNTC(IMAT,ICNT)-1)*KHKT(ISHELL)*MULT(ISTBAO(ISHELL))
C           IADRO = IADRO + 1
C 620    CONTINUE
C
         DO 630 ICNT = 1, NDIMC(IMAT)
            DO 700 IREPA = 0, MAXREP
               INDXBT(IADRO,IREPA) = IORBRP(IREPA)+(IPNTC(IMAT,ICNT)-1)
     &           *MLTCMP(NHKT(ISHELL),KHKT(ISHELL),ISTBAO(ISHELL),IREPA)
  700       CONTINUE
            IADRO = IADRO + 1
  630    CONTINUE
C
  500 CONTINUE
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
#include <mxorb.h>
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
#include <mxorb.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aobtch.h>
C
C     Sorts by...
C
C     ... angular momentum
C
      AOBTGT = NHKTBT(IAO) .GT. NHKTBT(JAO)
      IF (NHKTBT(IAO) .EQ. NHKTBT(JAO)) THEN
C
C        ... number of primitives
C
         AOBTGT = NPRFBT(IAO) .LT. NPRFBT(JAO)
         IF (NPRFBT(IAO) .EQ. NPRFBT(JAO)) THEN
C
C           ... number of contracted
C
            AOBTGT = NCTFBT(IAO) .GT. NCTFBT(IAO)
            IF (NCTFBT(IAO) .EQ. NCTFBT(JAO)) THEN
C
C              ... multiplicity
C
               AOBTGT = ISTBBT(IAO) .GT. ISTBBT(JAO)
               IF (ISTBBT(IAO) .EQ. ISTBBT(JAO)) THEN
                   AOBTGT = NCNTBT(IAO) .GT. NCNTBT(JAO)
               END IF
            END IF
         END IF
      END IF
C
      RETURN
      END
C  /* Deck mltcmp */
      FUNCTION MLTCMP(NHKTA,KHKTA,MULA,IREPA)
#include <implicit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
#include <symmet.h>
#include <ibtfun.h>
C
      MLTCMP = 0
      DO 100 ICMPA = 1, KHKTA
         IPARA = IBTXOR(IREPA,ISYMAO(NHKTA,ICMPA))
         IF (IBTAND(MULA,IPARA) .EQ. 0) MLTCMP = MLTCMP + 1
  100 CONTINUE
C
      RETURN
      END
