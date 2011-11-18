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
C  /* Deck fltinp */
      SUBROUTINE FLTINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 4)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <cbiflt.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.RESPON', '.STOP  '/
C
      NEWDEF = (WORD .EQ. '*FLOAT ')
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
     *            '" not recognized in FLTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FLTINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  RESPRI = .TRUE.
               GO TO 100
    4          CONTINUE
                  CUT    = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in FLTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in FLTINP.')
            END IF
      END IF
  300 CONTINUE
      IF (IPRINT .GT. 0) RESPRI = .TRUE.
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for FLTORB:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' FLTORB skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in FLTORB:',IPRINT
            END IF
            IF (RESPRI) THEN
               WRITE (LUPRI,'(A)') ' Orbital responses will be printed.'
            END IF
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after FLTORB.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck fltini */
      SUBROUTINE FLTINI
C
C     Initialize /CBIFLT/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbiflt.h>
C
      IPRINT = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      RESPRI = .FALSE.
      RETURN
      END
C  /* Deck fltorb */
      SUBROUTINE FLTORB(WORK,LWORK)
C
C     Jan 87 tuh
C
C     The purpose of this section is to calculate the correction to
C     second-order properties (Hessian matrix, dipole derivatives, and
C     polarizabilities) which arises from the coordinate dependence
C     of floating orbitals. We assume the the position of all centers
C     with zero charge have been completely optimized.
C
C     The final second-order properties are stored in
C     HESMOL(MXCOOR,MXCOOR), DIPFLT(3,MXCOOR) and POLFLT(3,3).
C     If no floating orbitals are used this arrays are set
C     equal to HESMOL(MXCOOR,MXCOOR), DIP1(3,MXCOOR) and
C     POLFLT(3,3)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <iratdef.h>
      LOGICAL POLDIP
      DIMENSION WORK(LWORK)
#include <abainf.h>
#include <cbiflt.h>
#include <taymol.h>
#include <moldip.h>
#include <nuclei.h>
      IF (SKIP) RETURN
      IF (.NOT. (MOLHES .OR. DIPDER .OR. POLAR)) RETURN
      CALL QENTER('FLTORB')
      IF (IPRINT.GT.5 .OR. (IPRINT.GT.0 .AND. NFLOAT.GT.0)) THEN
         CALL TITLER('Output from FLTORB','*',103)
         WRITE (LUPRI,'(A,I5)') ' Number of floating orbitals : ',NFLOAT
      END IF
C     IF (NFLOAT .GT. 0) THEN
      IF (.FALSE.) THEN
C
C        Number of independent and dependent coordinates:
C
         NIND = 3*NATOMS
         POLDIP = POLAR .OR. DIPDER
         IF (POLDIP) NIND = NIND + 3
         NDEP = 3*NFLOAT
         IF (IPRINT. GT. 1) THEN
            WRITE (LUPRI,'(A,I5)')
     *         ' Number of independent coordinates: ',NIND
            WRITE (LUPRI,'(A,I5)')
     *         ' Number of dependent coordinates  : ',NDEP
         END IF
C
C        Work space allocation
C
         KRED = 1
         KGII = KRED + NIND*NIND
         KGID = KGII + NIND*NIND
         KGDI = KGID + NIND*NDEP
         KGDD = KGDI + NDEP*NIND
         KPCN = KGDD + NDEP*NDEP
         KIND = KPCN + 6*NDEP
         KDEP = KIND + (NIND + IRAT - 1)/IRAT
         KPVT = KDEP + (NDEP + IRAT - 1)/IRAT
         KWRK = KPVT + (NDEP + IRAT - 1)/IRAT
         IF (KWRK.GT.LWORK) CALL STOPIT('FLTORB',' ',KWRK,LWORK)
C
C        Set up matrices GII, GID, GDD and GDI
C
         CALL FLTGMT(WORK(KGII),WORK(KGID),WORK(KGDI),WORK(KGDD),
     *               WORK(KIND),WORK(KDEP),NIND,NDEP)
C
C        Fold in GDD: GRED = GII - GID*(GDD-1)*GDI
C
         CALL FOLDIN(WORK(KRED),WORK(KGII),WORK(KGID),WORK(KGDI),
     *               WORK(KGDD),WORK(KPVT),NIND,NDEP)
C
C        Print responses
C
         IF (RESPRI) THEN
            CALL FLTGDI(WORK(KGDI),NIND,NDEP,POLDIP)
         END IF
         IF (POLAR) THEN
            CALL FLTPAN(WORK(KGID),WORK(KGDI),WORK(KPCN),NIND,NDEP)
         END IF
C
C        Pick out matrices HESFLT, DIPFLT, POLFLT i from GRED
C
         CALL FLTHDP(WORK(KRED),NIND)
      ELSE
         IF (MOLHES) THEN
            DO 200 I = 1, 3*NUCDEP
               DO 300 J = 1, I
                  HESMOL(J,I) = HESMOL(I,J)
  300          CONTINUE
  200       CONTINUE
         END IF
         IF (DIPDER) THEN
            DO 400 I = 1, 3*NUCDEP
               DO 500 J = 1, 3
                  DIPFLT(J,I) = DIP1(J,I)
  500          CONTINUE
  400       CONTINUE
         END IF
         IF (POLAR) THEN
             DO 600 I = 1, 3
                DO 700 J = 1, 3
                   POLFLT(J,I) = POLARS(J,I)
  700           CONTINUE
  600        CONTINUE
         END IF
      END IF
      IF (CUT) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Program stopped after FLTORB as requested.'
         CALL QUIT(' ***** End of ABACUS (in FLTORB) *****')
      END IF
      CALL QEXIT('FLTORB')
      RETURN
      END
C  /* Deck fltgmt */
      SUBROUTINE FLTGMT(GII,GID,GDI,GDD,IVIND,IVDEP,NIND,NDEP)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <cbiflt.h>
#include <abainf.h>
#include <taymol.h>
#include <nuclei.h>
#include <moldip.h>
      DIMENSION GII(NIND,NIND), GID(NIND,NDEP), GDI(NDEP,NIND),
     *          GDD(NDEP,NDEP), IVIND(NIND), IVDEP(NDEP)
      LOGICAL DEPEND
      IF (IPRINT .GE. 5) CALL TITLER('Output from FLTGMT','*',103)
C
C     Set up address vectors
C
      IIND = 0
      IDEP = 0
      DO 100 ICOOR = 1, 3*NUCDEP
         DEPEND = NINT(CHARGE((ICOOR + 2)/3)) .EQ. 0
         IF (DEPEND) THEN
            IDEP = IDEP + 1
            IVDEP(IDEP) = ICOOR
         ELSE
            IIND = IIND + 1
            IVIND(IIND) = ICOOR
         END IF
  100 CONTINUE
      IF (DIPDER .OR. POLAR) THEN
         IVIND(IIND + 1) = - 1
         IVIND(IIND + 2) = - 2
         IVIND(IIND + 3) = - 3
      END IF
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('IVIND vector',-1)
         WRITE (LUPRI,'(2X,10I5)') (IVIND(I),I=1,NIND)
         CALL HEADER('IVDEP vector',-1)
         WRITE (LUPRI,'(2X,10I5)') (IVDEP(I),I=1,NDEP)
      END IF
C
C     Construct GII
C
      DO 200 I = 1, NIND
         DO 210 J = 1, NIND
            II = IVIND(I)
            JJ = IVIND(J)
            IF (II.GT.0 .AND. JJ.GT.0) GII(I,J) = HESMOL( II, JJ)
            IF (II.GT.0 .AND. JJ.LT.0) GII(I,J) = DIP1  (-JJ, II)
            IF (II.LT.0 .AND. JJ.GT.0) GII(I,J) = DIP1  (-II, JJ)
            IF (II.LT.0 .AND. JJ.LT.0) GII(I,J) = - POLARS(-II,-JJ)
  210    CONTINUE
  200 CONTINUE
C
C     Construct GID and GDI
C
      DO 300 I = 1, NIND
         DO 310 J = 1, NDEP
            II = IVIND(I)
            JJ = IVDEP(J)
            IF (II.GT.0) ELEMNT = HESMOL( II,JJ)
            IF (II.LT.0) ELEMNT = DIP1  (-II,JJ)
            GID(I,J) = ELEMNT
            GDI(J,I) = ELEMNT
  310    CONTINUE
  300 CONTINUE
C
C     Construct GDD
C
      DO 400 I = 1, NDEP
         DO 410 J = 1, NDEP
            GDD(I,J) = HESMOL(IVDEP(I),IVDEP(J))
  410    CONTINUE
  400 CONTINUE
C
C     Print section
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('GII matrix',-1)
         CALL OUTPUT(GII,1,NIND,1,NIND,NIND,NIND,1,LUPRI)
         CALL HEADER('GDI matrix',-1)
         CALL OUTPUT(GDI,1,NDEP,1,NIND,NDEP,NIND,1,LUPRI)
         CALL HEADER('GID matrix',-1)
         CALL OUTPUT(GID,1,NIND,1,NDEP,NIND,NDEP,1,LUPRI)
         CALL HEADER('GDD matrix',-1)
         CALL OUTPUT(GDD,1,NDEP,1,NDEP,NDEP,NDEP,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck foldin */
      SUBROUTINE FOLDIN(GRED,GII,GID,GDI,GDD,KPVT,NIND,NDEP)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <cbiflt.h>
#include <abainf.h>
      DIMENSION GII(NIND,NIND), GID(NIND,NDEP), GDI(NDEP,NIND),
     *          GDD(NDEP,NDEP), KPVT(NDEP),GRED(NIND,NIND)
C
      IF (IPRINT .GE. 5) CALL TITLER('Output from FOLDINT','*',103)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('GDD matrix',-1)
         CALL OUTPUT(GDD,1,NDEP,1,NDEP,NDEP,NDEP,1,LUPRI)
         CALL HEADER('GDI matrix',-1)
         CALL OUTPUT(GDI,1,NDEP,1,NIND,NDEP,NIND,1,LUPRI)
      END IF
      CALL DGESOL(NIND,NDEP,NDEP,NDEP,GDD,GDI,KPVT,INFO)
      IF (INFO .NE. 0) THEN
         WRITE (LUPRI,'(//,A,I5,A,/)')
     *      ' ERROR (FOLDIN) INFO = ',INFO, ' FROM DGESOL '
         CALL QUIT('(ABACUS.FLTORB.FOLDIN) ERROR in DGESOL ')
      END IF
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('GDD-1*GDI matrix',-1)
         CALL OUTPUT(GDI,1,NDEP,1,NIND,NDEP,NIND,1,LUPRI)
      END IF
      CALL DGEMM('N','N',NIND,NIND,NDEP,1.D0,
     &           GID,NIND,
     &           GDI,NDEP,0.D0,
     &           GRED,NIND)
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('GID*GDD-1*GDI matrix',-1)
         CALL OUTPUT(GRED,1,NIND,1,NIND,NIND,NIND,1,LUPRI)
      END IF
      DO 100 I = 1, NIND
         DO 200 J = 1, NIND
            GRED(I,J) = GII(I,J) - GRED(I,J)
  200    CONTINUE
  100 CONTINUE
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('GRED matrix',-1)
         CALL OUTPUT(GRED,1,NIND,1,NIND,NIND,NIND,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck flthdp */
      SUBROUTINE FLTHDP(GRED,NIND)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <cbiflt.h>
#include <abainf.h>
#include <taymol.h>
#include <nuclei.h>
#include <moldip.h>
      DIMENSION GRED(NIND,NIND)
C
      IF (IPRINT. GE .5) CALL TITLER('Output from FLTHDP','*',103)
      NNUC = NIND
      IF (DIPDER .OR. POLAR) NNUC = NIND - 3
C
C     Hessian matrix
C
C      DO 100 I = 1, NNUC
C         DO 200 J = 1, NNUC
C            HESFLT(J,I) = GRED(J,I)
C  200    CONTINUE
C  100 CONTINUE
C
C     Dipole gradient
C
      IF (DIPDER) THEN
         DO 300 I = 1, NNUC
            DO 400 J = 1, 3
               DIPFLT(J,I) = GRED(NNUC + J,I)
  400       CONTINUE
  300    CONTINUE
      END IF
C
C     Polarizabilities
C
      IF (POLAR) THEN
         DO 500 I = 1, 3
            DO 600 J = 1, 3
               POLFLT(J,I) = - GRED(NNUC + J,NNUC + I)
  600       CONTINUE
  500    CONTINUE
      END IF
      RETURN
      END
C  /* Deck fltgdi */
      SUBROUTINE FLTGDI(RESPON,NIND,NDEP,POLDIP)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      LOGICAL POLDIP
      DIMENSION RESPON(NDEP,NIND)
      CHARACTER*6 LABEL(5), NAME
#include <nuclei.h>
C
      WRITE (LUPRI,'(//)')
      CALL HEADER('Orbital Center Responses',-1)
      IOFFNM = NIND
      IF (POLDIP) IOFFNM = NIND - 3
      ISTR = 1
      NBATCH = (NIND + 4)/5
      DO 100 IBATCH = 1, NBATCH
         IEND = MIN(ISTR + 4,NIND)
         NUMB = IEND - ISTR + 1
         DO 110 I = 1, 5
            IADR = ISTR - 1 + I
            NAME = NAMEX(IADR)
            IF (IADR .EQ. IOFFNM + 1) NAME = '  EX  '
            IF (IADR .EQ. IOFFNM + 2) NAME = '  EY  '
            IF (IADR .EQ. IOFFNM + 3) NAME = '  EZ  '
            LABEL(I) = NAME
  110    CONTINUE
         WRITE (LUPRI,'(/16X,5(A6,6X))')
     *      (LABEL(I), I = 1, NUMB)
         LENH = 10 + NUMB*12
         WRITE (LUPRI,'(2X,70A1)') ('-', II = 1,LENH)
         WRITE (LUPRI,'()')
         DO 200 ICOOR = 1, NDEP
            WRITE (LUPRI,1000) NAMEX(IOFFNM + ICOOR),
     *         (RESPON(ICOOR,II),II=ISTR,IEND)
            IF (MOD(ICOOR,3) .EQ. 0) WRITE (LUPRI,'()')
 200     CONTINUE
         ISTR = ISTR + 5
 100  CONTINUE
      RETURN
 1000 FORMAT (5X,A6,5F12.6)
      END
      SUBROUTINE FLTPAN(GID,GDI,PCNT,NIND,NDEP)
C
C     The purpose of this subroutine is to calculate the contributions
C     from each floating orbital to the polarizability corrections.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D00, D100 = 100.0D00, THRESH = 1.0D-04)
#include <nuclei.h>
      DIMENSION GID(NIND,NDEP), GDI(NDEP,NIND), PCNT(6,NDEP)
C
      WRITE (LUPRI,'()')
      CALL DZERO(PCNT,6*NDEP)
      NAMOFF = 3*NATOMS
      IJ = 0
      DO 100 I = 1, 3
         IADR = NIND - 3 + I
         DO 110 J = I, 3
            JADR = NIND - 3 + J
            IJ = IJ + 1
C
C           Calculate full polarizability correction
C
            POLCOR = D0
            DO 200 K = 1, NDEP
               POLCOR = POLCOR + GID(IADR,K)*GDI(K,JADR)
  200       CONTINUE
            POLCOR = POLCOR/D100
C
C           Calculate fraction of correction from each orbital center
C
            IF (ABS(POLCOR) .GT. THRESH) THEN
               DO 300 K = 1, NDEP
                  PCNT(IJ,K) = GID(IADR,K)*GDI(K,JADR)/POLCOR
  300          CONTINUE
            END IF
  110    CONTINUE
  100 CONTINUE
      CALL HEADER('Contributions To Polarizabilities '//
     *            'From Each Orbital Center',-1)
      CALL HEADER(
     *   '          XX      XY      XZ      YY      YZ      ZZ  ',-6)
      DO 400 I = 1, NDEP
         WRITE (LUPRI,'(10X,A,2X,6(F5.1,3X))') NAMEX(NAMOFF + I),
     *                                        (PCNT(K,I),K=1,6)
         IF (MOD(I,3) .EQ. 0) WRITE (LUPRI,'()')
  400 CONTINUE
      RETURN
      END
