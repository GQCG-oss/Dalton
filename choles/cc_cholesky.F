C  /* Deck cc_cholesky */
      SUBROUTINE CC_CHOLESKY(DIAG,DIAOLD,INDRED,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez
C
C     Calculation of Cholesky decomposition based on an
C     integral-direct two-index approach.
C
C     NB! a non-direct approach has not been implemented.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (ONTM10 = 1.0D-10, ONTM18 = 1.0D-18)
      PARAMETER (MXCHRD = 200)
      PARAMETER (MAXQUA = 75)
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
      LOGICAL SETUP, LAST, DECO
      LOGICAL NEWSYM
      DIMENSION DIAG(*),DIAOLD(*) 
      DIMENSION INDRED(*)
      DIMENSION DIAMIN(8)
      DIMENSION WORK(LWORK)
      DIMENSION INDEXA(MXCORB), INDEXB(MXCORB)
      DIMENSION IDNTC1(MAXQUA)
      DIMENSION IQUVAB(MAXQUA,8), IQUVEA(MAXQUA,8), IQUVEB(MAXQUA,8)
      DIMENSION NQUAL(8), ITMP(8)
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <distcl.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccfro.h>
#include <ccfop.h>
#include <ccsections.h>
#include <ccdeco.h>
#include <choles.h>
#include <choskp.h>
#include <chodbg2.h>
#include <symsq.h>
#include <chotim.h>
#include <ccdeco2.h>
#include <choio.h>
C
      CHARACTER*10 NAME
      CHARACTER*7  ROOT
      DATA ROOT    /'CHOLES_'/
C
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C
      DIMENSION XSEL0(8),YSEL0(8)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C     Initializations.
C     ----------------
C
      TIMTOT = SECOND()
      TIMINI = SECOND()
      TIMINT = ZERO
C
      IPRINT = IPRCHO
C
      ICHOV  = 0
      NCHOV2 = 0
      NINT2  = 0
C
      XDISC   = 0.0D0
      XNDISC  = 0.0D0
C
C
C
C     Open files for storing the Cholesky vectors.
C     --------------------------------------------
C
      IF (.NOT. RSTCHO)
     &    CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'INI')
C
      LUDIAG = -1
      CALL GPOPEN(LUDIAG,'CHODIAG','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUM,.FALSE.)
C
      NAME  = ROOT(1:7)//'SEL'
C
      LUSEL = 98
      CALL WOPEN2(LUSEL,NAME,64,0)
C
C     Initialize integral calculation.
C     --------------------------------
C
      NTOINT = 0
C
      T1 = SECOND()
      JPRINT = 0
      CALL HERDI1(WORK,LWORK,JPRINT)
      T1 = SECOND() - T1
      TIMINT = TIMINT + T1
      IF (IPRINT .GT. 0) 
     &    WRITE(LUPRI,'(A,I5,/)')  'Number of shells:',MAXSHL
C
C     Calculate diagonal elements.
C     ----------------------------
C
      TIMDIA = SECOND()
      REWIND(LUDIAG)
      IF (RSTDIA) THEN
         READ(LUDIAG) (DIAG(I),I = 1,NDIAG)
      ELSE
c        CALL CC_GABOLD(DIAG,WORK,LWORK)
         CALL CC_GAB(DIAG,WORK,LWORK)
         WRITE(LUDIAG) (DIAG(I),I = 1,NDIAG)
         RSTDIA = .TRUE.
         CALL FLSHFO(LUDIAG)
      END IF
      CALL DCOPY(NDIAG,DIAG,1,DIAOLD,1)
C
C     Analyse diagonal by Cauchy-Schwarz.
C     ------------------------------------
C
      IF (IPRINT .GE. 10) THEN
      CALL AROUND('Cauchy-Schwarz screening based on initial diagonal')
      XTHR = THRDEF**2
C
      DO ISYM = 1,NSYM
C
         XC = ZERO
         DO I = 1,NNBST(ISYM)
            KOFF1 = IDIAG(ISYM) + I
            DO J = 1,I
               KOFF2 = IDIAG(ISYM) + J
               IF (DIAG(KOFF1)*DIAG(KOFF2) .LT. XTHR) THEN
                  XC = XC + 1.0D0 
               ENDIF
            ENDDO
         ENDDO
C
         XDIM = DFLOAT(NNBST(ISYM))*DFLOAT(NNBST(ISYM)+1)/2.0D0
         WRITE(LUPRI,*) 'DIAGONAL SCREENING: ',ISYM,XC,XDIM,
     *                  (XC/XDIM)*100.0D0
         CALL FLSHFO(LUPRI)
C
      ENDDO
      END IF
C
C     Zero out converged elements.
C     ----------------------------
C
      IF (SCDIAG) THEN
         IZER = 0
         DO ISYM = 1,NSYM
C
            XM = ZERO
            DO I = 1,NNBST(ISYM)
               KOFF1 = IDIAG(ISYM) + I
               IF (XM .LT. ABS(DIAG(KOFF1))) XM = ABS(DIAG(KOFF1))
            END DO
C
            XSCREEN = THRDEF*THRDEF/XM/THINDI
C
            DO I = 1,NNBST(ISYM)
               KOFF1 = IDIAG(ISYM) + I
               IF (ABS(DIAG(KOFF1)) .LT. XSCREEN) THEN
                     DIAG(KOFF1) = ZERO
                     IZER = IZER + 1
               END IF
            END DO
C
         END DO
      END IF
C
C     Get reduced set.
C     ----------------
C
      IF (REDUCE) THEN
C
         LRDTOT = 0
         DO ISYM = 1,NSYM
            ICOUNT = 0
            DO I = 1,NNBST(ISYM)
               KOFF1 = IDIAG(ISYM) + I
               IF (DIAG(KOFF1) .NE. 0) THEN
                  ICOUNT = ICOUNT + 1
                  KOFF2 = LREDU*(ISYM-1) + ICOUNT
                  INDRED(KOFF2) = I
               END IF
            END DO
            NREDUC(ISYM) = ICOUNT
            LRDTOT = LRDTOT + NREDUC(ISYM)
         END DO
C
         IF (IPRINT .GT. 0) THEN
             WRITE(LUPRI,'(//,2(9X,A,/))') 'Diagonal information',
     &                                     '--------------------'
             WRITE(LUPRI,'(6X,A,8X,A,7X,A)') 'ISAB',' NNBST ','NREDUC'
             DO ISYM = 1,NSYM
                WRITE(LUPRI,'(7X,I2,2I16)') ISYM,NNBST(ISYM),
     &                                      NREDUC(ISYM)
             END DO
             WRITE(LUPRI,*)
             WRITE(LUPRI,'(4X,A5,2I16)') 'Total', NDIAG, LRDTOT
             CALL FLSHFO(LUPRI)
         END IF
         LRDTOT = LREDU * NSYM
C
      ELSE
         IF (IPRINT .GT. 0) THEN
             WRITE(LUPRI,'(//,2(9X,A,/))') 'Diagonal information',
     &                                     '--------------------'
             WRITE(LUPRI,'(6X,A,8X,A,7X,A)') 'ISAB',' NNBST '
             DO ISYM = 1,NSYM
                WRITE(LUPRI,'(7X,I2,2I16)') ISYM,NNBST(ISYM)
             END DO
             WRITE(LUPRI,*)
             WRITE(LUPRI,'(4X,A5,2I16)') 'Total', NDIAG
             CALL FLSHFO(LUPRI)
         END IF
         LRDTOT = 1
      END IF
C
C     Analyse the diagonal.
C     ---------------------
C
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,'(//,A,/)') '>>> Checking initial diagonal'
         CALL CC_ANADI(DIAG,.TRUE.)
      END IF
      TIMDIA = SECOND() - TIMDIA
C
C     Setup index information.
C     ------------------------
C
      ISHLB  = 1
      JPRINT = 0
      SETUP  = .TRUE.
      DO ISHLA = 1,MAXSHL
C
         T1 = SECOND()
         JPRINT = 0
         CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     *               NUMDIA,NUMDIB,JPRINT,SETUP,1)
         T1 = SECOND() - T1
         TIMINT = TIMINT + T1
C
         DO I = 1,NUMDIA
            INAOSH(INDEXA(I)) = ISHLA
         ENDDO
      ENDDO
C
C     Restart.
C     --------
C
      IF (RSTCHO) THEN
C
         CALL CHO_NAME(ISYM,IFILE,LUCHO,NAME,'RST')
         CALL CHO_RESTART(DIAG,INDRED,WORK,LWORK)
C
      END IF
C
C     Allocation stuff.
C     -----------------
C
      MXNNBS = 0
      MXREDU = 0
      DO I = 1,NSYM
         IF (NNBST(I)  .GT. MXNNBS) MXNNBS = NNBST(I)
         IF (NREDUC(I) .GT. MXREDU) MXREDU = NREDUC(I)
      END DO
C
      LWRKX = LWORK / 3
C
      LENAB = MXNNBS + (MXNNBS-1)/IRAT + 2
C
      IF (REDUCE) THEN
         LENAB = MXREDU
      ELSE IF (.NOT. COMP) THEN
         LENAB = MXNNBS
      ELSE
         LENAB = MXNNBS + (MXNNBS-1)/IRAT + 2
      END IF
C
      LREAD = (NBAST*NBAST + 1)/IRAT + 1 
     &      + (NIBUF*LBUF-1)/IRAT + 1 + LBUF
     &      + (2*MAXQUA-1)/IRAT + 1
C
      MINCHO = LENAB + MXNNBS 
      MINEED = LREAD + MXNNBS
C
      NCHORD = MIN((LWRKX-MINEED)/MINCHO,MXCHRD)
C
Casm  Don't read more than 2 Gb (268435456 dw)
C
      MCHRD2 = 268435456 / MINCHO
      NCHORD = MIN(NCHORD,MCHRD2)
C
      IF (NCHORD .LT. 1) THEN 
         WRITE(LUPRI,*) 'Memory problem. NCHORD :',NCHORD
         CALL QUIT('Not enough  memory space in CC_CHOLESKY')
      END IF
C
      LCHO2 = LENAB*NCHORD + MXNNBS
      LCHO1 = LCHO2 + MXNNBS*NCHORD
      LNEED = LCHO1 + LREAD
C
      MAXSEL = MIN((LWORK-LNEED)/MXNNBS,MAXQUA)
C
Casm  Don't read more than 2 Gb (268435456 dw)
C
      MXSEL2 = 268435456 / MXNNBS
      MAXSEL = MIN(MAXSEL,MXSEL2)
C
      IF (MAXSEL .LT. 1) THEN 
         WRITE(LUPRI,*) 'Memory problem. MAXSEL :',MAXSEL
         CALL QUIT('Not enough  memory space in CC_CHOLESKY')
      END IF
C
      IF (IPRINT .GT. 2) THEN
         WRITE(LUPRI,'(/,A,I4)') 
     &         'Maximum number of vectors read        :',NCHORD
         WRITE(LUPRI,'(A,I4)') 
     &         'Maximum number of diagonals qualified :',MAXSEL
C
      END IF
      IF (IPRINT .GT. 0)
     &   WRITE(LUPRI,'(/,A,//)') 'End of Cholesky initializations'
      IF (IPRINT .EQ. 1) 
     &   WRITE(LUPRI,'(//4X,A,10X,A,5X,A,6X,A)')
     &         'Vectors','Distributions','Qualified','Max. diagonal'
      CALL FLSHFO(LUPRI)
C
C-------------------------------------------
C     The Cholesky passes should start here.
C-------------------------------------------
C
  100 CONTINUE
C
C     Find largest diagonal element.
C     ------------------------------
C
      CALL CC_DIASCR(DIAG)
C
      XMAX2 = -1.0D0
      DO ISYMAB = 1,NSYM
         XMAX1 = -1.0d0
         DO I = 1,NNBST(ISYMAB)
             XXX = DIAG(IDIAG(ISYMAB)+I)
             IF (XXX .GT. XMAX1) XMAX1 = XXX
         END DO
         DIAMIN(ISYMAB) = XMAX1 * SPAN
         IF (DIAMIN(ISYMAB) .LT. THRCOM) DIAMIN(ISYMAB) = THRCOM
         IF (XMAX1 .GT. XMAX2) XMAX2 = XMAX1
      END DO
C
      IF (ABS(XMAX2) .LT. THRCOM) THEN
c        CALL AROUND('Cholesky elements to discard')
c        WRITE(LUPRI,*) 'Number of Elements:',
c    *                  XDISC,XNDISC,(XDISC/XNDISC)*100.0D0
         CALL CHO_END(DIAG,NTOINT,XSEL0,YSEL0)
         RETURN
      END IF
C
      DO I = 1,NSYM
         NQUAL(I) = 0
      END DO
C
      NUMINT = 0
C
  110 CONTINUE
C
      CALL CC_DIAGLR0(XC,ISYMXC,ISHLA,ISHLB)
C
      IF (XC .EQ. ZERO) GOTO 120
C
      DO I = 1,NSYM
         ITMP(I) = NQUAL(I) + 1
      END DO
C
      THRCOM = THRDEF
C
C     Overall convergence check.
C     --------------------------
C
C
      DECO = .TRUE.
      IF (ABS(XC) .GE. DIAMIN(ISYMXC)) THEN
         DECO = .FALSE.
         GOTO 222
      END IF
C
      GOTO 120
  222 CONTINUE
C
C     Calculate integral distributions.
C     ---------------------------------
C
      SETUP = .FALSE.
C
      NUMINT = NUMINT + 1
      NTOINT = NTOINT + 1
C
      IF (DIASC1(ISHLA,ISHLB) .EQ. ZERO) THEN
          NWARN = NWARN + 1
          WRITE(LUPRI,*) 'WARNING: Something is wrong with this couple'
          WRITE(LUPRI,*) 'ISHLA, ISHLB, DIASC1(ISHLA,ISHLB)',
     &                    ISHLA, ISHLB, DIASC1(ISHLA,ISHLB)
      ENDIF
C
c     DO ISHLC = 1,MAXSHL
c        DO ISHLD = 1,MAXSHL
c           IF (DIASC1(ISHELC,ISHELD) .EQ. ZERO) THEN
c               IALSKP = IALSKP + 1
c           ELSE
c               IALDON = IALDON + 1
c           ENDIF
c        END DO
c     END DO
c     IALTOT = IALSKP + IALDON
C
c     IALTOT = IALDON + IALSKP
      IF (IPRINT. GT. 2 ) THEN
         WRITE(LUPRI,'(2A,2I6)') 'Entering new distribution ',
     &            'corresponding to shells :',ISHLA,ISHLB
         CALL FLSHFO(LUPRI)
      END IF
      
c     WRITE(LUPRI,'(3(A,I8,X))') 'CD distributions. Total :',IALTOT,
c    &                       'Calculated :',IALDON,'Skipped :',IALSKP
C
      T1 = SECOND()
      CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     *            NUMDIA,NUMDIB,JPRINT,SETUP,1)
      T1 = SECOND() - T1
      TIMINT = TIMINT + T1
C
      NEWDIS = .TRUE.
C
C     Select diagonals in this distribution.
C     --------------------------------------
C
      IALSEL = 0
      DO 300 NA = 1,NUMDIA
         INUM = 1
         IF (ISHLA .EQ. ISHLB) INUM = NA
         DO 400 NB = INUM,NUMDIB
C
            JSYMA  = ISAO(INDEXA(NA))
            JSYMB  = ISAO(INDEXB(NB))
            JSYMAB = MULD2H(JSYMB,JSYMA)
C
            A = INDEXA(NA) - IBAS(JSYMA)
            B = INDEXB(NB) - IBAS(JSYMB) 
C
            IF (JSYMA .EQ. JSYMB) THEN
               NAB = INDEX(A,B)
            ELSE IF (JSYMB .GT.JSYMA) THEN
               NAB = NBAS(JSYMA)*(B - 1) + A
            ELSE
               NAB = NBAS(JSYMB)*(A - 1) + B
            END IF
C
            NAB = IDIAG(JSYMAB) + IAODPK(JSYMA,JSYMB) + NAB


            IF ((DIAG(NAB).GE.DIAMIN(JSYMAB)) .AND.
     &          (DIAG(NAB).GE.THRCOM)) THEN
C
               IF (NQUAL(JSYMAB) .GE. MAXSEL) THEN
                  DECO = .TRUE.
                  GOTO 400
               ENDIF
C
               IALSEL = 1
C
               NQUAL(JSYMAB) = NQUAL(JSYMAB) + 1
               IQUAL = NQUAL(JSYMAB)
               IQUVAB(IQUAL,JSYMAB) = NAB
               IQUVEA(IQUAL,JSYMAB) = INDEXA(NA)
               IQUVEB(IQUAL,JSYMAB) = INDEXB(NB)
C
            ENDIF
C
  400    CONTINUE
  300 CONTINUE
C

      IF (IALSEL .NE. 1) THEN
         DECO = .TRUE.
         GOTO 120
      END IF
C
C
C     Read integrals from disk. 
C     -------------------------
C
      DO 350 ISYMAB = 1,NSYM
C
         I     = ITMP(ISYMAB)
         NUMAB = NQUAL(ISYMAB) - ITMP(ISYMAB) + 1
C
         IF (NUMAB .EQ. 0) GOTO 350
C
         KINT1 = 1
         KIREC = KINT1 + NNBST(ISYMAB)*NUMAB
         KEND1 = KIREC + (NBUFX(0) - 1)/IRAT + 1
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) CALL QUIT('Not enough space for CC_RDAON')
C
         NEWDIS = .TRUE.
         CALL CC_RDAON(WORK(KINT1),IQUVEA(I,ISYMAB),IQUVEB(I,ISYMAB),
     *                 NUMAB,ISYMAB,WORK(KIREC),WORK(KEND1),LWRK1)
C
         LUNIT  = 98
         NAME   = ROOT(1:7)//'SEL'
         LENGTH = NNBST(ISYMAB) * NUMAB
C
         IF (LENGTH .NE. 0) THEN
            IADR = MXNNBS * MAXSEL * (ISYMAB-1)
     &            + NNBST(ISYMAB) * (ITMP(ISYMAB)-1) + 1
            CALL PUTWA2(LUNIT,NAME,WORK(KINT1),IADR,LENGTH)
         END IF
C
  350 CONTINUE
C
      IF ((.NOT. DECO) .AND. (XMAX2 .LT. 1.0D-6)) GOTO 110
C
      NSEL = 0
      DO I = 1,NSYM
         NSEL = NSEL + NQUAL(I)
      END DO
      IF (NSEL .LT. 50) GOTO 110
C
  120 CONTINUE
C
      NINT2  = NINT2 + NUMINT
C
      IF (IPRINT .GT. 1) THEN
          WRITE(LUPRI,'(/,A,2I6)')
     &     'Number of calculated distributions :',NUMINT,NINT2
         WRITE(LUPRI,'(A,2X,8I4)')
     &     'Number of qualified diagonals      :',(NQUAL(I),I=1,NSYM)
         WRITE(LUPRI,'(A,D18.8)') 
     &     'Maximum diagonal among them        : ',XMAX2
         IF (IPRINT. GT. 2)
     &       WRITE(LUPRI,'(/,A,15X,A,15X,2(4X,A),6X,A)')
     &                 'ISAB    Iter','Treated','Conv','Neg','Max'
      END IF
C
C     Loop over symmetries in file to decompose.
C     ------------------------------------------
C
      ICHOV2 = ICHOV
      DIMAX  = ZERO
      DO 200 ISYMAB = 1,NSYM
C
         IF (NNBST(ISYMAB) .EQ. 0) GOTO 200
         IF (NQUAL(ISYMAB) .EQ. 0) GOTO 200
C
C        Allocations.
C        ------------
C
         IF (REDUCE) THEN
            LENAB = NREDUC(ISYMAB)
            LEMAB = LENAB
         ELSE IF (.NOT. COMP) THEN
            LENAB = NNBST(ISYMAB)
            LEMAB = NNBST(ISYMAB)
         ELSE
            LENAB = NNBST(ISYMAB) + (NNBST(ISYMAB)-1)/IRAT + 2
            LEMAB = LENAB
         END IF
         LENABT = LENAB + NNBST(ISYMAB)
C
C
C        LWRKX = LWORK - space for integrals - additional work space
C
         LWRKX = LWORK - NNBST(ISYMAB)*NQUAL(ISYMAB) - NNBST(ISYMAB)
C
         NCHORD = LWRKX / LENABT
         NCHORD = MIN(NCHORD,NUMCHO(ISYMAB))
         IF (NCHORD .EQ. 0) THEN
            NCHORD = MIN(10, LWRKX/LENABT)
         END IF
         NVEFIL = MAXLEN / LEMAB
         NCHORD = MIN(NCHORD,NVEFIL)
         NCHORD = MIN(NCHORD,MCHRD2)
         NCHORD = MIN(NCHORD,MXCHRD)
C
         LCHO2 = LENAB * NCHORD + NNBST(ISYMAB)
C
         KINT1 = 1
         KCHO1 = KINT1 + NNBST(ISYMAB)*NQUAL(ISYMAB)
         KCHO2 = KCHO1 + NNBST(ISYMAB)*NCHORD
         KEND1 = KCHO2 + LCHO2
         LWRK1 = LWORK - KEND1
c
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Insufficient work space in CC_CHOLESKY.1'
            WRITE(LUPRI,*) 'Needed and available memory',KEND1,LWORK
            CALL QUIT('Insufficient work space in CC_CHOLESKY.1')
         ENDIF
C
C        Read integrals from disk. 
C        -------------------------
C
         LUNIT  = 98
         NAME   = ROOT(1:7)//'SEL'
         LENGTH = NNBST(ISYMAB) * NQUAL(ISYMAB)
C
         IF (LENGTH .NE. 0) THEN
            IADR = MXNNBS * MAXSEL * (ISYMAB-1) + 1
            CALL GETWA2(LUNIT,NAME,WORK(KINT1),IADR,LENGTH)
         END IF
C

C        Batch over previous cholesky vectors.
C        -------------------------------------
C
         NBATV = (NUMCHO(ISYMAB)-1)/NCHORD + 1
C
         TIMRD = 0.0D0
         TIMPR = 0.0D0
         TIMDG = 0.0D0
C
         IBATV2 = 0
         DO IBATV = 1,NBATV
C
            IBATV1 = IBATV2 + 1
            IBATV2 = IBATV2 + NCHORD
            IF (IBATV2 .GT. NUMCHO(ISYMAB)) IBATV2 = NUMCHO(ISYMAB)
            NUMVEC = IBATV2 - IBATV1 + 1
            IF (NUMVEC .EQ. 0) GOTO 444
C
            TIMX = SECOND()
            CALL CC_GETCHO(WORK(KCHO1),NUMVEC,ISYMAB,
     &                     IBATV1,INDRED,WORK(KCHO2),LCHO2)
            TIMX  = SECOND() - TIMX
            TIMRD = TIMRD + TIMX
C
C           Prepare factor array.
C           ---------------------
C
            TIMX = SECOND()
            DO J = 1,NQUAL(ISYMAB)
               DO JJ = 1,NUMVEC
                  KOFF1 = KCHO1 
     *                  + NNBST(ISYMAB)*(JJ-1) 
     *                  + IQUVAB(J,ISYMAB) - IDIAG(ISYMAB) - 1
                  KOFF2 = KCHO2 + NUMVEC*(J-1) + JJ - 1
                  WORK(KOFF2) = WORK(KOFF1)
               ENDDO
            ENDDO
            TIMX  = SECOND() - TIMX
            TIMPR = TIMPR + TIMX
C
C           Subtraction.
C           ------------
C
            NTOTAB = MAX(NNBST(ISYMAB),1)
C
            TIMX = SECOND()
            CALL DGEMM('N','N',NNBST(ISYMAB),NQUAL(ISYMAB),NUMVEC,
     *                 -ONE,WORK(KCHO1),NTOTAB,WORK(KCHO2),
     *                 NUMVEC,ONE,WORK(KINT1),NTOTAB) 
            TIMX  = SECOND() - TIMX
            TIMDG = TIMDG + TIMX
C
         END DO
C
         IF (IPRINT .GT. 5) WRITE(LUPRI,'(A,3F12.2)') 
     &     'Time in reading, preparing, multiplying :',TIMRD,TIMPR,TIMDG
  444    CONTINUE
C
C        Decompose the vectors in core.
C        ------------------------------
C
         IDUMP = 0
         LAST  = .FALSE.
         DO 600 ICHO = 1,NQUAL(ISYMAB)
C
            XC = -ONE
            DO I = 1,NQUAL(ISYMAB)
               II = IQUVAB(I,ISYMAB) 
               IF (DIAG(II) .GT. XC) THEN
                  XC   = DIAG(II)
                  ICAB = II
                  IC   = I
               ENDIF
            ENDDO
C
            IF ((XC .LT. DIAMIN(ISYMAB)) .OR. (XC .LT. THRCOM)) THEN
               LAST = .TRUE.
               GOTO 333
            ENDIF    
C
            XD   = ONE/SQRT(XC)
            KOFF = KINT1 + NNBST(ISYMAB)*(IC-1)
            CALL DSCAL(NNBST(ISYMAB),XD,WORK(KOFF),1)
C
C           Screen the vector
C
ctst        CHOMAX = ZERO
ctst        DO I = 1,NNBST(ISYMAB)
ctst           KOFF = KINT1 + NNBST(ISYMAB)*(IC-1) + I - 1
ctst           IF (ABS(WORK(KOFF)) .GT. CHOMAX) 
ctst &             CHOMAX = ABS(WORK(KOFF))
ctst        END DO
C 
ctst        DO I = 1,NNBST(ISYMAB)
ctst           KOFF = KINT1 + NNBST(ISYMAB)*(IC-1) + I - 1
ctst           XNDISC = XNDISC + 1.0D0
ctst           IF (ABS(WORK(KOFF))*CHOMAX .LT. THRDEF)
ctst &            XDISC = XDISC + 1.0D0
ctst        END DO
C
ctst        XNUM   = DFLOAT(20 * NBAST)
ctst        THRSEL = THRCOM/CHOMAX/XNUM
            thrsel = zero
C
            XSEL1 = ZERO
            YSEL1 = ZERO
            DO I = 1,NNBST(ISYMAB)
               II = IDIAG(ISYMAB) + I
               KOFF = KINT1 + NNBST(ISYMAB)*(IC-1) + I - 1
               IF (DIAG(II) .EQ. ZERO) WORK(KOFF) = ZERO
c              IF (WORK(KOFF) .EQ. ZERO) THEN
c                 YSEL1 = YSEL1 + 1.0D0
c              ELSE IF (ABS(WORK(KOFF)) .LT. THRSEL) THEN
c                 WORK(KOFF) = ZERO
c                 XSEL1 = XSEL1 + 1.0D0
c              END IF
            END DO
C
            XZERO = NNBST(ISYMAB) - NREDUC(ISYMAB)
            YSEL1 = YSEL1 - XZERO
C
            XSEL0(ISYMAB) = XSEL0(ISYMAB) + XSEL1
            YSEL0(ISYMAB) = YSEL0(ISYMAB) + YSEL1
C
c           IF (IPRINT .GT. 15) WRITE(LUPRI,'(A,3X,2(A,F8.0,2X),A,I8)')
c    &                              'Elementents in vector',
c    &                              'Real 0 :',YSEL1,
c    &                              'Zeroed :',XSEL1,
c    &                              'Length :',LENAB
C
C           Update diagonal
C
            XM = ZERO
            DO I = 1,NNBST(ISYMAB)
               II = IDIAG(ISYMAB) + I
               KOFF = KINT1 + NNBST(ISYMAB)*(IC-1) + I - 1
               DIAG(II) = DIAG(II) - WORK(KOFF)*WORK(KOFF)
               IF (DIAG(II) .GT. XM) XM = DIAG(II)
            ENDDO
C
            IF (XM .GT. DIMAX) DIMAX = XM
            DIAMIN(ISYMAB) = XM * SPAN
            IF (DIAMIN(ISYMAB) .LT. THRCOM) DIAMIN(ISYMAB) = THRCOM
C
C           Analyze updated diagonal.
C           -------------------------
C
            OLDIAG     = DIAG(ICAB)
            DIAG(ICAB) = ZERO
            CALL CHO_CHKDIA(DIAG,DIAOLD,OLDIAG,XC,ICAB,ICHOV,
     &                      XM,ISYMAB,INDRED)
C
            DO I = 1,NQUAL(ISYMAB)
               II = IQUVAB(I,ISYMAB)
               IF (DIAG(II) .NE. ZERO) THEN
C
                  KOFF1 = KINT1 + NNBST(ISYMAB)*(IC-1)
                  KOFF2 = KINT1 + NNBST(ISYMAB)*(I-1) 
C
                  FACTOR = -WORK(KOFF1+II-IDIAG(ISYMAB)-1)
C
                  CALL DAXPY(NNBST(ISYMAB),FACTOR,
     *                       WORK(KOFF1),1,WORK(KOFF2),1)
               ENDIF
            ENDDO
C
C           Write cholesky vectors to disk.
C           -------------------------------
C
            IDUMP         = IDUMP + 1
            IDNTC1(IDUMP) = ICAB - IDIAG(ISYMAB)
C
            KOFF1 = KINT1 + NNBST(ISYMAB)*(IC-1)
            KOFF2 = KCHO1 + NNBST(ISYMAB)*(IDUMP-1)
            CALL DCOPY(NNBST(ISYMAB),WORK(KOFF1),1,WORK(KOFF2),1)
C
  333       CONTINUE
C
            IF ((IDUMP.EQ.NCHORD).OR.(ICHO.EQ.NQUAL(ISYMAB))
     *                           .OR.LAST) THEN
               CALL CC_PUTCHO(WORK(KCHO1),IDNTC1,NCHORD,IDUMP,
     &                        ISYMAB,INDRED,WORK(KCHO2),LCHO2)
               IF (LAST) GOTO 200
               IDUMP = 0
            ENDIF
C
  600    CONTINUE 
  200 CONTINUE

C
C     Statistics after treating distribution(s).
C     ------------------------------------------
C
      IF (IPRINT .GE. 1) THEN
         ICHOV2 = ICHOV - ICHOV2
         NQTOT  = 0
         DO JSYM = 1,NSYM
            NQTOT = NQTOT + NQUAL(JSYM)
         END DO
         IF (IPRINT .EQ. 1) THEN
            WRITE(LUPRI,'(I4,I9,8X,I3,I7,I13,D24.9)') ICHOV2,ICHOV,
     &                                 NUMINT,NINT2,NQTOT,DIMAX
         ELSE
            WRITE(LUPRI,'(/,A,2I6)')
     &           'Number of Cholesky vectors so far  : ',ICHOV2,ICHOV
            WRITE(LUPRI,'(A,D18.8)')
     &           'Max. diagonal after decomposition  : ',DIMAX
            WRITE(LUPRI,'(/,A,I5,A)') 'Histogram of diagonal elements:'
            CALL CC_ANADI(DIAG,.FALSE.)
         END IF
      END IF
C
C     Check if we can decompose more in the last distribution.
C     --------------------------------------------------------
C
      DO ISYMAB = 1,NSYM
         NQUAL(ISYMAB) = 0
         ITMP(ISYMAB)  = 1
      END DO
C
      XMAX2 = -1.0D0
      DO ISYMAB = 1,NSYM
         XMAX1 = -1.0d0
         DO I = 1,NNBST(ISYMAB)
             XXX = DIAG(IDIAG(ISYMAB)+I)
             IF (XXX .GT. XMAX1) XMAX1 = XXX
         END DO
         DIAMIN(ISYMAB) = XMAX1 * SPAN
         IF (DIAMIN(ISYMAB) .LT. THRCOM) DIAMIN(ISYMAB) = THRCOM
         IF (XMAX1 .GT. XMAX2) XMAX2 = XMAX1
      END DO
C
      DO 301 NA = 1,NUMDIA
         INUM = 1
         IF (ISHLA .EQ. ISHLB) INUM = NA
         DO 401 NB = INUM,NUMDIB
C
            JSYMA  = ISAO(INDEXA(NA))
            JSYMB  = ISAO(INDEXB(NB))
            JSYMAB = MULD2H(JSYMB,JSYMA)
C
            A = INDEXA(NA) - IBAS(JSYMA)
            B = INDEXB(NB) - IBAS(JSYMB) 
C
            IF (JSYMA .EQ. JSYMB) THEN
               NAB = INDEX(A,B)
            ELSE IF (JSYMB .GT.JSYMA) THEN
               NAB = NBAS(JSYMA)*(B - 1) + A
            ELSE
               NAB = NBAS(JSYMB)*(A - 1) + B
            END IF
C
            NAB = IDIAG(JSYMAB) + IAODPK(JSYMA,JSYMB) + NAB
C
            IF ((DIAG(NAB).GE.DIAMIN(JSYMAB)) .AND.
     &          (DIAG(NAB).GE.THRCOM)) THEN
C
               IF (NQUAL(JSYMAB) .GE. MAXSEL) THEN
                  DECO = .TRUE.
                  GOTO 401
               ENDIF
C
               NQUAL(JSYMAB) = NQUAL(JSYMAB) + 1
               IQUAL = NQUAL(JSYMAB)
               IQUVAB(IQUAL,JSYMAB) = NAB
               IQUVEA(IQUAL,JSYMAB) = INDEXA(NA)
               IQUVEB(IQUAL,JSYMAB) = INDEXB(NB)
C
            ENDIF
C
  401    CONTINUE
  301 CONTINUE
C
C     Should we do a backstep?
C
      NSEL = 0
      DO I = 1,NSYM
         NSEL = NSEL + NQUAL(I)
      END DO
      IF (NSEL .LT. 30) THEN
         IF (IPRINT .GT. 5) THEN
            WRITE(LUPRI,*) 
     &            'Not enough diagonals qualified at second pass'
            WRITE(LUPRI,*) 'NQUAL : ', (NQUAL(I),I=1,NSYM)
         END IF
         GOTO 100      
      END IF
C
C     Yes, we should
C
      NUMINT = 0
C
C     Read integrals from disk. 
C     -------------------------
C
      DO 351 ISYMAB = 1,NSYM
C
         I     = ITMP(ISYMAB)
         NUMAB = NQUAL(ISYMAB) - ITMP(ISYMAB) + 1
C
         IF (NUMAB .EQ. 0) GOTO 351
C
         KINT1 = 1
         KIREC = KINT1 + NNBST(ISYMAB)*NUMAB
         KEND1 = KIREC + (NBUFX(0) - 1)/IRAT + 1
         LWRK1 = LWORK - KEND1
C
         NEWDIS = .TRUE.
         CALL CC_RDAON(WORK(KINT1),IQUVEA(I,ISYMAB),IQUVEB(I,ISYMAB),
     *                 NUMAB,ISYMAB,WORK(KIREC),WORK(KEND1),LWRK1)
C
         LUNIT  = 98
         NAME   = ROOT(1:7)//'SEL'
         LENGTH = NNBST(ISYMAB) * NUMAB
C
         IF (LENGTH .NE. 0) THEN
            IADR = MXNNBS * MAXSEL * (ISYMAB-1)
     &            + NNBST(ISYMAB) * (ITMP(ISYMAB)-1) + 1
            CALL PUTWA2(LUNIT,NAME,WORK(KINT1),IADR,LENGTH)
         END IF
C
  351 CONTINUE
C
      GOTO 120
C
      END
C
C  /* Deck cc_rdao */
      SUBROUTINE CC_RDAO(XINT,IDELTA,IGAM,WORK,LWORK,IRECNR)
C
C     Written by Henrik Koch 25-Sep-1993
C
C     Purpose: Read destribution of AO integrals.
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION XINT(*),WORK(LWORK)
      DIMENSION IRECNR(*)
C
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccpack.h>
C
      ISYMD  = ISAO(IDELTA)
      ISYMG  = ISAO(IGAM)
      ISYMGD = MULD2H(ISYMD,ISYMG)
C
C----------------------------
C     Construct index arrays.
C----------------------------
C
      KADR2 = 1 
      KEND1 = KADR2 + (NBAST*NBAST + 1)/IRAT + 1
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         write(lupri,*) 'Insufficient space for allocation in CC_RDAO',
     &      KEND1,' > ',LWORK
         CALL QUIT('Insufficient space for allocation in CC_RDAO')
      END IF
C
      CALL CC_INX(WORK(KADR2),ISYMGD)
C
C--------------------
C     Construct XINT.
C--------------------
C
      CALL DZERO(XINT,NNBST(ISYMGD))
C
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     Buffer allocation
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      LRINT  = IRAT*LBUF
C
      KIBUF = KEND1
      KRBUF = KIBUF + (NIBUF*LBUF-1)/IRAT + 1
      KEND2 = KRBUF + LBUF
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) THEN
         write(lupri,*) 'Insufficient space for allocation in CC_RDAO',
     &      KEND2,' > ',LWORK
         CALL QUIT('Insufficient work space in CC_RDAO')
      ENDIF
C
      CALL CC_RDA1(XINT,WORK(KIBUF),WORK(KRBUF),IDELTA,IGAM,
     *             WORK(KADR2),IBIT1,IRECNR)
C
      RETURN
      END
C  /* Deck cc_rda1 */
      SUBROUTINE CC_RDA1(XINT,IBUF,RBUF,IDELTA,IGAM,KADR2,
     *                  IBIT1,IRECNR)
C
C     Written by Henrik Koch 25-Sep-1993
C
#include <implicit.h>
#include <ibtpar.h>
      DIMENSION XINT(*)
      DIMENSION KADR2(NBAST,NBAST)
      DIMENSION IBUF(LBUF,NIBUF),RBUF(LBUF)
      DIMENSION IRECNR(*)
      CHARACTER*8 FAODER
      LOGICAL OLDDX
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
Casm
      LOGICAL MULFIL
      PARAMETER (MULFIL = .FALSE.)
Casm
#include <ccinftap.h>
#include <ccorb.h>
#include <mxcent.h>
#include <eribuf.h>
#include <nuclei.h>
#include <inftap.h>
#include <chrnos.h>
#include <ibtfun.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      IF (NEWDIS) THEN
C
         NEWDIS = .FALSE.
C
         IF (LUINTR .LE. 0) THEN
            CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',
     &        'UNFORMATTED',IDUMMY,.FALSE.)
         END IF
         REWIND (LUINTR)
C
         DO 50 I = 1,NBUFX(0)
            READ(LUINTR) IRECNR(I)
   50    CONTINUE
C
      ENDIF
C
      IF (LUAORC(0) .LE. 0) THEN
            LBFINP = LBUF
C    
#if defined (VAR_INT64)
            LRECL = LBFINP + LBFINP + 1
#elif defined (SYS_NEC)
            LRECL = LBFINP + LBFINP/IRAT + 1
#else
            IF (NBASIS .LE. 255) THEN
               LRECL = 2*LBFINP + LBFINP + 1
            ELSE
               LRECL = 2*LBFINP + 2*LBFINP + 1
            END IF
#endif
            FAODER = 'AO2DIS'//CHRNOS(0)//CHRNOS(0)
            CALL GPOPEN(LUAORC(0),FAODER,'UNKNOWN','DIRECT',
     &           'UNFORMATTED',LRECL,OLDDX)
      END IF
C
      ICOUNT = 0
C
      IF (NIBUF .EQ. 1) THEN
C
         DO 100 J = 1,NBUFX(0)
C
            IF (IRECNR(J) .EQ. IDELTA) THEN
               ICOUNT = ICOUNT + 1
               NREC   = J
               IF (MULFIL) NREC = ICOUNT
               READ(LUAORC(0),ERR=2000,REC=NREC) RBUF,IBUF,NINT
               DO 110 I = 1,NINT
                  IR = IBTAND(IBTSHR(IBUF(I,1),2*NBITS),IBIT1)
                  IF (IR .EQ. IGAM) THEN
                     IP = IBTAND(       IBUF(I,1)         ,IBIT1)
                     IQ = IBTAND(IBTSHR(IBUF(I,1),  NBITS),IBIT1)
                     IADR = KADR2(IP,IQ) + 1
                     XINT(IADR) = RBUF(I)
                  ENDIF
  110          CONTINUE
            ENDIF
C
  100    CONTINUE
C
      ELSE
C
         DO 120 J = 1,NBUFX(0)
C
            IF (IRECNR(J) .EQ. IDELTA) THEN
               ICOUNT = ICOUNT + 1
               NREC   = J
               IF (MULFIL) NREC = ICOUNT
               READ(LUAORC(0),ERR=2000,REC=NREC) RBUF,IBUF,NINT
               DO 130 I = 1,NINT
                  IR = IBTAND(       IBUF(I,1)       ,IBIT1)
                  IF (IR .EQ. IGAM) THEN
                     IP = IBTAND(       IBUF(I,2)       ,IBIT1)
                     IQ = IBTAND(IBTSHR(IBUF(I,2),NBITS),IBIT1)
                     IADR = KADR2(IP,IQ) + 1
                     XINT(IADR) = RBUF(I)
                  ENDIF
  130          CONTINUE
            ENDIF
C
  120    CONTINUE
C
      ENDIF
C
      CALL GPCLOSE(LUAORC(0),'KEEP')
C
      RETURN
 2000 CALL QUIT('Error reading AOTWODIS')
      END
C  /* Deck ccrd_init */
      SUBROUTINE CC_INX(KADR2,ISYMAB)
C
C     asm 22-aug-1994
C
C     Purpose: Construct index arrays for CC_RDAO
C
#include <implicit.h>
C
      DIMENSION KADR2(NBAST,NBAST)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
         ICOUN2 = 0
         DO 210 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            IF (ISYMB .GT. ISYMA) THEN

               DO 220 B = 1,NBAS(ISYMB)
                  NB = IBAS(ISYMB) + B
C
                  DO 230 A = 1,NBAS(ISYMA)
                     NA = IBAS(ISYMA) + A
C
                     KADR2(NA,NB) = ICOUN2
                     KADR2(NB,NA) = ICOUN2
C
                     ICOUN2 = ICOUN2 + 1
C
  230             CONTINUE
  220          CONTINUE
C
            ELSE IF (ISYMA .EQ. ISYMB) THEN
C
               DO 240 B = 1,NBAS(ISYMB)
                  NB = IBAS(ISYMB) + B
C
                  DO 250 A = 1,B
                     NA = IBAS(ISYMA) + A
C
                     KADR2(NA,NB) = ICOUN2
                     KADR2(NB,NA) = ICOUN2
C
                     ICOUN2 = ICOUN2 + 1
C
  250             CONTINUE
  240          CONTINUE
C
            END IF
C
  210    CONTINUE
C
      RETURN
      END
C
C  /* Deck cc_gab */
      SUBROUTINE CC_GAB(DIAG,WORK,LWORK)
C
C     Written by Henrik Koch 15-juni-2000
C
C     Calculate diagonal elements. 
C
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
      LOGICAL   SETUP
      DIMENSION DIAG(*), WORK(LWORK)
      DIMENSION INDEXA(MXCORB), INDEXB(MXCORB)
      DIMENSION JNDEXA(MXCORB), JNDEXB(MXCORB)
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <distcl.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccfro.h>
#include <ccfop.h>
#include <ccsections.h>
#include <symsq.h>
#include <ccdeco.h>
#include <choles.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      DIACAL = .TRUE.
C
c     WRITE(LUPRI,'(A,I5)') 'Number of shells: ', MAXSHL
      DO ISHLA = 1,MAXSHL
         DO ISHLB = 1,ISHLA
C
            SETUP  = .FALSE.
            JPRINT = 0
            CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     *                  NUMDIA,NUMDIB,JPRINT,SETUP,2)
C
            DO I = 1,NUMDIA
               JNDEXA(INDEXA(I)) = I
            ENDDO
            DO I = 1,NUMDIB
               JNDEXB(INDEXB(I)) = I
            ENDDO
C
            KDIAG = 1
            KIREC = KDIAG + NUMDIA*NUMDIB
            KEND1 = KIREC + (NBUFX(0) - 1)/IRAT + 1
            LWRK1 = LWORK - KEND1 + 1
            IF (LWRK1 .LT. 0) CALL QUIT('Not enough space in CC_GAB')
C
            NEWDIS = .TRUE.
            CALL CC_RDAOX(WORK(KDIAG),NUMDIA,NUMDIB,JNDEXA,JNDEXB,
     *                    ISHLA,ISHLB,WORK(KEND1),LWRK1,WORK(KIREC))
C
            DO JA = 1,NUMDIA
               DO JB = 1,NUMDIB
C
                  IA = INDEXA(JA) 
                  IB = INDEXB(JB)
C
                  ISYMA  = ISAO(IA)
                  ISYMB  = ISAO(IB)
                  ISYMAB = MULD2H(ISYMA,ISYMB)
C
                  A = IA - IBAS(ISYMA)
                  B = IB - IBAS(ISYMB) 
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB = INDEX(A,B)
                  ELSE IF (ISYMB .GT.ISYMA) THEN
                     NAB = NBAS(ISYMA)*(B - 1) + A
                  ELSE
                     NAB = NBAS(ISYMB)*(A - 1) + B
                  END IF
C
                  NAB = IAODPK(ISYMA,ISYMB) + NAB
C
                  KOFF1 = IDIAG(ISYMAB) + NAB
                  KOFF2 = KDIAG + NUMDIA*(JB-1) + JA - 1 
C
                  DIAG(KOFF1) = WORK(KOFF2)
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
      DIACAL = .FALSE.
C
      RETURN
      END
C
C  /* Deck cc_chodrv */
      SUBROUTINE CC_CHODRV(WORK,LWORK)
C
C     Written by Alfredo Sanchez   20 july 2000
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <maxorb.h>
#include <iratdef.h>
C
      PARAMETER (ONTM10 = 1.0D-10)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CHARACTER*10 NAME
C
      DIMENSION TITMOL(12,2)
      DIMENSION WORK(LWORK)
      DIMENSION MMBST(8)
      LOGICAL   LDUM, FEXIST
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <choio.h>
#include <chorst.h>
C
#include <mxcent.h>
#include <eribuf.h>
C
C
C------------------------------------------
C     Suggest to use 64 bit integer version
C------------------------------------------
C
#if !defined (VAR_INT64)
      NINFO = NINFO + 1
      WRITE(LUPRI,'(///,80A1,/)') ('=',I=1,80)
      WRITE(LUPRI,'(38X,A,/)') 'INFO'
      WRITE(LUPRI,'(2X,A,A,/,2X,A,/)')
     &            'Full advantage of Cholesky decomposition based-',
     &            'methods is accomplished only ',
     &            'if 64-bit integers are used.'
      WRITE(LUPRI,'(2X,A,A,/,2X,A,A)')
     &            'Therefore, if you plan an intensive use of ',
     &            'Cholesky routines, rebuild Dalton ',
     &            'code with -DVAR_INT64 directive and ',
     &            '-i8 option.'
      WRITE(LUPRI,'(2X,A,A,/)') '(Remember to link to the proper ',
     &                          '64-bit libraries)'
      WRITE(LUPRI,'(2X,A,A)')
     &            'Note, however, that some other ',
     &            'functionalities do not work in this mode.'
      WRITE(LUPRI,'(/,80A1,///)') ('=',I=1,80)
#endif
C
C------------------------------------------------------
C     Get basis information and construct index arrays.
C------------------------------------------------------
C
      CALL GPOPEN(LUONEL,'AOONEINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUONEL)
C
      READ (LUONEL) 
      READ (LUONEL) NSYM,(NBAS(ISYM),ISYM=1,NSYM),POTNUC
C
      CALL GPCLOSE(LUONEL,'KEEP')
      CALL CHO_INDEX
      CALL FLSHFO(LUPRI)
C
C----------------
C     Read input.
C----------------
C
      CALL CC_CHOINP
      CALL FLSHFO(LUPRI)
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      LREDU = 1
      IF (REDUCE) THEN
         DO I = 1,NSYM
            IF (LREDU .LT. NNBST(I)) LREDU = NNBST(I)
         END DO
      END IF
C
      KDIAG = 1
      KOLD  = KDIAG + NDIAG
      KIND  = KOLD  + NDIAG
      KEND  = KIND  + NSYM * ((LREDU-1)/IRAT + 1)
      LWRK  = LWORK - KEND
      IF (LWRK .LT. 0) CALL QUIT('Not enough space in CC_CHODRV')
C
      CALL DZERO(WORK(KDIAG),NDIAG)
C
C----------------------
C     Do Cholesky loop.
C----------------------
C
      LURST = -1
      CALL GPOPEN(LURST,'CHOLESKY.RST','UNKNOWN','SEQUENTIAL',
     &            'UNFORMATTED',IDUM,.FALSE.)
C
      LUSEC = -1
      CALL GPOPEN(LUSEC,'CHOLESKY.SEC','UNKNOWN','SEQUENTIAL',
     &            'UNFORMATTED',IDUM,.FALSE.)
C
C-------------------------
C     Restart information.
C-------------------------
C
      IF (RSTCHO .OR. SKIP) THEN
C
         IF (SKIP) THEN
            WRITE (LUPRI,'(//A/A//)') 
     &         ' Cholesky decomposition skipped;',
     &         ' it is assumed that it was completed previously.'
C
         ELSE
C
            WRITE (LUPRI,'(//A/A/A//)') 
     &         '        -------------------------',
     &         '          Restarted calculation',
     &         '        -------------------------'
C
         END IF
C
C
         REWIND(LURST)
         READ(LURST,ERR=666,END=666) MSYM,THRCOM,LRDTOT,
     &               (MMBST(ISYM),ISYM=1,MSYM),
     &               (NUMCHO(ISYM),ISYM=1,MSYM),
     &               (IADRTO(ISYM),ISYM=1,MSYM),
     &               (NUMFIL(ISYM),ISYM=1,MSYM)
C
         IF (MSYM .NE. NSYM) THEN
            WRITE(LUPRI,*) 'Error in symmetry checking at restart'
            WRITE(LUPRI,*) 'NSYM expected :',NSYM
            WRITE(LUPRI,*) 'NSYM read     :',MSYM
            GOTO 666
         END IF
         DO I = 1,NSYM
            IF (MMBST(I) .NE. NNBST(I)) THEN
               WRITE(LUPRI,*) 'Error in basis checking at restart'
               WRITE(LUPRI,'(A,I2,1X,A,I10)') 'Symmetry',I,
     &                         'NNBST expected :',NNBST(I)
               WRITE(LUPRI,'(A,I2,1X,A,I10)') 'Symmetry',I,
     &                         'NNBST read     :',MMBST(I)
               GOTO 666
            END IF
         END DO
C
         READ(LURST,ERR=666,END=666) IDUM
C
         DO ISYM = 1,NSYM
            READ(LURST,ERR=666,END=666) 
     &          (LENCHO(I,ISYM),I=1,NUMCHO(ISYM))
         END DO
         DO ISYM = 1,NSYM
            READ(LURST,ERR=666,END=666) 
     &          (IDNTCH(I,ISYM),I=1,NUMCHO(ISYM))
         END DO
         DO ISYM = 1,NSYM
            READ(LURST,ERR=666,END=666) 
     &          (IADRLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
         END DO
         DO ISYM = 1,NSYM
            READ(LURST,ERR=666,END=666) 
     &          (IDNTLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
         END DO
C
         GOTO 777
  666    CONTINUE
C
C        Use secure file
C
         WRITE(LUPRI,'(//,A,/,A,//)') 'Error reading CHOLESKY.RST',
     &                   'Trying to use CHOLESKY.SEC instead'
         REWIND(LUSEC)
         READ(LUSEC) MSYM,THRCOM,LRDTOT,
     &               (MMBST(ISYM),ISYM=1,MSYM),
     &               (NUMCHO(ISYM),ISYM=1,MSYM),
     &               (IADRTO(ISYM),ISYM=1,MSYM),
     &               (NUMFIL(ISYM),ISYM=1,MSYM)
C
         IF (MSYM .NE. NSYM) THEN
            WRITE(LUPRI,*) 'Error in symmetry checking at restart'
            WRITE(LUPRI,*) 'NSYM expected :',NSYM
            WRITE(LUPRI,*) 'NSYM read     :',MSYM
            CALL QUIT('Error in symmetry checking at restart')
         END IF
         DO I = 1,NSYM
            IF (MMBST(I) .NE. NNBST(I)) THEN
               WRITE(LUPRI,*) 'Error in basis checking at restart'
               WRITE(LUPRI,'(A,I2,1X,A,I10)') 'Symmetry',ISYM,
     &                         'NNBST expected :',NNBST(ISYM)
               WRITE(LUPRI,'(A,I2,1X,A,I10)') 'Symmetry',ISYM,
     &                         'NNBST read     :',MMBST(ISYM)
               CALL QUIT('Error in basis checking at restart')
            END IF
         END DO
C
         READ(LUSEC) IDUM
C
         DO ISYM = 1,NSYM
            READ(LUSEC) (LENCHO(I,ISYM),I=1,NUMCHO(ISYM))
         END DO
         DO ISYM = 1,NSYM
            READ(LUSEC) (IDNTCH(I,ISYM),I=1,NUMCHO(ISYM))
         END DO
         DO ISYM = 1,NSYM
            READ(LUSEC) (IADRLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
         END DO
C
      END IF
C
  777 CONTINUE
      IF (SKIP) RETURN
C
C
         CALL CC_CHOLESKY(WORK(KDIAG),WORK(KOLD),WORK(KIND),
     &                    WORK(KEND),LWRK)
C
C
C     Poor man check result.
C
      IF (IPRCHO .GE. 0) THEN
         WRITE(LUPRI,'(//A//)') 'Poor man check of decomposition'
         CALL CHO_RESTART(WORK(KDIAG),WORK(KIND),WORK(KEND),LWRK)
      END IF
C
C     Write restart files
C
      CALL CHO_SAVE(WORK(KIND))
C
C     Close files
C
c     DO ISYM = 1,NSYM
c        CALL CHO_NAME(ISYM,LUCHO,NAME,.FALSE.)
c        CALL WCLOSE2(LUCHO,NAME,'KEEP')
c     END DO
C
      CALL GPCLOSE(LUDIAG,'KEEP')
      CALL GPCLOSE(LURST,'KEEP')
      CALL GPCLOSE(LUSEC,'KEEP')
      CALL GPCLOSE(LUAORC(0),'DELETE')
C
      LUSEL = 98
      CALL WCLOSE2(LUSEL,'CHOLES_SEL','DELETE')
C
      CALL GPINQ('AOTWOINT','EXIST',FEXIST)
      IF (FEXIST) THEN
         CALL SYSTEM('rm -f AOTWOINT')
      END IF
C
      WRITE(LUPRI,'(//)')
C
      RETURN
      END
C
C  /* Deck cc_choinp */
      SUBROUTINE CC_CHOINP
C
#include <implicit.h>
#include <gnrinf.h>
#include <maxorb.h>
#include <priunit.h>
#include <inftap.h>
#include <dummy.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <mxcent.h>
#include <eribuf.h>
#include <choscf.h>
#include <chorst.h>
#include <choio.h>
C
      PARAMETER (NTABLE = 15)
C
      LOGICAL SET,THRCHA,CHLBUF,SPANDI,THSCH1,THSCH2
      SAVE SET
C
      CHARACTER*7 WORD
      CHARACTER*7 TABLE(NTABLE)
C
      DATA SET /.FALSE./
      DATA TABLE /'.THINDI','.THSUDI','.NOSCDI','.NOCOMP','.THRCOM',
     &            '.LENBUF','.RSTDIA','.RSTCHO','.SPANDI','.DENDEC',
     &            '.THRSTR','.REDUCE','.SKIP  ','.PRINT ','.COMPLE'/
C
C
      IF (SET) RETURN
      SET = .TRUE.
C
C---------------------
C     Initializations.
C---------------------
C
      IPRCHO = MAX(1,IPRUSR)
C
      THINDI = 1.0D0 
      THSUDI = 1.0D+3
      THRDEF = 1.0D-8
      SPAN   = 1.0D-3
C
      LBUF = 250000
C
      RSTDIA = .FALSE.
      RSTCHO = .FALSE.
      REDUCE = .TRUE.
      COMP   = .FALSE.
      SCDIAG = .TRUE.
C
      THRCHA = .FALSE.
      THSCH1 = .FALSE.
      THSCH2 = .FALSE.
      CHLBUF = .FALSE.
      SPANDI = .FALSE.
C
      SKIP   = .FALSE.
C
C--------------------------------------
C     Initialize info for Cholesky SCF.
C--------------------------------------
C
      CCMODSK = .TRUE.
      LUCCMO  = 79
      FCCMO   = 'CHOCMO'
      CALL IZERO(NDVCS,8)
      THRDC   = 1.00D-12
C
C-------------------
C     Procces input.
C-------------------
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCMD)
C
  100 CONTINUE
C
      READ(LUCMD,'(A7)') WORD
      IF (WORD(1:7) .EQ. '*END OF') THEN
         GOTO 300
      ELSE IF (WORD(1:7) .NE. '**CHOLE') THEN
         GOTO 100
      END IF
C
  200 CONTINUE
C
      READ (LUCMD,'(A7)') WORD
      IF (WORD(1:1) .EQ. '!' .OR. WORD(1:1) .EQ. '#') THEN
         GOTO 200
      ELSE IF (WORD(1:1) .EQ. '*') THEN
         GOTO 300
      ELSE IF (WORD(1:1) .EQ. '.') THEN
         DO I = 1,NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) I
            END IF
         END DO
         WRITE(LUPRI,*) 'Illegal keyword in **CHOLESKY input: ',WORD
         CALL QUIT('Illegal keyword in **CHOLESKY input')
      ELSE
         WRITE(LUPRI,*) 'Illegal prompt in **CHOLESKY input'
         CALL QUIT('Illegal prompt in **CHOLESKY input')
      END IF

C
C-----------------------
C     Succesful keyword.
C-----------------------
C
C     thindi
C
    1 CONTINUE
         THSCH1 = .TRUE.
         READ(LUCMD,*) THINDI
      GOTO 200
C
C     thsudi
C
    2 CONTINUE
         THSCH2 = .TRUE.
         READ(LUCMD,*) THSUDI
      GOTO 200
C
C     noscdi
C
    3 CONTINUE
         SCDIAG = .FALSE.
      GOTO 200
C
C     nocomp
C
    4 CONTINUE
         COMP = .FALSE.
      GOTO 200
C
C     thrdef
C
    5 CONTINUE
         THRCHA = .TRUE.
         READ(LUCMD,*) THRDEF
      GOTO 200
C
C     lenbuf
C
    6 CONTINUE
         CHLBUF = .TRUE.
         READ(LUCMD,*) LBFINP
         IF (LBFINP .LT. 1)  CHLBUF = .FALSE.
ctst     LBFINP = MAX(LBUF,LBFINP)
      GOTO 200
C
C     rstdia
C
    7 CONTINUE
         RSTDIA = .TRUE.
      GOTO 200
C
C     rstcho
C
    8 CONTINUE
         RSTCHO = .TRUE.
      GOTO 200
C
C     spandi
C
    9 CONTINUE
         SPANDI = .TRUE.
         READ(LUCMD,*) SPAN
      GOTO 200
C
C     .DENDEC => decompose density matrix in DIIS SCF.
C                Read threshold for that decomposition.
C
   10 CONTINUE
         CCMODSK = .FALSE.
         READ(LUCMD,*) THRDC
      GOTO 200
C
C     thrstr
C
   11 CONTINUE
         THRSTA = .TRUE.
         READ(LUCMD,*) THRSTR
      GOTO 200
C
C     reduce
C
   12 CONTINUE
         REDUCE = .TRUE.
      GOTO 200
C
C     reduce
C
   13 CONTINUE
         SKIP = .TRUE.
      GOTO 200
C
C     print
C
   14 CONTINUE
         READ(LUCMD,*) IPRCHO
      GOTO 200
C
C     comple  
C
   15 CONTINUE
         REDUCE = .FALSE.
      GOTO 200
C
C------------------------------------------------------
C     The Cholesky input should now be completely read.
C------------------------------------------------------
C
  300 CONTINUE
      CALL GPCLOSE(LUCMD,'KEEP')
C
C
C---------------------
C     Inizializations.
C---------------------
C
      IF (CHLBUF) LBUF   = LBFINP
      IF (REDUCE) THEN
c        COMP   = .TRUE.
         SCDIAG = .TRUE.
      END IF
C
C
      IF (IPRCHO .LE. 0) GOTO 999
C
      CALL AROUND('Output from **CHOLESKY input')
C
      IF (REDUCE) THEN
          WRITE(LUPRI,'(A)') 'Cholesky vectors written in reduced set'
      ELSE
          WRITE(LUPRI,'(A)') 'Cholesky vectors written in complete set'
      END IF
      IF (SCDIAG) THEN
         IF (THSCH1) THEN
            WRITE(LUPRI,'(A,T60,1P,D14.4)')
     &  'First screening threshold on diagonal defined by user',THINDI
         ELSE
            WRITE(LUPRI,'(A,T60,1P,D14.4)')
     &  'Default first screening threshold on diagonal',THINDI
         END IF
         IF (THSCH2) THEN
            WRITE(LUPRI,'(A,T60,1P,D14.4)')
     &  'Second screening threshold on diagonal defined by user',THSUDI
         ELSE
            WRITE(LUPRI,'(A,T60,1P,D14.4)')
     &  'Default second screening threshold on diagonal',THSUDI
         END IF
      ELSE
         WRITE(LUPRI,'(A)') 'No screening on diagonal'
      END IF
      IF (.NOT. COMP) THEN
          WRITE(LUPRI,'(A)') 'Cholesky vectors are not compressed'
      ELSE
          WRITE(LUPRI,'(A)') 'Cholesky vectors are compressed'
      END IF
      IF (THRCHA) THEN
         WRITE(LUPRI,'(A,T60,1P,D14.4)')
     &  'Convergence threshold defined by user',THRDEF
      ELSE
         WRITE(LUPRI,'(A,T60,1P,D14.4)')
     &  'Default convergence threshold',THRDEF
      END IF
      IF (RSTCHO) THEN
         WRITE(LUPRI,'(A)') 'This is a restarted calculation'
      END IF
      IF (RSTDIA) THEN
         WRITE(LUPRI,'(A)') 'This calculation uses previous diagonal'
      END IF
      IF (CHLBUF) THEN
         WRITE(LUPRI,'(A,I8)') 'Buffer length:',LBUF
      ELSE
         WRITE(LUPRI,'(A,I8)') 'Default buffer length:',LBUF
      END IF
      IF (THRSTA) THEN
         WRITE(LUPRI,'(A,T60,1P,D14.4)')
     &  'Special restart threshold', THRSTR
      END IF
      IF (SPANDI) THEN
         WRITE(LUPRI,'(A,T60,1P,D14.4)') 'Span factor: ',SPAN
      ELSE
         WRITE(LUPRI,'(A,T60,1P,D14.4)') 'Default span factor: ',SPAN
      END IF
      WRITE(LUPRI,'(A,I5,I20)') 'File size (Gb; words) :',IGB,MAXLEN
      IF (.NOT. CCMODSK) THEN
         WRITE(LUPRI,'(A,T60,1P,D15.6)')
     &   'SCF density will be decomposed; threshold: ',THRDC
      ENDIF
C
  999 CONTINUE
      RETURN
      END
C
C  /* Deck cho_index */
      SUBROUTINE CHO_INDEX
C
C     Alfredo Sanchez.       21-Jul-2000
C
C     Construct part of commons required for Cholesky stuff.
C
#include <implicit.h>
C
      EXTERNAL INIDAT
C
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
#include <ccisao.h>
#include <ccdeco.h>
#include <choio.h>
#include <inftap.h>
#include <priunit.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C
C-------------
C     ccorb.h
C-------------
C
      N2BAST = 0
      NBAST  = 0
      ICOUNT = 0
C
      DO ISYM = 1,NSYM
C
         NBAST = NBAST + NBAS(ISYM)
         N2BAST = N2BAST + NBAS(ISYM)*NBAS(ISYM)
C
         IBAS(ISYM) = ICOUNT
         ICOUNT = ICOUNT + NBAS(ISYM)
C
      END DO
C
      NNBASX = (NBAST*(NBAST+1))/2
      N2BASX = NBAST*NBAST
C
C--------------
C     ccisao.h
C--------------
C
      IOFF = 0
      DO ISYM = 1,NSYM
         DO I = 1,NBAS(ISYM)
C
             IOFF = IOFF + 1
             ISAO(IOFF) = ISYM
C
         END DO
      END DO
C
C
C---------------
C     ccsdsym.h
C---------------
C
      ICOUNT = 0
      DO ISYM = 1,NSYM
C
         NNBST(ISYM) = 0
         N2BST(ISYM) = 0
C
         DO ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYM)
C
            N2BST(ISYM) = N2BST(ISYM) + NBAS(ISYMA)*NBAS(ISYMB)
C
            IF (ISYMB .GT. ISYMA) THEN
               NNBST(ISYM) = NNBST(ISYM) + NBAS(ISYMA)*NBAS(ISYMB)
            ELSE IF (ISYMB .EQ. ISYMA) THEN
               NNBST(ISYM) = NNBST(ISYM) + NBAS(ISYMA)*(NBAS(ISYMA)+1)/2
            ENDIF
C
         END DO
C
         I2BST(ISYM) = ICOUNT
         ICOUNT = ICOUNT + N2BST(ISYM)
C
      END DO
C
      DO ISYMD = 1,NSYM
C
         NDISAO(ISYMD) = 0
         ICOUNT = 0
C
         DO ISYMG = 1,NSYM
C
            ISYMAB = MULD2H(ISYMG,ISYMD)
C
            NDISAO(ISYMD) = NDISAO(ISYMD) + NNBST(ISYMAB)*NBAS(ISYMG)
C
            IDSAOG(ISYMG,ISYMD) = ICOUNT
            ICOUNT = ICOUNT + NNBST(ISYMAB)*NBAS(ISYMG)
C
         END DO
      END DO
C
C--------------
C     ccdeco.h
C--------------
C
      MAXDIA = 0
      ICOUNT = 0
      DO ISYM = 1,NSYM
         IF (.NOT. RSTCHO) NUMCHO(ISYM) = 0
         IF (MAXDIA .LT. NNBST(ISYM)) MAXDIA = NNBST(ISYM)
         IDIAG(ISYM)  = ICOUNT
         ICOUNT = ICOUNT + NNBST(ISYM)
      ENDDO
      ICHOV  = 0
      NDIAG  = ICOUNT
C
#if defined (VAR_INT64)
      IGB = 128
#else
      IGB = 16
#endif
#if defined (VAR_CHO2)
      IGB = 2
#endif
C
C     maxlen (wrd) = N Gb * 1073741824 wrd/Gb
C     (minus security margin, probably not needed)
C
      MAXLEN = 1073741823
      MAXLEN = MAXLEN * IGB
      MAXLEN = MAXLEN - MAXDIA
C
C-------------
C     symsq.h
C-------------
C
      DO ISYMAB = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         DO ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            IAODIS(ISYMA,ISYMB) = ICOUN1
            IAODPK(ISYMA,ISYMB) = ICOUN2
            IAODPK(ISYMB,ISYMA) = ICOUN2
C
            ICOUN1 = ICOUN1 + NBAS(ISYMA)*NBAS(ISYMB)
            IF (ISYMB .GT. ISYMA) THEN
               ICOUN2 = ICOUN2 + NBAS(ISYMA)*NBAS(ISYMB)
            ELSE IF (ISYMAB .EQ. 1) THEN
               ICOUN2 = ICOUN2 + NBAS(ISYMB)*(NBAS(ISYMB)+1)/2
            ENDIF
C
         END DO
      END DO
C
      DO ISYMAB = 1,NSYM
         DO ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            DO B = 1,NBAS(ISYMB)
               DO A = 1,NBAS(ISYMA)
C
                  NABSQ = I2BST(ISYMAB) + IAODIS(ISYMA,ISYMB)
     *                  + NBAS(ISYMA)*(B - 1) + A
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB = INDEX(A,B)
                  ELSE IF (ISYMB .GT.ISYMA) THEN
                     NAB = NBAS(ISYMA)*(B - 1) + A
                  ELSE
                     NAB = NBAS(ISYMB)*(A - 1) + B
                  END IF
C
                  NABPK = IAODPK(ISYMA,ISYMB) + NAB
C
                  IADRPK(NABSQ) = NABPK
C
               END DO
            END DO
         END DO
      END DO
C
C
      IF (IPRCHO .GT. 5) THEN
C
         CALL AROUND('Information from CHO_INDEX')
C
	 WRITE(LUPRI,1) 'NBAST  :',NBAST
	 WRITE(LUPRI,1) 'N2BAST :',N2BAST
	 WRITE(LUPRI,1) 'N2BASX :',N2BASX
	 WRITE(LUPRI,1) 'NNBASX :',NNBASX
C
         WRITE(LUPRI,*)
         WRITE(LUPRI,1) 'NBAS   :',(NBAS(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'IBAS   :',(IBAS(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'NNBST  :',(NNBST(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'N2BST  :',(N2BST(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NDISAO :',(NDISAO(I), I=1,NSYM)
C
         WRITE(LUPRI,*)
         DO I = 1,NSYM
            WRITE(LUPRI,1) 'IDSAOG :',(IDSAOG(I,J), J=1,NSYM)
         END DO
C
         WRITE(LUPRI,*)
         DO I = 1,NSYM
            WRITE(LUPRI,1) 'IAODPK :',(IAODIS(I,J), J=1,NSYM)
         END DO
C
         WRITE(LUPRI,*)
         DO I = 1,NSYM
            WRITE(LUPRI,1) 'IAODIS :',(IAODIS(I,J), J=1,NSYM)
         END DO
C
      END IF
C
    1 FORMAT(3X,A8,8I8)
C
      RETURN
      END
C  /* Deck cc_diascr */
      SUBROUTINE CC_DIASCR(DIAG)
C
C     Written by Henrik Koch 15-juni-2000
C
C     Find largest diagonal element in a shell. 
C
C
#include <implicit.h>
#include <maxorb.h>
      DIMENSION DIAG(*)
#include <ccisao.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <symsq.h>
#include <blocks.h>
#include <ccdeco2.h>
#include <priunit.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      DO ISHLA = 1,MAXSHL
         DO ISHLB = 1,ISHLA
C
            XTMP = -1.0D0
            DO 100 IA = 1,NBAST
               DO 200 IB = 1,NBAST
C
                  IF (INAOSH(IA).NE.ISHLA) GOTO 100
                  IF (INAOSH(IB).NE.ISHLB) GOTO 200
C
                  ISYMA  = ISAO(IA)
                  ISYMB  = ISAO(IB)
                  ISYMAB = MULD2H(ISYMA,ISYMB)
C
                  A = IA - IBAS(ISYMA)
                  B = IB - IBAS(ISYMB) 
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB = INDEX(A,B)
                  ELSE IF (ISYMB .GT.ISYMA) THEN
                     NAB = NBAS(ISYMA)*(B - 1) + A
                  ELSE
                     NAB = NBAS(ISYMB)*(A - 1) + B
                  END IF
C
                  NAB = IDIAG(ISYMAB) + IAODPK(ISYMA,ISYMB) + NAB
C
                  IF (XTMP .LT. DIAG(NAB)) THEN
                     XTMP = DIAG(NAB)
                     ITMP = ISYMAB
                  END IF
C
  200          CONTINUE
  100       CONTINUE
C
            DIASCR(ISHLA,ISHLB) = XTMP
            DIASCR(ISHLB,ISHLA) = XTMP
            ISYSCR(ISHLA,ISHLB) = ITMP
            ISYSCR(ISHLB,ISHLA) = ITMP
C
            DIASC1(ISHLA,ISHLB) = XTMP
            DIASC1(ISHLB,ISHLA) = XTMP
C
         ENDDO
      ENDDO
C
c     WRITE(LUPRI,'(//,A,/)') 'DIASCR'
c     DO I = 1,MAXSHL
c        WRITE(LUPRI,'(4D18.6)') (DIASCR(I,J),J=1,MAXSHL)
c     END DO
C
      RETURN
      END
C  /* Deck cc_choscr */
      SUBROUTINE CC_CHOSCR(NDIM,CHOVEC,ISYMAB,isha,ishb)
C
C     Written by Alfredo Sanchez  27-september-2000
C
C     Zero out elements in Cholesky vector screened.
C
C
#include <implicit.h>
#include <maxorb.h>
#include <ccisao.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <symsq.h>
#include <blocks.h>
#include <choskp.h>
C
      PARAMETER (ZERO = 0.0D0)
C
      DIMENSION CHOVEC(NDIM)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C
      CALL QUIT('non-implemented choscr called')

#ifdef NONIMPLEMENTED
c
      ISYMOP = 1
      JSYMCD = MULD2H(ISYMOP,ISYMAB)
C
      write(6,*) 'choscr called',jsymcd
C
      DO 100 ISHC = 1,MAXSHL
         DO 200 ISHD = 1,MAXSHL
C
            IF (SCRCD(ISHC,ISHD) .NE. 1) GOTO 200
C
            DO 300 IC = 1,NBAST
C
               IF (INAOSH(IC) .NE. ISHC) GOTO 300
C
               DO 400 ID = 1,NBAST
C
                  IF (INAOSH(ID) .NE. ISHD) GOTO 400
C
                  ISYMC  = ISAO(IC)
                  ISYMD  = ISAO(ID)
                  ISYMCD = MULD2H(ISYMC,ISYMD)
C
c                 IF (ISYMCD .NE. JSYMCD) GOTO 400
C
                  C = IC - IBAS(ISYMC)
                  D = ID - IBAS(ISYMD) 
C
                  IF (ISYMC .EQ. ISYMD) THEN
                     NCD = INDEX(C,D)
                  ELSE IF (ISYMD .GT.ISYMC) THEN
                     NCD = NBAS(ISYMC)*(D - 1) + C
                  ELSE
                     NCD = NBAS(ISYMC)*(C - 1) + D
                  END IF
C
                  NCD = IAODPK(ISYMC,ISYMD) + NCD
C
c                 CHOVEC(NCD) = ZERO
C
                  write(6,99)'ishc,ishd,c,d,isymc,isymd,ncd,cho',
     &                        ishc,ishd,c,d,isymc,isymd,ncd,chovec(ncd)
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
C
   99 format(a,7i5,d20.8)
#endif
      RETURN
      END
C
C/* Deck cho_restart */
C
      SUBROUTINE CHO_RESTART(DIAG,INDRED,WORK,LWORK)
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <chorst.h>
#include <chodbg2.h>
#include <choio.h>
C
      DIMENSION DIAG(*), WORK(LWORK)
      DIMENSION INDRED(*)
C
C
      IF (IPRCHO .GT. 0) CALL AROUND('Summary of Cholesky restart')
C
      NTOT = 0
      NSAV = 0
C
C     Read and analyze diagonal.
C     --------------------------
C
      REWIND(LUDIAG)
      READ(LUDIAG) (DIAG(I),I = 1,NDIAG)
C
      IF (IPRCHO .GT. 2) THEN
         WRITE(LUPRI,'(//,A,I10)') 'Number of diagonal elements :',NDIAG
         CALL FLSHFO(LUPRI)
C
         WRITE(LUPRI,'(//,A,//)') 'Analysis of read diagonal'
         CALL CC_ANADI(DIAG,.TRUE.)
      END IF
C
      DO ISYMAB = 1,NSYM
C
C        Initializations.
C        ----------------
C
         NDIM = NNBST(ISYMAB)
         NVEC = NUMCHO(ISYMAB)
         IF (NVEC .EQ. 0) GOTO 100
C
C        Allocate memory.
C        ----------------
C
         KDIAG = 1
         KEND1 = KDIAG + NDIM
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0)
     &   CALL QUIT('Not enough space in CHO_RESTART.1')
C
         CALL DZERO(WORK(KDIAG),NDIM)
C
C        Read and substract previous vectors.
C        ------------------------------------
C
         IF (REDUCE) THEN
            LWRKX = LWRK1
            LENAB = NDIM + NREDUC(ISYMAB)
         ELSE IF (COMP) THEN
            LWRKX = LWRK1 - NDIM
            LENAB = NDIM  + (NDIM-1)/IRAT + 2 + NDIM
         ELSE
            LWRKX = LWRK1 - 1
            LENAB = NDIM
         END IF
C
Casm     Apparently, it is not possible to read more than 2 Gb (268435456 dw)
C
         NRDMX  = LWRKX / LENAB
         NRDMX1 = 268435456 / LENAB
         NRDMX  = MIN(NRDMX,NRDMX1)
         IF (NRDMX .LT. 1)
     &      CALL QUIT('Allocation failure in CHO_RESTART')
C
         NBAT = (NVEC-1)/NRDMX + 1
C
         ILAST = 0
         DO IBAT = 1,NBAT
C
            IFIRST = ILAST + 1
            ILAST  = ILAST + NRDMX
            IF (ILAST .GT. NVEC) ILAST = NVEC
            NREAD = ILAST - IFIRST + 1
C
C           Read them in, unpack and contruct new diagonal.
C           -----------------------------------------------
C
            LCHO1 = NREAD * NDIM
            IF (REDUCE) THEN
               LCHO2 = NREAD * NREDUC(ISYMAB)
            ELSE IF (COMP) THEN
               LCHO2 = NREAD * (NDIM  + (NDIM-1)/IRAT + 2) + NDIM
            ELSE
               LCHO2 = 1
            END IF
C
            KCHO1 = KEND1
            KCHO2 = KCHO1 + LCHO1
            KEND2 = KCHO2 + LCHO2
C
            LWRK2 = LWORK - KEND2
            IF (LWRK2 .LT. 0)
     &         CALL QUIT('Not enough space in CHO_RESTART.2')
C
            CALL CC_GETCHO(WORK(KCHO1),NREAD,ISYMAB,
     &                        IFIRST,INDRED,WORK(KCHO2),LCHO2)
C
            DO I = IFIRST,ILAST
               II = I - IFIRST + 1
               DO J = 1,NDIM
                  KOFF1 = KDIAG + J - 1
                  KOFF2 = KCHO1 + NDIM*(II-1) + J - 1
                  WORK(KOFF1) = WORK(KOFF1) + WORK(KOFF2)*WORK(KOFF2)
               END DO
            END DO
C
         END DO
C
C        Compare.
C        --------
C
         ERRMAX = ZERO
         ERRMX1 = ZERO
         DO I = 1,NDIM
            KOFF1 = IDIAG(ISYMAB) + I
            KOFF2 = KDIAG + I - 1
C
            ERRO1 = ABS(DIAG(KOFF1)-WORK(KOFF2))
            IF (ERRO1 .GT. ERRMX1) THEN 
               ERRMX1 = ERRO1
               EXACT1 = DIAG(KOFF1)
            END IF
C
         END DO
C
C        Make the substraction.
C        ----------------------
C
         KOFF = IDIAG(ISYMAB) + 1
         CALL DAXPY(NDIM,-ONE,WORK(KDIAG),1,DIAG(KOFF),1)
C
C        Check and zero out converged diagonal elements.
C        -----------------------------------------------
C
         IF (IPRCHO .GE. 0) WRITE(LUPRI,'(/,A,I3)') 
     &                      'Check information for symmetry',ISYMAB
C
         DO I = 1,NUMCHO(ISYMAB)
            KOFF = IDNTCH(I,ISYMAB) + IDIAG(ISYMAB)
            DIAG(KOFF) = ZERO
         END DO
C
         XM = ZERO
         YM = 1.0D20
         DO I = 1,NDIM
            KOFF = IDIAG(ISYMAB) + I
            IF (DIAG(KOFF) .GT. XM) XM = DIAG(KOFF)
            IF (DIAG(KOFF) .LT. YM) YM = DIAG(KOFF)
         END DO
C
         INEG = 0
         DO I = 1,NNBST(ISYMAB)
C
            KOFF1  = IDIAG(ISYMAB) + I
            DIATMP = DIAG(KOFF1)
C
            IF (DIATMP .LT. -1.0D-10) WRITE(LUPRI,'(A,I8,D15.5)')
     &         'Negative diagonal:',I,DIATMP
C
            IF (DIATMP .LT. -1.0d-13) THEN
               DIAG(KOFF1) = ZERO
               INEG = INEG + 1
            END IF
C
         END DO
C
         DO I = 1,NNBST(ISYMAB)
            KOFF1 = IDIAG(ISYMAB) + I
            DIAKO = SQRT(ABS(DIAG(KOFF1)) * XM)*THSUDI
            IF (DIAKO .LT. THRCOM) THEN
                  IF (SCDIAG) DIAG(KOFF1) = ZERO
            END IF
         END DO
C
         ICONV = 0
         DO I = 1,NNBST(ISYMAB)
            KOFF1 = IDIAG(ISYMAB) + I
            DIAKO = ABS(DIAG(KOFF1))
            IF (DIAKO .LT. THRCOM) ICONV = ICONV + 1
         END DO
         IUNCO = NDIM - ICONV
C
         IF (IPRCHO .GE. 0) THEN
            WRITE(LUPRI,'(/,A,I8,/)') 'Number of vectors      :',NVEC
            WRITE(LUPRI,'(2(A,D18.8,/))') 'Maximum diagonal       :',XM,
     &                             'Minimum diagonal       :',YM
            WRITE(LUPRI,'(A,2D18.8,/)')'Maximum absolute error :',
     &                             ERRMX1,EXACT1
            WRITE(LUPRI,'(3(A,I8,/))') 'Converged   diagonals :',ICONV,
     &                          'Zeroed neg  diagonals :',INEG,
     &                          'Unconverged diagonals :',IUNCO
            CALL FLSHFO(LUPRI)
         END IF
C
  100    CONTINUE
C
      END DO
C
      IF (IPRCHO .GT. 1) THEN
         WRITE(LUPRI,'(//,2(9X,A,/))') 'Histogram of diagonal elements',
     &                          '------------------------------'
         CALL CC_ANADI(DIAG,.TRUE.)
C
         WRITE(LUPRI,'(//,2(9X,A,/))') 'Number of vectors',
     &                          '-----------------'
         WRITE(LUPRI,'(6X,A,8X,A,7X,A)') 'ISAB','Maximum','Actual'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(7X,I2,2I16)') ISYM,NNBST(ISYM),NUMCHO(ISYM)
         END DO
C
         WRITE(LUPRI,*)
         CALL FLSHFO(LUPRI)
      END IF
C
      RETURN 
      END
C
      SUBROUTINE CHO_RST1(DIAG,NDIM,SCR,LSCR,NTOT,NSAV,SUMRST,
     &                    VECNRM)
C
C     Written by Alfredo Sanchez   8 October 2000
C
C     Update Cholesky vector with previous ones
C
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <chodbg2.h>
C
      DIMENSION DIAG(NDIM),SCR(LSCR)
C
C
      NUMEL = INT(SCR(1))
C
      KCHISE = 2
      KCHSAV = KCHISE + (NUMEL - 1)/IRAT + 1
C
      CALL CHO_RST2(DIAG,SCR(KCHISE),SCR(KCHSAV),NUMEL,NDIM,SUMRST,
     &              VECNRM)
C
      NTOT = NTOT + NDIM
      NSAV = NSAV + NUMEL
C
      RETURN
      END
C
      SUBROUTINE CHO_RST2(DIAG,ISELEC,CHOSEL,NUMEL,NDIM,SUMRST,VECNRM)
C
C     Written by Alfredo Sanchez   13 de julio del 2000
C
C     Select the relevant Cholesky elements
C
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <chodbg2.h>
C
      DIMENSION DIAG(NDIM),CHOSEL(NUMEL),ISELEC(NUMEL)
C
C
      SUMTMP = ZERO
      DO I = 1,NUMEL
         J = ISELEC(I)
         XTEMP = CHOSEL(I)*CHOSEL(I)
         DIAG(J) = DIAG(J) + XTEMP
         SUMTMP = SUMTMP + CHOSEL(I)
      END DO
      SUMRST = SUMRST + SUMTMP*SUMTMP
C
      VECNRM = DDOT(NUMEL,CHOSEL,1,CHOSEL,1)
C
      RETURN
      END
C
C/* Deck cho_name */
      SUBROUTINE CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,TASK)
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccdeco.h>
#include <choio.h>

C
      CHARACTER*10 NAME
      CHARACTER*7  ROOT
      CHARACTER*1  CSYM,CFIL,UNDERL
      DATA UNDERL  /'_'/
      DATA ROOT    /'CHOLES_'/
C
      CHARACTER*3 TASK
      INTEGER     ITEMP
C
      SAVE NTOTFI
C
C
      IF (TASK .EQ. 'INI') THEN
C
C        Initialize iadrlu. It is exactly recalculated when needed
C
         DO J = 0,9
            ITEMP = MAXLEN*J + 1
             DO I = 1,NSYM
                IADRLU(J,I) = ITEMP
             END DO
         END DO
C
         J = 0
         WRITE(CFIL,'(I1)') J
         NTOTFI = 0
C
         DO I = 1,NSYM
            WRITE(CSYM,'(I1)') I
            LUCHO = 70 + I
            NAME  = ROOT(1:7)//CSYM//UNDERL//CFIL
            CALL WOPEN2(LUCHO,NAME,64,0)
C
            IADRTO(I)   = 1
            IADRLU(J,I) = 1
            IDNTLU(J,I) = LUCHO
            NTOTFI      = NTOTFI + 1
            NUMFIL(I)   = 1
         END DO
C
      ELSE IF (TASK .EQ. 'OLD') THEN
C
         J = IFILE
         I = ISYMAB
         WRITE(CFIL,'(I1)') J
         WRITE(CSYM,'(I1)') I
C
         NAME  = ROOT(1:7)//CSYM//UNDERL//CFIL
         LUCHO = IDNTLU(J,I)
C
      ELSE IF (TASK .EQ. 'NEW') THEN
C
         J = IFILE
         I = ISYMAB
         WRITE(CFIL,'(I1)') J
         WRITE(CSYM,'(I1)') I
C
         NTOTFI = NTOTFI + 1
         LUCHO  = 70 + NTOTFI
         NUMFIL(I)   = NUMFIL(I) + 1
         IDNTLU(J,I) = LUCHO
C
         NAME  = ROOT(1:7)//CSYM//UNDERL//CFIL
         CALL WOPEN2(LUCHO,NAME,64,0)
C
      ELSE IF (TASK .EQ. 'RST') THEN
C
         NTOTFI = 0
         DO I = 1,NSYM
            J = -1
            DO K = 1,NUMFIL(I)
               J = J + 1
C
               WRITE(CFIL,'(I1)') J
               WRITE(CSYM,'(I1)') I
               NTOTFI = NTOTFI + 1
C
               NAME  = ROOT(1:7)//CSYM//UNDERL//CFIL
               LUCHO = IDNTLU(J,I)
               CALL WOPEN2(LUCHO,NAME,64,0)
C
            END DO
         END DO
C
      ELSE
C
         WRITE(LUPRI,'(//3A/A)')
     &      ' Task ',TASK,' not implemented in CHO_NAME',
     &      ' Program will stop.'
         CALL QUIT('Error in CHO_NAME')
C
      END IF
C
      RETURN
      END
C
C  /* Deck cho_read */
      SUBROUTINE CHO_READ(CHOVEC,ISTART,NR,IND1,IND2,ISYMAB,IOPT,
     &                    WORK,LWORK)
C
C     Written by Alfredo Sanchez   13 de julio del 2001
C     - IOPT=3 introduced, tbp, May 2002.
C     
C
C     Read NR Cholesky vectors
C
C     IOPT = 1 : CHOVEC     only non zero ab elements
C                IND1       alpha index
C                IND2       beta  index
C
C     IOPT = 2:  CHOVEC     squared (with zeros), (a,b) index
C                IND1,IND2  dummy
C
C     IOPT = 3:  CHOVEC     packed  (with zeros), a.le.b
C                IND1,IND2  dummy
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <chodbg2.h>
C
      CHARACTER*10 NAME
      LOGICAL LOG,test
      DATA LOG /.FALSE./
      save nalf
      data nalf /0/
C
Corig DIMENSION CHOVEC(N2BST(ISYMAB),NR),WORK(LWORK)
      DIMENSION CHOVEC(NR,*),WORK(LWORK)  ! just for it to compile...
      DIMENSION IND1(*),IND2(*)
C
      INTEGER MDIM,JSTART
C
      WRITE(lupri,*) 'Very old version of CHO_READ called'
      WRITE(lupri,*) 'It does not work any more'
      CALL QUIT('Old CHO_READ called')
C
      nalf = nalf + 1
      NDSQ = N2BST(ISYMAB)
      NDIM = NNBST(ISYMAB)
      CALL DZERO(CHOVEC,NDSQ*NR)
C
      ILAST = ISTART - 1
      IEND  = ISTART + NR - 1
C
      IF (IOPT .EQ. 1) THEN
         CALL QUIT('IOPT = 1 not yet implemented')
      END IF
C
      IF ((IOPT.EQ.2) .OR. (IOPT.EQ.3)) THEN
C
         LENADD = 0
C
C        Find out where the first vector is.
C
C_common DO I = 0,NUMFIL(ISYMAB)-1
C_common    IREF = ISTRFI(I,ISYMAB)
C_common    IF (ISTART .GE. IREF) IFILE = I
C_common END DO
C
         IR = 0
  200    CONTINUE
            IFIRST = ILAST + 1
C
C           Get unit number and name.
C
C_old       CALL CHONAMEOLD(ISYMAB,IFILE,LUCHO,NAME,LOG,LOG,LOG)
C_notsoold  CALL CHONAME(ISYMAB,LUCHO,NAME,LOG)
            CALL WOPEN2(LUCHO,NAME,64,0)
C
C_16Gb      IADR   = INDCHO(IFIRST,ISYMAB)
            IF (REDUCE) THEN
               MDIM   = NREDUC(ISYMAB)
               JSTART = IFIRST - 1
               IADR = MDIM * JSTART + 1
            ELSE IF (.NOT. COMP) THEN
               MDIM   = NNBST(ISYMAB)
               JSTART = IFIRST - 1
               IADR = MDIM * JSTART + 1
            ELSE
               IADR = 1
               DO I = 1,IFIRST-1
                  MDIM = LENCHO(I,ISYMAB)
                  IADR = IADR + MDIM
               END DO
            END IF
            LENTOT = LENCHO(IFIRST,ISYMAB)
            LENADD = LENADD + LENTOT
            IREAD  = IFIRST
C
            IF (LENTOT .GT. LWORK) THEN
               WRITE(LUPRI,*) 'Not enough space in cho_read'
               WRITE(LUPRI,'(A,I10)') 'Needed    :',LENTOT
               WRITE(LUPRI,'(A,I10)') 'Available :',LWORK
               CALL QUIT('Not enough space in cho_read')
            END IF
  210       CONTINUE
C
               IF (IREAD .EQ. IEND) GOTO 220
C
C              See how many vectors can be read.
C
C_16Gb         INCLEN = LENCHO(IREAD+1,ISYMAB)
C_16Gb         LETEMP = LENADD + INCLEN
C_16Gb         IF (LETEMP .GT. MAXLEN) THEN
C_16Gb            IFILE  = IFILE + 1
C_16Gb            LENADD = 0
C_16Gb            GOTO 220
C_16Gb         END IF
C
!              LENTMP = LENTOT + INCLEN
! INCLEN undefined
               IF (LENTMP .LT. LWORK) THEN
                  IREAD  = IREAD + 1
                  LENTOT = LENTMP
!                 LENADD = LETEMP
! LETEMP undefined
                  GOTO 210
               END IF
C
  220       CONTINUE
            ILAST = IREAD
C
C           Read the vectors and move them to CHOVEC
C
            CALL GETWA2(LUCHO,NAME,WORK(1),IADR,LENTOT)
            CALL WCLOSE2(LUCHO,NAME,'KEEP')
C
            KOFVEC = 1
            DO I = IFIRST,ILAST
C
               IR = IR + 1
C
               LENG  = LENCHO(I,ISYMAB)
               NUMEL = INT(WORK(KOFVEC))
C
               KOFIND = KOFVEC + 1
               KOFVE2 = KOFIND + (NUMEL - 1)/IRAT + 1
C
               xnrm = ddot(ndim,chovec(1,ir),1,chovec(1,ir),1)
               if (xnrm .ne. 0.0D0) then
                  write(lupri,*) 'Bug!!!'
                  write(lupri,*) 'Vector,isymab,norm :',ir,isymab,xnrm
               end if
               CALL CHO_READ1(CHOVEC(1,IR),NDIM,WORK(KOFVEC),LENG,
     &                        VECNRM)
C
               KOFVEC = KOFVEC + LENG
C
            END DO
C
         IF (ILAST .LT. IEND) GOTO 200
         IF (IR .NE. NR) CALL QUIT('CHO_READ: Error in reading')
C
C        IOPT=3: don't square!
C
         IF (IOPT .EQ. 3) GOTO 999
C
C        Finally, square up each vector
C        tbp - 11/8 2001: check memory!
C
         IF (LWORK .LT. N2BST(ISYMAB)) THEN
            WRITE(LUPRI,*)
     &      'Insufficient memory for squaring vectors in cho_read'
            WRITE(LUPRI,*) 'Symmetry of vectors     : ',ISYMAB
            WRITE(LUPRI,*) 'Memory needed (words)   : ',N2BST(ISYMAB)
            WRITE(LUPRI,*) 'Memory available (words): ',LWORK
            CALL QUIT(' Insufficient memory in cho_read ')
         ENDIF
C
         DO I = 1,NR
            CALL CCSD_SYMSQ(CHOVEC(1,I),ISYMAB,WORK)
            CALL DCOPY(N2BST(ISYMAB),WORK,1,CHOVEC(1,I),1)
         END DO
C
      END IF
C
  999 RETURN
      END
C
      SUBROUTINE CHO_READ1(CHOVEC,NDIM,SCR,LSCR,VECNRM)
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION CHOVEC(NDIM),SCR(LSCR)
C
C
      NUMEL = INT(SCR(1))
C
      KCHISE = 2
      KCHSAV = KCHISE + (NUMEL - 1)/IRAT + 1
C
      CALL CHO_READ2(CHOVEC,SCR(KCHISE),SCR(KCHSAV),NUMEL,NDIM,VECNRM)
C
      RETURN
      END
C
      SUBROUTINE CHO_READ2(CHOVEC,ISELEC,CHOSEL,NUMEL,NDIM,VECNRM)
C
C     Written by Alfredo Sanchez   13 de julio del 2001
C
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <chodbg2.h>
C
      DIMENSION CHOVEC(NDIM),CHOSEL(NUMEL),ISELEC(NUMEL)
      SAVE NCALL
      DATA NCALL /0/
C
C
      call dzero(chovec,ndim)
      DO I = 1,NUMEL
         J = ISELEC(I)
         CHOVEC(J) = CHOSEL(I)
      END DO
C
      VECNR2 = DDOT(NDIM,CHOVEC,1,CHOVEC,1)
      VECNRM = DDOT(NUMEL,CHOSEL,1,CHOSEL,1)
      VECDIF = ABS(VECNRM - VECNR2)
C
      RETURN
      END
C  /* Deck cc_rdaon */
      SUBROUTINE CC_RDAON(XINT,IDELTA,IGAM,NUMPR,ISYMGD,
     *                    IRECNR,WORK,LWORK)
C
C     Written by Henrik Koch 25-Sep-1993
C
C     Purpose: Read destribution of AO integrals.
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION XINT(*),WORK(LWORK)
      DIMENSION IRECNR(*)
      DIMENSION IDELTA(*),IGAM(*)
C
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccpack.h>
C
C----------------------------
C     Construct index arrays.
C----------------------------
C
      KADR2 = 1 
      KEND1 = KADR2 + (NBAST*NBAST + 1)/IRAT + 1
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in CC_RDAO')
      END IF
C
      CALL CC_INXN(WORK(KADR2),ISYMGD)
C
C--------------------
C     Construct XINT.
C--------------------
C
      CALL DZERO(XINT,NNBST(ISYMGD)*NUMPR)
C
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     Buffer allocation
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      LRINT  = IRAT*LBUF
C
      KIBUF = KEND1
      KRBUF = KIBUF + (NIBUF*LBUF-1)/IRAT + 1
      KGAM  = KRBUF + LBUF
      KEND2 = KGAM  + (2*NUMPR-1)/IRAT + 1 
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufiicient work space in CC_RDAO')
      ENDIF
C
      CALL CC_RDA1N(XINT,WORK(KIBUF),WORK(KRBUF),IDELTA,IGAM,
     *             NUMPR,ISYMGD,WORK(KADR2),IBIT1,WORK(KGAM),IRECNR)
C
      RETURN
      END
C  /* Deck cc_rda1n */
      SUBROUTINE CC_RDA1N(XINT,IBUF,RBUF,IDELTA,IGAM,NUMPR,ISYMGD,
     *                  KADR2,IBIT1,KGAM,IRECNR)
C
C     Written by Henrik Koch 25-Sep-1993
C
#include <implicit.h>
#include <ibtpar.h>
      LOGICAL   OPT1, OPT2
      DIMENSION XINT(*)
      DIMENSION KADR2(NBAST,NBAST)
      DIMENSION IBUF(LBUF,NIBUF),RBUF(LBUF)
      DIMENSION IDELTA(NUMPR), IGAM(NUMPR)
      DIMENSION KGAM(NUMPR,2)
      DIMENSION IRECNR(*)
      CHARACTER*8 FAODER
      LOGICAL OLDDX
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <mxcent.h>
#include <nuclei.h>
#include <chrnos.h>
#include <eribuf.h>
#include <inftap.h>
#include <ibtfun.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      OPT1 = .FALSE.
      OPT2 = .FALSE.
C
      IF (NEWDIS) THEN
C
         NEWDIS = .FALSE.
C
         IF (LUINTR .LE. 0) THEN
            CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',
     &        'UNFORMATTED',IDUMMY,.FALSE.)
         END IF
         REWIND (LUINTR)
         DO 50 I = 1,NBUFX(0)
            READ(LUINTR) IRECNR(I)
   50    CONTINUE
C
      ENDIF
C
      IF (LUAORC(0) .LE. 0) THEN
            LBFINP = LBUF
C
#if defined (VAR_INT64)
            LRECL = LBFINP + LBFINP + 1
#elif defined (SYS_NEC)
            LRECL = LBFINP + LBFINP/IRAT + 1
#else
            IF (NBASIS .LE. 255) THEN
               LRECL = 2*LBFINP + LBFINP + 1
            ELSE
               LRECL = 2*LBFINP + 2*LBFINP + 1
            END IF
#endif
            FAODER = 'AO2DIS'//CHRNOS(0)//CHRNOS(0)
            CALL GPOPEN(LUAORC(0),FAODER,'UNKNOWN','DIRECT',
     &           'UNFORMATTED',LRECL,OLDDX)
      END IF
C
      IF (NIBUF .EQ. 1) THEN
C
         DO 100 J = 1,NBUFX(0)
C
            ICOUNT = 0
            DO K = 1,NUMPR
               IF (IRECNR(J) .EQ. IDELTA(K)) THEN
                  ICOUNT = ICOUNT + 1
                  IDEL   = IRECNR(J)
                  KGAM(ICOUNT,1) = IGAM(K)
                  KGAM(ICOUNT,2) = K
               ENDIF
            ENDDO
            IF (ICOUNT .NE. 0) OPT1 = .TRUE.
C
            IF (OPT1) THEN
               NREC   = J
               READ(LUAORC(0),ERR=2000,REC=NREC) RBUF,IBUF,NINT
               DO 110 I = 1,NINT
                  IR = IBTAND(IBTSHR(IBUF(I,1),2*NBITS),IBIT1)
                  DO L = 1,ICOUNT
                     IF (IR .EQ. KGAM(L,1)) THEN
                        OPT2 = .TRUE.
                        LSV  = KGAM(L,2)
                     ENDIF
                  ENDDO 
                  IF (OPT2) THEN
                     IP = IBTAND(       IBUF(I,1)         ,IBIT1)
                     IQ = IBTAND(IBTSHR(IBUF(I,1),  NBITS),IBIT1)
                     IADR = NNBST(ISYMGD)*(LSV-1) + KADR2(IP,IQ) + 1
                     XINT(IADR) = RBUF(I)
                     OPT2 = .FALSE.
                  ENDIF
C 
  110          CONTINUE
C
            ENDIF
C
            OPT1 = .FALSE.
C
  100    CONTINUE
C
      ELSE
C
         DO 120 J = 1,NBUFX(0)
C
            ICOUNT = 0
            DO K = 1,NUMPR
               IF (IRECNR(J) .EQ. IDELTA(K)) THEN
                  ICOUNT = ICOUNT + 1
                  IDEL   = IRECNR(J)
                  KGAM(ICOUNT,1) = IGAM(K)
                  KGAM(ICOUNT,2) = K
               ENDIF
            ENDDO
            IF (ICOUNT .NE. 0) OPT1 = .TRUE.
C
            IF (OPT1) THEN
               NREC   = J
               READ(LUAORC(0),ERR=2000,REC=NREC) RBUF,IBUF,NINT
               DO 130 I = 1,NINT
                  IR = IBTAND(       IBUF(I,1)       ,IBIT1)
                  DO L = 1,ICOUNT
                     IF (IR .EQ. KGAM(L,1)) THEN
                        OPT2 = .TRUE.
                        LSV  = KGAM(L,2)
                     ENDIF
                  ENDDO 
                  IF (OPT2) THEN
                     IP = IBTAND(       IBUF(I,2)       ,IBIT1)
                     IQ = IBTAND(IBTSHR(IBUF(I,2),NBITS),IBIT1)
                     IADR = NNBST(ISYMGD)*(LSV-1) + KADR2(IP,IQ) + 1
                     XINT(IADR) = RBUF(I)
                     OPT2 = .FALSE.
                  ENDIF
  130          CONTINUE
            ENDIF
C
            OPT1 = .FALSE.
C
  120    CONTINUE
C
      ENDIF
C
      RETURN
 2000 CALL QUIT('Error reading AOTWODIS')
      END
C  /* Deck ccrd_inxn */
      SUBROUTINE CC_INXN(KADR2,ISYMAB)
C
C     asm 22-aug-1994
C
C     Purpose: Construct index arrays for CC_RDAO
C
#include <implicit.h>
C
      DIMENSION KADR2(NBAST,NBAST)
C
#include <ccorb.h>
#include <ccsdsym.h>
C
         ICOUN2 = 0
         DO 210 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            IF (ISYMB .GT. ISYMA) THEN

               DO 220 B = 1,NBAS(ISYMB)
                  NB = IBAS(ISYMB) + B
C
                  DO 230 A = 1,NBAS(ISYMA)
                     NA = IBAS(ISYMA) + A
C
                     KADR2(NA,NB) = ICOUN2
                     KADR2(NB,NA) = ICOUN2
C
                     ICOUN2 = ICOUN2 + 1
C
  230             CONTINUE
  220          CONTINUE
C
            ELSE IF (ISYMA .EQ. ISYMB) THEN
C
               DO 240 B = 1,NBAS(ISYMB)
                  NB = IBAS(ISYMB) + B
C
                  DO 250 A = 1,B
                     NA = IBAS(ISYMA) + A
C
                     KADR2(NA,NB) = ICOUN2
                     KADR2(NB,NA) = ICOUN2
C
                     ICOUN2 = ICOUN2 + 1
C
  250             CONTINUE
  240          CONTINUE
C
            END IF
C
  210    CONTINUE
C
      RETURN
      END
C
C/* Deck cc_putcho */
      SUBROUTINE CC_PUTCHO(CHOVEC,IDNTC1,NCHORD,NWRITE,ISYMAB,
     &                     INDRED,WORK,LWORK)
C
C     Written by Alfredo Sanchez   julio 2002
C
C     Write NW vectors at the end of corresponding file.
C
C     The compressed set is without zeros but with labels
C     The reduced set has zeros deleted according to indred
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <choio.h>
C
      DIMENSION CHOVEC(*),WORK(LWORK)
      DIMENSION INDRED(*)
C
      PARAMETER (MXCHRD = 200)
      DIMENSION IDNTC1(MXCHRD)
C
      CHARACTER*10 NAME
C
C
C     Initial stuff
C
      NDIM   = NNBST(ISYMAB)
C
      IFIRST = NUMCHO(ISYMAB) + 1
      ILAST  = IFIRST + NWRITE - 1
C
      NUMCH0 = NUMCHO(ISYMAB) + NWRITE
      IF (NUMCH0 .GT. MXCHVC) THEN
         WRITE(LUPRI,'(A,I3)') 'Too many vectors for symmetry',ISYMAB
         WRITE(LUPRI,'(2(A,I7,/))') 'Required :', NUMCH0,
     &                          'Maximum  :', MXCHVC
         CALL QUIT('Too many Cholesky vectors')
      ELSE
         NUMCHO(ISYMAB) = NUMCH0
      END IF
C
C     Get adress and starting file
C
      JADR = IADRTO(ISYMAB)
C
      DO J = 0,NUMFIL(ISYMAB)-1
         IF (IADRLU(J,ISYMAB) .GT. JADR) GOTO 50
         IFILE = J
      END DO
   50 CONTINUE
       
      LUCHO = IDNTLU(IFILE,ISYMAB)
      IADR  = JADR - IADRLU(IFILE,ISYMAB) + 1
C
      CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'OLD')
C
C     Get length and check if a new file is needed
C     The compressed packing needs special treatment
C
      IF (COMP) GOTO 100
C
      IF (REDUCE) THEN
         MDIM = NREDUC(ISYMAB)
      ELSE 
         MDIM = NNBST(ISYMAB)
      END IF
C
      NFIL1 = NWRITE
      NFIL2 = 0
C
      LENGTH = NWRITE * MDIM
      
      JFINIS = IADR + LENGTH - 1
      IF (JFINIS .GT. MAXLEN) THEN
         IROOM  = MAXLEN - IADR + 1
         NFIL1  = IROOM / MDIM
         NFIL2  = NWRITE
c        write(6,*) 'File splitted in writing'
c        write(6,*) 'isymab,ifirst,ilast :',isymab,ifirst,ilast
c        write(6,*) 'nwrite,nfil1,nfil2 :',nwrite,nfil1,nfil2
c        write(6,*) 'iadr,jadr :',iadr,jadr
      END IF
C
C     Write non compressed vectors.
C
      IF (.NOT. REDUCE) THEN
C
         LENTOT = NFIL1 * MDIM
C
         DO I = 1,NFIL1
            II = I + IFIRST - 1
            IDNTCH(II,ISYMAB) = IDNTC1(I)
            LENCHO(II,ISYMAB) = MDIM
            IADRTO(ISYMAB)    = IADRTO(ISYMAB) + MDIM
         END DO
C
         CALL PUTWA2(LUCHO,NAME,CHOVEC,IADR,LENTOT)
C
      ELSE
C
         LENTOT = NFIL1 * MDIM
C
         KWRIT = 1
         KEND  = KWRIT + LENTOT
         LWRK  = LWORK - KEND + 1
         IF (LWRK .LT. 0) CALL QUIT('Not enough space in CC_PUTCHO')
C
         KOFF1 = KWRIT
         KOFF2 = 1
         LENG1 = NREDUC(ISYMAB)
         LENG2 = NNBST(ISYMAB)
         DO J = 1,NFIL1
C
c           write(6,*) 'reducing vector : ',j
            JJ = J + IFIRST - 1
C
            IDNTCH(JJ,ISYMAB) = IDNTC1(J)
            LENCHO(JJ,ISYMAB) = MDIM
            IADRTO(ISYMAB)    = IADRTO(ISYMAB) + MDIM
C
            DO I = 1,NREDUC(ISYMAB)
               II = INDRED(LREDU*(ISYMAB-1)+I)
               IOFF1 = KOFF1 + I  - 1
               IOFF2 = KOFF2 + II - 1
               WORK(IOFF1) = CHOVEC(IOFF2)
            END DO
C
            KOFF1 = KOFF1 + LENG1
            KOFF2 = KOFF2 + LENG2
C
         END DO
C
c        write(6,*) 'lucho,name :',lucho,name
c        write(6,*) 'first,last :',work(kwrit),work(kwrit+lentot-1)
c        write(6,*) 'chovec :',chovec(1),chovec(ioff2)
         CALL PUTWA2(LUCHO,NAME,WORK(KWRIT),IADR,LENTOT)
C
      END IF
      IF (NFIL2 .EQ. 0) GOTO 999
C
C     New file is needed
C
      IADR  = 1
      IFILE = IFILE + 1
      CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'NEW')
      IADRLU(IFILE,ISYMAB) = IADRTO(ISYMAB)
C
      IF (.NOT. REDUCE) THEN
C
         LENTOT = NFIL2 * MDIM
C
         DO I = NFIL1+1,NFIL2
            II = I + IFIRST - 1
            IDNTCH(II,ISYMAB) = IDNTC1(I)
            LENCHO(II,ISYMAB) = MDIM
            IADRTO(ISYMAB)    = IADRTO(ISYMAB) + MDIM
         END DO
C
         KOFF = NFIL1*MDIM + 1
         CALL PUTWA2(LUCHO,NAME,CHOVEC(KOFF),IADR,LENTOT)
C
      ELSE
C
         LENTOT = NFIL2 * MDIM
C
         KWRIT = 1
         KEND  = KWRIT + LENTOT
         LWRK  = LWORK - KEND + 1
         IF (LWRK .LT. 0) CALL QUIT('Not enough space in CC_PUTCHO')
C
         KOFF1 = KWRIT
         KOFF2 = 1 + NFIL1 * NNBST(ISYMAB)
         LENG1 = NREDUC(ISYMAB)
         LENG2 = NNBST(ISYMAB)
         DO J = NFIL1+1,NFIL2
C
c           write(6,*) 'reducing vector : ',j
            JJ = J + IFIRST - 1
C
            IDNTCH(JJ,ISYMAB) = IDNTC1(J)
            LENCHO(JJ,ISYMAB) = MDIM
            IADRTO(ISYMAB)    = IADRTO(ISYMAB) + MDIM
C
            DO I = 1,NREDUC(ISYMAB)
               II = INDRED(LREDU*(ISYMAB-1)+I)
               IOFF1 = KOFF1 + I  - 1
               IOFF2 = KOFF2 + II - 1
               WORK(IOFF1) = CHOVEC(IOFF2)
            END DO
C
            KOFF1 = KOFF1 + LENG1
            KOFF2 = KOFF2 + LENG2
C
         END DO
C
c        write(6,*) 'lucho,name :',lucho,name
c        write(6,*) 'iadr,jadr :',iadr,jadr
c        write(6,*) 'first,last :',work(kwrit),work(kwrit+lentot-1)
c        write(6,*) 'chovec :',chovec(1 + NFIL1*MDIM),chovec(ioff2)
         CALL PUTWA2(LUCHO,NAME,WORK(KWRIT),IADR,LENTOT)
C
      END IF
      GOTO 999
C
  100 CONTINUE
C
C     Pack vectors
C
      KOFF = 1
      KEND = 1
      LWRK = LWORK
C
      LENTOT = 0
C
      I = IFIRST - 1
  200 CONTINUE
C
         I = I + 1
C
         THRSEL = 0.0D0
         CALL CHO_PACK(THRSEL,NDIM,CHOVEC(KOFF),WORK(KEND),
     &                 LWRK,LENGTH,NUMEL)
C
         II = I - IFIRST + 1
         IDNTCH(I,ISYMAB) = IDNTC1(II)
         LENCHO(I,ISYMAB) = LENGTH
         IADRTO(ISYMAB)   = IADRTO(ISYMAB) + LENGTH
C
         KOFF = KOFF + NDIM
         KEND = KEND + LENGTH
         LWRK = LWRK - LENGTH
C
         LENTOT = LENTOT + LENGTH
         JFINIS = IADR   + LENTOT - 1
C
         IF (JFINIS .GT. MAXLEN) THEN
            LENTOT = LENTOT - LENGTH
            CALL PUTWA2(LUCHO,NAME,WORK,IADR,LENTOT)
C
            IADR  = 1
            IFILE = IFILE + 1
            CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'NEW')
C
            CALL DCOPY(LENGTH,WORK(KEND),1,WORK(1),1)
C
            KEND = 1 + LENGTH
            LWRK = LWORK - LENGTH
         END IF
C
      IF (I .LT. ILAST) GOTO 200
C
      CALL PUTWA2(LUCHO,NAME,WORK,IADR,LENTOT)
C
  999 CONTINUE
C
C     Restart info
C
      REWIND(LURST)
      WRITE(LURST) NSYM,THRCOM,LRDTOT,
     &             (NNBST(ISYM),ISYM=1,NSYM),
     &             (NUMCHO(ISYM),ISYM=1,NSYM),
     &             (IADRTO(ISYM),ISYM=1,NSYM),
     &             (NUMFIL(ISYM),ISYM=1,NSYM)
C
      WRITE(LURST) (INDRED(I),I=1,LRDTOT)
C
      DO ISYM = 1,NSYM
         WRITE(LURST) (LENCHO(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         WRITE(LURST) (IDNTCH(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         WRITE(LURST) (IADRLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
      DO ISYM = 1,NSYM
         WRITE(LURST) (IDNTLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
      CALL FLSHFO(LURST)
C
      RETURN
      END
C
C
      SUBROUTINE CHO_PACK(THRSEL,NDIM,CHOVEC,WORK,LWORK,LENGTH,NUM)
C
C     Written by Alfredo Sanchez   julio 2002
C
C     Select the relevant Cholesky elements
C
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION CHOVEC(NDIM),WORK(LWORK)
C
C
C     Memory allocation.
C
      KCHLEN = 1
      KCHISE = KCHLEN + 1
      KCHSAV = KCHISE + (NDIM - 1)/IRAT + 1
      KEND   = KCHSAV + NDIM
      LWRK   = LWORK  - KEND + 1
      IF (LWRK .LT. 0) THEN
         WRITE(LUPRI,*) 'Not enough space in CHO_PACK.1'
         WRITE(LUPRI,'(A,2I10)')'Available,needed:',LWRK,KEND
         CALL QUIT('Not enough space in CHO_PACK')
      END IF
C
      CALL CHO_PACK1(THRSEL,NDIM,CHOVEC,WORK(KCHISE),WORK(KCHSAV),
     &               NUM,.FALSE.)
C
      IF (LWRK .LT. NUM) THEN
         WRITE(LUPRI,*) 'Not enough space in CHO_PACK.2'
         WRITE(LUPRI,'(A,2I10)')'Available,needed:',LWRK,NUM
         CALL QUIT('Not enough space in CHO_PACK')
      END IF
C
      WORK(KCHLEN) = DFLOAT(NUM)
C
      KOFF = KCHISE + (NUM - 1)/IRAT + 1
      CALL DCOPY(NUM,WORK(KCHSAV),1,WORK(KEND),1)
      CALL DCOPY(NUM,WORK(KEND),1,WORK(KOFF),1)
         
      LENGTH = (NUM - 1)/IRAT + NUM + 2
C
      RETURN
      END
C
      SUBROUTINE CHO_PACK1(THRSEL,NDIM,CHOVEC,ISELEC,CHOSEL,
     &                     NUMEL,INFO)
C
C     Written by Alfredo Sanchez   julio 2002
C
C     Select the relevant Cholesky elements
C
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
      DIMENSION CHOVEC(NDIM),CHOSEL(NDIM),ISELEC(NDIM)
C
      SAVE NTOT,NSAV
      DATA NTOT,NSAV /2*0/
C
      LOGICAL INFO
C
C
      IF (INFO) THEN
          WRITE(LUPRI,'(//,A,2I15,//)') 
     &         'Number of elements (total and saved):',NTOT,NSAV
          CALL FLSHFO(LUPRI)
          RETURN
      END IF
C
      ICOUNT = 0
C
      DO I = 1,NDIM
         IF (ABS(CHOVEC(I)) .GE. THRSEL) THEN
            ICOUNT = ICOUNT + 1
            ISELEC(ICOUNT) = I
            CHOSEL(ICOUNT) = CHOVEC(I)
         END IF
      END DO
C
      NUMEL = ICOUNT
C
      NTOT = NTOT + NDIM
      NSAV = NSAV + NUMEL
C
      RETURN
      END
C
C/* Deck cc_getcho */
      SUBROUTINE CC_GETCHO(CHOVEC,NREAD,ISYMAB,IFIRST,INDRED,WORK,LWORK)
C
C     Written by Alfredo Sanchez   julio 2002
C
C     Read NREAD vectors.
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <choio.h>
C
      DIMENSION CHOVEC(*),WORK(LWORK)
      DIMENSION INDRED(*)
C
      CHARACTER*10 NAME
      INTEGER MDIM,LONG1
C
C
C     Initial stuff
C
      ILAST  = IFIRST + NREAD - 1
      IF (ILAST .GT. NUMCHO(ISYMAB)) THEN
          ILAST = NUMCHO(ISYMAB)
          NREAD = ILAST - IFIRST + 1
      END IF
      IF (ILAST .EQ. 0) RETURN
C
C     Determine file name and address
C
      IF (COMP) THEN
         JADR = 1
         DO I = 1,IFIRST-1
            MDIM = LENCHO(I,ISYMAB)
            JADR = JADR + MDIM
         END DO
      ELSE
         IF (REDUCE) THEN
            MDIM   = NREDUC(ISYMAB)
         ELSE
            MDIM   = NNBST(ISYMAB)
         END IF
         LONG1 = IFIRST - 1
         JADR  = MDIM * LONG1 + 1
      END IF
C
      DO J = 0,NUMFIL(ISYMAB)-1
         IF (IADRLU(J,ISYMAB) .GT. JADR) GOTO 50
         IFILE = J
      END DO
   50 CONTINUE
C
      IADR  = JADR - IADRLU(IFILE,ISYMAB) + 1
C
      CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'OLD')
C
C     Read non compressed vectors.
C
      IF (COMP) GOTO 100
C
      NFIL1 = NREAD
      NFIL2 = 0
C
      LENTOT = MDIM * NREAD
      JFINIS = IADR + LENTOT - 1
      IF (JFINIS .GT. MAXLEN) THEN
         IROOM  = MAXLEN - IADR + 1
         NFIL1  = IROOM / MDIM
         NFIL2  = NREAD - NFIL1
c        write(6,*) 'File splitted in reading'
c        write(6,*) 'isymab,ifirst,ilast :',isymab,ifirst,ilast
c        write(6,*) 'nfil1,nfil2 :',nfil1,nfil2
c        write(6,*) 'jadr,iadr : ',jadr,iadr
      END IF
C
C     Read complete non compressed set
C
      IF (.NOT. REDUCE) THEN
          IF (NFIL2 .EQ. 0) THEN
C
             CALL GETWA2(LUCHO,NAME,CHOVEC,IADR,LENTOT)
C
          ELSE
C
             LENTOT = NFIL1 * MDIM            
             CALL GETWA2(LUCHO,NAME,CHOVEC,IADR,LENTOT)
C
             LENTOT = NFIL2 * MDIM
             IFILE  = IFILE + 1
             CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'OLD')
C
             KOFF = MDIM * NFIL1 + 1
             IADR = 1
             CALL GETWA2(LUCHO,NAME,CHOVEC(KOFF),IADR,LENTOT)
C
          END IF
          RETURN
      END IF
C
C
C     Read reduced set.
C
      IF (REDUCE) THEN
C
          LENTOT = NREDUC(ISYMAB) * NREAD
          LENGTH = NNBST(ISYMAB) * NREAD
C
          KREAD = 1
          KEND  = KREAD + LENTOT
          LWRK  = LWORK - KEND + 1
          IF (LWRK .LT. 0) CALL QUIT('Not enough space in CC_GETCHO')
C
          CALL DZERO(CHOVEC,LENGTH)
          IF (NFIL2 .EQ. 0) THEN
             CALL GETWA2(LUCHO,NAME,WORK(KREAD),IADR,LENTOT)
          ELSE
C
             LENTOT = NFIL1 * MDIM
c            write(6,*) 'Reading vectors from :', NAME,NFIL1
c            write(6,*) 'jadr,iadr : ',jadr,iadr
c            write(6,*) 'kread :',kread,kread+lentot-1
             CALL GETWA2(LUCHO,NAME,WORK(KREAD),IADR,LENTOT)
c            write(6,*) 'First and last :',work(kread),
c    &                   work(kread+lentot-1)
C
             LENTOT = NFIL2 * MDIM
             IFILE  = IFILE + 1
             CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'OLD')
C
             KOFF = KREAD + NFIL1*MDIM
             IADR = 1
c            write(6,*) 'Reading vectors from :', NAME,NFIL2
c            write(6,*) 'jadr,iadr : ',jadr,iadr
c            write(6,*) 'kread :',koff,kread+lentot-1
             CALL GETWA2(LUCHO,NAME,WORK(KOFF),IADR,LENTOT)
c            write(6,*) 'First and last :',work(koff),
c    &                   work(kread+lentot-1)
         END IF
C
C        Get complete set
C
         KOFF1 = KREAD
         KOFF2 = 1
         LENG1 = NREDUC(ISYMAB)
         LENG2 = NNBST(ISYMAB)
         DO J = IFIRST,ILAST
            DO I = 1,NREDUC(ISYMAB)
               II = INDRED(LREDU*(ISYMAB-1)+I)
               IOFF1 = KOFF1 + I  - 1
               IOFF2 = KOFF2 + II - 1
               CHOVEC(IOFF2) = WORK(IOFF1)
            END DO
            KOFF1 = KOFF1 + LENG1
            KOFF2 = KOFF2 + LENG2
         END DO
C
         RETURN
C
      END IF
C
  100 CONTINUE
C
C     Read compressed vectors.
C
      LENGT2 = 0
      LENTOT = 0
      DO I = IFIRST,ILAST
         LENTOT = LENTOT + LENCHO(I,ISYMAB)
         JFINIS = IADR   + LENTOT
         IF (JFINIS .GT. MAXLEN) 
     &      LENGT2 = LENGT2 + LENCHO(I,ISYMAB)
      END DO
C
      KREAD = 1
      KEND1 = KREAD + LENTOT
      LWRK1 = LWORK - KEND1 + 1
C
      IF (LWRK1 .LT. 0) CALL QUIT('Not enough space in cc_getcho')
C
      CALL DZERO(CHOVEC,LENGTH)
      IF (LENGT2 .EQ. 0) THEN
         CALL GETWA2(LUCHO,NAME,WORK(KREAD),IADR,LENTOT)
      ELSE
C
         LENGT1 = LENTOT - LENGT2
         CALL GETWA2(LUCHO,NAME,WORK(KREAD),IADR,LENGT1)
C
         IFILE = IFILE + 1
         CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'OLD')
C
         IADR = 1
         KOFF = KREAD + LENGT1
         CALL GETWA2(LUCHO,NAME,WORK(KOFF),IADR,LENGT2)
C
      END IF
C
C     Unpack vectors
C
      KOFF1 = KREAD
      KOFF2 = 1
      LENG2 = NNBST(ISYMAB)
      DO I = IFIRST,ILAST
C
         LENG1 = LENCHO(I,ISYMAB)
         CALL CHO_UNPACK(WORK(KOFF1),LENG1,CHOVEC(KOFF2),LENG2)
C
         KOFF1 = KOFF1 + LENG1
         KOFF2 = KOFF2 + LENG2
C
      END DO
C
      RETURN
      END
C
      SUBROUTINE CHO_UNPACK(SCR,LSCR,CHOVEC,NDIM)
C
C     Alfredo Sanchez.  Julio 2002
C
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION CHOVEC(NDIM),SCR(LSCR)
C
C
      NUMEL = INT(SCR(1))
C
      KCHISE = 2
      KCHSAV = KCHISE + (NUMEL - 1)/IRAT + 1
C
      CALL CHO_UNPACK1(CHOVEC,SCR(KCHISE),SCR(KCHSAV),
     &                 NUMEL,NDIM)
C
      RETURN
      END
C
      SUBROUTINE CHO_UNPACK1(CHOVEC,ISELEC,CHOSEL,NUMEL,NDIM)
C
C     Written by Alfredo Sanchez   julio 2002
C
C
#include <implicit.h>
C
      DIMENSION CHOVEC(NDIM),CHOSEL(NUMEL),ISELEC(NUMEL)
C
C
      DO I = 1,NUMEL
         J = ISELEC(I)
         CHOVEC(J) = CHOSEL(I)
      END DO
C
      RETURN
      END
C
C  /* Deck cc_anadi */
      SUBROUTINE CC_ANADI(DIAG,ALL)
C
C     Alfredo Sanchez   Julio 2002
C
C     Histogram on diagonal.
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DIAG(*)
      LOGICAL ALL
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <priunit.h>
C
C
      NUMH = 0
      NUM0 = 0
      NUM2 = 0
      NUM4 = 0
      NUM6 = 0
      NUM8 = 0
      NUM10 = 0
      NUM12 = 0
      NUM15 = 0
      NUM20 = 0
      NUMZ  = 0
      NUMN  = 0
      NUM11 = 0
C
      DO I = 1,NDIAG
         IF (DIAG(I) .LT. 0.0D0) NUMN = NUMN + 1
c
         IF (DIAG(I) .GT. 1.0d2) THEN
              NUMH = NUMH + 1
         ELSE IF (DIAG(I) .GT. 1.0d0) THEN
               NUM0 = NUM0 + 1
         ELSE IF (DIAG(I) .GT. 1.0d-2) THEN
               NUM2 = NUM2 + 1
         ELSE IF (DIAG(I) .GT. 1.0d-4) THEN
               NUM4 = NUM4 + 1
         ELSE IF (DIAG(I) .GT. 1.0d-6) THEN
               NUM6 = NUM6 + 1
         ELSE IF (DIAG(I) .GT. 1.0d-8) THEN
               NUM8 = NUM8 + 1
         ELSE IF (DIAG(I) .GT. 1.0d-10) THEN
               NUM10 = NUM10 + 1
         ELSE IF (DIAG(I) .GT. 1.0d-12) THEN
               NUM12 = NUM12 + 1
         ELSE IF (DIAG(I) .GT. 1.0d-15) THEN
               NUM15 = NUM15 + 1
         ELSE IF (DIAG(I) .eq. ZERO) THEN
               NUMZ = NUMZ + 1
         ELSE
               NUM20 = NUM20 + 1
         END IF
      END DO
      NUMS = NUM10 + NUM12 + NUM15 + NUM20
c
      NUMT = NUMH
      WRITE(LUPRI,'(A,2I10)') 'Larger than 1.0D+2          : ',NUMH,NUMT
      NUMT = NUMT + NUM0
      WRITE(LUPRI,'(A,2I10)') 'Between 1.0D+2 and 1.0D+0   : ',NUM0,NUMT
      NUMT = NUMT + NUM2
      WRITE(LUPRI,'(A,2I10)') 'Between 1.0D+0 and 1.0D-2   : ',NUM2,NUMT
      NUMT = NUMT + NUM4
      WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-2 and 1.0D-4   : ',NUM4,NUMT
      NUMT = NUMT + NUM6
      WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-4 and 1.0D-6   : ',NUM6,NUMT
      NUMT = NUMT + NUM8
      WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-6 and 1.0D-8   : ',NUM8,NUMT
      IF (ALL) THEN
         NUMT = NUMT + NUM10
         WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-8 and 1.0D-10  : ',
     &                       NUM10,NUMT
         NUMT = NUMT + NUM12
         WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-10 and 1.0D-12 : ',
     &                       NUM12,NUMT
         NUMT = NUMT + NUM15
         WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-12 and 1.0D-15 : ',
     &                       NUM15,NUMT
         NUMT = NUMT + NUM20
         WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-15 and 0.0D+00 : ',
     &                       NUM20,NUMT
      ELSE
         NUMT = NUMT + NUMS
         WRITE(LUPRI,'(A,2I10)') 'Between 1.0D-8 and 0.0D+00  : ',
     &                       NUMS,NUMT
      END IF
      NUMT = NUMT + NUMZ 
      WRITE(LUPRI,'(A,2I10)') 'Exactly 0.0d+00             : ',NUMZ,NUMT
      WRITE(LUPRI,'(A,I10,///)') 'Negative diagonals          : ',NUMN
C
      RETURN
      END
C
C
C  /* Deck cho_chkdia */
      SUBROUTINE CHO_CHKDIA(DIAG,DIAOLD,OLDIAG,XC,ICAB,ICHOV,
     &                     XMAX1,ISYMAB,INDRED)
C
C     Alfredo Sanchez   Julio 2002
C
C     Check diagonal after updating.
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DIAG(*),DIAOLD(*)
      DIMENSION INDRED(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <choio.h>
#include <priunit.h>
C
C
      ICHOV = ICHOV + 1
C
      XM = ZERO
      INEG = 0
      DO I = 1,NNBST(ISYMAB)
C
         KOFF1  = IDIAG(ISYMAB) + I
         DIATMP = DIAG(KOFF1)
C
         IF (XM .LT. ABS(DIATMP)) XM = ABS(DIATMP)
C
         IF (DIATMP .LT. -1.0D-13) THEN
C
            DIAG(KOFF1) = ZERO
            INEG = INEG + 1
C
            IF (DIATMP .LT. -1.0D-10) THEN
                WRITE(LUPRI,'(A,I12,2D18.8)') 'Negdiag.,init :',
     &                           I,DIATMP,DIAOLD(KOFF1)
            END IF
C
         END IF
C
      END DO
C
      ICONV = 0
      DO I = 1,NNBST(ISYMAB)
         KOFF1 = IDIAG(ISYMAB) + I
         DIAKO = SQRT(ABS(DIAG(KOFF1)) * XM)*THSUDI
         IF (DIAKO .LT. THRCOM) THEN
               IF (SCDIAG) DIAG(KOFF1) = ZERO
               ICONV = ICONV + 1
         END IF
      END DO
C
      IF (IPRCHO .GT. 2) THEN
          WRITE(LUPRI,'(I3,I8,I10,1P,2D14.5,I8,I6,D14.5)')
     &         ISYMAB,ICHOV,ICAB,XC,OLDIAG,ICONV,INEG,XMAX1
          CALL FLSHFO(LUPRI)
      END IF
C
C     Write secure restart file
C
      IF (MOD(ICHOV,500) .EQ. 0) THEN
C
         REWIND(LUSEC)
         WRITE(LUSEC) NSYM,THRCOM,LRDTOT,
     &             (NNBST(ISYM),ISYM=1,NSYM),
     &             (NUMCHO(ISYM),ISYM=1,NSYM),
     &             (IADRTO(ISYM),ISYM=1,NSYM),
     &             (NUMFIL(ISYM),ISYM=1,NSYM)
C
         WRITE(LUSEC) (INDRED(I),I=1,LRDTOT)
C
         DO ISYM = 1,NSYM
            WRITE(LUSEC) (LENCHO(I,ISYM),I=1,NUMCHO(ISYM))
         END DO
         DO ISYM = 1,NSYM
            WRITE(LUSEC) (IDNTCH(I,ISYM),I=1,NUMCHO(ISYM))
         END DO
         DO ISYM = 1,NSYM
            WRITE(LUSEC) (IADRLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
         END DO
         DO ISYM = 1,NSYM
            WRITE(LUSEC) (IDNTLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
         END DO
         CALL FLSHFO(LUSEC)
C
      END IF
C
      RETURN
      END
C
C
C  /* Deck cho_end */
      SUBROUTINE CHO_END(DIAG,NTOINT,XSEL0,YSEL0)
C
C     Alfredo Sanchez   Julio 2002
C
C     Final analysis after overall convergence.
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DIAG(*)
      DIMENSION XSEL0(8),YSEL0(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <chotim.h>
#include <blocks.h>
#include <priunit.h>
C
C
      WRITE(LUPRI,'(/)')
C
      CALL AROUND('Summary of Cholesky decomposition')
      YC = -1.0d0
      DO I = 1,NDIAG
         IF (ABS(DIAG(I)) .GT. YC) THEN
            YC = ABS(DIAG(I))
            IY = I
         ENDIF
         IF (DIAG(I) .LT. -THRCOM) WRITE(LUPRI,'(A,I8,D15.5)')
     &      'Negative diagonal:',I,DIAG(I)
      ENDDO
      WRITE(LUPRI,'(/,2(9X,A,/))') 'Number of vectors',
     &                          '-----------------'
      WRITE(LUPRI,'(6X,A,9X,A,10X,A)') 'ISAB','Maximum','Actual'
      DO ISYM = 1,NSYM
         WRITE(LUPRI,'(7X,I2,2I16)') ISYM,NNBST(ISYM),NUMCHO(ISYM)
      END DO
C
c     WRITE(LUPRI,'(//,2(9X,A,/))') 'Number of elements',
c    &                              '------------------'
c     WRITE(LUPRI,'(6X,A,8X,A,3(9X,A))') 'ISAB','Reduced','Real 0',
c    &                                   'Zeroed','Non zero'
c     WRITE(LUPRI,'(6X,A,7X,A,1X,2(9X,A))') 'ISAB','Reduced','Real 0',
c    &                                   'Non zero'
c     DO ISYM = 1,NSYM
c        XELEM = DFLOAT(NUMCHO(ISYM))*DFLOAT(NREDUC(ISYM))
c        YELEM = XELEM - YSEL0(ISYM) - XSEL0(ISYM)
c        WRITE(LUPRI,'(7X,I2,4F16.0)') ISYM,XELEM,YSEL0(ISYM),
c    &                                 XSEL0(ISYM),YELEM
c        WRITE(LUPRI,'(7X,I2,4F16.0)') ISYM,XELEM,YSEL0(ISYM),
c    &                                 YELEM
c     END DO
C
      WRITE(LUPRI,'(//,A,D15.5)') 
     *     'Largest absolute element in diagonal:',DIAG(IY)
C
      MAXINT = MAXSHL * (MAXSHL+1) / 2
      WRITE(LUPRI,'(A,2I10,//)')
     &         'Number of distributions (total and actual):',
     &          MAXINT,NTOINT
C
      TIMTOT = SECOND() - TIMTOT
      TIMDEC = TIMTOT - (TIMINI + TIMDIA + TIMINT)
C
      WRITE(LUPRI,'(//,2(9X,A,/))') 'Time statistics',
     &                              '---------------'
      WRITE(LUPRI,'(/,6X,A,F10.2)')  'Initialization time : ', TIMINI
      WRITE(LUPRI,'(6X,A,F10.2)')    'Diagonal time       : ', TIMDIA
      WRITE(LUPRI,'(6X,A,F10.2)')    'Integral time       : ', TIMINT
      WRITE(LUPRI,'(6X,A,F10.2)')    'Decomposition time  : ', TIMDEC
      WRITE(LUPRI,'(/,6X,A,F10.2)')  'Total time          : ', TIMTOT
C
      WRITE(LUPRI,'(//,A,//)') 'Decomposition completed' 
C
      RETURN
      END
C  /* Deck cc_diaglr0 */
      SUBROUTINE CC_DIAGLR0(XC,ISYMXC,ISHELA,ISHELB)
C
C     Written by Henrik Koch 15-juni-2000
C
C     Find largest diagonal element. 
C
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxorb.h>
#include <maxash.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <blocks.h>
#include <ccdeco.h>
C
      XC = -1.0d0
      DO ISHLA = 1,MAXSHL
         DO ISHLB = ISHLA,MAXSHL
            XXX = DIASCR(ISHLA,ISHLB)
            IF (XXX. GT. XC) THEN
               XC = XXX
               ISHELA = ISHLA
               ISHELB = ISHLB
            END IF
         END DO
      END DO
C
      ISYMXC = ISYSCR(ISHELA,ISHELB)
C
      DIASCR(ISHELA,ISHELB) = ZERO
      DIASCR(ISHELB,ISHELA) = ZERO
C
      RETURN
      END
C
C  /* Deck cho_save */
      SUBROUTINE CHO_SAVE(INDRED)
C
C     Write restar files at end of decomposition
C
C     Alfredo Sanchez de Meras    Julo 2002
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDRED(*)
#include <maxorb.h>
#include <maxash.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <blocks.h>
#include <ccdeco.h>
#include <choio.h>
#include <priunit.h>
C
C
      REWIND(LUSEC)
C
      WRITE(LUSEC) NSYM,THRCOM,LRDTOT,
     &             (NNBST(ISYM),ISYM=1,NSYM),
     &             (NUMCHO(ISYM),ISYM=1,NSYM),
     &             (IADRTO(ISYM),ISYM=1,NSYM),
     &             (NUMFIL(ISYM),ISYM=1,NSYM)
C
      WRITE(LUSEC) (INDRED(I), I=1,LRDTOT)
C
      DO ISYM = 1,NSYM
         WRITE(LUSEC) (LENCHO(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         WRITE(LUSEC) (IDNTCH(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         WRITE(LUSEC) (IADRLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
      DO ISYM = 1,NSYM
         WRITE(LUSEC) (IDNTLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
      CALL FLSHFO(LUSEC)
C
C
      REWIND(LURST)
C
      WRITE(LURST) NSYM,THRCOM,LRDTOT,
     &             (NNBST(ISYM),ISYM=1,NSYM),
     &             (NUMCHO(ISYM),ISYM=1,NSYM),
     &             (IADRTO(ISYM),ISYM=1,NSYM),
     &             (NUMFIL(ISYM),ISYM=1,NSYM)
C
      WRITE(LURST) (INDRED(I), I=1,LRDTOT)
C
      DO ISYM = 1,NSYM
         WRITE(LURST) (LENCHO(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         WRITE(LURST) (IDNTCH(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         WRITE(LURST) (IADRLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
      DO ISYM = 1,NSYM
         WRITE(LURST) (IDNTLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
      CALL FLSHFO(LURST)
C
      RETURN
      END
C  /* Deck cc_rdaox */
      SUBROUTINE CC_RDAOX(XINT,NUMDIA,NUMDIB,JNDEXA,JNDEXB,
     *                    ISHLA,ISHLB,WORK,LWORK,IRECNR)
C
C     Written by Henrik Koch 25-Sep-1993
C
C     Purpose: Read destribution of AO integrals.
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <iratdef.h>
C
      DIMENSION XINT(*),WORK(LWORK)
      DIMENSION JNDEXA(*), JNDEXB(*)
      DIMENSION IRECNR(*)
C
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccpack.h>
C
C--------------------
C     Construct XINT.
C--------------------
C
      CALL DZERO(XINT,NUMDIA*NUMDIB)
C
      IBIT1 = 2**NBITS     - 1
      IBIT2 = 2**(2*NBITS) - 1
C
C     Buffer allocation
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      LRINT  = IRAT*LBUF
C
      KIBUF = 1
      KRBUF = KIBUF + (NIBUF*LBUF-1)/IRAT + 1
      KEND2 = KRBUF + LBUF
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufiicient work space in CC_RDAO')
      ENDIF
C
      CALL CC_RDAOX1(XINT,WORK(KIBUF),WORK(KRBUF),NUMDIA,NUMDIB,
     *             JNDEXA,JNDEXB,ISHLA,ISHLB,IBIT1,IRECNR)
C
      RETURN
      END
C  /* Deck cc_rdaox1 */
      SUBROUTINE CC_RDAOX1(XINT,IBUF,RBUF,NUMDIA,NUMDIB,
     *                  JNDEXA,JNDEXB,ISHLA,ISHLB,IBIT1,IRECNR)
C
C     Written by Henrik Koch 25-Sep-1993
C
#include <implicit.h>
#include <ibtpar.h>
#include <priunit.h>
      DIMENSION XINT(*)
      DIMENSION JNDEXA(*), JNDEXB(*)
      DIMENSION IBUF(LBUF,NIBUF),RBUF(LBUF)
      DIMENSION IRECNR(*)
      CHARACTER*8 FAODER
      LOGICAL OLDDX
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
#include <ccinftap.h>
#include <ccorb.h>
#include <mxcent.h>
#include <eribuf.h>
#include <nuclei.h>
#include <inftap.h>
#include <chrnos.h>
#include <ibtfun.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      IF (NEWDIS) THEN
C
         NEWDIS = .FALSE.
C
         IF (LUINTR .LE. 0) THEN
            CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',
     &        'UNFORMATTED',IDUMMY,.FALSE.)
         END IF
         REWIND (LUINTR)
C
         DO 50 I = 1,NBUFX(0)
            READ(LUINTR) IRECNR(I)
   50    CONTINUE
C
      ENDIF
C
      IF (LUAORC(0) .LE. 0) THEN
            LBFINP = LBUF
C    
#if defined (VAR_INT64)
            LRECL = LBFINP + LBFINP + 1
#elif defined (SYS_NEC)
            LRECL = LBFINP + LBFINP/IRAT + 1
#else
            IF (NBASIS .LE. 255) THEN
               LRECL = 2*LBFINP + LBFINP + 1
            ELSE
               LRECL = 2*LBFINP + 2*LBFINP + 1
            END IF
#endif
            FAODER = 'AO2DIS'//CHRNOS(0)//CHRNOS(0)
            CALL GPOPEN(LUAORC(0),FAODER,'UNKNOWN','DIRECT',
     &           'UNFORMATTED',LRECL,OLDDX)
      END IF
C
      ICOUNT = 0
C
      IF (NIBUF .EQ. 1) THEN
C
         DO 100 J = 1,NBUFX(0)
C
            IF (IRECNR(J) .EQ. 1) THEN
               ICOUNT = ICOUNT + 1
               NREC   = J
               READ(LUAORC(0),ERR=2000,REC=NREC) RBUF,IBUF,NINT
               DO 110 I = 1,NINT
                  IA = IBTAND(IBTSHR(IBUF(I,1),3*NBITS),IBIT1)
                  IB = IBTAND(IBTSHR(IBUF(I,1),2*NBITS),IBIT1)
                  IC = IBTAND(IBTSHR(IBUF(I,1),  NBITS),IBIT1)
                  ID = IBTAND(       IBUF(I,1)         ,IBIT1)
                  IADR = NUMDIA*(JNDEXB(IB) - 1) + JNDEXA(IA) 
                  XINT(IADR) = RBUF(I)
                  IF (ISHLA .EQ. ISHLB) THEN
                     IADR = NUMDIA*(JNDEXB(IA) - 1) + JNDEXA(IB) 
                     XINT(IADR) = RBUF(I)
                  ENDIF
  110          CONTINUE
            ENDIF
C
  100    CONTINUE
C
      ELSE
C
         DO 120 J = 1,NBUFX(0)
C
            IF (IRECNR(J) .EQ. 1) THEN
               ICOUNT = ICOUNT + 1
               NREC   = J
               READ(LUAORC(0),ERR=2000,REC=NREC) RBUF,IBUF,NINT
               DO 130 I = 1,NINT
                  IA = IBTAND(IBTSHR(IBUF(I,1),NBITS),IBIT1)
                  IB = IBTAND(       IBUF(I,1)       ,IBIT1)
                  IC = IBTAND(IBTSHR(IBUF(I,2),NBITS),IBIT1)
                  ID = IBTAND(       IBUF(I,2)       ,IBIT1)
                  IADR = NUMDIA*(JNDEXB(IB) - 1) + JNDEXA(IA) 
                  XINT(IADR) = RBUF(I)
                  IF (ISHLA .EQ. ISHLB) THEN
                     IADR = NUMDIA*(JNDEXB(IA) - 1) + JNDEXA(IB) 
                     XINT(IADR) = RBUF(I)
                  ENDIF
  130          CONTINUE
            ENDIF
C
  120    CONTINUE
C
      ENDIF
C
      CALL GPCLOSE(LUAORC(0),'KEEP')
C
      RETURN
 2000 CALL QUIT('Error reading AOTWODIS')
      END
C
C  /* Deck cc_gabold */
      SUBROUTINE CC_GABOLD(DIAG,WORK,LWORK)
C
C     Written by Henrik Koch 15-juni-2000
C
C     Calculate diagonal elements. 
C
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
      LOGICAL   SETUP
      DIMENSION DIAG(*), WORK(LWORK)
      DIMENSION INDEXA(MXCORB), INDEXB(MXCORB)
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <distcl.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccfro.h>
#include <ccfop.h>
#include <ccsections.h>
#include <symsq.h>
#include <ccdeco.h>
#include <choles.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      WRITE(LUPRI,'(A,I5)') 'Number of shells: ', MAXSHL
      write(lupri,*) 'Using old diagonal subroutine'
      DO ISHLA = 1,MAXSHL
         DO ISHLB = 1,ISHLA
C
            SETUP = .FALSE.
            CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     *                  NUMDIA,NUMDIB,IPRINT,SETUP,2)
C
            DO JA = 1,NUMDIA
               DO JB = 1,NUMDIB
C
                  IA = INDEXA(JA) 
                  IB = INDEXB(JB)
C
                  ISYMA  = ISAO(IA)
                  ISYMB  = ISAO(IB)
                  ISYMAB = MULD2H(ISYMA,ISYMB)
C
                  KCHO1  = 1 
                  KIREC  = KCHO1 + NNBST(ISYMAB)
                  KEND2  = KIREC + (NBUFX(0) - 1)/IRAT + 1
                  LWORK2 = LWORK - KEND2
                  IF (LWORK2 .LT. 0)
     &            CALL QUIT('Not enough space in CC_GAB')
C
                  NEWDIS = .TRUE.
                  CALL CC_RDAO(WORK(KCHO1),IA,IB,WORK(KEND2),LWORK2,
     *                      WORK(KIREC))
C
                  A = IA - IBAS(ISYMA)
                  B = IB - IBAS(ISYMB) 
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB = INDEX(A,B)
                  ELSE IF (ISYMB .GT.ISYMA) THEN
                     NAB = NBAS(ISYMA)*(B - 1) + A
                  ELSE
                     NAB = NBAS(ISYMB)*(A - 1) + B
                  END IF
C
                  NAB = IAODPK(ISYMA,ISYMB) + NAB
C
                  KOFF1 = IDIAG(ISYMAB) + NAB
                  KOFF2 = KCHO1 + NAB - 1
C
                  DIAG(KOFF1) = WORK(KOFF2)
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
      RETURN
      END
C
C
C
Casm  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C     For some extrange and unknown reason, the calculation of the
C     exchange contribution to Fock matrix using Cholesky decomposed
C     integrals fails in the call to DGEMM just before the end of the 
C     789 loop in subroutine SCF_CHOFCK1 (line 2118 of sirfck.F). Such
C     bug has appeared when linking against ATLAS library using ifort
C     on a Linux box and when using gfortran and system libraries on
C     a Mac. However, in both cases the bug dissapeared if the blas
C     routines distributed in the release were employed. Moreover, no 
C     bug appeared using ifort and MKL library. 
C
C     Simply calling the DGEMM routine from pdpack/gp_dblas3.F at that
C     particular point of SCF_CHOFCK1, things seem to be OK and all 
C     Cholesky tests pass.
C
Casm  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C
      SUBROUTINE OTHER_DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, 
     $                   B, LDB, BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  OTHER_DGEMM  is an exact copy of DGEMM code from pdpack/gp_dblas3.F
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
