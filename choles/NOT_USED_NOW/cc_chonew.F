C  /* Deck cc_cholesky */
      SUBROUTINE CC_CHOLESKY(DIAG,DIAOLD,ICO2RE,IRE2CO,WORK,LWORK)
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
      DIMENSION ICO2RE(*),IRE2CO(*)
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
c
c      write(lupri,*) '  ire2co     ico2re'
c      do i = 1,ndiag
c         write(lupri,'(2I10)') ire2co(i),ico2re(i)
c      end do
C
C     Initializations.
C     ----------------
C
      TIMTOT = SECOND()
      TIMINI = SECOND()
      TIMINT = ZERO
      TSELTO = ZERO
      TRDTO  = ZERO
      TPRTO  = ZERO
      TDGTO  = ZERO
C
      ICHOV  = 0
      NTOINT = 0
C
C     Open files for storing the Cholesky vectors.
C     --------------------------------------------
C
      IF (.NOT. RSTCHO)
     &    CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'INI')
C
      NAME  = ROOT(1:7)//'SEL'
      CALL WOPEN2U(28,NAME,64,0)
C
C     Setup index information.
C     ------------------------
C
      ISHLB  = 1
      IPRINT = 0
      SETUP  = .TRUE.
      DO ISHLA = 1,MAXSHL
C
         T1 = SECOND()
         CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     *               NUMDIA,NUMDIB,IPRINT,SETUP,1)
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
         CALL CHO_RESTART(DIAG,WORK,LWORK)
         CALL CHO_REDUCE(DIAG,ICO2RE,IRE2CO)
C
      END IF
C
C     Allocation stuff.
C     -----------------
C
      LWRKX = LWORK / 3
C
      MXREDU = 0
      MXNNBS = 0
      DO I = 1,NSYM
         IF (NREDUC(I) .GT. MXREDU) MXREDU = NREDUC(I) + 1
         IF (NNBST(I)  .GT. MXNNBS) MXNNBS = NNBST(I)
      END DO
C
      LENAB = MXREDU + (MXREDU-1)/IRAT + 2
C
      LREAD = (NBAST*NBAST + 1)/IRAT + 1 
     &      + (NIBUF*LBUF-1)/IRAT + 1 + LBUF
     &      + (2*MAXQUA-1)/IRAT + 1
C
      MINCHO = LENAB + MXREDU + 1
      MINEED = LREAD + MXNNBS
C
      NCHORD = MIN((LWRKX-MINEED)/MINCHO,MXCHRD)
      IF (NCHORD .LT. 1) THEN 
         WRITE(LUPRI,*) 'Memory problem. NCHORD :',NCHORD
         STOP 'Not enough space in CC_CHOLESKY'
      END IF
      NCHRDM = NCHORD
C
      LCHO2 = LENAB*NCHORD + MXREDU+ 1
      LCHO1 = LCHO2 + MXNNBS*NCHORD+ 1
      LNEED = LCHO1 + LREAD
C
      MAXSEL = MIN((LWORK-LNEED)/MXREDU,MAXQUA)
      IF (MAXSEL .LT. 1) THEN 
         WRITE(LUPRI,*) 'Memory problem. MAXSEL :',MAXSEL
         STOP 'Not enough space in CC_CHOLESKY'
      END IF
C
      WRITE(LUPRI,'(/,A,I4)') 'Maximum number of vectors read        :',
     &                    NCHRDM
      WRITE(LUPRI,'(A,I4)') 'Maximum number of diagonals qualified :',
     &                  MAXSEL
C
      WRITE(LUPRI,'(/,A,//)') 'End of Cholesky initializations'
      CALL FLSHFO(LUPRI)
C
      TIMINI = SECOND() - TIMINI
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
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      CALL CC_DIAGLR0(XC,ISYMXC,ISHLA,ISHLB)
C
      IF (XC .EQ. ZERO) GOTO 120
C
      DO I = 1,NSYM
         ITMP(I) = NQUAL(I) + 1
      END DO
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
          WRITE(LUPRI,*) 'Something is wrong with this couple'
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
      WRITE(LUPRI,'(2A,2I6)') 'Entering new distribution ',
     &            'corresponding to shells :',ISHLA,ISHLB
      CALL FLSHFO(LUPRI)
c     WRITE(LUPRI,'(3(A,I8,X))') 'CD distributions. Total :',IALTOT,
c    &                       'Calculated :',IALDON,'Skipped :',IALSKP
C
      T1 = SECOND()
      CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     *            NUMDIA,NUMDIB,IPRINT,SETUP,1)
      T1 = SECOND() - T1
      TIMINT = TIMINT + T1
C
      NEWDIS = .TRUE.
C
C     Select diagonals in this distribution.
C     --------------------------------------
C
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
C
C     Read integrals from disk. 
C     -------------------------
C
      TIMSEL = SECOND()
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
         IF (LWRK1 .LT. 0) STOP 'Not enough space in CC_CHOLESKY'
C
         NEWDIS = .TRUE.
         CALL CC_RDAON(WORK(KINT1),IQUVEA(I,ISYMAB),IQUVEB(I,ISYMAB),
     *                 NUMAB,ISYMAB,WORK(KIREC),WORK(KEND1),LWRK1)
C
         LUNIT  = 28
         NAME   = ROOT(1:7)//'SEL'
         LENGTH = NNBST(ISYMAB) * NUMAB
C
         IF (LENGTH .NE. 0) THEN
c           write(6,*) 'mxnnbs :',mxnnbs
c           write(6,*) 'maxsel :',maxsel
c           write(6,*) 'isymab :',isymab
c           write(6,*) 'nnbst  :',nnbst(isymab) 
c           write(6,*) 'itmp   :',itmp(isymab) 
            IADR = MXNNBS * MAXSEL * (ISYMAB-1)
     &            + NNBST(ISYMAB) * (ITMP(ISYMAB)-1) + 1
c           write(6,*) 'iadr   :',iadr  
            CALL PUTWA2U(LUNIT,NAME,WORK(KINT1),IADR,LENGTH)
c           write(lupri,*) 'Norm of integrals',isymab,
c    &                      ddot(length,work(kint1),1,work(kint1),1)
         END IF
C
  350 CONTINUE
      TIMSEL = SECOND() - TIMSEL
      WRITE(LUPRI,'(A,F12.2)') 'Time selecting :',TIMSEL
      TSELTO = TSELTO + TIMSEL
      
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
      WRITE(LUPRI,'(////A,I4)')
     &     'Number of calculated distributions :',NUMINT
      WRITE(LUPRI,'(A,8I4)')
     &     'Number of qualified diagonals      :',(NQUAL(I),I=1,NSYM)
      WRITE(LUPRI,'(A,D17.8//)') 
     &     'Maximum diagonal among them        :',XMAX2
      WRITE(LUPRI,'(A,15X,A,12X,2(4X,A),9X,A)')
     &  'ISAB    Iter','Treated','Conv','Neg','Max'
C
C     Loop over symmetries in file to decompose.
C     ------------------------------------------
C
      DO 200 ISYMAB = 1,NSYM
C
         IF (NREDUC(ISYMAB) .EQ. 0) GOTO 200
         IF (NQUAL(ISYMAB) .EQ. 0)  GOTO 200
C
C        Allocations.
C        ------------
C
C        LENAB : Maximum length of compressed vector
C        LENABT: Plus length of actual reduced vector
C        NVEFIL: Maximum number of vectors in one file
C 
         LENAB  = MXNNBS + (MXNNBS-1)/IRAT + 2
         LENABT = LENAB + NREDUC(ISYMAB) + 1
C
         NVEFIL = MAXLEN / LENAB
C
C        LWRKX = LWORK - space for factor array
         LWRKX = LWORK - NCHRDM*NQUAL(ISYMAB)
C
C        LWRKX = LWRKX - space for vectors of selected diagonals
         LWRKX = LWRKX - NREDUC(ISYMAB)*NQUAL(ISYMAB)
C
         NCHORD = LWRKX / LENABT
C
         NCHORD = MIN(NCHORD,NUMCHO(ISYMAB))
         IF (NCHORD .EQ. 0) NCHORD = MIN(10, LWRKX/LENABT)
         NCHORD = MIN(NCHORD,NVEFIL)
C
         NCHORD = MIN(NCHORD,NCHRDM)
C
c        write(lupri,*) 'Vectors read:',nchord
C
         KINT1 = 1
         KFUL1 = KINT1 + NREDUC(ISYMAB)*NQUAL(ISYMAB) + 1
         KEND1 = KFUL1 + NNBST(ISYMAB)*NQUAL(ISYMAB)
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Insufficient work space in CC_CHOLESKY.1'
            WRITE(LUPRI,*) 'Needed and available memory',KEND1,LWORK
            STOP 'Insufficient work space in CC_CHOLESKY.1'
         ENDIF
C
C        Read integrals from disk.
C        -------------------------
C
         LUNIT  = 28
         NAME   = ROOT(1:7)//'SEL'
         LENGTH = NNBST(ISYMAB) * NQUAL(ISYMAB)
C
C        Full vectors read in KFUL1 and passed reduced to KINT1
C
         IF (LENGTH .NE. 0) THEN
            IADR = MXNNBS * MAXSEL * (ISYMAB-1) + 1
            CALL GETWA2U(LUNIT,NAME,WORK(KFUL1),IADR,LENGTH)
c           if (isymab .eq. 4) then
c           write(lupri,*) 'nqual :',nqual(4)
c           write(lupri,*) 'Norm of read vector :',
c    &         ddot(length,work(kful1),1,work(kful1),1)
c           call flshfo(lupri)
c           end if
            CALL CHO_REDVEC(WORK(KFUL1),WORK(KINT1),NQUAL(ISYMAB),
     &                      IRE2CO,ISYMAB)
c           if (isymab .eq. 4) then
c           write(lupri,*) 'Norm of reduced vector :',
c    &   ddot(NREDUC(ISYMAB)*NQUAL(ISYMAB),work(kint1),1,work(kint1),1)
c           call flshfo(lupri)
c           end if
         END IF
C
C        Reallocate.
C        KCHO1 is reused later on to keep new vectors before writing to disk
C
         KCHO1 = KFUL1 
         KEND2 = KCHO1 + NREDUC(ISYMAB)*NCHORD + 1 
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Insufficient work space in CC_CHOLESKY.2'
            WRITE(LUPRI,*) 'Needed and available memory',KEND2,LWORK
            STOP 'Insufficient work space in CC_CHOLESKY.2'
         ENDIF
C
C        Batch over previous cholesky vectors.
C        -------------------------------------
C
         NBATV = (NUMCHO(ISYMAB)-1)/NCHORD + 1
c        if (isymab .eq. 4) write(lupri,*) 'nbatv',nbatv
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
c
c           write(lupri,*) 'Before calling cc_getcho, ICO2RE: '
c           do kk = 1,nnbst(isymab)
c              kkk = kk + idiag(isymab)
c              write(lupri,'(2i6)') ico2re(kkk)
c           end do
c           write(lupri,*) 'Calling cc_getcho from/to',ibatv1,ibatv2
            CALL CC_GETCHO(WORK(KCHO1),NUMVEC,ISYMAB,
     &                     IBATV1,WORK(KEND2),LWRK2,ICO2RE)
            TIMX  = SECOND() - TIMX
            TIMRD = TIMRD + TIMX
            TRDTO = TRDTO + TIMX 
c
c           if (isymab .eq. 4) then
c           do n = 1,numvec
c           write(lupri,*) 'Vector read',n
c           do i = 1,nreduc(isymab)
c              koff = kcho1 + nreduc(isymab)*(n-1) + i - 1
c              write(lupri,'(i4,d20.10,i4)') i,work(koff),
c    &                        ire2co(idiag(isymab)+i)
c           end do
c              write(lupri,*)
c           end do
c           end if
C
C           Prepare factor array.
C           ---------------------
C
            KCHO2 = KEND2
            KEND3 = KCHO2 + NUMVEC*NQUAL(ISYMAB)
            LWRK3 = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
               WRITE(LUPRI,*) 'Insufficient work space in CC_CHOLESKY.3'
               WRITE(LUPRI,*) 'Needed and available memory',KEND3,LWORK
               STOP 'Insufficient work space in CC_CHOLESKY.3'
            END IF
C
            TIMX = SECOND()
            DO J = 1,NQUAL(ISYMAB)
               DO JJ = 1,NUMVEC
                  KOFF1 = KCHO1 
     *                  + NREDUC(ISYMAB)*(JJ-1) 
     *                  + ICO2RE(IQUVAB(J,ISYMAB)) - 1
                  KOFF2 = KCHO2 + NUMVEC*(J-1) + JJ - 1
                  WORK(KOFF2) = WORK(KOFF1)
c                 if (isymab .eq.4)
c    &            write(lupri,*) 'j,jj,work(koff2)',j,jj,work(koff2),
c    &                           isymab
               ENDDO
            ENDDO
            TIMX  = SECOND() - TIMX
            TIMPR = TIMPR + TIMX
            TPRTO = TPRTO + TIMX 
C
C           Subtraction.
C           ------------
C
            NTOTAB = MAX(NREDUC(ISYMAB),1)
C
            TIMX = SECOND()
            CALL DGEMM('N','N',NREDUC(ISYMAB),NQUAL(ISYMAB),NUMVEC,
     *                 -ONE,WORK(KCHO1),NTOTAB,WORK(KCHO2),
     *                 NUMVEC,ONE,WORK(KINT1),NTOTAB) 
            TIMX  = SECOND() - TIMX
            TIMDG = TIMDG + TIMX
            TDGTO = TDGTO + TIMX 
C
         END DO
C
         WRITE(LUPRI,'(2A,3F12.2)') 'Time in reading, preparing,',
     &                          ' multiplying :',TIMRD,TIMPR,TIMDG
c          write(lupri,*) 'Norm of cho1 :',
c    &         ddot(ntotab,work(kcho1),1,work(kcho1),1)
c           write(lupri,*) 'Norm of cho2 :',
c    &         ddot(nqual(isymab),work(kcho2),1,work(kcho2),1)
  444    CONTINUE
c         if (isymab .eq. 4) then
c        len8 = nreduc(isymab)*nqual(isymab)
c        write(lupri,*) 'Norm of read vector after substraction :',
c    &         ddot(len8,work(kint1),1,work(kint1),1)
c        end if
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
c           if (isymab .eq. 4) write(lupri,*) 'xd',xd
            KOFF = KINT1 + NREDUC(ISYMAB)*(IC-1)
            CALL DSCAL(NREDUC(ISYMAB),XD,WORK(KOFF),1)
            kofsav = koff
C
C           Screen the vector
C
ctst        CHOMAX = ZERO
ctst        DO I = 1,NREDUC(ISYMAB)
ctst           KOFF = KINT1 + NREDUC(ISYMAB)*(IC-1) + I - 1
ctst           IF (ABS(WORK(KOFF)) .GT. CHOMAX) 
ctst &             CHOMAX = ABS(WORK(KOFF))
ctst        END DO
ctst        XNUM   = DFLOAT(20 * NBAST)
ctst        THRSEL = THRCOM/CHOMAX/XNUM
            thrsel = zero
C
            XSEL1 = ZERO
            YSEL1 = ZERO
            DO I = 1,NREDUC(ISYMAB)
               II = IDIAG(ISYMAB) + IRE2CO(IDIAG(ISYMAB) + I)
               KOFF = KINT1 + NREDUC(ISYMAB)*(IC-1) + I - 1
               IF (DIAG(II) .EQ. ZERO) WORK(KOFF) = ZERO
ctst               IF (WORK(KOFF) .EQ. ZERO) THEN
ctst                  YSEL1 = YSEL1 + 1.0D0
ctst               ELSE IF (ABS(WORK(KOFF)) .LT. THRSEL) THEN
ctst                  WORK(KOFF) = ZERO
ctst                  XSEL1 = XSEL1 + 1.0D0
ctst           END IF
            END DO
C
ctst            XZERO = NNBST(ISYMAB) - NREDUC(ISYMAB)
ctst            YSEL1 = YSEL1 - XZERO
C
ctst            XSEL0(ISYMAB) = XSEL0(ISYMAB) + XSEL1
ctst            YSEL0(ISYMAB) = YSEL0(ISYMAB) + YSEL1
C
            IF (IPRINT .GT. 3) WRITE(LUPRI,'(A,3X,2(A,F8.0,2X),A,I8)')
     &                              'Elementents in vector',
     &                              'Real 0 :',YSEL1,
     &                              'Zeroed :',XSEL1,
     &                              'Length :',LENAB
c           if (isymab .eq. 4)
c    &      write(lupri,*) 'Norm of vector before updating :',
c    &         ddot(nreduc(isymab),work(kofsav),1,work(kofsav),1),
c    &         '    ',isymab
c           call flshfo(lupri)
C
C           Update diagonal
C
            XM = ZERO
            DO I = 1,NREDUC(ISYMAB)
               II = IDIAG(ISYMAB) + IRE2CO(IDIAG(ISYMAB) + I)
               KOFF = KINT1 + NREDUC(ISYMAB)*(IC-1) + I - 1
               DIAG(II) = DIAG(II) - WORK(KOFF)*WORK(KOFF)
               IF (DIAG(II) .GT. XM) XM = DIAG(II)
            ENDDO
C
            DIAMIN(ISYMAB) = XM * SPAN
            IF (DIAMIN(ISYMAB) .LT. THRCOM) DIAMIN(ISYMAB) = THRCOM
C
C           Analyze updated diagonal.
C           -------------------------
C
            OLDIAG     = DIAG(ICAB)
            DIAG(ICAB) = ZERO
            CALL CHO_CHKDIA(DIAG,DIAOLD,OLDIAG,XC,ICAB,ICHOV,
     &                      XM,ISYMAB)
c           write(lupri,*) 'Diag(10)',diag(10)
C
C           Substract this vector from the other qualified ones.
C           ----------------------------------------------------
C
            DO I = 1,NQUAL(ISYMAB)
               II = IQUVAB(I,ISYMAB)
c              write(lupri,*) 'ii,diag,i',ii,diag(ii),i
               IF (DIAG(II) .NE. ZERO) THEN
C
                  KOFF1 = KINT1 + NREDUC(ISYMAB) * (IC-1)
                  KOFF2 = KINT1 + NREDUC(ISYMAB) * (I-1) 
C
                  JJ = ICO2RE(II)
                  FACTOR = -WORK(KOFF1+JJ-1)
C
c               write(lupri,*) 'Before daxpy'
c               write(lupri,*) 'ic,factor1,work(koff2),koff2 :',ic,
c    &                              factor,work(koff2),koff2
c               do kk = 1,nreduc(isymab)
c                  write(lupri,'(i4,2d20.10)') kk,work(koff1+kk-1),
c    &                                         work(koff2+kk-1)
c               end do
C
                  CALL DAXPY(NREDUC(ISYMAB),FACTOR,
     *                       WORK(KOFF1),1,WORK(KOFF2),1)
C
c               if (isymab .eq. 4) then
c               write(lupri,*) 'After daxpy'
c                 write(lupri,*) 'Norm of vector ii',ii,
c    &            ddot(nreduc(isymab),work(koff2),1,work(koff2),1)
c               do kk = 1,nreduc(isymab)
c                  write(lupri,'(i4,d20.10)') kk,work(koff2+kk-1)
c               end do
c               end if
C
               ENDIF
            ENDDO
C
C           Write cholesky vectors to disk.
C           -------------------------------
C
            IDUMP         = IDUMP + 1
            IDNTC1(IDUMP) = ICAB - IDIAG(ISYMAB)
C
            KOFF1 = KINT1 + NREDUC(ISYMAB)*(IC-1)
            KOFF2 = KCHO1 + NREDUC(ISYMAB)*(IDUMP-1)
            CALL DCOPY(NREDUC(ISYMAB),WORK(KOFF1),1,WORK(KOFF2),1)
c
c           if (isymab .eq. 4) then
c           if (numcho(isymab) .ge. 12) then
c           write(lupri,*) 'Vector to write of diagonal',idntc1(idump)
c           do i = 1,nreduc(isymab)
c              write(lupri,'(i4,d20.10,i6)') i,work(koff1+i-1),
c    &                                       ire2co(idiag(isymab)+i)
c           end do
c           call flshfo(lupri)
c           end if
c           end if
C
  333       CONTINUE
C
            IF ((IDUMP.EQ.NCHORD).OR.(ICHO.EQ.NQUAL(ISYMAB))
     *                           .OR.LAST) THEN
c
c              write(lupri,*) 'Dumping vectors isym/idump',isymab,idump
               CALL CC_PUTCHO(WORK(KCHO1),IDNTC1,NCHORD,IDUMP,
     &                        ISYMAB,IRE2CO,WORK(KEND2),LWRK2)
c
c              if (isymab.eq.4) write(lupri,*)'Vectors dump',idump
               IF (LAST) GOTO 200
               IDUMP = 0
            ENDIF
C
  600    CONTINUE 
  200 CONTINUE
C
C     Check if we can decompose more in this distribution.
C     ----------------------------------------------------
C
      DO ISYMAB = 1,NSYM
         NQUAL(ISYMAB) = 0
         ITMP(ISYMAB)  = 1
      END DO
C
      XMAX2 = -1.0D0
      DO ISYMAB = 1,NSYM
         XMAX1 = -1.0d0
         DO I = 1,NREDUC(ISYMAB)
             XXX = DIAG(IDIAG(ISYMAB)+IRE2CO(IDIAG(ISYMAB)+I))
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
         WRITE(LUPRI,*) 'Not enough diagonals qualified at second pass'
         WRITE(LUPRI,*) 'NQUAL : ', (NQUAL(I),I=1,NSYM)
         GOTO 101      
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
         LUNIT  = 28
         NAME   = ROOT(1:7)//'SEL'
         LENGTH = NNBST(ISYMAB) * NUMAB
C
         IF (LENGTH .NE. 0) THEN
            IADR = MXNNBS * MAXSEL * (ISYMAB-1)
     &            + NNBST(ISYMAB) * (ITMP(ISYMAB)-1) + 1
            CALL PUTWA2U(LUNIT,NAME,WORK(KINT1),IADR,LENGTH)
         END IF
C
  351 CONTINUE
C
      GOTO 120
C
  101 CONTINUE
C
C     Statistics after treating distribution(s).
C     ------------------------------------------
C
      WRITE(LUPRI,'(//,A,/)') ' >>> Checking diagonal'
      CALL CC_ANADI(DIAG,.FALSE.,.FALSE.,.TRUE.)
c     do i = 1,ndiag
c        write(lupri,*) i,diag(i)
c     end do
C
C     Get new reduced set and go back to begin of loop
C     ------------------------------------------------
C
      CALL CHO_REDUCE(DIAG,ICO2RE,IRE2CO)
c
c     write(lupri,*) 'diagonal and mapping vectors'
c     do isym=1,nsym
c        do i = 1,nnbst(isym)
c           koff = idiag(isym)+i
c           write(lupri,'(a,i5,f15.10,2i6)') 'diagonal,ico2re,ire2co',
c    &            i,diag(koff),ico2re(koff),ire2co(koff)
c        end do
c        write(lupri,*)
c     end do
C
      GOTO 100
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
         STOP 'Insufficient space for allocation in CC_RDAO'
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
         STOP 'Insufiicient work space in CC_RDAO'
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
#if defined (SYS_CRAY) || defined (SYS_T90)
            LRECL = LBFINP + LBFINP + 1
#endif
#if defined (SYS_AIX) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_LINUX)
            IF (NBASIS .LE. 255) THEN
               LRECL = 2*LBFINP + LBFINP + 1
            ELSE
               LRECL = 2*LBFINP + 2*LBFINP + 1
            END IF
#endif
#if defined (SYS_IRIX) || defined (SYS_DEC)
            IF (NBASIS .LE. 255) THEN
               LRECL = 2*LBFINP + LBFINP + 1
            ELSE
               LRECL = 2*LBFINP + 2*LBFINP + 1
            END IF
#endif
#if defined (SYS_NEC)
            LRECL = LBFINP + LBFINP/IRAT + 1
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
 2000 STOP 'Error reading AOTWODIS'
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
                  IF (LWORK2 .LT. 0) STOP 'Not enough space in CC_GAB'
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
#include <mxcent.h>
#include <aovec.h>
#include <blocks.h>
#include <eribuf.h>
C
      PARAMETER (ONTM10 = 1.0D-10)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CHARACTER*10 NAME
C
      DIMENSION TITMOL(12,2)
      DIMENSION WORK(LWORK)
      DIMENSION MMBST(8)
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <ccdeco2.h>
#include <choio.h>
C
C
      CALL AROUND('Using new Cholesky code')
C
C------------------------------------------------------
C     Get basis information and construct index arrays.
C------------------------------------------------------
C
      CALL GPOPEN(LUONEL,'AOONEINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUONEL)
C
      READ (LUONEL) ((TITMOL(I,J),I=1,12),J=1,2),
     &               NSYM,(NBAS(ISYM),ISYM=1,NSYM),POTNUC
C
      CALL GPCLOSE(LUONEL,'KEEP')
      CALL CHO_INDEX
c     CALL FLSHFO(LUPRI)
C
C----------------
C     Read input.
C----------------
C
      CALL CC_CHOINP
      CALL FLSHFO(LUPRI)
C
      IF (SKPDEC) RETURN
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      LREDU = 0
      DO ISYM = 1,NSYM
         IF (LREDU .LT. NNBST(ISYM)) LREDU = NNBST(ISYM)
      END DO
C
      KINRE = 1
      KDIAG = KINRE + NSYM * ((LREDU-1)/IRAT + 1)
      KOLD  = KDIAG + NDIAG
      KEND1 = KOLD  + NDIAG
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) STOP 'Not enough space in CC_CHODRV.1'
C
      LWRTMP = LWORK - KOLD
C----------------------
C     Get the diagonal.
C----------------------
C
      CALL DZERO(WORK(KDIAG),NDIAG)
C
      LUDIAG = 0
      CALL GPOPEN(LUDIAG,'CHODIAG','UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUM,.FALSE.)
C
C     Initialize integral calculation.
C     --------------------------------
C
      IPRINT = 0 
      T1 = SECOND()
      CALL HERDI1(WORK,LWORK,IPRINT)
      T1 = SECOND() - T1
      TIMINT = TIMINT + T1
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(A,I5)')  'Number of shells:',MAXSHL
      WRITE(LUPRI,*)
C
C     Calculate diagonal elements.
C     ----------------------------
C
      TIMDIA = SECOND()
      REWIND(LUDIAG)
      IF (RSTDIA) THEN
         READ(LUDIAG) (WORK(KDIAG+I-1),I = 1,NDIAG)
      ELSE
         CALL CC_GAB(WORK(KDIAG),WORK(KOLD),LWRTMP)
         WRITE(LUDIAG) (WORK(KDIAG+I-1),I = 1,NDIAG)
         RSTDIA = .TRUE.
         CALL FLSHFO(LUDIAG)
      END IF
c     write(lupri,*) 'Diagonal :',ndiag
c     do i = 1,ndiag
c        write(lupri,*) i,'   ',work(kdiag+i-1)
c     end do
c     CALL DCOPY(NDIAG,WORK(KDIAG),1,WORK(KOLD),1)
C
C     Analyse the diagonal.
C     ---------------------
C
      WRITE(LUPRI,'(//,A,/)') ' >>> Histogram of  initial diagonal'
      CALL CC_ANADI(WORK(KDIAG),.TRUE.,.TRUE.,.FALSE.)
      TIMDIA = SECOND() - TIMDIA
C
C     Screen out the initial diagonal.
C     --------------------------------
C
      IZER = 0
      XM   = ZERO
      DO ISYM = 1,NSYM
C
         XMS = ZERO
         DO I = 1,NNBST(ISYM)
            KOFF1 = KDIAG + IDIAG(ISYM) + I - 1
            IF (XMS .LT. ABS(WORK(KOFF1))) XMS = ABS(WORK(KOFF1))
         END DO
C
         XSCREEN = THRDEF*THRDEF/XMS/THINDI
C
         DO I = 1,NNBST(ISYM)
            KOFF1 = KDIAG + IDIAG(ISYM) + I - 1
            IF (ABS(WORK(KOFF1)) .LT. XSCREEN) THEN
               WORK(KOFF1) = ZERO
               IZER = IZER + 1
            END IF
         END DO
C
         IF (XM .LT. XMS) XM = XMS
C
      END DO
C
C     Be sure that the threshold is not too crazy.
C     --------------------------------------------
C
      IF (XM .LT. THRDEF) THEN
         WRITE(LUPRI,'(//,2(A,D17.10,/),/)') 
     &               'Maximum diagonal        :', XM,
     &               'Decomposition threshold :', THRDEF
         WRITE(LUPRI,'(A,/,A)') 
     &               'The descomposition is initially converged',
     &               'There must be some error. The program will stop'
         STOP 'Threshold too large'
      END IF
C
      THRCOM = XM * 1.0D-15
      IF (THRCOM .LT. THRDEF) THEN
         THRCOM = THRDEF
      ELSE
         WRITE (LUPRI,'(A,/,A)') 'WARNING : ',
     &         'Threshold redefined because of too large diagonals'
      END IF
C
      WRITE(LUPRI,'(/,A,2D20.10,//)') 'Largest diagonal,thrcom :',
     &                                 XM,THRCOM
C
C     Get reference reduced set and reset memory.
C     -------------------------------------------
C
      CALL CHO_REDINI(WORK(KINRE),WORK(KDIAG),WORK(KOLD),
     &                WORK(KEND1),LWRK1)
C
      CALL DCOPY(LENRED,WORK(KOLD),1,WORK(KDIAG),1)
C
      KOLD  = KDIAG + LENRED
      KEND1 = KOLD  + LENRED
      LWRK1 = LWORK - KEND1
C
      CALL DCOPY(LENRED,WORK(KDIAG),1,WORK(KOLD),1)
C
C     Get first reduced set.
C     ----------------------
C
      KIND1 = KEND1
      KIND2 = KIND1 + (LENRED-1)/IRAT + 1
      KEND2 = KIND2 + (LENRED-1)/IRAT + 1
      LWRK2 = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Not enough memory in CC_CHODRV'
         WRITE(LUPRI,'(A,I10)') 'Needed    :', KEND2
         WRITE(LUPRI,'(A,I10)') 'Available :', LWRK2
         STOP 'Not enough space in CC_CHODRV.2'
      END IF
C
      ICOUNT = 0
      DO ISYM = 1,NSYM
         IREDIA(ISYM) = ICOUNT
         NREDUC(ISYM) = NREINI(ISYM)
         ICOUNT = ICOUNT + NREINI(ISYM)
      END DO
C
      CALL CHO_REDUCE(WORK(KDIAG),WORK(KIND1),WORK(KIND2))
C
C     Some info to output.
C
      WRITE(LUPRI,'(/,A,I10,//)') 'Zeroed initial diagonals : ',IZER
C
      WRITE(LUPRI,'(//,2(9X,A,/))') 
     &            'Diagonal information after screening',
     &            '------------------------------------'
      WRITE(LUPRI,'(6X,A,8X,A,3(7X,A))') 'ISAB',' NNBST ',' NREDUC',
     &                                   'IDIAG ','IREDIA'
C
      DO ISYM = 1,NSYM
         WRITE(LUPRI,'(7X,I2,4I16)') ISYM,NNBST(ISYM),NREDUC(ISYM),
     &                               IDIAG(ISYM),IREDIA(ISYM)
      END DO
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(5X,A,2I16)') 'Total',NDIAG,LENRED
      WRITE(LUPRI,*)
      CALL FLSHFO(LUPRI)
C
C----------------------
C     Open files
C----------------------
C
      LURST = 0
      CALL GPOPEN(LURST,'CHOLESKY.RST','UNKNOWN','SEQUENTIAL',
     &            'UNFORMATTED',IDUM,.FALSE.)
C
      LUSEC = 0
      CALL GPOPEN(LUSEC,'CHOLESKY.SEC','UNKNOWN','SEQUENTIAL',
     &            'UNFORMATTED',IDUM,.FALSE.)
C
C-------------------------
C     Restart information.
C-------------------------
C
      IF (RSTCHO) THEN
C
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) '        -------------------------'
         WRITE(LUPRI,*) '          Restarted calculation'
         WRITE(LUPRI,*) '        -------------------------'
         WRITE(LUPRI,*)
C
C        Use ldum instead of comp/reduce for compatibilty with Cholesky
C        files older than last splitting. 
C
         REWIND(LURST)
C_orig   READ(LURST,ERR=666,END=666) MSYM,THRCOM,COMP,REDUCE,LRDTOT,
         READ(LURST,ERR=666,END=666) MSYM,THRCOM,ldum,ldum,LRDTOT,
     &               (MMBST(ISYM),ISYM=1,MSYM),
     &               (NUMCHO(ISYM),ISYM=1,MSYM),
     &               (IADRTO(ISYM),ISYM=1,MSYM),
     &               (NUMFIL(ISYM),ISYM=1,MSYM)
C
         IF (MSYM .NE. NSYM) THEN
            WRITE(LUPRI,*) 'Error in symmetry checking at restart'
            GOTO 666
         END IF
         DO I = 1,NSYM
            IF (MMBST(I) .NE. NNBST(I)) THEN
               WRITE(LUPRI,*) 'Error in basis checking at restart'
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
         READ(LUSEC) MSYM,THRCOM,ldum,ldum,LRDTOT,
     &               (MMBST(ISYM),ISYM=1,MSYM),
     &               (NUMCHO(ISYM),ISYM=1,MSYM),
     &               (IADRTO(ISYM),ISYM=1,MSYM),
     &               (NUMFIL(ISYM),ISYM=1,MSYM)
C
         IF (MSYM .NE. NSYM) THEN
            WRITE(LUPRI,*) 'Error in symmetry checking at restart'
            STOP 'Error in symmetry checking at restart'
         END IF
         DO I = 1,NSYM
            IF (MMBST(I) .NE. NNBST(I)) THEN
               WRITE(LUPRI,*) 'Error in basis checking at restart'
               STOP 'Error in basis checking at restart'
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
         DO ISYM = 1,NSYM
            READ(LUSEC) (IDNTLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
         END DO
C
      END IF
C
  777 CONTINUE
C
C
      CALL CC_CHOLESKY(WORK(KDIAG),WORK(KOLD),WORK(KIND1),
     &                 WORK(KIND2),WORK(KEND),LWRK)
C
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'Integral decomposition completed'
      WRITE(LUPRI,*)
C
C     Poor man check result.
C
      IF (.NOT. SKPCHK) THEN
         WRITE(LUPRI,'(//,A,//)') 'Poor man check of decomposition'
         CALL CHO_RESTART(WORK(KDIAG),WORK(KEND),LWRK)
      END IF
C
C     Write restart files
C
      CALL CHO_SAVE(WORK(KIND1),WORK(KIND2))
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
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
C
      RETURN
      END
C
C  /* Deck cc_choinp */
      SUBROUTINE CC_CHOINP
C
#include <implicit.h>
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
      PARAMETER (NTABLE = 14)
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
     &            '.THRSTR','.REDUCE','.SKIP  ','.SKPCHK'/
C
C
      IF (SET) RETURN
      SET = .TRUE.
C
      CALL AROUND('Output from CHOLE_INP')
C
C---------------------
C     Initializations.
C---------------------
C
      THINDI = 1.0D+3
      THSUDI = 1.0D+3
      THRDEF = 1.0D-8
      SPAN   = 1.0D-3
C
      LBUF = 250000
C
      RSTDIA = .FALSE.
      RSTCHO = .FALSE.
      REDUCE = .FALSE.
      COMP   = .FALSE.
      SCDIAG = .TRUE.
C
      SKPDEC = .FALSE.
      SKPCHK = .FALSE.
C
      THRCHA = .FALSE.
      THSCH1 = .FALSE.
      THSCH2 = .FALSE.
      CHLBUF = .FALSE.
      SPANDI = .FALSE.
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
               GOTO(1,2,3,4,5,6,7,8,9,10,11,12) I
            END IF
         END DO
         WRITE(LUPRI,*) 'Illegal keyword in input: ',WORD
         STOP 'Illegal keyword in input'
      ELSE
         WRITE(LUPRI,*) 'Illegal prompt in input'
         STOP 'Illegal prompt in input'
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
         LBFINP = MAX(LBUF,LBFINP)
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
         SKPDEC = .TRUE.
      GOTO 200
C
C     reduce
C
   14 CONTINUE
         SKPCHK = .TRUE.
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
      IF (SKPDEC) THEN
         WRITE(LUPRI,'(A)') 'Cholesky decomposition skipped in this run'
         RETURN
      END IF
      IF (REDUCE) THEN
          WRITE(LUPRI,'(A)') 'Cholesky vectors written in reduced set'
      ELSE
          WRITE(LUPRI,'(A)') 'Cholesky vectors written in complete set'
      END IF
      IF (SCDIAG) THEN
         IF (THSCH1) THEN
            WRITE(LUPRI,'(2A,D14.4)') 'First screening threshold on ',
     &                            'diagonal defined by user',THINDI
         ELSE
            WRITE(LUPRI,'(2A,D14.4)') 'Default first screening ',
     &                            'threshold on diagonal',THINDI
         END IF
         IF (THSCH2) THEN
            WRITE(LUPRI,'(2A,D14.4)') 'Second screening threshold on ',
     &                            'diagonal defined by user',THSUDI
         ELSE
            WRITE(LUPRI,'(2A,D14.4)') 'Default second screening ',
     &                            'threshold on diagonal',THSUDI
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
         WRITE(LUPRI,'(2A,D14.4)') 'Convergence threshold defined ',
     &                         'by user',THRDEF
      ELSE
         WRITE(LUPRI,'(A,D14.4)') 'Default convergence threshold',
     &                         THRDEF
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
         WRITE(LUPRI,'(A,D14.4)') 'Special restart threshold', THRSTR
      END IF
      IF (SPANDI) THEN
         WRITE(LUPRI,'(A,D14.4)') 'Span factor: ',SPAN
      ELSE
         WRITE(LUPRI,'(A,D14.4)') 'Default span factor: ',SPAN
      END IF
      WRITE(LUPRI,'(A,I5,I20)') 'File size (Gb/words) :',IGB,MAXLEN
      IF (.NOT. CCMODSK) THEN
         WRITE(LUPRI,'(A,1P,D15.6)')
     &   'SCF density will be decomposed; threshold: ',THRDC
      ENDIF
C
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
#if defined (VAR_CHO2)
      IGB = 2
#else
#if defined (VAR_CHO4)
      IGB = 4
#else
#if defined(VAR_CHO8)
      IGB = 8
#else
#if defined(VAR_CHO32)
      IGB = 32
#else
#if defined(VAR_CHO64)
      IGB = 64
#else
#if defined(VAR_CHOTERA)
      IGB = 1024
#else
      IGB = 16
#endif
#endif
#endif
#endif
#endif
#endif
C
C     maxlen = N Gb / 8 bytes = N * 1024**3 / 8 
C     (minus security margin, probably not needed)
C
      MAXLEN = IGB * 1024*1024
      MAXLEN = MAXLEN * 1024 / 8 - MAXDIA
c     maxlen = maxlen / 64 / 16
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
c     IF (IPRINT .GT. 5) THEN
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
c     END IF
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
      DIMENSION CHOVEC(NDIM),scrcd(2,2)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C
      stop 'non-implemented choscr called'
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
      RETURN
      END
C
C/* Deck cho_restart */
C
      SUBROUTINE CHO_RESTART(DIAG,WORK,LWORK)
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
C
C
      CALL AROUND('Summary of Cholesky restart')
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
      WRITE(LUPRI,'(//,A,I10,//)') 'Number of diagonal elements :',NDIAG
      CALL FLSHFO(LUPRI)
C
      WRITE(LUPRI,'(//,A,//)') 'Analysis of read diagonal'
      CALL CC_ANADI(DIAG,.TRUE.,.FALSE.,.FALSE.)
C
      DO ISYMAB = 1,NSYM
C
C        Initializations.
C        ----------------
C
         NDIM = NNBST(ISYMAB)
         NVEC = NUMCHO(ISYMAB)
c        write(lupri,*) 'In restart, numcho:',nvec
         IF (NVEC .EQ. 0) GOTO 100
C
C        Allocate memory.
C        ----------------
C
         KDIAG = 1
         KEND1 = KDIAG + NDIM
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) STOP 'Not enough space in CHO_RESTART.1'
C
         CALL DZERO(WORK(KDIAG),NDIM)
C
C        Read and substract previous vectors.
C        ------------------------------------
C
         LWRKX = LWRK1 - NDIM
         LENAB = NDIM  + (NDIM-1)/IRAT + 2 + NDIM
C
         NRDMX = LWRKX / LENAB
c        nrdmx = min(300,nrdmx)
         IF (NRDMX .LT. 1) STOP 'Allocation failure in CHO_RESTART'
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
            LCHO2 = NREAD * (NDIM  + (NDIM-1)/IRAT + 2) + NDIM
C
            KCHO1 = KEND1
            KCHO2 = KCHO1 + LCHO1
            KEND2 = KCHO2 + LCHO2
C
            LWRK2 = LWORK - KEND2
            IF (LWRK2 .LT. 0) STOP 'Not enough space in CHO_RESTART.2'
C
            CALL CC_RDCHO(WORK(KCHO1),NREAD,ISYMAB,
     &                    IFIRST,WORK(KCHO2),LCHO2)
C
            DO I = IFIRST,ILAST
               II = I - IFIRST + 1
c              write(lupri,*)
c              write(lupri,*) 'In restart, vector no. ',ii
c              do kk = 1,ndim
c                 koff = kcho1 + ndim*(ii-1) + kk - 1
c                 write(lupri,*) kk,work(koff)
c              end do
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
         WRITE(LUPRI,'(/,A,I3)') 'Check information for symmetry',ISYMAB
         call flshfo(lupri)
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
         WRITE(LUPRI,'(/,A,I8,/)') 'Number of vectors      :',NVEC
         WRITE(LUPRI,'(2(A,D18.8,/))') 'Maximum diagonal       :',XM,
     &                             'Minimum diagonal       :',YM
         WRITE(LUPRI,'(A,2D18.8,/)')'Maximum absolute error :',
     &                             ERRMX1,EXACT1
         WRITE(LUPRI,'(3(A,I8,/))') 'Converged   diagonals :',ICONV,
     &                          'Zeroed neg  diagonals :',INEG,
     &                          'Unconverged diagonals :',IUNCO
         CALL FLSHFO(LUPRI)
C
  100    CONTINUE
C
      END DO
C
      WRITE(LUPRI,'(//,2(9X,A,/))') 'Histogram of diagonal elements',
     &                          '------------------------------'
      CALL CC_ANADI(DIAG,.TRUE.,.FALSE.,.FALSE.)
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
      INTEGER*8 ITEMP
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
            LUCHO = 30 + I
            NAME  = ROOT(1:7)//CSYM//UNDERL//CFIL
            CALL WOPEN2U(LUCHO,NAME,64,0)
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
         LUCHO  = 30 + NTOTFI
         NUMFIL(I)   = NUMFIL(I) + 1
         IDNTLU(J,I) = LUCHO
C
         NAME  = ROOT(1:7)//CSYM//UNDERL//CFIL
         CALL WOPEN2U(LUCHO,NAME,64,0)
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
               CALL WOPEN2U(LUCHO,NAME,64,0)
C
            END DO
         END DO
C
      ELSE
C
         WRITE(LUPRI,'(3A,/,A)') 'Task ',TASK,
     &                        ' not implemented in CHO_NAME',
     &                        'Program will stop'
         STOP 'Error in CHO_NAME'
C
      END IF
C
      RETURN
      END
C
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
         STOP 'Insufficient space for allocation in CC_RDAO'
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
         STOP 'Insufiicient work space in CC_RDAO'
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
#include <ibtfun.h>
#include <inftap.h>
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
#if defined (SYS_CRAY) || defined (SYS_T90)
            LRECL = LBFINP + LBFINP + 1
#endif
#if defined (SYS_AIX) || defined (SYS_SUN) || defined (SYS_HPUX) || defined (SYS_HAL) || defined (SYS_LINUX)
            IF (NBASIS .LE. 255) THEN
               LRECL = 2*LBFINP + LBFINP + 1
            ELSE
               LRECL = 2*LBFINP + 2*LBFINP + 1
            END IF
#endif
#if defined (SYS_IRIX) || defined (SYS_DEC)
            IF (NBASIS .LE. 255) THEN
               LRECL = 2*LBFINP + LBFINP + 1
            ELSE
               LRECL = 2*LBFINP + 2*LBFINP + 1
            END IF
#endif
#if defined (SYS_NEC)
            LRECL = LBFINP + LBFINP/IRAT + 1
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
 2000 STOP 'Error reading AOTWODIS'
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
     &                     IRE2CO,WORK,LWORK)
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
      DIMENSION IRE2CO(*)
C
      PARAMETER (MXCHRD = 200)
      DIMENSION IDNTC1(MXCHRD)
C
      CHARACTER*10 NAME
C
C
C     Initial stuff
C
      IF (NWRITE .EQ. 0) RETURN
C
      NDIM   = NREDUC(ISYMAB)
C
      IFIRST = NUMCHO(ISYMAB) + 1
      ILAST  = IFIRST + NWRITE - 1
C
      NUMCH0 = NUMCHO(ISYMAB) + NWRITE
      IF (NUMCH0 .GT. MXCHVC) THEN
         WRITE(LUPRI,'(A,I3)') 'Too many vectors for symmetry',ISYMAB
         WRITE(LUPRI,'(2(A,I7,/))') 'Required :', NUMCH0,
     &                          'Maximum  :', MXCHVC
         STOP 'Too many Cholesky vectors'
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
c        if (isymab .eq. 4) write(lupri,*) 'dumping vector ',i
         THRSEL = 0.0D0
         CALL CHO_PACK(THRSEL,NDIM,CHOVEC(KOFF),WORK(KEND),
     &                 LWRK,LENGTH,NUMEL,ISYMAB,IRE2CO)
c         write(lupri,*) 'isymab,i,numel,kend,work,lwrk',
c    &                    isymab,i,numel,kend,work(kend),lwrk
C
         II = I - IFIRST + 1
         IDNTCH(I,ISYMAB) = IDNTC1(II)
         LENCHO(I,ISYMAB) = LENGTH 
         IADRTO(ISYMAB)   = IADRTO(ISYMAB) + LENGTH
C
c        knum = int(work(kend))
c        kvec = kend + (knum-1)/2 + 2
c           if (isymab .eq. 4) then
c           if (numcho(isymab) .gt. 12) then
c           write(lupri,*) 'Vector in putcho of diagonal',idntc1(ii)
c           do kk = 1,knum
c              write(lupri,'(i4,d20.10,i6)') kk,work(kvec+kk-1)
c           end do
c           call flshfo(lupri)
c           end if
c           end if
c
         KOFF = KOFF + NDIM
         KEND = KEND + LENGTH
         LWRK = LWRK - LENGTH
C
         IF (LWRK .LT. 0) STOP 'Not enough space in CC_PUTCHO'
C
         LENTOT = LENTOT + LENGTH
         JFINIS = IADR   + LENTOT - 1
C
         IF (JFINIS .GT. MAXLEN) THEN
            LENTOT = LENTOT - LENGTH
            CALL PUTWA2U(LUCHO,NAME,WORK,IADR,LENTOT)
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
      CALL PUTWA2U(LUCHO,NAME,WORK,IADR,LENTOT)
c     write(lupri,*) 'Norm of compressed vector in putcho',
c    &      ddot(lentot,work,1,work,1)
C
C     Restart info
C
      REWIND(LURST)
      WRITE(LURST) NSYM,THRCOM,COMP,REDUCE,LRDTOT,
     &             (NNBST(ISYM),ISYM=1,NSYM),
     &             (NUMCHO(ISYM),ISYM=1,NSYM),
     &             (IADRTO(ISYM),ISYM=1,NSYM),
     &             (NUMFIL(ISYM),ISYM=1,NSYM)
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
C/* Deck cc_getcho */
      SUBROUTINE CC_GETCHO(CHOVEC,NREAD,ISYMAB,IFIRST,WORK,LWORK,ICO2RE)
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
      DIMENSION ICO2RE(*)
C
      CHARACTER*10 NAME
      INTEGER*8 MDIM,LONG1
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
      JADR = 1
      DO I = 1,IFIRST-1
         MDIM = LENCHO(I,ISYMAB)
         JADR = JADR + MDIM
      END DO
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
C     Read compressed vectors.
C
      LENGT2 = 0
      LENTOT = 0
c     write(lupri,*) 'Reading vector in cc_getcho. isymab :',
c    &                isymab
      DO I = IFIRST,ILAST
         LENTOT = LENTOT + LENCHO(I,ISYMAB)
         JFINIS = IADR   + LENTOT
         IF (JFINIS .GT. MAXLEN) 
     &      LENGT2 = LENGT2 + LENCHO(I,ISYMAB)
c        write(lupri,*) 'i,lencho(i,isymab),lentot:',
c    &                   i,lencho(i,isymab),lentot
      END DO
C
      KREAD = 1
      KEND1 = KREAD + LENTOT
      LWRK1 = LWORK - KEND1 + 1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Not enough space in CC_GETCHO.1'
         WRITE(LUPRI,'(A,I10)') 'Available memory :', LWORK
         WRITE(LUPRI,'(A,I10)') 'Needed    memory :', KEND1
         STOP 'Not enough space in cc_getcho'
      END IF
C
      LENGTH = NREDUC(ISYMAB)*NREAD
      CALL DZERO(CHOVEC,LENGTH)
C
      IF (LENGT2 .EQ. 0) THEN
         CALL GETWA2U(LUCHO,NAME,WORK(KREAD),IADR,LENTOT)
      ELSE
C
         LENGT1 = LENTOT - LENGT2
         CALL GETWA2U(LUCHO,NAME,WORK(KREAD),IADR,LENGT1)
C
         IFILE = IFILE + 1
         CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'OLD')
C
         IADR = 1
         KOFF = KREAD + LENGT1
         CALL GETWA2U(LUCHO,NAME,WORK(KOFF),IADR,LENGT2)
C
      END IF
c
c     write(lupri,*) 'Norm of compressed vector in getcho',
c    &               ddot(lentot,work(kread),1,work(kread),1)
C
C     Unpack vectors
C
      KOFF1 = KREAD
      KOFF2 = 1
      LENG2 = NREDUC(ISYMAB)
      DO I = IFIRST,ILAST
C
         LENG1 = LENCHO(I,ISYMAB)
c
c           if (isymab .eq. 4) then
c        knum = int(work(koff1))
c        kvec = koff1 + (knum-1)/2 + 2
c        write(lupri,*) 'ifirst,knum',ifirst,knum
c           if (i .gt. 10) then
c           write(lupri,*) 'Vector in getcho of diagonal',IDNTCH(I,4)
c           do kk = 1,knum
c              write(lupri,'(i4,d20.10,i6)') kk,work(kvec+kk-1)
c           end do
c           call flshfo(lupri)
c           end if
c           write(lupri,*) 'vec,work(kvec)',i,work(kvec),work(kvec-1),
c    &                      work(kvec+1)
c           end if
c
         CALL CHO_UNPACK(ISYMAB,WORK(KOFF1),LENG1,CHOVEC(KOFF2),
     &                   LENG2,ICO2RE)
C
c           if (isymab .eq. 4) then
c           if (numcho(isymab) .gt. 10) then
c           write(lupri,*) ' Vector unpacked',i
c           do iii = 1,nreduc(isymab)
c              kk = koff2 + iii - 1
c              write(lupri,*) iii, chovec(kk)
c           end do
c           end if
c           end if
         KOFF1 = KOFF1 + LENG1
         KOFF2 = KOFF2 + LENG2
C
      END DO
C
      RETURN
      END
C
C  /* Deck cc_anadi */
      SUBROUTINE CC_ANADI(DIAG,ALL,CAUCHY,MKRED)
C
C     Alfredo Sanchez   Julio 2002
C
C     Histogram on diagonal.
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DIAG(*)
      LOGICAL ALL,CAUCHY,MKRED
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
      WRITE(LUPRI,'(A,2I10)') 'Larger than 1.0D2           : ',NUMH,NUMT
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
      WRITE(LUPRI,'(/,A,2I10,//)') 'Negative diagonals          : ',NUMN
C
      IF (.NOT. CAUCHY) GOTO 888
C
C
C     Analyse diagonal by Cauchy-Schwarts.
C     ------------------------------------
C
      CALL AROUND('Cauchy-Schwarts screening based on initial diagonal')
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
C
  888 CONTINUE
      IF (.NOT. MKRED) GOTO 999
C
      DO ISYM = 1,NSYM
C
         NRE = 0
         DO I = 1,NNBST(ISYM)
            KOFF1 = IDIAG(ISYM) + I
            IF (DIAG(KOFF1) .NE. ZERO) NRE = NRE + 1
         END DO
C
         NREDUC(ISYM) = NRE
C
      END DO

   
 999  CONTINUE
      RETURN
      END
C
C
C  /* Deck cho_chkdia */
      SUBROUTINE CHO_CHKDIA(DIAG,DIAOLD,OLDIAG,XC,ICAB,ICHOV,
     &                     XMAX1,ISYMAB)
C
C     Alfredo Sanchez   Julio 2002
C
C     Check diagonal after updating.
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DIAG(*),DIAOLD(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <choio.h>
#include <priunit.h>
C
C
      ICHOV = ICHOV + 1
C
C     Loop over NNBST, Should be modified to NREDUC when tested
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
      WRITE(LUPRI,'(I3,I8,I10,2D14.5,I8,I6,D14.5)')
     &         ISYMAB,ICHOV,ICAB,XC,OLDIAG,ICONV,INEG,XMAX1
c     WRITE(LUPRI,'(I3,I8,I10,2D25.16,I8,I6,D25.16)')
c    &         ISYMAB,ICHOV,ICAB,XC,OLDIAG,ICONV,INEG,XMAX1
c     if (isymab .eq. 1)
c    &      write(lupri,'(a,2d25.16)')'chkdia 65453 65454',
c    &                             diag(65453),diag(65454)
      CALL FLSHFO(LUPRI)
C
C     Write secure restart file
C
      IF (MOD(ICHOV,500) .EQ. 0) THEN
C
         REWIND(LUSEC)
         WRITE(LUSEC) NSYM,THRCOM,COMP,REDUCE,LRDTOT,
     &             (NNBST(ISYM),ISYM=1,NSYM),
     &             (NUMCHO(ISYM),ISYM=1,NSYM),
     &             (IADRTO(ISYM),ISYM=1,NSYM),
     &             (NUMFIL(ISYM),ISYM=1,NSYM)
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
      WRITE(LUPRI,'(//,2(9X,A,/))') 'Number of vectors',
     &                          '-----------------'
      WRITE(LUPRI,'(6X,A,8X,A,7X,A)') 'ISAB','Maximum','Actual'
      DO ISYM = 1,NSYM
         WRITE(LUPRI,'(7X,I2,2I16)') ISYM,NNBST(ISYM),NUMCHO(ISYM)
      END DO
C
      WRITE(LUPRI,'(//,2(9X,A,/))') 'Number of elements',
     &                              '------------------'
      WRITE(LUPRI,'(6X,A,8X,A,3(7X,A))') 'ISAB','Reduced','Real 0',
     &                                   'Zeroed','Non zero'
      DO ISYM = 1,NSYM
         XELEM = DFLOAT(NUMCHO(ISYM))*DFLOAT(NREDUC(ISYM))
         YELEM = XELEM - YSEL0(ISYM) - XSEL0(ISYM)
         WRITE(LUPRI,'(7X,I2,4F16.0)') ISYM,XELEM,YSEL0(ISYM),
     &                                 XSEL0(ISYM),YELEM
      END DO
C
      WRITE(LUPRI,'(/,A,D15.5)') 
     *     'Largest absolute element in diagonal:',DIAG(IY)
C
      MAXINT = MAXSHL * (MAXSHL+1) / 2
      WRITE(LUPRI,'(/,A,2I10,//)')
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
      WRITE(LUPRI,'(6X,A,F10.2)')    '         Selecting  : ', TSELTO
      WRITE(LUPRI,'(6X,A,F10.2)')    '           Reading  : ', TRDTO
      WRITE(LUPRI,'(6X,A,F10.2)')    '         Preparing  : ', TPRTO
      WRITE(LUPRI,'(6X,A,F10.2)')    '       Multiplying  : ', TDGTO
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
C
C
      REWIND(LUSEC)
C
      WRITE(LUSEC) NSYM,THRCOM,COMP,REDUCE,LRDTOT,
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
      WRITE(LURST) NSYM,THRCOM,COMP,REDUCE,LRDTOT,
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
C
C  /* Deck cho_reduce */
      SUBROUTINE CHO_REDUCE(DIAG,ICO2RE,IRE2CO)
C
C     Get reduced to complete and complete to reduce set indices
C
C     ASM & VPM     April 2003
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0)
      DIMENSION DIAG(*),IRE2CO(*),ICO2RE(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <ccdeco2.h>
C
C       
      DO ISYM = 1,NSYM
C
         ICOUN1 = 0
         DO I = 1,NREINI(ISYM)
            KOFF = IREDIA(ISYM) + I
            IF (ABS(DIAG(KOFF)) .GT. ZERO) THEN
               ICOUN1 = ICOUN1 + 1
               ICO2RE(IREDIA(ISYM) + I) = ICOUN1
               IRE2CO(IREDIA(ISYM) + ICOUN1) = I
            ELSE
               ICO2RE(IDIAG(ISYM) + I) = NREDUC(ISYM) + 1
            END IF 
         END DO
      END DO
C
      RETURN
      END
C /* Deck cho_redvec */
      SUBROUTINE CHO_REDVEC(FULL,REDU,NUMVEC,IRE2CO,ISYM)
C
C     Get reduced vectors from full integral columns.
C
C     ASM & VPM     April 2003
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0)
      DIMENSION FULL(*),REDU(*),IRE2CO(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
C
C
      KOFF2 = 0
      DO IVEC = 1,NUMVEC
         KOFF1 = NNBST(ISYM) * (IVEC-1)
         DO I = 1,NREDUC(ISYM)
            II = IRE2CO(IDIAG(ISYM) + I)
            KOFF2 = KOFF2 + 1
            REDU(KOFF2) = FULL(KOFF1+II)
         END DO
      END DO
C
      RETURN
      END
C
C
C /* Deck cho_pack */
C
      SUBROUTINE CHO_PACK(THRSEL,NDIM,CHOVEC,WORK,LWORK,LENGTH,NUM,ISYM,
     &                    IRE2CO)
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
      DIMENSION IRE2CO(*)
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
         STOP 'Not enough space in CHO_PACK'
      END IF
C
      CALL CHO_PACK1(THRSEL,NDIM,CHOVEC,WORK(KCHISE),WORK(KCHSAV),
     &               NUM,.FALSE.,ISYM,IRE2CO)
C
      IF (LWRK .LT. NUM) THEN
         WRITE(LUPRI,*) 'Not enough space in CHO_PACK.2'
         WRITE(LUPRI,'(A,2I10)')'Available,needed:',LWRK,NUM
         STOP 'Not enough space in CHO_PACK'
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
     &                     NUMEL,INFO,ISYM,IRE2CO)
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
#include <ccdeco.h>
      DIMENSION CHOVEC(NDIM),CHOSEL(NDIM),ISELEC(NDIM)
      DIMENSION IRE2CO(*)
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
            ISELEC(ICOUNT) = IRE2CO(IDIAG(ISYM) + I)
            CHOSEL(ICOUNT) = CHOVEC(I)
c           write(lupri,*) 'In pack1, iselec,chosel',icount,
c    &                        iselec(icount),chosel(icount)
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
C /*Deck cho_unpack */
C
      SUBROUTINE CHO_UNPACK(ISYM,SCR,LSCR,CHOVEC,NDIM,ICO2RE)
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
      DIMENSION ICO2RE(*)
C
C
      NUMEL = INT(SCR(1))
C
      KCHISE = 2
      KCHSAV = KCHISE + (NUMEL - 1)/IRAT + 1
C
      CALL CHO_UNPACK1(ISYM,CHOVEC,SCR(KCHISE),SCR(KCHSAV),
     &                 NUMEL,NDIM,ICO2RE)
C
      RETURN
      END
C
      SUBROUTINE CHO_UNPACK1(ISYM,CHOVEC,ISELEC,CHOSEL,NUMEL,
     &                       NDIM,ICO2RE)
C
C     Written by Alfredo Sanchez   julio 2002
C
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccdeco.h>
C
      DIMENSION CHOVEC(NDIM),CHOSEL(NUMEL),ISELEC(NUMEL)
      DIMENSION ICO2RE(*)
C
C
c     write(lupri,*) 'Vector in unpack1',' red,val,full,chosel'
c     write(lupri,*) 'In cho_unpack1,numel',numel
      DO I = 1,NUMEL
         J = ICO2RE(IDIAG(ISYM) + ISELEC(I))
c        write(lupri,*)'i,iselec(i),ico2re(idiag(isym)+iselec(i))',
c    &                 'chosel(i) : ',i,iselec(i),
c    &                  ico2re(idiag(isym)+iselec(i)),chosel(i)
         CHOVEC(J) = CHOSEL(I)
c        write(lupri,'(2(i4,d20.10))') j,chovec(j),iselec(i),chosel(i)
      END DO
C
      RETURN
      END
C
C/* Deck cc_rdcho */
      SUBROUTINE CC_RDCHO(CHOVEC,NREAD,ISYMAB,IFIRST,WORK,LWORK)
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
C
      CHARACTER*10 NAME
      INTEGER*8 MDIM,LONG1
C
C
c     write(lupri,*) 'entering rdcho'
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
      JADR = 1
      DO I = 1,IFIRST-1
           MDIM = LENCHO(I,ISYMAB)
           JADR = JADR + MDIM
      END DO
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
C     Read compressed vectors.
C
      LENGT2 = 0
      LENTOT = 0
      DO I = IFIRST,ILAST
         LENTOT = LENTOT + LENCHO(I,ISYMAB)
c        write(lupri,*)'vector,lenght,lentot',i,lencho(i,isymab),lentot
c        call flshfo(lupri)
         JFINIS = IADR   + LENTOT
         IF (JFINIS .GT. MAXLEN) 
     &      LENGT2 = LENGT2 + LENCHO(I,ISYMAB)
      END DO
C
      KREAD = 1
      KEND1 = KREAD + LENTOT
      LWRK1 = LWORK - KEND1 + 1
C
      IF (LWRK1 .LT. 0) STOP 'Not enough space in cc_getcho'
C
      LENGTH = NREAD * NNBST(ISYMAB)
      CALL DZERO(CHOVEC,LENGTH)
C
      IF (LENGT2 .EQ. 0) THEN
         CALL GETWA2U(LUCHO,NAME,WORK(KREAD),IADR,LENTOT)
      ELSE
C
         LENGT1 = LENTOT - LENGT2
         CALL GETWA2U(LUCHO,NAME,WORK(KREAD),IADR,LENGT1)
C
         IFILE = IFILE + 1
         CALL CHO_NAME(ISYMAB,IFILE,LUCHO,NAME,'OLD')
C
         IADR = 1
         KOFF = KREAD + LENGT1
         CALL GETWA2U(LUCHO,NAME,WORK(KOFF),IADR,LENGT2)
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
         CALL CHO_UNCOMP(WORK(KOFF1),LENG1,CHOVEC(KOFF2),LENG2)
C
c        write(lupri,*) 'Vector uncompressed',koff2,chovec(koff2)
c        call flshfo(lupri)
         KOFF1 = KOFF1 + LENG1
         KOFF2 = KOFF2 + LENG2
C
      END DO
C
      RETURN
      END
C
      SUBROUTINE CHO_UNCOMP(SCR,LSCR,CHOVEC,NDIM)
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
      CALL CHO_UNCOMP1(CHOVEC,SCR(KCHISE),SCR(KCHSAV),
     &                 NUMEL,NDIM)
C
      RETURN
      END
C
      SUBROUTINE CHO_UNCOMP1(CHOVEC,ISELEC,CHOSEL,NUMEL,NDIM)
C
C     Written by Alfredo Sanchez   julio 2002
C
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION CHOVEC(NDIM),CHOSEL(NUMEL),ISELEC(NUMEL)
C
C
c     write(lupri,*) 'In uncomp1,numel:',numel
      call flshfo(lupri)
      DO I = 1,NUMEL
         J = ISELEC(I)
         CHOVEC(J) = CHOSEL(I)
      END DO
C
      RETURN
      END
C
C/* Deck cho_redini */
C
      SUBROUTINE CHO_REDINI(INDRED,DIAG,DIAOLD,WORK,LWORK)
C
C     LENRED : Real length of screened (reduced) diagonal
C     LRDTOT : Length allocated to INDRED 
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0)
      DIMENSION DIAG(*), DIAOLD(*), INDRED(*)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <ccdeco2.h>
C
      LENRED = 0
      KOFF3  = 0
      DO ISYM = 1,NSYM
         ICOUNT = 0
         DO I = 1,NNBST(ISYM)
            KOFF1 = IDIAG(ISYM) + I
            IF (DIAG(KOFF1) .NE. 0) THEN
C
               ICOUNT = ICOUNT + 1
C
               KOFF2 = LREDU*(ISYM-1) + ICOUNT
               INDRED(KOFF2) = I
C
               KOFF3 = KOFF3 + 1
               DIAOLD(KOFF3) = DIAG(KOFF1)
C
            END IF
         END DO
C
         NREINI(ISYM) = ICOUNT
         LENRED = LENRED + NREINI(ISYM)
C
      END DO
C
      WRITE(LUPRI,'(/,A,8I8)') 'NREDUC :',(NREDUC(I),I=1,NSYM)
      WRITE(LUPRI,'(A,I16,/)') 'LENRED :',LENRED
C
      LRDTOT = LREDU * NSYM
C
      RETURN
      END
