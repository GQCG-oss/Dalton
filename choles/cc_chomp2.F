C  /* Deck cc_chomp2 */
      SUBROUTINE CC_CHOMP2(WORK,LWORK,EMP2)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C
C        Calculate the MP2 energy correction (EMP2) using Cholesky
C        decomposed two-electron integrals.
C
C     Input variable IALMP2 (chomp2.h) decides the algorithm:
C
C        IALMP2 < 1: Same as IALMP2>3 below.
C
C        IALMP2 = 1: Calculate EMP2 directly from Cholesky vectors
C                    without permanently storing (ai|bj) integrals:
C                    batch over b-index only.
C
C        IALMP2 = 2: Calculate EMP2 through construction of (ai|bj)
C                    integrals (stored in memory).
C
C        IALMP2 = 3: Calculate EMP2 directly from Cholesky vectors
C                    without permanently storing (ai|bj) integrals:
C                    batch over a- and b-indices.
C
C        IALMP2 > 3: Use the IALMP2=2 algorithm if sufficient memory,
C                    otherwise use IALMP2=1 algorithm.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccorb.h>
#include <dccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <priunit.h>
#include <iratdef.h>

      INTEGER MINVEC(8), NVECTOT(8)

      PARAMETER (MINDEF = 5, INFO = 0, INFT = 5)
      PARAMETER (IOPEN = -1, IDEL = 0, IKEEP = 1, IMOAO = 1, IMOMO = 2)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)

      LOGICAL POSSBL

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOMP2')

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Start timing.
C     -------------

      TIMTOT = SECOND()

C     Check that this is a Cholesky run.
C     ----------------------------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,/)')
     &   'FATAL ERROR IN ',SECNAM,':',
     &   '- integrals must be Cholesky decomposed!'
         CALL QUIT('Fatal error in '//SECNAM)
      ENDIF

C     Read the MO coefficient matrix and orbital energies.
C     Reorder to appropriate CC ordering.
C     If requested, allocate space for amplitude diagonal.
C     NB!! Memory is reset below, so be careful if you
C          change anything in this allocation.
C     ----------------------------------------------------

      LCMO  = NLAMDS
      IF (SKIPTR) LCMO = 0
      NTOT1 = 0
      IF (MP2SAV) THEN
         DO ISYMAI = 1,NSYM
            NTOT1 = NTOT1 + NT1AM(ISYMAI)
         ENDDO
      ENDIF
      KEMO  = 1
      KTDIA = KEMO  + NORBTS
      KCMO  = KTDIA + MAX(NTOT1,1)
      KEND1 = KCMO  + LCMO
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,
     &   ' - Allocation: MO'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KEMO),WORK(KCMO),WORK(KEND1),LWRK1,
     &               .FALSE.,.TRUE.,(.NOT.SKIPTR))

C     Transform AO Cholesky vectors to MO (ai) basis.
C     Store resulting vectors on disk.
C     -----------------------------------------------

      TIMTRF = SECOND()
      CALL CHO_TRFAI(WORK(KCMO),WORK(KEND1),LWRK1)
      TIMTRF = SECOND() - TIMTRF

C     Reset memory; CMO no longer needed.
C     -----------------------------------

      KEND1 = KCMO
      LWRK1 = LWORK - KEND1 + 1

C     Decompose (ai|bj) integral matrix.
C     ----------------------------------

      IF (CHOMO) THEN
         TIMDEC = SECOND()
         CALL CC_DECOMP2(WORK(KEND1),LWRK1)
         TIMDEC = SECOND() - TIMDEC
      ENDIF

C     Get max. elements in Cholesky vectors.
C     --------------------------------------

      IF (SPRMP2) THEN
         NTOT = 0
         MXT1 = -1
         IF (CHOMO) THEN
            DO ISYM = 1,NSYM
               NTOT = NTOT + NCHOMO(ISYM)
               MXT1 = MAX(MXT1,NT1AM(ISYM))
            ENDDO
         ELSE
            DO ISYM = 1,NSYM
               NTOT = NTOT + NUMCHO(ISYM)
               MXT1 = MAX(MXT1,NT1AM(ISYM))
            ENDDO
         ENDIF
         KMAX   = KEND1
         KMAPAI = KMAX   + NTOT
         KMAI   = KMAPAI + (MXT1 - 1)/IRAT + 1
         KMAPBJ = KMAI   + (MXT1 - 1)/IRAT + 1
         KMBJ   = KMAPBJ + (MXT1 - 1)/IRAT + 1
         KEND1  = KMBJ   + (MXT1 - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) CALL QUIT('Insuf. mem. in '//SECNAM)
         IF (CHOMO) THEN
            CALL CC_CHOMAX(WORK(KMAX),WORK(KEND1),LWRK1,NCHOMO,2)
         ELSE
            CALL CC_CHOMAX(WORK(KMAX),WORK(KEND1),LWRK1,NUMCHO,1)
         ENDIF
         LALG  = 3
         KEND2 = KEND1
         LWRK2 = LWRK1
         GO TO 101
      ENDIF

C     Figure out which algorithm to use.
C     ----------------------------------

      IF (IALMP2 .EQ. 1) THEN

C        Force algorithm 1.
C        ------------------

         LALG = IALMP2

         KEND2 = KEND1
         LWRK2 = LWRK1

         XWRK2 = ONE*LWRK2
         DO ISYMB = 1,NSYM
            IF (XCKI(ISYMB) .GT. XWRK2) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory for MP2 algorithm 1 in ',SECNAM
               WRITE(LUPRI,'(5X,A,/)')
     &         '(IALMP2 = 1 was forced; try IALMP2=0 instead)'
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF
         ENDDO

      ELSE IF (IALMP2 .EQ. 2) THEN

C        Force algorithm 2.
C        ------------------

         LALG = IALMP2

         XLWRK = ONE*LWRK1
         IF (XT2SQ(1) .GT. XLWRK) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for IALMP2=2 in ',SECNAM
            WRITE(LUPRI,'(5X,A,F12.1,/,5X,A,F12.1)')
     &      'Need (more than): ',XT2SQ(1),
     &      'Available       : ',XLWRK
            WRITE(LUPRI,'(5X,A,/)')
     &      '(IALMP2=2 was forced; try IALMP2=0 instead)'
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         KXINT = KEND1
         KEND2 = KXINT + NT2SQ(1)
         LWRK2 = LWORK - KEND2 + 1

      ELSE IF (IALMP2 .EQ. 3) THEN

         LALG = IALMP2

         IF (MP2SAV) THEN
            WRITE(LUPRI,'(/,5X,A,/,5X,A,I2,/)')
     &      '***NOTICE*** MP2 amplitudes will not be saved',
     &      '- MP2SAV is incompatible with IALMP2 =',IALMP2
            KEND1 = KTDIA
            LWRK1 = LWORK - KEND1 + 1
            MP2SAV = .FALSE.
         ENDIF

         XLWRK = ONE*LWRK1
         XRHFT = ONE*NRHFT
         XNEED = ZERO
         DO ISYMB = 1,NSYM
            DO ISYMA = 1,NSYM
               DO ISYCHO = 1,NSYM
                  ISYMJ = MULD2H(ISYCHO,ISYMB)
                  ISYMI = MULD2H(ISYCHO,ISYMA)
                  XTST  = XT1AM(ISYCHO) + XRHF(ISYMI) + XRHF(ISYMJ)
                  XNEED = MAX(XNEED,XTST)
               ENDDO
            ENDDO
         ENDDO
         XNEED = XNEED + XRHFT*XRHFT

         IF (XNEED .GT. XLWRK) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for IALMP2=3 in ',SECNAM
            WRITE(LUPRI,'(5X,A,F12.1,/,5X,A,F12.1)')
     &      'Need (at least): ',XNEED,
     &      'Available      : ',XLWRK
            WRITE(LUPRI,'(5X,A,/)')
     &      '(IALMP2=3 was forced)'
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         KEND2 = KEND1
         LWRK2 = LWRK1

      ELSE

C        Find a suitable algorithm.
C        --------------------------

         CALL GET_MP2ALG(LWRK1,LALG)

C        Set up memory indices.
C        ----------------------

         KEND2 = KEND1
         LWRK2 = LWRK1
         IF (LALG .EQ. 2) THEN
            KXINT = KEND1
            KEND2 = KXINT + NT2SQ(1)
            LWRK2 = LWORK - KEND2 + 1
         ELSE IF ((LALG.EQ.3) .AND. MP2SAV) THEN
            WRITE(LUPRI,'(/,5X,A,/,5X,A,I2,/)')
     &      '***NOTICE*** MP2 amplitudes will not be saved',
     &      '- MP2SAV is incompatible with LALG =',LALG
            KEND2 = KTDIA
            LWRK2 = LWORK - KEND2 + 1
            MP2SAV = .FALSE.
         ENDIF

      ENDIF

C     Calculate EMP2 according to LALG.
C     ---------------------------------

  101 CONTINUE

      IF (LALG .EQ. 1) THEN

         TIMMP2 = SECOND()
         CALL MP2CHO_1(WORK(KTDIA),WORK(KEMO),WORK(KEND2),LWRK2,EMP2)
         TIMMP2 = SECOND() - TIMMP2
         
      ELSE IF (LALG .EQ. 2) THEN

         TIMMP2 = SECOND()
         CALL MP2CHO_2(WORK(KTDIA),WORK(KEMO),WORK(KXINT),
     &                 WORK(KEND2),LWRK2,EMP2)
         TIMMP2 = SECOND() - TIMMP2

      ELSE IF (LALG .EQ. 3) THEN

         IF (SPRMP2) THEN
            IF (SCRMP2 .LT. ZERO) THEN
               IF (CHOMO) THEN
                  SCRMP2 = THRMP2
               ELSE
                  SCRMP2 = THRCOM
               ENDIF
            ENDIF
            TIMMP2 = SECOND()
            CALL MP2CHO_4(WORK(KEMO),WORK(KMAX),WORK(KEND2),LWRK2,EMP2,
     &                    WORK(KMAPAI),WORK(KMAI),
     &                    WORK(KMAPBJ),WORK(KMBJ),SCRMP2)
            TIMMP2 = SECOND() - TIMMP2
         ELSE
            TIMMP2 = SECOND()
            CALL MP2CHO_3(WORK(KEMO),WORK(KEND2),LWRK2,EMP2)
            TIMMP2 = SECOND() - TIMMP2
         ENDIF

      ELSE

         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,I10,A,/,5X,A,I10,A,/)')
     &   'LOGICAL BUG IN ',SECNAM,':',
     &   ' LALG   = ',LALG,' should not be possible',
     &   ' IALMP2 = ',IALMP2,' (from input)'
         CALL QUIT('Logical bug in '//SECNAM)

      ENDIF

C     Print section.
C     --------------

      IF (MP2SAV .AND. (IPRINT.GE.INFT)) THEN

         MAXT1 = 0
         DO ISYMBJ = 1,NSYM
            MAXT1 = MAX(MAXT1,NT1AM(ISYMBJ))
         ENDDO
         XT1MX = ONE*MAXT1

         XLWRK = ONE*LWORK
         XNEED = XT1AM(1) + XT2AM(1) + XT1MX
         IF (XNEED .GT. XLWRK) GO TO 1000

         CALL HEADER(SECNAM//': MP2 Amplitudes',-1)

         KT1AM = 1
         KT2AM = KT1AM + NT1AM(1)
         KSCR1 = KT2AM + NT2AM(1)
         KENDT = KSCR1 + MAXT1
         LWRKT = LWORK - KENDT + 1

         IF (LWRKT .LT. 0) THEN
            WRITE(LUPRI,'(15X,A)')
     &      '*** NOTICE: Insufficient memory for print - skipping.'
            WRITE(LUPRI,'(15X,A,I10,/,15X,A,I10,//)')
     &      'Memory required for print: ',KENDT-1,
     &      'Total memory available   : ',LWORK
            GOTO 1000
         ENDIF

         CALL DZERO(WORK(KT1AM),NT1AM(1))

         DO ISYMBJ = 1,NSYM

            CALL OP_MP2TAM(IOPEN,ISYMBJ,LUTAM,1)

            DO NBJ = 1,NT1AM(ISYMBJ) 

               CALL RD_MP2TAM(WORK(KSCR1),NT1AM(ISYMBJ),NBJ,LUTAM)

               DO NAI = 1,NBJ

                  KOFFT = KT2AM + IT2AM(ISYMBJ,ISYMBJ)
     &                  + INDEX(NAI,NBJ) - 1
                  KOFFS = KSCR1 + NAI - 1

                  WORK(KOFFT) = WORK(KOFFS)

               ENDDO

            ENDDO

            CALL OP_MP2TAM(IKEEP,ISYMBJ,LUTAM,1)

         ENDDO

         CALL CC_PRAM(WORK(KT1AM),PT1,1)

      ENDIF

 1000 CONTINUE

      IF (IPRINT .GE. INFO) THEN

         TIMTOT = SECOND() - TIMTOT

         CALL AROUND('Output from '//SECNAM)

         CALL HEADER('Summary of MP2 Calculation',-1)
         IF (LALG .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '- (ai|#bj) integrals calculated on the fly'
         ELSE IF (LALG .EQ. 2) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '- full square (ai|bj) integral matrix held in core'
         ELSE IF (LALG .EQ. 3) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '- (#ai|#bj) integrals calculated on the fly'
         ELSE
            WRITE(LUPRI,'(5X,A)')
     &      '***WARNING: POSSIBLE LOGICAL BUG DETECTED: CHECK OUTPUT !!'
         ENDIF
         IF (IALMP2 .EQ. LALG) THEN
            WRITE(LUPRI,'(5X,A,I2)')
     &      '- algorithm forced by input option .ALGORI :',IALMP2
         ENDIF
         IF (CHOMO) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '- (ai|bj) integral matrix separately Cholesky decomposed'
         ENDIF
         IF (MP2SAV) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '- MP2 amplitudes (and diagonal) saved on disk'
         ENDIF
         WRITE(LUPRI,'(5X,A,F22.12,A,/)')
     &   'MP2 energy correction       : ',EMP2,' hartree'

         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for MO transformation    : ',TIMTRF,' seconds'
	 IF (CHOMO) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time for decomposing (ai|bj)  : ',TIMDEC,' seconds'
         ENDIF
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for MP2 energy correction: ',TIMMP2,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '--------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                    : ',TIMTOT,' seconds'

      ENDIF

      RETURN
      END
C  /* Deck cc_chomax */
      SUBROUTINE CC_CHOMAX(CHMAX,WORK,LWORK,NUMCHO,ITYP)
C
#include <implicit.h>
      DIMENSION CHMAX(*), WORK(LWORK)
      INTEGER NUMCHO(8)
#include <ccorb.h>
#include <ccsdsym.h>

      ICOUNT = 0

      DO ISYCHO = 1,NSYM

         CALL CHO_MOP(-1,ITYP,ISYCHO,LUCHO,1,1)

         KVEC  = 1
         KEND1 = KVEC  + NT1AM(ISYCHO)
         LWRK1 = LWORK - KEND1 + 1
         IF (LWRK1 .LT. 0) CALL QUIT('Insuf. mem. in CC_CHOMAX')

         DO JVEC = 1,NUMCHO(ISYCHO)

            ICOUNT = ICOUNT + 1

            CALL CHO_MOREAD(WORK(KVEC),NT1AM(ISYCHO),1,JVEC,LUCHO)

            CHMAX(ICOUNT) = ABS(WORK(KVEC))
            DO NAI = 2,NT1AM(ISYCHO)
               CHMAX(ICOUNT) = MAX(CHMAX(ICOUNT),ABS(WORK(KVEC+NAI-1)))
            ENDDO

         ENDDO

         CALL CHO_MOP(1,ITYP,ISYCHO,LUCHO,1,1)

      ENDDO


      RETURN
      END
C  /* Deck mp2cho_1 */
      SUBROUTINE MP2CHO_1(T2DIA,EMO,WORK,LWORK,EMP2)
C
C     Thomas Bondo Pedersen, May 2002.
C     - revised, tbp, July 2002.
C
C     Purpose:
C        Calculate the MP2 energy correction EMP2 directly
C        from Cholesky decomposed two-electron integrals without
C        permanent storage of (ai|bj) integral intermediates.
C
C     Notes:
C        This routine splits work space in two parts and performs
C        a double batch; one over b in (ai|bj) integrals, the other
C        over the MO Cholesky vectors. This implies that, in general, the
C        Cholesky vectors will be read several times.
C
C        This routine requires (at least) in the order of
C             VO^2 + VO + O
C        words of memory; it is thus designed for calculating EMP2
C        for rather large systems and/or basis sets.
C
C        If amplitudes are not to be saved (MP2SAV=.FALSE. from input),
C        array T2DIA may be dummy.
C
#include <implicit.h>
      DIMENSION T2DIA(*), EMO(*), WORK(LWORK)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <dccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccdeco.h>
#include <chomp2.h>
Casm
#include <chomp2_b.h>
Casm
#include <priunit.h>

      INTEGER BBEG, BEND
      INTEGER LVIR(8), IOFB1(8), IX1AM(8,8), NX1AM(8), IX2SQ(8)
      INTEGER IOFFD(8), NVECTOT(8)
      INTEGER LUTAM(8), LUCHMO(8), LSYM(8)

      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.00D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
Casm
      LOGICAL LDUM
Casm

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'MP2CHO_1')

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIME = ZERO
      TIMA = ZERO

C     Initial check of memory.
C     ------------------------

      XNEED = ZERO
      DO ISYMB = 1,NSYM
         XCHO = ZERO
         DO ISYCHO = 1,NSYM
            ISYMJ = MULD2H(ISYCHO,ISYMB)
            XTST  = XT1AM(ISYCHO) + XRHF(ISYMJ)
            XCHO  = MAX(XCHO,XTST)
         ENDDO
         XTOT  = XCKI(ISYMB) + XCHO
         XNEED = MAX(XNEED,XTOT)
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',INT(XNEED)
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      XLWORK = ONE*LWORK
      IF (XNEED .GT. XLWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',INT(XNEED),
     &   'Available memory       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initializations.
C     ----------------

      EMP2   = ZERO
      IF (RSTMP2) EMP2 = OLDEN2
      MXMEMT = 0
      MXMEML = 0
      XLWORK = ONE*LWORK

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

Casm
C     Restart info
C     ------------

      IF (RSTMP2) THEN
         WRITE(LUPRI,'(//,10X,A,/,10X,A,/,10X,A,/)') 
     &                               '-----------------------',
     &                               'Restarting Cholesky MP2',
     &                               '-----------------------'
         WRITE(LUPRI,'(A,I5,A,I2,/)') 'First orbital is nr',
     &                                IFVIRB,' of symmetry',IFSYMB
         IF (OLDKNO) THEN
            WRITE(LUPRI,'(A,F15.10)') 
     &           'Contribution from previous batches : ', OLDEN2
         ELSE
            WRITE(LUPRI,'(A,A)') 'Contribution from previous batches ',
     &            'not known. You must add it by yourself!!'
         END IF
         WRITE(LUPRI,*)
      END IF 
Casm

C     Set # vectors.
C     --------------

      IF (CHOMO) THEN
         ITYP = 2
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NCHOMO(ISYM)
         ENDDO
      ELSE
         ITYP = 1
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NUMCHO(ISYM)
         ENDDO
      ENDIF

C     Index array for T2DIA.
C     ----------------------

      ICOUNT = 0
      DO ISYMAI = 1,NSYM
         IOFFD(ISYMAI) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYMAI)
      ENDDO

C     Open files for MO Cholesky vectors.
C     -----------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IOPEN,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     If amplitudes are to be saved, open files.
C     ------------------------------------------

      IF (MP2SAV) THEN
         DTIME = SECOND()
         CALL OP_MP2TAM(IOPEN,LSYM,LUTAM,NSYM)
         DTIME = SECOND() - DTIME
         TIMA  = TIMA     + DTIME
      ENDIF

C     Debug: initialize norms.
C     ------------------------

      IF (LOCDBG) THEN
         XNRM = ZERO
         TNRM = ZERO
      ENDIF

C     Print EMO.
C     ----------

      IF (IPRINT .GT. 120) THEN
         CALL AROUND('EMO in '//SECNAM)
         DO ISYM = 1,NSYM
            IF ((NRHF(ISYM).GT.0) .OR. (NVIR(ISYM).GT.0)) THEN
               WRITE(LUPRI,'(//,29X,A,I1)')
     &         'Symmetry Block ',ISYM
            ELSE
               GO TO 112
            ENDIF
            IF (NRHF(ISYM) .EQ. 0) GO TO 111
            CALL HEADER('Occupied Part',-1)
            DO I = 1,NRHF(ISYM)
               KOFFI = IRHF(ISYM) + I
               WRITE(LUPRI,'(29X,F14.10)') EMO(KOFFI)
            ENDDO
  111       IF (NVIR(ISYM) .EQ. 0) GO TO 112
            CALL HEADER('Virtual Part',-1)
            DO A = 1,NVIR(ISYM)
               KOFFA = IVIR(ISYM) + A
               WRITE(LUPRI,'(29X,F14.10)') EMO(KOFFA)
            ENDDO
  112       CONTINUE
         ENDDO
      ENDIF

C     Determine memory split.
C     -----------------------

      CALL MP2_MEMSPL2(LWORK,NVECTOT,LWRK)
      LEFT = LWORK - LWRK

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         XMB = XLWORK*WTOMB
         WRITE(LUPRI,'(12X,A,/,12X,A,I10,A,F7.1,A)')
     &   '(ai|#bj) integrals will be calculated on the fly.',
     &   'Available memory: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(12X,A,1X,I10,A)')
     &   'Maximum memory for integral intermediates: ',LWRK,' words.'
         IF (MP2SAV) THEN
            WRITE(LUPRI,'(12X,A)')
     &      'MP2 amplitudes (and diagonal) will be saved on disk.'
         ENDIF
         WRITE(LUPRI,'(A)') ' '
         WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &   '   #b        First           Last         Mem. Usage',
     &   '        Time   ',
     &   '          (b, sym. b)     (b, sym. b)        (%)    ',
     &   '      (seconds)'
         WRITE(LUPRI,'(5X,A,A)')
     &   '----------------------------------------------------',
     &   '---------------'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Start of batched loop over b.
C     -----------------------------

Casm
C
C     BBEG  = 1
C
C     Restart MP2
C
      IF (RSTMP2) THEN
         BBEG = 0
         DO ISYM = 1,IFSYMB-1
            BBEG = BBEG + NVIR(ISYM)
         END DO
         BBEG = BBEG +  IFVIRB
      ELSE
         BBEG = 1
      END IF
      LMPRST = -1
      LDUM   = .FALSE.
      CALL GPOPEN(LMPRST,'CHOMP2_RST','UNKNOWN',' ','FORMATTED',
     &                  IDUM,LDUM)
      REWIND(LMPRST)
      WRITE(LMPRST,'(/,10X,A,/,10X,A,//)')
     &          'Cholesky MP2 restart information',
     &          '--------------------------------'
      WRITE(LMPRST,'(2X,3(A5,3X),5X,3(A5,3X),11X,A6)')
     &      'FstSy','FstOr','Bfrst','LstSy','LstOrb','Blast','Energy'
      CALL FLSHFO(LMPRST)
C
Casm
      NUMB  = 0
      NPASS = 0
  100 CONTINUE

         BBEG = BBEG + NUMB
         IF (BBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TBBAT = SECOND()

C        Find out how many b's can be treated.
C        -------------------------------------

         CALL GET_BTCH(LVIR,BBEG,NUMB,MEM,LWRK)

C        Check for errors.
C        -----------------

         BEND = BBEG + NUMB - 1
         IF (BEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First B           : ',BBEG,
     &      'Last  B           : ',BEND,
     &      'Number of virtuals: ',NVIRT
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Available memory  : ',LWRK,
     &      'Needed    memory  : ',MEM
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF
         
         IF (NUMB .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for b-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get symmetries of first and last b in batch.
C        --------------------------------------------

         ISBBEG = ISVI(BBEG)
         ISBEND = ISVI(BEND)

C        Complete allocation for integral intermediates.
C        -----------------------------------------------

         KXINT = 1
         KEND1 = KXINT + MEM
         LWRK1 = LWORK - KEND1 + 1

C        Initialize integral intermediates.
C        ----------------------------------

         CALL DZERO(WORK(KXINT),MEM)

C        Set up local index arrays.
C        --------------------------

         CALL IZERO(IOFB1,NSYM)

         IB1 = BBEG
         DO ISYMB = ISBBEG,ISBEND
            IOFB1(ISYMB) = IB1 + NRHFT - IVIR(ISYMB)
            IB1 = IB1 + LVIR(ISYMB)
         ENDDO

         ICOUN2 =  0
         DO ISYMBJ = 1,NSYM

            ICOUN1 = 0
            DO ISYMJ = 1,NSYM
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
               IX1AM(ISYMB,ISYMJ) = ICOUN1
               ICOUN1 = ICOUN1 + LVIR(ISYMB)*NRHF(ISYMJ)
            ENDDO
            NX1AM(ISYMBJ) = ICOUN1

            IX2SQ(ISYMBJ) = ICOUN2
            ICOUN2 = ICOUN2 + NT1AM(ISYMBJ)*NX1AM(ISYMBJ)

         ENDDO
         NX2SQ = ICOUN2

         IF (NX2SQ .NE. MEM) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Error detected in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &      'NX2SQ is calculated to be ',NX2SQ,
     &      'MEM determines alloc. as  ',MEM,
     &      'These numbers must be equal....'
            CALL QUIT('Error in '//SECNAM)
         ENDIF

C        Start loop over Cholesky symmetries.
C        ------------------------------------

         DO ISYCHO = 1,NSYM

            IF (NVECTOT(ISYCHO) .LE. 0) GO TO 999
            IF (NT1AM(ISYCHO)   .LE. 0) GO TO 999
            IF (NX1AM(ISYCHO)   .LE. 0) GO TO 999

C           Set up Cholesky batch in remaining memory.
C           ------------------------------------------

            MINMEM = NT1AM(ISYCHO) + NX1AM(ISYCHO)
            NVEC   = MIN(LWRK1/MINMEM,NVECTOT(ISYCHO))

            IF (NVEC .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory for Cholesky batch in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Memory left for Cholesky batch   : ',LWRK1,
     &         'Minimum required (pref. multiple): ',MINMEM
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            NBATCH = (NVECTOT(ISYCHO) - 1)/NVEC + 1

C           Start batch loop.
C           -----------------

            DO IBATCH = 1,NBATCH

               NUMV = NVEC
               IF (IBATCH .EQ. NBATCH) THEN
                  NUMV = NVECTOT(ISYCHO) - NVEC*(NBATCH - 1)
               ENDIF

               JVEC1 = NVEC*(IBATCH - 1) + 1

               KCHAI = KEND1
               KCHBJ = KCHAI + NT1AM(ISYCHO)*NUMV
               KEND2 = KCHBJ + NX1AM(ISYCHO)*NUMV
               LWRK2 = LWORK - KEND2 + 1

               KLAST  = KEND2 - 1
               MXMEML = MAX(MXMEML,KLAST)

C              Read MO Cholesky vectors.
C              -------------------------

               DTIME = SECOND()
               CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYCHO),NUMV,
     &                         JVEC1,LUCHMO(ISYCHO))
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

C              Get MO Cholesky subblock L(#bj,J).
C              ----------------------------------

               DTIME = SECOND()
               DO IVEC = 1,NUMV
                  DO ISYMJ = 1,NSYM

                     ISYMB = MULD2H(ISYMJ,ISYCHO)
                     IF (LVIR(ISYMB) .LE. 0) GOTO 998

                     DO J = 1,NRHF(ISYMJ)

                        KOFF1 = KCHAI + NT1AM(ISYCHO)*(IVEC - 1)
     &                        + IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                        + IOFB1(ISYMB) - 1
                        KOFF2 = KCHBJ + NX1AM(ISYCHO)*(IVEC - 1)
     &                        + IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1)

                        CALL DCOPY(LVIR(ISYMB),WORK(KOFF1),1,
     &                                         WORK(KOFF2),1)

                     ENDDO

  998                CONTINUE

                  ENDDO
               ENDDO
               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME

C              Calculate (ai|#bj) contribution.
C              --------------------------------

               DTIME = SECOND()

               NTOAI = MAX(NT1AM(ISYCHO),1)
               NTOBJ = MAX(NX1AM(ISYCHO),1)

               KOFFX = KXINT + IX2SQ(ISYCHO)

               CALL DGEMM('N','T',NT1AM(ISYCHO),NX1AM(ISYCHO),NUMV,
     &                    ONE,WORK(KCHAI),NTOAI,WORK(KCHBJ),NTOBJ,
     &                    ONE,WORK(KOFFX),NTOAI)

               DTIME = SECOND() - DTIME
               TIMC  = TIMC     + DTIME

            ENDDO

  999       CONTINUE

         ENDDO

C        Debug: Calculate norm of (ai|bj).
C        ---------------------------------

         IF (LOCDBG) THEN
            XNRM = XNRM + DDOT(NX2SQ,WORK(KXINT),1,WORK(KXINT),1)
         ENDIF

C        Calculate contribution to EMP2.
C        -------------------------------

         DTIME = SECOND()
         DO ISYMBJ = 1,NSYM

            ISYMAI = ISYMBJ

            KT2AM = KEND1
            KENDE = KT2AM + NT1AM(ISYMAI)
            LWRKE = LWORK - KENDE + 1

            KLAST  = KENDE - 1
            MXMEML = MAX(MXMEML,KLAST)

            IF (LWRKE .LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Insufficient memory in ',SECNAM,
     &         ' - Allocation: T2 subblock'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need     : ',KENDE-1,
     &         'Available: ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            DO ISYMJ = 1,NSYM

               ISYMB = MULD2H(ISYMJ,ISYMBJ)
               IF (LVIR(ISYMB) .LE. 0) GOTO 1000

               DO J = 1,NRHF(ISYMJ)

                  KOFFJ = IRHF(ISYMJ) + J

                  DO LB = 1,LVIR(ISYMB)

                     B     = IOFB1(ISYMB) + LB - 1
                     KOFFB = IVIR(ISYMB)  + B

                     LBJ = IX1AM(ISYMB,ISYMJ) + LVIR(ISYMB)*(J - 1)
     &                   + LB
                     NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                   + B

                     DO ISYMI = 1,NSYM

                        ISYMA  = MULD2H(ISYMI,ISYMAI)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = ISYMAJ

                        DO I = 1,NRHF(ISYMI)

                           KOFFI = IRHF(ISYMI) + I
                           LBI   = IX1AM(ISYMB,ISYMI)
     &                           + LVIR(ISYMB)*(I - 1) + LB

                           DO A = 1,NVIR(ISYMA)

                              KOFFA = IVIR(ISYMA) + A
                              NAI   = IT1AM(ISYMA,ISYMI)
     &                              + NVIR(ISYMA)*(I - 1) + A
                              NAJ   = IT1AM(ISYMA,ISYMJ)
     &                              + NVIR(ISYMA)*(J - 1) + A

                              NAIBJ = KXINT + IX2SQ(ISYMBJ)
     &                              + NT1AM(ISYMAI)*(LBJ - 1) + NAI - 1
                              NAJBI = KXINT + IX2SQ(ISYMAJ)
     &                              + NT1AM(ISYMAJ)*(LBI - 1) + NAJ - 1
                              KOFFT = KT2AM + NAI - 1

                              DENOM = EMO(KOFFA) + EMO(KOFFB)
     &                              - EMO(KOFFI) - EMO(KOFFJ)
                              WORK(KOFFT) = -WORK(NAIBJ)/DENOM
                              XAIBJ = TWO*WORK(NAIBJ) - WORK(NAJBI)
                              EMP2  = EMP2 + WORK(KOFFT)*XAIBJ

                           ENDDO

                        ENDDO

                     ENDDO

C                    Debug: calculate (packed) amplitude norm.
C                    -----------------------------------------

                     IF (LOCDBG) THEN
                        TNRM = TNRM + DDOT(NBJ,WORK(KT2AM),1,
     &                                         WORK(KT2AM),1)
                     ENDIF

C                    If requested, save (full square) amplitudes.
C                    Store diagonal element in T2DIA.
C                    --------------------------------------------

                     IF (MP2SAV) THEN
                        DTIM2 = SECOND()
                        CALL WR_MP2TAM(WORK(KT2AM),NT1AM(ISYMBJ),NBJ,
     &                                 LUTAM(ISYMBJ))
                        KOFF1 = KT2AM + NBJ - 1
                        KOFF2 = IOFFD(ISYMBJ) + NBJ
                        T2DIA(KOFF2) = WORK(KOFF1)
                        DTIM2 = SECOND() - DTIM2
                        TIMA  = TIMA     + DTIM2
                     ENDIF

                  ENDDO

               ENDDO

 1000          CONTINUE

            ENDDO

         ENDDO
         DTIME = SECOND() - DTIME
         TIME  = TIME     + DTIME

C        Print information from this b-batch.
C        ------------------------------------

         IF (IPRINT .GE. INFO) THEN
            TBBAT = SECOND() - TBBAT
            XMUSE = (D100*MXMEML)/XLWORK
            LBFST = IOFB1(ISBBEG)
            LBLST = IOFB1(ISBEND) + LVIR(ISBEND) - 1
         WRITE(LUPRI,'(5X,I6,4X,I6,1X,I1,8X,I6,1X,I1,9X,F6.2,8X,F10.2)')
     &      NUMB,LBFST,ISBBEG,LBLST,ISBEND,XMUSE,TBBAT
            CALL FLSHFO(LUPRI)
         ENDIF

C
C        Restart info
C
         WRITE(LMPRST,'(2(3X,I2,3X,I6,2X,I6,7X),F17.10)')
     &         ISBBEG,LBFST,BBEG,ISBEND,LBLST,BEND,EMP2
         CALL FLSHFO(LMPRST)

C        Update memory info.
C        -------------------

         MXMEMT = MAX(MXMEMT,MXMEML)
         MXMEML = 0

C        Update NPASS and go to next b-batch.
C        (Which will exit immediately if no more b's.)
C        ---------------------------------------------

         NPASS = NPASS + 1
         GO TO 100

C     Exit b-batch: calculation done.
C     -------------------------------

  200 CONTINUE

C     If requested, save amplitude diagonal.
C     Close file(s). Release units.
C     --------------------------------------

      IF (MP2SAV) THEN
         DTIME = SECOND()
         DO ISYMAI = 1,NSYM
            KOFFD = IOFFD(ISYMAI) + 1
            CALL WR_MP2TAM(T2DIA(KOFFD),NT1AM(ISYMAI),0,LUTAM(ISYMAI))
         ENDDO
         CALL OP_MP2TAM(IKEEP,LSYM,LUTAM,NSYM)
         DTIME = SECOND() - DTIME
         TIMA  = TIMA     + DTIME
      ENDIF

C     Close files for MO Cholesky vectors.
C     ------------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IKEEP,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Print section.
C     --------------

      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(5X,A,A,/)')
     &   '----------------------------------------------------',
     &   '---------------'
         CALL FLSHFO(LUPRI)
      ENDIF

      IF ((IPRINT.GT.0) .OR. LOCDBG) THEN
         IF (LOCDBG) THEN
            XNRM = DSQRT(XNRM)
            TNRM = DSQRT(TNRM)
            CALL HEADER('Debug info from '//SECNAM,-1)
            WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &      'Norm of (ai|bj) [full square]: ',XNRM
            WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &      'Norm of T2AM    [packed]     : ',TNRM
         ENDIF
         TIMT  = SECOND() - TIMT
         XMUSE = (D100*MXMEMT)/XLWORK
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &   'Memory reserved for integrals        : ',LWRK,
     &   'Memory reserved for Cholesky vectors : ',LEFT,
     &   'Total memory available               : ',LWORK
         WRITE(LUPRI,'(5X,A,4X,F6.2,A)')
     &   'Maximum memory usage                 : ',XMUSE,' %'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors: ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,I10)')
     &   ' - passes through Cholesky file(s)   : ',NPASS
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sort of MO vectors     : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating (ai|#bj)   : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating MP2 energy : ',TIME,' seconds'
         IF (MP2SAV) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used for saving MP2 amplitudes  : ',TIMA,' seconds'
         ENDIF
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used in total                   : ',TIMT,' seconds'
            WRITE(LUPRI,'(5X,A,/)')
     &      '(Debug calculations and printing included.)'
         ELSE
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Time used in total                   : ',TIMT,' seconds'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck get_btch */
      SUBROUTINE GET_BTCH(LVIR,BBEG,NUMB,MEM,LWRK)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Set up batch # 1 for MP2CHO_1.
C
C     LVIR is the local virtual-counter array (i.e. local version
C     og the global NVIR array).
C
#include <implicit.h>
      INTEGER LVIR(8)
      INTEGER BBEG
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CALL IZERO(LVIR,NSYM)

      B    = BBEG - 1
      MEM  = 0
      NUMB = 0
  100 CONTINUE

         B = B + 1
         IF (B .GT. NVIRT) GO TO 200

         ISYMB = ISVI(B)
         MEM   = MEM  + NCKI(ISYMB)
         NUMB  = NUMB + 1
         LVIR(ISYMB) = LVIR(ISYMB) + 1

         IF (MEM .LE. LWRK) THEN
            GO TO 100
         ELSE
            MEM  = MEM  - NCKI(ISYMB)
            NUMB = NUMB - 1
            LVIR(ISYMB) = LVIR(ISYMB) - 1
            GO TO 200
         ENDIF

  200 RETURN
      END
C  /* Deck mp2cho_2 */
      SUBROUTINE MP2CHO_2(T2DIA,EMO,XINT,WORK,LWORK,EMP2)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Calculate the MP2 energy correction EMP2 from Cholesky
C        decomposed two-electron integrals on disk. Calculate
C        the (ai|bj) integrals, stored as full square in XINT,
C        before final contraction to give EMP2.
C
C     Notes:
C        While expectedly faster than MP2CHO_1, this routine requires
C        substantially more memory, as the integrals must be held in
C        core (full square).
C
C        If amplitudes are not to be saved (MP2SAV=.FALSE. from input),
C        array T2DIA may be dummy.
C
C     NB: BE SURE TO UPDATE MP2CHO_2A IF YOU CHANGE ANY ALLOCATIONS IN
C         THIS ROUTINE.
C
#include <implicit.h>
      DIMENSION T2DIA(*), EMO(*), XINT(*), WORK(LWORK)
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <priunit.h>

      INTEGER IOFFD(8), NVECTOT(8)
      INTEGER LUTAM(8), LUCHMO(8), LSYM(8)

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0, TWO = 2.00D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.00D0)

      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IKEEP = 1)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'MP2CHO_2')

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMC = ZERO
      TIME = ZERO
      TIMA = ZERO

C     Initial check of memory.
C     ------------------------

      NEED = 0
      DO ISYCHO = 1,NSYM
         NEED = MAX(NEED,NT1AM(ISYCHO))
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initializations.
C     ----------------

      MXMEMT = 0
      XLWORK = ONE*LWORK

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

C     Set # vectors.
C     --------------

      IF (CHOMO) THEN
         ITYP = 2
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NCHOMO(ISYM)
         ENDDO
      ELSE
         ITYP = 1
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NUMCHO(ISYM)
         ENDDO
      ENDIF

C     Index array for T2DIA.
C     ----------------------

      ICOUNT = 0
      DO ISYMAI = 1,NSYM
         IOFFD(ISYMAI) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYMAI)
      ENDDO

C     Open files for MO Cholesky vectors.
C     -----------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IOPEN,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     If amplitudes are to be saved, open files.
C     ------------------------------------------

      IF (MP2SAV) THEN
         DTIME = SECOND()
         CALL OP_MP2TAM(IOPEN,LSYM,LUTAM,NSYM)
         DTIME = SECOND() - DTIME
         TIMA  = TIMA     + DTIME
      ENDIF

C     Print header.
C     -------------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         XMB = XLWORK*WTOMB
         WRITE(LUPRI,'(15X,A,/,15X,A,I10,A,F7.1,A)')
     &   'Full square (ai|bj) integral matrix held in core.',
     &   'Available memory: ',LWORK,' words (',XMB,' Mb).'
         IF (MP2SAV) THEN
            WRITE(LUPRI,'(15X,A)')
     &      'MP2 amplitudes (and diagonal) will be saved on disk.'
         ENDIF
         WRITE(LUPRI,'(A)') ' '
         CALL FLSHFO(LUPRI)
      ENDIF

C     Initialize integral array.
C     --------------------------

      CALL DZERO(XINT,NT2SQ(1))

C     Debug: initialize amplitude norm.
C     ---------------------------------

      IF (LOCDBG) THEN
         TNRM = ZERO
      ENDIF

C     Print EMO.
C     ----------

      IF (IPRINT .GT. 120) THEN
         CALL AROUND('EMO in '//SECNAM)
         DO ISYM = 1,NSYM
            IF ((NRHF(ISYM).GT.0) .OR. (NVIR(ISYM).GT.0)) THEN
               WRITE(LUPRI,'(//,29X,A,I1)')
     &         'Symmetry Block ',ISYM
            ELSE
               GO TO 112
            ENDIF
            IF (NRHF(ISYM) .EQ. 0) GO TO 111
            CALL HEADER('Occupied Part',-1)
            DO I = 1,NRHF(ISYM)
               KOFFI = IRHF(ISYM) + I
               WRITE(LUPRI,'(29X,F14.10)') EMO(KOFFI)
            ENDDO
  111       IF (NVIR(ISYM) .EQ. 0) GO TO 112
            CALL HEADER('Virtual Part',-1)
            DO A = 1,NVIR(ISYM)
               KOFFA = IVIR(ISYM) + A
               WRITE(LUPRI,'(29X,F14.10)') EMO(KOFFA)
            ENDDO
  112       CONTINUE
         ENDDO
      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(/,10X,A,/,10X,A)')
     &   'Cholesky    First      Last       Mem. Usage     Time   ',
     &   'Symmetry    Vector     Vector        (%)       (seconds)'
         WRITE(LUPRI,'(10X,A)')
     &   '--------------------------------------------------------'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NVECTOT(ISYCHO) .LE. 0) GO TO 999
         IF (NT1AM(ISYCHO)   .LE. 0) GO TO 999

C        Set up batch.
C        -------------

         NVEC = MIN(LWORK/NT1AM(ISYCHO),NVECTOT(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Memory left for Cholesky batch   : ',LWORK,
     &      'Minimum required (pref. multiple): ',NT1AM(ISYCHO)
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NVECTOT(ISYCHO) - 1)/NVEC + 1

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            TIBAT = SECOND()

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NVECTOT(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI = 1
            KBLST = KCHAI + NT1AM(ISYCHO)*NUMV - 1

            MXMEMT = MAX(MXMEMT,KBLST)

C           Read MO Cholesky vectors.
C           -------------------------

            DTIME = SECOND()
            CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYCHO),NUMV,
     &                      JVEC1,LUCHMO(ISYCHO))
            DTIME = SECOND() - DTIME
            TIMR  = TIMR     + DTIME

C           Calculate contribution to (ai|bj).
C           ----------------------------------

            DTIME = SECOND()

            NTOAI = NT1AM(ISYCHO)
            NTOBJ = NTOAI

            KOFF1 = IT2SQ(ISYCHO,ISYCHO) + 1

            CALL DGEMM('N','T',NT1AM(ISYCHO),NT1AM(ISYCHO),NUMV,
     &                 ONE,WORK(KCHAI),NTOAI,WORK(KCHAI),NTOBJ,
     &                 ONE,XINT(KOFF1),NTOAI)

            DTIME = SECOND() - DTIME
            TIMC  = TIMC     + DTIME

C           Print information from this batch.
C           ----------------------------------

            IF (IPRINT .GE. INFO) THEN
               TIBAT = SECOND() - TIBAT
               XMUSE = (D100*KBLST)/XLWORK
               JVECE = JVEC1 + NUMV - 1
               IF (JVEC1 .EQ. 1) THEN
                  WRITE(LUPRI,'(13X,I1,8X,I6,5X,I6,6X,F6.2,5X,F10.2)')
     &            ISYCHO,JVEC1,JVECE,XMUSE,TIBAT
               ELSE
                  WRITE(LUPRI,'(22X,I6,5X,I6,6X,F6.2,5X,F10.2)')
     &            JVEC1,JVECE,XMUSE,TIBAT
               ENDIF
               CALL FLSHFO(LUPRI)
            ENDIF

         ENDDO

  999    CONTINUE

      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         WRITE(LUPRI,'(10X,A)')
     &   '--------------------------------------------------------'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Calculate EMP2 from (ai|bj) integrals in XINT.
C     If requested, save amplitudes.
C     ----------------------------------------------

      DTIME = SECOND()
      EMP2  = ZERO
      DO ISYMBJ = 1,NSYM

         ISYMAI = ISYMBJ

         KT2AM = 1
         KENDE = KT2AM + NT1AM(ISYMAI)
         LWRKE = LWORK - KENDE + 1

         KLAST  = KENDE - 1
         MXMEMT = MAX(MXMEMT,KLAST)

         IF (LWRKE .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: T2 subblock'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KENDE-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         DO ISYMJ = 1,NSYM

            ISYMB = MULD2H(ISYMJ,ISYMBJ)

            DO J = 1,NRHF(ISYMJ)

               KOFFJ = IRHF(ISYMJ) + J

               DO B = 1,NVIR(ISYMB)

                  KOFFB = IVIR(ISYMB) + B
                  NBJ   = IT1AM(ISYMB,ISYMJ)
     &                  + NVIR(ISYMB)*(J - 1) + B

                  DO ISYMI = 1,NSYM

                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = ISYMAJ

                     DO I = 1,NRHF(ISYMI)

                        KOFFI = IRHF(ISYMI) + I
                        NBI   = IT1AM(ISYMB,ISYMI)
     &                        + NVIR(ISYMB)*(I - 1) + B

                        DO A = 1,NVIR(ISYMA)

                           KOFFA = IVIR(ISYMA) + A
                           NAI   = IT1AM(ISYMA,ISYMI)
     &                           + NVIR(ISYMA)*(I - 1) + A
                           NAJ   = IT1AM(ISYMA,ISYMJ)
     &                           + NVIR(ISYMA)*(J - 1) + A

                           NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     &                           + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           NAJBI = IT2SQ(ISYMAJ,ISYMBI)
     &                           + NT1AM(ISYMAJ)*(NBI - 1) + NAJ
                           KOFFT = KT2AM + NAI - 1

                           DENOM = EMO(KOFFA) + EMO(KOFFB)
     &                           - EMO(KOFFI) - EMO(KOFFJ)
                           WORK(KOFFT) = -XINT(NAIBJ)/DENOM
                           XAIBJ = TWO*XINT(NAIBJ) - XINT(NAJBI)
                           EMP2  = EMP2 + WORK(KOFFT)*XAIBJ

                        ENDDO

                     ENDDO

                  ENDDO

C                 Debug: calculate (packed) amplitude norm.
C                 -----------------------------------------

                  IF (LOCDBG) THEN
                     TNRM = TNRM + DDOT(NBJ,WORK(KT2AM),1,WORK(KT2AM),1)
                  ENDIF

C                 If requested, save (full square) amplitudes.
C                 Store diagonal element in T2DIA.
C                 --------------------------------------------

                  IF (MP2SAV) THEN
                     DTIM2 = SECOND()
                     CALL WR_MP2TAM(WORK(KT2AM),NT1AM(ISYMBJ),NBJ,
     &                              LUTAM(ISYMBJ))
                     KOFF1 = KT2AM + NBJ - 1
                     KOFF2 = IOFFD(ISYMBJ) + NBJ
                     T2DIA(KOFF2) = WORK(KOFF1)
                     DTIM2 = SECOND() - DTIM2
                     TIMA  = TIMA     + DTIM2
                  ENDIF

               ENDDO
               
            ENDDO

         ENDDO

      ENDDO
      DTIME = SECOND() - DTIME 
      TIME  = TIME     + DTIME

C     If requested, save amplitude diagonal.
C     Close file(s). Release units.
C     --------------------------------------

      IF (MP2SAV) THEN
         DTIME = SECOND()
         DO ISYMAI = 1,NSYM
            KOFFD = IOFFD(ISYMAI) + 1
            CALL WR_MP2TAM(T2DIA(KOFFD),NT1AM(ISYMAI),0,LUTAM(ISYMAI))
         ENDDO
         CALL OP_MP2TAM(IKEEP,LSYM,LUTAM,NSYM)
         DTIME = SECOND() - DTIME
         TIMA  = TIMA     + DTIME
      ENDIF

C     Close files for MO Cholesky vectors.
C     ------------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IKEEP,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Debug: print integral and amplitude norms.
C     ------------------------------------------

      IF (LOCDBG) THEN
         XNRM = DSQRT(DDOT(NT2SQ(1),XINT,1,XINT,1))
         TNRM = DSQRT(TNRM)
         CALL HEADER('Debug info from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &   'Norm of (ai|bj) [full square]: ',XNRM
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of T2AM    [packed]     : ',TNRM
      ENDIF

C     Print section.
C     --------------

      IF ((IPRINT.GT.0) .OR. LOCDBG) THEN
         TIMT  = SECOND() - TIMT
         XMUSE = (D100*MXMEMT)/XLWORK
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,4X,F6.2,A)')
     &   'Maximum memory usage                 : ',XMUSE,' %'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors: ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating (ai|bj)    : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating MP2 energy : ',TIME,' seconds'
         IF (MP2SAV) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used for saving MP2 amplitudes  : ',TIMA,' seconds'
         ENDIF
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used in total                   : ',TIMT,' seconds'
            WRITE(LUPRI,'(5X,A,/)')
     &      '(Debug calculations and printing included.)'
         ELSE
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Time used in total                   : ',TIMT,' seconds'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck mp2cho_2a */
      SUBROUTINE MP2CHO_2A(LWORK,MINVEC,POSSBL)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Run through the structure of MP2CHO_2 *without* doing any
C        calculations and check if sufficient memory is available.
C        Array MINVEC gives the minimum number of Cholesky vectors
C        (in each symmetry) that it must be possible to hold in core.
C
C     NB: BE SURE TO UPDATE THIS ROUTINE IF YOU CHANGE ANY ALLOCATIONS
C         IN THE "PARENT" ROUTINE MP2CHO_2.
C
#include <implicit.h>
      INTEGER MINVEC(8)
      LOGICAL POSSBL
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <priunit.h>

      INTEGER NVECTOT(8)

C     Set the default.
C     ----------------

      POSSBL = .TRUE.

C     MP2SAV-specific test.
C     ---------------------

      IF (MP2SAV) THEN
         MAXT1 = 0
         DO ISYMAI = 1,NSYM
            MAXT1 = MAX(MAXT1,NT1AM(ISYMAI))
         ENDDO
         IF (LWORK .LT. MAXT1) THEN
            POSSBL = .FALSE.
            RETURN
         ENDIF
      ENDIF

C     Run through "skeleton" of MP2CHO_2.
C     -----------------------------------

      DO ISYCHO = 1,NSYM

         IF (CHOMO) THEN
            NVECTOT(ISYCHO) = NCHOMO(ISYCHO)
         ELSE
            NVECTOT(ISYCHO) = NUMCHO(ISYCHO)
         ENDIF

         IF (NVECTOT(ISYCHO) .LE. 0) GO TO 999
         IF (NT1AM(ISYCHO)   .LE. 0) GO TO 999

C        Set up batch.
C        -------------

         NVEC = MIN(LWORK/NT1AM(ISYCHO),NVECTOT(ISYCHO))

         IF (NVEC .LT. MINVEC(ISYCHO)) THEN
            POSSBL = .FALSE.
            RETURN
         ENDIF

  999    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_choiajb */
      SUBROUTINE CC_CHOIAJB(XAIBJ,T1AM,WORK,LWORK)
C
C     Thomas Bondo Pedersen, August 2001.
C     - revised, tbp, May 2002.
C
C     Purpose:
C        Calculate the (ia|jb) integrals [sorted as ai,bj]
C        from Cholesky vectors. Lambda particle and hole
C        matrices are calculated in LAMMAT using T1AM, and
C        they are *not* assumed identical, though both are
C        assumed tot. sym.
C
C        For NOCCD, provided LAMMAT has been appropriately
C        modified, T1AM is dummy.
C
#include <implicit.h>
      DIMENSION XAIBJ(*), T1AM(*)
      DIMENSION WORK(LWORK)
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOIAJB')

C     Allocation for Lambda matrices.
C     -------------------------------

      KLAMDP = 1
      KLAMDH = KLAMDP + NGLMDT(1)
      KEND1  = KLAMDH + NGLMDT(1)
      LWRK1  = LWORK  - KEND1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - Allocation: Lambda'
         WRITE(LUPRI,'(5X,A,I10)')   'Need     : ',KEND1
         WRITE(LUPRI,'(5X,A,I10,/)') 'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate Lambda matrices.
C     --------------------------

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND1),LWRK1)

C     Initialize integral array XAIBJ.
C     --------------------------------

      CALL DZERO(XAIBJ,NT2AM(1))

C     Calculate (ia|jb) integrals.
C     ----------------------------

      CALL CC_CHOIAJB1(XAIBJ,WORK(KLAMDP),WORK(KLAMDH),
     &                 WORK(KEND1),LWRK1)

      RETURN
      END
C  /* Deck cc_choiajb1 */
      SUBROUTINE CC_CHOIAJB1(XAIBJ,XLAMDP,XLAMDH,WORK,LWORK)
C
C     Thomas Bondo Pedersen, August 2001.
C     - revised, tbp, May 2002.
C
C     Purpose: Calculate (ia|jb) integrals from Cholesky vectors
C              reading one vector at a time but keeping as many
C              MO (ia) vectors in memory as possible.
C
C     The Lambda matrices are assumed total symmetric; the integral
C     array XAIBJ must be initialized.
C
C     XAIBJ is packed and returned as ai<=bj.
C
#include <implicit.h>
      DIMENSION XAIBJ(*), XLAMDP(*), XLAMDH(*)
      DIMENSION WORK(LWORK)
      INTEGER AI, BJ, AIBJ
#include <iratdef.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOIAJB1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)
      PARAMETER (IOPTR = 2)

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMM = ZERO
      TIMS = ZERO
      TIMC = ZERO

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 999
         IF (N2BST(ISYCHO)  .LE. 0) GO TO 999
         IF (NT1AM(ISYCHO)  .LE. 0) GO TO 999

C        Allocate space for AO Cholesky vectors and
C        for integral intermediate.
C        ------------------------------------------

         MXAIB = 0
         MXB   = 0
         DO ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYCHO)
            NTST  = NT1AM(ISYCHO)*NVIR(ISYMB)
            MXAIB = MAX(MXAIB,NTST)
            MXB   = MAX(MXB,NVIR(ISYMB))
         ENDDO

         LREAD = MEMRD(1,ISYCHO,IOPTR)
         LSCR  = MAX(LREAD,NT1AO(ISYCHO))

         KCHOL = KEND0
         KSCR  = KCHOL + N2BST(ISYCHO)
         KXAIB = KSCR  + LSCR
         KEND1 = KXAIB + MXAIB
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &      'Need (more than)           : ',KEND1-1,
     &      'Available                  : ',LWORK
            MEMINC = KEND1 - LWORK + NT1AM(ISYCHO) + MXB - 1
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Increase memory by at least: ',MEMINC
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYCHO) + MXB
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory needed: ',MINMEM,
     &      'Available for batch  : ',LWRK1,
     &      'Available in total   : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI = KEND1
            KSCRB = KCHAI + NT1AM(ISYCHO)*NUMV

C           Read and transform Cholesky vectors in this batch.
C           --------------------------------------------------

            DO IVEC = 1,NUMV

               JVEC  = JVEC1 + IVEC - 1
               DTIME = SECOND()
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KSCR),LSCR)
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

               DTIME = SECOND()
               DO ISYMD = 1,NSYM

                  ISYMG = MULD2H(ISYMD,ISYCHO)
                  ISYMI = ISYMG
                  ISYMA = ISYMD

                  KOFFC = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFFP = IGLMRH(ISYMG,ISYMI) + 1
                  KOFFH = IGLMVI(ISYMD,ISYMA) + 1
                  KOFFM = KCHAI + NT1AM(ISYCHO)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  NTOTG = MAX(NBAS(ISYMG),1)
                  NTOTD = MAX(NBAS(ISYMD),1)
                  NTOTA = MAX(NVIR(ISYMA),1)

                  CALL DGEMM('T','N',
     &                       NBAS(ISYMD),NRHF(ISYMI),NBAS(ISYMG),
     &                       ONE,WORK(KOFFC),NTOTG,XLAMDP(KOFFP),NTOTG,
     &                       ZERO,WORK(KSCR),NTOTD)

                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     &                       ONE,XLAMDH(KOFFH),NTOTD,WORK(KSCR),NTOTD,
     &                       ZERO,WORK(KOFFM),NTOTA)

               ENDDO
               DTIME = SECOND() - DTIME
               TIMM  = TIMM     + DTIME

            ENDDO

C           Calculate (ia|jb) integral in loop over j.
C           ------------------------------------------

            DO ISYMJ = 1,NSYM

               ISYMB = MULD2H(ISYMJ,ISYCHO)

               IF (NVIR(ISYMB) .LE. 0) GO TO 998

               DO J = 1,NRHF(ISYMJ)

C                 Copy L(bj,J) --> L(b,J;j).
C                 --------------------------

                  KOFFB = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)

                  DTIME = SECOND()
                  DO IVEC = 1,NUMV
                     KOFF1 = KCHAI + NT1AM(ISYCHO)*(IVEC - 1) + KOFFB
                     KOFF2 = KSCRB + NVIR(ISYMB)*(IVEC - 1)
                     CALL DCOPY(NVIR(ISYMB),WORK(KOFF1),1,WORK(KOFF2),1)
                  ENDDO
                  DTIME = SECOND() - DTIME
                  TIMS  = TIMS     + DTIME

C                 Calculate intermediate result.
C                 ------------------------------

                  DTIME = SECOND()

                  CALL DGEMM('N','T',NT1AM(ISYCHO),NVIR(ISYMB),NUMV,
     &                       ONE,WORK(KCHAI),NT1AM(ISYCHO),
     &                       WORK(KSCRB),NVIR(ISYMB),
     &                       ZERO,WORK(KXAIB),NT1AM(ISYCHO))

C                 Add into packed integral array.
C                 -------------------------------

                  DO B = 1,NVIR(ISYMB)

                     BJ = KOFFB + B

                     DO AI = 1,BJ

                        KOFF = KXAIB + NT1AM(ISYCHO)*(B - 1) + AI - 1
                        AIBJ = IT2AM(ISYCHO,ISYCHO) + INDEX(AI,BJ)

                        XAIBJ(AIBJ) = XAIBJ(AIBJ) + WORK(KOFF)

                     ENDDO

                  ENDDO

                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME

               ENDDO

  998          CONTINUE

            ENDDO

         ENDDO

  999    CONTINUE

      ENDDO

C     Debug: print integral norm.
C     ---------------------------

      IF (LOCDBG) THEN
         XNRM = ZERO
         DO ISYM = 1,NSYM
            DO BJ = 1,NT1AM(ISYM)
               DO AI = 1,NT1AM(ISYM)
                  AIBJ = IT2AM(ISYM,ISYM) + INDEX(AI,BJ)
                  XNRM = XNRM + XAIBJ(AIBJ)*XAIBJ(AIBJ)
               ENDDO
            ENDDO
         ENDDO
         XNRM = DSQRT(XNRM)
         PNRM = DSQRT(DDOT(NT2AM(1),XAIBJ,1,XAIBJ,1))
         WRITE(LUPRI,'(/,5X,A,1P,D30.15)')
     &   'Norm of (ai|bj) [packed     ]: ',PNRM
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of (ai|bj) [full square]: ',XNRM
      ENDIF

C     Print section.
C     --------------

      IF ((IPRINT.GT.0) .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         CALL HEADER('Timing of '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors: ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for MO transformation      : ',TIMM,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sorting MO Cholesky    : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating (ai|bj)    : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used in total                   : ',TIMT,' seconds'
            WRITE(LUPRI,'(5X,A,/)')
     &      '(Debug calculations and printing included.)'
         ELSE
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Time used in total                   : ',TIMT,' seconds'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck diffiajb */
      SUBROUTINE DIFFIAJB(XAIBJ,T1AM,WORK,LWORK,CHOTIM,LSTOP)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Compare (ia|jb) integrals in XAIBJ to those obtained
C        from the integral-direct CCSD_IAJB routine.
C
C     Notes:
C        - The integral array XAIBJ is returned unchanged.
C        - CHOTIM is the time spend (outside) for calculating
C          the Cholesky integrals; for printing purposes. Printing
C          of this time may be suppressed by supplying a negative
C          CHOTIM.
C        - Logical variable LSTOP acts as input:
C          If, on input, LSTOP=.TRUE. then the execution is terminated
C          upon completion of this routine.
C
#include <implicit.h>
      DIMENSION XAIBJ(*), T1AM(*), WORK(LWORK)
      LOGICAL LSTOP, CHOINSV
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'DIFFIAJB')

      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0)

C     Allocation.
C     -----------

      KAIBJ = 1
      KEND1 = KAIBJ + NT2AM(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory for test in ',SECNAM,
     &   ' - Allocation: integrals'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK,
     &   ' - Test skipped, program continues....'
         GO TO 100
      ENDIF

C     Calculate integrals.
C     --------------------

      CHOINSV = CHOINT
      CHOINT  = .FALSE.

      TIMI = SECOND()
      CALL CCSD_IAJB(WORK(KAIBJ),T1AM,WORK(KEND1),LWRK1)
      TIMI = SECOND() - TIMI

      CHOINT = CHOINSV

C     Calculate difference and find max.
C     ----------------------------------

      CALL DAXPY(NT2AM(1),XMONE,XAIBJ,1,WORK(KAIBJ),1)
      CALL NOCC_T2LARG(WORK(KAIBJ),1,NAI,ISYMAI,NBJ,ISYMBJ,
     &                 NAIBJ,VAL)

C     Print.
C     ------

      DNRM = DSQRT(DDOT(NT2AM(1),WORK(KAIBJ),1,WORK(KAIBJ),1))
      XNRM = DSQRT(DDOT(NT2AM(1),XAIBJ,1,XAIBJ,1))
      IF (DNRM .GT. ZERO) THEN
         CTRB = 100.00D0*DABS(VAL)/DNRM
      ENDIF

      CALL HEADER(SECNAM//': Comparing (ia|jb) integrals',-1)
      IF (CHOTIM .GT. ZERO) THEN
         WRITE(LUPRI,'(5X,A,F15.2,A)')
     &   'Time for Cholesky    calculation: ',CHOTIM,' seconds'
      ENDIF
      WRITE(LUPRI,'(5X,A,F15.2,A)')
     & 'Time for int. direct calculation: ',TIMI,' seconds'
      WRITE(LUPRI,'(5X,A,1P,D15.6)')
     & 'Difference norm                 : ',DNRM
      WRITE(LUPRI,'(5X,A,1P,D15.6)')
     & '- norm of Cholesky integrals    : ',XNRM
      WRITE(LUPRI,'(5X,A,1P,D15.6)')
     & 'Largest element of difference   : ',VAL
      WRITE(LUPRI,'(5X,A,1P,D15.6)')
     & '- value of Cholesky integral    : ',XAIBJ(NAIBJ)
      IF (DNRM .GT. ZERO) THEN
         WRITE(LUPRI,'(5X,A,F15.2,A)')
     &    'Contribution to difference norm : ',CTRB,' %'
      ELSE
         WRITE(LUPRI,'(5X,A)')
     &    'Contribution to difference norm :  - '
      ENDIF
      WRITE(LUPRI,'(5X,A,5X,I10)')
     & 'Location; AIBJ                  : ',NAIBJ
      WRITE(LUPRI,'(5X,A,14X,I1)')
     & 'Location; symmetry block        : ',ISYMAI
      WRITE(LUPRI,'(5X,A,I7,1X,I7,/)')
     & 'Location; AI,BJ                 : ',NAI,NBJ

C     STOP, if requested.
C     -------------------

  100 IF (LSTOP) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   '- ',SECNAM,' stops execution on request...'
         CALL QUIT(' ***** Execution Killed on Request ***** ')
      ENDIF

      RETURN
      END
C  /* Deck cho_aibjd */
      SUBROUTINE CHO_AIBJD(DIAG,WORK,LWORK,ISYMAI,LUCHMO)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Calculate (ai|ai) integral diagonal of symmetry ISYMAI.
C
C     Notes:
C        MO Cholesky vectors, L(ai), are assumed available on
C        disk (file assumed open with unit LUCHMO).
C
#include <implicit.h>
      DIMENSION DIAG(*), WORK(LWORK) 
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <priunit.h>

      INTEGER AI

      PARAMETER (ZERO = 0.00D0)
      PARAMETER (INFO = 20)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_AIBJD')

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMC = ZERO

C     If nothing to do, return.
C     -------------------------

      IF ((NT1AM(ISYMAI).LE.0) .OR. (NUMCHO(ISYMAI).LE.0)) RETURN

C     Initialize diagonal.
C     --------------------

      CALL DZERO(DIAG,NT1AM(ISYMAI))

C     Set up batch.
C     -------------

      ISYCHO = ISYMAI
      MINMEM = NT1AM(ISYMAI)

      NVEC = MIN(LWORK/MINMEM,NUMCHO(ISYCHO))
      IF (NVEC .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory for batch in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',MINMEM,
     &   'Available memory       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C     Start batch loop.
C     -----------------

      DO IBATCH = 1,NBATCH

         NUMV = NVEC
         IF (IBATCH .EQ. NBATCH) THEN
            NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
         ENDIF

         JVEC1 = NVEC*(IBATCH - 1) + 1

         KCHAI = 1

C        Read MO vectors L(ai,J).
C        ------------------------

         DTIME = SECOND()
         CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYMAI),NUMV,JVEC1,LUCHMO)
         DTIME = SECOND() - DTIME
         TIMR  = TIMR     + DTIME

C        Calculate diagonal.
C        -------------------

         DTIME = SECOND()
         DO AI = 1,NT1AM(ISYMAI)
            KOFFC = KCHAI + AI - 1
            DIAG(AI) = DIAG(AI)
     &               + DDOT(NUMV,WORK(KOFFC),NT1AM(ISYMAI),
     &                           WORK(KOFFC),NT1AM(ISYMAI))
         ENDDO
         DTIME = SECOND() - DTIME
         TIMC  = TIMC     + DTIME

      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMT = SECOND() - TIMT
         XNRM = DSQRT(DDOT(NT1AM(ISYMAI),DIAG,1,DIAG,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   '(ai|ai) diagonal calculated, symmetry block: ',ISYMAI
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated diagonal: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for I/O of MO Cholesky vectors : ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for calculating diagonal       : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '--------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                          : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cho_aibjd_i */
      SUBROUTINE CHO_AIBJD_I(DIAG,CHAI,ISYMAI)
C
C     Thomas Bondo Pedersen, September 2002.
C
C     Purpose:
C        Calculate (ai|ai) integral diagonal of symmetry ISYMAI.
C
C     Notes:
C        MO Cholesky vectors, L(ai), are assumed available in
C        CHAI
C
#include <implicit.h>
      DIMENSION DIAG(*), CHAI(*)
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <priunit.h>

      INTEGER AI

      PARAMETER (INFO = 20)

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_AIBJD_I')

C     Initialize timing.
C     ------------------

      TIMT = SECOND()

C     If nothing to do, return.
C     -------------------------

      IF ((NT1AM(ISYMAI).LE.0) .OR. (NUMCHO(ISYMAI).LE.0)) RETURN

C     Initialize diagonal.
C     --------------------

      CALL DZERO(DIAG,NT1AM(ISYMAI))

C     Calculate diagonal.
C     -------------------

      DO AI = 1,NT1AM(ISYMAI)
         DIAG(AI) = DIAG(AI)
     &            + DDOT(NUMCHO(ISYMAI),CHAI(AI),NT1AM(ISYMAI),
     &                                  CHAI(AI),NT1AM(ISYMAI))
      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMT = SECOND() - TIMT
         XNRM = DSQRT(DDOT(NT1AM(ISYMAI),DIAG,1,DIAG,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   '(ai|ai) diagonal calculated, symmetry block: ',ISYMAI
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated diagonal: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time for calculating diagonal       : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cho_aibj */
      SUBROUTINE CHO_AIBJ(XINT,WORK,LWORK,LISTBJ,NUMBJ,ISYMBJ,LUCHMO)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Calculate (ai|bj) integrals for a list, LISTBJ, of NUMBJ column
C        (bj) indices of symmetry ISYMBJ.
C
C     Notes:
C        MO Cholesky vectors, L(ai), are assumed available on disk on
C        opened unit LUCHMO.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK) 
      INTEGER LISTBJ(*)
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <priunit.h>

      INTEGER BJ

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)
      PARAMETER (INFO = 20)

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_AIBJ')

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO

C     Return if nothing to do.
C     ------------------------

      IF ((NUMBJ.LE.0) .OR. (NT1AM(ISYMBJ).LE.0) .OR.
     &    (NUMCHO(ISYMBJ).LE.0)) RETURN

C     Test NUMBJ.
C     -----------

      IF (NUMBJ .GT. NT1AM(ISYMBJ)) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Error detected in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,A,I1,A,/,5X,A,I10,/)')
     &   'Number of bj requested: ',NUMBJ,' (sym. ',ISYMBJ,')',
     &   'Maximum value possible: ',NT1AM(ISYMBJ)
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Initialize integral array.
C     --------------------------

      CALL DZERO(XINT,NT1AM(ISYMBJ)*NUMBJ)

C     Set up batch.
C     -------------

      MINMEM = NT1AM(ISYMBJ) + NUMBJ

      NVEC = MIN(LWORK/MINMEM,NUMCHO(ISYMBJ))
      IF (NVEC .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory for batch in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',MINMEM,
     &   'Available memory       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      NBATCH = (NUMCHO(ISYMBJ) - 1)/NVEC + 1

C     Start batch loop.
C     -----------------

      DO IBATCH = 1,NBATCH

         NUMV = NVEC
         IF (IBATCH .EQ. NBATCH) THEN
            NUMV = NUMCHO(ISYMBJ) - NVEC*(NBATCH - 1)
         ENDIF

         JVEC1 = NVEC*(IBATCH - 1) + 1

         KCHAI = 1
         KCHBJ = KCHAI + NT1AM(ISYMBJ)*NUMV

C        Read MO Cholesky vectors L(ai,J).
C        ---------------------------------

         DTIME = SECOND()
         CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYMBJ),NUMV,JVEC1,LUCHMO)
         DTIME = SECOND() - DTIME
         TIMR  = TIMR     + DTIME

C        Calculate contribution.
C        -----------------------

         DTIME = SECOND()
         DO BJ = 1,NUMBJ
            KOFF1 = KCHAI + LISTBJ(BJ) - 1
            KOFF2 = KCHBJ + NUMV*(BJ - 1)
            CALL DCOPY(NUMV,WORK(KOFF1),NT1AM(ISYMBJ),WORK(KOFF2),1)
         ENDDO
         DTIME = SECOND() - DTIME
         TIMS  = TIMS     + DTIME

         DTIME = SECOND()
         CALL DGEMM('N','N',NT1AM(ISYMBJ),NUMBJ,NUMV,
     &              ONE,WORK(KCHAI),NT1AM(ISYMBJ),WORK(KCHBJ),NUMV,
     &              ONE,XINT,NT1AM(ISYMBJ))
         DTIME = SECOND() - DTIME
         TIMC  = TIMC     + DTIME

      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMT = SECOND() - TIMT
         XNRM = DSQRT(DDOT(NT1AM(ISYMBJ)*NUMBJ,XINT,1,XINT,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   'Integral symmetry block number ',ISYMBJ
         WRITE(LUPRI,'(5X,A,I8,A)')
     &   '(ai|bj) integrals calculated for the following ',NUMBJ,
     &   ' columns:'
         WRITE(LUPRI,'(10I8)') (LISTBJ(I),I=1,NUMBJ)
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated integral array: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for I/O of MO Cholesky vectors : ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for sorting MO Cholesky vectors: ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for calculating (ai|#[bj]) int.: ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '--------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                          : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cho_aibj_i */
      SUBROUTINE CHO_AIBJ_I(XINT,CHAI,CHBJ,LISTBJ,NUMBJ,ISYMBJ)
C
C     Thomas Bondo Pedersen, September 2002.
C
C     Purpose: Calculate (ai|#[bj]); as CHO_AIBJ except vectors
C              assumed in core.
C
#include <implicit.h>
      DIMENSION XINT(*), CHAI(*), CHBJ(*)
      INTEGER LISTBJ(*)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_AIBJ_I')

      PARAMETER (INFO = 20)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

C     Timing.
C     -------

      TIMT = SECOND()

      IF (LOCDBG) THEN
         LEN  = NT1AM(ISYMBJ)*NUMCHO(ISYMBJ)
         XNAI = DDOT(LEN,CHAI,1,CHAI,1)
         WRITE(LUPRI,'(A,A,1P,D15.7)')
     &   SECNAM,': Norm of CHAI before copy: ',XNAI
      ENDIF

C     Copy sub-block L(J,#[bj]).
C     --------------------------

      TIMS = SECOND()
      CALL CC_DCMCPAI(CHAI,CHBJ,LISTBJ,NUMBJ,ISYMBJ)
      TIMS = SECOND() - TIMS

      IF (LOCDBG) THEN
         LEN  = NT1AM(ISYMBJ)*NUMCHO(ISYMBJ)
         XNAI = DDOT(LEN,CHAI,1,CHAI,1)
         WRITE(LUPRI,'(A,A,1P,D15.7)')
     &   SECNAM,': Norm of CHAI after  copy: ',XNAI
         LEN  = NUMCHO(ISYMBJ)*NUMBJ
         XNBJ = DDOT(LEN,CHBJ,1,CHBJ,1)
         WRITE(LUPRI,'(A,A,1P,D15.7)')
     &   SECNAM,': Norm of CHBJ after  copy: ',XNBJ
      ENDIF

C     Calculate.
C     ----------

      TIMC   = SECOND()
      ISYMAI = ISYMBJ
      NTOAI  = MAX(NT1AM(ISYMAI),1)
      NTVEC  = MAX(NUMCHO(ISYMAI),1)
      CALL DGEMM('N','N',NT1AM(ISYMAI),NUMBJ,NUMCHO(ISYMAI),
     &           ONE,CHAI,NTOAI,CHBJ,NTVEC,
     &           ZERO,XINT,NTOAI)
      TIMC   = SECOND() - TIMC

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         TIMT = SECOND() - TIMT
         XNRM = DSQRT(DDOT(NT1AM(ISYMBJ)*NUMBJ,XINT,1,XINT,1))
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I1)')
     &   'Integral symmetry block number ',ISYMBJ
         WRITE(LUPRI,'(5X,A,I8,A)')
     &   '(ai|bj) integrals calculated for the following ',NUMBJ,
     &   ' columns:'
         WRITE(LUPRI,'(10I8)') (LISTBJ(I),I=1,NUMBJ)
         WRITE(LUPRI,'(5X,A,1P,D30.15,/)')
     &   'Norm of calculated integral array: ',XNRM
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for sorting MO Cholesky vectors: ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time for calculating (ai|#[bj]) int.: ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '--------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                          : ',TIMT,' seconds'
      ENDIF      

      RETURN
      END
C  /* Deck cho_trfai */
      SUBROUTINE CHO_TRFAI(CMO,WORK,LWORK)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Transform AO Cholesky vectors to MO (ai) basis
C        and write the resulting vectors to disk.
C
C     Control flag SKIPTR is set true at the end of this routine
C     to avoid future (re-)calculation.
C
#include <implicit.h>
      DIMENSION CMO(*), WORK(LWORK)
#include <iratdef.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CHO_TRFAI')

      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

      PARAMETER (IOPTR = 2, ITYP = 1, IOPEN = -1, IKEEP = 1, INFO = 3)

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMW = ZERO
      TIMM = ZERO

C     Skip if requested.
C     ------------------

      IF (SKIPTR) GOTO 1001

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Cholesky symmetry loop.
C     -----------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (N2BST(ISYCHO)  .LE. 0) GO TO 1000
         IF (NT1AM(ISYCHO)  .LE. 0) GO TO 1000

C        Open file for this symmetry.
C        ----------------------------

         CALL CHO_MOP(IOPEN,ITYP,ISYCHO,LUCHMO,1,1)

C        Allocation.
C        -----------

         ISYMAI = ISYCHO
         LREAD  = MEMRD(1,ISYCHO,IOPTR)
         LREAD  = MAX(LREAD,NT1AO(ISYCHO))

         KCHOL = KEND0
         KREAD = KCHOL + N2BST(ISYCHO)
         KEND1 = KREAD + LREAD
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - allocation: AO Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NT1AM(ISYMAI)

         NVEC = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))
         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required: ',MINMEM,
     &      'Available for batch    : ',LWRK1,
     &      'Available total memory : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JVEC1 = NVEC*(IBATCH - 1) + 1

            KCHAI = KEND1

C           Read AO vectors and transform to MO basis.
C           ------------------------------------------

            DO JVEC = 1,NUMV

               DTIME = SECOND()
               JCHOL = JVEC1 + JVEC - 1
               CALL CHO_READN(WORK(KCHOL),JCHOL,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KREAD),LREAD)
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

               DTIME = SECOND()
               DO ISYMD = 1,NSYM

                  ISYMG = MULD2H(ISYMD,ISYCHO)
                  ISYMA = ISYMG
                  ISYMI = ISYMD

                  NG    = NBAS(ISYMG)
                  NTOTG = MAX(NBAS(ISYMG),1)
                  ND    = NBAS(ISYMD)
                  NTOTD = MAX(NBAS(ISYMD),1)
                  NI    = NRHF(ISYMI)
                  NA    = NVIR(ISYMA)
                  NTOTA = MAX(NVIR(ISYMA),1)

                  KOFF1 = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFF2 = IGLMRH(ISYMD,ISYMI) + 1
                  KOFF3 = IGLMVI(ISYMG,ISYMA) + 1
                  KOFF4 = KCHAI + NT1AM(ISYCHO)*(JVEC - 1)
     &                  + IT1AM(ISYMA,ISYMI)

                  CALL DGEMM('N','N',NG,NI,ND,
     &                       ONE,WORK(KOFF1),NTOTG,CMO(KOFF2),NTOTD,
     &                       ZERO,WORK(KREAD),NTOTG)

                  CALL DGEMM('T','N',NA,NI,NG,
     &                       ONE,CMO(KOFF3),NTOTG,WORK(KREAD),NTOTG,
     &                       ZERO,WORK(KOFF4),NTOTA)

               ENDDO
               DTIME = SECOND() - DTIME
               TIMM  = TIMM     + DTIME

            ENDDO

C           Write MO Cholesky vectors to disk.
C           ----------------------------------

            DTIME = SECOND()
            CALL CHO_MOWRITE(WORK(KCHAI),NT1AM(ISYMAI),NUMV,JVEC1,
     &                       LUCHMO)
            DTIME = SECOND() - DTIME
            TIMW  = TIMW     + DTIME

         ENDDO

C        Close file for this symmetry.
C        -----------------------------

         CALL CHO_MOP(IKEEP,ITYP,ISYCHO,LUCHMO,1,1)

 1000    CONTINUE

      ENDDO

C     Print.
C     ------

 1001 IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         IF (SKIPTR) THEN
            WRITE(LUPRI,'(10X,A,/)')
     &      'Transformation of AO Cholesky vectors skipped.'
         ELSE
            TIMT = SECOND() - TIMT
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time for I/O of AO Cholesky vectors : ',TIMR,' seconds'
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time for I/O of MO Cholesky vectors : ',TIMW,' seconds'
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time for Cholesky MO transformation : ',TIMM,' seconds'
            WRITE(LUPRI,'(5X,A)')
     &      '--------------------------------------------------------'
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Total time                          : ',TIMT,' seconds'
         ENDIF
      ENDIF

C     MO transformation completed:
C     set SKIPTR flag to avoid future recalculations.
C     -----------------------------------------------

      SKIPTR = .TRUE.

      RETURN
      END
C  /* Deck cc_decomp2 */
      SUBROUTINE CC_DECOMP2(WORK,LWORK)
C
C     Thomas Bondo Pedersen, August 2002.
C
C     Purpose: Cholesky decomposition of (ia|jb) integral matrix
C              constructed from AO Cholesky vectors.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <iratdef.h>
#include <chomp2.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_DECOMP2')

      PARAMETER (ZERO = 0.0D0)

C     Set up info for decomposition.
C     ------------------------------

      IF (THRMP2 .LT. ZERO) THRMP2 = THRCOM
      IF (SPAMP2 .LT. ZERO) SPAMP2 = SPAN

      IMAT = 1

      KDIANL = 1
      KCOLUM = KDIANL + 3*NSYM
      KINDIA = KCOLUM + (NSYM - 1)/IRAT   + 1
      KEND1  = KINDIA + (MXDECM - 1)/IRAT + 1
      LWRK1  = LWORK  - KEND1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(6,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CC_DECMO(WORK(KDIANL),NCHOMO,WORK(KCOLUM),WORK(KINDIA),
     &              THRMP2,SPAMP2,
     &              FMP2S,SKIPCH,MXDECM,NCHORD,NT1AM,IPRINT,IMAT,
     &              NSYM,THZMP2,DUMMY,WORK(KEND1),LWRK1)

      RETURN
      END
C  /* Deck get_mp2alg */
      SUBROUTINE GET_MP2ALG(LWORK,LALG)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: Determine MP2 algorithm based on work space availability.
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccorb.h>
#include <dccsdsym.h>
#include <priunit.h>

      INTEGER VFRAC

      INTEGER NVECTOT(8), MINVEC(8)

      PARAMETER (NALG = 3)
      LOGICAL POSSBL(NALG)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'GET_MP2ALG')

      PARAMETER (MINDEF = 5, LFRAC = 10, VFRAC = 100)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Test LWORK.
C     -----------

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(//,5X,A,I10,/)')
     &   'Memory available, LWORK = ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      XWORK = ONE*LWORK

C     Initialize.
C     -----------

      DO IALG = 1,NALG
         POSSBL(IALG) = .FALSE.
      ENDDO

C     LALG=2: Full-square (ia|jb)
C           + 1/LFRAC of the vectors must fit in core.
C     ------------------------------------------------

      IF (XT2SQ(1) .GT. XWORK) THEN

         POSSBL(2) = .FALSE.

      ELSE

         LWRK1 = LWORK - NT2SQ(1)
         DO ISYCHO = 1,NSYM
            IF (CHOMO) THEN
               NVECTOT(ISYCHO) = NCHOMO(ISYCHO)
            ELSE
               NVECTOT(ISYCHO) = NUMCHO(ISYCHO)
            ENDIF
            IF (NVECTOT(ISYCHO) .GT. LFRAC) THEN
               MINVEC(ISYCHO) = NVECTOT(ISYCHO)/LFRAC
            ELSE
               MINVEC(ISYCHO) = MIN(NVECTOT(ISYCHO),MINDEF)
            ENDIF
         ENDDO
         CALL MP2CHO_2A(LWRK1,MINVEC,POSSBL(2))

      ENDIF

C     LALG=1: (ia|j#b) 
C           + L(ia,#J) + L(j#b,#J) must fit in core..
C     1/LFRAC Cholesky vectors must fit in core, and
C     1/VFRAC virtuals must fit in core.
C     -----------------------------------------------

      IF (.NOT. POSSBL(2)) THEN

         XMXCKI = ZERO
         ISYMB  = 0
         DO ISYM = 1,NSYM
            IF ((XCKI(ISYM).GT.XMXCKI) .AND. (NVIR(ISYM).GT.0)) THEN
               XMXCKI = XCKI(ISYM)
               ISYMB  = ISYM
            ENDIF
         ENDDO
         IF (NVIR(ISYMB) .GT. VFRAC) THEN
            NUMB = NVIR(ISYMB)/VFRAC
         ELSE
            NUMB = MIN(NVIR(ISYMB),MINDEF)
         ENDIF
         XMXCKI = XMXCKI*NUMB

         IF (XMXCKI .LE. ZERO) THEN

            POSSBL(1) = .FALSE.

         ELSE IF (XWORK .GT. XMXCKI) THEN

            POSSBL(1) = .TRUE.

            DO ISYCHO = 1,NSYM

               ISYMJ = MULD2H(ISYCHO,ISYMB)
               XNEED = XMXCKI
     &               + MINVEC(ISYCHO)*(XT1AM(ISYCHO) + NUMB*XRHF(ISYMJ))

               IF (XNEED .GT. XWORK) THEN
                  POSSBL(1) = .FALSE.
                  GO TO 1
               ENDIF

            ENDDO
    1       CONTINUE

         ELSE

            POSSBL(1) = .FALSE.

         ENDIF

      ENDIF

C     LALG=3: (i#a|j#b)
C           + L(kc,J) + L(i#a,J) + L(j#b,J) must fit in core.
C     #a, #b, and #J may be as small as one here.
C     -------------------------------------------------------

      IF ((.NOT.POSSBL(1)) .AND. (.NOT.POSSBL(2))) THEN

         POSSBL(3) = .TRUE.

         DO ISYMB = 1,NSYM
            DO ISYMA = 1,ISYMB
               ISYMIJ = MULD2H(ISYMA,ISYMB)
               XNEEDC = ZERO
               DO ISYCHO = 1,NSYM
                  ISYMI  = MULD2H(ISYCHO,ISYMA)
                  ISYMJ  = MULD2H(ISYCHO,ISYMB)
                  XNTST  = XT1AM(ISYCHO) + XRHF(ISYMI) + XRHF(ISYMJ)
                  IF (XNTST .GT. XNEEDC) XNEEDC = XNTST
               ENDDO
               XNEED = XMATIJ(ISYMIJ) + XNEEDC
               IF (XNEED .GT. XWORK) THEN
                  POSSBL(3) = .FALSE.
                  GO TO 2
               ENDIF
            ENDDO
         ENDDO
    2    CONTINUE

      ENDIF

C     Set LALG.
C     ---------

      IF (POSSBL(2)) THEN
         LALG = 2
      ELSE IF (POSSBL(1)) THEN
         LALG = 1
      ELSE IF (POSSBL(3)) THEN
         LALG = 3
      ELSE
         WRITE(LUPRI,'(//,5X,A,/,5X,A,A,A,I10,A,/)')
     &   'Insufficient memory for MP2 calculation:',
     &   'Work space passed to ',SECNAM,' is ',LWORK,' words'
         CALL QUIT('Insufficient memory for MP2 calculation')
      ENDIF

      RETURN
      END
C  /* Deck mp2cho_3 */
      SUBROUTINE MP2CHO_3(EMO,WORK,LWORK,EMP2)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose:
C        Calculate the MP2 energy correction EMP2 directly
C        from Cholesky decomposed two-electron integrals without
C        permanent storage of (ai|bj) integral intermediates.
C        The integrals are set up in batches of the two virtual
C        indices.
C
C     Notes:
C        This routine splits work space in two parts and performs
C        a triple batch; one over b and one over a in (ai|bj) integrals,
C        the third over the MO Cholesky vectors. This implies that, in
C        general, the Cholesky vectors will be read several times.
C
C        This routine requires at least in the order of
C             O*(N + 2)
C        words of memory.
C
#include <implicit.h>
      DIMENSION EMO(*), WORK(LWORK)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <priunit.h>

      INTEGER NVECTOT(8)
      INTEGER LUCHMO(8), LSYM(8)

      INTEGER IX2SQ(8)

      INTEGER BBEG, BEND
      INTEGER LVIRB(8), IOFB1(8), IX1AMB(8,8), NX1AMB(8)
      INTEGER LIBJ(8)
      INTEGER ABEG, AEND
      INTEGER LVIRA(8), IOFA1(8), IX1AMA(8,8), NX1AMA(8)

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.0D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
Casm
      LOGICAL LDUM
Casm
      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'MP2CHO_3')
Casm
#include <chomp2_b.h>
Casm

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIME = ZERO

C     Initial check of memory.
C     ------------------------

      NEED = 0
      DO ISYMB = 1,NSYM
         DO ISYMA = 1,ISYMB
            ISYMIJ = MULD2H(ISYMA,ISYMB)
            NEEDC  = 0
            DO ISYCHO = 1,NSYM
               ISYMI = MULD2H(ISYCHO,ISYMA)
               ISYMJ = MULD2H(ISYCHO,ISYMB)
               NTST  = NT1AM(ISYCHO) + NRHF(ISYMI) + NRHF(ISYMJ)
               NEEDC = MAX(NEEDC,NTST)
            ENDDO
            NTST = NMATIJ(ISYMIJ) + NEEDC
            NEED = MAX(NEED,NTST)
         ENDDO
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initializations.
C     ----------------

      EMP2   = ZERO
      IF (RSTMP2) EMP2 = OLDEN2
      MXMEMT = 0
      MXMEML = 0
      NPASS  = 0
      XLWORK = ONE*LWORK

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

Casm
C     Restart info
C     ------------

      IF (RSTMP2) THEN
         WRITE(LUPRI,'(//,A,/,A,/)') 'Restarting Cholesky MP2',
     &                               '-----------------------'
         WRITE(LUPRI,'(A,I5,A,I2,/)') 'First orbital is nr',
     &                                IFVIRB,' of symmetry',IFSYMB
         IF (OLDKNO) THEN
            WRITE(LUPRI,'(A,F15.10)')
     &           'Contribution from previous batches : ', OLDEN2
         ELSE
            WRITE(LUPRI,'(A,A)') 'Contribution from previous batches ',
     &            'not known. You must add it by yourself!!'
         END IF
         WRITE(LUPRI,*)
      END IF
Casm

C     Set # vectors.
C     --------------

      IF (CHOMO) THEN
         ITYP = 2
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NCHOMO(ISYM)
         ENDDO
      ELSE
         ITYP = 1
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NUMCHO(ISYM)
         ENDDO
      ENDIF

C     Open files for MO Cholesky vectors.
C     -----------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IOPEN,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Debug: initialize integral norm.
C     --------------------------------

      IF (LOCDBG) THEN
         XNRM = ZERO
         XDIM = ZERO
      ENDIF

C     Determine memory split.
C     -----------------------

      CALL MP2_MEMSPL3(LWORK,NVECTOT,LWRK)
      LEFT = LWORK - LWRK

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         XMB = XLWORK*WTOMB
         WRITE(LUPRI,'(12X,A,/,12X,A,I10,A,F7.1,A)')
     &   '(#ai|#bj) integrals will be calculated on the fly.',
     &   'Available memory: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(12X,A,1X,I10,A,/)')
     &   'Maximum memory for integral intermediates: ',LWRK,' words.'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Start batch loop over b.
C     ------------------------

Casm
C
C     Restart MP2
C
      IF (RSTMP2) THEN
         BBEG = 0
         DO ISYM = 1,IFSYMB-1
            BBEG = BBEG + NVIR(ISYM)
         END DO
         BBEG = BBEG +  IFVIRB
      ELSE
         BBEG = 1
      END IF
      LMPRST = -1
      LDUM   = .FALSE.
      CALL GPOPEN(LMPRST,'CHOMP2_RST','UNKNOWN','SEQUENTIAL',
     &            'FORMATTED',IDUM,LDUM)
      REWIND(LMPRST)
      WRITE(LMPRST,'(/,10X,A,/,10X,A,//)')
     &          'Cholesky MP2 restart information',
     &          '--------------------------------'
      WRITE(LMPRST,'(2X,3(A5,3X),5X,3(A5,3X),11X,A6)')
     &      'FstSy','FstOr','Bfrst','LstSy','LstOrb','Blast','Energy'
      CALL FLSHFO(LMPRST)
C
Casm
      IBATB = 0
Casm  BBEG  = 1
      NUMB  = 0
      NPASS = 0
  100 CONTINUE

         BBEG = BBEG + NUMB
         IF (BBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TBBAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATB = IBATB + 1

C        Find out how many b's can be treated,
C        estimating the number of a's along the way.
C        -------------------------------------------

         CALL GET_BTCHB(LVIRB,BBEG,NUMB,LWRK)

C        Check for errors.
C        -----------------

         BEND = BBEG + NUMB - 1
         IF (BEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'b-batch number    : ',IBATB
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First B           : ',BBEG,
     &      'Last  B           : ',BEND,
     &      'Number of virtuals: ',NVIRT
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory  : ',LWRK
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMB .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for b-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get symmetries of first and last b in batch.
C        --------------------------------------------

         ISBBEG = ISVI(BBEG)
         ISBEND = ISVI(BEND)

C        Set up index arrays.
C        --------------------

         CALL GET_INDVIR(BBEG,BEND,LVIRB,IOFB1,IX1AMB,NX1AMB)

         DO ISYIBJ = 1,NSYM
            LIBJ(ISYIBJ) = 0
            DO ISYMBJ = 1,NSYM
               ISYMI = MULD2H(ISYMBJ,ISYIBJ)
               LIBJ(ISYIBJ) = LIBJ(ISYIBJ)
     &                      + NRHF(ISYMI)*NX1AMB(ISYMBJ)
            ENDDO
         ENDDO

C        Print.
C        ------

         IF (LOCDBG) THEN
            WRITE(LUPRI,'(/,A,2I9)')
     &      'First and last b: ',BBEG,BEND
            WRITE(LUPRI,'(A,8I9)')
     &      'IOFB1  :',(IOFB1(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'LVIRB  :',(LVIRB(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NVIR   :',(NVIR(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NX1AMB :',(NX1AMB(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NT1AM  :',(NT1AM(I),I=1,NSYM)
            DO J = 1,NSYM
               WRITE(LUPRI,'(A,8I9)')
     &         'IX1AMB :',(IX1AMB(I,J),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'IT1AM  :',(IT1AM(I,J),I=1,NSYM)
            ENDDO
            WRITE(LUPRI,'(A,8I9)')
     &      'LIBJ   :',(LIBJ(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NCKI   :',(NCKI(I),I=1,NSYM)
            WRITE(LUPRI,*)
         ENDIF

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(20X,A,I6,A,/,20X,A)')
     &      'b-batch number ',IBATB,':',
     &      '======================'
            WRITE(LUPRI,'(20X,A,I6)')
     &      '#b   : ',NUMB
            WRITE(LUPRI,'(20X,A,I6,A,I1,/,20X,A,I6,A,I1,/)')
     &      'First: ',IOFB1(ISBBEG),'   symmetry: ',ISBBEG,
     &      'Last : ',IOFB1(ISBEND)+LVIRB(ISBEND)-1,'   symmetry: ',
     &                                                  ISBEND
            WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &      '   #a        First           Last         Mem. Usage',
     &      '        Time   ',
     &      '          (a, sym. a)     (a, sym. a)        (%)    ',
     &      '      (seconds)'
            WRITE(LUPRI,'(5X,A,A)')
     &      '----------------------------------------------------',
     &      '---------------'
            CALL FLSHFO(LUPRI)
         ENDIF

C        Start batch loop over a.
C        ------------------------

         IBATA = 0
         ABEG  = 1
         NUMA  = 0
  110    CONTINUE

            ABEG = ABEG + NUMA
            IF (ABEG .GT. NVIRT) GO TO 190

C           Time this batch.
C           ----------------

            TABAT = SECOND()

C           Update batch counter.
C           ---------------------

            IBATA = IBATA + 1

C           Find out how many a's can be treated.
C           -------------------------------------

            CALL GET_BTCHA(LIBJ,NX1AMB,LVIRA,ABEG,NUMA,MEM,LWRK,LWORK)

C           Check for errors.
C           -----------------

            AEND = ABEG + NUMA - 1
            IF (AEND .GT. NVIRT) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Batch error in ',SECNAM,' !'
               WRITE(LUPRI,'(5X,A)')
     &         '- dumping presumably most relevant info:'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'B-batch number    : ',IBATB,
     &         'A-batch number    : ',IBATA
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'First B           : ',BBEG,
     &         'Last  B           : ',BEND
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &         'First A           : ',ABEG,
     &         'Last  A           : ',AEND,
     &         'Number of virtuals: ',NVIRT
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Available memory  : ',LWRK,
     &         'Needed    memory  : ',MEM
               CALL QUIT('Batch error in '//SECNAM)
            ENDIF

            IF (NUMA .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory for a-batch in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10)')
     &         'Available memory: ',LWRK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

C           Get symmetries of first and last a in batch.
C           --------------------------------------------

            ISABEG = ISVI(ABEG)
            ISAEND = ISVI(AEND)

C           Set up index arrays.
C           --------------------

            CALL GET_INDVIR(ABEG,AEND,LVIRA,IOFA1,IX1AMA,NX1AMA)

C           Set up (#ai|#bj) index arrays.
C           ------------------------------

            ICOUNT = 0
            DO ISYMBJ = 1,NSYM
               ISYMAI = ISYMBJ
               IX2SQ(ISYMBJ) = ICOUNT
               ICOUNT = ICOUNT + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
            ENDDO
            NX2SQ = ICOUNT

            IF (NX2SQ .NE. MEM) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Error detected in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &         'NX2SQ is calculated to be ',NX2SQ,
     &         'MEM determines alloc. as  ',MEM,
     &         'These numbers must be equal....'
               CALL QUIT('Error in '//SECNAM)
            ENDIF

            IF (LOCDBG) THEN
               WRITE(LUPRI,'(/,A,2I9)')
     &         'First and last a: ',ABEG,AEND
               WRITE(LUPRI,'(A,8I9)')
     &         'IOFA1  :',(IOFA1(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'LVIRA  :',(LVIRA(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NVIR   :',(NVIR(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NX1AMA :',(NX1AMA(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NT1AM  :',(NT1AM(I),I=1,NSYM)
               DO J = 1,NSYM
                  WRITE(LUPRI,'(A,8I9)')
     &            'IX1AMA :',(IX1AMA(I,J),I=1,NSYM)
                  WRITE(LUPRI,'(A,8I9)')
     &            'IT1AM  :',(IT1AM(I,J),I=1,NSYM)
               ENDDO
               WRITE(LUPRI,'(A,8I9)')
     &         'IX2SQ  :',(IX2SQ(I),I=1,NSYM)
               WRITE(LUPRI,'(A,I9)')
     &         'NX2SQ  :',NX2SQ
               WRITE(LUPRI,'(A,I9,/)')
     &         'NT2SQ  :',NT2SQ(1)
            ENDIF

C           Complete allocation for integral intermediates.
C           -----------------------------------------------

            KXINT = 1
            KEND1 = KXINT + NX2SQ
            LWRK1 = LWORK - KEND1 + 1

C           Initialize integral intermediates.
C           ----------------------------------

            CALL DZERO(WORK(KXINT),NX2SQ)

C           Start loop over Cholesky symmetries.
C           ------------------------------------

            DO ISYCHO = 1,NSYM

C              If nothing to do here, cycle.
C              -----------------------------

               IF (NVECTOT(ISYCHO) .LE. 0) GO TO 999
               IF (NX1AMA(ISYCHO)  .LE. 0) GO TO 999
               IF (NX1AMB(ISYCHO)  .LE. 0) GO TO 999

C              Allocation for 1 full Cholesky vector L(ck).
C              --------------------------------------------

               KCHOL = KEND1
               KEND2 = KCHOL + NT1AM(ISYCHO)
               LWRK2 = LWORK - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Insufficient memory for full Cholesky vector in ',
     &            SECNAM
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'Need (more than) : ',KEND2-1,
     &            'Available        : ',LWORK
                  WRITE(LUPRI,'(5X,A,9X,I1,/,5X,A,I10)')
     &            'Cholesky symmetry: ',ISYCHO,
     &            'Size of vector   : ',NT1AM(ISYCHO)
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Size of (#ai|#bj): ',NX2SQ
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'b-batch number   : ',IBATB,
     &            'a-batch number   : ',IBATA
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'First B          : ',BBEG,
     &            'Last  B          : ',BEND
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'First A          : ',ABEG,
     &            'Last  A          : ',AEND
                  WRITE(LUPRI,'(5X,A,/)')
     &            'This error should not occur; may indicate a bug!'
                  CALL QUIT('Insufficient memory in '//SECNAM)
               ENDIF

C              Set up Cholesky batch in remaining memory.
C              ------------------------------------------

               MINMEM = NX1AMA(ISYCHO) + NX1AMB(ISYCHO)
               NVEC   = MIN(LWRK2/MINMEM,NVECTOT(ISYCHO))
               IF (NVEC .LE. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Insufficient memory for Cholesky batch in ',SECNAM
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Minimum memory required for vectors: ',MINMEM,
     &            'Memory available for batch         : ',LWRK2
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,A,I1,A)')
     &            'Memory available in total          : ',LWORK,
     &            'Total number of vectors            : ',
     &            NVECTOT(ISYCHO),' (symmetry ',ISYCHO,')'
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Size of (#ai|#bj)                  : ',NX2SQ
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Size of full vec.                  : ',NT1AM(ISYCHO)
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'b-batch number                     : ',IBATB,
     &            'a-batch number                     : ',IBATA
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'First B                            : ',BBEG,
     &            'Last  B                            : ',BEND
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &            'First A                            : ',ABEG,
     &            'Last  A                            : ',AEND
                  CALL QUIT('Insufficient memory in '//SECNAM)
               ENDIF
               NBATCH = (NVECTOT(ISYCHO) - 1)/NVEC + 1

               DO IBATCH = 1,NBATCH

                  NUMV = NVEC
                  IF (IBATCH .EQ. NBATCH) THEN
                     NUMV = NVECTOT(ISYCHO) - NVEC*(NBATCH - 1)
                  ENDIF

                  JVEC1 = NVEC*(IBATCH - 1) + 1

                  KCHBJ = KEND2
                  KCHAI = KCHBJ + NX1AMB(ISYCHO)*NUMV
                  KEND3 = KCHAI + NX1AMA(ISYCHO)*NUMV

                  KLAST  = KEND3 - 1
                  MXMEML = MAX(MXMEML,KLAST)

C                 Read vectors and copy out sub-blocks
C                 as L(#J,#bj) and L(#J,#ai).
C                 ------------------------------------

                  DO IVEC = 1,NUMV

                     DTIME = SECOND()
                     JVEC  = JVEC1 + IVEC - 1
                     CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYCHO),1,
     &                               JVEC,LUCHMO(ISYCHO))
                     DTIME = SECOND() - DTIME
                     TIMR  = TIMR     + DTIME

                     DTIME = SECOND()
                     DO ISYMJ = 1,NSYM

                        ISYMB = MULD2H(ISYMJ,ISYCHO)

                        IF (LVIRB(ISYMB) .GT. 0) THEN
                           DO J = 1,NRHF(ISYMJ)

                              LBJ1 = IX1AMB(ISYMB,ISYMJ)
     &                             + LVIRB(ISYMB)*(J - 1) + 1

                              KOFF1 = KCHOL + IT1AM(ISYMB,ISYMJ)
     &                              + NVIR(ISYMB)*(J - 1)
     &                              + IOFB1(ISYMB) - 1
                              KOFF2 = KCHBJ + NUMV*(LBJ1 - 1) + IVEC - 1

                              CALL DCOPY(LVIRB(ISYMB),WORK(KOFF1),1,
     &                                                WORK(KOFF2),NUMV)

                           ENDDO
                        ENDIF

                        ISYMI = ISYMJ
                        ISYMA = ISYMB

                        IF (LVIRA(ISYMA) .GT. 0) THEN
                           DO I = 1,NRHF(ISYMI)

                              LAI1 = IX1AMA(ISYMA,ISYMI)
     &                             + LVIRA(ISYMA)*(I - 1) + 1

                              KOFF1 = KCHOL + IT1AM(ISYMA,ISYMI)
     &                              + NVIR(ISYMA)*(I - 1)
     &                              + IOFA1(ISYMA) - 1
                              KOFF2 = KCHAI + NUMV*(LAI1 - 1) + IVEC - 1

                              CALL DCOPY(LVIRA(ISYMA),WORK(KOFF1),1,
     &                                                WORK(KOFF2),NUMV)

                           ENDDO
                        ENDIF

                     ENDDO
                     DTIME = SECOND() - DTIME
                     TIMS  = TIMS     + DTIME

                  ENDDO

C                 Calculate (#ai|#bj) contribution.
C                 ---------------------------------

                  DTIME = SECOND()

                  NTOAI = MAX(NX1AMA(ISYCHO),1)

                  KOFFX = KXINT + IX2SQ(ISYCHO)

                  CALL DGEMM('T','N',NX1AMA(ISYCHO),NX1AMB(ISYCHO),NUMV,
     &                       ONE,WORK(KCHAI),NUMV,WORK(KCHBJ),NUMV,
     &                       ONE,WORK(KOFFX),NTOAI)

                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME

               ENDDO

  999          CONTINUE

            ENDDO

C           Debug: Calculate norm of (ai|bj).
C           ---------------------------------

            IF (LOCDBG) THEN
               XNRM = XNRM + DDOT(NX2SQ,WORK(KXINT),1,WORK(KXINT),1)
               XDIM = XDIM + ONE*NX2SQ
            ENDIF

C           Calculate contribution to EMP2.
C           -------------------------------

            DTIME = SECOND()
            DO ISYMBJ = 1,NSYM

               ISYMAI = ISYMBJ

               DO ISYMJ = 1,NSYM

                  ISYMB = MULD2H(ISYMJ,ISYMBJ)
                  IF (LVIRB(ISYMB) .GT. 0) THEN

                     DO J = 1,NRHF(ISYMJ)

                        KOFFJ = IRHF(ISYMJ) + J

                        DO LB = 1,LVIRB(ISYMB)

                           B     = IOFB1(ISYMB) + LB - 1
                           KOFFB = IVIR(ISYMB)  + B

                           LBJ = IX1AMB(ISYMB,ISYMJ)
     &                         + LVIRB(ISYMB)*(J - 1) + LB

                           DO ISYMI = 1,NSYM

                              ISYMA = MULD2H(ISYMI,ISYMAI)
                              IF (LVIRA(ISYMA) .GT. 0) THEN

                                 ISYMAJ = MULD2H(ISYMA,ISYMJ)
                                 ISYMBI = ISYMAJ

                                 DO I = 1,NRHF(ISYMI)

                                    KOFFI = IRHF(ISYMI) + I

                                    LBI = IX1AMB(ISYMB,ISYMI)
     &                                  + LVIRB(ISYMB)*(I - 1) + LB

                                    DO LA = 1,LVIRA(ISYMA)

                                       A     = IOFA1(ISYMA) + LA - 1
                                       KOFFA = IVIR(ISYMA)  + A

                                       LAI = IX1AMA(ISYMA,ISYMI)
     &                                     + LVIRA(ISYMA)*(I - 1) + LA
                                       LAJ = IX1AMA(ISYMA,ISYMJ)
     &                                     + LVIRA(ISYMA)*(J - 1) + LA

                                       NAIBJ = KXINT + IX2SQ(ISYMBJ)
     &                                       + NX1AMA(ISYMAI)*(LBJ - 1)
     &                                       + LAI - 1
                                       NAJBI = KXINT + IX2SQ(ISYMBI)
     &                                       + NX1AMA(ISYMAJ)*(LBI - 1)
     &                                       + LAJ - 1

                                       DENOM = EMO(KOFFA) - EMO(KOFFI)
     &                                       + EMO(KOFFB) - EMO(KOFFJ)
                                       TAIBJ = -WORK(NAIBJ)/DENOM
                                       XAIBJ = TWO*WORK(NAIBJ)
     &                                       - WORK(NAJBI)

                                       EMP2 = EMP2 + TAIBJ*XAIBJ

                                    ENDDO

                                 ENDDO

                              ENDIF

                           ENDDO

                        ENDDO

                     ENDDO

                  ENDIF

               ENDDO

            ENDDO
            DTIME = SECOND() - DTIME
            TIME  = TIME     + DTIME

C           Print info from this a-batch.
C           -----------------------------

            IF (IPRINT .GE. INFO) THEN
               TABAT = SECOND() - TABAT
               XMUSE = (D100*MXMEML)/XLWORK
               LAFST = IOFA1(ISABEG)
               LALST = IOFA1(ISAEND) + LVIRA(ISAEND) - 1
         WRITE(LUPRI,'(I11,I10,1X,I1,8X,I6,1X,I1,9X,F6.2,F18.2,F15.8)')
     &         NUMA,LAFST,ISABEG,LALST,ISAEND,XMUSE,TABAT,EMP2
               CALL FLSHFO(LUPRI)
            ENDIF

C           Update memory info.
C           -------------------

            MXMEMT = MAX(MXMEMT,MXMEML)
            MXMEML = 0

C           Update NPASS and go to next a-batch.
C           (Which will exit immediately if no more a's.)
C           ---------------------------------------------

            NPASS = NPASS + 1
            GO TO 110

C        End of #a loop.
C        ---------------

  190    CONTINUE

C
C        Restart info
C
         WRITE(LMPRST,'(2(3X,I2,3X,I6,2X,I6,7X),F17.10)')
     &         ISBBEG,IOFB1(ISBBEG),BBEG,
     &         ISBEND,IOFB1(ISBEND)+LVIRB(ISBEND)-1,BEND,EMP2
         CALL FLSHFO(LMPRST)
 
C        Print info for this b-batch and got to next.
C        (Which will exit immediately if no more b's.)
C        ---------------------------------------------

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(5X,A,A,/)')
     &      '----------------------------------------------------',
     &      '---------------'
            TBBAT = SECOND() - TBBAT
            WRITE(LUPRI,'(20X,A,F10.2,A,/)')
     &      'Total time for this b-batch: ',TBBAT,' seconds'
            CALL FLSHFO(LUPRI)
         ENDIF

         GO TO 100

C     Exit b-batch: calculation done.
C     -------------------------------

  200 CONTINUE

C     Close files for MO Cholesky vectors.
C     ------------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IKEEP,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Print section.
C     --------------

      IF ((IPRINT.GT.0) .OR. LOCDBG) THEN
         IF (LOCDBG) THEN
            XNRM = DSQRT(XNRM)
            CALL HEADER('Debug info from '//SECNAM,-1)
            WRITE(LUPRI,'(5X,A,1P,D30.15,/,5X,A,1P,D30.15,/)')
     &      'Norm of (ai|bj) [full square]: ',XNRM,
     &      'Size of (ai|bj) [db. prec.]  : ',XDIM
         ENDIF
         TIMT  = SECOND() - TIMT
         XMUSE = (D100*MXMEMT)/XLWORK
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &   'Memory reserved for integrals        : ',LWRK,
     &   'Memory reserved for Cholesky vectors : ',LEFT,
     &   'Total memory available               : ',LWORK
         WRITE(LUPRI,'(5X,A,4X,F6.2,A)')
     &   'Maximum memory usage                 : ',XMUSE,' %'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors: ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,I10)')
     &   ' - passes through Cholesky file(s)   : ',NPASS
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sort of MO vectors     : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating (#ai|#bj)  : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating MP2 energy : ',TIME,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used in total                   : ',TIMT,' seconds'
            WRITE(LUPRI,'(5X,A,/)')
     &      '(Debug calculations and printing included.)'
         ELSE
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Time used in total                   : ',TIMT,' seconds'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck mp2cho_3o */
      SUBROUTINE MP2CHO_3O(EMO,WORK,LWORK,EMP2)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose:
C        Calculate the MP2 energy correction EMP2 directly
C        from Cholesky decomposed two-electron integrals without
C        permanent storage of (ai|bj) integral intermediates.
C        The integrals are set up in batches of the two virtual
C        indices.
C
C     Notes:
C        This routine splits work space in two parts and performs
C        a triple batch; one over b and one over a in (ai|bj) integrals,
C        the third over the MO Cholesky vectors. This implies that, in
C        general, the Cholesky vectors will be read several times.
C
C        This routine requires at least in the order of
C             O*(N + 2)
C        words of memory.
C
#include <implicit.h>
      DIMENSION EMO(*), WORK(LWORK)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <priunit.h>

      INTEGER NVECTOT(8)
      INTEGER LUCHMO(8), LSYM(8)

      INTEGER IX2SQ(8)

      INTEGER BBEG, BEND
      INTEGER LVIRB(8), IOFB1(8), IX1AMB(8,8), NX1AMB(8)
      INTEGER LIBJ(8)
      INTEGER ABEG, AEND
      INTEGER LVIRA(8), IOFA1(8), IX1AMA(8,8), NX1AMA(8)

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.0D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'MP2CHO_3O')

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIME = ZERO

C     Initial check of memory.
C     ------------------------

      NEED = 0
      DO ISYMB = 1,NSYM
         DO ISYMA = 1,ISYMB
            ISYMIJ = MULD2H(ISYMA,ISYMB)
            NEEDC  = 0
            DO ISYCHO = 1,NSYM
               ISYMI = MULD2H(ISYCHO,ISYMA)
               ISYMJ = MULD2H(ISYCHO,ISYMB)
               NTST  = NT1AM(ISYCHO) + NRHF(ISYMI) + NRHF(ISYMJ)
               NEEDC = MAX(NEEDC,NTST)
            ENDDO
            NTST = NMATIJ(ISYMIJ) + NEEDC
            NEED = MAX(NEED,NTST)
         ENDDO
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initializations.
C     ----------------

      EMP2   = ZERO
      MXMEMT = 0
      MXMEML = 0
      NPASS  = 0
      XLWORK = ONE*LWORK

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

C     Set # vectors.
C     --------------

      IF (CHOMO) THEN
         ITYP = 2
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NCHOMO(ISYM)
         ENDDO
      ELSE
         ITYP = 1
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NUMCHO(ISYM)
         ENDDO
      ENDIF

C     Open files for MO Cholesky vectors.
C     -----------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IOPEN,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Debug: initialize integral norm.
C     --------------------------------

      IF (LOCDBG) THEN
         XNRM = ZERO
         XDIM = ZERO
      ENDIF

C     Determine memory split.
C     -----------------------

      CALL MP2_MEMSPL3(LWORK,NVECTOT,LWRK)
      LEFT = LWORK - LWRK

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         XMB = XLWORK*WTOMB
         WRITE(LUPRI,'(12X,A,/,12X,A,I10,A,F7.1,A)')
     &   '(#ai|#bj) integrals will be calculated on the fly.',
     &   'Available memory: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(12X,A,1X,I10,A,/)')
     &   'Maximum memory for integral intermediates: ',LWRK,' words.'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Start batch loop over b.
C     ------------------------

      IBATB = 0
      BBEG  = 1
      NUMB  = 0
      NPASS = 0
  100 CONTINUE

         BBEG = BBEG + NUMB
         IF (BBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TBBAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATB = IBATB + 1

C        Find out how many b's can be treated,
C        estimating the number of a's along the way.
C        -------------------------------------------

         CALL GET_BTCHB(LVIRB,BBEG,NUMB,LWRK)

C        Check for errors.
C        -----------------

         BEND = BBEG + NUMB - 1
         IF (BEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'b-batch number    : ',IBATB
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First B           : ',BBEG,
     &      'Last  B           : ',BEND,
     &      'Number of virtuals: ',NVIRT
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory  : ',LWRK
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMB .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for b-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get symmetries of first and last b in batch.
C        --------------------------------------------

         ISBBEG = ISVI(BBEG)
         ISBEND = ISVI(BEND)

C        Set up index arrays.
C        --------------------

         CALL GET_INDVIR(BBEG,BEND,LVIRB,IOFB1,IX1AMB,NX1AMB)

         DO ISYIBJ = 1,NSYM
            LIBJ(ISYIBJ) = 0
            DO ISYMBJ = 1,NSYM
               ISYMI = MULD2H(ISYMBJ,ISYIBJ)
               LIBJ(ISYIBJ) = LIBJ(ISYIBJ)
     &                      + NRHF(ISYMI)*NX1AMB(ISYMBJ)
            ENDDO
         ENDDO

C        Print.
C        ------

         IF (LOCDBG) THEN
            WRITE(LUPRI,'(/,A,2I9)')
     &      'First and last b: ',BBEG,BEND
            WRITE(LUPRI,'(A,8I9)')
     &      'IOFB1  :',(IOFB1(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'LVIRB  :',(LVIRB(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NVIR   :',(NVIR(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NX1AMB :',(NX1AMB(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NT1AM  :',(NT1AM(I),I=1,NSYM)
            DO J = 1,NSYM
               WRITE(LUPRI,'(A,8I9)')
     &         'IX1AMB :',(IX1AMB(I,J),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'IT1AM  :',(IT1AM(I,J),I=1,NSYM)
            ENDDO
            WRITE(LUPRI,'(A,8I9)')
     &      'LIBJ   :',(LIBJ(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NCKI   :',(NCKI(I),I=1,NSYM)
            WRITE(LUPRI,*)
         ENDIF

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(20X,A,I6,A,/,20X,A)')
     &      'b-batch number ',IBATB,':',
     &      '======================'
            WRITE(LUPRI,'(20X,A,I6)')
     &      '#b   : ',NUMB
            WRITE(LUPRI,'(20X,A,I6,A,I1,/,20X,A,I6,A,I1,/)')
     &      'First: ',IOFB1(ISBBEG),'   symmetry: ',ISBBEG,
     &      'Last : ',IOFB1(ISBEND)+LVIRB(ISBEND)-1,'   symmetry: ',
     &                                                  ISBEND
            WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &      '   #a        First           Last         Mem. Usage',
     &      '        Time   ',
     &      '          (a, sym. a)     (a, sym. a)        (%)    ',
     &      '      (seconds)'
            WRITE(LUPRI,'(5X,A,A)')
     &      '----------------------------------------------------',
     &      '---------------'
            CALL FLSHFO(LUPRI)
         ENDIF

C        Start batch loop over a.
C        ------------------------

         IBATA = 0
         ABEG  = 1
         NUMA  = 0
  110    CONTINUE

            ABEG = ABEG + NUMA
            IF (ABEG .GT. NVIRT) GO TO 190

C           Time this batch.
C           ----------------

            TABAT = SECOND()

C           Update batch counter.
C           ---------------------

            IBATA = IBATA + 1

C           Find out how many a's can be treated.
C           -------------------------------------

            CALL GET_BTCHA(LIBJ,NX1AMB,LVIRA,ABEG,NUMA,MEM,LWRK,LWORK)

C           Check for errors.
C           -----------------

            AEND = ABEG + NUMA - 1
            IF (AEND .GT. NVIRT) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Batch error in ',SECNAM,' !'
               WRITE(LUPRI,'(5X,A)')
     &         '- dumping presumably most relevant info:'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'B-batch number    : ',IBATB,
     &         'A-batch number    : ',IBATA
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'First B           : ',BBEG,
     &         'Last  B           : ',BEND
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &         'First A           : ',ABEG,
     &         'Last  A           : ',AEND,
     &         'Number of virtuals: ',NVIRT
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Available memory  : ',LWRK,
     &         'Needed    memory  : ',MEM
               CALL QUIT('Batch error in '//SECNAM)
            ENDIF

            IF (NUMA .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory for a-batch in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10)')
     &         'Available memory: ',LWRK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

C           Get symmetries of first and last a in batch.
C           --------------------------------------------

            ISABEG = ISVI(ABEG)
            ISAEND = ISVI(AEND)

C           Set up index arrays.
C           --------------------

            CALL GET_INDVIR(ABEG,AEND,LVIRA,IOFA1,IX1AMA,NX1AMA)

C           Set up (#ai|#bj) index arrays.
C           ------------------------------

            ICOUNT = 0
            DO ISYMBJ = 1,NSYM
               ISYMAI = ISYMBJ
               IX2SQ(ISYMBJ) = ICOUNT
               ICOUNT = ICOUNT + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
            ENDDO
            NX2SQ = ICOUNT

            IF (NX2SQ .NE. MEM) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Error detected in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &         'NX2SQ is calculated to be ',NX2SQ,
     &         'MEM determines alloc. as  ',MEM,
     &         'These numbers must be equal....'
               CALL QUIT('Error in '//SECNAM)
            ENDIF

            IF (LOCDBG) THEN
               WRITE(LUPRI,'(/,A,2I9)')
     &         'First and last a: ',ABEG,AEND
               WRITE(LUPRI,'(A,8I9)')
     &         'IOFA1  :',(IOFA1(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'LVIRA  :',(LVIRA(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NVIR   :',(NVIR(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NX1AMA :',(NX1AMA(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NT1AM  :',(NT1AM(I),I=1,NSYM)
               DO J = 1,NSYM
                  WRITE(LUPRI,'(A,8I9)')
     &            'IX1AMA :',(IX1AMA(I,J),I=1,NSYM)
                  WRITE(LUPRI,'(A,8I9)')
     &            'IT1AM  :',(IT1AM(I,J),I=1,NSYM)
               ENDDO
               WRITE(LUPRI,'(A,8I9)')
     &         'IX2SQ  :',(IX2SQ(I),I=1,NSYM)
               WRITE(LUPRI,'(A,I9)')
     &         'NX2SQ  :',NX2SQ
               WRITE(LUPRI,'(A,I9,/)')
     &         'NT2SQ  :',NT2SQ(1)
            ENDIF

C           Complete allocation for integral intermediates.
C           -----------------------------------------------

            KXINT = 1
            KEND1 = KXINT + NX2SQ
            LWRK1 = LWORK - KEND1 + 1

C           Initialize integral intermediates.
C           ----------------------------------

            CALL DZERO(WORK(KXINT),NX2SQ)

C           Start loop over Cholesky symmetries.
C           ------------------------------------

            DO ISYCHO = 1,NSYM

C              If nothing to do here, cycle.
C              -----------------------------

               IF (NVECTOT(ISYCHO) .LE. 0) GO TO 999
               IF (NX1AMA(ISYCHO)  .LE. 0) GO TO 999
               IF (NX1AMB(ISYCHO)  .LE. 0) GO TO 999

C              Allocation for 1 full Cholesky vector L(ck).
C              --------------------------------------------

               KCHOL = KEND1
               KEND2 = KCHOL + NT1AM(ISYCHO)
               LWRK2 = LWORK - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Insufficient memory for full Cholesky vector in ',
     &            SECNAM
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'Need (more than) : ',KEND2-1,
     &            'Available        : ',LWORK
                  WRITE(LUPRI,'(5X,A,9X,I1,/,5X,A,I10)')
     &            'Cholesky symmetry: ',ISYCHO,
     &            'Size of vector   : ',NT1AM(ISYCHO)
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Size of (#ai|#bj): ',NX2SQ
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'b-batch number   : ',IBATB,
     &            'a-batch number   : ',IBATA
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'First B          : ',BBEG,
     &            'Last  B          : ',BEND
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'First A          : ',ABEG,
     &            'Last  A          : ',AEND
                  WRITE(LUPRI,'(5X,A,/)')
     &            'This error should not occur; may indicate a bug!'
                  CALL QUIT('Insufficient memory in '//SECNAM)
               ENDIF

C              Set up Cholesky batch in remaining memory.
C              ------------------------------------------

               MINMEM = NX1AMA(ISYCHO) + NX1AMB(ISYCHO)
               NVEC   = MIN(LWRK2/MINMEM,NVECTOT(ISYCHO))
               IF (NVEC .LE. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Insufficient memory for Cholesky batch in ',SECNAM
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Minimum memory required for vectors: ',MINMEM,
     &            'Memory available for batch         : ',LWRK2
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,A,I1,A)')
     &            'Memory available in total          : ',LWORK,
     &            'Total number of vectors            : ',
     &            NVECTOT(ISYCHO),' (symmetry ',ISYCHO,')'
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Size of (#ai|#bj)                  : ',NX2SQ
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Size of full vec.                  : ',NT1AM(ISYCHO)
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'b-batch number                     : ',IBATB,
     &            'a-batch number                     : ',IBATA
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &            'First B                            : ',BBEG,
     &            'Last  B                            : ',BEND
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &            'First A                            : ',ABEG,
     &            'Last  A                            : ',AEND
                  CALL QUIT('Insufficient memory in '//SECNAM)
               ENDIF
               NBATCH = (NVECTOT(ISYCHO) - 1)/NVEC + 1

               DO IBATCH = 1,NBATCH

                  NUMV = NVEC
                  IF (IBATCH .EQ. NBATCH) THEN
                     NUMV = NVECTOT(ISYCHO) - NVEC*(NBATCH - 1)
                  ENDIF

                  JVEC1 = NVEC*(IBATCH - 1) + 1

                  KCHBJ = KEND2
                  KCHAI = KCHBJ + NX1AMB(ISYCHO)*NUMV
                  KEND3 = KCHAI + NX1AMA(ISYCHO)*NUMV

                  KLAST  = KEND3 - 1
                  MXMEML = MAX(MXMEML,KLAST)

C                 Read vectors and copy out sub-blocks.
C                 -------------------------------------

                  DO IVEC = 1,NUMV

                     DTIME = SECOND()
                     JVEC  = JVEC1 + IVEC - 1
                     CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYCHO),1,
     &                               JVEC,LUCHMO(ISYCHO))
                     DTIME = SECOND() - DTIME
                     TIMR  = TIMR     + DTIME

                     DTIME = SECOND()
                     DO ISYMJ = 1,NSYM

                        ISYMB = MULD2H(ISYMJ,ISYCHO)

                        IF (LVIRB(ISYMB) .GT. 0) THEN
                           DO J = 1,NRHF(ISYMJ)

                              KOFF1 = KCHOL + IT1AM(ISYMB,ISYMJ)
     &                              + NVIR(ISYMB)*(J - 1)
     &                              + IOFB1(ISYMB) - 1
                              KOFF2 = KCHBJ + NX1AMB(ISYCHO)*(IVEC - 1)
     &                              + IX1AMB(ISYMB,ISYMJ)
     &                              + LVIRB(ISYMB)*(J - 1)

                              CALL DCOPY(LVIRB(ISYMB),WORK(KOFF1),1,
     &                                                WORK(KOFF2),1)

                           ENDDO
                        ENDIF

                        ISYMI = ISYMJ
                        ISYMA = ISYMB

                        IF (LVIRA(ISYMA) .GT. 0) THEN
                           DO I = 1,NRHF(ISYMI)

                              KOFF1 = KCHOL + IT1AM(ISYMA,ISYMI)
     &                              + NVIR(ISYMA)*(I - 1)
     &                              + IOFA1(ISYMA) - 1
                              KOFF2 = KCHAI + NX1AMA(ISYCHO)*(IVEC - 1)
     &                              + IX1AMA(ISYMA,ISYMI)
     &                              + LVIRA(ISYMA)*(I - 1)

                              CALL DCOPY(LVIRA(ISYMA),WORK(KOFF1),1,
     &                                                WORK(KOFF2),1)

                           ENDDO
                        ENDIF

                     ENDDO
                     DTIME = SECOND() - DTIME
                     TIMS  = TIMS     + DTIME

                  ENDDO

C                 Calculate (#ai|#bj) contribution.
C                 ---------------------------------

                  DTIME = SECOND()

                  NTOAI = MAX(NX1AMA(ISYCHO),1)
                  NTOBJ = MAX(NX1AMB(ISYCHO),1)

                  KOFFX = KXINT + IX2SQ(ISYCHO)

                  CALL DGEMM('N','T',NX1AMA(ISYCHO),NX1AMB(ISYCHO),NUMV,
     &                       ONE,WORK(KCHAI),NTOAI,WORK(KCHBJ),NTOBJ,
     &                       ONE,WORK(KOFFX),NTOAI)

                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME

               ENDDO

  999          CONTINUE

            ENDDO

C           Debug: Calculate norm of (ai|bj).
C           ---------------------------------

            IF (LOCDBG) THEN
               XNRM = XNRM + DDOT(NX2SQ,WORK(KXINT),1,WORK(KXINT),1)
               XDIM = XDIM + ONE*NX2SQ
            ENDIF

C           Calculate contribution to EMP2.
C           -------------------------------

            DTIME = SECOND()
            DO ISYMBJ = 1,NSYM

               ISYMAI = ISYMBJ

               DO ISYMJ = 1,NSYM

                  ISYMB = MULD2H(ISYMJ,ISYMBJ)
                  IF (LVIRB(ISYMB) .GT. 0) THEN

                     DO J = 1,NRHF(ISYMJ)

                        KOFFJ = IRHF(ISYMJ) + J

                        DO LB = 1,LVIRB(ISYMB)

                           B     = IOFB1(ISYMB) + LB - 1
                           KOFFB = IVIR(ISYMB)  + B

                           LBJ = IX1AMB(ISYMB,ISYMJ)
     &                         + LVIRB(ISYMB)*(J - 1) + LB

                           DO ISYMI = 1,NSYM

                              ISYMA = MULD2H(ISYMI,ISYMAI)
                              IF (LVIRA(ISYMA) .GT. 0) THEN

                                 ISYMAJ = MULD2H(ISYMA,ISYMJ)
                                 ISYMBI = ISYMAJ

                                 DO I = 1,NRHF(ISYMI)

                                    KOFFI = IRHF(ISYMI) + I

                                    LBI = IX1AMB(ISYMB,ISYMI)
     &                                  + LVIRB(ISYMB)*(I - 1) + LB

                                    DO LA = 1,LVIRA(ISYMA)

                                       A     = IOFA1(ISYMA) + LA - 1
                                       KOFFA = IVIR(ISYMA)  + A

                                       LAI = IX1AMA(ISYMA,ISYMI)
     &                                     + LVIRA(ISYMA)*(I - 1) + LA
                                       LAJ = IX1AMA(ISYMA,ISYMJ)
     &                                     + LVIRA(ISYMA)*(J - 1) + LA

                                       NAIBJ = KXINT + IX2SQ(ISYMBJ)
     &                                       + NX1AMA(ISYMAI)*(LBJ - 1)
     &                                       + LAI - 1
                                       NAJBI = KXINT + IX2SQ(ISYMBI)
     &                                       + NX1AMA(ISYMAJ)*(LBI - 1)
     &                                       + LAJ - 1

                                       DENOM = EMO(KOFFA) - EMO(KOFFI)
     &                                       + EMO(KOFFB) - EMO(KOFFJ)
                                       TAIBJ = -WORK(NAIBJ)/DENOM
                                       XAIBJ = TWO*WORK(NAIBJ)
     &                                       - WORK(NAJBI)

                                       EMP2 = EMP2 + TAIBJ*XAIBJ

                                    ENDDO

                                 ENDDO

                              ENDIF

                           ENDDO

                        ENDDO

                     ENDDO

                  ENDIF

               ENDDO

            ENDDO
            DTIME = SECOND() - DTIME
            TIME  = TIME     + DTIME

C           Print info from this a-batch.
C           -----------------------------

            IF (IPRINT .GE. INFO) THEN
               TABAT = SECOND() - TABAT
               XMUSE = (D100*MXMEML)/XLWORK
               LAFST = IOFA1(ISABEG)
               LALST = IOFA1(ISAEND) + LVIRA(ISAEND) - 1
         WRITE(LUPRI,'(5X,I6,4X,I6,1X,I1,8X,I6,1X,I1,9X,F6.2,8X,F10.2)')
     &         NUMA,LAFST,ISABEG,LALST,ISAEND,XMUSE,TABAT
               CALL FLSHFO(LUPRI)
            ENDIF

C           Update memory info.
C           -------------------

            MXMEMT = MAX(MXMEMT,MXMEML)
            MXMEML = 0

C           Update NPASS and go to next a-batch.
C           (Which will exit immediately if no more a's.)
C           ---------------------------------------------

            NPASS = NPASS + 1
            GO TO 110

C        End of #a loop.
C        ---------------

  190    CONTINUE

C        Print info for this b-batch and got to next.
C        (Which will exit immediately if no more b's.)
C        ---------------------------------------------

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(5X,A,A,/)')
     &      '----------------------------------------------------',
     &      '---------------'
            TBBAT = SECOND() - TBBAT
            WRITE(LUPRI,'(20X,A,F10.2,A,/)')
     &      'Total time for this b-batch: ',TBBAT,' seconds'
            CALL FLSHFO(LUPRI)
         ENDIF

         GO TO 100

C     Exit b-batch: calculation done.
C     -------------------------------

  200 CONTINUE

C     Close files for MO Cholesky vectors.
C     ------------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IKEEP,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Print section.
C     --------------

      IF ((IPRINT.GT.0) .OR. LOCDBG) THEN
         IF (LOCDBG) THEN
            XNRM = DSQRT(XNRM)
            CALL HEADER('Debug info from '//SECNAM,-1)
            WRITE(LUPRI,'(5X,A,1P,D30.15,/,5X,A,1P,D30.15,/)')
     &      'Norm of (ai|bj) [full square]: ',XNRM,
     &      'Size of (ai|bj) [db. prec.]  : ',XDIM
         ENDIF
         TIMT  = SECOND() - TIMT
         XMUSE = (D100*MXMEMT)/XLWORK
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &   'Memory reserved for integrals        : ',LWRK,
     &   'Memory reserved for Cholesky vectors : ',LEFT,
     &   'Total memory available               : ',LWORK
         WRITE(LUPRI,'(5X,A,4X,F6.2,A)')
     &   'Maximum memory usage                 : ',XMUSE,' %'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors: ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,I10)')
     &   ' - passes through Cholesky file(s)   : ',NPASS
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sort of MO vectors     : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating (#ai|#bj)  : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating MP2 energy : ',TIME,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used in total                   : ',TIMT,' seconds'
            WRITE(LUPRI,'(5X,A,/)')
     &      '(Debug calculations and printing included.)'
         ELSE
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Time used in total                   : ',TIMT,' seconds'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck get_indvir */
      SUBROUTINE GET_INDVIR(BBEG,BEND,LVIRB,IOFB1,IX1AMB,NX1AMB)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: Set up index arrays for MP2 calculation.
C
#include <implicit.h>
      INTEGER LVIRB(8), IOFB1(8), IX1AMB(8,8), NX1AMB(8)
      INTEGER BBEG, BEND
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>

      CALL IZERO(IOFB1,8)
      CALL IZERO(IX1AMB,64)
      CALL IZERO(NX1AMB,8)

      IB1 = BBEG
      DO ISYMB = ISVI(BBEG),ISVI(BEND)
         IOFB1(ISYMB) = IB1 + NRHFT - IVIR(ISYMB)
         IB1 = IB1 + LVIRB(ISYMB)
      ENDDO

      DO ISYMBJ = 1,NSYM
         ICOUN1 = 0
         DO ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            IX1AMB(ISYMB,ISYMJ) = ICOUN1
            ICOUN1 = ICOUN1 + LVIRB(ISYMB)*NRHF(ISYMJ)
         ENDDO
         NX1AMB(ISYMBJ) = ICOUN1
      ENDDO

      RETURN
      END
C  /* Deck get_btchb */
      SUBROUTINE GET_BTCHB(LVIRB,BBEG,NUMB,LWRK)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: Set up b-batch for MP2CHO_3:
C              (#ai|#bj) with #a estimated.
C
#include <implicit.h>
      INTEGER LVIRB(8)
      INTEGER BBEG
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER LVIRA(8)
      LOGICAL POSSBL(8)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'GET_BTCHB')

C     Check that first b is positive.
C     -------------------------------

      IF (BBEG .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,I10,/)')
     &   'Error in ',SECNAM,':',
     &   'NUMB = ',NUMB
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Check LWRK.
C     -----------

      IF (LWRK .LE. 0) THEN
         NUMB = 0
         RETURN
      ENDIF

      CALL IZERO(LVIRA,NSYM)
      CALL IZERO(LVIRB,NSYM)

      A    = 0
      B    = BBEG - 1
      NUMB = 0

C     Start recursive loop.
C     ---------------------

  100 CONTINUE

         A = A + 1
         B = B + 1
         IF (B .GT. NVIRT) THEN
            B = B - 1
            GO TO 200
         ENDIF

         ISYMA = ISVI(A)
         LVIRA(ISYMA) = LVIRA(ISYMA) + 1

         NUMB  = NUMB + 1
         ISYMB = ISVI(B)
         LVIRB(ISYMB) = LVIRB(ISYMB) + 1

         MEM = MEMAB(LVIRA,LVIRB)
         IF (MEM .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,I10,/,5X,A,/)')
     &      'Error in ',SECNAM,':',
     &      'MEM = ',MEM,
     &      '(The system is probably too large for this calculation)'
            CALL QUIT('Error in '//SECNAM)
         ELSE IF (MEM .LE. LWRK) THEN
            GO TO 100
         ELSE
            B = B - 1
            NUMB  = NUMB - 1
            LVIRB(ISYMB) = LVIRB(ISYMB) - 1
            GO TO 200
         ENDIF

C     Shrink b-batch to fit with at least one
C     a in all possible symmetries.
C     ---------------------------------------

  200 IF (NUMB .GT. 0) THEN
         DO JSYMA = 1,NSYM

            IF (NUMB .LE. 0) RETURN

            CALL IZERO(LVIRA,NSYM)
            LVIRA(JSYMA) = 1

            MEM = MEMAB(LVIRA,LVIRB)
            IF (MEM .GT. LWRK) THEN

  300          CONTINUE

                  IF (NUMB .LE. 0) RETURN

                  JSYMB = ISVI(B)
                  LVIRB(JSYMB) = LVIRB(JSYMB) - 1
                  B = B - 1
                  NUMB = NUMB - 1

                  MEM = MEMAB(LVIRA,LVIRB)
                  IF (MEM .GT. LWRK) GO TO 300

            ENDIF

         ENDDO
      ENDIF

      RETURN
      END
C  /* Deck memab */
      INTEGER FUNCTION MEMAB(LVIRA,LVIRB)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: calculate memory need.
C
#include <implicit.h>
      INTEGER LVIRA(8), LVIRB(8)
#include <ccorb.h>

      MEMAB = 0
      DO ISYMBJ = 1,NSYM
         ISYMAI = ISYMBJ
         ICOUBJ = 0
         ICOUAI = 0
         DO ISYM2 = 1,NSYM
            ISYM1  = MULD2H(ISYM2,ISYMBJ)
            ICOUBJ = ICOUBJ + LVIRB(ISYM1)*NRHF(ISYM2)
            ICOUAI = ICOUAI + LVIRA(ISYM1)*NRHF(ISYM2)
         ENDDO
         MEMAB = MEMAB + ICOUAI*ICOUBJ
      ENDDO

      RETURN
      END
C  /* Deck get_btcha */
      SUBROUTINE GET_BTCHA(LIBJ,NX1AMB,LVIRA,ABEG,NUMA,MEM,LWRK,LWORK)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: Set up a-batch for MP2CHO_3:
C              (#ai|#bj) with #b known.
C
#include <implicit.h>
      INTEGER LIBJ(8), NX1AMB(8), LVIRA(8)
      INTEGER ABEG
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'GET_BTCHA')

C     Check input memory.
C     -------------------

      IF ((LWRK.GT.LWORK) .OR. (LWRK.LT.0) .OR. (LWORK.LT.0)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Inconsistent input in ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'LWRK  = ',LWRK,
     &   'LWORK = ',LWORK
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Initialize.
C     -----------

      CALL IZERO(LVIRA,NSYM)

C     Start recursive loop.
C     ---------------------

      A    = ABEG - 1
      MEM  = 0
      NUMA = 0
  100 CONTINUE

         A = A + 1
         IF (A .GT. NVIRT) THEN
            A = A - 1
            GO TO 200
         ENDIF

         ISYMA = ISVI(A)
         MEM   = MEM  + LIBJ(ISYMA)
         NUMA  = NUMA + 1
         LVIRA(ISYMA) = LVIRA(ISYMA) + 1

         IF (MEM .LE. LWRK) THEN
            GO TO 100
         ELSE
            MEM  = MEM  - LIBJ(ISYMA)
            NUMA = NUMA - 1
            LVIRA(ISYMA) = LVIRA(ISYMA) - 1
            A = A - 1
            GO TO 200
         ENDIF

C     Shrink batch to make sure that at least one
C     Cholesky vector can be stored in memory alongside
C     the (#ai|#bj) integral batch.
C     --------------------------------------------


  200 IF (NUMA .GT. 0) THEN
         DO ISYCHO = 1,NSYM

            IF (NUMA .LE. 0) RETURN

            LEFT = LWORK - MEM
            NEED = MEMCHO(LVIRA,NX1AMB,ISYCHO)
            IF (NEED .GT. LEFT) THEN

  300          CONTINUE

                  IF (NUMA .LE. 0) RETURN

                  JSYMA = ISVI(A)
                  LVIRA(JSYMA) = LVIRA(JSYMA) - 1
                  A    = A - 1
                  NUMA = NUMA - 1
                  MEM  = MEM - LIBJ(JSYMA)

                  LEFT = LWORK - MEM
                  NEED = MEMCHO(LVIRA,NX1AMB,ISYCHO)
                  IF (NEED .GT. LEFT) GO TO 300

            ENDIF

         ENDDO
      ENDIF

      RETURN
      END
C  /* Deck memcho */
      INTEGER FUNCTION MEMCHO(LVIRA,NX1AMB,ISYCHO)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: Calculate memory need for holding one Cholesky
C              vector in MP2CHO_3.
C
#include <implicit.h>
      INTEGER LVIRA(8), NX1AMB(8)
#include <ccorb.h>
#include <ccsdsym.h>

      ICOUNT = 0
      DO ISYMA = 1,NSYM
         ISYMI  = MULD2H(ISYMA,ISYCHO)
         ICOUNT = ICOUNT + LVIRA(ISYMA)*NRHF(ISYMI)
      ENDDO

      MEMCHO = NT1AM(ISYCHO) + ICOUNT + NX1AMB(ISYCHO)

      RETURN
      END
C  /* Deck mp2_memspl2 */
      SUBROUTINE MP2_MEMSPL2(LWORK,NVECTOT,LWRK)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: Split memory for double batching in Cholesky MP2 (and CC2).
C
C     Input:
C        LWORK  : Work space available.
C        NVECTOT: Number of Cholesky vectors in each symmetry.
C
C     Output:
C        LWRK   : Dimension of work space reserved for
C                 (ai|#bj) integral batch.
C
#include <implicit.h>
      INTEGER NVECTOT(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccorb.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chomp2.h>
#include <chocc2.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'MP2_MEMSPL2')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Test LWORK.
C     -----------

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &   'Work space non-positive in ',SECNAM,': LWORK = ',LWORK
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Set Cholesky weight.
C     --------------------

      IF (MP2) THEN
         WCHOL = SPLITM
      ELSE IF (CC2) THEN
         WCHOL = SPLITC
      ELSE
         WCHOL = ONE
      ENDIF

      IF (WCHOL .LE. ZERO) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,1P,D22.15,/)')
     &   'Error in ',SECNAM,':',
     &   'Cholesky weight in memory split, WCHOL = ',WCHOL
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Calculate memory split IDIV according to relative
C     sizes of integrals and Cholesky vectors:
C     (ai|bj) integrals: XT2SQ(1)
C     Cholesky: 2*L(ai,J)
C     -------------------------------------------------

      XCHO = ZERO
      DO ISYCHO = 1,NSYM
         XTST = XT1AM(ISYCHO)*NVECTOT(ISYCHO)
         IF (XCHO .LT. XTST) XCHO = XTST
      ENDDO
      XCHO = TWO*XCHO

      IF ((XCHO.LE.ZERO) .OR. (XT2SQ(1).LE.ZERO)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,':'
         IF (XCHO .LE. ZERO) THEN
            WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &      'Calculated max. dimension of Cholesky part: ',XCHO
            WRITE(LUPRI,'(5X,A,/,8I10)')
     &      'Number of Cholesky vectors supplied:',
     &      (NVECTOT(ISYCHO),ISYCHO=1,NSYM)
         ENDIF
         IF (XT2SQ(1) .LE. ZERO) THEN
            WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &      'Integral dimension from dccsdsym.h: ',XT2SQ(1)
         ENDIF
         WRITE(LUPRI,*)
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      XCHOF = WCHOL*XCHO
      IF (XT2SQ(1) .GT. XCHOF) THEN
         ISPLIT = 1
         XFRAC  = XT2SQ(1)/XCHOF
      ELSE
         ISPLIT = 2
         XFRAC  = XCHOF/XT2SQ(1)
      ENDIF

      IFRAC = INT(XFRAC)
      IDIV  = MAX(IFRAC+1,2)

C     Split memory.
C     -------------

      LWRK = LWORK/IDIV
      IF (ISPLIT .EQ. 1) THEN
         LWRK = LWORK - LWRK
      ENDIF

C     Make sure that we do not reserve too much space for Cholesky vectors.
C     ---------------------------------------------------------------------

      LEFT = LWORK - LWRK
      XLFT = ONE*LEFT
      IF (XLFT .GT. XCHO) THEN
         LEFT = INT(XCHO)
         LWRK = LWORK - LEFT
      ENDIF

C     Make sure that (ai|#bj) integral batch can be held in core
C     for at least 1 virtual.
C     ----------------------------------------------------------

      XWRK  = ONE*LWRK
      XWORK = ONE*LWORK
      DO ISYMB = 1,NSYM
         IF (XWRK .LT. XCKI(ISYMB)) THEN
            XWRK = MIN(XWORK,XCKI(ISYMB))
            LWRK = INT(XWRK)
         ENDIF
      ENDDO

      RETURN
      END
C  /* Deck mp2_memspl3 */
      SUBROUTINE MP2_MEMSPL3(LWORK,NVECTOT,LWRK)
C
C     Thomas Bondo Pedersen, November 2002.
C
C     Purpose: Split memory for triple batching in Cholesky MP2.
C
C     Input:
C        LWORK  : Work space available.
C        NVECTOT: Number of Cholesky vectors in each symmetry.
C
C     Output:
C        LWRK   : Dimension of work space reserved for
C                 (#ai|#bj) integral batch.
C
#include <implicit.h>
      INTEGER NVECTOT(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccorb.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chomp2.h>
#include <chocc2.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'MP2_MEMSPL3')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)

C     Test LWORK.
C     -----------

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &   'Work space non-positive in ',SECNAM,': LWORK = ',LWORK
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Set Cholesky weight.
C     --------------------

      IF (MP2) THEN
         WCHOL = SPLITM
      ELSE IF (CC2) THEN
         WCHOL = SPLITC
      ELSE
         WCHOL = ONE
      ENDIF

      IF (WCHOL .LE. ZERO) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,1P,D22.15,/)')
     &   'Error in ',SECNAM,':',
     &   'Cholesky weight in memory split, WCHOL = ',WCHOL
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Calculate memory split IDIV according to relative
C     sizes of integrals and Cholesky vectors:
C     (ai|bj) integrals: XT2SQ(1)
C     Cholesky: 3*L(ai,J)
C     -------------------------------------------------

      XCHO = ZERO
      DO ISYCHO = 1,NSYM
         XTST = XT1AM(ISYCHO)*NVECTOT(ISYCHO)
         IF (XCHO .LT. XTST) XCHO = XTST
      ENDDO
      XCHO = THREE*XCHO

      IF ((XCHO.LE.ZERO) .OR. (XT2SQ(1).LE.ZERO)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,':'
         IF (XCHO .LE. ZERO) THEN
            WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &      'Calculated max. dimension of Cholesky part: ',XCHO
            WRITE(LUPRI,'(5X,A,/,8I10)')
     &      'Number of Cholesky vectors supplied:',
     &      (NVECTOT(ISYCHO),ISYCHO=1,NSYM)
         ENDIF
         IF (XT2SQ(1) .LE. ZERO) THEN
            WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &      'Integral dimension from dccsdsym.h: ',XT2SQ(1)
         ENDIF
         WRITE(LUPRI,*)
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      XCHOF = WCHOL*XCHO
      IF (XT2SQ(1) .GT. XCHOF) THEN
         ISPLIT = 1
         XFRAC  = XT2SQ(1)/XCHOF
      ELSE
         ISPLIT = 2
         XFRAC  = XCHOF/XT2SQ(1)
      ENDIF

      IFRAC = INT(XFRAC)
      IDIV  = MAX(IFRAC+1,2)

C     Split memory.
C     -------------

      LWRK = LWORK/IDIV
      IF (ISPLIT .EQ. 1) THEN
         LWRK = LWORK - LWRK
      ENDIF

C     Make sure that we do not reserve too much space for Cholesky vectors.
C     ---------------------------------------------------------------------

      LEFT = LWORK - LWRK
      XLFT = ONE*LEFT
      IF (XLFT .GT. XCHO) THEN
         LEFT = INT(XCHO)
         LWRK = LWORK - LEFT
      ENDIF

C     Make sure that (#ai|#bj) integral batch can be held in core
C     for at least 1 virtual.
C     -----------------------------------------------------------

      MAXIJ = 0
      MAXT1 = 0
      DO ISYM = 1,NSYM
         MAXIJ = MAX(MAXIJ,NMATIJ(ISYM))
         MAXT1 = MAX(MAXT1,NT1AM(ISYM))
      ENDDO

      IF (MAXIJ .GT. MAXT1) THEN
         DO ISYMB = 1,NSYM
            DO ISYMA = 1,ISYMB
               ISYMIJ = MULD2H(ISYMA,ISYMB)
               IF (LWRK .LT. NMATIJ(ISYMIJ)) THEN
                  LWRK = MIN(LWORK,NMATIJ(ISYMIJ))
               ENDIF
            ENDDO
         ENDDO
      ELSE
         LEFT = LWORK - LWRK
         DO ISYMI = 1,NSYM
            DO ISYMJ = 1,ISYMI
               NTST = MAXT1 + NRHF(ISYMI) + NRHF(ISYMJ)
               IF (LEFT .LT. NTST) THEN
                  LEFT = MIN(LWORK,NTST)
               ENDIF
            ENDDO
         ENDDO
         LWRK = LWORK - LEFT
      ENDIF

      RETURN
      END
C  /* Deck mp2cho_4 */
      SUBROUTINE MP2CHO_4(EMO,CHMAX,WORK,LWORK,EMP2,MAPAI,MAI,MAPBJ,MBJ,
     &                    THRSCR)
C
C     Henrik Koch and Thomas Bondo Pedersen, February 2003.
C
C     Purpose:
C        Calculate the MP2 energy correction EMP2 directly
C        from Cholesky decomposed two-electron integrals without
C        permanent storage of (ai|bj) integral intermediates.
C        The integrals are set up in batches of the two virtual
C        indices.
C
C        SPARSE VERSION!
C
C     Notes:
C        This routine splits work space in two parts and performs
C        a triple batch; one over b and one over a in (ai|bj) integrals,
C        the third over the MO Cholesky vectors. This implies that, in
C        general, the Cholesky vectors will be read several times.
C
C        This routine requires at least in the order of
C             O*(N + 2)
C        words of memory.
C
C     TODO/FIXME: There are some problems whenever the work space available
C                 is close to the actual minimum required. Might be fixed by
C                 ensuring that the initial b-batch will succeed for all
C                 subsequent a-batches.
C
#include <implicit.h>
      DIMENSION EMO(*), CHMAX(*), WORK(LWORK)
      INTEGER   MAPAI(*), MAI(*), MAPBJ(*), MBJ(*)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccdeco.h>
#include <chomp2.h>
#include <priunit.h>

      INTEGER NVECTOT(8)
      INTEGER LUCHMO(8), LSYM(8)

      INTEGER IX2SQ(8)

C-LBAT is the "batch of batches" length... to fit in cache...
      PARAMETER (LBAT = 500)

      INTEGER BBEG, BEND
      INTEGER LVIRB(8), IOFB1(8), IX1AMB(8,8), NX1AMB(8)
      INTEGER LIBJ(8)
      INTEGER ABEG, AEND
      INTEGER LVIRA(8), IOFA1(8), IX1AMA(8,8), NX1AMA(8)

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.0D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'MP2CHO_4')

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIME = ZERO

C     Initial check of memory.
C     ------------------------

      NEED = 0
      DO ISYMB = 1,NSYM
         DO ISYMA = 1,ISYMB
            ISYMIJ = MULD2H(ISYMA,ISYMB)
            NEEDC  = 0
            DO ISYCHO = 1,NSYM
               ISYMI = MULD2H(ISYCHO,ISYMA)
               ISYMJ = MULD2H(ISYCHO,ISYMB)
               NTST  = NT1AM(ISYCHO) + NRHF(ISYMI) + NRHF(ISYMJ)
               NEEDC = MAX(NEEDC,NTST)
            ENDDO
            NTST = NMATIJ(ISYMIJ) + NEEDC
            NEED = MAX(NEED,NTST)
         ENDDO
      ENDDO

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,A,I10)')
     &   'Work space supplied in ',SECNAM,' is non-positive: ',
     &   'LWORK = ',LWORK
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'Minimum memory required: ',NEED
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: initial'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Minimum memory required: ',NEED,
     &   'Total memory available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initializations.
C     ----------------

      EMP2   = ZERO
      MXMEMT = 0
      MXMEML = 0
      NPASS  = 0
      XOPC   = ZERO
      XLWORK = ONE*LWORK

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

C     Set # vectors.
C     --------------

      IF (CHOMO) THEN
         ITYP = 2
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NCHOMO(ISYM)
         ENDDO
      ELSE
         ITYP = 1
         DO ISYM = 1,NSYM
            NVECTOT(ISYM) = NUMCHO(ISYM)
         ENDDO
      ENDIF

C     Open files for MO Cholesky vectors.
C     -----------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IOPEN,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Debug: initialize integral norm.
C     --------------------------------

      IF (LOCDBG) THEN
         XNRM = ZERO
         XDIM = ZERO
      ENDIF

C     Determine memory split.
C     -----------------------

      CALL MP2_MEMSPL3(LWORK,NVECTOT,LWRK)
      LEFT = LWORK - LWRK

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN
         CALL AROUND('Output from '//SECNAM)
         XMB = XLWORK*WTOMB
         WRITE(LUPRI,'(12X,A,/,12X,A,I10,A,F7.1,A)')
     &   '(#ai|#bj) integrals will be calculated on the fly.',
     &   'Available memory: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(12X,A,1X,I10,A)')
     &   'Maximum memory for integral intermediates: ',LWRK,' words.'
         WRITE(LUPRI,'(12X,A,1P,D12.5)')
     &   'Cholesky screening threshold: ',THRSCR
         CALL FLSHFO(LUPRI)
      ENDIF

C     Start batch loop over b.
C     ------------------------

      IBATB = 0
      BBEG  = 1
      NUMB  = 0
      NPASS = 0
  100 CONTINUE

         BBEG = BBEG + NUMB
         IF (BBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TBBAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATB = IBATB + 1

C        Find out how many b's can be treated,
C        estimating the number of a's along the way.
C        -------------------------------------------

         CALL GET_BTCHB(LVIRB,BBEG,NUMB,LWRK)

C        Check for errors.
C        -----------------

         BEND = BBEG + NUMB - 1
         IF (BEND .GT. NVIRT) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'b-batch number    : ',IBATB
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First B           : ',BBEG,
     &      'Last  B           : ',BEND,
     &      'Number of virtuals: ',NVIRT
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory  : ',LWRK
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMB .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for b-batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get symmetries of first and last b in batch.
C        --------------------------------------------

         ISBBEG = ISVI(BBEG)
         ISBEND = ISVI(BEND)

C        Set up index arrays.
C        --------------------

         CALL GET_INDVIR(BBEG,BEND,LVIRB,IOFB1,IX1AMB,NX1AMB)

         DO ISYIBJ = 1,NSYM
            LIBJ(ISYIBJ) = 0
            DO ISYMBJ = 1,NSYM
               ISYMI = MULD2H(ISYMBJ,ISYIBJ)
               LIBJ(ISYIBJ) = LIBJ(ISYIBJ)
     &                      + NRHF(ISYMI)*NX1AMB(ISYMBJ)
            ENDDO
         ENDDO

C        Print.
C        ------

         IF (LOCDBG) THEN
            WRITE(LUPRI,'(/,A,2I9)')
     &      'First and last b: ',BBEG,BEND
            WRITE(LUPRI,'(A,8I9)')
     &      'IOFB1  :',(IOFB1(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'LVIRB  :',(LVIRB(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NVIR   :',(NVIR(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NX1AMB :',(NX1AMB(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NT1AM  :',(NT1AM(I),I=1,NSYM)
            DO J = 1,NSYM
               WRITE(LUPRI,'(A,8I9)')
     &         'IX1AMB :',(IX1AMB(I,J),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'IT1AM  :',(IT1AM(I,J),I=1,NSYM)
            ENDDO
            WRITE(LUPRI,'(A,8I9)')
     &      'LIBJ   :',(LIBJ(I),I=1,NSYM)
            WRITE(LUPRI,'(A,8I9)')
     &      'NCKI   :',(NCKI(I),I=1,NSYM)
            WRITE(LUPRI,*)
         ENDIF

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(20X,A,I6,A,/,20X,A)')
     &      'b-batch number ',IBATB,':',
     &      '======================'
            WRITE(LUPRI,'(20X,A,I6)')
     &      '#b   : ',NUMB
            WRITE(LUPRI,'(20X,A,I6,A,I1,/,20X,A,I6,A,I1,/)')
     &      'First: ',IOFB1(ISBBEG),'   symmetry: ',ISBBEG,
     &      'Last : ',IOFB1(ISBEND)+LVIRB(ISBEND)-1,'   symmetry: ',
     &                                                  ISBEND
            WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &      '   #a        First           Last         Mem. Usage',
     &      '        Time   ',
     &      '          (a, sym. a)     (a, sym. a)        (%)    ',
     &      '      (seconds)'
            WRITE(LUPRI,'(5X,A,A)')
     &      '----------------------------------------------------',
     &      '---------------'
            CALL FLSHFO(LUPRI)
         ENDIF

C        Start batch loop over a.
C        ------------------------

         IBATA = 0
         ABEG  = 1
         NUMA  = 0
  110    CONTINUE

            ABEG = ABEG + NUMA
            IF (ABEG .GT. NVIRT) GO TO 190

C           Time this batch.
C           ----------------

            TABAT = SECOND()

C           Update batch counter.
C           ---------------------

            IBATA = IBATA + 1

C           Find out how many a's can be treated.
C           -------------------------------------

            CALL GET_BTCHA(LIBJ,NX1AMB,LVIRA,ABEG,NUMA,MEM,LWRK,LWORK)

C           Check for errors.
C           -----------------

            AEND = ABEG + NUMA - 1
            IF (AEND .GT. NVIRT) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Batch error in ',SECNAM,' !'
               WRITE(LUPRI,'(5X,A)')
     &         '- dumping presumably most relevant info:'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'B-batch number    : ',IBATB,
     &         'A-batch number    : ',IBATA
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &         'First B           : ',BBEG,
     &         'Last  B           : ',BEND
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &         'First A           : ',ABEG,
     &         'Last  A           : ',AEND,
     &         'Number of virtuals: ',NVIRT
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Available memory  : ',LWRK,
     &         'Needed    memory  : ',MEM
               CALL QUIT('Batch error in '//SECNAM)
            ENDIF

            IF (NUMA .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory for a-batch in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10)')
     &         'Available memory: ',LWRK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

C           Get symmetries of first and last a in batch.
C           --------------------------------------------

            ISABEG = ISVI(ABEG)
            ISAEND = ISVI(AEND)

C           Set up index arrays.
C           --------------------

            CALL GET_INDVIR(ABEG,AEND,LVIRA,IOFA1,IX1AMA,NX1AMA)

C           Set up (#ai|#bj) index arrays.
C           ------------------------------

            ICOUNT = 0
            DO ISYMBJ = 1,NSYM
               ISYMAI = ISYMBJ
               IX2SQ(ISYMBJ) = ICOUNT
               ICOUNT = ICOUNT + NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)
            ENDDO
            NX2SQ = ICOUNT

            IF (NX2SQ .NE. MEM) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Error detected in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,/)')
     &         'NX2SQ is calculated to be ',NX2SQ,
     &         'MEM determines alloc. as  ',MEM,
     &         'These numbers must be equal....'
               CALL QUIT('Error in '//SECNAM)
            ENDIF

            IF (LOCDBG) THEN
               WRITE(LUPRI,'(/,A,2I9)')
     &         'First and last a: ',ABEG,AEND
               WRITE(LUPRI,'(A,8I9)')
     &         'IOFA1  :',(IOFA1(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'LVIRA  :',(LVIRA(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NVIR   :',(NVIR(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NX1AMA :',(NX1AMA(I),I=1,NSYM)
               WRITE(LUPRI,'(A,8I9)')
     &         'NT1AM  :',(NT1AM(I),I=1,NSYM)
               DO J = 1,NSYM
                  WRITE(LUPRI,'(A,8I9)')
     &            'IX1AMA :',(IX1AMA(I,J),I=1,NSYM)
                  WRITE(LUPRI,'(A,8I9)')
     &            'IT1AM  :',(IT1AM(I,J),I=1,NSYM)
               ENDDO
               WRITE(LUPRI,'(A,8I9)')
     &         'IX2SQ  :',(IX2SQ(I),I=1,NSYM)
               WRITE(LUPRI,'(A,I9)')
     &         'NX2SQ  :',NX2SQ
               WRITE(LUPRI,'(A,I9,/)')
     &         'NT2SQ  :',NT2SQ(1)
            ENDIF

C           Complete allocation for integral intermediates.
C           -----------------------------------------------

            KXINT = 1
            KEND1 = KXINT + NX2SQ
            LWRK1 = LWORK - KEND1 + 1

C           Initialize integral intermediates.
C           ----------------------------------

            CALL DZERO(WORK(KXINT),NX2SQ)

C           Start loop over Cholesky symmetries.
C           ------------------------------------

            ICHO = 0
            DO ISYCHO = 1,NSYM

C              If nothing to do here, cycle.
C              -----------------------------

               IF (NVECTOT(ISYCHO) .LE. 0) GO TO 999
               IF (NX1AMA(ISYCHO)  .LE. 0) GO TO 999
               IF (NX1AMB(ISYCHO)  .LE. 0) GO TO 999

C              Set up mappings.
C              ----------------

               ICOUAI = 0
               ICOUBJ = 0
               DO ISYMI = 1,NSYM

                  ISYMA = MULD2H(ISYMI,ISYCHO)

                  DO I = 1,NRHF(ISYMI)
                     DO LA = 1,LVIRA(ISYMA)
                        ICOUAI = ICOUAI + 1
                        NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1)
     &                      + IOFA1(ISYMA) + LA - 1
                        MAPAI(ICOUAI) = NAI
                     ENDDO
                  ENDDO

                  ISYMJ = ISYMI
                  ISYMB = ISYMA

                  DO J = 1,NRHF(ISYMJ)
                     DO LB = 1,LVIRB(ISYMB)
                        ICOUBJ = ICOUBJ + 1
                        NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                      + IOFB1(ISYMB) + LB - 1
                        MAPBJ(ICOUBJ) = NBJ
                     ENDDO
                  ENDDO

               ENDDO

C              Allocation.
C              -----------

               NBATAI = (NX1AMA(ISYCHO) - 1)/LBAT + 1
               NBATBJ = (NX1AMB(ISYCHO) - 1)/LBAT + 1

               KCHOL = KEND1
               KXAI  = KCHOL + NT1AM(ISYCHO)
               KXBJ  = KXAI  + NX1AMA(ISYCHO)
               KYAI  = KXBJ  + NX1AMB(ISYCHO)
               KYBJ  = KYAI  + LBAT
               KEND2 = KYBJ  + LBAT
               LWRK2 = LWORK - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient memory in '//SECNAM)
               ENDIF

               KLAST  = KEND2 - 1
               MXMEML = MAX(MXMEML,KLAST)

               DO JVEC = 1,NVECTOT(ISYCHO)

                  ICHO = ICHO + 1

                  DTIME = SECOND()
                  CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYCHO),1,JVEC,
     &                            LUCHMO(ISYCHO))
                  DTIME = SECOND() - DTIME
                  TIMR  = TIMR     + DTIME

C                 Copy subblocks.
C                 ---------------

                  CALL CC_CIACPS(WORK(KCHOL),WORK(KXAI),1,IOFA1,LVIRA,
     &                           NX1AMA,IX1AMA,ISYCHO)
                  CALL CC_CIACPS(WORK(KCHOL),WORK(KXBJ),1,IOFB1,LVIRB,
     &                           NX1AMB,IX1AMB,ISYCHO)

C                 Batch of batches.
C                 -----------------

                  ICOUSK = 0
                  DO IBATBJ = 1,NBATBJ

C                    Set up batch info.
C                    ------------------

                     NUMBJ = LBAT
                     IF (IBATBJ .EQ. NBATBJ) THEN
                        NUMBJ = NX1AMB(ISYCHO) - LBAT*(NBATBJ - 1)
                     ENDIF
                     LBJ1 = LBAT*(IBATBJ - 1) + 1

C                    Find max. within batch.
C                    -----------------------

                     XMXBJ = -1.0D8
                     DO LBJ = LBJ1,LBJ1+NUMBJ-1
                        KOFF  = KXBJ + LBJ - 1
                        XMXBJ = MAX(XMXBJ,ABS(WORK(KOFF)))
                     ENDDO

                     DO IBATAI = 1,NBATAI

C                       Set batch info.
C                       ---------------

                        NUMAI = LBAT
                        IF (IBATAI .EQ. NBATAI) THEN
                           NUMAI = NX1AMA(ISYCHO) - LBAT*(NBATAI - 1)
                        ENDIF
                        LAI1 = LBAT*(IBATAI - 1) + 1

C                       Find max. within batch.
C                       -----------------------

                        XMXAI = -1.0D8
                        DO LAI = LAI1,LAI1+NUMAI-1
                           KOFF  = KXAI + LAI - 1
                           XMXAI = MAX(XMXAI,ABS(WORK(KOFF)))
                        ENDDO

C                       Skip if possible.
C                       -----------------

                        IF (XMXAI*XMXBJ .LE. THRSCR) THEN
                           ICOUSK = ICOUSK + 1
                           GO TO 998
                        ENDIF

C                       Screening.
C                       ----------

c                       ICOUN1 = 0
c                       DO LAI = LAI1,LAI1+NUMAI-1
c                          KOFF = KXAI + LAI - 1
c                          XTST = ABS(WORK(KOFF))*XMXBJ
c                          IF (XTST .GT. THRSCR) THEN
c                             WORK(KYAI+ICOUN1) = WORK(KOFF)
c                             ICOUN1 = ICOUN1 + 1
c                             MAI(ICOUN1) = LAI
c                          ENDIF
c                       ENDDO
c                       IF (ICOUN1 .EQ. 0) GO TO 998

c                       ICOUN2 = 0
c                       DO LBJ = LBJ1,LBJ1+NUMBJ-1
c                          KOFF = KXBJ + LBJ - 1
c                          XTST = ABS(WORK(KOFF))*XMXAI
c                          IF (XTST .GT. THRSCR) THEN
c                             WORK(KYBJ+ICOUN2) = WORK(KOFF)
c                             ICOUN2 = ICOUN2 + 1
c                             MBJ(ICOUN2) = LBJ
c                          ENDIF
c                       ENDDO
c                       IF (ICOUN2 .EQ. 0) GO TO 998

c                       DTIME = SECOND()
c                       KOFFI = KXINT + IX2SQ(ISYCHO)
c                       CALL SPRMM(WORK(KOFFI),NX1AMA(ISYCHO),
c    &                             WORK(KYAI),ICOUN1,WORK(KYBJ),ICOUN2,
c    &                             MAI,MBJ)
c                       DTIME = SECOND() - DTIME
c                       TIMC  = TIMC     + DTIME
c                       XOPC  = XOPC + ICOUN1*ICOUN2

                        DO LLAI = 1,NUMAI
                           MAI(LLAI) = LAI1 + LLAI - 1
                        ENDDO

                        DO LLBJ = 1,NUMBJ
                           MBJ(LLBJ) = LBJ1 + LLBJ - 1
                        ENDDO

                        DTIME = SECOND()
                        KOFFI = KXINT + IX2SQ(ISYCHO)
                        KOFAI = KXAI + LAI1 - 1
                        KOFBJ = KXBJ + LBJ1 - 1
                        CALL SPRMM(WORK(KOFFI),NX1AMA(ISYCHO),
     &                             WORK(KOFAI),NUMAI,WORK(KOFBJ),NUMBJ,
     &                             MAI,MBJ)
                        DTIME = SECOND() - DTIME
                        TIMC  = TIMC     + DTIME
                        XOPC  = XOPC + NUMAI*NUMBJ

c                       DTIME = SECOND()
c                       KOFFI = KXINT + IX2SQ(ISYCHO)
c    &                        + NX1AMA(ISYCHO)*(LBJ1 - 1)
c    &                        + LAI1 - 1
c                       KOFAI = KXAI + LAI1 - 1
c                       KOFBJ = KXBJ + LBJ1 - 1
c                       CALL DGEMM('N','T',NUMAI,NUMBJ,1,
c    &                             ONE,WORK(KOFAI),NUMAI,
c    &                             WORK(KOFBJ),NUMBJ,
c    &                             ONE,WORK(KOFFI),NX1AMA(ISYCHO))
c                       DTIME = SECOND() - DTIME
c                       TIMC  = TIMC     + DTIME
c                       XOPC  = XOPC + NUMAI*NUMBJ

  998                   CONTINUE
   
                     ENDDO

                  ENDDO

      write(lupri,*) '***NUMBER OF XMXAI*XMXBJ SKIPS = ',ICOUSK,
     &              ' VEC: ',JVEC,' SYM.: ',ISYCHO
      write(lupri,*) '   NBATAI: ',NBATAI,' NBATBJ: ',NBATBJ,
     &               ' TOTAL: ',NBATAI*NBATBJ
      call flshfo(lupri)
      
               ENDDO

  999          CONTINUE

            ENDDO

C           Debug: Calculate norm of (ai|bj).
C           ---------------------------------

            IF (LOCDBG) THEN
               XNRM = XNRM + DDOT(NX2SQ,WORK(KXINT),1,WORK(KXINT),1)
               XDIM = XDIM + ONE*NX2SQ
            ENDIF

C           Calculate contribution to EMP2.
C           -------------------------------

            DTIME = SECOND()
            DO ISYMBJ = 1,NSYM

               ISYMAI = ISYMBJ

               DO ISYMJ = 1,NSYM

                  ISYMB = MULD2H(ISYMJ,ISYMBJ)
                  IF (LVIRB(ISYMB) .GT. 0) THEN

                     DO J = 1,NRHF(ISYMJ)

                        KOFFJ = IRHF(ISYMJ) + J

                        DO LB = 1,LVIRB(ISYMB)

                           B     = IOFB1(ISYMB) + LB - 1
                           KOFFB = IVIR(ISYMB)  + B

                           LBJ = IX1AMB(ISYMB,ISYMJ)
     &                         + LVIRB(ISYMB)*(J - 1) + LB

                           DO ISYMI = 1,NSYM

                              ISYMA = MULD2H(ISYMI,ISYMAI)
                              IF (LVIRA(ISYMA) .GT. 0) THEN

                                 ISYMAJ = MULD2H(ISYMA,ISYMJ)
                                 ISYMBI = ISYMAJ

                                 DO I = 1,NRHF(ISYMI)

                                    KOFFI = IRHF(ISYMI) + I

                                    LBI = IX1AMB(ISYMB,ISYMI)
     &                                  + LVIRB(ISYMB)*(I - 1) + LB

                                    DO LA = 1,LVIRA(ISYMA)

                                       A     = IOFA1(ISYMA) + LA - 1
                                       KOFFA = IVIR(ISYMA)  + A

                                       LAI = IX1AMA(ISYMA,ISYMI)
     &                                     + LVIRA(ISYMA)*(I - 1) + LA
                                       LAJ = IX1AMA(ISYMA,ISYMJ)
     &                                     + LVIRA(ISYMA)*(J - 1) + LA

                                       NAIBJ = KXINT + IX2SQ(ISYMBJ)
     &                                       + NX1AMA(ISYMAI)*(LBJ - 1)
     &                                       + LAI - 1
                                       NAJBI = KXINT + IX2SQ(ISYMBI)
     &                                       + NX1AMA(ISYMAJ)*(LBI - 1)
     &                                       + LAJ - 1

                                       DENOM = EMO(KOFFA) - EMO(KOFFI)
     &                                       + EMO(KOFFB) - EMO(KOFFJ)
                                       TAIBJ = -WORK(NAIBJ)/DENOM
                                       XAIBJ = TWO*WORK(NAIBJ)
     &                                       - WORK(NAJBI)

                                       EMP2 = EMP2 + TAIBJ*XAIBJ

                                    ENDDO

                                 ENDDO

                              ENDIF

                           ENDDO

                        ENDDO

                     ENDDO

                  ENDIF

               ENDDO

            ENDDO
            DTIME = SECOND() - DTIME
            TIME  = TIME     + DTIME

C           Print info from this a-batch.
C           -----------------------------

            IF (IPRINT .GE. INFO) THEN
               TABAT = SECOND() - TABAT
               XMUSE = (D100*MXMEML)/XLWORK
               LAFST = IOFA1(ISABEG)
               LALST = IOFA1(ISAEND) + LVIRA(ISAEND) - 1
         WRITE(LUPRI,'(5X,I6,4X,I6,1X,I1,8X,I6,1X,I1,9X,F6.2,8X,F10.2)')
     &         NUMA,LAFST,ISABEG,LALST,ISAEND,XMUSE,TABAT
               CALL FLSHFO(LUPRI)
            ENDIF

C           Update memory info.
C           -------------------

            MXMEMT = MAX(MXMEMT,MXMEML)
            MXMEML = 0

C           Update NPASS and go to next a-batch.
C           (Which will exit immediately if no more a's.)
C           ---------------------------------------------

            NPASS = NPASS + 1
            GO TO 110

C        End of #a loop.
C        ---------------

  190    CONTINUE

C        Print info for this b-batch and got to next.
C        (Which will exit immediately if no more b's.)
C        ---------------------------------------------

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(5X,A,A,/)')
     &      '----------------------------------------------------',
     &      '---------------'
            TBBAT = SECOND() - TBBAT
            WRITE(LUPRI,'(20X,A,F10.2,A,/)')
     &      'Total time for this b-batch: ',TBBAT,' seconds'
            CALL FLSHFO(LUPRI)
         ENDIF

         GO TO 100

C     Exit b-batch: calculation done.
C     -------------------------------

  200 CONTINUE

C     Close files for MO Cholesky vectors.
C     ------------------------------------

      DTIME = SECOND()
      CALL CHO_MOP(IKEEP,ITYP,LSYM,LUCHMO,NSYM,1)
      DTIME = SECOND() - DTIME
      TIMR  = TIMR     + DTIME

C     Print section.
C     --------------

      IF ((IPRINT.GT.0) .OR. LOCDBG) THEN
         IF (LOCDBG) THEN
            XNRM = DSQRT(XNRM)
            CALL HEADER('Debug info from '//SECNAM,-1)
            WRITE(LUPRI,'(5X,A,1P,D30.15,/,5X,A,1P,D30.15,/)')
     &      'Norm of (ai|bj) [full square]: ',XNRM,
     &      'Size of (ai|bj) [db. prec.]  : ',XDIM
         ENDIF
         XOPT = ZERO
         DO ISYCHO = 1,NSYM
            XLINT = XT1AM(ISYCHO)*XT1AM(ISYCHO)
            XOPT  = XOPT + XLINT*NVECTOT(ISYCHO)
         ENDDO
         TIMT  = SECOND() - TIMT
         XMUSE = (D100*MXMEMT)/XLWORK
         CALL HEADER('Information from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,1P,D12.5)')
     &   'Screening thr. on Cholesky vectors: ',THRSCR
         WRITE(LUPRI,'(5X,A,1P,D12.5)')
     &   'Actual  operation count           : ',XOPC
         WRITE(LUPRI,'(5X,A,1P,D12.5)')
     &   'Maximum operation count           : ',XOPT
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &   'Memory reserved for integrals        : ',LWRK,
     &   'Memory reserved for Cholesky vectors : ',LEFT,
     &   'Total memory available               : ',LWORK
         WRITE(LUPRI,'(5X,A,4X,F6.2,A)')
     &   'Maximum memory usage                 : ',XMUSE,' %'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O of Cholesky vectors: ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,I10)')
     &   ' - passes through Cholesky file(s)   : ',NPASS
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sort of MO vectors     : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating (#ai|#bj)  : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating MP2 energy : ',TIME,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '---------------------------------------------------------'
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Time used in total                   : ',TIMT,' seconds'
            WRITE(LUPRI,'(5X,A,/)')
     &      '(Debug calculations and printing included.)'
         ELSE
            WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &      'Time used in total                   : ',TIMT,' seconds'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck sprmm */
      SUBROUTINE SPRMM(XINT,NAI,XAI,LAI,XBJ,LBJ,MAI,MBJ)
C
#include <implicit.h>
      DIMENSION XINT(NAI,*), XAI(LAI), XBJ(LBJ)
      INTEGER   MAI(LAI), MBJ(LBJ)

      DO J = 1,LBJ
         DO I = 1,LAI
            XINT(MAI(I),MBJ(J)) = XINT(MAI(I),MBJ(J))
     &                          + XAI(I)*XBJ(J)
         ENDDO
      ENDDO

      RETURN
      END
