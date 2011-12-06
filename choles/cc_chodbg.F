C  /* Deck cc_chodbg */
      SUBROUTINE CC_CHODBG(WORK,LWORK)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose: Drive Cholesky debug calculations.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <chodbg.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHODBG')

C     Print a header.
C     ---------------

      IF (CHODBG) THEN
         WRITE(LUPRI,'(//,20X,A,A,A)')
     &   '>>>>> OUTPUT FROM ',SECNAM,' <<<<<'
      ELSE
         RETURN
      ENDIF

C     Check that this is a Cholesky run.
C     ----------------------------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(/,5X,A,/,5X,A,/,5X,A,A,/)')
     &   'The two-electron integrals *must* be Cholesky decomposed !',
     &   '(i.e. CHOINT flag must be set)',
     &   SECNAM,' returns immediately....'
         RETURN
      ENDIF

C     Check AO integrals.
C     ===================

      IF (DBAOIN) THEN

         CALL AROUND(SECNAM//': Test of AO Integrals')

         CALL DIFFAOIN(WORK,LWORK)

         CALL AROUND(SECNAM//': AO Integrals Tested')

         IF (STAOIN) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      '- ',SECNAM,' stops execution on request...'
            CALL QUIT(' ***** Execution Killed on Request *****')
         ENDIF

      ENDIF

  100 CONTINUE

C     Check (ia|jb) integrals.
C     ========================

      IF (DBIAJB) THEN

         CALL AROUND(SECNAM//': Test of (ia|jb) Integrals')

         KT2AM = 1
         KT1AM = KT2AM + NT2AM(1)
         KEND1 = KT1AM + NT1AM(1)
         LWRK1 = LWORK - KEND1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for (ia|jb) test in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &      'Need (more than): ',KEND1,
     &      'Available       : ',LWORK
            WRITE(LUPRI,'(5X,A,/)')
     &      '- Test skipped, program continues'
            GO TO 200
         ENDIF

C        Calculate Cholesky (ia|jb) integrals.
C        -------------------------------------

         CALL DZERO(WORK(KT1AM),NT1AM(1))

         CHOTIM = SECOND()
         CALL CC_CHOIAJB(WORK(KT2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
         CHOTIM = SECOND() - CHOTIM

C        Calculate (ia|jb) integral-direct and compare.
C        ----------------------------------------------

         CALL DZERO(WORK(KT1AM),NT1AM(1))
         CALL DIFFIAJB(WORK(KT2AM),WORK(KT1AM),WORK(KEND1),LWRK1,
     &                 CHOTIM,.FALSE.)

         CALL AROUND(SECNAM//': (ia|jb) Integrals Tested')

         IF (STIAJB) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      '- ',SECNAM,' stops execution on request...'
            CALL QUIT(' ***** Execution Killed on Request *****')
         ENDIF

      ENDIF

  200 CONTINUE

C     End of Cholesky debug section.
C     ------------------------------

      WRITE(LUPRI,'(//,20X,A,A,A,//)')
     & '>>>>> END OF ',SECNAM,' <<<<<'

C
      RETURN
      END
C  /* Deck diffaoin */
      SUBROUTINE DIFFAOIN(WORK,LWORK)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose: Compare Cholesky and standard integrals.
C
C     Version 2:
C        Calculate all distributions determined by the integral
C        program, batching over distributions as well as Cholesky
C        vectors.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccdeco.h>
#include <chodbg.h>
#include <priunit.h>

      DIMENSION INDEXA(MXCORB)
      DIMENSION ERRNRM(MXCORB), ERRMAX(MXCORB), XINSTD(MXCORB)

      PARAMETER (XMONE = -1.00D0, ZERO = 0.00D0)

      LOGICAL LOCDBG
c     PARAMETER (LOCDBG = .TRUE.)
      PARAMETER (LOCDBG = .false.)

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'DIFFAOIN')

C     Initialize timings.
C     -------------------

      TIMTOT  = SECOND()
      TIMHER1 = ZERO
      TIMHER2 = ZERO
      TIMCHOL = ZERO
      TIMRDAO = ZERO
      TIMCRAO = ZERO
      TIMDIS  = ZERO

      CALL AROUND('START OF '//SECNAM)

C     Check CHOINT flag.
C     ------------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(5X,A,A,A,/,5X,A,/)')
     &   'FATAL ERROR IN ',SECNAM,':',
     &   '- integrals *must* be Cholesky decomposed!'
         CALL QUIT('FATAL ERROR IN '//SECNAM)
      ENDIF

C     Initialize max. abs. error.
C     ---------------------------

      XMXERR = -1.0D10
      SUMCHO = ZERO
      SUMDIR = ZERO

C     Integrals are tot. sym.
C     -----------------------

      ISYMOP = 1

C     Open Cholesky distribution file.
C     --------------------------------

      DTIME   = SECOND()
      CALL CHRDAO(XDUM,IDUM1,IDUM2,0)
      DTIME   = SECOND() - DTIME
      TIMCRAO = TIMCRAO + DTIME

C     Start loop over integral distributions.
C     =======================================

      KEND1 = 1
      LWRK1 = LWORK

      DTIME  = SECOND()
      IF (HERDIR) THEN
         CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
      ELSE
         KCCFB1 = KEND1
         KINDXB = KCCFB1 + MXPRIM*MXCONT
         KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
         LWRK1  = LWORK  - KEND1
         CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &               KODPP1,KODPP2,KRDPP1,KRDPP2,
     &               KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &               WORK(KEND1),LWRK1,IPRERI)
         KEND1 = KFREE
         LWRK1 = LFREE
      ENDIF
      DTIME   = SECOND() - DTIME
      TIMHER1 = TIMHER1  + DTIME

      ICDEL1 = 0
      IF (HERDIR) THEN
         NTOT = MAXSHL
      ELSE
         NTOT = MXCALL
      ENDIF

      DO ILLL = 1,NTOT

         DTIME  = SECOND()
         IF (HERDIR) THEN
            CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                  IPRINT)
         ELSE
            CALL ERIDI2(ILLL,INDEXA,NUMDIS,
     &                  WORK(KODCL1),WORK(KODCL2),
     &                  WORK(KODBC1),WORK(KODBC2),
     &                  WORK(KRDBC1),WORK(KRDBC2),
     &                  WORK(KODPP1),WORK(KODPP2),
     &                  WORK(KRDPP1),WORK(KRDPP2),
     &                  WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1), LWRK1,IPRERI)
         ENDIF
         DTIME   = SECOND() - DTIME
         TIMHER2 = TIMHER2  + DTIME

C        Generate the same distributions from Cholesky vectors.
C        ------------------------------------------------------

         DTIME   = SECOND()
         CALL CCCHAO(WORK(KEND1),LWRK1,INDEXA,NUMDIS)
         DTIME   = SECOND() - DTIME
         TIMCHOL = TIMCHOL  + DTIME

C        Allocation needed for CCRDAO.
C        -----------------------------

         KRECNR = KEND1
         KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: rec.'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         END IF

C        Initialize address on Cholesky distr. file.
C        -------------------------------------------

         IADR = 1

C        Loop over distributions on disk.
C        --------------------------------

         DO IDEL2 = 1,NUMDIS

            IDEL   = INDEXA(IDEL2)
            ISYMD  = ISAO(IDEL)
            IDELT  = IDEL - IBAS(ISYMD)
            ISYDIS = MULD2H(ISYMD,ISYMOP)

            KXDIR = KEND1
            KXCHO = KXDIR + NDISAO(ISYDIS)
            KEND2 = KXCHO + NDISAO(ISYDIS)
            LWRK2 = LWORK - KEND2

            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Insufficient memory in ',SECNAM,
     &         ' - Allocation: Read distr.'
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need     : ',KEND2,
     &         'Available: ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

C           Read standard integral batch for this delta.
C           --------------------------------------------

            DTIME   = SECOND()
            CALL CCRDAO(WORK(KXDIR),IDEL,IDEL2,WORK(KEND2),LWRK2,
     &                  WORK(KRECNR),.TRUE.)
            DTIME   = SECOND() - DTIME
            TIMRDAO = TIMRDAO  + DTIME

C           Read Cholesky integral batch for this delta.
C           --------------------------------------------

            DTIME   = SECOND()
            CALL CHRDAO(WORK(KXCHO),NDISAO(ISYDIS),IADR,2)
            IADR    = IADR + NDISAO(ISYDIS)
            DTIME   = SECOND() - DTIME
            TIMCRAO = TIMCRAO  + DTIME

C           Debug: Calculate norms.
C           -----------------------

c          IF (LOCDBG) THEN
               CHONRM = DSQRT(DDOT(NDISAO(ISYDIS),WORK(KXCHO),1,
     &                                            WORK(KXCHO),1))
               DIRNRM = DSQRT(DDOT(NDISAO(ISYDIS),WORK(KXDIR),1,
     &                                            WORK(KXDIR),1))

c              WRITE(LUPRI,'(A,4I10,/,A,D25.15)')
c    &              'IDEL2,IDEL,IDELT,ISYMD: ',IDEL2,IDEL,IDELT,ISYMD,
c    &              'Norm of Cholesky integrals :',CHONRM 
c              WRITE(LUPRI,'(A,D25.15,/)')
c    &              'Norm of standard integrals :',DIRNRM
c          ENDIF

C
C           Sum test
C
            CALL SUMTST(WORK(KXDIR),WORK(KXCHO),ISYDIS,IDEL,
     &                  SUMCHO,SUMDIR)
C
C           Compare.
C           --------
C
            ialf2 = 0
            do ialf = 1,ndisao(isydis)
               koff1 = kxdir + ialf - 1
               koff2 = kxcho + ialf - 1
               if ((work(koff2) .eq. zero) .and.
     &             (work(koff1) .ne. zero)) ialf2 = ialf2 + 1
            end do
     
               WRITE(LUPRI,'(A,5I6,/,2(A,D20.10))')
     &  'IDEL2,IDEL,IDELT,ISYMD,ialf2: ',IDEL2,IDEL,IDELT,ISYMD,ialf2,
     &              '  Difference norm  :',CHONRM-DIRNRM,
     &              '  Norm of integral :',DIRNRM
C
c           call zercho(work(kxCHO),work(kxDIR),ISYDIS,IABGMX,
c    &                  IABMX,ISABMX,IGMX,ISGMX)
            CALL DAXPY(NDISAO(ISYDIS),XMONE,WORK(KXDIR),1,WORK(KXCHO),1)
            ERRNRM(IDEL) = DSQRT(DDOT(NDISAO(ISYDIS),WORK(KXCHO),1,
     &                                               WORK(KXCHO),1))

            DTIME  = SECOND()
            CALL DISMAX(WORK(KXCHO),ISYDIS,IABG,IAB,ISYMAB,IG,ISYMG,VAL)
            ERRMAX(IDEL) = VAL
            DTIME  = SECOND() - DTIME
            TIMDIS = TIMDIS   + DTIME

            KOFF = KXDIR + IABG - 1
            XINSTD(IDEL) = WORK(KOFF)

            IF (LOCDBG) THEN
               CALL HEADER('Error Report from '//SECNAM,-1)
               WRITE(LUPRI,'(5X,A,I10,4X,I1)')
     &         'Distribution (delta,sym. delta): ',IDELT,ISYMD
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         'Norm of standard distribution  : ',DIRNRM
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         'Norm of Cholesky distribution  : ',CHONRM
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         'Norm of difference distribution: ',ERRNRM(IDEL)
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         'Largest integral error         : ',VAL
               WRITE(LUPRI,'(5X,A,1P,D15.6)')
     &         ' - value of integral           : ',XINSTD(IDEL)
               WRITE(LUPRI,'(5X,A,I10,4X,I1)')
     &         'Location: ab, sym. ab          : ',IAB,ISYMAB
               WRITE(LUPRI,'(5X,A,I10,4X,I1,/)')
     &         'Location: gamma, sym. gamma    : ',IG,ISYMG
            ENDIF
c     CALL QUIT( 'End of test' )

            IF (DABS(VAL) .GT. XMXERR) THEN
               XMXERR = DABS(VAL)
               XMXINT = WORK(KOFF)
               IDELMX = IDEL
               ISDMX  = ISYMD
               IABGMX = IABG
               IABMX  = IAB
               ISABMX = ISYMAB
               IGMX   = IG
               ISGMX  = ISYMG
            ENDIF

         ENDDO

      ENDDO

C     Close and delete Cholesky distribution file.
C     --------------------------------------------

      DTIME   = SECOND()
      CALL CHRDAO(XDUM,IDUM1,IDUM2,4)
      DTIME   = SECOND() - DTIME
      TIMCRAO = TIMCRAO + DTIME

C     Print timings and error summary.
C     --------------------------------

      CALL HEADER('Error Summary from '//SECNAM,-1)
      WRITE(LUPRI,'(5X,A,/,5X,A)')
     & 'Delta Sym.     Error Norm      Max. Error   Std. Integral',
     & '---------------------------------------------------------'
      DO ISYMD = 1,NSYM
         DO D = 1,NBAS(ISYMD)
            IDEL = IBAS(ISYMD) + D
            IF (IDEL .EQ. IDELMX) THEN
               WRITE(LUPRI,1) D,ISYMD,
     &         ERRNRM(IDEL),ERRMAX(IDEL),XINSTD(IDEL)
            ELSE
               WRITE(LUPRI,2) D,ISYMD,
     &         ERRNRM(IDEL),ERRMAX(IDEL),XINSTD(IDEL)
            ENDIF
         ENDDO
      ENDDO
      WRITE(LUPRI,'(5X,A,//)')
     & '---------------------------------------------------------'
C
      WRITE(LUPRI,'(2(/A,D25.15)//)')' SUMDIR : ',SUMDIR,
     &                          ' SUMCHO : ',SUMCHO
C
      TIMTOT = SECOND() - TIMTOT
      CALL HEADER('Timing Report from '//SECNAM,-1)
      IF (HERDIR) THEN
         WRITE(LUPRI,'(20X,A,F10.2,A)')
     &   'Time used in HERDI1: ',TIMHER1,' seconds'
         WRITE(LUPRI,'(20X,A,F10.2,A)')
     &   'Time used in HERDI2: ',TIMHER2,' seconds'
      ELSE
         WRITE(LUPRI,'(20X,A,F10.2,A)')
     &   'Time used in ERIDI1: ',TIMHER1,' seconds'
         WRITE(LUPRI,'(20X,A,F10.2,A)')
     &   'Time used in ERIDI2: ',TIMHER2,' seconds'
      ENDIF
      WRITE(LUPRI,'(20X,A,F10.2,A)')
     & 'Time used in CCCHAO: ',TIMCHOL,' seconds'
      WRITE(LUPRI,'(20X,A,F10.2,A)')
     & 'Time used in CCRDAO: ',TIMRDAO,' seconds'
      WRITE(LUPRI,'(20X,A,F10.2,A)')
     & 'Time used in CHRDAO: ',TIMCRAO,' seconds'
      WRITE(LUPRI,'(20X,A,F10.2,A)')
     & 'Time used in DISMAX: ',TIMDIS,' seconds'
      WRITE(LUPRI,'(20X,A)')
     & '---------------------------------------'
      WRITE(LUPRI,'(20X,A,F10.2,A)')
     & 'Total time for test: ',TIMTOT,' seconds'

      CALL AROUND('END OF '//SECNAM)

      RETURN
    1 FORMAT(5X,I5,2X,I1,2X,1P,D15.6,1X,1P,D15.6,1X,1P,D15.6,
     &       ' < MAX. ERROR')
    2 FORMAT(5X,I5,2X,I1,2X,1P,D15.6,1X,1P,D15.6,1X,1P,D15.6)
      END
C  /* Deck dismax */
      SUBROUTINE DISMAX(XINT,ISYDIS,IABGMX,IABMX,ISABMX,IGMX,ISGMX,VAL)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Find max. abs. value in XINT integral distribution and
C        return its value and location.
C
#include <implicit.h>
      DIMENSION XINT(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER AB, ABG

      PARAMETER (ZERO = 0.00D0)

C     Initialize VAL.
C     ---------------

      VAL = ZERO

C     Find max.
C     ---------

      DO ISYMG = 1,NSYM

         ISYMAB = MULD2H(ISYMG,ISYDIS)

         DO G = 1,NBAS(ISYMG)
            DO AB = 1,NNBST(ISYMAB)

               ABG = IDSAOG(ISYMG,ISYDIS) + NNBST(ISYMAB)*(G - 1) + AB

               IF (DABS(XINT(ABG)) .GT. DABS(VAL)) THEN
                  VAL    = XINT(ABG)
                  IABGMX = ABG
                  IABMX  = AB
                  ISABMX = ISYMAB
                  IGMX   = G
                  ISGMX  = ISYMG
               ENDIF

            ENDDO
         ENDDO

      ENDDO

      RETURN
      END
C  /* Deck sumtst */
      SUBROUTINE SUMTST(XDIR,XCHO,ISYDIS,IDEL,SUMCHO,SUMDIR)
C
C     asm, May 2002.
C
C     Calculate  Sum (ab>=gd) (ab|gd) with Cholesky and standard integrals.
C
#include <implicit.h>
#include <maxorb.h>
      DIMENSION XDIR(*),XCHO(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.00D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C
      ISYMD = ISYDIS
      D = IDEL - IBAS(ISYMD)
C
      DO ISYMG = 1,NSYM
C
         ISYMAB = MULD2H(ISYMG,ISYDIS)
         ISYMGD = ISYMAB
C
         DO G = 1,NBAS(ISYMG)
C
            IF (ISYMG .EQ. ISYMD) THEN
               IF (D .GE. G) THEN
                  NGD = INDEX(G,D)
               ELSE
                  GOTO 100
               END IF
            ELSE IF (ISYMD .GT.ISYMG) THEN
               NGD = NBAS(ISYMG)*(D - 1) + G
            ELSE
               GOTO 100
            END IF
            NGD = IAODPK(ISYMG,ISYMD) + NGD
C
            DO NAB = 1,NNBST(ISYMAB)
C
               NABG = IDSAOG(ISYMG,ISYDIS) 
     &              + NNBST(ISYMAB)*(G - 1) + NAB
C
               SUMDIR = SUMDIR + XDIR(NABG)
               SUMCHO = SUMCHO + XCHO(NABG)
C
            END DO
  100       CONTINUE
         END DO
      END DO

      RETURN
      END
C  /* Deck ccchao */
      SUBROUTINE CCCHAO(WORK,LWORK,INDEXA,NUMDIS)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C
C        Calculate AO integral distributions from Cholesky
C        vectors and write them to disk.
C
C     Input:
C
C        INDEXA : delta indices (not sym. reduced!) of
C                 distributions to be calculated.
C
C        NUMDIS : Number of distributions to be calculated.
C
C     Notes:
C
C        The integral distributions are calculated in a double
C        batch; one over Cholesky vectors, and another over
C        requested distributions.
C
C        Currently, half of the available memory is allocated
C        for each batch. This is probably not optimal, but it
C        is simple.....
C
C        The distribution file *must* be open on calling this
C        routine. Furthermore, data present on that file will
C        be overwritten, and the file will *not* be closed here.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
      INTEGER INDEXA(*)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccisao.h>
#include <ccsdinp.h>
#include <symsq.h>
#include <chodbg2.h>
#include <priunit.h>

      INTEGER IPOINT(MXCORB)
      INTEGER IADR1

      PARAMETER (IDIV = 2, IOPTR = 3)
      PARAMETER (ZERO = 0.00D0, ONE = 1.00D0)

      CHARACTER*6 SECNAM
      PARAMETER (SECNAM = 'CCCHAO')

      LOGICAL LOCDBG
c     PARAMETER (LOCDBG = .TRUE.)
      PARAMETER (LOCDBG = .false.)

C     Statement function INDEX.
C     -------------------------

      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

C     Start timing.
C     -------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIMW = ZERO

C     If nothing to do, return.
C     -------------------------

      IF (NUMDIS .LE. 0) THEN
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(/,5X,A,A,A,I10,/,5X,A,A,/)')
     &      'Info from ',SECNAM,': NUMDIS = ',NUMDIS,
     &      SECNAM,' returns immediately...'
         ENDIF
         RETURN
      ENDIF

C     Only tot. sym. integrals can be treated.
C     ----------------------------------------

      IF (ISYMOP .NE. 1) THEN
         WRITE(LUPRI,'(//,5X,A,A,/,5X,A,I10,/)')
     &   'FATAL ERROR IN ',SECNAM,
     &   'Only tot. sym. integrals may be calculated; ISYMOP = ',ISYMOP
         CALL QUIT('FATAL ERROR IN '//SECNAM)
      ENDIF

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

C     Find largest NBAS.
C     ------------------

      NGMAX = 0
      DO ISYM = 1,NSYM
         NGMAX = MAX(NGMAX,NBAS(ISYM))
      ENDDO

C     Reduce memory for batching.
C     ---------------------------

      LWRK = LWRK0/IDIV

      IF (LWRK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,/,5X,A,I10,/,5X,A)')
     &   'Insufficient memory in ',SECNAM,
     &   'Available: ',LWORK,
     &   '(Increase memory significantly!)'
         WRITE(LUPRI,'(5X,A,I2,A,/)')
     &   'Memory was divided into ',IDIV,' parts'
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Start batch loop over distributions.
C     ------------------------------------

      IDIS1 = 1
      NUMD  = 0
      NPASS = 0
      IADR1 = 1
  100 CONTINUE

         IDIS1 = IDIS1 + NUMD
         IF (IDIS1 .GT. NUMDIS) GOTO 200

C        Figure out how many distributions can be treated.
C        -------------------------------------------------

         CALL GET_DBTCH(IDIS1,NUMD,NUMDIS,MEM,LWRK,INDEXA)

C        Check for errors.
C        -----------------

         IDIS2 = IDIS1 + NUMD - 1
         IF (IDIS2 .GT. NUMDIS) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Batch error in ',SECNAM,' !'
            WRITE(LUPRI,'(5X,A)')
     &      '- dumping presumably most relevant info:'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &      'First distribution: ',IDIS1,
     &      'Last  distribution: ',IDIS2,
     &      'Number of distr.  : ',NUMDIS
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Available memory for batch: ',LWRK,
     &      'Needed    memory for batch: ',MEM
            CALL QUIT('Batch error in '//SECNAM)
         ENDIF

         IF (NUMD .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available memory: ',LWRK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Complete allocation for distributions.
C        --------------------------------------

         KXINT = KEND0
         KEND1 = KXINT + MEM
         LWRK1 = LWORK - KEND1 + 1

C        Initialize integral intermediates.
C        ----------------------------------

         CALL DZERO(WORK(KXINT),MEM)

C        Set up distribution pointer array IPOINT.
C        -----------------------------------------

         ICOUNT = 0
         DO ID = 1,NUMD
            IDIS   = IDIS1 + ID - 1
            IDEL   = INDEXA(IDIS)
            ISYMD  = ISAO(IDEL)
            ISYDIS = MULD2H(ISYMD,ISYMOP)
            IPOINT(ID) = ICOUNT
            ICOUNT     = ICOUNT + NDISAO(ISYDIS)
         ENDDO

C        Start loop over Cholesky symmetries.
C        ------------------------------------

         DO ISYCHO = 1,NSYM

            IF (NUMCHO(ISYCHO) .LE. 0) GO TO 999
            IF (NNBST(ISYCHO)  .LE. 0) GO TO 999

C           Set up Cholesky batch in remaining memory.
C           ------------------------------------------

            LSCR1 = MEMRD(1,ISYCHO,IOPTR)
            LSCR  = MAX(LSCR1,NGMAX)

            MINMEM = NNBST(ISYCHO) + NGMAX
            NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

            IF (NVEC .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory for Cholesky batch in ',SECNAM
               MREQ = NNBST(ISYCHO) + LSCR
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Memory left for Cholesky batch   : ',LWRK1,
     &         'Minimum required (pref. multiple): ',MREQ
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            NVSAV = NVEC

            KCHOL = KEND1
  997       KSCRG = KCHOL + NNBST(ISYCHO)*NVEC
            LEFT  = LWORK - KSCRG + 1
            IF (LSCR1 .GT. LEFT) THEN
               NVEC = NVEC - 1
               IF (NVEC .LE. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A)')
     &            'Insufficient memory for Cholesky batch in ',SECNAM
                  MREQ = NNBST(ISYCHO) + LSCR
              WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/,5X,A,I10)')
     &            'Memory left for Cholesky batch   : ',LWRK1,
     &            'Minimum required (pref. multiple): ',MREQ,
     &            'Original number of vectors       : ',NVSAV,
     &            'Number of distr. held in memory  : ',NUMD
                  WRITE(LUPRI,'(5X,A,/)')
     &        '(Ocurred while shrinking batch to fit remaining memory.)'
                  CALL QUIT('Insufficient memory in '//SECNAM)
               ELSE
                  GO TO 997
               ENDIF
            ENDIF

            NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C           Start batch loop.
C           -----------------

            DO IBATCH = 1,NBATCH

               NUMV = NVEC
               IF (IBATCH .EQ. NBATCH) THEN
                  NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
               ENDIF

               JVEC1 = NVEC*(IBATCH - 1) + 1

               KCHOL = KEND1
               KSCRG = KCHOL + NNBST(ISYCHO)*NUMV
               KEND2 = KSCRG + MAX(LSCR1,NGMAX*NUMV)
               LWRK2 = LWORK - KEND2 + 1

C              Hopefully redundant test.
C              -------------------------

               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,'(//,5X,A,A,A)')
     &            'ALLOCATION BUG DETECTED IN ',SECNAM,':'
                  WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &            'LWORK = ',LWORK,
     &            'KEND2 = ',KEND2
                  CALL QUIT('Allocation bug in '//SECNAM)
               ENDIF

C              Debug: Print.
C              -------------

               IF (LOCDBG) THEN
c              IF (.FALSE.) THEN
                  KLAST = KSCRG + MAX(LSCR1,NGMAX*NUMV) - 1
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Cholesky batch information for pass    ',NPASS+1
                  WRITE(LUPRI,'(5X,A,I10)')
     &            '*************************************************'
                  WRITE(LUPRI,'(5X,A,4X,I10,A,3X,I10)')
     &            'This is batch ',IBATCH,' out of ',NBATCH
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Number of Chol. vec. in this batch :   ',NUMV
                  WRITE(LUPRI,'(5X,A,I6,A,I1)')
     &            ' - beginning with vector  number ',JVEC1,' of sym. ',
     &            ISYCHO
                  JVLST = JVEC1 + NUMV - 1
                  WRITE(LUPRI,'(5X,A,I6,A,I1)')
     &            ' - ending    with vector  number ',JVLST,' of sym. ',
     &            ISYCHO
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Number of distr. in this pass      :   ',NUMD
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'First distribution in this pass    :   ',IDIS1
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Last  distribution in this pass    :   ',IDIS2
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Memory for integral distributions  :   ',MEM
                  WRITE(LUPRI,'(5X,A,I10)')
     &            ' - last location for Chol. vecs.      ',KSCRG-1
                  WRITE(LUPRI,'(5X,A,I10)')
     &            'Last memory location to be used    :   ',KLAST
                  WRITE(LUPRI,'(5X,A,I10,/)')
     &            'Actual memory available            :   ',LWORK
                  IF (KLAST .GT. LWORK) THEN
                     WRITE(LUPRI,'(/,5X,A,A,A,/)')
     &               'BUG DETECTED IN ',SECNAM,' (SEE ABOVE)'
                     CALL QUIT('BUG DETECTED IN '//SECNAM)
                  ENDIF
               ENDIF

C              Read packed Cholesky vectors, incl. zeros.
C              ------------------------------------------

               DTIME = SECOND()
               tnorm = 0.0d0
               erralf = 0.0d0
               difalf = 0.0d0
               DO IVEC = 1,NUMV
                  KOFF = KCHOL + NNBST(ISYCHO)*(IVEC - 1)
                  JVEC = JVEC1 + IVEC - 1
                  CALL CHO_READN(WORK(KOFF),JVEC,1,WORK(KIND1),IDUM2,
     &                           ISYCHO,IOPTR,WORK(KSCRG),LSCR1)
                  xnorm = ddot(nnbst(isycho),work(koff),1,work(koff),1)
                  ynorm = abs(xnorm - vchnrm(jvec,isycho))
c        write(LUPRI,'(2(a,i5),2(a,d25.15))') 'Norm of read vector',
c    &            jvec,' of symmetry',isycho, ' :',
c    &            xnorm,'    difference :',ynorm
                  if (tnorm .lt. ynorm) tnorm = ynorm
c                 if (isycho .eq. 1)  then
c                    kalf = koff + 231
c                    write(LUPRI,'(a,d20.10)') 'Componente 232 :',
c    &                     work(kalf)  
c                 end if
               ENDDO
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME
               write(LUPRI,'(3(a,d15.5))') 'Maxerr in asum',erralf,
     &   ' in difnrm',difalf,' in rdnrm',tnorm

C              Loop over distributions in memory.
C              ----------------------------------

               DO ID = 1,NUMD

                  IDIS   = IDIS1 + ID - 1
                  IDEL   = INDEXA(IDIS)
                  ISYMD  = ISAO(IDEL)
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
                  ISYMG  = MULD2H(ISYDIS,ISYCHO)
                  D      = IDEL - IBAS(ISYMD)

                  IF (NBAS(ISYMG) .LE. 0) GO TO 998

C                 Copy gamma-columns.
C                 -------------------

                  DTIME = SECOND()

                  IF (ISYMG .EQ. ISYMD) THEN

                     DO IVEC = 1,NUMV
                        DO G = 1,NBAS(ISYMG)

                           KOFFC = KCHOL + NNBST(ISYCHO)*(IVEC - 1)
     &                           + IAODPK(ISYMG,ISYMD) + INDEX(G,D) - 1
                           KGD   = KSCRG + NBAS(ISYMG)*(IVEC - 1) + G
     &                           - 1

                           WORK(KGD) = WORK(KOFFC)

                        ENDDO
                     ENDDO

                  ELSE IF (ISYMG .LT. ISYMD) THEN

                     DO IVEC = 1,NUMV

                        KOFFC = KCHOL + NNBST(ISYCHO)*(IVEC - 1)
     &                        + IAODPK(ISYMG,ISYMD)
     &                        + NBAS(ISYMG)*(D - 1)
                        KOFFG = KSCRG + NBAS(ISYMG)*(IVEC - 1)

                        CALL DCOPY(NBAS(ISYMG),WORK(KOFFC),1,
     &                                         WORK(KOFFG),1)

                     ENDDO

                  ELSE

                     DO IVEC = 1,NUMV

                        KOFFC = KCHOL + NNBST(ISYCHO)*(IVEC - 1)
     &                        + IAODPK(ISYMD,ISYMG) + D - 1
                        KOFFG = KSCRG + NBAS(ISYMG)*(IVEC - 1)

                        CALL DCOPY(NBAS(ISYMG),WORK(KOFFC),NBAS(ISYMD),
     &                                         WORK(KOFFG),1)

                     ENDDO

                  ENDIF

                  DTIME = SECOND() - DTIME
                  TIMS  = TIMS     + DTIME

C                 Calculate distribution.
C                 -----------------------

                  DTIME = SECOND()

                  KOFFX = KXINT + IPOINT(ID) + IDSAOG(ISYMG,ISYMD)

                  NAB = NNBST(ISYCHO)
                  NG  = NBAS(ISYMG)

                  CALL DGEMM('N','T',NAB,NG,NUMV,
     &                       ONE,WORK(KCHOL),NAB,WORK(KSCRG),NG,
     &                       ONE,WORK(KOFFX),NAB)

                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME

  998             CONTINUE

               ENDDO

            ENDDO

  999       CONTINUE

         ENDDO

C        Write this batch of distributions to disk.
C        ------------------------------------------

         DTIME = SECOND()
         CALL CHRDAO(WORK(KXINT),MEM,IADR1,1)
         DTIME = SECOND() - DTIME
         TIMW  = TIMW     + DTIME

C        Debug: print norms.
C        -------------------

         IF (LOCDBG) THEN
            DO ID = 1,NUMD
               IDIS   = IDIS1 + ID - 1
               IDEL   = INDEXA(IDIS)
               ISYMD  = ISAO(IDEL)
               ISYDIS = MULD2H(ISYMD,ISYMOP)
               D      = IDEL - IBAS(ISYMD) 
               KOFFX  = KXINT + IPOINT(ID)
               DNORM  = DSQRT(DDOT(NDISAO(ISYDIS),WORK(KOFFX),1,
     &                                            WORK(KOFFX),1))
               WRITE(LUPRI,'(/,5X,A,A,3I10,1X,I1)')
     &         SECNAM,': IDIS, IDEL, D, ISYMD: ',IDIS,IDEL,D,ISYMD
               WRITE(LUPRI,'(5X,A,A,1P,D15.6,/)')
     &         SECNAM,': Distr. norm         : ',DNORM
            ENDDO
         ENDIF

C        Update info and go to next distribution batch.
C        (which will exit immediately, if no more delta's)
C        -------------------------------------------------

         IADR1 = IADR1 + MEM
         NPASS = NPASS + 1
         GO TO 100

C     Exit point for distribution batch.
C     ----------------------------------

  200 CONTINUE

C     Print info.
C     -----------

      IF (IPRINT .GT. 0) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,2X,A,A,I3,A,I3,A)')
     &   SECNAM,' generated ',NUMDIS,' distributions in ',
     &   NPASS,' passes through Cholesky files.'
         WRITE(LUPRI,'(2X,A,A,A,F10.2,A,/)')
     &   'Total time used in ',SECNAM,': ',TIMT,' seconds'
      ENDIF

      IF (LOCDBG) THEN
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Cholesky I/O     : ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Cholesky sorting : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for final contraction: ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used for distribution I/O : ',TIMW,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck get_dbtch */
      SUBROUTINE GET_DBTCH(IDIS1,NUMD,NUMDIS,MEM,LWRK,INDEXA)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Set up batch # 1 in CCCHAO.
C
#include <implicit.h>
      INTEGER INDEXA(*)
#include <maxorb.h>
#include <ccisao.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      IDIS = IDIS1 - 1
      MEM  = 0
      NUMD = 0
  100 CONTINUE

         IDIS = IDIS + 1
         IF (IDIS .GT. NUMDIS) GO TO 200

         IDEL   = INDEXA(IDIS)
         ISYMD  = ISAO(IDEL)
         ISYDIS = MULD2H(ISYMD,ISYMOP)
         MEM    = MEM  + NDISAO(ISYDIS)
         NUMD   = NUMD + 1

         IF (MEM .LE. LWRK) THEN
            GO TO 100
         ELSE
            MEM  = MEM  - NDISAO(ISYDIS)
            NUMD = NUMD - 1
            GO TO 200
         ENDIF

  200 RETURN
      END
C  /* Deck chrdao */
      SUBROUTINE CHRDAO(XINT,LEN,IADR1,IOPT)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        I/O handling of integral distributions generated from
C        Cholesky vectors. Uses crayio routines.
C
C     IOPT = 0 : Open file.
C
C     IOPT = 1 : Write to file, starting at IADR1.
C
C     IOPT = 2 : Read from file, starting at IADR1.
C
C     IOPT = 3 : Close file, keep it.
C
C     IOPT = 4 : Close file, delete it.
C
#include <implicit.h>
      DIMENSION XINT(*)
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccisao.h>
#include <priunit.h>

      CHARACTER*7 FCHODIS
      PARAMETER (LUCHODIS = 98, FCHODIS = 'CHODIST')

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHRDAO')

      INTEGER IADR1

      IF (IOPT .EQ. 0) THEN

         CALL WOPEN2(LUCHODIS,FCHODIS,64,0)

      ELSE IF (IOPT .EQ. 1) THEN

         IF (LEN .EQ. 0) RETURN
         CALL PUTWA2(LUCHODIS,FCHODIS,XINT,IADR1,LEN)

      ELSE IF (IOPT .EQ. 2) THEN

         IF (LEN .EQ. 0) RETURN
         CALL GETWA2(LUCHODIS,FCHODIS,XINT,IADR1,LEN)

      ELSE IF (IOPT .EQ. 3) THEN

         CALL WCLOSE2(LUCHODIS,FCHODIS,'KEEP')

      ELSE IF (IOPT .EQ. 4) THEN

         CALL WCLOSE2(LUCHODIS,FCHODIS,'DELETE')

      ELSE

         WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &   'FATAL ERROR IN ',SECNAM,': IOPT = ',IOPT
         CALL QUIT('FATAL ERROR IN '//SECNAM)

      ENDIF

      RETURN
      END
C
C  /* Deck zercho */
      SUBROUTINE ZERCHO(CHOINT,DIRINT,ISYDIS,IABGMX,IABMX,ISABMX,
     &                  IGMX,ISGMX)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Find max. abs. value in XINT integral distribution and
C        return its value and location.
C
#include <implicit.h>
      DIMENSION CHOINT(*),DIRINT(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER AB, ABG

      PARAMETER (ZERO = 0.00D0)

C     Find max.
C     ---------

      DO ISYMG = 1,NSYM

         ISYMAB = MULD2H(ISYMG,ISYDIS)

         DO G = 1,NBAS(ISYMG)
            DO AB = 1,NNBST(ISYMAB)

               ABG = IDSAOG(ISYMG,ISYDIS) + NNBST(ISYMAB)*(G - 1) + AB

               XCHO = CHOINT(ABG)
               XDIR = DIRINT(ABG)
               IF ((XCHO .EQ. ZERO) .AND. (XDIR .NE. ZERO)) THEN
                  IABGMX = ABG
                  IABMX  = AB
                  ISABMX = ISYMAB
                  IGMX   = G
                  ISGMX  = ISYMG
                  WRITE(LUPRI,'(A,5I10)') 
     &  'Zero choint, isymab,isymg,abg,ab,g',isymab,isymg,abg,ab,g
               ENDIF

            ENDDO
         ENDDO

      ENDDO

c
c     CALL QUIT( 'End of test' )
      RETURN
      END

