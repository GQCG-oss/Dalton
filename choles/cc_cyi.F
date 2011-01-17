C  /* Deck cc_cyi */
      SUBROUTINE CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &                  UVEC,ISYMU,VVEC,ISYMV,FAC2,NUMUV,
     &                  SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &                  JTYPZ,ISYMZ,NCHOZ,NUMZ,JTYPY,
     &                  WVEC,ISYMW,NUMW,SVEC,
     &                  WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2)
C
C     Coupled Cluster Cholesky Y Intermediates.
C     =========================================
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates according to
C
C        Y(ai,J) = sum(bj) X(ai,bj) * Z(bj,J)
C
C     where the amplitudes are generated in batches over virtual
C     indices from the general formula:
C
C        X(ai,bj) = (2-Pij) 1/[FREQ - e(ai,bj)] P(ai,bj)
C                 * [ FAC1 * sum(I=1,NUMP12) 
C                   * sum(k=1,NSYM) sum(J=1,NCHP12(k,I)) L1(ai,J,I) * L2(bj,J,I)
C                   + FAC2 * sum(I=1,NUMUV) UVEC(ai,I) * VVEC(bj,I)
C                   ]
C
C     (2-Pij) generates 2 Coulomb minus exchange, e(ai,bj) is the 
C     usual orbital energy denominator, and P(ai,bj) symmetrizes in
C     the ai,bj indices.
C
C     Additionally, the routine may update SVEC according to
C
C        SVEC(ai,I) = SVEC(ai,I) + sum(bj) X(ai,bj) * WMAT(bj,I)
C
C     for I = 1,NUMW.
C
C     Input:
C     ======
C
C     - specification of amplitude assembly:
C     JTYP1 : The "type" of the L1 Cholesky vectors used for file
C             opening in CHO_MOP, dimension: NUMP12.
C     ISYMP1: Symmetries of the L1-vectors, dimension: NUMP12.
C     JTYP2 : The "type" of the L2 Cholesky vectors used for file
C             opening in CHO_MOP, dimension: NUMP12.
C     ISYMP2: Symmetries of the L2-vectors, dimension: NUMP12.
C     NCHP12: # number of L1/L2-vectors, dimension: 8,NUMP12.
C     FAC1  : Scaling factor for L1/L2 part.
C     UVEC  : Array containing the U-vectors.
C     ISYMU : Symmetries of the U-vectors, dimension: NUMUV.
C     VVEC  : Array containing the V-vectors.
C     ISYMV : Symmetries of the V-vectors, dimension: NUMUV.
C     FAC2  : Scaling factor for U/V part.
C     SCD   : Diagonal scaling factor for each Cholesky integral assembly,
C             dimension: NUMP12.
C     SCDG  : Diagonal scaling factor of the final amplitudes.
C     IOPTDN: Include denominators (denominators are not needed if
C             the amplitudes have been separately decomposed).
C             IOPTDN = 1: include, else: exclude.
C     FOCKD : Canonical orbital energies (used only if IOPTDN = 1).
C     FREQ  : Frequency (used only if IOPTDN = 1).
C     IOPTCE: Set up 2 Coulomb minus exchange.
C             IOPTCE = 1: set up 2CME, else: don't.
C
C     - specification of Y intermediate calculation:
C     JTYPZ : The "type" of the Z Cholesky vectors used for file
C             opening in CHO_MOP, dimension: NUMZ.
C     ISYMZ : Symmetries of the Z-vectors, dimension: NUMZ.
C     NCHOZ :  # number of Y/Z-vectors, dimension: 8,NUMZ.
C     JTYPY : The "type" of the Z Cholesky vectors used for file
C             opening/read/write in CC_CYI i/o routines, dimension: NUMZ.
C
C     - specification of SVEC calculation:
C     WVEC  : Array containing the W vectors.
C     ISYMW : Symmetries of the W vectors, dimension: NUMW.
C
C     - work space:
C     WORK  : dimension: LWORK.
C
C     - file deleting at end of routine:
C     DELP1 : flags for deleting L1-vector files, dimension: NUMP12.
C     DELP2 : flags for deleting L2-vector files, dimension: NUMP12.
C
C
C     Output:
C     -------
C
C     The Y intermediates are returned on disk, as specified through JTYPY.
C
C     SVEC  : Array containing the S vectors (if NUMW > 0).
C     X2NRM : Norm of packed amplitude array, before 2CME.
C     X2CNM : Norm of packed amplitude array, after  2CME (if IOPTCE = 1).
C
C     Notes:
C     ------
C
C     - The Cholesky vectors must be available on disk, through CHO_MOP.
C     - NUMP12 <= 0 is treated as an error, i.e. the U*V part cannot be
C       done alone (yet). The problem is the structure of the memory split.
C     - There are several simplified interfaces to this routine, e.g.
C       CC_CYIV0 used for setting up the ground state vector function.
C
#include <implicit.h>
      DIMENSION UVEC(*), VVEC(*), WVEC(*), SVEC(*), FOCKD(*)
      DIMENSION SCD(NUMP12), WORK(LWORK)
      INTEGER JTYP1(NUMP12), ISYMP1(NUMP12)
      INTEGER JTYP2(NUMP12), ISYMP2(NUMP12)
      INTEGER NCHP12(8,NUMP12)
      INTEGER ISYMU(NUMUV), ISYMV(NUMUV)
      INTEGER JTYPZ(NUMZ), ISYMZ(NUMZ), NCHOZ(8,NUMZ), JTYPY(NUMZ)
      INTEGER ISYMW(NUMW)
      LOGICAL DELP1(NUMP12), DELP2(NUMP12)
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <ciarc.h>
#include <chocc2.h>
#include <priunit.h>
#include <cyit.h>     
#include <dummy.h>

      INTEGER ABEG, AEND, BBEG, BEND
      INTEGER LIBJ(8)

      LOGICAL LOCDBG, FULLIN
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*6 SECNAM
      PARAMETER (SECNAM = 'CC_CYI')

      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IDEL = 0, IKEEP = 1)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, HALF = 0.5D0)
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.0D0)
      PARAMETER (SMALL = 1.0D-4, TINY = 1.0D-14)

C     Initialize timings.
C     -------------------

      TIMCYI(1) = SECOND()
      DO ITMCYI = 2,NTMCYI
         TIMCYI(ITMCYI) = ZERO
      ENDDO
      DO ITCYIO = 1,NTCYIO
         TMCYIO(ITCYIO) = ZERO
      ENDDO
      DO ITCALC = 1,NTCALC
         TMCALC(ITCALC) = ZERO
      ENDDO

C     Check consistency of input.
C     ---------------------------

      CALL CC_CYINP1(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &               UVEC,ISYMU,VVEC,ISYMV,FAC2,NUMUV,
     &               SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &               JTYPZ,ISYMZ,NCHOZ,NUMZ,JTYPY,
     &               WVEC,ISYMW,NUMW,SVEC,
     &               WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2,IERR,LOCDBG,
     &               ISYINT)
      IF (IERR .NE. 0) CALL QUIT('Inconsistent input in '//SECNAM)

C     Issue warnings for small scaling factors, but do not stop.
C     ----------------------------------------------------------

      IF (NUMP12 .GT. 0) THEN
         IF (DABS(FAC1) .LE. SMALL) THEN
            WRITE(LUPRI,'(/,5X,A,A,/,5X,A,1P,D22.15)')
     &      '*** WARNING: Small scaling factor in ',SECNAM,
     &      '             FAC1 = ',FAC1
         ENDIF
      ENDIF
      IF (NUMUV .GT. 0) THEN
         IF (DABS(FAC2) .LE. SMALL) THEN
            WRITE(LUPRI,'(/,5X,A,A,/,5X,A,1P,D22.15)')
     &      '*** WARNING: Small scaling factor in ',SECNAM,
     &      '             FAC2 = ',FAC2
         ENDIF
      ENDIF

C     Global (static) configuration parameters of CC_CIA.
C     ---------------------------------------------------

      CIAMIO = .FALSE.
      GETMNM = .FALSE.

      IPRCIA = IPRINT - 30 + IPRLVL

C     Get the largest minimum memory need in CC_CIA (incl. integrals).
C     ----------------------------------------------------------------

      GETMNM = .TRUE.
      NEED   = 0
      IMAX   = 0
      DO I = 1,NUMP12

         ITYP1  = JTYP1(I)
         ITYP2  = JTYP2(I)
         ISYCH1 = ISYMP1(I)
         ISYCH2 = ISYMP2(I)
         DO ISYCHO = 1,NSYM
            NTOVEC(ISYCHO) = NCHP12(ISYCHO,I)
         ENDDO

         IPRCIA = IPRINT - 30 + IPRLVL

         LIDEN  = ITYP1 .EQ. ITYP2
         SYMTRZ = ITYP1 .NE. ITYP2

         CALL CC_CIA(DUMMY,DUMMY,IDUMMY,DUMMY,NEEDI,IDUMMY)

         IF (NEEDI .GT. NEED) THEN
            NEED = NEEDI
            IMAX = I
         ENDIF

      ENDDO
      GETMNM = .FALSE.

C     Initial check of memory.
C     ------------------------

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': Work space non-positive!'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &   'Min. memory required: ',NEED,
     &   'Work space supplied : ',LWORK
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
      MXMEML = 0
      NPASS  = 0
      NCALL  = 0
      NCALLY = 0
      XLWORK = ONE*LWORK
      X2NRM  = ZERO
      X2CNM  = ZERO

C     Determine memory split based on the assembly with largest
C     minimum need (IMAX).
C     ---------------------------------------------------------

      IF (IMAX .GT. 0) THEN

         ISYCH1 = ISYMP1(IMAX)
         ISYCH2 = ISYMP2(IMAX)
         ITYP1  = JTYP1(IMAX)
         ITYP2  = JTYP2(IMAX)
         DO ISYCHO = 1,NSYM
            NTOVEC(ISYCHO) = NCHP12(ISYCHO,IMAX)
         ENDDO

         LIDEN  = ITYP1 .EQ. ITYP2
         SYMTRZ = ITYP1 .NE. ITYP2
         CIAMIO = .FALSE.
         GETMNM = .FALSE.
         INXINT = IMAX .EQ. 1

         CALL CC_CYIMSP(LWORK,LWRK)
         LEFT = LWORK - LWRK

      ELSE

         WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &   'Error in ',SECNAM,': IMAX is non-positive: ',IMAX
         CALL QUIT('IMAX error in '//SECNAM)

      ENDIF

C     Check.
C     ------

      IF ((LWRK.LE.0) .OR. (LEFT.LE.0)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &   'Total memory available                  : ',LWORK,
     &   'Memory reserved for integrals/amplitudes: ',LWRK,
     &   'Memory reserved for Cholesky vectors    : ',LEFT
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN

         CALL AROUND('Output from '//SECNAM)

         WRITE(LUPRI,'(2X,A,I1,1X,1P,D22.15)')
     &   'Symmetry and dimension of integrals/amplitudes: ',
     &   ISYINT,XT2SQ(ISYINT)

         XMB  = XLWORK*WTOMB
         XPC1 = (D100*LWRK)/XLWORK
         XPC2 = (D100*LEFT)/XLWORK
         WRITE(LUPRI,'(2X,A,I10,A,F7.1,A)')
     &   'Total memory available: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(2X,A,I10,A,F7.1,A,/,2X,A,I10,A,F7.1,A,/)')
     &   'Maximum mem. allowed for integrals/amplitudes: ',LWRK,
     &   ' words (',XPC1,'%).',
     &   'Remaining memory used for Cholesky vectors   : ',LEFT,
     &   ' words (',XPC2,'%).'

         CALL FLSHFO(LUPRI)

      ENDIF

C     Start batch loop over b.
C     ------------------------

      IBATB = 0
      BBEG  = 1
      NUMB  = 0
  100 CONTINUE

         BBEG = BBEG + NUMB
         IF (BBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TBBAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATB = IBATB + 1

C        Initialize min. and max. Cholesky batch variables.
C        --------------------------------------------------

         MNCHBT = 230000000
         MXCHBT = -230000000

C        Find out how many b's can be treated, estimating
C        the number of a's along the way.
C        ------------------------------------------------

         CALL CC_CYIBTB(LVIRB,BBEG,NUMB,LWRK,ISYINT)

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

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(5X,A,I6,A,/,5X,A)')
     &      'b-batch number ',IBATB,':',
     &      '======================'
            WRITE(LUPRI,'(5X,A,I6)')
     &      '#b   : ',NUMB
            WRITE(LUPRI,'(5X,A,I6,A,I1,/,5X,A,I6,A,I1,/)')
     &      'First: ',IOFB1(ISBBEG),'   symmetry: ',ISBBEG,
     &      'Last : ',IOFB1(ISBEND)+LVIRB(ISBEND)-1,'   symmetry: ',
     &                                                  ISBEND
            WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &      '    #a   First     Last    Memory             CPU Time',
     &      ' (Seconds)',
     &      '       (a,sym.a) (a,sym.a) Usage (%)',
     &      '  Integrals  Y-interm.      Total'
            WRITE(LUPRI,'(5X,A,A)')
     &      '------------------------------------------------------',
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

            CALL CC_CYIBTA(LIBJ,NX1AMB,LVIRA,ABEG,NUMA,MEM,LWRK,ISYINT)

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

            CALL IZERO(IX2SQ,64)
            ICOUNT = 0
            DO ISYMBJ = 1,NSYM
               ISYMAI = MULD2H(ISYMBJ,ISYINT)
               IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
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

C           Allocation for integral intermediates.
C           --------------------------------------

            KXINT = 1
            KEND1 = KXINT + NX2SQ
            LWRK1 = LWORK - KEND1 + 1

C           Calculate integrals from Cholesky vectors.
C           ------------------------------------------

            TIMI = ZERO
            DO IP12 = 1,NUMP12

               ITYP1  = JTYP1(IP12)
               ITYP2  = JTYP2(IP12)
               ISYCH1 = ISYMP1(IP12)
               ISYCH2 = ISYMP2(IP12)
               DO ISYCHO = 1,NSYM
                  NTOVEC(ISYCHO) = NCHP12(ISYCHO,IP12)
               ENDDO

               LIDEN  = ITYP1 .EQ. ITYP2
               SYMTRZ = ITYP1 .NE. ITYP2
               INXINT = IP12  .EQ. 1

               SCDLOC = SCD(IP12)

               DTIME = SECOND()
               CALL CC_CIA(WORK(KXINT),WORK(KEND1),LWRK1,SCDLOC,KLAST,
     &                     NPASS)
               DTIME = SECOND() - DTIME
               TIMI  = TIMI     + DTIME

               NCALL  = NCALL + 1
               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)

               DO ISYCHO = 1,NSYM
                  MNCHBT = MIN(MNCHBT,NCHBAT(ISYCHO))
                  MXCHBT = MAX(MXCHBT,NCHBAT(ISYCHO))
               ENDDO

               IF (LOCDBG) THEN
                  IF (NX2SQ .EQ. NT2SQ(ISYINT)) THEN
                     CALL CC_CYICHKSYM(WORK(KXINT),ISYINT,ERRMAX)
                     WRITE(LUPRI,'(A,I3,A,1P,D22.15)')
     &               'Max. sym. error after call',NCALL,
     &               ' to CC_CIA: ',ERRMAX
                  ENDIF
               ENDIF

            ENDDO
            TIMCYI(2) = TIMCYI(2) + TIMI
            TMCALC(1) = TMCALC(1) + TIMCIA(4)

C           Scale integrals.
C           ----------------

            IF (FAC1 .NE. ONE) CALL DSCAL(NX2SQ,FAC1,WORK(KXINT),1)

C           Add FAC2*U*V terms (if any).
C           ----------------------------

            KOFFU = 1
            KOFFV = 1
            DO IUV = 1,NUMUV

               ISYMAI = ISYMU(IUV)
               ISYMBJ = ISYMV(IUV)

               DTIME = SECOND()
               CALL CC_CYILHX(UVEC(KOFFU),WORK(KXINT),VVEC(KOFFV),
     &                        WORK(KEND1),LWRK1,KLAST,ISYMAI,ISYMBJ,
     &                        FAC2)
               DTIME = SECOND() - DTIME
               TIMCYI(3) = TIMCYI(3) + DTIME

               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)

               KOFFU = KOFFU + NT1AM(ISYMAI)
               KOFFV = KOFFV + NT1AM(ISYMBJ)

               IF (LOCDBG) THEN
                  IF (NX2SQ .EQ. NT2SQ(ISYINT)) THEN
                     CALL CC_CYICHKSYM(WORK(KXINT),ISYINT,ERRMAX)
                     WRITE(LUPRI,'(A,I4,A,1P,D22.15)')
     &               'Max. sym. error after call',IUV,' to CC_CYILHX: ',
     &               ERRMAX
                  ENDIF
               ENDIF

            ENDDO

C           Scale diagonal (if needed).
C           ---------------------------

            IF ((SCDG.NE.ONE) .AND. (ISYINT.EQ.1)) THEN
               XTST   = ONE*NX2SQ
               DIFF   = DABS(XT2SQ(1) - XTST)
               FULLIN = DIFF .LE. TINY
               CALL CC_CIADSCL(WORK(KXINT),ISYINT,SCDG,FULLIN)
            ENDIF

C           Divide by orbital energies to obtain amplitudes,
C           if requested through IOPTDN.
C           ------------------------------------------------

            IF (IOPTDN .EQ. 1) THEN
               DTIME = SECOND()
               CALL CC_DNOM(FOCKD,WORK(KXINT),FREQ,ISYINT)
               DTIME = SECOND() - DTIME
               TIMCYI(4) = TIMCYI(4) + DTIME
            ENDIF

C           Calculate contribution to (packed) doubles norm.
C           ------------------------------------------------

            DTIME = SECOND()
            CALL CC_CYINRM(WORK(KXINT),ISYINT,X2NRM)
            DTIME = SECOND() - DTIME
            TIMCYI(5) = TIMCYI(5) + DTIME


C           Set up 2 Coulomb minus exchange (almost) in place, if
C           requested through IOPTCE, and compute contribution to
C           packed 2CME norm.
C           -----------------------------------------------------

            IF (IOPTCE .EQ. 1) THEN

               DTIME = SECOND()
               CALL CC_CYITCME(WORK(KXINT),WORK(KEND1),LWRK1,ISYINT,
     &                         KLAST)
               DTIME = SECOND() - DTIME
               TIMCYI(6) = TIMCYI(6) + DTIME

               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)

               DTIME = SECOND()
               CALL CC_CYINRM(WORK(KXINT),ISYINT,X2CNM)
               DTIME = SECOND() - DTIME
               TIMCYI(5) = TIMCYI(5) + DTIME

            ENDIF

C           Calculate SVECs (if any).
C           -------------------------

            KOFFS = 1
            KOFFW = 1
            DO IW = 1,NUMW

               ISYMBJ = ISYMW(IW)
               ISYMAI = MULD2H(ISYMBJ,ISYINT)

               DTIME = SECOND()
               CALL CC_CYICNI(SVEC(KOFFS),WVEC(KOFFW),WORK(KXINT),
     &                        WORK(KEND1),LWRK1,ISYMBJ,ISYINT,KLAST)
               DTIME = SECOND() - DTIME
               TIMCYI(7) = TIMCYI(7) + DTIME

               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)

               KOFFS = KOFFS + NT1AM(ISYMAI)
               KOFFW = KOFFW + NT1AM(ISYMBJ)

            ENDDO

C           Calculate Y intermediates.
C           --------------------------

            NCALLY = NCALLY + 1
            TIMY   = ZERO
            DO IZ = 1,NUMZ

               DTIME = SECOND()
               CALL CC_CYIM(WORK(KXINT),WORK(KEND1),LWRK1,ISYINT,
     &                      JTYPZ(IZ),ISYMZ(IZ),JTYPY(IZ),KLAST,
     &                      NCALLY,NCHOZ(1,IZ))
               DTIME = SECOND() - DTIME
               TIMY  = TIMY + DTIME

               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)

            ENDDO
            TIMCYI(8) = TIMCYI(8) + TIMY

C           Print information from this a batch.
C           ------------------------------------

            IF (IPRINT .GE. INFO) THEN
               TABAT = SECOND() - TABAT
               XMUSE = (D100*MXMEML)/XLWORK
               LAFST = IOFA1(ISABEG)
               LALST = IOFA1(ISAEND) + LVIRA(ISAEND) - 1
               WRITE(LUPRI,1) NUMA,LAFST,ISABEG,LALST,ISAEND,XMUSE,
     &                        TIMI,TIMY,TABAT
    1          FORMAT(5X,I6,1X,I6,1X,I1,2X,I6,1X,I1,2X,F7.2,3X,F10.2,1X,
     &                F10.2,1X,F10.2)
               CALL FLSHFO(LUPRI)
            ENDIF

C           Update memory info.
C           -------------------

            MXMEMT = MAX(MXMEMT,MXMEML)
            MXMEML = 0

C           Go to next a batch.
C           (Which will exit immediately if no more a's).
C           ---------------------------------------------

            GO TO 110

C        Exit point for a batch loop.
C        ----------------------------

  190    CONTINUE

C        Print info for this b-batch and got to next.
C        (Which will exit immediately if no more b's.)
C        ---------------------------------------------

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(5X,A,A,/)')
     &      '------------------------------------------------------',
     &      '---------------'
            TBBAT = SECOND() - TBBAT
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Total time for this b-batch: ',TBBAT,' seconds'
            WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &      'Accumulated amplitude 2-norm : ',DSQRT(X2NRM),
     &      'Accumulated 2CME am.  2-norm : ',DSQRT(X2CNM)
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum number of Cholesky batches: ',MNCHBT,
     &      'Maximum number of Cholesky batches: ',MXCHBT
            CALL FLSHFO(LUPRI)
         ENDIF

         GO TO 100

C     Exit point for b batch loop.
C     ----------------------------

  200 CONTINUE

C     If requested, delete Cholesky files.
C     ------------------------------------

      DO IP12 = 1,NUMP12
         IF (DELP1(IP12)) THEN
            DO ISYCHO = 1,NSYM
               CALL CHO_MOP(IOPEN,JTYP1(IP12),ISYCHO,LUCHO1,1,
     &                      ISYMP1(IP12))
               CALL CHO_MOP(IDEL,JTYP1(IP12),ISYCHO,LUCHO1,1,
     &                      ISYMP1(IP12))
            ENDDO
         ENDIF
         IF (DELP2(IP12)) THEN
            DO ISYCHO = 1,NSYM
               CALL CHO_MOP(IOPEN,JTYP2(IP12),ISYCHO,LUCHO2,1,
     &                      ISYMP2(IP12))
               CALL CHO_MOP(IDEL,JTYP2(IP12),ISYCHO,LUCHO2,1,
     &                      ISYMP2(IP12))
            ENDDO
         ENDIF
      ENDDO

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN

         TIMCYI(1) = SECOND() - TIMCYI(1)

         CALL HEADER(SECNAM//' Timing by Routine',-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CIA    : ',TIMCYI(2),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYILHX : ',TIMCYI(3),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_DNOM   : ',TIMCYI(4),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYINRM : ',TIMCYI(5),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYITCME: ',TIMCYI(6),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYICNI : ',TIMCYI(7),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYIM   : ',TIMCYI(8),' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time             : ',TIMCYI(1),' seconds'

         XTTOT = TMCYIO(1) + TMCYIO(2) + TMCALC(1) + TMCALC(2)
         XI1 = D100*TMCYIO(1)/TIMCYI(1)
         XI2 = D100*TMCYIO(2)/TIMCYI(1)
         XC1 = D100*TMCALC(1)/TIMCYI(1)
         XC2 = D100*TMCALC(2)/TIMCYI(1)
         XT  = D100*XTTOT/TIMCYI(1)
         CALL HEADER(SECNAM//' Timing by Task',-1)
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used for MO Cholesky I/O         : ',TMCYIO(1),
     &   ' seconds (',XI1,' %)'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used for Y intermediate I/O      : ',TMCYIO(2),
     &   ' seconds (',XI2,' %)'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used in DGEMM for integrals      : ',TMCALC(1),
     &   ' seconds (',XC1,' %)'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used in DGEMM for Y intermediates: ',TMCALC(2),
     &   ' seconds (',XC2,' %)'
         WRITE(LUPRI,'(5X,A,A)')
     &   '---------------------------------------------------',
     &   '-------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Total                                 : ',XTTOT,
     &   ' seconds (',XT,' %)'

         CALL FLSHFO(LUPRI)

      ENDIF

      RETURN
      END
C  /* Deck cc_cyinp1 */
      SUBROUTINE CC_CYINP1(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &                     UVEC,ISYMU,VVEC,ISYMV,FAC2,NUMUV,
     &                     SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &                     JTYPZ,ISYMZ,NCHOZ,NUMZ,JTYPY,
     &                     WVEC,ISYMW,NUMW,SVEC,
     &                     WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2,IERR,PRNT,
     &                     ISYINT)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Check input to CC_CYI1.
C              For practical purpose (allocation),
C              NUMP12 <= 0 is considered an error.
C
#include <implicit.h>
      DIMENSION UVEC(*), VVEC(*), WVEC(*), SVEC(*), FOCKD(*)
      DIMENSION SCD(NUMP12), WORK(LWORK)
      INTEGER JTYP1(NUMP12), ISYMP1(NUMP12)
      INTEGER JTYP2(NUMP12), ISYMP2(NUMP12)
      INTEGER NCHP12(8,NUMP12)
      INTEGER ISYMU(NUMUV), ISYMV(NUMUV)
      INTEGER JTYPZ(NUMZ), ISYMZ(NUMZ), NCHOZ(8,NUMZ), JTYPY(NUMZ)
      INTEGER ISYMW(NUMW)
      LOGICAL DELP1(NUMP12), DELP2(NUMP12)
      LOGICAL PRNT
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYINP1')

      PARAMETER (SMALL = 1.0D-4)

C     Set IERR.
C     ---------

      IERR = 0

C     Check symmetries.
C     -----------------

      IF (NUMP12 .GT. 0) THEN

         ISYINT = MULD2H(ISYMP1(1),ISYMP2(1))
         IF ((ISYINT.LT.1) .OR. (ISYINT.GT.NSYM)) THEN
            IERR = 1
            GO TO 100
         ENDIF

         DO I = 2,NUMP12
            ISYTST = MULD2H(ISYMP1(I),ISYMP2(I))
            IF (ISYTST .NE. ISYINT) THEN
               IERR = 2
               GO TO 100
            ENDIF
         ENDDO

         DO I = 1,NUMUV
            ISYTST = MULD2H(ISYMU(I),ISYMV(I))
            IF (ISYTST .NE. ISYINT) THEN
               IERR = 2
               GO TO 100
            ENDIF
         ENDDO

      ELSE

         IF (NUMUV .GT. 0) THEN

            ISYINT = MULD2H(ISYMU(1),ISYMV(1))
            IF ((ISYINT.LT.1) .OR. (ISYINT.GT.NSYM)) THEN
               IERR = 1
               GO TO 100
            ENDIF

            DO I = 2,NUMUV
               ISYTST = MULD2H(ISYMU(I),ISYMV(I))
               IF (ISYTST .NE. ISYINT) THEN
                  IERR = 2
                  GO TO 100
               ENDIF
            ENDDO

         ELSE

            IERR = -1
            GO TO 100

         ENDIF

C        For reasons of allocation, NUMP12 <= 0 must be
C        considered an error.
C        ----------------------------------------------

         IF (IERR .EQ. 0) THEN
            IERR = -3
            GO TO 100
         ENDIF

      ENDIF

C     Check that Y intermediates or U vectors are to be calculated.
C     -------------------------------------------------------------

      IF ((NUMZ.LE.0) .AND. (NUMW.LE.0)) THEN
         IERR = -2
         GO TO 100
      ENDIF

C     Print section.
C     --------------

  100 IF (PRNT .OR. (IERR.NE.0)) THEN

         WRITE(LUPRI,'(//,1X,A,A,A)')
     &   'Input to CC_CYI1 checked by ',SECNAM,':'
         WRITE(LUPRI,*) 'Error parameter IERR = ',IERR
         WRITE(LUPRI,*) 'Integral symmetry: ',ISYINT
         WRITE(LUPRI,*) 'NUMP12: ',NUMP12
         WRITE(LUPRI,*) 'JTYP1 : ',(JTYP1(I),I=1,NUMP12)
         WRITE(LUPRI,*) 'ISYMP1: ',(ISYMP1(I),I=1,NUMP12)
         WRITE(LUPRI,*) 'JTYP2 : ',(JTYP2(I),I=1,NUMP12)
         WRITE(LUPRI,*) 'ISYMP2: ',(ISYMP2(I),I=1,NUMP12)
         WRITE(LUPRI,*) 'SCD   : ',(SCD(I),I=1,NUMP12)
         WRITE(LUPRI,*) 'NCHP12:'
         DO I = 1,NUMP12
            WRITE(LUPRI,*) (NCHP12(J,I), J=1,NSYM)
         ENDDO
         WRITE(LUPRI,*) 'ISYMU : ',(ISYMU(I),I=1,NUMUV)
         WRITE(LUPRI,*) 'ISYMV : ',(ISYMV(I),I=1,NUMUV)
         WRITE(LUPRI,*) 'FAC1  : ',FAC1
         WRITE(LUPRI,*) 'FAC2  : ',FAC2
         WRITE(LUPRI,*) 'SCDG  : ',SCDG
         WRITE(LUPRI,*) 'FREQ  : ',FREQ
         WRITE(LUPRI,*) 'IOPTDN: ',IOPTDN
         WRITE(LUPRI,*) 'IOPTCE: ',IOPTCE
         WRITE(LUPRI,*) 'JTYPZ : ',(JTYPZ(I),I=1,NUMZ)
         WRITE(LUPRI,*) 'ISYMZ : ',(ISYMZ(I),I=1,NUMZ)
         WRITE(LUPRI,*) 'JTYPY : ',(JTYPY(I),I=1,NUMZ)
         WRITE(LUPRI,*) 'ISYMW : ',(ISYMW(I),I=1,NUMW)
         WRITE(LUPRI,*) 'LWORK : ',LWORK
         WRITE(LUPRI,*) 'DELP1 : ',(DELP1(I),I=1,NUMP12)
         WRITE(LUPRI,*) 'DELP2 : ',(DELP2(I),I=1,NUMP12)

         CALL FLSHFO(LUPRI)

      ENDIF

      RETURN
      END
C  /* Deck cc_dnom */
      SUBROUTINE CC_DNOM(FOCKD,XAM,FREQ,ISYMX)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate general doubles amplitudes as
C
C        XAM(#ai,#bj)
C         = - XAM(#ai,#bj)/[FOCKD(#a)-FOCKD(i)+FOCKD(#b)-FOCKD(j)-FREQ]
C
C              for batches of virtual orbitals.
C
C     Input:
C
C        XAM  : Integrals to be divided by orbital
C               energy denominators to obtain amplitude batch.
C        FOCKD: Canonical orbital energies.
C        FREQ : Perturbation frequency.
C               NOTE the sign in the formula above.
C        ISYMX: Symmetry of the amplitudes/integrals.
C
C     Information about the batches of virtual orbitals is assumed available
C     in the configuration file ciarc.h (see CC_CIA for details).
C
#include <implicit.h>
      DIMENSION FOCKD(*), XAM(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ciarc.h>
#include <priunit.h>

c      write(lupri,*) 'CC_DNOM: freq,isymx: ',freq,isymx
c      write(lupri,*) 'CC_DNOM: nx2sq     : ',nx2sq
c      write(lupri,*) 'CC_DNOM: koff1     : ',koff1
c      call flshfo(lupri)

      DO ISYMBJ = 1,NSYM

         ISYMAI = MULD2H(ISYMBJ,ISYMX)

c      write(lupri,*) 'CC_DNOM: isymai,isymbj: ',isymai,isymbj
c      write(lupri,*) 'CC_DNOM: nx1ama,nx1amb: ',nx1ama(isymai),
c     &                                          nx1amb(isymbj)
c      call flshfo(lupri)

         IF ((NX1AMB(ISYMBJ).GT.0).AND.(NX1AMA(ISYMAI).GT.0)) THEN

            DO ISYMJ = 1,NSYM

               ISYMB = MULD2H(ISYMJ,ISYMBJ)

               IF (LVIRB(ISYMB) .GT. 0) THEN

                  DO J = 1,NRHF(ISYMJ)

                     KOFFJ = IRHF(ISYMJ) + J

                     DO LB = 1,LVIRB(ISYMB)

                        B     = IOFB1(ISYMB) + LB - 1
                        KOFFB = IVIR(ISYMB)  + B

                        DBJ = FOCKD(KOFFB) - FOCKD(KOFFJ)

                        LBJ = IX1AMB(ISYMB,ISYMJ)
     &                      + LVIRB(ISYMB)*(J - 1) + LB

                        DO ISYMI = 1,NSYM

                           ISYMA = MULD2H(ISYMI,ISYMAI)

                           DO I = 1,NRHF(ISYMI)

                              KOFFI = IRHF(ISYMI) + I

                              DIBJ = DBJ - FOCKD(KOFFI)

                              DO LA = 1,LVIRA(ISYMA)

                                 A     = IOFA1(ISYMA) + LA - 1
                                 KOFFA = IVIR(ISYMA)  + A

                                 DAIBJ = FOCKD(KOFFA) + DIBJ
                                 DENOM = FREQ - DAIBJ

                                 LAI = IX1AMA(ISYMA,ISYMI)
     &                               + LVIRA(ISYMA)*(I - 1) + LA

                                 LAIBJ = IX2SQ(ISYMAI,ISYMBJ)
     &                                 + NX1AMA(ISYMAI)*(LBJ - 1)
     &                                 + LAI

c      write(lupri,*) 'aibj : ',laibj
c      write(lupri,*) 'denom: ',denom
c      call flshfo(lupri)

                                 XAM(LAIBJ) = XAM(LAIBJ)/DENOM

                              ENDDO

                           ENDDO

                        ENDDO

                     ENDDO

                  ENDDO

               ENDIF

            ENDDO

         ENDIF

      ENDDO

      RETURN
      END
C  /* Deck cc_cyimsp */
      SUBROUTINE CC_CYIMSP(LWORK,LWRK)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Split memory for triple batching in
C              Cholesky Y intermediate calculation based
C              on integral calculation in CC_CIA.
C
C     Input:
C
C        LWORK  : Work space available.
C
C        Logical control variables and # Cholesky vectors
C        in ciarc.h
C
C     Output:
C
C        LWRK   : Dimension of work space reserved for
C                 X(#ai,#bj) integral/amplitude batch.
C
C     Note that this routine makes assumptions about the memory
C     allocations in CC_CIA.
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccorb.h>
#include <dccsdsym.h>
#include <priunit.h>
#include <ciarc.h>
#include <chocc2.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIMSP')

      PARAMETER (MINDEF = 5, LFRAC = 5)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THREE = 3.0D0, FOUR = 4.0D0)

C     Test LWORK.
C     ===========

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &   'Work space non-positive in ',SECNAM,': LWORK = ',LWORK
         CALL QUIT('Error in '//SECNAM)
      ENDIF
      XWORK  = ONE*LWORK
      ISYINT = MULD2H(ISYCH1,ISYCH2)

C     Set Cholesky weight.
C     ====================

      WCHOL = SPLITC

      IF (WCHOL .LE. ZERO) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,1P,D22.15,/)')
     &   'Error in ',SECNAM,':',
     &   'Cholesky weight in memory split, WCHOL = ',WCHOL
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     See if we can hold the entire integral/amplitude matrix
C     in core, while at the same time storing at least 1/LFRAC
C     of the Cholesky vectors in each symmetry.
C     ========================================================

      XLFT = XWORK - XT2SQ(ISYINT)

      IF (XLFT .GT. ZERO) THEN

         LWRK = NT2SQ(ISYINT)
         LEFT = LWORK - LWRK

         IF (LIDEN) THEN

C           L1 = L2 (CC_CIA1).
C           ------------------

            DO ISYCHO = 1,NSYM
               IF (NTOVEC(ISYCHO) .GT. LFRAC) THEN
                  MINVEC = NTOVEC(ISYCHO)/LFRAC
               ELSE
                  MINVEC = MIN(NTOVEC(ISYCHO),MINDEF)
               ENDIF
               ISYMAI = MULD2H(ISYCHO,ISYCH1)
               NEED   = NT1AM(ISYMAI)*MINVEC
               IF (NEED .GT. LEFT) GO TO 100
            ENDDO

         ELSE

C           L1 != L2 (CC_CIA3).
C           -------------------

            DO ISYCHO = 1,NSYM
               IF (NTOVEC(ISYCHO) .GT. LFRAC) THEN
                  MINVEC = NTOVEC(ISYCHO)/LFRAC
               ELSE
                  MINVEC = MIN(NTOVEC(ISYCHO),MINDEF)
               ENDIF
               ISYMAI = MULD2H(ISYCHO,ISYCH1)
               ISYMBJ = MULD2H(ISYCHO,ISYCH2)
               NEED   = (NT1AM(ISYMAI) + NT1AM(ISYMBJ))*MINVEC
               IF (NEED .GT. LEFT) GO TO 100
            ENDDO

         ENDIF

C        If we reach this point, full calculation is possible.
C        -----------------------------------------------------

         RETURN

      ENDIF

C     Ressume processing here if batching is needed.
C     ----------------------------------------------

  100 CONTINUE

C     Calculate memory split IDIV according to relative
C     sizes of integrals and Cholesky vectors.
C     =================================================

      MAXIJ = 0
      DO ISYM = 1,NSYM
         MAXIJ = MAX(MAXIJ,NMATIJ(ISYM))
      ENDDO

      IF (LIDEN) THEN

C        L1 = L2 (CC_CIA2).
C        Integrals/amplitudes: XT2SQ(1)
C        MO Cholesky vectors : 2*L1(ai,J) + one additional vector.
C        ---------------------------------------------------------

         XCHO = ZERO
         DO ISYCHO = 1,NSYM
            ISYMAI = MULD2H(ISYCHO,ISYCH1)
            XTST   = XT1AM(ISYMAI)*(2*NTOVEC(ISYCHO) + 1)
            IF (XCHO .LT. XTST) THEN
               XCHO = XTST
            ENDIF
         ENDDO

         IF ((XCHO.LE.ZERO) .OR. (XT2SQ(1).LE.ZERO)) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Error in ',SECNAM,':'
            IF (XCHO .LE. ZERO) THEN
               WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &         'Calculated max. dimension of Cholesky part: ',XCHO
               WRITE(LUPRI,'(5X,A)')
     &         'Number of Cholesky vectors supplied:'
               WRITE(LUPRI,'(8I10)')
     &         (NTOVEC(ISYCHO),ISYCHO=1,NSYM)
            ENDIF
            IF (XT2SQ(1) .LE. ZERO) THEN
               WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &         'Integral dimension from dccsdsym.h: ',XT2SQ(1)
               WRITE(LUPRI,'(5X,A,I2)')
     &         'Integral symmetry:',ISYINT
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

C        Split memory.
C        -------------

         LWRK = LWORK/IDIV
         IF (ISPLIT .EQ. 1) THEN
            LWRK = LWORK - LWRK
         ENDIF

C        Make sure that we do not reserve too much space
C        for Cholesky vectors.
C        -----------------------------------------------

         LEFT = LWORK - LWRK
         XLFT = ONE*LEFT
         IF (XLFT .GT. XCHO) THEN
            LEFT = INT(XCHO)
            LWRK = LWORK - LEFT
         ENDIF

C        Make sure that integral/amplitude batch can be held in
C        core for at least 1 virtual. This section should be
C        relevant only for work spaces in the neighborhood of
C        the absolute minimal required.
C        ------------------------------------------------------

         IF (LWRK .LT. MAXIJ) THEN
            LWRK = MIN(LWORK,MAXIJ)
         ENDIF

         LEFT = LWORK - LWRK
         DO ISYMB = 1,NSYM
            DO ISYMA = 1,ISYMB
               DO ISYCHO = 1,NSYM
                  ISYMAI = MULD2H(ISYCHO,ISYCH1)
                  ISYMBJ = ISYMAI
                  ISYMI  = MULD2H(ISYMAI,ISYMA)
                  ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                  NTST   = NT1AM(ISYMAI) + NRHF(ISYMI) + NRHF(ISYMJ)
                  IF (LEFT .LT. NTST) THEN
                     LEFT = MIN(LWORK,NTST)
                     LWRK = LWORK - LEFT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

      ELSE

         IF (SYMTRZ .AND. CIAMIO) THEN

C           L1 != L2, symmetrization with minmum I/O (CC_CIA4_1).
C           Integrals/amplitudes: XT2SQ(ISYINT)
C           MO Cholesky vectors : 2*(L1(ai,J) + L2(bj,J)) + 2 additional.
C           -------------------------------------------------------------

            XCHO = ZERO
            DO ISYCHO = 1,NSYM
               ISYMAI = MULD2H(ISYCHO,ISYCH1)
               ISYMBJ = MULD2H(ISYCHO,ISYCH2)
               XT1TOT = XT1AM(ISYMAI) + XT1AM(ISYMBJ)
               XTST   = XT1TOT*(2*NTOVEC(ISYCHO) + 1)
               IF (XCHO .LT. XTST) XCHO = XTST
            ENDDO

            IF ((XCHO.LE.ZERO) .OR. (XT2SQ(ISYINT).LE.ZERO)) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Error in ',SECNAM,':'
               IF (XCHO .LE. ZERO) THEN
                  WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &            'Calculated max. dimension of Cholesky part: ',XCHO
                  WRITE(LUPRI,'(5X,A)')
     &            'Number of Cholesky vectors supplied:'
                  WRITE(LUPRI,'(8I10)')
     &            (NTOVEC(ISYCHO),ISYCHO=1,NSYM)
               ENDIF
               IF (XT2SQ(ISYINT) .LE. ZERO) THEN
                  WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &            'Integral dimension from dccsdsym.h: ',XT2SQ(ISYINT)
                  WRITE(LUPRI,'(5X,A,I2)')
     &            'Integral symmetry:',ISYINT
               ENDIF
               WRITE(LUPRI,*)
               CALL QUIT('Error in '//SECNAM)
            ENDIF

            XCHOF = WCHOL*XCHO
            IF (XT2SQ(ISYINT) .GT. XCHOF) THEN
               ISPLIT = 1
               XFRAC  = XT2SQ(ISYINT)/XCHOF
            ELSE
               ISPLIT = 2
               XFRAC  = XCHOF/XT2SQ(ISYINT)
            ENDIF

            IFRAC = INT(XFRAC)
            IDIV  = MAX(IFRAC+1,2)

C           Split memory.
C           -------------

            LWRK = LWORK/IDIV
            IF (ISPLIT .EQ. 1) THEN
               LWRK = LWORK - LWRK
            ENDIF

C           Make sure that we do not reserve too much space
C           for Cholesky vectors.
C           -----------------------------------------------

            LEFT = LWORK - LWRK
            XLFT = ONE*LEFT
            IF (XLFT .GT. XCHO) THEN
               LEFT = INT(XCHO)
               LWRK = LWORK - LEFT
            ENDIF

C           Make sure that integral/amplitude batch can be held in
C           core for at least 1 virtual. This section should be
C           relevant only for work spaces in the neighborhood of
C           the absolute minimal required.
C           ------------------------------------------------------

            IF (LWRK .LT. MAXIJ) THEN
               LWRK = MIN(LWORK,MAXIJ)
            ENDIF

            LEFT = LWORK - LWRK
            DO ISYMB = 1,NSYM
               DO ISYMA = 1,ISYMB
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  DO ISYCHO = 1,NSYM
                     ISYMAI = MULD2H(ISYCHO,ISYCH1)
                     ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                     ISYMI1 = MULD2H(ISYMAI,ISYMA)
                     ISYMJ1 = MULD2H(ISYMBJ,ISYMB)
                     ISYMI2 = MULD2H(ISYMBJ,ISYMA)
                     ISYMJ2 = MULD2H(ISYMAI,ISYMB)
                     NTST   = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
     &                      + NRHF(ISYMI1)  + NRHF(ISYMJ1)
     &                      + NRHF(ISYMI2)  + NRHF(ISYMJ2)
                     IF (LEFT .LT. NTST) THEN
                        LEFT = MIN(LWORK,NTST)
                        LWRK = LWORK - LEFT
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

         ELSE

C           L1 != L2, no symmetrization, or symmetrization by 2 passes
C                     through Cholesky files (CC_CIA4_3).
C           Integrals/amplitudes: XT2SQ(ISYINT)
C           MO Cholesky vectors : L1(ai,J) + L2(bj,J) + 2 additional.
C           ----------------------------------------------------------

            XCHO = ZERO
            DO ISYCHO = 1,NSYM
               ISYMAI = MULD2H(ISYCHO,ISYCH1)
               ISYMBJ = MULD2H(ISYCHO,ISYCH2)
               XT1TOT = XT1AM(ISYMAI) + XT1AM(ISYMBJ)
               XTST   = XT1TOT*(NTOVEC(ISYCHO) + 1)
               IF (XCHO .LT. XTST) XCHO = XTST
            ENDDO

            IF ((XCHO.LE.ZERO) .OR. (XT2SQ(ISYINT).LE.ZERO)) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Error in ',SECNAM,':'
               IF (XCHO .LE. ZERO) THEN
                  WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &            'Calculated max. dimension of Cholesky part: ',XCHO
                  WRITE(LUPRI,'(5X,A)')
     &            'Number of Cholesky vectors supplied:'
                  WRITE(LUPRI,'(8I10)')
     &            (NTOVEC(ISYCHO),ISYCHO=1,NSYM)
               ENDIF
               IF (XT2SQ(ISYINT) .LE. ZERO) THEN
                  WRITE(LUPRI,'(5X,A,1P,D30.15)')
     &            'Integral dimension from dccsdsym.h: ',XT2SQ(ISYINT)
                  WRITE(LUPRI,'(5X,A,I2)')
     &            'Integral symmetry:',ISYINT
               ENDIF
               WRITE(LUPRI,*)
               CALL QUIT('Error in '//SECNAM)
            ENDIF

            XCHOF = WCHOL*XCHO
            IF (XT2SQ(ISYINT) .GT. XCHOF) THEN
               ISPLIT = 1
               XFRAC  = XT2SQ(ISYINT)/XCHOF
            ELSE
               ISPLIT = 2
               XFRAC  = XCHOF/XT2SQ(ISYINT)
            ENDIF

            IFRAC = INT(XFRAC)
            IDIV  = MAX(IFRAC+1,2)

C           Split memory.
C           -------------

            LWRK = LWORK/IDIV
            IF (ISPLIT .EQ. 1) THEN
               LWRK = LWORK - LWRK
            ENDIF

C           Make sure that we do not reserve too much space
C           for Cholesky vectors.
C           -----------------------------------------------

            LEFT = LWORK - LWRK
            XLFT = ONE*LEFT
            IF (XLFT .GT. XCHO) THEN
               LEFT = INT(XCHO)
               LWRK = LWORK - LEFT
            ENDIF

C           Make sure that integral/amplitude batch can be held in
C           core for at least 1 virtual. This section should be
C           relevant only for work spaces in the neighborhood of
C           the absolute minimal required.
C           ------------------------------------------------------

            IF (LWRK .LT. MAXIJ) THEN
               LWRK = MIN(LWORK,MAXIJ)
            ENDIF

            LEFT = LWORK - LWRK
            DO ISYMB = 1,NSYM
               DO ISYMA = 1,ISYMB
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  DO ISYCHO = 1,NSYM
                     ISYMAI = MULD2H(ISYCHO,ISYCH1)
                     ISYMBJ = MULD2H(ISYCHO,ISYCH2)
                     ISYMI  = MULD2H(ISYMAI,ISYMA)
                     ISYMJ  = MULD2H(ISYMBJ,ISYMB)
                     NTST   = NT1AM(ISYMAI) + NT1AM(ISYMBJ)
     &                      + NRHF(ISYMI)   + NRHF(ISYMJ)
                     IF (LEFT .LT. NTST) THEN
                        LEFT = MIN(LWORK,NTST)
                        LWRK = LWORK - LEFT
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

         ENDIF

      ENDIF

      RETURN
      END
C  /* Deck cc_cyibtb */
      SUBROUTINE CC_CYIBTB(LVIRB,BBEG,NUMB,LWRK,ISYINT)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Set up b-batch for CC_CYI, estimating the number
C              of a's along the way.
C
C     Generalized version of GET_BTCHB in MP2 module.
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

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIBTB')

C     Check that first b is positive.
C     -------------------------------

      IF (BBEG .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,I10,/)')
     &   'Error in ',SECNAM,':',
     &   'BBEG = ',BBEG
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

      A = 0
      B = BBEG - 1

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

         MEM = MEMABG(LVIRA,LVIRB,ISYINT)
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
            NUMB = NUMB - 1
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

            MEM = MEMABG(LVIRA,LVIRB,ISYINT)
            IF (MEM .GT. LWRK) THEN

  300          CONTINUE

                  IF (NUMB .LE. 0) RETURN

                  JSYMB = ISVI(B)
                  LVIRB(JSYMB) = LVIRB(JSYMB) - 1
                  B = B - 1
                  NUMB = NUMB - 1

                  MEM = MEMABG(LVIRA,LVIRB,ISYINT)
                  IF (MEM .GT. LWRK) GO TO 300

            ENDIF

         ENDDO
      ENDIF

      RETURN
      END
C  /* Deck memabg */
      INTEGER FUNCTION MEMABG(LVIRA,LVIRB,ISYINT)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Evaluate memory required for (#ai|#bj)-type integrals
C              (of symmetry ISYINT) with virtual orbitals specified
C              LVIRA and LVIRB.
C
C     Generalized version of MEMAB in the MP2 module.
C
#include <implicit.h>
      INTEGER LVIRA(8), LVIRB(8)
#include <ccorb.h>

      MEMABG = 0
      DO ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYINT)
         ICOUBJ = 0
         ICOUAI = 0
         DO ISYMO = 1,NSYM
            ISYMB  = MULD2H(ISYMO,ISYMBJ)
            ISYMA  = MULD2H(ISYMO,ISYMAI)
            ICOUBJ = ICOUBJ + LVIRB(ISYMB)*NRHF(ISYMO)
            ICOUAI = ICOUAI + LVIRA(ISYMA)*NRHF(ISYMO)
         ENDDO
         MEMABG = MEMABG + ICOUAI*ICOUBJ
      ENDDO

      RETURN
      END
C  /* Deck cc_cyibta */
      SUBROUTINE CC_CYIBTA(LIBJ,NX1AMB,LVIRA,ABEG,NUMA,MEM,LWRK,ISYINT)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Set up a-batch for CC_CYI.
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
      PARAMETER (SECNAM = 'CC_CYIBTA')

C     Check that first a is positive.
C     -------------------------------

      IF (ABEG .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,I10,/)')
     &   'Error in ',SECNAM,':',
     &   'ABEG = ',ABEG
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

         ISYMA  = ISVI(A)
         ISYIBJ = MULD2H(ISYMA,ISYINT)
         MEM    = MEM  + LIBJ(ISYIBJ)
         NUMA   = NUMA + 1
         LVIRA(ISYMA) = LVIRA(ISYMA) + 1

         IF (MEM .LE. LWRK) THEN
            GO TO 100
         ELSE
            MEM  = MEM  - LIBJ(ISYIBJ)
            NUMA = NUMA - 1
            LVIRA(ISYMA) = LVIRA(ISYMA) - 1
            A = A - 1
            GO TO 200
         ENDIF

C     TODO: at this point we might want to check that
C           at least one Cholesky vector can be stored
C           in each symmetry alongside the integral batch.
C           For now, we simply trust the intial memory
C           split to have done this properly.
C     ----------------------------------------------------

  200 RETURN
      END
C  /* Deck cc_cyilhx */
      SUBROUTINE CC_CYILHX(X1AM,XINT,FIA,WORK,LWORK,MEM,ISYMX1,ISYMF,
     &                     FAC)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Add extra terms needed in CC_CYI for left-hand trf.
C
C              X(#ai,#bj) = X(#ai,#bj) + P(ai,bj) FAC * X1AM(ai) * FIA(jb)
C
#include <implicit.h>
      DIMENSION X1AM(*), XINT(*), FIA(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ciarc.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYILHX')

      LOGICAL FULLIN

      PARAMETER (TINY = 1.0D-14, ONE = 1.0D0)

C     Is this the full matrix?
C     ------------------------

      ISYINT = MULD2H(ISYMX1,ISYMF)
      XX2SQ  = ONE*NX2SQ
      DIFF   = XT2SQ(ISYINT) - XX2SQ
      FULLIN = DABS(DIFF) .LE. TINY

C     Calculate extra terms.
C     ----------------------

      IF (FULLIN) THEN

         MEM = 0

         KOFFX = IT2SQ(ISYMF,ISYMX1) + 1
         NTOAI = MAX(NT1AM(ISYMF),1)
         NTOBJ = MAX(NT1AM(ISYMX1),1)
         CALL DGEMM('N','T',NT1AM(ISYMF),NT1AM(ISYMX1),1,
     &              FAC,FIA,NTOAI,X1AM,NTOBJ,ONE,XINT(KOFFX),NTOAI)

         KOFFX = IT2SQ(ISYMX1,ISYMF) + 1
         NTMP  = NTOAI
         NTOAI = NTOBJ
         NTOBJ = NTMP
         CALL DGEMM('N','T',NT1AM(ISYMX1),NT1AM(ISYMF),1,
     &              FAC,X1AM,NTOAI,FIA,NTOBJ,ONE,XINT(KOFFX),NTOAI)

      ELSE

         KF1   = 1
         KX11  = KF1  + NX1AMA(ISYMF)
         KEND1 = KX11 + NX1AMB(ISYMX1)

         KX12  = 1
         KF2   = KX12 + NX1AMA(ISYMX1)
         KEND2 = KF2  + NX1AMB(ISYMF)

         MEM   = MAX(KEND1,KEND2) - 1
         LWRK  = LWORK - MEM

         IF (LWRK .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',MEM,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL CC_CIACPS(FIA,WORK(KF1),1,IOFA1,LVIRA,NX1AMA,IX1AMA,ISYMF)
         CALL CC_CIACPS(X1AM,WORK(KX11),1,IOFB1,LVIRB,NX1AMB,IX1AMB,
     &                  ISYMX1)

         KOFFX = IX2SQ(ISYMF,ISYMX1) + 1
         NTOAI = MAX(NX1AMA(ISYMF),1)
         NTOBJ = MAX(NX1AMB(ISYMX1),1)
         CALL DGEMM('N','T',NX1AMA(ISYMF),NX1AMB(ISYMX1),1,
     &              FAC,WORK(KF1),NTOAI,WORK(KX11),NTOBJ,
     &              ONE,XINT(KOFFX),NTOAI)

         CALL CC_CIACPS(X1AM,WORK(KX12),1,IOFA1,LVIRA,NX1AMA,IX1AMA,
     &                  ISYMX1)
         CALL CC_CIACPS(FIA,WORK(KF2),1,IOFB1,LVIRB,NX1AMB,IX1AMB,
     &                  ISYMF)

         KOFFX = IX2SQ(ISYMX1,ISYMF) + 1
         NTOAI = MAX(NX1AMA(ISYMX1),1)
         NTOBJ = MAX(NX1AMB(ISYMF),1)
         CALL DGEMM('N','T',NX1AMA(ISYMX1),NX1AMB(ISYMF),1,
     &              FAC,WORK(KX12),NTOAI,WORK(KF2),NTOBJ,
     &              ONE,XINT(KOFFX),NTOAI)

      ENDIF

      RETURN
      END
C  /* Deck cc_cyitcme */
      SUBROUTINE CC_CYITCME(X2AM,WORK,LWORK,ISYMX,MEM)
C
C     Thomas Bondo Pedersen, January 2003
C
C     Purpose: Set up Two Coulomb Minus Exchange of a batch of doubles
C              amplitudes X2AM(#ai,#bj) of symmetry ISYMX.
C              Work space needed is appr. 2*(#a); actual space used
C              is returned in MEM.
C
C     Information about the a- and b-batches is taken from ciarc.h
C
C     The overall structure is taken from the standard "exchange"
C     CCSD_T2TP routine by A. Sanchez and H. Koch.
C
#include <implicit.h>
      DIMENSION X2AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <ciarc.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CYITCME')

      PARAMETER (XMONE = -1.0D0, TWO = 2.0D0)

      MEM = 0

      DO ISYMJ = 1,NSYM
         DO J = 1,NRHF(ISYMJ)
            DO ISYMB = 1,NSYM

               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYMX)

               DO LB = 1,LVIRB(ISYMB)

                  LBJ = IX1AMB(ISYMB,ISYMJ) + LVIRB(ISYMB)*(J - 1) + LB

                  DO ISYMI = 1,ISYMJ

                     ISYMA  = MULD2H(ISYMI,ISYMAI)

                     IF (LVIRA(ISYMA) .GT. 0) THEN

                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)

                        KAIBJ = 1
                        KAJBI = KAIBJ + LVIRA(ISYMA)
                        KEND1 = KAJBI + LVIRA(ISYMA)
                        LWRK1 = LWORK - KEND1 + 1

                        KLAST = KEND1 - 1
                        MEM   = MAX(MEM,KLAST)

                        IF (LWRK1 .LT. 0) THEN
                      WRITE(LUPRI,'(//,5X,A,A,/,5X,A,I10,/,5X,A,I10,/)')
     &                     'Insufficient memory in ',SECNAM,
     &                     'Need     : ',KLAST,
     &                     'Available: ',LWORK
                           CALL QUIT('Insufficient memory in '//SECNAM)
                        ENDIF

                        IF (ISYMI .EQ. ISYMJ) THEN
                           NRHFI = J - 1
                        ELSE
                           NRHFI = NRHF(ISYMI)
                        ENDIF

                        DO I = 1,NRHFI

                           LBI = IX1AMB(ISYMB,ISYMI)
     &                         + LVIRB(ISYMB)*(I - 1) + LB

                           NAIBJ = IX2SQ(ISYMAI,ISYMBJ)
     &                           + NX1AMA(ISYMAI)*(LBJ - 1)
     &                           + IX1AMA(ISYMA,ISYMI)
     &                           + LVIRA(ISYMA)*(I - 1) + 1
                           NAJBI = IX2SQ(ISYMAJ,ISYMBI)
     &                           + NX1AMA(ISYMAJ)*(LBI - 1)
     &                           + IX1AMA(ISYMA,ISYMJ)
     &                           + LVIRA(ISYMA)*(J - 1) + 1

                           CALL DCOPY(LVIRA(ISYMA),X2AM(NAIBJ),1,
     &                                             WORK(KAIBJ),1)
                           CALL DCOPY(LVIRA(ISYMA),X2AM(NAJBI),1,
     &                                             WORK(KAJBI),1)
                           CALL DSCAL(LVIRA(ISYMA),TWO,X2AM(NAIBJ),1)
                           CALL DSCAL(LVIRA(ISYMA),TWO,X2AM(NAJBI),1)
                           CALL DAXPY(LVIRA(ISYMA),XMONE,WORK(KAJBI),1,
     &                                                   X2AM(NAIBJ),1)
                           CALL DAXPY(LVIRA(ISYMA),XMONE,WORK(KAIBJ),1,
     &                                                   X2AM(NAJBI),1)

                        ENDDO

                     ENDIF

                  ENDDO

               ENDDO

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_cyinrm */
      SUBROUTINE CC_CYINRM(X2AM,ISYMX,X2NRM)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate contribution to packed norm of
C              the amplitude batch X2AM(#ai,#bj).
C
#include <implicit.h>
      DIMENSION X2AM(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ciarc.h>
#include <priunit.h>

      INTEGER AI, BJ, AI1, AI2, BJ1, BJ2

      IF (ISYMX .EQ. 1) THEN

         DO ISYMBJ = 1,NSYM

            ISYMAI = ISYMBJ

            IF ((NX1AMA(ISYMAI).GT.0) .AND. (NX1AMB(ISYMBJ).GT.0)) THEN

               CALL CC_CYIFST(IOFA1,LVIRA,ISYMAI,AI1)
               CALL CC_CYILST(IOFB1,LVIRB,ISYMBJ,BJ2)

               IF (AI1 .GT. BJ2) GO TO 1000

               DO ISYMJ = 1,NSYM

                  ISYMB = MULD2H(ISYMJ,ISYMBJ)

                  IF ((LVIRB(ISYMB).GT.0) .AND. (NRHF(ISYMJ).GT.0)) THEN

                     DO ISYMI = 1,ISYMJ-1

                        ISYMA = MULD2H(ISYMI,ISYMAI)

                        IF ((LVIRA(ISYMA).GT.0) .AND.
     &                      (NRHF(ISYMI) .GT.0)) THEN

                           DO J = 1,NRHF(ISYMJ)
                              DO LB = 1,LVIRB(ISYMB)

                                 LBJ = IX1AMB(ISYMB,ISYMJ)
     &                               + LVIRB(ISYMB)*(J - 1) + LB

                                 DO I = 1,NRHF(ISYMI)
                                    DO LA = 1,LVIRA(ISYMA)

                                       LAI   = IX1AMA(ISYMA,ISYMI)
     &                                       + LVIRA(ISYMA)*(I - 1)
     &                                       + LA
                                       LAIBJ = IX2SQ(ISYMAI,ISYMBJ)
     &                                       + NX1AMA(ISYMAI)*(LBJ - 1)
     &                                       + LAI

                                       X2NRM = X2NRM
     &                                       + X2AM(LAIBJ)*X2AM(LAIBJ)

                                    ENDDO
                                 ENDDO

                              ENDDO
                           ENDDO

                        ENDIF

                     ENDDO

                     ISYMI = ISYMJ
                     ISYMA = ISYMB

                     IF (LVIRA(ISYMA) .GT. 0) THEN

                        DO J = 1,NRHF(ISYMJ)
                           DO LB = 1,LVIRB(ISYMB)

                              BJ  = IT1AM(ISYMB,ISYMJ)
     &                            + NVIR(ISYMB)*(J - 1)
     &                            + IOFB1(ISYMB) + LB - 1
                              LBJ = IX1AMB(ISYMB,ISYMJ)
     &                            + LVIRB(ISYMB)*(J - 1) + LB

                              DO I = 1,J
                                 DO LA = 1,LVIRA(ISYMA)

                                    AI  = IT1AM(ISYMA,ISYMI)
     &                                  + NVIR(ISYMA)*(I - 1)
     &                                  + IOFA1(ISYMA) + LA - 1

                                    IF (AI .GT. BJ) GO TO 900

                                    LAI   = IX1AMA(ISYMA,ISYMI)
     &                                    + LVIRA(ISYMA)*(I - 1) 
     &                                    + LA
                                    LAIBJ = IX2SQ(ISYMAI,ISYMBJ)
     &                                    + NX1AMA(ISYMAI)*(LBJ - 1)
     &                                    + LAI

                                    X2NRM = X2NRM 
     &                                    + X2AM(LAIBJ)*X2AM(LAIBJ)

                                 ENDDO
                              ENDDO                                 

  900                         CONTINUE

                           ENDDO
                        ENDDO

                     ENDIF

                  ENDIF

               ENDDO

            ENDIF

 1000       CONTINUE

         ENDDO

      ELSE

         DO ISYMBJ = 1,NSYM

            ISYMAI = MULD2H(ISYMBJ,ISYMX)

            IF (ISYMAI .LT. ISYMBJ) THEN

               KOFFX = IX2SQ(ISYMAI,ISYMBJ) + 1
               NAIBJ = NX1AMA(ISYMAI)*NX1AMB(ISYMBJ)

               X2NRM = X2NRM + DDOT(NAIBJ,X2AM(KOFFX),1,X2AM(KOFFX),1)

            ENDIF

         ENDDO

      ENDIF

      RETURN
      END
C  /* Deck cc_cyifst */
      SUBROUTINE CC_CYIFST(IOFA1,LVIRA,ISYMAI,AI1)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Find global index of first ai (sym. ISYMAI) in #ai batch.
C
#include <implicit.h>
      INTEGER IOFA1(8), LVIRA(8)
      INTEGER AI1
#include <ccorb.h>
#include <ccsdsym.h>

      AI1 = 0

      DO ISYMI = 1,NSYM
         ISYMA = MULD2H(ISYMI,ISYMAI)
         IF ((NRHF(ISYMI).GT.0) .AND. (LVIRA(ISYMA).GT.0)) THEN
            AI1 = IT1AM(ISYMA,ISYMI) + IOFA1(ISYMA)
            RETURN
         ENDIF
      ENDDO

      RETURN
      END
C  /* Deck cc_cyilst */
      SUBROUTINE CC_CYILST(IOFA1,LVIRA,ISYMAI,AI2)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Find global index of last ai (sym. ISYMAI) in #ai batch.
C
#include <implicit.h>
      INTEGER IOFA1(8), LVIRA(8)
      INTEGER AI2
#include <ccorb.h>
#include <ccsdsym.h>

      AI2 = 0

      DO ISYMI = NSYM,1,-1
         ISYMA = MULD2H(ISYMI,ISYMAI)
         IF ((NRHF(ISYMI).GT.0) .AND. (LVIRA(ISYMA).GT.0)) THEN
            AI2 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(NRHF(ISYMI) - 1)
     &          + IOFA1(ISYMA) + LVIRA(ISYMA) - 1
            RETURN
         ENDIF
      ENDDO

      RETURN
      END
C  /* Deck cc_cyicni */
      SUBROUTINE CC_CYICNI(OMEGA1,FOCK,X2AM,WORK,LWORK,ISYMF,ISYMX,MEM)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate contribution I,
C
C              OMEGA1(#bj) = OMEGA1(#bj) + sum(#ai) X2AM(#ai,#bj) * FOCK(i#a)
C
C     Information about the virtual batches is taken from ciarc.h.
C
#include <implicit.h>
      DIMENSION OMEGA1(*), FOCK(*), X2AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ciarc.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYICNI')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER (TINY = 1.0D-14)

C     Set symmetries.
C     ---------------

      ISYMAI = ISYMF
      ISYMBJ = MULD2H(ISYMAI,ISYMX)

C     Is this the full matrix?
C     ------------------------

      XX2SQ = ONE*NX2SQ
      DIFF  = XT2SQ(ISYMX) - XX2SQ

C     Calculate contribution.
C     -----------------------

      MEM = 0

      IF (DABS(DIFF) .LE. TINY) THEN


         NTOAI = MAX(NT1AM(ISYMAI),1)
         KOFFX = IT2SQ(ISYMAI,ISYMBJ) + 1

         CALL DGEMV('T',NT1AM(ISYMAI),NT1AM(ISYMBJ),
     &              ONE,X2AM(KOFFX),NTOAI,FOCK,1,
     &              ONE,OMEGA1,1)

      ELSE

         IF ((NX1AMA(ISYMAI).GT.0) .AND. (NX1AMB(ISYMBJ).GT.0)) THEN

            KFOCK = 1
            KOMEG = KFOCK + NX1AMA(ISYMAI)
            KEND1 = KOMEG + NX1AMB(ISYMBJ)

            MEM   = KEND1 - 1
            LWRK1 = LWORK - MEM

            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need     : ',MEM,
     &         'Available: ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            CALL CC_CIACPS(FOCK,WORK(KFOCK),1,IOFA1,LVIRA,NX1AMA,
     &                     IX1AMA,ISYMAI)

            KOFFX = IX2SQ(ISYMAI,ISYMBJ) + 1
            CALL DGEMV('T',NX1AMA(ISYMAI),NX1AMB(ISYMBJ),
     &                 ONE,X2AM(KOFFX),NX1AMA(ISYMAI),WORK(KFOCK),1,
     &                 ZERO,WORK(KOMEG),1)

            CALL CC_CYIADD(OMEGA1,WORK(KOMEG),ONE,IOFB1,LVIRB,NX1AMB,
     &                     IX1AMB,ISYMBJ)

         ENDIF

      ENDIF

      RETURN
      END
C  /* Deck cc_cyiadd */
      SUBROUTINE CC_CYIADD(X1AM,S1AM,ALPHA,IOFA1,LVIRA,NX1AMA,IX1AMA,
     &                     ISYMAI)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Add sub-matrix ALPHA*S1AM(#ai) into full matrix X1AM(ai).
C
#include <implicit.h>
      DIMENSION X1AM(*), S1AM(*)
      INTEGER   IOFA1(8), LVIRA(8), NX1AMA(8), IX1AMA(8,8)
#include <ccorb.h>
#include <ccsdsym.h>

      DO ISYMI = 1,NSYM
         ISYMA = MULD2H(ISYMI,ISYMAI)
         IF (LVIRA(ISYMA) .EQ. NVIR(ISYMA)) THEN
            NUMAI = NVIR(ISYMA)*NRHF(ISYMI)
            KOFFS = IX1AMA(ISYMA,ISYMI) + 1
            KOFFX = IT1AM(ISYMA,ISYMI)  + 1
            CALL DAXPY(NUMAI,ALPHA,S1AM(KOFFS),1,X1AM(KOFFX),1)
         ELSE IF (LVIRA(ISYMA) .GT. 0) THEN
            DO I = 1,NRHF(ISYMI)
               KOFFS = IX1AMA(ISYMA,ISYMI) + LVIRA(ISYMA)*(I - 1)
     &               + 1
               KOFFX = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1)
     &               + IOFA1(ISYMA)
               CALL DAXPY(LVIRA(ISYMA),ALPHA,S1AM(KOFFS),1,
     &                                       X1AM(KOFFX),1)
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END
C  /* Deck cc_cyim */
      SUBROUTINE CC_CYIM(X2AM,WORK,LWORK,ISYMX,ITYP3,ISYCH3,ISIDE,MEM,
     &                   NCALL,NUMCHO)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate contributions to Y-type intermediates,
C
C              Y(#bj,J) = sum(#ai) X2AM(#ai,#bj) * L(#ai,J)
C
C              where the information about the virtual batches is taken
C              from ciarc.h and the symmetry of X2AM is ISYMX. ITYP3 specifies
C              the Cholesky vectors of symmetry ISYCH3 L(#bj,J) via subroutine
C              CHO_MOP.
C
C              The file information for storing the Y-type intermediates
C              is done through subroutines CC_CYIO*, using a long-integer
C              version of the crayio routines.
C
C              For NCALL <= 1, the Y intermediates are initialized on disk
C              (if needed, i.e. if not full a batch). Thus, it is implicitly
C              assumed that the b batch loop is the outermost one in the
C              calling routine!
C
C     N.B!!!!! NUMCHO array (# Cholesky vectors) is passed as an argument
C              and *not* taken from ccdeco.h !!!!!!!!!!
C
C     TODO: The I/O of the Y-intermediates is probably crappy....
C
#include <implicit.h>
      DIMENSION X2AM(*), WORK(LWORK)
      INTEGER NUMCHO(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ciarc.h>
#include <priunit.h>
#include <cyit.h>

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CC_CYIM')

      LOGICAL FULLA

      PARAMETER (IOPEN = -1, IKEEP = 1)
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Initialize MEM usage.
C     ---------------------

      MEM = 0

C     Get result symmetry.
C     --------------------

      ISYMY = MULD2H(ISYMX,ISYCH3)

C     If first call and if this is batched in the a index,
C     initialize Y intermediates on disk.
C     ----------------------------------------------------

      FULLA = .TRUE.
      DO ISYMA = 1,NSYM
         IF (LVIRA(ISYMA) .NE. NVIR(ISYMA)) THEN
            FULLA = .FALSE.
            GO TO 99
         ENDIF
      ENDDO
  99  CONTINUE

      IF ((NCALL .LE. 1) .AND. (.NOT. FULLA)) THEN
         CALL CC_CYINI(WORK,LWORK,ISYMY,ISIDE,NUMCHO)
      ENDIF

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYCH3)
         ISYMBJ = MULD2H(ISYMAI,ISYMX)

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (NX1AMA(ISYMAI) .LE. 0) GO TO 1000
         IF (NX1AMB(ISYMBJ) .LE. 0) GO TO 1000

C        Open MO file for this symmetry.
C        -------------------------------

         CALL CHO_MOP(IOPEN,ITYP3,ISYCHO,LUCHAI,1,ISYCH3)

C        Open Y intermediate file for this symmetry and ISIDE.
C        -----------------------------------------------------

         CALL CC_CYIOP(IOPEN,ISYCHO,ISIDE)

         IF (NX1AMA(ISYMAI) .LT. NT1AM(ISYMAI)) THEN

C           Allocate space for one full L vector.
C           -------------------------------------

            KCHOL = 1
            KEND1 = KCHOL + NT1AM(ISYMAI)
            LWRK1 = LWORK - KEND1 + 1

            IF (LWRK1 .LE. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need (more than): ',KEND1,
     &         'Available       : ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

         ELSE

C           Don't allocate anything yet.
C           ----------------------------

            KEND1 = 1
            LWRK1 = LWORK

         ENDIF

C        Set up Cholesky batch.
C        ----------------------

         MINMEM = NX1AMA(ISYMAI) + NX1AMB(ISYMBJ)
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (pref. multiple) : ',MINMEM,
     &      'Total memory available: ',LWORK
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

C           Complete allocation.
C           --------------------

            KRESY = KEND1
            KCHAI = KRESY + NX1AMB(ISYMBJ)*NUMV
            KLAST = KCHAI + NX1AMA(ISYMAI)*NUMV - 1

            MEM = MAX(MEM,KLAST)

C           Read the Y intermediates from disk,
C           or, if full a batch, initialize.
C           -----------------------------------

            IF (FULLA) THEN
               CALL DZERO(WORK(KRESY),NX1AMB(ISYMBJ)*NUMV)
            ELSE
               CALL CC_CYIRD(WORK(KRESY),IOFB1,LVIRB,NX1AMB,IX1AMB,
     &                       NUMV,JVEC1,ISYCHO,ISYMY,ISIDE)
            ENDIF

C           Read Cholesky vectors.
C           ----------------------

            IF (NX1AMA(ISYMAI) .LT. NT1AM(ISYMAI)) THEN
               DO IVEC = 1,NUMV
                  JVEC = JVEC1 + IVEC - 1
                  CALL CHO_MOREAD(WORK(KCHOL),NT1AM(ISYMAI),1,JVEC,
     &                            LUCHAI)
                  CALL CC_CIACPS(WORK(KCHOL),WORK(KCHAI),IVEC,IOFA1,
     &                           LVIRA,NX1AMA,IX1AMA,ISYMAI)
               ENDDO
            ELSE
               CALL CHO_MOREAD(WORK(KCHAI),NX1AMA(ISYMAI),NUMV,JVEC1,
     &                         LUCHAI)
            ENDIF

C           Calculate Y intermediate.
C           -------------------------

            DTIME = SECOND()

            NAI = NX1AMA(ISYMAI)
            NBJ = NX1AMB(ISYMBJ)

            KOFFX = IX2SQ(ISYMAI,ISYMBJ) + 1

            CALL DGEMM('T','N',NBJ,NUMV,NAI,
     &                 ONE,X2AM(KOFFX),NAI,WORK(KCHAI),NAI,
     &                 ONE,WORK(KRESY),NBJ)

            DTIME = SECOND() - DTIME
            TMCALC(2) = TMCALC(2) + DTIME

C           Save Y intermediates on disk.
C           -----------------------------

            CALL CC_CYIWR(WORK(KRESY),IOFB1,LVIRB,NX1AMB,IX1AMB,
     &                    NUMV,JVEC1,ISYCHO,ISYMY,ISIDE)

         ENDDO

C        Close Y intermediate file for this symmetry and ISIDE.
C        ------------------------------------------------------

         CALL CC_CYIOP(IKEEP,ISYCHO,ISIDE)

C        Close MO file for this symmetry.
C        --------------------------------

         CALL CHO_MOP(IKEEP,ITYP3,ISYCHO,LUCHAI,1,ISYCH3)

C        Escape point if nothing to do in this symmetry.
C        -----------------------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_cyienr */
      SUBROUTINE CC_CYIENR(WORK,LWORK,ECC2)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Calculate Y intermediate contribution
C              to CC2 energy,
C
C              ECC2 = ECC2 + sum(aiJ) Y(ai,J) * L(ia,J)
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIENR')

C     Start loop over Cholesky symmetries.
C     ------------------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (NT1AM(ISYCHO)  .LE. 0) GO TO 1000

C        Open files with Cholesky vectors and Y intermediates.
C        -----------------------------------------------------

         CALL CHO_MOP(-1,1,ISYCHO,LUCHIA,1,1)
         CALL CC_CYIOP(-1,ISYCHO,0)

C        Set up batch.
C        -------------

         MINMEM = 2*NT1AM(ISYCHO)
         NVEC   = MIN(LWORK/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,A,I1,A,/,5X,A,I10,/)')
     &      'Minimum memory required   : ',MINMEM,
     &      ' (Cholesky symmetry ',ISYCHO,')',
     &      'Total work space available: ',LWORK
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
            NAIJ  = NT1AM(ISYCHO)*NUMV 

            KYIM  = 1
            KCHIA = KYIM  + NAIJ

C           Read Cholesky vectors and Y intermediates.
C           ------------------------------------------

            CALL CC_CYIRDF(WORK(KYIM),NUMV,JVEC1,ISYCHO,1,0)
            CALL CHO_MOREAD(WORK(KCHIA),NT1AM(ISYCHO),NUMV,JVEC1,LUCHIA)

C           Calculate energy contribution.
C           ------------------------------

            ECC2 = ECC2 + DDOT(NAIJ,WORK(KYIM),1,WORK(KCHIA),1)

         ENDDO

C        Close files with Cholesky vectors and Y intermediates.
C        ------------------------------------------------------

         CALL CHO_MOP(1,1,ISYCHO,LUCHIA,1,1)
         CALL CC_CYIOP(1,ISYCHO,0)

C        Escape point for empty symmetry block.
C        --------------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_cyichksym */
      SUBROUTINE CC_CYICHKSYM(X2AM,ISYMX,ERRMAX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Check symmetry of doubles array.
C
#include <implicit.h>
      DIMENSION X2AM(*)
#include <ccorb.h>
#include <ccsdsym.h>

      INTEGER AI, BJ, AIBJ, BJAI

      ERRMAX = -1.0D10

      IF (ISYMX .EQ. 1) THEN

         DO ISYMBJ = 1,NSYM
            ISYMAI = ISYMBJ
            DO BJ = 1,NT1AM(ISYMBJ)
               DO AI = 1,BJ-1
                  AIBJ = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(BJ - 1)
     &                 + AI
                  BJAI = IT2SQ(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)*(AI - 1)
     &                 + BJ
                  DIFF   = DABS(X2AM(AIBJ) - X2AM(BJAI))
                  ERRMAX = MAX(ERRMAX,DIFF)
               ENDDO
            ENDDO
         ENDDO

      ELSE

         DO ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYMX)
            IF (ISYMAI .LT. ISYMBJ) THEN
               DO BJ = 1,NT1AM(ISYMBJ)
                  DO AI = 1,NT1AM(ISYMAI)
                     AIBJ = IT2SQ(ISYMAI,ISYMBJ)
     &                    + NT1AM(ISYMAI)*(BJ - 1) + AI
                     BJAI = IT2SQ(ISYMBJ,ISYMAI)
     &                    + NT1AM(ISYMBJ)*(AI - 1) + BJ
                     DIFF   = DABS(X2AM(AIBJ) - X2AM(BJAI))
                     ERRMAX = MAX(ERRMAX,DIFF)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

      ENDIF

      RETURN
      END
C  /* Deck cc_cyicnf */
      SUBROUTINE CC_CYICNF(FOCKD,WORK,LWORK,ISYMTR,ISYFIA,ISIDE,
     &                     ITYP3,ISYCH3,NUMYVC,ITYP4,ISYCH4,NUMYV2,
     &                     SCD,SECNAM,IERR)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Set up basic configuration (ciarc.h) for CC_CYI.
C
#include <implicit.h>
      DIMENSION FOCKD(*), WORK(LWORK)
      INTEGER NUMYVC(8), NUMYV2(8)
      CHARACTER*(*) SECNAM
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <chocc2.h>
#include <ciarc.h>

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Set error variable.
C     -------------------

      IERR = 0

C     Go configure.
C     -------------

      IF (ISIDE .EQ. -3) THEN

C        Ground state one-electron density matrix.
C        Make sure that the ground state amplitudes
C        are decomposed at this stage.
C        ------------------------------------------

         IF (ISYMTR .NE. 1) THEN
            IERR = 3
            RETURN
         ENDIF

         IF (ISYFIA .NE. 1) THEN
            IERR = 4
            RETURN
         ENDIF

         IF (.NOT. CHOT2C) THEN
            CALL CC_CHOCIM1(FOCKD,DUM1,DUM2,WORK,LWORK,0,0)
         ENDIF

         ITYP1  = 1
         ITYP2  = 8
         ITYP3  = 6
         ISYCH1 = 1
         ISYCH2 = 1
         ISYCH3 = 1
         DO ISYCHO = 1,NSYM
            NTOVEC(ISYCHO) = NUMCHO(ISYCHO)
            NUMYVC(ISYCHO) = NCHMOC(ISYCHO)
         ENDDO

         IPRCIA = IPRINT - 30 + IPRLVL

         LIDEN  = .FALSE.
         SYMTRZ = .TRUE.
         CIAMIO = .FALSE.
         GETMNM = .FALSE.
         INXINT = .FALSE.

         SCD = ONE

      ELSE IF (ISIDE .EQ. -2) THEN

C        Right-hand side for 0th order multiplier equations.
C        ---------------------------------------------------

         IF (ISYMTR .NE. 1) THEN
            IERR = 3
            RETURN
         ENDIF

         ISYCH1 = 1
         ISYCH2 = ISYCH1
         ISYCH3 = 1
         ITYP1  = 1
         ITYP2  = ITYP1
         ITYP3  = 3
         DO ISYCHO = 1,NSYM
            NTOVEC(ISYCHO) = NUMCHO(ISYCHO)
            NUMYVC(ISYCHO) = NUMCHO(ISYCHO)
         ENDDO

         IPRCIA = IPRINT - 30 + IPRLVL

         LIDEN  = .TRUE.
         SYMTRZ = .FALSE.
         CIAMIO = .FALSE.
         GETMNM = .FALSE.
         INXINT = .TRUE.

         SCD = ONE

      ELSE IF (ISIDE .EQ. -1) THEN

C        LH Jacobian transformation.
C        ---------------------------

         IF (ISYFIA .NE. 1) THEN
            IERR = 4
            RETURN
         ENDIF

         ISYCH1 = 1
         ISYCH2 = ISYMTR
         ISYCH3 = 1
         ITYP1  = 1
         ITYP2  = 8
         ITYP3  = 3
         DO ISYCHO = 1,NSYM
            NTOVEC(ISYCHO) = NUMCHO(ISYCHO)
            NUMYVC(ISYCHO) = NUMCHO(ISYCHO)
         ENDDO

         IPRCIA = IPRINT - 30 + IPRLVL

         LIDEN  = .FALSE.
         SYMTRZ = .TRUE.
         CIAMIO = .FALSE.
         GETMNM = .FALSE.
         INXINT = .TRUE.

         SCD = ONE

      ELSE IF (ISIDE .EQ. 0) THEN

C        Ground state vector function.
C        -----------------------------

         IF (ISYMTR .NE. 1) THEN
            IERR = 3
            RETURN
         ENDIF

         IF (ISYFIA .NE. 1) THEN
            IERR = 4
            RETURN
         ENDIF

         ISYCH1 = 1
         ISYCH2 = ISYCH1
         ISYCH3 = 1
         IF (CHOT2) THEN
            ITYP1 = 6
            DO ISYM = 1,NSYM
               NTOVEC(ISYM) = NCHMO2(ISYM)
               NUMYVC(ISYM) = NUMCHO(ISYM)
            ENDDO
         ELSE IF (CHOMO2) THEN
            ITYP1 = 5
            DO ISYM = 1,NSYM
               NTOVEC(ISYM) = NCHMO2(ISYM)
               NUMYVC(ISYM) = NUMCHO(ISYM)
            ENDDO
         ELSE
            ITYP1 = 3
            DO ISYM = 1,NSYM
               NTOVEC(ISYM) = NUMCHO(ISYM)
               NUMYVC(ISYM) = NUMCHO(ISYM)
            ENDDO
         ENDIF
         ITYP2 = ITYP1
         ITYP3 = 1

         IPRCIA = IPRINT - 30 + IPRLVL

         LIDEN  = .TRUE.
         SYMTRZ = .FALSE.
         CIAMIO = .FALSE.
         GETMNM = .FALSE.
         INXINT = .TRUE.

         SCD = ONE

      ELSE IF (ISIDE .EQ. 1) THEN

C        RH Jacobian transformation.
C        ---------------------------

         IF (ISYFIA .NE. 1) THEN
            IERR = 4
            RETURN
         ENDIF

         ISYCH1 = 1
         ISYCH2 = ISYMTR
         ISYCH3 = 1
         ITYP1  = 3
         ITYP2  = 9
         ITYP3  = 1
         DO ISYCHO = 1,NSYM
            NTOVEC(ISYCHO) = NUMCHO(ISYCHO)
            NUMYVC(ISYCHO) = NUMCHO(ISYCHO)
         ENDDO

         IPRCIA = IPRINT - 30 + IPRLVL

         LIDEN  = .FALSE.
         SYMTRZ = .TRUE.
         CIAMIO = .FALSE.
         GETMNM = .FALSE.
         INXINT = .TRUE.

         SCD = ONE

      ELSE IF (ISIDE .EQ. 2) THEN

C        F-matrix transformation.
C        ------------------------

         IERR = 2

      ELSE IF (ISIDE .EQ. 3) THEN

C        Right-hand side for 1st order amplitude equations.
C        The ground state amplitudes must have been decomposed
C        at this stage.
C        -----------------------------------------------------

         IF (ISYFIA .NE. 1) THEN
            IERR = 4
            RETURN
         ENDIF

         IF (.NOT. CHOT2C) THEN
            IERR = 5
            RETURN
         ENDIF

         ISYCH1 = 1
         ISYCH2 = ISYMTR
         ISYCH3 = 1
         ITYP1  = 6
         ITYP2  = 9
         ITYP3  = 1
         DO ISYCHO = 1,NSYM
            NTOVEC(ISYCHO) = NCHMOC(ISYCHO)
            NUMYVC(ISYCHO) = NUMCHO(ISYCHO)
         ENDDO

         IPRCIA = IPRINT - 30 + IPRLVL

         LIDEN  = .FALSE.
         SYMTRZ = .TRUE.
         CIAMIO = .FALSE.
         GETMNM = .FALSE.
         INXINT = .TRUE.

         SCD = ONE

      ELSE

         IERR = 1

      ENDIF

      RETURN
      END
C  /* Deck cc_cyiajb */
      SUBROUTINE CC_CYIAJB(XINT,WORK,LWORK,FAC,XINI,KLAST,TIMI,NPASS,
     &                     NCALL)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate FAC*(i#a|j#b). The original configuration
C              file ciarc.h is restored on exit (the batch info is
C              unchanged, of course, in this routine). The result
C              is placed in XINT if XINI = .TRUE., else the result
C              is added to XINT.
C
#include <implicit.h>
      DIMENSION XINT(*), WORK(LWORK)
      LOGICAL XINI
#include <ciarc.h>
#include <cyit.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIAJB')

      LOGICAL LIDSAV, SMTSAV, MIOSAV, INXSAV, GTMSAV

      PARAMETER (SMALL = 1.0D-4, ONE = 1.0D0)

C     Save ciarc.h control variables.
C     -------------------------------

      IT1SAV = ITYP1
      IT2SAV = ITYP2
      IS1SAV = ISYCH1
      IS2SAV = ISYCH2
      LIDSAV = LIDEN
      SMTSAV = SYMTRZ
      MIOSAV = CIAMIO
      INXSAV = INXINT
      GTMSAV = GETMNM

C     Configure CC_CIA for (i#a|j#b) calculation.
C     -------------------------------------------

      ITYP1  = 1
      ITYP2  = ITYP1
      ISYCH1 = 1
      ISYCH2 = ISYCH1
      LIDEN  = .TRUE.
      SYMTRZ = .FALSE.
      CIAMIO = .FALSE.
      INXINT = XINI
      GETMNM = .FALSE.
      SCD    = ONE

C     Calculate integrals and scale.
C     ------------------------------

      IF (DABS(FAC) .LE. SMALL) THEN
         WRITE(LUPRI,'(/,5X,A,A,/,5X,A,1P,D22.15,/,5X,A,/)')
     &   '*** WARNING: Small integral scaling factor in ',SECNAM,
     &   '             FAC = ',FAC,
     &   '             Program continues nevertheless...'
      ENDIF

      IF (.NOT. XINI) THEN
         IF (FAC .NE. ONE) THEN
            FACI = ONE/FAC
            CALL DSCAL(NX2SQ,FACI,XINT,1)
         ENDIF
      ENDIF

      TIMI2 = SECOND()
      CALL CC_CIA(XINT,WORK,LWORK,SCD,KLAST,NPASS)
      TIMI2 = SECOND() - TIMI2

      IF (FAC .NE. ONE) THEN
         CALL DSCAL(NX2SQ,FAC,XINT,1)
      ENDIF

C     Update info.
C     ------------

      NCALL = NCALL + 1
      TIMI  = TIMI  + TIMI2

      TIMCYI(2) = TIMCYI(2) + TIMI2
      TMCALC(1) = TMCALC(1) + TIMCIA(4)

C     Restore ciarc.h
C     ---------------

      ITYP1  = IT1SAV
      ITYP2  = IT2SAV
      ISYCH1 = IS1SAV
      ISYCH2 = IS2SAV
      LIDEN  = LIDSAV
      SYMTRZ = SMTSAV
      CIAMIO = MIOSAV
      INXINT = INXSAV
      GETMNM = GTMSAV

      RETURN
      END
C  /* Deck cc_cyiold */
      SUBROUTINE CC_CYIOLD(FOCKD,OM1,X1AM,FIA,WORK,LWORK,
     &                     ISIDE,ISYMTR,ISYFIA,FREQ,X2NRM,X2CNM,DELCH2)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     OBSOLETE ROUTINE - KEPT ONLY FOR DEBUGGING PURPOSES!
C
C     Purpose: Calculate generalized Y intermediates
C
C              Y(ai,J) = sum(bj) [2*X(ai,bj)-X(aj,bi)] * Z(bj,J)
C
C              for use in Cholesky CC2 vector function calculation or
C              linear transformations without permanently storing doubles.
C              Only characteristic of the X amplitudes returned are the
C              packed (ai <= bj) norm of X [in X2NRM] and that of the
C              2CME array [in X2CNM].
C
C     Obligatory input:
C
C        FOCKD : Canonical orbital energies (CC ordering).
C        WORK  : Work array, dimension LWORK.
C        ISYMTR: Symmetry of trial vector.
C        ISIDE : Transformation type (see below) and specifier for disk
C                storage of the Y intermediates (see CC_CYIOP).
C
C     Additional input, depending on ISIDE:
C
C        ISIDE = -3 : Y intermediates are calculated for constructing the
C                     the CC2 ground state one-electron density matrix.
C                     The Y intermediates are returned on disk as specified
C                     in chocc2lhs.h. [I.e. using the same files as the Y
C                     intermediates in the left-hand side linear
C                     transformations].
C                     Additional input:
C                     OM1   : Dummy.
C                     X1AM  : Singles part of multipliers.
C                     FIA   : F(ia) matrix.
C                     ISYFIA: Symmetry of F(ia).
C                     FREQ  : Dummy (frequency is taken as zero).
C                     DELCH2: Logical flag for deleting files with perturbed
C                             Cholesky vectors.
C
C        ISIDE = -2 : Y intermediates are calculated for constructing the
C                     effective right-hand side for the 0th order multiplier
C                     equations. The Y intermediates are returned on disk
C                     as specified in chocc2lhs.h [i.e. using the same files
C                     as the Y intermediates in the left-hand side linear
C                     transformations].
C                     Additional input:
C                     OM1   : Dummy.
C                     X1AM  : Dummy.
C                     FIA   : Dummy.
C                     ISYFIA: Dummy.
C                     FREQ  : Dummy (frequency is taken as zero).
C                     DELCH2: Logical, ignored.
C
C        ISIDE = -1 : Y intermediates are calculated for left-hand side
C                     linear Jacobian transformations (incl. 0th order
C                     multipliers as a special case). The Y intermediates
C                     are returned on disk as specified in chocc2lhs.h.
C                     Additional input:
C                     OM1   : Dummy.
C                     X1AM  : Singles part of trial vector.
C                     FIA   : F(ia) matrix.
C                     ISYFIA: Symmetry of F(ia).
C                     FREQ  : Perturbation frequency.
C                     DELCH2: Logical flag for deleting files with perturbed
C                             Cholesky vectors, i.e. the ITYP2 files.
C
C        ISIDE =  0 : Y intermediates are calculated for the ground state
C                     CC2 vector function. Additionally, the I-term is
C                     evaluated. The Y intermediates are returned on disk
C                     as specified in chocc2.h.
C                     Additional input:
C                     OM1   : Result array (omega1) to which the I-term is
C                             added.
C                     X1AM  : Dummy.
C                     FIA   : F(ia) matrix.
C                     ISYFIA: Symmetry of F(ia).
C                     FREQ  : Dummy (frequency is taken as zero).
C                     DELCH2: Logical, ignored.
C
C        ISIDE =  1 : Y intermediates are calculated for right-hand side
C                     linear Jacobian transformations. Additionally, part
C                     of the I-term is calculated. The Y intermediates
C                     are returned on disk as specified in chocc2rhs.h.
C                     Additional input:
C                     OM1   : Result array (sigma1) to which the partial
C                             I-term is added.
C                     X1AM  : Dummy.
C                     FIA   : F(ia) matrix.
C                     ISYFIA: Symmetry of F(ia).
C                     FREQ  : Perturbation frequency.
C                     DELCH2: Logical flag for deleting files with perturbed
C                             Cholesky vectors, i.e. the ITYP2 files.
C
C        ISIDE =  2 : Not implemented.
C
C        ISIDE =  3 : Y intermediates are calculated for constructing the
C                     effective right-hand side for the 1st order amplitude
C                     equations. The Y intermediates are returned on disk
C                     as specified in chocc2rhs.h [i.e. using the same files
C                     as the Y intermediates in the right-hand side linear
C                     transformations].
C                     Additional input:
C                     OM1   : Result array (xi-X) to which the I-type term
C                             is added.
C                     X1AM  : Dummy.
C                     FIA   : F(ia) matrix.
C                     ISYFIA: Symmetry of F(ia).
C                     FREQ  : Perturbation frequency.
C                     DELCH2: Logical flag for deleting files with perturbed
C                             Cholesky vectors, i.e. the ITYP2 files.
C
C     Other assumptions:
C
C        - The relevant MO Cholesky vectors must be available on disk.
C          Routine CHO_MOP is used for opening the relevant files.
C
#include <implicit.h>
      DIMENSION FOCKD(*), OM1(*), X1AM(*), FIA(*)
      DIMENSION WORK(LWORK)
      LOGICAL DELCH2
#include <maxorb.h>
#include <ccisvi.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <ccdeco.h>
#include <ciarc.h>
#include <chocc2.h>
#include <priunit.h>
#include <cyit.h>

      INTEGER ABEG, AEND, BBEG, BEND
      INTEGER LIBJ(8), NUMYVC(8), NUMYV2(8)

      LOGICAL LOCDBG, CALYI2
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIOLD')

      PARAMETER (ISIMN = -3, ISIMX = 3)
      CHARACTER*46 ISITXT(ISIMN:ISIMX)

      PARAMETER (INFO = 3)
      PARAMETER (IOPEN = -1, IDEL = 0, IKEEP = 1)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, HALF = 0.5D0)
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.0D0)
      PARAMETER (TINY = 1.0D-14)

      DATA ISITXT / 'ground state one-electron density matrix      ',
     &              'right-hand side for 0th order multiplier eqs. ',
     &              'left-hand side linear Jacobian transformation ',
     &              'ground state vector function                  ',
     &              'right-hand side linear Jacobian transformation',
     &              'F-matrix transformation                       ',
     &              'right-hand side for 1st order amplitude eqs.  '/

C     Initialize timings.
C     -------------------

      TIMCYI(1) = SECOND()
      DO ITMCYI = 2,NTMCYI
         TIMCYI(ITMCYI) = ZERO
      ENDDO
      DO ITCYIO = 1,NTCYIO
         TMCYIO(ITCYIO) = ZERO
      ENDDO
      DO ITCALC = 1,NTCALC
         TMCALC(ITCALC) = ZERO
      ENDDO

C     Some debug stuff.
C     -----------------

      IF (LOCDBG) THEN
         X1NRM = DSQRT(DDOT(NT1AM(ISYMTR),X1AM,1,X1AM,1))
         FIANM = DSQRT(DDOT(NT1AM(ISYFIA),FIA,1,FIA,1))
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': norm of X1AM: ',X1NRM
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': norm of FIA : ',FIANM
      ENDIF

C     Overall configuration of the integral assembler CC_CIA
C     and subsequent Y intermediate calculation.
C     ------------------------------------------------------

      CALL CC_CYICNF(FOCKD,WORK,LWORK,ISYMTR,ISYFIA,ISIDE,ITYP3,ISYCH3,
     &               NUMYVC,ITYP4,ISYCH4,NUMYV2,SCD,SECNAM,IERR)

      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,I3,/,5X,A,I3)')
     &   'Configuration error in ',SECNAM,': CC_CYICNF returned IERR =',
     &   IERR,
     &   'Option ISIDE is ',ISIDE
         IF ((ISIDE.GE.ISIMN) .AND. (ISIDE.LE.ISIMX)) THEN
            WRITE(LUPRI,'(5X,A,/,5X,A,A)')
     &      'Apparant purpose of invocation:',
     &      'Y intermediates for ',ISITXT(ISIDE)
         ENDIF
         WRITE(LUPRI,'(5X,A,/)') 'Diagnostic of error:'
         IF (IERR .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A,/)')
     &      'Unknown ISIDE option'
            WRITE(LUPRI,'(5X,A,/,5X,A,/)')
     &      'Known ISIDE options:',
     &      '===================='
            WRITE(LUPRI,'(5X,A,/,5X,A)')
     &      'ISIDE                   Description',
     &      '-----------------------------------------------------'
            DO JSIDE = ISIMN,ISIMX
               WRITE(LUPRI,'(5X,I3,4X,A)') JSIDE,ISITXT(JSIDE)
            ENDDO
            WRITE(LUPRI,'(5X,A,/)')
     &      '-----------------------------------------------------'
         ELSE IF (IERR .EQ. 2) THEN
            WRITE(LUPRI,'(5X,A,/)')
     &      'Option recognized but not implemented yet'
         ELSE IF (IERR .EQ. 3) THEN
            WRITE(LUPRI,'(5X,A,I2,/)')
     &      'Symmetry error: ISYMTR =',ISYMTR
         ELSE IF (IERR .EQ. 4) THEN
            WRITE(LUPRI,'(5X,A,I2,/)')
     &      'Symmetry error: ISYFIA =',ISYFIA
         ELSE IF (IERR .EQ. 5) THEN
            WRITE(LUPRI,'(5X,A,/)')
     &      'Ground state amplitudes must be decomposed on entry'
         ELSE
            WRITE(LUPRI,'(5X,A,/)')
     &      'No diagnostic available!'
         ENDIF
         CALL QUIT('Configuration error in '//SECNAM)
      ENDIF

C     Initial check of memory.
C     ------------------------

      GETMNM = .TRUE.
      CALL CC_CIA(DUM1,DUM2,IDUM1,DUM3,NEED,IDUM2)
      GETMNM = .FALSE.

      IF (LWORK .LE. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': Work space non-positive!'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &   'Min. memory required: ',NEED,
     &   'Work space supplied : ',LWORK
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
      MXMEML = 0
      NPASS  = 0
      NCALL  = 0
      NCALLY = 0
      NCALY2 = 0
      ISYINT = MULD2H(ISYCH1,ISYCH2)
      XLWORK = ONE*LWORK
      X2NRM  = ZERO
      X2CNM  = ZERO

C     Determine memory split.
C     -----------------------

      CALL CC_CYIMSP(LWORK,LWRK)
      LEFT = LWORK - LWRK

C     Check.
C     ------

      IF ((LWRK.LE.0) .OR. (LEFT.LE.0)) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &   'Total memory available                  : ',LWORK,
     &   'Memory reserved for integrals/amplitudes: ',LWRK,
     &   'Memory reserved for Cholesky vectors    : ',LEFT
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN

         CALL AROUND('Output from '//SECNAM)

         WRITE(LUPRI,'(2X,A,A)')
     &   'Y intermediates for ',ISITXT(ISIDE)
         WRITE(LUPRI,'(2X,A,I1,1X,1P,D22.15)')
     &   'Symmetry and dimension of integrals/amplitudes: ',
     &   ISYINT,XT2SQ(ISYINT)

         XMB  = XLWORK*WTOMB
         XPC1 = (D100*LWRK)/XLWORK
         XPC2 = (D100*LEFT)/XLWORK
         WRITE(LUPRI,'(2X,A,I10,A,F7.1,A)')
     &   'Total memory available: ',LWORK,' words (',XMB,' Mb).'
         WRITE(LUPRI,'(2X,A,I10,A,F7.1,A,/,2X,A,I10,A,F7.1,A,/)')
     &   'Maximum mem. allowed for integrals/amplitudes: ',LWRK,
     &   ' words (',XPC1,'%).',
     &   'Remaining memory used for Cholesky vectors   : ',LEFT,
     &   ' words (',XPC2,'%).'

      ENDIF

C     Start batch loop over b.
C     ------------------------

      IBATB = 0
      BBEG  = 1
      NUMB  = 0
  100 CONTINUE

         BBEG = BBEG + NUMB
         IF (BBEG .GT. NVIRT) GO TO 200

C        Time this batch.
C        ----------------

         TBBAT = SECOND()

C        Update batch counter.
C        ---------------------

         IBATB = IBATB + 1

C        Initialize min. and max. Cholesky batch variables.
C        --------------------------------------------------

         MNCHBT = 230000000
         MXCHBT = -230000000

C        Find out how many b's can be treated, estimating
C        the number of a's along the way.
C        ------------------------------------------------

         CALL CC_CYIBTB(LVIRB,BBEG,NUMB,LWRK,ISYINT)

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

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(5X,A,I6,A,/,5X,A)')
     &      'b-batch number ',IBATB,':',
     &      '======================'
            WRITE(LUPRI,'(5X,A,I6)')
     &      '#b   : ',NUMB
            WRITE(LUPRI,'(5X,A,I6,A,I1,/,5X,A,I6,A,I1,/)')
     &      'First: ',IOFB1(ISBBEG),'   symmetry: ',ISBBEG,
     &      'Last : ',IOFB1(ISBEND)+LVIRB(ISBEND)-1,'   symmetry: ',
     &                                                  ISBEND
            WRITE(LUPRI,'(5X,A,A,/,5X,A,A)')
     &      '    #a   First     Last    Memory             CPU Time',
     &      ' (Seconds)',
     &      '       (a,sym.a) (a,sym.a) Usage (%)',
     &      '  Integrals  Y-interm.      Total'
            WRITE(LUPRI,'(5X,A,A)')
     &      '------------------------------------------------------',
     &      '---------------'
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

            CALL CC_CYIBTA(LIBJ,NX1AMB,LVIRA,ABEG,NUMA,MEM,LWRK,ISYINT)

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

            CALL IZERO(IX2SQ,64)
            ICOUNT = 0
            DO ISYMBJ = 1,NSYM
               ISYMAI = MULD2H(ISYMBJ,ISYINT)
               IX2SQ(ISYMAI,ISYMBJ) = ICOUNT
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

C           Allocation for integral intermediates.
C           --------------------------------------

            KXINT = 1
            KEND1 = KXINT + NX2SQ
            LWRK1 = LWORK - KEND1 + 1

C           2*(ia|jb) unperturbed integrals must be
C           added to obtain ground state multipliers.
C           -----------------------------------------

            IF (ISIDE .EQ. -3) THEN
               CALL CC_CYIAJB(WORK(KXINT),WORK(KEND1),LWRK1,TWO,.TRUE.,
     &                        KLAST,TIMI,NPASS,NCALL)
               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)
               DO ISYCHO = 1,NSYM
                  MNCHBT = MIN(MNCHBT,NCHBAT(ISYCHO))
                  MXCHBT = MAX(MXCHBT,NCHBAT(ISYCHO))
               ENDDO
               IF (LOCDBG) THEN
                  IF (NX2SQ .EQ. NT2SQ(ISYINT)) THEN
                     CALL CC_CYICHKSYM(WORK(KXINT),ISYINT,ERRMAX)
                     WRITE(LUPRI,'(A,I3,A,1P,D22.15)')
     &               'Max. sym. error after call',NCALL,
     &               ' to CC_CIA: ',ERRMAX
                  ENDIF
               ENDIF
            ENDIF

C           Calculate integrals.
C           --------------------

            TIMI = SECOND()
            CALL CC_CIA(WORK(KXINT),WORK(KEND1),LWRK1,SCD,KLAST,NPASS)
            TIMI = SECOND() - TIMI
            TIMCYI(2) = TIMCYI(2) + TIMI
            TMCALC(1) = TMCALC(1) + TIMCIA(4)

            NCALL  = NCALL + 1
            MUSED  = KEND1 - 1 + KLAST
            MXMEML = MAX(MXMEML,MUSED)

            DO ISYCHO = 1,NSYM
               MNCHBT = MIN(MNCHBT,NCHBAT(ISYCHO))
               MXCHBT = MAX(MXCHBT,NCHBAT(ISYCHO))
            ENDDO

            IF (LOCDBG) THEN
               IF (NX2SQ .EQ. NT2SQ(ISYINT)) THEN
                  CALL CC_CYICHKSYM(WORK(KXINT),ISYINT,ERRMAX)
                  WRITE(LUPRI,'(A,I3,A,1P,D22.15)')
     &            'Max. sym. error after call',NCALL,
     &            ' to CC_CIA: ',ERRMAX
               ENDIF
            ENDIF

C           Add extra terms,
C           X(#ai,#bj) = X(#ai,#bj) + P(ai,bj) F(ia) X1AM(bj).
C           --------------------------------------------------

            IF ((ISIDE.EQ.-1) .OR. (ISIDE.EQ.-3)) THEN
               DTIME = SECOND()
               CALL CC_CYILHX(X1AM,WORK(KXINT),FIA,
     &                        WORK(KEND1),LWRK1,KLAST,ISYMTR,ISYFIA,
     &                        ONE)
               DTIME = SECOND() - DTIME
               TIMCYI(3) = TIMCYI(3) + DTIME
               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)
               IF (LOCDBG) THEN
                  IF (NX2SQ .EQ. NT2SQ(ISYINT)) THEN
                     CALL CC_CYICHKSYM(WORK(KXINT),ISYINT,ERRMAX)
                     WRITE(LUPRI,'(A,1P,D22.15)')
     &               'Max. sym. error after call  1 to CC_CYILHX: ',
     &               ERRMAX
                  ENDIF
               ENDIF
            ENDIF

C           Divide by orbital energies to obtain amplitudes.
C           For left-hand transformations the sign convention
C           for frequency is different, but it is assumed
C           that the sign of FREQ has been taken care of in
C           calling routine! 
C           If ground state amplitudes were decomposed prior
C           to this routine, the denominators are already
C           included and we simply need to scale by -1 to 
C           get the right sign.
C           -------------------------------------------------

            IF ((ISIDE.EQ.-3) .OR. (ISIDE.EQ.-2)) THEN
               DTIME = SECOND()
               CALL CC_DNOM(FOCKD,WORK(KXINT),ZERO,ISYINT)
               DTIME = SECOND() - DTIME
               TIMCYI(4) = TIMCYI(4) + DTIME
c           ELSE IF (ISIDE .EQ. -1) THEN
c              FREQN = -FREQ
c              DTIME = SECOND()
c              CALL CC_DNOM(FOCKD,WORK(KXINT),FREQN,ISYINT)
c              DTIME = SECOND() - DTIME
c              TIMCYI(4) = TIMCYI(4) + DTIME
            ELSE IF (ISIDE .EQ. 0) THEN
               IF (CHOT2) THEN
                  CALL DSCAL(NX2SQ,XMONE,WORK(KXINT),1)
               ELSE
                  DTIME = SECOND()
                  CALL CC_DNOM(FOCKD,WORK(KXINT),ZERO,ISYINT)
                  DTIME = SECOND() - DTIME
                  TIMCYI(4) = TIMCYI(4) + DTIME
               ENDIF
            ELSE
               DTIME = SECOND()
               CALL CC_DNOM(FOCKD,WORK(KXINT),FREQ,ISYINT)
               DTIME = SECOND() - DTIME
               TIMCYI(4) = TIMCYI(4) + DTIME
            ENDIF

C           Calculate contribution to (packed) doubles norm.
C           ------------------------------------------------

            DTIME = SECOND()
            CALL CC_CYINRM(WORK(KXINT),ISYINT,X2NRM)
            DTIME = SECOND() - DTIME
            TIMCYI(5) = TIMCYI(5) + DTIME

C           Set up 2 Coulomb minus exchange (almost) in place.
C           --------------------------------------------------

            DTIME = SECOND()
            CALL CC_CYITCME(WORK(KXINT),WORK(KEND1),LWRK1,ISYINT,KLAST)
            DTIME = SECOND() - DTIME
            TIMCYI(6) = TIMCYI(6) + DTIME

            MUSED  = KEND1 - 1 + KLAST
            MXMEML = MAX(MXMEML,MUSED)

C           Calculate contribution to (packed) 2CME norm.
C           ---------------------------------------------

            DTIME = SECOND()
            CALL CC_CYINRM(WORK(KXINT),ISYINT,X2CNM)
            DTIME = SECOND() - DTIME
            TIMCYI(5) = TIMCYI(5) + DTIME

C           Calculate I-type terms.
C           -----------------------

            IF ((ISIDE.EQ.0) .OR. (ISIDE.EQ.1) .OR. (ISIDE.EQ.3)) THEN
               DTIME = SECOND()
               CALL CC_CYICNI(OM1,FIA,WORK(KXINT),WORK(KEND1),LWRK1,
     &                        ISYFIA,ISYINT,KLAST)
               DTIME = SECOND() - DTIME
               TIMCYI(7) = TIMCYI(7) + DTIME
               MUSED  = KEND1 - 1 + KLAST
               MXMEML = MAX(MXMEML,MUSED)
            ENDIF

C           Calculate Y intermediates.
C           --------------------------

            NCALLY = NCALLY + 1

            TIMY = SECOND()
            CALL CC_CYIM(WORK(KXINT),WORK(KEND1),LWRK1,ISYINT,ITYP3,
     &                   ISYCH3,ISIDE,KLAST,NCALLY,NUMYVC)
            TIMY = SECOND() - TIMY
            TIMCYI(8) = TIMCYI(8) + TIMY

            MUSED  = KEND1 - 1 + KLAST
            MXMEML = MAX(MXMEML,MUSED)

C           Print information from this a batch.
C           ------------------------------------

            IF (IPRINT .GE. INFO) THEN
               TABAT = SECOND() - TABAT
               XMUSE = (D100*MXMEML)/XLWORK
               LAFST = IOFA1(ISABEG)
               LALST = IOFA1(ISAEND) + LVIRA(ISAEND) - 1
               WRITE(LUPRI,1) NUMA,LAFST,ISABEG,LALST,ISAEND,XMUSE,
     &                        TIMI,TIMY,TABAT
    1          FORMAT(5X,I6,1X,I6,1X,I1,2X,I6,1X,I1,2X,F7.2,3X,F10.2,1X,
     &                F10.2,1X,F10.2)
            ENDIF

C           Update memory info.
C           -------------------

            MXMEMT = MAX(MXMEMT,MXMEML)
            MXMEML = 0

C           Go to next a batch.
C           (Which will exit immediately if no more a's).
C           ---------------------------------------------

            GO TO 110

C        Exit point for a batch loop.
C        ----------------------------

  190    CONTINUE

C        Print info for this b-batch and got to next.
C        (Which will exit immediately if no more b's.)
C        ---------------------------------------------

         IF (IPRINT .GE. INFO) THEN
            WRITE(LUPRI,'(5X,A,A,/)')
     &      '------------------------------------------------------',
     &      '---------------'
            TBBAT = SECOND() - TBBAT
            WRITE(LUPRI,'(5X,A,F10.2,A)')
     &      'Total time for this b-batch: ',TBBAT,' seconds'
            WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &      'Accumulated amplitude norm : ',DSQRT(X2NRM),
     &      'Accumulated 2CME am.  norm : ',DSQRT(X2CNM)
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum number of Cholesky batches: ',MNCHBT,
     &      'Maximum number of Cholesky batches: ',MXCHBT
            CALL FLSHFO(LUPRI)
         ENDIF

         GO TO 100

C     Exit point for b batch loop.
C     ----------------------------

  200 CONTINUE

C     If requested, delete files with perturbed Cholesky vectors.
C     -----------------------------------------------------------

      IF ((ISIDE.EQ.-3) .OR. (ISIDE.EQ.-1) .OR. (ISIDE.EQ.1) .OR.
     &    (ISIDE.EQ.3)) THEN
         IF (DELCH2) THEN
            DO ISYCHO = 1,NSYM
               CALL CHO_MOP(IOPEN,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)
               CALL CHO_MOP(IDEL,ITYP2,ISYCHO,LUCHO2,1,ISYCH2)
            ENDDO
            IF (LOCDBG) THEN
               WRITE(LUPRI,'(A,A)')
     &         SECNAM,': perturbed Cholesky vector files deleted...'
            ENDIF
         ENDIF
      ENDIF

C     Print.
C     ------

      IF (IPRINT .GE. INFO) THEN

         TIMCYI(1) = SECOND() - TIMCYI(1)

         CALL HEADER(SECNAM//' Timing by Routine',-1)
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CIA    : ',TIMCYI(2),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYILHX : ',TIMCYI(3),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_DNOM   : ',TIMCYI(4),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYINRM : ',TIMCYI(5),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYITCME: ',TIMCYI(6),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYICNI : ',TIMCYI(7),' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used in CC_CYIM   : ',TIMCYI(8),' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time             : ',TIMCYI(1),' seconds'

         XTTOT = TMCYIO(1) + TMCYIO(2) + TMCALC(1) + TMCALC(2)
         XI1 = D100*TMCYIO(1)/TIMCYI(1)
         XI2 = D100*TMCYIO(2)/TIMCYI(1)
         XC1 = D100*TMCALC(1)/TIMCYI(1)
         XC2 = D100*TMCALC(2)/TIMCYI(1)
         XT  = D100*XTTOT/TIMCYI(1)
         CALL HEADER(SECNAM//' Timing by Task',-1)
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used for MO Cholesky I/O         : ',TMCYIO(1),
     &   ' seconds (',XI1,' %)'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used for Y intermediate I/O      : ',TMCYIO(2),
     &   ' seconds (',XI2,' %)'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used in DGEMM for integrals      : ',TMCALC(1),
     &   ' seconds (',XC1,' %)'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Time used in DGEMM for Y intermediates: ',TMCALC(2),
     &   ' seconds (',XC2,' %)'
         WRITE(LUPRI,'(5X,A,A)')
     &   '---------------------------------------------------',
     &   '-------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,F7.2,A)')
     &   'Total                                 : ',XTTOT,
     &   ' seconds (',XT,' %)'

         CALL FLSHFO(LUPRI)

      ENDIF

      RETURN
      END
