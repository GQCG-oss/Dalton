C  /* Deck cc_caa */
      SUBROUTINE CC_CAA(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &                  UVEC,ISYMU,VVEC,ISYMV,FAC2,NUMUV,
     &                  SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &                  X2AM,WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2,MEMUSE)
C
C     Coupled Cluster Cholesky Amplitude Assembler.
C     =============================================
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate a batch of amplitudes according according to
C
C     X(#ai,#bj) = (2-Pij) 1/[FREQ - e(#ai,#bj)] P(#ai,#bj)
C                *[ FAC1 * sum(I=1,NUMP12)
C                 * sum(k=1,NSYM) sum(J=1,NCHP12(k,I)) L1(#ai,J,I) * L2(#bj,J,I)
C                 + FAC2 * sum(I=1,NUMUV) UVEC(#ai,I) * VVEC(#bj,I)
C                 ]
C
C     (2-Pij) generates 2 Coulomb minus exchange, e(ai,bj) is the
C     usual orbital energy denominator, and P(ai,bj) symmetrizes in
C     the ai,bj indices. Information about the batches of virtual orbitals
C     is taken from ciarc.h (see CC_CIA for details).
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
C     X2AM  : Array containing the requested amplitude batch.
C     X2NRM : Norm of packed amplitude array, before 2CME.
C     X2CNM : Norm of packed amplitude array, after  2CME (if IOPTCE = 1).
C     MEMUSE: Memory used in this routine.
C
C     Notes:
C     ------
C
C     - The Cholesky vectors must be available on disk, through CHO_MOP.
C
#include <implicit.h>
      DIMENSION UVEC(*), VVEC(*), FOCKD(*)
      DIMENSION SCD(NUMP12), X2AM(*), WORK(LWORK)
      INTEGER JTYP1(NUMP12), ISYMP1(NUMP12)
      INTEGER JTYP2(NUMP12), ISYMP2(NUMP12)
      INTEGER NCHP12(8,NUMP12)
      INTEGER ISYMU(NUMUV), ISYMV(NUMUV)
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
#include <dummy.h>

      CHARACTER*6 SECNAM
      PARAMETER (SECNAM = 'CC_CAA')

      LOGICAL LOCDBG, FULLIN
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (INFO = 20)
      PARAMETER (IOPEN = -1, IDEL = 0, IKEEP = 1)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, HALF = 0.5D0)
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (WTOMB = 7.62939453125D-6, D100 = 100.0D0)
      PARAMETER (SMALL = 1.0D-4, TINY = 1.0D-14)

C     Start timing.
C     -------------

      TIMT  = SECOND()
      TIMI  = ZERO
      TIMUV = ZERO

C     Return if nothing to do.
C     ------------------------

      IF ((NUMP12.LE.0) .AND. (NUMUV.LE.0)) THEN
         MEMUSE = 0
         RETURN
      ELSE
         IF (NX2SQ .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,I10,/)')
     &      SECNAM,': ERROR: NX2SQ = ',NX2SQ
            CALL QUIT('Error in '//SECNAM)
         ENDIF
      ENDIF

C     Debug print.
C     ------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) 'IOFA1:'
         WRITE(LUPRI,'(8I10)') (IOFA1(I),I=1,NSYM)
         WRITE(LUPRI,*) 'LVIRA:'
         WRITE(LUPRI,'(8I10)') (LVIRA(I),I=1,NSYM)
         WRITE(LUPRI,*) 'IOFB1:'
         WRITE(LUPRI,'(8I10)') (IOFB1(I),I=1,NSYM)
         WRITE(LUPRI,*) 'LVIRB:'
         WRITE(LUPRI,'(8I10)') (LVIRB(I),I=1,NSYM)
      ENDIF

C     Check input. Integral/amplitude symmetry ISYINT is set here.
C     ------------------------------------------------------------

      CALL CC_CAAIN1(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &               UVEC,ISYMU,VVEC,ISYMV,FAC2,NUMUV,
     &               SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &               X2AM,WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2,
     &               IERR,LOCDBG,ISYINT)
      IF (IERR .NE. 0) CALL QUIT('Inconsistent input in '//SECNAM)

C     Issue warnings for small scaling factors, but don't stop.
C     ---------------------------------------------------------

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

C     Initializations.
C     ----------------

      MEMUSE = 0
      X2NRM  = ZERO
      X2CNM  = ZERO

C     Calculate integrals from Cholesky vectors.
C     ------------------------------------------

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

         NPASS = 0

         DTIME = SECOND()
         CALL CC_CIA(X2AM,WORK,LWORK,SCDLOC,KLAST,
     &               NPASS)
         DTIME = SECOND() - DTIME
         TIMI  = TIMI     + DTIME

         MEMUSE = MAX(MEMUSE,KLAST)

      ENDDO

C     Scale integrals or initialize (if no integrals calculated).
C     -----------------------------------------------------------

      IF (NUMP12 .GT. 0) THEN
         IF (FAC1 .NE. ONE) CALL DSCAL(NX2SQ,FAC1,X2AM,1)
      ELSE
         CALL DZERO(X2AM,NX2SQ)
      ENDIF

C     Add FAC2*U*V terms (if any).
C     ----------------------------

      KOFFU = 1
      KOFFV = 1
      DO IUV = 1,NUMUV

         ISYMAI = ISYMU(IUV)
         ISYMBJ = ISYMV(IUV)

         DTIME = SECOND()
         CALL CC_CYILHX(UVEC(KOFFU),X2AM,VVEC(KOFFV),
     &                  WORK,LWORK,KLAST,ISYMAI,ISYMBJ,
     &                  FAC2)
         DTIME = SECOND() - DTIME
         TIMUV = TIMUV    + DTIME

         MEMUSE = MAX(MEMUSE,KLAST)

         KOFFU = KOFFU + NT1AM(ISYMAI)
         KOFFV = KOFFV + NT1AM(ISYMBJ)

      ENDDO

C     Scale diagonal (if needed).
C     ---------------------------

      IF ((SCDG.NE.ONE) .AND. (ISYINT.EQ.1)) THEN
         XTST   = ONE*NX2SQ
         DIFF   = DABS(XT2SQ(1) - XTST)
         FULLIN = DIFF .LE. TINY
         CALL CC_CIADSCL(X2AM,ISYINT,SCDG,FULLIN)
      ENDIF

C     Divide by orbital energies to obtain amplitudes,
C     if requested through IOPTDN.
C     ------------------------------------------------

      IF (IOPTDN .EQ. 1) THEN
         CALL CC_DNOM(FOCKD,X2AM,FREQ,ISYINT)
      ENDIF

C     Calculate contribution to (packed) doubles norm.
C     ------------------------------------------------

      CALL CC_CYINRM(X2AM,ISYINT,X2NRM)

C     Set up 2 Coulomb minus exchange (almost) in place, if
C     requested through IOPTCE, and compute contribution to
C     packed 2CME norm.
C     -----------------------------------------------------

      IF (IOPTCE .EQ. 1) THEN
         CALL CC_CYITCME(X2AM,WORK,LWORK,ISYINT,KLAST)
         CALL CC_CYINRM(X2AM,ISYINT,X2CNM)
         MEMUSE = MAX(MEMUSE,KLAST)
      ENDIF

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

      IF ((IPRINT.GE.INFO) .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         CALL HEADER('Output from '//SECNAM,-1)
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Memory available : ',LWORK,
     &   'Max. memory usage: ',MEMUSE
         WRITE(LUPRI,'(5X,A,I2)')
     &   'Amplitude symmetry:',ISYINT
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Packed amplitude norm: ',DSQRT(X2NRM)
         IF (IOPTCE .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &      'Packed 2CME am.  norm: ',DSQRT(X2CNM)
         ENDIF
         IF (IOPTDN .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &      'Denominator frequency: ',FREQ
         ENDIF
         WRITE(LUPRI,'(/,5X,A,I3,A,F10.2,A)')
     &   'Time used for',NUMP12,' integral calculations: ',TIMI,
     &   ' seconds.'
         WRITE(LUPRI,'(5X,A,I3,A,F10.2,A)')
     &   'Time used for',NUMUV,' U*V      calculations: ',TIMUV,
     &   ' seconds.'
         WRITE(LUPRI,'(5X,A)')
     &   '-----------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Time used in total                    : ',TIMT,' seconds.'
      ENDIF

      RETURN
      END
C  /* Deck cc_caain1 */
      SUBROUTINE CC_CAAIN1(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHP12,FAC1,NUMP12,
     &                     UVEC,ISYMU,VVEC,ISYMV,FAC2,NUMUV,
     &                     SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &                     X2AM,WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2,
     &                     IERR,PRNT,ISYINT)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Check input to CC_CAA.
C
#include <implicit.h>
      DIMENSION UVEC(*), VVEC(*), FOCKD(*)
      DIMENSION SCD(NUMP12), X2AM(*), WORK(LWORK)
      INTEGER JTYP1(NUMP12), ISYMP1(NUMP12)
      INTEGER JTYP2(NUMP12), ISYMP2(NUMP12)
      INTEGER NCHP12(8,NUMP12)
      INTEGER ISYMU(NUMUV), ISYMV(NUMUV)
      LOGICAL DELP1(NUMP12), DELP2(NUMP12)
      LOGICAL PRNT
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CAAIN1')

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

      ENDIF

C     Print section.
C     --------------

  100 IF (PRNT .OR. (IERR.NE.0)) THEN

         WRITE(LUPRI,'(//,1X,A,A,A)')
     &   'Input to CC_CAA checked by ',SECNAM,':'
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
         WRITE(LUPRI,*) 'LWORK : ',LWORK
         WRITE(LUPRI,*) 'DELP1 : ',(DELP1(I),I=1,NUMP12)
         WRITE(LUPRI,*) 'DELP2 : ',(DELP2(I),I=1,NUMP12)

         CALL FLSHFO(LUPRI)

      ENDIF

      RETURN
      END
