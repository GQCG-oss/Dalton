C  /* Deck cc_choatr */
      SUBROUTINE CC_CHOATR(SIGMA1,X1AM,FREQ,WORK,LWORK,ISYMX,NUMX,
     &                      ISIDE)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate NUMX linear transformations with the effective CC2
C              Jacobian from left (ISIDE=-1) or right (ISIDE=1).
C
C        ISIDE = -1: SIGMA1 = X1AM * A(FREQ)  [NOTE THE SIGN OF FREQ]
C
C        ISIDE =  1: SIGMA1 = A(FREQ) * X1AM
C
C     The (global) E intermediates must be available on disk (files opened
C     through subroutine CHO_IMOP).
C
#include <implicit.h>
      DIMENSION SIGMA1(*), X1AM(*), FREQ(NUMX), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOATR')

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,'(/,5X,A,A,I10)')
     &   SECNAM,': ISIDE = ',ISIDE
         WRITE(LUPRI,'(5X,A,A,I10)')
     &   SECNAM,': ISYMX = ',ISYMX
         WRITE(LUPRI,'(5X,A,A,I10)')
     &   SECNAM,': NUMX  = ',NUMX
         WRITE(LUPRI,'(5X,A,A)')
     &   SECNAM,': Frequencies:'
         WRITE(LUPRI,'(1P,D20.12)') (FREQ(I), I = 1,NUMX)
         WRITE(LUPRI,*)
      ENDIF

C     Calculate contribution from global E intermediates
C     (initialization of result vectors in SIGMA1).
C     ---------------------------------------------------

      CALL CC_CHOATR0(SIGMA1,X1AM,WORK,LWORK,ISYMX,NUMX,ISIDE)

      IF (LOCDBG) THEN
         DO I = 1,NUMX
            KOFF = NT1AM(ISYMX)*(I - 1) + 1
            SNRM = DSQRT(DDOT(NT1AM(ISYMX),SIGMA1(KOFF),1,
     &                                     SIGMA1(KOFF),1))
            WRITE(LUPRI,'(5X,A,I3,A,1P,D22.15)')
     &      'Norm of transformed vector',I,' after CC_CHOATR0: ',SNRM
         ENDDO
      ENDIF

C     Read the canonical orbital energies and calculate Lambda matrices.
C     ------------------------------------------------------------------

      KFOCKD = 1
      KLAMDP = KFOCKD + NORBTS
      KLAMDH = KLAMDP + NGLMDT(1)
      KT1AM  = KLAMDH + NGLMDT(1)
      KEND1  = KT1AM  + NT1AM(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: orb. en.'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &               .FALSE.,.TRUE.,.FALSE.)
      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     &            LWRK1)

C     Read unperturbed F(ia) into T1AM (no longer needed).
C     ----------------------------------------------------

      KFIA = KT1AM
      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(WORK(KFIA),NT1AM(1),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

C     Calculate C intermediates for ISIDE = -1, or
C     calculate the I2J contribution for ISIDE = 1.
C     ---------------------------------------------

      IF (ISIDE .EQ. -1) THEN
         CALL CC_CHOCIM(WORK(KFOCKD),X1AM,WORK(KEND1),LWRK1,ISYMX,NUMX)
      ELSE
         CALL CC_CHORTRI2J(WORK(KFOCKD),WORK(KLAMDP),WORK(KLAMDH),
     &                     SIGMA1,X1AM,WORK(KEND1),LWRK1,ISYMX,NUMX)
         IF (LOCDBG) THEN
            DO I = 1,NUMX
               KOFF = NT1AM(ISYMX)*(I - 1) + 1
               SNRM = DSQRT(DDOT(NT1AM(ISYMX),SIGMA1(KOFF),1,
     &                                        SIGMA1(KOFF),1))
               WRITE(LUPRI,'(5X,A,I3,A,1P,D22.15)')
     &         'Norm of transformed vector',I,' after CC_CHORTRI2J: ',
     &         SNRM
            ENDDO
         ENDIF
      ENDIF

C     Loop over trial vectors for calculating doubles-dependent terms.
C     ================================================================

      DO IX = 1,NUMX

C        Offset to vectors in SIGMA1 and X1AM.
C        -------------------------------------

         KOFIX = NT1AM(ISYMX)*(IX - 1) + 1

C        Calculate perturbed Cholesky vectors, store on disk.
C        ----------------------------------------------------

         CALL CC_CHOTG(WORK(KLAMDP),WORK(KLAMDH),X1AM(KOFIX),
     &                 WORK(KEND1),LWRK1,1,ISYMX,ISIDE)

C        Calculate Y intermediates. For ISIDE = 1, also
C        calculate the I1 contribution. Perturbed Cholesky
C        vector files are deleted at the end of CC_CYI.
C        -------------------------------------------------

c        CALL CC_CYI(WORK(KFOCKD),SIGMA1(KOFIX),X1AM(KOFIX),WORK(KFIA),
c    &               WORK(KEND1),LWRK1,ISIDE,ISYMX,1,FREQ(IX),X2NRM,
c    &               X2CNM,.TRUE.)

         IF (ISIDE .EQ. -1) THEN
            CALL CC_CYILTR(WORK(KFOCKD),X1AM(KOFIX),WORK(KFIA),
     &                     WORK(KEND1),LWRK1,ISYMX,FREQ(IX),X2NRM,
     &                     X2CNM,.FALSE.,.TRUE.)
         ELSE IF (ISIDE .EQ. 1) THEN
            CALL CC_CYIRTR(WORK(KFOCKD),SIGMA1(KOFIX),WORK(KFIA),
     &                     WORK(KEND1),LWRK1,ISYMX,FREQ(IX),X2NRM,
     &                     X2CNM,.FALSE.,.TRUE.)
         ELSE
            CALL QUIT('ISIDE error in '//SECNAM)
         ENDIF

         IF (LOCDBG) THEN
            SNRM = DSQRT(DDOT(NT1AM(ISYMX),SIGMA1(KOFIX),1,
     &                                     SIGMA1(KOFIX),1))
            WRITE(LUPRI,'(5X,A,I3,A,1P,D22.15)')
     &      'Norm of transformed vector',IX,' after CC_CYI: ',SNRM
         ENDIF

C        For ISIDE=-1: Calculate GIJ and H terms.
C        For ISIDE=+1: Calculate G   and H terms.
C        Delete Y intermediate files after use.
C        -----------------------------------------

         CALL CC_CHOATR1(WORK(KLAMDP),WORK(KLAMDH),SIGMA1(KOFIX),
     &                   X1AM(KOFIX),WORK(KEND1),LWRK1,ISYMX,ISIDE,
     &                   .TRUE.,IX)

         IF (LOCDBG) THEN
            SNRM = DSQRT(DDOT(NT1AM(ISYMX),SIGMA1(KOFIX),1,
     &                                     SIGMA1(KOFIX),1))
            WRITE(LUPRI,'(5X,A,I3,A,1P,D22.15)')
     &      'Norm of transformed vector',IX,' after CC_CHOATR1: ',SNRM
         ENDIF
 
      ENDDO

      RETURN
      END
C  /* Deck cc_cyirtr */
      SUBROUTINE CC_CYIRTR(FOCKD,SIGMA1,FIA,WORK,LWORK,ISYMX,FREQ,X2NRM,
     &                     X2CNM,DELP1,DELP2)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates for right-hand Jacobian
C              transformation.
C
#include <implicit.h>
      DIMENSION FOCKD(*), SIGMA1(*), FIA(*), WORK(LWORK)
      LOGICAL DELP1, DELP2
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>

C     Set Cholesky vector type for amplitude assembly.
C     ------------------------------------------------

      FAC1   = 1.0D0
      FAC2   = 0.0D0
      SCD    = 1.0D0
      SCDG   = 1.0D0

      NUMP12 = 1
      JTYP1  = 3
      JTYP2  = 9
      ISYMP1 = 1
      ISYMP2 = ISYMX

      IOPTDN = 1
      IOPTCE = 1

C     Set info for Y intermediates.
C     -----------------------------

      NUMZ  = 1
      JTYPZ = 1
      ISYMZ = 1
      JTYPY = 1

C     Calculate Y intermediates.
C     --------------------------

      CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NUMCHO,FAC1,NUMP12,
     &            DUMMY,IDUMMY,DUMMY,IDUMMY,FAC2,0,
     &            SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &            JTYPZ,ISYMZ,NUMCHO,NUMZ,JTYPY,
     &            FIA,1,1,SIGMA1,
     &            WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2)

      RETURN
      END
C  /* Deck cc_cyiltr */
      SUBROUTINE CC_CYILTR(FOCKD,X1AM,FIA,WORK,LWORK,ISYMX,FREQ,X2NRM,
     &                     X2CNM,DELP1,DELP2)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates for left-hand Jacobian
C              transformation.
C
#include <implicit.h>
      DIMENSION FOCKD(*), X1AM(*), FIA(*), WORK(LWORK)
      LOGICAL DELP1, DELP2
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>

C     Set Cholesky vector type for amplitude assembly.
C     ------------------------------------------------

      FAC1   = 1.0D0
      FAC2   = 1.0D0
      SCD    = 1.0D0
      SCDG   = 1.0D0

      NUMP12 = 1
      JTYP1  = 1
      JTYP2  = 8
      ISYMP1 = 1
      ISYMP2 = ISYMX

      IOPTDN = 1
      IOPTCE = 1

C     Set info for Y intermediates.
C     -----------------------------

      NUMZ  = 1
      JTYPZ = 3
      ISYMZ = 1
      JTYPY = -1

C     Calculate Y intermediates.
C     --------------------------

      CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NUMCHO,FAC1,NUMP12,
     &            X1AM,ISYMX,FIA,1,FAC2,1,
     &            SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &            JTYPZ,ISYMZ,NUMCHO,NUMZ,JTYPY,
     &            DUMMY,IDUMMY,0,DUMMY,
     &            WORK,LWORK,X2NRM,X2CNM,DELP1,DELP2)

      RETURN
      END
C  /* Deck cc_choatr1 */
      SUBROUTINE CC_CHOATR1(XLAMDP,XLAMDH,SIGMA1,X1AM,WORK,LWORK,ISYMX,
     &                      ISIDE,DELY,IX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the GIJ and H terms for ISIDE=-1, or
C              the G and H terms for ISIDE=1.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), SIGMA1(*), X1AM(*), WORK(LWORK)
      LOGICAL DELY
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOATR1')

C     Call appropriate routine.
C     -------------------------

      IF (ISIDE .EQ. -1) THEN

         KCIM  = 1
         KEND1 = KCIM  + NT1AM(ISYMX)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL CHO_IMOP(-1,3,LUCIM,ISYMX)
         CALL CHO_MOREAD(WORK(KCIM),NT1AM(ISYMX),1,IX,LUCIM)
         CALL CHO_IMOP(1,3,LUCIM,ISYMX)

         CALL CC_CHOLTRGIJH(XLAMDP,XLAMDH,SIGMA1,X1AM,WORK(KCIM),
     &                      WORK(KEND1),LWRK1,ISYMX)

      ELSE IF (ISIDE .EQ. 1) THEN

         CALL CC_CHORTRGH(XLAMDP,XLAMDH,SIGMA1,WORK,LWORK,ISYMX)

      ELSE

         WRITE(LUPRI,'(//,5X,A,A)')
     &   'ISIDE must be -1 or +1 in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'ISIDE = ',ISIDE
         CALL QUIT('Error in '//SECNAM)

      ENDIF

C     Delete Y intermediate files if requested.
C     -----------------------------------------

      IF (DELY) THEN
         DO ISYCHO = 1,NSYM
            CALL CC_CYIOP(-1,ISYCHO,ISIDE)
            CALL CC_CYIOP(0,ISYCHO,ISIDE)
         ENDDO
      ENDIF

      RETURN
      END
C  /* Deck cc_chortrgh */
      SUBROUTINE CC_CHORTRGH(XLAMDP,XLAMDH,SIGMA1,WORK,LWORK,ISYMY)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the G and H terms for right-hand Jacobian
C              transformations.
C
C     Formula:
C
C     SIGMA1(ai) = SIGMA1(ai)
C                + sum(g) LamP(g,a) sum(Jd) L(gd,J) sum(b) LamH(d,b) * Y(bi,J)
C                - sum(Jj) Y(aj,J) * L(ji,J)
C
C     where g,d refer to AO basis.
C     
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), SIGMA1(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER IOFFL(8), IOFFZ(8), IOFFY(8)

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHORTRGH')

      PARAMETER (IOPTR = 2)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

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

C     Allocation.
C     -----------

      KSIGMA = KEND0
      KEND1  = KSIGMA + NT1AO(ISYMY)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize AO G term.
C     ---------------------

      CALL DZERO(WORK(KSIGMA),NT1AO(ISYMY))

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000

         ISYMBI = MULD2H(ISYCHO,ISYMY)
         ISYMAJ = ISYMBI
         ISYMDI = ISYMBI

C        Allocation for one AO Cholesky vector.
C        --------------------------------------

         LREAD = MEMRD(1,ISYCHO,IOPTR)

         KCHOL = KEND1
         KREAD = KCHOL + N2BST(ISYCHO)
         KEND2 = KREAD + LREAD
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        SCRL: L(ji,#J),Y(b#J,i),L(g,d#J).
C        SCRZ: Y(bi,#J),Z(d#J,i).
C        ---------------------------------

         LENL = MAX(NMATIJ(ISYCHO),NT1AM(ISYMBI),N2BST(ISYCHO))
         LENZ = MAX(NT1AM(ISYMBI),NT1AO(ISYMDI))

         MINMEM = LENL + LENZ
         IF (MINMEM .GT. 0) THEN
            NVEC = MIN(LWRK2/MINMEM,NUMCHO(ISYCHO))
         ELSE IF (MINMEM .EQ. 0) THEN
            GO TO 1000
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &      'MINMEM is negative in ',SECNAM,': ',MINMEM
            CALL QUIT('Error in '//SECNAM)
         ENDIF

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I2,A)')
     &      '(Cholesky symmetry:',ISYCHO,')'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required for batch: ',MINMEM,
     &      'Memory available for batch       : ',LWRK2,
     &      'Total memory available           : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Open Y intermediate file.
C        -------------------------

         CALL CC_CYIOP(-1,ISYCHO,1)

C        Open occupied Cholesky vector file.
C        -----------------------------------

         CALL CHO_MOP(-1,4,ISYCHO,LUCHJI,1,1)

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF
            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Set up local index arrays.
C           --------------------------

            ICOUNL = 0
            ICOUNZ = 0
            ICOUNY = 0
            DO ISYMD = 1,NSYM
               ISYMI = MULD2H(ISYMD,ISYMDI)
               ISYMG = MULD2H(ISYMD,ISYCHO)
               ISYMB = ISYMD
               IOFFL(ISYMD) = ICOUNL
               IOFFZ(ISYMD) = ICOUNZ
               IOFFY(ISYMB) = ICOUNY
               LENDJ  = NBAS(ISYMD)*NUMV
               LENBJ  = NVIR(ISYMB)*NUMV
               ICOUNL = ICOUNL + NBAS(ISYMG)*LENDJ
               ICOUNZ = ICOUNZ + LENDJ*NRHF(ISYMI)
               ICOUNY = ICOUNY + LENBJ*NRHF(ISYMI)
            ENDDO

C           Complete allocation.
C           --------------------

            KSCRL = KEND2
            KSCRZ = KSCRL + LENL*NUMV

C           Read Y(bi,#J) in place of Z.
C           ----------------------------

            CALL CC_CYIRDF(WORK(KSCRZ),NUMV,JVEC1,ISYCHO,ISYMY,1)

C           Read L(ji,#J) in place of L(g,d#J).
C           -----------------------------------

            CALL CHO_MOREAD(WORK(KSCRL),NMATIJ(ISYCHO),NUMV,JVEC1,
     &                      LUCHJI)

C           Calculate H term.
C           -----------------

            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMI = MULD2H(ISYMJ,ISYCHO)
                  ISYMA = MULD2H(ISYMJ,ISYMAJ)

                  KOFF1 = KSCRZ + NT1AM(ISYMAJ)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMJ)
                  KOFF2 = KSCRL + NMATIJ(ISYCHO)*(IVEC - 1)
     &                  + IMATIJ(ISYMJ,ISYMI)
                  KOFF3 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTJ = MAX(NRHF(ISYMJ),1)

                  CALL DGEMM('N','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTJ,
     &                       ONE,SIGMA1(KOFF3),NTOTA)

               ENDDO
            ENDDO

C           Resort Y intermediate as Y(b#J,i), store in place of L(g,d#J).
C           --------------------------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYMBI)
                  LENBJ = NVIR(ISYMB)*NUMV

                  DO I = 1,NRHF(ISYMI)

                     KOFF1 = KSCRZ + NT1AM(ISYMBI)*(IVEC - 1)
     &                     + IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1)
                     KOFF2 = KSCRL + IOFFY(ISYMB) + LENBJ*(I - 1)
     &                     + NVIR(ISYMB)*(IVEC - 1)

                     CALL DCOPY(NVIR(ISYMB),WORK(KOFF1),1,WORK(KOFF2),1)

                  ENDDO

               ENDDO
            ENDDO

C           Backtransform the Y intermediates to get Z(d#J,i).
C           --------------------------------------------------

            DO ISYMB = 1,NSYM

               ISYMI = MULD2H(ISYMB,ISYMBI)
               ISYMD = ISYMB

               KOFF1 = IGLMVI(ISYMD,ISYMB) + 1
               KOFF2 = KSCRL + IOFFY(ISYMB)
               KOFF3 = KSCRZ + IOFFZ(ISYMD)

               NTOTD = MAX(NBAS(ISYMD),1)
               NTOTB = MAX(NVIR(ISYMB),1)

               CALL DGEMM('N','N',
     &                    NBAS(ISYMD),NUMV*NRHF(ISYMI),NVIR(ISYMB),
     &                    ONE,XLAMDH(KOFF1),NTOTD,WORK(KOFF2),NTOTB,
     &                    ZERO,WORK(KOFF3),NTOTD)

            ENDDO

C           Read AO Cholesky vectors and store as L(g,d#J).
C           -----------------------------------------------

            DO IVEC = 1,NUMV

               JVEC = JVEC1 + IVEC - 1
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KREAD),LREAD)

               DO ISYMD = 1,NSYM

                  ISYMG = MULD2H(ISYMD,ISYCHO)

                  LENGD = NBAS(ISYMG)*NBAS(ISYMD)

                  KOFF1 = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFF2 = KSCRL + IOFFL(ISYMD) + LENGD*(IVEC - 1)

                  CALL DCOPY(LENGD,WORK(KOFF1),1,WORK(KOFF2),1)

               ENDDO

            ENDDO

C           Calculate G term in AO basis.
C           -----------------------------

            DO ISYMD = 1,NSYM

               ISYMG = MULD2H(ISYMD,ISYCHO)
               ISYMI = MULD2H(ISYMD,ISYMDI)

               KOFF1 = KSCRL  + IOFFL(ISYMD)
               KOFF2 = KSCRZ  + IOFFZ(ISYMD)
               KOFF3 = KSIGMA + IT1AO(ISYMG,ISYMI)

               LENDJ = NBAS(ISYMD)*NUMV

               NTOTG = MAX(NBAS(ISYMG),1)
               NTODJ = MAX(LENDJ,1)

               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),LENDJ,
     &                    ONE,WORK(KOFF1),NTOTG,WORK(KOFF2),NTODJ,
     &                    ONE,WORK(KOFF3),NTOTG)

            ENDDO

         ENDDO

C        Close occupied Cholesky vector file.
C        ------------------------------------

         CALL CHO_MOP(1,4,ISYCHO,LUCHJI,1,1)

C        Close Y intermediate file.
C        --------------------------

         CALL CC_CYIOP(1,ISYCHO,1)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

C     Transform the G term to MO basis.
C     ---------------------------------

      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYMY)
         ISYMG = ISYMA

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFF1 = IGLMVI(ISYMG,ISYMA) + 1
         KOFF2 = KSIGMA + IT1AO(ISYMG,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFF1),NTOTG,WORK(KOFF2),NTOTG,
     &              ONE,SIGMA1(KOFF3),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choltrgijh */
      SUBROUTINE CC_CHOLTRGIJH(XLAMDP,XLAMDH,SIGMA1,X1AM,CIM,WORK,LWORK,
     &                         ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the GIJ and H terms for left-hand Jacobian
C              transformations.
C
C     Formula:
C
C     SIGMA1(ai) = SIGMA1(ai)
C                - sum(Jj) Y(aj,J) * L(ij,J)
C                + sum(g) LamH(g,a)
C                 * sum(Jd) L(gd,J)
C                   * [ sum(b) LamP(d,b) * {Y(bi,J) - sum(j) X(bj)*L(ij,J)}
C                     - sum(j) LamP(d,j) * sum(b) L(ib,J)*C(bj)
C                     + 2 LamP(d,i) * sum(bj) {L(jb,J)*C(bj) + L(bj,J)*X(bj)} ]
C
C     where g,d refer to AO basis, and L(gd,J) = L(dg,J) is used.
C     
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), SIGMA1(*), X1AM(*), CIM(*)
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER IOFFL(8), IOFFZ(8), IOFFY(8)

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CC_CHOLTRGIJH')

      PARAMETER (IOPTR = 2)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

c      snrm = dsqrt(ddot(nt1am(isymx),sigma1,1,sigma1,1))
c      xnrm = dsqrt(ddot(nt1am(isymx),x1am,1,x1am,1))
c      cnrm = dsqrt(ddot(nt1am(isymx),cim,1,cim,1))
c      write(lupri,'(A,A,1P,D22.15)')
c     & SECNAM,': norm of sigma1: ',snrm
c      write(lupri,'(A,A,1P,D22.15)')
c     & SECNAM,': norm of x1am  : ',xnrm
c      write(lupri,'(A,A,1P,D22.15)')
c     & SECNAM,': norm of cim   : ',cnrm

C     Set symmetries.
C     ---------------

      ISYCIM = ISYMX
      ISYMY  = ISYMX

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

C     Allocation.
C     -----------

      KSIGMA = KEND0
      KEND1  = KSIGMA + NT1AO(ISYMY)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize AO GIJ term.
C     -----------------------

      CALL DZERO(WORK(KSIGMA),NT1AO(ISYMY))

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000

         ISYMBI = MULD2H(ISYCHO,ISYMY)
         ISYMAJ = ISYMBI
         ISYMDI = ISYMBI
         ISYMKI = ISYMBI

C        Allocation for one AO Cholesky vector.
C        --------------------------------------

         LREAD = MEMRD(1,ISYCHO,IOPTR)
         MAXKI = 0
         DO ISYMI = 1,NSYM
            ISYMK = MULD2H(ISYMI,ISYMKI)
            NTEST = NRHF(ISYMK)*NRHF(ISYMI)
            MAXKI = MAX(MAXKI,NTEST)
         ENDDO
         LSCR  = MAX(LREAD,MAXKI)

         KCHOL = KEND1
         KREAD = KCHOL + N2BST(ISYCHO)
         KEND2 = KREAD + LSCR
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        SCRL: L(ij,#J),Y(b#J,i),L(ck,#J),L(ic,#J),L(g,d#J).
C        SCRZ: Y(bi,#J),Z(d#J,i).
C        ---------------------------------------------------

         LENL = MAX(NMATIJ(ISYCHO),NT1AM(ISYMBI),NT1AM(ISYCHO),
     &              N2BST(ISYCHO))
         LENZ = MAX(NT1AM(ISYMBI),NT1AO(ISYMDI))

         MINMEM = LENL + LENZ
         IF (ISYCHO .EQ. ISYMX) MINMEM = MINMEM + 1
         IF (MINMEM .GT. 0) THEN
            NVEC = MIN(LWRK2/MINMEM,NUMCHO(ISYCHO))
         ELSE IF (MINMEM .EQ. 0) THEN
            GO TO 1000
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,A,I10,/)')
     &      'MINMEM is negative in ',SECNAM,': ',MINMEM
            CALL QUIT('Error in '//SECNAM)
         ENDIF

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I2,A)')
     &      '(Cholesky symmetry:',ISYCHO,')'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required for batch: ',MINMEM,
     &      'Memory available for batch       : ',LWRK2,
     &      'Total memory available           : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Open Y intermediate file.
C        -------------------------

         CALL CC_CYIOP(-1,ISYCHO,-1)

C        Open MO Cholesky vector files.
C        ------------------------------

         CALL CHO_MOP(-1,1,ISYCHO,LUCHIA,1,1)
         IF (ISYCHO .EQ. ISYMX) CALL CHO_MOP(-1,3,ISYCHO,LUCHAI,1,1)
         CALL CHO_MOP(-1,4,ISYCHO,LUCHIJ,1,1)

C        Start batch loop.
C        -----------------

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF
            JVEC1 = NVEC*(IBATCH - 1) + 1

C           Set up local index arrays.
C           --------------------------

            ICOUNL = 0
            ICOUNZ = 0
            ICOUNY = 0
            DO ISYMD = 1,NSYM
               ISYMI = MULD2H(ISYMD,ISYMDI)
               ISYMG = MULD2H(ISYMD,ISYCHO)
               ISYMB = ISYMD
               IOFFL(ISYMD) = ICOUNL
               IOFFZ(ISYMD) = ICOUNZ
               IOFFY(ISYMB) = ICOUNY
               LENDJ  = NBAS(ISYMD)*NUMV
               LENBJ  = NVIR(ISYMB)*NUMV
               ICOUNL = ICOUNL + NBAS(ISYMG)*LENDJ
               ICOUNZ = ICOUNZ + LENDJ*NRHF(ISYMI)
               ICOUNY = ICOUNY + LENBJ*NRHF(ISYMI)
            ENDDO

C           Complete allocation.
C           --------------------

            KSCRL = KEND2
            KSCRZ = KSCRL + LENL*NUMV
            KUVEC = KSCRZ + LENZ*NUMV

C           Read Y(bi,#J) in place of Z.
C           ----------------------------

            CALL CC_CYIRDF(WORK(KSCRZ),NUMV,JVEC1,ISYCHO,ISYMY,-1)

c      ynrm1 = dsqrt(ddot(nt1am(isymbi)*numv,work(kscrz),1,
c     &                                      work(kscrz),1))
c      write(lupri,*) '...ynrm at read            : ',ynrm1

C           Read L(ij,#J) in place of L(g,d#J).
C           -----------------------------------

            CALL CHO_MOREAD(WORK(KSCRL),NMATIJ(ISYCHO),NUMV,JVEC1,
     &                      LUCHIJ)

C           Calculate H term.
C           -----------------

            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMI = MULD2H(ISYMJ,ISYCHO)
                  ISYMA = MULD2H(ISYMJ,ISYMAJ)

                  KOFF1 = KSCRZ + NT1AM(ISYMAJ)*(IVEC - 1)
     &                  + IT1AM(ISYMA,ISYMJ)
                  KOFF2 = KSCRL + NMATIJ(ISYCHO)*(IVEC - 1)
     &                  + IMATIJ(ISYMI,ISYMJ)
                  KOFF3 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTI,
     &                       ONE,SIGMA1(KOFF3),NTOTA)

               ENDDO
            ENDDO

C           Subtract sum(j) X(bj) * L(ij,#J) from Y(bi,#J).
C           -----------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMJ = 1,NSYM

                  ISYMI = MULD2H(ISYMJ,ISYCHO)
                  ISYMB = MULD2H(ISYMJ,ISYMX)

                  KOFF1 = IT1AM(ISYMB,ISYMJ) + 1
                  KOFF2 = KSCRL + NMATIJ(ISYCHO)*(IVEC - 1)
     &                  + IMATIJ(ISYMI,ISYMJ)
                  KOFF3 = KSCRZ + NT1AM(ISYMBI)*(IVEC - 1)
     &                  + IT1AM(ISYMB,ISYMI)

                  NTOTB = MAX(NVIR(ISYMB),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMB),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,X1AM(KOFF1),NTOTB,WORK(KOFF2),NTOTI,
     &                       ONE,WORK(KOFF3),NTOTB)

               ENDDO
            ENDDO

c      ynrm2 = dsqrt(ddot(nt1am(isymbi)*numv,work(kscrz),1,
c     &                                      work(kscrz),1))
c      write(lupri,*) '...ynrm after X-subtraction: ',ynrm2

C           Resort Y intermediate as Y(b#J,i), store in place of L(g,d#J).
C           --------------------------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYMBI)
                  LENBJ = NVIR(ISYMB)*NUMV

                  DO I = 1,NRHF(ISYMI)

                     KOFF1 = KSCRZ + NT1AM(ISYMBI)*(IVEC - 1)
     &                     + IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1)
                     KOFF2 = KSCRL + IOFFY(ISYMB) + LENBJ*(I - 1)
     &                     + NVIR(ISYMB)*(IVEC - 1)

                     CALL DCOPY(NVIR(ISYMB),WORK(KOFF1),1,WORK(KOFF2),1)

                  ENDDO

               ENDDO
            ENDDO

c      ynrm3 = dsqrt(ddot(nt1am(isymbi)*numv,work(kscrl),1,
c     &                                      work(kscrl),1))
c      write(lupri,*) '...ynrm after resort       : ',ynrm3

C           Backtransform the Y intermediates to get Z(d#J,i).
C           --------------------------------------------------

            DO ISYMB = 1,NSYM

               ISYMI = MULD2H(ISYMB,ISYMBI)
               ISYMD = ISYMB

               KOFF1 = IGLMVI(ISYMD,ISYMB) + 1
               KOFF2 = KSCRL + IOFFY(ISYMB)
               KOFF3 = KSCRZ + IOFFZ(ISYMD)

               NTOTD = MAX(NBAS(ISYMD),1)
               NTOTB = MAX(NVIR(ISYMB),1)

               CALL DGEMM('N','N',
     &                    NBAS(ISYMD),NUMV*NRHF(ISYMI),NVIR(ISYMB),
     &                    ONE,XLAMDP(KOFF1),NTOTD,WORK(KOFF2),NTOTB,
     &                    ZERO,WORK(KOFF3),NTOTD)

            ENDDO

            IF (ISYCHO .EQ. ISYMX) THEN

C              Read L(ck,#J) in place of L(g,d#J).
C              -----------------------------------

               CALL CHO_MOREAD(WORK(KSCRL),NT1AM(ISYCHO),NUMV,JVEC1,
     &                         LUCHAI)

C              Calculate U(#J) = 2 * sum(ck) L(ck,#J) * X(ck).
C              -----------------------------------------------

               NTOCK = MAX(NT1AM(ISYCHO),1)
               CALL DGEMV('T',NT1AM(ISYCHO),NUMV,
     &                    TWO,WORK(KSCRL),NTOCK,X1AM,1,
     &                    ZERO,WORK(KUVEC),1)

            ENDIF

C           Read L(ic,#J) in place of L(g,d#J).
C           -----------------------------------

            CALL CHO_MOREAD(WORK(KSCRL),NT1AM(ISYCHO),NUMV,JVEC1,LUCHIA)

            IF (ISYCHO .EQ. ISYCIM) THEN

C              Calculate U(#J) = U(#J) + 2 * sum(ck) L(kc,#J) * CIM(ck).
C              ---------------------------------------------------------

               NTOCK = MAX(NT1AM(ISYCHO),1)
               CALL DGEMV('T',NT1AM(ISYCHO),NUMV,
     &                    TWO,WORK(KSCRL),NTOCK,CIM,1,
     &                    ONE,WORK(KUVEC),1)

C              Set up Z contribution:
C              Z(d#J,i) = Z(d#J,i) + LamP(d,i) * U(#J)
C              ---------------------------------------

               DO ISYMD = 1,NSYM
                  ISYMI = ISYMD
                  DO I = 1,NRHF(ISYMI)
                     KOFF1 = IGLMRH(ISYMD,ISYMI)
     &                     + NBAS(ISYMD)*(I - 1) + 1
                     DO IVEC = 1,NUMV
                        FAC = WORK(KUVEC+IVEC-1)
                        KOFF2 = KSCRZ + IOFFZ(ISYMD)
     &                        + NBAS(ISYMD)*NUMV*(I - 1)
     &                        + NBAS(ISYMD)*(IVEC - 1)
                        CALL DAXPY(NBAS(ISYMD),FAC,XLAMDP(KOFF1),1,
     &                                             WORK(KOFF2),1)
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF

C           Transform: -sum(k) LamP(d,k) {sum(c) C(ck) * L(ic,#J)}
C           and add result into Z(d#J,i). Store intermediate result
C           in the scratch space for reading AO Cholesky vectors.
C           -------------------------------------------------------

            DO IVEC = 1,NUMV
               DO ISYMI = 1,NSYM

                  ISYMC = MULD2H(ISYMI,ISYCHO)
                  ISYMK = MULD2H(ISYMC,ISYCIM)

                  KOFF1 = IT1AM(ISYMC,ISYMK) + 1
                  KOFF2 = KSCRL + NT1AM(ISYCHO)*(IVEC - 1)
     &                  + IT1AM(ISYMC,ISYMI)

                  NTOTC = MAX(NVIR(ISYMC),1)
                  NTOTK = MAX(NRHF(ISYMK),1)

                  CALL DGEMM('T','N',
     &                       NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     &                       ONE,CIM(KOFF1),NTOTC,WORK(KOFF2),NTOTC,
     &                       ZERO,WORK(KREAD),NTOTK)

                  ISYMD = ISYMK
                  KOFF1 = IGLMRH(ISYMD,ISYMK) + 1
                  NTOTD = MAX(NBAS(ISYMD),1)

                  DO I = 1,NRHF(ISYMI)

                     KOFF2 = KREAD + NRHF(ISYMK)*(I - 1)
                     KOFF3 = KSCRZ + IOFFZ(ISYMD)
     &                     + NBAS(ISYMD)*NUMV*(I - 1)
     &                     + NBAS(ISYMD)*(IVEC - 1)

                     CALL DGEMV('N',NBAS(ISYMD),NRHF(ISYMK),
     &                          XMONE,XLAMDP(KOFF1),NTOTD,WORK(KOFF2),1,
     &                          ONE,WORK(KOFF3),1)

                  ENDDO

               ENDDO
            ENDDO

C           Read AO Cholesky vectors and store as L(g,d#J).
C           -----------------------------------------------

c      cnrm1 = 0.0d0

            DO IVEC = 1,NUMV

               JVEC = JVEC1 + IVEC - 1
               CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPTR,WORK(KREAD),LREAD)

c      cnrm1 = cnrm1 + ddot(n2bst(isycho),work(kchol),1,work(kchol),1)

               DO ISYMD = 1,NSYM

                  ISYMG = MULD2H(ISYMD,ISYCHO)

                  LENGD = NBAS(ISYMG)*NBAS(ISYMD)

                  KOFF1 = KCHOL + IAODIS(ISYMG,ISYMD)
                  KOFF2 = KSCRL + IOFFL(ISYMD) + LENGD*(IVEC - 1)

                  CALL DCOPY(LENGD,WORK(KOFF1),1,WORK(KOFF2),1)

               ENDDO

            ENDDO

c      cnrm1 = dsqrt(cnrm1)
c      cnrm2 = dsqrt(ddot(n2bst(isycho)*numv,work(kscrl),1,
c     &                                      work(kscrl),1))
c      write(lupri,*) '...Norm of Cholesky, read: ',cnrm1
c      write(lupri,*) '...Norm of Cholesky, sort: ',cnrm2

C           Calculate GIJ term in AO basis.
C           -------------------------------

            DO ISYMD = 1,NSYM

               ISYMG = MULD2H(ISYMD,ISYCHO)
               ISYMI = MULD2H(ISYMD,ISYMDI)

               KOFF1 = KSCRL  + IOFFL(ISYMD)
               KOFF2 = KSCRZ  + IOFFZ(ISYMD)
               KOFF3 = KSIGMA + IT1AO(ISYMG,ISYMI)

               LENDJ = NBAS(ISYMD)*NUMV

               NTOTG = MAX(NBAS(ISYMG),1)
               NTODJ = MAX(LENDJ,1)

               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),LENDJ,
     &                    ONE,WORK(KOFF1),NTOTG,WORK(KOFF2),NTODJ,
     &                    ONE,WORK(KOFF3),NTOTG)

            ENDDO

         ENDDO

C        Close MO Cholesky vector files.
C        -------------------------------

         CALL CHO_MOP(1,4,ISYCHO,LUCHIJ,1,1)
         IF (ISYCHO .EQ. ISYMX) CALL CHO_MOP(1,3,ISYCHO,LUCHAI,1,1)
         CALL CHO_MOP(1,1,ISYCHO,LUCHIA,1,1)

C        Close Y intermediate file.
C        --------------------------

         CALL CC_CYIOP(1,ISYCHO,-1)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

C     Transform the GIJ term to MO basis.
C     -----------------------------------

      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYMY)
         ISYMG = ISYMA

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFF1 = IGLMVI(ISYMG,ISYMA) + 1
         KOFF2 = KSIGMA + IT1AO(ISYMG,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDH(KOFF1),NTOTG,WORK(KOFF2),NTOTG,
     &              ONE,SIGMA1(KOFF3),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choatr0 */
      SUBROUTINE CC_CHOATR0(SIGMA1,X1AM,WORK,LWORK,ISYMX,NUMX,ISIDE)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate contributions from global (tot. sym.) E intermediates
C              to left-hand (ISIDE=-1) or right-hand (ISIDE=1) Jacobian
C              tranformations.
C
C     ISIDE = -1:
C
C        SIGMA1(ai) = sum(b) E(ba) * X(bi) - sum(j) X(aj) * E(ij)
C
C     ISIDE = 1:
C
C        SIGMA1(ai) = sum(b) E(ab) * X(bi) - sum(j) X(aj) * E(ji)
C
C     NOTICE: the content of SIGMA1 is overwritten!
C     
#include <implicit.h>
      DIMENSION SIGMA1(*), X1AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOATR0')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Read the E intermediates.
C     -------------------------

      KEIJ = 1
      KEAB = KEIJ  + NMATIJ(1)
      KEND = KEAB  + NMATAB(1)
      LWRK = LWORK - KEND + 1

      IF (LWRK .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CHO_IMOP(-1,1,LUE,1)
      CALL CHO_MOREAD(WORK(KEIJ),NMATIJ(1),1,1,LUE)
      CALL CHO_IMOP(1,1,LUE,1)

      CALL CHO_IMOP(-1,2,LUE,1)
      CALL CHO_MOREAD(WORK(KEAB),NMATAB(1),1,1,LUE)
      CALL CHO_IMOP(1,2,LUE,1)

      IF (LOCDBG) THEN
         EIJNRM = DSQRT(DDOT(NMATIJ(1),WORK(KEIJ),1,WORK(KEIJ),1))
         EABNRM = DSQRT(DDOT(NMATAB(1),WORK(KEAB),1,WORK(KEAB),1))
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': norm of EIJ: ',EIJNRM
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': norm of EAB: ',EABNRM
      ENDIF

C     Calculate contributions.
C     ------------------------

      IF (ISIDE .EQ. -1) THEN
         CALL CC_CHOATR0M1(SIGMA1,X1AM,WORK(KEIJ),WORK(KEAB),ISYMX,1,
     &                     NUMX)
      ELSE IF (ISIDE .EQ. 1) THEN
         CALL CC_CHOATR0P1(SIGMA1,X1AM,WORK(KEIJ),WORK(KEAB),ISYMX,1,
     &                     NUMX)
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,I10)')
     &   'ISIDE must be -1 or +1 in ',SECNAM,': ISIDE = ',ISIDE
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      RETURN
      END
C  /* Deck cc_choatr0m1 */
      SUBROUTINE CC_CHOATR0M1(SIGMA1,X1AM,EIJ,EAB,ISYMX,ISYME,NUMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate E-intermediate contributions, left-hand side style:
C
C              SIGMA1(ai) = sum(b) E(ba) * X(bi) - sum(j) X(aj) * E(ij)
C
C     NOTICE: the content of SIGMA1 is overwritten!
C
#include <implicit.h>
      DIMENSION SIGMA1(*), X1AM(*), EIJ(*), EAB(*)
#include <ccorb.h>
#include <ccsdsym.h>

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Get result symmetry.
C     --------------------

      ISYRES = MULD2H(ISYMX,ISYME)

      DO IX = 1,NUMX

         KOFX1 = NT1AM(ISYMX)*(IX - 1)  + 1
         KOFS1 = NT1AM(ISYRES)*(IX - 1) + 1

C        Virtual part.
C        -------------

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYMX)
            ISYMA = MULD2H(ISYMB,ISYME)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTB = MAX(NVIR(ISYMB),1)

            KOFFE = IMATAB(ISYMB,ISYMA) + 1
            KOFFX = KOFX1 + IT1AM(ISYMB,ISYMI)
            KOFFS = KOFS1 + IT1AM(ISYMA,ISYMI)

            CALL DGEMM('T','N',
     &                 NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                 ONE,EAB(KOFFE),NTOTB,X1AM(KOFFX),NTOTB,
     &                 ZERO,SIGMA1(KOFFS),NTOTA)

         ENDDO

C        Occupied part.
C        --------------

         DO ISYMJ = 1,NSYM

            ISYMA = MULD2H(ISYMJ,ISYMX)
            ISYMI = MULD2H(ISYMJ,ISYME)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTI = MAX(NRHF(ISYMI),1)

            KOFFX = KOFX1 + IT1AM(ISYMA,ISYMJ)
            KOFFE = IMATIJ(ISYMI,ISYMJ) + 1
            KOFFS = KOFS1 + IT1AM(ISYMA,ISYMI)

            CALL DGEMM('N','T',
     &                 NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                 XMONE,X1AM(KOFFX),NTOTA,EIJ(KOFFE),NTOTI,
     &                 ONE,SIGMA1(KOFFS),NTOTA)

         ENDDO

      ENDDO

      RETURN
      END
C  /* Deck cc_choatr0p1 */
      SUBROUTINE CC_CHOATR0P1(SIGMA1,X1AM,EIJ,EAB,ISYMX,ISYME,NUMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate E-intermediate contributions, right-hand side style:
C
C              SIGMA1(ai) = sum(b) E(ab) * X(bi) - sum(j) X(aj) * E(ji)
C
C     NOTICE: the content of SIGMA1 is overwritten!
C
#include <implicit.h>
      DIMENSION SIGMA1(*), X1AM(*), EIJ(*), EAB(*)
#include <ccorb.h>
#include <ccsdsym.h>

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Get result symmetry.
C     --------------------

      ISYRES = MULD2H(ISYMX,ISYME)

      DO IX = 1,NUMX

         KOFX1 = NT1AM(ISYMX)*(IX - 1)  + 1
         KOFS1 = NT1AM(ISYRES)*(IX - 1) + 1

C        Virtual part.
C        -------------

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYMX)
            ISYMA = MULD2H(ISYMB,ISYME)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTB = MAX(NVIR(ISYMB),1)

            KOFFE = IMATAB(ISYMA,ISYMB) + 1
            KOFFX = KOFX1 + IT1AM(ISYMB,ISYMI)
            KOFFS = KOFS1 + IT1AM(ISYMA,ISYMI)

            CALL DGEMM('N','N',
     &                 NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                 ONE,EAB(KOFFE),NTOTA,X1AM(KOFFX),NTOTB,
     &                 ZERO,SIGMA1(KOFFS),NTOTA)

         ENDDO

C        Occupied part.
C        --------------

         DO ISYMI = 1,NSYM

            ISYMJ = MULD2H(ISYMI,ISYME)
            ISYMA = MULD2H(ISYMJ,ISYMX)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTJ = MAX(NRHF(ISYMJ),1)

            KOFFX = KOFX1 + IT1AM(ISYMA,ISYMJ)
            KOFFE = IMATIJ(ISYMJ,ISYMI) + 1
            KOFFS = KOFS1 + IT1AM(ISYMA,ISYMI)

            CALL DGEMM('N','N',
     &                 NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                 XMONE,X1AM(KOFFX),NTOTA,EIJ(KOFFE),NTOTJ,
     &                 ONE,SIGMA1(KOFFS),NTOTA)

         ENDDO

      ENDDO

      RETURN
      END
C  /* Deck cc_chortri2j */
      SUBROUTINE CC_CHORTRI2J(FOCKD,XLAMDP,XLAMDH,SIGMA1,X1AM,
     &                        WORK,LWORK,ISYMX,NUMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the I2 and J terms for right-hand Jacobian
C              transformation.
C
C     The J term:
C
C        SIGMA1(ai) = SIGMA1(ai) + sum(ck) [2(ia|kc)-(ac|ka)] * X1AM(ck)
C
C     The I2 term:
C
C        SIGMA1(ai) = SIGMA1(ai)
C                   + sum(bj) [2t(ai,bj)-t(aj,bi)]
C                     * sum(ck) [2(ai|kc)-(ac|ki)] * X1AM(ck)
C
C     The calculation is done through construction of the Fock-matrix type
C     intermediate
C
C        W(al,be) = sum(ck) [2(al be|kc)-(al c|k be)] * X1AM(ck)
C
C     and, for the I2 term, through Cholesky decomposed 0th order doubles
C     amplitudes (as C intermediates in left-hand transformation). The
C     canonical orbital energies (FOCKD) are needed for this decomposition,
C     if it has not already been done.
C
#include <implicit.h>
      DIMENSION FOCKD(*), XLAMDP(*), XLAMDH(*), SIGMA1(*), X1AM(*)
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHORTRI2J')

C     Return if nothing to do.
C     ------------------------

      IF (NUMX .LE. 0) RETURN

C     Allocation.
C     -----------

      KLAMD2 = 1
      KWMAT  = KLAMD2 + MAX(NGLMRH(ISYMX),NT1AM(ISYMX))*NUMX
      KEND1  = KWMAT  + N2BST(ISYMX)*NUMX
      LWRK1  = LWORK  - KEND1 + 1

      KFJBP = KLAMD2

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: W'
         WRITE(LUPRI,'(5X,A,I10,A,I1,A)')
     &   'Number of trial vectors: ',NUMX,' (symmetry ',ISYMX,')'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Backtransform X1AM to obtain perturbed Lambda-hole matrices.
C     ------------------------------------------------------------

      CALL CC_CHOBCKTR(XLAMDH,X1AM,WORK(KLAMD2),ISYMX,NUMX)

C     Calculate perturbed Fock matrices (W intermediates).
C     ----------------------------------------------------

      CALL DZERO(WORK(KWMAT),N2BST(ISYMX)*NUMX)
      CALL CC_CHOFCKP(XLAMDP,WORK(KLAMD2),WORK(KWMAT),WORK(KEND1),LWRK1,
     &                1,ISYMX,NUMX)

C     Calculate J terms.
C     ------------------

      CALL CC_CHORTRJ(XLAMDP,XLAMDH,SIGMA1,WORK(KWMAT),WORK(KEND1),
     &                LWRK1,1,1,ISYMX,NUMX)

C     Calculate perturbed F(jb).
C     --------------------------

      CALL DZERO(WORK(KFJBP),NT1AM(ISYMX)*NUMX)
      CALL CC_CHORTRFJB(XLAMDP,XLAMDH,WORK(KFJBP),WORK(KWMAT),
     &                  WORK(KEND1),LWRK1,1,1,ISYMX,NUMX)

C     Calculate I2 terms.
C     -------------------

      KEND1 = KFJBP + NT1AM(ISYMX)*NUMX
      LWRK1 = LWORK - KEND1 + 1
      CALL CC_CHOCIM1(FOCKD,WORK(KFJBP),SIGMA1,WORK(KEND1),LWRK1,ISYMX,
     &                NUMX)

      RETURN
      END
C  /* Deck cc_chortrfjb */
      SUBROUTINE CC_CHORTRFJB(XLAMDP,XLAMDH,SIGMA1,WMAT,WORK,LWORK,
     &                        ISYMP,ISYMH,ISYMW,NUMW)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Transform NUMW AO matrices to MO basis as
C
C              SIGMA1(ia) = SIGMA1(ia)
C                         + sum(g,d) XLAMDH(g,a) * WMAT(d,g) * XLAMDP(d,i)
C
C              [with SIGMA1(ia) stored as ai]
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), SIGMA1(*), WMAT(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHORTRFJB')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Allocation.
C     -----------

      NEED = 0
      DO ISYMG = 1,NSYM
         ISYMD = MULD2H(ISYMG,ISYMW)
         ISYMI = MULD2H(ISYMD,ISYMP)
         NEEDS = NBAS(ISYMG)*NRHF(ISYMI)
         NEED  = MAX(NEED,NEEDS)
      ENDDO

      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',NEED,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate J terms.
C     ------------------

      ISYMIM = MULD2H(ISYMP,ISYMW)
      ISYMS  = MULD2H(ISYMH,ISYMIM)

      DO IW = 1,NUMW
         DO ISYMG = 1,NSYM

            ISYMD = MULD2H(ISYMG,ISYMW)
            ISYMI = MULD2H(ISYMD,ISYMP)
            ISYMA = MULD2H(ISYMG,ISYMH)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTD = MAX(NBAS(ISYMD),1)

            KOFFW = N2BST(ISYMW)*(IW - 1) + IAODIS(ISYMD,ISYMG) + 1
            KOFFP = IGLMRH(ISYMD,ISYMI)   + 1
            KOFFH = IGLMVI(ISYMG,ISYMA)   + 1
            KOFFS = NT1AM(ISYMS)*(IW - 1) + IT1AM(ISYMA,ISYMI)  + 1

            CALL DGEMM('T','N',NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMD),
     &                 ONE,WMAT(KOFFW),NTOTD,XLAMDP(KOFFP),NTOTD,
     &                 ZERO,WORK,NTOTG)
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &                 ONE,XLAMDH(KOFFH),NTOTG,WORK,NTOTG,
     &                 ONE,SIGMA1(KOFFS),NTOTA)

         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_chortrj */
      SUBROUTINE CC_CHORTRJ(XLAMDP,XLAMDH,SIGMA1,WMAT,WORK,LWORK,
     &                      ISYMP,ISYMH,ISYMW,NUMW)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate J term for right-hand Jacobian transformation
C              from the "perturbed" AO Fock matrix WMAT,
C
C              SIGMA1(ai) = SIGMA1(ai)
C                         + sum(g,d) XLAMDP(g,a) * WMAT(g,d) * XLAMDH(d,i)
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), SIGMA1(*), WMAT(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHORTRJ')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Allocation.
C     -----------

      NEED = 0
      DO ISYMD = 1,NSYM
         ISYMG = MULD2H(ISYMD,ISYMW)
         ISYMI = MULD2H(ISYMD,ISYMH)
         NEEDS = NBAS(ISYMG)*NRHF(ISYMI)
         NEED  = MAX(NEED,NEEDS)
      ENDDO

      IF (NEED .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',NEED,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate J terms.
C     ------------------

      ISYMIM = MULD2H(ISYMH,ISYMW)
      ISYMS  = MULD2H(ISYMP,ISYMIM)

      DO IW = 1,NUMW
         DO ISYMD = 1,NSYM

            ISYMG = MULD2H(ISYMD,ISYMW)
            ISYMI = MULD2H(ISYMD,ISYMH)
            ISYMA = MULD2H(ISYMG,ISYMP)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTD = MAX(NBAS(ISYMD),1)

            KOFFW = N2BST(ISYMW)*(IW - 1) + IAODIS(ISYMG,ISYMD) + 1
            KOFFH = IGLMRH(ISYMD,ISYMI)   + 1
            KOFFP = IGLMVI(ISYMG,ISYMA)   + 1
            KOFFS = NT1AM(ISYMS)*(IW - 1) + IT1AM(ISYMA,ISYMI)  + 1

            CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMD),
     &                 ONE,WMAT(KOFFW),NTOTG,XLAMDH(KOFFH),NTOTD,
     &                 ZERO,WORK,NTOTG)
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &                 ONE,XLAMDP(KOFFP),NTOTG,WORK,NTOTG,
     &                 ONE,SIGMA1(KOFFS),NTOTA)

         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_choaeffd */
      SUBROUTINE CC_CHOAEFFD(WORK,LWORK,MULSOL,TBAR)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Test RH and LH transformations by calculating
C              the effective CC2 Jacobian via transformations with
C              unit vectors.
C
C     NOTE: this routine requires a *lot* of memory!
C
#include <implicit.h>
      DIMENSION WORK(LWORK), TBAR(*)
      LOGICAL MULSOL
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOAEFFD')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      CALL AROUND(SECNAM//': Test of Jacobian tranformations')

C     Test memory.
C     ------------

      XLWORK = ONE*LWORK
      XNEED  = TWO*XT2SQ(1)
      IF (XNEED .GT. XLWORK) THEN
         WRITE(LUPRI,'(5X,A)')
     &   'Insufficient memory for test:'
         WRITE(LUPRI,'(5X,A,F15.1,A,/,5X,A,F15.1,A)')
     &   'Effective Jacobians alone require: ',XNEED,' words',
     &   'Total memory available is        : ',XLWORK,' words'
         WRITE(LUPRI,'(5X,A)')
     &   'Test is skipped!'
         GO TO 1000
      ENDIF

C     Allocation for effective Jacobians.
C     -----------------------------------

      KLEFT  = 1
      KRIGHT = KLEFT  + NT2SQ(1)
      KEND1  = KRIGHT + NT2SQ(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: Aeff'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Reduce printing.
C     ----------------

      IPRSAV = IPRINT
      IPRINT = -10

C     Left hand generation (row-wise).
C     --------------------------------

      ISIDE  = -1
      FREQ   = ZERO
      CALL CC_CHOAEFFD2(WORK(KLEFT),WORK(KEND1),LWRK1,FREQ,ISIDE)

c      CALL AROUND('Aeff from left-hand transformation (row-wise)')
c      CALL NOCC_PRT(WORK(KLEFT),1,'AIBJ')

C     Solve for multipliers.
C     ----------------------

      IF (MULSOL) THEN

         KETA  = KEND1
         KSCR  = KETA  + NT1AM(1)
         KTBAR = KSCR  + NT1AM(1)
         KEND  = KTBAR + NT1AM(1)
         LWRK  = LWORK - KEND + 1

         IF (LWRK .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for MULSOL in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Get right-hand side.
C        --------------------

         CALL CC_CHOETA(WORK(KETA),WORK(KEND),LWRK)
         CALL DCOPY(NT1AM(1),WORK(KETA),1,WORK(KTBAR),1)
         CALL DSCAL(NT1AM(1),XMONE,WORK(KTBAR),1)

C        Copy out sym. 1 block.
C        ----------------------

         KOFF1 = KLEFT + IT2SQ(1,1)
         LAIBJ = NT1AM(1)*NT1AM(1)
         CALL DCOPY(LAIBJ,WORK(KOFF1),1,WORK(KRIGHT),1)

C        Solve [Aeff]^T TBAR = -ETA.
C        ---------------------------

         INFO = -1
         CALL DGESOL(1,NT1AM(1),NT1AM(1),NT1AM(1),
     &               WORK(KRIGHT),WORK(KTBAR),WORK(KSCR),INFO)

         IF (INFO .NE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,I10,/)')
     &      SECNAM,': Error: DGESOL returned INFO = ',INFO
            CALL QUIT('DGESOL error in '//SECNAM)
         ENDIF

C        Check solution.
C        ---------------

         ETANRM = DSQRT(DDOT(NT1AM(1),WORK(KETA),1,WORK(KETA),1))
         TBRNRM = DSQRT(DDOT(NT1AM(1),WORK(KTBAR),1,WORK(KTBAR),1))

         WRITE(LUPRI,'(//,5X,A,A,1P,D22.15)')
     &   SECNAM,': norm of ETA : ',ETANRM
         WRITE(LUPRI,'(5X,A,A,1P,D22.15,/,5X,A)')
     &   SECNAM,': norm of TBAR: ',TBRNRM,
     &   '- going to calculate residual by trf. of solution vector...'

         CALL CC_CHORSDTST(WORK(KTBAR),WORK(KEND),LWRK)

         WRITE(LUPRI,'(5X,A)')
     &   '- going to calculate trf. vector by DGEMV using Aeff...'

         KOFF1 = KLEFT + IT2SQ(1,1)
         LAIBJ = NT1AM(1)*NT1AM(1)
         CALL DCOPY(LAIBJ,WORK(KOFF1),1,WORK(KRIGHT),1)
         CALL DGEMV('N',NT1AM(1),NT1AM(1),1.0D0,WORK(KRIGHT),NT1AM(1),
     &              WORK(KTBAR),1,0.0D0,WORK(KSCR),1)
         TRFNRM = DSQRT(DDOT(NT1AM(1),WORK(KSCR),1,WORK(KSCR),1))
         WRITE(LUPRI,'(5X,A,A,1P,D22.15,/)')
     &   SECNAM,': norm of trf. vector from DGEMV calc.: ',TRFNRM

         CALL DBGDIFAI(TBAR,WORK(KTBAR),TRGNRM,TSTNRM,ERRMIN,ERRMAX,1)
         WRITE(LUPRI,'(5X,A,A,1P,D22.15)')
     &   SECNAM,': norm of TBAR from conv. calc.: ',TRGNRM
         WRITE(LUPRI,'(5X,A,A,1P,D22.15)')
     &   SECNAM,': norm of TBAR from Chol. calc.: ',TBRNRM
         WRITE(LUPRI,'(5X,A,A,1P,D22.15)')
     &   SECNAM,': Min. difference              : ',ERRMIN
         WRITE(LUPRI,'(5X,A,A,1P,D22.15)')
     &   SECNAM,': Max. difference              : ',ERRMAX

C        Calculate density and properties.
C        ---------------------------------

         IF (NT2SQ(1) .GT. N2BST(1)) THEN

            IPRINT = IPRSAV

            KDEN = KRIGHT

            WRITE(LUPRI,'(/,5X,A,A)')
     &      SECNAM,': Calculating density matrix...'
            CALL CC_CHODEN(WORK(KTBAR),WORK(KDEN),WORK(KEND),LWRK)

            WRITE(LUPRI,'(5X,A,A)')
     &      SECNAM,': Calculating FOPs...'
            CALL CC_CHOFOP1(WORK(KDEN),WORK(KEND),LWRK,'CC2  ')

            IPRSAV = IPRINT
            IPRINT = -10

         ENDIF

      ENDIF

C     Transpose left-hand result.
C     ---------------------------

      CALL DCOPY(NT2SQ(1),WORK(KLEFT),1,WORK(KRIGHT),1)
      DO ISYMBJ = 1,NSYM
         ISYMAI = ISYMBJ
         DO NBJ = 1,NT1AM(ISYMBJ)
            KOFF1 = KRIGHT + IT2SQ(ISYMBJ,ISYMAI) + NBJ - 1
            KOFF2 = KLEFT  + IT2SQ(ISYMAI,ISYMBJ)
     &            + NT1AM(ISYMAI)*(NBJ - 1)
            CALL DCOPY(NT1AM(ISYMAI),WORK(KOFF1),NT1AM(ISYMBJ),
     &                               WORK(KOFF2),1)
         ENDDO
      ENDDO

c      CALL AROUND('Aeff from left-hand transformation (column-wise)')
c      CALL NOCC_PRT(WORK(KLEFT),1,'AIBJ')

C     Right hand generation (column-wise).
C     ------------------------------------

      ISIDE  = 1
      FREQ   = ZERO
      CALL CC_CHOAEFFD2(WORK(KRIGHT),WORK(KEND1),LWRK1,FREQ,ISIDE)

c      CALL AROUND('Aeff from right-hand transformation (column-wise)')
c      CALL NOCC_PRT(WORK(KRIGHT),1,'AIBJ')

C     Restore print level.
C     --------------------

      IPRINT = IPRSAV

C     Compare matrices obtained from LH and RH transformations.
C     ---------------------------------------------------------

      CALL DBGDIFAIBJ(WORK(KLEFT),WORK(KRIGHT),ALNRM,ARNRM,ERRMIN,
     &                ERRMAX,1,ISMNAI,ISMNBJ,MNAI,MNBJ,ISMXAI,ISMXBJ,
     &                MXAI,MXBJ,TARMIN,TARMAX)

      CALL DAXPY(NT2SQ(1),XMONE,WORK(KRIGHT),1,WORK(KLEFT),1)
      DIFNM2 = DDOT(NT2SQ(1),WORK(KLEFT),1,WORK(KLEFT),1)
      DIFNRM = DSQRT(DIFNM2)
      DIFRMS = DSQRT(DIFNM2/NT2SQ(1))

      CALL HEADER('Testing LH and RH eff. Jacobians',-1)
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Norm of left-hand  eff. Jacobian: ',ALNRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Norm of right-hand eff. Jacobian: ',ARNRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Difference of norms             : ',ALNRM-ARNRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Norm of difference              : ',DIFNRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Minimum absolute error          : ',ERRMIN
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Maximum absolute error          : ',ERRMAX
      WRITE(LUPRI,'(5X,A,1P,D22.15,/)')
     & 'RMS error                       : ',DIFRMS
      WRITE(LUPRI,'(5X,A)')
     & 'Location, min. abs. err.:'
      WRITE(LUPRI,'(5X,A,I10,A,I1)')
     & '   BJ = ',MNBJ,' of sym. ',ISMNBJ
      WRITE(LUPRI,'(5X,A,I10,A,I1)')
     & '   AI = ',MNAI,' of sym. ',ISMNAI
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & '   Value of LH Aeff: ',TARMIN
      WRITE(LUPRI,'(5X,A)')
     & 'Location, max. abs. err.:'
      WRITE(LUPRI,'(5X,A,I10,A,I1)')
     & '   BJ = ',MXBJ,' of sym. ',ISMXBJ
      WRITE(LUPRI,'(5X,A,I10,A,I1)')
     & '   AI = ',MXAI,' of sym. ',ISMXAI
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & '   Value of LH Aeff: ',TARMAX

c      CALL AROUND('L-R Aeff difference matrix')
c      CALL NOCC_PRT(WORK(KLEFT),1,'AIBJ')

 1000 CALL AROUND('End of '//SECNAM)

      RETURN
      END
C  /* Deck dbgdifaibj */
      SUBROUTINE DBGDIFAIBJ(TARGET,TEST,TRGNRM,TSTNRM,ERRMIN,ERRMAX,
     &                      ISYM,ISMNAI,ISMNBJ,MNAI,MNBJ,ISMXAI,ISMXBJ,
     &                      MXAI,MXBJ,TARMIN,TARMAX)
C
C     Thomas Bondo Pedersen, February 2003.
C
#include <implicit.h>
      DIMENSION TARGET(*), TEST(*)
#include <ccorb.h>
#include <ccsdsym.h>

      TRGNRM = DSQRT(DDOT(NT2SQ(1),TARGET,1,TARGET,1))
      TSTNRM = DSQRT(DDOT(NT2SQ(1),TEST,1,TEST,1))
      ERRMIN = 1.0D10
      ERRMAX = -1.0D10

      DO ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYM)
         DO NBJ = 1,NT1AM(ISYMBJ)
            DO NAI = 1,NT1AM(ISYMAI)
               NAIBJ = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ - 1)
     &               + NAI
               DIFF = DABS(TARGET(NAIBJ) - TEST(NAIBJ))
               IF (DIFF .LT. ERRMIN) THEN
                  ISMNAI = ISYMAI
                  ISMNBJ = ISYMBJ
                  MNAI   = NAI
                  MNBJ   = NBJ
                  ERRMIN = DIFF
                  TARMIN = TARGET(NAIBJ)
               ENDIF
               IF (DIFF .GT. ERRMAX) THEN
                  ISMXAI = ISYMAI
                  ISMXBJ = ISYMBJ
                  MXAI   = NAI
                  MXBJ   = NBJ
                  ERRMAX = DIFF
                  TARMAX = TARGET(NAIBJ)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_choaeffd2 */
      SUBROUTINE CC_CHOAEFFD2(AEFF,WORK,LWORK,FREQ,ISIDE)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Set up effective Jacobian by transformations
C              with unit vectors. For ISIDE = -1, AEFF is stored
C              transposed (i.e. row-wise).
C
#include <implicit.h>
      DIMENSION AEFF(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOAEFFD2')

      CHARACTER*5 ISITXT

      PARAMETER (ONE = 1.0D0)

      IF (ISIDE .EQ. -1) THEN
         ISITXT = 'left '
      ELSE IF (ISIDE .EQ. 1) THEN
         ISITXT = 'right'
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,I10,/)')
     &   'ISIDE must be -1 or +1 in ',SECNAM,':',
     &   'ISIDE = ',ISIDE
         CALL QUIT('ISIDE error in '//SECNAM)
      ENDIF

      WRITE(LUPRI,'(/,A,A,A,A,1P,D14.7,/)')
     & SECNAM,': eff. Jacobian transformations from ',ISITXT,
     & ' with FREQ=',FREQ

      DO ISYMTR = 1,NSYM

         KTRIAL = 1
         KEND1  = KTRIAL + NT1AM(ISYMTR)
         LWRK1  = LWORK  - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Trial'
            WRITE(LUPRI,'(5X,A,I1)')
     &      'Trial vector symmetry: ',ISYMTR
            WRITE(LUPRI,'(5X,A,I10)')
     &      'ISIDE = ',ISIDE
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         DO IVEC = 1,NT1AM(ISYMTR)

            WRITE(LUPRI,'(A,A,I10,A,I1,A,A)')
     &      SECNAM,': transforming unit vector number ',IVEC,
     &      ' of sym. ',ISYMTR,' from ',ISITXT

            CALL DZERO(WORK(KTRIAL),NT1AM(ISYMTR))
            WORK(KTRIAL+IVEC-1) = ONE
         
            KOFFA = IT2SQ(ISYMTR,ISYMTR) + NT1AM(ISYMTR)*(IVEC - 1) + 1
            CALL CC_CHOATR(AEFF(KOFFA),WORK(KTRIAL),FREQ,WORK(KEND1),
     &                     LWRK1,ISYMTR,1,ISIDE)

         ENDDO

      ENDDO

      RETURN
      END
