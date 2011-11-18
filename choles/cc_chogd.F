C  /* Deck cc_choeta */
      SUBROUTINE CC_CHOETA(ETA,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the effective right-hand side for
C              Cholesky CC2 0th order multiplier equation.
C              If the flag DSKETA is set (.TRUE.), then ETA
C              is assumed available on disk and will be read
C              instead. (The DSKETA flag is stored in chocc2.h)
C
C     Formula:
C
C        ETA(ai) = 2 * [ F(ia) + sum(Jb) L(ba,J) * Y(bi,J)
C                              - sum(Jj) Y(aj,J) * L(ij,J) ]
C
C     where
C
C        Y(ai,J) = sum(bj) tbar0(ai,bj) * L(bj,J)
C
C     and
C
C        tbar0(ai,bj) = - [2(ia|jb)-(ja|ib)]/[e(a)-e(i)+e(b)-e(j)]
C
C     (note the lacking factor of two in this definition of tbar0).
C
C     The Y intermediates are generated in CC_CYI and temporarily stored
C     in the files also used by the left-hand Jacobian transformations.
C
#include <implicit.h>
      DIMENSION ETA(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOETA')

      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0)

C     Set result symmetry.
C     --------------------

      ISYETA = 1

C     If requested, read ETA and return.
C     ----------------------------------

      IF (DSKETA) THEN

         CALL CHO_IMOP(-1,4,LUE,ISYETA)
         CALL CHO_MOREAD(ETA,NT1AM(ISYETA),1,1,LUE)
         CALL CHO_IMOP(1,4,LUE,ISYETA)

         RETURN

      ENDIF

C     Initialize ETA.
C     ---------------

      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(ETA,NT1AM(ISYETA),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

C     Read orbital energies.
C     ----------------------

      KFOCKD = 1
      KEND1  = KFOCKD + NORBTS
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: FOCK-diagonal'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &               .FALSE.,.TRUE.,.FALSE.)

C     Calculate Y intermediates.
C     --------------------------

c     ISIDE = -2
c     FREQ  = ZERO
c     CALL CC_CYI(WORK(KFOCKD),DUM1,DUM2,DUM3,WORK(KEND1),LWRK1,ISIDE,
c    &            ISYETA,1,FREQ,T2ANM,T2BNM,.FALSE.)

      FREQ = ZERO
      CALL CC_CYIETA(WORK(KFOCKD),WORK(KEND1),LWRK1,FREQ,T2ANM,T2BNM,
     &               .FALSE.,.FALSE.)

C     Calculate Lambda matrices.
C     --------------------------

      KLAMDP = 1
      KLAMDH = KLAMDP + NGLMDT(1)
      KT1AM  = KLAMDH + NGLMDT(1)
      KEND2  = KT1AM  + NT1AM(1)
      LWRK2  = LWORK  - KEND2 + 1

      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: Lambda'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND2-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND2),
     &            LWRK2)

C     Calculate contributions from Y intermediates, deleting files after use.
C     -----------------------------------------------------------------------

      KEND3 = KT1AM
      LWRK3 = LWORK - KEND3 + 1
      CALL CC_CHOETAY(WORK(KLAMDP),WORK(KLAMDH),ETA,WORK(KEND3),LWRK3,
     &                ISYETA,.TRUE.)

C     Scale by 2.
C     -----------

      CALL DSCAL(NT1AM(ISYETA),TWO,ETA,1)

C     Write to disk.
C     --------------

      CALL CHO_IMOP(-1,4,LUE,ISYETA)
      CALL CHO_MOWRITE(ETA,NT1AM(ISYETA),1,1,LUE)
      CALL CHO_IMOP(1,4,LUE,ISYETA)

C     Finally, set DSKETA flag.
C     -------------------------

      DSKETA = .TRUE.

      RETURN
      END
C  /* Deck cc_cyieta */
      SUBROUTINE CC_CYIETA(FOCKD,WORK,LWORK,FREQ,T2ANM,T2BNM,
     &                     DELP1,DELP2)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates for the effective right-hand
C              side for 0th order multipliers.
C
#include <implicit.h>
      DIMENSION FOCKD(*), WORK(LWORK)
      LOGICAL DELP1, DELP2
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>

C     Amplitude info.
C     ---------------

      FAC1   = 1.0D0
      FAC2   = 0.0D0
      SCD    = 1.0D0
      SCDG   = 1.0D0

      JTYP1  = 1
      JTYP2  = JTYP1
      ISYMP1 = 1
      ISYMP2 = ISYMP1
      NUMP12 = 1

      IOPTDN = 1
      IOPTCE = 1

C     Y info.
C     -------

      NUMZ  = 1
      JTYPZ = 3
      ISYMZ = 1
      JTYPY = -1

C     Calculate Y intermediates.
C     --------------------------

      CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NUMCHO,FAC1,NUMP12,
     &            DUMMY,IDUMMY,DUMMY,IDUMMY,FAC2,0,
     &            SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &            JTYPZ,ISYMZ,NUMCHO,NUMZ,JTYPY,
     &            DUMMY,IDUMMY,0,DUMMY,
     &            WORK,LWORK,T2ANM,T2BNM,DELP1,DELP2)

      RETURN
      END
C  /* Deck cc_choetay */
      SUBROUTINE CC_CHOETAY(XLAMDP,XLAMDH,ETA,WORK,LWORK,ISYMY,DELY)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate doubles contributions to the effective right-hand
C              side for 0th order multiplier equations.
C
C     Formula:
C
C     ETA(ai) = ETA(ai)
C             + sum(g) LamH(g,a) sum(Jd) L(gd,J) sum(b) LamP(d,b) * Y(bi,J)
C             - sum(Jj) Y(aj,J) * L(ij,J)
C
C     where g,d refer to AO basis and L(gd,J) = L(dg,J) has been used.
C
C     IF (DELY): Y intermediate files are deleted after use.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), ETA(*), WORK(LWORK)
      LOGICAL DELY
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER IOFFL(8), IOFFZ(8), IOFFY(8)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOETAY')

      PARAMETER (IOPTR = 2)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Set option for deleting or keeping Y intermediate files.
C     --------------------------------------------------------

      IOPTY = 1
      IF (DELY) IOPTY = 0

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

      KETA  = KEND0
      KEND1 = KETA  + NT1AO(ISYMY)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize.
C     -----------

      CALL DZERO(WORK(KETA),NT1AO(ISYMY))

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
C        SCRL: L(ij,#J),Y(b#J,i),L(g,d#J).
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

         CALL CC_CYIOP(-1,ISYCHO,-2)

C        Open occupied Cholesky vector file.
C        -----------------------------------

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

C           Read Y(bi,#J) in place of Z.
C           ----------------------------

            CALL CC_CYIRDF(WORK(KSCRZ),NUMV,JVEC1,ISYCHO,ISYMY,-2)

C           Read L(ij,#J) in place of L(g,d#J).
C           -----------------------------------

            CALL CHO_MOREAD(WORK(KSCRL),NMATIJ(ISYCHO),NUMV,JVEC1,
     &                      LUCHIJ)

C           Calculate occupied term.
C           ------------------------

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
     &                       ONE,ETA(KOFF3),NTOTA)

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
     &                    ONE,XLAMDP(KOFF1),NTOTD,WORK(KOFF2),NTOTB,
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

C           Calculate virtual term in AO basis.
C           -----------------------------------

            DO ISYMD = 1,NSYM

               ISYMG = MULD2H(ISYMD,ISYCHO)
               ISYMI = MULD2H(ISYMD,ISYMDI)

               KOFF1 = KSCRL + IOFFL(ISYMD)
               KOFF2 = KSCRZ + IOFFZ(ISYMD)
               KOFF3 = KETA  + IT1AO(ISYMG,ISYMI)

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

         CALL CHO_MOP(1,4,ISYCHO,LUCHIJ,1,1)

C        Close Y intermediate file, deleting if requested through IOPTY.
C        ---------------------------------------------------------------

         CALL CC_CYIOP(IOPTY,ISYCHO,-2)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

C     Transform the virtual term to MO basis.
C     ---------------------------------------

      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYMY)
         ISYMG = ISYMA

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFF1 = IGLMVI(ISYMG,ISYMA) + 1
         KOFF2 = KETA + IT1AO(ISYMG,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDH(KOFF1),NTOTG,WORK(KOFF2),NTOTG,
     &              ONE,ETA(KOFF3),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choxi0 */
      SUBROUTINE CC_CHOXI0(IXETRAN,NXETRAN,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the singles part of the conventional
C              right-hand side for the 1st order amplitude response
C              equations. Result is stored on disk, same place as
C              standard XI(0) vectors (i.e. using CC_WRRSP with
C              list 'O1 ').
C              This routine may be regarded as an "interface" to the
C              structure of the conventional CC code - it is called
C              from CCRHSVEC in place of CC_XIETA.
C
C     Input:
C
C     IXETRAN(1,*): operator indices in LBLOPR array.
C     IXETRAN(2,*): not used.
C     IXETRAN(3,*): indices for output XI vector on 'O1 ' list [FILXI list]
C     IXETRAN(4,*): not used.
C     IXETRAN(5,*): not used.
C     IXETRAN(6,*): not used.
C     IXETRAN(7,*): not used.
C     IXETRAN(8,*): not used.
C
C     NXETRAN : number of XI vectors/perturbation operators.
C
C     Formula:
C
C     XI(ai) = PROP(ai) + sum(bj) [2*t(ai,bj)-t(aj,bi)] * PROP(jb)
C
C     where the ground state doubles amplitudes are decomposed (if
C     it hasn't been done yet) as
C
C     t(ai,bj) = - sum(J) M(ai,J) * M(bj,J)
C
C     PROP designates the one-electron perturbation integrals, whose
C     precise nature (label,symmetry, etc.) is obtained through
C     IXETRAN(1,*).
C
#include <implicit.h>
#include <cclists.h>
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <chocc2.h>
#include <ccroper.h>
#include <dummy.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOXI0')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LPDBAS

      LOGICAL LONLYAO
      DATA LONLYAO /.TRUE./

      CHARACTER*8  LABEL
      CHARACTER*10 MODEL

      CHARACTER*3 FILXI
      DATA FILXI /'O1 '/

      INTEGER ILSTSYM  ! external function

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) ' Entering ',SECNAM
         WRITE(LUPRI,*) ' NXETRAN: ',NXETRAN
         WRITE(LUPRI,*) ' LWORK  : ',LWORK
         CALL FLSHFO(LUPRI)
      ENDIF

C     Initial tests.
C     --------------

      IF (.NOT. CHOINT) THEN
         WRITE(LUPRI,'(//,5X,A,A,/)')
     &   SECNAM,': Error: integrals not decomposed!'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      IF (CC2) THEN
         MODEL = 'CC2       '
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,/)')
     &   SECNAM,': Error: Only possible model is CC2'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Read orbital energies (if needed) and calculate Lambda matrices.
C     ----------------------------------------------------------------

      IF (CHOT2C) THEN
         LFOCKD = 0
      ELSE
         LFOCKD = NORBTS
      ENDIF

      KFOCKD = 1
      KLAMDP = KFOCKD + LFOCKD
      KLAMDH = KLAMDP + NGLMDT(1)
      KT1AM  = KLAMDH + NGLMDT(1)
      KEND   = KT1AM  + NT1AM(1)
      LWRK   = LWORK  - KEND + 1

      IF (LWRK .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: Lambda'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      IF (.NOT. CHOT2C) THEN
         CALL CHO_RDSIR(DUMMY,DUMMY,WORK(KFOCKD),DUMMY,WORK(KEND),LWRK,
     &                  .FALSE.,.TRUE.,.FALSE.)
      ENDIF
      IFAIL = -1
      CALL CHO_RDRST(DUMMY,WORK(KT1AM),DUMMY,.FALSE.,.TRUE.,.FALSE.,
     &               IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND),LWRK)

C     Reset memory (T1AM no longer needed).
C     -------------------------------------

      KEND1 = KT1AM
      LWRK1 = LWORK - KEND1 + 1

C     Loop over XI(0) vectors to compute.
C     -----------------------------------

      DO IVEC = 1,NXETRAN

C        Get operator and file indices.
C        ------------------------------

         IPROP = IXETRAN(1,IVEC)
         IFILE = IXETRAN(3,IVEC)

C        Extract and test info about operator.
C        -------------------------------------

         LABEL  = LBLOPR(IPROP)
         ISYMPR = ISYOPR(IPROP)
         IHERM  = ISYMAT(IPROP)
         LPDBAS = LPDBSOP(IPROP)
         ISYFIL = ILSTSYM(FILXI,IFILE)

         IF (LOCDBG) THEN
            WRITE(LUPRI,*) '   IVEC   = ',IVEC
            WRITE(LUPRI,*) '   IPROP  = ',IPROP
            WRITE(LUPRI,*) '   IFILE  = ',IFILE
            WRITE(LUPRI,*) '   LABEL  = ',LABEL
            WRITE(LUPRI,*) '   ISYMPR = ',ISYMPR
            WRITE(LUPRI,*) '   IHERM  = ',IHERM
            WRITE(LUPRI,*) '   LPDBAS = ',LPDBAS
            WRITE(LUPRI,*) '   ISYFIL = ',ISYFIL
            CALL FLSHFO(LUPRI)
         ENDIF

         IF (LPDBAS) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': WARNING:'
            WRITE(LUPRI,'(5X,A,A8,/,5X,A,//)')
     &      'It seems that the basis set depends on operator ',LABEL,
     &      '- this dependence will be ignored!'
         ENDIF

         IF ((ISYMPR.LT.1) .OR. (ISYMPR.GT.NSYM)) THEN
            WRITE(LUPRI,'(//,5X,A,A,A8,A,I2)')
     &      SECNAM,': Error: symmetry of operator ',LABEL,
     &      ' out of bounds:',ISYMPR
            CALL QUIT('Error in '//SECNAM)
         ENDIF

         IF (ISYMPR .NE. ISYFIL) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Symmetry mismatch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &      'Operator number  : ',IPROP,' (in LBLOPR list)',
     &      'Operator label   : ',LABEL
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'ISYOPR returned     : ',ISYMPR,
     &      '- expected (ILSTSYM): ',ISYFIL
            CALL QUIT('Symmetry mismatch in '//SECNAM)
         ENDIF

         IF ((IHERM.NE.-1) .AND. (IHERM.NE.1)) THEN
            WRITE(LUPRI,'(//,5X,A,A,A8,A,I2)')
     &      SECNAM,': Error: Hermiticity of operator ',LABEL,
     &      ' out of bounds:',IHERM
casm        CALL QUIT('Error in '//SECNAM)
         ENDIF

C        Allocation.
C        -----------

         KXI    = KEND1
         KPRPIA = KXI    + NT1AM(ISYMPR)
         KPROP  = KPRPIA + NT1AM(ISYMPR)
         KEND2  = KPROP  + N2BST(ISYMPR)
         LWRK2  = LWORK  - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: Xi'
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number  : ',IPROP,' (in LBLOPR list)',
     &      'Operator label   : ',LABEL,
     &      'Operator symmetry: ',ISYMPR
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read AO integrals and check info.
C        ---------------------------------

         IERR = -1
         CALL CCPRPAO(LABEL,LONLYAO,WORK(KPROP),IRREP,IREAL,IERR,
     &                WORK(KEND2),LWRK2)

         IF (IERR .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO'
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number  : ',IPROP,' (in LBLOPR list)',
     &      'Operator label   : ',LABEL,
     &      'Operator symmetry: ',ISYMPR
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERR
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERR .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': CCPRPAO failed to determine operator symmetry'
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number  : ',IPROP,' (in LBLOPR list)',
     &      'Operator label   : ',LABEL,
     &      'Operator symmetry: ',ISYMPR
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERR
            CALL QUIT('Irrep. error in CCPRPAO ('//SECNAM//')')
         ENDIF

         IF (IRREP .NE. ISYMPR) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Symmetry mismatch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &      'Operator number  : ',IPROP,' (in LBLOPR list)',
     &      'Operator label   : ',LABEL
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'CCPRPAO returned: ',IRREP,
     &      '- expected      : ',ISYMPR
            CALL QUIT('Symmetry mismatch in '//SECNAM)
         ENDIF

         IF (IREAL .NE. IHERM) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Hermiticity mismatch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number  : ',IPROP,' (in LBLOPR list)',
     &      'Operator label   : ',LABEL,
     &      'Operator symmetry: ',ISYMPR
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'CCPRPAO returned: ',IREAL,
     &      '- expected      : ',IHERM
casm        CALL QUIT('Hermiticity mismatch in '//SECNAM)
         ENDIF

C        Transform AO integrals to ai and ia.
C        (ai initializes Xi vector.)
C        ------------------------------------

         CALL CC_CHOTRFED(WORK(KLAMDP),WORK(KLAMDH),WORK(KPROP),
     &                    WORK(KXI),WORK(KPRPIA),
     &                    WORK(KEND2),LWRK2,1,1,ISYMPR)

C        Reset memory (AO property integrals no longer needed).
C        ------------------------------------------------------

         KEND3 = KPROP
         LWRK3 = LWORK - KEND3 + 1

C        Calculate the two first terms (i.e. the conventional
C        singles right-hand side). This is where the t-amplitudes
C        are decomposed, if it hasn't been done yet.
C        --------------------------------------------------------

         CALL CC_CHOCIM1(WORK(KFOCKD),WORK(KPRPIA),WORK(KXI),
     &                   WORK(KEND3),LWRK3,ISYMPR,1)

C        Write to disk.
C        --------------

         IOPT = 1
         CALL CC_WRRSP(FILXI,IFILE,ISYMPR,IOPT,MODEL,DUMMY,
     &                 WORK(KXI),DUMMY,WORK(KEND3),LWRK3)

         IF (LOCDBG) THEN
            XINRM = DSQRT(DDOT(NT1AM(ISYMPR),WORK(KXI),1,WORK(KXI),1))
            WRITE (LUPRI,*) 'Final result of ',SECNAM,':'
            WRITE (LUPRI,*) 'operator label:      ',LBLOPR(IPROP)
            WRITE (LUPRI,*) 'output file type:    ',FILXI(1:3)
            WRITE (LUPRI,*) 'index of output file:',IFILE
            WRITE (LUPRI,*) 'pert. dep. basis:    ',LPDBAS
            WRITE(LUPRI,*) 'Norm of XI: ',XINRM
c            CALL CC_PRP(WORK(KXI),DUMMY,ISYMPR,1,0)
            CALL NOCC_PRT(WORK(KXI),ISYMPR,'AI  ')
            CALL FLSHFO(LUPRI)
         END IF

      ENDDO

      RETURN
      END
C  /* Deck cc_choxi1 */
      SUBROUTINE CC_CHOXI1(XI,WORK,LWORK,ILSXI0,LABEL,FREQ,ISYMPR,IALG)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the Cholesky CC2 effective right-hand side for
C              1st order amplitude response equations for perturbation
C              operator LABEL and frequency FREQ.
C
C     Input:
C     WORK   is work space, dimension: LWORK.
C     ILSXI0 is the file list ('O1 ') number to read the conventional
C            singles part (calculated in CC_CHOXI0).
C     LABEL  is the perturbation operator label.
C     FREQ   is the perturbation frequency.
C     ISYMPR is the symmetry of the perturbation operator (and XI).
C     IALG   specifies the algorithm to use when calculating the 
C            double contribution:
C            IALG = 2: Force in-core storage of 2 full square doubles arrays.
C            IALG = 3: Force virtual-batching algorithm (using CC_CYI).
C            For any other value of IALG, the routine decides which algorithm
C            to use based on work space availability and the (very reasonable)
C            assumption that to store the doubles is fastest.
C
C     Output:
C     XI: effective right-hand side; need not be initialized on input
C         (whatever is stored in XI on input will be overwritten).
C
C
C     Formula (the conventional singles part of XI is initially
C     read from disk):
C
C     XI(ai) = XI(ai)
C
C            + sum(bj)  [2*XI2(ai,bj)-XI2(aj,bi)] * FIA(jb)
C
C            + sum(ckb) [2*XI2(ai,bj)-XI2(aj,bi)] * (kc|ab)
C
C            - sum(ckj) [2*XI2(ai,bj)-XI2(aj,bi)] * (kc|ji)
C
C     where [e(ai,bj) is the usual orbital energy denominator]
C
C     XI2(ai,bj) = - P(ai,bj) [ sum(c) t(ai,cj) * PROP(bc)
C                             - sum(k) t(ai,bk) * PROP(ji)]/[e(ai,bj) - FREQ]
C
C     and the ground state doubles amplitudes *must* have been decomposed as
C
C     t(ai,bj) = - sum(J) M(ai,J) * M(bj,J)
C
C     PROP designates the MO property integrals of operator LABEL
C     (calculated in this routine), and FIA is the unperturbed de-excitation
C     block of the (T1-transformed) Fock matrix.
C
C     TODO: the "doubles-direct" algorithm may be improved by doing things
C           much the same way as the in-core algorithm. I.e. set up T2AM,
C           transform with property integrals, calculate contributions.
C
#include <implicit.h>
      DIMENSION XI(*), WORK(LWORK)
      CHARACTER*8 LABEL
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chocc2.h>
#include <dummy.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOXI1')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LONLYAO
      DATA LONLYAO /.TRUE./

      CHARACTER*10 MODEL

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Set MODEL.
C     ----------

      IF (CHOINT .AND. CC2) THEN
         MODEL = 'CC2       '
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Error in ',SECNAM,': inconsistent flags:'
         WRITE(LUPRI,'(5X,A,L3)')   'CHOINT: ',CHOINT
         WRITE(LUPRI,'(5X,A,L3,/)') 'CC2   : ',CC2
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Debug print.
C     ------------

      IF (LOCDBG) THEN
         TIMT = SECOND()
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Entering ',SECNAM
         WRITE(LUPRI,*) 'LWORK  = ',LWORK
         WRITE(LUPRI,*) 'ILSXI0 = ',ILSXI0
         WRITE(LUPRI,*) 'LABEL  = ',LABEL
         WRITE(LUPRI,*) 'FREQ   = ',FREQ
         WRITE(LUPRI,*) 'ISYMPR = ',ISYMPR
         WRITE(LUPRI,*) 'IALG   = ',IALG
      ENDIF

C     Read the conventional XI from disk.
C     -----------------------------------

      IOPT = 1
      CALL CC_RDRSP('O1 ',ILSXI0,ISYMPR,IOPT,MODEL,XI,DUMMY)

      IF (LOCDBG) THEN
         XICNM = DSQRT(DDOT(NT1AM(ISYMPR),XI,1,XI,1))
         WRITE(LUPRI,*) 'XIcnm = ',XICNM
      ENDIF

C     Allocation.
C     -----------

      KLAMDP = 1
      KLAMDH = KLAMDP + NGLMDT(1)
      KPRPIJ = KLAMDH + NGLMDT(1)
      KPRPAB = KPRPIJ + NMATIJ(ISYMPR)
      KFIA   = KPRPAB + NMATAB(ISYMPR)
      KFOCKD = KFIA   + NT1AM(1)
      KEND1  = KFOCKD + NORBTS
      LWRK1  = LWORK  - KEND1 + 1

      KPRPAO = KFIA
      KEND0  = KPRPAO + N2BST(ISYMPR)
      LWRK0  = LWORK  - KEND0 + 1      ! work space for AO read/MO trf.

      KT1AM  = KFIA

      KEND = MAX(KEND0,KEND1)
      LWRK = LWORK - KEND + 1

      IF (LWRK .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: scr'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate Lambda matrices.
C     --------------------------

      IFAIL = -1
      CALL CHO_RDRST(DUMMY,WORK(KT1AM),DUMMY,.FALSE.,.TRUE.,.FALSE.,
     &               IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     &            LWRK1)

C     Read AO property integrals.
C     ---------------------------

      IERR = -1
      CALL CCPRPAO(LABEL,LONLYAO,WORK(KPRPAO),IRREP,IREAL,IERR,
     &             WORK(KEND0),LWRK0)
            
      IF (IERR .GT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': I/O error ocurred in CCPRPAO'
         WRITE(LUPRI,'(5X,A,A8)')
     &   'Operator: ',LABEL
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'CCPRPAO returned IERR = ',IERR 
         CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
      ELSE IF (IERR .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': CCPRPAO failed to determine operator symmetry'
         WRITE(LUPRI,'(5X,A,A8)')
     &   'Operator: ',LABEL
         WRITE(LUPRI,'(5X,A,I10,/)')
     &   'CCPRPAO returned IERR = ',IERR
         CALL QUIT('Irrep. error in CCPRPAO ('//SECNAM//')')
      ENDIF

      IF (IRREP .NE. ISYMPR) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Symmetry mismatch in ',SECNAM
         WRITE(LUPRI,'(5X,A,A8)')
     &   'Operator: ',LABEL
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'CCPRPAO returned: ',IRREP,
     &   '- expected      : ',ISYMPR
         CALL QUIT('Symmetry mismatch in '//SECNAM)
      ENDIF

C     Transform to MO basis (occupied and virtual blocks).
C     ----------------------------------------------------

      CALL CC_CHOTRFOV(WORK(KLAMDP),WORK(KLAMDH),WORK(KPRPAO),
     &                 WORK(KPRPIJ),WORK(KPRPAB),WORK(KEND0),LWRK0,
     &                 1,1,ISYMPR)

C     Read F(ia) and orbital energies.
C     --------------------------------

      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(WORK(KFIA),NT1AM(1),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

      CALL CHO_RDSIR(DUMMY,DUMMY,WORK(KFOCKD),DUMMY,WORK(KEND1),
     &               LWRK1,.FALSE.,.TRUE.,.FALSE.)

C     Calculate the frequency-dependent terms.
C     If algorithm is not forced from input, decide algorithm:
C     If memory allows, generate full-square t-amplitudes and transform
C     to obtain XI2. Else, transform amplitude Cholesky vectors and
C     use CC_CYI to obtain Y intermediates.
C     -----------------------------------------------------------------

      IF ((IALG .LE. 1) .OR. (IALG .GE. 4)) THEN
         XCHND = ZERO
         DO ISYCHO = 1,NSYM
            IF (NUMCHO(ISYCHO) .GT. 10) THEN
               MINVC = NUMCHO(ISYCHO)/10
            ELSE
               MINVC = MIN(NUMCHO(ISYCHO),5)
            ENDIF
            XMNVC = ONE*MINVC
            X2BST = ONE*N2BST(ISYCHO)
            XCHO  = XT1AO(ISYMPR) + X2BST
     &            + XMNVC*(XT1AO(ISYCHO) + X2BST)
            XCHND = MAX(XCHND,XCHO)
         ENDDO
         XWRK1 = ONE*LWRK1
         XNEED = MAX(XT2SQ(1),XT2SQ(ISYMPR)) + XT2SQ(ISYMPR) + XCHND
         IF (XNEED .LT. XWRK1) THEN
            LALG = 2
         ELSE
            LALG = 3
         ENDIF
      ELSE
         LALG = IALG
      ENDIF

      IF (LALG .EQ. 2) THEN
         CALL CC_CHOXI2(WORK(KFOCKD),WORK(KFIA),WORK(KLAMDP),
     &                  WORK(KLAMDH),XI,WORK(KPRPIJ),WORK(KPRPAB),
     &                  FREQ,WORK(KEND1),LWRK1,1,ISYMPR)
      ELSE IF (LALG .EQ. 3) THEN
         CALL CC_CHOXI3(WORK(KFOCKD),WORK(KFIA),WORK(KLAMDP),
     &                  WORK(KLAMDH),XI,WORK(KPRPIJ),WORK(KPRPAB),
     &                  FREQ,WORK(KEND1),LWRK1,1,ISYMPR)
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,I10)')
     &   'Logical error in ',SECNAM,': LALG = ',LALG
         CALL QUIT('Logical error in '//SECNAM)
      ENDIF

C     Debug print.
C     ------------

      IF (LOCDBG) THEN
         TIMT = SECOND() - TIMT
         XINRM = DSQRT(DDOT(NT1AM(ISYMPR),XI,1,XI,1))
         WRITE(LUPRI,*) 'Exiting ',SECNAM
         WRITE(LUPRI,*) 'LALG  = ',LALG
         WRITE(LUPRI,*) 'XInrm = ',XINRM
         WRITE(LUPRI,*) 'TIMT  = ',TIMT
      ENDIF

      RETURN
      END
C  /* Deck cc_choxi2 */
      SUBROUTINE CC_CHOXI2(FOCKD,FIA,XLAMDP,XLAMDH,XI,PROPIJ,PROPAB,
     &                     FREQ,WORK,LWORK,NFREQ,ISYMPR)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate frequency-dependent terms of the right-hand side
C              for 1st order amplitude equations.
C
C     See CC_CHOXI1 for the explicit formula.
C
C     WARNING: This routine allocates 2 full square doubles arrays. Use
C              CC_CHOXI3 for an implementation less demanding in terms of
C              memory (but more demanding in terms of CPU time).
C
#include <implicit.h>
      DIMENSION FOCKD(*), FIA(*), XLAMDP(*), XLAMDH(*), XI(*)
      DIMENSION PROPIJ(*), PROPAB(*)
      DIMENSION FREQ(NFREQ), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <ciarc.h>
#include <chocc2.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOXI2')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         FIANRM = DSQRT(DDOT(NT1AM(1),FIA,1,FIA,1))
         XLPNRM = DSQRT(DDOT(NGLMDT(1),XLAMDP,1,XLAMDP,1))
         XLHNRM = DSQRT(DDOT(NGLMDT(1),XLAMDH,1,XLAMDH,1))
         XINRM  = DSQRT(DDOT(NT1AM(ISYMPR)*NFREQ,XI,1,XI,1))
         PIJNRM = DSQRT(DDOT(NMATIJ(ISYMPR),PROPIJ,1,PROPIJ,1))
         PABNRM = DSQRT(DDOT(NMATAB(ISYMPR),PROPAB,1,PROPAB,1))
         CALL AROUND('Debug Info Entering '//SECNAM)
         WRITE(LUPRI,'(A,I1)')
     &   'Symmetry of property integrals: ',ISYMPR
         WRITE(LUPRI,'(A,I10,/,A)')
     &   'Number of frequencies: ',NFREQ,' - values:'
         WRITE(LUPRI,'(1P,4D20.10)') (FREQ(IFREQ), IFREQ = 1,NFREQ)
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of F(ia)  : ',FIANRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDP : ',XLPNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDH : ',XLHNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XI     : ',XINRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPIJ : ',PIJNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPAB : ',PABNRM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Allocation. T2AM MUST BE FIRST!
C     -------------------------------

      KT2AM = 1
      KX2AM = KT2AM + MAX(NT2SQ(1),NT2SQ(ISYMPR))
      KEND1 = KX2AM + NT2SQ(ISYMPR)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: doubles'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ELSE IF (LWRK1 .GT. LWORK) THEN
         WRITE(LUPRI,'(5X,A,A,A)')
     &   'Allocation error in ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10)')
     &   'Apparent need     : ',KEND1-1,
     &   'Available         : ',LWORK
         WRITE(LUPRI,'(5X,A,/)')
     &   '- problem size is probably too large: use CC_CHOXI3 instead!'
         CALL QUIT('Allocation error in '//SECNAM)
      ENDIF

C     Calculate doubles amplitudes (decomposed).
C     ==========================================

      IF (.NOT. CHOT2C) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/)')
     &   'Error in ',SECNAM,
     &   ': ground state amplitudes must be decomposed!'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      ISYCH1 = 1
      ISYCH2 = ISYCH1
      ITYP1  = 6
      ITYP2  = ITYP1
      DO ISYCHO = 1,NSYM
         NTOVEC(ISYCHO) = NCHMOC(ISYCHO)
      ENDDO
      DO ISYM = 1,NSYM
         IOFA1(ISYM)  = 1
         IOFB1(ISYM)  = 1
         LVIRA(ISYM)  = NVIR(ISYM)
         LVIRB(ISYM)  = NVIR(ISYM)
         NX1AMA(ISYM) = NT1AM(ISYM)
         NX1AMB(ISYM) = NT1AM(ISYM)
         DO ISYM1 = 1,NSYM
            IX1AMA(ISYM1,ISYM) = IT1AM(ISYM1,ISYM)
            IX1AMB(ISYM1,ISYM) = IT1AM(ISYM1,ISYM)
            IX2SQ(ISYM1,ISYM)  = IT2SQ(ISYM1,ISYM)
         ENDDO
      ENDDO
      NX2SQ  = NT2SQ(1)
      LIDEN  = .TRUE.
      SYMTRZ = .FALSE.
      CIAMIO = .FALSE.
      GETMNM = .FALSE.
      INXINT = .TRUE.
      IPRCIA = IPRINT - 30 + IPRLVL
      SCD    = ONE
      NPASS  = 0
      KLAST  = 0

      KEND = KX2AM
      LWRK = LWORK - KEND + 1
      CALL CC_CIA(WORK(KT2AM),WORK(KEND),LWRK,SCD,KLAST,NPASS)
      CALL DSCAL(NT2SQ(1),XMONE,WORK(KT2AM),1)

      IF (LOCDBG) THEN
         T2NRM = ZERO
         CALL CC_CYINRM(WORK(KT2AM),1,T2NRM)
         T2NRM = DSQRT(T2NRM)
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': Packed T2AM norm: ',T2NRM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Transform with property integrals to obtain standard rhs.
C     ---------------------------------------------------------

      CALL CC_CHOXI21(WORK(KT2AM),WORK(KX2AM),PROPIJ,PROPAB,1,ISYMPR)

C     Symmetrize and move result to T2AM (no longer needed).
C     ------------------------------------------------------

      CALL CC_CIASMT(WORK(KX2AM),ISYMPR)
      CALL DCOPY(NT2SQ(ISYMPR),WORK(KX2AM),1,WORK(KT2AM),1)

      IF (LOCDBG) THEN
         X2NRM = ZERO
         CALL CC_CYINRM(WORK(KX2AM),ISYMPR,X2NRM)
         X2NRM = DSQRT(X2NRM)
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': Packed X2i  norm: ',X2NRM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Re-order memory indices.
C     ------------------------

      KX2AM  = 1
      KXI2   = KX2AM  + NT2SQ(ISYMPR)
      KEND1  = KXI2   + NT2SQ(ISYMPR)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LOCDBG) THEN
         X2NRM = ZERO
         CALL CC_CYINRM(WORK(KX2AM),ISYMPR,X2NRM)
         X2NRM = DSQRT(X2NRM)
         WRITE(LUPRI,'(A,A,1P,D22.15,A)')
     &   SECNAM,': Packed X2i  norm: ',X2NRM,' - after mem. reorder'
         CALL FLSHFO(LUPRI)
      ENDIF

C     Calculate contributions in loop over frequencies.
C     =================================================

      DO IFREQ = 1,NFREQ

C        Offset to result vector.
C        ------------------------

         KOFXI = NT1AM(ISYMPR)*(IFREQ - 1) + 1

C        Copy X2AM to XI2.
C        -----------------

         CALL DCOPY(NT2SQ(ISYMPR),WORK(KX2AM),1,WORK(KXI2),1)

C        Divide by appropriate denominators.
C        N.B. CC_DNOM uses the index arrays of ciarc.h
C             - should not be a problem (if set correctly
C             above).
C        ------------------------------------------------

         CALL CC_DNOM(FOCKD,WORK(KXI2),FREQ(IFREQ),ISYMPR)

         IF (LOCDBG) THEN
            X2NRM = ZERO
            CALL CC_CYINRM(WORK(KXI2),ISYMPR,X2NRM)
            X2NRM = DSQRT(X2NRM)
            WRITE(LUPRI,'(A,A,1P,D22.15)')
     &      SECNAM,': Packed X2a  norm: ',X2NRM
            WRITE(LUPRI,'(A,I10)')
     &      ' - for frequency # ',IFREQ
            CALL FLSHFO(LUPRI)
         ENDIF

C        Set up 2CME amplitudes.
C        -----------------------

         CALL CC_CYITCME(WORK(KXI2),WORK(KEND1),LWRK1,ISYMPR,IDUM)

         IF (LOCDBG) THEN
            X2NRM = ZERO
            CALL CC_CYINRM(WORK(KXI2),ISYMPR,X2NRM)
            X2NRM = DSQRT(X2NRM)
            WRITE(LUPRI,'(A,A,1P,D22.15,A)')
     &      SECNAM,': Packed X2AM norm: ',X2NRM,' - i.e. 2CME'
            WRITE(LUPRI,'(A,I10)')
     &      ' - for frequency # ',IFREQ
            CALL FLSHFO(LUPRI)
         ENDIF

C        Calculate I-type term:
C        XI(ai) = XI(ai) + sum(bj) [2*XI2(ai,bj)-XI2(aj,bi)]*FIA(jb).
C        ------------------------------------------------------------

         ISYMBJ = 1
         ISYMAI = ISYMPR

         NTOAI = MAX(NT1AM(ISYMAI),1)

         CALL DGEMV('N',NT1AM(ISYMAI),NT1AM(ISYMBJ),
     &              ONE,WORK(KXI2),NTOAI,FIA,1,
     &              ONE,XI(KOFXI),1)

         IF (LOCDBG) THEN
            XINRM = DSQRT(DDOT(NT1AM(ISYMPR),XI(KOFXI),1,XI(KOFXI),1))
            WRITE(LUPRI,'(A,A,1P,D22.15)')
     &      SECNAM,': Norm XI (I-term): ',XINRM
            WRITE(LUPRI,'(A,I10)')
     &      ' - for frequency # ',IFREQ
            CALL FLSHFO(LUPRI)
         ENDIF

C        Calculate contributions from Y-type intermediates.
C        --------------------------------------------------

         CALL CC_CHOXI22(XLAMDP,XLAMDH,XI(KOFXI),WORK(KXI2),WORK(KEND1),
     &                   LWRK1,ISYMPR)

         IF (LOCDBG) THEN
            XINRM = DSQRT(DDOT(NT1AM(ISYMPR),XI(KOFXI),1,XI(KOFXI),1))
            WRITE(LUPRI,'(A,A,1P,D22.15)')
     &      SECNAM,': Norm XI (Y con.): ',XINRM
            WRITE(LUPRI,'(A,I10)')
     &      ' - for frequency # ',IFREQ
            CALL FLSHFO(LUPRI)
         ENDIF

      ENDDO

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         FIANRM = DSQRT(DDOT(NT1AM(1),FIA,1,FIA,1))
         XLPNRM = DSQRT(DDOT(NGLMDT(1),XLAMDP,1,XLAMDP,1))
         XLHNRM = DSQRT(DDOT(NGLMDT(1),XLAMDH,1,XLAMDH,1))
         XINRM  = DSQRT(DDOT(NT1AM(ISYMPR)*NFREQ,XI,1,XI,1))
         PIJNRM = DSQRT(DDOT(NMATIJ(ISYMPR),PROPIJ,1,PROPIJ,1))
         PABNRM = DSQRT(DDOT(NMATAB(ISYMPR),PROPAB,1,PROPAB,1))
         CALL AROUND('Debug Info Exiting '//SECNAM)
         WRITE(LUPRI,'(A,I1)')
     &   'Symmetry of property integrals: ',ISYMPR
         WRITE(LUPRI,'(A,I10,/,A)')
     &   'Number of frequencies: ',NFREQ,' - values:'
         WRITE(LUPRI,'(1P,4D20.10)') (FREQ(IFREQ), IFREQ = 1,NFREQ)
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of F(ia)  : ',FIANRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDP : ',XLPNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDH : ',XLHNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XI     : ',XINRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPIJ : ',PIJNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPAB : ',PABNRM
         CALL FLSHFO(LUPRI)
      ENDIF

      RETURN
      END
C  /* Deck cc_choxi3 */
      SUBROUTINE CC_CHOXI3(FOCKD,FIA,XLAMDP,XLAMDH,XI,PROPIJ,PROPAB,
     &                     FREQ,WORK,LWORK,NFREQ,ISYMPR)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate frequency-dependent terms of the right-hand side
C              for 1st order amplitude equations.
C
C     See CC_CHOXI1 for the explicit formula.
C
#include <implicit.h>
      DIMENSION FOCKD(*), FIA(*), XLAMDP(*), XLAMDH(*), XI(*)
      DIMENSION PROPIJ(*), PROPAB(*)
      DIMENSION FREQ(NFREQ), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <ciarc.h>
#include <chocc2.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOXI3')

      LOGICAL DELCH2, DELY

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         FIANRM = DSQRT(DDOT(NT1AM(1),FIA,1,FIA,1))
         XLPNRM = DSQRT(DDOT(NGLMDT(1),XLAMDP,1,XLAMDP,1))
         XLHNRM = DSQRT(DDOT(NGLMDT(1),XLAMDH,1,XLAMDH,1))
         XINRM  = DSQRT(DDOT(NT1AM(ISYMPR)*NFREQ,XI,1,XI,1))
         PIJNRM = DSQRT(DDOT(NMATIJ(ISYMPR),PROPIJ,1,PROPIJ,1))
         PABNRM = DSQRT(DDOT(NMATAB(ISYMPR),PROPAB,1,PROPAB,1))
         CALL AROUND('Debug Info Entering '//SECNAM)
         WRITE(LUPRI,'(A,I1)')
     &   'Symmetry of property integrals: ',ISYMPR
         WRITE(LUPRI,'(A,I10,/,A)')
     &   'Number of frequencies: ',NFREQ,' - values:'
         WRITE(LUPRI,'(1P,4D20.10)') (FREQ(IFREQ), IFREQ = 1,NFREQ)
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of F(ia)  : ',FIANRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDP : ',XLPNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDH : ',XLHNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XI     : ',XINRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPIJ : ',PIJNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPAB : ',PABNRM
      ENDIF

C     Set stuff for Y intermediate calculation.
C     Y intermediates are stored in same files
C     as the RH transformation Y intermediates.
C     -----------------------------------------

      ISIDE  = 3
      DELCH2 = .FALSE.

C     Set stuff for G- and H-type terms.
C     ----------------------------------

      JSIDE = 1
      DELY  = .FALSE.

C     Set type specifiers and symmetry for Cholesky vector transformation.
C     The transformed vectors are stored in the same files as the perturbed
C     vectors in the RH transformation.
C     ---------------------------------------------------------------------

      JTYP1  = 6
      JTYP2  = 9
      JSYCH1 = 1

C     Calculate contributions in loop over frequencies.
C     -------------------------------------------------

      DO IFREQ = 1,NFREQ

C        Offset to result vector for this frequency.
C        -------------------------------------------

         KOFXI = NT1AM(ISYMPR)*(IFREQ - 1) + 1

C        Transform Cholesky vectors from T2 decomposition
C        with property integrals. This is where the minus
C        sign from the decomposition is incorporated.
C        ------------------------------------------------

         CALL CC_CHOTRPRP(PROPIJ,PROPAB,WORK,LWORK,ISYMPR,JTYP1,JTYP2,
     &                    JSYCH1,NCHMOC)

C        Get Y intermediates, calculating the I-type contribution
C        as a byproduct. Delete transformed Cholesky vector files
C        after use in the last frequency.
C        --------------------------------------------------------

c        IF (IFREQ .EQ. NFREQ) DELCH2 = .TRUE.
c        CALL CC_CYI(FOCKD,XI(KOFXI),DUM1,FIA,WORK,LWORK,
c    &               ISIDE,ISYMPR,1,FREQ(IFREQ),XNRM,XI2NRM,DELCH2)

         IF (IFREQ .EQ. NFREQ) DELCH2 = .TRUE.
         CALL CC_CYIXI(FOCKD,XI(KOFXI),FIA,WORK,LWORK,FREQ(IFREQ),
     &                 ISYMPR,XNRM,XI2NRM,.FALSE.,DELCH2)

         IF (LOCDBG) THEN
            XINRM = DSQRT(DDOT(NT1AM(ISYMPR),XI(KOFXI),1,XI(KOFXI),1))
            WRITE(LUPRI,'(A,A,1P,D22.15)')
     &      SECNAM,': Norm XI (I-term): ',XINRM
            WRITE(LUPRI,'(A,I10)')
     &      ' - for frequency # ',IFREQ
         ENDIF

C        Calculate contributions from Y intermediates. Delete files
C        with Y intermediates after use in the last frequency.
C        ----------------------------------------------------------

         IF (IFREQ .EQ. NFREQ) DELY = .TRUE.
         CALL CC_CHOATR1(XLAMDP,XLAMDH,XI(KOFXI),DUM1,WORK,LWORK,ISYMPR,
     &                   JSIDE,DELY,IFREQ)

         IF (LOCDBG) THEN
            XINRM = DSQRT(DDOT(NT1AM(ISYMPR),XI(KOFXI),1,XI(KOFXI),1))
            WRITE(LUPRI,'(A,A,1P,D22.15)')
     &      SECNAM,': Norm XI (Y con.): ',XINRM
            WRITE(LUPRI,'(A,I10)')
     &      ' - for frequency # ',IFREQ
         ENDIF

      ENDDO

C     Debug: print.
C     -------------

      IF (LOCDBG) THEN
         FIANRM = DSQRT(DDOT(NT1AM(1),FIA,1,FIA,1))
         XLPNRM = DSQRT(DDOT(NGLMDT(1),XLAMDP,1,XLAMDP,1))
         XLHNRM = DSQRT(DDOT(NGLMDT(1),XLAMDH,1,XLAMDH,1))
         XINRM  = DSQRT(DDOT(NT1AM(ISYMPR)*NFREQ,XI,1,XI,1))
         PIJNRM = DSQRT(DDOT(NMATIJ(ISYMPR),PROPIJ,1,PROPIJ,1))
         PABNRM = DSQRT(DDOT(NMATAB(ISYMPR),PROPAB,1,PROPAB,1))
         CALL AROUND('Debug Info Exiting '//SECNAM)
         WRITE(LUPRI,'(A,I1)')
     &   'Symmetry of property integrals: ',ISYMPR
         WRITE(LUPRI,'(A,I10,/,A)')
     &   'Number of frequencies: ',NFREQ,' - values:'
         WRITE(LUPRI,'(1P,4D20.10)') (FREQ(IFREQ), IFREQ = 1,NFREQ)
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of F(ia)  : ',FIANRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDP : ',XLPNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XLAMDH : ',XLHNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of XI     : ',XINRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPIJ : ',PIJNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Norm of PROPAB : ',PABNRM
      ENDIF

      RETURN
      END
C  /* Deck cc_cyixi */
      SUBROUTINE CC_CYIXI(FOCKD,XI,FIA,WORK,LWORK,FREQ,ISYMXI,
     &                    XNRM,XI2NRM,DELP1,DELP2)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates for the effective right-hand
C              side for 1st order amplitude response.
C
#include <implicit.h>
      DIMENSION FOCKD(*), XI(*), FIA(*), WORK(LWORK)
      LOGICAL DELP1, DELP2
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>
#include <chocc2.h>

C     Amplitude info.
C     ---------------

      FAC1   = 1.0D0
      FAC2   = 0.0D0
      SCD    = 1.0D0
      SCDG   = 1.0D0

      JTYP1  = 6
      JTYP2  = 9
      ISYMP1 = 1
      ISYMP2 = ISYMXI
      NUMP12 = 1

      IOPTDN = 1
      IOPTCE = 1

C     Y info.
C     -------

      NUMZ  = 1
      JTYPZ = 1
      ISYMZ = 1
      JTYPY = 1

C     Calculate Y intermediates.
C     --------------------------

      CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NCHMOC,FAC1,NUMP12,
     &            DUMMY,IDUMMY,DUMMY,IDUMMY,FAC2,0,
     &            SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &            JTYPZ,ISYMZ,NUMCHO,NUMZ,JTYPY,
     &            FIA,1,1,XI,
     &            WORK,LWORK,XNRM,XI2NRM,DELP1,DELP2)

      RETURN
      END
C  /* Deck cc_choxi21 */
      SUBROUTINE CC_CHOXI21(T2AM,X2AM,PROPIJ,PROPAB,ISYMT,ISYMPR)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Transform T2AM. See CC_CHOXI1 and CC_CHOXI2 for details.
C
#include <implicit.h>
      DIMENSION T2AM(*), X2AM(*), PROPIJ(*), PROPAB(*)
#include <ccorb.h>
#include <ccsdsym.h>

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Get result symmetry.
C     --------------------

      ISYMX = MULD2H(ISYMT,ISYMPR)

C     Virtual part.
C     -------------

      DO ISYMBJ = 1,NSYM

         ISYMAI = MULD2H(ISYMBJ,ISYMX)
         ISYMCJ = MULD2H(ISYMAI,ISYMT)

         DO ISYMJ = 1,NSYM

            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            ISYMC = MULD2H(ISYMJ,ISYMCJ)

            DO J = 1,NRHF(ISYMJ)

               LBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
               LCJ = IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J - 1)

               KOFFT = IT2SQ(ISYMAI,ISYMCJ) + NT1AM(ISYMAI)*LCJ + 1
               KOFFP = IMATAB(ISYMB,ISYMC) + 1
               KOFFX = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*LBJ + 1

               NTOAI = MAX(NT1AM(ISYMAI),1)
               NTOTB = MAX(NVIR(ISYMB),1)

               CALL DGEMM('N','T',NT1AM(ISYMAI),NVIR(ISYMB),NVIR(ISYMC),
     &                    ONE,T2AM(KOFFT),NTOAI,PROPAB(KOFFP),NTOTB,
     &                    ZERO,X2AM(KOFFX),NTOAI)

            ENDDO

         ENDDO

      ENDDO

C     Occupied part.
C     --------------

      CALL DSCAL(NMATIJ(ISYMPR),XMONE,PROPIJ,1)
      DO ISYMBJ = 1,NSYM

         ISYMAI = MULD2H(ISYMBJ,ISYMX)
         ISYMBK = MULD2H(ISYMAI,ISYMT)

         DO ISYMJ = 1,NSYM

            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            ISYMK = MULD2H(ISYMB,ISYMBK)

            DO J = 1,NRHF(ISYMJ)

               LBJ   = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
               KOFFX = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*LBJ + 1

               DO K = 1,NRHF(ISYMK)

                  LBK = IT1AM(ISYMB,ISYMK)  + NVIR(ISYMB)*(K - 1)
                  KJ  = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(J - 1) + K

                  KOFFT = IT2SQ(ISYMAI,ISYMBK) + NT1AM(ISYMAI)*LBK + 1

                  LENAIB = NT1AM(ISYMAI)*NVIR(ISYMB)
                  FAC    = PROPIJ(KJ)

                  CALL DAXPY(LENAIB,FAC,T2AM(KOFFT),1,X2AM(KOFFX),1)

               ENDDO

            ENDDO

         ENDDO

      ENDDO
      CALL DSCAL(NMATIJ(ISYMPR),XMONE,PROPIJ,1)

      RETURN
      END
C  /* Deck cc_choxi22 */
      SUBROUTINE CC_CHOXI22(XLAMDP,XLAMDH,XI,X2AM,WORK,LWORK,ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the contributions
C
C        XI(ai) = XI(ai) + sum(ckb) X2AM(bi,ck) * (kc|ab)
C
C                        - sum(ckj) X2AM(aj,ck) * (kc|ji)
C
C     using Cholesky decomposed integrals. X2AM is assumed to be full-square.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), XI(*), X2AM(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      INTEGER IOFFL(8), IOFFZ(8), IOFFY(8)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOXI22')

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

C     Allocation for AO result.
C     -------------------------

      KXIAO = KEND0
      KEND1 = KXIAO + NT1AO(ISYMX)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: XIAO'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Initialize AO result.
C     ---------------------

      CALL DZERO(WORK(KXIAO),NT1AO(ISYMX))

C     Start Cholesky symmetry loop.
C     -----------------------------

      DO ISYCHO = 1,NSYM

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000

         ISYMBI = MULD2H(ISYCHO,ISYMX)
         ISYMAJ = ISYMBI
         ISYMDI = ISYMBI

C        Allocation for one Cholesky vector.
C        -----------------------------------

         LREAD = MEMRD(1,ISYCHO,IOPTR)

         KCHOL = KEND1
         KREAD = KCHOL + N2BST(ISYCHO)
         KEND2 = KREAD + LREAD
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: AO Chol.'
            WRITE(LUPRI,'(5X,A,I1)')
     &      'Cholesky symmetry: ',ISYCHO
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Set up batch.
C        SCRL: L(kc,#J), L(ji,#J), Y(b#J,i), L(g,d#J).
C        SCRZ: Y(bi,#J), Z(d#J,i).
C        ---------------------------------------------

         LENL = MAX(NT1AM(ISYCHO),NMATIJ(ISYCHO),NT1AM(ISYMBI),
     &              N2BST(ISYCHO))
         LENY = MAX(NT1AM(ISYMBI),NT1AO(ISYMDI))

         MINMEM = LENL + LENY
         IF (MINMEM .GT. 0) THEN
            NVEC = MIN(LWRK2/MINMEM,NUMCHO(ISYCHO))
         ELSE IF (MINMEM .EQ. 0) THEN
            GO TO 1000
         ELSE
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Error in ',SECNAM,':'
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'MINMEM is non-positive: ',MINMEM
            CALL QUIT('Error in '//SECNAM)
         ENDIF

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory for Cholesky batch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/,5X,A,I10,/)')
     &      'Minimum memory required   : ',MINMEM,
     &      'Memory available for batch: ',LWRK2,
     &      'Total memory available    : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

C        Open files.
C        -----------

         CALL CHO_MOP(-1,1,ISYCHO,LUCHKC,1,1)
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

C           Read L(kc,#J) into SCRL.
C           ------------------------

            CALL CHO_MOREAD(WORK(KSCRL),NT1AM(ISYCHO),NUMV,JVEC1,LUCHKC)

C           Calculate Y intermediates; store in SCRZ.
C           -----------------------------------------

            NTOBI = MAX(NT1AM(ISYMBI),1)
            NTOCK = MAX(NT1AM(ISYCHO),1)

            KOFFX = IT2SQ(ISYMBI,ISYCHO) + 1

            CALL DGEMM('N','N',NT1AM(ISYMBI),NUMV,NT1AM(ISYCHO),
     &                 ONE,X2AM(KOFFX),NTOBI,WORK(KSCRL),NTOCK,
     &                 ZERO,WORK(KSCRZ),NTOBI)

C           Read L(ji,#J) into SCRL.
C           ------------------------

            CALL CHO_MOREAD(WORK(KSCRL),NMATIJ(ISYCHO),NUMV,JVEC1,
     &                      LUCHJI)

C           Calculate occupied contribution.
C           --------------------------------

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
     &                       ONE,XI(KOFF3),NTOTA)

               ENDDO
            ENDDO

C           Resort Y intermediate as Y(b#J,i), store in SCRL.
C           -------------------------------------------------

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

C           Backtransform the Y intermediates to get Z(d#J,i) in SCRZ.
C           ----------------------------------------------------------

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

C           Read AO Cholesky vectors and store as L(g,d#J) in SCRL.
C           -------------------------------------------------------

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

C           Calculate virtual contribution in AO basis.
C           -------------------------------------------

            DO ISYMD = 1,NSYM

               ISYMG = MULD2H(ISYMD,ISYCHO)
               ISYMI = MULD2H(ISYMD,ISYMDI)

               KOFF1 = KSCRL + IOFFL(ISYMD)
               KOFF2 = KSCRZ + IOFFZ(ISYMD)
               KOFF3 = KXIAO + IT1AO(ISYMG,ISYMI)

               LENDJ = NBAS(ISYMD)*NUMV

               NTOTG = MAX(NBAS(ISYMG),1)
               NTODJ = MAX(LENDJ,1)

               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),LENDJ,
     &                    ONE,WORK(KOFF1),NTOTG,WORK(KOFF2),NTODJ,
     &                    ONE,WORK(KOFF3),NTOTG)

            ENDDO

         ENDDO

C        Close files.
C        ------------

         CALL CHO_MOP(1,4,ISYCHO,LUCHJI,1,1)
         CALL CHO_MOP(1,1,ISYCHO,LUCHKC,1,1)

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

C     Transform the virtual AO term to MO basis.
C     ------------------------------------------

      DO ISYMI = 1,NSYM

         ISYMA = MULD2H(ISYMI,ISYMX)
         ISYMG = ISYMA

         NTOTG = MAX(NBAS(ISYMG),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         KOFF1 = IGLMVI(ISYMG,ISYMA) + 1
         KOFF2 = KXIAO + IT1AO(ISYMG,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     &              ONE,XLAMDP(KOFF1),NTOTG,WORK(KOFF2),NTOTG,
     &              ONE,XI(KOFF3),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choxieta */
      SUBROUTINE CC_CHOXIETA(IXETRAN,NXETRAN,IOPTRES,IORDER,LISTL,
     &                       FILXI,IXDOTS,XCONS,FILETA,IEDOTS,ECONS,
     &                       MXVEC,WORK,LWORK)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Mimic the behaviour of CC_XIETA for Cholesky CC2.
C              Currently only used for Eta(A)*T(B) contributions to linear
C              response. See CC_XIETA for details (to the extent they are
C              applicable here).
C
#include <implicit.h>
#include <cclists.h>
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      INTEGER IXDOTS(MXVEC,NXETRAN), IEDOTS(MXVEC,NXETRAN)
      CHARACTER*(*) LISTL, FILXI, FILETA
      DIMENSION XCONS(MXVEC,NXETRAN), ECONS(MXVEC,NXETRAN)
      DIMENSION WORK(LWORK)
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOXIETA')

C     Return if nothing to do.
C     ------------------------

      IF (NXETRAN .LE. 0) RETURN

C     Check that the requested calculation has been implemented.
C     ----------------------------------------------------------

      JERR = 0

      IF (IOPTRES .NE. 5) THEN  ! Only dot products...
         WRITE(LUPRI,*) 'IOPTRES = ',IOPTRES,' not implemented!'
         JERR = JERR + 1
      ENDIF

      ICOUNT = 0
      DO ITRAN = 1,NXETRAN
         DO IVEC = 1,MXVEC
            ICOUNT = ICOUNT + IXDOTS(IVEC,ITRAN)
         ENDDO
      ENDDO
      IF (ICOUNT .NE. 0) THEN  ! No Xi vectors...
         WRITE(LUPRI,*) 'Xi vectors not implemented!'
         JERR = JERR + 1
      ENDIF

      IF (IORDER .NE. 1) THEN  ! Only 1st order implemented...
         WRITE(LUPRI,*) 'IORDER = ',IORDER,' not implemented!'
         JERR = JERR + 1
      ENDIF

      IF (LISTL(1:2) .NE. 'L0') THEN  ! Left vector must be 0th multipliers
         WRITE(LUPRI,*) 'LISTL = ',LISTL(1:3),' not implemented!'
         JERR = JERR + 1
      ENDIF

      IF (FILETA(1:2) .NE. 'R1') THEN  ! Right vector must be 1st am.
         WRITE(LUPRI,*) 'FILETA = ',FILETA(1:3),' not implemented!'
         JERR = JERR + 1
      ENDIF

      IF (JERR .NE. 0) THEN
         CALL FLSHFO(LUPRI)
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Do the calculation.
C     -------------------

      CALL CC_CHOETADOT(IXETRAN,NXETRAN,IEDOTS,ECONS,MXVEC,WORK,LWORK)

      RETURN
      END
C  /* Deck cc_choetadot */
      SUBROUTINE CC_CHOETADOT(IETRAN,NETRAN,IEDOTS,ECONS,MXDEN,
     &                        WORK,LWORK)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate EtaA*TB contributions to the linear response
C              function for Cholesky CC2. The contributions are added
C              into ECONS, mimicking the conventional code (see CC_XIETA).
C              The contributions are calculated from the 1st order right
C              densities, which must be available on disk, backtransformed
C              to the AO basis.
C
#include <implicit.h>
#include <cclists.h>
      INTEGER   IETRAN(MXDIM_XEVEC,NETRAN)
      INTEGER   IEDOTS(MXDEN,NETRAN)
      DIMENSION ECONS(MXDEN,NETRAN), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <ccroper.h>
#include <ccr1rsp.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOETADOT')

      LOGICAL LONLYAO
      DATA LONLYAO /.TRUE./

      CHARACTER*3 LISTD
      PARAMETER (LISTD = 'd01')

      CHARACTER*8 LABELA

      INTEGER ILSTSYM

C     Start loop through A operators (i.e. EtaA vectors).
C     ---------------------------------------------------

      DO IA = 1,NETRAN

         IOPERA = IETRAN(1,IA)
         LABELA = LBLOPR(IOPERA)
         ISYMA  = ISYOPR(IOPERA)
         IHERMA = ISYMAT(IOPERA)

C        Allocation.
C        -----------

         KOPERA = 1
         KDNBAO = KOPERA + N2BST(ISYMA)
         KEND1  = KDNBAO + N2BST(ISYMA)
         LWRK1  = LWORK  - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read integrals for A.
C        ---------------------

         IERR = -1
         CALL CCPRPAO(LABELA,LONLYAO,WORK(KOPERA),IRREP,IREAL,IERR,
     &                WORK(KEND1),LWRK1)

C        Check.
C        ------

         IF (IERR .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO'
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number     : ',IOPERA,' (in LBLOPR list)',
     &      'Operator label      : ',LABELA,
     &      'Operator symmetry   : ',ISYMA
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Operator Hermiticity: ',IHERMA
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERR
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERR .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': CCPRPAO failed to determine operator symmetry'
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8,/,5X,A,I10)')
     &      'Operator number     : ',IOPERA,' (in LBLOPR list)',
     &      'Operator label      : ',LABELA,
     &      'Operator symmetry   : ',ISYMA
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Operator Hermiticity: ',IHERMA
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERR
            CALL QUIT('Irrep. error in CCPRPAO ('//SECNAM//')')
         ENDIF

         IF (IRREP .NE. ISYMA) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Symmetry mismatch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &      'Operator number     : ',IOPERA,' (in LBLOPR list)',
     &      'Operator label      : ',LABELA
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Operator Hermiticity: ',IHERMA
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'CCPRPAO returned    : ',IRREP,
     &      '- expected          : ',ISYMA
            CALL QUIT('Symmetry mismatch in '//SECNAM)
         ENDIF

         IF (IREAL .NE. IHERMA) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Hermiticity mismatch in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,A,/,5X,A,2X,A8)')
     &      'Operator number  : ',IOPERA,' (in LBLOPR list)',
     &      'Operator label   : ',LABELA
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Operator symmetry: ',ISYMA
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'CCPRPAO returned : ',IREAL,
     &      '- expected       : ',IHERMA
casm        CALL QUIT('Hermiticity mismatch in '//SECNAM)
         ENDIF

C        Loop over densities (i.e. TB vectors).
C        --------------------------------------

         IDEN = 1
         DO WHILE ((IEDOTS(IDEN,IA).NE.0) .AND. IDEN.LE.MXDEN)

            ILSTDN = IEDOTS(IDEN,IA)
            ISYDEN = ILSTSYM(LISTD,ILSTDN)

            IF (ISYDEN .NE. ISYMA) THEN
               WRITE(LUPRI,'(//,5X,A,A,A)')
     &         'Symmetry mismatch in ',SECNAM,':'
               WRITE(LUPRI,'(5X,A,A)')
     &         'LISTD : ',LISTD
               WRITE(LUPRI,'(5X,A,I3)')
     &         'ISYDEN: ',ISYDEN
               WRITE(LUPRI,'(5X,A,I3,/)')
     &         'ISYMA : ',ISYMA
               CALL QUIT('Symmetry mismatch in '//SECNAM)
            ENDIF

            IOPT = 1
            CALL CC_RDRSPD(LISTD,ILSTDN,ISYDEN,IOPT,'CC2       ',
     &                     .FALSE.,WORK(KDNBAO))

            CON = DDOT(N2BST(ISYDEN),WORK(KOPERA),1,WORK(KDNBAO),1)
            ECONS(IDEN,IA) = ECONS(IDEN,IA) + CON

c-tbp:
c     write(lupri,*) 'ECONS(',IDEN,',',IA,') is:'
c     write(lupri,*) 'Eta(',LABELA,') * ',
c    &               'R1(',LRTLBL(IEDOTS(IDEN,IA)),',',
c    &                     FRQLRT(IEDOTS(IDEN,IA)),')'

            IDEN = IDEN + 1

         ENDDO

      ENDDO

      RETURN
      END
