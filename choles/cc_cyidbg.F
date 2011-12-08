C  /* Deck cc_cyidbg */
      SUBROUTINE CC_CYIDBG(WORK,LWORK)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Debug CC_CYI.
C
C     It is assumed that the CC2 energy code has completed at least
C     one iteration so that all MO Cholesky files are available.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <ccsdinp.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CYIDBG')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0)

C     Allocation.
C     -----------

      KFOCKD = 1
      KOMEG1 = KFOCKD + NORBTS
      KOMEGX = KOMEG1 + NT1AM(1)
      KFCKIA = KOMEGX + NT1AM(1)
      KEND1  = KFCKIA + NT1AM(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read canonical orbital energies.
C     --------------------------------

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &               .FALSE.,.TRUE.,.FALSE.)

C     Initialize.
C     -----------

      ECC2  = ZERO
      ECC2X = ZERO
      T2NRM = ZERO
      X2NRM = ZERO
      CALL DZERO(WORK(KOMEG1),NT1AM(1))
      CALL DZERO(WORK(KOMEGX),NT1AM(1))

      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(WORK(KFCKIA),NT1AM(1),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

C     Call debugged routine to get reference data.
C     --------------------------------------------

      CALL CC2_CHOYI1(WORK(KFOCKD),WORK(KOMEG1),WORK(KEND1),LWRK1,
     &                T2NRM,ECC2)

C     Make sure that target Y intermediates are removed from disk.
C     ------------------------------------------------------------

      DO ISYCHO = 1,NSYM
         CALL CC_CYIOP(-1,ISYCHO,0)
         CALL CC_CYIOP(0,ISYCHO,0)
      ENDDO

C     Call CC_CYI to do the same thing (more or less).
C     ------------------------------------------------

      LW1S  = LWRK1
      LWRK1 = MIN(LWRK1,MAX(NT2SQ(1)/10,5))
c     CALL CC_CYI(WORK(KFOCKD),WORK(KOMEGX),XDUM,WORK(KFCKIA),
c    &            WORK(KEND1),LWRK1,0,1,1,FDUM,X2NRM,X2CNM,.FALSE.)
      CALL CC_CYIV0(WORK(KFOCKD),WORK(KOMEGX),WORK(KFCKIA),
     &              WORK(KEND1),LWRK1,X2NRM,X2CNM,.FALSE.,.FALSE.)
      LWRK1 = LW1S

      CALL CC_CYIENR(WORK(KEND1),LWRK1,ECC2X)

C     Print result.
C     -------------

      CALL AROUND('Output from '//SECNAM)

      DIFFE = ECC2 - ECC2X
      WRITE(LUPRI,'(5X,A,F22.15,/,5X,A,F22.15,/,5X,A,1P,D22.15,/)')
     & 'Target CC2 energy: ',ECC2,
     & 'CC_CYI CC2 energy: ',ECC2X,
     & 'Difference       : ',DIFFE

      DIFFT = T2NRM - X2NRM
      WRITE(LUPRI,'(5X,A,F22.15,/,5X,A,F22.15,/,5X,A,1P,D22.15,/)')
     & 'Target T2 norm: ',T2NRM,
     & 'CC_CYI T2 norm: ',X2NRM,
     & 'Difference    : ',DIFFT

      ONRMT = DSQRT(DDOT(NT1AM(1),WORK(KOMEG1),1,WORK(KOMEG1),1))
      ONRMX = DSQRT(DDOT(NT1AM(1),WORK(KOMEGX),1,WORK(KOMEGX),1))
      DIFFN = ONRMT - ONRMX
      CALL DAXPY(NT1AM(1),XMONE,WORK(KOMEG1),1,WORK(KOMEGX),1)
      DNRM2 = DDOT(NT1AM(1),WORK(KOMEGX),1,WORK(KOMEGX),1)
      DRMS  = DSQRT(DNRM2/XT1AM(1))
      DNRM  = DSQRT(DNRM2)
      WRITE(LUPRI,'(5X,A,F22.15,/,5X,A,F22.15,/,5X,A,1P,D22.15)')
     & 'Target CC2 omega1 norm: ',ONRMT,
     & 'CC_CYI CC2 omega1 norm: ',ONRMX,
     & 'Difference of norms   : ',DIFFN
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15,/)')
     & 'Norm of difference    : ',DNRM,
     & 'RMS error             : ',DRMS

      RETURN
      END
C  /* Deck cc_choeimdbg */
      SUBROUTINE CC_CHOEIMDBG(WORK,LWORK)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Debug E intermediate routines.
C
C     Note: The Y intermediates and MO Cholesky vectors must be
C           available on disk.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOEIMDBG')

      LOGICAL EIJTHR, EABTHR

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, HALF = 0.5D0)

C     Make a pretty superficial test of the Y routines:
C     The Y-intermediate contributions to the CC2 energy
C     must be given by (tr[EIJ] - tr[EAB])/2.
C     --------------------------------------------------

      ETRG = ZERO
      CALL CC_CYIENR(WORK,LWORK,ETRG)

      KEIJ  = 1
      KEAB  = KEIJ  + NMATIJ(1)
      KEND1 = KEAB  + NMATAB(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL DZERO(WORK(KEIJ),NMATIJ(1))
      CALL DZERO(WORK(KEAB),NMATAB(1))

      ISYMY  = 1
      ISIDE  = 0
      ISYCHL = 1
      ITYP   = 1
      CALL CC_CHOEY(WORK(KEIJ),WORK(KEAB),WORK(KEND1),LWRK1,
     &              ISYMY,ISIDE,ISYCHL,ITYP,NUMCHO,.FALSE.)

      ECC2 = ZERO
      DO ISYM = 1,NSYM
         DO I = 1,NRHF(ISYM)
            KII  = KEIJ + IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1)
     &           + I - 1
            ECC2 = ECC2 + WORK(KII)
         ENDDO
         DO A = 1,NVIR(ISYM)
            KAA  = KEAB + IMATAB(ISYM,ISYM) + NVIR(ISYM)*(A - 1)
     &           + A - 1
            ECC2 = ECC2 - WORK(KAA)
         ENDDO
      ENDDO
      ECC2 = HALF*ECC2

      CALL HEADER('Comparing Y Contributions to CC2 Energy',-1)
      DIFF = ETRG - ECC2
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Target   CC2 energy contribution: ',ETRG
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'CC_CHOEY CC2 energy contribution: ',ECC2
      WRITE(LUPRI,'(5X,A,1P,D22.15,/)')
     & 'Difference                      : ',DIFF

C     If the E-intermediates are available on disk (from a
C     conventional run), read them and compare with Cholesky
C     results.
C     ------------------------------------------------------

      INQUIRE(FILE='CC_E2IM',EXIST=EIJTHR)
      INQUIRE(FILE='CC_E1IM',EXIST=EABTHR)

      IF (EIJTHR .AND. EABTHR) THEN

C        Allocation.
C        -----------

         KEIJTG = 1
         KEABTG = KEIJTG + NMATIJ(1)
         KEIJ   = KEABTG + NMATAB(1)
         KEAB   = KEIJ   + NMATIJ(1)
         KFOCKD = KEAB   + NMATAB(1)
         KT1AM  = KFOCKD + NORBTS
         KLAMDP = KT1AM  + NT1AM(1)
         KLAMDH = KLAMDP + NGLMDT(1)
         KEND1  = KLAMDH + NGLMDT(1)
         LWRK1  = LWORK  - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND1-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read target E-matrices.
C        -----------------------

         LUE = -1
         CALL GPOPEN(LUE,'CC_E2IM','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE)
         READ(LUE)(WORK(KEIJTG+I-1),I=1,NMATIJ(1))
         CALL GPCLOSE(LUE,'KEEP')

         LUE = -1
         CALL GPOPEN(LUE,'CC_E1IM','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE)
         READ(LUE)(WORK(KEABTG+I-1),I=1,NMATAB(1))
         CALL GPCLOSE(LUE,'KEEP')

C        Read orbital energies.
C        ----------------------

         CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &                  .FALSE.,.TRUE.,.FALSE.)

C        Read amplitudes.
C        ----------------

         IFAIL = -1
         CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,
     &                  IFAIL)

C        Calculate Lambda matrices.
C        --------------------------

         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &               WORK(KEND1),LWRK1)

C        Calculate Cholesky E-matrices.
C        ------------------------------

c         nwrk  = 3*n2bst(1) + 15*nt1amx
c         lw1s  = lwrk1
c         lwrk1 = min(lwrk1,nwrk)
c         write(lupri,*) '...note: calling CC_CHOEIM with LWRK1 = ',lwrk1

         CALL CC_CHOEIM0(WORK(KFOCKD),WORK(KT1AM),WORK(KLAMDP),
     &                   WORK(KLAMDH),WORK(KEIJ),WORK(KEAB),
     &                   WORK(KEND1),LWRK1,.FALSE.)

c         lwrk1 = lw1s

C        Compare.
C        --------

         CALL HEADER('Testing E-matrices',-1)

         TGNMIJ = DSQRT(DDOT(NMATIJ(1),WORK(KEIJTG),1,WORK(KEIJTG),1))
         TGNMAB = DSQRT(DDOT(NMATAB(1),WORK(KEABTG),1,WORK(KEABTG),1))
         CHNMIJ = DSQRT(DDOT(NMATIJ(1),WORK(KEIJ),1,WORK(KEIJ),1))
         CHNMAB = DSQRT(DDOT(NMATAB(1),WORK(KEAB),1,WORK(KEAB),1))

         CALL DAXPY(NMATIJ(1),XMONE,WORK(KEIJTG),1,WORK(KEIJ),1)
         CALL DAXPY(NMATAB(1),XMONE,WORK(KEABTG),1,WORK(KEAB),1)
         DIFNIJ = DDOT(NMATIJ(1),WORK(KEIJ),1,WORK(KEIJ),1)
         DIFRIJ = DSQRT(DIFNIJ/NMATIJ(1))
         DIFNIJ = DSQRT(DIFNIJ)
         DIFNAB = DDOT(NMATAB(1),WORK(KEAB),1,WORK(KEAB),1)
         DIFRAB = DSQRT(DIFNAB/NMATAB(1))
         DIFNAB = DSQRT(DIFNAB)

         DIFMIJ = -1.0D10
         DO ISYMJ = 1,NSYM
            ISYMI = ISYMJ
            DO J = 1,NRHF(ISYMJ)
               DO I = 1,NRHF(ISYMI)
                  KIJ = KEIJ + IMATIJ(ISYMI,ISYMJ)
     &                + NRHF(ISYMI)*(J - 1) + I - 1
                  DIFMIJ = MAX(DIFMIJ,DABS(WORK(KIJ)))
               ENDDO
            ENDDO
         ENDDO
         DIFMAB = -1.0D10
         DO ISYMB = 1,NSYM
            ISYMA = ISYMB
            DO B = 1,NVIR(ISYMB)
               DO A = 1,NVIR(ISYMA)
                  KAB = KEAB + IMATAB(ISYMA,ISYMB)
     &                + NVIR(ISYMA)*(B - 1) + A - 1
                  DIFMAB = MAX(DIFMAB,DABS(WORK(KAB)))
               ENDDO
            ENDDO
         ENDDO

         WRITE(LUPRI,'(5X,A,A,/)')
     &   'Note that the Cholesky E-intermediates are compared to the',
     &   ' *exact* ones!'

         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Target   EIJ norm : ',TGNMIJ
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Cholesky EIJ norm : ',CHNMIJ
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Difference        : ',TGNMIJ-CHNMIJ
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Norm of difference: ',DIFNIJ
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'RMS error         : ',DIFRIJ
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Max. abs. error   : ',DIFMIJ
         WRITE(LUPRI,'(5X,A,12X,I10,/)')
     &   'Number of elements: ',NMATIJ(1)

         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Target   EAB norm : ',TGNMAB
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Cholesky EAB norm : ',CHNMAB
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Difference        : ',TGNMAB-CHNMAB
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Norm of difference: ',DIFNAB
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'RMS error         : ',DIFRAB
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Max. abs. error   : ',DIFMAB
         WRITE(LUPRI,'(5X,A,12X,I10,/)')
     &   'Number of elements: ',NMATAB(1)

c         CALL HEADER('EIJ Difference Matrix',-1)
c         CALL NOCC_PRT(WORK(KEIJ),1,'IJ  ')
c         CALL HEADER('EAB Difference Matrix',-1)
c         CALL NOCC_PRT(WORK(KEAB),1,'AB  ')

      ELSE

         WRITE(LUPRI,'(/,5X,A,A,/)')
     &   SECNAM,': E-mat test not possible: conv. mat. not available!'

      ENDIF

      RETURN
      END
C  /* Deck cc_lamptdbg */
      SUBROUTINE CC_LAMPTDBG(WORK,LWORK)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Debug CC_LAMPT using CCLR_LAMTRA.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_LAMPTDBG')

C     Allocation.
C     -----------

      MXT1 = NT1AM(1)
      DO ISYM = 2,NSYM
         MXT1 = MAX(MXT1,NT1AM(ISYM))
      ENDDO

      MXLM = NGLMDT(1)
      DO ISYM = 2,NSYM
         MXLM = MAX(MXLM,NGLMDT(ISYM))
      ENDDO

      ISYLAM = 1

      KLAMDP = 1
      KLAMDH = KLAMDP + NGLMDT(ISYLAM)
      KX1AM  = KLAMDH + NGLMDT(ISYLAM)
      KLAMD2 = KX1AM  + MXT1
      KLAMPX = KLAMD2 + MXLM
      KLAMHX = KLAMPX + MXLM
      KEND1  = KLAMHX + MXLM
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate unperturbed Lambda matrices.
C     --------------------------------------

      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KX1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KX1AM),
     &            WORK(KEND1),LWRK1)

C     Tot. sym. check, ISIDE = 1.
C     ---------------------------

      ISIDE  = 1
      ISYMX  = 1
      CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPX),WORK(KLAMDH),
     &                 WORK(KLAMHX),WORK(KX1AM),ISYMX)
      CALL CC_LAMPT(WORK(KLAMDP),WORK(KLAMDH),WORK(KX1AM),WORK(KLAMD2),
     &              ISYLAM,ISYMX,ISIDE)
      CALL CC_LAMCHK(WORK(KLAMHX),WORK(KLAMD2),ISYMX,TRGNMO,TSTNMO,
     &               ERRMIO,ERRMXO,1)
      CALL CC_LAMCHK(WORK(KLAMPX),WORK(KLAMD2),ISYMX,TRGNMV,TSTNMV,
     &               ERRMIV,ERRMXV,2)
      CALL HEADER('Checking Perturbed Lambda Matrices',-1)
      WRITE(LUPRI,'(5X,A,I2,/,5X,A,I2)')
     & 'Symmetry of X1AM: ',ISYMX,
     & '           ISIDE: ',ISIDE
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Target norm, occupied part    : ',TRGNMO,
     & 'Test   norm, occupied part    : ',TSTNMO
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Min. abs. error, occupied part: ',ERRMIO,
     & 'Max. abs. error, occupied part: ',ERRMXO
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Target norm, virtual part     : ',TRGNMV,
     & 'Test   norm, virtual part     : ',TSTNMV
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Min. abs. error, virtual  part: ',ERRMIV,
     & 'Max. abs. error, virtual  part: ',ERRMXV

C     Tot. sym. check, ISIDE = -1.
C     ----------------------------

      ISIDE  = -1
      ISYMX  = 1
      CALL CCLR_LAMTRA(WORK(KLAMDH),WORK(KLAMHX),WORK(KLAMDP),
     &                 WORK(KLAMPX),WORK(KX1AM),ISYMX)
      CALL CC_LAMPT(WORK(KLAMDP),WORK(KLAMDH),WORK(KX1AM),WORK(KLAMD2),
     &              ISYLAM,ISYMX,ISIDE)
      CALL CC_LAMCHK(WORK(KLAMPX),WORK(KLAMD2),ISYMX,TRGNMO,TSTNMO,
     &               ERRMIO,ERRMXO,1)
      CALL CC_LAMCHK(WORK(KLAMHX),WORK(KLAMD2),ISYMX,TRGNMV,TSTNMV,
     &               ERRMIV,ERRMXV,2)
      CALL HEADER('Checking Perturbed Lambda Matrices',-1)
      WRITE(LUPRI,'(5X,A,I2,/,5X,A,I2)')
     & 'Symmetry of X1AM: ',ISYMX,
     & '           ISIDE: ',ISIDE
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Target norm, occupied part    : ',TRGNMO,
     & 'Test   norm, occupied part    : ',TSTNMO
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Min. abs. error, occupied part: ',ERRMIO,
     & 'Max. abs. error, occupied part: ',ERRMXO
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Target norm, virtual part     : ',TRGNMV,
     & 'Test   norm, virtual part     : ',TSTNMV
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Min. abs. error, virtual  part: ',ERRMIV,
     & 'Max. abs. error, virtual  part: ',ERRMXV

C     Loop over symmetries of X1AM.
C     -----------------------------

      DO ISYMX = 1,NSYM

C        Get a perturbed X1AM: read L(ia,1) of appropriate sym.
C        ------------------------------------------------------

         CALL CHO_MOP(-1,1,ISYMX,LUCHMO,1,1)
         CALL CHO_MOREAD(WORK(KX1AM),NT1AM(ISYMX),1,1,LUCHMO)
         CALL CHO_MOP(1,1,ISYMX,LUCHMO,1,1)

C        ISIDE = 1 test.
C        ---------------

         ISIDE = 1
         CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPX),WORK(KLAMDH),
     &                    WORK(KLAMHX),WORK(KX1AM),ISYMX)
         CALL CC_LAMPT(WORK(KLAMDP),WORK(KLAMDH),WORK(KX1AM),
     &                 WORK(KLAMD2),ISYLAM,ISYMX,ISIDE)
         CALL CC_LAMCHK(WORK(KLAMHX),WORK(KLAMD2),ISYMX,TRGNMO,TSTNMO,
     &                  ERRMIO,ERRMXO,1)
         CALL CC_LAMCHK(WORK(KLAMPX),WORK(KLAMD2),ISYMX,TRGNMV,TSTNMV,
     &                  ERRMIV,ERRMXV,2)
         CALL HEADER('Checking Perturbed Lambda Matrices',-1)
         WRITE(LUPRI,'(5X,A,I2,/,5X,A,I2)')
     &   'Symmetry of X1AM: ',ISYMX,
     &   '           ISIDE: ',ISIDE
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm, occupied part    : ',TRGNMO,
     &   'Test   norm, occupied part    : ',TSTNMO
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. error, occupied part: ',ERRMIO,
     &   'Max. abs. error, occupied part: ',ERRMXO
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm, virtual part     : ',TRGNMV,
     &   'Test   norm, virtual part     : ',TSTNMV
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. error, virtual  part: ',ERRMIV,
     &   'Max. abs. error, virtual  part: ',ERRMXV

C        ISIDE = -1 test.
C        ----------------

         ISIDE  = -1
         CALL CCLR_LAMTRA(WORK(KLAMDH),WORK(KLAMHX),WORK(KLAMDP),
     &                    WORK(KLAMPX),WORK(KX1AM),ISYMX)
         CALL CC_LAMPT(WORK(KLAMDP),WORK(KLAMDH),WORK(KX1AM),
     &                 WORK(KLAMD2),ISYLAM,ISYMX,ISIDE)
         CALL CC_LAMCHK(WORK(KLAMPX),WORK(KLAMD2),ISYMX,TRGNMO,TSTNMO,
     &                  ERRMIO,ERRMXO,1)
         CALL CC_LAMCHK(WORK(KLAMHX),WORK(KLAMD2),ISYMX,TRGNMV,TSTNMV,
     &                  ERRMIV,ERRMXV,2)
         CALL HEADER('Checking Perturbed Lambda Matrices',-1)
         WRITE(LUPRI,'(5X,A,I2,/,5X,A,I2)')
     &   'Symmetry of X1AM: ',ISYMX,
     &   '           ISIDE: ',ISIDE
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm, occupied part    : ',TRGNMO,
     &   'Test   norm, occupied part    : ',TSTNMO
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. error, occupied part: ',ERRMIO,
     &   'Max. abs. error, occupied part: ',ERRMXO
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm, virtual part     : ',TRGNMV,
     &   'Test   norm, virtual part     : ',TSTNMV
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. error, virtual  part: ',ERRMIV,
     &   'Max. abs. error, virtual  part: ',ERRMXV

      ENDDO

      RETURN
      END
C  /* Deck cc_lamchk */
      SUBROUTINE CC_LAMCHK(TARGET,TEST,ISYM,TRGNRM,TSTNRM,
     &                     ERRMIN,ERRMAX,IOPT)
C
#include <implicit.h>
      DIMENSION TARGET(*), TEST(*)
#include <ccorb.h>
#include <ccsdsym.h>

      ERRMIN = 1.0D10
      ERRMAX = -1.0D10
      TRGNRM = 0.0D0
      TSTNRM = 0.0D0

      IF (IOPT .EQ. 1) THEN

C        Occupied Lambda check.
C        ----------------------

         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMI,ISYM)
            DO I = 1,NRHF(ISYMI)
               DO A = 1,NBAS(ISYMA)
                  LAI = IGLMRH(ISYMA,ISYMI) + NBAS(ISYMA)*(I - 1) + A
                  ERR = DABS(TARGET(LAI)-TEST(LAI))
                  ERRMIN = MIN(ERRMIN,ERR)
                  ERRMAX = MAX(ERRMAX,ERR)
                  TRGNRM = TRGNRM + TARGET(LAI)*TARGET(LAI)
                  TSTNRM = TSTNRM + TEST(LAI)*TEST(LAI)
               ENDDO
            ENDDO
         ENDDO

      ELSE

C        Virtual Lambda check.
C        ---------------------

         DO ISYMA = 1,NSYM
            ISYMG = MULD2H(ISYMA,ISYM)
            DO A = 1,NVIR(ISYMA)
               DO G = 1,NBAS(ISYMG)
                  LGA = IGLMVI(ISYMG,ISYMA) + NBAS(ISYMG)*(A - 1) + G
                  ERR = DABS(TARGET(LGA)-TEST(LGA))
                  ERRMIN = MIN(ERRMIN,ERR)
                  ERRMAX = MAX(ERRMAX,ERR)
                  TRGNRM = TRGNRM + TARGET(LGA)*TARGET(LGA)
                  TSTNRM = TSTNRM + TEST(LGA)*TEST(LGA)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      TRGNRM = DSQRT(TRGNRM)
      TSTNRM = DSQRT(TSTNRM)

      RETURN
      END
C  /* Deck cc_chotgdbg */
      SUBROUTINE CC_CHOTGDBG(WORK,LWORK)
C
C     Thomas Bondo Pedersen, January 2003.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTGDBG')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Check perturbed Lambda matrices.
C     --------------------------------

      CALL CC_LAMPTDBG(WORK,LWORK)

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

C     Get unperturbed Lambda matrices.
C     --------------------------------

      ISYLAM = 1

      KLAMDP = KEND0
      KLAMDH = KLAMDP + NGLMDT(ISYLAM)
      KT1AM  = KLAMDH + NGLMDT(ISYLAM)
      KEND1  = KT1AM  + NT1AM(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

C     Reset memory.
C     -------------

      KEND1 = KT1AM
      LWRK1 = LWORK - KEND1 + 1

C     Loop over perturbation symmetries.
C     ----------------------------------

      DO ISYMX = 1,NSYM

C        Pick a vector to check.
C        -----------------------

         IF (ISYMX .LE. 2) THEN
            ISYCHO = ISYMX + 1
            ISYCHO = MIN(ISYCHO,NSYM)
         ELSE
            ISYCHO = ISYMX - 1
         ENDIF
         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000

         JVEC = MAX(NUMCHO(ISYCHO)/5,1)

C        Calculate the Cholesky vector L(ab,JVEC).
C        -----------------------------------------

         KLAB  = KEND1
         KCHOL = KLAB  + NMATAB(ISYCHO)
         KEND2 = KCHOL + N2BST(ISYCHO)
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                  ISYCHO,2,WORK(KEND2),LWRK2)

         DO ISYMB = 1,NSYM

            ISYMD = ISYMB
            ISYMG = MULD2H(ISYMD,ISYCHO)
            ISYMA = ISYMG

            NEED = NBAS(ISYMG)*NVIR(ISYMB)
            IF (NEED .GT. LWRK2) THEN
               WRITE(LUPRI,'(//,5X,A,A)')
     &         'Insufficient memory in ',SECNAM
               WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &         'Need (more than): ',KEND2+NEED-1,
     &         'Available       : ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            KOFFC = KCHOL  + IAODIS(ISYMG,ISYMD)
            KOFFH = KLAMDH + IGLMVI(ISYMD,ISYMB)
            KOFFP = KLAMDP + IGLMVI(ISYMG,ISYMA)
            KOFFR = KLAB   + IMATAB(ISYMA,ISYMB)

            NTOTG = MAX(NBAS(ISYMG),1)
            NTOTD = MAX(NBAS(ISYMD),1)
            NTOTA = MAX(NVIR(ISYMA),1)

            CALL DGEMM('N','N',NBAS(ISYMG),NVIR(ISYMB),NBAS(ISYMD),
     &                 ONE,WORK(KOFFC),NTOTG,WORK(KOFFH),NTOTD,
     &                 ZERO,WORK(KEND2),NTOTG)
            CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMG),
     &                 ONE,WORK(KOFFP),NTOTG,WORK(KEND2),NTOTG,
     &                 ZERO,WORK(KOFFR),NTOTA)

         ENDDO

C        Allocation.
C        -----------

         ISYMAI = MULD2H(ISYCHO,ISYMX)

         KX1AM  = KCHOL
         KLIJ   = KX1AM  + NT1AM(ISYMX)
         KTARGT = KLIJ   + NMATIJ(ISYCHO)
         KLAI   = KTARGT + NT1AM(ISYMAI)
         KEND3  = KLAI   + NT1AM(ISYMAI)
         LWRK3  = LWORK  - KEND3 + 1

         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND3-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read the first vector L(ia,1) that fits symmetry
C        into X1AM.
C        ------------------------------------------------

         CALL CHO_MOP(-1,1,ISYMX,LUCHMO,1,1)
         CALL CHO_MOREAD(WORK(KX1AM),NT1AM(ISYMX),1,1,LUCHMO)
         CALL CHO_MOP(1,1,ISYMX,LUCHMO,1,1)

C        Read the L(ij,JVEC).
C        --------------------

         CALL CHO_MOP(-1,4,ISYCHO,LUCHMO,1,1)
         CALL CHO_MOREAD(WORK(KLIJ),NMATIJ(ISYCHO),1,JVEC,LUCHMO)
         CALL CHO_MOP(1,4,ISYCHO,LUCHMO,1,1)

C        Test for ISIDE = -1.
C        ====================

c      nwrk  = 2*n2bst(isymx) + 25*nt1amx
c      lsav  = lwrk3
c      lwrk3 = min(lwrk3,nwrk)
c      write(lupri,*) '   ...calling cc_chotg with lwrk3 = ',lwrk3

         ISIDE = -1
         CALL CC_CHOTG(WORK(KLAMDP),WORK(KLAMDH),WORK(KX1AM),
     &                 WORK(KEND3),LWRK3,1,ISYMX,ISIDE)

c      lwrk3 = lsav

         CALL CHO_MOP(-1,8,ISYCHO,LUCHMO,1,ISYMX)
         CALL CHO_MOREAD(WORK(KLAI),NT1AM(ISYMAI),1,JVEC,LUCHMO)
         CALL CHO_MOP(1,8,ISYCHO,LUCHMO,1,ISYMX)

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYMX)
            ISYMA = MULD2H(ISYMB,ISYCHO)

            KOFF1 = KLAB   + IMATAB(ISYMB,ISYMA)
            KOFF2 = KX1AM  + IT1AM(ISYMB,ISYMI)
            KOFF3 = KTARGT + IT1AM(ISYMA,ISYMI)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTB = MAX(NVIR(ISYMB),1)

            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                 ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTB,
     &                 ZERO,WORK(KOFF3),NTOTA)

         ENDDO

         DO ISYMJ = 1,NSYM

            ISYMI = MULD2H(ISYMJ,ISYCHO)
            ISYMA = MULD2H(ISYMJ,ISYMX)

            KOFF1 = KX1AM  + IT1AM(ISYMA,ISYMJ)
            KOFF2 = KLIJ   + IMATIJ(ISYMI,ISYMJ)
            KOFF3 = KTARGT + IT1AM(ISYMA,ISYMI)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTI = MAX(NRHF(ISYMI),1)

            CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                 XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTI,
     &                 ONE,WORK(KOFF3),NTOTA)

         ENDDO

         CALL DBGDIFAI(WORK(KTARGT),WORK(KLAI),TRGNRM,TSTNRM,
     &                 ERRMIN,ERRMAX,ISYMAI)
         CALL HEADER('Test of CC_CHOTG',-1)
         WRITE(LUPRI,'(5X,A,I2,/,5X,A,I2,/,5X,A,I2,/,5X,A,I2)')
     &   'ISYMX  = ',ISYMX,
     &   'ISYCHO = ',ISYCHO,
     &   'ISYMAI = ',ISYMAI,
     &   'ISIDE  = ',ISIDE
         WRITE(LUPRI,'(5X,A,I10)')
     &   'Vector checked: ',JVEC
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm   : ',TRGNRM,
     &   'Test   norm   : ',TSTNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. err.: ',ERRMIN,
     &   'Max. abs. err.: ',ERRMAX

c         write(lupri,*) ' - difference matrix:'
c         call daxpy(nt1am(isymai),-1.0D0,work(ktargt),1,work(klai),1)
c         call nocc_prt(work(klai),isymai,'AI  ')

C        Test for ISIDE = +1.
C        ====================

c      nwrk  = 2*n2bst(isymx) + 25*nt1amx
c      lsav  = lwrk3
c      lwrk3 = min(lwrk3,nwrk)
c      write(lupri,*) '   ...calling cc_chotg with lwrk3 = ',lwrk3

         ISIDE = 1
         CALL CC_CHOTG(WORK(KLAMDP),WORK(KLAMDH),WORK(KX1AM),
     &                 WORK(KEND3),LWRK3,1,ISYMX,ISIDE)

c      lwrk3 = lsav
         
         CALL CHO_MOP(-1,9,ISYCHO,LUCHMO,1,ISYMX)
         CALL CHO_MOREAD(WORK(KLAI),NT1AM(ISYMAI),1,JVEC,LUCHMO)
         CALL CHO_MOP(1,9,ISYCHO,LUCHMO,1,ISYMX)

         DO ISYMI = 1,NSYM

            ISYMB = MULD2H(ISYMI,ISYMX)
            ISYMA = MULD2H(ISYMB,ISYCHO)
            
            KOFF1 = KLAB   + IMATAB(ISYMA,ISYMB)
            KOFF2 = KX1AM  + IT1AM(ISYMB,ISYMI)
            KOFF3 = KTARGT + IT1AM(ISYMA,ISYMI)
         
            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTB = MAX(NVIR(ISYMB),1) 

            CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                 ONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTB,
     &                 ZERO,WORK(KOFF3),NTOTA)

         ENDDO

         DO ISYMJ = 1,NSYM

            ISYMI = MULD2H(ISYMJ,ISYCHO)
            ISYMA = MULD2H(ISYMJ,ISYMX)
         
            KOFF1 = KX1AM  + IT1AM(ISYMA,ISYMJ)
            KOFF2 = KLIJ   + IMATIJ(ISYMJ,ISYMI)
            KOFF3 = KTARGT + IT1AM(ISYMA,ISYMI)

            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTJ = MAX(NRHF(ISYMJ),1)

            CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                 XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTJ,
     &                 ONE,WORK(KOFF3),NTOTA)

         ENDDO

         CALL DBGDIFAI(WORK(KTARGT),WORK(KLAI),TRGNRM,TSTNRM,
     &                 ERRMIN,ERRMAX,ISYMAI)
         CALL HEADER('Test of CC_CHOTG',-1)
         WRITE(LUPRI,'(5X,A,I2,/,5X,A,I2,/,5X,A,I2,/,5X,A,I2)')
     &   'ISYMX  = ',ISYMX,
     &   'ISYCHO = ',ISYCHO,
     &   'ISYMAI = ',ISYMAI,
     &   'ISIDE  = ',ISIDE
         WRITE(LUPRI,'(5X,A,I10)')
     &   'Vector checked: ',JVEC
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm   : ',TRGNRM,
     &   'Test   norm   : ',TSTNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. err.: ',ERRMIN,
     &   'Max. abs. err.: ',ERRMAX

 1000    CONTINUE

      ENDDO

C     Clean up.
C     ---------

      DO ISYCHO = 1,NSYM
         CALL CHO_MOP(-1,8,ISYCHO,LUCHMO,1,NSYM)
         CALL CHO_MOP(0,8,ISYCHO,LUCHMO,1,NSYM)
         CALL CHO_MOP(-1,9,ISYCHO,LUCHMO,1,NSYM)
         CALL CHO_MOP(0,9,ISYCHO,LUCHMO,1,NSYM)
      ENDDO

      RETURN
      END
C  /* Deck dbgdifai */
      SUBROUTINE DBGDIFAI(TARGET,TEST,TRGNRM,TSTNRM,ERRMIN,ERRMAX,
     &                    ISYMAI)
C
#include <implicit.h>
      DIMENSION TARGET(*), TEST(*)
#include <ccorb.h>
#include <ccsdsym.h>

      ERRMIN = 1.0D10
      ERRMAX = -1.0D10
      TRGNRM = DSQRT(DDOT(NT1AM(ISYMAI),TARGET,1,TARGET,1))
      TSTNRM = DSQRT(DDOT(NT1AM(ISYMAI),TEST,1,TEST,1))

      DO ISYMI = 1,NSYM
         ISYMA = MULD2H(ISYMI,ISYMAI)
         DO I = 1,NRHF(ISYMI)
            DO A = 1,NVIR(ISYMA)
               LAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
               ERR = DABS(TARGET(LAI) - TEST(LAI))
               ERRMIN = MIN(ERRMIN,ERR)
               ERRMAX = MAX(ERRMAX,ERR)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_chocimdbg */
      SUBROUTINE CC_CHOCIMDBG(WORK,LWORK)
C
C     Thomas Bondo Pedersen, January 2003.
C
C     Purpose: Debug CC_CHOCIM.
C
C     It is assumed that CC_CIA is bug-free!
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dccsdsym.h>
#include <priunit.h>
#include <ciarc.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOCIMDBG')

      INTEGER LSYM(8)

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      DO ISYM = 1,NSYM
         LSYM(ISYM) = ISYM
      ENDDO

C     Check memory.
C     -------------

      XWORK = ONE*LWORK
      IF (XT2SQ(1) .GT. LWORK) THEN
         WRITE(LUPRI,'(//,5X,A,A,//)')
     &   'Oooops. Not enough memory for test in ',SECNAM
         RETURN
      ENDIF

C     Allocation.
C     -----------

      KFOCKD = 1
      KT2AM  = KFOCKD + NORBTS
      KEND1  = KT2AM  + NT2SQ(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      write(lupri,*) 'lwork : ',lwork
      write(lupri,*) 'kfockd: ',kfockd
      write(lupri,*) 'kt2am : ',kt2am
      write(lupri,*) 'kend1 : ',kend1
      write(lupri,*) 'lwrk1 : ',lwrk1

C     Get Fock matrix diagonal.
C     -------------------------

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &               .FALSE.,.TRUE.,.FALSE.)

c      write(lupri,'(/,A)')
c     & 'Canonical orbital energies after read:'
c      do isym = 1,nsym
c         write(lupri,'(/,A,I1,A)') 'Occupied orbital energies, sym. ',
c     &                             isym,':'
c         write(lupri,'(A,/)')      '----------------------------------'
c         if (nrhf(isym) .gt. 0) then
c            do i = 1,nrhf(isym)
c               koffi = kfockd + irhf(isym) + i - 1
c               write(lupri,'(F15.8)') work(koffi)
c            enddo
c         else
c            write(lupri,'(A)') '- empty symmetry block.'
c         endif
c         write(lupri,'(/,A,I1,A)') 'Virtual orbital energies, sym. ',
c     &                             isym,':'
c         write(lupri,'(A,/)')      '---------------------------------'
c         if (nvir(isym) .gt. 0) then
c            do a = 1,nvir(isym)
c               koffa = kfockd + ivir(isym) + a - 1
c               write(lupri,'(F15.8)') work(koffa)
c            enddo
c         else
c            write(lupri,'(A)') '- empty symmetry block.'
c         endif
c         write(lupri,*)
c      enddo

C     Calculate 2CME amplitudes.
C     --------------------------

      IPRCIA = IPRLVL + 1

      LIDEN  = .TRUE.
      SYMTRZ = .FALSE.
      CIAMIO = .FALSE.
      INXINT = .TRUE.
      GETMNM = .FALSE.

      ISYCH1 = 1
      ISYCH2 = 1
      ITYP1  = 3
      ITYP2  = 3

      DO ISYM = 1,NSYM
         IOFA1(ISYM)  = 1
         LVIRA(ISYM)  = NVIR(ISYM)
         NX1AMA(ISYM) = NT1AM(ISYM)
         IOFB1(ISYM)  = 1
         LVIRB(ISYM)  = NVIR(ISYM)
         NX1AMB(ISYM) = NT1AM(ISYM)
         DO JSYM = 1,NSYM
            IX1AMA(JSYM,ISYM) = IT1AM(JSYM,ISYM)
            IX1AMB(JSYM,ISYM) = IT1AM(JSYM,ISYM)
            IX2SQ(JSYM,ISYM)  = IT2SQ(JSYM,ISYM)
         ENDDO
         NTOVEC(ISYM) = NUMCHO(ISYM)
      ENDDO
      NX2SQ = NT2SQ(1)

      SCD   = ONE
      NCALL = 0
      FREQ  = ZERO
      T2NRM = ZERO
      CALL CC_CIA(WORK(KT2AM),WORK(KEND1),LWRK1,SCD,KLAST,NCALL)
      CALL CC_DNOM(WORK(KFOCKD),WORK(KT2AM),FREQ,1)
      CALL CC_CYINRM(WORK(KT2AM),1,T2NRM)
      CALL CC_CYITCME(WORK(KT2AM),WORK(KEND1),LWRK1,1,KLAST)

      WRITE(LUPRI,'(A,1P,D22.15)')
     & 'Packed T2 norm in cc_chocimdbg: ',DSQRT(T2NRM)

C     Loop over perturbation symmetries.
C     ----------------------------------

      DO ISYMX = 1,NSYM

         IF (NUMCHO(ISYMX) .LE. 0) GO TO 1000
         IF (NT1AM(ISYMX)  .LE. 0) GO TO 1000

C        Pick a number of vectors.
C        -------------------------

         NUMX = MIN(NUMCHO(ISYMX),7)

C        Allocation.
C        -----------

         KCTARG = KEND1
         KCTEST = KCTARG + NT1AM(ISYMX)*NUMX
         KX1AM  = KCTEST + NT1AM(ISYMX)*NUMX
         KEND2  = KX1AM  + NT1AM(ISYMX)*NUMX
         LWRK2  = LWORK  - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Number of vectors: ',NUMX
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read perturbation vectors (i.e. Cholesky vectors).
C        --------------------------------------------------

         CALL CHO_MOP(-1,1,ISYMX,LUCHOX,1,1)
         CALL CHO_MOREAD(WORK(KX1AM),NT1AM(ISYMX),NUMX,1,LUCHOX)
         CALL CHO_MOP(1,1,ISYMX,LUCHOX,1,1)

C        Calculate target C intermediates.
C        ---------------------------------

         NTAI  = MAX(NT1AM(ISYMX),1)
         NTBJ  = NTAI
         KOFFT = KT2AM + IT2SQ(ISYMX,ISYMX)

         CALL DGEMM('N','N',NT1AM(ISYMX),NUMX,NT1AM(ISYMX),
     &              ONE,WORK(KOFFT),NTAI,WORK(KX1AM),NTBJ,
     &              ZERO,WORK(KCTARG),NTAI)

C        Calculate C intermediates using CC_CHOCIM.
C        ------------------------------------------

c         lwsav = lwrk2
c         nwrk  = (NUMX+6)*NT1AM(ISYMX) + (NUMX+8)*NMATIJ(1)
c         lwrk2 = min(lwrk2,nwrk)
c         write(lupri,*) '         calling cc_chocim with lwrk2=',lwrk2

         CALL CC_CHOCIM(WORK(KFOCKD),WORK(KX1AM),WORK(KEND2),LWRK2,
     &                  ISYMX,NUMX)

c         lwrk2 = lwsav

C        Read test results.
C        ------------------

         CALL CHO_IMOP(-1,3,LUCHOC,ISYMX)
         CALL CHO_MOREAD(WORK(KCTEST),NT1AM(ISYMX),NUMX,1,LUCHOC)
         CALL CHO_IMOP(0,3,LUCHOC,ISYMX)

C        Compare.
C        --------

         DO IX = 1,NUMX

            KOFF1 = KCTARG + NT1AM(ISYMX)*(IX - 1)
            KOFF2 = KCTEST + NT1AM(ISYMX)*(IX - 1)
            CALL DBGDIFAI(WORK(KOFF1),WORK(KOFF2),TRGNRM,TSTNRM,
     &                    ERRMIN,ERRMAX,ISYMX)

            CALL HEADER('Testing C Intermediates',-1)
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I3,A,I3)')
     &      'Symmetry of trial vectors: ',ISYMX,
     &      'Trial vector',IX,' of',NUMX
            WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &      'Target norm   : ',TRGNRM,
     &      'Test   norm   : ',TSTNRM
            WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &      'Min. abs. err.: ',ERRMIN,
     &      'Max. abs. err.: ',ERRMAX

         ENDDO

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck cc_choltrgijhd */
      SUBROUTINE CC_CHOLTRGIJHD(XLAMDP,XLAMDH,SIGMA1,X1AM,CIM,WORK,
     &                          LWORK,ISYMX,CALY,CALX,CALC)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the GIJ and H terms for left-hand Jacobian
C              transformations. DUMMY version.
C
C     Formula:
C
C     SIGMA1(ai) = SIGMA1(ai)
C                + sum(Jb) L(ba,J) * Y(bi,J)
C                - sum(Jj) Y(aj,J) * L(ij,J)
C                + 2 * sum(J) L(ia,J) sum(bj) L(jb,J) * C(bj)
C                -    sum(Jj) L(ja,J) sum(b)  L(ib,J) * C(bj)
C                + 2 * sum(J) L(ia,J) sum(bj) L(bj,J) * X(bj)
C                -    sum(Jj) L(ij,J) sum(b)  L(ba,J) * X(bj)
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), SIGMA1(*), X1AM(*), CIM(*)
      DIMENSION WORK(LWORK)
      LOGICAL CALY, CALX, CALC
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER* 14 SECNAM
      PARAMETER (SECNAM = 'CC_CHOLTRGIJHD')

      PARAMETER (IOPTR = 2)
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      snrm = dsqrt(ddot(nt1am(isymx),sigma1,1,sigma1,1))
      xnrm = dsqrt(ddot(nt1am(isymx),x1am,1,x1am,1))
      cnrm = dsqrt(ddot(nt1am(isymx),cim,1,cim,1))
      write(lupri,'(A,A,1P,D22.15)')
     & SECNAM,': norm of sigma1: ',snrm
      write(lupri,'(A,A,1P,D22.15)')
     & SECNAM,': norm of x1am  : ',xnrm
      write(lupri,'(A,A,1P,D22.15)')
     & SECNAM,': norm of cim   : ',cnrm
      write(lupri,*) SECNAM,': CALY, CALX, CALC: ',CALY,CALX,CALC
      write(lupri,*) 'lwork: ',lwork

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

      DO ISYCHO = 1,NSYM

         ISYMBI = MULD2H(ISYCHO,ISYMX)
         ISYMAJ = ISYMBI

C        Allocation.
C        -----------

         KCHOL = KEND0
         KCHIJ = KCHOL + MAX(NT1AM(ISYMBI),N2BST(ISYCHO))
         KCHAI = KCHIJ + NMATIJ(ISYCHO)
         KCHIA = KCHAI + NT1AM(ISYCHO)
         KCHAB = KCHIA + NT1AM(ISYCHO)
         KEND1 = KCHAB + NMATAB(ISYCHO)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insuf. mem. in '//SECNAM)
         ENDIF

         KYIM = KCHOL

C        Open files.
C        -----------

         CALL CC_CYIOP(-1,ISYCHO,-1)
         CALL CHO_MOP(-1,1,ISYCHO,LUCHIA,1,1)
         CALL CHO_MOP(-1,3,ISYCHO,LUCHAI,1,1)
         CALL CHO_MOP(-1,4,ISYCHO,LUCHIJ,1,1)

         DO JVEC = 1,NUMCHO(ISYCHO)

C           Read MO vectors.
C           ----------------

            CALL CHO_MOREAD(WORK(KCHIA),NT1AM(ISYCHO),1,JVEC,LUCHIA)
            CALL CHO_MOREAD(WORK(KCHAI),NT1AM(ISYCHO),1,JVEC,LUCHAI)
            CALL CHO_MOREAD(WORK(KCHIJ),NMATIJ(ISYCHO),1,JVEC,LUCHIJ)

C           Read AO vector.
C           ---------------

            CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                     ISYCHO,IOPTR,WORK(KEND1),LWRK1)

C           Transform to L(ab,J).
C           ---------------------

            DO ISYMB = 1,NSYM

               ISYMA = MULD2H(ISYMB,ISYCHO)
               ISYMD = ISYMB
               ISYMG = ISYMA

               NEED = NBAS(ISYMG)*NVIR(ISYMB)

               IF (NEED .GT. LWRK1) CALL QUIT(' WORK 1')

               NTOTG = MAX(NBAS(ISYMG),1)
               NTOTD = MAX(NBAS(ISYMD),1)

               KOFFC = KCHOL + IAODIS(ISYMG,ISYMD)
               KOFFH = IGLMVI(ISYMD,ISYMB) + 1

               CALL DGEMM('N','N',
     &                    NBAS(ISYMG),NVIR(ISYMB),NBAS(ISYMD),
     &                    ONE,WORK(KOFFC),NTOTG,XLAMDH(KOFFH),NTOTD,
     &                    ZERO,WORK(KEND1),NTOTG)

               NTOTA = MAX(NVIR(ISYMA),1)

               KOFFP = IGLMVI(ISYMG,ISYMA) + 1
               KOFFC = KCHAB + IMATAB(ISYMA,ISYMB)

               CALL DGEMM('T','N',
     &                    NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMG),
     &                    ONE,XLAMDP(KOFFP),NTOTG,WORK(KEND1),NTOTG,
     &                    ZERO,WORK(KOFFC),NTOTA)

            ENDDO

C           Contributions from X.
C           =====================

            IF (CALX) THEN

               IF (ISYCHO .EQ. ISYMX) THEN
                  FAC = TWO*DDOT(NT1AM(ISYMX),WORK(KCHAI),1,X1AM,1)
                  CALL DAXPY(NT1AM(ISYMX),FAC,WORK(KCHIA),1,SIGMA1,1)
               ENDIF

               DO ISYMJ = 1,NSYM

                  ISYMB = MULD2H(ISYMJ,ISYMX)
                  ISYMA = MULD2H(ISYMB,ISYCHO)
                  ISYMI = MULD2H(ISYMJ,ISYCHO)

                  NEED = NVIR(ISYMA)*NRHF(ISYMJ)

                  IF (NEED .GT. LWRK1) CALL QUIT(' WORK 2')

                  KOFF1 = KCHAB + IMATAB(ISYMB,ISYMA)
                  KOFF2 = IT1AM(ISYMB,ISYMJ) + 1

                  NTOTB = MAX(NVIR(ISYMB),1)
                  NTOTA = MAX(NVIR(ISYMA),1)

                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMJ),NVIR(ISYMB),
     &                       ONE,WORK(KOFF1),NTOTB,X1AM(KOFF2),NTOTB,
     &                       ZERO,WORK(KEND1),NTOTA)

                  KOFF1 = KCHIJ + IMATIJ(ISYMI,ISYMJ)
                  KOFF2 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KEND1),NTOTA,WORK(KOFF1),NTOTI,
     &                       ONE,SIGMA1(KOFF2),NTOTA)

               ENDDO

            ENDIF

C           Contributions from Y.
C           =====================

            IF (CALY) THEN

C              Read Y intermediate.
C              --------------------

               CALL CC_CYIRDF(WORK(KYIM),1,JVEC,ISYCHO,ISYMX,-1)

C              Calculate Y contributions.
C              --------------------------

               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYMBI)
                  ISYMA = MULD2H(ISYMB,ISYCHO)
   
                  KOFF1 = KCHAB + IMATAB(ISYMB,ISYMA)
                  KOFF2 = KYIM  + IT1AM(ISYMB,ISYMI)
                  KOFF3 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTB = MAX(NVIR(ISYMB),1)

                  CALL DGEMM('T','N',
     &                       NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &                       ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTB,
     &                       ONE,SIGMA1(KOFF3),NTOTA)

               ENDDO

               DO ISYMJ = 1,NSYM

                  ISYMA = MULD2H(ISYMJ,ISYMAJ)
                  ISYMI = MULD2H(ISYMJ,ISYCHO)

                  KOFF1 = KYIM  + IT1AM(ISYMA,ISYMJ)
                  KOFF2 = KCHIJ + IMATIJ(ISYMI,ISYMJ)
                  KOFF3 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTI,
     &                       ONE,SIGMA1(KOFF3),NTOTA)

               ENDDO

            ENDIF

C           Contributions from the C intermediate.
C           ======================================

            IF (CALC) THEN

               IF (ISYCHO .EQ. ISYMX) THEN
                  FAC = TWO*DDOT(NT1AM(ISYMX),WORK(KCHIA),1,CIM,1)
                  CALL DAXPY(NT1AM(ISYMX),FAC,WORK(KCHIA),1,SIGMA1,1)
               ENDIF

               DO ISYMI = 1,NSYM

                  ISYMB = MULD2H(ISYMI,ISYCHO)
                  ISYMJ = MULD2H(ISYMB,ISYMX)
                  ISYMA = MULD2H(ISYMJ,ISYCHO)

                  NEED = NRHF(ISYMI)*NRHF(ISYMJ)

                  IF (NEED .GT. LWRK1) CALL QUIT(' WORK 3')

                  KOFF1 = KCHIA + IT1AM(ISYMB,ISYMI)
                  KOFF2 = IT1AM(ISYMB,ISYMJ) + 1

                  NTOTB = MAX(NVIR(ISYMB),1)
                  NTOTI = MAX(NRHF(ISYMI),1)

                  CALL DGEMM('T','N',
     &                       NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),
     &                       ONE,WORK(KOFF1),NTOTB,CIM(KOFF2),NTOTB,
     &                       ZERO,WORK(KEND1),NTOTI)

                  KOFF1 = KCHIA + IT1AM(ISYMA,ISYMJ)
                  KOFF2 = IT1AM(ISYMA,ISYMI) + 1

                  NTOTA = MAX(NVIR(ISYMA),1)

                  CALL DGEMM('N','T',
     &                       NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &                       XMONE,WORK(KOFF1),NTOTA,WORK(KEND1),NTOTI,
     &                       ONE,SIGMA1(KOFF2),NTOTA)

               ENDDO

            ENDIF

         ENDDO

C        Close files.
C        ------------

         CALL CC_CYIOP(1,ISYCHO,-1)
         CALL CHO_MOP(1,4,ISYCHO,LUCHIJ,1,1)
         CALL CHO_MOP(1,3,ISYCHO,LUCHAI,1,1)
         CALL CHO_MOP(1,1,ISYCHO,LUCHIA,1,1)

      ENDDO

      RETURN
      END
C  /* Deck cc_choatrldbg */
      SUBROUTINE CC_CHOATRLHDBG(SIGMA1,X1AM,WORK,LWORK,ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Debug routines called by CC_CHOATR for ISIDE = -1.
C              CC_CYI, however, is assumed correct.....
C
#include <implicit.h>
      DIMENSION SIGMA1(*), X1AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <priunit.h>

      LOGICAL FALS, TRU
      PARAMETER (FALS = .FALSE., TRU = .TRUE.)

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CC_CHOATRLHDBG')

      PARAMETER (TINY = 1.0D-14)
      PARAMETER (ZERO = 0.0D0)

C     Allocation.
C     -----------

      KSIGM1 = 1
      KSIGM2 = KSIGM1 + NT1AM(ISYMX)
      KLAMDP = KSIGM2 + NT1AM(ISYMX)
      KLAMDH = KLAMDP + NGLMDT(1)
      KT1AM  = KLAMDH + NGLMDT(1)
      KFOCKD = KT1AM  + NT1AM(1)
      KX1AM  = KFOCKD + NORBTS
      KCIM1  = KX1AM  + NT1AM(ISYMX)
      KCIM2  = KCIM1  + NT1AM(ISYMX)
      KEND1  = KCIM2  + NT1AM(ISYMX)
      LWRK1  = LWORK  - KEND1 + 1

      KFIA = KT1AM

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insuf. mem. in '//SECNAM//' (1)')
      ENDIF

C     Read Fock diagonal, T1AM, and calculate Lambda matrices.
C     Replace T1AM with F(ia).
C     --------------------------------------------------------

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &               FALS,TRU,FALS)
      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,FALS,TRU,FALS,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     &            LWRK1)
      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(WORK(KFIA),NT1AM(1),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

c      lwsav = lwrk1
c      nwrk  = nmatij(1) + nmatab(1) + 25*nt1am(1) + n2bst(1)
c      lwrk1 = min(lwrk1,nwrk)
c      mxdsv = mxdecc
c      ncrsv = nchrdc
c      mxdecc = 5
c      nchrdc = 10
c      write(lupri,*) SECNAM,': reducing work space to ',lwrk1

C     Make sure the C intermediates are available.
C     --------------------------------------------

      CALL CC_CHOCIM(WORK(KFOCKD),X1AM,WORK(KEND1),LWRK1,ISYMX,1)

      CALL CHO_IMOP(-1,3,LUCIM,ISYMX)
      CALL CHO_MOREAD(WORK(KCIM1),NT1AM(ISYMX),1,1,LUCIM)
      CALL CHO_IMOP(1,3,LUCIM,ISYMX)

C     Make sure the Y intermediates are available.
C     --------------------------------------------

      FREQ = ZERO
      CALL DZERO(SIGMA1,NT1AM(ISYMX))
      CALL CC_CHOTG(WORK(KLAMDP),WORK(KLAMDH),X1AM,WORK(KEND1),LWRK1,1,
     &              ISYMX,-1)
c     CALL CC_CYI(WORK(KFOCKD),SIGMA1,X1AM,WORK(KFIA),
c    &            WORK(KEND1),LWRK1,-1,ISYMX,1,FREQ,X2NRM,
c    &            X2CNM,.TRUE.)
      CALL CC_CYILTR(WORK(KFOCKD),X1AM,WORK(KFIA),WORK(KEND1),LWRK1,
     &               ISYMX,FREQ,X2NRM,X2CNM,.FALSE.,.TRUE.)

      SNRM = DSQRT(DDOT(NT1AM(ISYMX),SIGMA1,1,SIGMA1,1))
      IF (SNRM .GE. TINY) THEN
         WRITE(LUPRI,'(A,A,1P,D22.15)')
     &   SECNAM,': error: sigma-norm after CC_CYI: ',SNRM
         WRITE(LUPRI,'(A)') ' - should be identically 0...'
         CALL QUIT('Error detected in '//SECNAM)
      ENDIF

C     Test CC_CHOATR0.
C     ----------------

      CALL CC_ECNLDBG(WORK(KSIGM1),X1AM,WORK(KEND1),LWRK1,ISYMX)
      CALL CC_CHOATR0(WORK(KSIGM2),X1AM,WORK(KEND1),LWRK1,ISYMX,1,-1)

      CALL DBGDIFAI(WORK(KSIGM1),WORK(KSIGM2),SI1NRM,SI2NRM,ERRMIN,
     &              ERRMAX,ISYMX)
      CALL HEADER('Testing CC_CHOATR0 (E contr.)',-1)
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Target  norm  : ',SI1NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'CHOATR0 norm  : ',SI2NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Min. abs. err.: ',ERRMIN
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Max. abs. err.: ',ERRMAX
      IF (ERRMAX .GT. TINY) THEN
         WRITE(LUPRI,'(5X,A,/)')
     &   '- error!'
         CALL QUIT('Error detected in '//SECNAM)
      ELSE
         WRITE(LUPRI,'(5X,A,/)')
     &   '- E contributions seem OK!'
      ENDIF

C     Test CC_CHOLTRGIJH: Y intermediate contributions only.
C     ------------------------------------------------------

      CALL DZERO(WORK(KSIGM1),NT1AM(ISYMX))
      CALL DZERO(WORK(KSIGM2),NT1AM(ISYMX))
      CALL DZERO(WORK(KCIM2),NT1AM(ISYMX))
      CALL DZERO(WORK(KX1AM),NT1AM(ISYMX))

      CALL CC_CHOLTRGIJHD(WORK(KLAMDP),WORK(KLAMDH),WORK(KSIGM1),
     &                    WORK(KX1AM),WORK(KCIM2),WORK(KEND1),LWRK1,
     &                    ISYMX,TRU,FALS,FALS)
      CALL CC_CHOLTRGIJH(WORK(KLAMDP),WORK(KLAMDH),WORK(KSIGM2),
     &                   WORK(KX1AM),WORK(KCIM2),WORK(KEND1),LWRK1,
     &                   ISYMX)

      CALL DBGDIFAI(WORK(KSIGM1),WORK(KSIGM2),SI1NRM,SI2NRM,ERRMIN,
     &              ERRMAX,ISYMX)
      CALL HEADER('Testing CC_CHOLTRGIJH: Y only',-1)
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Target  norm  : ',SI1NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'LTRGIJH norm  : ',SI2NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Min. abs. err.: ',ERRMIN
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Max. abs. err.: ',ERRMAX
      IF (ERRMAX .GT. TINY) THEN
         WRITE(LUPRI,'(5X,A,/)')
     &   '- error!'
         CALL QUIT('Error detected in '//SECNAM)
      ELSE
         WRITE(LUPRI,'(5X,A,/)')
     &   '- Y seems OK!'
      ENDIF

C     Test CC_CHOLTRGIJH: Y intermediate and X1AM contributions only.
C     ---------------------------------------------------------------

      CALL DZERO(WORK(KSIGM1),NT1AM(ISYMX))
      CALL DZERO(WORK(KSIGM2),NT1AM(ISYMX))
      CALL DZERO(WORK(KCIM2),NT1AM(ISYMX))
      CALL DCOPY(NT1AM(ISYMX),X1AM,1,WORK(KX1AM),1)

      CALL CC_CHOLTRGIJHD(WORK(KLAMDP),WORK(KLAMDH),WORK(KSIGM1),
     &                    WORK(KX1AM),WORK(KCIM2),WORK(KEND1),LWRK1,
     &                    ISYMX,TRU,TRU,FALS)
      CALL CC_CHOLTRGIJH(WORK(KLAMDP),WORK(KLAMDH),WORK(KSIGM2),
     &                   WORK(KX1AM),WORK(KCIM2),WORK(KEND1),LWRK1,
     &                   ISYMX)

      CALL DBGDIFAI(WORK(KSIGM1),WORK(KSIGM2),SI1NRM,SI2NRM,ERRMIN,
     &              ERRMAX,ISYMX)
      CALL HEADER('Testing CC_CHOLTRGIJH: Y and X only',-1)
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Target  norm  : ',SI1NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'LTRGIJH norm  : ',SI2NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Min. abs. err.: ',ERRMIN
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Max. abs. err.: ',ERRMAX
      IF (ERRMAX .GT. TINY) THEN
         WRITE(LUPRI,'(5X,A,/)')
     &   '- error!'
         CALL QUIT('Error detected in '//SECNAM)
      ELSE
         WRITE(LUPRI,'(5X,A,/)')
     &   '- Y and X seem OK!'
      ENDIF

C     Test CC_CHOLTRGIJH: All contributions.
C     --------------------------------------

      CALL DZERO(WORK(KSIGM1),NT1AM(ISYMX))
      CALL DZERO(WORK(KSIGM2),NT1AM(ISYMX))
      CALL DCOPY(NT1AM(ISYMX),WORK(KCIM1),1,WORK(KCIM2),1)
      CALL DCOPY(NT1AM(ISYMX),X1AM,1,WORK(KX1AM),1)

      CALL CC_CHOLTRGIJHD(WORK(KLAMDP),WORK(KLAMDH),WORK(KSIGM1),
     &                    WORK(KX1AM),WORK(KCIM2),WORK(KEND1),LWRK1,
     &                    ISYMX,TRU,TRU,TRU)
      CALL CC_CHOLTRGIJH(WORK(KLAMDP),WORK(KLAMDH),WORK(KSIGM2),
     &                   WORK(KX1AM),WORK(KCIM2),WORK(KEND1),LWRK1,
     &                   ISYMX)
      
      CALL DBGDIFAI(WORK(KSIGM1),WORK(KSIGM2),SI1NRM,SI2NRM,ERRMIN,
     &              ERRMAX,ISYMX)
      CALL HEADER('Testing CC_CHOLTRGIJH: all contributions',-1)
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Target  norm  : ',SI1NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'LTRGIJH norm  : ',SI2NRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Min. abs. err.: ',ERRMIN
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Max. abs. err.: ',ERRMAX
      IF (ERRMAX .GT. TINY) THEN
         WRITE(LUPRI,'(5X,A,/)')
     &   '- error!'
         CALL QUIT('Error detected in '//SECNAM)
      ELSE
         WRITE(LUPRI,'(5X,A,/)')
     &   '- All contributions seem OK!'
      ENDIF

c      lwrk1 = lwsav
c      mxdecc = mxdsv
c      nchrdc = ncrsv

      RETURN
      END
C  /* Deck cc_ecnldbg */
      SUBROUTINE CC_ECNLDBG(SIGMA1,X1AM,WORK,LWORK,ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Debug E contributions for LH trf.
C              SIGMA1 is initialized here...
C
#include <implicit.h>
      DIMENSION SIGMA1(*), X1AM(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_ECNLDBG')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Initialize.
C     -----------

      CALL DZERO(SIGMA1,NT1AM(ISYMX))

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

C     Calculate virtual contribution.
C     -------------------------------

      DO ISYMI = 1,NSYM

         ISYMB = MULD2H(ISYMI,ISYMX)
         ISYMA = ISYMB

         KOFFE = KEAB + IMATAB(ISYMB,ISYMA)
         KOFFX = IT1AM(ISYMB,ISYMI) + 1
         KOFFS = IT1AM(ISYMA,ISYMI) + 1

         NTOTB = MAX(NVIR(ISYMB),1)
         NTOTA = MAX(NVIR(ISYMA),1)

         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &              ONE,WORK(KOFFE),NTOTB,X1AM(KOFFX),NTOTB,
     &              ONE,SIGMA1(KOFFS),NTOTA)

      ENDDO

C     Calculate occupied contribution.
C     --------------------------------

      DO ISYMJ = 1,NSYM

         ISYMA = MULD2H(ISYMJ,ISYMX)
         ISYMI = ISYMJ

         KOFFX = IT1AM(ISYMA,ISYMJ) + 1
         KOFFE = KEIJ + IMATIJ(ISYMI,ISYMJ)
         KOFFS = IT1AM(ISYMA,ISYMI) + 1

         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTI = MAX(NRHF(ISYMI),1)

         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              XMONE,X1AM(KOFFX),NTOTA,WORK(KOFFE),NTOTI,
     &              ONE,SIGMA1(KOFFS),NTOTA)

      ENDDO

      RETURN
      END
C  /* Deck cc_choxidbg */
      SUBROUTINE CC_CHOXIDBG(WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Test right-hand side for 1st order amplitude response
C              equations.
C
C     N.B. most of this is internal consistency checking, not really
C          debugging....
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <chocc2.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CC_CHOXIDBG')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

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

      KFOCKD = KEND0
      KFIA   = KFOCKD + NORBTS
      KLAMDP = KFIA   + NT1AM(1)
      KLAMDH = KLAMDP + NGLMDT(1)
      KT1AM  = KLAMDH + NGLMDT(1)
      KEND1  = KT1AM  + NT1AM(1)
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read orbital energies.
C     ----------------------

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &               .FALSE.,.TRUE.,.FALSE.)

C     Read F(ia).
C     -----------

      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(WORK(KFIA),NT1AM(1),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

C     Get Lambda matrices.
C     --------------------

      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     &            LWRK1)

C     Reset memory.
C     -------------

      KEND1 = KT1AM
      LWRK1 = LWORK - KEND1 + 1

C     Decompose amplitudes.
C     ---------------------

      IF (.NOT. CHOT2C) THEN
         CALL CC_CHOCIM1(WORK(KFOCKD),DUM1,DUM2,WORK(KEND1),LWRK1,0,0)
      ENDIF

C     Loop over Cholesky symmetries.
C     ------------------------------

      DO ISYCHO = 1,NSYM

C        Pick a vector.
C        --------------

         IF (NUMCHO(ISYCHO) .GT. 10) THEN
            JVEC = MAX(NUMCHO(ISYCHO)/10,5)
         ELSE IF (NUMCHO(ISYCHO) .GT. 0) THEN
            JVEC = NUMCHO(ISYCHO)
         ELSE
            GO TO 1000
         ENDIF

C        Allocation.
C        -----------

         KCHOL = KEND1
         KCHIJ = KCHOL + N2BST(ISYCHO)
         KCHAI = KCHIJ + NMATIJ(ISYCHO)
         KCHIA = KCHAI + NT1AM(ISYCHO)
         KCHAB = KCHIA + NT1AM(ISYCHO)
         KTRIA = KCHAB + NMATAB(ISYCHO)
         KTRAI = KTRIA + NT1AM(ISYCHO)
         KTRIJ = KTRAI + NT1AM(ISYCHO)
         KXI1  = KTRIJ + NMATIJ(ISYCHO)
         KXI2  = KXI1  + NT1AM(ISYCHO)
         KEND2 = KXI2  + NT1AM(ISYCHO)
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,' - allocation: prop'
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need (more than): ',KEND2-1,
     &      'Available       : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Read vector.
C        ------------

         CALL CHO_READN(WORK(KCHOL),JVEC,1,WORK(KIND1),IDUM2,
     &                  ISYCHO,2,WORK(KEND2),LWRK2)

C        Transform.
C        ----------

c        CALL CC_CHOTRFMO(WORK(KLAMDP),WORK(KLAMDH),WORK(KCHOL),
c    &                    WORK(KCHIJ),WORK(KCHAI),
c    &                    WORK(KCHIA),WORK(KCHAB),
c    &                    WORK(KEND2),LWRK2,1,1,ISYCHO)

         CALL CC_CHOTRFED(WORK(KLAMDP),WORK(KLAMDH),WORK(KCHOL),
     &                    WORK(KCHAI),WORK(KCHIA),
     &                    WORK(KEND2),LWRK2,1,1,ISYCHO)

         CALL CC_CHOTRFOV(WORK(KLAMDP),WORK(KLAMDH),WORK(KCHOL),
     &                    WORK(KCHIJ),WORK(KCHAB),
     &                    WORK(KEND2),LWRK2,1,1,ISYCHO)

C        Read target vectors.
C        --------------------

         CALL CHO_MOP(-1,4,ISYCHO,LUCHIJ,1,1)
         CALL CHO_MOREAD(WORK(KTRIJ),NMATIJ(ISYCHO),1,JVEC,LUCHIJ)
         CALL CHO_MOP(1,4,ISYCHO,LUCHIJ,1,1)

         CALL CHO_MOP(-1,3,ISYCHO,LUCHAI,1,1)
         CALL CHO_MOREAD(WORK(KTRAI),NT1AM(ISYCHO),1,JVEC,LUCHAI)
         CALL CHO_MOP(1,3,ISYCHO,LUCHAI,1,1)

         CALL CHO_MOP(-1,1,ISYCHO,LUCHIA,1,1)
         CALL CHO_MOREAD(WORK(KTRIA),NT1AM(ISYCHO),1,JVEC,LUCHIA)
         CALL CHO_MOP(1,1,ISYCHO,LUCHIA,1,1)

C        Test (no test of virtual part!).
C        --------------------------------

         CALL DBGDIFIJ(WORK(KTRIJ),WORK(KCHIJ),TIJNRM,CIJNRM,ERMNIJ,
     &                 ERMXIJ,ISYCHO)
         CALL DBGDIFAI(WORK(KTRAI),WORK(KCHAI),TAINRM,CAINRM,ERMNAI,
     &                 ERMXAI,ISYCHO)
         CALL DBGDIFAI(WORK(KTRIA),WORK(KCHIA),TIANRM,CIANRM,ERMNIA,
     &                 ERMXIA,ISYCHO)

         CALL HEADER('Testing trf. IJ',-1)
         WRITE(LUPRI,'(A,I1,/,A,I10)')
     &   'Cholesky symmetry: ',ISYCHO,
     &   'Vector checked   : ',JVEC
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Target norm    : ',TIJNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'TRF    norm    : ',CIJNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Difference     : ',TIJNRM-CIJNRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Min. abs. err. : ',ERMNIJ
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Max. abs. err. : ',ERMXIJ

         CALL HEADER('Testing trf. AI',-1)
         WRITE(LUPRI,'(A,I1,/,A,I10)')
     &   'Cholesky symmetry: ',ISYCHO,
     &   'Vector checked   : ',JVEC
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Target norm    : ',TAINRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'TRF    norm    : ',CAINRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Difference     : ',TAINRM-CAINRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Min. abs. err. : ',ERMNAI
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Max. abs. err. : ',ERMXAI

         CALL HEADER('Testing trf. IA',-1)
         WRITE(LUPRI,'(A,I1,/,A,I10)')
     &   'Cholesky symmetry: ',ISYCHO,
     &   'Vector checked   : ',JVEC
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Target norm    : ',TIANRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'TRF    norm    : ',CIANRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Difference     : ',TIANRM-CIANRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Min. abs. err. : ',ERMNIA
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Max. abs. err. : ',ERMXIA

C        Calculate XI.
C        -------------

         FREQ = ZERO
         CALL DZERO(WORK(KXI1),NT1AM(ISYCHO))
         CALL DZERO(WORK(KXI2),NT1AM(ISYCHO))
         CALL CC_CHOXI2(WORK(KFOCKD),WORK(KFIA),WORK(KLAMDP),
     &                  WORK(KLAMDH),
     &                  WORK(KXI1),WORK(KCHIJ),WORK(KCHAB),FREQ,
     &                  WORK(KEND2),LWRK2,1,ISYCHO)
         CALL CC_CHOXI3(WORK(KFOCKD),WORK(KFIA),WORK(KLAMDP),
     &                  WORK(KLAMDH),
     &                  WORK(KXI2),WORK(KCHIJ),WORK(KCHAB),FREQ,
     &                  WORK(KEND2),LWRK2,1,ISYCHO)

         CALL DBGDIFAI(WORK(KXI1),WORK(KXI2),XI1NRM,XI2NRM,ERRMIN,
     &                 ERRMAX,ISYCHO)

         CALL HEADER('Testing XI',-1)
         WRITE(LUPRI,'(A,I1,/,A,I10)')
     &   'Cholesky symmetry: ',ISYCHO,
     &   'Vector used      : ',JVEC
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'XI1 norm (xi2) : ',XI1NRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'XI2 norm (xi3) : ',XI2NRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Difference     : ',XI1NRM-XI2NRM
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Min. abs. dif. : ',ERRMIN
         WRITE(LUPRI,'(A,1P,D22.15)')
     &   'Max. abs. dif. : ',ERRMAX

         IF (ISYCHO .EQ. 1) THEN

            CALL DZERO(WORK(KCHIJ),NMATIJ(1))
            DO ISYMI = 1,NSYM
               DO I = 1,NRHF(ISYMI)
                  II = KCHIJ + IMATIJ(ISYMI,ISYMI) + NRHF(ISYMI)*(I - 1)
     &               + I - 1
                  WORK(II) = ONE
               ENDDO
            ENDDO

            CALL DZERO(WORK(KCHAB),NMATAB(1))
            DO ISYMA = 1,NSYM
               DO A = 1,NVIR(ISYMA)
                  LA = KCHAB + IMATAB(ISYMA,ISYMA) + NVIR(ISYMA)*(A - 1)
     &               + A - 1
                  WORK(LA) = ONE
               ENDDO
            ENDDO

            FREQ = ZERO
            CALL DZERO(WORK(KXI1),NT1AM(ISYCHO))
            CALL DZERO(WORK(KXI2),NT1AM(ISYCHO))
            CALL CC_CHOXI2(WORK(KFOCKD),WORK(KFIA),WORK(KLAMDP),
     &                     WORK(KLAMDH),
     &                     WORK(KXI1),WORK(KCHIJ),WORK(KCHAB),FREQ,
     &                     WORK(KEND2),LWRK2,1,ISYCHO)
            CALL CC_CHOXI3(WORK(KFOCKD),WORK(KFIA),WORK(KLAMDP),
     &                     WORK(KLAMDH),
     &                     WORK(KXI2),WORK(KCHIJ),WORK(KCHAB),FREQ,
     &                     WORK(KEND2),LWRK2,1,ISYCHO)

            CALL DBGDIFAI(WORK(KXI1),WORK(KXI2),XI1NRM,XI2NRM,ERRMIN,
     &                    ERRMAX,ISYCHO)

            CALL HEADER('Testing XI, unit props',-1)
            WRITE(LUPRI,'(A,1P,D22.15)')
     &      'XI1 norm (xi2) : ',XI1NRM
            WRITE(LUPRI,'(A,1P,D22.15)')
     &      'XI2 norm (xi3) : ',XI2NRM
            WRITE(LUPRI,'(A,1P,D22.15)')
     &      'Difference     : ',XI1NRM-XI2NRM
            WRITE(LUPRI,'(A,1P,D22.15)')
     &      'Min. abs. dif. : ',ERRMIN
            WRITE(LUPRI,'(A,1P,D22.15)')
     &      'Min. abs. dif. : ',ERRMAX

         ENDIF

C        Escape point for empty symmetry.
C        --------------------------------

 1000    CONTINUE

      ENDDO

      RETURN
      END
C  /* Deck dbgdifij */
      SUBROUTINE DBGDIFIJ(TARGET,TEST,TRGNRM,TSTNRM,ERRMIN,ERRMAX,
     &                    ISYMIJ)
C
#include <implicit.h>
      DIMENSION TARGET(*), TEST(*)
#include <ccorb.h>
#include <ccsdsym.h>

      ERRMIN = 1.0D10
      ERRMAX = -1.0D10
      TRGNRM = DSQRT(DDOT(NMATIJ(ISYMIJ),TARGET,1,TARGET,1))
      TSTNRM = DSQRT(DDOT(NMATIJ(ISYMIJ),TEST,1,TEST,1))

      DO ISYMJ = 1,NSYM
         ISYMI = MULD2H(ISYMJ,ISYMIJ)
         DO J = 1,NRHF(ISYMJ)
            DO I = 1,NRHF(ISYMI)
               IJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
               ERR = DABS(TARGET(IJ) - TEST(IJ))
               ERRMIN = MIN(ERRMIN,ERR)
               ERRMAX = MAX(ERRMAX,ERR)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_chotg2dbg */
      SUBROUTINE CC_CHOTG2DBG(XLAMDP,XLAMDH,XLAMD2,X1AM,WORK,LWORK,
     &                        ISYMX)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: debug CC_CHOTG2.
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), XLAMD2(*), X1AM(*), WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOTG2DBG')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER (ITYPIA = 1, ITYPAI = 9, ITYPIJ = 10, ITYSCR = 11)

C     Move vectors to be tested to dummy file.
C     ========================================

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYMX)

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1000
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 1000

         MINMEM = NT1AM(ISYMAI)
         NVEC   = MIN(LWORK/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

         CALL CHO_MOP(-1,ITYPAI,ISYCHO,LUCHAI,1,ISYMX)
         CALL CHO_MOP(-1,ITYSCR,ISYCHO,LUSCR,1,ISYMX)

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF
            JVEC1 = NVEC*(IBATCH - 1) + 1

            CALL CHO_MOREAD(WORK,NT1AM(ISYMAI),NUMV,JVEC1,LUCHAI)
            CALL CHO_MOWRITE(WORK,NT1AM(ISYMAI),NUMV,JVEC1,LUSCR)

         ENDDO

         CALL CHO_MOP(0,ITYPAI,ISYCHO,LUCHAI,1,ISYMX)
         CALL CHO_MOP(1,ITYSCR,ISYCHO,LUSCR,1,ISYMX)

 1000    CONTINUE

      ENDDO

C     Now transform to obtain Lbar(ai,J) using chotg [which should work].
C     ===================================================================

      CALL CC_CHOTG(XLAMDP,XLAMDH,X1AM,WORK,LWORK,1,ISYMX,1)

C     Check Lbar(ai,J).
C     =================

      DO ISYCHO = 1,NSYM

         ISYMAI = MULD2H(ISYCHO,ISYMX)

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1001
         IF (NT1AM(ISYMAI)  .LE. 0) GO TO 1001

C        Allocation.
C        -----------

         KTARG = 1
         KTEST = KTARG + NT1AM(ISYMAI)
         KEND1 = KTEST + NT1AM(ISYMAI)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Pick a vector.
C        --------------

         JVEC = MIN(NUMCHO(ISYCHO),12)

         CALL CHO_MOP(-1,ITYPAI,ISYCHO,LUCHAI,1,ISYMX)
         CALL CHO_MOP(-1,ITYSCR,ISYCHO,LUSCR,1,ISYMX)
         CALL CHO_MOREAD(WORK(KTARG),NT1AM(ISYMAI),1,JVEC,LUCHAI)
         CALL CHO_MOREAD(WORK(KTEST),NT1AM(ISYMAI),1,JVEC,LUSCR)
         CALL CHO_MOP(1,ITYPAI,ISYCHO,LUCHAI,1,ISYMX)
         CALL CHO_MOP(0,ITYSCR,ISYCHO,LUSCR,1,ISYMX)

         CALL DBGDIFAI(WORK(KTARG),WORK(KTEST),TRGNRM,TSTNRM,
     &                 ERRMIN,ERRMAX,ISYMAI)

         CALL HEADER('Testing Lbar(ai,J)',-1)
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I3,A,I3)')
     &   'Symmetry of trial vector: ',ISYMX,
     &   'Picked Cholesky vector number',JVEC,' of sym.',ISYCHO
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm   : ',TRGNRM,
     &   'Test   norm   : ',TSTNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. err.: ',ERRMIN,
     &   'Max. abs. err.: ',ERRMAX

 1001    CONTINUE

      ENDDO

C     Check L(ij,J).
C     ==============

      DO ISYCHO = 1,NSYM

         ISYMIJ = MULD2H(ISYCHO,ISYMX)

         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 1002
         IF (NMATIJ(ISYMIJ) .LE. 0) GO TO 1002
         IF (NT1AM(ISYCHO)  .LE. 0) GO TO 1002

C        Allocation.
C        -----------

         KTARG = 1
         KTEST = KTARG + NMATIJ(ISYMIJ)
         KCHIA = KTEST + NMATIJ(ISYMIJ)
         KEND1 = KCHIA + NT1AM(ISYCHO)
         LWRK1 = LWORK - KEND1 + 1
         
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

C        Pick a vector.
C        --------------

         JVEC = MIN(NUMCHO(ISYCHO),21)

C        Calculate target vector.
C        ------------------------

         CALL CHO_MOP(-1,ITYPIA,ISYCHO,LUCHIA,1,1)
         CALL CHO_MOREAD(WORK(KCHIA),NT1AM(ISYCHO),1,JVEC,LUCHIA)
         CALL CHO_MOP(1,ITYPIA,ISYCHO,LUCHIA,1,1)

         DO ISYMJ = 1,NSYM

            ISYMB = MULD2H(ISYMJ,ISYMX)
            ISYMI = MULD2H(ISYMJ,ISYMIJ)

            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTI = MAX(NRHF(ISYMI),1)

            KOFFC = KCHIA + IT1AM(ISYMB,ISYMI)
            KOFFX = IT1AM(ISYMB,ISYMJ) + 1
            KOFFT = KTARG + IMATIJ(ISYMI,ISYMJ)

            CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMB),
     &                 1.0D0,WORK(KOFFC),NTOTB,X1AM(KOFFX),NTOTB,
     &                 0.0D0,WORK(KOFFT),NTOTI)

         ENDDO

C        Read test vector.
C        -----------------

         CALL CHO_MOP(-1,ITYPIJ,ISYCHO,LUCHIJ,1,ISYMX)
         CALL CHO_MOREAD(WORK(KTEST),NMATIJ(ISYMIJ),1,JVEC,LUCHIJ)
         CALL CHO_MOP(1,ITYPIJ,ISYCHO,LUCHIJ,1,ISYMX)

C        Compare.
C        --------

         CALL DBGDIFIJ(WORK(KTARG),WORK(KTEST),TRGNRM,TSTNRM,
     &                 ERRMIN,ERRMAX,ISYMIJ)

         CALL HEADER('Testing Lbar(ij,J)',-1)
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I3,A,I3)')
     &   'Symmetry of trial vector: ',ISYMX,
     &   'Picked Cholesky vector number',JVEC,' of sym.',ISYCHO
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Target norm   : ',TRGNRM,
     &   'Test   norm   : ',TSTNRM
         WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     &   'Min. abs. err.: ',ERRMIN,
     &   'Max. abs. err.: ',ERRMAX

 1002    CONTINUE

      ENDDO

      RETURN
      END
