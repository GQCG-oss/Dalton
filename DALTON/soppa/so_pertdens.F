C
      SUBROUTINE SO_PERTDENS(MODEL,CALCTYPE,NUMMAT,
     &                       FRVAL,NFRVAL,LABEL,
     &                       ISYMTR,IMAGPROP,FACTOR,
     &                       T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                       DENSAB,LDENSAB,DENSAI,LDENSAI,WORK,LWORK)
C
C     Rasmus Faber, November 2015
C
C     PURPOSE: Calculating the Response/perturbed density matrix
C              corresponding to the given solution/eigen vector.
C              This should ease the calculation of the final property.
C              In excitation energy calculations all eigevectors are
C              handled and both Excitation and De-Excitation parts are
C              written. For the calculation of response-functions only
C              one file is written.
C
C     INPUT:
C        MODEL    Determine which terms to include
C        CALCTYPE Signals which calculation type it is.
C        NUMMAT   Number of density matrices to calculate.
C        ISYMTR   Symmetry of the trial-vector
C        IMAGPROP Is the property imaginary
C        FACTOR   Factor with which to scale all contributions,
C                 use to correct for the fact that the solution vector
C                 is normalized in SO_OPTVEC (LR only)
C        T2AM     T2 amplitudes im complementary basis
C        DENSAB
C        DENSIJ   Second order densities (IA)
C        DENSAI
C
C     OUTPUT:
C        All written to disk
C
C
      use so_info, only: sop_linres, sop_excita, sop_stat_trh,
     &                   fn_rdens, fn_rdense, fn_rdensd,
     &                   so_has_doubles
C
      implicit none
C Get MULD2H (group multiplication table)
#include "ccorb.h"
C Symmetry-offsets in amplitudes
#include "ccsdsym.h"
C Symmetry-offsets in densities
#include "soppinf.h"
C Get irat
#include "iratdef.h"
C
C  Arguments
C
      CHARACTER*5, INTENT(IN) :: MODEL
      CHARACTER*8, INTENT(IN) :: LABEL
C
      LOGICAL,INTENT(IN) :: IMAGPROP
      INTEGER,INTENT(IN) ::  CALCTYPE, NUMMAT ,ISYMTR, NFRVAL,
     &                       LT2AM, LDENSIJ, LDENSAB, LDENSAI, LWORK
C
      DOUBLE PRECISION,INTENT(IN) :: T2AM(LT2AM),
     &                            DENSIJ(LDENSIJ), DENSAB(LDENSAB),
     &                            DENSAI(LDENSAI), FRVAL(NFRVAL),
     &                            FACTOR
C
      DOUBLE PRECISION, INTENT(INOUT) :: WORK(LWORK)
C

C
C     Local variables
      DOUBLE PRECISION EFACT, DFACT, DLABEL
      LOGICAL :: VALID, DOUBLES, STATIC
C
      INTEGER :: IDENS, IDUMMY
      INTEGER :: KSOLV1, KSOLV2, KRDENIJ, KRDENAB, ! Array positions
     &           KRDENAI
      INTEGER :: LSOLV1, LSOLV2, LRDENIJ, LRDENAB, LRDENAI,
     &           LRDENTOT                          ! Array lengths
      INTEGER :: KEND1, LWORK1
      INTEGER :: LURDENSE, LURDENSD     ! File handles
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
      DOUBLE PRECISION, PARAMETER :: SQRT2 = SQRT(2.0D0)
      LOGICAL :: DUMMY
      CALL QENTER('SO_PERTDENS')
C
C---------------------------------
C      Allocation of work space
C---------------------------------
C
      DOUBLES = SO_HAS_DOUBLES(MODEL)
      IF (.NOT.DOUBLES) THEN
         LSOLV1 = NT1AM(ISYMTR)
         LSOLV2 = 0
      ELSE
         LSOLV1 = NT1AM(ISYMTR)
         LSOLV2 = N2P2HOP(ISYMTR)
      ENDIF
      IF ( MODEL.EQ.'AORPA') THEN
         LRDENIJ = 0
         LRDENAB = 0
      ELSE
         LRDENIJ = NIJDEN(ISYMTR)
         LRDENAB = NABDEN(ISYMTR)
      ENDIF
      LRDENAI = NAIDEN(ISYMTR)
      LRDENTOT = LRDENIJ + LRDENAB + LRDENAI

      KRDENIJ = 1
      KRDENAB = KRDENIJ + LRDENIJ
      KRDENAI = KRDENAB + LRDENAB
      KSOLV1 = KRDENAI + LRDENAI
      KSOLV2 = KSOLV1 + LSOLV1
      KEND1  = KSOLV2 + LSOLV2
      LWORK1 = LWORK - KEND1

      CALL SO_MEMMAX ('SO_PERTDENS.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('SO_PERTDENS.1',' ',KEND1,LWORK)
C
      STATIC = (CALCTYPE.EQ.SOP_LINRES).AND.
     &         ABS(FRVAL(1)).LT.SOP_STAT_TRH
C
C--------------------------
C     Open the needed files
C--------------------------
C     Solution vectors
      CALL SO_OPEN(LUTR1E,FNTR1E,LSOLV1)
      IF (.NOT.STATIC) CALL SO_OPEN(LUTR1D,FNTR1D,LSOLV1)
      IF (DOUBLES) THEN
         CALL SO_OPEN(LUTR2E,FNTR2E,LSOLV2)
         IF(.NOT.STATIC) CALL SO_OPEN(LUTR2D,FNTR2D,LSOLV2)
      ENDIF
C     Density files
      LURDENSD = -1
      LURDENSE = -1
      SELECT CASE (CALCTYPE)
C        Linear response, only one file saved
         CASE ( SOP_LINRES )
            CALL GPOPEN(LURDENSD,FN_RDENS,'UNKNOWN','DIRECT',' ',
     &                  IRAT*LRDENTOT,DUMMY)
            LURDENSE = LURDENSD
C        Excitation energies, two files :
         CASE ( SOP_EXCITA )
            CALL GPOPEN(LURDENSE,FN_RDENSE,'UNKNOWN','DIRECT',' ',
     &                  IRAT*LRDENTOT,DUMMY)
            CALL GPOPEN(LURDENSD,FN_RDENSD,'UNKNOWN','DIRECT',' ',
     &                  IRAT*LRDENTOT,DUMMY)
C        Any other -- is an error
         CASE DEFAULT
            print *, 'Invalid CALCTYPE argument in so_pertdens'

      END SELECT
C
C------------------------------------
C     Setup factors for E and D parts
C------------------------------------
C     For static properties E and D contributions are the same:
C     Just calculate E with a factor of 2
      IF (STATIC) THEN
         EFACT = 2*FACTOR
      ELSE
         EFACT = FACTOR
      ENDIF
C
C     For real properties the D part has a minus, for excitation
C     energies, we keep both, so just keep a positive factor.
C
      IF ( IMAGPROP .OR. (CALCTYPE.EQ.SOP_EXCITA) ) THEN
         DFACT = FACTOR
      ELSE
         DFACT = -FACTOR
      ENDIF
C
      DO IDENS = 1, NUMMAT
         CALL DZERO( WORK(KRDENIJ), LRDENTOT)
C
C---------------------------------------
C        Process excitation contribution
C---------------------------------------
C
         CALL SO_READ(WORK(KSOLV1),LSOLV1,LUTR1E,FNTR1E,IDENS)
         IF (DOUBLES) CALL SO_READ(WORK(KSOLV2),LSOLV2,
     &                             LUTR2E,FNTR2E,IDENS)
C----------------------------
C        Singles contribution
C----------------------------
         CALL SO_PERTD1(MODEL,ISYMTR,EFACT,
     &               WORK(KRDENIJ),LRDENIJ,WORK(KRDENAB),LRDENAB,
     &               WORK(KRDENAI),LRDENAI,
     &               WORK(KSOLV1),LSOLV1,
     &               DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &               DENSAI,LDENSAI,WORK(KEND1),LWORK1)
C
C----------------------------------------------------
C        Doubles contributions to response densities.
C----------------------------------------------------
C
         IF (DOUBLES) THEN
C
            CALL SO_PERTD2(ISYMTR,EFACT,
     &                     WORK(KRDENIJ),LRDENIJ,WORK(KRDENAB),LRDENAB,
     &                     T2AM,LT2AM,WORK(KSOLV2),LSOLV2,
     &                     WORK(KEND1),LWORK1)
         END IF
C
         IF (CALCTYPE.EQ.SOP_EXCITA) THEN
            DLABEL = DBLE(IDENS)
            CALL SO_WRTVE(WORK(KRDENIJ),LRDENTOT,ISYMTR,LABEL,
     &                    DLABEL,LURDENSE)
            CALL DZERO( WORK(KRDENIJ), LRDENTOT)
         ENDIF
C
         IF (.NOT.STATIC) THEN
C
C------------------------------------------
C           Process de-excitation contribution
C------------------------------------------
C
            CALL SO_READ(WORK(KSOLV1),LSOLV1,LUTR1D,FNTR1D,IDENS)
            IF (DOUBLES) CALL SO_READ(WORK(KSOLV2),LSOLV2,
     &                                LUTR2D,FNTR2D,IDENS)
C-------------------------------
C           Singles contribution
C-------------------------------
            CALL SO_PERTD1(MODEL,ISYMTR,DFACT,
     &                     WORK(KRDENIJ),LRDENIJ,WORK(KRDENAB),LRDENAB,
     &                     WORK(KRDENAI),LRDENAI,
     &                     WORK(KSOLV1),LSOLV1,
     &                     DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                     DENSAI,LDENSAI,WORK(KEND1),LWORK1)
C
C----------------------------------------------------
C        Doubles contributions to response densities.
C----------------------------------------------------
C
            IF (DOUBLES) THEN
               CALL SO_PERTD2(ISYMTR,DFACT,
     &                        WORK(KRDENIJ),LRDENIJ,WORK(KRDENAB),
     &                        LRDENAB,T2AM,LT2AM,WORK(KSOLV2),LSOLV2,
     &                        WORK(KEND1),LWORK1)
            END IF
C
         END IF
C
C        For excitation calculations, label by number in irrep
C        (Since frequency based approach break down when we have
C         Degeneracies)
         IF (CALCTYPE.EQ.SOP_EXCITA) THEN
            DLABEL = DBLE(IDENS)
         ELSE
C           Else label by frequency of the field
            DLABEL = FRVAL(IDENS)
         ENDIF
C        Write on file
         CALL SO_WRTVE(WORK(KRDENIJ),LRDENTOT,ISYMTR,LABEL,
     &                  DLABEL,LURDENSD)
      END DO
C
C
C
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      IF(.NOT.STATIC) CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      IF(DOUBLES)THEN
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         IF(.NOT.STATIC) CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      ENDIF
C
      IF (CALCTYPE.EQ.SOP_EXCITA) CALL GPCLOSE(LURDENSE,'KEEP')
      CALL GPCLOSE(LURDENSD,'KEEP')

C
      CALL QEXIT('SO_PERTDENS')

      RETURN
      END SUBROUTINE
