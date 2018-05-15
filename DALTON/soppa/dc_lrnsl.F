C
C  /* Deck dc_lrnsl */
      SUBROUTINE DC_LRNSL(MODEL,ISYM,PROPAR,PROPAR1,PROPAR2,
     &                    DELETE_VECTORS,
     &                    FREQ,NFREQ,LABELS,
     &                    SOLV_LABELS,NLAB,
     &                    DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                    FOCKD,LFOCKD,T2AM,LT2AM,
     &                    IMAGINARY,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Anna Kristina Schnack-Petersen: April 2017
C
C     PURPOSE: Calculates the frequency dependent linear response
C              properties with an RPA(D) or HRPA(D) method with 
C              the atomic integral direct SOPPA program.
C
      use so_info, only: fn_rdens, sop_stat_trh, sop_dp,
     &                   sop_use_seller, so_has_doubles
      use so_data, only: fileinf
C
C#include "implicit.h"
      implicit none
#include "priunit.h"
C      
#include "soppinf.h"
#include "ccsdsym.h"
C REP <- Symbols of the ireducible representations.
C#include "pgroup.h"
#include "iratdef.h"
C

C
      integer, intent(in) :: isym, ! symmetry of this batch of perturbations
     &       LDENSIJ, LDENSAB, LDENSAI, LWORK, LFOCKD, LT2AM, ! array lengths
     &       NFREQ, NLAB
      REAL(sop_dp),INTENT(IN) :: DENSIJ(LDENSIJ), DENSAB(LDENSAB),
     &                     DENSAI(LDENSAI), FREQ(NFREQ),
     &                     FOCKD(LFOCKD), T2AM(LT2AM)
      REAL(sop_dp),INTENT(INOUT) :: WORK(LWORK)
      REAL(sop_dp),intent(inout) :: PROPAR(NFREQ,NLAB,NLAB),
     &                              PROPAR1(NFREQ,NLAB,NLAB),
     &                              PROPAR2(NFREQ,NLAB,NLAB)
      CHARACTER(LEN=8), intent(in) ::LABELS(NLAB)
      LOGICAL, intent(in) :: SOLV_LABELS(NLAB)
      CHARACTER(len=5), intent(in) :: MODEL
      LOGICAL, INTENT(IN) :: DELETE_VECTORS ! Whether to delete solution
                                            ! and residual vectors
      LOGICAL, INTENT(IN) :: IMAGINARY ! Wheter properties are imaginary
C
      CHARACTER*8 LABEL1, LABEL2
      CHARACTER*8 RTNLBL(2)
      REAL(sop_dp)  ::  PROPVAL, PROP1E, PROP1D, PROP2E, PROP2D,
     &                  DEN, XIDEN, Z1ABE,Z1BAE,Z1ABD,Z1BAD, GPFACTOR,
     &                  MFAC
      INTEGER :: LWORK1, LWORK2, LSOLV, LWORK3, L1VEC, L2VEC
C
      INTEGER :: KEND1, KEND2, KEND3, KSOLV, KGPVC1, KGPVC2,
     &           KTR1E, KTR1D, KSV1EA, KSV1DA, KSV1EB, KSV1DB,
     &           KRES1E, KRES1D, KRES2E, KRES2D, KRESO1E, KRESO1D,
     &           KRESO2E, KRESO2D, KEDIA2, KZ1E, KZ1D, KZ2E, KZ2D,
     &           KZ1EA, KZ1DA, KZ2EA, KZ2DA, KZ1EB, KZ1DB, KZ2EB, KZ2DB,
     &           KAssignedIndices
C
      INTEGER :: ILAB1, ILAB2, IFREQ, IDTYPE,
     &           maxnumjobs, LURDENS
C
      LOGICAL   IMAGPROP2, doubles, static, OLDDX
C
      REAL(sop_dp), PARAMETER :: HALF=0.5D0, ONE=1.0D0, TWO=2.0D0
C
      CHARACTER(len=*), parameter :: myname = 'DC_LRNSL'
      real(sop_dp) :: DDOT
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER(myname)
C
C     Initialize PROPAR
      PROPAR = 0.0D0
      PROPAR1 = 0.0D0
      PROPAR2 = 0.0D0
C
C---------------------------
C       Memory allocation 1
C---------------------------
C   As all singles vectors have the same length and all doubles vectors
C   likewise, these lengths are defined only once:
C
      L1VEC=NT1AM(ISYM)
      L2VEC=N2P2HOP(ISYM)
      KTR1E=1
      KTR1D=KTR1E+L1VEC
      KEND1=KTR1D+L1VEC
      LWORK1=LWORK-KEND1
C For MPI skal der allokeres plads til AssigendIndices
C---------------------------------
C       Memory allocation 2
C--------------------------------
      KGPVC1=KEND1
      KGPVC2=KGPVC1+L1VEC
      KRES1E=KGPVC2+L2VEC
      KRES1D=KRES1E+L1VEC
      KRES2E=KRES1D+L1VEC
      KRES2D=KRES2E+L2VEC
      KRESO1E=KRES2D+L2VEC
      KRESO1D=KRESO1E+L1VEC
      KRESO2E=KRESO1D+L1VEC
      KRESO2D=KRESO2E+L2VEC
      KEND2=KRESO2D+L2VEC
      LWORK2=LWORK-KEND2
C
      CALL SO_MEMMAX(myname//'.1',LWORK2)
      IF (LWORK2.LT.0) CALL STOPIT(myname//'.1',' ',KEND2,LWORK)
C-------------------------------
C   Begin loop over frequecies:
C-------------------------------
      DO IFREQ=1,NFREQ
C
         STATIC = ABS(FREQ(IFREQ)).LT.sop_stat_trh
         CALL SO_OPEN(LUTR1E,FNTR1E,L1VEC)
         IF (.NOT.STATIC) CALL SO_OPEN(LUTR1D,FNTR1D,L1VEC)
C
         LSOLV = L1VEC
         IF (.NOT.STATIC) LSOLV = 2*LSOLV
         DO ILAB1=1,NLAB
            LABEL1=LABELS(ILAB1)
            CALL SO_GETSOLV(WORK(KTR1E),LSOLV,LABEL1,FREQ(IFREQ),1,
     &                      STATIC,.FALSE.,ISYM)
C-------------------------------------------------------------------------
C   Save solution vectors in trial vector file, so Reduced E matrix is
C   calculated correctly.
C-------------------------------------------------------------------------
            CALL SO_WRITE(WORK(KTR1E),L1VEC,LUTR1E,FNTR1E,ILAB1)
            IF (.NOT.STATIC) CALL SO_WRITE(WORK(KTR1D),L1VEC,LUTR1D,
     &                       FNTR1D,ILAB1)
         END DO
C
         CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
         IF (.NOT.STATIC) CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
C
C--------------------------------------
C    Determine Linearly transformed -E
C--------------------------------------
         IF (STATIC)THEN
            IF (IMAGINARY) THEN
               IDTYPE=2
            ELSE
               IDTYPE=1
            END IF
         ELSE
            IDTYPE=0
         END IF
C
         CALL SO_ERES(MODEL, 0, NLAB, DENSIJ, LDENSIJ,
     &                DENSAB, LDENSAB,
     &                T2AM, LT2AM, FOCKD, LFOCKD, DENSAI, LDENSAI,
     &                 IFREQ, ISYM,IDTYPE,
#ifdef VAR_MPI
     &                WORK(KAssignedIndices), maxnumjobs,
#endif
     &                WORK(1),LWORK)
C
C----------------------------------------------------------------------------
C   Open result vector files and in the dynamic case also trial vector files
C----------------------------------------------------------------------------
C
         CALL SO_OPEN(LURS1E,FNRS1E,L1VEC)
         CALL SO_OPEN(LURS2E,FNRS2E,L2VEC)
         CALL SO_OPEN(LUTR1E,FNTR1E,L1VEC)
         IF (.NOT.STATIC) THEN
            CALL SO_OPEN(LURS1D,FNRS1D,L1VEC)
            CALL SO_OPEN(LURS2D,FNRS2D,L2VEC)
            CALL SO_OPEN(LUTR1D,FNTR1D,L1VEC)
         END IF
C
         DO ILAB1=1,NLAB

            IF (.NOT.SOLV_LABELS(ILAB1)) CYCLE
            LABEL1=LABELS(ILAB1)
C
C---------------------------------------
C   Read files into array
C---------------------------------------
C
            IF (MODEL .NE. 'DCHRP') THEN
C----------------------------------------------------------
C The singles result vectors are only used when not HRPA(D)
C----------------------------------------------------------
               CALL SO_READ(WORK(KRES1E),L1VEC,LURS1E,FNRS1E,ILAB1)
               IF (.NOT.STATIC) THEN
                  CALL SO_READ(WORK(KRES1D),L1VEC,LURS1D,FNRS1D,ILAB1)

               END IF
            END IF
            CALL SO_READ(WORK(KRES2E),L2VEC,LURS2E,FNRS2E,ILAB1)
            IF (.NOT.STATIC)THEN
               CALL SO_READ(WORK(KRES2D),L2VEC,LURS2D,FNRS2D,ILAB1)
            END IF
C----------------------------------------------------------
C    Determine 1p-1h and 1h-1p part of reduced S[2] matrix
C     - not necessary for HRPA(D) model
C----------------------------------------------------------
C
            IF ((.NOT.STATIC).AND.(MODEL.NE.'DCHRP')) THEN
               CALL SO_READ(WORK(KTR1E),L1VEC,LUTR1E,FNTR1E,ILAB1)
               CALL SO_READ(WORK(KTR1D),L1VEC,LUTR1D,FNTR1D,ILAB1)
               CALL SO_RES_O(WORK(KRESO1E),L1VEC,WORK(KRESO1D),L1VEC,
     &                      WORK(KTR1E),L1VEC,WORK(KTR1D),L1VEC,
     &                      DENSIJ,LDENSIJ,DENSAB,LDENSAB,ISYM,ISYM)
               CALL DSCAL(L1VEC,-ONE,WORK(KRESO1E),1)
               CALL DSCAL(L1VEC,-ONE,WORK(KRESO1D),1)
            END IF
C---------------------------------
C   Determine property gradients
C---------------------------------
            CALL SO_GETGP(WORK(KGPVC1),L1VEC,WORK(KGPVC2),L2VEC,
     &                    LABEL1,ISYM,IMAGPROP2,MODEL,T2AM,LT2AM,
     &                    DENSIJ,LDENSIJ,DENSAB,LDENSAB,
     &                    DENSAI,LDENSAI,WORK(KEND2),LWORK2)
            IF (.NOT.TRIPLET) THEN
               CALL SO_TMLTR(WORK(KGPVC2),HALF,ISYM)
            END IF
            IF (IMAGPROP2.NEQV.IMAGINARY)THEN
              CALL QUIT('Error in DC_LRNSL: Imaginary unclear')
            END IF
C--------------------------------------------------------------------------
C   Determine Y=(Omega*S-E)*X and Z for the singles part:
C   Z is the sum of the property gradient and -1/2*Y
C   - Not necessary for HRPA(D).
C--------------------------------------------------------------------------
            IF ((.NOT.STATIC).AND.(MODEL.NE.'DCHRP')) THEN
               CALL DAXPY(L1VEC,FREQ(IFREQ),WORK(KRESO1E),1,
     &                    WORK(KRES1E),1)
               CALL DAXPY(L1VEC,FREQ(IFREQ),WORK(KRESO1D),1,
     &                    WORK(KRES1D),1)
            END IF
C-------------------------------------------------
C   Determine gradient property deexcitation part 
C-------------------------------------------------
            IF (IMAGINARY) THEN
               GPFACTOR=ONE
            ELSE
               GPFACTOR=-ONE
            END IF
C---------------------------------------------
C   Determine singles Z vectors for the label
C-----------------------------------------------------------
C   For the HRPA(D) model we need only the property gradient
C------------------------------------------------------------
            IF (MODEL.EQ.'DCRPA')THEN
               MFAC = -HALF
            ELSE IF (MODEL.EQ.'DCHRP')THEN
               MFAC = 0.0D0
            END IF
            CALL DSCAL(L1VEC,MFAC,WORK(KRES1E),1)
            CALL DAXPY(L1VEC,ONE,WORK(KGPVC1),1,WORK(KRES1E),1)
            IF (.NOT.STATIC) THEN
               CALL DSCAL(L1VEC,MFAC,WORK(KRES1D),1)
               CALL DAXPY(L1VEC,GPFACTOR,WORK(KGPVC1),1,WORK(KRES1D),1)
            END IF
C
C-----------------------------------------------------
C   Write these new Z resultvectors to file:
C   The result vector singles files are used
C-----------------------------------------------------
C
            CALL SO_WRITE(WORK(KRES1E),L1VEC,LURS1E,FNRS1E,ILAB1)
            IF (.NOT.STATIC)THEN
               CALL SO_WRITE(WORK(KRES1D),L1VEC,LURS1D,FNRS1D,ILAB1)
            END IF
C
C-------------------------------------------
C   Determine Z vector for the doubles part
C-------------------------------------------
            CALL DSCAL(L2VEC,-ONE,WORK(KRES2E),1)
            CALL DAXPY(L2VEC,ONE,WORK(KGPVC2),1,WORK(KRES2E),1)
            IF (.NOT.STATIC) THEN
               CALL DSCAL(L2VEC,-ONE,WORK(KRES2D),1)
               CALL DAXPY(L2VEC,GPFACTOR,WORK(KGPVC2),1,WORK(KRES2D),1)
            END IF
C
C-----------------------------------------------------
C   Write these new Z resultvectors to file:
C   The result vector doubles files are used
C-----------------------------------------------------
C
            CALL SO_WRITE(WORK(KRES2E),L2VEC,LURS2E,FNRS2E,ILAB1)
            IF (.NOT.STATIC) THEN
               CALL SO_WRITE(WORK(KRES2D),L2VEC,LURS2D,FNRS2D,ILAB1)
            END IF
C
         END DO
C---------------------------------------------------
C   Calculate the property for the given frequency
C---------------------------------------------------
C   Memory allocation 3:
C------------------------
         KEDIA2=1
         KZ1EA=KEDIA2+L2VEC
         KZ1DA=KZ1EA+L1VEC
         KZ2EA=KZ1DA+L1VEC
         KZ2DA=KZ2EA+L2VEC
         KSV1EA=KZ2DA+L2VEC
         KSV1DA=KSV1EA+L1VEC
         KZ1EB=KSV1DA+L1VEC
         KZ1DB=KZ1EB+L1VEC
         KZ2EB=KZ1DB+L1VEC
         KZ2DB=KZ2EB+L2VEC
         KSV1EB=KZ2DB+L2VEC
         KSV1DB=KSV1EB+L1VEC
         KEND3=KSV1DB+L1VEC
         LWORK3=LWORK-KEND3
         CALL SO_MEMMAX(myname//'.2',LWORK3)
         IF (LWORK3 .LT. 0) CALL STOPIT(myname//'.2',' ',KEND3,LWORK)
C--------------------------------------------------------------------
C   Determine Diagonal E[2] elements for determining the D[0] matrix:
C--------------------------------------------------------------------
         IF (TRIPLET) THEN
            CALL SO_EDIAG2T(WORK(KEDIA2),L2VEC,FOCKD,LFOCKD,
     &                      ISYM,WORK(KEND3),LWORK3)
         ELSE
            CALL SO_EDIAG2(WORK(KEDIA2),L2VEC,FOCKD,LFOCKD,
     &                     ISYM,WORK(KEND3),LWORK3)
         END IF
C
C
C-------------------------------
C   Read the necessary vectors:
C------------------------------
         DO ILAB1=1,NLAB
            IF (.NOT.SOLV_LABELS(ILAB1)) CYCLE
            LABEL1=LABELS(ILAB1)
            CALL SO_READ(WORK(KZ1EA),L1VEC,LURS1E,FNRS1E,ILAB1)
C-------------------------------------------------------
C   Read different vectors for the two models:
C   RPA(D) requires solution vector and both Z-vectors
C   HRPA(D) requires only the doubles Z-vector 
C   (single Z-vector is property gradient).
C-------------------------------------------------------
            IF (MODEL.EQ.'DCRPA')THEN
               CALL SO_READ(WORK(KSV1EA),L1VEC,LUTR1E,FNTR1E,ILAB1)
               CALL SO_READ(WORK(KZ2EA),L2VEC,LURS2E,FNRS2E,ILAB1)
            ELSE IF (MODEL.EQ.'DCHRP')THEN
               CALL SO_READ(WORK(KZ2EA),L2VEC,LURS2E,FNRS2E,ILAB1)
            END IF
C---------------------------------------------------------------------------
C   In the case of the singlet state (one of) the doubles Z vectors must be
C   modified by spin adaption.
C---------------------------------------------------------------------------
            IF (.NOT.TRIPLET) CALL CCSD_TCMEPKX(WORK(KZ2EA),TWO,ISYM)
C
            IF (.NOT.STATIC) THEN
               CALL SO_READ(WORK(KZ1DA),L1VEC,LURS1D,FNRS1D,ILAB1)
               IF(MODEL.EQ.'DCRPA')THEN
                  CALL SO_READ(WORK(KSV1DA),L1VEC,LUTR1D,FNTR1D,ILAB1)
                  CALL SO_READ(WORK(KZ2DA),L2VEC,LURS2D,FNRS2D,ILAB1)
               ELSE IF (MODEL.EQ.'DCHRP')THEN
                  CALL SO_READ(WORK(KZ2DA),L2VEC,LURS2D,FNRS2D,ILAB1)
               END IF

C---------------------------------------------------------------------------
C   In the case of the singlet state (one of) the doubles Z vectors must be
C   modified by spin adaption.
C---------------------------------------------------------------------------
               IF (.NOT.TRIPLET) CALL CCSD_TCMEPKX(WORK(KZ2DA),TWO,ISYM)
            END IF
C
            DO ILAB2=ILAB1,NLAB
               IF (.NOT.SOLV_LABELS(ILAB2)) CYCLE
               LABEL2=LABELS(ILAB2)
               CALL SO_READ(WORK(KSV1EB),L1VEC,LUTR1E,FNTR1E,ILAB2)
C------------------------------------------------------------
C Both single excitation 'Z vectors' are not used in HRPA(D)
C------------------------------------------------------------
               IF (MODEL.NE.'DCHRP')THEN
                  CALL SO_READ(WORK(KZ1EB),L1VEC,LURS1E,FNRS1E,ILAB2)
               END IF
               CALL SO_READ(WORK(KZ2EB),L2VEC,LURS2E,FNRS2E,ILAB2)
               IF (.NOT.STATIC) THEN
                  CALL SO_READ(WORK(KSV1DB),L1VEC,LUTR1D,
     &                         FNTR1D,ILAB2)
                  IF (MODEL.NE.'DCHRP')THEN
                     CALL SO_READ(WORK(KZ1DB),L1VEC,LURS1D,FNRS1D,ILAB2)
                  END IF
                  CALL SO_READ(WORK(KZ2DB),L2VEC,LURS2D,FNRS2D,ILAB2)
               END IF
C------------------------------------------------------
C   Determine Singles part of polarization propagator:
C------------------------------------------------------
               IF (MODEL.EQ.'DCRPA')THEN
                  Z1ABE = DDOT(L1VEC,WORK(KSV1EA),1,WORK(KZ1EB),1)
                  Z1BAE = DDOT(L1VEC,WORK(KSV1EB),1,WORK(KZ1EA),1)
                  PROP1E=Z1ABE+Z1BAE
                  IF (.NOT.STATIC) THEN
                     Z1ABD = DDOT(L1VEC,WORK(KSV1DA),1,WORK(KZ1DB),1)
                     Z1BAD = DDOT(L1VEC,WORK(KSV1DB),1,WORK(KZ1DA),1)
                     PROP1D=Z1BAD+Z1ABD
                  ELSE
                     PROP1D=PROP1E
                  END IF
               ELSE IF (MODEL.EQ.'DCHRP')THEN
C------------------------------------------------------------------
C For HRPA(D) PROP1E is the dot-product of the solution vector and
C the single excitation part of the property gradient (saved in
C WORK(KZ1EA)
C-------------------------------------------------------------------
                  PROP1E = DDOT(L1VEC,WORK(KSV1EB),1,WORK(KZ1EA),1)
                  IF(.NOT.STATIC)THEN
                     PROP1D = DDOT(L1VEC,WORK(KSV1DB),1,WORK(KZ1DA),1)
                  ELSE
                     PROP1D = PROP1E
                  END IF
               END IF
C
C-----------------------------------------------------
C   Determine doubles part of polarization propagator:
C-----------------------------------------------------
               PROP2E=0.0D0
               DO I=1,L2VEC
                  DEN= (WORK(KEDIA2+I-1)-FREQ(IFREQ))
                  XIDEN=ONE/DEN
                  PROP2E=PROP2E+WORK(KZ2EA+I-1)*WORK(KZ2EB+I-1)*XIDEN
               END DO
               ! Do we need this?
               IF (.NOT.TRIPLET) THEN
                  PROP2E = HALF * PROP2E
               END IF
C
               IF(.NOT.STATIC)THEN
                  PROP2D=0.0D0
                  DO I=1,L2VEC
                     DEN= (WORK(KEDIA2+I-1)+FREQ(IFREQ))
                     XIDEN=ONE/DEN
                     PROP2D=PROP2D+WORK(KZ2DA+I-1)*WORK(KZ2DB+I-1)*XIDEN
                  END DO
                  IF(.NOT.TRIPLET) THEN
                     PROP2D = HALF * PROP2D
                  END IF
               ELSE
                  PROP2D=PROP2E
               END IF
C--------------------------------------------
C   Second order property:
C--------------------------------------------
               PROPVAL=PROP1E+PROP1D+PROP2E+PROP2D
               PROPAR(IFREQ,ILAB1,ILAB2)=PROPVAL
               PROPAR(IFREQ,ILAB2,ILAB1)=PROPVAL
               PROPAR1(IFREQ,ILAB1,ILAB2)=PROP1E+PROP1D
               PROPAR1(IFREQ,ILAB2,ILAB1)=PROP1E+PROP1D
               PROPAR2(IFREQ,ILAB1,ILAB2)=PROP2E+PROP2D
               PROPAR2(IFREQ,ILAB2,ILAB1)=PROP2E+PROP2D
            END DO !Over Label2
         END DO !Over Label1
C
C---------------------------------------------
C   Close result vector and trial vector files
C---------------------------------------------
C
         CALL SO_CLOSE(LURS1E,FNRS1E,L1VEC)
         CALL SO_CLOSE(LURS2E,FNRS2E,L2VEC)
         CALL SO_CLOSE(LUTR1E,FNRS1E,L1VEC)
         IF (.NOT.STATIC) THEN
            CALL SO_CLOSE(LURS1D,FNRS1D,L1VEC)
            CALL SO_CLOSE(LURS2D,FNRS2D,L2VEC)
            CALL SO_CLOSE(LUTR1D,FNTR1D,L1VEC)
         END IF
      END DO !Over frequencies
C
C---------------------------------------------------------
C     Delete files for this symmetry and dump fileinf list
C---------------------------------------------------------
C
      LURDENS = -1
      CALL GPOPEN(LURDENS, FN_RDENS, 'OLD','DIRECT','UNFORMATTED',
     &               IRAT*LDENSAI,OLDDX)
      IF (DELETE_VECTORS) THEN
         CALL GPCLOSE(LURDENS,'DELETE')
         STATIC = .TRUE.
         DO IFREQ = 1, NFREQ
            STATIC = STATIC.AND.(abs(freq(ifreq)).lt.sop_stat_trh)
         END DO
         CALL SO_DELVEC(STATIC,.FALSE.,ISYM)
         call fileinf%empty
      ELSE
         CALL GPCLOSE(LURDENS,'KEEP')
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT(myname)
      RETURN
      END

      SUBROUTINE DC_RSPOUT(PROPAR,PROPAR1,PROPAR2,ISYM,
     &                     FREQ,NFREQ,LABELS,NLAB)

      USE SO_INFO, ONLY: SOP_DP, SOP_STAT_TRH
      IMPLICIT NONE
C
C Get LUPRI, output unit (probably 6)      
#include "priunit.h"
C From codata we need varius conversion factors, XT*. 
#include "codata.h"
C inflr.h below requires MAXLBL from rspprp.h
#include "rspprp.h"
C We need 
C LBLLR <- Label of operators 
C NGPLR <- Number of operators of each symmetry
#include "inflr.h"
C TRIPLET <- Triplet property flag
#include "soppinf.h"

      REAL(SOP_DP), INTENT(IN) :: PROPAR(NFREQ,NLAB,NLAB),
     &                            FREQ(NFREQ),
     &                            PROPAR1(NFREQ,NLAB,NLAB),
     &                            PROPAR2(NFREQ,NLAB,NLAB)
      INTEGER, INTENT(IN) :: ISYM, NFREQ, NLAB
      CHARACTER(LEN=8) :: LABELS(NLAB)
      
      INTEGER :: ILAB1, ILAB2, IFREQ

      CALL HEADER('Final output of second order properties from'//
     &            ' linear response',-1)
      IF (TRIPLET) THEN
         WRITE (LUPRI,'(/A)') '@ Spin symmetry of operators: triplet'
      ELSE
         WRITE (LUPRI,'(/A)') '@ Spin symmetry of operators: singlet'
      END IF
      WRITE (LUPRI,'(/A/A)')
     &   ' Note that minus the linear response function:'//
     &   ' - << A; B >>(omega) is printed.'
      DO IFREQ = 1, NFREQ
         IF ( ABS(FREQ(IFREQ)) .LT. sop_stat_trh) THEN
            WRITE(LUPRI,'(/A/)')
     *      '@ FREQUENCY INDEPENDENT SECOND ORDER PROPERTIES'
         ELSE
            WRITE(LUPRI,'(/A/,5(/A,1P,D15.7),/)')
     * '@ FREQUENCY DEPENDENT SECOND ORDER PROPERTIES WITH FREQUENCY :',
     *      '@    a.u.:',FREQ(IFREQ),
     *      '@    cm-1:',XTKAYS*FREQ(IFREQ),
     *      '@    eV  :',XTEV*FREQ(IFREQ),
     *      '@  kJ/mol:',XKJMOL*FREQ(IFREQ),
     *      '@    nm  :',XTNM/FREQ(IFREQ)
         END IF
         DO ILAB1 = 1, NLAB
            DO ILAB2 = ILAB1, NLAB
               WRITE(LUPRI,'(A)')'Contribution from single excitations:'
               WRITE(LUPRI,'(5A,1P,D20.12)')
     *         '@ -<< ',LABELS(ILAB1),' ; ',LABELS(ILAB2),
     &         ' >> =',PROPAR1(IFREQ,ILAB1,ILAB2)
               WRITE(LUPRI,'(A)')'Contribution from double excitations:'
               WRITE(LUPRI,'(5A,1P,D20.12)')
     *         '@ -<< ',LABELS(ILAB1),' ; ',LABELS(ILAB2),
     &         ' >> =',PROPAR2(IFREQ,ILAB1,ILAB2)
               WRITE(LUPRI,'(A)')'Total second order property:'
               WRITE(LUPRI,'(5A,1P,D20.12)')
     *         '@ -<< ',LABELS(ILAB1),' ; ',LABELS(ILAB2),
     &         ' >> =',PROPAR(IFREQ,ILAB1,ILAB2)
            END DO
         END DO
      END DO
               
      END SUBROUTINE
