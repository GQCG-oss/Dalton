C
C   Triplet versions of Xi and Eta routines
C
c*DECK CC_ETAC13
      SUBROUTINE CC_ETAC13(ISYMC,LBLC,        !Operator stuff
     &                     ETAC,                !Result vector
     &                     LIST,ILSTNR,         !Left vector
     &                     LEOM,
     &                     XINT,WORK,LWORK)
C
C-----------------------------------------------------------------------
C     Purpose: Calculate the etaC(L) vector, when L is a triplet vector
C     and C is a singlet operator.
C     This means that the result vector will be a triplet vector.
C     (Modified version of the CCSD etaC(l0/l1) code
C
C     eta(tau_nu)= (<HF| + Sum(mu)L(1)<mu|)
C                         exp(-t)[C,tau_nu]exp(T)|HF>
C
C     Input:
C
C
C     LIST= 'L0' for zeroth order left amplitudes.
C                ISYML should be ISYMOP in this case.
C
C           'L1' for first order left amplitudes, read in from file
C                In this case the vector is found according to its list
C                number ILSTNR.
C
C                For L1 HF contribution is skipped.
C
C     C property integrals read according to LBLC
C
C     SLV98,OC: Allow for input of integrals if
C               LBLC.eq.'GIVE INT'
C
C     Sonia & Filip, Maj 2015
C
C-----------------------------------------------------------------------
C
      implicit none

      character(len=*), parameter :: myname = 'CC_ETAC13'
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccexci.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
!Sonia & Filip: find a better place
#include "ccsections.h"
!sonia
#include "second.h"

      DOUBLE PRECISION, PARAMETER :: TWO  = 2.0D0, HALF = 0.5D0,
     &                               ZERO = 0.0D0, ONEM = -1.0D0
      LOGICAL, PARAMETER :: LOCDBG = .false.
C
C     Input:
C
      INTEGER, INTENT(IN) :: ISYMC,  !  Symmetry of C
     &                       ILSTNR, ! List number for input
     &                       LWORK   !
      CHARACTER(LEN=*), INTENT(IN) ::  LBLC, !
     &                                 LIST  ! See above
      DOUBLE PRECISION, INTENT(OUT)::  ETAC(*), ! Output array
     &                                 WORK(LWORK)  ! Work array
      DOUBLE PRECISION, INTENT(INOUT) :: XINT(*)
      LOGICAL, INTENT(IN) :: LEOM ! Whether to include extra EOM-CC
                                  ! Contributions

      CHARACTER MODEL*10
      LOGICAL :: FCKCON, ETRAN
      INTEGER :: KT1AM, KT2AM, KL1AM, KL2AM, KCTMO, KLAMDP, KLAMDH,
     &           KEI1, KEI2, KEND1, KEND2, KEND21, KEND3, KEND4,
     &           KINTAI, KINTIA, KXMAT, KYMAT, KETAC, KCINT
      INTEGER :: KOFF1, KOFF2, KOFFP, KOFFM
      INTEGER :: ISYML, ISYRES, ISYMA, ISYMI
      INTEGER :: LEND1, LEND2, LEND21, LEND3, LEND4
      INTEGER :: IOPT
      DOUBLE PRECISION :: FF, XXI, XYI, ETA1, ETA2
      DOUBLE PRECISION :: TIMEC
C     EXTERNAL FUNCTIONS
      DOUBLE PRECISION :: DDOT
      INTEGER :: ILSTSYM
C
      CALL AROUND( 'Constructing ETA^{'// LBLC
     &              //'}('// LIST //') vector ')
C      IF ( (IPRINT .GT. 10).or.locdbg ) THEN
C         CALL AROUND( 'IN CCCI_ETAC: Constructing ETAC(LE) vector ')
C      ENDIF
C
C--------------------------------
C     find symmetry of D operator.
C--------------------------------
C
      ISYML = ILSTSYM(LIST,ILSTNR)
C
      ISYRES = MULD2H(ISYML,ISYMC)
      IF (( LIST .EQ. 'L0')) THEN
         CALL QUIT('Misuse of '//myname)
      ENDIF
C
      TIMEC = SECOND()
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
C--------------------
C     Allocate space.
C--------------------
C
      KCTMO  = 1
      KT1AM  = KCTMO  + N2BST(ISYMC)
      KLAMDP = KT1AM  + NT1AM(1)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
C
      LEND1  = LWORK  - KEND1
C
      IF ( .NOT. CCS) THEN

         KINTIA = KEND1
         KEND1  = KINTIA + NT1AM(ISYMC)
         LEND1  = LWORK - KEND1
         CALL DZERO(WORK(KintIA),NT1AM(isymc))
C
         KL1AM = KEND1
         KL2AM = KL1AM + NT1AM(ISYML)
         KEND2 = KL2AM + NT2SQ(ISYML)
         LEND2 = LWORK - KEND2
         KT2AM = KEND2
         KEND21= KT2AM + MAX(NT2AM(1),NT2AM(ISYML))
         LEND21= LWORK - KEND21
C
      ELSE
C
         KL1AM = KEND1
         KEND2 = KL1AM + NT1AM(ISYML)
         LEND2 = LEND1
         KEND21= KEND1
         LEND21= LEND1
C
      ENDIF
      KEI1   = KEND21
      KEI2   = KEI1   + NEMAT1(ISYMC)
      KEND3  = KEI2   + NMATIJ(ISYMC)
      LEND3  = LWORK  - KEND3
      IF (LEND3.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN '//myname)
C
C-----------------------
C     get T1 amplitudes.
C-----------------------
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      IF ( .NOT. CCS) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND21),LEND21)
C
C-------------------------------
C     get AO property integrals.
C-------------------------------
C
      CALL DZERO(WORK(KCTMO),N2BST(ISYMC))
      FF = 1.0D0
C SLV98,OC give integrals option
      IF (LBLC.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMC),XINT(1),1,WORK(KCTMO),1)
      ELSE
        FF = 1.0D0
        CALL CC_ONEP(WORK(KCTMO),WORK(KEND21),LEND21,FF,ISYMC,LBLC)
      ENDIF
C
C-----------------------------------------------
C     Make MO T1-transformed property integrals.
C-----------------------------------------------
C
      CALL CC_FCKMO(WORK(KCTMO),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND21),LEND21,ISYMC,1,1)
C
C----------------------------------------------------------
C     Extract Cia (stored ia) and reorder ai
C----------------------------------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMC)
C
         DO 110 A = 1,NVIR(ISYMA)
C
            DO 120 I = 1,NRHF(ISYMI)
C
               KOFF1 = KINTIA + IT1AM(ISYMA,ISYMI)
     &               + NVIR(ISYMA)*(I - 1) + A-1
               KOFF2 = KCTMO + IFCVIR(ISYMI,ISYMA)
     *               + NORB(ISYMI)*(A - 1) + I - 1
C
               WORK(KOFF1) = WORK(KOFF2)
C
  120       CONTINUE
  110    CONTINUE
C
  100 CONTINUE
C
C     Initialize ETAC
      CALL DZERO(ETAC,NT1AM(ISYRES)+2*NT2AM(ISYRES))

C----------------------------------------------
C     Read L1, and L2(+) multipliers/left vectors.
C----------------------------------------------
C
      IOPT = 65
      CALL CC_RDRSP(LIST,ILSTNR,ISYML,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KT2AM))
C
C--------------------------------
C     Put C into E matrix format.
C--------------------------------
C
      FCKCON = .TRUE.
      ETRAN  = .FALSE.
      CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                WORK(KCTMO),WORK(KEND3),LEND3,FCKCON,
     *                ETRAN,ISYMC)
C
C--------------------------------------------
C     etac1 =  sum(b)Lbi*Cba - sum(j)Laj*Cij.
C--------------------------------------------
C
      CALL CCLR_E1C1(ETAC,WORK(KL1AM),WORK(KEI1),WORK(KEI2),
     *               WORK(KEND3),LEND3,ISYML,ISYMC,'T')
C                   ~
C Square L2(+) onto L2
      IF (.NOT. CCS) THEN
         CALL CCRHS3_R2IJ(WORK(KT2AM),WORK(KEND3),LEND3,ISYML)
         CALL CC_T2SQ(WORK(KT2AM),WORK(KL2AM),ISYML)
C---------------------------------
C Calculate the doubles (+)G term
C---------------------------------
C Transpose i,j and a,b blocks of C
         CALL CC_EITR(WORK(KEI1),WORK(KEI2),WORK(KEND3),LEND3,ISYMC)
C Scale C by 1/2, as this term has a factor of 1/2
         CALL DSCAL(NMATAB(ISYMC),HALF,WORK(KEI1),1)
         CALL DSCAL(NMATIJ(ISYMC),HALF,WORK(KEI2),1)
C Calculate actual term
         KOFFP = NT1AM(ISYRES) + 1
         CALL CCRHS_E(ETAC(KOFFP),WORK(KL2AM),WORK(KEI1),WORK(KEI2),
     &                WORK(KEND3),LEND3,ISYML,ISYMC)
C Restore factor of C
         CALL DSCAL(NMATAB(ISYMC),TWO,WORK(KEI1),1)
         CALL DSCAL(NMATIJ(ISYMC),TWO,WORK(KEI2),1)
      END IF
C
C----------------------------------------------
C     Read L2(-) multipliers/left vector.
C----------------------------------------------
C     The routine needs to be modified to allow this
      IOPT = 128
      CALL CC_RDRSP(LIST,ILSTNR,ISYML,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KT2AM))
C
C We need L(+) - L(-) for X and Y in the following.
      CALL CC_T2SQ3A(WORK(KT2AM),WORK(KL2AM),ISYML,ONEM)
C
C------------------------------------------
C     Put T2 (packed) amplitudes in etac2.
C------------------------------------------
C
      IF (.NOT. CCS) THEN
         !read in T2AM (packed)
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
      ENDIF
C
C--------------------------------
C     Make X and Y intermediates.
C--------------------------------
CRF Triplet requires (L(+) - L(-))
      IF (.NOT. CCS) THEN
         KXMAT = KEND3
         KYMAT = KXMAT + NMATIJ(ISYML)
         KEND4 = KYMAT + NMATAB(ISYML)
         LEND4 = LWORK - KEND4
         IF (LEND4.LT. 0 )
     &        CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-2')
C
         IF ( DEBUG.or.LOCDBG ) THEN
            XYI   = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L1AM vector:              ',XYI
            XYI   = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L2AM vector:              ',XYI
            XXI   = DDOT(NT2AM(ISYMOP),WORK(KT2AM),1,WORK(KT2AM),1)
            WRITE(LUPRI,1) 'T2AM vector :                      ',XXI
         ENDIF
         CALL CC_XI(WORK(KXMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND4),LEND4)
         CALL CC_YI(WORK(KYMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND4),LEND4)
         IF ( DEBUG.or.LOCDBG ) THEN
            XYI   = DDOT(NMATAB(ISYML),WORK(KYMAT),1,WORK(KYMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: YI  intermediate is:      ',XYI
            XXI   = DDOT(NMATIJ(ISYML),WORK(KXMAT),1,WORK(KXMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: XI  intermediate is:      ',XXI
         ENDIF
      ELSE
         KEND4 = KEND3
         LEND4 = LEND3
      ENDIF
C
C----------------------------------------------
C     Calculate X and Y contributions to etac1.
C     etac1 = -sum(e)Cie*Yae - sum(l)Cla*Xli
C----------------------------------------------
C
      IF ( .NOT.CCS ) THEN
C
         CALL CC_21EFM(ETAC,WORK(KCTMO),ISYMC,WORK(KXMAT),
     *                 WORK(KYMAT),ISYML,WORK(KEND4),LEND4)
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            WRITE(LUPRI,1) 'Norm of eta1-after X&Y cont:       ',ETA1
         ENDIF
      ENDIF
C We have L(+)-L(-), we need L(+)+L(-)
      CALL CC_T2SQTRANSP(WORK(KL2AM),ISYML)
C
      IF (LEOM) THEN
C
C---------------------------------------------------------------
C     EOM contribution to ETAC:                     ~
C     etac = 2sum(ei)(L(+) + L(-))_{ckei}*(C_{ei} + t_{ei,fn}*C_{nf})
C---------------------------------------------------------------
         KCINT = KEND3
         KETAC = KCINT
         KINTAI = KCINT + MAX(NT1AM(ISYMC),NT1AM(ISYRES))
C        Extract C_{ai}
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMC)
            DO I = 1, NRHF(ISYMI)
               KOFF1 = KINTAI + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)
               KOFF2 = KCTMO + IFCRHF(ISYMA,ISYMI) +
     &                 NORB(ISYMA)*(I-1) + NRHF(ISYMA)
               CALL DCOPY(NVIR(ISYMA),WORK(KOFF2),1,WORK(KOFF1),1)
            END DO
         END DO
         !get ttilde (overwrites T2am)
         CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,1,1)
C            ~
C        Add t_{em,fn} *C_{nf} to C_{em}
         CALL CCG_LXD(WORK(KCINT),ISYMC,WORK(KINTIA),ISYMC,
     &                WORK(KT2AM),1,0)
         CALL DAXPY(NT1AM(ISYMC),1.D0,WORK(KCINT),1,WORK(KINTAI),1)
C        Calculate the term as L_{ai,em} * C'_{em}
         CALL CCG_LXD(WORK(KETAC),ISYRES,WORK(KINTAI),ISYMC,
     &                                   WORK(KL2AM),ISYML,1)
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),WORK(KETAC),1,WORK(KETAC),1)
            WRITE(LUPRI,1) 'Norm of alone Tbar_ck,ei*Cei: ',ETA1
         ENDIF
         !removed factor 2 to get FCI limit!
         CALL DAXPY(NT1AM(ISYRES),1.0D0,WORK(KETAC),1,ETAC(1),1)
      END IF
C
C---------------------------------
C Calculate the doubles (-)G term
C---------------------------------
C
      KOFFM = KOFFP + NT2AM(ISYRES)
      CALL CC_T2SQSYMSCAL(WORK(KL2AM),ISYML,ZERO)
C
      CALL CCRHS_E3(DUMMY,.FALSE.,WORK(KL2AM),WORK(KEI1),WORK(KEI2),
     &              WORK(KEND3),LEND3,ISYML,ISYMC,
     &              ETAC(KOFFM),.TRUE.)
C
C------------------------------------------------
C     Workspace for T2AM and X and Y is now free.
C     etac2 = P(ab,ij)(2l(ai)*Cjb - l(aj)*c(ib))
C------------------------------------------------
C
      IF (.NOT. CCS) THEN
C
         CALL CC_L1FCK3P(ETAC(1+NT1AM(ISYRES)),WORK(KL1AM),WORK(KINTIA),
     *                   ISYML,ISYMC)
C
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            ETA2 = DDOT(2*NT2AM(ISYRES),ETAC(1+NT1AM(ISYRES)),1,
     *                  ETAC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1-after L1c cont:       ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2-after L1c cont:       ',ETA2
         ENDIF
      ENDIF
C
C     Permute indices of plus vector
      CALL CCRHS3_IJ(ETAC(KOFFP),WORK(KEND2),LEND2,ISYRES)
C     For now, explicitly zero (-) diagonal!
      IF (ISYRES.EQ.1) CALL CCLR_DIASCL(ETAC(KOFFM),ZERO,ISYRES)
C
      IF (IPRINT .GT. 5 ) THEN
         TIMEC = SECOND() - TIMEC
         WRITE(LUPRI,9999) 'CCCI_ETA^C      ', TIMEC
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
9999  FORMAT(1x,'Time used in',2x,A18,2x,': ',f10.2,' seconds')
C
      END
C
      SUBROUTINE CC_XIC13(ISYMC,LBLC,        !Operator stuff
     &                    XIC,               !Result vector
     &                    LIST,ILSTNR,       !Left vector
     &                    LEOM,
     &                    XINT,WORK,LWORK)
C
C-----------------------------------------------------------------------
C     Purpose: Calculate the xiC(R) vector, when R is a triplet vector
C     and C is a singlet operator.
C     This means that the result vector will be a triplet vector.
C     (Modified version of the CCSD etaC(l0/l1) code
C
C     Input:
C
C     ISYMC   Symmetry of C
C     LBLC    Label of C integrals
C     XIC     Result vector (output)
C     LIST    List type of R vector
C     ILSTNR  List number of R vector
C     LEOM    Whether to include EOM disconnected terms (logical)
C     XINT    Can be used to pass C integrals in memory
C     WORK
C     LWORK
C
C     R. Faber 2017
C
C-----------------------------------------------------------------------
C
      implicit none

      character(len=*), parameter :: myname = 'CC_XI13'
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccexci.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
!Sonia & Filip: find a better place
#include "ccsections.h"
#include "second.h"

      DOUBLE PRECISION, PARAMETER :: TWO  = 2.0D0, HALF = 0.5D0,
     &                               ZERO = 0.0D0, ONEM = -1.0D0
      LOGICAL, PARAMETER :: LOCDBG = .false.
C
C     Input:
C
      INTEGER, INTENT(IN) :: ISYMC,  !  Symmetry of C
     &                       ILSTNR, ! List number for input
     &                       LWORK   !
      CHARACTER(LEN=*), INTENT(IN) ::  LBLC, !
     &                                 LIST  ! See above
      DOUBLE PRECISION, INTENT(OUT)::  XIC(*), ! Output array
     &                                 WORK(LWORK)  ! Work array
      DOUBLE PRECISION, INTENT(INOUT) :: XINT(*)
      LOGICAL, INTENT(IN) :: LEOM ! Whether to include extra EOM-CC
                                  ! Contributions

      CHARACTER MODEL*10
      LOGICAL :: FCKCON, ETRAN
      INTEGER :: KT1AM, KT2AM, KL1AM, KL2AM, KCTMO, KLAMDP, KLAMDH,
     &           KEI1, KEI2, KEND1, KEND2, KEND21, KEND3,
     &           KINTAI, KINTIA, KXMAT, KYMAT, KETAC, KCINT
      INTEGER :: KOFF1, KOFF2, KOFFP, KOFFM
      INTEGER :: ISYMR, ISYRES, ISYMA, ISYMI
      INTEGER :: LEND1, LEND2, LEND21, LEND3
      INTEGER :: IOPT
      DOUBLE PRECISION :: FF, XXI, XYI, ETA1, ETA2
      DOUBLE PRECISION :: TIMEC
C     EXTERNAL FUNCTIONS
      DOUBLE PRECISION :: DDOT
      INTEGER :: ILSTSYM
C
      CALL AROUND( 'Constructing XI^{'// LBLC
     &              //'}('// LIST //') vector ')
C      IF ( (IPRINT .GT. 10).or.locdbg ) THEN
C         CALL AROUND( 'IN CCCI_ETAC: Constructing ETAC(LE) vector ')
C      ENDIF
C
C--------------------------------
C     find symmetry of D operator.
C--------------------------------
C
      ISYMR = ILSTSYM(LIST,ILSTNR)
C
      ISYRES = MULD2H(ISYMR,ISYMC)
      IF (( LIST .EQ. 'R0')) THEN
         CALL QUIT('Misuse of '//myname)
      ENDIF
C
      TIMEC = SECOND()
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
C--------------------
C     Allocate space.
C--------------------
C
      KCTMO  = 1
      KT1AM  = KCTMO  + N2BST(ISYMC)
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
C
      LEND1  = LWORK  - KEND1
      IF ( .NOT. CCS) THEN

         KINTIA = KEND1
         KEND1  = KINTIA + NT1AM(ISYMC)
         LEND1  = LWORK - KEND1
         CALL DZERO(WORK(KintIA),NT1AM(isymc))
C
         KL1AM = KEND1
         KL2AM = KL1AM + NT1AM(ISYMR)
         KEND2 = KL2AM + NT2SQ(ISYMR)
         LEND2 = LWORK - KEND2
         KT2AM = KEND2
         KEND21= KT2AM + MAX(NT2AM(1),NT2AM(ISYMR))
         LEND21= LWORK - KEND21
C
      ELSE
C
         KL1AM = KEND1
         KEND2 = KL1AM + NT1AM(ISYMR)
         LEND2 = LEND1
         KEND21= KEND1
         LEND21= LEND1
C
      ENDIF
      KEI1   = KEND21
      KEI2   = KEI1   + NEMAT1(ISYMC)
      KEND3  = KEI2   + NMATIJ(ISYMC)
      LEND3  = LWORK  - KEND3
      IF (LEND3.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN '//myname)
      IF (LEND21.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN '//myname)
C
      IOPT = 3
C      CALL CC_RDRSP(LIST,ILSTNR,ISYMR,IOPT,MODEL,
C     *              ETAC,ETAC(1+NT1AM(ISYMR)))
C      ff = ddot(NT1AM(isyml)+2*nt2am(isyml),etac,1,etac,1)
C      write(lupri,*) 'Norm of total', ff

C
C-----------------------
C     get T1 amplitudes.
C-----------------------
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      IF ( .NOT. CCS) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND21),LEND21)
C
C-------------------------------
C     get AO property integrals.
C-------------------------------
C
      CALL DZERO(WORK(KCTMO),N2BST(ISYMC))
      FF = 1.0D0
C SLV98,OC give integrals option
      IF (LBLC.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMC),XINT(1),1,WORK(KCTMO),1)
      ELSE
        FF = 1.0D0
        CALL CC_ONEP(WORK(KCTMO),WORK(KEND21),LEND21,FF,ISYMC,LBLC)
      ENDIF
C
C-----------------------------------------------
C     Make MO T1-transformed property integrals.
C-----------------------------------------------
C
      CALL CC_FCKMO(WORK(KCTMO),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND21),LEND21,ISYMC,1,1)
C
C----------------------------------------------------------
C     Extract Cia (stored ia) and reorder ai
C----------------------------------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMC)
C
         DO 110 A = 1,NVIR(ISYMA)
C
            DO 120 I = 1,NRHF(ISYMI)
C
               KOFF1 = KINTIA + IT1AM(ISYMA,ISYMI)
     &               + NVIR(ISYMA)*(I - 1) + A-1
               KOFF2 = KCTMO + IFCVIR(ISYMI,ISYMA)
     *               + NORB(ISYMI)*(A - 1) + I - 1
C
               WORK(KOFF1) = WORK(KOFF2)
C
  120       CONTINUE
  110    CONTINUE
C
  100 CONTINUE
C
C     Initialize XIC
      CALL DZERO(XIC,NT1AM(ISYRES)+2*NT2AM(ISYRES))

C----------------------------------------------
C     Read R1, and R2(+) multipliers/left vectors.
C----------------------------------------------
C
      IOPT = 65
      CALL CC_RDRSP(LIST,ILSTNR,ISYMR,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KT2AM))
C
C--------------------------------
C     Put C into E matrix format.
C--------------------------------
C
      FCKCON = .TRUE.
      ETRAN  = .FALSE.
      CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                WORK(KCTMO),WORK(KEND3),LEND3,FCKCON,
     *                ETRAN,ISYMC)
C                   ~
C Square L2(+) onto L2
      IF (.NOT. CCS) THEN
         CALL CCRHS3_R2IJ(WORK(KT2AM),WORK(KEND3),LEND3,ISYMR)
         CALL CC_T2SQ(WORK(KT2AM),WORK(KL2AM),ISYMR)
C---------------------------------
C Calculate the doubles (+)G term
C---------------------------------
C Scale C by 1/2, as this term has a factor of 1/2
         CALL DSCAL(NMATAB(ISYMC),HALF,WORK(KEI1),1)
         CALL DSCAL(NMATIJ(ISYMC),HALF,WORK(KEI2),1)
C Calculate actual term
         KOFFP = NT1AM(ISYRES) + 1
         CALL CCRHS_E(XIC(KOFFP),WORK(KL2AM),WORK(KEI1),WORK(KEI2),
     &                WORK(KEND3),LEND3,ISYMR,ISYMC)
C Restore factor of C
         CALL DSCAL(NMATAB(ISYMC),TWO,WORK(KEI1),1)
         CALL DSCAL(NMATIJ(ISYMC),TWO,WORK(KEI2),1)
C
C----------------------------------------------
C     Read R2(-) multipliers/left vector.
C----------------------------------------------
C     The routine needs to be modified to allow this
         IOPT = 128
         CALL CC_RDRSP(LIST,ILSTNR,ISYMR,IOPT,MODEL,
     *                 WORK(KL1AM),WORK(KT2AM))
C
C We need R(+) + R(-) in the following.
         CALL CC_T2SQ3A(WORK(KT2AM),WORK(KL2AM),ISYMR,1.0D0)
C
C        Doubles contribution to Xi1:
C        Xi1 = sum_{em} C_{me} ((+)R+(-)R)_{ai,em}
         CALL CCG_LXD(XIC,ISYRES,WORK(KINTIA),ISYMC,
     &                WORK(KL2AM),ISYMR,1)
C
C---------------------------------
C Calculate the doubles (-)G term
C---------------------------------
C
         KOFFM = KOFFP + NT2AM(ISYRES)
         CALL CC_T2SQSYMSCAL(WORK(KL2AM),ISYMR,ZERO)
         CALL CCRHS_E3(DUMMY,.FALSE.,WORK(KL2AM),WORK(KEI1),WORK(KEI2),
     &                 WORK(KEND3),LEND3,ISYMR,ISYMC,
     &                 XIC(KOFFM),.TRUE.)
C
      END IF
C
C--------------------------------------------
C     xi1 =  sum(b)Lbi*Cba - sum(j)Laj*Cij.
C--------------------------------------------
C
      CALL CCLR_E1C1(XIC,WORK(KL1AM),WORK(KEI1),WORK(KEI2),
     *               WORK(KEND3),LEND3,ISYMR,ISYMC,'N')
C
C-----------------------------------------
C     Put T2 (packed) amplitudes in KT2AM.
C-----------------------------------------
C
      IF (.NOT. CCS) THEN
         !read in T2AM (packed)
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
      ENDIF
C
      IF (LEOM) THEN
C
C---------------------------------------------------------------
C     EOM contribution to XIC:      ~
C     (+/-) XiC2 = R_{ai}*(C_{bj} + t_{bj,fn}*C_{nf})
C---------------------------------------------------------------
         KCINT = KEND21
         KINTAI = KCINT + NT1AM(ISYMC)
C        Extract C_{ai}
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMC)
            DO I = 1, NRHF(ISYMI)
               KOFF1 = KINTAI + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)
               KOFF2 = KCTMO + IFCRHF(ISYMA,ISYMI) +
     &                 NORB(ISYMA)*(I-1) + NRHF(ISYMA)
               CALL DCOPY(NVIR(ISYMA),WORK(KOFF2),1,WORK(KOFF1),1)
            END DO
         END DO
         !get ttilde (overwrites T2am)
         CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,1,1)
C            ~
C        Add t_{em,fn} *C_{nf} to C_{em}
         CALL CCG_LXD(WORK(KCINT),ISYMC,WORK(KINTIA),ISYMC,
     &                WORK(KT2AM),1,0)
         CALL DAXPY(NT1AM(ISYMC),1.D0,WORK(KCINT),1,WORK(KINTAI),1)
C        Calculate the term as R_{ai} * C'_{bj}
         CALL CC_L1FCK3P(XIC(1+NT1AM(ISYRES)),WORK(KL1AM),WORK(KINTAI),
     *                   ISYMR,ISYMC)
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),WORK(KETAC),1,WORK(KETAC),1)
            WRITE(LUPRI,1) 'Norm of alone Tbar_ck,ei*Cei: ',ETA1
         ENDIF
      END IF
C
C------------------------------------------------
C     Workspace for T2AM and X and Y is now free.
C     etac2 = P(ab,ij)(2l(ai)*Cjb - l(aj)*c(ib))
C------------------------------------------------
C
      IF (.NOT. CCS) THEN
C
C
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),XIC(1),1,XIC(1),1)
            ETA2 = DDOT(2*NT2AM(ISYRES),XIC(1+NT1AM(ISYRES)),1,
     *                  XIC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1-after L1c cont:       ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2-after L1c cont:       ',ETA2
         ENDIF
      ENDIF
C
C     Permute indices of plus vector
      CALL CCRHS3_IJ(XIC(KOFFP),WORK(KEND2),LEND2,ISYRES)
C     For now, explicitly zero (-) diagonal!
      CALL CCLR_DIASCL(XIC(KOFFM),ZERO,ISYRES)


      IF (IPRINT .GT. 5 ) THEN
         TIMEC = SECOND() - TIMEC
         WRITE(LUPRI,9999) 'CCCI_ETA^C      ', TIMEC
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
9999  FORMAT(1x,'Time used in',2x,A18,2x,': ',f10.2,' seconds')
C
      END

