*=====================================================================*
      SUBROUTINE CC_RHTR_NODDY(IVEC,LUFC2,FC2AM,FREQC,
     &                         C1AM,FOCKC,ISYMC1,
     &                         FOCK0,XLAMDP0,XLAMDH0,
     &                         OMEGA1,OMEGA2,
     &                         WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to Jacobian transformation
*
*    Written by Christof Haettig, Februar 2003. 
*
*=====================================================================*
      IMPLICIT NONE  
#include <dummy.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0=1)

      CHARACTER FC2AM(*)
      INTEGER LWORK, ISYMC1, LUFC2, IVEC

#if defined (SYS_CRAY)
      REAL FREQC, DDOT, ONE, TWO, FF
      REAL WORK(LWORK), FOCK0(*), FOCKC(*)
      REAL XLAMDP0(*), XLAMDH0(*), C1AM(*)
      REAL OMEGA1(*), OMEGA2(*)
#else
      DOUBLE PRECISION FREQC, DDOT, ONE, TWO, FF
      DOUBLE PRECISION WORK(LWORK), FOCK0(*), FOCKC(*)
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), C1AM(*)
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*)
#endif
      PARAMETER ( ONE = 1.0D0, TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KSCR1, KFOCKD, KFOCK0, KFIELD, KFOCKC, KLAMPC, KLAMHC,
     &        KINT1S0, KINT2S0, KINT1T0, KINT2T0, KXIAJB, KYIAJB,
     &        KINT1SC, KINT2SC, KINT1TC, KINT2TC, KFLDC1, KFIELDAO,
     &        KT03AM, KT02AM, LUFOCK, LUSIFC,
     &        KEND1, LWRK1, KTC3AM, IF, IOPT, IJ, NIJ,
     &        KEND1A, LWRK1A, KTC2AM

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCRTRNOD')

      IF(DIRECT)CALL QUIT('DIRECT NOT IMPLEMENTED IN CC_RHTR_NODDY')

      IF (ISYMC1.NE.1 .OR. NSYM.NE.1)
     &  CALL QUIT('No symmetry yet implemented in CC_RHTR_NODDY!')

      IF (LOCDBG) THEN
        WRITE(LUPRI,*), 'CC_RHTR_NODDY> RESULT VECTOR ON INPUT:'
        WRITE(LUPRI,*), 'NORM^2 = ',
     &         DDOT(NT1AMX+NT2AMX,OMEGA1,1,OMEGA1,1)
        CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KEND1   = 1

      KFOCKD  = KEND1 
      KFOCK0  = KFOCKD + NORBT
      KFIELD  = KFOCK0 + NORBT*NORBT
      KEND1   = KFIELD + NORBT*NORBT

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB + NT1AMX*NT1AMX
      KEND1   = KYIAJB + NT1AMX*NT1AMX

      KTC3AM  = KEND1
      KEND1   = KTC3AM + NT1AMX*NT1AMX*NT1AMX

      ! what is above has to be kept until the end...
      ! everything below might be overwritten in CC_RHPART_NODDY
      KEND1A  = KEND1
      LWRK1A  = LWORK  - KEND1A
 
      KFOCKC  = KEND1
      KLAMPC  = KFOCKC + NORBT*NORBT
      KLAMHC  = KLAMPC + NLAMDT
      KFLDC1  = KLAMHC + NLAMDT
      KEND1   = KFLDC1 + NORBT*NORBT

      KINT1S0 = KEND1
      KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KINT1SC = KEND1
      KINT2SC = KINT1SC + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2SC + NRHFT*NRHFT*NT1AMX

      KINT1TC = KEND1
      KINT2TC = KINT1TC + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2TC + NRHFT*NRHFT*NT1AMX

      KT03AM  = KEND1
      KEND1   = KT03AM + NT1AMX*NT1AMX*NT1AMX

      KT02AM  = KEND1
      KTC2AM  = KT02AM + NT1AMX*NT1AMX*NT1AMX
      KEND1   = KTC2AM + NT1AMX*NT1AMX*NT1AMX

      KSCR1   = KEND1
      KEND1   = KSCR1  + NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_RHTR_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = 0
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = 0
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),XLAMDP0,XLAMDH0,
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*     Get response lambda matrices:
*---------------------------------------------------------------------*
      CALL CCLR_LAMTRA(XLAMDP0,WORK(KLAMPC),XLAMDH0,WORK(KLAMHC),
     &                 C1AM,ISYMC1) 

*---------------------------------------------------------------------*
*     Get the one electron integrals from the fields
*---------------------------------------------------------------------*
      IF (NONHF .AND. NFIELD .GT. 0) THEN
         KFIELDAO = KEND1
         KEND1A   = KFIELDAO + NORBT*NORBT
         LWRK1A   = LWORK  - KEND1A
         IF (LWRK1A .LT. 0) THEN
            CALL QUIT('Insufficient space in CC_RHTR_NODDY')
         ENDIF
         
         CALL DZERO(WORK(KFIELDAO),N2BST(ISYM0))
         DO I = 1, NFIELD
            FF = EFIELD(I)
            CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1A),LWRK1A,
     *                   FF,1,LFIELD(I))
         ENDDO
         CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)

         ! calculate external field in zero-order lambda basis
         CALL CC_FCKMO(WORK(KFIELD),XLAMDP0,XLAMDH0,
     *                 WORK(KEND1A),LWRK1A,1,1,1)

         ! one-index transformed field integrals for [V,R1]
         CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFLDC1),1)
         CALL CC_FCKMO(WORK(KFLDC1),XLAMDP0,WORK(KLAMHC),
     *                 WORK(KEND1A),LWRK1A,1,1,1)

         CALL CC_FCKMO(WORK(KFIELDAO),WORK(KLAMPC),XLAMDH0,
     *                 WORK(KEND1A),LWRK1A,1,1,1)
         CALL DAXPY(NORBT*NORBT,ONE,WORK(KFIELDAO),1,WORK(KFLDC1),1)

      ENDIF

*---------------------------------------------------------------------*
*     Compute some integrals:
*           XINT1S0 =  (CK|BD)
*           XINT2S0 =  (CK|LJ)
*           XINT1T0 =  (KC|BD)
*           XINT2T0 =  (KC|LJ)
*           XIAJB   = 2(IA|JB) - (IB|JA)
*           YIAJB   =  (IA|JB)
*---------------------------------------------------------------------*
      CALL CCSDT_INTS0_NODDY(.TRUE.,WORK(KXIAJB), WORK(KYIAJB),
     &                       .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .TRUE.,WORK(KINT1T0),WORK(KINT2T0),
     &                       XLAMDP0,XLAMDH0,
     &                       WORK(KEND1),LWRK1) 

*---------------------------------------------------------------------*
*     Loop over distributions of integrals and compute first-order
*     response versions XINT1SC=(ck|db)-bar and XINT2SC=(ck|lj)-bar
*---------------------------------------------------------------------*
      CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SC),WORK(KINT2SC),
     &                       .TRUE.,WORK(KINT1TC),WORK(KINT2TC),
     &                       XLAMDP0,XLAMDH0,
     &                       WORK(KLAMPC),WORK(KLAMHC),
     &                       WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Get zero-order amplitudes:
*---------------------------------------------------------------------*
      IF (LWRK1 .LT. NT2AMX) 
     &   CALL QUIT('Insufficient space in CC_RHTR_NODDY')

      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT02AM),ISYM0)
 
      IF ((NONHF) .AND. NFIELD.GT.0)  THEN
        IF (LWRK1 .LT. NT1AMX*NT1AMX*NT1AMX) 
     &    CALL QUIT('Out of memory for finite diff. in CC_RHTR_NODDY')
      END IF

      CALL CCSDT_T03AM(WORK(KT03AM),WORK(KINT1S0),WORK(KINT2S0),
     &                 WORK(KT02AM),WORK(KSCR1),WORK(KFOCKD),
     &                 WORK(KFIELD),WORK(KEND1))

*---------------------------------------------------------------------*
*     Add contribution from <mu_2|[[H,R1],T^0_3]|HF>
*---------------------------------------------------------------------*
      IF (LWRK1 .LT. NT1AMX*NT1AMX) 
     &   CALL QUIT('Insufficient space in CC_RHTR_NODDY')

      CALL DZERO(WORK(KEND1),NT1AMX*NT1AMX)

      CALL CCSDT_OMEGA2(WORK(KEND1),WORK(KINT1TC),WORK(KINT2TC),
     &                  WORK(KT03AM),FOCKC)
 
      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KEND1+IJ-1)
         END DO
      END DO

*---------------------------------------------------------------------*
*     Compute contribution from <mu_3|[H,T^B_2]|HF>/eps_3:
*---------------------------------------------------------------------*
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CC_RHTR_NODDY')
      ENDIF

      CALL DZERO(WORK(KTC3AM),NT1AMX*NT1AMX*NT1AMX)
 
      ! get doubles component of trial vector from file
      CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMC1),NT2AM(ISYMC1),
     *             IVEC,WORK(KEND1))

      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMC1)
      CALL CC_T2SQ(WORK(KEND1),WORK(KTC2AM),ISYMC1)
 
      CALL CCSDT_T3AM_R(WORK(KTC3AM),FREQC,
     &                  WORK(KINT1S0),WORK(KINT2S0),WORK(KTC2AM),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)

*---------------------------------------------------------------------*
*     Compute contribution from <mu_3|[H^B,T^0_2]|HF>/eps_3:
*---------------------------------------------------------------------*
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CC_RHTR_NODDY')
      ENDIF
 
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT02AM),ISYM0)
 
      CALL CCSDT_T3AM_R(WORK(KTC3AM),FREQC,
     &                  WORK(KINT1SC),WORK(KINT2SC),WORK(KT02AM),
     &                  WORK(KSCR1),WORK(KFOCKD),
     &                  .FALSE.,DUMMY,.FALSE.)
 
*---------------------------------------------------------------------*
*     Add finite difference contributions:
*---------------------------------------------------------------------*
      IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
        ! add <mu_3|[[V,R1],T^0_3]|HF>
        CALL CCSDT_XKSI3_2(WORK(KTC3AM),WORK(KFLDC1),WORK(KT03AM))

        ! add <mu_3|[[V,R2],T^0_2]|HF>
        CALL CCSDT_XKSI3_1(WORK(KTC3AM),WORK(KFIELD),
     &                     WORK(KT02AM),WORK(KTC2AM),ONE)
        CALL CCSDT_XKSI3_1(WORK(KTC3AM),WORK(KFIELD),
     &                     WORK(KTC2AM),WORK(KT02AM),ONE)
      END IF

*---------------------------------------------------------------------*
*     solve triples equations and add contribution to the
*     singles doubles result vectors:
*---------------------------------------------------------------------*
      CALL CC_RHPART_NODDY(OMEGA1,OMEGA2,WORK(KTC3AM),FREQC,
     &                     WORK(KFOCKD),WORK(KFOCK0),WORK(KFIELD),
     &                     WORK(KXIAJB),WORK(KINT1T0),WORK(KINT2T0),
     &                     WORK(KEND1A),LWRK1A)

*---------------------------------------------------------------------*
*     Print debug output and return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE(LUPRI,*),'CC_RHTR_NODDY> RESULT VECTOR ON OUTPUT:'
        WRITE(LUPRI,*),'NORM^2(OMEGA1)=',DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
        WRITE(LUPRI,*),'NORM^2(OMEGA2)=',DDOT(NT2AMX,OMEGA2,1,OMEGA2,1)
        CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF

      CALL QEXIT('CCRTRNOD')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_RHTR_NODDY                        *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CC_RHPART_NODDY(OMEGA1,OMEGA2,R3AM,FREQ,
     &                           FOCKD,FOCK0,FIELD,
     &                           XIAJB,XINT1T,XINT2T,
     &                           WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: solve 'right' triples equations and partition the
*             triples solution vector into an effective rhs vector
*
*    Written by Christof Haettig, Februar 2003. 
*
*=====================================================================*
      IMPLICIT NONE  
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      INTEGER LWORK
#if defined (SYS_CRAY)
      REAL FREQ, FIELD(*)
      REAL OMEGA1(*), OMEGA2(*), R3AM(*), FOCKD(*), FOCK0(*)
      REAL WORK(*), XIAJB(*), XINT1T(*), XINT2T(*)
#else
      DOUBLE PRECISION FREQ, FIELD(*)
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), R3AM(*), FOCKD(*), FOCK0(*)
      DOUBLE PRECISION WORK(*), XIAJB(*), XINT1T(*), XINT2T(*)
#endif
 
      LOGICAL TRANSPOSE
      INTEGER KSCR1, KEND1, LWRK1, KT3SCR, KOMEGA1, KOMEGA2, IJ, NIJ

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCRPARTN')

*---------------------------------------------------------------------*
*     Solve triples equation:
*---------------------------------------------------------------------*
      KSCR1 = 1
      KEND1 = KSCR1 + NT1AMX
      IF (NONHF) THEN
         KT3SCR = KEND1
         KEND1  = KT3SCR + NT1AMX*NT1AMX*NT1AMX
      END IF

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_RHPART_NODDY (1)')
      ENDIF

      TRANSPOSE = .FALSE.
      CALL CCSDT_3AM(R3AM,FREQ,WORK(KSCR1),FOCKD,
     &               NONHF,FIELD,TRANSPOSE,WORK(KT3SCR))

*---------------------------------------------------------------------*
*     Add contribution to singles doubles result vectors:
*---------------------------------------------------------------------*
      KOMEGA1 = 1
      KOMEGA2 = KOMEGA1 + NT1AMX
      KEND1   = KOMEGA2 + NT1AMX*NT1AMX

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_RHPART_NODDY (2)')
      ENDIF

      CALL DZERO(WORK(KOMEGA2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOMEGA1),NT1AMX)

      CALL CCSDT_OMEGA1(WORK(KOMEGA1),XIAJB,R3AM)
 
      CALL CCSDT_OMEGA2(WORK(KOMEGA2),XINT1T,XINT2T,R3AM,FOCK0)
                                                         
      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOMEGA1+I-1)
      END DO
 
      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOMEGA2+IJ-1)
         END DO
      END DO

      CALL QEXIT('CCRPARTN')
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_RHPART_NODDY                      *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCDOTRSP_NODDY(OMEGA1,OMEGA2,OMEGA3,SIGN,
     &                          ITRAN,LISTDP,IDOTS,DOTPROD,MXVEC,
     &                          XLAMDP,XLAMDH,FOCK0,FOCKD,
     &                          XIAJB,YIAJB,XINT1T,XINT2T,
     &                          XINT1S,XINT2S,
     &                          MODUL,LDBG,LDBGSD,SKIP_T3,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to dot products
*             similar as CCDOTRSP does it for singles and doubles
*
*    Written by Christof Haettig, Februar 2003. 
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>

      LOGICAL LDBG, LDBGSD, SKIP_T3
      CHARACTER*(*) LISTDP, MODUL
      INTEGER LWORK, MXVEC, ITRAN
      INTEGER IDOTS(MXVEC,*)

#if defined (SYS_CRAY)
      REAL OMEGA1(*), OMEGA2(*), OMEGA3(*), WORK(*)
      REAL XLAMDP(*), XLAMDH(*), FOCK0(*), FOCKD(*)
      REAL XIAJB(*), XINT1T(*), XINT2T(*)
      REAL YIAJB(*), XINT1S(*), XINT2S(*)
      REAL DOTPROD(MXVEC,*), DDOT, FREQC, SCON, DCON, TCON
      REAL SIXTH, TWO, SIGN
#else
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), OMEGA3(*), WORK(*)
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), FOCK0(*), FOCKD(*)
      DOUBLE PRECISION XIAJB(*), XINT1T(*), XINT2T(*)
      DOUBLE PRECISION YIAJB(*), XINT1S(*), XINT2S(*)
      DOUBLE PRECISION DOTPROD(MXVEC,*), DDOT, FREQC, SCON, DCON, TCON
      DOUBLE PRECISION SIXTH, TWO, SIGN
#endif
      PARAMETER ( SIXTH=1.0D0/6.0D0, TWO=2.0D0 )
 
      CHARACTER MODEL*(10)
      INTEGER INDEX, IVEC, IDLSTC, ISYMC, KSCR1, KLC3AM, KEND1, KLC1AM,
     &        KLC2AM, LWRK1, ILSTSYM, IOPT, KLAMPC, KLAMHC, KFOCKC,
     &        KINT1SC, KINT2SC, KEND2, LWRK2
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCDOTRSP')

*---------------------------------------------------------------------*
*     say hallo:
*---------------------------------------------------------------------*
      IF (LDBG) THEN
        WRITE(LUPRI,*) 'CCDOTRSP> call from ',MODUL
        WRITE(LUPRI,*) 'LISTDP = ',LISTDP
        WRITE(LUPRI,*) 'ITRAN, MAXVEC = ',ITRAN,MXVEC
        IF (LDBGSD) THEN
          WRITE(LUPRI,*) 'NORM^2(OMEGA1):',
     &        DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
          WRITE(LUPRI,*) 'NORM^2(OMEGA2):',
     &        DDOT(NT1AMX**2,OMEGA2,1,OMEGA2,1)
        END IF
        WRITE(LUPRI,*) 'NORM^2(OMEGA3):',
     &      DDOT(NT1AMX**3,OMEGA3,1,OMEGA3,1)
      END IF

*---------------------------------------------------------------------*
*     Memory allocation
*---------------------------------------------------------------------*
      KSCR1  = 1
      KLC3AM = KSCR1  + NT1AMX
      KEND1  = KLC3AM + NT1AMX*NT1AMX*NT1AMX

      IF (LDBGSD) THEN
        KLC1AM = KEND1
        KLC2AM = KLC1AM + NT1AMX
        KEND1  = KLC2AM + NT1AMX*NT1AMX
      END IF

      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCDOTRSP_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     loop over all vectors we need to dot on and compute the
*     vectors and then the triples contributions:
*---------------------------------------------------------------------*
      IVEC = 1
      DO WHILE (IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

        IDLSTC = IDOTS(IVEC,ITRAN)
        ISYMC  = ILSTSYM(LISTDP,IDLSTC)

        IF (.NOT.SKIP_T3) THEN
         IF     (LISTDP(1:3).EQ.'R1 '.OR.LISTDP(1:3).EQ.'RE '    ) THEN
           KLAMPC  = KEND1
           KLAMHC  = KLAMPC  + NLAMDT
           KFOCKC  = KLAMHC  + NLAMDT
           KINT1SC = KFOCKC  + NORBT*NORBT
           KINT2SC = KINT1SC + NT1AMX*NVIRT*NVIRT
           KEND2   = KINT2SC + NT1AMX*NRHFT*NRHFT
           LWRK2   = LWORK - KEND2

           IF (LWRK2 .LT. 0) THEN
             CALL QUIT('Insufficient space in CCDOTRSP_NODDY (T31)')
           ENDIF

           CALL CCSDT_T31_NODDY(WORK(KLC3AM),LISTDP,IDLSTC,FREQC,
     &                        .FALSE.,
     &                        .FALSE.,XINT1S,XINT2S,
     &                        .FALSE.,XINT1T,XINT2T,
     &                        .FALSE.,XIAJB, YIAJB,
     &                        WORK(KINT1SC),WORK(KINT2SC),
     &                        WORK(KLAMPC),WORK(KLAMHC),WORK(KFOCKC),
     &                        XLAMDP,XLAMDH,FOCK0,DUMMY,FOCKD,
     &                        WORK(KEND2),LWRK2)
 
         ELSEIF (LISTDP(1:3).EQ.'L1 '.OR.LISTDP(1:3).EQ.'LE '.OR.
     &           LISTDP(1:3).EQ.'M1 '.OR.LISTDP(1:3).EQ.'E0 '.OR.
     &           LISTDP(1:3).EQ.'N2 '                            ) THEN

           CALL CCSDT_TBAR31_NODDY(WORK(KLC3AM),FREQC,LISTDP,IDLSTC,
     &                             XLAMDP,XLAMDH,
     &                             FOCK0,FOCKD,WORK(KSCR1),
     &                             .FALSE.,XIAJB,DUMMY,
     &                             .FALSE.,XINT1T,XINT2T,
     &                             WORK(KEND1),LWRK1)
         ELSE
           CALL QUIT('Unknown or illegal list in CCDOTRSP_NODDY.')
         END IF
        END IF ! SKIP_T3

        IF (LDBG) THEN
          WRITE(LUPRI,*) 'CCDOTRSP_NODDY> IVEC,ITRAN:',IVEC,ITRAN
          WRITE(LUPRI,*) 'NORM^2(triples):',
     &      DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KLC3AM),1,WORK(KLC3AM),1)
        END IF

        IF (SKIP_T3) THEN
           TCON = 0.0D0
        ELSE
          TCON = SIGN * SIXTH *
     &      DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KLC3AM),1,OMEGA3,1)
        END IF

        DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TCON

        IF (LDBG) THEN
          IF (LDBGSD) THEN
            IOPT = 3
            Call CC_RDRSP(LISTDP,IDLSTC,ISYMC,IOPT,MODEL,
     &                    WORK(KLC1AM),WORK(KLC3AM))
            IF (LISTDP(1:3).EQ.'R1 '.OR.LISTDP(1:3).EQ.'RE '    ) THEN
               ! for right response vectors scale diagonal with 2
               Call CCLR_DIASCL(WORK(KLC3AM),TWO,ISYMC)
            ELSEIF (LISTDP(1:3).EQ.'L1 '.OR.LISTDP(1:3).EQ.'LE '.OR.
     &              LISTDP(1:3).EQ.'M1 '.OR.LISTDP(1:3).EQ.'E0 '.OR.
     &              LISTDP(1:3).EQ.'N2 '                        ) THEN
               ! no scaling for left response vectors
               CONTINUE
            ELSE
               ! print warning that debug output below might be wrong
               WRITE(LUPRI,*) 'do know how to scale diagonal...'
               WRITE(LUPRI,*) 'warning: DCON might be wrong!!!'
            END IF
            CALL CC_T2SQ(WORK(KLC3AM),WORK(KLC2AM),ISYMC)
            DCON = 0.5D0*DDOT(NT1AMX*NT1AMX,WORK(KLC2AM),1,OMEGA2,1)
            SCON = DDOT(NT1AMX,WORK(KLC1AM),1,OMEGA1,1) 
            WRITE(LUPRI,*) 'triples in singles vector:',SCON
            WRITE(LUPRI,*) 'triples in doubles vector:',DCON
          END IF
          WRITE(LUPRI,*) 'triples in triples vector:',TCON
        END IF

        IVEC = IVEC + 1

      END DO

      CALL QEXIT('CCDOTRSP')
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCDOTRSP_NODDY                       *
*---------------------------------------------------------------------*
