!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
*=====================================================================*
      SUBROUTINE CC_FDFBMAT(CMOQ,LISTL,IDLSTL,LISTR,IDLSTR,
     &                      RHODIF,WORK,LWORK)
C
C---------------------------------------------------------------------
C Test routine for calculating the orbital relaxation contribution
C to the CC F^B T^A transformation by finite difference on the usual
C F matrix transformation wrt the CMO coefficients.
C We pass the derivative C^(1) MO = CMOQ
C
C S. Coriani & Ch. Haettig, march 1999
C---------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
#include <ccinftap.h>
#include <dummy.h>
C
C------------------------------------------------------------
C     the displacement for the finite difference calculation:
C------------------------------------------------------------
      PARAMETER (DELTA = 1.0D-07, DELINV = 1.0D+7)
C------------------------------------------------------------
      DIMENSION WORK(LWORK), RHODIF(*), CMOQ(*)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, TWO = 2.0d0, XHALF = 0.5D0)
      CHARACTER MODEL*10, LISTL*(*), LISTR*(*)
      LOGICAL   LHTF
C
      IF (IPRINT.GT.10) THEN
        CALL AROUND('IN CC_FDFB : MAKING FIN. DIFF. CC F*T^A Vector')
      END IF
C
      IF (CCR12) CALL QUIT('Finite-difference of F*T^A Vector for '//
     &                     'CCR12 not adapted')
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      ISYMOP     = 1
C
      NTAMP      = NT1AMX + NT2AMX
      IF (CCS) NTAMP = NT1AMX
C
      KCMOREF    = 1
      KCMO       = KCMOREF  + NLAMDS
      KRHREF1    = KCMO     + NLAMDS
      KRHREF2    = KRHREF1  + NT1AMX
      KEND0      = KRHREF2  + NT2AMX
      LWRK0      = LWORK    - KEND0
C
      KT1AMP0    = KEND0                                 
      KOMEGA1    = KT1AMP0  + NT1AMX                     
      KOMEGA2    = KOMEGA1  + NT1AMX
      KT2AMP0    = KOMEGA2  + MAX(NT2AMX,2*NT2ORT(1),NT2AO(1))  
      KSCR2      = KT2AMP0  + NT2AMX
      KEND1      = KSCR2    + NT2AMX + NT1AMX
      LWRK1      = LWORK    - KEND1
C
      KRHO1      = KEND0
      KRHO2      = KRHO1    + NT1AMX
      KEND1A     = KRHO2    + NT2AMX
      LWRK1A     = LWORK    - KEND1A
C
      IF ( LWRK1.LT.0 .OR. LWRK1A.LT.0) THEN                   
         WRITE(LUPRI,*) 'Too little work space in CC_FDFBMAT '
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',MAX(KEND1,KEND1A)
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDFBMAT ')
      ENDIF
C
C---------------------
C     Initializations.
C---------------------
C
      CALL DZERO(RHODIF,NTAMP)
C
      IPRSAVE = IPRINT
C
C---------------------------------------------
C     Read the reference CMO vector from disk:
C---------------------------------------------
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (WORK(KCMOREF-1+I),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP') 
C
C-------------------------------------------------------------------
C     make sure that the correct response intermediates are on disc:
C-------------------------------------------------------------------
C
      CALL DZERO(WORK(KT1AMP0),NT1AMX)
      LHTF = .FALSE.
      CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),LHTF,
     &               .FALSE.,.FALSE.,WORK(KEND1),LWRK1)
      REWIND(LUIAJB)
      CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
      IOPT = 3
      CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2), 
     &            WORK(KT1AMP0),WORK(KT2AMP0),
     &            WORK(KEND1),LWRK1,'XXX')  
C
C-------------------------------------------------------------------
C     calculate the reference vector:
C-------------------------------------------------------------------
C use short-cut routine (for one single transformation)
C
      IPRINT = 0
      CALL CC_FTRAN(LISTL, IDLSTL, LISTR, IDLSTR, 
     &              WORK(KRHO1), LWRK0)
C
      IF (CCS) CALL DZERO(WORK(KRHO2),NT2AMX)
C
      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      IF (.NOT.CCS) RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)
C
      IF (IPRSAVE.GT.10) THEN
         WRITE (LUPRI,*) 'CC_FDFB: norm of reference RHO1:',RHO1N
         IF (.NOT.CCS) 
     *      WRITE (LUPRI,*) 'CC_FDFB: norm of reference RHO2:',RHO2N
      END IF
C
c  save result on Reference vector
C
      CALL DCOPY(NT1AMX,WORK(KRHO1),1,WORK(KRHREF1),1)          
      IF (.NOT.CCS) CALL DCOPY(NT2AMX,WORK(KRHO2),1,WORK(KRHREF2),1)
 
*----------------------------------------------------------------------*
*     Calculate the derivative of the vector function with respect
*     to the CMO vector by finite difference:
*----------------------------------------------------------------------*
 
      DO IDXDIF = 1, NLAMDS              
C
C        -------------------------------
C        add finite displacement to CMO:
C        -------------------------------
C
         CALL DCOPY(NLAMDS,WORK(KCMOREF),1,WORK(KCMO),1)
         WORK(KCMO-1 + IDXDIF) = WORK(KCMO-1 + IDXDIF) + DELTA  
C
         CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ(LUSIFC)
         READ(LUSIFC)
         WRITE(LUSIFC) (WORK(KCMO-1+I),I=1,NLAMDS)
         CALL GPCLOSE(LUSIFC,'KEEP')
C
C        -------------------------------------
C        calculate new response (global) intermediates:
C        -------------------------------------
C
         CALL DZERO(WORK(KT1AMP0),NT1AMX)                      
         LHTF = .FALSE.
         CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),
     &                  LHTF,.FALSE.,.FALSE.,WORK(KEND1),LWRK1)
         REWIND(LUIAJB)
         CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
         IOPT = 3
         CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))
C
         RSPIM = .TRUE.
         CALL CCRHSN(WORK(KRHO1),WORK(KRHO2),
     &               WORK(KT1AMP0),WORK(KT2AMP0),
     &               WORK(KEND1),LWRK1,'XXX')  
C
C        ---------------------------------
C        calculate the transformed vector:
C        ---------------------------------
C
         IPRINT = 0
         CALL CC_FTRAN(LISTL, IDLSTL, LISTR, IDLSTR,
     &                             WORK(KRHO1), LWRK0)
C
         IF (CCS) CALL DZERO(WORK(KRHO2),NT2AMX)
C
C        --------------------------------------------------
C        construct the row nb. IDXDIF of the result matrix:
C        --------------------------------------------------
C
         RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
         IF (.NOT.CCS) RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)
C
          IF (IPRSAVE.GT.10) THEN
            WRITE (LUPRI,*) 'CMO index:',IDXDIF
            WRITE (LUPRI,*) 'Norm of RHO1: ',RHO1N
            WRITE (LUPRI,*) 'Norm of RHO2: ',RHO2N
          END IF
C
c  Construct {[RES(delta)-RES(ref)]/delta}*C(1)
C
         CALL DAXPY(NT1AMX,-1.0D0,WORK(KRHREF1),1,WORK(KRHO1),1)       !
         IF (.NOT.CCS) 
     &     CALL DAXPY(NT2AMX,-1.0D0,WORK(KRHREF2),1,WORK(KRHO2),1)
         CALL DSCAL(NT1AMX,DELINV,WORK(KRHO1),1)
         IF (.NOT.CCS) 
     &     CALL DSCAL(NT2AMX,DELINV,WORK(KRHO2),1)
c
         CALL DAXPY(NT1AMX,CMOQ(IDXDIF),WORK(KRHO1),1,
     &                                     RHODIF(1),1)
         IF (.NOT.CCS) CALL DAXPY(NT2AMX,CMOQ(IDXDIF),WORK(KRHO2),1,
     &                                     RHODIF(1+NT1AMX),1)
C
         IF (IPRSAVE.GT.100) THEN
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),1,1,1)
            WRITE (LUPRI,*) 'CMOQ(index) = ', CMOQ(IDXDIF)
            WRITE (LUPRI,*) 'accumulated result vector:'
            CALL CC_PRP(RHODIF,RHODIF(NT1AMX+1),1,1,1)
         ENDIF
C
      END DO
C
C----------------------------------------------
C     restore the reference CMO vector on disk:
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      WRITE(LUSIFC) (WORK(KCMOREF-1+I),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C------------------------------------------------------------------
C     make sure that all intermediates on file are calculated with
C     the reference CMO vector:
C------------------------------------------------------------------
C
      CALL DZERO(WORK(KT1AMP0),NT1AMX)
      LHTF = .FALSE.
      CALL CCSD_IAJB(WORK(KT2AMP0),WORK(KT1AMP0),LHTF,
     &               .FALSE.,.FALSE.,WORK(KEND1),LWRK1)
      REWIND(LUIAJB)
      CALL WRITI(LUIAJB,IRAT*NT2AMX,WORK(KT2AMP0))
C
      IOPT = 3
      CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KT2AMP0))
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &            WORK(KT1AMP0),WORK(KT2AMP0),
     &            WORK(KEND1),LWRK1,'XXX')  
C
      IPRINT = IPRSAVE
C
      IF (IPRINT.GT.10) THEN
        CALL AROUND(' END OF CC_FDFBMAT ')
      END IF
C
      RETURN
      END
*=====================================================================*
      SUBROUTINE CC_FDFBMAT2(LISTR,IDLSTR,RHODIF1,RHODIF2,
     &                       LABEL,IRELAX,WORK,LWORK)
C
C---------------------------------------------------------------------
C Test routine for calculating a generalized CC F{O} matrix transformed 
C vector by finite difference on the Eta{O} vector
C
C Christof Haettig, march 1999
C---------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
#include <dummy.h>
#include <ccfro.h>
#include <ccroper.h>
#include <cclists.h>
C
C------------------------------------------------------------
C     the displacement for the finite difference calculation:
C------------------------------------------------------------
      PARAMETER (DELTA = 1.0D-07, DELINV = 1.0D+07)
C------------------------------------------------------------
C
      DIMENSION WORK(LWORK)
      DIMENSION RHODIF1(*), RHODIF2(*)
      CHARACTER*(3) LISTL, LISTR
      CHARACTER*(8) FILXI, FILETA
      CHARACTER MODEL*(10)
      INTEGER IXETRAN(MXDIM_XEVEC,3)
C     LOGICAL LORX 
C
      INTEGER IR1TAMP, IRHSR1, IETA1, IROPER, ILSTSYM
C
      PARAMETER (ONE=1.0d0, ZERO=0.0d0, TWO=2.0d0, HALF=0.5d0)
C
      IF (IPRINT.GT.5) THEN
         CALL AROUND('in CC_FDFBMAT2: making fini. diff. FB Matrix')
      ENDIF
C
C----------------------------------------------
C     set up IXETRAN array:
C----------------------------------------------
C
C     IRELAX = 0
C     IF (LORX) THEN
C       IRELAX = IR1TAMP(LABEL,LORX,FREQ,ISYHOP)
C     END IF

      ! set the IXETRAN array for one (XI,ETA) pair
      IXETRAN(1,1) = IROPER(LABEL,ISYHOP)
      IXETRAN(2,1) = 0
      IXETRAN(3,1) = -1
      IXETRAN(4,1) = 0
      IXETRAN(5,1) = IRELAX
C
      LISTL  = 'L0 '
      FILXI  = 'CCFDFBO1'
      FILETA = 'CCFDFBX1'
      IOPTRES = 0
      NXETRAN = 1
C
C----------------------------------------------
C     initializations & work space allocations.
C----------------------------------------------
C
      CALL DZERO(RHODIF1,NT1AMX)
      CALL DZERO(RHODIF2,NT2AMX)
C
      IPRSAV = IPRINT
      IPRINT = 0
C
      MODEL = 'UNKNOWN'
      IF (CCS)  MODEL = 'CCS'
      IF (CC2)  MODEL = 'CC2'
      IF (CCSD) MODEL = 'CCSD'
C
      KT1AMPSAV  = 1
      KT2AMPSAV  = KT1AMPSAV + NT1AMX
      KT1AMPA    = KT2AMPSAV + NT2AMX
      KT2AMPA    = KT1AMPA   + NT1AMX
      KEND0      = KT2AMPA   + NT2AMX
      LWRK0      = LWORK     - KEND0
C
      KT0AMP0    = KEND0
      KT1AMP0    = KT0AMP0   + 2*NALLAI(1)
      KOMEGA1    = KT1AMP0   + NT1AMX
      KOMEGA2    = KOMEGA1   + NT1AMX
      KT2AMP0    = KOMEGA2   + MAX(NT2AMX,2*NT2ORT(1),NT2AO(1))
      KSCR2      = KT2AMP0   + NT2AMX
      KEND1A     = KSCR2     + NT2AMX + NT1AMX
      LWRK1A     = LWORK     - KEND1A
C
      KRHO1      = KEND0
      KRHO2      = KRHO1     + NT1AMX
      KEND1B     = KRHO2     + NT2AMX
      LWRK1B     = LWORK     - KEND1B
C
      IF (LWRK1A.LT.0 .OR. LWRK1B.LT.0) THEN
         WRITE(LUPRI,*) 'TOO LITTLE WORK SPACE IN CC_FDFBMAT2::'
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',MAX(KEND1A,KEND1B)
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDFBMAT2: ')
      ENDIF

C     -------------------------------------------
C     Read the CC reference amplitudes from disk:
C     -------------------------------------------
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AMPSAV),
     *              WORK(KT2AMPSAV))
 
C     --------------------------------------------
C     Read the T^A reference amplitudes from disk:
C     --------------------------------------------
      IOPT  = 3
      ISYMR = ILSTSYM(LISTR,IDLSTR)
      CALL CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,
     &              WORK(KT1AMPA),WORK(KT2AMPA))
      IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,1)

*---------------------------------------------------------------------*
*     Add delta x t^A to cluster amplitudes and recalculate the 
*     response intermediates and the Eta^B vector:
*               Eta^B{t^0 + delta x t^A} 
*---------------------------------------------------------------------*

C     -------------------------------------------------------------
C     add finite displadement to t^0 and recalculate intermediates:
C     -------------------------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     &  CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
      CALL DAXPY(NT1AMX,DELTA,WORK(KT1AMPA),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     & CALL DAXPY(NT2AMX,DELTA,WORK(KT2AMPA),1,WORK(KT2AMP0),1)
 
      IOPT = 3
      CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,WORK(KT0AMP0),
     &              WORK(KT1AMP0),WORK(KT2AMP0),WORK(KEND1A),LWRK1A)
 
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),WORK(KT1AMP0),
     *            WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
 
C     ---------------------------------
C     calculate the transformed vector:
C     ---------------------------------
      IORDER = 1
      CALL CC_XIETA(IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL,
     &              FILXI,  IDUM, RDUM,
     &              FILETA, IDUM, RDUM,
     &              .FALSE.,0, WORK(KEND0), LWRK0 )

C     IOPT   = 3
C     IVEC   = IXETRAN(4,1)
C     ISYETA = ISYHOP
C     CALL CC_RDRSP('X1 ',IVEC,ISYETA,IOPT,MODEL,
C    &              WORK(KRHO1),WORK(KRHO2))

      LEN       = NT1AMX + NT2AMX
      IADRF_ETA = IXETRAN(4,1)
      LUETA = -1
      CALL WOPEN2(LUETA,FILETA,64,0)
      CALL GETWA2(LUETA,FILETA,WORK(KRHO1),IADRF_ETA,LEN)
      CALL WCLOSE2(LUETA,FILETA,'KEEP')

      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      IF (.NOT.CCS) RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'Norm of RHO1(t^0 + delta x t^C): ',RHO1N
         WRITE (LUPRI,*) 'Norm of RHO2(t^0 + delta x t^C): ',RHO2N
      END IF
 
      ! divide by 2*delta and copy to result vector:
      CALL DSCAL(NT1AMX,HALF*DELINV,WORK(KRHO1),1)
      IF (.NOT.CCS) CALL DSCAL(NT2AMX,HALF*DELINV,WORK(KRHO2),1)

      CALL DCOPY(NT1AMX,WORK(KRHO1),1,RHODIF1,1)
      IF (.NOT.CCS) CALL DCOPY(NT2AMX,WORK(KRHO2),1,RHODIF2,1)
 
*---------------------------------------------------------------------*
*     Substract delta x t^C to cluster amplitudes and recalculate the 
*     response intermediates and the F matrix transformed T^B vector:
*               F{t^0 - delta x t^C} t^B 
*---------------------------------------------------------------------*

C     -------------------------------------------------------------
C     add finite displadement to t^0 and recalculate intermediates:
C     -------------------------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     &  CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
      CALL DAXPY(NT1AMX,-DELTA,WORK(KT1AMPA),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     &  CALL DAXPY(NT2AMX,-DELTA,WORK(KT2AMPA),1,WORK(KT2AMP0),1)
 
      IOPT = 3
      CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,WORK(KT0AMP0),
     &              WORK(KT1AMP0),WORK(KT2AMP0),WORK(KEND1A),LWRK1A)
 
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &     WORK(KT1AMP0),WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
 
C     ---------------------------------
C     calculate the transformed vector:
C     ---------------------------------
      IORDER = 1
      CALL CC_XIETA(IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL,
     &              FILXI,  IDUM, RDUM,
     &              FILETA, IDUM, RDUM,
     &              .FALSE.,0, WORK(KEND0), LWRK0 )

C     IOPT   = 3
C     IVEC   = IXETRAN(4,1)
C     ISYETA = ISYHOP
C     CALL CC_RDRSP('X1 ',IVEC,ISYETA,IOPT,MODEL,
C    &              WORK(KRHO1),WORK(KRHO2))

      LEN       = NT1AMX + NT2AMX
      IADRF_ETA = IXETRAN(4,1)
      CALL WOPEN2(LUETA,FILETA,64,0)
      CALL GETWA2(LUETA,FILETA,WORK(KRHO1),IADRF_ETA,LEN)
      CALL WCLOSE2(LUETA,FILETA,'DELETE')

      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      IF (.NOT.CCS) RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'Norm of RHO1(t^0 + delta x t^C): ',RHO1N
         WRITE (LUPRI,*) 'Norm of RHO2(t^0 + delta x t^C): ',RHO2N
      END IF
 
      ! divide by 2*delta and substract from final result:
      CALL DAXPY(NT1AMX,-HALF*DELINV,WORK(KRHO1),1,RHODIF1,1)
      IF (.NOT.CCS) 
     & CALL DAXPY(NT2AMX,-HALF*DELINV,WORK(KRHO2),1,RHODIF2,1)

*---------------------------------------------------------------------*
*     fix the scale factor of the diagonal, print some output and
*     restore t^0 amplitudes and response intermediates on file:
*---------------------------------------------------------------------*

C     -----------------------------------------------------------------
C     scale diagonal with 1/2: (only for right vectors --> A,B,C,D mat)
C     -----------------------------------------------------------------
C     CALL CCLR_DIASCL(RHODIF2,TWO,1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'RESULT VECTOR FROM CC_FDFBMAT2:'
         CALL CC_PRP(RHODIF1,RHODIF2,1,1,1)
      ENDIF
 
C     --------------------------------------------
C     Restore the CC reference amplitudes on disk:
C     --------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     &  CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
 
      IOPT = 3
      CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,WORK(KT0AMP0),WORK(KT1AMP0),
     &              WORK(KT2AMP0),WORK(KEND1A),LWRK1A)
 
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &     WORK(KT1AMP0),WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
 
      IF (IPRSAV .GT. 5) THEN
         CALL AROUND(' END OF CC_FDFBMAT2:')
      ENDIF
 
      IPRINT = IPRSAV
 
      RETURN
      END
*=====================================================================*
*                      END OF SUBROUTINE CC_FDFBMAT2:                 *
*=====================================================================*

*=====================================================================*
      SUBROUTINE CC_FDFBMAT3(LISTR,IDLSTR,RHODIF1,RHODIF2,
     &                       LABEL,IRELAX,WORK,LWORK)
C
C---------------------------------------------------------------------
C Test routine for calculating the Eta{O} vector to test F{O}
C Sonia, dec 1999
C---------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
#include <ccfro.h>
#include <ccroper.h>
#include <cclists.h>
C
      DIMENSION WORK(LWORK)
      DIMENSION RHODIF1(*), RHODIF2(*)
      CHARACTER*(3) LISTL, LISTR
      CHARACTER*(8) FILXI, FILETA
      CHARACTER MODEL*(10)
      INTEGER IXETRAN(MXDIM_XEVEC,3)
C     LOGICAL LORX 
C
      INTEGER IR1TAMP, IRHSR1, IETA1, IROPER, ILSTSYM
C
      PARAMETER (ONE=1.0d0, ZERO=0.0d0, TWO=2.0d0, HALF=0.5d0)
      PARAMETER (DELTA = 1.0D-07, DELINV = 1.0D+7)
C
      IF (IPRINT.GT.5) THEN
         CALL AROUND('in CC_FDFBMAT3: calculate XIETA for ITEST 5')
      ENDIF
C
C----------------------------------------------
C     set up IXETRAN array:
C----------------------------------------------
C
C     IRELAX = 0
C     IF (LORX) THEN
C       IRELAX = IR1TAMP(LABEL,LORX,FREQ,ISYHOP)
C     END IF

      ! set the IXETRAN array for one (XI,ETA) pair
      IXETRAN(1,1) = IROPER(LABEL,ISYHOP)
      IXETRAN(2,1) = 0
      IXETRAN(3,1) = -1
      IXETRAN(4,1) = 0
      IXETRAN(5,1) = IRELAX
C
      LISTL  = 'L0 '
      FILXI  = 'CCFDFBO1'
      FILETA = 'CCFDFBX1'
      IOPTRES = 0
      NXETRAN = 1
C
C----------------------------------------------
C     initializations & work space allocations.
C----------------------------------------------
C
      CALL DZERO(RHO1,NT1AMX)
      CALL DZERO(RHO2,NT2AMX)
C
      IPRSAV = IPRINT
      IPRINT = 0
C
      MODEL = 'UNKNOWN'
      IF (CCS)  MODEL = 'CCS'
      IF (CC2)  MODEL = 'CC2'
      IF (CCSD) MODEL = 'CCSD'
C
      KT1AMPSAV  = 1
      KT2AMPSAV  = KT1AMPSAV + NT1AMX
      KT1AMPA    = KT2AMPSAV + NT2AMX
      KT2AMPA    = KT1AMPA   + NT1AMX
      KEND0      = KT2AMPA   + NT2AMX
      LWRK0      = LWORK     - KEND0
C
      KT0AMP0    = KEND0
      KT1AMP0    = KT0AMP0   + 2*NALLAI(1)
      KOMEGA1    = KT1AMP0   + NT1AMX
      KOMEGA2    = KOMEGA1   + NT1AMX
      KT2AMP0    = KOMEGA2   + MAX(NT2AMX,2*NT2ORT(1),NT2AO(1))
      KSCR2      = KT2AMP0   + NT2AMX
      KEND1A     = KSCR2     + NT2AMX + NT1AMX
      LWRK1A     = LWORK     - KEND1A
C
      KRHO1      = KEND0
      KRHO2      = KRHO1     + NT1AMX
      KEND1B     = KRHO2     + NT2AMX
      LWRK1B     = LWORK     - KEND1B
C
      IF (LWRK1A.LT.0 .OR. LWRK1B.LT.0) THEN
         WRITE(LUPRI,*) 'TOO LITTLE WORK SPACE IN CC_FDFBMAT3::'
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',MAX(KEND1A,KEND1B)
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDFBMAT2: ')
      ENDIF

C     -------------------------------------------
C     Read the CC reference amplitudes from disk:
C     -------------------------------------------
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AMPSAV),
     *              WORK(KT2AMPSAV))

C     --------------------------------------------
C     Read the T^A reference amplitudes from disk:
C     --------------------------------------------
      IOPT  = 3
      ISYMR = ILSTSYM(LISTR,IDLSTR)
      CALL CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,
     &              WORK(KT1AMPA),WORK(KT2AMPA))
      IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,1)

*---------------------------------------------------------------------*
*     Add delta x t^A to cluster amplitudes and recalculate the 
*     response intermediates and the Eta^B vector:
*               Eta^B{t^0 + delta x t^A} 
*---------------------------------------------------------------------*
C 
C     ---------------------------------
C     calculate the transformed vector:
C     ---------------------------------
      IORDER = 1
      CALL CC_XIETA(IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL,
     &              FILXI,  IDUM, RDUM,
     &              FILETA, IDUM, RDUM,
     &              .FALSE.,0, WORK(KEND0), LWRK0 )

C     IOPT   = 3
C     IVEC   = IXETRAN(4,1)
C     ISYETA = ISYHOP
C     CALL CC_RDRSP('X1 ',IVEC,ISYETA,IOPT,MODEL,
C    &              WORK(KRHO1),WORK(KRHO2))

      LEN       = NT1AMX + NT2AMX
      IADRF_ETA = IXETRAN(4,1)
      LUETA = -1
      CALL WOPEN2(LUETA,FILETA,64,0)
      CALL GETWA2(LUETA,FILETA,WORK(KRHO1),IADRF_ETA,LEN)
      CALL WCLOSE2(LUETA,FILETA,'KEEP')

      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      IF (.NOT.CCS) RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'Norm of RHO1(t^0 + delta x t^C): ',RHO1N
         WRITE (LUPRI,*) 'Norm of RHO2(t^0 + delta x t^C): ',RHO2N
      END IF
 
      ! divide by 2*delta and copy to result vector:
      CALL DSCAL(NT1AMX,HALF*DELINV,WORK(KRHO1),1)
      IF (.NOT.CCS) CALL DSCAL(NT2AMX,HALF*DELINV,WORK(KRHO2),1)

      CALL DCOPY(NT1AMX,WORK(KRHO1),1,RHODIF1,1)
      IF (.NOT.CCS) CALL DCOPY(NT2AMX,WORK(KRHO2),1,RHODIF2,1)
 
*---------------------------------------------------------------------*
*     Substract delta x t^C to cluster amplitudes and recalculate the 
*     response intermediates and the F matrix transformed T^B vector:
*               F{t^0 - delta x t^C} t^B 
*---------------------------------------------------------------------*

C     -------------------------------------------------------------
C     add finite displadement to t^0 and recalculate intermediates:
C     -------------------------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     &  CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
      CALL DAXPY(NT1AMX,-DELTA,WORK(KT1AMPA),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     &  CALL DAXPY(NT2AMX,-DELTA,WORK(KT2AMPA),1,WORK(KT2AMP0),1)
 
      IOPT = 3
      CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,WORK(KT0AMP0),
     &              WORK(KT1AMP0),WORK(KT2AMP0),WORK(KEND1A),LWRK1A)
 
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &     WORK(KT1AMP0),WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
 
C     ---------------------------------
C     calculate the transformed vector:
C     ---------------------------------
      IORDER = 1
      CALL CC_XIETA(IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL,
     &              FILXI,  IDUM, RDUM,
     &              FILETA, IDUM, RDUM,
     &              .FALSE.,0, WORK(KEND0), LWRK0 )

C     IOPT   = 3
C     IVEC   = IXETRAN(4,1)
C     ISYETA = ISYHOP
C     CALL CC_RDRSP('X1',IVEC,ISYETA,IOPT,MODEL,
C    &              WORK(KRHO1),WORK(KRHO2))

      LEN       = NT1AMX + NT2AMX
      IADRF_ETA = IXETRAN(4,1)
      CALL WOPEN2(LUETA,FILETA,64,0)
      CALL GETWA2(LUETA,FILETA,WORK(KRHO1),IADRF_ETA,LEN)
      CALL WCLOSE2(LUETA,FILETA,'DELETE')

      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      IF (.NOT.CCS) RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'Norm of RHO1(t^0 + delta x t^C): ',RHO1N
         WRITE (LUPRI,*) 'Norm of RHO2(t^0 + delta x t^C): ',RHO2N
      END IF
 
      ! divide by 2*delta and substract from final result:
      CALL DAXPY(NT1AMX,-HALF*DELINV,WORK(KRHO1),1,RHODIF1,1)
      IF (.NOT.CCS) 
     & CALL DAXPY(NT2AMX,-HALF*DELINV,WORK(KRHO2),1,RHODIF2,1)

*---------------------------------------------------------------------*
*     fix the scale factor of the diagonal, print some output and
*     restore t^0 amplitudes and response intermediates on file:
*---------------------------------------------------------------------*

C     -----------------------------------------------------------------
C     scale diagonal with 1/2: (only for right vectors --> A,B,C,D mat)
C     -----------------------------------------------------------------
C     CALL CCLR_DIASCL(RHODIF2,TWO,1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'RESULT VECTOR FROM CC_FDFBMAT2:'
         CALL CC_PRP(RHODIF1,RHODIF2,1,1,1)
      ENDIF
 
C     --------------------------------------------
C     Restore the CC reference amplitudes on disk:
C     --------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      IF (.NOT.CCS) 
     &  CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
 
      IOPT = 3
      CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,WORK(KT0AMP0),WORK(KT1AMP0),
     &              WORK(KT2AMP0),WORK(KEND1A),LWRK1A)
 
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),
     &     WORK(KT1AMP0),WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
 
      IF (IPRSAV .GT. 5) THEN
         CALL AROUND(' END OF CC_FDFBMAT2:')
      ENDIF
 
      IPRINT = IPRSAV
 
      RETURN
      END
*=====================================================================*
*                      END OF SUBROUTINE CC_FDFBMAT3:                 *
*=====================================================================*

