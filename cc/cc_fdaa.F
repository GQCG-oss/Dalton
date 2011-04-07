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
*---------------------------------------------------------------------*
c/* Deck CC_AATST */
*=====================================================================*
       SUBROUTINE CC_AATST(WORK,LWORK)
*---------------------------------------------------------------------*
*
* Purpose : provide some tests for the A{O} transformation
*           for more detailed information of the tests see below
*
*           noddy implementation for programmers use only...
*           ... to switch between the different test or for changing
*           the parameters the program must be recompiled...
*
* Christof Haettig, 1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "maxorb.h"
#include "ccroper.h"
#include "ccr1rsp.h"
#include "cco1rsp.h"
#include "ccx1rsp.h"
#include "ccfro.h"
#include "cclists.h"

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_AATST> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FREQ
      REAL ZERO, ONE, TWO, FOUR, FIVE, DUMMY
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION FREQ
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, FIVE, DUMMY
#endif 
      PARAMETER (ZERO = 0.0D0, ONE  = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0, FIVE = 5.0D0)

      CHARACTER*3 LISTR
      CHARACTER*8 LABEL
      LOGICAL LTWO, LORX
      INTEGER KFDIFF1, KFDIFF2, KNEW1, KNEW2, KOLD1, KOLD2, LWRK0
      INTEGER KEND1, LWRK1, IDLSTR, ISYHOP, IRELAX, N2VEC, ISYRES
      INTEGER ISYMR, ITEST, IORDER

* external functions:
      INTEGER IR1TAMP
      INTEGER IROPER
#if defined (SYS_CRAY) 
      REAL DSQRT, DDOT
#else 
      DOUBLE PRECISION DSQRT, DDOT
#endif 
    

*---------------------------------------------------------------------*
* set up information for rhs vectors for the different tests:
*---------------------------------------------------------------------*
*  number of simultaneous transformations:
c     NAATRAN = 1
*---------------------------------------------------------------------*
*  test 1: use zeroth-order Hamiltonian as perturbation.
*          (does not require that anything is done before CC_AATST)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST = 1
c     IORDER = 1
c     LABEL  = 'HAM0    '
c     LORX   = .TRUE.
c     LTWO   = .TRUE.
c     FREQ   = ZERO
c     ISYHOP = 1
c     LISTR  = 'R1 '
c     IDLSTR = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYMR)
c     ISYRES = MULD2H(ISYHOP,ISYMR)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 2: test the A{O} transformation for a non-relaxed one-elctron 
*          perturbation against the old implemenation.
*          (does only require that the integrals for the operator
*           are available on file)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ITEST = 2
      IORDER = 1
      LABEL  = 'ZDIPLEN '
      LORX   = .FALSE.
      LTWO   = .FALSE.
      FREQ   = ZERO
      ISYHOP = 1
      LISTR  = 'R1 '
      IDLSTR = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYMR)
      ISYRES = MULD2H(ISYHOP,ISYMR)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 3: set differentiated integrals to zero and test only the
*          orbital relaxation & reorthogonalization part of the 
*          A{O} matrix transformation
*          (requires that the CPHF equations for the `reference' 
*           operator, specified here, have been solved and the Kappa
*           vector is available on disc)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST = 3
c     IORDER = 1
c     LABEL  = 'ZDIPLEN '
c     LORX   = .TRUE.
c     LTWO   = .FALSE.
c     FREQ   = ZERO
c     ISYHOP = 1
c     LISTR  = 'R1 '
c     IDLSTR = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYMR)
c     ISYRES = MULD2H(ISYHOP,ISYMR)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 4: test the A{O} for a `orbital-relaxed'
*          one-electron perturbation.
*          (requires that the CPHF equations for the operator specified 
*           have been solved and the Kappa vector is available on disc)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST = 4
c     IORDER = 1
c     LABEL  = 'ZDIPLEN '
c     LORX   = .TRUE.
c     LTWO   = .FALSE.
c     FREQ   = ZERO
c     ISYHOP = 1
c     LISTR  = 'R1 '
c     IDLSTR = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYMR)
c     ISYRES = MULD2H(ISYHOP,ISYMR)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
*
*  test 5: test A{O} for a second-order perturbation
*          operator made up from a relaxed field with pert.-dep. basis
*          and an unrelaxed one-electron perturbation
*          (requires that the CPHF equations for the operator specified 
*           have been solved and the Kappa vector is available on disc)
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c     ITEST  = 5
c     IORDER = 2
c     LABEL1 = '1DHAM003'
c     LORX1  = .TRUE.
c     LTWO1  = .TRUE.
c     FREQ1  = ZERO
c     LABEL2 = 'ZDIPLEN '
c     LORX2  = .FALSE.
c     LTWO2  = .FALSE.
c     FREQ2  = ZERO
c     CALL CC_FIND_SO_OP(LABEL1,LABEL2,LABEL,ISYHOP,ISIGN,
c    &                   INUM,WORK,LWORK)
c     LORX   = ( LORX1 .AND. LORX2 )
c     LTWO   = ( LTWO1 .AND. LTWO2 )
c     FREQ   = FREQ1 + FREQ2
c     LISTR  = 'R1 '
c     IDLSTR = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYMR)
c     ISYRES = MULD2H(ISYHOP,ISYMR)
*---------------------------------------------------------------------*

      ! allow extensions in the vector lists for the next few lines
      LOPROPN = .TRUE.
      LO1OPN  = .TRUE.
      LX1OPN  = .TRUE.

      IRELAX = 0
      IF (LORX) THEN
        IRELAX = IR1TAMP(LABEL,LORX,FREQ,ISYHOP)
      END IF

      LPDBSOP(IROPER(LABEL,ISYHOP)) = LTWO

      ! disallow again extension in the vector lists
      LOPROPN = .FALSE.
      LO1OPN  = .FALSE.
      LX1OPN  = .FALSE.
C
      CALL AROUND('CC_AATST: test of A{O} transformation')
      WRITE (LUPRI,*) 'ITEST  =',ITEST
      WRITE (LUPRI,*) 'LABEL  =',LABEL
      WRITE (LUPRI,*) 'ISYHOP =',ISYHOP
      WRITE (LUPRI,*) 'LTWO   =',LTWO
      WRITE (LUPRI,*) 'LORX   =',LORX
      WRITE (LUPRI,*) 'FREQ   =',FREQ
      WRITE (LUPRI,*) 'IRELAX =',IRELAX
      WRITE (LUPRI,*) 'LISTR  =',LISTR 
      WRITE (LUPRI,*) 'IDLSTR =',IDLSTR
C 
      N2VEC = 1
      IF (CCS) N2VEC = 0

*---------------------------------------------------------------------*
*     allocate work space to keep the different result vectors:
*---------------------------------------------------------------------*
      KFDIFF1 = 1
      KFDIFF2 = KFDIFF1 + NT1AM(ISYRES)
      KNEW1   = KFDIFF2 + NT2AM(ISYRES)
      KNEW2   = KNEW1   + NT1AM(ISYRES)
      KEND1   = KNEW2   + NT2AM(ISYRES)
      LWRK1   = LWORK   - KEND1

*---------------------------------------------------------------------*
*     call new A{O} x T^B transformtion routine... (not yet available)
*---------------------------------------------------------------------*
c     CALL ....

*---------------------------------------------------------------------*
*     call the finite difference impementation based on CC_XIETA:
*---------------------------------------------------------------------*
      CALL CC_FDAAMAT(LISTR,IDLSTR,WORK(KFDIFF1),WORK(KFDIFF2),
     &                LABEL,IRELAX,WORK(KEND1),LWRK1)

      CALL AROUND('fin. diff. A{O} x T^B transformation result:')
      CALL CC_PRP(WORK(KFDIFF1),WORK(KFDIFF2),ISYRES,1,N2VEC)

*---------------------------------------------------------------------*
*     call old implementation for
*---------------------------------------------------------------------*
      IF (ITEST.EQ.2) THEN

        LWRK0 = LWRK1

        KOLD1 = KEND1
        KOLD2 = KOLD1 + NT1AM(ISYRES)
        KEND1 = KOLD2 + NT2AM(ISYRES)
        LWRK1 = LWORK - KEND1

        CALL CCCR_AA(LABEL,ISYHOP,LISTR,IDLSTR,DUMMY,WORK(KOLD1),LWRK0)

        CALL CCLR_DIASCL(WORK(KOLD2),TWO,ISYRES)

        CALL AROUND('old CCCR_AA one-electron A{O} implementation:')
        CALL CC_PRP(WORK(KOLD1),WORK(KOLD2),ISYRES,1,N2VEC)

        CALL DAXPY(NT1AM(ISYRES),-ONE,WORK(KFDIFF1),1,WORK(KOLD1),1)
        CALL DAXPY(NT2AM(ISYRES),-ONE,WORK(KFDIFF2),1,WORK(KOLD2),1)

        WRITE (LUPRI,*) 'Norm of difference between fin. diff. and old:'
        WRITE (LUPRI,*) ' single excitation part:',
     &    DSQRT(DDOT(NT1AM(ISYRES),WORK(KOLD1),1,WORK(KOLD1),1))
        IF (.NOT.CCS) WRITE (LUPRI,*) ' double excitation part:',
     &    DSQRT(DDOT(NT2AM(ISYRES),WORK(KOLD2),1,WORK(KOLD2),1))
        WRITE (LUPRI,*) ' difference vector:'
        CALL CC_PRP(WORK(KOLD1),WORK(KOLD2),ISYRES,1,N2VEC)

      END IF

      RETURN
      END 
*=====================================================================*

*=====================================================================*
      SUBROUTINE CC_FDAAMAT(LISTR,IDLSTR,RHODIF1,RHODIF2,
     &                      LABEL,IRELAX,WORK,LWORK)
C
C---------------------------------------------------------------------
C Test routine for calculating a generalized CC A{O} matrix transformed 
C vector by finite difference on the Xi{O} vector
C
C Christof Haettig, july 1999
C---------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
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
         CALL AROUND('in CC_FDAAMAT: making fini. diff. A{O} Matrix')
      ENDIF
C
C----------------------------------------------
C     set up IXETRAN array:
C----------------------------------------------
C
      ! set the IXETRAN array for one (XI,ETA) pair
      IXETRAN(1,1) = IROPER(LABEL,ISYHOP)
      IXETRAN(2,1) = 0
      IXETRAN(3,1) = 0
      IXETRAN(4,1) = -1
      IXETRAN(5,1) = IRELAX
C
      LISTL  = 'L0 '
      FILXI  = 'CCFDAAO1'
      FILETA = 'CCFDAAX1'
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
      IF (CCR12) CALL QUIT('Finite-difference A{O}-matrix for CCR12 '//
     &                     'not adapted')
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
         WRITE(LUPRI,*) 'TOO LITTLE WORK SPACE IN CC_FDAAMAT:'
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',MAX(KEND1A,KEND1B)
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDAAMAT: ')
      ENDIF

C     -------------------------------------------
C     Read the CC reference amplitudes from disk:
C     -------------------------------------------
      IOPT = 3
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMPSAV),
     &              WORK(KT2AMPSAV))
 
C     --------------------------------------------
C     Read the T^A reference amplitudes from disk:
C     --------------------------------------------
      IOPT  = 3
      ISYMR = ILSTSYM(LISTR,IDLSTR)
      CALL CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,
     &              WORK(KT1AMPA),WORK(KT2AMPA))
      CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,1)

*---------------------------------------------------------------------*
*     Add delta x t^A to cluster amplitudes and recalculate the 
*     response intermediates and the Xi^O vector:
*               Xi^O{t^0 + delta x t^A} 
*---------------------------------------------------------------------*

C     -------------------------------------------------------------
C     add finite displadement to t^0 and recalculate intermediates:
C     -------------------------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
      CALL DAXPY(NT1AMX,DELTA,WORK(KT1AMPA),1,WORK(KT1AMP0),1)
      CALL DAXPY(NT2AMX,DELTA,WORK(KT2AMPA),1,WORK(KT2AMP0),1)
C
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

      LEN       = NT1AMX + NT2AMX
      IADRF_XI  = IXETRAN(3,1)
      LUETA = -1
      CALL WOPEN2(LUETA,FILXI,64,0)
      CALL GETWA2(LUETA,FILXI,WORK(KRHO1),IADRF_XI,LEN)
      CALL WCLOSE2(LUETA,FILXI,'KEEP')

      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'Norm of RHO1(t^0 + delta x t^C): ',RHO1N
         WRITE (LUPRI,*) 'Norm of RHO2(t^0 + delta x t^C): ',RHO2N
      END IF
 
      ! divide by 2*delta and copy to result vector:
      CALL DSCAL(NT1AMX,HALF*DELINV,WORK(KRHO1),1)
      CALL DSCAL(NT2AMX,HALF*DELINV,WORK(KRHO2),1)

      CALL DCOPY(NT1AMX,WORK(KRHO1),1,RHODIF1,1)
      CALL DCOPY(NT2AMX,WORK(KRHO2),1,RHODIF2,1)
 
*---------------------------------------------------------------------*
*     Substract delta x t^C to cluster amplitudes and recalculate the 
*     response intermediates and Xi^O vector:
*               Xi^O{t^0 + delta x t^A} 
*---------------------------------------------------------------------*

C     -------------------------------------------------------------
C     add finite displadement to t^0 and recalculate intermediates:
C     -------------------------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
      CALL DAXPY(NT1AMX,-DELTA,WORK(KT1AMPA),1,WORK(KT1AMP0),1)
      CALL DAXPY(NT2AMX,-DELTA,WORK(KT2AMPA),1,WORK(KT2AMP0),1)
 
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

      LEN       = NT1AMX + NT2AMX
      IADRF_XI  = IXETRAN(4,1)
      CALL WOPEN2(LUETA,FILXI,64,0)
      CALL GETWA2(LUETA,FILXI,WORK(KRHO1),IADRF_XI,LEN)
      CALL WCLOSE2(LUETA,FILXI,'DELETE')

      RHO1N = DDOT(NT1AMX,WORK(KRHO1),1,WORK(KRHO1),1)
      RHO2N = DDOT(NT2AMX,WORK(KRHO2),1,WORK(KRHO2),1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'Norm of RHO1(t^0 + delta x t^C): ',RHO1N
         WRITE (LUPRI,*) 'Norm of RHO2(t^0 + delta x t^C): ',RHO2N
      END IF
 
      ! divide by 2*delta and substract from final result:
      CALL DAXPY(NT1AMX,-HALF*DELINV,WORK(KRHO1),1,RHODIF1,1)
      CALL DAXPY(NT2AMX,-HALF*DELINV,WORK(KRHO2),1,RHODIF2,1)

*---------------------------------------------------------------------*
*     fix the scale factor of the diagonal, print some output and
*     restore t^0 amplitudes and response intermediates on file:
*---------------------------------------------------------------------*

C     -----------------------------------------------------------------
C     scale diagonal with 1/2: (only for right vectors --> A,B,C,D mat)
C     -----------------------------------------------------------------
      CALL CCLR_DIASCL(RHODIF2,TWO,1)

      IF (IPRSAV.GT.10) THEN
         WRITE (LUPRI,*) 'RESULT VECTOR FROM CC_FDAAMAT:'
         CALL CC_PRP(RHODIF1,RHODIF2,1,1,1)
      ENDIF
 
C     --------------------------------------------
C     Restore the CC reference amplitudes on disk:
C     --------------------------------------------
      CALL DCOPY(NT1AMX,WORK(KT1AMPSAV),1,WORK(KT1AMP0),1)
      CALL DCOPY(NT2AMX,WORK(KT2AMPSAV),1,WORK(KT2AMP0),1)
 
      IOPT = 3
      CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,WORK(KT0AMP0),WORK(KT1AMP0),
     &              WORK(KT2AMP0),WORK(KEND1A),LWRK1A)
 
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KOMEGA1),WORK(KOMEGA2),WORK(KT1AMP0),
     *            WORK(KT2AMP0),WORK(KEND1A),LWRK1A,'XXX')
 
      IF (IPRSAV .GT. 5) THEN
         CALL AROUND(' END OF CC_FDAAMAT:')
      ENDIF
 
      IPRINT = IPRSAV
 
      RETURN
      END
*=====================================================================*
*                      END OF SUBROUTINE CC_FDAAMAT                   *
*=====================================================================*
