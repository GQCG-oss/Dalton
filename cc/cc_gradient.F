c /* deck cc_gradient */
*=====================================================================*
       SUBROUTINE CC_GRADIENT(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: calculation of the gradient
*
*             implemented wavefunctions models:
*                        CCSD
*
*     Written by Asger Halkier and Christof Haettig summer 1996.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <maxaqn_dec.h>
#  include <mxcent_dec.h>
#  include <nuclei_dec.h>
#  include <energy_dec.h>
#  include <taymol_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <nuclei.h>
#include <energy.h>
#include <taymol.h>
#include <cch2d.h>

* local parameters:

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0, GRDNRM
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0, GRDNRM
#endif

      LOGICAL SYM1ONLY
      INTEGER MAXDIF, IOPT

* external functions
#if defined (SYS_CRAY)
      REAL SECOND
      REAL DDOT
#else
      DOUBLE PRECISION SECOND
      DOUBLE PRECISION DDOT
#endif

*---------------------------------------------------------------------*
* print header for gradient section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<  OUTPUT FROM COUPLED CLUSTER ',
     &                                  'GRADIENT SECTION  >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*

* test for wavefunction model: 
      IF (.NOT. (CCS. OR. MP2. OR. cc2. OR. CCD .OR. CCSD) ) THEN
         CALL QUIT('CC_GRADIENT called for a model '//
     &             'not (yet) implemented.')
      END IF

* print available work space:
      IF (IPRINT .GT. 10) THEN
        WRITE(LUPRI,*) 'Available work space in CC_GRADIENT section:',
     &        LWORK
      END IF
  
* initialize timing:
      TIM0  = SECOND()

* initialize ABACUS common blocks DORPS and CBINUC: 
* set differentiation level to 1 and print level to IPRINT
* and variables on DORPS common block as needed gradient
C     MAXDIF = 1
C     CALL CCDORPSINI(MAXDIF,IPRINT,'GRADIENT')

      SYM1ONLY = .TRUE.
      CALL CC_SETDORPS('1DHAM   ',SYM1ONLY,IPRINT)

* initialize /CCH2D/ common block, controlling the calculation
* of derivative 2-el integrals for the CC program:
      TKTIME = .TRUE.
      RETUR  = .FALSE.
      NODV   = .FALSE.
      NOPV   = .FALSE.
      NOCONT = .FALSE.
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0

*---------------------------------------------------------------------*
* call ABACUS to calculate nuclear repulsion contribution to gradient:
* result is returned on GRADNN on ABACUS /ENERGY/ common block...
*---------------------------------------------------------------------*
      CALL NUCREP(WORK,WORK(MXCOOR*MXCOOR+1),WORK(2*MXCOOR*MXCOOR+1))

*---------------------------------------------------------------------*
* call CC_GRAD0 to calculate:
*        GRADKE -- contrib. from derivative of one-electron hamilt.
*        GRADEE -- contrib. from derivative of two-electron integr.
*        GRADFS -- contrib. from overlap derivative integr.
* arrays are allocated on ABACUS /ENERGY/ common block, GRADNA not used
*---------------------------------------------------------------------*
C
C==================================================================
C     Asger Halkier January 99. Change to two-step algorithm. First
C     call calculates everything else than GRADEE, which is instead
C     calculated in the second call to the new routine CC_GRAD2E.
C==================================================================
C     IOPT = 1
C     CALL CC_GRAD(GRADKE,GRADFS,WORK,LWORK,IOPT)
C     IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'Home from CC_GRAD.'
C     CALL CC_GRAD2E(GRADEE,WORK,LWORK)
C     IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'Home from CC_GRAD2E.'
C     CALL FLSHFO(LUPRI)
C==================================================================
C     SC & CH: Change back to one-step algorithm. Keep CC_GRAD in
C              two steps as used for DPT
C==================================================================
C
      CALL CC_GRAD0(GRADKE,GRADFS,GRADEE,WORK,LWORK)
      IF (IPRINT .GT. 10) 
     &    WRITE(LUPRI,*) ' Home from CC_GRAD0.'
      CALL FLSHFO(LUPRI)
C
*---------------------------------------------------------------------*
* zero nuclear attraction contribution since it is included in GRADKE
*---------------------------------------------------------------------*
      CALL DZERO(GRADNA,MXCOOR)

*---------------------------------------------------------------------*
* print gradient contributions (GRADNN is already printed by NUCREP)
*---------------------------------------------------------------------*
      IF (IPRINT .GT. 1) THEN
         CALL HEADER('One-electron integral gradient',-1)
         CALL PRIGRD(GRADKE,WORK,WORK(MXCOOR*MXCOOR+1))
         CALL HEADER('Two-electron integral gradient',-1)
         CALL PRIGRD(GRADEE,WORK,WORK(MXCOOR*MXCOOR+1))
         CALL HEADER('Reorthonormalization gradient',-1)
         CALL PRIGRD(GRADFS,WORK,WORK(MXCOOR*MXCOOR+1))
      END IF

*---------------------------------------------------------------------*
* calculated total gradient using ABACUS routines:
* (uses ABACUS common block /TAYMOL/)
*---------------------------------------------------------------------*
      CALL ZERGRD
      CALL ADDGRD(GRADNN)
      CALL ADDGRD(GRADEE)
      CALL ADDGRD(GRADFS)
      CALL ADDGRD(GRADKE)

      CALL HEADER('Molecular gradient',-1)
      CALL PRIGRD(GRDMOL,WORK,WORK(MXCOOR*MXCOOR+1))

      GRDNRM = DDOT(3*NUCDEP,GRDMOL,1,GRDMOL,1)

      WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *   'Molecular gradient norm:', GRDNRM

*---------------------------------------------------------------------*
* print timing & return:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,A,F12.2," seconds.")') '>>> Total time used ',
     &  'in gradient section:', SECOND() - TIM0

      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_GRADIENT                          *
*=====================================================================*
c /* deck ccdorpsini */
*=====================================================================*
       SUBROUTINE CCDORPSINI(MAXDIFCC,IPRINTCC,KEY)
*---------------------------------------------------------------------*
*
*  Purpose: initialize ABACUS common blocks /DOPRS/, /CBINUC/:
*
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <symmet.h>
#include <cbinuc.h>
#include <dorps.h>
#include <abainf.h>

       CHARACTER KEY*(*)
       INTEGER MAXDIFCC, IPRINTCC

*---------------------------------------------------------------------*
* copy input variables to /CBINUC/ common block:
*---------------------------------------------------------------------*
       MAXDIF = MAXDIFCC
       IPRINT = IPRINTCC

*---------------------------------------------------------------------*
* set variables on /DORPS/ common block:
*---------------------------------------------------------------------*
       IF ( KEY(1:8) .EQ. 'GRADIENT' ) THEN

C        ------------------------------------------------------------
C        for gradient we need only derivatives for totally symmetric
C        symmetry coordinates. note that in eri irepps go from 0 to 7
C        ------------------------------------------------------------
         MAXREP   = 7

         DOSYM(1) = .TRUE.
         DO IREP = 2, 8
           DOSYM(IREP) = .TRUE.
         END DO

         DOREPS(0) = .TRUE.
         DO IREP = 1, 7
           DOREPS(IREP) = .TRUE.
         END DO

         WRITE (LUPRI,*) 'DOREPS set.'
         CALL FLSHFO(LUPRI)

         DO I = 1, 3*MXCENT
           DOPERT(I,1) = .TRUE.
           DOPERT(I,2) = .FALSE.
         END DO

         WRITE (LUPRI,*) 'DOPERT set.'
         CALL FLSHFO(LUPRI)

         DO IATOM = 1, NUCIND
            DO IXYZ = 1, 3
               ICOOR = ICOOR + 1
               DOCOOR(IXYZ,IATOM) = .TRUE.
C              DO IREP = 0, 7
C                IF (DOREPS(IREP).AND.(IPTCNT(ICOOR,IREP,1).GT.0)) THEN
C                  DOCOOR(IXYZ,IATOM) = .TRUE.
C                END IF
C              END DO
            END DO
         END DO

         WRITE (LUPRI,*) 'DOCOOR set.'
         CALL FLSHFO(LUPRI)

       ELSE
         CALL QUIT('Unkown KEY in CCNUCINI.')
       END IF

       RETURN
       END
*---------------------------------------------------------------------*
C  /* Deck cc_geoflag */
      SUBROUTINE CC_GEOFLAG
C
C     Written by Asger Halkier & Christof Haettig 19/8 - 1998
C
C     Purpose: To set logical flag for integral file in case of
C              geometry optimization.
C              set RSPIM to enforce calculation of response 
C              intermediates together with the energy
C
#include <implicit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      KEEPAOTWO = MAX(KEEPAOTWO,1)
      RSPIM     = .TRUE.
      KEEPAOIN  = .TRUE.
C
      RETURN
      END
*---------------------------------------------------------------------*
