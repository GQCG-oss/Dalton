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
       SUBROUTINE CC_SETDORPS(LABELOP,SYM1ONLY,IPRINTCC)
*---------------------------------------------------------------------*
*
*    Purpose: set dorps.h common blocks for calculation of 
*             (derivative) two-electron integrals
*
*      LABELOP : operator label 
*         '1DHAMxxx'  -- geometric first derivatives
*         'dh/dBxxx'  -- magnetic first derivatives (for London orb.)
*
*      SYM1ONLY : if true, only totally symmetric integrals will
*                 be calculated (f.x. for gradient)
*
*      IPRINTCC : print level for integral program
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <nuclei.h>
#include <symmet.h>
#include <dorps.h>
#include <chrnos.h>
#include <exeinf.h>
#include <ccorb.h>
#include <cch2d.h>
#include <abainf.h>
#include <cbinuc.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_DERIV> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.) 

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      LOGICAL SYM1ONLY
      CHARACTER*8 LABELOP
      INTEGER IXYZ, IATOM, IPRINTCC, I

* analyze operator label:
      IF      ( LABELOP(1:5).EQ.'1DHAM' ) THEN

        MAXDIF = 1
        DO IATOM = 1, NUCIND
          DO IXYZ = 1, 3
            DOCOOR(IXYZ,IATOM) = .TRUE.
          END DO
        END DO

        DO I = 1, 3*MXCENT
          DOPERT(I,1) = .TRUE.
          DOPERT(I,2) = .FALSE.
        END DO

      ELSE IF ( LABELOP(1:5).EQ.'dh/dB' ) THEN

        MAXDIF = 1

      ELSE

        WRITE (LUPRI,*) 'Unknown operator label in CC_SETDORPS:',LABELOP
        CALL QUIT('Unknown operator label in CC_SETDORPS.')

      END IF

* initialize flags on common /CCH2D/:
      TKTIME = .TRUE.
      RETUR  = .FALSE.
      NODV   = .FALSE.
      NOPV   = .FALSE.
      NOCONT = .FALSE.
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0

* initialize print flag on common /CBINUC/:
      IPRINT = IPRINTCC

* initialize flags on common /ABAINF/:
      IF (SYM1ONLY) THEN
         DOSYM(1) = .TRUE.
         DO I = 2, 8
           DOSYM(I) = .FALSE.
         END DO
      ELSE
         DO I = 1, 8
           DOSYM(I) = .TRUE.
         END DO
      END IF

* initialize flags on common /DORPS/:
      IF (SYM1ONLY) THEN
         DOREPS(0) = .TRUE.
         DO I = 1, 7
           DOREPS(I) = .FALSE.
         END DO
      ELSE
         DO I = 0, 7
           DOREPS(I) = .TRUE.
         END DO
      END IF


      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_SETDORPS                          *
*=====================================================================*
