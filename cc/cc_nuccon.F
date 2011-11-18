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
      REAL*8  FUNCTION CC_NUCCON(LABEL,ISYM)
*---------------------------------------------------------------------*
*
*     Purpose: retrieve the nuclear contribution to the expectation
*              value of the operator LABEL with symmetry ISYM
*              
*     Christof Haettig, March 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <nuclei.h>
#include <dipole.h>
#include <symmet.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM
      CHARACTER LABEL*8

      REAL*8 ONE, ZERO

      PARAMETER(ONE=1.0D0, ZERO=0.0D0)

      INTEGER JSCOOR, JCOOR

*---------------------------------------------------------------------*Y
* if not total symmetric, return a zero:
*---------------------------------------------------------------------*Y
      IF (ISYM.NE.1) THEN

        CC_NUCCON = ZERO
         
        RETURN
         
      END IF

*---------------------------------------------------------------------*Y
* dipole gradient:
*---------------------------------------------------------------------*
      IF ( LABEL(4:6).EQ.'DPG') THEN

        READ(LABEL,'(I3)') JSCOOR
        IF (LABEL(8:8).EQ.'X') JCOOR = 1
        IF (LABEL(8:8).EQ.'Y') JCOOR = 2
        IF (LABEL(8:8).EQ.'Z') JCOOR = 3

        
        CC_NUCCON = DDIPN(IPTAX(JCOOR,1),JSCOOR)

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_NUCCON> LABEL, ISYM, JSCOOR, JCOOR:',
     &                         LABEL, ISYM, JSCOOR, JCOOR
           WRITE (LUPRI,*) 'CC_NUCCON> result:',CC_NUCCON
        END IF

*---------------------------------------------------------------------*Y
* gradient of the second moment of charge:
*---------------------------------------------------------------------*
      ELSE IF ( LABEL(3:5).EQ.'QDG') THEN

        CC_NUCCON = ZERO

*---------------------------------------------------------------------*Y
* gradient of the third moment of charge:
*---------------------------------------------------------------------*
      ELSE IF ( LABEL(3:5).EQ.'OCG') THEN

        CC_NUCCON = ZERO

*---------------------------------------------------------------------*Y
* shieldings:
*---------------------------------------------------------------------*
      ELSE IF ( LABEL(4:7).EQ.' NST') THEN

        CC_NUCCON = ZERO

*---------------------------------------------------------------------*Y
* gradient of the third moment of charge:
*---------------------------------------------------------------------*
      ELSE IF ( LABEL(2:6).EQ.'-CM1 ') THEN

        CC_NUCCON = ZERO

*---------------------------------------------------------------------*
* unknown operator: print error message and stop
*---------------------------------------------------------------------*
      ELSE
       WRITE (LUPRI,*) 'Error in CC_GET_NUCCON:'
       WRITE (LUPRI,*) 'No nuclear contrib. available for ',LABEL,
     &      'operator.'
       CALL QUIT('Unknown operator in CC_GET_NUCCON.')
      END IF

      RETURN
      END
*=====================================================================*
