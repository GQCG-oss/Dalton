*=====================================================================*
#if defined (SYS_CRAY)
      REAL FUNCTION CC_NUCCON(LABEL,ISYM)
#else
      DOUBLE PRECISION FUNCTION CC_NUCCON(LABEL,ISYM)
#endif
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
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <maxaqn_dec.h>
#  include <mxcent_dec.h>
#  include <mxorb_dec.h>
#  include <nuclei_dec.h>
#  include <symmet_dec.h>
#  include <dipole_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <nuclei.h>
#include <dipole.h>
#include <symmet.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM
      CHARACTER LABEL*8

#if defined (SYS_CRAY)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION ONE, ZERO
#endif
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
      ELSE IF ( LABEL(3:7).EQ.' NST ') THEN

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
       CALL QUIT('Unkown operator in CC_GET_NUCCON.')
      END IF

      RETURN
      END
*=====================================================================*
