C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
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
