C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
#include <single.h>
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
