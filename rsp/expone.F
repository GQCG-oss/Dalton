C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
C
C$Id: expone.F,v 1.1.1.1 2001-02-08 13:33:14 hjj Exp $
C
#include <single.h>
C  /* Deck aveprp */
      SUBROUTINE AVEPRP(CMO,UDV,WRK,LWRK)
C
C CALCULATE AVERAGE VALUE OF PROPERTIES
C
C
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 , D8 = 8.0D0, D0 = 0.0D0 )
      PARAMETER ( D1INF = 0.99999D0 , CKMXPR = 1.0D-6 )
      LOGICAL FNDLAB
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infave.h>
#include <inftap.h>
C
      DIMENSION CMO(*),UDV(*),WRK(*)
C
      IF (NAVTOT.LE.0) RETURN
C
C ALLOCATE WORK SPACE
C
      KMO    = 1
      KWRK1  = KMO    + NORBT*NORBT
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('AVEPRP',KWRK1-1,LWRK)
      IF (NAVTOT .GT. 0) THEN
         WRITE (LUPRI,'(/A//A/A)')
     &      ' *** Calculation of electronic expectation values ***',
     &      ' (note that to get e.g. a dipole moment you must multiply',
     &      '  this number by -1 and add the nuclear contribution.)'
      END IF
      CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      DO 100 IOP    = 1,NAVTOT
C
C Read AO property integrals
C
         IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                                  'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUPROP)
         IF ( FNDLAB(LBLAV(IOP),LUPROP) )  THEN
            CALL READT(LUPROP,NNBASX,WRK)
            IND   = IDAMAX(NNBASX,WRK,1)
            IF (ABS(WRK(IND)).GT.CKMXPR) THEN
               IIND  = INT( (SQRT(D1+D8*IND) - D1INF)/ D2 )
               IIMAX = IIND*(IIND+1)/2
               IF ( IIMAX .EQ. IND ) THEN
                  JIND = IIND
               ELSE
                  IIND = IIND + 1
                  JIND = IND  - IIMAX
               END IF
               KSYMPT = MULD2H( ISAO(IIND) , ISAO(JIND) )
               IF ( KSYMPT .EQ. 1 ) THEN
C
C Read AO property integrals and transform to MO basis.
C
                  KSYMOP = KSYMPT
                  CALL PRPGET (LBLAV(IOP),
     *                   CMO,WRK(KMO),ANTSYM,WRK(KWRK1),LWRK1)
C                 CALL PRPGET (WORD,CMO,PRPMO,WRK,LWRK)
C
C
C CALCULATE AVERAGE VALUE OF ONE ELECTRON OPERATOR
C
                  IPRONE = 0
                  IF (ANTSYM .GT. D0) THEN
                     WRITE(LUPRI,'(/3A)')
     *                    ' Property ',LBLAV(IOP),' is antisymmetric'//
     *                    ' and has therefore average value ZERO'
                  ELSE
                     CALL PRPONE(WRK(KMO),UDV,ONETOT,IPRONE,LBLAV(IOP))
                  END IF
               ELSE
                  WRITE(LUPRI,'(/3A,I3,A)')
     *            ' Property ',LBLAV(IOP),' has symmetry',KSYMPT,
     *            ' and has therefore average value ZERO'
               END IF
            ELSE
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/A/3A,1P,D15.7)')
     *        ' ** WARNING ** AVERAGE VALUE OF PROPERTY NOT CALCULATED',
     *        ' - ELEMENTS ON PROPERTY FILE  WITH LABEL: ',
     *           LBLAV(IOP),' ARE SMALLER THAN :',CKMXPR
            END IF
         ELSE
            NWARN  = NWARN + 1
            WRITE(LUPRI,'(/A/3A)')
     *        ' ** WARNING ** AVERAGE VALUE OF PROPERTY NOT CALCULATED',
     *        ' - LABEL ',LBLAV(IOP),' WAS NOT FOUND ON LUPROP.'
         END IF
 100  CONTINUE
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      RETURN
      END
