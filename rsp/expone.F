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
#include <single.h>
C
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
      PARAMETER ( D1INF = 0.99999D0 , CKMXPR = 1.0D-12 )
      LOGICAL FNDLB2
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
      CHARACTER*8 RTNLBL(2)
C
      CALL QENTER('AVEPRP')
      IF (NAVTOT.LE.0) GO TO 9999
C
C ALLOCATE WORK SPACE
C
      KMO    = 1
      KWRK1  = KMO    + NBAST*NBAST
C     ... sufficient space for N2BASX if "SQUARE" below
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('AVEPRP',KWRK1-1,LWRK)
C
      CALL HEADER('Calculation of electronic expectation values',0)
      WRITE (LUPRI,'(A/A)')
     &  ' (note that to get e.g. a dipole moment you must multiply the',
     &  '  electronic number by -1 and add the nuclear contribution.)'
      DO 100 IOP    = 1,NAVTOT
C
C Read AO property integrals and check for antisymmetry and for
C spatial symmetry (as identified by abs largest element).
C
         IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                                  'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUPROP)
         IF ( FNDLB2(LBLAV(IOP),RTNLBL,LUPROP)) THEN
            IF (RTNLBL(2).EQ.'SQUARE  ') THEN
               CALL READT(LUPROP,N2BASX,WRK)
               IND = IDAMAX(N2BASX,WRK,1)
            ELSE IF (RTNLBL(2).EQ.'SYMMETRI') THEN
               CALL READT(LUPROP,NNBASX,WRK)
               IND = IDAMAX(NNBASX,WRK,1)
            ELSE IF (RTNLBL(2).EQ.'ANTISYMM') THEN
               WRITE(LUPRI,'(/3A)')
     *            ' Property ',LBLAV(IOP),' is antisymmetric'//
     *            ' and has therefore average value ZERO'
               GO TO 100
            ELSE
               CALL QUIT('AVEPRP : No antisymmetry label on LUPROP')
            END IF
            IF (ABS(WRK(IND)).GT.CKMXPR) THEN
               IF (RTNLBL(2).NE.'SQUARE  ') THEN
C                 triangular matrix
                  IIND  = INT( (SQRT(D1+D8*IND) - D1INF)/ D2 )
                  IIMAX = IIND*(IIND+1)/2
                  IF ( IIMAX .EQ. IND ) THEN
                     JIND = IIND
                  ELSE
                     IIND = IIND + 1
                     JIND = IND  - IIMAX
                  END IF
               ELSE
C                 square matrix
                  IIND = 1 + (IND-1)/NBAST
                  JIND = IND - (IIND-1)*NBAST
               END IF
               KSYMPT = MULD2H( ISAO(IIND) , ISAO(JIND) )
               IF ( KSYMPT .EQ. 1 ) THEN
C
C Read AO property integrals and transform to MO basis.
C Calculate average value of one electron operator.
C
                  KSYMOP = KSYMPT
                  CALL PRPGET (LBLAV(IOP),
     *                   CMO,WRK(KMO),ANTSYM,WRK(KWRK1),LWRK1)
C                 CALL PRPGET (WORD,CMO,PRPMO,ANTSYM,WRK,LWRK)
C
                  IPRONE = 0
                  CALL PRPONE(WRK(KMO),UDV,ONETOT,IPRONE,LBLAV(IOP))
               ELSE
                  WRITE(LUPRI,'(/3A,I3,A)')
     *            ' Property ',LBLAV(IOP),' has symmetry',KSYMPT,
     *            ' and has therefore average value ZERO.'
               END IF
            ELSE
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/A/3A,1P,D10.2)')
     *        ' ** WARNING ** Average value of property not calculated',
     *        ' - all elements on property file AOPROPER with label: ',
     *           LBLAV(IOP),' are abs smaller than :',CKMXPR
            END IF
         ELSE
            NWARN  = NWARN + 1
            WRITE(LUPRI,'(/A/3A)')
     *        ' ** WARNING ** Average value of property not calculated',
     *        ' - Label ',LBLAV(IOP),' was not found on AOPROPER.'
         END IF
 100  CONTINUE
 9999 IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('AVEPRP')
      RETURN
      END
