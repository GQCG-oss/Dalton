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
#include "single.h"
C
C  /* Deck aveprp */
      SUBROUTINE AVEPRP(CMO,UDV,WRK,LWRK)
C
C Calculate average value (expectation value) of properties
C specified with .PROPAV under **RESPONS
C
C Revised Mar 2003 Hans Joergen Aa. Jensen
C
#include "implicit.h"
#include "priunit.h"
C
#include "infave.h"
C
      DIMENSION CMO(*),UDV(*),WRK(*)
C
      IF (NAVTOT.LE.0) RETURN
      CALL QENTER('AVEPRP')
C
      CALL HEADER('Calculation of electronic expectation values',0)
      WRITE (LUPRI,'(A/A)')
     &  ' (note that to get e.g. a dipole moment you must multiply the',
     &  '  electronic number by -1 and add the nuclear contribution.)'
C
      DO 100 IOP    = 1,NAVTOT
         CALL PRPAVE(LBLAV(IOP),ONETOT,CMO,UDV,WRK,LWRK)
 100  CONTINUE
C
      CALL QEXIT('AVEPRP')
      RETURN
      END
C  /* Deck prpave */
      SUBROUTINE PRPAVE(LBLAVE,AVEVAL,CMO,UDV,WRK,LWRK)
C
C Calculate average value (expectation value) of property
C with label 'LBLAVE' on AOPROPER
C
C Revised March 2003, Hans Joergen Aa. Jensen.
C
#include "implicit.h"
#include "dummy.h"
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 , D8 = 8.0D0, D0 = 0.0D0 )
      PARAMETER ( D1INF = 0.999999D0 , CKMXPR = 1.0D-12 )
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "inforb.h"
#include "infrsp.h"
#include "infind.h"
#include "wrkrsp.h"
#include "inftap.h"
C
      DIMENSION CMO(*),UDV(*),WRK(*)
      CHARACTER*8 LBLAVE, RTNLBL(2)
      LOGICAL FNDLB2
C
      CALL QENTER('PRPAVE')
C
C Allocate work spacE
C
      KMO    = 1
      KWRK1  = KMO    + NBAST*NBAST
C     ... sufficient space for N2BASX if "SQUARE" below
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('PRPAVE',KWRK1-1,LWRK)
C
C Read AO property integrals for LBLAVE and check for antisymmetry and
C for spatial symmetry (as identified by abs largest element).
C
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND (LUPROP)
      IF ( FNDLB2(LBLAVE,RTNLBL,LUPROP)) THEN
         IF (RTNLBL(2).EQ.'SQUARE  ') THEN
            CALL READT(LUPROP,N2BASX,WRK)
            IND = IDAMAX(N2BASX,WRK,1)
         ELSE IF (RTNLBL(2).EQ.'SYMMETRI') THEN
            CALL READT(LUPROP,NNBASX,WRK)
            IND = IDAMAX(NNBASX,WRK,1)
         ELSE IF (RTNLBL(2).EQ.'ANTISYMM') THEN
            WRITE(LUPRI,'(/3A)')
     *         ' Property ',LBLAVE,' is antisymmetric'//
     *         ' and has therefore average value ZERO'
            GO TO 100
         ELSE
            CALL QUIT('PRPAVE : No antisymmetry label on LUPROP')
         END IF
         IF (ABS(WRK(IND)).GT.CKMXPR) THEN
            IF (RTNLBL(2).NE.'SQUARE  ') THEN
C              triangular matrix
C              find solution to iind*(iind+1)/2 = ind :
               IIND  = INT( (SQRT(D1+D8*IND) - D1INF)/ D2 )
C              ... use D1INF to avoid round-off errors for diagonal elements
               IIMAX = IIND*(IIND+1)/2
               IF ( IIMAX .EQ. IND ) THEN
                  JIND = IIND
               ELSE
                  IIND = IIND + 1
                  JIND = IND  - IIMAX
               END IF
            ELSE
C              square matrix
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
               CALL PRPGET (LBLAVE,CMO,WRK(KMO),ANTSYM,WRK(KWRK1),LWRK1)
C              CALL PRPGET (WORD,CMO,PRPMO,ANTSYM,WRK,LWRK)
C
               IPRONE = 0
               CALL PRPONE(WRK(KMO),UDV,AVEVAL,IPRONE,LBLAVE)
            ELSE
               WRITE(LUPRI,'(/3A,I3,A)')
     *         ' Property ',LBLAVE,' has symmetry',KSYMPT,
     *         ' and has therefore average value ZERO.'
            END IF
         ELSE
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/A/3A,1P,D10.2)')
     *        ' ** WARNING ** Average value of property not calculated',
     *        ' - all elements on property file AOPROPER with label: ',
     *           LBLAVE,' are abs smaller than :',CKMXPR
            END IF
         ELSE
            NWARN  = NWARN + 1
            WRITE(LUPRI,'(/A/3A)')
     *        ' ** WARNING ** Average value of property not calculated',
     *        ' - Label ',LBLAVE,' was not found on AOPROPER.'
         END IF
 100  CONTINUE
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('PRPAVE')
      RETURN
      END
C  /* Deck prpone */
      SUBROUTINE PRPONE(PRPMO,UDV,ONETOT,IPRONE,PRPLBL)
C
C CALCULATE AVERAGE VALUE OF ONE ELECTRON OPERATOR in PRPMO
C
#include "implicit.h"
C
      CHARACTER*(*)PRPLBL
C
      DIMENSION PRPMO(NORBT,NORBT),UDV(NASHDI,NASHDI)
C
#include "priunit.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infdim.h"
#include "infpri.h"
C
      PARAMETER ( D2 = 2.0D0 , D0 = 0.0D0 )
      PARAMETER ( BIGLIM = 100000.0D0, SMLLIM = 0.01D0 )
C
      IPRORB = IPRONE
      IF (NOCCT .GT. 12) IPRORB = MAX(IPRORB,2)
C     hjaaj Mar 2003: default IPRRSP is 2, do not print individual
C     contributions by default if more than 12 occupied orbitals.
      IF (IPRRSP.GT.IPRORB) THEN
         WRITE(LUPRI,'(/A/3A/)')
     &      ' *** Individual non-zero orbital contributions',
     &      ' *** to the expectation value for property ',PRPLBL,' :'
      END IF
C
      ONEACT = D0
      ONEINA = D0
      ONESEC = D0
      DO 50 ISYM = 1,NSYM
         JSYM = MULD2H(KSYMOP,ISYM)
         NASHI = NASH(ISYM)
         NISHI = NISH(ISYM)
         IORBI = IORB(ISYM)
         IASHI = IASH(ISYM)
         NASHJ = NASH(JSYM)
         NISHJ = NISH(JSYM)
         IORBJ = IORB(JSYM)
         IASHJ = IASH(JSYM)
         IF ( .NOT. TRPLET ) THEN
            DO IINAC = 1,NISHI
               TMP = PRPMO(IORBI+IINAC,IORBI+IINAC) * D2
               IF (IPRRSP.GT.IPRORB .AND. TMP .NE. D0)
     &            WRITE(LUPRI,'(5X,A,2I3,A,I2,A,F15.8)')
     &            'Inactive     ',IINAC,IINAC,' in sym',ISYM,' :',TMP
               ONEINA = ONEINA + TMP
            END DO
         END IF
         DO 60 IA = 1,NASHI
            DO 70 JA = 1,NASHJ
               TMP = UDV(IASHI+IA,IASHJ+JA) *
     *               PRPMO(IORBI+NISHI+IA,IORBJ+NISHJ+JA)
               IF (IPRRSP.GT.IPRORB .AND. TMP .NE. D0)
     &            WRITE(LUPRI,'(5X,A,2I3,A,I2,A,F15.8)')
     &            'Active-active',IA,JA,' in sym',ISYM,' :',TMP
               ONEACT = ONEACT + TMP
 70         CONTINUE
 60      CONTINUE
 50   CONTINUE
C
      IF (SOPPA) CALL ONEMP2(PRPMO,UDV,ONESEC)
C
      ONETOT = ONEINA + ONEACT + ONESEC
      IF (IPRRSP.GE.IPRONE) THEN
         IF (ABS(ONETOT) .GT. SMLLIM .AND. ABS(ONETOT) .LT. BIGLIM)
     *                                                       THEN
            IF (.NOT. SOPPA) THEN
               WRITE(LUPRI,'(3(/5X,2A,F15.8))')
     *         PRPLBL,' inactive part:',ONEINA,
     *         PRPLBL,' active part  :',ONEACT,
     *         PRPLBL,' total        :',ONETOT
            ELSE
               WRITE(LUPRI,'(3(/5X,2A,F15.8))')
     *         PRPLBL,' Hartree-Fock part:',ONEINA,
     *         PRPLBL,' second order part:',ONESEC,
     *         PRPLBL,' total 2. order   :',ONETOT
            END IF
         ELSE
            IF (.NOT. SOPPA) THEN
               WRITE(LUPRI,'(3(/5X,2A,1P,D15.8))')
     *         PRPLBL,' inactive part:',ONEINA,
     *         PRPLBL,' active part  :',ONEACT,
     *         PRPLBL,' total        :',ONETOT
            ELSE
               WRITE(LUPRI,'(3(/5X,2A,1P,D15.8))')
     *         PRPLBL,' Hartree-Fock part:',ONEINA,
     *         PRPLBL,' second order part:',ONESEC,
     *         PRPLBL,' total 2. order   :',ONETOT
            END IF
         END IF
      END IF
C
C END OF PRPONE
C
      RETURN
      END
C  /* Deck mixs0 */
      SUBROUTINE MIXS0(CMO,UDV,WRK,LWRK)
C
C CALCULATE S(0) IN MIXED REPRESENTATION
C
C     NELECTRON = abs ( <0| [p,r] |0> )
C
#include "implicit.h"
C
#include "priunit.h"
#include "infs0.h"
#include "inforb.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
C
      DIMENSION CMO(*),UDV(*),WRK(*)
C
      IF (NOS0MX) RETURN
      CALL HEADER(' ** S(0) sum rule in mixed representation **',0)
C
C ALLOCATE WORK SPACE
C
      KDLEN  = 1
      KDVEL  = KDLEN  + NORBT*NORBT
      KWRK1  = KDVEL  + NORBT*NORBT
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('MIXS0',KWRK1-1,LWRK)
      DO 100 JSYMOP = 1,NSYM
         NOP = NGPS0(JSYMOP)
      DO 100 IOP    = 1,NOP
C
C Read AO property integrals and transform to MO basis.
C
         KSYMOP = JSYMOP
         CALL PRPGET (LBLS0(JSYMOP,IOP,1),
     *                CMO,WRK(KDLEN),ANTSYM,WRK(KWRK1),LWRK1)
         CALL PRPGET (LBLS0(JSYMOP,IOP,2),
     *                CMO,WRK(KDVEL),ANTSYM,WRK(KWRK1),LWRK1)
C        CALL PRPGET (WORD,CMO,PRPMO,WRK,LWRK)
C
C Calculate [ d/dx , x ]
C
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &              WRK(KDVEL),NORBT,
     &              WRK(KDLEN),NORBT,0.D0,
     &              WRK(KWRK1),NORBT)
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,-1.D0,
     &              WRK(KDLEN),NORBT,
     &              WRK(KDVEL),NORBT,1.D0,
     &              WRK(KWRK1),NORBT)
C
C Print atomic and molecular property integrals, if desired
C
         IF (IPRRSP.GT.125) THEN
            WRITE (LUPRI,'(/3A)')' [d/dq,q] PROPERTY INTEGRALS:',
     *         LBLS0(JSYMOP,IOP,2), LBLS0(JSYMOP,IOP,1)
            CALL OUTPUT(WRK(KWRK1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
C
C CALCULATE AVERAGE VALUE OF ONE ELECTRON OPERATOR
C
         WRITE(LUPRI,'(/5A,I2,/A)')
     *   ' Symmetry of ',
     *   LBLS0(JSYMOP,IOP,1),' and ', LBLS0(JSYMOP,IOP,2),
     *   ' is',JSYMOP,
     *   ' number of electrons from S(0) mixed :'
         KSYMOP = 1
         IPRONE = 0
         CALL PRPONE(WRK(KWRK1),UDV,ONETOT,IPRONE,'S(0) MIXED')
 100  CONTINUE
      RETURN
      END
