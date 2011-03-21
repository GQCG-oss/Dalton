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
#include "single.h"
C
C  /* Deck aveprp */
      SUBROUTINE AVEPRP(CMO,UDV,PVX,WRK,LWRK)
C
C Calculate average value (expectation value) of properties
C specified with .PROPAV under **RESPONS
C
C Revised Mar 2003 Hans Joergen Aa. Jensen
C
#include "implicit.h"
#include "priunit.h"
C
C infrsp.h: IPRRSP
#include "rspprp.h"
#include "infave.h"
#include "infrsp.h"
C
      DIMENSION CMO(*),UDV(*),PVX(*),WRK(*)
C
      CALL QENTER('AVEPRP')
C
      IF (N1AVTOT.GT.0) THEN
        CALL HEADER(
     &  'Calculation of electronic one-electron expectation values',0)
        WRITE (LUPRI,'(A/A)')
     &  ' (Note that to get e.g. a dipole moment you must multiply the',
     &  '  electronic number by -1 and add the nuclear contribution.)'
C
        DO IOP    = 1,N1AVTOT
          CALL PRP1AVE(LBL1AV(IOP),AVEVAL,CMO,UDV,WRK,LWRK,IPRRSP)
        END DO
      END IF
C
      IF (N2AVTOT.GT.0) THEN
        CALL HEADER(
     &  'Calculation of electronic two-electron expectation values',0)
        DO IOP    = 1,N2AVTOT
          write (lupri,'(2X,A)') (LBL2AV(i,IOP),i=1,2)
chj       CALL PRP2AVE(LBL2AV(1,IOP),AVEVAL,CMO,UDV,PVX,WRK,LWRK,IPRRSP)
        END DO
      END IF
C
      CALL QEXIT('AVEPRP')
      RETURN
      END
C  /* Deck prp1ave */
      SUBROUTINE PRP1AVE(LBLAVE,AVEVAL,CMO,UDV,WRK,LWRK,IPRINT)
C
C Calculate average value (expectation value) of property
C with label 'LBLAVE' on AOPROPER
C
C Revised March 2003, Hans Joergen Aa. Jensen.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER ( D0 = 0.0D0 ,  CKMXPR = 1.0D-12 )
C
#include "inforb.h"
#include "wrkrsp.h"
C
      CHARACTER*8 LBLAVE
      DIMENSION CMO(*),UDV(*),WRK(*)
C
      CALL QENTER('PRP1AVE')
C
C Allocate work space
C
      KPMO   = 1
      KWRK1  = KPMO  + N2ORBX
      LWRK1  = LWRK  - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('PRP1AVE',KWRK1-1,LWRK)
C
C Read AO property integrals for LBLAVE and check for antisymmetry
C and for spatial symmetry (as identified by abs largest element).
C
      KSYMP = -1
      CALL PRPGET (LBLAVE,CMO,WRK(KPMO),KSYMP,ANTSYM,
     &             WRK(KWRK1),LWRK1,IPRINT)
C     CALL PRPGET (WORD,CMO,PRPMO,KSYMP,ANTSYM,WRK,LWRK,IPRINT)
C
C Calculate average value of one electron operator.
C
      IF (ANTSYM .GT. D0) THEN
         AVEVAL = D0
         WRITE(LUPRI,'(/3A)')
     &      ' Property "',LBLAVE,'" is antisymmetric'//
     &      ' and has therefore average value ZERO'
      ELSE IF ( KSYMP .EQ. 1 ) THEN
          IPRONE = 0
          KSYMOP = 1
          CALL PRPONE(WRK(KPMO),UDV,AVEVAL,IPRONE,LBLAVE)
      ELSE
         AVEVAL = D0
         WRITE(LUPRI,'(/3A,I3,A)')
     &      ' Property "',LBLAVE,'" has symmetry',KSYMP,
     &      ' and has therefore average value ZERO.'
      END IF
      CALL QEXIT('PRP1AVE')
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
         DO 60 JA = 1,NASHJ
            DO 70 IA = 1,NASHI
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
      CALL QENTER('MIXS0 ')
      CALL HEADER(' ** S(0) sum rule in mixed representation **',0)
C
C ALLOCATE WORK SPACE
C
      KDLEN  = 1
      KDVEL  = KDLEN  + N2ORBX
      KWRK1  = KDVEL  + N2ORBX
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('MIXS0',KWRK1-1,LWRK)
      DO 100 JSYMOP = 1,NSYM
         NOP = NGPS0(JSYMOP)
      DO 100 IOP    = 1,NOP
C
C Read AO property integrals and transform to MO basis.
C
         KSYMP1 = -1
         KSYMP2 = -1
         CALL PRPGET (LBLS0(JSYMOP,IOP,1),
     *                CMO,WRK(KDLEN),KSYMP1,ASM,WRK(KWRK1),LWRK1,IPRRSP)
         CALL PRPGET (LBLS0(JSYMOP,IOP,2),
     *                CMO,WRK(KDVEL),KSYMP2,ASM,WRK(KWRK1),LWRK1,IPRRSP)
C        CALL PRPGET (WORD,CMO,PRPMO,KSYMP,ANTSYM,WRK,LWRK,IPRINT)
         IF (KSYMP1 .NE. JSYMOP .OR. KSYMP2 .NE. JSYMOP)
     &   CALL QUIT('MIXS0: unexpeted symmetry of property matrix')
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
      CALL QEXIT('MIXS0 ')
      RETURN
      END
