!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C infave.h: N1AVTOT,LBL1AV,AVE1VAL, ...
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
     & ' Calculation of electronic one-electron expectation values ',-1)
        WRITE (LUPRI,'(A/A)')
     &  ' (Note that to get e.g. a dipole moment you must multiply the',
     &  '  electronic number by -1 and add the nuclear contribution.)'
C
        DO IOP    = 1,N1AVTOT
          CALL PRP1AVE(LBL1AV(IOP),AVE1VAL(IOP),
     &       CMO,UDV,WRK,LWRK,IPRRSP)
        END DO
        CALL RSP_AVEDIP()
      END IF
C
      IF (N2AVTOT.GT.0) THEN
        CALL HEADER(
     & ' Calculation of electronic two-electron expectation values ',-1)
        DO IOP    = 1,N2AVTOT
          CALL PRP2AVE(LBL2AV(1,IOP),AVE2VAL(IOP),
     &       CMO,UDV,PVX,WRK,LWRK,IPRRSP)
        END DO
      END IF
C
      CALL QEXIT('AVEPRP')
      RETURN
      END
      SUBROUTINE RSP_AVEDIP()
C
C Print out of the dipole moment
C
#include "implicit.h"
#include "priunit.h"
C infave.h: LBL1AV(:), AVE1VAL(:)
C dipole.h: DIPMN(:)
#include "rspprp.h"
#include "infave.h"
#include "infrsp.h"
#include "mxcent.h"
#include "dipole.h"

      REAL*8 DIP(3), DIPMOM

      DO IOP = 1, N1AVTOT
        IF (LBL1AV(IOP).EQ.'XDIPLEN') THEN
          DIP(1) = DIPMN(1) - AVE1VAL(IOP)
        ELSE IF (LBL1AV(IOP).EQ.'YDIPLEN') THEN
          DIP(2) = DIPMN(2) - AVE1VAL(IOP)
        ELSE IF (LBL1AV(IOP).EQ.'ZDIPLEN') THEN
          DIP(3) = DIPMN(3) - AVE1VAL(IOP)
        END IF
      END DO
        CALL HEADER('Dipole moment - electronic and nuclear',-1)
        DIPMOM = SQRT(DIP(1)*DIP(1) + DIP(2)*DIP(2) +
     &       DIP(3)*DIP(3))
        WRITE (LUPRI,'(17X,A,15X,A,10X,A/3X,3F19.6)')
     *         'au','Debye','C m (/(10**-30)',
     *         DIPMOM, DEBYE*DIPMOM, DIPSI*DIPMOM
        CALL HEADER(
     $       'Dipole moment components',-1)
        CALL DP0PRI(DIP)

      I = 0
      DO IOP = 1, N1AVTOT
        IF (LBL1AV(IOP).EQ.'XLFDIPLN') THEN
          DIP(1) = DIPMN(1) - AVE1VAL(IOP)
          I = I + 1
        ELSE IF (LBL1AV(IOP).EQ.'YLFDIPLN') THEN
          DIP(2) = DIPMN(2) - AVE1VAL(IOP)
          I = I + 1
        ELSE IF (LBL1AV(IOP).EQ.'ZLFDIPLN') THEN
          DIP(3) = DIPMN(3) - AVE1VAL(IOP)
          I = I + 1
        END IF
      END DO
      IF (I.EQ.3) THEN
        CALL HEADER('Local-field corrected Dipole moment',-1)
        DIPMOM = SQRT(DIP(1)*DIP(1) + DIP(2)*DIP(2) +
     &       DIP(3)*DIP(3))
        WRITE (LUPRI,'(17X,A,15X,A,10X,A/3X,3F19.6)')
     *         'au','Debye','C m (/(10**-30)',
     *         DIPMOM, DEBYE*DIPMOM, DIPSI*DIPMOM
        CALL HEADER(
     $       'Local-field corrected dipole moment components',-1)
        CALL DP0PRI(DIP)
      END IF
      END
C  /* Deck prp1ave */
      SUBROUTINE PRP1AVE(LBL_1AVE,AVEVAL,CMO,UDV,WRK,LWRK,IPRINT)
C
C Calculate average value (expectation value) of property
C with label LBL_1AVE on AOPROPER
C
C Revised March 2003, Hans Joergen Aa. Jensen.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER ( D0 = 0.0D0 )
C
#include "inforb.h"
#include "wrkrsp.h"
C
      CHARACTER*8 LBL_1AVE
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
C Read AO property integrals for LBL_1AVE and check for antisymmetry
C and for spatial symmetry (as identified by abs largest element).
C
      KSYMP = -1
      CALL PRPGET (LBL_1AVE,CMO,WRK(KPMO),KSYMP,ANTSYM,
     &             WRK(KWRK1),LWRK1,IPRINT)
C     CALL PRPGET (WORD,CMO,PRPMO,KSYMP,ANTSYM,WRK,LWRK,IPRINT)
C
C Calculate average value of one electron operator.
C
      IF (ANTSYM .LT. D0) THEN
         AVEVAL = D0
         WRITE(LUPRI,'(/3A)')
     &      ' Property "',LBL_1AVE,'" is antisymmetric'//
     &      ' and has therefore average value ZERO'
      ELSE IF ( KSYMP .EQ. 1 ) THEN
         IPRONE = 0
         KSYMOP = 1
         CALL PRPONE(WRK(KPMO),UDV,AVEVAL,IPRONE,LBL_1AVE)
      ELSE
         AVEVAL = D0
         WRITE(LUPRI,'(/3A,I3,A)')
     &      ' Property "',LBL_1AVE,'" has symmetry',KSYMP,
     &      ' and has therefore average value ZERO.'
      END IF
      CALL QEXIT('PRP1AVE')
      RETURN
      END
C  /* Deck prp2ave */
      SUBROUTINE PRP2AVE(LBL_2AVE,AVE2VAL,CMO,UDV,PVX,WRK,LWRK,IPRINT)
C
C Calculate average value (expectation value) of 2-electron property
C built from two one-electron operators
C with labels LBL_2AVE(1:2) on AOPROPER
C
C Augus 2011, Hans Joergen Aa. Jensen.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      PARAMETER ( D0 = 0.0D0 )
C
#include "inforb.h"
#include "wrkrsp.h"
C
      CHARACTER*8 LBL_2AVE(2)
      DIMENSION   CMO(*),UDV(*),PVX(*),WRK(*)
      DIMENSION   ANTSYM(2), KSYMP(2)
C
      CALL QENTER('PRP2AVE')
C
C Allocate work space
C
      KFRSAV = 1
      KFREE  = 1
      LFREE  = LWRK
      CALL MEMGET('REAL',KP1MO,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KP2MO,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KP12MO,N2ORBX,WRK,KFREE,LFREE)
C
C Read AO property integrals for LBL_2AVE(1:2) and check for antisymmetry
C and for spatial symmetry (as identified by abs largest element).
C
      KSYMP(1:2) = -1
      CALL PRPGET (LBL_2AVE(1),CMO,WRK(KP1MO),KSYMP(1),ANTSYM(1),
     &             WRK(KFREE),LFREE,IPRINT)
      CALL PRPGET (LBL_2AVE(2),CMO,WRK(KP2MO),KSYMP(2),ANTSYM(2),
     &             WRK(KFREE),LFREE,IPRINT)
C     CALL PRPGET (WORD,CMO,PRPMO,KSYMP,ANTSYM,WRK,LWRK,IPRINT)

C     P12(p,q) = sum(r) P1(p,r) * P2(r,q)
      CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &           WRK(KP1MO),NORBT,
     &           WRK(KP2MO),NORBT,0.D0,
     &           WRK(KP12MO),NORBT)

      ANTSYM_tot = ANTSYM(1) * ANTSYM(2)
      KSYMP_tot  = MULD2H( KSYMP(1) , KSYMP(2) )

!     IF (IPRINT .GT. 0) THEN
         WRITE(LUPRI,'(/A/)')
     &      '============== Next operator product ====================='
         WRITE(LUPRI,*) 'Operator 1: ',LBL_2AVE(1), KSYMP(1), ANTSYM(1)
         WRITE(LUPRI,*) 'Operator 2: ',LBL_2AVE(2), KSYMP(2), ANTSYM(2)
         WRITE(LUPRI,*) 'Product   :         ', KSYMP_tot, ANTSYM_tot
!     END IF
C
C Calculate average value of one electron operator.
C
      KSYMOP = 1 ! symmetry of UDV and PVX
      IF (ANTSYM_tot .LT. D0) THEN
         AVE2VAL = D0
         WRITE(LUPRI,'(/5A)')
     &      ' Property "',LBL_2AVE(1),' x ',LBL_2AVE(2),
     &      '" is antisymmetric and has therefore average value ZERO'
      ELSE IF ( KSYMP_tot .NE. KSYMOP ) THEN
         AVE2VAL = D0
         WRITE(LUPRI,'(/5A,I2,A)')
     &      ' Property "',LBL_2AVE(1),' x ',LBL_2AVE(2),
     &      '" has symmetry',KSYMP_tot,
     &      ' and has therefore average value ZERO'
      ELSE
          IPRONE = 0
          CALL PRP2ONE(WRK(KP1MO),WRK(KP2MO),WRK(KP12MO),UDV,PVX,
     &                 AVE2VAL,IPRONE,LBL_2AVE)
      END IF
      CALL QEXIT('PRP2AVE')
      RETURN
      END
C  =====================================================================
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
      IF (SOPPA) THEN
         CALL ONEMP2(PRPMO,UDV,ONESEC)
      ELSE
         ONESEC = D0
      END IF
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
C  =====================================================================
      SUBROUTINE PRP2ONE(PRP1_MO,PRP2_MO,PRP12_MO,
     &                   UDV,PVX,AVE2VAL,IPR2ONE,PRPLBL)
C
C Calculate  AVERAGE VALUE of the product of the two ONE ELECTRON OPERATORS
C in PRP1_MO and PRP2_MO.
C Examples:
C   XDIPVEL and XDIPVEL for mass polarization
C   XANGMOM and XANGMOM for angular momentum expectation value (L**2)
C
C Written Aug. 2011 Hans Jørgen Aa. Jensen.
C
#include "implicit.h"
C
      CHARACTER*(*) PRPLBL(2)
C
      DIMENSION PRP1_MO(NORBT,NORBT),PRP2_MO(NORBT,NORBT)
      DIMENSION PRP12_MO(NORBT,NORBT)
      DIMENSION UDV(NASHDI,NASHDI), PVX(NASHDI,NASHDI,NASHDI,NASHDI)
      DIMENSION AVE1_PRP12(2), AVE1_PRP1(2), AVE1_PRP2(2) ! 1 - inactive, 2 - active
C
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
#include "infpri.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infind.h"
#include "infdim.h"
C
      PARAMETER ( D2 = 2.0D0 , D0 = 0.0D0 )
      PARAMETER ( BIGLIM = 100000.0D0, SMLLIM = 0.01D0 )
C
      IPRORB = IPR2ONE
      IF (NOCCT .GT. 12) IPRORB = MAX(IPRORB,2)
C     hjaaj Mar 2003: default IPRRSP is 2, do not print individual
C     contributions by default if more than 12 occupied orbitals.
!     IF (IPRRSP.GT.IPRORB) THEN
!        WRITE(LUPRI,'(/A/4A/)')
!    &      ' *** Individual non-zero orbital contributions',
!    &      ' *** to the expectation value for property ',PRPLBL(1:2),' :'
!     END IF
      IF (TRPLET) THEN
         WRITE(LUPRI,'(//A/A)')
     &      'WARNING: PRP2ONE called for triplet (not implemented)',
     &      'WARNING: AVE2VAL therefore arbitrarily set to zero'
         AVE2VAL = 0.0D0
         RETURN
      END IF
C
      IF (IPRRSP .GT. 25) THEN
         WRITE(LUPRI,'(/A)') ' Property matrix 1'
         CALL PRMGN(N2ORBX,PRP1_MO,1,4,LUPRI)
         CALL OUTPUT(PRP1_MO,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         WRITE(LUPRI,'(/A)') ' Property matrix 2'
         CALL PRMGN(N2ORBX,PRP2_MO,1,4,LUPRI)
         CALL OUTPUT(PRP2_MO,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         WRITE(LUPRI,'(/A)') ' Property matrix 1 * property matrix 2'
         CALL PRMGN(N2ORBX,PRP12_MO,1,4,LUPRI)
         CALL OUTPUT(PRP12_MO,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         WRITE(LUPRI,'(/A)') ' DV matrix'
         CALL PRMGN(N2ASHX,UDV,1,4,LUPRI)
         CALL OUTPUT(UDV,1,NASHT,1,NASHT,NASHT,NASHT,-1,LUPRI)
         WRITE(LUPRI,'(/A)') ' PV matrix'
         CALL PRMGN(N2ASHX*N2ASHX,PVX,1,4,LUPRI)
         CALL OUTPUT(PVX,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,-1,LUPRI)
      END IF
C
      AVE1_PRP12(1:2) = 0.0D0 !  sum(i) P12(i,i)*2; sum(u,v) P12(u,v)*D(u,v)
      AVE1_PRP1(1:2)  = 0.0D0 !  sum(i) P1(i,i)*2; sum(u,v) P1(u,v)*D(u,v)
      AVE1_PRP2(1:2)  = 0.0D0 !  sum(i) P2(i,i)*2; sum(u,v) P2(u,v)*D(u,v)
      AVE_PRP12_DC    = 0.0D0 ! -sum(i,j)   P1(j,i)*P2(i,j)*2
      AVE_PRP12_DV    = 0.0D0 ! -sum(i,u,v) P1(u,i)*P2(i,v)*DV(u,v)
      AVE_PRP12_PV    = 0.0D0 !  sum(u,v,x,y) PV(uv,xy) * P1(u,v) * P2(x,y)
      DO J = 1,NISHT
         J_INAC = ISX(J)
         AVE1_PRP1(1)  = AVE1_PRP1(1)  + PRP1_MO(J_INAC,J_INAC)*D2
         AVE1_PRP2(1)  = AVE1_PRP2(1)  + PRP2_MO(J_INAC,J_INAC)*D2
         AVE1_PRP12(1) = AVE1_PRP12(1) + PRP12_MO(J_INAC,J_INAC)*D2
         DO I = 1, NISHT
            I_INAC = ISX(I)
            AVE_PRP12_DC = AVE_PRP12_DC -
     &         PRP1_MO(I_INAC,J_INAC)*PRP2_MO(J_INAC,I_INAC)
         END DO
      END DO
      AVE_PRP12_DC = D2*AVE_PRP12_DC

      DO J = 1, NASHT
         J_ACT = ISX(NISHT+J)
         DO I = 1, NASHT
            I_ACT = ISX(NISHT+I)
            AVE1_PRP1(2) = AVE1_PRP1(2) +
     &         PRP1_MO(I_ACT,J_ACT) * UDV(I,J)
            AVE1_PRP2(2) = AVE1_PRP2(2) +
     &         PRP2_MO(I_ACT,J_ACT) * UDV(I,J)
            AVE1_PRP12(2) = AVE1_PRP12(2) +
     &         PRP12_MO(I_ACT,J_ACT) * UDV(I,J)
            DO L = 1, NISHT
               L_INAC = ISX(L)
               AVE_PRP12_DV = AVE_PRP12_DV -
     &         PRP1_MO(I_ACT,L_INAC)*PRP2_MO(L_INAC,J_ACT) * UDV(I,J)
            END DO

            AVE_PRP2_PV_IJ = 0.0D0
            DO L = 1, NASHT
               L_ACT = ISX(NISHT+L)
               DO K = 1, NASHT
                  K_ACT = ISX(NISHT+K)
                  AVE_PRP2_PV_IJ = AVE_PRP2_PV_IJ +
     &               PRP2_MO(K_ACT,L_ACT)*PVX(K,L,J,I)
               END DO
            END DO
            AVE_PRP12_PV = AVE_PRP12_PV +
     &           AVE_PRP2_PV_IJ*PRP1_MO(I_ACT,J_ACT)
         END DO
      END DO

      TERM_1 = AVE1_PRP1(1)*AVE1_PRP2(1)
      TERM_2 = AVE1_PRP1(1)*AVE1_PRP2(2) + AVE1_PRP1(2)*AVE1_PRP2(1)
      TERM_TOT = TERM_1 + TERM_2
     &         + AVE_PRP12_DC + AVE_PRP12_DV + AVE_PRP12_PV

!     IF (IPRINT .GT. ??) THEN
         WRITE(LUPRI,'(//A/)') ' --- individual contributions: --'
         WRITE(LUPRI,*) 'AVE1_PRP1  DC, DV, DTOT',AVE1_PRP1(1:2),
     &      AVE1_PRP1(1)+AVE1_PRP1(2)
         WRITE(LUPRI,*) 'AVE1_PRP2  DC, DV, DTOT',AVE1_PRP2(1:2),
     &      AVE1_PRP2(1)+AVE1_PRP2(2)
         WRITE(LUPRI,*) 'AVE1_PRP12 DC, DV, DTOT',AVE1_PRP12(1:2),
     &      AVE1_PRP12(1)+AVE1_PRP12(2)

         WRITE(LUPRI,*) 'AVE_PRP12_DC',AVE_PRP12_DC
         WRITE(LUPRI,*) 'AVE_PRP12_DV',AVE_PRP12_DV
         WRITE(LUPRI,*) 'AVE_PRP12_PV',AVE_PRP12_PV
         WRITE(LUPRI,*) 'AVE1_PRP1*AVE1_PRP2: DC*DC and DC*DV+DV*DC',
     &      TERM_1,TERM_2
         WRITE(LUPRI,*) 'PRP1*PRP2 (DC*DC - exch)',TERM_1 + AVE_PRP12_DC
!     END IF

      WRITE(LUPRI,'(/5A,T40,F20.10)')
     &   '@ < ',PRPLBL(1),'(1) * ',PRPLBL(2),'(2) > =',TERM_TOT
      WRITE(LUPRI,'( 5A,T40,F20.10)')
     &   '@ < ',PRPLBL(1),'(1) * ',PRPLBL(2),'(1) > =',
     &   AVE1_PRP12(1)+AVE1_PRP12(2)
      WRITE(LUPRI,'( 3A,T40,F20.10)')
     &   '@ < ',PRPLBL(1),'(1) > =',AVE1_PRP1(1) +AVE1_PRP1(2)
      WRITE(LUPRI,'( 3A,T40,F20.10)')
     &   '@ < ',PRPLBL(2),'(2) > =',AVE1_PRP2(1) +AVE1_PRP2(2)
      AVE2VAL = TERM_TOT
      RETURN
     
#ifdef MAYBE_USE_IN_THE_FUTURE
C
      IF (SOPPA) THEN

         ONESEC = 0.0D0
         WRITE (LUPRI,'(//A/A/A)')
     &      'WARNING, SOPPA not implemented here',
     &      'WARNING, second order part set to zero',
     &      'WARNING, expectation value is thus only '//
     &         'correct to first order'

         ONETOT = ONETOT + ONESEC

      END IF
C
      ONETOT = ONEINA + ONEACT
      IF (IPRRSP.GE.IPR2ONE) THEN
         IF (ABS(ONETOT) .GT. SMLLIM .AND. ABS(ONETOT) .LT. BIGLIM)
     *                                                       THEN
            IF (.NOT. SOPPA) THEN
               WRITE(LUPRI,'(3(/5X,3A,F15.8))')
     *         PRPLBL,' inactive part:',ONEINA,
     *         PRPLBL,' active part  :',ONEACT,
     *         PRPLBL,' total        :',ONETOT
            ELSE
               WRITE(LUPRI,'(3(/5X,3A,F15.8))')
     *         PRPLBL,' Hartree-Fock part:',ONEINA,
     *         PRPLBL,' second order part:',ONESEC,
     *         PRPLBL,' total 2. order   :',ONETOT
            END IF
         ELSE
            IF (.NOT. SOPPA) THEN
               WRITE(LUPRI,'(3(/5X,3A,1P,D15.8))')
     *         PRPLBL,' inactive part:',ONEINA,
     *         PRPLBL,' active part  :',ONEACT,
     *         PRPLBL,' total        :',ONETOT
            ELSE
               WRITE(LUPRI,'(3(/5X,3A,1P,D15.8))')
     *         PRPLBL,' Hartree-Fock part:',ONEINA,
     *         PRPLBL,' second order part:',ONESEC,
     *         PRPLBL,' total 2. order   :',ONETOT
            END IF
         END IF
      END IF
#endif
C
C END OF PRP2ONE
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
      CALL HEADER(' ** S(0) sum rule in mixed representation ** ',-1)
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
         IF (IPRRSP.GT.25) THEN
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
