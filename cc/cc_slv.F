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
C  /* Deck cc_slv */
       SUBROUTINE CC_SLV(AODEN,ETLM,DIELCONV,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster solvent effects.
C
C              CCS(CIS/HF)(nci), MP2(nci), CC2(nci), CCSD, CC3(nci), MCC2(nci)
C
C     SLV98,OC
C     Ove Christiansen, April 1998.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      LOGICAL   DIELCONV
      DIMENSION WORK(LWORK),AODEN(*),ETLM(NLMCU,2)
      CHARACTER MODEL*10
      CHARACTER MODELPRI*4
C
C-----------
C     Header 
C-----------
C
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*    Output from coupled cluster solvent program         '//
     *'         *'
      WRITE (LUPRI,'(1X,A,/)')
     *'*********************************************************'//
     *'**********'
C
          WRITE(LUPRI,'(/,1X,A)')
     *'+=====================================================+'
          WRITE(LUPRI,'(1X,A)')
     *'| Lmax     Cavity(a.u.)    Epsilon_st      Epsilon_op |'
          WRITE(LUPRI,'(1X,A)')
     *'+-----------------------------------------------------+'
          WRITE(LUPRI,'(1X,A,I3,4X,F13.8,2X,F13.8,3X,F13.8,A)')
     *    '| ',LMAXCU,RCAVCU,EPSTCU,EPOPCU,' |'
          WRITE(LUPRI,'(1X,A)')
     * '+=====================================================+'
C
      MODEL = 'CCSD'
      IF (CC2.AND.(.NOT.MCC2)) THEN
         CALL AROUND('Coupled Cluster model is: CC2')
         MODEL = 'CC2'
         MODELPRI = ' CC2'
      ENDIF
      IF (MCC2) THEN
         CALL AROUND('Coupled Cluster model is: MCC2')
         MODEL = 'MCC2'
         MODELPRI = 'MCC2'
      ENDIF
      IF (MP2) THEN
         CALL AROUND('Model is second order pert. theory: MP2 ')
         MODEL = 'MP2'
         MODELPRI = ' MP2'
      ENDIF
      IF (CCS.AND.(.NOT.CIS)) THEN
         CALL AROUND('Coupled Cluster model is: CCS')
         MODEL = 'CCS'
         MODELPRI = ' CCS'
      ENDIF
      IF (CCS.AND.CIS) THEN
         CALL AROUND('CIS model in use ')
         MODEL = 'CCS'
         MODELPRI = ' CIS'
      ENDIF
      IF (CCD) THEN
         CALL AROUND('Coupled Cluster model is: CCD')
         MODEL = 'CCD'
         MODELPRI = ' CCD'
      ENDIF
      IF (CC3  ) THEN
         CALL AROUND('Coupled Cluster model is: CC3')
         MODEL = 'CC3'
         MODELPRI = ' CC3'
         CALL QUIT('CC3 first order properties not implemented')
      ENDIF
      IF (CC1A) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1a')
         MODEL = 'CCSDT-1a'
         CALL QUIT( 'CCSDT-1a first order properties not implemented')
      ENDIF
      IF (CC1B) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1b')
         MODEL = 'CCSDT-1b'
         CALL QUIT( 'CCSDT-1b first order properties not implemented')
      ENDIF
      IF (CCSD) THEN
         CALL AROUND('Coupled Cluster model is: CCSD')
         MODEL = 'CCSD'
         MODELPRI = 'CCSD'
      ENDIF
C
c     IF ((.NOT. CCSD).OR.(.NOT.CC2)) THEN
c      CALL QUIT( ' Solvent not implemented for other than CC2 and CCSD')
c     ENDIF
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_SLV-1: Workspace:',LWORK
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SLV98,OC
C     Solvent section
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C------------------------------------------------------
C     Calculate the solvent contribution to the energy.
C------------------------------------------------------
C
      CALL CC_SLVE(ECCSLCON,AODEN,ETLM,WORK,LWORK)
C
C--------------------------------------------------------------
C     Calculate new solvent energy.
C--------------------------------------------------------------
C
      ECCCU = ECCGRS + ECCSLCON
      IF (ABS(ECCCU-ECCPR).LT.CVGESOL) LSLECVG = .TRUE.
      WRITE(LUPRI,'(1X,A,I3,A,F20.10)')
     *    'Solvent energy contribution in iteration',ICCSLIT,': ',
     *     ECCSLCON
      WRITE(LUPRI,'(1X,A,F20.10)')
     *    'CC energy in the  current solvent iteration: ',ECCCU  
      WRITE(LUPRI,'(1X,A,F20.10)')
     *    'CC energy in the previous solvent iteration: ',ECCPR  
      WRITE(LUPRI,*)'LSLECVG: ',LSLECVG
      WRITE(LUPRI,*)
     *' Change in Total energy in this solvent it.:',
     * ECCCU - ECCPR 
      ECCPR   = ECCCU  
      DIELCONV = .FALSE.
      IF (LSLECVG.AND.LSLTCVG.AND.LSLLCVG) THEN
        DIELCONV = .TRUE.
        WRITE(LUPRI,'(/,1X,A,I3,A)')
     *'Coupled cluster solvent equations are converged in ',ICCSLIT,
     *' solvent iterations'
        WRITE(LUPRI,'(/,1X,A8,A,F30.16,/)')
     *  MODEL,'converged energy in solvent:  ',ECCCU
        WRITE(LUPRI,'(/,1X,A8,A,F30.16,/)')
     *  MODEL,'solvation energy:             ',ECCSLCON
        WRITE(LURES,'(/,1X,A,I2,A,F8.4,A,F8.4,A,F8.4,A)')
     *   'Solvent: L_max=',LMAXCU,', R_cav=',RCAVCU,', Eps_st =',EPSTCU,
     *   ', Eps_op =',EPOPCU,': '
        WRITE(LURES,'(12X,A8,A,F20.10)')
     *  MODEL,' Total energy:           ',ECCCU
        WRITE(LURES,'(12X,A8,A,F20.10,/)')
     *  MODEL,' Solvation energy:       ',ECCSLCON
        ECCGRS = ECCCU    
      ELSE 
        ICCSLIT = ICCSLIT + 1
        IF (ICCSLIT.GT.MXCCSLIT) THEN
          WRITE(LUPRI,*) 'Maximum number of solvent iterations',
     *                MXCCSLIT,' is reached'
          CALL QUIT( 'Maximum number of solvent iterations reached')
        ENDIF
      ENDIF
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     SLV98,OC
C     End of solvent section.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*        End of coupled cluster solvent program          '//
     *'         *'
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
C
      END
C  /* Deck cc_slve */
       SUBROUTINE CC_SLVE(ESOLT,AODEN,ETLM,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Calculation of Coupled Cluster solvent energy. 
C              for given solvent defined by RCAVCU,EPSTCU,LMAXCU.
C
C              Based on SOLGRD from Sirius.
C
C     SLV98,OC
C     Ove Christiansen, April 1998.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
#include <thrzer.h>
C
      DIMENSION WORK(LWORK),AODEN(*),ETLM(NLMCU,2)
      CHARACTER MODEL*10
C
      LOGICAL     FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C
      DOUBLE PRECISION GL
      EXTERNAL GL
C
C-----------------------
C     Check g_l factors.
C-----------------------
C
      IF (IPRINT.GT.10) THEN
         WRITE(LUPRI,*) 'L  G_l(Eps) for Eps,Rcav=',EPSTCU,RCAVCU
         DO L = 0, LMAXCU
            WRITE(LUPRI,*) L,GL(L,EPSTCU,RCAVCU)
         ENDDO
      ENDIF
C
C-------------------------------------------------------
C     Check integrals and readin nuclear contributions.
C-------------------------------------------------------
C
      LUCSOL = -1
      CALL GPOPEN(LUCSOL,'AOSOLINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCSOL)
      CALL MOLLAB('SOLVRLM ',LUCSOL,LUERR)
C
      IF (FIRST) THEN
         READ (LUCSOL) LMAXSS, LMTOT, NNNBAS
         IF (LMAXSS .LT. LMAXCU) THEN
            WRITE (LUERR,'(//2A,2(/A,I5))') ' >>> CC_SLVE ERROR,',
     *      ' insufficient number of intgrals on LUCSOL',
     *      ' l max from CC     input :',LMAXCU,
     *      ' l max from LUCSOL  file  :',LMAXSS
            CALL QUIT('CC_SLVE: lmax on LUCSOL is too small')
         END IF
         IF ((LMAXSS+1)**2 .NE. LMTOT) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> CC_SLVE ERROR,',
     *      ' LUCSOL file info inconsistent',
     *      ' l_max               :',LMAXSS,
     *      ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *      ' LMTOT               :',LMTOT
            CALL QUIT('CC_SLVE: LUCSOL info not internally consistent')
         END IF
         IF (NNNBAS .NE. NBAST) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> CC_SLVE ERROR,',
     *      ' LUCSOL file info inconsistent with SIRIUS input',
     *      ' NBAST - LUCSOL       :',NNNBAS,
     *      ' NBAST - SIRIUS      :',NBAST
            CALL QUIT('CC_SLVE: LUCSOL info not '//
     &                'consistent with SIRIUS input.')
         END IF
         FIRST = .FALSE.
      ELSE
         READ (LUCSOL)
      END IF
      CALL READT(LUCSOL,NLMCU,ETLM(1,2))
C
C-------------------------------------
C     Print out nuclear contributions.
C-------------------------------------
C
      IF (IPRINT .GE. 10) THEN
         WRITE(LUPRI,'(/A/)')
     *      ' l, m, Tn(lm) - the nuclear contributions :'
         LM = 0
         DO 220 L = 0,LMAXCU 
            DO 210 M = -L,L
               LM = LM + 1
               WRITE(LUPRI,'(2I5,F15.10)') L,M,ETLM(LM,2)
  210       CONTINUE
            WRITE(LUPRI,'()')
  220    CONTINUE
      END IF
C
C-----------------------------------------
C     Calculate electronic conctributions.
C     1. Loop L
C     2. Get symmetries of Tlm
C     3. Loop M
C-----------------------------------------
C
      LM = 0
      DO 520 L = 0,LMAXCU 
         READ (LUCSOL) L1,(ISYTLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
           WRITE (LUERR,*) 
     &      'ERROR CC_SLVE: L from LUCSOL not as expected'
           WRITE (LUERR,*) 'L from 520 loop:',L
           WRITE (LUERR,*) 'L from LUCSOL   :',L1
           CALL QUIT('ERROR CC_SLVE: L from LUCSOL not as expected')
         END IF
C
        DO 500 M = -L,L
         LM = LM + 1
         IF (IPRINT .GE. 15) THEN
           WRITE(LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                               ' ===================='
           WRITE(LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYTLM(L+M+1))
         END IF
         IF (ISYTLM(L+M+1) .NE. 1) THEN
Chj         IF (ABS(ETLM(LM,2)) .GT. THRZER) THEN
Chj numerical errors have seen to give 1.0D-12,
Chj so we now use 1.0D-10 for test /Feb 2005
            IF (ABS(ETLM(LM,2)) .GT. 1.0D-10) THEN
              WRITE(LUPRI,*) 'ERROR CC_SLVE for l,m',L,M
              WRITE(LUPRI,*) 'Symmetry :',ISYTLM(L+M+1)
              WRITE(LUPRI,*) 'Tn(l,m) .ne. 0, but =',ETLM(LM,2)
              CALL QUIT( 'ERROR CC_SLVE: Tn(l,m) not 0 as expected')
            END IF
            ETLM(LM,2) = ZERO
C           ... to fix round-off errors in Tn(l,m) calculation
            IF (ISYTLM(L+M+1) .GT. 1) READ (LUCSOL)
            GO TO 500
         END IF
C
C--------------------------------
C        Read T(l,m) in AO basis.
C--------------------------------
C
         KTLMAO  = 1
         KWRK1   = KTLMAO  + N2BST(ISYMOP)
         LWRK1   = LWORK   - KWRK1          
         IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CC_SLVE, 1')
C
         CALL CC_SLVINT(WORK(KTLMAO),WORK(KWRK1),LWRK1,ISYMOP)
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('cc_slve: Tlm_ao matrix: - cc storage')
            CALL CC_PRFCKAO(WORK(KTLMAO),ISYMOP)
         ENDIF
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('One electron density in cc_slve')
            CALL CC_PRFCKAO(AODEN,ISYMOP)
         ENDIF
C
C------------------------------------------------------
C        Add electronic contribution TE(l,m) to T(l,m)
C        Contract CC density with integrals.
C        SLV98,OC
C------------------------------------------------------
C
         TELM     = DDOT(N2BST(ISYMOP),WORK(KTLMAO),1,AODEN,1)
C
         IF (IPRINT .GE. 6) THEN
            WRITE(LUPRI,'(A,2I5,/A,3F17.8)')
     *      ' >>> l, m :',L,M,
     *      '     Te(lm), Tn(lm), T(lm) :',
     *         TELM,ETLM(LM,2),ETLM(LM,2)-TELM
         END IF
C
         ETLM(LM,2) = ETLM(LM,2) - TELM
C
C To avoid numerical stability problems.  
C SLV98,OC Necessary???          
C
         IF (ABS(ETLM(LM,2)) .LE. THRZER) THEN
            ETLM(LM,2) = ZERO
            GO TO 500
         END IF
C 
  500   CONTINUE
  520 CONTINUE
C
C---------------------------------
C     Add up energy contributions.
C     Save <Tlm>'s.
C---------------------------------
C
      ESOLT = ZERO
      LM = 0
      DO 920 L = 0,LMAXCU
        DO 900 M = -L,L
          LM = LM + 1
          CCTLM(LM) = ETLM(LM,2)
          ETLM(LM,1) = GL(L,EPSTCU,RCAVCU)*ETLM(LM,2)*ETLM(LM,2)
          IF (IPRINT.GT.5) THEN
             WRITE(LUPRI,*) 'L,M,Energy cont.',L,M,ETLM(LM,1)
             WRITE(LUPRI,*) 'ETLM2,GL',ETLM(LM,2),GL(L,EPSTCU,RCAVCU)
          ENDIF
          ESOLT    = ESOLT  + ETLM(LM,1)
  900   CONTINUE
  920 CONTINUE
C
      CALL CCSL_TLMPUT
C
      CALL GPCLOSE(LUCSOL,'KEEP')
C
      END
C
      DOUBLE PRECISION FUNCTION GL(L,EPS,A)
C
C-----------------------------------------------------------------------------
C
C   Purpose: Calculate G_l(eps) factors for solvent calculations.
C            Cavity radius = A and dielectric constant = EPS
C            based on solfl
C     
C            Ove Christiansen, April 1998
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
C
C     Parameters: FLFAC = - 1 / (8 pi epsilon_null)   (S.I.)
C                       = - 1 / 2                     (a.u.)
C
      PARAMETER ( FLFAC = -0.5D0 , D1 = 1.0D0 )
C
      RL = L
      GL =   FLFAC * A**(-(2*L+1))* (RL + D1) * (EPS - D1)
     *       / (RL + EPS*(RL + D1))
      END
C
      DOUBLE PRECISION FUNCTION GL2(L,EPSST,EPSOP,A)
C
C-----------------------------------------------------------------------------
C
C   Purpose: Calculate G_l(epsst,epsop) factors for solvent calculations. 
C            Ove Christiansen, April 1998
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
C     
      EXTERNAL GL
      DOUBLE PRECISION GL
C
      GL2 = GL(L,EPSST,A) - GL(L,EPSOP,A)
C 
      END 
      SUBROUTINE CC_AOUP(AOUP,AOP,WORK,LWORK,ISYM)
C
C-----------------------------------------------------------------------------
C     Purpose: 
C              Transform from packed  to symmetry packed but unpacked form.
C
C     Ove Christiansen, April 1998.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      DIMENSION AOUP(*),AOP(*),WORK(LWORK)
C
      IF (LWORK.LT.(NBAST**2)) CALL QUIT( 'too little work in CC_AOUP')
C
      CALL DSPTGE(NBAST,AOP(1),WORK(1))
C
      IF (IPRINT.GT.50) THEN
        CALL AROUND('CC_AOUP: Integrals')
        CALL OUTPUT(WORK(1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
      DO ISYMA = 1, NSYM
        ISYMB = MULD2H(ISYM,ISYMA)
        DO A = 1, NBAS(ISYMA)
         DO B = 1, NBAS(ISYMB)
          IOLD = (IBAS(ISYMA)+A-1)*NBAST + (IBAS(ISYMB)+B)
          INEW = IAODIS(ISYMB,ISYMA) + NBAS(ISYMB)*(A-1) + B
          AOUP(INEW) = WORK(IOLD)
         END DO
        END DO
      END DO
      END
C  /* Deck cc_slvint */
      SUBROUTINE CC_SLVINT(TLMAO,WORK,LWORK,ISYMT)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Readin solvent integrals in coupled cluster format.
C
C     SLV98,OC
C     Ove Christiansen, April 1998.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),TLMAO(*)
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*) 'CC_SLVINT: Read in integrals'
        WRITE(LUPRI,*) 'Input symmetry claimed', isymt 
      ENDIF
C
      KTLMAOP = 1
      KWRK1   = KTLMAOP + NNBASX
      LWRK1   = LWORK   - KWRK1
      IF (LWRK1.LT.2*NNBASX) CALL QUIT( 'Too little work in CC_SLVINT')
C
      CALL READT(LUCSOL,NNBASX,WORK(KTLMAOP))
C
      IF (IPRINT .GE. 25) THEN
         CALL AROUND('CC_SLVINT: Tlm_ao matrix: - packed ')
         CALL OUTPAK(WORK(KTLMAOP),NBAST,1,LUPRI)
      END IF
C
      CALL CC_AOUP(TLMAO,WORK(KTLMAOP),WORK(KWRK1),LWRK1,ISYMT)
C
      END
C  /* Deck ccsl_rhstg */
      SUBROUTINE CCSL_RHSTG(FOCK,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              solvent effects.
C
C     SLV98,OC
C     Ove Christiansen, April 1998.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),FOCK(*)
C
      LOGICAL FIRST
      SAVE  FIRST
      DATA  FIRST/.TRUE./
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*) 
     *    'CCSL_RHSTG: Solvent contribution to CC equations.'
      ENDIF
C
C---------------------------------
C     Read in integrals from  file.
C---------------------------------
C
      LUCSOL = -1
      CALL GPOPEN(LUCSOL,'AOSOLINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCSOL)
      CALL MOLLAB('SOLVRLM ',LUCSOL,LUERR)
      IF (FIRST) THEN
         READ (LUCSOL) LMAXSS, LMTOT, NNNBAS
         IF (LMAXSS .LT. LMAXCU) THEN
            WRITE (LUERR,'(//2A,2(/A,I5))') ' >>> CCSL_RHSTG ERROR,',
     *      ' insufficient number of intgrals on LUCSOL',
     *      ' l max from CC     input :',LMAXCU,
     *      ' l max from LUCSOL  file  :',LMAXSS
            CALL QUIT('CCSL_RHSTG: lmax on LUCSOL is too small')
         END IF
         IF ((LMAXSS+1)**2 .NE. LMTOT) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> CCSL_RHSTG ERROR,',
     *      ' LUCSOL file info inconsistent',
     *      ' l_max               :',LMAXSS,
     *      ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *      ' LMTOT               :',LMTOT
          CALL QUIT('CCSL_RHSTG: LUCSOL info not internally consistent')
         END IF
         IF (NNNBAS .NE. NBAST) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> CCSL_RHSTG ERROR,',
     *      ' LUCSOL file info inconsistent with SIRIUS input',
     *      ' NBAST - LUCSOL       :',NNNBAS,
     *      ' NBAST - SIRIUS      :',NBAST
            CALL QUIT('CCSL_RHSTG: LUCSOL info not '//
     &                'consistent with SIRIUS input.')
         END IF
         FIRST = .FALSE.
      ELSE
         READ (LUCSOL)
      END IF
      CALL READT(LUCSOL,NLMCU,WORK(1))
      IF (IPRINT.GT.15) THEN
         WRITE(LUPRI,*) 'Common TLM'
         LM = 0
         DO L = 0,LMAXCU
           DO M = -L,L
             LM = LM + 1
             WRITE(LUPRI,*) 'LM,TLM:',LM,CCTLM(LM)
           ENDDO
         ENDDO 
      ENDIF
C
      CALL CCSL_TLMGET
C
      IF (IPRINT.GT.15) THEN
         WRITE(LUPRI,*) 'TLM from FIL'
         LM = 0
         DO L = 0,LMAXCU
           DO M = -L,L
             LM = LM + 1
             WRITE(LUPRI,*) 'LM,TLM:',LM,CCTLM(LM)
           ENDDO
         ENDDO 
      ENDIF
C
      LM = 0
      DO 520 L = 0,LMAXCU
         READ (LUCSOL) L1,(ISYTLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
           WRITE (LUERR,*) 
     *     'ERROR CCSL_RHSTG: L from LUCSOL not as expected'
           WRITE (LUERR,*) 'L from 520 loop:',L
           WRITE (LUERR,*) 'L from LUCSOL   :',L1
           CALL QUIT('ERROR CCSL_RHSTG: L from LUCSOL not as expected')
         END IF
C
        DO 500 M = -L,L
         LM = LM + 1
         IF (IPRINT .GE. 15) THEN
           WRITE(LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                               ' ===================='
           WRITE(LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYTLM(L+M+1))
         END IF
         IF (ISYTLM(L+M+1) .NE. 1) THEN
            READ (LUCSOL)
            GO TO 500
         ENDIF
C
C--------------------------------
C        Read T(l,m) in AO basis.
C--------------------------------
C
         ISYMTLM = ISYTLM(L+M+1)
         KTLMAO  = 1
         KWRK1   = KTLMAO  + N2BST(ISYMTLM)
         LWRK1   = LWORK   - KWRK1
         IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_RHSTG, 1')
C
         CALL CC_SLVINT(WORK(KTLMAO),WORK(KWRK1),LWRK1,ISYMTLM)
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('ccsl_rhstg: Tlm_ao matrix: - cc storage')
            CALL CC_PRFCKAO(WORK(KTLMAO),ISYMTLM)
         ENDIF
C
C SLV98,OC  Take care with sign of FACT!
C
         FACT =  -2.0D0*GL(L,EPSTCU,RCAVCU)*CCTLM(LM)
         IF (IPRINT.GT.10) THEN
            WRITE(LUPRI,*) 'CCSL_RHSTG: GL ',GL(L,EPSTCU,RCAVCU)
            WRITE(LUPRI,*) 'CCSL_RHSTG: TLM',CCTLM(LM)
            WRITE(LUPRI,*) 'CCSL_RHSTG: FAC',FACT
         ENDIF
C
C-----------------------------------------------------
C        Add contribution to effective AO fock matrix.
C-----------------------------------------------------
C
         CALL DAXPY(N2BST(ISYMTLM),FACT,WORK(KTLMAO),1,FOCK,1)
C
  500   CONTINUE
  520 CONTINUE
C
      CALL GPCLOSE(LUCSOL,'KEEP')
      END

C  /* Deck ccsl_ltrb */
      SUBROUTINE CCSL_LTRB(RHO1,RHO2,CTR1,CTR2,ISYMTR,LR,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Calculation of Coupled Cluster solvent T^gB contribution
C              to left and right Jacobian transformation.
C              <mu|exp(-T)T^g|CC> for LR = 0
C              F transformation for LR = F
C              P transformation for LR = P
C
C     LR = 'L','R','0','F','P'
C     SLV98,OC
C     Ove Christiansen, April/mai 1998.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*),CTR1(*),CTR2(*) 
C
      CHARACTER*8 LABEL,LIST*(2),LR*(1)
      CHARACTER*8 FILMME, FILMMX
      LOGICAL LEXIST
C
      DOUBLE PRECISION GL
      EXTERNAL GL
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCSL_LTRB: Solvent contribution to CC L. transf.'
        WRITE(LUPRI,*)'CCSL_LTRB: LWORK:', LWORK
        WRITE(LUPRI,*)'CCSL_LTRB: LR:', LR    
        WRITE(LUPRI,*)'CCSL_LTRB: ISYMTR:', ISYMTR
      ENDIF
      IF ( DEBUG .OR.(IPRINT.GT.10)) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
         RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
         WRITE(LUPRI,*) ' Norm af RHO1 in CCSL_LTGB on input:', RHO1N
         WRITE(LUPRI,*) ' Norm af RHO2 in CCSL_LTGB on input:', RHO2N
         IF (LR .NE. '0') THEN 
           RHO1N = DDOT(NT1AM(ISYMTR),CTR1,1,CTR1,1)
           RHO2N = DDOT(NT2AM(ISYMTR),CTR2,1,CTR2,1)
           WRITE(LUPRI,*) ' Norm af C1AM in CCSL_LTGB on input:', RHO1N
           WRITE(LUPRI,*) ' Norm af C2AM in CCSL_LTGB on input:', RHO2N
         END IF
      ENDIF
C
C     Note if CCSAV_LAM does not exist then
C     we have no contribution yet. 
C
c     INQUIRE(FILE='CCSAV_LAM',EXIST=LEXIST)
c     IF (.NOT.LEXIST) THEN
c        WRITE(LUPRI,*) ' CCSAV_LAM does not exits yet - no T^gB cont'
c        RETURN
c     ENDIF
C
C---------------------
C     Init parameters.
C---------------------
C
      NTAMP1  = NT1AM(ISYMTR)
      NTAMP2  = NT2AM(ISYMTR)
      IF (CCS)  NTAMP2  = 0
      NTAMP   = NTAMP1 + NTAMP2
C
      IF (DISCEX) THEN
        LUMMET = -1
        LUMMXI = -1
        FILMME = 'CCSL_ETA'
        FILMMX = 'CCSL__XI'
        CALL WOPEN2(LUMMET, FILMME, 64, 0)
        CALL WOPEN2(LUMMXI, FILMMX, 64, 0)
      END IF
C
C---------------------------------
C     Read in integrals from  file.
C---------------------------------
C
      LUCSOL = -1
      CALL GPOPEN(LUCSOL,'AOSOLINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCSOL)
      CALL MOLLAB('SOLVRLM ',LUCSOL,LUERR)
      READ (LUCSOL)
      CALL READT(LUCSOL,NLMCU,WORK(1))
C
      KTGB    = 1
      KWRK1   = KTGB    + N2BST(ISYMTR)
      LWRK1   = LWORK   - KWRK1
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_LTRB, 1')
      KWRK2   = KWRK1
      LWRK2   = LWRK1
C
      CALL DZERO(WORK(KTGB),N2BST(ISYMTR))
C
      LM = 0
      IADR1 = 1
      IADR2 = 1
C
      DO 600 L = 0,LMAXCU
         READ (LUCSOL) L1,(ISYTLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
           WRITE (LUERR,*) 
     *     'ERROR CCSL_LTRB: L from LUCSOL not as expected'
           WRITE (LUERR,*) 'L from 600 loop:',L
           WRITE (LUERR,*) 'L from LUCSOL   :',L1
           CALL QUIT('ERROR CCSL_LTRB: L from LUCSOL not as expected')
         END IF
C
        DO 500 M = -L,L
         LM = LM + 1
         IF (IPRINT .GE. 15) THEN
           WRITE(LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                               ' ===================='
           WRITE(LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYTLM(L+M+1))
         END IF
C
         ISYMNY = ABS(ISYTLM(L+M+1))
         NTA1  = NT1AM(ISYMNY)
         NTA2  = NT2AM(ISYMNY)
         NTA   = NTA1 + NTA2
C
C----------------------------------------------------------
C        Symmetry should be equal to trial vector symmetry.
C----------------------------------------------------------
C
         IF (ISYTLM(L+M+1) .NE. ISYMTR) THEN
            READ (LUCSOL)
            IADR1 = IADR1 + NTA
            IADR2 = IADR2 + NTA
            GO TO 500
         ENDIF
C
C-------------------------------------------------------------------
C        Calculate T^{gB}= -sum_lm 2G_l(ep)<B|T_lm|CC>T_lm operator.
C        1. Readin integrals again.
C        2. Calculate xi^T_lm  
C        3. Contract with B
C        4. scale integrals with 2G_l(ep)<B|T_lm|CC>
C        5. Add to T^{gB} integrals.
C-------------------------------------------------------------------
C
C
C--------------------------------
C        1. Readin integrals again.
C        Read T(l,m) in AO basis.
C--------------------------------
C
         KTLMAO  = KWRK1
         KWRK2   = KTLMAO  + N2BST(ISYMTR)
         LWRK2   = LWORK   - KWRK2
         IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_LTRB, 2')
C
         CALL CC_SLVINT(WORK(KTLMAO),WORK(KWRK2),LWRK2,ISYMTR)
C
         IF (IPRINT .GT. 25) THEN
            CALL AROUND('ccsl_ltrb: Tlm_ao matrix: - cc storage')
            CALL CC_PRFCKAO(WORK(KTLMAO),ISYMTR)
         ENDIF
C
C-----------------------------------------------------------
C        2. Calculate xi^T_lm  (for LR=R actually eta^T_lm )
C-----------------------------------------------------------
C
         KXI     = KWRK2
         KWRK3   = KXI     + NTAMP
         LWRK3   = LWORK   - KWRK3
         IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CCSL_LTRB, 3')
C
         IF (.NOT. (DISCEX)) THEN
           LABEL = 'GIVE INT'           
           IF ((LR.EQ.'L').OR.(LR.EQ.'0').OR.(LR.EQ.'P')) THEN
             CALL CC_XKSI(WORK(KXI),LABEL,ISYMTR,0,WORK(KTLMAO),
     *                    WORK(KWRK3),LWRK3)
           ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
             LIST  = 'L0'
             ILSTNR  = 1
             CALL CC_ETAC(ISYMTR,LABEL,WORK(KXI),
     *                    LIST,ILSTNR,0,WORK(KTLMAO),WORK(KWRK3),LWRK3)
           ENDIF
         ELSE
           IF ((LR.EQ.'L').OR.(LR.EQ.'P')) THEN
             CALL GETWA2(LUMMXI,FILMMX,WORK(KXI),IADR2,NTAMP)
             IADR2 = IADR2 + NTAMP
           ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
             LIST  = 'L0'
             ILSTNR  = 1
             CALL GETWA2(LUMMET,FILMME,WORK(KXI),IADR1,NTAMP)
             IADR1 = IADR1 + NTAMP
           ELSE IF (LR.EQ.'0') THEN
             LABEL = 'GIVE INT'
             CALL CC_XKSI(WORK(KXI),LABEL,ISYMTR,0,WORK(KTLMAO),
     *                    WORK(KWRK3),LWRK3)
           END IF
         END IF
C
C---------------------------
C        3. Contract with B
C---------------------------
C
         IF (LR.NE.'0') THEN
           KXI1 = KXI  
           KXI2 = KXI  + NTAMP1
C
           BXILMD1= DDOT(NTAMP1,CTR1,1,WORK(KXI1),1)
           BXILMD2= DDOT(NTAMP2,CTR2,1,WORK(KXI2),1)
           BXILMD = BXILMD1 + BXILMD2
C
C----------------------------------------------------
C           4. Find 2G_l(ep)<B|T_lm|CC> factor
C----------------------------------------------------
C
            GLFAC =  GL(L,EPOPCU,RCAVCU)
C
C SLV98,OC  
C
            FACT =   2.0D0*GLFAC*BXILMD   
C
            IF (IPRINT.GT.10) THEN
               WRITE(LUPRI,*) 'CCSL_LTRB: GL    ',GLFAC
               WRITE(LUPRI,*) 'CCSL_LTRB: BTLM1 ',BXILMD1  
               WRITE(LUPRI,*) 'CCSL_LTRB: BTLM2 ',BXILMD2  
               WRITE(LUPRI,*) 'CCSL_LTRB: BTLM  ',BXILMD   
               WRITE(LUPRI,*) 'CCSL_LTRB: FAC   ',FACT
            ENDIF         
C
C--------------------------------------------------------------
C           5. Add to T^{gB} integrals.(for LR=R actually T^gC)
C--------------------------------------------------------------
C
            CALL DAXPY(N2BST(ISYMTR),FACT,WORK(KTLMAO),1,
     *                 WORK(KTGB),1)

         ELSE IF (LR.EQ.'0') THEN
            KXI1 = KXI  
            KXI2 = KXI  + NTAMP1
            FACT =   -2.0D0*GL(L,EPSTCU,RCAVCU)*CCTLM(LM)
            CALL DAXPY(NT1AM(ISYMTR),FACT,WORK(KXI1),1,
     *                 RHO1,1) 
            CALL DAXPY(NT2AM(ISYMTR),FACT,WORK(KXI2),1,
     *                 RHO2,1) 
         ENDIF
C
  500   CONTINUE
  600 CONTINUE
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('ccsl_ltrb: T^gB_ao matrix: ')
         CALL CC_PRFCKAO(WORK(KTGB),ISYMTR)
      ENDIF
C
      IF (DISCEX) THEN
        CALL WCLOSE2(LUMMET, FILMME, 'KEEP')
        CALL WCLOSE2(LUMMXI, FILMMX, 'KEEP')
      END IF
C
      IF (LR.NE.'0') THEN
C
C-----------------------------------------------------
C       Calculate contribution from the T^gB operator.
C-----------------------------------------------------
C
        KETA    = KWRK2
        KWRK4   = KETA    + NTAMP
        LWRK4   = LWORK   - KWRK4
        IF (LWRK4.LT.0) CALL QUIT( 'Too little work in CCSL_LTRB 4')
        KETA1   = KETA
        KETA2   = KETA + NTAMP1
C
        IF ((LR.EQ.'L').OR.(LR.EQ.'F')) THEN
          LIST  = 'L0'
          LABEL = 'GIVE INT'           
          CALL CC_ETAC(ISYMTR,LABEL,WORK(KETA),
     *                 LIST,1,0,WORK(KTGB),WORK(KWRK4),LWRK4)
        ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'P')) THEN
          LABEL = 'GIVE INT'           
          CALL CC_XKSI(WORK(KETA),LABEL,ISYMTR,0,WORK(KTGB), 
     *                 WORK(KWRK4),LWRK4)
          IF (LR.EQ.'R')
     *       CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYMTR)
        ENDIF
C
        IF ( DEBUG .OR.(IPRINT.GT.10)) THEN
           RHO1N = DDOT(NT1AM(ISYMTR),WORK(KETA1),1,WORK(KETA1),1)
           RHO2N = DDOT(NT2AM(ISYMTR),WORK(KETA2),1,WORK(KETA2),1)
           WRITE(LUPRI,*) ' Norm af T^gB contribution to LHTR1:', RHO1N
           WRITE(LUPRI,*) ' Norm af T^gB contribution to LHTR2:', RHO2N
        ENDIF
C
        CALL DAXPY(NT1AM(ISYMTR),ONE,WORK(KETA1),1,RHO1,1)
        CALL DAXPY(NT2AM(ISYMTR),ONE,WORK(KETA2),1,RHO2,1)
        IF ( DEBUG .OR.(IPRINT.GT.10)) THEN
           RHO1N = DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
           RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
           WRITE(LUPRI,*) ' Norm af RHO1 in CCSL_LTGB:', RHO1N
           WRITE(LUPRI,*) ' Norm af RHO2 in CCSL_LTGB:', RHO2N
        ENDIF
C
      ENDIF
C
      CALL GPCLOSE(LUCSOL,'KEEP')
      END
C  /* Deck ccsl_tlmget*/
      SUBROUTINE CCSL_TLMGET
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      LUCTLM = -1
      CALL GPOPEN(LUCTLM,'CC_TL','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      LM = 0
      DO 820 L = 0,LMAXCU
        DO 800 M = -L,L
          LM = LM + 1
          READ(LUCTLM,'(I5,E25.15)',END=100,ERR=100) LM1,CCTLM(LM)
          IF (LM1.NE.LM) CALL QUIT( 'CCSL_TLMGET: Error on CC_TL')
          GOTO 800
  100     CCTLM(LM) = 0.0D0   
C
  800   CONTINUE
  820 CONTINUE
C
      CALL GPCLOSE(LUCTLM,'KEEP')
      RETURN
C
      END
C
C  /* Deck ccsl_tlmput*/
      SUBROUTINE CCSL_TLMPUT
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      LUCTLM = -1
      CALL GPOPEN(LUCTLM,'CC_TL','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
      LM = 0
      DO 920 L = 0,LMAXCU
        DO 900 M = -L,L
          LM = LM + 1
          WRITE(LUCTLM,'(I5,E25.15)') LM,CCTLM(LM) 
  900   CONTINUE
  920 CONTINUE
C
      CALL GPCLOSE(LUCTLM,'KEEP')
C
      END
****************************************************************
C  /* Deck ccsl_pbtr */
      SUBROUTINE CCSL_PBTR(RHO1,RHO2,ISYRES,
     *                     LISTL,IDLSTL,CTR1,ISYCTR,
     *                     LISTR,IDLSTR,BTR1,ISYBTR,
     *                     MODEL,WORK,LWORK) 
C
C---------------------------------------------------------------
C Calculates contributions from the T^gB , the ^CT^g 
C and ^CT^gB operators.
C JK+OC, 03
C---------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*),CTR1(*),BTR1(*)
C
      CHARACTER*(*) LISTL,LISTR,LIST*(2)
      CHARACTER*8 LABEL
      CHARACTER*10 MODEL
      LOGICAL LEXIST
      INTEGER IDLSTL, IDLSTR, ISYCTR, ISYBTR, ISYBC
C
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCSL_PBTR: ISYRES =',ISYRES
        WRITE(LUPRI,*)'CCSL_PBTR: ISYCTR =',ISYCTR
        WRITE(LUPRI,*)'CCSL_PBTR: ISYBTR =',ISYBTR
        WRITE(LUPRI,*)'CCSL_PBTR: LISTL =',LISTL
        WRITE(LUPRI,*)'CCSL_PBTR: LISTR =',LISTR
        WRITE(LUPRI,*)'CCSL_PBTR: IDLSTL =',IDLSTL
        WRITE(LUPRI,*)'CCSL_PBTR: IDLSTR =',IDLSTR
        CALL FLSHFO(LUPRI)
      ENDIF
C
C---------------------
C     Init parameters.
C---------------------
C
C     For the B (right) trial vector
C
      NAMP1B  = NT1AM(ISYBTR)
      NAMP2B  = NT2AM(ISYBTR)
      NAMPB   = NAMP1B + NAMP2B
C
C     For the C (left) trial vector
C
      NAMP1C  = NT1AM(ISYCTR)
      NAMP2C  = NT2AM(ISYCTR)
      NAMPC   = NAMP1C + NAMP2C
C
C     For the F = C X B vector 
C
      ISYBC = MULD2H(ISYBTR,ISYCTR)
      NAMP1F  = NT1AM(ISYBC)
      NAMP2F  = NT2AM(ISYBC)
      NAMPF   = NAMP1F + NAMP2F
C
C-----------------------------------------------------------
C       Calculate contribution from the effective operators.
C-----------------------------------------------------------
C
C     First we concentrate on the effective operator T^gB
C     -----------------------------------------------------
C
      KTGB    = 1
      KWRK1   = KTGB + N2BST(ISYBTR)
      LWRK1   = LWORK   - KWRK1
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_PBTR, 1')
C
      CALL DZERO(WORK(KTGB),N2BST(ISYBTR))
C
      CALL CCSL_TGB(BTR1,ISYBTR,LISTR,IDLSTR,WORK(KTGB),'ET',
     *              MODEL,WORK(KWRK1),LWRK1)
C
C     Symmetry: 
C     ---------
C
      ISYBC = MULD2H(ISYBTR,ISYCTR)
C
      IF (ISYBC .NE. ISYRES) THEN
        CALL QUIT( 'Symmetry problem in CCSL_PBTR')
      END IF
C
      KETA    = KWRK1
      KWRK2   = KETA + NAMPF 
      LWRK2   = LWORK   - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_PBTR, 2')
C
C     Note, LISTL .EQ. LE/L1 so the HF part of the following 
C     eta-transformation is skipped
C
      LABEL = 'GIVE INT'
      CALL CC_ETAC(ISYBTR,LABEL,WORK(KETA),
     *               LISTL,IDLSTL,0,WORK(KTGB),WORK(KWRK2),LWRK2)
C
      KETA1   = KETA
      KETA2   = KETA1 + NAMP1F 
C
      CALL DAXPY(NAMP1F,ONE,WORK(KETA1),1,RHO1,1)
      CALL DAXPY(NAMP2F,ONE,WORK(KETA2),1,RHO2,1)
C
C     Next, the effective operator ^CT^gB is considered
C----------------------------------------------------------
C
C     Symmetry of the operator:
C     -------------------------
C
      ISYBC = MULD2H(ISYBTR,ISYCTR)
C
      KTGB    = 1
      KWRK1   = KTGB + N2BST(ISYBC)
      LWRK1   = LWORK   - KWRK1
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_PBTR, 3')
C
      CALL DZERO(WORK(KTGB),N2BST(ISYBC))
C
      CALL CCSL_CTGB(LISTL,IDLSTL,CTR1,ISYCTR,
     *               LISTR,IDLSTR,BTR1,ISYBTR,
     *               WORK(KTGB),MODEL,WORK(KWRK1),LWRK1)
C
      KETA    = KWRK1
      KWRK2   = KETA + NAMPF
      LWRK2   = LWORK   - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_PBTR, 4')
C
      CALL DZERO(WORK(KETA),NAMPF)
C
      LABEL = 'GIVE INT'
      LIST  = 'L0'
      IDLINO = 1
C
      CALL CC_ETAC(ISYBC,LABEL,WORK(KETA),
     *             LIST,IDLINO,0,WORK(KTGB),WORK(KWRK2),LWRK2)
C
      KETA1   = KETA
      KETA2   = KETA1 + NAMP1F
C
      CALL DAXPY(NT1AM(ISYBC),ONE,WORK(KETA1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYBC),ONE,WORK(KETA2),1,RHO2,1)
C
C     Finally, we calculate the third contribution which is
C     a contribution from a T^gB operator but calculated
C     as a xi vector element.
C-----------------------------------------------------------
C
      KTGB    = 1
      KWRK1   = KTGB + N2BST(ISYCTR)
      LWRK1   = LWORK   - KWRK1
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_PBTR, 5')
C
      CALL DZERO(WORK(KTGB),N2BST(ISYCTR))
C
      CALL CCSL_TGB(CTR1,ISYCTR,LISTL,IDLSTL,WORK(KTGB),'XI',
     *              MODEL,WORK(KWRK1),LWRK1)
C
C     Symmetry:
C     ---------
C
      ISYBC = MULD2H(ISYBTR,ISYCTR)
C
      IF (ISYBC .NE. ISYRES) THEN
        CALL QUIT( 'Symmetry problem in CCSL_PBTR')
      END IF

      LABEL = 'GIVE INT'
      LIST  = 'L0'
      LISTNO = 1
C
C     (NB: Result vector from the following transformation is 
C     given at the beginning of WORK(KWRK1).)
C
      CALL CCLR_FA(LABEL,ISYCTR,LISTR,IDLSTR,
     &             LIST,LISTNO,WORK(KTGB),WORK(KWRK1),LWRK1)
C
C      Jacob, Husk at DDOTog DAXPY virker paa denne her maade !!  
C      TEMP1 = DDOT(NT1AM(ISYBC),WORK(KWRK1),1,WORK(KWRK1),1)
C      TEMP2 = DDOT(NT2AM(ISYBC),WORK(KWRK1+NT1AM(ISYBC)),1,
C     *                          WORK(KWRK1+NT1AM(ISYBC)),1)
C
      CALL DAXPY(NT1AM(ISYBC),ONE,WORK(KWRK1),1,RHO1,1) 
      CALL DAXPY(NT2AM(ISYBC),ONE,
     *           WORK(KWRK1+NT1AM(ISYBC)),1,RHO2,1)
C
      END
****************************************************************
C  /* Deck ccsl_ctgb */
      SUBROUTINE CCSL_CTGB(LISTL,IDLSTL,CTR1,ISYCTR,
     *                     LISTR,IDLSTR,BTR1,ISYBTR,
     *                     TGB,MODEL,WORK,LWORK) 
C
C JK+OC, 03
C-----------------------------------------------------------------------------
C
C   Construcs effective operator equal to
C   ^CT^gB = -2 Sum_lm g_l * Sum_sigma Sum_my
C            t^C_my <bar my|[T_lm,Tau_sigma]|CC> t^B_sigma T_lm
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),BTR1(*),CTR1(*)
      DIMENSION TGB(*) 
C
      CHARACTER*8 LABEL, LIST*(2)
      CHARACTER*(*) LISTL, LISTR
      CHARACTER*10 MODEL
      INTEGER IDLSTL, IDLSTR, ISYCTR, ISYBTR
      LOGICAL LEXIST
C
      DOUBLE PRECISION GL
      EXTERNAL GL
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCSL_PBTR: ISYCTR =',ISYCTR
        WRITE(LUPRI,*)'CCSL_PBTR: ISYBTR =',ISYBTR
        WRITE(LUPRI,*)'CCSL_PBTR: LISTL =',LISTL
        WRITE(LUPRI,*)'CCSL_PBTR: LISTR =',LISTR
        WRITE(LUPRI,*)'CCSL_PBTR: IDLSTL =',IDLSTL
        WRITE(LUPRI,*)'CCSL_PBTR: IDLSTR =',IDLSTR
        CALL FLSHFO(LUPRI)
      ENDIF
C
C---------------------
C     Init parameters.
C---------------------
C
C     For the B (right) trial vector
C 
      NAMP1B  = NT1AM(ISYBTR)
      NAMP2B  = NT2AM(ISYBTR)
      NAMPB   = NAMP1B + NAMP2B
C
C     For the C (left) trial vector
C
      NAMP1C  = NT1AM(ISYCTR)
      NAMP2C  = NT2AM(ISYCTR)
      NAMPC   = NAMP1C + NAMP2C
C
C     For the F = B x C vector
C
      ISYBC = MULD2H(ISYBTR,ISYCTR)
      NAMP1F  = NT1AM(ISYBC)
      NAMP2F  = NT2AM(ISYBC)
      NAMPF   = NAMP1F + NAMP2F
C
C-----------------------------------------------
C    Readin response vector (B, right) from file
C-----------------------------------------------
C
      KT2AMPA = 1
      KWRK1   = KT2AMPA +  NT2AM(ISYBTR)
      LWRK1   = LWORK - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insuff. work in CCSL_TGB')
      END IF
C
      IOPT = 2
      CALL CC_RDRSP(LISTR,IDLSTR,ISYBTR,IOPT,MODEL,
     *              WORK(KWRK1),WORK(KT2AMPA))
C
C     First we perform the contraction of of the eta'
C     element from the left hand side. The ' means that
C     the HF contribution is skipped.
C 
C---------------------------------
C     Read in integrals from  file.
C---------------------------------
C
      LUCSOL = -1
      CALL GPOPEN(LUCSOL,'AOSOLINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCSOL)
      CALL MOLLAB('SOLVRLM ',LUCSOL,LUERR)
      READ (LUCSOL)
      CALL READT(LUCSOL,NLMCU,WORK(KWRK1))
C
      LM = 0
      DO 600 L = 0,LMAXCU
         READ (LUCSOL) L1,(ISYTLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
           WRITE (LUERR,*)
     *     'ERROR CCSL_CTGB: L from LUCSOL not as expected'
           WRITE (LUERR,*) 'L from 600 loop:',L
           WRITE (LUERR,*) 'L from LUCSOL   :',L1
           CALL QUIT('ERROR CCSL_CTGB: L from LUCSOL not as expected')
         END IF
C
        DO 500 M = -L,L
         LM = LM + 1
         IF (IPRINT .GE. 15) THEN
           WRITE(LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                               ' ===================='
           WRITE(LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYTLM(L+M+1))
         END IF
C
C----------------------------------------------
C        Symmetry should be equal to symmetry
C        of the final contracted vector (ISYBC).
C----------------------------------------------
C
         ISYBC = MULD2H(ISYBTR,ISYCTR)
C
         IF (ISYTLM(L+M+1) .NE. ISYBC) THEN
           READ (LUCSOL)
           GO TO 500
         ENDIF
C
C--------------------------------
C        Read T(l,m) in AO basis.
C--------------------------------
C
         KTLMAO = KWRK1
         KWRK2  = KTLMAO + N2BST(ISYBC)
         LWRK2  = LWORK - KWRK2
         IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_CTGB, 2')
C
         CALL CC_SLVINT(WORK(KTLMAO),WORK(KWRK2),LWRK2,ISYBC)
C
         IF (IPRINT .GT. 25) THEN
            CALL AROUND('CCSL_CTGB : Tlm_ao matrix: - cc storage')
            CALL CC_PRFCKAO(WORK(KTLMAO),ISYBC)
         ENDIF
C
C---------------------------------------------------------------------
C        Calculate eta^T_lm
C        (LIST = 'LE' ensures that the HF part of eta^T_lm is skipped.
C---------------------------------------------------------------------
C
         KXI     = KWRK2
         KWRK3   = KXI     + NAMPB
         LWRK3   = LWORK   - KWRK3
         IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CCSL_CTGB, 3')
C
         LABEL = 'GIVE INT'
         CALL CC_ETAC(ISYBC,LABEL,WORK(KXI),
     *                LISTL,IDLSTL,0,WORK(KTLMAO),WORK(KWRK3),LWRK3)
C
C---------------------------------------------------------------
C        Contract with trial vector B (from the right hand side)
C---------------------------------------------------------------
C
         KXI1 = KXI
         KXI2 = KXI  + NT1AM(ISYBTR)
C
         BXILMD1= DDOT(NT1AM(ISYBTR),BTR1,1,WORK(KXI1),1)
         BXILMD2= DDOT(NT2AM(ISYBTR),WORK(KT2AMPA),1,WORK(KXI2),1)
         BXILMD = BXILMD1 + BXILMD2
C
C--------------------------------------------------------------
C        Find  factor from contracted vector and G_lm factor
C--------------------------------------------------------------
C
         GLFAC =  GL(L,EPOPCU,RCAVCU)
         FACT =   2.0D0*GLFAC*BXILMD
C
C--------------------------------------------------------------
C        Put factor on integrals and (for a given l,m) the
C        effective operator has now been constructed.
C--------------------------------------------------------------
C
           CALL DAXPY(N2BST(ISYBC),FACT,WORK(KTLMAO),1,
     *              TGB,1)
C---------------------------------------------------------------

  500   CONTINUE
  600 CONTINUE
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('ccsl_ctgb : T^gB_ao matrix: ')
         CALL CC_PRFCKAO(TGB,ISYBC)
      ENDIF
C
      CALL GPCLOSE(LUCSOL,'KEEP')
      END
*******************************************************************
C  /* Deck ccsl_tgb */
      SUBROUTINE CCSL_TGB(CTR1,ISYMTR,LISTIN,IDLIST,TGB,LR,
     *                    MODEL,WORK,LWORK)
C
C     JK+OC, 03
C-----------------------------------------------------------------------------
C
C   IF (LR.EQ.'ET') then the following effective operator is constructed
C   T^gB = -2 Sum_lm g_l * Sum_sigma
C             <Lambda|[T_lm,Tau_sigma]|CC> t^B_sigma T_lm
C
C   IF (LR.EQ.'XI') then the following effective operator is constructed
C   ^BT^g = -2 Sum_lm g_l * Sum_sigma
C             t^B_sigma <bar sigma|T_lm|CC> T_lm
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),CTR1(*)
      DIMENSION TGB(*)
C
      INTEGER ISYMTR, IDLIST
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      CHARACTER*(*) LISTIN
      CHARACTER*10 MODEL
C
      CHARACTER*8 LABEL,LR*(2),LIST
      CHARACTER*(8) FILMME, FILMMX
      LOGICAL LEXIST
C
      DOUBLE PRECISION GL
      EXTERNAL GL
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCSL_TGB : ISYMTR:', ISYMTR
        WRITE(LUPRI,*)'CCSL_TGB : LISTIN:', LISTIN
      ENDIF
C
C---------------------
C     Init parameters.
C---------------------
C
      NTAMP1  = NT1AM(ISYMTR)
      NTAMP2  = NT2AM(ISYMTR)
      IF (CCS)  NTAMP2  = 0
      NTAMP   = NTAMP1 + NTAMP2
C
C
      IF (DISCEX) THEN
        LUMMET = -1
        LUMMXI = -1
        FILMME = 'CCSL_ETA'
        FILMMX = 'CCSL__XI'
        CALL WOPEN2(LUMMET, FILMME, 64, 0)
        CALL WOPEN2(LUMMXI, FILMMX, 64, 0)
      END IF
C
C--------------------------------------
C     Readin trial vector from file
C--------------------------------------
C
      KT2AMPA = 1
      KWRK1   = KT2AMPA +  NT2AM(ISYMTR)
      LWRK1   = LWORK - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insuff. work in CCSL_TGB')
      END IF
C
      IOPT = 2
      CALL CC_RDRSP(LISTIN,IDLIST,ISYMTR,IOPT,MODEL,
     *              WORK(KDUM),WORK(KT2AMPA))
C
C---------------------------------
C     Read in integrals from  file.
C---------------------------------
C
      LUCSOL = -1
      CALL GPOPEN(LUCSOL,'AOSOLINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCSOL)
      CALL MOLLAB('SOLVRLM ',LUCSOL,LUERR)
      READ (LUCSOL)
      CALL READT(LUCSOL,NLMCU,WORK(KWRK1))
C
      TEMP = 0.00
      LM = 0
      IADR1 = 1
      IADR2 = 1
C
      DO 600 L = 0,LMAXCU
         READ (LUCSOL) L1,(ISYTLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
           WRITE (LUERR,*)
     *     'ERROR CCSL_TGB: L from LUCSOL not as expected'
           WRITE (LUERR,*) 'L from 600 loop:',L
           WRITE (LUERR,*) 'L from LUCSOL   :',L1
           CALL QUIT('ERROR CCSL_TGB: L from LUCSOL not as expected')
         END IF
C
        DO 500 M = -L,L
         LM = LM + 1
         IF (IPRINT .GE. 15) THEN
           WRITE(LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                               ' ===================='
           WRITE(LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYTLM(L+M+1))
         END IF
C
         ISYMNY = ABS(ISYTLM(L+M+1))
         NTA1  = NT1AM(ISYMNY)
         NTA2  = NT2AM(ISYMNY)
         NTA   = NTA1 + NTA2
C
C----------------------------------------------------------
C        Symmetry should be equal to trial vector symmetry.
C----------------------------------------------------------
C
         IF (ISYTLM(L+M+1) .NE. ISYMTR) THEN
            READ (LUCSOL)
            IADR1 = IADR1 + NTA
            IADR2 = IADR2 + NTA
            GO TO 500
         ENDIF
C
C--------------------------------
C        Read T(l,m) in AO basis.
C--------------------------------
C
         KTLMAO = KWRK1
         KWRK2  = KTLMAO + N2BST(ISYMTR)
         LWRK2  = LWORK - KWRK2
         IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_CTGB, 2')
C
         CALL CC_SLVINT(WORK(KTLMAO),WORK(KWRK2),LWRK2,ISYMTR)
C
         IF (IPRINT .GT. 25) THEN
            CALL AROUND('CCSL_TGB : Tlm_ao matrix: - cc storage')
            CALL CC_PRFCKAO(WORK(KTLMAO),ISYMTR)
         ENDIF
C
C---------------------------------------------------------------------
C        (LR.EQ.'ET') Calculate eta^T_lm
C        (LR.EQ.'XI') Calculate xi^T_lm
C---------------------------------------------------------------------
C
         KXI     = KWRK2
         KWRK3   = KXI     + NTAMP
         LWRK3   = LWORK   - KWRK3
         IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CCSL_TGB, 3')
C
         IF (.NOT. (DISCEX)) THEN
           IF (LR.EQ.'ET') THEN
             LABEL = 'GIVE INT'
             LIST  = 'L0'
             IDLINO = 1 
             CALL CC_ETAC(ISYMTR,LABEL,WORK(KXI),
     *                    LIST,IDLINO,0,WORK(KTLMAO),WORK(KWRK3),LWRK3)
C
           ELSE IF (LR.EQ.'XI') THEN
             LABEL = 'GIVE INT'
             CALL CC_XKSI(WORK(KXI),LABEL,ISYMTR,0,WORK(KTLMAO),
     *                   WORK(KWRK3),LWRK3)
           END IF
         ELSE
           IF (LR.EQ.'ET') THEN
             CALL GETWA2(LUMMET,FILMME,WORK(KXI),IADR1,NTAMP)
             IADR1 = IADR1 + NTAMP
           ELSE IF (LR.EQ.'XI') THEN
             CALL GETWA2(LUMMXI,FILMMX,WORK(KXI),IADR2,NTAMP)
             IADR2 = IADR2 + NTAMP
           END IF
         END IF
C
C--------------------------------------------------------------
C        Contract with trial vector: 
C
C        (If LR.EQ.'ET' then from the right hand side)
C        (If LR.EQ.'XI' then from the left hand side)
C
C--------------------------------------------------------------
C
         KXI1 = KXI
         KXI2 = KXI  + NTAMP1
C
C
         CXILMD1= DDOT(NTAMP1,CTR1,1,WORK(KXI1),1)
         CXILMD2= DDOT(NTAMP2,WORK(KT2AMPA),1,WORK(KXI2),1)
         CXILMD = CXILMD1 + CXILMD2
C
C--------------------------------------------------------------
C        Find  factor from contracted vector and G_lm factor
C--------------------------------------------------------------
C
         GLFAC =  GL(L,EPOPCU,RCAVCU)
         FACT =   2.0D0*GLFAC*CXILMD
C
C--------------------------------------------------------------
C        Put factor on integrals and (for a given l,m) the
C        effective operator has now been constructed.
C--------------------------------------------------------------
C
         CALL DAXPY(N2BST(ISYMTR),FACT,WORK(KTLMAO),1,
     *              TGB,1)
C---------------------------------------------------------------

  500   CONTINUE
  600 CONTINUE
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('ccsl_tgb : T^gB_ao matrix: ')
         CALL CC_PRFCKAO(TGB,ISYMTR)
      ENDIF
C
      CALL GPCLOSE(LUCSOL,'KEEP')
C
      IF (DISCEX) THEN
        CALL WCLOSE2(LUMMET, FILMME, 'KEEP')
        CALL WCLOSE2(LUMMXI, FILMMX, 'KEEP')
      END IF
C
      END
*******************************************************************
****************************************************************
C  /* Deck ccsl_gtr */
      SUBROUTINE CCSL_GTR(RHO1,RHO2,ISYRES,
     *                    LISTA,IDLSTA,ISYMTA,
     *                    LISTB,IDLSTB,ISYMTB,
     *                    LISTC,IDLSTC,ISYMTC,
     *                    MODEL,WORK,LWORK)  
C
C-----------------------------------------------------------------
C JK+OC, feb.03
C       Calculates the contribution to the CC/dielectric G matrx
C       transformations.
C
C       G_my,ny,sigma = 1/2 P(my,ny,sigma) 
C                       <Lambda|[[T^g,mu,tau_ny],tau_sigma]|CC >
C
C       The permutation creates 3 different contributions. 
C
C------------------------------------------------------------------
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*)
C
      CHARACTER*(3) LISTA, LISTB, LISTC
      CHARACTER*8 LABEL
      CHARACTER*10 MODEL
      LOGICAL LEXIST, LSAME
      INTEGER ISYCB,ISYRES,ISYMTA,ISYMTB,ISYMTC
      INTEGER IDLSTA,IDLSTB,IDLSTC
C
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCSL_GTR: MODEL =',MODEL(1:4)
        WRITE(LUPRI,*)'CCSL_GTR: ISYRES =',ISYRES
        WRITE(LUPRI,*)'CCSL_GTR: ISYMTA =',ISYMTA
        WRITE(LUPRI,*)'CCSL_GTR: ISYMTB =',ISYMTB
        WRITE(LUPRI,*)'CCSL_GTR: ISYMTC =',ISYMTC
        WRITE(LUPRI,*)'CCSL_GTR: LISTA =',LISTA
        WRITE(LUPRI,*)'CCSL_GTR: LISTB =',LISTB 
        WRITE(LUPRI,*)'CCSL_GTR: LISTC =',LISTC
        WRITE(LUPRI,*)'CCSL_GTR: IDLSTA =',IDLSTA
        WRITE(LUPRI,*)'CCSL_GTR: IDLSTB =',IDLSTB
        WRITE(LUPRI,*)'CCSL_GTR: IDLSTC =',IDLSTC
        CALL FLSHFO(LUPRI)
      ENDIF
C
C     Symmetry:
C     ---------
C
      ISYCB = MULD2H(ISYMTB,ISYMTC)
C
      IF (ISYCB .NE. ISYRES) THEN
        CALL QUIT( 'Symmetry problem in CCSL_GTR 1')
      END IF
C
C
C     First we concentrate on the effective operator T^gB
C     -----------------------------------------------------
C
      KT1AMPA = 1
      KTGB    = KT1AMPA + NT1AM(ISYMTB) 
      KWRK1   = KTGB + N2BST(ISYMTB)
      LWRK1   = LWORK   - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insuff. work in CCSL_GTR 1')
      END IF
C
      CALL DZERO(WORK(KT1AMPA),NT1AM(ISYMTB))
      CALL DZERO(WORK(KTGB),N2BST(ISYMTB))
C
C     Read one electron excitation part of response vector
C
      IOPT = 1
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     *              WORK(KT1AMPA),WORK(KDUM))
C
      CALL CCSL_TGB(WORK(KT1AMPA),ISYMTB,LISTB,IDLSTB,WORK(KTGB),'ET',
     *              MODEL,WORK(KWRK1),LWRK1)
C
      LABEL = 'GIVE INT'
      CALL CCLR_FA(LABEL,ISYMTB,LISTC,IDLSTC,
     &             LISTA,IDLSTA,WORK(KTGB),WORK(KWRK1),LWRK1)
C
C
      LSAME  =  (LISTC.EQ.LISTB  .AND. IDLSTC.EQ.IDLSTB)
      IF (LSAME) THEN
        FACTSLV = TWO
      ELSE
        FACTSLV = ONE
      END IF
C
      CALL DAXPY(NT1AM(ISYCB),FACTSLV,WORK(KWRK1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYCB),FACTSLV,
     *           WORK(KWRK1+NT1AM(ISYCB)),1,RHO2,1)
C
      IF (.NOT. (LSAME)) THEN
C
C       second we consider the effective operator T^gC
C       -----------------------------------------------------
C
        KT1AMPA = 1
        KTGC    = KT1AMPA + NT1AM(ISYMTC)
        KWRK1   = KTGC + N2BST(ISYMTC)
        LWRK1   = LWORK   - KWRK1
C
        IF (LWRK1 .LT. 0) CALL QUIT('Insuff. work in CCSL_GTR 2')
C
        CALL DZERO(WORK(KT1AMPA),NT1AM(ISYMTC))
        CALL DZERO(WORK(KTGC),N2BST(ISYMTC))
C
C       Read one electron excitation part of response vector
C
        IOPT = 1
        CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     *                WORK(KT1AMPA),WORK(KDUM))
C
        CALL CCSL_TGB(WORK(KT1AMPA),ISYMTC,LISTC,IDLSTC,WORK(KTGC),'ET',
     *                MODEL,WORK(KWRK1),LWRK1)
C
        LABEL = 'GIVE INT'
        CALL CCLR_FA(LABEL,ISYMTC,LISTB,IDLSTB,
     &               LISTA,IDLSTA,WORK(KTGC),WORK(KWRK1),LWRK1)
C
C
        CALL DAXPY(NT1AM(ISYCB),ONE,WORK(KWRK1),1,RHO1,1)
        CALL DAXPY(NT2AM(ISYCB),ONE,
     *             WORK(KWRK1+NT1AM(ISYCB)),1,RHO2,1)
C
      END IF
C
C     finally, we consider the effective operator T^g sigma
C     -----------------------------------------------------
C
      KTGCB = 1
      KWRK1 = KTGCB + N2BST(ISYCB)
      LWRK1 = LWORK  - KWRK1
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_GTR, 3')
C
      CALL DZERO(WORK(KTGCB),N2BST(ISYCB))
C
      CALL CCSL_TGCB(ISYRES,LISTA,LISTB,LISTC,
     *               IDLSTA,IDLSTB,IDLSTC,
     *               ISYMTA,ISYMTB,ISYMTC,WORK(KTGCB),
     *               MODEL,WORK(KWRK1),LWRK1) 
C
      NAMPF = NT1AM(ISYCB) + NT2AM(ISYCB)
C
      KETA    = KWRK1
      KWRK2   = KETA + NAMPF
      LWRK2   = LWORK  - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_GTR, 4')
C
      CALL DZERO(WORK(KETA),NAMPF)
C
      LABEL = 'GIVE INT'
      CALL CC_ETAC(ISYCB,LABEL,WORK(KETA),
     *             LISTA,IDLSTA,0,WORK(KTGCB),WORK(KWRK2),LWRK2)
C
      KETA1   = KETA
      KETA2   = KETA1 + NT1AM(ISYCB) 
C
      CALL DAXPY(NT1AM(ISYCB),ONE,WORK(KETA1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYCB),ONE,WORK(KETA2),1,RHO2,1)
C
      END
**************************************************************
C  /* Deck ccsl_tgcb */
      SUBROUTINE CCSL_TGCB(ISYRES,LISTA,LISTB,LISTC,
     *                     IDLSTA,IDLSTB,IDLSTC,
     *                     ISYMTA,ISYMTB,ISYMTC,TGCB,
     *                     MODEL,WORK,LWORK)
C
C    Constructs effective operator equal to
C    T^gCB = -2 sum_lm g_l <Lambda|[[T_lm,C],B]|CC>T_lm
C
C    JK+OC, 03
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK)
      DIMENSION TGCB(*)
C
      INTEGER KDUM
      INTEGER IDLSTA,IDLSTB,IDLSTC,ISYMTC,ISYMTA,ISYMTB,ISYRES
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      CHARACTER*(*) LISTA,LISTB,LISTC
      CHARACTER*10 MODEL
C
      CHARACTER*8 LABEL
      LOGICAL LEXIST
C
      DOUBLE PRECISION GL
      EXTERNAL GL
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCSL_TGCB: ISYRES =',ISYRES
        WRITE(LUPRI,*)'CCSL_TGCB: ISYTMA  =',ISYMTA
        WRITE(LUPRI,*)'CCSL_TGCB: ISYTMB  =',ISYMTB
        WRITE(LUPRI,*)'CCSL_TGCB: ISYTMC  =',ISYMTC
        WRITE(LUPRI,*)'CCSL_TGCB: LISTA  =',LISTA
        WRITE(LUPRI,*)'CCSL_TGCB: LISTB  =',LISTB
        WRITE(LUPRI,*)'CCSL_TGCB: LISTC  =',LISTC
        WRITE(LUPRI,*)'CCSL_TGCB: IDLSTA =',IDLSTA
        WRITE(LUPRI,*)'CCSL_TGCB: IDLSTB =',IDLSTB
        WRITE(LUPRI,*)'CCSL_TGCB: IDLSTC =',IDLSTC
        CALL FLSHFO(LUPRI)
      ENDIF
C--------------------------------------
C     Readin trial vector (B) from file
C--------------------------------------
C
      KT1AMB = 1
      KT2AMB = KT1AMB + NT1AM(ISYMTB)
      KWRK1  = KT2AMB + NT2AM(ISYMTB)
      LWRK1  = LWORK - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insuff. work in CCSL_TGCB 1')
      END IF
C
      IOPT = 3
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     *              WORK(KT1AMB),WORK(KT2AMB))
C
C------------------
C        Symmetry 1
C------------------
C
      ISYCB = MULD2H(ISYMTB,ISYMTC)
C
      IF (ISYCB .NE. ISYRES) THEN
        CALL QUIT( 'Symmetry problem in CCSL_TGCB')
      END IF
C
C---------------------------------
C     Read in integrals from  file.
C---------------------------------
C
      LUCSOL = -1
      CALL GPOPEN(LUCSOL,'AOSOLINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCSOL)
      CALL MOLLAB('SOLVRLM ',LUCSOL,LUERR)
      READ (LUCSOL)
      CALL READT(LUCSOL,NLMCU,WORK(KWRK1))
C
      TEMP = 0.00
      LM = 0
      DO 600 L = 0,LMAXCU
         READ (LUCSOL) L1,(ISYTLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
           WRITE (LUERR,*)
     *     'ERROR CCSL_TGCB: L from LUCSOL not as expected'
           WRITE (LUERR,*) 'L from 600 loop:',L
           WRITE (LUERR,*) 'L from LUCSOL   :',L1
           CALL QUIT('ERROR CCSL_TGCB: L from LUCSOL not as expected')
         END IF
C
        DO 500 M = -L,L
         LM = LM + 1
         IF (IPRINT .GE. 15) THEN
           WRITE(LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                               ' ===================='
           WRITE(LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYTLM(L+M+1))
         END IF
C
C------------------
C        Symmetry 2 
C------------------
C
         IF (ISYTLM(L+M+1) .NE. ISYCB) THEN
            READ (LUCSOL)
            GO TO 500
         ENDIF
C
C--------------------------------
C        Read T(l,m) in AO basis.
C--------------------------------
C
         KTLMAO = KWRK1
         KWRK2  = KTLMAO + N2BST(ISYCB)
         LWRK2  = LWORK - KWRK2
         IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_TGCB, 2')
C
         CALL CC_SLVINT(WORK(KTLMAO),WORK(KWRK2),LWRK2,ISYCB)
C
         IF (IPRINT .GT. 25) THEN
            CALL AROUND('CCSL_TGCB : Tlm_ao matrix: - cc storage')
            CALL CC_PRFCKAO(WORK(KTLMAO),ISYCB)
         ENDIF
C
C-------------------------------------------------------------------
C
        LABEL = 'GIVE INT'
        CALL CCLR_FA(LABEL,ISYCB,LISTC,IDLSTC,
     &               LISTA,IDLSTA,WORK(KTLMAO),WORK(KWRK2),LWRK2)
C
C
        CBDOT1 = DDOT(NT1AM(ISYMTB),WORK(KWRK2),1,WORK(KT1AMB),1)
        CBDOT2 = DDOT(NT2AM(ISYMTB),WORK(KWRK2+NT1AM(ISYMTB)),1,
     *                              WORK(KT2AMB),1)
C
        CBDOT = CBDOT1 + CBDOT2
C
C-------------------------
C        Find  G_lm factor
C-------------------------
C
         GLFAC =  GL(L,EPOPCU,RCAVCU)
         FACT =   2.0D0*GLFAC*CBDOT
C
C--------------------------------------------------------------
C        Put factor on integrals and (for a given l,m) the
C        effective operator has now been constructed.
C--------------------------------------------------------------
C
         CALL DAXPY(N2BST(ISYCB),FACT,WORK(KTLMAO),1,
     *              TGCB,1)
C---------------------------------------------------------------

  500   CONTINUE
  600 CONTINUE
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('ccsl_tgb : T^gCB_ao matrix: ')
         CALL CC_PRFCKAO(TGCB,ISYCB)
      ENDIF
C
      CALL GPCLOSE(LUCSOL,'KEEP')
C
      END
**********************************************************************
C  /* Deck ccsl_btr */
      SUBROUTINE CCSL_BTR(RHO1,RHO2,ISYMAB,
     *                    LISTA,IDLSTA,ISYMA,
     *                    LISTB,IDLSTB,ISYMB,
     *                    MODEL,WORK,LWORK)
C
C     JK+OC, 03
C---------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*)
C
      CHARACTER*(*) LISTA,LISTB
      CHARACTER*(3) LISTC
      CHARACTER*8 LABEL
      CHARACTER*10 MODEL
      LOGICAL LEXIST, LSAME
      INTEGER IDLSTA, IDLSTB, ISYMAB, ISYMA, ISYMB 
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCSL_BTR: ISYMAB =',ISYMAB
        WRITE(LUPRI,*)'CCSL_BTR: ISYMA  =',ISYMA
        WRITE(LUPRI,*)'CCSL_BTR: ISYMB  =',ISYMB
        WRITE(LUPRI,*)'CCSL_BTR: LISTA  =',LISTA
        WRITE(LUPRI,*)'CCSL_BTR: LISTB  =',LISTB
        WRITE(LUPRI,*)'CCSL_BTR: IDLSTA =',IDLSTA
        WRITE(LUPRI,*)'CCSL_BTR: IDLSTB =',IDLSTB
        CALL FLSHFO(LUPRI)
      ENDIF
C
C     First we concentrate on the effective operator T^gA
C     -----------------------------------------------------
C
      KT1AMPA = 1
      KTGA    = KT1AMPA + NT1AM(ISYMA)
      KWRK1   = KTGA + N2BST(ISYMA)
      LWRK1   = LWORK   - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insuff. work in CCSL_BTR 1')
      END IF
C
      CALL DZERO(WORK(KT1AMPA),NT1AM(ISYMA))
      CALL DZERO(WORK(KTGA),N2BST(ISYMA))
C
C     Read one electron excitation part of response vector
C
      IOPT = 1
      CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     *              WORK(KT1AMPA),WORK(KDUM))
C
      CALL CCSL_TGB(WORK(KT1AMPA),ISYMA,LISTA,IDLSTA,WORK(KTGA),'ET',
     *              MODEL,WORK(KWRK1),LWRK1)
C
      LABEL = 'GIVE INT'
      CALL CCCR_AA(LABEL,ISYMA,LISTB,IDLSTB,
     &            WORK(KTGA),WORK(KWRK1),LWRK1)
      CALL CCLR_DIASCL(WORK(KWRK1+NT1AM(ISYMAB)),2.0d0,ISYMAB)
C
C
      LSAME  =  (LISTA.EQ.LISTB  .AND. IDLSTA.EQ.IDLSTB)
      IF (LSAME) THEN
        FACTSLV = TWO
      ELSE
        FACTSLV = ONE
      END IF
C
      CALL DAXPY(NT1AM(ISYMAB),FACTSLV,WORK(KWRK1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYMAB),FACTSLV,
     *           WORK(KWRK1+NT1AM(ISYMAB)),1,RHO2,1)
C
C
      IF (.NOT.(LSAME)) THEN
C
C       Next we concentrate on the effective operator T^gB
C       -----------------------------------------------------
C
        KT1AMPB = 1
        KTGB    = KT1AMPB + NT1AM(ISYMB)
        KWRK1   = KTGB + N2BST(ISYMB)
        LWRK1   = LWORK   - KWRK1
C
        IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insuff. work in CCSL_BTR 2')
        END IF
C
        CALL DZERO(WORK(KT1AMPB),NT1AM(ISYMB))
        CALL DZERO(WORK(KTGB),N2BST(ISYMB))
C
C       Read one electron excitation part of response vector
C
        IOPT = 1
        CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     *                WORK(KT1AMPB),WORK(KDUM))

        CALL CCSL_TGB(WORK(KT1AMPB),ISYMB,LISTB,IDLSTB,WORK(KTGB),'ET',
     *                MODEL,WORK(KWRK1),LWRK1)
C
        LABEL = 'GIVE INT'
        CALL CCCR_AA(LABEL,ISYMB,LISTA,IDLSTA,
     &              WORK(KTGB),WORK(KWRK1),LWRK1)
        CALL CCLR_DIASCL(WORK(KWRK1+NT1AM(ISYMAB)),2.0d0,ISYMAB)
C
C
        CALL DAXPY(NT1AM(ISYMAB),ONE,WORK(KWRK1),1,RHO1,1)
        CALL DAXPY(NT2AM(ISYMAB),ONE,
     *             WORK(KWRK1+NT1AM(ISYMAB)),1,RHO2,1)
C 
      END IF
C
C     Finally, we concentrate on the effective operator T^gAB
C     -------------------------------------------------------
C
      KTGAB = 1
      KWRK1 = KTGAB + N2BST(ISYMAB)
      LWRK1 = LWORK  - KWRK1
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_BTR, 3')
C
      CALL DZERO(WORK(KTGAB),N2BST(ISYMAB))
C
      LISTC  = 'L0'
      IDLSTC =  1
      ISYMC  = ILSTSYM(LISTC,IDLSTC)
C
      CALL CCSL_TGCB(ISYMAB,LISTC,LISTA,LISTB,
     *               IDLSTC,IDLSTA,IDLSTB,
     *               ISYMC,ISYMA,ISYMB,WORK(KTGAB),
     *               MODEL,WORK(KWRK1),LWRK1)
C
      NAMPF = NT1AM(ISYMAB) + NT2AM(ISYMAB)
C
      KXI   = KWRK1
      KWRK2 = KXI + NAMPF
      LWRK2 = LWORK  - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_BTR, 4')
C
      CALL DZERO(WORK(KXI),NAMPF)
C
      LABEL = 'GIVE INT'
      CALL CC_XKSI(WORK(KXI),LABEL,ISYMAB,0,WORK(KTGAB),
     *             WORK(KWRK2),LWRK2)
C
      KXI1 = KXI
      KXI2 = KXI1 + NT1AM(ISYMAB)
C
      CALL CCLR_DIASCL(WORK(KXI2),2.0d0,ISYMAB)
C
      CALL DAXPY(NT1AM(ISYMAB),ONE,WORK(KXI1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYMAB),ONE,WORK(KXI2),1,RHO2,1)

      END
********************************************************************
C  /* Deck ccsl_xieta */
      SUBROUTINE CCSL_XIETA(WORK,LWORK)
C
C-----------------------------------------------------------------
C JK, april .03
C
C Calculates the CCSL ETA and XI vectors and write them to files.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK)
C
      INTEGER LUMMET, LUMMXI, ISYMTR
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      CHARACTER*(8) FILMME, FILMMX, LABEL
      CHARACTER*2 LIST
C
      LUMMET = -1
      LUMMXI = -1
      FILMME = 'CCSL_ETA'
      FILMMX = 'CCSL__XI'
C
      CALL WOPEN2(LUMMET, FILMME, 64, 0)
      CALL WOPEN2(LUMMXI, FILMMX, 64, 0)

C---------------------------------
C     Read in integrals from  file.
C---------------------------------
C
      LUCSOL = -1
      CALL GPOPEN(LUCSOL,'AOSOLINT','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCSOL)
      CALL MOLLAB('SOLVRLM ',LUCSOL,LUERR)
      READ (LUCSOL)
      CALL READT(LUCSOL,NLMCU,WORK(1))
C
      TEMP = 0.00
      LM = 0
      IADR1 = 1
      IADR2 = 1
C
      DO 600 L = 0,LMAXCU
         READ (LUCSOL) L1,(ISYTLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
           WRITE (LUERR,*)
     *     'ERROR CCSL_XIETA: L from LUCSOL not as expected'
           WRITE (LUERR,*) 'L from 600 loop:',L
           WRITE (LUERR,*) 'L from LUCSOL   :',L1
           CALL QUIT('ERROR CCSL_XIETA: L from LUCSOL not as expected')
         END IF
C
        DO 500 M = -L,L
          LM = LM + 1
C          IF (IPRINT .GE. 15) THEN
            WRITE(LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                                ' ===================='
            WRITE(LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYTLM(L+M+1))
C          END IF
C
C
C         Symmetry of the T_lm operator:
C         ------------------------------
C 
          ISYMTR = ABS(ISYTLM(L+M+1))
C
C         Number of amplitudes in the given symmetry: 
C         -------------------------------------------
C
          NTAMP1  = NT1AM(ISYMTR)
          NTAMP2  = NT2AM(ISYMTR)
          IF (CCS)  NTAMP2  = 0
          NTAMP   = NTAMP1 + NTAMP2
C
          KTLMAO = 1
          KWRK1  = KTLMAO + N2BST(ISYMTR)
          LWRK1  = LWORK - KWRK1
          IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCSL_XIETA, 1')
          CALL DZERO(WORK(KTLMAO),N2BST(ISYMTR))
C
          CALL CC_SLVINT(WORK(KTLMAO),WORK(KWRK1),LWRK1,ISYMTR)
C
          IF (IPRINT .GT. 25) THEN
             CALL AROUND('CCSL_XIETA : Tlm_ao matrix: - cc storage')
             CALL CC_PRFCKAO(WORK(KTLMAO),ISYMTR)
          ENDIF
C
C         Allocate space for the result vector:
C         -------------------------------------
C
          KXI     = KWRK1
          KWRK2   = KXI     + NTAMP
          LWRK2   = LWORK   - KWRK2
          IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCSL_XIETA, 2')
          CALL DZERO(WORK(KXI),NTAMP)
C
C         Perform the calculation
C         -----------------------
C
          LABEL = 'GIVE INT'
          LIST  = 'L0'
          IDLINO = 1
C
          CALL CC_ETAC(ISYMTR,LABEL,WORK(KXI),
     *                 LIST,IDLINO,0,WORK(KTLMAO),WORK(KWRK2),LWRK2)
          CALL PUTWA2(LUMMET,FILMME,WORK(KXI),IADR1,NTAMP)
          IADR1 = IADR1 + NTAMP
C
          CALL CC_XKSI(WORK(KXI),LABEL,ISYMTR,0,WORK(KTLMAO),
     *                 WORK(KWRK2),LWRK2)
          CALL PUTWA2(LUMMXI,FILMMX,WORK(KXI),IADR2,NTAMP)
          IADR2 = IADR2 + NTAMP
C
  500   CONTINUE
  600 CONTINUE
C
      CALL GPCLOSE(LUCSOL,'KEEP')
C
      CALL WCLOSE2(LUMMET, FILMME, 'KEEP')
      CALL WCLOSE2(LUMMXI, FILMMX, 'KEEP')
C
      END
C***********************************************************************
