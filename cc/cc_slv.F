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
C  /* Deck cc_slv */
       SUBROUTINE CC_SLV(AODEN,ETLM,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              solvent effects.
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
      IF (LSLECVG.AND.LSLTCVG.AND.LSLLCVG) THEN
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
            IF (ABS(ETLM(LM,2)) .GT. THRZER) THEN
              WRITE(LUPRI,*) 'ERROR CC_SLVE for l,m',L,M
              WRITE(LUPRI,*) 'Symmetry :',ISYTLM(L+M+1)
              WRITE(LUPRI,*) 'Tn(l,m) .ne. 0, but =',ETLM(LM,2)
C              CALL QUIT( 'ERROR CC_SLVE: Tn(l,m) not 0 as expected')
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
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*) 
     *    'CCSL_RHSTG: Solvent contribution to CC equations.'
      ENDIF
C
C---------------------------------
C     Readin integrals from  file.
C---------------------------------
C
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
         RHO1N = DDOT(NT1AM(ISYMTR),CTR1,1,CTR1,1)
         RHO2N = DDOT(NT2AM(ISYMTR),CTR2,1,CTR2,1)
         WRITE(LUPRI,*) ' Norm af C1AM in CCSL_LTGB on input:', RHO1N
         WRITE(LUPRI,*) ' Norm af C2AM in CCSL_LTGB on input:', RHO2N
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
C---------------------------------
C     Readin integrals from  file.
C---------------------------------
C
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
C----------------------------------------------------------
C        Symmetry should be equal to trial vector symmetry.
C----------------------------------------------------------
C
         IF (ISYTLM(L+M+1) .NE. ISYMTR) THEN
            READ (LUCSOL)
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
         LABEL = 'GIVE INT'           
         IF ((LR.EQ.'L').OR.(LR.EQ.'0').OR.(LR.EQ.'P')) THEN
            CALL CC_XKSI(WORK(KXI),LABEL,ISYMTR,0,WORK(KTLMAO),
     *                   WORK(KWRK3),LWRK3)
         ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
            LIST  = 'L0'
            ILSTNR  = 1
            CALL CC_ETAC(ISYMTR,LABEL,WORK(KXI),
     *                   LIST,ILSTNR,0,WORK(KTLMAO),WORK(KWRK3),LWRK3)
         ENDIF
C
C---------------------------
C        3. Contract with B
C---------------------------
C
         IF (LR.NE.'0') THEN
           KXI1 = KXI  
           KXI2 = KXI  + NTAMP1
C wrong:
c          IF (LR.EQ.'R')
c    *         CALL CCLR_DIASCL(CTR2,TWO,ISYMTR)
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
C  /* Deck ccsl_tlmput*/
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
      CALL GPOPEN(LUCTLM,'CC_TL','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      LM = 0
      DO 820 L = 0,LMAXCU
        DO 800 M = -L,L
          LM = LM + 1
          READ(LUCTLM,'(I5,3X,E22.15)',END=100) LM1,CCTLM(LM) 
c         write(lupri,*) 'lm,CCTLM(LM) read',lm,CCTLM(LM)
          IF (LM1.NE.LM) CALL QUIT( 'Error on CC_TL')
          GOTO 800
  100     CCTLM(LM) = 0.0D0   
c         write(lupri,*) 'lm,CCTLM(LM) set 0',lm,CCTLM(LM)
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
      CALL GPOPEN(LUCTLM,'CC_TL','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
C
      LM = 0
      DO 920 L = 0,LMAXCU
        DO 900 M = -L,L
          LM = LM + 1
          WRITE(LUCTLM,'(I5,3X,E22.15)') LM,CCTLM(LM) 
  900   CONTINUE
  920 CONTINUE
C
      CALL GPCLOSE(LUCTLM,'KEEP')
C
      END
