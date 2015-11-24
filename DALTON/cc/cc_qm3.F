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
C   /* Deck cc_qm3 */
       SUBROUTINE CC_QM3(AODEN,CONVERGED,WORK,LWORK)
C*********************************************************************
C
C      Purpose: Calculation of CC/MM wave-function and energy.
C      CCS(CIS/HF)(nci), MP2(nci), CC2(nci), CCSD, CC3(nci), MCC2(nci)
C      QM3, JK+AO+OC+KMI,01 
C
C*********************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
C
      LOGICAL   CONVERGED
      DIMENSION WORK(LWORK),AODEN(*),FFJ(3)
      CHARACTER MODEL*10
      CHARACTER MODELPRI*4
      CHARACTER LABEL*8
C
#if defined (SYS_CRAY)
      REAL              DEMMMM
#else
      DOUBLE PRECISION  DEMMMM
#endif
C
C--------------------------------------------------
C     Header for the CCMM output in each iteration:
C--------------------------------------------------
C
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*****************************************************'//
     *'************'
      WRITE (LUPRI,'(1X,A)')
     *'**** Output from Coupled Cluster/Molecular Mechanics program'//
     *' ****'
      WRITE (LUPRI,'(1X,A)')
     *'*****************************************************'//
     *'************'
C
C-----------------------------------------------
C     Testing if the CC model is CC2 or CCSD:
C-----------------------------------------------
C
      MODEL = 'CCSD'
C
      IF (CCSD) THEN
        CALL AROUND('Coupled Cluster model is: CCSD')
        MODEL = 'CCSD'
        MODELPRI = 'CCSD'
      ENDIF
C
      IF (CC2.AND.(.NOT.MCC2)) THEN
        CALL AROUND('Coupled Cluster model is: CC2')
        MODEL = 'CC2'
        MODELPRI = ' CC2'
      ENDIF
C
      IF (.NOT. (CCSD.OR.CC2)) THEN
        CALL QUIT('CCMM only implemented for the ' //
     &            'CC2 and CCSD parameterizations.') 
      ENDIF
C
      IF (IQM3PR .GT .10) 
     &WRITE(LUPRI,*) 'CC_MM-1: Workspace:',LWORK
C
C-----------------------------------------------
C     Calculate the E(QM/MM) part of the energy:
C-----------------------------------------------
C
      ECMCON = ZERO
         
      IF (.NOT. NYQMMM) THEN
        CALL CC_QM3E(ECMCON,EPOL,EEL,EELEL,EREP,AODEN,
     &             WORK,LWORK)
      ELSE IF (NYQMMM) THEN
        CALL CC_QMMME(ECMCON,EPOL,EEL,EELEL,EREP,AODEN,
     &             WORK,LWORK)
      END IF
C
C----------------------------
C     Calculate new E(QM/MM):
C----------------------------
C
      IF (.NOT. (OLDTG) .AND. .NOT. (NYQMMM)) THEN
        ECMCU = (ECCGRS - EELEL) + ECMCON 
        ECCGRS1 = (ECCGRS - EELEL)
      ELSE
        ECMCU = ECCGRS + ECMCON
        ECCGRS1 = ECCGRS
      END IF 
C
      IF (ABS(ECMCU-ECCPR).LT.CVGESOL) LSLECVG = .TRUE.
C
      WRITE(LUPRI,'(//1X,A,I3,A,F20.10)')
     *      'E(QM/MM) contribution in iteration',ICCSLIT,':    ',
     *       ECMCON
C
      WRITE(LUPRI,'(1X,A,F20.10)')
     *      'CC energy in the  current CCMM iteration: ',ECMCU  
C
      WRITE(LUPRI,'(1X,A,F20.10)')
     *      'CC energy in the previous CCMM iteration: ',ECCPR  
C
      WRITE(LUPRI,'(1X,A,E20.6//)')
     *      'Change in Total energy in this CCMM it.:  ',
     *       ECMCU - ECCPR 
C
      ECCPR   = ECMCU  
C
      CONVERGED = .FALSE.
      IF (LSLECVG.AND.LSLTCVG.AND.LSLLCVG) THEN
        CONVERGED = .TRUE.
        IF (LOEC3) THEN
C
          EPOL = ZERO
C
          WRITE(LUPRI,'(1X,A)')
     &    'Perturbative correction to polarization energy is calculated'
     &   ,'(Model SPC_EC3)'
C
          CALL EPERT2(EPOL,WORK,LWORK)
        END IF
C
C       Here we see if any static fields are to be added
C
        FFJ(1) = 0.0D0
        FFJ(2) = 0.0D0
        FFJ(3) = 0.0D0
C
        IF (RELMOM) THEN
          DO 330 IF =1, NFIELD
            IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
            IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
            IF (LFIELD(IF) .EQ. 'ZDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
  330     CONTINUE
        END IF
C
        CALL QM3_EMMMM(DEMMMM,TEMPOL,FFJ)
C
        IF (.NOT.(RELMOM)) PEDIP1 = 0.0D0
C
        EMMPOL = TEMPOL
        IF (.NOT.(RELMOM)) THEN
          EMM_MM = EMMELC + EMMVDW + EMMPOL
        ELSE
          EMM_MM = EMMELC + EMMVDW + EMMPOL + PEDIP1
        END IF
C
C---------------------------------------------------
C       Final output for the QM/MM results with the 
C       converged wave function:
C---------------------------------------------------
C
C       First the MM/MM energy with the QM induced 
C       dipole moments:
C-------------------------------------------------
C
        IF (.NOT. NYQMMM) THEN
         WRITE(LUPRI,'(//1X,72A)') ('-',J=1,72)
         WRITE(LUPRI,'(1X,A56,A16)')
     * '| --- The MM/MM classical interaction energy ',
     * '--- |'
         WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
         WRITE(LUPRI,'(1X,A1,A50,A7,F12.8,A2)') '|',
     *  '  Eelec = Sum_n,s[ (Q_n*Q_s)/|R_n - R_s| ]        ',
     *  '|      ',EMMELC,' |'
         IF (RELMOM) THEN
          WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
          WRITE(LUPRI,'(1X,A1,A50,A7,F12.8,A2)') '|',
     *  '  Epol  = - 1/2*Sum_a[ MYind_a*(E^s_a + E^ext)]         ',
     *  '|      ',EMMPOL,' |'
          WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
          WRITE(LUPRI,'(1X,A1,A50,A7,F12.8,A2)') '|',
     *  '  E_field = - Mu_perm* E^ext                            ',
     *  '|      ',PEDIP1,' |'
         ELSE
          WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
          WRITE(LUPRI,'(1X,A1,A50,A7,F12.8,A2)') '|',
     *  '  Epol  = - 1/2*Sum_a[ MYind_a*E^site_a ]         ',
     *  '|      ',EMMPOL,' |'
         END IF 
         WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
         WRITE(LUPRI,'(1X,A1,A50,A7,F12.8,A2)') '|',
     *  '  Evdw  = Sum_a[ A_ma/|R_ma|^12 - B_ma/|R_ma|^6 ] ',
     *  '|      ',EMMVDW,' |'
         WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
         WRITE(LUPRI,'(1X,A1,50A,A1,A20)')
     *              '|',('*',J=1,50),'|',
     *              '*******************|'
         WRITE(LUPRI,'(1X,72A)') ('-',J=1,72)
         IF (.NOT.(RELMOM)) THEN
          WRITE(LUPRI,'(1X,A1,A50,A7,F12.8,A2)') '|',
     *    '  E(MM/MM) = Eelec + Epol + Evdw                  ',
     *    '|      ',EMM_MM,' |'
          WRITE(LUPRI,'(1X,72A,//)') ('-',J=1,72)
         ELSE
          WRITE(LUPRI,'(1X,A1,A50,A7,F12.8,A2)') '|',
     *    '  E(MM/MM) = Eelec + Epol + E_field + Evdw         ',
     *    '|      ',EMM_MM,' |'
          WRITE(LUPRI,'(1X,72A,//)') ('-',J=1,72)
         END IF
        END IF
C
C--------------------------------------------
C  Second the CC/MM interaction energy terms:
C--------------------------------------------
C
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('=',J=1,70),'+'
C
        IF (MODEL .EQ. 'CC2')THEN
          WRITE(LUPRI,'(1X,A1,A15,A18,A3,A19,A15,A1)')
     *         '|','--------------- ',
     *         'Final output from ',MODEL,'/MM energy program ',
     *         '---------------','|'
        ELSE IF (MODEL .EQ. 'CCSD') THEN
          WRITE(LUPRI,'(1X,A1,A15,A18,A4,A19,A14,A1)')
     *         '|','-------------- ',
     *         'Final output from ',MODEL,'/MM energy program ',
     *         '--------------','|'
        END IF
C
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('=',J=1,70),'+'
C
        WRITE(LUPRI,'(1X,A2,A15,A3,A15,A3,A15,A3,A14,A2)')
     *        '| ','     Eelec     ',' | ','     Epol      ',
     *        ' | ','     Evdw      ',' | ','   E(QM/MM)   ',
     *        ' |'
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('-',J=1,70),'+'
        WRITE(LUPRI,'(1X,A2,F15.10,A3,F15.10,A3,F15.10,
     *                A3,F14.10,A2)')
     *        '| ',EEL,' | ',EPOL,' | ',ECLVDW,' | ',ECMCON,
     *        ' |'
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('-',J=1,70),'+'
C
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('=',J=1,70),'+'
C
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('-',J=1,70),'+'
        WRITE(LUPRI,'(1X,A2,A15,A3,A15,A3,A15,A3,A14,A2)')
     *        '| ',' <L|H(vac)|CC> ',' | ','<H(qm)+H(qmmm)>',
     *        ' | ','Delta E(mm/mm) ',' | ','E_rep         ',
     *        ' |'
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('-',J=1,70),'+'
        WRITE(LUPRI,'(1X,A2,F15.10,A3,F15.10,A3,F15.10,
     *                A3,F14.10,A2)')
     *        '| ',ECCGRS1,' | ',ECMCU,' | ',DEMMMM,' | ',EREP,
     *        ' |'
        WRITE(LUPRI,'(1X,A1,70A,A1//)') '+',('=',J=1,70),'+'
C
        IF (RELMOM) THEN
          WRITE(LUPRI,*)' '
          WRITE(LUPRI,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(LUPRI,*)' '
          WRITE(LUPRI,*)'Pind_a is determined without external field'
          WRITE(LUPRI,*)'Delta E(mm/mm) is WRONG !!!                '
          WRITE(LUPRI,*)' '
          WRITE(LUPRI,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(LUPRI,*)' '
        END IF
C------------------------------------------------------
C       Third information concerning the convergence is 
C       written: 
C------------------------------------------------------
C
        IF (LOITER) THEN
          WRITE(LUPRI,'(1X,A,I3,A)')
     *          'Maximum inner iterations for t set to ',
     *           MXTINIT,' in each outer iteration'
          WRITE(LUPRI,'(1X,A,I3,A/)')
     *          'Maximum inner iterations for t-bar set to ',
     *           MXLINIT,' in each outer iteration'
        END IF
C
        WRITE(LUPRI,'(/,1X,A,I3,A)')
     *  ' CCMM equations are converged in ',ICCSLIT,
     *  ' outer iterations'
        WRITE(LUPRI,'(1X,A,I3,A)')
     *  ' CCMM equations are converged in ',NSLVINIT,
     *  ' inner iterations'
C
        WRITE(LURES,'(12X,A4,A,F20.10)')
     *  MODELPRI,' Total  energy:             ',ECMCU
C
        WRITE(LURES,'(12X,A4,A,F20.10,/)')
     *  MODELPRI,' E(QM/MM)     :             ',ECMCON
C
        ECCGRS1 = ECMCU    
        LABEL = MODELPRI//'/MM '
        CALL WRIPRO(ECCGRS1,MODEL,0,
     *              LABEL,LABEL,LABEL,LABEL,
     *              ZERO,ZERO,ZERO,1,0,0,0)
C
        LABEL = 'E_QM/MM '
        CALL WRIPRO(ECMCON,MODEL,0,
     *              LABEL,LABEL,LABEL,LABEL,
     *              ZERO,ZERO,ZERO,1,0,0,0)

C------------------------------------------------------
C       Print information concerning the reaction field
C------------------------------------------------------

        CALL RFQMMM()

C------------------------------------------------------

      ELSE 
C
        ICCSLIT = ICCSLIT + 1
C
        IF (ICCSLIT.GT.MXCCSLIT) THEN
          WRITE(LUPRI,*) 'Maximum number of CCMM iterations',
     *                MXCCSLIT,' is reached'
          CALL QUIT( 'Maximum number of CCMM iterations reached')
        ENDIF
C
      ENDIF
C
C--------------------------------------
C     Ending header for each iteration:
C--------------------------------------
C
      WRITE (LUPRI,'(1X,A)')
     *'*****************************************************'//
     *'************'
      WRITE (LUPRI,'(1X,A)')
     *'******* End of Coupled Cluster/Molecular Mechanics program'//
     *' ******'
      WRITE (LUPRI,'(1X,A)')
     *'*****************************************************'//
     *'************'
C
       END
C
C**************************************************************
C  /* Deck cc_qm3e */
      SUBROUTINE CC_QM3E(ECCMM,EPOL,EEL,EELEL,EREP,AODEN,WORK,LWORK)
C**************************************************************
C
C     QM3: JK+AO+OC+KMI,01
C     Purpose: Only calculation of E(QM/MM)
C     QMMM: sneskov, 09
C     Purpose: Updates the induced dipoles besides calculating E(QM/MM)
C
C     In new QMMM implementation (NYQMMM) the strategy is as follows
C     1) Calculate F^sta=F^el+F^mul+F^nuc
C     2) Calculate induced dipole by transforming with relay 
C     matrix(MMMAT) or by iterative means (MMITER)
C     3) Put induced moments to file
C     4) Calculate Epol=-1/2sum_a\mu^ind_aF^sta_a
C     5) Calculate Eelesta=-sum_s<N_s>
C
C     Also the LGFILE keyword is reset such that a new G matrix is
C     written to file when TGT is called next time.
C
C**************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "ccinftap.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "thrzer.h"

C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 1.0D0/2.0D0)
      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)
C
      DIMENSION WORK(LWORK),AODEN(*)
      CHARACTER LABEL*8

      DIMENSION CCNS(MXQM3),ALPHAMM(MXQM3),FFJ(3)
      INTEGER KS
C
      CALL QENTER('CC_QM3E')
C------------------------------------------------------------------
C       Init parameters
C------------------------------------------------------------------

      ISYMOP = 1
C
C------------------------------------------------------------------
C       Dynamical allocation for CCMM 
C------------------------------------------------------------------
C

C------------------------------------------------------------------
C      Dynamical allocation. If all models are SPC (LOSPC = .TRUE.)
C      the electric field from the electrons are not needed. Then 
C      only space for Ns in AO basis is allocated.
C------------------------------------------------------------------
C

         IF ( .NOT. (LOSPC) ) THEN
C
            KRAAOx  = 1
            KRAAOy = KRAAOx + N2BST(ISYMOP)
            KRAAOz = KRAAOy + N2BST(ISYMOP)
            KOMMSx = KRAAOz + N2BST(ISYMOP)
            KOMMSy = KOMMSx + NCOMS
            KOMMSz = KOMMSy + NCOMS
            KNSAO =  KOMMSz + NCOMS
            KRAx = KNSAO + N2BST(ISYMOP)
            KRAy = KRAx + NCOMS
            KRAz = KRAy + NCOMS
            KWRK1 = KRAz + NCOMS
            LWRK1 = LWORK - KWRK1
C
            CALL DZERO(WORK(KRAAOx),3*N2BST(ISYMOP))
            CALL DZERO(WORK(KOMMSx),3*NCOMS)
            CALL DZERO(WORK(KNSAO),N2BST(ISYMOP))
            CALL DZERO(WORK(KRAx),3*NCOMS)
C
         ELSE
C
            KNSAO = 1
            KWRK1 = KNSAO + N2BST(ISYMOP)
            LWRK1 = LWORK - KWRK1
C
            CALL DZERO(WORK(KNSAO),N2BST(ISYMOP))
C
         END IF
C
         IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CC_QM3RAINT')
C 
C----------------------------
C     Read Rr(a) in AO basis:
C----------------------------
C
C     First we see if any static fields are to be added
C
         FFJ(1) = 0.0D0
         FFJ(2) = 0.0D0
         FFJ(3) = 0.0D0
C
         IF (RELMOM) THEN
            DO 330 IF =1, NFIELD
               IF (LFIELD(IF).EQ. 'XDIPLEN') FFJ(1)=FFJ(1)+EFIELD(IF)
               IF (LFIELD(IF).EQ. 'YDIPLEN') FFJ(2)=FFJ(2)+EFIELD(IF)
               IF (LFIELD(IF).EQ. 'ZDIPLEN') FFJ(3)=FFJ(3)+EFIELD(IF)
  330       CONTINUE
         END IF
C
         IF (.NOT. LOSPC) THEN
            CALL CC_QM3RAINT(WORK(KRAAOx),WORK(KRAAOy),WORK(KRAAOz),
     &                   AODEN,WORK(KRAx),WORK(KRAy),WORK(KRAz),
     &                   WORK(KWRK1),LWRK1,ISYMOP)
         END IF
C
C------------------------------
C        Read N(s) in AO basis:
C------------------------------
C
         CALL CC_QM3NSINT(WORK(KNSAO),AODEN,CCNS,
     &                 WORK(KWRK1),LWRK1,ISYMOP)
C
C-------------------------------------
C     Add up the energy contributions:
C     1) Epol:
C-------------------------------------
C 
         ECCMM = ZERO
C
         IF (.NOT. LOSPC) THEN
            CALL CC_GET31('OBARFILE',NCOMS,
     *                     WORK(KOMMSx),WORK(KOMMSy),WORK(KOMMSz)) 
C
            KS = 0

            DO 110 I = 1, ISYTP
               IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
                  DO 120 J = NSYSBG(I), NSYSED(I)
                     DO 121 K = 1, NUALIS(I)
                        KS = KS+1
                        RAT = ZERO
                        RAT =  0.5D0 * ALPIMM(I,K)*WORK(KRAx+KS-1) *
     &                  (WORK(KRAx + KS -1) + WORK(KOMMSx + KS - 1))
     &                  + 0.5D0 * ALPIMM(I,K) * WORK(KRAy + KS -1) *
     &                  (WORK(KRAy + KS -1) + WORK(KOMMSy + KS - 1))
     &                  + 0.5D0 * ALPIMM(I,K) * WORK(KRAz + KS -1) *
     &                  (WORK(KRAz + KS -1) + WORK(KOMMSz + KS - 1))
                        ECCMM = ECCMM - RAT
  121                CONTINUE 
  120             CONTINUE
               END IF
  110       CONTINUE
C
            TEMP1 = ECCMM
            CALL QM3_OTILDE(OTILDE,FFJ)
         ELSE
            OTILDE = 0.0D0
            TEMP1 = 0.0D0
         END IF
C
         ECCMM = ECCMM + OTILDE
         EPOL = ECCMM
C
C--------
C 2) Eel:
C--------
C
         ETEMP = 0.0D0
         L = 0
C
         DO 130 I = 1, ISYTP
            IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
               DO 140 J = NSYSBG(I), NSYSED(I) 
                  DO 150 K = 1, NSISY(I)
                     L = L + 1
                     ETEMP = ETEMP - CCNS(L) 
  150             CONTINUE  
  140          CONTINUE
            END IF
  130    CONTINUE
C
         EELEL = ETEMP

         ECCMM = ECCMM + ENUQM3 + EELEL
C
         EEL = ENUQM3 + EELEL
C
C-------------
C 3) E(VDW):
C-------------
C
         ECCMM = ECCMM + ECLVDW
C
C-----------------------------------------------
C 4) E_repulsion
C
         CALL CCMM_REP2(EREP,AODEN,WORK(KWRK1),LWRK1)
         ECCMM = ECCMM + EREP
C
C-----------------------------------------------
C   Testing the energy contributions one by one:
C-----------------------------------------------
C
         IF (IQM3PR .GT. 5) THEN
            WRITE(LUPRI,'(//,A)')
            WRITE(LUPRI,'(/A)')
     *     ' +====================================='
     *     //'======================================'
     *     //'==================================+'
            WRITE(LUPRI,'(1X,A)')
     *    '|--------- '
     *     //'The different energy contributions outlined'
     *     //' ---------|'
            WRITE(LUPRI,'(A)')
     *    ' +====================================='
     *     //'======================================'
     *     //'==================================+'
            WRITE(LUPRI,'(1X,A)')
     *    '| Evdw                | E-nuclear           |'
     *     //' -Sum_s <N_s>        |'
     *     //'| -alpha*Sum_a[<Rra> | O-tilde             |'
            WRITE(LUPRI,'(1X,A)')
     *    '| contribution        | contribution        |'
     *     //'                     |'
     *     //'*{<Rra>+OmmS}]       | contribution        |'
            WRITE(LUPRI,'(A)')
     *    ' +-------------------------------------'
     *     //'--------------------------------------'
     *     //'----------------------------------+'
            WRITE(LUPRI,'(1X,A,F20.16,A,F20.16,A,F20.16,A,
     &             F20.16,A,F20.16,A)')
     *    '|',ECLVDW,' |',ENUQM3,' |',EELEL,' |',TEMP1,' |',
     *   OTILDE,' |'
            WRITE(LUPRI,'(A)')
     *    ' +====================================='
     *     //'======================================'
     *     //'==================================+'
            WRITE(LUPRI,'(//,A)')
         END IF
      
      CALL QEXIT('CC_QM3E')
C
      END
C
C**************************************************************
C  /* Deck cc_qm3raint */
      SUBROUTINE CC_QM3RAINT(RAAOx,RAAOy,RAAOz,AODEN,
     &                       RAx,RAy,RAz,WORK,LWORK,
     &                       ISYMT)
C**************************************************************
C
C     Readin mm-gradient integrals in coupled cluster format
C     and evaluate the expectation value of Rra. Calculates
C     the induced dipole moment at the center of mass a and
C     evaluates the O(mmss) factor at the center of mass a
C
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "qm3.h"
#include "ccslvinf.h"
#include "ccinftap.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
C
      DIMENSION WORK(LWORK)
      DIMENSION RAAOx(*),RAAOy(*),RAAOz(*)
      DIMENSION AODEN(*)
      DIMENSION RAx(*)
      DIMENSION RAy(*)
      DIMENSION RAz(*)
      DIMENSION EELEC(3,MXQM3)
      DIMENSION FFJ(3)
      CHARACTER*8 LABEL
      INTEGER   KL
      LOGICAL   LOINDM
C
C--------------------------------------------------------
C     Initializing the electronic electric field to zero:
C--------------------------------------------------------
C
      DO 879 I = 1, MXQM3
        DO 880 J = 1, 3
          EELEC(J,I) = 0.0D0 
  880   CONTINUE
  879 CONTINUE
C
      IF (IQM3PR .GT. 10) THEN
        WRITE(LUPRI,*) 'CC_QM3RAINT: Read in integrals'
        WRITE(LUPRI,*) 'Input symmetry claimed', isymt 
      ENDIF
C
      LUCCEF = -1
      CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCEF)
C
      KL = 0
C
      DO 950 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 951 J = NSYSBG(I), NSYSED(I)
            DO 952 K = 1, NUALIS(I)
C
              KL = KL +1
C
              LABEL = 'READ INT'
              CALL CCMMINT(LABEL,LUCCEF,RAAOx,DUMMY,IRREP,
     *                   ISYM,IERR,WORK,LWORK)
              RAx(KL) = DDOT(N2BST(ISYMOP),RAAOx,1,AODEN,1)
C
              LABEL = 'READ INT'
              CALL CCMMINT(LABEL,LUCCEF,RAAOy,DUMMY,IRREP,
     *                   ISYM,IERR,WORK,LWORK)
              RAy(KL) = DDOT(N2BST(ISYMOP),RAAOy,1,AODEN,1)
C
              LABEL = 'READ INT'
              CALL CCMMINT(LABEL,LUCCEF,RAAOz,DUMMY,IRREP,
     *                   ISYM,IERR,WORK,LWORK)
              RAz(KL) = DDOT(N2BST(ISYMOP),RAAOz,1,AODEN,1)
C
              IF (IQM3PR .GT. 10) THEN 
                WRITE(LUPRI,*)'RAx(KL) =',RAx(KL)
                WRITE(LUPRI,*)'RAy(KL) =',RAy(KL)
                WRITE(LUPRI,*)'RAz(KL) =',RAz(KL)
              END IF
C
  952       CONTINUE
  951     CONTINUE
        END IF 
  950 CONTINUE
C
      IF (KL .NE. NCOMS) THEN
        CALL QUIT('Error in no. of center of masses in CC_QM3RAINT')
      END IF
C
      DO 534 LM = 1,NCOMS
        EELEC(1,LM) = RAx(LM)
        EELEC(2,LM) = RAy(LM)
        EELEC(3,LM) = RAz(LM)
  534 CONTINUE
C
      CALL GPCLOSE(LUCCEF,'KEEP')
C
C     If RELMOM is true we want to include the external field(s) in
C     the determination of the induced dipole moments
C
        FFJ(1) = 0.0D0
        FFJ(2) = 0.0D0
        FFJ(3) = 0.0D0
C
      IF (RELMOM) THEN
        DO 330 IF =1, NFIELD
          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'ZDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
  330   CONTINUE
      END IF 
C
      IF (FIXMOM) THEN
        WRITE(LUPRI,'(/A)')'FIXMOM: NO ITER. DETERM. OF MYIND'
      ELSE IF (LOSPC) THEN
        WRITE(LUPRI,'(/A)')'All MM models are SPC, INDMOM not called'
      ELSE 
        LOINDM = .FALSE.
        CALL INDMOM(EELEC,LOINDM,FFJ)
      END IF
C
      CALL QM3_OBAR(FFJ)
      CALL CC_PUT31('CC_RA',NCOMS,RAx,RAy,RAz)
C
      END
C
C**************************************************************
C  /* Deck cc_qm3nsint */
      SUBROUTINE CC_QM3NSINT(NSAO,AODEN,CCNS,WORK,LWORK,ISYMT)
C**************************************************************
C
C     Readin mm-potentiale energy integrals in coupled cluster
C     format and calculate the expectation value of Ns.
C
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccinftap.h"
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
C
      DIMENSION WORK(LWORK),NSAO(*)
      DIMENSION AODEN(*)
      DIMENSION CCNS(*)
C
      CHARACTER*8 LABEL
C
      LUCCPO = -1
      CALL GPOPEN(LUCCPO,'POTMM','UNKNOWN',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCPO)
C
      LM = 0
C
      DO 940 I= 1, ISYTP
        IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
          DO 941 J = NSYSBG(I), NSYSED(I)
            DO 942 K = 1, NSISY(I)
              LM = LM + 1
              LABEL = 'READ INT'
              CALL CCMMINT(LABEL,LUCCPO,NSAO,DUMMY,IRREP,
     &                     ISYM,IERR,WORK,LWORK)
              CCNS(LM) = DDOT(N2BST(ISYMOP),NSAO,1,AODEN,1)
  942       CONTINUE         
  941     CONTINUE
        END IF
  940 CONTINUE
C
      CALL GPCLOSE(LUCCPO,'KEEP')
C
C-------------------
C     Print section:
C-------------------
C
      IF (IQM3PR .GT. 5) THEN
        WRITE(LUPRI,'(//,A)')
     *        ' +======================+'
        WRITE(LUPRI,'(A)')
     *        ' |Site| <N_s>           |' 
        WRITE(LUPRI,'(1X,A)')
     *        '+======================+'
C
        LS = 0
C
        DO 215 I = 1, ISYTP
          IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO 216 J = NSYSBG(I), NSYSED(I)
              DO 217 K = 1, NSISY(I)
C
                LS = LS + 1
C
                WRITE(LUPRI,'(1X,A,I3,A,F16.10,A)')
     *                      '| ', LS,'|', CCNS(LS),' |'
                WRITE(LUPRI,'(1X,A)')
     *                      '+----------------------+'
  217         CONTINUE
  216       CONTINUE
          END IF
  215   CONTINUE
        WRITE(LUPRI,'(//,A)')
      END IF
C
      END
C
C**************************************************************
C  /* Deck ccmm_rhstg */
      SUBROUTINE CCMM_RHSTG(FOCK,WORK,LWORK)
C**************************************************************
C
C     Direct calculation of CC/MM solvent effects. If
C     OLDTG = .TRUE. the Ns part is included in the T^g
C     operator. However, if OLDTG = .FALSE. the contrinution is 
C     allready in the h1 "vacuum" operator.  
C
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "qm3.h"
#include "ccinftap.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
C
      DIMENSION WORK(LWORK),FOCK(*)
      DIMENSION ALPHAMM(MXQM3)
C
      CHARACTER*8 LABEL
      CHARACTER*9 FILMMR
C
      LOGICAL FIRST
      SAVE  FIRST
C
C------------------------------------------------------------------
C     If a system is model SPC_ECX (X=1,2,3,4) the polarization of
C     the system does not contribute to the optimization of the 
C     wave function. Also, if OLDTG = .FALSE. the partial charges 
C     are included in h1. So if all models are SPC -> RETURN from 
C     CCMM_RHSTG.
C------------------------------------------------------------------
C
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) RETURN
C
C----------------------------
C       Dynamical allocation:
C----------------------------
C
      IF ( (.NOT. (OLDTG)) .AND. (.NOT. (LOSPC)) ) THEN
C
C--------------------------
C       No space needed for 
C       Ns in AO basis:
C--------------------------
C
        KTAO  =  1 
        KRAAOx = KTAO   + NNBASX
        KRAAOy = KRAAOx + NNBASX
        KRAAOz = KRAAOy + NNBASX
        KENSAx = KRAAOz + NNBASX
        KENSAy = KENSAx + NCOMS
        KENSAz = KENSAy + NCOMS
        KRAx = KENSAz + NCOMS 
        KRAy = KRAx + NCOMS
        KRAz = KRAy + NCOMS
        KWRK1 = KRAz + NCOMS 
        LWRK1   = LWORK   - KWRK1
C
        CALL DZERO(WORK(KTAO),4*NNBASX)
        CALL DZERO(WORK(KENSAx),3*NCOMS)
        CALL DZERO(WORK(KRAx),3*NCOMS)
C
      ELSE IF ( (OLDTG) .AND. (.NOT. (LOSPC)) ) THEN
C
C--------------------------
C       Space needed for 
C       Ns in AO basis:
C--------------------------
C
        KTAO  =  1
        KRAAOx = KTAO   + NNBASX
        KRAAOy = KRAAOx + NNBASX
        KRAAOz = KRAAOy + NNBASX
        KENSAx = KRAAOz + NNBASX
        KENSAy = KENSAx + NCOMS
        KENSAz = KENSAy + NCOMS
        KNSAO = KENSAz + NCOMS
        KRAx = KNSAO + NNBASX 
        KRAy = KRAx + NCOMS
        KRAz = KRAy + NCOMS
        KWRK1 = KRAz + NCOMS
        LWRK1   = LWORK   - KWRK1
C
        CALL DZERO(WORK(KTAO),4*NNBASX)
        CALL DZERO(WORK(KENSAx),3*NCOMS)
        CALL DZERO(WORK(KNSAO),NNBASX)
        CALL DZERO(WORK(KRAx),3*NCOMS)
C
      ELSE IF ( (OLDTG) .AND. (LOSPC) ) THEN
C
C--------------------------
C       No space needed for 
C       alpha contributers:
C--------------------------
C
        KTAO  = 1
        KNSAO = KTAO + NNBASX
        KWRK1 = KNSAO + NNBASX
        LWRK1   = LWORK   - KWRK1
C
        CALL DZERO(WORK(KTAO),2*NNBASX)
C
      END IF 
C
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_RHSTG' )
C
C-------------------------------------------------
C     If OLDTG and LOSPC are .TRUE. -> GOTO 888
C     which is after the alpha part in CCMM_RHSTG.
C-------------------------------------------------
C
      IF ( (OLDTG) .AND. (LOSPC) ) GOTO 888 
C
      CALL CC_GET31('CC_RA',NCOMS,
     *               WORK(KRAx),WORK(KRAy),WORK(KRAz))
C
C-------------------
C     Print section:
C-------------------
C
      IF (IQM3PR .GT. 10) THEN
        WRITE(LUPRI,'(//,A)')
     *        ' +============================================'
     *      //'==============+'
        WRITE(LUPRI,'(A)')
     *        ' | COM| <Rra>_x         | <Rra>_y         |'
     *      //' <Rra>_z         |'
        WRITE(LUPRI,'(1X,A)')
     *        '+============================================'
     *      //'==============+'
C
         NUM = 0
C
         DO 215 I = 1, ISYTP
           IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
             DO 216 J = NSYSBG(I), NSYSED(I)
               DO 217 K=1,NUALIS(I)
C
                 NUM = NUM + 1
C
                 WRITE(LUPRI,'(1X,A,I3,A,F16.10,A,
     &                         F16.10,A,F16.10,A)')
     *           '| ', J,'|', WORK(KRAx + NUM - 1),' |',
     *           WORK(KRAy + NUM - 1),' |', 
     *           WORK(KRAz + NUM - 1),' |'
                 WRITE(LUPRI,'(1X,A)')
     *           '+---------------------------------------------'
     *           //'-------------+'
  217          CONTINUE 
  216        CONTINUE
           END IF
  215    CONTINUE 
         WRITE(LUPRI,'(//,A)')
      END IF
C
      CALL CC_GET31('ENSAFILE',NCOMS,
     *               WORK(KENSAx),WORK(KENSAy),WORK(KENSAz))

C
      IF (IQM3PR .GT. 10) THEN
        WRITE(LUPRI,'(//,A)')
     *        ' +============================================'
     *      //'==============+'
        WRITE(LUPRI,'(A)')
     *        ' | COM| ENSA_x          | ENSA_y          |'
     *      //' ENSA_z          |'
        WRITE(LUPRI,'(1X,A)')
     *        '+============================================'
     *      //'==============+'
C
         NUM = 0
C
         DO 415 I = 1, ISYTP
           IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
C
             NUM = NUM + 1
C
             DO 416 J = NSYSBG(I), NSYSED(I)
               DO 417 K=1,NUALIS(I)
                 WRITE(LUPRI,'(1X,A,I3,A,F16.10,A,
     &                         F16.10,A,F16.10,A)')
     *           '| ', J,'|', WORK(KENSAx + NUM - 1),' |',
     *           WORK(KENSAy + NUM - 1),' |', 
     *           WORK(KENSAz + NUM - 1),' |'
                 WRITE(LUPRI,'(1X,A)')
     *           '+---------------------------------------------'
     *           //'-------------+'
  417          CONTINUE
  416        CONTINUE
           END IF
  415    CONTINUE
         WRITE(LUPRI,'(//,A)')
      END IF
C
C---------------------------------------
C Read Rr(s) in ao basis and add to tao:
C---------------------------------------
C
      LUCCEF = -1
      CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCEF)
C
      LM = 0
C
      DO 520 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 521 J = NSYSBG(I), NSYSED(I)
            DO 522 K = 1, NUALIS(I)
              LM = LM + 1
C
              CALL READT(LUCCEF,NNBASX,WORK(KRAAOx))
C
              IF (IQM3PR .GE. 15) THEN
                WRITE (LUPRI,'(/A)') ' Rra_ao x matrix:'
                CALL OUTPAK(WORK(KRAAOx),NBAST,1,LUPRI)
              END IF
C
              CALL READT(LUCCEF,NNBASX,WORK(KRAAOy))
C
              IF (IQM3PR .GE. 15) THEN
                WRITE (LUPRI,'(/A)') ' Rra_ao y matrix:'
                CALL OUTPAK(WORK(KRAAOy),NBAST,1,LUPRI)
              END IF
C
              CALL READT(LUCCEF,NNBASX,WORK(KRAAOz))
C
              IF (IQM3PR .GE. 15) THEN
                WRITE (LUPRI,'(/A)') ' Rra_ao z matrix:'
                CALL OUTPAK(WORK(KRAAOz),NBAST,1,LUPRI)
              END IF
C
              IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
                FACx =  - ALPIMM(I,K) * 
     &                  (WORK(KRAx + LM - 1)  
     &                + 0.5D0 * WORK(KENSAx + LM - 1)) 
                FACy =  - ALPIMM(I,K) * 
     &                  (WORK(KRAy + LM - 1)  
     &                + 0.5D0 * WORK(KENSAy + LM - 1))
                FACz =  - ALPIMM(I,K) * 
     &                  (WORK(KRAz + LM - 1)  
     &                + 0.5D0 * WORK(KENSAz + LM - 1))
C
                CALL DAXPY(NNBASX,FACx,WORK(KRAAOx),
     *                   1,WORK(KTAO),1)
C
                CALL DAXPY(NNBASX,FACy,WORK(KRAAOy),
     *                   1,WORK(KTAO),1)
C
                CALL DAXPY(NNBASX,FACz,WORK(KRAAOz),
     *                   1,WORK(KTAO),1)
              END IF
  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE
C
      IF (IQM3PR.GT.14) THEN
        WRITE (LUPRI,*) ' NORM of TAO matrix /alpha contr.: ',
     *  DDOT(NNBASX,WORK(KTAO),1,WORK(KTAO),1)
        WRITE (LUPRI,'(/A)') ' TAO matrix: '
        CALL OUTPAK(WORK(KTAO),NBAST,1,LUPRI)
      ENDIF
C
      CALL GPCLOSE(LUCCEF,'KEEP')
C 
C--------------------------------------
C     End of alpha part in CCMM_RHSTG:
C--------------------------------------
C
  888 CONTINUE
C
C-----------------------------------------------------------
C     If OLDTG we want to include the partial MM charges in
C     the T^g operator in the opt. of the wave function.
C-----------------------------------------------------------
C
      IF (OLDTG) THEN
C
        LUCCPO = -1
        CALL GPOPEN(LUCCPO,'POTMM','UNKNOWN',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND (LUCCPO)
C
        FAC1 = -1.0D0
C
        L = 0
C
        DO 525 I = 1, ISYTP
          IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO 526 J = NSYSBG(I), NSYSED(I)
              DO 527 K = 1,NSISY(I)
C
                L = L +1
C
                CALL READT(LUCCPO,NNBASX,WORK(KNSAO))
C
                IF (IQM3PR .GE. 4) THEN
                  WRITE (LUPRI,'(/A,I3,A)') 
     *                ' N(',L,')_ao matrix: '
                  CALL OUTPAK(WORK(KNSAO),NBAST,1,LUPRI)
                END IF
C
                CALL DAXPY(NNBASX,FAC1,WORK(KNSAO),
     *                     1,WORK(KTAO),1)
  527         CONTINUE
  526       CONTINUE
          END IF
  525   CONTINUE
C
        IF (IQM3PR.GT.14) THEN
           WRITE (LUPRI,*) ' NORM of TAO matrix /nsao: ',
     *     DDOT(NNBASX,WORK(KTAO),1,WORK(KTAO),1)
           WRITE (LUPRI,'(/A)') ' TAO matrix: '
           CALL OUTPAK(WORK(KTAO),NBAST,1,LUPRI)
        ENDIF
C
        CALL GPCLOSE(LUCCPO,'KEEP')
C
      END IF
C
C--------------------------------------------------
C     Add contribution to effective AO fock matrix:
C--------------------------------------------------
C
      FACT = 1.0D0
C
C--------------------------
C     Dynamical allocation:
C--------------------------
C
      KINT = KWRK1
      KWRK2 = KINT + N2BST(ISYMOP)
      LWRK2 = LWORK - KWRK2 
C
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_RHSTG' )
C
      CALL DZERO(WORK(KINT),N2BST(ISYMOP))
C
      IF (IQM3PR .GT. 14) THEN
        WRITE (LUPRI,*) 'NORM of TAO matrix. CCMMINT: ',
     *  DDOT(NNBASX,WORK(KTAO),1,WORK(KTAO),1)
        WRITE (LUPRI,'(/A)') ' TAO matrix: '
        CALL OUTPAK(WORK(KTAO),NBAST,1,LUPRI)
      ENDIF
C
      LABEL= 'GIVE INT'
C
      CALL CCMMINT(LABEL,LUCCEF,WORK(KINT),WORK(KTAO),
     &             IRREP,ISYM,IERR,WORK(KWRK2),LWRK2)
C
      IF (IQM3PR .GT.14) THEN
        CALL AROUND( 'CCMM cont. to AO matrix' )
        CALL CC_PRFCKAO(WORK(KINT),ISYMOP)
      ENDIF
C
      IF (IQM3PR .GT.14) THEN
        CALL AROUND( 'Usual Fock AO matrix' )
        CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
      CALL DAXPY(N2BST(ISYMOP),FACT,WORK(KINT),1,FOCK,1)
C
      IF (IQM3PR .GT.14) THEN
        CALL AROUND( 'Modified Fock AO matrix before rep.')
        CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
      IF (LOSHAW) THEN
        CALL DZERO(WORK(KINT),N2BST(ISYMOP))
C
        FILMMR  = 'MMREPOP_0'
        LUMMRE  = -1
        IADRFIL = 1
        NDATA   = N2BST(ISYMOP)
C
        CALL WOPEN2(LUMMRE,FILMMR,64,0)
C
        CALL GETWA2(LUMMRE,FILMMR,WORK(KINT),IADRFIL,NDATA)
C
        IF (REPTST) THEN
          NBASLO = MAX(NBAST,1)
C
          KINT1 = KWRK2
          KINT2 = KINT1 + N2BST(ISYMOP)
          KWRK3 = KINT2 + N2BST(ISYMOP)
          LWRK3 = LWORK - KWRK3
C
          IF (LWRK3 .LT. 0) THEN
            CALL QUIT( 'Too litle work in CCMM_RHSTG Rep. 1')
          END IF
C
          CALL DZERO(WORK(KINT1),2*N2BST(ISYMOP))
          CALL DCOPY(N2BST(ISYMOP),WORK(KINT),1,WORK(KINT1),1)
C
          DO 668 KR=1,NREPMT
            CALL DGEMM('N','N',NBAST,NBAST,NBAST,
     *                 ONE,WORK(KINT),NBASLO,
     *                 WORK(KINT1),NBASLO,
     *                 ZERO,WORK(KINT2),NBASLO)
C
            CALL DCOPY(N2BST(ISYMOP),WORK(KINT2),1,WORK(KINT),1) 
  668     CONTINUE
        END IF        
C
          TAL1 = DDOT(N2BST(ISYMOP),WORK(KINT),1,WORK(KINT),1)
          TAL2 = SQRT(TAL1)
C
        CALL DAXPY(N2BST(ISYMOP),FACT,WORK(KINT),1,FOCK,1)
C
        IF (IQM3PR .GT.14) THEN
          CALL AROUND( 'Modified Fock AO matrix after rep.')
          CALL CC_PRFCKAO(FOCK,ISYMOP)
        ENDIF
C
        CALL WCLOSE2(LUMMRE,FILMMR, 'KEEP')
C
      END IF !LOSHAW
C
      END
C
C**************************************************************
C  /* Deck ccmm_ltrb */
      SUBROUTINE CCMM_LTRB(RHO1,RHO2,CTR1,CTR2,
     &                     ISYMTR,LR,WORK,LWORK)
C**************************************************************
C
C     Calculation of CCMM T^gB contribution to left and right
C     Jacobian transformation:
C
C     <mu|exp(-T)T^g|CC> for LR = 0
C     F transformation for LR = F
C     P transformation for LR = P
C
C     LR = 'L','R','0','F','P'
C
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "qm3.h"
#include "ccinftap.h"
#include "ccslvinf.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*),CTR1(*),CTR2(*) 
      DIMENSION ALPHAMM(MXQM3)
C
      CHARACTER*8 LABEL,LABEL1,LIST*(2),LR*(1)
      CHARACTER*(8) FILMME, FILMMX
      CHARACTER*9 FILMMR
      LOGICAL LEXIST
C
C----------------------------------------------
C     If all models are SPC and OLDTG = .FALSE.
C     -> RETURN from CCMM_LTRB:
C----------------------------------------------
C
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) RETURN
C
C-----------------------------------------------------
C     Also, if all models are SPC and OLDTG are .TRUE. 
C     only enter CCMM_LTRB if LR .NE.'0':
C-----------------------------------------------------
C
      IF ( (LOSPC) .AND. (OLDTG) .AND. (LR .NE.'0') ) RETURN
C
      IF (IQM3PR .GT. 10) THEN
        WRITE(LUPRI,*)'CCMM_LTRB: CC/MM contribution to CC L. transf.'
        WRITE(LUPRI,*)'CCMM_LTRB: LWORK:', LWORK
        WRITE(LUPRI,*)'CCMM_LTRB: LR:', LR    
      ENDIF
C
      IF ( DEBUG .OR.(IQM3PR .GT. 10)) THEN
        RHO1N = DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
        RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
        WRITE(LUPRI,*) ' Norm of RHO1 in CCMM_LTGB on input:', RHO1N
        WRITE(LUPRI,*) ' Norm of RHO2 in CCMM_LTGB on input:', RHO2N
        RHO1N = DDOT(NT1AM(ISYMTR),CTR1,1,CTR1,1)
        RHO2N = DDOT(NT2AM(ISYMTR),CTR2,1,CTR2,1)
        WRITE(LUPRI,*) ' Norm af C1AM in CCMM_LTGB on input:', RHO1N
        WRITE(LUPRI,*) ' Norm af C2AM in CCMM_LTGB on input:', RHO2N
      ENDIF
C
C---------------------
C     Init parameters:
C---------------------
C
      NTAMP1  = NT1AM(ISYMTR)
      NTAMP2  = NT2AM(ISYMTR)
      NTAMP   = NTAMP1 + NTAMP2
C
      IF (CCS)  NT2AM(ISYMTR) = IZERO
C------------------------------
C
      IF (DISCEX) THEN
        LUMMET = -1
        LUMMXI = -1
        FILMME = 'CCMM_ETA'
        FILMMX = 'CCMM__XI'
        CALL WOPEN2(LUMMET, FILMME, 64, 0)
        CALL WOPEN2(LUMMXI, FILMMX, 64, 0)
      END IF
C
C-----------------------------------------
C
      KTGB = 1
      KWRK1 = KTGB + N2BST(ISYMTR)
      LWRK1  = LWORK - KWRK1
C
      IF (LWRK1 .LT. 0) CALL QUIT( ' Too litle work in CCMM_LTRB 1' )
C
      CALL DZERO(WORK(KTGB),N2BST(ISYMTR))
C
C--------------------------------------------------------
C     If OLDTG, LOSPC = .TRUE. and (LR .EQ.'0') only the
C     contribution due to the partial charges in the T^g
C     operator -> GOTO 888 (after the alpha part).
C--------------------------------------------------------
C
      IF ( (LOSPC) .AND. (OLDTG) .AND. (LR .EQ.'0') ) GOTO 888
C
C--------------------------------------------------------
C  Calculate T^{gB}= - alpha sum_a <B|Rra|CC>Rra
C  and contributions from this term. For LR = 0
C  calculates <mu|T^g|CC>.
C--------------------------------------------------------
C
C--------------------------
C     Dynamical allocation:
C--------------------------
C
      KRAAOx = KWRK1
      KRAAOy = KRAAOx + N2BST(ISYMTR)
      KRAAOz = KRAAOy + N2BST(ISYMTR)
      KENSAx = KRAAOz + N2BST(ISYMTR)
      KENSAy = KENSAx + NCOMS
      KENSAz = KENSAy + NCOMS
      KRAx   = KENSAz + NCOMS
      KRAy   = KRAx + NCOMS
      KRAz   = KRAy + NCOMS
      KWRK2  = KRAz + NCOMS
      LWRK2  = LWORK   - KWRK2
C
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_LTRB, 2')
C
      KXIx  = KWRK2
      KXIy  = KXIx + NTAMP
      KXIz  = KXIy + NTAMP
      KWRK3 = KXIz  + NTAMP
      LWRK3   = LWORK   - KWRK3
C
      IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CCMM_LTRB, 3')
C
      CALL DZERO(WORK(KRAAOx),3*N2BST(ISYMTR))
      CALL DZERO(WORK(KENSAx),3*NCOMS)
      CALL DZERO(WORK(KRAx),3*NCOMS)
      CALL DZERO(WORK(KXIx),3*NTAMP)
C 
      CALL CC_GET31('ENSAFILE',NCOMS,
     *               WORK(KENSAx),WORK(KENSAy),WORK(KENSAz))
      CALL CC_GET31('CC_RA',NCOMS,
     *               WORK(KRAx),WORK(KRAy),WORK(KRAz))
C
      LUCCEF = -1
      CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCEF)
C
C----------------------------
C     Readin integrals again,
C     Read Rra in AO basis:
C----------------------------
C
      LM = 0
      IADR1 = 1
      IADR2 = 1
      LEN = NTAMP
C
      DO 700 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 701 J = NSYSBG(I), NSYSED(I)
            DO 702 K = 1, NUALIS(I)
C
              LM = LM + 1
              LABEL = 'READ INT'
C
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOx),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOy),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOz),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
C
C------------------------------------------------------------
C     Calculate xi^Rra  (for LR=R actually eta^Rra )
C
C     Only if MDLWRD(I) = 'SPC_E01' are we going to calculate
C     contributions from the T^g operator to the response 
C     equations.
C------------------------------------------------------------
C
              IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
                IF (.NOT. (DISCEX)) THEN
C
                  LABEL = 'GIVE INT'           
                  IF ( (LR.EQ.'L').OR.(LR.EQ.'0')
     &              .OR.(LR.EQ.'P') ) THEN
                    CALL CC_XKSI(WORK(KXIx),
     *                   LABEL,ISYMTR,0,WORK(KRAAOx),
     *                   WORK(KWRK3),LWRK3)
                  ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
                    LIST  = 'L0'
                    ILSTNR  = 1
                    CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIx),
     *                            LIST,ILSTNR,0,WORK(KRAAOx),
     *                            WORK(KWRK3),LWRK3)
                  END IF
C
                  IF ( (LR.EQ.'L').OR.(LR.EQ.'0')
     &              .OR.(LR.EQ.'P') ) THEN
                    CALL CC_XKSI(WORK(KXIy),LABEL,ISYMTR,
     *                           0,WORK(KRAAOy),
     *                           WORK(KWRK3),LWRK3)
                  ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
                    LIST  = 'L0'
                    ILSTNR  = 1
                    CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIy),
     *                           LIST,ILSTNR,0,WORK(KRAAOy),
     *                           WORK(KWRK3),LWRK3)
                  END IF
C
                  IF ( (LR.EQ.'L').OR.(LR.EQ.'0')
     &              .OR.(LR.EQ.'P') ) THEN
                    CALL CC_XKSI(WORK(KXIz),LABEL,ISYMTR,
     *                           0,WORK(KRAAOz),
     *                   WORK(KWRK3),LWRK3)
                  ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
                    LIST  = 'L0'
                    ILSTNR  = 1
                    CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIz),
     *              LIST,ILSTNR,0,WORK(KRAAOz),WORK(KWRK3),LWRK3)
                  END IF
C
                ELSE
C
                  IF ( (LR.EQ.'L').OR.(LR.EQ.'P') ) THEN
                    CALL GETWA2(LUMMXI,FILMMX,WORK(KXIx),IADR2,LEN)
                    IADR2 = IADR2 + LEN
                  ELSE IF ( (LR.EQ.'R').OR.(LR.EQ.'F') ) THEN
                    CALL GETWA2(LUMMET,FILMME,WORK(KXIx),IADR1,LEN)
                    IADR1 = IADR1 + LEN
                  ELSE IF (LR.EQ.'0') THEN
                    LABEL = 'GIVE INT' 
                    CALL CC_XKSI(WORK(KXIx),LABEL,ISYMTR,0,WORK(KRAAOx),
     *                           WORK(KWRK3),LWRK3)
                  END IF
C
                  IF ( (LR.EQ.'L').OR.(LR.EQ.'P') ) THEN
                    CALL GETWA2(LUMMXI,FILMMX,WORK(KXIy),IADR2,LEN)
                    IADR2 = IADR2 + LEN
                  ELSE IF ( (LR.EQ.'R').OR.(LR.EQ.'F') ) THEN
                    CALL GETWA2(LUMMET,FILMME,WORK(KXIy),IADR1,LEN)
                    IADR1 = IADR1 + LEN
                  ELSE IF (LR.EQ.'0') THEN
                    LABEL = 'GIVE INT'
                    CALL CC_XKSI(WORK(KXIy),LABEL,ISYMTR,0,WORK(KRAAOy),
     *                           WORK(KWRK3),LWRK3)
                  END IF
C               
                  IF ( (LR.EQ.'L').OR.(LR.EQ.'P') ) THEN
                    CALL GETWA2(LUMMXI,FILMMX,WORK(KXIz),IADR2,LEN)
                    IADR2 = IADR2 + LEN
                  ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
                    CALL GETWA2(LUMMET,FILMME,WORK(KXIz),IADR1,LEN)
                    IADR1 = IADR1 + LEN
                  ELSE IF (LR.EQ.'0') THEN
                    LABEL = 'GIVE INT'
                    CALL CC_XKSI(WORK(KXIz),LABEL,ISYMTR,0,WORK(KRAAOz),
     *                           WORK(KWRK3),LWRK3)
                  END IF
                END IF
C
C-------------------------------
C               Contract with B:
C-------------------------------
C
                IF (LR.NE.'0') THEN
C
                  KXI1x = KXIx  
                  KXI2x = KXIx + NTAMP1
                  BXILMD1x= DDOT(NTAMP1,CTR1,1,WORK(KXI1x),1)
                  BXILMD2x= DDOT(NTAMP2,CTR2,1,WORK(KXI2x),1)
                  BXILMDx = BXILMD1x + BXILMD2x
C
                  KXI1y = KXIy 
                  KXI2y = KXIy + NTAMP1
                  BXILMD1y = DDOT(NTAMP1,CTR1,1,WORK(KXI1y),1)
                  BXILMD2y = DDOT(NTAMP2,CTR2,1,WORK(KXI2y),1)
                  BXILMDy = BXILMD1y + BXILMD2y
C
                  KXI1z = KXIz 
                  KXI2z = KXIz  + NTAMP1
                  BXILMD1z = DDOT(NTAMP1,CTR1,1,WORK(KXI1z),1)
                  BXILMD2z = DDOT(NTAMP2,CTR2,1,WORK(KXI2z),1)
                  BXILMDz  = BXILMD1z + BXILMD2z
C

                  FACx =  - ALPIMM(I,K) * (BXILMDx)
                  FACy =  - ALPIMM(I,K) * (BXILMDy)
                  FACz =  - ALPIMM(I,K) * (BXILMDz)
C
                  IF (IQM3PR .GT. 5) THEN
                    WRITE(LUPRI,*)'---------------------------' //
     *                            '-------------------'
                    WRITE(LUPRI,*)'The factors calculated at cent.' //
     *                            ' of mass no.',LM
                    WRITE(LUPRI,*)'---------------------------' //
     *                            '-------------------'
                    WRITE(LUPRI,*)'FACx=',FACx
                    WRITE(LUPRI,*)'FACy=',FACy
                    WRITE(LUPRI,*)'FACz=',FACz
                    WRITE(LUPRI,*)'---------------------------' //
     *                            '-------------------'
                  END IF
C
C------------------------------------------
C                 Add to T^{qB} integrals.
C                 (for LR=R actually T^gC):
C------------------------------------------
C
                  CALL DAXPY(N2BST(ISYMTR),FACx,WORK(KRAAOx),1,
     *                       WORK(KTGB),1)
C
                  CALL DAXPY(N2BST(ISYMTR),FACy,WORK(KRAAOy),1,
     *                       WORK(KTGB),1)
C
                  CALL DAXPY(N2BST(ISYMTR),FACz,WORK(KRAAOz),1,
     *                       WORK(KTGB),1)
C
                END IF
                IF (LR.EQ.'0') THEN
C
                  KXI1x = KXIx  
                  KXI2x = KXIx  + NTAMP1
C
                  FACx =  - ALPIMM(I,K) * (WORK(KRAx + LM - 1) 
     &                    + 0.5D0 * WORK(KENSAx + LM - 1))
C
                  CALL DAXPY(NT1AM(ISYMTR),FACx,WORK(KXI1x),1,
     *                       RHO1,1)
                  CALL DAXPY(NT2AM(ISYMTR),FACx,WORK(KXI2x),1,
     *                       RHO2,1)
C
                  KXI1y= KXIy
                  KXI2y = KXIy  + NTAMP1
C
                  FACy =  - ALPIMM(I,K) * (WORK(KRAy + LM - 1) 
     &                    + 0.5D0 * WORK(KENSAy + LM - 1))
C
                  CALL DAXPY(NT1AM(ISYMTR),FACy,WORK(KXI1y),1,
     *                       RHO1,1)
                  CALL DAXPY(NT2AM(ISYMTR),FACy,WORK(KXI2y),1,
     *                       RHO2,1)
C
                  KXI1z = KXIz
                  KXI2Z = KXIz  + NTAMP1
C
                  FACz =  - ALPIMM(I,K) * (WORK(KRAz + LM - 1) 
     &                    + 0.5D0 * WORK(KENSAz + LM - 1))
C
                  CALL DAXPY(NT1AM(ISYMTR),FACz,WORK(KXI1z),1,
     *                       RHO1,1)
                  CALL DAXPY(NT2AM(ISYMTR),FACz,WORK(KXI2z),1,
     *                       RHO2,1)
C
                END IF
              END IF
  702       CONTINUE
  701     CONTINUE
        END IF
  700 CONTINUE
C
      CALL GPCLOSE(LUCCEF,'KEEP')
C
      IF (DISCEX) THEN
        CALL WCLOSE2(LUMMET, FILMME, 'KEEP')
        CALL WCLOSE2(LUMMXI, FILMMX, 'KEEP')
      END IF
C
C-----------------------
C     End of alpha part:
C-----------------------
C
  888 CONTINUE
C
C------------------------------
C     If OLDTG = .TRUE. 
C     N_s part -> T^g operator:
C------------------------------
C  
      IF (OLDTG) THEN
        IF (LR.EQ.'0') THEN
C
          KNSAO = KWRK1
          KWRK2 = KNSAO + N2BST(ISYMTR)
          LWRK2 = LWORK - KWRK2
C
          CALL DZERO(WORK(KNSAO),N2BST(ISYMTR))
C
          IF (LWRK2 .LT. 0) 
     &        CALL QUIT( ' Too litle work in CCMM_LTRB 4b' )
C
          KXI   = KWRK2
          KWRK3 = KXI + NTAMP
          LWRK3 = LWORK - KWRK3
C
          CALL DZERO(WORK(KXI),NTAMP)
C
          IF (LWRK3.LT.0) 
     &        CALL QUIT( 'Too little work in CCMM_LTRB, 5')
C
          LUCCPO = -1
          CALL GPOPEN(LUCCPO,'POTMM','OLD',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
          REWIND (LUCCPO)
C
          FAC1 = -1.0D0
          LM   =  0
C
          DO 800 I = 1, ISYTP
            IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
              DO 801 J = NSYSBG(I), NSYSED(I)
                DO 802 K = 1,NSISY(I)
C
                  LM = LM +1
                  LABEL = 'READ INT'
C
                  CALL CCMMINT(LABEL,LUCCPO,WORK(KNSAO),DUMMY,IRREP,
     &                         ISYM,IERR,WORK(KWRK3),LWRK3)
C
C---------------------------------
C                 Calculate xi^Ns:
C---------------------------------
C
                  LABEL = 'GIVE INT'
C
                  CALL CC_XKSI(WORK(KXI),LABEL,ISYMTR,0,
     *                         WORK(KNSAO),WORK(KWRK3),LWRK3)
C
                  KXI1 = KXI
                  KXI2 = KXI  + NTAMP1
C
                  CALL DAXPY(NT1AM(ISYMTR),FAC1,WORK(KXI1),1,
     *                       RHO1,1)
                  CALL DAXPY(NT2AM(ISYMTR),FAC1,WORK(KXI2),1,
     *                       RHO2,1)
  802           CONTINUE
  801         CONTINUE
            END IF
  800     CONTINUE
        IF ( DEBUG .OR.(IQM3PR .GT. 10)) THEN
          RHO1N = DDOT(NT1AM(ISYMTR),WORK(KXI1),1,WORK(KXI1),1)
          RHO2N = DDOT(NT2AM(ISYMTR),WORK(KXI2),1,WORK(KXI2),1)
          WRITE(LUPRI,*) ' Norm af elstat contribution to LHTR1:', RHO1N
          WRITE(LUPRI,*) ' Norm af elstat contribution to LHTR2:', RHO2N
        END IF
C
          CALL GPCLOSE(LUCCPO,'KEEP')
C
        END IF
      END IF
C
C-------------------------------------------------------
C
C     Introduce repulsion component into T^g
C
C-------------------------------------------------------
C
      IF (LOSHAW) THEN 
        IF (LR.EQ.'0') THEN
C
        KINT  = KWRK1
        KWRK2 = KINT  + N2BST(ISYMOP)
        LWRK2 = LWORK - KWRK2
C
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT( 'Too litle work in CCMM_LTRB Rep. 1')
        END IF
C
        CALL DZERO(WORK(KINT),N2BST(ISYMOP))
C
        FILMMR  = 'MMREPOP_0'
        LUMMRE  = -1
        IADRFIL = 1
        NDATA   = N2BST(ISYMOP)
C
        CALL WOPEN2(LUMMRE,FILMMR,64,0)
        CALL GETWA2(LUMMRE,FILMMR,WORK(KINT),IADRFIL,NDATA)
C
        IF (REPTST) THEN
          NBASLO = MAX(NBAST,1)
C
          KINT1 = KWRK2
          KINT2 = KINT1 + N2BST(ISYMOP)
          KWRK3 = KINT2 + N2BST(ISYMOP)
          LWRK3 = LWORK - KWRK3
C
          IF (LWRK3.LT.0) THEN
            CALL QUIT( 'Too little work in CCMM_LTRB Rep. 1a')
          END IF
C
          CALL DZERO(WORK(KINT1),2*N2BST(ISYMOP))
          CALL DCOPY(N2BST(ISYMOP),WORK(KINT),1,WORK(KINT1),1)
C
          DO 668 KR=1,NREPMT
            CALL DGEMM('N','N',NBAST,NBAST,NBAST,
     *                 ONE,WORK(KINT),NBASLO,
     *                 WORK(KINT1),NBASLO,
     *                 ZERO,WORK(KINT2),NBASLO)
C
            CALL DCOPY(N2BST(ISYMOP),WORK(KINT2),1,WORK(KINT),1)
  668     CONTINUE
        END IF
C
        IF (IQM3PR .GT.14) THEN
          TAL1 = DDOT(N2BST(ISYMOP),WORK(KINT),1,WORK(KINT),1)
          TAL2 = SQRT(TAL1)
          WRITE(lupri,*)'Norm of rep. operator:',TAL2
        END IF
C
        KXI   = KWRK2
        KWRK3 = KXI + NTAMP
        LWRK3 = LWORK - KWRK3
C
        IF (LWRK3.LT.0) THEN
          CALL QUIT( 'Too little work in CCMM_LTRB Rep. 2')
        END IF
C
        CALL DZERO(WORK(KXI),NTAMP)
C
        LABEL = 'GIVE INT'
        CALL CC_XKSI(WORK(KXI),LABEL,ISYMTR,0,
     *               WORK(KINT),WORK(KWRK3),LWRK3)
C
        KXI1 = KXI
        KXI2 = KXI  + NTAMP1
C
          TAL3 = DDOT(NT1AM(ISYMTR),WORK(KXI1),1,WORK(KXI1),1)
          TAL4 = DDOT(NT2AM(ISYMTR),WORK(KXI2),1,WORK(KXI2),1)
C
        FAC2 = 1.0D0
        CALL DAXPY(NT1AM(ISYMTR),FAC2,WORK(KXI1),1,
     *             RHO1,1)
        CALL DAXPY(NT2AM(ISYMTR),FAC2,WORK(KXI2),1,
     *                   RHO2,1)
C
        CALL WCLOSE2(LUMMRE,FILMMR, 'KEEP')
C
C test test test
        END IF
C test slut
      END IF
C
C---------------------------------------------------
C     Calculate contribution from the T^gB operator:
C---------------------------------------------------
C
      IF (LR.NE.'0') THEN
C
        KETA    = KWRK1
        KWRK2   = KETA    + NTAMP
        LWRK2   = LWORK   - KWRK2
C
        IF (LWRK2 .LT. 0) 
     &    CALL QUIT( ' Too litle work in CCMM_LTRB 6' )
C
        CALL DZERO(WORK(KETA),NTAMP)
C
        KETA1   = KETA
        KETA2   = KETA + NTAMP1
C
        IF ((LR.EQ.'L').OR.(LR.EQ.'F')) THEN
C
          LIST  = 'L0'
          LABEL = 'GIVE INT'           
C
          CALL CC_ETAC(ISYMTR,LABEL,WORK(KETA),
     *                  LIST,1,0,WORK(KTGB),WORK(KWRK2),LWRK2)
        ELSE IF ((LR.EQ.'R').OR.(LR.EQ.'P')) THEN
C
          LABEL = 'GIVE INT'           
C
          CALL CC_XKSI(WORK(KETA),LABEL,ISYMTR,0,WORK(KTGB), 
     *               WORK(KWRK2),LWRK2)
          IF (LR.EQ.'R')
     *      CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYMTR)
        END IF
C
        IF ( DEBUG .OR.(IQM3PR .GT. 10)) THEN
          RHO1N = DDOT(NT1AM(ISYMTR),WORK(KETA1),1,WORK(KETA1),1)
          RHO2N = DDOT(NT2AM(ISYMTR),WORK(KETA2),1,WORK(KETA2),1)
          WRITE(LUPRI,*) ' Norm af T^gB contribution to LHTR1:', RHO1N
          WRITE(LUPRI,*) ' Norm af T^gB contribution to LHTR2:', RHO2N
        END IF
C
        CALL DAXPY(NT1AM(ISYMTR),ONE,WORK(KETA1),1,RHO1,1)
        CALL DAXPY(NT2AM(ISYMTR),ONE,WORK(KETA2),1,RHO2,1)
C
        IF ( DEBUG .OR.(IQM3PR .GT. 10)) THEN
          RHO1N = DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
          RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
          WRITE(LUPRI,*) ' Norm af RHO1 in CCMM_LTGB:', RHO1N
          WRITE(LUPRI,*) ' Norm af RHO2 in CCMM_LTGB:', RHO2N
        END IF
      END IF
C
      END
C
C********************************************************************
C  /* Deck ccmmint */
      SUBROUTINE CCMMINT(LABEL,LU,RAAOx,Xint,IRREP,
     *                   ISYM,IERR,WORK,LWORK)
C********************************************************************
C
C  Purpose: read property one-electron AO integrals from file or copy
C  from array and transform to CC format.   
C
C       input:   LU    -- Logical unit number for file to be read.
C
C       output:           property AO integrals in coupled cluster
C                         storage scheme (dimension of ouput vector
C                         is N2BST(IRREP)
C                IRREP -- irrep symmetry of the operator
C                ISYM  -- +1 for a symmetric operator
C                         -1 for an antisymmetric operator
C                          0 if unknown
C
C********************************************************************
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccorb.h"
#include "inftap.h"
#include "ccisao.h"
#include "dummy.h"
#include "qm3.h"
C
C local parameters:
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
#if defined (SYS_CRAY)
      REAL CKMXPR
#else
      DOUBLE PRECISION CKMXPR
#endif
      PARAMETER (CKMXPR = 1.0d-12)
C
C input:
C
      CHARACTER*8 LABEL
      CHARACTER*8 RRAO
      INTEGER IRREP, ISYM, IERR, LWORK, LU
C
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL RAAOx(*),Xint(*)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION RAAOx(*),Xint(*)
#endif
C
      LOGICAL LOPENED
      INTEGER IDX, IDXI, IDXJ
      INTEGER KEND0, LEND0, KINTEG, ISYMA, ISYMB, IOLD, INEW
C
C functions:
C
      INTEGER IDAMAX
C
C--------------------------
C     Dynamical allocation:
C-------------------------- 
C
      KINTEG = 1
      KEND0  = KINTEG + N2BASX
      LEND0  = LWORK - KEND0
C
      IF (LEND0 .LT. 0) CALL QUIT('Insufficient memory in CCRRA.')
C 
      CALL DZERO(WORK(KINTEG),N2BASX)
C
      ISYM = +1
C
      IF (LABEL.EQ.'GIVE INT') THEN
         IF (IQM3PR .GT. 15) WRITE (LUPRI,*) ' Starting to copy!' 
         CALL DCOPY(NNBASX,XINT,1,WORK(KEND0),1)
      ELSE
        CALL READT(LU,NNBASX,WORK(KEND0))
      END IF
C
      CALL DSPTGE(NBAST,WORK(KEND0),WORK(KINTEG))
C
      IF (IQM3PR .GT. 99 .OR. LOCDBG) THEN
        CALL AROUND('AOQM3INT: -integrals')
        CALL OUTPUT(WORK(KINTEG),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C-------------------------------------
C Find irrep symmetry of the operator:
C-------------------------------------
C
      IDX = IDAMAX(N2BASX,WORK(KINTEG),1)
C
      IF ( ABS(WORK(KINTEG-1+IDX)) .GT. CKMXPR ) THEN
        IDXI  = (IDX+NBAST-1) / NBAST
        IDXJ  = IDX - (IDXI-1)*NBAST
        IRREP = MULD2H(ISAO(IDXI),ISAO(IDXJ))
      ELSE
        IRREP = 0
        IERR  = -1
C
        IF (IQM3PR .GT. 5) THEN        
          WRITE (LUPRI,*)
     &         'WARNING: irrep symmetry can not be determined.'
          WRITE (LUPRI,*)
     &         'Irrep set to 1!!!! '
        END IF
C
        IRREP = 1
        CALL FLSHFO(LUPRI)
      END IF
C
C---------------------------------------------
C Resort integrals to coupled cluster storage:
C---------------------------------------------
C
      DO ISYMA = 1, NSYM
        ISYMB = MULD2H(IRREP,ISYMA)
        DO A = 1, NBAS(ISYMA)
        DO B = 1, NBAS(ISYMB)
          IOLD = (IBAS(ISYMA)+A-1)*NBAST + (IBAS(ISYMB)+B)
          INEW = IAODIS(ISYMB,ISYMA) + NBAS(ISYMB)*(A-1) + B
          RAAOx(INEW) = WORK(KINTEG-1+IOLD)
        END DO
        END DO
      END DO
C
      RETURN
      END
C
C**************************************************************
C  /* Deck hfint */
      SUBROUTINE HFINT(EINTF,WORK,LWORK)
C**************************************************************
C
C This routine calculates the electrostatic HF/MM electronic 
C interaction energy. Note that the interaction energy
C is calculated using the wave function optimized with
C the induced dipoles determined at the CC level of theory.
C Thus, this is only a pseudo HF result. 
C
C**************************************************************
C      
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccinftap.h"
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
C
      DIMENSION WORK(LWORK)
C
      CHARACTER*8 LABEL
C
C--------------------------
C     Dynamical allocation:
C--------------------------
C
      KDENS = 1
      KNSAO = KDENS + N2BST(ISYMOP)
      KHFNS = KNSAO + N2BST(ISYMOP)
      KWRK2 = KHFNS + NUSITE
      LWRK2 = LWORK - KWRK2
C
      IF (LWRK2 .LT. 0) THEN
        WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
        CALL QUIT('Insufficient memory for CCS AO-density in '//
     &            'HFINT')
      ENDIF
C
      CALL DZERO(WORK(KDENS),N2BST(ISYMOP))
      CALL DZERO(WORK(KNSAO),N2BST(ISYMOP))
      CALL DZERO(WORK(KHFNS),NUSITE)
C
C------------------------------------------
C     First we calculate the HF AO density:
C------------------------------------------
C
      CALL CCS_D1AO(WORK(KDENS),WORK(KWRK2),LWRK2)
C
      IF (IQM3PR .GT. 50) THEN
        CALL AROUND('CCS One electron density in HFINT')
        CALL CC_PRFCKAO(WORK(KDENS),1)
      ENDIF
C
C-----------------------------------------
C     Read integrals from file
C     and perform the calculation of HFNS: 
C-----------------------------------------
C
      LUCCPO = -1
      CALL GPOPEN(LUCCPO,'POTMM','UNKNOWN',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCPO)
C
      LM = 0
      DO 950 I = 1, ISYTP
        IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
          DO 951 J = NSYSBG(I), NSYSED(I)
            DO 952 K = 1,NSISY(I)
              LM = LM +1
              LABEL = 'READ INT'
              CALL CCMMINT(LABEL,LUCCPO,WORK(KNSAO),DUMMY,
     *                     IRREP, ISYM,IERR,
     *                     WORK(KWRK2),LWRK2)
              WORK(KHFNS+LM-1) = DDOT(N2BST(ISYMOP),
     *                           WORK(KNSAO),1,
     *                           WORK(KDENS),1)
  952       CONTINUE 
  951     CONTINUE
        END IF 
  950 CONTINUE
C
      CALL GPCLOSE(LUCCPO,'KEEP')
C
C-------------------
C     Print section:
C-------------------
C
      IF (IQM3PR .GT. 5) THEN
         WRITE(LUPRI,'(//,A)')
     *        ' +======================+'
         WRITE(LUPRI,'(A)')
     *        ' |Site| <N_s>           |'
         WRITE(LUPRI,'(1X,A)')
     *        '+======================+'
         LS = 0
         DO 215 I = 1, ISYTP
           IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
             DO 216 J = NSYSBG(I), NSYSED(I)
               DO 217 K = 1, NSISY(I)
                 LS = LS + 1
                 WRITE(LUPRI,'(1X,A,I3,A,F16.10,A)')
     *                       '| ', LS,'|', WORK(KHFNS+LS-1),' |'
                 WRITE(LUPRI,'(1X,A)')
     *                       '+----------------------+'
  217          CONTINUE
  216        CONTINUE
           END IF
  215    CONTINUE
         WRITE(LUPRI,'(//,A)')
      END IF
C
C-------------------------------
C     Adding each contribution :
C-------------------------------
C
      EINT = 0.0D0
      L = 0
C
      DO 130 I = 1, ISYTP
        IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
          DO 140 J = NSYSBG(I), NSYSED(I)
            DO 150 K = 1, NSISY(I)
              L = L + 1
              EINT = EINT - WORK(KHFNS+L-1)
  150       CONTINUE
  140     CONTINUE
        END IF
  130 CONTINUE
C
      EINTF = EINT
C
      END
C
C**************************************************************
C  /* Deck epert2 */
      SUBROUTINE EPERT2(EPOL,WORK,LWORK)
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccinftap.h"
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
C
      DIMENSION WORK(LWORK)
      DIMENSION EELEC(3,MXQM3)
      DIMENSION FFJ(3)
      LOGICAL   LOINDM
C
      DO 879 I = 1, MXQM3
         DO 880 J = 1, 3
            EELEC(J,I) = 0.0D0
  880    CONTINUE
  879 CONTINUE
C
C--------------------------
C     Dynamical allocation:
C--------------------------
C
      KRAx = 1
      KRAy = KRAx + NCOMS
      KRAz = KRAy + NCOMS
      KOMMSx = KRAz + NCOMS
      KOMMSy = KOMMSx + NCOMS
      KOMMSz = KOMMSy + NCOMS
      KWRK1 = KOMMSz + NCOMS
      LWRK1 = LWORK - KWRK1
C
      CALL DZERO(WORK(KRAx),6*NCOMS)
C
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work space in CC_QM3')
C
      CALL CC_GET31('CC_RA',NCOMS,
     *               WORK(KRAx),WORK(KRAy),WORK(KRAz))
C
      DO 111 I = 1,NCOMS
        EELEC(1,I) = WORK(KRAx + I - 1)
        EELEC(2,I) = WORK(KRAy + I - 1)
        EELEC(3,I) = WORK(KRAz + I - 1)
  111 CONTINUE
C
        FFJ(1) = 0.0D0
        FFJ(2) = 0.0D0
        FFJ(3) = 0.0D0
C
      IF (RELMOM) THEN
        DO 330 IF =1, NFIELD
          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
  330   CONTINUE
C
        WRITE(LUPRI,*)'STATIC ELECTRIC FIELD ADDED (x,y,z)'
        WRITE(LUPRI,*) FFJ(1),FFJ(2),FFJ(3)
      END IF
C
      LOINDM = .TRUE.
      CALL INDMOM(EELEC,LOINDM,FFJ)
      CALL QM3_OBAR(FFJ)
      CALL QM3_OTILDE(OTILDE,FFJ)
C
      CALL CC_GET31('OBARFILE',NCOMS,
     *                 WORK(KOMMSx),WORK(KOMMSy),WORK(KOMMSz))
C
      EQM3T = ZERO
      KS = 0
C
      DO 110 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 120 J = NSYSBG(I), NSYSED(I)
            DO 134 K = 1, NUALIS(I)
              KS = KS+1
              RAT = ZERO
              RAT =  0.5D0 * ALPIMM(I,K) * WORK(KRAx + KS -1) *
     &        (WORK(KRAx + KS -1) + WORK(KOMMSx + KS - 1))
     &        + 0.5D0 * ALPIMM(I,K) * WORK(KRAy + KS -1) *
     &        (WORK(KRAy + KS -1) + WORK(KOMMSy + KS - 1))
     &        + 0.5D0 * ALPIMM(I,K) * WORK(KRAz + KS -1) *
     &        (WORK(KRAz + KS -1) + WORK(KOMMSz + KS - 1))
              EQM3T = EQM3T - RAT
  134       CONTINUE
  120     CONTINUE
        END IF
  110 CONTINUE
C
       IF (IQM3PR .GT. 10) THEN
         DO 777 I=1,NCOMS
           WRITE(LUPRI,*)'WORK(KOMMSx) =',WORK(KOMMSx + I - 1)
           WRITE(LUPRI,*)'WORK(KOMMSy) =',WORK(KOMMSy + I - 1)
           WRITE(LUPRI,*)'WORK(KOMMSz) =',WORK(KOMMSz + I - 1)
  777    CONTINUE
       END IF

      EPOL = EQM3T
      EPOL = EPOL + OTILDE
C
      END
C
C**************************************************************
C  /* Deck cc_put31 */
      SUBROUTINE CC_PUT31(FLNAME,NULOOP,OMMSx,OMMSy,OMMSz)
C**************************************************************
C
#include "implicit.h"
#include "dummy.h"
C
      CHARACTER*(*) FLNAME
      INTEGER   NMBU,NULOOP
      DIMENSION OMMSx(*) , OMMSy(*) , OMMSz(*)
C
      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'UNKNOWN',' ','FORMATTED',IDUMMY,.FALSE.)
C
      REWIND (NMBU)
      LM = 0
C
      DO 820 L = 1,NULOOP
        LM = LM + 1
        WRITE(NMBU,'(I5,3E25.15)') LM,OMMSx(LM),OMMSy(LM),OMMSz(LM)
  820 CONTINUE
C
      REWIND (NMBU)
      CALL GPCLOSE(NMBU,'KEEP')
C
      END
C
C**************************************************************
C  /* Deck cc_get31 */
      SUBROUTINE CC_GET31(FLNAME,NULOOP,OMMSx,OMMSy,OMMSz)
C**************************************************************
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
C
C
      INTEGER   NMBU,NULOOP
      CHARACTER*(*) FLNAME
      DIMENSION OMMSx(*) , OMMSy(*) , OMMSz(*)
C
      LOGICAL FILE_EXIST
C
      INQUIRE(FILE=FLNAME,EXIST=FILE_EXIST)
      IF (.NOT. FILE_EXIST) THEN
         DO LM = 1,NULOOP
            OMMSx(LM) = 0.0D0
            OMMSy(LM) = 0.0D0
            OMMSz(LM) = 0.0D0
         END DO
         GO TO 9000
      END IF
C

      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'OLD',' ','FORMATTED',IDUMMY,.FALSE.)
 
      DO 820 LM = 1,NULOOP

        READ(NMBU,'(I5,3E25.15)',END=700,ERR=700)
     &      LM1, OMMSx(LM), OMMSy(LM), OMMSz(LM)
        IF (LM.NE.LM1) THEN
           write (lupri,*) 'Error in CC_GET31 on ',FLNAME
           write (lupri,*) 'LM, LM1, NULOOP',LM,LM1,NULOOP
           CALL QUIT( 'Error in CC_GET31')
        END IF
  700   CONTINUE
  820 CONTINUE
 
      CALL GPCLOSE(NMBU,'KEEP')
C
 9000 RETURN
C
      END
C
C**************************************************************
C
C  /* Deck ccmm_pbtr */
      SUBROUTINE CCMM_PBTR(RHO1,RHO2,ISYRES,
     *                     LISTL,IDLSTL,CTR1,ISYCTR,
     *                     LISTR,IDLSTR,BTR1,ISYBTR,
     *                     MODEL,WORK,LWORK)
C
C---------------------------------------------------------------
C Calculates contributions from the T^gB , the ^CT^g
C and ^CT^gB operators.
C---------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*),CTR1(*),BTR1(*)
C
      CHARACTER*(*) LISTL,LISTR,LIST*(2)
      CHARACTER*8 LABEL
      CHARACTER MODEL*10
      LOGICAL LEXIST
      INTEGER IDLSTL, IDLSTR, ISYCTR, ISYBTR, ISYBC
C
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCMM_PBTR: ISYRES =',ISYRES
        WRITE(LUPRI,*)'CCMM_PBTR: ISYCTR =',ISYCTR
        WRITE(LUPRI,*)'CCMM_PBTR: ISYBTR =',ISYBTR
        WRITE(LUPRI,*)'CCMM_PBTR: LISTL =',LISTL
        WRITE(LUPRI,*)'CCMM_PBTR: LISTR =',LISTR
        WRITE(LUPRI,*)'CCMM_PBTR: IDLSTL =',IDLSTL
        WRITE(LUPRI,*)'CCMM_PBTR: IDLSTR =',IDLSTR
        CALL FLSHFO(LUPRI)
      ENDIF
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CCMM_TGB:
C----------------------------------------------
C
      IF (LOSPC) RETURN
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
      IF (CCS) THEN
        NT2AM(ISYBTR)  = IZERO
        NT2AM(ISYCTR)  = IZERO
        NT2AM(ISYBC)   = IZERO
      END IF
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
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_PBTR, 1')
C
      CALL DZERO(WORK(KTGB),N2BST(ISYBTR))
C
      CALL CCMM_TGB(BTR1,ISYBTR,LISTR,IDLSTR,WORK(KTGB),'ET',
     *              MODEL,WORK(KWRK1),LWRK1)
C
C     Symmetry:
C     ---------
C
      ISYBC = MULD2H(ISYBTR,ISYCTR)
C
      IF (ISYBC .NE. ISYRES) THEN
        CALL QUIT( 'Symmetry problem in CCMM_PBTR')
      END IF
C
      KETA    = KWRK1
      KWRK2   = KETA + NAMPF
      LWRK2   = LWORK   - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_PBTR, 2')
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
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_PBTR, 3')
C
      CALL DZERO(WORK(KTGB),N2BST(ISYBC))
C
      CALL CCMM_CTGB(LISTL,IDLSTL,CTR1,ISYCTR,
     *               LISTR,IDLSTR,BTR1,ISYBTR,
     *               WORK(KTGB),MODEL,WORK(KWRK1),LWRK1)
C
      KETA    = KWRK1
      KWRK2   = KETA + NAMPF
      LWRK2   = LWORK   - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_PBTR, 4')
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
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_PBTR, 5')
C
      CALL DZERO(WORK(KTGB),N2BST(ISYCTR))
C
      CALL CCMM_TGB(CTR1,ISYCTR,LISTL,IDLSTL,WORK(KTGB),'XI',
     *              MODEL,WORK(KWRK1),LWRK1)
C
C     Symmetry:
C     ---------
C
      ISYBC = MULD2H(ISYBTR,ISYCTR)
C
      IF (ISYBC .NE. ISYRES) THEN
        CALL QUIT( 'Symmetry problem in CCMM_PBTR')
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
      CALL DAXPY(NT1AM(ISYBC),ONE,WORK(KWRK1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYBC),ONE,
     *           WORK(KWRK1+NT1AM(ISYBC)),1,RHO2,1)
C
      END
*******************************************************************
C  /* Deck ccmm_tgb */
      SUBROUTINE CCMM_TGB(CTR1,ISYMTR,LISTIN,IDLIST,TGB,LR,
     *                    MODEL,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C   IF (LR.EQ.'ET') then the following effective operator is constructed
C   T^gB = - Sum_a * Sum_sigma
C             <Lambda|[Rr_a,Tau_sigma]|CC> t^B_sigma Rr_a
C
C   IF (LR.EQ.'XI') then the following effective operator is constructed
C   ^BT^g = - Sum_a * Sum_sigma
C             t^B_sigma <bar sigma|Rr_a|CC> Rr_a 
C
C   JK+OC, marts 03
C   KS, modifed to allow for J term needed from the left. I.e when cG operator is needed
C   in the NYQMMM formulation. This routine is not used for calculating Gc term (LR=ET) - see
C   instead CCMM_TRANSFORMER in cc_qmmm
C-----------------------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
C
      DIMENSION WORK(LWORK),CTR1(*)
      DIMENSION TGB(*)
C
      INTEGER ISYMTR, IDLIST
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      CHARACTER*(*) LISTIN
      CHARACTER MODEL*10
C
      CHARACTER*8 LABEL,LR*(2),LIST
      CHARACTER*(8) FILMME, FILMMX
      LOGICAL LEXIST, LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
C
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCMM_TGB : ISYMTR:', ISYMTR
      ENDIF
C
C----------------------------------------------
C     If all models are SPC 
C     -> RETURN from CCMM_TGB:
C----------------------------------------------
C
      IF (LOSPC) RETURN
C
C---------------------
C     Init parameters.
C---------------------
C
      NTAMP1  = NT1AM(ISYMTR)
      NTAMP2  = NT2AM(ISYMTR)
      NTAMP   = NTAMP1 + NTAMP2
C
      IF (CCS)  NT2AM(ISYMTR) = IZERO
C
      IF (DISCEX) THEN
        LUMMET = -1
        LUMMXI = -1
        FILMME = 'CCMM_ETA'
        FILMMX = 'CCMM__XI'
        CALL WOPEN2(LUMMET, FILMME, 64, 0)
        CALL WOPEN2(LUMMXI, FILMMX, 64, 0)
      END IF
C
C----------------------------------------
C     Readin trial vector from file.
C----------------------------------------
C
      KT2AMPA = 1
      KWRK1   = KT2AMPA +  NT2AM(ISYMTR)
      LWRK1   = LWORK - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insuff. work in CCMM_TGB 1')
      END IF
C
      IOPT = 2
      CALL CC_RDRSP(LISTIN,IDLIST,ISYMTR,IOPT,MODEL,
     *              WORK(KDUM),WORK(KT2AMPA))

      IF (LOCDEB) THEN
        RHO2N = DDOT(NT2AM(ISYMTR),WORK(KT2AMPA),1,WORK(KT2AMPA),1)
        WRITE(LUPRI,*)'Norm af B2AM in CCMM_TGB on input:,1', RHO2N
      END IF
C
      IF (.NOT. NYQMMM) THEN
C------------------------------------
C     Dynamical allocation for CCMM :
C------------------------------------
C
        KRAAOx = KWRK1
        KRAAOy = KRAAOx + N2BST(ISYMTR)
        KRAAOz = KRAAOy + N2BST(ISYMTR)
        KWRK2  = KRAAOz + N2BST(ISYMTR)
        LWRK2  = LWORK   - KWRK2
C  
        IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_TGB, 2')
C  
        KXIx  = KWRK2
        KXIy  = KXIx + NTAMP
        KXIz  = KXIy + NTAMP
        KWRK3 = KXIz  + NTAMP
        LWRK3   = LWORK   - KWRK3
C  
        IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CCMM_TGB, 3')
C  
        CALL DZERO(WORK(KRAAOx),3*N2BST(ISYMTR))
        CALL DZERO(WORK(KXIx),3*NTAMP)
C  
C  ----------------------------
C       Read Rra in AO basis:
C  ----------------------------
C  
        LUCCEF = -1
        CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND(LUCCEF)
C  
        LM = 0
        IADR1 = 1
        IADR2 = 1
        LEN = NTAMP
C  
        DO 700 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 701 J = NSYSBG(I), NSYSED(I)
              DO 702 K = 1, NUALIS(I)
C  
                LM = LM + 1
                LABEL = 'READ INT'
C  
                CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOx),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
                CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOy),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
                CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOz),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
C  
C---------------------------------------------------------------------
C               (LR.EQ.'ET') Calculate eta^T_lm
C               (LR.EQ.'XI') Calculate xi^T_lm
C               Only if MDLWRD(I) = 'SPC_E01' are we going to calculate
C               contributions from the T^g operator to the response
C               equations.
C---------------------------------------------------------------------
C  
                IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
                  IF (.NOT. (DISCEX)) THEN
                    IF (LR.EQ.'ET') THEN
                      LABEL = 'GIVE INT'
                      LIST  = 'L0'
                      IDLINO = 1
C  
                      CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIx),
     *                   LIST,IDLINO,0,WORK(KRAAOx),WORK(KWRK3),LWRK3)
                      CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIy),
     *                   LIST,IDLINO,0,WORK(KRAAOy),WORK(KWRK3),LWRK3)
                      CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIz),
     *                   LIST,IDLINO,0,WORK(KRAAOz),WORK(KWRK3),LWRK3)
C  
                    ELSE IF (LR.EQ.'XI') THEN 
                      LABEL = 'GIVE INT'
C  
                    CALL CC_XKSI(WORK(KXIx),LABEL,ISYMTR,0,WORK(KRAAOx),
     *                           WORK(KWRK3),LWRK3)
                    CALL CC_XKSI(WORK(KXIy),LABEL,ISYMTR,0,WORK(KRAAOy),
     *                           WORK(KWRK3),LWRK3)
                    CALL CC_XKSI(WORK(KXIz),LABEL,ISYMTR,0,WORK(KRAAOz),
     *                           WORK(KWRK3),LWRK3) 
C  
                    END IF
                  ELSE 
                    IF (LR.EQ.'ET') THEN
                      CALL GETWA2(LUMMET,FILMME,WORK(KXIx),IADR1,LEN)
                      IADR1 = IADR1 + LEN
                      CALL GETWA2(LUMMET,FILMME,WORK(KXIy),IADR1,LEN)
                      IADR1 = IADR1 + LEN
                      CALL GETWA2(LUMMET,FILMME,WORK(KXIz),IADR1,LEN)
                      IADR1 = IADR1 + LEN
                    ELSE IF (LR.EQ.'XI') THEN
                      CALL GETWA2(LUMMXI,FILMMX,WORK(KXIx),IADR2,LEN)
                      IADR2 = IADR2 + LEN
                      CALL GETWA2(LUMMXI,FILMMX,WORK(KXIy),IADR2,LEN)
                      IADR2 = IADR2 + LEN
                      CALL GETWA2(LUMMXI,FILMMX,WORK(KXIz),IADR2,LEN)
                      IADR2 = IADR2 + LEN
                    END IF
                  END IF
C  
C--------------------------------------------------------------
C                 Contract with trial vector:
C                 (If LR.EQ.'ET' then from the right hand side)
C                 (If LR.EQ.'XI' then from the left hand side)
C--------------------------------------------------------------
C  
                  KXI1x = KXIx
                  KXI2x = KXIx + NTAMP1
                  BXILMD1x= DDOT(NTAMP1,CTR1,1,WORK(KXI1x),1)
                  BXILMD2x= DDOT(NTAMP2,WORK(KT2AMPA),1,WORK(KXI2x),1)
                  BXILMDx = BXILMD1x + BXILMD2x
C  
                  KXI1y = KXIy
                  KXI2y = KXIy + NTAMP1
                  BXILMD1y = DDOT(NTAMP1,CTR1,1,WORK(KXI1y),1)
                  BXILMD2y = DDOT(NTAMP2,WORK(KT2AMPA),1,WORK(KXI2y),1)
                  BXILMDy = BXILMD1y + BXILMD2y
C  
                  KXI1z = KXIz
                  KXI2z = KXIz  + NTAMP1
                  BXILMD1z = DDOT(NTAMP1,CTR1,1,WORK(KXI1z),1)
                  BXILMD2z = DDOT(NTAMP2,WORK(KT2AMPA),1,WORK(KXI2z),1)
                  BXILMDz  = BXILMD1z + BXILMD2z
C  
                  FACx =  - ALPIMM(I,K) * (BXILMDx)
                  FACy =  - ALPIMM(I,K) * (BXILMDy)
                  FACz =  - ALPIMM(I,K) * (BXILMDz)
C  
                  IF (IQM3PR .GT. 5) THEN
                    WRITE(LUPRI,*)'---------------------------' //
     *                          '-------------------'
                    WRITE(LUPRI,*)'The factors calculated at cent.' //
     *                          ' of mass no.',LM
                    WRITE(LUPRI,*)'---------------------------' //
     *                          '-------------------'
                    WRITE(LUPRI,*)'FACx=',FACx
                    WRITE(LUPRI,*)'FACy=',FACy
                    WRITE(LUPRI,*)'FACz=',FACz
                    WRITE(LUPRI,*)'---------------------------' //
     *                          '-------------------'
                 END IF
C  

C  --------------------------------------------------------------
C          Put factor on integrals 
C  --------------------------------------------------------------
C  
                    CALL DAXPY(N2BST(ISYMTR),FACx,WORK(KRAAOx),1,
     *                       TGB,1)
                    CALL DAXPY(N2BST(ISYMTR),FACy,WORK(KRAAOy),1,
     *                       TGB,1)
                    CALL DAXPY(N2BST(ISYMTR),FACz,WORK(KRAAOz),1,
     *                       TGB,1)
C  
                END IF
  702       CONTINUE
  701     CONTINUE
        END IF
  700 CONTINUE
C
      CALL GPCLOSE(LUCCEF,'KEEP')
C
      ELSE IF (NYQMMM) THEN

        IF (LR .NE. 'XI' ) CALL QUIT('CCMM_TGB xi flag not allowed')
        KDENS = KWRK1
        KWRK2 = KDENS + N2BST(ISYMTR)
        LWRK2 = LWORK - KWRK2

        CALL DZERO(WORK(KDENS),N2BST(ISYMTR))

C----------------------------
C     Construct auxiliary densities
C----------------------------
        CALL CCMM_D1AO(WORK(KDENS),CTR1,WORK(KT2AMPA),MODEL,'L',
     *                 LISDUM,IDLDUM,ISYDUM,
     *                 WORK(KWRK2),LWRK2)
C----------------------------
C     Construct effective G operator
C----------------------------

        CALL CCMM_GEFF(WORK(KDENS),TGB(1),WORK(KWRK2),LWRK2)

      END IF
      IF (DISCEX) THEN
        CALL WCLOSE2(LUMMET, FILMME, 'KEEP')
        CALL WCLOSE2(LUMMXI, FILMMX, 'KEEP')
      END IF
C
      END
************************************************************
C  /* Deck ccmm_ctgb */
      SUBROUTINE CCMM_CTGB(LISTL,IDLSTL,CTR1,ISYCTR,
     *                     LISTR,IDLSTR,BTR1,ISYBTR,
     *                     TGB,MODEL,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C   Construcs effective operator equal to
C   ^CT^gB = - Sum_a * Sum_sigma Sum_my
C            t^C_my <bar my|[Rr_a,Tau_sigma]|CC> t^B_sigma * Rr_a
C
C   JK+OC, marts 03
C-----------------------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
C
      DIMENSION WORK(LWORK),BTR1(*),CTR1(*)
      DIMENSION TGB(*)
C
      CHARACTER*8 LABEL, LIST*(2)
      CHARACTER*(*) LISTL, LISTR
      CHARACTER MODEL*10
      INTEGER IDLSTL, IDLSTR, ISYCTR, ISYBTR
      LOGICAL LEXIST
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCMM_CTGB: ISYCTR =',ISYCTR
        WRITE(LUPRI,*)'CCMM_CTGB: ISYBTR =',ISYBTR
        WRITE(LUPRI,*)'CCMM_CTGB: LISTL =',LISTL
        WRITE(LUPRI,*)'CCMM_CTGB: LISTR =',LISTR
        WRITE(LUPRI,*)'CCMM_CTGB: IDLSTL =',IDLSTL
        WRITE(LUPRI,*)'CCMM_CTGB: IDLSTR =',IDLSTR
        CALL FLSHFO(LUPRI)
      ENDIF
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CCMM_TGB:
C----------------------------------------------
C
      IF (LOSPC) RETURN
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
      IF (CCS) THEN
        NT2AM(ISYBTR) = IZERO
        NT2AM(ISYCTR) = IZERO
        NT2AM(ISYBC)  = IZERO
      END IF 

C-----------------------------------------------
C    Readin response vector (B, right) from file
C-----------------------------------------------
C
      KT2AMPA = 1
      KWRK1   = KT2AMPA +  NT2AM(ISYBTR)
      LWRK1   = LWORK - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insuff. work in CCMM_TGB')
      END IF
C
      IOPT = 2
      CALL CC_RDRSP(LISTR,IDLSTR,ISYBTR,IOPT,MODEL,
     *              WORK(KWRK1),WORK(KT2AMPA))
C
C------------------------------------
C     Dynamical allocation for CCMM :
C------------------------------------
C
      KRAAOx = KWRK1
      KRAAOy = KRAAOx + N2BST(ISYBC)
      KRAAOz = KRAAOy + N2BST(ISYBC)
      KWRK2  = KRAAOz + N2BST(ISYBC)
      LWRK2  = LWORK   - KWRK2
C
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_CTGB, 2')
C
      KXIx  = KWRK2
      KXIy  = KXIx + NAMPB
      KXIz  = KXIy + NAMPB
      KWRK3 = KXIz  + NAMPB
      LWRK3   = LWORK   - KWRK3
C
      IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CCMM_CTGB, 3')
C
      CALL DZERO(WORK(KRAAOx),3*N2BST(ISYBC))
      CALL DZERO(WORK(KXIx),3*NAMPB)
C
C----------------------------
C     Read Rra in AO basis:
C----------------------------
C
      LUCCPO = -1
      CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCEF)
C
      LM = 0
C
      DO 700 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 701 J = NSYSBG(I), NSYSED(I)
            DO 702 K = 1, NUALIS(I)
C
              LM = LM + 1
              LABEL = 'READ INT'
C
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOx),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOy),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOz),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK3),LWRK3)
C
C---------------------------------------------------------------------
C        Calculate eta^Rr_a
C        (LIST = 'LE' ensures that the HF part of eta^Rr_a is skipped.
C---------------------------------------------------------------------
C
              IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
                LABEL = 'GIVE INT'
                CALL CC_ETAC(ISYBC,LABEL,WORK(KXIx),
     *               LISTL,IDLSTL,0,WORK(KRAAOx),WORK(KWRK3),LWRK3)
                CALL CC_ETAC(ISYBC,LABEL,WORK(KXIy),
     *               LISTL,IDLSTL,0,WORK(KRAAOy),WORK(KWRK3),LWRK3)                
                CALL CC_ETAC(ISYBC,LABEL,WORK(KXIz),
     *               LISTL,IDLSTL,0,WORK(KRAAOz),WORK(KWRK3),LWRK3)
C
C---------------------------------------------------------------
C        Contract with trial vector B (from the right hand side)
C---------------------------------------------------------------
C
                KXI1x = KXIx
                KXI2x = KXIx  + NT1AM(ISYBTR)
                BXILMD1x= DDOT(NT1AM(ISYBTR),BTR1,1,WORK(KXI1x),1)
                BXILMD2x= DDOT(NT2AM(ISYBTR),WORK(KT2AMPA),
     *                         1,WORK(KXI2x),1)
                BXILMDx = BXILMD1x + BXILMD2x
C
                KXI1y = KXIy
                KXI2y = KXIy  + NT1AM(ISYBTR)
                BXILMD1y= DDOT(NT1AM(ISYBTR),BTR1,1,WORK(KXI1y),1)
                BXILMD2y= DDOT(NT2AM(ISYBTR),WORK(KT2AMPA),
     *                         1,WORK(KXI2y),1)
                BXILMDy = BXILMD1y + BXILMD2y
C
                KXI1z = KXIz
                KXI2z = KXIz  + NT1AM(ISYBTR)
                BXILMD1z= DDOT(NT1AM(ISYBTR),BTR1,1,WORK(KXI1z),1)
                BXILMD2z= DDOT(NT2AM(ISYBTR),WORK(KT2AMPA),
     *                         1,WORK(KXI2z),1)
                BXILMDz = BXILMD1z + BXILMD2z
C
                FACx =  - ALPIMM(I,K) * (BXILMDx)
                FACy =  - ALPIMM(I,K) * (BXILMDy)
                FACz =  - ALPIMM(I,K) * (BXILMDz)
C
                IF (IQM3PR .GT. 5) THEN
                  WRITE(LUPRI,*)'---------------------------' //
     *                          '-------------------'
                  WRITE(LUPRI,*)'The factors calculated at cent.' //
     *                          ' of mass no.',LM
                  WRITE(LUPRI,*)'---------------------------' //
     *                          '-------------------'
                  WRITE(LUPRI,*)'FACx=',FACx
                  WRITE(LUPRI,*)'FACy=',FACy
                  WRITE(LUPRI,*)'FACz=',FACz
                  WRITE(LUPRI,*)'---------------------------' //
     *                          '-------------------'
                END IF
C
C--------------------------------------------------------------
C               Put factor on integrals 
C--------------------------------------------------------------
C
                CALL DAXPY(N2BST(ISYBC),FACx,WORK(KRAAOx),1,
     *                     TGB,1)
                CALL DAXPY(N2BST(ISYBC),FACy,WORK(KRAAOy),1,
     *                     TGB,1)
                CALL DAXPY(N2BST(ISYBC),FACz,WORK(KRAAOz),1,
     *                     TGB,1)
C
C---------------------------------------------------------------
C
              END IF
  702       CONTINUE
  701     CONTINUE
        END IF
  700 CONTINUE
C
      CALL GPCLOSE(LUCCEF,'KEEP')
      END
C*******************************************************************
C  /* Deck ccmm_btr */
      SUBROUTINE CCMM_BTR(RHO1,RHO2,ISYMAB,
     *                    LISTA,IDLSTA,ISYMA,
     *                    LISTB,IDLSTB,ISYMB,
     *                    MODEL,WORK,LWORK)
C
C     Calculates the CC/MM contribution to the B-matrix 
C     transformation.
C     JK+OC, marts 03
C---------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*)
C
      CHARACTER*(*) LISTA,LISTB
      CHARACTER*(3) LISTC
      CHARACTER*8 LABEL
      CHARACTER MODEL*10
      LOGICAL LEXIST, LSAME, LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
      INTEGER IDLSTA, IDLSTB, ISYMAB, ISYMA, ISYMB
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CCMM_TGB:
C----------------------------------------------
C
      IF (LOSPC) RETURN
C
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCMM_BTR: ISYMAB =',ISYMAB
        WRITE(LUPRI,*)'CCMM_BTR: ISYMA  =',ISYMA
        WRITE(LUPRI,*)'CCMM_BTR: ISYMB  =',ISYMB
        WRITE(LUPRI,*)'CCMM_BTR: LISTA  =',LISTA
        WRITE(LUPRI,*)'CCMM_BTR: LISTB  =',LISTB
        WRITE(LUPRI,*)'CCMM_BTR: IDLSTA =',IDLSTA
        WRITE(LUPRI,*)'CCMM_BTR: IDLSTB =',IDLSTB
        CALL FLSHFO(LUPRI)
      ENDIF
C
      IF (CCS) THEN
        NT2AM(ISYMAB) = IZERO
        NT2AM(ISYMA)  = IZERO
        NT2AM(ISYMB)  = IZERO
      END IF
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
         CALL QUIT('Insuff. work in CCMM_BTR 1')
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

      IF (LOCDEB) THEN
        RHO1N = DDOT(NT1AM(ISYMAB),RHO1,1,RHO1,1)
        RHO2N = DDOT(NT2AM(ISYMAB),RHO2,1,RHO2,1)
        WRITE(LUPRI,*)'Norm of RHO1 in CCMM_QRTRANS on input:,1', RHO1N
        WRITE(LUPRI,*)'Norm of RHO2 in CCMM_QRTRANS on input:,1', RHO2N
        RHO1N = DDOT(NT1AM(ISYMA),WORK(KT1AMPA),1,WORK(KT1AMPA),1)
        WRITE(LUPRI,*)'Norm af B1AM in CCMM_QRTRANS on input:,1', RHO1N
      END IF
C
      CALL CCMM_TGB(WORK(KT1AMPA),ISYMA,LISTA,IDLSTA,WORK(KTGA),'ET',
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

      IF (LOCDEB) THEN
         TAL1= DDOT(NT1AM(ISYMAB),WORK(KWRK1),1,WORK(KWRK1),1)
         TAL2= DDOT(NT2AM(ISYMAB),WORK(KWRK1+NT1AM(ISYMAB)),1,
     *              WORK(KWRK1+NT1AM(ISYMAB)),1)
         WRITE(LUPRI,*) 'Printing first contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
         TAL1= DDOT(NT1AM(ISYMAB),RHO1,1,RHO1,1)
         TAL2= DDOT(NT2AM(ISYMAB),RHO2,1,RHO2,1)
         WRITE(LUPRI,*) 'Printing RHO total. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF
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
           CALL QUIT('Insuff. work in CCMM_BTR 2')
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

        IF (LOCDEB) THEN
         RHO1N = DDOT(NT1AM(ISYMAB),RHO1,1,RHO1,1)
         RHO2N = DDOT(NT2AM(ISYMAB),RHO2,1,RHO2,1)
         WRITE(LUPRI,*)'Norm of RHO1 in CCMM_QRTRANS on input:,1', RHO1N
         WRITE(LUPRI,*)'Norm of RHO2 in CCMM_QRTRANS on input:,1', RHO2N
         RHO1N = DDOT(NT1AM(ISYMB),WORK(KT1AMPA),1,WORK(KT1AMPA),1)
         WRITE(LUPRI,*)'Norm af B1AM in CCMM_QRTRANS on input:,1', RHO1N
        END IF

        CALL CCMM_TGB(WORK(KT1AMPB),ISYMB,LISTB,IDLSTB,WORK(KTGB),'ET',
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

        IF (LOCDEB) THEN
         TAL1= DDOT(NT1AM(ISYMAB),WORK(KWRK1),1,WORK(KWRK1),1)
         TAL2= DDOT(NT2AM(ISYMAB),WORK(KWRK1+NT1AM(ISYMAB)),1,
     *              WORK(KWRK1+NT1AM(ISYMAB)),1)
         WRITE(LUPRI,*) 'Printing second contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2

         TAL1= DDOT(NT1AM(ISYMAB),RHO1,1,RHO1,1)
         TAL2= DDOT(NT2AM(ISYMAB),RHO2,1,RHO2,1)
         WRITE(LUPRI,*) 'Printing RHO total. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
        END IF
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
C      TAL2= DDOT(N2BST(ISYMAB),WORK(KTGAB),1,WORK(KTGAB),1)
C      WRITE(LUPRI,*) 'Printing TGAB before build. 
C     &                     Norm2 of operator: ', TAL2
      CALL CCMM_TGCB(ISYMAB,LISTC,LISTA,LISTB,
     *               IDLSTC,IDLSTA,IDLSTB,
     *               ISYMC,ISYMA,ISYMB,WORK(KTGAB),
     *               MODEL,WORK(KWRK1),LWRK1)
C
C      TAL2= DDOT(N2BST(ISYMAB),WORK(KTGAB),1,WORK(KTGAB),1)
C      WRITE(LUPRI,*) 'Printing TGAB after build. 
C     &                     Norm2 of operator: ', TAL2

      NAMPF = NT1AM(ISYMAB) + NT2AM(ISYMAB)
C
      KXI   = KWRK1
      KWRK2 = KXI + NAMPF
      LWRK2 = LWORK  - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_BTR, 4')
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

      IF (LOCDEB) THEN
         TAL1= DDOT(NT1AM(ISYMAB),WORK(KXI1),1,WORK(KXI1),1)
         TAL2= DDOT(NT2AM(ISYMAB),WORK(KXI2),1,WORK(KXI2),1)
         WRITE(LUPRI,*) 'Printing XI^TGCB contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF
C
      CALL DAXPY(NT1AM(ISYMAB),ONE,WORK(KXI1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYMAB),ONE,WORK(KXI2),1,RHO2,1)

      END
***********************************************************
C  /* Deck ccmm_tgcb */
      SUBROUTINE CCMM_TGCB(ISYRES,LISTA,LISTB,LISTC,
     *                     IDLSTA,IDLSTB,IDLSTC,
     *                     ISYMTA,ISYMTB,ISYMTC,TGCB,
     *                     MODEL,WORK,LWORK)
C
C     Constructs effective operator equal to
C     T^gCB = - sum_a <Lambda|[[Rr_a,C],B]|CC> Rr_a
C     JK+OC, marts 03
C-----------------------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
C
      DIMENSION WORK(LWORK)
      DIMENSION TGCB(*)
C
      INTEGER KDUM
      INTEGER IDLSTA,IDLSTB,IDLSTC,ISYMTC,ISYMTA,ISYMTB,ISYRES
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
      CHARACTER*(*) LISTA,LISTB,LISTC
      CHARACTER MODEL*10
C
      CHARACTER*8 LABEL
      LOGICAL LEXIST
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CCMM_TGB:
C----------------------------------------------
C
      IF (LOSPC) RETURN
C
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
        CALL QUIT('Insuff. work in CCMM_TGCB 1')
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
        CALL QUIT( 'Symmetry problem in CCMM_TGCB')
      END IF
C
C------------------------------------
C     Dynamical allocation for CCMM :
C------------------------------------
C
      KRAAOx = KWRK1
      KRAAOy = KRAAOx + N2BST(ISYCB)
      KRAAOz = KRAAOy + N2BST(ISYCB)
      KWRK2  = KRAAOz + N2BST(ISYCB)
      LWRK2  = LWORK   - KWRK2
C
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_CTGB, 2')
C
      CALL DZERO(WORK(KRAAOx),3*N2BST(ISYCB))
C
C----------------------------
C     Read Rra in AO basis:
C----------------------------
C
      LUCCPO = -1
      CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCEF)
C
      LM = 0
C
      DO 700 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 701 J = NSYSBG(I), NSYSED(I)
            DO 702 K = 1, NUALIS(I)
C
              LM = LM + 1
              LABEL = 'READ INT'
C
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOx),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK2),LWRK2)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOy),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK2),LWRK2)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOz),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK2),LWRK2)
C
C-------------------------------------------------------------------
C
              IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
C
                LABEL = 'GIVE INT'
                CALL CCLR_FA(LABEL,ISYCB,LISTC,IDLSTC,
     &               LISTA,IDLSTA,WORK(KRAAOx),WORK(KWRK2),LWRK2)
C
                CBDOT1x = DDOT(NT1AM(ISYMTB),WORK(KWRK2),
     &                         1,WORK(KT1AMB),1)
                CBDOT2x = DDOT(NT2AM(ISYMTB),WORK(KWRK2+NT1AM(ISYMTB)),
     &                         1,WORK(KT2AMB),1)
                CBDOTx = CBDOT1x + CBDOT2x
                FACx =  - ALPIMM(I,K) * (CBDOTx)
                CALL DAXPY(N2BST(ISYCB),FACx,WORK(KRAAOx),1,
     *                     TGCB,1)
C
                LABEL = 'GIVE INT'
                CALL CCLR_FA(LABEL,ISYCB,LISTC,IDLSTC,
     &               LISTA,IDLSTA,WORK(KRAAOy),WORK(KWRK2),LWRK2)
C
                CBDOT1y = DDOT(NT1AM(ISYMTB),WORK(KWRK2),
     &                         1,WORK(KT1AMB),1)
                CBDOT2y = DDOT(NT2AM(ISYMTB),WORK(KWRK2+NT1AM(ISYMTB)),
     &                         1,WORK(KT2AMB),1)
                CBDOTy = CBDOT1y + CBDOT2y
                FACy =  - ALPIMM(I,K) * (CBDOTy)
                CALL DAXPY(N2BST(ISYCB),FACy,WORK(KRAAOy),1,
     *                     TGCB,1)
C
                LABEL = 'GIVE INT'
                CALL CCLR_FA(LABEL,ISYCB,LISTC,IDLSTC,
     &               LISTA,IDLSTA,WORK(KRAAOz),WORK(KWRK2),LWRK2)
C
                CBDOT1z = DDOT(NT1AM(ISYMTB),WORK(KWRK2),
     &                         1,WORK(KT1AMB),1)
                CBDOT2z = DDOT(NT2AM(ISYMTB),WORK(KWRK2+NT1AM(ISYMTB)),
     &                         1,WORK(KT2AMB),1)
                CBDOTz = CBDOT1z + CBDOT2z
                FACz =  - ALPIMM(I,K) * (CBDOTz)
                CALL DAXPY(N2BST(ISYCB),FACz,WORK(KRAAOz),1,
     *                     TGCB,1)
C
              END IF
  702       CONTINUE
  701     CONTINUE
        END IF
  700 CONTINUE
C
      CALL GPCLOSE(LUCCEF,'KEEP')
      END
C****************************************************************
C  /* Deck ccsl_gtr */
      SUBROUTINE CCMM_GTR(RHO1,RHO2,ISYRES,
     *                    LISTA,IDLSTA,ISYMTA,
     *                    LISTB,IDLSTB,ISYMTB,
     *                    LISTC,IDLSTC,ISYMTC,
     *                    MODEL,WORK,LWORK)
C
C-----------------------------------------------------------------
C JK+OC, marts.03
C       Calculates the contribution to the CC/MM G matrx
C       transformations.
C
C       G_my,ny,sigma = 1/2 P(my,ny,sigma)
C                       <Lambda|[[T^g,mu,tau_ny],tau_sigma]|CC >
C
C       The permutation creates 3 different contributions.
C
C------------------------------------------------------------------
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*)
C
      CHARACTER*(*) LISTA, LISTB, LISTC
      CHARACTER*8 LABEL
      CHARACTER MODEL*10
      LOGICAL LEXIST, LSAME, LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
      INTEGER ISYCB,ISYRES,ISYMTA,ISYMTB,ISYMTC
      INTEGER IDLSTA,IDLSTB,IDLSTC
      INTEGER FACTSL
C
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CCMM_TGB:
C----------------------------------------------
C
      IF (LOSPC) RETURN
C
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*)'CCMM_GTR: ISYRES =',ISYRES
        WRITE(LUPRI,*)'CCMM_GTR: ISYMTA =',ISYMTA
        WRITE(LUPRI,*)'CCMM_GTR: ISYMTB =',ISYMTB
        WRITE(LUPRI,*)'CCMM_GTR: ISYMTC =',ISYMTC
        WRITE(LUPRI,*)'CCMM_GTR: LISTA =',LISTA
        WRITE(LUPRI,*)'CCMM_GTR: LISTB =',LISTB
        WRITE(LUPRI,*)'CCMM_GTR: LISTC =',LISTC
        WRITE(LUPRI,*)'CCMM_GTR: IDLSTA =',IDLSTA
        WRITE(LUPRI,*)'CCMM_GTR: IDLSTB =',IDLSTB
        WRITE(LUPRI,*)'CCMM_GTR: IDLSTC =',IDLSTC
        CALL FLSHFO(LUPRI)
      ENDIF
C
C     Symmetry:
C     ---------
C
      ISYCB = MULD2H(ISYMTB,ISYMTC)
C
      IF (ISYCB .NE. ISYRES) THEN
        CALL QUIT( 'Symmetry problem in CCMM_GTR 1')
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
         CALL QUIT('Insuff. work in CCMM_GTR 1')
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

      IF (LOCDEB) THEN
         RHO1N = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
         RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
         WRITE(LUPRI,*)'Norm of RHO1 in CCMM_QRTRANS on input:,1', RHO1N
         WRITE(LUPRI,*)'Norm of RHO2 in CCMM_QRTRANS on input:,1', RHO2N
         RHO1N = DDOT(NT1AM(ISYMTB),WORK(KT1AMPA),1,WORK(KT1AMPA),1)
         WRITE(LUPRI,*)'Norm af B1AM in CCMM_QRTRANS on input:,1', RHO1N
      END IF

      CALL CCMM_TGB(WORK(KT1AMPA),ISYMTB,LISTB,IDLSTB,WORK(KTGB),'ET',
     *              MODEL,WORK(KWRK1),LWRK1)
C
      LABEL = 'GIVE INT'
      CALL CCLR_FA(LABEL,ISYMTB,LISTC,IDLSTC,
     &             LISTA,IDLSTA,WORK(KTGB),WORK(KWRK1),LWRK1)
C
      LSAME  =  (LISTC.EQ.LISTB  .AND. IDLSTC.EQ.IDLSTB)
      IF (LSAME) THEN
        FACTSLV = TWO
      ELSE
        FACTSLV = ONE
      END IF
C
      IF (LOCDEB) THEN
         TAL1= DDOT(NT1AM(ISYCB),WORK(KWRK1),1,WORK(KWRK1),1)
         TAL2= DDOT(NT2AM(ISYCB),WORK(KWRK1+NT1AM(ISYCB)),1,
     *              WORK(KWRK1+NT1AM(ISYCB)),1)
         WRITE(LUPRI,*) 'Printing transformed first contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF
C        
      CALL DAXPY(NT1AM(ISYCB),FACTSLV,WORK(KWRK1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYCB),FACTSLV,
     *           WORK(KWRK1+NT1AM(ISYCB)),1,RHO2,1)

      IF (LOCDEB) THEN
         TAL1= DDOT(NT1AM(ISYCB),RHO1,1,RHO1,1)
         TAL2= DDOT(NT2AM(ISYCB),RHO2,1,RHO2,1)
         WRITE(LUPRI,*) 'Printing RHO total. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF
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
        IF (LWRK1 .LT. 0) CALL QUIT('Insuff. work in CCMM_GTR 2')
C
        CALL DZERO(WORK(KT1AMPA),NT1AM(ISYMTC))
        CALL DZERO(WORK(KTGC),N2BST(ISYMTC))
C
C       Read one electron excitation part of response vector
C
        IOPT = 1
        CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     *                WORK(KT1AMPA),WORK(KDUM))

        IF (LOCDEB) THEN
         RHO1N = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
         RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
         WRITE(LUPRI,*)'Norm of RHO1 in CCMM_QRTRANS on input:,1', RHO1N
         WRITE(LUPRI,*)'Norm of RHO2 in CCMM_QRTRANS on input:,1', RHO2N
         RHO1N = DDOT(NT1AM(ISYMTB),WORK(KT1AMPA),1,WORK(KT1AMPA),1)
         WRITE(LUPRI,*)'Norm af B1AM in CCMM_QRTRANS on input:,1', RHO1N
        END IF

        CALL CCMM_TGB(WORK(KT1AMPA),ISYMTC,LISTC,IDLSTC,WORK(KTGC),'ET',
     *                MODEL,WORK(KWRK1),LWRK1)
C
        LABEL = 'GIVE INT'
        CALL CCLR_FA(LABEL,ISYMTC,LISTB,IDLSTB,
     &               LISTA,IDLSTA,WORK(KTGC),WORK(KWRK1),LWRK1)
C
        IF (LOCDEB) THEN
          TAL1= DDOT(NT1AM(ISYCB),WORK(KWRK1),1,WORK(KWRK1),1)
          TAL2= DDOT(NT2AM(ISYCB),WORK(KWRK1+NT1AM(ISYCB)),1,
     *              WORK(KWRK1+NT1AM(ISYCB)),1)
          WRITE(LUPRI,*) 'Printing transformed second contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
        END IF
CC
        CALL DAXPY(NT1AM(ISYCB),ONE,WORK(KWRK1),1,RHO1,1)
        CALL DAXPY(NT2AM(ISYCB),ONE,
     *             WORK(KWRK1+NT1AM(ISYCB)),1,RHO2,1)

        IF (LOCDEB) THEN
          TAL1= DDOT(NT1AM(ISYCB),RHO1,1,RHO1,1)
          TAL2= DDOT(NT2AM(ISYCB),RHO2,1,RHO2,1)
          WRITE(LUPRI,*) 'Printing RHO total. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
        END IF
C
      END IF
C     finally, we consider the effective operator T^g sigma
C     -----------------------------------------------------
C
      KTGCB = 1
      KWRK1 = KTGCB + N2BST(ISYCB)
      LWRK1 = LWORK  - KWRK1
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_GTR, 3')
C
      CALL DZERO(WORK(KTGCB),N2BST(ISYCB))
C
      CALL CCMM_TGCB(ISYRES,LISTA,LISTB,LISTC,
     *               IDLSTA,IDLSTB,IDLSTC,
     *               ISYMTA,ISYMTB,ISYMTC,WORK(KTGCB),
     *               MODEL,WORK(KWRK1),LWRK1)
C
C      TAL2= DDOT(N2BST(ISYRES),WORK(KTGCB),1,WORK(KTGCB),1)
C      WRITE(LUPRI,*) 'Printing TGAB after build. 
C     &                     Norm2 of operator: ', TAL2
      NAMPF = NT1AM(ISYCB) + NT2AM(ISYCB)
C
      KETA    = KWRK1
      KWRK2   = KETA + NAMPF
      LWRK2   = LWORK  - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_GTR, 4')
C
      CALL DZERO(WORK(KETA),NAMPF)
C
      LABEL = 'GIVE INT'
      CALL CC_ETAC(ISYCB,LABEL,WORK(KETA),
     *             LISTA,IDLSTA,0,WORK(KTGCB),WORK(KWRK2),LWRK2)
C
      KETA1   = KETA
      KETA2   = KETA1 + NT1AM(ISYCB)

      IF (LOCDEB) THEN
        TAL1= DDOT(NT1AM(ISYCB),WORK(KETA1),1,WORK(KETA1),1)
        TAL2= DDOT(NT2AM(ISYCB),WORK(KETA2),1,
     *              WORK(KETA2),1)
        WRITE(LUPRI,*) 'Printing third GTR contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF
C
C
      CALL DAXPY(NT1AM(ISYCB),ONE,WORK(KETA1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYCB),ONE,WORK(KETA2),1,RHO2,1)
C
      END
************************************************************
C  /* Deck ccsl_gtr */
      SUBROUTINE CCMM_XIETA(WORK,LWORK)
C
C-----------------------------------------------------------------
C JK, april .03
C
C Calculates the CCMM ETA and XI vectors and write them to files.
C
C------------------------------------------------------------------
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
C
      DIMENSION WORK(LWORK)
C
      CHARACTER*(8) FILMME, FILMMX, LABEL
      CHARACTER*2 LIST
      INTEGER LUMMET, LUMMXI, ISYMTR
      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
C
C---------------------
C     Init parameters.
C---------------------
C
      ISYMTR = ISYMOP
      NTAMP1  = NT1AM(ISYMTR)
      NTAMP2  = NT2AM(ISYMTR)
      IF (CCS)  NTAMP2  = 0
      NTAMP   = NTAMP1 + NTAMP2
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CCMM_XIETA
C----------------------------------------------
C
      IF (LOSPC) RETURN
C
      IF (IPRINT.GT.10) THEN
        CALL FLSHFO(LUPRI)
      ENDIF
C
      LUMMET = -1
      LUMMXI = -1
      FILMME = 'CCMM_ETA'
      FILMMX = 'CCMM__XI'
C
      CALL WOPEN2(LUMMET, FILMME, 64, 0)
      CALL WOPEN2(LUMMXI, FILMMX, 64, 0)
C
C------------------------------------
C     Dynamical allocation for CCMM :
C------------------------------------
C
      KRAAOx = 1
      KRAAOy = KRAAOx + N2BST(ISYMTR)
      KRAAOz = KRAAOy + N2BST(ISYMTR)
      KWRK1  = KRAAOz + N2BST(ISYMTR)
      LWRK1  = LWORK   - KWRK1
C
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_XIETA, 1')
C
      KXIx  = KWRK1
      KXIy  = KXIx + NTAMP
      KXIz  = KXIy + NTAMP
      KWRK2 = KXIz  + NTAMP
      LWRK2   = LWORK   - KWRK2
C
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_XIETA, 2')
C
      CALL DZERO(WORK(KRAAOx),3*N2BST(ISYMTR))
      CALL DZERO(WORK(KXIx),3*NTAMP)
C
C----------------------------
C     Read Rra in AO basis:
C----------------------------
C
      LUCCPO = -1
      CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUCCEF)
C
      LM = 0
      IADR1 = 1 
      IADR2 = 1
      LEN = NTAMP
C
      DO 700 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 701 J = NSYSBG(I), NSYSED(I)
            DO 702 K = 1, NUALIS(I)
C
              LM = LM + 1
              LABEL = 'READ INT'
C
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOx),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK2),LWRK2)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOy),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK2),LWRK2)
              CALL CCMMINT(LABEL,LUCCEF,WORK(KRAAOz),DUMMY,IRREP,
     &             ISYM,IERR,WORK(KWRK2),LWRK2)
C
              IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
                LABEL = 'GIVE INT'
                LIST  = 'L0'
                IDLINO = 1
C
                CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIx),
     *               LIST,IDLINO,0,WORK(KRAAOx),WORK(KWRK2),LWRK2)
                CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIy),
     *               LIST,IDLINO,0,WORK(KRAAOy),WORK(KWRK2),LWRK2)
                CALL CC_ETAC(ISYMTR,LABEL,WORK(KXIz),
     *               LIST,IDLINO,0,WORK(KRAAOz),WORK(KWRK2),LWRK2)
C
                CALL PUTWA2(LUMMET,FILMME,WORK(KXIx),IADR1,LEN)
                IADR1 = IADR1 + LEN
                CALL PUTWA2(LUMMET,FILMME,WORK(KXIy),IADR1,LEN) 
                IADR1 = IADR1 + LEN
                CALL PUTWA2(LUMMET,FILMME,WORK(KXIz),IADR1,LEN)
                IADR1 = IADR1 + LEN
C
C
                CALL CC_XKSI(WORK(KXIx),LABEL,ISYMTR,0,WORK(KRAAOx),
     *                       WORK(KWRK2),LWRK2)
                CALL CC_XKSI(WORK(KXIy),LABEL,ISYMTR,0,WORK(KRAAOy),
     *                       WORK(KWRK2),LWRK2)
                CALL CC_XKSI(WORK(KXIz),LABEL,ISYMTR,0,WORK(KRAAOz),
     *                       WORK(KWRK2),LWRK2)
C
                CALL PUTWA2(LUMMXI,FILMMX,WORK(KXIx),IADR2,LEN)
                IADR2 = IADR2 + LEN
                CALL PUTWA2(LUMMXI,FILMMX,WORK(KXIy),IADR2,LEN)
                IADR2 = IADR2 + LEN
                CALL PUTWA2(LUMMXI,FILMMX,WORK(KXIz),IADR2,LEN)
                IADR2 = IADR2 + LEN
C
              END IF
  702       CONTINUE
  701     CONTINUE
        END IF
  700 CONTINUE
C
      CALL GPCLOSE(LUCCEF,'KEEP')
C
      CALL WCLOSE2(LUMMET, FILMME, 'KEEP')
      CALL WCLOSE2(LUMMXI, FILMMX, 'KEEP')
C
      END
****************************************************************
C  /* deck ccsl_gtr */
      SUBROUTINE CCMM_REP1(WORK,LWORK)
C
C----------------------------------------------------------------
C JK+OC, mai .03
C
C
C-----------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "inftap.h"
#include "dummy.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "ccslvinf.h"
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      DIMENSION WORK(LWORK)
C
      CHARACTER*5 FILMM
C
      KCMO  = 1
      KEND1 = KCMO  + NLAMDS
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
        WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
        CALL QUIT('Insufficient memory for allocation in CCMM_REP1 1')
      ENDIF
C
      IF (NSYM .NE. 1) THEN
       WRITE(LUPRI,*)'ONLYMO is not implemented for other than NSYM = 1'
       CALL QUIT('Symmetry problem in CCMM_REP1')
      END IF
C
C     First we write to the INFO file
C
      NTAL = NRHFS(1)*NBAS(1)
C
      LUMMMO = -1
      CALL GPOPEN(LUMMMO,'MMMO_INFO','UNKNOWN',' ',
     &                       'FORMATTED',IDUMMY,.FALSE.)
C
      WRITE(LUMMMO,'(I10,5X,I10,5X,I10)') NTAL, NRHFS(1), NBAS(1)
      CALL GPCLOSE(LUMMMO,'KEEP')


      CALL CC_GET_CMO(WORK(KCMO))
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
      FILMM   = 'MMMOS'
      LUMMMO  = -1
      IADRFIL = 1
      NDATA   = NTAL
C
      CALL WOPEN2(LUMMMO,FILMM,64,0)
      CALL PUTWA2(LUMMMO,FILMM,WORK(KCMO),
     *                      IADRFIL,NDATA)
      CALL WCLOSE2(LUMMMO,FILMM,'KEEP')
C
      KORBE = 1
      KEND1 = KORBE + NORBTS
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
        WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
        CALL QUIT('Insufficient memory for allocation in CCMM_REP1 2')
      ENDIF
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &           .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KORBE+I-1), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      IF (IPRINT .GT. 20) THEN
            CALL AROUND('HF - Orbital energies. ')
            WRITE(LUPRI,*) (WORK(KORBE+I-1), I=1,NORBTS)
      ENDIF
C
      IF (IPRINT .GT. 5) THEN
            CALL AROUND('HF - Occupied Orbital energies. ')
            WRITE(LUPRI,*) (WORK(KORBE+I-1), I=1,NRHFS(1))
      ENDIF
C
      FILMM   = 'MMEPS'
      LUMMEP  = -1
      IADRFIL = 1
      NDATA   = NRHFS(1)
C
      CALL WOPEN2(LUMMEP,FILMM,64,0)
      CALL PUTWA2(LUMMEP,FILMM,WORK(KORBE),
     *                       IADRFIL,NDATA)
      CALL WCLOSE2(LUMMEP,FILMM,'KEEP')
C
C
      END
C**************************************************************
C  /* Deck cc_qm3nsint */
      SUBROUTINE CCMM_REP2(EREP,AODEN,WORK,LWORK)
C**************************************************************
C
C     Calculates the non-coulombic repulsive contribution to the 
C     QM/MM interaction energy.
C     JK+OC, mai 03
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccinftap.h"
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
C
      
      DIMENSION WORK(LWORK)
      DIMENSION AODEN(*)
C
      CHARACTER*9 FILMMR
C
C-----------------------------------
C
      EREP = 0.0D0
C
      IF (.NOT. (LOSHAW)) RETURN
C
      KINT  = 1
      KWRK1 = KINT + N2BST(ISYMOP)
      LWRK1 = LWORK - KWRK1
C
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_REP2 1' )
C
      CALL DZERO(WORK(KINT),N2BST(ISYMOP))
C
      FILMMR  = 'MMREPOP_0'
      LUMMRE  = -1
      IADRFIL = 1
      NDATA   = N2BST(ISYMOP)
C
      CALL WOPEN2(LUMMRE,FILMMR,64,0)
C
      CALL GETWA2(LUMMRE,FILMMR,WORK(KINT),IADRFIL,NDATA)
C
      IF (REPTST) THEN
        NBASLO = MAX(NBAST,1)
C
        KINT1 = KWRK1
        KINT2 = KINT1 + N2BST(ISYMOP)
        KWRK2 = KINT2 + N2BST(ISYMOP)
        LWRK2 = LWORK - KWRK2
C
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT( 'Too litle work in CCMM_REP2 Rep. 2')
        END IF
C
        CALL DZERO(WORK(KINT1),2*N2BST(ISYMOP))
        CALL DCOPY(N2BST(ISYMOP),WORK(KINT),1,WORK(KINT1),1)
C
        DO 668 KR=1,NREPMT
          CALL DGEMM('N','N',NBAST,NBAST,NBAST,
     *               ONE,WORK(KINT),NBASLO,
     *               WORK(KINT1),NBASLO,
     *               ZERO,WORK(KINT2),NBASLO)
C
          CALL DCOPY(N2BST(ISYMOP),WORK(KINT2),1,WORK(KINT),1)
  668   CONTINUE
      END IF
C
      IF (IQM3PR .GT.14) THEN
        TAL1 = DDOT(N2BST(ISYMOP),WORK(KINT),1,WORK(KINT),1)
        TAL2 = SQRT(TAL1)
      END IF
C
      EREP = DDOT(N2BST(ISYMOP),WORK(KINT),1,AODEN,1)
C
      CALL WCLOSE2(LUMMRE,FILMMR, 'KEEP')
C
      END
