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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
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
#include <qm3.h>
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
         
      CALL CC_QM3E(ECMCON,EPOL,EEL,EELEL,EREP,AODEN,
     &             WORK,LWORK)
C
C----------------------------
C     Calculate new E(QM/MM):
C----------------------------
C
      IF (.NOT. (OLDTG)) THEN
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
        WRITE(LUPRI,'(//1X,72A)') ('-',J=1,72)
        WRITE(LUPRI,'(1X,A56,A16)')
     * '|<<<<<<<<<<<<<<< The MM/MM classical interaction energy ',
     * '>>>>>>>>>>>>>>>|'
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
C
C--------------------------------------------
C  Second the CC/MM interaction energy terms:
C--------------------------------------------
C
        WRITE(LUPRI,'(1X,A1,70A,A1)') '+',('=',J=1,70),'+'
C
        IF (MODEL .EQ. 'CC2')THEN
          WRITE(LUPRI,'(1X,A1,A15,A18,A3,A19,A15,A1)')
     *         '|','<<<<<<<<<<<<<< ',
     *         'Final output from ',MODEL,'/MM energy program ',
     *         '>>>>>>>>>>>>>>>','|'
        ELSE IF (MODEL .EQ. 'CCSD') THEN
          WRITE(LUPRI,'(1X,A1,A15,A18,A4,A19,A14,A1)')
     *         '|','<<<<<<<<<<<<<< ',
     *         'Final output from ',MODEL,'/MM energy program ',
     *         '>>>>>>>>>>>>>>','|'
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
     *  MODELPRI,' Total energy:           ',ECMCU
C
        WRITE(LURES,'(12X,A4,A,F20.10,/)')
     *  MODELPRI,' E(QM/MM):       ',ECMCON
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
C     Purpose: Calculation of E(QM/MM). 
C     QM3, JK+AO+OC+KMI,01
C     Modified by KS, nov 09
C     Purpose: Updates the induced dipoles besides calculating E(QM/MM)
C
C**************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxorb.h>
#include <orgcom.h>
#include <nuclei.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <exeinf.h>
#include <ccfdgeo.h>
#include <thrzer.h>

C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 1.0D0/2.0D0)
C
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB, FNDLAB, EXCENT
      DIMENSION WORK(LWORK),AODEN(*)
      CHARACTER LABEL*8
      INTEGER NTOTI

      DIMENSION CCNS(MXQM3),ALPHAMM(MXQM3),FFJ(3)
      INTEGER KS
C
      CALL QENTER('CC_QM3E')
C this assignment is temporary
C------------------------------------------------------------------
C       Init parameters
C------------------------------------------------------------------

      LOCDEB=.FALSE.
      NDIM = 3*NNZAL
      ISYMOP = 1
C
C------------------------------------------------------------------
C       Dynamical allocation for CCMM 
C------------------------------------------------------------------
C

      IF (.NOT. NYQMMM) THEN
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
     *    '|<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
     *     //'The different energy contributions outlined'
     *     //' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|'
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
      
      ELSE IF (NYQMMM) THEN

         KRAAO  = 1
         KQMELF = KRAAO +3*NNBASX
         KELF = KQMELF + NDIM
         KRAAOx = KELF + NDIM
         KRAAOy = KRAAOx + N2BST(ISYMOP)
         KRAAOz = KRAAOy + N2BST(ISYMOP)
         KINDIP = KRAAOz + N2BST(ISYMOP)
         KRELAY = KINDIP + NDIM
         KMULTI = KRELAY + NDIM*NDIM
         KNSAO = KMULTI + NNBASX
         KINT = KNSAO + NNBASX
         KWRK1 = KINT + N2BST(ISYMOP)
         LWRK1 = LWORK - KWRK1
C
         CALL DZERO(WORK(KRAAO),3*NNBASX)
         CALL DZERO(WORK(KQMELF),NDIM)
         CALL DZERO(WORK(KELF),NDIM)
         CALL DZERO(WORK(KRAAOx),3*N2BST(ISYMOP))
         CALL DZERO(WORK(KINDIP),NDIM)
         CALL DZERO(WORK(KRELAY),NDIM*NDIM)
C
C
         IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CC_QM3E')
C 
C sne implement for finite field? 
C     First we see if any static fields are to be added
C
C        FFJ(1) = 0.0D0
C        FFJ(2) = 0.0D0
C        FFJ(3) = 0.0D0
C
C sne 
C      IF (RELMOM) THEN
C        DO 330 IF =1, NFIELD
C          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
C          IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
C          IF (LFIELD(IF) .EQ. 'ZDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
C  330   CONTINUE
C      END IF
C
C First we find the field due to the distributed multipoles
         IF (IPOLTP .GT. 0) THEN
            LRI = 0 ! Row-index in large matrix

            DO 100 I = 1,MMCENT

               IF(ZEROAL(I).EQ. -1) THEN
                  GOTO 100
               END IF
 
               LCI = 0 ! Column index in the large matrix

               DO 200 J=1,MMCENT

                  IF (J .NE. I) THEN     
                     DIST2 = (MMCORD(1,I)-MMCORD(1,J))**2 +
     &                      (MMCORD(2,I)-MMCORD(2,J))**2 +
     &                      (MMCORD(3,I)-MMCORD(3,J))**2
                     DIST = SQRT(DIST2)

                     IF(DIST .GT. RCUTMM) THEN
                        LCI = LCI +3
                        IF (LOCDEB) THEN
                           WRITE(LUPRI,*) 'Skipping element 
     &                                     in T2 CC', I,J
                        ENDIF
                        GOTO 200
                     ENDIF
               
                     EXCENT=.FALSE.
               
                     DO 250 L=1,NEXLST
                        IF(EXLIST(L,I).EQ.EXLIST(1,J)) EXCENT=.TRUE.
 250                 CONTINUE

C Form F vector due to the permament moments

C              A) The point-charge contribution
                     IF (.NOT. (EXCENT) ) THEN
                        IF((ABS(MUL0MM(J)) .GT. THRMM) 
     &                             .AND. (NMULT.GE.0) ) THEN
        
                           ELFLDX = 0.0D0 
                           ELFLDY = 0.0D0
                           ELFLDZ = 0.0D0 

                           CALL GET_CHARGE_ELFLD(MUL0MM(J),
     *                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     *                          ELFLDX,ELFLDY,ELFLDZ)
                           WORK(KELF+LRI+0) = WORK(KELF+LRI+0)+ELFLDX
                           WORK(KELF+LRI+1) = WORK(KELF+LRI+1)+ELFLDY
                           WORK(KELF+LRI+2) = WORK(KELF+LRI+2)+ELFLDZ

                        ENDIF
            
C              B) The dipole contribution
                        IF (NMULT .GE. 1) THEN
                           ELFLDX = 0.0D0
                           ELFLDY = 0.0D0
                           ELFLDZ = 0.0D0
                           CALL GET_DIPOLE_ELFLD(MUL1MM(1,J),
     *                              MUL1MM(2,J),MUL1MM(3,J),
     *                              MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                              MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     *                              ELFLDX,ELFLDY,ELFLDZ)
                           WORK(KELF+LRI+0) = WORK(KELF+LRI+0)+ELFLDX
                           WORK(KELF+LRI+1) = WORK(KELF+LRI+1)+ELFLDY
                           WORK(KELF+LRI+2) = WORK(KELF+LRI+2)+ELFLDZ
                        ENDIF

C              C) The quadrupole contribution
                        IF (NMULT .GE. 2) THEN
                           ELFLDX = 0.0D0
                           ELFLDY = 0.0D0
                           ELFLDZ = 0.0D0
                           CALL GET_QUADRUPOLE_ELFLD(
     *                             MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     *                             MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     *                             MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                             MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     *                             ELFLDX,ELFLDY,ELFLDZ)
                           WORK(KELF+LRI+0) = WORK(KELF+LRI+0)+ELFLDX
                           WORK(KELF+LRI+1) = WORK(KELF+LRI+1)+ELFLDY
                           WORK(KELF+LRI+2) = WORK(KELF+LRI+2)+ELFLDZ
                        ENDIF 
            
C              D) The octopole contribution
                        IF (NMULT .GE. 3) THEN
                           ELFLDX = 0.0D0
                           ELFLDY = 0.0D0
                           ELFLDZ = 0.0D0
                           CALL GET_OCTUPOLE_ELFLD(
     *                          MUL3MM(1,J),MUL3MM(2,J),MUL3MM(3,J),
     *                          MUL3MM(4,J),MUL3MM(5,J),MUL3MM(6,J),
     *                          MUL3MM(7,J),MUL3MM(8,J),MUL3MM(9,J),
     *                          MUL3MM(10,J),
     *                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     *                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     *                          ELFLDX,ELFLDY,ELFLDZ)
                           WORK(KELF+LRI+0) = WORK(KELF+LRI+0)+ELFLDX
                           WORK(KELF+LRI+1) = WORK(KELF+LRI+1)+ELFLDY
                           WORK(KELF+LRI+2) = WORK(KELF+LRI+2)+ELFLDZ
                        ENDIF
                     ENDIF
     
                     IF (ZEROAL(J) .EQ. -1) GOTO 200

                     LCI = LCI + 3
                  END IF
 200           CONTINUE

C Then we calculate the field due to the QM electrons
               KPATOM = 0
               NOSIM  = 3
               TOFILE = .FALSE.
               TRIMAT = .TRUE.
               EXP1VL = .FALSE. 
      
               DIPORG(1) = MMCORD(1,I)
               DIPORG(2) = MMCORD(2,I)
               DIPORG(3) = MMCORD(3,I) 

               DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *               (MMCORD(2,I)-QMCOM(2))**2 + 
     *               (MMCORD(3,I)-QMCOM(3))**2
               DIST = SQRT(DIST2)
               IF (DIST .GT. RCUTMM) THEN
                  LRI = LRI + 3
                  IF (LOCDEB) THEN
                     WRITE(LUPRI,*) 'Skipping site in F vector ', I
                  ENDIF
                  GOTO 100
               ENDIF

               CALL DZERO(WORK(KRAAO),3*NNBASX)
               RUNQM3 = .TRUE.
               CALL GET1IN(WORK(KRAAO),'NEFIELD',NOSIM,WORK(KWRK1),
     *                  LWRK1,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     *                  TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)
               RUNQM3 = .FALSE.

               IF ( (IPQMMM.GT.15) .OR. (LOCDEB) ) THEN
                  WRITE(LUPRI,'(/A)') 'Packed E_x matrix in CC_QM3E:'
                  CALL OUTPAK(WORK(KRAAO),NBAST,1,LUPRI)

                  WRITE(LUPRI,'(/A)') 'Packed E_y matrix in CC_QM3E:'
                  CALL OUTPAK(WORK(KRAAO+NNBASX),NBAST,1,LUPRI)

                  WRITE(LUPRI,'(/A)') 'Packed  E_z matrix in CC_QM3E:'
                  CALL OUTPAK(WORK(KRAAO+2*NNBASX),NBAST,1,LUPRI)
               ENDIF

C Unpack the integrals
               CALL DZERO(WORK(KRAAOx),3*N2BST(ISYMOP))
               LABEL = 'GIVE INT'
               CALL CCMMINT(LABEL,LUQMMM,WORK(KRAAOx),WORK(KRAAO),
     *                     IRREP,ISYM,IERR,WORK(KWRK1),LWRK1)
               CALL CCMMINT(LABEL,LUQMMM,WORK(KRAAOy),
     *                      WORK(KRAAO+NNBASX),IRREP,ISYM,IERR,
     *                      WORK(KWRK1),LWRK1)
               CALL CCMMINT(LABEL,LUQMMM,WORK(KRAAOz),
     *                      WORK(KRAAO+2*NNBASX),IRREP,ISYM,IERR,
     *                      WORK(KWRK1),LWRK1)
C Dot with the AO density

               EXELCO = DDOT(N2BST(ISYMOP),AODEN,1,WORK(KRAAOx),1)
               EYELCO = DDOT(N2BST(ISYMOP),AODEN,1,WORK(KRAAOy),1)
               EZELCO = DDOT(N2BST(ISYMOP),AODEN,1,WORK(KRAAOz),1)

               IF (LOCDEB) THEN
                  WRITE(LUPRI,*) 'QM,el field: CC',EXELCO,EYELCO,EZELCO
               END IF
 
               WORK(KQMELF+LRI+0) = WORK(KQMELF+LRI+0)+EXELCO
               WORK(KQMELF+LRI+1) = WORK(KQMELF+LRI+1)+EYELCO
               WORK(KQMELF+LRI+2) = WORK(KQMELF+LRI+2)+EZELCO
C Add to the electric field from the multipoles 
               WORK(KELF+LRI+0) = WORK(KELF+LRI+0) + WORK(KQMELF+LRI+0)
               WORK(KELF+LRI+1) = WORK(KELF+LRI+1) + WORK(KQMELF+LRI+1)
               WORK(KELF+LRI+2) = WORK(KELF+LRI+2) + WORK(KQMELF+LRI+2)

C Now the field from the QM nuclei

               TJKX = 0.0D0
               TJKY = 0.0D0
               TJKZ = 0.0D0

               DO 145 J=1,NUCIND ! Dalton keyword for nr of QM nuclei
                  CALL GET_CHARGE_ELFLD(CHARGE(J),
     *                   CORD(1,J),CORD(2,J),CORD(3,J),
     *                   MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     *                   ELFLDX,ELFLDY,ELFLDZ)

                  WORK(KELF+LRI+0)=WORK(KELF+LRI+0)+ELFLDX
                  WORK(KELF+LRI+1)=WORK(KELF+LRI+1)+ELFLDY
                  WORK(KELF+LRI+2)=WORK(KELF+LRI+2)+ELFLDZ
            
                  TJKX = TJKX + ELFLDX
                  TJKY = TJKY + ELFLDY
                  TJKZ = TJKZ + ELFLDZ
 145           CONTINUE

               IF (LOCDEB) THEN
                  WRITE(LUPRI,*) 'Nuclear field: CC ', TJKX,TJKY,TJKZ
               ENDIF
               LRI = LRI + 3

 100        CONTINUE

C     Read the relay matric from file and construct the induced moments


            LUQMMM = -1
            IF(LUQMMM .LT. 0) THEN
               CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &                       'UNFORMATTED',IDUMMY,.FALSE.)
            ENDIF
            REWIND(LUQMMM)
      
            IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
               WRITE(LUPRI,*) 'reading relay of length: ', NDIM*NDIM
               CALL READT(LUQMMM,NDIM*NDIM,WORK(KRELAY))
            ELSE
               CALL QUIT('Problem reading the relay 
     &                    matrix from QMMMIM file')
            ENDIF

            CALL GPCLOSE(LUQMMM,'KEEP')

            IF (IPQMMM .GE. 5) THEN
               WRITE(LUPRI,*) 'The relay matrix is read from file'
               CALL OUTPUT(WORK(KRELAY),1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
            ENDIF

            IF (IPRINT .GT. 1) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,1051)
               WRITE(LUPRI,1050)
               WRITE(LUPRI,1051)
               WRITE(LUPRI,*)
            ENDIF

            IF(LOCDEB) THEN
               WRITE(LUPRI,*) 'F vector: '
                DO I=1,NDIM
                   WRITE(LUPRI,*) WORK(KELF+I-1)
                END DO
            ENDIF

C calculate the corresponding induced moments
            NTOTI = MAX(NDIM,1)
            CALL DGEMV('N',NDIM,NDIM,ONE,WORK(KRELAY),NTOTI,
     *                 WORK(KELF),1,ZERO,WORK(KINDIP),1)
 
            IF (IPQMMM .GT. 1) THEN
               WRITE(LUPRI,*)
               WRITE(LUPRI,1030)
               WRITE(LUPRI,*)
               WRITE(LUPRI,1000)
               WRITE(LUPRI,1010)
               WRITE(LUPRI,1000)
            ENDIF

            IINIM = 0

            DO I=1,MMCENT
               IF (ZEROAL(I) .EQ. -1) THEN
                  DIPX = 0.0D0
                  DIPY = 0.0D0
                  DIPZ = 0.0D0
               ELSE
                  DIPX = WORK(KINDIP + IINIM )
                  DIPY = WORK(KINDIP + IINIM +1)
                  DIPZ = WORK(KINDIP + IINIM +2)
                  IINIM = IINIM +3
               ENDIF
               IF (IPQMMM .GT. 1) THEN
                  WRITE(LUPRI,1020) I, DIPX, DIPY, DIPZ
               END IF
            END DO

            IF(IPRINT .GT. 1) THEN
               WRITE(LUPRI,1000)
               WRITE(LUPRI,*)
            ENDIF

C Finally write the updated induced moments to file 
            CALL PUT_TO_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDIP))

         END IF
C Finally construct the effective Tg operator and subsequently calculate
C the polarization energy 
C-------------------------------------
C     Add up the energy contributions:
C     1) Epol:
C-------------------------------------
         ECCMM = ZERO
         EPOL = ZERO
         ETEMP = ZERO

C calculate the polarization energy as -1/2 sum_a mu_a^ind*F^sta
         ETEMP = DDOT(NDIM,WORK(KINDIP),1,WORK(KELF),1)
      
C Remember the factor of -1/2
         EPOL = -HALF*ETEMP

C Now add the electrostatic energy contribution
C we do this by getting the multipole integrals from file
C and taking the expectation value. The integrals are at our disposal
C from the preceeding HF/MM calculation handled by sirius


C--------
C 2) Eel:
C--------
C
C the charge correction
         LABEL = 'GIVE INT'
         CALL DZERO(WORK(KNSAO),NNBASX)
         CALL DZERO(WORK(KMULTI),NNBASX)
         IF (NMULT .GE. 0) THEN

            LUQMMM = -1
            CALL GPOPEN(LUQMMM,'MU0INT','OLD',' ',
     &           'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND(LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KMULTI))

C sne: Minus carried in the integrals cf sirqm3
            CALL DAXPY(NNBASX,ONE,WORK(KMULTI),1,WORK(KNSAO),1)
            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF
     

C the dipole correction

         IF (NMULT .GE. 1) THEN
            CALL DZERO(WORK(KMULTI),NNBASX)
            CALL GPOPEN(LUQMMM,'MU1INT','OLD',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND(LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KMULTI))

C sne: Minus correct in daxpy cf sirqm3
            CALL DAXPY(NNBASX,-ONE,WORK(KMULTI),1,WORK(KNSAO),1)
            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF

C the quadrupole correction

         IF (NMULT .GE. 2) THEN
            CALL DZERO(WORK(KMULTI),NNBASX)
            CALL GPOPEN(LUQMMM,'MU2INT','OLD',' ',
     &                  'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND(LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KMULTI))

C sne: Minus carried in the integrals cf sirqm3
            CALL DAXPY(NNBASX,ONE,WORK(KMULTI),1,WORK(KNSAO),1)
            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF
      
C the octopole correction
         IF (NMULT .GE. 3) THEN
            CALL DZERO(WORK(KMULTI),NNBASX)
            CALL GPOPEN(LUQMMM,'MU3INT','OLD',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND(LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KMULTI))

C sne: Minus carried in the integrals cf sirqm3
            CALL DAXPY(NNBASX,ONE,WORK(KMULTI),1,WORK(KNSAO),1)
            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF

C Unpack the multipole integrals and dot with the density to
C obtain the expectation value

         LABEL = 'GIVE INT'
         CALL CCMMINT(LABEL,LUQMMM,WORK(KINT),WORK(KNSAO),
     *               IRREP,ISYM,IERR,WORK(KWRK1),LWRK1)
 
         EEL = DDOT(N2BST(ISYMOP),WORK(KINT),1,AODEN,1)

C Finally we add the energy between QM nuclei and MM multipoles. This 
C was found in the previous SCF/MM calculation and is readily available 
C in the COMMON block. 
         EEL = EEL + ENUMUL

         ECCMM = EPOL + EEL
         IF (LOCDEB .OR. IPQMMM .GT. 4) THEN
            WRITE(LUPRI,*) 'EEL: ', EEL
            WRITE(LUPRI,*) 'EPOL: ', EPOL
            WRITE(LUPRI,*) 'ECCMM: ', ECCMM
            WRITE(LUPRI,*) 'The nuclear contribution: ', ENUMUL
         END IF
C
C
C-------------
C 3) Edisp
C-------------
C
C sne - New qmmm not implemented for dispersion yet
C      ECCMM = ECCMM + ECLVDW
C
C-----------------------------------------------
C sne New qmmm not implemented for explicit repulsion yet
C 4) E_repulsion
C
C      CALL CCMM_REP2(EREP,AODEN,WORK(KWRK1),LWRK1)
C      ECCMM = ECCMM + EREP
C
C-----------------------------------------------
C   Writing out the final energy contributions one by one:
C-----------------------------------------------
C
 1050 FORMAT('  Induced dipole moments   ')
 1051 FORMAT(2X,'=',22('-'),'=',2X)
 1000 FORMAT(1X,51('='))
 1020 FORMAT(1X,I6,3(4X,F10.6))
 1010 FORMAT(' | Site |     X     |     Y     |     Z     |')
 1030 FORMAT(' Total induced dipole moments: ')

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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
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
#include <qm3.h>
#include <ccslvinf.h>
#include <ccinftap.h>
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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
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
#include <ccinftap.h>
#include <qm3.h>
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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
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
#include <qm3.h>
#include <ccinftap.h>
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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
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
#include <qm3.h>
#include <ccinftap.h>
#include <ccslvinf.h>
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
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>
#include <inftap.h>
#include <ccisao.h>
#include <dummy.h>
#include <qm3.h>
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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
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
#include <ccinftap.h>
#include <qm3.h>
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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
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
#include <ccinftap.h>
#include <qm3.h>
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
#include <implicit.h>
#include <dummy.h>
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
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
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
#include <qm3.h>
#include <ccinftap.h>
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
#include <qm3.h>
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
      CHARACTER MODEL*10
C
      CHARACTER*8 LABEL,LR*(2),LIST
      CHARACTER*(8) FILMME, FILMMX
      LOGICAL LEXIST
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
C
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
C----------------------------
C     Read Rra in AO basis:
C----------------------------
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
C             (LR.EQ.'ET') Calculate eta^T_lm
C             (LR.EQ.'XI') Calculate xi^T_lm
C             Only if MDLWRD(I) = 'SPC_E01' are we going to calculate
C             contributions from the T^g operator to the response
C             equations.
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
C               Contract with trial vector:
C               (If LR.EQ.'ET' then from the right hand side)
C               (If LR.EQ.'XI' then from the left hand side)
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

C--------------------------------------------------------------
C        Put factor on integrals 
C--------------------------------------------------------------
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
#include <qm3.h>
#include <ccinftap.h>
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
#include <qm3.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*)
C
      CHARACTER*(*) LISTA,LISTB
      CHARACTER*(3) LISTC
      CHARACTER*8 LABEL
      CHARACTER MODEL*10
      LOGICAL LEXIST, LSAME
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
      CALL CCMM_TGCB(ISYMAB,LISTC,LISTA,LISTB,
     *               IDLSTC,IDLSTA,IDLSTB,
     *               ISYMC,ISYMA,ISYMB,WORK(KTGAB),
     *               MODEL,WORK(KWRK1),LWRK1)
C
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
#include <qm3.h>
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
#include <qm3.h>
#include <ccinftap.h>
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*)
C
      CHARACTER*(*) LISTA, LISTB, LISTC
      CHARACTER*8 LABEL
      CHARACTER MODEL*10
      LOGICAL LEXIST, LSAME
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
      CALL DAXPY(NT1AM(ISYCB),FACTSLV,WORK(KWRK1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYCB),FACTSLV,
     *           WORK(KWRK1+NT1AM(ISYCB)),1,RHO2,1)
C
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

        CALL CCMM_TGB(WORK(KT1AMPA),ISYMTC,LISTC,IDLSTC,WORK(KTGC),'ET',
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
#include <qm3.h>
#include <ccinftap.h>
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
************************************************************
C  /* Deck quitmm */
      SUBROUTINE QUITMM(TEXT)
#include <implicit.h>
#include <priunit.h>
      CHARACTER TEXT*(*)
C
C     Stamp date and time and hostname to output
C
      CALL TSTAMP('  --- PROGRAM STOPS ---',
     &   LUPRI)
      WRITE (LUPRI,'(/2A/)') ' Reason: ',TEXT
      CALL GETTIM(CTOT,WTOT)
      CALL TIMTXT('>>>> Total CPU  time used in DALTON:',CTOT,LUPRI)
      CALL TIMTXT('>>>> Total wall time used in DALTON:',WTOT,LUPRI)
C
      STOP 100
      END
****************************************************************
C
C  /* deck ccsl_gtr */
      SUBROUTINE CCMM_REP1(WORK,LWORK)
C
C----------------------------------------------------------------
C JK+OC, mai .03
C
C
C-----------------------------------------------------------------
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <dummy.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccslvinf.h>
#include <qm3.h>
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
C******************************************************************
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
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
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
#include <ccinftap.h>
#include <qm3.h>
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
C***********************************************************************
C  /* Deck ccmm_tgt */
      SUBROUTINE CCMM_TGT(FOCK,WORK,LWORK)
C**************************************************************
C
C     Purpose: Update t parameters
C     In the new QMMM implementation it replaces CCMM_RHSTG
C
C     sneskov oct 09
C**************************************************************
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <qmmm.h>
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
#include <qm3.h>
#include <ccinftap.h>
#include <nuclei.h>
#include <orgcom.h>
C
      PARAMETER (ZERO = 0.0D0, MONE=-1.0D0 ,ONE = 1.0D0,
     &             IZERO = 0 , TWO = 2.0D0)
      DIMENSION WORK(LWORK),FOCK(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB
C
      CHARACTER*8 LABEL
C
C
C-------------------------
C       Init parameters
C-------------------------
      NDIM=3*NNZAL
C--------------------------
C       Dynamical allocation
C--------------------------
C
      KTAO  =  1 
      KRAAO = KTAO   + NNBASX
      KINDIP = KRAAO + 3*NNBASX
      KNSAO = KINDIP + NDIM
      KINT = KNSAO + NNBASX
      KWRK1 = KINT + N2BST(ISYMOP)
      LWRK1   = LWORK   - KWRK1

      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_TGT' )
        
      CALL DZERO(WORK(KTAO),NNBASX)
      CALL DZERO(WORK(KRAAO),3*NNBASX)
      CALL DZERO(WORK(KINDIP),NDIM)
      CALL DZERO(WORK(KNSAO),NNBASX)
      CALL DZERO(WORK(KINT),N2BST(ISYMOP))
C
C Backup dipole origin
      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)
C
C     First we construct the polarization contribution to Tg.
C
C     The first time we enter this routine we read the converged 
C     induced dipoles from the preceeding HF/MM calculation.
C     In later iterations we read the induced dipoles written in 
C     CC_QM3E
      
     
      IF (IPOLTP .GT. 0) THEN 
         WRITE(LUPRI,*) 'Reading the induced dipoles from file.'
      
         CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDIP))
C
C-------------------
C     Print section:
C-------------------
C
      
         IF (IPQMMM .GT. 10) THEN
            WRITE(LUPRI,'(//,A)')
     *        ' +============================================'
     *       //'==============+'
           WRITE(LUPRI,'(A)')
     *        ' | MMcent | <mu^ind>_x  | <mu^ind>_y      |'
     *       //' <mu^ind>_z      |'
           WRITE(LUPRI,'(1X,A)')
     *        '+============================================'
     *       //'==============+'
C
           LRI = 0
C
           DO 100 I = 1, MMCENT

              IF (ZEROAL(I) .EQ. -1) GOTO 100
C
C
              WRITE(LUPRI,'(1X,A,I3,A,F16.10,A,
     &                         F16.10,A,F16.10,A)')
     *           '| ', I,'|', WORK(KINDIP + LRI+0 ),' |',
     *              WORK(KINDIP + LRI+1 ),' |', 
     *              WORK(KINDIP + LRI+2 ),' |'
                    WRITE(LUPRI,'(1X,A)')
     *           '+---------------------------------------------'
     *           //'-------------+'
              LRI = LRI + 3
 100       CONTINUE 
           WRITE(LUPRI,'(//,A)')
         END IF
C
C
C---------------------------------------
C Calculate Rr(s) in ao basis and add to tao:
C Note we are doing everything integral-direct to avoid
C a possible storage-problem when working with large systems.
C---------------------------------------
C
C
         LRI = 0
C
         DO 200 I = 1, MMCENT

            IF (ZEROAL(I) .EQ. -1) GOTO 200

            DIPORG(1) = MMCORD(1,I)
            DIPORG(2) = MMCORD(2,I)
            DIPORG(3) = MMCORD(3,I)
 
            DIST2 = (MMCORD(1,I)-QMCOM(1))**2 + 
     &            (MMCORD(2,I)-QMCOM(2))**2 + 
     &            (MMCORD(3,I)-QMCOM(3))**2  
            DIST = SQRT(DIST2)

            IF (DIST .GT. RCUTMM) THEN
               LRI = LRI+3
               IF (LOCDEB) THEN
                  WRITE(LUPRI,*) 'Skipping site in mu^ind vector 
     &                        in CCMM_TGT', I
               ENDIF
               GOTO 200
            ENDIF
      
            KPATOM = 0
            NOCOMP = 3
            TOFILE = .FALSE.
            TRIMAT = .TRUE.
            EXP1VL = .FALSE.   

            RUNQM3 = .TRUE.
            CALL DZERO(WORK(KRAAO),3*NNBASX)
            CALL GET1IN(WORK(KRAAO),'NEFIELD',NOCOMP,WORK(KWRK1),
     *              LWRK1,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     *              TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)
            RUNQM3 = .FALSE.

            FACx =  - (WORK(KINDIP + LRI +0 ))  
            FACy =  - (WORK(KINDIP + LRI +1 ))  
            FACz =  - (WORK(KINDIP + LRI +2 ))  
C
            CALL DAXPY(NNBASX,FACx,WORK(KRAAO),
     *                      1,WORK(KTAO),1)
C
            CALL DAXPY(NNBASX,FACy,WORK(KRAAO+NNBASX),
     *                        1,WORK(KTAO),1)
C 
            CALL DAXPY(NNBASX,FACz,WORK(KRAAO+2*NNBASX),
     *                     1,WORK(KTAO),1)

            LRI = LRI + 3
 200     CONTINUE
C
         IF (IPQMMM.GT.14) THEN
            WRITE (LUPRI,*) ' NORM of TAO matrix after pol. contr.: ',
     *      DDOT(NNBASX,WORK(KTAO),1,WORK(KTAO),1)
            WRITE (LUPRI,'(/A)') ' TAO matrix: '
            CALL OUTPAK(WORK(KTAO),NBAST,1,LUPRI)
         ENDIF
      END IF
C
C--------------------------------------
C     Now the electrostatic contribution to Tg. 
C     We get the multipole integrals from the HF/MM
C     calculation available from file.
C--------------------------------------
C
C We do all the corrections simultaneosly. We will analyze
C the HF/MM results if we want to consider each multipole contribution
C separately.

C 1) The charge contribution
C
C
      IF (NMULT .GE. 0) THEN
C
         LUQMMM = -1
         CALL GPOPEN(LUQMMM,'MU0INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQMMM)
         CALL READT(LUQMMM,NNBASX,WORK(KNSAO))
         CALL GPCLOSE(LUQMMM,'KEEP')
         CALL OUTPAK(WORK(KNSAO),NBAST,1,LUPRI)

C sne the minus is in the integrals so add to t
         CALL DAXPY(NNBASX,ONE,WORK(KNSAO),
     *                   1,WORK(KTAO),1)

      END IF
C 2) The dipole contribution
      IF (NMULT .GE. 1) THEN
         CALL DZERO(WORK(KNSAO),NNBASX)
         CALL GPOPEN(LUQMMM,'MU1INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQMMM)
         CALL READT(LUQMMM,NNBASX,WORK(KNSAO))

C sne: Keep minus in daxpy cf sirm3
         CALL DAXPY(NNBASX,-ONE,WORK(KNSAO),
     *                   1,WORK(KTAO),1)

         CALL GPCLOSE(LUQMMM,'KEEP')
      END IF
C 3) The quadrupole contribution
      IF (NMULT .GE. 2) THEN
         CALL DZERO(WORK(KNSAO),NNBASX)
         CALL GPOPEN(LUQMMM,'MU2INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQMMM)
         CALL READT(LUQMMM,NNBASX,WORK(KNSAO))

C Sne: the minus is carried in the integrals 
         CALL DAXPY(NNBASX,ONE,WORK(KNSAO),
     *                   1,WORK(KTAO),1)

         CALL GPCLOSE(LUQMMM,'KEEP')
      END IF
C 4) The octopole contribution
      IF (NMULT .GE. 3) THEN
         CALL DZERO(WORK(KNSAO),NNBASX)
         CALL GPOPEN(LUQMMM,'MU3INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQMMM)
         CALL READT(LUQMMM,NNBASX,WORK(KNSAO))

C sne: the minus is carried in the integrals
         CALL DAXPY(NNBASX,ONE,WORK(KNSAO),
     *                   1,WORK(KTAO),1)

         CALL GPCLOSE(LUQMMM,'KEEP')
      END IF
C
C--------------------------------------------------
C     Add total contribution to effective AO fock matrix:
C--------------------------------------------------
C
      IF (IPQMMM .GT. 14) THEN
        WRITE(LUPRI,*) 'NORM of TAO matrix after electrostatic cont: ',
     *  DDOT(NNBASX,WORK(KTAO),1,WORK(KTAO),1)
        WRITE (LUPRI,*) ' TAO matrix: '
        CALL OUTPAK(WORK(KTAO),NBAST,1,LUPRI)
      ENDIF
C
C
C--------------------------
C     Dynamical allocation:
C--------------------------


      LABEL= 'GIVE INT'
      WRITE(LUPRI,*) 'Print IRREP, ISYM, IERR', IRREP, ISYM, IERR
      CALL CCMMINT(LABEL,-1,WORK(KINT),WORK(KTAO),
     &             IRREP,ISYM,IERR,WORK(KWRK1),LWRK1)
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'CCMM cont. to AO matrix' )
         CALL OUTPUT(WORK(KINT),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
      CALL DAXPY(N2BST(ISYMOP),ONE,WORK(KINT),1,FOCK,1)
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Modified Fock AO matrix when leaving TGT')
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
C Restore dipole origin
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      END
C
C***********************************************************************
C
C  /* Deck ccmm_tgtbar */
      SUBROUTINE CCMM_TGTBAR(RHO1,RHO2,CTR1,CTR2,MODEL,
     *                          ISYMTR,LR,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C   Construct effective operators from contracted densities
C   Find the polarization part of Tg by applying the relay matrix
C
C   IF (LR.EQ.'R').OR.(LR.EQ.'F') then the following effective 
C   operator is constructed 
C   T^gC = - Sum_a alpha_a*(tr(Rr_a*D^C)Rr_a)
C   Note we denote this operator G^C in the general notation
C
C   IF (LR.EQ.'L').OR.(LR.EQ.'P') then the following effective operator 
C   is constructed
C   ^BT^g = - Sum_a * alpha_a Sum_sigma (tr(Rr_a*D^B)Rr_a)
C   Note we also denote this operator G^B in the general notation
C
C   Finally, if LR.EQ.0 then construct
C
C   Kristian Sneskov, oct 09 - alternative to CCMM_LTRB routine chosen
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <qmmm.h>
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
#include <qm3.h>
#include <ccinftap.h>
#include <nuclei.h>
#include <orgcom.h>

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
      INTEGER NDIM, NTOTI 
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*),CTR1(*),CTR2(*)
      DIMENSION ALPHAMM(MXQM3)
C
      CHARACTER MODEL*10 
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB, FNDLAB
C
      CHARACTER*8 LABEL,LR*(1), LIST*(2)
C
C
      IF (IPQMMM.GT.10) THEN
        WRITE(LUPRI,*)'CCMM_TGTBAR : ISYMTR:', ISYMTR,
     &               ' and LR: ', LR
      ENDIF
C
C
C---------------------
C     Init parameters.
C---------------------
C
      NTAMP1  = NT1AM(ISYMTR)
      NTAMP2  = NT2AM(ISYMTR)
      NTAMP   = NTAMP1 + NTAMP2
      NDIM = 3*NNZAL
C
      IF (CCS.OR.CC2)  CALL QUIT('CCMM_TGTBAR: Not implemented for CCS
     &                       or CC2')
C
      IF (DISCEX) CALL QUIT('CCMM_TGTBAR: Not implemented for DISCEX')
C
C------------------------------------
C     Dynamical allocation for CCMM :
C------------------------------------
C
      KRAAO  = 1
      KRAAOx = KRAAO + 3*NNBASX
      KRAAOy = KRAAOx + N2BST(ISYMTR)
      KRAAOz = KRAAOy + N2BST(ISYMTR)
      KWRK1  = KRAAOz + N2BST(ISYMTR)
      LWRK1  = LWORK   - KWRK1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insuff. work in CCMM_TGTBAR, 1')
      END IF
CC
      KDENS = KWRK1
      KWRK2 = KDENS + N2BST(ISYMTR)
      LWRK2 = LWORK - KWRK2
      IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CCMM_TGTBAR, 2')
C
      KETA = KWRK2
      KWRK3 = KETA + NTAMP
      LWRK3 = LWORK - KWRK3
      IF (LWRK3 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK3
         CALL QUIT( 'Too little work in CCMM_TGTBAR, 2')
      END IF

      KTG = KWRK3 
      KWRK4 = KTG + N2BST(ISYMTR)
      LWRK4 = LWORK - KWRK4
      IF (LWRK4.LT.0) CALL QUIT( 'Too little work in CCMM_TGTBAR, 4')

      KRELAY = KWRK4
      KEFIEL = KRELAY + NDIM*NDIM
      KINDIP = KEFIEL + NDIM
      KNSAO = KINDIP + NDIM
      KWRK5 = KNSAO + NNBASX
      LWRK5 = LWORK - KWRK5
      IF (LWRK5.LT.0) CALL QUIT('Too little work in CCMM_TGTBAR, 5')

C Initialize the working matrices
      CALL DZERO(WORK(KRAAO),3*NNBASX)
      CALL DZERO(WORK(KRAAOx),3*N2BST(ISYMTR))
      CALL DZERO(WORK(KDENS),N2BST(ISYMTR))
      CALL DZERO(WORK(KETA),NTAMP) 
      CALL DZERO(WORK(KTG),N2BST(ISYMTR))
      CALL DZERO(WORK(KRELAY),NDIM*NDIM)
      CALL DZERO(WORK(KEFIEL),NDIM)
C
C
      IF ( IPQMMM .GT. 10 .OR. DEBUG) THEN
        RHO1N = DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
        RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
        WRITE(LUPRI,*) ' Norm of RHO1 in CCMM_LTGB on input:', RHO1N
        WRITE(LUPRI,*) ' Norm of RHO2 in CCMM_LTGB on input:', RHO2N
        RHO1N = DDOT(NT1AM(ISYMTR),CTR1,1,CTR1,1)
        RHO2N = DDOT(NT2AM(ISYMTR),CTR2,1,CTR2,1)
        WRITE(LUPRI,*) ' Norm af C1AM in CCMM_LTGB on input:', RHO1N
        WRITE(LUPRI,*) ' Norm af C2AM in CCMM_LTGB on input:', RHO2N
      ENDIF
C Backup dipole origin
      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      IF (LR.NE.'0' .AND. IPOLTP .GT. 0) THEN
C----------------------------
C     Read the relay matrix from file
C----------------------------

         LUQMMM = -1
         CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &              'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND(LUQMMM)

         IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
            CALL READT(LUQMMM,NDIM*NDIM,WORK(KRELAY))
         ELSE
            CALL QUIT('Problem reading the relay matrix from the
     &                 QMMMIM file')
         END IF
   
         CALL GPCLOSE(LUQMMM,'KEEP')

         IF (IPQMMM .GE. 5) THEN
            WRITE(LUPRI,*) 'The relay matrix is read from file'
            CALL OUTPUT(WORK(KRELAY),1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
         END IF

C----------------------------
C     Construct auxiliary densities
C----------------------------
         CALL CCMM_D1AO(WORK(KDENS),CTR1,CTR2,MODEL,LR,
     *                 WORK(KWRK5),LWRK5)
C Calculate corresponding induced dipoles 
C        IF LR.EQ.0 THEN LOAD DIPOLES INSTEAD
C
C
         LRI=0
C
         DO 200 I = 1,MMCENT

            IF(ZEROAL(I) .EQ. -1) GOTO 200

            DIPORG(1) = MMCORD(1,I)
            DIPORG(2) = MMCORD(2,I)    
            DIPORG(3) = MMCORD(3,I)

            DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
            DIST = SQRT(DIST2)

            IF (DIST .GT. RCUTMM) THEN
                LRI = LRI +1
               IF (LOCDEB) THEN 
                  WRITE(LUPRI,*) 'Skipping site in mu^ind vector
     &                              in CCMM_TGTBAR', I
               ENDIF
               GOTO 200
            ENDIF
            
            KPATOM = 0
            NOCOMP = 3
            TOFILE = .FALSE.
            TRIMAT = .TRUE.
            EXP1VL = .FALSE.
C----------------------------
C     Calculate Rra directly in AO basis:
C----------------------------
            RUNQM3 = .TRUE.
            CALL GET1IN(WORK(KRAAO),'NEFIELD',NOCOMP,WORK(KWRK5),
     *              LWRK5,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,
     *              TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)
            RUNQM3 = .FALSE.

            LABEL = 'GIVE INT'

            CALL CCMMINT(LABEL,LUQMMM,WORK(KRAAOx),WORK(KRAAO),
     &             IRREP, ISYM,IERR,WORK(KWRK5),LWRK5)
            CALL CCMMINT(LABEL,LUQMMM,WORK(KRAAOy),WORK(KRAAO+NNBASX),
     &             IRREP, ISYM,IERR,WORK(KWRK5),LWRK5)
            CALL CCMMINT(LABEL,LUQMMM,WORK(KRAAOz),WORK(KRAAO+2*NNBASX),
     &             IRREP, ISYM,IERR,WORK(KWRK5),LWRK5)
C
C 
C Do the trace operation of the product matrices
C
            PRODx = DDOT(N2BST(ISYMTR),WORK(KRAAOx),1,WORK(KDENS),1)
            PRODy = DDOT(N2BST(ISYMTR),WORK(KRAAOy),1,WORK(KDENS),1)
            PRODz = DDOT(N2BST(ISYMTR),WORK(KRAAOz),1,WORK(KDENS),1)

            WORK(KEFIEL + LRI + 0) = PRODx
            WORK(KEFIEL + LRI + 1) = PRODy
            WORK(KEFIEL + LRI + 2) = PRODz 

 200     CONTINUE


         IF (IPQMMM .GE. 5) THEN
            WRITE(LUPRI,*) 'Transformed electric fields'
            CALL OUTPUT(WORK(KEFIEL),1,NDIM,1,1,NDIM,1,1,LUPRI)
            WRITE(LUPRI,*)
         END IF

C Note that we don't add the fields from the QM nuclei and the multipoles
C since we are only interested in the dipoles from the auxiliary densities.

C Now calculate the corresponding induced dipoles
         NTOTI = MAX(NDIM,1)

         CALL DGEMV('N',NDIM,NDIM,ONE,WORK(KRELAY),NTOTI,
     &            WORK(KEFIEL),1,ZERO,WORK(KINDIP),1)
        
         IF (IPQMMM .GE. 5) THEN
            WRITE(LUPRI,*) 'Transformed induced dipole moments'
            CALL OUTPUT(WORK(KINDIP),1,NDIM,1,1,NDIM,1,1,LUPRI)
            WRITE(LUPRI,*)
         END IF

      ELSE IF ( LR .EQ. '0' .AND. IPOLTP .GT. 0) THEN
         CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDIP))
C-------------------
C     Print section:
C-------------------
C
         IF (IPQMMM .GT. 10) THEN
            WRITE(LUPRI,'(//,A)')
     *        ' +============================================'
     *      //'==============+'
            WRITE(LUPRI,'(A)')
     *        ' | MMcent | <mu^ind>_x         | <mu^ind>_y         |'
     *      //' <mu^ind>_z         |'
            WRITE(LUPRI,'(1X,A)')
     *        '+============================================'
     *      //'==============+'
C
            LRI = 0
C
            DO 250 I = 1, MMCENT

               IF (ZEROAL(I) .EQ. -1) GOTO 250
C
C
               WRITE(LUPRI,'(1X,A,I3,A,F16.10,A,
     &                         F16.10,A,F16.10,A)')
     *           '| ', I,'|', WORK(KINDIP + LRI ),' |',
     *           WORK(KINDIP + LRI +1 ),' |', 
     *           WORK(KINDIP + LRI +2),' |'
               WRITE(LUPRI,'(1X,A)')
     *           '+---------------------------------------------'
     *           //'-------------+'
               LRI = LRI + 3
 250        CONTINUE 
            WRITE(LUPRI,'(//,A)')
         END IF
      END IF
     

C To complete the contruction of the Tg operator we need to
C construct the E integrals again - note that  if
C LR==0 so it is actually the first time here. 
C Anyway, we deliberately do this integral-direct
C since we want to avoid a possible storage bottleneck for large systems.


      IF ((LR .EQ. '0' .AND. IPOLTP .GT. 0) .OR. LR.NE.'0') THEN
         LRI = 0

         DO 300 I=1,MMCENT

            IF (ZEROAL(I) .EQ. -1) GOTO 300
         
            DIPORG(1) = MMCORD(1,I)
            DIPORG(2) = MMCORD(2,I)
            DIPORG(3) = MMCORD(3,I)

            DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *              (MMCORD(2,I)-QMCOM(2))**2 +
     *              (MMCORD(3,I)-QMCOM(3))**2
            DIST  = SQRT(DIST2)

            IF (DIST .GT. RCUTMM) THEN
               LRI = LRI+3
               IF (LOCDEB) THEN
                  WRITE(LUPRI,*) 'Skipping site in mu^ind vector
     &                           in CCMM_TGTBAR 2', I
               END IF
               GOTO 300
            END IF

            CALL DZERO(WORK(KRAAO),3*NNBASX)
            CALL DZERO(WORK(KNSAO),NNBASX)

            KPATOM = 0
            NOCOMP = 3
            TOFILE = .FALSE.
            TRIMAT = .TRUE.
            EXP1VL = .FALSE.

            RUNQM3 = .TRUE.
            CALL GET1IN(WORK(KRAAO),'NEFIELD',NOCOMP,WORK(KWRK5),
     &               LWRK5,LABINT,INTREP,INTADR,I,TOFILE,
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)
            RUNQM3 = .FALSE.
       
            FACx = -WORK(KINDIP + LRI + 0)
 
            FACy = -WORK(KINDIP + LRI + 1)

            FACz = -WORK(KINDIP + LRI + 2)
        
             
            CALL DAXPY(NNBASX,FACx,WORK(KRAAO),
     *               1,WORK(KNSAO),1)

            CALL DAXPY(NNBASX,FACy,WORK(KRAAO+NNBASX),
     *               1,WORK(KNSAO),1)

            CALL DAXPY(NNBASX,FACz,WORK(KRAAO+2*NNBASX),
     *               1,WORK(KNSAO),1)

            LRI = LRI + 3
 300     CONTINUE

C   sne: This workspace allocation should be temporary. There
C   must be tons of room earlier in WORK.
         KINT = KWRK5
         KWRK6 = KINT + N2BST(ISYMTR)
         LWRK6 = LWORK - KWRK6
         IF (LWRK6.LT.0) CALL QUIT( 'Too little work in CCMM_TGTBAR' )
         LABEL = 'GIVE INT'
         CALL DZERO(WORK(KINT),N2BST(ISYMTR))

         CALL CCMMINT(LABEL,LUQMMM,WORK(KINT),WORK(KNSAO),
     &           IRREP, ISYM,IERR,WORK(KWRK6),LWRK6)

         CALL DAXPY(N2BST(ISYMTR),ONE,WORK(KINT),1,WORK(KTG),1)

      END IF
C Now the electrostatic contribution. This was once called OLDTG
C--------------------------------------
C     Now the electrostatic contribution to Tg. 
C     We get the multipole integrals from the HF/MM
C     calculation available from file.
C     We do all the corrections simultaneosly. We will analyze
C     the HF/MM results if we want to consider each multipole contribution
C     separately.

C 1) The charge contribution
C
C--------------------------
C     Dynamical allocation:
C--------------------------
C
      IF ( LR .EQ. '0') THEN
         IF (NMULT .GE. 0) THEN

C
            LABEL = 'GIVE INT'
            CALL DZERO(WORK(KINT),N2BST(ISYMTR))
            CALL DZERO(WORK(KNSAO),NNBASX)

            LUQMMM = -1
            CALL GPOPEN(LUQMMM,'MU0INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND (LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KNSAO))
     
            CALL CCMMINT(LABEL,LUQMMM,WORK(KINT),WORK(KNSAO),
     &             IRREP, ISYM,IERR,WORK(KWRK6),LWRK6)

C sne: Minus should already be carried in the integrals cf sirqm3
            CALL DAXPY(N2BST(ISYMTR),ONE,WORK(KINT),
     *                   1,WORK(KTG),1)

            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF
C 2) The dipole contribution
         IF (NMULT .GE. 1) THEN

            CALL DZERO(WORK(KNSAO),NNBASX)
            CALL DZERO(WORK(KINT),N2BST(ISYMTR))

            CALL GPOPEN(LUQMMM,'MU1INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND (LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KNSAO))
            CALL CCMMINT(LABEL,LUQMMM,WORK(KINT),WORK(KNSAO),
     &             IRREP, ISYM,IERR,WORK(KWRK6),LWRK6)

C sne: Minus should be correct in dipole case cf sirqm3
            CALL DAXPY(N2BST(ISYMTR),-ONE,WORK(KINT),
     *                   1,WORK(KTG),1)
C
            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF

C 3) The quadrupole contribution
         IF (NMULT .GE. 2) THEN 
            CALL DZERO(WORK(KNSAO),NNBASX)
            CALL DZERO(WORK(KINT),N2BST(ISYMTR))

            CALL GPOPEN(LUQMMM,'MU2INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND (LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KNSAO))

C sne: Minus should already be carried in the integrals cf sirqm3
            CALL DAXPY(N2BST(ISYMTR),ONE,WORK(KINT),
     *                   1,WORK(KTG),1)

            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF

C 4) The octopole contribution

         IF (NMULT .GE. 3) THEN
            CALL DZERO(WORK(KNSAO),NNBASX)
            CALL DZERO(WORK(KINT),N2BST(ISYMTR))

            CALL GPOPEN(LUQMMM,'MU3INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND (LUQMMM)
            CALL READT(LUQMMM,NNBASX,WORK(KNSAO))

C sne: Minus should already be carried in the integrals cf sirqm3
            CALL DAXPY(N2BST(ISYMTR),ONE,WORK(KINT),
     *                   1,WORK(KTG),1)

            CALL GPCLOSE(LUQMMM,'KEEP')
         END IF
      END IF
   


C---------------------------------------------------------------------
C             (LR.EQ.'L').OR.(LR.EQ.'F')                Calculate eta^Geff 
C             (LR.EQ.'R').OR.(LR.EQ.'P').OR.(LR.EQ.'0') Calculate xi^Geff
C---------------------------------------------------------------------
C
      KETA1 = KETA
      KETA2 = KETA1 + NTAMP1
C      

      IF ( (LR.EQ.'L') .OR. (LR.EQ.'F') ) THEN
         LABEL = 'GIVE INT'
         LIST  = 'L0'
         IDLINO = 1
C
         CALL CC_ETAC(ISYMTR,LABEL,WORK(KETA),
     *             LIST,IDLINO,0,WORK(KTG),WORK(KWRK6),LWRK6)
C
C if LR=0 also do the XKSI transformation of KTG
      ELSE IF ( (LR.EQ.'R') .OR. (LR.EQ.'P') .OR. (LR.EQ.'0')) THEN 
         LABEL = 'GIVE INT'
C
         CALL CC_XKSI(WORK(KETA),LABEL,ISYMTR,0,WORK(KTG),
     *                 WORK(KWRK6),LWRK6) 
C
         IF (LR .EQ. 'R') THEN
            CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYMTR)
         END IF
      END IF
C

      IF(DEBUG .OR. IPQMMM .GT. 14) THEN
         TAL1= DDOT(NT1AM(ISYMTR),WORK(KETA1),1,WORK(KETA1),1)
         TAL2= DDOT(NT2AM(ISYMTR),WORK(KETA2),1,WORK(KETA2),1)
         WRITE(LUPRI,*) 'Printing TGTBAR contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF

         CALL DAXPY(NT1AM(ISYMTR),ONE,WORK(KETA1),1,RHO1,1)
         CALL DAXPY(NT2AM(ISYMTR),ONE,WORK(KETA2),1,RHO2,1)
C
      IF (DEBUG .OR. IPQMMM .GT. 14) THEN
         TAL1= DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
         TAL2= DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
         WRITE(LUPRI,*) 'Printing RHO:
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF
C Restore dipole origin
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
         
      END
C
C***********************************************************************
C  /* Deck ccmm_d1ao */
      SUBROUTINE CCMM_D1AO(AODEN,CTR1,CTR2,MODEL,LR,
     *                    WORK,LWORK)
C
C     Kristian Sneskov July 2009 inspired by CCX_D1AO
C
C     Purpose: To calculate contributions to the one electron 
C              auxiliary density matrix appearing in TINE CCMM_LTRB 
C              and return it backtransformed to AO-basis in AODEN.
C
C     IF(LR=R.OR.LR=F) -> Calc D^C
C     IF(LR=L.OR.LR=P) -> Calc D^B
C     Current models: CCS, CC2 and CCSD    
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccexci.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <cclr.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK), CTR1(*), CTR2(*)

      CHARACTER MODEL*10,MODDUM*10
      CHARACTER LR*1
C
      CALL QENTER('CCMM_D1AO')
C
C---------------------------------------------
C     Find symmetry of left and right vectors. 
C---------------------------------------------
C
      ISYML = ILSTSYM('L0',1) ! Symmetry of multipliers and b-trial
      ISYMR = ILSTSYM('R0',1) ! symmetry of (t2)amplitudes
      IF (ISYMR .NE. ISYML) 
     &     CALL QUIT('CCMM_D1AO: Density not total sym.')
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
      KONEAI = 1
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KL1AM  = KONEIA + NT1AMX
      KL2AM  = KL1AM  + NT1AM(ISYML) ! stores t2 when D^B and tbar2 when D^C
      KT1AM  = KL2AM  + NT2SQ(ISYML)
      KR2AM  = KT1AM  + NT1AM(ISYMR)
      KEND0  = KR2AM  + NT2AM(ISYMTR)
C
      KR2AMT = KEND0
      KEND1  = KR2AMT  + NT2AM(ISYMTR)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('CCMM_D1AO: Not enough memory')
      ENDIF
C
C--------------------------------------------
C     Initialize one electron density arrays.
C--------------------------------------------
C
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      CALL DZERO(WORK(KONEAI),NT1AMX)
      CALL DZERO(WORK(KONEIA),NT1AMX)
C--------------------------------------------
C     Initialize working arrays.
C     NOTE if D^B -> KL1AM is just a dummy array
C--------------------------------------------
      CALL DZERO(WORK(KL1AM),NT1AM(ISYML))
      CALL DZERO(WORK(KL2AM),NT2SQ(ISYML))
      CALL DZERO(WORK(KT1AM),NT1AM(ISYMR))
      CALL DZERO(WORK(KR2AM),NT2AM(ISYMTR))
      CALL DZERO(WORK(KR2AMT),NT2AM(ISYMTR))
C
C----------------------
C     Read left vector.
C----------------------
C
C
      IF((LR.EQ.'R').OR.(LR.EQ.'F')) THEN
         IOPT = 3
         CALL CC_RDRSP('L0',1,ISYML,IOPT,MODDUM,
     *             WORK(KL1AM),WORK(KR2AMT))
C Read in t1 amplitudes
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
C     Copy CTR2 to R2AM since we will scale diagonal
         CALL DCOPY(NT2AM(ISYMR),CTR2,1,WORK(KR2AM),1) 

      IF (ISYMR.EQ.1) CALL CCLR_DIASCL(WORK(KR2AM),TWO,ISYMR)
         IF (DEBUG) THEN
            TAL1 = DDOT(NT1AMX,WORK(KL1AM),1,WORK(KL1AM),1)
            TAL2 = DDOT(NT2AMX,WORK(KR2AMT),1,WORK(KR2AMT),1)
            WRITE(LUPRI,*) 'CCMM_TGDEN: Norm of L1: ', TAL1, 
     &                     'Norm of L2: ', TAL2
         END IF
      ELSE IF ( (LR .EQ. 'L') .OR. (LR .EQ. 'P') ) THEN
         CALL DCOPY(NT2AM(ISYMTR),CTR2,1,WORK(KR2AMT),1)
         IOPT = 3
         CALL CC_RDRSP('R0',1,ISYMR,IOPT,MODDUM,
     *              WORK(KT1AM),WORK(KR2AM)) 
         IF (DEBUG) THEN
            TAL1 = DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
            TAL2 = DDOT(NT2AMX,WORK(KR2AM),1,WORK(KR2AM),1)
            WRITE(LUPRI,*) 'CCMM_TGDEN: Norm of T1: ', TAL1, 
     &                     'Norm of T2: ', TAL2
         END IF
      END IF
C
C----------------------
C     Square up tbar_2 multipliers or b_2 elements
C----------------------
C
      CALL CC_T2SQ(WORK(KR2AMT),WORK(KL2AM),ISYML)
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
C Note that diascl amplitudes are expected 
      IF (.NOT. MP2) THEN
         CALL DCOPY(NT2AM(ISYMR),WORK(KR2AM),1,WORK(KR2AMT),1) 
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KR2AMT),1.0D0,ISYMR,IOPTTCME)
      ENDIF
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ia) is stored transposed
C-------------------------------
C
      KXMAT  = KEND1
      KYMAT  = KXMAT  + NMATIJ(1)
      KEND2  = KYMAT  + NMATAB(1)
      LWRK2  = LWORK  - KEND1
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 2. alloc. in CCMM_D1AO')
      ENDIF
      CALL DZERO(WORK(KXMAT),NMATIJ(1))
      CALL DZERO(WORK(KYMAT),NMATAB(1))
C
C-----------------------------------------------------------
C     Calculate X-intermediate of L2AM- and R2AM-amplitudes.
C-----------------------------------------------------------
C Note that KL2AM must be squared when input in these routines
C
      CALL CC_XI(WORK(KXMAT),WORK(KL2AM),ISYML,WORK(KR2AM),ISYMR,
     *           WORK(KEND2),LWRK2) 
C
C-----------------------------------------------------------
C     Calculate Y-intermediate of L2AM- and R2AM-amplitudes.
C-----------------------------------------------------------
C
      CALL CC_YI(WORK(KYMAT),WORK(KL2AM),ISYML,WORK(KR2AM),ISYMR,
     *           WORK(KEND2),LWRK2) 
C
C--------------------------------------------------------------
C     Construct <L2|[Emn,R2]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND2),LWRK2,1)
      CALL DAXPY(NMATIJ(1),-ONE,WORK(KXMAT),1,WORK(KONEIJ),1)
C
      IF ( (LR .EQ. 'R') .OR. (LR .EQ. 'F') ) THEN
C
C--------------------------------------------------------------
C     Construct <HF|[Emn,R1]|HF> contribution.
C--------------------------------------------------------------
C
         CALL DAXPY(NT1AMX,TWO,CTR1,1,WORK(KONEIA),1)
C
C
C--------------------------------------------------------------
C     Construct <L1|[Emn,R1]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
         CALL CC_DXIJ(WORK(KONEIJ),WORK(KL1AM),ISYML,CTR1,ISYMR)
         CALL CC_DXAB(WORK(KONEAB),WORK(KL1AM),ISYML,CTR1,ISYMR)
C
C--------------------------------------------------------------
C     Construct <L1|[Eia,R2]|HF> contribution to DXia
C--------------------------------------------------------------
C
         CALL CC_DXIA12(WORK(KONEIA),WORK(KL1AM),ISYML,
     *         WORK(KR2AMT),ISYMR)
C
C----------------------------------------------------------
C     Construct <L2|[[Eia,R1],T2]|HF> contribution to DXia.
C----------------------------------------------------------
         KT2AM  = KEND0
         KEND3  = KT2AM  + NT2AM(ISYMOP) 
         LWRK3  = LWORK  - KEND3 
C
         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for 3. alloc. in CCMM_D1AO')
         ENDIF
         CALL CC_DXIA21(WORK(KONEIA),WORK(KL2AM),ISYML,
     *               CTR1,ISYMR,WORK(KT2AM),ISYMR,
     *               WORK(KEND3),LWRK3)

      ELSE IF ((LR .EQ. 'L') .OR. (LR .EQ. 'P')) THEN
C
C----------------------------------------------------------
C     Construct <B1|[Eai|HF> contribution to DXai.
C----------------------------------------------------------
         CALL DAXPY(NT1AMX,ONE,CTR1,1,WORK(KONEAI),1)

C--------------------------------------------------------------
C     Construct <B1|[Eia,T2]|HF> contribution to DXia
C--------------------------------------------------------------
         CALL CC_DXIA12(WORK(KONEIA),CTR1,ISYMTR,WORK(KR2AMT),ISYMR)
      END IF 
C
C--------------------------
C     Write out test norms.
C--------------------------
C
      IF (DEBUG) THEN
         XIJ = DDOT(NMATIJ(1),WORK(KONEIJ),1,WORK(KONEIJ),1)
         XAB = DDOT(NMATAB(1),WORK(KONEAB),1,WORK(KONEAB),1)
         XAI = DDOT(NT1AMX,WORK(KONEAI),1,WORK(KONEAI),1)
         XIA = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
         WRITE(LUPRI,*) 'Norms: DXIJ',XIJ
         WRITE(LUPRI,*) 'Norms: DXAB',XAB
         WRITE(LUPRI,*) 'Norms: DXAI',XAI
         WRITE(LUPRI,*) 'Norms: DXIA',XIA
      ENDIF
C
C----------------------------------
C     Calculate the lambda matrices.
C----------------------------------
C
      KLAMDP = KEND0
      KLAMDH = KLAMDP + NLAMDT
      KEND4  = KLAMDH + NLAMDT
      LWRK4  = LWORK  - KEND4 
      IF (LWRK4 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 4. alloc. in CCMM_D1AO')
      ENDIF
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND4),
     *            LWRK4)
C
C--------------------------------------------------------
C     Add the one electron density in the AO-basis.
C--------------------------------------------------------
CC
      ISDEN = 1
      CALL CC_DENAO(AODEN,ISDEN,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KEND4),LWRK4)
C
      CALL QEXIT('CCMM_D1AO')
C
   1  FORMAT(1x,A35,1X,E20.10)
      RETURN
      END
C
