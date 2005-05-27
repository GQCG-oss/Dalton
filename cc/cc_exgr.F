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
c*DECK CC_EXGR
       SUBROUTINE CC_EXGR(WORK,LWORK)
C
C----------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster analytical
C              excited state first order properties and gradients.
C
C              CIS, CCS, CC2, CCSD
C
C     Solves for excited state t-bar amplitudes 
C                           = Lagrangian multipliers.
C     Calculates first order properties: dipole moment, 
C     quadrupole moment, electric field gradients, 
C     relativistic corrections, electronic moments.
C
C     Written by Ove Christiansen April 1997.
C     Restructured for incorporating CC3 by Christof Haettig, 2003.
C---------------------------------------------------------------------
      IMPLICIT NONE
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <codata.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccfop.h>
#include <ccexgr.h>
#include <ccexci.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccrspprp.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccfield.h>
#include <ccroper.h>
#include <ccfro.h>
#include <exeinf.h>
#include <infvar.h>
#include <dipole.h>
#include <quadru.h>
#include <nqcc.h>
#include <ccnoddy.h>

      INTEGER IZERO, LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK), ELSEMO(3,3), SKODE(3,3), SKODN(3,3)
      REAL ZERO, ONE, TWO, XD, XN, EXC, FF, DUM, DDOT
      REAL DARW, VELO, PROP
#else
      DOUBLE PRECISION WORK(LWORK), ELSEMO(3,3), SKODE(3,3), SKODN(3,3)
      DOUBLE PRECISION ZERO, ONE, TWO, XD, XN, EXC, FF, DUM, DDOT
      DOUBLE PRECISION DARW, VELO, PROP
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)

      CHARACTER MODEL*10, MODDUM*10
      CHARACTER MODELPRI*4, MODELPRI2*30
      CHARACTER LABEL*8
      CHARACTER LABR12*3

      CHARACTER*(8) FNE0DIA, FNE0DIA2, FNE0DAB, FNE0DIJ
      PARAMETER (FNE0DIA='CC_E0DIA',FNE0DIA2='CCE0DIA2',
     &           FNE0DAB='CC_E0DAB',FNE0DIJ ='CC_E0DIJ')

      INTEGER IX, IEX, ISYMX, KDENS, KWRK2, LWRK2, IXNR, ILLNR, KDENS2,
     &        KWRK3, LWRK3, ILRNR, IASGER, IDIP, KONEP, KWRK4, LWRK4,
     &        ISY, IOPT, IJ, LENGTH, KDIAG, KAXIS, ICCPRI, IRC, IOP
C
C#include <leinf.h>
C
      CALL QENTER('CC_EXGR')
C
C------------------------------------
C     Header of Property calculation.
C------------------------------------
C
C     CALL TIMER('START ',TIMEIN,TIMOUT)
C
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A)')
     *'*<<<<<<<<<<<<< OUTPUT FROM COUPLED CLUSTER RESPONSE  >>>>'//
     *'>>>>>>>>>*'
      IF ( CCFOP  ) THEN
         WRITE (LUPRI,'(1X,A)')
     *   '*                                                        '//
     *   '         *'
         WRITE (LUPRI,'(1X,A)')
     *   '*<<<<<  CALCULATION OF EXCITED STATE FIRST ORDER PROPERTI'//
     *   'ES >>>>>>*'
      ENDIF
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A,/)')
     *'*********************************************************'//
     *'**********'
C
      LABR12 = ' '
      IF (CCR12) LABR12 = '-R12'
C
      MODEL = 'CCSD'
      IF (CC2) THEN
         CALL AROUND('Coupled Cluster model is: CC2'//LABR12)
         MODEL = 'CC2'
         MODELPRI = ' CC2'
      ENDIF
      IF (CCS.AND.(.NOT.CIS)) THEN
         CALL AROUND('Coupled Cluster model is: CCS')
         MODEL = 'CCS'
         MODELPRI = ' CCS'
      ENDIF
      IF (CCS.AND.CIS) THEN
         CALL AROUND('Model is not CC but CIS crap.')
         MODEL = 'CCS'
         MODELPRI = ' CIS'
      ENDIF
      IF (CCD) THEN
         CALL AROUND('Coupled Cluster model is: CCD'//LABR12)
         MODEL = 'CCD'
         MODELPRI = ' CCD'
      ENDIF
      IF (CC3  ) THEN
         CALL AROUND('Coupled Cluster model is: CC3'//LABR12)
         MODEL = 'CC3'
         MODELPRI = ' CC3'
      ENDIF
      IF (CC1A) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1a'//LABR12)
         MODEL = 'CCSDT-1a'
         WRITE(LUPRI,*) 
     *    'CCSDT-1a X-state first order properties not implemented'
         CALL QEXIT('CC_EXGR')
         RETURN
      ENDIF
      IF (CC1B) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1b'//LABR12)
         MODEL = 'CCSDT-1b'
         WRITE(LUPRI,*) 
     *    'CCSDT-1b X-state first order properties not implemented'
         CALL QEXIT('CC_EXGR')
         RETURN
      ENDIF
      IF (CCSD) THEN
         CALL AROUND('Coupled Cluster model is: CCSD'//LABR12)
         MODEL = 'CCSD'
         MODELPRI = 'CCSD'
      ENDIF
      MODDUM=MODEL
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_EXGR-1: Workspace:',LWORK

*----------------------------------------------------------------------*
* Loop over (selected) excited states:
*----------------------------------------------------------------------*
      DO 3000 IX = 1, NXGRST
         IEX   = IXGRST(IX)
         ISYMX = ISYEXC(IEX)
         IXNR  = IEX - ISYOFE(ISYMX)
         WRITE(LUPRI,*) 'Total, actual nr.,loop index',NXGRST,IEX,IX
         WRITE(LUPRI,*) 'Sym, reduced nr.',ISYMX,IXNR

*======================================================================*
* Section 1: compute the excited state 1e^- AO density matrix
*======================================================================*
         KDENS = 1    
         KWRK2 = KDENS  + N2BST(ISYMOP)
         LWRK2 = LWORK  - KWRK2

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
            CALL QUIT('Insufficient memory for one '//
     &           'e-density in CC_EXGR')
         ENDIF

*----------------------------------------------------------------------*
* Calculate contribution from <HF + Tbar^f|E_pq|CC>:
* (this is done using the ground state routines, but interchanging
*  Tbar^0 ("L0 ") by Tbar^f ("E0 "). Note that for CCS here only the
*  ground state SCF density is calculated, for the contribution of
*  Tbar^f see below)  
*----------------------------------------------------------------------*
         IF (CCS) THEN
           CALL CCS_D1AO(WORK(KDENS),WORK(KWRK2),LWRK2)
         ELSE
           ILLNR = IX

           IF (CC3) THEN
             ! calculate CC3 contributions in MO basis and pass them 
             ! on the direct access files FNE0D* to CC_D1AO
             IF (NODDY_DEN) THEN
               CALL CCSDT_XI_CONT_NODDY('E0 ',DUMMY,ILLNR,1,
     &                                  IDUMMY,IDUMMY,0,0,
     &                                  .TRUE.,FNE0DIA,FNE0DIA2,
     &                                   FNE0DAB,FNE0DIJ,
     &                                  WORK(KWRK2),LWRK2) 
             ELSE

               CALL QUIT('CC_EXGR> Real code for CC3 Xi den. missing!')
     
             END IF
           END IF

           CALL CC_D1AO(WORK(KDENS),DUM,WORK(KWRK2),LWRK2,MODDUM,
     *                  'E0 ',ILLNR,.FALSE.,
     &                  FNE0DIA,FNE0DIA2,FNE0DAB,FNE0DIJ)
         END IF

         IF (FROIMP .OR. FROEXP) THEN
           CALL CC_FCD1AO(WORK(KDENS),WORK(KWRK2),LWRK2,MODDUM)
         ENDIF

         IF (DEBUG) THEN
            XD = DDOT(N2BST(ISYMOP),WORK(KDENS),1,WORK(KDENS),1)
            WRITE(LUPRI,*) 'Norm of dens-1: ',XD
         ENDIF

*----------------------------------------------------------------------*
* Calculate contribution from <Ebar^f|[E_pq,E^f]|CC>, where Ebar^f
* and E^f are the left and right eigenvectors (on lists "LE " and "RE ")
* Note that for CCS the contribution from <Tbar^f|E_pq|CC> is included
* here
*----------------------------------------------------------------------*
         KDENS2 = KWRK2
         KWRK3  = KDENS2 + N2BST(ISYMOP)
         LWRK3  = LWORK - KWRK3
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
            CALL QUIT('Insufficient memory for one '//
     &           'e-density in CC_EXGR')
         ENDIF

         CALL DZERO(WORK(KDENS2),N2BST(ISYMOP))
         ILLNR = IEX
         ILRNR = ILLNR

         IF (CCS) THEN
           CALL CCSX_D1AO(WORK(KDENS2),WORK(KWRK3),LWRK3,
     *                   'LE ',ILLNR,'RE ',ILRNR,'E0 ',IX)
         ELSE
           CALL CCX_D1AO(WORK(KDENS2),WORK(KWRK3),LWRK3,MODDUM,
     *                  'LE ',ILLNR,'RE ',ILRNR)
         END IF
         CALL DAXPY(N2BST(ISYMOP),ONE,WORK(KDENS2),1,WORK(KDENS),1)

         IF (DEBUG) THEN
            XD = DDOT(N2BST(ISYMOP),WORK(KDENS),1,WORK(KDENS),1)
            WRITE(LUPRI,*) 'Norm of dens-2: ',XD
         ENDIF

         IF (IPRINT .GT. 50) THEN
            CALL AROUND('One electron unrelaxed density in cc_exgr')
            CALL CC_PRFCKAO(WORK(KDENS),1)
         ENDIF

         KWRK3 = KWRK2
         LWRK3 = LWRK2

*======================================================================*
* Section 2: evaluate different properties and print:
*======================================================================*
         MODELPRI2 = 'Unrelaxed excited state '//MODELPRI
         CALL AROUND(MODELPRI2//' First-order properties: ')

         ISYMX = ISYEXC(IEX)
         IXNR  = IEX - ISYOFE(ISYMX)
         EXC = EIGVAL(IEX) 
         WRITE(LUPRI,'(/,1X,A,/,1X,A,I5,/,1X,A,I5,/,1X,A,F10.6,F12.6)') 
     *    'Excited state properties for ',
     *    'Excited state nr.:           ',IXNR,
     *    'Excited state sym:           ',ISYMX,
     *    'Excitation energy (au.,eV):  ',EXC,EXC*XTEV
         CALL FLSHFO(LUPRI)
C
C=======================================
C     Calculate molecular dipole moment.
C=======================================
C
      IF (XDIPMO) THEN
C
         CALL AROUND(' Electric Dipole Moment ')
C
C-------------------------------------------
C        Calculate the nuclear contribution.
C-------------------------------------------
C
         IASGER = IPRINT - 5
         CALL DIPNUC(WORK(KWRK3),WORK(KWRK3),IASGER,.FALSE.)
C
         DO 100 IDIP = 1,3
C
            IF (IDIP .EQ. 1) LABEL = 'XDIPLEN '
            IF (IDIP .EQ. 2) LABEL = 'YDIPLEN '
            IF (IDIP .EQ. 3) LABEL = 'ZDIPLEN '
C
C----------------------------------
C           get property integrals.
C----------------------------------
C
            KONEP  = KWRK3
            KWRK4  = KONEP  + N2BST(ISYMOP)
            LWRK4  = LWORK  - KWRK4
C
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
               CALL QUIT('Insufficient memory for '//
     &              'DIPLEN-int. in CC_EXGR')
            ENDIF
C
            CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
            FF = 1.0D0
            ISY = -1
            CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
            IF (IPRINT .GT. 50) THEN
               CALL AROUND('One electron property integrals in cc_fop')
               CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
            ENDIF
C
C----------------------------------------------
C        Calculate the electronic contribution.
C----------------------------------------------
C
            IF (ISY .EQ. 1 ) THEN
               DIPME(IDIP) = -DDOT(N2BST(ISYMOP),WORK(KONEP),1,
     *                             WORK(KDENS),1)
            ELSE
               DIPME(IDIP) = 0
            ENDIF
            DIPMN(IDIP) = DIPMN(IDIP) + DIPME(IDIP)
C
C--------------------------------
C           Store on prpc list.
C--------------------------------
C
            CALL CC_PRPC(DIPMN(IDIP),MODEL,-11,LABEL,LABEL,LABEL,LABEL,
     *                   EXC,DUMMY,DUMMY,ISY,ISYMX,1,IXNR)
C
  100    CONTINUE
C
C---------------------
C        Print result.
C---------------------
C
         IF (IASGER .GT. 0) THEN
            CALL HEADER('Electronic contribution to dipole moment',-1)
            CALL DP0PRI(DIPME)
         ENDIF
         CALL HEADER('Total Molecular Dipole Moment',-1)
         CALL DP0PRI(DIPMN)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
C===========================================
C     Calculate molecular quadrupole moment.
C===========================================
C
      IF (XQUADR) THEN
C
         CALL AROUND(' Electric Quadrupole Moment ')
C
C-------------------------------------------
C        Calculate the nuclear contribution.
C-------------------------------------------
C
         IOPT   = 1
         IASGER = -1
         CALL CCNUCQUA(WORK(KWRK3),LWRK3,IOPT,IASGER)
         CALL DZERO(QDREL,3*3)
C
         IJ = 0
         DO 110 I = 1,3
            DO 120 J = I,3
               IJ = IJ + 1
C
               IF (IJ .EQ. 1) LABEL = 'XXTHETA '
               IF (IJ .EQ. 2) LABEL = 'XYTHETA '
               IF (IJ .EQ. 3) LABEL = 'XZTHETA '
               IF (IJ .EQ. 4) LABEL = 'YYTHETA '
               IF (IJ .EQ. 5) LABEL = 'YZTHETA '
               IF (IJ .EQ. 6) LABEL = 'ZZTHETA '
C
C-------------------------------------
C              get property integrals.
C-------------------------------------
C
               KONEP  = KWRK3
               KWRK4  = KONEP  + N2BST(ISYMOP)
               LWRK4  = LWORK  - KWRK4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
                  CALL QUIT('Insufficient memory for THETA-int. '//
     &                 'in CC_EXGR')
               ENDIF
C
               CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
               FF = 1.0D0
               ISY = -1
               CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
               IF (IPRINT .GT. 50) THEN
                  CALL AROUND('One electron property int. in cc_fop')
                  CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
               ENDIF
C
C-------------------------------------------------
C           Calculate the electronic contribution.
C-------------------------------------------------
C
               LENGTH = N2BST(ISYMOP)
C
               IF (ISY .EQ. 1) THEN
                 CALL CCELQUA(WORK(KONEP),WORK(KDENS),LENGTH,I,J,QDREL)
               ENDIF
C
  120       CONTINUE
  110    CONTINUE
C
C------------------------
C        Reorder storing.
C------------------------
C
         CALL CC_QUAREO(QDREL,SKODE)
         CALL CC_QUAREO(QDRNUC,SKODN)
C
C---------------------
C        Print result.
C---------------------
C
         IF (IPRINT .GT. 9) THEN
            CALL HEADER('Nuclear contr. to quadrupole moment',-1)
            WRITE(LUPRI,474) 'X','Y','Z'
            CALL OUTPUT(SKODN,1,3,1,3,3,3,1,LUPRI)
            CALL HEADER('Electronic contr. to quadrupole moment',-1)
            WRITE(LUPRI,474) 'X','Y','Z'
            CALL OUTPUT(SKODE,1,3,1,3,3,3,1,LUPRI)
         ENDIF
         CALL DAXPY(9,-ONE,SKODE,1,SKODN,1)
         CALL HEADER('Total Molecular quadrupole moment',-1)
         WRITE(LUPRI,474) 'X','Y','Z'
         CALL OUTPUT(SKODN,1,3,1,3,3,3,1,LUPRI)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
C==================================================
C     Calculate electronic second moment of charge.
C==================================================
C
      IF (XSECMO) THEN
C
         CALL AROUND(' Electronic second moment of charge ')
C
         CALL DZERO(ELSEMO,9)
C
         IJ = 0
         DO 115 I = 1,3
            DO 125 J = I,3
               IJ = IJ + 1
C
               IF (IJ .EQ. 1) LABEL = 'XXSECMOM'
               IF (IJ .EQ. 2) LABEL = 'XYSECMOM'
               IF (IJ .EQ. 3) LABEL = 'XZSECMOM'
               IF (IJ .EQ. 4) LABEL = 'YYSECMOM'
               IF (IJ .EQ. 5) LABEL = 'YZSECMOM'
               IF (IJ .EQ. 6) LABEL = 'ZZSECMOM'
C
C-------------------------------------
C              get property integrals.
C-------------------------------------
C
               KONEP  = KWRK3
               KWRK4  = KONEP  + N2BST(ISYMOP)
               LWRK4  = LWORK  - KWRK4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
                  CALL QUIT('Insufficient memory for '//
     &                 'SECMOM-int. in CC_EXGR')
               ENDIF
C
               CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
               FF = 1.0D0
               ISY = -1
               CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
               IF (IPRINT .GT. 50) THEN
                  CALL AROUND('One electron property int. in cc_fop')
                  CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
               ENDIF
C
C-------------------------------------------------
C           Calculate the electronic contribution.
C-------------------------------------------------
C
               LENGTH = N2BST(ISYMOP)
C
               IF (ISY.EQ.1) THEN
                 CALL CCELQUA(WORK(KONEP),WORK(KDENS),LENGTH,I,J,ELSEMO)
               ENDIF
C
C--------------------------------
C           Store on prpc list.
C--------------------------------
C
            CALL CC_PRPC(ELSEMO(I,J),MODEL,-11,LABEL,LABEL,LABEL,LABEL,
     *                   EXC,DUMMY,DUMMY,ISY,ISYMX,1,IXNR)
C
  125       CONTINUE
  115    CONTINUE
C
C------------------------
C        Reorder storing.
C------------------------
C
         CALL CC_QUAREO(ELSEMO,SKODE)
C
C---------------------
C        Print result.
C---------------------
C
         WRITE(LUPRI,474) 'X','Y','Z'
         CALL OUTPUT(SKODE,1,3,1,3,3,3,1,LUPRI)
         CALL CC_TNSRAN(SKODE,WORK(KWRK3),LWRK3)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
  474 FORMAT(20X,A1,14X,A1,14X,A1)
C
C=======================================
C     Calculate electric field gradient.
C=======================================
C
      IF (XNQCC) THEN
C
         CALL AROUND(' Electric Field Gradients ')
C
C-------------------------------------------
C        Calculate the nuclear contribution.
C-------------------------------------------
C
         IOPT   = 2
         IASGER = IPRINT - 5
         CALL CCNUCQUA(WORK(KWRK3),LWRK3,IOPT,IASGER)
C
C----------------------------------------------
C        Calculate the electronic contribution.
C----------------------------------------------
C
         LENGTH = N2BST(ISYMOP)
         CALL CCELEFG(WORK(KDENS),LENGTH,WORK(KWRK3),LWRK3,IASGER)
C
C---------------------
C        Print result.
C---------------------
C
         KDIAG = KWRK3
         KAXIS = KDIAG + 3*MXCENT
         KWRK4 = KAXIS + 9*MXCENT
         LWRK4 = LWORK - KWRK4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
            CALL QUIT('Insufficient memory for EFG-results in CC_EXGR')
         ENDIF
C
         IASGER = 2
         ICCPRI = 2
         CALL NQCRES(IASGER,WORK(KDIAG),WORK(KAXIS),ICCPRI)
C
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
C=========================================================
C     Relativistic corrections to the ground-state energy.
C=========================================================
C
      IF (XRELCO) THEN
C
         CALL AROUND(' Relativistic corrections to the ground-state'
     *               //' energy ')
C
         DO 130 IRC = 1,2
C
            IF (IRC .EQ. 1) LABEL = 'DARWIN  '
            IF (IRC .EQ. 2) LABEL = 'MASSVELO'
C
C-----------------------------
C           get the integrals.
C-----------------------------
C
            KONEP  = KWRK3
            KWRK4  = KONEP  + N2BST(ISYMOP)
            LWRK4  = LWORK  - KWRK4
C
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
               CALL QUIT('Insufficient memory for '//
     &              'Darwin-int. in CC_EXGR')
            ENDIF
C
            CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
            FF = 1.0D0
            ISY = 1
            CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
            IF (IPRINT .GT. 50) THEN
               CALL AROUND('Relativistic integrals in cc_fop')
               CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
            ENDIF
C
C-------------------------------------
C           Calculate the corrections.
C-------------------------------------
C
            IF (IRC .EQ. 1) THEN
               DARW = DDOT(N2BST(ISYMOP),WORK(KONEP),1,WORK(KDENS),1)
            ELSE IF (IRC .EQ. 2) THEN
               VELO = DDOT(N2BST(ISYMOP),WORK(KONEP),1,WORK(KDENS),1)
            ENDIF
C
  130    CONTINUE
C
C----------------------
C     Write out result.
C----------------------
C
      WRITE(LUPRI,*) ' '
      WRITE(LUPRI,131) 'Darwin term:', DARW, 'Mass-Velocity term:', VELO
      WRITE(LUPRI,132) '----------- ',       '------------------ '
      WRITE(LUPRI,*) ' '
      WRITE(LUPRI,*) ' '
      WRITE(LUPRI,133) 'Total relativistic correction:', DARW + VELO
      WRITE(LUPRI,134) '----------------------------- '
      WRITE(LUPRI,*) ' '
      WRITE(LUPRI,*) ' '
C
  131 FORMAT(9X,A12,1X,F9.6,7X,A19,1X,F9.6)
  132 FORMAT(9X,A12,17X,A19)
  133 FORMAT(19X,A30,1X,F9.6)
  134 FORMAT(19X,A30)
C
      ENDIF
C
C--------------------------------------------------------------
C     Section for general operator APROP represented by LABEL.
C     Note that only the electronic contribution is calculated.
C--------------------------------------------------------------
C
      DO 140 IOP = 1, NAXGRO
C
         LABEL = LBLOPR(IAXGRO(IOP))
C
         IF (IOP .EQ. 1) CALL AROUND( 
     *               ' Electronic contribution to operator ')
C
C--------------------------
C        get the integrals.
C--------------------------
C
         KONEP  = KWRK3
         KWRK4  = KONEP  + N2BST(ISYMOP)
         LWRK4  = LWORK  - KWRK4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
            CALL QUIT('Insufficient memory for property '//
     &           'integrals in CC_EXGR')
         ENDIF
C
         CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
         FF = 1.0D0
         ISY = -1
         CALL CC_ONEP(WORK(KONEP),WORK(KWRK4),LWRK4,FF,ISY,LABEL)
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND('Property integrals in cc_fop')
            CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
         ENDIF
C
C--------------------------------------------------------------------
C        Calculate the electronic contribution to the given property.
C--------------------------------------------------------------------
C
         PROP = DDOT(N2BST(ISYMOP),WORK(KONEP),1,WORK(KDENS),1)
         CALL CC_PRPC(PROP,MODEL,-11,LABEL,LABEL,LABEL,LABEL,
     *                EXC,DUMMY,DUMMY,ISY,ISYMX,1,IXNR)

C
C----------------------
C        Write out result.
C----------------------
C
         WRITE(LUPRI,*) ' '
         IF (ISY.EQ.1) WRITE(LUPRI,141) LABEL//':', PROP
         IF (ISY.NE.1) WRITE(LUPRI,142) LABEL//':', 'zero by symmetry'
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) ' '
C
  141    FORMAT(20X,A9,1X,F10.6)
  142    FORMAT(20X,A9,1X,A)
C
  140 CONTINUE
C
C------------------------------------------------------------------------
C     End of loop for write out of excited state one-electron properties.
C------------------------------------------------------------------------
C
 3000 CONTINUE
C
      CALL QEXIT('CC_EXGR')
      RETURN
      END
C  /* Deck ccx_d1ao */
      SUBROUTINE CCX_D1AO(AODEN,WORK,LWORK,MODEL,
     *                    LLIST,ILLNR,RLIST,ILRNR)
C
C     Ove Christiansen April 1997 inspired by CC_D1AO
C
C     Purpose: To calculate contributions to the excited state
C              one electron density matrix and return it backtransformed
C              to AO-basis in AODEN.
C
C     Current models: CCS, CC2 and CCSD    
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK)
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
#include <ccnoddy.h>
C
      CHARACTER*8 FNTOC
      CHARACTER*4 FNDKBC
      CHARACTER*5 FNDKBC3
      CHARACTER*6 FNDELD, FNCKJD, FN3VI, FN3FOPX
      CHARACTER*7 FN3FOP2X
 
      PARAMETER (FNDELD  = 'CKDELD'  , FNCKJD  = 'CKJDEL'  )
      PARAMETER (FNDKBC  = 'DKBC'    , FNTOC   = 'CCSDT_OC')
      PARAMETER (FN3VI   = 'CC3_VI'  )
      PARAMETER (FN3FOPX = 'PTFOPX'  , FN3FOP2X= 'PTFOP2X' )
      PARAMETER (FNDKBC3 = 'DKBC3'   )                      
C
      PARAMETER ( ISYM0 = 1 )
      CHARACTER MODEL*10,MODDUM*10
      CHARACTER LLIST*(*),RLIST*(*)
C
      CALL QENTER('CCX_D1AO')
C
C---------------------------------------------
C     Find symmetry of left and right vectors. 
C---------------------------------------------
C
      ISYMR = ISYEXC(ILRNR)
      ISYML = ISYEXC(ILLNR)
      IF (ISYMR .NE. ISYML) 
     &     CALL QUIT('CCX_D1AO: Density not total sym.')
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
      KL2AM  = KL1AM  + NT1AM(ISYML)
      KT1AM  = KL2AM  + NT2SQ(ISYML)
      KR1AM  = KT1AM  + NT1AM(ISYMOP)
      KEND0  = KR1AM  + NT1AM(ISYMR)
C
      KR2AM  = KEND0
      KR2AMT = KR2AM  + NT2AM(ISYMR)
      KEND1  = KR2AMT + NT2AM(ISYMR)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
        WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
        CALL QUIT('Insufficient memory for initial '//
     &       'allocation in CCX_D1AO')
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
C
C----------------------
C     Read left vector.
C----------------------
C
      IOPT = 3
      CALL CC_RDRSP(LLIST,ILLNR,ISYML,IOPT,MODDUM,
     *              WORK(KL1AM),WORK(KR2AM))
C
C--------------------------------
C     Square up zeta2 amplitudes.
C--------------------------------
C
      CALL CC_T2SQ(WORK(KR2AM),WORK(KL2AM),ISYML)
C
C----------------------
C     Read rigth vector.
C----------------------
C
      IOPT = 3
      CALL CC_RDRSP(RLIST,ILRNR,ISYMR,IOPT,MODDUM,
     *              WORK(KR1AM),WORK(KR2AM))
      IF (ISYMR.EQ.1) CALL CCLR_DIASCL(WORK(KR2AM),TWO,ISYMR)
C
C---------------------------------------------------
C     Read zero'th order cluster singles amplitudes.
C---------------------------------------------------
C
      IOPT = 1
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
      IF (.NOT. MP2) THEN
         CALL DCOPY(NT2AM(ISYMR),WORK(KR2AM),1,WORK(KR2AMT),1)
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KR2AMT),1.0D0,ISYMR,IOPTTCME)
      ENDIF
C
C---------------------
C     Test amplitudes.
C---------------------
C
      IF ( DEBUG ) THEN
         XLV1 = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
         XLV2 = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
         WRITE(LUPRI,1) 'Norm of Response vector: L1AM    ',XLV1
         WRITE(LUPRI,1) 'Norm of Response vector: L2AM    ',XLV2
         XLV1 = DDOT(NT1AM(ISYML),WORK(KR1AM),1,WORK(KR1AM),1)
         XLV2 = DDOT(NT2AM(ISYML),WORK(KR2AM),1,WORK(KR2AM),1)
         WRITE(LUPRI,1) 'Norm of Response vector: R1AM    ',XLV1
         WRITE(LUPRI,1) 'Norm of Response vector: R2AM    ',XLV2
         XLV2 = DDOT(NT2AM(ISYML),WORK(KR2AMT),1,WORK(KR2AMT),1)
         WRITE(LUPRI,1) 'Norm of Response vector: R2AM-Tr.',XLV2
         XLV1 = DDOT(NT1AM(ISYMOP),WORK(KT1AM),1,WORK(KT1AM),1)
         WRITE(LUPRI,1) 'Norm of referencevector: T1AM    ',XLV1
      ENDIF
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ai)=ZETA(ai) &
C     D(ia) is stored transposed
C-------------------------------
C
      KXMAT  = KEND1
      KYMAT  = KXMAT  + NMATIJ(1)
      KEND2  = KYMAT  + NMATAB(1)
      LWRK2  = LWORK  - KEND1
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insufficient memory for 2. alloc. in CCX_D1AO')
      ENDIF
C
C-----------------------------------------------------------
C     Calculate X-intermediate of L2AM- and R2AM-amplitudes.
C-----------------------------------------------------------
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
C     DXai is zero.
C--------------------------------------------------------------
C     Construct <L2|[Emn,R2]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND2),LWRK2,1)
      CALL DAXPY(NMATIJ(1),-ONE,WORK(KXMAT),1,WORK(KONEIJ),1)
C
C--------------------------------------------------------------
C     Construct <L1|[Emn,R1]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
C
      CALL CC_DXIJ(WORK(KONEIJ),WORK(KL1AM),ISYML,WORK(KR1AM),ISYMR)
      CALL CC_DXAB(WORK(KONEAB),WORK(KL1AM),ISYML,WORK(KR1AM),ISYMR)
C
C--------------------------------------------------------------
C     Construct <L1|[Eia,R2]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
C
      CALL CC_DXIA12(WORK(KONEIA),WORK(KL1AM),ISYML,WORK(KR2AMT),ISYMR)
C
C----------------------------------------------------------
C     Construct <L2|[[Eia,R1],T2]|HF> contribution to DXia.
C----------------------------------------------------------
C
      KT2AM  = KEND0
      KEND3  = KT2AM  + NT2AM(ISYMOP) 
      LWRK3  = LWORK  - KEND3 
C
      IF (LWRK3 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
         CALL QUIT('Insufficient memory for 3. alloc. in CCX_D1AO')
      ENDIF
      CALL CC_DXIA21(WORK(KONEIA),WORK(KL2AM),ISYML,
     *               WORK(KR1AM),ISYMR,WORK(KT2AM),ISYMOP,
     *               WORK(KEND3),LWRK3)
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
C     Calculate the lamda matrices.
C----------------------------------
C
      KLAMDP = KEND0
      KLAMDH = KLAMDP + NLAMDT
      KEND4  = KLAMDH + NLAMDT
      LWRK4  = LWORK  - KEND4 
      IF (LWRK4 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND4
         CALL QUIT('Insufficient memory for 4. alloc. in CCX_D1AO')
      ENDIF
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND4),
     *            LWRK4)

*======================================================================*
*     Get the triples contributions and add them:
*======================================================================*
 
      IF (CCSDT) THEN

        KFOCK0 = KEND4
        KEND5  = KFOCK0 + N2BST(ISYM0)

        KDAB   = KEND5
        KDIJ   = KDAB   + NMATAB(ISYM0)
        KDIA   = KDIJ   + NMATIJ(ISYM0)
        KEND5  = KDIA   + NT1AM(ISYM0)  

C       --------------------------------
C       Get SCF Fock matrix in MO basis:
C       --------------------------------
        LWRK5  = LWORK  - KEND5 
        IF (LWRK5 .LT. 0) THEN
          WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND5
          CALL QUIT('Insufficient memory for 5. alloc. in CCX_D1AO')
        ENDIF

        ! read zeroth-order AO Fock matrix from file:
        LUFOCK = 0
        CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        REWIND(LUFOCK)
        READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
        CALL GPCLOSE(LUFOCK,'KEEP')
 
        CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMDP),WORK(KLAMDH),
     &                WORK(KEND5),LWRK5,1,1,1) 

C       -------------------------
C       Open some integral files:
C       -------------------------
        LUDELD  = 0
        LUCKJD  = 0
        LUDKBC  = 0
        LUTOC   = 0
        LU3VI   = 0
        LUDKBC3 = 0
        LU3FOPX = 0
        LU3FOP2X= 0
   
        CALL WOPEN2(LUDELD,FNDELD,64,0)
        CALL WOPEN2(LUCKJD,FNCKJD,64,0)
        CALL WOPEN2(LUDKBC,FNDKBC,64,0)
        CALL WOPEN2(LUTOC,FNTOC,64,0)
        CALL WOPEN2(LU3VI,FN3VI,64,0)
        CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
        CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
        CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)

C       ----------------------------------------------------------
C       Get the triples contributions:
C       ----------------------------------------------------------
        IF (NODDY_DEN) THEN

          CALL CCSDT_ADEN_NODDY(LLIST,ILLNR,RLIST,ILRNR,
     *                    WORK(KLAMDP),WORK(KLAMDH),WORK(KFOCK0),
     *                    WORK(KDIJ),WORK(KDAB),.TRUE.,WORK(KDIA),
     *                    .FALSE.,DUMMY,WORK(KEND5),LWRK5,
     *                    LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                    FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                    LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                    LU3FOP2X,FN3FOP2X)

        ELSE

          CALL QUIT('call to CC3_ADEN(?) missing in CCX_D1AO.')

        END IF

C       ---------------------------------------
C       Add triples contributions to densities:
C       ---------------------------------------
        CALL DAXPY(NMATAB(ISYM0),ONE,WORK(KDAB),1,WORK(KONEAB),1)
        CALL DAXPY(NMATIJ(ISYM0),ONE,WORK(KDIJ),1,WORK(KONEIJ),1)
        CALL DAXPY(NT1AM(ISYM0), ONE,WORK(KDIA),1,WORK(KONEIA),1)

C       ---------------------
C       Close integral files:
C       ---------------------
        CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
        CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
        CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
        CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
        CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
        CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
        CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
        CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')  

      END IF
C
C--------------------------------------------------------
C     Add the one electron density in the AO-basis.
C--------------------------------------------------------
C
      ISDEN = 1
      CALL CC_DENAO(AODEN,ISDEN,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KEND4),LWRK4)
C
      CALL QEXIT('CCX_D1AO')
C
   1  FORMAT(1x,A35,1X,E20.10)
      RETURN
      END
C  /* Deck cc_dxij */
      SUBROUTINE CC_DXIJ(DIJ,XL1AM,ISYML,XR1AM,ISYMR) 
C
C     Written by Ove Christiansen April 1997
C
C     Purpose: To add contributions from the L1 and R1 -amplitudes to
C              the excited state CC one electron density in MO-basis.
C              <L1|[Eij,R1]|HF>  
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, XMONE = -1.0D0, ONE = 1.0D00)
      DIMENSION DIJ(*), XL1AM(*), XR1AM(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_DXIJ')
C
C---------------------------------------------------------
C     Add contribution.
C     Assume ISYMR = ISYML and density is total symmetric.
C---------------------------------------------------------
C
      ISYRES = MULD2H(ISYMR,ISYML)
      IF (ISYRES .NE. ISYMOP) CALL QUIT('CC_DXIJ require tot.sym Dens.')
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ  = MULD2H(ISYRES,ISYMI)
         ISYMC  = MULD2H(ISYMR,ISYMI)
C
         KOFF1  = IT1AM(ISYMC,ISYMI)  + 1
         KOFF2  = IT1AM(ISYMC,ISYMJ)  + 1
         KOFF3  = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),
     *              XMONE,XR1AM(KOFF1),NTOTC,XL1AM(KOFF2),NTOTC,ONE,
     *              DIJ(KOFF3),NTOTI)
C
  100 CONTINUE
C
      CALL QEXIT('CC_DXIJ')
C
      END 
C  /* Deck cc_dxab */
      SUBROUTINE CC_DXAB(DAB,XL1AM,ISYML,XR1AM,ISYMR)
C
C     Written by Ove Christiansen April 1997
C
C     Purpose: To add contributions from the L1 and R1 -amplitudes to
C              the excited state CC one electron density in MO-basis.
C              <L1|[Eab,R1]|HF>
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, XMONE = -1.0D0, ONE = 1.0D00)
      DIMENSION DAB(*), XL1AM(*), XR1AM(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_DXAB')
C
C---------------------------------------------------------
C     Add contribution.
C     Assume ISYMR = ISYML and density is total symmetric.
C---------------------------------------------------------
C
      ISYRES = MULD2H(ISYMR,ISYML)
      IF (ISYRES .NE. ISYMOP) CALL QUIT('CC_DXAB require tot.sym Dens.')
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMB  = MULD2H(ISYMA,ISYRES)
         ISYMK  = MULD2H(ISYMA,ISYMR)
C
         KOFF1 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF2 = IT1AM(ISYMB,ISYMK)  + 1
         KOFF3 = IMATAB(ISYMA,ISYMB) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),
     *              ONE,XL1AM(KOFF1),NTOTA,XR1AM(KOFF2),NTOTB,ONE, 
     *              DAB(KOFF3),NTOTA)
C
  100 CONTINUE
C
      CALL QEXIT('CC_DXAB')
C
      END
C  /* Deck cc_dxai12 */
      SUBROUTINE CC_DXIA12(DIA,XL1AM,ISYML,XR2AM,ISYMR)
C
C     Written by Ove Christiansen April 1997
C
C     Purpose: To add contributions from the L1 and R2 -amplitudes to
C              the excited state CC one electron density in MO-basis.
C              <L1|[Eia,R2]|HF>
C              NB - 2c-E XR2AM assumed 
C              and the Dia block is stored transposed, i.e. like a t1-amplitude!
C
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
      PARAMETER(ZERO = 0.0D0, XMONE = -1.0D0, ONE = 1.0D00)
      DIMENSION DIA(*), XL1AM(*), XR2AM(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CC_DXIA12')
C
C---------------------------------------------------------
C     Assume ISYMR = ISYML and density is total symmetric.
C---------------------------------------------------------
C
      ISYRES = MULD2H(ISYMR,ISYML)
      IF (ISYRES .NE. ISYMOP) CALL QUIT('CC_DXIA12 require '//
     &     'tot.sym Dens.')
C
      ISYMAI = ISYRES
      ISYMDK = MULD2H(ISYMR,ISYMAI)
C
C--------------------------
C     Set up density block.
C--------------------------
C
      DO 100 NAI = 1,NT1AM(ISYMAI)
         DO 110 NDK = 1,NT1AM(ISYMDK)
C
            IF (ISYMAI .EQ. ISYMDK) THEN
               NDKAI = IT2AM(ISYMDK,ISYMAI) + INDEX(NDK,NAI)
            ELSE IF (ISYMDK .LT. ISYMAI) THEN
               NDKAI = IT2AM(ISYMDK,ISYMAI) 
     *               + NT1AM(ISYMDK)*(NAI-1) + NDK
            ELSE IF (ISYMDK .GT. ISYMAI) THEN
               NDKAI = IT2AM(ISYMDK,ISYMAI) 
     *               + NT1AM(ISYMAI)*(NDK-1) + NAI
            ENDIF
            DIA(NAI) = DIA(NAI) + XR2AM(NDKAI)*XL1AM(NDK)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_DXIA12')
C
      RETURN
      END
C  /* Deck cc_dxia21 */
      SUBROUTINE CC_DXIA21(DIA,XL2AM,ISYML,XR1AM,ISYMR,
     *                     T2AM,ISYMT2,WORK,LWORK)
C
C     Written by Ove Christiansen April 1997
C
C     Purpose: Construct <L2|[[Eia,R1],T2]|HF> contribution to DXia, 
C              the excited state CC one electron density in MO-basis.
C              The Dia block is stored transposed, i.e. like a t1-amplitude.
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      PARAMETER(ZERO = 0.0D0, XMONE = -1.0D0, ONE = 1.0D00)
      DIMENSION DIA(*), XL2AM(*), XR1AM(*), T2AM(*),WORK(LWORK)
      CHARACTER*10 MODEL
C
      CALL QENTER('CC_DXIA21')
C
C---------------------------------------------------------
C     Assume ISYMR = ISYML and density is total symmetric.
C---------------------------------------------------------
C
      ISYRES = MULD2H(ISYMR,ISYML)
      IF (ISYRES .NE. ISYMOP)
     &     CALL QUIT('CC_DXIA21 require tot.sym Dens.')
      IF (ISYRES .NE. ISYMT2)
     &     CALL QUIT('symmetry mismatch in CC_DXIA21.')
C
C---------------------------------------------------
C     Read zero'th order cluster doubles amplitudes.
C---------------------------------------------------
C
      IOPT = 2
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,T2AM)
C
C----------------------------------------------------------
C     Construct <L2|[[Eia,R1],T2]|HF> contribution to DXia.
C----------------------------------------------------------
C
      KONEAB = 1
      KONEIJ = KONEAB + NMATAB(ISYML)
      KXMAT  = KONEIJ + NMATIJ(ISYML)
      KYMAT  = KXMAT  + NMATIJ(ISYML)
      KEND1  = KYMAT  + NMATAB(ISYML)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 1. alloc. in CC_DXIA21.')
      ENDIF
C
C-----------------------------------------------------------
C     Calculate X-intermediate of L2AM- and T2AM-amplitudes.
C-----------------------------------------------------------
C
         XLV2 = DDOT(NT2AM(ISYMT2),T2AM,1,T2AM,1)
         WRITE(LUPRI,*) 'Norm of Response vector: T2AM   .',XLV2
         CALL CC_XI(WORK(KXMAT),XL2AM,ISYML,T2AM,ISYMT2,
     *           WORK(KEND1),LWRK1)
         XIJ = DDOT(NMATIJ(ISYML),WORK(KXMAT),1,WORK(KXMAT),1)
         WRITE(LUPRI,*) 'Norms: XMAT',XIJ
C
C-----------------------------------------------------------
C     Calculate Y-intermediate of L2AM- and T2AM-amplitudes.
C-----------------------------------------------------------
C
         CALL CC_YI(WORK(KYMAT),XL2AM,ISYML,T2AM,ISYMT2,
     *           WORK(KEND1),LWRK1)
C
      CALL DZERO(WORK(KONEIJ),NMATIJ(ISYML))
      CALL DCOPY(NMATAB(ISYML),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND1),LWRK1,ISYML)
      CALL DAXPY(NMATIJ(ISYML),XMONE,WORK(KXMAT),1,WORK(KONEIJ),1)
      IF (DEBUG) THEN
         XIJ = DDOT(NMATIJ(ISYML),WORK(KONEIJ),1,WORK(KONEIJ),1)
         XAB = DDOT(NMATAB(ISYML),WORK(KONEAB),1,WORK(KONEAB),1)
         WRITE(LUPRI,*) 'Norms: DXIJ',XIJ
         WRITE(LUPRI,*) 'Norms: DXAB',XAB
      ENDIF
C
      CALL CC_IA21(DIA,WORK(KONEAB),WORK(KONEIJ),ISYML,XR1AM,ISYMR,
     *             WORK(KEND1),LWRK1)
C
      CALL QEXIT('CC_DXIA21')
C
      END
      SUBROUTINE CC_IA21(DIA,DAB,DIJ,ISYMD,XR1AM,ISYMR,WORK,LWORK)
C
C     Written by Ove Christiansen April 1997 
C
C     Version: 1.0
C
C     Purpose: Contract DAB and DIJ with R1AM into DIA
C              The Dia block is stored transposed, i.e. like a t1-amplitude!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, XMONE = -1.0D0)
      DIMENSION DIA(*),DAB(*),DIJ(*),XR1AM(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_IA21')
C
C-----------------------------------------------------
C     Assume ISYMR = ISYMD and Dia is total symmetric.
C-----------------------------------------------------
C
      ISYRES = MULD2H(ISYMR,ISYMD)
      IF (ISYRES .NE. ISYMOP) CALL QUIT('CC_IA21 require tot.sym Dens.')
      ISYMAI = ISYRES
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMA,ISYMAI)
         ISYMK = MULD2H(ISYMA,ISYMR)
C
         KOFF1 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF2 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTI = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *              ONE,XR1AM(KOFF1),NTOTA,DIJ(KOFF2),NTOTI,ONE,
     *              DIA(KOFF3),NTOTA)
C
         ISYMC = MULD2H(ISYMI,ISYMR)
C
         KOFF1 = IMATAB(ISYMC,ISYMA) + 1
         KOFF2 = IT1AM(ISYMC,ISYMI)  + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTI = MAX(NRHF(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),
     *              XMONE,DAB(KOFF1),NTOTC,XR1AM(KOFF2),NTOTC,ONE,
     *              DIA(KOFF3),NTOTA)
C
  100 CONTINUE
C
      CALL QEXIT('CC_IA21')
C
      END
      SUBROUTINE CCSX_D1AO(AODEN,WORK,LWORK, 
     *                     LLIST,ILLNR,RLIST,ILRNR,L0LIST,ILNRL0)
C
C     Ove Christiansen April 1997 inspired by CC_D1AO
C
C     Purpose: To calculate contributions to the excited state
C              one electron density matrix and return it backtransformed
C              to AO-basis in AODEN.
C           <LE1|[Emn,RE1]|HF> contribution. 
C           For CCS but not CIS also;
C           <L1|Emn|HF> 
C
C     Current models: CCS, CC2 and CCSD
C
#include <implicit.h>      
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK)
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
      CHARACTER MODEL*5,MODDUM*10
      CHARACTER LLIST*(*),RLIST*(*)
C
      CALL QENTER('CCSX_D1AO')
C
C---------------------------------------------
C     Find symmetry of left and right vectors.
C---------------------------------------------
C
      ISYMR = ISYEXC(ILRNR)
      ISYML = ISYEXC(ILLNR)
      IF (ISYMR .NE. ISYML)
     &     CALL QUIT('CCSX_D1AO: Density not total sym.')
C
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
      KR1AM  = KL1AM  + NT1AM(ISYML)
      KT1AM  = KR1AM  + NT1AM(ISYMR)
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
      LWRK1  = LWORK  - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
        WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
        CALL QUIT('Insufficient memory for initial '//
     &       'allocation in CCSX_D1AO')
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
C
C-----------------------------
C     Read Left  eigen-vector.
C-----------------------------
C
      IOPT = 1
      CALL CC_RDRSP(LLIST,ILLNR,ISYML,IOPT,MODDUM,
     *              WORK(KL1AM),WORK(KEND1))
C
C-----------------------------
C     Read rigth eigen-vector.
C-----------------------------
C
      IOPT = 1
      CALL CC_RDRSP(RLIST,ILRNR,ISYMR,IOPT,MODDUM,
     *              WORK(KR1AM),WORK(KEND1))
C
C--------------------------------------------------------------
C     Construct <L1|[Emn,R1]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
C
      CALL CC_DXIJ(WORK(KONEIJ),WORK(KL1AM),ISYML,WORK(KR1AM),ISYMR)
      CALL CC_DXAB(WORK(KONEAB),WORK(KL1AM),ISYML,WORK(KR1AM),ISYMR)
C
C-------------------------------
C     <L1|Emn|HF>  Contribution
C-------------------------------
C
      IF (CCS.AND.(.NOT.CIS)) THEN
         IOPT = 1
         CALL CC_RDRSP(L0LIST,ILNRL0,ISYMOP,IOPT,MODDUM,
     *                 WORK(KONEAI),WORK(KEND1))
      ENDIF 
C
C-------------------------------
C     Get MO coefficient matrix.
C-------------------------------
C
      CALL DZERO(WORK(KT1AM),NT1AMX)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
C--------------------------------------------------------
C     Add the one electron density in the AO-basis.
C--------------------------------------------------------
C
      CALL DZERO(AODEN,N2BST(1))
      ISDEN = 1
      CALL CC_DENAO(AODEN,ISDEN,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
      CALL QEXIT('CCSX_D1AO')
C
      END
c*DECK TNSRAN
      SUBROUTINE CC_TNSRAN(TNSR,WORK,LWORK)
C
C------------------------------------------------------------------------
C
C     Call TNSRAN and write out selected info.
C
C------------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <ccsdinp.h>
C
      PARAMETER (THR = 1.0D-08)
      DIMENSION TNSR(3,3),PVAL(3),PAXIS(3,3)
C
      CALL QENTER('CC_TNSRAN')
C
      CALL TNSRAN(TNSR,PVAL,PAXIS,ALFSQ,BETSQ,ITST,ITST2,
     *            APAR,APEN,XKAPPA,IPAR)
C
      WRITE(LUPRI,'(/,1X,A38,F14.6)')
     *              'Alfa**2 Invariant:            '
     *            //'            ',ALFSQ
      WRITE(LUPRI,'(1X,A38,F14.6)')
     *           'Beta**2 Invariant:            '
     *            //'            ',BETSQ
      SHPAL = SQRT(ALFSQ)
      ANINV = SQRT(BETSQ)
      WRITE(LUPRI,'(/,1X,A42,F10.6,A)') 'Isotropic Property:       '
     *         //'                 ',SHPAL,' a.u.'
      WRITE(LUPRI,'(1X,A42,F10.6,A)') 'Property anisotropy invariant:'
     *      //'            ',ANINV,' a.u.'

      CALL QEXIT('CC_TNSRAN')
C
      END 

