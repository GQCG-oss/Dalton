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
C
C---------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
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
C
      LOGICAL LINQCC, LPROJECT, TRIPLET, LREDS,B0SKIP
      DIMENSION WORK(LWORK), ELSEMO(3,3), SKODE(3,3), SKODN(3,3)
      CHARACTER*8 FC1AM,FC2AM,FRHO1,FRHO2,FRHO12,FC12AM,FS12AM,FS2AM
      CHARACTER*3 APROXR12
      PARAMETER (FC1AM ='CCR_C1AM',FC2AM ='CCR_C2AM')
      PARAMETER (FRHO1 ='CCR_RHO1',FRHO2 ='CCR_RHO2')
      PARAMETER (FRHO12='CCR_RH12',FC12AM='CCR_C12M',FS12AM='CCR_S12M')
      PARAMETER (FS2AM ='CCR_S2DM')
      CHARACTER MODEL*10
      CHARACTER MODELPRI*4, MODELPRI2*30
      CHARACTER LABEL*8, LABR12*4
      CHARACTER*2 LIST
      CHARACTER*5 FNDPTIA, FNDPTAB, FNDPTIJ
      CHARACTER*6 FNDPTIA2
      PARAMETER(FNDPTIA='DPTIA', FNDPTIA2 = 'DPTIA2',
     *          FNDPTAB='DPTAB' ,FNDPTIJ  = 'DPTIJ'  ) 
C
#include <leinf.h>
C
      CALL QENTER('CC_EXGR')
C
C------------------------------------
C     Header of Property calculation.
C------------------------------------
C
C     CALL TIMER('START ',TIMEIN,TIMOUT)
C
      LUFC1  = -1
      LUFC2  = -1
      LUFC12 = -1
      LUFR1  = -1
      LUFR2  = -1
      LUFR12 = -1
      LUFS12 = -1
      LUFS2  = -1

      IF (CCR12) THEN
         LABR12 = '-R12'
      ELSE
         LABR12 = '    '
      END IF
C
      IF (CCR12) CALL QUIT('R12 not yet implemented in CC_EXGR.')
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
      MODEL = 'CCSD      '
      IF (CC2) THEN
         CALL AROUND('Coupled Cluster model is: CC2'//LABR12)
         MODEL = 'CC2       '
         MODELPRI = ' CC2'
      ENDIF
      IF (CCS.AND.(.NOT.CIS)) THEN
         CALL AROUND('Coupled Cluster model is: CCS')
         MODEL = 'CCS       '
         MODELPRI = ' CCS'
      ENDIF
      IF (CCS.AND.CIS) THEN
         CALL AROUND('Model is not CC but CIS crap.')
         MODEL = 'CCS       '
         MODELPRI = ' CIS'
      ENDIF
      IF (CCD) THEN
         CALL AROUND('Coupled Cluster model is: CCD'//LABR12)
         MODEL = 'CCD       '
         MODELPRI = ' CCD'
      ENDIF
      IF (CC3  ) THEN
         CALL AROUND('Coupled Cluster model is: CC3'//LABR12)
         MODEL = 'CC3       '
         MODELPRI = ' CC3'
         WRITE(LUPRI,*)
     *    'CC3 X-state first order properties not implemented'
         RETURN
      ENDIF
      IF (CC1A) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1a'//LABR12)
         MODEL = 'CCSDT-1a  '
         WRITE(LUPRI,*)
     *    'CCSDT-1a X-state first order properties not implemented'
         RETURN
      ENDIF
      IF (CC1B) THEN
         CALL AROUND('Coupled Cluster model is: CCSDT-1b'//LABR12)
         MODEL = 'CCSDT-1b  '
         WRITE(LUPRI,*)
     *    'CCSDT-1b X-state first order properties not implemented'
         RETURN
      ENDIF
      IF (CCSD) THEN
         CALL AROUND('Coupled Cluster model is: CCSD'//LABR12)
         MODEL = 'CCSD      '
         MODELPRI = 'CCSD'
      ENDIF
C
      NSIDIN = NSIDE
C
      LREDS = CCR12
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_EXGR-1: Workspace:',LWORK
C
C-----------------------------
C     Initialize Variables.
C-----------------------------
C
      ISYM = ISYMOP
      ISYMTR = ISYMOP
      IF (CCS) THEN
        NCCVAR = NT1AM(ISYM)
      ELSE      
        NCCVAR = NT1AM(ISYM) + NT2AM(ISYM)
        IF (CCR12) NCCVAR = NCCVAR + NTR12AM(ISYM)
      END IF
      TRIPLET = .FALSE.
      CALL CCLR_LEINFI(TRIPLET)
      THRLE = THRLEQ
      LINQCC = .TRUE.
      NLOAD  = 1
      ISIDE  = -1
      LIST  = 'E0'
      IF (NCCVAR .LE. 0) THEN
         WRITE(LUPRI,*) 'There are no amplitudes of this symmetry'
         CALL QUIT('There are no amplitudes of this symmetry')
      ENDIF
C
C-------------------------------------------
C     For CIS jump solution for multipliers.
C-------------------------------------------
C
      IF (CCS.AND.CIS) GOTO 47
C
C-----------------------------------------------
C     Calculate appropriate F/B-transformations.
C-----------------------------------------------
C
c     IF ((.NOT.E0SKIP).AND.(.NOT.B0SKIP)) THEN
c hjaaj jan-07: B0SKIP is not defined
      IF (.NOT.E0SKIP) THEN
         CALL CC_FRE(WORK,LWORK)
      ENDIF
C
C-----------------------------------------------------------------
C     Start loop over lists over right hand lists to be solved.
C     Find out how many from common block info.
C-----------------------------------------------------------------
C
      IBOFF = 0
      IBEND = NXGRST
      NEQSYM  = IBEND - IBOFF
      IF (E0SKIP) THEN
         NEQSYM = 0
         WRITE(LUPRI,*) 'Skipping solving for excited state t(0)'
      ENDIF
C
C------------------------------------------------------------------
C     Make bathching over nr. of simulatneous vectors to be solved.
C     Find out how many simultaneously;
C     all or take the number from input.
C------------------------------------------------------------------
C
      IF (NEQSYM .EQ. 0 ) THEN
         NSIM = 0
         NBAT = 0
      ELSE
         IF ( NSIMLE .EQ. 0 ) THEN
            NSIM  = NEQSYM
         ELSE
            NSIM  = NSIMLE
         ENDIF
         NBAT  = (NEQSYM-1)/NSIM + 1
      ENDIF
C
      IF (DEBUG ) THEN
         WRITE(LUPRI,*) '         '
         WRITE(LUPRI,*) 'CC_EXGR      NEQSYM = ',NEQSYM
         WRITE(LUPRI,*) 'CC_EXGR      IBOFF  = ',IBOFF
         WRITE(LUPRI,*) 'CC_EXGR      IBEND  = ',IBEND
         WRITE(LUPRI,*) 'CC_EXGR      NSIM   = ',NSIM
         WRITE(LUPRI,*) 'CC_EXGR      NBAT   = ',NBAT
      ENDIF
C
      DO 2000 IRB = 1, NBAT
C
C-------------------------------------------------
C        Find start number for this batch on list.
C        and the number of equations.
C-------------------------------------------------
C
         IRST  = IBOFF + (IRB-1)*NSIM + 1
         NSIMR = MIN(NSIM,NEQSYM - (IRB-1)*NSIM)
         IREND = IRST  + NSIMR - 1
C
         IF (DEBUG ) THEN
            WRITE(LUPRI,*) 'CC_EXGR      IRST   = ',IRST
            WRITE(LUPRI,*) 'CC_EXGR      IREND  = ',IREND
            WRITE(LUPRI,*) 'CC_EXGR      NSIMR  = ',NSIMR
         ENDIF
C
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(/,1x,A,I2,A)')
     *        'Solving ',NSIMR,' response equations for states;'
            DO 2005 IR = IRST,IREND
                 IEX = IXGRST(IR)
                 WRITE(LUPRI,'(1X,A,F10.6,A,I3,A,I3)')
     *             'Energy:',EIGVAL(IEX),
     *           ' and nr.: ',IEX-ISYOFE(ISYEXC(IEX)),' of symmetry: ',
     *           ISYEXC(IEX)
 2005       CONTINUE
         ENDIF
C
         IF (IPRINT.GT.10)
     *         WRITE(LUPRI,*) 'CC_EXGR      Workspace:',LWORK
C
         CALL FLSHFO(LUPRI)
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         KIPLAC = 1
         KREDH  = KIPLAC + MAXRED
         KREDGD = KREDH  + MAXRED*MAXRED
         KEIVAL = KREDGD + MAXRED
         KSOLEQ = KEIVAL + MAXRED
         KWRK1  = KSOLEQ + MAXRED*MAXRED
         IF (LREDS) THEN
           KREDS = KWRK1
           KWRK1 = KREDS + MAXRED*MAXRED
         END IF
         LWRK1  = LWORK  - KWRK1
C
         IF (LWRK1.LT. 0 )
     *         CALL QUIT(' TOO LITTLE WORKSPACE IN CC_RSPSOL')
C
C-------------------------
C        Zero frequencies.
C-------------------------
C
         CALL DZERO(WORK(KEIVAL),MAXRED)
C
C-------------------
C        Open files.
C-------------------
C
         CALL CC_FILOP(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2)
C
C-----------------------------
C        Create start vectors.
C-----------------------------
C
         TRIPLET  = .FALSE.
         LPROJECT = .FALSE.
         CALL CCEQ_STR(FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2,LPROJECT,ISTATPRJ,
     *                 TRIPLET,ISIDE,IRST,NSIMR,NUPVEC,
     *                 NREDH,WORK(KEIVAL),WORK(KIPLAC),
     *                 WORK(KWRK1),LWRK1,LIST)
C
C------------------------------------------
C        Solve equations by call to solver.
C------------------------------------------
C
         write(LUPRI,*) 'exgr, irst:',irst
         ECURR = 0.0D0
         CALL CCEQ_SOL(LIST,LPROJECT,ISTATPRJ,ECURR,
     *                 FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 LREDS,WORK(KREDS),FS12AM,LUFS12,FS2AM,LUFS2,
     *                 LINQCC,TRIPLET,ISIDE,IRST,NSIMR,NUPVEC,
     *                 NREDH,WORK(KREDH),WORK(KEIVAL),
     *                 WORK(KSOLEQ),WORK(KWRK1),LWRK1,APROXR12)
C
C----------------------------------
C        Analysis of solution vectors.
C----------------------------------
C
         NVARPT= LETOT + 2*NALLAI(ISYM)
         KWRK2 = KWRK1 + NVARPT
         LWRK2 = LWORK - KWRK2
C
         THRESH = 0.05
         MAXLIN = 100
         NSIMUL = MIN(NSIM,LWRK2/LETOT )
         NBATCH = (NSIM-1)/NSIMUL + 1
         IOFF1  = 1
         ICOUNT = 0
         ILSTNR = IRST
C
         DO 500 I = 1,NBATCH
            IOFF2 = MIN(NSIMUL,NSIM - (I-1)*NSIMUL)
            CALL CCCONV(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *                  TRIPLET,CCR12,NREDH,IOFF1,IOFF2,WORK(KSOLEQ),
     *                  WORK(KWRK2),WORK(KWRK1))
C
            IF ( IPRINT .GT. 10 ) THEN
               RHO1N = DDOT(NT1AM(ISYM),WORK(KWRK2),1,WORK(KWRK2),1)
               RHO2N = DDOT(NT2AM(ISYM),WORK(KWRK2+NT1AMX),1,
     *                      WORK(KWRK2+NT1AMX),1)
               WRITE(LUPRI,*) 'Norm of Lambda vector :',RHO1N
               WRITE(LUPRI,*) 'Norm of Lambda vector :',RHO2N
            ENDIF
C
            IF ( IPRINT .GT. 30 ) THEN
               CALL AROUND('CC_EXGR: E0 vector in mo basis' )
               CALL OUTPUT(WORK(KWRK2),1,LETOT,1,NSIM,
     *                     LETOT,NSIM,1,LUPRI)
            ENDIF
C
            DO 510 J = 1,IOFF2
               ICOUNT = ICOUNT + 1
               IF (IPRINT .GT. 1) THEN
                 WRITE(LUPRI,'(//1X,A)')
     *'Analysis of the X-st 0-th order Lagrangian multipliers: '
         WRITE(LUPRI,'(1X,A)')
     *'--------------------------------------------------------'
                 CALL CC_PRAM(WORK(KWRK2 + (ICOUNT-1)*LETOT),
     *                        PT1LOCAL,ISYM,.FALSE.)
               ENDIF
               IF (CCSTST) THEN
                  CALL DZERO(WORK(KWRK2+(ICOUNT-1)*LETOT+
     *                       NT1AM(ISYM)),NT2AMX)
               ENDIF
C
C-----------------------------------------
C              Save response vectors on file.
C-----------------------------------------
C
               KT1    = KWRK2 + (ICOUNT-1)*LETOT
               KT2    = KWRK2 + (ICOUNT-1)*LETOT + NT1AM(ISYM)
               IOPT   = 3
               CALL CC_WRRSP('E0',ILSTNR,ISYM,IOPT,MODEL,DUMMY,
     *                        WORK(KT1),WORK(KT2),WORK(KWRK1),NVARPT)
               ILSTNR = ILSTNR + 1
C
  510       CONTINUE
            IOFF1 = IOFF1 + NSIMUL
  500    CONTINUE
C
C-------------------------------
C        Close and delete files.
C-------------------------------
C
         CALL CC_FILCL(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2)
         CALL FLSHFO(LUPRI)
C
 2000 CONTINUE
C
C--------------------------------
C     Spaghetti goto end for CIS.
C--------------------------------
C
  47  CONTINUE
C
C-----------------------------------------------------------------------
C     Calculate one electron AO-density and thereby CC nat.occ.num.
C-----------------------------------------------------------------------
C
      DO 3000 IX = 1, NXGRST
C
         IEX   = IXGRST(IX)
         ISYMX = ISYEXC(IEX)
         IXNR  = IEX - ISYOFE(ISYMX)
         WRITE(LUPRI,*) 'Total, actual nr.,loop index',NXGRST,IEX,IX
         WRITE(LUPRI,*) 'Sym, reduced nr.',ISYMX,IXNR
         IF (.NOT.CCS) THEN
            KDENS = 1    
            KWRK2 = KDENS  + N2BST(ISYMOP)
            LWRK2 = LWORK  - KWRK2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
               CALL QUIT('Insufficient memory for one '//
     &              'e-density in CC_EXGR')
            ENDIF
C
            ILLNR = IX
C
Casm        ifort does not like LIST. Use 'E0' instead.
C
            CALL CC_D1AO(IDUMMY,.FALSE.,WORK(KDENS),DUM,WORK(KWRK2),
     &                   LWRK2,MODEL,'E0',ILLNR,.FALSE.,
     &                   FNDPTIA, FNDPTIA2, FNDPTAB, FNDPTIJ)
casm        CALL CC_D1AO(IDUMMY,.FALSE.,WORK(KDENS),DUM,WORK(KWRK2),
casm &                   LWRK2,MODEL,LIST,ILLNR,.FALSE.,
casm &                   FNDPTIA, FNDPTIA2, FNDPTAB, FNDPTIJ)
            IF (FROIMP .OR. FROEXP) THEN
              CALL CC_FCD1AO(WORK(KDENS),WORK(KWRK2),LWRK2,MODELPRI)
            ENDIF
C
            IF (DEBUG) THEN
               XD = DDOT(N2BST(ISYMOP),WORK(KDENS),1,WORK(KDENS),1)
               WRITE(LUPRI,*) 'Norm of dens-1: ',XD
            ENDIF
C
            KDENS2 = KWRK2
            KWRK3  = KDENS2 + N2BST(ISYMOP)
            LWRK3  = LWORK - KWRK3
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
               CALL QUIT('Insufficient memory for one '//
     &              'e-density in CC_EXGR')
            ENDIF
            CALL DZERO(WORK(KDENS2),N2BST(ISYMOP))
            ILLNR = IEX
            ILRNR = ILLNR
            CALL CCX_D1AO(WORK(KDENS2),WORK(KWRK3),LWRK3,MODELPRI,
     *                   'LE',ILLNR,'RE',ILRNR)
            CALL DAXPY(N2BST(ISYMOP),ONE,WORK(KDENS2),1,WORK(KDENS),1)
C
            IF (DEBUG) THEN
               XD = DDOT(N2BST(ISYMOP),WORK(KDENS),1,WORK(KDENS),1)
               WRITE(LUPRI,*) 'Norm of dens-2: ',XD
            ENDIF
C
            IF (IPRINT .GT. 50) THEN
               CALL AROUND('One electron unrelaxed density in cc_exgr')
               CALL CC_PRFCKAO(WORK(KDENS),1)
            ENDIF
C
            CALL FLSHFO(LUPRI)
         ENDIF
C
C------------------------------------------------------------------------
C        Calculate the simple one electron AO-density in CCS calculation.
C------------------------------------------------------------------------
C
         IF (CCS) THEN
C
            KDENS = 1
            KWRK3 = KDENS + N2BST(ISYMOP)
            LWRK3 = LWORK - KWRK3
C
            KDENS2= KWRK3 
            KWRK4 = KDENS2+ N2BST(ISYMOP)
            LWRK4 = LWORK - KWRK3
C
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK4
               CALL QUIT('Insufficient memory for CCS '//
     &              'AO-density in CC_EXGR')
            ENDIF
C
C------------------------------------------
C           Usual <HF|Emn|HF> contribution. 
C------------------------------------------
C
            CALL CCS_D1AO(WORK(KDENS),WORK(KWRK3),LWRK3)
            IF (FROIMP .OR. FROEXP) THEN
              CALL CC_FCD1AO(WORK(KDENS),WORK(KWRK3),LWRK3,MODELPRI)
            ENDIF
            XN = DDOT(N2BST(ISYMOP),WORK(KDENS),1,WORK(KDENS),1)
            WRITE(LUPRI,*) 'Norm of AO density-1' ,XN
C
C------------------------------------------------
C           <LE1|[Emn,RE1]|HF> contribution. 
C           For CCS but not CIS also; <L1|Emn|HF> 
C------------------------------------------------
C
            ILLNR = IEX
            ILRNR = ILLNR
            CALL CCSX_D1AO(WORK(KDENS2),WORK(KWRK4),LWRK4,
     *                    'LE',ILLNR,'RE',ILRNR,'E0',IX)
            XN = DDOT(N2BST(ISYMOP),WORK(KDENS2),1,WORK(KDENS2),1)
            WRITE(LUPRI,*) 'Norm of AO density2' ,XN
            CALL DAXPY(N2BST(ISYMOP),ONE,WORK(KDENS2),1,WORK(KDENS),1)
         ELSE
            KWRK3 = KWRK2
            LWRK3 = LWRK2
         ENDIF
C
         MODELPRI2 = 'Unrelaxed excited state '//MODELPRI
         CALL AROUND(MODELPRI2//' First-order properties: ')
C
         ISYMX = ISYEXC(IEX)
         IXNR  = IEX - ISYOFE(ISYMX)
         EXC = EIGVAL(IEX) 
         WRITE(LUPRI,'(/,1X,A,/,1X,A,I5,/,1X,A,I5,/,1X,A,F10.6,F12.6)') 
     *    'Excited state properties for ',
     *    'Excited state nr.:           ',IXNR,
     *    'Excited state sym:           ',ISYMX,
     *    'Excitation energy (au.,eV):  ',EXC,EXC*XTEV
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
         IASGER = IPRINT - 9
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
            CALL WRIPRO(DIPMN(IDIP),MODEL,-11,LABEL,LABEL,LABEL,LABEL,
     *                  EXC,DUMMY,DUMMY,ISY,ISYMX,1,IXNR)
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
            CALL WRIPRO(ELSEMO(I,J),MODEL,-11,LABEL,LABEL,LABEL,LABEL,
     *                  EXC,DUMMY,DUMMY,ISY,ISYMX,1,IXNR)
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
         CALL WRIPRO(PROP,MODEL,-11,LABEL,LABEL,LABEL,LABEL,
     *               EXC,DUMMY,DUMMY,ISY,ISYMX,1,IXNR)

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
C------------------------------------
C     Restore NSIDE to input/default.
C------------------------------------
C
      IF (.NOT. CCS) NSIDE = NSIDIN
C
      CALL QEXIT('CC_EXGR')
C
      RETURN
      END
C  /* Deck ccx_d1ao */
      SUBROUTINE CCX_D1AO(AODEN,WORK,LWORK,MODEL,
     *                    LLIST,ILLNR,RLIST,ILRNR)
C
C     Ove Christiansen April 1997 inspired by CC_D1AO
C     KS Aug 2010: minor modifications to share routines
C     with CCMM_D2AO
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
C
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
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
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
c test
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

C KS: Read amplitudes outside DXIA21 routine
      IOPT = 2
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
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
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 4. alloc. in CCX_D1AO')
      ENDIF
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND4),
     *            LWRK4)
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
C     KS, Aug 10:
C     Modified to calculate general <L2|[[Eia,R1],R2]|HF>
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
      LOGICAL LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
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
      IF (LOCDEB) THEN
        TAL1 = DDOT(NT2AMX,XL2AM,1,XL2AM,1)
        WRITE(LUPRI,*) 'CC_EXGR: Norm of L0_2: ', TAL1 
        TAL1 = DDOT(NT1AMX,XR1AM,1,XR1AM,1)
        WRITE(LUPRI,*) 'CC_EXGR: Norm of B_1: ', TAL1 
        TAL1 = DDOT(NT2AMX,T2AM,1,T2AM,1)
        WRITE(LUPRI,*) 'CC_EXGR: Norm of C_2: ', TAL1 
      END IF
C KS schedule delete next 6 lines
C since they are now read outside routine and carried
C in memory at this time
C---------------------------------------------------
C     Read zero'th order cluster doubles amplitudes.
C---------------------------------------------------
C
C      IOPT = 2
C      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,T2AM)
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
         CALL CC_XI(WORK(KXMAT),XL2AM,ISYML,T2AM,ISYMT2,
     *           WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------
C     Calculate Y-intermediate of L2AM- and T2AM-amplitudes.
C-----------------------------------------------------------
C
         CALL CC_YI(WORK(KYMAT),XL2AM,ISYML,T2AM,ISYMT2,
     *           WORK(KEND1),LWRK1)
C         TAL1 = DDOT(NMATIJ(ISYML),WORK(KXMAT),1,WORK(KXMAT),1)
C         TAL2 = DDOT(NMATAB(ISYML),WORK(KYMAT),1,WORK(KYMAT),1)
C         WRITE(LUPRI,*) 'CC_EXGR: X intermediate: ',TAL1
C         WRITE(LUPRI,*) 'CC_EXGR: Y intermediate: ',TAL2
C
      CALL DZERO(WORK(KONEIJ),NMATIJ(ISYML))
      CALL DCOPY(NMATAB(ISYML),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND1),LWRK1,ISYML)
      CALL DAXPY(NMATIJ(ISYML),XMONE,WORK(KXMAT),1,WORK(KONEIJ),1)
      IF (LOCDEB) THEN
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
C  /* Deck ccmm_dia */
      SUBROUTINE CCMM_DIA(DIA,DIJ,ISYML,XR1AM,ISYMR)
C
C     Purpose: Calculates DIA block of quadratic density to be use
C     for a quadratic response CC/MM calculation. It takes as input 
C     a pseudo-density of product coefficients \sum_a t^a_i b_p^a = Dij
C     calculated in the CC_DXIJ routine.
C             
C     KS, Aug 2010       
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, XMONE = -1.0D0, ONE = 1.0D00)
      DIMENSION DIA(*), DIJ(*), XR1AM(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCMM_DIA')
C
C---------------------------------------------------------
C     Add contribution.
C     Assume ISYMR = ISYML and density is total symmetric.
C---------------------------------------------------------
C
      ISYRES = MULD2H(ISYMR,ISYML)
      IF (ISYRES.NE.ISYMOP) CALL QUIT('CCMM_DIA require tot.sym Dens.')
C

      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMA,ISYRES)
         ISYMK = MULD2H(ISYMA,ISYMR)
C
         KOFF1 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF2 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTI = MAX(NRHF(ISYMI),1)

         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *              ONE,XR1AM(KOFF1),NTOTA,DIJ(KOFF2),NTOTI,ONE,
     *              DIA(KOFF3),NTOTA)

  100 CONTINUE
C
      CALL QEXIT('CCMM_DIA')
C
      END

