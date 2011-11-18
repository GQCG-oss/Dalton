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
c*DECK CC_QR2RSD
       SUBROUTINE CC_QR2RSD(WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              quadratic response second residue calculation.
C
C              CIS, CCS, CC2, CCSD
C
C     Ove Christiansen April 1997.
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
      PARAMETER (TOLFRQ=1.0D-08,ONE=1.0D0,XMONE=-1.0D0,THR=1.0D-08)
C
#include <iratdef.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccinftap.h>
#include <ccsdinp.h>
#include <cclrinf.h>
#include <ccexci.h>
#include <cclres.h>
#include <ccroper.h>
#include <ccqr2r.h>
C
      LOGICAL LTST,LCALC,LDIP
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*10,MODELP*10,MODEL1*10,CHSYM*2,FILEX*10,FILOLD*10
      CHARACTER LABELA*8,LABELB*8
C
#include <leinf.h>
C
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
C
C------------------------------------
C     Header of Property calculation.

C
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A)')
     *'*<<<<<<<<< OUTPUT FROM COUPLED CLUSTER QUADRATIC RESPONSE >'//
     *'>>>>>>>>*'
      IF ( XOSCST ) THEN
         WRITE (LUPRI,'(1X,A)')
     *   '*                                                        '//
     *   '         *'
         WRITE (LUPRI,'(1X,A)')
     *   '*<<<<< CALCULATION OF CC EXCITED STATE TRANSITION STRENGT'//
     *   'HS >>>>>>*'
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
         MODEL = 'CC2'
      ENDIF
      IF (CCS) THEN
         MODEL = 'CCS'
      ENDIF
      IF (CC3  ) THEN
         MODEL = 'CC3'
         WRITE(LUPRI,'(/,1x,A)') 
     *    'CC3 Oscillator strengths not implemented yet'
         RETURN
      ENDIF
      IF (CC1A) THEN
         MODEL = 'CCSDT-1a'
         WRITE(LUPRI,'(/,1x,A)') 
     *    'CC1A Oscillator strengths not implemented yet'
         RETURN
      ENDIF
      IF (CCSD) THEN
         MODEL = 'CCSD'
      ENDIF
C
      IF (CIS) THEN
         MODELP = 'CIS'
      ELSE
         MODELP = MODEL
      ENDIF
C
      CALL AROUND( 'Calculation of '//MODELP// ' residues ')
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_QR2RSD Workspace:',LWORK
C
C-----------------------
C     Calculate property
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      NAQR2   = NQR2OP*NXQR2ST
      NBQR2   = NQR2OP*NXQR2ST
C
      KOSCS    = 1
      KOSCSF   = KOSCS  + NAQR2  
      KEND1    = KOSCSF + NBQR2  
      LEND1    = LWORK  - KEND1
      CALL DZERO(WORK(KOSCS),NAQR2)
      CALL DZERO(WORK(KOSCSF),NBQR2)
C
C----------------------------------------
C     Calculate linear response residues. 
C----------------------------------------
C
      DO 1000 IQR2 = 1, NXQR2ST
        ISTI = IQR2STI(IQR2)
        ISTF = IQR2STF(IQR2)
        ISYMI = ISYEXC(ISTI)   
        ISYMF = ISYEXC(ISTF) 
        ISTSYI = ISTI - ISYOFE(ISYMI)
        ISTSYF = ISTF - ISYOFE(ISYMF)
        EIGVI  = EIGVAL(ISTI)
        EIGVF  = EIGVAL(ISTF)
        IF (IPRINT .GT. 5) THEN
          WRITE(LUPRI,'(/,1x,A,2(/,1X,A,I3,1X,A,I3,A,F16.8))')
     *    'Calculating quadratic response 2. residues for: ',
     *     'State nr. ',ISTSYI,'of symmetry ',ISYMI,
     *     ' and eigenvalue: ',EIGVI,
     *     'State nr. ',ISTSYF,'of symmetry ',ISYMF,
     *     ' and eigenvalue: ',EIGVF 
        ENDIF
C
        DO 2000 IOPER = 1, NQR2OP
          ISYMA  = ISYOPR(IAQR2OP(IOPER))
          ISYMB  = ISYOPR(IBQR2OP(IOPER))
          ISYMAI = MULD2H(ISYMA,ISYMI)
          ISYMBF = MULD2H(ISYMB,ISYMF)
C
C------------------------------
C           Calculate residues.
C------------------------------
C
          LABELA = LBLOPR(IAQR2OP(IOPER))
          LABELB = LBLOPR(IBQR2OP(IOPER))
          IF ((ISYMF.EQ.ISYMAI).AND.(ISYMBF.EQ.ISYMI)) THEN
            KOSCSOF  = NQR2OP*(IQR2-1) + IOPER + KOSCS - 1
            CALL CC_QR2(LABELA,ISYMA,
     *                  ISTI,ISYMI,ISTF,ISYMF,   
     *                  WORK(KOSCSOF),WORK(KEND1),LEND1)
            KOSCSOF2 = NQR2OP*(IQR2-1) + IOPER + KOSCSF - 1
            CALL CC_QR2(LABELB,ISYMB,
     *                  ISTF,ISYMF,ISTI,ISYMI,   
     *                  WORK(KOSCSOF2),WORK(KEND1),LEND1)
          ENDIF
 2000   CONTINUE
 1000 CONTINUE
C
C---------------------------------------------------
C     Output quadratic response residue properties.
C     IF XOSCST put into oscillator strength tensor.
C---------------------------------------------------
C
      KOSCS2 = KEND1
      KTRS   = KOSCS2 + NEXCI*NEXCI*3*3
      KVELST = KTRS   + NEXCI*NEXCI*3*3
      KVELST2= KVELST + NEXCI*NEXCI*3*3
      KEND2  = KVELST2+ NEXCI*NEXCI*3*3
      LEND2  = LWORK  - KEND2
      CALL DZERO(WORK(KOSCS2),NEXCI*NEXCI*3*3)
      CALL DZERO(WORK(KTRS),NEXCI*NEXCI*3*3)
      CALL DZERO(WORK(KVELST),NEXCI*NEXCI*3*3)
      CALL DZERO(WORK(KVELST2),NEXCI*NEXCI*3*3)
C
      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     *  'Transition moments between excited states in atomic units(L):'
      WRITE(LUPRI,'(1X,A6,A,/)') '------',
     *  '-------------------------------------------------------------'
C
      DO IOPER = 1, NQR2OP
        ISYMA  = ISYOPR(IAQR2OP(IOPER))
        LABELA = LBLOPR(IAQR2OP(IOPER))
        ISYMB  = ISYOPR(IBQR2OP(IOPER))
        LABELB = LBLOPR(IBQR2OP(IOPER))
        DO IQR2  = 1, NXQR2ST
          ISTI     = IQR2STI(IQR2)
          ISTF     = IQR2STF(IQR2)
          ISYMI  = ISYEXC(ISTI)
          ISYMF  = ISYEXC(ISTF)
          ISYMAI = MULD2H(ISYMA,ISYMI)
          ISYMBF = MULD2H(ISYMB,ISYMF)
          ISTSYI = ISTI - ISYOFE(ISYMI)
          ISTSYF = ISTF - ISYOFE(ISYMF)
          EIGVI  = EIGVAL(ISTI)    
          EIGVF  = EIGVAL(ISTF)    
          IF ((ISYMF.EQ.ISYMAI).AND.(ISYMBF.EQ.ISYMI)) THEN
            K1     = NQR2OP*(IQR2-1) + IOPER + KOSCS - 1
            WRITE(LUPRI,'(1X,I2,1X,I2,1X,2(I2,F10.6),2X,A3,
     &           A8,A6,1X,F15.8)')
     *      IOPER,IQR2,ISTI,EIGVI,ISTF,EIGVF,
     *      '<i|',LABELA,'|f> = ',WORK(K1)    
          ENDIF
        END DO  
      END DO
C
      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     *  'Transition moments between excited states in atomic units(R):'
      WRITE(LUPRI,'(1X,A6,A,/)') '------',
     *  '-------------------------------------------------------------'
C
      DO IOPER = 1, NQR2OP
        ISYMA  = ISYOPR(IAQR2OP(IOPER))
        LABELA = LBLOPR(IAQR2OP(IOPER))
        ISYMB  = ISYOPR(IBQR2OP(IOPER))
        LABELB = LBLOPR(IBQR2OP(IOPER))
        DO IQR2  = 1, NXQR2ST
          ISTI     = IQR2STI(IQR2)
          ISTF     = IQR2STF(IQR2)
          ISYMI  = ISYEXC(ISTI)
          ISYMF  = ISYEXC(ISTF)
          ISYMBF = MULD2H(ISYMB,ISYMF)
          ISYMAI = MULD2H(ISYMA,ISYMI)
          ISTSYI = ISTI - ISYOFE(ISYMI)
          ISTSYF = ISTF - ISYOFE(ISYMF)
          EIGVI  = EIGVAL(ISTI)    
          EIGVF  = EIGVAL(ISTF)    
          IF ((ISYMF.EQ.ISYMAI).AND.(ISYMBF.EQ.ISYMI)) THEN
            K1     = NQR2OP*(IQR2-1) + IOPER + KOSCSF - 1
            WRITE(LUPRI,'(1X,I2,1X,I2,1X,2(I2,F10.6),2X,A3,A8,
     &           A6,1X,F15.8)')
     *      IOPER,IQR2,ISTF,EIGVF,ISTI,EIGVI,
     *      '<f|',LABELB,'|i> = ',WORK(K1)    
          ENDIF
        END DO  
      END DO
C
      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     *  'transition strength between excited states in atomic units:'
C
      WRITE(LUPRI,'(1X,A6,A,/)') '------',
     *  '-----------------------------------------------------------'
C
      DO IOPER = 1, NQR2OP
        ISYMA  = ISYOPR(IAQR2OP(IOPER))
        ISYMB  = ISYOPR(IBQR2OP(IOPER))
        LABELA = LBLOPR(IAQR2OP(IOPER))
        LABELB = LBLOPR(IBQR2OP(IOPER))
        DO IQR2  = 1, NXQR2ST
          ISTI     = IQR2STI(IQR2)
          ISTF     = IQR2STF(IQR2)
          ISYMI  = ISYEXC(ISTI)
          ISYMF  = ISYEXC(ISTF)
          ISYMBF = MULD2H(ISYMB,ISYMF)
          ISYMAI = MULD2H(ISYMA,ISYMI)
          ISTSYI = ISTI - ISYOFE(ISYMI)
          ISTSYF = ISTF - ISYOFE(ISYMF)
          EIGVI  = EIGVAL(ISTI)    
          EIGVF  = EIGVAL(ISTF)    
          IF ((ISYMF.EQ.ISYMAI).AND.(ISYMBF.EQ.ISYMI)) THEN
            K1     = NQR2OP*(IQR2-1) + IOPER + KOSCS - 1
            K2     = NQR2OP*(IQR2-1) + IOPER + KOSCSF - 1 
            RESIDUE = WORK(K1)*WORK(K2)
            IF (RESIDUE.GE.0.0D0) THEN
              SQRRES=SQRT(RESIDUE)
            ELSE 
              SQRRES=-SQRT(-RESIDUE)
            ENDIF
            WRITE(LUPRI,'(1X,A,A8,A1,A8,A2,F9.6,A,F9.6,A,F14.8,
     &           A,F12.8,A)')
     *      'S{',LABELA,',',LABELB,'}(',EIGVI,' -> ',EIGVF,') =',
     *      RESIDUE,' ( ',SQRRES,')'
            IF (XOSCST) THEN
              IADR1 = 0
              IADR2 = 0
              IF (LABELA(1:5).EQ.'XDIPL') IADR1 = 1
              IF (LABELA(1:5).EQ.'YDIPL') IADR1 = 2
              IF (LABELA(1:5).EQ.'ZDIPL') IADR1 = 3
              IF (LABELB(1:5).EQ.'XDIPL') IADR2 = 1
              IF (LABELB(1:5).EQ.'YDIPL') IADR2 = 2
              IF (LABELB(1:5).EQ.'ZDIPL') IADR2 = 3
              IF ((IADR1+IADR2).GE.2) THEN
                IEXAD = NEXCI*(ISTF - 1) + ISTI
                IOSCS2 = 3*3*(IEXAD-1)+3*(IADR2-1)+IADR1+KOSCS2-1
                WORK(IOSCS2) = RESIDUE     
              ENDIF
            ENDIF
            IF (XVELST) THEN
              IADR1 = 0
              IADR2 = 0
              IF (LABELA(1:5).EQ.'XDIPV') IADR1 = 1
              IF (LABELA(1:5).EQ.'YDIPV') IADR1 = 2
              IF (LABELA(1:5).EQ.'ZDIPV') IADR1 = 3
              IF (LABELB(1:5).EQ.'XDIPV') IADR2 = 1
              IF (LABELB(1:5).EQ.'YDIPV') IADR2 = 2
              IF (LABELB(1:5).EQ.'ZDIPV') IADR2 = 3
              IF ((IADR1+IADR2).GE.2) THEN
                IEXAD = NEXCI*(ISTF - 1) + ISTI
                IOSCS2 = 3*3*(IEXAD-1)+3*(IADR2-1)+IADR1+KVELST-1
                WORK(IOSCS2) = RESIDUE     
              ENDIF
            ENDIF
          ELSE
            SQRRES  = 0.0D0
            RESIDUE = 0.0D0
          ENDIF
          IF (LABELA.EQ.LABELB) THEN
             CALL WRIPRO(SQRRES,MODELP,-2,
     *                   LABELA,LABELB,LABELA,LABELB,
     *                   EIGVI,EIGVF,EIGVF-EIGVI,MULD2H(ISYMAI,ISYMF),
     *                   0,0,0)
             OSCCON = (EIGVF-EIGVI)*SQRRES*SQRRES 
             CALL WRIPRO(OSCCON,MODEL,-22,
     &                   LABELA,LABELB,LABELA,LABELB,
     &                   EIGVI,EIGVF,EIGVF-EIGVI,MULD2H(ISYMAI,ISYMF),
     &                   ISYME,1,ISTATE)
          ENDIF
        END DO  
c       WRITE(LUPRI,*) ' ' 
      END DO   

C
C-------------------------------------------
C     Perform analysis for dipole strenghts.
C-------------------------------------------
C
      CALL DCOPY(3*3*NEXCI*NEXCI,WORK(KOSCS2),1,WORK(KTRS),1)
      CALL DCOPY(3*3*NEXCI*NEXCI,WORK(KVELST),1,WORK(KVELST2),1)
C
C-----------------------------------------------
C     Write out strength for CCS, CC2, and CCSD.
C-----------------------------------------------
C
      LUOSC = LURES
      IF (XOSCST.AND. (CCS.OR.CC2.OR.CCSD)) THEN
C
         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '| Sym.I/F | I / F |       ',MODELP,' excited sta'
     *    //'te transition properties      |'
         WRITE(LUOSC,'(A)')
     *     ' |         |       +-----------------------------'
     *    //'------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '|         |       | Dipole Strength(a.u.) | Oscillator stre'
     *    //'ngth  | Direction  |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
         DO 9001 ISYMF = 1, NSYM  
          DO 9002 ISYMI = 1, ISYMF
           DO 9003 IEXF  = 1, NCCEXCI(ISYMF,1)
            DO 9004 IEXI  = 1, NCCEXCI(ISYMI,1)

              ISTI   = ISYOFE(ISYMI) + IEXI
              ISTF   = ISYOFE(ISYMF) + IEXF
              IEXAD  = NEXCI*(ISTF - 1) + ISTI
              IF ((.NOT.SELQR2).AND.(ISTI.GE.ISTF)) GO TO 9004

              LCALC  = .FALSE.
              LDIP   = .TRUE. 
              DO IQR2  = 1, NXQR2ST
                ISTII  = IQR2STI(IQR2)
                ISTIF  = IQR2STF(IQR2)
                IF ((ISTII.EQ.ISTI).AND.(ISTIF.EQ.ISTF)) LCALC = .TRUE.
              END DO
              KOSCSI = KOSCS2 + 3*3*(IEXAD-1)
              KTRSI  = KTRS   + 3*3*(IEXAD-1)
              CALL CC_XOSCPRI(WORK(KTRSI),WORK(KOSCSI),
     *                        EIGVAL(ISTI),EIGVAL(ISTF),
     *                        IEXI,ISYMI,IEXF,ISYMF,
     *                        WORK(KEND2),LEND2,MODELP,LCALC,
     *                        LDIP,LUOSC)
 9004       CONTINUE
 9003      CONTINUE
 
             IF (.NOT.((ISYMI.EQ.NSYM).OR.
     *           (NCCEXCI(ISYMI,1).EQ.0).OR.
     *           (NCCEXCI(ISYMF,1).EQ.0))) THEN
                NREST = 0
                DO 9013 ISYM2 = ISYMI+1,NSYM
                   NREST = NREST + NCCEXCI(ISYM2,1)
 9013           CONTINUE
                IF (NREST.EQ.0) GOTO 9002
                WRITE(LUOSC,'(A)')
     *          ' +----------------------------------------------'
     *         //'-------------------------------+'
             ENDIF
 9002     CONTINUE
 9001    CONTINUE
C
         WRITE(LUOSC,'(A)')
     *    ' +=============================================='
     *   //'===============================+'
C
      ENDIF
C
      LUOSC = LURES
      IF (XVELST.AND. (CCS.OR.CC2.OR.CCSD)) THEN
C
         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '| Sym.I/F | I / F |       ',MODELP,' excited sta'
     *    //'te transition properties      |'
         WRITE(LUOSC,'(A)')
     *     ' |         |       +-----------------------------'
     *    //'------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '|         |       | Veloc. Strength(a.u.) | Oscillator stre'
     *    //'ngth  | Direction  |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
         DO 9006 ISYMF = 1, NSYM  
          DO 9007 ISYMI = 1, ISYMF
           DO 9008 IEXF  = 1, NCCEXCI(ISYMF,1)
            DO 9009 IEXI  = 1, NCCEXCI(ISYMI,1)

              ISTI   = ISYOFE(ISYMI) + IEXI
              ISTF   = ISYOFE(ISYMF) + IEXF
              IEXAD  = NEXCI*(ISTF - 1) + ISTI
              IF ((.NOT.SELQR2).AND.(ISTI.GE.ISTF)) GO TO 9009

              LCALC  = .FALSE.
              LDIP   = .FALSE.
              DO IQR2  = 1, NXQR2ST
                ISTII  = IQR2STI(IQR2)
                ISTIF  = IQR2STF(IQR2)
                IF ((ISTII.EQ.ISTI).AND.(ISTIF.EQ.ISTF)) LCALC = .TRUE.
              END DO
              KOSCSI = KVELST + 3*3*(IEXAD-1)
              KTRSI  = KVELST2+ 3*3*(IEXAD-1)
              CALL CC_XOSCPRI(WORK(KTRSI),WORK(KOSCSI),
     *                        EIGVAL(ISTI),EIGVAL(ISTF),
     *                        IEXI,ISYMI,IEXF,ISYMF,
     *                        WORK(KEND2),LEND2,MODELP,LCALC,
     *                        LDIP,LUOSC)
 9009       CONTINUE
 9008      CONTINUE
 
             IF (.NOT.((ISYMI.EQ.NSYM).OR.
     *           (NCCEXCI(ISYMI,1).EQ.0).OR.
     *           (NCCEXCI(ISYMF,1).EQ.0))) THEN
                NREST = 0
                DO 9018 ISYM2 = ISYMI+1,NSYM
                   NREST = NREST + NCCEXCI(ISYM2,1)
 9018           CONTINUE
                IF (NREST.EQ.0) GOTO 9007
                WRITE(LUOSC,'(A)')
     *          ' +----------------------------------------------'
     *         //'-------------------------------+'
             ENDIF
 9007     CONTINUE
 9006    CONTINUE
C
         WRITE(LUOSC,'(A)')
     *    ' +=============================================='
     *   //'===============================+'
C
      ENDIF
C
      END
c*DECK CC_QR2
      SUBROUTINE CC_QR2(LABELA,ISYMA,
     *                  ISTI,ISYMI,ISTF,ISYMF,
     *                  RES1,WORK,LWORK)
C
C------------------------------------------------------------------------
C
C     Purpose: Calculate second residue of the quadratic response function.
C
C     Written by Ove Christiansen 26-4-1996
C
C------------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccexci.h>
#include <ccqr2r.h>
#include <dummy.h>
C
      PARAMETER( TWO = 2.0D00,TOLFRQ=1.0D-08 )

      DIMENSION WORK(LWORK)
      CHARACTER*8 LABELA,MODEL*10
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_QR2: Calculating residues   ')
      ENDIF
C
C------------------------
C     Allocate workspace.
C------------------------
C
      IF (ISYMA.NE.MULD2H(ISYMI,ISYMF)) 
     *                       CALL QUIT('Symmetry mismatch-2 in CC_QR2R')
C
      ISYMAI = MULD2H(ISYMA,ISYMI)
      NTAMPAI = NT1AM(ISYMAI) + NT2AM(ISYMAI)
      IF ( CCS ) NTAMPAI = NT1AM(ISYMAI)
C
      KETA  = 1
      KEND1 = KETA  + NTAMPAI
      LEND1 = LWORK - KEND1
      IF (LEND1 .LT. 0)
     *      CALL QUIT('Insufficient space for allocation in CC_QR2-1')
      KR1  = KEND1
      KR11 = KEND1
      KR12 = KEND1 + NT1AM(ISYMAI)
      KEND2 = KR1  + NTAMPAI
      LEND2 = LWORK - KEND2
      IF (LEND2 .LT. 0)
     *      CALL QUIT('Insufficient space for allocation in CC_QR2-2')
C
C--------------------------------------
C     Calculate Aa-matrix contribution.
C--------------------------------------
C
      CALL CC_ETAC(ISYMA,LABELA,WORK(KETA),'LE',ISTI,0,
     *             DUMMY,WORK(KEND1),LEND1)
C
      KR11 = KR1
      KR12 = KR1 + NT1AM(ISYMF)
      IOPT   = 3
      CALL CC_RDRSP('RE',ISTF,ISYMF,IOPT,MODEL,WORK(KR11),
     *              WORK(KR12))
C
      EATBCN = DDOT(NTAMPAI,WORK(KETA),1,WORK(KR1),1)
C
      IF ( IPRINT .GT. 9 ) THEN
          WRITE(LUPRI,*) ' Singles contribution:',
     *       DDOT(NT1AM(ISYMAI),WORK(KETA),1,WORK(KR1),1)
          IF (.NOT. CCS) WRITE(LUPRI,*) ' Doubles contribution:',
     *       DDOT(NT2AM(ISYMAI),WORK(KETA+NT1AM(ISYMAI)),1,
     *       WORK(KR1+NT1AM(ISYMAI)),1)
      ENDIF
C
C------------------------------------
C     Add to response function array.
C------------------------------------
C
      IF (IPRINT .GT. 2 ) THEN
          WRITE(LUPRI,'(1X,A3,A8,A3,A,F10.6)')
     *    '<i|',LABELA,'|f>',' LEi*A*REf cont. = ',EATBCN
      ENDIF
      RES1       = EATBCN  + RES1      
C
      NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
      IF ( CCS ) NTAMPA = NT1AM(ISYMA)
      KETA  = 1
      KEND1 = KETA  + NTAMPA
      LEND1 = LWORK - KEND1
      KR1  = KEND1
      KR11 = KEND1
      KR12 = KEND1 + NT1AM(ISYMA)
      KEND2 = KETA  + NTAMPA
      LEND2 = LWORK - KEND2
      IF (LEND2 .LT. 0)
     *      CALL QUIT('Insufficient space for allocation in CC_QR2-2')
C
C-------------------------------------
C     Calculate B-matrix contribution.
C-------------------------------------
C
      IF ((.NOT. CIS).AND.(.NOT.QR22N1)) THEN
         WRITE(LUPRI,*) 'Have not been programmed this stupid way. '
      ENDIF
      IF ((.NOT. CIS).AND.QR22N1) THEN
        CALL CC_XKSI(WORK(KETA),LABELA,ISYMA,0,DUMMY,WORK(KEND1),LEND1)
        ILSTNR = IN2AMP(ISTI,-EIGVAL(ISTI),ISYMI,
     *                  ISTF,EIGVAL(ISTF),ISYMF)
        CALL CC_RDRSP('N2',ILSTNR,ISYMA,IOPT,MODEL,WORK(KR11),
     *               WORK(KR12))
      ENDIF
C
      EATBCN = DDOT(NTAMPA,WORK(KETA),1,WORK(KR1),1)
C
      IF (IPRINT .GT. 9 .AND. (.NOT.CIS)) THEN
          WRITE(LUPRI,*) ' Singles contribution:',
     *       DDOT(NT1AM(ISYMA),WORK(KETA),1,WORK(KR1),1)
          IF (.NOT. CCS) WRITE(LUPRI,*) ' Doubles contribution:',
     *       DDOT(NT2AM(ISYMA),WORK(KETA+NT1AM(ISYMA)),1,
     *       WORK(KR1+NT1AM(ISYMA)),1)
      ENDIF
C
C------------------------------------
C     Add to response function array.
C------------------------------------
C
      IF (IPRINT .GT. 2 .AND.(.NOT.CIS)) THEN
        IF (.NOT.QR22N1) THEN
          WRITE(LUPRI,'(1X,A3,A8,A3,A,F10.6)')
     *    '<i|',LABELA,'|f>',' LEi*B*REf*t  = ',EATBCN
        ELSE
          WRITE(LUPRI,'(1X,A3,A8,A3,A,F10.6)')
     *    '<i|',LABELA,'|f>',' Nif*KsiAcont. = ',EATBCN
        ENDIF
      ENDIF
      IF (.NOT.CIS) RES1       = EATBCN  + RES1      
C
      RETURN
      END
      SUBROUTINE CC_XOSCPRI(TRS,OSC,EIGVI,EIGVF,IEXI,ISYMI,IEXF,ISYMF,
     *                      WORK,LWORK,MODEL,LCALC,LDIP,LUOSC)
C
C------------------------------------------------------------------------------
C
C              Write out transition/oscillator strength between excited states.
C
C     Written by Ove Christiansen April 1997 based on CC_OSCPRI
C
C------------------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
      PARAMETER (TOLFRQ = 1.0D-08,ONE= 1.0D0,THR = 1.0D-08)
C
#include <iratdef.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccqr2r.h>
C
      DIMENSION OSC(*),PVAL(3),PAXIS(3,3)
      CHARACTER MODEL*10,CDIP*7
      LOGICAL LCALC,LDIP
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_XOSCPRI: Output transition properties ' )
      ENDIF
C
C------------------------------------------
C     write out transition strength matrix.
C------------------------------------------
C
      IF (LCALC) THEN 
      EIGVFI = EIGVF - EIGVI 
      WRITE(LUPRI,'(//,1X,A6,A,I2,A,I3,A,I2,A,I3,A,/,A,F11.8,1X,
     *          F11.8,A,F11.8)') 
     *    MODEL(1:6),
     *    'Transition strength matrix for (',ISYMI,',',IEXI,') -> (', 
     *    ISYMF,',',IEXF,') transition',
     *    ' Excitation energies:',EIGVI,EIGVF,' and f-i energy ',EIGVFI
      IF (LDIP) THEN
         WRITE(LUPRI,'(1X,A)') 'Dipole gauge'
      ELSE
         WRITE(LUPRI,'(1X,A)') 'Velocity gauge'
      ENDIF
      CALL OUTPUT(TRS,1,3,1,3,3,3,1,LUPRI)
C
      CALL  TNSRAN(TRS,PVAL,PAXIS,
     *             ALFSQ,BETSQ,ITST,ITST2,
     *             APAR1,APEN1,XKAPPA,IPAR)
      WRITE(LUPRI,'(/,1X,A,/)')
     *    'Principal values of diagonalized transition strength matrix:'
      WRITE(LUPRI,'(1X,A)') '            a.u.               '
      WRITE(LUPRI,'(1X,A,F12.8)') '1     ',PVAL(1)
      WRITE(LUPRI,'(1X,A,F12.8)') '2     ',PVAL(2)
      WRITE(LUPRI,'(1X,A,F12.8)') '3     ',PVAL(3)
      WRITE(LUPRI,'(/,1X,A,/)')
     *    'Principal axis of diagonalized transition strength matrix:'
      CALL OUTPUT(PAXIS,1,3,1,3,3,3,1,LUPRI)
      TRA = PVAL(1)+PVAL(2)+PVAL(3)
C
      IF (IPAR .EQ.1) CDIP = '   X   '
      IF (IPAR .EQ.2) CDIP = '   Y   '
      IF (IPAR .EQ.3) CDIP = '   Z   '
      IF (IPAR .EQ.4) CDIP = ' (X,Y) '
      IF (IPAR .EQ.5) CDIP = ' (X,Z) '
      IF (IPAR .EQ.6) CDIP = ' (Y,Z) '
      IF (IPAR .EQ.7) CDIP = '(X,Y,Z)'
      IF (IPAR .EQ.8) CDIP = '   -   '
C
C------------------------------------------
C     First scale it - then
C     write out oscillator strength matrix.
C------------------------------------------
C
      IF (LDIP) THEN
         FACT = EIGVFI*2.0D0/3.0D0
      ELSE
         FACT = 2.0D0/(3.0D0*EIGVFI)
      ENDIF
      CALL DSCAL(3*3,FACT,OSC,1)
      WRITE(LUPRI,'(//,1X,A6,A,I2,A,I3,A,I2,A,I3,A,/,A,F11.8,1X,
     *          F11.8,A,F11.8/)') 
     *    MODEL(1:6),
     *    'Oscillator strength matrix for (',ISYMI,',',IEXI,') -> (', 
     *    ISYMF,',',IEXF,') transition',
     *    ' Excitation energies:',EIGVI,EIGVF,' and f-i energy ',EIGVFI
      CALL OUTPUT(OSC,1,3,1,3,3,3,1,LUPRI)
      CALL TNSRAN(OSC,PVAL,PAXIS,
     *            ALFSQ,BETSQ,ITST,ITST2,
     *            APAR2,APEN2,XKAPPA,IPAR)
      WRITE(LUPRI,'(/,1X,A,/)')
     *    'Principal values of diagonalized oscillator strength matrix:'
      WRITE(LUPRI,'(1X,A)') '            a.u.               '
      WRITE(LUPRI,'(1X,A,F12.8)') '1     ',PVAL(1)
      WRITE(LUPRI,'(1X,A,F12.8)') '2     ',PVAL(2)
      WRITE(LUPRI,'(1X,A,F12.8)') '3     ',PVAL(3)
      WRITE(LUPRI,'(/,1X,A,/)')
     *    'Principal axis of diagonalized oscillator strength matrix:'
      CALL OUTPUT(PAXIS,1,3,1,3,3,3,1,LUPRI)
      OSCS = PVAL(1)+PVAL(2)+PVAL(3)
c
      NR = 2
      IF ((.NOT.SELQR2).AND.(ISYMI.EQ.ISYMF)) NR = 3
c     IF ((IEXF+IEXI).EQ.NR) THEN
         WRITE(LUOSC,9988) ISYMI,ISYMF,IEXI,IEXF,TRA,OSCS,CDIP
c     ELSE
c        WRITE(LUOSC,9989) IEXI,IEXF,TRA,OSCS,CDIP
c     ENDIF
C
      ELSE IF (.NOT.LCALC) THEN
        CDIP = '   ?   '
c       IF ((IEXF+IEXI).EQ.NR) THEN
           WRITE(LUOSC,9986) ISYMI,ISYMF,IEXI,IEXF,'Not calculated',
     *                       'Not calculated',CDIP
c       ELSE
c          WRITE(LUOSC,9987) IEXI,IEXF,'Not calculated',
c    *                       'Not calculated',CDIP
c       ENDIF
      ENDIF
C
 9986 FORMAT(1X,'|',I3,I3,'   |',I3,I3,' | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,' | ')
 9987 FORMAT(1X,'|         |',I3,I3,' | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,' | ')
 9988 FORMAT(1X,'|',I3,I3,'   |',I3,I3,' | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,' | ')
 9989 FORMAT(1X,'|         |',I3,I3,' | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,' | ')
C
      END
