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
c*DECK CC_LR
       SUBROUTINE CC_LR(WORK,LWORK)
C
C----------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              polarizabilities.
C              (without orbital relaxation)
C
C              CIS, CCS, CC2, CCSD
C
C     Written by Ove Christiansen februar 1996.
C     Modified version for general linear response properties
C     Ove Christiansen November 1996.
C     New loop structure for general prop. Ove Christiansen April 1997.
C     SCF model added. Christof Haettig November 1998.
C     1/2 C^{+/-w} symmetrization introduced. Ch. Haettig, March 1999.
C
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (TOLFRQ=1.0D-08,ONE=1.0D0,XMONE=-1.0D0,THR=1.0D-08)
      INTEGER LUFCK
      PARAMETER (HALF = 0.5D0, ZERO = 0.0D0)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
C
#include <iratdef.h>
#include <inftap.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cclrinf.h>
#include <ccroper.h>
#include <ccr1rsp.h>
#include <ccrspprp.h>
#include <ccexpfck.h>
#include <ccfro.h>
#include <leinf.h>
#include <symmet.h>
#include <codata.h>
C
      LOGICAL FTSAV,LRLXA,LRLXB,LPDBSA,LPDBSB,LPRTSCF,OPTST,NOKAPPA
      LOGICAL SHIELD
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*10,MODELP*10
      CHARACTER LABELA*8, LABELB*8, LABSOP*8
      SAVE LPRTSCF
      DATA LPRTSCF /.TRUE./
      PARAMETER ( TWO = 2.0D0 )
C
C
C
C------------------------------------
C     Header of Property calculation.
C------------------------------------
C
      CALL QENTER('CC_LR')
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A)')
     *'*<<<<<<<<<< OUTPUT FROM COUPLED CLUSTER LINEAR RESPONSE >'//
     *'>>>>>>>>>*'
      IF ( DIPPOL ) THEN
         WRITE (LUPRI,'(1X,A)')
     *   '*                                                        '//
     *   '         *'
         WRITE (LUPRI,'(1X,A)')
     *   '*<<<<<<<<<<      CALCULATION OF CC POLARIZABILITIES     >'//
     *   '>>>>>>>>>*'
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
         MODEL = 'CC2       '
      ENDIF
      IF (CCS) THEN
         MODEL = 'CCS       '
      ENDIF
      IF (CC3  ) THEN
         MODEL = 'CC3       '
         WRITE(LUPRI,'(/,1x,A)') 'CC3 Polari not implemented yet'
         CALL QEXIT('CC_LR')
         RETURN
      ENDIF
      IF (CC1A) THEN
         MODEL = 'CCSDT-1a  '
         WRITE(LUPRI,'(/,1x,A)') 'CC1A Polari not implemented yet'
         CALL QEXIT('CC_LR')
         RETURN
      ENDIF
      IF (CC1B) THEN
         MODEL = 'CCSDT-1b  '
         WRITE(LUPRI,'(/,1x,A)') 'CC1B Polari not implemented yet'
         CALL QEXIT('CC_LR')
         RETURN
      ENDIF
      IF (CCSD) THEN
         MODEL = 'CCSD      '
      ENDIF
C
      IF (CIS) THEN
         MODELP = 'CIS       '
      ELSE
         MODELP = MODEL
      ENDIF
C
      CALL AROUND( 'Calculation of '//MODELP//
     *             ' linear response properties ')
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_LR Workspace:',LWORK
C
      CALL FLSHFO(LUPRI)
C
      NLRPRP  = NLROP*NBLRFR
C
C     --------------------------------------------------------------
C     open AOPROPER file for GETGPV routine of the RSP program...
C     --------------------------------------------------------------
C
      CALL CC_SIRINF(NCMOT,NASHT,N2ASHX,LCINDX)
C
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
C
C     -------------------------------
C     allocate workspace for results:
C     -------------------------------
C
      KCMO    = 1
      KUDV    = KCMO    + NCMOT
      KXINDX  = KUDV    + N2ASHX
      KR2EFF  = KXINDX  + LCINDX
      KFOCK0  = KR2EFF  + N2BST(1)
      KOVERLP = KFOCK0  + N2BST(1)
      KEND1   = KOVERLP + N2BST(1)

      KPOL    = KEND1
      KPOLF   = KPOL    + 2*NLRPRP
      KPOLSCF = KPOLF   + 2*NLRPRP
      KEND1   = KPOLSCF + 2*NLRPRP

      LEND1   = LWORK    - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_LR.')
      END IF

      CALL DZERO(WORK(KPOL),2*NLRPRP)
      CALL DZERO(WORK(KPOLF),2*NLRPRP)
      CALL DZERO(WORK(KPOLSCF),2*NLRPRP)
C
C     ------------------------------
C     read MO coefficient from file:
C     ------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NCMOT,WORK(KCMO))
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C     ------------------------------------
C     loop over operators and frequencies:
C     ------------------------------------
C
      NSCF = 0
C
      DO 1000 IOPER = 1, NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        LRLXA  = LALORX(IOPER)
        LRLXB  = LBLORX(IOPER)
        ISAMA  = ISYMAT(IOPERA)
        ISAMB  = ISYMAT(IOPERB)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        LPDBSA = LPDBSOP(IOPERA)
        LPDBSB = LPDBSOP(IOPERB)

        ISAPROP = ISAMA * ISAMB

        IF (ISYMA.EQ.ISYMB) THEN
          DO IFREQ = 1, NBLRFR
          DO ISIGN = +1, -1, -2
C
            IOFSGN = ((-ISIGN+1)/2) * NLRPRP
C
            SIGN   = DBLE(ISIGN)
            FREQA  = SIGN * ALRFR(IFREQ)
            FREQB  = SIGN * BLRFR(IFREQ)
C
            IF (IPRINT .GT. 5 .OR. LOCDBG) THEN
              WRITE(LUPRI,'(/,1x,A,F16.8,/,A,I2,/,3A,/,A,2L3,/,A,2L3)')
     *         'Calculating response property with frequency',FREQB,
     *         ' Operator symmetry = ',ISYMB,
     *         ' Labels = ',LABELA, LABELB,
     *         ' orbital relaxation flags = ',LRLXA, LRLXB,
     *         ' pert.-dep. basis set flags = ',LPDBSA, LPDBSB
            ENDIF
C
            KPRP1 = KPOL  + IOFSGN + NBLRFR*(IOPER-1) + IFREQ - 1
            KPRP2 = KPOLF + IOFSGN + NBLRFR*(IOPER-1) + IFREQ - 1
C
C-------------------------------------------
C           The etaA*tB(omeg) contributions.
C-------------------------------------------
C
            CALL CC_EATB(LABELA,ISYMA,FREQA,LRLXA,LPDBSA,
     *                   LABELB,ISYMB,FREQB,LRLXB,LPDBSB,
     *                   WORK(KPRP1),WORK(KEND1),LEND1)
C
            IF ( .NOT. ASYMSD) THEN
C
C-------------------------------------------------------
C             IF ordinatry form the calculate EATB form.
C-------------------------------------------------------
C
              CALL CC_EATB(LABELB,ISYMB,FREQB,LRLXB,LPDBSB,
     *                     LABELA,ISYMA,FREQA,LRLXA,LPDBSA,
     *                     WORK(KPRP1),WORK(KEND1),LEND1)
C
C--------------------------------------------------
C             The FtA(-omeg)*tB(omeg) contribution.
C--------------------------------------------------
C
              IF (.NOT.CIS) THEN
               CALL CC_FABCON(LABELA,ISYMA,FREQA,LRLXA,
     *                        LABELB,ISYMB,FREQB,LRLXB,
     *                        WORK(KPRP2),WORK(KEND1),LEND1)
              ENDIF
C
C---------------------------------------------------
C             The Pt-barA(-omeg)*t-barB(omeg)
C             contribution for solvent calculations.
C---------------------------------------------------
C
              IF (CCSLV) THEN
                CALL CC_PABCON(LABELA,ISYMA,FREQA,LRLXA,
     *                         LABELB,ISYMB,FREQB,LRLXB,
     *                         WORK(KPRP2),WORK(KEND1),LEND1)
              ENDIF
C
            ELSE
C
C-------------------------------------------------------------------
C           Use asymmetric form for calculating polarizability.
C           Calculate trivial LAKSIB contribution to polarizability.
C-------------------------------------------------------------------
C
               CALL CC_LAKSIB(LABELA,ISYMA,FREQA,LRLXA,
     *                        LABELB,ISYMB,FREQB,LRLXB,
     *                        WORK(KPRP2),WORK(KEND1),LEND1)
C
            ENDIF

C-------------------------------------------------------------------
C             construct the zeroth-order eff. CC Fock matrix in MO
C             and initialize the 'eff.' sec. order connection matrix
C-------------------------------------------------------------------
              IF ( (LRLXA .OR. LPDBSA) .AND. (LRLXB .OR. LPDBSB) ) THEN
                IFOCK = IEFFFOCK('HAM0    ',ISYM0,1)
                IADRF = IADRFCK(1,IFOCK)

                LUFCK = -1
                CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
                CALL GETWA2(LUFCK,FILFCKEFF,WORK(KFOCK0),
     &                      IADRF,N2BST(ISYM0))
                CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')

                CALL RDONEL('OVERLAP ',.TRUE.,WORK(KEND1),NBAST)
                CALL CCSD_SYMSQ(WORK(KEND1),ISYM0,WORK(KOVERLP))

                CALL CC_EFFCKMO(WORK(KFOCK0),ISYM0,WORK(KCMO),
     &                          WORK(KOVERLP),WORK(KEND1),LEND1)

              ELSE
                CALL DZERO(WORK(KFOCK0),N2BST(1))
              END IF

              CALL DZERO(WORK(KR2EFF),N2BST(1))

C-------------------------------------------------------------------
C             construct the X^(1) interm. for the A perturbation and
C             calculate its contribution to the response function:
C-------------------------------------------------------------------
              RLXBCON = ZERO

              IF (LRLXB.OR.LPDBSB) THEN

                 KXIMA  = KEND1
                 KAPB   = KXIMA  + N2BST(ISYMA)
                 KQMATH = KAPB   + 2*NALLAI(ISYMB)
                 KQMATP = KQMATH + MAX(N2BST(ISYMB),N2BST(ISYMA))
                 KRMAT  = KQMATP + MAX(N2BST(ISYMB),N2BST(ISYMA))
                 KAPBSQ = KRMAT  + MAX(N2BST(ISYMB),N2BST(ISYMA))
                 KQTRP  = KAPBSQ + N2BST(ISYMB)
                 KEND2  = KQTRP  + MAX(N2BST(ISYMB),N2BST(ISYMA))
                 LWRK2  = LWORK  - KEND2
                 IF (LWRK2 .LT. 0) THEN
                   CALL QUIT('Insufficient memory in CC_LR.')
                 END IF


                 CALL CCRLXXIM(WORK(KXIMA),ISYMA,LABELA,LRLXA,LPDBSA,
     &                         FREQA,WORK(KCMO),WORK(KEND2),LWRK2)

                 IF (LRLXB) THEN
                    IKAPPA = IR1KAPPA(LABELB,FREQB,ISYMB)
                    CALL CC_RDHFRSP('R1 ',IKAPPA,ISYMB,WORK(KAPB))
                 ELSE
                    CALL DZERO(WORK(KAPB),2*NALLAI(ISYMB))
                 END IF

                 CALL CC_GET_RMAT(WORK(KRMAT),IOPERB,1,ISYMB,
     &                            WORK(KEND2),LWRK2)
                 NOKAPPA = .FALSE.
                 CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),
     &                        WORK(KRMAT),WORK(KAPB),
     &                        ISAMB,ISYMB,NOKAPPA,WORK(KCMO),
     &                        WORK(KEND2),LWRK2)

                 DO ISYM1 = 1, NSYM
                    ISYM2 = MULD2H(ISYM1,ISYMB)
                    KOFF1 = KQMATH + IAODIS(ISYM1,ISYM2)
                    KOFF2 = KQTRP  + IAODIS(ISYM2,ISYM1)
                    CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                          WORK(KOFF1),WORK(KOFF2))
                 END DO
                 CALL DCOPY(N2BST(ISYMB),WORK(KQTRP),1,WORK(KQMATH),1)
                 CALL DSCAL(N2BST(ISYMB),-HALF,WORK(KQMATH),1)

                 DO ISYM1 = 1, NSYM
                    ISYM2 = MULD2H(ISYM1,ISYMB)
                    KOFF1 = KQMATP + IAODIS(ISYM1,ISYM2)
                    KOFF2 = KQTRP  + IAODIS(ISYM2,ISYM1)
                    CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                          WORK(KOFF1),WORK(KOFF2))
                 END DO
                 CALL DCOPY(N2BST(ISYMB),WORK(KQTRP),1,WORK(KQMATP),1)
                 CALL DSCAL(N2BST(ISYMB),-HALF,WORK(KQMATP),1)

                 RLXBCON =
     &                - DDOT(N2BST(ISYMA),WORK(KQMATH),1,WORK(KXIMA),1)
     &                - DBLE(ISAMA) *
     &                  DDOT(N2BST(ISYMA),WORK(KQMATP),1,WORK(KXIMA),1)

                 IF (LOCDBG) THEN
                   WRITE(LUPRI,*) 'XIMA for RLXBCON:'
                   CALL CC_PRONELAO(WORK(KXIMA),ISYMA)
                   WRITE(LUPRI,*) 'transpose QMATH:'
                   CALL CC_PRONELAO(WORK(kqtrp),ISYMB)
                   WRITE(LUPRI,*) 'RLXBCON:',RLXBCON
                 END IF

                 WORK(KPRP1) = WORK(KPRP1) + RLXBCON

                 CALL CCKAPPASQ(WORK(KAPBSQ),WORK(KAPB),ISYMB,'N')

                 CALL CC_GET_RMAT(WORK(KRMAT),IOPERA,1,ISYMA,
     &                            WORK(KEND2),LWRK2)

                 NOKAPPA = .TRUE.
                 CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),
     &                        WORK(KRMAT),DUMMY,
     &                        ISAMA,ISYMA,NOKAPPA,WORK(KCMO),
     &                        WORK(KEND2),LWRK2)

                 CALL CC_MMOMMO('N','N',+1.0D0,WORK(KAPBSQ),ISYMB,
     &                          WORK(KQMATH),ISYMA,1.0D0,WORK(KR2EFF),1)
                 CALL CC_MMOMMO('N','N',-1.0D0,WORK(KQMATH),ISYMA,
     &                          WORK(KAPBSQ),ISYMB,1.0D0,WORK(KR2EFF),1)

                 IF (LOCDBG .OR. IPRINT.GT.1) THEN
                    WRITE (LUPRI,*) 'CC_LR> RLXBCON = ',RLXBCON
                    WRITE (LUPRI,*) 'CC_LR> PRP1    = ',WORK(KPRP1)
                 END IF
                 IF (LOCDBG) THEN
                    WRITE (LUPRI,*) 'RMAT A: AO'
                    CALL CC_PRONELAO(WORK(KRMAT),ISYMA)
                    WRITE (LUPRI,*) 'RMAT A: MO'
                    CALL CC_PRONELAO(WORK(KQMATH),ISYMA)
                    WRITE (LUPRI,*) 'KAPPA B:'
                    CALL CC_PRONELAO(WORK(KAPBSQ),ISYMB)
                    WRITE (LUPRI,*) 'KR2EFF:'
                    CALL CC_PRONELAO(WORK(KR2EFF),ISYM0)
                 END IF
              END IF

C-------------------------------------------------------------------
C             construct the X^(1) interm. for the B perturbation and
C             calculate its contribution to the response function:
C-------------------------------------------------------------------
              RLXACON = ZERO

              IF (LRLXA .OR. LPDBSA) THEN

                 KXIMB  = KEND1
                 KAPA   = KXIMB  + N2BST(ISYMB)
                 KAPASQ = KAPA   + 2*NALLAI(ISYMA)
                 KRMAT  = KAPASQ + N2BST(ISYMA)
                 KQMATH = KRMAT  + MAX(N2BST(ISYMA),N2BST(ISYMB))
                 KQMATP = KQMATH + MAX(N2BST(ISYMA),N2BST(ISYMB))
                 KQTRP  = KQMATP + MAX(N2BST(ISYMA),N2BST(ISYMB))
                 KEND2  = KQTRP  + MAX(N2BST(ISYMA),N2BST(ISYMB))
                 LWRK2  = LWORK  - KEND2
                 IF (LWRK2 .LT. 0) THEN
                   CALL QUIT('Insufficient memory in CC_LR.')
                 END IF


                 CALL CCRLXXIM(WORK(KXIMB),ISYMB,LABELB,LRLXB,LPDBSA,
     &                         FREQB,WORK(KCMO),WORK(KEND2),LWRK2)

                 IF (LRLXA) THEN
                    IKAPPA = IR1KAPPA(LABELA,FREQA,ISYMA)
                    CALL CC_RDHFRSP('R1 ',IKAPPA,ISYMA,WORK(KAPA))
                 ELSE
                    CALL DZERO(WORK(KAPA),2*NALLAI(ISYMA))
                 END IF

                 CALL CC_GET_RMAT(WORK(KRMAT),IOPERA,1,ISYMA,
     &                            WORK(KEND2),LWRK2)

                 NOKAPPA = .FALSE.
                 CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),
     &                        WORK(KRMAT),WORK(KAPA),
     &                        ISAMA,ISYMA,NOKAPPA,WORK(KCMO),
     &                        WORK(KEND2),LWRK2)

                 DO ISYM1 = 1, NSYM
                    ISYM2 = MULD2H(ISYM1,ISYMB)
                    KOFF1 = KQMATH + IAODIS(ISYM1,ISYM2)
                    KOFF2 = KQTRP  + IAODIS(ISYM2,ISYM1)
                    CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                          WORK(KOFF1),WORK(KOFF2))
                 END DO
                 CALL DCOPY(N2BST(ISYMB),WORK(KQTRP),1,WORK(KQMATH),1)
                 CALL DSCAL(N2BST(ISYMB),-HALF,WORK(KQMATH),1)

                 DO ISYM1 = 1, NSYM
                    ISYM2 = MULD2H(ISYM1,ISYMB)
                    KOFF1 = KQMATP + IAODIS(ISYM1,ISYM2)
                    KOFF2 = KQTRP  + IAODIS(ISYM2,ISYM1)
                    CALL TRSREC(NBAS(ISYM1),NBAS(ISYM2),
     &                          WORK(KOFF1),WORK(KOFF2))
                 END DO
                 CALL DCOPY(N2BST(ISYMB),WORK(KQTRP),1,WORK(KQMATP),1)
                 CALL DSCAL(N2BST(ISYMB),-HALF,WORK(KQMATP),1)

                 RLXACON =
     &                - DDOT(N2BST(ISYMB),WORK(KQMATH),1,WORK(KXIMB),1)
     &                - DBLE(ISAMB) *
     &                  DDOT(N2BST(ISYMB),WORK(KQMATP),1,WORK(KXIMB),1)

                 if (locdbg) then
                   WRITE(LUPRI,*) 'XIMB for RLXACON:'
                   call cc_pronelao(work(kximb),isymb)
                   WRITE(LUPRI,*) 'transpose QMATH:'
                   call cc_pronelao(work(KQTRP),isymb)
                   WRITE(LUPRI,*) 'RLXACON:',RLXACON
                 end if

                 WORK(KPRP1) = WORK(KPRP1) + RLXACON


                 CALL CCKAPPASQ(WORK(KAPASQ),WORK(KAPA),ISYMA,'N')

                 CALL CC_GET_RMAT(WORK(KRMAT),IOPERB,1,ISYMB,
     &                            WORK(KEND2),LWRK2)

                 NOKAPPA = .TRUE.
                 CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),
     &                        WORK(KRMAT),DUMMY,
     &                        ISAMB,ISYMB,NOKAPPA,WORK(KCMO),
     &                        WORK(KEND2),LWRK2)


                 CALL CC_MMOMMO('N','N',+1.0D0,WORK(KAPASQ),ISYMA,
     &                          WORK(KQMATH),ISYMB,1.0D0,WORK(KR2EFF),1)
                 CALL CC_MMOMMO('N','N',-1.0D0,WORK(KQMATH),ISYMB,
     &                          WORK(KAPASQ),ISYMA,1.0D0,WORK(KR2EFF),1)

                 IF (LOCDBG .OR. IPRINT.GT.1) THEN
                    WRITE (LUPRI,*) 'CC_LR> RLXACON = ',RLXACON
                    WRITE (LUPRI,*) 'CC_LR> PRP1    = ',WORK(KPRP1)
                 END IF
                 IF (LOCDBG) THEN
                    WRITE (LUPRI,*) 'RMAT B:'
                    CALL CC_PRONELAO(WORK(KQMATH),ISYMB)
                    WRITE (LUPRI,*) 'KAPPA A:'
                    CALL CC_PRONELAO(WORK(KAPASQ),ISYMA)
                    WRITE (LUPRI,*) 'KR2EFF:'
                    CALL CC_PRONELAO(WORK(KR2EFF),ISYM0)
                 END IF
              END IF

              IF (LPDBSA .OR. LPDBSB) THEN
                CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,ISGNSOP,
     *                             INUM,WORK(KEND1),LEND1)
                IF (INUM.LT.0) CALL QUIT('Operator error in CC_LR.')
                IEXPV = IEXPECT(LABSOP,ISYSOP)
                XSTAT = DBLE(ISGNSOP) *
     *                    ( EXPVALUE(1,IEXPV) + EXPVALUE(2,IEXPV) )
                XNUCL = CC_NUCCON(LABSOP,ISYSOP)

                XREO = TWO*DDOT(N2BST(1),WORK(KR2EFF),1,WORK(KFOCK0),1)

                IF (LOCDBG .OR. IPRINT.GT.1) THEN
                   WRITE (LUPRI,*) LABSOP,
     *                    EXPVALUE(1,IEXPV),EXPVALUE(2,IEXPV)
                   WRITE (LUPRI,*)
     *                    'CC_LR>  contrib. of Fock^(eff,0) :',XREO
                END IF
                IF (LOCDBG) THEN
                   WRITE (LUPRI,*)
     *                    'CC_LR>  [K^(A),R^(B)]+[K^(B),R^(A)] :'
                   CALL CC_PRONELAO(WORK(KR2EFF),1)
                   WRITE (LUPRI,*) 'CC_LR>  Fock^(eff,0) :'
                   CALL CC_PRONELAO(WORK(KFOCK0),1)
                END IF
              ELSE
                XSTAT = ZERO
                XNUCL = ZERO
                XREO  = ZERO
              END IF

              WORK(KPRP1) = WORK(KPRP1) + XREO + XSTAT - XNUCL

              IF (LOCDBG .OR. IPRINT.GT.10) THEN
                 WRITE (LUPRI,*) 'CC_LR> RLXACON   = ',RLXACON
                 WRITE (LUPRI,*) 'CC_LR> RLXBCON   = ',RLXBCON
                 WRITE (LUPRI,*) 'CC_LR> XSTAT(CC) = ',XSTAT
                 WRITE (LUPRI,*) 'CC_LR> XNUCL     = ',XNUCL
                 WRITE (LUPRI,*) 'CC_LR> XREO      = ',XREO
                 WRITE (LUPRI,*) 'CC_LR> PRP1      = ',WORK(KPRP1)
                 WRITE (LUPRI,*) 'CC_LR> PRP2      = ',WORK(KPRP2)
              END IF
C
C--------------------------------------------------------------
C             in relaxed case calculate SCF result if possible:
C--------------------------------------------------------------
C
              IF (LRLXA.AND.LRLXB) THEN

               IF (LEND1 .LT. 4*NALLAI(ISYMA)) THEN
                 CALL QUIT('Insufficient memory in CC_LR.')
               END IF

               KG1    = KEND1
               LWRKG1 = LWORK - KG1

               KG2    = KG1    + NALLAI(ISYMA)
               KAPPA1 = KG2    + NALLAI(ISYMA)
               KAPPA2 = KAPPA1 + NALLAI(ISYMA)

               NSCF = NSCF + 1
               KPRP = KPOLSCF + IOFSGN + NBLRFR*(IOPER-1) + IFREQ - 1

               IDXR = IR1KAPPA(LABELA,+FREQA,ISYMA)
               CALL CC_GETHFGD(IDXR,'R1 ',LRTHFLBL,IDUM,IDUM,RDUM,
     *                         ISYLRTHF,FRQLRTHF,IDUM,NLRTHFLBL,
     *                         MAXTLBL,IREAL,WORK(KCMO),WORK(KUDV),
     *                         WORK(KXINDX),FRVAL,WORK(KG1),LWRKG1)

               IDXR = IR1KAPPA(LABELB,+FREQB,ISYMB)
               CALL CC_RDHFRSP('R1 ',IDXR,ISYMB,WORK(KAPPA1))

               XRLXAB=DDOT(2*NALLAI(ISYMB),WORK(KAPPA1),1,WORK(KG1),1)


               IDXR = IR1KAPPA(LABELB,+FREQB,ISYMB)
               CALL CC_GETHFGD(IDXR,'R1 ',LRTHFLBL,IDUM,IDUM,RDUM,
     *                         ISYLRTHF,FRQLRTHF,IDUM,NLRTHFLBL,
     *                         MAXTLBL,IREAL,WORK(KCMO),WORK(KUDV),
     *                         WORK(KXINDX),FRVAL,WORK(KG1),LWRKG1)

               IDXR = IR1KAPPA(LABELA,+FREQA,ISYMA)
               CALL CC_RDHFRSP('R1 ',IDXR,ISYMA,WORK(KAPPA1))

               XRLXBA=DDOT(2*NALLAI(ISYMA),WORK(KAPPA1),1,WORK(KG1),1)

               WORK(KPRP) = XRLXAB

               ERROR = XRLXBA - DBLE(ISAPROP) * XRLXAB

               IF (LOCDBG.OR.DABS(ERROR).GT.THRLEQ.OR.IPRINT.GT.1) THEN
                  WRITE (LUPRI,*)'CC_LR>', LABELA,FREQA,LABELB,FREQB
                  WRITE (LUPRI,*)'CC_LR> ',XRLXAB,XRLXBA,ERROR,THRLEQ
                  IF (ERROR.GT.THRLEQ) THEN
                     WRITE (LUPRI,*)
     *                     'Warning: large errors in SCF second-',
     *                       'order property encountered!!!'
                  END IF
               END IF

               KFOCK1 = KEND1
               KR1DEN = KFOCK1 + N2BST(ISYMA)
               KEND2  = KR1DEN + N2BST(ISYMB)
               LWRK2  = LWORK  - KEND2

               LUFCK = -1
               IFOCK = IEFFFOCK(LABELA,ISYM,1)
               IADRF = IADRFCK(2,IFOCK)
               CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
               CALL GETWA2(LUFCK,FILFCKEFF,WORK(KFOCK1),
     &                     IADRF,N2BST(ISYMA))
               CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')

               CALL CC_HFR1DEN(WORK(KR1DEN),IOPERB,1,ISYMB,
     &                         WORK(KEND2),LWRK2)

               XREOB = -TWO * DDOT(N2BST(ISYMA),WORK(KFOCK1),1,
     &                                          WORK(KR1DEN),1)
               IF (LOCDBG .OR. IPRINT.GT.1) THEN
                  WRITE (LUPRI,*) 'CC_LR> XREOB = ',XREOB
               END IF

               KFOCK1 = KEND1
               KR1DEN = KFOCK1 + N2BST(ISYMB)
               KEND2  = KR1DEN + N2BST(ISYMA)
               LWRK2  = LWORK  - KEND2

               LUFCK = -1
               IFOCK = IEFFFOCK(LABELB,ISYM,1)
               IADRF = IADRFCK(2,IFOCK)
               CALL WOPEN2(LUFCK,FILFCKEFF,64,0)
               CALL GETWA2(LUFCK,FILFCKEFF,WORK(KFOCK1),
     &                     IADRF,N2BST(ISYMB))
               CALL WCLOSE2(LUFCK,FILFCKEFF,'KEEP')

               CALL CC_HFR1DEN(WORK(KR1DEN),IOPERA,1,ISYMA,
     &                         WORK(KEND2),LWRK2)

               XREOA = - TWO * DDOT(N2BST(ISYMA),WORK(KFOCK1),1,
     &                                           WORK(KR1DEN),1)
               IF (LOCDBG .OR. IPRINT.GT.1) THEN
                 WRITE (LUPRI,*) 'CC_LR> XREOA = ',XREOA
               END IF

               IF (LPDBSA .OR. LPDBSB) THEN
                 CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYSOP,
     *                              ISGNSOP,INUM,WORK(KEND1),LEND1)
                 IF (INUM.LT.0) CALL QUIT('Operator error in CC_LR.')
                 IEXPV = IEXPECT(LABSOP,ISYSOP)
                 XSTAT = EXPVALUE(3,IEXPV) + EXPVALUE(4,IEXPV)
                 XNUCL = CC_NUCCON(LABSOP,ISYSOP)
               ELSE
                 XSTAT = ZERO
                 XNUCL = ZERO
               END IF

               WORK(KPRP) = WORK(KPRP) + XREOA+XREOB+XNUCL+XSTAT

               IF (LOCDBG .OR. IPRINT.GT.1) THEN
                 WRITE (LUPRI,*) 'SCF <<',LABELA,';',LABELB,'>> : '
                 WRITE (LUPRI,*) 'relaxation contribution:',XRLXAB
                 WRITE (LUPRI,*) 'reorthog.  contribution:',XREOA+XREOB
                 WRITE (LUPRI,*) 'static electronic cont.:',XSTAT
                 WRITE (LUPRI,*) 'nuclear    contribution:',XNUCL
                 WRITE (LUPRI,*) 'total result           :',WORK(KPRP)
               END IF

              END IF
C
          END DO
          END DO
C
        ENDIF
 1000 CONTINUE
C
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
C
C------------------------------------
C     Output SCF response properties:
C------------------------------------
C
      IF ( NSCF.GT.1 .AND. (LPRTSCF.OR.LOCDBG) ) THEN
C
        WRITE(LUPRI,'(//,1X,A)')
     *    'SCF linear response properties in atomic units:'
        WRITE(LUPRI,'(1X,A,/)')
     *    '-----------------------------------------------'
C
        DO IOPER  = 1,NLROP
          IOPERA = IALROP(IOPER)
          IOPERB = IBLROP(IOPER)
          LRLXA  = LALORX(IOPER)
          LRLXB  = LBLORX(IOPER)
          ISYMA  = ISYOPR(IOPERA)
          ISYMB  = ISYOPR(IOPERB)
          LABELA = LBLOPR(IOPERA)
          LABELB = LBLOPR(IOPERB)
          LPDBSA = LPDBSOP(IOPERA)
          LPDBSB = LPDBSOP(IOPERB)
          IF(LRLXA.AND.LRLXB)THEN
           DO IFREQ = 1, NBLRFR
             KPRP1 = KPOLSCF + NBLRFR*(IOPER-1) + IFREQ - 1
             IF (ISYMA.EQ.ISYMB) THEN
               WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8)') '<<',
     *           LABELA,',',LABELB,'>>(',BLRFR(IFREQ),') =',WORK(KPRP1)
              ELSE
               WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8,A)') '<<',
     *           LABELA,',',LABELB,'>>(',BLRFR(IFREQ),') =',WORK(KPRP1),
     *           ' BY SYMMETRY !'
              ENDIF
           END DO
          END IF
        END DO
C
        LPRTSCF = .FALSE.
C
      END IF
C
C-------------------------------------------------
C        Output Linear response properties.
C        IF DIPPOL put into polarizability tensor.
C-------------------------------------------------
C
      KPOL2 = KEND1
      KEND2 = KPOL2 + NBLRFR*3*3
      LEND2 = LWORK - KEND2
C
      CALL DZERO(WORK(KPOL2),3*3*NBLRFR)
C
      CALL DAXPY(2*NLRPRP,ONE,WORK(KPOLF),1,WORK(KPOL),1)
C
      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     *  'linear response properties in atomic units:'
      WRITE(LUPRI,'(1X,A,/)')
     *  '-------------------------------------------------'
C
      DO 4000 IOPER  = 1,NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        ISYMAB = MULD2H(ISYMA,ISYMB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        ISAMA  = ISYMAT(IOPERA)
        ISAMB  = ISYMAT(IOPERB)

        ISAPROP = ISAMA * ISAMB
        SIGN    = DBLE(ISAPROP)

        IF ((LABELA(1:5).EQ.'dh/dB'.AND.LABELB(1:4).EQ.'PSO ').OR.
     *      (LABELB(1:5).EQ.'dh/dB'.AND.LABELA(1:4).EQ.'PSO ')     )THEN
         SHIELD = .TRUE.
         FACTOR = 1.0D06 * ALPHA2 ! conversion to ppm
        ELSE
         SHIELD = .FALSE.
         FACTOR = 1.0D0
        END IF

        DO IFREQ = 1, NBLRFR
            KPRP1P = KPOL +          NBLRFR*(IOPER-1) + IFREQ - 1
            KPRP1M = KPOL + NLRPRP + NBLRFR*(IOPER-1) + IFREQ - 1

            RESULT = HALF*( WORK(KPRP1P) + SIGN * WORK(KPRP1M) )
            ERROR  = HALF*( WORK(KPRP1P) - SIGN * WORK(KPRP1M) )

            IF (IPRINT.GT.11 .OR. ISAPROP.EQ.0) THEN

              IF (ISAPROP .EQ. 0) THEN
                WRITE(LUPRI,'(/1X,A,/1X,A)')
     *           'Cannot determine if real or imaginary property...',
     *           'the non-symmetrized results for +/- w are:'
              ELSE
                 WRITE(LUPRI,'(/1X,A)') 'non-symmetrized '//
     &                'results for +/-w:'
              ENDIF

              WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8)')
     *        '<<',LABELA,',',LABELB,
     *        '>>(',BLRFR(IFREQ),') =',WORK(KPRP1P)
              WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8)')
     *        '<<',LABELA,',',LABELB,
     *        '>>(',-BLRFR(IFREQ),') =',WORK(KPRP1M)

              WRITE(LUPRI,'(1X,A)')
     &             'symmetric/antisymmetric contributions:'
              WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8,A,F15.8)')
     *        '<<',LABELA,',',LABELB, '>>(',-BLRFR(IFREQ),') =',
     *            HALF*(WORK(KPRP1P)+WORK(KPRP1M)),' / ',
     *            HALF*(WORK(KPRP1P)-WORK(KPRP1M))

              IF      (ISAPROP .EQ. +1) THEN
                WRITE(LUPRI,'(1X,2A,/1X,2A)')
     *            'the symmetric contribution corresponds to ',
     *            'the (real) physical result,',
     *            'the antisymmetric contribution is an artifact of ',
     *            'the non-symmetric CC parametrization.'
              ELSE IF (ISAPROP .EQ. -1) THEN
                WRITE(LUPRI,'(1X,2A,/1X,2A)')
     *            'the antisymmetric contribution corresponds to ',
     *            'the imaginary part of the physical result,',
     *            'the symmetric contribution is an artifact of ',
     *            'the non-symmetric CC parametrization.'
              ENDIF

            ELSE

              IF (SHIELD) THEN
                WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8,2X,F15.8)')
     *          '<<',LABELA,',',LABELB,
     *          '>>(',BLRFR(IFREQ),') =',RESULT,FACTOR*RESULT
              ELSE
                IF (ISYMA.EQ.ISYMB) THEN
                  WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8)')
     *            '<<',LABELA,',',LABELB,
     *            '>>(',BLRFR(IFREQ),') =',RESULT
                ELSE
                  WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F9.6,A,F15.8,A)')
     *            '<<',LABELA,',',LABELB,
     *            '>>(',BLRFR(IFREQ),') =',RESULT,' BY SYMMETRY !'
                ENDIF
              END IF

            ENDIF
            CALL WRIPRO(RESULT,MODELP,2,
     *                  LABELA,LABELB,LABELA,LABELB,
     *                  BLRFR(IFREQ),BLRFR(IFREQ),BLRFR(IFREQ),ISYMAB,
     *                  0,0,0)
         END DO

         IF (DIPPOL.AND.(LABELA(2:7).EQ.'DIPLEN')
     *             .AND.(LABELB(2:7).EQ.'DIPLEN')) THEN
           DO 6000 IFREQ = 1, NBLRFR
             KPRP1P = KPOL +          NBLRFR*(IOPER-1) + IFREQ - 1
             KPRP1M = KPOL + NLRPRP + NBLRFR*(IOPER-1) + IFREQ - 1

             RESULT = HALF*( WORK(KPRP1P) + SIGN * WORK(KPRP1M) )
             ERROR  = HALF*( WORK(KPRP1P) - SIGN * WORK(KPRP1M) )

             KPOLOF = KPOL2 + 3*3*(IFREQ-1) - 1

             IF (LABELA(1:2).EQ.'XD') IADR1 = 1
             IF (LABELA(1:2).EQ.'YD') IADR1 = 2
             IF (LABELA(1:2).EQ.'ZD') IADR1 = 3
             IF (LABELB(1:2).EQ.'XD') IADR2 = 1
             IF (LABELB(1:2).EQ.'YD') IADR2 = 2
             IF (LABELB(1:2).EQ.'ZD') IADR2 = 3
             IPOL = KPOLOF + 3*(IADR2-1) + IADR1
             WORK(IPOL) = RESULT
 6000      CONTINUE
         ENDIF
 4000 CONTINUE
C
C---------------------------------
C     Perform analysis for DIPPOL.
C---------------------------------
C
      IF (DIPPOL) THEN
         DO 9000 IFREQ = 1, NBLRFR
            KPOLI = KPOL2 + 3*3*(IFREQ-1)
            CALL DSCAL(9,XMONE,WORK(KPOLI),1)
            CALL CC_POLPRI(WORK(KPOLI),BLRFR(IFREQ))
 9000    CONTINUE
      ENDIF
C
C-------------
      CALL QEXIT('CC_LR')
      RETURN
      END
c*DECK CC_EATB
      SUBROUTINE CC_EATB(LABELA,ISYMA,FREQA,LRLXA,LPDBSA,
     *                   LABELB,ISYMB,FREQB,LRLXB,LPDBSB,
     *                   PRP,WORK,LWORK)
C
C----------------------------------------------------------------------
C
C   Purpose: Calculate etaA*tB contribution to second order properties.
C
C
C   Written by Ove Christiansen 21-6-1996
C   New version november 1996.
C
C----------------------------------------------------------------------
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
#include <dummy.h>
C
      PARAMETER( TWO = 2.0D00,TOLFRQ=1.0D-08 )
      DIMENSION WORK(LWORK)
      CHARACTER LABELA*8,LABELB*8,MODEL*10
      LOGICAL LRLXA, LRLXB, LPDBSA, LPDBSB
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_EATB: Calculating polarizabilty ')
      ENDIF
C
C------------------------
C     Allocate workspace.
C------------------------
C
      IF (ISYMA .NE. ISYMB ) CALL QUIT('Symmetry mismatch in CC_EATB')
      NTAMPB = NT1AM(ISYMB) + NT2AM(ISYMB)
      IF ( CCS ) NTAMPB = NT1AM(ISYMB)
      NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
      IF ( CCS ) NTAMPA = NT1AM(ISYMA)
C
      KETA  = 1
      KEND1 = KETA  + NTAMPA
      LEND1 = LWORK - KEND1

      KETA1 = KETA
      KETA2 = KETA1 + NT1AM(ISYMA)
C
      KR1   = KEND1
      KEND2 = KR1   + NTAMPB
      LEND2 = LWORK - KEND2
C
      IF (LEND2 .LT. 0)
     *      CALL QUIT('Insufficient space for allocation in CC_EATB')
C
C----------------------------------------------
C     Calculate contribution to polarizability.
C----------------------------------------------
C
      IF (LRLXA .OR. LPDBSA) THEN
         ILSTETA = IETA1(LABELA,LRLXA,FREQA,ISYMA)
         IOPT    = 3
         CALL CC_RDRSP('X1 ',ILSTETA,ISYMA,IOPT,MODEL,
     *                 WORK(KETA1),WORK(KETA2))
         IF (DEBUG) THEN
            WRITE (LUPRI,*) 'IETA1:',ILSTETA
            WRITE (LUPRI,*) 'norm(eta1):',
     *         DDOT(NT1AM(ISYMA),WORK(KETA1),1,WORK(KETA1),1)
            WRITE (LUPRI,*) 'norm(eta2):',
     *         DDOT(NT2AM(ISYMA),WORK(KETA2),1,WORK(KETA2),1)
         END IF
      ELSE
         CALL CC_ETAC(ISYMA,LABELA,WORK(KETA),'L0',1,0,
     *                DUMMY,WORK(KEND1),LEND1)
      END IF
C
      KR11 = KR1
      KR12 = KR1 + NT1AM(ISYMB)
      ILSTNR = IR1TAMP(LABELB,LRLXB,FREQB,ISYMB)
      IOPT   = 3
      CALL CC_RDRSP('R1 ',ILSTNR,ISYMB,IOPT,MODEL,WORK(KR11),
     *              WORK(KR12))
      IF (IPRINT .GT. 40 ) THEN
         CALL AROUND( 'In CC_EATB:  RSP vector ' )
         CALL CC_PRP(WORK(KR1),WORK(KR1+NT1AM(ISYMB)),ISYMB,1,1)
      ENDIF
      EATBCN = DDOT(NTAMPA,WORK(KETA),1,WORK(KR1),1)
C
      IF ( IPRINT .GT. 9 ) THEN
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
      IF (IPRINT .GT. 2 ) THEN
          WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F10.6,A,F14.10)')
     *    '<<',LABELA,',',LABELB,'>>(',
     *    FREQB,') EtaA*tB cont. = ',EATBCN
      ENDIF
      PRP  = EATBCN + PRP
C
      RETURN
      END
c*DECK CC_FABCON
      SUBROUTINE CC_FABCON(LABELA,ISYMA,FREQA,LRLXA,
     *                     LABELB,ISYMB,FREQB,LRLXB,
     *                     PRP,WORK,LWORK)
C
C----------------------------------------------------------------------
C
C     Purpose: Calculate F*TA(-omeg)*TB(omeg)
C
C     Written by Ove Christiansen 21-6-1996
C     New version 7-11-1996
C
C----------------------------------------------------------------------
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
#include <leinf.h>
C
      PARAMETER( TWO = 2.0D00,HALF=0.5D00,TOLFRQ=1.0D-08 )
      DIMENSION WORK(LWORK)
      CHARACTER LABELA*8,LABELB*8,MODEL*10
      LOGICAL LRLXA,LRLXB
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_FABCON: Calculating polarizabilty F-cont.')
      ENDIF
C
      NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
      IF ( CCS ) NTAMPA = NT1AM(ISYMA)
      NTAMPB = NT1AM(ISYMB) + NT2AM(ISYMB)
      IF ( CCS ) NTAMPB = NT1AM(ISYMB)
      IF (ISYMA .NE. ISYMB ) CALL QUIT('Symmetry mismatch in CC_FABCON')
C
C-----------------------------------------------
C     Loop perturbations of this symmetry class.
C-----------------------------------------------
C
      KR1   = 1
      KEND1 = KR1 + NTAMPB
      LEND1 = LWORK - KEND1
      IF (LEND1.LT. 0 )
     &     CALL QUIT(' TOO LITTLE WORKSPACE IN CC_FABCON-1 ')
C
C------------------------------
C     Get F-transformed vector.
C------------------------------
C
      KR11 = KR1
      KR12 = KR1 + NT1AM(ISYMB)
      ILSTNR = IR1TAMP(LABELB,LRLXB,FREQB,ISYMB)
      IOPT   = 3
      CALL CC_RDRSP('F1',ILSTNR,ISYMB,IOPT,MODEL,WORK(KR11),
     *              WORK(KR12))
      IF (IPRINT .GT. 40 ) THEN
         CALL AROUND( 'In CC_EATB:  F*RSP vector ' )
         CALL CC_PRP(WORK(KR1),WORK(KR1+NT1AM(ISYMB)),ISYMB,1,1)
      ENDIF
C
      IF ( DEBUG ) THEN
         XLV  = DDOT(NTAMPB, WORK(KR1),1,WORK(KR1),1)
         WRITE(LUPRI,1) 'Norm of F_Response vector:         ',XLV
      ENDIF
C
      KR2   = KEND1
      KEND2 = KR2 + NTAMPA
      LEND2 = LWORK - KEND2
      IF (LEND2.LT. 0 )
     &     CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ABFCON-2 ')
C
C-----------------------------------------------------------
C     Get response vectors and do the dot with the F*vector.
C-----------------------------------------------------------
C
      KR21 = KR2
      KR22 = KR2 + NT1AM(ISYMA)
      ILSTNR = IR1TAMP(LABELA,LRLXA,FREQA,ISYMA)
      IOPT   = 3
      CALL CC_RDRSP('R1',ILSTNR,ISYMA,IOPT,MODEL,WORK(KR21),
     *              WORK(KR22))
      IF ( DEBUG ) THEN
         XLV  = DDOT(NTAMPA, WORK(KR2),1,WORK(KR2),1)
         WRITE(LUPRI,1) 'Norm of Response vector:         ',XLV
      ENDIF
C
      FABCON = DDOT(NTAMPA,WORK(KR1),1,WORK(KR2),1)
      IF ( IPRINT .GT. 9 ) THEN
         WRITE(LUPRI,*) ' Singles contribution:',
     *      DDOT(NT1AM(ISYMA),WORK(KR1),1,WORK(KR2),1)
         IF (.NOT. CCS) WRITE(LUPRI,*) ' Doubles contribution:',
     *      DDOT(NT2AM(ISYMA),WORK(KR1+NT1AM(ISYMA)),1,
     *      WORK(KR2+NT1AM(ISYMA)),1)
      ENDIF
      IF (IPRINT .GT. 2 ) THEN
         WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F10.6,A,F14.10)')
     *   '<<',LABELA,',',LABELB,'>>(',
     *   FREQB,') F*tA*tB cont. = ',FABCON
      ENDIF
      PRP       = PRP       + FABCON
C
   1  FORMAT(1x,A35,1X,E20.10)
      RETURN
      END
c*DECK CC_LAKSIB
      SUBROUTINE CC_LAKSIB(LABELA,ISYMA,FREQA,LRLXA,
     *                     LABELB,ISYMB,FREQB,LRLXB,
     *                     PRP,WORK,LWORK)
C
C----------------------------------------------------------------------
C
C   Purpose: Calculate LD*ksiC contribution to second order properties.
C            For use in calculation of molecular properties from
C            Asymmetric formulaes not in accordance with 2n+2 rule for
C            the multipliers, left vector, t-bar, lamdas, zeta or
C            whatever your preferred choice is today.
C
C     Written by Ove Christiansen 17-10-1996/7-11-1996
C
C----------------------------------------------------------------------
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
C
      PARAMETER( TWO = 2.0D00,TOLFRQ=1.0D-08 )
      DIMENSION WORK(LWORK)
      CHARACTER LABELA*8,LABELB*8,MODEL*10
      LOGICAL LRLXA, LRLXB
C
      IF ( IPRINT .GT. 5 ) THEN
         CALL AROUND( 'IN CC_LAKSIB: Calculating polarizabilty '
     *                 //'contribution')
         WRITE(LUPRI,'(/,1x,A,F16.8,/,A,I2,/,3A,/,A,2L3)')
     *   'Calculating response property with frequency',FREQB,
     *   ' Operator symmetry = ',ISYMB,
     *   ' Labels = ',LABELA, LABELB,
     *   ' orbital relaxation flags = ',LRLXA, LRLXB
      ENDIF
C
C------------------------
C     Allocate workspace.
C------------------------
C
      NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
      IF ( CCS ) NTAMPA = NT1AM(ISYMA)
      NTAMPB = NT1AM(ISYMB) + NT2AM(ISYMB)
      IF ( CCS ) NTAMPB = NT1AM(ISYMB)
      IF (ISYMA .NE. ISYMB ) CALL QUIT('Symmetry mismatch in CC_LAKSIB')
C
      KKSI = 1
      KEND1 = KKSI + NTAMPA
      LEND1 = LWORK - KEND1

      KKSI1 = KKSI
      KKSI2 = KKSI1 + NT1AM(ISYMA)
C
      KR1   = KEND1
      KEND2 = KR1   + NTAMPB
      LEND2 = LWORK - KEND2
C
      IF (LEND2 .LT. 0)
     *      CALL QUIT('Insufficient space for allocation in CC_LAKSIB')
C
C----------------------------------------------
C     Calculate contribution to polarizability.
C----------------------------------------------
C
      IF (LRLXA) THEN
         ILSTRHS = IRHSR1(LABELA,LRLXA,FREQA,ISYMA)
         IOPT    = 3
         CALL CC_RDRSP('O1 ',ILSTRHS,ISYMA,IOPT,MODEL,
     *                 WORK(KKSI1),WORK(KKSI2))
         IF (DEBUG) THEN
            WRITE (LUPRI,*) 'IRHSR1:',ILSTRHS
            WRITE (LUPRI,*) 'norm(xksi1):',
     *         DDOT(NT1AM(ISYMA),WORK(KKSI1),1,WORK(KKSI1),1)
            WRITE (LUPRI,*) 'norm(xksi2):',
     *         DDOT(NT2AM(ISYMA),WORK(KKSI2),1,WORK(KKSI2),1)
            call cc_prp(work(kksi1),work(kksi2),isyma,1,1)
         END IF
      ELSE
         CALL CC_XKSI(WORK(KKSI),LABELA,ISYMA,0,DUMMY,WORK(KEND1),LEND1)
      END IF
C
      KR11 = KR1
      KR12 = KR1 + NT1AM(ISYMB)
      ILSTNR = IL1ZETA(LABELB,LRLXB,FREQB,ISYMB)
      IOPT   = 3
      CALL CC_RDRSP('L1',ILSTNR,ISYMB,IOPT,MODEL,WORK(KR11),
     *              WORK(KR12))
      ABCON = DDOT(NTAMPA,WORK(KR1),1,WORK(KKSI),1)
      IF ( DEBUG ) THEN
         XLV  = DDOT(NTAMPB, WORK(KR1),1,WORK(KR1),1)
         WRITE(LUPRI,1) 'Norm of Response vector:         ',XLV
      ENDIF
C
      IF ( IPRINT .GT. 9 ) THEN
          WRITE(LUPRI,*) ' Singles contribution:',
     *       DDOT(NT1AM(ISYMA),WORK(KKSI),1,WORK(KR1),1)
          IF (.NOT. CCS) WRITE(LUPRI,*) ' Doubles contribution:',
     *       DDOT(NT2AM(ISYMA),WORK(KKSI+NT1AM(ISYMA)),1,
     *       WORK(KR1+NT1AM(ISYMA)),1)
      ENDIF
      IF (IPRINT .GT. 2 ) THEN
         WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F10.6,A,F10.6)')
     *   '<<',LABELA,',',LABELB,'>>(',
     *   FREQB,') LB*ksiA cont. = ',ABCON
      ENDIF
      PRP  = PRP + ABCON
C
   1  FORMAT(1x,A35,1X,E20.10)
      RETURN
      END
c*DECK CC_POLPRI
      SUBROUTINE CC_POLPRI(POL,FRQ)
C
C----------------------------------------------------------------------
C
C   Purpose: Calculate LD*ksiC contribution to second order properties.
C            For use in calculation of molecular properties from
C            Asymmetric formulaes not in accordance with 2n+2 rule for
C            the multipliers, left vector, t-bar, lamdas, zeta or
C            whatever your preferred choice is today.
C
C     Written by Ove Christiansen 17-10-1996/7-11-1996
C
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
      PARAMETER (TOLFRQ = 1.0D-08,ONE= 1.0D0,THR = 1.0D-08)
      PARAMETER (DPOLAUTSI = 1.648778D-41, QPOLAUTSI = 4.617048 D-62 )
C
C DPOL C2m2J-1, QPOL C2m4J-1
C
#include <iratdef.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <cclrinf.h>
#include <ccrspprp.h>
C
      DIMENSION POL(*),PVAL(3),PAXIS(3,3)
      CHARACTER MODEL*10
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_POLPRI: Output polarizabilities   ' )
      ENDIF
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CIS) MODEL = 'CIS       '
      IF (CC2) MODEL = 'CC2       '
C
      IF (.NOT.(CCS.OR.CC2.OR.CCSD)) THEN
         WRITE(LUPRI,'(A)')
     &        ' CC_POLPRI: Do not want to calculate anything'
     *                  //' else than CCS, CC2 and CCSD properties '
         CALL QUIT('Model not CCS, CC2, or CCSD in CC_POLPRI')
      ENDIF
C
C--------------------------------------
C        Find the frequency components.
C--------------------------------------
C
         WRITE(LUPRI,'(//,1X,A6,A,F10.6,/)') MODEL(1:6),
     *       'polarizability for frequency: ',FRQ
         CALL OUTPUT(POL,1,3,1,3,3,3,1,LUPRI)
cmbh: print polarizability for MidasCpp
         call wripro(POL(1),'  '//MODEL(1:6)//'  ',2,
     *               'X_DIPLEN','X_DIPLEN','X_DIPLEN','X_DIPLEN',
     *               FRQ,FRQ,FRQ,1,0,0,0) 
         call wripro(POL(2),'  '//MODEL(1:6)//'  ',2,
     *               'X_DIPLEN','Y_DIPLEN','X_DIPLEN','Y_DIPLEN',
     *               FRQ,FRQ,FRQ,1,0,0,0) 
         call wripro(POL(3),'  '//MODEL(1:6)//'  ',2,
     *               'X_DIPLEN','Z_DIPLEN','X_DIPLEN','Z_DIPLEN',
     *               FRQ,FRQ,FRQ,1,0,0,0) 
         call wripro(POL(5),'  '//MODEL(1:6)//'  ',2,
     *               'Y_DIPLEN','Y_DIPLEN','Y_DIPLEN','Y_DIPLEN',
     *               FRQ,FRQ,FRQ,1,0,0,0) 
         call wripro(POL(6),'  '//MODEL(1:6)//'  ',2,
     *               'Y_DIPLEN','Z_DIPLEN','Y_DIPLEN','Z_DIPLEN',
     *               FRQ,FRQ,FRQ,1,0,0,0) 
         call wripro(POL(9),'  '//MODEL(1:6)//'  ',2,
     *               'Z_DIPLEN','Z_DIPLEN','Z_DIPLEN','Z_DIPLEN',
     *               FRQ,FRQ,FRQ,1,0,0,0) 
cmbh end
C
         CALL  TNSRAN(POL,PVAL,PAXIS,
     *                ALFSQ,BETSQ,ITST,ITST2,
     *                APAR,APEN,XKAPPA,IPAR)
         WRITE(LUPRI,'(/,1X,A38,F14.6)')
     *              'Alfa**2 Invariant:            '
     *            //'            ',ALFSQ
         WRITE(LUPRI,'(1X,A38,F14.6)')
     *           'Beta**2 Invariant:            '
     *            //'            ',BETSQ
         SHPAL = SQRT(ALFSQ)
         ANINV = SQRT(BETSQ)
         WRITE(LUPRI,'(/,1X,A42,F10.6,A)') 'Isotropic Polarizability: '
     *         //'                 ',SHPAL,' a.u.'
         WRITE(LUPRI,'(1X,A42,F10.6,A)') 'Polarizability anisotropy '
     *      //'invariant:      ',ANINV,' a.u.'
         IF (ITST .EQ. 0) THEN
          IF (ITST2 .EQ. 3) THEN
             WRITE(LUPRI,'(/,1X,A)')
     *           'Polarizability has spherical symmetry:'
             WRITE(LUPRI,'(1X,A,F10.6,A,3X,E15.6,A)')
     *   'Isotropic polarizabilty: ',APAR,' a.u.',APAR*DPOLAUTSI,' S.I.'
          ELSE IF (ITST2 .EQ. 1) THEN
             WRITE(LUPRI,'(/,1X,A,/)')
     *   'Polarizability has cylinder symmetry: '
             WRITE(LUPRI,'(1X,A,F10.6,A,3X,E15.6,A)')
     *   'Parallel component:      ',APAR,' a.u.',APAR*DPOLAUTSI,' S.I.'
             WRITE(LUPRI,'(1X,A,F10.6,A,3X,E15.6,A)')
     *   'Perpendicular component: ',APEN,' a.u.',APEN*DPOLAUTSI,' S.I.'
             WRITE(LUPRI,'(/,1X,A42,F12.6)')
     *   'Dimensionless polarizability anisotropy:  ',XKAPPA
          ELSE IF (ITST2. EQ. 0) THEN
             WRITE(LUPRI,'(/,1X,A,/)')
     *          'Polarizability is diagonal with diagonal values:   '
                WRITE(LUPRI,'(1X,A)')
     *      '        a.u.          S.I. '
                WRITE(LUPRI,'(1X,A,F10.6,3X,E15.6)')
     *          'XX  ',PVAL(1),PVAL(1)*DPOLAUTSI
                WRITE(LUPRI,'(1X,A,F10.6,3X,E15.6)')
     *          'YY  ',PVAL(2),PVAL(2)*DPOLAUTSI
                WRITE(LUPRI,'(1X,A,F10.6,3X,E15.6)')
     *          'ZZ  ',PVAL(3),PVAL(3)*DPOLAUTSI
          ENDIF
         ELSE
             WRITE(LUPRI,'(/,1X,A,/)')
     *           'Principal values of diagonalized Polarizability:'
             WRITE(LUPRI,'(1X,A)')
     *      '        a.u.          S.I. '
             WRITE(LUPRI,'(1X,A,F10.6,3X,E15.6)')
     *          '1     ',PVAL(1),PVAL(1)*DPOLAUTSI
             WRITE(LUPRI,'(1X,A,F10.6,3X,E15.6)')
     *          '2     ',PVAL(2),PVAL(2)*DPOLAUTSI
             WRITE(LUPRI,'(1X,A,F10.6,3X,E15.6)')
     *          '3     ',PVAL(3),PVAL(3)*DPOLAUTSI
             WRITE(LUPRI,'(/,1X,A,/)')
     *           'Principal axis of diagonalized Polarizability:'
             CALL OUTPUT(POL,1,3,1,3,3,3,1,LUPRI)
         ENDIF
         WRITE(LUPRI,'(/,1X,A,E18.8,A,/)')
     *      'Conversion factor (a.u. - S.I.):',DPOLAUTSI,' (C^2m^2J^-1)'
C
            CALL WRIPRO(SHPAL,MODEL,2,
     *                  'isoalpha','isoalpha','isoalpha','isoalpha',
     *                  FRQ,DUMMY,DUMMY,1,0,0,0)
C
            CALL WRIPRO(ANINV,MODEL,2,
     *                  'anis_inv','anis_inv','anis_inv','anis_inv',
     *                  FRQ,DUMMY,DUMMY,1,0,0,0)
C
      END
c*DECK CC_LRESID
       SUBROUTINE CC_LRESID(WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              linear response residue calculation.
C
C              CCS, CC2, CCSD
C
C     Modified version for general linear response properties
C     Ove Christiansen November 1996.
C
C     Symmetrization (C+/-w operator)
C     Thomas Bondo Pedersen, January 2005.
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
C
      LOGICAL LCALC
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*10,MODELP*10
      CHARACTER LABELA*8,LABELB*8
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
#include <leinf.h>
C
#include <mxcent.h>
#include <maxaqn.h>
#include <symmet.h>
#include <codata.h>
C
      PARAMETER (RAUSI = FPEPS0*(HBAR**3)*1.0D55/(EMASS**2))
      PARAMETER (RAUCGS = ECHARGE*ECHARGE*HBAR*XTANG*CCM*1.0D36/EMASS)
C
      TIMTOT = SECOND()
      NTOT   = 0
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
     *'*<<<<<<<<<< OUTPUT FROM COUPLED CLUSTER LINEAR RESPONSE >'//
     *'>>>>>>>>>*'
      IF ( OSCSTR ) THEN
         WRITE (LUPRI,'(1X,A)')
     *   '*                                                        '//
     *   '         *'
         WRITE (LUPRI,'(1X,A)')
     *   '*<<<<<<<<<<      CALCULATION OF CC OSCILLATOR STRENGTHS  >'//
     *   '>>>>>>>>>*'
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
         MODEL = 'CC2       '
      ENDIF
      IF (CCS) THEN
         MODEL = 'CCS       '
      ENDIF
      IF (CC3  ) THEN
         MODEL = 'CC3       '
         WRITE(LUPRI,'(/,1x,A)')
     *    'CC3 Oscillator strengths not implemented yet'
         RETURN
      ENDIF
      IF (CC1A) THEN
         MODEL = 'CCSDT-1a  '
         WRITE(LUPRI,'(/,1x,A)')
     *    'CC1A Oscillator strengths not implemented yet'
         RETURN
      ENDIF
      IF (CCSD) THEN
         MODEL = 'CCSD      '
      ENDIF
C
      IF (CIS) THEN
         MODELP = 'CIS       '
      ELSE
         MODELP = MODEL
      ENDIF
C
      CALL AROUND( 'Calculation of '//MODELP// ' residues ')
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_LRESID Workspace:',LWORK
C
C-------------------------------------------------------------------------
C     Calculate polarizabilities in loops over symmetries and frequencies.
C-------------------------------------------------------------------------
C
      CALL FLSHFO(LUPRI)
C
      NALRPRP = NLRSOP*NXLRSST
      NBLRPRP = NLRSOP*NXLRSST
C
      KOSCS    = 1
      KOSCSF   = KOSCS  + NALRPRP
      KSYMB    = KOSCSF + NBLRPRP
      KSYMA    = KSYMB  + NBLRPRP
      KEND1    = KSYMA  + NALRPRP
      LEND1    = LWORK  - KEND1
C
      IF (LEND1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CC_LRESID [1]')
      END IF
C
      CALL DZERO(WORK(KOSCS),NALRPRP)
      CALL DZERO(WORK(KOSCSF),NBLRPRP)
      CALL DZERO(WORK(KSYMB),NBLRPRP)
      CALL DZERO(WORK(KSYMA),NALRPRP)
C
C----------------------------------------------
C     Loop over states and operators requested.
C----------------------------------------------
C
      DO 1000 IRSD  = 1, NXLRSST
        ISTATE = ILRSST(IRSD)
        ISYME  = ISYEXC(ISTATE)
        ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
        EIGV   = EIGVAL(ISTATE)
        IF (IPRINT .GT. 5) THEN
          WRITE(LUPRI,'(/,1x,A,I3,/1X,A,I3,A,F16.8)')
     *    'Calculating linear response residues for state',ISTSY,
     *    'of symmetry ',ISYME,' and with eigenvalue: ',EIGV
        ENDIF
C
        DO 2000 IOPER = 1, NLRSOP
          ISYMA  = ISYOPR(IALRSOP(IOPER))
          ISYMB  = ISYOPR(IBLRSOP(IOPER))

          IF ((ISYME.EQ.ISYMA).AND.(ISYME.EQ.ISYMB)) THEN

            LABELA = LBLOPR(IALRSOP(IOPER))
            LABELB = LBLOPR(IBLRSOP(IOPER))
C
C----------------------------------------
C           Calculate transition moments.
C----------------------------------------
C
            KRES1 = KOSCS  + NLRSOP*(IRSD-1) + IOPER - 1
            KRES2 = KOSCSF + NLRSOP*(IRSD-1) + IOPER - 1
            CALL CC_LRSD(LABELA,ISYMA,
     *                   LABELB,ISYMB,
     *                   ISTATE,WORK(KRES1),WORK(KRES2),
     *                   WORK(KEND1),LEND1)
C
            KRES3 = KSYMB + NLRSOP*(IRSD-1) + IOPER - 1
            KRES4 = KSYMA + NLRSOP*(IRSD-1) + IOPER - 1
            IF (LABELA .EQ. LABELB) THEN
               WORK(KRES3) = WORK(KRES1)
               WORK(KRES4) = WORK(KRES2)
            ELSE
               CALL CC_LRSD(LABELB,ISYMB,
     *                      LABELA,ISYMA,
     *                      ISTATE,WORK(KRES3),WORK(KRES4),
     *                      WORK(KEND1),LEND1)
            END IF
            IF (LOCDBG) THEN
               WRITE(LUPRI,*) ' Residue symmetrization:'
               WRITE(LUPRI,*) '   Exc. state: ',ISTSY,' of sym. ',
     &                            ISYME,':'
               WRITE(LUPRI,*) '   T(0f,',LABELA,') = ',WORK(KRES1)
               WRITE(LUPRI,*) '   T(f0,',LABELB,') = ',WORK(KRES2)
               WRITE(LUPRI,*) '   T(0f,',LABELB,') = ',WORK(KRES3)
               WRITE(LUPRI,*) '   T(f0,',LABELA,') = ',WORK(KRES4)
               CALL FLSHFO(LUPRI)
            END IF
C
          ENDIF
 2000   CONTINUE
 1000 CONTINUE
C
C-----------------------------------------
C     Output Linear response properties.
C     Save requested transition strengths.
C-----------------------------------------
C
      IF (OSCSTR) THEN
         LOSCIL = NEXCI*3*3
      ELSE
         LOSCIL = 0
      END IF
C
      IF (VELSTR) THEN
         LOSCIV = NEXCI*3*3
      ELSE
         LOSCIV = 0
      END IF
C
      IF (MIXSTR) THEN
         LOSCIM = NEXCI*3*3
      ELSE
         LOSCIM = 0
      END IF
C
      IF (ROTLEN) THEN
         LROTL = NEXCI*3
         LCHKL = NEXCI
      ELSE
         LROTL = 0
         LCHKL = 0
      ENDIF
C
      IF (ROTVEL) THEN
         LROTV = NEXCI*3
         LCHKV = NEXCI
      ELSE
         LROTV = 0
         LCHKV = 0
      ENDIF
C
      IF (RTNLEN) THEN
         LRQL = NEXCI*3*9
         LRML = NEXCI*3*3
         NWRL = 0
      ELSE
         LRQL = 0
         LRML = 0
      ENDIF
C
      IF (RTNVEL) THEN
         LRQV = NEXCI*3*9
         LRMV = NEXCI*3*3
         NWRV = 0
      ELSE
         LRQV = 0
         LRMV = 0
      ENDIF
C
      KOSCS2 = KEND1
      KTRS   = KOSCS2  + LOSCIL
      KVELST = KTRS    + LOSCIL
      KVELST2= KVELST  + LOSCIV
      KMIXST = KVELST2 + LOSCIV
      KMIXST2= KMIXST  + LOSCIM
      KROTL  = KMIXST2 + LOSCIM
      KROTV  = KROTL   + LROTL
      KRQL   = KROTV   + LROTV
      KRML   = KRQL    + LRQL
      KRQL2  = KRML    + LRML
      KRML2  = KRQL2   + LRML
      KRQV   = KRML2   + LRML
      KRMV   = KRQV    + LRQV
      KRQV2  = KRMV    + LRMV
      KRMV2  = KRQV2   + LRMV
      KCHKL  = KRMV2   + LRMV
      KCHKV  = KCHKL   + LCHKL
      KEND2  = KCHKV   + LCHKV
      LEND2  = LWORK   - KEND2
C
      IF (LEND2 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CC_LRESID [2]')
      END IF
C
      IF (OSCSTR) THEN
         CALL DZERO(WORK(KOSCS2),LOSCIL)
         CALL DZERO(WORK(KTRS),LOSCIL)
      END IF
      IF (VELSTR) THEN
         CALL DZERO(WORK(KVELST),LOSCIV)
         CALL DZERO(WORK(KVELST2),LOSCIV)
      END IF
      IF (MIXSTR) THEN
         CALL DZERO(WORK(KMIXST),LOSCIM)
         CALL DZERO(WORK(KMIXST2),LOSCIM)
      END IF
      IF (ROTLEN) THEN
         CALL DZERO(WORK(KROTL),LROTL)
         CALL DZERO(WORK(KROTL),LROTL)
         CALL DZERO(WORK(KCHKL),LCHKL)
      END IF
      IF (ROTVEL) THEN
         CALL DZERO(WORK(KROTV),LROTV)
         CALL DZERO(WORK(KROTV),LROTV)
         CALL DZERO(WORK(KCHKV),LCHKV)
      END IF
      IF (RTNLEN) THEN
         CALL DZERO(WORK(KRQL),LRQL)
         CALL DZERO(WORK(KRML),LRML)
         CALL DZERO(WORK(KRQL2),LRML)
         CALL DZERO(WORK(KRML2),LRML)
      END IF
      IF (RTNVEL) THEN
         CALL DZERO(WORK(KRQV),LRQV)
         CALL DZERO(WORK(KRMV),LRMV)
         CALL DZERO(WORK(KRQV2),LRMV)
         CALL DZERO(WORK(KRMV2),LRMV)
      END IF
C
      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     *  'Right transition moments in atomic units:'
      WRITE(LUPRI,'(1X,A,/)')
     *  '-----------------------------------------------'
C
      DO IOPER = 1, NLRSOP
        ISYMA  = ISYOPR(IALRSOP(IOPER))
        ISYMB  = ISYOPR(IBLRSOP(IOPER))
        LABELA = LBLOPR(IALRSOP(IOPER))
        LABELB = LBLOPR(IBLRSOP(IOPER))
        DO IRSD  = 1, NXLRSST
          ISTATE = ILRSST(IRSD)
          ISYME  = ISYEXC(ISTATE)
          ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
          EIGV   = EIGVAL(ISTATE)
          IF ((ISYME.EQ.ISYMA).AND.(ISYME.EQ.ISYMB)) THEN
            K1     = NLRSOP*(IRSD-1) + IOPER + KOSCS - 1
            WRITE(LUPRI,'(1X,I2,F15.6,2X,A1,A8,A6,1X,F15.8)')
     *      ISTATE,EIGV,'<',LABELA,'|f> = ',WORK(K1)
          ENDIF
        END DO
      END DO
C
      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     *  'Left  transition moments in atomic units:'
      WRITE(LUPRI,'(1X,A,/)')
     *  '-----------------------------------------------'
C
      DO IOPER = 1, NLRSOP
        ISYMA  = ISYOPR(IALRSOP(IOPER))
        ISYMB  = ISYOPR(IBLRSOP(IOPER))
        LABELA = LBLOPR(IALRSOP(IOPER))
        LABELB = LBLOPR(IBLRSOP(IOPER))
        DO IRSD  = 1, NXLRSST
          ISTATE = ILRSST(IRSD)
          ISYME  = ISYEXC(ISTATE)
          ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
          EIGV   = EIGVAL(ISTATE)
          IF ((ISYME.EQ.ISYMA).AND.(ISYME.EQ.ISYMB)) THEN
            K1     = NLRSOP*(IRSD-1) + IOPER + KOSCSF - 1
            WRITE(LUPRI,'(1X,I2,F15.6,2X,A3,A8,A4,1X,F15.8)')
     *         ISTATE,EIGV,'<f|',LABELB,'> = ',WORK(K1)
          ENDIF
        END DO
      END DO
C
      CALL FLSHFO(LUPRI)
C
C----------------------------------------------------------------
C     Calculate linear response residues from transition moments,
C     incl. symmetrization.
C----------------------------------------------------------------
C
      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     *  'linear response residue property in atomic units:'
C
      WRITE(LUPRI,'(1X,A,/)')
     *  '-------------------------------------------------------'
C
      DO IOPER = 1, NLRSOP
        ISYMA  = ISYOPR(IALRSOP(IOPER))
        ISYMB  = ISYOPR(IBLRSOP(IOPER))
        LABELA = LBLOPR(IALRSOP(IOPER))
        LABELB = LBLOPR(IBLRSOP(IOPER))
        DO IRSD  = 1, NXLRSST
          ISTATE = ILRSST(IRSD)
          ISYME  = ISYEXC(ISTATE)
          ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
          EIGV   = EIGVAL(ISTATE)
          ISYMEA = MULD2H(ISYME,ISYMA)
          IF ((ISYME.EQ.ISYMA).AND.(ISYME.EQ.ISYMB)) THEN
            NTOT = NTOT + 1
            K1 = NLRSOP*(IRSD-1) + IOPER + KOSCS - 1
            K2 = NLRSOP*(IRSD-1) + IOPER + KOSCSF - 1
            K3 = NLRSOP*(IRSD-1) + IOPER + KSYMB  - 1
            K4 = NLRSOP*(IRSD-1) + IOPER + KSYMA  - 1
            IHERMA = ISYMAT(IALRSOP(IOPER))
            IHERMB = ISYMAT(IBLRSOP(IOPER))
            ISASB  = IHERMA*IHERMB
            IF (ISASB .EQ. 0) THEN
               WRITE(LUPRI,*) ' WARNING: operators ',LABELA,LABELB,
     &                        ' have undefined hermiticities: ',
     &                          IHERMA,IHERMB
               WRITE(LUPRI,*) ' Residue not appropriately symmetrized..'
               CALL FLSHFO(LUPRI)
               SIGN = 1.0D0
            ELSE
               SIGN = DBLE(ISASB)
            ENDIF
            RESIDAB = WORK(K1)*WORK(K2)
            RESIDBA = WORK(K3)*WORK(K4)
            RESIDUE = 0.5D0*(RESIDAB + SIGN*RESIDBA)
            IF (RESIDUE.GE.0.0D0) THEN
              SQRRES=SQRT(RESIDUE)
            ELSE
              SQRRES=-SQRT(-RESIDUE)
            ENDIF
            WRITE(LUPRI,'(1X,A6,A8,A1,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     *      'RES{<<',LABELA,',',LABELB,'>>(',EIGV,')} =',
     *      RESIDUE,' ( ',SQRRES,')'
            IF (LOCDBG) THEN
               WRITE(LUPRI,*) '   A,B: ',RESIDAB,
     &                        ' B,A: ',RESIDBA,
     &                        ' combination: ',ISASB
            END IF
            IF (OSCSTR) THEN ! length gauge oscillator strength
              IF (LABELA(2:7).EQ.'DIPLEN' .AND.
     &            LABELB(2:7).EQ.'DIPLEN') THEN
                 IF (LABELA(1:1).EQ.'X') IADR1 = 1
                 IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                 IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                 IF (LABELB(1:1).EQ.'X') IADR2 = 1
                 IF (LABELB(1:1).EQ.'Y') IADR2 = 2
                 IF (LABELB(1:1).EQ.'Z') IADR2 = 3
                 IF ((IADR1+IADR2).GE.2) THEN
                   IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KOSCS2-1
                   WORK(IOSCS2) = RESIDUE
                 END IF
              END IF
            END IF
            IF (VELSTR) THEN ! velocity gauge oscillator strength
              IF (LABELA(2:7).EQ.'DIPVEL' .AND.
     &            LABELB(2:7).EQ.'DIPVEL') THEN
                 IF (LABELA(1:1).EQ.'X') IADR1 = 1
                 IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                 IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                 IF (LABELB(1:1).EQ.'X') IADR2 = 1
                 IF (LABELB(1:1).EQ.'Y') IADR2 = 2
                 IF (LABELB(1:1).EQ.'Z') IADR2 = 3
                 IF ((IADR1+IADR2).GE.2) THEN
                   IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KVELST-1
                   WORK(IOSCS2) = RESIDUE
                 END IF
              END IF
            END IF
            IF (MIXSTR) THEN ! Mixed gauge oscillator strength
              IF (LABELA(2:7).EQ.'DIPLEN' .AND.
     &            LABELB(2:7).EQ.'DIPVEL') THEN
                 IF (LABELA(1:1).EQ.'X') IADR1 = 1
                 IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                 IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                 IF (LABELB(1:1).EQ.'X') IADR2 = 1
                 IF (LABELB(1:1).EQ.'Y') IADR2 = 2
                 IF (LABELB(1:1).EQ.'Z') IADR2 = 3
                 IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KMIXST-1
                 WORK(IOSCS2) = RESIDUE
               END IF
            END IF
            IF (ROTLEN) THEN ! Length gauge rotatory strength
               IF (LABELA(2:7) .EQ. 'DIPLEN' .AND.
     &             LABELB(2:7) .EQ. 'ANGMOM') THEN
                 IF (LABELA(1:1).EQ.'X') IADR1 = 1
                 IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                 IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                 IF (LABELB(1:1).EQ.'X') IADR2 = 1
                 IF (LABELB(1:1).EQ.'Y') IADR2 = 2
                 IF (LABELB(1:1).EQ.'Z') IADR2 = 3
                 IF (IADR1 .EQ. IADR2) THEN
                    IROTST = KROTL + 3*(ISTATE-1) + IADR1 - 1
                    WORK(IROTST) = RESIDUE
                 END IF
               END IF
            END IF
            IF (ROTVEL) THEN ! Velocity gauge rotatory strength
               IF (LABELA(2:7) .EQ. 'DIPVEL' .AND.
     &             LABELB(2:7) .EQ. 'ANGMOM') THEN
                 IF (LABELA(1:1).EQ.'X') IADR1 = 1
                 IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                 IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                 IF (LABELB(1:1).EQ.'X') IADR2 = 1
                 IF (LABELB(1:1).EQ.'Y') IADR2 = 2
                 IF (LABELB(1:1).EQ.'Z') IADR2 = 3
                 IF (IADR1 .EQ. IADR2) THEN
                    IROTST = KROTV + 3*(ISTATE-1) + IADR1 - 1
                    WORK(IROTST) = RESIDUE
                 END IF
               END IF
            END IF
            IF (RTNLEN) THEN
               IF (LABELA(2:7) .EQ. 'DIPLEN') THEN
                  IF (LABELB(3:8) .EQ. 'SECMOM') THEN
                     IF (LABELA(1:1).EQ.'X') IADR1 = 1
                     IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                     IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                     IF (LABELB(1:2).EQ.'XX') THEN
                        IADR23 = 1
                        IADR32 = 0
                     ELSE IF (LABELB(1:2).EQ.'XY') THEN
                        IADR23 = 4
                        IADR32 = 2
                     ELSE IF (LABELB(1:2).EQ.'XZ') THEN
                        IADR23 = 7
                        IADR32 = 3
                     ELSE IF (LABELB(1:2).EQ.'YY') THEN
                        IADR23 = 5
                        IADR32 = 0
                     ELSE IF (LABELB(1:2).EQ.'YZ') THEN
                        IADR23 = 8
                        IADR32 = 6
                     ELSE IF (LABELB(1:2).EQ.'ZZ') THEN
                        IADR23 = 9
                        IADR32 = 0
                     END IF
                     IRTEN = KRQL + 3*9*(ISTATE-1)
     &                     + 3*(IADR23-1) + IADR1 - 1
                     WORK(IRTEN) = RESIDUE
                     IF (IADR32 .NE. 0) THEN
                        IRTEN = KRQL + 3*9*(ISTATE-1)
     &                        + 3*(IADR32-1) + IADR1 - 1
                        WORK(IRTEN) = RESIDUE
                     END IF
                  ELSE IF (LABELB(2:7) .EQ. 'ANGMOM') THEN
                     IF (LABELA(1:1).EQ.'X') IADR1 = 1
                     IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                     IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                     IF (LABELB(1:1).EQ.'X') IADR2 = 1
                     IF (LABELB(1:1).EQ.'Y') IADR2 = 2
                     IF (LABELB(1:1).EQ.'Z') IADR2 = 3
                     IRTEN = KRML + 3*3*(ISTATE-1)
     &                     + 3*(IADR2-1) + IADR1 - 1
                     WORK(IRTEN) = RESIDUE
                  END IF
               END IF
            END IF
            IF (RTNVEL) THEN
               IF (LABELA(2:7) .EQ. 'DIPVEL') THEN
                  IF (LABELB(3:8) .EQ. 'ROTSTR') THEN
                     IF (LABELA(1:1).EQ.'X') IADR1 = 1
                     IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                     IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                     IF (LABELB(1:2).EQ.'XX') THEN
                        IADR23 = 1
                        IADR32 = 0
                     ELSE IF (LABELB(1:2).EQ.'XY') THEN
                        IADR23 = 4
                        IADR32 = 2
                     ELSE IF (LABELB(1:2).EQ.'XZ') THEN
                        IADR23 = 7
                        IADR32 = 3
                     ELSE IF (LABELB(1:2).EQ.'YY') THEN
                        IADR23 = 5
                        IADR32 = 0
                     ELSE IF (LABELB(1:2).EQ.'YZ') THEN
                        IADR23 = 8
                        IADR32 = 6
                     ELSE IF (LABELB(1:2).EQ.'ZZ') THEN
                        IADR23 = 9
                        IADR32 = 0
                     END IF
                     IRTEN = KRQV + 3*9*(ISTATE-1)
     &                     + 3*(IADR23-1) + IADR1 - 1
                     WORK(IRTEN) = RESIDUE
                     IF (IADR32 .NE. 0) THEN
                        IRTEN = KRQV + 3*9*(ISTATE-1)
     &                        + 3*(IADR32-1) + IADR1 - 1
                        WORK(IRTEN) = RESIDUE
                     END IF
                  ELSE IF (LABELB(2:7) .EQ. 'ANGMOM') THEN
                     IF (LABELA(1:1).EQ.'X') IADR1 = 1
                     IF (LABELA(1:1).EQ.'Y') IADR1 = 2
                     IF (LABELA(1:1).EQ.'Z') IADR1 = 3
                     IF (LABELB(1:1).EQ.'X') IADR2 = 1
                     IF (LABELB(1:1).EQ.'Y') IADR2 = 2
                     IF (LABELB(1:1).EQ.'Z') IADR2 = 3
                     IRTEN = KRMV + 3*3*(ISTATE-1)
     &                     + 3*(IADR2-1) + IADR1 - 1
                     WORK(IRTEN) = RESIDUE
                  END IF
               END IF
            END IF
          ELSE
            RESIDUE = 0.0D0
            SQRRES  = 0.0D0
          ENDIF
          IF (LABELA.EQ.LABELB) THEN
             CALL WRIPRO(SQRRES,MODEL,-1,
     *                   LABELA,LABELB,LABELA,LABELB,
     *                   EIGV,EIGV,EIGV,ISYMEA,ISYME,1,ISTATE)
          ENDIF
        END DO
      END DO

C
C-----------------------------------------------
C     Perform analysis for oscillator strengths.
C-----------------------------------------------
C
      IF (OSCSTR) CALL DCOPY(LOSCIL,WORK(KOSCS2),1,WORK(KTRS),1)
      IF (VELSTR) CALL DCOPY(LOSCIV,WORK(KVELST),1,WORK(KVELST2),1)
      IF (MIXSTR) CALL DCOPY(LOSCIM,WORK(KMIXST),1,WORK(KMIXST2),1)
C
C-------------------------------------------------------------
C     Write out strength for CCS, CC2, and CCSD on unit LUOSC.
C-------------------------------------------------------------
C
      LUOSC = LURES
      IF (OSCSTR .AND. (CCS.OR.CC2.OR.CCSD)) THEN
C
         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Transition properti'
     *     //'es                    |'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        | Dipole Strength(a.u.) | Oscillator stre'
     *    //'ngth  | Direction   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
         DO 9001 ISYM  = 1, NSYM
          DO 9002 IEX   = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOSCSI = KOSCS2 + 3*3*(ISTATE-1)
           KTRSI  = KTRS   + 3*3*(ISTATE-1)
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_OSCPRI(WORK(KTRSI),WORK(KOSCSI),EIGV,
     *                    IEX,ISYM,WORK(KEND2),LEND2,MODELP,LCALC,
     *                    LDIP,LUOSC)
 9002     CONTINUE

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO 9003 ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
 9003        CONTINUE
             IF (NREST.EQ.0) GOTO 9001
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          ENDIF
 9001    CONTINUE
C
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      ENDIF
C
      LUOSC = LURES
      IF (VELSTR .AND. (CCS.OR.CC2.OR.CCSD)) THEN
C
         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Transition properti'
     *     //'es                    |'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        | Veloc. Strength(a.u.) | Oscillator stre'
     *    //'ngth  | Direction   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
         DO 9005 ISYM  = 1, NSYM
          DO 9006 IEX   = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOSCSI = KVELST + 3*3*(ISTATE-1)
           KTRSI  = KVELST2+ 3*3*(ISTATE-1)
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_OSCPRI(WORK(KTRSI),WORK(KOSCSI),EIGV,
     *                    IEX,ISYM,WORK(KEND2),LEND2,MODELP,LCALC,
     *                    LDIP,LUOSC)
 9006     CONTINUE

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO 9007 ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
 9007        CONTINUE
             IF (NREST.EQ.0) GOTO 9005
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          ENDIF
 9005    CONTINUE
C
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      ENDIF
C
      IF (MIXSTR .AND. (CCS.OR.CC2.OR.CCSD)) THEN
C
         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Mixed    Gauge Osci'
     *     //'llator Strength       |'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        | Dipole Strength(a.u.) | Oscillator stre'
     *    //'ngth  | Direction   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
         DO ISYM  = 1, NSYM
          DO IEX   = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOSCSI = KMIXST + 3*3*(ISTATE-1)
           KTRSI  = KMIXST2+ 3*3*(ISTATE-1)
           LCALC  = .FALSE.
           LDIP   = 3
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_OSCPRI(WORK(KTRSI),WORK(KOSCSI),EIGV,
     *                    IEX,ISYM,WORK(KEND2),LEND2,MODELP,LCALC,
     *                    LDIP,LUOSC)
          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9008
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9008     CONTINUE
         END DO
C
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      END IF
C
      LUOSC = LURES
      IF (ROTLEN .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rota'
     *     //'tory Strength         |'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Direction   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KTRSI  = KROTL + 3*(ISTATE-1)
           KSTREN = KCHKL + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_ROTPRI(WORK(KTRSI),WORK(KSTREN),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9009
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9009     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

      END IF
C
      LUOSC = LURES
      IF (ROTVEL .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rota'
     *     //'tory Strength         |'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Direction   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KTRSI  = KROTV + 3*(ISTATE-1)
           KSTREN = KCHKV + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_ROTPRI(WORK(KTRSI),WORK(KSTREN),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9010
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9010     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

      END IF

      LUOSC = LURES
      IF (RTNLEN .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rot.'
     *     //'Str. Tensor, El. Quad.|'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFQ  = KRQL  + 3*9*(ISTATE-1)
           KOFQ2  = KRQL2 + 3*3*(ISTATE-1)
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTQPRI(WORK(KOFFQ),WORK(KOFQ2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,NWRL)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9011
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9011     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rot.'
     *     //'Str. Tensor, Mag. Dip.|'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFM  = KRML  + 3*3*(ISTATE-1)
           KOFM2  = KRML2 + 3*3*(ISTATE-1)
           KSTREN = KCHKL + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTMPRI(WORK(KOFFM),WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRL)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9012
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9012     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         CALL DAXPY(LRML,1.0D0,WORK(KRQL2),1,WORK(KRML2),1)  ! Get total tensor (in KRML2)

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rot.'
     *     //'Str. Tensor, Total    |'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFM2  = KRML2 + 3*3*(ISTATE-1)
           KSTREN = KCHKL + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTTPRI(WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRL)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9013
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9013     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         IF (NWRL .NE. 0) THEN
            WRITE(LUOSC,'(//,1X,A,I4,A)')
     &      '***NOTICE:',NWRL,' warnings issued for Rot. Str. Tensors.'
            WRITE(LUOSC,'(1X,A)')
     &      '           Length gauge tensors are wrong!'
         END IF

      END IF

      LUOSC = LURES
      IF (RTNVEL .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rot.'
     *     //'Str. Tensor, El. Quad.|'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFQ  = KRQV  + 3*9*(ISTATE-1)
           KOFQ2  = KRQV2 + 3*3*(ISTATE-1)
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTQPRI(WORK(KOFFQ),WORK(KOFQ2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,NWRV)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9014
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9014     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rot.'
     *     //'Str. Tensor, Mag. Dip.|'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFM  = KRMV  + 3*3*(ISTATE-1)
           KOFM2  = KRMV2 + 3*3*(ISTATE-1)
           KSTREN = KCHKV + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTMPRI(WORK(KOFFM),WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRV)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9015
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9015     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         CALL DAXPY(LRMV,1.0D0,WORK(KRQV2),1,WORK(KRMV2),1)  ! Get total tensor (in KRMV2)

         WRITE(LUOSC,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rot.'
     *     //'Str. Tensor, Total    |'
         WRITE(LUOSC,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     *     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     *    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFM2  = KRMV2 + 3*3*(ISTATE-1)
           KSTREN = KCHKV + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTTPRI(WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRV)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9016
             WRITE(LUOSC,'(A)')
     *       ' +----------------------------------------------'
     *      //'-------------------------------+'
          END IF
 9016     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     *     ' +=============================================='
     *    //'===============================+'

         IF (NWRV .NE. 0) THEN
            WRITE(LUOSC,'(//,1X,A,I4,A)')
     &      '***NOTICE:',NWRV,' warnings issued for Rot. Str. Tensors.'
            WRITE(LUOSC,'(1X,A)')
     &      '           Velocity gauge tensors are wrong!'
         END IF

      END IF

      LUOSC = LURES
      IF (ROTLEN .OR. ROTVEL .OR. RTNLEN .OR. RTNVEL) THEN
         WRITE(LUOSC,'(/,1X,A)')
     &   'Conversion factors for rotatory strengths:'
         WRITE(LUOSC,'(3X,A,F15.10,A)')
     &   'SI  units:   1 a.u. = ',RAUSI,'D-55 A^2 m^3 s'
         WRITE(LUOSC,'(3X,A,F15.10,A)')
     &   'cgs units:   1 a.u. = ',RAUCGS,'D-40 cm^5 g s^-2'
      END IF
C
      TIMTOT = SECOND() - TIMTOT
      WRITE(LUPRI,'(/,1X,A,I5,A,F10.2,A)')
     & '>>> Time for',NTOT,' linear response residues: ',
     & TIMTOT,' seconds.'
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
c*DECK CC_LRSD
      SUBROUTINE CC_LRSD(LABELA,ISYMA,
     *                   LABELB,ISYMB,
     *                   ISTATE,RES1,RES2,WORK,LWORK)
C
C------------------------------------------------------------------------
C
C     Purpose: Calculate etaA*tB contribution to second order properties.
C
C
C     Written by Ove Christiansen 21-6-1996
C     New version november 1996.
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
#include <cclres.h>
#include <dummy.h>
C
      PARAMETER( TWO = 2.0D00,TOLFRQ=1.0D-08 )

      DIMENSION WORK(LWORK)
      CHARACTER LABELA*8,LABELB*8,MODEL*10
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_LRSD: Calculating residues   ')
      ENDIF
C
C------------------------
C     Allocate workspace.
C------------------------
C
      IF (ISYMA .NE. ISYMB ) CALL QUIT('Symmetry mismatch in CC_EATB')
      NTAMPB = NT1AM(ISYMB) + NT2AM(ISYMB)
      IF ( CCS ) NTAMPB = NT1AM(ISYMB)
      NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
      IF ( CCS ) NTAMPA = NT1AM(ISYMA)
C
      KETA  = 1
      KEND1 = KETA  + NTAMPA
      LEND1 = LWORK - KEND1
C
      KKSI  = KETA
C
      KR1   = KEND1
      KEND2 = KR1   + NTAMPB
      LEND2 = LWORK - KEND2
C
      IF (LEND2 .LT. 0)
     *      CALL QUIT('Insufficient space for allocation in CC_EATB')
C
C---------------------------------------------
C     Calculate first contribution to residue.
C---------------------------------------------
C
      CALL CC_ETAC(ISYMA,LABELA,WORK(KETA),'L0',1,0,
     *             DUMMY,WORK(KEND1),LEND1)
C
      KR11 = KR1
      KR12 = KR1 + NT1AM(ISYMB)
      IOPT   = 3
      CALL CC_RDRSP('RE',ISTATE,ISYMA,IOPT,MODEL,WORK(KR11),
     *              WORK(KR12))
C
      EATBCN = DDOT(NTAMPA,WORK(KETA),1,WORK(KR1),1)
C
      IF ( IPRINT .GT. 9 ) THEN
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
      IF (IPRINT .GT. 2 ) THEN
          WRITE(LUPRI,'(1X,A1,A8,A3,A,F10.6)')
     *    '<',LABELA,'|f>',' EtaA*RE cont. = ',EATBCN
      ENDIF
      RES1       = EATBCN  + RES1
C
C-------------------------------------
C     Calculate F-matrix contribution.
C-------------------------------------
C
      IF ((.NOT. CIS).AND.(.NOT.LRS2N1)) THEN
        IOPT   = 3
        KF11   = KETA
        KF12   = KETA + NT1AM(ISYMA)
        ILSTNR = IR1TAMP(LABELA,.FALSE.,-EIGVAL(ISTATE),ISYMA)
        CALL CC_RDRSP('F1',ILSTNR,ISYMB,IOPT,MODEL,WORK(KF11),
     *               WORK(KF12))
        IF (IPRINT .GT. 40 ) THEN
          CALL AROUND( 'In CC_LRSD:  F-transformed resp. vector ' )
          CALL CC_PRP(WORK(KF11),WORK(KF12),ISYMB,1,1)
        ENDIF
      ENDIF
      IF ((.NOT. CIS).AND.LRS2N1) THEN
        CALL CC_XKSI(WORK(KETA),LABELA,ISYMA,0,DUMMY,WORK(KEND1),LEND1)
        ILSTNR = ILRMAMP(ISTATE,EIGVAL(ISTATE),ISYMA)
        CALL CC_RDRSP('M1',ILSTNR,ISYMA,IOPT,MODEL,WORK(KR11),
     *               WORK(KR12))
      ENDIF
C
      EATBCN = DDOT(NTAMPA,WORK(KETA),1,WORK(KR1),1)
C
      IF ( IPRINT .GT. 9 ) THEN
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
      IF ((IPRINT.GT.2).AND.(.NOT. CIS)) THEN
        IF (.NOT.LRS2N1) THEN
          WRITE(LUPRI,'(1X,A1,A8,A3,A,F10.6)')
     *    '<',LABELA,'|f>',' F*taA*RE cont. = ',EATBCN
        ELSE
          WRITE(LUPRI,'(1X,A1,A8,A3,A,F10.6)')
     *    '<',LABELA,'|f>',' Mf*KsiA  cont. = ',EATBCN
        ENDIF
      ENDIF
C
      IF (.NOT.CIS) RES1       = EATBCN  + RES1
C
C---------------------------------------------
C     Calculate second contribution to residue.
C---------------------------------------------
C
      CALL CC_XKSI(WORK(KETA),LABELB,ISYMB,0,DUMMY,WORK(KEND1),LEND1)
C
      KR11   = KR1
      KR12   = KR1 + NT1AM(ISYMB)

      CALL CC_RDRSP('LE',ISTATE,ISYMB,IOPT,MODEL,WORK(KR11),
     *              WORK(KR12))
      IF (IPRINT .GT. 40 ) THEN
         CALL AROUND( 'In CC_LRSD:  Left Eigen vector ' )
         CALL CC_PRP(WORK(KR1),WORK(KR1+NT1AM(ISYMB)),ISYMB,1,1)
      ENDIF
C
      EATBCN = DDOT(NTAMPA,WORK(KETA),1,WORK(KR1),1)
C
      IF ( IPRINT .GT. 9 ) THEN
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
      IF (IPRINT .GT. 2 ) THEN
          WRITE(LUPRI,'(1X,A3,A8,A1,A,F10.6)')
     *    '<f|',LABELB,'>',' LE*XksiB cont. = ',EATBCN
      ENDIF
      RES2       = EATBCN  + RES2
C
      RETURN
      END
c*DECK CC_OSCPRI
      SUBROUTINE CC_OSCPRI(TRS,OSC,EIGV,IEX,ISYM,WORK,LWORK,MODEL,LCALC,
     *                     LDIP,LUOSC)
C
C------------------------------------------------------------------------
C
C     Purpose: Calculate LD*ksiC contribution to second order properties.
C              For use in calculation of molecular properties from
C              Asymmetric formulaes not in accordance with 2n+2 rule for
C              the multipliers, left vector, t-bar, lamdas, zeta or
C              whatever your preferred choice is today.
C
C     Written by Ove Christiansen 17-10-1996/7-11-1996
C
C------------------------------------------------------------------------
C
#include <implicit.h>
#include <pgroup.h>
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
C
      DIMENSION OSC(*),PVAL(3),PAXIS(3,3)
      CHARACTER MODEL*10,CDIP*7
      LOGICAL LCALC
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_OSCPRI: Output transition properties ' )
      ENDIF
C
C------------------------------------------
C     write out transition strength matrix.
C------------------------------------------
C

      IMULT = 1  ! force singlet spin symmetry...

      IF (LCALC) THEN
C
      WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     *    'Transition strength matrix for state nr.',IEX,
     *     ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
      IF (LDIP .EQ. 1) THEN
        WRITE(LUPRI,'(1X,A)') 'Gauge: length'
      ELSE IF (LDIP .EQ. 2) THEN
        WRITE(LUPRI,'(1X,A)') 'Gauge: velocity'
      ELSE IF (LDIP .EQ. 3) THEN
        WRITE(LUPRI,'(1X,A)') 'Gauge: mixed length/velocity'
      ELSE
        WRITE(LUPRI,'(1X,A)') 'Gauge: UNKNOWN'
        WRITE(LUPRI,'(1X,A)') '- scaling factors will be incorrect!'
      ENDIF
      CALL OUTPUT(TRS,1,3,1,3,3,3,1,LUPRI)
C
      CALL TNSRAN(TRS,PVAL,PAXIS,
     *            ALFSQ,BETSQ,ITST,ITST2,
     *            APAR1,APEN1,XKAPPA,IPAR)
      WRITE(LUPRI,'(/,1X,A,/)')
     *    'Principal values of diagonalized transition strength matrix:'
      WRITE(LUPRI,'(1X,A)') '            a.u.               '
      WRITE(LUPRI,'(1X,A,F16.8)') '1 ',PVAL(1)
      WRITE(LUPRI,'(1X,A,F16.8)') '2 ',PVAL(2)
      WRITE(LUPRI,'(1X,A,F16.8)') '3 ',PVAL(3)
      WRITE(LUPRI,'(/,1X,A,/)')
     *    'Principal axis of diagonalized transition strength matrix:'
      CALL OUTPUT(PAXIS,1,3,1,3,3,3,1,LUPRI)
      TRA = PVAL(1)+PVAL(2)+PVAL(3)
C
C------------------------------------------
C     First scale it - then
C     write out oscillator strength matrix.
C------------------------------------------
C
      IF (LDIP .EQ. 1) THEN
         FACT = EIGV*2.0D0/3.0D0
      ELSE IF (LDIP .EQ. 2) THEN
         FACT = -2.0D0/(3.0D0*EIGV)
      ELSE IF (LDIP .EQ. 3) THEN
         FACT = -2.0D0/3.0D0
      ELSE
         FACT = 1.0D0
      ENDIF
      CALL DSCAL(3*3,FACT,OSC,1)
      WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     *    ' oscillator strength matrix for state nr.',IEX,
     *    ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
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

      CALL WRIPRO(OSCS,MODEL,400,
     &            "OSCI-LEN","OSCI-LEN","OSCI-LEN","OSCI-LEN",
     &            EIGV,EIGV,EIGV,ISYM,ISYM,1,IEX)

      CDIP = 'unknown'
      IF (IPAR .EQ.1) CDIP = '   X   '
      IF (IPAR .EQ.2) CDIP = '   Y   '
      IF (IPAR .EQ.3) CDIP = '   Z   '
      IF (IPAR .EQ.4) CDIP = ' (X,Y) '
      IF (IPAR .EQ.5) CDIP = ' (X,Z) '
      IF (IPAR .EQ.6) CDIP = ' (Y,Z) '
      IF (IPAR .EQ.7) CDIP = '(X,Y,Z)'
      IF (IPAR .EQ.8) CDIP = '   -   '
c
c     IF ( IEX .EQ. 1) THEN
C IMULT = 1 is hardwired in since for linear response residues
C only singlet states have a non-vanishing oscillator strength anyway
         WRITE(LUOSC,9988) IMULT,REP(ISYM-1),IEX,TRA,OSCS,CDIP
c     ELSE
c        WRITE(LUOSC,9989) IEX,TRA,OSCS,CDIP
c     ENDIF
C
      ELSE IF (.NOT.LCALC) THEN
         CDIP = '   ?   '
c        IF ( IEX .EQ. 1) THEN
           WRITE(LUOSC,9986) IMULT,REP(ISYM-1),IEX,'Not calculated',
     *                       'Not calculated',CDIP
c        ELSE
c          WRITE(LUOSC,9987) IEX,'Not calculated','Not calculated',CDIP
c        ENDIF
      ENDIF
C
 9986 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,'  |')
 9987 FORMAT(1X,'|       | ',I4,'   | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,'  |')
 9988 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9989 FORMAT(1X,'|       | ',I4,'   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
C
      END
      SUBROUTINE CC_TSTAV(ILSTNR,VEC,WORK,LWORK,IOPTTST)
C
C----------------------------------------------------------------------
C
C     Purpose: Calculate first order property from first order response
C              vectors to test these.
C              NOT MEANT to advocate this way of calculating
C              expectation values.
C
C     Written by Ove Christiansen 10-5-1996 / 2.0: 13-3-1997
C
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <iratdef.h>
C
      LOGICAL LOCDBG
      PARAMETER( LOCDBG = .FALSE. )
      PARAMETER( TWO = 2.0D00,XHALF = 0.5D00 )
      DIMENSION WORK(LWORK),VEC(*)
      CHARACTER*10 MODEL
C
#include <ccorb.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccroper.h>
#include <ccr1rsp.h>
#include <ccx1rsp.h>
#include <leinf.h>
C
C-------------------------------------------------------------
C        Calculate response contribution to expectation value.
C-------------------------------------------------------------
C
      IF (.NOT.CCS) THEN
         NVAR   = NT1AM(ISYMOP) + NT2AM(ISYMOP)
         IF (CCR12) THEN
           NVAR = NVAR + NTR12AM(ISYMOP)
         ENDIF
         KETA   = 1
         KWRK1  = KETA  + NVAR
         LWRK1  = LWORK - KWRK1
         IF (LWRK1 .LT. 0 )
     &        CALL QUIT('Too little workspace in CC_TSTAV-1')
         IF      (IOPTTST.EQ.0) THEN
            CALL CC_ETA(WORK(KETA),WORK(KWRK1),LWRK1)
         ELSE IF (IOPTTST.EQ.1) THEN
            IOPT  = 3
            CALL CC_RDRSP('L0 ',0,ISYMOP,IOPT,MODEL,WORK(KETA),
     *                    WORK(KETA+NT1AM(ISYMOP)))
            IF (CCR12) THEN
              IOPT = 32
              CALL CC_RDRSP('L0 ',0,ISYMOP,IOPT,MODEL,DUMMY,
     *                      WORK(KETA+NT1AM(ISYMOP)+NT2AM(ISYMOP)))
            ENDIF  
         ELSE
            WRITE(LUPRI,*) 'IOPTTST = ',IOPTTST
            CALL QUIT('ILLEGAL VALUE FOR IOPTTST IN CC_TSTAV.')
         END IF
         PROPRSP = DDOT(NVAR,WORK(KETA),1,VEC,1)

         IF (LOCDBG) THEN
           write(lupri,*) 'Input vector:'
           call cc_prp(vec,vec(nt1am(isymop)+1),isymop,1,1)
           if (CCR12) call cc_prpr12(vec(1+nt1am(isymop)+nt2am(isymop)),
     *                               isymop,1,.true.)
           write(lupri,*) 'L0/X0 vector:'
           call cc_prp(work(keta),work(keta+nt1am(isymop)),isymop,1,1)
           if (CCR12) call cc_prpr12(work(keta+nt1am(isymop)+
     *                               nt2am(isymop)),isymop,1,.true.)
           write(lupri,*) 'PROPRSP:',PROPRSP
         END IF
      ELSE
         PROPRSP = 0.0D0
         KWRK1   = 1
         LWRK1   = LWORK
      ENDIF
C
C------------------------------------------
C     Calculate average value contribution.
C------------------------------------------
C
      ! find operator index
      ISYM  = 1
      IOPER = IROPER(LRTLBL(ILSTNR),ISYM)
C
      IF ( LORXLRT(ILSTNR) .OR. LPDBSOP(IOPER) ) THEN
        ! if the orbitals are allowed to relax in the field or if the
        ! basis set depends on the perturbation, read the average
        ! value contribution from the ccx1rsp.h common blocks
        ILSTETA = IETA1(LRTLBL(ILSTNR),LORXLRT(ILSTNR),
     &                  FRQLRT(ILSTNR),ISYM)
        PROPAVE = AVEX1(ILSTETA)
      ELSE
        ! if it is a simple unrelaxed one-electron perturbation
        ! calculate the average value contribution in CC_AVE
        FF = 1.0D00
        CALL CC_AVE(PROPAVE,LRTLBL(ILSTNR),WORK(KWRK1),LWRK1,FF)
      END IF
C
      WRITE(LUPRI,'(1X,A,A)') 'Operator: :   ',LRTLBL(ILSTNR)
      WRITE(LUPRI,'(1X,A,F16.10)') 'Average contribution:   ',
     *                         PROPAVE
      WRITE(LUPRI,'(1X,A,F16.10)') 'Response contribution:  ',
     *                         PROPRSP
      WRITE(LUPRI,'(1X,A,F16.10)') 'Total expectation value:',
     *                         PROPAVE + PROPRSP
C
      CALL FLSHFO(LUPRI)
      END
      SUBROUTINE CC_AVE(XVALUE,LBL,WORK,LWORK,FF)
C
C-----------------------------------------------------------------------
C
C     Purpose: Calculate <HF|C|CC> contribution to first order property.
C              C is assumed to be a one-electron operator.
C
C     Written by Ove Christiansen 10-5-1996
C
C     Bug-Fix for frozen-core calculations: Chr. Neiss  22-04-2005
C
C-----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <iratdef.h>
C
      CHARACTER LBL*(*), MODEL*10
      DIMENSION WORK(LWORK)
      INTEGER ICMO(8,8), NCMO(8), IGLMRHS(8,8), NGLMDS(8), IGLMVIS(8,8)
C
#include <ccorb.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <leinf.h>
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_AVE ')
      ENDIF
C
      DO ISYM = 1, NSYM
         ICOUNT  = 0
         ICOUNT2 = 0
         ICOUNT3 = 0
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            ICMO(ISYM1,ISYM2)    = ICOUNT
            IGLMRHS(ISYM1,ISYM2) = ICOUNT2
            ICOUNT  = ICOUNT  + NBAS(ISYM1)*NORBS(ISYM2)
            ICOUNT2 = ICOUNT2 + NBAS(ISYM1)*NRHFS(ISYM2)
            ICOUNT3 = ICOUNT3 + NBAS(ISYM1)*NRHFS(ISYM2)
         END DO
         NCMO(ISYM)   = ICOUNT
         NGLMDS(ISYM) = ICOUNT2
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            IGLMVIS(ISYM1,ISYM2) = ICOUNT3
            ICOUNT3 = ICOUNT3 + NBAS(ISYM1)*NVIRS(ISYM2)
         END DO
      END DO
C
      KONEP  = 1
      KT1AM  = KONEP  + N2BST(ISYMOP)
      KLAMDPS= KT1AM  + NT1AMX
      KLAMDHS= KLAMDPS+ NGLMDS(1)
      KEND1  = KLAMDHS+ NGLMDS(1)
      LWRK1  = LWORK  - KEND1
      IF ( LWRK1 .LT. 0 )
     *     CALL QUIT(' Too little workspace in CC_AVE')
C
      CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
      CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,1,LBL)
C
      IF (.NOT.CCS) THEN
         IOPT = 1
         CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      ELSE IF (CCS ) THEN
         CALL DZERO(WORK(KT1AM),NT1AMX)
      ENDIF
C
      CALL LAMMATS(WORK(KLAMDPS),WORK(KLAMDHS),WORK(KT1AM),ISYMOP,
     *             .FALSE.,.FALSE.,NGLMDS,IGLMRHS,IGLMVIS,ICMO,
     *             WORK(KEND1),LWRK1)
C
      XVALUE = 0.0D0
C 
      DO ISYM = 1, NSYM

        KSCR1 = KEND1
        KEND2 = KSCR1 + NBAS(ISYM) * NRHFS(ISYM)
        LWRK2 = LWORK  - KEND2
        IF ( LWRK2 .LT. 0 ) THEN
          WRITE (LUPRI,*) 'LWORK, LWRK2: ',WORK, LWRK2
          CALL QUIT('Too little workspace in CC_AVE')
        END IF

        NBAS1 = MAX(NBAS(ISYM),1)
        KOFF1 = KONEP   + IAODIS(ISYM,ISYM)
        KOFF2 = KLAMDHS + IGLMRHS(ISYM,ISYM)

        CALL DGEMM('N','N',NBAS(ISYM),NRHFS(ISYM),NBAS(ISYM),
     *             1.0D0,WORK(KOFF1),NBAS1,WORK(KOFF2),NBAS1,
     *             0.0D0,WORK(KSCR1),NBAS1)

        KOFF2 = KLAMDPS + IGLMRHS(ISYM,ISYM)

        TRACE = DDOT(NBAS(ISYM)*NRHFS(ISYM),
     &                 WORK(KOFF2),1,WORK(KSCR1),1)
        XVALUE = XVALUE + 2.0D0 * TRACE
      END DO
C
      END
c*DECK CC_XKSI
      SUBROUTINE CC_XKSI(XKSI,LBPERT,ISYMPT,IOPTCC2,XINT,WORK,LWORK)
C
C----------------------------------------------------------------------
C
C     Purpose: Calculate XKSI vector.
C
C     IOPTCC2 = 1 -- use for CC2 the CMO vector instead of the lambda
C                    matrices to transform the Fock mat. in the E-term
C
C     SLV98,OC: Allow for input of integrals if
C               LBPERT.eq.'GIVE INT'
C
C     Written by Ove Christiansen 16 februar 1996
C
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <iratdef.h>
C
      PARAMETER( TWO = 2.0D00,XHALF = 0.5D00 )
      LOGICAL FCKCON,ETRAN
      INTEGER IOPTCC2
      CHARACTER LBPERT*(*), MODEL*10
      DIMENSION XKSI(*),WORK(LWORK),XINT(*)
C
#include <ccorb.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <leinf.h>
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_XKSI: Constructing XKSI vector ')
      ENDIF
C
C-------------------------------------------------------------------
C     Read in AO property integrals and transform to T1 transformed
C     MO basis.
C-------------------------------------------------------------------
C
      KFOCK  = 1
      KT1AM  = KFOCK  + N2BST(ISYMPT)
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        KCMO   = KEND1
        KFCKHF = KCMO    + NLAMDT
        KEND1  = KFCKHF  + N2BST(ISYMPT)
      END IF
C
      LEND1  = LWORK  - KEND1
C
      IF ( .NOT. CCS) THEN
C
         KT2AM = KEND1
         KEND2 = KT2AM + NT2SQ(1)
         LEND2 = LWORK - KEND2
C
         KT2PK = KEND2
         KEND3 = KT2PK + NT2AMX
         LEND3 = LWORK - KEND3
C
      ELSE
C
         KEND2 = KEND1
         LEND2 = LEND1
         KEND3 = KEND1
         LEND3 = LEND1
C
      ENDIF
C
      IF (LEND3 .LT. 0 ) THEN
         WRITE(LUPRI,*) 'Requested workspace, available workspace =',
     *               KEND3,LWORK
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_XKSI-1')
      ENDIF
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
C
      IF (.NOT.(CCS.OR.CCP2)) THEN
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2PK))
         CALL CC_T2SQ(WORK(KT2PK),WORK(KT2AM),1)
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND2),LEND2)
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        REWIND(LUSIFC)
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ(LUSIFC)
        READ(LUSIFC)
        READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
        CALL GPCLOSE(LUSIFC,'KEEP')
        CALL CMO_REORDER(WORK(KCMO),WORK(KEND2),LEND2)
      END IF
C
      CALL DZERO(WORK(KFOCK),N2BST(ISYMPT))
C
C SLV98,OC if option for solvent
C
      IF (LBPERT.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMPT),XINT,1,WORK(KFOCK),1)
      ELSE
        FF = 1.0D0
        CALL CC_ONEP(WORK(KFOCK),WORK(KEND2),LEND2,FF,ISYMPT,LBPERT)
      ENDIF
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        CALL DCOPY(N2BST(ISYMPT),WORK(KFOCK),1,WORK(KFCKHF),1)
        CALL CC_FCKMO(WORK(KFCKHF),WORK(KCMO),WORK(KCMO),
     *                WORK(KEND2),LEND2,ISYMPT,1,1)
      END IF
C
      CALL CC_FCKMO(WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND2),LEND2,ISYMPT,1,1)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC_XKSI: MO^(t1) property matrix' )
         CALL CC_PRFCKMO(WORK(KFOCK),ISYMPT)
      ENDIF
C
C------------------------------
C     Contract into ksi vector.
C     first zero result.
C------------------------------
C
      CALL DZERO(XKSI(1),NT1AM(ISYMPT))
      IF (.NOT. CCS) CALL DZERO(XKSI(1+NT1AM(ISYMPT)),NT2AM(ISYMPT))
C
C----------------------
C     Calculate J-term.
C----------------------
C
      CALL CCRHS_J(XKSI(1),ISYMPT,WORK(KFOCK))
C
      IF (.NOT. CCS) THEN
C
C----------------------------------
C        Calculate E contributions.
C----------------------------------
C
         KEI1  = KEND2
         KEI2  = KEI1 + NEMAT1(ISYMPT)
         KEND3 = KEI2 + NMATIJ(ISYMPT)
         LEND3 = LWORK - KEND3
C
         IF (LEND3.LT. 0 )
     &        CALL QUIT(' TOO LITTLE WORKSPACE IN CC_XKSI-2')
C
         FCKCON = .TRUE.
         ETRAN  = .FALSE.
C
         IF (CC2 .AND. IOPTCC2.EQ.1) THEN
           CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KCMO),
     *                     WORK(KFCKHF),WORK(KEND3),LEND3,FCKCON,
     *                     ETRAN,ISYMPT)
         ELSE
           CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                     WORK(KFOCK),WORK(KEND3),LEND3,FCKCON,
     *                     ETRAN,ISYMPT)
         END IF
C
         CALL CCRHS_E(XKSI(1+NT1AM(ISYMPT)),WORK(KT2AM),WORK(KEI1),
     *                WORK(KEI2),WORK(KEND3),LEND3,ISYMOP,ISYMPT)
C
         CALL CCLR_DIASCL(XKSI(1+NT1AM(ISYMPT)),XHALF,ISYMPT)
C
C-------------------------
C        Calculate I-term.
C-------------------------
C
         CALL CCRHS_T2TR(WORK(KT2AM),WORK(KEND2),LEND2,1)
C
         CALL CCRHS_I(XKSI(1),WORK(KT2AM),WORK(KFOCK),
     *             WORK(KEND2),LEND2,ISYMOP,ISYMPT)
C
      ENDIF
C
      IF (IPRINT .GT. 40 ) THEN
         NC2 = 1
         IF ( CCS ) NC2 = 0
         CALL AROUND( 'In CC_XKSI:  XKSI vector ' )
         CALL CC_PRP(XKSI(1),XKSI(1+NT1AM(ISYMPT)),ISYMPT,1,NC2)
      ENDIF
C
      IF ( IPRINT .GT. 10 ) THEN
         XKSI1 = DDOT(NT1AM(ISYMPT),XKSI(1),1,XKSI(1),1)
         WRITE(LUPRI,*) 'Norm of XKSI1: ',XKSI1
         IF ( .NOT. CCS ) THEN
            XKSI2 = DDOT(NT2AM(ISYMPT),XKSI(1+NT1AM(ISYMPT)),
     *               1,XKSI(1+NT1AM(ISYMPT)),1)
            WRITE(LUPRI,*) 'Norm of XKSI2: ',XKSI2
         ENDIF
         CALL AROUND( 'END OF CC_XKSI ')
      ENDIF
C
      END
c*DECK CC_ETAC
      SUBROUTINE CC_ETAC(ISYMC,LBLC,ETAC,LIST,ILSTNR,IOPTCC2,
     *                   XINT,WORK,LWORK)
C
C-----------------------------------------------------------------------
C
C     Purpose: Calculate ETAC vector.
C
C     Important note: Requires work space of dimension of
C             NT2AM + NT2SQ in addition to ETAC, so please take care.
C
C     eta(tau,nu)= (<HF| + Sum(mu)L(0 or 1)<mu|)
C                         exp(-t)[C,tau,nu]exp(T)|HF>
C
C     LIST= 'L0' for zeroth order left amplitudes.
C                ISYML should be ISYMOP in this case.
C           'L1' for first order left amplitudes, read in from file
C                In this case the vector is found according to its list
C                number ILSTNR.
C
C                For L1 HF contribution is skipped.
C
C     IOPTCC2 = 1 -- transform for CC2 the Fock matrix entering the
C                    E term contribution with CMO vector instead with
C                    Lambda matrices
C
C     C property integrals read according to LBLC
C
C     SLV98,OC: Allow for input of integrals if
C               LBLC.eq.'GIVE INT'
C
C
C     Written by Ove Christiansen 20-6-1996/1-11-1996
C
C
C-----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <cclr.h>
#include <ccexci.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
C
      PARAMETER( TWO = 2.0D00, XHALF = 0.5D00 )
      DIMENSION ETAC(*),WORK(LWORK),XINT(*)
      CHARACTER LBLC*(*),LIST*(*),MODEL*10
      INTEGER IOPTCC2
      LOGICAL  FCKCON,ETRAN
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_ETAC: Constructing EtaC vector ')
      ENDIF
C
C--------------------------------
C     find symmetry of D operator.
C--------------------------------
C
      ISYML = ILSTSYM(LIST,ILSTNR)
C
      ISYRES = MULD2H(ISYML,ISYMC)
      IF (( LIST .EQ. 'L0').AND.(ISYML.NE.1)) THEN
         CALL QUIT('Misuse of CC_ETAC')
      ENDIF
C
      TIMEC = SECOND()
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
C--------------------
C     Allocate space.
C--------------------
C
      KCTMO  = 1
      KT1AM  = KCTMO  + N2BST(ISYMC)
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        KCMO   = KEND1
        KFCKHF = KCMO   + NLAMDT
        KEND1  = KFCKHF + N2BST(ISYMC)
      END IF
C
      LEND1  = LWORK  - KEND1
C
      IF ( .NOT. CCS) THEN
C
         KL1AM = KEND1
         KL2AM = KL1AM + NT1AM(ISYML)
         KEND2 = KL2AM + NT2SQ(ISYML)
         LEND2 = LWORK - KEND2
         KT2AM = KEND2
         KEND21= KT2AM + MAX(NT2AM(ISYML),NT2AM(1))
         LEND21= LWORK - KEND2
C
      ELSE
C
         KL1AM = KEND1
         KEND2 = KL1AM + NT1AM(ISYML)
         LEND2 = LEND1
         KEND21= KEND1
         LEND21= LEND1
C
      ENDIF
C
      IF (LEND21.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-1')
C
C-----------------------
C     get T1 amplitudes.
C-----------------------
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      IF ( .NOT. CCS) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND21),LEND21)
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        REWIND(LUSIFC)
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ(LUSIFC)
        READ(LUSIFC)
        READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
        CALL GPCLOSE(LUSIFC,'KEEP')
        CALL CMO_REORDER(WORK(KCMO),WORK(KEND21),LEND21)
      END IF
C
C-------------------------------
C     get AO property integrals.
C-------------------------------
C
      CALL DZERO(WORK(KCTMO),N2BST(ISYMC))
      FF = 1.0D0
C SLV98,OC give integrals option
      IF (LBLC.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMC),XINT(1),1,WORK(KCTMO),1)
      ELSE
        FF = 1.0D0
        CALL CC_ONEP(WORK(KCTMO),WORK(KEND21),LEND21,FF,ISYMC,LBLC)
      ENDIF
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        CALL DCOPY(N2BST(ISYMC),WORK(KCTMO),1,WORK(KFCKHF),1)
        CALL CC_FCKMO(WORK(KFCKHF),WORK(KCMO),WORK(KCMO),
     *                WORK(KEND21),LEND21,ISYMC,1,1)
      END IF
C
C-----------------------------------------------
C     Make MO T1-transformed property integrals.
C-----------------------------------------------
C
      CALL CC_FCKMO(WORK(KCTMO),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND21),LEND21,ISYMC,1,1)
C
C----------------------------------------------
C     Calculate 2Cai Hartree-Fock contribution.
C----------------------------------------------
C
      CALL DZERO(ETAC,NT1AM(ISYRES))
C
      IF (LIST .EQ. 'L0') THEN
         DO 100 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMC)
C
            DO 110 A = 1,NVIR(ISYMA)
C
               DO 120 I = 1,NRHF(ISYMI)
C
                  KOFF1 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                  KOFF2 = KCTMO + IFCVIR(ISYMI,ISYMA)
     *                  + NORB(ISYMI)*(A - 1) + I - 1
C
                  ETAC(KOFF1) = TWO*WORK(KOFF2)
C
  120          CONTINUE
  110       CONTINUE
C
  100    CONTINUE
C
      ENDIF
C
      IF ( DEBUG ) THEN
         ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,1) 'Norm of ETAC - First contribution:',ETA1
      ENDIF
C
C------------------------
C     IF CCS then return.
C------------------------
C
      IF ( CCS .AND. (LIST .EQ. 'L0')) RETURN
C
C----------------------------------------------
C     Read zero'th order amplitude multipliers.
C----------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP(LIST,ILSTNR,ISYML,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KT2AM))
      IF (.NOT. CCS) CALL CC_T2SQ(WORK(KT2AM),WORK(KL2AM),ISYML)
C
C--------------------------------
C     Put T2 amplitudes in etac2.
C--------------------------------
C
      IF (.NOT. CCS) THEN
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
      ENDIF
C
C--------------------------------
C     Make X and Y intermediates.
C--------------------------------
C
      IF (.NOT. CCS) THEN
         KXMAT = KEND21
         KYMAT = KXMAT + NMATIJ(ISYML)
         KEND3 = KYMAT + NMATAB(ISYML)
         LEND3 = LWORK - KEND3
         IF (LEND3.LT. 0 )
     &        CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-2')
C
         IF ( DEBUG ) THEN
            XYI   = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L1AM vector:              ',XYI
            XYI   = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L2AM vector:              ',XYI
            XXI   = DDOT(NT2AM(ISYMOP),WORK(KT2AM),1,WORK(KT2AM),1)
            WRITE(LUPRI,1) 'T2AM vector :                      ',XXI
         ENDIF
         CALL CC_XI(WORK(KXMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND3),LEND3)
         CALL CC_YI(WORK(KYMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND3),LEND3)
         IF ( DEBUG ) THEN
            XYI   = DDOT(NMATAB(ISYML),WORK(KYMAT),1,WORK(KYMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: YI  intermediate is:      ',XYI
            XXI   = DDOT(NMATIJ(ISYML),WORK(KXMAT),1,WORK(KXMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: XI  intermediate is:      ',XXI
         ENDIF
      ELSE
         KEND3 = KEND2
         LEND3 = LEND2
      ENDIF
C
C----------------------------------------------
C     Calculate X and Y contributions to etac1.
C     etac1 = -sum(e)Cie*Yae - sum(l)Cla*Xli
C----------------------------------------------
C
      IF ( (.NOT.CCS) .AND. (.NOT.(CC2.AND.IOPTCC2.EQ.1)) ) THEN
         CALL CC_21EFM(ETAC,WORK(KCTMO),ISYMC,WORK(KXMAT),
     *                 WORK(KYMAT),ISYML,WORK(KEND3),LEND3)
C
         IF ( DEBUG ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            WRITE(LUPRI,1) 'Norm of eta1-after X&Y cont:       ',ETA1
         ENDIF
      ENDIF
C
C------------------------------------------------
C     Workspace for T2AM and X and Y is now free.
C     etac2 = P(ab,ij)(2l(ai)*Cjb - l(aj)*c(ib))
C------------------------------------------------
C
      IF (.NOT. CCS) THEN
         CALL DZERO(ETAC(1+NT1AM(ISYRES)),NT2AM(ISYRES))
         CALL CC_L1FCK(ETAC(1+NT1AM(ISYRES)),WORK(KL1AM),WORK(KCTMO),
     *                 ISYML,ISYMC,WORK(KEND2),LEND2)
C
         IF ( DEBUG ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            ETA2 = DDOT(NT2AM(ISYRES),ETAC(1+NT1AM(ISYRES)),1,
     *                  ETAC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1-after L1c cont:       ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2-after L1c cont:       ',ETA2
         ENDIF
      ENDIF
C
      KEI1   = KEND2
      KEI2   = KEI1   + NEMAT1(ISYMC)
      KEND3  = KEI2   + NMATIJ(ISYMC)
      LEND3  = LWORK  - KEND3
      IF (LEND3.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-3')
C
C--------------------------------
C     Put A into E matrix format.
C--------------------------------
C
      FCKCON = .TRUE.
      ETRAN  = .FALSE.
      CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                WORK(KCTMO),WORK(KEND3),LEND3,FCKCON,
     *                ETRAN,ISYMC)
C
C--------------------------------------------
C     etac1 =  sum(b)Lbi*Cba - sum(j)Laj*Cij.
C--------------------------------------------
C
      IF ( DEBUG ) THEN
         XE1 = DDOT(NMATAB(ISYMC),WORK(KEI1),1,WORK(KEI1),1)
         XE2 = DDOT(NMATIJ(ISYMC),WORK(KEI2),1,WORK(KEI2),1)
         WRITE(LUPRI,1) 'Norm of EI1  -after EFCK:          ',XE1
         WRITE(LUPRI,1) 'Norm of EI2  -after EFCK:          ',XE2
         ETA1 = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
         WRITE(LUPRI,1) 'Norm of L1AM before  CCLR_E1C1:    ',ETA1
      ENDIF
C
c test
c     kei11= kend3
c     kei21= kei11+ NMATAB(ISYMC)
c     kend3 = kei21+ NMATIJ(ISYMC)
c     lend3 = lwork -kend3
c     call dzero(work(kei11),NMATAB(ISYMC))
c     call dzero(work(kei21),NMATIJ(ISYMC))
c     call dcopy(NMATAB(ISYMC),work(kei1),1,work(kei11),1)
c     call dcopy(NMATIJ(ISYMC),work(kei2),1,work(kei21),1)
c     CALL CCLR_E1C1(ETAC,WORK(KL1AM),WORK(KEI11),WORK(KEI21),
c    *               WORK(KEND3),LEND3,ISYML,ISYMC,'T')
c test
C
      CALL CCLR_E1C1(ETAC,WORK(KL1AM),WORK(KEI1),WORK(KEI2),
     *               WORK(KEND3),LEND3,ISYML,ISYMC,'T')
C
      IF (DEBUG ) THEN
         ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
         WRITE(LUPRI,1) 'Norm of eta1 - after CCLR_E1C1:    ',ETA1
      ENDIF
C
C---------------------------------------------------------------
C     etac2 = P(ab,ij)(sum(e)2L(aiej)*Ceb - sym(k)L(aibk)*c(jk))
C---------------------------------------------------------------
C
      IF (.NOT. CCS) THEN
C
         IF (CC2 .AND. IOPTCC2.EQ.1) THEN
           FCKCON = .TRUE.
           ETRAN  = .FALSE.
           CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KCMO),
     *             WORK(KFCKHF),WORK(KEND3),LEND3,FCKCON,ETRAN,ISYMC)
         END IF

         CALL CC_EITR(WORK(KEI1),WORK(KEI2),WORK(KEND3),LEND3,
     *                ISYMC)
C
         CALL CCRHS_E(ETAC(1+NT1AM(ISYRES)),WORK(KL2AM),
     *                WORK(KEI1),WORK(KEI2),WORK(KEND3),
     *                LEND3,ISYML,ISYMC)
C
         IF (IPRINT .GT. 40 ) THEN
            CALL AROUND( 'In CC_ETAC:  EtaC vector ' )
            CALL CC_PRP(ETAC(1),ETAC(1+NT1AM(ISYRES)),ISYMC,1,1)
         ENDIF
C
         IF (DEBUG .OR. ( IPRINT .GT. 20 )) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            ETA2 = DDOT(NT2AM(ISYRES),ETAC(1+NT1AM(ISYRES)),1,
     *                  ETAC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1 - end of CC_ETAC:     ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2 - end of CC_ETAC:     ',ETA2
            CALL AROUND( 'END OF CC_ETAC ')
         ENDIF
      ENDIF
C
      IF (IPRINT .GT. 5 ) THEN
         TIMEC = SECOND() - TIMEC
         WRITE(LUPRI,9999) 'CC_ETA          ', TIMEC
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
9999  FORMAT(1x,'Time used in',2x,A18,2x,': ',f10.2,' seconds')
C
      END
c /* Deck polsym */
      SUBROUTINE POLSYM(A,FACT)
C
#include <implicit.h>
C
      DIMENSION A(3,3)
C
      DO 10 I = 1, 3
        DO 20 J = 1, I -1
           A(J,I) = (A(J,I) + A(I,J))*FACT
           A(I,J) = A(J,I)
  20    CONTINUE
        A(I,I) = 2.0D00*A(I,I)*FACT
  10  CONTINUE
C
      RETURN
      END
c*DECK TNSRAN
      SUBROUTINE TNSRAN(TNSR,PVAL,PAXIS,ALFSQ,BETSQ,ITST,ITST2,
     *                  APAR,APEN,XKAPPA,IPAR)
C
C------------------------------------------------------------------------
C
C     Purpose: Analyse 3 by 3 tensot and
C
C              1. calculate rotatinal invariants
C                 alfa**2 = ((TNSRxx+TNSRyy+TNSRzz)**2)/9
C                 beta**2 = [(TNSRxx-TNSRyy)**2 +
C                           (TNSRxx-TNSRzz)**2 +
C                           (TNSRyy-TNSRzz)**2 +
C                           +3(TNSRxy**2+TNSRxy**2+TNSRxy**2+
C                            TNSRxy**2+TNSRxy**2+TNSRxy**2)]/2
C
C              2. Diagonal, block-diagonal, all elements differ.
C                 itst = 0,   2,                 6
C                 itst = nr. of non-zero out of diagonal elements.
C              3. If diagonal then a. no symmetry.       itst2 = 3
C                                  b. cylinder symmetry. itst2 = 1
C                                  c. Spherical symmetry.itst2 = 0
C
C              4. If not diagonal then diagonalize
C
C
C     Written by Ove Christiansen 18-10-1996
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
      DIMENSION AMAT(3,3),WI(3),V1(3),FV1(3)
      LOGICAL D12,D13,D23,D1122,D1133,D2233,LBD
C
      APAR   = 0.0D0
      APEN   = 0.0D0
      XKAPPA = 0.0D0
      XX = TNSR(1,1)
      YX = TNSR(2,1)
      ZX = TNSR(3,1)
      XY = TNSR(1,2)
      YY = TNSR(2,2)
      ZY = TNSR(3,2)
      XZ = TNSR(1,3)
      YZ = TNSR(2,3)
      ZZ = TNSR(3,3)
C
      ALFSQ  = (XX+YY+ZZ)**2/9.0D0
      BETSQ  = ((XX-YY)**2+(XX-ZZ)**2+(YY-ZZ)**2 +
     *         3*(XY**2+YX**2+XZ**2+XZ**2+YZ**2+ZY**2))/2.0D0
C
      IF ((ABS(XY-YX)+ABS(XZ-ZX)+ABS(YZ-ZY)).GT.THR) THEN
C
         WRITE(LUPRI,'(/,1X,A)')
     *          'Tensor is not symmetric on input in TNSRAN'
         WRITE(LUPRI,'(1X,A,/,1X,A)')
     * 'I will symmetrice it for you to get a real symmmetric ',
     * 'tensor according to:  2*AlfaXY(om) = <<X,Y>>(om)+<<X,Y>>(-om)'
         WRITE(LUPRI,'(1X,A)')
     * '                                   = <<X,Y>>(om)+<<Y,X>>(om) '
         CALL POLSYM(TNSR,0.5D0)
         WRITE(LUPRI,'(1X,A)') 'Tensor is now: '
         CALL OUTPUT(TNSR,1,3,1,3,3,3,1,LUPRI)
C
      ENDIF
C
      D12 = (ABS(XY) .GT. THR )
      D13 = (ABS(XZ) .GT. THR )
      D23 = (ABS(YZ) .GT. THR )
C
      ITST = 0
      IF (D12) ITST = ITST + 1
      IF (D13) ITST = ITST + 1
      IF (D23) ITST = ITST + 1
C
      ITST = ITST*2
C
      IF ( ITST .EQ. 0 ) THEN
C
C------------------------------------
C        Section for diagonal tensor.
C------------------------------------
C
         PVAL(1) = TNSR(1,1)
         PVAL(2) = TNSR(2,2)
         PVAL(3) = TNSR(3,3)
         CALL DUNIT(PAXIS,3)
C
C------------------------------------------------
C        determine number of equivalent elements.
C------------------------------------------------
C
         D1122 = (ABS(XX-YY) .LT. THR )
         D1133 = (ABS(XX-ZZ) .LT. THR )
         D2233 = (ABS(YY-ZZ) .LT. THR )
C
         ITST2 = 0
         IF (D1122) ITST2 = ITST2 + 1
         IF (D1133) ITST2 = ITST2 + 1
         IF (D2233) ITST2 = ITST2 + 1
C
         IF (ITST2 .EQ. 3) THEN
C
            IF (IPRINT .GT. 10) WRITE(LUPRI,'(/,1X,A,/)')
     *           'TNSRAN: Tensor is spherical symmetric.'
            APAR = ZZ
            APEN = XX
            IPAR = 8
C
         ELSE IF (ITST2 .EQ. 1) THEN
C
            IF (IPRINT .GT. 10) WRITE(LUPRI,'(/,1X,A,/)')
     *           'TNSRAN: Tensor has cylinder symmetry.'
            IF ( D1122 ) THEN
               APAR = ZZ
               APEN = XX
               IPAR = 3
            ENDIF
            IF ( D1133 ) THEN
               APAR = YY
               APEN = XX
               IPAR = 2
            ENDIF
            IF ( D2233 ) THEN
               APAR = XX
               APEN = YY
               IPAR = 1
            ENDIF
            XKAPPA = (APAR - APEN)/(3*SQRT(ALFSQ))
C
         ELSE IF (ITST2 .EQ. 0) THEN
C
            IF (IPRINT .GT. 10) WRITE(LUPRI,'(/,1X,A,/)')
     *         'TNSRAN: Tensor is a diagonal asym. top.'
            IF (ABS(ZZ).LT.THR) IPAR=4
            IF (ABS(YY).LT.THR) IPAR=5
            IF (ABS(XX).LT.THR) IPAR=6
C
         ENDIF
C
      ELSE
         IF (IPRINT .GT. 10) WRITE(LUPRI,'(/,1X,A,I2,A,/)')
     *        'TNSRAN: Tensor has ',ITST,
     *         ' out of diagonal elements'
C
         IF (ITST .EQ. 2) THEN
            LBD = .TRUE.
            IF (D12) ISPAC = 3
            IF (D13) ISPAC = 2
            IF (D23) ISPAC = 1
            IF (D12) IPAR  = 4
            IF (D13) IPAR  = 5
            IF (D23) IPAR  = 6
         ENDIF
C
         MATZ = 1
         CALL DCOPY(3*3,TNSR,1,AMAT,1)
         CALL RG(3,3,AMAT,PVAL,WI,MATZ,PAXIS,V1,FV1,IERR)
         CALL RGORD(3,3,PVAL,WI,PAXIS,.FALSE.)
      ENDIF
C
C------------------------------------
C     A little Self consistency test.
C------------------------------------
C
      XX = PVAL(1)
      YY = PVAL(2)
      ZZ = PVAL(3)
      ALFSQ2  = (XX+YY+ZZ)**2/9.0D0
      BETSQ2  = (((XX-YY)**2+(XX-ZZ)**2+(YY-ZZ)**2)/2.0D0)
C
      IF ((ABS(ALFSQ-ALFSQ2).GT.THR).OR.(ABS(ALFSQ-ALFSQ2).GT.THR))
     *        THEN
         WRITE(LUPRI,'(/,1X,A)') 'Rotational invariants before '
     *      //'and after diagonalization is '
         WRITE(LUPRI,'(1X,A,2F15.10)') 'Alfa**2',ALFSQ,ALFSQ2
         WRITE(LUPRI,'(1X,A,2F15.10)') 'Beta**2',BETSQ,BETSQ2
         WRITE(LUPRI,'(1X,A)') 'Check the diagonalization'
      ENDIF
C
      IF ((ABS(XX).GT.THR).AND.(ABS(YY).GT.THR).AND.(ABS(ZZ).GT.THR))
     *  IPAR = 7
C
      END
c*DECK CC_PABCON
      SUBROUTINE CC_PABCON(LABELA,ISYMA,FREQA,LRLXA,
     *                     LABELB,ISYMB,FREQB,LRLXB,
     *                     PRP,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Calculate T-barA(-omeg)*Tbar-B(omeg)*P contribution to LRF.
C
C     Written by Ove Christiansen May 1998 - based on CC_FABCON
C     (for that reason somethings are called R that really are L
C      and F instead of P)
C
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <leinf.h>
C
      PARAMETER( TWO = 2.0D00,HALF=0.5D00,TOLFRQ=1.0D-08 )
      DIMENSION WORK(LWORK)
      CHARACTER LABELA*8,LABELB*8,MODEL*10
      LOGICAL LRLXA,LRLXB
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_FABCON: Calculating polarizabilty P-cont.')
      ENDIF
C
      NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
      IF ( CCS ) NTAMPA = NT1AM(ISYMA)
      NTAMPB = NT1AM(ISYMB) + NT2AM(ISYMB)
      IF ( CCS ) NTAMPB = NT1AM(ISYMB)
      IF (ISYMA .NE. ISYMB ) CALL QUIT('Symmetry mismatch in CC_PABCON')
C
C-----------------------------------------------
C     Loop perturbations of this symmetry class.
C-----------------------------------------------
C
      KR1   = 1
      KEND1 = KR1 + NTAMPB
      LEND1 = LWORK - KEND1
C
C------------------------------
C     Get P-transformed vector.
C------------------------------
C
      KR11 = KR1
      KR12 = KR1 + NT1AM(ISYMB)
C
      CALL DZERO(WORK(KR1),NTAMPB)
      CALL CC_PTB(WORK(KR1),LABELB,ISYMB,FREQB,LRLXB,WORK(KEND1),LEND1)
C
      IF (IPRINT .GT. 40 ) THEN
         CALL AROUND( 'In CC_EATB:  P*RSP vector ' )
         CALL CC_PRP(WORK(KR1),WORK(KR1+NT1AM(ISYMB)),ISYMB,1,1)
      ENDIF
C
      IF ( DEBUG ) THEN
         XLV  = DDOT(NTAMPB, WORK(KR1),1,WORK(KR1),1)
         WRITE(LUPRI,1) 'Norm of P*Response vector:         ',XLV
      ENDIF
C
      KR2   = KEND1
      KEND2 = KR2 + NTAMPA
      LEND2 = LWORK - KEND2
      IF (LEND2.LT.0) CALL QUIT('TOO LITTLE WORKSPACE IN CC_ABFCON-2')
C
C-----------------------------------------------------------
C     Get response vectors and do the dot with the P*vector.
C-----------------------------------------------------------
C
      KR21   = KR2
      KR22   = KR2 + NT1AM(ISYMA)
      ILSTNR = IL1ZETA(LABELA,LRLXA,FREQA,ISYMA)
      IOPT   = 3
      CALL CC_RDRSP('L1',ILSTNR,ISYMA,IOPT,MODEL,WORK(KR21),
     *              WORK(KR22))
      IF ( DEBUG ) THEN
         XLV  = DDOT(NTAMPA, WORK(KR2),1,WORK(KR2),1)
         WRITE(LUPRI,1) 'Norm of Response vector:         ',XLV
      ENDIF
C
      FABCON = DDOT(NTAMPA,WORK(KR1),1,WORK(KR2),1)
      IF ( IPRINT .GT. 9 ) THEN
         WRITE(LUPRI,*) ' Singles contribution:',
     *      DDOT(NT1AM(ISYMA),WORK(KR1),1,WORK(KR2),1)
         IF (.NOT. CCS) WRITE(LUPRI,*) ' Doubles contribution:',
     *      DDOT(NT2AM(ISYMA),WORK(KR1+NT1AM(ISYMA)),1,
     *      WORK(KR2+NT1AM(ISYMA)),1)
      ENDIF
      IF (IPRINT .GT. 2 ) THEN
         WRITE(LUPRI,'(1X,A2,A8,A1,A8,A3,F10.6,A,F10.6)')
     *   '<<',LABELA,',',LABELB,'>>(',
     *   FREQB,') LA*LB*P cont. = ',FABCON
      ENDIF
      PRP       = PRP       - FABCON
C
   1  FORMAT(1x,A35,1X,E20.10)
      RETURN
      END
c*DECK CC_PRPC
       SUBROUTINE CC_PRPC(PROP,LABEL,NORD,LABX,LABY,LABZ,LABU,
     *                   FRQY,FRQZ,FRQU,ISYMIN,ISYMEX,ISPINEX,IEX)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Add response property to list of property information to be
C              passed to numerical differentiantion/averaging.
C
C     Ove Christiansen August 1999.
C
C     NORD = 1    exp. value
C            2    Linear response function
C            3    Quadratic response function
C            4    Cubic response function
C           -1    ground - excited  transition matrix element
C           -2    excited - excited transition matrix element (not implemented yet)
C           -3    ground - excited transition strength
C           -4    excited - excited transition strength (not implemented yet)
C           -11    First order excited state property
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
C
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <prpc.h>
#include <ccinftap.h>
C
      LOGICAL EXIST,L1,L2,L3,L4,LI1,LI2
      PARAMETER (TOLFRQ=1.0D-08,ONE=1.0D0,XMONE=-1.0D0,TOLEXCI =1.0D-02)
C
      CHARACTER LABEL*10, LABX*8, LABY*8, LABZ*8, LABU*8
C
C--------------------------------------------------
C     Test if this property is already on the list.
C     In that case find address else update NPRPC
C--------------------------------------------------
C
C
      IF (NOEONL .AND. (NORD.EQ.0)) THEN
C         if energy and NOEONList = true then skip addition to list.
        RETURN
      ELSE
        EXIST = .FALSE.
        IF (EXIST) THEN
c          IPRPC = IHIT
        ELSE
           NPRPC = NPRPC + 1
           IPRPC = NPRPC
        ENDIF
C
        WRITE(LUPRPC,
     *   '(I5,I3,I4,1X,A10,E23.16,4(1X,A8),3E23.16,3I4)')
     *   IPRPC,ISYMIN,NORD,LABEL,PROP,
     *   LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYMEX,ISPINEX,IEX
      ENDIF
C
      END
*---------------------------------------------------------------------*
      SUBROUTINE CC_AVE2(VALUE,IDLSTX,IDLSTY,WORK,LWORK)
C-----------------------------------------------------------------------
C     Purpose: Calculate <HF|[[H,T^x],T^y]+[X,T^y]+[Y,T^x]|CC> 
C              contribution to second order property.
C              IDLSTX,IDLSTY - indeces of first-order amplitudes
C     Written by Christof Haettig, Mai 2003
C-----------------------------------------------------------------------
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccr1rsp.h>

      INTEGER ISYM0
      PARAMETER ( ISYM0 = 1 )
 
      CHARACTER LISTR1*3, LISTR2*3, MODEL*10, LABELX*8, LABELY*8
      INTEGER IDLSTX, IDLSTY, LWORK

#if defined (SYS_CRAY)  
      REAL WORK(LWORK), VALUE, DDOT
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION WORK(LWORK), VALUE, DDOT
      DOUBLE PRECISION ZERO, ONE, TWO
#endif
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
  
      INTEGER ISYMX, ISYMY, ISYMXY, KT1AM0, KLAMP0, KLAMH0, KT1AMX,
     &        KT1AMY, KEND1, LWRK1, IOPT, KFOCKX, KXIA, KEND2, LWRK2,
     &        KFOCKY, KYIA, KXIAJB, KT1AM, IRREP, ISYMM, IERR


      VALUE = ZERO

      ISYMX  = ISYLRT(IDLSTX)
      ISYMY  = ISYLRT(IDLSTY)
      ISYMXY = MULD2H(ISYMX,ISYMY)

      LABELX = LRTLBL(IDLSTX)
      LABELY = LRTLBL(IDLSTY)

      IF (ISYMXY.NE.1) RETURN
C
      KT1AM0 = 1
      KLAMP0 = KT1AM0 + NT1AM(ISYM0)
      KLAMH0 = KLAMP0 + NLAMDT
      KT1AMX = KLAMH0 + NLAMDT
      KT1AMY = KT1AMX + NT1AM(ISYMX)
      KEND1  = KT1AMY + NT1AM(ISYMY)
      LWRK1  = LWORK  - KEND1
      IF (LWRK1.LT.0)CALL QUIT(' Too little workspace in CC_AVE2')

C     -----------------------------------------------------------
C     read amplitudes:
C     -----------------------------------------------------------
      IOPT = 1
      CALL CC_RDRSP('R1',IDLSTX,ISYMX,IOPT,MODEL,WORK(KT1AMX),DUMMY)
      CALL CC_RDRSP('R1',IDLSTY,ISYMY,IOPT,MODEL,WORK(KT1AMY),DUMMY)

      ! read zeroth-order singles amplitudes and compute Lambda
      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AM0),DUMMY)
      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AM0),
     *            WORK(KEND1),LWRK1)
 
C     -----------------------------------------------------------
C     compute <HF|[X,T^Y]|HF>
C     -----------------------------------------------------------
      KFOCKX = KEND1
      KXIA   = KFOCKX + N2BST(ISYMX)
      KEND2  = KXIA   + NT1AM(ISYMX)
      LWRK2  = LWORK  - KEND2
      IF (LWRK2.LT.0)CALL QUIT(' Too little workspace in CC_AVE2')

      ! get X integrals:
      CALL CCPRPAO(LABELX,.TRUE.,WORK(KFOCKX),IRREP,ISYMM,IERR,
     &              WORK(KEND2),LWRK2)
      IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMX)) THEN
        CALL QUIT('CC_AVE2: error reading operator '//LABELX)
      ELSE IF (IERR.LT.0) THEN
        CALL DZERO(WORK(KFOCKX),N2BST(ISYMX))
      END IF
      CALL CC_FCKMO(WORK(KFOCKX),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND2),LWRK2,ISYMX,1,1)
      CALL CC_FOCK_RESORT(DUMMY,.FALSE.,WORK(KXIA),.TRUE.,
     &      DUMMY,.FALSE.,DUMMY,.FALSE.,WORK(KFOCKX),ISYMX) 
 

      VALUE = VALUE + TWO * 
     &   DDOT(NT1AM(ISYMX),WORK(KXIA),1,WORK(KT1AMY),1)

C     -----------------------------------------------------------
C     compute <HF|[Y,T^X]|HF>
C     -----------------------------------------------------------
      KFOCKY = KEND1
      KYIA   = KFOCKY + N2BST(ISYMX)
      KEND2  = KYIA   + NT1AM(ISYMX)
      LWRK2  = LWORK  - KEND2
      IF (LWRK2.LT.0)CALL QUIT(' Too little workspace in CC_AVE2')

      ! get Y integrals:
      CALL CCPRPAO(LABELY,.TRUE.,WORK(KFOCKY),IRREP,ISYMM,IERR,
     &              WORK(KEND2),LWRK2)
      IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMY)) THEN
        CALL QUIT('CC_AVE2: error reading operator '//LABELY)
      ELSE IF (IERR.LT.0) THEN
        CALL DZERO(WORK(KFOCKY),N2BST(ISYMY))
      END IF
      CALL CC_FCKMO(WORK(KFOCKY),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND2),LWRK2,ISYMY,1,1)
      CALL CC_FOCK_RESORT(DUMMY,.FALSE.,WORK(KYIA),.TRUE.,
     &      DUMMY,.FALSE.,DUMMY,.FALSE.,WORK(KFOCKY),ISYMY) 
 
      VALUE = VALUE + TWO * 
     &   DDOT(NT1AM(ISYMX),WORK(KYIA),1,WORK(KT1AMX),1)

C     -----------------------------------------------------------
C     get packed L(ia,jb) integrals and evaluate the 
C     projection contribution <HF|[[H,T^X],T^Y]|CC>
C     -----------------------------------------------------------
      KXIAJB = KEND1
      KXIA   = KXIAJB + NT2AM(ISYM0)
      KEND2  = KXIA   + NT1AM(ISYMX)
      LWRK2  = LWORK  - KEND2
      IF (LWRK2.LT.0)CALL QUIT(' Too little workspace in CC_AVE2')

      CALL CCG_RDIAJB(WORK(KXIAJB),NT2AM(ISYM0))

      IOPT = 1
      Call CCSD_TCMEPK(WORK(KXIAJB),ONE,ISYM0,IOPT) 

      IOPT = 0
      CALL DZERO(WORK(KXIA),NT1AM(ISYMX))
      CALL CCG_LXD(WORK(KXIA),ISYMX,WORK(KT1AMX),ISYMX,
     &             WORK(KXIAJB),ISYM0,IOPT)

      VALUE = VALUE + TWO * 
     &   DDOT(NT1AM(ISYMX),WORK(KXIA),1,WORK(KT1AMY),1)
 
      RETURN
      END
*---------------------------------------------------------------------*
      SUBROUTINE CC_TSTAV2(IDLSTR2,VEC,WORK,LWORK,IOPTTST)
C----------------------------------------------------------------------
C     Purpose: Calculate second-order properties from the second-order
C              amplitude response to test these
C     Written by Christof Haettig, May 2003
C----------------------------------------------------------------------
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccr1rsp.h>
#include <ccr2rsp.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IDLSTR2, LWORK, IOPTTST, ISYM0
      PARAMETER (ISYM0 = 1)

#if defined (SYS_CRAY)  
      REAL PROPAVE, PROPRSP, WORK(*), VEC(*), DDOT
#else
      DOUBLE PRECISION PROPAVE, PROPRSP, WORK(*), VEC(*), DDOT
#endif

      LOGICAL LORX
      CHARACTER*10 MODEL
      INTEGER KETA, KEND1, LWRK1, IOPT, IDLSTX, IDLSTY,
     &        IR1TAMP, NVAR

      IF (CCS) THEN
         PROPRSP = 0.0D0
      ELSE
         NVAR   = NT1AM(ISYM0) + NT2AM(ISYM0)
         IF (CCR12) THEN
           NVAR = NVAR + NTR12AM(ISYM0)
         ENDIF
         KETA   = 1 
         KEND1  = KETA  + NVAR
         LWRK1  = LWORK - KEND1 
         IF (LWRK1.LT.0) CALL QUIT('Too little workspace in CC_TSTAV2')
         IF      (IOPTTST.EQ.0) THEN
            CALL CC_ETA(WORK(KETA),WORK(KEND1),LWRK1)
         ELSE IF (IOPTTST.EQ.1) THEN
            IOPT  = 3
            CALL CC_RDRSP('L0 ',0,ISYM0,IOPT,MODEL,WORK(KETA),
     *                    WORK(KETA+NT1AM(ISYMOP)))
            IF (CCR12) THEN
              IOPT = 32
              CALL CC_RDRSP('L0 ',0,ISYMOP,IOPT,MODEL,DUMMY,
     *                      WORK(KETA+NT1AM(ISYMOP)+NT2AM(ISYMOP)))
            ENDIF
         ELSE
            WRITE(LUPRI,*) 'IOPTTST = ',IOPTTST
            CALL QUIT('ILLEGAL VALUE FOR IOPTTST IN CC_TSTAV2.')
         END IF
         PROPRSP = DDOT(NVAR,WORK(KETA),1,VEC,1)
 
         IF (LOCDBG) THEN
           write(lupri,*) 'Input vector:'
           call cc_prp(vec,vec(nt1am(isymop)+1),isymop,1,1)
           if (CCR12) call cc_prpr12(vec(1+nt1am(isymop)+nt2am(isymop)),
     *                               isymop,1,.false.)
           write(lupri,*) 'L0/X0 vector:'
           call cc_prp(work(keta),work(keta+nt1am(isymop)),isymop,1,1)
           if (CCR12) call cc_prpr12(work(keta+nt1am(isymop)+
     *                               nt2am(isymop)),isymop,1,.false.)
           write(lupri,*) 'PROPRSP:',PROPRSP
         END IF
      ENDIF

      IDLSTX = IR1TAMP(LBLR2T(IDLSTR2,1),LORXR2T(IDLSTR2,1),
     &                 FRQR2T(IDLSTR2,1), ISYR2T(IDLSTR2,1))
      IDLSTY = IR1TAMP(LBLR2T(IDLSTR2,2),LORXR2T(IDLSTR2,2),
     &                 FRQR2T(IDLSTR2,2), ISYR2T(IDLSTR2,2))

      LORX = LORXR2T(IDLSTR2,1) .OR. LORXR2T(IDLSTR2,2)

      IF ( LORX ) THEN
        CALL QUIT('No relaxation implemented in CC_TSTAV2.')
      ELSE
        ! if it is a simple unrelaxed one-electron perturbation
        ! calculate the average value contribution in CC_AVE
        CALL CC_AVE2(PROPAVE,IDLSTX,IDLSTY,WORK,LWORK)
      END IF

      WRITE(LUPRI,'(1X,3A)') 'Operators   : ',
     *  LBLR2T(IDLSTR2,1),LBLR2T(IDLSTR2,2)
      WRITE(LUPRI,'(1X,A,2F16.4)') 'Frequencies : ',
     *  FRQR2T(IDLSTR2,1),FRQR2T(IDLSTR2,2)
      WRITE(LUPRI,'(1X,A,F16.10)') 'Average contribution:   ',
     *                         PROPAVE
      WRITE(LUPRI,'(1X,A,F16.10)') 'Response contribution:  ',
     *                         PROPRSP
      WRITE(LUPRI,'(1X,A,F16.10)') 'Total second-order property:',
     *                         PROPAVE + PROPRSP

      RETURN
      END 
*---------------------------------------------------------------------*
C  /* Deck cc_rotpri */
      SUBROUTINE CC_ROTPRI(RIN,STREN,EIGV,IEX,ISYM,MODEL,LCALC,LDIP,
     &                     LUOSC)
C
C     Thomas Bondo Pedersen, January 2005.
C     - based on CC_OSCPRI by Ove Christiansen.
C
C     Purpose: Print rotatory strengths.
C
#include <implicit.h>
      DIMENSION    RIN(3)
      CHARACTER*10 MODEL
      LOGICAL      LCALC
#include <priunit.h>
#include <pgroup.h>
#include <codata.h>
#include <ccsdinp.h>

      PARAMETER (RAUSI  = FPEPS0*(HBAR**3)*1.0D55/(EMASS**2))
      PARAMETER (RAUCGS = ECHARGE*ECHARGE*HBAR*XTANG*CCM*1.0D36/EMASS)

      DIMENSION   ROT(3), STR(3)
      INTEGER     POL(3)
      CHARACTER*7 CDIP

      PARAMETER (THRPOL = 1.0D-8)  ! Same threshold as in TNSRAN for polarization...

      DATA POL /1,10,100/

      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_ROTPRI: Output Rotatory Strengths ' )
      END IF

      IMULT = 1  ! force singlet spin symmetry...

      IF (LCALC) THEN

C-tbp: ANGMOM sign fixed here:
         CALL DSCAL(3,-1.0D0,RIN,1)

         CALL DCOPY(3,RIN,1,ROT,1)

         WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     &   'Rotatory strength for state nr.',IEX,
     &   ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
         IF (LDIP .EQ. 1) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: length'
            FACT = -0.5D0
         ELSE IF (LDIP .EQ. 2) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: velocity'
            IF (ABS(EIGV) .LT. 1.0D-8) THEN
               FACT = -1.0D16
            ELSE
               FACT = -1.0D0/(2.0D0*EIGV)
            END IF
         ELSE
            WRITE(LUPRI,'(3X,A)') 'Gauge: UNKNOWN'
            WRITE(LUPRI,'(3X,A)') '- scaling factors will be incorrect!'
            FACT = 1.0D0
         ENDIF
         CALL DZERO(STR,3)
         CALL DSCAL(3,FACT,ROT,1)
         DO I = 1,3
            STR(1) = STR(1) + ROT(I)
         END DO
         STR(2) = RAUSI*STR(1)
         STR(3) = RAUCGS*STR(1)
         WRITE(LUPRI,'(/,3X,A)') 'Rotatory strength components (a.u.):'
         WRITE(LUPRI,'(10X,A1,15X,A1,15X,A1)') 'X','Y','Z'
         WRITE(LUPRI,'(3X,F15.10,1X,F15.10,1X,F15.10,/)')
     &   ROT(1),ROT(2),ROT(3)
         WRITE(LUPRI,'(3X,A,F15.7,/,3X,A,F15.7,/,3X,A,F15.7)')
     &   'Total Rotatory Strength in Atomic Units      : ',STR(1),
     &   'Total Rotatory Strength in 10-55   A^2 m^3 s : ',STR(2),
     &   'Total Rotatory Strength in 10-40 cm^5 g s^-2 : ',STR(3)

         STREN = STR(1)

         IPOL = 0
         DO I = 1,3
            IF (ABS(ROT(I)) .GT. THRPOL) IPOL = IPOL + POL(I)
         END DO
         IF (IPOL .EQ.   1) THEN
            CDIP = '   X   '
         ELSE IF (IPOL .EQ.  10) THEN
            CDIP = '   Y   '
         ELSE IF (IPOL .EQ. 100) THEN
            CDIP = '   Z   '
         ELSE IF (IPOL .EQ.  11) THEN
            CDIP = ' (X,Y) '
         ELSE IF (IPOL .EQ. 101) THEN
            CDIP = ' (X,Z) '
         ELSE IF (IPOL .EQ. 110) THEN
            CDIP = ' (Y,Z) '
         ELSE IF (IPOL .EQ. 111) THEN
            CDIP = '(X,Y,Z)'
         ELSE
            CDIP = '   -   '
         ENDIF

         WRITE(LUOSC,9988) IMULT,REP(ISYM-1),IEX,STR(2),STR(3),CDIP

         CALL FLSHFO(LUPRI)

      ELSE

         CDIP = '   ?   '
         WRITE(LUOSC,9986) IMULT,REP(ISYM-1),IEX,'Not calculated',
     &                     'Not calculated',CDIP

      END IF

      RETURN

 9986 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,'  |')
 9987 FORMAT(1X,'|       | ',I4,'   | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,'  |')
 9988 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9989 FORMAT(1X,'|       | ',I4,'   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')

      END
C  /* Deck cc_rtqpri */
      SUBROUTINE CC_RTQPRI(RQIN,RQOUT,EIGV,IEX,ISYM,MODEL,LCALC,LDIP,
     &                     LUOSC,NWAR)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Print rotatory strength tensors, el. quadrupole contribution.
C
#include <implicit.h>
      DIMENSION    RQIN(3,9), RQOUT(3,3)
      CHARACTER*10 MODEL
      LOGICAL      LCALC
#include <priunit.h>
#include <pgroup.h>
#include <codata.h>
#include <ccsdinp.h>

      PARAMETER (RAUSI  = FPEPS0*(HBAR**3)*1.0D55/(EMASS**2))
      PARAMETER (RAUCGS = ECHARGE*ECHARGE*HBAR*XTANG*CCM*1.0D36/EMASS)

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_RTQPRI')

      DIMENSION   SQ(3,3,3), AVE(3)
      DIMENSION   RQ(3,3,3)
      CHARACTER*7 CDIP
      LOGICAL     WARN
#if defined (SYS_CRAY)
      REAL              LEVICI(3,3,3)
#else
      DOUBLE PRECISION  LEVICI(3,3,3)
#endif

      PARAMETER (TINY = 1.0D-12)

      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND('IN CC_RTQPRI: El. Quadr. Rotatory'
     &               //' Strength Tensors')
      END IF

      IMULT = 1  ! force singlet spin symmetry...

      IF (LCALC) THEN

         CALL DCOPY(3*9,RQIN,1,SQ,1)

         CALL DZERO(LEVICI,3*3*3)
         LEVICI(1,2,3) = 1.0D0
         LEVICI(2,1,3) = -1.0D0
         LEVICI(3,1,2) = 1.0D0
         LEVICI(1,3,2) = -1.0D0
         LEVICI(2,3,1) = 1.0D0
         LEVICI(3,2,1) = -1.0D0

         WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     &   'El. quadr. rotatory strength tensor for state nr.',IEX,
     &   ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
         IF (LDIP .EQ. 1) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: length'
            FACT = -3.0D0*EIGV/4.0D0
         ELSE IF (LDIP .EQ. 2) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: velocity'
            IF (ABS(EIGV) .LT. 1.0D-8) THEN
               FACT = -1.0D16
            ELSE
               FACT = 3.0D0/(4.0D0*EIGV)
            END IF
         ELSE
            WRITE(LUPRI,'(3X,A)') 'Gauge: UNKNOWN'
            WRITE(LUPRI,'(3X,A)') '- scaling factors will be incorrect!'
            FACT = 1.0D0
         ENDIF

         IERR = 0
         DO I = 1,3
            JERR = 0
            DO J = 1,3
               DO K = 1,J
                  JK = 3*(K - 1) + J
                  KJ = 3*(J - 1) + K
                  DIFF = ABS(RQIN(I,JK) - RQIN(I,KJ))
                  IF (DIFF .GT. 1.0D-14) JERR = JERR + 1
               END DO
            END DO
            IERR = IERR + JERR
         END DO
         IF (IERR .NE. 0) THEN
            WRITE(LUPRI,*) SECNAM,': non-symmetric rank-3 tensor',
     &                     ' on entry'
            WRITE(LUPRI,*) 'This will lead to non-zero average!!!'
            WRITE(LUPRI,*) 'Residues from input:'
            WRITE(LUPRI,'(1X,A,F12.8)') 'X,XX: ',RQIN(1,1)
            WRITE(LUPRI,'(1X,A,F12.8)') 'X,XY: ',RQIN(1,4)
            WRITE(LUPRI,'(1X,A,F12.8)') 'X,XZ: ',RQIN(1,7)
            WRITE(LUPRI,'(1X,A,F12.8)') 'X,YY: ',RQIN(1,5)
            WRITE(LUPRI,'(1X,A,F12.8)') 'X,YZ: ',RQIN(1,8)
            WRITE(LUPRI,'(1X,A,F12.8)') 'X,ZZ: ',RQIN(1,9)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Y,XX: ',RQIN(2,1)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Y,XY: ',RQIN(2,4)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Y,XZ: ',RQIN(2,7)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Y,YY: ',RQIN(2,5)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Y,YZ: ',RQIN(2,8)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Y,ZZ: ',RQIN(2,9)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Z,XX: ',RQIN(3,1)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Z,XY: ',RQIN(3,4)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Z,XZ: ',RQIN(3,7)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Z,YY: ',RQIN(3,5)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Z,YZ: ',RQIN(3,8)
            WRITE(LUPRI,'(1X,A,F12.8)') 'Z,ZZ: ',RQIN(3,9)
            CALL QUIT('Error in '//SECNAM)
         END IF

         CALL DSCAL(3*3*3,FACT,SQ,1)
         CALL DZERO(RQ,3*3*3)
         DO K = 1,3
            DO J = 1,3
               DO M = 1,3
                  DO L = 1,3
                     RQ(J,K,1) = RQ(J,K,1)
     &                         + LEVICI(L,M,J)*SQ(L,M,K)
                  END DO
               END DO
            END DO
         END DO
         CALL POLSYM(RQ(1,1,1),0.5D0)
         CALL DAXPY(3*3,RAUSI,RQ(1,1,1),1,RQ(1,1,2),1)
         CALL DAXPY(3*3,RAUCGS,RQ(1,1,1),1,RQ(1,1,3),1)
         CALL DZERO(AVE,3)
         DO I = 1,3
            DO J = 1,3
               AVE(I) = AVE(I) + RQ(J,J,I)
            END DO
            AVE(I) = AVE(I)/3.0D0
         END DO
         WRITE(LUPRI,'(/,3X,A)')
     & 'Electric quadrupole rotatory strength tensor components (a.u.):'
         CALL OUTPUT(RQ(1,1,1),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)')
     &   'Orientational average: ',AVE(1)
         WRITE(LUPRI,'(/,3X,A,A)')
     &   'Electric quadrupole rotatory strength tensor components ',
     &   '(D-55 SI):'
         CALL OUTPUT(RQ(1,1,2),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)')
     &   'Orientational average: ',AVE(2)
         WRITE(LUPRI,'(/,3X,A,A)')
     &   'Electric quadrupole rotatory strength tensor components ',
     &   '(D-40 cgs):'
         CALL OUTPUT(RQ(1,1,3),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)')
     &   'Orientational average: ',AVE(3)
         DIFF = AVE(1)
         WARN = ABS(DIFF) .GT. TINY
         IF (WARN) THEN
            WRITE(LUPRI,9990)
            NWAR = NWAR + 1
         END IF

         DO J = 1,3
            DO K = J,3
               CDIP = '   ?   '
               IF ((J.EQ.1) .AND. (K.EQ.1)) THEN
                  CDIP = '  XX   '
                  WRITE(LUOSC,9988) IMULT,REP(ISYM-1),IEX,
     &                              RQ(J,K,2),RQ(J,K,3),CDIP
               ELSE
                  IF (J .EQ. 1) THEN
                     IF (K .EQ. 2) THEN
                        CDIP = '  XY   '
                     ELSE IF (K .EQ. 3) THEN
                        CDIP = '  XZ   '
                     END IF
                  ELSE IF (J .EQ. 2) THEN
                     IF (K .EQ. 2) THEN
                        CDIP = '  YY   '
                     ELSE IF (K .EQ. 3) THEN
                        CDIP = '  YZ   '
                     END IF
                  ELSE IF (J .EQ. 3) THEN
                     IF (K .EQ. 3) THEN
                        CDIP = '  ZZ   '
                     END IF
                  END IF
                  WRITE(LUOSC,9987) RQ(J,K,2),RQ(J,K,3),CDIP
               END IF
            END DO
         END DO

         CALL DCOPY(3*3,RQ(1,1,1),1,RQOUT,1)

      ELSE

         CDIP = '   ?   '
         WRITE(LUOSC,9986) IMULT,REP(ISYM-1),IEX,'Not calculated',
     &                     'Not calculated',CDIP

      END IF

      RETURN

 9986 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,'  |')
 9987 FORMAT(1X,'|      ',' |     ','   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9988 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9990 FORMAT(1X,'***WARNING*** Incorrect average!!!')

      END
C  /* Deck cc_rtmpri */
      SUBROUTINE CC_RTMPRI(RMIN,RMOUT,EIGV,IEX,ISYM,MODEL,LCALC,LDIP,
     &                     LUOSC,CHKSTR,NWAR)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Print rotatory strength tensors, magn. dipole contribution.
C
#include <implicit.h>
      DIMENSION    RMIN(3,3), RMOUT(3,3)
      CHARACTER*10 MODEL
      LOGICAL      LCALC
#include <priunit.h>
#include <pgroup.h>
#include <codata.h>
#include <ccsdinp.h>

      PARAMETER (RAUSI  = FPEPS0*(HBAR**3)*1.0D55/(EMASS**2))
      PARAMETER (RAUCGS = ECHARGE*ECHARGE*HBAR*XTANG*CCM*1.0D36/EMASS)

      DIMENSION   SM(3,3), AVE(3)
      DIMENSION   RM(3,3,3)
      CHARACTER*7 CDIP
      LOGICAL     WARN

      PARAMETER (TINY = 1.0D-12)

      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND('IN CC_RTMPRI: Magn. Dip. Rotatory'
     &               //' Strength Tensors')
      END IF

      IMULT = 1  ! force singlet spin symmetry...

      IF (LCALC) THEN

C-tbp: ANGMOM sign fixed here:
         CALL DSCAL(3*3,-1.0D0,RMIN,1)

         TRA = 0.0D0
         DO K = 1,3
            DO J = 1,3
               SM(J,K) = -RMIN(K,J)
            END DO
            TRA = TRA + RMIN(K,K)
         END DO

         WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     &   'Magn. dip. rotatory strength tensor for state nr.',IEX,
     &   ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
         IF (LDIP .EQ. 1) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: length'
            FACT = -0.75D0
         ELSE IF (LDIP .EQ. 2) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: velocity'
            IF (ABS(EIGV) .LT. 1.0D-8) THEN
               FACT = -1.0D16
            ELSE
               FACT = -3.0D0/(4.0D0*EIGV)
            END IF
         ELSE
            WRITE(LUPRI,'(3X,A)') 'Gauge: UNKNOWN'
            WRITE(LUPRI,'(3X,A)') '- scaling factors will be incorrect!'
            FACT = 1.0D0
         ENDIF
         TRA = TRA*FACT
         CALL DSCAL(3*3,FACT,SM,1)
         CALL DZERO(RM,3*3*3)
         DO K = 1,3
            DO J = 1,3
               RM(J,K,1) = SM(J,K)
            END DO
            RM(K,K,1) = RM(K,K,1) + TRA
         END DO
         CALL POLSYM(RM(1,1,1),0.5D0)
         CALL DAXPY(3*3,RAUSI,RM(1,1,1),1,RM(1,1,2),1)
         CALL DAXPY(3*3,RAUCGS,RM(1,1,1),1,RM(1,1,3),1)
         CALL DZERO(AVE,3)
         DO I = 1,3
            DO J = 1,3
               AVE(I) = AVE(I) + RM(J,J,I)
            END DO
            AVE(I) = AVE(I)/3.0D0
         END DO
         WRITE(LUPRI,'(/,3X,A)')
     &   'Magnetic dipole rotatory strength tensor components (a.u.):'
         CALL OUTPUT(RM(1,1,1),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)')
     &   'Orientational average: ',AVE(1)
         WRITE(LUPRI,'(/,3X,A,A)')
     &   'Magnetic dipole rotatory strength tensor components ',
     &   '(D-55 SI):'
         CALL OUTPUT(RM(1,1,2),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)')
     &   'Orientational average: ',AVE(2)
         WRITE(LUPRI,'(/,3X,A,A)')
     &   'Magnetic dipole rotatory strength tensor components ',
     &   '(D-40 cgs):'
         CALL OUTPUT(RM(1,1,3),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)')
     &   'Orientational average: ',AVE(3)
         DIFF = AVE(1) - CHKSTR
         WARN = ABS(DIFF) .GT. TINY
         IF (WARN) THEN
            WRITE(LUPRI,9990)
            NWAR = NWAR + 1
         END IF

         DO J = 1,3
            DO K = J,3
               CDIP = '   ?   '
               IF ((J.EQ.1) .AND. (K.EQ.1)) THEN
                  CDIP = '  XX   '
                  WRITE(LUOSC,9988) IMULT,REP(ISYM-1),IEX,
     &                              RM(J,K,2),RM(J,K,3),CDIP
               ELSE
                  IF (J .EQ. 1) THEN
                     IF (K .EQ. 2) THEN
                        CDIP = '  XY   '
                     ELSE IF (K .EQ. 3) THEN
                        CDIP = '  XZ   '
                     END IF
                  ELSE IF (J .EQ. 2) THEN
                     IF (K .EQ. 2) THEN
                        CDIP = '  YY   '
                     ELSE IF (K .EQ. 3) THEN
                        CDIP = '  YZ   '
                     END IF
                  ELSE IF (J .EQ. 3) THEN
                     IF (K .EQ. 3) THEN
                        CDIP = '  ZZ   '
                     END IF
                  END IF
                  WRITE(LUOSC,9987) RM(J,K,2),RM(J,K,3),CDIP
               END IF
            END DO
         END DO

         CALL DCOPY(3*3,RM(1,1,1),1,RMOUT,1)

      ELSE

         CDIP = '   ?   '
         WRITE(LUOSC,9986) IMULT,REP(ISYM-1),IEX,'Not calculated',
     &                     'Not calculated',CDIP

      END IF

      RETURN

 9986 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,'  |')
 9987 FORMAT(1X,'|      ',' |     ','   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9988 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9990 FORMAT(1X,'***WARNING*** Incorrect average!!!')

      END
C  /* Deck cc_rttpri */
      SUBROUTINE CC_RTTPRI(RTIN,EIGV,IEX,ISYM,MODEL,LCALC,LDIP,
     &                     LUOSC,CHKSTR,NWAR)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Print rotatory strength tensors, total.
C
#include <implicit.h>
      DIMENSION    RTIN(3,3)
      CHARACTER*10 MODEL
      LOGICAL      LCALC
#include <priunit.h>
#include <pgroup.h>
#include <codata.h>
#include <ccsdinp.h>

      PARAMETER (RAUSI  = FPEPS0*(HBAR**3)*1.0D55/(EMASS**2))
      PARAMETER (RAUCGS = ECHARGE*ECHARGE*HBAR*XTANG*CCM*1.0D36/EMASS)

      DIMENSION   RTOT(3,3,3), AVE(3)
      CHARACTER*7 CDIP
      LOGICAL     WARN

      PARAMETER (TINY = 1.0D-12)

      IMULT = 1

      IF (LCALC) THEN

         CALL DCOPY(3*3,RTIN,1,RTOT(1,1,1),1)
         CALL DZERO(RTOT(1,1,2),3*3)
         CALL DAXPY(3*3,RAUSI,RTOT(1,1,1),1,RTOT(1,1,2),1)
         CALL DZERO(RTOT(1,1,3),3*3)
         CALL DAXPY(3*3,RAUCGS,RTOT(1,1,1),1,RTOT(1,1,3),1)

         CALL DZERO(AVE,3)
         DO I = 1,3
            DO J = 1,3
               AVE(I) = AVE(I) + RTOT(J,J,I)
            END DO
            AVE(I) = AVE(I)/3.0D0
         END DO

         WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     &   'Total rotatory strength tensor for state nr.',IEX,
     &   ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
         IF (LDIP .EQ. 1) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: length'
         ELSE IF (LDIP .EQ. 2) THEN
            WRITE(LUPRI,'(3X,A)') 'Gauge: velocity'
         ELSE
            WRITE(LUPRI,'(3X,A)') 'Gauge: UNKNOWN'
         ENDIF
         WRITE(LUPRI,'(/,3X,A)')
     &   'Total rotatory strength tensor components (a.u.):'
         CALL OUTPUT(RTOT(1,1,1),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)') 'Scalar strength: ',AVE(1)
         WRITE(LUPRI,'(/,3X,A,A)')
     &   'Total rotatory strength tensor components ',
     &   '(D-55 SI):'
         CALL OUTPUT(RTOT(1,1,2),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)') 'Scalar strength: ',AVE(2)
         WRITE(LUPRI,'(/,3X,A,A)')
     &   'Total rotatory strength tensor components ',
     &   '(D-40 cgs):'
         CALL OUTPUT(RTOT(1,1,3),1,3,1,3,3,3,1,LUPRI)
         WRITE(LUPRI,'(3X,A,1P,D17.10)') 'Scalar strength: ',AVE(3)
         DIFF = AVE(1) - CHKSTR
         WARN = ABS(DIFF) .GT. TINY
         IF (WARN) THEN
            WRITE(LUPRI,9990)
            NWAR = NWAR + 1
         END IF

         DO J = 1,3
            DO K = J,3
               CDIP = '   ?   '
               IF ((J.EQ.1) .AND. (K.EQ.1)) THEN
                  CDIP = '  XX   '
                  WRITE(LUOSC,9988) IMULT,REP(ISYM-1),IEX,
     &                              RTOT(J,K,2),RTOT(J,K,3),CDIP
               ELSE
                  IF (J .EQ. 1) THEN
                     IF (K .EQ. 2) THEN
                        CDIP = '  XY   '
                     ELSE IF (K .EQ. 3) THEN
                        CDIP = '  XZ   '
                     END IF
                  ELSE IF (J .EQ. 2) THEN
                     IF (K .EQ. 2) THEN
                        CDIP = '  YY   '
                     ELSE IF (K .EQ. 3) THEN
                        CDIP = '  YZ   '
                     END IF
                  ELSE IF (J .EQ. 3) THEN
                     IF (K .EQ. 3) THEN
                        CDIP = '  ZZ   '
                     END IF
                  END IF
                  WRITE(LUOSC,9987) RTOT(J,K,2),RTOT(J,K,3),CDIP
               END IF
            END DO
         END DO

      ELSE

         CDIP = '   ?   '
         WRITE(LUOSC,9986) IMULT,REP(ISYM-1),IEX,'Not calculated',
     &                     'Not calculated',CDIP

      END IF

      RETURN

 9986 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',A16,4X,
     *       '  |',A15,5X,'  | ',A7,'  ',1X,'  |')
 9987 FORMAT(1X,'|      ',' |     ','   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9988 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',F16.7,4X,
     *       '  |',F15.7,5X,'  | ',A7,'  ',1X,'  |')
 9990 FORMAT(1X,'***WARNING*** Incorrect scalar strength!!!')

      END
C  /* Deck cc_sopr */
      SUBROUTINE CC_SOPR(WORK,LWORK)
C
C     Thomas Bondo Pedersen, January 2005.
C     - based on CC_LRESID by Ove Christiansen.
C
C     Purpose: Calculate linear response residues.
C              The Eta and Ksi vectors are calculated only once.
C
C     NOTE: it is probably better to use *CCOPA ....
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <codata.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclres.h>
#include <ccrspprp.h>
#include <ccroper.h>
#include <ccexci.h>
#include <ccexcinf.h>
#include <dummy.h>
#include <ccinftap.h>
#include <ccsections.h>

      PARAMETER (RAUSI  = FPEPS0*(HBAR**3)*1.0D55/(EMASS**2))
      PARAMETER (RAUCGS = ECHARGE*ECHARGE*HBAR*XTANG*CCM*1.0D36/EMASS)

      CHARACTER*7 SECNAM
      PARAMETER (SECNAM = 'CC_SOPR')

      LOGICAL LOCDBG,LCALC
      PARAMETER (LOCDBG = .FALSE.)
      CHARACTER*16 DBGMSG
      PARAMETER (DBGMSG = 'CC_SOPR[debug]: ')

      CHARACTER*8  LABELA, LABELB
      CHARACTER*10 MODEL, MODELP

      INTEGER NLOCS(8)

      INTEGER ILRES

      CALL QENTER(SECNAM)

C     Start timing.
C     -------------

      TIMTOT = SECOND()

C     Initialize counter (# residues).
C     --------------------------------

      NTOT = 0

C     Print header.
C     -------------

      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<  OUTPUT FROM COUPLED CLUST',
     &                               'ER LINEAR RESPONSE   >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<      CALCULATION OF SECOND',
     &                               ' ORDER RESIDUES      >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************'

      MODEL = 'CCSD      '
      IF (CC2) THEN
         MODEL = 'CC2       '
      END IF
      IF (MCC2) THEN
         MODEL = 'MCC2      '
      END IF
      IF (CCS) THEN
         MODEL = 'CCS       '
      END IF
      IF (CC3  ) THEN
         MODEL = 'CC3       '
         WRITE(LUPRI,'(/,1X,A)')
     *    'CC3 linear response residues not implemented yet'
         RETURN
      END IF
      IF (CC1A) THEN
         MODEL = 'CCSDT-1a  '
         WRITE(LUPRI,'(/,1X,A)')
     *    'CC1A linear response residues not implemented yet'
         RETURN
      END IF
      IF (CCSD) THEN
         MODEL = 'CCSD      '
      END IF

      IF (CIS) THEN
         MODELP = 'CIS      '
      ELSE
         MODELP = MODEL
      END IF

      CALL AROUND(SECNAM//': Calculation of '//MODELP//' Residues')
      IF (IPRINT .GT. 10) THEN 
         WRITE(LUPRI,*) SECNAM,': LWORK = ',LWORK
      END IF
      CALL FLSHFO(LUPRI)

C     Count number of selected states in each symmetry.
C     -------------------------------------------------

      CALL IZERO(NLOCS,NSYM)
      DO IRSD = 1,NXLRSST
         ISTATE = ILRSST(IRSD)
         ISYME  = ISYEXC(ISTATE)
         NLOCS(ISYME) = NLOCS(ISYME) + 1
      END DO

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) DBGMSG,'NLOCS: ',(NLOCS(I),I=1,NSYM)
         CALL FLSHFO(LUPRI)
      END IF

C     Check that any residues requested.
C     ----------------------------------

      NTEST = NXLRSST*NLRSOP
      IF (NTEST .LE. 0) THEN
         WRITE(LUPRI,'(/,1X,A,A)')
     &   SECNAM,': No residues requested.'
         WRITE(LUPRI,'(1X,A,I10,/,1X,A,I10,/)')
     &   'Number of selected  excited  states :',NXLRSST,
     &   'Number of requested operator doubles:',NLRSOP
         GO TO 999
      END IF

C     Allocation 1.
C     -------------

      NTRMOM = NXLRSST*NPRLBL

      KRIGHT = 1
      KLEFT  = KRIGHT + NTRMOM
      KEND1  = KLEFT  + NTRMOM
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [1]')
      END IF

C     Initialize transition moment arrays.
C     ------------------------------------

      CALL DZERO(WORK(KRIGHT),NTRMOM)
      CALL DZERO(WORK(KLEFT),NTRMOM)

C     Loop through operators in PRPLBL.
C     ---------------------------------

      DO IPRLBL = 1,NPRLBL

         LABELA = PRPLBL(IPRLBL)

C        Check that the operator enters in at least 1 residue requested.
C        ---------------------------------------------------------------

         IAB   = 1
         IOPER = ILRES(LABELA,'A')
         IF (IOPER .LE. 0) THEN
            IAB   = 2
            IOPER = ILRES(LABELA,'B')
         END IF

         IF (IOPER .GT. 0) THEN

            IF (IAB .EQ. 1) THEN
               ISYMA = ISYOPR(IALRSOP(IOPER))
            ELSE IF (IAB .EQ. 2) THEN
               ISYMA = ISYOPR(IBLRSOP(IOPER))
            ELSE
               WRITE(LUPRI,*) SECNAM,': lllegal IAB: ',IAB
               CALL QUIT('Internal error in '//SECNAM//' [IAB 1]')
            END IF

            IF (NLOCS(ISYMA) .GT. 0) THEN

C              Allocation 2.
C              -------------

               NTAMP = NT1AM(ISYMA)
               IF (.NOT. CCS) NTAMP = NTAMP + NT2AM(ISYMA)

               KETA  = KEND1
               KEND2 = KETA  + NTAMP
               LWRK2 = LWORK - KEND2 + 1

               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient memory in '//SECNAM//' [2]')
               END IF

C              Offsets to right and left moments.
C              ----------------------------------

               KOFFR = KRIGHT + NXLRSST*(IPRLBL - 1)
               KOFFL = KLEFT  + NXLRSST*(IPRLBL - 1)

C              Calculate etaA vector.
C              ----------------------

               CALL CC_ETAC(ISYMA,LABELA,WORK(KETA),'L0',1,0,DUMMY,
     &                      WORK(KEND2),LWRK2)

C              Calculate contribution to right transition moment:
C              etaA*RE for all excited states of matching symmetry.
C              ----------------------------------------------------

               CALL CC_TRRETA(ISYMA,LABELA,WORK(KOFFR),WORK(KETA),
     &                        WORK(KEND2),LWRK2,MODEL)

C              Calculate contribution to right transition moment:
C              [F*tA(-wf)]*RE for all excited states of matching symmetry,
C              if requested.
C              -----------------------------------------------------------

               IF ((.NOT.CIS) .AND. (.NOT.LRS2N1)) THEN
                  CALL CC_TRRFTA(ISYMA,LABELA,WORK(KOFFR),
     &                           WORK(KEND1),LWRK1,MODEL)
               END IF

C              Calculate ksiA vector.
C              ----------------------

               KKSI = KETA
               CALL CC_XKSI(WORK(KKSI),LABELA,ISYMA,0,DUMMY,
     &                      WORK(KEND2),LWRK2)

C              Calculate contribution to right transition moment:
C              Mf(wf)*ksiA for all excited states of matching symmetry,
C              if requested.
C              --------------------------------------------------------

               IF ((.NOT.CIS) .AND. LRS2N1) THEN
                  CALL CC_TRRKSI(ISYMA,LABELA,WORK(KOFFR),WORK(KKSI),
     &                           WORK(KEND2),LWRK2,MODEL)
               END IF

C              Calculate left transition moment:
C              LE*ksiA for all excited states of matching symmetry.
C              ----------------------------------------------------

               CALL CC_TRLKSI(ISYMA,LABELA,WORK(KOFFL),WORK(KKSI),
     &                        WORK(KEND2),LWRK2,MODEL)

            END IF

         END IF

      END DO

C     Print right transition moments.
C     -------------------------------

      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     &  'Right transition moments in atomic units:'
      WRITE(LUPRI,'(1X,A,/)')
     &  '-----------------------------------------------'
 
      DO IPRLBL = 1,NPRLBL

         LABELA = PRPLBL(IPRLBL)

         IAB   = 1
         IOPER = ILRES(LABELA,'A')
         IF (IOPER .LE. 0) THEN
            IAB   = 2
            IOPER = ILRES(LABELA,'B')
         END IF

         IF (IOPER .GT. 0) THEN

            IF (IAB .EQ. 1) THEN
               ISYMA = ISYOPR(IALRSOP(IOPER))
            ELSE IF (IAB .EQ. 2) THEN
               ISYMA = ISYOPR(IBLRSOP(IOPER))
            ELSE
               WRITE(LUPRI,*) SECNAM,': lllegal IAB: ',IAB
               CALL QUIT('Internal error in '//SECNAM//' [IAB 2]')
            END IF

            IF (NLOCS(ISYMA) .GT. 0) THEN
               DO IRSD = 1,NXLRSST
                  ISTATE = ILRSST(IRSD)
                  ISYME  = ISYEXC(ISTATE)
                  ISTSY  = ISTATE - ISYOFE(ISYME)
                  EIGV   = EIGVAL(ISTATE)
                  IF (ISYME .EQ. ISYMA) THEN
                    KOFF = KRIGHT + NXLRSST*(IPRLBL - 1) + IRSD - 1
                    WRITE(LUPRI,'(1X,I2,F15.6,2X,A1,A8,A6,1X,F15.8)')
     &              ISTATE,EIGV,'<',LABELA,'|f> = ',WORK(KOFF)
                  END IF
               END DO
            END IF

         END IF

      END DO

C     Print left transition moments.
C     ------------------------------

      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     &  'Left  transition moments in atomic units:'
      WRITE(LUPRI,'(1X,A,/)')
     &  '-----------------------------------------------'

      DO IPRLBL = 1,NPRLBL

         LABELA = PRPLBL(IPRLBL)

         IAB   = 1
         IOPER = ILRES(LABELA,'A')
         IF (IOPER .LE. 0) THEN
            IAB   = 2
            IOPER = ILRES(LABELA,'B')
         END IF

         IF (IOPER .GT. 0) THEN

            IF (IAB .EQ. 1) THEN
               ISYMA = ISYOPR(IALRSOP(IOPER))
            ELSE IF (IAB .EQ. 2) THEN
               ISYMA = ISYOPR(IBLRSOP(IOPER))
            ELSE
               WRITE(LUPRI,*) SECNAM,': lllegal IAB: ',IAB
               CALL QUIT('Internal error in '//SECNAM//' [IAB 3]')
            END IF

            IF (NLOCS(ISYMA) .GT. 0) THEN
               DO IRSD = 1,NXLRSST
                  ISTATE = ILRSST(IRSD)
                  ISYME  = ISYEXC(ISTATE)
                  ISTSY  = ISTATE - ISYOFE(ISYME)
                  EIGV   = EIGVAL(ISTATE)
                  IF (ISYME .EQ. ISYMA) THEN
                    KOFF = KLEFT + NXLRSST*(IPRLBL - 1) + IRSD - 1
                    WRITE(LUPRI,'(1X,I2,F15.6,2X,A3,A8,A4,1X,F15.8)')
     &              ISTATE,EIGV,'<f|',LABELA,'> = ',WORK(KOFF)
                  END IF
               END DO
            END IF

         END IF

      END DO

      CALL FLSHFO(LUPRI)

C     Allocation 3.
C     -------------

      IF (OSCSTR) THEN
         LOSCIL = NEXCI*3*3
      ELSE
         LOSCIL = 0
      END IF

      IF (VELSTR) THEN
         LOSCIV = NEXCI*3*3
      ELSE
         LOSCIV = 0
      END IF

      IF (MIXSTR) THEN
         LOSCIM = NEXCI*3*3
      ELSE
         LOSCIM = 0
      END IF

      IF (ROTLEN) THEN
         LROTL = NEXCI*3
         LCHKL = NEXCI
      ELSE
         LROTL = 0
         LCHKL = 0
      ENDIF

      IF (ROTVEL) THEN
         LROTV = NEXCI*3
         LCHKV = NEXCI
      ELSE
         LROTV = 0
         LCHKV = 0
      ENDIF

      IF (RTNLEN) THEN
         LRQL = NEXCI*3*9
         LRML = NEXCI*3*3
         NWRL = 0
      ELSE
         LRQL = 0
         LRML = 0
      ENDIF

      IF (RTNVEL) THEN
         LRQV = NEXCI*3*9
         LRMV = NEXCI*3*3
         NWRV = 0
      ELSE
         LRQV = 0
         LRMV = 0
      ENDIF

      KOSCS2 = KEND1
      KTRS   = KOSCS2  + LOSCIL
      KVELST = KTRS    + LOSCIL
      KVELST2= KVELST  + LOSCIV
      KMIXST = KVELST2 + LOSCIV
      KMIXST2= KMIXST  + LOSCIM
      KROTL  = KMIXST2 + LOSCIM
      KROTV  = KROTL   + LROTL
      KRQL   = KROTV   + LROTV
      KRML   = KRQL    + LRQL
      KRQL2  = KRML    + LRML
      KRML2  = KRQL2   + LRML
      KRQV   = KRML2   + LRML
      KRMV   = KRQV    + LRQV
      KRQV2  = KRMV    + LRMV
      KRMV2  = KRQV2   + LRMV
      KCHKL  = KRMV2   + LRMV
      KCHKV  = KCHKL   + LCHKL
      KEND3  = KCHKV   + LCHKV
      LWRK3  = LWORK   - KEND3

      IF (LWRK3 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [3]')
      END IF

      IF (OSCSTR) THEN
         CALL DZERO(WORK(KOSCS2),LOSCIL)
         CALL DZERO(WORK(KTRS),LOSCIL)
      END IF
      IF (VELSTR) THEN
         CALL DZERO(WORK(KVELST),LOSCIV)
         CALL DZERO(WORK(KVELST2),LOSCIV)
      END IF
      IF (MIXSTR) THEN
         CALL DZERO(WORK(KMIXST),LOSCIM)
         CALL DZERO(WORK(KMIXST2),LOSCIM)
      END IF
      IF (ROTLEN) THEN
         CALL DZERO(WORK(KROTL),LROTL)
         CALL DZERO(WORK(KROTL),LROTL)
         CALL DZERO(WORK(KCHKL),LCHKL)
      END IF
      IF (ROTVEL) THEN
         CALL DZERO(WORK(KROTV),LROTV)
         CALL DZERO(WORK(KROTV),LROTV)
         CALL DZERO(WORK(KCHKV),LCHKV)
      END IF
      IF (RTNLEN) THEN
         CALL DZERO(WORK(KRQL),LRQL)
         CALL DZERO(WORK(KRML),LRML)
         CALL DZERO(WORK(KRQL2),LRML)
         CALL DZERO(WORK(KRML2),LRML)
      END IF
      IF (RTNVEL) THEN
         CALL DZERO(WORK(KRQV),LRQV)
         CALL DZERO(WORK(KRMV),LRMV)
         CALL DZERO(WORK(KRQV2),LRMV)
         CALL DZERO(WORK(KRMV2),LRMV)
      END IF

C     Calculate linear response residues from transition moments,
C     incl. symmetrization.
C     -----------------------------------------------------------

      WRITE(LUPRI,'(//,1X,A6,A)') MODELP(1:6),
     &  'linear response residue property in atomic units:'
      WRITE(LUPRI,'(1X,A,/)')
     &  '-------------------------------------------------------'

      DO IOPER = 1,NLRSOP

         ISYMA  = ISYOPR(IALRSOP(IOPER))
         ISYMB  = ISYOPR(IBLRSOP(IOPER))
         LABELA = LBLOPR(IALRSOP(IOPER))
         LABELB = LBLOPR(IBLRSOP(IOPER))

         IPROPA = INDPRPCC(LABELA)
         IPROPB = INDPRPCC(LABELB)

         DO IRSD = 1,NXLRSST

            ISTATE = ILRSST(IRSD)
            ISYME  = ISYEXC(ISTATE)
            ISTSY  = ISTATE - ISYOFE(ISYME)
            EIGV   = EIGVAL(ISTATE)
            ISYMEA = MULD2H(ISYME,ISYMA)

            IF ((ISYME.EQ.ISYMA) .AND. (ISYME.EQ.ISYMB)) THEN

               NTOT = NTOT + 1

               KA = NXLRSST*(IPROPA - 1) + IRSD
               KB = NXLRSST*(IPROPB - 1) + IRSD

               K1 = KRIGHT + KA - 1
               K2 = KLEFT  + KB - 1
               K3 = KRIGHT + KB - 1
               K4 = KLEFT  + KA - 1

               IHERMA = ISYMAT(IALRSOP(IOPER))
               IHERMB = ISYMAT(IBLRSOP(IOPER))
               ISASB  = IHERMA*IHERMB

               IF (ISASB .EQ. 0) THEN
                  WRITE(LUPRI,*) ' WARNING: operators ',LABELA,LABELB,
     &                           ' have undefined hermiticities: ',
     &                             IHERMA,IHERMB
                  WRITE(LUPRI,*)
     &            ' Residue not appropriately symmetrized..'
                  SIGN = 1.0D0
               ELSE
                  SIGN = DBLE(ISASB)
               ENDIF
               RESIDAB = WORK(K1)*WORK(K2)
               RESIDBA = WORK(K3)*WORK(K4)
               RESIDUE = 0.5D0*(RESIDAB + SIGN*RESIDBA)
               IF (RESIDUE.GE.0.0D0) THEN
                 SQRRES=SQRT(RESIDUE)
               ELSE
                 SQRRES=-SQRT(-RESIDUE)
               ENDIF
               WRITE(LUPRI,'(1X,A6,A8,A1,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     &         'RES{<<',LABELA,',',LABELB,'>>(',EIGV,')} =',
     &         RESIDUE,' ( ',SQRRES,')'
               IF (BOTHLRS) THEN
                  WRITE(LUPRI,'(1X,A,F12.8,A)')
     &            '  (Unsymmetrized residue: ',RESIDAB,')'
               END IF

               IF (OSCSTR) THEN ! Length gauge oscillator strength
                  IF (LABELA(2:7).EQ.'DIPLEN' .AND.
     &                LABELB(2:7).EQ.'DIPLEN') THEN
                     IF (LABELA(1:5).EQ.'XDIPL') IADR1 = 1
                     IF (LABELA(1:5).EQ.'YDIPL') IADR1 = 2
                     IF (LABELA(1:5).EQ.'ZDIPL') IADR1 = 3
                     IF (LABELB(1:5).EQ.'XDIPL') IADR2 = 1
                     IF (LABELB(1:5).EQ.'YDIPL') IADR2 = 2
                     IF (LABELB(1:5).EQ.'ZDIPL') IADR2 = 3
                     IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KOSCS2-1
                     WORK(IOSCS2) = RESIDUE
                  END IF
               END IF
               IF (VELSTR) THEN ! Velocity gauge oscillator strength
                  IF (LABELA(2:7).EQ.'DIPVEL' .AND.
     &                LABELB(2:7).EQ.'DIPVEL') THEN
                     IF (LABELA(1:5).EQ.'XDIPV') IADR1 = 1
                     IF (LABELA(1:5).EQ.'YDIPV') IADR1 = 2
                     IF (LABELA(1:5).EQ.'ZDIPV') IADR1 = 3
                     IF (LABELB(1:5).EQ.'XDIPV') IADR2 = 1
                     IF (LABELB(1:5).EQ.'YDIPV') IADR2 = 2
                     IF (LABELB(1:5).EQ.'ZDIPV') IADR2 = 3
                     IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KVELST-1
                     WORK(IOSCS2) = RESIDUE
                  END IF
               END IF
               IF (MIXSTR) THEN ! Mixed gauge oscillator strength
                  IF (LABELA(2:7).EQ.'DIPLEN' .AND.
     &                LABELB(2:7).EQ.'DIPVEL') THEN
                     IF (LABELA(1:5).EQ.'XDIPL') IADR1 = 1
                     IF (LABELA(1:5).EQ.'YDIPL') IADR1 = 2
                     IF (LABELA(1:5).EQ.'ZDIPL') IADR1 = 3
                     IF (LABELB(1:5).EQ.'XDIPV') IADR2 = 1
                     IF (LABELB(1:5).EQ.'YDIPV') IADR2 = 2
                     IF (LABELB(1:5).EQ.'ZDIPV') IADR2 = 3
                     IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KMIXST-1
                     WORK(IOSCS2) = RESIDUE
                  END IF
               END IF
               IF (ROTLEN) THEN ! Length gauge rotatory strength
                  IF (LABELA(2:7) .EQ. 'DIPLEN' .AND.
     &                LABELB(2:7) .EQ. 'ANGMOM') THEN
                     IF (LABELA(1:5).EQ.'XDIPL') IADR1 = 1
                     IF (LABELA(1:5).EQ.'YDIPL') IADR1 = 2
                     IF (LABELA(1:5).EQ.'ZDIPL') IADR1 = 3
                     IF (LABELB(1:5).EQ.'XANGM') IADR2 = 1
                     IF (LABELB(1:5).EQ.'YANGM') IADR2 = 2
                     IF (LABELB(1:5).EQ.'ZANGM') IADR2 = 3
                     IF (IADR1 .EQ. IADR2) THEN
                        IROTST = KROTL + 3*(ISTATE-1) + IADR1 - 1
                        WORK(IROTST) = RESIDUE
                     END IF
                  END IF
               END IF
               IF (ROTVEL) THEN ! Velocity gauge rotatory strength
                  IF (LABELA(2:7) .EQ. 'DIPVEL' .AND.
     &                LABELB(2:7) .EQ. 'ANGMOM') THEN
                     IF (LABELA(1:5).EQ.'XDIPV') IADR1 = 1
                     IF (LABELA(1:5).EQ.'YDIPV') IADR1 = 2
                     IF (LABELA(1:5).EQ.'ZDIPV') IADR1 = 3
                     IF (LABELB(1:5).EQ.'XANGM') IADR2 = 1
                     IF (LABELB(1:5).EQ.'YANGM') IADR2 = 2
                     IF (LABELB(1:5).EQ.'ZANGM') IADR2 = 3
                     IF (IADR1 .EQ. IADR2) THEN
                        IROTST = KROTV + 3*(ISTATE-1) + IADR1 - 1
                        WORK(IROTST) = RESIDUE
                     END IF
                  END IF
               END IF
               IF (RTNLEN) THEN
                  IF (LABELA(2:7) .EQ. 'DIPLEN') THEN
                     IF (LABELB(3:8) .EQ. 'SECMOM') THEN
                        IADR1 = -999999
                        IF (LABELA(1:5).EQ.'XDIPL') IADR1 = 1
                        IF (LABELA(1:5).EQ.'YDIPL') IADR1 = 2
                        IF (LABELA(1:5).EQ.'ZDIPL') IADR1 = 3
                        IADR23 = -999999
                        IADR32 = -999999
                        IF (LABELB(1:5).EQ.'XXSEC') THEN
                           IADR23 = 1
                           IADR32 = 0
                        ELSE IF (LABELB(1:5).EQ.'XYSEC') THEN
                           IADR23 = 4
                           IADR32 = 2
                        ELSE IF (LABELB(1:5).EQ.'XZSEC') THEN
                           IADR23 = 7
                           IADR32 = 3
                        ELSE IF (LABELB(1:5).EQ.'YYSEC') THEN
                           IADR23 = 5
                           IADR32 = 0
                        ELSE IF (LABELB(1:5).EQ.'YZSEC') THEN
                           IADR23 = 8
                           IADR32 = 6
                        ELSE IF (LABELB(1:5).EQ.'ZZSEC') THEN
                           IADR23 = 9
                           IADR32 = 0
                        END IF
                        IF ((IADR1.LT.0) .OR. (IADR23.LT.0) .OR.
     &                      (IADR32.LT.0)) THEN
                           CALL QUIT('RQL error in '//SECNAM)
                        END IF
                        IRTEN = KRQL + 3*9*(ISTATE-1)
     &                        + 3*(IADR23-1) + IADR1 - 1
                        WORK(IRTEN) = RESIDUE
                        IF (IADR32 .NE. 0) THEN
                           IRTEN = KRQL + 3*9*(ISTATE-1)
     &                           + 3*(IADR32-1) + IADR1 - 1
                           WORK(IRTEN) = RESIDUE
                        END IF
                     ELSE IF (LABELB(2:7) .EQ. 'ANGMOM') THEN
                        IF (LABELA(1:5).EQ.'XDIPL') IADR1 = 1
                        IF (LABELA(1:5).EQ.'YDIPL') IADR1 = 2
                        IF (LABELA(1:5).EQ.'ZDIPL') IADR1 = 3
                        IF (LABELB(1:5).EQ.'XANGM') IADR2 = 1
                        IF (LABELB(1:5).EQ.'YANGM') IADR2 = 2
                        IF (LABELB(1:5).EQ.'ZANGM') IADR2 = 3
                        IRTEN = KRML + 3*3*(ISTATE-1)
     &                        + 3*(IADR2-1) + IADR1 - 1
                        WORK(IRTEN) = RESIDUE
                     END IF
                  END IF
               END IF
               IF (RTNVEL) THEN
                  IF (LABELA(2:7) .EQ. 'DIPVEL') THEN
                     IF (LABELB(3:8) .EQ. 'ROTSTR') THEN
                        IF (LABELA(1:5).EQ.'XDIPV') IADR1 = 1
                        IF (LABELA(1:5).EQ.'YDIPV') IADR1 = 2
                        IF (LABELA(1:5).EQ.'ZDIPV') IADR1 = 3
                        IF (LABELB(1:5).EQ.'XXROT') THEN
                           IADR23 = 1
                           IADR32 = 0
                        ELSE IF (LABELB(1:5).EQ.'XYROT') THEN
                           IADR23 = 4
                           IADR32 = 2
                        ELSE IF (LABELB(1:5).EQ.'XZROT') THEN
                           IADR23 = 7
                           IADR32 = 3
                        ELSE IF (LABELB(1:5).EQ.'YYROT') THEN
                           IADR23 = 5
                           IADR32 = 0
                        ELSE IF (LABELB(1:5).EQ.'YZROT') THEN
                           IADR23 = 8
                           IADR32 = 6
                        ELSE IF (LABELB(1:5).EQ.'ZZROT') THEN
                           IADR23 = 9
                           IADR32 = 0
                        END IF
                        IRTEN = KRQV + 3*9*(ISTATE-1)
     &                        + 3*(IADR23-1) + IADR1 - 1
                        WORK(IRTEN) = RESIDUE
                        IF (IADR32 .NE. 0) THEN
                           IRTEN = KRQV + 3*9*(ISTATE-1)
     &                           + 3*(IADR32-1) + IADR1 - 1
                           WORK(IRTEN) = RESIDUE
                        END IF
                     ELSE IF (LABELB(2:7) .EQ. 'ANGMOM') THEN
                        IF (LABELA(1:5).EQ.'XDIPV') IADR1 = 1
                        IF (LABELA(1:5).EQ.'YDIPV') IADR1 = 2
                        IF (LABELA(1:5).EQ.'ZDIPV') IADR1 = 3
                        IF (LABELB(1:5).EQ.'XANGM') IADR2 = 1
                        IF (LABELB(1:5).EQ.'YANGM') IADR2 = 2
                        IF (LABELB(1:5).EQ.'ZANGM') IADR2 = 3
                        IRTEN = KRMV + 3*3*(ISTATE-1)
     &                        + 3*(IADR2-1) + IADR1 - 1
                        WORK(IRTEN) = RESIDUE
                     END IF
                  END IF
               END IF

            ELSE

               RESIDUE = 0.0D0
               SQRRES  = 0.0D0

            END IF

            IF (LABELA .EQ. LABELB) THEN
               CALL WRIPRO(SQRRES,MODEL,-1,
     &                     LABELA,LABELB,LABELA,LABELB,
     &                     EIGV,EIGV,EIGV,ISYMEA,ISYME,1,ISTATE)
               OSCCON = EIGV*SQRRES*SQRRES 
               CALL WRIPRO(OSCCON,MODEL,-21,
     &                     LABELA,LABELB,LABELA,LABELB,
     &                     EIGV,EIGV,EIGV,ISYMEA,ISYME,1,ISTATE)
            END IF

         END DO

      END DO

      CALL FLSHFO(LUPRI)

C     Print summary on unit LURES.
C     ----------------------------

      LUOSC = LURES

      IF (OSCSTR) CALL DCOPY(LOSCIL,WORK(KOSCS2),1,WORK(KTRS),1)
      IF (VELSTR) CALL DCOPY(LOSCIV,WORK(KVELST),1,WORK(KVELST2),1)
      IF (MIXSTR) CALL DCOPY(LOSCIM,WORK(KMIXST),1,WORK(KMIXST2),1)

      IF (OSCSTR .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Osci'
     &     //'llator Strength       |'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        | Dipole Strength(a.u.) | Oscillator stre'
     &    //'ngth  | Direction   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM  = 1, NSYM
            DO IEX   = 1, NCCEXCI(ISYM,1)
               ISTATE = ISYOFE(ISYM) + IEX
               EIGV   = EIGVAL(ISTATE)
               KOSCSI = KOSCS2 + 3*3*(ISTATE-1)
               KTRSI  = KTRS   + 3*3*(ISTATE-1)
               LCALC  = .FALSE.
               LDIP   = 1
               DO IRSD  = 1, NXLRSST
                 ISTATE = ILRSST(IRSD)
                 ISYME  = ISYEXC(ISTATE)
                 ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
                 IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
               END DO
               CALL CC_OSCPRI(WORK(KTRSI),WORK(KOSCSI),EIGV,
     &                        IEX,ISYM,WORK(KEND2),LEND2,MODELP,LCALC,
     &                        LDIP,LUOSC)
            END DO
            IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
               NREST = 0
               DO ISYM2 = ISYM+1,NSYM
                  NREST = NREST + NCCEXCI(ISYM2,1)
               END DO
               IF (NREST.EQ.0) GOTO 9001
               WRITE(LUOSC,'(A)')
     &         ' +----------------------------------------------'
     &        //'-------------------------------+'
            END IF
 9001       CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

      END IF

      IF (VELSTR .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Osci'
     &     //'llator Strength       |'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        | Dipole Strength(a.u.) | Oscillator stre'
     &    //'ngth  | Direction   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM  = 1, NSYM
            DO IEX   = 1, NCCEXCI(ISYM,1)
               ISTATE = ISYOFE(ISYM) + IEX
               EIGV   = EIGVAL(ISTATE)
               KOSCSI = KVELST + 3*3*(ISTATE-1)
               KTRSI  = KVELST2+ 3*3*(ISTATE-1)
               LCALC  = .FALSE.
               LDIP   = 2
               DO IRSD  = 1, NXLRSST
                  ISTATE = ILRSST(IRSD)
                  ISYME  = ISYEXC(ISTATE)
                  ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
                  IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
               END DO
               CALL CC_OSCPRI(WORK(KTRSI),WORK(KOSCSI),EIGV,
     &                        IEX,ISYM,WORK(KEND2),LEND2,MODELP,LCALC,
     &                        LDIP,LUOSC)
            END DO
            IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
               NREST = 0
               DO ISYM2 = ISYM+1,NSYM
                  NREST = NREST + NCCEXCI(ISYM2,1)
               END DO
               IF (NREST.EQ.0) GOTO 9005
               WRITE(LUOSC,'(A)')
     &         ' +----------------------------------------------'
     &         //'-------------------------------+'
            END IF
 9005       CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

      END IF

      IF (MIXSTR .AND. (CCS.OR.CC2.OR.CCSD)) THEN
 
         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Mixed    Gauge Osci'
     &     //'llator Strength       |'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        | Dipole Strength(a.u.) | Oscillator stre'
     &    //'ngth  | Direction   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'
 
         DO ISYM  = 1, NSYM
            DO IEX   = 1, NCCEXCI(ISYM,1)
               ISTATE = ISYOFE(ISYM) + IEX
               EIGV   = EIGVAL(ISTATE)
               KOSCSI = KMIXST + 3*3*(ISTATE-1)
               KTRSI  = KMIXST2+ 3*3*(ISTATE-1)
               LCALC  = .FALSE.
               LDIP   = 3
               DO IRSD  = 1, NXLRSST
                  ISTATE = ILRSST(IRSD)
                  ISYME  = ISYEXC(ISTATE)
                  ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
                  IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
               END DO
               CALL CC_OSCPRI(WORK(KTRSI),WORK(KOSCSI),EIGV,
     &                        IEX,ISYM,WORK(KEND2),LEND2,MODELP,LCALC,
     &                        LDIP,LUOSC)
            END DO
            IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
               NREST = 0
               DO ISYM2 = ISYM+1,NSYM
                  NREST = NREST + NCCEXCI(ISYM2,1)
               END DO
               IF (NREST.EQ.0) GOTO 9008
               WRITE(LUOSC,'(A)')
     &         ' +----------------------------------------------'
     &        //'-------------------------------+'
            END IF
 9008       CONTINUE
         END DO
 
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'
 
      END IF

      IF (ROTLEN .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rota'
     &     //'tory Strength         |'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Direction   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KTRSI  = KROTL + 3*(ISTATE-1)
           KSTREN = KCHKL + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_ROTPRI(WORK(KTRSI),WORK(KSTREN),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9009
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9009     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

      END IF

      IF (ROTVEL .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rota'
     &     //'tory Strength         |'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Direction   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KTRSI  = KROTV + 3*(ISTATE-1)
           KSTREN = KCHKV + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_ROTPRI(WORK(KTRSI),WORK(KSTREN),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9010
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9010     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

      END IF

      IF (RTNLEN .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rot.'
     &     //'Str. Tensor, El. Quad.|'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFQ  = KRQL  + 3*9*(ISTATE-1)
           KOFQ2  = KRQL2 + 3*3*(ISTATE-1)
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTQPRI(WORK(KOFFQ),WORK(KOFQ2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,NWRL)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9011
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9011     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rot.'
     &     //'Str. Tensor, Mag. Dip.|'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFM  = KRML  + 3*3*(ISTATE-1)
           KOFM2  = KRML2 + 3*3*(ISTATE-1)
           KSTREN = KCHKL + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTMPRI(WORK(KOFFM),WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRL)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9012
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9012     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         CALL DAXPY(LRML,1.0D0,WORK(KRQL2),1,WORK(KRML2),1)  ! Get total tensor (in KRML2)

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Length   Gauge Rot.'
     &     //'Str. Tensor, Total    |'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFM2  = KRML2 + 3*3*(ISTATE-1)
           KSTREN = KCHKL + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 1
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTTPRI(WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRL)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9013
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9013     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         IF (NWRL .NE. 0) THEN
            WRITE(LUOSC,'(//,1X,A,I4,A)')
     &      '***NOTICE:',NWRL,' warnings issued for Rot. Str. Tensors.'
            WRITE(LUOSC,'(1X,A)')
     &      '           Length gauge tensors are wrong!'
         END IF

      END IF

      IF (RTNVEL .AND. (CCS.OR.CC2.OR.CCSD)) THEN

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rot.'
     &     //'Str. Tensor, El. Quad.|'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFQ  = KRQV  + 3*9*(ISTATE-1)
           KOFQ2  = KRQV2 + 3*3*(ISTATE-1)
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTQPRI(WORK(KOFFQ),WORK(KOFQ2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,NWRV)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9014
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9014     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rot.'
     &     //'Str. Tensor, Mag. Dip.|'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFFM  = KRMV  + 3*3*(ISTATE-1)
           KOFM2  = KRMV2 + 3*3*(ISTATE-1)
           KSTREN = KCHKV + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTMPRI(WORK(KOFFM),WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRV)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9015
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9015     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         CALL DAXPY(LRMV,1.0D0,WORK(KRQV2),1,WORK(KRMV2),1)  ! Get total tensor (in KRMV2)

         WRITE(LUOSC,'(//A)')
     &     ' +=============================================='
     &    //'===============================+'
         WRITE(LUOSC,'(1X,A26,A10,A)')
     &     '|  sym. | Exci.  |        ',MODELP,' Velocity Gauge Rot.'
     &     //'Str. Tensor, Total    |'
         WRITE(LUOSC,'(A)')
     &     ' |(spin, |        +-----------------------------'
     &    //'-------------------------------+'
         WRITE(LUOSC,'(1X,A)')
     &     '| spat) |        |        D-55 SI        |      D-40 cgs  '
     &    //'      | Component   |'
         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         DO ISYM = 1, NSYM
          DO IEX = 1, NCCEXCI(ISYM,1)
           ISTATE = ISYOFE(ISYM) + IEX
           EIGV   = EIGVAL(ISTATE)
           KOFM2  = KRMV2 + 3*3*(ISTATE-1)
           KSTREN = KCHKV + ISTATE - 1
           LCALC  = .FALSE.
           LDIP   = 2
           DO IRSD  = 1, NXLRSST
             ISTATE = ILRSST(IRSD)
             ISYME  = ISYEXC(ISTATE)
             ISTSY  = ILRSST(IRSD)-ISYOFE(ISYME)
             IF ((ISYME.EQ.ISYM).AND.(IEX.EQ.ISTSY)) LCALC =.TRUE.
           END DO
           CALL CC_RTTPRI(WORK(KOFM2),EIGV,IEX,ISYM,MODELP,
     &                    LCALC,LDIP,LUOSC,WORK(KSTREN),NWRV)

          END DO

          IF (.NOT.((ISYM .EQ. NSYM).OR.(NCCEXCI(ISYM,1).EQ.0))) THEN
             NREST = 0
             DO ISYM2 = ISYM+1,NSYM
                NREST = NREST + NCCEXCI(ISYM2,1)
             END DO
             IF (NREST.EQ.0) GOTO 9016
             WRITE(LUOSC,'(A)')
     &       ' +----------------------------------------------'
     &      //'-------------------------------+'
          END IF
 9016     CONTINUE
         END DO

         WRITE(LUOSC,'(A)')
     &     ' +=============================================='
     &    //'===============================+'

         IF (NWRV .NE. 0) THEN
            WRITE(LUOSC,'(//,1X,A,I4,A)')
     &      '***NOTICE:',NWRV,' warnings issued for Rot. Str. Tensors.'
            WRITE(LUOSC,'(1X,A)')
     &      '           Velocity gauge tensors are wrong!'
         END IF

      END IF

      IF (ROTLEN .OR. ROTVEL .OR. RTNLEN .OR. RTNVEL) THEN
         WRITE(LUOSC,'(/,1X,A)')
     &   'Conversion factors for rotatory strengths:'
         WRITE(LUOSC,'(3X,A,F15.10,A)')
     &   'SI  units:   1 a.u. = ',RAUSI,'D-55 A^2 m^3 s'
         WRITE(LUOSC,'(3X,A,F15.10,A)')
     &   'cgs units:   1 a.u. = ',RAUCGS,'D-40 cm^5 g s^-2'
      END IF

C     Print timings and exit.
C     -----------------------

 999  TIMTOT = SECOND() - TIMTOT
      WRITE(LUPRI,'(/,1X,A,I7,A,F10.2,A)')
     & '>>> Total time for',NTOT,' linear response residues: ',
     & TIMTOT,' seconds.'
      CALL FLSHFO(LUPRI)

      CALL QEXIT(SECNAM)

      RETURN
      END
C  /* Deck ilres */
      INTEGER FUNCTION ILRES(LABEL,LIST)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Find the first index of the operator LABEL on the residue
C              list indicated by LIST ('A' or 'B') for which the total
C              residue is symmetry--allowed.
C              If LABEL is not on the list, ILRES = -1
C              If LIST is illegal,          ILRES = -2.
C
#include <implicit.h>
      CHARACTER*8 LABEL
      CHARACTER*1 LIST
#include <cclres.h>
#include <ccroper.h>

      CHARACTER*8 LOCLAB

      ILRES = -1

      IF (LIST .EQ. 'A') THEN

         DO IOPER = 1,NLRSOP
            LOCLAB = LBLOPR(IALRSOP(IOPER))
            IF (LABEL(1:8) .EQ. LOCLAB(1:8)) THEN
               ISYMA = ISYOPR(IALRSOP(IOPER))
               ISYMB = ISYOPR(IBLRSOP(IOPER))
               IF (ISYMA .EQ. ISYMB) THEN
                  ILRES = IOPER
                  RETURN
               END IF
            END IF
         END DO

      ELSE IF (LIST .EQ. 'B') THEN

         DO IOPER = 1,NLRSOP
            LOCLAB = LBLOPR(IBLRSOP(IOPER))
            IF (LABEL(1:8) .EQ. LOCLAB(1:8)) THEN
               ISYMA = ISYOPR(IALRSOP(IOPER))
               ISYMB = ISYOPR(IBLRSOP(IOPER))
               IF (ISYMA .EQ. ISYMB) THEN
                  ILRES = IOPER
                  RETURN
               END IF
            END IF
         END DO

      ELSE

         ILRES = -2

      END IF

      RETURN
      END
C  /* Deck cc_trreta */
      SUBROUTINE CC_TRRETA(ISYMA,LABELA,TRRMOM,ETA,WORK,LWORK,MODEL)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Calculate etaA*RE contributions to right ground-excited
C              state transition moments for all excited states of matching
C              symmetry.
C
#include <implicit.h>
      DIMENSION TRRMOM(*), ETA(*), WORK(LWORK)
      CHARACTER*8  LABELA
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <cclres.h>
#include <ccexci.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_TRRETA')

C     Allocation.
C     -----------

      NTAMP = NT1AM(ISYMA)
      IF (.NOT. CCS) NTAMP = NTAMP + NT2AM(ISYMA)

      KRE1 = 1
      KRE2 = KRE1 + NT1AM(ISYMA)
      IF (CCS) THEN
         KEND = KRE2
      ELSE
         KEND = KRE2 + NT2AM(ISYMA)
      END IF
      LWRK = LWORK - KEND + 1

      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM)
      END IF

C     Loop over requested excited states.
C     -----------------------------------

      DO IRSD = 1,NXLRSST

         ISTATE = ILRSST(IRSD)
         ISYME  = ISYEXC(ISTATE)

         IF (ISYME .EQ. ISYMA) THEN

C           Calculate contribution.
C           -----------------------

            IOPT = 3
            CALL CC_RDRSP('RE',ISTATE,ISYMA,IOPT,MODEL,WORK(KRE1),
     &                    WORK(KRE2))

            CONTR = DDOT(NTAMP,ETA,1,WORK(KRE1),1)

            TRRMOM(IRSD) = TRRMOM(IRSD) + CONTR

            IF (IPRINT .GT. 2) THEN
               ISTSY = ISTATE - ISYOFE(ISYME)
               WRITE(LUPRI,'(1X,A1,A8,A3,A,F12.6,A,I3,A,I2,A)')
     &         '<',LABELA,'|f>',' EtaA*RE cont. = ',CONTR,
     &         '  (f:',ISTSY,' of sym.',ISYME,')'
            END IF

         END IF

      END DO

      RETURN
      END
C  /* Deck cc_trlksi */
      SUBROUTINE CC_TRLKSI(ISYMA,LABELA,TRLMOM,XKSI,WORK,LWORK,MODEL)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Calculate LE*ksiA contributions to left ground-excited
C              state transition moments for all excited states of matching
C              symmetry.
C
#include <implicit.h>
      DIMENSION TRLMOM(*), XKSI(*), WORK(LWORK)
      CHARACTER*8  LABELA
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <cclres.h>
#include <ccexci.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_TRLKSI')

C     Allocation.
C     -----------

      NTAMP = NT1AM(ISYMA)
      IF (.NOT. CCS) NTAMP = NTAMP + NT2AM(ISYMA)

      KLE1 = 1
      KLE2 = KLE1 + NT1AM(ISYMA)
      IF (CCS) THEN
         KEND = KLE2
      ELSE
         KEND = KLE2 + NT2AM(ISYMA)
      END IF
      LWRK = LWORK - KEND + 1

      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM)
      END IF

C     Loop over requested excited states.
C     -----------------------------------

      DO IRSD = 1,NXLRSST

         ISTATE = ILRSST(IRSD)
         ISYME  = ISYEXC(ISTATE)

         IF (ISYME .EQ. ISYMA) THEN

C           Calculate contribution.
C           -----------------------

            IOPT = 3
            CALL CC_RDRSP('LE',ISTATE,ISYMA,IOPT,MODEL,WORK(KLE1),
     &                    WORK(KLE2))

            CONTR = DDOT(NTAMP,WORK(KLE1),1,XKSI,1)

            TRLMOM(IRSD) = TRLMOM(IRSD) + CONTR

            IF (IPRINT .GT. 2) THEN
               ISTSY = ISTATE - ISYOFE(ISYME)
               WRITE(LUPRI,'(1X,A3,A8,A1,A,F12.6,A,I3,A,I2,A)')
     &         '<f|',LABELA,'>',' LE*ksiA cont. = ',CONTR,
     &         '  (f:',ISTSY,' of sym.',ISYME,')'
            END IF

         END IF

      END DO

      RETURN
      END
C  /* Deck cc_trrksi */
      SUBROUTINE CC_TRRKSI(ISYMA,LABELA,TRRMOM,XKSI,WORK,LWORK,MODEL)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Calculate Mf*ksiA contributions to right ground-excited
C              state transition moments for all excited states of matching
C              symmetry.
C
#include <implicit.h>
      DIMENSION TRRMOM(*), XKSI(*), WORK(LWORK)
      CHARACTER*8  LABELA
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <cclres.h>
#include <ccexci.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_TRRKSI')

C     Allocation.
C     -----------

      NTAMP = NT1AM(ISYMA)
      IF (.NOT. CCS) NTAMP = NTAMP + NT2AM(ISYMA)

      KMF1 = 1
      KMF2 = KMF1 + NT1AM(ISYMA)
      IF (CCS) THEN
         KEND = KMF2
      ELSE
         KEND = KMF2 + NT2AM(ISYMA)
      END IF
      LWRK = LWORK - KEND + 1

      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM)
      END IF

C     Loop over requested excited states.
C     -----------------------------------

      DO IRSD = 1,NXLRSST

         ISTATE = ILRSST(IRSD)
         ISYME  = ISYEXC(ISTATE)

         IF (ISYME .EQ. ISYMA) THEN

C           Calculate contribution.
C           -----------------------

            IOPT   = 3
            ILSTNR = ILRMAMP(ISTATE,EIGVAL(ISTATE),ISYMA)
            CALL CC_RDRSP('M1',ILSTNR,ISYMA,IOPT,MODEL,WORK(KMF1),
     &                    WORK(KMF2))

            CONTR = DDOT(NTAMP,WORK(KMF1),1,XKSI,1)

            TRRMOM(IRSD) = TRRMOM(IRSD) + CONTR

            IF (IPRINT .GT. 2) THEN
               ISTSY = ISTATE - ISYOFE(ISYME)
               WRITE(LUPRI,'(1X,A1,A8,A3,A,F12.6,A,I3,A,I2,A)')
     &         '<',LABELA,'|f>',' Mf*ksiA cont. = ',CONTR,
     &         '  (f:',ISTSY,' of sym.',ISYME,')'
            END IF

         END IF

      END DO

      RETURN
      END
C  /* Deck cc_trrfta */
      SUBROUTINE CC_TRRFTA(ISYMA,LABELA,TRRMOM,WORK,LWORK,MODEL)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Purpose: Calculate [F*tA(-wf)]*RE contributions to right ground-excited
C              state transition moments for all excited states of matching
C              symmetry.
C
#include <implicit.h>
      DIMENSION TRRMOM(*), WORK(LWORK)
      CHARACTER*8  LABELA
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <cclres.h>
#include <ccexci.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_TRRFTA')

C     Allocation.
C     -----------

      NTAMP = NT1AM(ISYMA)
      IF (.NOT. CCS) NTAMP = NTAMP + NT2AM(ISYMA)

      KFTA = 1
      KRE  = KFTA  + NTAMP
      KEND = KRE   + NTAMP
      LWRK = LWORK - KEND + 1

      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM)
      END IF

      KFTA1 = KFTA
      KFTA2 = KFTA + NT1AM(ISYMA)
      KRE1  = KRE
      KRE2  = KRE  + NT1AM(ISYMA)

C     Loop over requested excited states.
C     -----------------------------------

      DO IRSD = 1,NXLRSST

         ISTATE = ILRSST(IRSD)
         ISYME  = ISYEXC(ISTATE)

         IF (ISYME .EQ. ISYMA) THEN

C           Calculate contribution.
C           -----------------------

            IOPT   = 3
            ILSTNR = IR1TAMP(LABELA,.FALSE.,-EIGVAL(ISTATE),ISYMA)
            CALL CC_RDRSP('F1',ILSTNR,ISYMA,IOPT,MODEL,WORK(KFTA1),
     &                    WORK(KFTA2))
            IOPT   = 3
            CALL CC_RDRSP('RE',ISTATE,ISYMA,IOPT,MODEL,WORK(KRE1),
     &                    WORK(KRE2))

            CONTR = DDOT(NTAMP,WORK(KFTA),1,WORK(KRE),1)

            TRRMOM(IRSD) = TRRMOM(IRSD) + CONTR

            IF (IPRINT .GT. 2) THEN
               ISTSY = ISTATE - ISYOFE(ISYME)
               WRITE(LUPRI,'(1X,A1,A8,A3,A,F12.6,A,I3,A,I2,A)')
     &         '<',LABELA,'|f>',' [F*tA(-wf)]*RE cont. = ',CONTR,
     &         '  (f:',ISTSY,' of sym.',ISYME,')'
            END IF

         END IF

      END DO

      RETURN
      END
c*DECK WRIPRO
       SUBROUTINE WRIPRO(PROP,LABEL,NORD,LABX,LABY,LABZ,LABU,
     *                   FRQY,FRQZ,FRQU,ISYMIN,ISYMEX,ISPINEX,IEX)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Add response property to list of property information to be
C              passed to numerical differentiantion/averaging.
C
C     Ove Christiansen August 1999.
C
C     NORD = 0    energy (ground or excited)
C            1    exp. value
C            2    Linear response function
C            3    Quadratic response function
C            4    Cubic response function
C           -1    ground - excited  transition matrix element, <0|x|i>
C           -2    excited - excited transition matrix element, |<i|x|f>|
C           -11   First order excited state property, <i|x|i>
C           -20   <0|x|i><i|y|0>
C           -21   w*<0|x|i><i|y|0>
C           -22   (w_f - w_i)*|<i|x|f>|^2
C           -30   D_pa
C           -31   D_pe
C           -32   D_pc
C           -33   w1w2D_pa
C           -34   w1w2D_pe
C           -35   w1w2D_pc
C           -400  oscillator strength
C            401  chemical shielding isotropic
C            402  chemical shielding tensor
C-----------------------------------------------------------------------------
C
#include <implicit.h>
#include <maxorb.h>
C
#include <dummy.h>
#include <iratdef.h>
#include <priunit.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <prpc.h>
#include <inftap.h>
C
      LOGICAL EXIST,L1,L2,L3,L4,LI1,LI2
      PARAMETER (TOLFRQ=1.0D-08,ONE=1.0D0,XMONE=-1.0D0,TOLEXCI =1.0D-02)
C
      CHARACTER LABEL*10, LABX*8, LABY*8, LABZ*8, LABU*8
C
C--------------------------------------------------
C
C
      IF (NOEONL .AND. (NORD.EQ.0)) THEN
C         if energy and NOEONList = true then skip addition to list.
        RETURN
      ELSE
        EXIST = .FALSE.
        IF (EXIST) THEN
c          IPRMI = IHIT
        ELSE
           NPRMI = NPRMI + 1
           IPRMI = NPRMI
        ENDIF
C
        WRITE(LUNDPF,
     *   '(I5,I3,I4,1X,A10,E23.16,4(1X,A8),3E23.16,3I4)')
     *   IPRMI,ISYMIN,NORD,LABEL,PROP,
     *   LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYMEX,ISPINEX,IEX
        WRITE(LUNMPF,
     *   '(I5,I3,I4,1X,A10,E23.16,4(1X,A8),3E23.16,3I4)')
     *   IPRMI,ISYMIN,NORD,LABEL,PROP,
     *   LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYMEX,ISPINEX,IEX
      ENDIF
C
      END
      SUBROUTINE stripblanks(tobestripped)
#include <priunit.h>
c
c     mbh: transform nuclei string, e.g. 'C0   1' into
c          'C0_1   '
c
      character tobestripped*8,helper*8,lastchar*1
      integer i,j,idx
c
      idx=1
      helper='        '
      write(lupri,*)'String on input: "',tobestripped,'"'
c
c strip all blanks from beginning
c
      do 10 i=1,8
         if(tobestripped(i:i).ne.' ') goto 15
10    continue
15    continue
c
c add to helper until we hit a blank again
c
      do 20 j=i,8
         if(tobestripped(j:j).ne.' ') then
            helper(idx:idx)=tobestripped(j:j)
         else 
            goto 25
         endif
         idx=idx+1
20    continue
25    continue
      helper(idx:idx)='_'
      lastchar='_'
      idx=idx+1
c
c strip all blanks from here to next 'item'
c
      do 30 i=j,8
         if(tobestripped(i:i).ne.' ') then
            helper(idx:idx)=tobestripped(i:i)
            lastchar=tobestripped(i:i)
            idx=idx+1
         endif
30    continue
      if(lastchar.eq.'_') helper(idx-1:idx-1)=' '
      write(lupri,*)'String on input: "',helper,'"'
      tobestripped=helper 
      end
C--------------------------------------------------------------
