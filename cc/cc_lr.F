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
c*DECK CC_LR    
       SUBROUTINE CC_LR(WORK,LWORK)
C
C----------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              polarizabilities.
C              (without orbital relaxation)
C
C              CIS, CCS, CC2, CCSD, MCC2
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
#include <maxorb.h>
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
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <symmet.h>
#include <alphac.h>
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
         MODEL = 'CC2'
      ENDIF
      IF (MCC2) THEN
         MODEL = 'MCC2'
      ENDIF
      IF (CCS) THEN
         MODEL = 'CCS'
      ENDIF
      IF (CC3  ) THEN
         MODEL = 'CC3'
         WRITE(LUPRI,'(/,1x,A)') 'CC3 Polari not implemented yet'
         CALL QEXIT('CC_LR')
         RETURN
      ENDIF
      IF (CC1A) THEN
         MODEL = 'CCSDT-1a'
         WRITE(LUPRI,'(/,1x,A)') 'CC1A Polari not implemented yet'
         CALL QEXIT('CC_LR')
         RETURN
      ENDIF
      IF (CC1B) THEN
         MODEL = 'CCSDT-1b'
         WRITE(LUPRI,'(/,1x,A)') 'CC1B Polari not implemented yet'
         CALL QEXIT('CC_LR')
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
C     open AOPROPER file for GETGP routine of the RSP program...
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
     *                        
C
            ENDIF
 
C-------------------------------------------------------------------
C             construct the zeroth-order eff. CC Fock matrix in MO
C             and initialize the 'eff.' sec. order connection matrix
C-------------------------------------------------------------------
              IF ( (LRLXA .OR. LPDBSA) .AND. (LRLXB .OR. LPDBSB) ) THEN
                IFOCK = IEFFFOCK('HAM0    ',ISYM0,1)
                IADRF = IADRFCK(1,IFOCK)

                LUFCK = 0
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
                    IOPT = 1
                    IKAPPA = IR1TAMP(LABELB,LRLXB,FREQB,ISYMB)
                    CALL CC_RDHFRSP('R1',IKAPPA,ISYMB,IOPT,MODEL,
     &                              WORK(KAPB))
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
                    IOPT = 1
                    IKAPPA = IR1TAMP(LABELA,LRLXA,FREQA,ISYMA)
                    CALL CC_RDHFRSP('R1',IKAPPA,ISYMA,IOPT,MODEL,
     &                              WORK(KAPA))
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

               IDXR = IR1TAMP(LABELA,LRLXA,+FREQA,ISYMA)
               CALL CC_GETHFGD(IDXR,'R1 ',LRTLBL,IDUM,IDUM,RDUM,
     *                         ISYLRT,FRQLRT,LORXLRT,IDUM,NLRTLBL,
     *                         MAXTLBL,IREAL,WORK(KCMO),WORK(KUDV),
     *                         WORK(KXINDX),FRVAL,WORK(KG1),LWRKG1)

               IOPT = 1
               IDXR = IR1TAMP(LABELB,LRLXB,+FREQB,ISYMB)
               CALL CC_RDHFRSP('R1',IDXR,ISYMB,IOPT,MODEL,WORK(KAPPA1))

               XRLXAB=DDOT(2*NALLAI(ISYMB),WORK(KAPPA1),1,WORK(KG1),1)


               IDXR = IR1TAMP(LABELB,LRLXB,+FREQB,ISYMB)
               CALL CC_GETHFGD(IDXR,'R1 ',LRTLBL,IDUM,IDUM,RDUM,
     *                         ISYLRT,FRQLRT,LORXLRT,IDUM,NLRTLBL,
     *                         MAXTLBL,IREAL,WORK(KCMO),WORK(KUDV),
     *                         WORK(KXINDX),FRVAL,WORK(KG1),LWRKG1)

               IOPT = 1
               IDXR = IR1TAMP(LABELA,LRLXA,+FREQA,ISYMA)
               CALL CC_RDHFRSP('R1',IDXR,ISYMA,IOPT,MODEL,WORK(KAPPA1))

               XRLXBA=DDOT(2*NALLAI(ISYMA),WORK(KAPPA1),1,WORK(KG1),1)

               WORK(KPRP) = XRLXAB
      
               ERROR = XRLXBA - DBLE(ISAPROP) * XRLXAB

               IF (LOCDBG.OR.DABS(ERROR).GT.THRLEQ.OR.IPRINT.GT.1) THEN
                  WRITE (LUPRI,*) 'CC_LR>', LABELA,FREQA,LABELB,FREQB
                  WRITE (LUPRI,*) 'CC_LR> ',XRLXAB ,XRLXBA, ERROR, THRLEQ
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

               LUFCK = 0
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

               LUFCK = 0
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
     *            'the anitsymmetric contribution corresponds to ',
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
            CALL CC_PRPC(RESULT,MODELP,2,
     *                   LABELA,LABELB,LABELA,LABELB,
     *                   BLRFR(IFREQ),BLRFR(IFREQ),BLRFR(IFREQ),ISYMAB)
         END DO

         IF (DIPPOL.AND.(LABELA(2:4).EQ.'DIP')
     *             .AND.(LABELB(2:4).EQ.'DIP')) THEN
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
      IF (MCC2) MODEL = 'MCC2      '
C
      IF (.NOT.(CCS.OR.CC2.OR.CCSD)) THEN
         WRITE(LUPRI,'(A)')
     &        ' CC_SOLRSP: Do not want to calculate anything'
     *                  //' else than CCS, CC2 and CCSD properties '
         CALL QUIT('Model not CCS, CC2, MCC2 or CCSD in CC_SOLRSP')
      ENDIF
C
C--------------------------------------
C        Find the frequency components.
C--------------------------------------
C
         WRITE(LUPRI,'(//,1X,A6,A,F10.6,/)') MODEL(1:6),
     *       'polarizability for frequency: ',FRQ
         CALL OUTPUT(POL,1,3,1,3,3,3,1,LUPRI)
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
      END
c*DECK CC_LRESID
       SUBROUTINE CC_LRESID(WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              linear response residue calculation.
C
C              CCS, CC2, CCSD, MCC2
C
C     Modified version for general linear response properties
C     Ove Christiansen November 1996.
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
      LOGICAL LCALC,LDIP
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*10,MODELP*10
      CHARACTER LABELA*8,LABELB*8
C
#include <leinf.h>
C
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
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
         MODEL = 'CC2'
      ENDIF
      IF (MCC2) THEN
         MODEL = 'MCC2'
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
      KEND1    = KOSCSF + NBLRPRP 
      LEND1    = LWORK  - KEND1
      CALL DZERO(WORK(KOSCS),NALRPRP)
      CALL DZERO(WORK(KOSCSF),NBLRPRP)
C
C----------------------------------------
C     Calculate linear response residues. 
C----------------------------------------
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
C------------------------------
C           Calculate residues.
C------------------------------
C
            KOSCSOF  = NLRSOP*(IRSD-1) + IOPER + KOSCS - 1
            KOSCSOF2 = NLRSOP*(IRSD-1) + IOPER + KOSCSF - 1
            CALL CC_LRSD(LABELA,ISYMA,
     *                   LABELB,ISYMB,
     *                   ISTATE,WORK(KOSCSOF),WORK(KOSCSOF2),
     *                   WORK(KEND1),LEND1)
C
          ENDIF
 2000   CONTINUE
 1000 CONTINUE
C
C-----------------------------------------------
C     Output Linear response properties.
C     IF DIPOSCS put into polarizability tensor.
C-----------------------------------------------
C
      KOSCS2 = KEND1
      KTRS   = KOSCS2 + NEXCI*3*3
      KVELST = KTRS   + NEXCI*3*3
      KVELST2= KVELST + NEXCI*3*3
      KEND2  = KVELST2+ NEXCI*3*3
      LEND2  = LWORK  - KEND2
      CALL DZERO(WORK(KOSCS2),NEXCI*3*3)
      CALL DZERO(WORK(KTRS),NEXCI*3*3)
      CALL DZERO(WORK(KVELST),NEXCI*3*3)
      CALL DZERO(WORK(KVELST2),NEXCI*3*3)
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
            K1     = NLRSOP*(IRSD-1) + IOPER + KOSCS - 1
            K2     = NLRSOP*(IRSD-1) + IOPER + KOSCSF - 1
            RESIDUE = WORK(K1)*WORK(K2)
            IF (RESIDUE.GE.0.0D0) THEN
              SQRRES=SQRT(RESIDUE)
            ELSE
              SQRRES=-SQRT(-RESIDUE)
            ENDIF
            WRITE(LUPRI,'(1X,A6,A8,A1,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     *      'RES{<<',LABELA,',',LABELB,'>>(',EIGV,')} =',
     *      RESIDUE,' ( ',SQRRES,')'
            IF (OSCSTR) THEN
              IADR1 = 0
              IADR2 = 0
              IF (LABELA(1:5).EQ.'XDIPL') IADR1 = 1
              IF (LABELA(1:5).EQ.'YDIPL') IADR1 = 2
              IF (LABELA(1:5).EQ.'ZDIPL') IADR1 = 3
              IF (LABELB(1:5).EQ.'XDIPL') IADR2 = 1
              IF (LABELB(1:5).EQ.'YDIPL') IADR2 = 2
              IF (LABELB(1:5).EQ.'ZDIPL') IADR2 = 3
              IF ((IADR1+IADR2).GE.2) THEN
                IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KOSCS2-1
                WORK(IOSCS2) = RESIDUE     
              ENDIF
            ENDIF
            IF (VELSTR) THEN
              IADR1 = 0
              IADR2 = 0
              IF (LABELA(1:5).EQ.'XDIPV') IADR1 = 1
              IF (LABELA(1:5).EQ.'YDIPV') IADR1 = 2
              IF (LABELA(1:5).EQ.'ZDIPV') IADR1 = 3
              IF (LABELB(1:5).EQ.'XDIPV') IADR2 = 1
              IF (LABELB(1:5).EQ.'YDIPV') IADR2 = 2
              IF (LABELB(1:5).EQ.'ZDIPV') IADR2 = 3
              IF ((IADR1+IADR2).GE.2) THEN
                IOSCS2 = 3*3*(ISTATE-1)+3*(IADR2-1)+IADR1+KVELST-1
                WORK(IOSCS2) = RESIDUE     
              ENDIF
            ENDIF
          ELSE
            RESIDUE = 0.0D0
            SQRRES  = 0.0D0
          ENDIF
          IF (LABELA.EQ.LABELB) THEN
             CALL CC_PRPC(SQRRES,MODELP,-1,
     *                    LABELA,LABELB,LABELA,LABELB,
     *                    EIGV,EIGV,EIGV,ISYMEA)
          ENDIF
        END DO  
      END DO   

C
C--------------------------------------------------------
C     Perform analysis for dipole and velocity strenghts.
C--------------------------------------------------------
C
      CALL DCOPY(3*3*NEXCI,WORK(KOSCS2),1,WORK(KTRS),1)
      CALL DCOPY(3*3*NEXCI,WORK(KVELST),1,WORK(KVELST2),1)
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
           LDIP   = .TRUE. 
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
           LDIP   = .FALSE.
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
      LOGICAL LCALC,LDIP
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_OSCPRI: Output transition properties ' )
      ENDIF
C
C------------------------------------------
C     write out transition strength matrix.
C------------------------------------------
C
      IF (LCALC) THEN 
C
      WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     *    'Transition strength matrix for state nr.',IEX,
     *     ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
      IF (LDIP) THEN
        WRITE(LUPRI,'(1X,A)') 'Dipole gauge'
      ELSE 
        WRITE(LUPRI,'(1X,A)') 'Dipole gauge'
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
C------------------------------------------
C     First scale it - then
C     write out oscillator strength matrix.
C------------------------------------------
C
      IF (LDIP) THEN
         FACT = EIGV*2.0D0/3.0D0
      ELSE
         FACT = 2.0D0/(3.0D0*EIGV)
      ENDIF
      CALL DSCAL(3*3,FACT,OSC,1)
      WRITE(LUPRI,'(//,1X,A6,A,I3,A,I2,/,A7,A,F11.8,/)') MODEL(1:6),
     *    'oscillator strength matrix for state nr.',IEX,
     *     ' of symmetry',ISYM,MODEL(1:6),'excitation energy:',EIGV
      CALL OUTPUT(OSC,1,3,1,3,3,3,1,LUPRI)
      CALL  TNSRAN(OSC,PVAL,PAXIS,
     *             ALFSQ,BETSQ,ITST,ITST2,
     *             APAR2,APEN2,XKAPPA,IPAR)
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
C
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
         IMULT = 1
         WRITE(LUOSC,9988) IMULT,REP(ISYM-1),IEX,TRA,OSCS,CDIP
c     ELSE
c        WRITE(LUOSC,9989) IEX,TRA,OSCS,CDIP
c     ENDIF
C
      ELSE IF (.NOT.LCALC) THEN
        CDIP = '   ?   '
c       IF ( IEX .EQ. 1) THEN
           WRITE(LUOSC,9986) IMULT,REP(ISYM-1),IEX,'Not calculated',
     *                       'Not calculated',CDIP
c       ELSE
c          WRITE(LUOSC,9987) IEX,'Not calculated','Not calculated',CDIP
c       ENDIF
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
      SUBROUTINE CC_TSTAV(ILSTNR,VEC,WORK,LWORK)
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
      PARAMETER( TWO = 2.0D00,XHALF = 0.5D00 )
      DIMENSION WORK(LWORK),VEC(*)
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
C--------------------------
C     Calculate eta vector.
C--------------------------
C
      IF (.NOT.CCS) THEN
         KETA   = 1 
         KWRK1  = KETA  + NT1AM(ISYMOP) + NT2AM(ISYMOP)
         LWRK1  = LWORK - KWRK1
         IF (LWRK1 .LT. 0 )
     &        CALL QUIT('Too little workspace in CC_TSTAV-1')
         CALL CC_ETA(WORK(KETA),WORK(KWRK1),LWRK1)
C
C-------------------------------------------------------------
C        Calculate response contribution to expectation value.
C-------------------------------------------------------------
C
         PROPRSP = DDOT(NT1AM(ISYMOP)+NT2AM(ISYMOP),WORK(KETA),1,VEC,1)
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
      KONEP  = 1
      KT1AM  = KONEP  + N2BST(ISYMOP)
      KLAMDP = KT1AM  + NT1AMX
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
      LWRK1  = LWORK  - KEND1
      IF ( LWRK1 .LT. 0 )
     *     CALL QUIT(' Too little workspace in CC_AVE')
C
      CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
      CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,1,LBL)
C
      IF (.NOT.CCS) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
      ELSE IF (CCS ) THEN
         CALL DZERO(WORK(KT1AM),NT1AMX)
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND1),LWRK1)
C
      ISYFAO = 1
      ISYMPA = 1
      ISYMHO = 1
C
      CALL CC_FCKMO(WORK(KONEP),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC_AVE: Pert. matrix' )
         CALL CC_PRFCKMO(WORK(KONEP),ISYMOP)
      ENDIF
C
      XVALUE = 0.0D0
C
      DO 200 ISYMI = 1, NSYM
C
         DO 210 I = 1, NRHF(ISYMI)
C
            KII  = KONEP - 1 + IFCRHF(ISYMI,ISYMI)
     *           + NORB(ISYMI)*(I - 1) + I
C
            XVALUE = XVALUE + 2.0D00*WORK(KII)
C
  210    CONTINUE
  200 CONTINUE
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
        LUSIFC = 0
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
      IF (MCC2) MODEL = 'MCC2      '
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
        LUSIFC = 0
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
         CALL RGORD(3,3,PVAL,WI,PAXIS)
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
     *                   FRQY,FRQZ,FRQU,ISYMIN)
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
C
      LOGICAL EXIST,L1,L2,L3,L4,LI1,LI2
      PARAMETER (TOLFRQ=1.0D-08,ONE=1.0D0,XMONE=-1.0D0,TOLEXCI =1.0D-02)
      CHARACTER LABEL*10, LABX*8, LABY*8, LABZ*8, LABU*8
C
C--------------------------------------------------
C     Test if this property is already on the list.
C     In that case find address else update NPRPC
C--------------------------------------------------
C
      EXIST = .FALSE.
C
C     New decision....
C
c     DO IPRPC = 1, NPRPC
c       IF ((CPRPC(IPRPC).EQ.LABEL).AND.(NORD.EQ.IOPRPC(IPRPC))) THEN
c        L1 = .TRUE. 
c        L2 = .TRUE. 
c        L3 = .TRUE. 
c        L4 = .TRUE. 
c        LI1 = .TRUE.
c        LI2 = .TRUE.
c
c        IF ((NORD.GE.1).AND.(LABX.EQ.LAB0(IOPRPC(IPRPC)))) THEN
c          L1 = .TRUE.
c        ELSE 
c          L1 = .FALSE.
c        ENDIF
c        IF ((NORD.GE.2).AND.(LABY.EQ.LAB1(IOPRPC(IPRPC))).AND.
c    *       (ABS(FRQY-FRQ1(IPRPC)).LT.TOLFRQ)) THEN
c          L2 = .TRUE.
c        ELSE 
c          L2 = .FALSE.
c        ENDIF
c        IF ((NORD.GE.3).AND.(LABZ.EQ.LAB2(IOPRPC(IPRPC))).AND.
c    *       (ABS(FRQZ-FRQ2(IPRPC)).LT.TOLFRQ)) THEN
c          L3 = .TRUE.
c        ELSE 
c          L3 = .FALSE.
c        ENDIF
c        IF ((NORD.GE.4).AND.(LABU.EQ.LAB3(IOPRPC(IPRPC))).AND.
c    *       (ABS(FRQU-FRQ3(IPRPC)).LT.TOLFRQ)) THEN
c          L4 = .TRUE.
c        ELSE 
c          L4 = .FALSE.
c        ENDIF
c        IF ((NORD.EQ.-1).AND.
c    *       (LABX.EQ.LAB0(IOPRPC(IPRPC))).AND.
c    *       (LABY.EQ.LAB1(IOPRPC(IPRPC))).AND.
c    *       (ABS(FRQY-FRQ1(IPRPC)).LT.TOLEXCI)) THEN
c          LI1 = .TRUE.
c        ELSE
c          LI1 = .FALSE.
c        ENDIF
c        IF ((NORD.EQ.-2).AND.
c    *       (LABX.EQ.LAB0(IOPRPC(IPRPC))).AND.
c    *       (LABY.EQ.LAB1(IOPRPC(IPRPC))).AND.
c    *       (ABS(FRQY-FRQ1(IPRPC)).LT.TOLEXCI).AND.
c    *       (ABS(FRQZ-FRQ2(IPRPC)).LT.TOLEXCI)) THEN
c          LI2 = .TRUE.
c        ELSE
c          LI2 = .FALSE.
c        ENDIF
c       ENDIF
c       EXIST = L1.AND.L2.AND.L3.AND.L4.AND.LI1.AND.LI2
c       IF (EXIST) THEN
c         IHIT = IPRPC
c         GOTO 100
c       ENDIF
c     ENDDO
c
c100  CONTINUE
C
      IF (EXIST) THEN
         IPRPC = IHIT
      ELSE
         NPRPC = NPRPC + 1
         IPRPC = NPRPC
      ENDIF 
      IF (NPRPC.GT.MXPRPC) THEN
        WRITE(LUPRI,*) ' Too many properties on PRPC list - ',
     *    ' I will not add more!'
      ENDIF
C
C--------------------------------------------------
C     Overwritehis property is already on the list.
C--------------------------------------------------
C
      IOPRPC(IPRPC) = NORD
      CPRPC(IPRPC)  = LABEL
C
      LAB0(IPRPC) = LABX
      LAB1(IPRPC) = LABY
      LAB2(IPRPC) = LABZ
      LAB3(IPRPC) = LABU
C
      FRQ1(IPRPC) = FRQY
      FRQ2(IPRPC) = FRQZ
      FRQ3(IPRPC) = FRQU
C
      ISAVSY(IPRPC) = ISYMIN
      PRPC(IPRPC) = PROP
C
      END
