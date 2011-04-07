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
c  /* Deck cc_drv */
      SUBROUTINE CC_DRV(WORK,LWORK)
C
C
C     Ove Christiansen 14-11-1995
C
C
C     Purpose: Driver routine for coupled cluster modules in
C              calculation of
C              ground state energies, excitation energies,
C              oscilator strenghts,
C              frequency dependent polarizabilities and first order
C              properties without orbital relaxation.
C
C              Pass control to;
C
C              CCSD_ENERGY for (integral direct) calculation of
C                          groundstate energies and amplitudes.
C
C              CC_EXCI     for (integral direct) calculation of
C                          excitation energies and amplitudes.
C
C              CC_FOP      for integral direct calculation of first
C                          order properties.
C                          Also calculates the left amplitudes.
C
C              CC_EXGR     for (integral direct) calculation of excited
C                          state gradient properties.
C                          Currently only first order properties.
C
C              CC_LR &     Linear Response or Second-Order Properties:
C              CC_SOP      polarizabilities, property gradients etc.
C
C              CC_LRESID & Control calculation of residues and
C              CC_SOPR     oscillator strenths.
C
C              CC_HYPPOL   calculation of quadratic response properties.
C
C              CC_2HYP     calculation of cubic response properties.
C
C              CC_3HYP     calculation of quartic response properties.
C
C              CC_TPA      Two-Photon Absorption strengths and moments.
C
C              CC_TMCAL    calculation of third order moments
C
C              CC_MCDCAL   calculation of B terms in MCD
C
C     Loop over all models specified. The order is supposed
C     to be the optimal one from a computational point of view.
C     Energy calculations are restarted from the preceding
C     calculation. Excitation vectors is also saved
C     and the excitation energy calculation is restarted
C     from preceding model solution as well. Similar for
C     response vectors in property calculations.
C
C     (ground state) geometry optimizations using the
C     analytic gradient are done for the first (i.e. lowest
C     order) model specified in the input.
C
C     Order: CIS,CCS,MP2,CC(2)(=CIS(D)),CC2,CCD,CCSD,
C            CCSD(T),CCSDR(T),CCSDR(1a),CCSDR(1b),
C            CC(3),CCSDR(3), CC3,CC1A,CC1B
C            (see list below for name of logical flags.)
C
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <iratdef.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccinftap.h>
#include <inftap.h>
#include <ccgr.h>
#include <cclrinf.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccfdgeo.h>
#include <ccslvinf.h>
#include <ccexci.h>
#include <prpc.h>
#include <ccfop.h>
#include <ccfro.h>
#include <eribuf.h>
#include <cbirea.h>
#include <r12int.h>
#include <grdccpt.h>
#include <ccfield.h>
#include <qm3.h>
C
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
C
      PARAMETER(NMODEL = 18)
      LOGICAL   LCCSAV(0:NMODEL), LTCEXC, LPTEXC
      LOGICAL   LTBAR0, LRSPSYM, LLEFTEIG
      LOGICAL   MOLGRD, LHTF, NOAUXBSAVE
      LOGICAL   CC3RSP, CCR12RSP, CCR12LIM
      PARAMETER (CC3RSP = .TRUE.)
      DIMENSION WORK(LWORK)
      CHARACTER FILPQIM*(6), FILXYM*(6)
      PARAMETER (FILPQIM = 'CCPQ00', FILXYM = 'CC_XYM')
      CHARACTER*3 APROXR12
      CHARACTER*8 NAME(8)
      DATA NAME  /'CCAOIN_1','CCAOIN_2','CCAOIN_3','CCAOIN_4',
     *            'CCAOIN_5','CCAOIN_6','CCAOIN_7','CCAOIN_8'/
      COMMON/SORTIO/LUAOIN(8)

      INTEGER NCCEXCI2(8,3),ISYOFE2(8),ITROFE2(8)
C
      INTEGER ISYEXC2(MAXEXCI),IMULTE2(MAXEXCI)
C
      CALL QENTER('CC_DRV')
      CALL GETTIM(TSTART,WSTART)
C
      IF (.NOT. (ONLYMO)) THEN
      CALL CCTAPINI
      WRITE (LUPRI,'(A,/)') '  '
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                    START OF COUPLED CLUSTER CALCULATION  '//
     *'                   *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A,/)')
     *'*********************************************************'//
     *'**********************'
C
      IDUM = 0
C
      IF (IPRINT.GE.5) WRITE(LUPRI,'(/1X,A,I15,/)') 'LWORK = ',LWORK
      END IF
C
      IF (ONLYMO) THEN
        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
        REWIND LUSIFC
C
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ (LUSIFC) NSYM, NORBTS, NBAST, NLAMDS, (NRHFS(I),I=1,NSYM),
     &                (NORBS(I),I=1,NSYM), (NBAS(I),I=1,NSYM), DUM, DUM
C
        CALL GPCLOSE(LUSIFC,'KEEP')
C
        CALL CCMM_REP1(WORK,LWORK)
C
C       Unit no. for info file set to zero
C
        LUMMMO = -1
C
        CALL QUITMM('SCF orbitals and energies calculated')
C
      END IF
C----------------------------------------------------
C     read some flags from /ABAINF/ common block
C     and close SIRIFC in case of gradient requested.
C----------------------------------------------------
C
      CALL CCRDABAINF(MOLGRD)
      IF (LUSIFC .GT. 0) CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------------------
C     initialize common block for Cray I/O:
C----------------------------------------------------------------------
C
      CALL INITWIO()
C
C----------------------------------------------------------------------
C     Open files needed open throughout entire CC calculation.
C----------------------------------------------------------------------
C
      LURES = -1
      CALL GPOPEN(LURES,'CC_RES','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LURES)

      LUPRPC = -1
      CALL GPOPEN(LUPRPC,'CC_PRPC','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUPRPC)

      LUIAJB = -1
      CALL GPOPEN(LUIAJB,'CCSD_IAJB','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUIAJB)
C
C---------------------------------
C     Initialize energy variables.
C---------------------------------
C
      ECCGRS = 0.0D0
      OMECCX = 0.0D0
      ECCXST = 0.0D0
C
C----------------------------------------------------------------------*
C     Save information on the number of states etc. for subsequent
C     calls in numerical differentiations.
C----------------------------------------------------------------------*
C
      NEXCI2 = NEXCI
      DO IMULT = 1, 3, 2
       DO ISYM =1, 8
         NCCEXCI2(ISYM,IMULT) = NCCEXCI(ISYM,IMULT)
         IF (IMULT.EQ.1) ISYOFE2(ISYM)  = ISYOFE(ISYM)
         IF (IMULT.EQ.3) ITROFE2(ISYM)  = ITROFE(ISYM)
       ENDDO
      ENDDO
      DO IEX=1,MAXEXCI
         ISYEXC2(IEX) = ISYEXC(IEX)
         IMULTE2(IEX) = IMULTE(IEX)
      ENDDO
      NPRPC = 0
      NPRMI = 0
C
C---------------------------------------------------
C     Set some logicals needed for R12 calculations:
C---------------------------------------------------
C
      R12CAL = R12CAL .OR. LMULBS

      ! request that some intermediates needed for CC2-R12
      ! are calculated during MP2-R12 calculation
      CC2R12INT = R12CAL .AND. CC2
      CCSDR12INT = R12CAL .AND. (CCSD.OR.CCPT.OR.CCP3.OR.CC3)

      ! switch on R12 flag used in CC routines
      CCR12 = R12CAL

      ! CCR12 not possible as CCS:
      IF (CCR12.AND.CCS) CALL QUIT('CCS not reasonable with R12')
C
C-------------------------------------
C     Set logicals: 4,6,8 is obsolete.
C-------------------------------------
C
      LCCSAV(0)  = CIS
      LCCSAV(1)  = CCS
      LCCSAV(2)  = MP2
      LCCSAV(3)  = CCP2
      LCCSAV(4)  = .FALSE.
      LCCSAV(5)  = CC2
      LCCSAV(6)  = .FALSE.
      LCCSAV(7)  = CCD
      LCCSAV(8)  = .FALSE.
      LCCSAV(9)  = CCSD
      LCCSAV(10) = CCPT
      LCCSAV(11) = CCRT
      LCCSAV(12) = CCR1A
      LCCSAV(13) = CCR1B
      LCCSAV(14) = CCP3
      LCCSAV(15) = CCR3
      LCCSAV(16) = CC3
      LCCSAV(17) = CC1A
      LCCSAV(18) = CC1B
      APROXR12   = '   '
C
      LTCEXC = .FALSE.
C
      MMOD = 0
      DO I = 0, NMODEL
        IF (LCCSAV(I)) MMOD = MMOD + 1
      END DO
C
C------------------------------------------
C     Call the CCSD initialization routine.
C------------------------------------------
C
      ISYMOP = 1
      LABEL  = 'TRCCINT '
      IF (LMULBS) THEN
         NOAUXBSAVE = NOAUXB
         NOAUXB = .TRUE.
      END IF
      CALL CCSD_INIT1(WORK,LWORK)
C
C-----------------------------------------------------------------
C     flag for ground state zeroth-order Lagrange multiplier:
C-----------------------------------------------------------------
      LTBAR0 = (CCLRSD.OR.CCLR.OR.CCQR.OR.CCCR.OR.CC4R.OR.
     &          CC5R.OR.CCQR2R.OR.CCEXGR.OR.CCTM.OR.
     &          CCMCD.OR.CCEXLR.OR.MOLGRD.OR.CCSLV.OR.CCDERI.OR.
     &          CCOPA.OR.CCTPA.OR.CCXOPA)

C-----------------------------------------------------------------
C     flag for kappa-bar-0 Lagrange multiplier:
C-----------------------------------------------------------------
      RELORB = MOLGRD .OR. CCDERI .OR. RELORB

C-----------------------------------------------------------------
C     consistency check: no gradient for solvent...
C-----------------------------------------------------------------
      IF ((MOLGRD.OR.CCDERI) .AND. CCSLV) THEN
        WRITE(LUPRI,*) 'No analytic derivatives in solvent available!'
        IF (MOLGRD) THEN
          CALL QUIT("No optimiz. with anal. gradient in solvent.")
        END IF
        IF (CCDERI) THEN
          WRITE(LUPRI,*) 'CCDERI flag will be ignored.'
          CCDERI = .FALSE.
        END IF
      END IF

C-----------------------------------------------------------------
C     flag for deleting AOTWOINT integral file:
C         KEEPAOTWO  = 0  --> delete AOTWOINT in CCSD_SORTAO
C         KEEPAOTWO  < 2  --> delete AOTWOINT in CC_GRAD
C         KEEPAOTWO >= 2  --> keep AOTWOINT file until the end
C-----------------------------------------------------------------
      IF (RELORB .AND. MMOD.LE.1) KEEPAOTWO = MAX(KEEPAOTWO,1)
      IF (RELORB .AND. MMOD.GT.1) KEEPAOTWO = MAX(KEEPAOTWO,2)

C-----------------------------------------------------------------
C     flag for response calculation (-> execute CCRSPSYM routine):
C-----------------------------------------------------------------
      LRSPSYM = (CCLRSD.OR.CCLR.OR.CCQR.OR.CCCR.OR.CC4R.OR.
     &           CC5R.OR.CCQR2R.OR.CCEXGR.OR.CCTM.OR.
     &           CCMCD.OR.CCEXLR.OR.CCOPA.OR.CCTPA.OR.CCXOPA)

      RSPIM   = RSPIM .OR. LRSPSYM .OR. LTBAR0

C-----------------------------------------------------------------
C     set flag for CCR12 response beyond excitation energies
C     (controls the calculation of certain integrals)
C-----------------------------------------------------------------
      CCR12RSP = CCR12 .AND. (LRSPSYM .OR. NONHF .OR. CCFOP) .AND.
     &           .NOT. R12PRP
      IF (CCR12RSP.AND.(IANCC2.NE.1)) THEN
        WRITE(LUPRI,*) 'CCR12RSP = ',CCR12RSP
        WRITE(LUPRI,*) 'IANR12   = ',IANCC2
        CALL QUIT('Sorry, only Ansatz 1 implemented '//
     &            'for CC response')
      END IF

C-----------------------------------------------------------------
C     flag for left eigenvector calculation:
C-----------------------------------------------------------------
      LLEFTEIG = CCOPA .OR. CCTPA .OR. CCXOPA .OR.
     &           CCMCD .OR. CCTM .OR. CCEXLR

C-----------------------------------------------------------------
C     LPTEXC  - to control if any perturbative triples corrections
C     to excitation energies and thus calculated left vectors.
C-----------------------------------------------------------------
C
      LPTEXC = .FALSE.
      IF (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B) LPTEXC = .TRUE.
C
C----------------------------------------------------------------
C     CCR12LIM - controls if intermediates for CCR12 left transf.
C                are needed
C----------------------------------------------------------------
C
      CCR12LIM = CCR12.AND.(CCR12RSP.OR.LLEFTEIG.OR.LHTR)
C
C------------------------------------------
C     Sort AO-integrals into distributions.
C------------------------------------------
C
      IF ((.NOT. CHOINT) .AND. (.NOT. DIRECT)) THEN
         CALL CCSD_SORTAO(WORK,LWORK)
      ENDIF
C
C-----------------------------------------------------------
C     Calculate the ( ia | jb ) integrals and write to disk.
C-----------------------------------------------------------
C
      IF ((.NOT. CHOINT) .AND. (.NOT.LISKIP)) THEN
         IF (LMULBS.AND.MMOD.LE.1.AND.MP2.AND.(.NOT.NATVIR).AND.
     *       (.NOT.R12PRP)) THEN
           ! if only MP2-R12 the (IA|JB) integrals are not needed
           ! on file LUIAJB
           CONTINUE
         ELSE IF (LMULBS.AND.HERDIR) THEN
           CALL QUIT('CC-R12 with HERDIR not implemented')
         ELSE
           KT1AM = 1
           KT2AM = KT1AM + NT1AMX
           KEND  = KT2AM + NT2AMX
           LWRK  = LWORK - KEND
           IF (LWRK .LT. 0) CALL QUIT ('Insufficient memory in CC_DRV')
           CALL DZERO(WORK(KT1AM),NT1AMX)
           LHTF = .FALSE.
           CALL CCSD_IAJB(WORK(KT2AM),WORK(KT1AM),LHTF,
     *                         .FALSE.,.FALSE.,WORK(KEND),LWRK)
           REWIND(LUIAJB)
           CALL WRITI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KT2AM))
         ENDIF
      ENDIF
C
C--------------------------------------------------------------------
C     Calculate semi-natural virtual orbitals from MP2-Guess for R12:
C--------------------------------------------------------------------
C
      IF ((.NOT.LISKIP).AND.NATVIR) THEN
         CALL CC_R12NO(WORK(KT1AM),WORK(KT2AM),WORK(KEND),LWRK)
      END IF
C
      IF (LMULBS) NOAUXB = NOAUXBSAVE
C
C----------------------
C     Loop over models.
C----------------------
C
      INTTR  = .FALSE.
      IMOD   = 0
C
      DO 100 I = 0, NMODEL
C
         CALL CC_FALSE()
C
         IF (.NOT.LCCSAV(I)) THEN
C
            IF ((I .EQ. NMODEL).AND.(IMOD.EQ.0)) THEN
               WRITE(LUPRI,'(1x,A)') 'No CC model specified - '//
     *                            'default is CCSD'
               CCSD = .TRUE.
            ELSE
               GOTO 100
            ENDIF
C
         ENDIF
C
         IF (I.EQ.0)  CIS    = LCCSAV(0)
         IF (I.EQ.0)  CCS    = LCCSAV(0)
         IF (I.EQ.1)  CCS    = LCCSAV(1)
         IF (I.EQ.2)  MP2    = LCCSAV(2)
         IF (I.EQ.3)  CCP2   = LCCSAV(3)
         IF (I.EQ.4)  GOTO 100
         IF (I.EQ.4)  CC2    = LCCSAV(4)
         IF (I.EQ.5)  CC2    = LCCSAV(5)
         IF (I.EQ.6)  GOTO 100
         IF (I.EQ.7)  CCD    = LCCSAV(7)
         IF (I.EQ.8)  GOTO 100
         IF (I.EQ.9)  CCSD   = LCCSAV(9)
         IF (I.EQ.10) CCPT   = LCCSAV(10)
         IF (I.EQ.11) CCRT   = LCCSAV(11)
         IF (I.EQ.12) CCR1A  = LCCSAV(12)
         IF (I.EQ.13) CCR1B  = LCCSAV(13)
         IF (I.EQ.14) CCP3   = LCCSAV(14)
         IF (I.EQ.15) CCR3   = LCCSAV(15)
         IF (I.EQ.16) CC3    = LCCSAV(16)
         IF (I.EQ.17) CC1A   = LCCSAV(17)
         IF (I.EQ.18) CC1B   = LCCSAV(18)
C
         IF (CCS  .AND. LCCSAV(3)) GOTO 100
         IF (MP2  .AND. LCCSAV(3)) GOTO 100
         IF (CCP2 ) MP2 = .TRUE.
C        IF (CCSD .AND. LCCSAV(10)) GOTO 100
         IF (CCSD .AND. LCCSAV(11)) GOTO 100
         IF (CCSD .AND. LCCSAV(12)) GOTO 100
         IF (CCSD .AND. LCCSAV(13)) GOTO 100
         IF (CCSD .AND. LCCSAV(14) .AND.(.NOT.LCCSAV(15))) GOTO 100
         IF (CCP3 .AND. LCCSAV(15)) GOTO 100
C
         IMOD = IMOD + 1
C
C        -----------------------------------------------------------
C        Save integral transformation (if not triples) and
C        restart from previous cc amplitudes and excitation vectors.
C        -----------------------------------------------------------
         IF (IMOD.GT.1) THEN
            STOLD = .TRUE.
            IF ((I.GT.3).AND..NOT.(LCCSAV(1).AND.(IMOD.EQ.2)))
     *           CCRSTR  = .TRUE.
         ENDIF

         ! in geometry optimizations we can restart from the
         ! amplitudes of the previous points
         IF (MOLGRD) CCRSTR = .TRUE.

         IF (I .GT. 9) INTTR = .TRUE.
C
         IF ((I.GT.15).AND.(.NOT.CCSD)) THEN
            CCSDT  =  .TRUE.
         ENDIF
C
         IF (CCR12 .AND. .NOT.(CCS .OR. MP2 .OR. CC2.OR.
     &       CCSD.OR.CCPT.OR.CCP3.OR.CC3)) THEN
           CALL QUIT('R12 not yet implemented for this CC model')
         END IF

C        -----------------------------------------------------------
C        in geometry optimization runs compute the gradient only
C        for the first model
C        -----------------------------------------------------------
         IF (IMOD.GT.1 .AND. MOLGRD) GOTO 100
C
C----------------------------------------------
C        **************************************
C        ***** Zeroth-order CC amplitudes *****
C        **************************************
C----------------------------------------------
C
         IF (I .GT. 13) LTCEXC = .FALSE.
         IF (LTCEXC) GO TO 100
         IF (CCP2) CCS = .FALSE.
         IF (CCR3) CCP3 = .TRUE.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        Start loop over R12 ansaetze and approximations
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
         IANR12 = IANCC2
         IAPR12 = IAPCC2
C
          WRITE(LUPRI,'(//A,I3)') ' CCR12 ANSATZ = ',IANR12
          WRITE(LUPRI,'(/A,I3/)') ' CCR12 APPROX = ',IAPR12

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        1. Start loop over different solvents.
C        2. Start loop for self-consistent CC solvent
C           solution.
C        Do at least one(which could be vacuum).
C        SLV98,OC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
         NSLV = MAX(1,NCCSLV)
C
         DO 9999 ISLV = 1, NSLV
C
          LMAXCU = LMAXCC(ISLV)
          NLMCU  = (LMAXCU+1)*(LMAXCU+1)
          EPSTCU = EPSTCC(ISLV)
          EPOPCU = EPOPCC(ISLV)
          RCAVCU = RCAVCC(ISLV)
          LSTBTR = .FALSE.
          LSLECVG= .FALSE.
          LSLTCVG= .FALSE.
          LSLLCVG= .FALSE.
          ICCSLIT= 0
C
C         Spaghetti GOTO for finding CCSCRF
C
   31     CONTINUE
C
C          ! no energy in gradient step for geometry optimization
           IF (.NOT.MOLGRD) THEN

              CALL CCSD_ENERGY(WORK,LWORK,APROXR12,CCR12RSP,CCR12LIM)
c
c              write(lupri,*) 'after CCSD_ENERGY: APROXR12 = ',APROXR12
c

           END IF
C
C          switch off flags for sort and transformation
C          done in CCSD_ENERGY
C
           INTTR = .FALSE.
C
           IF (CCP2) CCS = .TRUE.
C
C--------------------------------------
C          ****************************
C          ***** CC Lambda vector *****
C          ****************************
C--------------------------------------
C

           IF ( (CCFOP.OR.LTBAR0) .AND. 
     &          (.NOT. (CCR3.OR.CCR1A.OR.CCR1B.OR.CCRT.OR.CCR3))) THEN
C
             IF (CCR12.AND.MP2.AND..NOT.R12PRP) THEN
                 NWARN = NWARN + 1
                 WRITE(LUPRI,'(//A/A)')
     &              'WARNING: Please specify R12PRP in your input...',
     &              'WARNING: CCFOP will be ignored...'
             ELSEIF (.NOT. (CCSDT .AND. (.NOT. CC3))) THEN
                 !...................................
                 !SONIA: I am not sure but let us try
                 !...................................
                 !IF ( I .EQ. 10 ) CALL CC3_FILOP()
                 !...................................
                 IF (R12PRP) THEN
                    DO IPDD = 1,5
                       IF (IPDD .EQ.2.OR.IPDD.EQ.3.OR.IPDD.EQ.5)THEN
                           CALL CC_FOP(IPDD,WORK,LWORK,APROXR12)
                       ENDIF
                    ENDDO
                 ELSE
                    CALL CC_FOP(IDUMMY,WORK,LWORK,APROXR12)
                 ENDIF
             ENDIF
C
           ENDIF
C
C------------------------------------------------------------
C           SLV98,OC
C           End loop for self-consistent CC solvent solution.
C           Force restart after first run.
C           After CC SCRF convergence then left
C           transformations are no longer static.
C------------------------------------------------------------
C
         IF (CCSLV) THEN
            CCRSTR  = .TRUE.
            IF (LSLECVG.AND.LSLTCVG.AND.LSLLCVG) THEN
              LSTBTR = .TRUE.
            ELSE
              GOTO 31
            ENDIF
         ENDIF
C
         IF (CCSLV .AND. DISCEX) THEN
C           IF (CCEXCI .OR. CCEXGR .OR. LRSPSYM .OR. CCLRSD .OR.
C     *         CCQR2R .OR. CCSM .OR. CCLR .OR. CCQR ) THEN
             IF (.NOT. CCMM) THEN
               CALL CCSL_XIETA(WORK,LWORK)
             ELSE IF (CCMM) THEN
               CALL CCMM_XIETA(WORK,LWORK)
             END IF
C           END IF
         END IF
C
C----------------------------------------------------------------------*
C
C        CCDERI -- gradient calculation forced in CC input
C
C        MOLGRD -- gradient required from DALTON (via /ABAINF/)
C                  some short cuts are made in the latter case:
C                  1) the energy calculation and all property/excitation
C                     energy calculations are saved because they were
C                     done in a previous call to the CC program
C                  2) we do the optimization (and thus the gradient)
C                     for the first (i.e. lowest order) model
c
C----------------------------------------------------------------------*
C
         IF (CCDERI .OR. MOLGRD) THEN
            CALL CC_GRADIENT(WORK,LWORK)
            IGRDCCPT = IGRDCCPT + 1

            IF (MOLGRD) GOTO 100
         END IF
C
C---------------------------------------------------------------------
C        If triples corrections have been calculated goto end of loop.
C        or do all triple corrections in one call.
C        But if both ccsd(t) and cc(3) ground state energies are wanted
C        wait to next time.
C---------------------------------------------------------------------
C
         IF (LTCEXC) GO TO 90
C
         IF ((I.GE.10).AND.(I.LE.13)) THEN
C
            LTCEXC = .TRUE.
            CCRT   = LCCSAV(11)
            CCR1A  = LCCSAV(12)
            CCR1B  = LCCSAV(13)
C
         ENDIF
C
         IF ( CCR3 ) CCT = .TRUE.
C
C---------------------------------------------
C        ***********************************
C        ***** CC Excitation energies  *****
C        ***********************************
C---------------------------------------------
C
         IF (CCEXCI) THEN
C
            IF (((.NOT. MP2).OR.CCP2).AND.(.NOT.(CCD.OR.CCPT))) THEN
C
               IF ((.NOT. LHTR).AND.(.NOT.RESKIP)) THEN
                  CALL CC_EXCI(WORK,LWORK,'RE ',APROXR12)
               ENDIF
C
               IF (CCSPIC.AND.CCS) CALL CC_REDEIG(WORK,LWORK,OMPCCS)
               IF (CC2PIC.AND.CC2) CALL CC_REDEIG(WORK,LWORK,OMPCC2)
               IF (CCSDPI.AND.CCSD) CALL CC_REDEIG(WORK,LWORK,OMPCCSD)
C
               IF (LHTR.OR.CCLRSD.OR.CCQR2R.OR.LPTEXC.OR.CCP2.OR.
     *             LLEFTEIG.OR.CCEXGR) THEN
                 IF ((I .LT. 17).AND.(.NOT.LESKIP)) THEN
                   CALL CC_EXCI(WORK,LWORK,'LE ',APROXR12)
                 ELSE
                   WRITE(LUPRI,*) 'No triples left excitation energies '
                 ENDIF
               ENDIF
            ELSE
               WRITE(LUPRI,'(/,1X,A)')
     &              ' No MP2 and CCD excitation energies'
            ENDIF
C
         ENDIF
C
C------------------------------------------------------------
C        Settings for numerical gradient (driven from dalton)
C------------------------------------------------------------
C
         CALL CC_FDGD
C
C------------------------------------------------
C        The rest is not implemented for triples.
C------------------------------------------------
C
         IF (I.GT.9 .AND. I.NE.16) GOTO 90
C
C------------------------------------------------------------
C        ****************************************************
C        ***** general set up for response calculations *****
C        ****************************************************
C------------------------------------------------------------
C
         IF (LRSPSYM) THEN
           Call CCRSPSYM(MOLGRD,WORK,LWORK)
         END IF
C
C---------------------------------------------------------------
C        precalculate expectatation values & eff. Fock matrices:
C---------------------------------------------------------------
C
          IF (CCS .or. CC2 .or. CCSD) THEN
            I1DXORD = 0
            IOPREL  = 0
            CALL CC_GRAD2(I1DXORD,WORK,LWORK)
C           CALL CCEFFFOCK(I1DXORD,IOPREL,WORK,LWORK)
          END IF

C
C-----------------------------------------------
C        ***************************************
C        ***** Solve CC-response equations.*****
C        ***************************************
C-----------------------------------------------
C
         IF (LRSPSYM) THEN
           CALL CC_SOLRSP(WORK,LWORK,APROXR12)
         END IF
C
C-----------------------------------------------------
C        *********************************************
C        ***** Excited state gradient properties *****
C        *********************************************
C-----------------------------------------------------
C
         IF (CCEXGR) THEN
           IF (.NOT. CCSLV) THEN
            IF ((CCS.or.CC2.or.CCSD.or.CC3).and.(.not.CCR12)) THEN
             CALL CC_EXGR(WORK,LWORK)
            ELSE
             WRITE(LUPRI,*) 'Excited state gradient properties'//
     &          'are not implemented for the requested model.'
            END IF
           ELSE 
            WRITE(LUPRI,*) 'Solvent model not implemented for'//
     &                     'excited state gradients.'
           ENDIF
         ENDIF 
C
C-----------------------------------------------
C        ***************************************
C        ***** CC linear response residues *****
C        ***************************************
C-----------------------------------------------
C
         IF (CCLRSD) THEN
           IF ((CCS .or. CC2 .or. CCSD) .and. (.not. CCR12)) THEN
            IF (OLDLRS) THEN
               CALL CC_LRESID(WORK,LWORK)
            ELSE
               CALL CC_SOPR(WORK,LWORK)
            END IF
           ELSE
            WRITE(LUPRI,*) 'Linear response residues'//
     &         ' are not implemented for the requested model.'
           END IF
         END IF

         IF (CCOPA) THEN 
           IF ((CCS .or. CC2 .or. CCSD .or. (CC3.AND.CC3RSP)) 
     &          .and. (.not. CCR12)) THEN
            CALL CC_OPA(WORK,LWORK)  
           ELSE
            WRITE(LUPRI,*) 'One-photon absorption strengths'//
     &         ' are not implemented for the requested model.'
           END IF
         END IF
C
C--------------------------------------------------------
C        ************************************************
C        ***** CC quadratic reponse second residues *****
C        ************************************************
C--------------------------------------------------------
C
         IF (CCQR2R) THEN
           IF ((CCS .or. CC2 .or. CCSD) .and. (.not. CCR12)) THEN
            CALL CC_QR2RSD(WORK,LWORK)
           ELSE
            WRITE(LUPRI,*) 'Quadratic response residues'//
     &         ' are not implemented for the requested model.'
           END IF
         END IF

         IF (CCXOPA) THEN
           IF ((CCS .or. CC2 .or. CCSD .or. (CC3.AND.CC3RSP))
     &          .and. (.not. CCR12)) THEN
            CALL CC_XOPA(WORK,LWORK)
           ELSE
            WRITE(LUPRI,*) 'One-photon absorption strengths'//
     &         ' are not implemented for the requested model.'
           END IF
         END IF
C
C-------------------------------------------------------
C        ***********************************************
C        *****  CC two-photon transition moments   *****
C        ***********************************************
C-------------------------------------------------------
C
         IF (CCTPA) THEN 
           IF ((CCS .or. CC2 .or. CCSD .or. (CC3.AND.CC3RSP))
     &          .and. (.not. CCR12)) THEN
            CALL CC_TPA(WORK,LWORK)  
           ELSE
            WRITE(LUPRI,*) 'two-photon absorption strengths'//
     &         ' are not implemented for the requested model.'
           END IF
         ENDIF
C
C------------------------------------------------------
C        **********************************************
C        ***** CC third-order transition moments  *****
C        **********************************************
C------------------------------------------------------
C
         IF (CCTM) THEN
           IF ((CCS .or. CC2 .or. CCSD) .and. (.not. CCR12) 
     &          .AND. (.NOT. CCSLV)) THEN
            CALL CC_TMCAL(WORK,LWORK)
           ELSE
            WRITE(LUPRI,*) 'third-order transition moments'//
     &         ' are not implemented for the requested model.'
           END IF
         ENDIF
C
C-------------------------------------------------------
C        ***********************************************
C        ***** CC magnetic circular dichroism      *****
C        ***********************************************
C-------------------------------------------------------
C
         IF (CCMCD) THEN
           IF ((CCS .or. CC2 .or. CCSD) .and. (.not. CCR12)
     &          .AND. (.NOT. CCSLV)) THEN
            CALL CC_MCDCAL(WORK,LWORK)
           ELSE
            WRITE(LUPRI,*) 'magnetic circular dichroism'//
     &         ' is not implemented for the requested model.'
           END IF
         ENDIF
C
C-----------------------------------------
C        *******************************
C        ***** CC Polarizabilities *****
C        *******************************
C-----------------------------------------
C
         IF ( CCLR .AND. NBLRFR.GT.0) THEN
C
           IF (CIS .OR. CCS .OR. CC2 .OR. CCSD .OR. (CC3.AND.CC3RSP))
     &     THEN
             IF ((OLDLR .AND. .NOT.CC3) .OR. CIS) THEN
               CALL CC_LR(WORK,LWORK)
             ELSE
               CALL CC_SOP(WORK,LWORK)
             END IF
           ELSE
              WRITE(LUPRI,'(/,1X,A)')
     &             ' No polarizabilities for present model'
           ENDIF
C
         ENDIF
C
C-----------------------------------------
C        *****************************
C        ***** CC Cauchy Moments *****
C        *****************************
C-----------------------------------------
C
         IF ( CCLR .AND. CAUCHY) THEN
           IF ((CCS .or. CC2 .or. CCSD .OR. (CC3.AND.CC3RSP)) 
     &         .AND. (.NOT.CCR12) .AND. (.NOT. CCSLV) ) THEN
             CALL CC_CAUCHY(WORK,LWORK)
           ELSE
             WRITE(LUPRI,'(/,1X,2A)') 'Cauchy moments',
     &        ' are not yet implemented for the requested model'
           ENDIF
         ENDIF
C
C----------------------------------------------------
C        ********************************************
C        ***** CC excited state linear response *****
C        ********************************************
C----------------------------------------------------
C
         IF (CCEXLR) THEN 
           IF ((CCS .or. CC2 .or. CCSD .or. (CC3.AND.CC3RSP))
     &          .and. (.not. CCR12) .AND. (.NOT. CCSLV)) THEN
            CALL CC_EXLR(WORK,LWORK)  
           ELSE
            WRITE(LUPRI,*) 'Excited state linear response'//
     &         ' is not implemented for the requested model.'
           END IF
         END IF
C
C--------------------------------------------------
C        ******************************************
C        ***** CC first Hyperpolarizabilities *****
C        ******************************************
C--------------------------------------------------
C
         IF ( CCQR ) THEN
C
            IF (CCS .or. CC2 .or. CCSD .or. (CC3.AND.CC3RSP)) THEN
               CALL CC_HYPPOL(WORK,LWORK)
            ELSE
             WRITE(LUPRI,'(/,1X,2A)') 'first hyperpolarizabilities',
     &        ' are not yet implemented for the requested model'
            ENDIF
C
         ENDIF
C
C---------------------------------------------------
C        *******************************************
C        ***** CC second Hyperpolarizabilities *****
C        *******************************************
C---------------------------------------------------
C
         IF ( CCCR ) THEN
C
            IF ( (CCS .or. CC2 .or. CCSD .or. (CC3.AND.CC3RSP)) 
     &         .AND. (.NOT. CCSLV) ) THEN
               CALL CC_2HYP(WORK,LWORK)
            ELSE
              WRITE(LUPRI,'(/,1X,A)')
     &              'requested model not yet implemented'
            ENDIF
C
         ENDIF
C
C--------------------------------------------------
C        ******************************************
C        ***** CC third Hyperpolarizabilities *****
C        ******************************************
C--------------------------------------------------
C
         IF ( CC4R ) THEN
C
            IF ( (CCS .or. CC2 .or. CCSD) .AND. (.NOT. CCSLV) ) THEN
               CALL CC_3HYP(WORK,LWORK)
            ELSE
              WRITE(LUPRI,'(/,1X,A)')
     &              'requested model not yet implemented'
            ENDIF
C
         ENDIF
C
C---------------------------------------------------
C        *******************************************
C        ***** CC fourth Hyperpolarizabilities *****
C        *******************************************
C---------------------------------------------------
C
         IF ( CC5R ) THEN
C
            IF ( (CCS .or. CC2 .or. CCSD) .AND. (.NOT. CCSLV) ) THEN
               CALL CC_4HYP(WORK,LWORK)
            ELSE
              WRITE(LUPRI,'(/,1X,A)')
     &              'requested model not yet implemented'
            ENDIF
C
         ENDIF
C
C--------------------------------------------------
C        ******************************************
C        *****  B,F,C,D... matrix test calls  *****
C        ******************************************
C--------------------------------------------------
C
         IF ( .TRUE. ) THEN
C
            IF (CCS .or. CC2 .or. CCSD) THEN
C              !NOTE: Only CC_FTSTNEW, CC_BTST are compatible with CC-R12!
C              CALL CC_BTST(WORK,LWORK,APROXR12)
C              CALL CC_CTST(WORK,LWORK)
C              CALL CC_DTST(WORK,LWORK)
C              CALL CC_FTSTNEW(WORK,LWORK,APROXR12)
C              CALL CC_GTST(WORK,LWORK)
C              CALL CC_HTST(WORK,LWORK)
C              CALL CC_XETST(WORK,LWORK)
C              CALL CC_FBTST(WORK,LWORK)
C              CALL CC_AATST(WORK,LWORK)
C              CALL CC_BATST(WORK,LWORK)
            ELSE
              WRITE(LUPRI,'(/,1X,A)')
     &              'requested model not yet implemented'
            ENDIF
C
         ENDIF
C
  90     CONTINUE
C
C SLV98,OC End of solvent loop
C
 9999    CONTINUE
C
C-------------------------------------------------
C        Close files used in triples calculations.
C-------------------------------------------------
C
!SONIA SONIA SONIA
!SONIA SONIA SONIA: take care. Closed files we need for geoopt????
!SONIA SONIA SONIA

      IF ( I .GT. 9 ) CALL CC3_FILCL()
C
 100  CONTINUE
C
C-------------------------------------
C     Make summary of CC calculations.
C-------------------------------------
C
      CALL CC_RESUME(WORK,LWORK)
C
      CCS    = LCCSAV(1)
      MP2    = LCCSAV(2)
      CCP2   = LCCSAV(3)
      CC2    = LCCSAV(5)
      CCD    = LCCSAV(7)
      CCSD   = LCCSAV(9)
      CCPT   = LCCSAV(10)
      CCRT   = LCCSAV(11)
      CCR1A  = LCCSAV(12)
      CCR1B  = LCCSAV(13)
      CCP3   = LCCSAV(14)
      CCR3   = LCCSAV(15)
      CC3    = LCCSAV(16)
      CC1A   = LCCSAV(17)
      CC1B   = LCCSAV(18)
C
C     Close files.
C
C
      CALL GPCLOSE(LURES,'KEEP')
      CALL GPCLOSE(LUPRPC,'KEEP')
      CALL GPCLOSE(LUIAJB,'KEEP')
      IF (LUAORC(0).GT.0) CALL GPCLOSE(LUAORC(0),'DELETE')
      IF (LUINTR.GT.0)    CALL GPCLOSE(LUINTR,'DELETE')
      IF (LUINTA.GT.0) THEN
        IF (MOLGRD.AND.(.NOT.DIRECT)) THEN
          CALL GPCLOSE(LUINTA,'DELETE')
        ELSE
          CALL GPCLOSE(LUINTA,'KEEP')
        END IF
      END IF
C
      DO I = 1, 8
        IF (LUAOIN(I) .GT. 0) THEN
          IF ((.NOT.KEEPAOIN).OR.MOLGRD) THEN
            CALL WCLOSE2(LUAOIN(I),NAME(I),'DELETE')
          ELSE
            ! needed for gradient calculation...
            CALL WCLOSE2(LUAOIN(I),NAME(I),'KEEP')
          END IF
        END IF
      END DO

      CALL WCLOSEALL()
C
C------------------------------------------------
C     Restore information on number of roots etc.
C     if there is an additional call to cc_drv
C------------------------------------------------
C
      NEXCI = NEXCI2
      DO IMULT = 1, 3, 2
       DO ISYM =1, 8
         NCCEXCI(ISYM,IMULT) = NCCEXCI2(ISYM,IMULT)
         IF (IMULT.EQ.1) ISYOFE(ISYM)  = ISYOFE2(ISYM)
         IF (IMULT.EQ.3) ITROFE(ISYM)  = ITROFE2(ISYM)
       ENDDO
      ENDDO
      DO IEX=1,MAXEXCI
         ISYEXC(IEX) = ISYEXC2(IEX)
         IMULTE(IEX) = IMULTE2(IEX)
      ENDDO
C
      CALL GETTIM(TEND,WEND)
      TUSED = TEND - TSTART
      WUSED = WEND - WSTART
      IF (IPRSTAT .GT. 0) THEN
         WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for CC :',TUSED,WUSED
         CALL FLSHFO(LUSTAT)
      END IF
      WRITE (LUPRI,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for CC :',TUSED,WUSED
      CALL TSTAMP(' ',LUPRI)
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('CC_DRV')
C
      END
C  /* Deck cc_resume */
      SUBROUTINE CC_RESUME(WORK,LWORK)
C
C     Ove Christiansen 14-11-1995
C
C     Purpose: Resume all results obtained in this call to coupled
C              cluster routines: ground state energies, excitation energies etc.
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <iratdef.h>
#include <dummy.h>
      DIMENSION WORK(LWORK)
      CHARACTER*80  STR
      CHARACTER STHELP*10
      CHARACTER LABEL*10, LABEL1*8
      
#include <codata.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cclr.h>
#include <prpc.h>
#include <ccpack.h>
#include <inftap.h>
#include <ccinftap.h>
#include <numder.h>
C
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./
C
      CALL QENTER('CC_RESUME')
C
      REWIND(LURES)
C
      WRITE (LUPRI,'(A,/)') '  '
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                   SUMMARY OF COUPLED CLUSTER CALCULATION '//
     *'                   *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A,/)')
     *'*********************************************************'//
     *'**********************'
C
 100  READ(LURES,'(A80)',END=200) STR
C
         WRITE(LUPRI,'(A80)') STR
         GOTO 100
C
 200  CONTINUE
C
      IF ( CCEXCI ) THEN
         WRITE(LUPRI,'(/A11,F15.5,A6)') ' 1 a.u. = ',XTEV,' eV.  '
         WRITE(LUPRI,'(A11,F15.5,A6)') ' 1 a.u. = ',XTKAYS,' cm-1.'
      ENDIF
C
      IF (FIRSTCALL) THEN
         FIRSTCALL = .FALSE.
         CALL CC_WPRPC
C
C-----------------------------------------------
C           Save the list ("old prpc") on file.
C           Count the number of properties in
C           each symmetry class
C-----------------------------------------------
C
         DO ISYM=1,8
            NPRPINSYM(ISYM) = 0
         ENDDO
         LUPRPCO = -1
         CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',' ','FORMATTED',
     *        IDUMMY,.FALSE.)
         REWIND(LUPRPCO)
         REWIND(LUPRPC)

         NPRPCO = NPRPC
         DO I = 1, NPRPC
            READ(LUPRPC,
     *         '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *         IPRPC,ISYMIN,NORD,LABEL,PROP,
     *         LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYEX,ISPEX,IEX
            WRITE(LUPRPCO,
     *         '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *         IPRPC,ISYMIN,NORD,LABEL,PROP,
     *         LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYEX,ISPEX,IEX
            NPRPINSYM(ISYMIN) = NPRPINSYM(ISYMIN) + 1
         END DO
         CALL GPCLOSE(LUPRPCO,'KEEP')
C
      ELSE
C
C        Reorder the list of properties to match the firstcall order
C
         CALL CC_PRPCREORDER
C
C        *** For numerical differentiation of cc-properties, write the ***
C        *** property to file.                                         ***
C
         IF (NMRDRP.GT.0) THEN
            KPRPC = 1
            KEND  = KPRPC + NPRPC
            LEND  = LWORK - KEND
            IF (LEND.LE.0)CALL QUIT('Too little workspace in cc_resume')
            REWIND(LUPRPC)

            DO I = 1, NPRPC
               READ(LUPRPC,
     *          '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *          IPRPC,ISYMIN,NORD,LABEL,WORK(KPRPC+I-1),
     *          LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYEX,ISPEX,IEX
            ENDDO
            CALL WRAVFL(WORK(KPRPC),NPRPC,1,1,'CC PROPER',IPRINT)
         END IF
      ENDIF
C
C     Writeout the list of properties to summary section of output file.
C

      CALL FLSHFO(LUPRI)
      IF ((NPRPC.GT.0).AND.(IPRINT.GT.100)) THEN
        CALL PRPRPC(LUPRPC,1,DUMMY,NPRPC)
      ENDIF
C
      IF ( LPACKINT ) THEN
         WRITE (LUPRI,'(//10X,A,G8.2)')
     &     'Threshold used for packing of integrals:',THRPCKINT
         WRITE (LUPRI,'(10X,A,F8.2,A)')
     &     'Time used for packing and unpacking:    ',
     &             PCKTIME, ' seconds'
      END IF
      WRITE (LUPRI,'(A,/)') '  '
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                      END OF COUPLED CLUSTER CALCULATION  '//
     *'                   *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'                     *'
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********************'
      WRITE (LUPRI,'(1x,A,/)')
     *'*********************************************************'//
     *'**********************'
C
      LABEL1 = 'ALL_DONE'
      STHELP = 'THE_END   '
      ZERO=0.0D0
      CALL CC_PRPC(ZERO,STHELP,666,
     *             LABEL1,LABEL1,LABEL1,LABEL1,
     *             ZERO,ZERO,ZERO,1,0,0,0)

      CALL QEXIT('CC_RESUME')
C
      END
C  /* Deck cc_false */
      SUBROUTINE CC_FALSE()
C
C     Ove Christiansen 14-11-1995
C
C     Purpose: set all model flags to false.
C
#include <implicit.h>
#include <ccsdinp.h>
#include <cclr.h>
C
      CALL QENTER('CC_FALSE')
C
      CIS    = .FALSE.
      CCS    = .FALSE.
      MP2    = .FALSE.
      CCP2   = .FALSE.
      CC2    = .FALSE.
      CCD    = .FALSE.
      CCSD   = .FALSE.
      CCPT   = .FALSE.
      CCP3   = .FALSE.
      CCRT   = .FALSE.
      CCR3   = .FALSE.
      CCR1A  = .FALSE.
      CCR1B  = .FALSE.
      CC3    = .FALSE.
      CC1A   = .FALSE.
      CC1B   = .FALSE.
      CCSDT  = .FALSE.
C
      CCT    = .FALSE.
      STCCS  = .FALSE.
C
      CALL QEXIT('CC_FALSE')
C
      END
C  /* Deck cc_false */
      SUBROUTINE CC_FDGD
C
C     Ove Christiansen 13-2-1997
C
C     Purpose: Set energy for finite difference gradient calculation.
C
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccgr.h>
#include <ccfdgeo.h>
#include <infopt.h>
C
      CALL QENTER('CC_FDGD')
C
      IF (IPRINT.GT.2) THEN
        WRITE (LUPRI,'(/,1X,A)')
     *  '*********************************************************'//
     *  '**********'
      END IF
      IF ((IXSTSY.GT.0).AND.(IXSTAT.GT.0)) THEN
         ECORR = ECCXST
         IF (IPRINT.GT.2) THEN
           WRITE(LUPRI,'(A,I3,A,I3)') ' Excited state nr.',IXSTAT,
     *      ' of symmetry ',IXSTSY
           WRITE(LUPRI,'(A,F20.10)') ' Ground state energy: ', ECCGRS
           WRITE(LUPRI,'(A,F20.10)') ' Excitation energy:   ', OMECCX
           WRITE(LUPRI,'(A,F20.10)') ' Excited state energy:',ECCXST
         ENDIF
      ELSE
         ECORR = ECCGRS
         IF (IPRINT.GT.2) THEN
           WRITE(LUPRI,'(/,1X,A)') 'Ground state energy is target '
         ENDIF
      ENDIF
C
      IF (IPRINT.GT.2) THEN
       WRITE(LUPRI,'(/,A,F30.20)')    ' Ecorr =              ',Ecorr
       WRITE(LUPRI,'(/,1X,A)')
     * '*********************************************************'//
     * '**********'
      END IF
C
      CALL QEXIT('CC_FDGD')
C
      END
*=====================================================================*
       SUBROUTINE CCRDABAINF(MOLGRDCC)
*---------------------------------------------------------------------*
*
*    Purpose: read flags from /ABAINF/ common block used in CC
*
*=====================================================================*
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>

      LOGICAL MOLGRDCC

      CALL QENTER('CCRDABAINF')

      MOLGRDCC = MOLGRD

      CALL QEXIT('CCRDABAINF')

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CCRDABAINF                           *
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCTAPINI
*---------------------------------------------------------------------*
*
*     Purpose: Initialize coupled cluster unit numbers and file names
*
*=====================================================================*
#include <implicit.h>
#include <ccinftap.h>
C
      CALL QENTER('CCTAPINI')
C
      LUBFDN  = -1
      LURES   = -1
      LUCSOL  = -1
C
      CALL QEXIT('CCTAPINI')
C
      RETURN
      END

*=====================================================================*
      SUBROUTINE CC_WPRPC
C
C     Ove Christiansen Nov. 1999
C
C     Purpose: Write out operators.
C
#include <implicit.h>
#include <ccsdinp.h>
#include <ccroper.h>
#include <dummy.h>
C
      CALL QENTER('CC_WPRPC  ')
C
      LUPRP = -1
      CALL GPOPEN(LUPRP,'PRPCOP','UNKNOWN',' ','FORMATTED',
     *            IDUMMY,.FALSE.)
C
      DO IROPER = 1,NRSOLBL
         WRITE(LUPRP,'(I3,A8,I3)') IROPER,LBLOPR(IROPER),ISYOPR(IROPER)
      END DO
C
      CALL GPCLOSE(LUPRP,'KEEP')
C
      CALL QEXIT('CC_WPRPC  ')
C
      END

*=====================================================================*
       SUBROUTINE CC_PRPCREORDER
C
C     Ove Christiansen Nov. 1999
C
C     Purpose: Reorder PRPC list to fit old PRPC2 lists.
C              Property list for numerical differentiation.
C
#include <implicit.h>
#include <prpc.h>
#include <priunit.h>
#include <ccinftap.h>
      CHARACTER LABELI*10, LABXI*8, LABYI*8, LABZI*8, LABUI*8
      CHARACTER LABELJ*10, LABXJ*8, LABYJ*8, LABZJ*8, LABUJ*8
C
c     DIMENSION FRQ12(MXPRPC),FRQ22(MXPRPC),FRQ32(MXPRPC)
c     DIMENSION PRPC2(MXPRPC)
c     CHARACTER*8  LAB02(MXPRPC), LAB12(MXPRPC), LAB22(MXPRPC),
c    *             LAB32(MXPRPC)
c     CHARACTER*10 CPRPC2(MXPRPC)
c     DIMENSION IOPRPC2(MXPRPC), ISAVSY2(MXPRPC)
      PARAMETER (TOLPRPFD=1.0D-10,TOLE=1.0D-02)
C
      CALL QENTER('CC_WPRPCREORDER')
C
C  TEST TEST
C     IF (NPRPC.NE.NPRPCO) CALL QUIT('Mistmatch in WPRPC')
C  TEST TEST
C
C-----------------------------------------------
C     Save the current list on ("prpc_tmp") on file.
C-----------------------------------------------
C
      LUPRPCTMP = -1
      CALL GPOPEN(LUPRPCTMP,'CC_PRPC_TMP','UNKNOWN',' ','FORMATTED',
     *           IDUMMY,.FALSE.)
      REWIND(LUPRPCTMP)
      REWIND(LUPRPC)
       DO I = 1, NPRPC
         READ(LUPRPC,
     *     '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *     IPRPCI,ISYMINI,NORDI,LABELI,PROPI,
     *     LABXI,LABYI,LABZI,LABUI,FRQYI,FRQZI,FRQUI,
     *     ISYEXI,ISPEXI,IEXI
         WRITE(LUPRPCTMP,
     *     '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *     IPRPCI,ISYMINI,NORDI,LABELI,PROPI,
     *     LABXI,LABYI,LABZI,LABUI,FRQYI,FRQZI,FRQUI,
     *     ISYEXI,ISPEXI,IEXI
      ENDDO
C
C-----------------------------------------------
C     Open the saved list ("prpc_o") on file.
C-----------------------------------------------
C
      LUPRPCO   = -1
      CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',' ','FORMATTED',
     *            IDUMMY,.FALSE.)
C
C     Rewind luprpc to make ready for the ordered set of data.
C
      REWIND(LUPRPC)
C
C     Rewind luprpco to prepare the reference set of data.
C
      REWIND(LUPRPCO)
      DO 300 I = 1, NPRPC
       READ(LUPRPCO,
     *    '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *    IPRPCI,ISYMINI,NORDI,LABELI,PROPI,
     *    LABXI,LABYI,LABZI,LABUI,FRQYI,FRQZI,FRQUI,
     *    ISYEXI,ISPEXI,IEXI

C
C     Rewind luprpc_tmp to prepare the data that is reordered
C
       REWIND(LUPRPCTMP)
       DO 400 J = 1, NPRPC
        READ(LUPRPCTMP,
     *     '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *     IPRPCJ,ISYMINJ,NORDJ,LABELJ,PROPJ,
     *     LABXJ,LABYJ,LABZJ,LABUJ,FRQYJ,FRQZJ,FRQUJ,
     *     ISYEXJ,ISPEXJ,IEXJ
        IF (LABELJ.EQ.LABELI) THEN
         IF (NORDI.EQ.NORDJ) THEN
          IF ((LABXJ.EQ.LABXI).AND.
     *        (LABYJ.EQ.LABYI).AND.
     *        (LABZJ.EQ.LABZI).AND.
     *        (LABUJ.EQ.LABUI))  THEN
C Take care with excitation energies !!!!
           IF (
     *         ((NORDI.GT.0).AND.
     *          ((ABS(FRQYJ-FRQYI).LT.TOLPRPFD).AND.
     *           (ABS(FRQZJ-FRQZI).LT.TOLPRPFD).AND.
     *           (ABS(FRQUJ-FRQUI).LT.TOLPRPFD))).OR.
     *         ((NORDI.LE.0).AND.
     *          ((ABS(FRQYJ-FRQYI).LT.TOLE).AND.
     *           (ABS(FRQZJ-FRQZI).LT.TOLE).AND.
     *           (ABS(FRQUJ-FRQUI).LT.TOLE)))) THEN
              WRITE(LUPRPC,
     *          '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *          IPRPCI,ISYMINJ,NORDJ,LABELJ,PROPJ,
     *          LABXJ,LABYJ,LABZJ,LABUJ,FRQYJ,FRQZJ,FRQUJ,
     *          ISYEXJ,ISPEXJ,IEXJ
              GOTO 300
           ENDIF
          ENDIF
         ENDIF
        ENDIF
 400   CONTINUE
 300  CONTINUE

      CALL GPCLOSE(LUPRPCTMP,'KEEP')
      CALL GPCLOSE(LUPRPCO,'KEEP')
      CALL QEXIT('CC_WPRPCREORDER')
C
      END
C
C
      SUBROUTINE PRPRPC(LU,IOPT,PROPIN,NUOFPR)
C
C     Outputs the prpc to lupri the list of properties
C     that is stored on the file opened with unit number LU
C     NB: File should be open on entry
C     IF IOPT = 1 the properties on file is output.
C     IF IOPT = 2 the properties in vector propin is output with
C                 labels etc. comming from the prpc list.
C
#include <implicit.h>
#include <priunit.h>
#include <prpc.h>
      CHARACTER STHELP*10
      CHARACTER LABEL*10, LABX*8, LABY*8, LABZ*8, LABU*8
      DIMENSION PROPIN(*)

      WRITE (LUPRI,'(A,/)') '  '
      WRITE (LUPRI,'(A)')
     *' +---------------------------------------------------------+ '
      WRITE (LUPRI,'(A)')
     *' | List of selected calculated molecular properties (a.u.) | '
      WRITE (LUPRI,'(A,/)')
     *' +---------------------------------------------------------+ '
      STHELP = 'Dummy   '
      WRITE (LUPRI,*) ' Properties are read from unit number ', LU
      REWIND(LU)
      DO I = 1, NUOFPR
       READ(LU,
     *      '(I5,I3,I4,1X,A10,1X,E22.16,4(1X,A8),3(1X,E22.16),3I4)')
     *      IPRPC,ISYMIN,NORD,LABEL,PROP,
     *      LABX,LABY,LABZ,LABU,FRQY,FRQZ,FRQU,ISYEX,ISPEX,IEX
       IF (LABEL .EQ. 'THE_END   ') THEN
         WRITE(LUPRI,'(A)')
         WRITE(LUPRI,'(1X,A)') LABEL
         return
       END IF
       IF (I.NE.IPRPC) WRITE(LUPRI,*) ' PROBLEM OM PROPERTY LIST?'
       IF (IOPT.EQ.2) PROP = PROPIN(I)

       IF (STHELP.NE.LABEL) WRITE(LUPRI,'(/,1X,A10,A)') LABEL,
     *   ' properties :'
       STHELP = LABEL
       IF (NORD.EQ.0) THEN
         IF (LABX.EQ."Excita  ") THEN
          IF (IOPT.NE.2) THEN
            WRITE(LUPRI,'(1X,A,I3,I2,2(A,I1),A,I3,3(A,F18.10))')
     *      '#',I,ISYMIN,
     *       ' Sy: ',ISYEX,' Sp: ',ISPEX,' Exc:',IEX,
     *       ' E_tot:',FRQU,' w:',PROP,' E_grs:',FRQY
          ELSE
            WRITE(LUPRI,'(1X,A,I3,I2,2(A,I1),A,I3,1(A,F18.10))')
     *      '#',I,ISYMIN,
     *       ' Sy: ',ISYEX,' Sp: ',ISPEX,' Exc:',IEX,
     *       ' E_tot:',PROP
          ENDIF
         ELSE IF (LABX.EQ."Exctot  ") THEN
            WRITE(LUPRI,'(1X,A,I3,I2,A,A8,A,F18.10)')
     *      '#',I,ISYMIN,' ',LABX,' Excited state Energy = ',PROP
         ELSE
            WRITE(LUPRI,'(1X,A,I3,I2,A,A8,A,F18.10)')
     *      '#',I,ISYMIN,' ',LABX,' Ground state Energy = ',PROP
         ENDIF
       ENDIF
       IF (NORD.EQ.1) WRITE(LUPRI,'(1X,A,I3,I2,A,A8,A3,F20.12)')
     *  '#',I,ISYMIN,' <',LABX,'> = ',PROP
       IF (NORD.EQ.-11)
     *   WRITE(LUPRI,'(1X,A,I3,I2,A,A8,A,I1,I2,I3,F9.6,A,F20.12)')
     *  '#',I,ISYMIN,' <',LABX,'>(',ISYEX,ISPEX,IEX,FRQY,') = ',PROP
       IF (NORD.EQ.2) WRITE(LUPRI,
     * '(1X,A,I3,I2,A,A8,A1,A8,A3,F9.6,A,F20.12)')
     *  '#',I,ISYMIN,' <<',LABX,',',LABY,
     *  '>>(',FRQY,') =',PROP
       IF (NORD.EQ.3) WRITE(LUPRI,
     * '(1X,A,I3,I2,A,A8,A1,A8,A1,A8,A3,F9.6,A1,F9.6,A,F20.12)')
     *  '#',I,ISYMIN,' <<',LABX,',',LABY,',',LABZ,
     *  '>>(',FRQY,',',FRQZ,') =',PROP
       IF (NORD.EQ.4) WRITE(LUPRI,
     * '(1X,A,I3,I2,A,A8,A1,A8,A1,A8,A1,A8,A3,F9.6,A1,F9.6,A1,F9.6,
     *    A,F20.12)')
     *  '#',I,ISYMIN,' <<',LABX,',',LABY,',',LABZ,
     *   ',',LABU,
     *  '>>(',FRQY,',',FRQZ,',',FRQU,') =',PROP
       IF (NORD.EQ.-1) WRITE(LUPRI,
     *   '(1X,A,I3,I2,A,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     *   '#',I,ISYMIN,' |<O|',LABX,'|i(',FRQY,')>| =',PROP
       IF (NORD.EQ.-2) WRITE(LUPRI,
     *   '(1X,A,I3,I2,A,F9.6,A,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     *   '#',I,ISYMIN,' |<i(',FRQY,')|',LABX,'|f(',FRQZ,
     *   ')>| =',PROP
       IF (NORD.EQ.-20) WRITE(LUPRI,
     *   '(1X,A,I3,I2,A,A8,A1,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     *   '#',I,ISYMIN,' RES{<<',LABX,',',LABY,
     *   '>>(',FRQY,')} =',
     *   PROP,' ( ',SQRT(ABS(PROP)),')'
      ENDDO
      END
