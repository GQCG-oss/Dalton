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
C              CC_LR       control the calculation of polarizabilities.
C                          Frequency dependent, without orbital 
C                          relaxation.
C
C              CC_LRESID   Control calculation of residues and 
C                          oscillator strenths.
C
C              CC_HYPPOL   calculation of quadratic response properties.
C
C              CC_2HYP     calculation of cubic response properties.
C
C              CC_3HYP     calculation of quartic response properties.
C
C              CC_SMCAL    calculation of second order moments
C
C              CC_TMCAL    calculation of third order moments
C
C              CC_MCDCAL   calculation of B terms in MCD
C
C     Loop over all models specified. The order is supposed 
C     to be the optimal one from a computational point of view.
C     Energy calculations are restarted from the preceeding
C     calculation. Excitation vectors is also saved
C     and the excitation energy calculation is restarted
C     from preceeding model solution as well. Similar for
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
C
      PARAMETER(NMODEL = 18)
      LOGICAL   LCCSAV(0:NMODEL), LTCEXC, LPTEXC
      LOGICAL   LTBAR0, LRSPSYM, LLEFTEIG
      LOGICAL   MOLGRD, LSORT
      DIMENSION WORK(LWORK)
      CHARACTER FILPQIM*(6), FILXYM*(6)
      PARAMETER (FILPQIM = 'CCPQ00', FILXYM = 'CC_XYM')
      CHARACTER*8 NAME(8)
      DATA NAME  /'CCAOIN_1','CCAOIN_2','CCAOIN_3','CCAOIN_4',
     *            'CCAOIN_5','CCAOIN_6','CCAOIN_7','CCAOIN_8'/
      COMMON/SORTIO/LUAOIN(8)

      INTEGER NCCEXCI2(8,3),ISYOFE2(8),ITROFE2(8)
C
      INTEGER ISYEXC2(MAXEXCI),IMULTE2(MAXEXCI)
C
      CALL QENTER('CC_DRV')
C
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
C
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
      CALL GPOPEN(LURES,'CC_RES','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LURES)
C
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
      DO IPRPC = 1, MXPRPC
         PRPC(IPRPC)   = 1000000.0
         CPRPC(IPRPC)  = 'DDDDDDDDDD'
         IOPRPC(IPRPC) = 0
         LAB0(IPRPC)   = 'DDDDDDDD'
         LAB1(IPRPC)   = 'DDDDDDDD'
         LAB2(IPRPC)   = 'DDDDDDDD'
         LAB3(IPRPC)   = 'DDDDDDDD'
         FRQ1(IPRPC)   = 1000000.0
         FRQ2(IPRPC)   = 1000000.0
         FRQ3(IPRPC)   = 1000000.0
         ISAVSY(IPRPC) = 0
      ENDDO
C
C-------------------------------------
C     Set logicals: 4,6,8 is obsolete.
C-------------------------------------
C
      LCCSAV(0)  = CIS
      LCCSAV(1)  = CCS
      LCCSAV(2)  = MP2
      LCCSAV(3)  = CCP2
      LCCSAV(4)  = MCC2
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
C
      IMOD   = 0
      LTCEXC = .FALSE.
C
      MMOD = 0 
      DO I = 0, NMODEL
        IF (LCCSAV(I)) MMOD = MMOD + 1
      END DO
 
C-----------------------------------------------------------------
C     set flags controlling AO integral sort for non-direct runs
C     and AO -> MO transformation to (ia|jb) etc. integrals 
C-----------------------------------------------------------------
      LSORT  = .TRUE.
      INTTR  = .TRUE.

C-----------------------------------------------------------------
C     flag for ground state zeroth-order Lagrange multiplier:
C-----------------------------------------------------------------
      LTBAR0 = (CCLRSD.OR.CCLR.OR.CCQR.OR.CCCR.OR.CC4R.OR.
     &          CC5R.OR.CCQR2R.OR.CCEXGR.OR.CCSM.OR.CCTM.OR.
     &          CCMCD.OR.CCEXLR.OR.MOLGRD.OR.CCSLV.OR.CCDERI)

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
     &           CC5R.OR.CCQR2R.OR.CCEXGR.OR.CCSM.OR.CCTM.OR.
     &           CCMCD.OR.CCEXLR)

      RSPIM   = RSPIM .OR. LRSPSYM .OR. LTBAR0

C-----------------------------------------------------------------
C     flag for left eigenvector calculation:
C-----------------------------------------------------------------
      LLEFTEIG = CCMCD .OR. CCSM .OR. CCTM .OR. CCEXLR

C-----------------------------------------------------------------
C     LPTEXC  - to control if any perturbative triples corrections 
C     to excitation energies and thus calculated left vectors.
C-----------------------------------------------------------------
C 
      LPTEXC = .FALSE.
      IF (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B) LPTEXC = .TRUE.
C
C----------------------
C     Loop over models.
C----------------------
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
         IF (I.EQ.4)  MCC2   = LCCSAV(4)
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
C        1. Start loop over different solvents.
C        2. Start loop for self-consistent CC solvent
C           solution.
C        Do at least one(which could be vacuum).
C        SLV98,OC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
         NSLV = MAX(1,NCCSLV)
C
C     Dummy assignment to avoid weird bug on SUN.
C     This should only be a temporary solution - VB
C
         NN = NSLV + NCCSLV
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

              CALL CCSD_ENERGY(LSORT,WORK,LWORK)

           END IF
C
C          switch off flags for sort and transformation
C          done in CCSD_ENERGY
C
           INTTR = .FALSE.
           LSORT = .FALSE.
C
           IF (CCP2) CCS = .TRUE.
C
C--------------------------------------
C          ****************************
C          ***** CC Lambda vector *****
C          ****************************
C--------------------------------------
C
           IF ((CCFOP.OR.LTBAR0)
     *       .AND.(.NOT.
     *       (CCSDT.OR.CCR3.OR.CCR1A.OR.CCR1B.OR.CCRT.OR.CCR3))) THEN
C
              CALL CC_FOP(WORK,LWORK)
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
                  CALL CC_EXCI(WORK,LWORK,'RE ')
               ENDIF
C
               IF (CCSPIC.AND.CCS) CALL CC_REDEIG(WORK,LWORK,OMPCCS)
               IF (CC2PIC.AND.CC2) CALL CC_REDEIG(WORK,LWORK,OMPCC2)
               IF (CCSDPI.AND.CCSD) CALL CC_REDEIG(WORK,LWORK,OMPCCSD)
C
               IF (LHTR.OR.CCLRSD.OR.CCQR2R.OR.LPTEXC.OR.CCP2.OR.
     *             LLEFTEIG.OR.CCEXGR) THEN
                 IF ((I .LT. 16).AND.(.NOT.LESKIP)) THEN
                   CALL CC_EXCI(WORK,LWORK,'LE ')
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
         IF (I.GT.9) GOTO 90 
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
C-----------------------------------------------------
C        *********************************************
C        ***** Excited state gradient properties *****
C        *********************************************
C-----------------------------------------------------
C
         IF (CCEXGR) THEN
            CALL CC_EXGR(WORK,LWORK)
         ENDIF 
C
C---------------------------------------------------------------
C        precalculate expectatation values & eff. Fock matrices:
C---------------------------------------------------------------
C
          I1DXORD = 0
          IOPREL  = 0
          CALL CC_GRAD2(I1DXORD,WORK,LWORK)
C         CALL CCEFFFOCK(I1DXORD,IOPREL,WORK,LWORK)

C
C-----------------------------------------------
C        ***************************************
C        ***** Solve CC-response equations.*****
C        ***************************************
C-----------------------------------------------
C
         IF (LRSPSYM) THEN 
            CALL CC_SOLRSP(WORK,LWORK)  
         END IF
C
C-----------------------------------------------
C        ***************************************
C        ***** CC linear response residues *****
C        ***************************************
C-----------------------------------------------
C
         IF (CCLRSD) THEN 
            CALL CC_LRESID(WORK,LWORK)  
         END IF
C
C--------------------------------------------------------
C        ************************************************
C        ***** CC quadratic reponse second residues *****
C        ************************************************
C--------------------------------------------------------
C
         IF (CCQR2R) THEN
            CALL CC_QR2RSD(WORK,LWORK)
         END IF
C
C-------------------------------------------------------
C        ***********************************************
C        ***** CC second-order transition moments  *****
C        ***********************************************
C-------------------------------------------------------
C
         IF (CCSM) THEN 
            CALL CC_SMCAL(WORK,LWORK)  
         ENDIF
C
C------------------------------------------------------
C        **********************************************
C        ***** CC third-order transition moments  *****
C        **********************************************
C------------------------------------------------------
C
         IF (CCTM) THEN 
            CALL CC_TMCAL(WORK,LWORK)  
         ENDIF
C
C-------------------------------------------------------
C        ***********************************************
C        ***** CC magnetic circular dichroism      *****
C        ***********************************************
C-------------------------------------------------------
C
         IF (CCMCD) THEN
            CALL CC_MCDCAL(WORK,LWORK)
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
            IF (((.NOT. MP2).OR.CCP2).AND.(.NOT.CCD)) THEN
               IF (OLDLR) THEN
                 CALL CC_LR(WORK,LWORK)
               ELSE
                 CALL CC_SOP(WORK,LWORK)
               END IF
            ELSE
               WRITE(LUPRI,'(/,1X,A)') 
     &              ' No MP2 and CCD Polarizabilities '
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
           IF (CCS .or. CC2 .or. CCSD) THEN
             CALL CC_CAUCHY(WORK,LWORK)
           ELSE
             WRITE(LUPRI,'(/,1X,2A)') 'Cauchy moments',
     &        ' are not yet implemented for the requested model'
           ENDIF
         ENDIF
C
C--------------------------------------------
C        ******************************************
C        ***** CC first Hyperpolarizabilities *****
C        ******************************************
C--------------------------------------------
C
         IF ( CCQR ) THEN
C
            IF (CCS .or. CC2 .or. CCSD) THEN
               CALL CC_HYPPOL(WORK,LWORK)
            ELSE
              WRITE(LUPRI,'(/,1X,A)') 'requested model not '//
     &              'yet implemented'
            ENDIF
C
         ENDIF
C
C--------------------------------------------
C        *******************************************
C        ***** CC second Hyperpolarizabilities *****
C        *******************************************
C--------------------------------------------
C
         IF ( CCCR ) THEN
C
            IF (CCS .or. CC2 .or. CCSD) THEN
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
            IF (CCS .or. CC2 .or. CCSD) THEN
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
            IF (CCS .or. CC2 .or. CCSD) THEN
               CALL CC_4HYP(WORK,LWORK)
            ELSE
              WRITE(LUPRI,'(/,1X,A)') 
     &              'requested model not yet implemented'
            ENDIF
C
         ENDIF
C
C----------------------------------------------------
C        ********************************************
C        ***** CC excited state linear response *****
C        ********************************************
C----------------------------------------------------
C
         IF (CCEXLR) THEN 
            CALL CC_EXLR(WORK,LWORK)  
         END IF
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
C              CALL CC_BTST(WORK,LWORK)
C              CALL CC_CTST(WORK,LWORK)
C              CALL CC_DTST(WORK,LWORK)
C
C              ! the following calls to CC_PQIM and CC_XYIM should be
C              ! moved inside to the CC_FTST routine...
C              IF (.NOT.(CCS.OR.MP2.OR.CC2)) THEN
C                 CALL CC_PQIM('L0',0,'R0',0,FILPQIM,WORK,LWORK)
C              END IF
C              IF (.NOT.CCS) THEN
C                 CALL CC_XYMIM('L0',0,'R0',0,FILXYM,WORK,LWORK)
C              END IF
C              CALL CC_FTST(WORK,LWORK)
C
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
      MCC2   = LCCSAV(4)
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
      CALL GPCLOSE(LUIAJB,'KEEP')
      IF (LUAORC(0) .GT. 0) CALL GPCLOSE(LUAORC,'DELETE')
      IF (LUINTR.GT.0) CALL GPCLOSE(LUINTR,'DELETE')
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
      CHARACTER*80  STR,STHELP*10
C
#include <codata.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cclr.h>
#include <prpc.h>
#include <ccpack.h>
#include <ccinftap.h>
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
C
C        Save the list of properties.
C
         NPRPCO = NPRPC
         DO I = 1, NPRPC
            PRPCO(I)  = PRPC(I)
            CPRPCO(I) = CPRPC(I)
            IOPRPCO(I)= IOPRPC(I)
            LAB0O(I)  = LAB0(I)
            LAB1O(I)  = LAB1(I)
            LAB2O(I)  = LAB2(I)
            LAB3O(I)  = LAB3(I)
            FRQ1O(I)  = FRQ1(I)
            FRQ2O(I)  = FRQ2(I)
            FRQ3O(I)  = FRQ3(I)
            ISAVSYO(I)= ISAVSY(I)
         ENDDO
         FIRSTCALL = .FALSE.
         CALL CC_WPRPC
      ELSE
C
C        Reorder the list of properties to match the firstcall order
C
         CALL CC_PRPCREORDER
      ENDIF
C
      IF (NPRPC.GT.0) THEN
       WRITE (LUPRI,'(A,/)') '  '
       WRITE (LUPRI,'(A)')
     * ' +---------------------------------------------------------+ '
       WRITE (LUPRI,'(A)')
     * ' | List of selected calculated molecular properties (a.u.) | '
       WRITE (LUPRI,'(A,/)')
     * ' +---------------------------------------------------------+ '
       STHELP = 'Dummy'
       DO I = 1, NPRPC
        NORD = IOPRPC(I)
        IF (STHELP.NE.CPRPC(I)) WRITE(LUPRI,'(/,1X,A10,A)') CPRPC(I),
     *    ' properties :'
        STHELP = CPRPC(I)
        IF (NORD.EQ.1) WRITE(LUPRI,'(1X,A,I3,I2,A,A8,A3,F20.12)')
     *   '#',I,ISAVSY(I),' <',LAB0(I),'> = ',PRPC(I)
        IF (NORD.EQ.2) WRITE(LUPRI,
     *  '(1X,A,I3,I2,A,A8,A1,A8,A3,F9.6,A,F20.12)')
     *   '#',I,ISAVSY(I),' <<',LAB0(I),',',LAB1(I),
     *   '>>(',FRQ1(I),') =',PRPC(I)
        IF (NORD.EQ.-1) WRITE(LUPRI,
     *    '(1X,A,I3,I2,A,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     *    '#',I,ISAVSY(I),' |<O|',LAB0(I),'|i(',FRQ1(I),')>| =',PRPC(I)
        IF (NORD.EQ.-3) WRITE(LUPRI,
     *    '(1X,A,I3,I2,A,A8,A1,A8,A3,F9.6,A,F15.8,A,F12.8,A)')
     *    '#',I,ISAVSY(I),' RES{<<',LAB0(I),',',LAB1(I),
     *    '>>(',FRQ1(I),')} =',
     *    PRPC(I),' ( ',SQRT(ABS(PRPC(I))),')'
       ENDDO
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
      MCC2   = .FALSE.
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
      LUTOC   = 0
      LU3VI   = 0
      LU3VI2  = 0
      LU3VI3  = 0
      LU3VI4  = 0
      LU3SRT  = 0
      LUDELD  = 0
      LUDKBC  = 0
      LUCKJD  = 0
      LUIAJB  = 0
      LUBFDN  = 0
      LURES   = 0
      LUCSOL  = 0
      LUCTLM  = 0
      LU3SRT2 = 0
      LU3SRT3 = 0
      LUCKJD2 = 0
      LUCKJD3 = 0
      LUDELD2 = 0
      LUDELD3 = 0
      LUDELD4 = 0
      LUDELD5 = 0
      LUDKBC2 = 0
      LUDKBC3 = 0
      LUDKBC4 = 0
      LUDKBC5 = 0
C
      FNTOC   = 'CCSDT_OC'
      FN3VI   = 'CC3_VI'
      FN3VI2  = 'CC3_VI12'
      FN3VI3  = 'CC3_VI13'
      FN3VI4  = 'CC3_VI14'
      FN3SRT  = 'CC3_SORT'
      FNDELD  = 'CKDELD'
      FNDKBC  = 'DKBC'
      FNCKJD  = 'CKJDEL'
      FN3SRT2 = 'CC3_SORT2'
      FN3SRT3 = 'CC3_SORT3'
      FNCKJD2 = 'CKJDEL2'
      FNCKJD3 = 'CKJDEL3'
      FNDELD2 = 'CKDELD2'
      FNDELD3 = 'CKDELD3'
      FNDELD4 = 'CKDELD4'
      FNDELD5 = 'CKDELD5'
      FNDKBC2 = 'DKBC2'
      FNDKBC3 = 'DKBC3'
      FNDKBC4 = 'DKBC4'
      FNDKBC5 = 'DKBC5'
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
      LUPRP = 0
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
C
      DIMENSION FRQ12(MXPRPC),FRQ22(MXPRPC),FRQ32(MXPRPC)
      DIMENSION PRPC2(MXPRPC)
      CHARACTER*8  LAB02(MXPRPC), LAB12(MXPRPC), LAB22(MXPRPC),       
     *             LAB32(MXPRPC)
      CHARACTER*10 CPRPC2(MXPRPC)
      DIMENSION IOPRPC2(MXPRPC), ISAVSY2(MXPRPC) 
      PARAMETER (TOLPRPFD=1.0D-10,TOLEX=1.0D-02)
C
      CALL QENTER('CC_WPRPCREORDER')
C
      IF (NPRPC.NE.NPRPCO) CALL QUIT('Mistmatch in WPRPC')
      DO I = 1, NPRPC
         PRPC2(I)  = PRPC(I)
         CPRPC2(I) = CPRPC(I)
         IOPRPC2(I)= IOPRPC(I)
         LAB02(I)  = LAB0(I)
         LAB12(I)  = LAB1(I)
         LAB22(I)  = LAB2(I)
         LAB32(I)  = LAB3(I)
         FRQ12(I)  = FRQ1(I)
         FRQ22(I)  = FRQ2(I)
         FRQ32(I)  = FRQ3(I)
         ISAVSY2(I)= ISAVSY(I) 
      ENDDO
C
      DO 100 I = 1, NPRPC
       DO 200 J = 1, NPRPC
        IF (CPRPC2(J).EQ.CPRPCO(I))  THEN
         IF (IOPRPC2(J).EQ.IOPRPCO(I)) THEN 
          IF ((LAB02(J).EQ.LAB0O(I)).AND.
     *        (LAB12(J).EQ.LAB1O(I)).AND.
     *        (LAB22(J).EQ.LAB2O(I)).AND.
     *        (LAB32(J).EQ.LAB3O(I)))  THEN
C Take care with excitation energies !!!!
           IF (((IOPRPCO(I).GT.0).AND.
     *          ((ABS(FRQ12(J)-FRQ1O(I)).LT.TOLPRPFD).AND.
     *           (ABS(FRQ22(J)-FRQ2O(I)).LT.TOLPRPFD).AND.
     *           (ABS(FRQ32(J)-FRQ3O(I)).LT.TOLPRPFD))).OR.
     *         ((IOPRPCO(I).LT.0).AND.     
     *          ((ABS(FRQ12(J)-FRQ1O(I)).LT.TOLEX).AND.
     *           (ABS(FRQ22(J)-FRQ2O(I)).LT.TOLEX).AND.
     *           (ABS(FRQ32(J)-FRQ3O(I)).LT.TOLEX)))) THEN
              PRPC(I)  = PRPC2(J)
              CPRPC(I) = CPRPC2(J)
              IOPRPC(I)= IOPRPC2(J)
              LAB0(I)  = LAB02(J)
              LAB1(I)  = LAB12(J)
              LAB2(I)  = LAB22(J)
              LAB3(I)  = LAB32(J)
              FRQ1(I)  = FRQ12(J)
              FRQ2(I)  = FRQ22(J)
              FRQ3(I)  = FRQ32(J)
              ISAVSY(I) = ISAVSY2(J) 
              GOTO 100
           ENDIF
          ENDIF
         ENDIF
        ENDIF
 200   CONTINUE
 100  CONTINUE
C
C Loop old list
C
      CALL QEXIT('CC_WPRPCREORDER')
C
      END

