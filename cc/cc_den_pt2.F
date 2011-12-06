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
C  /* Deck cc_den_pt2 */
      SUBROUTINE CC_DEN_PT2(ETAAI,ETAIJ,ETAAB,WORK,LWORK,
     &                      IOPT,LTSTEN,en2pt)
C
C     Written by S. Coriani, based on CC_DEN_PT
C     January 2002
C
C     Version: 1.0
C
C     Purpose: 
C     drive the calculation of the "pure d_pqrs(T)" contributions to the 
C     ^kappabar-eta_pq RHS of the orbital multipliers
C     LTSTEN = true, test densities via energy calculation
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION ETAAI(*), ETAIJ(*), ETAAB(*), WORK(LWORK)
      LOGICAL LTSTEN
#include <ccorb.h>
#include <ccisao.h>
#include <r12int.h>
#include <inftap.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccinftap.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ccfro.h>
C
      CHARACTER MODEL*10
      CHARACTER NAME1*8
      CHARACTER NAME2*8

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
C
      CALL QENTER('CC_DEN_PT2')
C
      CALL HEADER('Construct part of rhs for CCSD(T)-kappa-0',-1)
C
C-----------------------------------------
C     Initialization of timing parameters.
C-----------------------------------------
C
      TIMTOT = ZERO
      TIMTOT = SECOND()
      TIMDEN = ZERO
      TIMRES = ZERO
      TIRDAO = ZERO
      TIMHE2 = ZERO
      TIMONE = ZERO
      TIMONE = SECOND()


      IF (LTSTEN) EN2PT=ZERO
C
C----------------------------------------------------
C     Both zeta- and t-vectors are totally symmetric.
C----------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C----------------------------------------
C     Get CMO coefficients
C----------------------------------------
C
      KT1AM_PT = 1
      KEND0 = KT1AM_PT + NT1AM(ISYMOP)
      CALL DZERO(WORK(KT1AM_PT),NT1AM(ISYMOP))

      !KEND0 = 1

      KCMO  = KEND0
      KEND1 = KCMO  + NLAMDS
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation 1 CC_DEN_PT2')
      ENDIF
C
!      IF (FROIMP) THEN
C
C-------------------------------------------
C        Get the FULL MO coefficient matrix.
C-------------------------------------------
C
!         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
!      ELSE
!         !Sonia: get coefficients for (T) part and reorder
!         CALL CC_GET_CMO(WORK(KCMO))
!         CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
!      ENDIF
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
 
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE (LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C-------------------------------------
C  Two-electron part starts here.....
C-------------------------------------
C
      TIMONE = SECOND() - TIMONE
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0


      xnaigd = zero
      xnijgd = zero

      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRINT)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                        WORK(KODCL1),WORK(KODCL2),
     *                        WORK(KODBC1),WORK(KODBC2),
     *                        WORK(KRDBC1),WORK(KRDBC2),
     *                        WORK(KODPP1),WORK(KODPP2),
     *                        WORK(KRDPP1),WORK(KRDPP2),
     *                        WORK(KCCFB1),WORK(KINDXB),
     *                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMHE2 = TIMHE2   + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                CALL QUIT('Insufficient core in CC_DEN_PT2')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
C---------------------------------------------------------
C              Sonia
C              Work space allocation for the (T) densities
C              with third index backtransformed to gamma
C              All gammas together
C---------------------------------------------------------
C
               ISYDEN = ISYMD
C
               KD2IJG_PT = KEND1
               KD2AIG_PT = KD2IJG_PT + ND2IJG(ISYDEN)
               KD2IAG_PT = KD2AIG_PT + ND2AIG(ISYDEN)
               KD2ABG_PT = KD2IAG_PT + ND2AIG(ISYDEN)
               KEND2     = KD2ABG_PT + ND2ABG(ISYDEN)
               LWRK2     = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
                  CALL QUIT('Insufficient space for allocation '//
     &                      '2and1/2 in CC_DEN_PT2')
               ENDIF
C
C-------------------------------------------------------
C              Initialize 4 two electron density arrays.
C-------------------------------------------------------
C
               CALL DZERO(WORK(KD2IJG_PT),ND2IJG(ISYDEN))
               CALL DZERO(WORK(KD2AIG_PT),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2IAG_PT),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2ABG_PT),ND2ABG(ISYDEN))
C
C-------------------------------------------------------------------
C              Calculate the two electron density d(pq,gamma;delta).
C-------------------------------------------------------------------
C
               AUTIME = SECOND()
C
               if (.true.) then
               CALL CC_DEN2_PT(WORK(KD2IJG_PT),WORK(KD2AIG_PT),
     &                         WORK(KD2IAG_PT),WORK(KD2ABG_PT),
     &                         WORK(KCMO),1,
     &                         WORK(KEND2),LWRK2,
     &                         IDEL,ISYMD)
               end if
C

               if (.false.) then
                  xtest = ddot(ND2IJG(ISYDEN),WORK(KD2IJG_PT),1,
     &                                     WORK(KD2IJG_PT),1)
               write(lupri,*)'norm of ij,gamma;delta', 
     &                        xtest, isymd, idel
                  xnijgd = xnijgd + xtest
                  xtest = ddot(ND2AIG(ISYDEN),WORK(KD2AIG_PT),1,
     &                                     WORK(KD2AIG_PT),1)
               write(lupri,*)'norm of ai,gamma;delta', 
     &                        xtest, isymd, idel
                  xtest = ddot(ND2AIG(ISYDEN),WORK(KD2IAG_PT),1,
     &                                     WORK(KD2IAG_PT),1)
               write(lupri,*)'norm of ia,gamma;delta', 
     &                        xtest, isymd, idel

                  xnaigd = xnaigd + xtest

                  xtest = ddot(ND2ABG(ISYDEN),WORK(KD2ABG_PT),1,
     &                                     WORK(KD2ABG_PT),1)
               write(lupri,*)'norm of ab;gamma;delta', 
     &                        xtest, isymd, idel
               end if

               AUTIME = SECOND() - AUTIME
               TIMDEN = TIMDEN + AUTIME
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND2
               KEND3  = KXINT  + NDISAO(ISYDIS)
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND3
                  CALL QUIT('Insufficient space for allocation '//
     &                      '3 in CC_DEN_PT2')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               AUTIME = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     *                     WORK(KRECNR),DIRECT)
               AUTIME = SECOND() - AUTIME
               TIRDAO = TIRDAO + AUTIME
C
C------------------------------------------------------
C              Start loop over second AO-index (gamma).
C------------------------------------------------------
C
               AUTIME = SECOND()
C
               DO 130 ISYMG = 1,NSYM
C
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
                  DO 140 G = 1,NBAS(ISYMG)
C
                     KD2GIJ = KD2IJG_PT + ID2IJG(ISYMPQ,ISYMG)
     *                      + NMATIJ(ISYMPQ)*(G - 1) 
                     KD2GAI = KD2AIG_PT + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG_PT + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GAB = KD2ABG_PT + ID2ABG(ISYMPQ,ISYMG)
     *                      + NMATAB(ISYMPQ)*(G - 1)
C
C-----------------------------------------------
C                    Work space allocation four.
C-----------------------------------------------
C
                     KINTAO = KEND3
                     KD2AOB = KINTAO + N2BST(ISYMPQ)
                     KEND4  = KD2AOB + N2BST(ISYMPQ)
                     LWRK4  = LWORK  - KEND4
C
                     IF (LWRK4 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND4
                        CALL QUIT('Insufficient space in CC_DEN_PT2')
                     ENDIF
C
                     CALL DZERO(WORK(KD2AOB),N2BST(ISYMPQ))
C
C-------------------------------------------------------
C                    Square up AO-integral distribution.
C-------------------------------------------------------
C
                     KOFFIN = KXINT + IDSAOG(ISYMG,ISYDIS) 
     *                      + NNBST(ISYMPQ)*(G - 1) 
C
                     CALL CCSD_SYMSQ(WORK(KOFFIN),ISYMPQ,
     *                               WORK(KINTAO))

C
C---------------------------------------------------------
C
C---------------------------------------------------------
C
                     if (ltsten) then

                        CALL CC_DENAO(WORK(KD2AOB),ISYMPQ,
     *                                WORK(KD2GAI),WORK(KD2GAB),
     *                                WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                                WORK(KCMO),1,WORK(KCMO),1,
     *                                WORK(KEND4),LWRK4)
C
C---------------------------------------------------------------------
C                    Add relaxation terms to set up effective density.
C---------------------------------------------------------------------
C
!                        IF (IOPT .EQ. 3) THEN
C
!                           ICON = 1
!                           CALL CC_D2EFF(WORK(KD2AOB),G,ISYMG,IDEL,
!     *                                   ISYMD,WORK(KKABAO),
!     *                                   WORK(KDHFAO),ICON)
C
!                        ENDIF
C
C----------------------------------------------------------------------
C                    Calculate the 2 e- density contribution to E-ccsd.
C----------------------------------------------------------------------
C
                        EN2PT = EN2PT + HALF*DDOT(N2BST(ISYMPQ),
     *                                    WORK(KD2AOB),1,WORK(KINTAO),1)
C
                     end if
C
C-----------------------------------------------
C                    Work space allocation five.
C-----------------------------------------------
C
                        KIJINT = KEND4
                        KAIINT = KIJINT + NMATIJ(ISYMPQ)
                        KIAINT = KAIINT + NT1AM(ISYMPQ)
                        KABINT = KIAINT + NT1AM(ISYMPQ)
                        KEND5  = KABINT + NMATAB(ISYMPQ)
                        LWRK5  = LWORK  - KEND5
C
                        IF (LWRK5 .LT. 0) THEN
                           WRITE(LUPRI,*) 'Available:', LWORK
                           WRITE(LUPRI,*) 'Needed:', KEND5
                           CALL QUIT('Insufficient work space '//
     &                               'in CC_DEN_PT2')
                        ENDIF
C
C----------------------------------------------------------------
C                       Transform 2 integral indices to MO-basis.
C----------------------------------------------------------------
C
                        ISYM = ISYMPQ
                        CALL CCDINTMO(WORK(KIJINT),WORK(KIAINT),
     *                                WORK(KABINT),WORK(KAIINT),
     *                                WORK(KINTAO),WORK(KCMO),
     *                                WORK(KCMO),WORK(KEND5),
     *                                LWRK5,ISYM)
C
C-------------------------------------------------------------------
C                       Calculate 2 e- contribution to Zeta-Kappa-0.
C-------------------------------------------------------------------
C
                        ISYM = ISYMPQ

                        IF (IOPT.EQ.2) THEN

                          CALL CCPT_ETARS_2E(ETAIJ,ETAAB,
     &                                  WORK(KIJINT),WORK(KAIINT),
     &                                  WORK(KIAINT),WORK(KABINT),
     &                                  WORK(KD2GIJ),WORK(KD2GAI),
     &                                  WORK(KD2GIA),WORK(KD2GAB),
     &                                  WORK(KEND5),LWRK5,ISYM)
                          if (.false.) then
                          XETAIJ = DDOT(NMATIJ(1),ETAIJ(1),1,
     &                             ETAIJ(1),1)
               WRITE(LUPRI,*) 'Norm of eta_ij (2e) after (T)', XETAIJ
                          XETAAB = DDOT(NMATAB(1),ETAAB(1),1,
     &                             ETAAB(1),1)
               WRITE(LUPRI,*) 'Norm of eta_ab (2e) after (T)', XETAAB
                          end if

                        END IF

                        CALL CCPT_ETAAI_2E(ETAAI,
     &                                WORK(KIJINT),WORK(KAIINT),
     &                                WORK(KIAINT),WORK(KABINT),
     &                                WORK(KD2GIJ),WORK(KD2GAI),
     &                                WORK(KD2GIA),WORK(KD2GAB),
     &                                WORK(KEND5),LWRK5,
     &                                ISYM)
                          if (.false.) then
                         XETAAI = DDOT(NALLAI(1),ETAAI(1),1,
     &                                 ETAAI(1),1)
               WRITE(LUPRI,*) 'Norm of eta_ai (2e) after (T)', XETAAI
                          end if
C

  140                CONTINUE
  130             CONTINUE
C
                  AUTIME = SECOND() - AUTIME
                  TIMRES = TIMRES + AUTIME
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE

C
C-----------------------
C     Regain work space.
C-----------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C------------------------
C
C------------------------
C
      if (ltsten) then
         write(lupri,*)'CC_DEN_PT2--> EN2PT: ', EN2PT
      end if
C
C-----------------------
C     Write out timings.
C-----------------------
C
  99  TIMTOT = SECOND() - TIMTOT
C
      IF (IPRINT .GT. 3) THEN
         WRITE (LUPRI,*) ' '
         WRITE (LUPRI,*) 'Requested density dependent '//
     &        'quantities calculated'
         WRITE (LUPRI,*) 'Total time used in CC_DEN_PT2:', TIMTOT
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE (LUPRI,*) 'Time used for setting up 2 e- density:',TIMDEN
         WRITE (LUPRI,*) 'Time used for contraction with integrals:',
     &        TIMRES
         WRITE (LUPRI,*) 'Time used for reading 2 e- AO-integrals:',
     &        TIRDAO
         WRITE (LUPRI,*) 'Time used for calculating 2 e- AO-integrals:',
     *              TIMHE2
         WRITE (LUPRI,*) 'Time used for 1 e- density & intermediates:',
     *              TIMONE
      ENDIF
C
      CALL QEXIT('CC_DEN_PT2')
      RETURN
      END
