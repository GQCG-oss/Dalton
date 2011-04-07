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
*---------------------------------------------------------------------*
c/* Deck CC_CMAT */
*=====================================================================*
      SUBROUTINE CC_CMAT( ICTRAN,  NCTRAN, LISTA, LISTB, LISTC,
     &                    IOPTRES, FILCMA, ICDOTS,CCONS, MXVEC,
     &                    WORK,    LWORK                       )
*---------------------------------------------------------------------*
*
*    Purpose: AO-direct calculation of a linear transformation of three
*             CC amplitude vectors, T^A, T^B and T^C, with the coupled
*             cluster C matrix (second derivative of the CC jacobian 
*             with respect to the amplitudes)
*          
*             The linear transformations are calculated for a list
*             of T^A, T^B and T^C vectors: 
*
*                LISTA       -- type of T^A vectors
*                LISTB       -- type of T^B vectors
*                LISTC       -- type of T^C vectors
*                ICTRAN(1,*) -- indeces of T^A vectors
*                ICTRAN(2,*) -- indeces of T^B vectors
*                ICTRAN(3,*) -- indeces of T^C vectors
*                ICTRAN(4,*) -- indeces or addresses of result vectors
*                NCTRAN      -- number of requested transformations
*                FILCMA      -- file name / list type of result vectors
*                               or list type of vectors to be dotted on
*                ICDOTS      -- indeces of vectors to be dotted on
*                CCONS       -- contains the dot products on return
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  all result vectors are written to a direct
*                          access file, FILCMA is used as file name
*                          the start addresses of the vectors are
*                          returned in ICTRAN(4,*)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          ICTRAN(4,*). N.B.: if WORK is not large
*                          enough iopt is automatically reset to 0!!!
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, FILCMA is used
*                          as list type and ICTRAN(4,*) as list index
*                          NOTE that ICTRAN(4,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WRRSP, FILCMA is used
*                          as list type and ICTRAN(4,*) as list index
*                          NOTE that ICTRAN(4,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by FILCMA and the indeces from ICDOTS
*                          the result of the dot products is returned
*                          in the CCONS array
*
*
*     Written by Christof Haettig, april-june 1997.
*
*     included CC-R12: Christian Neiss, nov. 2005
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "mxcent.h"
#include "ccorb.h"
#include "cciccset.h"
#include "cbieri.h"
#include "distcl.h"
#include "iratdef.h"
#include "eribuf.h"
#include "ccisao.h"
#include "ccfield.h"
#include "aovec.h"
#include "blocks.h"
#include "r12int.h"
#include "ccr12int.h"

* local parameters:
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CC_CMAT> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
      INTEGER ISYM0
      PARAMETER( ISYM0 = 1 ) ! symmetry of the reference state
      INTEGER ISYOVOV
      PARAMETER( ISYOVOV = 1 ) ! symmetry of (ia|jb) integrals 
      
      INTEGER LUCMAT

      CHARACTER*(*) LISTA, LISTB, LISTC, FILCMA
      INTEGER IOPTRES
      INTEGER NCTRAN, MXVEC, LWORK
      INTEGER ICTRAN(4,NCTRAN)
      INTEGER ICDOTS(MXVEC,NCTRAN)

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL ZERO, ONE, TWO, HALF
      REAL CCONS(MXVEC,NCTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      DOUBLE PRECISION CCONS(MXVEC,NCTRAN)
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, HALF = 0.5d0)


      CHARACTER*(10) MODEL, MODELW
      CHARACTER*(3) LISTXA, LISTXB, LISTXC
      CHARACTER*(8) CBAFIL, DBAFIL
      INTEGER INDEXA(MXCORB)
      INTEGER IOPTW, IOPT, ITRAN, IADRTH
      INTEGER LENALL, LEN, IOFFCD, ICYCLE, ILLL, IERROR
      INTEGER KEND0, LEND0, KENDE2, LENDE2, KENDE3, LENDE3, KEND, LEND
      INTEGER KLIAJB, KT2AMPA, KT2AMPC, KT2AMPB, KTHETA1, KTHETA2
      INTEGER KT1AMPA, KT1AMPB, KT1AMPC, KFOCKA, KFOCKB, KFOCKC
      INTEGER KFABVV, KFABOO, KFBCVV, KFBCOO, KFCAVV, KFCAOO, KTHETA0
      INTEGER KX4O, KXIAJB, KKINTC, KCBAR, KDBAR, LUCBAR, LUDBAR
      INTEGER IDLSTA, IDLSTB, IDLSTC, ISYRES, ISYMAB, ISYCBAR, ISYDBAR
      INTEGER ISYMTA, ISYMTB, ISYMTC, ISYMFA, ISYMFB, ISYMFC, ISYX4O
      INTEGER ISYFAB, ISYFBC, ISYFCA, ISYMD1, NTOT, ISYMKC
      INTEGER NTOSYM, KCCFB1, KINDXB, KFREE, LFREE, KENDSV, LWRKSV
      INTEGER KODCL1, KODBC1, KRDBC1, KODPP1, KRDPP1, KRECNR
      INTEGER KODCL2, KODBC2, KRDBC2, KODPP2, KRDPP2, NUMDIS
      INTEGER IDEL2, ISYDEL, IDEL, KXINT
      INTEGER KT1AMP0, KLAMDH0, KLAMDP0, KLAMDHA, KLAMDPA
      INTEGER KLAMDHB, KLAMDPB, KLAMDPC, KLAMDHC
      INTEGER KENDF1, LENDF1, KENDF2, LENDF2, KEND1, LEND1
      INTEGER KEND2, LEND2, KEND3, LEND3, KENDA3, LENDA3
      INTEGER KENDB3, LENDB3, IOPTTCME
      INTEGER IDUM, IAN, KTR12AM, KVINT, LUNIT
      
#if defined (SYS_CRAY) 
      REAL XNORM
#else 
      DOUBLE PRECISION XNORM
#endif 




* external functions:
      INTEGER ILSTSYM

#if defined (SYS_CRAY) 
      REAL DDOT
#else 
      DOUBLE PRECISION DDOT
#endif 
  
      CALL QENTER('CC_CMAT')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_CMAT')
        IF (DIRECT) WRITE(LUPRI,'(/1X,A)') 'AO direct transformation'
        WRITE (LUPRI,*) MSGDBG,'LISTA : ',LISTA(1:2)
        WRITE (LUPRI,*) MSGDBG,'LISTB : ',LISTB(1:2)
        WRITE (LUPRI,*) MSGDBG,'LISTC : ',LISTC(1:2)
        WRITE (LUPRI,*) MSGDBG,'FILCMA: ',FILCMA
        WRITE (LUPRI,*) MSGDBG,'NCTRAN: ',NCTRAN
        WRITE (LUPRI,*) MSGDBG,'IOPTRES:',IOPTRES
        WRITE (LUPRI,*) MSGDBG,'MXVEC:',MXVEC
        WRITE (LUPRI,*) MSGDBG,'LWORK:',LWORK
        CALL FLSHFO(LUPRI)
      END IF
      
      IF (CCSDT) THEN
        WRITE(LUPRI,'(/1x,a)') 'C matrix transformations not '
     &          //'implemented for triples yet...'
        CALL QUIT('Triples not implemented for C '//
     &            'matrix transformations')
      END IF

      IF ( .not. (CCS .or. CC2 .or. CCSD) ) THEN
        WRITE(LUPRI,'(/1x,a)') 'CC_CMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_CMAT...'
        CALL QUIT('Unknown CC method in CC_CMAT.')
      END IF

      IF (    LISTA(1:1).NE.'R' 
     &   .OR. LISTB(1:1).NE.'R'
     &   .OR. LISTC(1:1).NE.'R' ) THEN
        WRITE(LUPRI,*) 'LISTA, LISTB and LISTC must refer to',
     &                    ' t-amplituded vectors in CC_CMAT.'
        CALL QUIT('Illegal LISTA or LISTB or LISTC in CC_CMAT.')
      END IF

      IF (ISYMOP .NE. 1) THEN
        WRITE(LUPRI,*) 'ISYMOP = ',ISYMOP
        WRITE(LUPRI,*) 'CC_CMAT is not implemented for ISYMOP.NE.1'
        CALL QUIT('CC_CMAT is not implemented for ISYMOP.NE.1')
      END IF

C     IF (NCTRAN .GT. MAXSIM) THEN
C       WRITE(LUPRI,*) 'NCTRAN = ', NCTRAN
C       WRITE(LUPRI,*) 'MAXSIM = ', MAXSIM
C       WRITE(LUPRI,*) 'number of requested transformation is larger'
C       WRITE(LUPRI,*) 'than the maximum number of allowed ',
C    &                 'simultaneous transformation.'
C       CALL QUIT('Error in CC_CMAT: NCTRAN is larger than MAXSIM.')
C     END IF

* check return option for the result vectors:
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN

        LUCMAT = -1
        CALL WOPEN2(LUCMAT,FILCMA,64,0)

      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
        CONTINUE
      ELSE IF (IOPTRES.EQ.5) THEN
        IF (MXVEC*NCTRAN.NE.0) CALL DZERO(CCONS,MXVEC*NCTRAN)  
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_CMAT.')
      END IF

* construct 'MODEL' string for CC_WRRSP routine and set write option:
      IF (CCS) THEN
         MODELW = 'CCS       '
         IOPTW  = 1
      ELSE IF (CC2) THEN
         MODELW = 'CC2       '
         IOPTW  = 3
      ELSE IF (CCSD) THEN
         MODELW = 'CCSD      '
         IOPTW  = 3
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_CMAT.')
      END IF



* start of scratch space for the following:
      KEND0 = 1
      LEND0 = LWORK

*=====================================================================*
* Calculate J term and N^5 terms E1 and E2 in a loop over all required
* C matrix transformations
*
* the N^5 terms H, G, I vanish for the C matrix 
*
* All N^5 terms drop out for CCS, the E1, E2 vanish also for CC2.
*=====================================================================*

* memory allocation:
      KLIAJB = KEND0
      KENDE2 = KLIAJB + NT2AM(ISYOVOV)
      LENDE2 = LWORK - KENDE2

      IF (LENDE2.LT.0) THEN
        CALL QUIT('Insufficient memory in CC_CMAT.(1)')
      END IF

* read packed (ia|jb) integrals and calculate L(ia|jb) in place: 
      Call CCG_RDIAJB(WORK(KLIAJB),NT2AM(ISYOVOV))
      
      IOPTTCME = 1
      Call CCSD_TCMEPK(WORK(KLIAJB),ONE,ISYOVOV,IOPTTCME)

*---------------------------------------------------------------------*
* start loop over all C matrix transformations
*---------------------------------------------------------------------*
      IADRTH = 1

      DO ITRAN = 1, NCTRAN
        IDLSTA = ICTRAN(1,ITRAN)
        IDLSTB = ICTRAN(2,ITRAN)
        IDLSTC = ICTRAN(3,ITRAN)

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG, 'ITRAN :',ITRAN
          WRITE (LUPRI,*) MSGDBG, 'LISTA/IDLSTA:',LISTA(1:2),IDLSTA
          WRITE (LUPRI,*) MSGDBG, 'LISTB/IDLSTB:',LISTB(1:2),IDLSTB
          WRITE (LUPRI,*) MSGDBG, 'LISTC/IDLSTC:',LISTC(1:2),IDLSTC
          CALL FLSHFO(LUPRI)
        END IF
*---------------------------------------------------------------------*
* set symmetries and allocate memory:
*---------------------------------------------------------------------*
        ISYMTA = ILSTSYM(LISTA,IDLSTA)
        ISYMTB = ILSTSYM(LISTB,IDLSTB)
        ISYMTC = ILSTSYM(LISTC,IDLSTC)

        ISYMFA = MULD2H(ISYOVOV,ISYMTA)
        ISYMFB = MULD2H(ISYOVOV,ISYMTB)
        ISYMFC = MULD2H(ISYOVOV,ISYMTC)

        ISYFAB = MULD2H(ISYMFA,ISYMTB)
        ISYFBC = MULD2H(ISYMFB,ISYMTC)
        ISYFCA = MULD2H(ISYMFC,ISYMTA)

        ISYRES = MULD2H(ISYFAB,ISYMTC)

* allocate memory for double excitation response vectors:
* (overlaps with memory for two-electron integrals)
        KENDE3  = KLIAJB + NT2AM(ISYOVOV)
        IF (.NOT.CCS) THEN
          KT2AMPA = KLIAJB
          KT2AMPB = KLIAJB
          KT2AMPC = KLIAJB

          KENDE3  = MAX(KENDE3,KLIAJB + NT2SQ(ISYMTA))
          KENDE3  = MAX(KENDE3,KLIAJB + NT2SQ(ISYMTB))
          KENDE3  = MAX(KENDE3,KLIAJB + NT2SQ(ISYMTC))
        END IF

* allocate memory for the result vector:
        KTHETA1 = KENDE3
        KENDE3  = KTHETA1 + NT1AM(ISYRES)

        IF (.NOT.CCS) THEN
          KTHETA2 = KENDE3
          KENDE3  = KTHETA2 + NT2AM(ISYRES)
        END IF

        KT1AMPA = KENDE3
        KT1AMPB = KT1AMPA + NT1AM(ISYMTA)
        KT1AMPC = KT1AMPB + NT1AM(ISYMTB)
        KFOCKA  = KT1AMPC + NT1AM(ISYMTC)
        KFOCKB  = KFOCKA  + NT1AM(ISYMFA)
        KFOCKC  = KFOCKB  + NT1AM(ISYMFB)
        KFABVV  = KFOCKC  + NT1AM(ISYMFC)
        KFABOO  = KFABVV  + NMATAB(ISYFAB)
        KFBCVV  = KFABOO  + NMATIJ(ISYFAB)
        KFBCOO  = KFBCVV  + NMATAB(ISYFBC)
        KFCAVV  = KFBCOO  + NMATIJ(ISYFBC)
        KFCAOO  = KFCAVV  + NMATAB(ISYFCA)
        KENDE3  = KFCAOO  + NMATIJ(ISYFCA)
        LENDE3  = LWORK - KENDE3

        IF (LENDE3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_CMAT.(2)')
        END IF

*---------------------------------------------------------------------*
* read single excitation part of the response vectors:
*---------------------------------------------------------------------*
* read A response amplitudes:
        IOPT = 1
        Call CC_RDRSP(LISTA,IDLSTA,ISYMTA,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KDUM))

* read B response amplitudes:
        IOPT = 1
        Call CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     &                WORK(KT1AMPB),WORK(KDUM))

* read C response amplitudes:
        IOPT = 1
        Call CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                WORK(KT1AMPC),WORK(KDUM))

*---------------------------------------------------------------------*
* calculate occupied/virtual part of first order Fock matrices:
*---------------------------------------------------------------------*
* calculate tF^A_ia:
        Call CCG_LXD(WORK(KFOCKA), ISYMFA, WORK(KT1AMPA), ISYMTA, 
     &               WORK(KLIAJB), ISYOVOV, 0)
      
* calculate tF^B_ia:
        Call CCG_LXD(WORK(KFOCKB), ISYMFB, WORK(KT1AMPB), ISYMTB,
     &               WORK(KLIAJB), ISYOVOV, 0)

* calculate tF^C_ia:
        Call CCG_LXD(WORK(KFOCKC), ISYMFC, WORK(KT1AMPC), ISYMTC, 
     &               WORK(KLIAJB), ISYOVOV, 0)
      
*---------------------------------------------------------------------*
* calculate double one-index transformed Fock matrices:
*---------------------------------------------------------------------*
* calculate occ/occ block for FAB:
        Call CCG_1ITROO(WORK(KFABOO),ISYFAB,
     &                  WORK(KFOCKB),ISYMFB,
     &                  WORK(KT1AMPA),ISYMTA    )

        IF (LENDE3.LT.NMATIJ(ISYFAB)) 
     &       CALL QUIT('Out of memory in CC_CMAT.')

        Call CCG_1ITROO(WORK(KENDE3),ISYFAB,
     &                  WORK(KFOCKA),ISYMFA,
     &                  WORK(KT1AMPB),ISYMTB    )

        Call DAXPY(NMATIJ(ISYFAB),ONE,WORK(KENDE3),1,WORK(KFABOO),1)

* calculate vir/vir block for FAB:
        Call CCG_1ITRVV(WORK(KFABVV),ISYFAB,
     &                  WORK(KFOCKB),ISYMFB,
     &                  WORK(KT1AMPA),ISYMTA  )

        IF (LENDE3.LT.NMATAB(ISYFAB)) 
     &       CALL QUIT('Out of memory in CC_CMAT.')

        Call CCG_1ITRVV(WORK(KENDE3),ISYFAB,
     &                  WORK(KFOCKA),ISYMFA,
     &                  WORK(KT1AMPB),ISYMTB  )

        Call DAXPY(NMATAB(ISYFAB),ONE,WORK(KENDE3),1,WORK(KFABVV),1)

* calculate occ/occ block for FBC:
        Call CCG_1ITROO(WORK(KFBCOO),ISYFBC,
     &                  WORK(KFOCKB),ISYMFB,
     &                  WORK(KT1AMPC),ISYMTC    )

        IF (LENDE3.LT.NMATIJ(ISYFBC)) 
     &       CALL QUIT('Out of memory in CC_CMAT.')

        Call CCG_1ITROO(WORK(KENDE3),ISYFBC,
     &                  WORK(KFOCKC),ISYMFC,
     &                  WORK(KT1AMPB),ISYMTB    )

        Call DAXPY(NMATIJ(ISYFBC),ONE,WORK(KENDE3),1,WORK(KFBCOO),1)

* calculate vir/vir block for FBC:
        Call CCG_1ITRVV(WORK(KFBCVV),ISYFBC,
     &                  WORK(KFOCKB),ISYMFB,
     &                  WORK(KT1AMPC),ISYMTC  )

        IF (LENDE3.LT.NMATAB(ISYFBC)) 
     &       CALL QUIT('Out of memory in CC_CMAT.')

        Call CCG_1ITRVV(WORK(KENDE3),ISYFBC,
     &                  WORK(KFOCKC),ISYMFC,
     &                  WORK(KT1AMPB),ISYMTB  )

        Call DAXPY(NMATAB(ISYFBC),ONE,WORK(KENDE3),1,WORK(KFBCVV),1)

* calculate occ/occ block for FCA:
        Call CCG_1ITROO(WORK(KFCAOO),ISYFCA,
     &                  WORK(KFOCKC),ISYMFC,
     &                  WORK(KT1AMPA),ISYMTA    )

        IF (LENDE3.LT.NMATIJ(ISYFCA))
     &       CALL QUIT('Out of memory in CC_CMAT.')

        Call CCG_1ITROO(WORK(KENDE3),ISYFCA,
     &                  WORK(KFOCKA),ISYMFA,
     &                  WORK(KT1AMPC),ISYMTC    )

        Call DAXPY(NMATIJ(ISYFCA),ONE,WORK(KENDE3),1,WORK(KFCAOO),1)

* calculate vir/vir block for FCA:
        Call CCG_1ITRVV(WORK(KFCAVV),ISYFCA,
     &                  WORK(KFOCKC),ISYMFC,
     &                  WORK(KT1AMPA),ISYMTA  )

        IF (LENDE3.LT.NMATAB(ISYFCA))
     &       CALL QUIT('Out of memory in CC_CMAT.')

        Call CCG_1ITRVV(WORK(KENDE3),ISYFCA,
     &                  WORK(KFOCKA),ISYMFA,
     &                  WORK(KT1AMPC),ISYMTC  )

        Call DAXPY(NMATAB(ISYFCA),ONE,WORK(KENDE3),1,WORK(KFCAVV),1)

*---------------------------------------------------------------------*
* J term
*---------------------------------------------------------------------*
        IF (LENDE3.LT.NT1AM(ISYRES))
     &       CALL QUIT('Out of memory in CC_CMAT.')

* initialize single-excitation part of the result vector:
        CALL DZERO(WORK(KTHETA1),NT1AM(ISYRES))

* calculate triple one-index transformed contribution T^C x F^AB:
        IOPT  = 1
        Call CCG_1ITRVO( WORK(KENDE3), ISYRES,
     &                   WORK(KFABOO), WORK(KFABVV), ISYFAB,
     &                   WORK(KT1AMPC), ISYMTC, IOPT          )

        CALL DAXPY(NT1AM(ISYRES),HALF,WORK(KENDE3),1,WORK(KTHETA1),1)


* calculate triple one-index transformed contribution T^A x F^BC:
        IOPT  = 1
        Call CCG_1ITRVO( WORK(KENDE3), ISYRES,
     &                   WORK(KFBCOO), WORK(KFBCVV), ISYFBC,
     &                   WORK(KT1AMPA), ISYMTA, IOPT          )

        CALL DAXPY(NT1AM(ISYRES),HALF,WORK(KENDE3),1,WORK(KTHETA1),1)


* calculate triple one-index transformed contribution T^B x F^CA:
        IOPT  = 1
        Call CCG_1ITRVO( WORK(KENDE3), ISYRES,
     &                   WORK(KFCAOO), WORK(KFCAVV), ISYFCA,
     &                   WORK(KT1AMPB), ISYMTB, IOPT          )

        CALL DAXPY(NT1AM(ISYRES),HALF,WORK(KENDE3),1,WORK(KTHETA1),1)

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG,'THETA1 after J term:'
          CALL CC_PRP(WORK(KTHETA1),WORK(KDUM),ISYRES,1,0)
          CALL FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
* initialize double-excitation part of the result vector:
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN
          CALL DZERO(WORK(KTHETA2),NT2AM(ISYRES))
        END IF
 
*---------------------------------------------------------------------*
* E1 and E2 terms 
*---------------------------------------------------------------------*
        IF (.NOT. (CCS.OR.CC2.OR.CCSTST) ) THEN

* check available scratch space:
          IF (     LENDE3 .LT. NT2AM(ISYMTA)
     &        .OR. LENDE3 .LT. NT2AM(ISYMTB) 
     &        .OR. LENDE3 .LT. NT2AM(ISYMTC) ) THEN
            CALL QUIT('Insufficient work space in CC_BMAT.(3)')
          END IF

* read double excitation part of T^C vector and unpack to full matrix:
          IOPT = 2
          CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KENDE3))
          CALL CCLR_DIASCL(WORK(KENDE3),TWO,ISYMTC)
          CALL CC_T2SQ(WORK(KENDE3),WORK(KT2AMPC),ISYMTC)

* calculate the contribution to THETA2:
          CALL CCRHS_E(WORK(KTHETA2),WORK(KT2AMPC),WORK(KFABVV),
     &                 WORK(KFABOO),WORK(KENDE3),LENDE3,ISYMTC,ISYFAB)


* read double excitation part of T^A vector and unpack to full matrix:
          IOPT = 2
          CALL CC_RDRSP(LISTA,IDLSTA,ISYMTA,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KENDE3))
          CALL CCLR_DIASCL(WORK(KENDE3),TWO,ISYMTA)
          CALL CC_T2SQ(WORK(KENDE3),WORK(KT2AMPA),ISYMTA)

* calculate the contribution to THETA2:
          CALL CCRHS_E(WORK(KTHETA2),WORK(KT2AMPA),WORK(KFBCVV),
     &                 WORK(KFBCOO),WORK(KENDE3),LENDE3,ISYMTA,ISYFBC)


* read double excitation part of T^B vector and unpack to full matrix:
          IOPT = 2
          CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KENDE3))
          CALL CCLR_DIASCL(WORK(KENDE3),TWO,ISYMTB)
          CALL CC_T2SQ(WORK(KENDE3),WORK(KT2AMPB),ISYMTB)

* calculate the contribution to THETA2:
          CALL CCRHS_E(WORK(KTHETA2),WORK(KT2AMPB),WORK(KFCAVV),
     &                 WORK(KFCAOO),WORK(KENDE3),LENDE3,ISYMTB,ISYFCA)


* restore L(ia|jb) integrals: 
          IF (ITRAN.LT.NCTRAN) THEN
            Call CCG_RDIAJB(WORK(KLIAJB),NT2AM(ISYOVOV))
            IOPTTCME = 1
            Call CCSD_TCMEPK(WORK(KLIAJB),ONE,ISYOVOV,IOPTTCME)
          END IF

        END IF

*---------------------------------------------------------------------*
* save result vector:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'THETA after E term:'
        CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
        WRITE (LUPRI,*) MSGDBG,'IOPTRES:',IOPTRES
        WRITE (LUPRI,*) MSGDBG,'FILCMA:',FILCMA
        WRITE (LUPRI,*) MSGDBG,'LUCMAT:',LUCMAT
        WRITE (LUPRI,*) MSGDBG,'ICTRAN(4,ITRAN):',ICTRAN(4,ITRAN)
        CALL FLSHFO(LUPRI)
      END IF

      KTHETA0 = -999999


      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN

*       write to a common direct acces file,
*       store start address in ICTRAN(4,ITRAN)

        ICTRAN(4,ITRAN) = IADRTH

        CALL PUTWA2(LUCMAT,FILCMA,WORK(KTHETA1),IADRTH,NT1AM(ISYRES))
        IADRTH = IADRTH + NT1AM(ISYRES)

        IF (.NOT.CCS) THEN
          CALL PUTWA2(LUCMAT,FILCMA,WORK(KTHETA2),IADRTH,NT2AM(ISYRES))
          IADRTH = IADRTH + NT2AM(ISYRES)
        END IF


      ELSE IF (IOPTRES .EQ. 3) THEN

*       write to a sequential file by call to CC_WRRSP,
*       use FILCMA as LIST type and ICTRAN(4,ITRAN) as index

        CALL CC_WRRSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
     &                WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                WORK(KENDE3),LENDE3)

      ELSE IF (IOPTRES .EQ. 4) THEN

*       add to a vector on a sequential file by call to CC_WARSP,
*       use FILCMA as LIST type and ICTRAN(4,ITRAN) as index

        CALL CC_WARSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
     &                WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                WORK(KENDE3),LENDE3)

      ELSE IF (IOPTRES .EQ. 5) THEN
       
*       calculate required dot products

        CALL CCDOTRSP(ICDOTS,CCONS,IOPTW,FILCMA,ITRAN,NCTRAN,MXVEC,
     &                WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &                WORK(KENDE3),LENDE3)

      ELSE
        CALL QUIT('illegal value for IOPT in CC_CMAT.')
      END IF

* end of loop over C matrix transformations:
      END DO

*=====================================================================*
* end of J, E1 and E2 terms. 
*=====================================================================*

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'J and E (CC2/CCSD) term finished...'
        CALL FLSHFO(LUPRI)
      END IF

* that's it for CCS:
      IF (CCS.OR.CCSTST) GOTO 8888

*=====================================================================*
* F term: requires AO integrals...
*
*         Loop over AO integral shells 
*            Loop over C matrix transformations 
*               Loop over AO integral distributions
*  
*                  Calculation of F term contributions
*
*               End loop
*            End loop
*         End loop
*
*  F term drops out for CCS.
*
*=====================================================================*
      IF (.NOT. (CCS.OR.CCSTST) ) THEN


*---------------------------------------------------------------------*
* initialize integral calculation
*---------------------------------------------------------------------*

        KEND = KEND0
        LEND = LEND0

        IF (DIRECT) THEN
           NTOSYM = 1

           IF (HERDIR) THEN
             CALL HERDI1(WORK(KEND),LEND,IPRERI)
           ELSE
             KCCFB1 = KEND
             KINDXB = KCCFB1 + MXPRIM*MXCONT
             KEND   = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
             LEND   = LWORK  - KEND
             CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                   KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                   KFREE,LFREE,KEND,WORK(KCCFB1),WORK(KINDXB),
     *                   WORK(KEND),LEND,IPRERI)
 
             KEND = KFREE
             LEND = LFREE
           END IF

           KENDSV = KEND
           LWRKSV = LEND
        ELSE
           NTOSYM = NSYM
        END IF

*---------------------------------------------------------------------*
* start loop over AO integrals shells:
*---------------------------------------------------------------------*
      DO ISYMD1 = 1, NTOSYM

        IF (DIRECT) THEN
          IF (HERDIR) THEN
             NTOT = MAXSHL
          ELSE
             NTOT = MXCALL
          ENDIF                   
        ELSE
          NTOT = NBAS(ISYMD1)
        END IF

        DO ILLL = 1, NTOT

          IF (DIRECT) THEN
            KEND = KENDSV
            LEND = LWRKSV

            IF (HERDIR) THEN
               CALL HERDI2(WORK(KEND),LEND,INDEXA,ILLL,NUMDIS, 
     &                     IPRINT) 
            ELSE
              CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                    WORK(KODCL1),WORK(KODCL2),
     *                    WORK(KODBC1),WORK(KODBC2),
     *                    WORK(KRDBC1),WORK(KRDBC2),
     *                    WORK(KODPP1),WORK(KODPP2),
     *                    WORK(KRDPP1),WORK(KRDPP2),
     *                    WORK(KCCFB1),WORK(KINDXB),
     *                    WORK(KEND), LEND,IPRERI)
            END IF

            KRECNR = KEND
            KEND   = KRECNR + (NBUFX(0) - 1)/IRAT + 1
            LEND   = LWORK - KEND
 
            IF (LEND .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_CMAT. (1a)')
            END IF

          ELSE
            NUMDIS = 1
          END IF


*---------------------------------------------------------------------*
*         start loop over C matrix transformations:
*---------------------------------------------------------------------*
          DO ITRAN = 1, NCTRAN

            IDLSTA = ICTRAN(1,ITRAN)
            IDLSTB = ICTRAN(2,ITRAN)
            IDLSTC = ICTRAN(3,ITRAN)

            ISYMTA = ILSTSYM(LISTA,IDLSTA)
            ISYMTB = ILSTSYM(LISTB,IDLSTB)
            ISYMTC = ILSTSYM(LISTC,IDLSTC)

            ISYRES = MULD2H(MULD2H(ISYMTA,ISYMTB),ISYMTC)

* single excitation parts of coupled cluster vectors:
            KT1AMP0 = KEND
            KT1AMPA = KT1AMP0 + NT1AM(ISYM0)
            KT1AMPB = KT1AMPA + NT1AM(ISYMTA)
            KT1AMPC = KT1AMPB + NT1AM(ISYMTB)
            KENDF1  = KT1AMPC + NT1AM(ISYMTC)

* Lambda-hole matrices:
            KLAMDH0 = KENDF1
            KLAMDHA = KLAMDH0 + NGLMDT(ISYM0)
            KLAMDHB = KLAMDHA + NGLMDT(ISYMTA)
            KLAMDHC = KLAMDHB + NGLMDT(ISYMTB)
            KENDF1  = KLAMDHC + NGLMDT(ISYMTC)

* Lambda-particle matrices:
            KLAMDP0 = KENDF1
            KLAMDPA = KLAMDP0 + NGLMDT(ISYM0)
            KLAMDPB = KLAMDPA + NGLMDT(ISYMTA)
            KLAMDPC = KLAMDPB + NGLMDT(ISYMTB)
            KENDF1  = KLAMDPC + NGLMDT(ISYMTC)

* the result vector:
            KTHETA1 = KENDF1
            KTHETA2 = KTHETA1 + NT1AM(ISYRES)
            KENDF1  = KTHETA2 + NT2AM(ISYRES)
            LENDF1  = LWORK - KENDF1

            IF (LENDF1 .LT. 0) THEN
              CALL QUIT('Insufficient memory in CC_CMAT.(1F)')
            END IF

  
* read coupled cluster vectors:
            IOPT = 1
            CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,
     &                    WORK(KT1AMP0),WORK(KDUM))

            IOPT = 1
            CALL CC_RDRSP(LISTA,IDLSTA,ISYMTA,IOPT,MODEL,
     &                    WORK(KT1AMPA),WORK(KDUM))

            IOPT = 1
            CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     &                    WORK(KT1AMPB),WORK(KDUM))

            IOPT = 1
            CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                    WORK(KT1AMPC),WORK(KDUM))

* calculate unperturbed Lambda matrices:
            CALL LAMMAT(WORK(KLAMDP0),WORK(KLAMDH0),WORK(KT1AMP0),
     &                  WORK(KENDF1),LENDF1)

* calculate response Lambda matrices:
            CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPA),WORK(KLAMDH0),
     &                       WORK(KLAMDHA),WORK(KT1AMPA),ISYMTA)

            CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPB),WORK(KLAMDH0),
     &                       WORK(KLAMDHB),WORK(KT1AMPB),ISYMTB)

            CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPC),WORK(KLAMDH0),
     &                       WORK(KLAMDHC),WORK(KT1AMPC),ISYMTC)


* read result vector:
            IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
              CALL GETWA2(LUCMAT,FILCMA,WORK(KTHETA1),
     &                    ICTRAN(4,ITRAN),NT1AM(ISYRES)+NT2AM(ISYRES))
            ELSE IF (IOPTRES.EQ.3) THEN
              IOPT = 3
              CALL CC_RDRSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPT,MODEL,
     &                      WORK(KTHETA1),WORK(KTHETA2))
            ELSE IF (IOPTRES.EQ.4) THEN
              CALL DZERO( WORK(KTHETA1), NT1AM(ISYRES) )
              CALL DZERO( WORK(KTHETA2), NT2AM(ISYRES) )
            ELSE IF (IOPTRES.EQ.5) THEN
              CALL DZERO( WORK(KTHETA1), NT1AM(ISYRES) )
              CALL DZERO( WORK(KTHETA2), NT2AM(ISYRES) )
            ELSE
              CALL QUIT('Error in CC_CMAT.')
            END IF

*---------------------------------------------------------------------*
*        loop over number of distributions on the disk:
*---------------------------------------------------------------------*
          DO IDEL2  = 1, NUMDIS

            IF (DIRECT) THEN
              IDEL   = INDEXA(IDEL2)
              IF (NOAUXB) THEN
                IDUM = 1
                CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
              END IF
              ISYDEL = ISAO(IDEL)
            ELSE
              IDEL   = IBAS(ISYMD1) + ILLL
              ISYDEL = ISYMD1
            END IF

*           read AO integral distribution and calculate integrals with
*           one index transformed to occupied MO (particle):

            KXINT  = KENDF1
            KENDF2 = KXINT  + NDISAO(ISYDEL)
            LENDF2 = LWORK - KENDF2

            IF (LENDF2 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_CMAT. (3)')
            END IF

            CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KENDF2),LENDF2,
     &                  WORK(KRECNR),DIRECT)
*.....................................................................*

* set option for CC_MOFCON routine:
            IOPT = 3

*.....................................................................*
            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPA),WORK(KLAMDHA),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPC),WORK(KLAMDH0),
     &                      ISYMTA,ISYMTB,ISYMTC,ISYM0,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPA),WORK(KLAMDHA),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDP0),WORK(KLAMDHC),
     &                      ISYMTA,ISYMTB,ISYM0,ISYMTC,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

*.....................................................................*
            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPC),WORK(KLAMDHC),
     &                      WORK(KLAMDP0),WORK(KLAMDHA),
     &                      ISYMTB,ISYMTC,ISYM0,ISYMTA,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDP0),WORK(KLAMDH0),
     &                      WORK(KLAMDPC),WORK(KLAMDHA),
     &                      ISYMTB,ISYM0,ISYMTC,ISYMTA,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

*.....................................................................*
            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPC),WORK(KLAMDHC),
     &                      WORK(KLAMDPA),WORK(KLAMDH0),
     &                      ISYMTB,ISYMTC,ISYMTA,ISYM0,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)
 
            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDP0),WORK(KLAMDH0),
     &                      WORK(KLAMDPA),WORK(KLAMDHC),
     &                      ISYMTB,ISYM0,ISYMTA,ISYMTC,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)
 
*.....................................................................*


          END DO ! IDEL2
*---------------------------------------------------------------------*
*         end of the loop over integral distributions:
*---------------------------------------------------------------------*
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,'THETA after F term:'
            WRITE (LUPRI,*) MSGDBG, 'ITRAN :',ITRAN
            WRITE (LUPRI,*) MSGDBG, 'LISTA/IDLSTA:',LISTA(1:2),
     &           ICTRAN(1,ITRAN)
            WRITE (LUPRI,*) MSGDBG, 'LISTB/IDLSTB:',LISTB(1:2),
     &           ICTRAN(2,ITRAN)
            WRITE (LUPRI,*) MSGDBG, 'LISTC/IDLSTC:',LISTC(1:2),
     &           ICTRAN(3,ITRAN)
            CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
            WRITE (LUPRI,*) MSGDBG,'IOPTRES:',IOPTRES
            WRITE (LUPRI,*) MSGDBG,'FILCMA:',FILCMA
            WRITE (LUPRI,*) MSGDBG,'LUCMAT:',LUCMAT
            WRITE (LUPRI,*) MSGDBG,'ICTRAN(4,ITRAN):',ICTRAN(4,ITRAN)
            CALL FLSHFO(LUPRI)
          END IF

          KTHETA0 = -999999

          IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
            CALL PUTWA2(LUCMAT,FILCMA,WORK(KTHETA1),
     &                  ICTRAN(4,ITRAN),NT1AM(ISYRES)+NT2AM(ISYRES))
c         ELSE IF (IOPTRES.EQ.3) THEN
c           CALL CC_WRRSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
c    &                    WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
c    &                    WORK(KENDF1),LENDF1)
          ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
            CALL CC_WARSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
     &                    WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                    WORK(KENDF1),LENDF1)
          ELSE IF (IOPTRES.EQ.5) THEN
            CALL CCDOTRSP(ICDOTS,CCONS,IOPTW,FILCMA,ITRAN,NCTRAN,MXVEC,
     &                    WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &                    WORK(KENDF1),LENDF1)
          ELSE
            CALL QUIT('Error in CC_CMAT.')
          END IF

*---------------------------------------------------------------------*
*         end of the loop over C matrix transformations:
*---------------------------------------------------------------------*
          END DO ! ITRAN
       END DO ! ILLL
      END DO ! ISYMD1
*=====================================================================*
* End of Loop over AO-integrals
*=====================================================================*

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG,'F term section finished...'
          CALL FLSHFO(LUPRI)
        END IF

      END IF ! (.NOT.CCS)
*=====================================================================*
* end of F term section
*=====================================================================*

* that's it for CC2 
      IF (CC2) GOTO 8888


*=====================================================================*
* Calculate the N^6 terms A, B, C and D in a loop over all required
* C matrix transformations.
*
* All N^6 terms drop out for CCS and CC2.
*=====================================================================*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSTST) ) THEN


*---------------------------------------------------------------------*
* loop over all C matrix transformations:
*---------------------------------------------------------------------*
      DO ITRAN = 1, NCTRAN
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG, 'N^6 term section:'
          WRITE (LUPRI,*) MSGDBG, 'ITRAN :',ITRAN
          WRITE (LUPRI,*) MSGDBG, 'LISTA/IDLSTA:',LISTA(1:2),
     &         ICTRAN(1,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTB/IDLSTB:',LISTB(1:2),
     &         ICTRAN(2,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTC/IDLSTC:',LISTC(1:2),
     &         ICTRAN(3,ITRAN)
          CALL FLSHFO(LUPRI)
        END IF

        ISYMTA = ILSTSYM(LISTA,ICTRAN(1,ITRAN))
        ISYMTB = ILSTSYM(LISTC,ICTRAN(2,ITRAN))
        ISYMTC = ILSTSYM(LISTB,ICTRAN(3,ITRAN))
        ISYRES = MULD2H(MULD2H(ISYMTA,ISYMTB),MULD2H(ISYMTC,ISYOVOV))

* read result vector:
        KTHETA1 = KEND0
        KTHETA2 = KTHETA1 + NT1AM(ISYRES)
        KEND1   = KTHETA2 + NT2AM(ISYRES)
        LEND1   = LWORK - KEND1

        IF (LEND1.LT.0) THEN
          CALL QUIT('Insufficient memory in CC_CMAT. (1b)')
        END IF

        IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
          CALL GETWA2(LUCMAT,FILCMA,WORK(KTHETA1),
     &                ICTRAN(4,ITRAN),NT1AM(ISYRES)+NT2AM(ISYRES))
        ELSE IF (IOPTRES.EQ.3) THEN
          IOPT = 3
          CALL CC_RDRSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPT,MODEL,
     &                  WORK(KTHETA1),WORK(KTHETA2))
        ELSE IF (IOPTRES.EQ.4) THEN
          CALL DZERO( WORK(KTHETA1), NT1AM(ISYRES) )
          CALL DZERO( WORK(KTHETA2), NT2AM(ISYRES) )
        ELSE IF (IOPTRES.EQ.5) THEN
          CALL DZERO( WORK(KTHETA1), NT1AM(ISYRES) )
          CALL DZERO( WORK(KTHETA2), NT2AM(ISYRES) )
        ELSE
          CALL QUIT('Error in CC_CMAT.')
        END IF

* loop over cyclic permutations of A, B and C:
        DO ICYCLE = 1, 3
          IF (ICYCLE.EQ.1) THEN
             IDLSTA = ICTRAN(1,ITRAN)
             IDLSTB = ICTRAN(2,ITRAN)
             IDLSTC = ICTRAN(3,ITRAN)
           
             LISTXA = LISTA
             LISTXB = LISTB
             LISTXC = LISTC
          ELSE IF (ICYCLE.EQ.2) THEN
             IDLSTB = ICTRAN(1,ITRAN)
             IDLSTC = ICTRAN(2,ITRAN)
             IDLSTA = ICTRAN(3,ITRAN)
             
             LISTXB = LISTA
             LISTXC = LISTB
             LISTXA = LISTC
          ELSE IF (ICYCLE.EQ.3) THEN
             IDLSTC = ICTRAN(1,ITRAN)
             IDLSTA = ICTRAN(2,ITRAN)
             IDLSTB = ICTRAN(3,ITRAN)
             
             LISTXC = LISTA
             LISTXA = LISTB
             LISTXB = LISTC
          ELSE
             CALL QUIT('Error in CC_CMAT.')
          END IF

          ISYMTA = ILSTSYM(LISTXA,IDLSTA)
          ISYMTB = ILSTSYM(LISTXB,IDLSTB)
          ISYMTC = ILSTSYM(LISTXC,IDLSTC)
          ISYMAB = MULD2H(ISYMTA,ISYMTB)


*---------------------------------------------------------------------*
* read single excitation parts of T^A and T^B and keep them in 
* memory during the loop:
*---------------------------------------------------------------------*
        KT1AMPA  = KEND1
        KT1AMPB  = KT1AMPA + NT1AM(ISYMTA)
        KEND2    = KT1AMPB + NT1AM(ISYMTB)
        LEND2    = LWORK - KEND2

        IF (LEND2 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_CMAT. (2b)')
        END IF
       
* read single excitation part of T^A vector:
        IOPT = 1
        CALL CC_RDRSP(LISTXA,IDLSTA,ISYMTA,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KDUM))

* read single excitation part of T^B vector:
        IOPT = 1
        CALL CC_RDRSP(LISTXB,IDLSTB,ISYMTB,IOPT,MODEL,
     &                WORK(KT1AMPB),WORK(KDUM))


*=====================================================================*
* A term:
*=====================================================================*

* calculate Gamma^AB: intermediate:
        ISYX4O = MULD2H(ISYOVOV,ISYMAB)

        KX4O    = KEND2
        KXIAJB  = KX4O    + NGAMMA(ISYX4O)
        KENDA3  = KXIAJB  + NT2AM(ISYOVOV)
        LENDA3  = LWORK   - KENDA3

        IF (LENDA3 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_CMAT. (A)')
        END IF

* read (ia|jb) integrals from file:
        Call CCG_RDIAJB (WORK(KXIAJB),NT2AM(ISYOVOV))

* calculate double one-index transformed (ik|jl) integrals:
        IOPT = 1
        CALL CCG_4O(WORK(KX4O),ISYX4O,WORK(KXIAJB),ISYOVOV,
     &              WORK(KT1AMPA),ISYMTA,WORK(KT1AMPB),ISYMTB,
     &              WORK(KENDA3),LENDA3,IOPT)

* read double excitation part of T^C vector and unpack to full matrix:
        KT2AMPC = KX4O    + NGAMMA(ISYX4O)
        KENDA3  = KT2AMPC + NT2SQ(ISYMTC)
        LENDA3  = LWORK   - KENDA3

        IF (LENDA3 .LE. NT2AM(ISYMTC)) THEN
          CALL QUIT('Insufficient work space in CC_CMAT. (A2)')
        END IF

        IOPT = 2
        CALL CC_RDRSP(LISTXC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                WORK(KDUM),WORK(KENDA3))

        CALL CCLR_DIASCL(WORK(KENDA3),TWO,ISYMTC)

        CALL CC_T2SQ(WORK(KENDA3),WORK(KT2AMPC),ISYMTC)


* contract T^C with Gamma^AB to A term contribution:
        IOPT = 1
        CALL CCRHS_A(WORK(KTHETA2),WORK(KT2AMPC),WORK(KX4O),
     &               WORK(KENDA3),LENDA3,ISYX4O,ISYMTC,IOPT)

        IF (LOCDBG) THEN
          XNORM=DDOT(NGAMMA(ISYX4O),WORK(KX4O),1,WORK(KX4O),1)
          WRITE (LUPRI,*) 'Norm of X4O intermediate:',XNORM
          WRITE (LUPRI,*) MSGDBG,'THETA after A term:'
          WRITE (LUPRI,*) MSGDBG, 'ITRAN/ICYCLE :',ITRAN,ICYCLE
          WRITE (LUPRI,*) MSGDBG, 'LISTA/IDLSTA:',LISTA(1:2),
     &         ICTRAN(1,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTB/IDLSTB:',LISTB(1:2),
     &         ICTRAN(2,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTC/IDLSTC:',LISTC(1:2),
     &         ICTRAN(3,ITRAN)
          XNORM=DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) MSGDBG, 'Norm of THETA2:',XNORM
          CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
          CALL FLSHFO(LUPRI)
        END IF


*=====================================================================*
* B term:
*=====================================================================*
        ISYMKC  = MULD2H(ISYOVOV,ISYMTC)

        KKINTC  = KEND2
        KXIAJB  = KKINTC  + N3ORHF(ISYMKC)
        KT2AMPC = KXIAJB  + NT2SQ(ISYOVOV)
        KENDB3  = KT2AMPC + MAX(NT2AM(ISYMTC),NT2AM(ISYOVOV))
        LENDB3  = LWORK - KENDB3

        IF (LENDB3 .LT. 0) THEN
          CALL QUIT('Insufficient memory in CC_CMAT. (B)')
        END IF

* read (ia|jb) integrals from file and unpack them to a full matrix:
        Call CCG_RDIAJB (WORK(KT2AMPC),NT2AM(ISYOVOV))

        CALL CC_T2SQ(WORK(KT2AMPC),WORK(KXIAJB),ISYOVOV)


* read (ia|jb) integrals from file:
        IOPT = 2
        CALL CC_RDRSP(LISTXC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                WORK(KDUM),WORK(KT2AMPC))
        CALL CCLR_DIASCL(WORK(KT2AMPC),TWO,ISYMTC)

* construct PHI^C intermediate:
        CALL CC_MI(WORK(KKINTC),WORK(KXIAJB),ISYOVOV,
     &             WORK(KT2AMPC),ISYMTC,WORK(KENDB3),LENDB3)

* for CCSD(R12) add correction term:
        IF (CCR12.AND..NOT.(CCS.OR.CC2)) THEN
          KTR12AM = KXIAJB
          KVINT   = KTR12AM + NTR12AM(ISYMTC)
          KENDB3  = KVINT   + NTR12AM(1)
          LENDB3  = LWORK - KENDB3
          IF (LENDB3.LT.0) THEN
            CALL QUIT('Insufficient work space in CC_CMAT (B-R12)')
          END IF

          IOPT = 32
          CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,WORK(KDUM),
     &                  WORK(KTR12AM))

          LUNIT = -1
          CALL GPOPEN(LUNIT,FCCR12V,'OLD',' ','UNFORMATTED',
     &                KDUM,.FALSE.)
 9999     READ(LUNIT) IAN
          READ(LUNIT) (WORK(KVINT-1+I), I=1, NTR12AM(1))
          IF (IAN.NE.IANR12) GOTO 9999
          CALL GPCLOSE(LUNIT,'KEEP')

          IOPT = 2
          CALL CC_R12CV(WORK(KKINTC),WORK(KTR12AM),ISYMTC,WORK(KVINT),
     &                  1,IOPT,WORK(KENDB3),LENDB3)

        END IF

* double oneindex-transform PHI^C intermediate to result vector:
        CALL CCC_OVOV(WORK(KTHETA2),ISYRES,WORK(KKINTC),ISYMKC,
     &                WORK(KT1AMPA),ISYMTA,WORK(KT1AMPB),ISYMTB,
     &                WORK(KENDB3),LENDB3)



        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG,'THETA after B term:'
          WRITE (LUPRI,*) MSGDBG, 'ITRAN/ICYCLE :',ITRAN,ICYCLE
          WRITE (LUPRI,*) MSGDBG, 'LISTA/IDLSTA:',LISTA(1:2),
     &         ICTRAN(1,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTB/IDLSTB:',LISTB(1:2),
     &         ICTRAN(2,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTC/IDLSTC:',LISTC(1:2),
     &         ICTRAN(3,ITRAN)
          XNORM=DDOT(N3ORHF(ISYMKC),WORK(KKINTC),1,WORK(KKINTC),1)
          WRITE (LUPRI,*) MSGDBG, 'Norm of KINTC:',XNORM
          XNORM=DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) MSGDBG, 'Norm of THETA2:',XNORM
          CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
          CALL FLSHFO(LUPRI)
        END IF

*=====================================================================*
* C term: requires 5 x 1/2 V^2 O^2  !!!
* D term: requires 5 x 1/2 V^2 O^2  !!!
*=====================================================================*
        ISYCBAR = MULD2H(ISYMTC,ISYOVOV)
        ISYDBAR = MULD2H(ISYMTC,ISYOVOV)

        KXIAJB  = KEND2
        KT2AMPC = KXIAJB  + NT2SQ(ISYOVOV)
        KCBAR   = KT2AMPC + MAX(NT2AM(ISYMTC),NT2AM(ISYOVOV))
        KEND3   = KCBAR   + NT2SQ(ISYCBAR)
        LEND3   = LWORK - KEND3

        KDBAR   = KCBAR

        IF (LEND3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (C)')
        END IF

* read (ia|jb) and square them:
        Call CCG_RDIAJB (WORK(KT2AMPC),NT2AM(ISYOVOV))

        Call CC_T2SQ (WORK(KT2AMPC), WORK(KXIAJB), ISYOVOV)

* read double excitation part of T^C vector:
        IOPT = 2
        CALL CC_RDRSP(LISTXC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                WORK(KDUM),WORK(KT2AMPC)          )

        CALL CCLR_DIASCL(WORK(KT2AMPC),TWO,ISYMTC)

* calculate CBAR^C intermediate:
        IOPT   = 1
        CBAFIL = '--------'
        LUCBAR = -1
        IOFFCD = -1
        CALL CCB_CDBAR('C',WORK(KXIAJB),ISYOVOV, WORK(KT2AMPC),ISYMTC,
     &                     WORK(KCBAR), ISYCBAR, WORK(KEND3),LEND3,
     &                     CBAFIL, LUCBAR, IOFFCD, IOPT)

* double transform CBAR^C intermediate to C term of the result vector:
        CALL CCB_22CD(WORK(KTHETA2),ISYRES,WORK(KCBAR),ISYCBAR, 
     &                WORK(KT1AMPA),ISYMTA,WORK(KT1AMPB),ISYMTB,
     &                'C', WORK(KEND3),LEND3)


* calculate DBAR^D intermediate:
        IOPT   = 1
        DBAFIL = '--------'
        LUDBAR = -1
        IOFFCD = -1
        CALL CCB_CDBAR('D',WORK(KXIAJB),ISYOVOV, WORK(KT2AMPC),ISYMTC,
     &                     WORK(KDBAR), ISYDBAR, WORK(KEND3),LEND3,
     &                     DBAFIL, LUDBAR, IOFFCD, IOPT)

* double transform DBAR^C intermediate to D term of the result vector:
        CALL CCB_22CD(WORK(KTHETA2),ISYRES,WORK(KDBAR),ISYDBAR, 
     &                WORK(KT1AMPA),ISYMTA,WORK(KT1AMPB),ISYMTB,
     &                'D', WORK(KEND3),LEND3)


        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG,'THETA after C and D terms:'
          WRITE (LUPRI,*) MSGDBG, 'ITRAN/ICYCLE :',ITRAN,ICYCLE
          WRITE (LUPRI,*) MSGDBG, 'LISTA/IDLSTA:',LISTA(1:2),
     &         ICTRAN(1,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTB/IDLSTB:',LISTB(1:2),
     &         ICTRAN(2,ITRAN)
          WRITE (LUPRI,*) MSGDBG, 'LISTC/IDLSTC:',LISTC(1:2),
     &         ICTRAN(3,ITRAN)
          XNORM=DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) MSGDBG, 'Norm of THETA2:',XNORM
          CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
          CALL FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
* End loop over cyclic permutations, save result vector
*---------------------------------------------------------------------*
        END DO ! ICYCLE

        KTHETA0 = -999999

        IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
          CALL PUTWA2(LUCMAT,FILCMA,WORK(KTHETA1),
     &                ICTRAN(4,ITRAN),NT1AM(ISYRES)+NT2AM(ISYRES))
c       ELSE IF (IOPTRES.EQ.3) THEN
c         CALL CC_WRRSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
c    &                  WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
c    &                  WORK(KEND1),LEND1)
        ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
          CALL CC_WARSP(FILCMA,ICTRAN(4,ITRAN),ISYRES,IOPTW,MODELW,
     &                  WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                  WORK(KEND1),LEND1)
        ELSE IF (IOPTRES.EQ.5) THEN
          CALL CCDOTRSP(ICDOTS,CCONS,IOPTW,FILCMA,ITRAN,NCTRAN,MXVEC,
     &                  WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &                  WORK(KEND1),LEND1)
        ELSE
          CALL QUIT('Error in CC_CMAT.')
        END IF
*---------------------------------------------------------------------*
* End loop over C matrix transformations
*---------------------------------------------------------------------*
      END DO

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG,'N^6 term section finished...'
          CALL FLSHFO(LUPRI)
        END IF

      END IF ! (.NOT. (CCS .OR. CC2))
*=====================================================================*
* End of the N^6 term section
*=====================================================================*

*=====================================================================*
* restore result vectors and clean up and return:
*=====================================================================*
8888   CONTINUE

*---------------------------------------------------------------------*
* if IOPTRES=1 and enough work space available, read result
* vectors back into memory:
*---------------------------------------------------------------------*

* check size of work space:
      IF (IOPTRES .EQ. 1) THEN
        LENALL = IADRTH-1
        IF (LENALL .GT. LWORK) IOPTRES = 0
      END IF

* read the result vectors back into memory:
      IF (IOPTRES .EQ. 1) THEN

        CALL GETWA2(LUCMAT,FILCMA,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
          DO ITRAN = 1, NCTRAN
            IF (ITRAN.LT.NCTRAN) THEN
              LEN     = ICTRAN(4,ITRAN+1)-ICTRAN(4,ITRAN)
            ELSE
              LEN     = IADRTH-ICTRAN(4,NCTRAN)
            END IF
            KTHETA1 = ICTRAN(4,ITRAN)
            XNORM   = DDOT(LEN, WORK(KTHETA1),1, WORK(KTHETA1),1)
            WRITE (LUPRI,*) 'Read C matrix transformation nb. ',ITRAN
            WRITE (LUPRI,*) 'Adress, length, NORM:',ICTRAN(4,NCTRAN),
     &           LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

*---------------------------------------------------------------------*
* close C matrix file & return
*---------------------------------------------------------------------*
* check return option for the result vectors:
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
        CALL WCLOSE2(LUCMAT,FILCMA,'KEEP')

      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_CMAT.')
      END IF


*=====================================================================*

      CALL QEXIT('CC_CMAT')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CC_CMAT
*=====================================================================*

*---------------------------------------------------------------------*
c/* Deck CC_CTST */
*=====================================================================*
       SUBROUTINE CC_CTST(WORK,LWORK)
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_CTST> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER MXCTRAN
      PARAMETER (MXCTRAN = 2)

      INTEGER LWORK
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL DDOT, RDUM
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION DDOT, RDUM
#endif 

      CHARACTER*(3) LISTA, LISTB, LISTC, LISTD
      CHARACTER*(8) FILCMA, LABELA
      CHARACTER*(10) MODEL
      INTEGER IOPTRES, IDUM
      INTEGER ICTRAN(4,MXCTRAN), NCTRAN
      INTEGER IDLSTA, IDLSTB, IDLSTC, ISYMA, ISYMB, ISYMC, ISYMABC
      INTEGER KTHETA1, KTHETA2, KT1AMPC, KT2AMPC, KRESLT1, KRESLT2
      INTEGER KEND1, LEND1, IOPT, ISYMD, KT1AMPD, KT2AMPD, IDLSTD

* external function:
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER ILSTSYM


      CALL QENTER('CC_CTST')


*---------------------------------------------------------------------*
* call C matrix transformation:
*---------------------------------------------------------------------*
      LISTA = 'R1'
      LISTB = 'R1'
      LISTC = 'R1'
      LISTD = 'L1'
      IDLSTA = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTB = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTC = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTD = IL1ZETA('ZDIPLEN ',.FALSE.,0.0D0,1)
      ICTRAN(1,1) = IDLSTA
      ICTRAN(2,1) = IDLSTB
      ICTRAN(3,1) = IDLSTC
      NCTRAN = 1

      IOPTRES = 1
      FILCMA  = 'CCCMAT'

      WRITE (LUPRI,*) 'CCTST: call CC_CMAT:'

      CALL CC_CMAT(ICTRAN,  NCTRAN,  LISTA,  LISTB, LISTC,
     &             IOPTRES, FILCMA, IDUM, RDUM, 0, WORK, LWORK )


      ISYMA  = ILSTSYM(LISTA,IDLSTA)
      ISYMB  = ILSTSYM(LISTB,IDLSTB)
      ISYMC  = ILSTSYM(LISTC,IDLSTC)
      ISYMD  = ILSTSYM(LISTD,IDLSTD)
      ISYMABC = MULD2H(MULD2H(ISYMA,ISYMB),ISYMC)

      KTHETA1 = ICTRAN(4,1)
      KTHETA2 = KTHETA1 + NT1AM(ISYMABC)
      KEND1   = KTHETA2 + NT2AM(ISYMABC)
      LEND1   = LWORK - KEND1

      IF (NSYM.EQ.1 .AND. LOCDBG) THEN
        KT1AMPC = KTHETA2 + NT2AM(ISYMABC)
        KT2AMPC = KT1AMPC + NT1AM(ISYMC)
        KRESLT1 = KT2AMPC + NT2AM(ISYMC)
        KRESLT2 = KRESLT1 + NT1AM(ISYMABC)
        KEND1   = KRESLT2 + NT2AM(ISYMABC)
        LEND1   = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_CTST.')
        END IF

        IOPT = 3
        Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &                WORK(KT1AMPC),WORK(KT2AMPC))

        ! zero singles or doubles C vector:
C       CALL DZERO(WORK(KT1AMPC),NT1AM(ISYMC))
C       CALL DZERO(WORK(KT2AMPC),NT2AM(ISYMC))
        CALL DZERO(WORK(KRESLT1),NT1AM(ISYMABC)+NT2AM(ISYMABC))
        IPRINT  = 5

        CALL CC_FDC(NT1AM(ISYMABC),NT2AM(ISYMABC),
     >              LISTA,IDLSTA,LISTB,IDLSTB,  
     >              WORK(KT1AMPC), WORK(KRESLT1),
     >              WORK(KEND1), LEND1)

        IPRINT  = 0

        IF (.TRUE.) THEN
          WRITE (LUPRI,*) 
     *          'LISTA, IDLSTA, ISYMA:',LISTA(1:2),IDLSTA,ISYMA
          WRITE (LUPRI,*) 
     *          'LISTB, IDLSTB, ISYMB:',LISTB(1:2),IDLSTB,ISYMB
          WRITE (LUPRI,*) 
     *          'LISTC, IDLSTC, ISYMC:',LISTC(1:2),IDLSTC,ISYMC
          WRITE (LUPRI,*) 'ISYMABC:',ISYMABC
          WRITE (LUPRI,*)
          WRITE (LUPRI,*) 'finite difference Theta vector:'
          Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMABC,1,1)
          WRITE (LUPRI,*) 'analytical Theta vector:'
          Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMABC,1,1)
        END IF

        KT1AMPD = KEND1
        KT2AMPD = KT1AMPD + NT1AM(ISYMABC)
        KEND1   = KT2AMPD + NT2AM(ISYMABC)
        LEND1   = LWORK - KEND1
        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_CTST.')
        END IF

        IOPT = 3
        Call CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,
     &                WORK(KT1AMPD),WORK(KT2AMPD))
        CALL CCLR_DIASCL(WORK(KT2AMPD),0.5d0,ISYMD)

        WRITE (LUPRI,*) 'Dot product with T^D for analytical C matrix:',
     >   DDOT(NT1AM(ISYMD),WORK(KT1AMPD),1,WORK(KTHETA1),1)+
     >   DDOT(NT2AM(ISYMD),WORK(KT2AMPD),1,WORK(KTHETA2),1) 

        WRITE (LUPRI,*) 'Dot product with T^D for numerical C matrix:',
     >   DDOT(NT1AM(ISYMD),WORK(KT1AMPD),1,WORK(KRESLT1),1)+
     >   DDOT(NT2AM(ISYMD),WORK(KT2AMPD),1,WORK(KRESLT2),1) 



        Call DAXPY(NT1AM(ISYMABC),-1.0d0,WORK(KTHETA1),1,
     &                                  WORK(KRESLT1),1)
        IF (.NOT.CCS) THEN
          Call DAXPY(NT2AM(ISYMABC),-1.0d0,WORK(KTHETA2),1,
     &                                    WORK(KRESLT2),1)
        ELSE
          Call DZERO(WORK(KRESLT2),NT2AM(ISYMABC))
        END IF

        WRITE (LUPRI,*) 'Norm of difference between analytical THETA '
     >           // 'vector and the numerical result:'
        WRITE (LUPRI,*) 'singles excitation part:',
     >   DSQRT(DDOT(NT1AM(ISYMA),WORK(KRESLT1),1,WORK(KRESLT1),1))
        WRITE (LUPRI,*) 'double excitation part: ',
     >   DSQRT(DDOT(NT2AM(ISYMA),WORK(KRESLT2),1,WORK(KRESLT2),1))

        WRITE (LUPRI,*) 'difference vector:'
        Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMABC,1,1)

        CALL FLSHFO(LUPRI)


      ELSE IF (NSYM.NE.1 .AND. LOCDBG) THEN
        WRITE (LUPRI,*) 'analytical Theta vector:'
        Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMABC,1,1)

        WRITE (LUPRI,*) 'CC_CTST> can not calculate finite '//
     &       'difference C matrix'
        WRITE (LUPRI,*) 'CC_CTST> with symmetry.'

        KT1AMPD = KEND1
        KT2AMPD = KT1AMPD + NT1AM(ISYMABC)
        KEND1   = KT2AMPD + NT2AM(ISYMABC)
        LEND1   = LWORK - KEND1
        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_CTST.')
        END IF

        IOPT = 3
        Call CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,
     &                WORK(KT1AMPD),WORK(KT2AMPD))
        CALL CCLR_DIASCL(WORK(KT2AMPD),0.5d0,ISYMD)

        WRITE (LUPRI,*) 'Dot product with T^D for analytical C matrix:',
     >   DDOT(NT1AM(ISYMD),WORK(KT1AMPD),1,WORK(KTHETA1),1)+
     >   DDOT(NT2AM(ISYMD),WORK(KT2AMPD),1,WORK(KTHETA2),1) 

      END IF

      CALL QEXIT('CC_CTST')

      RETURN
      END 
*=====================================================================*
*---------------------------------------------------------------------*
c/*  Deck CCC_OVOV */
*=====================================================================*
      SUBROUTINE CCC_OVOV(THETA2, ISYRES, XKINT,  ISYXKI, 
     &                    T1AMPA, ISYMTA, T1AMPB, ISYMTB,
     &                    WORK, LWORK )
*---------------------------------------------------------------------*
*
*     Purpose: double one-index transform XKINT_{iljk} intermediate 
*              with T^A and T^B single-excitation amplitudes 
*
*              THETA2(ai,bj) += P^{AB}  t^A_{ak} t^B_{bl} XKINT_{iljk}
*
*              needed for CCSD part of C matrix transformation
*
*     Christof Haettig, Maj 1997
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include <implicit.h>
#endif
 
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
 
#if defined (SYS_CRAY)
      REAL ZERO, ONE
#else
      DOUBLE PRECISION ZERO, ONE, TWO
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMTA, ISYMTB, ISYXKI, ISYRES
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL THETA2(*)     ! dimension (NT2AM(ISYRES))
      REAL XKINT(*)      ! dimension (N3ORHF(ISYXKI))
      REAL T1AMPA(*)     ! dimension (NT1AM(ISYMTA))
      REAL T1AMPB(*)     ! dimension (NT1AM(ISYMTB))
      REAL WORK(LWORK)
      REAL SUM1, DDOT, SUM2
#else
      DOUBLE PRECISION THETA2(*)     ! dimension (NT2AM(ISYRES))
      DOUBLE PRECISION XKINT(*)      ! dimension (N3ORHF(ISYXKI))
      DOUBLE PRECISION T1AMPA(*)     ! dimension (NT1AM(ISYMTA))
      DOUBLE PRECISION T1AMPB(*)     ! dimension (NT1AM(ISYMTB))
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION SUM1, DDOT, SUM2
#endif

      INTEGER ISYMA, ISYMI, ISYMB, ISYMJ, ISYMK, ISYML
      INTEGER ISYMAI, ISYMBJ, ISYMJL, ISYJLI, NAIBJ
      INTEGER KSCR1, KSCR2, KEND1, LEND1, KEND2, LEND2, KOFF1, KOFF2
      INTEGER KOFFX, KOFFT, NTOTJLI, NTOTA, NTOTB, NTOTJ, NAI, NBJ


* statement function:
      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCC_OVOV')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
* check symmetries:
      IF (MULD2H(ISYMTA,ISYMTB) .NE. MULD2H(ISYRES,ISYXKI)) THEN
        CALL QUIT('Symmetry mismatch in CCC_OVOV.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CCC_OVOV> ISYRES:', ISYRES
        WRITE (LUPRI,*) 'CCC_OVOV> ISYXKI:', ISYXKI
        WRITE (LUPRI,*) 'CCC_OVOV> ISYMTA:', ISYMTA
        WRITE (LUPRI,*) 'CCC_OVOV> ISYMTB:', ISYMTB
        WRITE (LUPRI,*) 'CCC_OVOV> XKINT:'
C       WRITE (LUPRI,'(5F14.6)') XKINT
      END IF

      SUM1 = 0.0d0

*---------------------------------------------------------------------*
* loop over A indeces:
*---------------------------------------------------------------------*
      DO ISYMA = 1, NSYM
        ISYMK  = MULD2H(ISYMA,ISYMTA)
        ISYJLI = MULD2H(ISYXKI,ISYMK)

      IF (NRHF(ISYMK).NE.0) THEN


        KSCR1 = 1
        KEND1 = KSCR1 + NMAIJK(ISYJLI)
        LEND1 = LWORK - KEND1
    
        IF (LEND1.LT.0) THEN
          CALL QUIT('Insufficient memory in CCC_OVOV.')
        END IF

      DO A = 1, NVIR(ISYMA)

*---------------------------------------------------------------------*
* transform K index:  sum_k XKINT_{jli;k} * t^A(ka)
*---------------------------------------------------------------------*
          KOFFX = I3ORHF(ISYJLI,ISYMK) + 1
          KOFFT = IT1AM(ISYMA,ISYMK) + A

          NTOTJLI = MAX(1,NMAIJK(ISYJLI))
          NTOTA   = MAX(1,NVIR(ISYMA))

          Call DGEMV('N', NMAIJK(ISYJLI), NRHF(ISYMK),
     &                ONE,  XKINT(KOFFX),  NTOTJLI, 
     &                      T1AMPA(KOFFT), NTOTA,
     &                ZERO, WORK(KSCR1),   1          )

*---------------------------------------------------------------------*
* loop over I indeces:
*---------------------------------------------------------------------*
      DO ISYMI = 1, NSYM
        ISYMJL = MULD2H(ISYJLI,ISYMI)
        ISYMAI = MULD2H(ISYMA,ISYMI)
        ISYMBJ = MULD2H(ISYRES,ISYMAI)

        KSCR2 = KEND1
        KEND2 = KSCR2 + NT1AM(ISYMBJ)
        LEND2 = LWORK - KEND2
   
        IF (LEND2.LT.0) THEN
          CALL QUIT('Insufficient memory in CCC_OVOV.')
        END IF


      DO I = 1, NRHF(ISYMI)

*---------------------------------------------------------------------*
* transform L index:  sum_l t^B(bl) * SCR_{jl;i}
*---------------------------------------------------------------------*
      DO ISYMB = 1, NSYM
        ISYMJ = MULD2H(ISYMBJ,ISYMB)
        ISYML = MULD2H(ISYMTB,ISYMB)
        IF (MULD2H(ISYMJ,ISYML).NE.ISYMJL)
     &       CALL QUIT('Error in CCC_OVOV (1).')
        
        KOFF1 = KSCR1 + IMAIJK(ISYMJL,ISYMI) +
     &            NMATIJ(ISYMJL)*(I-1) + IMATIJ(ISYMJ,ISYML)
        KOFF2 = KSCR2 + IT1AM(ISYMB,ISYMJ)
        KOFFT = IT1AM(ISYMB,ISYML) + 1

        NTOTJ = MAX(1,NRHF(ISYMJ))
        NTOTB = MAX(1,NVIR(ISYMB))

        Call DGEMM('N','T',NVIR(ISYMB),NRHF(ISYMJ),NRHF(ISYML),
     &              ONE, T1AMPB(KOFFT), NTOTB, WORK(KOFF1), NTOTJ,
     &             ZERO, WORK(KOFF2),   NTOTB )
      END DO

*---------------------------------------------------------------------*
* store result in THETA2 vector:
*---------------------------------------------------------------------*
      NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A

      IF      (ISYMAI .EQ. ISYMBJ) THEN

        WORK(KSCR2-1+NAI) = TWO * WORK(KSCR2-1+NAI)

        DO NBJ = 1, NT1AM(ISYMBJ)
          NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
          THETA2(NAIBJ) = THETA2(NAIBJ) + WORK(KSCR2-1+NBJ)
        END DO

      ELSE IF (ISYMAI .LT. ISYMBJ) THEN

        DO NBJ = 1, NT1AM(ISYMBJ)
          NAIBJ = IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)+NAI
          THETA2(NAIBJ) = THETA2(NAIBJ) + WORK(KSCR2-1+NBJ)
        END DO

      ELSE IF (ISYMBJ .LT. ISYMAI) THEN

        DO NBJ = 1, NT1AM(ISYMBJ)
          NAIBJ = IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMBJ)*(NAI-1)+NBJ
          THETA2(NAIBJ) = THETA2(NAIBJ) + WORK(KSCR2-1+NBJ)
        END DO

      ELSE
        CALL QUIT('Error in CCC_OVOV (2).')
      END IF

      

*---------------------------------------------------------------------*
* end loops over I and A:
*---------------------------------------------------------------------*
      END DO ! I
      END DO ! ISYMI

      END DO ! A
      END IF ! NRHF(ISYMK).NE.0
      END DO ! ISYMA

      CALL QEXIT('CCC_OVOV')

      RETURN
      END
*---------------------------------------------------------------------*
*                   END OF SUBROUTINE CCC_OVOV                        *
*---------------------------------------------------------------------*
