
!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
*---------------------------------------------------------------------*
c/* Deck CC_KMATRIX */
*=====================================================================*
      SUBROUTINE CC_KMATRIX(IKTRAN, NKTRAN, LISTA, LISTB, IOPTRES,
     &                    FILKMA, IKDOTS, KCONS, MXVEC, WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over K matrix transformations
*             (needed if the number of transformations exceeds the
*              limit MAXSIM defined on ccsdio.h )
*        
*    CCMM JK+OC, modyfied CC_BMATRIX
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "maxorb.h"
#include "ccsdio.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*(*) LISTA, LISTB, FILKMA
      INTEGER IOPTRES
      INTEGER NKTRAN, MXVEC, LWORK
      INTEGER IKTRAN(3,NKTRAN)
      INTEGER IKDOTS(MXVEC,NKTRAN)
      
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL KCONS(MXVEC,NKTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION KCONS(MXVEC,NKTRAN)
#endif 

      INTEGER MAXKTRAN, NTRAN, ISTART, IBATCH, NBATCH

      MAXKTRAN = MAXSIM

      NBATCH = (NKTRAN+MAXKTRAN-1)/MAXKTRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Batching over K matrix transformations:'
        WRITE (LUPRI,*) 'nb. of batches needed:', NBATCH
      END IF
  
      DO IBATCH = 1, NBATCH
        ISTART = (IBATCH-1) * MAXKTRAN + 1
        NTRAN  = MIN(NKTRAN-(ISTART-1),MAXKTRAN)

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',ISTART
          WRITE (LUPRI,*) '# transf.:',NTRAN
        END IF

        CALL CC_KMAT(IKTRAN(1,ISTART), NTRAN,
     &                LISTA, LISTB, IOPTRES, FILKMA, 
     &                IKDOTS(1,ISTART), KCONS(1,ISTART), 
     &                MXVEC, WORK,LWORK)

      END DO

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_KMATRIX                           *
*---------------------------------------------------------------------*

*---------------------------------------------------------------------*
c/* Deck CC_BMAT */
*=====================================================================*
      SUBROUTINE CC_KMAT(IBTRAN, NBTRAN, LISTA, LISTB, IOPTRES,
     &                    FILBMA, IBDOTS, BCONS, MXVEC, WORK, LWORK )
*---------------------------------------------------------------------*
*
*             The linear transformations are calculated for a list
*             of bar{T^A} vectors and a list of bar{T^B} vectors: 
*
*                LISTA       -- type of bar{T^A} vectors
*                LISTB       -- type of bar{T^B} vectors
*                IBTRAN(1,*) -- indeces of bar{T^A} vectors
*                IBTRAN(2,*) -- indeces of bar{T^B} vectors
*                IBTRAN(3,*) -- indeces or addresses of result vectors
*                NBTRAN      -- number of requested transformations
*                FILBMA      -- file name / list type of result vectors
*                               or list type of vectors to be dotted on
*                IBDOTS      -- indeces of vectors to be dotted on
*                BCONS       -- contains the dot products on return
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  all result vectors are written to a direct
*                          access file, FILBMA is used as file name
*                          the start addresses of the vectors are
*                          returned in IBTRAN(3,*)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          IBTRAN(3,*). N.B.: if WORK is not large
*                          enough IOPTRES is automatically reset to 0!!
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, FILBMA is used
*                          as list type and IBTRAN(3,*) as list index
*                          NOTE that IBTRAN(3,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WARSP, FIBBMA is used
*                          as list type and IBTRAN(3,*) as list index
*                          NOTE that IBTRAN(3,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by FILBMA and the indeces from IBDOTS
*                          the result of the dot products is returned
*                          in the BCONS array
*
*
*           CCMM JK+OC, modyfied version of CC_BMAT
*
*=====================================================================*
      USE PELIB_INTERFACE, ONLY: USE_PELIB, PELIB_IFC_QRTRANSFORMER
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsections.h"
#include "maxorb.h"
#include "mxcent.h"
#include "ccsdio.h"
#include "ccorb.h"
#include "iratdef.h"
#include "eribuf.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "second.h"

* local parameters:
      CHARACTER MSGDBG*(16)
      PARAMETER (MSGDBG='[debug] CC_KMAT> ')

      LOGICAL LOCDBG, LSAME
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space

      
      INTEGER LUBMAT

      CHARACTER*(*) LISTA, LISTB, FILBMA

      INTEGER IOPTRES
      INTEGER NBTRAN, MXVEC, LWORK
      INTEGER IBTRAN(3,NBTRAN)
      INTEGER IBDOTS(MXVEC,NBTRAN)

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL ZERO, ONE, TWO
      REAL DUM, XNORM, DUMMY
      REAL BCONS(MXVEC,NBTRAN)
      REAL FACTSLV
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ZERO, ONE, TWO
      DOUBLE PRECISION DUM, XNORM, DUMMY
      DOUBLE PRECISION BCONS(MXVEC,NBTRAN)
      DOUBLE PRECISION FACTSLV
      DOUBLE PRECISION TAL1, TAL2, RHO1N, RHO2N
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0)

      CHARACTER*(10) MODEL, MODELW, CDUMMY
      CHARACTER*8 LABEL
      CHARACTER RSPTYP*(1)
      INTEGER NBATCH


      INTEGER ITRAN, IDLSTA, IDLSTB, IOPT
      INTEGER ISYMA, ISYMB, ISYMAB
      INTEGER IBATCH,IADRTH
      INTEGER KEND1, LEN, LENALL
      INTEGER KEND2, LWRK2, KEND3, LWRK3, LWRK1
      INTEGER IVEC
      INTEGER KT2AMPA, KTHETA0
      INTEGER KTHETA1, KTHETA2, KT1AMPA, KT1AMPB 
      INTEGER IOPTW, IDUMMY
      INTEGER KTGB, KETA, KETA1, KETA2, NAMPF
      INTEGER KTGA

* external functions:
      INTEGER ILSTSYM

#if defined (SYS_CRAY) 
      REAL DDOT, DTIME, TIMALL, TIMTRN
#else 
      DOUBLE PRECISION DDOT, DTIME, TIMALL, TIMTRN
#endif 
  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_KMAT')
        WRITE (LUPRI,*) 'LISTA : ',LISTA
        WRITE (LUPRI,*) 'LISTB : ',LISTB
        WRITE (LUPRI,*) 'FILBMA: ',FILBMA
        WRITE (LUPRI,*) 'NKTRAN: ',NBTRAN
        WRITE (LUPRI,*) 'IOPTRES:',IOPTRES
        CALL FLSHFO(LUPRI)
      END IF
      
      IF (CCSDT) THEN
        WRITE(LUPRI,'(/1x,a)') 'K matrix transformations not '
     &          //'implemented for triples yet...'
        CALL QUIT('Triples not implemented for K matrix '//
     &            'transformations')
      END IF

      IF ( .not. (CCS .or. CC2 .or. CCSD) ) THEN
        WRITE(LUPRI,'(/1x,a)') 'CC_KMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_KMAT...'
        CALL QUIT('Unknown CC method in CC_KMAT.')
      END IF

      IF (.NOT. DUMPCD) THEN
        WRITE(LUPRI,*) 'DUMPCD = ',DUMPCD
        WRITE(LUPRI,*) 'CC_KMAT requires DUMPCD=.TRUE.'
        CALL QUIT('DUMPCD=.FALSE. , CC_KMAT requires DUMPCD=.TRUE.')
      END IF

      IF (.NOT. RSPIM) THEN
        WRITE(LUPRI,*) 'RSPIM = ',RSPIM
        WRITE(LUPRI,*) 'CC_BMAT requires RSPIM=.TRUE.'
        CALL QUIT('RSPIM=.FALSE. , CC_KMAT requires RSPIM=.TRUE.')
      END IF

      IF (ISYMOP .NE. 1) THEN
        WRITE(LUPRI,*) 'ISYMOP = ',ISYMOP
        WRITE(LUPRI,*) 'CC_KMAT is not implemented for ISYMOP.NE.1'
        CALL QUIT('CC_KMAT is not implemented for ISYMOP.NE.1')
      END IF

      IF (NBTRAN .GT. MAXSIM) THEN
        WRITE(LUPRI,*) 'NBTRAN = ', NBTRAN
        WRITE(LUPRI,*) 'MAXSIM = ', MAXSIM
        WRITE(LUPRI,*) 'number of requested transformation is larger'
        WRITE(LUPRI,*) 'than the maximum number of allowed ',
     &                 'simultaneous transformation.'
        WRITE(LUPRI,*) 'Error in CC_KMAT: NBTRAN is larger than MAXSIM.'
        CALL QUIT('Error in CC_KMAT: NBTRAN is larger than MAXSIM.')
      END IF

      IF (IPRINT.GT.0) THEN
 
         WRITE (LUPRI,'(//1X,A1,50("="),A1)')'+','+'

         WRITE (LUPRI,'(1x,A52)')
     &         '|        K MATRIX TRANSFORMATION SECTION           |'

         IF (IOPTRES.EQ.3) THEN
            WRITE (LUPRI,'(1X,A52)')
     &         '|          (result is written to file)             |'
         ELSE IF (IOPTRES.EQ.4) THEN
            WRITE (LUPRI,'(1X,A52)')
     &         '|     (result is added to a vector on file)        |'
         ELSE IF (IOPTRES.EQ.5) THEN
            WRITE (LUPRI,'(1X,A52)')
     &         '|    (result used to calculate dot products)       |'
         END IF
        
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'

      END IF

* initialize timings:
      TIMALL  = SECOND()

* set option and model to write vectors to file:
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
         CALL QUIT('Unknown coupled cluster model in CC_KMAT.')
      END IF


* check return option for the result vectors:
      LUBMAT = -1
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
         CALL WOPEN2(LUBMAT, FILBMA, 64, 0)
      ELSE IF (IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4) THEN
         CONTINUE
      ELSE IF (IOPTRES .EQ. 5) THEN
         IF (MXVEC*NBTRAN.NE.0) CALL DZERO(BCONS,MXVEC*NBTRAN)
      ELSE
         CALL QUIT('Illegal value of IOPTRES in CC_KMAT.')
      END IF
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CCMM_TGB:
C----------------------------------------------
C
      IF (LOSPC) RETURN
C
*=====================================================================*
* calculate K matrix transformations:
*=====================================================================*
      IADRTH = 1
      DO ITRAN = 1, NBTRAN

        IDLSTA = IBTRAN(1,ITRAN)
        IDLSTB = IBTRAN(2,ITRAN)

        ISYMA  = ILSTSYM(LISTA,IDLSTA)
        ISYMB  = ILSTSYM(LISTB,IDLSTB)
        ISYMAB = MULD2H(ISYMA,ISYMB)

        TIMTRN = SECOND()

*---------------------------------------------------------------------*
* allocate work space for the result vector:
*---------------------------------------------------------------------*
        IF (CCS) THEN
          KTHETA1 = 1
          KTHETA2 = KDUM
          KEND1   = KTHETA1 + NT1AM(ISYMAB)
          LWRK1 = LWORK - KEND1
          CALL DZERO(WORK(KTHETA1),NT1AM(ISYMAB))
        ELSE 
          KTHETA1 = 1 
          KTHETA2 = KTHETA1 + NT1AM(ISYMAB)
          KEND1   = KTHETA2 + NT2AM(ISYMAB)
          LWRK1 = LWORK - KEND1
          CALL DZERO(WORK(KTHETA1),NT1AM(ISYMAB))
          CALL DZERO(WORK(KTHETA2),NT2AM(ISYMAB))
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'K matrix transformation for ITRAN,',ITRAN
         WRITE (LUPRI,*) 'IADRTH:',IADRTH
         WRITE (LUPRI,*) 'LISTA,IDLSTA:',LISTA,IDLSTA
         WRITE (LUPRI,*) 'LISTB,IDLSTB:',LISTB,IDLSTB
         WRITE (LUPRI,*) 'ISYMA,ISYMB,ISYMAB:',ISYMA,ISYMB,ISYMAB
         CALL FLSHFO(LUPRI)
        END IF
C
C 
        IF (NYQMMM) THEN
          
          RSPTYP='K'
          CALL CCMM_QRTRANSFORMER(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,
     *                         LISTB,IDLSTB,ISYMB,
     *                         LISTA,IDLSTA,ISYMA,
     *                         MODELW,RSPTYP,WORK(KEND1),LWRK1)
C
        ELSE IF ((.NOT. NYQMMM) .AND. (.NOT. USE_PELIB())) THEN
          KTGB  = KEND1 
          KEND2 = KTGB + N2BST(ISYMB)
          LWRK2 = LWORK   - KEND2 
          IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CC_KMAT, 1')
C
          CALL DZERO(WORK(KTGB),N2BST(ISYMB))
C
C------------------------------------------------
C       Trial vector (B left) one excitation part
C------------------------------------------------
C
          KT1AMPB = KEND2
          KEND3   = KT1AMPB +  NT1AM(ISYMB)
          LWRK3   = LWORK -  KEND3 
          IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insuff. work in CC_KMAT 2')
          END IF
          CALL DZERO(WORK(KT1AMPB),NT1AM(ISYMB))
C
          IOPT = 1
          CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     *                WORK(KT1AMPB),WORK(KDUM))
C
          IF (.NOT. CCMM ) CALL CCSL_TGB(WORK(KT1AMPB),ISYMB,
     *                                 LISTB,IDLSTB,WORK(KTGB),
     *                                 'XI',MODEL,
     *                                 WORK(KEND3),LWRK3)
    
C
          IF (CCMM) CALL CCMM_TGB(WORK(KT1AMPB),ISYMB,
     *                          LISTB,IDLSTB,WORK(KTGB),
     *                         'XI',MODEL,
     *                          WORK(KEND3),LWRK3)

          NAMPF   = NT1AM(ISYMAB) + NT2AM(ISYMAB)
C
          KETA    = KEND2
          KEND3   = KETA + NAMPF
          LWRK3   = LWORK   - KEND3 
          IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CC_KMAT, 3')
          CALL DZERO(WORK(KETA),NAMPF)
C
          IF ( LOCDBG ) THEN
            TAL1= DDOT(N2BST(ISYMAB),WORK(KTGB),1,WORK(KTGB),1)
            WRITE(LUPRI,*) 'Printing norm2 TGCMO in KMAT: ', TAL1
          END IF

          LABEL = 'GIVE INT'
          CALL CC_ETAC(ISYMB,LABEL,WORK(KETA),
     *                 LISTA,IDLSTA,0,WORK(KTGB),WORK(KEND3),LWRK3)

C
          KETA1   = KETA
          KETA2   = KETA + NT1AM(ISYMAB)

          IF (LOCDBG) THEN
            TAL1= DDOT(NT1AM(ISYMAB),WORK(KETA1),1,WORK(KETA1),1)
            TAL2= DDOT(NT2AM(ISYMAB),WORK(KETA2),1,
     *              WORK(KETA2),1)
            WRITE(LUPRI,*) 'Printing first contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
          END IF
C  
          LSAME  =  (LISTA.EQ.LISTB  .AND. IDLSTA.EQ.IDLSTB)
          IF (LSAME) THEN
            FACTSLV = TWO
          ELSE
            FACTSLV = ONE
          END IF
C
          CALL DAXPY(NT1AM(ISYMAB),FACTSLV,WORK(KETA1),1,
     *               WORK(KTHETA1),1)
          CALL DAXPY(NT2AM(ISYMAB),FACTSLV,WORK(KETA2),1,
     *               WORK(KTHETA2),1)
C
C------------------------------------------------------------------
C
          IF (.NOT. (LSAME)) THEN
C
            KTGA  = KEND1
            KEND2 = KTGA + N2BST(ISYMA)
            LWRK2 = LWORK   - KEND2
            IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CC_KMAT, 4')
            CALL DZERO(WORK(KTGA),N2BST(ISYMA))
C
C------------------------------------------------
C         Trial vector (A left) one excitation part
C------------------------------------------------
C
            KT1AMPA = KEND2
            KEND3   = KT1AMPA +  NT1AM(ISYMA)
            LWRK3   = LWORK -  KEND3
            IF (LWRK3 .LT. 0) THEN
              CALL QUIT('Insuff. work in CC_KMAT 5')
            END IF
            CALL DZERO(WORK(KT1AMPA),NT1AM(ISYMA))
C
            IOPT = 1
            CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     *                  WORK(KT1AMPA),WORK(KDUM))
C
            IF (.NOT. CCMM) CALL CCSL_TGB(WORK(KT1AMPA),ISYMA,LISTA,
     *                    IDLSTA,WORK(KTGA),'XI',
     *                    MODEL,WORK(KEND3),LWRK3)
C
            IF (CCMM) CALL CCMM_TGB(WORK(KT1AMPA),ISYMA,LISTA,
     *                   IDLSTA,WORK(KTGA),'XI',
     *                   MODEL,WORK(KEND3),LWRK3)
C
            NAMPF   = NT1AM(ISYMAB) + NT2AM(ISYMAB)
C
            KETA    = KEND2
            KEND3   = KETA + NAMPF
            LWRK3   = LWORK   - KEND3
            IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CC_KMAT, 6')
            CALL DZERO(WORK(KETA),NAMPF)
C
            LABEL = 'GIVE INT'
            CALL CC_ETAC(ISYMA,LABEL,WORK(KETA),
     *                   LISTB,IDLSTB,0,WORK(KTGA),WORK(KEND3),LWRK3)
C
            KETA1   = KETA
            KETA2   = KETA + NT1AM(ISYMAB)

            IF (LOCDBG) THEN
              TAL1= DDOT(NT1AM(ISYMAB),WORK(KETA1),1,WORK(KETA1),1)
              TAL2= DDOT(NT2AM(ISYMAB),WORK(KETA2),1,
     *              WORK(KETA2),1)
              WRITE(LUPRI,*) 'Printing second contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
            END IF
C
            CALL DAXPY(NT1AM(ISYMAB),ONE,WORK(KETA1),1,WORK(KTHETA1),1)
            CALL DAXPY(NT2AM(ISYMAB),ONE,WORK(KETA2),1,WORK(KTHETA2),1)
C
          END IF

        END IF ! NYQMMM
C
        IF (USE_PELIB()) THEN
          RSPTYP='K'
          CALL PELIB_IFC_QRTRANSFORMER(WORK(KTHETA1),WORK(KTHETA2),
     &                   ISYMAB,LISTB,IDLSTB,ISYMB,LISTA,IDLSTA,ISYMA,
     &                   MODELW,RSPTYP,WORK(KEND1),LWRK1)
        END IF

*---------------------------------------------------------------------*
* write result vector to output:
*---------------------------------------------------------------------*

        IF (IOPTRES .EQ. 0  .OR. IOPTRES .EQ. 1) THEN

*       write to a common direct access file, 
*       store start address in IBTRAN(3,ITRAN)

        IBTRAN(3,ITRAN) = IADRTH

        CALL PUTWA2(LUBMAT,FILBMA,WORK(KTHETA1),IADRTH,NT1AM(ISYMAB))
        IADRTH = IADRTH + NT1AM(ISYMAB)

        IF (.NOT.CCS) THEN
          CALL PUTWA2(LUBMAT,FILBMA,WORK(KTHETA2),IADRTH,NT2AM(ISYMAB))
          IADRTH = IADRTH + NT2AM(ISYMAB)
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'K matrix transformation nb. ',ITRAN,
     &          ' saved on file.'
         WRITE (LUPRI,*) 'ADRESS, LENGTH:',
     &        IBTRAN(3,ITRAN),IADRTH-IBTRAN(3,ITRAN)
         XNORM = DDOT(NT1AM(ISYMAB),WORK(KTHETA1),1,WORK(KTHETA1),1)
         IF (.NOT.CCS) XNORM = XNORM +
     &           DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
         WRITE (LUPRI,*) 'Norm:', XNORM

         Call AROUND('K matrix transformation written to file:')
         Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,1,1)
        END IF

      ELSE IF ( IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4 ) THEN

*        write to a sequential file by a call to CC_WRRSP/CC_WARSP,
*        use FILBMA as LIST type and IBTRAN(3,ITRAN) as index
         KTHETA0 = -999999
         IF (IOPTRES.EQ.3) THEN
           CALL CC_WRRSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND1),LWRK1)
         ELSE IF (IOPTRES.EQ.4) THEN
           CALL CC_WARSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND1),LWRK1)
         END IF

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Write K * ',LISTA,' * ',LISTB,
     &              ' transformation',
     &              ' as ',FILBMA,' type vector to file.'
           WRITE (LUPRI,*) 'index of inp. A vector:',IBTRAN(1,ITRAN)
           WRITE (LUPRI,*) 'index of inp. B vector:',IBTRAN(2,ITRAN)
           WRITE (LUPRI,*) 'index of result vector:',IBTRAN(3,ITRAN)
           LEN = NT1AM(ISYMAB) + NT2AM(ISYMAB)
           IF (CCS) LEN = NT1AM(ISYMAB)
           XNORM = DDOT(LEN,WORK(KTHETA1),1,WORK(KTHETA1),1)
           WRITE (LUPRI,*) 'norm^2 of result vector:',XNORM
         END IF
      ELSE IF (IOPTRES.EQ.5) THEN
         IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KTHETA2),TWO,ISYMAB)
         CALL CCDOTRSP(IBDOTS,BCONS,IOPTW,FILBMA,ITRAN,NBTRAN,MXVEC,
     &                 WORK(KTHETA1),WORK(KTHETA2),ISYMAB,
     &                 WORK(KEND1),LWRK1)
      ELSE
        CALL QUIT('Illegal value for IOPTRES in CC_KMAT.')
      END IF

      TIMTRN = SECOND() - TIMTRN
      
      IF (IPRINT.GT.0) THEN

         IF (IOPTRES.EQ.5) THEN
            IVEC = 1
            DO WHILE (IBDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
               IVEC = IVEC + 1
            END DO    
            WRITE (LUPRI,'(1X,2(A,I5),A,I6,A,F10.2,A)')'| ',IDLSTA,
     &        '    | ',IDLSTB,'    | ',IVEC-1,'       | ',TIMTRN,'  |'
         ELSE
            WRITE (LUPRI,'(1X,2(A,I5),A,I6,A,F10.2,A)') '| ',IDLSTA, 
     &           '    | ',
     &        IDLSTB,'    | ',IBTRAN(3,ITRAN),'       | ',TIMTRN,'  |'
         END IF 

      END IF

*---------------------------------------------------------------------*
* End of loop over K matrix transformations
*---------------------------------------------------------------------*
      END DO
      WRITE (LUPRI,'(1X,A1,50("="),A1,//)') '+','+'

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

        CALL GETWA2(LUBMAT,FILBMA,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
          DO ITRAN = 1, NBTRAN
            IF (ITRAN.LT.NBTRAN) THEN
              LEN     = IBTRAN(3,ITRAN+1)-IBTRAN(3,ITRAN)
            ELSE
              LEN     = IADRTH-IBTRAN(3,NBTRAN)
            END IF
            KTHETA1 = IBTRAN(3,ITRAN)
            XNORM   = DDOT(LEN, WORK(KTHETA1),1, WORK(KTHETA1),1)
            WRITE (LUPRI,*) 'Read K matrix transformation nb. ',NBTRAN
            WRITE (LUPRI,*) 'Adress, length, NORM:',IBTRAN(3,NBTRAN),
     &                      LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

*---------------------------------------------------------------------*
* close K matrix file, print timings & return
*---------------------------------------------------------------------*

      IF (IOPTRES.EQ.0 ) THEN
        CALL WCLOSE2(LUBMAT, FILBMA, 'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
        CALL WCLOSE2(LUBMAT, FILBMA, 'DELETE')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_KMAT.')
      END IF

*=====================================================================*

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CC_KMAT
*=====================================================================*

