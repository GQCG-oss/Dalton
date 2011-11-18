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
*=====================================================================*
      SUBROUTINE CC_JMATRIX(IJTRAN, NJTRAN, LISTL, IOPTRES,
     &                      FILJMA, IJDOTS, CJCON, MXVEC, WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over J matrix transformations
*             (needed if the number of transformations exceeds the
*              limit MAXSIM defined on ccsdio.h )
*
*        
*     Written by Christof Haettig, November 1998, based on CC_BMATRIX.
*     modified by JK in order to perform J matrix transformations (03)
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include "maxorb.h"
#include "ccsdio.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*(*) LISTL, FILJMA
      INTEGER IOPTRES
      INTEGER NJTRAN, MXVEC, LWORK
      INTEGER IJTRAN(3,NJTRAN)
      INTEGER IJDOTS(MXVEC,NJTRAN)
      
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL CJCON(MXVEC,NJTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION CJCON(MXVEC,NJTRAN)
#endif 

      INTEGER MAXJTRAN, NTRAN, ISTART, IBATCH, NBATCH

      MAXJTRAN = MAXSIM

      NBATCH = (NJTRAN+MAXJTRAN-1)/MAXJTRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Batching over J matrix transformations:'
        WRITE (LUPRI,*) 'nb. of batches needed:', NBATCH
      END IF
  
      DO IBATCH = 1, NBATCH
        ISTART = (IBATCH-1) * MAXJTRAN + 1
        NTRAN  = MIN(NJTRAN-(ISTART-1),MAXJTRAN)

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',ISTART
          WRITE (LUPRI,*) '# transf.:',NTRAN
        END IF

        CALL CC_JMAT(IJTRAN(1,ISTART), NTRAN,
     &               LISTL, IOPTRES, FILJMA, 
     &               IJDOTS(1,ISTART), CJCON(1,ISTART), 
     &               MXVEC, WORK, LWORK)

      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_JMATRIX                           *
*---------------------------------------------------------------------*
c/* Deck CC_JMAT */
*=====================================================================*
      SUBROUTINE CC_JMAT(IJTRAN, NJTRAN, LISTA, IOPTRES,
     &                    FILJMA, IJDOTS, CJCON, MXVEC, WORK, LWORK )
*---------------------------------------------------------------------*
*             The linear transformations are calculated for a list
*             of bar{T^A} vectors.  
*
*                LISTA       -- type of bar{T^A} vectors
*                IJTRAN(1,*) -- indeces of {T^B} vectors (0)
*                IJTRAN(2,*) -- indeces of \bar{T^A} vectors
*                IJTRAN(3,*) -- indeces or addresses of result vectors
*                NJTRAN      -- number of requested transformations
*                FILJMA      -- file name / list type of result vectors
*                               or list type of vectors to be dotted on
*                IJDOTS      -- indeces of vectors to be dotted on
*                CJCON       -- contains the dot products on return
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  all result vectors are written to a direct
*                          access file, FILJMA is used as file name
*                          the start addresses of the vectors are
*                          returned in IJTRAN(3,*)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          IJTRAN(3,*). N.B.: if WORK is not large
*                          enough IOPTRES is automatically reset to 0!!
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, FILJMA is used
*                          as list type and IJTRAN(3,*) as list index
*                          NOTE that IJTRAN(3,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WARSP, FILJMA is used
*                          as list type and IJTRAN(3,*) as list index
*                          NOTE that IJTRAN(3,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by FILJMA and the indeces from IJDOTS
*                          the result of the dot products is returned
*                          in the CJCON array
*
*
*           CCMM JK, modyfied version of CC_FMAT
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
#include "ccsections.h"
#include "maxorb.h"
#include "mxcent.h"
#include "ccsdio.h"
#include "ccorb.h"
#include "iratdef.h"
#include "eribuf.h"
#include "ccslvinf.h"
#include "second.h"
#include "qm3.h"

* local parameters:
      CHARACTER MSGDBG*(16)
      PARAMETER (MSGDBG='[debug] CC_KMAT> ')

      LOGICAL LOCDBG, LSAME
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space

      
      INTEGER LUBMAT

      CHARACTER*(*) LISTA, FILJMA

      INTEGER IOPTRES
      INTEGER NJTRAN, MXVEC, LWORK
      INTEGER IJTRAN(3,NJTRAN)
      INTEGER IJDOTS(MXVEC,NJTRAN)

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL ZERO, ONE, TWO
      REAL DUM, XNORM, DUMMY
      REAL CJCON(MXVEC,NJTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ZERO, ONE, TWO
      DOUBLE PRECISION DUM, XNORM, DUMMY
      DOUBLE PRECISION CJCON(MXVEC,NJTRAN)
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0)

      CHARACTER*(10) MODEL, MODELW
      CHARACTER*8 LABEL
      INTEGER NBATCH


      INTEGER ITRAN, IDLSTA, IOPT
      INTEGER ISYMA
      INTEGER IBATCH,IADRTH
      INTEGER KEND1, LEN, LENALL
      INTEGER KEND2, LWRK2, KEND3, LWRK3, LWRK1
      INTEGER IVEC
      INTEGER KT2AMPA, KTHETA0
      INTEGER KTHETA1, KTHETA2, KT1AMPA, KT1AMPB 
      INTEGER IOPTW, IDUMMY
      INTEGER KTGB, KXI, KXI1, KXI2, NAMPF

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
        Call AROUND('ENTERED CC_JMAT')
        WRITE (LUPRI,*) 'LISTA : ',LISTA
        WRITE (LUPRI,*) 'FILJMA: ',FILJMA
        WRITE (LUPRI,*) 'NJTRAN: ',NJTRAN
        WRITE (LUPRI,*) 'IOPTRES:',IOPTRES
        CALL FLSHFO(LUPRI)
      END IF
      
      IF (CCSDT) THEN
        WRITE(LUPRI,'(/1x,a)') 'J matrix transformations not '
     &          //'implemented for triples yet...'
        CALL QUIT('Triples not implemented for J matrix '//
     &            'transformations')
      END IF

      IF ( .not. (CCS .or. CC2 .or. CCSD) ) THEN
        WRITE(LUPRI,'(/1x,a)') 'CC_JMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_JMAT...'
        CALL QUIT('Unknown CC method in CC_JMAT.')
      END IF

      IF (.NOT. DUMPCD) THEN
        WRITE(LUPRI,*) 'DUMPCD = ',DUMPCD
        WRITE(LUPRI,*) 'CC_JMAT requires DUMPCD=.TRUE.'
        CALL QUIT('DUMPCD=.FALSE. , CC_JMAT requires DUMPCD=.TRUE.')
      END IF

      IF (.NOT. RSPIM) THEN
        WRITE(LUPRI,*) 'RSPIM = ',RSPIM
        WRITE(LUPRI,*) 'CC_JMAT requires RSPIM=.TRUE.'
        CALL QUIT('RSPIM=.FALSE. , CC_JMAT requires RSPIM=.TRUE.')
      END IF

      IF (ISYMOP .NE. 1) THEN
        WRITE(LUPRI,*) 'ISYMOP = ',ISYMOP
        WRITE(LUPRI,*) 'CC_JMAT is not implemented for ISYMOP.NE.1'
        CALL QUIT('CC_JMAT is not implemented for ISYMOP.NE.1')
      END IF

      IF (NJTRAN .GT. MAXSIM) THEN
        WRITE(LUPRI,*) 'NJTRAN = ', NJTRAN
        WRITE(LUPRI,*) 'MAXSIM = ', MAXSIM
        WRITE(LUPRI,*) 'number of requested transformation is larger'
        WRITE(LUPRI,*) 'than the maximum number of allowed ',
     &                 'simultaneous transformation.'
        WRITE(LUPRI,*) 'Error in CC_JMAT: NJTRAN is larger than MAXSIM.'
        CALL QUIT('Error in CC_JMAT: NJTRAN is larger than MAXSIM.')
      END IF

      IF (IPRINT.GT.0) THEN
 
         WRITE (LUPRI,'(//1X,A1,50("="),A1)')'+','+'

         WRITE (LUPRI,'(1x,A52)')
     &         '|        J MATRIX TRANSFORMATION SECTION           |'

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
         CALL QUIT('Unknown coupled cluster model in CC_JMAT.')
      END IF


* check return option for the result vectors:
      LUBMAT = -1
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
         CALL WOPEN2(LUBMAT, FILJMA, 64, 0)
      ELSE IF (IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4) THEN
         CONTINUE
      ELSE IF (IOPTRES .EQ. 5) THEN
         IF (MXVEC*NJTRAN.NE.0) CALL DZERO(CJCON,MXVEC*NJTRAN)
      ELSE
         CALL QUIT('Illegal value of IOPTRES in CC_JMAT.')
      END IF
C
C----------------------------------------------
C     If all models are SPC
C     -> RETURN from CC_JMAT:
C----------------------------------------------
C
      IF (LOSPC) RETURN
C
*=====================================================================*
* calculate J matrix transformations:
*=====================================================================*
      IADRTH = 1
      DO ITRAN = 1, NJTRAN

        IDLSTA = IJTRAN(2,ITRAN)
        ISYMA  = ILSTSYM(LISTA,IDLSTA)

        TIMTRN = SECOND()

*---------------------------------------------------------------------*
* allocate work space for the result vector:
*---------------------------------------------------------------------*
        IF (CCS) THEN
          KTHETA1 = 1
          KTHETA2 = KDUM
          KEND1   = KTHETA1 + NT1AM(ISYMA)
          LWRK1 = LWORK - KEND1
          CALL DZERO(WORK(KTHETA1),NT1AM(ISYMA))
        ELSE 
          KTHETA1 = 1 
          KTHETA2 = KTHETA1 + NT1AM(ISYMA)
          KEND1   = KTHETA2 + NT2AM(ISYMA)
          LWRK1 = LWORK - KEND1
          CALL DZERO(WORK(KTHETA1),NT1AM(ISYMA))
          CALL DZERO(WORK(KTHETA2),NT2AM(ISYMA))
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'J matrix transformation for ITRAN,',ITRAN
         WRITE (LUPRI,*) 'IADRTH:',IADRTH
         WRITE (LUPRI,*) 'LISTA,IDLSTA:',LISTA,IDLSTA
         WRITE (LUPRI,*) 'ISYMA:',ISYMA
         CALL FLSHFO(LUPRI)
        END IF
C
        KTGB  = KEND1 
        KEND2 = KTGB + N2BST(ISYMA)
        LWRK2 = LWORK   - KEND2 
        IF (LWRK2.LT.0) CALL QUIT( 'Too little work in CC_JMAT, 1')
C
        CALL DZERO(WORK(KTGB),N2BST(ISYMA))
C
C------------------------------------------------
C       Trial vector (A left) one excitation part
C------------------------------------------------
C
        KT1AMPA = KEND2
        KEND3   = KT1AMPA +  NT1AM(ISYMA)
        LWRK3   = LWORK -  KEND3 
        IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insuff. work in CC_JMAT 2')
        END IF
        CALL DZERO(WORK(KT1AMPA),NT1AM(ISYMA))
C
        IOPT = 1
        CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     *                WORK(KT1AMPA),WORK(KDUM))
C
        IF (.NOT. CCMM ) CALL CCSL_TGB(WORK(KT1AMPA),ISYMA,
     *                                 LISTA,IDLSTA,WORK(KTGB),
     *                                 'XI',MODEL,
     *                                 WORK(KEND3),LWRK3)
    
C
        IF (CCMM) CALL CCMM_TGB(WORK(KT1AMPA),ISYMA,
     *                          LISTA,IDLSTA,WORK(KTGB),
     *                         'XI',MODEL,
     *                          WORK(KEND3),LWRK3)

        NAMPF   = NT1AM(ISYMA) + NT2AM(ISYMA)
C
        KXI    = KEND2
        KEND3  = KXI + NAMPF
        LWRK3  = LWORK   - KEND3 
        IF (LWRK3.LT.0) CALL QUIT( 'Too little work in CC_JMAT, 3')
        CALL DZERO(WORK(KXI),NAMPF)
C
        LABEL = 'GIVE INT'
        CALL CC_XKSI(WORK(KXI),LABEL,ISYMA,0,WORK(KTGB),
     *               WORK(KEND3),LWRK3)
C
        KXI1   = KXI
        KXI2   = KXI + NT1AM(ISYMA)
C
        CALL DAXPY(NT1AM(ISYMA),ONE,WORK(KXI1),1,WORK(KTHETA1),1)
        CALL DAXPY(NT2AM(ISYMA),ONE,WORK(KXI2),1,WORK(KTHETA2),1)
C
*---------------------------------------------------------------------*
* write result vector to output:
*---------------------------------------------------------------------*

      IF (IOPTRES .EQ. 0  .OR. IOPTRES .EQ. 1) THEN

*       write to a common direct access file, 
*       store start address in IJTRAN(3,ITRAN)

        IJTRAN(3,ITRAN) = IADRTH

        CALL PUTWA2(LUBMAT,FILJMA,WORK(KTHETA1),IADRTH,NT1AM(ISYMA))
        IADRTH = IADRTH + NT1AM(ISYMA)

        IF (.NOT.CCS) THEN
          CALL PUTWA2(LUBMAT,FILJMA,WORK(KTHETA2),IADRTH,NT2AM(ISYMA))
          IADRTH = IADRTH + NT2AM(ISYMA)
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'J matrix transformation nb. ',ITRAN,
     &          ' saved on file.'
         WRITE (LUPRI,*) 'ADRESS, LENGTH:',
     &        IJTRAN(3,ITRAN),IADRTH-IJTRAN(3,ITRAN)
         XNORM = DDOT(NT1AM(ISYMA),WORK(KTHETA1),1,WORK(KTHETA1),1)
         IF (.NOT.CCS) XNORM = XNORM +
     &           DDOT(NT2AM(ISYMA),WORK(KTHETA2),1,WORK(KTHETA2),1)
         WRITE (LUPRI,*) 'Norm:', XNORM

         Call AROUND('J matrix transformation written to file:')
         Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMA,1,1)
        END IF

      ELSE IF ( IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4 ) THEN

*        write to a sequential file by a call to CC_WRRSP/CC_WARSP,
*        use FILJMA as LIST type and IJTRAN(3,ITRAN) as index
         KTHETA0 = -999999
         IF (IOPTRES.EQ.3) THEN
           CALL CC_WRRSP(FILJMA,IJTRAN(3,ITRAN),ISYMA,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND1),LWRK1)
         ELSE IF (IOPTRES.EQ.4) THEN
           CALL CC_WARSP(FILJMA,IJTRAN(3,ITRAN),ISYMA,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND1),LWRK1)
         END IF

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Write J ',LISTA,
     &              ' transformation',
     &              ' as ',FILJMA,' type vector to file.'
           WRITE (LUPRI,*) 'index of inp. A vector:',IJTRAN(2,ITRAN)
           WRITE (LUPRI,*) 'index of result vector:',IJTRAN(3,ITRAN)
           LEN = NT1AM(ISYMA) + NT2AM(ISYMA)
           IF (CCS) LEN = NT1AM(ISYMA)
           XNORM = DDOT(LEN,WORK(KTHETA1),1,WORK(KTHETA1),1)
           WRITE (LUPRI,*) 'norm^2 of result vector:',XNORM
         END IF
      ELSE IF (IOPTRES.EQ.5) THEN
         IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KTHETA2),TWO,ISYMA)
         CALL CCDOTRSP(IJDOTS,CJCON,IOPTW,FILJMA,ITRAN,NJTRAN,MXVEC,
     &                 WORK(KTHETA1),WORK(KTHETA2),ISYMA,
     &                 WORK(KEND1),LWRK1)
      ELSE
        CALL QUIT('Illegal value for IOPTRES in CC_JMAT.')
      END IF

      TIMTRN = SECOND() - TIMTRN
      
      IF (IPRINT.GT.0) THEN

         IF (IOPTRES.EQ.5) THEN
            IVEC = 1
            DO WHILE (IJDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
               IVEC = IVEC + 1
            END DO    
            WRITE (LUPRI,'(1X,2(A,I5),A,F10.2,A)')'| ',IDLSTA,
     &        '    |  ',IVEC-1,'       | ',TIMTRN,'  |'
         ELSE
            WRITE (LUPRI,'(1X,2(A,I5),A,F10.2,A)') '| ',IDLSTA, 
     &           '    | ',
     &        '    | ',IJTRAN(3,ITRAN),'       | ',TIMTRN,'  |'
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

        CALL GETWA2(LUBMAT,FILJMA,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
          DO ITRAN = 1, NJTRAN
            IF (ITRAN.LT.NJTRAN) THEN
              LEN     = IJTRAN(3,ITRAN+1)-IJTRAN(3,ITRAN)
            ELSE
              LEN     = IADRTH-IJTRAN(3,NJTRAN)
            END IF
            KTHETA1 = IJTRAN(3,ITRAN)
            XNORM   = DDOT(LEN, WORK(KTHETA1),1, WORK(KTHETA1),1)
            WRITE (LUPRI,*) 'Read J matrix transformation nb. ',NJTRAN
            WRITE (LUPRI,*) 'Adress, length, NORM:',IJTRAN(3,NJTRAN),
     &                      LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

*---------------------------------------------------------------------*
* close J matrix file, print timings & return
*---------------------------------------------------------------------*

      IF (IOPTRES.EQ.0 ) THEN
        CALL WCLOSE2(LUBMAT, FILJMA, 'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
        CALL WCLOSE2(LUBMAT, FILJMA, 'DELETE')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_JMAT.')
      END IF

*=====================================================================*

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CC_JMAT
*=====================================================================*

