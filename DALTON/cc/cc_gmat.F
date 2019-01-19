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
c/* Deck CC_GMATRIX */
*=====================================================================*
       SUBROUTINE CC_GMATRIX( LISTA,  ! inp: A Zeta vector list
     &                        LISTB,  ! inp: B amplitude list
     &                        LISTC,  ! inp: C amplitude list
     &                        LISTD,  ! inp: D amplitude list
     &                        NGTRAN, ! inp: nb. of G matrix transf.
     &                        MXVEC,  ! inp: max nb. of dot products
     &                        IGTRAN, ! inp: indeces of G mat. transf.
     &                        IGDOTS, ! inp: indeces of dot products
     &                        GCON,   ! inp: array for dot products
     &                        WORK,   ! scr: work space
     &                        LWORK,  ! inp: length of work space
     &                        IOPTRES)! inp: output option     
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over G matrix transformations
*
*             for more detailed documentation see: CC_GMAT
*        
*     Written by Christof Haettig, April 1999, based on CC_FMATRIX.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "maxorb.h"
#include "ccslvinf.h"
#include "qm3.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MAXGTRAN
      PARAMETER (MAXGTRAN = 100)

      CHARACTER*(*) LISTA, LISTB, LISTC, LISTD
      INTEGER IOPTRES, NGTRAN, MXVEC, LWORK
      INTEGER IGTRAN(4,NGTRAN)
      INTEGER IGDOTS(MXVEC,NGTRAN)
      
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL GCON(MXVEC,NGTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION GCON(MXVEC,NGTRAN)
#endif 

      INTEGER NTRAN, ISTART, IBATCH, NBATCH
      INTEGER KTIMTRN, KDELTA1, KT1AMPB, KLAMDHB, KLAMDPB, KT1AMPC  
      INTEGER KLAMDHC, KLAMDPC, KCTR1, KCTR2, KFOCKB, KFOCKC  
      INTEGER KFBCOO, KFBCVV  
      INTEGER KEND, LEND

      CALL QENTER('CC_GMATRIX')

      NBATCH = (NGTRAN+MAXGTRAN-1)/MAXGTRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Batching over G matrix transformations:'
        WRITE (LUPRI,*) 'nb. of batches needed:', NBATCH
      END IF
  
      DO IBATCH = 1, NBATCH
        ISTART = (IBATCH-1) * MAXGTRAN + 1
        NTRAN  = MIN(NGTRAN-(ISTART-1),MAXGTRAN)

        KTIMTRN  = 1
        KDELTA1  = KTIMTRN  + NTRAN
        KT1AMPB  = KDELTA1  + NTRAN
        KLAMDHB  = KT1AMPB  + NTRAN
        KLAMDPB  = KLAMDHB  + NTRAN
        KT1AMPC  = KLAMDPB  + NTRAN
        KLAMDHC  = KT1AMPC  + NTRAN
        KLAMDPC  = KLAMDHC  + NTRAN
        KCTR1    = KLAMDPC  + NTRAN
        KCTR2    = KCTR1    + NTRAN
        KFOCKB   = KCTR2    + NTRAN
        KFOCKC   = KFOCKB   + NTRAN
        KFBCOO   = KFOCKC   + NTRAN
        KFBCVV   = KFBCOO   + NTRAN
        KEND     = KFBCVV   + NTRAN
        LEND     = LWORK    - KEND

        IF (LEND .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient work space in CC_GMATRIX.'
           WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
           WRITE (LUPRI,*) 'Need at least:',KEND, ' words.'
           CALL QUIT('Insufficient work space in CC_GMATRIX.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',ISTART
          WRITE (LUPRI,*) '# transf.:',NTRAN
        END IF

        CALL CC_GMAT( LISTA,  LISTB,  LISTC,  LISTD,  NTRAN, MXVEC,  
     &                IGTRAN(1,ISTART),IGDOTS(1,ISTART),GCON(1,ISTART),
     &                WORK(KTIMTRN), WORK(KDELTA1), WORK(KT1AMPB),
     &                WORK(KLAMDHB), WORK(KLAMDPB), WORK(KT1AMPC),
     &                WORK(KLAMDHC), WORK(KLAMDPC), WORK(KCTR1),
     &                WORK(KCTR2),   WORK(KFOCKB),  WORK(KFOCKC),
     &                WORK(KFBCOO),  WORK(KFBCVV),
     &                WORK(KEND),   LEND,  IOPTRES)

      END DO

      CALL QEXIT('CC_GMATRIX')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_GMATRIX                           *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
c/* Deck CC_GMAT */
*=====================================================================*
       SUBROUTINE CC_GMAT( LISTA,  ! inp: A Zeta vector list
     &                     LISTB,  ! inp: B amplitude list
     &                     LISTC,  ! inp: C amplitude list
     &                     LISTD,  ! inp: D amplitude list
     &                     NGTRAN, ! inp: nb. of G matrix transf.
     &                     MXVEC,  ! inp: max nb. of dot products
     &                     IGTRAN, ! inp: index array of G mat. transf.
     &                     IGDOTS, ! inp: index array of dot products
     &                     GCON,   ! inp: array for dot products
     &                     TIMTRN, ! scr: 
     &                     KDELTA1,! scr: 
     &                     KT1AMPB,! scr:
     &                     KLAMDHB,! scr: 
     &                     KLAMDPB,! scr: 
     &                     KT1AMPC,! scr:
     &                     KLAMDHC,! scr: 
     &                     KLAMDPC,! scr: 
     &                     KCTR1,  ! scr:
     &                     KCTR2,  ! scr: 
     &                     KFOCKB, ! scr: 
     &                     KFOCKC, ! scr:
     &                     KFBCOO, ! scr: 
     &                     KFBCVV, ! scr:      
     &                     WORK,   ! scr: work space
     &                     LWORK,  ! inp: length of work space
     &                     IOPTRES)! inp: output option
*---------------------------------------------------------------------*
*
*    Purpose: AO-direct calculation of a linear transformation of two
*             CC amplitude vectors, T^B and T^C, with the CC G matrix
*             (second partial derivative of the CC lagrangian)
*
*             index combinations for L^A, T^B, T^C passed in IGTRAN,
*             indeces for T^D in IGDOTS
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  all result vectors are written to a direct
*                          access file, LISTD is used as file name,
*                          the start addresses of the vectors are
*                          returned in IGTRAN(4,*)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          IGTRAN(4,*). N.B.: if WORK is not large
*                          enough, IOPTRES is automatically reset to 0
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, LISTD is used
*                          as list type and IGTRAN(4,*) as list index
*                          NOTE that IGTRAN(4,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WARSP, LISTD is used
*                          as list type and IGTRAN(4,*) as list index
*                          NOTE that IGTRAN(4,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by LISTD and the indeces from IGDOTS
*                          the result of the dot products is returned
*                          in the GCON array
*
*
*    symmetries/variables:
*              
*           ISYRES : result vector 
*           ISYCTR : lagrangian multipliers 
*           ISYMTB : B response amplitudes 
*           ISYMTC : C response amplitudes
*
*           DELTA1 : single excitation part (total sum)
*
*    Note: it is assumed, that the integrals (ia|jb) are on disk
*
*          used vectors of the size (1/2) V^2 O^2 :
*            CCS:  (ia|jb)
*            CC2:  (ia|jb), CTR2
*            CCSD: (ia|jb), CTR2, TAmpB2, TAmpC2, DELTA2
*
*     Written by Christof Haettig, Aug. 1996; restructed, Sept. 1998.
*
*     included CC-R12: Christian Neiss, november 2005
* 
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "iratdef.h"
#include "cbieri.h"
#include "mxcent.h"
#include "eritap.h"
#include "maxorb.h"
#include "distcl.h"
#include "ccorb.h"
#include "ccisao.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "aovec.h"
#include "blocks.h"
#include "second.h"
#include "ccnoddy.h"
#include "r12int.h"
#include "ccr12int.h"
#include "ccsections.h"
#include "ccslvinf.h"
#include "qm3.h"

* local parameters:
      CHARACTER MSGDBG*(16)
      PARAMETER (MSGDBG='[debug] CC_GMAT>')
#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO, THREE
#else
      DOUBLE PRECISION ZERO, ONE, TWO, THREE
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, THREE = 3.0d0)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      LOGICAL LDOUBL
      PARAMETER (LDOUBL = .TRUE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
      INTEGER LUGMAT

      CHARACTER*(*) LISTA, LISTB, LISTC, LISTD
      CHARACTER*6 FILGMA
      INTEGER MXVEC, NGTRAN, IOPTRES, LWORK
      INTEGER IGTRAN(4,NGTRAN)
      INTEGER IGDOTS(MXVEC,NGTRAN)

      CHARACTER*(1) RSPTYP
      CHARACTER*(8) FNTOC
      CHARACTER*(7) FN3FOP2X
      CHARACTER*(6) FNCKJD, FN3VI, FN3FOPX
      CHARACTER*(5) FNDKBC3
      INTEGER LUTOC,LU3FOP2X,LUCKJD, LU3VI, LU3FOPX,LUDKBC3
      PARAMETER (FNTOC ='CCSDT_OC',FN3FOP2X='PTFOP2X',
     *           FNCKJD='CKJDEL',FN3VI ='CC3_VI'  ,
     *           FN3FOPX='PTFOPX' ,FNDKBC3='DKBC3')


#if defined (SYS_CRAY)
      REAL GCON(MXVEC,NGTRAN)
      REAL WORK(LWORK)
      REAL DDOT, FACT, XNORM, TIMTRN(NGTRAN)
      REAL DTIME, TIMALL, TIMFCK, TIMIO, TIMCCS, TIMRDAO
      REAL TIM21CD, TIM21SD, TIMDBL, TIMINT, TIMTRBT,CONVRT
      REAL TIMTRIP
#else
      DOUBLE PRECISION GCON(MXVEC,NGTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION DDOT, FACT, XNORM, TIMTRN(NGTRAN)
      DOUBLE PRECISION DTIME, TIMALL, TIMFCK, TIMIO, TIMCCS, TIMRDAO
      DOUBLE PRECISION TIM21CD, TIM21SD, TIMDBL, TIMINT, TIMTRBT,CONVRT
      DOUBLE PRECISION TIMTRIP
#endif

      INTEGER KT1AMP0, KT1AMPB(NGTRAN), KT1AMPC(NGTRAN), KCTR1(NGTRAN) 
      INTEGER KLAMDP0, KLAMDPB(NGTRAN), KLAMDPC(NGTRAN)
      INTEGER KLAMDH0, KLAMDHB(NGTRAN), KLAMDHC(NGTRAN)
      INTEGER KFOCKB(NGTRAN), KFOCKC(NGTRAN)
      INTEGER KFBCOO(NGTRAN), KFBCVV(NGTRAN)
      INTEGER KDELTA1(NGTRAN), KCTR2(NGTRAN)

      CHARACTER MODEL*(10), MODELW*(10), APROXR12*(3)
      LOGICAL LSAME, LLSAME, CCSDR12
      INTEGER INDEXA(MXCORB_CC)
      INTEGER IDLSTA, IDLSTB, IDLSTC, IDLSTD
      INTEGER ISYRES, ISYCTR, ISYMTB, ISYMTC, ISYMTD
      INTEGER KEND0, KEND1, KEND1A, KEND2, KFREE, KENDSV, KENDE3
      INTEGER LEND0, LEND1, LEND1A, LEND2, LFREE, LWRKSV, LENDE3
      INTEGER KENDB2, KENDB3, KENDA2, KENDA3, KEND3, KXIAJB,KC4O,KX4O
      INTEGER LENDB2, LENDB3, LENDA2, LENDA3, LEND3, KDELTA2, KCDTERM
      INTEGER ISYTCTB, ISYOVOV, ISYMFB, ISYMFC, ISYFBC, ISYDEL
      INTEGER KLIAJB, KSCR, LSCR, LWRK0, ITRAN, ISTRT, IEND, IFILE
      INTEGER KTTZET, KTZBOO, KTZBVV, KXINT, KDUM, KZETA2, IOPTW
      INTEGER KT2AMPB, KT2AMPC, KEMAT1, KEMAT2
      INTEGER KDELTBC, KDELTCB, KINDXB, KCCFB1
      INTEGER KODCL1, KODCL2, KODBC1, KODBC2, KRDBC1, KRDBC2
      INTEGER KODPP1, KODPP2, KRDPP1, KRDPP2, KRECNR
      INTEGER NTOSYM,ISYMD1,NTOT,ILLL,NUMDIS,IDEL,IDEL2, LENALL,IVEC
      INTEGER ISYDIS,ISYC4O,ISYX4O,IOPT,ADR,IERR,IADRDEL,LEN1,LEN2
      INTEGER ISYMI, ISYMJ, ISYMA, ISYMB, NAB, NBA, NIJ, NJI
      INTEGER KDSRHF0,KDSRHFBC,ISYMXB,ISYMXC,ISYCTB,ISYCTC,IOPTWE
      INTEGER ITAIKJ(8,8), NTAIKJ(8), ICOUNT, ISYM, ISYMAIK
      INTEGER KXBIAJK, KEBIAJK, KZBIAJK, KCBIAJK, KXLMDHB, KXLMDPB
      INTEGER KXCIAJK, KECIAJK, KZCIAJK, KCCIAJK, KXLMDHC, KXLMDPC
      INTEGER MDISAO,MDSRHF,M2BST,MT2BGD,MT2SQ,MT2AM,MSCRATCH,MFREE

      INTEGER ILSTSYM, IOPTTCME, IOPTWR12, LENMOD
      INTEGER IDUM
      INTEGER KDELTAR12, KR12SCR, KVABKL, LUNIT, LENR12

      CALL QENTER('CC_GMAT')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_GMAT')
        IF (DIRECT) WRITE(LUPRI,'(/1X,A)') 'AO direct transformation'
        CALL FLSHFO(LUPRI)
      END IF
      
      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
        WRITE(LUPRI,'(/1x,a)') 'CC_GMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_GMAT...'
        CALL QUIT('Unknown CC method in CC_GMAT.')
      END IF

      IF (LISTB(1:1).NE.'R' .OR. LISTC(1:1).NE.'R') THEN
        WRITE(LUPRI,'(2A)') 'LISTB and LISTC must refer to',
     &                    ' t-amplitude vectors in CC_GMAT.'
        CALL QUIT('Illegal LISTB or LISTC in CC_GMAT.')
      END IF

      !set CCSDR12
      CCSDR12 = CCR12 .AND. .NOT.(CC2.OR.CCS)

      IF (IPRINT.GT.0) THEN
 
         WRITE (LUPRI,'(//1X,A1,50("="),A1)')'+','+'

         IF (LISTA(1:2).EQ.'L0') THEN
            WRITE (LUPRI,'(1x,A52)')
     &         '|        G MATRIX TRANSFORMATION SECTION           |'
         ELSE
            WRITE (LUPRI,'(1X,A52)')
     &         '|    GENERALIZED G MATRIX TRANSFORMATION SECTION   |'
         END IF

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
        
         IF (IOPTRES.EQ.5) THEN
            WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
            WRITE (LUPRI,'(1X,A52)')
     &         '| L vec. | R vec. | R vec. |  # prod. |            |'
            WRITE (LUPRI,'(1X,4(A,A3),A)') '| ',LISTA(1:3), ' No.| ',
     &         LISTB(1:3),' No.| ',LISTC(1:3),' No.| with ',
     &         LISTD(1:3),' | time/secs  |'
         ELSE 
            WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
            WRITE (LUPRI,'(1X,A52)')
     &         '| L vec. | R vec. | R vec. |  result  |            |'
            WRITE (LUPRI,'(1X,A2,A3,3(A,A3),A)') '| ',LISTA(1:3),
     &         ' No.| ',LISTB(1:3),' No.| ',LISTC(1:3),' No.|  ',
     &         LISTD(1:3),' No. | time/secs  |'
         END IF

         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'

      END IF
      
*---------------------------------------------------------------------*
* initialize timings:
*---------------------------------------------------------------------*
      TIMALL  = SECOND()
      TIMIO   = ZERO
      TIMFCK  = ZERO
      TIMCCS  = ZERO
      TIMINT  = ZERO
      TIMRDAO = ZERO
      TIM21CD = ZERO
      TIM21SD = ZERO
      TIMDBL  = ZERO
      TIMTRBT = ZERO
      TIMTRIP = ZERO

*---------------------------------------------------------------------*
* flush print unit, set dummy work space address to minus a large numb.
*---------------------------------------------------------------------*
      Call FLSHFO(LUPRI)

      KDUM  = -99 999 999  ! dummy address

      IF (LOCDBG) THEN
        WRITE(LUPRI,'(/1x,a,i15)') 'work space in CC_GMAT:',LWORK
        Call FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
* if IOPTRES=0,1 open direct access file for result vectors:
* if IOPTRES=5   initialize GCON array
*---------------------------------------------------------------------*
      IF ( IOPTRES.EQ.0  .OR.  IOPTRES.EQ.1 ) THEN

         FILGMA = LISTD(1:6)

         LUGMAT = -1
         CALL WOPEN2(LUGMAT,FILGMA,64,0)
         IADRDEL = 1

      ELSE IF ( IOPTRES.EQ.3  .OR.  IOPTRES.EQ.4 ) THEN
         CONTINUE
      ELSE IF ( IOPTRES.EQ.5 ) THEN
         IF (MXVEC*NGTRAN.NE.0) CALL DZERO(GCON,MXVEC*NGTRAN)
      ELSE
         CALL QUIT('Illegal value of IOPTRES in CC_GMAT.')
      END IF

*---------------------------------------------------------------------*
* set special symmetry array for one-index transformed integrals
* and several maximum dimensions:
*---------------------------------------------------------------------*
      MDISAO = 0
      MDSRHF = 0
      M2BST  = 0
      MT2BGD = 0
      MT2SQ  = 0
      MT2AM  = 0
      DO  ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO 
        NTAIKJ(ISYM) = ICOUNT
        MDISAO = MAX(MDISAO,NDISAO(ISYM))
        MDSRHF = MAX(MDSRHF,NDSRHF(ISYM))
        M2BST  = MAX(M2BST ,N2BST(ISYM) )
        MT2BGD = MAX(MT2BGD,NT2BGD(ISYM))
        MT2SQ  = MAX(MT2SQ ,NT2SQ(ISYM) )
        MT2AM  = MAX(MT2AM ,NT2AM(ISYM) )
      END DO

*---------------------------------------------------------------------*
* allocate work space for (ia|jb) at the very end of the work space
* and read them in and calculate L(ia|jb) in place:
*---------------------------------------------------------------------*
      ISYOVOV = ISYM0
      KLIAJB  = LWORK  - NT2AM(ISYM0)
      LWRK0   = KLIAJB - 1

      IF (LWRK0 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_GMAT. (0)')
      END IF

      Call CCG_RDIAJB(WORK(KLIAJB),NT2AM(ISYOVOV))
 
      IOPTTCME = 1
      Call CCSD_TCMEPK(WORK(KLIAJB),ONE,ISYOVOV,IOPTTCME)

*---------------------------------------------------------------------*
* allocate work space for zeroth-order amplitudes and Lambda matrices,
* read amplitudes and calculate Lambda matrices:
*---------------------------------------------------------------------*
      KT1AMP0 = 1
      KLAMDP0 = KT1AMP0 + NT1AM(ISYM0)
      KLAMDH0 = KLAMDP0 + NGLMDT(ISYM0)
      KEND1   = KLAMDH0 + NGLMDT(ISYM0)
      LEND1   = LWRK0   - KEND1
      
      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_GMAT. (1a)')
      END IF

* read zeroth-order amplitudes:
      IOPT = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

* get zeroth-order Lambda matrices:
      Call LAMMAT(WORK(KLAMDP0),WORK(KLAMDH0),WORK(KT1AMP0),
     &            WORK(KEND1),LEND1)

*=====================================================================*
* global intermediate & CCS section: 
*     allocate work space for global two-index intermediates
*     and calculate them in a loop over transformations,
*     then calculate the `CCS' contributions to the result vector
*=====================================================================*

      DO ITRAN = 1, NGTRAN
       
         TIMTRN(ITRAN) = - SECOND()
C
C        ------------------------------------------
C        set vector indeces and several symmetries:
C        ------------------------------------------
C
         IDLSTA = IGTRAN(1,ITRAN)
         IDLSTB = IGTRAN(2,ITRAN)
         IDLSTC = IGTRAN(3,ITRAN)

         ISYCTR = ILSTSYM(LISTA,IDLSTA)
         ISYMTB = ILSTSYM(LISTB,IDLSTB)
         ISYMTC = ILSTSYM(LISTC,IDLSTC)

         ISYTCTB = MULD2H(ISYMTB,ISYMTC)   ! TAmpB x TAmpC
         ISYRES  = MULD2H(ISYTCTB,ISYCTR)  ! result vector

         ISYMFB = MULD2H(ISYMTB,ISYOVOV)
         ISYMFC = MULD2H(ISYMTC,ISYOVOV)
         ISYFBC = MULD2H(ISYTCTB,ISYOVOV)
C
         IF (LOCDBG) THEN
           WRITE (LUPRI,*) LISTA,LISTB,LISTC,LISTD
           WRITE (LUPRI,*) IDLSTA,IDLSTB,IDLSTC,ITRAN
           WRITE (LUPRI,*) ISYCTR,ISYMTB,ISYMTC,ISYRES,ISYTCTB
         END IF
C
C        ------------------------------------------------------------
C        allocate single parts of the vectors and Fock intermediates:
C        ------------------------------------------------------------
C
         KDELTA1(ITRAN) = KEND1
         KT1AMPB(ITRAN) = KDELTA1(ITRAN) + NT1AM(ISYRES)
         KLAMDHB(ITRAN) = KT1AMPB(ITRAN) + NT1AM(ISYMTB)
         KLAMDPB(ITRAN) = KLAMDHB(ITRAN) + NGLMDT(ISYMTB)
         KT1AMPC(ITRAN) = KLAMDPB(ITRAN) + NGLMDT(ISYMTB)
         KLAMDHC(ITRAN) = KT1AMPC(ITRAN) + NT1AM(ISYMTC)
         KLAMDPC(ITRAN) = KLAMDHC(ITRAN) + NGLMDT(ISYMTC)
         KCTR1(ITRAN)   = KLAMDPC(ITRAN) + NGLMDT(ISYMTC)
         KFOCKB(ITRAN)  = KCTR1(ITRAN)   + NT1AM(ISYCTR)
         KFOCKC(ITRAN)  = KFOCKB(ITRAN)  + NT1AM(ISYMFB)
         KFBCOO(ITRAN)  = KFOCKC(ITRAN)  + NT1AM(ISYMFC)
         KFBCVV(ITRAN)  = KFBCOO(ITRAN)  + NMATIJ(ISYFBC)
         KEND1          = KFBCVV(ITRAN)  + NMATAB(ISYFBC)

         LEND1 = LWRK0 - KEND1

*        allocate scratch vector needed in Fock matrices calculations
*        and for the first CCS term:
         KSCR   = KEND1
         LSCR   = MAX(NMATIJ(ISYFBC),NMATAB(ISYFBC),NT1AM(ISYRES))

         IF (LEND1 .LT. LSCR) THEN
           CALL QUIT('Insufficient work space in CC_GMAT. (1b)')
         END IF
C
C        ---------------------------------------------
C        initialize single parts of the result vector:
C        ---------------------------------------------
C
*        initialize singles excitation result vector:
         Call DZERO(WORK(KDELTA1(ITRAN)),NT1AM(ISYRES))
 
         DTIME = SECOND()
C
C        --------------------------------------------------------------
C        for CCS and zeroth order zeta vector all contributions vanish:
C        ---> we skip the rest of the loop
C        --------------------------------------------------------------
C
         IF (CCS .AND. LISTA(1:2).EQ.'L0') GOTO 666
C
*        read A response zeta vector:
         IOPT = 1
         Call CC_RDRSP(LISTA,IDLSTA,ISYCTR,IOPT,MODEL,
     &                 WORK(KCTR1(ITRAN)),WORK(KDUM))

*        read B response amplitudes:
         IOPT = 1
         Call CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     &                 WORK(KT1AMPB(ITRAN)),WORK(KDUM))

*        read C response amplitudes:
         IOPT = 1
         Call CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                 WORK(KT1AMPC(ITRAN)),WORK(KDUM))

         TIMIO = TIMIO + SECOND() - DTIME

         IF (LOCDBG) THEN
           Call AROUND('zeroth order T1 amplitudes:')
           Call CC_PRP(WORK(KT1AMP0),WORK(KDUM),ISYM0,1,0)

           Call AROUND('ZETA1 multipliers:')
           WRITE (LUPRI,*) 
     *                'List, Index, Symmetry:',LISTA, IDLSTA, ISYCTR
           Call CC_PRP(WORK(KCTR1(ITRAN)),WORK(KDUM),ISYCTR,1,0)

           Call AROUND('response T1 amplitudes B:')
           WRITE (LUPRI,*) 
     *                'List, Index, Symmetry:',LISTB, IDLSTB, ISYMTB
           Call CC_PRP(WORK(KT1AMPB(ITRAN)),WORK(KDUM),ISYMTB,1,0)

           Call AROUND('response T1 amplitudes C:')
           WRITE (LUPRI,*) 
     *                'List, Index, Symmetry:',LISTC, IDLSTC, ISYMTC
           Call CC_PRP(WORK(KT1AMPC(ITRAN)),WORK(KDUM),ISYMTC,1,0)

           Call FLSHFO(LUPRI)
         END IF
C
C        -----------------------------------
C        calculate response Lambda matrices:
C        -----------------------------------
C
*        calculate B response Lambda matrices:
         Call CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPB(ITRAN)),
     &                    WORK(KLAMDH0),WORK(KLAMDHB(ITRAN)),
     &                    WORK(KT1AMPB(ITRAN)),ISYMTB)

*        calculate C response Lambda matrices:
         Call CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPC(ITRAN)),
     &                    WORK(KLAMDH0),WORK(KLAMDHC(ITRAN)),
     &                    WORK(KT1AMPC(ITRAN)),ISYMTC)
C
C        -------------------------------------------------------------
C        calculate occupied/virtual part of first order Fock matrices:
C        -------------------------------------------------------------
C
         DTIME  = SECOND()

*        calculate tF^B_ia:
         Call CCG_LXD(WORK(KFOCKB(ITRAN)), ISYMFB,
     &                WORK(KT1AMPB(ITRAN)),ISYMTB,
     &                WORK(KLIAJB),        ISYOVOV,0)
      
*        calculate tF^C_ia:
         Call CCG_LXD(WORK(KFOCKC(ITRAN)), ISYMFC, 
     &                WORK(KT1AMPC(ITRAN)),ISYMTC, 
     &                WORK(KLIAJB),        ISYOVOV,0)
      
         IF (LOCDBG) THEN
           Call AROUND(' one-index B transformed Fock marix ')
           Call CC_PRP(WORK(KFOCKB(ITRAN)),WORK(KDUM),ISYMFB,1,0)

           Call AROUND(' one-index C transformed Fock marix ')
           Call CC_PRP(WORK(KFOCKC(ITRAN)),WORK(KDUM),ISYMFC,1,0)
         END IF
C
C        ------------------------------------------------
C        double one-index transformed Fock matrix ttF^BC: 
C        ------------------------------------------------
C
*        calculate occ/occ block of double one-index 
*        transformed Fock matrix:
         Call CCG_1ITROO(WORK(KFBCOO(ITRAN)), ISYFBC,
     &                   WORK(KFOCKB(ITRAN)), ISYMFB,
     &                   WORK(KT1AMPC(ITRAN)),ISYMTC)

         Call CCG_1ITROO(WORK(KSCR),          ISYFBC,
     &                   WORK(KFOCKC(ITRAN)), ISYMFC,
     &                   WORK(KT1AMPB(ITRAN)),ISYMTB)

         Call DAXPY(NMATIJ(ISYFBC),ONE,WORK(KSCR),1,
     &              WORK(KFBCOO(ITRAN)),1)

*        calculate vir/vir block of double one-index 
*        transformed Fock matrix:
         Call CCG_1ITRVV(WORK(KFBCVV(ITRAN)), ISYFBC,
     &                   WORK(KFOCKB(ITRAN)), ISYMFB,
     &                   WORK(KT1AMPC(ITRAN)),ISYMTC  )

         Call CCG_1ITRVV(WORK(KSCR),          ISYFBC,
     &                   WORK(KFOCKC(ITRAN)), ISYMFC,
     &                   WORK(KT1AMPB(ITRAN)),ISYMTB  )

         Call DAXPY(NMATAB(ISYFBC),ONE,WORK(KSCR),1,
     &              WORK(KFBCVV(ITRAN)),1)

         TIMFCK  = TIMFCK + SECOND() - DTIME 
C
C        -------------------------------------------------
C          CCS part:  < Zeta_1 | [ttH^BC, tau_1] | HF>
C        -------------------------------------------------
C
         DTIME = SECOND()

*        do one-index transformation with Zeta vector:
         IOPT  = 2
         Call CCG_1ITRVO(WORK(KSCR),ISYRES,
     &                   WORK(KFBCOO(ITRAN)),WORK(KFBCVV(ITRAN)),
     &                   ISYTCTB,WORK(KCTR1(ITRAN)),ISYCTR,IOPT  )

*        add contribution to result vector:
         Call DAXPY(NT1AM(ISYRES),ONE,WORK(KSCR),1,
     &             WORK(KDELTA1(ITRAN)),1)

         IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG, 'ISYCTR,-RES,-TCTB:',ISYCTR,ISYRES,
     &           ISYTCTB
          WRITE (LUPRI,*) MSGDBG, 'CTR1:'
          Call CC_PRP(WORK(KCTR1(ITRAN)),WORK,ISYCTR,1,0)
          WRITE (LUPRI,*) MSGDBG, 
     &         'double one-index trans. Fock mat. (o/o):'
          WRITE (LUPRI,'(5D21.14)')
     &         (WORK(KFBCOO(ITRAN)-1+I),I=1,NMATIJ(ISYFBC))
          WRITE (LUPRI,*) MSGDBG, 
     &         'double one-index trans. Fock mat. (v/v):'
          WRITE (LUPRI,'(5D21.14)') 
     &         (WORK(KFBCVV(ITRAN)-1+I),I=1,NMATAB(ISYFBC))
          WRITE (LUPRI,*) MSGDBG, 
     &         'contribution of double trans. F to DELTA:'
          Call CC_PRP(WORK(KSCR),WORK,ISYRES,1,0)
         END IF
C
C        -------------------------------------------------------
C        second CCS contribution (involving L(ia,jb) integrals):
C        -------------------------------------------------------
C
         ISYCTB = MULD2H(ISYCTR,ISYMTB)

*        allocate temporay memory for tZeta^B, ttZeta^BC matrices 
*        and a scratch vector: 
         KTTZET = KEND1
         KTZBOO = KTTZET + NT1AM(ISYRES)
         KTZBVV = KTZBOO + NMATIJ(ISYCTB)
         KSCR   = KTZBVV + NMATAB(ISYCTB)
         KEND2  = KSCR   + NT1AM(ISYRES)
         LEND2  = LWORK  - KEND2

         IF (LEND2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_GMAT. (2)')
         END IF

*        calculate one-index transformed Zeta matrix:
         Call CCG_1ITROO(WORK(KTZBOO),         ISYCTB,
     &                   WORK(KCTR1(ITRAN)),   ISYCTR,
     &                   WORK(KT1AMPB(ITRAN)), ISYMTB )

         Call CCG_1ITRVV(WORK(KTZBVV),         ISYCTB,
     &                   WORK(KCTR1(ITRAN)),   ISYCTR,
     &                   WORK(KT1AMPB(ITRAN)), ISYMTB )

* calculate double one-index transformed Zeta matrix:
         IOPT  = 1
         Call CCG_1ITRVO(WORK(KTTZET),ISYRES,
     &                   WORK(KTZBOO),WORK(KTZBVV),ISYCTB,
     &                   WORK(KT1AMPC(ITRAN)),ISYMTC,IOPT )

* contract with L(ia|jb):
         Call CCG_LXD(WORK(KSCR), ISYRES, WORK(KTTZET), 
     &                ISYRES, WORK(KLIAJB), ISYOVOV,0)

         IF (LOCDBG) THEN
           Call AROUND('double one-index transformed Zeta matrix:')
           Call CC_PRP(WORK(KTTZET),WORK(KDUM),ISYRES,1,0)
           Call AROUND('Contribution of L(ia|jb) to CCS part:')
           Call CC_PRP(WORK(KSCR),WORK(KDUM),ISYRES,1,0)
         END IF

* add result to delta vector: 
         Call DAXPY(NT1AM(ISYRES),ONE,WORK(KSCR),1,
     &              WORK(KDELTA1(ITRAN)),1)

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG, 
     &           'Delta vector at the end of the CCS part:'
           Call CC_PRP(WORK(KDELTA1(ITRAN)),WORK(KDUM),ISYRES,1,0)
           WRITE (LUPRI,*) MSGDBG, 'Norm of Delta1 after CCS part:',
     &       DSQRT(DDOT(NT1AM(ISYRES),WORK(KDELTA1(ITRAN)),1,
     &                                WORK(KDELTA1(ITRAN)),1))
         END IF
C
 666     CONTINUE
C
         TIMCCS = TIMCCS + SECOND() - DTIME

         TIMTRN(ITRAN) = TIMTRN(ITRAN) + SECOND()

      END DO
*---------------------------------------------------------------------*
* end of CCS part; regain work space from L(ia|jb)
* in memory now for each transformation:
*
*     KT1AMP0: singles cluster amplitudes 
*     KLAMDH0: zeroth-order Lambda hole matrix 
*     KLAMDP0: zeroth-order Lambda particle matrix 
*     KDELTA1: singles result vector
*     KT1AMPB: singles response amplitudes B
*     KLAMDHB: response Lambda hole matrix for B
*     KLAMDPB: response LAmbda particle matrix for B
*     KT1AMPC: singles response amplitudes C
*     KLAMDHC: response Lambda hole matrix for C
*     KLAMDPC: response LAmbda particle matrix for C
*     KCTR1  : singles response multipliers A
*     KFOCKB : Fock matrix tF^B(i a)  
*     KFOCKC : Fock matrix tF^C(i a)
*     KFBCOO : Fock matrix ttF^BC(i j)
*     KFBCVV : Fock matrix ttF^BC(a b)
*
*     KEND1  : free work space
*
*---------------------------------------------------------------------*

      LEND1 = LWORK - KEND1
      
      IF (CCS) GOTO 777

*=====================================================================*
*
* CC2 part: C and D terms: <Zeta_2| [ttH^BC, tau_1] |HF>
*
* this part needs the AO integrals and the double excitation part
* of the Lagrangian multipliers (Zeta_2) --> we make a compromise
* between memory and CPU requirements:
*
* -- outermost batched loop over transformations with different Zeta_2
*    read all Zeta_2 vectors for this batch
*    -- loop over integral distributions
*       -- loop over transformations within this batch
*          calculate CCG_21CD contributions
*       -- end loop
*    -- end loop
* -- end loop
*
* >> within the loop over AO integrals no I/O over vectors
* >> integral recalculation scales only with nb. of Zeta_2 vectors
* >> only one integral calculation for real G matrix: Zeta = 'L0'
*
*=====================================================================*

*---------------------------------------------------------------------*
* reserve memory for AO integrals and work space ERI & CCG_21CD:
*---------------------------------------------------------------------*
      MSCRATCH = MDISAO + 2*MDSRHF + 10*M2BST + 2*MT2BGD
     &           + MXPRIM*MXCONT   + (8*MXSHEL*MXCONT + 1)/IRAT
      MSCRATCH = MAX(MSCRATCH,MT2AM)

      IF ( LEND1 .LT. MSCRATCH+MT2SQ ) THEN
        CALL QUIT('Insufficient work space in CC_GMAT. (CC2)')
      END IF

*---------------------------------------------------------------------*
* start batched loop over transformations with variable batch lengths:
*---------------------------------------------------------------------*
      IEND = 0

      DO WHILE ( IEND.LT.NGTRAN )
C
C       --------------------------------------------------
C       determine length of next batch and read the double
C       excitation parts of the Lagrangian multipliers:
C       (read as packed matrix, which is then squared up)
C       --------------------------------------------------
C
        ISTRT  = IEND + 1

        IDLSTA = -999
        KEND2  = KEND1
        LEND2  = LEND1
        MFREE  = LEND1 - MSCRATCH

        DO WHILE (MFREE.GT.MT2SQ .AND. IEND.LT.NGTRAN)
          IEND = IEND + 1 
          IF ( IGTRAN(1,IEND) .EQ. IDLSTA ) THEN
             KCTR2(IEND) = KCTR2(IEND-1)
          ELSE
             KCTR2(IEND) = KEND2
             IDLSTA = IGTRAN(1,IEND)
             ISYCTR = ILSTSYM(LISTA,IDLSTA)
             KEND2  = KCTR2(IEND) + NT2SQ(ISYCTR)
             LEND2  = LWORK - KEND2
             MFREE  = LEND2 - MSCRATCH

             IOPT = 2
             DTIME = SECOND()
             Call CC_RDRSP(LISTA,IDLSTA,ISYCTR,IOPT,MODEL,
     &                         WORK(KDUM),WORK(KEND2))
             TIMIO = TIMIO + SECOND() - DTIME

             Call CC_T2SQ (WORK(KEND2), WORK(KCTR2(IEND)), ISYCTR)
             
          END IF
        END DO

*---------------------------------------------------------------------*
* initialize integral calculation and start loop over calls to ERI:
*---------------------------------------------------------------------*
        DTIME = SECOND()
       
        IF (DIRECT) THEN
          NTOSYM = 1

          IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND2),LEND2,IPRERI)
          ELSE
            KCCFB1 = KEND2
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND2  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LEND2  = LWORK  - KEND2
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND2,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND2),LEND2,IPRERI)
            KEND2 = KFREE
            LEND2 = LFREE
          END IF
        ELSE
          NTOSYM = NSYM
        END IF

        TIMINT = TIMINT + SECOND() - DTIME

        KENDSV = KEND2
        LWRKSV = LEND2

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
  
          DTIME = SECOND()

          IF (DIRECT) THEN

            KEND2 = KENDSV
            LEND2 = LWRKSV

            IF (HERDIR) THEN
              CALL HERDI2(WORK(KEND2),LEND2,INDEXA,ILLL,NUMDIS,
     &                    IPRINT)          
            ELSE
              CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0D0,
     *                    WORK(KODCL1),WORK(KODCL2),
     *                    WORK(KODBC1),WORK(KODBC2),
     *                    WORK(KRDBC1),WORK(KRDBC2),
     *                    WORK(KODPP1),WORK(KODPP2),
     *                    WORK(KRDPP1),WORK(KRDPP2),
     *                    WORK(KCCFB1),WORK(KINDXB),
     *                    WORK(KEND2),LEND2,IPRERI)
            END IF

            KRECNR = KEND2
            KEND2  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
            LEND2  = LWORK - KEND2
          
            IF (LEND2 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_GMAT.')
            END IF

          ELSE
            NUMDIS = 1
          END IF

          TIMINT = TIMINT + SECOND() - DTIME

*---------------------------------------------------------------------*
*       loop over number of distributions on the disk
*---------------------------------------------------------------------*
          DO IDEL2 = 1, NUMDIS

            IF(DIRECT) THEN
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

            ISYDIS = MULD2H(ISYDEL,ISYMOP)
C
C           ---------------------------------------------------
C           allocate work space for AO integrals, read them in
C           and transform gamma index with zeroth-order XLAMDH0
C           ---------------------------------------------------
C
            KXINT    = KEND2
            KDSRHF0  = KXINT    + NDISAO(ISYDIS)
            KDSRHFBC = KDSRHF0  + NDSRHF(ISYDIS)
            KEND3    = KDSRHFBC + MDSRHF
            LEND3    = LWORK - KEND3

            IF (LEND3 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_GMAT. (3)')
            END IF

            DTIME = SECOND()
            CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LEND3,
     &                  WORK(KRECNR),DIRECT)
            TIMRDAO = TIMRDAO + SECOND() - DTIME
           
            DTIME = SECOND()
            CALL CCTRBT(WORK(KXINT),WORK(KDSRHF0),WORK(KLAMDH0),ISYM0,
     &                  WORK(KEND3),LEND3,ISYDIS)
            TIMTRBT = TIMTRBT + SECOND() - DTIME

C
C           ---------------------------------------
C           loop over transformation in this batch:
C           ---------------------------------------
C
            DO ITRAN = ISTRT, IEND

               TIMTRN(ITRAN) = TIMTRN(ITRAN) - SECOND()

               IDLSTA = IGTRAN(1,ITRAN)
               IDLSTB = IGTRAN(2,ITRAN)
               IDLSTC = IGTRAN(3,ITRAN)

               ISYCTR = ILSTSYM(LISTA,IDLSTA)
               ISYMTB = ILSTSYM(LISTB,IDLSTB)
               ISYMTC = ILSTSYM(LISTC,IDLSTC)

               ISYMXB = MULD2H(ISYDIS,ISYMTB)
               ISYMXC = MULD2H(ISYDIS,ISYMTC)

               ISYTCTB = MULD2H(ISYMTB,ISYMTC)   
               ISYRES  = MULD2H(ISYTCTB,ISYCTR)  
 
               KXLMDHB = KLAMDHB(ITRAN)
               KXLMDHC = KLAMDHC(ITRAN)
               KXLMDPB = KLAMDPB(ITRAN)
               KXLMDPC = KLAMDPC(ITRAN)

               LSAME  =  (LISTC.EQ.LISTB  .AND. IDLSTC.EQ.IDLSTB) 

C              -----------------------------------------------------
C              calculate the first 21CD term (symmetric in B,C):
C              -----------------------------------------------------
               LLSAME = LSAME
               FACT   = ONE
               DTIME  = SECOND()
               Call CCG_21CD( WORK(KDELTA1(ITRAN)),ISYRES, 
     &                        WORK(KCTR2(ITRAN)),  ISYCTR, 
     &                        WORK(KLAMDH0), WORK(KLAMDP0),
     &                        WORK(KXLMDHB), WORK(KXLMDPB), ISYMTB,
     &                        WORK(KXLMDHC), WORK(KXLMDPC), ISYMTC,
     &                        WORK(KDSRHF0), ISYDIS,  IDEL, ISYDEL, 
     &                        LLSAME, FACT,  WORK(KEND3), LEND3 )
               TIM21CD = TIM21CD + SECOND() - DTIME
 

C              ------------------------------------------------------
C              transform gamma to occupied, using the XLAMDHB matrix 
C              and calculate the second 21CD term (symmetric in C,0):
C              ------------------------------------------------------
               DTIME = SECOND()
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHFBC),WORK(KXLMDHB),
     &                     ISYMTB,WORK(KEND3),LEND3,ISYDIS)
               TIMTRBT = TIMTRBT + SECOND() - DTIME
          
               FACT   = ONE
               IF (LSAME) FACT = TWO
               LLSAME = .FALSE.
               DTIME  = SECOND()
               Call CCG_21CD( WORK(KDELTA1(ITRAN)),ISYRES, 
     &                        WORK(KCTR2(ITRAN)),  ISYCTR, 
     &                        WORK(KLAMDH0), WORK(KLAMDP0),
     &                        WORK(KXLMDHC), WORK(KXLMDPC), ISYMTC,
     &                        WORK(KLAMDH0), WORK(KLAMDP0), ISYM0,
     &                        WORK(KDSRHFBC),ISYMXB, IDEL,  ISYDEL, 
     &                        LLSAME, FACT,  WORK(KEND3), LEND3 )
               TIM21CD = TIM21CD + SECOND() - DTIME


C              -----------------------------------------------------
C              transform gamma to occupied, using the XLAMDHC matrix
C              and calculate the third 21CD term (symmetric in B,0):
C              -----------------------------------------------------
               IF (.NOT. LSAME) THEN
                  DTIME  = SECOND()
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHFBC),WORK(KXLMDHC),
     &                        ISYMTC,WORK(KEND3),LEND3,ISYDIS)
                  TIMTRBT = TIMTRBT + SECOND() - DTIME

                  FACT   = ONE
                  LLSAME = .FALSE.
                  DTIME  = SECOND()
                  Call CCG_21CD( WORK(KDELTA1(ITRAN)),ISYRES, 
     &                           WORK(KCTR2(ITRAN)),  ISYCTR, 
     &                           WORK(KLAMDH0), WORK(KLAMDP0),
     &                           WORK(KLAMDH0), WORK(KLAMDP0), ISYM0,
     &                           WORK(KXLMDHB), WORK(KXLMDPB), ISYMTB,
     &                           WORK(KDSRHFBC),ISYMXC, IDEL,  ISYDEL, 
     &                           LLSAME, FACT,  WORK(KEND3), LEND3 )
                  TIM21CD = TIM21CD + SECOND() - DTIME
               END IF

               TIMTRN(ITRAN) = TIMTRN(ITRAN) + SECOND()

            END DO ! loop over transformations for this batch

          END DO ! IDEL2 = 1, NUMDIS
          END DO ! ILLL = 1, NTOT
        END DO ! ISYMD1 = 1, NTOSYM
      END DO ! loop over batches of transformations

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,
     &        'Delta vectors at the end of the CC2 part:'
        DO ITRAN = 1, NGTRAN
          IDLSTA  = IGTRAN(1,ITRAN)
          IDLSTB  = IGTRAN(2,ITRAN)
          IDLSTC  = IGTRAN(3,ITRAN)
          ISYTCTB = MULD2H(ILSTSYM(LISTB,IDLSTB),ILSTSYM(LISTC,IDLSTC))
          ISYRES  = MULD2H(ISYTCTB,ILSTSYM(LISTA,IDLSTA))  
          WRITE (LUPRI,*) MSGDBG, 'ITRAN  = ',ITRAN
          WRITE (LUPRI,*) MSGDBG, 'IDLSTA,IDLSTB,IDLSTC=',IDLSTA,
     &         IDLSTB,IDLSTC
          Call CC_PRP(WORK(KDELTA1(ITRAN)),WORK(KDUM),ISYRES,1,0)
          WRITE (LUPRI,*) MSGDBG, 'Norm of Delta1 after CC2 part:',
     &           DSQRT(DDOT(NT1AM(ISYRES),WORK(KDELTA1(ITRAN)),1,
     &                                    WORK(KDELTA1(ITRAN)),1))
        END DO
      END IF
*---------------------------------------------------------------------*
* end of CC2 part; recover memory form CTR2
* in memory now for each transformation:
*
*    KT1AMP0: singles cluster amplitudes 
*    KLAMDH0: zeroth-order Lambda hole matrix 
*    KLAMDP0: zeroth-order LAmbda particle matrix 
*    KDELTA1: singles excitation part of the result vector
*    KT1AMPB: single response amplitudes B
*    KLAMDHB: response Lambda hole matrix for B
*    KLAMDPB: response LAmbda particle matrix for B
*    KT1AMPC: single response amplitudes C
*    KLAMDHC: response Lambda hole matrix for C
*    KLAMDPC: response LAmbda particle matrix for C
*    KCTR1  : single response multipliers A
*    KFOCKB : Fock matrix tF^B(i a) 
*    KFOCKC : Fock matrix tF^C(i a)
*    KFBCOO : Fock matrix ttF^BC(i j)
*    KFBCVV : Fock matrix ttF^BC(a b)
*
*    KEND1  : free work space
*
*---------------------------------------------------------------------*

C     ---------------------------
C     here we continue for CCS...
C     ---------------------------
 777  CONTINUE

      KEND0 = KEND1

*=====================================================================*
* the remaining parts are more memory intensive, but to not require
* the AO integrals --> we complete the transformations one by one
* and save the result one file, or calculate the required dot products
*=====================================================================*

      DO ITRAN = 1, NGTRAN

         TIMTRN(ITRAN) = TIMTRN(ITRAN) - SECOND()

         IDLSTA = IGTRAN(1,ITRAN)
         IDLSTB = IGTRAN(2,ITRAN)
         IDLSTC = IGTRAN(3,ITRAN)

         ISYCTR = ILSTSYM(LISTA,IDLSTA)
         ISYMTB = ILSTSYM(LISTB,IDLSTB)
         ISYMTC = ILSTSYM(LISTC,IDLSTC)

         ISYCTB = MULD2H(ISYCTR,ISYMTB)
         ISYCTC = MULD2H(ISYCTR,ISYMTC)

         ISYTCTB = MULD2H(ISYMTB,ISYMTC)   
         ISYRES  = MULD2H(ISYTCTB,ISYCTR)  

         ISYMFB = MULD2H(ISYMTB,ISYOVOV)
         ISYMFC = MULD2H(ISYMTC,ISYOVOV)
         ISYFBC = MULD2H(ISYTCTB,ISYOVOV)
C
         IF (CCS. OR. CC2) GOTO 888
C
C        -------------------------------------------------------
C        allocate work space for one-index transformed integrals
C        and Lagrangian multipliers:
C        -------------------------------------------------------
C
         KXBIAJK = KEND0
         KEBIAJK = KXBIAJK + NTAIKJ(ISYMTB)
         KZBIAJK = KEBIAJK + NTAIKJ(ISYMTB)
         KCBIAJK = KZBIAJK + NTAIKJ(ISYCTB)
         KEND1   = KCBIAJK + NTAIKJ(ISYCTB)

         KXCIAJK = KEND1
         KECIAJK = KXCIAJK + NTAIKJ(ISYMTC)
         KZCIAJK = KECIAJK + NTAIKJ(ISYMTC)
         KCCIAJK = KZCIAJK + NTAIKJ(ISYCTC)
         KEND1   = KCCIAJK + NTAIKJ(ISYCTC)

         LEND1   = LWORK   - KEND1

*=====================================================================*
* CCSD part for DELTA1
*=====================================================================*
 
         DTIME = SECOND()
C
C        ------------------------------------------
C        DELTBC = <Zeta_2| [[tH^B,T2C], tau_1] |HF>
C        ------------------------------------------
C
*        read double excitation part of response vector TAmpC:
         KT2AMPC = KEND1
         KDELTBC = KT2AMPC + NT2AM(ISYMTC)
         KEND2   = KDELTBC + NT1AM(ISYRES)
         LEND2   = LWORK - KEND2

         IF (LEND2 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_GMAT.')
         END IF

         IOPT = 2
         Call CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     &                 WORK(KDUM),WORK(KT2AMPC))

         Call CCLR_DIASCL(WORK(KT2AMPC),TWO,ISYMTC)

         Call DZERO(WORK(KDELTBC),NT1AM(ISYRES))

         Call CCG_21CCSD(WORK(KDELTBC),ISYRES,LISTA,IDLSTA, 
     &                   WORK(KFOCKB(ITRAN)), ISYMFB,
     &                   WORK(KT1AMPB(ITRAN)),ISYMTB, 
     &                   WORK(KT2AMPC), ISYMTC,
     &                   LISTC, IDLSTC, 
     &                   WORK(KXBIAJK), WORK(KEBIAJK),        
     &                   WORK(KZBIAJK), WORK(KCBIAJK),        
     &                   WORK(KEND2),   LEND2                 )

         LSAME  =  (LISTC.EQ.LISTB  .AND. IDLSTC.EQ.IDLSTB) 

         IF (.NOT. LSAME) THEN
            CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KDELTBC),1,
     &                                   WORK(KDELTA1(ITRAN)),1)
         ELSE
            CALL DAXPY(NT1AM(ISYRES),TWO,WORK(KDELTBC),1,
     &                                   WORK(KDELTA1(ITRAN)),1)

            CALL DCOPY(NTAIKJ(ISYMTB),WORK(KXBIAJK),1,WORK(KXCIAJK),1)
            CALL DCOPY(NTAIKJ(ISYMTB),WORK(KEBIAJK),1,WORK(KECIAJK),1)
            CALL DCOPY(NTAIKJ(ISYCTB),WORK(KZBIAJK),1,WORK(KZCIAJK),1)
            CALL DCOPY(NTAIKJ(ISYCTB),WORK(KCBIAJK),1,WORK(KCCIAJK),1)
         END IF


         IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG, '1. CCG_21CCSD contrib. to Delta1:'
           WRITE (LUPRI,*) MSGDBG, 'Norm of DELTBC:',
     &       DSQRT(DDOT(NT1AM(ISYRES),WORK(KDELTBC),1,WORK(KDELTBC),1))
           XNORM = ZERO
           IF (ISYMTC.EQ.ISYRES) XNORM  = 
     &       DDOT(NT1AM(ISYRES),WORK(KDELTBC),1,WORK(KT1AMPC(ITRAN)),1)
           WRITE(LUPRI,*) MSGDBG, 'DELTBC dotted with C amplitudes:',
     &          XNORM
           XNORM = ZERO
           IF (ISYMTB.EQ.ISYRES) XNORM  = 
     &       DDOT(NT1AM(ISYRES),WORK(KDELTBC),1,WORK(KT1AMPB(ITRAN)),1)
           WRITE (LUPRI,*) MSGDBG, 'DELTBC dotted with B amplitudes:',
     &          XNORM
           WRITE (LUPRI,*) MSGDBG, 'DELTBC:'
           CALL CC_PRP(WORK(KDELTBC),WORK(KDUM),ISYRES,1,0)
           CALL FLSHFO(LUPRI)
         END IF
C
C        ------------------------------------------
C        DELTCB: <Zeta_2| [[tH^C,T2B], tau_1] |HF>
C        ------------------------------------------
C
         IF (.NOT. LSAME) THEN

*            read double excitation part of response vector TAmpB:
             KT2AMPB = KEND1
             KDELTCB = KT2AMPB + NT2AM(ISYMTB)
             KEND2   = KDELTCB + NT1AM(ISYRES)
             LEND2   = LWORK - KEND2

             IF (LEND2 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_GMAT.')
             END IF

             IOPT = 2
             Call CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     &                     WORK(KDUM),WORK(KT2AMPB))

             Call CCLR_DIASCL(WORK(KT2AMPB),TWO,ISYMTB)

             Call DZERO(WORK(KDELTCB),NT1AM(ISYRES))

             Call CCG_21CCSD (WORK(KDELTCB),ISYRES,LISTA,IDLSTA, 
     &                        WORK(KFOCKC(ITRAN)),  ISYMFC, 
     &                        WORK(KT1AMPC(ITRAN)), ISYMTC, 
     &                        WORK(KT2AMPB), ISYMTB,
     &                        LISTB, IDLSTB, 
     &                        WORK(KXCIAJK), WORK(KECIAJK),         
     &                        WORK(KZCIAJK), WORK(KCCIAJK),         
     &                        WORK(KEND2),   LEND2  )

             Call DAXPY(NT1AM(ISYRES),ONE,WORK(KDELTCB),1,
     &                                    WORK(KDELTA1(ITRAN)),1)

         END IF

         IF (LOCDBG .AND. .NOT. LSAME) THEN
           WRITE (LUPRI,*) MSGDBG, '2. CCG_21CCSD contrib. to Delta1:'
           WRITE (LUPRI,*) MSGDBG, 'Norm of DELTCB:',
     &       DSQRT(DDOT(NT1AM(ISYRES),WORK(KDELTCB),1,WORK(KDELTCB),1))
           WRITE (LUPRI,*) MSGDBG, 'DELTCB dotted with C amplitudes:',
     &       DDOT(NT1AM(ISYRES),WORK(KDELTCB),1,WORK(KT1AMPC(ITRAN)),1)
           WRITE (LUPRI,*) MSGDBG, 'DELTCB dotted with B amplitudes:',
     &       DDOT(NT1AM(ISYRES),WORK(KDELTCB),1,WORK(KT1AMPB(ITRAN)),1)
           WRITE (LUPRI,*) MSGDBG, 'DELTCB:'
           Call CC_PRP(WORK(KDELTCB),WORK(KDUM),ISYRES,1,0)
           Call FLSHFO(LUPRI)
         END IF

         TIM21SD = TIM21SD + SECOND() - DTIME

         IF (LDOUBL) THEN

*=====================================================================*
* CCSD part for DELTA2
*=====================================================================*

         DTIME = SECOND()
C
C        --------------------------------------------------------
C        allocate work space for double excitation result vector:
C        --------------------------------------------------------
C
         KDELTA2 = KEND1
         KEND1   = KDELTA2 + NT2AM(ISYRES)
         LEND1   = LWORK   - KEND1

         IF (LEND1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_GMAT. '//
     &           '(doubles part)')
         END IF

*        initialize DELTA2 array:
         Call DZERO(WORK(KDELTA2), NT2AM(ISYRES))
C
C        -------------------------------------
C        B term: requires packed multipliers 
C                and squared (ia|jb) integrals
C        -------------------------------------
C
           ISYTCTB = MULD2H(ISYMTC,ISYMTB)
           ISYC4O  = MULD2H(ISYTCTB,ISYCTR)

*          memory allocation: 
*          1) squared integrals & scratch space for packed integrals
           KXIAJB = KEND1
           KSCR   = KXIAJB + NT2SQ(ISYOVOV)
           KENDB2 = KSCR   + NT2AM(ISYOVOV)
           LENDB2 = LWORK  - KENDB2

*          2) packed multipliers & double one-index transf. multipl.:
           KZETA2 = KSCR
           KC4O   = KZETA2 + NT2AM(ISYCTR)
           KENDB3 = KC4O  + NGAMMA(ISYC4O)
           LENDB3 = LWORK - KENDB3

           IF ( LENDB2 .LT. 0 .OR. LENDB3 .LT. 0 ) THEN
             CALL QUIT('Insufficient work space in CC_GMAT. (B2/B3)')
           END IF

*          read packed (ia|jb) integrals and square them up 
           Call CCG_RDIAJB(WORK(KSCR),NT2AM(ISYOVOV))

           Call CC_T2SQ (WORK(KSCR), WORK(KXIAJB), ISYOVOV)

*          read packed multipliers:
           KDUM = - 99 999 999  ! dummy address
           IOPT = 2
           Call CC_RDRSP(LISTA,IDLSTA,ISYCTR,IOPT,MODEL,
     &                   WORK(KDUM),WORK(KZETA2))

*          calculate double one-index transformed lagr. multipliers:
           IOPT = 1
           Call CCG_4O(WORK(KC4O),ISYC4O, WORK(KZETA2),ISYCTR,
     &                 WORK(KT1AMPB(ITRAN)),ISYMTB,
     &                 WORK(KT1AMPC(ITRAN)),ISYMTC,
     &                 WORK(KENDB3),LENDB3,IOPT)

*          contract them with the (ia|jb) integrals:
           IOPT = 2
           Call CCRHS_A(WORK(KDELTA2),WORK(KXIAJB),WORK(KC4O), 
     &                  WORK(KENDB3), LENDB3, ISYC4O, ISYOVOV, IOPT)

C
C        -----------------------------------------
C        A term: requires packed (ia|jb) integrals 
C                and squared multipliers
C        -----------------------------------------
C
*        set symmetry for double transformed (ik|jl) integrals
         ISYTCTB = MULD2H(ISYMTC,ISYMTB)
         ISYX4O  = MULD2H(ISYTCTB,ISYOVOV)

*        memory allocation: 
*        1) squared multipliers & scratch space for packed multipliers
         KZETA2 = KEND1
         KSCR   = KZETA2 + NT2SQ(ISYCTR)
         KENDA2 = KSCR  + NT2AM(ISYCTR)
         LENDA2 = LWORK - KENDA2

*        2) packed (ia|jb) integrals & (ik|jl) integrals
         KXIAJB = KSCR
         KX4O   = KXIAJB + NT2AM(ISYOVOV)
         KENDA3 = KX4O   + NGAMMA(ISYX4O)
         LENDA3 = LWORK  - KENDA3

         IF ( LENDA2 .LT. 0 .OR. LENDA3 .LT. 0 ) THEN
           CALL QUIT('Insufficient work space in CC_GMAT. (A2/A3)')
         END IF

*        read packed lagrangian multipliers and square them up 
         IOPT = 2
         Call CC_RDRSP(LISTA,IDLSTA,ISYCTR,IOPT,MODEL,
     &                     WORK(KDUM),WORK(KSCR))

         Call CC_T2SQ (WORK(KSCR), WORK(KZETA2), ISYCTR)

*        read (ia|jb) integrals:
         Call CCG_RDIAJB(WORK(KXIAJB),NT2AM(ISYOVOV))


*        calculate double one-index transformed (ik|jl) integrals:
         IOPT = 1
         Call CCG_4O(WORK(KX4O),ISYX4O, WORK(KXIAJB),ISYOVOV,
     &               WORK(KT1AMPB(ITRAN)),ISYMTB,
     &               WORK(KT1AMPC(ITRAN)),ISYMTC,
     &               WORK(KENDA3),LENDA3,IOPT)


*        contract them with the multipliers:
                IOPT = 2
                Call CCRHS_A(WORK(KDELTA2),WORK(KZETA2),WORK(KX4O),
     &                       WORK(KENDA3),LENDA3,ISYX4O,ISYCTR,IOPT)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG, 'finished A term.'
            Call CC_PRP(WORK(KDELTA1(ITRAN)),WORK(KDELTA2),ISYRES,1,1)
            Call FLSHFO(LUPRI)
         END IF
C
C        -------------------------------------
C        E term: requires squared multipliers
C                (still in memory form A term)
C        -------------------------------------
C
         KZETA2 = KEND1
         KEMAT1 = KZETA2 + NT2SQ(ISYCTR)
         KEMAT2 = KEMAT1 + NMATAB(ISYFBC)
         KENDE3 = KEMAT2 + NMATIJ(ISYFBC)
         LENDE3 = LWORK - KENDE3
 
         IF ( LENDE3 .LT. 0 ) THEN
           CALL QUIT('Insufficient work space in CC_GMAT. (E3)')
         END IF

*        transpose ttF^BC(a b) --> EMAT1(b a)
         DO ISYMA = 1, NSYM
           ISYMB = MULD2H(ISYMA,ISYFBC)
           DO A = 1, NVIR(ISYMA)
           DO B = 1, NVIR(ISYMB)
             NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B-1) + A
             NBA = IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(A-1) + B
         
             WORK(KEMAT1 - 1 + NBA) = WORK(KFBCVV(ITRAN) - 1 + NAB)
           END DO
           END DO
         END DO


*        transpose ttF^BC(i j) --> EMAT2(j i)
         DO ISYMI = 1, NSYM
           ISYMJ = MULD2H(ISYMI,ISYFBC)
           DO I = 1, NRHF(ISYMI)
           DO J = 1, NRHF(ISYMJ)
             NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + I
             NJI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I-1) + J
         
             WORK(KEMAT2 - 1 + NJI) = WORK(KFBCOO(ITRAN) - 1 + NIJ)
           END DO
           END DO
         END DO

*        combine EMAT1/EMAT2 with lagrangian multipliers:
         Call CCRHS_E(WORK(KDELTA2),WORK(KZETA2),WORK(KEMAT1), 
     &                WORK(KEMAT2), WORK(KENDE3),LENDE3,ISYCTR,ISYFBC)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG, 'finished E term.'
            Call CC_PRP(WORK(KDELTA1(ITRAN)),WORK(KDELTA2),ISYRES,1,1)
            Call FLSHFO(LUPRI)
         END IF

*---------------------------------------------------------------------*
* D term & C term:
*---------------------------------------------------------------------*

* recover work space from integrals and Lagrangian multipliers,
* and allocate work space for C & D term as square matrix:
            KCDTERM = KZETA2
            KEND2   = KCDTERM +  NT2SQ(ISYRES)
            LEND2   = LWORK - KEND2
          
            IF ( LEND2 .LT. 0 ) THEN
              CALL QUIT('Insufficient work space in CC_GMAT. (C & D)')
            END IF

            Call CCG_22CD(WORK(KDELTA2),WORK(KCDTERM),ISYRES,
     &                    WORK(KXBIAJK),WORK(KEBIAJK),ISYMTB,         
     &                    WORK(KZBIAJK),WORK(KCBIAJK),ISYCTB,         
     &                    WORK(KXCIAJK),WORK(KECIAJK),ISYMTC,         
     &                    WORK(KZCIAJK),WORK(KCCIAJK),ISYCTC,         
     &                    LSAME, WORK(KEND2), LEND2)


         TIMDBL = TIMDBL + SECOND() - DTIME

         END IF ! (LDOUBL)

*=====================================================================*
* add CC3 corrections:
*=====================================================================*
         DTIME = SECOND()

         IF (CCSDT) THEN
            IF (NODDY_GMAT) THEN
               CALL CCSDT_GMAT_NODDY(LISTA,IDLSTA,LISTB,IDLSTB,
     &                               LISTC,IDLSTC,
     &                               WORK(KDELTA1(ITRAN)),WORK(KDELTA2),
     &                               WORK(KEND1),LEND1)
            ELSE

             LUCKJD   = -1
             LUTOC    = -1
             LU3VI    = -1
             LUDKBC3  = -1
             LU3FOPX  = -1
             LU3FOP2X = -1
C
             CALL WOPEN2(LUCKJD,FNCKJD,64,0)
             CALL WOPEN2(LUTOC,FNTOC,64,0)
             CALL WOPEN2(LU3VI,FN3VI,64,0)
             CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
             CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
             CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)

             CALL CC3_GMATNEW(LISTA,IDLSTA,LISTB,IDLSTB,
     *                        LISTC,IDLSTC,
     *                        WORK(KDELTA1(ITRAN)),WORK(KDELTA2),
     *                        ISYRES,
     *                        WORK(KEND1),LEND1,
     *                        LUTOC,FNTOC,
     *                        LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                        LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)

             CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
             CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
             CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
             CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
             CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
             CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')


            END IF
         END IF

         TIMTRIP = TIMTRIP + SECOND() - DTIME

C
C       CCSLV part, JK+OC
C       The result vectors are ind WORK(KDELTA1(ITRAN)) and WORK(KDELTA2)
C       Futhermore, we need to know the list numbers, their ID's and symmetries
C       For the zeta (A) vector these are: LISTA,IDLSTA,ISYCTR
C       For the B trial vector these are : LISTB,IDLSTB,ISYMTB
C       For the C trial vector these are : LISTC,IDLSTC,ISYMTC
C       Also, we have as arguments the symmetry of B x C : ISYTCTB
C       and the symmetry of A x B X C : ISYRES
C       work arrays used are WORK(KEND2) and LEND2
C
        IF (CCSLV) THEN
          IF (.NOT. CCMM) CALL CCSL_GTR(WORK(KDELTA1(ITRAN)),
     *                    WORK(KDELTA2),ISYRES,
     *                    LISTA,IDLSTA,ISYCTR,
     *                    LISTB,IDLSTB,ISYMTB,
     *                    LISTC,IDLSTC,ISYMTC,
     *                    MODEL,WORK(KEND2),LEND2)
C
          IF (CCMM) THEN
            IF (.NOT. NYQMMM) THEN
              CALL CCMM_GTR(WORK(KDELTA1(ITRAN)),
     *                   WORK(KDELTA2),ISYRES,
     *                   LISTA,IDLSTA,ISYCTR,
     *                   LISTB,IDLSTB,ISYMTB,
     *                   LISTC,IDLSTC,ISYMTC,
     *                   MODEL,WORK(KEND2),LEND2)
C 
            ELSE IF (NYQMMM) THEN
              RSPTYP='G'
              CALL CCMM_QRTRANSFORMER(WORK(KDELTA1(ITRAN)),
     *                   WORK(KDELTA2),ISYRES,
     *                   LISTB,IDLSTB,ISYMTB,
     *                   LISTC,IDLSTC,ISYMTC,
     *                   MODEL,RSPTYP,WORK(KEND2),LEND2)
            END IF
          END IF
        END IF

*=====================================================================*
* CCSD part for CC-R12
*=====================================================================*

         IF (CCSDR12) THEN
C
C          ------------------------------------------------------------
C          allocate work space for r12 double excitation result vector
C          right after conv. doubles and reset KEND1:
C          ------------------------------------------------------------
C
           KDELTAR12 = KEND1
           KEND1     = KDELTAR12 + NTR12AM(ISYRES)
           LEND1     = LWORK   - KEND1

           IF (LEND1 .LT. 0) THEN
             CALL QUIT('Insufficient work space in CC_GMAT. '//
     &                 '(r12 doubles part)')
           END IF

*          initialize DELTAR12 array:
           Call DZERO(WORK(KDELTAR12), NTR12AM(ISYRES))
C
C          add B' term
C
           !read V_(alpha beta)^(kl) from disk
           KVABKL = KEND1
           KEND2  = KVABKL + NVABKL(1)
           LEND2  = LWORK - KEND2
           IF (LEND2 .LT. 0) THEN
             CALL QUIT('Insuff. work space for V^(albe)_kl in CC_GMAT')
           END IF
           lunit = -1
           call gpopen(lunit,fvabkl,'unknown',' ','unformatted',
     &                 idum,.false.)
           read (lunit)(work(kvabkl+i-1),i=1,nvabkl(1))
           call gpclose(lunit,'KEEP')
C
           IOPT = 2
           CALL CC_R12MKVIRT(WORK(KVABKL),WORK(KLAMDPB(ITRAN)),ISYMTB,
     &                       WORK(KLAMDPC(ITRAN)),ISYMTC,
     &                       'R12VCBDBKL',IOPT,WORK(KEND2),LEND2)
C
           !read doubles Lagrangian multipliers
           KZETA2 = KEND1
           KR12SCR= KZETA2 + NT2AM(ISYCTR)
           KEND2  = KR12SCR+ NTR12SQ(ISYRES)
           LEND2  = LWORK - KEND2
           IF (LEND2 .LT. 0) THEN
             CALL QUIT('Insuff. work space in CC_GMAT')
           END IF
C
           CALL DZERO(WORK(KR12SCR),NTR12SQ(ISYRES))
C
           IOPT = 2
           CALL CC_RDRSP(LISTA,IDLSTA,ISYCTR,IOPT,MODEL,
     &                   WORK(KDUM),WORK(KZETA2))
C
           !calculate contribution
           CALL CCRHS_BPP(WORK(KR12SCR),WORK(KZETA2),ISYCTR,.TRUE.,
     &                 'R12VCBDBKL',ISYTCTB,WORK(KEND2),LEND2)
C
C          resort result into a symmetry packed triangular matrix
C
           IOPT = 1
           CALL CCR12PCK2(WORK(KDELTAR12),ISYRES,.FALSE.,
     &                    WORK(KR12SCR),'T',IOPT)
           CALL CCLR_DIASCLR12(WORK(KDELTAR12),0.5D0*KETSCL,ISYRES)
C
         END IF !(CCSDR12)

*=====================================================================*
* end CCSD part for DELTA2 / DELTAR12; save result vector on file or
* calculate the required dot products:
*=====================================================================*
C
C        ---------------------------------
C        here we continue for CCS and CC2:
C        ---------------------------------
 888     CONTINUE
         IF (CC2) THEN
            KDELTA2 = KEND0
            KEND1   = KDELTA2 + NT2AM(ISYRES)
            LEND1   = LWORK   - KEND1

            IF (LEND1 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_GMAT. '//
     &                  '(save section)')
            END IF

            CALL DZERO(WORK(KDELTA2),NT2AM(ISYRES))
C
            IF (CCSLV) THEN
              IF (.NOT. CCMM) CALL CCSL_GTR(WORK(KDELTA1(ITRAN)),
     *                        WORK(KDELTA2),ISYRES,
     *                        LISTA,IDLSTA,ISYCTR,
     *                        LISTB,IDLSTB,ISYMTB,
     *                        LISTC,IDLSTC,ISYMTC,
     *                        MODEL,WORK(KEND1),LEND1)
C
              IF (CCMM) THEN 
                 IF (.NOT.NYQMMM)THEN
                    CALL CCMM_GTR(WORK(KDELTA1(ITRAN)),
     *                            WORK(KDELTA2),ISYRES,
     *                            LISTA,IDLSTA,ISYCTR,
     *                            LISTB,IDLSTB,ISYMTB,
     *                            LISTC,IDLSTC,ISYMTC,
     *                            MODEL,WORK(KEND1),LEND1)
                 ELSE IF (NYQMMM) THEN
                    RSPTYP = 'G'
                    CALL CCMM_QRTRANSFORMER(WORK(KDELTA1(ITRAN)),
     *                            WORK(KDELTA2),ISYRES,
     *                            LISTB,IDLSTB,ISYMTB,
     *                            LISTC,IDLSTC,ISYMTC,
     *                            MODEL,RSPTYP,WORK(KEND1),LEND1)
                 END IF
              END IF
            END IF
         END IF
C
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 
     &           'Final result for G matrix transformation ',ITRAN
            Call CC_PRP(WORK(KDELTA1(ITRAN)),WORK(KDELTA2),ISYRES,1,1)
            IF (CCSDR12) CALL CC_PRPR12(WORK(KDELTAR12),ISYRES,1,.TRUE.)
         END IF

         IFILE  = IGTRAN(4,ITRAN)
         IDLSTD = IGTRAN(4,ITRAN)
         LEN1   = NT1AM(ISYRES)
         LEN2   = NT2AM(ISYRES)
         LENR12 = NTR12AM(ISYRES)
         
         IF (CCS) THEN
           IOPTW  = 1
           MODELW = 'CCS       '
         ELSE IF (CC2) THEN
           IOPTW  = 3
           MODELW = 'CC2       '
         ELSE IF (CCSD) THEN
           IOPTW  = 3
           MODELW = 'CCSD      '
         ELSE IF (CC3) THEN
           IOPTW  = 3
           IOPTWE = 24
           MODELW = 'CC3       '
         ELSE
           CALL QUIT('Unknown coupled cluster model in CC_GMAT.')
         END IF
         IF (CCR12) THEN
           APROXR12 = '   '
           CALL CCSD_MODEL(MODELW,LENMOD,10,MODELW,10,APROXR12)
           IOPTWR12 = 32
         END IF

         IF ( IOPTRES.EQ.0 .OR. IOPTRES.EQ.1 ) THEN
C  
C          ----------------------------------
C          save result on direct access file:
C          ----------------------------------
C
           IGTRAN(4,ITRAN) = IADRDEL

           DTIME = SECOND()

           CALL PUTWA2(LUGMAT,FILGMA,WORK(KDELTA1(ITRAN)),IADRDEL,LEN1)
           IADRDEL = IADRDEL + LEN1

           IF (.NOT.CCS) THEN
              CALL PUTWA2(LUGMAT,FILGMA,WORK(KDELTA2),IADRDEL,LEN2)
              IADRDEL = IADRDEL + LEN2
           END IF

           IF (CCSDR12) THEN
              CALL PUTWA2(LUGMAT,FILGMA,WORK(KDELTAR12),IADRDEL,LENR12)
              IADRDEL = IADRDEL + LENR12
           END IF

           TIMIO = TIMIO + SECOND() - DTIME

         ELSE IF (IOPTRES.EQ.3) THEN
C  
C          -----------------------------
C          write to file using CC_WRRSP:
C          -----------------------------
C
           DTIME = SECOND()
           CALL CC_WRRSP(LISTD,IFILE,ISYRES,IOPTW,MODELW,WORK(KDUM),
     &           WORK(KDELTA1(ITRAN)),WORK(KDELTA2),WORK(KEND1),LEND1)
           IF (CCSDR12) THEN
             CALL CC_WRRSP(LISTD,IFILE,ISYRES,IOPTWR12,MODELW,
     &                     WORK(KDUM),WORK(KDUM),WORK(KDELTAR12),
     &                     WORK(KEND1),LEND1)
           END IF
           IF (CCSDT) THEN
             CALL CC_WRRSP(LISTD,IFILE,ISYRES,IOPTWE,MODELW,WORK(KDUM),
     &             WORK(KDELTA1(ITRAN)),WORK(KDELTA2),WORK(KEND1),LEND1)
           END IF
           TIMIO = TIMIO + SECOND() - DTIME

         ELSE IF (IOPTRES.EQ.4) THEN
C  
C          -----------------------------
C          write to file using CC_WARSP:
C          -----------------------------
C
           DTIME = SECOND()
           CALL CC_WARSP(LISTD,IFILE,ISYRES,IOPTW,MODELW,WORK(KDUM),
     &           WORK(KDELTA1(ITRAN)),WORK(KDELTA2),WORK(KEND1),LEND1)
           IF (CCSDR12) THEN
             CALL CC_WARSP(LISTD,IFILE,ISYRES,IOPTWR12,MODELW,
     &                     WORK(KDUM),WORK(KDUM),WORK(KDELTAR12),
     &                     WORK(KEND1),LEND1)
           END IF
           IF (CCSDT) THEN
             CALL CC_WARSP(LISTD,IFILE,ISYRES,IOPTWE,MODELW,WORK(KDUM),
     &             WORK(KDELTA1(ITRAN)),WORK(KDELTA2),WORK(KEND1),LEND1)
           END IF
           TIMIO = TIMIO + SECOND() - DTIME

         ELSE IF (IOPTRES.EQ.5) THEN
C  
C          --------------------------------
C          calculate required dot products:
C          --------------------------------
C
           IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KDELTA2),TWO,ISYRES)
           CALL CCDOTRSP(IGDOTS,GCON,IOPTW,LISTD,ITRAN,NGTRAN,MXVEC,
     &                   WORK(KDELTA1(ITRAN)),WORK(KDELTA2),ISYRES,
     &                   WORK(KEND1),LEND1)
           IF (CCSDR12) THEN
             CALL CCDOTRSP(IGDOTS,GCON,IOPTWR12,LISTD,ITRAN,NGTRAN,
     &                     MXVEC,WORK(KDUM),WORK(KDELTAR12),ISYRES,
     &                     WORK(KEND1),LEND1)
           END IF

         ELSE
            CALL QUIT('Illegal value for IOPTRES in CC_GMAT.')
         END IF

         TIMTRN(ITRAN) = TIMTRN(ITRAN) + SECOND()

         IF (IPRINT.GT.0) THEN
            IF (IOPTRES.EQ.5) THEN
               IVEC = 1
               DO  WHILE (IGDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                  IVEC = IVEC + 1
               END DO
               WRITE (LUPRI,'(1X,3(A,I5),A,I6,A,F10.2,A)') '| ',IDLSTA,
     &           '  | ',IDLSTB,'  | ',IDLSTC,'  | ',IVEC-1,
     &           '  | ',TIMTRN(ITRAN),'  |'
            ELSE
               WRITE (LUPRI,'(1X,3(A,I5),A,I6,A,F10.2,A)') '| ',IDLSTA,
     &           '  | ',IDLSTB,'  | ',IDLSTC,'  | ',IDLSTD,
     &           '  | ',TIMTRN(ITRAN),'  |'
            END IF
         END IF

      END DO ! ITRAN

*=====================================================================*
* finished with all G matrix transformations for this time:
* for IOPTRES = 1 and enough memory read result back into memory,
* in any case close G matrix file & return
*=====================================================================*

*     check size of work space:
      IF (IOPTRES .EQ. 1) THEN
         LENALL = IADRDEL - 1
         IF (LENALL .GT. LWORK) IOPTRES = 0
      END IF

*     read result vectors back into memory:
      IF (IOPTRES .EQ. 1) THEN
        DTIME = SECOND()
        CALL GETWA2(LUGMAT,FILGMA,WORK,1,LENALL)
        TIMIO = TIMIO + SECOND() - DTIME
      END IF

*     close G matrix file, if open:
      IF ( IOPTRES.EQ.0  .OR.  IOPTRES.EQ.1 ) THEN

         CALL WCLOSE2(LUGMAT,FILGMA,'KEEP')

      END IF

*---------------------------------------------------------------------*
* print timings and return:
*---------------------------------------------------------------------*
      TIMALL = SECOND() - TIMALL

      IF (IPRINT.GE.0) THEN
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
         WRITE (LUPRI,'(1X,A,I4,A,F10.2,A)')
     &     '| total time for',NGTRAN,' G transforms.:',TIMALL,' secs.|'
         IF (TIMALL.GT.1.0D0) THEN
          WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
          CONVRT = 100.0D0/TIMALL
          WRITE (LUPRI,
     &        '(1X,"|  % of time used in ",A,":",F10.2,"      |")')
     &      'ERI          ', TIMINT*CONVRT,
     &      'CCRDAO       ', TIMRDAO*CONVRT,
     &      'CCTRBT       ', TIMTRBT*CONVRT,
     &      'CCG_21CD     ', TIM21CD*CONVRT,
     &      'Fock interm. ', TIMFCK*CONVRT,
     &      'CCS part     ', TIMCCS*CONVRT,
     &      '21 CCSD part ', TIM21SD*CONVRT,
     &      'doubles part ', TIMDBL*CONVRT,
     &      'triples part ', TIMTRIP*CONVRT,
     &      'I/O          ', TIMIO*CONVRT
          END IF

          WRITE (LUPRI,'(1X,A1,50("="),A1,//)') '+','+'
      END IF

      CALL QEXIT('CC_GMAT')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_1ITROO */
*=====================================================================*
      SUBROUTINE CCG_1ITROO(TFMAT,ISYMTF,FMAT,ISYMFM,T1,ISYMT1)
*---------------------------------------------------------------------*
*
*     Purpose: calculate occupied/occupied block of one-index
*              transformed (Fock) matrix F:
*
*              tF(ij)  = sum(b) F(ib) * T1(bj) = [F,T1]_ij
*
*              needed, e.g., for double one-index transformed Fock mat.
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
   
      INTEGER ISYMTF, ISYMT1, ISYMFM
      
#if defined (SYS_CRAY)
      REAL TFMAT(*)  ! dimension (NMATIJ(ISYMTF))
      REAL FMAT(*)   ! dimension (NT1AM(ISYMFM))
      REAL T1(*)     ! dimension (NT1AM(ISYMT1))
#else
      DOUBLE PRECISION TFMAT(*)  ! dimension (NMATIJ(ISYMTF))
      DOUBLE PRECISION FMAT(*)   ! dimension (NT1AM(ISYMFM))
      DOUBLE PRECISION T1(*)     ! dimension (NT1AM(ISYMT1))
#endif

      INTEGER ISYMI, ISYMJ, ISYMB, NIJ, NBI, NBJ

      CALL QENTER('CCG_1ITROO')

* check symmetries:
      IF (ISYMTF .NE. MULD2H(ISYMFM,ISYMT1)) THEN
        CALL QUIT('Symmetry mismatch in CCG_1ITROO.')
      END IF

* initialize resulting matrix:
      CALL DZERO(TFMAT, NMATIJ(ISYMTF))

      DO ISYMI = 1, NSYM
        ISYMJ = MULD2H(ISYMI,ISYMTF)
        ISYMB = MULD2H(ISYMI,ISYMFM)

        DO I = 1, NRHF(ISYMI)
        DO J = 1, NRHF(ISYMJ)
          NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I

          DO B = 1, NVIR(ISYMB)
            NBI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1) + B
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
            TFMAT(NIJ) = TFMAT(NIJ) + FMAT(NBI) * T1(NBJ)
          END DO

        END DO
        END DO

      END DO
 
      CALL QEXIT('CCG_1ITROO')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_1ITRVO */
*=====================================================================*
      SUBROUTINE CCG_1ITRVO(tF,ISYMTF,FOO,FVV,ISYMFM,T1,ISYMT1,IOPT)
*---------------------------------------------------------------------*
*
*     Purpose: calculate vir/occ or occ/vir block of one-index
*              transformed (Fock) matrix F:
*
*  iopt=1 :    tF(ai)  = + sum_b F(ab) * T1(bi) - sum_j T1(aj) * F(ji)
*
*  iopt=2 :    tF(ia)  = + sum_b T1(ib) * F(ba) - sum_j F(ij) * T1(ja)
*
*              needed, e.g., for Delta1 vector in CCQR_G
*
*  iopt=1 referes to a one-index transformation of F with a vir/occ
*         matrix (e.g. single excitation amplitudes)
*
*  iopt=2 referes to a one-index transformation of F with a occ/vir
*         matrix (e.g. single excitation part of zeta vector)
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
   
      INTEGER ISYMTF, ISYMT1, ISYMFM, IOPT
      
#if defined (SYS_CRAY)
      REAL TF(*)      ! dimension (NT1AM(ISYMTF))
      REAL FOO(*)     ! dimension (NMATIJ(ISYMFM))
      REAL FVV(*)     ! dimension (NMATAB(ISYMFM))
      REAL T1(*)      ! dimension (NT1AM(ISYMT1))
#else
      DOUBLE PRECISION TF(*)    ! dimension (NT1AM(ISYMTF))
      DOUBLE PRECISION FOO(*)   ! dimension (NMATIJ(ISYMFM))
      DOUBLE PRECISION FVV(*)   ! dimension (NMATAB(ISYMFM))
      DOUBLE PRECISION T1(*)    ! dimension (NT1AM(ISYMT1))
#endif

      INTEGER ISYMI,ISYMJ,ISYMA,ISYMB,NAI,NAB,NBA,NBI,NAJ,NJI,NIJ

      CALL QENTER('CCG_1ITRVO')

* check symmetries:
      IF (ISYMTF .NE. MULD2H(ISYMFM,ISYMT1)) THEN
        CALL QUIT('Symmetry mismatch in CCG_1ITRVO.')
      END IF

* initialize resulting matrix:
      CALL DZERO(tF, NT1AM(ISYMTF))

      DO ISYMI = 1, NSYM
        ISYMA = MULD2H(ISYMI,ISYMTF)
        ISYMB = MULD2H(ISYMI,ISYMT1)
        ISYMJ = MULD2H(ISYMA,ISYMT1)

      IF (IOPT .EQ. 1) THEN

        DO I = 1, NRHF(ISYMI)
        DO A = 1, NVIR(ISYMA)
          NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A

          DO B = 1, NVIR(ISYMB)
            NBI = IT1AM(ISYMB,ISYMI)  + NVIR(ISYMB)*(I - 1) + B
            NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + A
            tF(NAI) = tF(NAI) + FVV(NAB) * T1(NBI)
          END DO

          DO J = 1, NRHF(ISYMJ)
            NAJ = IT1AM(ISYMA,ISYMJ)  + NVIR(ISYMA)*(J - 1) + A
            NJI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I - 1) + J
            tF(NAI) = tF(NAI) - T1(NAJ) * FOO(NJI)
          END DO

        END DO
        END DO

      ELSE IF (IOPT .EQ. 2) THEN

        DO I = 1, NRHF(ISYMI)
        DO A = 1, NVIR(ISYMA)
          NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A

          DO B = 1, NVIR(ISYMB)
            NBI = IT1AM(ISYMB,ISYMI)  + NVIR(ISYMB)*(I - 1) + B
            NBA = IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(A - 1) + B
            tF(NAI) = tF(NAI) + T1(NBI) * FVV(NBA)
          END DO

          DO J = 1, NRHF(ISYMJ)
            NAJ = IT1AM(ISYMA,ISYMJ)  + NVIR(ISYMA)*(J - 1) + A
            NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
            tF(NAI) = tF(NAI) - FOO(NIJ) * T1(NAJ)
          END DO

        END DO
        END DO

      ELSE
        CALL QUIT('CCG_1ITRVO called with illegal value for IOPT.')
      ENDIF

      END DO
 
      CALL QEXIT('CCG_1ITRVO')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_1ITRVV */
*=====================================================================*
      SUBROUTINE CCG_1ITRVV(TFMAT,ISYMTF,FMAT,ISYMFM,T1,ISYMT1)
*---------------------------------------------------------------------*
*
*     Purpose: calculate virtual/virtual block of one-index
*              transformed (Fock) matrix F:
*
*              tF(ab)  = - sum(i) T1(ai) * F(ib)  = [F,T1]_ab
*
*              needed, e.g., for double one-index transformed Fock mat.
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
   
      INTEGER ISYMTF, ISYMT1, ISYMFM
      
#if defined (SYS_CRAY)
      REAL TFMAT(*)
      REAL FMAT(*), T1(*)
#else
      DOUBLE PRECISION TFMAT(*)
      DOUBLE PRECISION FMAT(*), T1(*)
#endif

      INTEGER ISYMA, ISYMI, ISYMB, NAI, NBI, NAB

      CALL QENTER('CCG_1ITRVV')

* check symmetries:
      IF (ISYMTF .NE. MULD2H(ISYMFM,ISYMT1)) THEN
        CALL QUIT('Symmetry mismatch in CCG_1ITRVV.')
      END IF

* initialize resulting matrix:
      CALL DZERO(TFMAT, NMATAB(ISYMTF))

      DO ISYMA = 1, NSYM
        ISYMB = MULD2H(ISYMA,ISYMTF)
        ISYMI = MULD2H(ISYMA,ISYMT1)

        DO A = 1, NVIR(ISYMA)
        DO B = 1, NVIR(ISYMB)
          NAB  = IMATAB(ISYMA,ISYMB)  + NVIR(ISYMA)*(B - 1) + A

          DO I = 1, NRHF(ISYMI)
            NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
            NBI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1) + B
            TFMAT(NAB) = TFMAT(NAB) - T1(NAI) * FMAT(NBI) 
          END DO

        END DO
        END DO

      END DO
 
      CALL QEXIT('CCG_1ITRVV')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_21CCSD */
*=====================================================================*
       SUBROUTINE CCG_21CCSD (DELTA1, ISYRES,  ! out: Delta vector
     &                        LISTA,  IZETVA,  ! inp: A resp. Zeta vec. 
     &                        FockB,  ISYMFB,  ! inp: B transf. Fock Ma
     &                        T1AMPB, ISYMTB,  ! inp: resp. B singles
     &                        T2AMPC, ISYMTC,  ! inp: resp. C doubles
     &                        LISTC,  IDLSTC,  ! inp: resp. C 
     &                        XIAJK,  EIAJK,   ! out: transf. integrals
     &                        ZIAJK,  CIAJK,   ! out: transf. multipl.
     &                        WORK,   LWORK   )! work space
*---------------------------------------------------------------------*
*
*    Purpose: calculate CCSD contribution to DELTA1:
*
*        DELTA1 = DELTA1 + <Zeta_2| [[tH^B,T2C], tau_1] |HF>
*
*
*    symmetries/variables:
*              
*             ISYRES : result vector DELTA1
*             ISYCTR : lagrangian multipliers (zeta vector) CTR1, CTR2
*             ISYMTB : response amplitudes T1AMPB, (T2AmpB)
*             ISYMTC : response amplitudes T2AMPC, (T1AMPC)
*
*    Note: the double excitation response amplitudes T2AMPC are
*          expected in packed lower triangular storage 
*
*          the double excitation part of lagrangian multiplier vector
*          CTR2 and the (ia|jb) MO integrals will be read from disk
*
*          required vectors of the size (1/2) V^2 O^2 :
*            CCSD: (ia|jb), CTR2, T2AMPC
*
*          (ia|jb) integrals are assumed to have the symmetry ISYMOP,
*          which is transferred via COMMON /CCSDSYM/
*
*     Written by Christof Haettig, August 1996.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "r12int.h"
#include "ccr12int.h"

* local parameters:
      CHARACTER MSGDBG*(19)
      PARAMETER (MSGDBG='[debug] CC_21CCSD> ')
#if defined (SYS_CRAY)
      REAL ONE, ZERO, HALF, TWO
#else
      DOUBLE PRECISION ONE, ZERO, HALF, TWO
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, HALF = 0.5d0, TWO = 2.0d0)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
    

* input variables:
      CHARACTER*(*) LISTA, LISTC
      CHARACTER*(10) MODEL
      INTEGER ISYRES, ISYCTR, ISYMFB, ISYMTB, ISYMTC, LWORK, IZETVA
      INTEGER IDLSTC

#if defined (SYS_CRAY)
      REAL DELTA1(*), FOCKB(*), T1AMPB(*), T2AMPC(*)
      REAL XIAJK(*), EIAJK(*), ZIAJK(*), CIAJK(*)
      REAL WORK(LWORK)     
      REAL SWAP
#else
      DOUBLE PRECISION DELTA1(*), FOCKB(*), T1AMPB(*), T2AMPC(*)   
      DOUBLE PRECISION XIAJK(*), EIAJK(*), ZIAJK(*), CIAJK(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION SWAP
#endif

* local variables:
      LOGICAL LMOINT
      INTEGER ITAIKJ(8,8), NTAIKJ(8)
      INTEGER KEND0, KEND1, KEND1B, KEND2, KEND4
      INTEGER LEND0, LEND1, LEND1B, LEND2, LEND4
      INTEGER KXCMAT, KYCMAT, KDXYMAT, KXIAJB, KCTR2
      INTEGER KSCR, KDUM, KMCINT, KKCINT
      INTEGER ISYTCTB, ISYOVOV, ISYMKC, ISYMMC, ISYMXY, ISYDXY
      INTEGER ISYTCIN, ISYMC, ISYMJKL, NCI, KOFF, ISYMS, ISYMKLJ
      INTEGER ISYMI, ISYMJ, MAXJ, NIJ, NJI, IOPT, IPOS
      INTEGER ISYM, ICOUNT, ISYMAIK, ISYMB, ISYCTB, NTOTC
      INTEGER ISYMK, ISYIAJ, ISYMAI, ISYMBJ, NBJ, KOFF1, KOFF2, KOFF3
      INTEGER KT2AMP, ISYMDLK, KOFFAIK, KOFFDKL, ISYMAL, ISYML, ISYMA
      INTEGER KCINT, ISYMDKL, KEND1C, LEND1C, NTOTAI, NTOTB, NTOTJKL
      INTEGER KWTINT, KVTINT, KWHINT, KVHINT, KEXC, NTOTDKL, NTOTA
      INTEGER IOPTTCME, KTR12AM, KVINT, IAN, LUNIT 
   
#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

      INTEGER ILSTSYM

      CALL QENTER('CCG_21CCSD ')

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'just entered this routine...'
        WRITE (LUPRI,*) MSGDBG, 'LWORK:',LWORK
        Call FLSHFO(LUPRI)
      END IF

      DO  ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO 
        NTAIKJ(ISYM) = ICOUNT
      END DO

*---------------------------------------------------------------------*
* begin: set start of work space and set & check symmetries
*---------------------------------------------------------------------*
      ISYOVOV = ISYMOP                  ! 2-el integrals
      ISYCTR  = ILSTSYM(LISTA,IZETVA)   ! lagrangian multipliers CTR2
      ISYTCTB = MULD2H(ISYMTC,ISYMTB)   ! Wtilde, Vtilde
      ISYCTB  = MULD2H(ISYCTR,ISYMTB)   ! one-index transformed CTR2
      ISYMKC  = MULD2H(ISYOVOV,ISYMTC)  ! K^C intermediate
      ISYMMC  = MULD2H(ISYCTR,ISYMTC)   ! M^C intermediate

      IF (MULD2H(ISYRES,ISYOVOV) .NE. MULD2H(ISYCTR,ISYTCTB)) THEN
        CALL QUIT('Symmetry mismatch in CCG_21CCSD.')
      END IF

      KWTINT = 1
      KVTINT = KWTINT + NTAIKJ(ISYTCTB)
      KWHINT = KVTINT + NTAIKJ(ISYTCTB)
      KVHINT = KWHINT + NTAIKJ(ISYRES)
      KEND0  = KVHINT + NTAIKJ(ISYRES)

      KKCINT = KEND0
      KMCINT = KKCINT + N3ORHF(ISYMKC)
      KEND0  = KMCINT + N3ORHF(ISYMMC)

      LEND0  = LWORK - KEND0
      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCG_21CCSD. (0)')
      END IF


*---------------------------------------------------------------------*
* read (ia|jb) integrals from disk and square up:
*---------------------------------------------------------------------*
      KXIAJB = KEND0
      KEND1  = KXIAJB + NT2SQ(ISYOVOV)
      LEND1  = LWORK  - KEND1

      If (LEND1 .LT. NT2AM(ISYOVOV)) THEN
        CALL QUIT('Insufficient work space in CCG_21CCSD. (1a)')
      END IF

* read & unpack (ia|jb) integrals:
      Call CCG_RDIAJB(WORK(KEND1), NT2AM(ISYOVOV))

CCN
C     WRITE(LUPRI,*) '(ia|jb) integrals:'
C     CALL CC_PRP(WORK(KDUM),WORK(KEND1),ISYOVOV,0,1)
C     WRITE(LUPRI,*) 'T^C(el,dj) response amplitudes:'
C     CALL CC_PRP(WORK(KDUM),T2AMPC,ISYMTC,0,1)
CCN

      Call CC_T2SQ(WORK(KEND1), WORK(KXIAJB), ISYOVOV)

*---------------------------------------------------------------------*
* precalculate K^C intermediate: K(i,jkl) = sum_ed (ie|kd) T^C(el,dj)
* this and the (ia|jk) integrals are the only intermediates, which 
* require the (ia|jb) integrals as full square matrix:
*---------------------------------------------------------------------*
      ISYMKC  = MULD2H(ISYOVOV,ISYMTC)  ! K^C intermediate

      Call CC_MI(WORK(KKCINT),WORK(KXIAJB), ISYOVOV,
     &           T2AMPC,ISYMTC,WORK(KEND1),LEND1)

*---------------------------------------------------------------------*
* for CCSD(R12) (and higher) add Sum_mn (V^\dagger)_(im,kn) T^C(ml,nj)
*---------------------------------------------------------------------*
      IF (CCR12.AND..NOT.(CCS.OR.CC2)) THEN
        KTR12AM = KEND1
        KVINT   = KTR12AM + NTR12AM(ISYMTC)
        KEND1   = KVINT   + NTR12AM(1)
        LEND1   = LWORK - KEND1
        IF (LEND1.LT.0) THEN
          CALL QUIT('Insufficient work space in CCG_21CCSD. (1b)')
        END IF

        KDUM = - 99 999 999 ! dummy address
        IOPT = 32
        CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,WORK(KDUM),
     &                WORK(KTR12AM))

        LUNIT = -1
        CALL GPOPEN(LUNIT,FCCR12V,'OLD',' ','UNFORMATTED',
     &              KDUM,.FALSE.)
 9999   READ(LUNIT) IAN
        READ(LUNIT) (WORK(KVINT-1+I), I=1, NTR12AM(1))
        IF (IAN.NE.IANR12) GOTO 9999
        CALL GPCLOSE(LUNIT,'KEEP')

        IOPT = 2
        CALL CC_R12CV(WORK(KKCINT),WORK(KTR12AM),ISYMTC,WORK(KVINT),
     &                1,IOPT,WORK(KEND1),LEND1) 

      END IF

      IF (.FALSE.) THEN
        WRITE (LUPRI,*) MSGDBG,'response K intermediate:'
        WRITE(*,'(5f12.6)') (WORK(KKCINT-1+I),I=1,N3ORHF(ISYMKC))
      END IF

*---------------------------------------------------------------------*
* calculate (ia|jk) integrals, stored as I^j(ai,k)
*       and (ja|ik) integrals, stored as I^j(ai,k)
* finally replace (ia|jk) by (ia|jk) - 0.5*(ja|ik)
*---------------------------------------------------------------------*
      CALL CCG_TRANS4(XIAJK,ISYMTB,WORK(KXIAJB),ISYOVOV,T1AMPB,ISYMTB)

      CALL CCG_SORT1(XIAJK,EIAJK,ISYMTB,0)

      CALL DAXPY(NTAIKJ(ISYMTB),-HALF,EIAJK,1,XIAJK,1)

*---------------------------------------------------------------------*
* read & unpack lagrangian multipliers and re-read packed integrals:
*---------------------------------------------------------------------*
* 1) squared multipliers & scratch space for packed multipliers:
      KCTR2  = KEND0
      KSCR   = KCTR2  + NT2SQ(ISYCTR)
      KEND1B = KSCR   + NT2AM(ISYCTR)
      LEND1B = LWORK - KEND1B

* 2) squared multipliers & packed integrals:
      KXIAJB = KCTR2  + NT2SQ(ISYCTR)
      KEND1  = KXIAJB + NT2AM(ISYOVOV)
      LEND1  = LWORK - KEND1

      If (LEND1B .LT. 0 .OR. LEND1 .LT.0) THEN
        CALL QUIT('Insufficient work space in CCG_21CCSD. (1b/1)')
      END IF

* read packed lagrangian multipliers and square them up
      KDUM = - 99 999 999 ! dummy address

      IOPT = 2
      Call CC_RDRSP(LISTA,IZETVA,ISYCTR,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KSCR))

      Call CC_T2SQ(WORK(KSCR), WORK(KCTR2), ISYCTR)

* read packed integrals (ia|jb) and calculate in place L(iajb):
      CALL CCG_RDIAJB(WORK(KXIAJB), NT2AM(ISYOVOV))
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KXIAJB),ONE,ISYOVOV,IOPTTCME)

*---------------------------------------------------------------------*
* calculate X^C and Y^C intermediates and evaluate A+B and E terms:
*---------------------------------------------------------------------*
      ISYMXY  = MULD2H(ISYCTR,ISYMTC)
      ISYDXY  = MULD2H(ISYCTR,ISYTCTB)

      KXCMAT  = KEND1
      KYCMAT  = KXCMAT  + NMATIJ(ISYMXY)
      KDXYMAT = KYCMAT  + NMATAB(ISYMXY)
      KSCR    = KDXYMAT + NT1AM(ISYDXY)
      KEND2   = KSCR    + NT1AM(ISYRES)
      LEND2   = LWORK - KEND2

      IF (LEND2 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCG_21CCSD. (A+B)')
      END IF

* calculate X^C & Y^C intermediate:
      Call CC_XI(WORK(KXCMAT),WORK(KCTR2), ISYCTR,
     &                        T2AMPC,ISYMTC,WORK(KEND2),LEND2)

      Call CC_YI(WORK(KYCMAT),WORK(KCTR2), ISYCTR,
     &                        T2AMPC,ISYMTC,WORK(KEND2),LEND2)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'response X intermediate:'
        WRITE (LUPRI,'(5f12.6)') (WORK(KXCMAT-1+I),I=1,NMATIJ(ISYMXY))
        WRITE (LUPRI,*) MSGDBG,'response Y intermediate:'
        WRITE (LUPRI,'(2f12.6)') (WORK(KYCMAT-1+I),I=1,NMATAB(ISYMXY))
      END IF

* calculate XY^C:  XY^C_ij = X^C_ji,  XY^C_bd = -Y^C_bd
* 1.) transpose X^C intermediate
      DO ISYMI = 1, NSYM
        ISYMJ = MULD2H(ISYMI,ISYMXY)
        IF (ISYMJ .LE. ISYMI) THEN
          DO I = 1, NRHF(ISYMI)
            MAXJ =  NRHF(ISYMJ)
            IF (ISYMJ .EQ. ISYMI) MAXJ = I-1
          DO J = 1, MAXJ
            NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + I
            NJI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I-1) + J
            SWAP = WORK(KXCMAT-1+NIJ)
            WORK(KXCMAT-1+NIJ) = WORK(KXCMAT-1+NJI)
            WORK(KXCMAT-1+NJI) = SWAP
          END DO
          END DO
        END IF
      END DO

* 2.) multiply Y^C intermediate with -1:
      Call DSCAL(NMATAB(ISYMXY), -ONE, WORK(KYCMAT), 1)

*.....................................................................*
* A+B term: 
*.....................................................................*
* calculate effective density matrix DXY = [XY^C,T1^B]
      IOPT  = 1
      Call CCG_1ITRVO(WORK(KDXYMAT),ISYDXY,WORK(KXCMAT),WORK(KYCMAT),
     &                ISYMXY,T1AMPB,ISYMTB,IOPT)

* contract with L(ia|jb):
      Call CCG_LXD(WORK(KSCR), ISYRES, WORK(KDXYMAT), ISYDXY, 
     &             WORK(KXIAJB),ISYOVOV,0)

* add contribution to DELTA1
      Call DAXPY (NT1AM(ISYRES), ONE, WORK(KSCR),1, DELTA1, 1)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'T1AMPB:'
        WRITE (LUPRI,'(5f12.6)') (T1AMPB(I),I=1,NT1AM(ISYMTB))
        WRITE (LUPRI,*) MSGDBG,'DXY intermediate:'
        WRITE (LUPRI,'(5f12.6)') (WORK(KDXYMAT-1+I),I=1,NT1AM(ISYRES))
        WRITE (LUPRI,*) MSGDBG, 'A+B term:'
        CALL CC_PRP(WORK(KSCR),WORK,ISYRES,1,0)
      END IF

*.....................................................................*
* E term: 
*.....................................................................*

* do one-index transformation of XY^C with F^B:
      IOPT  = 2
      Call CCG_1ITRVO(WORK(KSCR),ISYRES,
     &                WORK(KXCMAT),WORK(KYCMAT),ISYMXY,
     &                FockB,       ISYMFB,      IOPT    )

* add contribution to DELTA1
      Call DAXPY (NT1AM(ISYRES), ONE, WORK(KSCR),1, DELTA1, 1)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'FockB:'
        WRITE (LUPRI,'(5f12.6)') (FockB(I),I=1,NT1AM(ISYMFB))
        WRITE (LUPRI,*) MSGDBG, 'E term:'
        CALL CC_PRP(WORK(KSCR),WORK,ISYRES,1,0)
        WRITE (LUPRI,*) MSGDBG, 'Delta1 now:'
        CALL CC_PRP(DELTA1,WORK,ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*

*---------------------------------------------------------------------*
* precalculate M^C intermediate: M(i,jkl) = sum_ed Z(ei,dk) T^C(el,dj)
* (overwrites X^C and Y^C intermediates)
*---------------------------------------------------------------------*
      ISYMMC = MULD2H(ISYCTR,ISYMTC)

      Call CC_MI(WORK(KMCINT),WORK(KCTR2),ISYCTR, 
     &           T2AMPC,ISYMTC,WORK(KEND1),LEND1)

      IF (.FALSE.) THEN
        WRITE (LUPRI,*) MSGDBG,'response M intermediate:'
        WRITE(LUPRI,'(5f12.6)') (WORK(KMCINT-1+I),I=1,N3ORHF(ISYMMC))
      END IF

*---------------------------------------------------------------------*
* calculate Zeta(ia|jk), stored as Z^j(ai,k)
*       and Zeta(ja|ik), stored as C^j(ai,k)
* finally replace Zeta(ja|ik) by 2 Zeta(ja|ik) + Zeta(ia|jk)
*---------------------------------------------------------------------*
      CALL CCG_TRANS4(ZIAJK,ISYCTB,WORK(KCTR2),ISYCTR,T1AMPB,ISYMTB)

      CALL CCG_SORT1(ZIAJK,CIAJK,ISYCTB,0)

      CALL DAXPY(NTAIKJ(ISYCTB),HALF,ZIAJK,1,CIAJK,1)

*---------------------------------------------------------------------*
*  calculate Wtilde and What intermediates:
*    Wtilde(dl,k;j) = sum_ia [2(ia|jk) - (ja|ik)] [2T(ai,dl)-T(al,di)]
*    What(dl,k;j)   = sum_ia Z(ia|jk) [2T(ai,dl)-T(al,di)]
*  calculate Vtilde and Vhat intermediates:
*    Vhat(dl,k;j)   = sum_ia [2Z(ja|ik)+Z(ia|jk)] T(al,di)
*    Vtilde(dl,k;j) = sum_ia (ja|ik) T(al,di)
*---------------------------------------------------------------------*
      KT2AMP = KEND0
      KEND1  = KT2AMP + NT2SQ(ISYMTC)
      LEND1  = LWORK  - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCG_21CCSD. (1b)')
      END IF

* square up T2 amplitudes and calculate in place 2 Cou - Exc comb.:
      CALL CC_T2SQ(T2AMPC,WORK(KT2AMP),ISYMTC)
      CALL CCRHS_T2TR(WORK(KT2AMP),WORK(KEND1),LEND1,ISYMTC)

* Wtilde intermediate:
      DO ISYMJ = 1, NSYM
         ISYMAIK = MULD2H(ISYMTB,ISYMJ)
         ISYMDLK = MULD2H(ISYMAIK,ISYMTC)

         DO J = 1, NRHF(ISYMJ)
           KOFFAIK = 1+ITAIKJ(ISYMAIK,ISYMJ)+NT2BCD(ISYMAIK)*(J-1)
           KOFFDKL = ITAIKJ(ISYMDLK,ISYMJ) + NT2BCD(ISYMDLK)*(J-1)

           IOPT = 1
           CALL CC_ZWVI(WORK(KWTINT+KOFFDKL),WORK(KT2AMP),ISYMTC,
     &                  XIAJK(KOFFAIK),ISYMAIK,WORK(KEND1),LEND1,IOPT)
         END DO
      END DO

* What intermediate:
      DO ISYMJ = 1, NSYM
         ISYMAIK = MULD2H(ISYCTB,ISYMJ)
         ISYMDLK = MULD2H(ISYMAIK,ISYMTC)

         DO J = 1, NRHF(ISYMJ)
           KOFFAIK = 1+ITAIKJ(ISYMAIK,ISYMJ)+NT2BCD(ISYMAIK)*(J-1)
           KOFFDKL = ITAIKJ(ISYMDLK,ISYMJ) + NT2BCD(ISYMDLK)*(J-1)

           IOPT = 1
           CALL CC_ZWVI(WORK(KWHINT+KOFFDKL),WORK(KT2AMP),ISYMTC,
     &                  ZIAJK(KOFFAIK),ISYMAIK,WORK(KEND1),LEND1,IOPT)
         END DO
      END DO


* recover squared T2 amplitudes and transpose occupied indeces:
      CALL CC_T2SQ(T2AMPC,WORK(KT2AMP),ISYMTC)
      CALL CCSD_T2TP(WORK(KT2AMP),WORK(KEND1),LEND1,ISYMTC)

* Vtilde intermediate:
      DO ISYMJ = 1, NSYM
         ISYMAIK = MULD2H(ISYMTB,ISYMJ)
         ISYMDLK = MULD2H(ISYMAIK,ISYMTC)

         DO J = 1, NRHF(ISYMJ)
           KOFFAIK = 1+ITAIKJ(ISYMAIK,ISYMJ)+NT2BCD(ISYMAIK)*(J-1)
           KOFFDKL = ITAIKJ(ISYMDLK,ISYMJ) + NT2BCD(ISYMDLK)*(J-1)

           IOPT = 1
           CALL CC_ZWVI(WORK(KVTINT+KOFFDKL), WORK(KT2AMP), ISYMTC,
     &                  EIAJK(KOFFAIK),ISYMAIK,WORK(KEND1),LEND1,IOPT)
         END DO
      END DO

* Vhat intermediate:
      DO ISYMJ = 1, NSYM
         ISYMAIK = MULD2H(ISYCTB,ISYMJ)
         ISYMDLK = MULD2H(ISYMAIK,ISYMTC)

         DO J = 1, NRHF(ISYMJ)
           KOFFAIK = 1+ITAIKJ(ISYMAIK,ISYMJ)+NT2BCD(ISYMAIK)*(J-1)
           KOFFDKL = ITAIKJ(ISYMDLK,ISYMJ) + NT2BCD(ISYMDLK)*(J-1)

           IOPT = 1
           CALL CC_ZWVI(WORK(KVHINT+KOFFDKL), WORK(KT2AMP), ISYMTC,
     &                  CIAJK(KOFFAIK),ISYMAIK,WORK(KEND1),LEND1,IOPT)
         END DO
      END DO

*---------------------------------------------------------------------*
* reorganization of work space:
*  KKCINT :  K^C intermediate
*  KMCINT :  M^C intermediate
*  KCTR2  :  zeta vector (packed)
*  KXIAJB :  (ia|jb) integrals (packed)
*  KEND2  :  free work space
*---------------------------------------------------------------------*
      KCTR2  = KEND0
      KXIAJB = KCTR2  + NT2AM(ISYCTR)
      KEND2  = KXIAJB + NT2AM(ISYOVOV)
      LEND2  = LWORK  - KEND2

      IF (LEND2 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCG_21CCSD.')
      END IF

* read (ia|jb) integrals:
      Call CCG_RDIAJB(WORK(KXIAJB), NT2AM(ISYOVOV))

* read packed lagrangian multipliers and square them up
      IOPT = 2
      Call CC_RDRSP(LISTA,IZETVA,ISYCTR,IOPT,MODEL,
     &              WORK(KDUM),WORK(KCTR2))

*---------------------------------------------------------------------*
* H term: contract Wtilde/Vtilde interm. with the Lagrangian multipl.:
*---------------------------------------------------------------------*
      DO ISYMA = 1, NSYM

         ISYMDKL = MULD2H(ISYCTR,ISYMA)
         ISYMI   = MULD2H(ISYTCTB,ISYMDKL)

         KCINT  = KEND2
         KEXC   = KCINT + NT2BCD(ISYMDKL)
         KEND1C = KEXC  + NT2BCD(ISYMDKL)
         LEND1C = LWORK - KEND1C
        
         IF (LEND1C .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCG_21CCSD. (1c)')
         END IF

         DO A = 1, NVIR(ISYMA)
           
            CALL CCG_TI(WORK(KCINT),ISYMDKL,WORK(KCTR2),ISYCTR,A,ISYMA)

            CALL DCOPY(NT2BCD(ISYMDKL),WORK(KCINT),1,WORK(KEXC),1)
            CALL CCLT_P21I(WORK(KEXC),ISYMDKL,WORK(KEND1C),LEND1C,
     &                     IT2BCD, NT2BCD, IT1AM, NT1AM, NVIR)
            CALL DAXPY(NT2BCD(ISYMDKL),HALF,WORK(KCINT),1,WORK(KEXC),1)

            KOFF1   = ITAIKJ(ISYMDKL,ISYMI) 
            KOFF2   = IT1AM(ISYMA,ISYMI) + A 
            NTOTDKL = MAX(NT2BCD(ISYMDKL),1)
            NTOTA   = MAX(NVIR(ISYMA),1)

            CALL DGEMV('T',NT2BCD(ISYMDKL),NRHF(ISYMI),
     &                 -ONE,WORK(KWTINT+KOFF1),NTOTDKL,WORK(KCINT),1,
     &                  ONE,DELTA1(KOFF2),NTOTA)

            CALL DGEMV('T',NT2BCD(ISYMDKL),NRHF(ISYMI),
     &                  ONE,WORK(KVTINT+KOFF1),NTOTDKL,WORK(KEXC),1,
     &                  ONE,DELTA1(KOFF2),NTOTA)

         END DO
      END DO


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'Delta1 after H term in CCG_21CCSD:'
        Call CC_PRP(DELTA1,WORK,ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*
* I term: contract What/Vhat intermediates with the integrals:
*---------------------------------------------------------------------*
      DO ISYMA = 1, NSYM

         ISYMDKL = MULD2H(ISYOVOV,ISYMA)
         ISYMI   = MULD2H(ISYRES,ISYMDKL)

         KCINT  = KEND2
         KEXC   = KCINT + NT2BCD(ISYMDKL)
         KEND1C = KEXC  + NT2BCD(ISYMDKL)
         LEND1C = LWORK - KEND1C
        
         IF (LEND1C .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCG_21CCSD. (1c)')
         END IF

         DO A = 1, NVIR(ISYMA)
           
           CALL CCG_TI(WORK(KCINT),ISYMDKL,WORK(KXIAJB),ISYOVOV,A,ISYMA)

           CALL DCOPY(NT2BCD(ISYMDKL),WORK(KCINT),1,WORK(KEXC),1)
           CALL CCLT_P21I(WORK(KEXC),ISYMDKL,WORK(KEND1C),LEND1C,
     &                    IT2BCD, NT2BCD, IT1AM, NT1AM, NVIR)
           CALL DAXPY(NT2BCD(ISYMDKL),-HALF,WORK(KEXC),1,WORK(KCINT),1)

           KOFF1   = ITAIKJ(ISYMDKL,ISYMI) 
           KOFF2   = IT1AM(ISYMA,ISYMI) + A
           NTOTDKL = MAX(NT2BCD(ISYMDKL),1)
           NTOTA   = MAX(NVIR(ISYMA),1)

           CALL DGEMV('T',NT2BCD(ISYMDKL),NRHF(ISYMI),
     &                -ONE,WORK(KWHINT+KOFF1),NTOTDKL,WORK(KCINT),1,
     &                 ONE,DELTA1(KOFF2),NTOTA)

           CALL DGEMV('T',NT2BCD(ISYMDKL),NRHF(ISYMI),
     &                 ONE,WORK(KVHINT+KOFF1),NTOTDKL,WORK(KEXC),1,
     &                 ONE,DELTA1(KOFF2),NTOTA)

         END DO
      END DO

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'Delta1 after I term in CCG_21CCSD:'
        Call CC_PRP(DELTA1,WORK,ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*
* resort transformed integrals and multipliers for G and F term:
*---------------------------------------------------------------------*

      CALL DAXPY(NTAIKJ(ISYMTB),HALF,EIAJK,1,XIAJK,1)
      CALL CCG_SORT1(XIAJK,EIAJK,ISYMTB,4)

      CALL CCG_SORT1(ZIAJK,CIAJK,ISYCTB,4)

*---------------------------------------------------------------------*
* loop over symmetry blocks and calculate G and F term contributions:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'Delta1 before G and F term:'
        Call CC_PRP(DELTA1,WORK,ISYRES,1,0)
      END IF
      DO ISYMC = 1, NSYM

          ISYMI   = MULD2H(ISYMC,ISYRES)
          NCI     = IT1AM(ISYMC,ISYMI) + 1
          NTOTC   = MAX(NVIR(ISYMC),1)
  
C         --------------------------------------
C         G term: contract (k l^B|j c) integrals 
C                with M^C(i,klj) intermediate:
C         --------------------------------------
          ISYMJKL = MULD2H(ISYMI,ISYMMC)

          KOFF1   = KMCINT + I3ORHF(ISYMJKL,ISYMI) 
          KOFF2   = ISJIKA(ISYMJKL,ISYMC) + 1
          NTOTJKL = MAX(NMAIJK(ISYMJKL),1)

          CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMI),NMAIJK(ISYMJKL),
     &               ONE,EIAJK(KOFF2),NTOTC,WORK(KOFF1),NTOTJKL,
     &               ONE,DELTA1(NCI), NTOTC)
 
C         ----------------------------------------
C         F term: contract S^{B,c}(klj) multiplier
C                 with K^C(i,klj) intermediate
C         ----------------------------------------
          ISYMS   = MULD2H(ISYCTR,ISYMTB)
          ISYMKLJ = MULD2H(ISYMC,ISYMS)
  
          KOFF1   = KKCINT + I3ORHF(ISYMKLJ,ISYMI)  
          KOFF2   = ISJIKA(ISYMKLJ,ISYMC) + 1
          NTOTJKL = MAX(NMAIJK(ISYMKLJ),1)

          CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMI),NMAIJK(ISYMKLJ),
     &               ONE,CIAJK(KOFF2),NTOTC,WORK(KOFF1),NTOTJKL,
     &               ONE,DELTA1(NCI), NTOTC)
CCN
C       WRITE(LUPRI,*) 'K^C(i,klj):'
C       CALL OUTPUT(WORK(KOFF1),1,NTOTJKL,1,NRHF(ISYMI),
C    &              NTOTJKL,NRHF(ISYMI),1,LUPRI)
C       WRITE(LUPRI,*) 'CIAJK:'
C       CALL OUTPUT(CIAJK(KOFF2),1,NVIR(ISYMC),
C    &              1,NMAIJK(ISYMKLJ),NVIR(ISYMC),NMAIJK(ISYMKLJ),
C    &              1,LUPRI)
C       WRITE (LUPRI,*) MSGDBG, 'after Symmetry block:',ISYMC
C       Call CC_PRP(DELTA1,WORK,ISYRES,1,0)
CCN
      END DO ! ISYMC 

CCN
C     WRITE(LUPRI,*) 'S^{B,c}(klj): '
C     DO ISYMC = 1, NSYM
C       ISYMS   = MULD2H(ISYCTR,ISYMTB)
C       ISYMKLJ = MULD2H(ISYMC,ISYMS)
C       WRITE(LUPRI,*) 'Symmetry block: ',ISYMKLJ,ISYMC
C       CALL OUTPUT(CIAJK(1+ISJIKA(ISYMKLJ,ISYMC)),1,NVIR(ISYMC),
C    &              1,NMAIJK(ISYMKLJ),NVIR(ISYMC),NMAIJK(ISYMKLJ),
C    &              1,LUPRI)
C     END DO
CCN

*---------------------------------------------------------------------*
* print debug output and return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'Delta1 at the end of CCG_21CCSD:'
        Call CC_PRP(DELTA1,WORK,ISYRES,1,0)
      END IF

      CALL QEXIT('CCG_21CCSD ')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_21CD */
*=====================================================================*
      SUBROUTINE CCG_21CD( RESVEC,ISYRES,CTR2,ISYCTR,XLAMDH0,XLAMDP0,
     &                     XLAMDHA,XLAMDPA,ISYMLA,
     &                     XLAMDHB,XLAMDPB,ISYMLB,
     &                     DSRHFC,ISYMXC,IDEL,ISYDEL,
     &                     LSAME, FACTOR,WORK,LWORK)
*---------------------------------------------------------------------*
*
*  Purpose: to calculate G and K matrix contributions which are
*           derivatives of the 21C and 21D contribution to RHO1
*           (the left hand side transformed trial vector).
*           for that this routine has to be called 3 times with
*           appropriate permutations of the XLAMDH/XLAMDP matrices
* 
*
*
*  21C contribution:
*
*  RESVEC(a i) =  RESVEC(a i) - FACTOR * XLAMDP0(del i) * sum_{k j c} 
*           [(del j^C|c^B k^A) + (del j^C|c^A k^B)] * CTR2(a j c k)
*
*
*  21D contribution:
*  
*  RESVEC(a i) = RESVEC(a i) + sum_{bet} FACTOR * XLAMDH0(bet a) * 
*                 sum_{j alp} (del j^C|alp bet) * 
*                    [ CTR2(alp^B i j del^A) + CTR2(alp^A i j delB)]
*  
*   
*  For LSAME=.true. XLAMDHA & XLAMDHB and XLAMDPA & XLAMDPB are assumed
*  to be the same --> some intermediates are multiplied by 2 instead
*  of calculating them again with A and B exchanged...
*
*   symmetries & variables:   
*
*            ISYRES : result vector RESVEC
*            ISYCTR : lagrangian multipliers (zeta vector) CTR2
*            ISYMXC : symmetry of DSRHFC
*            ISYDEL : symmetry of SAO basis function IDEL
*            ISYMLA : Lambda matrices XLAMDHA  XLAMDPA
*            ISYMLB : Lambda matrices XLAMDHB, XLAMDPB
*                     Lambda matrices XLAMDH0, XLAMDP0 assumed tot.sym.
*
*   used for CC2 part of CCQR_G
*
*   Christof Haettig, August 1996, revised in August 1998
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "dummy.h"
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
 

#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO, FACT, FACTOR
#else
      DOUBLE PRECISION ZERO, ONE, TWO, FACT, FACTOR
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0)

      LOGICAL LSAME
      INTEGER ISYRES, ISYCTR, ISYDEL, ISYMLA, ISYMLB, ISYMXC
      INTEGER IDEL, LWORK

#if defined (SYS_CRAY)
      REAL RESVEC(*)    ! dimension (NT1AM(ISYRES))
      REAL CTR2(*)      ! dimension (NT2AM(ISYCTR))
      REAL DSRHFC(*)
      REAL XLAMDH0(*)   ! dimension (NLAMDT)
      REAL XLAMDP0(*)   ! dimension (NLAMDT)
      REAL XLAMDHA(*)   ! dimension (NGLMDT(ISYMLA))
      REAL XLAMDPA(*)   ! dimension (NGLMDT(ISYMLA))
      REAL XLAMDHB(*)   ! dimension (NGLMDT(ISYMLB))
      REAL XLAMDPB(*)   ! dimension (NGLMDT(ISYMLB))
      REAL WORK(LWORK)
      REAL DDOT
#else
      DOUBLE PRECISION RESVEC(*)    ! dimension (NT1AM(ISYRES))
      DOUBLE PRECISION CTR2(*)      ! dimension (NT2AM(ISYCTR))
      DOUBLE PRECISION DSRHFC(*)
      DOUBLE PRECISION XLAMDH0(*)   ! dimension (NLAMDT)
      DOUBLE PRECISION XLAMDP0(*)   ! dimension (NLAMDT)
      DOUBLE PRECISION XLAMDHA(*)   ! dimension (NGLMDT(ISYMLA))
      DOUBLE PRECISION XLAMDPA(*)   ! dimension (NGLMDT(ISYMLA))
      DOUBLE PRECISION XLAMDHB(*)   ! dimension (NGLMDT(ISYMLB))
      DOUBLE PRECISION XLAMDPB(*)   ! dimension (NGLMDT(ISYMLB))
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION DDOT
#endif

      INTEGER ISYLALB, ISYMAO, ISYMMO
      INTEGER KDRES, KCRES, KEND0, LEND0, KEND1, LEND1
      INTEGER ISYCTA, KCTR2A, IOPT, KLAM0, KOFFR, ISYMBET, NTOTA
      INTEGER ISYMJ, ISYALBE, KSQC, KEND2, LEND2
      INTEGER KOFFC, KOFFB, ISYMALP, ISYMC, ISYMK, KHALF, KMO
      INTEGER KOFFIC, KOFFLB, KOFFLA
      INTEGER ISYMCI, ISYMI, KOFFZ, KOFFSV, NTOTC, NTOTBE, NTOTAL
      INTEGER ISYMCK, ISYMAJ, ISYMA, NAJ, KCTR2, ISYMD, KRESV
      INTEGER NTOTB, ISYMAI, KOFFZC, KOFFZB, ISYCTAB
      INTEGER KCTR2AB, NHALF, ISYCTB
      INTEGER KMOINT


C      INTEGER INDEX
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_21CD')

* set & check symmetries:
      ISYLALB = MULD2H(ISYMLA,ISYMLB)    ! Lambda C x Lambda B
      ISYMAO  = MULD2H(ISYMXC,ISYDEL)    ! AO integrals
      ISYMMO  = MULD2H(ISYMAO,ISYLALB)   ! MO integrals

      IF (ISYRES .NE. MULD2H(ISYMMO,ISYCTR)) THEN
        CALL QUIT('Symmetry mismatch in CCG_21CD')
      END IF

* allocate scratch arrays for half-transformed result vectors:
      KDRES = 1
      KCRES = KDRES + NT1AO(ISYRES)
      KEND1 = KCRES + NVIRT
      LEND1 = LWORK - KEND1

      Call DZERO(WORK(KDRES),NT1AO(ISYRES))
      Call DZERO(WORK(KCRES),NVIRT)

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCG_21CD (0).')
      END IF

*---------------------------------------------------------------------*
* calculate one-index back-transformed CTR2 vector using XLAMDPA:
*---------------------------------------------------------------------*
      ISYCTA  = MULD2H(ISYDEL,MULD2H(ISYCTR,ISYMLA))
      ISYCTB  = MULD2H(ISYDEL,MULD2H(ISYCTR,ISYMLB))
      ISYCTAB = MULD2H(ISYCTA,ISYMLB)
      
      KCTR2A   = KEND1
      KCTR2AB  = KCTR2A  + MAX(NT2BCD(ISYCTA),NT2BCD(ISYCTB))
      KEND1    = KCTR2AB + NT2BGD(ISYCTAB)
      LEND1    = LWORK - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCG_21CD (1b).')
      END IF

      D = IDEL - IBAS(ISYDEL)

* back transform the virtual indeces of CTR2 with XLAMDPA/XLAMDPB:
      IOPT = 1 
      Call CC_T2AO(CTR2, XLAMDPA, ISYMLA, WORK(KCTR2A),
     &             WORK(KEND1), LEND1, IDEL, ISYDEL, ISYCTR, IOPT)

      IOPT = 0
      CALL CC_BFAO(WORK(KCTR2AB),WORK(KCTR2A),ISYCTA,XLAMDPB,ISYMLB,
     &             D,ISYDEL,DUMMY,IDUMMY,DUMMY,IDUMMY,ZERO,IOPT)

      IF (LSAME) THEN
         CALL DSCAL(NT2BGD(ISYCTAB),TWO,WORK(KCTR2AB),1)
      ELSE 
         IOPT = 1 
         Call CC_T2AO(CTR2, XLAMDPB, ISYMLB, WORK(KCTR2A),
     &                WORK(KEND1), LEND1, IDEL, ISYDEL, ISYCTR, IOPT)
         IOPT = 0
         CALL CC_BFAO(WORK(KCTR2AB),WORK(KCTR2A),ISYCTB,XLAMDPA,ISYMLA,
     &                D,ISYDEL,DUMMY,IDUMMY,DUMMY,IDUMMY,ONE,IOPT)
      END IF

*---------------------------------------------------------------------*
* start loop over occupied index j:
*---------------------------------------------------------------------*
      DO ISYMJ = 1, NSYM
        ISYALBE = MULD2H(ISYMJ,ISYMXC)
        ISYMCK  = MULD2H(ISYALBE,MULD2H(ISYMLA,ISYMLB))

        KSQC   = KEND1
        KMOINT = KSQC + N2BST(ISYALBE)
        KEND2  = KMOINT + NT1AM(ISYMCK)
        LEND2  = LWORK - KEND2

        IF (LEND2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCG_21CD (2).')
        END IF


      DO J = 1, NRHF(ISYMJ)

        CALL DZERO(WORK(KMOINT),NT1AM(ISYMCK))

C       ---------------------------------------------
C       unpack integral distribution (al be | j del):
C       ---------------------------------------------
        KOFFC = IDSRHF(ISYALBE,ISYMJ) + NNBST(ISYALBE)*(J-1) + 1

        Call CCSD_SYMSQ (DSRHFC(KOFFC), ISYALBE, WORK(KSQC))


C       ----------------------------------------------------
C       loop over symmetry blocks of unpacked distributions:
C       ----------------------------------------------------
        DO ISYMBET = 1, NSYM
          ISYMALP = MULD2H(ISYMBET,ISYALBE)

*---------------------------------------------------------------------*
* calculate (c^B k^A |j^C del) + (c^A k^B |j^c del) integrals:
*---------------------------------------------------------------------*
          ISYMC  = MULD2H(ISYMALP,ISYMLB)
          ISYMK  = MULD2H(ISYMBET,ISYMLA)

          NHALF = NBAS(ISYMALP) * NRHF(ISYMK)

          IF (LEND2 .LT. NHALF) THEN
            CALL QUIT('Insufficient work space in CCG_21CD (3).')
          END IF

* first XLAMDPB(alp c) * XLAMDHA(bet k) * (alp bet|j^C del):
          KOFFIC = KSQC + IAODIS(ISYMALP,ISYMBET)
          KOFFLA = IGLMRH(ISYMBET,ISYMK) + 1
          KOFFLB = IGLMVI(ISYMALP,ISYMC) + 1
          KMO    = KMOINT + IT1AM(ISYMC,ISYMK)
          KHALF  = KEND2
        
          NTOTAL = MAX(NBAS(ISYMALP),1)
          NTOTBE = MAX(NBAS(ISYMBET),1)
          NTOTC  = MAX(NVIR(ISYMC),1)
             
          Call DGEMM('N','N',NBAS(ISYMALP),NRHF(ISYMK),NBAS(ISYMBET),
     &               ONE,  WORK(KOFFIC),NTOTAL, XLAMDHA(KOFFLA),NTOTBE,
     &               ZERO, WORK(KHALF),NTOTAL)

          Call DGEMM('T','N',NVIR(ISYMC), NRHF(ISYMK), NBAS(ISYMALP),
     &               ONE, XLAMDPB(KOFFLB),NTOTAL, WORK(KHALF),NTOTAL,
     &               ONE, WORK(KMO),NTOTC )


* add XLAMDPA(alp c) * XLAMDHB(bet k) * (alp bet|j^C del):
         IF (.NOT. LSAME) THEN

           ISYMC  = MULD2H(ISYMALP,ISYMLA)
           ISYMK  = MULD2H(ISYMBET,ISYMLB)

           NHALF = NBAS(ISYMALP) * NRHF(ISYMK)

           IF (LEND2 .LT. NHALF) THEN
             CALL QUIT('Insufficient work space in CCG_21CD (3).')
           END IF

           KOFFIC = KSQC + IAODIS(ISYMALP,ISYMBET)
           KOFFLB = IGLMRH(ISYMBET,ISYMK) + 1
           KOFFLA = IGLMVI(ISYMALP,ISYMC) + 1
           KMO    = KMOINT + IT1AM(ISYMC,ISYMK)
           KHALF  = KEND2
        
           NTOTAL = MAX(NBAS(ISYMALP),1)
           NTOTBE = MAX(NBAS(ISYMBET),1)
           NTOTC  = MAX(NVIR(ISYMC),1)

           Call DGEMM('N','N',NBAS(ISYMALP),NRHF(ISYMK),NBAS(ISYMBET),
     &                ONE,  WORK(KOFFIC),NTOTAL,XLAMDHB(KOFFLB),NTOTBE,
     &                ZERO, WORK(KHALF),NTOTAL)

           Call DGEMM('T','N',NVIR(ISYMC), NRHF(ISYMK), NBAS(ISYMALP),
     &                ONE, XLAMDPA(KOFFLA),NTOTAL,WORK(KHALF),NTOTAL,
     &                ONE, WORK(KMO),NTOTC )
         END IF

*---------------------------------------------------------------------*
* calculate the D contribution: 
*  FACTOR * sum_alpha I(alpha beta;j^C del) * CTR2(alpha i, j; del)^AB 
*---------------------------------------------------------------------*
         ISYMAI = MULD2H(ISYMJ,ISYCTAB) 
         ISYMI  = MULD2H(ISYMALP,ISYMAI)

         KOFFSV = KDRES + IT1AO(ISYMBET,ISYMI) 
         KOFFZB = KCTR2AB + IT2BGD(ISYMAI,ISYMJ) 
     &                    + NT1AO(ISYMAI)*(J-1) + IT1AO(ISYMALP,ISYMI)
         KOFFIC = KSQC + IAODIS(ISYMALP,ISYMBET)

         NTOTBE = MAX(NBAS(ISYMBET),1)

         NTOTAL = MAX(NBAS(ISYMALP),1)
         Call DGEMM('T','N',NBAS(ISYMBET),NRHF(ISYMI),NBAS(ISYMALP),
     &              FACTOR, WORK(KOFFIC),NTOTAL,WORK(KOFFZB),NTOTAL,
     &              ONE, WORK(KOFFSV),NTOTBE )
           
*---------------------------------------------------------------------*
* end loop over blocks of unpacked AO integrals
*---------------------------------------------------------------------*
        END DO ! ISYMBET

*---------------------------------------------------------------------*
* calculate the C contribution: -sum_ck I(c k;j del) * CTR2(c k;a j) 
*---------------------------------------------------------------------*
        ISYMAJ = MULD2H(ISYMCK,ISYCTR)
        ISYMA  = MULD2H(ISYMJ,ISYMAJ)

        IF (MULD2H(ISYMA,ISYDEL) .NE. ISYRES) THEN
           CALL QUIT('Symmetry mismatch in CCG_21CD.')
        END IF

        FACT = ONE * FACTOR
        IF (LSAME) FACT = TWO * FACTOR

        DO A = 1, NVIR(ISYMA)
           NAJ   = IT1AM(ISYMA,ISYMJ)   + NVIR(ISYMA)*(J-1)    + A
           KCTR2 = IT2SQ(ISYMCK,ISYMAJ) + NT1AM(ISYMCK)*(NAJ-1) + 1
           KOFFR = KCRES-1 + IVIR(ISYMA)-NRHFT + A
 
           WORK(KOFFR) = WORK(KOFFR) -
     &        FACT * DDOT(NT1AM(ISYMCK),WORK(KMOINT),1,CTR2(KCTR2),1)
        END DO

*---------------------------------------------------------------------*
* end loop over symmetry and index of J
*---------------------------------------------------------------------*
      END DO ! J
      END DO ! ISYMJ

*---------------------------------------------------------------------*
* scale C contribution with XLAMDP0
*---------------------------------------------------------------------*
      DO ISYMI = 1, NSYM
        ISYMA = MULD2H(ISYMI,ISYRES)
        ISYMD = ISYMI

      DO I = 1, NRHF(ISYMI)
        KOFFR = KCRES + IVIR(ISYMA)-NRHFT
        KRESV = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + 1
        KLAM0 = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I-1) 
     &            + IDEL-IBAS(ISYDEL)
        
        Call DAXPY(NVIR(ISYMA),XLAMDP0(KLAM0),WORK(KOFFR),1,
     &             RESVEC(KRESV),1)

      END DO
      END DO

*---------------------------------------------------------------------*
* transform D contribution to MO respresentation using XLAMDH0
*---------------------------------------------------------------------*
      DO ISYMA = 1, NSYM
        ISYMI  = MULD2H(ISYMA,ISYRES)
        ISYMBET = ISYMA

        KOFFSV = KDRES + IT1AO(ISYMBET,ISYMI) 
        KLAM0  = IGLMVI(ISYMBET,ISYMA) + 1
        KRESV  = IT1AM(ISYMA,ISYMI) + 1

        NTOTBE = MAX(NBAS(ISYMBET),1)
        NTOTA  = MAX(NVIR(ISYMA),1)

        Call DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMBET),
     &             ONE, XLAMDH0(KLAM0),NTOTBE, WORK(KOFFSV),NTOTBE,
     &             ONE, RESVEC(KRESV),NTOTA)

      END DO

*---------------------------------------------------------------------*
* that's it! return.
*---------------------------------------------------------------------*

      CALL QEXIT('CCG_21CD')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_22CD */
*=====================================================================*
      SUBROUTINE CCG_22CD(DELTA2,CDTERM,ISYRES,
     &                    XBIAJK,EBIAJK,ISYMXB,
     &                    ZBIAJK,CBIAJK,ISYMZB,         
     &                    XCIAJK,ECIAJK,ISYMXC,
     &                    ZCIAJK,CCIAJK,ISYMZC,         
     &                    LSAME, WORK, LWORK)
*---------------------------------------------------------------------*
*
*  Purpose: to calculate the contribution to the G matrix which
*           are analog to the 22C/D contribution to the left transf.
*
*  symmetries & variables:   
*
*            ISYRES  : result vector DELTA2, scratch vector CDTERM
*            ISYMXB  : one-index transf. integrals XBIAJK, EBIAJK
*            ISYMZB  : one-index transf. multipliers ZBIAJK, CBIAJK
*            ISYMXC  : one-index transf. integrals XCIAJK, ECIAJK
*            ISYMZC  : one-index transf. multipliers ZCIAJK, CCIAJK
*
*  if LSAME = .TRUE. C and B perturbations are assumed to be identical
*
*  Christof Haettig, September 1998
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
      LOGICAL LSAME
      INTEGER ISYRES, ISYMXB, ISYMZB, ISYMXC, ISYMZC, LWORK

#if defined (SYS_CRAY)
      REAL ZERO, HALF, ONE, TWO, THREE, FACT, DDOT
      REAL DELTA2(*), CDTERM(*), WORK(LWORK)
      REAL XBIAJK(*),EBIAJK(*),ZBIAJK(*),CBIAJK(*)
      REAL XCIAJK(*),ECIAJK(*),ZCIAJK(*),CCIAJK(*)
#else
      DOUBLE PRECISION DELTA2(*), CDTERM(*), WORK(LWORK)
      DOUBLE PRECISION XBIAJK(*),EBIAJK(*),ZBIAJK(*),CBIAJK(*)         
      DOUBLE PRECISION XCIAJK(*),ECIAJK(*),ZCIAJK(*),CCIAJK(*)         
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, THREE, FACT, DDOT
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0)
      PARAMETER (THREE = 3.0d0, HALF = 0.5d0)

      INTEGER ITAIKJ(8,8), NTAIKJ(8), ISYMAIK, ISYM, ISYMJ, ICOUNT
      INTEGER ISYMAJ, ISYMBI, ISYMAI, ISYMBJ, ISYMKL
      INTEGER NTOTAJ, NTOTBI, NBJ, NJ, NAI, NAIBJ
      INTEGER KSCR, KOFF1, KOFF2, KOFF3

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_22CD')

* check symmetries:
      IF (MULD2H(ISYMXB,ISYMZC).NE.ISYRES .OR.
     &    MULD2H(ISYMXC,ISYMZB).NE.ISYRES      ) THEN
        CALL QUIT('Symmetry mismatch in CCG_22CD.')
      END IF

      DO  ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO 
        NTAIKJ(ISYM) = ICOUNT
      END DO

      KSCR  = 1

      IF (LWORK .LT. MAX(NTAIKJ(ISYMXB),NTAIKJ(ISYMXC)) ) THEN
        CALL QUIT('Insufficient memory in CCG_22CD.')
      END IF

C     CALL DZERO(CDTERM,NT2SQ(ISYRES))

      FACT = HALF
      IF (LSAME) FACT = ONE

C     ----------------------------------------------------------
C     calculate 2 x Cou - Exc of one-index transformed integrals
C     and resort to L(ib,kl) storage for DGEMM
C     ----------------------------------------------------------
      CALL CCG_SORT1(XBIAJK,WORK(KSCR),ISYMXB,0)
      CALL DSCAL(NTAIKJ(ISYMXB),-ONE,WORK(KSCR),1)
      CALL DAXPY(NTAIKJ(ISYMXB),TWO,XBIAJK,1,WORK(KSCR),1)
      CALL CCG_SORT1(WORK(KSCR),EBIAJK,ISYMXB,1)

C     ----------------------------------------------------
C     resort one-index transformed Lagrangian multipliers:
C     ----------------------------------------------------
      CALL CCG_SORT1(ZCIAJK,CCIAJK,ISYMZC,2)

C     ---------------------------------------------
C     contract L^C with I^B to D term contribution:
C     ---------------------------------------------
      DO ISYMBI = 1, NSYM
         ISYMAJ = MULD2H(ISYRES,ISYMBI)
         ISYMKL = MULD2H(ISYMXB,ISYMBI)

         KOFF1 = ISAIKL(ISYMAJ,ISYMKL) + 1
         KOFF2 = ISAIKL(ISYMBI,ISYMKL) + 1
         KOFF3 = IT2SQ(ISYMAJ,ISYMBI)  + 1

         NTOTAJ = MAX(NT1AM(ISYMAJ),1)
         NTOTBI = MAX(NT1AM(ISYMBI),1)

         CALL DGEMM('N','T',NT1AM(ISYMAJ),NT1AM(ISYMBI),NMATIJ(ISYMKL),
     &              -FACT,CCIAJK(KOFF1),NTOTAJ,EBIAJK(KOFF2),NTOTBI,
     &               ZERO,CDTERM(KOFF3),NTOTAJ)
      END DO

C     WRITE (LUPRI,*) 'NORM^2(EBIAJK):',DDOT(NTAIKJ(ISYMXB),EBIAJK,1,EBIAJK,1)
C     WRITE (LUPRI,*) 'NORM^2(CCIAJK):',DDOT(NTAIKJ(ISYMZC),CCIAJK,1,CCIAJK,1)

C     WRITE (LUPRI,*) 'CDTERM after 1. contribution:'
C     CALL CC_PRSQ(WORK,CDTERM,ISYRES,0,1)

      IF (.NOT. LSAME) THEN

C       ----------------------------------------------------------
C       calculate 2 x Cou - Exc of one-index transformed integrals
C       and resort to L(ib,kl) storage for DGEMM
C       ----------------------------------------------------------
        CALL CCG_SORT1(XCIAJK,WORK(KSCR),ISYMXC,0)
        CALL DSCAL(NTAIKJ(ISYMXC),-ONE,WORK(KSCR),1)
        CALL DAXPY(NTAIKJ(ISYMXC),TWO,XCIAJK,1,WORK(KSCR),1)
        CALL CCG_SORT1(WORK(KSCR),ECIAJK,ISYMXC,1)

C       ----------------------------------------------------
C       resort one-index transformed Lagrangian multipliers:
C       ----------------------------------------------------
        CALL CCG_SORT1(ZBIAJK,CBIAJK,ISYMZB,2)

C     WRITE (LUPRI,*) 'NORM^2(ECIAJK):',DDOT(NTAIKJ(ISYMXC),ECIAJK,1,ECIAJK,1)
C     WRITE (LUPRI,*) 'NORM^2(CBIAJK):',DDOT(NTAIKJ(ISYMZB),CBIAJK,1,CBIAJK,1)

C       ---------------------------------------------
C       contract L^B with I^C to D term contribution:
C       ---------------------------------------------
        DO ISYMBI = 1, NSYM
         ISYMAJ = MULD2H(ISYRES,ISYMBI)
         ISYMKL = MULD2H(ISYMXC,ISYMBI)

         KOFF1 = ISAIKL(ISYMAJ,ISYMKL) + 1
         KOFF2 = ISAIKL(ISYMBI,ISYMKL) + 1
         KOFF3 = IT2SQ(ISYMAJ,ISYMBI)  + 1

         NTOTAJ = MAX(NT1AM(ISYMAJ),1)
         NTOTBI = MAX(NT1AM(ISYMBI),1)

         CALL DGEMM('N','T',NT1AM(ISYMAJ),NT1AM(ISYMBI),NMATIJ(ISYMKL),
     &              -FACT,CBIAJK(KOFF1),NTOTAJ,ECIAJK(KOFF2),NTOTBI,
     &                ONE,CDTERM(KOFF3),NTOTAJ)
        END DO

      END IF

C     WRITE (LUPRI,*) 'CDTERM after 2. contribution:'
C     CALL CC_PRSQ(WORK,CDTERM,ISYRES,0,1)

C     --------------------------------------
C     apply (2 - Pij) to the D contribution:
C     --------------------------------------
      CALL CCRHS_T2TR(CDTERM,WORK,LWORK,ISYRES)

C     ------------------------------------------------------------
C     transpose occupied indeces before adding the C contribution:
C     ------------------------------------------------------------
      CALL CCSD_T2TP(CDTERM,WORK,LWORK,ISYRES)

C     WRITE (LUPRI,*) 'CDTERM before 3. contribution:'
C     CALL CC_PRSQ(WORK,CDTERM,ISYRES,0,1)

C     WRITE (LUPRI,*) 'ISYRES:',ISYRES
C     WRITE (LUPRI,*) 'ISYMXB,ISYMXC:',ISYMXB,ISYMXC
C     WRITE (LUPRI,*) 'ISYMZB,ISYMZC:',ISYMZB,ISYMZC

C     WRITE (LUPRI,*) 'NORM^2(XBIAJK):',DDOT(NTAIKJ(ISYMXB),XBIAJK,1,XBIAJK,1)
C     WRITE (LUPRI,*) 'NORM^2(ZCIAJK):',DDOT(NTAIKJ(ISYMZC),ZCIAJK,1,ZCIAJK,1)

C     --------------------------------------------------
C     resort one-index transformed integrals to (lb,ik):
C     --------------------------------------------------
      CALL CCG_SORT1(XBIAJK,EBIAJK,ISYMXB,0)
C     WRITE (LUPRI,*) 'NORM^2(EBIAJK):',DDOT(NTAIKJ(ISYMXB),EBIAJK,1,EBIAJK,1)
      CALL CCG_SORT1(EBIAJK,XBIAJK,ISYMXB,1)

C     -------------------------------------------------------------
C     calculate 2 x Exc + Cou of one-index transformed multipliers:
C     and resort to Zeta(ja,kl) storage for DGEMM
C     -------------------------------------------------------------
      CALL CCG_SORT1(ZCIAJK,CCIAJK,ISYMZC,0)
C     WRITE (LUPRI,*) 'NORM^2(CCIAJK):',DDOT(NTAIKJ(ISYMZC),CCIAJK,1,CCIAJK,1)
      CALL DAXPY(NTAIKJ(ISYMZC),TWO,CCIAJK,1,ZCIAJK,1)
C     WRITE (LUPRI,*) 'NORM^2(ZCIAJK):',DDOT(NTAIKJ(ISYMZC),ZCIAJK,1,ZCIAJK,1)
      CALL CCG_SORT1(ZCIAJK,CCIAJK,ISYMZC,2)

C     WRITE (LUPRI,*) 'NORM^2(XBIAJK):',DDOT(NTAIKJ(ISYMXB),XBIAJK,1,XBIAJK,1)
C     WRITE (LUPRI,*) 'NORM^2(CCIAJK):',DDOT(NTAIKJ(ISYMZC),CCIAJK,1,CCIAJK,1)

C     ---------------------------------------------
C     contract L^C with I^B to C term contribution:
C     ---------------------------------------------
      DO ISYMBI = 1, NSYM
         ISYMAJ = MULD2H(ISYRES,ISYMBI)
         ISYMKL = MULD2H(ISYMXB,ISYMBI)

         KOFF1 = ISAIKL(ISYMAJ,ISYMKL) + 1
         KOFF2 = ISAIKL(ISYMBI,ISYMKL) + 1
         KOFF3 = IT2SQ(ISYMAJ,ISYMBI)  + 1

         NTOTAJ = MAX(NT1AM(ISYMAJ),1)
         NTOTBI = MAX(NT1AM(ISYMBI),1)

         CALL DGEMM('N','T',NT1AM(ISYMAJ),NT1AM(ISYMBI),NMATIJ(ISYMKL),
     &               FACT,CCIAJK(KOFF1),NTOTAJ,XBIAJK(KOFF2),NTOTBI,
     &                ONE,CDTERM(KOFF3),NTOTAJ)
      END DO

C     WRITE (LUPRI,*) 'CDTERM after 3. contribution:'
C     CALL CC_PRSQ(WORK,CDTERM,ISYRES,0,1)

      IF (.NOT. LSAME) THEN

C       --------------------------------------------------
C       resort one-index transformed integrals to (lb,ik):
C       --------------------------------------------------
        CALL CCG_SORT1(XCIAJK,ECIAJK,ISYMXC,0)
        CALL CCG_SORT1(ECIAJK,XCIAJK,ISYMXC,1)

C       -------------------------------------------------------------
C       calculate 2 x Exc + Cou of one-index transformed multipliers:
C       and resort to Zeta(ja,kl) storage for DGEMM
C       -------------------------------------------------------------
        CALL CCG_SORT1(ZBIAJK,CBIAJK,ISYMZB,0)
        CALL DAXPY(NTAIKJ(ISYMZB),TWO,CBIAJK,1,ZBIAJK,1)
        CALL CCG_SORT1(ZBIAJK,CBIAJK,ISYMZB,2)

C       ---------------------------------------------
C       contract L^B with I^C to C term contribution:
C       ---------------------------------------------
        DO ISYMBI = 1, NSYM
         ISYMAJ = MULD2H(ISYRES,ISYMBI)
         ISYMKL = MULD2H(ISYMXC,ISYMBI)

         KOFF1 = ISAIKL(ISYMAJ,ISYMKL) + 1
         KOFF2 = ISAIKL(ISYMBI,ISYMKL) + 1
         KOFF3 = IT2SQ(ISYMAJ,ISYMBI)  + 1

         NTOTAJ = MAX(NT1AM(ISYMAJ),1)
         NTOTBI = MAX(NT1AM(ISYMBI),1)

         CALL DGEMM('N','T',NT1AM(ISYMAJ),NT1AM(ISYMBI),NMATIJ(ISYMKL),
     &               FACT,CBIAJK(KOFF1),NTOTAJ,XCIAJK(KOFF2),NTOTBI,
     &                ONE,CDTERM(KOFF3),NTOTAJ)
        END DO

      END IF

C     WRITE (LUPRI,*) 'CDTERM after 4. contribution:'
C     CALL CC_PRSQ(WORK,CDTERM,ISYRES,0,1)

C     ---------------------------------------------
C     'back'-transposition of the occupied indeces:
C     ---------------------------------------------
      CALL CCSD_T2TP(CDTERM,WORK,LWORK,ISYRES)

*---------------------------------------------------------------------*
*     storage in result vector:
*---------------------------------------------------------------------*
      DO ISYMBJ = 1, NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYRES)

         IF (ISYMAI .EQ. ISYMBJ) THEN

            DO NBJ = 1, NT1AM(ISYMBJ)
               NJ = IT2SQ(ISYMAI,ISYMBJ)+NT1AM(ISYMAI)*(NBJ-1)
                 
               DO NAI = 1, NT1AM(ISYMAI)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  DELTA2(NAIBJ) = DELTA2(NAIBJ) + CDTERM(NJ + NAI)
               END DO

               NAIBJ = IT2AM(ISYMBJ,ISYMBJ) + INDEX(NBJ,NBJ)
               DELTA2(NAIBJ) = DELTA2(NAIBJ) + CDTERM(NJ + NBJ)

            END DO

         ELSE IF (ISYMAI .LT. ISYMBJ) THEN

            NJ    = IT2SQ(ISYMAI,ISYMBJ) + 1
            NAIBJ = IT2AM(ISYMAI,ISYMBJ) + 1
            CALL DAXPY(NT1AM(ISYMAI)*NT1AM(ISYMBJ),ONE,CDTERM(NJ),1,
     &                                              DELTA2(NAIBJ),1)
 
         ELSE IF (ISYMAI .GT. ISYMBJ) THEN

            DO NBJ = 1, NT1AM(ISYMBJ)
               NJ    = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1) + 1
               NAIBJ = IT2AM(ISYMBJ,ISYMAI) + NBJ
      
               Call DAXPY(NT1AM(ISYMAI),ONE,CDTERM(NJ),1,
     &                    DELTA2(NAIBJ),NT1AM(ISYMBJ))
            END DO 

         END IF

      END DO 

      CALL QEXIT('CCG_22CD')

      RETURN 
      END
*---------------------------------------------------------------------*
c/* Deck CCG_4O */
*=====================================================================*
      SUBROUTINE CCG_4O(X4O, ISYM4O, XOVOV, ISYOVOV, T1AMPB, ISYMTB, 
     &                  T1AMPC, ISYMTC, WORK, LWORK, IOPT)
*---------------------------------------------------------------------*
*
*   Purpose: transform (ia|jb) integrals with B1 and C1 amplitudes to 
*
*            X4O(ik,jl) = (i k^B|j l^C) + (i k^C|j l^B)
*                       =  sum(ab) (ia|jb) * T1AMPB(a k) * T1AMPC(b l) 
*                         +sum(ab) (ia|jb) * T1AMPC(a k) * T1AMPB(b l) 
*
*   IOPT = 1  :  packed integrals (jb|id), output in NGAMMA format
*   IOPT = 2  :  full square of integrals, output in NGAMMA format
*   IOPT = 3  :  packed integrals (jb|id), output in N3ORHF format
*   IOPT = 4  :  full square of integrals, output in N3ORHF format
* 
*   needed for CCSD part of CCQR_G (22A/22B contributions)
*   and for H matrix (CC_HMAT)
*
*   Christof Haettig, August 1996, changes for H matrix February 1998
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
      CHARACTER MSGDBG*(16)
      PARAMETER (MSGDBG='[debug] CCG_4O> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM4O, ISYOVOV, ISYMTB, ISYMTC, IOPT
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL X4O(*)
      REAL XOVOV(*)
      REAL T1AMPB(*)      ! dimension (NT1AM(ISYMTB))
      REAL T1AMPC(*)      ! dimension (NT1AM(ISYMTC))
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION X4O(*)
      DOUBLE PRECISION XOVOV(*)
      DOUBLE PRECISION T1AMPB(*)      ! dimension (NT1AM(ISYMTB))
      DOUBLE PRECISION T1AMPC(*)      ! dimension (NT1AM(ISYMTC))
      DOUBLE PRECISION WORK(LWORK)
#endif

      INTEGER KEND0, ISYMB, ISYMIAJ, ISYCIKJ, ISYBIKJ, KXCIKJ, KXBIKJ
      INTEGER KEND1, LEND1, ISYTCTB, IOPT2

C      INTEGER INDEX
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_4O')

* begin:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'just entered this routine...'
      END IF

* check symmetries:
      ISYTCTB = MULD2H(ISYMTC,ISYMTB)
      IF (ISYM4O .NE. MULD2H(ISYTCTB,ISYOVOV)) THEN
        CALL QUIT('Symmetry mismatch in CCG_4O.')
      END IF

* check IOPT and initialize result array:
      If (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
        CALL DZERO( X4O, NGAMMA(ISYM4O) )
        IOPT2 = IOPT
      ELSE IF (IOPT.EQ.3 .OR. IOPT.EQ.4) THEN
        CALL DZERO( X4O, N3ORHF(ISYM4O) )
        IOPT2 = IOPT - 2
      ELSE
        CALL QUIT('Illegal value for IOPT in CCG_4O.') 
      ENDIF

* set start of work space:
      KEND0 = 1


* start outer loop over symmetry blocks:
      DO ISYMB = 1, NSYM

*---------------------------------------------------------------------*
*       memory allocation:
*---------------------------------------------------------------------*
        ISYMIAJ = MULD2H(ISYMB,ISYOVOV)
        ISYCIKJ = MULD2H(ISYMIAJ,ISYMTC)
        ISYBIKJ = MULD2H(ISYMIAJ,ISYMTB)
         
        KXCIKJ = KEND0
        KXBIKJ = KEND0
        KEND1  = KEND0 + MAX( NMAIJK(ISYCIKJ) , NMAIJK(ISYBIKJ) )
        LEND1  = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCG_4O.')
        END IF
        
        DO B = 1, NVIR(ISYMB)

*---------------------------------------------------------------------*
* calculate one-index transformed (i k^C|j b) distributions
* and update X4O(ik,jl) by (i k^C|j b) * B1(b l) for all i,k,j,l
*---------------------------------------------------------------------*
          Call CCG_TRANS2(WORK(KXCIKJ), ISYCIKJ, XOVOV, ISYOVOV,
     &                    T1AMPC, ISYMTC, B, ISYMB, IOPT2)

          Call CCG_4OS(X4O, ISYM4O, WORK(KXCIKJ), ISYCIKJ,
     &                 T1AMPB, ISYMTB, B, ISYMB, IOPT)
*---------------------------------------------------------------------*
* calculate one-index transformed (i k^B|j b) distributions
* and update X4O(ik,jl) by (i k^B|j b) * C1(b l) for all i,k,j,l
*---------------------------------------------------------------------*
          Call CCG_TRANS2(WORK(KXBIKJ), ISYBIKJ, XOVOV, ISYOVOV,
     &                    T1AMPB, ISYMTB, B, ISYMB, IOPT2)

          Call CCG_4OS(X4O, ISYM4O, WORK(KXBIKJ), ISYBIKJ,
     &                 T1AMPC, ISYMTC, B, ISYMB, IOPT)
        
*---------------------------------------------------------------------*
        END DO ! B
      END DO ! ISYMB

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'leaving this routine...'
      END IF

      CALL QEXIT('CCG_4O')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_4OS */
*=====================================================================*
      SUBROUTINE CCG_4OS(X4O,ISYM4O, XIKJB,ISYIKJ, 
     &                   T1,ISYMT1, B,ISYMB, IOPT)
*---------------------------------------------------------------------*
*
*     Purpose: update I(ik,jl) integrals by:
*
*     X4O(ik,jl) = X4O(ik,jl) + (i k|j b) * T1(b l)   for all i,k,j,l
*
*     called by CCG_4O
*
*     Christof Haettig, August 1996, changes for H matrix: Februar 1998
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CCG_4OS> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM4O, ISYIKJ, ISYMT1, ISYMB, IOPT

#if defined (SYS_CRAY)
      REAL X4O(*)
      REAL XIKJB(*)   ! dimension (NMAIJK(ISYIKJ))
      REAL T1(*)      ! dimension (NT1AM(ISYMT1))
#else
      DOUBLE PRECISION X4O(*)
      DOUBLE PRECISION XIKJB(*)   ! dimension (NMAIJK(ISYIKJ))
      DOUBLE PRECISION T1(*)      ! dimension (NT1AM(ISYMT1))
#endif

      INTEGER ISYML, NBL, ISYMJ, ISYMJL, ISYMIK, NJL, NIK, NIKJ, NIKJL
      INTEGER NLJ, NLJKI, OFFIK, OFFI, ISYLJK, ISYMI, ISYMLJ, ISYMK

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_4OS')

* begin:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'just entered this routine...'
      END IF

* check symmetries:
      IF (MULD2H(ISYM4O,ISYMT1) .NE. MULD2H(ISYIKJ,ISYMB)) THEN
        CALL QUIT('Symmetry mismatch in CCG_4OS.')
      END IF

*=====================================================================*
* store result in NGAMMA format:
*=====================================================================*
      IF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN

*---------------------------------------------------------------------*
* start loop over occupied index l:
*---------------------------------------------------------------------*
      ISYML = MULD2H(ISYMB,ISYMT1)

      DO L = 1, NRHF(ISYML)

        NBL = IT1AM(ISYMB,ISYML) + NVIR(ISYMB)*(L-1) + B

*---------------------------------------------------------------------*
* add contribution to X4O:
*---------------------------------------------------------------------*
        DO ISYMJ = 1, NSYM

          ISYMJL = MULD2H(ISYMJ,ISYML)
          ISYMIK = MULD2H(ISYMJL,ISYM4O)
              
          IF (ISYMIK .LT. ISYMJL) THEN
            DO J = 1, NRHF(ISYMJ)
              NJL = IMATIJ(ISYMJ,ISYML) + NRHF(ISYMJ)*(L-1) + J

              DO NIK = 1, NMATIJ(ISYMIK)
                NIKJ  = IMAIJK(ISYMIK,ISYMJ) 
     &                     + NMATIJ(ISYMIK)*(J-1)   + NIK
                NIKJL = IGAMMA(ISYMIK,ISYMJL)
     &                     + NMATIJ(ISYMIK)*(NJL-1) + NIK

                X4O(NIKJL) = X4O(NIKJL) + XIKJB(NIKJ) * T1(NBL)
              END DO
            END DO

          ELSE IF (ISYMIK .EQ. ISYMJL) THEN

            DO J = 1, NRHF(ISYMJ)
              NJL = IMATIJ(ISYMJ,ISYML) + NRHF(ISYMJ)*(L-1) + J

              DO NIK = 1, NJL
                NIKJ  = IMAIJK(ISYMIK,ISYMJ) 
     &                     + NMATIJ(ISYMIK)*(J-1) + NIK
                NIKJL = IGAMMA(ISYMIK,ISYMJL) + INDEX(NIK,NJL)

                X4O(NIKJL) = X4O(NIKJL) + XIKJB(NIKJ) * T1(NBL)
              END DO
            END DO

          END IF

        END DO ! ISYMJ

      END DO ! L 

*=====================================================================*
* store result in N3ORHF format:
*=====================================================================*
      ELSE IF (IOPT.EQ.3 .OR. IOPT.EQ.4) THEN

      ISYML = MULD2H(ISYMB,ISYMT1)

      DO ISYMJ = 1, NSYM
      DO ISYMI = 1, NSYM

        ISYMLJ = MULD2H(ISYML,ISYMJ)
        ISYMIK = MULD2H(ISYMLJ,ISYM4O)
        ISYMK  = MULD2H(ISYMIK,ISYMI)
        ISYLJK = MULD2H(ISYMLJ,ISYMK)

      DO I = 1, NRHF(ISYMI)
        OFFI  = I3ORHF(ISYLJK,ISYMI) + NMAIJK(ISYLJK)*(I-1)

        DO K = 1, NRHF(ISYMK)

          OFFIK = OFFI + IMAIJK(ISYMLJ,ISYMK) + NMATIJ(ISYMLJ)*(K-1)
          NIK   = IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K-1) + I

          DO J = 1, NRHF(ISYMJ)

            NIKJ  = IMAIJK(ISYMIK,ISYMJ) + NMATIJ(ISYMIK)*(J-1) + NIK

            DO L = 1, NRHF(ISYML)

              NBL   = IT1AM(ISYMB,ISYML) + NVIR(ISYMB)*(L-1) + B
              NLJ   = IMATIJ(ISYML,ISYMJ) + NRHF(ISYML)*(J-1) + L
              NLJKI = OFFIK + NLJ
  
              X4O(NLJKI) = X4O(NLJKI) + XIKJB(NIKJ) * T1(NBL)
  
            END DO ! L 

          END DO ! J

        END DO ! K

      END DO ! I  

      END DO ! ISYMI
      END DO ! ISYMJ

*=====================================================================*
* return:
*=====================================================================*
      ELSE
         CALL QUIT('Illegal value for IOPT in CCG_4OS.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'leaving this routine...'
      END IF

      CALL QEXIT('CCG_4OS')

      RETURN
      END
*---------------------------------------------------------------------*
c/*  Deck CCG_OOVV */
*=====================================================================*
      SUBROUTINE CCG_OOVV(Jint, Kint, ISYJKINT, XOVOV, ISYOVOV, 
     &                    B1AMP, ISYMTB, C1AMP, ISYMTC,
     &                    WORK, LWORK, D, ISYMD, IOPT,IPCK)
*---------------------------------------------------------------------*
*
*     Purpose: transform a three-index batch of (jb|id) integrals
*              (with fixed index d) with B1 and C1 amplitudes to 
*
*              Jint(el,j) = (e^C l^B|j d) + (e^B l^C|j d)
*  and
*              Kint(el,j) = (j l^B|e^C d) + (j l^C|e^B d)
*
*
*              IOPT = 1  : Jint only
*              IOPT = 2  : Kint only
*              IOPT = 3  : Jint & Kint 
*
*              IPCK = 1  :  packed integrals (jb|id)
*              IPCK = 2  :  full square of integrals
* 
*              needed for CCSD part of CCQR_G (CCG_22C contribiution)
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
#if defined (SYS_CRAY)
      REAL ZERO, ONE
#else
      DOUBLE PRECISION ZERO, ONE
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0)

      INTEGER ISYJKINT, ISYOVOV, ISYMTB, ISYMTC, ISYMD
      INTEGER IPCK, LWORK, IOPT

#if defined (SYS_CRAY)
      REAL Kint(*)       ! dimension (NT2BCD(ISYJKINT))
      REAL Jint(*)       ! dimension (NT2BCD(ISYJKINT))
      REAL XOVOV(*)
      REAL B1AMP(*)       ! dimension (NT1AM(ISYMTB))
      REAL C1AMP(*)       ! dimension (NT1AM(ISYMTC))
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION Kint(*)       ! dimension (NT2BCD(ISYJKINT))
      DOUBLE PRECISION Jint(*)       ! dimension (NT2BCD(ISYJKINT))
      DOUBLE PRECISION XOVOV(*)
      DOUBLE PRECISION B1AMP(*)       ! dimension (NT1AM(ISYMTB))
      DOUBLE PRECISION C1AMP(*)       ! dimension (NT1AM(ISYMTC))
      DOUBLE PRECISION WORK(LWORK)
#endif

      INTEGER ISYTCTB, ISYMIAJ, ISYCJLI, ISYBJLI, ISYJLE, ISYMEL
      INTEGER KXCJLI, KXBJLI, JXCJLI, JXBJLI, KEND1, LEND1, ISYML
      INTEGER ISYMJL, ISYME, KKSCR, KJSCR, KEND2, LEND2, ISYMI, NELJ
      INTEGER KOFFK, KOFFJ, KT1, NTOTJL, NTOTE, ISYMJ, NJL, NJLE, NEL

C      INTEGER INDEX
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_OOVV')

* check symmetries:
      ISYTCTB = MULD2H(ISYMTC,ISYMTB)

      IF (MULD2H(ISYOVOV,ISYTCTB) .NE. MULD2H(ISYJKINT,ISYMD)) THEN
        CALL QUIT('Symmetry mismatch in CCG_OOVV.')
      END IF

*---------------------------------------------------------------------*
* memory allocation
*---------------------------------------------------------------------*
      ISYMIAJ = MULD2H(ISYMD,ISYOVOV)
      ISYCJLI = MULD2H(ISYMIAJ,ISYMTC)
      ISYBJLI = MULD2H(ISYMIAJ,ISYMTB)
      ISYJLE  = MULD2H(ISYMIAJ,ISYTCTB)

      KXCjli = 1
      KXBjli = KXCjli + NMAIJK(ISYCJLI)
      JXCjli = KXBjli + NMAIJK(ISYBJLI)
      JXBjli = JXCjli + NMAIJK(ISYCJLI)
      KEND1  = JXBjli + NMAIJK(ISYBJLI)
      LEND1  = LWORK - KEND1

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCG_OOVV.')
      END IF

*---------------------------------------------------------------------*
* calculate one-index transformed (j l^C|i d) distribution
* and contruct distribution with j and i indeces exchanged
*---------------------------------------------------------------------*
      Call CCG_TRANS2(WORK(KXCjli),ISYCJLI,XOVOV,ISYOVOV,C1AMP,ISYMTC,
     &                  D,ISYMD,IPCK)

      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 3) THEN
        Call CCLT_PI3O(WORK(JXCjli),WORK(KXCjli),ISYCJLI)
      END IF

*---------------------------------------------------------------------*
* calculate one-index transformed (j l^B|i d) distribution
* and contruct distribution with j and i indeces exchanged
*---------------------------------------------------------------------*
      Call CCG_TRANS2(WORK(KXBjli),ISYBJLI,XOVOV,ISYOVOV,B1AMP,ISYMTB,
     &                  D,ISYMD,IPCK)

      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 3) THEN
        Call CCLT_PI3O(WORK(JXBjli),WORK(KXBjli),ISYBJLI)
      END IF

*---------------------------------------------------------------------*
* calculate  I^d(el,j) = (j l^C|e^B d) + (j l^B|e^C d):
*---------------------------------------------------------------------*
      DO ISYMJL = 1, NSYM
        ISYME   = MULD2H(ISYMJL,ISYJLE)

        KKSCR = KEND1
        KJSCR = KKSCR + NMATIJ(ISYMJL) * NVIR(ISYME)
        KEND2 = KJSCR + NMATIJ(ISYMJL) * NVIR(ISYME)
        LEND2 = LWORK - KEND2

        IF (LEND2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCG_22C.')
        END IF

*.....................................................................*
* first contribution Kint: -sum_i (j l^C|i d) * T1B(e i)
*.....................................................................*
        ISYMI = MULD2H(ISYMJL,ISYCJLI)

        KOFFK = KXCjli + IMAIJK(ISYMJL,ISYMI) 
        KOFFJ = JXCjli + IMAIJK(ISYMJL,ISYMI) 
        KT1   = IT1AM(ISYME,ISYMI) + 1

        NTOTJL = MAX(NMATIJ(ISYMJL),1)
        NTOTE  = MAX(NVIR(ISYME),1)

        IF (IOPT .EQ. 2 .OR. IOPT .EQ. 3) THEN
          Call DGEMM('N','T',NMATIJ(ISYMJL),NVIR(ISYME),NRHF(ISYMI),
     &               -ONE, WORK(KOFFK), NTOTJL, B1AMP(KT1), NTOTE,
     &               ZERO, WORK(KKSCR), NTOTJL )
        END IF

        IF (IOPT .EQ. 1 .OR. IOPT .EQ. 3) THEN
*.....................................................................*
* analogously for Jint: -sum_i (i l^C|j d) * T1B(e i)
*.....................................................................*
          Call DGEMM('N','T',NMATIJ(ISYMJL),NVIR(ISYME),NRHF(ISYMI),
     &               -ONE, WORK(KOFFJ), NTOTJL, B1AMP(KT1), NTOTE,
     &               ZERO, WORK(KJSCR), NTOTJL )
        END IF


*.....................................................................*
* add second contribution: -sum_i (j l^C|i d) * T1C(e i)
*.....................................................................*
        ISYMI = MULD2H(ISYMJL,ISYBJLI)

        KOFFK = KXBjli + IMAIJK(ISYMJL,ISYMI)
        KOFFJ = JXBjli + IMAIJK(ISYMJL,ISYMI)
        KT1   = IT1AM(ISYME,ISYMI) + 1

        NTOTJL = MAX(NMATIJ(ISYMJL),1)
        NTOTE  = MAX(NVIR(ISYME),1)

        IF (IOPT .EQ. 2 .OR. IOPT .EQ. 3) THEN
          Call DGEMM('N','T',NMATIJ(ISYMJL),NVIR(ISYME),NRHF(ISYMI),
     &               -ONE, WORK(KOFFK), NTOTJL, C1AMP(KT1), NTOTE,
     &               +ONE, WORK(KKSCR), NTOTJL )
        END IF

        IF (IOPT .EQ. 1 .OR. IOPT .EQ. 3) THEN
*.....................................................................*
* analogously for Jint: -sum_i (i l^C|j d) * T1C(e i)
*.....................................................................*
          Call DGEMM('N','T',NMATIJ(ISYMJL),NVIR(ISYME),NRHF(ISYMI),
     &               -ONE, WORK(KOFFJ), NTOTJL, C1AMP(KT1), NTOTE,
     &               +ONE, WORK(KJSCR), NTOTJL )
        END IF


*---------------------------------------------------------------------*
* resort result from I(j l;e d) to I(e l;j d):
*---------------------------------------------------------------------*
        DO ISYMJ = 1, NSYM
          ISYML  = MULD2H(ISYMJ,ISYMJL)
          ISYMEL =  MULD2H(ISYME,ISYML)

          IF (IOPT .EQ. 2 .OR. IOPT .EQ. 3) THEN
            DO J = 1, NRHF(ISYMJ)
            DO L = 1, NRHF(ISYML)
            DO E = 1, NVIR(ISYME)
              NJL  = IMATIJ(ISYMJ,ISYML) + NRHF(ISYMJ)*(L-1) + J
              NJLE = NMATIJ(ISYMJL)*(E-1) + NJL

              NEL  = IT1AM(ISYME,ISYML) + NVIR(ISYME)*(L-1) + E
              NELJ = IT2BCD(ISYMEL,ISYMJ) + NT1AM(ISYMEL)*(J-1) + NEL

              Kint(NELJ) = WORK(KKSCR-1 + NJLE)

            END DO
            END DO
            END DO
          END IF

          IF (IOPT .EQ. 1 .OR. IOPT .EQ. 3) THEN
            DO J = 1, NRHF(ISYMJ)
            DO L = 1, NRHF(ISYML)
            DO E = 1, NVIR(ISYME)
              NJL  = IMATIJ(ISYMJ,ISYML) + NRHF(ISYMJ)*(L-1) + J
              NJLE = NMATIJ(ISYMJL)*(E-1) + NJL

              NEL  = IT1AM(ISYME,ISYML) + NVIR(ISYME)*(L-1) + E
              NELJ = IT2BCD(ISYMEL,ISYMJ) + NT1AM(ISYMEL)*(J-1) + NEL

              Jint(NELJ) = WORK(KJSCR-1 + NJLE)

            END DO
            END DO
            END DO
          END IF

        END DO ! ISYMJ
          
      END DO ! ISYMJL

      CALL QEXIT('CCG_OOVV')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_RDIAJB */
*=====================================================================*
      SUBROUTINE CCG_RDIAJB(XIAJB,LIAJB)
*---------------------------------------------------------------------*
*     Purpose: read packed (ia|jb) integrals into memory
*              integrals read from unit 52, which is assumed open
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "iratdef.h"
#include "ccinftap.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LIAJB
      
#if defined (SYS_CRAY)
      REAL XIAJB(LIAJB)
#else
      DOUBLE PRECISION XIAJB(LIAJB)
#endif

      CALL QENTER('CCG_RDIAJB')

* check dimension:
      IF (LIAJB .LT. NT2AM(ISYMOP)) THEN
        CALL QUIT(
     *       'Insufficient memory for (ia|jb) integrals in CCG_RDIAJB.')
      END IF 

* rewind integral file:
      REWIND(LUIAJB)

* read (ia|jb) integrals:
      CALL READI(LUIAJB,IRAT*NT2AM(1),XIAJB)

      IF (LOCDBG) THEN
         CALL AROUND( 'CCG_RDIAJB:  Integrals (ia|jb) ' )
         CALL CC_PRP(0.00D0,XIAJB,ISYMOP,0,1)
      ENDIF

      CALL QEXIT('CCG_RDIAJB')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_TI */
*=====================================================================*
      SUBROUTINE CCG_TI(Tint,ISYMTI,T2amp,ISYMTAM,C,ISYMC)
*---------------------------------------------------------------------*
*
*     Purpose: to get a three-index submatrix T^{c}(bk,j) out of the
*              packed matrix of T(dk,cj) amplitudes
*
*              T2amp -- packed matrix of T2 amplitudes
*              Tint  -- three-index submatrix of T2 with fixed c
*
*              needed for CCSD part of CCQR_G
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
 
#if defined (SYS_CRAY)
      REAL T2amp(*)       ! dimension (NT2AM(ISYMTAM))
      REAL Tint(*)       ! dimension (NT2BCD(ISYMTI))
#else
      DOUBLE PRECISION T2amp(*)       ! dimension (NT2AM(ISYMTAM))
      DOUBLE PRECISION Tint(*)       ! dimension (NT2BCD(ISYMTI))
#endif

      INTEGER ISYMTI, ISYMTAM, ISYMC, ISYMDK, ISYMCJ, ISYMJ
      INTEGER NCJ, NDK, NDKCJ, NDKJ

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_TI')

* check symmetries:
      IF (ISYMTI .NE. MULD2H(ISYMTAM,ISYMC)) THEN
        CALL QUIT('SYMMETRY MISMATCH IN CCG_TI.')
      END IF


      DO ISYMDK = 1, NSYM
        ISYMCJ = MULD2H(ISYMDK,ISYMTAM)
        ISYMJ  = MULD2H(ISYMCJ,ISYMC)

      DO J = 1, NRHF(ISYMJ)

        NCJ = IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J-1) + C
 
        IF (ISYMDK. EQ. ISYMCJ) THEN

          DO NDK = 1, NT1AM(ISYMDK)
            NDKCJ = IT2AM(ISYMDK,ISYMCJ) + INDEX(NDK,NCJ)
            NDKJ  = IT2BCD(ISYMDK,ISYMJ) + NT1AM(ISYMDK)*(J-1) + NDK
            Tint(NDKJ) = T2amp(NDKCJ)
          END DO

        ELSE IF (ISYMDK. LT. ISYMCJ) THEN

          NDKCJ = IT2AM(ISYMDK,ISYMCJ) + NT1AM(ISYMDK)*(NCJ-1) + 1
          NDKJ  = IT2BCD(ISYMDK,ISYMJ) + NT1AM(ISYMDK)*(J-1)   + 1
          Call DCOPY(NT1AM(ISYMDK),T2amp(NDKCJ),1,Tint(NDKJ),1)

        ELSE IF (ISYMDK. GT. ISYMCJ) THEN

          NDKCJ = IT2AM(ISYMCJ,ISYMDK) + NCJ
          NDKJ  = IT2BCD(ISYMDK,ISYMJ) + NT1AM(ISYMDK)*(J-1) + 1
          Call DCOPY(NT1AM(ISYMDK),T2amp(NDKCJ),NT1AM(ISYMCJ),
     &                             Tint(NDKJ),1)

        END IF

      END DO ! J
      END DO ! ISYMDK

      CALL QEXIT('CCG_TI')
C
      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_TRANS13 */
*=====================================================================*
      SUBROUTINE CCG_TRANS13(Xiaec,ISYMIAE,XOVOV,ISYOVOV,
     &                       B1AMP,ISYMB1,C,ISYMC,IPOS,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: transform the first (IPOS=1) or the third (IPOS=3) index
*              of a three-index batch of (i a|k c) integrals (with 
*              fixed index c) with the single excitation amplitudes 
*              stored in B1AMP. 
*
*     NOTE the minus sign in the transformation!
*
*
*     IPOS = 1 :    (e^B a|i c) = sum_k -B1(ek) * (k a|i c)
*
*     IPOS = 3 :    (i a|e^B c) = sum_k -B1(ek) * (i a|k c)
*
*     the result is in both cases returned in Xiaec stored as 
*     (a i|e^B;c) using the dimension/adress arrays NCKATR & ICKATR
*
*     (needed for CCSD part of CCQR_G)
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
      CHARACTER MSGDBG*(21)
      PARAMETER (MSGDBG='[debug] CCG_TRANS13> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMIAE, ISYOVOV, ISYMB1, ISYMC, LWORK, IPOS

#if defined (SYS_CRAY)
      REAL XIAEC(*)       ! dimension (NCKATR(ISYMIAE))
      REAL XOVOV(*)       ! dimension (NT2AM(ISYOVOV))
      REAL B1AMP(*)       ! dimension (NT1AM(ISYMB1))
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XIAEC(*)       ! dimension (NCKATR(ISYMIAE))
      DOUBLE PRECISION XOVOV(*)       ! dimension (NT2AM(ISYOVOV))
      DOUBLE PRECISION B1AMP(*)       ! dimension (NT1AM(ISYMB1))
      DOUBLE PRECISION WORK(LWORK)
#endif

      INTEGER ISYMAIK, KAIK, KFREE, LFREE, ISYMAI, ISYMK, ISYME
      INTEGER NEK, NAI, NAIE, NAIK

C      INTEGER INDEX
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_TRANS13')

* check symmetries and IPOS option:
      IF (MULD2H(ISYMIAE,ISYOVOV) .NE. MULD2H(ISYMB1,ISYMC)) THEN
        CALL QUIT('Symmetry mismatch in CCG_TRANS13.')
      END IF
      IF (IPOS .NE. 1 .AND. IPOS .NE. 3) THEN
        CALL QUIT('CCG_TRANS13 called with illegal value for IPOS.')
      END IF

* initialize result:
      Call DZERO (Xiaec, NCKATR(ISYMIAE))

* get submatrix I(ai,k) of integrals:
      ISYMAIK = MULD2H(ISYMC,ISYOVOV)
      KAIK    = 1
      KFREE   = KAIK  + NT2BCD(ISYMAIK)
      LFREE   = LWORK - NT2BCD(ISYMAIK)
      
      IF (LFREE .lt. 0) THEN
        CALL QUIT('Insufficient work space in CCG_TRANS13.')
      END IF

      Call CCG_TI (WORK(KAIK),ISYMAIK,XOVOV,ISYOVOV,C,ISYMC)
      
* if the first index should be transformed, exchange i and k:
      IF (IPOS .eq. 1) THEN
        Call CCLT_P21I(WORK(KAIK),ISYMAIK,WORK(KFREE),LFREE,
     &                 IT2BCD, NT2BCD, IT1AM, NT1AM, NVIR)
      END IF

* do the transformation:
      DO ISYMAI = 1, NSYM
        ISYMK = MULD2H(ISYMAI,ISYMAIK)
        ISYME = MULD2H(ISYMK,ISYMB1)

        DO K = 1, NRHF(ISYMK)
        DO E = 1, NVIR(ISYME)
          NEK  = IT1AM(ISYME,ISYMK) + NVIR(ISYME)*(K-1) + E
        DO NAI = 1, NT1AM(ISYMAI)
          NAIE = ICKATR(ISYMAI,ISYME) + NT1AM(ISYMAI)*(E-1) + NAI
          NAIK = IT2BCD(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1) + NAI
          Xiaec(NAIE) = Xiaec(NAIE) - B1AMP(NEK) * WORK(KAIK-1+NAIK)
        END DO
        END DO
        END DO

        IF (LOCDBG) THEN
          DO E = 1, NVIR(ISYME)
            WRITE (LUPRI,*) MSGDBG, 'ISYMAI,E:',ISYMAI,E
            NAIE = ICKATR(ISYMAI,ISYME) + NT1AM(ISYMAI)*(E-1) + 1
            Call CC_PRP(Xiaec(NAIE),WORK,ISYMAI,1,0)
          END DO
        END IF
     
      END DO


      CALL QEXIT('CCG_TRANS13')
C
      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_TRANS2 */
*=====================================================================*
      SUBROUTINE CCG_TRANS2(Xjlic,ISYMJLI,XOVOV,ISYOVOV,B1AMP,ISYMB1,
     &                      C, ISYMC, IOPT)
*---------------------------------------------------------------------*
*
*     Purpose: transform a three-index batch of (jd|ic) integrals
*              (with fixed index c) with B1 amplitudes to (j l^B|i c)
*              (one-index transformed integrals with second index
*              transformed)
*
*              I^c(jl;i) = sum_a  (j a|i c) * B1AMP(a l)
*
*         IOPT = 1 : XOVOV contains (ia|jb) integrals packed
*                2 : full square matrix of (ia|jb) in XOVOV
*                    --> transform virtual index of the leading pair
*
*              needed for CCSD part of CCQR_G
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
      INTEGER ISYMJLI, ISYMB1, ISYOVOV, ISYMC, IOPT

#if defined (SYS_CRAY)
      REAL XJLIC(*)       ! dimension (NMAIJK(ISYMJLI))
      REAL XOVOV(*)
      REAL B1AMP(*)       ! dimension (NT1AM(ISYMB1))
#else
      DOUBLE PRECISION XJLIC(*)       ! dimension (NMAIJK(ISYMJLI))
      DOUBLE PRECISION XOVOV(*)
      DOUBLE PRECISION B1AMP(*)       ! dimension (NT1AM(ISYMB1))
#endif

      INTEGER ISYMDJ, ISYMCI, ISYMJL, ISYMI, ISYMJ, ISYML, ISYMD
      INTEGER NCI, NDJ, NDL, NJL, NJLI, NDJCI

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_TRANS2')

* check symmetries:
      IF (MULD2H(ISYMJLI,ISYMC) .NE. MULD2H(ISYMB1,ISYOVOV)) 
     &  CALL QUIT('SYMMETRY MISMATCH IN CCG_TRANS2.')

* initialize result:
      Call DZERO (Xjlic, NMAIJK(ISYMJLI))

      DO ISYMDJ = 1, NSYM
        ISYMCI = MULD2H(ISYMDJ,ISYOVOV)
        ISYMI  = MULD2H(ISYMCI,ISYMC)

      DO I = 1, NRHF(ISYMI)

        NCI = IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I-1) + C
 
        DO ISYMD = 1, NSYM
          ISYMJ = MULD2H(ISYMD,ISYMDJ)
          ISYML = MULD2H(ISYMD,ISYMB1)
          ISYMJL = MULD2H(ISYMJ,ISYML)

        IF (IOPT .EQ. 1) THEN

          IF (ISYMDJ .EQ. ISYMCI) THEN

            DO L = 1, NRHF(ISYML)
            DO J = 1, NRHF(ISYMJ)
              NJL   = IMATIJ(ISYMJ,ISYML)  + NRHF(ISYMJ)*(L-1)    + J
              NJLI  = IMAIJK(ISYMJL,ISYMI) + NMATIJ(ISYMJL)*(I-1) + NJL
            DO D = 1, NVIR(ISYMD)
              NDJ   = IT1AM(ISYMD,ISYMJ)   + NVIR(ISYMD)*(J-1)    + D
              NDL   = IT1AM(ISYMD,ISYML)   + NVIR(ISYMD)*(L-1)    + D
              NDJCI = IT2AM(ISYMDJ,ISYMCI) + INDEX(NDJ,NCI)
              Xjlic(NJLI) = Xjlic(NJLI) + XOVOV(NDJCI) * B1AMP(NDL)
            END DO
            END DO
            END DO

          ELSE IF (ISYMDJ .LT. ISYMCI) THEN

            DO L = 1, NRHF(ISYML)
            DO J = 1, NRHF(ISYMJ)
              NJL   = IMATIJ(ISYMJ,ISYML)  + NRHF(ISYMJ)*(L-1)    + J
              NJLI  = IMAIJK(ISYMJL,ISYMI) + NMATIJ(ISYMJL)*(I-1) + NJL
            DO D = 1, NVIR(ISYMD)
              NDJ   = IT1AM(ISYMD,ISYMJ)   + NVIR(ISYMD)*(J-1)    + D
              NDL   = IT1AM(ISYMD,ISYML)   + NVIR(ISYMD)*(L-1)    + D
              NDJCI = IT2AM(ISYMDJ,ISYMCI) + NT1AM(ISYMDJ)*(NCI-1)+NDJ
              Xjlic(NJLI) = Xjlic(NJLI) + XOVOV(NDJCI) * B1AMP(NDL)
            END DO
            END DO
            END DO

          ELSE IF (ISYMDJ .GT. ISYMCI) THEN

            DO L = 1, NRHF(ISYML)
            DO J = 1, NRHF(ISYMJ)
              NJL   = IMATIJ(ISYMJ,ISYML)  + NRHF(ISYMJ)*(L-1)    + J
              NJLI  = IMAIJK(ISYMJL,ISYMI) + NMATIJ(ISYMJL)*(I-1) + NJL
            DO D = 1, NVIR(ISYMD)
              NDJ   = IT1AM(ISYMD,ISYMJ)   + NVIR(ISYMD)*(J-1)    + D
              NDL   = IT1AM(ISYMD,ISYML)   + NVIR(ISYMD)*(L-1)    + D
              NDJCI = IT2AM(ISYMCI,ISYMDJ) + NT1AM(ISYMCI)*(NDJ-1)+NCI
              Xjlic(NJLI) = Xjlic(NJLI) + XOVOV(NDJCI) * B1AMP(NDL)
            END DO
            END DO
            END DO

          END IF

        ELSE IF (IOPT .EQ. 2) THEN

          DO L = 1, NRHF(ISYML)
          DO J = 1, NRHF(ISYMJ)
            NJL   = IMATIJ(ISYMJ,ISYML)  + NRHF(ISYMJ)*(L-1)    + J
            NJLI  = IMAIJK(ISYMJL,ISYMI) + NMATIJ(ISYMJL)*(I-1) + NJL
          DO D = 1, NVIR(ISYMD)
            NDJ   = IT1AM(ISYMD,ISYMJ)   + NVIR(ISYMD)*(J-1)    + D
            NDL   = IT1AM(ISYMD,ISYML)   + NVIR(ISYMD)*(L-1)    + D
            NDJCI = IT2SQ(ISYMDJ,ISYMCI) + NT1AM(ISYMDJ)*(NCI-1)+ NDJ
            Xjlic(NJLI) = Xjlic(NJLI) + XOVOV(NDJCI) * B1AMP(NDL)
          END DO
          END DO
          END DO

        ELSE
          CALL QUIT('Illegal value for parameter IOPT in CCG_TRANS2.')
        END IF

        END DO ! ISYMD

      END DO ! J
      END DO ! ISYMDK

      CALL QEXIT('CCG_TRANS2')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_TRANS4 */
*=====================================================================*
      SUBROUTINE CCG_TRANS4(XIAJK,ISYMIAJK,XOVOV,ISYOVOV,T1AMP,ISYMT1)
*---------------------------------------------------------------------*
*
*     Purpose: transform (ia|jb) integrals to (ia|jk) using T1(bk)
*
*              I^j(ai;k) = sum_b  (i a|j b) * T1AMP(b k)
*
*              XOVOV contains (ia|jb) integrals as full matrix
*
*     needed for CCSD part of CCQR_G
*
*     Christof Haettig, August 1998
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
      INTEGER ISYMIAJK, ISYMT1, ISYOVOV

#if defined (SYS_CRAY)
      REAL XIAJK(*), XOVOV(*), T1AMP(*)    
      REAL ONE, ZERO
#else
      DOUBLE PRECISION XIAJK(*), XOVOV(*), T1AMP(*) 
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0)

      INTEGER ITAIKJ(8,8), NTAIKJ(8)
      INTEGER ISYMB, ISYMK, ISYMAI, ISYMBJ, ISYMAIK, ISYMJ
      INTEGER NBJ, KOFF1, KOFF2, KOFF3, ICOUNT, ISYM, NTOTB, NTOTAI

C      INTEGER INDEX
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_TRANS4')

*---------------------------------------------------------------------*
* check symmetries & and set index arrays:
*---------------------------------------------------------------------*
      IF (ISYMIAJK .NE. MULD2H(ISYMT1,ISYOVOV))  THEN
         CALL QUIT('SYMMETRY MISMATCH IN CCG_TRANS2.')
      END IF

      DO  ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO 
        NTAIKJ(ISYM) = ICOUNT
      END DO

*---------------------------------------------------------------------*
* calculate (ia|jk) integrals and store as I^j(ai,k)
*---------------------------------------------------------------------*
      DO ISYMBJ = 1, NSYM
         ISYMAI  = MULD2H(ISYOVOV,ISYMBJ)
            
         DO ISYMB  = 1, NSYM
            ISYMK   = MULD2H(ISYMB,ISYMT1)
            ISYMAIK = MULD2H(ISYMAI,ISYMK)
            ISYMJ   = MULD2H(ISYMBJ,ISYMB)
 
            DO J = 1, NRHF(ISYMJ)

              NBJ    = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + 1
              KOFF1  = IT1AM(ISYMB,ISYMK) + 1
              KOFF2  = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1) + 1
              KOFF3  = ITAIKJ(ISYMAIK,ISYMJ) + NT2BCD(ISYMAIK)*(J-1) 
     &                                       + IT2BCD(ISYMAI,ISYMK) + 1

              NTOTB  = MAX(NVIR(ISYMB),1)
              NTOTAI = MAX(NT1AM(ISYMAI),1)

              CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),NVIR(ISYMB),
     &                   ONE,XOVOV(KOFF2),NTOTAI,T1AMP(KOFF1),NTOTB,
     &                   ZERO,XIAJK(KOFF3),NTOTAI)

            END DO

         END DO
      END DO

      CALL QEXIT('CCG_TRANS4')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_SORT1 */
*=====================================================================*
      SUBROUTINE CCG_SORT1(X,Y,ISYMX,IOPT)
*---------------------------------------------------------------------*
*
*     Purpose: resort integrals stored as X(ai,k,j) to
*
*              IOPT = 0    --    Y(ai,j,k) 
*                     1    --    Y(ai,jk) 
*                     2    --    Y(ai,kj) 
*                     4    --    Y(a,jki) 
*                     5    --    Y(jki,a)
*
*     Christof Haettig, August 1998
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
 
      INTEGER ISYMX, IOPT

#if defined (SYS_CRAY)
      REAL X(*), Y(*)
#else
      DOUBLE PRECISION X(*), Y(*)
#endif
      INTEGER ITAIKJ(8,8), NTAIKJ(8)
      INTEGER ISYM,ICOUNT,ISYMJ,ISYMAIK,ISYMAI,ISYMK,KOFF1,KOFF2,NJK
      INTEGER ISYMJK,ISYMAJ,ISYMI,ISYMA,ISYMAJK,ISYJIK,NAI,NJKI

      CALL QENTER('CCG_SORT1')

*---------------------------------------------------------------------*
* set index arrays for X integrals:
*---------------------------------------------------------------------*
      DO  ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO 
        NTAIKJ(ISYM) = ICOUNT
      END DO

      IF (IOPT.EQ.0) THEN
*---------------------------------------------------------------------*
* resort X(ai,k,j) to Y(aj,k,i)
*---------------------------------------------------------------------*
        DO ISYMJ = 1, NSYM
        DO ISYMK = 1, NSYM

           ISYMJK  = MULD2H(ISYMJ,ISYMK)
           ISYMAI  = MULD2H(ISYMJK,ISYMX)
           ISYMAIK = MULD2H(ISYMAI,ISYMK)
 
           DO J = 1, NRHF(ISYMJ)
           DO K = 1, NRHF(ISYMK)

              DO ISYMI = 1, NSYM

                 ISYMA   = MULD2H(ISYMI,ISYMAI)
                 ISYMAJ  = MULD2H(ISYMA,ISYMJ)
                 ISYMAJK = MULD2H(ISYMK,ISYMAJ)

                 KOFF1 = ITAIKJ(ISYMAIK,ISYMJ) + NT2BCD(ISYMAIK)*(J-1)
     &                 + IT2BCD(ISYMAI,ISYMK)  + NT1AM(ISYMAI)*(K-1) 
     &                 + IT1AM(ISYMA,ISYMI)    + 1

                 KOFF2 = ITAIKJ(ISYMAJK,ISYMI) 
     &                 + IT2BCD(ISYMAJ,ISYMK)  + NT1AM(ISYMAJ)*(K-1) 
     &                 + IT1AM(ISYMA,ISYMJ)    + NVIR(ISYMA)*(J-1) + 1

              DO I = 1, NRHF(ISYMI)

                 CALL DCOPY(NVIR(ISYMA),X(KOFF1),1,Y(KOFF2),1)

                 KOFF2 = KOFF2 + NT2BCD(ISYMAJK)
                 KOFF1 = KOFF1 + NVIR(ISYMA)
 
              END DO
              END DO
           END DO
           END DO

        END DO
        END DO

      ELSE IF (IOPT.EQ.1 .OR. IOPT.EQ.2) THEN
*---------------------------------------------------------------------*
* resort X(ai,k,j) to Y(ai,jk) / Y(ai,kj)
*---------------------------------------------------------------------*
        DO ISYMJ = 1, NSYM
        DO ISYMK = 1, NSYM

           ISYMJK  = MULD2H(ISYMJ,ISYMK)
           ISYMAI  = MULD2H(ISYMJK,ISYMX)
           ISYMAIK = MULD2H(ISYMAI,ISYMK)
 
           DO J = 1, NRHF(ISYMJ)
           DO K = 1, NRHF(ISYMK)

              IF (IOPT.EQ.1) THEN
                 NJK   = IMATIJ(ISYMJ,ISYMK) + NRHF(ISYMJ)*(K-1) + J
              ELSE
                 NJK   = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(J-1) + K
              END IF

              KOFF1 = ITAIKJ(ISYMAIK,ISYMJ) + NT2BCD(ISYMAIK)*(J-1)
     &              + IT2BCD(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1) + 1

              KOFF2 = ISAIKL(ISYMAI,ISYMJK) + NT1AM(ISYMAI)*(NJK-1) + 1

              CALL DCOPY(NT1AM(ISYMAI),X(KOFF1),1,Y(KOFF2),1)

           END DO
           END DO

        END DO
        END DO

      ELSE IF (IOPT.EQ.4) THEN
*---------------------------------------------------------------------*
* resort X(ai,k,j) to Y(a,jki)
*---------------------------------------------------------------------*
        DO ISYMJ = 1, NSYM
        DO ISYMK = 1, NSYM
        DO ISYMI = 1, NSYM

           ISYMJK  = MULD2H(ISYMJ,ISYMK)
           ISYJIK  = MULD2H(ISYMJK,ISYMI)
           ISYMAI  = MULD2H(ISYMJK,ISYMX)
           ISYMAIK = MULD2H(ISYMAI,ISYMK)
           ISYMA   = MULD2H(ISYMI,ISYMAI)
 
           DO J = 1, NRHF(ISYMJ)
           DO K = 1, NRHF(ISYMK)
           DO I = 1, NRHF(ISYMI)

             NJK   = IMATIJ(ISYMJ,ISYMK)  + NRHF(ISYMJ)*(K-1) + J
             NJKI  = IMAIJK(ISYMJK,ISYMI) + NMATIJ(ISYMJK)*(I-1) + NJK
             NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I-1) + 1

             KOFF1 = ITAIKJ(ISYMAIK,ISYMJ) + NT2BCD(ISYMAIK)*(J-1)
     &             + IT2BCD(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1) + NAI

             KOFF2 = ISJIKA(ISYJIK,ISYMA) + NVIR(ISYMA)*(NJKI-1) + 1

             CALL DCOPY(NVIR(ISYMA),X(KOFF1),1,Y(KOFF2),1)

           END DO
           END DO
           END DO

        END DO
        END DO
        END DO

      ELSE IF (IOPT.EQ.5) THEN
*---------------------------------------------------------------------*
* resort X(ai,k,j) to Y(jki,a)
*---------------------------------------------------------------------*
        DO ISYMJ = 1, NSYM
        DO ISYMK = 1, NSYM
        DO ISYMI = 1, NSYM

           ISYMJK  = MULD2H(ISYMJ,ISYMK)
           ISYJIK  = MULD2H(ISYMJK,ISYMI)
           ISYMAI  = MULD2H(ISYMJK,ISYMX)
           ISYMAIK = MULD2H(ISYMAI,ISYMK)
           ISYMA   = MULD2H(ISYMI,ISYMAI)
 
           DO J = 1, NRHF(ISYMJ)
           DO K = 1, NRHF(ISYMK)
           DO I = 1, NRHF(ISYMI)

             NJK   = IMATIJ(ISYMJ,ISYMK)  + NRHF(ISYMJ)*(K-1) + J
             NJKI  = IMAIJK(ISYMJK,ISYMI) + NMATIJ(ISYMJK)*(I-1) + NJK
             NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I-1) + 1

             KOFF1 = ITAIKJ(ISYMAIK,ISYMJ) + NT2BCD(ISYMAIK)*(J-1)
     &             + IT2BCD(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1) + NAI

             KOFF2 = ISJIKA(ISYJIK,ISYMA) + NJKI

             CALL DCOPY(NVIR(ISYMA),X(KOFF1),1,Y(KOFF2),NMAIJK(ISYJIK))

           END DO
           END DO
           END DO

        END DO
        END DO
        END DO

*---------------------------------------------------------------------*
      ELSE
        CALL QUIT('Illegal option in CCG_SORT1.')
      END IF

      CALL QEXIT('CCG_SORT1')

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CCG_LXD */
*=====================================================================*
      SUBROUTINE CCG_LXD(FTWO,ISYFCK,DENS,ISYDNS,XLIAJB,ISYOVOV,IOPT)
*---------------------------------------------------------------------*
*
*     Purpose: calculate contraction of packed (squared)
*               L(ia,jb) integrals with a density matrix D(jb)
*
*     FTWO(ai)  = sum(bj) L(ia,jb) * DENS(jb)
*
*     IOPT = 0 : packed integrals
*     IOPT = 1 : squared integrals
*
*     needed, e.g., for the two electron contribution to 
*     the first order responses of the Fock matrix
*
*     Christof Haettig, August 1996
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"     
#include "ccsdsym.h"

      INTEGER ISYFCK, ISYDNS, ISYOVOV, IOPT, KOFF1
      
#if defined (SYS_CRAY)
      REAL FTWO(*), DENS(*), XLIAJB(*)
      REAL ZERO, ONE
#else
      DOUBLE PRECISION FTWO(*), DENS(*), XLIAJB(*)
      DOUBLE PRECISION ZERO, ONE
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      INTEGER ISYMAI, ISYMBJ, NAI, NBJ, NAIBJ, NTOTAI, NTOTBJ

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCG_LXD')

* check symmetries:
      IF (ISYFCK .NE. MULD2H(ISYDNS,ISYOVOV)) THEN
        CALL QUIT('Symmetry mismatch in CCG_LXD.')
      END IF

* initialize result vector:
      CALL DZERO(FTWO, NT1AM(ISYFCK))

* calculate contraction FTWO(ia) = sum_bj L(ia,jb) * D(jb):
      ISYMAI = ISYFCK
      ISYMBJ = ISYDNS

      IF (IOPT.EQ.0 .AND. ISYMAI. EQ. ISYMBJ) THEN

        DO NAI = 1, NT1AM(ISYMAI)
          DO NBJ = 1, NT1AM(ISYMBJ)
 
            NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ) 
            FTWO(NAI) = FTWO(NAI) + XLIAJB(NAIBJ) * DENS(NBJ)
 
          END DO
        END DO

      ELSE IF (IOPT.EQ.0 .AND. ISYMAI .GT. ISYMBJ) THEN

        KOFF1  = IT2AM(ISYMBJ,ISYMAI) + 1
        NTOTBJ = MAX(1,NT1AM(ISYMBJ))

        CALL DGEMV('T',NT1AM(ISYMBJ),NT1AM(ISYMAI),ONE,XLIAJB(KOFF1),
     &             NTOTBJ,DENS,1,ZERO,FTWO,1)
       
      ELSE IF (IOPT.EQ.0 .AND. ISYMAI .LT. ISYMBJ) THEN

        KOFF1  = IT2AM(ISYMAI,ISYMBJ) + 1
        NTOTAI = MAX(1,NT1AM(ISYMAI))

        CALL DGEMV('N',NT1AM(ISYMAI),NT1AM(ISYMBJ),ONE,XLIAJB(KOFF1),
     &             NTOTAI,DENS,1,ZERO,FTWO,1)
       
      ELSE IF (IOPT.EQ.1) THEN

        KOFF1  = IT2SQ(ISYMAI,ISYMBJ) + 1
        NTOTAI = MAX(1,NT1AM(ISYMAI))

        CALL DGEMV('N',NT1AM(ISYMAI),NT1AM(ISYMBJ),ONE,XLIAJB(KOFF1),
     &             NTOTAI,DENS,1,ZERO,FTWO,1)

      END IF
 
      CALL QEXIT('CCG_LXD')

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_GMAT_NODDY(LISTL,IDLSTL,LISTB,IDLSTB,
     &                            LISTC,IDLSTC,
     &                            OMEGA1,OMEGA2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to G matrix transformation
*
*  (G T^B T^C)^eff_1,2 = (G T^B T^C)_1,2(CCSD) 
*                            + (G T^B T^C)_1,2(L_3)
*        
*     Written by Christof Haettig, April 2002 
*     based on CCSDT_FMAT_NODDY
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "dummy.h"
#include "priunit.h"
#include "ccsdinp.h"
#include "maxorb.h"
#include "ccsdsym.h"
#include "ccfield.h"
#include "ccorb.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER*3 LISTL, LISTB, LISTC
      INTEGER LWORK, IDLSTL, IDLSTB, IDLSTC

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL OMEGA1(NT1AMX), OMEGA2(NT2AMX)
      REAL DDOT, TWO, XNORMVAL, FREQL
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION OMEGA1(NT1AMX), OMEGA2(NT2AMX)
      DOUBLE PRECISION DDOT, TWO, XNORMVAL, FREQL
#endif
      PARAMETER (TWO = 2.0D0)

      CHARACTER*10 MODEL
      INTEGER KXINT, KXIAJB, KYIAJB, KT1AMP0, KLAMP0, KLAMH0, KEND1,
     &        KFOCK0, LWRK1, KLAMPB, KLAMHB, KLAMPC, KLAMHC, KT1AMPB,
     &        KT1AMPC, KINT1T0, KINT2T0, KINT1SBC, KINT2SBC,
     &        KBIOVVO, KBIOOVV, KBIOOOO, KBIVVVV,
     &        KCIOVVO, KCIOOVV, KCIOOOO, KCIVVVV,
     &        KBCIOVVO, KBCIOOVV, KBCIOOOO, KBCIVVVV,
     &        KOME1, KOME2, KL1AM, KL2AM, KL3AM, KT2AM, KFOCKD,
     &        KSCR1, KDUM, KFIELD, KFIELDAO
      INTEGER ISYMD, ILLL, IDEL, ISYDIS, IJ, NIJ, LUSIFC, INDEX,
     &        ISYMC, ISYMB, LUFOCK, KEND2, LWRK2, IOPT, ILSTSYM, ISYML

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      KDUM = 1
      CALL QENTER('CCSDT_GMAT_NODDY')

      IF (DIRECT) CALL QUIT('CCSDT_GMAT_NODDY: DIRECT NOT IMPLEMENTED')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KEND1   = KFOCKD + NORBT

      KFOCK0  = KEND1
      KEND1   = KFOCK0  + NORBT*NORBT

      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF

      KT1AMP0 = KEND1
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KT1AMPB = KEND1
      KLAMPB  = KT1AMPB + NT1AMX
      KLAMHB  = KLAMPB  + NLAMDT
      KEND1   = KLAMHB  + NLAMDT

      KT1AMPC = KEND1
      KLAMPC  = KT1AMPC + NT1AMX
      KLAMHC  = KLAMPC  + NLAMDT
      KEND1   = KLAMHC  + NLAMDT

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KBIOVVO = KEND1
      KBIOOVV = KBIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KBIOOOO = KBIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KBIVVVV = KBIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1   = KBIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KCIOVVO = KEND1
      KCIOOVV = KCIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KCIOOOO = KCIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KCIVVVV = KCIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1   = KCIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KBCIOVVO = KEND1
      KBCIOOVV = KBCIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KBCIOOOO = KBCIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KBCIVVVV = KBCIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1    = KBCIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KINT1SBC = KEND1
      KINT2SBC = KINT1SBC + NT1AMX*NVIRT*NVIRT
      KEND1    = KINT2SBC + NRHFT*NRHFT*NT1AMX

      KOME1   = KEND1
      KOME2   = KOME1  + NT1AMX
      KEND1   = KOME2  + NT1AMX*NT1AMX

      KL1AM   = KEND1
      KL2AM   = KL1AM  + NT1AMX
      KL3AM   = KL2AM  + NT1AMX*NT1AMX
      KEND1   = KL3AM  + NT1AMX*NT1AMX*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Get response  Lambda matrices:
*---------------------------------------------------------------------*
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 1
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KT1AMPB),WORK(KDUM))

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPB),
     &                 WORK(KLAMH0),WORK(KLAMHB),WORK(KT1AMPB),ISYMB)

      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT  = 1
      Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KT1AMPC),WORK(KDUM))

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPC),
     &                 WORK(KLAMH0),WORK(KLAMHC),WORK(KT1AMPC),ISYMC)

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)      

*---------------------------------------------------------------------*
*     Get the one electron integrals from the fields
*---------------------------------------------------------------------*
      IF (NONHF) THEN
         CALL DZERO(WORK(KFIELDAO),NORBT*NORBT)
         DO I = 1, NFIELD
            CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1),LWRK1,
     *                   EFIELD(I),1,LFIELD(I))
         ENDDO
         CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)

         ! calculate external field in zero-order lambda basis
         CALL CC_FCKMO(WORK(KFIELD),WORK(KLAMP0),WORK(KLAMH0),
     *                 WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)
      ENDIF

*---------------------------------------------------------------------*
*     Compute integrals needed for the following contributions:
*---------------------------------------------------------------------*

      CALL DZERO(WORK(KBIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KBIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KCIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KCIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KCIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KCIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KBCIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBCIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBCIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KBCIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KXIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)

      CALL DZERO(WORK(KINT1T0),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T0),NRHFT*NRHFT*NT1AMX)

      CALL DZERO(WORK(KINT1SBC),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2SBC),NRHFT*NRHFT*NT1AMX)

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  WORK(KDUM),DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1(WORK(KINT1T0),WORK(KINT2T0),
     &                       WORK(KLAMP0),WORK(KLAMH0),
     &                       WORK(KXINT),IDEL)

            CALL CC3_TRAN2(WORK(KXIAJB),WORK(KYIAJB),
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barB K-barC|B D)
            ! XINT2S = XINT2S + (C-barB K-barC|L J)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barB K|B-barC D)
            ! XINT2S = XINT2S + (C-barB K|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barC K-barB|B D)
            ! XINT2S = XINT2S + (C-barC K-barB|L J)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHB),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C K-barB|B-barC D)
            ! XINT2S = XINT2S + (C K-barB|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMHB),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barC K|B-barB D)
            ! XINT2S = XINT2S + (C-barC K|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C K-barC|B-barB D)
            ! XINT2S = XINT2S + (C K-barC|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBC),WORK(KINT2SBC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMHC),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C             (OV|VO)-B  (OO|VV)-B  (OO|OO)-B  (VV|VV)-B
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KBIOVVO),WORK(KBIOOVV),
     &                         WORK(KBIOOOO),WORK(KBIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KBIOVVO),WORK(KBIOOVV),
     &                         WORK(KBIOOOO),WORK(KBIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C             (OV|VO)-C  (OO|VV)-C  (OO|OO)-C  (VV|VV)-C
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KCIOVVO),WORK(KCIOOVV),
     &                         WORK(KCIOOOO),WORK(KCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KCIOVVO),WORK(KCIOOVV),
     &                         WORK(KCIOOOO),WORK(KCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C           (OV|VO)-BC  (OO|VV)-BC  (OO|OO)-BC  (VV|VV)-BC
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KBCIOVVO),WORK(KBCIOOVV),
     &                         WORK(KBCIOOOO),WORK(KBCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KBCIOVVO),WORK(KBCIOOVV),
     &                         WORK(KBCIOOOO),WORK(KBCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

      if (locdbg) then
        write(lupri,*) 'norm^2(OVVO-B):',
     &    ddot(NRHFT*NVIRT*NVIRT*NRHFT,WORK(KBIOVVO),1,WORK(KBIOVVO),1)
        write(lupri,*) 'norm^2(OOVV-B):',
     &    ddot(NRHFT*NVIRT*NVIRT*NRHFT,WORK(KBIOOVV),1,WORK(KBIOOVV),1)
        write(lupri,*) 'norm^2(OOOO-B):',
     &    ddot(NRHFT*NRHFT*NRHFT*NRHFT,WORK(KBIOOOO),1,WORK(KBIOOOO),1)
        write(lupri,*) 'norm^2(VVVV-B):',
     &    ddot(NVIRT*NVIRT*NVIRT*NVIRT,WORK(KBIVVVV),1,WORK(KBIVVVV),1)
        
        write(lupri,*) 'norm^2(OVVO-C):',
     &    ddot(NRHFT*NVIRT*NVIRT*NRHFT,WORK(KCIOVVO),1,WORK(KCIOVVO),1)
        write(lupri,*) 'norm^2(OOVV-C):',
     &    ddot(NRHFT*NVIRT*NVIRT*NRHFT,WORK(KCIOOVV),1,WORK(KCIOOVV),1)
        write(lupri,*) 'norm^2(OOOO-C):',
     &    ddot(NRHFT*NRHFT*NRHFT*NRHFT,WORK(KCIOOOO),1,WORK(KCIOOOO),1)
        write(lupri,*) 'norm^2(VVVV-C):',
     &    ddot(NVIRT*NVIRT*NVIRT*NVIRT,WORK(KCIVVVV),1,WORK(KCIVVVV),1)
        
        write(lupri,*) 'norm^2(OVVO-BC):',
     &   ddot(NRHFT*NVIRT*NVIRT*NRHFT,WORK(KBCIOVVO),1,WORK(KBCIOVVO),1)
        write(lupri,*) 'norm^2(OOVV-BC):',
     &   ddot(NRHFT*NVIRT*NVIRT*NRHFT,WORK(KBCIOOVV),1,WORK(KBCIOOVV),1)
        write(lupri,*) 'norm^2(OOOO-BC):',
     &   ddot(NRHFT*NRHFT*NRHFT*NRHFT,WORK(KBCIOOOO),1,WORK(KBCIOOOO),1)
        write(lupri,*) 'norm^2(VVVV-BC):',
     &   ddot(NVIRT*NVIRT*NVIRT*NVIRT,WORK(KBCIVVVV),1,WORK(KBCIVVVV),1)
        
        write(lupri,*) 'norm^2(XIAJB):',
     &      ddot(NT1AMX*NT1AMX,WORK(KXIAJB),1,WORK(KXIAJB),1)
        write(lupri,*) 'norm^2(YIAJB):',
     &      ddot(NT1AMX*NT1AMX,WORK(KYIAJB),1,WORK(KYIAJB),1)
        
        write(lupri,*) 'norm^2(INT1T0):',
     &      ddot(NT1AMX*NVIRT*NVIRT,WORK(KINT1T0),1,WORK(KINT1T0),1)
        write(lupri,*) 'norm^2(INT2T0):',
     &      ddot(NRHFT*NRHFT*NT1AMX,WORK(KINT2T0),1,WORK(KINT2T0),1)
        
        write(lupri,*) 'norm^2(INT1SBC):',
     &      ddot(NT1AMX*NVIRT*NVIRT,WORK(KINT1SBC),1,WORK(KINT1SBC),1)
        write(lupri,*) 'norm^2(INT2SBC):',
     &      ddot(NRHFT*NRHFT*NT1AMX,WORK(KINT2SBC),1,WORK(KINT2SBC),1)
      end if

*---------------------------------------------------------------------*
*     Compute L^0_3 multipliers:
*---------------------------------------------------------------------*
      ISYML = ILSTSYM(LISTL,IDLSTL)

      CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)

      IF (LISTL(1:3).EQ.'L0 ') THEN
        ! read left amplitudes from file and square up the doubles part
        IF (LWRK1 .LT. NT2AMX) 
     *     CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')
        IOPT = 3
        Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     *                WORK(KL1AM),WORK(KEND1))
        CALL CC_T2SQ(WORK(KEND1),WORK(KL2AM),ISYML)

        IF (NONHF .AND. LWRK1.LT.NT1AMX*NT1AMX*NT1AMX) 
     *     CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')

        CALL CCSDT_L03AM(WORK(KL3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                   WORK(KXIAJB),WORK(KFOCK0),WORK(KL1AM),
     *                   WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD),
     *                   WORK(KFIELD),WORK(KEND1))

      ELSE IF (LISTL(1:3).EQ.'L1 ' .OR. LISTL(1:3).EQ.'LE ' .OR.
     &         LISTL(1:3).EQ.'M1 ' .OR. LISTL(1:3).EQ.'N2 ' .OR.
     &         LISTL(1:3).EQ.'E0 '                              ) THEN

        CALL CCSDT_TBAR31_NODDY(WORK(KL3AM),FREQL,LISTL,IDLSTL,
     &                        WORK(KLAMP0),WORK(KLAMH0),
     &                        WORK(KFOCK0),WORK(KFOCKD),WORK(KSCR1),
     &                        WORK(KXIAJB),WORK(KINT1T0),WORK(KINT2T0),
     &                        WORK(KEND1),LWRK1)

      ELSE

        CALL QUIT('CCSDT_GMAT_NODDY> LISTL NOT AVAILABLE:'//LISTL)
      
      END IF

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)
  
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'NORM^2(L3AM):',
     *    DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KL3AM),1,WORK(KL3AM),1)
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[[H^BC,T^0_2],\tau_nu_1]|HF>
*                          and   <L_3|[[H^B, T^C_2],\tau_nu_1]|HF>
*                          and   <L_3|[[H^C, T^B_2],\tau_nu_1]|HF>
*---------------------------------------------------------------------*
      KT2AM  = KEND1
      KEND1  = KT2AM + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_GMAT_NODDY')
      ENDIF

      CALL DZERO(WORK(KOME1),NT1AMX)

      ! read T^0 doubles amplitudes from file and square up 
      IOPT   = 2
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYM0)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KBCIOOOO),WORK(KBCIOVVO),
     &                         WORK(KBCIOOVV),WORK(KBCIVVVV),
     &                         WORK(KT2AM))

      if (locdbg) then
        write(lupri,*) 'norm^2(t2am)-0:',
     &   ddot(nt1amx*nt1amx,work(kt2am),1,work(kt2am),1)
        write(lupri,*) 'norm^2(l3am):',
     *    ddot(nt1amx*nt1amx*nt1amx,work(kl3am),1,work(kl3am),1)
        write(lupri,*) 'norm^2(ome1) after <l3|[[Hbc,T0],tau]|0>:',
     &   ddot(nt1amx,work(kome1),1,work(kome1),1)
      end if


      ! read T^C doubles amplitudes from file and square up 
      IOPT   = 2
      Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMC)
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMC)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KBIOOOO),WORK(KBIOVVO),
     &                         WORK(KBIOOVV),WORK(KBIVVVV),
     &                         WORK(KT2AM))

      if (locdbg) then
        write(lupri,*) 'norm^2(t2am)-C:',
     &   ddot(nt1amx*nt1amx,work(kt2am),1,work(kt2am),1)
        write(lupri,*) 'norm^2(l3am):',
     *    ddot(nt1amx*nt1amx*nt1amx,work(kl3am),1,work(kl3am),1)
        write(lupri,*) 'norm^2(ome1) after <l3|[[Hb,T0],tau]|0>:',
     &   ddot(nt1amx,work(kome1),1,work(kome1),1)
      end if


      ! read T^B doubles amplitudes from file and square up 
      IOPT   = 2
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMB)
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMB)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KCIOOOO),WORK(KCIOVVO),
     &                         WORK(KCIOOVV),WORK(KCIVVVV),
     &                         WORK(KT2AM))

      if (locdbg) then
        write(lupri,*) 'norm^2(t2am)-B:',
     &   ddot(nt1amx*nt1amx,work(kt2am),1,work(kt2am),1)
        write(lupri,*) 'norm^2(l3am):',
     *    ddot(nt1amx*nt1amx*nt1amx,work(kl3am),1,work(kl3am),1)
        write(lupri,*) 'norm^2(ome1) after <l3|[[Hc,T0],tau]|0>:',
     &   ddot(nt1amx,work(kome1),1,work(kome1),1)
      end if

      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      END DO

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[H^BC,\tau_nu_2]|HF>
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)

      CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),WORK(KL3AM),
     *                         WORK(KINT1SBC),WORK(KINT2SBC))
      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         END DO
      END DO

      IF (LOCDBG) THEN
       WRITE(LUPRI,*)'G matrix transformation after triples contrib:'
        Call CC_PRP(OMEGA1,OMEGA2,1,1,1)
      END IF

      CALL QEXIT('CCSDT_GMAT_NODDY')

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_GMAT_NODDY                     *
*---------------------------------------------------------------------*
