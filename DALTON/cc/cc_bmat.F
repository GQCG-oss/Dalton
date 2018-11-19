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
c/* Deck CC_BMATRIX */
*=====================================================================*
      SUBROUTINE CC_BMATRIX( IBTRAN, NBTRAN, LISTA, LISTB, IOPTRES,
     &                       FILBMA, IBDOTS, BCONS, MXVEC, 
     &                       DO_O2,  WORK,   LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over B matrix transformations
*             (needed if the number of transformations exceeds the
*              limit MAXSIM defined on ccsdio.h )
*        
*     Written by Christof Haettig, March 1998.
*
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

      LOGICAL DO_O2
      CHARACTER*(*) LISTA, LISTB, FILBMA
      INTEGER IOPTRES
      INTEGER NBTRAN, MXVEC, LWORK
      INTEGER IBTRAN(3,NBTRAN)
      INTEGER IBDOTS(MXVEC,NBTRAN)
      
      REAL*8 WORK(LWORK) 
      REAL*8 BCONS(MXVEC,NBTRAN)

      INTEGER MAXBTRAN, NTRAN, ISTART, IBATCH, NBATCH

      CALL QENTER('CC_BMATRIX')
C
      MAXBTRAN = MAXSIM

      NBATCH = (NBTRAN+MAXBTRAN-1)/MAXBTRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Batching over B matrix transformations:'
        WRITE (LUPRI,*) 'nb. of batches needed:', NBATCH
      END IF
  
      DO IBATCH = 1, NBATCH
        ISTART = (IBATCH-1) * MAXBTRAN + 1
        NTRAN  = MIN(NBTRAN-(ISTART-1),MAXBTRAN)

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',ISTART
          WRITE (LUPRI,*) '# transf.:',NTRAN
        END IF

        CALL CC_BMAT( IBTRAN(1,ISTART), NTRAN,
     &                LISTA, LISTB, IOPTRES, FILBMA, 
     &                IBDOTS(1,ISTART), BCONS(1,ISTART), 
     &                MXVEC, DO_O2, WORK, LWORK )

      END DO

      CALL QEXIT('CC_BMATRIX')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_BMATRIX                           *
*---------------------------------------------------------------------*

*---------------------------------------------------------------------*
c/* Deck CC_BMAT */
*=====================================================================*
      SUBROUTINE CC_BMAT( IBTRAN, NBTRAN, LISTA, LISTB, IOPTRES,
     &                    FILBMA, IBDOTS, BCONS, MXVEC, DO_O2,
     &                    WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: AO-direct calculation of a linear transformation of two
*             CC amplitude vectors, T^A and T^B, with the CC B matrix
*             (derivatives of the CC lagrangian with respect to t)
*          
*             The linear transformations are calculated for a list
*             of T^A vectors and a list of T^B vectors: 
*
*                LISTA       -- type of T^A vectors
*                LISTB       -- type of T^B vectors
*                IBTRAN(1,*) -- indeces of T^A vectors
*                IBTRAN(2,*) -- indeces of T^B vectors
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
*                          file by a call to CC_WARSP, FILBMA is used
*                          as list type and IBTRAN(3,*) as list index
*                          NOTE that IBTRAN(3,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by FILBMA and the indeces from IBDOTS
*                          the result of the dot products is returned
*                          in the BCONS array
*
*     Written by Christof Haettig, Januar/Februar 1997.
*     BF terms rewritten in October 1998, Christof Haettig
*     CC3 noddy version, April 2002, Christof Haettig
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "mxcent.h"
#include "ccsdio.h"
#include "ccorb.h"
#include "cciccset.h"
#include "cbieri.h"
#include "distcl.h"
#include "iratdef.h"
#include "eritap.h"
#include "ccisao.h"
#include "ccfield.h"
#include "aovec.h"
#include "blocks.h"
#include "second.h"
#include "ccnoddy.h"
#include "ccr1rsp.h"
#include "r12int.h"
#include "ccsections.h"
#include "ccslvinf.h"
#include "qm3.h"
!#include "qmmm.h"

* local parameters:
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CC_BMAT> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL APPEND, NOAPPEND
      PARAMETER (APPEND = .TRUE., NOAPPEND = .FALSE.)

      INTEGER KDUM, IDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
      INTEGER ISYM0
      PARAMETER( ISYM0 = 1 ) ! symmetry of the reference state
      INTEGER ISYOVOV
      PARAMETER( ISYOVOV = 1 ) ! symmetry of (ia|jb) integrals 
      
      INTEGER LUBF, LUBFD, LUC, LUD, LUF, LUFK, LUR
      INTEGER LUAIBJ, LUCBAR, LUDBAR, LUBMAT
      CHARACTER*(8) BFFIL, CBAFIL, DBAFIL, CTFIL, DTFIL, RFIL
      CHARACTER*(8) FFIL, FKFIL, FNBFD, FNAIBJ
      PARAMETER (BFFIL ='CCCR_BFI', FNBFD ='CCBFDENS',
     &           CBAFIL='CCCR_CBA', DBAFIL='CCCR_DBA',
     &           CTFIL ='CCCR_CIM', DTFIL ='CCCR_DIM', 
     &           FFIL  ='CCCR_FIM', FKFIL ='CCCR_FKI',
     &           FNAIBJ='CCB_AIBJ', RFIL  ='CCCR_RIM')


      CHARACTER*(1) RSPTYP
      CHARACTER*(*) LISTA, LISTB, FILBMA
      LOGICAL DO_O2
      INTEGER IOPTRES
      INTEGER NBTRAN, MXVEC, LWORK
      INTEGER IBTRAN(3,NBTRAN)
      INTEGER IBDOTS(MXVEC,NBTRAN)

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL ZERO, ONE, TWO, FREQ
      REAL DUM, XNORM, FF, DUMMY
      REAL BCONS(MXVEC,NBTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ZERO, ONE, TWO, FREQ
      DOUBLE PRECISION DUM, XNORM, FF, DUMMY
      DOUBLE PRECISION BCONS(MXVEC,NBTRAN)
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0)

      CHARACTER*(3) LIST, LIST2A, LIST2B
      CHARACTER*(10) MODEL, MODELW, CDUMMY
      INTEGER INDEXA(MXCORB_CC)
      INTEGER INTMED1(2,2*MAXSIM), NINT1
      INTEGER INTMEDA(2,MAXSIM), NINTA
      INTEGER INTMED2(4,MAXSIM), NINT2
      INTEGER I1HGH(0:MAXSIM), I2HGH(0:MAXSIM), NBATCH
      INTEGER IOFFCD(0:MAXSIM+1)
      INTEGER IADRD(MXCORB_CC,2*MAXSIM) ! big static array :-(
      INTEGER IT2F(MXCORB_CC,2*MAXSIM)  ! big static array :-(
      INTEGER KLAMP(2*MAXSIM), KLAMH(2*MAXSIM), KDENS(2*MAXSIM)
      INTEGER KFOCK(2*MAXSIM), KRHO2(2*MAXSIM), KRIM(2*MAXSIM)
      INTEGER KFOCKOO(2*MAXSIM), KFOCKOV(2*MAXSIM), KFOCKVV(2*MAXSIM)
      INTEGER KXBAR(2*MAXSIM), KYBAR(2*MAXSIM)
      INTEGER KOMEGA2(MAXSIM) 
      INTEGER KLAMPA(MAXSIM),KLAMHA(MAXSIM)
      INTEGER KLAMPB(MAXSIM),KLAMHB(MAXSIM)

      LOGICAL NEWFTERM
      PARAMETER (NEWFTERM = .TRUE.)

      LOGICAL LGAMMA, LO3BF, OSQSAV, OORSAV
      INTEGER ITRAN, ISYM, IDLST, IDLSTA, IDLSTB, IOPT, ICORE, ICON, IF
      INTEGER ISYMA, ISYMB, ISYMAB, ISYMD1, NTOSYM, IDEL, ISYDEL
      INTEGER IINT1, IINT2, ISYM1, ISYCDBAR, IDXA, ISYX4O, IOPTG
      INTEGER IINT1A, IINT1B, IINTA, ICDEL2, NTOT, ILLL, NUMDIS
      INTEGER IBATCH, IDEL2, IADRTH, IERR, IOFFCDB, IOPTB, IADRBFD
      INTEGER MT2BGD, MDISAO, MDSRHF, KINDXB, KCCFB1, NVEC2
      INTEGER MSCRATCH, MEMAVAIL, NNWORK, NWORK, NSECMAR
      INTEGER KFOCK0, KDENS0, KT1AMP0, KLAMP0, KLAMH0, KEND0, LWRK0
      INTEGER KEND, LWRK, KENDSV, LWRKSV, KFREE, LFREE, JEND1, KEND1
      INTEGER KEND2, LWRK2, JEND2, KEND3, LWRK3, KEND4, LWRK4, LWRK1
      INTEGER KODCL1, KODCL2, KODBC1, KODBC2, KRDBC1, KRDBC2
      INTEGER KODPP1, KODPP2, KRDPP1, KRDPP2, KRECNR, KWRKSV
      INTEGER KXINT, KDSRHF, KLIAJB, KFOCK0OO, KFOCK0OV, KFOCK0VV
      INTEGER LEN, LENR, LENBF, LENF, LENFK, LENALL, IADRF, IVEC
      INTEGER KXIAJB, KT2AMP0, KT2AMPA, KCDBAR, KTHETA0, KFCKC0
      INTEGER KTHETA1, KTHETA2, KT1AMPA, KT1AMPB, KXLAMPA, KXLAMHA
      INTEGER KFCKAOO, KFCKAVV, KFCKBOO, KFCKBVV, KDNSC0
      INTEGER KFCKABOO, KFCKABOV, KFCKABVV, KXAIBJ, KBDRHF, KDCRHF
      INTEGER KBF0, LUBF0, KCBAR0, KDBAR0, KX4O, KSCR, KSCR2, KSCR1
      INTEGER KLAMDPB, KLAMDHB, KLAMDPA, KLAMDHA, IOPTW, IDUMMY
      INTEGER NBSRHF(8), IBSRHF(8,8), ICOUNT, ISYMAK, ISYBET
      INTEGER IOPTTCME, IOPTWE, KTHETA1EFF, KTHETA2EFF, KATRAN2
      INTEGER IOPTWR12,LENMOD,KTHETAR12,KATRANR12,IAMP
      CHARACTER APROXR12*3

* external functions:
      INTEGER ICCSET1
      INTEGER ICCSET2
      INTEGER ILSTSYM

#if defined (SYS_CRAY) 
      REAL DTIME, CONVRT, TIMALL, TIMTRN, TIMIO, TIMPRE
      REAL TIMA, TIMBF, TIMF, TIME, TIMI, TIMC, TIMD, TIMIM0
      REAL TIMINT, TIMRDAO, TIMTRBT, TIMIMA, TIMIMAB, TIMFCK
      REAL DDOT, FREQLST
#else 
      DOUBLE PRECISION DDOT, FREQLST
      DOUBLE PRECISION DTIME, CONVRT, TIMALL, TIMTRN, TIMIO, TIMPRE
      DOUBLE PRECISION TIMA,TIMBF,TIMF,TIME,TIMI,TIMC,TIMD,TIMIM0
      DOUBLE PRECISION TIMINT, TIMRDAO, TIMTRBT, TIMIMA, TIMIMAB,TIMFCK
#endif 
  
      CALL QENTER('CC_BMAT')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_BMAT')
        IF (DIRECT) WRITE(LUPRI,'(/1X,A)') 'AO direct transformation'
        WRITE (LUPRI,*) 'LISTA : ',LISTA
        WRITE (LUPRI,*) 'LISTB : ',LISTB
        WRITE (LUPRI,*) 'FILBMA: ',FILBMA
        WRITE (LUPRI,*) 'NBTRAN: ',NBTRAN
        WRITE (LUPRI,*) 'IOPTRES:',IOPTRES
        CALL FLSHFO(LUPRI)
      END IF
      
      IF ( .not. (CCS .or. CC2 .or. CCSD .or. CC3) ) THEN
        WRITE(LUPRI,'(/1x,a)') 'CC_BMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_BMAT...'
        CALL QUIT('Unknown CC method in CC_BMAT.')
      END IF

      IF (LISTA(1:1).NE.'R' .OR. LISTB(1:1).NE.'R') THEN
        WRITE(LUPRI,*) 'LISTA and LISTB must refer to t-amplitude',
     &                    ' vectors in CC_BMAT.'
        CALL QUIT('Illegal LISTA or LISTB in CC_BMAT.')
      END IF

      IF (.NOT. DUMPCD) THEN
        WRITE(LUPRI,*) 'DUMPCD = ',DUMPCD
        WRITE(LUPRI,*) 'CC_BMAT requires DUMPCD=.TRUE.'
        CALL QUIT('DUMPCD=.FALSE. , CC_BMAT requires DUMPCD=.TRUE.')
      END IF

      IF (ISYMOP .NE. 1) THEN
        WRITE(LUPRI,*) 'ISYMOP = ',ISYMOP
        WRITE(LUPRI,*) 'CC_BMAT is not implemented for ISYMOP.NE.1'
        CALL QUIT('CC_BMAT is not implemented for ISYMOP.NE.1')
      END IF

      IF (NBTRAN .GT. MAXSIM) THEN
        WRITE(LUPRI,*) 'NBTRAN = ', NBTRAN
        WRITE(LUPRI,*) 'MAXSIM = ', MAXSIM
        WRITE(LUPRI,*) 'number of requested transformation is larger'
        WRITE(LUPRI,*) 'than the maximum number of allowed ',
     &                 'simultaneous transformation.'
        WRITE(LUPRI,*) 'Error in CC_BMAT: NBTRAN is larger than MAXSIM.'
        CALL QUIT('Error in CC_BMAT: NBTRAN is larger than MAXSIM.')
      END IF

      IF (IPRINT.GT.0) THEN
 
         WRITE (LUPRI,'(//1X,A1,50("="),A1)')'+','+'

         WRITE (LUPRI,'(1x,A52)')
     &         '|        B MATRIX TRANSFORMATION SECTION           |'

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
      TIMIO   = ZERO
      TIMPRE  = ZERO
      TIMFCK  = ZERO
      TIMF    = ZERO
      TIMA    = ZERO
      TIMBF   = ZERO
      TIME    = ZERO
      TIMI    = ZERO
      TIMC    = ZERO
      TIMD    = ZERO
      TIMINT  = ZERO
      TIMRDAO = ZERO
      TIMTRBT = ZERO
      TIMIM0  = ZERO
      TIMIMA  = ZERO
      TIMIMAB = ZERO

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
      ELSE IF (CC3) THEN
         MODELW = 'CC3       '
         IOPTW  = 3
         IOPTWE = 24
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_BMAT.')
      END IF
      IF (CCR12) THEN
        APROXR12 = '   '
        CALL CCSD_MODEL(MODELW,LENMOD,10,MODELW,10,APROXR12)
        IOPTWR12 = 32
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
         CALL QUIT('Illegal value of IOPTRES in CC_BMAT.')
      END IF
 
* precalculate symmetry array for BSRHF:
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO

*=====================================================================*
* build nonredundant arrays of response vectors and pairs of them
* for which intermediates have to be calculated
*=====================================================================*
      DTIME = SECOND()

* array for intermediates that depend on one response vector
      NINT1 = 0
      DO ITRAN = 1, NBTRAN
        I=ICCSET1(INTMED1,LISTA,IBTRAN(1,ITRAN),NINT1,2*MAXSIM,APPEND)
        I=ICCSET1(INTMED1,LISTB,IBTRAN(2,ITRAN),NINT1,2*MAXSIM,APPEND)
      END DO 

* array for intermediates that are only required for the A vectors:
      NINTA = 0
      DO ITRAN = 1, NBTRAN
        I=ICCSET1(INTMEDA,LISTA,IBTRAN(1,ITRAN),NINTA,MAXSIM,APPEND)
      END DO 

* array for intermediates that depend on two response vectors
      NINT2 = 0
      DO ITRAN = 1, NBTRAN
        I=ICCSET2(INTMED2,LISTA,IBTRAN(1,ITRAN), 
     &                    LISTB,IBTRAN(2,ITRAN),NINT2,MAXSIM,APPEND)
      END DO 


      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/A)')'List of response vector for '//
     &        'AO intermediates:'
        WRITE (LUPRI,'((/5X,2I5))') ((INTMED1(I,J),I=1,2),J=1,NINT1)
        WRITE (LUPRI,'(/A)') 'List of response vector for '//
     &       'MO intermediates:'
        WRITE (LUPRI,'((/5X,2I5))') ((INTMEDA(I,J),I=1,2),J=1,NINTA)
        WRITE (LUPRI,'(/A)') 'List of vector pairs for '//
     &       'AO F intermediates:'
        WRITE (LUPRI,'((/5X,4I5))') ((INTMED2(I,J),I=1,4),J=1,NINT2)
      END IF
 
      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* estimate scratch space requirements
*---------------------------------------------------------------------*
      DTIME = SECOND()

      MT2BGD = 0
      MDISAO = 0
      MDSRHF = 0
      DO ISYM = 1, NSYM
        MT2BGD = MAX(MT2BGD,NT2BGD(ISYM))
        MDISAO = MAX(MDISAO,NDISAO(ISYM))
        MDSRHF = MAX(MDSRHF,NDSRHF(ISYM))
      END DO

*     5 x a NT2BGD type intermediate  
*     + integral arrays + some reserve

      MSCRATCH = 5*MT2BGD + MDISAO + 10*N2BASX  
      IF (CCSD.OR.CCSDT) MSCRATCH = MAX(MSCRATCH,MDISAO+5*MDSRHF) 

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_BMAT> scratch space estimate MSCRATCH:',
     &                  MSCRATCH
        CALL FLSHFO(LUPRI)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* estimate memory for 'in core' version and batched versions:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      MEMAVAIL = LWORK - MSCRATCH

      NWORK  = 0
      NBATCH = 1
      IF (CCS) THEN
        NSECMAR = 10 * N2BASX
      ELSE IF (CC2 .OR. CCSD .OR. CCSDT) THEN
        NSECMAR = 10 * MT2BGD
      ELSE
          CALL QUIT('Unknown CC model in CC_BMAT.')
      END IF

      I1HGH(0) = 0
      I2HGH(0) = 0

* intermediates that dependent on one response vector:
* (see routine ccbpre1 for details)
      DO IINT1 = 1, NINT1
        LIST  = VTABLE(INTMED1(2,IINT1))
        IDLST = INTMED1(1,IINT1)
        ISYM  = ILSTSYM(LIST,IDLST)
 
        NNWORK = 2*NGLMDT(ISYM) + 2*N2BST(ISYM)
        IF (CCSD.OR.CCSDT) THEN
          NNWORK = 2*NGLMDT(ISYM)+NT2AOIJ(ISYM)+NEMAT1(ISYM)
        END IF

        IF( (NWORK+NNWORK+NSECMAR).GT.MEMAVAIL ) THEN
          I1HGH(NBATCH) = IINT1 - 1
          I2HGH(NBATCH) = 0
       
          NBATCH = NBATCH + 1
          NWORK  = 0
        END IF
        NWORK = NWORK + NNWORK
        IF (NWORK .GT. LWORK) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CC_BMAT. (01)'
          WRITE (LUPRI,*) 'Need at least:',NNWORK, ' words.'
          CALL FLSHFO(LUPRI)
          CALL QUIT('Insufficient work space in CC_BMAT. (01)')
        END IF
      END DO

* intermediates that dependent on two response vectors:
* (see routine ccbpre2 for details)
      DO IINT2 = 1, NINT2
        LIST2A = VTABLE(INTMED2(2,IINT2))
        LIST2B = VTABLE(INTMED2(4,IINT2))
        IDLSTA = INTMED2(1,IINT2)
        IDLSTB = INTMED2(3,IINT2)
        ISYMA  = ILSTSYM(LIST2A,IDLSTA)
        ISYMB  = ILSTSYM(LIST2B,IDLSTB)
        ISYMAB = MULD2H(ISYMA,ISYMB)
 
        IF (CCS) THEN
          NNWORK = 0
        ELSE IF (CC2) THEN
          NNWORK = NT2AM(ISYMAB) + 2*NGLMDT(ISYMA) + 2*NGLMDT(ISYMB)
        ELSE IF (CCSD.OR.CCSDT) THEN
          NNWORK = 2*NGLMDT(ISYMA) + 2*NGLMDT(ISYMB)
        ELSE
          CALL QUIT('Unknown CC model in CC_BMAT.')
        END IF

        IF( (NWORK+NNWORK+NSECMAR).GT.MEMAVAIL ) THEN
          I1HGH(NBATCH) = NINT1
          I2HGH(NBATCH) = IINT2 - 1
       
          NBATCH = NBATCH + 1
          NWORK  = 0
        END IF
        NWORK = NWORK + NNWORK
        IF (NWORK .GT. LWORK) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CC_BMAT. (02)'
          WRITE (LUPRI,*) 'Need at least:',NNWORK,' words.'
          CALL FLSHFO(LUPRI)
          CALL QUIT('Insufficient work space in CC_BMAT. (02)')
        END IF
      END DO
 
      I1HGH(NBATCH) = NINT1
      I2HGH(NBATCH) = NINT2

      IF   (LOCDBG .AND. (NBATCH.EQ.1)) THEN
        WRITE (LUPRI,*) 'CC_BMAT> one batch only... '//
     &                  'will be done in core.'
        WRITE (LUPRI,*) 'CC_BMAT> memory for intermediates: ', NWORK
        WRITE (LUPRI,*) 'CC_BMAT> remaining scratch space: ',LWORK-NWORK
        CALL FLSHFO(LUPRI)
      ELSE IF (LOCDBG .AND. (NBATCH.GT.1)) THEN
        WRITE (LUPRI,*) 'CC_BMAT> more than one batch... '//
     &                  'choose I/O algorithm.'
        WRITE (LUPRI,*) 'CC_BMAT> max. memory for intermediates: ',
     &                  MEMAVAIL
        WRITE (LUPRI,*) 'CC_BMAT> number of batches: ',NBATCH
        CALL FLSHFO(LUPRI)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* read zeroth-order singles amplitudes, allocate space for Fock matrix,
* and prepare zeroth-order lambda matrices and density: 
*---------------------------------------------------------------------*
      DTIME = SECOND()

      KFOCK0   = 1
      KFOCK0OO = KFOCK0   + N2BAST
      KFOCK0OV = KFOCK0OO + NMATIJ(ISYM0)
      KFOCK0VV = KFOCK0OV + NT1AMX
      KDENS0   = KFOCK0VV + NMATAB(ISYM0)
      KT1AMP0  = KDENS0   + N2BAST
      KLAMP0   = KT1AMP0  + NT1AMX
      KLAMH0   = KLAMP0   + NLAMDT
      KEND0    = KLAMH0   + NLAMDT

      IF (FROIMP.OR.FROEXP) THEN
        KFCKC0 = KEND0
        KDNSC0 = KFCKC0 + N2BAST
        KEND0  = KDNSC0 + N2BAST
      ELSE
        KFCKC0 = KDUM
        KDNSC0 = KDUM
      END IF

      LWRK0    = LWORK - KEND0
      IF (LWRK0 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_BMAT. (0)')
      END IF

* read zeroth order amplitudes:
      IOPT   = 1
      Call CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

* get unperturbed Lambda matrices:
      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND0),LWRK0)

* calculate the density matrix:
      ICORE = 1
      CALL CC_AODENS(WORK(KLAMP0),WORK(KLAMH0),WORK(KDENS0),
     &               ISYM0,ICORE, WORK(KEND0),LWRK0)

* calculate pure core contribution to the density matrix,
* and initialize core contribution to Fock matrix with zeros
      IF (FROIMP.OR.FROEXP) THEN
        ICORE = 0 ! exclude core contribution
        CALL CC_AODENS(WORK(KLAMP0),WORK(KLAMH0),WORK(KDNSC0),
     &                 ISYM0,ICORE, WORK(KEND0),LWRK0)
        CALL DSCAL(N2BAST,-ONE,WORK(KDNSC0),1)
        CALL DAXPY(N2BAST,+ONE,WORK(KDENS0),1,WORK(KDNSC0),1)
        CALL DZERO(WORK(KFCKC0),N2BAST)
      END IF 

* initialize Fock matrix with the one-electron integrals:
      CALL CCRHS_ONEAO(WORK(KFOCK0),WORK(KEND0),LWRK0)
      DO IF= 1, NFIELD
        FF = EFIELD(IF)
        CALL CC_ONEP(WORK(KFOCK0),WORK(KEND0),LWRK0,FF,1,LFIELD(IF) )
      END DO
C
C------------------------------------------------------------------------
C     CCMM, 03 JK+OC
C     Solvent/QMMM  contribution to one-electron integrals.
C     T^g contribution to transformation.
C------------------------------------------------------------------------
C
      IF (CCSLV) THEN
         IF (.NOT.CCMM) CALL CCSL_RHSTG(WORK(KFOCK0),WORK(KEND0),LWRK0)
         IF (CCMM) THEN
            IF (.NOT. NYQMMM) THEN
            CALL CCMM_RHSTG(WORK(KFOCK0),WORK(KEND0),LWRK0)
            ELSE IF (NYQMMM) THEN
              IF (HFFLD) THEN
                WRITE(LUPRI,*) 'Is it justified to do B transfor with
     &  with a HFFLD?'
                CALL QUIT('HFFLD not implemented for QR')
              ELSE
                CALL CCMM_ADDG(WORK(KFOCK0),WORK(KEND0),LWRK0)
              END IF
            END IF
         ENDIF
      ENDIF
C
C------------------------------------------------------------------------
C
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'norm of T1AMP0:',
     &        DDOT(NT1AMX,WORK(KT1AMP0),1,WORK(KT1AMP0),1)
        WRITE (LUPRI,*) 'norm of XLAMP0:',
     &        DDOT(NLAMDT,WORK(KLAMP0),1,WORK(KLAMP0),1)
        WRITE (LUPRI,*) 'norm of XLAMH0:',
     &        DDOT(NLAMDT,WORK(KLAMH0),1,WORK(KLAMH0),1)
        WRITE (LUPRI,*) 'norm of DENS0:',
     &        DDOT(N2BAST,WORK(KDENS0),1,WORK(KDENS0),1)
        WRITE (LUPRI,*) 'norm of FOCK0:',
     &        DDOT(N2BAST,WORK(KFOCK0),1,WORK(KFOCK0),1)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME
*---------------------------------------------------------------------*
* open files for BF, C, D, F and Fock matrix intermediates:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      CALL CCBOPEN(LUBF,LUCBAR,LUDBAR,LUC,LUD,LUF,LUFK,LUR,
     &             BFFIL,CBAFIL,DBAFIL,CTFIL,DTFIL,FFIL,FKFIL,RFIL,
     &             LENBF, LENF, LENFK, LENR,
     &             NINT1, NINT2, WORK(KEND0), LWRK0 )

* open file for effective densities in BF term:
      LUBFD  = -1
      LUAIBJ = -1
      IF (.NOT.(CCS.OR.CC2)) CALL WOPEN2(LUBFD,  FNBFD,  64, 0)
      IF (.NOT.CCS)          CALL WOPEN2(LUAIBJ, FNAIBJ, 64, 0)

* initialize offsets for C & D intermediates:
      ICDEL2 = 0

* initialize offset for F term integrals:
      IADRF = 1

      TIMPRE = TIMPRE + SECOND() - DTIME

*---------------------------------------------------------------------*
* precalculate effective densities for BF intermediates:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IADRBFD = 1

      IF (.NOT. (CCS .OR. CC2)) THEN
         DO IINT1 = 1, NINT1
            LIST   = VTABLE(INTMED1(2,IINT1))
            IDLST  = INTMED1(1,IINT1)
            ISYMA  = ILSTSYM(LIST,IDLST)

            KT1AMPA = KEND0
            KT2AMPA = KT1AMPA + NT1AM(ISYMA)
            KXLAMHA = KT2AMPA + NT2SQ(ISYMA)
            KXLAMPA = KXLAMHA + NGLMDT(ISYMA)
            KEND1   = KXLAMPA + NGLMDT(ISYMA)
            LWRK1   = LWORK   - KEND1

            IF (LWRK1 .LT. NT2AM(ISYMA)) THEN
               CALL QUIT('Insufficient work space in CC_BMAT.(CCBFDEN)')
            END IF

C           ------------------------------------------------------
C           read response amplitudes, scale and square T2 part
C           ------------------------------------------------------
            IOPT = 3
            CALL CC_RDRSP(LIST,IDLST,ISYMA,IOPT,MODEL,
     *                    WORK(KT1AMPA),WORK(KEND1))
            CALL CCLR_DIASCL(WORK(KEND1),TWO,ISYMA)
            CALL CC_T2SQ(WORK(KEND1),WORK(KT2AMPA),ISYMA)
      
C           ------------------------------------------------------
C           calculate response lambda matrices:
C           ------------------------------------------------------
            CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KXLAMPA),
     *                       WORK(KLAMH0),WORK(KXLAMHA),
     *                       WORK(KT1AMPA),ISYMA)

C           ---------------------------------------------------------
C           calculate effective density for BF term and store on disk
C           ---------------------------------------------------------
            IOPT = 3
            CALL CC_BFDEN(WORK(KT2AMPA),ISYMA, DUMMY, IDUMMY,
     *                    WORK(KLAMH0), ISYM0, WORK(KLAMH0),ISYM0,
     *                    WORK(KXLAMHA),ISYMA, DUMMY,  IDUMMY,
     *                    FNBFD,  LUBFD,IADRD, IADRBFD,
     *                    IINT1,  IOPT, .FALSE., WORK(KEND1),LWRK1)

         END DO

      END IF

      TIMBF  = TIMBF  + SECOND() - DTIME
      TIMIMA = TIMIMA + SECOND() - DTIME

*---------------------------------------------------------------------*
* if all vectors and intermediates fit into the memory, read all
* response vectors before the loop over AO integral shells:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (NBATCH .EQ. 1) THEN

            CALL CCBPRE1(INTMED1, 1, NINT1,
     &                   KRHO2, KLAMP, KLAMH, KDENS, KFOCK, KRIM,
     &                   LUBF,BFFIL,LENBF,LUFK,FKFIL,LENFK,
     &                   LUR,RFIL,LENR,
     &                   WORK(KLAMP0), WORK(KLAMH0),
     &                   WORK, LWORK, KEND0, JEND1 )
            KEND1 = JEND1

            CALL CCBPRE2(INTMED2,1,NINT2,LUF,FFIL,LENF,
     &                   KOMEGA2,KLAMPA,KLAMHA,KLAMPB,KLAMHB,
     &                   WORK(KLAMP0),WORK(KLAMH0),
     &                   WORK, LWORK, KEND1, JEND1           )
            KEND1 = JEND1

            IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'CC_BMAT> allocated work '//
     &              'space for intermediates:'
             WRITE (LUPRI,*) 'CC_BMAT> KRHO2 :',(KRHO2(I),I=1,NINT1)
             WRITE (LUPRI,*) 'CC_BMAT> KLAMP :',(KLAMP(I),I=1,NINT1)
             WRITE (LUPRI,*) 'CC_BMAT> KLAMH :',(KLAMH(I),I=1,NINT1)
             WRITE (LUPRI,*) 'CC_BMAT> KDENS :',(KDENS(I),I=1,NINT1)
             WRITE (LUPRI,*) 'CC_BMAT> KFOCK :',(KFOCK(I),I=1,NINT1)
             WRITE (LUPRI,*) 'CC_BMAT> KOMEGA2:',(KOMEGA2(I),I=1,NINT2)
             WRITE (LUPRI,*) 'CC_BMAT> KLAMPA:',(KLAMPA(I),I=1,NINT2)
             WRITE (LUPRI,*) 'CC_BMAT> KLAMHA:',(KLAMHA(I),I=1,NINT2)
             WRITE (LUPRI,*) 'CC_BMAT> KLAMPB:',(KLAMPB(I),I=1,NINT2)
             WRITE (LUPRI,*) 'CC_BMAT> KLAMHB:',(KLAMHB(I),I=1,NINT2)
             WRITE (LUPRI,*) 'CC_BMAT> KEND1:',KEND1
             CALL FLSHFO(LUPRI)
            END IF

      ELSE
        KEND1 = KEND0
      END IF

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_BMAT. (1)')
      END IF
 
      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* initialize integral calculation
*---------------------------------------------------------------------*
      DTIME = SECOND()

      KEND = KEND1
      LWRK = LWRK1

      IF (DIRECT) THEN
         NTOSYM = 1

         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND),LWRK,IPRERI)
         ELSE
           KCCFB1 = KEND
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND   = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK   = LWORK  - KEND
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND),LWRK,IPRERI)
   
           KEND = KFREE
           LWRK = LFREE
         END IF

         KENDSV = KEND
         LWRKSV = LWRK
      ELSE
         NTOSYM = NSYM
      END IF

      TIMINT = TIMINT + SECOND() - DTIME
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

          DTIME = SECOND()
 
          IF (DIRECT) THEN
            KEND = KENDSV
            LWRK = LWRKSV
C 
            IF (HERDIR) THEN
              CALL HERDI2(WORK(KEND),LWRK,INDEXA,ILLL,NUMDIS,
     &                    IPRINT)
            ELSE
              CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                    WORK(KODCL1),WORK(KODCL2),
     *                    WORK(KODBC1),WORK(KODBC2),
     *                    WORK(KRDBC1),WORK(KRDBC2),
     *                    WORK(KODPP1),WORK(KODPP2),
     *                    WORK(KRDPP1),WORK(KRDPP2),
     *                    WORK(KCCFB1),WORK(KINDXB),
     *                    WORK(KEND), LWRK,IPRERI)
            END IF
C
            KRECNR = KEND
            KEND   = KRECNR + (NBUFX(0) - 1)/IRAT + 1
            LWRK   = LWORK - KEND
  
            IF (LWRK .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_BMAT. (1a)')
            END IF
 
          ELSE
            NUMDIS = 1
          END IF
 
          TIMINT = TIMINT + SECOND() - DTIME
*---------------------------------------------------------------------*
*        if out of core: allocate memory and get response vectors:
*---------------------------------------------------------------------*
          DO IBATCH = 1, NBATCH
             KEND2 = KEND ! reset memory for each batch
     
             IF (LOCDBG) THEN
               WRITE (LUPRI,*) MSGDBG, 
     &               IBATCH,'-th. batch out of ',NBATCH
               WRITE (LUPRI,*) MSGDBG, 'I1:',I1HGH(IBATCH-1)+1,' -- ',
     &                                  I1HGH(IBATCH)
               WRITE (LUPRI,*) MSGDBG, 'I2:',I2HGH(IBATCH-1)+1,' -- ',
     &                                  I2HGH(IBATCH)
             END IF
 
             IF (NBATCH.GT.1) THEN
 
               DTIME = SECOND()

               CALL CCBPRE1(INTMED1,I1HGH(IBATCH-1)+1,I1HGH(IBATCH),
     &                      KRHO2, KLAMP, KLAMH, KDENS, KFOCK, KRIM,
     &                      LUBF,BFFIL,LENBF,LUFK,FKFIL,LENFK,
     &                      LUR,RFIL,LENR,
     &                      WORK(KLAMP0), WORK(KLAMH0),
     &                      WORK, LWORK, KEND2, JEND2 )
               KEND2 = JEND2
 
               CALL CCBPRE2(INTMED2,I2HGH(IBATCH-1)+1,I2HGH(IBATCH),
     &                      LUF,FFIL,LENF,
     &                      KOMEGA2,KLAMPA,KLAMHA,KLAMPB,KLAMHB,
     &                      WORK(KLAMP0),WORK(KLAMH0),
     &                      WORK, LWORK, KEND2, JEND2              )
               KEND2 = JEND2

               TIMPRE = TIMPRE + SECOND() - DTIME
 
             END IF

             LWRK2 = LWORK - KEND2
             IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_BMAT. (2)')
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
            
            KXINT  = KEND2
            KEND3  = KXINT  + NDISAO(ISYDEL)

            IF (CCSD.OR.CCSDT) THEN
               KBDRHF = KEND3
               KDCRHF = KBDRHF + NBSRHF(ISYDEL)
               KDSRHF = KDCRHF + NBSRHF(ISYDEL)
               KEND3  = KDSRHF + NDSRHF(ISYDEL)
            END IF

            LWRK3  = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_BMAT. (3)')
            END IF
 
            DTIME = SECOND()
            CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     &                  WORK(KRECNR),DIRECT)
            TIMRDAO = TIMRDAO + SECOND() - DTIME

            IF (CCSD.OR.CCSDT) THEN

               DTIME = SECOND()

               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMP0),ISYM0,
     &                     WORK(KEND3),LWRK3,ISYDEL)

               CALL CC_BFBSORT(WORK(KDSRHF),WORK(KBDRHF),ISYDEL)
               
               CALL CCB_CDSORT(WORK(KXINT),ISYDEL,WORK(KDCRHF),
     *                         WORK(KLAMP0),ISYM0,WORK(KEND3),LWRK3)

               TIMTRBT = TIMTRBT + SECOND() - DTIME
               TIMIM0  = TIMIM0  + SECOND() - DTIME

               KEND3 = KDSRHF
               LWRK3 = LWORK  - KEND3

            END IF

*           calculate zeroth-order Fock matrix (Fhat)
            IF (IBATCH .EQ. 1) THEN
              DTIME = SECOND()

              CALL CC_AOFOCK(WORK(KXINT), WORK(KDENS0),
     *                       WORK(KFOCK0),WORK(KEND3),
     *                       LWRK3,IDEL,ISYDEL,.FALSE.,DUMMY,ISYM0)

              IF (FROIMP.OR.FROEXP) THEN
                 CALL CC_AOFOCK(WORK(KXINT), WORK(KDNSC0),
     *                          WORK(KFCKC0),WORK(KEND3),
     *                          LWRK3,IDEL,ISYDEL,.FALSE.,DUMMY,ISYM0)
              END IF                           

              TIMFCK = TIMFCK + SECOND() - DTIME
              TIMIM0 = TIMIM0 + SECOND() - DTIME
            END IF
 
*           calculate intermediates that depend on one response vector:
 
            DO IINT1 = I1HGH(IBATCH-1)+1, I1HGH(IBATCH)
              LIST   = VTABLE(INTMED1(2,IINT1))
              IDLST  = INTMED1(1,IINT1)
              ISYM1  = ILSTSYM(LIST,IDLST)

*             calculate addresses for C & D intermediates:
              IT2DLR(IDEL,IINT1) = ICDEL2
              ICDEL2 = ICDEL2 + NT2BCD(MULD2H(ISYDEL,ISYM1))

              DTIME = SECOND()
              CALL CCBINT1(WORK(KXINT), WORK(KBDRHF), WORK(KDCRHF),
     &                     IDEL, ISYDEL,  WORK(KRHO2(IINT1)),  
     &                     WORK(KLAMP0),  WORK(KLAMH0), 
     &                     WORK(KLAMP(IINT1)),  WORK(KLAMH(IINT1)),
     &                     ISYM1, IINT1,
     &                     WORK(KDENS(IINT1)),  WORK(KFOCK(IINT1)),  
     &                     WORK(KRIM(IINT1)),
     &                     LUC, CTFIL, LUD, DTFIL, 
     &                     LUBFD, FNBFD, IADRD(1,IINT1), 
     &                     WORK(KEND3),   LWRK3,          
     &                     TIMFCK, TIMBF, TIMC, TIMD  )
              TIMIMA = TIMIMA + SECOND() - DTIME

              
               IF (LOCDBG) THEN
                 WRITE (LUPRI,'(2A,3i5)') 
     &                 ' CC_BMAT> returned form CCBINT1 for',
     &                 ' IDEL,ISYDEL,IINT1=',IDEL,ISYDEL,IINT1
                 IF (.NOT.(CCS.OR.CC2)) THEN
                   XNORM = DDOT(NT2AOIJ(ISYM1),WORK(KRHO2(IINT1)),1,
     &                               WORK(KRHO2(IINT1)),1)
                   WRITE (LUPRI,*) 'CC_BMAT> norm of BF int.:',XNORM
                 END IF
                 IF (.NOT.(CCSD.OR.CCSDT)) THEN
                   XNORM = DDOT(N2BST(ISYM1),WORK(KFOCK(IINT1)),1,
     &                                       WORK(KFOCK(IINT1)),1)
                   WRITE (LUPRI,*) 'CC_BMAT> norm of FOCK int.:',XNORM
                 END IF
                 CALL FLSHFO(LUPRI)
               END IF

            END DO
 
 
*           calculate A & B vector dependend intermediates:
 
            DO IINT2 = I2HGH(IBATCH-1)+1, I2HGH(IBATCH)
              LIST2A = VTABLE(INTMED2(2,IINT2))
              LIST2B = VTABLE(INTMED2(4,IINT2))
              IDLSTA = INTMED2(1,IINT2)
              IDLSTB = INTMED2(3,IINT2)
              ISYMA  = ILSTSYM(LIST2A,IDLSTA)
              ISYMB  = ILSTSYM(LIST2B,IDLSTB)
              ISYMAB = MULD2H(ISYMA,ISYMB)

              DTIME = SECOND()
              CALL CCBINT2(WORK(KXINT), IDEL, ISYDEL, 
     &                WORK(KOMEGA2(IINT2)), ISYMAB,
     &                LUAIBJ, FNAIBJ, IT2F(1,IINT2), IADRF, NEWFTERM,
     &                WORK(KLAMPA(IINT2)),WORK(KLAMHA(IINT2)),ISYMA,
     &                WORK(KLAMPB(IINT2)),WORK(KLAMHB(IINT2)),ISYMB,
     &                WORK(KLAMP0), WORK(KLAMH0), 
     &                WORK(KEND3),  LWRK3                           )
              TIMF    = TIMF    + SECOND() - DTIME
              TIMIMAB = TIMIMAB + SECOND() - DTIME

               IF (LOCDBG) THEN
                 WRITE (LUPRI,*) 
     &                 'CC_BMAT> returned form CCBINT2 for IDEL,',
     &                 'ISYDEL,IINT2=',IDEL,ISYDEL,IINT2
                 IF (CC2) THEN
                   XNORM = DDOT(NT2AM(ISYMAB),WORK(KOMEGA2(IINT2)),1,
     &                                        WORK(KOMEGA2(IINT2)),1)
                   WRITE (LUPRI,*) 'CC_BMAT> norm of F int.:',XNORM
                 END IF
                 CALL FLSHFO(LUPRI)
               END IF

            END DO
 
          END DO ! IDEL2
*---------------------------------------------------------------------*
*         end of the loop over integral distributions:
*         if batched I/O algorithm used, save result on disc:
*---------------------------------------------------------------------*
          IF (NBATCH.GT.1) THEN
            DTIME = SECOND()
            CALL CCBSAVE(IBATCH,  I1HGH, I2HGH, INTMED1, INTMED2,
     &                   KRHO2,   LUBF,  BFFIL, LENBF,
     &                   KOMEGA2, LUF,   FFIL,  LENF,
     &                   KFOCK,   LUFK,  FKFIL, LENFK,
     &                   KRIM,    LUR,   RFIL,  LENR,
     &                   NINT1,   NINT2, WORK,  LWORK )
            TIMIO = TIMIO + SECOND() - DTIME
          END IF


        END DO ! IBATCH
       END DO ! ILLL
      END DO ! ISYMD1
*=====================================================================*
* End of Loop over AO-integrals
*=====================================================================*

*---------------------------------------------------------------------*
* if in-core algorithm used, save results now on disc:
*---------------------------------------------------------------------*
      IF (NBATCH.EQ.1) THEN
        DTIME = SECOND()
        CALL CCBSAVE(1,       I1HGH, I2HGH, INTMED1, INTMED2,
     &               KRHO2,   LUBF,  BFFIL, LENBF,
     &               KOMEGA2, LUF,   FFIL,  LENF,
     &               KFOCK,   LUFK,  FKFIL, LENFK,
     &               KRIM,    LUR,   RFIL,  LENR,
     &               NINT1,   NINT2, WORK,  LWORK )
        TIMIO = TIMIO + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
* recover workspace:
*---------------------------------------------------------------------*
      KEND1 = KEND0
      LWRK1 = LWRK0


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Loop over AO-integrals completed ',
     &           ' & AO intermediates saved on file.'
        WRITE (LUPRI,*) 'recover work space: KEND1,LWRK1=',KEND1,LWRK1
        WRITE (LUPRI,*) 'norm of XLAMH0:',
     &        DDOT(NLAMDT,WORK(KLAMH0),1,WORK(KLAMH0),1)
        CALL FLSHFO(LUPRI)
      END IF

*=====================================================================*
* calculate CBAR and DBAR intermediates:
*=====================================================================*
      IF (.NOT. (CCS .OR. CC2)) THEN

* initialize offset:
      IOFFCD(0) = 0

* set KCDBAR to a dummy address:
      KCDBAR = KDUM

* read (ia|jb) and square them:
      KXIAJB = KEND1
      KEND2  = KXIAJB + NT2SQ(ISYOVOV)
      LWRK2  = LWORK - KEND2

      IF (LWRK2 .LT. NT2AM(ISYOVOV)) THEN
        CALL QUIT('Insufficient work space in CC_BMAT. (4)')
      END IF

      DTIME = SECOND()

      Call CCG_RDIAJB (WORK(KEND2),NT2AM(ISYOVOV))

      Call CC_T2SQ (WORK(KEND2), WORK(KXIAJB), ISYOVOV)

      TIMIO = TIMIO + SECOND() - DTIME

*-----------------------------
* zeroth order intermediates:
*-----------------------------
      ISYCDBAR = MULD2H(ISYM0,ISYOVOV)

* read vector:
      KT2AMP0 = KEND2
      KEND3   = KT2AMP0 + NT2AM(ISYM0)
      LWRK3   = LWORK - KEND3
   
      IF (LWRK3 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_BMAT. (5)')
      END IF

      IOPT = 2
      DTIME = SECOND()
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KT2AMP0))
      TIMIO = TIMIO + SECOND() - DTIME

* CBAR intermediate:
      IOPT = 2
      DTIME = SECOND()
      CALL CCB_CDBAR('C', WORK(KXIAJB),ISYOVOV, WORK(KT2AMP0),ISYM0,
     &                    WORK(KCDBAR),ISYCDBAR, WORK(KEND3),LWRK3,
     &                    CBAFIL, LUCBAR, IOFFCD(0), IOPT)
      TIMC   = TIMC   + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME

* DBAR intermediate:
      IOPT = 2
      DTIME = SECOND()
      CALL CCB_CDBAR('D', WORK(KXIAJB),ISYOVOV, WORK(KT2AMP0),ISYM0,
     &                    WORK(KCDBAR),ISYCDBAR, WORK(KEND3),LWRK3,
     &                    DBAFIL, LUDBAR, IOFFCD(0), IOPT)
      TIMD   = TIMD   + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME

* increment offset:
      IOFFCD(1) = IOFFCD(0) + NT2SQ(ISYCDBAR)


*---------------------------------------------
* calculate intermediates for all A responses:
*---------------------------------------------
      DO IDXA = 1, NINTA
        IDLSTA   = INTMEDA(1,IDXA)
        ISYMA    = ILSTSYM(LISTA,IDLSTA)
        ISYCDBAR = MULD2H(ISYMA,ISYOVOV)

* read vector:
        KT2AMPA = KEND2
        KEND3   = KT2AMPA + NT2AM(ISYMA)
        LWRK3   = LWORK - KEND3
   
        IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (6)')
        END IF

        DTIME = SECOND()

        IOPT = 2
        CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &                WORK(KDUM),WORK(KT2AMPA)          )
        CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,ISYMA)

        TIMIO = TIMIO + SECOND() - DTIME

* CBAR intermediate:
        IF (LOCDBG) WRITE (LUPRI,*) 'CBAR', IDXA, ':'
        IOPT = 2
        DTIME = SECOND()
        CALL CCB_CDBAR('C', WORK(KXIAJB),ISYOVOV, WORK(KT2AMPA),ISYMA,
     &                      WORK(KCDBAR),ISYCDBAR, WORK(KEND3),LWRK3,
     &                      CBAFIL, LUCBAR, IOFFCD(IDXA), IOPT)
        TIMC   = TIMC + SECOND() - DTIME
        TIMIMA = TIMIMA + SECOND() - DTIME

* DBAR intermediate:
        IF (LOCDBG) WRITE (LUPRI,*) 'DBAR', IDXA, ':'
        IOPT = 2
        DTIME = SECOND()
        CALL CCB_CDBAR('D', WORK(KXIAJB),ISYOVOV, WORK(KT2AMPA),ISYMA,
     &                      WORK(KCDBAR),ISYCDBAR, WORK(KEND3),LWRK3,
     &                      DBAFIL, LUDBAR, IOFFCD(IDXA), IOPT)
        TIMD   = TIMD + SECOND() - DTIME
        TIMIMA = TIMIMA + SECOND() - DTIME

* increment offset:
        IOFFCD(IDXA+1) = IOFFCD(IDXA) + NT2SQ(ISYCDBAR)

      END DO

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'all CBAR & DBAR intermediates calculated...'
        CALL FLSHFO(LUPRI)
      END IF


* correct response BF^A/B intermediates for CCSD(R12) and higher:

      IF (CCR12) THEN
        DO IINT1 = 1, NINT1
          IDLST = INTMED1(1,IINT1)
          LIST  = VTABLE(INTMED1(2,IINT1))
          ISYM = ILSTSYM(LIST,IDLST) 

          !allocate scratch memory
          KSCR1 = KEND1
          KEND2 = KSCR1 + NT2AO(ISYM) 
          LWRK2 = LWORK - KEND2
          IF (LWRK2 .LE. 0) THEN
            CALL QUIT('Insufficient work space in CC_BMAT.')
          END IF
          
          CALL DZERO(WORK(KSCR1),NT2AO(ISYM))

          !calculate contribution:
          IOPT = 1
          IAMP = 2
          CALL CCRHS_BP(WORK(KSCR1),ISYM,IOPT,IAMP,DUMMY,IDUMMY,
     &                  IDUMMY,LIST,IDLST,KETSCL,WORK(KEND2),LWRK2)

          !read in response BF intermediate
          KSCR2 = KEND2 
          KEND2 = KSCR2 + NT2AOIJ(ISYM)
          LWRK2 = LWORK - KEND2
          IF (LWRK2 .LE. 0) THEN
            CALL QUIT('Insufficient work space in CC_BMAT.')
          END IF
          DTIME = SECOND()
          CALL CC_RVEC(LUBF,BFFIL,LENBF,NT2AOIJ(ISYM),IINT1,WORK(KSCR2))
          TIMIO = TIMIO + SECOND() - DTIME

          !transform beta index to occupied and add on BF intermediate:
          OSQSAV = OMEGSQ
          OORSAV = OMEGOR 
          OMEGSQ = .FALSE.
          OMEGOR = .FALSE.
          ICON = 4
          CALL CC_T2MO(DUMMY,DUMMY,1,WORK(KSCR1),DUMMY,WORK(KSCR2),
     &                 WORK(KLAMP0),WORK(KLAMP0),1,WORK(KEND2),LWRK2,
     &                 ISYM,ICON)
          OMEGSQ = OSQSAV
          OMEGOR = OORSAV

          !write out response BF intermediate
          DTIME = SECOND()
          CALL CC_WVEC(LUBF,BFFIL,LENBF,NT2AOIJ(ISYM),IINT1,
     &                 WORK(KSCR2)) 
          TIMIO = TIMIO + SECOND() - DTIME

        END DO
      END IF  ! (CCR12)
 
      END IF  ! (.NOT. (CCS .OR. CC2))

*=====================================================================*
* transform zeroth-order Fock matrix to the MO representation:
*=====================================================================*
      DTIME = SECOND()
      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)

      CALL CC_GATHEROO(WORK(KFOCK0),WORK(KFOCK0OO),ISYM0)
      CALL CC_GATHEROV(WORK(KFOCK0),WORK(KFOCK0OV),ISYM0)
      CALL CC_GATHERVV(WORK(KFOCK0),WORK(KFOCK0VV),ISYM0)

      TIMFCK = TIMFCK + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME
*=====================================================================*
* transform the response Fock^{*} matrices to MO representation
* and calculate the XBAR and YBAR intermediates:
*=====================================================================*
* read (ia|jb) integrals, calculate L(ia|jb) in place and square up:
* (stored and the upper end of the work space,
*  to keep the lower end free for intermediates)
      KLIAJB  = LWORK - NT2SQ(ISYOVOV)
      LFREE   = KLIAJB-1 

      LWRK1   = LWORK - KEND1
      IF ( LWRK1 .LT. (NT2AM(ISYOVOV)+NT2SQ(ISYOVOV)) ) THEN
        CALL QUIT('Insufficient work space in CC_BMAT. (6c)')
      END IF

      DTIME = SECOND()

      CALL CCG_RDIAJB(WORK(KEND1), NT2AM(ISYOVOV))

      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KEND1),ONE,ISYOVOV,IOPTTCME)

      CALL CC_T2SQ(WORK(KEND1), WORK(KLIAJB), ISYOVOV)

      TIMIO = TIMIO + SECOND() - DTIME


      DTIME = SECOND()
   
      DO IINT1 = 1, NINT1
        LIST  = VTABLE(INTMED1(2,IINT1))
        IDLST = INTMED1(1,IINT1)
        ISYM  = ILSTSYM(LIST,IDLST)
        CALL CCBINT3(LIST, IDLST,
     &               LUFK, FKFIL, LENFK, IINT1, KFOCK(IINT1),
     &               KFOCKOO(IINT1), KFOCKOV(IINT1), KFOCKVV(IINT1),
     &               KXBAR(IINT1), KYBAR(IINT1),
     &               WORK(KLIAJB), ISYOVOV,
     &               WORK(KLAMP0), WORK(KLAMH0),
     &               WORK, LFREE, KEND1, JEND1,
     &               TIMFCK,TIMIO,TIME)
        KEND1 = JEND1
      END DO

      TIMIMA = TIMIMA + SECOND() - DTIME

      
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LE. 0) THEN
        CALL QUIT('Insufficient work space in CC_BMAT. (7)')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'all FOCK intermediates transformed '//
     &                  'to MO basis...'
        CALL FLSHFO(LUPRI)
      END IF

      IF (IPRINT.GT.0) THEN

         WRITE (LUPRI,'(1X,A,F10.2,A)')
     &    '| time for zero order intermediat.:',TIMIM0 ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINT1,' sets of 1. ord. int.:',TIMIMA ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINT2,' sets of 2. ord. int.:',TIMIMAB,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A)')
     &    '| intermediates calculated in ',NBATCH,' batches          |'

         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'

         IF (IOPTRES.EQ.5) THEN
            WRITE (LUPRI,'(1X,A)')
     &         '| R vector | R vector |  # products  |             |'
            WRITE (LUPRI,'(1X,3(A,A3),A)') '|  ',LISTA(1:3), ' No. |  ',
     &       LISTB(1:3),' No. |   with ', FILBMA(1:3),
     &       '   |  time/secs  |'
         ELSE 
            WRITE (LUPRI,'(1X,A)')
     &         '| R vector | R vector |    result    |             |'
            WRITE (LUPRI,'(1X,A2,A3,2(A,A3),A)') '|  ',LISTA(1:3), 
     &        '  No. |  ', LISTB(1:3),' No. |   ', FILBMA(1:3),
     &        '  No.   |  time/secs  |'
         END IF

         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
      END IF

*=====================================================================*
* calculate B matrix transformations:
*=====================================================================*
      IADRTH = 1
      DO ITRAN = 1, NBTRAN

        IDLSTA = IBTRAN(1,ITRAN)
        IDLSTB = IBTRAN(2,ITRAN)

        ISYMA  = ILSTSYM(LISTA,IDLSTA)
        ISYMB  = ILSTSYM(LISTB,IDLSTB)
        ISYMAB = MULD2H(ISYMA,ISYMB)

        IINT1A = ICCSET1(INTMED1,LISTA,IDLSTA,NINT1,2*MAXSIM,NOAPPEND)
        IINT1B = ICCSET1(INTMED1,LISTB,IDLSTB,NINT1,2*MAXSIM,NOAPPEND)
        IINTA  = ICCSET1(INTMEDA,LISTA,IDLSTA,NINTA,MAXSIM,NOAPPEND)
        IINT2  = ICCSET2(INTMED2,LISTA,IDLSTA,
     &                           LISTB,IDLSTB,NINT2,MAXSIM,NOAPPEND)

        TIMTRN = SECOND()

*---------------------------------------------------------------------*
* allocate work space for the result vector:
*---------------------------------------------------------------------*
        IF (CCS) THEN
          KTHETA1 = KEND1
          KTHETA2 = KDUM
          KEND2   = KTHETA1 + NT1AM(ISYMAB)
        ELSE 
          KTHETA1 = KEND1
          KTHETA2 = KTHETA1 + NT1AM(ISYMAB)
          KEND2   = KTHETA2 + NT2AM(ISYMAB)
          IF (CCR12) THEN
            KTHETAR12 = KTHETA2 + NT2AM(ISYMAB)
            KEND2     = KTHETAR12 + NTR12AM(ISYMAB)
          END IF
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'B matrix transformation for ITRAN,',ITRAN
         WRITE (LUPRI,*) 'IADRTH:',IADRTH
         WRITE (LUPRI,*) 'LISTA,IDLSTA:',LISTA,IDLSTA
         WRITE (LUPRI,*) 'LISTB,IDLSTB:',LISTB,IDLSTB
         WRITE (LUPRI,*) 'ISYMA,ISYMB,ISYMAB:',ISYMA,ISYMB,ISYMAB
         WRITE (LUPRI,*) 'IINT1A,IINT1B,IINTA,IINT2:',IINT1A,IINT1B,
     &                    IINTA,IINT2
         CALL FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
* read the single excitation part of the response vectors and 
* calculate the response Lambda matrices:
*---------------------------------------------------------------------*
        DTIME = SECOND()

        KT1AMPA = KEND2
        KT1AMPB = KT1AMPA + NT1AM(ISYMA)
        KLAMDPB = KT1AMPB + NT1AM(ISYMB)
        KLAMDHB = KLAMDPB + NGLMDT(ISYMB)
        KLAMDPA = KLAMDHB + NGLMDT(ISYMB)
        KLAMDHA = KLAMDPA + NGLMDT(ISYMA)
        KEND2   = KLAMDHA + NGLMDT(ISYMA)
        LWRK2   = LWORK - KEND2

        IF (LWRK2 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (8)')
        END IF
      
        IOPT = 1

        CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KDUM))

        CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &                WORK(KT1AMPB),WORK(KDUM))

        CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMDPA), WORK(KLAMH0),
     &                   WORK(KLAMDHA),WORK(KT1AMPA),ISYMA)

        CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMDPB), WORK(KLAMH0),
     &                   WORK(KLAMDHB),WORK(KT1AMPB),ISYMB)

        TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* calculate Ftilde^{A} = [Fhat, T1^A] + Ftilde^{A,*} and Ftilde^{B}
*---------------------------------------------------------------------*
        DTIME = SECOND()

        KFCKAOO = KEND2
        KFCKAVV = KFCKAOO + NMATIJ(ISYMA)
        KFCKBOO = KFCKAVV + NMATAB(ISYMA)
        KFCKBVV = KFCKBOO + NMATIJ(ISYMB)
        KEND2   = KFCKBVV + NMATAB(ISYMB)
        LWRK2   = LWORK - KEND2

        IF (LWRK2 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (9)')
        END IF

        IF (.NOT.(CCSD.OR.CCSDT)) THEN

*         Ftilde^{A}, occupied/occupied blocks:
          CALL CCG_1ITROO(WORK(KFCKAOO),ISYMA, 
     &                    WORK(KFOCK0OV), ISYM0, WORK(KT1AMPA),ISYMA )

          CALL DAXPY(NMATIJ(ISYMA), ONE, 
     &               WORK(KFOCKOO(IINT1A)),1, WORK(KFCKAOO), 1)

*         Ftilde^{B}, occupied/occupied blocks:
          CALL CCG_1ITROO(WORK(KFCKBOO),ISYMB,
     &                    WORK(KFOCK0OV), ISYM0, WORK(KT1AMPB),ISYMB )

          CALL DAXPY(NMATIJ(ISYMB), ONE, 
     &               WORK(KFOCKOO(IINT1B)),1, WORK(KFCKBOO), 1)

*         Ftilde^{A}, virtual/virtual blocks:
          CALL CCG_1ITRVV(WORK(KFCKAVV),ISYMA, 
     &                    WORK(KFOCK0OV), ISYM0, WORK(KT1AMPA),ISYMA  )

          CALL DAXPY(NMATAB(ISYMA), ONE, 
     &               WORK(KFOCKVV(IINT1A)),1, WORK(KFCKAVV), 1)

*         Ftilde^{B}, virtual/virtual blocks:
          CALL CCG_1ITRVV(WORK(KFCKBVV),ISYMB, 
     &                    WORK(KFOCK0OV), ISYM0, WORK(KT1AMPB),ISYMB  )

          CALL DAXPY(NMATAB(ISYMB), ONE, 
     &               WORK(KFOCKVV(IINT1B)),1, WORK(KFCKBVV), 1)

          IF (LOCDBG) THEN
            XNORM=DDOT(NMATIJ(ISYMA),WORK(KFCKAOO),1,WORK(KFCKAOO),1)
            WRITE (LUPRI,*) 'Norm of FCKAOO:',XNORM
            XNORM=DDOT(NMATAB(ISYMA),WORK(KFCKAVV),1,WORK(KFCKAVV),1)
            WRITE (LUPRI,*) 'Norm of FCKAVV:',XNORM
            XNORM=DDOT(NMATIJ(ISYMB),WORK(KFCKBOO),1,WORK(KFCKBOO),1)
            WRITE (LUPRI,*) 'Norm of FCKBOO:',XNORM
            XNORM=DDOT(NMATAB(ISYMB),WORK(KFCKBVV),1,WORK(KFCKBVV),1)
            WRITE (LUPRI,*) 'Norm of FCKBVV:',XNORM
            CALL FLSHFO(LUPRI)
          END IF
        
        END IF

        TIMFCK = TIMFCK + SECOND() - DTIME
*---------------------------------------------------------------------*
* initialize the singles part of the result vector THETA:
*---------------------------------------------------------------------*
        CALL DZERO(WORK(KTHETA1),NT1AM(ISYMAB))

*---------------------------------------------------------------------*
* contributions which do not require any double amplitudes:
*---------------------------------------------------------------------*

*------------------------------------------------------------
* F contribution: transform (a i~|delta j~) integrals to MO
* for CC2 add here also remaining part of the F contribution:
*------------------------------------------------------------
      DTIME = SECOND()

      IF (CCSTST.AND.(.NOT.CCS)) CALL DZERO(WORK(KTHETA2),NT2AM(ISYMAB))

      IF ( CCSD .OR. CCSDT) THEN

         KXAIBJ = KEND2
         KEND3  = KXAIBJ + NT2SQ(ISYMAB)
         LWRK3  = LWORK  - KEND3
 
         IF (LWRK3 .LE. 0) THEN
            CALL QUIT('Insufficient work space in CC_BMAT. (CC_IAJB2)')
         END IF

         CALL DZERO(WORK(KXAIBJ),NT2SQ(ISYMAB))

         IOPT = 1
         CALL CC_IAJB2(WORK(KXAIBJ),ISYMAB,IOPT,.FALSE.,.FALSE.,.FALSE.,
     &                 LUAIBJ,FNAIBJ,IT2F(1,IINT2),WORK(KLAMP0),ISYM0,
     &                 IDUMMY,CDUMMY,IDUMMY,DUMMY,IDUMMY,
     &                 WORK(KEND3),LWRK3)

         IOPT = 1
         CALL CC_T2PK(WORK(KTHETA2),WORK(KXAIBJ),ISYMAB,IOPT)

      ELSE IF ( .NOT. (CCS .OR. CCSTST) ) THEN
         
         LEN = NT2AM(ISYMAB)
         CALL CC_RVEC(LUF,FFIL,LENF,LEN,IINT2,WORK(KTHETA2))

      END IF

      IF ( LOCDBG .AND. .NOT.(CCS.OR.CCSTST) ) THEN
         WRITE (LUPRI,*) 'read F contribution from file:'
         XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
         WRITE (LUPRI,*) 'Norm of THETA2 after F contribution:',XNORM
         CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,1,1)
         CALL FLSHFO(LUPRI)
      END IF

      TIMF = TIMF + SECOND() - DTIME
*---------------------------------------
* BF contribution, LAM^A x LAM^B x BF^0:
*---------------------------------------
      IF (.NOT. (CCS .OR. CC2 .OR. CCSTST) ) THEN

        KBF0  = KEND2
        KEND3 = KBF0 + 2*NT2ORT(ISYM0)
        LWRK3 = LWORK - KEND3

        IF (LWRK3 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (12)')
        END IF

* read zeroth-order BF intermediate:
        LUBF0 = -1
        DTIME = SECOND()
        CALL GPOPEN(LUBF0,'CC_BFIM','OLD',' ','UNFORMATTED',KDUM,
     &              .FALSE.)
        READ(LUBF0) (WORK(KBF0-1+I),I=1,2*NT2ORT(1))
        CALL GPCLOSE(LUBF0,'KEEP')
        TIMIO = TIMIO + SECOND() - DTIME

* transform to MO representation using two response Lambda matrices:
* (skip the calculation of the Gamma intermediate)
        ICON   = 2
        IOPTG  = 0
        LGAMMA = .FALSE.
        LO3BF  = .FALSE.
        DTIME  = SECOND()
        CALL CC_T2MO3(DUM,DUM,1,WORK(KBF0),
     *                WORK(KTHETA2),DUM,DUM,DUM, 
     *                WORK(KLAMDPA),ISYMA,WORK(KLAMDPB),ISYMB,
     *                WORK(KEND3),LWRK3,ISYM0,ICON,
     *                LGAMMA,IOPTG,LO3BF,.FALSE.)
        TIMBF = TIMBF + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(2*NT2ORT(ISYM0),WORK(KBF0),1,WORK(KBF0),1)
          WRITE (LUPRI,*) 'read BF(0) intermediate from file, norm is:',
     &                    XNORM
          XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after BF contribution:',XNORM
          CALL FLSHFO(LUPRI)
        END IF

      END IF

*---------------------------------------
* C contribution, CBAR^0 x T1^A x T1^B:
*---------------------------------------
      IF (.NOT. (CCS .OR. CC2 .OR. CCSTST )) THEN
        KCBAR0 = KEND2
        KEND3  = KCBAR0 + NT2SQ(ISYM0)
        LWRK3 = LWORK - KEND3

        IF (LWRK3 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (13)')
        END IF

        DTIME = SECOND()

* read in CBAR^0 intermediate:
        CALL GETWA2(LUCBAR,CBAFIL,WORK(KCBAR0),IOFFCD(0)+1,NT2SQ(ISYM0))

        TIMIO = TIMIO + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(NT2SQ(ISYM0),WORK(KCBAR0),1,WORK(KCBAR0),1)
          WRITE (LUPRI,*) 'read CBAR0 intermediate from file, norm is:',
     &                     XNORM
          CALL FLSHFO(LUPRI)
        END IF


        DTIME = SECOND()
        CALL CCB_22CD(WORK(KTHETA2),ISYMAB,WORK(KCBAR0),ISYM0, 
     &                WORK(KT1AMPA),ISYMA, WORK(KT1AMPB),ISYMB,
     &                'C', WORK(KEND3),LWRK3)
        TIMC = TIMC + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after C contribution:',XNORM
        END IF

      END IF

*---------------------------------------
* D contribution, DBAR^0 x T1^A x T1^B:
*---------------------------------------
      IF (.NOT. (CCS .OR. CC2 .OR. CCSTST) ) THEN
        KDBAR0 = KEND2
        KEND3  = KDBAR0 + NT2SQ(ISYM0)
        LWRK3  = LWORK - KEND3

        IF (LWRK3 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (14)')
        END IF

        DTIME = SECOND()

* read in CBAR^0 intermediate:
        CALL GETWA2(LUDBAR,DBAFIL,WORK(KDBAR0),IOFFCD(0)+1,NT2SQ(ISYM0))

        TIMIO = TIMIO + SECOND() - DTIME

        DTIME = SECOND()
        CALL CCB_22CD(WORK(KTHETA2),ISYMAB,WORK(KDBAR0),ISYM0, 
     &                WORK(KT1AMPA),ISYMA, WORK(KT1AMPB),ISYMB,
     &                'D', WORK(KEND3),LWRK3)
        TIMD = TIMD + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(NT2SQ(ISYM0),WORK(KDBAR0),1,WORK(KDBAR0),1)
          WRITE (LUPRI,*) 'read DBAR0 intermediate from file, norm is:',
     &                    XNORM
          XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after D contribution:',XNORM
          CALL FLSHFO(LUPRI)
        END IF

      END IF


*---------------------------------------------------------------------*
* CCSD contributions that require the zeroth order double amplitudes:
*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2) ) THEN
        KT2AMP0 = KEND2
        KEND3   = KT2AMP0 + NT2SQ(ISYM0)
        LWRK3   = LWORK - KEND3 

        IF (LWRK3 .LT. NT2AM(ISYM0)) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (15)')
        END IF

        IOPT = 2
        DTIME = SECOND()
        CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KEND3))

        CALL CC_T2SQ(WORK(KEND3),WORK(KT2AMP0),ISYM0)
        TIMIO = TIMIO + SECOND() - DTIME

      END IF

*----------------------------------
* A contribution, T^0 x (ki|lj)^AB:
*----------------------------------
      IF (.NOT. (CCS .OR. CC2 .OR. CCSTST) ) THEN

        ISYX4O = MULD2H(ISYOVOV,ISYMAB)
        IF (MULD2H(ISYX4O,ISYM0).NE.ISYMAB) THEN
          CALL QUIT('Symmetry mismatch in CC_BMAT, A term '//
     &              'T^0 x (ki|lj)^AB.')
        END IF

        KXIAJB = KEND3
        KX4O   = KXIAJB + NT2AM(ISYOVOV)
        KEND4  = KX4O   + NGAMMA(ISYX4O)
        LWRK4  = LWORK - KEND4

        IF (LWRK4 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (16)')
        END IF

* read (ia|jb) integrals:
        DTIME = SECOND()
        Call CCG_RDIAJB (WORK(KXIAJB),NT2AM(ISYOVOV))
        TIMIO = TIMIO + SECOND() - DTIME

* calculate double one-index transformed (ik|jl) integrals:
        DTIME = SECOND()

        IOPT = 1
        CALL CCG_4O(WORK(KX4O),ISYX4O,WORK(KXIAJB),ISYOVOV,
     &              WORK(KT1AMPA),ISYMA,WORK(KT1AMPB),ISYMB,
     &              WORK(KEND4),LWRK4,IOPT)

* calculate the contribution to THETA2:
        IOPT = 1
        CALL CCRHS_A(WORK(KTHETA2),WORK(KT2AMP0),WORK(KX4O),
     &               WORK(KEND4),LWRK4,ISYX4O,ISYM0,IOPT)

        TIMA = TIMA + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(NGAMMA(ISYX4O),WORK(KX4O),1,WORK(KX4O),1)
          WRITE (LUPRI,*) 'Norm of X4O intermediate:',XNORM
          XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after A contribution:',XNORM
          CALL FLSHFO(LUPRI)
        END IF

      END IF ! (.NOT. (CCS .OR. CC2))

*------------------------------------
* E1 & E2 contributions, T^0 x F^AB:
*------------------------------------
      IF (.NOT. (CCS .OR. CC2 .OR. CCSTST) ) THEN
        KFCKABOO = KEND3
        KFCKABVV = KFCKABOO + NMATIJ(ISYMAB)
        KSCR     = KFCKABVV + NMATAB(ISYMAB)
        KEND4    = KSCR + MAX(NMATIJ(ISYMAB),NMATAB(ISYMAB))
        LWRK4    = LWORK - KEND4

        IF (LWRK4 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_BMAT. (17)')
        END IF

        DTIME = SECOND()

* calculate occ/occ block of double one-index transformed Fock matrix:
        Call CCG_1ITROO(WORK(KFCKABOO),       ISYMAB,
     &                  WORK(KFOCKOV(IINT1A)),ISYMA,
     &                  WORK(KT1AMPB),        ISYMB  )

        Call CCG_1ITROO(WORK(KSCR),           ISYMAB,
     &                  WORK(KFOCKOV(IINT1B)),ISYMB,
     &                  WORK(KT1AMPA),        ISYMA  )

        Call DAXPY(NMATIJ(ISYMAB),ONE, WORK(KSCR),1, WORK(KFCKABOO),1)

* calculate vir/vir block of double one-index transformed Fock matrix:
        Call CCG_1ITRVV(WORK(KFCKABVV),       ISYMAB,
     &                  WORK(KFOCKOV(IINT1A)),ISYMA,
     &                  WORK(KT1AMPB),        ISYMB  )

        Call CCG_1ITRVV(WORK(KSCR),           ISYMAB,
     &                  WORK(KFOCKOV(IINT1B)),ISYMB,
     &                  WORK(KT1AMPA),        ISYMA  )

        Call DAXPY(NMATAB(ISYMAB),ONE, WORK(KSCR),1, WORK(KFCKABVV),1)

* calculate the contribution to THETA2:
        CALL CCRHS_E(WORK(KTHETA2),WORK(KT2AMP0),WORK(KFCKABVV),
     &               WORK(KFCKABOO),WORK(KEND4),LWRK4,ISYM0,ISYMAB)

        TIME = TIME + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(NMATIJ(ISYMAB),WORK(KFCKABOO),1,WORK(KFCKABOO),1)
          WRITE (LUPRI,*) 'Norm of KFCKABOO intermediate:',XNORM
          XNORM=DDOT(NMATAB(ISYMAB),WORK(KFCKABVV),1,WORK(KFCKABVV),1)
          WRITE (LUPRI,*) 'Norm of KFCKABVV intermediate:',XNORM
          XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after E contribution:',XNORM
          WRITE (LUPRI,*) 'THETA2 after E contribution:'
          Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,1,1)
          CALL FLSHFO(LUPRI)
        END IF


      END IF ! (.NOT. (CCS .OR. CC2 .OR. CCSTST) ) 

*---------------------------------------------------------------------*
* CC2/CCSD contributions that require the response B double amplitudes,
* and/or the response A BF and GAMMA intermediates, add here
* contributions from the CBAR and DBAR intermediates
*---------------------------------------------------------------------*
      IOFFCDB = IOFFCD(IINTA)
      IOPTB   = 1

      CALL CCBMAT2(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,IOPTB,IOFFCDB,
     &             LISTB,  IDLSTB, IINT1A, ISYMA, ISYMB,
     &             WORK(KFOCKOV(IINT1A)), WORK(KFOCKVV(IINT1A)), 
     &             WORK(KFCKAVV), WORK(KFCKAOO), WORK(KFOCK0OV),
     &             WORK(KFCKC0),
     &             WORK(KYBAR(IINT1A)),  WORK(KXBAR(IINT1A)),
     &             BFFIL,  CTFIL,  DTFIL,   DBAFIL, CBAFIL, RFIL,
     &             LUBF,   LENBF,  LUC,LUD, LUCBAR, LUDBAR, LUR,LENR,
     &             WORK(KLAMDPB), WORK(KLAMDPA), WORK(KLAMDHA),
     &             WORK(KLAMP0),  WORK(KLAMH0),
     &             WORK(KEND2), LWRK2, LISTA, IDLSTA,
     &             TIMIO,TIMA,TIMBF,TIME,TIMC,TIMD,TIMI)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'returned from CCBMAT2 (first call).'
        IF (.NOT. (CCS .OR. CCSTST) ) THEN
          XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after these contributions:',
     &                    XNORM
        END IF
        CALL FLSHFO(LUPRI)
      END IF


*---------------------------------------------------------------------*
* CC2/CCSD contributions that require the response A double amplitudes,
* and/or the response B BF and GAMMA intermediates, skip here
* contributions from the CBAR and DBAR intermediates
*---------------------------------------------------------------------*
      IOFFCDB = -99 999 999 ! dummy address
      IOPTB   = 0

      CALL CCBMAT2(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,IOPTB,IOFFCDB,
     &             LISTA,  IDLSTA, IINT1B, ISYMB, ISYMA,
     &             WORK(KFOCKOV(IINT1B)), WORK(KFOCKVV(IINT1B)), 
     &             WORK(KFCKBVV), WORK(KFCKBOO), WORK(KFOCK0OV), 
     &             WORK(KFCKC0),
     &             WORK(KYBAR(IINT1B)),  WORK(KXBAR(IINT1B)),
     &             BFFIL,  CTFIL,  DTFIL,   DBAFIL, CBAFIL, RFIL,
     &             LUBF,   LENBF,  LUC,LUD, LUCBAR, LUDBAR, LUR, LENR,
     &             WORK(KLAMDPA), WORK(KLAMDPB), WORK(KLAMDHB),
     &             WORK(KLAMP0),  WORK(KLAMH0),
     &             WORK(KEND2), LWRK2, LISTB, IDLSTB, 
     &             TIMIO,TIMA,TIMBF,TIME,TIMC,TIMD,TIMI)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'returned from CCBMAT2 (second call).'
        IF (.NOT. (CCS .OR. CCSTST) ) THEN
          XNORM=DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
          WRITE (LUPRI,*) 'Norm of THETA2 after these contributions:',
     &                    XNORM
        END IF
        CALL FLSHFO(LUPRI)
      END IF

C     initialize R12 part of result:
      IF (CCR12) CALL DZERO(WORK(KTHETAR12),NTR12AM(ISYMAB))

      IF (CCSLV) THEN
        IF (.NOT. CCMM) CALL CCSL_BTR(WORK(KTHETA1),WORK(KTHETA2),
     *                                ISYMAB,LISTA,IDLSTA,ISYMA,
     *                                LISTB,IDLSTB,ISYMB,
     *                                MODEL,WORK(KEND2),LWRK2)
C
        IF (CCMM) THEN
          IF (.NOT. NYQMMM) THEN
            CALL CCMM_BTR(WORK(KTHETA1),WORK(KTHETA2),
     *                          ISYMAB,LISTA,IDLSTA,ISYMA,
     *                          LISTB,IDLSTB,ISYMB,
     *                          MODEL,WORK(KEND2),LWRK2)
          ELSE IF (NYQMMM) THEN
            RSPTYP='B'
            CALL CCMM_QRTRANSFORMER(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,
     *                          LISTA,IDLSTA,ISYMA,
     *                          LISTB,IDLSTB,ISYMB,
     *                          MODEL,RSPTYP,WORK(KEND2),LWRK2)
          END IF
        END IF
      END IF
*---------------------------------------------------------------------*
* if DO_O2 flag is set include A{x} t^y + A{y} t^x contribution
* to the second-order rhs vector (here we do the CCS/CC2/CCSD parts):
*---------------------------------------------------------------------*
      IF (DO_O2) THEN
         IF ((FILBMA(1:3).NE.'O2 ' .AND. FILBMA(1:3).NE.'EO1') .OR.
     &       IOPTRES.GE.5) THEN
           CALL QUIT('Illegal result type for DO_O2 flag in CC_BMAT')
         END IF

         KATRAN2 = KEND2 + NT1AM(ISYMAB)
         KATRANR12 = KATRAN2 + NT2AM(ISYMAB) 

         IF (LISTA.EQ.'R1 ') THEN
            CALL CCCR_AA(LRTLBL(IDLSTA),ISYMA,LISTB,IDLSTB,DUMMY,
     &                   WORK(KEND2),LWRK2)
            CALL DAXPY(NT1AM(ISYMAB),ONE,WORK(KEND2),1,WORK(KTHETA1),1)
            IF (.NOT.CCS) THEN
              CALL CCLR_DIASCL(WORK(KATRAN2),2.0D0,ISYMAB)
              CALL DAXPY(NT2AM(ISYMAB),ONE,WORK(KATRAN2),1,
     &                                     WORK(KTHETA2),1)
            END IF
            IF (CCR12) THEN
              CALL DAXPY(NTR12AM(ISYMAB),ONE,WORK(KATRANR12),1,
     &                                      WORK(KTHETAR12),1)
            END IF
         END IF

         IF (LISTB.EQ.'R1 ') THEN
            CALL CCCR_AA(LRTLBL(IDLSTB),ISYMB,LISTA,IDLSTA,DUMMY,
     &                   WORK(KEND2),LWRK2)
            CALL DAXPY(NT1AM(ISYMAB),ONE,WORK(KEND2),1,WORK(KTHETA1),1)
            IF (.NOT.CCS) THEN
              CALL CCLR_DIASCL(WORK(KATRAN2),2.0D0,ISYMAB)
              CALL DAXPY(NT2AM(ISYMAB),ONE,WORK(KATRAN2),1,
     &                                     WORK(KTHETA2),1)
            END IF
            IF (CCR12) THEN
              CALL DAXPY(NTR12AM(ISYMAB),ONE,WORK(KATRANR12),1,
     &                                      WORK(KTHETAR12),1)
            END IF
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate R12 contributions:
*
* C. Neiss,  June 2005
*---------------------------------------------------------------------* 
      IF (CCR12) THEN
        CALL CC_R12BMAT(WORK(KTHETA1),WORK(KTHETA2),
     &                  WORK(KTHETAR12),ISYMAB,
     &                  LISTA,IDLSTA,WORK(KT1AMPA),ISYMA,
     &                  LISTB,IDLSTB,WORK(KT1AMPB),ISYMB,
     &                  WORK(KLAMDPA),WORK(KLAMDHA),
     &                  WORK(KLAMDPB),WORK(KLAMDHB),
     &                  WORK(KLAMP0),WORK(KLAMH0),
     &                  WORK(KEND2),LWRK2)
      END IF

*---------------------------------------------------------------------*
* add CC3 contribution
*---------------------------------------------------------------------*
      IF (CCSDT) THEN
        IF (IOPTRES.LT.5) THEN
           KTHETA1EFF = KEND2
           KTHETA2EFF = KTHETA1EFF + NT1AM(ISYMAB)
           KEND2      = KTHETA2EFF + NT2AM(ISYMAB)
           LWRK2      = LWORK - KEND2
C
           CALL DZERO(WORK(KTHETA1EFF),NT1AM(ISYMAB))
           CALL DZERO(WORK(KTHETA2EFF),NT2AM(ISYMAB))
        END IF

        IF (NODDY_BMAT) THEN
         
           ! Old style noddy code:
           CALL CCSDT_BMAT_NODDY(LISTA,IDLSTA,LISTB,IDLSTB,IOPTRES,
     &                           WORK(KLAMP0),WORK(KLAMH0),
     &                           WORK(KTHETA1),WORK(KTHETA2),
     &                           WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     &                           IBDOTS,BCONS,FILBMA,ITRAN,
     &                           NBTRAN,MXVEC,WORK(KEND2),LWRK2)
    
           IF (DO_O2) THEN
             FREQ = FREQLST(FILBMA,IBTRAN(3,ITRAN))
             IF (LISTA.EQ.'R1 ') THEN
               CALL CCSDT_AAMAT_NODDY(IOPTRES,FREQ,LRTLBL(IDLSTA),ISYMA,
     &                                LISTB,IDLSTB,.FALSE.,
     &                                WORK(KTHETA1),WORK(KTHETA2),
     &                                WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     &                                IBDOTS,BCONS,FILBMA,ITRAN,
     &                                NBTRAN,MXVEC,WORK(KEND2),LWRK2)
             END IF
             IF (LISTB.EQ.'R1 ') THEN
               CALL CCSDT_AAMAT_NODDY(IOPTRES,FREQ,LRTLBL(IDLSTB),ISYMB,
     &                                LISTA,IDLSTA,.FALSE.,
     &                                WORK(KTHETA1),WORK(KTHETA2),
     &                                WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     &                                IBDOTS,BCONS,FILBMA,ITRAN,
     &                                NBTRAN,MXVEC,WORK(KEND2),LWRK2)
             END IF
C            ! noddy code based on similar intermediates as real code:
C            ! is here called as dummy routine (which means, that it
C            ! actually doesn't change the result vectors theta*)
C            write(lupri,*) 'call now ccsdt_o32_noddy...'
C            CALL CCSDT_O32_NODDY(LISTA,IDLSTA,LISTB,IDLSTB,
C    *                           FILBMA,IBTRAN(3,ITRAN),
C    *                           WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
C    *                           WORK(KEND2),LWRK2)
           END IF 

        ELSE

            IF (IOPTRES .LT. 5) THEN
               IF (DO_O2) THEN
                  ! OMEGAEFF = OMEGAEFF + contribution from aamat
                  ! We assume that the singles and doubles contirbutions
                  !(sitting in KTHETA1 and KTHETA2) are added in CC3_BMATSD

                  !Project triples part of AAMAT into singles and 
                  !doubles space
                  CALL CC3_AAMATSD(LISTA,IDLSTA,
     *                                LISTB,IDLSTB,
     &                                DUMMY,DUMMY, 
     &                                WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     &                                ISYMAB,WORK(KEND2),LWRK2)
C
                  !Calculate triples contributions that enter directly
                  !doubles space simultanously for AAMAT and BMAT
                  CALL CC3_AABMAT_DOUB(WORK(KTHETA2),
     *                               LISTA,IDLSTA,LISTB,IDLSTB,
     *                               WORK(KEND2),LWRK2)

c
c
                  ! OMEGAEFF = OMEGAEFF + contribution from bmatsd
                  ! At the end: OMEGAEFF = OMEGAEFF + OMEGA
c
                  CALL CC3_BMATSD(WORK(KTHETA1),WORK(KTHETA2),
     &                            WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     &                            ISYMAB,
     *                            LISTA,IDLSTA,LISTB,IDLSTB,
     *                            WORK(KEND2),LWRK2)
               ELSE
                 WRITE(LUPRI,*)'Second-order rhs equatioons only'
                 WRITE(LUPRI,*)'implemented for TXY and EfX.'
                 CALL QUIT('Wrong DO_O2 value in CC_BMAT ')
               END IF
            ELSE 
              WRITE(LUPRI,*)'IOPTRES = ',IOPTRES
              CALL QUIT('Illegal IOPTRES in CC_BMAT (real CC3 code)')
            END IF


        END IF

      END IF

*---------------------------------------------------------------------*
* write result vector to output:
*---------------------------------------------------------------------*
      DTIME = SECOND()

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
        IF (CCR12) THEN
          CALL PUTWA2(LUBMAT,FILBMA,WORK(KTHETAR12),IADRTH,
     &                NTR12AM(ISYMAB))
          IADRTH = IADRTH + NTR12AM(ISYMAB)
        END IF

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'B matrix transformation nb. ',ITRAN,
     &          ' saved on file.'
         WRITE (LUPRI,*) 'ADRESS, LENGTH:',
     &        IBTRAN(3,ITRAN),IADRTH-IBTRAN(3,ITRAN)
         XNORM = DDOT(NT1AM(ISYMAB),WORK(KTHETA1),1,WORK(KTHETA1),1)
         IF (.NOT.CCS) XNORM = XNORM +
     &        DDOT(NT2AM(ISYMAB),WORK(KTHETA2),1,WORK(KTHETA2),1)
         IF (CCR12) XNORM = XNORM +
     &        DDOT(NTR12AM(ISYMAB),WORK(KTHETAR12),1,WORK(KTHETAR12),1)
         WRITE (LUPRI,*) 'Norm:', XNORM

         Call AROUND('B matrix transformation written to file:')
         Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,1,1)
         IF (CCR12) CALL CC_PRPR12(WORK(KTHETAR12),ISYMAB,1,.TRUE.)
        END IF

      ELSE IF ( IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4 ) THEN

*        write to a sequential file by a call to CC_WRRSP/CC_WARSP,
*        use FILBMA as LIST type and IBTRAN(3,ITRAN) as index
         KTHETA0 = -999999
         IF (IOPTRES.EQ.3) THEN
           CALL CC_WRRSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LWRK2)
           IF (CCR12) THEN
             CALL CC_WRRSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTWR12,
     &                     MODELW,DUMMY,DUMMY,WORK(KTHETAR12),
     &                     WORK(KEND2),LWRK2)
           END IF
           IF (CCSDT) THEN
             CALL CC_WRRSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTWE,MODELW,
     &                     WORK(KTHETA0),WORK(KTHETA1EFF),
     &                     WORK(KTHETA2EFF),WORK(KEND2),LWRK2)
           END IF
         ELSE IF (IOPTRES.EQ.4) THEN
           IF (CCSDT) THEN
             CALL CC_WARSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTWE,MODELW,
     &                     WORK(KTHETA0),WORK(KTHETA1EFF),
     &                     WORK(KTHETA2EFF),WORK(KEND2),LWRK2)
           END IF
           CALL CC_WARSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LWRK2)
           IF (CCR12) THEN
             CALL CC_WARSP(FILBMA,IBTRAN(3,ITRAN),ISYMAB,IOPTWR12,
     &                     MODELW,DUMMY,DUMMY,WORK(KTHETAR12),
     &                     WORK(KEND2),LWRK2)
           END IF
         END IF

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Write B * ',LISTA,' * ',LISTB,
     &              ' transformation',
     &              ' as ',FILBMA,' type vector to file.'
           WRITE (LUPRI,*) 'index of inp. A vector:',IBTRAN(1,ITRAN)
           WRITE (LUPRI,*) 'index of inp. B vector:',IBTRAN(2,ITRAN)
           WRITE (LUPRI,*) 'index of result vector:',IBTRAN(3,ITRAN)
           NVEC2 = 1
           LEN   = NT1AM(ISYMAB) + NT2AM(ISYMAB)
           IF (CCS) THEN
             NVEC2 = 0
             LEN   = NT1AM(ISYMAB)
           END IF
           IF (CCR12) LEN = LEN + NTR12AM(ISYMAB)
           XNORM = DDOT(LEN,WORK(KTHETA1),1,WORK(KTHETA1),1)
           WRITE (LUPRI,*) 'norm^2 of result vector:',XNORM
           WRITE (LUPRI,*) 'Listing of the result vector:'
           CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,1,NVEC2)
           IF (CCR12) CALL CC_PRPR12(WORK(KTHETAR12),ISYMAB,1,.TRUE.)
           IF (CCSDT) THEN
             XNORM = DDOT(LEN,WORK(KTHETA1EFF),1,WORK(KTHETA1EFF),1)
             WRITE (LUPRI,*) 'norm^2 of eff. result vector:',XNORM
             WRITE (LUPRI,*) 'Listing of the eff. result vector:'
             CALL CC_PRP(WORK(KTHETA1EFF),WORK(KTHETA2EFF),ISYMAB,1,1)
           END IF
         END IF
      ELSE IF (IOPTRES.EQ.5) THEN
         CALL CCDOTRSP(IBDOTS,BCONS,IOPTW,FILBMA,ITRAN,NBTRAN,MXVEC,
     &                 WORK(KTHETA1),WORK(KTHETA2),ISYMAB,
     &                 WORK(KEND2),LWRK2)
         IF (CCR12) THEN
           CALL CCDOTRSP(IBDOTS,BCONS,IOPTWR12,FILBMA,ITRAN,NBTRAN,
     &                   MXVEC,DUMMY,WORK(KTHETAR12),ISYMAB,
     &                   WORK(KEND2),LWRK2)
         END IF
      ELSE
        CALL QUIT('Illegal value for IOPTRES in CC_BMAT.')
      END IF
      
      TIMIO = TIMIO + SECOND() - DTIME

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
* End of loop over B matrix transformations
*---------------------------------------------------------------------*
      END DO

*---------------------------------------------------------------------*
* close & remove scratch files:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (.NOT. (CCS.OR.CC2)) THEN
         CALL WCLOSE2(LUBF,   BFFIL,  'DELETE')
         CALL WCLOSE2(LUC,    CTFIL,  'DELETE')
         CALL WCLOSE2(LUD,    DTFIL,  'DELETE')
         CALL WCLOSE2(LUCBAR, CBAFIL, 'DELETE')
         CALL WCLOSE2(LUDBAR, DBAFIL, 'DELETE')
         CALL WCLOSE2(LUR,    RFIL,   'DELETE')
      END IF

      IF (.NOT. CCS) CALL WCLOSE2(LUAIBJ,FNAIBJ,'DELETE') 
      IF (.NOT. CCS) CALL WCLOSE2(LUF,   FFIL,  'DELETE')
      IF (CCSD.OR.CCSDT)      CALL WCLOSE2(LUBFD, FNBFD, 'DELETE')

      CALL WCLOSE2(LUFK, FKFIL, 'DELETE')

      TIMIO = TIMIO + SECOND() - DTIME

*---------------------------------------------------------------------*
* if IOPTRES=1 and enough work space available, read result
* vectors back into memory:
*---------------------------------------------------------------------*
      DTIME = SECOND()

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
            WRITE (LUPRI,*) 'Read B matrix transformation nb. ',NBTRAN
            WRITE (LUPRI,*) 'Adress, length, NORM:',IBTRAN(3,NBTRAN),
     &                      LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

      TIMIO = TIMIO + SECOND() - DTIME
*---------------------------------------------------------------------*
* close B matrix file, print timings & return
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (IOPTRES.EQ.0 ) THEN
        CALL WCLOSE2(LUBMAT, FILBMA, 'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
        CALL WCLOSE2(LUBMAT, FILBMA, 'DELETE')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_BMAT.')
      END IF

      TIMIO  = TIMIO + SECOND() - DTIME
      TIMALL = SECOND() - TIMALL

      IF (IPRINT.GE.0) THEN
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
         WRITE (LUPRI,'(1X,A,I4,A,F10.2,A)') 
     &     '| total time for',NBTRAN,' B transforms.:',TIMALL,' secs.|'
      IF (TIMALL .GE. 1.0D0) THEN
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
         CONVRT = 100.0D0/TIMALL
         WRITE (LUPRI,'(1X,"|  % of time used in ",A,":",
     &        F10.2,"      |")')
     &      'start up org.', TIMPRE*CONVRT,
     &      'Fock interm. ', TIMFCK*CONVRT,
     &      'ERI          ', TIMINT*CONVRT,
     &      'CCRDAO       ', TIMRDAO*CONVRT,
     &      '(**|k del)   ', TIMTRBT*CONVRT,
     &      'singles part ', TIMI*CONVRT,
     &      'A term       ', TIMA*CONVRT,
     &      'BF term      ', TIMBF*CONVRT,
     &      'F term       ', TIMF*CONVRT,
     &      'E term       ', TIME*CONVRT,
     &      'C term       ', TIMC*CONVRT,
     &      'D term       ', TIMD*CONVRT,
     &      'I/O          ', TIMIO*CONVRT
      END IF

         WRITE (LUPRI,'(1X,A1,50("="),A1,//)') '+','+'
      END IF
*=====================================================================*

      CALL QEXIT('CC_BMAT')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CC_BMAT
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBMAT2 */
*=====================================================================*
      SUBROUTINE CCBMAT2(THETA1, THETA2, ISYRES, IOPTB,  IOFFCDB,
     &                   LISTB,  IDLSTB, IINT1A,  ISYMA, ISYMB,
     &                   FOCKA,  FOCKVV, FCKAVV, FCKAOO, FCK0OV, FCKC0,
     &                   YBARA,  XBARA,  BFIL,   CTFIL,  DTFIL,   
     &                   DBAFIL, CBAFIL, RFIL,   LUBF,   LENBF, 
     &                   LUC,LUD, LUCBAR, LUDBAR,LUR,    LENR,
     &                   XLAMPB, XLAMPA, XLAMHA, XLAMP0, XLAMH0,  
     &                   WORK,   LWORK,  LISTA,  IDLSTA, TIMIO,   
     &                   TIMA,   TIMBF, TIME, TIMC, TIMD,  TIMI )
*---------------------------------------------------------------------*
*    Purpose: calculate CC2 & CCSD contributions to the B matrix
*             transformation that require the response B double
*             amplitudes and/or the response A BF and GAMMA 
*             intermediates
*
*    Written by Christof Haettig, Januar 1997.
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccorb.h"
#include "ccfield.h"
#include "second.h"
#include "ccsections.h"
#include "ccslvinf.h"
#include "qm3.h"
!#include "qmmm.h"


* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0, KDUM
      PARAMETER (ISYM0 = 1)
      PARAMETER (KDUM = +99 999 999)

#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO, HALF, FACB
#else
      DOUBLE PRECISION ZERO, ONE, TWO, HALF, FACB
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, HALF = 0.5d0)

      LOGICAL LGAMMA, LO3BF, LRCON, LGCON, FCKCON
      CHARACTER*(*) BFIL, CTFIL, DTFIL, DBAFIL, CBAFIL, RFIL
      CHARACTER*(*) LISTB, LISTA
      CHARACTER*(10) MODEL
      INTEGER LWORK, IOFFCDB
      INTEGER IOPTB, IOPTE
      INTEGER LUBF, LENBF, LUC, LUD, LUCBAR, LUDBAR, LUR, LENR
      INTEGER ISYMA, ISYMB, ISYRES, ISYMAB
      INTEGER IINT1A, IDLSTB, IDLSTA
      INTEGER KEND0, KGAMMA, KBFA, KEND1, LWRK1, LEN, ICON, IDUM
      INTEGER KT2AMPB, KGAMMAX, LWRK0, KEND2, LWRK2, KSCR, IOPT, IOPTG
      INTEGER KEMAT1A, KEMAT2A, NRHO, KRIM, KGIM, KCON, KFCKA
      INTEGER IF, KT1AMPA, KOV, KPERT, KT1AMPB, KONEHG, KONEHR

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL THETA1(*), THETA2(*), XLAMP0(*), XLAMH0(*)
      REAL XLAMPB(*), XLAMPA(*), XLAMHA(*)
      REAL FOCKA(*), FOCKVV(*), FCKAVV(*), FCKAOO(*) 
      REAL YBARA(*), XBARA(*), FCK0OV(*), FCKC0(*) 
      REAL DUM, XNORM, DTIME
      REAL TIMIO, TIMBF, TIMA, TIME, TIMI, TIMC, TIMD
      REAL DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION THETA1(*), THETA2(*), XLAMP0(*), XLAMH0(*)
      DOUBLE PRECISION XLAMPB(*), XLAMPA(*), XLAMHA(*)
      DOUBLE PRECISION FOCKA(*), FOCKVV(*), FCKAVV(*), FCKAOO(*)  
      DOUBLE PRECISION YBARA(*), XBARA(*), FCK0OV(*), FCKC0(*)
      DOUBLE PRECISION FF, DUM, XNORM, DTIME
      DOUBLE PRECISION TIMIO, TIMBF, TIMA, TIME, TIMI, TIMC, TIMD
      DOUBLE PRECISION DDOT
#endif

      CALL QENTER('CCBMAT2')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Entered CCBMAT2.'
        WRITE (LUPRI,*) 'norm of XLAMH0:',
     &       DDOT(NLAMDT,XLAMH0,1,XLAMH0,1)
        CALL FLSHFO(LUPRI)
      END IF

      ISYMAB = MULD2H(ISYMA,ISYMB)

      IF (ISYRES .NE. ISYMAB) THEN
         CALL QUIT('Symmetry mismatch in CCBMAT2.')
      END IF

*-----------------------------------------------
* allocate work space and read T1AMB amplitudes:
*-----------------------------------------------
      KT1AMPB = 1
      KCON    = KT1AMPB + NT1AM(ISYMB)
      KEND0   = KCON    + NT1AM(ISYMAB)
      LWRK0   = LWORK   - KEND0

      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CCBMAT2. (0)')
      END IF

      IOPT = 1
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KT1AMPB),WORK(KDUM))

*--------------------------------------------------
* for CCS calculate here J contribution and return:
*--------------------------------------------------
      IF (CCS .OR. CCSTST) THEN

         IOPT = 1
         CALL CCG_1ITRVO(WORK(KCON),ISYMAB,FCKAOO,FOCKVV,ISYMA,
     &                   WORK(KT1AMPB),ISYMB,IOPT)

         CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KCON),1,THETA1,1)

         CALL QEXIT('CCBMAT2')

         RETURN

      END IF
      

      KEMAT1A = KEND0
      KEMAT2A = KEMAT1A + NMATAB(ISYMA)
      KRIM    = KEMAT2A + NMATIJ(ISYMA)
      KEND0   = KRIM    + NEMAT1(ISYMA)
      LWRK0   = LWORK   - KEND0

*----------------------------------------
* BF contribution: LAM^B x LAM^0 x BF^A:
*----------------------------------------
      IF (.NOT. CC2) THEN

        NRHO = NT2AOIJ(ISYMA)

        KGAMMA  = KEND0
        KGIM    = KGAMMA  + NGAMMA(ISYMA)
        KBFA    = KGIM    + NT1AO(ISYMA)
        KEND1   = KBFA    + NRHO
        LWRK1   = LWORK   - KEND1

        IF (LWRK1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCBMAT2. (1)')
        END IF

* read BF^A intermediate:
        DTIME = SECOND()
        CALL CC_RVEC(LUBF,BFIL,LENBF,NRHO,IINT1A,WORK(KBFA))
        TIMIO = TIMIO + SECOND() - DTIME

* initialize GAMMA^A and G^A intermediates:
        CALL DZERO(WORK(KGAMMA),NGAMMA(ISYMA))
        CALL DZERO(WORK(KGIM),  NT1AO(ISYMA))

* transform to MO representation using one response Lambda matrix (B)
* and calculate the Gamma intermediate using XLAMP0
        ICON    = 2
        IOPTG   = 2
        LGAMMA  = .TRUE.
        LO3BF   = .TRUE.
        DTIME = SECOND()
        CALL CC_T2MO3(DUM,DUM,1,WORK(KBFA),
     *                THETA2,WORK(KGAMMA),WORK(KGIM),DUM, 
     *                XLAMP0,ISYM0,WORK(KT1AMPB),ISYMB,
     *                WORK(KEND1),LWRK1,ISYMA,ICON,
     *                LGAMMA,IOPTG,LO3BF,.FALSE.)
        TIMBF = TIMBF + SECOND() - DTIME


        IF (LOCDBG) THEN
          XNORM=DDOT(NGAMMA(ISYMA),WORK(KGAMMA),1,WORK(KGAMMA),1)
          WRITE (LUPRI,*) 'Norm of response GAMMA intermediate:',XNORM
          XNORM=DDOT(NRHO,WORK(KBFA),1,WORK(KBFA),1)
          WRITE (LUPRI,*) 'Norm of response BF intermediate:',XNORM
          XNORM=DDOT(NT2AM(ISYMAB),THETA2,1,THETA2,1)
          WRITE (LUPRI,*) 'Norm of THETA2 after BF contribution:',XNORM
          WRITE (LUPRI,*) 'THETA after BF contribution:'
          Call CC_PRP(THETA1,THETA2,ISYMAB,1,1)
          CALL FLSHFO(LUPRI)
        END IF
 
      END IF

*-------------------------------------------------------------------
* for CCSD calculate here E intermed. from G, R, YBAR intermediates:
* This might miss a frozen core term!!!!!!!!
*-------------------------------------------------------------------
      IF (CCSD.OR.CCSDT) THEN
         KONEHR = KEND1
         KONEHG = KONEHR + MAX(N2BST(ISYM0),N2BST(ISYMA))
         KEND1  = KONEHG + MAX(N2BST(ISYM0),N2BST(ISYMA))
         LWRK1  = LWORK  - KEND1

         IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCBMAT2. (0)')
         END IF

         DTIME = SECOND()

         CALL CCRHS_ONEAO(WORK(KONEHR),WORK(KEND1),LWRK1)
         DO IF= 1, NFIELD
           FF = EFIELD(IF)
           CALL CC_ONEP(WORK(KONEHR),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
         END DO
C
C------------------------------------------------------------------------
C     CCMM, 03 JK+OC
C     Solvent/QMMM  contribution to one-electron integrals.
C     T^g contribution to transformation.
C------------------------------------------------------------------------
C
      IF (CCSLV) THEN
         IF (.NOT.CCMM) CALL CCSL_RHSTG(WORK(KONEHR),WORK(KEND1),LWRK1)
         IF (CCMM) THEN
             IF (.NOT. NYQMMM) THEN
               CALL CCMM_RHSTG(WORK(KONEHR),WORK(KEND1),LWRK1)
             ELSE IF (NYQMMM) THEN
               IF (HFFLD) THEN
                 WRITE(LUPRI,*) 'Is it justified to do B transfor with
     &  with a HFFLD?'
                 CALL QUIT('HFFLD not implemented for QR')
               ELSE
                 CALL CCMM_ADDG(WORK(KONEHR),WORK(KEND1),LWRK1)
               END IF
             END IF
         END IF
      ENDIF
C
C------------------------------------------------------------------------
C
         IF (FROIMP.OR.FROEXP) THEN
           CALL DAXPY(N2BST(ISYM0),ONE,FCKC0,1,WORK(KONEHR),1)
         END IF
         CALL DCOPY(N2BST(ISYM0),WORK(KONEHR),1,WORK(KONEHG),1)

         CALL CC_FCKMO(WORK(KONEHR),XLAMPA,XLAMH0,WORK(KEND1),LWRK1,
     *                 ISYM0,ISYMA,ISYM0)
         CALL CC_FCKMO(WORK(KONEHG),XLAMP0,XLAMHA,WORK(KEND1),LWRK1,
     *                 ISYM0,ISYM0,ISYMA)

         TIME = TIME + SECOND() - DTIME

*        read R^A intermediate:
         DTIME = SECOND()
         CALL CC_RVEC(LUR,RFIL,LENR,NEMAT1(ISYMA),IINT1A,WORK(KRIM))
         TIMIO = TIMIO + SECOND() - DTIME

         DTIME = SECOND()

*        transform AO indizes of R^A and G^A intermediates and add
*        one-electron hamiltonian contributions:
         LRCON  = .TRUE.
         LGCON  = .TRUE.
         FCKCON = .TRUE.
         IOPT   = 1
         CALL CC_EIM(WORK(KEMAT1A),WORK(KEMAT2A),
     *               WORK(KRIM),DUM,WORK(KGIM),DUM,
     *               WORK(KONEHR),WORK(KONEHG),
     *               XLAMH0,XLAMP0,ISYM0,DUM,DUM,IDUM,
     *               FCKCON,LRCON,LGCON,.FALSE.,IOPT,ISYMA)

*        add T2^A contribution to E1 intermediate:
         CALL DAXPY(NMATAB(ISYMA),ONE,YBARA,1,WORK(KEMAT1A),1)

         TIME = TIME + SECOND() - DTIME

         IF (LOCDBG) THEN
            XNORM = DDOT(NEMAT1(ISYMA),WORK(KRIM),1,WORK(KRIM),1)
            WRITE (LUPRI,*) 'Norm^2 of RIM:',XNORM
            CALL AROUND( 'CCSD implementation of E2 intermediate:')
            CALL CC_PREI(WORK(KEMAT1A),WORK(KEMAT2A),ISYMA,1)
            CALL AROUND( 'ONEHR intermediate:')
            CALL CC_PRFCKMO(WORK(KONEHR),ISYMA)
            CALL AROUND( 'ONEHG intermediate:')
            CALL CC_PRFCKMO(WORK(KONEHG),ISYMA)
         END IF

      END IF

*-------------------------------------------------------------------
* read double excitation response T^B amplitudes and square them up:
*-------------------------------------------------------------------
      KT2AMPB = KEND0
      KEND0   = KT2AMPB + NT2SQ(ISYMB)
      LWRK0   = LWORK - KEND0

      KGAMMAX = KGAMMA

      KGAMMA  = KEND0
      KEND1   = KGAMMA  + NGAMMA(ISYMA)
      LWRK1   = LWORK - KEND1

      KSCR = KEND0
      IF (.NOT. CC2) THEN ! take care of Gamma intermediate,
         KSCR = KEND1     ! shift it behind the T^B amplitudes
         DO I = NGAMMA(ISYMA), 1, -1
           WORK(KGAMMA-1+I) = WORK(KGAMMAX-1+I)
         END DO
      END IF
      KEND2 = KSCR + NT2AM(ISYMB)
      LWRK2 = LWORK - KEND2

      IF (LWRK2 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCBMAT2. (3)')
      END IF

      DTIME = SECOND()

      IOPT = 2
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KSCR))

      CALL CCLR_DIASCL(WORK(KSCR),TWO,ISYMB)

      CALL CC_T2SQ(WORK(KSCR),WORK(KT2AMPB),ISYMB)

      TIMIO = TIMIO + SECOND() - DTIME

*------------------------------------------------------
* A contribution: T^B x GAMMA^A, (requires squared T2^B)
*------------------------------------------------------
      IF (.NOT. CC2) THEN
        IOPT = 1
        DTIME = SECOND()

*       add E2 intermediate contribution to diagonal of GAMMA:
        CALL CC_GAMMA2(WORK(KGAMMA),WORK(KEMAT2A),ISYMA)

*       calculate A term contribution:
        CALL CCRHS_A(THETA2,WORK(KT2AMPB),WORK(KGAMMA),
     &               WORK(KEND1),LWRK1,ISYMA,ISYMB,IOPT)

        TIMA = TIMA + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(NGAMMA(ISYMA),WORK(KGAMMA),1,WORK(KGAMMA),1)
          WRITE (LUPRI,*) 'Norm of response GAMMA intermediated:',XNORM
          XNORM=DDOT(NT2AM(ISYMAB),THETA2,1,THETA2,1)
          WRITE (LUPRI,*) 'Norm of THETA2 after A contribution:',XNORM
          CALL FLSHFO(LUPRI)
        END IF

      END IF

*--------------------------------------------------------------
* E1 & E2 contributions, T^B x E1/E2^A: (requires squared T2^B)
*--------------------------------------------------------------

      DTIME = SECOND()

      IF(CC2.AND.(.NOT.CCSTST).AND.
     &   ((NFIELD.GT.0) .OR. CCSLV)) THEN
 
         KPERT   = KEND0
         KT1AMPA = KPERT   + N2BST(ISYM0)
         KOV     = KT1AMPA + NT1AM(ISYMA)
         KEND1   = KOV     + NT1AM(ISYM0)
         LWRK1   = LWORK   - KEND1

         IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCBMAT2. (4)')
         END IF
        
         ! read finite field perturb. integrals and transform to MO
         CALL DZERO(WORK(KPERT),N2BST(ISYM0))
         DO IF = 1, NFIELD
            FF = EFIELD(IF)
            CALL CC_ONEP(WORK(KPERT),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
         END DO
C
C----------------------------------------------------------------------
C        CCMM, 03 JK+OC
C        Solvent/QMMM  contribution to one-electron integrals.
C        T^g contribution to transformation.
C----------------------------------------------------------------------
C
         IF (CCSLV) THEN
           IF (.NOT.CCMM) CALL CCSL_RHSTG(WORK(KPERT),WORK(KEND1),LWRK1)
           IF (CCMM) THEN
               IF(.NOT. NYQMMM) THEN
              CALL CCMM_RHSTG(WORK(KPERT),WORK(KEND1),LWRK1)
               ELSE IF (NYQMMM) THEN
                 IF (HFFLD) THEN
                   WRITE(LUPRI,*) 'Is it justified to do B transfor with
     &  with a HFFLD?'
                   CALL QUIT('HFFLD not implemented for QR')
                 ELSE
                   CALL CCMM_ADDG(WORK(KPERT),WORK(KEND1),LWRK1)
                 END IF
               END IF  
           END IF
         ENDIF
C
C------------------------------------------------------------------------
C
         CALL CC_FCKMO(WORK(KPERT),XLAMP0,XLAMH0,WORK(KEND1),LWRK1,
     &                 ISYM0,1,1)

         ! gather occupied/virtual block and calculate one-index 
         ! transformed V^A (occ/occ and vir/vir blocks)
         CALL CC_GATHEROV(WORK(KPERT),WORK(KOV),ISYM0)
         IOPT = 1
         CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &                 WORK(KT1AMPA),WORK(KDUM))
         CALL CCG_1ITROO(WORK(KEMAT2A),ISYMA,
     &                   WORK(KOV), ISYM0, WORK(KT1AMPA),ISYMA)
         CALL CCG_1ITRVV(WORK(KEMAT1A),ISYMA, 
     &                   WORK(KOV), ISYM0, WORK(KT1AMPA),ISYMA)
       
*        calculate the contribution to THETA2:
         CALL CCRHS_E(THETA2,WORK(KT2AMPB),WORK(KEMAT1A),
     &                WORK(KEMAT2A),WORK(KEND1),LWRK1,ISYMB,ISYMA)
 
      END IF

      TIME = TIME + SECOND() - DTIME

      IF (LOCDBG .AND. .NOT.(CCSD.OR.CCSDT)) THEN
        XNORM=DDOT(NT2SQ(ISYMB),WORK(KT2AMPB),1,WORK(KT2AMPB),1)
        WRITE (LUPRI,*) 'Norm of response T2AMPB amplitudes:',XNORM
        XNORM=DDOT(NMATAB(ISYMA),WORK(KEMAT1A),1,WORK(KEMAT1A),1)
        WRITE (LUPRI,*) 'Norm of response EMAT1A intermediated:',XNORM
        XNORM=DDOT(NMATIJ(ISYMA),WORK(KEMAT2A),1,WORK(KEMAT2A),1)
        WRITE (LUPRI,*) 'Norm of response EMAT2A intermediated:',XNORM
        XNORM=DDOT(NT2AM(ISYMAB),THETA2,1,THETA2,1)
        WRITE (LUPRI,*) 'Norm of THETA2 after E contribution:',XNORM
        WRITE (LUPRI,*) 'THETA2 after E contribution:'
        Call CC_PRP(THETA1,THETA2,ISYMAB,1,1)
        WRITE (LUPRI,*) 'ISYMA, ISYMB, ISYMAB:',ISYMA,ISYMB,ISYMAB
        CALL AROUND( 'response E-intermediates (1-el part)')
        CALL CC_PREI(FCKAVV,FCKAOO,ISYMA,1)
        CALL AROUND( 'response E-intermediates')
        CALL CC_PREI(WORK(KEMAT1A),WORK(KEMAT2A),ISYMA,1)
        CALL FLSHFO(LUPRI)
      END IF

*----------------------------------------------------------------
* C contribution: T^B x (CTILDE^A + CBAR^A)  or with CTILDE only:
*----------------------------------------------------------------
      IF (.NOT. CC2) THEN

* transpose occupied indeces of the amplitudes:
        CALL CCSD_T2TP(WORK(KT2AMPB),WORK(KEND0),LWRK0,ISYMB)

        IOPT  = 2
        IOPTE = 1
        FACB  = ONE
        DTIME = SECOND()
        CALL CCRHS_CIO2(THETA2,WORK(KT2AMPB),XLAMH0,
     &                  WORK(KEND0),LWRK0,ISYMB,ISYMA,
     &                  LUC,CTFIL,IINT1A,IOPT,
     &                  IOPTB,LUCBAR,CBAFIL,IOFFCDB,FACB,
     &                  IOPTE,WORK(KEMAT1A),.FALSE.) 
        TIMC = TIMC + SECOND() - DTIME
 
* restore original amplitudes:
        CALL CCSD_T2TP(WORK(KT2AMPB),WORK(KEND0),LWRK0,ISYMB)
 
        IF (LOCDBG) THEN
          XNORM=DDOT(NT2AM(ISYMAB),THETA2,1,THETA2,1)
          WRITE (LUPRI,*) 'Norm of THETA2 after C contribution:',XNORM
          CALL FLSHFO(LUPRI)
        END IF

      END IF


*----------------------------------------------------------------
* calculate  2*t^B(iajb) - t^B(ibja) in place of the T2^B vector:
*----------------------------------------------------------------
      CALL CCRHS_T2TR(WORK(KT2AMPB),WORK(KEND0),LWRK0,ISYMB)

*----------------------------------
* I contribution, (2T^B-T^B) x F^A
*----------------------------------
      DTIME = SECOND()

      IOPT = 1
      CALL CCG_LXD(WORK(KCON),ISYMAB,FOCKA,ISYMA,
     &             WORK(KT2AMPB),ISYMB,IOPT)
      CALL DAXPY(NT1AM(ISYMAB),ONE,WORK(KCON),1,THETA1,1)

      TIMI = TIMI + SECOND() - DTIME

      IF (LOCDBG) THEN
        XNORM=DDOT(NT1AM(ISYMAB),THETA1,1,THETA1,1)
        WRITE (LUPRI,*) 'Norm of THETA1 after I contribution:',XNORM
        WRITE (LUPRI,*) 'THETA1 after I contribution:'
        Call CC_PRP(THETA1,THETA2,ISYMAB,1,0)
        CALL FLSHFO(LUPRI)
      END IF

*----------------------------------------------------------------------
* D contribution: (2T^B-T^B) x (DTILDE^A + DBAR^A) or with DTILDE only:
*----------------------------------------------------------------------
      IF (.NOT. CC2) THEN

        IOPT  = 2
        IOPTE = 1
        FACB  = ONE
        DTIME = SECOND()
        CALL CCRHS_DIO2(THETA2,WORK(KT2AMPB),XLAMH0,
     &                  WORK(KEND0),LWRK0,ISYMB,ISYMA,
     &                  LUD,DTFIL,LUC,CTFIL,IINT1A,IOPT,
     &                  IOPTB,LUDBAR,DBAFIL,IOFFCDB,FACB, 
     &                  IOPTE,WORK(KEMAT1A),.FALSE.) 
        TIMD = TIMD + SECOND() - DTIME

        IF (LOCDBG) THEN
          XNORM=DDOT(NT2AM(ISYMAB),THETA2,1,THETA2,1)
          WRITE (LUPRI,*) 'Norm of THETA2 after D contribution:',XNORM
          CALL FLSHFO(LUPRI)
        END IF

      END IF

*--------------------------
* J, G and H contributions:
*--------------------------
      IF (CCSD.OR.CCSDT) THEN
         KT1AMPA = KEND0
         KFCKA   = KT1AMPA + NT1AM(ISYMA)
         KEND0   = KFCKA   + NMATAB(ISYMA)
         LWRK0   = LWORK   - KEND0

         IF (LWRK0 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCBMAT2. (J)')
         END IF

         DTIME = SECOND()
         
         IOPT = 1
         CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &                 WORK(KT1AMPA),WORK(KDUM))

         CALL CCG_1ITRVV(WORK(KFCKA),ISYMA,FCK0OV,ISYM0,
     &                   WORK(KT1AMPA),ISYMA)

         CALL DAXPY(NMATAB(ISYMA),-ONE,WORK(KFCKA),1,WORK(KEMAT1A),1)

      ELSE
         CALL DCOPY(NMATAB(ISYMA),    FOCKVV,1,WORK(KEMAT1A),1)
         CALL DAXPY(NMATAB(ISYMA),ONE,YBARA, 1,WORK(KEMAT1A),1)

         CALL DCOPY(NMATIJ(ISYMA),    FCKAOO,1,WORK(KEMAT2A),1)
         CALL DAXPY(NMATIJ(ISYMA),ONE,XBARA, 1,WORK(KEMAT2A),1)
      END IF

      IOPT = 1
      CALL CCG_1ITRVO(WORK(KCON),ISYMAB,
     &                WORK(KEMAT2A),WORK(KEMAT1A),ISYMA,
     &                WORK(KT1AMPB),ISYMB,IOPT)

      CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KCON),1,THETA1,1)

      TIMI = TIMI + SECOND() - DTIME

      IF (LOCDBG) THEN
         CALL AROUND( 'Intermediates for J, G & H terms:')
         CALL CC_PREI(WORK(KEMAT1A),WORK(KEMAT2A),ISYMA,1)
         WRITE (LUPRI,*) 'THETA1 after J, G & H contributions:'
         Call CC_PRP(THETA1,THETA2,ISYMAB,1,0)
      END IF

*----------------------------------------------------------------------

      CALL QEXIT('CCBMAT2')

      RETURN

      END

*=====================================================================*
*            END OF SUBROUTINE CCBMAT2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBPRE1 */
*=====================================================================*
      SUBROUTINE CCBPRE1(INTMED1,ISTART,IEND,
     &                   KRHO2,KLAMP,KLAMH,KDENS,KFOCK,KRIM,
     &                   LUBF,BFFIL,LENBF,LUFK,FKFIL,LENFK,
     &                   LUR,RFIL,LENR,
     &                   XLAMDP,XLAMDH,WORK,LWORK,KENDIN,KENDOUT)
*---------------------------------------------------------------------*
*    Purpose: prepare for calculation of intermediates that depend
*             on the AO integrals and one response vector
*
*    N.B.: this routine allocates work space for CC_BMAT 
*             INPUT  end of used space: KENDIN
*             OUTPUT end of used space: KENDOUT
*      
*     Written by Christof Haettig, Januar/Februar 1997.
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccorb.h"
#include "cciccset.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER KDUM 
      PARAMETER (KDUM = -99 999 999)  ! dummy address

      INTEGER LWORK, KENDIN, KENDOUT
      INTEGER ISTART, IEND
      INTEGER LUBF, LENBF, LUFK, LENFK, LUR, LENR
      CHARACTER*(*) BFFIL, FKFIL, RFIL
      INTEGER INTMED1(2,IEND)
      INTEGER KLAMP(IEND), KLAMH(IEND)
      INTEGER KRHO2(IEND), KDENS(IEND), KFOCK(IEND), KRIM(IEND)

      CHARACTER*(10) MODEL
      CHARACTER*(3) LIST
      INTEGER IOPT, IC, ISYM, IDX, IDLST
      INTEGER LEN, NRHO, KT1AMP
      
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XLAMDP(NLAMDT), XLAMDH(NLAMDT)
      REAL XNORM, DDOT
      REAL TWO
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XLAMDP(NLAMDT), XLAMDH(NLAMDT)
      DOUBLE PRECISION XNORM, DDOT
      DOUBLE PRECISION TWO
#endif
      PARAMETER (TWO = 2.0d0)

* external functions:
      INTEGER ILSTSYM

      CALL QENTER('CCBPRE1')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      KENDOUT = KENDIN
           
      DO IDX = ISTART, IEND
        LIST  = VTABLE(INTMED1(2,IDX))
        IDLST = INTMED1(1,IDX)
        ISYM  = ILSTSYM(LIST,IDLST)
        NRHO  = NT2AOIJ(ISYM)

* memory allocation:
        IF (CCS .OR. CC2) THEN
          KRIM(IDX)   = KDUM
          KRHO2(IDX)  = KDUM
          KLAMP(IDX)  = KENDOUT
          KLAMH(IDX)  = KLAMP(IDX) + NGLMDT(ISYM)
          KDENS(IDX)  = KLAMH(IDX) + NGLMDT(ISYM)
          KFOCK(IDX)  = KDENS(IDX) + N2BST(ISYM)
          KENDOUT     = KFOCK(IDX) + N2BST(ISYM)
        ELSE IF (CCSD.OR.CCSDT) THEN
          KRIM(IDX)   = KENDOUT
          KRHO2(IDX)  = KRIM(IDX)  + NEMAT1(ISYM)
          KLAMP(IDX)  = KRHO2(IDX) + NRHO
          KLAMH(IDX)  = KLAMP(IDX) + NGLMDT(ISYM)
          KENDOUT     = KLAMH(IDX) + NGLMDT(ISYM)
          KDENS(IDX)  = KDUM
          KFOCK(IDX)  = KDUM
        ELSE
          CALL QUIT('Unknown CC model in CCBPRE1.')
        END IF

        IF ( (LWORK-KENDOUT) .LE. NT1AM(ISYM) ) THEN
          CALL QUIT('Insufficient work space in CCBPRE1.')
        END IF

* read singles part of the response vector and
* calculate response Lambda matrices:
        KT1AMP = KENDOUT
        IOPT   = 1 
        CALL CC_RDRSP(LIST,IDLST,ISYM,IOPT,MODEL,
     &                WORK(KT1AMP),WORK(KENDOUT))

        CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMP(IDX)),
     &                   XLAMDH,WORK(KLAMH(IDX)),
     &                   WORK(KT1AMP),ISYM)

* calculate response density matrix:
        IF (.NOT.(CCSD.OR.CCSDT)) THEN
          IC = 0 ! no core contribution
          CALL CC_AODENS(XLAMDP,WORK(KLAMH(IDX)),WORK(KDENS(IDX)),ISYM,
     &                 IC,WORK(KENDOUT),LWORK-KENDOUT)
        END IF

* recover the R and BF intermediates:
        IF (CCSD.OR.CCSDT) THEN

          CALL CC_RVEC(LUR,RFIL,LENR,NEMAT1(ISYM),IDX,WORK(KRIM(IDX)))

          CALL CC_RVEC(LUBF,BFFIL,LENBF,NRHO,IDX,WORK(KRHO2(IDX)))

        END IF

* recover the response Fock matrix:
        IF (.NOT.(CCSD.OR.CCSDT)) THEN
           LEN = N2BST(ISYM)
           CALL CC_RVEC(LUFK,FKFIL,LENFK,LEN,IDX,WORK(KFOCK(IDX)))
        END IF

        IF (LOCDBG) THEN
          XNORM = DDOT(NGLMDT(ISYM),
     &                 WORK(KLAMP(IDX)),1,WORK(KLAMP(IDX)),1)
          WRITE (LUPRI,*) 'Norm of response LAMDP nb. ', IDX, ' is ',
     &         XNORM
          XNORM = DDOT(NGLMDT(ISYM),
     &                 WORK(KLAMH(IDX)),1,WORK(KLAMH(IDX)),1)
          WRITE (LUPRI,*) 'Norm of response LAMDH nb. ', IDX, ' is ', 
     &         XNORM
          IF (.NOT.(CCSD.OR.CCSDT)) THEN
            XNORM = DDOT(LEN,WORK(KDENS(IDX)),1,WORK(KDENS(IDX)),1)
            WRITE (LUPRI,*) 'Norm of response DENSITY nb. ', IDX,
     &           ' is ',XNORM
            XNORM = DDOT(LEN,WORK(KFOCK(IDX)),1,WORK(KFOCK(IDX)),1)
            WRITE (LUPRI,*) 'Norm of recovered FOCK nb. ', IDX, ' is ',
     &           XNORM
          ELSE IF (CCSD.OR.CCSDT) THEN
            XNORM = DDOT(NRHO,WORK(KRHO2(IDX)),1,WORK(KRHO2(IDX)),1)
            WRITE (LUPRI,*) 'Norm of recovered BF nb. ', IDX, ' is ',
     &           XNORM
          END IF
        END IF

      END DO

      CALL QEXIT('CCBPRE1')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCBPRE1
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBPRE2 */
*=====================================================================*
      SUBROUTINE CCBPRE2(INTMED2,ISTART,IEND,LUF,FFIL,LENF,
     &                   KOMEGA2,KLAMPA,KLAMHA,KLAMPB,KLAMHB,
     &                   XLAMDP,XLAMDH,WORK,LWORK,KENDIN,KENDOUT )
*---------------------------------------------------------------------*
*    Purpose: prepare for calculation of intermediates that depend
*             on the AO integrals and two response vector
*
*    N.B.: this routine allocates work space for CC_BMAT 
*             INPUT  end of used space: KENDIN
*             OUTPUT end of used space: KENDOUT
*      
*     Written by Christof Haettig, Januar/Februar 1997.
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "cciccset.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER KDUM
      PARAMETER (KDUM = +99 999 999) ! dummy address on work space

      INTEGER LWORK, KENDIN, KENDOUT
      INTEGER ISTART, IEND
      INTEGER LUF, LENF
      INTEGER INTMED2(4,IEND)
      INTEGER KLAMPA(IEND), KLAMHA(IEND)
      INTEGER KLAMPB(IEND), KLAMHB(IEND)
      INTEGER KOMEGA2(IEND)

      CHARACTER*(*) FFIL
      CHARACTER*(3) LISTA, LISTB
      CHARACTER*(10) MODEL
      INTEGER KT1AMPA, KT1AMPB
      INTEGER LEN, KEND1, IOPT
      INTEGER IDLSTA, IDLSTB, ISYMA, ISYMB, ISYMAB, IINT2
      
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XLAMDP(NLAMDT), XLAMDH(NLAMDT)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XLAMDP(NLAMDT), XLAMDH(NLAMDT)
#endif

* external functions:
      INTEGER ILSTSYM

      CALL QENTER('CCBPRE2')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      KENDOUT = KENDIN
           
      DO IINT2 = ISTART, IEND
        LISTA  = VTABLE(INTMED2(2,IINT2))
        LISTB  = VTABLE(INTMED2(4,IINT2))
        IDLSTA = INTMED2(1,IINT2)
        IDLSTB = INTMED2(3,IINT2)
        ISYMA  = ILSTSYM(LISTA,IDLSTA)
        ISYMB  = ILSTSYM(LISTB,IDLSTB)
        ISYMAB = MULD2H(ISYMA,ISYMB)

        KOMEGA2(IINT2)  = KDUM

        IF (CCS) THEN
          KLAMPA(IINT2)   = KDUM
          KLAMHA(IINT2)   = KDUM
          KLAMPB(IINT2)   = KDUM
          KLAMHB(IINT2)   = KDUM
          KENDOUT         = KENDOUT
        ELSE 
          IF (CC2) THEN
           KOMEGA2(IINT2) = KENDOUT
           KENDOUT        = KOMEGA2(IINT2)  + NT2AM(ISYMAB)
          END IF

          KLAMPA(IINT2)   = KENDOUT
          KLAMHA(IINT2)   = KLAMPA(IINT2)   + NGLMDT(ISYMA)
          KLAMPB(IINT2)   = KLAMHA(IINT2)   + NGLMDT(ISYMA)
          KLAMHB(IINT2)   = KLAMPB(IINT2)   + NGLMDT(ISYMB)
          KENDOUT         = KLAMHB(IINT2)   + NGLMDT(ISYMB)

          KT1AMPA = KENDOUT 
          KT1AMPB = KT1AMPA + NT1AM(ISYMA)
          KEND1   = KT1AMPB + NT1AM(ISYMB)

          IF ( (LWORK-KEND1) .LE. 0 ) THEN
            CALL QUIT('Insufficient work space in CCBPRE2.')
          END IF

* recover F intermediate:
          IF (CC2) THEN
            LEN = NT2AM(ISYMAB)
            CALL CC_RVEC(LUF,FFIL,LENF,LEN,IINT2,WORK(KOMEGA2(IINT2)))
          END IF

* A response Lambda matrices:
          IOPT = 1 ! read singles response vector
          CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &                  WORK(KT1AMPA),WORK(KDUM) )

          ! calculate response Lambda matrices:
          CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPA(IINT2)),
     &                     XLAMDH,WORK(KLAMHA(IINT2)),
     &                     WORK(KT1AMPA),ISYMA)

* B response Lambda matrices:
          IOPT = 1 ! read singles response vector
          CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &                  WORK(KT1AMPB),WORK(KDUM) )

          ! calculate response Lambda matrices:
          CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPB(IINT2)),
     &                     XLAMDH,WORK(KLAMHB(IINT2)),
     &                     WORK(KT1AMPB),ISYMB)

        END IF
      END DO

      CALL QEXIT('CCBPRE2')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCBPRE2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBINT1 */
*=====================================================================*
       SUBROUTINE CCBINT1(XINT,  BDRHF, DCRHF,
     &                    IDEL, ISYMD, RHO2A,
     &                    XLAMP0, XLAMH0, XLAMPA, XLAMHA,
     &                    ISYMA,  IVECA,  DENSA,  FOCKA, RIMA,
     &                    LUC,    CTFIL,  LUD,    DTFIL, 
     &                    LUBFD,  FNBFD,  IADRBFD, WORK,   LWORK,
     &                    TIMFCK, TIMBF, TIMC, TIMD  )
*---------------------------------------------------------------------*
*
*    Purpose: calculate intermediates for B matrix transformation
*             which require AO integrals and depend on one response
*             vector:
*
*               RHO^{BF,A}, Ctilde^A, Dtilde^A, Ftilde^{A,*}, R^A
*
*     input:  XLAMPA - response A Lambda particle matrix
*             XLAMHA - response A Lambda hole matrix
*             XLAMP0 - zeroth order Lambda particle matrix
*             XLAMH0 - zeroth order Lambda hole matrix
*             DENSA  - A response of the density matrix
*
*    output:  RHO2A - updated RHO^{BF,A} intermediate
*             FOCKA - updated Ftilde^{A,*} intermediate
*             RIMA  - updated R^{A,*} intermediate
*
*    written to file:  contributions to Ctilde^A and Dtilde^A intermed.
*
*
*    Written by Christof Haettig, Januar/Februar 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "second.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
      
      CHARACTER*(*) CTFIL, DTFIL, FNBFD
      INTEGER LWORK, IDEL, ISYMD, ISYMA
      INTEGER LUC, LUD, LUBFD, IVECA
      INTEGER IADRBFD(*)

      INTEGER KEND0, KEND1, LWRK0, LWRK1, KSCRCM, IOPT, ISYMM
      INTEGER IADR, KMGD, NMGD, NRHO, IOPTR

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XINT(*), BDRHF(*), DCRHF(*)
      REAL XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)   
      REAL RHO2A(*), DENSA(*), FOCKA(*), RIMA(*)      

      REAL FACTC, FACTD, DUMMY, XNORM
      REAL TIMFCK, TIMBF, TIMC, TIMD, DTIME
      REAL DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XINT(*), BDRHF(*), DCRHF(*)
      DOUBLE PRECISION XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)  
      DOUBLE PRECISION RHO2A(*), DENSA(*), FOCKA(*), RIMA(*)   

      DOUBLE PRECISION FACTC, FACTR, DUMMY, XNORM
      DOUBLE PRECISION TIMFCK, TIMBF, TIMC, TIMD, DTIME
      DOUBLE PRECISION DDOT
#endif


      CALL QENTER('CCBINT1')

  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
       KEND0 = 1
       LWRK0 = LWORK

*---------------------------------------------------------------------*
* Ftilde^{A,*}:
*---------------------------------------------------------------------*
       IF (.NOT. (CCSD.OR.CCSDT)) THEN

cch
c            XNORM=DDOT(N2BST(ISYMA),FOCKA,1,FOCKA,1)
c       WRITE (LUPRI,*) 'CCBINT1> norm of FOCKA matrix (before):',XNORM
cch
          DTIME = SECOND()
          CALL CC_AOFOCK( XINT, DENSA, FOCKA,
     *                    WORK(KEND0),LWRK0,IDEL,ISYMD,.FALSE.,DUMMY,
     *                    ISYMA )
          TIMFCK = TIMFCK + SECOND() - DTIME
 
cch       IF (LOCDBG) THEN
c            WRITE(LUPRI,*) 'ISYMD, ISYMA:',ISYMD, ISYMA
c            XNORM = DDOT(NDISAO(ISYMD),XINT,1,XINT,1)
c            WRITE (LUPRI,*) 'CCBINT1> norm of XINT matrix:',XNORM
c            XNORM=DDOT(N2BST(ISYMA),DENSA,1,DENSA,1)
c            WRITE (LUPRI,*) 'CCBINT1> norm of DENSA matrix:',XNORM
c            XNORM=DDOT(N2BST(ISYMA),FOCKA,1,FOCKA,1)
c            WRITE (LUPRI,*) 'CCBINT1> norm of FOCKA matrix:',XNORM
cch       END IF
 
       END IF

*---------------------------------------------------------------------*
* RHO^{BF,A}:
*---------------------------------------------------------------------*
       IF (.NOT. (CCS .OR. CC2)) THEN
         DTIME = SECOND()

         ISYMM = MULD2H(ISYMD,ISYMA)
         NMGD  = NT2BGD(ISYMM)

         KMGD  = KEND0
         KEND1 = KMGD  + NMGD
         LWRK1 = LWORK - KEND1

         IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCBINT1.')
         END IF

*        read delta batch of the effective density:
         IADR = IADRBFD(IDEL)
         CALL GETWA2(LUBFD,FNBFD,WORK(KMGD),IADR,NMGD)

*        update BF intermediate:
         CALL CC_BFIB(RHO2A,BDRHF,ISYMD,WORK(KMGD),ISYMM,
     *                WORK(KEND1),LWRK1)

         TIMBF = TIMBF + SECOND() - DTIME

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCBINT1> IDEL, ISYMD:',IDEL,ISYMD
            XNORM=DDOT(NT2AOIJ(ISYMA),RHO2A,1,RHO2A,1)
            WRITE (LUPRI,*) 'CCBINT1> norm of RHO2A:',XNORM
         END IF

       END IF

*---------------------------------------------------------------------*
* Ctilde^A, Dtilde^A, and R^A:
*---------------------------------------------------------------------*
       IF (.NOT.(CCS .OR. CC2)) THEN
 
          IOPTR = 1
          FACTR = -2.0D0
          DTIME = SECOND()
          CALL CC_CDB(DCRHF, ISYMD, IDEL, ISYMD, LUC, CTFIL, IVECA,
     *                XLAMP0, XLAMH0, XLAMPA, XLAMHA, ISYMA,
     *                IOPTR, FACTR, RIMA, WORK(KEND0), LWRK0 )
          TIMC = TIMC + SECOND() - DTIME


          IOPTR = 1
          FACTR = 1.0D0
          DTIME = SECOND()
          CALL CC_CDB(BDRHF, ISYMD, IDEL, ISYMD, LUD, DTFIL, IVECA,
     *                XLAMP0, XLAMH0, XLAMPA, XLAMHA, ISYMA,
     *                IOPTR, FACTR, RIMA, WORK(KEND0), LWRK0 )
          TIMD = TIMD + SECOND() - DTIME

       ENDIF

*---------------------------------------------------------------------*

       CALL QEXIT('CCBINT1')

       RETURN
       END
*=====================================================================*
*            END OF SUBROUTINE CCBINT1
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBINT2 */
*=====================================================================*
       SUBROUTINE CCBINT2(XINT, IDEL, ISYMD, OMEGA2, ISYOMEG,
     &                    LUAIBJ, FNAIBJ, IT2F, IADRF, NEWFTERM,
     &                    XLAMPA, XLAMHA, ISYMA,
     &                    XLAMPB, XLAMHB, ISYMB,
     &                    XLAMP0, XLAMH0, WORK,   LWORK   )
*---------------------------------------------------------------------*
*
*    Purpose: calculate intermediates for B matrix transformation
*             which require AO integrals and depend on two response
*             vectors:
*
*               
*     input:  XINT           - AO integral distribution
*             IDEL           - delta index of XINT
*             XLAMP0, XLAMH0 - ordinary zeroth order Lambda matrices
*             XLAMPA, XLAMHA - response A Lambda matrices
*             XLAMPB, XLAMHB - response B Lambda matrices
*
*     output: F term contribution to B matrix in MO representation
*
*
*     Written by Christof Haettig, Januar/Februar 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
      
      LOGICAL NEWFTERM
      CHARACTER*(*) FNAIBJ
      INTEGER IT2F(*)
      INTEGER LUAIBJ, LWORK, IDEL, ISYMD, ISYMA, ISYMB, ISYOMEG
      INTEGER IOPT, ISYALBE, IGAM, KXAIBJ, LEN, ISYIAJ, IADRF
      INTEGER ISYGAM, KEND1, LWRK1, IDUMMY, KOFF


#if defined (SYS_CRAY)
      REAL WORK(LWORK), XINT(*), OMEGA2(*)
      REAL XLAMP0(*), XLAMH0(*)
      REAL XLAMPA(*), XLAMHA(*), XLAMPB(*), XLAMHB(*)    
      REAL DUMMY
#else
      DOUBLE PRECISION WORK(LWORK), XINT(*), OMEGA2(*)
      DOUBLE PRECISION XLAMP0(*), XLAMH0(*)
      DOUBLE PRECISION XLAMPA(*), XLAMHA(*), XLAMPB(*), XLAMHB(*)    
      DOUBLE PRECISION DUMMY
#endif


      CALL QENTER('CCBINT2')

  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (CCS) THEN
        CONTINUE

      ELSE IF (CC2) THEN
*---------------------------------------------------------------------*
* for CC2 calculate the complete F term:
*---------------------------------------------------------------------*
        IOPT = 3

        CALL CC_MOFCON2(XINT,OMEGA2,
     &                  XLAMPA,XLAMHA,XLAMPB,XLAMHB,
     &                  XLAMP0,XLAMH0,ISYMA,ISYMB,ISYM0,ISYM0,
     &                  WORK,LWORK,IDEL,ISYMD,ISYOMEG,ISYM0,IOPT)

        CALL CC_MOFCON2(XINT,OMEGA2,
     &                  XLAMPB,XLAMHB,XLAMP0,XLAMH0,
     &                  XLAMP0,XLAMHA,ISYMB,ISYM0,ISYM0,ISYMA,
     &                  WORK,LWORK,IDEL,ISYMD,ISYOMEG,ISYM0,IOPT)

        CALL CC_MOFCON2(XINT,OMEGA2,
     &                  XLAMPB,XLAMHB,XLAMP0,XLAMH0,
     &                  XLAMPA,XLAMH0,ISYMB,ISYM0,ISYMA,ISYM0,
     &                  WORK,LWORK,IDEL,ISYMD,ISYOMEG,ISYM0,IOPT)

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'DEBUG_CCBINT2> used CC2 version of F term.'
        END IF

      ELSE IF (CCSD.OR.CCSDT) THEN
*---------------------------------------------------------------------*
* for CCSD calculate only (a i^A | b j^B) + (a i^B | b j^A) :
*---------------------------------------------------------------------*
        ISYIAJ = MULD2H(ISYMD,MULD2H(ISYMA,ISYMB))
        LEN    = NT2BCD(ISYIAJ)

        KXAIBJ = 1
        KEND1  = KXAIBJ + LEN
        LWRK1  = LWORK  - KEND1

        IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCBINT2.')
        END IF

        CALL DZERO(WORK(KXAIBJ),LEN)
           
        DO ISYGAM = 1, NSYM
        DO G = 1, NBAS(ISYGAM)

           ISYALBE = MULD2H(ISYMD,ISYGAM)

           IGAM = G + IBAS(ISYGAM)

           KOFF = IDSAOG(ISYGAM,ISYMD) + NNBST(ISYALBE)*(G-1) + 1

           IOPT = 0
           CALL CC_AIBJ( XINT(KOFF), ISYALBE, DUMMY, IDUMMY,
     &                   IDEL, IGAM,    WORK(KXAIBJ), DUMMY, 
     &                   XLAMHA,ISYMA,  XLAMHB,ISYMB,
     &                   XLAMP0,ISYM0,  WORK(KEND1),  LWRK1,
     &                   IOPT,  .FALSE., .FALSE.             )

        END DO
        END DO

        CALL PUTWA2(LUAIBJ, FNAIBJ, WORK(KXAIBJ), IADRF, LEN)

        IT2F(IDEL) = IADRF
        IADRF      = IADRF + LEN

*---------------------------------------------------------------------*
      ELSE 
        CALL QUIT('Unknown Coupled Cluster model in CCBINT2.')
      END IF

      CALL QEXIT('CCBINT2')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCBINT2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBINT3 */
*=====================================================================*
       SUBROUTINE CCBINT3(LIST,IDLST,LUFK,FKFIL,LENFK,IDXFK,
     &                    KFOCK,KFOCKOO,KFOCKOV,KFOCKVV,KXBAR,KYBAR,
     &                    XLIAJB,ISYOVOV,XLAMP0,XLAMH0,
     &                    WORK,LWORK,KENDIN,KENDOUT, 
     &                    TIMFCK,TIMIO,TIME)
*---------------------------------------------------------------------*
*
*    Purpose: calculate some intermediates for B matrix transformation
*             which do NOT require AO integrals and depend on one 
*             response vector:
*
*              Ftilde^{A,*} (o/o, o/v and v/v blocks), Xbar^A, Ybar^A
*
*    N.B.: this routine allocates work space for CC_BMAT 
*             INPUT  end of used space: KENDIN
*             OUTPUT end of used space: KENDOUT
*
*     Written by Christof Haettig, Januar 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "second.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER KDUM, ISYM0
      PARAMETER ( KDUM = +99 999 999 )  ! dummy address
      PARAMETER (ISYM0 = 1) ! reference state symmetry
      
      CHARACTER*(*) FKFIL, LIST
      INTEGER LWORK, KENDIN, KENDOUT
      INTEGER LUFK, LENFK, IDXFK, IDLST, ISYOVOV

      CHARACTER*(10) MODEL
      INTEGER KFOCK, KFOCKOO, KFOCKOV, KFOCKVV, KXBAR, KYBAR
      INTEGER KEND1, LWRK1, KT1AMP, KT2AMP, LEN, ISYMA, IOPT

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XLAMP0(NLAMDT), XLAMH0(NLAMDT), XLIAJB(*)
      REAL TWO, XNORM, TIMFCK, TIMIO, TIME, DTIME
      REAL DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XLAMP0(NLAMDT), XLAMH0(NLAMDT), XLIAJB(*)
      DOUBLE PRECISION TWO, XNORM, TIMFCK, TIMIO, TIME, DTIME
      DOUBLE PRECISION DDOT
#endif
      PARAMETER (TWO = 2.0d0)

* external functions:
      INTEGER ILSTSYM

      CALL QENTER('CCBINT3')

  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYMA = ILSTSYM(LIST,IDLST)

      KFOCKOV = KENDIN
      KENDOUT = KFOCKOV + NT1AM(ISYMA)

      IF (.NOT.(CCSD.OR.CCSDT)) THEN
         KFOCK   = KENDOUT
         KFOCKOO = KFOCK   + N2BST(ISYMA)
         KFOCKVV = KFOCKOO + NMATIJ(ISYMA)
         KENDOUT = KFOCKVV + NMATAB(ISYMA)
      ELSE
         KFOCK   = KDUM
         KFOCKOO = KDUM
         KFOCKVV = KDUM
      END IF
      
      IF (CCS) THEN
         KXBAR   = KDUM
         KYBAR   = KDUM
      ELSE 
         KXBAR   = KENDOUT
         KYBAR   = KXBAR   + NMATIJ(ISYMA)
         KENDOUT = KYBAR   + NMATAB(ISYMA)
      END IF

      KT1AMP = KENDOUT
      KEND1  = KT1AMP + NT1AM(ISYMA)

      IF (.NOT. CCS) THEN
        KT2AMP  = KEND1
        KEND1   = KT2AMP   + NT2AM(ISYMA)
      END IF

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCBINT3.')
      END IF


      IF (.NOT.(CCSD.OR.CCSDT)) THEN
         DTIME  = SECOND()

*        read AO Ftilde^{A,*} matrix:
         LEN = N2BST(ISYMA)
         CALL CC_RVEC(LUFK,FKFIL,LENFK,LEN,IDXFK,WORK(KFOCK))

*        transform to MO representation:
         CALL CC_FCKMO(WORK(KFOCK),XLAMP0,XLAMH0,
     &                 WORK(KEND1),LWRK1,ISYMA,ISYM0,ISYM0)

*        resort occ/occ and vir/vir blocks:
         CALL CC_GATHEROO(WORK(KFOCK),WORK(KFOCKOO),ISYMA)
         CALL CC_GATHERVV(WORK(KFOCK),WORK(KFOCKVV),ISYMA)

         TIMFCK = TIMFCK  + SECOND() - DTIME
      END IF

* read the response A amplitudes:
      DTIME = SECOND()
      IOPT  = 1
      IF (.NOT.CCS) IOPT = 3
      CALL CC_RDRSP(LIST,IDLST,ISYMA,IOPT,MODEL,
     &              WORK(KT1AMP),WORK(KT2AMP))
      IF (.NOT.CCS) Call CCLR_DIASCL(WORK(KT2AMP),TWO,ISYMA)
      TIMIO = TIMIO  + SECOND() - DTIME

* calculate the occ/vir block of the one-index transformed Fock matrix:
      IOPT = 1
      DTIME  = SECOND()
      CALL CCG_LXD(WORK(KFOCKOV),ISYMA,WORK(KT1AMP),ISYMA,
     &             XLIAJB,ISYM0,IOPT)
      TIMFCK = TIMFCK  + SECOND() - DTIME

* calculate the XBAR and YBAR intermediates:
      IF (.NOT.CCS) THEN
        DTIME  = SECOND()
        CALL CC_XBAR(WORK(KXBAR),XLIAJB,ISYOVOV,
     &               WORK(KT2AMP),ISYMA, WORK(KEND1),LWRK1 )

        CALL CC_YBAR(WORK(KYBAR),XLIAJB,ISYOVOV,
     &               WORK(KT2AMP),ISYMA, WORK(KEND1),LWRK1 )
        TIME = TIME  + SECOND() - DTIME
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG CCBINT3>  IDLST = ',IDLST
        XNORM = DDOT(N2BST(ISYMA),WORK(KFOCK),1,WORK(KFOCK),1)
        WRITE (LUPRI,*) 'DEBUG CCBINT3>  NORM^2 of FOCK = ',XNORM
        XNORM = DDOT(NMATIJ(ISYMA),WORK(KFOCKOO),1,WORK(KFOCKOO),1)
        WRITE (LUPRI,*) 'DEBUG CCBINT3>  NORM^2 of FOCKOO = ',XNORM
        XNORM = DDOT(NT1AM(ISYMA),WORK(KFOCKOV),1,WORK(KFOCKOV),1)
        WRITE (LUPRI,*) 'DEBUG CCBINT3>  NORM^2 of FOCKOV = ',XNORM
        XNORM = DDOT(NMATAB(ISYMA),WORK(KFOCKVV),1,WORK(KFOCKVV),1)
        WRITE (LUPRI,*) 'DEBUG CCBINT3>  NORM^2 of FOCKVV = ',XNORM
      END IF

      CALL QEXIT('CCBINT3')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCBINT3
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBOPEN*/
*=====================================================================*
      SUBROUTINE CCBOPEN(LUBF,LUCBAR,LUDBAR,LUC,LUD,LUF,LUFK,LUR,
     &                BFFIL,CBAFIL,DBAFIL,CTFIL,DTFIL,FFIL,FKFIL,RFIL,
     &                   LENBF, LENF, LENFK, LENR, 
     &                   NINT1, NINT2, WORK, LWORK      )
*---------------------------------------------------------------------*
*    Purpose: open files for intermediates in B matrix transformation
*      
*     Written by Christof Haettig, Januar/Februar 1997.
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccorb.h"
#include "ccsdsym.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LUBF, LUCBAR, LUDBAR, LUC, LUD, LUF, LUFK, LUR
      INTEGER LWORK, LENBF, LENF, LENFK, LENR, NINT1, NINT2
      CHARACTER*(*) BFFIL,CBAFIL,DBAFIL,CTFIL,DTFIL,FFIL,FKFIL,RFIL
      INTEGER LEN, IINT1, IINT2, ISYM
      
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION WORK(LWORK)
#endif

      CALL QENTER('CCBOPEN')

*---------------------------------------------------------------------*
* open files for local intermediates:
*---------------------------------------------------------------------*
      LUBF   = -1
      LUC    = -1
      LUD    = -1
      LUCBAR = -1
      LUDBAR = -1
      LUR    = -1
      LUF    = -1
      LUFK   = -1
      IF (.NOT. (CCS.OR.CC2)) THEN
         CALL WOPEN2(LUBF,   BFFIL,  64, 0)
         CALL WOPEN2(LUC,    CTFIL,  64, 0)
         CALL WOPEN2(LUD,    DTFIL,  64, 0)
         CALL WOPEN2(LUCBAR, CBAFIL, 64, 0)
         CALL WOPEN2(LUDBAR, DBAFIL, 64, 0)
         CALL WOPEN2(LUR,    RFIL,   64, 0)
      END IF

      IF (.NOT.CCS) THEN
         CALL WOPEN2(LUF,  FFIL, 64, 0)
      END IF

      CALL WOPEN2(LUFK, FKFIL,64, 0)

*---------------------------------------------------------------------*
* calculate a common vector length for BF intermediates and
* initialize them with zeros:
*---------------------------------------------------------------------*
      IF (.NOT. (CCS.OR.CC2)) THEN
         LENBF = 0
         DO ISYM = 1, NSYM
            LENBF = MAX(LENBF,NT2AOIJ(ISYM))
         END DO

         IF (LWORK .LT. LENBF) THEN
           CALL QUIT('OUT OF MEMORY IN CCBOPEN.')
         END IF
      
         CALL DZERO(WORK,LENBF)

         DO IINT1 = 1, NINT1
           LEN = LENBF
           CALL CC_WVEC(LUBF,BFFIL,LENBF,LEN,IINT1,WORK)
         END DO

      END IF

*---------------------------------------------------------------------*
* calculate a common vector length for R intermediates and
* initialize them with zeros:
*---------------------------------------------------------------------*
      IF (.NOT. (CCS.OR.CC2)) THEN
         LENR = 0
         DO ISYM = 1, NSYM
            LENR = MAX(LENR,NEMAT1(ISYM))
         END DO

         IF (LWORK .LT. LENR) THEN
           CALL QUIT('OUT OF MEMORY IN CCBOPEN.')
         END IF
      
         CALL DZERO(WORK,LENR)

         DO IINT1 = 1, NINT1
           LEN = LENR
           CALL CC_WVEC(LUR,RFIL,LENR,LEN,IINT1,WORK)
         END DO

      END IF

*---------------------------------------------------------------------*
* calculate a common vector length for the F intermediates and
* initialize them with zeros:
*---------------------------------------------------------------------*
      IF (CC2) THEN
        LENF = 0
        DO ISYM = 1, NSYM
          LENF = MAX(LENF,NT2AM(ISYM))
        END DO

        IF (LWORK .LT. LENF) THEN
          CALL QUIT('OUT OF MEMORY IN CCBOPEN.')
        END IF
   
        CALL DZERO(WORK,LENF)

        DO IINT2 = 1, NINT2
          LEN = LENF
          CALL CC_WVEC(LUF,FFIL,LENF,LEN,IINT2,WORK)
        END DO
      END IF

      
*---------------------------------------------------------------------*
* calculate a common vector length for the response Fock matrices and
* initialize them with zeros:
*---------------------------------------------------------------------*
      LENFK = 0
      DO ISYM = 1, NSYM
        LENFK = MAX(LENFK,N2BST(ISYM))
      END DO

      IF (LWORK .LT. LENFK) THEN
        CALL QUIT('OUT OF MEMORY IN CCBOPEN.')
      END IF
   
      CALL DZERO(WORK,LENFK)

      DO IINT1 = 1, NINT1
        LEN = LENFK
        CALL CC_WVEC(LUFK,FKFIL,LENFK,LEN,IINT1,WORK)
      END DO

      CALL QEXIT('CCBOPEN')

      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CCBOPEN                          *
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCBSAVE*/
*=====================================================================*
      SUBROUTINE CCBSAVE(IBATCH, I1HGH, I2HGH, INTMED1, INTMED2,
     &                   KRHO2,  LUBF, BFFIL, LENBF,
     &                   KOMEG,  LUF,  FFIL,  LENF,
     &                   KFOCK,  LUFK, FKFIL, LENFK,
     &                   KRIM,   LUR,  RFIL,  LENR,
     &                   NINT1,  NINT2,WORK,  LWORK )
*---------------------------------------------------------------------*
*    Purpose: save intermediates for B matrix transformation on file
*      
*     Written by Christof Haettig, Januar/Februar 1997.
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "cciccset.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LUBF, LUF, LUFK, LUR
      INTEGER NINT1, NINT2, LENBF, LENF, LENFK, LENR, LWORK, IBATCH
      CHARACTER*(*) BFFIL, FFIL, FKFIL, RFIL
      INTEGER I1HGH(0:IBATCH), I2HGH(0:IBATCH)
      INTEGER INTMED1(2,NINT1), INTMED2(4,NINT2)
      INTEGER KRHO2(NINT1), KFOCK(NINT1), KOMEG(NINT2), KRIM(NINT2)

      CHARACTER*(3) LIST, LISTA, LISTB
      INTEGER IDLST, IDLSTA, IDLSTB, ISYM, ISYMA, ISYMB, ISYMAB, LEN
      INTEGER IINT1, IINT2
      
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XNORM
      REAL DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XNORM
      DOUBLE PRECISION DDOT
#endif

* external function:
      INTEGER ILSTSYM

      CALL QENTER('CCBSAVE')

*---------------------------------------------------------------------*
* Fock, BF and R intermediates:
*---------------------------------------------------------------------*
      DO IINT1 = I1HGH(IBATCH-1)+1, I1HGH(IBATCH)
        LIST   = VTABLE(INTMED1(2,IINT1))
        IDLST  = INTMED1(1,IINT1)
        ISYM   = ILSTSYM(LIST,IDLST)

* BF intermediate:
        IF (.NOT. (CCS .OR. CC2)) THEN
          LEN = NT2AOIJ(ISYM)
          CALL CC_WVEC(LUBF,BFFIL, LENBF,LEN,IINT1,WORK(KRHO2(IINT1)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KRHO2(IINT1)),1,WORK(KRHO2(IINT1)),1)
            WRITE (LUPRI,*) 'Norm of saved BF intermediate nb. ',
     &                IINT1, ' is ', XNORM
          END IF
        END IF

* R intermediate:
        IF (.NOT. (CCS .OR. CC2)) THEN
          LEN = NEMAT1(ISYM)
          CALL CC_WVEC(LUR,RFIL,LENR,LEN,IINT1,WORK(KRIM(IINT1)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KRIM(IINT1)),1,WORK(KRIM(IINT1)),1)
            WRITE (LUPRI,*) 'Norm of saved R intermediate nb. ',
     &                IINT1, ' is ', XNORM
          END IF
        END IF

* Fock intermediate:
        IF (.NOT.(CCSD.OR.CCSDT)) THEN
          LEN = N2BST(ISYM)
          CALL  CC_WVEC (LUFK,FKFIL,LENFK,LEN,IINT1,WORK(KFOCK(IINT1)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KFOCK(IINT1)),1,WORK(KFOCK(IINT1)),1)
            WRITE (LUPRI,*) 'Norm of saved FOCK intermediate nb. ',
     &                IINT1, ' is ', XNORM
          END IF
        END IF


      END DO

*---------------------------------------------------------------------*
* F term:
*---------------------------------------------------------------------*
      IF (CC2) THEN
        DO IINT2 = I2HGH(IBATCH-1)+1, I2HGH(IBATCH)
          LISTA  = VTABLE(INTMED2(2,IINT2))
          LISTB  = VTABLE(INTMED2(4,IINT2))
          IDLSTA = INTMED2(1,IINT2)
          IDLSTB = INTMED2(3,IINT2)
          ISYMA  = ILSTSYM(LISTA,IDLSTA)
          ISYMB  = ILSTSYM(LISTB,IDLSTB)
          ISYMAB = MULD2H(ISYMA,ISYMB)
          LEN    = NT2AM(ISYMAB)
          CALL CC_WVEC(LUF,FFIL,LENF,LEN,IINT2,WORK(KOMEG(IINT2)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KOMEG(IINT2)),1,WORK(KOMEG(IINT2)),1)
            WRITE (LUPRI,*) 'Norm of saved F intermediate nb. ',
     &                IINT2, ' is ', XNORM
          END IF
        END DO
      END IF

      CALL QEXIT('CCBSAVE')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCBSAVE
*=====================================================================*

*---------------------------------------------------------------------*
c/* Deck CCB_22CD */
*=====================================================================*
      SUBROUTINE CCB_22CD(THETA2,ISYRES,CDBAR,ISYMCD, 
     &                    T1AMPA,ISYMTA,T1AMPB,ISYMTB,TERM, WORK,LWORK)
*---------------------------------------------------------------------*
*
*  Purpose: to calculate the contribution to the B matrix which
*           are analog to the 22C/D contribution to the right transf.
*
*  assumes:   result vector      THETA2  packed
*             intermediate       CDBAR   squared
*
*  TERM = 'C'  : calculate 22C contribution 
*
*  TERM = 'D'  : calculate 22D contribution 
*               
* 
*  symmetries & variables:   
*
*            ISYRES : result vector THETA2
*            ISYMCD : CDBAR intermediate
*            ISYMTA : response vector A
*            ISYMTB : response vector B
*
*  Christof Haettig, January 1997, based on CCG_22CD
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
      REAL ZERO, HALF, ONE, TWO
#else
      DOUBLE PRECISION ZERO, HALF, ONE, TWO
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, HALF = 0.5d0)

      CHARACTER TERM*(1)

      INTEGER ISYRES, ISYMCD, ISYMTB, ISYMTA
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL THETA2(*)     ! dimension (NT2AM(ISYRES))
      REAL CDBAR(*)      ! dimension (NT2SQ(ISYMCD))
      REAL T1AMPA(*)     ! dimension (NT1AM(ISYMTA))
      REAL T1AMPB(*)     ! dimension (NT1AM(ISYMTB))
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION THETA2(*)     ! dimension (NT2AM(ISYRES))
      DOUBLE PRECISION CDBAR(*)      ! dimension (NT2SQ(ISYMCD))
      DOUBLE PRECISION T1AMPA(*)     ! dimension (NT1AM(ISYMTA))
      DOUBLE PRECISION T1AMPB(*)     ! dimension (NT1AM(ISYMTB))
      DOUBLE PRECISION WORK(LWORK)
#endif

      INTEGER ISYMB, ISYMAIJ, ISYMCKJ, ISYTATB, ISYMAI, ISYMBJ, ISYMJ
      INTEGER KJINT, KSCRT, KEND2, LEND2, KDUM
      INTEGER IOPT, IPCK, NAI, NJ, NBJ, NAIBJ

      INTEGER INDEX
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      CALL QENTER('CCB_22CD')

* check symmetries:
      ISYTATB = MULD2H(ISYMTA,ISYMTB)

      IF (ISYRES .NE. MULD2H(ISYTATB,ISYMCD)) THEN
        CALL QUIT('Symmetry mismatch in CCB_22CD.')
      END IF

* check TERM option:
      IF (TERM .NE. 'C' .AND.  TERM .NE. 'D') THEN
        CALL QUIT('CCB_22CD CALLed with illegal TERM option.')
      END IF

      DO ISYMB = 1, NSYM
        ISYMAIJ  = MULD2H(ISYMB,ISYRES)   ! b batch of result vector
        ISYMCKJ  = MULD2H(ISYMB,ISYMCD)  ! batch of CDBAR^{b}_{ld|i}

        IF (ISYMAIJ .NE. MULD2H(ISYMCKJ,ISYTATB)) THEN
          CALL QUIT('Symmetry mismatch in CCB_22CD.')
        END IF

        KJINT = 1
        KSCRT = KJINT + NT2BCD(ISYMAIJ)
        KEND2 = KSCRT + NT2BCD(ISYMAIJ)
        LEND2 = LWORK - KEND2

        IF (LEND2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCB_22CD.')
        END IF
      
        DO B = 1, NVIR(ISYMB)

*---------------------------------------------------------------------*
*         calculate double transformed CDBAR intermediate 
*          CDBARtt(ai;j) = CDBAR_{a^A i^B, b j} + CDBAR_{a^B i^A, b j}
*---------------------------------------------------------------------*

          IOPT = 1  ! coulomb type result (no exchange type)
          IPCK = 2  ! CDBAR intermediate is stored as a squared vector
          KDUM = -99 999 999 ! dummy address
          CALL CCG_OOVV(WORK(KJINT),WORK(KDUM),ISYMAIJ,CDBAR,ISYMCD,
     &                  T1AMPA, ISYMTA, T1AMPB, ISYMTB, 
     &                  WORK(KEND2), LEND2, B, ISYMB, IOPT, IPCK)

*---------------------------------------------------------------------*
*         for 22D contribution scale with +1/2
*---------------------------------------------------------------------*
          IF (TERM .EQ. 'D') THEN
            CALL DSCAL(NT2BCD(ISYMAIJ),HALF,WORK(KJINT),1)
          END IF

*---------------------------------------------------------------------*
*         for the 22C contribution apply +1/2 * (1 + 2 * Pij)
*---------------------------------------------------------------------*
          IF (TERM .EQ. 'C') THEN
           CALL DCOPY(NT2BCD(ISYMAIJ), WORK(KJINT),1, WORK(KSCRT),1)

           CALL CCLT_P21I(WORK(KSCRT), ISYMAIJ, WORK(KEND2), LEND2,
     &                    IT2BCD, NT2BCD, IT1AM, NT1AM, NVIR)

           CALL DAXPY(NT2BCD(ISYMAIJ),TWO,WORK(KSCRT),1,WORK(KJINT),1)

           CALL DSCAL(NT2BCD(ISYMAIJ),HALF,WORK(KJINT),1)
          END IF
   
*---------------------------------------------------------------------*
*         storage in result vector:
*---------------------------------------------------------------------*
          DO ISYMJ = 1, NSYM
            ISYMAI = MULD2H(ISYMJ,ISYMAIJ)
            ISYMBJ = MULD2H(ISYMJ,ISYMB)
            
            IF (MULD2H(ISYMAI,ISYMBJ) .NE. ISYRES) THEN
              CALL QUIT('Symmetry mismatch in CCB_22CD.')
            END IF

            DO J = 1, NRHF(ISYMJ)
         
              NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
              NJ  = KJINT-1+ IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1)

              IF (ISYMAI .EQ. ISYMBJ) THEN
                 
                DO NAI = 1, NT1AM(ISYMAI)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                
                  THETA2(NAIBJ) = THETA2(NAIBJ) + WORK(NJ + NAI)

                  IF (NAI .EQ. NBJ) THEN
                    THETA2(NAIBJ) = THETA2(NAIBJ) + WORK(NJ + NAI)
                  END IF
                END DO

              ELSE IF (ISYMAI .LT. ISYMBJ) THEN

               NAIBJ = IT2AM(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)+1
      
               CALL DAXPY (NT1AM(ISYMAI), ONE, WORK(NJ+1), 1,
     &                                    THETA2(NAIBJ), 1)

              ELSE IF (ISYMAI .GT. ISYMBJ) THEN

               NAIBJ = IT2AM(ISYMBJ,ISYMAI) + NBJ
      
               CALL DAXPY (NT1AM(ISYMAI), ONE, WORK(NJ+1), 1,
     &                                    THETA2(NAIBJ), NT1AM(ISYMBJ))

              END IF

            END DO ! J
          END DO ! ISYMJ

        END DO ! B
      END DO ! ISYMB

      CALL QEXIT('CCB_22CD')

      RETURN 
      END
*---------------------------------------------------------------------*
*                   END OF ROUTINE CCB_2CD                            *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
c/* Deck CCB_CDBAR */
*=====================================================================*
      SUBROUTINE CCB_CDBAR(TYPE, XIAJB, ISYINT, T2AMP, ISYTAM,
     &                     CDBAR, ISYCDBAR, WORK, LWORK,
     &                     FILE, LUNIT, IOFFSET, IOPT)
*---------------------------------------------------------------------*
*    Purpose: calculate CBAR/DBAR intermediates
*
*    TYPE='C' : calculate CBAR intermediate
*    TYPE='D' : calculate DBAR intermediate
*
*    IOPT=1,3 : store intermediate in CDBAR
*    IOPT=2,3 : write CDBAR intermediate to FILE (LUNIT), starting
*               at position IOFFSET+1
*
*
*    N.B. for TYPE='D' the amplitudes T2AMP will be overwritten by
*         2*t(ai|bj) - t(aj|bi)
*      
*    Written by Christof Haettig, Januar/Februar 1997.
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*(*) TYPE, FILE
      INTEGER LWORK, ISYTAM, ISYINT, ISYCDBAR, LUNIT, IOFFSET, IOPT
      INTEGER ISYMA, ISYMI, ISYMCK, ISYTINT, ISYCINT, ISYMAI
      INTEGER KTINT, KCINT, KEND1, LWRK1, KOFF1, KOFF2, LEN, NAI, IERR

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL T2AMP(*)      ! dimension (NT2AM(ISYTAM))
      REAL XIAJB(*)      ! dimension (NT2SQ(ISYINT))
      REAL CDBAR(*)      ! dimension (NT2SQ(ISYCDBAR))
      REAL ONE, XNORM, DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION T2AMP(*)      ! dimension (NT2AM(ISYTAM))
      DOUBLE PRECISION XIAJB(*)      ! dimension (NT2SQ(ISYINT))
      DOUBLE PRECISION CDBAR(*)      ! dimension (NT2SQ(ISYCDBAR))
      DOUBLE PRECISION ONE, XNORM, DDOT
#endif
      PARAMETER (ONE = 1.0d0)

      INTEGER IOPTZWVI, ISYTIN, ISYCIN, IOPTTCME

      CALL QENTER('CCB_CDBAR')

*---------------------------------------------------------------------*
* check symmetries:
*---------------------------------------------------------------------*
      IF (MULD2H(ISYTAM,ISYINT) .NE. ISYCDBAR) THEN
        WRITE (LUPRI,*) 'ERROR> SYMMETRY MISMATCH IN CCB_CDBAR.'
        CALL QUIT('SYMMETRY MISMATCH IN CCB_CDBAR.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Entered CDBAR: TYPE, LWORK:',TYPE,LWORK
        XNORM =  0.0d0
      END IF

*---------------------------------------------------------------------*
* prepare (ia|jb) integrals, and amplitudes for contraction:
*---------------------------------------------------------------------*
      IF (TYPE(1:1).EQ.'C') THEN

*       for CBAR intermediate transpose (ia|jb) to (ja|ib):

        CALL CCSD_T2TP(XIAJB,WORK,LWORK,ISYINT)

        IOPTZWVI = 2

      ELSE IF (TYPE(1:1).EQ.'D') THEN

*       for DBAR intermediate calculate L(ia|jb) and 
*       2T(ia|jb)-T(ib|ja) in place:
      
        CALL CCRHS_T2TR(XIAJB,WORK,LWORK,ISYINT)
        IOPTTCME = 1
        CALL CCSD_TCMEPK(T2AMP,ONE,ISYTAM,IOPTTCME)

        IOPTZWVI = 1

      ELSE
        CALL QUIT('ILLEGAL OPTION IN CCB_CDBAR.')
      END IF

*---------------------------------------------------------------------*
* start loop over virtual orbital index a:
*---------------------------------------------------------------------*
      DO ISYMA = 1, NSYM
        ISYTIN = MULD2H(ISYTAM,ISYMA)
        ISYCIN = MULD2H(ISYTIN,ISYINT)

        KTINT = 1
        KCINT = KTINT + NT2BCD(ISYTIN)
        KEND1 = KCINT + NT2BCD(ISYCIN)

        LWRK1 = LWORK - KEND1
        IF (LWRK1 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CCB_CDBAR.')
        END IF

      DO A = 1, NVIR(ISYMA)
        
*---------------------------------------------------------------------*
* get t^{a}(bj,i) submatrix of the t amplitudes
*---------------------------------------------------------------------*
        CALL CCG_TI(WORK(KTINT),ISYTIN,T2AMP,ISYTAM,A,ISYMA)

*---------------------------------------------------------------------*
* calculate CBAR^{a}(ck,i) intermediate
*---------------------------------------------------------------------*
        Call CC_ZWVI(WORK(KCINT), XIAJB, ISYINT, WORK(KTINT),
     &               ISYTIN, WORK(KEND1), LWRK1, IOPTZWVI)

*---------------------------------------------------------------------*
* resort to CDBAR(ck,ai) and write to output variable/file:
* (for CBAR intermediate scale with -1)
*---------------------------------------------------------------------*
        DO ISYMI = 1, NSYM
          ISYMCK = MULD2H(ISYCIN,ISYMI)
          ISYMAI = MULD2H(ISYMA,ISYMI)
          LEN    = NT1AM(ISYMCK)

        IF (LEN.GT.0) THEN

        DO I = 1, NRHF(ISYMI)
          NAI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
          KOFF1 = KCINT + IT2BCD(ISYMCK,ISYMI) + NT1AM(ISYMCK)*(I-1) 
          KOFF2 = IT2SQ(ISYMCK,ISYMAI) + NT1AM(ISYMCK)*(NAI-1) + 1

          IF (TYPE(1:1).EQ.'C') CALL DSCAL(LEN,-ONE,WORK(KOFF1),1)
      
          IF (IOPT.EQ.1 .OR. IOPT.EQ.3) THEN
            CALL DCOPY(LEN,WORK(KOFF1),1,CDBAR(KOFF2),1)
          END IF
          IF (LOCDBG) THEN
            XNORM = XNORM + DDOT(LEN,WORK(KOFF1),1,WORK(KOFF1),1)
          END IF

          IF (IOPT.EQ.2 .OR. IOPT.EQ.3) THEN
            CALL PUTWA2(LUNIT,FILE,WORK(KOFF1),IOFFSET+KOFF2,LEN)
          END IF

        END DO

        END IF

        END DO


      END DO
      END DO

*---------------------------------------------------------------------*
* reconstruct (ia|jb) integrals:
*---------------------------------------------------------------------*
      IF (TYPE(1:1).EQ.'C') THEN
        CALL CCSD_T2TP(XIAJB,WORK,LWORK,ISYINT)
      ELSE IF (TYPE(1:1).EQ.'D') THEN
        CALL CCRHS_T2BT(XIAJB,WORK,LWORK,ISYINT)
      ELSE
        CALL QUIT('ILLEGAL OPTION IN CCB_CDBAR.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Norm of ',TYPE,'BAR intermediate is ',XNORM
      END IF

      CALL QEXIT('CCB_CDBAR')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCB_CDBAR
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_FDB(NC1VEC,NC2VEC,NCR12VEC,TXAM,TYAM,RESULT,
     &                  WORK,LWORK,APROXR12)
*---------------------------------------------------------------------
* Test routine for calculating the CC B matrix by finite difference 
* on the right hand Jacobian transformation.
* Ch. Haettig, februar 1997
*
* adapted for CC-R12
* Christian Neiss, november 2005
*---------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "iratdef.h"
#include "ccorb.h"
#include "aovec.h"
#include "ccsdinp.h"
#include "cclr.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "leinf.h"
#include "r12int.h"
#include "ccr12int.h"
C
      DIMENSION WORK(LWORK),ITADR(2),RESULT(*)
      PARAMETER (XHALF = 0.5D00,XMTWO = -2.0D00, DELTA = 1.0D-07)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, TWO = 2.0d0)
      CHARACTER MODEL*10, APROXR12*3
      LOGICAL L1TST,L2TST, LETST
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CC_FDB')
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
      IF (CC3) MODEL = 'CC3       '
C
      IF (CCR12) CALL CCSD_MODEL(MODEL,LENMOD,24,MODEL,10,APROXR12)
C
      IF (IPRINT.GT.5) THEN
         CALL AROUND( 'IN CC_FDB  : MAKING FINITE DIFF. CC B Matrix')
      ENDIF
C
C----------------------
C     Set Test options.
C----------------------
C
      L1TST = .FALSE.
      L2TST = .FALSE.
      LETST = .TRUE.
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      ISYMOP     = 1
C
      IF (CCR12) THEN
        NTAMR12  = NTR12AM(ISYMTR)
      ELSE
        NTAMR12  = 0
      END IF
C
      NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR) + NTAMR12
      NTAMP2     = NTAMP*(NC1VEC + NC2VEC + NCR12VEC)
      KF         = 1
      KRHO1      = KF       + NTAMP2
      KRHO2      = KRHO1    + NT1AMX
      KRHO12     = KRHO2    + MAX(NT2AMX,NT2AM(ISYMTR))
C     KC1AM      = KRHO12   + NTAMR12
C     KC2AM      = KC1AM    + NT1AM(ISYMTR)
C     KC12AM     = KC2AM    + NT2AM(ISYMTR)
C     KEND1      = KC2AM   
C    *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2SQ(ISYMTR),
C    *                 NT2R12(ISYMTR)) + NTAMR12
C     KEND1      = KC12AM   + NTAMR12
C     LWRK1      = LWORK    - KEND1
C
C     KRHO1D     = KEND1
C     KRHO2D     = KRHO1D   + NT1AMX
C     KRHO12D    = KRHO2D   + NT2AMX
C     KEND2      = KRHO2D   
C    *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2AO(ISYMTR),
C    *                 2*NT2ORT(ISYMTR)) + NTAMR12
C     LWRK2      = LWORK    - KEND2
C
      KC1AM      = KRHO12   + NTAMR12
      KEND1      = KC1AM    + MAX(NT1AM(ISYMTR),NTAMP)
      LWRK1      = LWORK    - KEND1
C
      KRHO1D     = KEND1
      KRHO2D     = KRHO1D   + NT1AMX
      KRHO12D    = KRHO2D   + NT2AMX
      KC2AM      = KRHO2D
     *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2AO(ISYMTR),
     *                 2*NT2ORT(ISYMTR)) + NTAMR12 
      KC12AM     = KC2AM    + NT2AM(ISYMTR)
      KEND2      = KC2AM
     *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2SQ(ISYMTR),
     *                 NT2R12(ISYMTR)) + NTAMR12 
      LWRK2      = LWORK    - KEND2     
C
      IF (IPRINT .GT. 100 ) THEN
         WRITE(LUPRI,*) ' IN CC_FDB: KF      =  ',KF     
         WRITE(LUPRI,*) ' IN CC_FDB: KRHO1   =  ',KRHO1
         WRITE(LUPRI,*) ' IN CC_FDB: KRHO2   =  ',KRHO2
         WRITE(LUPRI,*) ' IN CC_FDB: KC1AM   =  ',KC1AM
         WRITE(LUPRI,*) ' IN CC_FDB: KC2AM   =  ',KC2AM
         WRITE(LUPRI,*) ' IN CC_FDB: KRHO1D  =  ',KRHO1D
         WRITE(LUPRI,*) ' IN CC_FDB: KRHO2D  =  ',KRHO2D
         WRITE(LUPRI,*) ' IN CC_FDB: KEND2   =  ',KEND2
         WRITE(LUPRI,*) ' IN CC_FDB: LWRK2   =  ',LWRK2
      ENDIF
      IF (LWRK2.LT.0 ) THEN
         WRITE(LUPRI,*) 'Too little work space in CC_FDB '
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',KEND2
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDB ')
      ENDIF
      KF2   = KF      + NC1VEC*NTAMP
      KFR   = KF      + (NC1VEC+NC2VEC)*NTAMP
C
C---------------------
C     Initializations.
C---------------------
C
      CALL DZERO(WORK(KC1AM),NT1AMX)
      CALL DZERO(WORK(KC2AM),NT2AMX)
      CALL DZERO(WORK(KC12AM),NTAMR12)
      CALL DZERO(WORK(KF),NTAMP2)
      IF (ABS(DELTA) .GT. 1.0D-15 ) THEN 
         DELTAI = 1.0D00/DELTA
      ELSE
         CALL QUIT('DELTA too small in CC_FDB')
C        DELTAI = 1
      ENDIF
      X11 = 0.0D00
      X12 = 0.0D00
      X21 = 0.0D00
      X22 = 0.0D00
      XNJ = 0.0D00
      XR1 = 0.0D00
      X1R = 0.0D00
      XR2 = 0.0D00
      X2R = 0.0D00
      XRR = 0.0D00
C
C------------------------------------------------
C     Read the CC reference amplitudes From disk.
C------------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KC1AM),WORK(KC2AM))
C
      IF (CCR12) THEN
        IOPT = 32
        CALL CC_RDRSP('R0 ',0,1,IOPT,MODEL,DUMMY,WORK(KC12AM))
      END IF
C
C----------------------------------------------
C     Save the CC reference amplitudes on disk.
C----------------------------------------------
C
      LUTAM = -1
      CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUTAM)
      WRITE(LUTAM) (WORK(KC1AM + I -1 ), I = 1, NT1AMX)
      WRITE(LUTAM) (WORK(KC2AM + I -1 ), I = 1, NT2AMX)
      WRITE(LUTAM) (WORK(KC12AM+ I -1 ), I = 1, NTAMR12) 
      CALL GPCLOSE(LUTAM,'KEEP')
C
      IF (IPRINT.GT.125) THEN
         RHO1N = DDOT(NT1AMX,WORK(KC1AM),1,WORK(KC1AM),1)
         RHO2N = DDOT(NT2AMX,WORK(KC2AM),1,WORK(KC2AM),1)
         RHO12N = DDOT(NTAMR12,WORK(KC12AM),1,WORK(KC12AM),1)
         WRITE(LUPRI,*) 'Norm of T1AM: ',RHO1N
         WRITE(LUPRI,*) 'Norm of T2AM: ',RHO2N
         IF (CCR12) WRITE(LUPRI,*) 'Norm of R12 amplitudes: ',RHO12N
         CALL CC_PRP(WORK(KC1AM),WORK(KC2AM),1,1,1)
         IF (CCR12) CALL CC_PRPR12(WORK(KC12AM),1,1,.TRUE.)
      ENDIF
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *            WORK(KC2AM),WORK(KEND2),LWRK2,APROXR12)
C
C------------------------------------------
C     Calculate reference A*T vector.
C------------------------------------------
C
      CALL DCOPY(NTAMP,TXAM,1,WORK(KRHO1D),1)
      ISIDE = +1
      CALL CC_ATRR(0.0D0,1,ISIDE,WORK(KRHO1D),LWRK1,.FALSE.,DUMMY,
     &             APROXR12)
C
C-------------------------
C     Zero out components.
C-------------------------
C
      IF (LCOR .OR. LSEC) THEN
C
         CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
C
      ENDIF
C
      IF (IPRINT.GT.2) THEN
         RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
         RHO12N = DDOT(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12D),1)
         WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'ref'
         WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'ref'
         IF (CCR12) WRITE(LUPRI,*) 'Norm of RHO-R12: ',RHO12N,'ref'
      ENDIF
      IF (IPRINT.GT.125) THEN
         CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
         IF (CCR12) CALL CC_PRPR12(WORK(KRHO12D),1,1,.TRUE.)
      ENDIF

      CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1),1)
      CALL DCOPY(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2),1)
      CALL DCOPY(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12),1)
C
C==================================================
C calculate intermediates for response vector TXAM:
C==================================================
C
      IF (.FALSE.) THEN
         IOPT = 3
         CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2)
C
         IF (CCR12) THEN
           IOPT = 32
           CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,DUMMY,DUMMY,
     &                   WORK(KC12AM),WORK(KEND2),LWRK2)
         END IF
C
         WRITE (LUPRI,*) 'NTAMP:',NTAMP
         WRITE (LUPRI,*) 'NORM TXAM:',DDOT(NTAMP,TXAM,1,TXAM,1)
C
         RSPIM = .TRUE.
         CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND2),LWRK2,APROXR12) 
      END IF
C-------------------------------------
C read E intermediates:
C-------------------------------------
      IF (LETST) THEN
         KEI1  = KEND2
         KEI2  = KEI1 + NMATAB(1)
         KEND3 = KEI2 + NMATIJ(1)
         LWRK3 = LWORK - KEND3

         IF (LWRK3.LT.0 ) THEN
           CALL QUIT('Insufficient memory in CC_FDB.')
         END IF
C
         LUE1 = -1
         CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE1)
         READ (LUE1)(WORK(KEI1+ J-1),J = 1,NMATAB(ISYMOP))
         CALL GPCLOSE(LUE1,'KEEP' )
C
         LUE2 = -1
         CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUE2)
         READ (LUE2) (WORK(KEI2+ J-1),J = 1,NMATIJ(ISYMOP))
         CALL GPCLOSE(LUE2,'KEEP' )
C
         CALL AROUND( 'E^X-intermediates read from disk ')
         CALL CC_PREI(WORK(KEI1),WORK(KEI2),ISYMOP,1)
      END IF

C
C=============================================
C     Calculate B-matrix by finite difference.
C=============================================
C
      DO 100 I = 1, NC1VEC
         WRITE (LUPRI,*) 'singles index:',I
C
C----------------------------------------
C        Add finite displadement to t and 
C        calculate new intermediates.
C----------------------------------------
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         REWIND(LUTAM)
         READ(LUTAM) (WORK(KC1AM + J -1 ) , J = 1, NT1AMX)
         READ(LUTAM) (WORK(KC2AM + J -1 ) , J = 1, NT2AMX)
         READ(LUTAM) (WORK(KC12AM+ J -1 ) , J = 1, NTAMR12)
         CALL GPCLOSE(LUTAM,'KEEP')
C
         TI   = SECOND()
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1 ) + DELTA
         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
         ENDIF
C
         IOPT = 3
         CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2)
C
         IF (CCR12) THEN
           IOPT = 32
           CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,DUMMY,DUMMY,
     &                   WORK(KC12AM),WORK(KEND2),LWRK2)
         END IF
C
         RSPIM = .TRUE.
         CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND2),LWRK2,APROXR12) 
C
C---------------------------------------------
C        Get the CC response vector again.
C---------------------------------------------
C
         CALL DCOPY(NTAMP,TXAM,1,WORK(KC1AM),1)
C
C---------------------------------------
C        For Test zero part of T vector.
C---------------------------------------
C     
         IF ( L1TST ) THEN
C           CALL DZERO(WORK(KC2AM),NT2AMX)
C           CALL DZERO(WORK(KC12AM),NTAMR12)
            CALL DZERO(WORK(KC1AM+NT1AMX),NT2AMX)
            CALL DZERO(WORK(KC1AM+NT1AMX+NT2AMX),NTAMR12)
         ENDIF
         IF ( L2TST ) THEN
            CALL DZERO(WORK(KC1AM),NT1AMX)
C           CALL DZERO(WORK(KC12AM),NTAMR12)
            CALL DZERO(WORK(KC1AM+NT1AMX+NT2AMX),NTAMR12)
         ENDIF
C
C------------------
C        Transform.
C------------------
C
         CALL DCOPY(NTAMP,WORK(KC1AM),1,WORK(KRHO1D),1)

         ISIDE = +1
         CALL CC_ATRR(0.0D0,1,ISIDE,WORK(KRHO1D),LWRK1,.FALSE.,DUMMY,
     &             APROXR12)
C
         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
         ENDIF
C
         IF (IPRINT.GT.2) THEN
            RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
            RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
            RHO12N = DDOT(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12D),1)
            WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'ai=',I
            WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'ai=',I
            IF (CCR12) WRITE(LUPRI,*) 'Norm of RHO-R12: ',RHO12N,'ai=',I
         ENDIF
         IF (IPRINT.GT.125) THEN
            CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
            IF (CCR12) CALL CC_PRPR12(WORK(KRHO12D),1,1,.TRUE.)
         ENDIF

         CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
         CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
         CALL DAXPY(NTAMR12,-1.0D00,WORK(KRHO12),1,WORK(KRHO12D),1)
         CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
         CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
         CALL DSCAL(NTAMR12,DELTAI,WORK(KRHO12D),1)
         CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *              WORK(KF+NTAMP*(I-1)),1)
         CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *              WORK(KF+NTAMP*(I-1)+NT1AMX),1)
         CALL DCOPY(NTAMR12,WORK(KRHO12D),1,
     &              WORK(KF+NTAMP*(I-1)+NT1AMX+NT2AMX),1)
         X11 = X11 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         X21 = X21 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
         XR1 = XR1 + DDOT(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12D),1)
C
         TI   = SECOND() - TI
         IF (IPRINT.GT.5 ) THEN
            WRITE(LUPRI,*) '  '
            WRITE(LUPRI,*) 'FDB ROW NR. ',I,' DONE IN ',TI,' SEC.'
         ENDIF
C
 100  CONTINUE
C
C----------------------------------------------------------------
C     Loop over T2 amplitudes. Take care of diagonal t2 elements
C     is in a different convention in the energy code.
C     Factor 1/2 from right , and factor 2 from left.
C----------------------------------------------------------------
C
      IF (.NOT. (CCS .OR. CCSTST)) THEN
      DO 200 NAI = 1, NT1AMX
        DO 300 NBJ = 1, NAI
         I = INDEX(NAI,NBJ)
C
         IF (I.LE.NC2VEC) THEN
           WRITE (LUPRI,*) 'doubles index:',I
C
C--------------------------------------------
C          Add finite displacement to t and
C          calculate new intermediates.
C-------------------------------------------
C
           LUTAM = -1
           CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           REWIND(LUTAM)
           READ(LUTAM) (WORK(KC1AM + J -1 ) , J = 1, NT1AMX)
           READ(LUTAM) (WORK(KC2AM + J -1 ) , J = 1, NT2AMX)
           READ(LUTAM) (WORK(KC12AM+ J -1 ) , J = 1, NTAMR12)
           CALL GPCLOSE(LUTAM,'KEEP')
C
           TI   = SECOND()
           DELT = DELTA
           IF (NAI.EQ.NBJ) DELT = 2*DELTA
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) + DELT
           IF (LCOR .OR. LSEC) THEN
             CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
           ENDIF
C
           IOPT = 3
           CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                   WORK(KC2AM),WORK(KEND2),LWRK2)
C
           IF (CCR12) THEN
             IOPT = 32
             CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,DUMMY,DUMMY,
     &                     WORK(KC12AM),WORK(KEND2),LWRK2)
           END IF
C
           RSPIM = .TRUE.
           CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2,APROXR12) 
C
C-----------------------------------------------
C          Get the CC response vector again.
C-----------------------------------------------
C
           CALL DCOPY(NTAMP,TXAM,1,WORK(KC1AM),1)
C
C-----------------------------------------
C          For Test zero part of T vector.
C-----------------------------------------
C     
           IF ( L1TST ) THEN
C             CALL DZERO(WORK(KC2AM),NT2AMX)
C             CALL DZERO(WORK(KC12AM),NTAMR12)
              CALL DZERO(WORK(KC1AM+NT1AMX),NT2AMX)
              CALL DZERO(WORK(KC1AM+NT1AMX+NT2AMX),NTAMR12)
           ENDIF
           IF ( L2TST ) THEN
              CALL DZERO(WORK(KC1AM),NT1AMX)
C             CALL DZERO(WORK(KC12AM),NTAMR12)
              CALL DZERO(WORK(KC1AM+NT1AMX+NT2AMX),NTAMR12)
           ENDIF
C
           RHO1N = DDOT(NT1AMX,WORK(KC1AM),1,WORK(KC1AM),1)
           RHO2N = DDOT(NT2AMX,WORK(KC2AM),1,WORK(KC2AM),1)
           RHO12N = DDOT(NTAMR12,WORK(KC12AM),1,WORK(KC12AM),1)
           IF ( DEBUG ) THEN
              WRITE(LUPRI,*) 'Norm of L1AM-inp: ',RHO1N
              WRITE(LUPRI,*) 'Norm of L2AM-inp: ',RHO2N
              IF (CCR12) WRITE(LUPRI,*) 'Norm of LR12AM-inp: ',RHO12N
           ENDIF
C
C--------------------
C          Transform.
C--------------------
C
           CALL DCOPY(NTAMP,WORK(KC1AM),1,WORK(KRHO1D),1)

           ISIDE = +1
           CALL CC_ATRR(0.0D0,1,ISIDE,WORK(KRHO1D),LWRK1,.FALSE.,DUMMY,
     &             APROXR12)
C
           IF (LCOR .OR. LSEC) THEN
              CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
           ENDIF
C
           IF (IPRINT.GT.2) THEN
             RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
             RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
             RHO12N = DDOT(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12D),1)
             WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'aibj=',I
             WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'aibj=',I
             IF (CCR12) WRITE(LUPRI,*) 'Norm of RHO-R12: ',RHO12N,
     &                                 'aibj=',I
           ENDIF
           IF (IPRINT.GT.125) THEN
            CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
            IF (CCR12) CALL CC_PRPR12(WORK(KRHO12D),1,1,.TRUE.)
           ENDIF
C
           CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
           CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
           CALL DAXPY(NTAMR12,-1.0D00,WORK(KRHO12),1,WORK(KRHO12D),1)
           CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
           CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
           CALL DSCAL(NTAMR12,DELTAI,WORK(KRHO12D),1)
           CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KF2+NTAMP*(I-1)),1)
           CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *              WORK(KF2+NTAMP*(I-1)+NT1AMX),1)
           CALL DCOPY(NTAMR12,WORK(KRHO12D),1,
     &              WORK(KF2+NTAMP*(I-1)+NT1AMX+NT2AMX),1)
C
           X12 = X12 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
           X22 = X22 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
           XR2 = XR2 + DDOT(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12D),1)
           TI   = SECOND() - TI
           IF (IPRINT.GT.5 ) THEN
              WRITE(LUPRI,*) '  '
              WRITE(LUPRI,*) 'FDB ROW NR. ',I+NT1AMX,
     *                  ' DONE IN ',TI,' SEC.'
           ENDIF
C
         ENDIF
C
 300    CONTINUE
 200  CONTINUE
      END IF
C
C----------------------------------------------------------------
C     Loop over R12 amplitudes.
C----------------------------------------------------------------
C
      IF (CCR12) THEN
      DO NKI = 1, NMATKI(1)
        DO NLJ = 1, NKI
         I = INDEX(NKI,NLJ)
C
         IF (I.LE.NCR12VEC) THEN
           WRITE (LUPRI,*) 'R12 doubles index:',I
C
C--------------------------------------------
C          Add finite displacement to t and
C          calculate new intermediates.
C-------------------------------------------
C
           LUTAM = -1
           CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           REWIND(LUTAM)
           READ(LUTAM) (WORK(KC1AM + J -1 ) , J = 1, NT1AMX)
           READ(LUTAM) (WORK(KC2AM + J -1 ) , J = 1, NT2AMX)
           READ(LUTAM) (WORK(KC12AM+ J -1 ) , J = 1, NTAMR12)
           CALL GPCLOSE(LUTAM,'KEEP')
C
           TI   = SECOND()
           DELT = DELTA
           IF (NKI.EQ.NLJ) DELT = KETSCL*DELTA
           WORK(KC12AM + I -1) = WORK(KC12AM+I -1) + DELT
           IF (LCOR .OR. LSEC) THEN
             CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
           ENDIF
C
           IOPT = 3
           CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                   WORK(KC2AM),WORK(KEND2),LWRK2)
C
           IF (CCR12) THEN
             IOPT = 32
             CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,DUMMY,DUMMY,
     &                     WORK(KC12AM),WORK(KEND2),LWRK2)
           END IF
C
           RSPIM = .TRUE.
           CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2,APROXR12) 
C
C-----------------------------------------------
C          Get the CC response vector again.
C-----------------------------------------------
C
           CALL DCOPY(NTAMP,TXAM,1,WORK(KC1AM),1)
C
C-----------------------------------------
C          For Test zero part of T vector.
C-----------------------------------------
C     
           IF ( L1TST ) THEN
C             CALL DZERO(WORK(KC2AM),NT2AMX)
C             CALL DZERO(WORK(KC12AM),NTAMR12)
              CALL DZERO(WORK(KC1AM+NT1AMX),NT2AMX)
              CALL DZERO(WORK(KC1AM+NT1AMX+NT2AMX),NTAMR12)
           ENDIF
           IF ( L2TST ) THEN
              CALL DZERO(WORK(KC1AM),NT1AMX)
C             CALL DZERO(WORK(KC12AM),NTAMR12)
              CALL DZERO(WORK(KC1AM+NT1AMX+NT2AMX),NTAMR12)
           ENDIF
C
           RHO1N = DDOT(NT1AMX,WORK(KC1AM),1,WORK(KC1AM),1)
           RHO2N = DDOT(NT2AMX,WORK(KC2AM),1,WORK(KC2AM),1)
           RHO12N = DDOT(NTAMR12,WORK(KC12AM),1,WORK(KC12AM),1)
           IF ( DEBUG ) THEN
              WRITE(LUPRI,*) 'Norm of L1AM-inp: ',RHO1N
              WRITE(LUPRI,*) 'Norm of L2AM-inp: ',RHO2N
              IF (CCR12) WRITE(LUPRI,*) 'Norm of LR12AM-inp: ',RHO12N
           ENDIF
C
C--------------------
C          Transform.
C--------------------
C
           CALL DCOPY(NTAMP,WORK(KC1AM),1,WORK(KRHO1D),1)

           ISIDE = +1
           CALL CC_ATRR(0.0D0,1,ISIDE,WORK(KRHO1D),LWRK1,.FALSE.,DUMMY,
     &             APROXR12)
C
           IF (LCOR .OR. LSEC) THEN
              CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
           ENDIF
C
           IF (IPRINT.GT.2) THEN
             RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
             RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
             RHO12N = DDOT(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12D),1)
             WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'kilj=',I
             WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'kilj=',I
             IF (CCR12) WRITE(LUPRI,*) 'Norm of RHO-R12: ',RHO12N,
     &                                 'kilj=',I
           ENDIF
           IF (IPRINT.GT.125) THEN
            CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
            IF (CCR12) CALL CC_PRPR12(WORK(KRHO12D),1,1,.TRUE.)
           ENDIF
C
           CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
           CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
           CALL DAXPY(NTAMR12,-1.0D00,WORK(KRHO12),1,WORK(KRHO12D),1)
           CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
           CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
           CALL DSCAL(NTAMR12,DELTAI,WORK(KRHO12D),1)
           CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KFR+NTAMP*(I-1)),1)
           CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *              WORK(KFR+NTAMP*(I-1)+NT1AMX),1)
           CALL DCOPY(NTAMR12,WORK(KRHO12D),1,
     &              WORK(KFR+NTAMP*(I-1)+NT1AMX+NT2AMX),1)
C
           X1R = X1R + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
           X2R = X2R + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
           XRR = XRR + DDOT(NTAMR12,WORK(KRHO12D),1,WORK(KRHO12D),1)
           TI   = SECOND() - TI
           IF (IPRINT.GT.5 ) THEN
              WRITE(LUPRI,*) '  '
              WRITE(LUPRI,*) 'FDB ROW NR. ',I+NT1AMX+NT2AMX,
     *                  ' DONE IN ',TI,' SEC.'
           ENDIF
C
         ENDIF
C
        END DO
      END DO
      END IF
C
      WRITE(LUPRI,*) '    '
      WRITE(LUPRI,*) '**  FINITE DIFF WITH DELTA ',DELTA, '**'
      WRITE(LUPRI,*) '    '
      IF ( IPRINT .GT. 4 ) THEN
         CALL AROUND('FINITE DIFF. CC B*Tx-Matrix - 11 & 21 (& R12_1)'//
     &               ' PART')
         CALL OUTPUT(WORK(KF),1,NTAMP,1,NC1VEC,NTAMP,NC1VEC,1,LUPRI)
         CALL AROUND('FINITE DIFF. CC B*Tx-Matrix - 12 & 22 (& R12_2)'//
     &               ' PART')
         CALL OUTPUT(WORK(KF+NTAMP*NC1VEC),1,NTAMP,1,NC2VEC,
     *               NTAMP,NC2VEC,1,LUPRI)
         IF (CCR12) THEN
           CALL AROUND('FINITE DIFF. CC B*Tx-Matrix - 1R12 & 2R12 & '//
     &                 'R12R12 PART')
           CALL OUTPUT(WORK(KF+NTAMP*(NC1VEC+NC2VEC)),1,NTAMP,
     &                 1,NCR12VEC,NTAMP,NCR12VEC,1,LUPRI)
         END IF
      ENDIF

      XNJ = X11 + X12 + X21 + X22 + XR1 + XR2 + X1R + X2R + XRR
      WRITE(LUPRI,*) '  '
      WRITE(LUPRI,*) ' NORM OF FIN. DIFF. B*Tx-Matrix.', SQRT(XNJ)
      WRITE(LUPRI,*) '  '
      WRITE(LUPRI,*) ' NORM OF 11 PART OF FD. B*Tx-mat.: ', SQRT(X11)
      IF (.NOT.(CCS.OR.CCSTST)) THEN
         WRITE(LUPRI,*) ' NORM OF 21 PART OF FD. B*Tx-mat.: ', SQRT(X21)
         WRITE(LUPRI,*) ' NORM OF 12 PART OF FD. B*Tx-mat.: ', SQRT(X12)
         WRITE(LUPRI,*) ' NORM OF 22 PART OF FD. B*Tx-mat.: ', SQRT(X22)
      ENDIF
      IF (CCR12) THEN
         WRITE(LUPRI,*) ' NORM OF R12_1 PART OF FD. B*Tx-mat.: ',
     &                  SQRT(XR1)
         WRITE(LUPRI,*) ' NORM OF R12_2 PART OF FD. B*Tx-mat.: ',
     &                  SQRT(XR2)
         WRITE(LUPRI,*) ' NORM OF 1R12 PART OF FD. B*Tx-mat.: ',
     &                  SQRT(X1R)
         WRITE(LUPRI,*) ' NORM OF 2R12 PART OF FD. B*Tx-mat.: ',
     &                  SQRT(X2R)
         WRITE(LUPRI,*) ' NORM OF R12R12 PART OF FD. B*Tx-mat.: ',
     &                  SQRT(XRR)
      END IF
C
C--------------------------------------
C     Calculate Matrix times Ty vector.
C--------------------------------------
C
      CALL DGEMV('N',NTAMP,NTAMP,ONE,WORK(KF),NTAMP,TYAM,1,
     *           ZERO,RESULT,1)

      IF (CCS.OR.CCSTST) THEN
         CALL DZERO(RESULT(NT1AM(1)+1),NT2AM(1))
      END IF
C-----------------------------------------------------------
C     scale diagonal with 1/2 and print norm of the vectors:
C-----------------------------------------------------------
      CALL CCLR_DIASCL(RESULT(NT1AM(1)+1),TWO,1)

      WRITE (LUPRI,*) 'NTAMP:',NTAMP
      WRITE (LUPRI,*) 'NORM^2 OF TX VEC.:',
     &     DDOT(NT1AM(1)+NT2AM(1),TXAM,1,TXAM,1)
      WRITE (LUPRI,*) 'single-excitation part:',
     &     DDOT(NT1AM(1),TXAM,1,TXAM,1)
      WRITE (LUPRI,*) 'NORM^2 OF TY VEC.:',
     &     DDOT(NT1AM(1)+NT2AM(1),TYAM,1,TYAM,1)
      WRITE (LUPRI,*) 'single-excitation part:',
     &     DDOT(NT1AM(1),TYAM,1,TYAM,1)
      WRITE (LUPRI,*) 'NORM^2 OF RESULT VECTOR:',
     &     DDOT(NTAMP,RESULT,1,RESULT,1)
      WRITE (LUPRI,*) 'single-excitation part:',
     &     DDOT(NT1AM(1),RESULT,1,RESULT,1)
C
C-------------------------------------------------
C     Restore the CC reference amplitudes on disk.
C-------------------------------------------------
C
      LUTAM = -1
      CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUTAM)
      READ(LUTAM) (WORK(KC1AM + I -1 ) , I = 1, NT1AMX)
      READ(LUTAM) (WORK(KC2AM + I -1 ) , I = 1, NT2AMX)
      READ(LUTAM) (WORK(KC12AM+ J -1 ) , J = 1, NTAMR12)
      CALL GPCLOSE(LUTAM,'DELETE')
C
      IOPT = 3
      CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *              WORK(KC2AM),WORK(KEND2),LWRK2)
C
      IF (CCR12) THEN
        IOPT = 32
        CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,DUMMY,DUMMY,
     &                WORK(KC12AM),WORK(KEND2),LWRK2)
      END IF
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),WORK(KC2AM),
     *            WORK(KEND2),LWRK2,APROXR12) 
C
      CALL AROUND(' END OF CC_FDB ')
C
      CALL QEXIT('CC_FDB')
C
      RETURN
      END
*=====================================================================*
c /* deck CCDOTRSP */
*=====================================================================*
       SUBROUTINE CCDOTRSP(IDOTS,DOTPROD,IOPT,TYPE,ITRAN,MXTRAN,MXVEC,
     &                     VEC1,VEC2,ISYVEC,WORK,LWORK                )
*---------------------------------------------------------------------*
*
*    Purpose: dot vector on a batch of other vectors
*             which are read from file
*
*             IDOTS   --  arrays of vectors to be dotted on
*             DOTPROD --  array for resulting dot products
*
*             IOPT    --  1 : use only singles
*                         2 : use only doubles
*                         3 : use singles and doubles
*
*             TYPE    --  type of vectors to be dotted on
*             ITRAN   --  index of present vector in the
*                         matrices IDOTS, DOTPROD
*
*             VEC1,VEC2 -- singles and doubles part of vector
*             ISYVEC    -- symmetry of vector
*
*             MXTRAN,MXVEC -- dimensions for IDOTS, DOTPROD
*
*    Written by Christof Haettig, june 1997.
*
*    adapted for R12 by Christian Neiss,  Feb. 2005
*
*=====================================================================*
      IMPLICIT NONE  
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER TYPE*(*)
      INTEGER IOPT, ITRAN, MXVEC, MXTRAN, ISYVEC, LWORK
      INTEGER IDOTS(MXVEC,MXTRAN)

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,MXTRAN)
      REAL VEC1(*), VEC2(*)
      REAL WORK(LWORK)
      REAL CON1, CON2
      REAL ZERO, XNORM1,XNORM2
      REAL DDOT
#else
      DOUBLE PRECISION DOTPROD(MXVEC,MXTRAN)
      DOUBLE PRECISION VEC1(*), VEC2(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION CON1, CON2
      DOUBLE PRECISION ZERO, XNORM1,XNORM2
      DOUBLE PRECISION DDOT
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER MODEL*(10)
      INTEGER KZETA1, KZETA2, KEND, LEND, IVEC, ISYCTR, IZETAV

* external functions:
      INTEGER ILSTSYM


      CALL QENTER('CCDOTRSP')

*---------------------------------------------------------------------*
* allocate memory for vectors to be read from file:
*---------------------------------------------------------------------*
      KZETA1 = 1
      KEND   = KZETA1 + NT1AM(ISYVEC)
      IF (IOPT.EQ.32) THEN
        KZETA2 = KEND
        KEND   = KZETA2 + NTR12AM(ISYVEC)
      ELSE IF (IOPT.GT.1) THEN
        KZETA2 = KEND
        KEND   = KZETA2 + NT2AM(ISYVEC)
      END IF
      LEND   = LWORK - KEND

      IF (LEND .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCDOTRSP.')
      END IF

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CCDOTRSP:',IOPT
         XNORM1 = 0.0d0
         XNORM2 = 0.0d0
         IF (IOPT.EQ.1 .OR. IOPT.EQ.3) 
     &    XNORM1=DDOT(NT1AM(ISYVEC),VEC1,1,VEC1,1)
         IF (IOPT.EQ.2 .OR. IOPT.EQ.3) 
     &    XNORM2=DDOT(NT2AM(ISYVEC),VEC2,1,VEC2,1)
         IF (IOPT.EQ.32) 
     &    XNORM2=DDOT(NTR12AM(ISYVEC),VEC2,1,VEC2,1)
         IF (IOPT.NE.32) THEN
           WRITE (LUPRI,*) 'XNORM1/XNORM2:',XNORM1,XNORM2
           CALL CC_PRP(VEC1,VEC2,ISYVEC,1,1)
         ELSE
           WRITE (LUPRI,*) 'XNORM_R12:',XNORM2
           CALL CC_PRPR12(VEC2,ISYVEC,1,.true.)
         END IF
      END IF

*---------------------------------------------------------------------*
* loop over all vectors to be dotted on:
*---------------------------------------------------------------------*
      IVEC = 1
      DO WHILE (IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

        IZETAV = IDOTS(IVEC,ITRAN)
        ISYCTR = ILSTSYM(TYPE,IZETAV)

        IF (ISYCTR.NE.ISYVEC) THEN
          WRITE (LUPRI,*) 'symmetry mismatch in CCDOTRSP:'
          WRITE (LUPRI,*) 'TYPE:',TYPE
          WRITE (LUPRI,*) 'ISYCTR:',ISYCTR
          WRITE (LUPRI,*) 'ISYVEC:',ISYVEC
          CALL QUIT('symmetry mismatch in CCDOTRSP.')
        END IF

        CALL CC_RDRSP(TYPE,IZETAV,ISYCTR,IOPT,MODEL,
     &                WORK(KZETA1),WORK(KZETA2))
C
        IF (IOPT.EQ.1 .OR. IOPT.EQ.3) THEN
          CON1 = DDOT(NT1AM(ISYCTR),WORK(KZETA1),1,VEC1,1)
        ELSE
          CON1 = ZERO
        END IF

        IF (IOPT.EQ.2 .OR. IOPT.EQ.3) THEN
          IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KZETA2),0.5d0,ISYCTR)
          CON2 = DDOT(NT2AM(ISYCTR),WORK(KZETA2),1,VEC2,1)
        ELSE IF (IOPT.EQ.32) THEN
          CON2 = DDOT(NTR12AM(ISYCTR),WORK(KZETA2),1,VEC2,1)
        ELSE
          CON2 = ZERO
        END IF

        DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + CON1 + CON2

        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CCDOTRSP:',IOPT,DOTPROD(IVEC,ITRAN),
     &          CON1, CON2
         XNORM1 = 0.0d0
         XNORM2 = 0.0d0
         IF (IOPT.EQ.1 .OR. IOPT.EQ.3) 
     &    XNORM1=DDOT(NT1AM(ISYCTR),WORK(KZETA1),1,WORK(KZETA1),1)
         IF (IOPT.EQ.2 .OR. IOPT.EQ.3) THEN 
           XNORM2=DDOT(NT2AM(ISYCTR),WORK(KZETA2),1,WORK(KZETA2),1)
         ELSE IF(IOPT.EQ.32) THEN
           XNORM2=DDOT(NTR12AM(ISYCTR),WORK(KZETA2),1,WORK(KZETA2),1)
         END IF
         WRITE (LUPRI,*) 'TYPE,IZETAV,XNORM:',TYPE(1:2),IZETAV,XNORM1,
     &        XNORM2
        END IF

        IVEC = IVEC + 1

      END DO
      
      CALL QEXIT('CCDOTRSP')

      RETURN
      END 

*---------------------------------------------------------------------*
*               END OF SUBROUTINE CCDOTRSP                            *
*---------------------------------------------------------------------*

*---------------------------------------------------------------------*
c/* Deck CC_BTST */
*=====================================================================*
       SUBROUTINE CC_BTST(WORK,LWORK,APROXR12)
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "dummy.h"
#include "r12int.h"
#include "ccr12int.h"

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_BTST> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .TRUE.)
      INTEGER MXBTRAN
      PARAMETER (MXBTRAN = 2)

      INTEGER LWORK
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL DDOT
      REAL AATRAN1, EATRAN1, AATRAN2, EATRAN2,RDUM(2)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION AATRAN1, EATRAN1, AATRAN2, EATRAN2, RDUM(2)
#endif 

      CHARACTER*(3) LISTA, LISTB, LISTC, APROXR12
      CHARACTER*(8) FILBMA, LABELA
      CHARACTER*(10) MODEL
      INTEGER IOPTRES
      INTEGER IBTRAN(3,MXBTRAN), NBTRAN
      INTEGER IDLSTA, IDLSTB, ISYMA, ISYMB, ISYMAB, IOPT
      INTEGER KRESLT1, KRESLT2, KT1AMPA, KT1AMPB, KT2AMPA, KT2AMPB
      INTEGER KTHETA1, KTHETA2, KEND1, LWRK1, LEND1, LEND2, LEND3
      INTEGER KEND2, KEND3, KETA1, KETA2, KT1AMPC, KT2AMPC, NTAMP
      INTEGER ISYMAC, ISYMC, IDLSTC, IDUM(2)
      INTEGER KTHETAR12, KT12AMPB, KT12AMPA, KRESLTR12
      INTEGER NTR12AB, NTR12A, NTR12B

* external function:
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER ILSTSYM


      CALL QENTER('CC_BTST')


*---------------------------------------------------------------------*
* call B matrix transformation:
*---------------------------------------------------------------------*
      LISTA = 'R1 '
      LISTB = 'R1 '
      IDLSTA = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYMA)
      IDLSTB = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYMB)
      IBTRAN(1,1) = IDLSTA
      IBTRAN(2,1) = IDLSTB
      NBTRAN = 1

      ISYMAB = MULD2H(ISYMA,ISYMB)

      IOPTRES = 1
      FILBMA  = 'CC__BMAT'

      CALL CC_BMAT(IBTRAN,  NBTRAN, LISTA,  LISTB, IOPTRES, 
     &             FILBMA, IDUM, RDUM, 0, .FALSE., WORK, LWORK )

      IF (CCR12) THEN
        NTR12AB = NTR12AM(ISYMAB)
        NTR12A  = NTR12AM(ISYMA)
        NTR12B  = NTR12AM(ISYMB)
      ELSE
        NTR12AB = 0
        NTR12A  = 0
        NTR12B  = 0
      END IF

      KTHETA1 = IBTRAN(3,1)
      KTHETA2 = KTHETA1 + NT1AM(ISYMAB)
      KTHETAR12 = KTHETA2 + NT2AM(ISYMAB)
      KEND1     = KTHETAR12 + NTR12AB

      IF (NSYM.EQ.1 .AND. LOCDBG) THEN
        KT1AMPB = KEND1
        KT2AMPB = KT1AMPB + NT1AM(ISYMB)
        KT12AMPB = KT2AMPB + NT2AM(ISYMB)
        KT1AMPA  = KT12AMPB + NTR12B
        KT2AMPA = KT1AMPA + NT1AM(ISYMA)
        KT12AMPA = KT2AMPA + NT2AM(ISYMA)
        KRESLT1  = KT12AMPA + NTR12A
        KRESLT2 = KRESLT1 + NT1AM(ISYMAB)
        KRESLTR12 = KRESLT2 + NT2AM(ISYMAB)
        KEND1     = KRESLTR12 + NTR12AB
        LEND1   = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_BTST.')
        END IF

        IOPT = 3
        Call CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KT2AMPA))
        Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &                WORK(KT1AMPB),WORK(KT2AMPB))

        IF (CCR12) THEN
          IOPT = 32
          CALL CC_RDRSP(LISTA,IDLSTA,ISYMA,IOPT,MODEL,
     &                  DUMMY,WORK(KT12AMPA))
          Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &                  DUMMY,WORK(KT12AMPB))
        END IF

        ! zero doubles of B and/or C vector:
C       CALL DZERO(WORK(KT1AMPA),NT1AM(ISYMA))
C       CALL DZERO(WORK(KT1AMPB),NT1AM(ISYMB))
C       CALL DZERO(WORK(KT2AMPA),NT2AM(ISYMA))
C       CALL DZERO(WORK(KT2AMPB),NT2AM(ISYMB))
        CALL DZERO(WORK(KRESLT1),NT1AM(ISYMAB)+NT2AM(ISYMAB)+
     &               NTR12AB)
C       IPRINT  = 5

        CALL CC_FDB(NT1AM(ISYMAB),NT2AM(ISYMAB),NTR12AB,
     >              WORK(KT1AMPA), WORK(KT1AMPB), WORK(KRESLT1),
     >              WORK(KEND1), LEND1, APROXR12)

C       IPRINT  = 0

        IF (.TRUE.) THEN
          WRITE (LUPRI,*) 'LISTA, IDLSTA, ISYMA:',LISTA,IDLSTA,ISYMA
          WRITE (LUPRI,*) 'LISTB, IDLSTB, ISYMB:',LISTB,IDLSTB,ISYMB
          WRITE (LUPRI,*) 'ISYMAB:',ISYMAB
          WRITE (LUPRI,*)
          WRITE (LUPRI,*) 'finite difference Theta vector:'
          Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMAB,1,1)
          IF (CCR12) THEN
            Call CC_PRPR12(WORK(KRESLTR12),ISYMAB,1,.TRUE.)
          END IF
          WRITE (LUPRI,*) 'analytical Theta vector:'
          Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,1,1)
          IF (CCR12) THEN
            Call CC_PRPR12(WORK(KTHETAR12),ISYMAB,1,.TRUE.)
          END IF
          WRITE(LUPRI,*) 'Norm of analytical Theta vector: ',
     &      DSQRT(DDOT(NT1AM(ISYMAB)+NT2AM(ISYMAB)+NTR12AB,
     &                 WORK(KTHETA1),1,WORK(KTHETA1),1))
        END IF

        Call DAXPY(NT1AM(ISYMAB),-1.0d0,WORK(KTHETA1),1,
     &                                  WORK(KRESLT1),1)
        IF (.NOT.(CCS.OR.CCSTST)) THEN
          Call DAXPY(NT2AM(ISYMAB),-1.0d0,WORK(KTHETA2),1,
     &                                    WORK(KRESLT2),1)
        ELSE
          Call DZERO(WORK(KRESLT2),NT2AM(ISYMAB))
        END IF

        IF (CCR12) THEN
          CALL DAXPY(NTR12AB,-1.0d0,WORK(KTHETAR12),1,
     &                              WORK(KRESLTR12),1)
        END IF

        WRITE (LUPRI,*) 'Norm of difference between analytical THETA '
     >           // 'vector and the numerical result:'
        WRITE (LUPRI,*) 'singles excitation part:',
     >   DSQRT(DDOT(NT1AM(ISYMAB),WORK(KRESLT1),1,WORK(KRESLT1),1))
        WRITE (LUPRI,*) 'double excitation part: ',
     >   DSQRT(DDOT(NT2AM(ISYMAB),WORK(KRESLT2),1,WORK(KRESLT2),1))
        IF (CCR12) THEN
          WRITE (LUPRI,*) 'R12 double excitation part: ',
     &      DSQRT(DDOT(NTR12AB,WORK(KRESLTR12),1,
     &                         WORK(KRESLTR12),1))
        END IF

        WRITE (LUPRI,*) 'difference vector:'
        Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMAB,1,1)
        IF (CCR12) THEN
          Call CC_PRPR12(WORK(KRESLTR12),ISYMAB,1,.TRUE.)
        END IF

        CALL FLSHFO(LUPRI)


      ELSE IF (NSYM.NE.1 .AND. LOCDBG) THEN

        WRITE (LUPRI,*) 
     &        'CC_BTST> can not calculate finite difference B matrix'
        WRITE (LUPRI,*) 'CC_BTST> with symmetry.'

        WRITE (LUPRI,*) 'analytical Theta vector:'
        Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMAB,1,1)
        IF (CCR12) THEN
          Call CC_PRPR12(WORK(KTHETAR12),ISYMAB,1,.TRUE.)
        END IF
        WRITE(LUPRI,*) 'Norm of analytical Theta vector: ',
     &    DSQRT(DDOT(NT1AM(ISYMAB)+NT2AM(ISYMAB)+NTR12AB,
     &               WORK(KTHETA1),1,WORK(KTHETA1),1))
        CALL FLSHFO(LUPRI)

      END IF

*=====================================================================*
* test A{A} transformation:
*=====================================================================*
C     This part is NOT adapted for CC-R12 yet!
      IF (.FALSE.) THEN
      WRITE (LUPRI,'(/,/5X,A)') 'TEST A{A} TRANSFORMATION:'
      IF (CCR12) CALL QUIT('Not adapted for CC-R12')

      LABELA = 'ZDIPLEN '
      LISTC  = 'L1 '
      IDLSTC = IL1ZETA('ZDIPLEN ',.FALSE.,0.0D0,1)
      ISYMC  = ILSTSYM(LISTC,IDLSTC)

      KRESLT1 = 1
      KRESLT2 = KRESLT1 + NT1AM(ISYMAB)
      KT1AMPC = KRESLT2 + NT2AM(ISYMAB)
      KT2AMPC = KT1AMPC + NT1AM(ISYMC)
      KEND1   = KT2AMPC + NT2AM(ISYMC)
      LEND1   = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_BTST.')
        END IF

      CALL CCCR_AA(LABELA,ISYMA,LISTB,IDLSTB,DUMMY,WORK,LWORK)


      IOPT = 3
      Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KT1AMPC),WORK(KT2AMPC))

      IF (ISYMC.NE.ISYMAB) THEN
        CALL QUIT('SYMMETRY MISMATCH IN CC_BTST.')
      END IF
      
      AATRAN1 = DDOT(NT1AM(ISYMC),WORK(KRESLT1),1,WORK(KT1AMPC),1)
      IF (.NOT.CCS) THEN
        AATRAN2 = DDOT(NT2AM(ISYMC),WORK(KRESLT2),1,WORK(KT2AMPC),1)
      END IF


      ISYMAC = MULD2H(ISYMA,ISYMC)

      KETA1 = KEND1
      KETA2 = KETA1 + NT1AM(ISYMAC)
      KEND2 = KETA2 + NT2AM(ISYMAC)
      LEND2 = LWORK - KEND2

      KT1AMPB = KEND2
      KT2AMPB = KT1AMPB + NT1AM(ISYMB)
      KEND3   = KT2AMPB + NT2AM(ISYMB)
      LEND3   = LWORK - KEND3

        IF (LEND3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_BTST.')
        END IF

      CALL CC_ETAC(ISYMA,LABELA,WORK(KEND1),
     &             LISTC,IDLSTC,0,DUMMY,WORK(KEND2),LEND2)

      IOPT = 3
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KT1AMPB),WORK(KT2AMPB))
      
      EATRAN1 = DDOT(NT1AM(ISYMB),WORK(KETA1),1,WORK(KT1AMPB),1)
      IF (.NOT.CCS) THEN
        EATRAN2 = DDOT(NT2AM(ISYMB),WORK(KETA2),1,WORK(KT2AMPB),1)
      END IF

      WRITE (LUPRI,*) 'comparison of the results:'
      WRITE (LUPRI,*) 'C x AATRAN(A,B):', AATRAN1+AATRAN2, AATRAN1,
     &     AATRAN2
      WRITE (LUPRI,*) 'EATRAN(C,A) x B:', EATRAN1+EATRAN2, EATRAN1,
     &     EATRAN2

      END IF
    
      CALL QEXIT('CC_BTST')

      RETURN
      END 
*=====================================================================*
*=====================================================================*
C  /* Deck ccb_cdsort */
      SUBROUTINE CCB_CDSORT(XINT,ISYDIS,DDRHF,XLAMDP,ISYXLP,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculated presorted one-index transformed integrals
*              for D intermediate:
*
*           DDRHF(k,alp bet) = 2 (alp bet|k del) - (k bet|alp del)
*
*     Written by Christof Haettig November 1998
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"
#include "symsq.h"

      INTEGER ISYALP, ISYBET, ISYGAM, ISYRHF, ISYDIS, ISYXLP, LWORK
      INTEGER ISYMAB, ISYMBK, ISYMGK, NRHFK, NBASA, KSCR1, KSCR2
      INTEGER KOFF0, KOFF1, KOFF2, KOFF3, KOFF4, KOFF5, KEND1, LWRK1
      INTEGER ICOUNT, NBSRHF(8), IBSRHF(8,8), ISYM, ISYMAK, ISYMK
     

#if defined (SYS_CRAY)
      REAL XINT(*), DDRHF(*), XLAMDP(*), WORK(LWORK)
      REAL TWO, ONE, ZERO
#else
      DOUBLE PRECISION XINT(*), DDRHF(*), XLAMDP(*), WORK(LWORK)
      DOUBLE PRECISION TWO, ONE, ZERO
#endif
      PARAMETER (TWO = 2.0D0, ONE = 1.0D0, ZERO = 0.0D0)
C
      CALL QENTER('CCB_CDSORT')
C
C     --------------------------------------
C     precalculate symmetry array for DDRHF:
C     --------------------------------------
C
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO
C
      ISYRHF = MULD2H(ISYDIS,ISYXLP)
C
      CALL DZERO(DDRHF,NBSRHF(ISYRHF))
C
      DO ISYGAM = 1, NSYM
C 
         ISYMAB = MULD2H(ISYGAM,ISYDIS)
         ISYMBK = MULD2H(ISYMAB,ISYXLP)
C
         KSCR1 = 1
         KSCR2 = KSCR1 + N2BST(ISYMAB)
         KEND1 = KSCR2 + NT1AO(ISYMBK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in CCB_CDSORT.')
         END IF
C
         DO G = 1, NBAS(ISYGAM)
C
           KOFF0 = IDSAOG(ISYGAM,ISYDIS) + NNBST(ISYMAB)*(G-1) + 1
           CALL CCSD_SYMSQ(XINT(KOFF0),ISYMAB,WORK(KSCR1))
C
           DO ISYALP = 1, NSYM

             ISYBET = MULD2H(ISYMAB,ISYALP) 
             ISYMK  = MULD2H(ISYALP,ISYXLP)
             ISYMGK = MULD2H(ISYGAM,ISYMK)
C
C            -------------------------------
C            transform the alpha index to k:
C            -------------------------------
C
             KOFF1 = IGLMRH(ISYALP,ISYMK) + 1
             KOFF2 = KSCR1 + IAODIS(ISYALP,ISYBET)
             KOFF3 = KSCR2 + IT1AO(ISYBET,ISYMK)
             NBASA = MAX(NBAS(ISYALP),1)
             NRHFK = MAX(NRHF(ISYMK),1)

             CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYBET),NBAS(ISYALP),
     &                  ONE,XLAMDP(KOFF1),NBASA,WORK(KOFF2),NBASA,
     &                  ZERO,WORK(KOFF3),NRHFK)
C
C            --------------------------
C            store as DDRHF(gam k;bet):
C            --------------------------
C
             DO B = 1, NBAS(ISYBET)

                KOFF4 = KSCR2 + IT1AO(ISYBET,ISYMK) + NRHF(ISYMK)*(B-1)
                KOFF5 = IBSRHF(ISYMGK,ISYBET) + NT1AO(ISYMGK)*(B-1)
     &                + IT1AO(ISYGAM,ISYMK) + G

                CALL DCOPY(NRHF(ISYMK),WORK(KOFF4),1,
     &                                 DDRHF(KOFF5),NBAS(ISYGAM))
             END DO
C
           END DO
C
         END DO
C
      END DO
C
      CALL QEXIT('CCB_CDSORT')
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck cc_cdb */
      SUBROUTINE CC_CDB(DDRHF, ISYRHF, IDEL, ISYDEL, LUD, DFIL, IV,
     &                  XLAMDP, XLAMDH, XLAMPC, XLAMHC, ISYXLC, 
     &                  IOPTR, FACTR, RIM, WORK, LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: calculate the D intermediate in the B matrix transf.
*
*     Written by Christof Haettig November 1998
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"

      CHARACTER*(*) DFIL
      INTEGER LWORK, ISYXLC, ISYRHF, LUD, IV, IOPTR, IDEL, ISYDEL

#if defined (SYS_CRAY)
      REAL DDRHF(*), RIM(*), WORK(LWORK)
      REAL XLAMDP(*), XLAMDH(*), XLAMPC(*), XLAMHC(*)
      REAL FACTR, ONE, ZERO, DDOT, XNORM
#else
      DOUBLE PRECISION DDRHF(*), RIM(*), WORK(LWORK)
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), XLAMPC(*), XLAMHC(*)
      DOUBLE PRECISION FACTR, ONE, ZERO, DDOT, XNORM
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INTEGER NBSRHF(8), IBSRHF(8,8), ISYM, ISYALK, ISYBET, ICOUNT
      INTEGER ISYMI, ISYMA, ISYMAI, KOFF1, KOFF2, KOFF3, KOFF4, KOFF5
      INTEGER NAI, ISYAIK, KSCR1, KSCR2, KEND1, LWRK1, KOFF6, ISYALP
      INTEGER ISYMK, NT1AK, NBASB, NBASA, NRHFK, IADR
C
      CALL QENTER('CC_CDB')
C
C     --------------------------------------
C     precalculate symmetry array for BSRHF:
C     --------------------------------------
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYALK = 1, NSYM
           ISYBET = MULD2H(ISYALK,ISYM)
           IBSRHF(ISYALK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYALK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO
C
      ISYAIK = MULD2H(ISYRHF,ISYXLC)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + MAX(NT2BGD(ISYAIK),NT2BGD(ISYRHF))
      KEND1 = KSCR2 + NT2BCD(ISYAIK)
      LWRK1 = LWORK - KEND1
C
      CALL DZERO(WORK(KSCR2),NT2BCD(ISYAIK))
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CC_CDB.')
      END IF
C
      DO ISYALK = 1, NSYM

         ISYBET = MULD2H(ISYALK,ISYRHF)
C
C        -------------------------------------------------
C        transform beta index to i and alpha index to a^C:
C        -------------------------------------------------
C
         ISYMI  = ISYBET
C
         KOFF1 = IBSRHF(ISYALK,ISYBET) + 1
         KOFF2 = IGLMRH(ISYBET,ISYMI)  + 1
         KOFF3 = KSCR1 + IT2BGT(ISYMI,ISYALK) 

         NT1AK = MAX(NT1AO(ISYALK),1)
         NBASB = MAX(NBAS(ISYBET),1)

         CALL DGEMM('N','N',NT1AO(ISYALK),NRHF(ISYMI),NBAS(ISYBET),
     &              ONE,DDRHF(KOFF1),NT1AK,XLAMDH(KOFF2),NBASB,
     &              ZERO,WORK(KOFF3),NT1AK)

         DO I = 1, NRHF(ISYMI)

           DO ISYALP = 1, NSYM

             ISYMK  = MULD2H(ISYALK,ISYALP)
             ISYMA  = MULD2H(ISYALP,ISYXLC)
             ISYMAI = MULD2H(ISYMA,ISYMI)
 
             NAI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + 1

             KOFF5 = IGLMVI(ISYALP,ISYMA) + 1
             KOFF6 = KSCR2 + IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI-1)
             KOFF4 = KSCR1 + IT2BGT(ISYMI,ISYALK) + 
     &                       NT1AO(ISYALK)*(I-1) + IT1AO(ISYALP,ISYMK)

             NBASA = MAX(NBAS(ISYALP),1)
             NRHFK = MAX(NRHF(ISYMK),1)

             CALL DGEMM('T','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYALP),
     &                  ONE,WORK(KOFF4),NBASA,XLAMPC(KOFF5),NBASA,
     &                  ONE,WORK(KOFF6),NRHFK)

           END DO

         END DO

      END DO
C
      DO ISYALK = 1, NSYM

         ISYBET = MULD2H(ISYALK,ISYRHF)
C
C        -------------------------------------------------
C        transform beta index to i^C and alpha index to a:
C        -------------------------------------------------
C
         ISYMI  = MULD2H(ISYBET,ISYXLC)
C
         KOFF1 = IBSRHF(ISYALK,ISYBET) + 1
         KOFF2 = IGLMRH(ISYBET,ISYMI)  + 1
         KOFF3 = KSCR1 + IT2BGT(ISYMI,ISYALK) 

         NT1AK = MAX(NT1AO(ISYALK),1)
         NBASB = MAX(NBAS(ISYBET),1)

         CALL DGEMM('N','N',NT1AO(ISYALK),NRHF(ISYMI),NBAS(ISYBET),
     &              ONE,DDRHF(KOFF1),NT1AK,XLAMHC(KOFF2),NBASB,
     &              ZERO,WORK(KOFF3),NT1AK)

         DO I = 1, NRHF(ISYMI)

           DO ISYALP = 1, NSYM

             ISYMK  = MULD2H(ISYALK,ISYALP)
             ISYMA  = ISYALP
             ISYMAI = MULD2H(ISYMA,ISYMI)

             NAI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + 1

             KOFF5 = IGLMVI(ISYALP,ISYMA) + 1
             KOFF6 = KSCR2 + IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI-1)
             KOFF4 = KSCR1 + IT2BGT(ISYMI,ISYALK) + 
     &                       NT1AO(ISYALK)*(I-1) + IT1AO(ISYALP,ISYMK)

             NBASA = MAX(NBAS(ISYALP),1)
             NRHFK = MAX(NRHF(ISYMK),1)

             CALL DGEMM('T','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYALP),
     &                  ONE,WORK(KOFF4),NBASA,XLAMDP(KOFF5),NBASA,
     &                  ONE,WORK(KOFF6),NRHFK)

           END DO

         END DO
 
      END DO
C
C     -------------------------------
C     write the intermediate to disk:
C     -------------------------------
C
      IADR   = IT2DLR(IDEL,IV) + 1
      ISYAIK = MULD2H(ISYRHF,ISYXLC)
C
C     XNORM = DDOT(NT2BCD(ISYAIK),WORK(KSCR2),1,WORK(KSCR2),1)
C     WRITE (LUPRI,*) 'IDEL,XNORM:',IDEL,XNORM
C
      CALL PUTWA2(LUD,DFIL,WORK(KSCR2),IADR,NT2BCD(ISYAIK))
C       
C     --------------------------------------------------
C     calculate contribution to R intermediate as trace:
C     --------------------------------------------------
C
      IF (IOPTR.EQ.1 .AND. NT2BCD(ISYAIK).GT.0 ) THEN

         D = IDEL - IBAS(ISYDEL)

         DO ISYMAI = 1, NSYM

            ISYMK = MULD2H(ISYMAI,ISYAIK)
            ISYMI = ISYMK
            ISYMA = MULD2H(ISYMAI,ISYMI)

            DO I = 1, NRHF(ISYMI)
            DO A = 1, NVIR(ISYMA)

               NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I-1)   + A
               KOFF1 = IEMAT1(ISYMA,ISYDEL) + NVIR(ISYMA)*(D-1)   + A
               KOFF3 = IT2BCT(ISYMI,ISYMAI) + NRHF(ISYMI)*(NAI-1) + I
               
               RIM(KOFF1) = RIM(KOFF1) + FACTR * WORK(KSCR2-1+KOFF3)

            END DO
            END DO

         END DO 

      END IF
C
      CALL QEXIT('CC_CDB')
C
      RETURN
      END
*=====================================================================*
c /* deck cc_aibj */
*=====================================================================*
       SUBROUTINE CC_AIBJ( X0INT,   ISY0ALBE, X1INT,   ISY1ALBE, 
     &                     IDEL,    IGAM,     X0AIBJ,  X1AIBJ, 
     &                     XLAMDHA, ISYXLA,   XLAMDHB, ISYXLB,
     &                     XLAMDP0, ISYXL0,   WORK,    LWORK,   
     &                     IOPT,    LDERIV,   LRELAX             )
*---------------------------------------------------------------------*
*
*   Purpose: generalized transformation to (ai|bj) integrals
*            for the two-index (**|gam del) approach
*            assumes three-index array XAIBJ in core
*
*            this routine transforms the indeces ia and j, the 
*            transformation of the delta index to b has to be done
*            from the outside.  
*
*            The (ai|bj) integrals are calculated by 
*            transforming the gamma index with XLAMDH matrices, 
*            this assumes symmetric AO integrals 
*                   --> factor (-1) for London orbitals i guess...
*
*            X0INT, X1IAJB : usual integrals
*            X1INT, X1IAJB : derivative integrals
*
*
*            IOPT=0: (a i^A|del j^B) as for F term in energy code
*
*            IOPT=1: not used
*
*            IF LDERIV=.TRUE. transform also the derivative integrals:
*               not yet implemented...
*
*            IF LRELAX=.TRUE. include relaxtion contribution to the
*            derivative integrals:
*               not yet implemented...
*
*                  i^A  transform. with XLAMDHA with sym. ISYXLA
*                  j^B  transform. with XLAMDHB with sym. ISYXLB
*                  a    transform. with XLAMDP0 with sym. ISYXL0
*             
*    Written by Christof Haettig, October 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccisao.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)


#if defined (SYS_CRAY)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      LOGICAL LDERIV, LRELAX
      INTEGER IDEL, IGAM, ISY0ALBE, ISY1ALBE, LWORK, IOPT, IDUMMY
      INTEGER ISYXL0, ISYXLA, ISYXLB
      
#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDHA(*), XLAMDHB(*)
      REAL X0INT(*),   X1INT(*),   X0AIBJ(*),  X1AIBJ(*) 
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDHA(*), XLAMDHB(*)
      DOUBLE PRECISION X0INT(*),   X1INT(*),   X0AIBJ(*),  X1AIBJ(*)
      DOUBLE PRECISION WORK(LWORK)
#endif

      INTEGER ISYMAI, ISYGAM, ISYALP, ISYBET
      INTEGER KSCR1, KSCR2, KSCR4, KEND1, LWRK1
      INTEGER KOFF1, KOFF2, KOFF4, KLAMD
      INTEGER NBASA, NBASB, NVIRA, ISYMA, ISYMI

      CALL QENTER('CC_AIBJ')

*---------------------------------------------------------------------*
*     work space allocation:
*
*     KSCR1 --  I^{del,gam}(alp bet)
*     KSCR2 --  I^{del,gam}(i bet)
*     KSCR4 --  I^{del,gam}(i a)
*
*---------------------------------------------------------------------*
      KSCR1   = 1
      KSCR2   = KSCR1 + N2BST(ISY0ALBE)
      KSCR4   = KSCR2 + NBAST*NRHFT
      KEND1   = KSCR4 + NVIRT*NRHFT

      LWRK1   = LWORK - KEND1

      IF ( LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_AIBJ.')
      END IF

*---------------------------------------------------------------------*
*     square integral matrix up 
*---------------------------------------------------------------------*

      CALL CCSD_SYMSQ(X0INT,ISY0ALBE,WORK(KSCR1))

*---------------------------------------------------------------------*
*     transform alpha index to I using XLAMDHA
*      -- store (i bet|gam del) in SCR2 
*---------------------------------------------------------------------*
      KOFF2 = KSCR2

      DO ISYMI = 1, NSYM
              
        ISYALP = MULD2H(ISYXLA,ISYMI)
        ISYBET = MULD2H(ISYALP,ISY0ALBE)

        KOFF1 = KSCR1 + IAODIS(ISYALP,ISYBET) 
        KLAMD = IGLMRH(ISYALP,ISYMI) + 1

        NBASA = MAX(NBAS(ISYALP),1)
        NBASB = MAX(NBAS(ISYBET),1)

        CALL DGEMM('T','N',NBAS(ISYBET),NRHF(ISYMI),NBAS(ISYALP),
     *             ONE,WORK(KOFF1),NBASA,XLAMDHA(KLAMD),
     *             NBASA,ZERO,WORK(KOFF2),NBASB)

        KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)
          
      END DO

*---------------------------------------------------------------------*
*     transform beta index to a using XLAMDP0 
*      -- store (i a|gam del) in SCR4 
*---------------------------------------------------------------------*
      KOFF2 = KSCR2

      DO ISYMI = 1, NSYM
              
        ISYALP = MULD2H(ISYXLA,ISYMI)
        ISYBET = MULD2H(ISYALP,ISY0ALBE)
        ISYMA  = MULD2H(ISYXL0,ISYBET)

        KLAMD = IGLMVI(ISYBET,ISYMA) + 1
        KOFF4 = KSCR4 + IT1AM(ISYMA,ISYMI)

        NBASB = MAX(NBAS(ISYBET),1)
        NVIRA = MAX(NVIR(ISYMA),1)

        CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBET),
     *             ONE,XLAMDP0(KLAMD),NBASB,WORK(KOFF2),
     *             NBASB,ZERO,WORK(KOFF4),NVIRA)

        KOFF2 = KOFF2 + NBAS(ISYBET)*NRHF(ISYMI)

      END DO

*---------------------------------------------------------------------*
*     Add the contribution to the result XIAJB vector
*     transform thereby gamma to j using XLAMDH1 and XLAMDH2
*---------------------------------------------------------------------*
      ISYGAM  = ISAO(IGAM)

C     --------------------------
C     add (i a|j del) to X1IAJB:
C     --------------------------
      ISYMAI = MULD2H(ISY0ALBE,MULD2H(ISYXL0,ISYXLA))

      CALL CC_IAJB1(IGAM, WORK(KSCR4), ISYMAI, ISYGAM,
     &              XLAMDHB, ISYXLB, X0AIBJ, .FALSE., IDUMMY)

      CALL QEXIT('CC_AIBJ')

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_AIBJ                           *
*=====================================================================*
*=====================================================================*
C  /* Deck ccbpre1int */
      SUBROUTINE CCBPRE1INT(INTMED1,NINT1,IOFFCD,IADRBFD,
     &                      CBAFIL,DBAFIL,FNBFD,
     &                      XLAMDP0,XLAMDH0,
     &                      TIMIO,TIMC,TIMD,TIMBF,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: precalculate some intermediates for B matrix transform.
*              which depend only on one amplitude response vector:
*
*              CBAR, DBAR, and the BF density 
*
*              the results are written to direct acces files
*
*     Written by Christof Haettig November 1998
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "cciccset.h"
#include "second.h"

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER LUCBAR, LUDBAR, LUBFD

      CHARACTER*(*) CBAFIL, DBAFIL, FNBFD
      INTEGER LWORK, IDLSTR, NINT1
      INTEGER IADRBFD(*), ISTARTPQ, IOFFCD(0:NINT1), INTMED1(2,NINT1)

#if defined (SYS_CRAY)
      REAL WORK(LWORK), XLAMDP0(*), XLAMDH0(*)
      REAL TIMC, TIMD, TIMBF, TIMIO, DTIME
      REAL TWO, ONE, ZERO, DDOT, DUMMY(2)
#else
      DOUBLE PRECISION WORK(LWORK), XLAMDP0(*), XLAMDH0(*)
      DOUBLE PRECISION TIMC, TIMD, TIMBF, TIMIO, DTIME
      DOUBLE PRECISION TWO, ONE, ZERO, DDOT, DUMMY(2)
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CHARACTER MODEL*(10), LIST*(3)
      INTEGER IOPT,IINT1,ISTARTBFD,IDLST,ISYMA,KT2AMSQ,KXPACK,IDUMMY
      INTEGER KT1AMPA, KT2AMPA, KXLAMHA, KXLAMPA, KXIAJB, KEND1, LWRK1
      INTEGER ILSTSYM

      CALL QENTER('CCBPRE1INT')

*---------------------------------------------------------------------*
* test CC model, open files and do some initializations:
*---------------------------------------------------------------------*
      IF (CCS .OR. CC2) THEN
         CALL QEXIT('CCBPRE1INT')
         RETURN
      ENDIF

      LUBFD  = -1
      LUCBAR = -1
      LUDBAR = -1
      CALL WOPEN2(LUCBAR,CBAFIL,64,0)
      CALL WOPEN2(LUDBAR,DBAFIL,64,0)
      CALL WOPEN2(LUBFD, FNBFD, 64,0)

      ISTARTBFD = 1
      IOFFCD(1) = 0

*---------------------------------------------------------------------*
* loop over all amplitude vectors and compute the intermediates:
*---------------------------------------------------------------------*
      DO IINT1 = 1, NINT1
         LIST  = VTABLE(INTMED1(2,IINT1))
         IDLST = INTMED1(1,IINT1)
         ISYMA = ILSTSYM(LIST,IDLST)

         KT1AMPA = 1
         KT2AMPA = KT1AMPA + NT1AM(ISYMA)
         KXLAMHA = KT2AMPA + MAX(NT2AM(ISYMA),NT2AM(ISYM0))
         KXLAMPA = KXLAMHA + NGLMDT(ISYMA)
         KXIAJB  = KXLAMPA + NGLMDT(ISYMA)
         KEND1   = KXIAJB  + MAX(NT2SQ(ISYM0),NT2SQ(ISYMA))
         LWRK1   = LWORK   - KEND1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCBPRE1INT.')
         END IF

* reuse integral array for squared amplitudes and 
* amplitude array for packed integrals:
         KT2AMSQ = KXIAJB
         KXPACK  = KT2AMPA

C        -----------------------------------------------------------
C        read packed integrals, square up and read packed amplitudes
C        -----------------------------------------------------------
         DTIME = SECOND()

         CALL CCG_RDIAJB(WORK(KXPACK),NT2AM(ISYM0))
         CALL CC_T2SQ(WORK(KXPACK),WORK(KXIAJB),ISYM0)

         IOPT = 3
         CALL CC_RDRSP(LIST,IDLST,ISYMA,IOPT,MODEL,
     *                 WORK(KT1AMPA),WORK(KT2AMPA))
         CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,ISYMA)

         TIMIO = TIMIO + SECOND() - DTIME

C        -----------------------------------------------------------
C        calculate CBAR intermediate:
C        -----------------------------------------------------------
         DTIME = SECOND()
         IOPT  = 2
         CALL CCB_CDBAR('C',WORK(KXIAJB),ISYM0,WORK(KT2AMPA),ISYMA,
     &                  DUMMY,ISYMA, WORK(KEND1),LWRK1,
     &                  CBAFIL,LUCBAR,IOFFCD(IINT1),IOPT)
         TIMC = TIMC + SECOND() - DTIME
            
C        -----------------------------------------------------------
C        calculate DBAR intermediate:
C        -----------------------------------------------------------
         DTIME = SECOND()
         IOPT  = 2
         CALL CCB_CDBAR('D',WORK(KXIAJB),ISYM0,WORK(KT2AMPA),ISYMA,
     &                  DUMMY,ISYMA, WORK(KEND1),LWRK1,
     &                  DBAFIL,LUDBAR,IOFFCD(IINT1),IOPT)
         TIMD = TIMD + SECOND() - DTIME

C        ---------------------------------------------------------
C        increment offset for CBAR & DBAR intermediate:
C        ---------------------------------------------------------
         IF (IINT1.LT.NINT1) THEN
            IOFFCD(IINT1+1) = IOFFCD(IINT1) + NT2SQ(ISYMA)
         END IF

C        ------------------------------------------------------
C        reread response amplitudes, scale and square T2 part
C        ------------------------------------------------------
         DTIME = SECOND()

         IOPT  = 3
         CALL CC_RDRSP(LIST,IDLST,ISYMA,IOPT,MODEL,
     *                 WORK(KT1AMPA),WORK(KT2AMPA))
         CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,ISYMA)
         CALL CC_T2SQ(WORK(KT2AMPA),WORK(KT2AMSQ),ISYMA)

         TIMIO = TIMIO + SECOND() - DTIME
     
C        ------------------------------------------------------
C        calculate response lambda matrices:
C        ------------------------------------------------------
         CALL CCLR_LAMTRA(XLAMDP0,WORK(KXLAMPA),
     *                    XLAMDH0,WORK(KXLAMHA),
     *                    WORK(KT1AMPA),ISYMA)

C        ---------------------------------------------------------
C        calculate effective density for BF term and store on disk
C        ---------------------------------------------------------
         DTIME = SECOND()
         IOPT  = 3
         CALL CC_BFDEN(WORK(KT2AMSQ),ISYMA, DUMMY, IDUMMY,
     *                 XLAMDH0,      ISYM0, XLAMDH0,ISYM0,
     *                 WORK(KXLAMHA),ISYMA, DUMMY,  IDUMMY,
     *                 FNBFD,  LUBFD,IADRBFD, ISTARTBFD,
     *                 IINT1,  IOPT, .FALSE., WORK(KEND1),LWRK1)
         TIMBF  = TIMBF  + SECOND() - DTIME

      END DO

*---------------------------------------------------------------------*
*     that's it; close files and return:
*---------------------------------------------------------------------*
      CALL WCLOSE2(LUCBAR,CBAFIL,'KEEP')
      CALL WCLOSE2(LUDBAR,DBAFIL,'KEEP')
      CALL WCLOSE2(LUBFD, FNBFD, 'KEEP')

      CALL QEXIT('CCBPRE1INT')

      RETURN
      END
*=====================================================================*
c/* Deck CC_R12BMAT */
*=====================================================================*
       SUBROUTINE CC_R12BMAT(THETA1, THETA2, THETAR12, ISYRES,
     &                       LISTA, IDLSTA, T1AMPA, ISYMA,
     &                       LISTB, IDLSTB, T1AMPB, ISYMB,
     &                       LAMDPA, LAMDHA, LAMDPB, LAMDHB,
     &                       LAMP0, LAMH0, WORK, LWRK)
*---------------------------------------------------------------------*
*
*    Purpose: calculate R12 contributions for B-matrix transformations
*
*    C. Neiss  june 2005
*
*=====================================================================*
       implicit none
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "dummy.h"
#include "r12int.h"
#include "ccr12int.h"
#include "ccfield.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LV, LVAJKL, LVIJKL
      INTEGER LWRK, ISYRES, IDLSTA, ISYMA, IDLSTB, ISYMB, KEND1, LWRK1
      INTEGER ISYM1, ISYM2, IDLST1, KEIM, KSCR, KT1AMP, KVAJKL, KVIJKL
      INTEGER KEND0, KT12AMP, KXINTTRI, KXINTSQ, KFIELDAO, IFLD
      INTEGER LUNIT, IAN, IOPT 
      INTEGER KVABKL, KEND2, LWRK2
      CHARACTER LISTA*3, LISTB*3, CDUMMY*3, LIST1*3

#if defined (SYS_CRAY)
      REAL WORK(LWRK), THETA1(*), THETA2(*), THETAR12(*),
     &     T1AMPA(*),T1AMPB(*),
     &     LAMDPA(*), LAMDHA(*), LAMDPB(*), LAMDHB(*),
     &     LAMP0(*), LAMH0(*) 
      REAL TIM1, TIM2, TIM3
      REAL ONE
#else
      DOUBLE PRECISION WORK(LWRK), THETA1(*), THETA2(*), THETAR12(*),
     &                 T1AMPA(*),T1AMPB(*),
     &                 LAMDPA(*), LAMDHA(*), LAMDPB(*), LAMDHB(*),
     &                 LAMP0(*), LAMH0(*)
      DOUBLE PRECISION TIM1, TIM2, TIM3 
      DOUBLE PRECISION ONE
#endif      
      PARAMETER (ONE = 1.0D0)

      CALL QENTER('CC_R12BMAT')
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'Entered CC_R12BMAT'
        CALL FLSHFO(LUPRI)
      ENDIF
C
      IF (ISYRES.NE.MULD2H(ISYMA,ISYMB)) THEN
        CALL QUIT('Symmetry error in CC_R12BMAT')
      ENDIF
C
      IF (CC2) THEN
C
C     --------------------------------------------------------------
C     calculate G'-Term Singles contributions:
C     do first E-intermediate calculation, then contract with
C     singles Lagrangian multipliers and add to conventional term
C     --------------------------------------------------------------
C
      KEIM  = 1
      KSCR  = KEIM + MAX(NMATIJ(ISYMA),NMATIJ(ISYMB))
      KT1AMP= KSCR + MAX(NMATAB(ISYMA),NMATAB(ISYMB))
      KEND1 = KT1AMP + MAX(NT1AM(ISYMA),NT1AM(ISYMB))
      LWRK1 = LWRK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_R12BMAT')
      END IF
C
      DO I = 1, 2
C       E(R12)-Intermediate:
C
        IF (I.EQ.1) THEN
          ISYM1 = ISYMA
          ISYM2 = ISYMB
          LIST1 = LISTA
          IDLST1 = IDLSTA
          CALL DCOPY(NT1AM(ISYMB),T1AMPB,1,WORK(KT1AMP),1)
        ELSE IF (I.EQ.2) THEN
          ISYM1 = ISYMB
          ISYM2 = ISYMA
          LIST1 = LISTB
          IDLST1 = IDLSTB
          CALL DCOPY(NT1AM(ISYMA),T1AMPA,1,WORK(KT1AMP),1)
        END IF
C
        CALL DZERO(WORK(KEIM),NMATIJ(ISYM1))
        CALL CCRHS_EINTP(WORK(KEIM),LAMP0,WORK(KEND1),LWRK1,
     &                   2,ISYM1,IDUMMY,IDUMMY,IDUMMY,LIST1,IDLST1)
C
        CALL DZERO(WORK(KSCR),NMATAB(ISYM1))
        CALL CCLR_E1C1(THETA1,WORK(KT1AMP),WORK(KSCR),WORK(KEIM),
     &                 WORK(KEND1),LWRK1,ISYM2,ISYM1,'N')
C
        IF (LOCDBG) THEN
          WRITE(LUPRI,*) "E Intermediates in CC_R12BMAT:"
          CALL CC_PREI(WORK(KSCR),WORK(KEIM),ISYM1,1)
        END IF
C
      END DO
C
      END IF
C
C     --------------------------------------------------------------
C     calculate F'-Term r12 doubles contribution
C     --------------------------------------------------------------
C
      KVIJKL = 1
      KEND0  = KVIJKL + NTR12SQ(ISYRES)
      KVAJKL = KEND0
      KSCR   = KVAJKL + NVAJKL(ISYMA)
      KEND1  = KSCR + NTR12AM(ISYRES)
      LWRK1 = LWRK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_R12BMAT')
      END IF
C
      CALL DZERO(WORK(KVAJKL),NVAJKL(ISYMA))
      CALL DZERO(WORK(KVIJKL),NTR12SQ(ISYRES))
C
      IF (.NOT.USEVABKL) THEN
C
        IOPT = 1
        CALL CC_R12MKVAMKL0(WORK(KVAJKL),NVAJKL(ISYMA),IOPT,LAMDHA,
     &                      ISYMA,WORK(KEND1),LWRK1)
        CALL CC_MOFCONR12(LAMDHA,ISYMA,IDUMMY,IDUMMY,IDUMMY,
     &                    IDUMMY,DUMMY,DUMMY,WORK(KVAJKL),IDUMMY,
     &                    .FALSE.,.TRUE.,.FALSE.,2,
     &                    TIM1,TIM2,TIM3,
     &                    DUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     &                    WORK(KEND1),LWRK1)
C
      ELSE
C
        KVABKL = KEND1
        KEND2  = KVABKL + NVABKL(1)
        LWRK2  = LWRK - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_R12BMAT')
        END IF
C
        LV = .TRUE.
        LVAJKL = .FALSE.
        LVIJKL = .FALSE.
        CALL CC_R12MKVTF(WORK(KVABKL),WORK(KVAJKL),DUMMY,
     &                   LAMDHA,ISYMA,
     &                   LV,LVIJKL,LVAJKL,CDUMMY,WORK(KEND2),LWRK2)
C
      END IF
C
      CALL CC_R12MKVIJKL(WORK(KVAJKL),ISYMA,LAMDHB,ISYMB,
     &                   WORK(KEND1),LWRK1,.TRUE.,ONE,WORK(KVIJKL))
      CALL CC_R12PKLIJ(WORK(KVIJKL),ISYRES,'T',WORK(KEND1),LWRK1)
C
C     --------------------------------------------------------------
C     add finite field contributions
C     --------------------------------------------------------------
      IF (NONHF) THEN
        !allocate memory:
        KT12AMP  = KEND0
        KXINTTRI = KT12AMP + MAX(NTR12SQ(ISYMA),NTR12SQ(ISYMB))
        KXINTSQ  = KXINTTRI + MAX(NR12R12P(1),NTR12SQ(ISYRES))
        KFIELDAO = KXINTSQ + NR12R12SQ(1)
        KSCR     = KFIELDAO + N2BST(1)
        KEND1    = KSCR + NTR12SQ(ISYRES)
        LWRK1    = LWRK - KEND1
        IF (LWRK1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_R12BMAT')
        END IF
C
        !initialize fields:
        CALL DZERO(WORK(KFIELDAO),N2BST(1))
C
        IF (ISYMOP.NE.1) CALL QUIT('ISYMOP .NE. 1 not implemented')
C
        !sum up fields:
        DO IFLD = 1, NFIELD
          IF ( NHFFIELD(IFLD) ) THEN
            CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1),LWRK1,
     *                   EFIELD(IFLD),1,LFIELD(IFLD))
          ELSE IF (.NOT. NHFFIELD(IFLD)) THEN
            CALL QUIT('CCR12 response can only handle '//
     &                'unrelaxed orbitals (w.r.t. the perturbation)')
          END IF
        END DO
C
        !read r12 overlap matrix
        LUNIT = -1
        CALL GPOPEN(LUNIT,FCCR12X,'OLD',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND(LUNIT)
 8888   READ(LUNIT) IAN
        READ(LUNIT) (WORK(KXINTTRI-1+I), I=1, NR12R12P(1))
        IF (IAN.NE.IANR12) GOTO 8888
        CALL GPCLOSE(LUNIT,'KEEP')
        IOPT = 2
        CALL CCR12UNPCK2(WORK(KXINTTRI),1,WORK(KXINTSQ),'N',IOPT)
C
        DO I = 1, 2
          IF (I.EQ.1) THEN
            ISYM1 = ISYMA
            ISYM2 = ISYMB
            LIST1 = LISTA
            IDLST1 = IDLSTA
          ELSE IF (I.EQ.2) THEN
            ISYM1 = ISYMB
            ISYM2 = ISYMA
            LIST1 = LISTB
            IDLST1 = IDLSTB
          END IF
          !read R12 response amplitudes
          CALL CC_R12GETCT(WORK(KT12AMP),ISYM1,2,KETSCL,.FALSE.,'N',
     &                 DUMMY,DUMMY,DUMMY,LIST1,IDLST1,WORK(KEND1),LWRK1)
          !calculate....
          IF (I.EQ.1) THEN
            CALL CC_R12XI2A(WORK(KSCR),ISYRES,WORK(KT12AMP),ISYM1,
     &                      WORK(KFIELDAO),1,LAMP0,LAMDHB,ISYM2,'N',
     &                      WORK(KEND1),LWRK1)
            CALL DCOPY(NTR12SQ(ISYRES),WORK(KSCR),1,WORK(KXINTTRI),1)
          ELSE IF (I.EQ.2) THEN
            CALL CC_R12XI2A(WORK(KSCR),ISYRES,WORK(KT12AMP),ISYM1,
     &                      WORK(KFIELDAO),1,LAMP0,LAMDHA,ISYM2,'N',
     &                      WORK(KEND1),LWRK1)
          END IF
        END DO
C
        CALL DAXPY(NTR12SQ(ISYRES),ONE,WORK(KXINTTRI),1,WORK(KSCR),1)
        CALL CC_R12XI2B(WORK(KVIJKL),'T',WORK(KXINTSQ),1,'N',
     &                  WORK(KSCR),ISYRES,'N',-ONE)
C
      END IF
C 
      !pack to triangular format:
      IOPT = 1
      CALL CCR12PCK2(WORK(KSCR),ISYRES,.FALSE.,WORK(KVIJKL),'T',IOPT)
      CALL CCLR_DIASCLR12(WORK(KSCR),BRASCL,ISYRES)
      !add to result:
      CALL DAXPY(NTR12AM(ISYRES),ONE,WORK(KSCR),1,THETAR12,1)
C
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) "Theta at end of CC_R12BMAT:"
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
         CALL CC_PRPR12(THETAR12,ISYRES,1,.TRUE.)
         WRITE(LUPRI,*) 'Leaving CC_R12BMAT'
      END IF
C
      CALL QEXIT('CC_R12BMAT')
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CC_R12BMAT                       *
*=====================================================================*
