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
c/* Deck CC_FMATRIX */
*=====================================================================*
      SUBROUTINE CC_FMATRIX(IFTRAN, NFTRAN, LISTL, LISTR, IOPTRES,
     &                      FILFMA, IFDOTS, FCONS, MXVEC, WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over F matrix transformations
*             (needed if the number of transformations exceeds the
*              limit MAXSIM defined on ccsdio.h )
*
*             for more detailed documentation see: CC_FMAT
*        
*     Written by Christof Haettig, November 1998, based on CC_BMATRIX.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "maxorb.h"
Cholesky
#include "ccdeco.h"
Cholesky
#include "ccsdio.h"
#include "ccnoddy.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*(*) LISTL, LISTR, FILFMA
      CHARACTER*8 FN3VI2, FNTOC
      CHARACTER*4 FNDKBC
      CHARACTER*6 FNDELD, FNCKJD, FN3VI
      INTEGER IOPTRES
      INTEGER NFTRAN, MXVEC, LWORK
      INTEGER IFTRAN(3,NFTRAN)
      INTEGER IFDOTS(MXVEC,NFTRAN)
      
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FCONS(MXVEC,NFTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION FCONS(MXVEC,NFTRAN)
COMMENT COMMENT
      DOUBLE PRECISION DDOT
COMMENT COMMENT
#endif 

      LOGICAL LADD
      INTEGER MAXFTRAN, NTRAN, ISTART, IBATCH, NBATCH
      INTEGER KIADRBFD, KIADRZ0, KIADRPQ0, KIADRPQR, KIADRPQMO   
      INTEGER KIADRBFI, KRHO1, KRHO2, KFIM, KMGDL, KRIM, KZDEN, KGIM 
      INTEGER KLAMPA, KLAMHA, KCTR2, KLAMP, KLAMH, KDENS, KFOCK
      INTEGER KEND, LEND
      INTEGER MXRAMP, MXLAMP, MXBTRAN, KB1TRAN, KB1DOTS, KB1CON,
     &        KB2TRAN, KB2DOTS, KB2CON, NB1TRAN, NB2TRAN

      CALL QENTER('CC_FMATRIX')

      IF (IOPTRES.EQ.1) 
     &      CALL QUIT('IOPTRES cannnot be used with CC_FMATRIX')

*---------------------------------------------------------------------*
* Main section:   CC_FMATNEW, driven by a loop over transformed vectors
* -------------
*   singles and doubles models:
*               calculation of the single and double excitation
*               parts of the transformed vectors and respective
*               dot products (IOPTRES=5).
*
*   triple models: 
*               calculation of the effective transformed vectors,
*               for IOPTRES=5 only singles and doubles contributions
*               and contributions from <L3|...|HF> are computed here
*               the remaining contributions <L2|[...,R3]|HF> are
*               calculated in CCSDT_FMAT_CONT (see below).
*---------------------------------------------------------------------*
*
Cholesky
C
      IF (CHOINT .AND. CC2) THEN

         CALL CC_CHOFMAT(IFTRAN, NFTRAN, LISTL, LISTR, IOPTRES,
     &                   FILFMA, IFDOTS, FCONS, MXVEC, WORK, LWORK)

         GOTO 9999
      END IF
C
Cholesky
*
      MAXFTRAN = MAXSIM

      NBATCH = (NFTRAN+MAXFTRAN-1)/MAXFTRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Batching over F matrix transformations:'
        WRITE (LUPRI,*) 'nb. of batches needed:', NBATCH
      END IF
  
      DO IBATCH = 1, NBATCH
        ISTART = (IBATCH-1) * MAXFTRAN + 1
        NTRAN  = MIN(NFTRAN-(ISTART-1),MAXFTRAN)

        KIADRBFD    = 1
        KIADRZ0     = KIADRBFD    + MXCORB_CC*NTRAN
        KIADRPQ0    = KIADRZ0     + MXCORB_CC*NTRAN
        KIADRPQR    = KIADRPQ0    + MXCORB_CC*NTRAN
        KIADRPQMO   = KIADRPQR    + MXCORB_CC*NTRAN
        KIADRBFI    = KIADRPQMO   + MXCORB_CC*NTRAN
        KRHO1       = KIADRBFI    + MXCORB_CC*NTRAN
        KRHO2       = KRHO1       + NTRAN
        KFIM        = KRHO2       + NTRAN
        KMGDL       = KFIM        + NTRAN
        KRIM        = KMGDL       + NTRAN
        KZDEN       = KRIM        + NTRAN
        KGIM        = KZDEN       + NTRAN
        KLAMPA      = KGIM        + NTRAN
        KLAMHA      = KLAMPA      + NTRAN
        KCTR2       = KLAMHA      + NTRAN
        KLAMP       = KCTR2       + NTRAN
        KLAMH       = KLAMP       + NTRAN
        KDENS       = KLAMH       + NTRAN
        KFOCK       = KDENS       + NTRAN
        KEND        = KFOCK       + NTRAN
        LEND        = LWORK       - KEND

        IF (LEND .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient work space in CC_FMATRIX.'
           WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
           WRITE (LUPRI,*) 'Need at least:',KEND, ' words.'
           CALL QUIT('Insufficient work space in CC_FMATRIX.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',ISTART
          WRITE (LUPRI,*) '# transf.:',NTRAN
        END IF

        CALL CC_FMATNEW(IFTRAN(1,ISTART), NTRAN,
     &                  LISTL, LISTR, IOPTRES, FILFMA, 
     &                  IFDOTS(1,ISTART), FCONS(1,ISTART), 
     &                  WORK(KIADRBFD), WORK(KIADRZ0), WORK(KIADRPQ0),
     &                  WORK(KIADRPQR),WORK(KIADRPQMO),WORK(KIADRBFI),
     &                  WORK(KRHO1),WORK(KRHO2),WORK(KFIM),WORK(KMGDL),
     &                  WORK(KRIM),WORK(KZDEN),WORK(KGIM),WORK(KLAMPA),
     &                  WORK(KLAMHA), WORK(KCTR2), WORK(KLAMP),
     &                  WORK(KLAMH), WORK(KDENS), WORK(KFOCK),
     &                  MXVEC, WORK(KEND), LEND )

      END DO
*---------------------------------------------------------------------*
* special triples section:
* ------------------------
*               compute the contributions <L2|[...,R3]|HF> to
*               F matrix contractions 
*---------------------------------------------------------------------*
      IF (IOPTRES.EQ.5 .AND. CCSDT .AND. CCSDT_F_ALTER) THEN
         MXRAMP  = MAX(NFTRAN,MXVEC)
         MXLAMP  = NFTRAN
         MXBTRAN = MXRAMP*MXRAMP

         KB1TRAN = 1
         KB1DOTS = KB1TRAN + MXBTRAN * 3
         KB1CON  = KB1DOTS + MXBTRAN * MXLAMP
         KEND    = KB1CON  + MXBTRAN * MXLAMP

         IF (LISTR(1:3).NE.FILFMA(1:3)) THEN
           KB2TRAN = KEND
           KB2DOTS = KB2TRAN + MXBTRAN * 3
           KB2CON  = KB2DOTS + MXBTRAN * MXLAMP
           KEND    = KB2CON  + MXBTRAN * MXLAMP
         ELSE
           KB2TRAN = -999999
           KB2DOTS = -999999
           KB2CON  = -999999
         END IF

         LEND = LWORK - KEND
         IF (LEND.LE.0) CALL QUIT('Out of memory in CC_FMATRIX.')
         
         CALL DZERO(WORK(KB1TRAN),MXBTRAN*3)
         CALL DZERO(WORK(KB1DOTS),MXBTRAN*MXLAMP)
         CALL DZERO(WORK(KB1CON), MXBTRAN*MXLAMP)
         IF (LISTR(1:3).NE.FILFMA(1:3)) THEN
           CALL DZERO(WORK(KB2TRAN),MXBTRAN*3)
           CALL DZERO(WORK(KB2DOTS),MXBTRAN*MXLAMP)
           CALL DZERO(WORK(KB2CON), MXBTRAN*MXLAMP)
         END IF

         LADD = .FALSE.
         CALL CCSDT_F2B_SETUP(IFTRAN,IFDOTS,FCONS,NFTRAN,MXVEC,LADD,
     &                        WORK(KB1TRAN),WORK(KB1DOTS),
     &                        WORK(KB1CON),NB1TRAN,MXLAMP,
     &                        WORK(KB2TRAN),WORK(KB2DOTS),
     &                        WORK(KB2CON),NB2TRAN,MXLAMP,
     &                        LISTL,LISTR,FILFMA,MXBTRAN,MXBTRAN)

         IF (NODDY_FMAT) THEN
             CALL CCSDT_FBC_NODDY(CCSDT_F_ALTER,
     &                            WORK(KB1TRAN),WORK(KB1DOTS),
     &                            WORK(KB1CON),NB1TRAN,MXLAMP,
     &                            LISTL,LISTR,FILFMA,WORK(KEND),LEND) 
         ELSE
             FNCKJD = 'CKJDEL'
             FNDKBC = 'DKBC'
             FNDELD = 'CKDELD'
             FNTOC  = 'CCSDT_OC'
             FN3VI  = 'CC3_VI'
             FN3VI2 = 'CC3_VI12'

             IF (LISTR(1:3).EQ.'R1 '.OR.LISTR(1:3).EQ.'RE '
     *           .OR.LISTR(1:3).EQ.'RC ') THEN
                CALL CCSDT_FBMAT(WORK(KB1TRAN),WORK(KB1DOTS),
     *                        WORK(KB1CON),NB1TRAN,MXLAMP,
     *                        LISTL,LISTR,FILFMA,FNCKJD,FNDKBC,
     *                        FNDELD,FNTOC,FN3VI,FN3VI2,
     *                        WORK(KEND),LEND) 
             ELSE IF (LISTR(1:3).EQ.'R2 '.OR.LISTR(1:3).EQ.'ER1') THEN
                CALL CC3_FBMATT3ZU(WORK(KB1TRAN),WORK(KB1DOTS),
     *                        WORK(KB1CON),NB1TRAN,MXLAMP,
     *                        LISTL,LISTR,FILFMA,
     *                        WORK(KEND),LEND) 
             ELSE
                CALL QUIT('Unknown LISTR in CC_FMAT.')
             END IF

         END IF

         IF (LISTR(1:3).NE.FILFMA(1:3)) THEN
           IF (NODDY_FMAT) THEN
             CALL CCSDT_FBC_NODDY(CCSDT_F_ALTER,
     &                            WORK(KB2TRAN),WORK(KB2DOTS),
     &                            WORK(KB2CON),NB2TRAN,MXLAMP,
     &                            LISTL,FILFMA,LISTR,WORK(KEND),LEND) 
           ELSE
             FNCKJD = 'CKJDEL'
             FNDKBC = 'DKBC'
             FNDELD = 'CKDELD'
             FNTOC  = 'CCSDT_OC'
             FN3VI  = 'CC3_VI'
             FN3VI2 = 'CC3_VI12'

             IF (FILFMA(1:3).EQ.'R1 '.OR.FILFMA(1:3).EQ.'RE '
     *           .OR.FILFMA(1:3).EQ.'RC ') THEN
                CALL CCSDT_FBMAT(WORK(KB2TRAN),WORK(KB2DOTS),
     *                        WORK(KB2CON),NB2TRAN,MXLAMP,
     *                        LISTL,FILFMA,LISTR,FNCKJD,FNDKBC,
     *                        FNDELD,FNTOC,FN3VI,FN3VI2,
     *                        WORK(KEND),LEND)
             ELSE IF (FILFMA(1:3).EQ.'R2 '.OR.FILFMA(1:3).EQ.'ER1') THEN
                CALL CC3_FBMATT3ZU(WORK(KB2TRAN),WORK(KB2DOTS),
     *                        WORK(KB2CON),NB2TRAN,MXLAMP,
     *                        LISTL,FILFMA,LISTR,
     *                        WORK(KEND),LEND)
             ELSE
                CALL QUIT('Unknown FILFMA in CC_FMAT.')
             END IF

           END IF
         END IF

         LADD = .TRUE.
         CALL CCSDT_F2B_SETUP(IFTRAN,IFDOTS,FCONS,NFTRAN,MXVEC,LADD,
     &                        WORK(KB1TRAN),WORK(KB1DOTS),
     &                        WORK(KB1CON),NB1TRAN,MXLAMP,
     &                        WORK(KB2TRAN),WORK(KB2DOTS),
     &                        WORK(KB2CON),NB2TRAN,MXLAMP,
     &                        LISTL,LISTR,FILFMA,MXBTRAN,MXBTRAN)

      END IF

      IF (LOCDBG) WRITE(LUPRI,*) 'Leave CC_FMATRIX...'

 9999 CONTINUE
      CALL QEXIT('CC_FMATRIX')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_FMATRIX                           *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
c/* Deck CC_FTRAN */
*=====================================================================*
      SUBROUTINE CC_FTRAN(LISTL, IDLSTL, LISTR, IDLSTR, WORK, LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: shortcut to F matrix transformation for a single
*             transformation (should be avoided and maybe completely
*             eliminated)
*        
*     Written by Christof Haettig, Januar 1999, based on CC_FMATRIX.
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

      INTEGER LUFMAT

      CHARACTER*(*) LISTL, LISTR
      CHARACTER*(8) FILFMA
      INTEGER IOPTRES, IDLSTL, IDLSTR, IDUM, NFTRAN, LWORK
      INTEGER IFTRAN(3,1)
      INTEGER KIADRBFD, KIADRZ0, KIADRPQ0, KIADRPQR, KIADRPQMO   
      INTEGER KIADRBFI, KRHO1, KRHO2, KFIM, KMGDL, KRIM, KZDEN, KGIM 
      INTEGER KLAMPA, KLAMHA, KCTR2, KLAMP, KLAMH, KDENS, KFOCK
      INTEGER KEND, LEND, ISYML, ISYMR, ISYRES, LENALL
      
      INTEGER ILSTSYM

      REAL*8  WORK(LWORK), RDUM

      CALL QENTER('CC_FTRAN')

      IF (CCSDT.OR.CC3) CALL
     &  QUIT('Triples not carried through in CC_FTRAN!')     


      NFTRAN      = 1
      IFTRAN(1,1) = IDLSTL
      IFTRAN(2,1) = IDLSTR
      IOPTRES     = 0
      FILFMA      = 'CC__FMAT'

      KIADRBFD    = 1
      KIADRZ0     = KIADRBFD    + MXCORB_CC*NFTRAN
      KIADRPQ0    = KIADRZ0     + MXCORB_CC*NFTRAN
      KIADRPQR    = KIADRPQ0    + MXCORB_CC*NFTRAN
      KIADRPQMO   = KIADRPQR    + MXCORB_CC*NFTRAN
      KIADRBFI    = KIADRPQMO   + MXCORB_CC*NFTRAN
      KRHO1       = KIADRBFI    + MXCORB_CC*NFTRAN
      KRHO2       = KRHO1       + NFTRAN
      KFIM        = KRHO2       + NFTRAN
      KMGDL       = KFIM        + NFTRAN
      KRIM        = KMGDL       + NFTRAN
      KZDEN       = KRIM        + NFTRAN
      KGIM        = KZDEN       + NFTRAN
      KLAMPA      = KGIM        + NFTRAN
      KLAMHA      = KLAMPA      + NFTRAN
      KCTR2       = KLAMHA      + NFTRAN
      KLAMP       = KCTR2       + NFTRAN
      KLAMH       = KLAMP       + NFTRAN
      KDENS       = KLAMH       + NFTRAN
      KFOCK       = KDENS       + NFTRAN
      KEND        = KFOCK       + NFTRAN
      LEND        = LWORK       - KEND

      IF (LEND .LT. 0) THEN
         WRITE (LUPRI,*) 'Insufficient work space in CC_FTRAN.'
         WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
         WRITE (LUPRI,*) 'Need at least:',KEND, ' words.'
         CALL QUIT('Insufficient work space in CC_FTRAN.')
      END IF

      CALL CC_FMATNEW( IFTRAN, NFTRAN, LISTL, LISTR, IOPTRES, FILFMA, 
     &                 IDUM, RDUM, 
     &                 WORK(KIADRBFD), WORK(KIADRZ0), WORK(KIADRPQ0),
     &                 WORK(KIADRPQR),WORK(KIADRPQMO),WORK(KIADRBFI),
     &                 WORK(KRHO1),WORK(KRHO2),WORK(KFIM),WORK(KMGDL),
     &                 WORK(KRIM),WORK(KZDEN),WORK(KGIM),WORK(KLAMPA),
     &                 WORK(KLAMHA), WORK(KCTR2), WORK(KLAMP),
     &                 WORK(KLAMH), WORK(KDENS), WORK(KFOCK),
     &                 0, WORK(KEND), LEND )

      ISYML = ILSTSYM(LISTL,IDLSTL)
      ISYMR = ILSTSYM(LISTR,IDLSTR)
      ISYRES = MULD2H(ISYML,ISYMR)

      LENALL = NT1AM(ISYRES) + NT2AM(ISYRES)
      IF (CCS) LENALL = NT1AM(ISYRES)
      IF (LENALL.GT.LWORK) THEN
         WRITE (LUPRI,*) 'Insufficient work space in CC_FTRAN.'
         WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
         WRITE (LUPRI,*) 'Need at least:',LENALL, ' words.'
         CALL QUIT('Insufficient work space in CC_FTRAN.')
      END IF

      LUFMAT = -1
      CALL WOPEN2(LUFMAT, FILFMA, 64, 0)
      CALL GETWA2(LUFMAT,FILFMA,WORK(1),1,LENALL) 
      CALL WCLOSE2(LUFMAT, FILFMA, 'DELETE')

      CALL QEXIT('CC_FTRAN')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_FTRAN                             *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
c/* Deck CC_FMATNEW */
*=====================================================================*
      SUBROUTINE CC_FMATNEW(IFTRAN, NFTRAN, LISTL, LISTR, IOPTRES,
     &                      FILFMA, IFDOTS, FCONS, 
     &                      IADRBFD, IADRZ0, IADRPQ0, IADRPQR,
     &                      IADRPQMO, IADRBFI, KRHO1, KRHO2, KFIM, 
     &                      KMGDL, KRIM, KZDEN, KGIM, KLAMPA,
     &                      KLAMHA, KCTR2, KLAMP, KLAMH,
     &                      KDENS, KFOCK, MXVEC, WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: AO-direct calculation of a linear transformation of two
*             CC amplitude vectors, L and R, with the CC F matrix
*             (or B matrix for L different from zeroth-order multipl.)
*
*                L.EQ.'L0'  --->      F * R       
*                L.NE.'L0'  --->  L * B * R       
*          
*             the linear transformations are calculated for a list
*             of L vectors and a list of R vectors: 
*
*                LISTL       -- type of L vectors
*                LISTR       -- type of R vectors
*                IFTRAN(1,*) -- indeces of L vectors
*                IFTRAN(2,*) -- indeces of R vectors
*                IFTRAN(3,*) -- indeces or addresses of result vectors
*                NFTRAN      -- number of requested transformations
*                FILFMA      -- file name / list type of result vectors
*                               or list type of vectors to be dotted on
*                IFDOTS      -- indeces of vectors to be dotted on
*                FCONS       -- contains the dot products on return
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  all result vectors are written to a direct
*                          access file, FILFMA is used as file name
*                          the start addresses of the vectors are
*                          returned in IFTRAN(3,*)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          IFTRAN(3,*). N.B.: if WORK is not large
*                          enough IOPTRES is automatically reset to 0!!
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, FILFMA is used
*                          as list type and IFTRAN(3,*) as list index
*                          NOTE that IFTRAN(3,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WARSP, FILFMA is used
*                          as list type and IFTRAN(3,*) as list index
*                          NOTE that IFTRAN(3,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by FILFMA and the indeces from IFDOTS
*                          the result of the dot products is returned
*                          in the FCONS array
*
*     Written by Christof Haettig, November 1998, based on CC_FMAT.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsections.h"
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
#include "ccslvinf.h"
#include "ccnoddy.h"
#include "second.h"
#include "dummy.h"
#include "r12int.h"
#include "qm3.h"
!#include "qmmm.h"

* local parameters:
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CC_FMAT> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL APPEND, NOAPPEND
      PARAMETER (APPEND = .TRUE., NOAPPEND = .FALSE.)

      INTEGER KDUM
      PARAMETER( KDUM = +99 999 999 ) ! dummy address for work space
      INTEGER ISYM0
      PARAMETER( ISYM0 = 1 ) ! symmetry of the reference state
      
      INTEGER LUBF,   LUBFD,  LUC,    LUD,    LUFK,   LURIM, LUBFI
      INTEGER LURHO,  LUFIM,  LUMGD,  LUXIM,  LUYIM,  LUBDZ0, LUO3
      INTEGER LUCBAR, LUDBAR, LUPQRR, LUPQR0, LUPQMO, LUFMAT, LUZDN
      INTEGER LUGIM,  LUMIM
      INTEGER LUCKJD, LUDKBC, LUTOC, LU3VI, LUDKBC3, LU3FOP, LU3FOP2
      INTEGER LU3VI2, LU3FOPX, LU3FOP2X

      CHARACTER*(8) BFFIL,  FNBFD,  CTFIL, DTFIL, FILFCK, FILRIM, FILO3
      CHARACTER*(8) RHOFIL,FIMFIL,FILMGD,FILXIM,FILYIM,FNBFI,FILZDN
      CHARACTER*(8) CBAFIL, DBAFIL, FNPQRR, FNPQR0, FNPQMO, FNBDZ0
      CHARACTER*(8) GIMFIL, FNTOC, FN3VI2, FILMIM
      CHARACTER*(7) FN3FOP2X
      CHARACTER*(6) FNCKJD, FN3VI, FN3FOP2, FNDELD, FN3FOPX
      CHARACTER*(5) FNDKBC3, FN3FOP
      CHARACTER*(4) FNDKBC
      CHARACTER RSPTYP(1)
      PARAMETER (BFFIL ='CCFM_BFI', FNBFD ='CCBFDENS',
     &           CBAFIL='CCFM_CBA', DBAFIL='CCFM_DBA',
     &           CTFIL ='CCFM_CIM', DTFIL ='CCFM_DIM', 
     &           FIMFIL='CCFM_FIM', RHOFIL='CCFM_RHO',
     &           FILRIM='CCFM_RIM', FILMGD='CCFM_MGD',
     &           FILXIM='CCFM_XIM', FILYIM='CCFM_YIM',
     &           FILMIM='CCFM_MIM',
     &           FILFCK='CCFM_FCK', FNPQRR='CCFMPQRR',
     &           FNPQR0='CCFMPQR0', FNPQMO='CCFMPQMO',
     &           FNBDZ0='CCFMBDZ0', FNBFI ='CCFM_BZI',
     &           FILO3 ='CCFM_O3I', FILZDN='CCFM_ZDN',
     &           GIMFIL='CCFM_GIM', FNCKJD='CKJDEL',
     &           FNDKBC='DKBC'    , FNTOC ='CCSDT_OC',
     &           FN3VI ='CC3_VI'  , FNDKBC3='DKBC3',
     &           FN3FOP='PTFOP'   , FN3FOP2='PTFOP2',
     &           FN3FOPX='PTFOPX' , FN3FOP2X='PTFOP2X',
     &           FNDELD='CKDELD'  , FN3VI2 = 'CC3_VI12' )


      CHARACTER*(*) LISTL, LISTR, FILFMA
      INTEGER IOPTRES, IOPTTCME, IDUM
      INTEGER NFTRAN, MXVEC, LWORK
      INTEGER IFTRAN(3,NFTRAN)
      INTEGER IFDOTS(MXVEC,NFTRAN)


#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FCONS(MXVEC,NFTRAN)
      REAL DUM, XNORM, FF, DUMMY
      REAL DTIME, CONVRT, TIMALL
      REAL TIMFCK, TIMBF, TIMC, TIMD, TIMIMR, TIMXYM, TIMIO
      REAL TIMPQ, TIMIML, TIMIM2, TIMPRE, TIMINT, TIMRDAO
      REAL TIMTRBT, TIMIM0, TIMTRN, TIMGZ, TIMFG
      REAL ZERO, ONE, TWO, HALF, FACT
C     REAL ECURRSAV
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION FCONS(MXVEC,NFTRAN)
      DOUBLE PRECISION XNORM, FF
      DOUBLE PRECISION DTIME, CONVRT, TIMALL
      DOUBLE PRECISION TIMFCK, TIMBF, TIMC, TIMD, TIMIMR, TIMXYM, TIMIO
      DOUBLE PRECISION TIMPQ, TIMIML, TIMIM2, TIMPRE, TIMINT, TIMRDAO
      DOUBLE PRECISION TIMTRBT, TIMIM0, TIMTRN, TIMGZ, TIMFG
      DOUBLE PRECISION ZERO, ONE, TWO, HALF, FACT
C     DOUBLE PRECISION ECURRSAV
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, HALF = 0.5d0)

      CHARACTER*(10) MODEL, MODELW, CDUMMY
      CHARACTER*(3) APROXR12
      LOGICAL LLSAME, CCSDR12, OSQSAV, OORSAV
      INTEGER INDEXA(MXCORB_CC)
      INTEGER INTMEDR(2,MAXSIM), NINTR, IINTR
      INTEGER INTMEDL(4,MAXSIM), NINTL, IINTL
      INTEGER INTMED2(4,MAXSIM), NINT2, IINT2
      INTEGER IRHGH(0:MAXSIM), I2HGH(0:MAXSIM), IOFFCD(0:MAXSIM+1)
      INTEGER IADRBFD(MXCORB_CC,NFTRAN), IADRZ0(MXCORB_CC,NFTRAN)
      INTEGER IADRPQ0(MXCORB_CC,NFTRAN), IADRPQR(MXCORB_CC,NFTRAN)
      INTEGER IADRPQMO(MXCORB_CC,NFTRAN),IADRBFI(MXCORB_CC,NFTRAN)
      INTEGER KRHO1(NFTRAN),KFIM(NFTRAN),KMGDL(NFTRAN),KRIM(NFTRAN)
      INTEGER KZDEN(NFTRAN),KGIM(NFTRAN)
      INTEGER KLAMPA(NFTRAN),KLAMHA(NFTRAN),KRHO2(NFTRAN),KCTR2(NFTRAN)
      INTEGER KLAMP(NFTRAN),KLAMH(NFTRAN),KDENS(NFTRAN),KFOCK(NFTRAN)
      INTEGER LENX,LENY,LENM,LENMGD,LENR,LENFIM,LENRHO,LENFK,LENBF,LEN
      INTEGER MT2BGD,MDISAO,MDSRHF,MSCRATCH,MEMAVAIL,NWORK,NSECMAR
      INTEGER NNWORK, ICDEL2, IBATCH, KEND2, JEND2, LWRK2
      INTEGER IOPTW, ITRAN, IOPT, ICORE, IF, IOPTRSP, ISTARTBFI, IVEC
      INTEGER ISYM, ICOUNT, ISYMAK, ISYBET, IBSRHF(8,8), NBSRHF(8)
      INTEGER ISYML, ISYMR, ISYRES, IDLSTL, IDLSTR, KEND1, JEND1
      INTEGER KFOCK0, KDENS0, KT1AMP0, KDSRHFA, ISYMXA, KFCKC0, KDNSC0
      INTEGER KLAMP0, KLAMH0, KTHETA1, KTHETA2, KZETA1, KT1AMPA
      INTEGER KCCFB1, KINDXB, KODCL1, KODCL2, KODBC1, KODBC2, NBATCH
      INTEGER KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, KRDPP2, KFREE
      INTEGER LFREE, KEND, LWRK, KENDSV, LWRKSV, KEND0, LWRK0, LWRK1
      INTEGER NTOSYM, NTOT, KRECNR, ISYMD1, ILLL, NUMDIS, IDEL2, IDEL
      INTEGER ISYDEL, KXINT, KEND3, K3OINT, KBFRHF, KDCRHF, KDSRHF
      INTEGER KFINT, KRIMA, KFOCKA, KTHETA0, KYINT, KYBARA, LENZDN
      INTEGER LWRK3, IADRTH, KLAMDPA, KLAMDHA, LENALL, KT2AMPA
      INTEGER ISYCTR, ISYAMP, KLIAJB, KA2IM, KXIDJL, KXJLID
      INTEGER KXIAJB, NDBLE, LENGIM, KGZETA, KFOCK0AO, KXBARA
      INTEGER KTHETA1EFF, KTHETA2EFF, IOPTWE
      INTEGER KTHETAR12,LENMOD,IOPTWR12,KSCR1,KSCR2,ICON,IAMP
      INTEGER KMINT

* external functions:
      INTEGER ICCSET1
      INTEGER ICCSET2
      INTEGER ILSTSYM

#if defined (SYS_CRAY) 
      REAL DDOT
#else 
      DOUBLE PRECISION DDOT
#endif 
  
      CALL QENTER('CC_FMATNEW')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_FMAT')
        IF (DIRECT) WRITE(LUPRI,'(/1X,A)') 'AO direct transformation'
        WRITE (LUPRI,*) MSGDBG, 'LISTL : ',LISTL(1:3)
        WRITE (LUPRI,*) MSGDBG, 'LISTR : ',LISTR(1:3)
        WRITE (LUPRI,*) MSGDBG, 'FILFMA: ',FILFMA(1:3)
        WRITE (LUPRI,*) MSGDBG, 'NFTRAN: ',NFTRAN
        WRITE (LUPRI,*) MSGDBG, 'IOPTRES:',IOPTRES
        CALL FLSHFO(LUPRI)
      END IF
      
      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
        WRITE (LUPRI,'(/1x,2a)')'CC_FMAT called for a Coupled Cluster ',
     &          'method not implemented in CC_FMAT...'
        CALL QUIT('Unknown CC method in CC_FMAT.')
      END IF

      IF (.NOT. DUMPCD) THEN
        WRITE (LUPRI,*) 'DUMPCD = ',DUMPCD
        WRITE (LUPRI,*) 'CC_FMAT requires DUMPCD=.TRUE.'
        CALL QUIT('DUMPCD=.FALSE. , CC_FMAT requires DUMPCD=.TRUE.')
      END IF

      IF (NFTRAN .GT. MAXSIM) THEN
        WRITE (LUPRI,*) 
     *      'Error in CC_FMAT: NFTRAN is larger than MAXSIM.'
        CALL QUIT('Error in CC_FMAT: NFTRAN is larger than MAXSIM.')
      END IF

      IF (IPRINT.GT.2) THEN
 
         WRITE (LUPRI,'(//1X,A1,50("="),A1)') '+','+'

         WRITE (LUPRI,'(1x,A52)')
     &         '|        F MATRIX TRANSFORMATION SECTION           |'

         IF (IOPTRES.EQ.0) THEN
            WRITE (LUPRI,'(1X,A52)')
     &         '|  (results are written to a direct access file)   |'
         ELSE IF (IOPTRES.EQ.1) THEN
            WRITE (LUPRI,'(1X,A52)')
     &         '|        (results are returned in memory)          |'
         ELSE IF (IOPTRES.EQ.3) THEN
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
      TIMIM0  = ZERO
      TIMIML  = ZERO
      TIMIMR  = ZERO
      TIMIM2  = ZERO
      TIMPRE  = ZERO
      TIMXYM  = ZERO
      TIMPQ   = ZERO
      TIMFCK  = ZERO
      TIMBF   = ZERO
      TIMC    = ZERO
      TIMD    = ZERO
      TIMFG   = ZERO
      TIMGZ   = ZERO
      TIMINT  = ZERO
      TIMRDAO = ZERO
      TIMTRBT = ZERO
      TIMIO   = ZERO

* set option and model to write vectors to file:
      IOPTW    = 0
      IOPTWR12 = 0
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
         CALL QUIT('Unknown coupled cluster model in CC_FMAT.')
      END IF
      IF (CCR12) THEN
        APROXR12 = '   '
        CALL CCSD_MODEL(MODELW,LENMOD,10,MODELW,10,APROXR12)
        IOPTWR12 = 32
      END IF

      !set CCSDR12
      CCSDR12 = CCR12 .AND. .NOT.(CC2.OR.CCS)

* check return option for the result vectors:
      LUFMAT = -1
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
         CALL WOPEN2(LUFMAT, FILFMA, 64, 0)
      ELSE IF (IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4) THEN
         CONTINUE
      ELSE IF (IOPTRES .EQ. 5) THEN
         IF (MXVEC*NFTRAN.NE.0) CALL DZERO(FCONS,MXVEC*NFTRAN)
      ELSE
         CALL QUIT('Illegal value of IOPTRES in CC_FMAT.')
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

* array for intermediates which depend on the left vectors:
      NINTL = 0
      DO ITRAN = 1, NFTRAN
        I=ICCSET2(INTMEDL,LISTL,IFTRAN(1,ITRAN),
     &                    'R0 ',0,              NINTL,MAXSIM,APPEND)
      END DO 

* array for intermediates which depend on the right vectors:
      NINTR = 0
      DO ITRAN = 1, NFTRAN
        I=ICCSET1(INTMEDR,LISTR,IFTRAN(2,ITRAN),NINTR,MAXSIM,APPEND)
      END DO 

* array for intermediates that depend on left and right vectors: 
      NINT2 = 0
      DO ITRAN = 1, NFTRAN
        I=ICCSET2(INTMED2,LISTL,IFTRAN(1,ITRAN), 
     &                    LISTR,IFTRAN(2,ITRAN),NINT2,MAXSIM,APPEND)
      END DO 


      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/A)')
     &        'List of response vector for left intermediates:'
        WRITE (LUPRI,'((/5X,2I5))') ((INTMEDL(I,J),I=1,2),J=1,NINTL)
        WRITE (LUPRI,'(/A)')
     &       'List of response vector for right intermediates:'
        WRITE (LUPRI,'((/5X,2I5))') ((INTMEDR(I,J),I=1,2),J=1,NINTR)
        WRITE (LUPRI,'(/A)')
     &       'List of vector pairs for 2. order intermediates:'
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

      MSCRATCH = 5*MT2BGD + MDISAO + MDSRHF + 10*N2BASX  
      IF (CC2)               MSCRATCH = MSCRATCH+MDISAO+2*MDSRHF
      IF (.NOT.(CCS.OR.CC2)) MSCRATCH = MSCRATCH+MDISAO+4*MDSRHF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'scratch space estimate MSCRATCH:',
     &        MSCRATCH
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
          CALL QUIT('Unknown CC model in CC_FMAT.')
      END IF

      IRHGH(0) = 0
      I2HGH(0) = 0

* intermediates that dependent on right response vectors:
* (see routine CCBPRE1 for details)
      DO IINTR = 1, NINTR
        IDLSTR = INTMEDR(1,IINTR)
        ISYM   = ILSTSYM(LISTR,IDLSTR)
 
        NNWORK = 2*NGLMDT(ISYM) + 2*N2BST(ISYM)
        IF (CCSD .OR. CCSDT) THEN
          NNWORK = 2*NGLMDT(ISYM)+NT2AOIJ(ISYM)+NEMAT1(ISYM)
        END IF

        IF( (NWORK+NNWORK+NSECMAR).GT.MEMAVAIL ) THEN
          IRHGH(NBATCH) = IINTR - 1
          I2HGH(NBATCH) = 0
       
          NBATCH = NBATCH + 1
          NWORK  = 0
        END IF
        NWORK = NWORK + NNWORK
        IF (NWORK .GT. LWORK) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CC_FMAT. (01)'
          WRITE (LUPRI,*) 'Need at least:',NNWORK, ' words.'
          CALL FLSHFO(LUPRI)
          CALL QUIT('Insufficient work space in CC_FMAT. (01)')
        END IF
      END DO

* intermediates that dependent on left and right vectors:
* (see routine CCFPRE2 for details)
      DO IINT2 = 1, NINT2
        IDLSTL = INTMED2(1,IINT2)
        IDLSTR = INTMED2(3,IINT2)
        ISYML  = ILSTSYM(LISTL,IDLSTL)
        ISYMR  = ILSTSYM(LISTR,IDLSTR)
        ISYRES = MULD2H(ISYML,ISYMR)
 
        IF (CCS) THEN
          NNWORK = 2*N2BST(ISYRES)
        ELSE IF (CC2) THEN
          NNWORK = 2*N2BST(ISYRES) + 2*NGLMDT(ISYMR) +
     &               NT2SQ(ISYML)  + NT1AM(ISYRES)
        ELSE IF (CCSD .OR. CCSDT) THEN
          NNWORK = 2*NGLMDT(ISYMR) + 
     &               NT1AO(ISYRES) + NT1AM(ISYRES) + NT2AOIJ(ISYRES)
        ELSE
          CALL QUIT('Unknown CC model in CC_FMAT.')
        END IF

        IF( (NWORK+NNWORK+NSECMAR).GT.MEMAVAIL ) THEN
          IRHGH(NBATCH) = NINTR
          I2HGH(NBATCH) = IINT2 - 1
       
          NBATCH = NBATCH + 1
          NWORK  = 0
        END IF
        NWORK = NWORK + NNWORK
        IF (NWORK .GT. LWORK) THEN
          WRITE (LUPRI,*) 'Insufficient work space in CC_FMAT. (02)'
          WRITE (LUPRI,*) 'Need at least:',NNWORK,' words.'
          CALL FLSHFO(LUPRI)
          CALL QUIT('Insufficient work space in CC_FMAT. (02)')
        END IF
      END DO
 
      IRHGH(NBATCH) = NINTR
      I2HGH(NBATCH) = NINT2

      IF   (LOCDBG .AND. (NBATCH.EQ.1)) THEN
        WRITE (LUPRI,*) MSGDBG,'one batch only... will be done in core.'
        WRITE (LUPRI,*) MSGDBG, 'memory for intermediates: ', NWORK
        WRITE (LUPRI,*) MSGDBG, 'remaining scratch space: ', LWORK-NWORK
        CALL FLSHFO(LUPRI)
      ELSE IF (LOCDBG .AND. (NBATCH.GT.1)) THEN
        WRITE (LUPRI,*) MSGDBG, 
     &        'more than one batch... choose I/O algorithm.'
        WRITE (LUPRI,*) MSGDBG, 'max. memory for intermediates: ',
     &       MEMAVAIL
        WRITE (LUPRI,*) MSGDBG, 'number of batches: ',NBATCH
        CALL FLSHFO(LUPRI)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* read zeroth-order singles amplitudes, allocate space for Fock matrix,
* and prepare zeroth-order lambda matrices and density: 
*---------------------------------------------------------------------*
      DTIME = SECOND()

      KFOCK0AO = 1
      KFOCK0   = KFOCK0AO + N2BAST

      KDENS0   = KFOCK0   + N2BAST
      KT1AMP0  = KDENS0   + N2BAST
      KLAMP0   = KT1AMP0  + NT1AMX
      KLAMH0   = KLAMP0   + NLAMDT
      KEND0    = KLAMH0   + NLAMDT

      IF (FROIMP.OR.FROEXP) THEN
        KFCKC0   = KEND0
        KDNSC0   = KFCKC0 + N2BAST
        KEND0    = KDNSC0 + N2BAST
      END IF

      LWRK0    = LWORK - KEND0

      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_FMAT. (0)')
      END IF

* read zeroth order amplitudes:
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

* get unperturbed Lambda matrices:
      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND0),LWRK0)

* calculate the density matrix:
      ICORE = 1 ! include core contribution
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
 
* Solvent contribution to one-electron integrals. 
* T^g contribution to transformation. SLV98,OC, CCMM JK, NYQMMM KS
      IF (CCSLV) THEN
         IF (.NOT.CCMM) CALL CCSL_RHSTG(WORK(KFOCK0),WORK(KEND0),LWRK0)
         IF (CCMM) THEN
            IF (.NOT. NYQMMM) THEN
               CALL CCMM_RHSTG(WORK(KFOCK0),WORK(KEND0),LWRK0)
            ELSE IF (NYQMMM) THEN
              IF (HFFLD) THEN
                CALL CCMM_ADDGHF(WORK(KFOCK0),WORK(KEND0),LWRK0)
              ELSE
                CALL CCMM_ADDG(WORK(KFOCK0),WORK(KEND0),LWRK0)
              END IF
            END IF
         END IF
      ENDIF
C
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG, 'norm of T1AMP0:',
     &        DDOT(NT1AMX,WORK(KT1AMP0),1,WORK(KT1AMP0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of XLAMP0:',
     &        DDOT(NLAMDT,WORK(KLAMP0),1,WORK(KLAMP0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of XLAMH0:',
     &        DDOT(NLAMDT,WORK(KLAMH0),1,WORK(KLAMH0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of DENS0:',
     &        DDOT(N2BAST,WORK(KDENS0),1,WORK(KDENS0),1)
        WRITE (LUPRI,*) MSGDBG, 'norm of FOCK0:',
     &        DDOT(N2BAST,WORK(KFOCK0),1,WORK(KFOCK0),1)
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME

*---------------------------------------------------------------------*
* precalculate some intermediates which depend only on the left 
* vectors: X, Y, M, Chi, Yps, P, Q, and the BFZeta density
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IOPTRSP = 1
      CALL CCFPREINT(INTMEDL,NINTL,WORK(KLAMP0),WORK(KLAMH0),
     &             CDUMMY,LENX,CDUMMY,LENY,CDUMMY,LENM,
     &             CDUMMY,LENMGD,CDUMMY,LENZDN,
     &             FNPQR0,IADRPQ0,FNPQMO,IADRPQMO,FNBDZ0,IADRZ0,
     &             TIMXYM,TIMBF,TIMIO,TIMPQ,IOPTRSP,WORK(KEND0),LWRK0)

      TIMIML = TIMIML + SECOND() - DTIME

*---------------------------------------------------------------------*
* precalculate some intermediates which depend only on the right
* vectors: CBAR, DBAR, and the BF density
*---------------------------------------------------------------------*
      DTIME = SECOND()

      CALL CCBPRE1INT(INTMEDR,NINTR,IOFFCD,IADRBFD,
     &                CBAFIL,DBAFIL,FNBFD,
     &                WORK(KLAMP0),WORK(KLAMH0),
     &                TIMIO,TIMC,TIMD,TIMBF,WORK(KEND0),LWRK0)

      TIMIMR = TIMIMR + SECOND() - DTIME

*---------------------------------------------------------------------*
* precalculate some intermediates which depend on left and right
* vectors: response X, Y, M, Chi, Yps, P, Q, and BFZeta density
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IOPTRSP = 2
      CALL CCFPREINT(INTMED2,NINT2,WORK(KLAMP0),WORK(KLAMH0),
     &             FILXIM,LENX,FILYIM,LENY,FILMIM,LENM,
     &             FILMGD,LENMGD,FILZDN,LENZDN,
     &             FNPQRR,IADRPQR,CDUMMY,IDUMMY,CDUMMY,IDUMMY,
     &             TIMXYM,TIMBF,TIMIO,TIMPQ,IOPTRSP,WORK(KEND0),LWRK0)

      TIMIM2 = TIMIM2 + SECOND() - DTIME

*---------------------------------------------------------------------*
* open files for local intermediates generated in the loop over the
* AO integrals and which need to be initialized here:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      CALL CCFOPEN(LUBF,  LUFK,   LURIM,  LUFIM,  LUGIM,  LURHO,
     &             BFFIL, FILFCK, FILRIM, FIMFIL, GIMFIL, RHOFIL,
     &             LENBF, LENFK,  LENR,   LENFIM, LENGIM, LENRHO,
     &             NINTR, NINT2,  WORK(KEND0),    LWRK0           )

* open some other files needed in the loop over integrals, but need
* not to be initialized:
      LUBFD  = -1
      LUZDN  = -1
      LUC    = -1
      LUD    = -1
      LUMGD  = -1
      LUBFI  = -1
      LUBFD  = -1
      LUPQR0 = -1
      LUPQRR = -1
      LUBDZ0 = -1
      LUO3   = -1
      IF (CCS.OR.CC2) THEN
         CALL WOPEN2(LUZDN, FILZDN, 64,0)
      ELSE IF (.NOT.(CCS.OR.CC2)) THEN
         CALL WOPEN2(LUC,   CTFIL, 64,0)
         CALL WOPEN2(LUD,   DTFIL, 64,0)
         CALL WOPEN2(LUMGD, FILMGD,64,0)
         CALL WOPEN2(LUBFI, FNBFI, 64,0)
         CALL WOPEN2(LUBFD, FNBFD, 64,0)
         CALL WOPEN2(LUPQR0,FNPQR0,64,0)
         CALL WOPEN2(LUPQRR,FNPQRR,64,0)
         CALL WOPEN2(LUBDZ0,FNBDZ0,64,0)
         CALL WOPEN2(LUO3 , FILO3, 64,0)
      END IF

* initialize offsets C and D intermediates:
      ICDEL2    = 0
      
* initialize start address for BFI intermediates:
      ISTARTBFI = 1

      TIMPRE = TIMPRE + SECOND() - DTIME
*---------------------------------------------------------------------*
* if all vectors and intermediates fit into the memory, read all
* response vectors before the loop over AO integral shells:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (NBATCH .EQ. 1) THEN

            CALL CCBPRE1(INTMEDR, 1, NINTR,
     &                   KRHO2, KLAMP, KLAMH, KDENS, KFOCK, KRIM,
     &                   LUBF,BFFIL,LENBF,LUFK,FILFCK,LENFK,
     &                   LURIM,FILRIM,LENR,
     &                   WORK(KLAMP0), WORK(KLAMH0),
     &                   WORK, LWORK, KEND0, JEND1 )
            KEND1 = JEND1

            CALL CCFPRE2(INTMED2,1,NINT2,
     &                   KFIM, LUFIM,FIMFIL,LENFIM,
     &                   KGIM, LUGIM,GIMFIL,LENGIM,
     &                   KRHO1,LURHO,RHOFIL,LENRHO,
     &                   KMGDL,LUMGD,FILMGD,LENMGD,
     &                   KZDEN,LUZDN,FILZDN,LENZDN,
     &                   KCTR2,KLAMPA,KLAMHA,
     &                   WORK(KLAMP0),WORK(KLAMH0),
     &                   WORK, LWORK, KEND1, JEND1  )
            KEND1 = JEND1

            IF (LOCDBG) THEN
             WRITE (LUPRI,*) MSGDBG,
     &              'allocated work space for intermediates:'
             WRITE (LUPRI,*) MSGDBG,'KRHO2 :',(KRHO2(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KLAMP :',(KLAMP(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KLAMH :',(KLAMH(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KDENS :',(KDENS(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KFOCK :',(KFOCK(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KRIM  :',(KRIM(I),I=1,NINTR)
             WRITE (LUPRI,*) MSGDBG,'KRHO1 :',(KRHO1(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KZDEN :',(KZDEN(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KFIM  :',(KFIM(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KGIM  :',(KGIM(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KMGDL :',(KMGDL(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KLAMPA:',(KLAMPA(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KLAMHA:',(KLAMHA(I),I=1,NINT2)
             WRITE (LUPRI,*) MSGDBG,'KEND1:',KEND1
             CALL FLSHFO(LUPRI)
            END IF

      ELSE
        KEND1 = KEND0
      END IF

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_FMAT. (1)')
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
* start loop over AO integral shells:
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
 
            IF (HERDIR) THEN
               CALL HERDI2(WORK(KEND),LWRK,INDEXA,ILLL,NUMDIS,
     &                     IPRINT)
            ELSE
               CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                     WORK(KODCL1),WORK(KODCL2),
     *                     WORK(KODBC1),WORK(KODBC2),
     *                     WORK(KRDBC1),WORK(KRDBC2),
     *                     WORK(KODPP1),WORK(KODPP2),
     *                     WORK(KRDPP1),WORK(KRDPP2),
     *                     WORK(KCCFB1),WORK(KINDXB),
     *                     WORK(KEND),LWRK,IPRERI)
            END IF

            KRECNR = KEND
            KEND   = KRECNR + (NBUFX(0) - 1)/IRAT + 1
            LWRK   = LWORK - KEND
  
            IF (LWRK .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_FMAT. (1a)')
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
               WRITE (LUPRI,*) MSGDBG, IBATCH,
     &                         '-th. batch out of ',NBATCH
               WRITE (LUPRI,*) MSGDBG, 'IR:',IRHGH(IBATCH-1)+1,' -- ',
     &                                  IRHGH(IBATCH)
               WRITE (LUPRI,*) MSGDBG, 'I2:',I2HGH(IBATCH-1)+1,' -- ',
     &                                  I2HGH(IBATCH)
             END IF
 
             IF (NBATCH.GT.1) THEN
 
               DTIME = SECOND()

               CALL CCBPRE1(INTMEDR,IRHGH(IBATCH-1)+1,IRHGH(IBATCH),
     &                      KRHO2, KLAMP, KLAMH, KDENS, KFOCK, KRIM,
     &                      LUBF,BFFIL,LENBF,LUFK,FILFCK,LENFK,
     &                      LURIM,FILRIM,LENR,
     &                      WORK(KLAMP0), WORK(KLAMH0),
     &                      WORK, LWORK, KEND2, JEND2 )
               KEND2 = JEND2
 
               CALL CCFPRE2(INTMED2,I2HGH(IBATCH-1)+1,I2HGH(IBATCH),
     &                      KFIM, LUFIM,FIMFIL,LENFIM,
     &                      KGIM, LUGIM,GIMFIL,LENGIM,
     &                      KRHO1,LURHO,RHOFIL,LENRHO,
     &                      KMGDL,LUMGD,FILMGD,LENMGD,
     &                      KZDEN,LUZDN,FILZDN,LENZDN,
     &                      KCTR2,KLAMPA,KLAMHA,
     &                      WORK(KLAMP0),WORK(KLAMH0),
     &                      WORK, LWORK, KEND2, JEND2 )
               KEND2 = JEND2

               TIMPRE = TIMPRE + SECOND() - DTIME
 
             END IF

             LWRK2 = LWORK - KEND2
             IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_FMAT. (2)')
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
cch
c              write(lupri,*) 'NOAUXB,IDEL2,IDEL,ISYDEL:',
c    &                         NOAUXB,IDEL2,IDEL,ISYDEL
cch
 
*           read AO integral distribution and calculate integrals with
*           one index transformed to occupied MO (particle):
            
            KXINT  = KEND2
            KEND3  = KXINT  + NDISAO(ISYDEL)

            IF (CC2) THEN
               KDSRHF  = KEND3
               KDSRHFA = KDSRHF  + NDSRHF(ISYDEL)
               KEND3   = KDSRHFA + MDSRHF
            ELSE IF (CCSD .OR. CCSDT) THEN
               K3OINT  = KEND3
               KBFRHF  = K3OINT + NMAIJK(ISYDEL)
               KDCRHF  = KBFRHF + NBSRHF(ISYDEL)
               KDSRHF  = KDCRHF + NBSRHF(ISYDEL)
               KEND3   = KDSRHF + NDSRHF(ISYDEL)
            END IF

            LWRK3  = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_FMAT. (3)')
            END IF
 
            DTIME = SECOND()
            CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     &                  WORK(KRECNR),DIRECT)
            TIMRDAO = TIMRDAO + SECOND() - DTIME

C
            IF (LOCDBG) THEN
              XNORM = DDOT(NDISAO(ISYDEL),WORK(KXINT),1,WORK(KXINT),1)
              WRITE (LUPRI,*) 'IDEL,XNORM:',IDEL,XNORM
            END IF
C

            IF (CCSD .OR. CCSDT) THEN
C              -----------------------------------------------------
C              some integral transformations and presortings for
C              BF, C, and D intermediate and 21G term:
C                   DSRHF  --  standard (**|kdel) integrals
C                   BFRHF  --  (**|kdel) presorted for B & D interm.
C                   DCRHF  --  (**|kdel) presorted for C & D interm.
C                   3OINT  --  (ij|kdel) integrals for 21G term
C              -----------------------------------------------------
               DTIME = SECOND()

               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMP0),ISYM0,
     &                     WORK(KEND3),LWRK3,ISYDEL)

               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMH0),
     &                       ISYM0,WORK(KLAMP0),WORK(KEND3),
     &                       LWRK3,IDEL,ISYDEL,LUO3,FILO3)

               CALL CC_BFBSORT(WORK(KDSRHF),WORK(KBFRHF),ISYDEL)
               
               CALL CCB_CDSORT(WORK(KXINT),ISYDEL,WORK(KDCRHF),
     &                         WORK(KLAMP0),ISYM0,WORK(KEND3),LWRK3)

               TIMTRBT = TIMTRBT + SECOND() - DTIME
               TIMIM0  = TIMIM0  + SECOND() - DTIME

               KEND3  = KDSRHF
               LWRK3  = LWORK  - KEND3

               KDSRHF = KDUM

            END IF


C           ----------------------------------------------------------
C           calculate zeroth-order Fock matrix (Fhat), as well as the
C           pure core contribution to the zeroth-order Fock matrix
C           ----------------------------------------------------------
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
 

C           ----------------------------------------------------------
C           calculate intermediates that depend only on right vectors:
C           ----------------------------------------------------------
            DO IINTR = IRHGH(IBATCH-1)+1, IRHGH(IBATCH)
              IDLSTR = INTMEDR(1,IINTR)
              ISYMR  = ILSTSYM(LISTR,IDLSTR)

*             calculate addresses for C & D intermediates:
              IT2DLR(IDEL,IINTR) = ICDEL2
              ICDEL2 = ICDEL2 + NT2BCD(MULD2H(ISYDEL,ISYMR))

              DTIME = SECOND()
              CALL CCBINT1(WORK(KXINT), WORK(KBFRHF), WORK(KDCRHF),
     &                     IDEL, ISYDEL,  WORK(KRHO2(IINTR)),  
     &                     WORK(KLAMP0),  WORK(KLAMH0), 
     &                     WORK(KLAMP(IINTR)),  WORK(KLAMH(IINTR)),
     &                     ISYMR, IINTR,
     &                     WORK(KDENS(IINTR)),  WORK(KFOCK(IINTR)),  
     &                     WORK(KRIM(IINTR)),
     &                     LUC, CTFIL, LUD, DTFIL, 
     &                     LUBFD, FNBFD, IADRBFD(1,IINTR), 
     &                     WORK(KEND3),   LWRK3,          
     &                     TIMFCK, TIMBF, TIMC, TIMD  )
              TIMIMR = TIMIMR + SECOND() - DTIME

            END DO
 
 
C           ----------------------------------------------------------
C           calculate intermediates that depend on both the left and
C           the right vectors:
C           ----------------------------------------------------------
            DO IINT2 = I2HGH(IBATCH-1)+1, I2HGH(IBATCH)
              IDLSTL = INTMED2(1,IINT2)
              IDLSTR = INTMED2(3,IINT2)
              ISYML  = ILSTSYM(LISTL,IDLSTL)
              ISYMR  = ILSTSYM(LISTR,IDLSTR)
              ISYRES = MULD2H(ISYML,ISYMR)

              IINTL  = ICCSET2(INTMEDL,LISTL,IDLSTL,
     &                         'R0 ',0,NINTL,MAXSIM,NOAPPEND)

              DTIME = SECOND()
              CALL CCFINT2(IDEL,ISYDEL,ISYML,WORK(KXINT),
     &                   WORK(KBFRHF),WORK(K3OINT),
     &                   WORK(KMGDL(IINT2)),WORK(KZDEN(IINT2)),
     &                   WORK(KFIM(IINT2)), WORK(KGIM(IINT2)),
     &                   WORK(KRHO1(IINT2)),
     &                   LUBFI,FNBFI,IADRBFI(1,IINT2),ISTARTBFI,
     &                   LUPQR0,FNPQR0,IADRPQ0(1,IINTL),
     &                   LUPQRR,FNPQRR,IADRPQR(1,IINT2),
     &                   LUBDZ0,FNBDZ0,IADRZ0(1,IINTL),
     &                   WORK(KLAMPA(IINT2)),WORK(KLAMHA(IINT2)),ISYMR,
     &                   WORK(KLAMP0), WORK(KLAMH0), 
     &                   WORK(KEND3),  LWRK3,
     &                   TIMFCK, TIMBF, TIMGZ, TIMIO, TIMFG      )
              TIMIM2 = TIMIM2 + SECOND() - DTIME

            END DO
 
C           ----------------------------------------------------------
C           for CC2 calculate the 21CD term contributions.
C           ----------------------------------------------------------
            IF (CC2) THEN
C            ........................................................
C            DSRHF  --  (**|kdel) integrals, transformed with XLAMH0:
C            ........................................................
             DTIME = SECOND()
             CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMH0),ISYM0,
     &                   WORK(KEND3),LWRK3,ISYDEL)
             TIMTRBT = TIMTRBT + SECOND() - DTIME
             TIMIM0  = TIMIM0  + SECOND() - DTIME

             DO IINT2 = I2HGH(IBATCH-1)+1, I2HGH(IBATCH)
              IDLSTL = INTMED2(1,IINT2)
              IDLSTR = INTMED2(3,IINT2)
              ISYML  = ILSTSYM(LISTL,IDLSTL)
              ISYMR  = ILSTSYM(LISTR,IDLSTR)
              ISYRES = MULD2H(ISYML,ISYMR)

              IINTL  = ICCSET2(INTMEDL,LISTL,IDLSTL,
     &                         'R0 ',0,NINTL,MAXSIM,NOAPPEND)

C             .............................................
C             calculate the first term depending on DSRHF0:
C             .............................................
              LLSAME = .FALSE.
              FACT   = ONE
              DTIME = SECOND()
              CALL CCG_21CD( WORK(KRHO1(IINT2)),ISYRES,
     &                   WORK(KCTR2(IINT2)),ISYML,
     &                   WORK(KLAMH0), WORK(KLAMP0), 
     &                   WORK(KLAMHA(IINT2)),WORK(KLAMPA(IINT2)),ISYMR,
     &                   WORK(KLAMH0), WORK(KLAMP0), ISYM0,
     &                   WORK(KDSRHF), ISYDEL, IDEL, ISYDEL,
     &                   LLSAME, FACT, WORK(KEND3),  LWRK3 )

C             ..............................................
C             calculate the second term depending on DSRHFA:
C             ..............................................
              DTIME = SECOND()
              CALL CCTRBT(WORK(KXINT),WORK(KDSRHFA),
     &                    WORK(KLAMHA(IINT2)),ISYMR,
     &                    WORK(KEND3),LWRK3,ISYDEL)
              TIMTRBT = TIMTRBT + SECOND() - DTIME


              LLSAME = .TRUE.
              FACT   = HALF
              ISYMXA = MULD2H(ISYDEL,ISYMR)
              DTIME = SECOND()
              CALL CCG_21CD( WORK(KRHO1(IINT2)),ISYRES,
     &                   WORK(KCTR2(IINT2)),ISYML,
     &                   WORK(KLAMH0), WORK(KLAMP0), 
     &                   WORK(KLAMH0), WORK(KLAMP0), ISYM0,
     &                   WORK(KLAMH0), WORK(KLAMP0), ISYM0,
     &                   WORK(KDSRHFA), ISYMXA, IDEL, ISYDEL,
     &                   LLSAME, FACT, WORK(KEND3),  LWRK3 )
              TIMIM2 = TIMIM2 + SECOND() - DTIME

             END DO
            END IF

          END DO ! IDEL2
*---------------------------------------------------------------------*
*         end of the loop over integral distributions:
*         if batched I/O algorithm used, save result on disc:
*---------------------------------------------------------------------*
          IF (NBATCH.GT.1) THEN
            DTIME = SECOND()
            CALL CCFSAVE(IBATCH,  IRHGH, I2HGH, INTMEDR, INTMED2,
     &                   KRHO2,   LUBF,  BFFIL, LENBF,
     &                   KFOCK,   LUFK,  FILFCK,LENFK,
     &                   KRIM,    LURIM, FILRIM,LENR,
     &                   KFIM,    LUFIM, FIMFIL,LENFIM,
     &                   KGIM,    LUGIM, GIMFIL,LENGIM,
     &                   KRHO1,   LURHO, RHOFIL,LENRHO,
     &                   NINTR,   NINT2, WORK,  LWORK )
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
        CALL CCFSAVE(1,       IRHGH, I2HGH, INTMEDR, INTMED2,
     &               KRHO2,   LUBF,  BFFIL, LENBF,
     &               KFOCK,   LUFK,  FILFCK,LENFK,
     &               KRIM,    LURIM, FILRIM,LENR,
     &               KFIM,    LUFIM, FIMFIL,LENFIM,
     &               KGIM,    LUGIM, GIMFIL,LENGIM,
     &               KRHO1,   LURHO, RHOFIL,LENRHO,
     &               NINTR,   NINT2, WORK,  LWORK )
        TIMIO = TIMIO + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
* close & delete some files which are no longer needed:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (CCS.OR.CC2) THEN
         CALL WCLOSE2(LUZDN, FILZDN,'DELETE')
      ELSE IF (.NOT.(CCS.OR.CC2)) THEN
         CALL WCLOSE2(LUMGD, FILMGD,'DELETE')
         CALL WCLOSE2(LUPQR0,FNPQR0,'DELETE')
         CALL WCLOSE2(LUPQRR,FNPQRR,'DELETE')
         CALL WCLOSE2(LUBDZ0,FNBDZ0,'DELETE')
         CALL WCLOSE2(LUO3 , FILO3, 'DELETE')
      END IF

      TIMPRE = TIMPRE + SECOND() - DTIME

*---------------------------------------------------------------------*
* recover workspace:
*---------------------------------------------------------------------*
      KEND1 = KEND0
      LWRK1 = LWRK0


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) MSGDBG,'Loop over AO-integrals completed ',
     &           ' & AO intermediates saved on file.'
        WRITE (LUPRI,*) MSGDBG,'recover work space: KEND1,LWRK1=',
     &       KEND1,LWRK1
        WRITE (LUPRI,*) MSGDBG,'norm of XLAMH0:',
     &        DDOT(NLAMDT,WORK(KLAMH0),1,WORK(KLAMH0),1)
        CALL FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
* correct response BF intermediate for CCSD(R12) and higher           *
*---------------------------------------------------------------------*
      IF (CCSDR12) THEN
        DO IINTR = 1, NINTR
          IDLSTR = INTMEDR(1,IINTR)
          ISYM   = ILSTSYM(LISTR,IDLSTR)
          !allocate scratch memory
          KSCR1 = KEND1
          KEND2 = KSCR1 + NT2AO(ISYM) 
C         KEND2 = KSCR1 + 2*NT2ORT(ISYAMP) 
          LWRK2 = LWORK - KEND2
          IF (LWRK2 .LE. 0) THEN
            CALL QUIT('Insufficient work space in CC_FMAT. (7)')
          END IF
          
          CALL DZERO(WORK(KSCR1),NT2AO(ISYM))
C         CALL DZERO(WORK(KSCR1),2*NT2ORT(ISYAMP))

          !calculate contribution:
          IOPT = 1
C         IOPT = 0
          IAMP = 2
          CALL CCRHS_BP(WORK(KSCR1),ISYM,IOPT,IAMP,DUMMY,IDUMMY,
     &                  IDUMMY,LISTR,IDLSTR,KETSCL,WORK(KEND2),LWRK2)

          !read in response BF intermediate
          KSCR2 = KEND2 
          KEND2 = KSCR2 + NT2AOIJ(ISYM)
          LWRK2 = LWORK - KEND2
          IF (LWRK2 .LE. 0) THEN
            CALL QUIT('Insufficient work space in CC_FMAT. (7)')
          END IF
          DTIME = SECOND()
          CALL CC_RVEC(LUBF,BFFIL,LENBF,NT2AOIJ(ISYM),IINTR,
     &                 WORK(KSCR2))
          TIMIO = TIMIO + SECOND() - DTIME

          !transform beta index to occupied and add on BF intermediate:
          OSQSAV = OMEGSQ
          OORSAV = OMEGOR 
          OMEGSQ = .FALSE.
          OMEGOR = .FALSE.
C         OMEGOR = .TRUE.
          ICON = 4
          CALL CC_T2MO(DUMMY,DUMMY,1,WORK(KSCR1),DUMMY,WORK(KSCR2),
     &                 WORK(KLAMP0),WORK(KLAMP0),1,WORK(KEND2),LWRK2,
     &                 ISYM,ICON)
          OMEGSQ = OSQSAV
          OMEGOR = OORSAV

          !write out response BF intermediate
          DTIME = SECOND()
          CALL CC_WVEC(LUBF,BFFIL,LENBF,NT2AOIJ(ISYM),IINTR,
     &                 WORK(KSCR2)) 
          TIMIO = TIMIO + SECOND() - DTIME

        END DO
      END IF      
      
*---------------------------------------------------------------------*
* open files for CBAR, DBAR, X & Y intermediates:
*---------------------------------------------------------------------*
      LUCBAR = -1
      LUDBAR = -1
      LUPQMO = -1
      LUXIM  = -1
      LUYIM  = -1
      LUMIM  = -1
      IF (.NOT.(CCS.OR.CC2)) THEN
         CALL WOPEN2(LUCBAR,CBAFIL, 64,0)
         CALL WOPEN2(LUDBAR,DBAFIL, 64,0)
         CALL WOPEN2(LUPQMO,FNPQMO, 64,0)
      END IF

      IF (.NOT.CCS) THEN
         CALL WOPEN2(LUXIM, FILXIM, 64,0)
         CALL WOPEN2(LUYIM, FILYIM, 64,0)
      END IF

      IF (CCSDR12) CALL WOPEN2(LUMIM, FILMIM, 64,0)

*---------------------------------------------------------------------*
* save a copy of the zeroth-order AO Fock matrix:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      CALL DCOPY(N2BAST,WORK(KFOCK0),1,WORK(KFOCK0AO),1)

      TIMFCK = TIMFCK + SECOND() - DTIME
      TIMIM0 = TIMIM0 + SECOND() - DTIME

*---------------------------------------------------------------------*
* do some printout:
*---------------------------------------------------------------------*
      IF (IPRINT.GT.2) THEN

         WRITE (LUPRI,'(1X,A,F10.2,A)')
     &    '| time for zero order intermediat.:',TIMIM0 ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINTL,' sets of left interm.:',TIMIML ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINTR,' sets of right inter.:',TIMIMR ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A,F10.2,A)')
     &    '| time for',NINT2,' sets of 2. ord. int.:',TIMIM2 ,' secs.|'
         WRITE (LUPRI,'(1X,A,I3,A)')
     &    '| intermediates calculated in ',NBATCH,' batches          |'

         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'

         IF (IOPTRES.EQ.5) THEN
            WRITE (LUPRI,'(1X,A)')
     &         '| L vector | R vector |  # products  |             |'
            WRITE (LUPRI,'(1X,3(A,A3),A)') '|  ',LISTL(1:3), ' No. |  ',
     &       LISTR(1:3),' No. |   with ', FILFMA(1:3),
     &       '   |  time/secs  |'
         ELSE 
            WRITE (LUPRI,'(1X,A)')
     &         '| L vector | R vector |    result    |             |'
            WRITE (LUPRI,'(1X,A2,A3,2(A,A3),A)') '|  ',LISTL(1:3), 
     &        '  No. |  ', LISTR(1:3),' No. |   ', FILFMA(1:3),
     &        '  No.   |  time/secs  |'
         END IF

         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
      END IF

*=====================================================================*
* calculate F matrix transformations:
*=====================================================================*
      IADRTH = 1
      DO ITRAN = 1, NFTRAN

        IDLSTL = IFTRAN(1,ITRAN)
        IDLSTR = IFTRAN(2,ITRAN)

        ISYCTR = ILSTSYM(LISTL,IDLSTL)
        ISYAMP = ILSTSYM(LISTR,IDLSTR)
        ISYRES = MULD2H(ISYAMP,ISYCTR)

        IINTL  = ICCSET2(INTMEDL,LISTL,IDLSTL,
     &                          'R0 ',0,      NINTL,MAXSIM,NOAPPEND)
        IINTR  = ICCSET1(INTMEDR,LISTR,IDLSTR,NINTR,MAXSIM,NOAPPEND)
        IINT2  = ICCSET2(INTMED2,LISTL,IDLSTL,
     &                           LISTR,IDLSTR,NINT2,MAXSIM,NOAPPEND)

        TIMTRN = SECOND()

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG,'F matrix transformation for ITRAN,',
     &          ITRAN
           WRITE (LUPRI,*) MSGDBG,'IADRTH:',IADRTH
           WRITE (LUPRI,*) MSGDBG,'LISTL,IDLSTL:',LISTL(1:3),IDLSTL
           WRITE (LUPRI,*) MSGDBG,'LISTR,IDLSTR:',LISTR(1:3),IDLSTR
           WRITE (LUPRI,*) MSGDBG,'ISYCTR,ISYAMP,ISYRES:',ISYCTR,ISYAMP,
     &          ISYRES
           WRITE (LUPRI,*) MSGDBG,'IINTL,IINTR,IINT2:',IINTL,IINTR,IINT2
        END IF

*---------------------------------------------------------------------*
* read the single excitation part of the response vectors and 
* lagrangian multipliers and calculate the response Lambda matrices:
*---------------------------------------------------------------------*
        KTHETA1 = KEND1
        KEND2   = KTHETA1 + NT1AM(ISYRES)
        KTHETA2 = KDUM

        CALL DZERO(WORK(KTHETA1),NT1AM(ISYRES))

        DTIME = SECOND()

        KZETA1  = KEND2
        KT1AMPA = KZETA1  + NT1AM(ISYCTR)
        KLAMDPA = KT1AMPA + NT1AM(ISYAMP)
        KLAMDHA = KLAMDPA + NGLMDT(ISYAMP)
        KEND2   = KLAMDHA + NGLMDT(ISYAMP)
        LWRK2   = LWORK - KEND2

        IF (LWRK2 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CC_FMAT. (8)')
        END IF
      

        IOPT = 1

        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                WORK(KZETA1),WORK(KDUM))

        CALL CC_RDRSP(LISTR,IDLSTR,ISYAMP,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KDUM))

        CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMDPA), WORK(KLAMH0),
     &                   WORK(KLAMDHA),WORK(KT1AMPA),ISYAMP)

        TIMPRE = TIMPRE + SECOND() - DTIME

*---------------------------------------------------------------------*
*       do some memory allocation and read (ia|jb) integrals:
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN
           KXBARA = KEND2
           KYBARA = KXBARA + NMATIJ(ISYAMP)
           KA2IM  = KYBARA + NMATAB(ISYAMP)
           KEND2  = KA2IM  + NT1AM(ISYRES)
  
           KLIAJB = KEND2
           KXIAJB = KLIAJB + NT2SQ(ISYM0)
           KEND3  = KXIAJB + NT2AM(ISYM0)
        ELSE IF (LISTL(1:2).EQ.'L0') THEN
           KXIAJB = KEND2
           KEND3  = KXIAJB + NT2AM(ISYM0)
        ELSE
           KEND3  = KEND2
        END IF

        LWRK3  = LWORK - KEND3

        IF (LWRK3 .LE. 0) THEN
           CALL QUIT('Insufficient work space in CC_FMAT.')
        END IF

C       -----------------------
C       read (ia|jb) integrals:
C       -----------------------
        IF ((.NOT.CCS) .OR. LISTL(1:2).EQ.'L0') THEN
           CALL CCG_RDIAJB(WORK(KXIAJB),NT2AM(ISYM0))
        END IF

C       --------------------------------------------------------
C       for CCSD resort (ia|jb) integrals to full square matrix:
C       --------------------------------------------------------
        IF (.NOT.CCS) THEN
           CALL CC_T2SQ(WORK(KXIAJB),WORK(KLIAJB),ISYM0)
        END IF

*---------------------------------------------------------------------*
*       calculate the projection on <HF|:
*---------------------------------------------------------------------*
        IF (LISTL(1:2).EQ.'L0') THEN

C          -----------------------------------------------------------
C          calculate packed L(ia,jb) integrals and evaluate the 
C          projection contribution <HF|[[H,T],t_mu1]|CC>
C          -----------------------------------------------------------
           IOPTTCME = 1
           CALL CCSD_TCMEPK(WORK(KXIAJB),ONE,ISYM0,IOPTTCME)

           IOPT = 0
           CALL CCG_LXD(WORK(KTHETA1),ISYRES,WORK(KT1AMPA),ISYAMP,
     &                  WORK(KXIAJB),ISYM0,IOPT)
           CALL DSCAL(NT1AM(ISYRES),TWO,WORK(KTHETA1),1)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) MSGDBG,
     &             'norm of THETA1 after <HF| contrib.:',
     &          DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
C              WRITE (LUPRI,*) 
C     &           MSGDBG,'result vector after <HF| contribution:'
C              CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)
           END IF
        END IF

*---------------------------------------------------------------------*
*       for CCSD calculate second 21G contribution
*---------------------------------------------------------------------*
        IF (CCSD .OR. CCSDT) THEN

C          ----------------------------------------------------------
C          memory allocation, reusing space for packed (ia|jb) array:
C          ----------------------------------------------------------
           KXIDJL  = KLIAJB  + NT2SQ(ISYM0)
           KXJLID  = KXIDJL  + NTRAOC(ISYAMP)
           KEND3   = KXJLID  + NTRAOC(ISYAMP)
           LWRK3   = LWORK - KEND3

           IF (LWRK3 .LE. 0) THEN
             CALL QUIT('Insufficient work space in CC_FMAT. (21G)')
           END IF

C          ----------------------------------------------------------
C          calculate g_iljd = (ia|jd) * T_al and compute 21G contrib.
C          ----------------------------------------------------------
           CALL CCG_TRANS4(WORK(KXIDJL),ISYAMP,WORK(KLIAJB),ISYM0,
     &                     WORK(KT1AMPA),ISYAMP)

           IOPT = 5
           CALL CCG_SORT1(WORK(KXIDJL),WORK(KXJLID),ISYAMP,IOPT)

           CALL CC_21GMO(WORK(KTHETA1),ISYRES,WORK(KXJLID),ISYAMP,
     &                   LUPQMO,FNPQMO,IADRPQMO(1,IINTL),ISYCTR,
     &                   WORK(KEND3),LWRK3)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'Theta1 after 21GMO contribution:'
              CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)
           END IF
        END IF

*---------------------------------------------------------------------*
*       for CC2 & CCSD precalculate the YBARA and A2 intermediates
*       for CC2 precalculate also the XBARA intermediate
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN

           KT2AMPA = KLIAJB  + NT2SQ(ISYM0)
           KEND3   = KT2AMPA + MAX(NT2AM(ISYAMP),NT2AM(ISYM0))
           LWRK3   = LWORK - KEND3

           IF (LWRK3 .LE. 0) THEN
             CALL QUIT('Insufficient work space in CC_FMAT. '//
     &                 '(YBARA/A2IM)')
           END IF

C          ----------------------------------------------
C          calculate Liajb and read response amplitudes:
C          ----------------------------------------------
           CALL CCRHS_T2TR(WORK(KLIAJB),WORK(KEND3),LWRK3,ISYM0)
         
           IOPT = 2
           CALL CC_RDRSP(LISTR,IDLSTR,ISYAMP,IOPT,MODEL,
     &                   WORK(KDUM),WORK(KT2AMPA))
           CALL CCLR_DIASCL(WORK(KT2AMPA),TWO,ISYAMP)

C          ----------------------------------------------
C          XBARA_kj = sum_dem T_djem L_mekd
C          ----------------------------------------------
           IF (CC2) THEN
              CALL CC_XBAR(WORK(KXBARA),WORK(KLIAJB),ISYM0,
     &                     WORK(KT2AMPA),ISYAMP,WORK(KEND3),LWRK3)
           END IF

C          ----------------------------------------------
C          YBARA_ba = sum_dlk T_dlbk L_ldka
C          ----------------------------------------------
           CALL CC_YBAR(WORK(KYBARA),WORK(KLIAJB),ISYM0,
     &                  WORK(KT2AMPA),ISYAMP,WORK(KEND3),LWRK3)

C          ----------------------------------------------
C          A2IM = sum_bj (2T_aibj - T_ajbi) Zeta_bj
C          ----------------------------------------------
           IOPTTCME = 1
           CALL CCSD_TCMEPK(WORK(KT2AMPA),ONE,ISYAMP,IOPTTCME)
           CALL CCG_LXD(WORK(KA2IM),ISYRES,WORK(KZETA1),ISYCTR,
     &                  WORK(KT2AMPA),ISYAMP,0)

        END IF

*---------------------------------------------------------------------*
* allocate & initialize doubles result vector:
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN
          KTHETA2 = KEND2
          KEND2   = KTHETA2 + NT2AM(ISYRES)
          CALL DZERO(WORK(KTHETA2),NT2AM(ISYRES))
        END IF

*---------------------------------------------------------------------*
* allocate work space some `small' intermediates:
*---------------------------------------------------------------------*
        IF (CCS.OR.CC2) THEN
          KFOCKA  = KEND2
          KGZETA  = KFOCKA + N2BST(ISYAMP)
          KEND2   = KGZETA + N2BST(ISYRES)
        END IF

        IF (.NOT.CCS) THEN
          KRIMA   = KEND2
          KXINT   = KRIMA + NEMAT1(ISYAMP)
          KYINT   = KXINT + NMATIJ(ISYRES)
          KFINT   = KYINT + NMATAB(ISYRES)
          KEND2   = KFINT + NT1AO(ISYRES) 
        END IF

        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LE. NT1AM(ISYRES)) THEN
          CALL QUIT('Insufficient work space in CC_FMAT. (8)')
        END IF

*---------------------------------------------------------------------*
* for CCS/CC2 read the response Fock matrix and the Fock matrix like
* CCS/CC2 GZeta intermediate and calculate the CC_11A contribution
*---------------------------------------------------------------------*
        IF (CCS.OR.CC2) THEN
           IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'Norm^2 of THETA1 before GZETA contrib.: ',
     &          DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
           END IF
       
           CALL CC_RVEC(LUFK,FILFCK,LENFK,N2BST(ISYAMP),IINTR,
     &                  WORK(KFOCKA))
           CALL CC_FCKMO(WORK(KFOCKA),WORK(KLAMP0),WORK(KLAMH0),
     &                   WORK(KEND2),LWRK2,ISYAMP,ISYM0,ISYM0)

           CALL CC_RVEC(LUGIM,GIMFIL,LENGIM,N2BST(ISYRES),IINT2,
     &                  WORK(KGZETA))
           CALL CC_11A(WORK(KTHETA1),WORK(KGZETA),ISYRES,WORK(KLAMH0),
     &                 WORK(KLAMP0),WORK(KEND2),LWRK2)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'CCS GZETA intermediate:'
              CALL CC_PRFCKAO(WORK(KGZETA),ISYRES)
              WRITE (LUPRI,*) 'Norm^2 of THETA1 after GZETA contrib.: ',
     &          DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
           END IF
        END IF
      
*---------------------------------------------------------------------*
* for CC2 and CCSD read precalculated X and Y intermediates from file,
* and add the precalculated contributions from Rho1 to result vector:
*---------------------------------------------------------------------*
        IF (.NOT.CCS) THEN
           CALL CC_RVEC(LUXIM,FILXIM,LENX,NMATIJ(ISYRES),IINT2,
     &                  WORK(KXINT))

           CALL CC_RVEC(LUYIM,FILYIM,LENY,NMATAB(ISYRES),IINT2,
     &                  WORK(KYINT))
           
           CALL CC_RVEC(LURHO,RHOFIL,LENRHO,NT1AM(ISYRES),IINT2,
     &                  WORK(KEND2))

           CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KEND2),1,WORK(KTHETA1),1)

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'restored X and Y intermediates:'
              CALL CC_PREI(WORK(KYINT),WORK(KXINT),ISYRES,1)
              WRITE (LUPRI,*) LUXIM,FILXIM,LENX,NMATIJ(ISYRES),IINT2
              WRITE (LUPRI,*) LUYIM,FILYIM,LENY,NMATAB(ISYRES),IINT2

              WRITE (LUPRI,*) 'Rho1 read from file:'
              CALL CC_PRP(WORK(KEND2),WORK(KEND2),ISYRES,1,0)
              WRITE (LUPRI,*) 'result vector after Rho1 was added:'
              CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)

              CALL FLSHFO(LUPRI) 
           END IF

        END IF

*---------------------------------------------------------------------*
* for CCSD read also the F and R intermediates;
* the F intermediate is transformed to MO and added to result vector
*---------------------------------------------------------------------*
        IF (.NOT.(CCS.OR.CC2)) THEN
           CALL CC_RVEC(LUFIM,FIMFIL,LENFIM,NT1AO(ISYRES),IINT2,
     &                  WORK(KFINT))
           CALL CC_T1AM(WORK(KTHETA1),ISYRES,WORK(KFINT),ISYRES,
     &                  WORK(KLAMH0),ISYM0,ONE)

           CALL CC_RVEC(LURIM,FILRIM,LENR,NEMAT1(ISYAMP),IINTR,
     &                  WORK(KRIMA))
        END IF

*---------------------------------------------------------------------*
* calculate the remaining contributions:
*---------------------------------------------------------------------*

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG,'norm of THETA1 before CCFMAT2:',
     &       DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
           IF (.NOT.CCS) 
     &      WRITE (LUPRI,*) MSGDBG,'norm of THETA2 before CCFMAT2:',
     &       DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
           WRITE (LUPRI,*) 'result vector before entering CCFMAT2:'
           CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,0)
        END IF

        CALL CCFMAT2(WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &               LISTL,IDLSTL, WORK(KZETA1), ISYCTR,
     &               LISTR,IDLSTR, WORK(KT1AMPA),ISYAMP,
     &               WORK(KLAMDPA),WORK(KLAMDHA),
     &               WORK(KLAMP0), WORK(KLAMH0),
     &               WORK(KXINT),  WORK(KYINT),
     &               WORK(KXBARA), WORK(KYBARA), 
     &               WORK(KRIMA),  WORK(KA2IM),
     &               WORK(KFOCKA), WORK(KFOCK0AO), WORK(KFCKC0),
     &               LUBFI,FNBFI,IADRBFI(1,IINT2),
     &               LUBF, BFFIL,LENBF,IINTR,
     &               LUC,  CTFIL, LUCBAR, CBAFIL,
     &               LUD,  DTFIL, LUDBAR, DBAFIL,
     &               IOFFCD(IINTR), WORK(KEND2),LWRK2)

*---------------------------------------------------------------------*
* allocate memory for R12-part of result vector (stored as symmetry 
* packed triangular matrix):
*---------------------------------------------------------------------*
        IF (CCR12) THEN
          KTHETAR12 = KEND2
          KEND2    = KTHETAR12 + NTR12AM(ISYRES)
          LWRK2    = LWORK - KEND2
          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient memory for R12-FMATRIX!')
          END IF
          CALL DZERO(WORK(KTHETAR12),NTR12AM(ISYRES))

        END IF
C
*---------------------------------------------------------------------*
* calculate R12 contributions:
* 
* C. Neiss,  Autumn 2004
*---------------------------------------------------------------------*
C
       IF (CCR12) THEN
         IF (LOCDBG) THEN
           CALL AROUND('Call R12 section for FMATRIX')
         END IF
C        
         CALL CC_R12FMAT(WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KTHETAR12), ISYRES,
     &                   LISTL,IDLSTL, WORK(KZETA1), ISYCTR,
     &                   LISTR,IDLSTR, ISYAMP,
     &                   WORK(KLAMDPA),WORK(KLAMDHA),
     &                   WORK(KLAMP0), WORK(KLAMH0),
     &                   WORK(KXINT),LUMIM,FILMIM,LENM,IINT2,
     &                   WORK(KEND2),LWRK2)
         IF (LOCDBG) THEN
           CALL AROUND('Returned from R12 section for FMATRIX')
         END IF
C
       ENDIF
C
C-------------------------------------------
C     SLV98,OC, CCMM JK, NYQMMM KS 10
C     Calculate T^{gB} solvent contribution.
C     Unsquared T2 is required!
C-------------------------------------------
C
      IF (CCSLV.AND.LSTBTR) THEN
C
         KT2AMPA = KEND2
         KEND3   = KT2AMPA +  NT2AM(ISYAMP)
         LWRK3   = LWORK - KEND3
         IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insuff. work in CC_FMAT. in CCSLV part')
         END IF
C
C
         IF (LISTL(1:2).EQ.'L0') THEN
           IOPT = 2
           CALL CC_RDRSP(LISTR,IDLSTR,ISYAMP,IOPT,MODEL,
     *                   WORK(KDUM),WORK(KT2AMPA))

           IF (.NOT. CCMM) CALL CCSL_LTRB(WORK(KTHETA1),WORK(KTHETA2),
     *                    WORK(KT1AMPA),WORK(KT2AMPA),ISYAMP,'F',
     *                    WORK(KEND3),LWRK3)
           IF (CCMM) THEN
              IF (.NOT. NYQMMM) THEN
                 CALL CCMM_LTRB(WORK(KTHETA1),WORK(KTHETA2),
     *                    WORK(KT1AMPA),WORK(KT2AMPA),ISYAMP,'F',
     *                    WORK(KEND3),LWRK3)
              ELSE IF (NYQMMM .AND. (.NOT. HFFLD)) THEN
                 CALL CCMM_TRANSFORMER(WORK(KTHETA1),WORK(KTHETA2),
     *                         WORK(KT1AMPA),WORK(KT2AMPA),MODEL,ISYAMP,
     *                         'F',WORK(KEND3),LWRK3)
              END IF
           END IF
         ELSE IF (LISTL(1:2).NE.'L0') THEN
           IF (.NOT. CCMM) CALL CCSL_PBTR(WORK(KTHETA1),WORK(KTHETA2),
     *                                    ISYRES,
     *                                    LISTL,IDLSTL, WORK(KZETA1),
     *                                    ISYCTR,
     *                                    LISTR,IDLSTR, WORK(KT1AMPA),
     *                                    ISYAMP,MODEL,
     *                                    WORK(KEND2),LWRK2)
           IF (CCMM) THEN
             IF (.NOT. NYQMMM) THEN
               CALL CCMM_PBTR(WORK(KTHETA1),WORK(KTHETA2),
     *                              ISYRES,
     *                              LISTL,IDLSTL, WORK(KZETA1),
     *                              ISYCTR,
     *                              LISTR,IDLSTR, WORK(KT1AMPA),
     *                              ISYAMP,MODEL,
     *                              WORK(KEND2),LWRK2)
             ELSE IF (NYQMMM .AND. (.NOT. HFFLD ) ) THEN
               RSPTYP='F'
               CALL CCMM_QRTRANSFORMER(WORK(KTHETA1),WORK(KTHETA2),
     *                              ISYRES,
     *                              LISTL,IDLSTL,ISYCTR,
     *                              LISTR,IDLSTR,ISYAMP,
     *                              MODEL,RSPTYP,WORK(KEND2),LWRK2)
             END IF
           END IF 
         END IF
      ENDIF
C
C----------------------------------------------------------
C
        
        IF (CCSDT) THEN

           ! allocate memory for effective rhs vector
           IF ( IOPTRES.GE.0 .AND. IOPTRES.LE.4 ) THEN
             KTHETA1EFF = KEND2
             KTHETA2EFF = KTHETA1EFF + NT1AM(ISYRES)
             KEND2      = KTHETA2EFF + NT2AM(ISYRES)
             LWRK2      = LWORK - KEND2
             IF (LWRK2 .LE. NT1AM(ISYRES)) THEN
               CALL QUIT('Insufficient work space in CC_FMAT. (9)')
             END IF
           END IF
 
           IF ( (.NOT. NODDY_FMAT)) THEN 
 
C            ----------------------------------------
C            Calculate <L_2|[[H,T^B_3],\tau_nu_1|HF>
C            ----------------------------------------
             IF (IOPTRES.LT.5) THEN

               CALL CC3_FT3B(LISTL,IDLSTL,LISTR,IDLSTR,IOPTRES,
     *                       WORK(KTHETA1),ISYRES,WORK(KEND2),LWRK2,
     *                       FNCKJD,FNDKBC,FNDELD,FNTOC)

             END IF
 
             LUCKJD   = -1
             LUDKBC   = -1
             LUTOC    = -1
             LU3VI    = -1
             LU3VI2   = -1
             LUDKBC3  = -1
             LU3FOP   = -1
             LU3FOP2  = -1
             LU3FOPX  = -1
             LU3FOP2X = -1
 
             CALL WOPEN2(LUCKJD,FNCKJD,64,0)
             CALL WOPEN2(LUDKBC,FNDKBC,64,0)
             CALL WOPEN2(LUTOC,FNTOC,64,0)
             CALL WOPEN2(LU3VI,FN3VI,64,0)
             CALL WOPEN2(LU3VI2,FN3VI2,64,0)
             CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
             CALL WOPEN2(LU3FOP,FN3FOP,64,0)
             CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
             CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
             CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)

C
C            ---------------------------------------------
C             Calculate   <L_3|[[H^B,T^0_2],\tau_nu_1]|HF> 
C                       + <L_3|[[H^0,T^B_2],\tau_nu_1]|HF> 
C             and 
C                         <L_3|[H^B,\tau_nu_2]|HF>
C            ---------------------------------------------
             IF (IOPTRES.EQ.5) THEN
               ! in  CC3_BFMAT Omegaeff is used as scratch array.
               KTHETA1EFF = KEND2
               KTHETA2EFF = KTHETA1EFF + NT1AM(ISYRES)
               KEND3      = KTHETA2EFF + NT2AM(ISYRES)
               LWRK3      = LWORK - KEND3
   
               IF (LWRK3 .LT. 0) THEN
                 CALL QUIT('Insufficient work space in CC_FMAT. (9b)')
               END IF

             ELSE
               KEND3 = KEND2
               LWRK3 = LWRK2
             END IF
 
             CALL DZERO(WORK(KTHETA1EFF),NT1AM(ISYRES))
             CALL DZERO(WORK(KTHETA2EFF),NT2AM(ISYRES))

             CALL CC3_BFMAT(LISTL,IDLSTL,LISTR,IDLSTR,
     *                     WORK(KTHETA1),   WORK(KTHETA2),
     *                     WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     *                     ISYRES,
     *                     WORK(KEND3),LWRK3,
     *                     LUTOC,FNTOC,
     *                     LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                     LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                     LU3VI2,FN3VI2,LU3FOP,FN3FOP,LU3FOP2,FN3FOP2)

C
C            ------------------------------------------
C             Calculate <F_3|[[H,T^0_2],\tau_nu_1]|HF>
C             and       <F_3|[H,\tau_nu_2]|HF>
C            ------------------------------------------
             IF (IOPTRES.LT.5) THEN
 
               CALL DZERO(WORK(KTHETA1EFF),NT1AM(ISYRES))
               CALL DZERO(WORK(KTHETA2EFF),NT2AM(ISYRES))

               CALL CC3_FMATSD(LISTL,IDLSTL,LISTR,IDLSTR,
     *                     WORK(KTHETA1),WORK(KTHETA2),WORK(KTHETA1EFF),
     *                     WORK(KTHETA2EFF),ISYRES,WORK(KEND2),LWRK2,
     *                     LUDKBC,FNDKBC,LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                     LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                     LU3FOP,FN3FOP,LU3FOP2,FN3FOP2)
             ENDIF

             CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
             CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
             CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
             CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
             CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
             CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
             CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
             CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')
             CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
             CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')
 
           ELSE
 
             IF ( IOPTRES.GE.0 .AND. IOPTRES.LE.4 ) THEN
               CALL DZERO(WORK(KTHETA1EFF),NT1AM(ISYRES))
               CALL DZERO(WORK(KTHETA2EFF),NT2AM(ISYRES))
             END IF
 
             CALL CCSDT_FMAT_NODDY(LISTL,IDLSTL,LISTR,IDLSTR,
     &                             IOPTRES,CCSDT_F_ALTER,
     &                             WORK(KTHETA1),   WORK(KTHETA2),
     &                             WORK(KTHETA1EFF),WORK(KTHETA2EFF),
     &                             IFDOTS,FCONS,FILFMA,ITRAN,
     &                             NFTRAN,MXVEC,WORK(KEND2),LWRK2)
 
           ENDIF
 
           IF (IPRINT .GT. 110) THEN
              CALL AROUND('Omega1 and Omega2 after  F(cc3) : ')
              CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),1,1,1)
           ENDIF

        END IF



*---------------------------------------------------------------------*
* write result vector to output:
*---------------------------------------------------------------------*
      DTIME = SECOND()

      KTHETA0 = -999999

      IF (IOPTRES .EQ. 0  .OR. IOPTRES .EQ. 1) THEN

*       write to a common direct access file, 
*       store start address in IFTRAN(3,ITRAN)

        IFTRAN(3,ITRAN) = IADRTH

        CALL PUTWA2(LUFMAT,FILFMA,WORK(KTHETA1),IADRTH,NT1AM(ISYRES))
        IADRTH = IADRTH + NT1AM(ISYRES)

        IF (.NOT.CCS) THEN
          CALL PUTWA2(LUFMAT,FILFMA,WORK(KTHETA2),IADRTH,NT2AM(ISYRES))
          IADRTH = IADRTH + NT2AM(ISYRES)
        END IF

        IF (CCR12) THEN
          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'THETAR12 in CC_FMATNEW:'
            WRITE(LUPRI,*) 'Will write to file ', FILFMA, LUFMAT
            CALL OUTPAK(WORK(KTHETAR12),NMATKI(ISYRES),1,LUPRI)
          END IF
          CALL PUTWA2(LUFMAT,FILFMA,WORK(KTHETAR12),IADRTH,
     &                NTR12AM(ISYRES))
          IADRTH = IADRTH + NTR12AM(ISYRES)
        END IF

C       some debug output:
        IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'F matrix transformation nb. ',ITRAN,
     &          ' saved on file.'
         WRITE (LUPRI,*) 'ADRESS, LENGTH:',
     &        IFTRAN(3,ITRAN),IADRTH-IFTRAN(3,ITRAN)
         XNORM = DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
         IF (.NOT.CCS) XNORM = XNORM +
     &           DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
         IF (CCR12) XNORM = XNORM + 
     &           DDOT(NTR12AM(ISYRES),WORK(KTHETAR12),1,
     &                WORK(KTHETAR12),1)
         WRITE (LUPRI,*) 'Norm:', XNORM

         Call AROUND('F matrix transformation written to file:')
         NDBLE = 1
         IF (CCS) NDBLE = 0
         Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,NDBLE)
         IF (CCR12) THEN
           write(lupri,*) 'THETA_R12 in CC_FMATNEW:'
           CALL OUTPAK(WORK(KTHETAR12),NMATKI(ISYRES),1,LUPRI)
         END IF
        END IF
C

      ELSE IF ( IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4 ) THEN
         ! write to a sequential file by a call to CC_WRRSP/CC_WARSP,
         ! use FILFMA as LIST type and IFTRAN(3,ITRAN) as index
         IF (IOPTRES.EQ.3) THEN
           CALL CC_WRRSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LWRK2)
           IF (CCR12) THEN
             CALL CC_WRRSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTWR12,
     &                     MODELW,DUMMY,DUMMY,WORK(KTHETAR12),
     &                     WORK(KEND2),LWRK2)
           END IF
           IF (CCSDT) THEN
             CALL CC_WRRSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTWE,MODELW,
     &                     WORK(KTHETA0),WORK(KTHETA1EFF),
     &                     WORK(KTHETA2EFF),WORK(KEND2),LWRK2)
           END IF
         ELSE IF (IOPTRES.EQ.4) THEN
           CALL CC_WARSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTW,MODELW,
     &                   WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                   WORK(KEND2),LWRK2)
           IF (CCR12) THEN
             CALL CC_WARSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTWR12,
     &                     MODELW,DUMMY,DUMMY,WORK(KTHETAR12),
     &                     WORK(KEND2),LWRK2)
           END IF
           IF (CCSDT) THEN
             CALL CC_WARSP(FILFMA,IFTRAN(3,ITRAN),ISYRES,IOPTWE,MODELW,
     &                     WORK(KTHETA0),WORK(KTHETA1EFF),
     &                     WORK(KTHETA2EFF),WORK(KEND2),LWRK2)
           END IF
         END IF

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Write ',LISTL(1:3),' * B * ',LISTR(1:3),
     &      ' transformation',' as ',FILFMA(1:3),' type vector to file.'
           WRITE (LUPRI,*) 'index of inp. left  vector:',IFTRAN(1,ITRAN)
           WRITE (LUPRI,*) 'index of inp. right vector:',IFTRAN(2,ITRAN)
           WRITE (LUPRI,*) 'index of     result vector:',IFTRAN(3,ITRAN)
           XNORM = DDOT(NT1AM(ISYRES),WORK(KTHETA1),1,WORK(KTHETA1),1)
           IF (.NOT.CCS) XNORM = XNORM +
     &       DDOT(NT2AM(ISYRES),WORK(KTHETA2),1,WORK(KTHETA2),1)
           IF (CCR12) THEN
             XNORM = XNORM + DDOT(NTR12AM(ISYRES),WORK(KTHETAR12),1, 
     &               WORK(KTHETAR12),1)
           END IF 
           WRITE (LUPRI,*) 'norm^2 of result vector:',XNORM
           CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
           IF (CCR12) CALL CC_PRPR12(WORK(KTHETAR12),ISYRES,1,.TRUE.)
           IF (CCSDT) THEN
             WRITE(LUPRI,*) 'CCSDT effective vector:'
             XNORM = DDOT(NT1AM(ISYRES),WORK(KTHETA1EFF),1,
     &         WORK(KTHETA1EFF),1) + DDOT(NT2AM(ISYRES),
     &           WORK(KTHETA2EFF),1,WORK(KTHETA2EFF),1)
             WRITE (LUPRI,*) 'norm^2:',XNORM
             CALL CC_PRP(WORK(KTHETA1EFF),WORK(KTHETA2EFF),ISYRES,1,1)
           END IF
         END IF
      ELSE IF (IOPTRES.EQ.5) THEN
         IF (LOCDBG.AND.CCSDT) THEN
           WRITE(LUPRI,*) 'FCONS TRIPLES CONTRIBUTION:'
           IVEC = 1
           DO WHILE (IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            WRITE(LUPRI,*)'FCONS:',IVEC,ITRAN,FCONS(IVEC,ITRAN),IOPTW
            IVEC = IVEC + 1
           END DO
         END IF

         IF (CCR12) THEN
           CALL CCDOTRSP(IFDOTS,FCONS,IOPTWR12,FILFMA,ITRAN,NFTRAN,
     &                   MXVEC,DUMMY,WORK(KTHETAR12),ISYRES,
     &                   WORK(KEND2),LWRK2)
         END IF

         IF (LOCDBG) THEN
           write(lupri,*)'FCONS R12 doubles contributions:'
           IVEC = 1
           DO WHILE (IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            WRITE(LUPRI,*)'FCONS:',IVEC,ITRAN,FCONS(IVEC,ITRAN),IOPTW
            IVEC = IVEC + 1
           END DO
         END IF

         IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KTHETA2),TWO,ISYRES)
         CALL CCDOTRSP(IFDOTS,FCONS,IOPTW,FILFMA,ITRAN,NFTRAN,MXVEC,
     &                 WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &                 WORK(KEND2),LWRK2)
         IF (LOCDBG) THEN
           IVEC = 1
           DO WHILE (IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            WRITE(LUPRI,*)'FCONS:',IVEC,ITRAN,FCONS(IVEC,ITRAN),IOPTW
            IVEC = IVEC + 1
           END DO
         END IF
      ELSE
        CALL QUIT('Illegal value for IOPTRES in CC_FMAT.')
      END IF
      
      TIMIO = TIMIO + SECOND() - DTIME

      TIMTRN = SECOND() - TIMTRN

      IF (IPRINT.GT.2) THEN

         IF (IOPTRES.EQ.5) THEN
            IVEC = 1
            DO WHILE (IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
               IVEC = IVEC + 1
            END DO    
            WRITE (LUPRI,'(1X,2(A,I5),A,I6,A,F10.2,A)')'| ',IDLSTL,
     &        '    | ',IDLSTR,'    | ',IVEC-1,'       | ',TIMTRN,'  |'
         ELSE
            WRITE (LUPRI,'(1X,2(A,I5),A,I6,A,F10.2,A)') '| ',IDLSTL,
     &        '    | ',IDLSTR,'    | ',IFTRAN(3,ITRAN),'       | ',
     &           TIMTRN,'  |'
         END IF 

      END IF

*---------------------------------------------------------------------*
* End of loop over F matrix transformations
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
         CALL WCLOSE2(LURIM,  FILRIM, 'DELETE')
         CALL WCLOSE2(LUFIM,  FIMFIL, 'DELETE')
         CALL WCLOSE2(LUPQMO, FNPQMO, 'DELETE')
         CALL WCLOSE2(LUBFI,  FNBFI,  'DELETE')
      ELSE
         CALL WCLOSE2(LUGIM,  GIMFIL, 'DELETE')
         CALL WCLOSE2(LUFK,   FILFCK, 'DELETE')
      END IF

      IF (.NOT.CCS) THEN
         CALL WCLOSE2(LUXIM,FILXIM, 'DELETE')
         CALL WCLOSE2(LUYIM,FILYIM, 'DELETE')
         CALL WCLOSE2(LURHO,RHOFIL, 'DELETE')
      END IF
      IF (CCSD .OR. CCSDT)      CALL WCLOSE2(LUBFD, FNBFD, 'DELETE')

      IF (CCSDR12) CALL WCLOSE2(LUMIM,FILMIM,'DELETE')


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

        CALL GETWA2(LUFMAT,FILFMA,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
          DO ITRAN = 1, NFTRAN
            IF (ITRAN.LT.NFTRAN) THEN
              LEN     = IFTRAN(3,ITRAN+1)-IFTRAN(3,ITRAN)
            ELSE
              LEN     = IADRTH-IFTRAN(3,NFTRAN)
            END IF
            KTHETA1 = IFTRAN(3,ITRAN)
            XNORM   = DDOT(LEN, WORK(KTHETA1),1, WORK(KTHETA1),1)
            WRITE (LUPRI,*) 'Read F matrix transformation nb. ',NFTRAN
            WRITE (LUPRI,*) 'Adress, length, NORM:',IFTRAN(3,NFTRAN),
     &           LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

      TIMIO = TIMIO + SECOND() - DTIME
*---------------------------------------------------------------------*
* close F matrix file, print timings & return
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (IOPTRES.EQ.0 ) THEN
        CALL WCLOSE2(LUFMAT, FILFMA, 'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
        CALL WCLOSE2(LUFMAT, FILFMA, 'DELETE')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_FMAT.')
      END IF

      TIMIO  = TIMIO + SECOND() - DTIME
      TIMALL = SECOND() - TIMALL

      IF (IPRINT.GT.2) THEN
         WRITE (LUPRI,'(1X,A1,50("-"),A1)') '+','+'
         WRITE (LUPRI,'(1X,A,I4,A,F10.2,A)')
     &     '| total time for',NFTRAN,' F transforms.:',TIMALL,' secs.|'
        IF (TIMALL.GE.1.0D0) THEN
         CONVRT = 100.0D0/TIMALL
         WRITE (LUPRI,'(1X,A1,50("-"),A1)')'+','+'
         WRITE
     &    (LUPRI,'(1X,"|  % of time used in ",A,":",F10.2,"      |")')
     &      'start up org.', TIMPRE*CONVRT,
     &      'Fock interm. ', TIMFCK*CONVRT,
     &      'ERI          ', TIMINT*CONVRT,
     &      'CCRDAO       ', TIMRDAO*CONVRT,
     &      '(**|k del)   ', TIMTRBT*CONVRT,
     &      'BF term      ', TIMBF*CONVRT,
     &      'C term       ', TIMC*CONVRT,
     &      'GZ terms     ', TIMGZ*CONVRT,
     &      'FG terms     ', TIMFG*CONVRT,
     &      'D term       ', TIMD*CONVRT,
     &      'I/O          ', TIMIO*CONVRT
        END IF
        WRITE (LUPRI,'(1X,A1,50("="),A1,//)') '+','+'
      END IF
*---------------------------------------------------------------------*
* that's it; return:
*---------------------------------------------------------------------*

      CALL QEXIT('CC_FMATNEW')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CC_FMAT
*=====================================================================*
*=====================================================================*
C  /* Deck ccfpreint */
      SUBROUTINE CCFPREINT(INTMED2,NINT2,XLAMDP,XLAMDH,
     &                     FILXIM,LENX,
     &                     FILYIM,LENY,
     &                     FILMIM,LENM,
     &                     FILMGD,LENMGD,
     &                     FILZDN,LENZDN,
     &                     FILPQ,  IADRPQ,
     &                     FILPQMO,IADRPQMO,
     &                     FNBDZ0, IADRZ0,
     &                     TIMXYM,TIMBF,TIMIO,TIMPQ,
     &                     IOPTRSP,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: precalculate some intermediates for the F matrix 
*              transformation which depend on left and right vectors:
*
*              X, Y, M, Chi, Yps, P, Q, and the BZeta density 
*
*              the results are written to direct acces files
*
*     IOPTRSP = 1 computes zeroth-order version of BZeta density 
*                 and P and Q intermediates in MO and backtransformed
*                 X, Y are only used local
*
*     IOPTRSP = 2 computes F matrix version of BZeta density,
*                 and only backtransformed P and Q intermediates 
*                 X, Y are written to file
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
#include "ccsdinp.h"
#include "cciccset.h"
#include "ccfield.h"
#include "second.h"

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER LUXIM, LUYIM, LUBDZ0, LUPQ, LUPQMO, LUMGD, LUZDN, LUMIM

      CHARACTER*(*) FILYIM,FILXIM,FILMGD,FILPQ,FILPQMO,FNBDZ0,FILZDN,
     &              FILMIM
      INTEGER LWORK,IDLSTL,IDLSTR,NINT2,LENX,LENY,LENMGD,LENM,IOPTRSP
      INTEGER IADRPQ(MXCORB_CC,*), IADRZ0(MXCORB_CC,*), INTMED2(4,*)
      INTEGER IADRPQMO(MXCORB_CC,*), LENZDN

#if defined (SYS_CRAY)
      REAL WORK(LWORK), XLAMDP(*), XLAMDH(*)
      REAL TIMXYM, TIMBF, TIMIO, TIMPQ, DTIME
      REAL TWO, ONE, ZERO, DDOT, DUMMY, XNORM
#else
      DOUBLE PRECISION WORK(LWORK), XLAMDP(*), XLAMDH(*)
      DOUBLE PRECISION TIMXYM, TIMBF, TIMIO, TIMPQ, DTIME
      DOUBLE PRECISION TWO, ONE, ZERO, DUMMY
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
 
      CHARACTER*(10) MODEL
      CHARACTER*(3) LISTL, LISTR, LISTLOLD, LISTROLD
      LOGICAL NEWLVEC, NEWRVEC
      INTEGER IOPT, IOPTY, ISYML, ISYMR, ISYINT
      INTEGER KXINT, KYINT, KMINT, KMGDL, KCHI, KEND1, LWRK1, ISYM
      INTEGER ISYMK, ISYMC, ISYMI, KOFF1, KOFF2, KOFF3, NRHFK, NVIRC
      INTEGER KZETA1, KZETA2, KT1AMP, KT2AMP, IINT2, ISTARTPQ
      INTEGER MT1AM, MT2AM, MT2SQ, M3ORHF, ILASTL, ILASTR
      INTEGER ILSTSYM, MGLMRH, ISTARTZ0, ISTARTPQMO, IDUMMY
      INTEGER MGLMDT, KLAMPA, KLAMHA, KZDENS

      CALL QENTER('CCFPREINT')

*---------------------------------------------------------------------*
* set some symmetries & dimensions and do some initializations:
*---------------------------------------------------------------------*
      IF (CCS .AND. IOPTRSP.EQ.1) THEN
         CALL QEXIT('CCFPREINT')
         RETURN
      END IF

      LENX    = 0
      LENY    = 0
      LENMGD  = 0
      LENZDN  = 0
      MT1AM   = 0
      MT2AM   = 0
      MT2SQ   = 0
      MGLMRH  = 0
      MGLMDT  = 0
      M3ORHF  = 0
      DO ISYM = 1, NSYM
        LENX    = MAX(LENX  ,NMATIJ(ISYM))
        LENY    = MAX(LENY  ,NMATAB(ISYM))
        LENMGD  = MAX(LENMGD,NT2AOIJ(ISYM))
        LENZDN  = MAX(LENZDN,N2BST(ISYM))
        MT1AM   = MAX(MT1AM, NT1AM(ISYM))
        MT2AM   = MAX(MT2AM, NT2AM(ISYM))
        MT2SQ   = MAX(MT2SQ, NT2SQ(ISYM))
        MGLMRH  = MAX(MGLMRH,NGLMRH(ISYM))
        MGLMDT  = MAX(MGLMDT,NGLMDT(ISYM))
        M3ORHF  = MAX(M3ORHF,N3ORHF(ISYM))
      END DO
      LENM    = M3ORHF

      ISTARTZ0   = 1
      ISTARTPQ   = 1
      ISTARTPQMO = 1

*---------------------------------------------------------------------*
* open direct access files:
*---------------------------------------------------------------------*

      LUPQ   = -1
      LUZDN  = -1
      LUXIM  = -1
      LUYIM  = -1
      LUMGD  = -1
      LUBDZ0 = -1
      LUPQMO = -1
      LUMIM  = -1
C
      IF (CCSD .OR. CCSDT) THEN
        CALL WOPEN2(LUPQ, FILPQ, 64,0)
      END IF
      
      IF      (IOPTRSP.EQ.2 .AND. (CCS.OR.CC2)) THEN
        CALL WOPEN2(LUZDN,FILZDN,64,0)
      END IF

      IF (IOPTRSP.EQ.2 .AND. (.NOT.CCS)) THEN
        CALL WOPEN2(LUXIM,FILXIM,64,0)
        CALL WOPEN2(LUYIM,FILYIM,64,0)
        IF (.NOT.CC2) CALL WOPEN2(LUMGD,FILMGD,64,0)
      END IF

      IF (IOPTRSP.EQ.1 .AND. (.NOT.CCS) .AND. (.NOT.CC2)) THEN
        CALL WOPEN2(LUBDZ0,FNBDZ0, 64,0)
        CALL WOPEN2(LUPQMO,FILPQMO,64,0)
      END IF

      IF (IOPTRSP.EQ.2 .AND. (.NOT.CCS) .AND. (.NOT.CC2) 
     &     .AND. CCR12) THEN
        CALL WOPEN2(LUMIM,FILMIM,64,0)
      END IF

*---------------------------------------------------------------------*
* allocate work space:
*---------------------------------------------------------------------*
      KZETA1 = 1
      KT1AMP = KZETA1 + MT1AM
      KEND1  = KT1AMP + MT1AM

      IF (CCS.OR.CC2) THEN
         KLAMHA = KEND1
         KLAMPA = KLAMHA + MGLMDT
         KZDENS = KLAMPA + MGLMDT
         KEND1  = KZDENS + LENZDN
      END IF

      IF (.NOT.CCS) THEN
         KZETA2 = KEND1
         KT2AMP = KZETA2 + MT2SQ
         KXINT  = KT2AMP + MT2AM
         KYINT  = KXINT  + LENX
         KEND1  = KYINT  + LENY
      END IF

      IF (CCSD .OR. CCSDT) THEN
         KCHI   = KEND1
         KMINT  = KCHI  + MAX(LENX,MGLMRH)
         KMGDL  = KMINT + M3ORHF
         KEND1  = KMGDL + LENMGD 
      END IF

      LWRK1 = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CCFPRE2INT.')
      END IF

*---------------------------------------------------------------------*
* loop over pairs of left and right vectors:
*---------------------------------------------------------------------*
      LISTLOLD = '   '
      LISTROLD = '   '
      ILASTL = 0
      ILASTR = 0
      DO IINT2 = 1, NINT2
         LISTL  = VTABLE(INTMED2(2,IINT2))
         IDLSTL = INTMED2(1,IINT2)
         LISTR  = VTABLE(INTMED2(4,IINT2))
         IDLSTR = INTMED2(3,IINT2)
         ISYML  = ILSTSYM(LISTL,IDLSTL)
         ISYMR  = ILSTSYM(LISTR,IDLSTR)
         ISYINT = MULD2H(ISYML,ISYMR)

         NEWLVEC = (IINT2.EQ.1) .OR. (LISTL.NE.LISTLOLD) 
     &                          .OR. (IDLSTL.NE.ILASTL)

         NEWRVEC = (IINT2.EQ.1) .OR. (LISTR.NE.LISTROLD) 
     &                          .OR. (IDLSTR.NE.ILASTR)

*---------------------------------------------------------------------*
* if needed read new lagrangian multiplier and/or amplitude vector:
* Note that if new lagrangian multipliers are read, the amplitudes
* are overwritten and must then also be read from file.
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF (NEWLVEC) THEN
         IOPT = 1
         IF ( (.NOT.(CCS.OR.CC2)) .OR. (CC2.AND.NONHF) ) IOPT = 3
         CALL CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     &                 WORK(KZETA1),WORK(KT2AMP))
         IF (IOPT.EQ.3) CALL CC_T2SQ(WORK(KT2AMP),WORK(KZETA2),ISYML)

         LISTLOLD = LISTL
         ILASTL   = IDLSTL
      END IF

      IF (NEWRVEC .OR. NEWLVEC) THEN
         IOPT = 1
         IF ( (.NOT.(CCS.OR.CC2)) .OR. (CC2.AND.NONHF) ) IOPT = 3
         CALL CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,
     &                 WORK(KT1AMP),WORK(KT2AMP))
         IF ( (IOPT.EQ.3) .AND. (.NOT. LISTR(1:2).EQ.'R0') ) THEN
            CALL CCLR_DIASCL(WORK(KT2AMP),TWO,ISYMR)
         END IF

         LISTROLD = LISTR
         ILASTR   = IDLSTR
      END IF

      TIMIO = TIMIO + SECOND() - DTIME

*---------------------------------------------------------------------*
*     for CCS & CC2  calculate the double AO backtransformed Zeta1 
*     vector (the `Zeta' density):
*---------------------------------------------------------------------*
      IF ((CCS.OR.CC2) .AND. IOPTRSP.EQ.2) THEN
     
         CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPA),XLAMDH,WORK(KLAMHA),
     &                    WORK(KT1AMP),ISYMR)

         CALL DZERO(WORK(KZDENS),LENZDN)

         IOPT = 2
         CALL CC_ZETADEN(WORK(KZDENS),WORK(KZETA1),ISYML,
     &                   XLAMDP,XLAMDH,ISYM0,
     &                   WORK(KLAMPA),WORK(KLAMHA),ISYMR,
     &                   IOPT,WORK(KEND1),LWRK1)

         CALL CC_WVEC(LUZDN,FILZDN,LENZDN,N2BST(ISYINT),IINT2,
     &                WORK(KZDENS))

      END IF
*---------------------------------------------------------------------*
*     calculate X, Y and (MO) Chi intermediate:
*     save X and Y on file, Chi is only needed for BFZeta density
*---------------------------------------------------------------------*
      DTIME = SECOND()

      IF ( (.NOT.(CCS.OR.CC2)) .OR. (CC2.AND.NONHF) ) THEN
        CALL CC_XI(WORK(KXINT),WORK(KZETA2),ISYML,WORK(KT2AMP),ISYMR,
     &             WORK(KEND1),LWRK1)

        CALL CC_YI(WORK(KYINT),WORK(KZETA2),ISYML,WORK(KT2AMP),ISYMR,
     &             WORK(KEND1),LWRK1)
      END IF
      
      IF (CCSD .OR. CCSDT) THEN

         IF (IOPTRSP.EQ.1) THEN

           CALL CCLT_CHI(WORK(KZETA1),WORK(KXINT),ISYML,
     &                   XLAMDP,ISYM0,WORK(KCHI),1)

         ELSE IF (IOPTRSP.EQ.2) THEN

           DO ISYMK = 1, NSYM
              ISYMC = MULD2H(ISYMK,ISYMR)
              ISYMI = MULD2H(ISYMC,ISYML)

              KOFF1 = KT1AMP + IT1AM(ISYMC,ISYMK) 
              KOFF2 = KZETA1 + IT1AM(ISYMC,ISYMI)
              KOFF3 = KCHI   + IMATIJ(ISYMK,ISYMI)

              NRHFK = MAX(NRHF(ISYMK),1)
              NVIRC = MAX(NVIR(ISYMC),1)

              CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     *                   -ONE,WORK(KOFF1),NVIRC,WORK(KOFF2),NVIRC,
     *                   ZERO,WORK(KOFF3),NRHFK)

           END DO
           CALL DAXPY(NMATIJ(ISYINT),-ONE,WORK(KXINT),1,WORK(KCHI),1)

         END IF

      END IF

      TIMXYM = TIMXYM + SECOND() - DTIME

      IF (IOPTRSP.EQ.2 .AND. (.NOT.CCS)) THEN
       DTIME = SECOND()
       CALL CC_WVEC(LUXIM,FILXIM,LENX,NMATIJ(ISYINT),IINT2,WORK(KXINT))
       CALL CC_WVEC(LUYIM,FILYIM,LENY,NMATAB(ISYINT),IINT2,WORK(KYINT))
       TIMIO = TIMIO + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
*     for CCSD calculate also the M intermediates
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
         DTIME = SECOND()
         CALL CC_MI(WORK(KMINT),WORK(KZETA2),ISYML,WORK(KT2AMP),ISYMR,
     &              WORK(KEND1),LWRK1)
         TIMXYM = TIMXYM + SECOND() - DTIME

         IF (IOPTRSP.EQ.2 .AND. CCR12) THEN
           DTIME = SECOND()
           CALL CC_WVEC(LUMIM,FILMIM,LENM,N3ORHF(ISYINT),IINT2,
     &                  WORK(KMINT))
           TIMIO = TIMIO + SECOND() - DTIME
         END IF
      END IF

*---------------------------------------------------------------------*
*     calculate the P and Q intermediates and write them to file:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
        DTIME = SECOND()

        IF (IOPTRSP.EQ.1) THEN

          IOPTY = 1
          CALL CC_PQIMO(WORK(KZETA2),ISYML,WORK(KT2AMP),ISYMR,
     &                  WORK(KYINT),IOPTY,FILPQMO,LUPQMO,
     &                  IADRPQMO,ISTARTPQMO,IINT2,
     &                  WORK(KEND1),LWRK1)

          CALL CC_PQIAO(FILPQMO,LUPQMO,IADRPQMO,ISYINT,
     &                  FILPQ,LUPQ,IADRPQ,ISTARTPQ,IINT2,
     &                  XLAMDH,ISYM0,WORK(KEND1),LWRK1)

        ELSE IF (IOPTRSP.EQ.2) THEN

          IOPTY = 1
          CALL CC_PQI2(WORK(KZETA2),ISYML,WORK(KT2AMP),ISYMR,
     &                 WORK(KYINT),IOPTY,FILPQ,LUPQ,
     &                 IADRPQ,ISTARTPQ,IINT2,WORK(KEND1),LWRK1)

        END IF

        TIMPQ = TIMPQ + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
*     calculate effective density for left BFZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN

         IF (IOPTRSP.EQ.1) THEN

           IOPT = 8
           CALL CC_BFDEN(WORK(KZETA2),ISYML,WORK(KMINT),ISYINT,
     &                   XLAMDP,    ISYM0,  XLAMDP,ISYM0,
     &                   WORK(KCHI),ISYINT, DUMMY, IDUMMY,
     &                   FNBDZ0,LUBDZ0,IADRZ0,ISTARTZ0,
     &                   IINT2,IOPT,.FALSE.,WORK(KEND1),LWRK1)

         ELSE IF (IOPTRSP.EQ.2) THEN

           DTIME = SECOND()
           CALL CC_BFDENF(WORK(KZETA2), ISYML, WORK(KMINT),ISYINT,
     &                    XLAMDP,ISYM0, WORK(KCHI), ISYINT,
     &                    WORK(KT1AMP), ISYMR, WORK(KMGDL),
     &                    WORK(KEND1),LWRK1)
  
           CALL CC_WVEC(LUMGD,FILMGD,LENMGD,NT2AOIJ(ISYINT),
     &                  IINT2,WORK(KMGDL))

         END IF
         TIMBF = TIMBF + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
*     and loop
*---------------------------------------------------------------------*
      END DO

*---------------------------------------------------------------------*
* close direct access files and return:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
        CALL WCLOSE2(LUPQ, FILPQ, 'KEEP')
      END IF

      IF      (IOPTRSP.EQ.2 .AND. (CCS.OR.CC2)) THEN
        CALL WCLOSE2(LUZDN,FILZDN,'KEEP')
      END IF

      IF (IOPTRSP.EQ.2 .AND. (.NOT.CCS)) THEN
        CALL WCLOSE2(LUXIM,FILXIM,'KEEP')
        CALL WCLOSE2(LUYIM,FILYIM,'KEEP')
        IF (.NOT.CC2) CALL WCLOSE2(LUMGD,FILMGD,'KEEP')
      END IF

      IF (IOPTRSP.EQ.1 .AND. (.NOT.CCS) .AND. (.NOT.CC2)) THEN
        CALL WCLOSE2(LUBDZ0,FNBDZ0, 'KEEP')
        CALL WCLOSE2(LUPQMO,FILPQMO,'KEEP')
      END IF

      IF (IOPTRSP.EQ.2 .AND. (.NOT.CCS) .AND. (.NOT.CC2)
     &     .AND. CCR12) THEN
        CALL WCLOSE2(LUMIM,FILMIM,'KEEP')
      END IF

      CALL QEXIT('CCFPREINT')

      RETURN
      END
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCFOPEN*/
*=====================================================================*
      SUBROUTINE CCFOPEN( LUBF,  LUFK,  LUR,  LUFIM,  LUGIM,  LURHO,  
     &                    BFFIL, FKFIL, RFIL, FIMFIL, GIMFIL, RHOFIL, 
     &                    LENBF, LENFK, LENR, LENFIM, LENGIM, LENRHO, 
     &                    NINT1, NINT2, WORK, LWORK                  )
*---------------------------------------------------------------------*
*    Purpose: open files for intermediates in F matrix transformation
*             which are used under the loop over AO integrals and
*             need to be initialized at the beginning.
*      
*     Written by Christof Haettig, November 1998.
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

      CHARACTER*(*) BFFIL, FKFIL, RFIL, FIMFIL, GIMFIL, RHOFIL
      INTEGER       LUBF,  LUFK,  LUR,  LUFIM,  LUGIM,  LURHO
      INTEGER       LENBF, LENFK, LENR, LENFIM, LENGIM, LENRHO
      INTEGER       NINT1, NINT2, LWORK
      INTEGER       IINT1, IINT2, ISYM
      
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION WORK(LWORK)
#endif

      CALL QENTER('CCFOPEN')

*---------------------------------------------------------------------*
* open files for local intermediates:
*---------------------------------------------------------------------*
      LUBF  = -1
      LUR   = -1
      LUFIM = -1
      LUGIM = -1
      LUFK  = -1
      LURHO = -1
      IF (.NOT. (CCS.OR.CC2)) THEN
         CALL WOPEN2(LUBF, BFFIL, 64,0)
         CALL WOPEN2(LUR,  RFIL,  64,0)
         CALL WOPEN2(LUFIM,FIMFIL,64,0)
      ELSE IF (CCS.OR.CC2) THEN
         CALL WOPEN2(LUGIM, GIMFIL,64, 0)
         CALL WOPEN2(LUFK,  FKFIL, 64, 0)
      END IF

      IF (.NOT.CCS) THEN
         CALL WOPEN2(LURHO,RHOFIL,64,0)
      END IF

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
           CALL QUIT('OUT OF MEMORY IN CCFOPEN.')
         END IF
      
         CALL DZERO(WORK,LENBF)

         DO IINT1 = 1, NINT1
           CALL CC_WVEC(LUBF,BFFIL,LENBF,LENBF,IINT1,WORK)
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
           CALL QUIT('OUT OF MEMORY IN CCFOPEN.')
         END IF
      
         CALL DZERO(WORK,LENR)

         DO IINT1 = 1, NINT1
           CALL CC_WVEC(LUR,RFIL,LENR,LENR,IINT1,WORK)
         END DO

      END IF

*---------------------------------------------------------------------*
* calculate a common vector length for the response Fock matrices and
* initialize them with zeros:
*---------------------------------------------------------------------*
      IF (CCS.OR.CC2) THEN
         LENFK = 0
         DO ISYM = 1, NSYM
           LENFK = MAX(LENFK,N2BST(ISYM))
         END DO
         LENGIM = LENFK

         IF (LWORK .LT. LENFK) THEN
           CALL QUIT('OUT OF MEMORY IN CCFOPEN.')
         END IF
   
         CALL DZERO(WORK,LENFK)

         DO IINT1 = 1, NINT1
           CALL CC_WVEC(LUFK,FKFIL,LENFK,LENFK,IINT1,WORK)
         END DO

         DO IINT2 = 1, NINT2
           CALL CC_WVEC(LUGIM,GIMFIL,LENGIM,LENGIM,IINT2,WORK)
         END DO
      END IF

*---------------------------------------------------------------------*
* calculate a common vector length for the RHO1 vectors and
* initialize them with zeros:
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN
         LENRHO = 0
         DO ISYM = 1, NSYM
            LENRHO = MAX(LENRHO,NT1AM(ISYM))
         END DO

         IF (LWORK .LT. LENRHO) THEN
           CALL QUIT('OUT OF MEMORY IN CCFOPEN.')
         END IF
      
         CALL DZERO(WORK,LENRHO)

         DO IINT2 = 1, NINT2
           CALL CC_WVEC(LURHO,RHOFIL,LENRHO,LENRHO,IINT2,WORK)
         END DO
      END IF

*---------------------------------------------------------------------*
* calculate a common vector length for the F intermediates and
* initialize them with zeros:
*---------------------------------------------------------------------*
      IF (.NOT. (CCS.OR.CC2)) THEN
         LENFIM = 0
         DO ISYM = 1, NSYM
            LENFIM = MAX(LENFIM,NT1AO(ISYM))
         END DO

         IF (LWORK .LT. LENFIM) THEN
           CALL QUIT('OUT OF MEMORY IN CCFOPEN.')
         END IF
      
         CALL DZERO(WORK,LENFIM)

         DO IINT2 = 1, NINT2
           CALL CC_WVEC(LUFIM,FIMFIL,LENFIM,LENFIM,IINT2,WORK)
         END DO

      END IF

      CALL QEXIT('CCFOPEN')

      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CCFOPEN                          *
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCFSAVE*/
*=====================================================================*
      SUBROUTINE CCFSAVE(IBATCH, I1HGH, I2HGH, INTMED1, INTMED2,
     &                   KRHO2,  LUBF, BFFIL, LENBF,
     &                   KFOCK,  LUFK, FKFIL, LENFK,
     &                   KRIM,   LUR,  RFIL,  LENR,
     &                   KFIM,   LUFIM,FIMFIL,LENFIM,
     &                   KGIM,   LUGIM,GIMFIL,LENGIM,
     &                   KRHO1,  LURHO,RHOFIL,LENRHO,
     &                   NINT1,  NINT2,WORK,  LWORK )
*---------------------------------------------------------------------*
*    Purpose: save intermediates for F matrix transformation on file
*      
*     Written by Christof Haettig, November 1998.
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

      INTEGER LUBF, LUFK, LUR, LUFIM, LUGIM, LURHO
      INTEGER LENBF,LENFK,LENR,LENFIM,LENGIM, LENRHO
      INTEGER NINT1,NINT2,LWORK,IBATCH
      CHARACTER*(*) BFFIL, FIMFIL, GIMFIL, FKFIL, RFIL, RHOFIL
      INTEGER I1HGH(0:IBATCH), I2HGH(0:IBATCH)
      INTEGER INTMED1(2,NINT1), INTMED2(4,NINT2)
      INTEGER KRHO2(NINT1), KFOCK(NINT1), KRIM(NINT1)
      INTEGER KFIM(NINT2), KRHO1(NINT2), KGIM(NINT2)

      CHARACTER*(3) LIST, LISTL, LISTR
      INTEGER IDLST, IDLSTL, IDLSTR, ISYM, ISYML, ISYMR, ISYRES, LEN
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

      CALL QENTER('CCFSAVE')

*---------------------------------------------------------------------*
* Fock, BF, and R intermediates:
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
            WRITE (LUPRI,*) 'Norm of saved BF        nb. ',IINT1,' is ',
     &           XNORM
          END IF
        END IF

* R intermediate:
        IF (.NOT. (CCS .OR. CC2)) THEN
          LEN = NEMAT1(ISYM)
          CALL CC_WVEC(LUR,RFIL,LENR,LEN,IINT1,WORK(KRIM(IINT1)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KRIM(IINT1)),1,WORK(KRIM(IINT1)),1)
            WRITE (LUPRI,*) 'Norm of saved RIM       nb. ',IINT1,' is ',
     &           XNORM
          END IF
        END IF

* Fock intermediate:
        IF (CCS .OR. CC2) THEN
          LEN = N2BST(ISYM)
          CALL  CC_WVEC (LUFK,FKFIL,LENFK,LEN,IINT1,WORK(KFOCK(IINT1)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KFOCK(IINT1)),1,WORK(KFOCK(IINT1)),1)
            WRITE (LUPRI,*) 'Norm of saved FOCK mat. nb. ',IINT1,' is ',
     &           XNORM
          END IF
        END IF

      END DO

*---------------------------------------------------------------------*
* RHO1, F and the CCS G intermediate:
*---------------------------------------------------------------------*
      DO IINT2 = I2HGH(IBATCH-1)+1, I2HGH(IBATCH)
        LISTL  = VTABLE(INTMED2(2,IINT2))
        LISTR  = VTABLE(INTMED2(4,IINT2))
        IDLSTL = INTMED2(1,IINT2)
        IDLSTR = INTMED2(3,IINT2)
        ISYML  = ILSTSYM(LISTL,IDLSTL)
        ISYMR  = ILSTSYM(LISTR,IDLSTR)
        ISYRES = MULD2H(ISYML,ISYMR)


* Rho1:
        IF (.NOT.CCS) THEN
          LEN    = NT1AM(ISYRES)
          CALL CC_WVEC(LURHO,RHOFIL,LENRHO,LEN,IINT2,WORK(KRHO1(IINT2)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KRHO1(IINT2)),1,WORK(KRHO1(IINT2)),1)
            WRITE (LUPRI,*) 'Norm of saved Rho1      nb. ',IINT2,' is ',
     &           XNORM
          END IF
        END IF

* F intermediate:
        IF (.NOT.(CCS.OR.CC2)) THEN
          LEN    = NT1AO(ISYRES)
          CALL CC_WVEC(LUFIM,FIMFIL,LENFIM,LEN,IINT2,WORK(KFIM(IINT2)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KFIM(IINT2)),1,WORK(KFIM(IINT2)),1)
            WRITE (LUPRI,*) 'Norm of saved FIM       nb. ',IINT2,' is ',
     &           XNORM
          END IF
        END IF

* the CCS/CC2 G intermediate:
        IF (CCS.OR.CC2) THEN
          LEN    = N2BST(ISYRES)
          CALL CC_WVEC(LUGIM,GIMFIL,LENGIM,LEN,IINT2,WORK(KGIM(IINT2)))
          IF (LOCDBG) THEN
            XNORM = DDOT(LEN,WORK(KGIM(IINT2)),1,WORK(KGIM(IINT2)),1)
            WRITE (LUPRI,*) 'Norm of saved GIM       nb. ',IINT2,' is ',
     &           XNORM
          END IF
        END IF

      END DO


      CALL QEXIT('CCFSAVE')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCFSAVE
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCFPRE2 */
*=====================================================================*
      SUBROUTINE CCFPRE2(INTMED2,ISTART,IEND,
     &                   KFIM, LUFIM,FIMFIL,LENFIM,
     &                   KGIM, LUGIM,GIMFIL,LENGIM,
     &                   KRHO1,LURHO,RHOFIL,LENRHO,
     &                   KMGDL,LUMGD,MGDFIL,LENMGD,
     &                   KZDEN,LUZDN,FILZDN,LENZDN,
     &                   KCTR2,KLAMPA,KLAMHA,XLAMDP,XLAMDH,
     &                   WORK,LWORK,KENDIN,KENDOUT )
*---------------------------------------------------------------------*
*    Purpose: prepare for calculation of intermediates that depend
*             on the AO integrals and left and right vectors
*
*    N.B.: this routine allocates work space for CC_FMAT 
*             INPUT  end of used space: KENDIN
*             OUTPUT end of used space: KENDOUT
*      
*     Written by Christof Haettig, November 1998.
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
      INTEGER LUFIM,LENFIM, LURHO,LENRHO, LUMGD,LENMGD
      INTEGER LUZDN,LENZDN, LUGIM,LENGIM
      INTEGER INTMED2(4,IEND)
      INTEGER KLAMPA(IEND), KLAMHA(IEND), KZDEN(IEND), KCTR2(IEND)
      INTEGER KFIM(IEND), KRHO1(IEND), KMGDL(IEND), KGIM(IEND)

      CHARACTER*(*) FIMFIL, RHOFIL, MGDFIL, FILZDN, GIMFIL
      CHARACTER*(3) LISTL, LISTR
      CHARACTER*(10) MODEL
      INTEGER KT1AMPA, KZETA2
      INTEGER LEN, KEND1, IOPT
      INTEGER IDLSTL, IDLSTR, ISYML, ISYMR, ISYRES, IINT
      
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL XLAMDP(NLAMDT), XLAMDH(NLAMDT), ddot
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XLAMDP(NLAMDT), XLAMDH(NLAMDT), ddot
#endif

* external functions:
      INTEGER ILSTSYM

      CALL QENTER('CCFPRE2')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      KENDOUT = KENDIN
           
      DO IINT = ISTART, IEND
        LISTL  = VTABLE(INTMED2(2,IINT))
        LISTR  = VTABLE(INTMED2(4,IINT))
        IDLSTL = INTMED2(1,IINT)
        IDLSTR = INTMED2(3,IINT)
        ISYML  = ILSTSYM(LISTL,IDLSTL)
        ISYMR  = ILSTSYM(LISTR,IDLSTR)
        ISYRES = MULD2H(ISYML,ISYMR)

        IF (CCS) THEN
          KCTR2(IINT)  = KDUM
          KLAMPA(IINT) = KDUM
          KLAMHA(IINT) = KDUM
          KRHO1(IINT)  = KDUM
          KMGDL(IINT)  = KDUM
          KFIM(IINT)   = KDUM
          KGIM(IINT)   = KENDOUT
          KZDEN(IINT)  = KGIM(IINT)     + N2BST(ISYRES)
          KENDOUT      = KZDEN(IINT)    + N2BST(ISYRES)
        ELSE IF (CC2) THEN
          KMGDL(IINT)  = KDUM
          KFIM(IINT)   = KDUM
          KLAMPA(IINT) = KENDOUT
          KLAMHA(IINT) = KLAMPA(IINT)   + NGLMDT(ISYMR)
          KGIM(IINT)   = KLAMHA(IINT)   + NGLMDT(ISYMR)
          KZDEN(IINT)  = KGIM(IINT)     + N2BST(ISYRES)
          KRHO1(IINT)  = KZDEN(IINT)    + N2BST(ISYRES)
          KCTR2(IINT)  = KRHO1(IINT)    + NT1AM(ISYRES)
          KENDOUT      = KCTR2(IINT)    + NT2SQ(ISYML)
        ELSE 
          KCTR2(IINT)  = KDUM
          KGIM(IINT)   = KDUM
          KZDEN(IINT)  = KDUM
          KLAMPA(IINT) = KENDOUT
          KLAMHA(IINT) = KLAMPA(IINT)   + NGLMDT(ISYMR)
          KRHO1(IINT)  = KLAMHA(IINT)   + NGLMDT(ISYMR)
          KMGDL(IINT)  = KRHO1(IINT)    + NT1AM(ISYRES)
          KFIM(IINT)   = KMGDL(IINT)    + NT2AOIJ(ISYRES)
          KENDOUT      = KFIM(IINT)     + NT1AO(ISYRES)
        END IF

        KT1AMPA = KENDOUT 
        KEND1   = KT1AMPA + NT1AM(ISYMR)

        IF (CC2) THEN
          KZETA2 = KEND1
          KEND1  = KZETA2 + NT2AM(ISYML)
        END IF

        IF ( (LWORK-KEND1) .LE. 0 ) THEN
          CALL QUIT('Insufficient work space in CCFPRE2.')
        END IF


*       recover Zeta density for CCS/CC2:
        IF (CCS.OR.CC2) THEN
          LEN = N2BST(ISYRES)
          CALL CC_RVEC(LUZDN,FILZDN,LENZDN,LEN,IINT,WORK(KZDEN(IINT)))
        END IF
          
*       recover G intermediate for CCS/CC2:
        IF (CCS.OR.CC2) THEN
          LEN = N2BST(ISYRES)
          CALL CC_RVEC(LUGIM,GIMFIL,LENGIM,LEN,IINT,WORK(KGIM(IINT)))
        END IF

       
*       recover RHO1:
        IF (.NOT.CCS) THEN
          LEN = NT1AM(ISYRES)
          CALL CC_RVEC(LURHO,RHOFIL,LENRHO,LEN,IINT,WORK(KRHO1(IINT)))
        END IF

*       read singles response vector and calculate response Lambda:
        IF (.NOT.CCS) THEN
          IOPT = 1 ! read singles response vector
          CALL CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,
     &                  WORK(KT1AMPA),WORK(KDUM) )

          ! calculate response Lambda matrices:
          CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPA(IINT)),
     &                     XLAMDH,WORK(KLAMHA(IINT)),
     &                     WORK(KT1AMPA),ISYMR)
        END IF

*       read doubles Zeta response vector:
        IF (CC2) THEN
          IOPT = 2 
          CALL CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     &                  WORK(KDUM),WORK(KZETA2))
          CALL CC_T2SQ(WORK(KZETA2),WORK(KCTR2(IINT)),ISYML)
        END IF

*       recover FIM intermediate:
        IF (CCSD .OR. CCSDT) THEN
          LEN = NT1AO(ISYRES)
          CALL CC_RVEC(LUFIM,FIMFIL,LENFIM,LEN,IINT,WORK(KFIM(IINT)))
        END IF

*       recover XMGDL density:
        IF (CCSD .OR. CCSDT) THEN
          LEN = NT2AOIJ(ISYRES)
          CALL CC_RVEC(LUMGD,MGDFIL,LENMGD,LEN,IINT,WORK(KMGDL(IINT)))
        END IF

      END DO


      CALL QEXIT('CCFPRE2')

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CCFPRE2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCFINT2 */
*=====================================================================*
       SUBROUTINE CCFINT2(IDEL,   ISYDEL,  ISYCTR,
     &                    XINT,   BFRHF,   X3OINT,  XMGDL, ZDEN,  
     &                    FIM,    GIM,     RHO1,
     &                    LUBFI,  FNBFI,   IADRBFI, ISTARTBFI, 
     &                    LUPQR0, FNPQR0,  IADRPQ0,
     &                    LUPQRR, FNPQRR,  IADRPQR,
     &                    LUBDZ0, FNBDZ0,  IADRZ0,
     &                    XLAMPA, XLAMHA,  ISYMA,
     &                    XLAMP0, XLAMH0,  WORK,   LWORK,
     &                    TIMFCK, TIMBF, TIMGZ, TIMIO, TIMFG      )
*---------------------------------------------------------------------*
*
*    Purpose: calculate intermediates for F matrix transformation
*             which require AO integrals and depend on both the
*             response amplitude and the Lagrangian multiplier vectors
*
*             BFZeta   --   written to file 
*             FIM      --   intermediate for 21F term
*             GIM      --   intermediate for CCS case
*             RHO1     --   contribution to  21G term
*            (GZIM     --   intermediate for E   term, added to FIM)
*               
*     input:  XINT           - AO integral distribution
*             BFRHF          - (**|kdel) integrals sorted for CC_BFIF
*             X3OINT         - (ij|kdel) integrals for CC_21H
*             XMGDL          - effective density for CC_BFIF
*             ZDEN           - effective density for GIM intermediate
*             XLAMPA, XLAMHA - response A Lambda matrices
*             XLAMP0, XLAMH0 - ordinary zeroth order Lambda matrices
*
*     Written by Christof Haettig, November 1998.
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
#include "ccorb.h"
#include "second.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
      
      CHARACTER*(*) FNBFI, FNPQR0, FNPQRR, FNBDZ0
      INTEGER LUBFI, LUPQR0, LUPQRR, LUBDZ0, ISTARTBFI
      INTEGER IDEL, ISYDEL, ISYMA, LWORK, ISYRES
      INTEGER ISYCTR, ISYMGD, ISYABK, ISYZT0, ISYZT1, LEN0, LEN1
      INTEGER KDSRHA, KMGD, KEND1, LWRK1, IADR, IOPT
      INTEGER KPINT0, KQINT0, KPINT1, KQINT1, KWINT1
      INTEGER IADRPQ0(MXCORB_CC),IADRPQR(MXCORB_CC)
      INTEGER IADRBFI(MXCORB_CC),IADRZ0(MXCORB_CC)


#if defined (SYS_CRAY)
      REAL WORK(LWORK), XINT(*), BFRHF(*), X3OINT(*)
      REAL XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)
      REAL XMGDL(*), ZDEN(*), FIM(*), GIM(*), RHO1(*)
      REAL TIMBF, TIMGZ, TIMIO, TIMFG, TIMFCK
      REAL DUMMY, DTIME
#else
      DOUBLE PRECISION WORK(LWORK), XINT(*), BFRHF(*), X3OINT(*)
      DOUBLE PRECISION XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)
      DOUBLE PRECISION XMGDL(*), ZDEN(*), FIM(*), GIM(*), RHO1(*)
      DOUBLE PRECISION TIMBF, TIMGZ, TIMIO, TIMFG, TIMFCK
      DOUBLE PRECISION DUMMY, DTIME
#endif


      CALL QENTER('CCFINT2')

  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYRES = MULD2H(ISYCTR,ISYMA)

*---------------------------------------------------------------------*
* calculate the CCS/CC2 GZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCS.OR.CC2) THEN
         DTIME = SECOND()
         CALL CC_AOFOCK(XINT,ZDEN,GIM,WORK,LWORK,IDEL,ISYDEL,.FALSE.,
     &                  DUMMY,ISYRES)
         TIMFCK = TIMFCK + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
* calculate BZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
         DTIME = SECOND()
         CALL CC_BFIF(BFRHF,ISYDEL,XMGDL,ISYRES,LUBFI,FNBFI,
     &                IADRBFI,ISTARTBFI,IDEL,WORK,LWORK)
         TIMBF = TIMBF + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
* calculate the GZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
         DTIME = SECOND()

         ISYABK = MULD2H(ISYDEL,ISYMA)
         ISYMGD = MULD2H(ISYDEL,ISYCTR)

         KDSRHA = 1
         KMGD   = KDSRHA + NDSRHF(ISYABK)
         KEND1  = KMGD   + NT2BGD(ISYMGD)
         LWRK1  = LWORK  - KEND1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in CCFINT2. (21F & 21G)')
         END IF
        
         CALL CCTRBT(XINT,WORK(KDSRHA),XLAMHA,ISYMA,
     &               WORK(KEND1),LWRK1,ISYDEL)

         IADR = IADRZ0(IDEL)
         CALL GETWA2(LUBDZ0,FNBDZ0,WORK(KMGD),IADR,NT2BGD(ISYMGD))

         IOPT = 0
         CALL CC_GIM(WORK(KDSRHA),ISYABK,WORK(KMGD),ISYMGD,FIM,IOPT,
     &               WORK(KEND1),LWRK1)
      
         TIMGZ = TIMGZ + SECOND() - DTIME
      END IF

*---------------------------------------------------------------------*
* calculate 21F term and one part of the 21G term:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN

         ISYZT0 = MULD2H(ISYDEL,ISYCTR)
         ISYZT1 = MULD2H(ISYZT0,ISYMA)
         LEN0   = NT2BCD(ISYZT0)
         LEN1   = NT2BCD(ISYZT1)

         KPINT0 = 1
         KQINT0 = KPINT0 + LEN0
         KPINT1 = KQINT0 + LEN0
         KQINT1 = KPINT1 + LEN1
         KWINT1 = KQINT1 + LEN1
         KEND1  = KWINT1 + LEN1
         LWRK1  = LWORK  - KEND1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in CCFINT2. (21F & 21G)')
         END IF
        
         DTIME = SECOND()

         IADR = IADRPQ0(IDEL)
         CALL GETWA2(LUPQR0,FNPQR0,WORK(KPINT0),IADR,LEN0)
         IADR = IADRPQ0(IDEL) + LEN0
         CALL GETWA2(LUPQR0,FNPQR0,WORK(KQINT0),IADR,LEN0)
      
         IADR = IADRPQR(IDEL)
         CALL GETWA2(LUPQRR,FNPQRR,WORK(KPINT1),IADR,LEN1)
         IADR = IADRPQR(IDEL) + LEN1
         CALL GETWA2(LUPQRR,FNPQRR,WORK(KQINT1),IADR,LEN1)

         TIMIO = TIMIO + SECOND() - DTIME

         IOPT = 2
         DTIME = SECOND()
         CALL CC_21I2(FIM,XINT,ISYDEL,DUMMY,0,
     *                WORK(KPINT0),WORK(KQINT0),ISYZT0,
     *                WORK(KPINT1),WORK(KQINT1),ISYZT1,
     *                XLAMP0,XLAMH0,ISYM0,XLAMPA,ISYMA,
     *                WORK(KEND1),LWRK1,IOPT,
     *                .TRUE.,.FALSE.,.FALSE.)
         TIMFG = TIMFG + SECOND() - DTIME 

         CALL DZERO(WORK(KWINT1),LEN1)

         DTIME  = SECOND()
         CALL CC_21H(RHO1,ISYRES, WORK(KQINT1),
     *               WORK(KWINT1),WORK(KPINT1),ISYZT1,
     *               X3OINT,ISYM0,
     *               WORK(KEND1),LWRK1,ISYDEL)
         TIMFG = TIMFG + SECOND() - DTIME

      END IF
      
*---------------------------------------------------------------------*
* that's it, return:
*---------------------------------------------------------------------*

      CALL QEXIT('CCFINT2')
C
      RETURN

      END
*=====================================================================*
*            END OF SUBROUTINE CCFINT2
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CCFMAT2 */
*=====================================================================*
       SUBROUTINE CCFMAT2(THETA1, THETA2,  ISYRES,
     &                    LISTL,  IDLSTL,  ZETA1,  ISYCTR,  
     &                    LISTR,  IDLSTR,  T1AMPA, ISYAMP,   
     &                    XLAMPA, XLAMHA,  XLAMP0, XLAMH0,  
     &                    XINT,   YINT,    XBARA,  YBARA, 
     &                    RIMA,   A2IM,    FOCKA,  FOCK0AO, FCKC0,
     &                    LUBFI,  FNBFI,   IADRBFI,
     &                    LUBF,   BFFIL,   LENBF,  IINTR,
     &                    LUC,    CTFIL,   LUCBAR, CBAFIL,
     &                    LUD,    DTFIL,   LUDBAR, DBAFIL,
     &                    IOFFCD, WORK,    LWORK                    )
*---------------------------------------------------------------------*
*
*    Purpose: calculate final contributions to F matrix transformation
*             from precalculated intermediates
*
*     Written by Christof Haettig, November 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "dummy.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsections.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccfield.h"
#include "ccslvinf.h"
#include "qm3.h"
!#include "qmmm.h"

* local parameters:
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CCFMAT2> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0, LUBF0
      PARAMETER (ISYM0 = 1)
      
      CHARACTER*(*) LISTL, LISTR, FNBFI, BFFIL
      CHARACTER*(*) CTFIL, CBAFIL, DTFIL, DBAFIL
      INTEGER IDLSTL, IDLSTR, ISYRES, ISYCTR, ISYAMP, LWORK
      INTEGER LUBFI, IADRBFI(MXCORB_CC), IOFFCD
      INTEGER LUC, LUCBAR, LUD, LUDBAR, LUBF, LENBF, IINTR
      INTEGER IOPTTCME
C
C
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL THETA1(*), THETA2(*), ZETA1(*), T1AMPA(*)
      REAL XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)
      REAL XINT(*), YINT(*), XBARA(*), YBARA(*)
      REAL RIMA(*), A2IM(*), FOCKA(*), FOCK0AO(*), FCKC0(*)
      REAL XNORM, DDOT, ZERO, ONE, HALF, FACB, FF
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION THETA1(*), THETA2(*), ZETA1(*), T1AMPA(*)
      DOUBLE PRECISION XLAMP0(*), XLAMH0(*), XLAMPA(*), XLAMHA(*)
      DOUBLE PRECISION XINT(*), YINT(*), XBARA(*), YBARA(*), RIMA(*)
      DOUBLE PRECISION A2IM(*), FOCKA(*), FOCK0AO(*), FCKC0(*)
      DOUBLE PRECISION XNORM, DDOT, ZERO, ONE, HALF, FACB, FF
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0)

      LOGICAL LGAMMA, LBFZETA, LO3BF, LRCON, LGCON, FCKCON
      CHARACTER*(10) MODEL
      INTEGER KGZETA, KGIMA, KEND0, LWRK0, ISYMB, ISYMJ, KFOCKA
      INTEGER IOPTG, ICON, KEND1, LWRK1, KBZETA, KFOCK
      INTEGER IOPT, KGAMMA, KBF, KZETA2, KOFF1, KOFF2, KLIAJB
      INTEGER KYPS, KCHI, ISYMA, ISYMD, ISYMK, NVIRA, NVIRD, KOFF3
      INTEGER IOPTE, IOPTB, KCDBAR, KEMAT1, KEMAT2
      INTEGER KONEHR, KONEHG, KONEH, KA2CON
      INTEGER NRHFK, NVIRC, ISYMC, ISYMI, IF, IOPTR12
      INTEGER KFIELD, KFIELDR, KFIELDG, ISYLMA

  
      CALL QENTER('CCFMAT2')

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (CCS) THEN
         KONEHR  = 1
         KONEHG  = KONEHR + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEMAT1  = KONEHG + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEMAT2  = KEMAT1 + NMATAB(ISYAMP)
         KEND0   = KEMAT2 + NMATIJ(ISYAMP)
      ELSE 
         KGZETA  = 1
         KGIMA   = KGZETA + NT1AO(ISYRES)
         KFOCKA  = KGIMA  + NT1AO(ISYAMP)
         KYPS    = KFOCKA + NT1AM(ISYAMP)
         KEMAT1  = KYPS   + NMATAB(ISYRES)
         KEMAT2  = KEMAT1 + MAX(NMATAB(ISYAMP),NEMAT1(ISYAMP))
         KONEH   = KEMAT2 + NMATIJ(ISYAMP)
         KONEHR  = KONEH  + N2BST(ISYM0)
         KONEHG  = KONEHR + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KA2CON  = KONEHG + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEND0   = KA2CON + NT1AM(ISYRES)
      END IF

      IF (CC2.AND.NONHF) THEN
         KFIELD  = KEND0
         KFIELDR = KFIELD  + N2BST(ISYM0)
         KFIELDG = KFIELDR + MAX(N2BST(ISYM0),N2BST(ISYAMP))
         KEND0   = KFIELDG + MAX(N2BST(ISYM0),N2BST(ISYAMP))
      END IF

      LWRK0 = LWORK - KEND0
      IF (LWRK0.LT.0) THEN
         CALL QUIT('Insufficient memory in CCFMAT2. (0)')
      END IF

*---------------------------------------------------------------------*
* for CCS complete here the calculation of the FOCK^A matrix and
* put it into ONEHG and ONEHR arrays used for construction of E1 & E2:
*---------------------------------------------------------------------*
      IF (CCS.OR.CC2) THEN

         CALL DCOPY(N2BST(ISYM0),FOCK0AO,1,WORK(KONEHR),1)
         CALL DCOPY(N2BST(ISYM0),FOCK0AO,1,WORK(KONEHG),1)

         CALL CC_FCKMO(WORK(KONEHR),XLAMPA,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYAMP,ISYM0)
         
         CALL CC_FCKMO(WORK(KONEHG),XLAMP0,XLAMHA,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYAMP)

         CALL DAXPY(N2BST(ISYAMP),ONE,WORK(KONEHG),1,WORK(KONEHR),1)
         CALL DAXPY(N2BST(ISYAMP),ONE,FOCKA,       1,WORK(KONEHR),1)

         CALL DCOPY(N2BST(ISYAMP),WORK(KONEHR),1,WORK(KONEHG),1)
         
      END IF

*---------------------------------------------------------------------*
* get one-electron hamiltonian and transform to MO:
* if frozen core used, add core-only zeroth-order Fock matrix
*---------------------------------------------------------------------*
      IF (.NOT.(CCS.OR.CC2)) THEN
         CALL CCRHS_ONEAO(WORK(KONEH),WORK(KEND0),LWRK0)
         DO IF = 1, NFIELD
           FF = EFIELD(IF)
           CALL CC_ONEP(WORK(KONEH),WORK(KEND0),LWRK0,FF,1,LFIELD(IF))
         END DO
C
C---------------------------------------------------------------------------
C     Solvent contribution to one-electron integrals. 
C     T^g contribution to transformation.
C    SLV98,OC, NYQMMM10 KS
C------------------------------------------------------------------------------
C
         IF (CCSLV .AND. (.NOT. CCMM)) THEN             
            CALL CCSL_RHSTG(WORK(KONEH),WORK(KEND0),LWRK0)
         ENDIF
C
         IF (CCMM) THEN
            IF (.NOT. NYQMMM) THEN
               CALL CCMM_RHSTG(WORK(KONEH),WORK(KEND0),LWRK0)
            ELSE IF (NYQMMM) THEN
              IF (HFFLD) THEN
                 CALL CCMM_ADDGHF(WORK(KONEH),WORK(KEND0),LWRK0)
              ELSE
                CALL CCMM_ADDG(WORK(KONEH),WORK(KEND0),LWRK0)
              END IF
            END IF
         ENDIF
C-----------------------------------------------------------------

         IF (FROIMP.OR.FROEXP) THEN
           CALL DAXPY(N2BST(ISYM0),ONE,FCKC0,1,WORK(KONEH),1)
         END IF

         CALL DCOPY(N2BST(ISYM0),WORK(KONEH),1,WORK(KONEHR),1)
         CALL DCOPY(N2BST(ISYM0),WORK(KONEH),1,WORK(KONEHG),1)
    
         CALL CC_FCKMO(WORK(KONEH),XLAMP0,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYM0)
         CALL CC_FCKMO(WORK(KONEHR),XLAMPA,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYAMP,ISYM0)
         CALL CC_FCKMO(WORK(KONEHG),XLAMP0,XLAMHA,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYAMP)
      END IF

*---------------------------------------------------------------------*
* get external fields added after the SCF and transform to MO:
*---------------------------------------------------------------------*
      IF (CC2.AND.NONHF) THEN
         CALL DZERO(WORK(KFIELD),N2BST(ISYM0))
         DO IF = 1, NFIELD
           FF = EFIELD(IF)
           CALL CC_ONEP(WORK(KFIELD),WORK(KEND0),LWRK0,FF,1,LFIELD(IF))
         END DO

C
C-----------------------------------------------------------------------------------------
C     Solvent contribution to one-electron integrals. 
C     T^g contribution to transformation.
C     SLV98,OC, NYQMMM10 KS
C------------------------------------------------------------------------------
C
         IF (CCSLV .AND. (.NOT. CCMM)) THEN
            CALL CCSL_RHSTG(WORK(KFIELD),WORK(KEND0),LWRK0)
         ENDIF
C
         IF (CCMM) THEN
            IF (.NOT. NYQMMM) THEN
               CALL CCMM_RHSTG(WORK(KFIELD),WORK(KEND0),LWRK0)
            ELSE IF (NYQMMM) THEN
              IF (HFFLD) THEN
                CALL CCMM_ADDGHF(WORK(KFIELD),WORK(KEND0),LWRK0)
              ELSE
                CALL CCMM_ADDG(WORK(KFIELD),WORK(KEND0),LWRK0)
              END IF
            END IF
         ENDIF
C------------------------------------------------------------------------------
C
         CALL DCOPY(N2BST(ISYM0),WORK(KFIELD),1,WORK(KFIELDR),1)
         CALL DCOPY(N2BST(ISYM0),WORK(KFIELD),1,WORK(KFIELDG),1)
    
         CALL CC_FCKMO(WORK(KFIELD),XLAMP0,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYM0)
         CALL CC_FCKMO(WORK(KFIELDR),XLAMPA,XLAMH0,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYAMP,ISYM0)
         CALL CC_FCKMO(WORK(KFIELDG),XLAMP0,XLAMHA,WORK(KEND0),LWRK0,
     &                 ISYM0,ISYM0,ISYAMP)
      END IF

*---------------------------------------------------------------------*
* calculate BZeta contribution and GZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
         KBZETA = KEND0
         KEND1  = KBZETA + NT2AO(ISYRES)
         LWRK1  = LWORK  - KEND1
      
         IF (LWRK1.LT.0) THEN
            CALL QUIT('Insufficient memory in CCFMAT2.')
         END IF
 
         CALL CC_BFIFSORT(WORK(KBZETA),ISYRES,LUBFI,FNBFI,IADRBFI,
     &                    WORK(KEND1),LWRK1)

         CALL DZERO(WORK(KGZETA),NT1AO(ISYRES))

         ICON    = 1
         IOPTG   = 1
         LGAMMA  = .FALSE.
         LO3BF   = .FALSE.
         LBFZETA = .TRUE.
         CALL CC_T2MO3(DUMMY,DUMMY,1,WORK(KBZETA),
     &                 THETA2,DUMMY,WORK(KGZETA),DUMMY,
     &                 XLAMH0,ISYM0,XLAMH0,ISYM0,
     &                 WORK(KEND1),LWRK1,ISYRES,
     &                 ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
      
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,'norm of THETA1 after B+E term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,'norm of THETA2 after B+E term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C           CALL AROUND('THETA2 after B+E term:')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate the E1' contribution from the GZeta intermediate:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN

         CALL CC_T1AM(THETA1,ISYRES,WORK(KGZETA),ISYRES,
     &                XLAMH0,ISYM0,ONE)

         IF (LOCDBG) THEN
            CALL AROUND("THETA1 after E1' term:")
            CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate D2 contribution to Theta1 results vector: 
*---------------------------------------------------------------------*
      IF (.NOT.(CCS.OR.CC2)) THEN

         CALL CC_21EFM(THETA1,WORK(KONEH),ISYM0,XINT,YINT,ISYRES,
     &                 WORK(KEND0),LWRK0)

      END IF

      IF (LOCDBG .AND. .NOT.CCS) THEN
         CALL AROUND('THETA after singles excit. D2 term:')
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*
* for CC2/NONHF and CCSD read double excitation part of the Lagrangian 
* multipliers and expand to a full square matrix:
*---------------------------------------------------------------------*
      IF ((CC2.AND.NONHF) .OR. CCSD .OR. CCSDT) THEN
        KZETA2 = KEND0
        KEND1  = KZETA2 + NT2SQ(ISYCTR)
        LWRK1  = LWORK  - KEND1
      
        IF (LWRK1.LT.NT2AM(ISYCTR)) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        IOPT = 2
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,DUMMY,WORK(KEND1))
        CALL CC_T2SQ(WORK(KEND1),WORK(KZETA2),ISYCTR)

      ELSE
        KEND1 = KEND0
        LWRK1 = LWORK  - KEND1
      END IF

*---------------------------------------------------------------------*
* calculate the first E2 term contribution from the zero-order BF:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
        KZETA2 = KEND0
        KBF    = KZETA2 + NT2SQ(ISYCTR)
        KEND1  = KBF    + 2 * NT2ORT(ISYM0)
        LWRK1  = LWORK  - KEND1
     
        IF (LWRK1.LT.0) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        LUBF0 = -1
        CALL GPOPEN(LUBF0,'CC_BFIM','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        READ(LUBF0) (WORK(KBF-1+I),I=1,2*NT2ORT(ISYM0))
        CALL GPCLOSE(LUBF0,'KEEP')


        ICON    = 3 
        IOPTG   = 0
        LGAMMA  = .FALSE.
        LO3BF   = .FALSE.
        LBFZETA = .FALSE.
        CALL CC_T2MO3(THETA1,WORK(KZETA2),ISYCTR,WORK(KBF),
     &                DUMMY,DUMMY,DUMMY,DUMMY,
     &                XLAMP0,ISYM0,XLAMPA,ISYAMP,
     &                WORK(KEND1),LWRK1,ISYM0,
     &                ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
     
        IF (LOCDBG) THEN
           WRITE (LUPRI,*) MSGDBG,
     &          'norm(THETA1) after first E2 (LT21C) term:',
     &       DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
           IF (.NOT.CCS) 
     &      WRITE (LUPRI,*) MSGDBG,
     &          'norm(THETA2) after first E2 (LT21C) term:',
     &       DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C          CALL AROUND('THETA after first E2 (LT21C) term:')
C          CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
        END IF

      END IF

*---------------------------------------------------------------------*
* calculate second E2 term contribution and the GAMMA and G 
* intermediates from the response BF intermediate: 
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
        KZETA2 = KEND0
        KBF    = KZETA2 + NT2SQ(ISYCTR)
        KGAMMA = KBF    + NT2AOIJ(ISYAMP)
        KEND1  = KGAMMA + NGAMMA(ISYAMP)
        LWRK1  = LWORK  - KEND1
      
        IF (LWRK1.LT.NT2AM(ISYCTR)) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        CALL CC_RVEC(LUBF,BFFIL,LENBF,NT2AOIJ(ISYAMP),IINTR,WORK(KBF))

        CALL DZERO(WORK(KGAMMA),NGAMMA(ISYAMP))
        CALL DZERO(WORK(KGIMA), NT1AO(ISYAMP))


        ICON    = 3
        IOPTG   = 2
        LGAMMA  = .TRUE.
        LO3BF   = .TRUE.
        LBFZETA = .FALSE.
        CALL CC_T2MO3(THETA1,WORK(KZETA2),ISYCTR,WORK(KBF),
     &                DUMMY,WORK(KGAMMA),WORK(KGIMA),DUMMY,
     &                XLAMP0,ISYM0,XLAMP0,ISYM0,
     &                WORK(KEND1),LWRK1,ISYAMP,
     &                ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
      
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG,
     &          'norm(THETA1) after second E2 (LT21C) term:',
     &      DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
          WRITE (LUPRI,*) MSGDBG,
     &         'norm(THETA2) after second E2 (LT21C) term:',
     &      DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
          WRITE (LUPRI,*) MSGDBG,'the GIMA intermediate:'
          WRITE (LUPRI,'(5F12.6)') (WORK(KGIMA-1+I),I=1,NT1AO(ISYAMP))
C         CALL AROUND('THETA after second E2 (LT21C) term:')
C         CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
          XNORM = DDOT(NGAMMA(ISYAMP),WORK(KGAMMA),1,WORK(KGAMMA),1)
          WRITE (LUPRI,*) 'Norm(GAMMA) = ',XNORM
        END IF

      END IF

*---------------------------------------------------------------------*
* for CCS/CCSD calculate here E1/E1* and E2 intermediates as nedded 
* for the B2 contribution from the G and YBAR intermediates:
*---------------------------------------------------------------------*
      LRCON  = .FALSE.
      LGCON  = .FALSE.
      FCKCON = .TRUE.
      IOPT   = 1

      IF (CCSD .OR. CCSDT) LGCON  = .TRUE.

      CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     &            DUMMY,DUMMY,WORK(KGIMA),DUMMY,
     &            WORK(KONEHR),WORK(KONEHG),
     &            XLAMH0,XLAMP0,ISYM0,DUMMY,DUMMY,IDUMMY,
     &            FCKCON,LRCON,LGCON,.FALSE.,IOPT,ISYAMP)

      IF (CC2) THEN
         CALL DAXPY(NMATAB(ISYAMP),ONE,YBARA,1,WORK(KEMAT1),1)
         CALL DAXPY(NMATIJ(ISYAMP),ONE,XBARA,1,WORK(KEMAT2),1)
      END IF

      IF (CCSD .OR. CCSDT) THEN
         CALL DAXPY(NMATAB(ISYAMP),ONE,YBARA,1,WORK(KEMAT1),1)
      END IF

      IF (LOCDBG) THEN
         CALL AROUND('E-intermediates out from CC_EIM:')
         CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYAMP,1)
      END IF

*---------------------------------------------------------------------*
* calculate B2 contribution to Theta1 results vector: 
*---------------------------------------------------------------------*
      CALL CCLR_E1C1(THETA1,ZETA1,WORK(KEMAT1),WORK(KEMAT2),
     &               WORK(KEND1),LWRK1,ISYCTR,ISYAMP,'T')

      IF (LOCDBG) THEN
         CALL AROUND('THETA1 B2 term:')
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
      END IF

*---------------------------------------------------------------------*
* for CCSD calculate here the E intermed. from G, R and YBAR intermed.:
*---------------------------------------------------------------------*
      IF (.NOT.(CCS.OR.CC2)) THEN

         LRCON  = .TRUE.
         LGCON  = .TRUE.
         FCKCON = .TRUE.
         IOPT   = 1
         CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     &               RIMA,DUMMY,WORK(KGIMA),DUMMY,
     &               WORK(KONEHR),WORK(KONEHG),
     &               XLAMH0,XLAMP0,ISYM0,DUMMY,DUMMY,IDUMMY,
     &               FCKCON,LRCON,LGCON,.FALSE.,IOPT,ISYAMP)

         CALL DAXPY(NMATAB(ISYAMP),ONE,YBARA,1,WORK(KEMAT1),1)

         IF (LOCDBG) THEN
            CALL AROUND('E-intermediates out from CC_EIM:')
            CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYAMP,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate the A term and add to Theta2 result vector:
*---------------------------------------------------------------------*
      IF (CCSD .OR. CCSDT) THEN
         CALL CC_GAMMA2(WORK(KGAMMA),WORK(KEMAT2),ISYAMP)

         IOPT = 2
         CALL CCRHS_A(THETA2,WORK(KZETA2),WORK(KGAMMA),
     &                WORK(KEND1),LWRK1,ISYAMP,ISYCTR,IOPT)
      
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA1) after doubles excit A term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA2) after doubles excit A term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C           CALL AROUND('THETA after doubles excit. A term:')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* read (ia|jb) integrals and calculate L(ia|jb) in place:
* (needed for H term and for A2 term)
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN
         KLIAJB = KEND0
         KEND1  = KLIAJB + NT2AM(ISYM0)
         LWRK1  = LWORK  - KEND1
      
         IF (LWRK1.LT.0) THEN
            CALL QUIT('Insufficient memory in CCFMAT2.')
         END IF
      
         CALL CCG_RDIAJB(WORK(KLIAJB),NT2AM(ISYM0))
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KLIAJB),ONE,ISYM0,IOPTTCME)
      END IF

*---------------------------------------------------------------------*
* calculate A2 contribution:
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN
     
         CALL CCG_LXD(WORK(KA2CON),ISYRES,A2IM,ISYRES,
     &                WORK(KLIAJB),ISYM0,0)

         CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KA2CON),1,THETA1,1)

         IF (LOCDBG) THEN
            CALL AROUND('THETA1 A2 term:')
            CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
            CALL AROUND('A2IM intermediate:')
            CALL CC_PRP(A2IM,THETA2,ISYRES,1,0)
            CALL AROUND('A2CON contribution:')
            CALL CC_PRP(WORK(KA2CON),THETA2,ISYRES,1,0)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate H term:
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN
         KFOCK  = KEND1
         KCHI   = KFOCK  + N2BST(ISYAMP)
         KEND1  = KCHI   + NMATIJ(ISYRES)
         LWRK1  = LWORK  - KEND1
      
         IF (LWRK1.LT.0) THEN
            CALL QUIT('Insufficient memory in CCFMAT2.')
         END IF
      
C        ------------------------
C        calculate F^A_jb matrix:
C        ------------------------
         IOPT = 0
         CALL CCG_LXD(WORK(KFOCKA),ISYAMP,T1AMPA,ISYAMP,
     &                WORK(KLIAJB),ISYM0,IOPT)

         IF (LOCDBG) THEN
            CALL AROUND('FOCKA_jb matrix:')
            CALL CC_PRP(WORK(KFOCKA),WORK,ISYRES,1,0)
         END IF

CTST
         CALL DZERO(WORK(KFOCK),N2BST(ISYAMP))
CTST

C        -----------------------------------------
C        resort Fock^A_jb into a full Fock matrix:
C        -----------------------------------------
         DO ISYMB = 1, NSYM
            ISYMJ = MULD2H(ISYAMP,ISYMB)

            DO J = 1, NRHF(ISYMJ)
            DO B = 1, NVIR(ISYMB)
               KOFF1 = IFCVIR(ISYMJ,ISYMB) + NORB(ISYMJ)*(B-1) + J
               KOFF2 = IT1AM(ISYMB,ISYMJ)  + NVIR(ISYMB)*(J-1) + B
               WORK(KFOCK-1+KOFF1) = WORK(KFOCKA-1+KOFF2)
            END DO
            END DO

         END DO

C        ---------------------------------------
C        calculate the simple Fock contribution:
C        ---------------------------------------
         CALL CC_L1FCK(THETA2,ZETA1,WORK(KFOCK),ISYCTR,ISYAMP,
     &                 WORK(KEND1),LWRK1)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA1) after doubles excit H term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA2) after doubles excit H term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C          CALL AROUND('THETA after doubles exci. H term (Fock part):')
C          CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
           XNORM = DDOT(NT1AM(ISYCTR),ZETA1,1,ZETA1,1)
           WRITE (LUPRI,*) 'Norm of ZETA1:',XNORM
         END IF

C        ----------------------------------------------
C        calculate the transposed Ypsilon intermediate:
C        (for CC2 exclude the contribution from YINT)
C        ----------------------------------------------
         IF (CC2) THEN
            CALL DZERO(WORK(KYPS),NMATAB(ISYRES))
         ELSE
            CALL DCOPY(NMATAB(ISYRES),YINT,1,WORK(KYPS),1)
         END IF

         DO ISYMA = 1, NSYM

            ISYMD = MULD2H(ISYMA,ISYRES)
            ISYMK = MULD2H(ISYMA,ISYCTR)

            KOFF1 = IT1AM(ISYMD,ISYMK) + 1
            KOFF2 = IT1AM(ISYMA,ISYMK) + 1
            KOFF3 = KYPS + IMATAB(ISYMD,ISYMA) 

            NVIRD = MAX(NVIR(ISYMD),1)
            NVIRA = MAX(NVIR(ISYMA),1)

            CALL DGEMM('N','T',NVIR(ISYMD),NVIR(ISYMA),NRHF(ISYMK),
     &                 ONE,T1AMPA(KOFF1),NVIRD,ZETA1(KOFF2),NVIRA,
     &                 ONE,WORK(KOFF3),NVIRD)
         END DO    
      
         CALL DSCAL(NMATAB(ISYRES),HALF,WORK(KYPS),1)

C        ------------------------------------------------------
C        calculate the Chi intermediate:
C        (for CCSD this contribution is included in the B term)
C        ------------------------------------------------------
         CALL DZERO(WORK(KCHI),NMATIJ(ISYRES))

         IF (CC2) THEN

           DO ISYMK = 1, NSYM
              ISYMC = MULD2H(ISYMK,ISYAMP)
              ISYMI = MULD2H(ISYMC,ISYCTR)

              KOFF1 = IT1AM(ISYMC,ISYMK) + 1
              KOFF2 = IT1AM(ISYMC,ISYMI) + 1
              KOFF3 = KCHI + IMATIJ(ISYMK,ISYMI)

              NRHFK = MAX(NRHF(ISYMK),1)
              NVIRC = MAX(NVIR(ISYMC),1)

              CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     *                   ONE,T1AMPA(KOFF1),NVIRC,ZETA1(KOFF2),NVIRC,
     *                   ONE,WORK(KOFF3),NRHFK)

           END DO

           CALL DSCAL(NMATIJ(ISYRES),HALF,WORK(KCHI),1)

         END IF
         
C        --------------------------------
C        contract with (jb|id) integrals:
C        --------------------------------
         CALL CC_22EC(THETA2,WORK(KLIAJB),WORK(KCHI),WORK(KYPS),
     &                ISYRES,WORK(KEND1),LWRK1)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA1) after doubles excit H term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            IF (.NOT.CCS) 
     &       WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA2) after doubles excit H term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
C           CALL AROUND('THETA after doubles excit. H term (2. part):')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
            CALL AROUND('the Ypsilon/Chi intermediates:')
            CALL CC_PREI(WORK(KYPS),WORK(KCHI),ISYRES,1)
            CALL AROUND('the Y and X intermediates:')
            CALL CC_PREI(YINT,XINT,ISYRES,1)
C           CALL AROUND('the integrals:')
C           CALL CC_PRSQ(WORK,WORK(KLIAJB),ISYM0,0,1)
         END IF

      END IF

*---------------------------------------------------------------------*
* calculate the D & C term contributions:
*---------------------------------------------------------------------*
      IF (.NOT. (CCS.OR.CC2)) THEN

        KCDBAR = KEND0
        KZETA2 = KCDBAR  + NT2SQ(ISYAMP)
        KEND1  = KZETA2  + NT2AM(ISYCTR)
        LWRK1  = LWORK   - KEND1

        IF (LWRK1 .LT. NT2AM(ISYCTR)) THEN
           CALL QUIT('Insufficient memory in CCFMAT2.')
        END IF

        IOPT = 2
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                DUMMY,WORK(KZETA2))

C       -----------------------------------------------
C       complete the calculation of the D intermediate:
C       (including the E1 contribution to the diagonal)
C       -----------------------------------------------
        IOPT  = 3
        IOPTB = 1
        IOPTE = 1
        FACB  = ONE
        CALL CCRHS_DIO3(WORK(KCDBAR),DUMMY,XLAMH0,
     &                  WORK(KEND1),LWRK1,1,ISYAMP,
     &                  LUD,DTFIL,LUC,CTFIL,IINTR,IOPT,
     &                  IOPTB,LUDBAR,DBAFIL,IOFFCD,FACB,
     &                  IOPTE,WORK(KEMAT1))

C       -----------------------------------------------
C       contract with the lagrangian multiplier vector:
C       -----------------------------------------------
        ioptr12 = 0
        CALL CC_CD('D',-1,ioptr12,THETA2,ISYRES,WORK(KZETA2),ISYCTR,
     &             WORK(KCDBAR),ISYAMP,WORK(KEND1),LWRK1)

        IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &          'norm(THETA1) after doubles excit D term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA2) after doubles excit D term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
            WRITE (LUPRI,*) MSGDBG,
     &           'norm(D int.) after doubles excit D term:',
     &        DDOT(NT2SQ(ISYAMP),WORK(KCDBAR),1,WORK(KCDBAR),1)
C           CALL AROUND('THETA2 after D term:')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
        END IF

C       -----------------------------------------------
C       complete the calculation of the C intermediate:
C       (including the E1 contribution to the diagonal)
C       -----------------------------------------------
        IOPT  = 3
        IOPTB = 1
        IOPTE = 1
        FACB  = ONE
        CALL CCRHS_CIO3(WORK(KCDBAR),DUMMY,XLAMH0,
     &                  WORK(KEND1),LWRK1,1,ISYAMP,
     &                  LUC,CTFIL,IINTR,IOPT,
     &                  IOPTB,LUCBAR,CBAFIL,IOFFCD,FACB,
     &                  IOPTE,WORK(KEMAT1))
       
C       -----------------------------------------------
C       contract with the lagrangian multiplier vector:
C       -----------------------------------------------
        ioptr12 = 0
        CALL CC_CD('C',-1,ioptr12,THETA2,ISYRES,WORK(KZETA2),ISYCTR,
     &             WORK(KCDBAR),ISYAMP,WORK(KEND1),LWRK1)

        IF (LOCDBG) THEN
            WRITE (LUPRI,*) MSGDBG,
     &          'norm(THETA1) after doubles excit C term:',
     &        DDOT(NT1AM(ISYRES),THETA1,1,THETA1,1)
            WRITE (LUPRI,*) MSGDBG,
     &           'norm(THETA2) after doubles excit C term:',
     &        DDOT(NT2AM(ISYRES),THETA2,1,THETA2,1)
            WRITE (LUPRI,*) MSGDBG,
     &           'norm(C int.) after doubles excit C term:',
     &        DDOT(NT2SQ(ISYAMP),WORK(KCDBAR),1,WORK(KCDBAR),1)
C           CALL AROUND('THETA2 after C term:')
C           CALL CC_PRP(THETA1,THETA2,ISYRES,0,1)
        END IF

      END IF

*---------------------------------------------------------------------*
* that's it, return:
*---------------------------------------------------------------------*

      CALL QEXIT('CCFMAT2')

      RETURN

      END
*=====================================================================*
*            END OF SUBROUTINE CCFMAT2
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_21GMO(THETA1,ISYRES,XILJD,ISYINT,
     &                    LUPQMO,FILPQMO,IADRPQMO,ISYPQI,
     &                    WORK,LWORK)
*---------------------------------------------------------------------*
*
* Purpose: Calculate a 21G like contribution to the F matrix transf.
*          from (il|jd) integrals and P & Q intermediates in MO
*
*     Theta1 = Theta1 - P^d_alj g^X_iljd + 1/2 (P+Q)^d_alj g^X_jlid
*
*     symmetries: ISYRES - THETA1 result vector
*                 ISYINT - (il|jd) integrals
*                 ISYPQI - P and Q intermediates on file
*
* Christof Haettig, November 1998
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"

      CHARACTER*(*) FILPQMO
      INTEGER LWORK, ISYRES, ISYINT, ISYPQI, LUPQMO, IADRPQMO(*)

#if defined (SYS_CRAY)
      REAL XILJD(*), THETA1(*), WORK(LWORK)
#else
      DOUBLE PRECISION XILJD(*), THETA1(*), WORK(LWORK)
#endif

      INTEGER IVIRD, ISYALJ, ISYJLI, ISYMD, IADR, LEN, KOFF1
      INTEGER KPINT, KWINT, KQINT, KEND1, LWRK1
 
      CALL QENTER('CC_21GMO')

*---------------------------------------------------------------------*
*     check symmetries:
*---------------------------------------------------------------------*
      IF (MULD2H(ISYINT,ISYPQI).NE.ISYRES) THEN
         CALL QUIT('Symmetry mismatch in CC_21GMO.')
      END IF

*---------------------------------------------------------------------*
*     contract P and Q intermediates with intgrals in a loop over the
*     virtual index D
*---------------------------------------------------------------------*
      DO ISYMD = 1, NSYM
         ISYALJ = MULD2H(ISYPQI,ISYMD)
         ISYJLI = MULD2H(ISYINT,ISYMD)
 
         LEN   = NT2BCD(ISYALJ)

         KPINT = 1
         KQINT = KPINT + LEN
         KWINT = KQINT + LEN
         KEND1 = KWINT + LEN
         LWRK1 = LWORK - KEND1

         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient memory in CC_21GMO.')
         END IF

         CALL DZERO(WORK(KWINT),LEN)

         DO D = 1, NVIR(ISYMD)
            IVIRD = IVIR(ISYMD) + D
            IADR  = IADRPQMO(IVIRD)
            CALL GETWA2(LUPQMO,FILPQMO,WORK(KPINT),IADR,    LEN)
            CALL GETWA2(LUPQMO,FILPQMO,WORK(KQINT),IADR+LEN,LEN)

            KOFF1 = ISJIKA(ISYJLI,ISYMD) + NMAIJK(ISYJLI)*(D-1) + 1
            
            CALL CC_21H(THETA1,ISYRES,WORK(KQINT),WORK(KWINT),
     &                  WORK(KPINT),ISYALJ,XILJD(KOFF1),ISYINT,
     &                  WORK(KEND1),LWRK1,ISYMD)
         END DO

      END DO

      CALL QEXIT('CC_21GMO')

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_21GMO                          *
*=====================================================================*
*---------------------------------------------------------------------*
c/* Deck CC_ZETADEN */
*=====================================================================*
      SUBROUTINE CC_ZETADEN(ZDEN,   ZETA1,  ISYCTR, 
     &                      XLAMP1, XLAMH1, ISYXL1,
     &                      XLAMP2, XLAMH2, ISYXL2,
     &                      IOPT,   WORK,   LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: calculate `Zeta'-density for F-matrix: a double AO
*             backtransformed Zeta1 vector
*
*    IOPT = 1  -  Jacobian: XLAMP1,XLAMH1 same as XLAMP2/XLAMH2
*    IOPT = 2  -  F matrix: XLAMP1,XLAMH1 different from XLAMP2/XLAMH2
*        
*     Written by Christof Haettig, November 1998
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IOPT, ISYCTR, ISYXL1, ISYXL2, LWORK
      
#if defined (SYS_CRAY) 
      REAL ZDEN(*), ZETA1(*), WORK(LWORK
      REAL XLAMP1(*), XLAMH1(*), XLAMP2(*), XLAMH2(*)
      REAL ONE, ZERO
#else 
      DOUBLE PRECISION ZDEN(*), ZETA1(*), WORK(LWORK) 
      DOUBLE PRECISION XLAMP1(*), XLAMH1(*), XLAMP2(*), XLAMH2(*)
      DOUBLE PRECISION ONE, ZERO
#endif 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

      INTEGER ISYM12, ISALBE, ISYMAL, ISYMBE, ISYMB, ISYMK
      INTEGER KOFF1, KOFF2, KOFF3, NTOTBE, NTOTAL, NTOTB

      CALL QENTER('CC_ZETADEN')

*---------------------------------------------------------------------*
*     set symmetries:
*---------------------------------------------------------------------*
 
      ISYM12 = MULD2H(ISYXL1,ISYXL2)
      ISALBE = MULD2H(ISYCTR,ISYM12)

*---------------------------------------------------------------------*
*     LambdaP2 x LambdaH1 x Zeta
*---------------------------------------------------------------------*
 
      DO ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMAL,ISALBE)
         ISYMB  = MULD2H(ISYMAL,ISYXL2)
         ISYMK  = MULD2H(ISYMBE,ISYXL1)
C
         IF (LWORK .LT. NBAS(ISYMAL)*NRHF(ISYMK) ) THEN
            CALL QUIT('Insufficient memory in CC_ZETADEN.')
         ENDIF
C
         KOFF1  = IGLMVI(ISYMAL,ISYMB) + 1
         KOFF2  = IT1AM(ISYMB,ISYMK) + 1
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),NVIR(ISYMB),
     *              ONE,XLAMP2(KOFF1),NTOTAL,ZETA1(KOFF2),NTOTB,
     *              ZERO,WORK,NTOTAL)
C
         KOFF2  = IGLMRH(ISYMBE,ISYMK) + 1
         KOFF3  = IAODIS(ISYMAL,ISYMBE) + 1
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHF(ISYMK),
     *              ONE,WORK,NTOTAL,XLAMH1(KOFF2),NTOTBE,
     *              ONE,ZDEN(KOFF3),NTOTAL)

      END DO
 
*---------------------------------------------------------------------*
*     LambdaP1 x LambdaH2 x Zeta
*---------------------------------------------------------------------*
 
      IF (IOPT .EQ. 2) THEN
C
         DO ISYMAL = 1,NSYM
C
            ISYMBE = MULD2H(ISYMAL,ISALBE)
            ISYMB  = MULD2H(ISYMAL,ISYXL1)
            ISYMK  = MULD2H(ISYMBE,ISYXL2)
C
            IF (LWORK .LT. NBAS(ISYMAL)*NRHF(ISYMK) ) THEN
               CALL QUIT('Insufficient memory in CC_ZETADEN.')
            ENDIF
C
            KOFF1  = IGLMVI(ISYMAL,ISYMB) + 1
            KOFF2  = IT1AM(ISYMB,ISYMK) + 1
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTB  = MAX(NVIR(ISYMB),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMK),NVIR(ISYMB),
     *                 ONE,XLAMP1(KOFF1),NTOTAL,ZETA1(KOFF2),NTOTB,
     *                 ZERO,WORK,NTOTAL)
C
            KOFF2  = IGLMRH(ISYMBE,ISYMK) + 1
            KOFF3  = IAODIS(ISYMAL,ISYMBE) + 1
            NTOTAL = MAX(NBAS(ISYMAL),1)
            NTOTBE = MAX(NBAS(ISYMBE),1)
C
            CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHF(ISYMK),
     *                 ONE,WORK,NTOTAL,XLAMH2(KOFF2),NTOTBE,
     *                 ONE,ZDEN(KOFF3),NTOTAL)
C
         END DO
C
      ENDIF
C
      CALL QEXIT('CC_ZETADEN')
C
      RETURN
      END
*=====================================================================*
C  /* Deck cclt_chi */
      SUBROUTINE CCLT_CHI(CTR1,XI,ISYCTR,XLAMDP,ISYLAM,CHI,IOPTX)
C
C     Purpose: To calculate the Chi intermediate:
C
C     Chi(alpha i)  =   sum_c XLAMDP(alpha c) CTR1(c i) 
C                     - sum_k XLAMDP(alpha k)   XI(k i) 
C
C     ISYCTR : symmetry of CTR1, XI
C     ISYLAM : symmetry of XLAMDP
C
C     if IOPTX <> 1 skip contribution from XI intermediate
C
C     Christof Haettig, September 1998
C     based on Asger Halkiers CCLT_CT1AO routine 
C
C
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION CTR1(*),XLAMDP(*),CHI(*),XI(*)
C
      CALL QENTER('CCLT_CHI')
C
C---------------------------------------------
C     Half-transformation to AO-basis of CTR1.
C---------------------------------------------
C
      ISYMAO = MULD2H(ISYCTR,ISYLAM)
      ISYMCJ = ISYCTR 
C
      CALL DZERO(CHI,NGLMDT(ISYMAO))
C
      DO ISYMAL = 1,NSYM
C
         ISYMC = MULD2H(ISYMAL,ISYLAM)
         ISYMJ = MULD2H(ISYMC,ISYMCJ)
C
         KOFF1 = IGLMVI(ISYMAL,ISYMC) + 1
         KOFF2 = IT1AM(ISYMC,ISYMJ)   + 1
         KOFF3 = IGLMRH(ISYMAL,ISYMJ) + 1
C
         NTOTBA = MAX(NBAS(ISYMAL),1)
         NTOTVI = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NVIR(ISYMC),
     *               ONE,XLAMDP(KOFF1),NTOTBA,CTR1(KOFF2),NTOTVI,
     *               ONE,CHI(KOFF3),NTOTBA)
C
         IF (IOPTX.EQ.1) THEN
C
           KOFF4 = IMATIJ(ISYMC,ISYMJ)  + 1
           KOFF5 = IGLMRH(ISYMAL,ISYMC) + 1
C
           NTOTRH = MAX(NRHF(ISYMC),1)
C
           CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NRHF(ISYMC),
     *                -ONE,XLAMDP(KOFF5),NTOTBA,XI(KOFF4),NTOTRH,
     *                 ONE,CHI(KOFF3),NTOTBA)
C
         END IF
C
      END DO
C
      CALL QEXIT('CCLT_CHI')
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck cc_bfif */
      SUBROUTINE CC_BFIF(BFRHF,ISYRHF,XMGD,ISYMGD,
     &                   LUBFI,FNBFI,IADRBF,IADR,IDEL,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: contract effective density with (**|k delta) integrals
*              to the BF intermediate in the F matrix transformation
*              (special version of the CC_BFI routine for F matrix)
*
*     BFRHF  : (**|k delta) integral distribution (symmetry ISYRHF)
*     XMGD   : effective density for BF term (symmetry ISYMGD)
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

      CHARACTER*(*) FNBFI
      INTEGER LWORK, ISYMGD, ISYRHF, IADR, IDEL, LUBFI, IADRBF(*)

#if defined (SYS_CRAY)
      REAL BFRHF(*), XMGD(*), WORK(LWORK)
      REAL ZERO, ONE, HALF
#else
      DOUBLE PRECISION BFRHF(*), XMGD(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, HALF
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
C
      INTEGER NBFRHF(8), IBFRHF(8,8), ISYM, ICOUNT, ISYMAK, ISYBET
      INTEGER ISYMIJ, ISYMIJB, KOFF1, KOFF2, KOFF3, NBASB, NTOTAK
C
      CALL QENTER('CC_BFIF')
C
C     --------------------------------------
C     precalculate symmetry array for BFRHF:
C     --------------------------------------
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBFRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBFRHF(ISYM) = ICOUNT
      END DO
C
      ISYMIJB = MULD2H(ISYRHF,ISYMGD)
C
      IF (LWORK .LT. ND2IJG(ISYMIJB)) THEN
         WRITE (LUPRI,*) 'LWORK =',LWORK
         WRITE (LUPRI,*) 'need ',ND2IJG(ISYMIJB)
         CALL QUIT('Insufficient memory in CC_BFIF.')
      END IF
C
      DO ISYMIJ = 1, NSYM
C
         ISYMAK = MULD2H(ISYMGD,ISYMIJ)
         ISYBET = MULD2H(ISYMAK,ISYRHF)
C
         KOFF1  = IT2AOIJ(ISYMAK,ISYMIJ) + 1
         KOFF3  = IBFRHF(ISYMAK,ISYBET)  + 1
         KOFF2  = ID2IJG(ISYMIJ,ISYBET)  + 1
         NTOTAK = MAX(NT1AO(ISYMAK),1)
         NBASB  = MAX(NBAS(ISYBET),1)
C
         CALL DGEMM('T','N',NBAS(ISYBET),NMATIJ(ISYMIJ),NT1AO(ISYMAK),
     &              ONE, BFRHF(KOFF3),NTOTAK, XMGD(KOFF1),NTOTAK,
     &              ZERO,WORK(KOFF2),NBASB)
C
      END DO
C
      IADRBF(IDEL) = IADR
C
      CALL PUTWA2(LUBFI,FNBFI,WORK,IADR,ND2IJG(ISYMIJB))
C
      IADR = IADR + ND2IJG(ISYMIJB)
C
      CALL QEXIT('CC_BFIF')
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck cc_bfifsort */
      SUBROUTINE CC_BFIFSORT(BFI,ISYRES,LUBFI,FNBFI,IADRBF,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: read and sort the BFZeta intermediate in F mat. transf.
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

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*(*) FNBFI
      INTEGER LWORK, ISYRES, IADR, LUBFI, IADRBF(*)

#if defined (SYS_CRAY)
      REAL BFI(*), WORK(LWORK)
      REAL ONE, TWO
#else
      DOUBLE PRECISION BFI(*), WORK(LWORK)
      DOUBLE PRECISION ONE, TWO
#endif
      PARAMETER(TWO = 2.0D0, ONE = 1.0D0)
C
      INTEGER ISYDEL,ISYMIJB,ISYBET,ISYMI,ISYMJ,ISYMIJ,ISYBEJ,ISYDEI
      INTEGER NDI, NBJ, NIJ, NDIBJ, KOFF1, IDEL, NDB, INDEX
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('CC_BFIFSORT')
C
      CALL DZERO(BFI,NT2AO(ISYRES))
C
      DO ISYDEL = 1, NSYM

         ISYMIJB = MULD2H(ISYDEL,ISYRES)

         IF (LWORK .LT. ND2IJG(ISYMIJB) ) THEN
            CALL QUIT('Insufficient memory in CC_BFIFSORT.')
         END IF

      DO D = 1, NBAS(ISYDEL)
    
         IDEL = D + IBAS(ISYDEL)
         IADR = IADRBF(IDEL)

         CALL GETWA2(LUBFI,FNBFI,WORK,IADR,ND2IJG(ISYMIJB))

         DO ISYMJ = 1, NSYM
         DO ISYMI = 1, NSYM

           ISYMIJ = MULD2H(ISYMI,ISYMJ)
           ISYBET = MULD2H(ISYMIJB,ISYMIJ)
           ISYBEJ = MULD2H(ISYBET,ISYMJ)
           ISYDEI = MULD2H(ISYDEL,ISYMI)
 
           DO J = 1, NRHF(ISYMJ)
           DO I = 1, NRHF(ISYMI)
           DO B = 1, NBAS(ISYBET)

              NIJ   = IMATIJ(ISYMI,ISYMJ)   + NRHF(ISYMI)*(J-1)    + I
              KOFF1 = ID2IJG(ISYMIJ,ISYBET) + NBAS(ISYBET)*(NIJ-1) + B
               
              NBJ   = IT1AO(ISYBET,ISYMJ) + NBAS(ISYBET)*(J-1) + B
              NDI   = IT1AO(ISYDEL,ISYMI) + NBAS(ISYDEL)*(I-1) + D

              IF      (ISYDEI .EQ. ISYBEJ) THEN
                NDIBJ = IT2AO(ISYDEI,ISYBEJ) + INDEX(NDI,NBJ)
                IF (NDI.EQ.NBJ) WORK(KOFF1) = TWO * WORK(KOFF1)
              ELSE IF (ISYDEI .LT. ISYBEJ) THEN
                NDIBJ = IT2AO(ISYDEI,ISYBEJ)+NT1AO(ISYDEI)*(NBJ-1)+NDI
              ELSE IF (ISYDEI .GT. ISYBEJ) THEN
                NDIBJ = IT2AO(ISYDEI,ISYBEJ)+NT1AO(ISYBEJ)*(NDI-1)+NBJ
              END IF

              BFI(NDIBJ) = BFI(NDIBJ) + WORK(KOFF1)

            END DO
            END DO
            END DO

         END DO
         END DO

      END DO
      END DO
C
      IF (LOCDBG) THEN
         CALL AROUND("CC_BFIFSORT: The BF intermediate:")

         WRITE (LUPRI,*) 'START OF BFI:',(BFI(I),I=1,10)

         DO ISYMJ  = 1, NSYM
         DO ISYMI  = 1, NSYM
         DO ISYBET = 1, NSYM
          
           ISYMIJ = MULD2H(ISYMI,ISYMJ)
           ISYBEJ = MULD2H(ISYMJ,ISYBET)
           ISYDEI = MULD2H(ISYBEJ,ISYRES)
           ISYDEL = MULD2H(ISYDEI,ISYMI)

           IF (LWORK .LT. NBAS(ISYBET)*NBAS(ISYDEL) ) THEN
             CALL QUIT('Insufficient memory in CC_BFIFSORT.')
           END IF

           DO J = 1, NRHF(ISYMJ)
           DO I = 1, NRHF(ISYMI)
             
             DO B = 1, NBAS(ISYBET)
             DO D = 1, NBAS(ISYDEL)

               NBJ = IT1AO(ISYBET,ISYMJ) + NBAS(ISYBET)*(J-1) + B
               NDI = IT1AO(ISYDEL,ISYMI) + NBAS(ISYDEL)*(I-1) + D
               NDB = NBAS(ISYDEL)*(B-1) + D

               IF     (ISYDEI .EQ. ISYBEJ) THEN
                 NDIBJ = IT2AO(ISYDEI,ISYBEJ) + INDEX(NDI,NBJ)
               ELSEIF (ISYDEI .LT. ISYBEJ) THEN
                 NDIBJ = IT2AO(ISYDEI,ISYBEJ)+NT1AO(ISYDEI)*(NBJ-1)+NDI
               ELSEIF (ISYDEI .GT. ISYBEJ) THEN
                 NDIBJ = IT2AO(ISYDEI,ISYBEJ)+NT1AO(ISYBEJ)*(NDI-1)+NBJ
               ENDIF
                  
               WORK(NDB) = BFI(NDIBJ)

             END DO
             END DO

             WRITE (LUPRI,*) 'ISYMJ,ISYMI,ISYBET,ISYDEL,J,I:',
     &                 ISYMJ,ISYMI,ISYBET,ISYDEL,J,I

             CALL OUTPUT(WORK,1,NBAS(ISYDEL),1,NBAS(ISYBET),
     &                   NBAS(ISYDEL),NBAS(ISYBET),1,LUPRI)

           END DO
           END DO

         END DO
         END DO
         END DO

      END IF
C
      CALL QEXIT('CC_BFIFSORT')
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_GIM(DSRHF,ISYRHF,XMGD,ISYMGD,GIM,IOPT,WORK,LWORK)
*---------------------------------------------------------------------*
*
* Purpose: Calculate G intermediate from effective BF density and
*          (**|k del) integrals
*
*    IOPT = 0  :  compute left-hand-side GZeta intermediate
*
*    IOPT = 1  :  compute right-hand-side G intermediate, using
*                 2 Cou - Exc combination of XMGD 
*                 Note: this overwrites XMGD!
*
*    symmetries:     ISYRHF  --  DSRHF
*                    ISYMGD  --  XMGD
*
* Christof Haettig, November 1998
*---------------------------------------------------------------------*
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      DIMENSION DSRHF(*),XMGD(*),GIM(*),WORK(LWORK)
C
      CALL QENTER('CC_GIM')
C
C     ----------------------------------------------------------
C     if right-hand-side G requested replace XMGD by 2Cou - Exc:
C     ----------------------------------------------------------
      IF (IOPT.EQ.1) THEN
         CALL CC_MGDTCME(XMGD,ISYMGD,WORK,LWORK)
      END IF
C
      DO ISYMK = 1,NSYM
C
         ISYMAG = MULD2H(ISYMK,ISYRHF)
         ISYMGI = MULD2H(ISYMK,ISYMGD)
C
         IF (LWORK .LT. N2BST(ISYMAG)) THEN
            CALL QUIT('Insufficient space in CC_GIM')
         ENDIF
C
         DO K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYMAG,ISYMK) + NNBST(ISYMAG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK)
C
            DO ISYMI = 1,NSYM
C
               ISYMG = MULD2H(ISYMI,ISYMGI)
               ISYMA = MULD2H(ISYMG,ISYMAG)
C
               KOFF4 = IAODIS(ISYMA,ISYMG) + 1
               KOFF5 = IT2BGD(ISYMGI,ISYMK) + NT1AO(ISYMGI)*(K - 1) 
     *               + IT1AO(ISYMG,ISYMI) + 1
               KOFF6 = IT1AO(ISYMA,ISYMI) + 1
C
               NBASG = MAX(NBAS(ISYMG),1)
               NBASA = MAX(NBAS(ISYMA),1)
C
               CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF4),NBASA,XMGD(KOFF5),NBASG,
     *                    ONE,GIM(KOFF6), NBASA)
C
            END DO
         END DO
      END DO
 
      CALL QEXIT('CC_GIM')
C
      RETURN
      END
*=====================================================================*
*                   END OF SUBROUTINE CC_GIM                          *
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_MGDTCME(XMGD,ISYMGD,WORK,LWORK)
*---------------------------------------------------------------------*
*
* Purpose: Calculate 2 Cou - Exc combination of two-index AO back-
*          transformed amplitude batch XMGD(gam i;j)^del
*
* Christof Haettig, November 1998
*---------------------------------------------------------------------*
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      PARAMETER(TWO=2.0D0,ONE=1.0D0)
      DIMENSION XMGD(*),WORK(LWORK)
C
      CALL QENTER('CC_MGDTCME')
C
      DO ISYMJ = 1,NSYM
C
        ISYMGI = MULD2H(ISYMJ,ISYMGD)
C
        DO J = 1,NRHF(ISYMJ)
C
          DO ISYMI = 1,ISYMJ
C
            ISYMG  = MULD2H(ISYMI,ISYMGI)
            ISYMGJ = MULD2H(ISYMG,ISYMJ)
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NBAS(ISYMG)
            KEND1 = KSCR2 + NBAS(ISYMG)
            LWRK1 = LWORK - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space in CC_MGDTCME')
            ENDIF
C
            IF (ISYMI .EQ. ISYMJ) THEN
               NRHFI = J - 1
            ELSE
               NRHFI = NRHF(ISYMI)
            END IF
C
            DO I = 1,NRHFI
C
              NGIJ = IT2BGD(ISYMGI,ISYMJ) + NT1AO(ISYMGI)*(J-1)
     *             + IT1AO(ISYMG,ISYMI) + NBAS(ISYMG)*(I-1) + 1
              NGJI = IT2BGD(ISYMGJ,ISYMI) + NT1AO(ISYMGJ)*(I-1)
     *             + IT1AO(ISYMG,ISYMJ) + NBAS(ISYMG)*(J-1) + 1
C
              CALL DCOPY(NBAS(ISYMG),XMGD(NGIJ),1,WORK(KSCR1),1)
              CALL DCOPY(NBAS(ISYMG),XMGD(NGJI),1,WORK(KSCR2),1)
              CALL DSCAL(NBAS(ISYMG),TWO,XMGD(NGIJ),1)
              CALL DSCAL(NBAS(ISYMG),TWO,XMGD(NGJI),1)
              CALL DAXPY(NBAS(ISYMG),-ONE,WORK(KSCR2),1,XMGD(NGIJ),1)
              CALL DAXPY(NBAS(ISYMG),-ONE,WORK(KSCR1),1,XMGD(NGJI),1)
C
            END DO
          END DO
        END DO
      END DO
C
      CALL QEXIT('CC_MGDTCME')
C
      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CC_MGDTCME                       *
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCSDT_F2B_SETUP(IFTRAN,IFDOTS,FCON,NFTRAN,MXVEC,LADD,
     &                           IB1TRAN,IB1DOTS,B1CON,NB1TRAN,MXB1VEC,
     &                           IB2TRAN,IB2DOTS,B2CON,NB2TRAN,MXB2VEC,
     &                           LISTL,LISTB,LISTC,MXB1TRAN,MXB2TRAN)
*---------------------------------------------------------------------*
*
* Purpose: transform from F matrix like loop structure to a B matrix
*          like loop structure needed for CCSDT_FBMAT_CONT
*
* Christof Haettig, May 2002
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"

      CHARACTER*(25) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCSDT_F2B_SETUP> ')  
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
*--------
* input:
*--------
      CHARACTER LISTL*(3), LISTC*(3), LISTB*(3)
      LOGICAL LADD
      INTEGER NFTRAN, MXVEC, MXB1VEC, MXB1TRAN, MXB2VEC, MXB2TRAN
      INTEGER IFTRAN(3,NFTRAN), IFDOTS(MXVEC,NFTRAN)
#if defined (SYS_CRAY) 
      REAL FCON(MXVEC,NFTRAN)
#else
      DOUBLE PRECISION FCON(MXVEC,NFTRAN)
#endif
      
*--------
* output:
*--------
      INTEGER NB1TRAN, NB2TRAN
      INTEGER IB1TRAN(3,MXB1TRAN), IB1DOTS(MXB1VEC,MXB1TRAN)
      INTEGER IB2TRAN(3,MXB2TRAN), IB2DOTS(MXB2VEC,MXB2TRAN)
#if defined (SYS_CRAY) 
      REAL B1CON(MXB1VEC,MXB1TRAN), B2CON(MXB2VEC,MXB2TRAN)
#else
      DOUBLE PRECISION B1CON(MXB1VEC,MXB1TRAN), B2CON(MXB2VEC,MXB2TRAN)
#endif

* local:
      INTEGER MXVAB1, MXVAB2, ITRAN, IZETA, ITAMPB, ITAMPC, IVEC,
     *        JVEC, JTRAN
#if defined (SYS_CRAY) 
      REAL CONT1, CONT2
#else
      DOUBLE PRECISION CONT1, CONT2
#endif

      CALL QENTER('CCFCSU')

      IF (.NOT. LADD) THEN
        NB1TRAN = 0
        NB2TRAN = 0
      END IF
      MXVAB1 = 0
      MXVAB2 = 0

      DO ITRAN = 1, NFTRAN

        IZETA  = IFTRAN(1,ITRAN)
        IF      (IZETA.EQ.0 .AND. LISTL(1:3).EQ.'L0 ') THEN
          IZETA = 1
        ELSE IF (IZETA.EQ.0) THEN
          CALL QUIT('Illegal index for left vector...')
        END IF

        ITAMPB = IFTRAN(2,ITRAN)

        IVEC = 1
        DO WHILE(IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
          ITAMPC = IFDOTS(IVEC,ITRAN)

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) MSGDBG,'IZETA,ITAMPB,ITAMPC:',
     &                             IZETA,ITAMPB,ITAMPC
          END IF

          ! R3=ITAMPB, R1=ITAMPC
          CALL CC_SETF12(IB1TRAN,IB1DOTS,MXB1TRAN,MXB1VEC,
     &                   ITAMPB,ITAMPC,IZETA,JTRAN,JVEC)
          NB1TRAN = MAX(NB1TRAN,JTRAN)
          MXVAB1  = MAX(MXVAB1,JVEC)
          CONT1   = B1CON(JVEC,JTRAN)

          IF (LISTB(1:3).EQ.LISTC(1:3)) THEN
           ! R3=ITAMPC, R1=ITAMPB on the same list
           CALL CC_SETF12(IB1TRAN,IB1DOTS,MXB1TRAN,MXB1VEC,
     &                    ITAMPC,ITAMPB,IZETA,JTRAN,JVEC)
           NB1TRAN = MAX(NB1TRAN,JTRAN)
           MXVAB1  = MAX(MXVAB1,JVEC)
           CONT2   = B1CON(JVEC,JTRAN)
          ELSE
           ! R3=ITAMPC, R1=ITAMPB on a second list
           CALL CC_SETF12(IB2TRAN,IB2DOTS,MXB2TRAN,MXB2VEC,
     &                    ITAMPC,ITAMPB,IZETA,JTRAN,JVEC)
           NB2TRAN = MAX(NB2TRAN,JTRAN)
           MXVAB2  = MAX(MXVAB2,JVEC)
           CONT2   = B2CON(JVEC,JTRAN)
          END IF

          IF (LADD) THEN
            FCON(IVEC,ITRAN) = FCON(IVEC,ITRAN) + CONT1 + CONT2
         
            IF (LOCDBG) THEN
              WRITE(LUPRI,*)'ITRAN,IVEC:',ITRAN,IVEC
              WRITE(LUPRI,*)'FCON(before):',FCON(IVEC,ITRAN)-CONT1-CONT2
              WRITE(LUPRI,'(A,3I5)') 'IZETA,ITAMPB,ITAMPC:',
     &                                IZETA,ITAMPB,ITAMPC
              WRITE(LUPRI,*) '<L|[[H,RC_1],RB_3]|HF> = ',CONT1
              WRITE(LUPRI,*) '<L|[[H,RB_1],RC_3]|HF> = ',CONT2
              WRITE(LUPRI,*) 'FCON(after):',FCON(IVEC,ITRAN)
            END IF
          END IF

          IVEC = IVEC + 1
        END DO

      END DO
C
      IF (LOCDBG) THEN
        WRITE (LUPRI,*)
        WRITE (LUPRI,*) 'List of <L2|[[H,RC_1],RB_3]|HF>:'
        DO ITRAN = 1, NB1TRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IB1TRAN(I,ITRAN),I=1,2),(IB1DOTS(I,ITRAN),I=1,MXVAB1)
        END DO
        WRITE (LUPRI,*)
        IF (LISTB(1:3).NE.LISTC(1:3)) THEN
          WRITE (LUPRI,*) 'List of <L2|[[H,RB_1],RC_3]|HF>:'
          DO ITRAN = 1, NB2TRAN
            WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &       (IB2TRAN(I,ITRAN),I=1,2),(IB2DOTS(I,ITRAN),I=1,MXVAB2)
          END DO
        WRITE (LUPRI,*)
        END IF
      END IF
C
      CALL QEXIT('CCFCSU')
C
      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CCSDT_F2B_SETUP                  *
*=====================================================================*

*---------------------------------------------------------------------*
c/* Deck CC_R12FMAT */
*=====================================================================*
       SUBROUTINE CC_R12FMAT(THETA1, THETA2, THETAR12, 
     &                       ISYRES, LISTL, IDLSTL,
     &                       ZETA1, ISYCTR, LISTR, IDLSTR,  
     &                       ISYAMP, LAMDPA, LAMDHA, LAMP0, LAMH0,
     &                       XINT,LUMIM,FILMIM,LENM,IINT2,WORK,LWRK)  
*---------------------------------------------------------------------*
*
*    Purpose: calculate R12 contributions for F-matrix
*
*    C. Neiss, C. Haettig  autumn 2004
*    Added finite fields, Chr. Neiss, June 2005
*    Added CCSD(R12), Chr. Neiss, Nov. 2005 
*
*    THETA1        singles part of result vector
*    THETA2        doubles part of result vector
*    THETAR12      R12-doubles part of result vector
*    ISYRES        symmetry number of result vector
*    LISTL         type of L-vector; see CC_RDRSP.F
*    IDLSTL        index of L-vector; see CC_RDRSP.F
*    ZETA1         singles part of Lagrangian multipliers
*    ISYCTR        symmetry of Lagrangian multipliers
*    LISTR         type of right response vector
*    IDLSTR        index of right response vector
*    ISYAMP        symmetry of right response vector
*    LAMDPA        Lambda^p-Matrix contracted with singles part of 
*                  response vector
*    LAMDHA        Lambda^h-Matrix contracted with singles part of
*                  response vector
*    LAMP0         Lambda^p-Matrix
*    LAMH0         Lambda^h-Matrix
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
    
      LOGICAL LV, LVIJKL, LVAJKL

      INTEGER LWRK, ISYRES, IDLSTL, ISYCTR, IDLSTR, ISYAMP 
      INTEGER KEIM, ISYME, LWRK1, LWRK2, KEND1, KEND2, KEND3, LWRK3,
     &        KR12SCR, NR12SCR, KVAJKL,
     &        KZETA12, KCTR2, IOPT
      INTEGER ISYM1,ISYM2,LUNIT
      INTEGER LUMIM,LENM,IINT2,KMINT
      INTEGER KVABKL, KVINT
      INTEGER KT12AMP, KXINTTRI, KXINTSQ, KSCR, KFIELDAO, IFLD, IAN
      CHARACTER LISTL*3, LISTR*3, CDUMMY*8, MODEL*10
      CHARACTER*(*) FILMIM

#if defined (SYS_CRAY)
      REAL WORK(LWRK), THETA1(*), THETA2(*), THETAR12(*), ZETA1(*), 
     &     LAMDPA(*), LAMDHA(*), LAMP0(*), LAMH0(*), XINT(*)
      REAL ZERO, ONE, TWO, HALF, DDOT
      REAL TIM1, TIM2, TIM3
#else
      DOUBLE PRECISION WORK(LWRK), THETA1(*), THETA2(*), THETAR12(*),
     &                 ZETA1(*), XINT(*), 
     &                 LAMDPA(*), LAMDHA(*), LAMP0(*), LAMH0(*)
      DOUBLE PRECISION ZERO, ONE, TWO, HALF, DDOT
      DOUBLE PRECISION TIM1, TIM2, TIM3
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)

  
      CALL QENTER('CC_R12FMAT')

      IF (IANR12 .NE. 1) THEN
        CALL QUIT('Only available for CC-R12/ANSATZ 1')
      END IF
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) "Theta1 and Theta2 before R12 section:"
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
      END IF

      IF (CC2) THEN
C     --------------------------------------------------------------
C     calculate G'-Term Singles contribution:
C     do first E-intermediate calculation, then contract with 
C     singles Lagrangian multipliers and add to conventional term
C     --------------------------------------------------------------
C
C     E-Intermediate:
C
      ISYME = ISYAMP
      IF (ISYRES .NE. MULD2H(ISYME,ISYCTR)) THEN
        WRITE(LUPRI,*) 'isyres, isyamp, isyctr: ',isyres, isyamp, isyctr
        CALL QUIT('Symmetry mismatch in cc_r12fmat!')
      END IF       

      KEIM  = 1 
      KEND1 = KEIM + NMATIJ(ISYME)
      LWRK1 = LWRK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in cc_r12fmat')
      END IF 
  
      CALL DZERO(WORK(KEIM),NMATIJ(ISYME))
      CALL CCRHS_EINTP(WORK(KEIM),LAMP0,WORK(KEND1),LWRK1,
     &                    2,ISYAMP,CDUMMY,IDUMMY,IDUMMY,LISTR,IDLSTR)
C
C     contraction with singles Lagrange-multipliers:
C
      KEND2 = KEND1 + NMATAB(ISYME)
      LWRK2 = LWRK - KEND2
      IF (LWRK2 .LT. 0) THEN
        CALL QUIT('Insufficient work space in cc_r12fmat')
      END IF

      CALL DZERO(WORK(KEND1),NMATAB(ISYME))

      CALL CCLR_E1C1(THETA1,ZETA1,WORK(KEND1),WORK(KEIM),WORK(KEND2),
     &               LWRK2,ISYCTR,ISYME,'T')          
C      
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) "ZETA1:"
         CALL CC_PRP(ZETA1,DUMMY,ISYCTR,1,0)
         WRITE(LUPRI,*) "E Intermediates:"
         CALL CC_PREI(WORK(KEND1),WORK(KEIM),ISYME,1)
         write (lupri,*) "Norm^2 (E Intermediates)",
     &      ddot(nmatab(isyme),work(kend1),1,work(kend1),1),
     &      ddot(nmatij(isyme),work(keim),1,work(keim),1)
         write(lupri,*) "Theta1 after G' term:"
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,0)
      END IF
C
      END IF !CC2
C
C     --------------------------------------------------------------
C     calculate F'-Term Singles contribution
C     --------------------------------------------------------------
C
C     reallocate memory
C
      KVAJKL = 1
      KEND1 = KVAJKL + NVAJKL(ISYAMP)
      LWRK1 = LWRK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in cc_r12fmat')
      END IF
C
C     calculate V_{kl}^{alpha jbar}-Intermediate
C
      IF (.NOT.USEVABKL) THEN
        IOPT = 1 
        CALL CC_R12MKVAMKL0(WORK(KVAJKL),NVAJKL(ISYAMP),IOPT,LAMDHA,
     &                      ISYAMP,WORK(KEND1),LWRK1)
        CALL CC_MOFCONR12(LAMDHA,ISYAMP,IDUMMY,IDUMMY,IDUMMY,
     &                    IDUMMY,DUMMY,DUMMY,WORK(KVAJKL),IDUMMY,
     &                    .FALSE.,.TRUE.,.FALSE.,2,
     &                    TIM1,TIM2,TIM3,
     &                    IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     &                    WORK(KEND1),LWRK1)
C
      ELSE
C
        KVABKL = KEND1
        KEND2  = KVABKL + NVABKL(1)
        LWRK2  = LWRK - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insuff. work space for V^(albe)_kl in cc_r12fmat')
        END IF
C
        LV = .TRUE.
        LVAJKL = .FALSE.
        LVIJKL = .FALSE.
        CALL CC_R12MKVTF(WORK(KVABKL),WORK(KVAJKL),DUMMY,
     &                   LAMDHA,ISYAMP,
     &                   LV,LVIJKL,LVAJKL,CDUMMY,WORK(KEND2),LWRK2)
C
      END IF
C
      if (locdbg) then
        write(lupri,*) 'Norm^2 of VAJbarKL in CC_R12FMAT: ',
     &             ddot(nvajkl(isyamp),work(kvajkl),1,work(kvajkl),1)
        write(lupri,*) 'VAJbarKL in CC_R12FMAT:'
        do isym2 = 1, nsym
          isym1 = muld2h(isyamp,isym2)
          write(lupri,*) 'Block isymaj, isymkl: ',isym1,isym2
          if ((nt1ao(isym1).eq.0).or.(nmatkl(isym2).eq.0)) then
           write(lupri,*) 'This block is empty'
          else 
           call output(work(kvajkl+ivajkl(isym1,isym2)),1,nt1ao(isym1),
     &                 1,nmatkl(isym2),nt1ao(isym1),nmatkl(isym2),1,
     &                 lupri)
          end if
        end do
      end if
C
C     read in r12 doubles Lagrangian multipliers...
C 
      KZETA12 = KEND1
      KEND1   = KZETA12 + NTR12SQ(ISYCTR)
      LWRK1   = LWRK - KEND1
      
      CALL CC_R12GETCT(WORK(KZETA12),ISYCTR,2,2.0D0*BRASCL,.FALSE.,
     &                 'T',DUMMY,CDUMMY,DUMMY,LISTL,IDLSTL,
     &                 WORK(KEND1),LWRK1)
C
C     ...and contract:
C
      CALL CCRHS_GP0(THETA1,ISYRES,WORK(KVAJKL),ISYAMP,LAMH0,
     &               1,WORK(KZETA12),ISYCTR,.FALSE.,DUMMY,LOCDBG,
     &               WORK(KEND1),LWRK1)
C
C     ---------------------------------------------------------------
C     allocate memory for R12-doubles part of result vector (squared)
C     ---------------------------------------------------------------
C
      KR12SCR = 1
      KEND1   = KR12SCR + NTR12SQ(ISYRES)
      LWRK1   = LWRK - KEND1
      IF (LWRK1 .LT. 0) THEN
        CALL QUIT('Insufficient work space in cc_r12fmat')
      END IF
      CALL DZERO(WORK(KR12SCR),NTR12SQ(ISYRES))
C
C     --------------------------------------------------------------
C     calculate G'-Term R12-doubles contribution:
C     --------------------------------------------------------------
C
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Call now CCLHTR_GP!'
      END IF
      CALL CCLHTR_GP(ZETA1,ISYCTR,LAMDPA,ISYAMP,WORK(KR12SCR),ISYRES,
     &               WORK(KEND1),LWRK1) 
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'KR12SCR after CCLHTR_GP'
        call cc_prsqr12(WORK(KR12SCR),isyres,'T',1,.false.)
      END IF
C
C     ------------------
C     add terms for CCSD
C     ------------------
C
      IF (CCSD .OR. CCSDT) THEN
C
C       add B' term
C
C       !read V_(alpha beta)^(kl) from disk
        KVABKL = KEND1
        KEND2  = KVABKL + NVABKL(1)
        LWRK2  = LWRK - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insuff. work space for V^(albe)_kl in CC_R12FMAT')
        END IF
        lunit = -1
        call gpopen(lunit,fvabkl,'unknown',' ','unformatted',
     &            idummy,.false.)
        read (lunit)(work(kvabkl+i-1),i=1,nvabkl(1))
        call gpclose(lunit,'KEEP')
C
        IOPT = 2
        CALL CC_R12MKVIRT(WORK(KVABKL),LAMDPA,ISYAMP,LAMP0,1,
     &                    'R12VCBDTKL',IOPT,WORK(KEND2),LWRK2)
C
        !read doubles Lagrangian multipliers
        KCTR2 = KEND1
        KEND2 = KCTR2 + NT2AM(ISYCTR)
        LWRK2  = LWRK - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insuff. work space in CC_R12FMAT')
        END IF
        IOPT = 2
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,DUMMY,WORK(KCTR2))
C
        !calculate contribution
        CALL CCRHS_BPP(WORK(KR12SCR),WORK(KCTR2),ISYCTR,.TRUE.,
     &                 'R12VCBDTKL',ISYAMP,WORK(KEND2),LWRK2)
C
C       add A' term    
C
        !read V-intermediate from disk:
        KVINT = KEND1
        KEND2 = KVINT + NTR12AM(1)
        LWRK2 = LWRK - KEND2
        IF (LWRK2 .LT. 0) THEN
          WRITE(LUPRI,*) 'Available:', LWRK, 'Needed:', KEND2
          CALL QUIT('Insufficient memory for V-intermediate '//
     &              'in CC_R12FMAT')
        ENDIF
        lunit = -1
        call gpopen(lunit,fccr12v,'old',' ','unformatted',
     &                  idummy,.FALSE.)
 9999   read(lunit) ian
        read(lunit) (work(kvint+i), i=0, ntr12am(1)-1)
        if (ian.ne.ianr12) goto 9999
        call gpclose(lunit,'KEEP')
C
        !read M-intermediate from disk:
        KMINT = KEND2
        KEND3 = KMINT + N3ORHF(ISYRES)
        LWRK3 = LWRK - KEND3
        IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insuff. work space in CC_R12FMAT')
        END IF
        CALL CC_RVEC(LUMIM,FILMIM,LENM,N3ORHF(ISYRES),IINT2,
     &               WORK(KMINT))
C
        CALL CCLHTR_AP(WORK(KR12SCR),WORK(KVINT),WORK(KMINT),
     &                 ISYRES,WORK(KEND3),LWRK3)
C
C       add E' term
C            
        CALL CCLHTR_EP(WORK(KR12SCR),WORK(KVINT),XINT,
     &                 ISYRES,WORK(KEND2),LWRK2)
      END IF
C
C     resort result into a symmetry packed triangular matrix and 
C     apply the projection operator 0.5*P_(ij)^(kl):
C
      IOPT = 1
      CALL CCR12PCK2(THETAR12,ISYRES,.TRUE.,WORK(KR12SCR),'T',
     &               IOPT)
      CALL CCLR_DIASCLR12(THETAR12,0.5D0*KETSCL,ISYRES) 
      CALL DSCAL(NTR12AM(ISYRES),0.5D0,THETAR12,1)

      if (locdbg) then 
        write(lupri,*) 'theta_r12 after packing'
        call cc_prpr12(thetar12,isyres,1,.false.)
      end if

C     ---------------------------------------------------
C     add finite field contributions:
C     ---------------------------------------------------
      IF (NONHF) THEN
C       allocate memory:
        KZETA12 = KEND1
        KT12AMP = KZETA12 + NTR12SQ(ISYCTR)
        KXINTTRI= KT12AMP + NTR12SQ(ISYAMP)
        KXINTSQ = KXINTTRI+ NR12R12P(1)
        KSCR    = KXINTSQ + NR12R12SQ(1)
        KFIELDAO= KSCR    + NTR12SQ(ISYRES)
        KEND2   = KFIELDAO+ N2BST(1)
        LWRK2   = LWRK - KEND2
        IF (LWRK2.LT.0) THEN
          CALL QUIT('Insufficient work space for R12 in CC_R12FMAT')
        END IF

C       initialize fields:
        CALL DZERO(WORK(KFIELDAO),N2BST(1))
        CALL DZERO(WORK(KSCR),NTR12SQ(ISYRES))
        CALL DZERO(WORK(KR12SCR),NTR12SQ(ISYRES))
C
        IF (ISYMOP.NE.1) CALL QUIT('ISYMOP .NE. 1 not implemented')
       
C       sum up fields:
        DO IFLD = 1, NFIELD
          IF ( NHFFIELD(IFLD) ) THEN
            CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND2),LWRK2,
     *                   EFIELD(IFLD),1,LFIELD(IFLD))
          ELSE IF (.NOT. NHFFIELD(IFLD)) THEN
            CALL QUIT('CCR12 response can only handle '//
     &                'unrelaxed orbitals (w.r.t. the perturbation)')
          END IF
        END DO

C       read r12 Lagrangian multipliers:
        CALL CC_R12GETCT(WORK(KZETA12),ISYCTR,2,2.0D0*BRASCL,.FALSE.,
     &                   'N',DUMMY,CDUMMY,DUMMY,LISTL,IDLSTL,
     &                   WORK(KEND2),LWRK2)

C       read r12 response amplitudes:
        CALL CC_R12GETCT(WORK(KT12AMP),ISYAMP,2,KETSCL,.FALSE.,'N',
     &               DUMMY,CDUMMY,DUMMY,LISTR,IDLSTR,WORK(KEND2),LWRK2)

C       read r12 overlap matrix:
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
C       calculate singles contribution:
C
        CALL CC_R12ETAA(THETA1,ISYRES,WORK(KZETA12),ISYCTR,
     &                  WORK(KT12AMP),ISYAMP,WORK(KXINTSQ),
     &                  WORK(KFIELDAO),1,LAMP0,LAMH0,.FALSE.,
     &                  WORK(KEND2),LWRK2)
C
C       calculate r12 doubles contribution:
C 
        CALL CC_R12XI2A(WORK(KSCR),ISYRES,WORK(KZETA12),ISYCTR,
     &                  WORK(KFIELDAO),1,LAMP0,LAMDHA,
     &                  ISYAMP,'T',WORK(KEND2),LWRK2)
        CALL CC_R12XI2B(WORK(KR12SCR),'N',WORK(KXINTSQ),1,'N',
     &                  WORK(KSCR),ISYRES,'N',-ONE)
C
        IOPT = 1
        CALL CCR12PCK2(WORK(KSCR),ISYRES,.FALSE.,WORK(KR12SCR),
     &                 'N',IOPT)
        CALL CCLR_DIASCLR12(WORK(KSCR),0.5D0*KETSCL,ISYRES)
        CALL DAXPY(NTR12AM(ISYRES),ONE,WORK(KSCR),1,THETAR12,1)
        
      END IF

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) "Theta at end of CC_R12FMAT:"
         CALL CC_PRP(THETA1,THETA2,ISYRES,1,1)
         CALL CC_PRPR12(THETAR12,ISYRES,1,.TRUE.)
      END IF

      CALL QEXIT('CC_R12FMAT')
      
      RETURN
      END

*=====================================================================*
*                  END OF SUBROUTINE CC_R12FMAT                  *
*=====================================================================*
