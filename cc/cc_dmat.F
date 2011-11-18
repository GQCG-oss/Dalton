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
c/* Deck CC_DMAT */
*=====================================================================*
      SUBROUTINE CC_DMAT( IDTRAN,  NDTRAN, LISTA, LISTB, LISTC, LISTD,
     &                    IOPTRES, FILDMA, IDDOTS,DCONS, MXVEC,
     &                    WORK,    LWORK                              )
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
*                LISTD       -- type of T^D vectors
*                IDTRAN(1,*) -- indeces of T^A vectors
*                IDTRAN(2,*) -- indeces of T^B vectors
*                IDTRAN(3,*) -- indeces of T^C vectors
*                IDTRAN(4,*) -- indeces of T^D vectors
*                IDTRAN(5,*) -- indeces or addresses of result vectors
*                NDTRAN      -- number of requested transformations
*                FILDMA      -- file name / list type of result vectors
*                               or list type of vectors to be dotted on
*                IDDOTS      -- indeces of vectors to be dotted on
*                DCONS       -- contains the dot products on return
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  all result vectors are written to a direct
*                          access file, FILDMA is used as file name
*                          the start addresses of the vectors are
*                          returned in IDTRAN(5,*)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          IDTRAN(3,5). N.B.: if WORK is not large
*                          enough iopt is automatically reset to 0!!!
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, FILDMA is used
*                          as list type and IDTRAN(5,*) as list index
*                          NOTE that IDTRAN(5,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WRRSP, FILCMA is used
*                          as list type and ICTRAN(4,*) as list index
*                          NOTE that ICTRAN(4,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by FILDMA and the indeces from IDDOTS
*                          the result of the dot products is returned
*                          in the DCONS array
*
*
*     Written by Christof Haettig, may 1997.
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

* local parameters:
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CC_DMAT> ')

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
      
      INTEGER LUDMAT

      CHARACTER*(*) LISTA, LISTB, LISTC, LISTD, FILDMA
      INTEGER IOPTRES
      INTEGER NDTRAN, MXVEC, LWORK
      INTEGER IDTRAN(5,NDTRAN)
      INTEGER IDDOTS(MXVEC,NDTRAN)

#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL ZERO, ONE, TWO, HALF
      REAL DCONS(MXVEC,NDTRAN)
#else 
      DOUBLE PRECISION DCONS(MXVEC,NDTRAN)
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
#endif 
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, HALF = 0.5d0)


      CHARACTER*(10) MODEL, MODELW
      INTEGER INDEXA(MXCORB)
      INTEGER IOPTW, IOPT, ITRAN, IADRTH, IOPTRW
      INTEGER LENALL, LEN, IOFFCD, ICYCLE, ILLL, IERROR
      INTEGER KTHETA0, KTHETA1, KTHETA2
      INTEGER KT1AMPA, KT1AMPB, KT1AMPC, KT1AMPD
      INTEGER IDLSTA, IDLSTB, IDLSTC, IDLSTD
      INTEGER ISYMTA, ISYMTB, ISYMTC, ISYMTD, ISYRES
      INTEGER NTOSYM, KCCFB1, KINDXB, KFREE, LFREE, KENDSV, LWRKSV
      INTEGER KODCL1, KODBC1, KRDBC1, KODPP1, KRDPP1, KRECNR
      INTEGER KODCL2, KODBC2, KRDBC2, KODPP2, KRDPP2, NUMDIS
      INTEGER IDEL2, ISYDEL, IDEL, KXINT, KDSRHF
      INTEGER KT1AMP0, KLAMDH0, KLAMDP0, KLAMDHA, KLAMDPA
      INTEGER KLAMDHB, KLAMDPB, KLAMDPC, KLAMDHC, KLAMDHD, KLAMDPD
      INTEGER NTAMP, NTOT, ISYMD1
      INTEGER KENDF1, LENDF1, KENDF2, LENDF2, KEND, LEND, KEND0, LEND0
      
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
  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        Call AROUND('ENTERED CC_DMAT')
        IF (DIRECT) WRITE(LUPRI,'(/1X,A)') 'AO direct transformation'
        WRITE (LUPRI,*) MSGDBG,'LISTA : ',LISTA
        WRITE (LUPRI,*) MSGDBG,'LISTB : ',LISTB
        WRITE (LUPRI,*) MSGDBG,'LISTC : ',LISTC
        WRITE (LUPRI,*) MSGDBG,'LISTD : ',LISTD
        WRITE (LUPRI,*) MSGDBG,'FILDMA: ',FILDMA
        WRITE (LUPRI,*) MSGDBG,'NDTRAN: ',NDTRAN
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
        WRITE(LUPRI,'(/1x,a)') 'CC_DMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_DMAT...'
        CALL QUIT('Unknown CC method in CC_DMAT.')
      END IF

      IF (    LISTA(1:1).NE.'R' 
     &   .OR. LISTB(1:1).NE.'R'
     &   .OR. LISTC(1:1).NE.'R'
     &   .OR. LISTD(1:1).NE.'R' ) THEN
        WRITE(LUPRI,*)'LISTA, LISTB, LISTC or LISTD must refer to',
     &                    ' t-amplituded vectors in CC_DMAT.'
        CALL QUIT('Illegal LISTA or LISTB or LISTC in CC_DMAT.')
      END IF

      IF (ISYMOP .NE. 1) THEN
        WRITE(LUPRI,*) 'ISYMOP = ',ISYMOP
        WRITE(LUPRI,*) 'CC_DMAT is not implemented for ISYMOP.NE.1'
        CALL QUIT('CC_DMAT is not implemented for ISYMOP.NE.1')
      END IF

C     IF (NDTRAN .GT. MAXSIM) THEN
C       WRITE(LUPRI,*) 'NDTRAN = ', NDTRAN
C       WRITE(LUPRI,*) 'MAXSIM = ', MAXSIM
C       WRITE(LUPRI,*) 'number of requested transformation is larger'
C       WRITE(LUPRI,*) 'than the maximum number of allowed ',
C    &                 'simultaneous transformation.'
C       CALL QUIT('Error in CC_DMAT: NDTRAN is larger than MAXSIM.')
C     END IF

* check return option for the result vectors:
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
        LUDMAT = -1
        CALL WOPEN2(LUDMAT,FILDMA,64,0)

      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
        CONTINUE
      ELSE IF (IOPTRES.EQ.5) THEN
        IF (MXVEC*NDTRAN.NE.0) CALL DZERO(DCONS,MXVEC*NDTRAN)
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_DMAT.')
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
         CALL QUIT('Unknown coupled cluster model in CC_DMAT.')
      END IF

* start of scratch space for the following:
      KEND0 = 1
      LEND0 = LWORK

*=====================================================================*
* initialize result vectors with zeros:
*=====================================================================*
      IADRTH = 1
      DO ITRAN = 1, NDTRAN

        ISYMTA = ILSTSYM(LISTA,IDTRAN(1,ITRAN))
        ISYMTB = ILSTSYM(LISTB,IDTRAN(2,ITRAN))
        ISYMTC = ILSTSYM(LISTC,IDTRAN(3,ITRAN))
        ISYMTD = ILSTSYM(LISTD,IDTRAN(4,ITRAN))

        ISYRES = MULD2H(MULD2H(ISYMTA,ISYMTB),MULD2H(ISYMTC,ISYMTD))

        IF (CCS) THEN
          NTAMP = NT1AM(ISYRES)
        ELSE
          NTAMP = NT1AM(ISYRES) + NT2AM(ISYRES)
        END IF

        IF (LWORK .LT. NTAMP) THEN
          CALL QUIT('Insufficient memory in CC_DMAT (1).')
        END IF

        CALL DZERO(WORK,NTAMP)


        IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN

          IDTRAN(5,ITRAN) = IADRTH
          CALL PUTWA2(LUDMAT,FILDMA,WORK,IADRTH,NTAMP)
          IADRTH = IADRTH + NTAMP

        ELSE IF (IOPTRES.EQ.3) THEN

          CALL CC_WRRSP(FILDMA,IDTRAN(5,ITRAN),ISYRES,IOPTW,MODELW,
     &                  WORK,WORK,WORK,WORK(NTAMP),LWORK-NTAMP)

        ELSE IF (IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN

          CONTINUE

        ELSE
          WRITE (LUPRI,*) 'Error in CC_DMAT: illegal value for IOPTRES.'
          CALL QUIT('Error in CC_DMAT.')
        END IF

      END DO

*=====================================================================*

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
      IF (.NOT. (CCS .OR. CCSTST) ) THEN


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
           END IF
 
           KEND = KFREE
           LEND = LFREE

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
          END IF
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
     *                     WORK(KODCL1),WORK(KODCL2),
     *                     WORK(KODBC1),WORK(KODBC2),
     *                     WORK(KRDBC1),WORK(KRDBC2),
     *                     WORK(KODPP1),WORK(KODPP2),
     *                     WORK(KRDPP1),WORK(KRDPP2),
     *                     WORK(KCCFB1),WORK(KINDXB),
     *                     WORK(KEND), LEND,IPRERI)
            END IF

            KRECNR = KEND
            KEND   = KRECNR + (NBUFX(0) - 1)/IRAT + 1
            LEND   = LWORK - KEND
 
            IF (LEND .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_DMAT. (2)')
            END IF

          ELSE
            NUMDIS = 1
          END IF


*---------------------------------------------------------------------*
*         start loop over C matrix transformations:
*---------------------------------------------------------------------*
          DO ITRAN = 1, NDTRAN

            IDLSTA = IDTRAN(1,ITRAN)
            IDLSTB = IDTRAN(2,ITRAN)
            IDLSTC = IDTRAN(3,ITRAN)
            IDLSTD = IDTRAN(4,ITRAN)

            ISYMTA = ILSTSYM(LISTA,IDLSTA)
            ISYMTB = ILSTSYM(LISTB,IDLSTB)
            ISYMTC = ILSTSYM(LISTC,IDLSTC)
            ISYMTD = ILSTSYM(LISTD,IDLSTD)

            ISYRES =MULD2H(MULD2H(ISYMTA,ISYMTB),MULD2H(ISYMTC,ISYMTD))

* single excitation parts of coupled cluster vectors:
            KT1AMP0 = KEND
            KT1AMPA = KT1AMP0 + NT1AM(ISYM0)
            KT1AMPB = KT1AMPA + NT1AM(ISYMTA)
            KT1AMPC = KT1AMPB + NT1AM(ISYMTB)
            KT1AMPD = KT1AMPC + NT1AM(ISYMTC)
            KENDF1  = KT1AMPD + NT1AM(ISYMTD)

* Lambda-hole matrices:
            KLAMDH0 = KENDF1
            KLAMDHA = KLAMDH0 + NGLMDT(ISYM0)
            KLAMDHB = KLAMDHA + NGLMDT(ISYMTA)
            KLAMDHC = KLAMDHB + NGLMDT(ISYMTB)
            KLAMDHD = KLAMDHC + NGLMDT(ISYMTC)
            KENDF1  = KLAMDHD + NGLMDT(ISYMTD)

* Lambda-particle matrices:
            KLAMDP0 = KENDF1
            KLAMDPA = KLAMDP0 + NGLMDT(ISYM0)
            KLAMDPB = KLAMDPA + NGLMDT(ISYMTA)
            KLAMDPC = KLAMDPB + NGLMDT(ISYMTB)
            KLAMDPD = KLAMDPC + NGLMDT(ISYMTC)
            KENDF1  = KLAMDPD + NGLMDT(ISYMTD)

* the result vector:
            KTHETA1 = KENDF1
            KTHETA2 = KTHETA1 + NT1AM(ISYRES)
            KENDF1  = KTHETA2 + NT2AM(ISYRES)
            LENDF1  = LWORK - KENDF1

            IF (LENDF1 .LT. 0) THEN
              CALL QUIT('Insufficient memory in CC_DMAT. (3)')
            END IF

  
* read coupled cluster vectors:
            IOPTRW = 1
            CALL CC_RDRSP('R0',0,ISYM0,IOPTRW,MODEL,
     &                    WORK(KT1AMP0),WORK(KDUM))

            IOPTRW = 1
            CALL CC_RDRSP(LISTA,IDLSTA,ISYMTA,IOPTRW,MODEL,
     &                    WORK(KT1AMPA),WORK(KDUM))

            IOPTRW = 1
            CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPTRW,MODEL,
     &                    WORK(KT1AMPB),WORK(KDUM))

            IOPTRW = 1
            CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPTRW,MODEL,
     &                    WORK(KT1AMPC),WORK(KDUM))

            IOPTRW = 1
            CALL CC_RDRSP(LISTD,IDLSTD,ISYMTD,IOPTRW,MODEL,
     &                    WORK(KT1AMPD),WORK(KDUM))

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

            CALL CCLR_LAMTRA(WORK(KLAMDP0),WORK(KLAMDPD),WORK(KLAMDH0),
     &                       WORK(KLAMDHD),WORK(KT1AMPD),ISYMTD)


* read result vector:
            IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
              CALL GETWA2(LUDMAT,FILDMA,WORK(KTHETA1),
     &                    IDTRAN(5,ITRAN),NT1AM(ISYRES)+NT2AM(ISYRES))
            ELSE IF (IOPTRES.EQ.3) THEN
              IOPTRW = 3
              CALL CC_RDRSP(FILDMA,IDTRAN(5,ITRAN),ISYRES,IOPTRW,
     &                      MODEL,WORK(KTHETA1),WORK(KTHETA2))
            ELSE IF (IOPTRES.EQ.4) THEN
              CALL DZERO( WORK(KTHETA1), NT1AM(ISYRES) )
              CALL DZERO( WORK(KTHETA2), NT2AM(ISYRES) )
            ELSE IF (IOPTRES.EQ.5) THEN
              CALL DZERO( WORK(KTHETA1), NT1AM(ISYRES) )
              CALL DZERO( WORK(KTHETA2), NT2AM(ISYRES) )
            ELSE
              CALL QUIT('Error in CC_DMAT.')
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
              CALL QUIT('Insufficient work space in CC_DMAT. (4)')
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
     &                      WORK(KLAMDPC),WORK(KLAMDHD),
     &                      ISYMTA,ISYMTB,ISYMTC,ISYMTD,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPA),WORK(KLAMDHA),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPD),WORK(KLAMDHC),
     &                      ISYMTA,ISYMTB,ISYMTD,ISYMTC,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

*.....................................................................*
            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPC),WORK(KLAMDHC),
     &                      WORK(KLAMDPD),WORK(KLAMDHA),
     &                      ISYMTB,ISYMTC,ISYMTD,ISYMTA,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPD),WORK(KLAMDHD),
     &                      WORK(KLAMDPC),WORK(KLAMDHA),
     &                      ISYMTB,ISYMTD,ISYMTC,ISYMTA,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)

*.....................................................................*
            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPC),WORK(KLAMDHC),
     &                      WORK(KLAMDPA),WORK(KLAMDHD),
     &                      ISYMTB,ISYMTC,ISYMTA,ISYMTD,
     &                      WORK(KENDF2),LENDF2,IDEL,
     &                      ISYDEL,ISYRES,ISYM0,IOPT)
 
            CALL CC_MOFCON2(WORK(KXINT),WORK(KTHETA2),
     &                      WORK(KLAMDPB),WORK(KLAMDHB),
     &                      WORK(KLAMDPD),WORK(KLAMDHD),
     &                      WORK(KLAMDPA),WORK(KLAMDHC),
     &                      ISYMTB,ISYMTD,ISYMTA,ISYMTC,
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
            WRITE (LUPRI,*) MSGDBG, 'LISTA/IDLSTA:',LISTA,
     &           IDTRAN(1,ITRAN)
            WRITE (LUPRI,*) MSGDBG, 'LISTB/IDLSTB:',LISTB,
     &           IDTRAN(2,ITRAN)
            WRITE (LUPRI,*) MSGDBG, 'LISTC/IDLSTC:',LISTC,
     &           IDTRAN(3,ITRAN)
            WRITE (LUPRI,*) MSGDBG, 'LISTD/IDLSTD:',LISTD,
     &           IDTRAN(4,ITRAN)
            CALL CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYRES,1,1)
            WRITE (LUPRI,*) MSGDBG,'IOPTRES:',IOPTRES
            WRITE (LUPRI,*) MSGDBG,'FILDMA:',FILDMA
            WRITE (LUPRI,*) MSGDBG,'LUDMAT:',LUDMAT
            WRITE (LUPRI,*) MSGDBG,'IDTRAN(5,ITRAN):',IDTRAN(5,ITRAN)
            CALL FLSHFO(LUPRI)
C           IOPTRW = 3
C           CALL CC_RDRSP('L2',1,ISYRES,IOPTRW,MODEL,
C    &                    WORK(KENDF1),WORK(KENDF1+NT1AM(ISYRES)))
C           WRITE (LUPRI,*) MSGDBG,'L2 x THETA:',
C    &        DDOT(NT1AM(ISYRES)+NT2AM(ISYRES),WORK(KENDF1),1,
C    &                WORK(KTHETA1),1)
          END IF

          KTHETA0 = -999999

          IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
            CALL PUTWA2(LUDMAT,FILDMA,WORK(KTHETA1),
     &                  IDTRAN(5,ITRAN),NT1AM(ISYRES)+NT2AM(ISYRES))
c         ELSE IF (IOPTRES.EQ.3) THEN
c           CALL CC_WRRSP(FILDMA,IDTRAN(5,ITRAN),ISYRES,IOPTW,MODELW,
c    &                    WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
c    &                    WORK(KENDF1),LENDF1)
          ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
            CALL CC_WARSP(FILDMA,IDTRAN(5,ITRAN),ISYRES,IOPTW,MODELW,
     &                    WORK(KTHETA0),WORK(KTHETA1),WORK(KTHETA2),
     &                    WORK(KENDF1),LENDF1)
          ELSE IF (IOPTRES.EQ.5) THEN
            IOPTW = 2
            CALL CCDOTRSP(IDDOTS,DCONS,IOPTW,FILDMA,ITRAN,NDTRAN,MXVEC,
     &                    WORK(KTHETA1),WORK(KTHETA2),ISYRES,
     &                    WORK(KENDF1),LENDF1)
          ELSE
            CALL QUIT('Error in CC_DMAT.')
          END IF

*---------------------------------------------------------------------*
*         end of the loop over D matrix transformations:
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

      END IF
*=====================================================================*
* end of F term section
*=====================================================================*

*=====================================================================*
* restore result vectors and clean up and return:
*=====================================================================*
8888  CONTINUE

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

        CALL GETWA2(LUDMAT,FILDMA,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
          DO ITRAN = 1, NDTRAN
            IF (ITRAN.LT.NDTRAN) THEN
              LEN     = IDTRAN(5,ITRAN+1)-IDTRAN(5,ITRAN)
            ELSE
              LEN     = IADRTH-IDTRAN(5,NDTRAN)
            END IF
            KTHETA1 = IDTRAN(5,ITRAN)
            XNORM   = DDOT(LEN, WORK(KTHETA1),1, WORK(KTHETA1),1)
            WRITE (LUPRI,*) 'Read D matrix transformation nb. ',ITRAN
            WRITE (LUPRI,*) 'Adress, length, NORM:',IDTRAN(5,NDTRAN),
     &           LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF 

*---------------------------------------------------------------------*
* close D matrix file & return
*---------------------------------------------------------------------*
* check return option for the result vectors:
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN

        CALL WCLOSE2(LUDMAT,FILDMA,'KEEP')

      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_DMAT.')
      END IF


*=====================================================================*

      RETURN
      END
*=====================================================================*
*            END OF SUBROUTINE CC_DMAT
*=====================================================================*

*---------------------------------------------------------------------*
c/* Deck CC_DTST */
*=====================================================================*
       SUBROUTINE CC_DTST(WORK,LWORK)
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
      PARAMETER (MSGDBG='[debug] CC_DTST> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER MXCTRAN
      PARAMETER (MXCTRAN = 2)

      INTEGER LWORK
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL DDOT, RDUM, D1, D2
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION DDOT, RDUM(2), D1, D2
#endif 

      CHARACTER*(3) LISTA, LISTB, LISTC, LISTD, LISTL
      CHARACTER*(8) FILDMA, LABELA
      CHARACTER*(10) MODEL
      INTEGER IOPTRES
      INTEGER IDTRAN(5,MXCTRAN), NDTRAN
      INTEGER IDLSTA, IDLSTB, IDLSTC, IDLSTD, IDUM(2), IZETAV
      INTEGER ISYMA, ISYMB, ISYMC, ISYMD, ISYMABCD
      INTEGER KTHETA1, KTHETA2, KZETA1, KZETA2
      INTEGER KT1AMPD, KT2AMPD, KRESLT1, KRESLT2
      INTEGER KEND1, LEND1, IOPT, IOPTRW

* external function:
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER ILSTSYM



*---------------------------------------------------------------------*
* call C matrix transformation:
*---------------------------------------------------------------------*
      LISTL = 'L1'
      LISTA = 'R1'
      LISTB = 'R1'
      LISTC = 'R1'
      LISTD = 'R1'
C     IZETAV = IL2ZETA('ZDIPLEN ',0.0D0,1,'ZDIPLEN ',0.0D0,1)
      IZETAV = IL1ZETA('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTA = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTB = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTC = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTD = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDTRAN(1,1) = IDLSTA
      IDTRAN(2,1) = IDLSTB
      IDTRAN(3,1) = IDLSTC
      IDTRAN(4,1) = IDLSTD
      NDTRAN = 1

      IOPTRES = 1
      FILDMA  = 'CCDMAT'

      CALL CC_DMAT(IDTRAN,  NDTRAN,
     &             LISTA,  LISTB, LISTC, LISTD,
     &             IOPTRES, FILDMA, IDUM, RDUM, 0, WORK, LWORK )


      ISYMA  = ILSTSYM(LISTA,IDLSTA)
      ISYMB  = ILSTSYM(LISTB,IDLSTB)
      ISYMC  = ILSTSYM(LISTC,IDLSTC)
      ISYMD  = ILSTSYM(LISTD,IDLSTD)
      ISYMABCD = MULD2H(MULD2H(ISYMA,ISYMB),MULD2H(ISYMC,ISYMD))

      KTHETA1 = IDTRAN(5,1)
      KTHETA2 = KTHETA1 + NT1AM(ISYMABCD)

      IF (.TRUE.) THEN
        KZETA1  = KTHETA2 + NT2AM(ISYMABCD)
        KZETA2  = KZETA1  + NT1AM(ISYMABCD)
        KEND1   = KZETA2  + NT2AM(ISYMABCD)
        LEND1   = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_DTST.')
        END IF

        IOPTRW = 3
        Call CC_RDRSP(LISTL,IZETAV,ISYMABCD,IOPTRW,MODEL,
     &                WORK(KZETA1),WORK(KZETA2))

        D1=DDOT(NT1AM(ISYMABCD),WORK(KZETA1),1,WORK(KTHETA1),1)
        D2=DDOT(NT2AM(ISYMABCD),WORK(KZETA2),1,WORK(KTHETA2),1)

        WRITE (LUPRI,*) 'Dot product with left vector:',D1+D2
        WRITE (LUPRI,*) 'singles excitation part:', D1
        WRITE (LUPRI,*) 'double excitation part: ', D2
        WRITE (LUPRI,*) LISTL, IZETAV, ISYMABCD
        Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMABCD,1,1)
      END IF

      IF (NSYM.EQ.1 .AND. LOCDBG) THEN
        KT1AMPD = KTHETA2 + NT2AM(ISYMABCD)
        KT2AMPD = KT1AMPD + NT1AM(ISYMD)
        KRESLT1 = KT2AMPD + NT2AM(ISYMD)
        KRESLT2 = KRESLT1 + NT1AM(ISYMABCD)
        KEND1   = KRESLT2 + NT2AM(ISYMABCD)
        LEND1   = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_DTST.')
        END IF

        IOPTRW = 3
        Call CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPTRW,MODEL,
     &                WORK(KT1AMPD),WORK(KT2AMPD))

        ! zero singles or doubles C vector:
C       CALL DZERO(WORK(KT1AMPD),NT1AM(ISYMD))
C       CALL DZERO(WORK(KT2AMPD),NT2AM(ISYMD))
        CALL DZERO(WORK(KRESLT1),NT1AM(ISYMABCD)+NT2AM(ISYMABCD))
        IPRINT  = 5

        CALL CC_FDD(NT1AM(ISYMABCD),NT2AM(ISYMABCD),
     >              LISTA,IDLSTA,LISTB,IDLSTB,LISTC,IDLSTC,
     >              WORK(KT1AMPD), WORK(KRESLT1),
     >              WORK(KEND1), LEND1)

        IPRINT  = 0

        IF (.TRUE.) THEN
          WRITE (LUPRI,*) 'LISTA, IDLSTA, ISYMA:',LISTA,IDLSTA,ISYMA
          WRITE (LUPRI,*) 'LISTB, IDLSTB, ISYMB:',LISTB,IDLSTB,ISYMB
          WRITE (LUPRI,*) 'LISTC, IDLSTC, ISYMC:',LISTC,IDLSTC,ISYMC
          WRITE (LUPRI,*) 'LISTD, IDLSTD, ISYMD:',LISTD,IDLSTD,ISYMD
          WRITE (LUPRI,*) 'ISYMABCD:',ISYMABCD
          WRITE (LUPRI,*)
          WRITE (LUPRI,*) 'finite difference Theta vector:'
          Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMABCD,1,1)
          WRITE (LUPRI,*) 'analytical Theta vector:'
          Call CC_PRP(WORK(KTHETA1),WORK(KTHETA2),ISYMABCD,1,1)
        END IF

        Call DAXPY(NT1AM(ISYMABCD),-1.0d0,WORK(KTHETA1),1,
     &                                  WORK(KRESLT1),1)
        IF (.NOT.CCS) THEN
          Call DAXPY(NT2AM(ISYMABCD),-1.0d0,WORK(KTHETA2),1,
     &                                    WORK(KRESLT2),1)
        ELSE
          Call DZERO(WORK(KRESLT2),NT2AM(ISYMABCD))
        END IF

        WRITE (LUPRI,*) 'Norm of difference between analytical THETA '
     >           // 'vector and the numerical result:'
        WRITE (LUPRI,*) 'singles excitation part:',
     >   DSQRT(DDOT(NT1AM(ISYMABCD),WORK(KRESLT1),1,WORK(KRESLT1),1))
        WRITE (LUPRI,*) 'double excitation part: ',
     >   DSQRT(DDOT(NT2AM(ISYMABCD),WORK(KRESLT2),1,WORK(KRESLT2),1))

        WRITE (LUPRI,*) 'difference vector:'
        Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYMABCD,1,1)

        CALL FLSHFO(LUPRI)


      ELSE IF (NSYM.NE.1 .AND. LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_DTST> can not calculate finite '//
     &        'difference D matrix'
        WRITE (LUPRI,*) 'CC_DTST> with symmetry.'
      END IF

      RETURN
      END 
*=====================================================================*
