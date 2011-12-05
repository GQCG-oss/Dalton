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
      SUBROUTINE CCSDT_ETA_CONT(LISTL,LISTR,NODDY,
     &                          IDOTS,DOTPROD,
     &                          NXETRAN,IXETRAN,MXVEC,
     &                          NLEFT,IEASTEND,
     &                          NDENS,IDXLVEC,IDXRVEC,
     &                          WORK,LWORK,
     &                          FNDELD,FNCKJD,FNDKBC,FNTOC,
     &                          FN3VI,FNDKBC3,FN3FOPX,FN3FOP2X)
*---------------------------------------------------------------------*
*
*    Purpose: just a little driver to call the real routines...
*
*    Written by Christof Haettig, Mai 2002, Friedrichstal
*
*=====================================================================*
      IMPLICIT NONE
#include <ccsdsym.h>
#include <priunit.h>
#include <ccorb.h>
#include <dummy.h>
#include <cclists.h>
#include <ccroper.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*12 FNDEFF
      PARAMETER( FNDEFF = 'TMP-XIETADEN')

      INTEGER ISYM0
      PARAMETER( ISYM0 = 1 )

      CHARACTER*3 LISTL, LISTR
      LOGICAL NODDY
      INTEGER LWORK, MXVEC, NXETRAN, NLEFT, NDENS
      INTEGER IDOTS(MXVEC,NXETRAN), IEASTEND(2,NLEFT),
     &        IDXLVEC(NDENS), IDXRVEC(NDENS), 
     &        IXETRAN(MXDIM_XEVEC,NXETRAN)
      CHARACTER*(*) FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3
      CHARACTER*(*) FN3FOPX, FN3FOP2X

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NXETRAN), TRIPCON, CC_XIETA_DENCON
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NXETRAN), TRIPCON, CC_XIETA_DENCON
      DOUBLE PRECISION WORK(LWORK)
#endif

      CHARACTER MODEL*(10), LABELH*(8)
      LOGICAL SKIPXI, SKIPETA
      INTEGER LUDELD, LUCKJD, LUDKBC, LUTOC, LU3VI, LUDKBC3, LU3FOPX,
     *        LU3FOP2X, LUDEFF, LUFOCK,
     *        KLAMP0, KLAMH0, KFOCK0, KEND0, LWRK0, KT1AMP0, KEND1,
     *        LWRK1, IOPT, ILEFT, IDLSTL, IDENS, IDLSTR, ISYMR, ISYML,
     *        ISYDEN, KDAB, KDIJ, KDIA, IADRIJ, IADRAB, IADRIA, IOPER,
     *        ISYHOP, ISYCTR, ISYETA, IVEC, IDX, M2BST, ITRAN, ISYM
C

      INTEGER ILSTSYM

*---------------------------------------------------------------------*
*     begin
*---------------------------------------------------------------------*
      CALL QENTER('CCSDT_ETA_CONT')

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'entered CCSDT_ETA_CONT... NODDY=',NODDY
        WRITE(LUPRI,*) 'LISTL,LISTR:',LISTL,LISTR
      END IF

      LUDEFF  = -1
      CALL WOPEN2(LUDEFF,FNDEFF,64,0)

      M2BST = 0
      DO ISYM = 1, NSYM
        M2BST = MAX(M2BST,N2BST(ISYM))
      END DO
*---------------------------------------------------------------------*
*     initialize 0.th-order Lambda and Fock matrices:
*---------------------------------------------------------------------*
      KLAMP0  = 1
      KLAMH0  = KLAMP0  + NLAMDT
      KFOCK0  = KLAMH0  + NLAMDT
      KEND0   = KFOCK0  + N2BST(ISYM0)
      LWRK0   = LWORK   - KEND0

      KT1AMP0 = KEND0
      KEND1   = KT1AMP0 + NT1AMX
      LWRK1   = LWORK   - KEND1

      IF (LWRK1.LT.0) CALL QUIT('Out of memory in CCSDT_ETA_CONT')

      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),DUMMY)

      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*     read zeroth-order AO Fock matrix from file: 
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,1,1,1)

*---------------------------------------------------------------------*
*     precompute effective densities and store on file:
*---------------------------------------------------------------------*
      LUDELD  = -1
      LUCKJD  = -1
      LUDKBC  = -1
      LUTOC   = -1
      LU3VI   = -1
      LUDKBC3 = -1
      LU3FOPX = -1
      LU3FOP2X= -1
 
      CALL WOPEN2(LUDELD,FNDELD,64,0)
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
      CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
      CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)

      DO ILEFT = 1, NLEFT
        IDLSTL = IDXLVEC(IEASTEND(1,ILEFT))

        DO IDENS = IEASTEND(1,ILEFT), IEASTEND(2,ILEFT)

         IF (IDXLVEC(IDENS).NE.IDLSTL) 
     *     CALL QUIT('Loop error in CCSDT_ETA_CONT.')

         IDLSTR = IDXRVEC(IDENS)
         ISYMR  = ILSTSYM(LISTR,IDLSTR)
         ISYML  = ILSTSYM(LISTL,IDLSTL)
         ISYDEN = MULD2H(ISYMR,ISYML)
         
         KDAB   = KEND0
         KDIJ   = KDAB   + NMATAB(ISYDEN)
         KDIA   = KDIJ   + NMATIJ(ISYDEN)
         KEND1  = KDIA   + NT1AM(ISYDEN)
         LWRK1  = LWORK  - KEND1

         IF (LWRK1.LT.0) CALL QUIT('Out of memory in CCSDT_ETA_CONT')

         IF (LISTL.EQ.'L0 ') THEN
C         --------------------
C         Compute the density:
C         --------------------
          CALL CCSDT_ETA_DEN(LISTL,IDLSTL,LISTR,IDLSTR,
     *                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     *                       WORK(KDIJ),WORK(KDAB),WORK(KDIA),
     *                       WORK(KEND1),LWRK1,
     *                       LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                       FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                       LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                       LU3FOP2X,FN3FOP2X)
         ELSE

           IF (NODDY) THEN
             CALL CCSDT_ADEN_NODDY(LISTL,IDLSTL,LISTR,IDLSTR,
     *                       WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     *                       WORK(KDIJ),WORK(KDAB),.TRUE.,WORK(KDIA),
     *                       .FALSE.,DUMMY,WORK(KEND1),LWRK1,
     *                       LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                       FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                       LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                       LU3FOP2X,FN3FOP2X)
           ELSE
              IF (LISTR(1:3).EQ.'R1 '.OR.LISTR(1:3).EQ.'RE ') THEN
                CALL CC3_ADEN(LISTL,IDLSTL,LISTR,IDLSTR,
     *                         WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     *                         WORK(KDIJ),WORK(KDAB),
     *                         .TRUE.,WORK(KDIA),ISYDEN,
     *                         .FALSE.,DUMMY,
     *                         WORK(KEND1),LWRK1,
     *                         LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                         FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                         LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                         LU3FOP2X,FN3FOP2X)
              ELSEIF (LISTR(1:3).EQ.'R2 '.OR.LISTR(1:3).EQ.'ER1') THEN
                 CALL CC3_ADEN_CUB(LISTL,IDLSTL,LISTR,IDLSTR,
     *                         WORK(KLAMP0),WORK(KLAMH0),WORK(KFOCK0),
     *                         WORK(KDIJ),WORK(KDAB),WORK(KDIA),ISYDEN,
     *                         WORK(KEND1),LWRK1,
     *                         LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                         FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                         LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                         LU3FOP2X,FN3FOP2X)
              ELSE
                 CALL QUIT('Unknown LISTR in CCSDT_ETA_CONT.')
              ENDIF

           END IF

         END IF


         IADRIJ = 1 + M2BST*(IDENS-1)
         CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDIJ),IADRIJ,NMATIJ(ISYDEN))

         IADRAB = IADRIJ + NMATIJ(ISYDEN)
         CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDAB),IADRAB,NMATAB(ISYDEN))

         IADRIA = IADRAB + NMATAB(ISYDEN)
         CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDIA),IADRIA,NT1AM(ISYDEN))

        END DO ! IDENS

      END DO ! ILEFT
         
      CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
      CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
      CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')

*---------------------------------------------------------------------*
*     compute from the densities the polarizability contributions:
*---------------------------------------------------------------------*
      DO ITRAN = 1, NXETRAN
        IOPER  = IXETRAN(1,ITRAN)  ! operator index
        IDLSTL = IXETRAN(2,ITRAN)  ! Left vector index

        ISYHOP = ISYOPR(IOPER)     ! symmetry of O operator
        LABELH = LBLOPR(IOPER)     ! operator label

        SKIPXI  = ( IXETRAN(3,ITRAN) .EQ. -1 )
        SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )

        IF (.NOT.SKIPETA) THEN
          ISYCTR = ILSTSYM(LISTL,IDLSTL) ! sym. of Left (ZETA) vector
          ISYETA = MULD2H(ISYHOP,ISYCTR) ! sym. of ETA result vect.

          IVEC = 1
          DO WHILE(IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            IDLSTR = IDOTS(IVEC,ITRAN)
            ISYMR  = ILSTSYM(LISTR,IDLSTR)

            IF (ISYMR.EQ.ISYETA) THEN
              ! find index of density
              IDENS = 0
              DO IDX = 1, NDENS
                IF (IDLSTL.EQ.IDXLVEC(IDX) .AND.
     &              IDLSTR.EQ.IDXRVEC(IDX)       ) THEN
                  IDENS = IDX
                END IF 
              END DO
              
              IF (IDENS.EQ.0) 
     &          CALL QUIT('Density not found in ccsdt_eta_cont.')
     
              TRIPCON = CC_XIETA_DENCON(IDENS,LABELH,ISYHOP,
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     LUDEFF,FNDEFF,M2BST,WORK(KEND0),LWRK0)

              IF (LOCDBG) THEN
                WRITE(LUPRI,*) 'IVEC,ITRAN,TRIPCON:',
     &                          IVEC,ITRAN,TRIPCON
                WRITE(LUPRI,*) 'DOTPROD before,after:',
     &            DOTPROD(IVEC,ITRAN),DOTPROD(IVEC,ITRAN)+TRIPCON
              END IF

              DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TRIPCON

            END IF 

            IVEC = IVEC + 1
          END DO 

        END IF
      END DO ! ITRAN


      CALL WCLOSE2(LUDEFF,FNDEFF,'DELETE')

      CALL QEXIT('CCSDT_ETA_CONT')
      RETURN
      END 
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCSDT_ETA_DEN(LISTL,IDLSTL,LISTR,IDLSTR,
     *                         XLAMDP,XLAMDH,FOCK0,
     *                         DENIJ,DENAB,DENIA,
     *                         WORK,LWORK,
     *                         LUDELD,FNDELD,LUCKJD, FNCKJD,
     *                         LUDKBC,FNDKBC,LUTOC,  FNTOC,
     *                         LU3VI, FN3VI, LUDKBC3,FNDKBC3,
     *                         LU3FOP,FN3FOP,LU3FOP2,FN3FOP2)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to a contraction of 
*             an Eta vector with some 'right' vector 
*             
*             Calculation is done via a Eta-density:
*
*             Eta(X) x T^Y = sum_pq D^eta_pq X_pq
*
*    Written by Christof Haettig, Mai 2002, Aarhus
*
***********************************************************************
*
*    Spring 2004, Filip Pawlowski, Aarhus: Modified to treat also the 
*                                          Cauchy moments.
*            
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccr1rsp.h>
#include <ccexci.h>
#include <ccrc1rsp.h>
C
      LOGICAL LOCDBG
      LOGICAL LCAUCHY
      PARAMETER(LOCDBG = .FALSE.)
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
C
      CHARACTER*10 FNT2Y
      PARAMETER(FNT2Y = 'CCT2Y__TMP' )
C
      CHARACTER FNDPTIJ*5, FNDPTAB*5
      PARAMETER(FNDPTIJ = 'DPTIJ', FNDPTAB='DPTAB')
C
      CHARACTER LISTL*3, LISTR*3
      CHARACTER*(*) FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3
      CHARACTER*(*) FN3FOP, FN3FOP2
      
      INTEGER LUDELD, LUCKJD, LUDKBC, LUTOC, LU3VI, LUDKBC3, LU3FOP,
     *        LU3FOP2, IDLSTR, IDLSTL

      INTEGER ISYCTR, ISTAMY, KT2Y, ISYDEN, KDAB, KDIJ, KDIA1,
     *        KDIA2, KDIA3, KDIA4, KMMAT, KYMMAT, KT2TP, KL1AM, KL2TP
      INTEGER LUPTIJ, LUPTAB, LUT2Y

      CHARACTER LABELY*8, MODEL*10, CDUMMY*1
      INTEGER ISYM, ILEN, ISYMIM, NIMFN(8), IOPT, KEND0, LWRK0,
     &        KD0AB, KD0IJ, KT1AMY, KEND2, LWRK2, IOPTT2, ISYMA,
     &        ISYMI, ISYMB, ISYMJ, KOFFDBA, KOFFDIJ, KOFFDIA, KOFFTBI,
     &        KOFFTAJ, NVIRA, NVIRB, NRHFI, ISYMFN, LWORK
C
      INTEGER NCAU

#if defined (SYS_CRAY)
      REAL XLAMDP(*), XLAMDH(*), FOCK0(*)
      REAL PRECISION FREQY, WORK(*), DDOT
      REAL ZERO, ONE, TWO, HALF
      REAL DENIJ(*),DENAB(*),DENIA(*)
#else
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), FOCK0(*)
      DOUBLE PRECISION FREQY, WORK(*), DDOT
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      DOUBLE PRECISION DENIJ(*),DENAB(*),DENIA(*)
#endif
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INTEGER ILSTSYM
C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
      CALL QENTER('CTETADEN')
 
      IF (LISTL(1:3).EQ.'L0 ') THEN
        ISYCTR = ISYM0
      ELSE
        ! ups, probably higher-order response, not yet implemented here
        CALL QUIT('Unknown type of left vector in CCSDT_ETA_DEN.')
      END IF
C
      !initialize LCAUCHY and NCAU
      LCAUCHY = .FALSE.
      NCAU = 0

      IF (LISTR(1:3).EQ.'R1 ') THEN
        ISTAMY = ISYLRT(IDLSTR)
        FREQY  = FRQLRT(IDLSTR)
        LABELY = LRTLBL(IDLSTR)

        IF (LORXLRT(IDLSTR)) 
     &   CALL QUIT('NO ORBITAL RELAXED PERTURBATION IN CCSDT_ETA_DEN')

      ELSE IF (LISTR(1:3).EQ.'RE ') THEN
        ISTAMY  = ILSTSYM(LISTR,IDLSTR)
        FREQY  = EIGVAL(IDLSTR)
        LABELY = '- none -'

      ELSE IF (LISTR(1:3).EQ.'RC ') THEN
        ISTAMY  = ISYLRC(IDLSTR)
        FREQY   = ZERO
        LABELY  = LRCLBL(IDLSTR)
        NCAU    = ILRCAU(IDLSTR)
        LCAUCHY = .TRUE.

      ELSE
        ! ups, probably higher-order response, not yet implemented here
        WRITE(LUPRI,*)'LISTR = ',LISTR
        CALL QUIT('Unknown type of right vector in CCSDT_ETA_DEN.')
      END IF
C
C------------------------------------------------------------
C     compute dimensions of M matrices NIMFN:
C------------------------------------------------------------

      DO ISYM = 1, NSYM
        ILEN = 0
        DO ISYMFN = 1, NSYM
          ISYMIM = MULD2H(ISYM,ISYMFN)
          ILEN   = ILEN + NMATIJ(ISYMIM)*NT1AM(ISYMFN)
        END DO
        NIMFN(ISYM) = ILEN
      END DO
C
C------------------------------------------------------------------
C     resort response amplitudes and dump them on a temporary file:
C------------------------------------------------------------------
C
      KT2Y  = 1
      KEND0 = KT2Y  + NT2SQ(ISTAMY)
      LWRK0 = LWORK - KEND0

      IF (LWRK0 .LT. NT2AM(ISTAMY)) THEN
       CALL QUIT('Not enough memory for squared T2Y (CCSDT_ETA_DEN)')
      ENDIF

      IOPT = 2
      CALL CC_RDRSP(LISTR,IDLSTR,ISTAMY,IOPT,MODEL,DUMMY,WORK(KEND0))

      Call CCLR_DIASCL(WORK(KEND0),TWO,ISTAMY)

      CALL CC_T2SQ(WORK(KEND0),WORK(KT2Y),ISTAMY)
 
      LUT2Y = -1
      CALL WOPEN2(LUT2Y,FNT2Y,64,0)
      CALL PUTWA2(LUT2Y,FNT2Y,WORK(KT2Y),1,NT2SQ(ISTAMY))
C
C-------------------------------------------------------
C     initial allocations, prepare T2 and L2 amplitudes:
C-------------------------------------------------------
C
      ISYDEN = MULD2H(ISYCTR,ISTAMY)

      KEND0  = 1

      KDAB   = KEND0
      KDIJ   = KDAB   + NMATAB(ISYDEN)
      KDIA1  = KDIJ   + NMATIJ(ISYDEN)
      KDIA2  = KDIA1  + NT1AM(ISYDEN)
      KDIA3  = KDIA2  + NT1AM(ISYDEN)
      KDIA4  = KDIA3  + NT1AM(ISYDEN)
      KEND0  = KDIA4  + NT1AM(ISYDEN)

      KMMAT  = KEND0
      KYMMAT = KMMAT  + NIMFN(ISYCTR)
      KEND0  = KYMMAT + NIMFN(ISYDEN)

      KT2TP  = KEND0
      KL1AM  = KT2TP + NT2SQ(ISYM0)
      KL2TP  = KL1AM + NT1AM(ISYCTR)
      KEND0  = KL2TP + NT2SQ(ISYCTR)

      LWRK0  = LWORK - KEND0
 
      IF (LWRK0.LT.NT2SQ(ISYM0) .OR. LWRK0.LT.NT2SQ(ISYCTR)) THEN
       CALL QUIT('Not enough memory to construct T2TP (CCSDT_ETA_DEN)')
      ENDIF

 
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KT2TP))

      CALL CC_T2SQ(WORK(KT2TP),WORK(KEND0),ISYM0)

      CALL CC3_T2TP(WORK(KT2TP),WORK(KEND0),ISYM0)
 
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
     *    DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)
 
      IOPT = 3
      CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KL2TP))

      CALL CC_T2SQ(WORK(KL2TP),WORK(KEND0),ISYCTR)

      CALL CC3_T2TP(WORK(KL2TP),WORK(KEND0),ISYCTR)

      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
     *    DDOT(NT2SQ(ISYCTR),WORK(KL2TP),1,WORK(KL2TP),1)

      ! initialize density matrices
      CALL DZERO(WORK(KDAB), NMATAB(ISYDEN))
      CALL DZERO(WORK(KDIJ), NMATIJ(ISYDEN))
      CALL DZERO(WORK(KDIA1),NT1AM(ISYDEN))
      CALL DZERO(WORK(KDIA2),NT1AM(ISYDEN))
      CALL DZERO(WORK(KDIA3),NT1AM(ISYDEN))
      CALL DZERO(WORK(KDIA4),NT1AM(ISYDEN))

      ! initialize M-matrix intermediates
      CALL DZERO(WORK(KMMAT), NIMFN(ISYCTR))
      CALL DZERO(WORK(KYMMAT),NIMFN(ISYDEN))

      !frequency of LISTL is handled inside
      CALL CCSDT_ETA_FA_DEN(LISTL,IDLSTL,ISYCTR,
     *                      LISTR,IDLSTR,ISTAMY,FREQY,LABELY,
     *                      LCAUCHY,NCAU,
     *                      XLAMDP,XLAMDH,FOCK0,
     *                      WORK(KT2TP),WORK(KL2TP),WORK(KL1AM),
     *                      ISYDEN,WORK(KDAB),WORK(KDIJ),
     *                      .TRUE.,WORK(KDIA3),
     *                      .TRUE.,WORK(KMMAT),.TRUE.,WORK(KYMMAT),
     *                      WORK(KEND0),LWRK0,
     *                      LUDELD,FNDELD,LUCKJD, FNCKJD,
     *                      LUDKBC,FNDKBC,LUTOC,  FNTOC,
     *                      LU3VI, FN3VI, LUDKBC3,FNDKBC3,
     *                      LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                      LUT2Y,FNT2Y)

C     ---------------------------------------------
C     D(ia) <-- D(ia) + sum_fnm y^t(am,fn) M(imfn):
C     ---------------------------------------------
      IOPTT2 = 1
      CALL CCSDT_ETA_TM2(WORK(KDIA2),ISYDEN,WORK(KMMAT),ISYCTR,
     &                   DUMMY,LUT2Y,FNT2Y,ISTAMY,IOPTT2,
     &                   WORK(KEND0),LWRK0)

C     ---------------------------------------------
C     D(ia) <-- D(ia) + sum_fnm t(am,fn) y^M(imfn):
C     ---------------------------------------------
      IOPTT2 = 0
      CALL CCSDT_ETA_TM2(WORK(KDIA4),ISYDEN,WORK(KYMMAT),ISYDEN,
     &                   WORK(KT2TP),IDUMMY,CDUMMY,ISYM0,IOPTT2,
     &                   WORK(KEND0),LWRK0)

C
C     ----------------------------------------------------------------
C     D(ia) <-- D(ia) - sum_b y^t(bi) D^0(ba) + sum_j y^t(aj) D^0(ij):
C     ----------------------------------------------------------------
      ! read D^0(ab) 
      KD0AB  = KEND0
      KD0IJ  = KD0AB  + NMATAB(ISYM0)
      KT1AMY = KD0IJ  + NMATIJ(ISYM0)
      KEND2  = KT1AMY + NT1AM(ISTAMY)
      LWRK2  = LWORK  - KEND2
 
      IF (LWRK2 .LT. 0) THEN
       CALL QUIT('Out of memory in CCSDT_ETA_DEN.')
      ENDIF

      IOPT = 1
      CALL CC_RDRSP(LISTR,IDLSTR,ISTAMY,IOPT,MODEL,WORK(KT1AMY),DUMMY)

      LUPTIJ = -1
      LUPTAB = -1

      ! read intermediates from ground state density from file...
      CALL WOPEN2( LUPTIJ,FNDPTIJ,64,0)
      CALL GETWA2( LUPTIJ,FNDPTIJ,WORK(KD0IJ),1,NMATIJ(ISYM0))
      CALL WCLOSE2(LUPTIJ,FNDPTIJ,'KEEP')
      
      CALL WOPEN2( LUPTAB,FNDPTAB,64,0)
      CALL GETWA2( LUPTAB,FNDPTAB,WORK(KD0AB),1,NMATAB(ISYM0))
      CALL WCLOSE2(LUPTAB,FNDPTAB,'KEEP')
      

      DO ISYMA = 1, NSYM
        ISYMI = MULD2H(ISYDEN,ISYMA)
        ISYMB = MULD2H(ISYCTR,ISYMA)
        ISYMJ = MULD2H(ISYCTR,ISYMI)

        KOFFDBA = KD0AB + IMATAB(ISYMB,ISYMA)
        KOFFDIJ = KD0IJ + IMATIJ(ISYMI,ISYMJ)
        KOFFDIA = KDIA1 + IT1AM(ISYMA,ISYMI)
        KOFFTBI = KT1AMY+ IT1AM(ISYMB,ISYMI)
        KOFFTAJ = KT1AMY+ IT1AM(ISYMA,ISYMJ)

        NVIRA   = MAX(NVIR(ISYMA),1)
        NVIRB   = MAX(NVIR(ISYMB),1)
        NRHFI   = MAX(NRHF(ISYMI),1)

        CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &             -ONE,WORK(KOFFDBA),NVIRB,WORK(KOFFTBI),NVIRB,
     &              ONE,WORK(KOFFDIA),NVIRA)

        CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              ONE,WORK(KOFFTAJ),NVIRA,WORK(KOFFDIJ),NRHFI,
     &              ONE,WORK(KOFFDIA),NVIRA)

      END DO
C
C-----------------------------------------
C     copy/add densities to result arrays:
C-----------------------------------------
C
      CALL DCOPY(NMATIJ(ISYDEN),WORK(KDIJ),1,DENIJ,1)
      CALL DCOPY(NMATAB(ISYDEN),WORK(KDAB),1,DENAB,1)

      CALL DCOPY(NT1AM(ISYDEN),WORK(KDIA1),1,DENIA,1)
      CALL DAXPY(NT1AM(ISYDEN),+1.0D0,WORK(KDIA2),1,DENIA,1)
      CALL DAXPY(NT1AM(ISYDEN),-1.0D0,WORK(KDIA3),1,DENIA,1)
      CALL DAXPY(NT1AM(ISYDEN),+1.0D0,WORK(KDIA4),1,DENIA,1)


C-----------------------------------------
C     close/deleta scratch files:
C-----------------------------------------
      CALL WCLOSE2(LUT2Y,FNT2Y,'DELETE')
C
      CALL QEXIT('CTETADEN')
      RETURN
      END 
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCSDT_ETA_FA_DEN(LISTL,IDLSTL,ISYMBL,
     *                            LISTR,IDLSTR,ISTAMY,FREQY,LABELY,
     *                            LCAUCHY,NCAU,
     *                            XLAMDP,XLAMDH,FOCK0,
     *                            T2TP,ZL2TP,ZL1AM,
     *                            ISYDEN,DAB,DIJ,DO_DIA,DIA,
     *                            DO_MMAT,XMMAT,DO_YMMAT,YMMAT,
     *                            WORK,LWORK,
     *                            LUDELD,FNDELD,LUCKJD, FNCKJD,
     *                            LUDKBC,FNDKBC,LUTOC,  FNTOC,
     *                            LU3VI, FN3VI, LUDKBC3,FNDKBC3,
     *                            LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                            LUT2Y,FNT2Y)
*---------------------------------------------------------------------*
*
*    Purpose: compute contractions needed for triples contributions
*             to the Eta and F{A} densities
*             
*
*    Written by Christof Haettig, Mai 2002, Aarhus
*
***********************************************************************
*
*     Spring 2004, Filip Pawlowski, Aarhus: Modified to treat also the 
*                                           Cauchy moments.
*            
*=====================================================================*
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccinftap.h>
#include <inftap.h>
#include <ccexci.h>
C
      LOGICAL LOCDBG, TOT_T3Y
      LOGICAL LCAUCHY
      PARAMETER(LOCDBG = .FALSE., TOT_T3Y = .TRUE.)
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
C
      CHARACTER*10 FNTBAR, FNWMAT
      PARAMETER(FNTBAR = 'CCTBAR_TMP', FNWMAT = 'CCWMAT_TMP')
C
      CHARACTER LISTL*3, LISTR*3
      CHARACTER*12 FN3SRTR, FNDELDR
      CHARACTER*16 FNCKJDR, FNDKBCR
      CHARACTER*(*) FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3
      CHARACTER*(*) FN3FOP, FN3FOP2, FNT2Y
C
      PARAMETER(FN3SRTR  = 'CCSDT_ETAFD1', FNDELDR  = 'CCSDT_ETAFD3')

      CHARACTER MODEL*10, LABELY*8, CDUMMY*1

      LOGICAL DO_DIA, DO_MMAT, DO_YMMAT
 
      INTEGER LWORK, IDLSTL, IDLSTR

      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR
      INTEGER KFOCKD, KEND0, LWRK0,
     *        KTROC0, KXIAJB, KEND1, LWRK1, KINTOC, KEND2, LWRK2,
     *        LENGTH, IOFF, ISYMD, ISYCKB, ISTAMY,
     *        KTRVI1, KTRVI2, KTRVI0,
     *        KTRVI3, KEND3, LWRK3, KINTVI, KEND4, LWRK4, 
     *        ISALJB, ISYMBD, ISCKIJ, KBLSMATBD, KSMATBD,
     *        KDIAG, KDIAGW, ISYMC, ISYMK, KOFF1, KOFF2, KOFF3,
     *        KINDSQ, KINDEXB, KTMAT, LENSQ, KINDSQW,
     *        KFCKBA, KTRVI4, KTRVI5, KTRVI6,
     *        KUMATBD,  KBLUMATBD, 
     *        KTROC02, KTRVI7, KTRVI8, KTRVI9, KTRVI10,
     *        KTRVI11, KTRVI12, KTRVI13,
     *        ISYCKD, KSMATDB, KUMATDB, KEND5, LWRK5,
     *        ISALJD, KINDEXD, KBLSMATDB, KBLUMATDB,
     *        KTRVI14, KTRVI15, KTRVI16, KTRVI17, KTRVI18, KTRVI19,
     *        KTRVI20, ISYTMP, KTROC01, KTROC21, KTROC03, KTROC23,
     *        LUDELD, LUCKJD, LUDKBC, LUTOC, LU3VI, LUDKBC3, LU3FOP,
     *        LU3FOP2, LUTBAR, LUWMAT, LUT2Y, 
     *        ISYMAI, ISYML, ISYAIL, ISAIBJ, ISYMDL, ISYMBJ,
     *        KOFFA, IOPT, NBJ, IADR, ISYMJ, ISYMB, ISYMI, NRHFI,
     *        ISYMBL, ISAIB, ISYM, ILEN, ISYMFN, ISYMIM, ISYDEN,
     *        ISWMAT, KFOCKY, IRREP, NVIRB,
     *        IERR, KWMAT, LENSQW, ISTBAR, NDL, ISYMN, ISYEMF, KTBAR,
     *        ISYMBN, ISYMEM, ISYMF, KT2Y, KOFFD, NTOTEM, NNVIRF,
     *        ISYMFI, NNEMF, ISYMM, ISYME, ISYMEI, ISYEIL, 
     *        KFN, KOFFM, NVIRE, NRFHI, ISYMDN, KOFFT, KBN, KOFFW, 
     *        ISYEMN, ISYMEN, ISYENL, ISYEML,
     *        KT1B, KT2B, KTROC0Y, KTRVI0Y,
     *        IMAT, ISCKBY, KTRVI8Y,
     *        ISCKDY
      INTEGER ISWTL(8,8), ISTLN(8,8), NIMFN(8)

      INTEGER ISBLALJB,ISBLALJD,ISBLCKIJ,KBLINDEXB,KBLINDEXD
      INTEGER KBLDIAG,LENSQBL,KBLINDSQ
C
      INTEGER NCAU,MCAU,IDLSTRM,KOFFOCC,KOFFINTD,KOFFINTB
C
      !functions:
      INTEGER ILSTSYM,ILRCAMP
C
      integer kx3am
 
#if defined (SYS_CRAY)
      REAL XLAMDP(*), XLAMDH(*), FOCK0(*)
      REAL WORK(LWORK)
      REAL FREQY, DDOT
      REAL ZERO, ONE, TWO, HALF
      REAL DAB(*),DIJ(*),DIA(*)
      REAL XMMAT(*),YMMAT(*)
      REAL T2TP(*),ZL2TP(*),ZL1AM(*)
      REAL FREQL
#else
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), FOCK0(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQY, DDOT
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      DOUBLE PRECISION DAB(*),DIJ(*),DIA(*)
      DOUBLE PRECISION XMMAT(*),YMMAT(*)
      DOUBLE PRECISION T2TP(*),ZL2TP(*),ZL1AM(*)
      DOUBLE PRECISION FREQL
#endif
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
C------------------------------------------------------------
C     some initializations:
C     set offset arrays ISWTL and ISTLN and dimensions NIMFN:
C------------------------------------------------------------
      CALL QENTER('CTETAFAD')
C
      !read frequency for listl
      IF (LISTL(1:3) .EQ. 'L0 ') THEN
         FREQL = 0.0D0
      ELSE IF (LISTL(1:3) .EQ. 'LE ') THEN
         FREQL = EIGVAL(IDLSTL)
      ELSE
        WRITE(LUPRI,*) 'LISTL = ', LISTL(1:3) 
        WRITE(LUPRI,*) 'CCSDT_ETA_FA_DEN handles only L0 or LE as LISTL'
        CALL QUIT('Unknown left list in CCSDT_ETA_FA_DEN')
      END IF
C
      DO ISYMDL = 1, NSYM
        IOFF = 0        
        DO ISYML = 1, NSYM
          ISWMAT = MULD2H(ISYMDL,ISYML)
          ISWTL(ISWMAT,ISYML) = IOFF
          IOFF = IOFF + NT2SQ(ISWMAT)*NRHF(ISYML)
        END DO 
      END DO 

      DO ISAIBJ = 1, NSYM
        IOFF = 0
        DO ISYMJ = 1, NSYM
          ISAIB = MULD2H(ISAIBJ,ISYMJ)
          ISTLN(ISAIB,ISYMJ) = IOFF
          IOFF = IOFF + NCKATR(ISAIB)*NRHF(ISYMJ)
        END DO
      END DO

      DO ISYM = 1, NSYM
        ILEN = 0
        DO ISYMFN = 1, NSYM
          ISYMIM = MULD2H(ISYM,ISYMFN)
          ILEN   = ILEN + NMATIJ(ISYMIM)*NT1AM(ISYMFN)
        END DO
        NIMFN(ISYM) = ILEN
      END DO
C
C-------------------------------------------------------
C     initial allocations, prepare T2B amplitudes:
C-------------------------------------------------------
C
      KEND0  = 1

      IF ( (LISTR(1:3).EQ.'R1 ') .OR. (LISTR(1:3).EQ.'RC ') ) THEN
        KFOCKY = KEND0  
        KEND0  = KFOCKY + N2BST(ISTAMY)
      END IF

      KFOCKD = KEND0
      KFCKBA = KFOCKD + NORBTS
      KEND0  = KFCKBA + NT1AMX

      IF (TOT_T3Y) THEN
         KT1B   = KEND0
         KT2B   = KT1B + (NCAU+1)*NT1AM(ISTAMY)
         KEND0  = KT2B + (NCAU+1)*NT2SQ(ISTAMY)
      END IF

      LWRK0  = LWORK - KEND0
 
      IF (LWRK0.LT.0) THEN
       CALL QUIT('Not enough memory in CCSDT_ETA_FA_DEN')
      ENDIF

C
C     ---------------------------------------------------------------------
C     Get KT1B and KT2B vectors
C
C     For non-Cacuhy calculations (LCAUCHY = .FALSE., NCAU = 0) the
C     loop below is dúmmy.
C
C     For Cacuhy calculations (LCAUCHY = .TRUE., NCAU >= 0) 
C     KT1B and KT2B correspond to singles and doubles part of Cauchy
C     vectors.
C     ---------------------------------------------------------------------
C
      IF (TOT_T3Y) THEN
C
         IF (LWRK0.LT.NT2SQ(ISTAMY)) THEN
          CALL QUIT('Not enough memory for T2Y (CCSDT_ETA_FA_DEN)')
         ENDIF
C
         DO MCAU = 0, NCAU
C
            IF (LCAUCHY) THEN
              !Get the list for current MCAU
              IDLSTRM = ILRCAMP(LABELY,MCAU,ISTAMY)
              !Consistency check
              IF (MCAU.EQ.NCAU .AND. IDLSTRM.NE.IDLSTR) THEN
               CALL QUIT('Sth wrong in Cauchy loop in CCSDT_ETA_FA_DEN')
              END IF
            ELSE
              IDLSTRM = IDLSTR
            END IF
C
            IOPT = 3
            KOFF1 = KT1B + MCAU*NT1AM(ISTAMY)
            KOFF2 = KT2B + MCAU*NT2SQ(ISTAMY)
            CALL CC_RDRSP(LISTR,IDLSTRM,ISTAMY,IOPT,MODEL,
     *                    WORK(KOFF1),WORK(KOFF2))
            CALL CCLR_DIASCL(WORK(KOFF2),TWO,ISTAMY)
            CALL CC_T2SQ(WORK(KOFF2),WORK(KEND0),ISTAMY)
            CALL CC3_T2TP(WORK(KOFF2),WORK(KEND0),ISTAMY)
C
         END DO !MCAU
C
      END IF
C
C     ----------------------------------------------------------------------
C     Calculate (ck|de)-tilde(Y) and (ck|lm)-tilde(Y) needed in construction
C     of W3 ( <mu3|[[H^0,T1B],T2^0]|HF> term).
C     ----------------------------------------------------------------------
C
      IF (TOT_T3Y) THEN
C
         DO MCAU = 0, NCAU
C
            !Open temporary files
            LU3SRTR  = -1
            LUDELDR  = -1
            CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
            CALL WOPEN2(LUDELDR,FNDELDR,64,0)
C
            !Generate file names for T1-transformed integrals
            IF (LCAUCHY) THEN
               CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
            ELSE
               FNCKJDR  = 'CCSDT_____ETAFD2'
               FNDKBCR  = 'CCSDT_____ETAFD4'
            END IF
C
            !Open the files for T1-transformed integrals
            LUCKJDR  = -1
            LUDKBCR  = -1
            CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
            CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
C
            !Get offset for the T1Y amplitudes
            KOFF1 = KT1B + MCAU*NT1AM(ISTAMY)
C
            !Construct T1-transformed integrals and store on file
            CALL CC3_BARINT(WORK(KOFF1),ISTAMY,XLAMDP,
     *                      XLAMDH,WORK(KEND0),LWRK0,
     *                      LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
            CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISTAMY,LU3SRTR,FN3SRTR,
     *                     LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                     IDUMMY,CDUMMY)
C
            CALL CC3_SINT(XLAMDH,WORK(KEND0),LWRK0,ISTAMY,
     *                    LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
            !Close the files for T1-transformed integrals
            CALL WCLOSE2(LUCKJDR,FNCKJDR,'KEEP')
            CALL WCLOSE2(LUDKBCR,FNDKBCR,'KEEP')
C
            !Close and delete temporary files
            CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
            CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
C
         END DO !MCAU
C
      ENDIF
c
      if (.false.) then
         kx3am  = kend0
         kend0 = kx3am + nt1amx*nt1amx*nt1amx
         call dzero(work(kx3am),nt1amx*nt1amx*nt1amx)
         lwrk0 = lwork - kend0
         if (lwrk0 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend0
            call quit('Insufficient space (T3) in CCSDT_ETA_FA_DEN')
         END IF
      endif
C
C-------------------------------------
C     Read canonical orbital energies:
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
 
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
 
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------------------------
C     Delete frozen orbitals in Fock diagonal.
C---------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
C-----------------------------------------------------
C     Sort the fock matrix
C-----------------------------------------------------
C
      DO ISYMC = 1,NSYM

         ISYMK = MULD2H(ISYMC,ISYM0)

         DO K = 1,NRHF(ISYMK)

            DO C = 1,NVIR(ISYMC)

               KOFF1 = IFCVIR(ISYMK,ISYMC) + 
     *                 NORB(ISYMK)*(C - 1) + K
               KOFF2 = IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C

               WORK(KFCKBA-1+KOFF2) = FOCK0(KOFF1)

            ENDDO
         ENDDO
      ENDDO

      IF (LOCDBG) THEN
       CALL AROUND('In CCSDT_ETA_FA_DEN: Triples Fock MO matrix (sort)')
       CALL CC_PRFCKMO(WORK(KFCKBA),ISYM0)
      ENDIF
C
C-----------------------------------------------------------
C     Prepare one-electron operators needed to compute the
C     amplitude response vectors:
C-----------------------------------------------------------
C
      IF ( (LISTR(1:3).EQ.'R1 ') .OR. (LISTR(1:3).EQ.'RC ') ) THEN
        ! read property integrals from file:
        CALL CCPRPAO(LABELY,.TRUE.,WORK(KFOCKY),IRREP,IMAT,IERR,
     *               WORK(KEND0),LWRK0)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISTAMY)) THEN
          WRITE(LUPRI,*) 'ISTAMY:',ISTAMY
          WRITE(LUPRI,*) 'IRREP :',IRREP
          WRITE(LUPRI,*) 'IERR  :',IERR
          CALL QUIT('CCSDT_ETA_FA_DEN: error reading operator '//LABELY)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(WORK(KFOCKY),N2BST(ISTAMY))
        END IF
 
        ! transform property integrals to Lambda-MO basis
        CALL CC_FCKMO(WORK(KFOCKY),XLAMDP,XLAMDH,WORK(KEND0),LWRK0,
     *                ISTAMY,1,1)
      END IF

C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
      KTROC0  = KEND0
      KTROC02 = KTROC0  + NTRAOC(ISYM0)
      KXIAJB  = KTROC02 + NTRAOC(ISYM0)
      KEND1   = KXIAJB  + NT2AM(ISYM0)
      IF (TOT_T3Y) THEN
         KTROC0Y = KEND1
         KEND1   = KTROC0Y + (NCAU+1)*NTRAOC(ISTAMY)
      ENDIF

      KTROC01 = KEND1
      KTROC21 = KTROC01 + NTRAOC(ISYM0)
      KTROC03 = KTROC21 + NTRAOC(ISYM0)
      KTROC23 = KTROC03 + NTRAOC(ISYM0)
      KEND1   = KTROC23 + NTRAOC(ISYM0)
      LWRK1   = LWORK   - KEND1

      KINTOC  = KEND1
      KEND2   = KINTOC + MAX(NTOTOC(ISYM0),NTOTOC(ISTAMY))
      LWRK2   = LWORK  - KEND2

      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CCSDT_ETA_FA_DEN')
      END IF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      LENGTH = IRAT*NT2AM(ISYM0)

      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))

      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYM0,1)
C
C------------------------
C     Occupied integrals.
C------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYM0) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF

      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of OCC-INT ',
     *    DDOT(NTOTOC(ISYM0),WORK(KINTOC),1,WORK(KINTOC),1)
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),XLAMDP,
     *               WORK(KEND2),LWRK2,ISYM0)

C
C--------------------------------------------------------------------
C     Read in T1Y-transformed occupied integrals used to construct W3
C--------------------------------------------------------------------
C
      IF (TOT_T3Y) THEN
C
         DO MCAU = 0,NCAU
C
            !Generate file name for T1Y-transformed integrals
            IF (LCAUCHY) THEN
               !(FNDKBCR is not needed here)
               CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
            ELSE
               FNCKJDR  = 'CCSDT_____ETAFD2'
            END IF
C
            !Open the file for T1Y-transformed integrals
            LUCKJDR  = -1
            CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
C
            !Get the offset for the integrals depending on MCAU
            KOFF1 = KTROC0Y + MCAU*NTRAOC(ISTAMY)
C
            IOFF = 1
            IF (NTOTOC(ISTAMY) .GT. 0) THEN
               CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,
     *                     NTOTOC(ISTAMY))
            ENDIF
C
            IF (LOCDBG)
     *         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (Y transformed) ',
     *           DDOT(NTOTOC(ISTAMY),WORK(KINTOC),1,WORK(KINTOC),1)
C
            CALL CC3_TROCC(WORK(KINTOC),WORK(KOFF1),XLAMDP,
     *                     WORK(KEND2),LWRK2,ISTAMY)
C
            !Close and delete the file for T1Y-transformed occupied integrals
            CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
C
         END DO !MCAU
C
      ENDIF
C
C-----------------------------------------------------------
C     Construct 2*C-E of the integrals.
C     Have integral for both (ij,k,a) and (a,k,j,i)
C-----------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYM0) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF
 
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC01),XLAMDH,
     *               WORK(KEND2),LWRK2,ISYM0)
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',
     *    DDOT(NTOTOC(ISYM0),WORK(KINTOC),1,WORK(KINTOC),1)
 
      CALL CCSDT_TCMEOCC(WORK(KTROC01),WORK(KTROC21),ISYM0)
 
      CALL CCFOP_SORT(WORK(KTROC01),WORK(KTROC03),ISYM0,1)
 
      CALL CCFOP_SORT(WORK(KTROC21),WORK(KTROC23),ISYM0,1)
 
      CALL CCFOP_SORT(WORK(KTROC0),WORK(KTROC02),ISYM0,1)
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of TROC0 ',
     *     DDOT(NTRAOC(ISYM0),WORK(KTROC0),1,WORK(KTROC0),1)
C
C---------------------------------------------
C     Open files for Tbar and W intermediates:
C---------------------------------------------
C
      LUTBAR = -1
      LUWMAT = -1
  
      CALL WOPEN2(LUTBAR,FNTBAR,64,0)
      CALL WOPEN2(LUWMAT,FNWMAT,64,0)

C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKB = MULD2H(ISYMD,ISYM0)
         ISCKBY = MULD2H(ISTAMY,ISYMD)
         IF (LOCDBG) 
     *      WRITE(LUPRI,*) 'In CCSDT_ETA_FA_DEN: ISYCKB :',ISYCKB

         DO D = 1,NVIR(ISYMD)
C
C           ------------------
C           Memory allocation.
C           ------------------
            KTRVI0  = KEND1
            KTRVI1  = KTRVI0 + NCKATR(ISYCKB)
            KTRVI2  = KTRVI1 + NCKATR(ISYCKB)
            KTRVI3  = KTRVI2 + NCKATR(ISYCKB)
            KTRVI4  = KTRVI3 + NCKATR(ISYCKB)
            KTRVI5  = KTRVI4 + NCKATR(ISYCKB)
            KTRVI6  = KTRVI5 + NCKATR(ISYCKB)
            KTRVI7  = KTRVI6 + NCKATR(ISYCKB)
            KEND3   = KTRVI7 + NCKATR(ISYCKB)
            LWRK3   = LWORK  - KEND3
           
            KTRVI14 = KEND3
            KTRVI15 = KTRVI14 + NCKATR(ISYCKB)
            KTRVI18 = KTRVI15 + NCKATR(ISYCKB)
            KTRVI19 = KTRVI18 + NCKATR(ISYCKB)
            KEND3   = KTRVI19 + NCKATR(ISYCKB)
            LWRK3   = LWORK  - KEND3

            IF (TOT_T3Y) THEN
               KTRVI0Y  = KEND3
               KEND3    = KTRVI0Y  + (NCAU+1)*NCKATR(ISCKBY)
               LWRK3    = LWORK    - KEND3
            ENDIF
           
            KINTVI = KEND3
            KEND4  = KINTVI + MAX(NCKA(ISYMD),NCKA(ISYCKB))
            LWRK4  = LWORK  - KEND4
           
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CCSDT_ETA_FA_DEN(99)')
            END IF
C
C           ------------------------------------
C           Integrals used in s3am.
C           ------------------------------------
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                     NCKATR(ISYCKB))
            ENDIF

C
C           -------------------------------------------------------------
C           Read T1Y-transformed virt ints (fixed D) used to construct W3
C           -------------------------------------------------------------
C
            IF (TOT_T3Y) THEN
C
               DO MCAU = 0, NCAU
C
                  !Generate file names for T1Y-transformed integrals
                  IF (LCAUCHY) THEN
                     !(FNCKJDR is not needed here)
                     CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
                  ELSE
                     FNDKBCR  = 'CCSDT_____ETAFD4'
                  END IF
C
                  !Open the files for T1Y-transformed integrals
                  LUDKBCR  = -1
                  CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
C
                  !Get the offset for a given MCAU
                  KOFF1 = KTRVI0Y + MCAU*NCKATR(ISCKBY)
C
                  !Read in the integrals
                  IOFF = ICKBD(ISCKBY,ISYMD) + NCKATR(ISCKBY)*(D-1) + 1
                  IF (NCKATR(ISCKBY) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KOFF1),IOFF,
     &                           NCKATR(ISCKBY))
                  ENDIF
C
                  !Close the file for T1B-transformed virtual integrals
                  !(and keep it: it will be needed in B-loop)
                  CALL WCLOSE2(LUDKBCR,FNDKBCR,'KEEP')
C
               END DO !MCAU
C
            ENDIF
C
C           -------------------------------------------
C           Read 2*C-E of integral used for t3-bar
C           -------------------------------------------
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KTRVI4),IOFF,
     &                     NCKATR(ISYCKB))
            ENDIF
C
C           -------------------------------------------------
C           Integrals used for t3-bar for cc3
C           -------------------------------------------------
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI14),IOFF,
     &                     NCKATR(ISYCKB))
            ENDIF
            CALL CCSDT_SRVIR3(WORK(KTRVI14),WORK(KEND4),
     *                        ISYMD,D,ISYM0)
            CALL CCSDT_SRTVIR(WORK(KTRVI14),WORK(KTRVI15),WORK(KEND4)
     *                        ,LWRK4,ISYMD,ISYM0)
C
C           ------------------------------------------------
C           Sort the integrals for s3am and for t3-bar
C           ------------------------------------------------
C
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISYM0)

            CALL CCSDT_SRTVIR(WORK(KTRVI4),WORK(KTRVI5),WORK(KEND4),
     *                        LWRK4,ISYMD,ISYM0)
C
C           -------------------------------------------
C           Read virtual integrals used in contraction.
C           -------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF

            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMDH,
     *                       ISYMD,D,ISYM0,WORK(KEND4),LWRK4)
C
C           ---------------------------------------------
C           Calculate virtual integrals used in q3am.
C           ---------------------------------------------
C
            CALL DCOPY(NCKATR(ISYCKB),WORK(KTRVI1),1,WORK(KTRVI3),1)

            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDT_ETA_FA_DEN (1)')
            END IF

            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND4),ISYMD,D,ISYM0)
C
C           ----------------------------------------------------
C           Read virtual integrals used in q3am/u3am for t3-bar.
C           ----------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF

            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI19),XLAMDP,
     *                       ISYMD,D,ISYM0,WORK(KEND4),LWRK4)

            IF (LWRK4 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     *                   'CCSDT_ETA_FA_DEN  (CC3 TRVI)')
            END IF

            CALL CCSDT_SRTVIR(WORK(KTRVI19),WORK(KTRVI18),WORK(KEND4)
     *                        ,LWRK4,ISYMD,ISYM0)

            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOP,FN3FOP,WORK(KTRVI6),IOFF,
     *                     NCKATR(ISYCKB))
            ENDIF

            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDT_ETA_FA_DEN (2)')
            END IF

            CALL DCOPY(NCKATR(ISYCKB),WORK(KTRVI6),1,WORK(KTRVI7),1)

            CALL CCSDT_SRVIR3(WORK(KTRVI6),WORK(KEND4),ISYMD,D,ISYM0)
C
C
            DO ISYMB = 1,NSYM

               ISYMBD  = MULD2H(ISYMB,ISYMD)

               ISALJB  = MULD2H(ISYMB,ISYM0)
               ISALJD = MULD2H(ISYMD,ISYM0)
               ISBLALJB  = MULD2H(ISYMB,ISYMBL)
               ISBLALJD = MULD2H(ISYMD,ISYMBL)
C
               ISCKIJ  = MULD2H(ISYMBD,ISYM0)
               ISBLCKIJ  = MULD2H(ISYMBD,ISYMBL)
C
               ISYCKD  = MULD2H(ISYM0,ISYMB)
               ISCKDY  = MULD2H(ISTAMY,ISYMB)
               ISWMAT  = MULD2H(ISCKIJ,ISTAMY)

               IF (LOCDBG) THEN
                  WRITE(LUPRI,*) 'In CCSDT_ETA_FA_DEN: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CCSDT_ETA_FA_DEN: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CCSDT_ETA_FA_DEN: ISALJB:',ISALJB
                  WRITE(LUPRI,*) 'In CCSDT_ETA_FA_DEN: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CCSDT_ETA_FA_DEN: ISCKIJ:',ISCKIJ
               ENDIF
C
C              Can use kend3 since we do not need the integrals anymore.
               KSMATBD   = KEND3
               KBLSMATBD  = KSMATBD   + NCKIJ(ISCKIJ)
               KSMATDB  = KBLSMATBD  + NCKIJ(ISBLCKIJ)
               KUMATBD   = KSMATDB  + NCKIJ(ISCKIJ)
               KBLUMATBD  = KUMATBD   + NCKIJ(ISCKIJ)
               KUMATDB  = KBLUMATBD  + NCKIJ(ISBLCKIJ)
               KDIAG   = KUMATDB  + NCKIJ(ISCKIJ)
               KDIAGW  = KDIAG   + NCKIJ(ISCKIJ)
               KINDSQW = KDIAGW  + NCKIJ(ISWMAT)
               KINDSQ  = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEXB  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEXD = KINDEXB  + (NCKI(ISALJB)  - 1)/IRAT + 1
               KTMAT   = KINDEXD + (NCKI(ISALJD) - 1)/IRAT + 1
               KTRVI8  = KTMAT   + MAX(NCKIJ(ISCKIJ),NCKIJ(ISBLCKIJ))
               KTRVI9  = KTRVI8  + NCKATR(ISYCKD)
               KTRVI10 = KTRVI9  + NCKATR(ISYCKD)
               KEND4   = KTRVI10 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4
C
               KBLINDEXB  = KEND4
               KBLINDEXD = KBLINDEXB  + (NCKI(ISBLALJB)  - 1)/IRAT + 1
               KEND4   = KBLINDEXD + (NCKI(ISBLALJD) - 1)/IRAT + 1
               LWRK4   = LWORK   - KEND4
C
               KBLINDSQ  = KEND4
               KEND4  = KBLINDSQ  + (6*NCKIJ(ISBLCKIJ) - 1)/IRAT + 1
               LWRK4   = LWORK   - KEND4
C
               KBLDIAG   = KEND4
               KEND4  = KBLDIAG   + NCKIJ(ISBLCKIJ)
               LWRK4   = LWORK   - KEND4
C
               KTRVI16 = KEND4
               KTRVI17 = KTRVI16 + NCKATR(ISYCKD)
               KTRVI20 = KTRVI17 + NCKATR(ISYCKD)
               KEND4   = KTRVI20 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4
               IF (TOT_T3Y) THEN
                  KTRVI8Y = KEND4
                  KEND4   = KTRVI8Y + (NCAU+1)*NCKATR(ISCKDY)
               ENDIF

               KBLSMATDB  = KEND4
               KBLUMATDB  = KBLSMATDB  + NCKIJ(ISBLCKIJ)
               KTRVI11 = KBLUMATDB  + NCKIJ(ISBLCKIJ)
               KTRVI12 = KTRVI11 + NCKATR(ISYCKD)
               KTRVI13 = KTRVI12 + NCKATR(ISYCKD)
               KEND4   = KTRVI13 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4

               KWMAT   = KEND4
               KEND4   = KWMAT   + NCKIJ(ISWMAT)
               LWRK4   = LWORK   - KEND4

               KINTVI  = KEND4
               KEND5   = KINTVI  + MAX(NCKA(ISYMB),NCKA(ISYCKD))
               LWRK5   = LWORK   - KEND5

               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CCSDT_ETA_FA_DEN')
               END IF
C
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KBLDIAG),WORK(KFOCKD),ISBLCKIJ)
C
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)

               IF (LOCDBG) THEN
                 WRITE(LUPRI,*) 'Norm of DIA  ',
     *              DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,WORK(KDIAG),1)
                 WRITE(LUPRI,*) 'Norm of DIA_W',
     *              DDOT(NCKIJ(ISWMAT),WORK(KDIAGW),1,WORK(KDIAGW),1)
               END IF
C
C              -----------------------
C              Construct index arrays.
C              -----------------------
C
               LENSQ  = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               LENSQBL  = NCKIJ(ISBLCKIJ)
               CALL CC3_INDSQ(WORK(KBLINDSQ),LENSQBL,ISBLCKIJ)
C
               LENSQW = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT) 
C
               CALL CC3_INDEX(WORK(KBLINDEXB),ISBLALJB)
               CALL CC3_INDEX(WORK(KBLINDEXD),ISBLALJD)
               CALL CC3_INDEX(WORK(KINDEXB),ISALJB)
               CALL CC3_INDEX(WORK(KINDEXD),ISALJD)

               DO B = 1,NVIR(ISYMB)
                  IF (LOCDBG) write(lupri,*) 'For B, D : ',B,D
C
C                 ----------------------------
C                 Read and transform integrals 
C                 ----------------------------
                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B - 1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI16),IOFF,
     *                           NCKATR(ISYCKD))
                  ENDIF
                  CALL CCSDT_SRVIR3(WORK(KTRVI16),WORK(KEND5),
     *                              ISYMB,B,ISYM0)
                  CALL CCSDT_SRTVIR(WORK(KTRVI16),WORK(KTRVI17),
     *                              WORK(KEND4),LWRK4,ISYMB,ISYM0)
C
C                 ------------------------------------
C                 Read and transform integrals used in 
C                 second S-bar and U-bar
C                 ------------------------------------
                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B-1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KTRVI11),
     *                           IOFF,NCKATR(ISYCKD))
                  ENDIF

                  CALL CCSDT_SRTVIR(WORK(KTRVI11),WORK(KTRVI12),
     *                              WORK(KEND5),LWRK5,ISYMB,
     *                              ISYM0)

                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B - 1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LU3FOP,FN3FOP,WORK(KTRVI13),IOFF,
     *                           NCKATR(ISYCKD))
                  ENDIF

                  IOFF = ICKAD(ISYCKD,ISYMB) 
     *                 + NCKA(ISYCKD)*(B - 1) + 1
                  IF (NCKA(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                           NCKA(ISYCKD))
                  ENDIF

                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI20),
     *                             XLAMDP,ISYMB,B,ISYM0,
     *                             WORK(KEND4),LWRK4)
C
C                 ----------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for B,D for T3-BAR:
C                 ----------------------------------------------------
                  CALL DZERO(WORK(KBLSMATBD),NCKIJ(ISBLCKIJ))
 
                  CALL CCFOP_SMAT(FREQL,ZL1AM,ISYMBL,ZL2TP,
     *                            ISYMBL,WORK(KTMAT),
     *                            WORK(KFCKBA),WORK(KXIAJB),ISYM0,
     *                            WORK(KTRVI14),WORK(KTRVI15),
     *                            WORK(KTRVI4),WORK(KTRVI5),
     *                            WORK(KTROC01),WORK(KTROC21),
     *                            ISYM0,WORK(KFOCKD),WORK(KBLDIAG),
     *                            WORK(KBLSMATBD),WORK(KEND4),LWRK4,
     *                            WORK(KBLINDEXB),
     *                            WORK(KBLINDSQ),LENSQBL,
     *                            ISYMB,B,ISYMD,D)

                  CALL DSCAL(NCKIJ(ISBLCKIJ),HALF,WORK(KBLSMATBD),1)

                  CALL T3_FORBIDDEN(WORK(KBLSMATBD),ISYMBL,ISYMB,B,
     *                              ISYMD,D)
C
C                 ----------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for D,B for T3-BAR:
C                 ----------------------------------------------------
                  CALL DZERO(WORK(KBLSMATDB),NCKIJ(ISBLCKIJ))

                  CALL CCFOP_SMAT(FREQL,ZL1AM,ISYMBL,ZL2TP,
     *                            ISYMBL,WORK(KTMAT),WORK(KFCKBA),
     *                            WORK(KXIAJB),ISYM0,
     *                            WORK(KTRVI16),WORK(KTRVI17),
     *                            WORK(KTRVI11),WORK(KTRVI12),
     *                            WORK(KTROC01),WORK(KTROC21),
     *                            ISYM0,WORK(KFOCKD),WORK(KBLDIAG),
     *                            WORK(KBLSMATDB),WORK(KEND4),LWRK4,
     *                            WORK(KBLINDEXD),WORK(KBLINDSQ),
     *                            LENSQBL,ISYMD,D,ISYMB,B)

                  CALL DSCAL(NCKIJ(ISBLCKIJ),HALF,WORK(KBLSMATDB),1)
C
C
                  CALL T3_FORBIDDEN(WORK(KBLSMATDB),ISYMBL,ISYMD,D,
     *                              ISYMB,B)
C
C                 ------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d for t3-bar.
C                 ------------------------------------------------
                  CALL DZERO(WORK(KBLUMATBD),NCKIJ(ISBLCKIJ))

                  CALL CCFOP_UMAT(FREQL,ZL1AM,ISYMBL,ZL2TP,
     *                            ISYMBL,
     *                            WORK(KXIAJB),ISYM0,WORK(KFCKBA),
     *                            WORK(KTRVI19),WORK(KTRVI7),
     *                            WORK(KTROC03),WORK(KTROC23),ISYM0,
     *                            WORK(KFOCKD),WORK(KBLDIAG),
     *                            WORK(KBLUMATBD),
     *                            WORK(KTMAT),WORK(KEND4),LWRK4,
     *                            WORK(KBLINDSQ),LENSQBL,
     *                            ISYMB,B,ISYMD,D)

                  CALL DSCAL(NCKIJ(ISBLCKIJ),HALF,WORK(KBLUMATBD),1)
C
C
                  CALL T3_FORBIDDEN(WORK(KBLUMATBD),ISYMBL,ISYMB,B,
     *                              ISYMD,D)
C
C                 ------------------------------------------------
C                 Calculate U(ci,jk) for fixed d,b for t3-bar.
C                 ------------------------------------------------
                  CALL DZERO(WORK(KBLUMATDB),NCKIJ(ISBLCKIJ))

                  CALL CCFOP_UMAT(FREQL,ZL1AM,ISYMBL,ZL2TP,
     *                            ISYMBL,WORK(KXIAJB),ISYM0,
     *                            WORK(KFCKBA),WORK(KTRVI20),
     *                            WORK(KTRVI13),WORK(KTROC03),
     *                            WORK(KTROC23),ISYM0,
     *                            WORK(KFOCKD),WORK(KBLDIAG),
     *                            WORK(KBLUMATDB),WORK(KTMAT),
     *                            WORK(KEND4),LWRK4,WORK(KBLINDSQ),
     *                            LENSQBL,ISYMD,D,ISYMB,B)

                  CALL DSCAL(NCKIJ(ISBLCKIJ),HALF,WORK(KBLUMATDB),1)

                  CALL T3_FORBIDDEN(WORK(KBLUMATDB),ISYMBL,ISYMD,D,
     *                              ISYMB,B)
C
C                 -------------------------------------------
C                 Sum up the S-bar and U-bar to get a real T3-bar
C                 -------------------------------------------
                  CALL CC3_T3BD(ISBLCKIJ,WORK(KBLSMATBD),
     *                                 WORK(KBLSMATDB),
     *                                 WORK(KBLUMATBD),WORK(KBLUMATDB),
     *                                 WORK(KTMAT),WORK(KBLINDSQ),
     *                                 LENSQBL)
C
C                 -------------------------------------
C                 Write T3-bar as T3^D(ai,bj,l) to disc
C                 -------------------------------------
                  DO ISYML = 1, NSYM
                   ISYMDL = MULD2H(ISYMD,ISYML)
                   ISAIBJ = MULD2H(ISYMBL,ISYMDL)
                   DO L = 1, NRHF(ISYML)
                    DO ISYMJ = 1, NSYM
                     ISYMBJ = MULD2H(ISYMJ,ISYMB)
                     ISYMAI = MULD2H(ISAIBJ,ISYMBJ)
                     ISYAIL = MULD2H(ISYMAI,ISYML)
                     ISAIB  = MULD2H(ISYMAI,ISYMB)
                     DO J = 1, NRHF(ISYMJ)
    
                      KOFFA = ISAIKJ(ISYAIL,ISYMJ)+ NCKI(ISYAIL)*(J-1)
     *                      + ISAIK(ISYMAI, ISYML)+NT1AM(ISYMAI)*(L-1)

                      IADR = ISWTL(ISAIBJ,ISYML) + NT2SQ(ISAIBJ)*(L-1)+
     *                       ISTLN(ISAIB,ISYMJ)  + NCKATR(ISAIB)*(J-1)+
     *                       ICKATR(ISYMAI,ISYMB)+ NT1AM(ISYMAI)*(B-1)+1

                      CALL PUTWA2(LUTBAR,FNTBAR,WORK(KTMAT+KOFFA),
     *                            IADR,NT1AM(ISYMAI))
                     END DO
                    END DO
                   END DO
                  END DO
C
C                 ---------------------------------------------
C                 Initialize WMAT:
C                 ---------------------------------------------
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))

C                 ---------------------------------------------
C                 Read and transform integrals used in second S
C                 ---------------------------------------------
                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B - 1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI8),IOFF,
     *                           NCKATR(ISYCKD))
                  ENDIF


                  CALL CCSDT_SRTVIR(WORK(KTRVI8),WORK(KTRVI9),
     *                              WORK(KEND4),LWRK4,ISYMB,ISYM0)
C
C                 --------------------------------------------------
C                 Read T1Y-transformed vir int (fixed B) used for W3
C                 --------------------------------------------------
C
                  IF (TOT_T3Y) THEN
C
                     DO MCAU = 0, NCAU
C
                        !Generate file names for T1Y-transformed integrals
                        IF (LCAUCHY) THEN
                           !(FNCKJDR is not needed here)
                           CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
                        ELSE
                           FNDKBCR  = 'CCSDT_____ETAFD4'
                        END IF
C
                        !Open the files for T1Y-transformed integrals
                        LUDKBCR  = -1
                        CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
C
                        !Get the offset for a given MCAU
                        KOFF1 = KTRVI8Y + MCAU*NCKATR(ISCKDY)
C
                        !Read in the integrals
                        IOFF = ICKBD(ISCKDY,ISYMB) + 
     &                         NCKATR(ISCKDY)*(B - 1) + 1
                        IF (NCKATR(ISCKDY) .GT. 0) THEN
                         CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KOFF1),IOFF,
     &                               NCKATR(ISCKDY))
                        ENDIF
C
                        !Close and keep the file for T1-transformed virt int
                        !...or delte it if you don't need it anymore

                        IF (      (ISYMD.EQ.NSYM) .AND. (ISYMB.EQ.NSYM)
     *                      .AND. (D.EQ.NVIR(ISYMD))
     *                      .AND. (B.EQ.NVIR(ISYMB))) THEN
C
                               CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
                        ELSE
                               CALL WCLOSE2(LUDKBCR,FNDKBCR,'KEEP')
                        END IF
C
                     END DO !MCAU
C
                  ENDIF
C
C                 -----------------------------------------
C                 Read virtual integrals used in second U
C                 -----------------------------------------
                  IOFF = ICKAD(ISYCKD,ISYMB) 
     *                 + NCKA(ISYCKD)*(B - 1) + 1
                  IF (NCKA(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                           NCKA(ISYCKD))
                  ENDIF

                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI10),
     *                             XLAMDH,ISYMB,B,ISYM0,
     *                             WORK(KEND5),LWRK5)
C
C                 --------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for B,D.
C                 --------------------------------------------------
                  CALL CC3_SMAT(0.0D0,T2TP,ISYM0,WORK(KTMAT),
     *                          WORK(KTRVI0),WORK(KTRVI2),
     *                          WORK(KTROC0),ISYM0,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMATBD),WORK(KEND4),LWRK4,
     *                          WORK(KINDEXB),WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)

                  CALL T3_FORBIDDEN(WORK(KSMATBD),ISYM0,ISYMB,B,ISYMD,D)
C
                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of SMAT :',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KSMATBD),1,WORK(KSMATBD),1)
C
C                 --------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for D,B.
C                 --------------------------------------------------
                  CALL CC3_SMAT(0.0D0,T2TP,ISYM0,WORK(KTMAT),
     *                          WORK(KTRVI8),WORK(KTRVI9),
     *                          WORK(KTROC0),ISYM0,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMATDB),WORK(KEND4),LWRK4,
     *                          WORK(KINDEXD),WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)

                  CALL T3_FORBIDDEN(WORK(KSMATDB),ISYM0,ISYMD,D,ISYMB,B)
C
                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of SMAT3:',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KSMATDB),1,WORK(KSMATDB),1)
C
C                 ---------------------------------
C                 Calculate U(ci,jk) for fixed b,d.
C                 ---------------------------------
                  CALL CC3_UMAT(0.0D0,T2TP,ISYM0,WORK(KTRVI1),
     *                          WORK(KTROC02),ISYM0,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMATBD),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)

                  CALL T3_FORBIDDEN(WORK(KUMATBD),ISYM0,ISYMB,B,ISYMD,D)
C
                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of UMAT :',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KUMATBD),1,WORK(KUMATBD),1)
C
C                 ---------------------------------
C                 Calculate U(ci,jk) for fixed d,b.
C                 ---------------------------------
                  CALL CC3_UMAT(0.0D0,T2TP,ISYM0,WORK(KTRVI10),
     *                          WORK(KTROC02),ISYM0,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMATDB),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)

                  CALL T3_FORBIDDEN(WORK(KUMATDB),ISYM0,ISYMD,D,ISYMB,B)

                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of UMAT3:',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KUMATDB),1,WORK(KUMATDB),1)
C
C                 -----------------------------------
C                 Sum up the S and U to get T30
C                 -----------------------------------
                  CALL CC3_T3BD(ISCKIJ,WORK(KSMATBD),WORK(KSMATDB),
     *                                 WORK(KUMATBD),WORK(KUMATDB),
     *                                 WORK(KTMAT),WORK(KINDSQ),LENSQ)
C
C                 -----------------------------------
C                 Use the T30 in TMAT to calculate W3
C                 -----------------------------------
                  IF ( (LISTR(1:3).EQ.'R1 ') .OR. 
     *                 (LISTR(1:3).EQ.'RC ') ) THEN
C
                     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of TMAT:',
     *                  DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,WORK(KTMAT),1)
C
C                    --------------------------------------------------------
C                    Calculate the  term <mu3|[Y,T3]|HF> virtual contribution 
C                    added in W^BD(KWMAT)
C                    --------------------------------------------------------
C
                     CALL WBD_V(WORK(KTMAT),ISCKIJ,WORK(KFOCKY),ISTAMY,
     *                          WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)

                     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_V):',
     *                  DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)

C
C                    ---------------------------------------------------------
C                    Calculate the  term <mu3|[Y,T3]|HF> occupied contribution 
C                    added in W^BD(KWMAT)
C                    ---------------------------------------------------------
C
                     CALL WBD_O(WORK(KTMAT),ISCKIJ,WORK(KFOCKY),ISTAMY,
     *                          WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)

                     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_O)',
     *                  DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)

C
C                    ---------------------------------------------------------
C                    Calculate the  term <mu3|[[Y,T2],T2]|HF>
C                    added in W^BD(KWMAT)
C                    ---------------------------------------------------------
C
                     CALL WBD_T2(.FALSE.,B,ISYMB,D,ISYMD,
     *                           T2TP,ISYM0,T2TP,ISYM0,
     *                           WORK(KFOCKY),ISTAMY,
     *                           WORK(KINDSQW),LENSQW,
     *                           WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)

                     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_T2)',
     *                  DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)
c
*                 call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *                        1,work(kx3am),4)
C
                  END IF
C
C                 --------------------------------------------------------------
C                 Calculate the terms <mu3|[H^0,T2Y]|HF> and <mu3|[H^Y,T2^0]|HF>
C                 added in W^BD(KWMAT)
C                 --------------------------------------------------------------
C
                  DO MCAU = 0, NCAU !for non-Cauchy calc this is dummy loop
C
                     IF (TOT_T3Y) THEN
C
                        !Get offset for KT2B vec for a given MCAU
                        KOFF2 = KT2B + MCAU*NT2SQ(ISTAMY)

                        !Calculate the term <mu3|[H^0,T2Y]|HF> 


                        CALL WBD_GROUND(WORK(KOFF2),ISTAMY,WORK(KTMAT),
     *                                  WORK(KTRVI0),WORK(KTRVI8),
     *                                  WORK(KTROC0),1,
     *                                  WORK(KWMAT),
     *                                  WORK(KEND4),LWRK4,
     *                                  WORK(KINDSQW),LENSQW,
     *                                  ISYMB,B,ISYMD,D)
C
                        !Get offset for T1Y-trans occ int for a given MCAU
                        KOFFOCC  = KTROC0Y + MCAU*NTRAOC(ISTAMY)
C
                        !Get offset for T1Y-trans virt int with fixed D for MCAU
                        KOFFINTD = KTRVI0Y + MCAU*NCKATR(ISCKBY)
                        !Get offset for T1Y-trans virt int with fixed B for MCAU
                        KOFFINTB = KTRVI8Y + MCAU*NCKATR(ISCKDY)
C
                        !Calculate the term <mu3|[H^Y,T2^0]|HF>
                        CALL WBD_GROUND(T2TP,ISYM0,WORK(KTMAT),
     *                                  WORK(KOFFINTD),WORK(KOFFINTB),
     *                                  WORK(KOFFOCC),ISTAMY,
     *                                  WORK(KWMAT),
     *                                  WORK(KEND4),LWRK4,
     *                                  WORK(KINDSQW),LENSQW,
     *                                  ISYMB,B,ISYMD,D)
                     ENDIF

                     CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQY, 
     *                            ISWMAT,WORK(KWMAT), 
     *                            WORK(KDIAGW),WORK(KFOCKD)) 

                     CALL T3_FORBIDDEN(WORK(KWMAT),ISTAMY,
     *                                 ISYMB,B,ISYMD,D) 
C
                  END DO !MCAU
C
                  !AFTER ALL THE SIGN SHOULD NOT BE TURNED !!!!!

*                 IF (LCAUCHY) THEN
*                    !trun the sign C_T  <-  -C_T
*                    CALL DSCAL(NCKIJ(ISTAMY),-1.0D0,WORK(KWMAT),1)
*                 END IF !LCAUCHY

*                 call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *                        1,work(kx3am),4)

                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_DIA)',
     *               DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)
C
C                 -------------------------------------
C                 Write WMAT as WMAT^D(ai,bj,l) to disc
C                 -------------------------------------
                  DO ISYML = 1, NSYM
                   ISYMDL = MULD2H(ISYMD,ISYML)
                   ISAIBJ = MULD2H(ISTAMY,ISYMDL)
                   DO L = 1, NRHF(ISYML)
                    DO ISYMJ = 1, NSYM
                     ISYMBJ = MULD2H(ISYMJ,ISYMB)
                     ISYMAI = MULD2H(ISAIBJ,ISYMBJ)
                     ISYAIL = MULD2H(ISYMAI,ISYML)
                     DO J = 1, NRHF(ISYMJ)
    
                      KOFFA = ISAIKJ(ISYAIL,ISYMJ)+ NCKI(ISYAIL)*(J-1)
     *                      + ISAIK(ISYMAI, ISYML)+NT1AM(ISYMAI)*(L-1)

                      NBJ  = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
                      IADR = ISWTL(ISAIBJ,ISYML)+NT2SQ(ISAIBJ)*(L-1)+
     *                    IT2SQ(ISYMAI,ISYMBJ)+NT1AM(ISYMAI)*(NBJ-1)+1

                      CALL PUTWA2(LUWMAT,FNWMAT,WORK(KWMAT+KOFFA),
     *                            IADR,NT1AM(ISYMAI))
                     END DO
                    END DO
                   END DO
                  END DO

               ENDDO   ! B
            ENDDO      ! ISYMB

            DO ISYML = 1, NSYM

             ISYMDL = MULD2H(ISYMD,ISYML)
             ISWMAT = MULD2H(ISYMDL,ISTAMY)
             ISTBAR = MULD2H(ISYMBL,ISYMDL)

             KWMAT  = KEND1
             KT2Y   = KWMAT + NT2SQ(ISWMAT)
             KEND2  = KT2Y  + NT1AM(ISWMAT)
             LWRK2  = LWORK - KEND2

             IF ( LWRK2 .LT. 0 ) THEN
               CALL QUIT('Out of memory in CCSDT_ETA_FA_DEN (x)')
             ENDIF

             DO L = 1, NRHF(ISYML)

C              --------------------------------------------
C              Read WMAT from file and generate WMAT-tilde:
C              --------------------------------------------
               IADR = ISWTL(ISWMAT,ISYML) + NT2SQ(ISWMAT)*(L-1) + 1
               CALL GETWA2(LUWMAT,FNWMAT,WORK(KWMAT),IADR,NT2SQ(ISWMAT))

               CALL CC_T2MOD(WORK(KWMAT),ISWMAT,HALF)

C              --------------------------------------------
C              Read response doubles amplitudes T2^y,DL:
C              --------------------------------------------
               NDL  = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
               IADR = IT2SQ(ISWMAT,ISYMDL) + NT1AM(ISWMAT)*(NDL-1) + 1

               CALL GETWA2(LUT2Y,FNT2Y,WORK(KT2Y),IADR,NT1AM(ISWMAT))


C              -----------------------------------
C              Loop over outermost occupied index:
C              -----------------------------------
               DO ISYMN = 1, NSYM
                ISYEMF = MULD2H(ISTBAR,ISYMN)

                KTBAR  = KEND2
                KEND3  = KTBAR + NCKATR(ISYEMF)
                LWRK3  = LWORK - KEND2

                IF ( LWRK3 .LT. 0 ) THEN
                  CALL QUIT('Out of memory in CCSDT_ETA_FA_DEN (x2)')
                ENDIF

                DO N = 1, NRHF(ISYMN)

C                ----------------------
C                Read T3-BAR from file:
C                ----------------------
                 IADR = ISWTL(ISTBAR,ISYML) + NT2SQ(ISTBAR)*(L-1) + 
     *                  ISTLN(ISYEMF,ISYMN) + NCKATR(ISYEMF)*(N-1)+1
                 CALL GETWA2(LUTBAR,FNTBAR,WORK(KTBAR),IADR,
     *                       NCKATR(ISYEMF))
                 
C                -------------------------------------------------------
C                D(fb) <- D(fb)+ sum_em tbar^DLN(em,f) Wtilde^DL(em,bN):
C                -------------------------------------------------------
                 DO ISYMBN = 1, NSYM
                   ISYMEM = MULD2H(ISWMAT,ISYMBN)
                   ISYMB  = MULD2H(ISYMBN,ISYMN)
                   ISYMF  = MULD2H(ISYEMF,ISYMEM)

                   KOFFT  = KTBAR + ICKATR(ISYMEM,ISYMF)
                   KBN    = IT1AM(ISYMB,ISYMN)+NVIR(ISYMB)*(N-1)+1
                   KOFFW  = KWMAT + IT2SQ(ISYMEM,ISYMBN) +
     *                              NT1AM(ISYMEM)*(KBN-1)
                   KOFFD  = 1  + IMATAB(ISYMF,ISYMB)
                 
                   NTOTEM = MAX(NT1AM(ISYMEM),1)
                   NNVIRF = MAX(NVIR(ISYMF),1)
                    
                   CALL DGEMM('T','N',NVIR(ISYMF),NVIR(ISYMB),
     *                        NT1AM(ISYMEM),
     *                        ONE,WORK(KOFFT),NTOTEM,WORK(KOFFW),NTOTEM,
     *                        ONE,DAB(KOFFD),NNVIRF)
                 
                 END DO

C                -------------------------------------------------------
C                D(iN) <- D(iN)- sum_em tbar^DLN(em,f) Wtilde^DL(em,fi):
C                -------------------------------------------------------
                 DO ISYMFI = 1, NSYM
                   ISYMEM = MULD2H(ISWMAT,ISYMFI)
                   ISYMF  = MULD2H(ISYEMF,ISYMEM)
                   ISYMI  = MULD2H(ISYMFI,ISYMF)
                   
                   KOFFT  = KTBAR + ICKATR(ISYMEM,ISYMF)
                   KOFFW  = KWMAT + IT2SQ(ISYMEM,ISYMFI) +
     *                              NT1AM(ISYMEM)*IT1AM(ISYMF,ISYMI)
                   KOFFD  = 1  + IMATIJ(ISYMI,ISYMN) + 
     *                              NRHF(ISYMI)*(N-1)

                   NNEMF  = MAX(NT1AM(ISYMEM)*NVIR(ISYMF),1)

                   CALL DGEMV('T',NT1AM(ISYMEM)*NVIR(ISYMF),NRHF(ISYMI),
     *                        -ONE,WORK(KOFFW),NNEMF,WORK(KOFFT),1,
     *                         ONE,DIJ(KOFFD),1)
                 
                 END DO

C                -------------------------------------------------------
C                M(imfN)   <-- M(imfN) + sum_em t^DL(ei) tbar^DLN(em,f):
C                -------------------------------------------------------
                 IF (DO_MMAT) THEN
                  DO ISYMF = 1, NSYM
                   ISYMFN = MULD2H(ISYMF,ISYMN)
                   ISYMEM = MULD2H(ISYEMF,ISYMF)
                   DO ISYMM = 1, NSYM
                     ISYME  = MULD2H(ISYMEM,ISYMM)
                     ISYMI  = MULD2H(ISYMDL,ISYME) 
                     ISYMIM = MULD2H(ISYMM, ISYMI)
                     ISYMEI = MULD2H(ISYME, ISYMI)
                     ISYEIL = MULD2H(ISYMEI,ISYML)

                     DO F = 1, NVIR(ISYMF)

                       KOFFA  = 1 + IT2SP(ISYEIL,ISYMD) +
     *                    NCKI(ISYEIL) *(D-1) + ISAIK(ISYMEI,ISYML) +
     *                    NT1AM(ISYMEI)*(L-1) + IT1AM(ISYME,ISYMI)

                       KOFFT  = KTBAR + ICKATR(ISYMEM,ISYMF) +
     *                    NT1AM(ISYMEM)*(F-1) + IT1AM(ISYME,ISYMM)

                       KFN = IT1AM(ISYMF,ISYMN) + NVIR(ISYMF)*(N-1) + F
                       KOFFM  = 1 + ISAIKL(ISYMFN,ISYMIM) +
     *                    NMATIJ(ISYMIM)*(KFN-1) + IMATIJ(ISYMI,ISYMM)

                       NVIRE  = MAX(NVIR(ISYME),1)
                       NRHFI  = MAX(NRHF(ISYMI),1)

                       CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMM),
     *                         NVIR(ISYME),
     *                         ONE,T2TP(KOFFA),NVIRE,WORK(KOFFT),NVIRE,
     *                         ONE,XMMAT(KOFFM),NRHFI)

                     END DO
                   END DO
                  END DO
                 END IF

C                -------------------------------------------------------
C                y^M(imfN) <- M(imfN)+ sum_em y^t^DL(ei) tbar^DLN(em,f):
C                -------------------------------------------------------
                 IF (DO_YMMAT) THEN
                  DO ISYMF = 1, NSYM
                   ISYMFN = MULD2H(ISYMF,ISYMN)
                   ISYMEM = MULD2H(ISYEMF,ISYMF)
                   DO ISYMM = 1, NSYM
                     ISYME  = MULD2H(ISYMEM,ISYMM)
                     ISYMI  = MULD2H(MULD2H(ISYMDL,ISYME),ISTAMY)
                     ISYMIM = MULD2H(ISYMM, ISYMI)
                     ISYMEI = MULD2H(ISYME, ISYMI)
                     ISYEIL = MULD2H(ISYMEI,ISYML)

                     DO F = 1, NVIR(ISYMF)

                       KOFFA  = KT2Y  + IT1AM(ISYME,ISYMI)

                       KOFFT  = KTBAR + ICKATR(ISYMEM,ISYMF) +
     *                    NT1AM(ISYMEM)*(F-1) + IT1AM(ISYME,ISYMM)

                       KFN = IT1AM(ISYMF,ISYMN) + NVIR(ISYMF)*(N-1) + F
                       KOFFM  = 1 + ISAIKL(ISYMFN,ISYMIM) +
     *                    NMATIJ(ISYMIM)*(KFN-1) + IMATIJ(ISYMI,ISYMM)

                       NVIRE  = MAX(NVIR(ISYME),1)
                       NRHFI  = MAX(NRHF(ISYMI),1)

                       CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMM),
     *                        NVIR(ISYME),
     *                        ONE,WORK(KOFFA),NVIRE,WORK(KOFFT),NVIRE,
     *                        ONE,YMMAT(KOFFM),NRHFI)

                     END DO
                   END DO
                  END DO
                 END IF

                IF (DO_DIA) THEN
C                ------------------------------------------------------
C                D(Lb) <-- D(Lb) - sum_em tbar^DN(em) Wtilde^DL(em,bN):
C                ------------------------------------------------------
                 ISYMDN = MULD2H(ISYMD,ISYMN)
                 ISYMEM = MULD2H(ISYMBL,ISYMDN)
                 ISYMBN = MULD2H(ISWMAT,ISYMEM)
                 ISYMB  = MULD2H(ISYMBN,ISYMN)
                 ISYEMN = MULD2H(ISYMEM,ISYMN)

                 KOFFT  = 1 + IT2SP(ISYEMN,ISYMD) + 
     *                    NCKI(ISYEMN) * (D-1) + ISAIK(ISYMEM,ISYMN)+
     *                    NT1AM(ISYMEM)* (N-1)

                 KBN    = IT1AM(ISYMB,ISYMN)+NVIR(ISYMB)*(N-1)+1
                 KOFFW  = KWMAT  + IT2SQ(ISYMEM,ISYMBN) +
     *                             NT1AM(ISYMEM)*(KBN-1)

                 KOFFD  = 1 + IT1AM(ISYMB,ISYML)+NVIR(ISYMB)*(L-1)

                 NTOTEM = MAX(NT1AM(ISYMEM),1)

                 CALL DGEMV('T',NT1AM(ISYMEM),NVIR(ISYMB),
     *                      -ONE,WORK(KOFFW),NTOTEM,ZL2TP(KOFFT),1,
     *                       ONE,DIA(KOFFD),1)
                   

C                -----------------------------------------------------
C                D(mb) <-- D(mb) - sum_e tbar^DL(eN) Wtilde^DL(em,bN):
C                -----------------------------------------------------
                 DO ISYMBN = 1, NSYM
                   ISYMEM = MULD2H(ISWMAT,ISYMBN)
                   ISYMB  = MULD2H(ISYMBN,ISYMN)
                   ISYMEN = MULD2H(ISYMBL,ISYMDL)
                   ISYME  = MULD2H(ISYMEN,ISYMN)
                   ISYMM  = MULD2H(ISYMEM,ISYME)
                   ISYENL = MULD2H(ISYMEN,ISYML)

                   NVIRE = MAX(NVIR(ISYME),1)
                   NVIRB = MAX(NVIR(ISYMB),1)

                   DO B = 1, NVIR(ISYMB)

                     KOFFT  = 1 + IT2SP(ISYENL,ISYMD) + 
     *                        NCKI(ISYENL) * (D-1)+ISAIK(ISYMEN,ISYML)+
     *                        NT1AM(ISYMEN)* (L-1)+IT1AM(ISYME,ISYMN) +
     *                        NVIR(ISYME)  * (N-1)

                     KBN    = IT1AM(ISYMB,ISYMN)+NVIR(ISYMB)*(N-1)+B
                     KOFFW  = KWMAT + IT2SQ(ISYMEM,ISYMBN) +
     *                        NT1AM(ISYMEM)*(KBN-1)+IT1AM(ISYME,ISYMM)

                     KOFFD  = IT1AM(ISYMB,ISYMM) + B

                     CALL DGEMV('T',NVIR(ISYME),NRHF(ISYMM),
     *                          -ONE,WORK(KOFFW),NVIRE,ZL2TP(KOFFT),1,
     *                           ONE,DIA(KOFFD),NVIRB)
                   END DO
                   
                 END DO
                END IF ! DO_DIA

                END DO ! N
               END DO ! ISYMN

              IF (DO_DIA) THEN
C              --------------------------------------------------------
C              D(nb) <-- D(nb) + 2 sum_em tbar^DL(em) Wtilde^DL(em,bN):
C              --------------------------------------------------------
               ISYMEM = MULD2H(ISYMBL,ISYMDL)
               ISYMBN = MULD2H(ISWMAT,ISYMEM)
               ISYEML = MULD2H(ISYMEM,ISYML)

               KOFFT  = 1 + IT2SP(ISYEML,ISYMD) + 
     *                  NCKI(ISYEML) * (D-1) + ISAIK(ISYMEM,ISYML) +
     *                  NT1AM(ISYMEM)* (L-1)

               KOFFW  = KWMAT + IT2SQ(ISYMEM,ISYMBN)

               NTOTEM = MAX(NT1AM(ISYMEM),1)
   
               CALL DGEMV('T',NT1AM(ISYMEM),NT1AM(ISYMBN),
     *                    TWO,WORK(KOFFW),NTOTEM,ZL2TP(KOFFT),1,
     *                    ONE,DIA,1)
              END IF ! DO_DIA


             END DO ! L
            END DO ! ISYML
 
         ENDDO       ! D
      ENDDO          ! ISYMD
 
C
C-------------
C     End
C-------------
C
      CALL WCLOSE2(LUTBAR,FNTBAR,'DELETE')
      CALL WCLOSE2(LUWMAT,FNWMAT,'DELETE')
C
C---------------------------------------------------------------------
C     It might have happened that 'CC3_CAUINT_V*' or 'CCSDT_____ETAFD4'
C     files have not been deleted up there. Do it now!
C---------------------------------------------------------------------
C     
      IF (LCAUCHY) THEN
         CALL DELETE_FILES('CC3_CAUINT_V*')
      ELSE
         CALL DELETE_FILES('CCSDT_____ETAFD4')
      END IF

C
*     call print_pt3(work(kx3am),1,4)
C
      CALL QEXIT('CTETAFAD')
C
      RETURN
      END
*=====================================================================*
      SUBROUTINE CC_T2MOD(T2AMP,ISYAMP,FAC)
*---------------------------------------------------------------------*
*  
*     Purpose: construct for a squared vector the following 
*              combination in place:
*
*              T2-new(ai,bj) = T2-old(ai,bj) + FAC * T2-old(bj,ai)
*
*    Written by Christof Haettig, Mai 2002, Aarhus
*             
*=====================================================================*
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL T2AMP(*), XAIBJ, XBJAI, FAC
#else
      DOUBLE PRECISION T2AMP(*), XAIBJ, XBJAI, FAC
#endif

      INTEGER ISYMAI, ISYMBJ, NAI, NBJ, NAIBJ, NBJAI, NAIAI, ISYAMP

      CALL QENTER('CC_T2MOD')

      DO ISYMAI = 1, NSYM
        ISYMBJ = MULD2H(ISYMAI,ISYAMP)

        IF (ISYMAI.GT.ISYMBJ) THEN

          DO NAI = 1, NT1AM(ISYMAI)
            DO NBJ = 1, NT1AM(ISYMBJ)
   
              NAIBJ = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)+NAI
              NBJAI = IT2SQ(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)*(NAI-1)+NBJ

              XAIBJ = T2AMP(NAIBJ)
              XBJAI = T2AMP(NBJAI)

              T2AMP(NAIBJ) = XAIBJ + FAC * XBJAI
              T2AMP(NBJAI) = XBJAI + FAC * XAIBJ

            END DO
          END DO

        ELSE IF (ISYMAI.EQ.ISYMBJ) THEN

          DO NAI = 1, NT1AM(ISYMAI)
            DO NBJ = 1, NAI-1
              
              NAIBJ = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1)+NAI
              NBJAI = IT2SQ(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)*(NAI-1)+NBJ

              XAIBJ = T2AMP(NAIBJ)
              XBJAI = T2AMP(NBJAI)

              T2AMP(NAIBJ) = XAIBJ + FAC * XBJAI
              T2AMP(NBJAI) = XBJAI + FAC * XAIBJ

            END DO

            NAIAI = IT2SQ(ISYMAI,ISYMAI) + NT1AM(ISYMAI)*(NAI-1)+NAI
            T2AMP(NAIAI) = (1.0D0 + FAC) * T2AMP(NAIAI) 

          END DO

        END IF

      END DO

      CALL QEXIT('CC_T2MOD')
      RETURN
      END 
*=====================================================================* 
*=====================================================================*
      SUBROUTINE CCSDT_ETA_TM2(DIA,ISYDIA,XMMAT,ISMMAT,T2TP,
     &                         LUT2,FNT2,ISYMT2,IOPTT2,WORK,LWORK)
*---------------------------------------------------------------------*
*  
*     Purpose: compute density matrix contribution
*
*              D(ia) <-- D(ia) + sum_fnm T(am,fn) M(imfn):
*
*              Note: D(ia) is stored as D(ai) using NT1AM, IT1AM
*
*              IOPTT2 = 0 : take T2 from T2TP
*                     = 1 : read T2 from LUT2
*
*    Written by Christof Haettig, Mai 2002, Aarhus
*             
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL DIA(*), XMMAT(*), T2TP(*), WORK(*)
      REAL ONE
#else
      DOUBLE PRECISION DIA(*), XMMAT(*), T2TP(*), WORK(*)
      DOUBLE PRECISION ONE
#endif
      PARAMETER(ONE=1.0D0)

      CHARACTER*(*) FNT2
      INTEGER LUT2, ISYDIA, ISMMAT, ISYMT2, LWORK, IOPTT2

      INTEGER ISYMFN, ISYMAM, ISYMIM, KTAM, KEND1, LWRK1, ISYMN,
     *        ISYMF, ISYAMN, NFN, IADR, ISYMM, ISYMI, ISYMA, KOFFM,
     *        KOFFT, KOFFW, NVIRA, NRHFI, KOFFD

      CALL QENTER('CCSDT_ETA_TM')

      IF (MULD2H(ISYMT2,ISMMAT).NE.ISYDIA)
     *  CALL QUIT('Symmetry mismatch in CCSDT_ETA_TM2')

      DO ISYMFN = 1, NSYM
        ISYMAM = MULD2H(ISYMFN,ISYMT2)
        ISYMIM = MULD2H(ISYMFN,ISMMAT)

        KTAM  = 1
        KEND1 = KTAM  + NT1AM(ISYMAM)
        LWRK1 = LWORK - KEND1

        IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CCSDT_ETA_TM2')
        ENDIF

        DO ISYMN = 1, NSYM
         ISYMF  = MULD2H(ISYMFN,ISYMN)
         ISYAMN = MULD2H(ISYMAM,ISYMN)

         DO N = 1, NRHF(ISYMN)
         DO F = 1, NVIR(ISYMF) 
          NFN = IT1AM(ISYMF,ISYMN) + NVIR(ISYMF) * (N-1) + F
        
          IF      (IOPTT2.EQ.0) THEN
            IADR  = 1 + IT2SP(ISYAMN,ISYMF) + NCKI(ISYAMN) *(F-1) + 
     *                  ISAIK(ISYMAM,ISYMN) + NT1AM(ISYMAM)*(N-1) 
            CALL DCOPY(NT1AM(ISYMAM),T2TP(IADR),1,WORK(KTAM),1)
          ELSE IF (IOPTT2.EQ.1) THEN
            IADR = IT2SQ(ISYMAM,ISYMFN) + NT1AM(ISYMAM)*(NFN-1) + 1
            CALL GETWA2(LUT2,FNT2,WORK(KTAM),IADR,NT1AM(ISYMAM))
          ELSE
            CALL QUIT('Illegal case IOPTT2 in CCSDT_ETA_TM')
          END IF
        
          DO ISYMM = 1, NSYM
            ISYMI = MULD2H(ISYMIM,ISYMM)
            ISYMA = MULD2H(ISYMAM,ISYMM)
        
            KOFFM = 1 + ISAIKL(ISYMFN,ISYMIM) +
     *              NMATIJ(ISYMIM)*(NFN-1) + IMATIJ(ISYMI,ISYMM)
        
            KOFFT = KTAM + IT1AM(ISYMA,ISYMM)
        
            KOFFD = 1 + IT1AM(ISYMA,ISYMI)
        
            NVIRA = MAX(NVIR(ISYMA),1)
            NRHFI = MAX(NRHF(ISYMI),1)
        
            CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMM),
     *                 ONE,WORK(KOFFT),NVIRA,XMMAT(KOFFM),NRHFI,
     *                 ONE,DIA(KOFFD),NVIRA)
        
          END DO

         END DO ! N 
         END DO ! F
        
        END DO ! ISYMN 
      END DO ! ISYMFN

      CALL QEXIT('CCSDT_ETA_TM2')
      RETURN
      END 
*=====================================================================* 
*=====================================================================*
      SUBROUTINE CC_FOCK_RESORT(FIJ,LOO,FIA,LOV,FAI,LVO,FAB,LVV,
     &                          FOCK,ISYFCK)
*---------------------------------------------------------------------*
*  
*     Purpose: grep occupied/occupied, occupied/virtual etc. blocks
*              out of the complete MO Fock matrix
*
*              if LOO.eq.true  --> return occupied/occupied block FIJ
*              if LOV.eq.true  --> return occupied/virtual  block FIA
*              if LVO.eq.true  --> return virtual/occupied  block FAI
*              if LVV.eq.true  --> return virtual/virtual   block FAB
*
*
*    Written by Christof Haettig, Mai 2002, Aarhus
*             
*=====================================================================*
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>

#if defined (SYS_CRAY)
      REAL FIA(*), FIJ(*), FAI(*), FAB(*), FOCK(*)
#else
      DOUBLE PRECISION FIA(*), FIJ(*), FAI(*), FAB(*), FOCK(*)
#endif

      LOGICAL LOV, LOO, LVV, LVO
      INTEGER ISYFCK, ISYMK, ISYMC, ISYMI, ISYMJ, ISYMA, ISYMB,
     &        KOFF1, KOFF2

      CALL QENTER('CC_FOCK_RESORT')


      IF (LOV) THEN
        DO ISYMC = 1,NSYM
          ISYMK = MULD2H(ISYMC,ISYFCK)
          DO K = 1,NRHF(ISYMK)
            DO C = 1,NVIR(ISYMC)
              KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
              KOFF2 = IT1AM(ISYMC,ISYMK)  + NVIR(ISYMC)*(K - 1) + C
              FIA(KOFF2) = FOCK(KOFF1)
            ENDDO
          ENDDO
        ENDDO
      END IF

      IF (LVO) THEN
        DO ISYMC = 1,NSYM
          ISYMK = MULD2H(ISYMC,ISYFCK)
          DO K = 1,NRHF(ISYMK)
            KOFF1 = 1+IFCRHF(ISYMK,ISYMC)+NORB(ISYMC)*(K-1)+NRHF(ISYMK)
            KOFF2 = 1+IT1AM(ISYMC,ISYMK) +NVIR(ISYMC)*(K-1)
            CALL DCOPY(NVIR(ISYMC),FOCK(KOFF1),1,FAI(KOFF2),1)
          ENDDO
        ENDDO
      END IF


      IF (LVV) THEN
       DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYFCK)
         DO B = 1,NVIR(ISYMB)
           KOFF1 = 1+IFCVIR(ISYMA,ISYMB)+NORB(ISYMA)*(B-1)+NRHF(ISYMA) 
           KOFF2 = 1+IMATAB(ISYMA,ISYMB)+NVIR(ISYMA)*(B-1) 
           CALL DCOPY(NVIR(ISYMA),FOCK(KOFF1),1,FAB(KOFF2),1)
         END DO
       END DO
      END IF

      IF (LOO) THEN
        DO ISYMJ = 1,NSYM
          ISYMI = MULD2H(ISYMJ,ISYFCK)
          DO J = 1,NRHF(ISYMJ)
            KOFF1 = 1 + IFCRHF(ISYMI,ISYMJ) + NORB(ISYMI)*(J - 1)
            KOFF2 = 1 + IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1)
            CALL DCOPY(NRHF(ISYMI),FOCK(KOFF1),1,FIJ(KOFF2),1) 
          END DO
        END DO
      END IF

      CALL QEXIT('CC_FOCK_RESORT')
      RETURN
      END 
*=====================================================================* 
*=====================================================================*
      SUBROUTINE CC_XIETA_DENPREP(IXETRAN,NXETRAN,MXVEC,
     &                            IXDOTS,LISTDPX,
     &                            IEDOTS,LISTDPE,LISTL,
     &                            IDXL_XIDEN,NXIDEN,MXXIDEN,
     &                            IDXL_EADEN,IDXR_EADEN,NEADEN,MXEADEN,
     &                            IEASTEND,NEALEFT,MXEALEFT)
*---------------------------------------------------------------------*
*  
*     Purpose: set up loops for the calculation of effective Xi/Eta
*              densities
*
*    Written by Christof Haettig, Mai 2002, Aarhus
*             
*=====================================================================*
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccroper.h>
#include <cclists.h>

* input:
      CHARACTER*3 LISTL, LISTDPE, LISTDPX
      INTEGER NXETRAN, IORDER, MXVEC
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      INTEGER IXDOTS(MXVEC,NXETRAN), IEDOTS(MXVEC,NXETRAN)

* output for Xi densities:
      INTEGER NXIDEN, MXXIDEN
      INTEGER IDXL_XIDEN(MXXIDEN)

* output for Eta densities:
      INTEGER NEADEN, MXEADEN, NEALEFT, MXEALEFT
      INTEGER IDXL_EADEN(MXEADEN)
      INTEGER IDXR_EADEN(MXEADEN)
      INTEGER IEASTEND(2,MXEALEFT)

* local variables:
      CHARACTER*8 LABELH
      LOGICAL LTWOEL, LRELAX, SKIPXI, SKIPETA, CHANGES
      INTEGER IVEC, ITRAN, IOPER, IDLSTL, IRELAX, ISYHOP, ISYCTR,
     &        ISYETA, IDLSTR, ISYMR, ISYML, IDX, IEADEN, IXIDEN

* external functions:
      INTEGER ILSTSYM

      CALL QENTER('CXIETADP')


      NEADEN = 0
      NXIDEN = 0

      DO ITRAN = 1, NXETRAN
        IOPER  = IXETRAN(1,ITRAN)  ! operator index
        IDLSTL = IXETRAN(2,ITRAN)  ! Left vector index
        IRELAX = IXETRAN(5,ITRAN)  ! flag for relax. contrib.

        ISYHOP = ISYOPR(IOPER)     ! symmetry of O operator
        LABELH = LBLOPR(IOPER)     ! operator label
        LTWOEL = LPDBSOP(IOPER)    ! two-electron contribution to O?

        SKIPXI  = ( IXETRAN(3,ITRAN) .EQ. -1 )
        SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )

        LRELAX  = LTWOEL .OR. (IRELAX.GE.1)
        IF (LRELAX) CALL QUIT(
     &    'Relaxed perturbations not yet implemented CC_XIETA_DENPREP')
 

C       ---------------------------------------------------------------
C       set up list of effecitive Eta{O} x R or L x A{O} x R densities:
C       ---------------------------------------------------------------
        IF (.NOT.SKIPETA) THEN
          ISYCTR = ILSTSYM(LISTL,IDLSTL) ! sym. of Left (ZETA) vector
          ISYETA = MULD2H(ISYHOP,ISYCTR) ! sym. of ETA result vect.

          IVEC = 1
          DO WHILE(IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            IDLSTR = IEDOTS(IVEC,ITRAN)
            ISYMR  = ILSTSYM(LISTDPE,IDLSTR)

            IF (ISYMR.EQ.ISYETA) THEN

              ! check if respective density already on our to-do list
              IEADEN = 0
              DO IDX = 1, NEADEN
                IF (IDLSTL.EQ.IDXL_EADEN(IDX) .AND.
     &              IDLSTR.EQ.IDXR_EADEN(IDX)       ) THEN
                  IEADEN = IDX
                END IF 
              END DO
              
              ! not found, put on to-do list for densities
              IF (IEADEN.EQ.0) THEN
                NEADEN = NEADEN + 1
                IF (NEADEN.GT.MXEADEN) CALL QUIT('NEADEN out of range')
                IDXL_EADEN(NEADEN) = IDLSTL
                IDXR_EADEN(NEADEN) = IDLSTR
              END IF

            END IF 

            IVEC = IVEC + 1
          END DO
        END IF

C       --------------------------------------------
C       set up list of effecitive L Xi{O} densities:
C       --------------------------------------------
        IF (.NOT.SKIPXI) THEN
          IVEC = 1
          DO WHILE(IXDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            IDLSTL = IXDOTS(IVEC,ITRAN)
            ISYML  = ILSTSYM(LISTDPX,IDLSTL)

            IF (ISYML.EQ.ISYHOP) THEN

              ! check if respective density already on our to-do list
              IXIDEN = 0
              DO IDX = 1, NXIDEN
                IF (IDLSTL.EQ.IDXL_XIDEN(IDX)) THEN
                  IXIDEN = IDX
                END IF 
              END DO
              
              ! not found, put on to-do list for densities
              IF (IXIDEN.EQ.0) THEN
                NXIDEN = NXIDEN + 1
                IF (NXIDEN.GT.MXXIDEN) CALL QUIT('NXIDEN out of range')
                IDXL_XIDEN(NXIDEN) = IDLSTL
              END IF

            END IF 

            IVEC = IVEC + 1
          END DO
        END IF

      END DO

C     -------------------------------------------------------------
C     sort list of effecitive Eta{O} x R or L x A{O} x R densities:
C     -------------------------------------------------------------
      CHANGES = .TRUE.
      DO WHILE(CHANGES)
        CHANGES = .FALSE.
        DO IEADEN = 2, NEADEN
          IF (IDXL_EADEN(IEADEN-1).GT.IDXL_EADEN(IEADEN)) THEN
            CHANGES = .TRUE.
            IDLSTL = IDXL_EADEN(IEADEN)
            IDLSTR = IDXR_EADEN(IEADEN)
            IDXL_EADEN(IEADEN) = IDXL_EADEN(IEADEN-1)
            IDXR_EADEN(IEADEN) = IDXR_EADEN(IEADEN-1)
            IDXL_EADEN(IEADEN-1) = IDLSTL
            IDXR_EADEN(IEADEN-1) = IDLSTR
          END IF
        END DO
      END DO

C     -------------------------------------------------------
C     set start and end indeces for loop over left vectors in
C     calculation of effective Eta densities:
C     -------------------------------------------------------
      IF (NEADEN.GE.1) THEN
        NEALEFT       = 1
        IEASTEND(1,1) = 1
        IEASTEND(2,1) = 1
        DO IEADEN = 2, NEADEN
          IF (IDXL_EADEN(IEADEN-1).NE.IDXL_EADEN(IEADEN)) THEN
            NEALEFT = NEALEFT + 1
            IEASTEND(1,NEALEFT) = IEADEN
          END IF
          IEASTEND(2,NEALEFT) = IEADEN
        END DO
      ELSE
        NEALEFT       = 0
        IEASTEND(1,1) = 0
        IEASTEND(2,1) = 0
      END IF

      CALL QEXIT('CXIETADP')
      RETURN
      END 
*=====================================================================* 
*=====================================================================*
#if defined (SYS_CRAY)
      REAL
#else
      DOUBLE PRECISION 
#endif
     & FUNCTION CC_XIETA_DENCON(IDENS,LABELH,ISYHOP,XLAMDP,XLAMDH,
     &                          LUDEFF,FNDEFF,M2BST,WORK,LWORK)
*---------------------------------------------------------------------*
*  
*     Purpose: contract effective density with one-electron integrals
*
*     Written by Christof Haettig, Mai 2002, Friedrichstal
*             
*=====================================================================*
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER LABELH*(*), FNDEFF*(*)
      INTEGER IDENS, ISYHOP, LUDEFF, M2BST, LWORK

#if defined (SYS_CRAY)
      REAL 
#else
      DOUBLE PRECISION
#endif
     &  WORK(LWORK), XLAMDP(*), XLAMDH(*), DDOT

      INTEGER KDAB,KDIJ,KDIA,KEND1,KFOCK,KFOCKIJ,KFOCKIA,KFOCKAB,
     &        LWRK1, IADRIJ, IADRAB, IADRIA, IRREP, IMAT, IERR

      CALL QENTER('CXIETDC')

      KDAB    = 1
      KDIJ    = KDAB    + NMATAB(ISYHOP)
      KDIA    = KDIJ    + NMATIJ(ISYHOP)
      KEND1   = KDIA    + NT1AM(ISYHOP)

      KFOCK   = KEND1
      KEND1   = KFOCK   + N2BST(ISYHOP)

      KFOCKIJ = KEND1
      KFOCKIA = KFOCKIJ + NMATIJ(ISYHOP)
      KFOCKAB = KFOCKIA +  NT1AM(ISYHOP)
      KEND1   = KFOCKAB + NMATAB(ISYHOP)

      LWRK1   = LWORK   - KEND1
      IF (LWRK1.LT.0) CALL QUIT('Out of memory in CC_XIETA_DENCON')

      IADRIJ = 1 + M2BST*(IDENS-1)
      CALL GETWA2(LUDEFF,FNDEFF,WORK(KDIJ),IADRIJ,NMATIJ(ISYHOP))

      IADRAB = IADRIJ + NMATIJ(ISYHOP)
      CALL GETWA2(LUDEFF,FNDEFF,WORK(KDAB),IADRAB,NMATAB(ISYHOP))

      IADRIA = IADRAB + NMATAB(ISYHOP)
      CALL GETWA2(LUDEFF,FNDEFF,WORK(KDIA),IADRIA,NT1AM(ISYHOP))

      CALL CCPRPAO(LABELH,.TRUE.,WORK(KFOCK),IRREP,IMAT,IERR,
     *             WORK(KEND1),LWRK1)
      IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYHOP)) THEN
        WRITE(LUPRI,*) 'ISYHOP:',ISYHOP
        WRITE(LUPRI,*) 'IRREP :',IRREP
        WRITE(LUPRI,*) 'IERR  :',IERR
        WRITE(LUPRI,*) 'LABEL :',LABELH
        CALL QUIT('CC_XIETA_DENCON: error reading operator ')
      ELSE IF (IERR.LT.0) THEN
        CALL DZERO(WORK(KFOCK),N2BST(ISYHOP))
      END IF

      ! transform property integrals to Lambda-MO basis
      CALL CC_FCKMO(WORK(KFOCK),XLAMDP,XLAMDH,WORK(KEND1),LWRK1,
     *              ISYHOP,1,1)

      CALL CC_FOCK_RESORT(WORK(KFOCKIJ),.TRUE.,WORK(KFOCKIA),.TRUE.,
     &                    DUMMY,.FALSE.,       WORK(KFOCKAB),.TRUE.,
     &                    WORK(KFOCK),ISYHOP)

      CC_XIETA_DENCON = 
     &    DDOT( NT1AM(ISYHOP),WORK(KDIA),1,WORK(KFOCKIA),1)  +
     &    DDOT(NMATIJ(ISYHOP),WORK(KDIJ),1,WORK(KFOCKIJ),1)  +
     &    DDOT(NMATAB(ISYHOP),WORK(KDAB),1,WORK(KFOCKAB),1)

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CC_XIETA_DENCON>', ISYHOP,LABELH
        WRITE(LUPRI,*) 'NORM^2(DIJ):',
     &    DDOT(NMATIJ(ISYHOP),WORK(KDIJ),1,WORK(KDIJ),1)
        WRITE(LUPRI,*) 'NORM^2(DAB):',
     &    DDOT(NMATAB(ISYHOP),WORK(KDAB),1,WORK(KDAB),1)
        WRITE(LUPRI,*) 'NORM^2(DIA):',
     &    DDOT(NT1AM(ISYHOP),WORK(KDIA),1,WORK(KDIA),1)
        WRITE(LUPRI,*) 'NORM^2(FOCKIJ):',
     &    DDOT(NMATIJ(ISYHOP),WORK(KFOCKIJ),1,WORK(KFOCKIJ),1)
        WRITE(LUPRI,*) 'NORM^2(KFOCKAB):',
     &    DDOT(NMATAB(ISYHOP),WORK(KFOCKAB),1,WORK(KFOCKAB),1)
        WRITE(LUPRI,*) 'NORM^2(FOCKIA):',
     &    DDOT(NT1AM(ISYHOP),WORK(KFOCKIA),1,WORK(KFOCKIA),1)
        WRITE(LUPRI,*) 'DIA x FOCKIA:',
     &    DDOT(NT1AM(ISYHOP),WORK(KDIA),1,WORK(KFOCKIA),1)
        WRITE(LUPRI,*) 'DIJ x FOCKIJ:',
     &    DDOT(NMATIJ(ISYHOP),WORK(KDIJ),1,WORK(KFOCKIJ),1)
        WRITE(LUPRI,*) 'DAB x FOCKAB:',
     &    DDOT(NMATAB(ISYHOP),WORK(KDAB),1,WORK(KFOCKAB),1)
      END IF

      CALL QEXIT('CXIETDC')
      RETURN
      END 
*=====================================================================* 
