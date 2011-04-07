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
C  /* Deck cc3_xi_cont */
      SUBROUTINE CC3_XI_CONT(LISTL,DOTPROD,
     &                        IDXL_XIDEN,NXIDEN,
     &                        IXDOTS,IXETRAN,MXVEC,NXETRAN,
     &                        WORK,LWORK)
*---------------------------------------------------------------------*
*     Purpose: compute density corresponding to a contraction of      *
*              the triples parts of a Xi vector and (first-order)     *
*              lagrangian multipliers - needed for assymetric form    *
*              of polarizability                                      *
*                                                                     *
*              t3bar^Y xi3^X                                          *
*                            = <L3Y|[X,T3]|HF>                        *
*                            + <L3Y|[[X,T2],T2]|HF>                   *
*                            + <L2Y|[X,T3]|HF>                        *
*                                                                     *
* Written by Poul Jorgensen and Filip Pawlowski, Fall 2002, Aarhus    *
*---------------------------------------------------------------------*
      IMPLICIT NONE  

#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <cclists.h>
#include <ccroper.h>

      CHARACTER*8 FNTOC, FN3VI2
      CHARACTER*6 FN3VI, FN3FOPX, FN3FOP2, FNCKJD, FNDELD
      CHARACTER*5 FNDKBC3, FN3FOP
      CHARACTER*7 FN3FOP2X
      CHARACTER*4 FNDKBC

      PARAMETER (FNTOC   = 'CCSDT_OC' , FN3VI   = 'CC3_VI')
      PARAMETER (FN3VI2  = 'CC3_VI12' , FNDKBC3 = 'DKBC3' )
      PARAMETER (FN3FOPX  = 'PTFOPX'   , FN3FOP2X = 'PTFOP2X'  )
      PARAMETER (FN3FOP  = 'PTFOP'   , FN3FOP2 = 'PTFOP2'  )
      PARAMETER (FNCKJD  = 'CKJDEL'  , FNDKBC  = 'DKBC')
      PARAMETER (FNDELD  = 'CKDELD')

      CHARACTER*12 FNDEFF
      PARAMETER( FNDEFF = 'TMP-XIETADEN' )

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*3 LISTL
      INTEGER LWORK

      INTEGER NXETRAN, MXVEC, NXIDEN
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      INTEGER IXDOTS(MXVEC,NXETRAN)
      INTEGER IDXL_XIDEN(NXIDEN)

#if defined (SYS_CRAY)
      REAL WORK(LWORK), DUMMY, TRIPCON
      REAL DOTPROD(MXVEC,NXETRAN),ONE
#else
      DOUBLE PRECISION WORK(LWORK), DUMMY, TRIPCON
      DOUBLE PRECISION DOTPROD(MXVEC,NXETRAN),ONE
#endif

      CHARACTER MODEL*10, LABELH*8
      LOGICAL LTWOEL, LRELAX, SKIPXI
      INTEGER LUTOC,LU3VI,LU3VI2,LUDKBC3,LU3FOPX,LU3FOP2X,LU3FOP  
      INTEGER LU3FOP2,LUCKJD,LUDKBC,LUDELD  
      INTEGER LUDEFF
      INTEGER M2BST,ISYM,IDENS,ISYML
      INTEGER KLAMP0,KLAMH0,KEND1,LWRK1,KDIJ,KDAB,KDIA,KT1AMP
      INTEGER IADRIJ,IADRAB,IADRIA,IOPT,ITRAN,IOPER,IRELAX,ISYHOP
      INTEGER IVEC,IDLSTL,IDX
C
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
C
      PARAMETER (ONE = 1.0D0)
C
      INTEGER ILSTSYM
#if defined (SYS_CRAY)
      REAL CC_XIETA_DENCON
#else
      DOUBLE PRECISION CC_XIETA_DENCON
#endif

      
*---------------------------------------------------------------------*
*     some initializations:
*---------------------------------------------------------------------*
      CALL QENTER('CC3XICNT')

      LUDEFF = -1
      CALL WOPEN2(LUDEFF,FNDEFF,64,0)

      M2BST = 0
      DO ISYM = 1, NSYM
        M2BST = MAX(M2BST,N2BST(ISYM))
      END DO 
C
*---------------------------------------------------------------------*
*     precalculate density matrices:
*---------------------------------------------------------------------*
   
      DO IDENS = 1, NXIDEN
     
        IDLSTL = IDXL_XIDEN(IDENS)
        ISYML  = ILSTSYM(LISTL,IDLSTL)
C
        KDIJ   = 1
        KDAB   = KDIJ   + NMATIJ(ISYML)
        KDIA   = KDAB   + NMATAB(ISYML)
        KEND1  = KDIA   + NT1AM(ISYML)
        LWRK1  = LWORK - KEND1
C
        IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient space in CC3_XI_CONT (2)')
        ENDIF
C
        CALL DZERO(WORK(KDIJ),NMATIJ(ISYML))
        CALL DZERO(WORK(KDAB),NMATAB(ISYML))
        CALL DZERO(WORK(KDIA),NT1AM(ISYML))
C
             LUTOC    = -1
             LU3VI    = -1
             LU3VI2   = -1
             LUDKBC3  = -1
             LU3FOPX  = -1
             LU3FOP2X = -1
             LU3FOP   = -1
             LU3FOP2  = -1
             LUCKJD   = -1
             LUDKBC   = -1
             LUDELD   = -1
C
             CALL WOPEN2(LUTOC,FNTOC,64,0)
             CALL WOPEN2(LU3VI,FN3VI,64,0)
             CALL WOPEN2(LU3VI2,FN3VI2,64,0)
             CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
             CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
             CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)
             CALL WOPEN2(LU3FOP,FN3FOP,64,0)
             CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
             CALL WOPEN2(LUCKJD,FNCKJD,64,0)
             CALL WOPEN2(LUDKBC,FNDKBC,64,0)
             CALL WOPEN2(LUDELD,FNDELD,64,0)
C
        
        CALL CC3_XI_DEN(WORK(KDIJ),WORK(KDAB),WORK(KDIA),LISTL,IDLSTL,
     *                   LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                   LUDKBC3,FNDKBC3,
     *                   LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                   LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                   LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                   WORK(KEND1),LWRK1)

        IF (IPRINT .GT. 55) THEN
          WRITE(LUPRI,*)'DAB density after CC3_XI_DEN '
          CALL PRINT_MATAB(WORK(KDAB),ISYML)
          WRITE(LUPRI,*)'DIJ density after CC3_XI_DEN '
          CALL PRINT_MATIJ(WORK(KDIJ),ISYML)
          WRITE(LUPRI,*)'DIA density after CC3_XI_DEN '
          CALL PRINT_MATAI(WORK(KDIA),ISYML)
        END IF


C
             CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
             CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
             CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
             CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
             CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
             CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')
             CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
             CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')
             CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
             CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
             CALL WCLOSE2(LUDELD,FNDELD,'KEEP')



        IADRIJ = 1 + M2BST*(IDENS-1)
C
C       WRITE(LUPRI,*)'Density D(ij) written on file in CC3_XI_CONT'
C       CALL PRINT_MATIJ(WORK(KDIJ),ISYML)
C
        CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDIJ),IADRIJ,NMATIJ(ISYML))

        IADRAB = IADRIJ + NMATIJ(ISYML)
C
C       WRITE(LUPRI,*)'Density D(ab) written on file in CC3_XI_CONT'
C       CALL PRINT_MATAB(WORK(KDAB),ISYML)
C
        CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDAB),IADRAB,NMATAB(ISYML))


        IADRIA = IADRAB + NMATAB(ISYML)
        CALL DSCAL(NT1AM(ISYML),-ONE,WORK(KDIA),1)
C
C       WRITE(LUPRI,*)'Density D(ia) after scaling written on file'
C    *                //' in CC3_XI_CONT'
C       CALL PRINT_MATAI(WORK(KDIA),ISYML)
C
        CALL PUTWA2(LUDEFF,FNDEFF,WORK(KDIA),IADRIA,NT1AM(ISYML))

      END DO

*---------------------------------------------------------------------*
*     calculate triples contributions to L x Xi{O} contractions:
*---------------------------------------------------------------------*
C
*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KT1AMP = 1
      KLAMP0 = KT1AMP + NT1AM(ISYM0)
      KLAMH0 = KLAMP0 + NLAMDT
      KEND1 = KLAMH0 + NLAMDT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC3_XI_CONT (1)')
      ENDIF
C
*---------------------------------------------------------------------*
*     initialize 0.th-order Lambda:
*---------------------------------------------------------------------*
      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP),DUMMY)

      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP),
     &            WORK(KEND1),LWRK1)
C
      DO ITRAN = 1, NXETRAN
        IOPER  = IXETRAN(1,ITRAN)  ! operator index
        IRELAX = IXETRAN(5,ITRAN)  ! flag for relax. contrib.

        ISYHOP = ISYOPR(IOPER)     ! symmetry of O operator
        LABELH = LBLOPR(IOPER)     ! operator label
        LTWOEL = LPDBSOP(IOPER)    ! two-electron contribution to O?

        SKIPXI = ( IXETRAN(3,ITRAN) .EQ. -1 )

        LRELAX = LTWOEL .OR. (IRELAX.GE.1)
        IF (LRELAX) CALL QUIT(
     &    'Relaxed contributions not yet implemented in CC3_XI_CONT')

        IF (.NOT.SKIPXI) THEN
          IVEC = 1
          DO WHILE(IXDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            IDLSTL = IXDOTS(IVEC,ITRAN)
            ISYML  = ILSTSYM(LISTL,IDLSTL)

            IF (ISYML.EQ.ISYHOP) THEN
                
              ! find index of density
              IDENS = 0
              DO IDX = 1, NXIDEN
                IF (IDLSTL.EQ.IDXL_XIDEN(IDX)) THEN
                  IDENS = IDX
                END IF
              END DO

              IF (IDENS.EQ.0) 
     &          CALL QUIT('Density not found in CC3_XI_CONT')

              TRIPCON = CC_XIETA_DENCON(IDENS,LABELH,ISYHOP,
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     LUDEFF,FNDEFF,M2BST,WORK(KEND1),LWRK1)

              IF (LOCDBG) THEN
                WRITE(LUPRI,*) 'CC3_XI_CONT>'
                WRITE(LUPRI,*) 'IVEC,ITRAN,LABELH:',IVEC,ITRAN,LABELH
                WRITE(LUPRI,*) 'tbar3 x xi3 :',TRIPCON 
              END IF

              DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TRIPCON

            END IF 

            IVEC = IVEC + 1
          END DO

        END IF
      END DO ! ITRAN

      CALL WCLOSE2(LUDEFF,FNDEFF,'DELETE')

      CALL QEXIT('CC3XICNT')
      RETURN
      END
*---------------------------------------------------------------------*


C  /* Deck cc3_xi_den */
      SUBROUTINE CC3_XI_DEN(DIJ,DAB,DIA,LISTBY,IDLSTBY,
     *                   LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                   LUDKBC3,FNDKBC3,
     *                   LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                   LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                   LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                   WORK,LWORK)

*---------------------------------------------------------------------*
*
*    Purpose: compute triples component of eta vector 
*             projected into the singles and doubles space
*             
*    W3BAR^Y = ( <L2|[Y,tau3]|HF> + <L3|[Y^,tau3]|HF>
*
*            +   <L2|[H^Y,tau3]|HF> 
*
*            +    <L1^Y + L2^Y|[H^,tau3]|HF> )/ (w+epsiln(tau3))
*
*
*    Written by Poul Jorgensen and Filip Pawlowski, Fall 2002, Aarhus
*            
*=====================================================================*
C
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
#include <cc3t3d.h>
#include <ccl1rsp.h>
#include <ccr1rsp.h>
#include <cclrmrsp.h>
#include <ccexci.h>
#include <ccn2rsp.h>

C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
      CHARACTER LISTL*3, LISTBY*3, LISTY*3, MODEL*10, LABELB*8
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X 
      CHARACTER*(*) FNDKBC,FNDELD,FN3VI2,FN3FOP,FN3FOP2,FNCKJD
C
      CHARACTER*10 FNT3
      CHARACTER*13 FNWBMAT
      PARAMETER(FNT3 = 'CC3_T3_TMP', FNWBMAT = 'CC3_WBMAT_TMP')
C

      LOGICAL   LOCDBG,LORXB
      PARAMETER (LOCDBG = .FALSE.)
C
      INTEGER   IDLSTL, LWORK
      INTEGER   ISYRES
      INTEGER   LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER   KLAMP0,KLAMH0,KFOCKD,KFCKBA,KT2TP,KL1AM,KL2TP
      INTEGER   KEND0,LWRK0,KT1AMP0,KEND1,LWRK1
      INTEGER   ISYCTR, ISYMC, ISYMK, KOFF1, KOFF2
      INTEGER   KXIAJB, KT3BOG1, KT3BOL1, KT3BOG2, KT3BOL2
      INTEGER   KEND2, LWRK2
      INTEGER   IOPT, LENGTH, ISYCKB, ISYMD, ISYMBD
      INTEGER   KT3BVDL1, KT3BVDL2, KT3BVDL3, KEND3, LWRK3
      INTEGER   KT3BVDG1, KT3BVDG2, KT3BVDG3, KINTVI, KTRVI6
      INTEGER   KEND4, LWRK4, IOFF, ISYMB, ISYALJ, ISYALJ2, ISCKIJ
      INTEGER   ISYCKD, LENSQ
      INTEGER   KSMAT2, KUMAT2, KDIAG, KINDSQ, KINDEX, KINDEX2, KTMAT
      INTEGER   KT3BVBG1, KT3BVBG2, KT3BVBG3, KSMAT4, KUMAT4, KT3BVBL1
      INTEGER   KT3BVBL2, KT3BVBL3, KEND5, LWRK5
      INTEGER   KWMAT, ISWMAT, ISYMIM
      INTEGER   KINDSQW, LENSQW
      INTEGER   ISCKB2, KDIAGW
      INTEGER   ISINT1, ISINT2, ISCKB1
      INTEGER   ISYML,ISYMDL,ISAIBJ,ISYMJ,ISYMBJ,ISYMA,ISYAIL,ISAIB
      INTEGER   ISYMAI,IADR,NBJ
      INTEGER   LUDKBC,LUDELD,LU3VI2,LU3FOP,LU3FOP2,LUCKJD
      INTEGER   LUT3,LUWBMAT,LUFOCK
      INTEGER   KS3MAT,KU3MAT,KS3MAT3,KU3MAT3,KT3VBG1,KT3VBG2,KT3VBG3
      INTEGER   KW3BXVDGX1,KW3BXVDGX2,KW3BXVDLX1,KW3BXVDLX2
      INTEGER   ISYALJBY,ISYALJDY,ISINT2Y,ISYFCKY,ISINT1Y,ISYCKBY
      INTEGER   KINDEXBY,KINDEXDY,KT3MAT,KT3OG1,KT3OG2,KLAMPY
      INTEGER   KLAMHY,KFCKYCK,KW3BXOG1,KW3BXOL1,KW3BXOGX1
      INTEGER   KW3BXOLX1,KT1Y,KT3VDG1,KT3VDG2,KT3VDG3,KW3BXVDG1
      INTEGER   KW3BXVDG2,KW3BXVDL1,KW3BXVDL2,KL1BY,KL2TPBY,KFOCK0
      INTEGER   KFOCKY
      INTEGER   ISYMBY,IDLSTBY,IDLSTY
      INTEGER   ISYMY
      INTEGER   IRREP,IERR
      INTEGER   IR1TAMP
      INTEGER   ISYM
      INTEGER ISYALJLE, ISYALJ2LE, KINDEXLE, KINDEX2LE
C
C
COMMENT COMENT
COMMENT COMENT
      INTEGER kx3am
COMMENT COMENT
C
COMMENT COMENT
C

#if defined (SYS_CRAY)
      REAL      FREQ, FREQY, HALF, DDOT
      REAL      WORK(LWORK)
      REAL      XNORMVAL
      REAL      DAB(*),DIJ(*),DIA(*)
#else
      DOUBLE PRECISION      FREQ, FREQY, HALF, DDOT
      DOUBLE PRECISION      WORK(LWORK)
      DOUBLE PRECISION      XNORMVAL
      DOUBLE PRECISION      DAB(*),DIJ(*),DIA(*)
#endif
C
      PARAMETER(HALF = 0.5D0)
C
      INTEGER ILSTSYM
C
      CALL QENTER('CC3_XI_DEN')
C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
*     LISTL = 'L0 '
*     IDLSTL = 0 
*     ISYCTR = ISYM0

      IF (LISTBY(1:3).EQ.'L1 ') THEN

         ! get symmetry, frequency and integral label from common blocks
         ! defined in ccl1rsp.h
         ISYMBY  = ISYLRZ(IDLSTBY)
         FREQ  = FRQLRZ(IDLSTBY)
         LABELB = LRZLBL(IDLSTBY)
         LORXB  = LORXLRZ(IDLSTBY)

         IF (LORXB) CALL QUIT('NO ORBITAL RELAX. IN CC3_XI_DEN')

        ! LISTY has to be R1 when LISTBY is L1 !!!
        LISTY  = 'R1 '
        IDLSTY = IR1TAMP(LABELB,LORXB,FREQ,ISYMBY)
        ! get symmetry and frequency from common blocks
        ! defined in ccl1rsp.h
        ISYMY  = ISYLRT(IDLSTY)
        FREQY  = FRQLRT(IDLSTY)
    
         !get the associated L0 list
         LISTL = 'L0 '
         IDLSTL = 0 
         ISYCTR = ISYM0
         IF (FREQ .NE. FREQY) THEN
            WRITE(LUPRI,*) 'FREQ = ', FREQ
            WRITE(LUPRI,*) 'FREQY = ', FREQY
            CALL QUIT('Frequency mismatch in CC3_XI_DEN(0)')
         END IF
C
      ELSE IF (LISTBY(1:3).EQ.'M1 ') THEN

        ISYMBY = ISYLRM(IDLSTBY)
Cold
C        ISYMBY = ILSTSYM(LISTBY,IDLSTBY) 
Cold
        FREQ  = FRQLRM(IDLSTBY)
        LABELB = '- none -'
        ! LISTY has to be RE when LISTBY is M1 !!!
        LISTY  = 'RE '
        IDLSTY = ILRM(IDLSTBY)
        ! get symmetry and frequency from common blocks
        ISYMY   = ILSTSYM(LISTY,IDLSTY)
        FREQY  = FREQ

         !get the associated L0 list
         LISTL = 'L0 '
         IDLSTL = 0 
         ISYCTR = ISYM0
C
         IF (FREQ .NE. FREQY) THEN
            WRITE(LUPRI,*) 'FREQ = ', FREQ
            WRITE(LUPRI,*) 'FREQY = ', FREQY
            CALL QUIT('Frequency mismatch in CC3_XI_DEN(1)')
         END IF

      ELSE IF (LISTBY(1:3).EQ.'LE ') THEN
        ISYMBY = ILSTSYM(LISTBY,IDLSTBY)
        FREQ  = -EIGVAL(IDLSTBY)
        LABELB = '- none -'
        ! LISTY is empty when LISTBY is LE !!!
        LISTY  = '---'
        IDLSTY = -99
        ! get symmetry and frequency from common blocks
        ISYMY   = 1
        FREQY  = FREQ
    
         !get the associated L0 list WHICH IS NOT REALLY USED IN THIS CASE
         LISTL = 'L0 '
         IDLSTL = 0
         ISYCTR = ISYM0
***************************
*       IMPORTANT !!!
***************************
        ! The term  <L2|[H^Y,tau3]|HF> which corresponds to 
        ! LISTY is absent for LISTBY = 'LE '.
        ! To avoid allocation problems WE HAVE INITIALIZED all
        ! the arrays refering to an artificial ISYMY = 1
        ! and JUMPED places where files are read in using ISYMY
        ! or LISTY.
***************************
         IF (FREQ .NE. FREQY) THEN
            WRITE(LUPRI,*) 'FREQ = ', FREQ
            WRITE(LUPRI,*) 'FREQY = ', FREQY
            CALL QUIT('Frequency mismatch in CC3_XI_DEN(2)')
         END IF
      ELSE IF (LISTBY(1:3).EQ.'N2 ') THEN
         ISYMBY = ILSTSYM(LISTBY,IDLSTBY)
         FREQ   = FRQIN2(IDLSTBY) + FRQFN2(IDLSTBY)
         LABELB = '- none -'

        ! LISTY has to be RE when LISTBY is N2 !!!
        LISTY  = 'RE '
        IDLSTY = IFN2(IDLSTBY)
        ! get symmetry and frequency from common blocks
        ISYMY   = ILSTSYM(LISTY,IDLSTY)
        FREQY  = FRQFN2(IDLSTBY)

       !get the associated zero-order left list 
        LISTL = 'LE '
        IDLSTL = IIN2(IDLSTBY)
        ISYCTR = ILSTSYM(LISTL,IDLSTL)

      ELSE
        CALL QUIT('Unknown list in CC3_XI_DEN')
      END IF
C
c     IF (FREQ .NE. FREQY) THEN
c        WRITE(LUPRI,*) 'FREQ = ', FREQ
c        WRITE(LUPRI,*) 'FREQY = ', FREQY
c        CALL QUIT('Frequency mismatch in CC3_XI_DEN')
c     END IF
C
C-------------------------------------------------------
C     initial allocations, orbital energy, fock matrix and T2 and L2 :
C-------------------------------------------------------
C
C     Symmetry of integrals in contraction:
C
      ISINT1 = 1
      ISINT2 = 1
C
      KLAMP0 = 1
      KLAMH0  = KLAMP0  + NLAMDT
      KFOCKD  = KLAMH0  + NLAMDT
      KFCKBA  = KFOCKD  + NORBTS
      KT2TP   = KFCKBA  + NT1AMX 
      KL1AM   = KT2TP   + NT2SQ(ISYM0)
      KL2TP   = KL1AM   + NT1AM(ISYCTR)
      KEND0   = KL2TP   + NT2SQ(ISYCTR)
      LWRK0   = LWORK   - KEND0
C
      KL1BY   = KEND0
      KL2TPBY   = KL1BY   + NT1AM(ISYMBY)
      KFOCK0  = KL2TPBY   + NT2SQ(ISYMBY)
      KFOCKY  = KFOCK0    + N2BST(ISYM0)
      KEND0   = KFOCKY    + N2BST(ISYMBY)
C
      KT1AMP0 = KEND0
      KEND1   = KT1AMP0 + NT1AMX
      LWRK1   = LWORK   - KEND1

C
C-------------------------------------
C     Read in lamdap and lamdh
C-------------------------------------
C
      CALL GET_LAMBDA0(WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),LWRK1)
C
C---------------------------------------------------------------------
C     Read zeroth-order AO Fock matrix from file and trasform it to 
C     lambda basis
C---------------------------------------------------------------------
C
      CALL GET_FOCK0(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),
     *               LWRK1)
C
C---------------------------------------------------------------------
C     Read the matrix the property integrals and trasform it to lambda 
C     basis
C---------------------------------------------------------------------
C
      IF (LISTBY(1:3).EQ.'L1 ') THEN
         CALL GET_FOCKX(WORK(KFOCKY),LABELB,IDLSTBY,ISYMBY,WORK(KLAMP0),
     *                  ISYM0,WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
      END IF
C
C-------------------------------------
C     Read T2 amplitudes 
C-------------------------------------
C
      IOPT = 2
      CALL GET_T1_T2(IOPT,.FALSE.,DUMMY,WORK(KT2TP),'R0',0,ISYM0,
     *                WORK(KEND1),LWRK1)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
     *    DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)
C
C-------------------------------------
C     Read L1 and L2 amplitudes 
C-------------------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1AM),WORK(KL2TP),LISTL,IDLSTL,
     *               ISYCTR,WORK(KEND1),LWRK1)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
     *    DDOT(NT2SQ(ISYCTR),WORK(KL2TP),1,WORK(KL2TP),1)

C
C-------------------------------------
C     Read L1BY and L2BY multipliers 
C-------------------------------------
C
      IOPT  = 3
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1BY),WORK(KL2TPBY),LISTBY,
     *               IDLSTBY,ISYMBY,WORK(KEND1),LWRK1)

*     write(lupri,*)'nrm KL1BY, LISTBY,IDLSTBY ',
*    * LISTBY,IDLSTBY,
*    * ddot(NT1AM(ISYMBY),WORK(KL1BY),1,WORK(KL1BY),1)

*     write(lupri,*)'nrm KL2TPBY, LISTBY,IDLSTBY ',
*    * LISTBY,IDLSTBY,
*    * ddot(NT2SQ(ISYMBY),WORK(KL2TPBY),1,WORK(KL2TPBY),1)
C
C---------------------------------------------------------------
C     Read canonical orbital energies and delete frozen orbitals 
C     in Fock diagonal, if required
C---------------------------------------------------------------
C
      CALL GET_ORBEN(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C--------------------------------------------
C     Sort the Fock matrix to get F(ck) block
C--------------------------------------------
C
      CALL SORT_FOCKCK(WORK(KFCKBA),WORK(KFOCK0),ISYM0)
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT If we want to sum the T3 amplitudes
COMMENT COMMENT COMMENT
C
      if (.false.) then
         kx3am  = kend1
         kend1 = kx3am + nrhft*nrhft*nrhft*nvirt*nvirt*nvirt
         call dzero(work(kx3am),nrhft*nrhft*nrhft*nvirt*nvirt*nvirt)
         lwrk0 = lwork - kend1
         if (lwrk0 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend1
            call quit('Insufficient space (T3) in CC3_XI_DEN')
         END IF
      endif
COMMENT COMMENT COMMENT
C      write(lupri,*) 'WMAT after dzero'
C      call print_pt3(work(kx3am),ISYMBY,4)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
      ISINT2Y = MULD2H(ISYMY,ISINT2)
      ISYFCKY = MULD2H(ISYMOP,ISYMY)

      KXIAJB   = KEND1
      KEND1   = KXIAJB  + NT2AM(ISYM0)

      KT3BOG1 = KEND1
      KT3BOL1 = KT3BOG1 + NTRAOC(ISYM0)
      KT3BOG2 = KT3BOL1 + NTRAOC(ISYM0)
      KT3BOL2 = KT3BOG2 + NTRAOC(ISYM0)
      KT3OG1  = KT3BOL2 + NTRAOC(ISYM0)
      KT3OG2 = KT3OG1  + NTRAOC(ISINT2)
      KLAMPY  = KT3OG2 + NTRAOC(ISINT2)
      KLAMHY  = KLAMPY  + NLAMDT
      KEND1   = KLAMHY  + NLAMDT
      LWRK1   = LWORK   - KEND1
C
      CALL DZERO(WORK(KLAMPY),NLAMDT)
      CALL DZERO(WORK(KLAMHY),NLAMDT)
C
      KFCKYCK    = KEND1
      KW3BXOG1   = KFCKYCK    + NT1AM(ISYFCKY)
      KW3BXOL1   = KW3BXOG1   + NTRAOC(ISYM0)
      KW3BXOGX1   = KW3BXOL1   + NTRAOC(ISYM0)
      KW3BXOLX1   = KW3BXOGX1   + NTRAOC(ISINT2Y)
      KEND1      = KW3BXOLX1   + NTRAOC(ISINT2Y)
      LWRK1      = LWORK      - KEND1
C
      CALL DZERO(WORK(KW3BXOGX1),NTRAOC(ISINT2Y))
      CALL DZERO(WORK(KW3BXOLX1),NTRAOC(ISINT2Y))
      CALL DZERO(WORK(KFCKYCK),NT1AM(ISYFCKY))
C
      KT1Y  = KEND1
      KEND2  = KT1Y + NT1AM(ISYMY)
      LWRK2   = LWORK  - KEND2
C
      CALL DZERO(WORK(KT1Y),NT1AM(ISYMY))
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_XI_DEN')
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
C---------------------------------------------------
C     Prepare to construct the occupied integrals...
C---------------------------------------------------
C
C        isint1  - symmetry of integrals in standard H, transformed
C                  with LambdaH_0
C        isint1y - symmetry of integrals in standard H, transformed
C                  with LambdaH_Y

      ISINT1  = 1
      ISINT1Y = MULD2H(ISINT1,ISYMY)
C
C--------------------------
C     Get Lambda_Y matrices
C--------------------------
C
      IF (.NOT.(LISTBY(1:3).EQ.'LE ') ) THEN
         CALL GET_LAMBDAX(WORK(KLAMPY),WORK(KLAMHY),LISTY,IDLSTY,ISYMY,
     *                    WORK(KLAMP0),WORK(KLAMH0),WORK(KEND2),LWRK2)
C
C------------------------------------------------------------------
C        Calculate the F^Y matrix (kc elements evaluated and stored 
C        as ck)
C------------------------------------------------------------------
C
         IOPT = 1
         CALL GET_T1_T2(IOPT,.FALSE.,WORK(KT1Y),DUMMY,LISTY,IDLSTY,
     *                  ISYMY,WORK(KEND2),LWRK2)
         CALL CC3LR_MFOCK(WORK(KFCKYCK),WORK(KT1Y),WORK(KXIAJB),ISYFCKY)
C
         ! From now on WORK(KEND1) is used again, since we do not need
         ! KT1Y amplitudes any more...
      END IF
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_0 multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BAR0(LUTOC,FNTOC,WORK(KLAMH0),ISYM0,WORK(KT3BOG1),
     *                   WORK(KT3BOL1),WORK(KT3BOG2),WORK(KT3BOL2),
     *                   WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3_0 amplitudes
C-----------------------------------------------------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,WORK(KLAMP0),ISINT2,WORK(KT3OG1),
     *                WORK(KT3OG2),WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_Y multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BARX(.FALSE.,
     *                   LUTOC,FNTOC,ISYMOP,WORK(KLAMH0),ISYM0,ISINT1,
     *                   WORK(KLAMHY),ISYMY,ISINT1Y,WORK(KW3BXOG1),
     *                   WORK(KW3BXOL1),WORK(KW3BXOGX1),WORK(KW3BXOLX1),
     *                   WORK(KEND1),LWRK1)
C
C---------------------------------------------
C     Open files for Tbar and W intermediates:
C---------------------------------------------
C
      LUT3    = -1
      LUWBMAT = -1

      CALL WOPEN2(LUT3,FNT3,64,0)
      CALL WOPEN2(LUWBMAT,FNWBMAT,64,0)
C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKB  = MULD2H(ISYMD,ISYM0)
         ISYCKBY  = MULD2H(ISYMD,ISYMY)
         ISCKB1 = MULD2H(ISINT1,ISYMD)
         ISCKB2 = MULD2H(ISYMD,ISYM0)
C
         DO D = 1,NVIR(ISYMD)
C
C           ------------------
C           Memory allocation.
C           ------------------
            KT3VDG1  = KEND1
            KT3VDG2  = KT3VDG1  + NCKATR(ISCKB2)
            KT3VDG3   = KT3VDG2  + NCKATR(ISCKB2)
            KEND1   = KT3VDG3 + NCKATR(ISCKB2)
C
            KT3BVDL1  = KEND1
            KT3BVDL2  = KT3BVDL1 + NCKATR(ISYCKB)
            KT3BVDL3  = KT3BVDL2 + NCKATR(ISYCKB)
            KEND3   = KT3BVDL3 + NCKATR(ISYCKB)
            LWRK3   = LWORK  - KEND3
           
            KT3BVDG1 = KEND3
            KT3BVDG2 = KT3BVDG1 + NCKATR(ISYCKB)
            KT3BVDG3 = KT3BVDG2 + NCKATR(ISYCKB)
            KEND3   = KT3BVDG3 + NCKATR(ISYCKB)
            LWRK3   = LWORK  - KEND3

            KW3BXVDG1  = KEND3
            KW3BXVDG2  = KW3BXVDG1  + NCKATR(ISYCKB)
            KW3BXVDL1  = KW3BXVDG2  + NCKATR(ISYCKB)
            KW3BXVDL2  = KW3BXVDL1  + NCKATR(ISYCKB)
            KEND3     = KW3BXVDL2  + NCKATR(ISYCKB)
            LWRK3     = LWORK     - KEND3

            KW3BXVDGX1  = KEND3
            KW3BXVDGX2  = KW3BXVDGX1  + NCKATR(ISYCKBY)
            KW3BXVDLX1  = KW3BXVDGX2  + NCKATR(ISYCKBY)
            KW3BXVDLX2  = KW3BXVDLX1  + NCKATR(ISYCKBY)
            KEND3     = KW3BXVDLX2  + NCKATR(ISYCKBY)
            LWRK3     = LWORK     - KEND3
C
            KINTVI = KEND3
            KTRVI6 = KINTVI + MAX(NCKA(ISYCKB),NCKA(ISYCKBY))
            KEND4  = KTRVI6 + NCKATR(ISYCKB) 
            LWRK4  = LWORK  - KEND4
           
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CC3_XI_DEN')
            END IF
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed D) which are required 
C           to calculate t3_0 amplitudes
C-----------------------------------------------------------------------
C
            CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                        WORK(KT3VDG1),WORK(KT3VDG2),WORK(KT3VDG3),
     *                        WORK(KLAMH0),ISYMD,D,WORK(KEND4),LWRK4)
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed D) which are required 
C           to calculate t3bar_0 multipliers
C-----------------------------------------------------------------------
C
            CALL INTVIR_T3BAR0_D(LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                           LUDKBC3,FNDKBC3,LU3VI,FN3VI,ISYM0,
     *                           WORK(KT3BVDL1),WORK(KT3BVDG1),
     *                           WORK(KT3BVDG2),WORK(KT3BVDL2),
     *                           WORK(KT3BVDG3),WORK(KT3BVDL3),
     *                           WORK(KLAMP0),ISYMD,D,WORK(KEND4),LWRK4)
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed D) which are required 
C           to calculate t3bar_X multipliers
C-----------------------------------------------------------------------
C
            CALL INTVIR_T3BARX_D(.FALSE.,
     *                           ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                           LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                           WORK(KW3BXVDGX1),WORK(KW3BXVDG1),
     *                           WORK(KW3BXVDGX2),WORK(KW3BXVDG2),
     *                           WORK(KW3BXVDLX1),WORK(KW3BXVDL1),
     *                           WORK(KW3BXVDLX2),WORK(KW3BXVDL2),
     *                           WORK(KLAMPY),ISYMY,WORK(KLAMP0),
     *                           ISYM0,ISYMD,D,WORK(KEND4),LWRK4)
C
            DO ISYMB = 1,NSYM

               ISYALJ  = MULD2H(ISYMB,ISYM0)
               ISYALJ2 = MULD2H(ISYMD,ISYM0)
C
               IF (LISTBY(1:3).EQ.'N2 ') THEN
                  ISYALJLE  = MULD2H(ISYMB,ISYCTR)
                  ISYALJ2LE = MULD2H(ISYMD,ISYCTR)
               END IF
C
               ISYALJBY  = MULD2H(ISYMB,ISYMBY)
               ISYALJDY = MULD2H(ISYMD,ISYMBY)
               ISYMBD  = MULD2H(ISYMD,ISYMB)
c              ISCKIJ  = MULD2H(ISYMBD,ISYCTR)
               ISCKIJ  = MULD2H(ISYMBD,ISYM0)
               ISWMAT  = MULD2H(ISCKIJ,ISYMBY)
               ISYCKD  = MULD2H(ISYM0,ISYMB)

C              Can use kend3 since we do not need the integrals anymore.
               KSMAT2     = KEND3
               KUMAT2     = KSMAT2    + NCKIJ(ISCKIJ)
               KDIAG      = KUMAT2    + NCKIJ(ISCKIJ)
               KDIAGW     = KDIAG     + NCKIJ(ISCKIJ)
               KINDSQ     = KDIAGW    + NCKIJ(ISWMAT)
               KINDSQW    = KINDSQ    + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX     = KINDSQW   + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEX2    = KINDEX    + (NCKI(ISYALJ)  - 1)/IRAT + 1
               KINDEXBY   = KINDEX2   + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KINDEXDY   = KINDEXBY  + (NCKI(ISYALJBY)  - 1)/IRAT + 1
               KTMAT      = KINDEXDY  + (NCKI(ISYALJDY) - 1)/IRAT + 1
               KT3MAT     = KTMAT     + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWMAT))
               KWMAT      = KT3MAT    + NCKIJ(ISCKIJ)
               KEND4      = KWMAT     + NCKIJ(ISWMAT)
               LWRK4      = LWORK     - KEND4
C
               IF (LISTBY(1:3).EQ.'N2 ') THEN
                  KINDEXLE = KEND4
                  KINDEX2LE = KINDEXLE + (NCKI(ISYALJLE)  - 1)/IRAT + 1
                  KEND4     = KINDEX2LE + (NCKI(ISYALJ2LE) - 1)/IRAT + 1
                  LWRK4      = LWORK     - KEND4
               ELSE 
                  KINDEXLE = KINDEX
                  KINDEX2LE = KINDEX2
               END IF

               KS3MAT   = KEND4
               KU3MAT   = KS3MAT  + NCKIJ(ISCKIJ)
               KS3MAT3  = KU3MAT  + NCKIJ(ISCKIJ)
               KU3MAT3  = KS3MAT3 + NCKIJ(ISCKIJ)
               KEND4    = KU3MAT3 + NCKIJ(ISCKIJ)

               KT3VBG1  = KEND4
               KT3VBG2  = KT3VBG1  + NCKATR(ISYCKD)
               KT3VBG3   = KT3VBG2  + NCKATR(ISYCKD)
               KEND4   = KT3VBG3 + NCKATR(ISYCKD)

               KT3BVBG1 = KEND4
               KT3BVBG2 = KT3BVBG1 + NCKATR(ISYCKD)
               KT3BVBG3 = KT3BVBG2 + NCKATR(ISYCKD)
               KEND4   = KT3BVBG3 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4

               KSMAT4  = KEND4
               KUMAT4  = KSMAT4  + NCKIJ(ISCKIJ)
               KT3BVBL1 = KUMAT4  + NCKIJ(ISCKIJ)
               KT3BVBL2 = KT3BVBL1 + NCKATR(ISYCKD)
               KT3BVBL3 = KT3BVBL2 + NCKATR(ISYCKD)
               KEND4   = KT3BVBL3 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4


               KINTVI  = KEND4
               KEND5   = KINTVI  + MAX(NCKA(ISYMB),NCKA(ISYCKD))
               LWRK5   = LWORK   - KEND5

               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CC3_XI_DEN')
               END IF
C
C
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG), WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)

C
C              -----------------------
C              Construct index arrays.
C              -----------------------
C
               LENSQ  = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               LENSQW  = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT)
               CALL CC3_INDEX(WORK(KINDEX),ISYALJ)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJ2)
               CALL CC3_INDEX(WORK(KINDEXBY),ISYALJBY)
               CALL CC3_INDEX(WORK(KINDEXDY),ISYALJDY)
C
               IF (LISTBY(1:3).EQ.'N2 ') THEN
                  CALL CC3_INDEX(WORK(KINDEXLE),ISYALJLE)
                  CALL CC3_INDEX(WORK(KINDEX2LE),ISYALJ2LE)
               END IF
C
               DO B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed B) which are required 
C           to calculate t3_0 amplitudes
C           (the same routine as in d-loop is used)
C-----------------------------------------------------------------------
C
                 CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                             WORK(KT3VBG1),WORK(KT3VBG2),
     *                             WORK(KT3VBG3),WORK(KLAMH0),ISYMB,B,
     *                             WORK(KEND5),LWRK5)

C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed B) which are required 
C           to calculate t3bar_0 multipliers
C           (the same routine as in d-loop is used)
C-----------------------------------------------------------------------
C
                  CALL INTVIR_T3BAR0_D(LU3FOPX,FN3FOPX,LU3FOP2X,
     *                                 FN3FOP2X,LUDKBC3,FNDKBC3,
     *                                 LU3VI,FN3VI,ISYM0,WORK(KT3BVBL1),
     *                                 WORK(KT3BVBG1),WORK(KT3BVBG2),
     *                                 WORK(KT3BVBL2),WORK(KT3BVBG3),
     *                                 WORK(KT3BVBL3),WORK(KLAMP0),
     *                                 ISYMB,B,WORK(KEND5),LWRK5)
C
C-----------------------------------------------------
C                 Get T3_BD amplitudes (using S and U)
C-----------------------------------------------------
C
                  CALL GET_T30_BD(ISYM0,ISINT2,WORK(KT2TP),ISYM0,
     *                            WORK(KT3MAT),WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KINDSQ),LENSQ,WORK(KS3MAT),
     *                            WORK(KT3VDG1),WORK(KT3VDG2),
     *                            WORK(KT3OG1),WORK(KINDEX),
     *                            WORK(KS3MAT3),WORK(KT3VBG1),
     *                            WORK(KT3VBG2),WORK(KINDEX2),
     *                            WORK(KU3MAT),WORK(KT3VDG3),
     *                            WORK(KT3OG2),WORK(KU3MAT3),
     *                            WORK(KT3VBG3),ISYMB,B,ISYMD,D,ISCKIJ,
     *                            WORK(KEND5),LWRK5)
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
c       call sum_pt3(work(KT3MAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),4)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C
C                 -------------------------------------
C                 Write T3 amplitudes as T3^D(ai,bj,l) to disc
C                 -------------------------------------
                  CALL WRITE_T3_DL(LUT3,FNT3,WORK(KT3MAT),ISYM0,
     *                             ISYMD,ISYMB,B)
C
C---------------------------------------------------------
C                 Get T3bar_BD multipliers (using S and U)
C---------------------------------------------------------
C
                  IF (LISTBY(1:3).EQ.'L1 ') THEN
                   CALL GET_T3BAR0_BD(ISYCTR,WORK(KL1AM),ISYCTR,
     *                                WORK(KL2TP),ISYCTR,WORK(KTMAT),
     *                                WORK(KFCKBA),WORK(KFOCKD),
     *                                WORK(KDIAG),WORK(KXIAJB),ISYM0,
     *                                ISYM0,WORK(KINDSQ),LENSQ,
     *                                WORK(KSMAT2),WORK(KT3BVDG1),
     *                                WORK(KT3BVDG2),WORK(KT3BVDL1),
     *                                WORK(KT3BVDL2),WORK(KT3BOG1),
     *                                WORK(KT3BOL1),WORK(KINDEX),
     *                                WORK(KSMAT4),WORK(KT3BVBG1),
     *                                WORK(KT3BVBG2),WORK(KT3BVBL1),
     *                                WORK(KT3BVBL2),WORK(KINDEX2),
     *                                WORK(KUMAT2),WORK(KT3BVDG3),
     *                                WORK(KT3BVDL3),WORK(KT3BOG2),
     *                                WORK(KT3BOL2),WORK(KUMAT4),
     *                                WORK(KT3BVBG3),WORK(KT3BVBL3),
     *                                ISYMB,B,ISYMD,D,ISCKIJ,
     *                                WORK(KEND5),LWRK5)
C
c      call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),4)
C
                  END IF
C
C----------------------------------------------------
C                 Get T3barX_BD multipliers (using W)
C----------------------------------------------------
C
                  IF (LISTBY(1:3).EQ.'L1 ') THEN
c     call dzero(WORK(KFCKBA),NT1AMX)
c
c     call dzero(WORK(KW3BXVDL2),NCKATR(ISYCKB))
c     call dzero(WORK(KW3BXVDL1),NCKATR(ISYCKB))
c     call dzero(WORK(KW3BXVDG2),NCKATR(ISYCKB))
c     call dzero(WORK(KW3BXVDG1),NCKATR(ISYCKB))
c
c     call dzero(WORK(KW3BXOL1),NTRAOC(ISYM0))
c     call dzero(WORK(KW3BXOG1),NTRAOC(ISYM0))
                  CALL GET_T3BARX_BD(.FALSE.,
     *                               WORK(KTMAT),ISCKIJ,WORK(KFOCKY),
     *                               ISYMBY,WORK(KWMAT),ISWMAT,
     *                               WORK(KL2TP),ISYCTR,WORK(KFCKYCK),
     *                               ISYFCKY,WORK(KW3BXVDLX2),
     *                               WORK(KW3BXVDLX1),WORK(KW3BXVDGX2),
     *                               WORK(KW3BXVDGX1),WORK(KW3BXOLX1),
     *                               WORK(KW3BXOGX1),ISINT2Y,
     *                               WORK(KINDEX),WORK(KINDEX2),
     *                               WORK(KINDSQW),LENSQW,WORK(KL2TPBY),
     *                               ISYMBY,WORK(KFCKBA),ISYM0,
     *                               WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                               WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                               WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                               ISINT2,WORK(KINDEXBY),
     *                               WORK(KINDEXDY),WORK(KL1BY),ISYMBY,
     *                               WORK(KXIAJB),ISINT1,-FREQ,
     *                               WORK(KDIAGW),WORK(KFOCKD),B,ISYMB,
     *                               D,ISYMD,ISYMBY,WORK(KEND5),LWRK5)
                 END IF
C
                 IF ( (LISTBY(1:3).EQ.'M1 ') 
     *               .OR. (LISTBY(1:3).EQ.'N2 ') ) THEN
                 CALL GET_M3BAR_BD(WORK(KTMAT),WORK(KWMAT),ISWMAT,
     *                             WORK(KL2TP),ISYCTR,WORK(KFCKYCK),
     *                             ISYFCKY,WORK(KW3BXVDLX2),
     *                             WORK(KW3BXVDLX1),WORK(KW3BXVDGX2),
     *                             WORK(KW3BXVDGX1),WORK(KW3BXOLX1),
     *                             WORK(KW3BXOGX1),ISINT2Y,
     *                             WORK(KINDEXLE),WORK(KINDEX2LE),
     *                             WORK(KINDSQW),LENSQW,WORK(KL2TPBY),
     *                             ISYMBY,WORK(KFCKBA),ISYM0,
     *                             WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                             WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                             WORK(KW3BXOL1),WORK(KW3BXOG1),ISINT2,
     *                             WORK(KINDEXBY),WORK(KINDEXDY),
     *                             WORK(KL1BY),ISYMBY,WORK(KXIAJB),
     *                             ISINT1,-FREQ,WORK(KDIAGW),
     *                             WORK(KFOCKD),B,ISYMB,D,ISYMD,ISYMBY,
     *                             WORK(KEND5),LWRK5)


*       call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *             ISWMAT,work(kx3am),4)
                 END IF 
C
                 IF (LISTBY(1:3).EQ.'LE ') THEN
                 CALL GET_L3BAR_BD(WORK(KTMAT),WORK(KWMAT),ISWMAT,
     *                             WORK(KINDSQW),LENSQW,WORK(KL2TPBY),
     *                             ISYMBY,WORK(KFCKBA),ISYM0,
     *                             WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                             WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                             WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                             ISINT2,WORK(KINDEXBY),WORK(KINDEXDY),
     *                             WORK(KL1BY),ISYMBY,WORK(KXIAJB),
     *                             ISINT1,-FREQ,WORK(KDIAGW),
     *                             WORK(KFOCKD),B,ISYMB,D,ISYMD,ISYMBY,
     *                             WORK(KEND5),LWRK5)
c       call sum_pt3(work(kwmat),isymb,b,isymd,d,
c    *             ISWMAT,work(kx3am),4)
                 END IF
C
C
C
C                 -------------------------------------
C                 Write WMAT as WMAT^D(ai,bj,l) to disc
C                 -------------------------------------
                  CALL WRITE_T3_DL(LUWBMAT,FNWBMAT,WORK(KWMAT),ISYMBY,
     *                             ISYMD,ISYMB,B)
C
C------------------------------------------------
C    Calculate density
C        D(ia) = <L3Y|[[E_ia,T2],T2]|HF>
C------------------------------------------------
C
                  CALL CC3_XI_DEN_IA(DIA,WORK(KWMAT),ISWMAT,WORK(KT2TP),
     *                               ISYM0,WORK(KINDSQW),LENSQW,B,ISYMB,
     *                               D,ISYMD,WORK(KEND5),LWRK5)
C
C
C------------------------------------------------
C    Calculate density
C        D(ia) = <L2Y|[E_ia,T3]|HF>
C------------------------------------------------
C
                  CALL CCFOP_ONED(DIA,WORK(KL2TPBY),ISYMBY,WORK(KS3MAT),
     *                          WORK(KT3MAT),ISYM0,WORK(KINDSQ),
     *                          LENSQ,WORK(KEND5),LWRK5,ISYMB,B,ISYMD,D)
C

               ENDDO   ! B
            ENDDO      ! ISYMB
C
C------------------------------------------------
C    Calculate densities
C         D(ab) = <L3Y|[E_ab,T3]|HF>
C         D(ij) = <L3Y|[E_ij,T3]|HF>
C------------------------------------------------
C
            CALL CC3_XI_DEN_ABIJ(DAB,DIJ,ISYM0,ISYMBY,LUT3,FNT3,
     *                           LUWBMAT,FNWBMAT,WORK(KEND5),LWRK5,
     *                           ISYMD)
C

         ENDDO       ! D
      ENDDO          ! ISYMD 
C
COMMENT COMMENT
c      write(lupri,*) 'cont to t3bx in CC3_XI_CONT'
c      call print_pt3(work(kx3am),ISYM0,4)
COMMENT COMMENT
C
C---------------------------------
C     Close the files
C---------------------------------
C
      CALL WCLOSE2(LUT3,FNT3,'DELETE')
      CALL WCLOSE2(LUWBMAT,FNWBMAT,'DELETE')
C

C
C-------------
C     End
C-------------
C
C
      CALL QEXIT('CC3_XI_DEN')
C
      RETURN
      END
C  /* Deck wbarbd_l1 */
      SUBROUTINE WBARBD_L1(T1AM,ISYMT1,TMAT,XIAJB,
     *                      ISINT1,
     *                      WMAT,WORK,LWORK,
     *                      INDSQ,LENSQ,ISYMB,B,ISYMC,C)
*---------------------------------------------------------------------*
*
*    Purpose: compute Tbar1^Y contribution to triples component of 
*    first-order multipliers vector:
*
*    <Tbar1^Y|[H_0^,tau3]|HF> = P^(abc)_(ijk) (  t1bar^Y(ai)*L(jbkc) 
*                                              - t1bar^Y(ak)*L(jbic)  )
*
*             Use W intermmediates: 
*
*    WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(jBkC)             
*                                  - T1(ak)*L(jBiC)
*                                  + T1(ai)*L(kCjB)
*                                  - T1(aj)*L(kCiB)
*
*    Written by Filip Pawlowski, Fall 2002, Aarhus
*            
*=====================================================================*
C

      IMPLICIT NONE
C
      INTEGER ISYMT1, ISINT1, LENSQ, ISYMB, ISYMC
      INTEGER ISYMBC, ISYRES, JSAIKJ, LENGTH, ISYMK, ISYMJ 
      INTEGER ISYMAI, ISYAIK, ISYMJK, ISYMCK, NBJ, NCK, ISYMBJ
      INTEGER NCKBJ, NBJCK, NAI, NAIKJ 
      INTEGER INDEX, INDSQ(LENSQ,6)
      INTEGER LWORK
C
#if defined (SYS_CRAY)
      REAL T1AM(*), TMAT(*), XIAJB(*) 
      REAL WMAT(*),WORK(LWORK)
      real xnormval,ddot
#else
      DOUBLE PRECISION T1AM(*), TMAT(*), XIAJB(*) 
      DOUBLE PRECISION WMAT(*),WORK(LWORK)  
      double precision xnormval,ddot
#endif
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('WBARBD_L1')
C
      ISYRES  = MULD2H(ISYMT1,ISINT1)
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYRES,ISYMBC)
      LENGTH = NCKIJ(JSAIKJ)
C
C--------------------------------------------
C     First contribution from both T1 terms
C
C    WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(jBkC)             
C                                  - T1(ak)*L(jBiC)
C
C--------------------------------------------
C
      ISYMJK = MULD2H(ISYMBC,ISINT1)
C
C---------------------------------------
C     Contract the integrals with T1.
C---------------------------------------
C
      CALL DZERO(TMAT,LENGTH)
C
      ISYMAI = ISYMT1
      DO ISYMJ = 1, NSYM
         ISYMK  = MULD2H(ISYMJK,ISYMJ)
         ISYAIK = MULD2H(ISYMK,ISYMAI)
         ISYMCK = MULD2H(ISYMC,ISYMK)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
         DO J = 1, NRHF(ISYMJ)
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            DO K = 1, NRHF(ISYMK)
C
               NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               NCKBJ = IT2AM(ISYMCK,ISYMBJ) + INDEX(NCK,NBJ)
C
               DO NAI = 1, NT1AM(ISYMAI)
C
                  NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK)
     *                  + NT1AM(ISYMAI)*(K - 1) + NAI
C
                  TMAT(NAIKJ) = T1AM(NAI)*XIAJB(NCKBJ)
C
               ENDDO
            ENDDO
         ENDDO
C
      ENDDO
C
C----------------------------------------
C     Sum the result into WMAT.
C----------------------------------------
C
      JSAIKJ = MULD2H(ISYMAI,ISYMJK)
      DO I = 1, NCKIJ(JSAIKJ)
C         First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(jBkC)             
          WMAT(I) = WMAT(I) + TMAT(I)
C         Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T1(ak)*L(jBiC)
          WMAT(I) = WMAT(I) - TMAT(INDSQ(I,1))
      ENDDO
C
C--------------------------------------------
C    Second contribution from both T1 terms
C
C    WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(kCjB)
C                                  - T1(aj)*L(kCiB)
C
C
C--------------------------------------------
C
      ISYMJK = MULD2H(ISYMBC,ISINT1)
C
C---------------------------------------
C     Contract the integrals with T1.
C---------------------------------------
C
      CALL DZERO(TMAT,LENGTH)
C
      ISYMAI = ISYMT1
      DO ISYMK = 1, NSYM
         ISYAIK = MULD2H(ISYMK,ISYMAI)
         ISYMJ  = MULD2H(ISYMJK,ISYMK)
         ISYMCK = MULD2H(ISYMC,ISYMK)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
         DO K = 1, NRHF(ISYMK)
C
            NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
            DO J = 1, NRHF(ISYMJ)
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               NBJCK = IT2AM(ISYMBJ,ISYMCK) + INDEX(NBJ,NCK)
C
               DO NAI = 1, NT1AM(ISYMAI)
C
                  NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK)
     *                  + NT1AM(ISYMAI)*(K-1) + NAI
C
                  TMAT(NAIKJ) = T1AM(NAI)*XIAJB(NBJCK)
C
               ENDDO
            ENDDO
         ENDDO
C
      ENDDO
c
C
C----------------------------------------
C     Sum the result into WMAT.
C----------------------------------------
C
      JSAIKJ = MULD2H(ISYMAI,ISYMJK)
      DO I = 1, NCKIJ(JSAIKJ)
C         First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(kCjB)
          WMAT(I) = WMAT(I) + TMAT(I)
C         Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T1(aj)*L(kCiB)
          WMAT(I) = WMAT(I) - TMAT(INDSQ(I,5))
      ENDDO
C
      CALL QEXIT('WBARBD_L1')
C
      RETURN
      END
C
C  /* Deck cc3_xi_den_abij */
      SUBROUTINE CC3_XI_DEN_ABIJ(DAB,DIJ,ISYM0,ISYFKY,LUT3,FNT3,LUWBMAT,
     *                      FNWBMAT,WORK,LWORK,ISYMD)
C
C=========================================================================
C    Calculate the densities 
C
C     (1)    D(ab) = <L3Y|[E_ab,T3]|HF> 
C
C            D(ab) = 1/2 tbarY^(dea)_(lmn) t^(deb)_(lmn) 
C
C                  = ( 1/2 W^de(anml) + W^da(emnl) ) t^(deb)_(lmn)
C                  =       Wtilde_bar^DL(em,aN)      T3^DL(em,bN)
C
C
C     (2)    D(ij) = <L3Y|[E_ij,T3]|HF> 
C
C            D(ij) = - 1/2 tbarY^(def)_(lmj) t^(def)_(lmi) 
C
C                  = ( 1/2 W^de(fjml) + W^df(emjl) ) t^(deb)_(lmi)
C                  =       Wtilde_bar^DL(em,aJ)      T3^DL(em,fi)
C
C
C     Common notation for (1) and (2):
C
C                   a -> f     j -> n
C                   ======     ======
C
C     Outside the densities D(ab) snd D(ij) are to be contracted 
C     with X_(ab) and X_(ij), respectively
C
C
C    Written by Filip Pawlowski, Fall 2002, Aarhus
C=========================================================================
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cc3t3d.h>
C
      INTEGER ISYM0,ISYFKY,LUT3,LUWBMAT,LWORK,ISYMD
      INTEGER ISYML,ISYMDL,ISWMAT,ISYMT3,ISYMN,ISYEMF,ISYMBN,ISYMEM
      INTEGER ISYMB,ISYMF,ISYMFI,ISYMI,ISYEMB,ISYMFN
      INTEGER KT3,KWMAT,KEND1,LWRK1
      INTEGER KOFF1,KOFF2,KOFF3,KBN,KFN
      INTEGER NTOTEM,NTOTF,NNEMF
      INTEGER IADR
C
      CHARACTER*(*) FNT3,FNWBMAT
C
#if defined (SYS_CRAY)
      REAL DAB(*),DIJ(*),WORK(LWORK),ONE,HALF
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION DAB(*),DIJ(*),WORK(LWORK),ONE,HALF
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C
      PARAMETER(ONE = 1.0D0, HALF = 0.5D0)
C
      CALL QENTER('CC3_XI_DEN_ABIJ')
C
      DO ISYML = 1,NSYM
C
         ISYMDL = MULD2H(ISYMD,ISYML)
         ISWMAT = MULD2H(ISYFKY,ISYMDL)
         ISYMT3 = MULD2H(ISYM0,ISYMDL)
C
         KT3  = 1
         KWMAT  = KT3 + NT2SQ(ISYMT3)
         KEND1 = KWMAT + NT2SQ(ISWMAT)
         LWRK1  = LWORK - KEND1

C
         IF ( LWRK1 .LT. 0 ) THEN
           CALL QUIT('Out of memory in CC3_XI_DEN_ABIJ (x)')
         ENDIF
C
         DO L = 1, NRHF(ISYML)
C
C           --------------------------------------------
C           Read T3 amplitudes from file:
C           --------------------------------------------
C
            IADR = ISWTL(ISYMT3,ISYML) + NT2SQ(ISYMT3)*(L-1) + 1
            CALL GETWA2(LUT3,FNT3,WORK(KT3),IADR,NT2SQ(ISYMT3))
C
C           --------------------------------------------
C           Read WMAT_bar from file and generate WMAT-tilde:
C           --------------------------------------------
C
            IADR = ISWTL(ISWMAT,ISYML) + NT2SQ(ISWMAT)*(L-1) + 1
            CALL GETWA2(LUWBMAT,FNWBMAT,WORK(KWMAT),IADR,NT2SQ(ISWMAT))

            CALL CC_T2MOD(WORK(KWMAT),ISWMAT,HALF)

C
C           -----------------------------------
C           Loop over outermost occupied index:
C           -----------------------------------
            DO ISYMN = 1, NSYM
               ISYEMF = MULD2H(ISWMAT,ISYMN)
               ISYEMB = MULD2H(ISYMT3,ISYMN)
C
               DO N = 1, NRHF(ISYMN)
C
C                 -------------------------------------------------------
C                 D(fb) <- D(fb)+ sum_em Wtilde_bar^DL(em,fN) T3^DL(em,bN):
C                 -------------------------------------------------------
                  DO ISYMEM = 1, NSYM
                     ISYMB  = MULD2H(ISYEMB,ISYMEM)
                     ISYMF  = MULD2H(ISYEMF,ISYMEM)
                     ISYMFN = MULD2H(ISYMF,ISYMN)
                     ISYMBN = MULD2H(ISYMB,ISYMN)

                     KFN    = IT1AM(ISYMF,ISYMN)+NVIR(ISYMF)*(N-1)+1
                     KOFF1  = KWMAT + IT2SQ(ISYMEM,ISYMFN)
     *                              + NT1AM(ISYMEM)*(KFN-1)
                     KBN    = IT1AM(ISYMB,ISYMN)+NVIR(ISYMB)*(N-1)+1
                     KOFF2  = KT3   + IT2SQ(ISYMEM,ISYMBN)
     *                              + NT1AM(ISYMEM)*(KBN-1)
                     KOFF3  = IMATAB(ISYMF,ISYMB) + 1


                     NTOTEM = MAX(NT1AM(ISYMEM),1)
                     NTOTF  = MAX(NVIR(ISYMF),1)

                     CALL DGEMM('T','N',NVIR(ISYMF),NVIR(ISYMB),
     *                          NT1AM(ISYMEM),ONE,WORK(KOFF1),NTOTEM,
     *                          WORK(KOFF2),NTOTEM,ONE,DAB(KOFF3),NTOTF)

                  END DO ! ISYMEM
C                   -------------------------------------------------------
C                   D(iN) <- D(iN)- sum_emf Wtilde_bar^DL(em,fN) t^DL(em,fi):
C                   -------------------------------------------------------
                  DO ISYMEM = 1, NSYM
                     ISYMFI = MULD2H(ISYMT3,ISYMEM)
                     ISYMF  = MULD2H(ISYEMF,ISYMEM)
                     ISYMI  = MULD2H(ISYMFI,ISYMF)
                     ISYMFN = MULD2H(ISYMF,ISYMN)

C
                     KOFF1  = KT3   + IT2SQ(ISYMEM,ISYMFI) 
     *                              + NT1AM(ISYMEM)*IT1AM(ISYMF,ISYMI)

                     KFN    = IT1AM(ISYMF,ISYMN)+NVIR(ISYMF)*(N-1)+1
                     KOFF2  = KWMAT + IT2SQ(ISYMEM,ISYMFN)
     *                              + NT1AM(ISYMEM)*(KFN-1)

                     KOFF3  = IMATIJ(ISYMI,ISYMN) 
     *                              + NRHF(ISYMI)*(N-1) + 1

                     NNEMF  = MAX(NT1AM(ISYMEM)*NVIR(ISYMF),1)

                     CALL DGEMV('T',NT1AM(ISYMEM)*NVIR(ISYMF),
     *                         NRHF(ISYMI),-ONE,WORK(KOFF1),NNEMF,
     *                         WORK(KOFF2),1,ONE,DIJ(KOFF3),1)

                  END DO ! ISYMFI

               END DO ! N 
            END DO    ! ISYMN 
         END DO       ! L 
      END DO          ! ISYML 
C
      CALL QEXIT('CC3_XI_DEN_ABIJ')
C
      RETURN
      END
C
C  /* Deck cc3_xi_den_ia */
      SUBROUTINE CC3_XI_DEN_IA(DIA,WMAT,ISWMAT,T2TP,ISYMT2,INDSQ,LENSQ,
     *                              IB,ISYMIB,ID,ISYMID,WORK,LWORK)
C
C=========================================================================
C    Calculate the density 
C
C            D(ia) = <L3Y|[[E_ia,T2],T2]|HF> 
C
C            D(ia) = - tbarY^(dbc)_(jlk) t^(db)_(ji) t^(ca)_(kl) 
C
C                  = ( - W^BC(djkl) T2TP(djiB) T2TP(alkC)
C                      - W^DB(cklj) T2TP(DjiB) T2TP(ckla)
C                      - W^CD(bljk) T2TP(bijD) T2TP(alkC) )  
C
C            Density stored as AI:  DS(ai) = D(ia)
C
C            Outside the density DS(ai) is to be contracted with X_(ia)
C
C
C    Written by Filip Pawlowski, Fall 2002, Aarhus
C=========================================================================
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISWMAT,ISYMT2,IB,ISYMIB,ID,ISYMID,LWORK
      INTEGER ISYMB,ISYMC,ISYDJI,ISYALK,ISYKLI,ISYMDJ,ISYMKL,ISYMI
      INTEGER ISYMK,ISYMAL,ISYMA,ISYML,ISYMD,ISYMBD,ISYMJI,ISYCKLI
      INTEGER ISYMJ,ISYCKL,ISYBIJ,ISYLKI,ISYMBI,ISYMBJ,ISYMLK,ISYMAK
      INTEGER LENGTH,KTMAT,KEND1,LWRK1,KGKLI,KT2KLA
      INTEGER KT2JI,KGCKLI,KGLKI,KT2BJI,KEND2,LWRK2,KT2ALK
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER NTOTDJ,NTOTKL,NTOTA,NTOTJ,NTOTCKL,NTOTBJ,NTOTLK
      INTEGER LENSQ,INDSQ(LENSQ,6)
      INTEGER ISYMIBID,ISYDEN
C
#if defined (SYS_CRAY)
      REAL DIA(*),WMAT(*),T2TP(*),WORK(LWORK),ONE,TWO
#else
      DOUBLE PRECISION DIA(*),WMAT(*),T2TP(*),WORK(*),ONE,TWO
#endif
C
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CC3_XI_DEN_IA')
C
C
      LENGTH = NCKIJ(ISWMAT)
C
      IF (LENSQ .NE. LENGTH) THEN
         WRITE(LUPRI,*) 'LENSQ = ', LENSQ
         WRITE(LUPRI,*) 'LENGTH = ', LENGTH
         CALL QUIT('Inconsistency in length in CC3_XI_DEN_IA')
      END IF
C
      KTMAT  = 1
      KEND1  = KTMAT + NCKIJ(ISWMAT)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_XI_DEN_IA (0)')
      ENDIF
C
      ISYMIBID = MULD2H(ISYMIB,ISYMID)
      ISYDEN   = MULD2H(ISWMAT,ISYMIBID)
C
C
C=========================================================================
C    calculate  first contribution: 
C                 DS(ai) = DS(ai) - W^BC(djkl) T2TP(djiB) T2TP(alkC)
C                                                         
C
C                                          G(kl,i)        T^C(kl,a)
C=========================================================================
C
      B = IB
      C = ID
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
C------------------------
C Set symmetry flags
C------------------------
C
      ISYDJI = MULD2H(ISYMT2,ISYMB)
      ISYALK = MULD2H(ISYMT2,ISYMC)
      ISYKLI = MULD2H(ISWMAT,ISYDJI)
C
C-----------------------
C Memory allocation
C-----------------------
C
      KGKLI  = KEND1
      KT2KLA = KGKLI  + NMAIJK(ISYKLI)
      KEND2  = KT2KLA + NCKI(ISYALK)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_XI_DEN_IA (1)')
      ENDIF
C
C---------------------
C Initialize
C---------------------
C
      CALL DZERO(WORK(KTMAT),LENGTH)
      CALL DZERO(WORK(KT2KLA),NCKI(ISYALK))
      CALL DZERO(WORK(KGKLI),NMAIJK(ISYKLI))

C
C-------------------------------------------------------------------------
C     If symmetry, sort W^BC(djkl) as W^BC(dj,kl)
C-------------------------------------------------------------------------
C
      DO I = 1, LENGTH
         WORK(KTMAT-1+I) = WMAT(I)
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK2 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_XI_DEN_IA'
     *                 //' (CC_GATHER-1)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND2),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND2),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------------------------------------------
C     Calculate G(kl,i) = - W^BC(dj,kl) T2TP^B(dj,i)
C-------------------------------------------------------------------------
C
      DO ISYMDJ = 1,NSYM
C
         ISYMKL = MULD2H(ISWMAT,ISYMDJ)
         ISYMI  = MULD2H(ISYDJI,ISYMDJ)
C
         KOFF1  = KTMAT + ISAIKL(ISYMDJ,ISYMKL)  
         KOFF2  = IT2SP(ISYDJI,ISYMB)
     *          + NCKI(ISYDJI)*(B-1)
     *          + ISAIK(ISYMDJ,ISYMI)
     *          + 1
         KOFF3 = KGKLI  + IMAIJK(ISYMKL,ISYMI)
C
         NTOTDJ = MAX(NT1AM(ISYMDJ),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         CALL DGEMM('T','N',NMATIJ(ISYMKL),NRHF(ISYMI),
     *              NT1AM(ISYMDJ),-ONE,WORK(KOFF1),NTOTDJ,
     *              T2TP(KOFF2),NTOTDJ,ONE,WORK(KOFF3),NTOTKL)
C
      END DO ! ISYMDJ
C
C-------------------------------------
C     Sort T2   T^C(kl,a) = T2TP(alkC)
C-------------------------------------
C
      DO ISYMK = 1, NSYM
         ISYMAL = MULD2H(ISYALK,ISYMK)
         DO ISYMA = 1, NSYM
            ISYML = MULD2H(ISYMAL,ISYMA)
            ISYMKL = MULD2H(ISYMK,ISYML)
C
            DO K = 1, NRHF(ISYMK)
               DO L = 1, NRHF(ISYML)
                  KOFF1 = IT2SP(ISYALK,ISYMC)
     *                  + NCKI(ISYALK)*(C-1)
     *                  + ICKI(ISYMAL,ISYMK)
     *                  + NT1AM(ISYMAL)*(K-1)
     *                  + IT1AM(ISYMA,ISYML)
     *                  + NVIR(ISYMA)*(L-1)
     *                  + 1
                  KOFF2 = KT2KLA - 1           
     *                  + IMAIJA(ISYMKL,ISYMA)
     *                  + IMATIJ(ISYMK,ISYML)
     *                  + NRHF(ISYMK)*(L-1)
     *                  + K
C
                  CALL DCOPY(NVIR(ISYMA),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMKL))
               ENDDO ! L
            ENDDO    ! K
         ENDDO       ! ISYMA
      ENDDO          ! ISYMK
C
C-------------------------------------------------------------------------
C     Calculate D(ia) = T^C(kl,a)  G(kl,i) 
C-------------------------------------------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMKL = MULD2H(ISYALK,ISYMA)
         ISYMI  = MULD2H(ISYKLI,ISYMKL)
C
         KOFF1 = KT2KLA + IMAIJA(ISYMKL,ISYMA)
         KOFF2 = KGKLI  + IMAIJK(ISYMKL,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMKL),
     *              ONE,WORK(KOFF1),NTOTKL,
     *              WORK(KOFF2),NTOTKL,ONE,DIA(KOFF3),NTOTA)
C
      END DO ! ISYMA

C
C=========================================================================
C    calculate  second contribution: 
C                 DS(ai) = DS(ai) - W^DB(cklj) T2TP(DjiB) T2TP(ckla)
C                                                         
C                                          G(ckl,i)        T(ckl,a)
C=========================================================================
C
      D = IB
      B = ID
      ISYMD = ISYMIB
      ISYMB = ISYMID
C
C------------------------
C Set symmetry flags
C------------------------
C
      ISYMBD  = MULD2H(ISYMB,ISYMD)
      ISYMJI  = MULD2H(ISYMBD,ISYMT2)
      ISYCKLI = MULD2H(ISWMAT,ISYMJI)
C
C-----------------------
C Memory allocation
C-----------------------
C
      KT2JI  = KEND1
      KGCKLI   = KT2JI + NMATIJ(ISYMJI)
      KEND2   = KGCKLI  + NCKIJ(ISYCKLI)
      LWRK2   = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_XI_DEN_IA (2)')
      ENDIF
C
C---------------------
C Initialize
C---------------------
C
      CALL DZERO(WORK(KT2JI),NMATIJ(ISYMJI))
      CALL DZERO(WORK(KGCKLI),NCKIJ(ISYCKLI))
C
C ---------------------------------
C     Sort T2TP(DjiB) as T^{DB}_{ji} 
C ---------------------------------
C
      ISYDJI = MULD2H(ISYMJI,ISYMD)
      DO ISYMI = 1,NSYM
         ISYMDJ = MULD2H(ISYDJI,ISYMI)
         ISYMJ = MULD2H(ISYMJI,ISYMI)
         DO I = 1,NRHF(ISYMI)
            DO J = 1,NRHF(ISYMJ)
               KOFF1 = IT2SP(ISYDJI,ISYMB)
     *                  + NCKI(ISYDJI)*(B-1)
     *                  + ICKI(ISYMDJ,ISYMI)
     *                  + NT1AM(ISYMDJ)*(I-1)
     *                  + IT1AM(ISYMD,ISYMJ)
     *                  + NVIR(ISYMD)*(J-1)
     *                  + D
               KOFF2 = IMATIJ(ISYMJ,ISYMI)
     *                  + NRHF(ISYMJ)*(I-1)
     *                  + J
               WORK(KT2JI-1+KOFF2) = T2TP(KOFF1)
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------------------------------------------------
C     Calculate G(ckl,i) = - W^DB(cklj) T^{DB}_{ji}
C-------------------------------------------------------------------------
C
      DO ISYMI = 1,NSYM
C
         ISYMJ   = MULD2H(ISYMJI,ISYMI)
         ISYCKL  = MULD2H(ISWMAT,ISYMJ)
C
         KOFF1  = ISAIKJ(ISYCKL,ISYMJ) + 1
         KOFF2  = IMATIJ(ISYMJ,ISYMI)  + KT2JI
         KOFF3  = ISAIKJ(ISYCKL,ISYMI) + KGCKLI
C
         NTOTJ      =  MAX(1,NRHF(ISYMJ))
         NTOTCKL    =  MAX(1,NCKI(ISYCKL))
C
         CALL DGEMM('N','N',NCKI(ISYCKL),NRHF(ISYMI),
     *               NRHF(ISYMJ),-ONE,WMAT(KOFF1),NTOTCKL,
     *               WORK(KOFF2),NTOTJ,ONE,WORK(KOFF3),NTOTCKL)

      END DO ! ISYMI
C
C-------------------------------------------------------------------------
C       DS(ai) = DS(ai) + T(ckl,a)  G(ckl,i) 
C-------------------------------------------------------------------------
C
      DO ISYMA = 1,NSYM
C
         ISYCKL  = MULD2H(ISYMT2,ISYMA)
         ISYMI  = MULD2H(ISYCKLI,ISYCKL)
C
         KOFF1  = IT2SP(ISYCKL,ISYMA) + 1
         KOFF2  = ISAIKJ(ISYCKL,ISYMI) + KGCKLI
         KOFF3  = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTCKL   = MAX(1,NCKI(ISYCKL))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *               NCKI(ISYCKL),ONE,T2TP(KOFF1),NTOTCKL,
     *               WORK(KOFF2),NTOTCKL,ONE,DIA(KOFF3),NTOTA)
      END DO ! ISYMA

C
C========================================================================= C    calculate  third contribution: 
C                 DS(ai) = DS(ai) - W^CD(bljk) T2TP(bijD) T2TP(alkC)
C
C                                          G(lk,i)        T^C(lk,a)
C=========================================================================
C
      C = IB
      D = ID
      ISYMC = ISYMIB
      ISYMD = ISYMID
C
C------------------------
C Set symmetry flags
C------------------------
C
      ISYBIJ = MULD2H(ISYMT2,ISYMD)
      ISYALK = MULD2H(ISYMT2,ISYMC)
      ISYLKI = MULD2H(ISWMAT,ISYBIJ)
C
C-----------------------
C Memory allocation
C-----------------------
C
      KGLKI  = KEND1
      KEND2  = KGLKI  + NMAIJK(ISYLKI)
C
      KT2BJI = KEND2
      KEND2  = KT2BJI + NCKI(ISYBIJ)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_XI_DEN_IA (3)')
      ENDIF
C
C---------------------
C Initialize
C---------------------
C
      CALL DZERO(WORK(KTMAT),LENGTH)
      CALL DZERO(WORK(KT2BJI),NCKI(ISYBIJ))
      CALL DZERO(WORK(KGLKI),NMAIJK(ISYLKI))
C
C-------------------------------------------------------------------------
C     Sort W^CD(bljk) as W^CD(bjlk)
C-------------------------------------------------------------------------
C
      DO I = 1, LENGTH
         WORK(KTMAT-1+I) = WMAT(INDSQ(I,1))
      ENDDO
C
C-------------------------------------------------------------------------
C     If symmetry, sort W^CD(bjlk) as W^CD(bj,lk)
C-------------------------------------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK2 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_XI_DEN_IA'
     *                 //' (CC_GATHER-2)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND2),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND2),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Sort T2   T^D(bj,i) = T2TP(bijD)
C-------------------------------------
C
      DO ISYMJ = 1, NSYM
         ISYMBI = MULD2H(ISYBIJ,ISYMJ)
         DO ISYMB = 1, NSYM
            ISYMI = MULD2H(ISYMBI,ISYMB)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO J = 1, NRHF(ISYMJ)
               DO I = 1, NRHF(ISYMI)
                  KOFF1 = IT2SP(ISYBIJ,ISYMD)
     *                  + NCKI(ISYBIJ)*(D-1)
     *                  + ICKI(ISYMBI,ISYMJ)
     *                  + NT1AM(ISYMBI)*(J-1)
     *                  + IT1AM(ISYMB,ISYMI)
     *                  + NVIR(ISYMB)*(I-1)
     *                  + 1
                  KOFF2 = KT2BJI        
     *                  + ICKI(ISYMBJ,ISYMI)
     *                  + NT1AM(ISYMBJ)*(I-1)
     *                  + IT1AM(ISYMB,ISYMJ)
     *                  + NVIR(ISYMB)*(J-1)
C
                  CALL DCOPY(NVIR(ISYMB),T2TP(KOFF1),1,
     *                       WORK(KOFF2),1) 
               ENDDO ! I
            ENDDO    ! J
         ENDDO       ! ISYMB
      ENDDO          ! ISYMJ
C    
C-------------------------------------------------------------------------
C     Calculate G(lk,i) = - W^CD(bj,lk) T^D(bj,i)     
C-------------------------------------------------------------------------
C
      DO ISYMBJ = 1,NSYM
C
         ISYMLK = MULD2H(ISWMAT,ISYMBJ)
         ISYMI  = MULD2H(ISYBIJ,ISYMBJ)
C
         KOFF1  = KTMAT  + ISAIKL(ISYMBJ,ISYMLK) 
         KOFF2  = KT2BJI + ICKI(ISYMBJ,ISYMI)
         KOFF3  = KGLKI  + IMAIJK(ISYMLK,ISYMI)
C
         NTOTBJ = MAX(NT1AM(ISYMBJ),1)
         NTOTLK = MAX(NMATIJ(ISYMLK),1)
C
         CALL DGEMM('T','N',NMATIJ(ISYMLK),NRHF(ISYMI),
     *              NT1AM(ISYMBJ),-ONE,WORK(KOFF1),NTOTBJ,
     *              WORK(KOFF2),NTOTBJ,ONE,WORK(KOFF3),NTOTLK)
C

      END DO ! ISYMBJ
C
C-------------------------------------
C     Sort T2   T^C(lk,a) = T2TP(alkC)
C-------------------------------------
C
      KT2ALK = KEND2
      KEND2  = KT2ALK + NCKI(ISYALK)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_XI_DEN_IA (4)')
      ENDIF
C
      CALL DZERO(WORK(KT2ALK),NCKI(ISYALK))
C
      DO ISYML = 1, NSYM
         ISYMAK = MULD2H(ISYALK,ISYML)
         DO ISYMA = 1, NSYM
            ISYMK = MULD2H(ISYMAK,ISYMA)
            ISYMLK = MULD2H(ISYML,ISYMK)
            ISYMAL = MULD2H(ISYML,ISYMA)
C
            DO L = 1, NRHF(ISYML)
               DO K = 1, NRHF(ISYMK)
                  KOFF1 = IT2SP(ISYALK,ISYMC)
     *                  + NCKI(ISYALK)*(C-1)
     *                  + ICKI(ISYMAL,ISYMK)
     *                  + NT1AM(ISYMAL)*(K-1)
     *                  + IT1AM(ISYMA,ISYML)
     *                  + NVIR(ISYMA)*(L-1)
     *                  + 1
                  KOFF2 = KT2ALK - 1
     *                  + IMAIJA(ISYMLK,ISYMA)
     *                  + IMATIJ(ISYML,ISYMK)
     *                  + NRHF(ISYML)*(K-1)
     *                  + L
C
                  CALL DCOPY(NVIR(ISYMA),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMLK))
               ENDDO ! K
            ENDDO    ! L
         ENDDO       ! ISYMA
      ENDDO          ! ISYML
C
C-------------------------------------------------------------------------
C     Calculate DS(ai) = T^C(lk,a)  G(lk,i) 
C-------------------------------------------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMLK = MULD2H(ISYALK,ISYMA)
         ISYMI  = MULD2H(ISYLKI,ISYMLK)
C
         KOFF1 = KT2ALK + IMAIJA(ISYMLK,ISYMA)
         KOFF2 = KGLKI  + IMAIJK(ISYMLK,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTLK = MAX(NMATIJ(ISYMLK),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMLK),
     *              ONE,WORK(KOFF1),NTOTLK,
     *              WORK(KOFF2),NTOTLK,ONE,DIA(KOFF3),NTOTA)
C
      END DO ! ISYMA

      CALL QEXIT('CC3_XI_DEN_IA')
C
      RETURN
      END

C  /* Deck print_matij */
      SUBROUTINE PRINT_MATIJ(MATIJ,ISYMIJ)
**************************************************************
*                                                            *
*     The printing routine, mainly for debugging purposes.   *
*                                                            *
*     Print M(i,j) matrix symmetry-wise                      *
*                                                            *
*     Filip Pawlowski, 04-Dec-2002, Aarhus                   *
*                                                            *
**************************************************************
C
      IMPLICIT NONE 
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      INTEGER ISYMIJ,ISYMJ,ISYMI
      INTEGER KOFF1
C
#if defined (SYS_CRAY)
      REAL MATIJ(*)
#else
      DOUBLE PRECISION MATIJ(*)
#endif
C
      CALL QENTER('PRINT_MATIJ')
C
      DO ISYMJ = 1, NSYM
         ISYMI = MULD2H(ISYMIJ,ISYMJ)
C
         KOFF1 = IMATIJ(ISYMI,ISYMJ) + 1
C
         WRITE(LUPRI,*) 'Symmetry block number(i,j): ',ISYMI,ISYMJ
C
         CALL OUTPUT_E(MATIJ(KOFF1),1,NRHF(ISYMI),1,NRHF(ISYMJ),
     *               NRHF(ISYMI),NRHF(ISYMJ),1,LUPRI)
C
         WRITE(LUPRI,*) ' '
C
      END DO ! ISYMJ
C
      CALL QEXIT('PRINT_MATIJ')
C
      RETURN
      END

C  /* Deck print_matab */
      SUBROUTINE PRINT_MATAB(MATAB,ISYMAB)
**************************************************************
*                                                            *
*     The printing routine, mainly for debugging purposes.   *
*                                                            *
*     Print M(a,b) matrix symmetry-wise                      *
*                                                            *
*     Filip Pawlowski, 04-Dec-2002, Aarhus                   *
*                                                            *
**************************************************************
C
      IMPLICIT NONE 
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

C
      INTEGER ISYMAB,ISYMB,ISYMA
      INTEGER KOFF1
C
#if defined (SYS_CRAY)
      REAL MATAB(*)
#else
      DOUBLE PRECISION MATAB(*)
#endif
C
      CALL QENTER('PRINT_MATAB')
C
      DO ISYMB = 1, NSYM
         ISYMA = MULD2H(ISYMAB,ISYMB)
C
         KOFF1 = IMATAB(ISYMA,ISYMB) + 1
C
         WRITE(LUPRI,*) 'Symmetry block number(a,b): ',ISYMA,ISYMB
C
         CALL OUTPUT_E(MATAB(KOFF1),1,NVIR(ISYMA),1,NVIR(ISYMB),
     *               NVIR(ISYMA),NVIR(ISYMB),1,LUPRI)
C
         WRITE(LUPRI,*) ' '
C
      END DO ! ISYMB
C
      CALL QEXIT('PRINT_MATAB')
C
      RETURN
      END

C  /* Deck print_matia */
      SUBROUTINE PRINT_MATIA(MATIA,ISYMIA)
**************************************************************
*                                                            *
*     The printing routine, mainly for debugging purposes.   *
*                                                            *
*     Print M(i,a) matrix symmetry-wise                      *
*                                                            *
*     Filip Pawlowski, 04-Dec-2002, Aarhus                   *
*                                                            *
**************************************************************
C
      IMPLICIT NONE 
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

C
      INTEGER ISYMIA,ISYMA,ISYMI
      INTEGER KOFF1
C
#if defined (SYS_CRAY)
      REAL MATIA(*)
#else
      DOUBLE PRECISION MATIA(*)
#endif
C
      CALL QENTER('PRINT_MATIA')
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYMIA,ISYMA)
C
         KOFF1 = IT1AMT(ISYMI,ISYMA) + 1
C
         WRITE(LUPRI,*) 'Symmetry block number(i,a): ',ISYMI,ISYMA
C
         CALL OUTPUT_E(MATIA(KOFF1),1,NRHF(ISYMI),1,NVIR(ISYMA),
     *               NRHF(ISYMI),NVIR(ISYMA),1,LUPRI)
C
         WRITE(LUPRI,*) ' '
C
      END DO ! ISYMA
C
      CALL QEXIT('PRINT_MATIA')
C
      RETURN
      END

C  /* Deck print_matai */
      SUBROUTINE PRINT_MATAI(MATAI,ISYMAI)
**************************************************************
*                                                            *
*     The printing routine, mainly for debugging purposes.   *
*                                                            *
*     Print M(a,i) matrix symmetry-wise                      *
*                                                            *
*     Filip Pawlowski, 04-Dec-2002, Aarhus                   *
*                                                            *
**************************************************************
C
      IMPLICIT NONE 
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

C
      INTEGER ISYMAI,ISYMI,ISYMA
      INTEGER KOFF1
C
#if defined (SYS_CRAY)
      REAL MATAI(*)
#else
      DOUBLE PRECISION MATAI(*)
#endif
C
      CALL QENTER('PRINT_MATAI')
C
      DO ISYMI = 1, NSYM
         ISYMA = MULD2H(ISYMAI,ISYMI)
C
         KOFF1 = IT1AM(ISYMA,ISYMI) + 1
C
         WRITE(LUPRI,*) 'Symmetry block number(a,i): ',ISYMA,ISYMI
C
         CALL OUTPUT_E(MATAI(KOFF1),1,NVIR(ISYMA),1,NRHF(ISYMI),
     *               NVIR(ISYMA),NRHF(ISYMI),1,LUPRI)
C
         WRITE(LUPRI,*) ' '
C
      END DO ! ISYMI
C
      CALL QEXIT('PRINT_MATAI')
C
      RETURN
      END
C  /* Deck output_e */
      SUBROUTINE OUTPUT_E(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,
     *                    COLDIM,NCTL,LUPRI)
C.......................................................................
C Revised 15-Dec-1983 by Hans Jorgen Aa. Jensen.
C         16-Jun-1986 hjaaj ( removed Hollerith )
C
C OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
C AND COLUMNS.  THE INPUT IS AS FOLLOWS;
C
C        AMATRX(',').........MATRIX TO BE OUTPUT
C
C        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
C
C        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
C
C        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
C
C        ROWDIM..............ROW DIMENSION OF AMATRX(',')
C
C        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
C
C        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
C                                                   2 FOR DOUBLE SPACE
C                                                   3 FOR TRIPLE SPACE
C
C THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
C PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
C FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
C CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
C OF COLUMNS.
C
C AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
C          FLORIDA, GAINESVILLE
C REVISED; FEBRUARY 26, 1971
C
C.......................................................................
C
#include <implicit.h>
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
C
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
C
      AMAX = ZERO
      DO 10 J = COLLOW,COLHI
         DO 10 I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
   10 CONTINUE
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
C        use F output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      ELSE
C        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
C
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
C
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.ZERO) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,COLHI)
    3 RETURN
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
C2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
C2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
