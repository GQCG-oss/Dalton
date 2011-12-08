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
C  /* Deck cc3_hmat */
      SUBROUTINE CC3_HMAT(LISTL,IDLSTL,LISTB,IDLSTB,
     *                    LISTC,IDLSTC,
     *                    LISTD,IDLSTD,
     *                    OMEGA1,OMEGA2,
     *                    ISYRES,
     *                    WORK,LWORK)
*
**********************************************************************
*
* Calculate following contributions to hmat:
*
*  <L3|[H^BD,T^{C}_{2}],tau_{1}]|HF>
*  <L3|[H^CD,T^{B}_{2}],tau_{1}]|HF>
*  <L3|[H^BC,T^{D}_{2}],tau_{1}]|HF>
*  <L3|[H^BCD,\tau_nu_2]|HF>
*
* L0 or L1 list may be used for the multipliers; thus hmat or dmat 
* contributions may be calculated, respectively.
*
*---------------------------------------------------------------------
*
* HOWEVER... L1 list would mean D matrix and thus quartic response
*            calculations, which is, in general, NOT implemented
*            in CC3. Therefore, though this routine is general enough
*            to handle L1 list, we are quitting for this list !!!!!!!
*
*---------------------------------------------------------------------
*
* If we have L0 then WB3X = .false.  and  SKIPGEI = .true.
*
* If we have L1 (or other similar lists---like LE,M1,N2---that require
*                W intermediate)  then WB3X = .true.  and  SKIPGEI = .false.
*
*---------------------------------------------------------------------
*
* F. Pawlowski and P. Jorgensen, Spring 2003.
*
* (based on the cc3_bfmat routine)
*
*---------------------------------------------------------------------
*
* April-2004, Aarhus, FP: VVVV integrals removed, flags LVVVV and SKIPGEI 
*                         introduced.
*
**********************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccl1rsp.h>
#include <ccsdsym.h>
#include <ccr1rsp.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccinftap.h>
#include <ccsdinp.h>
C
      LOGICAL SKIPGEI
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
C
      CHARACTER LISTL0*3, LISTL*3,LISTB*3,LISTC*3,LISTD*3,LISTL1R*3,
     *          LABELL1*8
C
      CHARACTER*(8) FN3VI2,FNTOC
      CHARACTER*(7) FN3FOP2X
      CHARACTER*(6) FN3FOP2,FN3VI,FN3FOPX
      CHARACTER*(5) FN3FOP,FNDKBC3
C
      PARAMETER (FN3VI2 = 'CC3_VI12', FN3FOP='PTFOP',  FN3FOP2='PTFOP2')
      PARAMETER (FNTOC  = 'CCSDT_OC', FN3VI ='CC3_VI', FNDKBC3='DKBC3')
      PARAMETER (FN3FOPX= 'PTFOPX'  , FN3FOP2X='PTFOP2X')
C
      CHARACTER*13 FNDELDR,FNDKBCR,FNDKBCR4,FNCKJDR,FN4VBD,FN4VCD,
     *             FN4VBC,FN3SRTR
C
      CHARACTER*7 FNGEIB,FNFEIB
      CHARACTER*6 FNN1B
      PARAMETER(FNGEIB='N1_GEIB' , FNFEIB='N1_FEIB' , FNN1B='N1MATB' )
      INTEGER LUGEIB,LUFEIB,LUN1B
C
      CHARACTER*7 FNGEIC,FNFEIC
      CHARACTER*6 FNN1C
      PARAMETER(FNGEIC='N1_GEIC' , FNFEIC='N1_FEIC' , FNN1C='N1MATC' )
      INTEGER LUGEIC,LUFEIC,LUN1C
C
      CHARACTER*7 FNGEID,FNFEID
      CHARACTER*6 FNN1D
      PARAMETER(FNGEID='N1_GEID' , FNFEID='N1_FEID' , FNN1D='N1MATD' )
      INTEGER LUGEID,LUFEID,LUN1D
C
      INTEGER IDLSTL0,IDLSTL,IDLSTB,IDLSTC,IDLSTD,IDLSTL1R
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER LU3VI2,LU3FOP,LU3FOP2
      INTEGER LUDELDR,LUDKBCR,LUDKBCR4,LUCKJDR,LU4VBD,LU4VCD,LU4VBC
      INTEGER LU3SRTR
C
      CHARACTER CDUMMY*1
C
      LOGICAL   LOCDBG,LORXL1
      PARAMETER (LOCDBG = .FALSE.)
C
      INTEGER ISYRES,LWORK
      INTEGER ISYML1,ISYML1R,ISYMRB
      INTEGER ISINT2,KRBJIA,KLAMP0,KLAMH0,KFOCKD,KFOCK0CK,KT2TP,KL1AM
      INTEGER KL2TP,KEND1,LWRK1
      INTEGER KL1L1,KL2L1,KFOCKL1
      INTEGER KT2RB,KT2RC,KFOCK0
      INTEGER IOPT
      INTEGER ISINT1RB,ISINT2RB,ISINT1BC,ISINT1RBD,ISINT1BCD
      INTEGER ISINT2RBD
      INTEGER KBDIOOOO,KBDIOVVO,KBDIOOVV,KBCIOOOO,KBCIOVVO,KBCIOOVV
      INTEGER ISINT2L1R,ISYFCKL1R,KXIAJB,KT3BOG1,KT3BOL1,KT3BOG2
      INTEGER KT3BOL2,KLAMPL1R,KLAMHL1R
      INTEGER KW3XOGX1,KFOCKL1RCK,KW3BXOGX1,KW3BXOLX1
      INTEGER KW3BXOG1,KW3BXOL1,KT1L1R
      INTEGER LENGTH
      INTEGER ISINT1,ISINT1L1R
      INTEGER IOFF
      INTEGER ISYMD,ISYCKBD0,ISYCKBDL1R,ISYCKBDRB,ISYCKBDRBC
      INTEGER ISYCKBDRBD
      INTEGER KVVVVBD,KVVVVBC
      INTEGER KT3BVDL1,KT3BVDL2,KT3BVDL3,KEND3,LWRK3
      INTEGER KT3BVDG1,KT3BVDG2,KT3BVDG3
      INTEGER KW3BXVDG1,KW3BXVDG2,KW3BXVDL1,KW3BXVDL2,KW3BXVDGX1
      INTEGER KW3BXVDGX2,KW3BXVDLX1,KW3BXVDLX2
      INTEGER KTRVIR,KTRVIR1,KEND4,LWRK4
      INTEGER ISYMB,ISYALJB0,ISYALJD0,ISYALJBL1,ISYALJDL1,ISYMBD
      INTEGER ISCKIJ,ISWBMAT,ISYCKD
      INTEGER KSMAT2,KUMAT2,KDIAG,KINDSQ
      INTEGER KDIAGWB,KINDSQWB
      INTEGER KINDEX,KINDEX2
      INTEGER KINDEXBL1,KINDEXDL1,KTMAT,KW3BMAT
      INTEGER KT3BVBG1,KT3BVBG2,KT3BVBG3,KSMAT4,KUMAT4,KT3BVBL1
      INTEGER KT3BVBL2,KT3BVBL3,KEND5,LWRK5
      INTEGER KINTOC,KTROCR,KTROCR1 

      INTEGER LENSQ,LENSQWB
      INTEGER ISTBD,KTB,LENSQTB,KINDSQTB,ISTB

      INTEGER ISYMRC,ISYMBC,ISYMRBD,ISYMCD,ISINT2RC,KCDIOOOO,KCDIOVVO,
     *        KCDIOOVV
      INTEGER ISYCKBDRC,KVVVVCD,ISINT1RC,ISINT1RCD
      INTEGER ISYCKBDRCD

      INTEGER ISYMRD,ISYMBCD,KT2RD,ISINT2RCD,ISYCKBDRBCD
C
      INTEGER IR1TAMP
      INTEGER ILSTSYM
C
      INTEGER KEND0,LWRK0
C
      INTEGER ISYMN1B,ISYMN2B,ISYMN1C,ISYMN2C,ISYMN1D,ISYMN2D
      INTEGER KN2MATD,KN2MATB,KN2MATC,KINDSQND,KINDSQNB,KINDSQNC
      INTEGER LENSQND,LENSQNB,LENSQNC,ISGEID,ISFEID,ISGEIB,ISFEIB
      INTEGER ISGEIC,ISFEIC,LENGTHGEID,LENGTHGEIB,LENGTHGEIC
      INTEGER KGEID,KFEID,KGEIB,KFEIB,KGEIC,KFEIC
      INTEGER IADR,KLAMPB,KLAMHB,KLAMPC,KLAMHC,KLAMPD,KLAMHD
c
      integer kx3am
c
      LOGICAL WB3X
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2(*)
      REAL WORK(LWORK)
      REAL FREQL1,FREQL1R
      REAL DDOT,XNORMVAL
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQL1,FREQL1R
      DOUBLE PRECISION DDOT,XNORMVAL
#endif
C
      CALL QENTER('HMAT')
C
c     write(lupri,*)'BEFORE '
c     write(lupri,*)'omega1 before CC3_HMAT, LISTL  ', LISTL
c     call PRINT_MATAI(OMEGA1,ISYRES)
c     xnormval = ddot(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
c     write(lupri,*)'norm omega2 before CC3_HMAT ', xnormval
c     CALL OUTPAK(OMEGA2,NT1AMX,1,LUPRI)
C
C----------------------------------------------------
C     Initialise character strings and open files
C----------------------------------------------------
C    
      CDUMMY = ' '
C
      SKIPGEI = .FALSE.
C
      LU3VI2   = -1
      LU3FOP   = -1
      LU3FOP2  = -1
      LUTOC    = -1
      LU3VI    = -1
      LUDKBC3  = -1
      LU3FOPX  = -1
      LU3FOP2X = -1


C
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
      CALL WOPEN2(LU3FOP,FN3FOP,64,0)
      CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
      CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
      CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)
C
C-------------------------
C     Open temporary files
C-------------------------
C
      LU4VBD   = -1
      LU4VCD   = -1
      LU4VBC   = -1
      LU3SRTR  = -1
      LUCKJDR  = -1
      LUDELDR  = -1
      LUDKBCR  = -1
      LUDKBCR4 = -1
C     
      FN4VBD     = 'CC3_FMAT_TMP1'
      FN4VCD     = 'CC3_FMAT_TMP2'
      FN4VBC    = 'CC3_FMAT_TMP3'
      FN3SRTR  = 'CC3_FMAT_TMP4'
      FNCKJDR  = 'CC3_FMAT_TMP5'
      FNDELDR  = 'CC3_FMAT_TMP6'
      FNDKBCR  = 'CC3_FMAT_TMP7'
      FNDKBCR4 = 'CC3_FMAT_TMP8'
C
      IF (LVVVV) THEN
         CALL WOPEN2(LU4VBD,FN4VBD,64,0)
         CALL WOPEN2(LU4VCD,FN4VCD,64,0)
         CALL WOPEN2(LU4VBC,FN4VBC,64,0)
      END IF
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
      CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
      CALL WOPEN2(LUDKBCR4,FNDKBCR4,64,0)
C
C------------------------------------------------------------
C     lists handling
C------------------------------------------------------------
C
      !The routine is general enough to handle LISTL = 'L1 '
      !(but it is NOT debugged for this list!!!),
      !so in principle it can also calculate D matrix contribution 
      !to quartic response function (e.g., third hyoerpolarizability).
      !However, the quartic response function has not been implemented
      !in general for CC3, so we should quit (at the current stage
      !of development) for LISTL = 'L1 '
      IF (LISTL(1:3).EQ.'L1 ') THEN
         WRITE(LUPRI,*)'LISTL = ',LISTL(1:3)
         WRITE(LUPRI,*)'This means quartic response calculations in CC3'
         WRITE(LUPRI,*)'It is not implemented yet... So sth. is wrong!'
         CALL QUIT('Illegal left list in CC3_HMAT')
      END IF
C
      LISTL0 = 'L0 '
      IDLSTL0 = 0
C
      IF (LISTL(1:3).EQ.'L1 ') THEN

         ! get symmetry, frequency and integral label from common blocks
         ! defined in ccl1rsp.h
         ISYML1  = ISYLRZ(IDLSTL)
         FREQL1  = FRQLRZ(IDLSTL)
         LABELL1 = LRZLBL(IDLSTL)
         LORXL1  = LORXLRZ(IDLSTL)
c
c        write(lupri,*)'LABELL1 in CC3_HMAT ', LABELL1
c

         IF (LORXL1) CALL QUIT('NO ORBITAL RELAX. IN CC3_HMAT')

        LISTL1R  = 'R1 '
        IDLSTL1R = IR1TAMP(LABELL1,LORXL1,FREQL1,ISYML1)
        ! get symmetry and frequency from common blocks
        ! defined in ccl1rsp.h
        ISYML1R  = ISYLRT(IDLSTL1R)
        FREQL1R  = FRQLRT(IDLSTL1R)
C
        IF (ISYML1 .NE. ISYML1R) THEN
           WRITE(LUPRI,*)'ISYML1: ', ISYML1
           WRITE(LUPRI,*)'ISYML1R: ', ISYML1R
           CALL QUIT('Symmetry mismatch in CC3_HMAT')
        END IF
C
        IF (FREQL1R .NE. FREQL1) THEN
           WRITE(LUPRI,*)'FREQL1R: ', FREQL1R
           WRITE(LUPRI,*)'FREQL1: ', FREQL1
           CALL QUIT('Frequency mismatch in CC3_HMAT')
        END IF
C
      ELSE IF (LISTL(1:3).EQ.'L0 ') THEN
*        CONTINUE
         SKIPGEI = .TRUE.
      ELSE
         CALL QUIT('Unknown left list in CC3_HMAT')
      END IF
C
      IF (.NOT.LVVVV) THEN
         !Open files for N1MAT^B intermediates
         LUGEIB = -1
         LUFEIB = -1
         LUN1B  = -1
         CALL WOPEN2(LUFEIB,FNFEIB,64,0)
         CALL WOPEN2(LUN1B,FNN1B,64,0)
         IF (.NOT.SKIPGEI) THEN
            CALL WOPEN2(LUGEIB,FNGEIB,64,0)
         END IF
         !Open files for N1MAT^C intermediates
         LUGEIC = -1
         LUFEIC = -1
         LUN1C  = -1
         CALL WOPEN2(LUFEIC,FNFEIC,64,0)
         CALL WOPEN2(LUN1C,FNN1C,64,0)
         IF (.NOT.SKIPGEI) THEN
            CALL WOPEN2(LUGEIC,FNGEIC,64,0)
         END IF
         !Open files for N1MAT^D intermediates
         LUGEID = -1
         LUFEID = -1
         LUN1D  = -1
         CALL WOPEN2(LUFEID,FNFEID,64,0)
         CALL WOPEN2(LUN1D,FNN1D,64,0)
         IF (.NOT.SKIPGEI) THEN
            CALL WOPEN2(LUGEID,FNGEID,64,0)
         END IF
      END IF


c     IF (LISTB(1:3).EQ.'R1 ') THEN
c        ! get symmetry, frequency and integral label for right list 
c        ! from common blocks defined in ccr1rsp.h
c       ISYMRB  = ISYLRT(IDLSTB)
c     ELSE
c        CALL QUIT('Unknown right list in CC3_HMAT')
c     END IF
      ISYMRB = ILSTSYM(LISTB,IDLSTB)
      ISYMRC = ILSTSYM(LISTC,IDLSTC)
      ISYMRD = ILSTSYM(LISTD,IDLSTD)
C
      ISYMBC = MULD2H(ISYMRB,ISYMRC)
      ISYMRBD = MULD2H(ISYMRB,ISYMRD)
      ISYMCD = MULD2H(ISYMRC,ISYMRD)
C
      ISYMBCD = MULD2H(ISYMBC,ISYMRD)
C
      ISINT1 = ISYM0
      ISINT2 = ISYM0
C
      IF (.NOT.LVVVV) THEN
        IF ( LISTL(1:3).EQ.'L1 ') THEN
           !Symmetries for N1^B and N2^B intermediates
           ISYMN1B = MULD2H(ISYML1,ISYMRB)
           ISYMN2B = MULD2H(ISYML1,ISYMRB)
           !Symmetries for N1^C and N2^C intermediates
           ISYMN1C = MULD2H(ISYML1,ISYMRC)
           ISYMN2C = MULD2H(ISYML1,ISYMRC)
           !Symmetries for N1^D and N2^D intermediates
           ISYMN1D = MULD2H(ISYML1,ISYMRD)
           ISYMN2D = MULD2H(ISYML1,ISYMRD)
        ELSE IF (LISTL(1:3).EQ.'L0 ') THEN
           !Symmetries for N1^B and N2^B intermediates
           ISYMN1B = MULD2H(ISYM0,ISYMRB)
           ISYMN2B = MULD2H(ISYM0,ISYMRB)
           !Symmetries for N1^C and N2^C intermediates
           ISYMN1C = MULD2H(ISYM0,ISYMRC)
           ISYMN2C = MULD2H(ISYM0,ISYMRC)
           !Symmetries for N1^D and N2^D intermediates
           ISYMN1D = MULD2H(ISYM0,ISYMRD)
           ISYMN2D = MULD2H(ISYM0,ISYMRD)
        ELSE
           CALL QUIT('Unknown left list in CC3_HMAT(x)')
        END IF
      END IF
C
      IF (LVVVV) THEN
         KRBJIA = 1
      ELSE
         KN2MATD = 1
         KN2MATB = KN2MATD + NCKIJ(ISYMN2D)
         KN2MATC = KN2MATB + NCKIJ(ISYMN2B)
         KRBJIA  = KN2MATC + NCKIJ(ISYMN2C)
      END IF
      KLAMP0 = KRBJIA   + NT2SQ(ISYRES)
      KLAMH0  = KLAMP0  + NLAMDT
      KFOCKD  = KLAMH0  + NLAMDT
      KFOCK0CK  = KFOCKD  + NORBTS
      KT2TP   = KFOCK0CK  + NT1AMX
      KL1AM   = KT2TP   + NT2SQ(ISYM0)
      KL2TP   = KL1AM   + NT1AM(ISYM0)
      KEND0   = KL2TP   + NT2SQ(ISYM0)
      LWRK0   = LWORK   - KEND0
C
      IF (LISTL(1:3).EQ.'L1 ') THEN 
         KL1L1   = KEND0
         KL2L1   = KL1L1   + NT1AM(ISYML1)
         KFOCKL1   = KL2L1   + NT2SQ(ISYML1)
         KEND0   = KFOCKL1    + N2BST(ISYML1)
         LWRK0    = LWORK - KEND0
      END IF
C
      KT2RB   = KEND0
      KFOCK0  = KT2RB   + NT2SQ(ISYMRB)
      KEND0  = KFOCK0    + N2BST(ISYM0)
      LWRK0    = LWORK - KEND0
C
      KT2RC   = KEND0
      KEND0  = KT2RC   + NT2SQ(ISYMRC)
      LWRK0    = LWORK - KEND0
C
      KT2RD   = KEND0
      KEND0  = KT2RD   + NT2SQ(ISYMRD)
      LWRK0    = LWORK - KEND0
C
      IF (.NOT.LVVVV) THEN
         KINDSQND = KEND0
         KINDSQNB = KINDSQND + (6*NCKIJ(ISYMN2D) - 1)/IRAT + 1
         KINDSQNC = KINDSQNB + (6*NCKIJ(ISYMN2B) - 1)/IRAT + 1
         KEND0    = KINDSQNC + (6*NCKIJ(ISYMN2C) - 1)/IRAT + 1
         LWRK0    = LWORK    - KEND0
      END IF
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_HMAT')
      END IF
C
      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
C
      IF (.NOT.LVVVV) THEN
         CALL DZERO(WORK(KN2MATD),NCKIJ(ISYMN2D))
         CALL DZERO(WORK(KN2MATB),NCKIJ(ISYMN2B))
         CALL DZERO(WORK(KN2MATC),NCKIJ(ISYMN2C))
C
         !index array for N2^D
         LENSQND = NCKIJ(ISYMN2D)
         CALL CC3_INDSQ(WORK(KINDSQND),LENSQND,ISYMN2D)
         !index array for N2^B
         LENSQNB = NCKIJ(ISYMN2B)
         CALL CC3_INDSQ(WORK(KINDSQNB),LENSQNB,ISYMN2B)
         !index array for N2^C
         LENSQNC = NCKIJ(ISYMN2C)
         CALL CC3_INDSQ(WORK(KINDSQNC),LENSQNC,ISYMN2C)
      END IF
C
C-------------------------------------
C     Read in lamdap and lamdh
C-------------------------------------
C
      CALL GET_LAMBDA0(WORK(KLAMP0),WORK(KLAMH0),WORK(KEND0),LWRK0)
C
C---------------------------------------------------------------------
C     Read zeroth-order AO Fock matrix from file and trasform it to 
C     lambda basis
C---------------------------------------------------------------------
C
      CALL GET_FOCK0(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),WORK(KEND0),
     *               LWRK0)
C
C---------------------------------------------------------------------
C     Read the matrix the property integrals and trasform it to lambda 
C     basis for L1 list 
C---------------------------------------------------------------------
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         CALL GET_FOCKX(WORK(KFOCKL1),LABELL1,IDLSTL,ISYML1,
     *                  WORK(KLAMP0),ISYM0,
     *                  WORK(KLAMH0),ISYM0,WORK(KEND0),LWRK0)
      END IF
C
C-------------------------------------
C     Read T2 zeroth-order amplitudes 
C-------------------------------------
C
      IOPT = 2
      CALL GET_T1_T2(IOPT,.FALSE.,DUMMY,WORK(KT2TP),'R0',0,ISYM0,
     *                WORK(KEND0),LWRK0)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
     *    DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)
C
C--------------------------------------------
C     Read L1 and L2 zeroth-order multipliers 
C--------------------------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1AM),WORK(KL2TP),LISTL0,
     *               IDLSTL0,ISYM0,WORK(KEND0),LWRK0)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
     *    DDOT(NT2SQ(ISYM0),WORK(KL2TP),1,WORK(KL2TP),1)
C
C---------------------------------------------
C     Read L1L1 and L2L1 multipliers (L1 list)
C---------------------------------------------
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         IOPT  = 3
         CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1L1),WORK(KL2L1),LISTL,
     *                  IDLSTL,ISYML1,WORK(KEND0),LWRK0)
      END IF
C
C-------------------------------------
C     Read T1R1 and L2R1 amplitudes 
C-------------------------------------
C
      IOPT  = 2
      CALL GET_T1_T2(IOPT,.TRUE.,DUMMY,WORK(KT2RB),LISTB,
     *               IDLSTB,ISYMRB,WORK(KEND0),LWRK0)
C
C-------------------------------------
C     Read T1R1 and L2R1 amplitudes 
C-------------------------------------
C
      IOPT  = 2
      CALL GET_T1_T2(IOPT,.TRUE.,DUMMY,WORK(KT2RC),LISTC,
     *               IDLSTC,ISYMRC,WORK(KEND0),LWRK0)
C
C-------------------------------------
C     Read T1R1 and L2R1 amplitudes 
C-------------------------------------
C
      IOPT  = 2
      CALL GET_T1_T2(IOPT,.TRUE.,DUMMY,WORK(KT2RD),LISTD,
     *               IDLSTD,ISYMRD,WORK(KEND0),LWRK0)
C
C----------------------------------------
C     Integrals [H,T1Y] where Y is LISTB
C----------------------------------------
C
      ISINT1RB = MULD2H(ISINT1,ISYMRB)
      ISINT1RC = MULD2H(ISINT1,ISYMRC)
      ISINT2RB = MULD2H(ISINT2,ISYMRB)
      ISINT2RC = MULD2H(ISINT2,ISYMRC)
C
      ISINT1RBD = MULD2H(ISINT1,ISYMRBD)
      ISINT2RBD = MULD2H(ISINT2,ISYMRBD)
      ISINT1RCD = MULD2H(ISINT1,ISYMCD)
      ISINT2RCD = MULD2H(ISINT2,ISYMCD)
C
      ISINT1BC = MULD2H(ISINT1,ISYMBC)
      ISINT1BCD = MULD2H(ISINT1,ISYMBCD)
C
      CALL CC3_3BARINT(ISYMRB,LISTB,IDLSTB,ISYMRC,LISTC,IDLSTC,
     *                 ISYMRD,LISTD,IDLSTD,.TRUE.,
     *                 WORK(KLAMP0),WORK(KLAMH0),WORK(KEND0),LWRK0,
     *                 LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISINT1BCD,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMH0),WORK(KEND0),LWRK0,ISINT1BCD,
     *              LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
      !with option IOPT = 2 it just calculates LUDKBCR4 needed 
      !in contractions
      IOPT = 2
      CALL CC3_TCME(DUMMY,ISINT1BCD,WORK(KEND0),LWRK0,IDUMMY,CDUMMY,
     *              LUDKBCR,FNDKBCR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *              IDUMMY,CDUMMY,LUDKBCR4,FNDKBCR4,IOPT)
C
C---------------------------------------------------------------
C     Read canonical orbital energies and delete frozen orbitals 
C     in Fock diagonal, if required
C---------------------------------------------------------------
C
      CALL GET_ORBEN(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
C--------------------------------------------
C     Sort the Fock matrix to get F(ck) block
C--------------------------------------------
C
      CALL SORT_FOCKCK(WORK(KFOCK0CK),WORK(KFOCK0),ISYM0)
C
C----------------------------------------
C     If we want to sum the T3 amplitudes
C----------------------------------------
C
      if (.false.) then
         kx3am  = kend0
         kend0 = kx3am + nrhft*nrhft*nrhft*nvirt*nvirt*nvirt
         call dzero(work(kx3am),nrhft*nrhft*nrhft*nvirt*nvirt*nvirt)
         lwrk0 = lwork - kend0
         if (lwrk0 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend0
            call quit('Insufficient space (T3) in CC3_HMAT')
         END IF
      endif
C
C      write(lupri,*) 'WBMAT after dzero'
C      call print_pt3(work(kx3am),ISYML1,4)
C
C--------------------------------------------------------------
C     Calculate the normal g^0 integrals for
C     OOOO, OVVO and OOVV integrals used in contraction
C     (VVVV is stored on file LU4VBD and read in in D-loop)
C--------------------------------------------------------------
C
      KBDIOOOO = KEND0
      KBDIOVVO = KBDIOOOO + N3ORHF(ISINT1RBD)
      KBDIOOVV = KBDIOVVO + NT2SQ(ISINT1RBD)
      KEND0   = KBDIOOVV + NT2SQ(ISINT1RBD)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_HMAT (g^0[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KBDIOVVO),NT2SQ(ISINT1RBD))
      CALL DZERO(WORK(KBDIOOVV),NT2SQ(ISINT1RBD))
C
      CALL CC3_INT(WORK(KBDIOOOO),WORK(KBDIOOVV),WORK(KBDIOVVO),
     *             ISINT1RBD,LU4VBD,FN4VBD,
     *             .TRUE.,LISTB,IDLSTB,ISYMRB,
     *             .TRUE.,LISTD,IDLSTD,ISYMRD,
     *             WORK(KEND0),LWRK0)
C
C
C--------------------------------------------------------------
C     Calculate the normal g^0 integrals for
C     OOOO, OVVO and OOVV integrals used in contraction
C     (VVVV is stored on file LU4VCD and read in in D-loop)
C--------------------------------------------------------------
C
      KCDIOOOO = KEND0
      KCDIOVVO = KCDIOOOO + N3ORHF(ISINT1RCD)
      KCDIOOVV = KCDIOVVO + NT2SQ(ISINT1RCD)
      KEND0   = KCDIOOVV + NT2SQ(ISINT1RCD)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_HMAT (g^0[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KCDIOVVO),NT2SQ(ISINT1RCD))
      CALL DZERO(WORK(KCDIOOVV),NT2SQ(ISINT1RCD))
C
      CALL CC3_INT(WORK(KCDIOOOO),WORK(KCDIOOVV),WORK(KCDIOVVO),
     *             ISINT1RCD,LU4VCD,FN4VCD,
     *             .TRUE.,LISTC,IDLSTC,ISYMRC,
     *             .TRUE.,LISTD,IDLSTD,ISYMRD,
     *             WORK(KEND0),LWRK0)
C
C---------------------------------------------------------------
C     Calculate the g^B integrals for
C     OOOO^B, OVVO^B and OOVV^B integrals used in contraction
C     (VVVV^B is stored on file LU4VBC and read in in D-loop)
C---------------------------------------------------------------
C
      KBCIOOOO = KEND0
      KBCIOVVO = KBCIOOOO + N3ORHF(ISINT1BC)
      KBCIOOVV = KBCIOVVO + NT2SQ(ISINT1BC)
      KEND0   = KBCIOOVV + NT2SQ(ISINT1BC)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_HMAT (g^B[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KBCIOOOO),N3ORHF(ISINT1BC))
      CALL DZERO(WORK(KBCIOVVO),NT2SQ(ISINT1BC))
      CALL DZERO(WORK(KBCIOOVV),NT2SQ(ISINT1BC))
C
      CALL CC3_INT(WORK(KBCIOOOO),WORK(KBCIOOVV),WORK(KBCIOVVO),
     *             ISINT1BC,LU4VBC,FN4VBC,
     *             .TRUE.,LISTB,IDLSTB,ISYMRB,
     *             .TRUE.,LISTC,IDLSTC,ISYMRC,
     *             WORK(KEND0),LWRK0)
C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         ISINT2L1R = MULD2H(ISYML1R,ISINT2)
         ISYFCKL1R = MULD2H(ISYMOP,ISYML1R)
      END IF
C
      KXIAJB   = KEND0
      KEND0   = KXIAJB  + NT2AM(ISYM0)

      KT3BOG1 = KEND0
      KT3BOL1 = KT3BOG1 + NTRAOC(ISYM0)
      KT3BOG2 = KT3BOL1 + NTRAOC(ISYM0)
      KT3BOL2 = KT3BOG2 + NTRAOC(ISYM0)
      KEND0  = KT3BOL2 + NTRAOC(ISYM0)
      LWRK0   = LWORK   - KEND0
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         KW3XOGX1 = KEND0
         KFOCKL1RCK   = KW3XOGX1 + NTRAOC(ISINT2L1R)
         KW3BXOGX1    = KFOCKL1RCK    + NT1AM(ISYFCKL1R)
         KW3BXOLX1   = KW3BXOGX1   + NTRAOC(ISINT2L1R)
         KEND0      = KW3BXOLX1   + NTRAOC(ISINT2L1R)
         LWRK0      = LWORK      - KEND0
C
         KW3BXOG1   = KEND0
         KW3BXOL1   = KW3BXOG1   + NTRAOC(ISYM0)
         KEND0   = KW3BXOL1   + NTRAOC(ISYM0)
         LWRK0      = LWORK      - KEND0
C
         KLAMPL1R = KEND0
         KLAMHL1R  = KLAMPL1R  + NLAMDT
         KEND0   = KLAMHL1R  + NLAMDT
         LWRK0   = LWORK   - KEND0
C
      END IF
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         KT1L1R  = KEND0
         KEND0  = KT1L1R + NT1AM(ISYML1R)
         LWRK0   = LWORK  - KEND0
      END IF
C
      KINTOC = KEND0
      KTROCR = KINTOC + NTOTOC(ISINT1BCD)
      KTROCR1 = KTROCR + NTRAOC(ISINT1BCD)
      KEND0   = KTROCR1 + NTRAOC(ISINT1BCD)
      LWRK0   = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_HMAT(x)')
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
C        ISINT1L1R - symmetry of integrals in standard H, transformed
C                  with LambdaH_L1R

      ISINT1  = 1
C
C--------------------------
C     Get Lambda for right list depended on left LISTL list
C--------------------------
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
C
         ISINT1L1R = MULD2H(ISINT1,ISYML1R)
C
         CALL GET_LAMBDAX(WORK(KLAMPL1R),WORK(KLAMHL1R),LISTL1R,
     *                    IDLSTL1R,
     *                    ISYML1R,
     *                    WORK(KLAMP0),WORK(KLAMH0),WORK(KEND0),LWRK0)

C
C------------------------------------------------------------------
C        Calculate the F^L1R matrix (kc elements evaluated and stored 
C        as ck)
C------------------------------------------------------------------
C
         IOPT = 1
         CALL GET_T1_T2(IOPT,.FALSE.,WORK(KT1L1R),DUMMY,LISTL1R,
     *                  IDLSTL1R,ISYML1R,WORK(KEND0),LWRK0)
         CALL CC3LR_MFOCK(WORK(KFOCKL1RCK),WORK(KT1L1R),WORK(KXIAJB),
     *                    ISYFCKL1R)
C
      END IF
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_0 multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BAR0(LUTOC,FNTOC,WORK(KLAMH0),ISYM0,WORK(KT3BOG1),
     *                   WORK(KT3BOL1),WORK(KT3BOG2),WORK(KT3BOL2),
     *                   WORK(KEND0),LWRK0)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_Y multipliers                                             
C-----------------------------------------------------------------
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         CALL INTOCC_T3BARX(.FALSE.,
     *                      LUTOC,FNTOC,ISYMOP,WORK(KLAMH0),ISYM0,
     *                      ISINT1,
     *                      WORK(KLAMHL1R),ISYML1R,ISINT1L1R,
     *                      WORK(KW3BXOG1),
     *                      WORK(KW3BXOL1),WORK(KW3BXOGX1),
     *                      WORK(KW3BXOLX1),
     *                      WORK(KEND0),LWRK0)
      END IF
C
C--------------------------------------------------------------------
C     Read in R1-transformed occupied integrals used in contractions, 
C     transform with lambda_P^0 and sort
C--------------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT1BCD) .GT. 0) THEN
         CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,
     *               NTOTOC(ISINT1BCD))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROCR),WORK(KLAMP0),
     *                    WORK(KEND0),LWRK0,ISINT1BCD)
C
      CALL CCFOP_SORT(WORK(KTROCR),WORK(KTROCR1),ISINT1BCD,1)
C
      CALL CC3_LSORT1(WORK(KTROCR),ISINT1BCD,WORK(KEND0),LWRK0,5)
C
      DO ISYMD = 1,NSYM
C
         ISYCKBD0  = MULD2H(ISYMD,ISYM0)
C
         ISYCKBDRB  = MULD2H(ISYMD,ISINT2RB)
         ISYCKBDRC  = MULD2H(ISYMD,ISINT2RC)
C
         ISYCKBDRBD  = MULD2H(ISYMD,ISINT2RBD)
         ISYCKBDRCD  = MULD2H(ISYMD,ISINT2RCD)
         ISYCKBDRBC  = MULD2H(ISYMD,ISINT1BC)
C
         ISYCKBDRBCD  = MULD2H(ISYMD,ISINT1BCD)
C
        IF (.NOT.LVVVV) THEN
            !Symmetry of arrays needed to construct N1MAT^D
            ISGEID  = MULD2H(ISYMN1D,ISYMD)
            ISFEID  = MULD2H(ISYMN1D,ISYMD)
            !Symmetry of arrays needed to construct N1MAT^B
            ISGEIB  = MULD2H(ISYMN1B,ISYMD)
            ISFEIB  = MULD2H(ISYMN1B,ISYMD)
            !Symmetry of arrays needed to construct N1MAT^C
            ISGEIC  = MULD2H(ISYMN1C,ISYMD)
            ISFEIC  = MULD2H(ISYMN1C,ISYMD)
         END IF
C
         IF (LISTL(1:3).EQ.'L1 ') THEN
            ISYCKBDL1R  = MULD2H(ISYMD,ISYML1R)
         END IF
C
         IF (LVVVV) THEN
            KVVVVBD  = KEND0
            KVVVVCD  = KVVVVBD  + NMAABC(ISYCKBDRBD)
            KVVVVBC  = KVVVVCD  + NMAABC(ISYCKBDRCD)
            KEND1    = KVVVVBC  + NMAABC(ISYCKBDRBC)
            LWRK1    = LWORK    - KEND1
         ELSE
            IF (SKIPGEI) THEN !we want to have the dummy allocations
                              !to skip GEI, but keep the code simple.
               LENGTHGEID = 1
               LENGTHGEIB = 1
               LENGTHGEIC = 1
            ELSE
               LENGTHGEID = NCKATR(ISGEID)
               LENGTHGEIB = NCKATR(ISGEIB)
               LENGTHGEIC = NCKATR(ISGEIC)
            END IF
            !Arrays needed to construct N1MAT^D
            KGEID = KEND0
            KFEID = KGEID + LENGTHGEID
            KEND1 = KFEID + NCKATR(ISFEID)
            LWRK1 = LWORK - KEND1
C
            !Arrays needed to construct N1MAT^B
            KGEIB  = KEND1
            KFEIB  = KGEIB  + LENGTHGEIB
            KEND1  = KFEIB  + NCKATR(ISFEIB)
            LWRK1  = LWORK  - KEND1
C
            !Arrays needed to construct N1MAT^C
            KGEIC  = KEND1
            KFEIC  = KGEIC  + LENGTHGEIC
            KEND1  = KFEIC  + NCKATR(ISFEIC)
            LWRK1  = LWORK  - KEND1
C
         END IF

C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space (isymd-loop) in CC3_HMAT')
         END IF
C
         DO D = 1,NVIR(ISYMD)
C
            KT3BVDL1  = KEND1
            KT3BVDL2  = KT3BVDL1 + NCKATR(ISYCKBD0)
            KT3BVDL3  = KT3BVDL2 + NCKATR(ISYCKBD0)
            KEND3   = KT3BVDL3 + NCKATR(ISYCKBD0)
            LWRK3   = LWORK  - KEND3

            KT3BVDG1 = KEND3
            KT3BVDG2 = KT3BVDG1 + NCKATR(ISYCKBD0)
            KT3BVDG3 = KT3BVDG2 + NCKATR(ISYCKBD0)
            KEND3   = KT3BVDG3 + NCKATR(ISYCKBD0)
            LWRK3   = LWORK  - KEND3

            IF (LISTL(1:3).EQ.'L1 ') THEN
               KW3BXVDG1  = KEND3
               KW3BXVDG2  = KW3BXVDG1  + NCKATR(ISYCKBD0)
               KW3BXVDL1  = KW3BXVDG2  + NCKATR(ISYCKBD0)
               KW3BXVDL2  = KW3BXVDL1  + NCKATR(ISYCKBD0)
               KEND3     = KW3BXVDL2  + NCKATR(ISYCKBD0)
               LWRK3     = LWORK     - KEND3
C
               KW3BXVDGX1  = KEND3
               KW3BXVDGX2  = KW3BXVDGX1  + NCKATR(ISYCKBDL1R)
               KW3BXVDLX1  = KW3BXVDGX2  + NCKATR(ISYCKBDL1R)
               KW3BXVDLX2  = KW3BXVDLX1  + NCKATR(ISYCKBDL1R)
               KEND3     = KW3BXVDLX2  + NCKATR(ISYCKBDL1R)
               LWRK3     = LWORK     - KEND3
            END IF
C
            KTRVIR = KEND3
            KTRVIR1 = KTRVIR + NCKATR(ISYCKBDRBCD)
            KEND4 = KTRVIR1 + NCKATR(ISYCKBDRBCD)
            LWRK4  = LWORK  - KEND4

            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space (d-loop) in CC3_HMAT')
            END IF
C
            IF (.NOT.LVVVV) THEN
               CALL DZERO(WORK(KFEID),NCKATR(ISFEID))
               CALL DZERO(WORK(KFEIB),NCKATR(ISFEIB))
               CALL DZERO(WORK(KFEIC),NCKATR(ISFEIC))
               CALL DZERO(WORK(KGEID),LENGTHGEID)
               CALL DZERO(WORK(KGEIB),LENGTHGEIB)
               CALL DZERO(WORK(KGEIC),LENGTHGEIC)
            END IF
C
            IF (LVVVV) THEN
C
C---------------------------------------------------------------------
C              Read in g^0_{vvvv} (used in contraction) for a given D
C---------------------------------------------------------------------
C
               IF (NMAABC(ISYCKBDRBD) .GT. 0) THEN
                  IOFF = I3VVIR(ISYCKBDRBD,ISYMD)
     *                 + NMAABC(ISYCKBDRBD)*(D-1)
     *                 + 1
                  CALL GETWA2(LU4VBD,FN4VBD,WORK(KVVVVBD),IOFF,
     *                        NMAABC(ISYCKBDRBD))
               ENDIF
C
C---------------------------------------------------------------------
C              Read in g^0_{vvvv} (used in contraction) for a given D
C---------------------------------------------------------------------
C
               IF (NMAABC(ISYCKBDRCD) .GT. 0) THEN
                  IOFF = I3VVIR(ISYCKBDRCD,ISYMD)
     *                 + NMAABC(ISYCKBDRCD)*(D-1)
     *                 + 1
                  CALL GETWA2(LU4VCD,FN4VCD,WORK(KVVVVCD),IOFF,
     *                        NMAABC(ISYCKBDRCD))
               ENDIF
C
C---------------------------------------------------------------------
C              Read in g^B_{vvvv} (used in contraction) for a given D
C---------------------------------------------------------------------
C
               IF (NMAABC(ISYCKBDRBC) .GT. 0) THEN
                  IOFF = I3VVIR(ISYCKBDRBC,ISYMD)
     *                 + NMAABC(ISYCKBDRBC)*(D-1)
     *                 + 1
                  CALL GETWA2(LU4VBC,FN4VBC,WORK(KVVVVBC),IOFF,
     *                        NMAABC(ISYCKBDRBC))
               ENDIF
C
            END IF !LVVVV
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
            IF (LISTL(1:3).EQ.'L1 ') THEN
               CALL INTVIR_T3BARX_D(.FALSE.,
     *                           ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                           LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                           WORK(KW3BXVDGX1),WORK(KW3BXVDG1),
     *                           WORK(KW3BXVDGX2),WORK(KW3BXVDG2),
     *                           WORK(KW3BXVDLX1),WORK(KW3BXVDL1),
     *                           WORK(KW3BXVDLX2),WORK(KW3BXVDL2),
     *                           WORK(KLAMPL1R),ISYML1R,WORK(KLAMP0),
     *                           ISYM0,ISYMD,D,WORK(KEND4),LWRK4)
            END IF
C
C---------------------------------------------------------------------
C           Read and sort R1-transformed integrals used in contraction.
C---------------------------------------------------------------------
C
            IF (NCKATR(ISYCKBDRBCD) .GT. 0) THEN
               IOFF = ICKBD(ISYCKBDRBCD,ISYMD)+NCKATR(ISYCKBDRBCD)*(D-1)
     *              + 1
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVIR),IOFF,
     &                     NCKATR(ISYCKBDRBCD))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISYCKBDRBCD)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_HMAT (TRVI)')
            END IF
C
            CALL CCSDT_SRVIR3(WORK(KTRVIR),WORK(KEND4),ISYMD,D,
     *                        ISINT1BCD)
C
C
            IF (NCKATR(ISYCKBDRBCD) .GT. 0) THEN
               IOFF = ICKBD(ISYCKBDRBCD,ISYMD)+NCKATR(ISYCKBDRBCD)*(D-1)
     *              + 1
               CALL GETWA2(LUDKBCR4,FNDKBCR4,WORK(KTRVIR1),IOFF,
     &                     NCKATR(ISYCKBDRBCD))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISYCKBDRBCD)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_HMAT (TRVI1)')
            END IF
C
            CALL CCSDT_SRVIR3(WORK(KTRVIR1),WORK(KEND4),ISYMD,D,
     *                        ISINT1BCD)
C

C
            DO ISYMB = 1,NSYM
C
               ISYALJB0  = MULD2H(ISYMB,ISYM0)
               ISYALJD0 = MULD2H(ISYMD,ISYM0)
               ISYMBD  = MULD2H(ISYMD,ISYMB)
               ISCKIJ  = MULD2H(ISYMBD,ISYM0)
               ISYCKD  = MULD2H(ISYM0,ISYMB)
               IF (LISTL(1:3).EQ.'L1 ') THEN
                  ISYALJBL1  = MULD2H(ISYMB,ISYML1)
                  ISYALJDL1 = MULD2H(ISYMD,ISYML1)
                  ISWBMAT  = MULD2H(ISCKIJ,ISYML1)
               END IF
C
               KSMAT2     = KEND4
               KUMAT2     = KSMAT2    + NCKIJ(ISCKIJ)
               KDIAG      = KUMAT2    + NCKIJ(ISCKIJ)
               KINDSQ     = KDIAG     + NCKIJ(ISCKIJ)
               KEND5      = KINDSQ    + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1

               IF (LISTL(1:3).EQ.'L1 ') THEN
                  KDIAGWB     = KEND5
                  KINDSQWB     = KDIAGWB    + NCKIJ(ISWBMAT)
                  KEND5    = KINDSQWB  + (6*NCKIJ(ISWBMAT) - 1)/IRAT + 1
               END IF

               KINDEX     = KEND5
               KINDEX2    = KINDEX    + (NCKI(ISYALJB0)  - 1)/IRAT + 1
               KEND5      = KINDEX2   + (NCKI(ISYALJD0) - 1)/IRAT + 1
C
               IF (LISTL(1:3).EQ.'L1 ') THEN
                  KINDEXBL1   = KEND5
                  KINDEXDL1  = KINDEXBL1 + (NCKI(ISYALJBL1)-1)/IRAT + 1
                  KTMAT      = KINDEXDL1  + (NCKI(ISYALJDL1)-1)/IRAT + 1
                  KW3BMAT    = KTMAT + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWBMAT))
                  KEND5      = KW3BMAT     + NCKIJ(ISWBMAT)
                  LWRK5      = LWORK     - KEND5
               ELSE
                  KTMAT      = KEND5
                  KEND5      = KTMAT + NCKIJ(ISCKIJ)
               END IF

               KT3BVBG1 = KEND5
               KT3BVBG2 = KT3BVBG1 + NCKATR(ISYCKD)
               KT3BVBG3 = KT3BVBG2 + NCKATR(ISYCKD)
               KEND5   = KT3BVBG3 + NCKATR(ISYCKD)
               LWRK5   = LWORK   - KEND5

               KSMAT4  = KEND5
               KUMAT4  = KSMAT4  + NCKIJ(ISCKIJ)
               KT3BVBL1 = KUMAT4  + NCKIJ(ISCKIJ)
               KT3BVBL2 = KT3BVBL1 + NCKATR(ISYCKD)
               KT3BVBL3 = KT3BVBL2 + NCKATR(ISYCKD)
               KEND5   = KT3BVBL3 + NCKATR(ISYCKD)
               LWRK5   = LWORK   - KEND5
C
               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space(isymb-loop) 
     *                       in CC3_HMAT')
               END IF
C
C--------------------------------------------
C              Construct part of the diagonal
C--------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
               IF (LISTL(1:3).EQ.'L1 ') THEN
                  CALL CC3_DIAG(WORK(KDIAGWB),WORK(KFOCKD),ISWBMAT)
               END IF
C
C------------------------------------
C              Construct index arrays
C------------------------------------
C
               LENSQ  = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               IF (LISTL(1:3).EQ.'L1 ') THEN
                  LENSQWB  = NCKIJ(ISWBMAT)
                  CALL CC3_INDSQ(WORK(KINDSQWB),LENSQWB,ISWBMAT)
               END IF
C
               CALL CC3_INDEX(WORK(KINDEX),ISYALJB0)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJD0)
               IF (LISTL(1:3).EQ.'L1 ') THEN
                  CALL CC3_INDEX(WORK(KINDEXBL1),ISYALJBL1)
                  CALL CC3_INDEX(WORK(KINDEXDL1),ISYALJDL1)
               END IF
C
               DO B = 1,NVIR(ISYMB)
C
                  IF (LISTL(1:3).EQ.'L1 ') THEN
                     CALL DZERO(WORK(KW3BMAT),NCKIJ(ISWBMAT))
                  END IF
C
C-----------------------------------------------------------------------
C                 Construct virtual integrals (for fixed B) which are 
C                 required to calculate t3bar_0 multipliers
C                 (the same routine as in d-loop is used)
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
C---------------------------------------------------------
C                 Get T3bar_BD multipliers (using S and U)
C---------------------------------------------------------
C
                  CALL GET_T3BAR0_BD(ISYM0,WORK(KL1AM),ISYM0,
     *                               WORK(KL2TP),ISYM0,WORK(KTMAT),
     *                               WORK(KFOCK0CK),WORK(KFOCKD),
     *                               WORK(KDIAG),WORK(KXIAJB),ISYM0,
     *                               ISYM0,WORK(KINDSQ),LENSQ,
     *                               WORK(KSMAT2),WORK(KT3BVDG1),
     *                               WORK(KT3BVDG2),WORK(KT3BVDL1),
     *                               WORK(KT3BVDL2),WORK(KT3BOG1),
     *                               WORK(KT3BOL1),WORK(KINDEX),
     *                               WORK(KSMAT4),WORK(KT3BVBG1),
     *                               WORK(KT3BVBG2),WORK(KT3BVBL1),
     *                               WORK(KT3BVBL2),WORK(KINDEX2),
     *                               WORK(KUMAT2),WORK(KT3BVDG3),
     *                               WORK(KT3BVDL3),WORK(KT3BOG2),
     *                               WORK(KT3BOL2),WORK(KUMAT4),
     *                               WORK(KT3BVBG3),WORK(KT3BVBL3),
     *                               ISYMB,B,ISYMD,D,ISCKIJ,
     *                               WORK(KEND5),LWRK5)
c
c       call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),5)
C
C----------------------------------------------------
C                 Get T3barX_BD multipliers (using W)
C----------------------------------------------------
C
                  IF (LISTL(1:3).EQ.'L1 ') THEN
                    CALL GET_T3BARX_BD(.FALSE.,
     *                               WORK(KTMAT),ISCKIJ,WORK(KFOCKL1),
     *                               ISYML1,WORK(KW3BMAT),ISWBMAT,
     *                               WORK(KL2TP),ISYM0,WORK(KFOCKL1RCK),
     *                               ISYFCKL1R,WORK(KW3BXVDLX2),
     *                               WORK(KW3BXVDLX1),WORK(KW3BXVDGX2),
     *                               WORK(KW3BXVDGX1),WORK(KW3BXOLX1),
     *                               WORK(KW3BXOGX1),ISINT2L1R,
     *                               WORK(KINDEX),WORK(KINDEX2),
     *                               WORK(KINDSQWB),LENSQWB,WORK(KL2L1),
     *                               ISYML1,WORK(KFOCK0CK),ISYM0,
     *                               WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                               WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                               WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                               ISINT2,WORK(KINDEXBL1),
     *                               WORK(KINDEXDL1),WORK(KL1L1),ISYML1,
     *                               WORK(KXIAJB),ISINT1,-FREQL1,
     *                               WORK(KDIAGWB),WORK(KFOCKD),B,ISYMB,
     *                               D,ISYMD,ISYML1,WORK(KEND5),LWRK5)
c
c       call sum_pt3(work(KW3BMAT),isymb,b,isymd,d,
c    *             ISWBMAT,work(kx3am),4)
                  END IF
C
C
C-------------------------------------------------------------------
C                 Set up variables depending on the LISTL (L1 or L0)
C-------------------------------------------------------------------
C
C If WB3X = .TRUE. (i.e, LISTL = L1), then KTB contains wbar_X
C else
C we get tbar_0 in KTB 
C
C Note Wbar_X and Tbar_0 is stored in the same way for fixed B and D
C
C Inside the routines CC3_W3_OMEGA1, CC3_W3_CY2V, CC3_W3_CY2O terms 
C are selected depending on WB3X

                  IF (LISTL(1:3).EQ.'L1 ') THEN
                     WB3X = .TRUE.
                     ISTBD = ISWBMAT
                     KTB   = KW3BMAT
                     LENSQTB = LENSQWB
                     KINDSQTB = KINDSQWB
                  ELSE
                     WB3X = .FALSE.
                     ISTBD = ISCKIJ
                     !KTMAT will be destroyed so an extra array for
                     !zeroth-order multipliers is needed
                     KTB   = KEND5
                     KEND5 = KTB + NCKIJ(ISCKIJ)
                     LWRK5 = LWORK - KEND5
                     IF (LWRK5 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Memory available : ',LWORK
                        WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                        CALL QUIT('Too little space(setup) 
     *                              in CC3_HMAT')
                     END IF
                     CALL DCOPY(NCKIJ(ISCKIJ),WORK(KTMAT),1,WORK(KTB),1)

                     LENSQTB = LENSQ
                     KINDSQTB = KINDSQ
                  END IF
                  ISTB = MULD2H(ISTBD,ISYMBD)
C
                  !Check the consistency of logical flags
                  IF (WB3X.AND.(.NOT.SKIPGEI)) THEN
                    CONTINUE
                  ELSE IF ((.NOT.WB3X).AND.SKIPGEI) THEN
                    CONTINUE
                  ELSE
                    WRITE(LUPRI,*)'WB3X    = ',WB3X
                    WRITE(LUPRI,*)'SKIPGEI = ',SKIPGEI
                    WRITE(LUPRI,*)'WB3X and SKIPGEI should be opposite'
                    CALL QUIT('Inconsistent flags in CC3_GMATNEW')
                  END IF
C
C-------------------------------------------------
C                 Calculate <L3|[H^BCD,\tau_nu_2]|HF>
C-------------------------------------------------
C
                  CALL CC3_W3_CY2V(OMEGA2,ISYRES,WORK(KRBJIA),
     *                             WORK(KTB),ISTBD,WORK(KTMAT),
     *                             WORK(KTRVIR),WORK(KTRVIR1),ISINT1BCD,
     *                             WORK(KEND5),LWRK5,WORK(KINDSQTB),
     *                             LENSQTB,ISYMB,B,ISYMD,D,WB3X)
C
                  CALL CC3_W3_CY2O(OMEGA2,ISYRES,WORK(KTB),
     *                             ISTBD,
     *                             WORK(KTMAT),WORK(KTROCR),
     *                             WORK(KTROCR1),ISINT1BCD,
     *                             WORK(KEND5),LWRK5,WORK(KINDSQTB),
     *                             LENSQTB,ISYMB,B,ISYMD,D,WB3X)
C
                  IF (LVVVV) THEN
C
                     !Calculate <L3|[H^BD,T^{C}_{2}],tau_{1}]|HF>
                     CALL CC3_W3_OMEGA1(OMEGA1,ISYRES,WORK(KTB),
     *                                  WORK(KTMAT),ISTB,
     *                                  WORK(KBDIOOOO),WORK(KBDIOVVO),
     *                                  WORK(KBDIOOVV),WORK(KVVVVBD),
     *                                  ISINT1RBD,WORK(KT2RC),ISYMRC,
     *                                  WORK(KEND5),LWRK5,
     *                                  LENSQTB,WORK(KINDSQTB),
     *                                  ISYMB,B,ISYMD,D,WB3X)
C
                     !Calculate <L3|[H^CD,T^{B}_{2}],tau_{1}]|HF>
                     CALL CC3_W3_OMEGA1(OMEGA1,ISYRES,WORK(KTB),
     *                                  WORK(KTMAT),ISTB,
     *                                  WORK(KCDIOOOO),WORK(KCDIOVVO),
     *                                  WORK(KCDIOOVV),WORK(KVVVVCD),
     *                                  ISINT1RCD,WORK(KT2RB),ISYMRB,
     *                                  WORK(KEND5),LWRK5,
     *                                  LENSQTB,WORK(KINDSQTB),
     *                                  ISYMB,B,ISYMD,D,WB3X)
C
                     !Calculate <L3|[H^BC,T^{D}_{2}],tau_{1}]|HF>
                     CALL CC3_W3_OMEGA1(OMEGA1,ISYRES,WORK(KTB),
     *                                  WORK(KTMAT),ISTB,
     *                                  WORK(KBCIOOOO),WORK(KBCIOVVO),
     *                                  WORK(KBCIOOVV),WORK(KVVVVBC),
     *                                  ISINT1BC,WORK(KT2RD),ISYMRD,
     *                                  WORK(KEND5),LWRK5,
     *                                  LENSQTB,WORK(KINDSQTB),
     *                                  ISYMB,B,ISYMD,D,WB3X)
C
                  ELSE
C
                     CALL DSCAL(NCKIJ(ISTBD),-1.0D0,WORK(KTB),1)

                     !<L3|[H^BD,T^{C}_{2}],tau_{1}]|HF> in terms of N1 and N2 
                     CALL WT2_N1N2(WORK(KTB),ISTB,
     *                       WORK(KT2RC),ISYMRC,
     *                       WORK(KGEIC),WORK(KFEIC),
     *                       ISYMN1C,
     *                       WORK(KN2MATC),ISYMN2C,
     *                       B,ISYMB,D,ISYMD,
     *                       WORK(KINDSQTB),LENSQTB,
     *                       WORK(KINDSQNC),LENSQNC,
     *                       WORK(KEND5),LWRK5,
     *                       WB3X)
C
                     !<L3|[H^CD,T^{B}_{2}],tau_{1}]|HF>in terms of N1 and N2 
                     CALL WT2_N1N2(WORK(KTB),ISTB,
     *                       WORK(KT2RB),ISYMRB,
     *                       WORK(KGEIB),WORK(KFEIB),
     *                       ISYMN1B,
     *                       WORK(KN2MATB),ISYMN2B,
     *                       B,ISYMB,D,ISYMD,
     *                       WORK(KINDSQTB),LENSQTB,
     *                       WORK(KINDSQNB),LENSQNB,
     *                       WORK(KEND5),LWRK5,
     *                       WB3X)
C
                     !<L3|[H^BC,T^{D}_{2}],tau_{1}]|HF>in terms of N1 and N2 
                     CALL WT2_N1N2(WORK(KTB),ISTB,
     *                       WORK(KT2RD),ISYMRD,
     *                       WORK(KGEID),WORK(KFEID),
     *                       ISYMN1D,
     *                       WORK(KN2MATD),ISYMN2D,
     *                       B,ISYMB,D,ISYMD,
     *                       WORK(KINDSQTB),LENSQTB,
     *                       WORK(KINDSQND),LENSQND,
     *                       WORK(KEND5),LWRK5,
     *                       WB3X)

                   END IF

               ENDDO   ! B
            ENDDO      ! ISYMB
C
            IF (.NOT.LVVVV) THEN
C
C            ----------------------------------------------------------
C            Put KGEI(ge,i)^F and KFEI(fe,i)^G (which are intermediates
C            for N1MAT(fge,i) ) to files (for fixed F=D and G=D).
C            ----------------------------------------------------------

             IF (.NOT.SKIPGEI) THEN
              !Put KGEID to file as (gei,F)   (fixed F corresponds to D)
              IADR = ICKBD(ISGEID,ISYMD) + NCKATR(ISGEID)*(D-1) + 1
              CALL PUTWA2(LUGEID,FNGEID,WORK(KGEID),IADR,NCKATR(ISGEID))
             END IF
C
             !Put KFEID to file as (fei,G)   (fixed G corresponds to D)
             IADR = ICKBD(ISFEID,ISYMD) + NCKATR(ISFEID)*(D-1) + 1
             CALL PUTWA2(LUFEID,FNFEID,WORK(KFEID),IADR,NCKATR(ISFEID))
C
             IF (.NOT.SKIPGEI) THEN
              !Put KGEIB to file as (gei,F)   (fixed F corresponds to D)
              IADR = ICKBD(ISGEIB,ISYMD) + NCKATR(ISGEIB)*(D-1) + 1
              CALL PUTWA2(LUGEIB,FNGEIB,WORK(KGEIB),IADR,NCKATR(ISGEIB))
             END IF
C
             !Put KFEIB to file as (fei,G)   (fixed G corresponds to D)
             IADR = ICKBD(ISFEIB,ISYMD) + NCKATR(ISFEIB)*(D-1) + 1
             CALL PUTWA2(LUFEIB,FNFEIB,WORK(KFEIB),IADR,NCKATR(ISFEIB))
C
             IF (.NOT.SKIPGEI) THEN
              !Put KGEIC to file as (gei,F)   (fixed F corresponds to D)
              IADR = ICKBD(ISGEIC,ISYMD) + NCKATR(ISGEIC)*(D-1) + 1
              CALL PUTWA2(LUGEIC,FNGEIC,WORK(KGEIC),IADR,NCKATR(ISGEIC))
             END IF
C
             !Put KFEIC to file as (fei,G)   (fixed G corresponds to D)
             IADR = ICKBD(ISFEIC,ISYMD) + NCKATR(ISFEIC)*(D-1) + 1
             CALL PUTWA2(LUFEIC,FNFEIC,WORK(KFEIC),IADR,NCKATR(ISFEIC))
C
            END IF
C
         ENDDO       ! D
      ENDDO          ! ISYMD 

C
C------------------------------------------------------
C     Accumulate RBJIA from the virtual contribution
C     in OMEGA2
C------------------------------------------------------
C
      CALL CC3_RBJIA(OMEGA2,ISYRES,WORK(KRBJIA))
C
      IF (.NOT.LVVVV) THEN
C
         KLAMPB = KEND0
         KLAMHB = KLAMPB + NGLMDT(ISYMRB)
         KLAMPC = KLAMHB + NGLMDT(ISYMRB)
         KLAMHC = KLAMPC + NGLMDT(ISYMRC)
         KLAMPD = KLAMHC + NGLMDT(ISYMRC)
         KLAMHD = KLAMPD + NGLMDT(ISYMRD)
         KEND0  = KLAMHD + NGLMDT(ISYMRD)
         LWRK0  = LWORK  - KEND0
         IF (LWRK0 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND0
            CALL QUIT('Insufficient space in CC3_GMATNEW(final)')
         END IF
C
         !Lambda^B (particle) needed in backtransformation of N1^C
         CALL GET_LAMBDAX(WORK(KLAMPB),WORK(KLAMHB),LISTB,IDLSTB,ISYMRB,
     *                    WORK(KLAMP0),WORK(KLAMH0),
     *                    WORK(KEND0),LWRK0)
C
         !Lambda^C (particle) needed in backtransformation of N1^C
         CALL GET_LAMBDAX(WORK(KLAMPC),WORK(KLAMHC),LISTC,IDLSTC,ISYMRC,
     *                    WORK(KLAMP0),WORK(KLAMH0),
     *                    WORK(KEND0),LWRK0)
C
         !Lambda^D (particle) needed in backtransformation of N1^D
         CALL GET_LAMBDAX(WORK(KLAMPD),WORK(KLAMHD),LISTD,IDLSTD,ISYMRD,
     *                    WORK(KLAMP0),WORK(KLAMH0),
     *                    WORK(KEND0),LWRK0)
C
         !Read (gei,F) and (fei,G) intermediates from files
         !add them and put the result to a file as (fge,I)
         CALL N1_RESORT(ISYMN1C,LUN1C,FNN1C,LUGEIC,FNGEIC,LUFEIC,FNFEIC,
     *                  WORK(KEND0),LWRK0,SKIPGEI)
C
         !Calculate <T3|[[H^BD,T2^C],tau_ai]|HF> except VVVV contribution
         CALL N1N2_G(LUN1C,FNN1C,
     *                     ISYMN1C,
     *                     WORK(KN2MATC),ISYMN2C,
     *                     WORK(KBDIOVVO),WORK(KBDIOOVV),
     *                     WORK(KBDIOOOO),ISINT1RBD,
     *                     OMEGA1,ISYRES,
     *                     WORK(KINDSQNC),LENSQNC,
     *                     WORK(KEND0),LWRK0)
C
         !Calculate VVVV contribution to <T3|[[H^BD,T2^C],tau_ai]|HF>
         IOPT = 1 !(double) T1 one-index backtransformation
         CALL  N1_GV4(IOPT,
     *                LUN1C,FNN1C,
     *                ISYMN1C,
     *                WORK(KLAMPB),ISYMRB,
     *                WORK(KLAMPD),ISYMRD,
     *                WORK(KLAMH0),1,
     *                WORK(KLAMH0),1,
     *                OMEGA1,ISYRES,
     *                WORK(KEND0),LWRK0)
C
C         REPEAT THINGS FOR N1MATB
C
         !Read (gei,F) and (fei,G) intermediates from files
         !add them and put the result to a file as (fge,I)
         CALL N1_RESORT(ISYMN1B,LUN1B,FNN1B,LUGEIB,FNGEIB,LUFEIB,FNFEIB,
     *                  WORK(KEND0),LWRK0,SKIPGEI)
C
         !Calculate <T3|[[H^CD,T2^B],tau_ai]|HF> except VVVV contribution
         CALL N1N2_G(LUN1B,FNN1B,
     *                     ISYMN1B,
     *                     WORK(KN2MATB),ISYMN2B,
     *                     WORK(KCDIOVVO),WORK(KCDIOOVV),
     *                     WORK(KCDIOOOO),ISINT1RCD,
     *                     OMEGA1,ISYRES,
     *                     WORK(KINDSQNB),LENSQNB,
     *                     WORK(KEND0),LWRK0)
C
         !Calculate VVVV contribution to <T3|[[H^CD,T2^B],tau_ai]|HF>
         IOPT = 1 !(double) T1 one-index backtransformation
         CALL  N1_GV4(IOPT,
     *                LUN1B,FNN1B,
     *                ISYMN1B,
     *                WORK(KLAMPC),ISYMRC,
     *                WORK(KLAMPD),ISYMRD,
     *                WORK(KLAMH0),1,
     *                WORK(KLAMH0),1,
     *                OMEGA1,ISYRES,
     *                WORK(KEND0),LWRK0)
C
C         REPEAT THINGS FOR N1MAT^D
C
C
         !Read (gei,F) and (fei,G) intermediates from files
         !add them and put the result to a file as (fge,I)
         CALL N1_RESORT(ISYMN1D,LUN1D,FNN1D,LUGEID,FNGEID,LUFEID,FNFEID,
     *                  WORK(KEND0),LWRK0,SKIPGEI)
C
         !Calculate <T3|[[H^BC,T2^D],tau_ai]|HF> except VVVV contribution
         CALL N1N2_G(LUN1D,FNN1D,
     *                     ISYMN1D,
     *                     WORK(KN2MATD),ISYMN2D,
     *                     WORK(KBCIOVVO),WORK(KBCIOOVV),
     *                     WORK(KBCIOOOO),ISINT1BC,
     *                     OMEGA1,ISYRES,
     *                     WORK(KINDSQND),LENSQND,
     *                     WORK(KEND0),LWRK0)
C
         !Calculate VVVV contribution to <T3|[[H^BC,T2^D],tau_ai]|HF>
         IOPT = 1 !(double) T1 one-index backtransformation
         CALL  N1_GV4(IOPT,
     *                LUN1D,FNN1D,
     *                ISYMN1D,
     *                WORK(KLAMPB),ISYMRB,
     *                WORK(KLAMPC),ISYMRC,
     *                WORK(KLAMH0),1,
     *                WORK(KLAMH0),1,
     *                OMEGA1,ISYRES,
     *                WORK(KEND0),LWRK0)
C
      END IF
C

C
C
c     write(lupri,*)'AFTER '
c     write(lupri,*)'omega1 after CC3_HMAT,LISTL  ', LISTL
c     call PRINT_MATAI(OMEGA1,ISYRES)
*     xnormval = ddot(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
*     write(lupri,*)'norm omega1 after CC3_HMAT ', xnormval
*     xnormval = ddot(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
*     write(lupri,*)'norm omega2 after CC3_HMAT ', xnormval
c     CALL OUTPAK(OMEGA2,NT1AMX,1,LUPRI)
COMMENT COMMENT
c      write(lupri,*) 't3barx  in CC3_HMAT'
c      call print_pt3(work(kx3am),ISYM0,4)
COMMENT COMMENT
C
C
C-----------------
C     Close files
C-----------------
C
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
      CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
      CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
      CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')
C
C-------------------------------------
C     Close and delete temporary files
C-------------------------------------
C
      IF (LVVVV) THEN
         CALL WCLOSE2(LU4VBD,FN4VBD,'DELETE')
         CALL WCLOSE2(LU4VCD,FN4VCD,'DELETE')
         CALL WCLOSE2(LU4VBC,FN4VBC,'DELETE')
      END IF
C
      IF (.NOT.LVVVV) THEN
         !Close files for N1MATB intermediates
         CALL WCLOSE2(LUFEIB,FNFEIB,'DELETE')
         CALL WCLOSE2(LUN1B,FNN1B,'DELETE')
         IF (.NOT.SKIPGEI) THEN
            CALL WCLOSE2(LUGEIB,FNGEIB,'DELETE')
         END IF
         !Close files for N1MATC intermediates
         CALL WCLOSE2(LUFEIC,FNFEIC,'DELETE')
         CALL WCLOSE2(LUN1C,FNN1C,'DELETE')
         IF (.NOT.SKIPGEI) THEN
            CALL WCLOSE2(LUGEIC,FNGEIC,'DELETE')
         END IF
         !Close files for N1MATD intermediates
         CALL WCLOSE2(LUFEID,FNFEID,'DELETE')
         CALL WCLOSE2(LUN1D,FNN1D,'DELETE')
         IF (.NOT.SKIPGEI) THEN
            CALL WCLOSE2(LUGEID,FNGEID,'DELETE')
         END IF
      END IF
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
      CALL WCLOSE2(LUDKBCR4,FNDKBCR4,'DELETE')
     
C-------------
C     End
C-------------
C
      CALL QEXIT('HMAT')
C
      RETURN
      END


