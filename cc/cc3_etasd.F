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
      SUBROUTINE CC3_ETASD(LISTL,IDLSTL,FOCKY,ISYFKY,
     *                     FREQ,FOCK0, ETA1,ETA2,ETA1EFF,
     *                     ETA2EFF,ISYRES,WORK,LWORK,LUDKBC,FNDKBC,
     *                     LUCKJD,FNCKJD, 
     *                     LUTOC,FNTOC,LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                     LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples component of eta vector 
*             projected into the singles and doubles space
*             
*    ETA3 = ( <L2|[Y,tau3]|HF> + <L3|[Y^,tau3]|HF>)/ (w+epsiln(tau3))
*
*    ETA1EFF(ai) = ETA1(ai) + < ETA3|[[H^,T2],tau(ai)]|HF>     
*
*    ETA2EFF(aibj) = ETA2(aibj) + < ETA3|[H^,tau(aibj)]|HF>     
*
*
*    To these we add
*
*    ETA1EFF(ai) = ETA1EFF(ai) + <L3|[[X,tau_ai],T_3]|HF>
*
*    ETA2EFF(aibj) = ETA2EFF(aibj) + <L3|[[X,E_aiE_bj],T_2]|HF>
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
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
      CHARACTER LISTL*3, MODEL*10
C
      CHARACTER*6 FNGEI,FNFEI
      CHARACTER*5 FNN1
      PARAMETER( FNGEI = 'N1_GEI' , FNFEI = 'N1_FEI' , FNN1 = 'N1MAT' )
      INTEGER LUGEI,LUFEI,LUN1
C
      CHARACTER*(*) FNCKJD, FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X 
      CHARACTER*(*)  FNDKBC
      CHARACTER*13 FN4V
      CHARACTER*11 FNDKBC4
      CHARACTER*1 CDUMMY
      LOGICAL   LOCDBG
      INTEGER   LU4V
      INTEGER   IDLSTL, LWORK
      INTEGER   ISYFKY, ISYRES
      INTEGER   LUCKJD, LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER   KLAMP0,KLAMH0,KFOCKD,KFCKBA,KT2TP,KL1AM,KL2TP
      INTEGER   KEND0,LWRK0,KT1AMP0,KEND1,LWRK1
      INTEGER   ISYCTR, ISYMC, ISYMK, KOFF1, KOFF2
      INTEGER   KXIAJB, KTROC01, KTROC21, KTROC03, KTROC23
      INTEGER   KINTOC, KEND2, LWRK2
      INTEGER   IOPT, LENGTH, ISYCKB, ISYMD, ISYMBD
      INTEGER   KTRVI4, KTRVI5, KTRVI7, KEND3, LWRK3
      INTEGER   KTRVI14, KTRVI15, KTRVI18, KTRVI19, KINTVI, KTRVI6
      INTEGER   KEND4, LWRK4, IOFF, ISYMB, ISYALJ, ISYALJ2, ISCKIJ
      INTEGER   ISYCKD, LENSQ
      INTEGER   KSMAT2, KUMAT2, KDIAG, KINDSQ, KINDEX, KINDEX2, KTMAT
      INTEGER   KTRVI16, KTRVI17, KTRVI20, KSMAT4, KUMAT4, KTRVI11
      INTEGER   KTRVI12, KTRVI13, KEND5, LWRK5
      INTEGER   KWMAT, ISWMAT, ISYMIM
      INTEGER   KINDSQW, LENSQW
      INTEGER   KOIOOOO, KOIOVVO, KOIOOVV, KVVVV, ISCKB2, KDIAGW
      INTEGER   LUDKBC4, LUDKBC
      INTEGER   ISINT1, ISINT2, KRBJIA, KTROC, KTROC1, ISCKB1
      INTEGER   KTRVI, KTRVI1
      INTEGER   KRAIJB,KFCKYCK,ISYMT2
C
      INTEGER ISYMN1,ISYMN2,KN2MAT,KINDSQN,LENSQN,ISGEI,ISFEI,KGEI,KFEI
      INTEGER IADR
C
      INTEGER kx3am
C

#if defined (SYS_CRAY)
      REAL      FREQ, HALF, DDOT
      REAL      FOCKY(*), FOCK0(*), WORK(LWORK)
      REAL      ETA1(*), ETA2(*), ETA1EFF(*), ETA2EFF(*)
      REAL      XNORM
#else
      DOUBLE PRECISION      FREQ, HALF, DDOT
      DOUBLE PRECISION      FOCKY(*), FOCK0(*), WORK(LWORK)
      DOUBLE PRECISION      ETA1(*), ETA2(*), ETA1EFF(*), ETA2EFF(*)
      DOUBLE PRECISION      XNORM
#endif
C
      PARAMETER(HALF = 0.5D0)
C
      CALL QENTER('CC3_ETASD')
C
C
C     CALL DCOPY(NT1AM(ISYRES),ETA1,1,ETA1EFF,1)
C     CALL DCOPY(NT2AM(ISYRES),ETA2,1,ETA2EFF,1)
C
C----------------------------------------------------
C     Initialise character strings and open files
C----------------------------------------------------
C
      CDUMMY = ' '

      LU4V     = -1
C
      FN4V     = 'CC3_ETASD_TMP'
C
      IF (LVVVV) THEN
         CALL WOPEN2(LU4V,FN4V,64,0)
      END IF
C
      LUDKBC4 = -1
      FNDKBC4 = 'CC3_L3_TMP1'
C
      CALL WOPEN2(LUDKBC4,FNDKBC4,64,0)
C
      IF (.NOT.LVVVV) THEN
         !Open files for N1MAT intermediates
         LUGEI = -1
         LUFEI = -1
         LUN1  = -1
         CALL WOPEN2(LUGEI,FNGEI,64,0)
         CALL WOPEN2(LUFEI,FNFEI,64,0)
         CALL WOPEN2(LUN1,FNN1,64,0)
      END IF



C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
 
      IF (LISTL(1:3).EQ.'L0 ') THEN
        ISYCTR = ISYM0
      ELSE
        ! ups, probably higher-order response, not yet implemented here
        CALL QUIT('Unknown type of left vector in CC3_ETASD.')
      END IF

C
C-------------------------------------------------------
C     initial allocations, orbital energy, fock matrix and T2 and L2 :
C-------------------------------------------------------
C
C     Symmetry of integrals in contraction:
C
      ISINT1 = 1
      ISINT2 = 1
      ISYMT2 = ISYM0
      ISYMIM = MULD2H(ISYCTR,ISYFKY)
      IF (.NOT.LVVVV) THEN
         !Symmetries for N1 and N2 intermediates
         ISYMN1 = MULD2H(ISYMIM,ISYMT2) 
         ISYMN2 = MULD2H(ISYMIM,ISYMT2) 
      END IF
C    
      IF (LVVVV) THEN 
         KRBJIA = 1
      ELSE 
         KN2MAT = 1
         KRBJIA = KN2MAT + NCKIJ(ISYMN2)
      END IF
C
      KRAIJB = KRBJIA   + NT2SQ(ISYRES)
      KLAMP0 = KRAIJB   + NT2SQ(ISYRES)
      KLAMH0  = KLAMP0  + NLAMDT
      KFOCKD  = KLAMH0  + NLAMDT
      KFCKBA  = KFOCKD  + NORBTS
      KFCKYCK = KFCKBA  + NT1AMX 
      KT2TP   = KFCKYCK + NT1AM(ISYFKY)
      KL1AM   = KT2TP   + NT2SQ(ISYM0)
      KL2TP   = KL1AM   + NT1AM(ISYCTR)
      KEND0   = KL2TP   + NT2SQ(ISYCTR)
      LWRK0   = LWORK   - KEND0
C
      IF (.NOT.LVVVV) THEN
         KINDSQN = KEND0
         KEND0  = KINDSQN + (6*NCKIJ(ISYMN2) - 1)/IRAT + 1
         LWRK0  = LWORK - KEND0 
      END IF
C
      KT1AMP0 = KEND0
      KEND1   = KT1AMP0 + NT1AMX
      LWRK1   = LWORK   - KEND1

      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
      CALL DZERO(WORK(KRAIJB),NT2SQ(ISYRES))
C
      IF (.NOT.LVVVV) THEN
         CALL DZERO(WORK(KN2MAT),NCKIJ(ISYMN2))
C
         !index array for N2
         LENSQN = NCKIJ(ISYMN2)
         CALL CC3_INDSQ(WORK(KINDSQN),LENSQN,ISYMN2)
      END IF
C
C-------------------------------------
C     Read in lamdap and lamdh
C-------------------------------------
C
      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),DUMMY)

      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

C
C-----------------------------------------------------------
C     Calculate 2*C-E and store 
C     FNDKBC3, FN3FOP and FN3FOP2 for f.o.p. later. 
C     That's option IOPT = 1
C
C     With option IOPT = 2 it just calculates LUDKBC4.
C-----------------------------------------------------------
C
      CALL CC3_TCME(WORK(KLAMP0),ISINT1,WORK(KEND1),LWRK1,IDUMMY,CDUMMY,
     *              LUDKBC,FNDKBC,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *              IDUMMY,CDUMMY,LUDKBC4,FNDKBC4,2)


C
C-------------------------------------
C     Read T2 amplitudes 
C-------------------------------------
C
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KT2TP))

      CALL CC_T2SQ(WORK(KT2TP),WORK(KEND0),ISYM0)

      CALL CC3_T2TP(WORK(KT2TP),WORK(KEND0),ISYM0)
C
C     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
C    *    DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)
C
C-------------------------------------
C     Read L2 amplitudes 
C-------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KL2TP))

      CALL CC_T2SQ(WORK(KL2TP),WORK(KEND0),ISYCTR)

      CALL CC3_T2TP(WORK(KL2TP),WORK(KEND0),ISYCTR)

C     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
C    *    DDOT(NT2SQ(ISYCTR),WORK(KL2TP),1,WORK(KL2TP),1)

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
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT If we want to sum the T3 amplitudes
COMMENT COMMENT COMMENT
      if (.false.) then
         kx3am  = kend0
         kend0 = kx3am + nrhft*nrhft*nrhft*nvirt*nvirt*nvirt
         call dzero(work(kx3am),nrhft*nrhft*nrhft*nvirt*nvirt*nvirt)
         lwrk0 = lwork - kend0
         if (lwrk0 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend0
            call quit('Insufficient space (T3) in CC3_etasd')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT

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

C     IF (LOCDBG) THEN
C        CALL AROUND('In CC3_ETASD: Fock MO matrix (sort)')
C        CALL CC_PRFCKMO(WORK(KFCKBA),ISYM0) 
C     ENDIF

C-----------------------------------------------------
C     Sort the FOCKY matrix
C-----------------------------------------------------
C
      DO ISYMC = 1,NSYM

         ISYMK = MULD2H(ISYMC,ISYFKY)

         DO K = 1,NRHF(ISYMK)

            DO C = 1,NVIR(ISYMC)

               KOFF1 = IFCVIR(ISYMK,ISYMC) +
     *                 NORB(ISYMK)*(C - 1) + K
               KOFF2 = IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C

               WORK(KFCKYCK-1+KOFF2) = FOCKY(KOFF1)

            ENDDO
         ENDDO
      ENDDO

*     IF (LOCDBG) THEN
*        CALL AROUND('In CC3_ETASD: FOCKY MO matrix (sort)')
*        CALL CC_PRFCKMO(WORK(KFCKYCK),ISYFKY)
*     ENDIF


C
C--------------------------------------------------------------
C     Calculate the normal g^0 integrals for
C     OOOO, OVVO and OOVV integrals. VVVV is stored on file
C--------------------------------------------------------------
C
      KOIOOOO = KEND0
      KOIOVVO = KOIOOOO + N3ORHF(ISYM0)
      KOIOOVV = KOIOVVO + NT2SQ(ISYM0)
      KEND0   = KOIOOVV + NT2SQ(ISYM0)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_FMAT (g^0[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KOIOVVO),NT2SQ(ISYMOP))
      CALL DZERO(WORK(KOIOOVV),NT2SQ(ISYMOP))
C
      CALL CC3_INT(WORK(KOIOOOO),WORK(KOIOOVV),WORK(KOIOVVO),
     *             ISYMOP,LU4V,FN4V,
     *             .FALSE.,'DUMMY',IDUMMY,IDUMMY,
     *             .FALSE.,'DUMMY',IDUMMY,IDUMMY,
     *             WORK(KEND0),LWRK0)
C
C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
      KTROC   = KEND0
      KTROC1  = KTROC   + NTRAOC(ISINT2)
      KXIAJB  = KTROC1  + NTRAOC(ISINT2)
      KEND0   = KXIAJB  + NT2AM(ISYM0)

      KTROC01 = KEND0
      KTROC21 = KTROC01 + NTRAOC(ISYM0)
      KTROC03 = KTROC21 + NTRAOC(ISYM0)
      KTROC23 = KTROC03 + NTRAOC(ISYM0)
      KEND0   = KTROC23 + NTRAOC(ISYM0)
      LWRK0   = LWORK   - KEND0

      KINTOC  = KEND0
      KEND2   = KINTOC + MAX(NTOTOC(ISYM0),NTOTOC(ISINT2))
      LWRK2   = LWORK  - KEND2

      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_ETASD')
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
C------------------------------------------------------------
C     Read in integrals used in contractions and transform.
C------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC),WORK(KLAMP0),
     *                    WORK(KEND2),LWRK2,ISINT2)
C
      CALL CCFOP_SORT(WORK(KTROC),WORK(KTROC1),ISINT2,1)
C
      CALL CC3_LSORT1(WORK(KTROC),ISINT2,WORK(KEND2),LWRK2,5)

C
C------------------------
C     Occupied integrals.
C------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYM0) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF

*     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of OCC-INT ',
*    *    DDOT(NTOTOC(ISYM0),WORK(KINTOC),1,WORK(KINTOC),1)
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
 
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC01),WORK(KLAMH0),
     *               WORK(KEND2),LWRK2,ISYM0)
*     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',
*    *    DDOT(NTOTOC(ISYM0),WORK(KINTOC),1,WORK(KINTOC),1)
 
      CALL CCSDT_TCMEOCC(WORK(KTROC01),WORK(KTROC21),ISYM0)
 
      CALL CCFOP_SORT(WORK(KTROC01),WORK(KTROC03),ISYM0,1)
 
      CALL CCFOP_SORT(WORK(KTROC21),WORK(KTROC23),ISYM0,1)

C
C-----------------------------------------------------
C     Calculate <L3|[[X,tau_ai],T_3]|HF>
C-----------------------------------------------------     
C
      CALL CC3_ETA_1(WORK(KFCKYCK),ISYFKY,ETA1,ISYRES,
     *                     WORK(KEND2),LWRK2)

C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKB  = MULD2H(ISYMD,ISYM0)
         ISCKB1 = MULD2H(ISINT1,ISYMD)
         ISCKB2 = MULD2H(ISYMD,ISYM0)
C
         IF (.NOT.LVVVV) THEN
            !Symmetry of arrays needed to construct N1MAT
            ISGEI  = MULD2H(ISYMN1,ISYMD)
            ISFEI  = MULD2H(ISYMN1,ISYMD)
         END IF
C
C        IF (LOCDBG) WRITE(LUPRI,*) 'In CC3_ETASDPJ: ISYCKB :',ISYCKB
C
         KTRVI  = KEND0
         KTRVI1 = KTRVI  + NCKATR(ISCKB1)
         KEND1  = KTRVI1 + NCKATR(ISCKB1)
         LWRK1  = LWORK - KEND1
C
         IF (LVVVV) THEN
            KVVVV = KEND1
            KEND1 = KVVVV +  NMAABC(ISCKB2)
            LWRK1 = LWORK - KEND1
         END IF

C
         IF (.NOT.LVVVV) THEN
            !Arrays needed to construct N1MAT
            KGEI  = KEND1
            KFEI  = KGEI  + NCKATR(ISGEI)
            KEND1 = KFEI  + NCKATR(ISFEI)
            LWRK1 = LWORK - KEND1
         END IF
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in CC3_ETASD')
         END IF
C



         DO D = 1,NVIR(ISYMD)
C
            IF (.NOT.LVVVV) THEN
               CALL DZERO(WORK(KGEI),NCKATR(ISGEI))
               CALL DZERO(WORK(KFEI),NCKATR(ISFEI))
            END IF

C
            IF (LVVVV) THEN
C--------------------------------------------
C              Read in g^0_{vvvv} for a given D
C--------------------------------------------
C
               IF (NMAABC(ISCKB2) .GT. 0) THEN
                  IOFF = I3VVIR(ISCKB2,ISYMD)
     *                 + NMAABC(ISCKB2)*(D-1)
     *                 + 1
                  CALL GETWA2(LU4V,FN4V,WORK(KVVVV),IOFF,NMAABC(ISCKB2))
               ENDIF
C
            END IF

C
C------------------------------------------------------------
C           Read and transform integrals used in contraction.
C------------------------------------------------------------
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
               IOFF = ICKBD(ISCKB1,ISYMD) + NCKATR(ISCKB1)*(D - 1) + 1
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI),IOFF,
     &                     NCKATR(ISCKB1))
            ENDIF
C
            IF (LWRK1 .LT. NCKATR(ISCKB1)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_ETASD (TRVI)')
            END IF
C
            CALL CCSDT_SRVIR3(WORK(KTRVI),WORK(KEND1),ISYMD,D,ISINT1)
C
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
               IOFF = ICKBD(ISCKB1,ISYMD) + NCKATR(ISCKB1)*(D - 1) + 1
               CALL GETWA2(LUDKBC4,FNDKBC4,WORK(KTRVI1),IOFF,
     &                     NCKATR(ISCKB1))
            ENDIF
C
            IF (LWRK1 .LT. NCKATR(ISCKB1)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_ETASD (TRVI1)')
            END IF
C
            CALL CCSDT_SRVIR3(WORK(KTRVI1),WORK(KEND1),ISYMD,D,ISINT1)
C

C
C
C           ------------------
C           Memory allocation.
C           ------------------
            KTRVI4  = KEND1
            KTRVI5  = KTRVI4 + NCKATR(ISYCKB)
            KTRVI7  = KTRVI5 + NCKATR(ISYCKB)
            KEND3   = KTRVI7 + NCKATR(ISYCKB)
            LWRK3   = LWORK  - KEND3
           
            KTRVI14 = KEND3
            KTRVI15 = KTRVI14 + NCKATR(ISYCKB)
            KTRVI18 = KTRVI15 + NCKATR(ISYCKB)
            KTRVI19 = KTRVI18 + NCKATR(ISYCKB)
            KEND3   = KTRVI19 + NCKATR(ISYCKB)
            LWRK3   = LWORK  - KEND3

           
            KINTVI = KEND3
            KTRVI6 = KINTVI + MAX(NCKA(ISYMD),NCKA(ISYCKB))
            KEND4  = KTRVI6 + NCKATR(ISYCKB) 
            LWRK4  = LWORK  - KEND4
           
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CC3_ETASD')
            END IF
C
C           -------------------------------------------
C           Read 2*C-E of integral used for t3-bar
C           -------------------------------------------
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI4),IOFF,
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
C           Sort the integrals for t3-bar
C           ------------------------------------------------
C

            CALL CCSDT_SRTVIR(WORK(KTRVI4),WORK(KTRVI5),WORK(KEND4),
     *                        LWRK4,ISYMD,ISYM0)
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

            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI19),WORK(KLAMP0),
     *                       ISYMD,D,ISYM0,WORK(KEND4),LWRK4)

            IF (LWRK4 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     *                   'CC3_ETASD  (CC3 TRVI)')
            END IF

            CALL CCSDT_SRTVIR(WORK(KTRVI19),WORK(KTRVI18),WORK(KEND4)
     *                        ,LWRK4,ISYMD,ISYM0)

            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI6),IOFF,
     *                     NCKATR(ISYCKB))
            ENDIF

            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_ETASD (2)')
            END IF

            CALL DCOPY(NCKATR(ISYCKB),WORK(KTRVI6),1,WORK(KTRVI7),1)
C
C
            DO ISYMB = 1,NSYM

               ISYALJ  = MULD2H(ISYMB,ISYM0)
               ISYALJ2 = MULD2H(ISYMD,ISYM0)
               ISYMBD  = MULD2H(ISYMD,ISYMB)
               ISCKIJ  = MULD2H(ISYMBD,ISYCTR)
               ISWMAT  = MULD2H(ISCKIJ,ISYFKY)
               ISYCKD  = MULD2H(ISYM0,ISYMB)

C              Can use kend3 since we do not need the integrals anymore.
               KSMAT2  = KEND3
               KUMAT2  = KSMAT2  + NCKIJ(ISCKIJ)
               KDIAG   = KUMAT2  + NCKIJ(ISCKIJ)
               KDIAGW  = KDIAG   + NCKIJ(ISCKIJ)
               KINDSQ  = KDIAGW  + NCKIJ(ISWMAT)
               KINDSQW = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX  = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ)  - 1)/IRAT + 1
               KTMAT   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KWMAT   = KTMAT   + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWMAT))
               KEND4   = KWMAT   + NCKIJ(ISWMAT)
               LWRK4   = LWORK   - KEND4

               KTRVI16 = KEND4
               KTRVI17 = KTRVI16 + NCKATR(ISYCKD)
               KTRVI20 = KTRVI17 + NCKATR(ISYCKD)
               KEND4   = KTRVI20 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4

               KSMAT4  = KEND4
               KUMAT4  = KSMAT4  + NCKIJ(ISCKIJ)
               KTRVI11 = KUMAT4  + NCKIJ(ISCKIJ)
               KTRVI12 = KTRVI11 + NCKATR(ISYCKD)
               KTRVI13 = KTRVI12 + NCKATR(ISYCKD)
               KEND4   = KTRVI13 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4


               KINTVI  = KEND4
               KEND5   = KINTVI  + MAX(NCKA(ISYMB),NCKA(ISYCKD))
               LWRK5   = LWORK   - KEND5

               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CC3_ETASD')
               END IF
C
C
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG), WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)

C              IF (LOCDBG) THEN
C                WRITE(LUPRI,*) 'Norm of DIA  ',
C    *              DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,WORK(KDIAG),1)
C              END IF
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

               DO B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
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
                     CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI11),
     *                           IOFF,NCKATR(ISYCKD))
                  ENDIF

                  CALL CCSDT_SRTVIR(WORK(KTRVI11),WORK(KTRVI12),
     *                              WORK(KEND5),LWRK5,ISYMB,
     *                              ISYM0)

                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B - 1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI13),IOFF,
     *                           NCKATR(ISYCKD))
                  ENDIF

                  IOFF = ICKAD(ISYCKD,ISYMB) 
     *                 + NCKA(ISYCKD)*(B - 1) + 1
                  IF (NCKA(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                           NCKA(ISYCKD))
                  ENDIF

                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI20),
     *                             WORK(KLAMP0),ISYMB,B,ISYM0,
     *                             WORK(KEND4),LWRK4)
C
C                 ----------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for B,D for T3-BAR:
C                 ----------------------------------------------------
                  CALL DZERO(WORK(KSMAT2),NCKIJ(ISCKIJ))
 
                  CALL CCFOP_SMAT(0.0D0,WORK(KL1AM),ISYCTR,WORK(KL2TP),
     *                            ISYCTR,WORK(KTMAT),
     *                            WORK(KFCKBA),WORK(KXIAJB),ISYM0,
     *                            WORK(KTRVI14),WORK(KTRVI15),
     *                            WORK(KTRVI4),WORK(KTRVI5),
     *                            WORK(KTROC01),WORK(KTROC21),
     *                            ISYM0,WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KSMAT2),WORK(KEND4),LWRK4,
     *                            WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                            ISYMB,B,ISYMD,D)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT2),1)

                  CALL T3_FORBIDDEN(WORK(KSMAT2),ISYCTR,ISYMB,B,ISYMD,D)
C
C                 ----------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for D,B for T3-BAR:
C                 ----------------------------------------------------
                  CALL DZERO(WORK(KSMAT4),NCKIJ(ISCKIJ))

                  CALL CCFOP_SMAT(0.0D0,WORK(KL1AM),ISYCTR,WORK(KL2TP),
     *                            ISYCTR,WORK(KTMAT),WORK(KFCKBA),
     *                            WORK(KXIAJB),ISYM0,
     *                            WORK(KTRVI16),WORK(KTRVI17),
     *                            WORK(KTRVI11),WORK(KTRVI12),
     *                            WORK(KTROC01),WORK(KTROC21),
     *                            ISYM0,WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KSMAT4),WORK(KEND4),LWRK4,
     *                            WORK(KINDEX2),WORK(KINDSQ),
     *                            LENSQ,ISYMD,D,ISYMB,B)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT4),1)
C
C
                  CALL T3_FORBIDDEN(WORK(KSMAT4),ISYCTR,ISYMD,D,ISYMB,B)
C
C                 ------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d for t3-bar.
C                 ------------------------------------------------
                  CALL DZERO(WORK(KUMAT2),NCKIJ(ISCKIJ))

                  CALL CCFOP_UMAT(0.0D0,WORK(KL1AM),ISYCTR,WORK(KL2TP),
     *                            ISYCTR,
     *                            WORK(KXIAJB),ISYM0,WORK(KFCKBA),
     *                            WORK(KTRVI19),WORK(KTRVI7),
     *                            WORK(KTROC03),WORK(KTROC23),ISYM0,
     *                            WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KUMAT2),
     *                            WORK(KTMAT),WORK(KEND4),LWRK4,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KUMAT2),1)
C
C
                  CALL T3_FORBIDDEN(WORK(KUMAT2),ISYCTR,ISYMB,B,ISYMD,D)
C
C                 ------------------------------------------------
C                 Calculate U(ci,jk) for fixed d,b for t3-bar.
C                 ------------------------------------------------
                  CALL DZERO(WORK(KUMAT4),NCKIJ(ISCKIJ))

                  CALL CCFOP_UMAT(0.0D0,WORK(KL1AM),ISYCTR,WORK(KL2TP),
     *                            ISYCTR,WORK(KXIAJB),ISYM0,
     *                            WORK(KFCKBA),WORK(KTRVI20),
     *                            WORK(KTRVI13),WORK(KTROC03),
     *                            WORK(KTROC23),ISYM0,
     *                            WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KUMAT4),WORK(KTMAT),
     *                            WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                            LENSQ,ISYMD,D,ISYMB,B)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KUMAT4),1)

                  CALL T3_FORBIDDEN(WORK(KUMAT4),ISYCTR,ISYMD,D,ISYMB,B)
C
C                 -------------------------------------------
C                 Sum up the S-bar and U-bar to get a real T3-bar
C                 -------------------------------------------
                  CALL CC3_T3BD(ISCKIJ,WORK(KSMAT2),WORK(KSMAT4),
     *                                 WORK(KUMAT2),WORK(KUMAT4),
     *                                 WORK(KTMAT),WORK(KINDSQ),LENSQ)
C
C                 -------------------------------------
C
C----------------------------------------------------------
C                 Calculate <L3|[[X,E_aiE_bj],T_2]|HF>
C----------------------------------------------------------
C
*                 ISYMT2 = ISYM0
C
                  CALL CC3_ETA_2(WORK(KTMAT),ISCKIJ,WORK(KFCKYCK),
     *                           ISYFKY,WORK(KT2TP),ISYMT2,ETA2,
     *                           ISYRES,WORK(KRAIJB),WORK(KINDSQ),
     *                           LENSQ,WORK(KINDEX2),ISYALJ2,
     *                           ISYMB,B,ISYMD,D,WORK(KEND4),LWRK4)
C
C    calculate     <L3|[Y^,tau3]|HF>
C
                  CALL WBARBD_V(WORK(KTMAT),ISCKIJ,FOCKY,ISYFKY,
     *                 WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)
C
                  CALL WBARBD_O(WORK(KTMAT),ISCKIJ,FOCKY,ISYFKY,
     *                 WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)
 
C    calculate     <L2|[Y,tau3]|HF>
C
                  CALL WBARBD_T2(B,ISYMB,D,ISYMD,WORK(KL2TP),ISYCTR,
     *                           FOCKY,ISYFKY,WORK(KWMAT),ISWMAT)
C
COMMENT COMMENT COMMENT
C       call sum_pt3(work(kwmat),isymb,b,isymd,d,
C    *             iswmat,work(kx3am),4)
COMMENT COMMENT COMMENT
C
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,-FREQ,
     *                         ISWMAT,WORK(KWMAT),
     *                         WORK(KDIAGW),WORK(KFOCKD))
C
                  CALL T3_FORBIDDEN(WORK(KWMAT),ISYFKY,ISYMB,B,ISYMD,D)
C

COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C       call sum_pt3(work(kwmat),isymb,b,isymd,d,
C    *             iswmat,work(kx3am),4)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
*                 ISYMIM = MULD2H(ISWMAT,ISYMBD)
C
C-----------------------------------------
C Calculate <E_3|[H,\tau_nu_2]|HF>
C-----------------------------------------
C
                  CALL CC3_W3_CY2V(ETA2EFF,ISYRES,WORK(KRBJIA),
     *                             WORK(KWMAT),ISWMAT,
     *                             WORK(KTMAT),WORK(KTRVI),WORK(KTRVI1),
     *                             ISINT1,WORK(KEND4),LWRK4,
     *                             WORK(KINDSQW),LENSQW,
     *                             ISYMB,B,ISYMD,D,.TRUE.)

                  CALL CC3_W3_CY2O(ETA2EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                             WORK(KTMAT),WORK(KTROC),WORK(KTROC1),
     *                             ISINT1,WORK(KEND4),LWRK4,
     *                             WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D,
     *                             .TRUE.)

C
C-----------------------------------------
C Calculate <E_3|[[H,T^0_2],\tau_nu_1]|HF>
C-----------------------------------------
C
                  IF (LVVVV) THEN
                    CALL CC3_W3_OMEGA1(ETA1EFF,ISYRES,WORK(KWMAT),
     *                                 WORK(KTMAT),ISYMIM,
     *                                 WORK(KOIOOOO),WORK(KOIOVVO),
     *                                 WORK(KOIOOVV),WORK(KVVVV),ISYM0,
     *                                 WORK(KT2TP),ISYM0,
     *                                 WORK(KEND4),LWRK4,
     *                                 LENSQW,WORK(KINDSQW),
     *                                 ISYMB,B,ISYMD,D,.TRUE.)
                   ELSE  
                     CALL DSCAL(NCKIJ(ISWMAT),-1.0D0,WORK(KWMAT),1)

                     !Construct N1 and N2 intermediates
                     CALL WT2_N1N2(WORK(KWMAT),ISYMIM,
     *                       WORK(KT2TP),ISYM0,
     *                       WORK(KGEI),WORK(KFEI),
     *                       ISYMN1,
     *                       WORK(KN2MAT),ISYMN2,
     *                       B,ISYMB,D,ISYMD,
     *                       WORK(KINDSQW),LENSQW,
     *                       WORK(KINDSQN),LENSQN,
     *                       WORK(KEND4),LWRK4,
     *                       .TRUE.)

                   END IF

               ENDDO   ! B
            ENDDO      ! ISYMB

            IF (.NOT.LVVVV) THEN
C
C              ----------------------------------------------------------
C              Put KGEI(ge,i)^F and KFEI(fe,i)^G (which are intermediates
C              for N1MAT(fge,i) ) to files (for fixed F=D and G=D).
C              ----------------------------------------------------------

               !Put KGEI to file as (gei,F)   (fixed F corresponds to D)
               IADR = ICKBD(ISGEI,ISYMD) + NCKATR(ISGEI)*(D-1) + 1
               CALL PUTWA2(LUGEI,FNGEI,WORK(KGEI),IADR,NCKATR(ISGEI))
C
               !Put KFEI to file as (fei,G)   (fixed G corresponds to D)
               IADR = ICKBD(ISFEI,ISYMD) + NCKATR(ISFEI)*(D-1) + 1
               CALL PUTWA2(LUFEI,FNFEI,WORK(KFEI),IADR,NCKATR(ISFEI))
C
            END IF


         ENDDO       ! D
      ENDDO          ! ISYMD 
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Vccupied  cont ) 
C     in ETA2EFF 
C
C     Accumulate RAIJB from <L3|[[X,E_aiE_bj],T_2]|HF> 
C     in ETA2EFF
C------------------------------------------------------
C
      CALL CC3_RBJIA(ETA2EFF,ISYRES,WORK(KRBJIA))
      CALL CC3_RBJIA(ETA2EFF,ISYRES,WORK(KRAIJB))

C
      IF (.NOT.LVVVV) THEN
C
         !Read (gei,F) and (fei,G) intermediates from files
         !add them and put the result to a file as (fge,I)
         CALL N1_RESORT(ISYMN1,LUN1,FNN1,LUGEI,FNGEI,LUFEI,FNFEI,
     *                  WORK(KEND0),LWRK0,.FALSE.)
C
         !Calculate <T3|[[H,T2],tau_ai]|HF> except VVVV contribution
         CALL N1N2_G(LUN1,FNN1,
     *                     ISYMN1,
     *                     WORK(KN2MAT),ISYMN2,
     *                     WORK(KOIOVVO),WORK(KOIOOVV),
     *                     WORK(KOIOOOO),ISYM0,
     *                     ETA1EFF,ISYRES,
     *                     WORK(KINDSQN),LENSQN,
     *                     WORK(KEND0),LWRK0)
C
         !Calculate VVVV contribution to <T3|[[H,T2],tau_ai]|HF>
         IOPT = 0 !normal Lambda matrices used in backtransformation
         CALL  N1_GV4(IOPT,
     *                LUN1,FNN1,
     *                ISYMN1,
     *                WORK(KLAMP0),1,
     *                WORK(KLAMP0),1,
     *                WORK(KLAMH0),1,
     *                WORK(KLAMH0),1,
     *                ETA1EFF,ISYRES,
     *                WORK(KEND0),LWRK0)
C
      END IF
C
COMMENT COMMENT
C      write(lupri,*) 'WMAT in CC3_ETASD : '
C      call print_pt3(work(kx3am),ISYFKY,4)
COMMENT COMMENT
C
C-------------------------------
C     Close and delete files
C-------------------------------
C
      IF (LVVVV) THEN
         CALL WCLOSE2(LU4V,FN4V,'DELETE')
      END IF
C
      CALL WCLOSE2(LUDKBC4,FNDKBC4,'DELETE')
C
      IF (.NOT.LVVVV) THEN
         !Close files for N1MAT intermediates
         CALL WCLOSE2(LUGEI,FNGEI,'DELETE')
         CALL WCLOSE2(LUFEI,FNFEI,'DELETE')
         CALL WCLOSE2(LUN1,FNN1,'DELETE')
      END IF
C
C------------------------------------------
C     Accumulate CCSD and CC3 contributions
C------------------------------------------
C
      DO I = 1,NT2AM(ISYRES)
         ETA2EFF(I) = ETA2EFF(I) + ETA2(I)
      END DO
C
      DO I = 1,NT1AM(ISYRES)
         ETA1EFF(I) = ETA1EFF(I) + ETA1(I)
      END DO
C

C
C-------------
C     End
C-------------
C
C
      CALL QEXIT('CC3_ETASD')
C
      RETURN
      END

C  /* Deck wbarbd_v */
      SUBROUTINE WBARBD_V(TMAT,ISTMAT,FOCKY,ISYFKY, 
     *                 WMAT,ISWMAT,WRK,LWRK)
C 
C WBD(a,i,k,j) = WBD(a,i,k,j) + sum (f)  focky(f,a)*tmatBD(f,i,k,j) 
C
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C

C
      IMPLICIT NONE
C
      INTEGER LWRK, KFCAF, KEND0, LWRK0, KOFF1, KOFF2 
      INTEGER NF, KOFFY, KOFFT, KOFFW 
      INTEGER ISTMAT, ISYFKY, ISWMAT, ISFIKJ, ISYFIK 
      INTEGER ISYMA, ISYAI, ISYAIK, NA
      INTEGER ISYMF, ISYMJ, ISYMK, ISYMI, ISYFI
#if defined (SYS_CRAY)
      REAL TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      REAL HALF, ONE
#else
      DOUBLE PRECISION TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      DOUBLE PRECISION HALF, ONE
#endif
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('WBARBD_V')
C
C
C RESORT VIR-VIR  FOCKY ELEMENTS (A,F)
C
C
      KFCAF  = 1
      KEND0  = KFCAF + NMATAB(ISYFKY)
      LWRK0  = LWRK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WBARBD_V')
      END IF
C
      DO ISYMF = 1,NSYM
         ISYMA = MULD2H(ISYMF,ISYFKY)
         DO F = 1,NVIR(ISYMF)
            KOFF1 = IFCVIR(ISYMA,ISYMF) + NORB(ISYMA)*(F - 1)
     *                                  + NRHF(ISYMA) + 1
            KOFF2 = KFCAF + IMATAB(ISYMA,ISYMF) + NVIR(ISYMA)*(F - 1) 
            CALL DCOPY(NVIR(ISYMA),FOCKY(KOFF1),1,WRK(KOFF2),1)
         END DO
      END DO
C
C CARRY OUT MATRIX MULTIPLICATION
C WBD(a,i,k,j) = WBD(a,i,k,j) + sum (f)  focky(f,a)*tmatBD(f,i,k,j) 
C 
      ISFIKJ = ISTMAT
      DO ISYMJ = 1,NSYM
         ISYFIK =MULD2H(ISYMJ,ISFIKJ)
         DO J   = 1,NRHF(ISYMJ)
            DO ISYMK = 1,NSYM
               ISYFI = MULD2H(ISYMK,ISYFIK)
               DO K  = 1,NRHF(ISYMK)
                  DO ISYMI = 1,NSYM
                     ISYMF = MULD2H(ISYFI,ISYMI)   
                     ISYMA = MULD2H(ISYFKY,ISYMF)
                     ISYAIK = MULD2H(ISWMAT,ISYMJ)
                     ISYAI = MULD2H(ISYAIK,ISYMK)
                     NA    = MAX(1,NVIR(ISYMA))
                     NF    = MAX(1,NVIR(ISYMF))
                     KOFFY = KFCAF + IMATAB(ISYMF,ISYMA) 
                     KOFFT = ISAIKJ(ISYFIK,ISYMJ)
     *                      + NCKI(ISYFIK)*(J-1)
     *                      + ISAIK(ISYFI,ISYMK) 
     *                      + NT1AM(ISYFI)*(K-1)
     *                      + IT1AM(ISYMF,ISYMI) + 1
                     KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *                      + NCKI(ISYAIK)*(J-1)
     *                      + ISAIK(ISYAI,ISYMK) 
     *                      + NT1AM(ISYAI)*(K-1)
     *                      + IT1AM(ISYMA,ISYMI) + 1
C
C SYMMETRY BETWEEN BJ AND CK INTRODUCE A FACTOR TWO
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
C WBD(a,i,k,j) = WBD(a,i,k,j) + sum (f)  focky(f,a)*tmatBD(f,i,k,j) 
C
                     CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI), 
     *                          NVIR(ISYMF),-ONE,WRK(KOFFY),NF,
     *                          TMAT(KOFFT),NF,ONE,WMAT(KOFFW),NA) 
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('WBARBD_V')
C
      RETURN
      END
C
C  /* Deck wbarbd_o */
      SUBROUTINE WBARBD_O(TMAT,ISTMAT,FOCKY,ISYFKY,
     *                 WMAT,ISWMAT,WRK,LWRK)
C 
C WBD(a,i,k,j) = WBD(a,i,k,j) - sum (l) tmatBD(a,l,k,j)*focky(i,l)
C
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C

      IMPLICIT NONE
C
      INTEGER LWRK, KFCLI, KEND0, LWRK0, KOFF1, KOFF2 
      INTEGER NI, KOFFY, KOFFT, KOFFW 
      INTEGER ISTMAT, ISYFKY, ISWMAT, ISALKJ
      INTEGER ISYMA, ISYAI, ISYAIK, ISYALK, ISYAL, NA
      INTEGER ISYMJ, ISYMK, ISYMI, ISYML, ISYFI
C
#if defined (SYS_CRAY)
      REAL TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      REAL MHALF, ONE 
#else
      DOUBLE PRECISION TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      DOUBLE PRECISION MHALF, ONE
#endif
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER ( ONE = 1.0D0)
C
      CALL QENTER('WBARBD_O')
C

C
C RESORT OCC-OCC  FOCKY ELEMENTS (L,I)
C
C
      KFCLI  = 1
      KEND0  = KFCLI + NMATIJ(ISYFKY)
      LWRK0  = LWRK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WBARBD_O')
      END IF
C   
      DO ISYMI = 1,NSYM
         ISYML = MULD2H(ISYMI,ISYFKY)
         DO I = 1,NRHF(ISYMI)
             KOFF1 = IFCRHF(ISYML,ISYMI) + NORB(ISYML)*(I - 1) + 1
             KOFF2 = KFCLI + IMATIJ(ISYML,ISYMI) + NRHF(ISYML)*(I - 1)
             CALL DCOPY(NRHF(ISYML),FOCKY(KOFF1),1,WRK(KOFF2),1) 
         END DO
      END DO
C
C CARRY OUT MATRIX MULTIPLICATION
C WBD(a,i,k,j) = WBD(a,i,k,j) - sum (l) tmatBD(a,l,k,j)*focky(i,l)
C 
      ISALKJ = ISTMAT
      DO ISYMJ = 1,NSYM
         ISYALK =MULD2H(ISYMJ,ISALKJ)
         DO J = 1,NRHF(ISYMJ)
            DO ISYMK = 1,NSYM
               ISYAL = MULD2H(ISYMK,ISYALK)
               DO K  = 1,NRHF(ISYMK)
                  DO ISYML = 1,NSYM
                     ISYMA = MULD2H(ISYAL,ISYML)   
                     ISYMI = MULD2H(ISYFKY,ISYML)
                     ISYAIK = MULD2H(ISWMAT,ISYMJ)
                     ISYAI = MULD2H(ISYAIK,ISYMK)
                     NA    = MAX(1,NVIR(ISYMA))
                     NI    = MAX(1,NRHF(ISYMI))
                     KOFFY = KFCLI + IMATIJ(ISYMI,ISYML) 
                     KOFFT = ISAIKJ(ISYALK,ISYMJ)
     *                      + NCKI(ISYALK)*(J-1)
     *                      + ISAIK(ISYAL,ISYMK) 
     *                      + NT1AM(ISYAL)*(K-1)
     *                      + IT1AM(ISYMA,ISYML) + 1
                     KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *                      + NCKI(ISYAIK)*(J-1)
     *                      + ISAIK(ISYAI,ISYMK) 
     *                      + NT1AM(ISYAI)*(K-1)
     *                      + IT1AM(ISYMA,ISYMI) + 1
C
C SYMMETRY BETWEEN BJ AND CK INTRODUCE A FACTOR TWO
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
                     CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),
     *                          NRHF(ISYML),ONE,TMAT(KOFFT),NA,
     *                          WRK(KOFFY),NI,ONE,WMAT(KOFFW),NA)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('WBARBD_O')
C
      RETURN
      END
C
C  /* Deck wbarbd_t2 */
      SUBROUTINE WBARBD_T2(B,ISYMB,D,ISYMD,T2TP,ISYMT2,FOCKY,ISYFKY,
     *                 WMAT,ISWMAT)
C 
C WBD(a,i,k,j) = WBD(a,i,k,j) +
C        
C       focky(j,B)*t2(ai,Dk) - focky(k,B)*t2(ai,Dj) 
C       focky(k,D)*t2(ai,Bj) - focky(j,D)*t2(ai,Bk) 
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C

      IMPLICIT NONE
C
      INTEGER ISYMB, ISYMD, ISYMT2, ISYFKY, ISWMAT 
      INTEGER ISYMJ, KJB, KJD, ISYMK, KKB, KKD, ISYMI, ISYIJ, ISYIK
      INTEGER ISYMA, ISYAI, ISYAIK, ISYAIJ, KAIKD, KAIJD, KAIJB
      INTEGER KAIKB, KAIKJ
C
#if defined (SYS_CRAY)
      REAL T2TP(*), FOCKY(*), WMAT(*)
#else
      DOUBLE PRECISION T2TP(*), FOCKY(*), WMAT(*)
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      CALL QENTER('WBARBD_T2')
C
C
C       focky(j,B)*t2(ai,Dk) - focky(k,B)*t2(ai,Dj) 
C       focky(k,D)*t2(ai,Bj) - focky(j,D)*t2(ai,Bk) 
C
      
C
C (1)   wmat(aikj) = wmat(aikj) +  focky(j,B)*t2(ai,Dk) 
C
      ISYMJ = MULD2H(ISYFKY,ISYMB)
      ISYAIK = MULD2H(ISYMT2,ISYMD)
      DO ISYMK = 1,NSYM
         ISYAI = MULD2H(ISYAIK,ISYMK)
         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYAI,ISYMI)
            DO J = 1,NRHF(ISYMJ)
               KJB = IFCVIR(ISYMJ,ISYMB) + NORB(ISYMJ)*(B - 1) + J
               DO K = 1,NRHF(ISYMK)
                  DO I = 1,NRHF(ISYMI)
                     DO A = 1,NVIR(ISYMA)
                        KAIKD = IT2SP(ISYAIK,ISYMD)
     *                        + NCKI(ISYAIK)*(D-1)
     *                        + ISAIK(ISYAI,ISYMK)
     *                        + NT1AM(ISYAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J-1)
     *                        + ISAIK(ISYAI,ISYMK)
     *                        + NT1AM(ISYAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        WMAT(KAIKJ) = WMAT(KAIKJ) 
     *                              + FOCKY(KJB)*T2TP(KAIKD)    
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
          
C
C (2)  wmat(aikj) = wmat(aikj) - focky(k,B)*t2(ai,Dj) 
C
      ISYMK = MULD2H(ISYFKY,ISYMB)
      ISYAIJ = MULD2H(ISYMT2,ISYMD)
      DO ISYMJ = 1,NSYM
         ISYAI = MULD2H(ISYAIJ,ISYMJ)
         ISYAIK = MULD2H(ISYAI,ISYMK)
         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYAI,ISYMI)
            DO J = 1,NRHF(ISYMJ)
               DO K = 1,NRHF(ISYMK)
               KKB = IFCVIR(ISYMK,ISYMB) + NORB(ISYMK)*(B - 1) + K
                  DO I = 1,NRHF(ISYMI)
                     DO A = 1,NVIR(ISYMA)

                        KAIJD = IT2SP(ISYAIJ,ISYMD)
     *                        + NCKI(ISYAIJ)*(D-1)
     *                        + ISAIK(ISYAI,ISYMJ)
     *                        + NT1AM(ISYAI)*(J-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J-1)
     *                        + ISAIK(ISYAI,ISYMK)
     *                        + NT1AM(ISYAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        WMAT(KAIKJ) = WMAT(KAIKJ) 
     *                              - FOCKY(KKB)*T2TP(KAIJD)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      
C
C (3)  wmat(aikj) = wmat(aikj) + focky(k,D)*t2(ai,Bj) 
C
      ISYMK = MULD2H(ISYFKY,ISYMD)
      ISYAIJ = MULD2H(ISYMT2,ISYMB)
      DO ISYMJ = 1,NSYM
         ISYAI = MULD2H(ISYAIJ,ISYMJ)
         ISYAIK = MULD2H(ISYAI,ISYMK)
         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYAI,ISYMI)
            DO J = 1,NRHF(ISYMJ)
               DO K = 1,NRHF(ISYMK)
                  KKD = IFCVIR(ISYMK,ISYMD) + NORB(ISYMK)*(D - 1) + K
                  DO I = 1,NRHF(ISYMI)
                     DO A = 1,NVIR(ISYMA)

                        KAIJB = IT2SP(ISYAIJ,ISYMB)
     *                        + NCKI(ISYAIJ)*(B-1)
     *                        + ISAIK(ISYAI,ISYMJ)
     *                        + NT1AM(ISYAI)*(J-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J-1)
     *                        + ISAIK(ISYAI,ISYMK)
     *                        + NT1AM(ISYAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        WMAT(KAIKJ) = WMAT(KAIKJ) 
     *                              + FOCKY(KKD)*T2TP(KAIJB)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

C
C (4)  wmat(aikj) = wmat(aikj) -  focky(j,D)*t2(ai,Bk) 
C
      ISYMJ = MULD2H(ISYFKY,ISYMD)
      ISYAIK = MULD2H(ISYMT2,ISYMB)
      DO ISYMK = 1,NSYM
         ISYAI = MULD2H(ISYAIK,ISYMK)
         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYAI,ISYMI)
            DO J = 1,NRHF(ISYMJ)
               KJD = IFCVIR(ISYMJ,ISYMD) + NORB(ISYMJ)*(D - 1) + J
               DO K = 1,NRHF(ISYMK)
                  DO I = 1,NRHF(ISYMI)
                     DO A = 1,NVIR(ISYMA)

                        KAIKB = IT2SP(ISYAIK,ISYMB)
     *                        + NCKI(ISYAIK)*(B-1)
     *                        + ISAIK(ISYAI,ISYMK)
     *                        + NT1AM(ISYAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J-1)
     *                        + ISAIK(ISYAI,ISYMK)
     *                        + NT1AM(ISYAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A

                        WMAT(KAIKJ) = WMAT(KAIKJ)
     *                              - FOCKY(KJD)*T2TP(KAIKB)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL QEXIT('WBARBD_T2')
C
      RETURN
      END
C  /* Deck cc3_w3_omega1 */
      SUBROUTINE CC3_W3_OMEGA1(OMEGA1,ISYRES,WMAT,TMAT,ISYMIM,
     *                         XOOOO,XOVVO,XOOVV,XVVVV,ISYINT,
     *                         T2TP,ISYMT2,WORK,LWORK,LENSQ,INDSQ,
     *                         ISYMIB,IB,ISYMID,ID,W3X)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate the L3 contributions to omega1.
C
      IMPLICIT NONE
C
      LOGICAL W3X
C
      INTEGER ISYMIM, ISYINT, ISYMT2, LWORK, ISYMIB, IB, ISYMID, ID
      INTEGER LENSQ, INDSQ(LENSQ,6)
      INTEGER ISYMBD, ISCKIJ, ISYCKM, LENGTH, ISYTMP, KSCR1
      INTEGER KEND1, LWRK1, KEND2, LWRK2, ISYMCK, ISYMIJ, ISYMDM
      INTEGER ISYMM, KOFF1, KOFF2, KOFF3, NTOTIJ, NTOTCK
      INTEGER KT2TMP, ISYMCM, ISYMK, ISYMC
      INTEGER ISYRES, ISYMI, ISYOOO, NTOIJK, NTOTB, NTOTI, NBI
      INTEGER ISYVVV, ISYEIJ, ISYMKM, ISYME, KSCR2, ISYMEK
      INTEGER ISYMCE, ISYMAC, ISYMA, NTOTCE, NTOTA, ISYENI, ISYMEN
      INTEGER ISYDLM, ISYMN, NTODLM, NTOTE, ISYMDN, ISYDNI
      INTEGER NTOTEN, ISYENF, ISYELM, ISYMLM, ISYML, ISYMFN, ISYFNI
      INTEGER ISYMEL, ISYLMI, NTOTFN, NTOTLM, ISYMEI, KSCR3, ISYMF
      INTEGER ISYMFI, ISYMDL, ISYMIN, NTOTDL, NTOTIN, ISYMD, ISYTMP2
      INTEGER ISYMMN, ISYAMN, ISYDMN, NTOTMN, ISYBMN, ISYMBN, ISYMDI
      INTEGER KMIMAT, KRMAT, KSORT, ISWMAT, ISYAON, ISAONM, KOFFTM
      INTEGER KOFFTP, KOFFRE, NL, NAON, ISONM, KOFFI, KOFFAI, NONM, NA 
      INTEGER ISRMAT, ISYMBA, ISYMDLM, NTOTDLM, ISYMT2W, ISYMBAE
      INTEGER ISYMAE, KT2W
      INTEGER ISYMMLE, ISYMENI, ISYMLNI, KGD, KT2B, KMLNI, ISYMDNI
      INTEGER ISYMNI, ISYMMLN, ISYMML, NTOTML, NTOTN, NTOTMLN
      INTEGER ISYMEML, ISYMEM, ISYMBIN, ISYMBI
      INTEGER ISYMBLM, ISYFINM, ISYMMI, ISYBMI, ISYMBM, ISYMFIN
      INTEGER NTOTFIN, NTOTFNK, NTOTL, IOPT, ISYMFMN, NTOTFMN, ISYMB
      INTEGER ISYFNIM, ISYMFNI, NTOTFNI, ISYMFNM, NTOTFNM, KSCR4
C
#if defined (SYS_CRAY)
      REAL ZERO, ONE
      REAL OMEGA1(*), WMAT(*), TMAT(*)
      REAL XOOOO(*), XOVVO(*), XOOVV(*), XVVVV(*)
      REAL T2TP(*), WORK(LWORK)
      REAL ZERO, ONE, DDOT, XNORM
#else
      DOUBLE PRECISION OMEGA1(*), WMAT(*), TMAT(*)
      DOUBLE PRECISION XOOOO(*), XOVVO(*), XOOVV(*), XVVVV(*)
      DOUBLE PRECISION T2TP(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, DDOT, XNORM
#endif
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO= 0.0D0, ONE= 1.0D0)
C
C If W3X = .TRUE., then WMAT contains wbar_X and all terms in 
C 1-6 are included
C else
C we assume that we get tbar_0 in as WMAT and calculate only
C those terms in 1-6 which are marked with star (*).
C
C Note Wbar_X and Tbar_0 is stored in the same way for fixed IB and ID
C
C
C CALCULATE THE TERMS 1-6 USING W INTERMEDIATE
C
C 1.
C     + L^{dae}_{lno} t^{de}_{lm} g_{inmo}
C 2.
C     + L^{dfg}_{lim} t^{de}_{lm} g_{fage}
C 3.
C     - L^{daf}_{lmn} t^{de}_{lm} g_{iefn}
C 4.
C     - L^{daf}_{lnm} t^{de}_{lm} g_{infe}
C 5.
C     - L^{def}_{lin} t^{de}_{lm} g_{mafn}
C 6.
C     - L^{def}_{lni} t^{de}_{lm} g_{mnfa}
C
C 1.  Calculate contribution from g_{oooo}
C ( Wae(dlon)(*) + Wad(eoln) + Wed(anlo) )* t(dl,em) * g(inmo)
C 
C 2.  Calculate contribution from g_{vvvv}
C ( Wdg(fiml) + Wdf(gmil) + Wfg(dlmi)(*) )* t(dl,em) * g(fage)
C 
C 3.  Calculate contributions from g_{ovvo}
C -( Waf(dlnm)(*) + Wad(fnlm) + Wdf(amnl) )* t(dl,em) * g(iefn)
C
C 4.  Calculate contributions from g_{oovv}
C -( Waf(dlmn)(*) + Wad(fmln) + Wfd(anlm) )* t(dl,em) * g(infe)
C
C 5.  Calculate contributions from g_{ovvo}
C -( Wfe(dlin)(*) + Wfd(eiln) + Wde(fnil) )* t(dl,em) * g(mafn)
C
C 6.  Calculate contributions from g_{oovv}
C -( Wfe(dlni)(*) + Wfd(enli) + Wde(finl) )* t(dl,em) * g(mnfa)
C
      CALL QENTER('CC3_W3_OMEGA1')
C
      ISYTMP  = MULD2H(ISYMIM,ISYMT2)
      ISYTMP2 = MULD2H(ISYTMP,ISYINT)
      IF (ISYRES .NE. ISYTMP2) THEN
         CALL QUIT('Symmetry mimatch in CC3_W3_OMEGA1')
      ENDIF
C
      ISYMBD = MULD2H(ISYMIB,ISYMID)
      ISCKIJ = MULD2H(ISYMIM,ISYMBD)
      ISWMAT = ISCKIJ
C
      LENGTH = NCKIJ(ISCKIJ)
      IF (LENSQ.NE.LENGTH) THEN
         WRITE(LUPRI,*)'CC3_W3_OMEGA1, SYMMETRY MISMATCH : LENGTH
     *   NE LENSQ'
         STOP ' CC3_W3_OMEGA1, LENGTH NE LENSQ'
      END IF
C
C================================================
C 1.    Calculate contribution from g_{oooo}
C     + L^{dae}_{lno} t^{de}_{lm} g_{inmo} =
C ( Wae(dlon) + Wad(eoln) + Wed(anlo) )* t(dl,em) * g(inmo)
C================================================
C
      ISYCKM = MULD2H(ISYMT2,ISYMID)
      ISYTMP = MULD2H(ISCKIJ,ISYCKM) ! Symmetry of intermediate
C
      KSCR1 = 1
      KEND1 = KSCR1 + NMAIJK(ISYTMP)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_W3_OMEGA1')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NMAIJK(ISYTMP))
C
C--------------------------------------------
C     First contribution to intermediate
C--------------------------------------------
C
C
C  Wae(dlon) * t(dl,em) * g(inmo)
C
C  WAD(ckij) * t2tp(ckmD) * g(Ijmi)
C     GAD(ij,m)           
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(I)
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-1)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      DO ISYMCK = 1, NSYM
C
         ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
         ISYMM  = MULD2H(ISYCKM,ISYMCK)
C
         KOFF1  = ISAIKL(ISYMCK,ISYMIJ)
     *          + 1
         KOFF2  = IT2SP(ISYCKM,ISYMID)
     *          + NCKI(ISYCKM)*(ID-1)
     *          + ICKI(ISYMCK,ISYMM)
     *          + 1
         KOFF3  = KSCR1
     *          + IMAIJK(ISYMIJ,ISYMM)
C
         NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
         NTOTCK = MAX(1,NT1AM(ISYMCK))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIJ),NRHF(ISYMM),NT1AM(ISYMCK),
     *              ONE,TMAT(KOFF1),NTOTCK,T2TP(KOFF2),NTOTCK,
     *              ONE,WORK(KOFF3),NTOTIJ)
C
      ENDDO
c
*     write(lupri,*)'kscr1(on,m) in cc3_w3_lhtr ib,id ',ib,id
*     call output(work(kscr1),1,nmatij(1),1,nrhf(1),
*    *            nmatij(1),nrhf(1),
*    *            1,lupri)
C
C--------------------------------------------
C     Second contribution to intermediate
C--------------------------------------------
C
C    Wad(eoln)  t(dl,em) * g(inmo)
C
C    WAD(ckij * tD(ckm)  * g(Ijmi)
C     GAD(ij,m)          
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,1))
         ENDDO
C
         IF (NSYM .GT. 1) THEN
            IF (LWRK1 .LT. LENGTH) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-2)')
            ENDIF
            CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
         ENDIF
C
         KT2TMP = KEND1
         KEND2  = KT2TMP + NCKI(ISYCKM)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Memory exceeded in CC3_W3_OMEGA1')
         ENDIF
C
C        sort t2tp(ckmD) as TD(cmk)
C
         DO ISYMK = 1, NSYM
            ISYMCM = MULD2H(ISYCKM,ISYMK)
            DO ISYMC = 1, NSYM
               ISYMM  = MULD2H(ISYMCM,ISYMC)
               ISYMCK = MULD2H(ISYMC,ISYMK)
C
               DO K = 1, NRHF(ISYMK)
                  DO M = 1, NRHF(ISYMM)
C
                     KOFF1 = IT2SP(ISYCKM,ISYMID)
     *                     + NCKI(ISYCKM)*(ID-1)
     *                     + ICKI(ISYMCK,ISYMM)
     *                     + NT1AM(ISYMCK)*(M-1)
     *                     + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
     *                     + 1
                     KOFF2 = KT2TMP
     *                     + ICKI(ISYMCM,ISYMK)
     *                     + NT1AM(ISYMCM)*(K-1)
     *                     + IT1AM(ISYMC,ISYMM)
     *                     + NVIR(ISYMC)*(M-1)
C
                     CALL DCOPY(NVIR(ISYMC),T2TP(KOFF1),1,WORK(KOFF2),1)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         DO ISYMCK = 1, NSYM
C
            ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
            ISYMM  = MULD2H(ISYCKM,ISYMCK)
C
            KOFF1  = ISAIKL(ISYMCK,ISYMIJ)
     *             + 1
            KOFF2  = KT2TMP
     *             + ICKI(ISYMCK,ISYMM)
            KOFF3  = KSCR1
     *             + IMAIJK(ISYMIJ,ISYMM)
C
            NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
            NTOTCK = MAX(1,NT1AM(ISYMCK))
C
            CALL DGEMM('T','N',NMATIJ(ISYMIJ),NRHF(ISYMM),NT1AM(ISYMCK),
     *                 ONE,TMAT(KOFF1),NTOTCK,WORK(KOFF2),NTOTCK,
     *                 ONE,WORK(KOFF3),NTOTIJ)
C
         ENDDO
C
      END IF
c
*     write(lupri,*)'kscr1(on,m) in cc3_w3_lhtr ib,id ',ib,id
*     call output(work(kscr1),1,nmatij(1),1,nrhf(1),
*    *            nmatij(1),nrhf(1),
*    *            1,lupri)
C
C------------------------------------------------
C     Contract the intermediate with g_{oooo}
C------------------------------------------------
C
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYOOO = MULD2H(ISYINT,ISYMI)
      DO I = 1, NRHF(ISYMI)
         NBI = IT1AM(ISYMIB,ISYMI) + NVIR(ISYMIB)*(I-1) + IB
         KOFF1 = I3ORHF(ISYOOO,ISYMI)
     *         + NMAIJK(ISYOOO)*(I-1)
     *         + 1
C
C
C 1.1 (+ 1.2) omega contribution addomega
C
         OMEGA1(NBI) = OMEGA1(NBI) 
     *               - DDOT(NMAIJK(ISYOOO),XOOOO(KOFF1),1,WORK(KSCR1),1)
      ENDDO
C
C-----------------------------------------
C  Wed(anlo) * t(dl,em) * g(inmo)
C
C  WBD(anlo) * t2t2(DlmB) * g(inmo)
C
C  TMAT(aonl) * T(lm) * G(onmi)
C       
C       R(aonm)
C-----------------------------------------
C
C  TMAT(aonl) = WBD(anlo)
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,2))
         END DO
C
C -------------------------------
C        Sort T^{DB}_{lm} as T_{lm}
C -------------------------------
C
         ISYMBD = MULD2H(ISYMIB,ISYMID)
         ISYMLM = MULD2H(ISYMBD,ISYMT2)
         ISYDLM = MULD2H(ISYMLM,ISYMID)
         ISRMAT = MULD2H(ISWMAT,ISYMLM)
C   
         KMIMAT  = 1 
         KRMAT   = KMIMAT + NMATIJ(ISYMLM)
         KSORT   = KRMAT  + N3VOOO(ISRMAT)
         KEND1   = KSORT  + N3VOOO(ISRMAT)
         LWRK1   = LWORK  - KEND1
C
         IF (LWRK1 .LT.0 ) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (sort)')
         ENDIF
C
         DO ISYMM = 1,NSYM
            ISYML = MULD2H(ISYMLM,ISYMM)
            ISYMDL = MULD2H(ISYDLM,ISYMM)
            DO M = 1,NRHF(ISYMM)
               DO L = 1,NRHF(ISYML)
                  KOFF1 = IT2SP(ISYDLM,ISYMIB)
     *                        + NCKI(ISYDLM)*(IB-1)
     *                        + ICKI(ISYMDL,ISYMM)
     *                        + NT1AM(ISYMDL)*(M-1)
     *                        + IT1AM(ISYMID,ISYML)
     *                        + NVIR(ISYMID)*(L-1)
     *                        + ID
                  KOFF2 = IMATIJ(ISYML,ISYMM)
     *                        + NRHF(ISYML)*(M-1)
     *                        + L
                  WORK(KMIMAT-1+KOFF2) = T2TP(KOFF1)
               ENDDO
            ENDDO
         ENDDO
C
C        R(aonm) = sum_l ( TMAT(aonl) * T(lm) )
C
         DO ISYMM = 1,NSYM
            ISYML   = MULD2H(ISYMLM,ISYMM)
            ISYAON  = MULD2H(ISWMAT,ISYML)
            ISAONM  = MULD2H(ISYAON,ISYMM)
            KOFFTM  = ISAIKJ(ISYAON,ISYML) + 1
            KOFFTP  = IMATIJ(ISYML,ISYMM)  + KMIMAT
            KOFFRE  = ISAIKJ(ISYAON,ISYMM) + KRMAT
            NL      =  MAX(1,NRHF(ISYML))
            NAON    =  MAX(1,NCKI(ISYAON))

C
            CALL DGEMM('N','N',NCKI(ISYAON),NRHF(ISYMM),
     *                 NRHF(ISYML),ONE,TMAT(KOFFTM),NAON,
     *                 WORK(KOFFTP),NL,ZERO,WORK(KOFFRE),NAON) 
C
         END DO
*     write(lupri,*)'krmat(aon,m) in cc3_w3_lhtr ib,id ',ib,id
*     call output(work(krmat),1,ncki(1),1,nrhf(1),
*    *            ncki(1),nrhf(1),
*    *            1,lupri)

C
C        omega(ai) = sum_onm ( R(aonm) * G(onmi) )
C       
C
         IF (NSYM .GT. 1) THEN
           CALL CC3_SRTVOOO(WORK(KSORT),WORK(KRMAT),ISRMAT)
           CALL DCOPY(N3VOOO(ISRMAT),WORK(KSORT),1,WORK(KRMAT),1)
         ENDIF
C

         DO ISYMI = 1,NSYM
            ISYMA  = MULD2H(ISYRES,ISYMI)
            ISONM  = MULD2H(ISYINT,ISYMI)
            KOFFRE = I3VOOO(ISYMA,ISONM) + KRMAT 
            KOFFI  = I3ORHF(ISONM,ISYMI) + 1
            KOFFAI = IT1AM(ISYMA,ISYMI)  + 1
            NONM   = MAX(1,NMAIJK(ISONM))
            NA  = MAX(1,NVIR(ISYMA))

C
C           1.3 omega contribution addomega
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                  NMAIJK(ISONM),-ONE,WORK(KOFFRE),NA,
     *                  XOOOO(KOFFI),NONM,ONE,OMEGA1(KOFFAI),NA)
         END DO
C
      END IF
C
C=============================================
C 2.     Calculate contribution from g_{vvvv}
C     + L^{dfg}_{lim} t^{de}_{lm} g_{fage} =
C ( Wdg(fiml) + Wdf(gmil) + Wfg(dlmi) )* t(dl,em) * g(fage)
C=============================================
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C
         ISYCKM = MULD2H(ISYMT2,ISYMIB)
         ISYEIJ = ISYCKM
         ISYTMP = MULD2H(ISCKIJ,ISYEIJ)
C
         KSCR1  = 1
         KT2TMP = KSCR1  + NCKATR(ISYTMP)
         KEND1  = KT2TMP + NCKI(ISYEIJ)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (VVVV T2-sort)')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),NCKATR(ISYTMP))
C
C----------------
C        Sort T2   tB(km,c) = t2tp(ckmB) 
C----------------
C
         DO ISYMK = 1, NSYM
            ISYMCM = MULD2H(ISYCKM,ISYMK)
            DO ISYMC = 1, NSYM
               ISYMM = MULD2H(ISYMCM,ISYMC)
               ISYMKM = MULD2H(ISYMK,ISYMM)
               ISYMCK = MULD2H(ISYMK,ISYMC)
C
               DO K = 1, NRHF(ISYMK)
                  DO M = 1, NRHF(ISYMM)
                     KOFF1 = IT2SP(ISYCKM,ISYMIB)
     *                     + NCKI(ISYCKM)*(IB-1)
     *                     + ICKI(ISYMCK,ISYMM)
     *                     + NT1AM(ISYMCK)*(M-1)
     *                     + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
     *                     + 1
                     KOFF2 = KT2TMP - 1
     *                     + IMAIJA(ISYMKM,ISYMC)
     *                     + IMATIJ(ISYMK,ISYMM)
     *                     + NRHF(ISYMK)*(M-1)
     *                     + K
C
                     CALL DCOPY(NVIR(ISYMC),T2TP(KOFF1),1,
     *                          WORK(KOFF2),NMATIJ(ISYMKM))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C--------------------------------
C        First intermediate
C--------------------------------
C
C  Wdg(fiml) * t(dl,em) * g(fage)
C
C  WBD(ck,ij) * tB(ij,e) * g(caDe)      
C  
C  TMAT(ck,ij) * tB(ij,e) * gD(cae)
C       G(ck,e) * GD(ce,a)
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(I)
         ENDDO
C
         IF (NSYM .GT. 1) THEN
            IF (LWRK1 .LT. LENGTH) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-3)')
            ENDIF
            CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
         ENDIF
C
         DO ISYMCK = 1, NSYM
            ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
            ISYME  = MULD2H(ISYEIJ,ISYMIJ)
C
            KOFF1 = ISAIKL(ISYMCK,ISYMIJ)
     *            + 1
            KOFF2 = KT2TMP
     *            + IMAIJA(ISYMIJ,ISYME)
            KOFF3 = KSCR1
     *            + ICKATR(ISYMCK,ISYME)
C
            NTOTCK = MAX(1,NT1AM(ISYMCK))
            NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
C
            CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYME),
     *                 NMATIJ(ISYMIJ),ONE,TMAT(KOFF1),NTOTCK,
     *                 WORK(KOFF2),NTOTIJ,ONE,WORK(KOFF3),
     *                 NTOTCK)
         ENDDO
C
C---------------------
C        Sort result.
C---------------------
C
         KSCR2  = KEND1
         KEND2  = KSCR2  + NCKATR(ISYTMP)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (Sorting-1)')
         ENDIF
C
         DO ISYMC = 1, NSYM
            ISYMEK = MULD2H(ISYMC,ISYTMP)
            DO ISYMK = 1, NSYM
               ISYME  = MULD2H(ISYMK,ISYMEK)
               ISYMCE = MULD2H(ISYMC,ISYME)
               ISYMCK = MULD2H(ISYMC,ISYMK)
C
               DO K = 1, NRHF(ISYMK)
                  DO E = 1, NVIR(ISYME)
C
                     KOFF1 = KSCR1
     *                     + ICKATR(ISYMCK,ISYME)
     *                     + NT1AM(ISYMCK)*(E-1)
     *                     + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
                     KOFF2 = KSCR2
     *                     + ICKASR(ISYMCE,ISYMK)
     *                     + NMATAB(ISYMCE)*(K-1)
     *                     + IMATAB(ISYMC,ISYME)
     *                     + NVIR(ISYMC)*(E-1)
C
                     CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),1,
     *                          WORK(KOFF2),1)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         CALL DCOPY(NCKATR(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
C----------------------------------------
C        Sort and contract with integral.
C----------------------------------------
C
         ISYVVV = MULD2H(ISYINT,ISYMID)
C
         KSCR2 = KEND1
         KEND2 = KSCR2 + NMAABC(ISYVVV)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (Sort/Contract)')
         ENDIF
C
         DO ISYME = 1, NSYM
            ISYMAC = MULD2H(ISYVVV,ISYME)
            DO ISYMC = 1, NSYM
               ISYMA  = MULD2H(ISYMAC,ISYMC)
               ISYMCE = MULD2H(ISYMC,ISYME)
C
               DO A = 1, NVIR(ISYMA)
                  DO E = 1, NVIR(ISYME)
C
                     KOFF1 = IMAABC(ISYMAC,ISYME)
     *                     + NMATAB(ISYMAC)*(E-1)
     *                     + IMATAB(ISYMC,ISYMA)
     *                     + NVIR(ISYMC)*(A-1)
     *                     + 1
                     KOFF2 = KSCR2
     *                     + IMAABC(ISYMCE,ISYMA)
     *                     + NMATAB(ISYMCE)*(A-1)
     *                     + IMATAB(ISYMC,ISYME)
     *                     + NVIR(ISYMC)*(E-1)
C
                     CALL DCOPY(NVIR(ISYMC),XVVVV(KOFF1),1,
     *                          WORK(KOFF2),1)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         DO ISYMA = 1, NSYM
            ISYMI = MULD2H(ISYMA,ISYRES)
            ISYMCE = MULD2H(ISYMA,ISYVVV)
C
            KOFF1 = KSCR2
     *            + IMAABC(ISYMCE,ISYMA)
            KOFF2 = KSCR1
     *            + ICKASR(ISYMCE,ISYMI)
            KOFF3 = IT1AM(ISYMA,ISYMI)
     *            + 1
C
            NTOTCE = MAX(1,NMATAB(ISYMCE))
            NTOTA  = MAX(1,NVIR(ISYMA))
C
C
C
C        2.1 omega contribution addomega
C
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATAB(ISYMCE),
     *                 -ONE,WORK(KOFF1),NTOTCE,WORK(KOFF2),NTOTCE,
     *                 ONE,OMEGA1(KOFF3),NTOTA)
         ENDDO
C
C---------------------------------------
C        Second contribution.
C---------------------------------------
C
C  Wdf(gmil) * t(dl,em) * g(fage)
C
C  WBD(gmil) * t(emlB) * g(geDa)
C
C  TMAT(giml) * tB(mle) * gD(gea)
C      
C  TMAT(ckij) * tB(ij,e) * gD(ce,a)
C           G(ck,e)
C  sort     G(ci,e) as M(ce,i)
C 
C  eta(ai) = eta(ai) + gD(ce,a) * M(ce,i) 
C
         ISYCKM = MULD2H(ISYMT2,ISYMIB)
         ISYEIJ = ISYCKM
         ISYTMP = MULD2H(ISCKIJ,ISYEIJ)
C
         KSCR1  = 1
         KT2TMP = KSCR1  + NCKATR(ISYTMP)
         KEND1  = KT2TMP + NCKI(ISYCKM)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (T2-sort-2)')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),NCKATR(ISYTMP))
C
C----------------
C        Sort T2   tB(km,c) = t2tp(ckmB)
C----------------
C
         DO ISYMK = 1, NSYM
            ISYMCM = MULD2H(ISYCKM,ISYMK)
            DO ISYMC = 1, NSYM
               ISYMM = MULD2H(ISYMCM,ISYMC)
               ISYMKM = MULD2H(ISYMK,ISYMM)
               ISYMCK = MULD2H(ISYMK,ISYMC)
C
               DO K = 1, NRHF(ISYMK)
                  DO M = 1, NRHF(ISYMM)
                     KOFF1 = IT2SP(ISYCKM,ISYMIB)
     *                     + NCKI(ISYCKM)*(IB-1)
     *                     + ICKI(ISYMCK,ISYMM)
     *                     + NT1AM(ISYMCK)*(M-1)
     *                     + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
     *                     + 1
                     KOFF2 = KT2TMP - 1
     *                     + IMAIJA(ISYMKM,ISYMC)
     *                     + IMATIJ(ISYMK,ISYMM)
     *                     + NRHF(ISYMK)*(M-1)
     *                     + K
C
                     CALL DCOPY(NVIR(ISYMC),T2TP(KOFF1),1,
     *                          WORK(KOFF2),NMATIJ(ISYMKM))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C--------------------------------
C     Second intermediate
C--------------------------------
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,1))
         ENDDO
C
         IF (NSYM .GT. 1) THEN
            IF (LWRK1 .LT. LENGTH) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-3)')
            ENDIF
            CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
         ENDIF
C
         DO ISYMCK = 1, NSYM
            ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
            ISYME  = MULD2H(ISYEIJ,ISYMIJ)
C
            KOFF1 = ISAIKL(ISYMCK,ISYMIJ)
     *         + 1
            KOFF2 = KT2TMP
     *         + IMAIJA(ISYMIJ,ISYME)
            KOFF3 = KSCR1
     *         + ICKATR(ISYMCK,ISYME)
C
            NTOTCK = MAX(1,NT1AM(ISYMCK))
            NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
C
            CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYME),
     *              NMATIJ(ISYMIJ),ONE,TMAT(KOFF1),NTOTCK,
     *              WORK(KOFF2),NTOTIJ,ONE,WORK(KOFF3),
     *              NTOTCK)
         ENDDO
C
C---------------------
C        Sort result.
C---------------------
C
         KSCR2  = KEND1
         KEND2  = KSCR2  + NCKATR(ISYTMP)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (Sorting-1)')
         ENDIF
C
         DO ISYMC = 1, NSYM
            ISYMEK = MULD2H(ISYMC,ISYTMP)
            DO ISYMK = 1, NSYM
               ISYME  = MULD2H(ISYMK,ISYMEK)
               ISYMCE = MULD2H(ISYMC,ISYME)
               ISYMCK = MULD2H(ISYMC,ISYMK)
C
               DO K = 1, NRHF(ISYMK)
                  DO E = 1, NVIR(ISYME)
C
                     KOFF1 = KSCR1 - 1
     *                     + ICKATR(ISYMCK,ISYME)
     *                     + NT1AM(ISYMCK)*(E-1)
     *                     + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
     *                     + 1
                     KOFF2 = KSCR2 - 1
     *                     + ICKASR(ISYMCE,ISYMK)
     *                     + NMATAB(ISYMCE)*(K-1)
     *                     + IMATAB(ISYMC,ISYME)
     *                     + NVIR(ISYMC)*(E-1)
     *                     + 1
C
                     CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),1,
     *                          WORK(KOFF2),1)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         CALL DCOPY(NCKATR(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
C----------------------------------------
C        Contract with integral.
C----------------------------------------
C
         ISYVVV = MULD2H(ISYINT,ISYMID)
C
         DO ISYMA = 1, NSYM
            ISYMI = MULD2H(ISYMA,ISYRES)
            ISYMCE = MULD2H(ISYMA,ISYVVV)
C
            KOFF1 = IMAABC(ISYMCE,ISYMA)
     *            + 1
            KOFF2 = KSCR1
     *            + ICKASR(ISYMCE,ISYMI)
            KOFF3 = IT1AM(ISYMA,ISYMI)
     *            + 1
C
            NTOTCE = MAX(1,NMATAB(ISYMCE))
            NTOTA  = MAX(1,NVIR(ISYMA))
C
C        2.2 omega contribution addomega
C
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATAB(ISYMCE),
     *                 -ONE,XVVVV(KOFF1),NTOTCE,WORK(KOFF2),NTOTCE,
     *                 ONE,OMEGA1(KOFF3),NTOTA)
         ENDDO
C
      END IF
C
C third contribution
C
C-----------------------------------------
C  Wfg(dlmi) * t(dl,em) * g(fage)
C
C  WBD(dlmi) * t2tp(dlme) * g(BaDe)
C
C   t2tp(dlme) *  WBD(dlmi) * gD(Bae)
C
C   eta (ai) = eta(ai) + gDB(ae) * t2w(ei)  
C----------------------------------------
C
C sort gD(Bae) as gDB(ae)
C
      ISWMAT = ISCKIJ
      ISYMT2W = MULD2H(ISWMAT,ISYMT2)
      ISYMBAE = MULD2H(ISYINT,ISYMID)
      ISYMAE = MULD2H(ISYMBAE,ISYMIB)
C
      KSCR2 = 1
      KT2W  = KSCR2 + NMATAB(ISYMAE)
      KEND2 = KT2W  + NT1AM(ISYMT2W)
      LWRK2 = LWORK - KEND2
C
      CALL DZERO(WORK(KT2W),NT1AM(ISYMT2W))
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_W3_OMEGA1 (Sort/Contract)')
      ENDIF
C
      DO ISYME = 1, NSYM
         ISYMBA = MULD2H(ISYMBAE,ISYME)
         ISYMA  = MULD2H(ISYMBA,ISYMIB)
C
         DO E = 1, NVIR(ISYME)
            DO A = 1, NVIR(ISYMA)
C
               KOFF1 = IMAABC(ISYMBA,ISYME)
     *               + NMATAB(ISYMBA)*(E-1)
     *               + IMATAB(ISYMIB,ISYMA)
     *               + NVIR(ISYMIB)*(A-1)
     *               + IB
               KOFF2 = KSCR2
     *               + IMATAB(ISYMA,ISYME)
     *               + NVIR(ISYMA)*(E-1)
     *               + A -1
C
               WORK(KOFF2) = XVVVV(KOFF1) 
C
            ENDDO
         ENDDO
      ENDDO
C
C  t2w(ei) = t2tp(dlme) *  WBD(dlmi) 
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(I)
      ENDDO
C
      DO ISYMI = 1,NSYM
         ISYMDLM = MULD2H(ISWMAT,ISYMI)
         ISYME   = MULD2H(ISYMDLM,ISYMT2)
C
         KOFF1 = IT2SP(ISYMDLM,ISYME) + 1
         KOFF2 = ISAIKJ(ISYMDLM,ISYMI) + 1
         KOFF3 = KT2W + IT1AM(ISYME,ISYMI)
C
         NTOTDLM = MAX(1,NCKI(ISYMDLM))
         NTOTE   = MAX(1,NVIR(ISYME))
C
         CALL DGEMM('T','N',NVIR(ISYME),NRHF(ISYMI),
     *              NCKI(ISYMDLM),ONE,T2TP(KOFF1),NTOTDLM,
     *              TMAT(KOFF2),NTOTDLM,ONE,WORK(KOFF3),
     *              NTOTE)

C
      END DO
C
C   eta (ai) = eta(ai) + gDB(ae) * t2w(ei)  
C
      DO ISYMI = 1,NSYM
         ISYME = MULD2H(ISYMT2W,ISYMI)
         ISYMA   = MULD2H(ISYMAE,ISYME)
C
         KOFF1 = KSCR2 + IMATAB(ISYMA,ISYME)
         KOFF2 = KT2W  + IT1AM(ISYME,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA = MAX(1,NVIR(ISYMA))
         NTOTE = MAX(1,NVIR(ISYME))
C
C 2.3 omega contribution addomega
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *              NVIR(ISYME),-ONE,WORK(KOFF1),NTOTA,
     *              WORK(KOFF2),NTOTE,ONE,OMEGA1(KOFF3),
     *              NTOTA)
      END DO
C
C================================================
C 3.  Calculate contributions from g_{ovvo}
C     - L^{daf}_{lmn} t^{de}_{lm} g_{iefn} =
C -( Waf(dlnm) + Wad(fnlm) + Wdf(amnl) )* t(dl,em) * g(iefn)
C
C================================================
C
C--------------------------------
C     First contribution
C--------------------------------
C - Waf(dlnm) * t(dl,em) * g(iefn)
C
C   t2tp(dlme) * Tmat(dlmn)   xovvo(Dnie)
C
C  omega(ib,i) = omega(ib,i) +  t2T(en) * gD(eni)
C
      ISYMEN = MULD2H(ISCKIJ,ISYMT2)
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYENI = MULD2H(ISYMEN,ISYMI)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NT1AM(ISYMEN)
      KEND1 = KSCR2 + NCKI(ISYENI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_W3_OMEGA1 (OVVO-1)')
      ENDIF
C
C t2T(en) =   t2tp(dlme) * Tmat(dlmn) 
C
      CALL DZERO(WORK(KSCR1),NT1AM(ISYMEN))
      CALL DZERO(WORK(KSCR2),NCKI(ISYENI))
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(INDSQ(I,3))
      ENDDO
C
      DO ISYME = 1, NSYM
         ISYDLM = MULD2H(ISYME,ISYMT2)
         ISYMN  = MULD2H(ISYMEN,ISYME)
C
         KOFF1 = IT2SP(ISYDLM,ISYME)
     *         + 1
         KOFF2 = ISAIKJ(ISYDLM,ISYMN)
     *         + 1
         KOFF3 = KSCR1
     *         + IT1AM(ISYME,ISYMN)
C
         NTODLM = MAX(1,NCKI(ISYDLM))
         NTOTE  = MAX(1,NVIR(ISYME))
C
         CALL DGEMM('T','N',NVIR(ISYME),NRHF(ISYMN),NCKI(ISYDLM),
     *              -ONE,T2TP(KOFF1),NTODLM,TMAT(KOFF2),NTODLM,
     *              ONE,WORK(KOFF3),NTOTE)
C
      ENDDO
C
C sort integrals xovvo(Dnie) as gD(eni) 
C
      DO ISYME = 1, NSYM
         ISYMN  = MULD2H(ISYMEN,ISYME)
         ISYMDN = MULD2H(ISYMN,ISYMID)
         ISYDNI = MULD2H(ISYMDN,ISYMI)
C
         DO E = 1, NVIR(ISYME)
            DO N = 1, NRHF(ISYMN)
C
               KOFF1 = IT2SP(ISYDNI,ISYME)
     *               + NCKI(ISYDNI)*(E-1)
     *               + ICKI(ISYMDN,ISYMI)
     *               + IT1AM(ISYMID,ISYMN)
     *               + NVIR(ISYMID)*(N-1)
     *               + ID
C
               KOFF2 = KSCR2 - 1
     *               + ICKI(ISYMEN,ISYMI)
     *               + IT1AM(ISYME,ISYMN)
     *               + NVIR(ISYME)*(N-1)
     *               + E
C
               CALL DCOPY(NRHF(ISYMI),XOVVO(KOFF1),NT1AM(ISYMDN),
     *                    WORK(KOFF2),NT1AM(ISYMEN))
C
            ENDDO
         ENDDO
      ENDDO
C
C  omega(ib,i) = omega(ib,i) +  t2T(en) * gD(eni)
C
      KOFF1  = KSCR2
     *       + ICKI(ISYMEN,ISYMI)
      KOFF3  = IT1AM(ISYMIB,ISYMI)
     *       + IB
C
      NTOTEN = MAX(1,NT1AM(ISYMEN))
      NTOTB  = MAX(1,NVIR(ISYMIB))
C
C 3.1 omega contribution addomega
C
      CALL DGEMV('T',NT1AM(ISYMEN),NRHF(ISYMI),-ONE,WORK(KOFF1),
     *           NTOTEN,WORK(KSCR1),1,ONE,OMEGA1(KOFF3),NTOTB)
C
C
C--------------------------------
C     Second contribution
C
C - Wad(fnlm) * t(dl,em) * g(iefn)
C
C   TMAT(fnlm) * xovvo(fnie)    t2tp(emlD)
C
C   omega1(Bi) = omega(Bi) - Tx(lmie) * t2D(lme) 
C--------------------------------
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C

         ISYTMP = MULD2H(ISCKIJ,ISYINT)
         ISYMI  = MULD2H(ISYRES,ISYMIB)
         ISYENF = MULD2H(ISYINT,ISYMI)
         ISYELM = MULD2H(ISYMT2,ISYMID)
C
         KSCR1 = 1
         KEND1 = KSCR1 + NMAIJA(ISYELM)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (OVVO-2)')
         ENDIF
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(I)
         ENDDO
C
         IF (NSYM .GT. 1) THEN
            IF (LWORK .LT. LENGTH) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-4)')
            ENDIF
            CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK,1,TMAT,1)
         ENDIF
C
C sort t2 amplitudes  t2tp(elmD) as t2D(mle)
C
         DO ISYME = 1, NSYM
            ISYMLM = MULD2H(ISYELM,ISYME)
            DO ISYML = 1, NSYM
               ISYMM  = MULD2H(ISYMLM,ISYML)
               ISYMEL = MULD2H(ISYME,ISYML)
C
               DO L = 1, NRHF(ISYML)
                  DO M = 1, NRHF(ISYMM)
C
                     KOFF1 = IT2SP(ISYELM,ISYMID)
     *                  + NCKI(ISYELM)*(ID-1)
     *                  + ICKI(ISYMEL,ISYMM)
     *                  + NT1AM(ISYMEL)*(M-1)
     *                  + IT1AM(ISYME,ISYML)
     *                  + NVIR(ISYME)*(L-1)
     *                  + 1
C
                     KOFF2 = KSCR1 - 1
     *                  + IMAIJA(ISYMLM,ISYME)
     *                  + IMATIJ(ISYMM,ISYML)
     *                  + NRHF(ISYMM)*(L-1)
     *                  + M
C
                     CALL DCOPY(NVIR(ISYME),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMLM))
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         DO ISYME = 1, NSYM
C
            ISYMFN = MULD2H(ISYME,ISYENF)
            ISYMLM = MULD2H(ISCKIJ,ISYMFN)
            ISYFNI = MULD2H(ISYMI,ISYMFN)
            ISYLMI = MULD2H(ISYMI,ISYMLM)
C
            KSCR2 = KEND1
            KEND2 = KSCR2 + NMAIJK(ISYLMI)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (OVVO-3)')
            ENDIF
C
            DO E = 1, NVIR(ISYME)
C
               KOFF1 = ISAIKL(ISYMFN,ISYMLM)
     *            + 1
               KOFF2 = IT2SP(ISYFNI,ISYME)
     *            + NCKI(ISYFNI)*(E-1)
     *            + ICKI(ISYMFN,ISYMI)
     *            + 1
               KOFF3 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
C
               NTOTFN = MAX(1,NT1AM(ISYMFN))
               NTOTLM = MAX(1,NMATIJ(ISYMLM))
C
C   TMAT(fnlm) * xovvo(fnie)  
C
C
               CALL DGEMM('T','N',NMATIJ(ISYMLM),NRHF(ISYMI),
     *                  NT1AM(ISYMFN),
     *                 -ONE,TMAT(KOFF1),NTOTFN,XOVVO(KOFF2),NTOTFN,
     *                 ZERO,WORK(KOFF3),NTOTLM)
C
               KOFF1 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
               KOFF2 = KSCR1
     *            + IMAIJA(ISYMLM,ISYME)
     *            + NMATIJ(ISYMLM)*(E-1)
               KOFF3 = IT1AM(ISYMIB,ISYMI)
     *            + IB
C
               NTOTB  = MAX(1,NVIR(ISYMIB))
               NTOTIJ = MAX(1,NMATIJ(ISYMLM))
C
C   omega1(Bi) = omega(Bi) - Tx(lmie) * t2D(lme) 
C
C
C 3.2 omega contribution addomega
C
               CALL DGEMV('T',NMATIJ(ISYMLM),NRHF(ISYMI),-ONE, 
     *                 WORK(KOFF1),
     *                 NTOTIJ,WORK(KOFF2),1,ONE,OMEGA1(KOFF3),NTOTB)
C
C
             ENDDO
         ENDDO
C
C----------------------------------------------
C third contribution
C
C - Wdf(amnl) * t(dl,em) * g(iefn)
C
C - WBD(amnl) * t2tp(emlB) * xovvo(Dnie)
C
C - TMAT(amln)   t2B(mle)  * gD(eni)
C
C omega1(ai) = omega1(ai) - TMAT(amln) * t2g(mlni)
C----------------------------------------------
C
         ISYMMLE = MULD2H(ISYMT2,ISYMIB)
         ISYMENI = MULD2H(ISYINT,ISYMID)
         ISYMLNI = MULD2H(ISYMMLE,ISYMENI)
         KGD   = 1
         KT2B  = KGD  + NCKI(ISYMENI)
         KMLNI = KT2B + NCKI(ISYMMLE) 
         KSCR1 = KMLNI + N3ORHF(ISYMLNI)
         KEND1 = KSCR1 + NCKIJ(ISWMAT)
         LWRK1 = LWORK - KEND1
C
         CALL DZERO(WORK(KMLNI),N3ORHF(ISYMLNI))
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,3))
         ENDDO
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 3.3')
         ENDIF
C
         IF (NSYM .GT. 1) THEN
           CALL CC3_SRTVOOO(WORK(KSCR1),TMAT,ISWMAT)
           CALL DCOPY(LENGTH,WORK(KSCR1),1,TMAT,1)
         END IF
C

C
C----------------
C     Sort T2   tB(kmc) = t2tp(ckmB)
C----------------
C
         ISYCKM = MULD2H(ISYMT2,ISYMIB)
         DO ISYMK = 1, NSYM
            ISYMCM = MULD2H(ISYCKM,ISYMK)
            DO ISYMC = 1, NSYM
               ISYMM = MULD2H(ISYMCM,ISYMC)
               ISYMKM = MULD2H(ISYMK,ISYMM)
               ISYMCK = MULD2H(ISYMK,ISYMC)
C
               DO K = 1, NRHF(ISYMK)
                  DO M = 1, NRHF(ISYMM)
                     KOFF1 = IT2SP(ISYCKM,ISYMIB)
     *                  + NCKI(ISYCKM)*(IB-1)
     *                  + ICKI(ISYMCK,ISYMM)
     *                  + NT1AM(ISYMCK)*(M-1)
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K-1)
     *                  + 1
                     KOFF2 = KT2B - 1
     *                  + IMAIJA(ISYMKM,ISYMC)
     *                  + IMATIJ(ISYMK,ISYMM)
     *                  + NRHF(ISYMK)*(M-1)
     *                  + K
C
                     CALL DCOPY(NVIR(ISYMC),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMKM))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C---------------------------------------
C sort xovvo(Dnie) as gD(eni) 
C---------------------------------------
C 
         DO ISYME = 1,NSYM 
            ISYMDNI = MULD2H(ISYINT,ISYME) 
            ISYMNI  = MULD2H(ISYMDNI,ISYMID)
            DO ISYMI = 1,NSYM 
               ISYMN = MULD2H(ISYMNI,ISYMI) 
               ISYMEN = MULD2H(ISYMENI,ISYMI) 
               ISYMDN = MULD2H(ISYMDNI,ISYMI)
               DO E = 1,NVIR(ISYME)
                  DO I = 1,NRHF(ISYMI)
                     DO N = 1,NRHF(ISYMN) 
                        KOFF1 = IT2SP(ISYMDNI,ISYME)
     *                     + NCKI(ISYMDNI)*(E-1)
     *                     + ICKI(ISYMDN,ISYMI)
     *                     + NT1AM(ISYMDN)*(I-1)
     *                     + IT1AM(ISYMID,ISYMN)
     *                     + NVIR(ISYMID)*(N-1)
     *                     + ID
                        KOFF2 = KGD - 1
     *                     + ICKI(ISYMEN,ISYMI)
     *                     + NT1AM(ISYMEN)*(I-1)
     *                     + IT1AM(ISYME,ISYMN)
     *                     + NVIR(ISYME)*(N-1)
     *                     + E
                        WORK(KOFF2) = XOVVO(KOFF1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C - WBD(amnl) * t2tp(emlB) * xovvo(Dnie)
C
C - TMAT(amln)   t2B(mle)  * gD(eni)
C
         DO ISYMI = 1,NSYM
            ISYMMLN = MULD2H(ISYMLNI,ISYMI)
            ISYMEN = MULD2H(ISYMENI,ISYMI)
            DO ISYMN = 1,NSYM
               ISYME = MULD2H(ISYMEN,ISYMN)
               ISYMML = MULD2H(ISYMMLE,ISYME) 
               DO I = 1,NRHF(ISYMI)
C
                  KOFF1 = KT2B
     *               + IMAIJA(ISYMML,ISYME) 
                  KOFF2 = KGD 
     *               + ICKI(ISYMEN,ISYMI)
     *               + NT1AM(ISYMEN)*(I-1)
     *               + IT1AM(ISYME,ISYMN)

                  KOFF3 = KMLNI
     *               + I3ORHF(ISYMMLN,ISYMI) 
     *               + NMAIJK(ISYMMLN)*(I-1)
     *               + IMAIJK(ISYMML,ISYMN)
C
                  NTOTML = MAX(1,NMATIJ(ISYMML))
                  NTOTE  = MAX(1,NVIR(ISYME))
C          
C
                  CALL DGEMM('N','N',NMATIJ(ISYMML),NRHF(ISYMN),
     *              NVIR(ISYME),
     *              ONE,WORK(KOFF1),NTOTML,WORK(KOFF2),NTOTE,
     *              ONE,WORK(KOFF3),NTOTML)
C

               END DO
            END DO
         END DO
C
C omega1(ai) = omega1(ai) - TMAT(amln) * t2Bgd(mlni)  
C
C                                      t2B(mle) * gD(eni)
C
         DO ISYMI = 1,NSYM
            ISYMMLN = MULD2H(ISYMLNI,ISYMI)
            ISYMA   = MULD2H(ISWMAT,ISYMMLN)
C
            KOFF1 = I3VOOO(ISYMA,ISYMMLN) + 1
            KOFF2 = KMLNI + I3ORHF(ISYMMLN,ISYMI)
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
            NTOTA = MAX(1,NVIR(ISYMA))
            NTOTMLN = MAX(1,NMAIJK(ISYMMLN))
C
C 3.3 omega contribution addomega
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NMAIJK(ISYMMLN),
     *              ONE,TMAT(KOFF1),NTOTA,WORK(KOFF2),NTOTMLN,
     *              ONE,OMEGA1(KOFF3),NTOTA)
         END DO
C
      END IF
C
C================================================
C 4.  Calculate contributions from g_{oovv}
C     - L^{daf}_{lnm} t^{de}_{lm} g_{infe} =
C -( Waf(dlmn) + Wad(fmln) + Wfd(anlm) )* t(dl,em) * g(infe)
C================================================
C
C-------------------------------- 
C     First contribution
C--------------------------------
C - Waf(dlmn) * t(dl,em) * g(infe)
C
C t2T(en) = t2tp(glme) * TMAT(glmn) 
C
C g_{inDe} = xoovv(Dine) = gD(eni)
C  
C wmega1(Ai) = wmega1(Ai) + gD(eni) * t2T(en) 
C
      ISYMEN = MULD2H(ISCKIJ,ISYMT2)
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYENI = MULD2H(ISYMEN,ISYMI)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NT1AM(ISYMEN)
      KEND1 = KSCR2 + NCKI(ISYENI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_W3_OMEGA1 (OOVV-1)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NT1AM(ISYMEN))
      CALL DZERO(WORK(KSCR2),NCKI(ISYENI))
C
C
C t2T(en) = t2tp(glme) * TMAT(glmn) 
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(I)
      ENDDO
C
      DO ISYME = 1, NSYM
         ISYDLM = MULD2H(ISYME,ISYMT2)
         ISYMN  = MULD2H(ISYMEN,ISYME)
C
         KOFF1 = IT2SP(ISYDLM,ISYME)
     *         + 1
         KOFF2 = ISAIKJ(ISYDLM,ISYMN)
     *         + 1
         KOFF3 = KSCR1
     *         + IT1AM(ISYME,ISYMN)
C
         NTODLM = MAX(1,NCKI(ISYDLM))
         NTOTE  = MAX(1,NVIR(ISYME))
C
         CALL DGEMM('T','N',NVIR(ISYME),NRHF(ISYMN),NCKI(ISYDLM),
     *              -ONE,T2TP(KOFF1),NTODLM,TMAT(KOFF2),NTODLM,
     *              ONE,WORK(KOFF3),NTOTE)
C
      ENDDO
C
C g_{inDe} = xoovv(Dine) = gD(eni)
C  
      DO ISYME = 1, NSYM
         ISYMN  = MULD2H(ISYMEN,ISYME)
         ISYMEI = MULD2H(ISYME,ISYMI)
         ISYMDN = MULD2H(ISYMN,ISYMID)
         ISYMDI = MULD2H(ISYMI,ISYMID)
         ISYDNI = MULD2H(ISYMDN,ISYMI)
C
         DO E = 1, NVIR(ISYME)
            DO N = 1, NRHF(ISYMN)
C
               KOFF1 = IT2SP(ISYDNI,ISYME)
     *               + NCKI(ISYDNI)*(E-1)
     *               + ICKI(ISYMDI,ISYMN)
     *               + NT1AM(ISYMDI)*(N-1)
     *               + IT1AM(ISYMID,ISYMI)
     *               + ID
C
               KOFF2 = KSCR2 - 1
     *               + ICKI(ISYMEN,ISYMI)
     *               + IT1AM(ISYME,ISYMN)
     *               + NVIR(ISYME)*(N-1)
     *               + E
CC
               CALL DCOPY(NRHF(ISYMI),XOOVV(KOFF1),NVIR(ISYMID),
     *                    WORK(KOFF2),NT1AM(ISYMEN))
C
            ENDDO
         ENDDO
      ENDDO
C  
C wmega1(Ai) = wmega1(Ai) + gD(eni) * t2T(en) 
C
      KOFF1  = KSCR2
     *       + ICKI(ISYMEN,ISYMI)
      KOFF3  = IT1AM(ISYMIB,ISYMI)
     *       + IB
C
      NTOTEN = MAX(1,NT1AM(ISYMEN))
      NTOTB  = MAX(1,NVIR(ISYMIB))
C
C 4.1 omega contribution addomega
C
      CALL DGEMV('T',NT1AM(ISYMEN),NRHF(ISYMI),-ONE,WORK(KOFF1),
     *           NTOTEN,WORK(KSCR1),1,ONE,OMEGA1(KOFF3),NTOTB)
C
C--------------------------------
C     Second contribution
C--------------------------------
C
C - Wad(fmln) * t(dl,em) * g(infe)
C
C g_{infe} = xoovv(fine) = ge(fni)
C
C t2tp(elmD) = tD(mle)
C
C  TX(lmie) = TMAT(fnlm) * ge(fni)        tD(mle) 
C  
C wmega1(Ai) = wmega1(Ai) -  TX(lmie) *  tD(mle) 
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C

         ISYTMP = MULD2H(ISCKIJ,ISYINT)
         ISYMI  = MULD2H(ISYRES,ISYMIB)
         ISYENF = MULD2H(ISYINT,ISYMI)
         ISYELM = MULD2H(ISYMT2,ISYMID)
C
         KSCR1 = 1
         KEND1 = KSCR1 + NMAIJA(ISYELM)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (OOVV-2)')
         ENDIF
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,5))
         ENDDO
C
         IF (NSYM .GT. 1) THEN
            IF (LWORK .LT. LENGTH) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-5)')
            ENDIF
            CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK,1,TMAT,1)
         ENDIF
C
C t2tp(elmD) = tD(mle)
C
         DO ISYME = 1, NSYM
            ISYMLM = MULD2H(ISYELM,ISYME)
            DO ISYML = 1, NSYM
               ISYMM  = MULD2H(ISYMLM,ISYML)
               ISYMEL = MULD2H(ISYME,ISYML)
C
               DO L = 1, NRHF(ISYML)
                  DO M = 1, NRHF(ISYMM)
C
                     KOFF1 = IT2SP(ISYELM,ISYMID)
     *                  + NCKI(ISYELM)*(ID-1)
     *                  + ICKI(ISYMEL,ISYMM)
     *                  + NT1AM(ISYMEL)*(M-1)
     *                  + IT1AM(ISYME,ISYML)
     *                  + NVIR(ISYME)*(L-1)
     *                  + 1
C
                     KOFF2 = KSCR1 - 1
     *                  + IMAIJA(ISYMLM,ISYME)
     *                  + IMATIJ(ISYMM,ISYML)
     *                  + NRHF(ISYMM)*(L-1)
     *                  + M
C
                     CALL DCOPY(NVIR(ISYME),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMLM))
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         DO ISYME = 1, NSYM
C
            ISYMFN = MULD2H(ISYME,ISYENF)
            ISYMLM = MULD2H(ISCKIJ,ISYMFN)
            ISYFNI = MULD2H(ISYMI,ISYMFN)
            ISYLMI = MULD2H(ISYMI,ISYMLM)
C
            KSCR2 = KEND1
            KSCR3 = KSCR2 + NMAIJK(ISYLMI)
            KEND2 = KSCR3 + NCKI(ISYFNI)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (OOVV-3)')
            ENDIF
C
C g_{infe} = xoovv(fine) = ge(fni)
C
            DO E = 1, NVIR(ISYME)
C
               DO ISYMF = 1, NSYM
                  ISYMN  = MULD2H(ISYMFN,ISYMF)
                  ISYMFI = MULD2H(ISYMI,ISYMF)
                  DO F = 1, NVIR(ISYMF)
                     DO N = 1, NRHF(ISYMN)
C
                        KOFF1 = IT2SP(ISYFNI,ISYME)
     *                     + NCKI(ISYFNI)*(E-1)
     *                     + ICKI(ISYMFI,ISYMN)
     *                     + NT1AM(ISYMFI)*(N-1)
     *                     + IT1AM(ISYMF,ISYMI)
     *                     + F
C
                        KOFF2 = KSCR3 - 1
     *                     + ICKI(ISYMFN,ISYMI)
     *                     + IT1AM(ISYMF,ISYMN)
     *                     + NVIR(ISYMF)*(N-1)
     *                     + F
C
                        CALL DCOPY(NRHF(ISYMI),XOOVV(KOFF1),NVIR(ISYMF),
     *                          WORK(KOFF2),NT1AM(ISYMFN))
                     ENDDO
                  ENDDO
               ENDDO
C
               KOFF1 = ISAIKL(ISYMFN,ISYMLM)
     *            + 1
               KOFF2 = KSCR3
     *            + ICKI(ISYMFN,ISYMI)
               KOFF3 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
C
               NTOTFN = MAX(1,NT1AM(ISYMFN))
               NTOTLM = MAX(1,NMATIJ(ISYMLM))
C
C  TX(lmie) = TMAT(fnlm) * ge(fni)   
C  
               CALL DGEMM('T','N',NMATIJ(ISYMLM),NRHF(ISYMI),
     *                  NT1AM(ISYMFN),
     *                 -ONE,TMAT(KOFF1),NTOTFN,WORK(KOFF2),NTOTFN,
     *                 ZERO,WORK(KOFF3),NTOTLM)
C
               KOFF1 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
               KOFF2 = KSCR1
     *            + IMAIJA(ISYMLM,ISYME)
     *            + NMATIJ(ISYMLM)*(E-1)
               KOFF3 = IT1AM(ISYMIB,ISYMI)
     *            + IB
C
               NTOTB  = MAX(1,NVIR(ISYMIB))
               NTOTIJ = MAX(1,NMATIJ(ISYMLM))
C  
C wmega1(Ai) = wmega1(Ai) -  TX(lmie) *  tD(mle) 
C
C
C 4.2 omega contribution addomega
C
               CALL DGEMV('T',NMATIJ(ISYMLM),NRHF(ISYMI),-ONE,
     *                 WORK(KOFF1),
     *                 NTOTIJ,WORK(KOFF2),1,ONE,OMEGA1(KOFF3),NTOTB)
C

            ENDDO
         ENDDO
C
C--------------------------------
C     Third  contribution
C--------------------------------
C
C - Wfd(anlm) * t(dl,em) * g(infe)
C
C - WBD(anlm) * t2tp(emlD) * xoovv(Bine)
C
C      tD(mle) * xB(eni) = tdxB(mlni)
C
C omega1(ai) = omega1(ai) +  TMAT(amln) * tdxB(mlni)
C
C
         ISYMMLE  = MULD2H(ISYMT2,ISYMID)
         ISYMENI  = MULD2H(ISYINT,ISYMIB)
         ISYMLNI = MULD2H(ISYMMLE,ISYMENI)
C
         KSCR1 = 1
         KSCR2 = KSCR1 + NMAIJA(ISYMMLE)
         KSCR3 = KSCR2 + NMAIJA(ISYMENI) 
         KSCR4 = KSCR3 + N3ORHF(ISYMLNI) 
         KEND1 = KSCR4 + NCKIJ(ISWMAT) 
         LWRK1 = LWORK - KEND1
C
         CALL DZERO(WORK(KSCR3),N3ORHF(ISYMLNI))
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 4.3')
         ENDIF
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,5))
         ENDDO
C
         IF (LWORK .LT. NCKIJ(ISCKIJ)) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 4.3.1')
         ENDIF
C
         IF (NSYM .GT. 1) THEN
           CALL CC3_SRTVOOO(WORK(KSCR4),TMAT,ISWMAT)
           CALL DCOPY(LENGTH,WORK(KSCR4),1,TMAT,1)
         END IF
C
C
C t2tp(emlD) = tD(mle) 
C
         ISYMEML = MULD2H(ISYMT2,ISYMID)
         DO ISYME = 1, NSYM
            ISYMML = MULD2H(ISYMEML,ISYME)
            DO ISYML = 1, NSYM
               ISYMM  = MULD2H(ISYMML,ISYML)
               ISYMEM = MULD2H(ISYME,ISYMM)
C
               DO E = 1,NVIR(ISYME)
                  DO L = 1, NRHF(ISYML)
                     DO M = 1, NRHF(ISYMM)
C
                        KOFF1 = IT2SP(ISYMEML,ISYMID)
     *                     + NCKI(ISYMEML)*(ID-1)
     *                     + ICKI(ISYMEM,ISYML)
     *                     + NT1AM(ISYMEM)*(L-1)
     *                     + IT1AM(ISYME,ISYMM)
     *                     + NVIR(ISYME)*(M-1)
     *                     + E
C
                        KOFF2 = KSCR1 - 1
     *                     + IMAIJA(ISYMML,ISYME)
     *                     + NMATIJ(ISYMML)*(E-1)
     *                     + IMATIJ(ISYMM,ISYML)
     *                     + NRHF(ISYMM)*(L-1)
     *                     + M
C
                        WORK(KOFF2) = T2TP(KOFF1)
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C
C  xoovv(Bine) =  xB(eni)
C
C
         DO ISYME = 1, NSYM
            ISYMBIN = MULD2H(ISYINT,ISYME)
            DO ISYMN = 1,NSYM
               ISYMEN = MULD2H(ISYME,ISYMN)
               ISYMBI = MULD2H(ISYMBIN,ISYMN)
               ISYMI  = MULD2H(ISYMBI,ISYMIB)
               DO E = 1, NVIR(ISYME)
                  DO N = 1, NRHF(ISYMN)
                     DO I = 1, NRHF(ISYMI)
C
                        KOFF1 = IT2SP(ISYMBIN,ISYME)
     *                     + NCKI(ISYMBIN)*(E-1)
     *                     + ICKI(ISYMBI,ISYMN)
     *                     + NT1AM(ISYMBI)*(N-1)
     *                     + IT1AM(ISYMIB,ISYMI)
     *                     + NVIR(ISYMIB)*(I-1)
     *                     + IB
C
                        KOFF2 = KSCR2 - 1
     *                     + ICKI(ISYMEN,ISYMI)
     *                     + NT1AM(ISYMEN)*(I-1)
     *                     + IT1AM(ISYME,ISYMN)
     *                     + NVIR(ISYME)*(N-1)
     *                     + E
C
                        WORK(KOFF2) = XOOVV(KOFF1)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C                tD(mle) * xB(eni) = tdxB(mlni)
C
         DO ISYMI = 1,NSYM
            ISYMEN = MULD2H(ISYMENI,ISYMI)
            DO ISYMN = 1,NSYM 
               ISYME = MULD2H(ISYMEN,ISYMN)
               ISYMML = MULD2H(ISYMMLE,ISYME)
               ISYMMLN = MULD2H(ISYMML,ISYMN)
               DO I = 1,NRHF(ISYMI)
 
                  KOFF1 = KSCR1 
     *               + IMAIJA(ISYMML,ISYME)
                  KOFF2 = KSCR2 
     *               + ICKI(ISYMEN,ISYMI) 
     *               + NT1AM(ISYMEN)*(I-1)
     *               + IT1AM(ISYME,ISYMN)
                  KOFF3 = KSCR3
     *               + I3ORHF(ISYMMLN,ISYMI)
     *               + NMAIJK(ISYMMLN) * (I-1) 
     *               + IMAIJK(ISYMML,ISYMN)
C
                  NTOTML = MAX(1,NMATIJ(ISYMML))
                  NTOTE  = MAX(1,NVIR(ISYME))
C
                  CALL DGEMM('N','N',NMATIJ(ISYMML),NRHF(ISYMN),
     *                     NVIR(ISYME),-ONE,WORK(KOFF1),
     *                     NTOTML,WORK(KOFF2),NTOTE,
     *                     ONE,WORK(KOFF3),NTOTML)
C
               END DO
            END DO
         END DO
C
C omega1(ai) = omega1(ai) +  TMAT(amln) * tdxB(mlni)
C
         DO ISYMI = 1,NSYM
            ISYMMLN = MULD2H(ISYMLNI,ISYMI)
            ISYMA   = MULD2H(ISYRES,ISYMI)
C
            KOFF1 = I3VOOO(ISYMA,ISYMMLN) + 1
            KOFF2 = KSCR3  + I3ORHF(ISYMMLN,ISYMI)
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
            NTOTMLN = MAX(1,NMAIJK(ISYMMLN))
            NTOTA   = MAX(1,NVIR(ISYMA))
C
C 4.3 omega contribution addomega
C
            CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *               NMAIJK(ISYMMLN),-ONE,TMAT(KOFF1),
     *               NTOTA,WORK(KOFF2),NTOTMLN,
     *               ONE,OMEGA1(KOFF3),NTOTA)
C
         END DO
C
      END IF
C
C================================================
C 5.  Calculate contributions from g_{ovvo}
C     - L^{def}_{lin} t^{de}_{lm} g_{mafn}
C -( Wfe(dlin) + Wfd(eiln) + Wde(fnil) )* t(dl,em) * g(mafn)
C================================================
C
      ISYDLM = MULD2H(ISYMT2,ISYMID)
      ISYTMP = MULD2H(ISYDLM,ISCKIJ)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NMAIJK(ISYTMP)
      KEND1 = KSCR2 + NCKI(ISYDLM)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_W3_OMEGA1 (OVVO-4)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NMAIJK(ISYTMP))
C
C----------------------------------------------
C -( Wfe(dlin) + Wfd(eiln) )* t(dl,em) * g(mafn)
C----------------------------------------------
C    
C - Wfe(dlin) * t(dl,em) * g(mafn)
C
C  Tt2(inm) = TMAT(dlin) * t2tp(dlmD)
C
C - Wfd(eiln) * t(dl,em) * g(mafn)
C
C Tt2(inm) = Tt2(inm) + TMAT(elin) * t2tp(emlD)
C                                    tD(elm)
C g(maBn) = xovvo(Bnma) = IB(mna)
C 
C Tt2(inm) = T(mni)
C
C omega(ai) = omega(ai) + IB(mna) * T(mni) 
C
C-----------------------------------------------
C     First contribution to intermediate
C-----------------------------------------------
C
C - Wfe(dlin) * t(dl,em) * g(mafn)
      DO I = 1, LENGTH
         TMAT(I) = WMAT(I)
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-6)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
C  Tt2(inm) = TMAT(dlin) * t2tp(dlmD)
C
      DO ISYMDL = 1, NSYM
         ISYMM  = MULD2H(ISYDLM,ISYMDL)
         ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
         KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
         KOFF2 = IT2SP(ISYDLM,ISYMID)
     *         + NCKI(ISYDLM)*(ID-1)
     *         + ICKI(ISYMDL,ISYMM)
     *         + 1
         KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
         NTOTDL = MAX(1,NT1AM(ISYMDL))
         NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,T2TP(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
      ENDDO
c
*     write(lupri,*)'kscr1(on,m) in cc3_w3_lhtr ib,id ',ib,id
*     call output(work(kscr1),1,nmatij(1),1,nrhf(1),
*    *            nmatij(1),nrhf(1),
*    *            1,lupri)

C
C-----------------------------------------------
C     Second contribution to intermediate
C-----------------------------------------------
C
C Tt2(inm) = Tt2(inm) + TMAT(elin) * t2tp(emlD)
C                                    tD(elm)
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C           
         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,1))
         ENDDO
C
         IF (NSYM .GT. 1) THEN
            IF (LWRK1 .LT. LENGTH) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-6)')
            ENDIF
            CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
         ENDIF
C
C t2tp(emlD) = tD(elm)
C
         DO ISYMM = 1, NSYM
            ISYMDL = MULD2H(ISYDLM,ISYMM)
            DO ISYMD = 1, NSYM
               ISYML  = MULD2H(ISYMDL,ISYMD)
               ISYMDM = MULD2H(ISYMD,ISYMM)
               DO M = 1, NRHF(ISYMM)
                  DO L = 1, NRHF(ISYML)
C
                     KOFF1 = IT2SP(ISYDLM,ISYMID)
     *                  + NCKI(ISYDLM)*(ID-1)
     *                  + ICKI(ISYMDL,ISYMM)
     *                  + NT1AM(ISYMDL)*(M-1)
     *                  + IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L-1)
     *                  + 1
C
                     KOFF2 = KSCR2
     *                  + ICKI(ISYMDM,ISYML)
     *                  + NT1AM(ISYMDM)*(L-1)
     *                  + IT1AM(ISYMD,ISYMM)
     *                  + NVIR(ISYMD)*(M-1)
C
                     CALL DCOPY(NVIR(ISYMD),T2TP(KOFF1),1,
     *                       WORK(KOFF2),1)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C Tt2(inm) = Tt2(inm) + TMAT(elin) * t2tp(emlD)
C                                    tD(elm)
         DO ISYMDL = 1, NSYM
            ISYMM  = MULD2H(ISYDLM,ISYMDL)
            ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
            KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
            KOFF2 = KSCR2
     *         + ICKI(ISYMDL,ISYMM)
            KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
            NTOTDL = MAX(1,NT1AM(ISYMDL))
            NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
            CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,WORK(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
         ENDDO
c
*     write(lupri,*)'kscr1(on,m) in cc3_w3_lhtr ib,id ',ib,id
*     call output(work(kscr1),1,nmatij(1),1,nrhf(1),
*    *            nmatij(1),nrhf(1),
*    *            1,lupri)

C
      END IF
C
C-------------------------------------------------------
C     Sort intermediate and integrals and contract
C-------------------------------------------------------
C
      KEND1 = KSCR2 + NMAIJK(ISYTMP)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_W3_OMEGA1 (OVVO-5)')
      ENDIF
C
C Tt2(inm) = T(mni)
C
      DO ISYMM = 1, NSYM
         ISYMIN = MULD2H(ISYTMP,ISYMM)
         DO ISYMN = 1, NSYM
            ISYMI  = MULD2H(ISYMIN,ISYMN)
            ISYMMN = MULD2H(ISYMM,ISYMN)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = KSCR1
     *                  + IMAIJK(ISYMIN,ISYMM)
     *                  + NMATIJ(ISYMIN)*(M-1)
     *                  + IMATIJ(ISYMI,ISYMN)
     *                  + NRHF(ISYMI)*(N-1)
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJK(ISYMMN,ISYMI)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M
C
                  CALL DCOPY(NRHF(ISYMI),WORK(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL DCOPY(NMAIJK(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
      ISYAMN = MULD2H(ISYINT,ISYMIB)
C
C g(maBn) = xovvo(Bnma) = IB(mna)
C
      DO ISYMA = 1, NSYM
         ISYMMN = MULD2H(ISYMA,ISYAMN)
         ISYBMN = MULD2H(ISYINT,ISYMA)
         DO ISYMM = 1, NSYM
            ISYMN  = MULD2H(ISYMMN,ISYMM)
            ISYMBN = MULD2H(ISYBMN,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = IT2SP(ISYBMN,ISYMA)
     *                  + ICKI(ISYMBN,ISYMM)
     *                  + NT1AM(ISYMBN)*(M-1)
     *                  + IT1AM(ISYMIB,ISYMN)
     *                  + NVIR(ISYMIB)*(N-1)
     *                  + IB
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJA(ISYMMN,ISYMA)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M
C
                  CALL DCOPY(NVIR(ISYMA),XOVVO(KOFF1),NCKI(ISYBMN),
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C omega(ai) = omega(ai) + IB(mna) * T(mni)
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMMN = MULD2H(ISYAMN,ISYMA)
C
         KOFF1 = KSCR2
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KSCR1
     *         + IMAIJK(ISYMMN,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NTOTMN = MAX(1,NMATIJ(ISYMMN))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
C 5.1 + 5.2 omega contribution addomega
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              -ONE,WORK(KOFF1),NTOTMN,WORK(KOFF2),NTOTMN,
     *              ONE,OMEGA1(KOFF3),NTOTA)
      ENDDO
C
C-------------------------------------
C     Third contribution
C
C - Wde(fnil) * t(dl,em) * g(mafn)
C-------------------------------------
C
C   TMAT(fnil) *  t2tp(BlmD       xovvo(fnma)
C
C   Tt2(fnim)  * xovvo(fnma)
C   tt2(fnmi)
C
C   omega(ai) = omega(ai) +  xoovv(fnma) * tt2(fnmi)
C
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C           

         ISYMBLM = MULD2H(ISYMT2,ISYMID)
         ISYMLM  = MULD2H(ISYMBLM,ISYMIB)
         ISYFNIM = MULD2H(ISYMLM,ISCKIJ)
C
         KSCR1 = 1
         KSCR2 = KSCR1 + NMATIJ(ISYMLM)
         KEND1 = KSCR2 + NCKIJ(ISYFNIM)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory CC3_W3_OMEG 6.)')
         ENDIF
C
C
         CALL DZERO(WORK(KSCR2),NCKIJ(ISYFNIM))
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(I)
         ENDDO

C
C -------------------------------
C     Sort T^{BD}_{mi} as T_{mi}
C -------------------------------
C
         ISYMBD = MULD2H(ISYMIB,ISYMID)
         ISYMMI = MULD2H(ISYMBD,ISYMT2)
         ISYBMI = MULD2H(ISYMMI,ISYMIB)
C
         DO ISYMI = 1,NSYM
            ISYMM = MULD2H(ISYMMI,ISYMI)
            ISYMBM = MULD2H(ISYBMI,ISYMI)
            DO I = 1,NRHF(ISYMI)
               DO M = 1,NRHF(ISYMM)
                  KOFF1 = IT2SP(ISYBMI,ISYMID)
     *                  + NCKI(ISYBMI)*(ID-1)
     *                  + ICKI(ISYMBM,ISYMI)
     *                  + NT1AM(ISYMBM)*(I-1)
     *                  + IT1AM(ISYMIB,ISYMM)
     *                  + NVIR(ISYMIB)*(M-1)
     *                  + IB
                  KOFF2 = IMATIJ(ISYMM,ISYMI)
     *                  + NRHF(ISYMM)*(I-1)
     *                  + M
                  WORK(KSCR1-1+KOFF2) = T2TP(KOFF1)
               ENDDO
            ENDDO
C
            KOFF3 = KSCR1 + IMATIJ(ISYMI,ISYMM)
C
         ENDDO
C
C - TMAT(fnil) *  t2tp(BlmD       xovvo(fnma)
C
         DO ISYML = 1,NSYM
            ISYMFNI = MULD2H(ISCKIJ,ISYML)
            ISYMM = MULD2H(ISYMLM,ISYML)
C
            KOFF1 = ISAIKJ(ISYMFNI,ISYML) + 1
            KOFF2 = KSCR1 + IMATIJ(ISYML,ISYMM)
            KOFF3 = KSCR2 + ISAIKJ(ISYMFNI,ISYMM)
C
            NTOTL   = MAX(1,NRHF(ISYML))
            NTOTFNI = MAX(1,NCKI(ISYMFNI))
C
            CALL DGEMM('N','N',NCKI(ISYMFNI),NRHF(ISYMM),NRHF(ISYML),
     *             -ONE,TMAT(KOFF1),NTOTFNI,WORK(KOFF2),NTOTL,
     *              ONE,WORK(KOFF3),NTOTFNI)
C
         END DO
C
C sort indices of result vedtor FNIM as FNMI
C
         IOPT = 2
         CALL CC3_LSORT2(WORK(KSCR2),ISYFNIM,WORK(KEND1),LWRK1,IOPT)
C
C   omega(ai) = omega(ai) +  xoovv(fnma) * tt2(fnmi)
C

         DO ISYMI = 1,NSYM
            ISYMFNM = MULD2H(ISYFNIM,ISYMI)
            ISYMA   = MULD2H(ISYRES,ISYMI)
C
            KOFF1 = IT2SP(ISYMFNM,ISYMA) + 1
            KOFF2 = KSCR2 + ISAIKJ(ISYMFNM,ISYMI)
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
            NTOTFNM = MAX(1,NCKI(ISYMFNM))
            NTOTA   = MAX(1,NVIR(ISYMA))
C
C 5.3 omega contribution addomega
C
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISYMFNM),
     *              -ONE,XOVVO(KOFF1),NTOTFNM,WORK(KOFF2),NTOTFNM,
     *              ONE,OMEGA1(KOFF3),NTOTA)
C
         END DO
C
      END IF
C
C================================================
C 6.  Calculate contributions from g_{oovv}
C     - L^{def}_{lni} t^{de}_{lm} g_{mnfa}
C -( Wfe(dlni) + Wfd(enli) + Wde(finl) )* t(dl,em) * g(mnfa)
C================================================
C
      ISYDLM = MULD2H(ISYMT2,ISYMID)
      ISYTMP = MULD2H(ISYDLM,ISCKIJ)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NMAIJK(ISYTMP)
      KEND1 = KSCR2 + NCKI(ISYDLM)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_W3_OMEGA1 (OVVO-4)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NMAIJK(ISYTMP))
C
C-----------------------------------------------
C     Calculate first and second contribution
C -( Wfe(dlni) + Wfd(enli) )* t(dl,em) * g(mnfa)
C-----------------------------------------------
C
C - Wfe(dlni) * t(dl,em) * g(mnfa)
C
C  Tt2(inm) = TMAT(dlin) * t2tp(dlmD)
C
C - Wfd(enli) * t(dl,em) * g(mnfa)
C
C  Tt2(inm) = Tt2(inm) + TMAT(elin) * t2tp(emlD)
C                                    tD(eml)
C  Tt2(inm) = T(mni)
C
C  g(mnBa) = xoovv(Bmna) = IB(mna)
C
C omega(ai) = omega(ai) + IB(mna) * T(mni)
C
C-----------------------------------------------
C     First contribution to intermediate
C-----------------------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(INDSQ(I,3))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-6)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
C  Tt2(inm) = TMAT(dlin) * t2tp(dlmD)
C
      DO ISYMDL = 1, NSYM
         ISYMM  = MULD2H(ISYDLM,ISYMDL)
         ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
         KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
         KOFF2 = IT2SP(ISYDLM,ISYMID)
     *         + NCKI(ISYDLM)*(ID-1)
     *         + ICKI(ISYMDL,ISYMM)
     *         + 1
         KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
         NTOTDL = MAX(1,NT1AM(ISYMDL))
         NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,T2TP(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
      ENDDO
C
C-----------------------------------------------
C     Second contribution to intermediate
C-----------------------------------------------
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C           

         DO I = 1, LENGTH
            TMAT(I) = WMAT(INDSQ(I,4))
         ENDDO
C
C  Tt2(inm) = Tt2(inm) + TMAT(elin) * t2tp(emlD)
C                                    tD(eml)
         IF (NSYM .GT. 1) THEN
            IF (LWRK1 .LT. LENGTH) THEN
               CALL QUIT('Out of memory in CC3_W3_OMEGA1 (CC_GATHER-6)')
            ENDIF
            CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
         ENDIF
C
         DO ISYMM = 1, NSYM
            ISYMDL = MULD2H(ISYDLM,ISYMM)
            DO ISYMD = 1, NSYM
               ISYML  = MULD2H(ISYMDL,ISYMD)
               ISYMDM = MULD2H(ISYMD,ISYMM)
               DO M = 1, NRHF(ISYMM)
                  DO L = 1, NRHF(ISYML)
C
                     KOFF1 = IT2SP(ISYDLM,ISYMID)
     *                  + NCKI(ISYDLM)*(ID-1)
     *                  + ICKI(ISYMDL,ISYMM)
     *                  + NT1AM(ISYMDL)*(M-1)
     *                  + IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L-1)
     *                  + 1
C
                     KOFF2 = KSCR2
     *                  + ICKI(ISYMDM,ISYML)
     *                  + NT1AM(ISYMDM)*(L-1)
     *                  + IT1AM(ISYMD,ISYMM)
     *                  + NVIR(ISYMD)*(M-1)
C
                     CALL DCOPY(NVIR(ISYMD),T2TP(KOFF1),1,
     *                       WORK(KOFF2),1)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
         DO ISYMDL = 1, NSYM
            ISYMM  = MULD2H(ISYDLM,ISYMDL)
            ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
            KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
            KOFF2 = KSCR2
     *         + ICKI(ISYMDL,ISYMM)
            KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
            NTOTDL = MAX(1,NT1AM(ISYMDL))
            NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
            CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,WORK(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
         ENDDO
C
      END IF
C
C-------------------------------------------------------
C     Sort intermediate and integrals and contract
C-------------------------------------------------------
C
      KEND1 = KSCR2 + NMAIJK(ISYTMP)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_W3_OMEGA1 (OVVO-5)')
      ENDIF
C
C Tt2(inm) = T(mni)
C
      DO ISYMM = 1, NSYM
         ISYMIN = MULD2H(ISYTMP,ISYMM)
         DO ISYMN = 1, NSYM
            ISYMI  = MULD2H(ISYMIN,ISYMN)
            ISYMMN = MULD2H(ISYMM,ISYMN)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = KSCR1
     *                  + IMAIJK(ISYMIN,ISYMM)
     *                  + NMATIJ(ISYMIN)*(M-1)
     *                  + IMATIJ(ISYMI,ISYMN)
     *                  + NRHF(ISYMI)*(N-1)
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJK(ISYMMN,ISYMI)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M
C
                  CALL DCOPY(NRHF(ISYMI),WORK(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL DCOPY(NMAIJK(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
      ISYAMN = MULD2H(ISYINT,ISYMIB)
C
C  g(mnBa) = xoovv(Bmna) = IB(mna)
C
      DO ISYMA = 1, NSYM
         ISYMMN = MULD2H(ISYMA,ISYAMN)
         ISYBMN = MULD2H(ISYINT,ISYMA)
         DO ISYMM = 1, NSYM
            ISYMN  = MULD2H(ISYMMN,ISYMM)
            ISYMBN = MULD2H(ISYBMN,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = IT2SP(ISYBMN,ISYMA)
     *                  + ICKI(ISYMBN,ISYMM)
     *                  + NT1AM(ISYMBN)*(M-1)
     *                  + IT1AM(ISYMIB,ISYMN)
     *                  + NVIR(ISYMIB)*(N-1)
     *                  + IB
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJA(ISYMMN,ISYMA)
     *                  + IMATIJ(ISYMN,ISYMM)
     *                  + NRHF(ISYMN)*(M-1)
     *                  + N
C
                  CALL DCOPY(NVIR(ISYMA),XOOVV(KOFF1),NCKI(ISYBMN),
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C omega(ai) = omega(ai) + IB(mna) * T(mni)
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMMN = MULD2H(ISYAMN,ISYMA)
C
         KOFF1 = KSCR2
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KSCR1
     *         + IMAIJK(ISYMMN,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NTOTMN = MAX(1,NMATIJ(ISYMMN))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
C 6.1 + 6.2 omega contribution addomega
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              -ONE,WORK(KOFF1),NTOTMN,WORK(KOFF2),NTOTMN,
     *              ONE,OMEGA1(KOFF3),NTOTA)
      ENDDO
C-----------------------------------------------
C     Calculate third contribution
C - Wde(finl) * t(dl,em) * g(mnfa)
C----------------------------------------------
C
C - TMAT(finl) * t2tp(BlmD)    xoovv(fmna)
C
C   Tt2(finm)  *  xoovv(fmna)
C   tt2(fmni)
C
C   omega(ai) = omega(ai) +  xoovv(fmna) * tt2(fmni)
C
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C           

         ISYMBLM = MULD2H(ISYMT2,ISYMID)
         ISYMLM  = MULD2H(ISYMBLM,ISYMIB)
         ISYFINM = MULD2H(ISYMLM,ISCKIJ)
C
         KSCR1 = 1
         KSCR2 = KSCR1 + NMATIJ(ISYMLM)
         KEND1 = KSCR2 + NCKIJ(ISYFINM)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory CC3_W3_OMEG 6.)')
         ENDIF
C
C
         CALL DZERO(WORK(KSCR2),NCKIJ(ISYFINM))
C
         DO I = 1, LENGTH
            TMAT(I) = WMAT(I)
         ENDDO
C
C -------------------------------
C     Sort T^{BD}_{mi} as T_{mi}
C -------------------------------
C
         ISYMMI = MULD2H(ISYMBD,ISYMT2)
         ISYBMI = MULD2H(ISYMMI,ISYMIB)
C
         DO ISYMI = 1,NSYM
            ISYMM = MULD2H(ISYMMI,ISYMI)
            ISYMBM = MULD2H(ISYBMI,ISYMI)
            DO I = 1,NRHF(ISYMI)
               DO M = 1,NRHF(ISYMM)
                  KOFF1 = IT2SP(ISYBMI,ISYMID)
     *                  + NCKI(ISYBMI)*(ID-1)
     *                  + ICKI(ISYMBM,ISYMI)
     *                  + NT1AM(ISYMBM)*(I-1)
     *                  + IT1AM(ISYMIB,ISYMM)
     *                  + NVIR(ISYMIB)*(M-1)
     *                  + IB
                  KOFF2 = IMATIJ(ISYMM,ISYMI)
     *                  + NRHF(ISYMM)*(I-1)
     *                  + M
                  WORK(KSCR1-1+KOFF2) = T2TP(KOFF1)
               ENDDO
            ENDDO
         ENDDO
C
C - TMAT(finl) * t2tp(BlmD)    xoovv(fmna)
C
         ISYMLM = MULD2H(ISYMBD,ISYMT2)
         DO ISYML = 1,NSYM
            ISYMFIN = MULD2H(ISCKIJ,ISYML)
            ISYMM = MULD2H(ISYMLM,ISYML)
C
            KOFF1 = ISAIKJ(ISYMFIN,ISYML) + 1
            KOFF2 = KSCR1 + IMATIJ(ISYML,ISYMM)
            KOFF3 = KSCR2 + ISAIKJ(ISYMFIN,ISYMM)
C
            NTOTL   = MAX(1,NRHF(ISYML))
            NTOTFIN = MAX(1,NCKI(ISYMFIN))
C
            CALL DGEMM('N','N',NCKI(ISYMFIN),NRHF(ISYMM),NRHF(ISYML),
     *             -ONE,TMAT(KOFF1),NTOTFIN,WORK(KOFF2),NTOTL,
     *              ONE,WORK(KOFF3),NTOTFIN)
         END DO
C
C sort indices of result vedtor FINM as FMNI
C
         IOPT = 4 
         CALL CC3_LSORT2(WORK(KSCR2),ISYFINM,WORK(KEND1),LWRK1,IOPT)
C
C   omega(ai) = omega(ai) +  xoovv(fmna) * tt2(fmni)
C
         DO ISYMI = 1,NSYM
            ISYMFMN = MULD2H(ISYFINM,ISYMI)
            ISYMA   = MULD2H(ISYRES,ISYMI)
C
            KOFF1 = IT2SP(ISYMFMN,ISYMA) + 1   
            KOFF2 = KSCR2 + ISAIKJ(ISYMFMN,ISYMI)
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
            NTOTFMN = MAX(1,NCKI(ISYMFMN))
            NTOTA   = MAX(1,NVIR(ISYMA))
C
C 6.3 omega contribution addomega
C
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISYMFMN),
     *              -ONE,XOOVV(KOFF1),NTOTFMN,WORK(KOFF2),NTOTFMN,
     *              ONE,OMEGA1(KOFF3),NTOTA)
         END DO
C
      END IF
C
C-------------------
C     End.
C-------------------
C
      CALL QEXIT('CC3_W3_OMEGA1')
C
      RETURN
      END
C  /* Deck cc3_srtvooo */
      SUBROUTINE CC3_SRTVOOO(GVPOOO,GVOOO,ISCKIJ)
C
C   Sort matrix  of symmetry ISCKIJ stored in GVOOO(V,O,O,O) 
C   as GVPOOO(V,OOO) 
C
C
      IMPLICIT NONE
C

      INTEGER ISCKIJ, ISYMJ, ISYCKI, ISYMI, ISYMCK, ISYMIJ, ISYMC, ISYMK
      INTEGER ISYMKI, ISYKIJ, NKIJ, KVPOOO, KVOOO 
C
#if defined (SYS_CRAY)
      REAL GVPOOO(*),GVOOO(*)
#else
      DOUBLE PRECISION GVPOOO(*),GVOOO(*)
#endif
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

C
      CALL QENTER('CC3_SRTVOOO')
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYCKI = MULD2H(ISYMJ,ISCKIJ)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMCK = MULD2H(ISYMI,ISYCKI)
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
            DO 120 ISYMC = 1,NSYM
C
               ISYMK  = MULD2H(ISYMCK,ISYMC) 
               ISYMKI = MULD2H(ISYMK,ISYMI)
               ISYKIJ = MULD2H(ISYMIJ,ISYMK)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                   DO 140 I = 1,NRHF(ISYMI)
C
                      DO 150 K = 1,NRHF(ISYMK)
C
                         NKIJ = IMAIJK(ISYMKI,ISYMJ)
     *                        + NMATIJ(ISYMKI)*(J-1)
     *                        + IMATIJ(ISYMK,ISYMI)
     *                        + NRHF(ISYMK)*(I - 1) + K
C
                         DO 160 C = 1,NVIR(ISYMC)
C
                            KVOOO  = ISAIKJ(ISYCKI,ISYMJ)
     *                             + NCKI(ISYCKI)*(J-1)
     *                             + ISAIK(ISYMCK,ISYMI)
     *                             + NT1AM(ISYMCK)*(I - 1)
     *                             + IT1AM(ISYMC,ISYMK)
     *                             + NVIR(ISYMC)*(K - 1)
     *                             + C
C
                            KVPOOO = I3VOOO(ISYMC,ISYKIJ)
     *                             + NVIR(ISYMC)*(NKIJ - 1)
     *                             + C
C
                            GVPOOO(KVPOOO) =  GVOOO(KVOOO)
C
  160                    CONTINUE
  150                 CONTINUE
  140              CONTINUE
  130           CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_SRTVOOO')
C
      RETURN
      END

C  /* Deck cc3_w3_cy2v */
      SUBROUTINE CC3_W3_CY2V(XI2,ISYRES,RBJIA,WMAT,ISWMAT,TMAT,TRVIR,
     *                      TRVIR1,ISYINT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMIB,IB,ISYMID,ID,W3X)
C
C    
C Modified based on CC3_CY2V routine where TWO*WMAT is replaced by WMAT
C and the other terms in TMAT are removed.
C
C General symmetry: ISWMAT is symmetry of WMAT and TMAT intermediates.
C                   (exclusive isymib,isymid)
C                   ISYINT is symmetry of integrals in TRVIR and TROVIR1.
C                   ISYRES = ISWMAT*ISYINT*ISYMIB*ISYMID
C                   RBJIA intermediate result vector
C
C
C If W3X = .TRUE., then WMAT contains wbar_X and all terms 
C 1-3 are included
C else
C we assume that we get tbar_0 in as WMAT and calculate only
C term 3. 
C
C Note Wbar_X and Tbar_0 is stored in the same way for fixed IB and ID
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{cdk} (2W^{bcd}_{jik}-W^{bcd}_{kij}-W^{bcd}_{jki}) * (ac|kd)
C (1):
C     = sum_{ck} (2W^{DB}_{cijk}-W^{DB}_{cikj}-W^{DB}_{ckji}) * (ac|kD)
C (2):
C     + sum_{dk} (2W^{BC}_{dkij}-W^{BC}_{djik}-W^{BC}_{dikj}) * (aC|kd) 
C (3):
C    + sum_{k} (2W^{CD}_{bjki}-W^{CD}_{bkji}-W^{CD}_{bjik}) * (aC|kD)
C
      IMPLICIT NONE
C
      LOGICAL W3X
C
      INTEGER ISYRES, ISWMAT, ISYINT , LENSQ, LWORK
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYMIB, IB, ISYMID, ID
      INTEGER INDEX ,LENGTH, ISCKIJ, ISYMB, ISYAIJ, KRMAT, KEND, LEND
      INTEGER ISYMJ ,ISYBJ, ISYAI, ISYCKI, ISYCK, ISYMA, NTOTCK, NVIRA
      INTEGER KOFF1 , KOFF2, KOFF3, ISYMD, ISDKIJ, ISYDKI
      INTEGER ISYDK, NTOTDK, ISBJIK, ISYCKA
      INTEGER ISYKA, KINT, ISYBJI, KOFFT, KOFFM, KOFFR
      INTEGER NTOK, NTOBJI, ISYMK, ISYMC, ISYMI
#if defined (SYS_CRAY)
      REAL XI2(*), RBJIA(*), WMAT(*), TMAT(*), TRVIR(*), TRVIR1(*) 
      REAL WORK(*)
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION XI2(*), RBJIA(*), WMAT(*), TMAT(*) 
      DOUBLE PRECISION TRVIR(*), TRVIR1(*), WORK(*)
      DOUBLE PRECISION ZERO, ONE, TWO
#endif

C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_W3_CY2V')
      LENGTH = NCKIJ(ISWMAT)
C
C------------------------
C (1):
C     = sum_{ck} (2W^{DB}_{cijk}-W^{DB}_{cikj}-W^{DB}_{ckji}) * (ac|kD)
C
C------------------------
C
C   TMAT(ckij) = 2W^{DB}_{cijk}-W^{DB}_{cikj}-W^{DB}_{ckji}
C
C   Use here the symmetry between BJ, DK !!
C
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C           

         ISCKIJ = ISWMAT 
C
         ISYMB = ISYMIB
         ISYMD = ISYMID
         B = IB
         D = ID
C
         ISYAIJ = MULD2H(ISYMB,ISYRES)
         KRMAT  = 1
         KEND   = KRMAT + NCKI(ISYAIJ)
         LEND   = LWORK - KEND
         IF (LEND .LT. 0)THEN
            CALL QUIT('1. Insufficient core in CC3_W3_CY2V')
         ENDIF
         CALL DZERO(WORK(KRMAT), NCKI(ISYAIJ))
         DO I = 1,LENGTH
C
            TMAT(I) =  WMAT(INDSQ(I,1))
C
         ENDDO
C
C---------------------------
C     (ac|kD) = TRVIR^D_{ck,a}
C
C     RMAT^B_(aij) =  TRVIR^D_{ck,a} )^T *  TMAT(ckij)
C---------------------------
C
         DO 200 ISYMJ = 1,NSYM
C
            ISYBJ = MULD2H(ISYMB,ISYMJ)
            ISYAI = MULD2H(ISYBJ,ISYRES)
            ISYCKI = MULD2H(ISCKIJ,ISYMJ)
C
            DO 210 J = 1,NRHF(ISYMJ)
C
               DO 220 ISYMI = 1,NSYM
C
                  ISYCK = MULD2H(ISYCKI,ISYMI)
                  ISYMA  = MULD2H(ISYAI,ISYMI)
C 
                  NTOTCK = MAX(NT1AM(ISYCK),1)
                  NVIRA  = MAX(NVIR(ISYMA),1)
C
                  KOFF1 = ICKATR(ISYCK,ISYMA) + 1
                  KOFF2 = ISAIKJ(ISYCKI,ISYMJ)
     *               + NCKI(ISYCKI)*(J - 1)
     *               + ISAIK(ISYCK,ISYMI)  + 1
                  KOFF3 = KRMAT
     *               + ISAIK(ISYAI,ISYMJ)
     *               + NT1AM(ISYAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI)
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                    NT1AM(ISYCK),
     *                   -ONE,TRVIR(KOFF1),NTOTCK,TMAT(KOFF2),NTOTCK,
     *                    ONE,WORK(KOFF3),NVIRA) 
C
  220          CONTINUE
C
  210       CONTINUE
C
  200    CONTINUE
C
         ISYAIJ = MULD2H(ISYMB,ISYRES)
         IF ( NCKI(ISYAIJ).GT.0 ) THEN
           CALL CC3_RACC(XI2,WORK(KRMAT),ISYMB,B,ISYRES)
         END IF
C
C
C (2):
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{dk} (2W^{BC}_{dkij}-W^{BC}_{djik}-W^{BC}_{dikj}) * (aC|kd)
C
C  TMAT(dkij) = 2W^{BC}_{dkij}-W^{BC}_{djik}-W^{BC}_{dikj}
C
         ISYMC = ISYMID
         ISYMB = ISYMIB
         B  =  IB
         C  =  ID 
C
         ISDKIJ = ISWMAT
         ISYAIJ = MULD2H(ISYMB,ISYRES)
         KRMAT  = 1
         KEND   = KRMAT + NCKI(ISYAIJ)
         LEND   = LWORK - KEND
         IF (LEND .LT. 0)THEN
            CALL QUIT('2. Insufficient core in CC3_W3_CY2V')
         ENDIF
         CALL DZERO(WORK(KRMAT), NCKI(ISYAIJ))
C 
         DO I = 1,LENGTH
C
            TMAT(I) =  WMAT(I) 
C
         ENDDO
C
C     (ad|kC) = TRVIR1^C(dka)
C
C  RMAT^B_(aij) = ( TRVIR1^C(dka) )^T * TMAT(dkij)
C
         DO 300 ISYMJ = 1,NSYM
C
            ISYBJ = MULD2H(ISYMB,ISYMJ)
            ISYAI = MULD2H(ISYBJ,ISYRES)
            ISYDKI = MULD2H(ISDKIJ,ISYMJ)
C
            DO 310 J = 1,NRHF(ISYMJ)
C
               DO 320 ISYMI = 1,NSYM
C
                  ISYDK = MULD2H(ISYDKI,ISYMI)
                  ISYMA  = MULD2H(ISYAI,ISYMI)
C
                  NTOTDK = MAX(NT1AM(ISYDK),1)
                  NVIRA  = MAX(NVIR(ISYMA),1)
C
                  KOFF1 = ICKATR(ISYDK,ISYMA) + 1
                  KOFF2 = ISAIKJ(ISYDKI,ISYMJ)
     *               + NCKI(ISYDKI)*(J - 1)
     *               + ISAIK(ISYDK,ISYMI)  + 1
                  KOFF3 = KRMAT + ISAIK(ISYAI,ISYMJ)
     *               + NT1AM(ISYAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI) 
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                    NT1AM(ISYDK),
     *                   -ONE,TRVIR1(KOFF1),NTOTDK,TMAT(KOFF2),NTOTDK,
     *                    ONE,WORK(KOFF3),NVIRA)
C
  320          CONTINUE
  310       CONTINUE
  300    CONTINUE
C
         ISYAIJ = MULD2H(ISYMB,ISYRES)
         IF ( NCKI(ISYAIJ).GT.0 ) THEN
           CALL CC3_RACC(XI2,WORK(KRMAT),ISYMB,B,ISYRES)
         END IF
C
      END IF
C
C-----------------------------------
C (3):
C    + sum_{k} (2W^{CD}_{bjki}-W^{CD}_{bkji}-W^{CD}_{bjik}) * (ac|kD)
C
C-----------------------------------
C
      ISYMC = ISYMIB
      ISYMD = ISYMID
      C = IB
      D = ID
C
C     TMAT(bjik) = 2W^{CD}_{bjki}-W^{CD}_{bkji}-W^{CD}_{bjik}
C
      ISBJIK = ISWMAT
C
      DO I = 1,LENGTH
C
         TMAT(I) =  WMAT(INDSQ(I,3))
C
      ENDDO
C
C     (ac|kD) = TRVIR^D_{Cka} = M^DC_{ka}
C
c
      ISYCKA = MULD2H(ISYINT,ISYMD)
      ISYKA  = MULD2H(ISYCKA,ISYMC)
C
      DO  ISYMA = 1,NSYM
         ISYCK  = MULD2H(ISYCKA,ISYMA)
         ISYMK  = MULD2H(ISYMC,ISYCK)
         DO A = 1,NVIR(ISYMA)
            DO K = 1,NRHF(ISYMK) 
               KOFFM = IT1AMT(ISYMK,ISYMA)
     *               + NRHF(ISYMK)*(A-1) + K 
               KINT  = ICKATR(ISYCK,ISYMA)
     *               + NT1AM(ISYCK)*(A-1) 
     *               + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K-1) + C 
               WORK(KOFFM) = TRVIR(KINT)
            END DO
         END DO
      END DO
C
C  RBJIA(bjia) = RBJIA(bjia) + TMAT(bjik) *  M^DC_{ka}
C
      DO ISYMA = 1,NSYM
         ISYMK = MULD2H(ISYMA,ISYKA)
         ISYBJI = MULD2H(ISBJIK,ISYMK)
         KOFFT  = ISAIKJ(ISYBJI,ISYMK) + 1 
         KOFFM  = IT1AMT(ISYMK,ISYMA)  + 1
         KOFFR  = IT2SP(ISYBJI,ISYMA)  + 1
         NTOK = MAX(NRHF(ISYMK),1)
         NTOBJI = MAX(NCKI(ISYBJI),1)
C
         CALL DGEMM('N','N',NCKI(ISYBJI),NVIR(ISYMA),
     *               NRHF(ISYMK),-ONE,TMAT(KOFFT),NTOBJI, 
     *               WORK(KOFFM),NTOK,ONE,RBJIA(KOFFR),NTOBJI)
C
      END DO
C
      CALL QEXIT('CC3_W3_CY2V')
      RETURN
      END

C  /* Deck cc3_w3_cy2o */
      SUBROUTINE CC3_W3_CY2O(XI2,ISYRES,WMAT,ISWMAT,TMAT,TROCC,
     *                      TROCC1,ISYINT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMIB,IB,ISYMID,ID,W3X)
C
C
C    
C Modified based on CC3_CY2O routine where TWO*WMAT is replaced by WMAT
C and the other terms in TMAT are removed.
C
C     General symmetry: ISWMAT is symmetry of WMAT and TMAT intermediates.
C                       (exclusive isymib,isymid)
C                       ISYINT is symmetry of integrals in TROCC and TROCC1.
C                       ISYRES = ISWMAT*ISYINT*ISYMIB*ISYMID
C
C If W3X = .TRUE., then WMAT contains wbar_X and all terms
C 1-3 are included
C else
C we assume that we get tbar_0 in as WMAT and calculate only
C term 1.
C
C Note Wbar_X and Tbar_0 is stored in the same way for fixed IB and ID

C
C   XI2(aibj) =  XI2(aibj) 
C
C     + sum_{ckl} (2W^{abc}_{ikl}-W^{abc}_{lki}-W^{abc}_{ilk}) * (lc|kj)
C (1):
C     = sum_{kl} (2W^{BC}_{ailk}-W^{BC}_{alik}-W^{BC}_{aikl}) * (lc|kj) 
C (2):
C     + sum_{kl} (2W^{CA}_{bkil}-W^{CA}_{bkli}-W^{CA}_{blik}) * (lc|kj) 
C (3):
C    + sum_{ckl} (2W^{AB}_{clki}-W^{AB}_{cikl}-W^{AB}_{ckli}) * (lc|kj) 
C
      IMPLICIT NONE
C
      LOGICAL W3X
C
      INTEGER ISYRES, ISWMAT, ISYINT, LWORK, LENSQ 
      INTEGER ISYMIB, IB, ISYMID, ID 
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYMC, ISYMB, ISYBC, JSAIKL, LENGTH, ISYMJ, ISYBJ
      INTEGER ISYAI, ISYKL, NTOTAI 
      INTEGER NTOTKL, KOFF1, KOFF2, KOFF3 
      INTEGER ISYMA, ISYAB, JSCKLI, JSBIKL, ISYMI
      INTEGER ISYCKL, NTOCKL, NRHFI 
      INTEGER KRMAT1, KRMAT2, KEND, LEND, INDEX
      INTEGER ISYCK, ISYCJ, ISYAC, ISYBI, ISYKLJ, NTOTBI
      INTEGER ISWB, ISWD, NCKIMX
      INTEGER ISYAIJ, ISYBJI, ISYIJ 
      INTEGER KTMATX
C
#if defined (SYS_CRAY) 
      REAL XI2(*), WMAT(*), TMAT(*), TROCC(*), TROCC1(*)
      REAL WORK(*) 
      REAL TWO, ONE, ZERO, XTMAT, XRMAT, DDOT
#else
      DOUBLE PRECISION XI2(*), WMAT(*), TMAT(*)
      DOUBLE PRECISION TROCC(*), TROCC1(*) 
      DOUBLE PRECISION WORK(*) 
      DOUBLE PRECISION TWO, ONE, ZERO, XTMAT, XRMAT, DDOT
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_W3_CY2O')
C
      
      ISWB   = MULD2H(ISWMAT,ISYMIB)
      ISWD   = MULD2H(ISWMAT,ISYMID)
      ISYAIJ = MULD2H(ISYMIB,ISYRES)
      ISYBJI = MULD2H(ISYMID,ISYRES)
      NCKIMX = MAX(NCKI(ISWB),NCKI(ISWD),NCKI(ISYAIJ),NCKI(ISYBJI)) 
C
      KRMAT1 = 1
      KRMAT2 = KRMAT1 + NCKIMX
      KTMATX = KRMAT2 + NCKIMX
      KEND   = KTMATX + NCKIJ(ISWMAT)
      LEND   = LWORK  - KEND
      IF (LWORK .LT. KEND) THEN
         CALL QUIT('Insufficient space in CY2O')
      ENDIF
C
C-------------------------
C     First occupied term.
C-------------------------
C (1):
C             XI2(aibj) =  XI2(aibj) 
C
C     + sum_{kl} (2W^{BC}_{ailk}-W^{BC}_{alik}-W^{BC}_{aikl}) * (lc|kj)

C
      B = IB
      C = ID
C
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
      IF (NCKI(ISYAIJ).GT.0 ) THEN 
        CALL DZERO(WORK(KRMAT1),NCKI(ISYAIJ))
C
        ISYBC = MULD2H(ISYMB,ISYMC)
        JSAIKL = ISWMAT 
C
        LENGTH = NCKIJ(JSAIKL)
C
C---------------------------------------
C
C     TMAT(aikl) = (2W^{BC}_{ailk}-W^{BC}_{alik}-W^{BC}_{aikl}) 
C
C---------------------------------------
C
        DO I = 1,LENGTH
C
           TMAT(I) =   WMAT(INDSQ(I,3))
C
        ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
      
        IF (NSYM .GT. 1) THEN
           CALL CC_GATHER(LENGTH,WORK(KTMATX),TMAT,INDSQ(1,6))
           CALL DCOPY(LENGTH,WORK(KTMATX),1,TMAT,1)
        ENDIF
C
        IF (IPRINT .GT. 55) THEN
           XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
           WRITE(LUPRI,*) 'In CC3_W3_CY2O: 1. Norm of TMAT = ',XTMAT
        ENDIF
C
C-----------------------
C     First contraction.
C-----------------------
C
C         TROCC(kljC) = (lc|kj)
C
C   XI2(aibj) =  XI2(aibj) 
C             + sum(kl) TMAT(aikl)* TROCC(kljC) 
C
        DO 200 ISYMJ = 1,NSYM
C
           ISYBJ = MULD2H(ISYMB,ISYMJ)
           ISYAI = MULD2H(ISYBJ,ISYRES)
           ISYKL = MULD2H(JSAIKL,ISYAI)
           ISYKLJ = MULD2H(ISYKL,ISYMJ)
C
           NTOTAI = MAX(NT1AM(ISYAI),1)
           NTOTKL = MAX(NMATIJ(ISYKL),1)
C
           KOFF1  = ISAIKL(ISYAI,ISYKL) + 1
           KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *          + NMAJIK(ISYKLJ)*(C - 1)
     *          + ISJIK(ISYKL,ISYMJ) + 1
           KOFF3  = KRMAT1 + ISAIK(ISYAI,ISYMJ) 
C
           CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMJ),NMATIJ(ISYKL),
     *                ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *                ONE,WORK(KOFF3),NTOTAI) 
C
  200   CONTINUE
C
        IF (IPRINT .GT. 55) THEN
           XRMAT = DDOT(NCKI(ISYAIJ),WORK(KRMAT1),1,WORK(KRMAT1),1)
           WRITE(LUPRI,*) '1. In CC3_W3_CY2O: Norm of RMAT1 =  ',XRMAT
        ENDIF
C
        CALL CC3_RACC(XI2,WORK(KRMAT1),ISYMB,B,ISYRES)
C
      END IF
C
C--------------------------
C     Second occupied term.
C--------------------------
C
C (2):
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{kl} (2W^{CA}_{bkil}-W^{CA}_{bkli}-W^{CA}_{blik}) * (lc|kj)
C 
      ! Calculate only when contracting with first-order multipliers
      IF (W3X) THEN
C           

         A = ID
         C = IB
C
         ISYMA = ISYMID
         ISYMC = ISYMIB
C
         ISYBJI = MULD2H(ISYMA,ISYRES)
         IF (NCKI(ISYBJI).GT.0) THEN
           CALL DZERO(WORK(KRMAT1),NCKI(ISYBJI))
C
           ISYAC = MULD2H(ISYMA,ISYMC)
           JSBIKL = ISWMAT 
C
C---------------------------------------
C
C     TMAT(bikl) = (2W^{CA}_{bkil}-W^{CA}_{bkli}-W^{CA}_{blik}) 
C
C---------------------------------------
C
           LENGTH = NCKIJ(JSBIKL)
C
           DO I = 1,LENGTH
C
              TMAT(I) =   WMAT(INDSQ(I,1)) 
C
           ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
           IF (NSYM .GT. 1) THEN
              CALL CC_GATHER(LENGTH,WORK(KTMATX),TMAT,INDSQ(1,6))
              CALL DCOPY(LENGTH,WORK(KTMATX),1,TMAT,1)
           ENDIF
C
           IF (IPRINT .GT. 55) THEN
              XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
              WRITE(LUPRI,*) 'In CC3_W3_CY2O: 2. Norm of TMAT = ',XTMAT
           ENDIF
C
C------------------------
C     Second contraction.
C------------------------
C
C         TROCC(kljC) = (lc|kj)
C
C      M^A_(bij) =   sum(kl) TMAT(bikl)* TTROCC(kljC)
C

           DO 400 ISYMJ = 1,NSYM
C
              ISYCJ = MULD2H(ISYMC,ISYMJ)
              ISYKL = MULD2H(ISYINT,ISYCJ)
              ISYBI = MULD2H(ISYKL,ISWMAT)
              ISYKLJ = MULD2H(ISYKL,ISYMJ)
C
              NTOTBI = MAX(NT1AM(ISYBI),1)
              NTOTKL = MAX(NMATIJ(ISYKL),1)
C
              KOFF1  = ISAIKL(ISYBI,ISYKL) + 1
              KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *            + NMAJIK(ISYKLJ)*(C - 1)
     *            + ISJIK(ISYKL,ISYMJ) + 1
              KOFF3  = KRMAT1 + ISAIK(ISYBI,ISYMJ) 
C

              CALL DGEMM('N','N',NT1AM(ISYBI),NRHF(ISYMJ),NMATIJ(ISYKL),
     *                ONE,TMAT(KOFF1),NTOTBI,TROCC(KOFF2),NTOTKL,
     *                ONE,WORK(KOFF3),NTOTBI)
C
  400      CONTINUE
C

           IF (IPRINT .GT. 55) THEN
              XRMAT = DDOT(NCKI(ISYBJI),WORK(KRMAT1),1,WORK(KRMAT1),1)
              WRITE(LUPRI,*) '2.In CC3_W3_CY2O: Norm of RMAT1 =  ',XRMAT 
           ENDIF
C 
C reorder result vector M^A_(bij) as M^A_(bji)
C
           CALL CC3_MAJI(WORK(KRMAT1),WORK(KRMAT2),ISYMA,A,ISYRES) 
C
           IF (IPRINT .GT. 55) THEN
              XRMAT = DDOT(NCKI(ISYBJI),WORK(KRMAT2),1,WORK(KRMAT2),1)
              WRITE(LUPRI,*) 
     *         'In CC3_W3_CY2O: Reorder :Norm of RMAT2 =  ',XRMAT 
           ENDIF
C
C add vector to XI2
C
           CALL CC3_RACC(XI2,WORK(KRMAT2),ISYMA,A,ISYRES) 
C
         END IF
C
C-------------------------
C     Third occupied term.
C-------------------------
C (3):
C             XI2(aibj) =  XI2(aibj) 
C
C    + sum_{ckl} (2W^{AB}_{clki}-W^{AB}_{cikl}-W^{AB}_{ckli}) * (lc|kj) 
C
         B = ID
         A = IB
C
         ISYMB = ISYMID
         ISYMA = ISYMIB
C
         ISYAB = MULD2H(ISYMA,ISYMB)
         ISYIJ = MULD2H(ISYAB,ISYRES)
         IF (NMATIJ(ISYIJ).GT.0) THEN
           CALL DZERO(WORK,NMATIJ(ISYIJ))
C
           JSCKLI = ISWMAT 
C
           LENGTH = NCKIJ(JSCKLI)
C
C----------------------------------
C
C   TMAT(ckli) = 2W^{AB}_{clki}-W^{AB}_{cikl}-W^{AB}_{ckli}
C
C----------------------------------
C
           DO I = 1,LENGTH
C
              TMAT(I) =  WMAT(INDSQ(I,1))    
C
           ENDDO
C
C
           IF (IPRINT .GT. 55) THEN
              XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
              WRITE(LUPRI,*) 'In CC3_W3_CY2O: 3. Norm of TMAT = ',XTMAT
           ENDIF
C
C-----------------------
C     Third contraction.
C-----------------------
C
C
C         TROCC1(cklj) = (lc|kj)
C
C   XI2(aibj) =  XI2(aibj)
C             + sum(kl) TMAT(ckli)* TROCC1(cklj) 
C
C
           DO 600 ISYMJ = 1,NSYM
C
              ISYBJ = MULD2H(ISYMB,ISYMJ)
              ISYAI = MULD2H(ISYBJ,ISYRES)
              ISYMI  = MULD2H(ISYAI,ISYMA)
              ISYCKL = MULD2H(ISYMI,JSCKLI)
C
              IF (LWORK .LT. NRHF(ISYMI)*NRHF(ISYMJ)) THEN
                 CALL QUIT('Insufficient memory in CCSDT_CY2O')
              END IF
C
              NTOCKL = MAX(NCKI(ISYCKL),1)
              NRHFI  = MAX(NRHF(ISYMI),1)
C
              KOFF1  = ISAIKJ(ISYCKL,ISYMI) + 1
              KOFF2  = ISAIKJ(ISYCKL,ISYMJ) + 1
              KOFF3  = IMATIJ(ISYMI,ISYMJ) + 1
C

              CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NCKI(ISYCKL),
     *                ONE,TMAT(KOFF1),NTOCKL,TROCC1(KOFF2),NTOCKL,
     *                ONE,WORK(KOFF3),NRHFI)
C
C

  600      CONTINUE
C
           IF (IPRINT .GT. 55) THEN
              XRMAT = DDOT(NMATIJ(ISYIJ),WORK,1,WORK,1)
              WRITE(LUPRI,*) '3.In CC3_W3_CY2O: Norm of RABIJ =  ',XRMAT 
           ENDIF
C
           CALL CC3_RABIJ(XI2,WORK,ISYMA,A,ISYMB,B,ISYRES)
C
         END IF
C
      END IF
C
      CALL QEXIT('CC3_W3_CY2O')
C
      RETURN
      END
C
C  /* Deck cc3_eta_2 */
      SUBROUTINE CC3_ETA_2(TMAT,ISTMAT,FOCKYCK,ISYFKY,   
     *                 T2TP,ISYMT2,ETA2EFF,ISYRES,RAIJB,
     *                 INDSQ,LENSQ,INDAJLC,ISYAJL,
     *                 ISYMIB,IB,ISYMID,ID,WRK,LWRK)
C 
C     Calculate:
C
C        tbar_mu3  <mu_3|[[X,E_aiE_bj],T_2]|HF> =
C
C         - sum_ckdl ( tbar(aiblck) t(ckdl) X(jd) + tbar(aidjck) t(ckdl) X(lb) )
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C
      IMPLICIT NONE
C
      INTEGER ISTMAT,ISYFKY,ISYMT2,ISYRES,ISYAJL,ISYMIB,ISYMID 
      INTEGER ISYMB,ISYMC,ISYMDKL,ISYMKLJ,ISYRMAT,ISYMJ,ISYMD
      INTEGER ISYMKL,ISYML,ISYMK,ISYMDK,ISYMLJ,ISYMKJ,ISYMBJ
      INTEGER ISYMAI,ISYMDLK,ISYMDC,ISYMKB,ISYMDL,ISYAIJ
      INTEGER IB,ID
      INTEGER NTOTD,NTOTK,NTOTAI,NTOTKL,NTOTB,NTOTAIJ
      INTEGER LENSQ,INDSQ(LENSQ,6),INDAJLC(*)
      INTEGER KOFF,KOFF1,KOFF2,KOFF3
      INTEGER KDKL,KKJL,KKLJ,KTMAT,KRMAT,KEND1,KT2KL,KKB
      INTEGER LWRK,LWRK1

#if defined (SYS_CRAY)
      REAL TMAT(*),FOCKYCK(*),T2TP(*),ETA2EFF(*),RAIJB(*)
      REAL WRK(*)
      REAL ONE,XRMAT,DDOT
#else
      DOUBLE PRECISION TMAT(*),FOCKYCK(*),T2TP(*),ETA2EFF(*),RAIJB(*)
      DOUBLE PRECISION WRK(*)
      DOUBLE PRECISION ONE,XRMAT,DDOT
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('CC3_ETA_2')
C
C--------------------------
C First contribution
C--------------------------
C
C       -  tbar(aiblck) t(ckdl) X(jd) 
C
C       TMAT^BC(aikl) T2TP(dlkC) FOCKYCK(dj)
C
C                     T^C(dkl) FOCKYCK(dj)
C
C                        TF^C(kjl)
C
C       TMAT^BC(aikl) G^C(klj)                 
C
      B = IB
      ISYMB = ISYMIB
      C = ID
      ISYMC = ISYMID
C
      ISYMDKL = MULD2H(ISYMT2,ISYMC)
      ISYMKLJ = MULD2H(ISYMDKL,ISYFKY)
      ISYRMAT = MULD2H(ISYRES,ISYMB)
C
      KDKL    = 1
      KKJL    = KDKL   + NCKI(ISYMDKL)
      KKLJ    = KKJL   + NMAJIK(ISYMKLJ)
      KTMAT   = KKLJ   + NMAJIK(ISYMKLJ)
      KRMAT   = KTMAT  + NCKIJ(ISTMAT)
      KEND1   = KRMAT  + NCKI(ISYRMAT) 
      LWRK1   = LWRK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Not enough space in cc3_eta_2 (1. contribution)')
      END IF
C
      CALL DZERO(WRK(KKJL),NMAJIK(ISYMKLJ))
      CALL DZERO(WRK(KRMAT),NCKI(ISYRMAT))
C
C     T2TP(dlkC) put in WORK(dkl)
C
      KOFF = IT2SP(ISYMDKL,ISYMC) + NCKI(ISYMDKL)*(C - 1) + 1
      CALL CC_GATHER(NCKI(ISYMDKL),WRK(KDKL),T2TP(KOFF),INDAJLC)
C
C         TF^C(kjl) =  T^C(dkl) FOCKYCK(dj)
C
      DO ISYMJ = 1,NSYM
         ISYMD = MULD2H(ISYMJ,ISYFKY)
         ISYMKL = MULD2H(ISYMDKL,ISYMD)
         DO ISYML = 1,NSYM
            ISYMK = MULD2H(ISYMKL,ISYML)
            ISYMDK = MULD2H(ISYMD,ISYMK)
            ISYMKJ = MULD2H(ISYMJ,ISYMK)
            DO L = 1,NRHF(ISYML)
C
               KOFF1 = KDKL + ICKI(ISYMDK,ISYML)
     *                      + NT1AM(ISYMDK)*(L-1)
     *                      + IT1AM(ISYMD,ISYMK)
               KOFF2 = IT1AM(ISYMD,ISYMJ) + 1     
               KOFF3 = KKJL + ISJIK(ISYMKJ,ISYML)
     *                      + NMATIJ(ISYMKJ)*(L-1)
     *                      + IMATIJ(ISYMK,ISYMJ)
C
               NTOTD = MAX(1,NVIR(ISYMD))
               NTOTK = MAX(1,NRHF(ISYMK))
C
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),  
     *                    NVIR(ISYMD),ONE,WRK(KOFF1),NTOTD,
     *                    FOCKYCK(KOFF2),NTOTD,ONE,WRK(KOFF3),
     *                    NTOTK)
            END DO                    
         END DO                       
      END DO
C
C       G^C(klj) = TF^C(kjl)
C
      DO ISYMJ = 1,NSYM
         DO ISYML = 1,NSYM
            ISYMKJ = MULD2H(ISYMKLJ,ISYML)
            ISYMK = MULD2H(ISYMKJ,ISYMJ)
            ISYMKL = MULD2H(ISYMK,ISYML)
            DO J = 1,NRHF(ISYMJ)
               DO L = 1,NRHF(ISYML)
                  DO K = 1,NRHF(ISYMK) 
C
                     KOFF1 = KKJL - 1 
     *                            + ISJIK(ISYMKJ,ISYML)
     *                            + NMATIJ(ISYMKJ)*(L-1)
     *                            + IMATIJ(ISYMK,ISYMJ)
     *                            + NRHF(ISYMK)*(J-1)
     *                            + K 
                     KOFF2 = KKLJ - 1
     *                            + ISJIK(ISYMKL,ISYMJ)
     *                            + NMATIJ(ISYMKL)*(J-1)
     *                            + IMATIJ(ISYMK,ISYML)
     *                            + NRHF(ISYMK)*(L-1)
     *                            + K
C
                    WRK(KOFF2) = WRK(KOFF1)
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL DCOPY(LENSQ,TMAT,1,WRK(KTMAT),1)
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
      IF (NSYM .GT. 1) THEN
        CALL CC_GATHER(LENSQ,WRK(KTMAT),TMAT,INDSQ(1,6))
      ENDIF
C
C----------------------------------
C   ETA2EFF(aibj) =  ETA2EFF(aibj) +
C
C       RAIJB1 = TMAT^BC(aikl) G^C(klj)                 
C----------------------------------
C
      DO  ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES) 
         ISYMKL = MULD2H(ISTMAT,ISYMAI)
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)

         KOFF1  = KTMAT + ISAIKL(ISYMAI,ISYMKL) 
         KOFF2  = KKLJ  + ISJIK(ISYMKL,ISYMJ) 
         KOFF3  = KRMAT + ISAIK(ISYMAI,ISYMJ)
C
         CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *              ONE,WRK(KOFF1),NTOTAI,WRK(KOFF2),NTOTKL,
     *              ONE,WRK(KOFF3),NTOTAI)
C
      END DO            
C
        IF (IPRINT .GT. 55) THEN
           XRMAT = DDOT(NCKI(ISYRMAT),WRK(KRMAT),1,WRK(KRMAT),1)
           WRITE(LUPRI,*) '1. In CC3_ETA_2: Norm of RMAT =  ',XRMAT
        ENDIF
C
        CALL CC3_RACC(ETA2EFF,WRK(KRMAT),ISYMB,B,ISYRES)
C
C----------------------------------
C Second contribution.
C----------------------------------
C
C       TMAT^DC(aikj) T2TP(DlkC) FOCKYCK(bl)
C
C                     T^DC(kl) FOCKYCK(bl)
C
C       WORK^DC(aijk) G^DC(kb)                 
C
      D = IB
      ISYMD = ISYMIB
      C = ID
      ISYMC = ISYMID
C
      ISYMDLK = MULD2H(ISYMT2,ISYMC)
      ISYMDC = MULD2H(ISYMD,ISYMC)
      ISYMKL = MULD2H(ISYMDC,ISYMT2)
      ISYMKB = MULD2H(ISYMKL,ISYFKY)
C
      KT2KL  = 1
      KKB   = KT2KL + NMATIJ(ISYMKL)
      KTMAT = KKB   + NT1AM(ISYMKB)
      KEND1 = KTMAT + NCKIJ(ISTMAT)
      LWRK1 = LWRK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Not enough space in cc3_eta_2 (2. contribution)')
      END IF
C
      CALL DZERO(WRK(KKB),NT1AM(ISYMKB))
C
C -------------------------------
C      T2TP(DlkC)  sorted as T^{DC}_{kl} equal T_{kl}
C -------------------------------
C
      DO ISYML = 1,NSYM
         ISYMK = MULD2H(ISYMKL,ISYML)
         ISYMDL = MULD2H(ISYMDLK,ISYMK)
         DO L = 1,NRHF(ISYML)
            DO K = 1,NRHF(ISYMK)
               KOFF1 = IT2SP(ISYMDLK,ISYMC)
     *                  + NCKI(ISYMDLK)*(C-1)
     *                  + ICKI(ISYMDL,ISYMK)
     *                  + NT1AM(ISYMDL)*(K-1)
     *                  + IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L-1)
     *                  + D
               KOFF2 = IMATIJ(ISYMK,ISYML)
     *                  + NRHF(ISYMK)*(L-1)
     *                  + K
               WRK(KT2KL-1+KOFF2) = T2TP(KOFF1)
            ENDDO
         ENDDO
      ENDDO
C
C         TF^DC(kb) =  T^DC(kl) FOCKYCK(bl)
C
      DO ISYML = 1,NSYM
         ISYMK = MULD2H(ISYMKL,ISYML)
         ISYMB = MULD2H(ISYFKY,ISYML)
C
         KOFF1 = KT2KL + IMATIJ(ISYMK,ISYML)
         KOFF2 = IT1AM(ISYMB,ISYML) + 1     
         KOFF3 = KKB + IT1AMT(ISYMK,ISYMB)
C
         NTOTK = MAX(1,NRHF(ISYMK))
         NTOTB = MAX(1,NVIR(ISYMB))
C
         CALL DGEMM('N','T',NRHF(ISYMK),NVIR(ISYMB),  
     *              NRHF(ISYML),ONE,WRK(KOFF1),NTOTK,
     *              FOCKYCK(KOFF2),NTOTB,ONE,WRK(KOFF3),
     *              NTOTK)
C
      END DO                       
C
C     Sort TMAT^DC(aikj) as WORK^DC(aijk)
C
      DO I = 1, NCKIJ(ISTMAT)
         WRK(KTMAT-1+I) = TMAT(INDSQ(I,3))
      ENDDO
C
C----------------------------------
C
C   RAIJB1 = WORK^DC(aijk) TF^DC(kb)
C
C----------------------------------
C
      DO  ISYMB = 1,NSYM
C
         ISYMK = MULD2H(ISYMB,ISYMKB)
         ISYAIJ = MULD2H(ISTMAT,ISYMK) 
C
         KOFF1  = KTMAT + ISAIKJ(ISYAIJ,ISYMK) 
         KOFF2  = KKB  + IT1AMT(ISYMK,ISYMB) 
         KOFF3  = IT2SP(ISYAIJ,ISYMB) + 1
C
         NTOTAIJ = MAX(NCKI(ISYAIJ),1)
         NTOTK = MAX(NRHF(ISYMK),1)
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NVIR(ISYMB),NRHF(ISYMK),
     *              ONE,WRK(KOFF1),NTOTAIJ,WRK(KOFF2),NTOTK,
     *              ONE,RAIJB(KOFF3),NTOTAIJ)
C
      END DO            
C
        IF (IPRINT .GT. 55) THEN
           XRMAT = DDOT(NT2SQ(ISYRES),RAIJB,1,RAIJB,1)
           WRITE(LUPRI,*) '2. In CC3_ETA_2: Norm of RAIJB =  ',XRMAT
        ENDIF
C
      CALL QEXIT('CC3_ETA_2')
C
      RETURN
      END
C
C  /* Deck cc3_eta_1 */
      SUBROUTINE CC3_ETA_1(FOCKYCK,ISYFKY,ETA1EFF,ISYRES,
     *                     WRK,LWRK)
C
C---------------------------------------------------------------------
C
C     Calculate:
C
C     <L3|[[X,tau_ai],T_3]|HF> =
C
C       sum_l X(la) D^0(li) - sum_d X(id) D^0(ad)
C
C      where
C      D^0(pq) = <L3|[E_pq,T_3]|HF>  ; p,q = a,b or i,j
C
C
      IMPLICIT NONE
C
      CHARACTER FNDPTIJ*5, FNDPTAB*5
      PARAMETER(FNDPTIJ = 'DPTIJ', FNDPTAB='DPTAB')
C
      INTEGER ISYFKY,ISYRES,ISYM0,ISYMI,ISYMA,ISYML,ISYMD
      INTEGER KD0AB,KD0IJ,KEND1
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER NTOTL,NTOTA,NTOTD
      INTEGER LWRK,LWRK1
      INTEGER LUPTIJ, LUPTAB
C
#if defined (SYS_CRAY)
      REAL FOCKYCK(*),ETA1EFF(*),WRK(*),ONE
#else
      DOUBLE PRECISION FOCKYCK(*),ETA1EFF(*),WRK(*),ONE
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      PARAMETER (ONE = 1.0D0)
C 
      CALL QENTER('CC3_ETA_1')
C
      ISYM0 = 1
C
      KD0AB  = 1
      KD0IJ  = KD0AB  + NMATAB(ISYM0)
      KEND1  = KD0IJ  + NMATIJ(ISYM0)
      LWRK1  = LWRK   - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Not enough space in cc3_eta_1')
      END IF
C
      LUPTIJ = -1
      LUPTAB = -1
C
      ! read intermediates from ground state density from file...
C
      CALL WOPEN2(LUPTIJ,FNDPTIJ,64,0)
      CALL GETWA2(LUPTIJ,FNDPTIJ,WRK(KD0IJ),1,NMATIJ(ISYM0))
      CALL WCLOSE2(LUPTIJ,FNDPTIJ,'KEEP')

      CALL WOPEN2(LUPTAB,FNDPTAB,64,0)
      CALL GETWA2(LUPTAB,FNDPTAB,WRK(KD0AB),1,NMATAB(ISYM0))
      CALL WCLOSE2(LUPTAB,FNDPTAB,'KEEP')
C
C     ----------------------------------------------------------------
C       sum_l     X(la)      D^0(li) 
C             FOCKYCK(al) 
C     ----------------------------------------------------------------
C
      DO ISYMI = 1, NSYM
         ISYMA = MULD2H(ISYRES,ISYMI)
         ISYML = MULD2H(ISYFKY,ISYMA)
C
         KOFF1   = IT1AM(ISYMA,ISYML) + 1
         KOFF2   = KD0IJ + IMATIJ(ISYML,ISYMI)
         KOFF3   = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTL   = MAX(NRHF(ISYML),1)
         NTOTA   = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYML),
     *               ONE,FOCKYCK(KOFF1),NTOTA,WRK(KOFF2),NTOTL,
     *               ONE,ETA1EFF(KOFF3),NTOTA)
C
      END DO
C
C
C     ----------------------------------------------------------------
C     - sum_d     D^0(ad)    X(id)      
C                          FOCKYCK(di) 
C     ----------------------------------------------------------------
C
      DO ISYMI = 1, NSYM
         ISYMA = MULD2H(ISYRES,ISYMI)
         ISYMD = MULD2H(ISYFKY,ISYMI)
C 
         KOFF1   = KD0AB + IMATAB(ISYMA,ISYMD)
         KOFF2   = IT1AM(ISYMD,ISYMI) + 1
         KOFF3   = IT1AM(ISYMA,ISYMI) + 1
C 
         NTOTA   = MAX(NVIR(ISYMA),1)
         NTOTD   = MAX(NVIR(ISYMD),1)
C        
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMD),
     *               -ONE,WRK(KOFF1),NTOTA,FOCKYCK(KOFF2),NTOTD,
     *               ONE,ETA1EFF(KOFF3),NTOTA)
C
      END DO
C
      CALL QEXIT('CC3_ETA_1')
C
      RETURN
      END
C

