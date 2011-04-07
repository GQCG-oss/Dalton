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
      SUBROUTINE CC3_FMATSD(LISTL,IDLSTL,LISTY,IDLSTY,
     *                     OMEGA1,OMEGA2,OMEGA1EFF,
     *                     OMEGA2EFF,ISYRES,WORK,LWORK,LUDKBC,FNDKBC,
     *                     LUCKJD,FNCKJD, 
     *                     LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                     LU3FOP,FN3FOP,LU3FOP2,FN3FOP2)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples component of F*T^Y vector 
*             projected into the singles and doubles space
*
*             Use W intermmediate as these are required for 
*             ETA3 projection and therefore are needed to 
*             construct t3bar^Y(omega_Y).
*             
*    F*T^Y3 =  <L2|[H^Y,tau3]|HF>/ (w+epsiln(tau3))
*
*    F*T^Y1(ai) = F*T^Y1(ai) + < F*T^Y3|[[H^,T2],tau(ai)]|HF>     
*
*    F*T^Y2(aibj) = F*T^Y2(aibj) + < F*T^Y3|[H^,tau(aibj)]|HF>     
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
#include <ccr1rsp.h>
#include <ccexci.h>
C

      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)

      INTEGER   ISYRES,ISYCTR,ISYMY,ISINT1,ISINT2,ISINT1Y,ISINT2Y
      INTEGER   ISYFCKY,ISYMD,ISYCKB,ISCKB1,ISCKB2,ISCKB2Y,ISYMB
      INTEGER   ISYMBD,ISWMAT,ISYALJB,ISYALJC,ISYMIM
      INTEGER   ISYMT2
      INTEGER   IDLSTL,IDLSTY
      INTEGER   KRBJIA,KLAMP0,KLAMH0,KFOCKD,KT2TP,KL2TP
      INTEGER   KEND0,KT1AMP0,KEND1,KOIOOOO,KOIOVVO,KOIOOVV
      INTEGER   KTROC,KTROC1,KXIAJB,KTROCCG,KTROCCL,KFCKYCK
      INTEGER   KLAMPY,KLAMHY,KINTOC,KT1Y,KVVVV,KTRVI,KTRVI1
      INTEGER   KVGDKCB,KVGDKBC,KVLDKCB,KVLDKBC,KINTVI
      INTEGER   KDIAGW,KINDSQW,KINDEXB,KINDEXC,KTMAT,KWMAT
      INTEGER   LWORK,LWRK0,LWRK1,LENGTH,LENSQW
      INTEGER   IOPT
      INTEGER   IOFF
      INTEGER   LU4V
      INTEGER   LUCKJD, LUTOC, LU3VI, LU3VI2, LU3FOP, LU3FOP2
      INTEGER   LUDKBC, LUDKBC4
C
      INTEGER ISYMN1,ISYMN2,KN2MAT,KINDSQN,LENSQN,ISGEI,ISFEI,KGEI,KFEI
      INTEGER IADR
      INTEGER KEND2,LWRK2
C
      INTEGER kx3am
C

      CHARACTER LISTL*3,LISTY*3,MODEL*10
      CHARACTER*(*) FNCKJD, FNTOC, FN3VI, FN3VI2, FN3FOP, FN3FOP2 
      CHARACTER*(*)  FNDKBC
      CHARACTER*14 FN4V
      CHARACTER*11 FNDKBC4
      CHARACTER*1 CDUMMY
C
      CHARACTER*6 FNGEI,FNFEI
      CHARACTER*5 FNN1
      PARAMETER( FNGEI = 'N1_GEI' , FNFEI = 'N1_FEI' , FNN1 = 'N1MAT' )
      INTEGER LUGEI,LUFEI,LUN1
C
      LOGICAL   LOCDBG

#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      REAL WORK(*)
      REAL FREQY,FREQL,FREQYML
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION  OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      DOUBLE PRECISION  WORK(*)
      DOUBLE PRECISION  FREQY,FREQL,FREQYML
      DOUBLE PRECISION  XNORMVAL,DDOT
#endif
C
      INTEGER ILSTSYM
C
      CALL QENTER('CC3FMATSD')
C
C----------------------------------------------------
C     Initialise character strings and open files
C----------------------------------------------------
C
      CDUMMY = ' '
C
      LU4V     = -1
C
      FN4V     = 'CC3_FMATSD_TMP'
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
        FREQL  = 0.0D0
      ELSE IF (LISTL(1:3).EQ.'LE ') THEN
        ISYCTR = ILSTSYM(LISTL,IDLSTL)
        FREQL  = EIGVAL(IDLSTL)
      ELSE
        ! ups, probably higher-order response, not yet implemented here
        CALL QUIT('Unknown type of left vector in CC3_FMATSD.')
      END IF

      ISYMT2 = ISYCTR

      IF (LISTY(1:3).EQ.'R1 ') THEN

        ! get symmetry, frequency and integral label from common blocks
        ! defined in ccl1rsp.h
        ISYMY  = ISYLRT(IDLSTY)
        FREQY  = FRQLRT(IDLSTY)
      ELSE IF (LISTY(1:3).EQ.'RE ') THEN
C
        ISYMY  = ILSTSYM(LISTY,IDLSTY)
        FREQY  = EIGVAL(IDLSTY)
C
      ELSE
        CALL QUIT('Unknown/Illegal list in CC3_FMATSD')
      END IF
C
      !frequency asscociated with the L1 and M1 and N2, respectivly.
      FREQYML = FREQY - FREQL
C
C-------------------------------------------------------
C     initial allocations, orbital energy, fock matrix and T2 and L2 :
C-------------------------------------------------------
C
C     Symmetry of integrals in contraction:
C
      ISINT1 = 1
      ISINT2 = 1
      ISYFCKY = MULD2H(ISYMOP,ISYMY)
      ISYMIM = MULD2H(ISYFCKY,ISYMT2)
      IF (.NOT.LVVVV) THEN
         !Symmetries for N1 and N2 intermediates
         ISYMN1 = MULD2H(ISYMIM,ISYM0)
         ISYMN2 = MULD2H(ISYMIM,ISYM0)
      END IF
C
      IF (LVVVV) THEN
         KRBJIA = 1
      ELSE
         KN2MAT = 1
         KRBJIA = KN2MAT + NCKIJ(ISYMN2)
      END IF
      KLAMP0 = KRBJIA   + NT2SQ(ISYRES)
      KLAMH0  = KLAMP0  + NLAMDT
      KFOCKD  = KLAMH0  + NLAMDT
      KT2TP  = KFOCKD  + NORBTS
      KL2TP   = KT2TP   + NT2SQ(ISYM0)
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
C----------------------------------------------------------------------------
C     With option IOPT = 2 it just calculates LUDKBC4i needed in contractions.
C----------------------------------------------------------------------------
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

C     IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
C    *    DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)
C
C-------------------------------------
C     Read L2 amplitudes 
C-------------------------------------
C
      IOPT = 2
      CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     *              DUMMY,WORK(KL2TP))

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
            call quit('Insufficient space (T3) in CC3_FMATSD')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT

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
      ISINT1Y = MULD2H(ISYMY,ISINT1)
      ISINT2Y = MULD2H(ISYMY,ISINT2)

      KTROC    = KEND0
      KTROC1   = KTROC    + NTRAOC(ISINT2)
      KXIAJB   = KTROC1   + NTRAOC(ISINT2)
      KTROCCG   = KXIAJB   + NT2AM(ISYM0)
      KTROCCL  = KTROCCG   + NTRAOC(ISINT2Y)
      KEND0    = KTROCCL  + NTRAOC(ISINT2Y)

      KFCKYCK = KEND0
      KLAMPY  = KFCKYCK + NT1AM(ISYFCKY)
      KLAMHY  = KLAMPY  + NLAMDT
      KEND0   = KLAMHY  + NLAMDT
      LWRK0   = LWORK   - KEND0

      KINTOC = KEND0
      KT1Y   = KINTOC + MAX(NTOTOC(ISINT2Y),NTOTOC(ISINT2))
      KEND1  = KT1Y + NT1AM(ISYMY)
      LWRK1  = LWORK  - KEND1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_FMATSD')
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
C---------------------------------------------
C        Regain work space and create lam^C
C---------------------------------------------
C

C        isint1 - symmetry of integrals in standard H

      ISINT1  = 1
      ISINT1Y = MULD2H(ISINT1,ISYMY)
C

      IOPT = 1
C
      CALL CC_RDRSP(LISTY,IDLSTY,ISYMY,IOPT,MODEL,
     *              WORK(KT1Y),DUMMY)
C
      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPY),WORK(KLAMH0),
     *                 WORK(KLAMHY),WORK(KT1Y),ISYMY)
C
C------------------------------------------------------------------
C        Calculate the F^Y matrix (kc elements evaluated and stored 
C        as ck)
C------------------------------------------------------------------
C
      CALL CC3LR_MFOCK(WORK(KFCKYCK),WORK(KT1Y),WORK(KXIAJB),ISYFCKY)
C
C------------------------
C     Occupied integrals.
C
C     Read in integrals for WMAT etc.
C-----------------------
C
      IOFF = 1
      IF (NTOTOC(ISYMOP) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
      ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XNORMVAL
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k-) and sort  G(j,i,k-,a)
C----------------------------------------------------------------------
C
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROCCG),
     *                 WORK(KLAMHY),ISYMY,
     *                 WORK(KEND1),LWRK1)
C
C----------------------------------------------------------------------
C integrals g(ia|j k-) sorted as G(j,i,k-,a) resorted as G(k-,i,j,a)
C                       on input KTROCCG   on output     KTROCCG
C           L(ia|j k-)                                   KTROCCL
C----------------------------------------------------------------------
C
      CALL CC3_SRTOCC(WORK(KTROCCG),WORK(KTROCCL),ISINT1Y)
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
     *                    WORK(KEND1),LWRK1,ISINT2)
C
      CALL CCFOP_SORT(WORK(KTROC),WORK(KTROC1),ISINT2,1)
C
      CALL CC3_LSORT1(WORK(KTROC),ISINT2,WORK(KEND1),LWRK1,5)

C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKB  = MULD2H(ISYMD,ISYM0)
         ISCKB1 = MULD2H(ISINT1,ISYMD)
         ISCKB2 = MULD2H(ISYMD,ISYM0)
         ISCKB2Y = MULD2H(ISYMD,ISINT2Y)
C
         IF (.NOT.LVVVV) THEN
            !Symmetry of arrays needed to construct N1MAT
            ISGEI  = MULD2H(ISYMN1,ISYMD)
            ISFEI  = MULD2H(ISYMN1,ISYMD)
         END IF
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
         KVGDKCB  = KEND1
         KVGDKBC  = KVGDKCB  + NCKATR(ISCKB2Y)
         KVLDKCB  = KVGDKBC  + NCKATR(ISCKB2Y)
         KVLDKBC  = KVLDKCB  + NCKATR(ISCKB2Y)
         KEND1    = KVLDKBC  + NCKATR(ISCKB2Y)
         LWRK1    = LWORK    - KEND1

         KINTVI   = KEND1
         KEND1    = KINTVI + NCKATR(ISCKB2Y)
         LWRK1    = LWORK  - KEND1
C
         IF (LWRK1 .LT. NCKATR(ISCKB2Y)) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK1
            WRITE(LUPRI,*) 'Memory needed    : ',NCKATR(ISCKB2Y)
            CALL QUIT('Insufficient space in CC3_FMATSD')
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

C
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
C           Read in, transform and sort integrals used in 
C           construction of W intermediate.
C------------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
C
C           Read in g(c1^h k1^p | delta b2^h) integrals
C
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF
C
C          Transform g(c1^h k1^p | delta b2^h) integrals
C          to g(c1^h k1^p | d2^pY b2^h) 
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKCB),
     *                      WORK(KLAMPY),ISYMY,
     *                      ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C
C          Sort g(c1^h k1^p | d2^pY b2^h) as 
C          g(d2^pY k1^p | c1^h b2^h)
C
            CALL CCSDT_SRVIR3(WORK(KVGDKCB),WORK(KEND1),ISYMD,D,ISYMY)

C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
C
C           Read in g(b2^h k1^p | delta c1^h) integrals and 
C           transform and sort --- see above
C
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF

            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKBC),
     *                          WORK(KLAMPY),ISYMY,
     *                          ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C

            CALL CCSDT_SRVIR3(WORK(KVGDKBC),WORK(KEND1),ISYMD,D,ISYMY)

C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
C
C           Read in L(c1^h k1^p | delta b2^h) integrals
C
               CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF
C
C          Transform L(c1^h k1^p | delta b2^h) integrals
C          to L(c1^h k1^p | d2^pY b2^h) 
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKCB),
     *                          WORK(KLAMPY),ISYMY,
     *                          ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C
C          Sort L(c1^h k1^p | d2^pY b2^h) as 
C          L(d2^pY k1^p | c1^h b2^h)
C
            CALL CCSDT_SRVIR3(WORK(KVLDKCB),WORK(KEND1),ISYMD,D,ISYMY)
C
C           Read in L(b2^h k1^p | delta c1^h) integrals and 
C           transform and sort --- see above
C
            CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKBC),
     *                          WORK(KLAMPY),ISYMY,
     *                          ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C
            CALL CCSDT_SRVIR3(WORK(KVLDKBC),WORK(KEND1),ISYMD,D,ISYMY)
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
     &                   'CC3_FMATSD (TRVI)')
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
     &                   'CC3_FMATSD (TRVI1)')
            END IF
C
            CALL CCSDT_SRVIR3(WORK(KTRVI1),WORK(KEND1),ISYMD,D,ISINT1)
C

C
            DO ISYMB = 1,NSYM

               ISYMBD  = MULD2H(ISYMD,ISYMB)
               ISWMAT  = MULD2H(ISYMIM,ISYMBD)
               ISYALJB = MULD2H(ISYMT2,ISYMB)
               ISYALJC = MULD2H(ISYMT2,ISYMD)

C              Can use kend3 since we do not need the integrals anymore.
               KDIAGW   = KEND1
               KINDSQW  = KDIAGW  + NCKIJ(ISWMAT)
               KINDEXB  = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEXC  = KINDEXB + (NCKI(ISYALJB) - 1)/IRAT + 1
               KTMAT    = KINDEXC + (NCKI(ISYALJC) - 1)/IRAT + 1
               KWMAT    = KTMAT   + NCKIJ(ISWMAT)
               KEND2    = KWMAT   + NCKIJ(ISWMAT)
               LWRK2    = LWORK   - KEND2

               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND2
                  CALL QUIT('Insufficient space in CC3_FMATSD')
               END IF
C
C
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
C
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)

C              IF (LOCDBG) THEN
C                WRITE(LUPRI,*) 'Norm of DIA  ',
C    *              DDOT(NCKIJ(ISWMAT),WORK(KDIAGW),1,WORK(KDIAGW),1)
C              END IF
C
C              -----------------------
C              Construct index arrays.
C              -----------------------
C
               LENSQW  = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT)
               CALL CC3_INDEX(WORK(KINDEXB),ISYALJB)
               CALL CC3_INDEX(WORK(KINDEXC),ISYALJC)


               DO B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
C
C    calculate     <L2|[H^Y,tau3]|HF>
C
                  CALL WBARBD_TMAT(WORK(KL2TP),ISYCTR,WORK(KWMAT),
     *                             WORK(KTMAT),ISWMAT,WORK(KFCKYCK),
     *                             ISYFCKY,WORK(KVLDKBC),WORK(KVLDKCB),
     *                             WORK(KVGDKBC),WORK(KVGDKCB),
     *                             WORK(KTROCCL),WORK(KTROCCG),
     *                             ISINT2Y,WORK(KEND2),LWRK2,
     *                             WORK(KINDEXB),WORK(KINDEXC),
     *                             WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)
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
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,-FREQYML,
     *                         ISWMAT,WORK(KWMAT),
     *                         WORK(KDIAGW),WORK(KFOCKD))
C
                  CALL T3_FORBIDDEN(WORK(KWMAT),ISYMIM,ISYMB,B,ISYMD,D)
C
C
                  CALL CC3_W3_CY2V(OMEGA2EFF,ISYRES,WORK(KRBJIA),
     *                             WORK(KWMAT),ISWMAT,
     *                             WORK(KTMAT),WORK(KTRVI),WORK(KTRVI1),
     *                             ISINT1,WORK(KEND2),LWRK2,
     *                             WORK(KINDSQW),LENSQW,
     *                             ISYMB,B,ISYMD,D,.TRUE.)
C
                  CALL CC3_W3_CY2O(OMEGA2EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                             WORK(KTMAT),WORK(KTROC),WORK(KTROC1),
     *                             ISINT1,WORK(KEND2),LWRK2,
     *                             WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D,
     *                             .TRUE.)
C
                  IF (LVVVV) THEN
                    CALL CC3_W3_OMEGA1(OMEGA1EFF,ISYRES,WORK(KWMAT),
     *                                 WORK(KTMAT),ISYMIM,
     *                                 WORK(KOIOOOO),WORK(KOIOVVO),
     *                                 WORK(KOIOOVV),WORK(KVVVV),ISYM0,
     *                                 WORK(KT2TP),ISYM0,
     *                                 WORK(KEND2),LWRK2,
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
     *                       WORK(KEND2),LWRK2,
     *                       .TRUE.)

                   END IF


               ENDDO   ! B
            ENDDO      ! ISYMB
C
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
C
COMMENT COMMENT
C      write(lupri,*) 'WMAT in CC3_FMATSD : '
C      call print_pt3(work(kx3am),ISYFCKY,4)
COMMENT COMMENT
C


C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Vccupied  cont ) 
C     in XI2EFF 
C------------------------------------------------------
C
      CALL CC3_RBJIA(OMEGA2EFF,ISYRES,WORK(KRBJIA))
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
     *                     OMEGA1EFF,ISYRES,
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
     *                OMEGA1EFF,ISYRES,
     *                WORK(KEND0),LWRK0)
C
      END IF

C
C
C       write(lupri,*) 'OMEGA1EFF inside cc3_fmatsd'
C       do i = 1,nt1am(isyres)
C           write(lupri,*) i,OMEGA1EFF(i)
C       end do
C
C     write(lupri,*) 'OMEGA2EFF inside cc3_fmatsd'
C     XNORMVAL = ddot(nt2am(isyres),OMEGA2EFF,1,OMEGA2EFF,1)
C     write(lupri,*) 'OMEGA2EFF Norm =', XNORMVAL
C     call outpak(OMEGA2EFF,NT1AM(isyres),1,lupri)
C
C
      DO I = 1,NT1AM(ISYRES)
         OMEGA1EFF(I) = OMEGA1EFF(I) + OMEGA1(I)
      END DO
C
      DO I = 1,NT2AM(ISYRES)
         OMEGA2EFF(I) = OMEGA2EFF(I) + OMEGA2(I)
      END DO
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
C-------------
C     End
C-------------
C
C
      CALL QEXIT('CC3FMATSD')
C
      RETURN
      END



C  /* Deck wbarbd_tmat */
      SUBROUTINE WBARBD_TMAT(T2TP,ISYMT2,WMAT,TMAT,ISWMAT,FOCK,
     *                      ISYFOCK,VLDKBC,VLDKCB,VGDKBC,VGDKCB,TROCCL,
     *                      TROCCG,ISINT2,WORK,LWORK,INDAJLB,
     *                      INDAJLC,INDSQ,LENSQ,ISYMB,B,ISYMC,C)
C
C     Written by Kasper Hald, Fall 2001.
C
C     General symmetry: ISINT2 is symmetry of integrals 
C                       ISYMT2 is symmetry of T2TP
C
C     Virtual integrals stored as:
C          L(kcd^b) -> IC(d^kB):  VLDKBC
C          L(kcd^b) -> IB(d^kC):  VLDKCB
C          g(kcd^b) -> IC(d^kB):  VGDKBC
C          g(kcd^b) -> IB(d^kC):  VGDKCB

C     Occupied integrals stored as:
C          L(ia|j k-) -> I(k-,i,j,a): TROCCL 
C          g(ia|j k-) -> I(k-,i,j,a): TROCCG 
C
C
C     WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijB)*F(kC) 
C                                   - T2TP(aikB)*F(jC)
C                                   + T2TP(aikC)*F(jB)
C                                   - T2TP(aijC)*F(kB)
C
C                   + T2TP(aijd)*L(d^BkC)
C                   - T2TP(ajkd)*g(iBd^C)
C                   + T2TP(aikd)*L(d^CjB)
C                   - T2TP(akjd)*g(iCd^B)
C
C                   + T2TP(ailB)*L(jl^kC)
C                   - T2TP(alkB)*g(il^jC)
C                   + T2TP(ailC)*L(kl^jB)
C                   - T2TP(aljC)*g(il^kB)
C

      IMPLICIT NONE
C
      INTEGER ISYMT2,ISWMAT,ISINT2,ISYMB,ISYMC,ISYRES,ISYMBC
      INTEGER JSAIKJ,ISYMK,ISYAIJ,ISYMJ,ISYMAI,ISYAIK,ISYMI
      INTEGER ISYMA,ISYMDK,ISYMD,ISYMDI,ISYAJK,ISYMDJ,ISYAIL
      INTEGER ISYLKJ,ISYMLK,ISYML,ISYALK,ISYLJI,ISYAKJ,ISYMLJ
      INTEGER ISYMAK,ISYAJL,ISYLKI,ISYMAJ,ISYFOCK
      INTEGER NAI,NAIJB,NCK,NAIKJ,NCJ,NAIKB,NBJ,NAIKC,NAIJC,NBK
      INTEGER NTOAIJ,NVIRD,NTOAJK,NTOAIK,NTOAKJ,NTOTAI,NRHFL
      INTEGER NTOTAK,NTOTAJ
      INTEGER INDAJLB,INDAJLC,LENSQ,INDSQ(LENSQ,6),INDEX
      INTEGER KOFF1,KOFF2,KOFF3,KALK,KEND1,KALJ,KOFF
      INTEGER LWORK,LENGTH,LWRK1

#if defined (SYS_CRAY)
      REAL T2TP(*),WMAT(*),TMAT(*),FOCK(*)
      REAL VLDKBC(*),VLDKCB(*),VGDKBC(*),VGDKCB(*),TROCCL(*),TROCCG(*)
      REAL WORK(*)
      REAL XWMAT,ONE,DDOT
#else
      DOUBLE PRECISION T2TP(*),WMAT(*),TMAT(*),FOCK(*)
      DOUBLE PRECISION VLDKBC(*),VLDKCB(*),VGDKBC(*),VGDKCB(*),TROCCL(*)
      DOUBLE PRECISION TROCCG(*),WORK(*)
      DOUBLE PRECISION XWMAT,ONE,DDOT
#endif
C
      PARAMETER(ONE = 1.0D0)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('WBARBD_TMAT')
C
      ISYRES = MULD2H(ISYMT2,ISINT2)
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYRES,ISYMBC)
      LENGTH = NCKIJ(JSAIKJ)
C
C-----------------------------------------------------------------------
C    
C     WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijB)*F(kC) 
C                                   - T2TP(aikB)*F(jC)
C                                   + T2TP(aikC)*F(jB)
C                                   - T2TP(aijC)*F(kB)
C
C-----------------------------------------------------------------------
C     Contribution from both Fock terms:
C-----------------------------------------------------------------------
C
C
C     WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijB)*F(kC) 
C
      ISYMK  = MULD2H(ISYFOCK,ISYMC)
      ISYAIJ = MULD2H(ISYMT2,ISYMB)
C
      DO ISYMJ = 1, NSYM
         ISYMAI = MULD2H(ISYAIJ,ISYMJ)
         ISYAIK = MULD2H(ISYMK,ISYMAI)
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
C
            DO J = 1, NRHF(ISYMJ)
C
               DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                 Index for sorted T2 amplitudes.
C
                  NAIJB = IT2SP(ISYAIJ,ISYMB)
     *                  + NCKI(ISYAIJ)*(B - 1)
     *                  + ICKI(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + NAI
C
                  DO K = 1, NRHF(ISYMK)
C
                     NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) + C 
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K-1)
     *                     + NAI
                            
C
C Fock 1.0 contribution addWMAT
C

                     WMAT(NAIKJ) = WMAT(NAIKJ) + T2TP(NAIJB)*FOCK(NCK)
C
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(aikB)*F(jC)
C
      ISYMJ  = MULD2H(ISYFOCK,ISYMC)
      ISYAIK = MULD2H(ISYMT2,ISYMB)
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
C
            DO J = 1, NRHF(ISYMJ)
               NCJ = IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J-1) + C
C
               DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                 Index for sorted T2 amplitudes.
C
                  DO K = 1, NRHF(ISYMK)
C
                     NAIKB = IT2SP(ISYAIK,ISYMB)
     *                     + NCKI(ISYAIK)*(B - 1)
     *                     + ICKI(ISYMAI,ISYMK) 
     *                     + NT1AM(ISYMAI)*(K - 1) + NAI
C
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K-1)
     *                     + NAI

C
C Fock 2.0 contribution addWMAT
C
                     WMAT(NAIKJ) = WMAT(NAIKJ) - T2TP(NAIKB)*FOCK(NCJ)
C
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C
C     WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aikC)*F(jB)
C                                  
      ISYMJ  = MULD2H(ISYFOCK,ISYMB)
      ISYAIK = MULD2H(ISYMT2,ISYMC)
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
C
            DO J = 1, NRHF(ISYMJ)
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
               DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
C 
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                 Index for sorted T2 amplitudes.
C 
                  DO K = 1, NRHF(ISYMK)
C                    
                     NAIKC = IT2SP(ISYAIK,ISYMC)
     *                     + NCKI(ISYAIK)*(C - 1)
     *                     + ICKI(ISYMAI,ISYMK) 
     *                     + NT1AM(ISYMAI)*(K - 1) + NAI
C
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K-1)
     *                     + NAI

C                    
C Fock 3.0 contribution addWMAT
C
                     WMAT(NAIKJ) = WMAT(NAIKJ) + T2TP(NAIKC)*FOCK(NBJ)
C                 
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C
C     WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(aijC)*F(kB)
C
      ISYMK  = MULD2H(ISYFOCK,ISYMB)
      ISYAIJ = MULD2H(ISYMT2,ISYMC)
C
      DO ISYMJ = 1, NSYM
         ISYMAI = MULD2H(ISYAIJ,ISYMJ)
         ISYAIK = MULD2H(ISYMK,ISYMAI)
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
C
            DO J = 1, NRHF(ISYMJ)
C
               DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                 Index for sorted T2 amplitudes.
C
                  NAIJC = IT2SP(ISYAIJ,ISYMC)
     *                  + NCKI(ISYAIJ)*(C - 1)
     *                  + ICKI(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + NAI
C
                  DO K = 1, NRHF(ISYMK)
C
                     NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K-1) + B
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K-1)
     *                     + NAI

C
C Fock 4.0 contribution addWMAT
C
                     WMAT(NAIKJ) = WMAT(NAIKJ) - T2TP(NAIJC)*FOCK(NBK)
C
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------------------------------------------
C     First virtual contribution of L term.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijd)*L(d^BkC)
C---------------------------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
      ISYMDK = MULD2H(ISYMBC,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      CALL DZERO(TMAT,LENGTH)
C
      DO ISYMK = 1,NSYM
C
         ISYMD  = MULD2H(ISYMK,ISYMDK)
         ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIJ,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDK,ISYMB) + NT1AM(ISYMDK)*(B - 1)
     *         + IT1AM(ISYMD,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
         NTOAIJ = MAX(1,NCKI(ISYAIJ))
         NVIRD  = MAX(NVIR(ISYMD),1)
C                    
C Virtual-L 1.0 contribution addWMAT
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *              NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIJ,
     *              VLDKBC(KOFF2),NVIRD,ONE,
     *              TMAT(KOFF3),NTOAIJ)
C
      ENDDO
C
C      CALL CC_GATHER(LENGTH,SMAT,WORK,INDSQ(1,3))
      DO I = 1,LENGTH
         WMAT(I) = WMAT(I) + TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: vir-L 1. Norm of WMAT ',XWMAT
      ENDIF


C---------------------------------------------------------
C     First virtual contribution of g term.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(ajkd)*g(iBd^C)
C---------------------------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
      ISYMDI = MULD2H(ISYMBC,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      CALL DZERO(TMAT,LENGTH)
C
      DO ISYMI = 1,NSYM
C
         ISYMD  = MULD2H(ISYMI,ISYMDI)
         ISYAJK = MULD2H(ISYMI,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAJK,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDI,ISYMB) + NT1AM(ISYMDI)*(B - 1)
     *         + IT1AM(ISYMD,ISYMI)   + 1
         KOFF3 = ISAIKJ(ISYAJK,ISYMI) + 1
C
         NTOAJK = MAX(1,NCKI(ISYAJK))
         NVIRD  = MAX(NVIR(ISYMD),1)
C
C Virtual-g 1.0 contribution addWMAT
C
         CALL DGEMM('N','N',NCKI(ISYAJK),NRHF(ISYMI),
     *              NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAJK,
     *              VGDKCB(KOFF2),NVIRD,ONE,
     *              TMAT(KOFF3),NTOAJK)
C
      ENDDO
C     
      DO I = 1,LENGTH
         WMAT(I) = WMAT(I) - TMAT(INDSQ(I,5))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: vir-g 1. Norm of WMAT ',XWMAT
      ENDIF
C
C---------------------------------------------------------
C     Second virtual contribution of L term.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aikd)*L(d^CjB)
C---------------------------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
      ISYMDJ = MULD2H(ISYMBC,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      CALL DZERO(TMAT,LENGTH)
C
      DO ISYMJ = 1,NSYM
C
         ISYMD  = MULD2H(ISYMJ,ISYMDJ)
         ISYAIK = MULD2H(ISYMJ,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIK,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDJ,ISYMB) + NT1AM(ISYMDJ)*(B - 1)
     *         + IT1AM(ISYMD,ISYMJ)   + 1
         KOFF3 = ISAIKJ(ISYAIK,ISYMJ) + 1
C
         NTOAIK = MAX(1,NCKI(ISYAIK))
         NVIRD  = MAX(NVIR(ISYMD),1)
C
C Virtual-L 2.0 contribution addWMAT
C
         CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *              NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIK,
     *              VLDKCB(KOFF2),NVIRD,ONE,
     *              TMAT(KOFF3),NTOAIK)
C
      ENDDO
C
      DO I = 1,LENGTH
         WMAT(I) = WMAT(I) + TMAT(I)
      ENDDO
C     
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: vir-L 2. Norm of WMAT ',XWMAT
      ENDIF
C
C---------------------------------------------------------
C     Second virtual contribution of g term.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(akjd)*g(iCd^B) 
C---------------------------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
      ISYMDI = MULD2H(ISYMBC,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      CALL DZERO(TMAT,LENGTH)
C
      DO ISYMI = 1,NSYM
C
         ISYMD  = MULD2H(ISYMI,ISYMDI)
         ISYAKJ = MULD2H(ISYMI,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAKJ,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDI,ISYMB) + NT1AM(ISYMDI)*(B - 1)
     *         + IT1AM(ISYMD,ISYMI)   + 1
         KOFF3 = ISAIKJ(ISYAKJ,ISYMI) + 1
C
         NTOAKJ = MAX(1,NCKI(ISYAKJ))
         NVIRD  = MAX(NVIR(ISYMD),1)
C
C Virtual-g 2.0 contribution addWMAT
C
         CALL DGEMM('N','N',NCKI(ISYAKJ),NRHF(ISYMI),
     *              NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAKJ,
     *              VGDKBC(KOFF2),NVIRD,ONE,
     *              TMAT(KOFF3),NTOAKJ)
C
      ENDDO
C
C     
      DO I = 1,LENGTH
         WMAT(I) = WMAT(I) - TMAT(INDSQ(I,2))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: vir-g 2. Norm of WMAT ',XWMAT
      ENDIF
C
C----------------------------------------
C     First occupied L contribution.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj)  
C                    + T2TP(ailB)*L(jl^kC) 
C                                      
C                      TB(ail)*LC(l^kj) = R(aikj)
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
      ISYAIL = MULD2H(ISYMB,ISYMT2)
      ISYLKJ = MULD2H(ISYMC,ISINT2)
C
      CALL DZERO(TMAT,LENGTH)
C
      DO ISYMJ = 1,NSYM
C
         ISYMLK = MULD2H(ISYMJ,ISYLKJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMK = 1,NSYM
C
               ISYML  = MULD2H(ISYMK,ISYMLK)
               ISYMAI = MULD2H(ISYAIL,ISYML)
               ISYAIK = MULD2H(ISYMAI,ISYMK)
C
               KOFF1 = IT2SP(ISYAIL,ISYMB)
     *               + NCKI(ISYAIL)*(B - 1)
     *               + ICKI(ISYMAI,ISYML) + 1
               KOFF2 = ISJIKA(ISYLKJ,ISYMC)
     *               + NMAJIK(ISYLKJ)*(C - 1)
     *               + ISJIK(ISYMLK,ISYMJ)
     *               + NMATIJ(ISYMLK)*(J - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
               KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *               + NCKI(ISYAIK)*(J - 1)
     *               + ICKI(ISYMAI,ISYMK) + 1
C
               NTOTAI = MAX(1,NT1AM(ISYMAI))
               NRHFL  = MAX(1,NRHF(ISYML))
C
C Occupied-L 1.0 contribution addWMAT
C
               CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                    NRHF(ISYML),ONE,T2TP(KOFF1),NTOTAI,
     *                    TROCCL(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,NCKIJ(JSAIKJ)
         WMAT(I) = WMAT(I) - TMAT(I)
      ENDDO
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: occ-L 1. Norm of WMAT ',XWMAT
      ENDIF
C
C----------------------------------------
C     Second occupied L contribution.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj)  
C                    + T2TP(ailC)*L(kl^jB)
C                                      
C                      TC(ail)*LB(l^jk) = R(aijk)
C
C----------------------------------------
C

      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
      ISYAIL = MULD2H(ISYMC,ISYMT2)
      ISYLKJ = MULD2H(ISYMB,ISINT2)
C
      CALL DZERO(TMAT,LENGTH)
C
      DO ISYMJ = 1,NSYM
C
         ISYMLK = MULD2H(ISYMJ,ISYLKJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMK = 1,NSYM
C
               ISYML  = MULD2H(ISYMK,ISYMLK)
               ISYMAI = MULD2H(ISYAIL,ISYML)
               ISYAIK = MULD2H(ISYMAI,ISYMK)
C
               KOFF1 = IT2SP(ISYAIL,ISYMC)
     *               + NCKI(ISYAIL)*(C - 1)
     *               + ICKI(ISYMAI,ISYML) + 1
               KOFF2 = ISJIKA(ISYLKJ,ISYMB)
     *               + NMAJIK(ISYLKJ)*(B - 1)
     *               + ISJIK(ISYMLK,ISYMJ)
     *               + NMATIJ(ISYMLK)*(J - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
               KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *               + NCKI(ISYAIK)*(J - 1)
     *               + ICKI(ISYMAI,ISYMK) + 1
C
               NTOTAI = MAX(1,NT1AM(ISYMAI))
               NRHFL  = MAX(1,NRHF(ISYML))
C
C Occupied-L 2.0 contribution addWMAT
C
               CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                    NRHF(ISYML),ONE,T2TP(KOFF1),NTOTAI,
     *                    TROCCL(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,NCKIJ(JSAIKJ)
         WMAT(I) = WMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: occ-L 2. Norm of WMAT ',XWMAT
      ENDIF

C
C
C----------------------------------------
C     First occupied g contribution.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj)  
C                    - T2TP(alkB)*g(il^jC)
C                                      
C                      TB(akl)*gC(l^ji) = R(akji)
C
C----------------------------------------
C
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
      ISYALK = MULD2H(ISYMB,ISYMT2)
      ISYLJI = MULD2H(ISYMC,ISINT2)
C
      KALK = 1
      KEND1  = KALK   + NCKI(ISYALK)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Not enough space in WBARBD_TMAT')
      END IF
C
      CALL DZERO(TMAT,NCKIJ(JSAIKJ))
C
C
C     T2TP(alkB) put in WORK(akl)
C
      KOFF = IT2SP(ISYALK,ISYMB) + NCKI(ISYALK)*(B - 1) + 1


      CALL CC_GATHER(NCKI(ISYALK),WORK(KALK),T2TP(KOFF),INDAJLB)
C
      DO ISYMI = 1,NSYM
C
         ISYAKJ = MULD2H(JSAIKJ,ISYMI)
         ISYMLJ = MULD2H(ISYLJI,ISYMI)
         DO I = 1,NRHF(ISYMI)
C
            DO ISYML = 1,NSYM
C
               ISYMAK = MULD2H(ISYALK,ISYML)
               ISYMJ  = MULD2H(ISYMLJ,ISYML)
C
               KOFF1 = KALK
     *               + ICKI(ISYMAK,ISYML) 
               KOFF2 = ISJIKA(ISYLJI,ISYMC)
     *               + NMAJIK(ISYLJI)*(C - 1)
     *               + ISJIK(ISYMLJ,ISYMI)
     *               + NMATIJ(ISYMLJ)*(I - 1)
     *               + IMATIJ(ISYML,ISYMJ) + 1
               KOFF3 = ISAIKJ(ISYAKJ,ISYMI)
     *               + NCKI(ISYAKJ)*(I - 1)
     *               + ICKI(ISYMAK,ISYMJ) + 1
C
               NTOTAK = MAX(1,NT1AM(ISYMAK))
               NRHFL  = MAX(1,NRHF(ISYML))
C
C Occupied-g 1.0 contribution addWMAT
C
               CALL DGEMM('N','N',NT1AM(ISYMAK),NRHF(ISYMJ),
     *                    NRHF(ISYML),ONE,WORK(KOFF1),NTOTAK,
     *                    TROCCG(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                    NTOTAK)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,NCKIJ(JSAIKJ)
         WMAT(I) = WMAT(I) + TMAT(INDSQ(I,2))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: occ-g 1. Norm of WMAT ',XWMAT
      ENDIF
C

C----------------------------------------
C     Second occupied g contribution.
C
C      WMAT^BC(aikj) = WMAT^BC(aikj)  
C                    - T2TP(aljC)*g(il^kB)
C                                      
C                      TC(ajl)*gB(l^ki) = R(ajki)
C----------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
      ISYAJL = MULD2H(ISYMC,ISYMT2)
      ISYLKI = MULD2H(ISYMB,ISINT2)
C
      KALJ = 1
      KEND1  = KALJ   + NCKI(ISYAJL)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Not enough space in WBARBD_TMAT')
      END IF
C
      CALL DZERO(TMAT,NCKIJ(JSAIKJ))
C
C
C     T2TP(aljC) put in WORK(ajl)
C
      KOFF = IT2SP(ISYAJL,ISYMC) + NCKI(ISYAJL)*(C - 1) + 1
      CALL CC_GATHER(NCKI(ISYAJL),WORK(KALJ),T2TP(KOFF),INDAJLC)
C
      DO ISYMI = 1,NSYM
C
         ISYAJK = MULD2H(JSAIKJ,ISYMI)
         ISYMLK = MULD2H(ISYLKI,ISYMI)
         DO I = 1,NRHF(ISYMI)
C
            DO ISYML = 1,NSYM
C
               ISYMAJ = MULD2H(ISYAJL,ISYML)
               ISYMK  = MULD2H(ISYMLK,ISYML)
C
               KOFF1 = KALJ
     *               + ICKI(ISYMAJ,ISYML) 
               KOFF2 = ISJIKA(ISYLKI,ISYMB)
     *               + NMAJIK(ISYLKI)*(B - 1)
     *               + ISJIK(ISYMLK,ISYMI)
     *               + NMATIJ(ISYMLK)*(I - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
               KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *               + NCKI(ISYAJK)*(I - 1)
     *               + ICKI(ISYMAJ,ISYMK) + 1
C
               NTOTAJ = MAX(1,NT1AM(ISYMAJ))
               NRHFL  = MAX(1,NRHF(ISYML))
C
C Occupied-g 2.0 contribution addWMAT
C
               CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMK),
     *                    NRHF(ISYML),ONE,WORK(KOFF1),NTOTAJ,
     *                    TROCCG(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                    NTOTAJ)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,NCKIJ(JSAIKJ)
         WMAT(I) = WMAT(I) + TMAT(INDSQ(I,5))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBARBD_TMAT: occ-g 2. Norm of WMAT ',XWMAT
      ENDIF
C
      CALL QEXIT('WBARBD_TMAT')
C
      RETURN
      END
C  /* Deck cc3_srtocc */
      SUBROUTINE CC3_SRTOCC(XINT,XLINT,ISYINT)
C
C----------------------------------------------------------------------
C integrals g(ia|j k-) sorted as G(j,i,k-,a) resorted as G(k-,i,j,a)
C                   input XINT  G(j,i,k-,a)  output  XINT  G(k-,i,j,a) 
C           L(ia|j k-)  output                       XLINT L(k-,i,j,a) 
C
C----------------------------------------------------------------------
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMA, ISYINT, ISYMK, ISYMJIK, ISYMJI, ISYMJ, ISYMKIJ
      INTEGER ISYMKI, ISYMKJI, ISYMKJ, ISYMJIKA, ISYMI
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL XINT(*), XLINT(*)
      REAL D2
#else
      DOUBLE PRECISION XINT(*), XLINT(*)
      DOUBLE PRECISION D2 
#endif
C
      PARAMETER ( D2 = 2.0D0)
C
C
      CALL QENTER('CC3_SRTOCC')
C
C----------------------------------------------------------------------
C integrals g(ia|j k-) sorted as XINT(j,i,k-,a) resorted  XINT(k-,i,j,a)
C----------------------------------------------------------------------
C
      DO ISYMA  = 1,NSYM
C
         ISYMJIK = MULD2H(ISYINT,ISYMA)
         ISYMKIJ = MULD2H(ISYINT,ISYMA)
         DO ISYMK = 1,NSYM
C
            ISYMJI = MULD2H(ISYMJIK,ISYMK)
            DO ISYMI = 1,NSYM
C
               ISYMJ = MULD2H(ISYMJI,ISYMI)
               ISYMKI = MULD2H(ISYMKIJ,ISYMJ)
C
               DO A = 1,NVIR(ISYMA)
C
                  DO K = 1,NRHF(ISYMK)
C
                     DO J = 1,NRHF(ISYMJ)
C
                        DO I = 1,NRHF(ISYMI)
C
                           KOFF1 = ISJIKA(ISYMJIK,ISYMA)
     *                        + NMAJIK(ISYMJIK)*(A-1)
     *                        + ISJIK(ISYMJI,ISYMK)
     *                        + NMATIJ(ISYMJI)*(K - 1)
     *                        + IMATIJ(ISYMJ,ISYMI)
     *                        + NRHF(ISYMJ)*(I - 1) + J
C
                           KOFF2 = ISJIKA(ISYMKIJ,ISYMA)
     *                        + NMAJIK(ISYMKIJ)*(A-1)
     *                        + ISJIK(ISYMKI,ISYMJ)
     *                        + NMATIJ(ISYMKI)*(J - 1)
     *                        + IMATIJ(ISYMK,ISYMI)
     *                        + NRHF(ISYMK)*(I - 1) + K

                           XLINT(KOFF2) = XINT(KOFF1)
C
                        END DO! J
                     ENDDO    ! I
                  ENDDO       ! K
               ENDDO          ! A
            ENDDO             ! ISYMI
         ENDDO                ! ISYMK
      ENDDO                   ! ISYMA
      
C
      CALL DCOPY(NCKIJ(ISYINT),XLINT,1,XINT,1)
C
C----------------------------------------------------------------------
C           L(ia|j k-) = 2.0 * XINT(k-,i,j,a) -  XINT(k-,i,j,a) 
C           L(ia|j k-) given as  XLINT(k-,i,j,a)
C----------------------------------------------------------------------
C
      DO ISYMA  = 1,NSYM
C
         ISYMKIJ = MULD2H(ISYINT,ISYMA)
         DO ISYMJ = 1,NSYM
C
            ISYMKI = MULD2H(ISYMKIJ,ISYMJ)
            DO ISYMI = 1,NSYM
C
               ISYMK   = MULD2H(ISYMKI,ISYMI)
               ISYMKJ  = MULD2H(ISYMK,ISYMJ)
               ISYMKJI = MULD2H(ISYMKI,ISYMJ)
C
               DO A = 1,NVIR(ISYMA)
C
                  DO J = 1,NRHF(ISYMJ)
C
                     DO I = 1,NRHF(ISYMI)
C
                        DO K = 1,NRHF(ISYMK)
                           KOFF1 = ISJIKA(ISYMKIJ,ISYMA)
     *                        + NMAJIK(ISYMKIJ)*(A-1)
     *                        + ISJIK(ISYMKI,ISYMJ)
     *                        + NMATIJ(ISYMKI)*(J - 1)
     *                        + IMATIJ(ISYMK,ISYMI)
     *                        + NRHF(ISYMK)*(I - 1) + K
C
                           KOFF2 = ISJIKA(ISYMKJI,ISYMA)
     *                        + NMAJIK(ISYMKJI)*(A-1)
     *                        + ISJIK(ISYMKJ,ISYMI)
     *                        + NMATIJ(ISYMKJ)*(I - 1)
     *                        + IMATIJ(ISYMK,ISYMJ) 
     *                        + NRHF(ISYMK)*(J - 1) + K

                           XLINT(KOFF1) = D2 * XINT(KOFF1) - XINT(KOFF2)
C
                        END DO! K
                     ENDDO    ! I
                  ENDDO       ! J
               ENDDO          ! A
            ENDDO             ! ISYMI
         ENDDO                ! ISYMJ
      ENDDO                   ! ISYMA

      CALL QEXIT('CC3_SRTOCC')
C
      RETURN
      END
