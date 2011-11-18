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
C  /* Deck cc_yi */
      SUBROUTINE CC_YI(YMAT,CLTR2,ISYMLV,T2AM,ISYMRV,WORK,LWORK)
C
C     Written by Asger Halkier 16/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the Y-intermediat entering the E'-term
C              of both the 2.1-block and 2.2-block of Jacobian,
C              and in the F-matrix and etaC.
C
C     It is assumed that CLTR2 is squared up, and that T2AM is not.
C
C     Ove Christiansen 20-6-1996:
C
C     General symmetry of both CLTR2 and T2AM for F-mat. and etaC.
C
#include <implicit.h>
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION YMAT(*),CLTR2(*),T2AM(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYRES = MULD2H(ISYMLV,ISYMRV)
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      CALL DZERO(YMAT,NMATAB(ISYRES))
C
      DO 100 ISYMCM = 1,NSYM
C
         ISYMFN = MULD2H(ISYMCM,ISYMLV)
         ISYMEN = MULD2H(ISYMCM,ISYMRV)
C
         DO 110 ISYMN = 1,NSYM
C
            ISYME = MULD2H(ISYMN,ISYMEN)
            ISYMF = MULD2H(ISYMN,ISYMFN)
C
C-------------------------------------
C           Work space allocation one.
C-------------------------------------
C
            KT2SM = 1
            KC2SM = KT2SM + NT1AM(ISYMCM)*NVIR(ISYME)
            KEND1 = KC2SM + NT1AM(ISYMCM)*NVIR(ISYMF)
            LWRK1 = LWORK - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_YI')
            ENDIF
C
            CALL DZERO(WORK,KEND1)
C
            DO 120 N = 1,NRHF(ISYMN)
C
C------------------------------------------------
C              Copy submatrix out of packed T2AM.
C------------------------------------------------
C
               DO 130 E = 1,NVIR(ISYME)
C
                  NEN = IT1AM(ISYME,ISYMN) + NVIR(ISYME)*(N - 1) + E
C
                  IF (ISYMCM .EQ. ISYMEN) THEN
C
                     DO 140 NCM = 1,NT1AM(ISYMCM)
C
                        NCMEN = IT2AM(ISYMCM,ISYMEN) + INDEX(NCM,NEN)
                        NECM  = KT2SM + NVIR(ISYME)*(NCM - 1) + E - 1
C
                        WORK(NECM) = T2AM(NCMEN)
C
  140                CONTINUE
C
                  ELSE IF (ISYMCM .LT. ISYMEN) THEN
C
                     DO 150 NCM = 1,NT1AM(ISYMCM)
C
                        NCMEN = IT2AM(ISYMCM,ISYMEN)
     &                        + NT1AM(ISYMCM)*(NEN - 1) + NCM
                        NECM  = KT2SM + NVIR(ISYME)*(NCM - 1) + E - 1
C
                        WORK(NECM) = T2AM(NCMEN)
C
  150                CONTINUE
C
                  ELSE IF (ISYMCM .GT. ISYMEN) THEN
C
                     DO 160 NCM = 1,NT1AM(ISYMCM)
C
                        NCMEN = IT2AM(ISYMEN,ISYMCM)
     &                        + NT1AM(ISYMEN)*(NCM - 1) + NEN
                        NECM  = KT2SM + NVIR(ISYME)*(NCM - 1) + E - 1
C
                        WORK(NECM) = T2AM(NCMEN)
C
  160                CONTINUE
C
                  ENDIF
C
  130          CONTINUE
C
C-----------------------------------------
C              Copy submatrix out of CLTR2.
C-----------------------------------------
C
               DO 170 F = 1,NVIR(ISYMF)
C
                  NFN = IT1AM(ISYMF,ISYMN) + NVIR(ISYMF)*(N - 1) + F
C
                  DO 180 NCM = 1,NT1AM(ISYMCM)
C
                     NCMFN = IT2SQ(ISYMCM,ISYMFN)
     &                     + NT1AM(ISYMCM)*(NFN - 1) + NCM
                     NCMF  = KC2SM + NT1AM(ISYMCM)*(F - 1) + NCM - 1
C
                     WORK(NCMF) = CLTR2(NCMFN)
C
  180             CONTINUE
  170          CONTINUE
C
C-------------------------------------------
C              Contract the two submatrices.
C-------------------------------------------
C
               KOFF1  = IMATAB(ISYME,ISYMF) + 1
C
               NTOTE  = MAX(NVIR(ISYME),1)
               NTOTCM = MAX(NT1AM(ISYMCM),1)
C
               CALL DGEMM('N','N',NVIR(ISYME),NVIR(ISYMF),NT1AM(ISYMCM),
     &                    ONE,WORK(KT2SM),NTOTE,WORK(KC2SM),NTOTCM,
     &                    ONE,YMAT(KOFF1),NTOTE)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) '>>>>>CC_YI: Norm of Y-intermediate:',
     &    DDOT(NMATAB(ISYRES),YMAT,1,YMAT,1)
      END IF
C
      RETURN
      END
C  /* Deck cc_xi */
      SUBROUTINE CC_XI(XMAT,CLTR2,ISYMLV,T2AM,ISYMRV,WORK,LWORK)
C
C     Written by Asger Halkier 16/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the X-intermediat entering the E'-term
C              of both the 2.1-block and 2.2-block of Jacobian,
C              and in the F-matrix and etaC.
C
C     Ove Christiansen 20-6-1996:
C
C     General symmetry of both CLTR2 and T2AM for F-mat. and etaC.
C
C
#include <implicit.h>
      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XMAT(*),CLTR2(*),T2AM(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      ISYRES = MULD2H(ISYMLV,ISYMRV)
C
      CALL DZERO(XMAT,NMATIJ(ISYRES))
C
      DO 100 ISYMEM = 1,NSYM
C
         ISYMFL = MULD2H(ISYMEM,ISYMRV)
         ISYMFJ = MULD2H(ISYMEM,ISYMLV)
C
         DO 110 ISYMF = 1,NSYM
C
            ISYML = MULD2H(ISYMF,ISYMFL)
            ISYMJ = MULD2H(ISYMF,ISYMFJ)
C
C-------------------------------------
C           Work space allocation one.
C-------------------------------------
C
            KT2SM = 1
            KC2SM = KT2SM + NT1AM(ISYMEM)*NRHF(ISYML)
            KEND1 = KC2SM + NT1AM(ISYMEM)*NRHF(ISYMJ)
            LWRK1 = LWORK - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_XI')
            ENDIF
C
            CALL DZERO(WORK,KEND1)
C
            DO 120 F = 1,NVIR(ISYMF)
C
C--------------------------------------
C           Copy submatrix of T2AM out.
C--------------------------------------
C
               DO 130 L = 1,NRHF(ISYML)
C
                  NFL = IT1AM(ISYMF,ISYML) + NVIR(ISYMF)*(L - 1) + F
C
                  IF (ISYMEM .EQ. ISYMFL) THEN
C
                     DO 140 NEM = 1,NT1AM(ISYMEM)
C
                        NEMFL = IT2AM(ISYMEM,ISYMFL) + INDEX(NEM,NFL)
                        NLEM  = KT2SM + NRHF(ISYML)*(NEM - 1) + L - 1
C
                        WORK(NLEM) = T2AM(NEMFL)
C
  140                CONTINUE
C
                  ELSE IF (ISYMEM .LT. ISYMFL) THEN
C
                     DO 150 NEM = 1,NT1AM(ISYMEM)
C
                        NEMFL = IT2AM(ISYMEM,ISYMFL)
     &                        + NT1AM(ISYMEM)*(NFL - 1) + NEM
                        NLEM  = KT2SM + NRHF(ISYML)*(NEM - 1) + L - 1
C
                        WORK(NLEM) = T2AM(NEMFL)
C
  150                CONTINUE
C
                  ELSE IF (ISYMEM .GT. ISYMFL) THEN
C
                     DO 160 NEM = 1,NT1AM(ISYMEM)
C
                        NEMFL = IT2AM(ISYMFL,ISYMEM)
     &                        + NT1AM(ISYMFL)*(NEM - 1) + NFL
                        NLEM  = KT2SM + NRHF(ISYML)*(NEM - 1) + L - 1
C
                        WORK(NLEM) = T2AM(NEMFL)
C
  160                CONTINUE
C
                  ENDIF
C
  130          CONTINUE
C
C-----------------------------------------
C              Copy submatrix of CLTR2 out.
C-----------------------------------------
C
               DO 170 J = 1,NRHF(ISYMJ)
C
                  NFJ = IT1AM(ISYMF,ISYMJ) + NVIR(ISYMF)*(J - 1) + F
C
                  DO 180 NEM = 1,NT1AM(ISYMEM)
C
                     NEMFJ = IT2SQ(ISYMEM,ISYMFJ)
     &                     + NT1AM(ISYMEM)*(NFJ - 1) + NEM
                     NEMJ  = KC2SM + NT1AM(ISYMEM)*(J - 1) + NEM - 1
C
                     WORK(NEMJ) = CLTR2(NEMFJ)
C
  180             CONTINUE
  170          CONTINUE
C
C-------------------------------------------
C              Contract the two submatrices.
C-------------------------------------------
C
               KOFF1  = IMATIJ(ISYML,ISYMJ) + 1
C
               NTOTL  = MAX(NRHF(ISYML),1)
               NTOTEM = MAX(NT1AM(ISYMEM),1)
C
               CALL DGEMM('N','N',NRHF(ISYML),NRHF(ISYMJ),NT1AM(ISYMEM),
     &                    ONE,WORK(KT2SM),NTOTL,WORK(KC2SM),NTOTEM,ONE,
     &                    XMAT(KOFF1),NTOTL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF (LOCDBG) THEN 
        WRITE(LUPRI,*) '>>>>>CC_XI: Norm of X-intermediate:', 
     &    DDOT(NMATIJ(ISYRES),XMAT,1,XMAT,1)
      END IF
C
      RETURN
      END
C  /* Deck cc_21efm */
      SUBROUTINE CC_21EFM(RHO1,FCKMO,ISYFCK,XMAT,YMAT,ISYMIM,WORK,LWORK)
C
C     Written by Asger Halkier 28/7 - 1995.
C
C     Version: 1.0
C
C     Purpose: To calculate the Fock contributions to the E' term
C              of the 1.2-block.
C
C     Ove Christiansen 20-6-1996 
C              General symmetry of fckmo(ISYFCK) and X and Y 
C              intermediates ISYMIM
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION RHO1(*),FCKMO(*),XMAT(*),YMAT(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      IF (LWORK .LT. NT1AM(ISYFCK)) THEN
         CALL QUIT('Insufficient work-space-area in CC_12EMF')
      ENDIF
C
      ISYRES = MULD2H(ISYFCK,ISYMIM)
C
C-----------------------------------
C     Extract the fock matrix F(ie).
C-----------------------------------
C
      DO 50 ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISYFCK)
C
         DO 60 K = 1,NRHF(ISYMK)
C
            DO 70 C = 1,NVIR(ISYMC)
C
               KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
               KOFF2 = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               WORK(KOFF2) = FCKMO(KOFF1)
C
   70       CONTINUE
   60    CONTINUE
   50 CONTINUE
C
C-----------------------------------------------
C     Calculate contribution from the XMAT-part.
C-----------------------------------------------
C
      DO 80 ISYMA = 1,NSYM
C
         ISYML = MULD2H(ISYMA,ISYFCK)
         ISYMI = MULD2H(ISYMA,ISYRES)
C
         KOFF1 = IT1AM(ISYMA,ISYML) + 1
         KOFF2 = IMATIJ(ISYML,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTL = MAX(NRHF(ISYML),1)
         NTOTA = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYML),-ONE,
     &              WORK(KOFF1),NTOTA,XMAT(KOFF2),NTOTL,ONE,
     &              RHO1(KOFF3),NTOTA)
C
  80  CONTINUE
C
C-----------------------------------------------
C     Calculate contribution from the YMAT-part.
C-----------------------------------------------
C
      DO 90 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMA,ISYRES)
         ISYME = MULD2H(ISYMI,ISYFCK)
C
         KOFF1 = IMATAB(ISYME,ISYMA) + 1
         KOFF2 = IT1AM(ISYME,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTE = MAX(NVIR(ISYME),1)
         NTOTA = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYME),-ONE,
     &              YMAT(KOFF1),NTOTE,WORK(KOFF2),NTOTE,ONE,
     &              RHO1(KOFF3),NTOTA)
C
  90  CONTINUE
C
      RETURN
      END
C  /* Deck cc_22ec */
      SUBROUTINE CC_22EC(RHO2,T2AM,XMAT,YMAT,ISYMXY,WORK,LWORK)
C
C     Written by Asger Halkier 21/11 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the coulomb part of the 22E-prime
C              contribution to RHO2.
C
C     Ove Christiansen 17-9-1996:
C              Gen. sym. input of XMAT and YMAT: ISYMXY
C              (ia|jb) quantity is total symmetric and packed
C              as ai,bj.
C
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO2(*), T2AM(*), XMAT(*), YMAT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYRES = MULD2H(ISYMOP,ISYMXY)
      ISYMAT = ISYRES
C
      DO 100 ISYMJ = 1,NSYM
C
         IF (NRHF(ISYMJ) .EQ. 0) GOTO 100
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYRES)
               ISYMAN = MULD2H(ISYMBJ,ISYMOP)
               ISYMFI = ISYMAN
C
               IF (NVIR(ISYMB) .EQ. 0) GOTO 120
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KSCR1 = 1
               KSCR2 = KSCR1 + NT1AM(ISYMAI)
               KEND1 = KSCR2 + NT1AM(ISYMAN)
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_22EC')
               ENDIF
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KSCR1),NT1AM(ISYMAI))
C
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
C-------------------------------------------------
C                 Copy submatrix out of integrals.
C-------------------------------------------------
C
                  DO 140 NAN = 1,NT1AM(ISYMAN)
C
                     NANBJ = IT2AM(ISYMAN,ISYMBJ) + INDEX(NAN,NBJ)
C
                     WORK(KSCR2 + NAN - 1) = T2AM(NANBJ)
C
  140             CONTINUE
C
C-----------------------------------------------------------------
C                 Contraction of integrals with X- and Y-matrices.
C-----------------------------------------------------------------
C
                  DO 150 ISYMAX = 1,NSYM
C
                     ISYMNX = MULD2H(ISYMAX,ISYMAN)
                     ISYMIX = MULD2H(ISYMNX,ISYMAT)
                     ISYMFY = ISYMAX
                     ISYMIY = MULD2H(ISYMFY,ISYMFI)
                     ISYMAY = MULD2H(ISYMFY,ISYMAT)
C
                     KOFF1 = KSCR2 + IT1AM(ISYMAX,ISYMNX)
                     KOFF2 = IMATIJ(ISYMNX,ISYMIX) + 1
                     KOFF3 = KSCR1 + IT1AM(ISYMAX,ISYMIX)
C
                     NTOTA = MAX(NVIR(ISYMAX),1)
                     NTOTN = MAX(NRHF(ISYMNX),1)
C
                     CALL DGEMM('N','N',NVIR(ISYMAX),NRHF(ISYMIX),
     &                          NRHF(ISYMNX),ONE,WORK(KOFF1),NTOTA,
     &                          XMAT(KOFF2),NTOTN,ONE,WORK(KOFF3),
     &                          NTOTA)
C
                     KOFF4 = IMATAB(ISYMFY,ISYMAY) + 1
                     KOFF5 = KSCR2 + IT1AM(ISYMFY,ISYMIY)
                     KOFF6 = KSCR1 + IT1AM(ISYMAY,ISYMIY)
C
                     NTOTF = MAX(NVIR(ISYMFY),1)
                     NTOTA = MAX(NVIR(ISYMAY),1)
C
                     CALL DGEMM('T','N',NVIR(ISYMAY),NRHF(ISYMIY),
     &                          NVIR(ISYMFY),ONE,YMAT(KOFF4),NTOTF,
     &                          WORK(KOFF5),NTOTF,ONE,WORK(KOFF6),
     &                          NTOTA)
C
  150             CONTINUE
C
C-----------------------------------------------
C                 Storage in RHO2 result vector.
C-----------------------------------------------
C
                  IF (ISYMAI .EQ. ISYMBJ) THEN
C
                     WORK(KSCR1 + NBJ - 1) = TWO*WORK(KSCR1 + NBJ - 1)
C
                     DO 160 NAI = 1,NT1AM(ISYMAI)
C
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                        + INDEX(NAI,NBJ)
C
                        RHO2(NAIBJ) = RHO2(NAIBJ)
     &                              - TWO*WORK(KSCR1 + NAI - 1)
C
  160                CONTINUE
C
                  ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
                     KOFF7 = IT2AM(ISYMAI,ISYMBJ)
     &                     + NT1AM(ISYMAI)*(NBJ - 1) + 1
C
                     CALL DAXPY(NT1AM(ISYMAI),-TWO,WORK(KSCR1),1,
     &                          RHO2(KOFF7),1)
C
                  ELSE IF (ISYMAI .GT. ISYMBJ) THEN
C
                     KOFF8 = IT2AM(ISYMBJ,ISYMAI) + NBJ
C
                     CALL DAXPY(NT1AM(ISYMAI),-TWO,WORK(KSCR1),1,
     &                          RHO2(KOFF8),NT1AM(ISYMBJ))
C
                  ENDIF
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_22ee */
      SUBROUTINE CC_22EE(RHO2,T2AM,XMAT,YMAT,ISYMXY,WORK,LWORK)
C
C     Written by Asger Halkier 21/11 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the exchange part of the 22E-prime
C              contribution to RHO2.
C
C     Ove Christiansen 17-9-1996:
C              Gen. sym. input of XMAT and YMAT: ISYMXY
C              (ia|jb) quantity is total symmetric and packed
C              as ai,bj.
C
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO2(*), T2AM(*), XMAT(*), YMAT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYRES = MULD2H(ISYMOP,ISYMXY)
      ISYMAT = ISYRES
C
      DO 100 ISYMI = 1,NSYM
C
         IF (NRHF(ISYMI) .EQ. 0) GOTO 100
C
         DO 110 I = 1,NRHF(ISYMI)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMAJ = MULD2H(ISYMBI,ISYRES)
               ISYMAN = MULD2H(ISYMBI,ISYMOP)
               ISYMFJ = ISYMAN
C
               IF (NVIR(ISYMB) .EQ. 0) GOTO 120
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KSCR1 = 1
               KSCR2 = KSCR1 + NT1AM(ISYMAJ)
               KEND1 = KSCR2 + NT1AM(ISYMAN)
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_22EE')
               ENDIF
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KSCR1),NT1AM(ISYMAJ))
C
                  NBI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1) + B
C
C-------------------------------------------------
C                 Copy submatrix out of integrals.
C-------------------------------------------------
C
                  DO 140 NAN = 1,NT1AM(ISYMAN)
C
                     NANBI = IT2AM(ISYMAN,ISYMBI) + INDEX(NAN,NBI)
C
                     WORK(KSCR2 + NAN - 1) = T2AM(NANBI)
C
  140             CONTINUE
C
C-----------------------------------------------------------------
C                 Contraction of integrals with X- and Y-matrices.
C-----------------------------------------------------------------
C
                  DO 150 ISYMAX = 1,NSYM
C
                     ISYMNX = MULD2H(ISYMAX,ISYMAN)
                     ISYMJX = MULD2H(ISYMNX,ISYMAT)
                     ISYMFY = ISYMAX
                     ISYMJY = MULD2H(ISYMFY,ISYMFJ)
                     ISYMAY = MULD2H(ISYMFY,ISYMAT)
C
                     KOFF1 = KSCR2 + IT1AM(ISYMAX,ISYMNX)
                     KOFF2 = IMATIJ(ISYMNX,ISYMJX) + 1
                     KOFF3 = KSCR1 + IT1AM(ISYMAX,ISYMJX)
C
                     NTOTA = MAX(NVIR(ISYMAX),1)
                     NTOTN = MAX(NRHF(ISYMNX),1)
C
                     CALL DGEMM('N','N',NVIR(ISYMAX),NRHF(ISYMJX),
     &                          NRHF(ISYMNX),ONE,WORK(KOFF1),NTOTA,
     &                          XMAT(KOFF2),NTOTN,ONE,WORK(KOFF3),
     &                          NTOTA)
C
                     KOFF4 = IMATAB(ISYMFY,ISYMAY) + 1
                     KOFF5 = KSCR2 + IT1AM(ISYMFY,ISYMJY)
                     KOFF6 = KSCR1 + IT1AM(ISYMAY,ISYMJY)
C
                     NTOTF = MAX(NVIR(ISYMFY),1)
                     NTOTA = MAX(NVIR(ISYMAY),1)
C
                     CALL DGEMM('T','N',NVIR(ISYMAY),NRHF(ISYMJY),
     &                          NVIR(ISYMFY),ONE,YMAT(KOFF4),NTOTF,
     &                          WORK(KOFF5),NTOTF,ONE,WORK(KOFF6),
     &                          NTOTA)
C
  150             CONTINUE
C
C-----------------------------------------------
C                 Storage in RHO2 result vector.
C-----------------------------------------------
C
                  DO 160 ISYMJ = 1,NSYM
C
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
                     ISYMAI = MULD2H(ISYMBJ,ISYRES)
                     ISYMA  = MULD2H(ISYMJ,ISYMAJ)
C
                     DO 170 J = 1,NRHF(ISYMJ)
C
                        NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1)
     &                      + B
C
                        IF (ISYMAI .EQ. ISYMBJ) THEN
C
                           DO 180 A = 1,NVIR(ISYMA)
C
                              NAJ   = IT1AM(ISYMA,ISYMJ)
     &                              + NVIR(ISYMA)*(J - 1) + A
                              NAI   = IT1AM(ISYMA,ISYMI)
     &                              + NVIR(ISYMA)*(I - 1) + A
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                              + INDEX(NAI,NBJ)
C
                              RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                    + WORK(KSCR1 + NAJ - 1)
C
                              IF (NAI .EQ. NBJ) THEN
C
                                 RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                       + WORK(KSCR1 + NAJ - 1)
C
                              ENDIF
C
  180                      CONTINUE
C
                        ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
                           KOFF7 = KSCR1 + IT1AM(ISYMA,ISYMJ)
     &                           + NVIR(ISYMA)*(J - 1)
                           KOFF8 = IT2AM(ISYMAI,ISYMBJ)
     &                           + NT1AM(ISYMAI)*(NBJ - 1)
     &                           + IT1AM(ISYMA,ISYMI)
     &                           + NVIR(ISYMA)*(I - 1) + 1
C
                           CALL DAXPY(NVIR(ISYMA),ONE,WORK(KOFF7),1,
     &                                RHO2(KOFF8),1)
C
                        ELSE IF (ISYMAI .GT. ISYMBJ) THEN
C
                           DO 190 A = 1,NVIR(ISYMA)
C
                              NAJ   = IT1AM(ISYMA,ISYMJ)
     &                              + NVIR(ISYMA)*(J - 1) + A
                              NAI   = IT1AM(ISYMA,ISYMI)
     &                              + NVIR(ISYMA)*(I - 1) + A
                              NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     &                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
C
                              RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                    + WORK(KSCR1 + NAJ - 1)
C
  190                      CONTINUE
C
                        ENDIF
C
  170                CONTINUE
  160             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_22c */
      SUBROUTINE CC_22C(RHO2,CTR2,ISYVEC,XLAMDH,WORK,LWORK,
     *                  ISYCIM,LUC,CFIL,IV,IOPT)
C
C     Written by Asger Halkier 27/11 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the 22C contribution to RHO2.
C
C     Ove Christiansen 17-9-1996: 
C
C              modified to account for general
C              non. total symmetric vectors (ISYVEC) and
C              intermediates (ISYCIM). LUC and CFIL is
C              used to control from which file the
C              intermediate is obtained.
C
C              if iopt = 1 the C intermediate is assumed
C              to be as in energy calc.
C
C              if iopt ne. 1 we use the intermediate
C              on lud with address given according to
C              transformed vector nr iv (iv is not 1
C              if several vectors are transformed
C              simultaneously.)
C
C              in ordinary Lambda equations: iv=1,iopt=1
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO2(*), CTR2(*), XLAMDH(*), WORK(LWORK)
      CHARACTER CFIL*(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <maxorb.h>
#include <ccsdio.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYRES = MULD2H(ISYVEC,ISYCIM)
C
      DO 100 ISYMD = 1,NSYM
C
         ISYMB  = ISYMD
         ISYIEM = MULD2H(ISYMD,ISYCIM)
         ISYAJI = MULD2H(ISYMB,ISYRES)
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         KPINT = 1
         KSCR1 = KPINT + NT2BCD(ISYIEM)
         KEND1 = KSCR1 + NT2BCD(ISYAJI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_22C')
         ENDIF
C
         DO 110 IDEL = 1,NBAS(ISYMD)
C
C---------------------------------------------------
C           Read P-intermediate submatrix from disc.
C---------------------------------------------------
C
            ID = IDEL + IBAS(ISYMD)
C
            IF ( IOPT .EQ. 1) THEN
               IOFF = IT2DEL(ID) + 1
            ELSE
               IOFF = IT2DLR(ID,IV) + 1
            ENDIF
C
            LEN = NT2BCD(ISYIEM)
C
            IF (LEN .GT. 0) THEN
               CALL GETWA2(LUC,CFIL,WORK(KPINT),IOFF,LEN)
            ENDIF
C
C--------------------------------------------------------
C           Contraction of intermediate and trial vector.
C--------------------------------------------------------
C
            DO 120 ISYMAJ = 1,NSYM
C
               ISYMEM = MULD2H(ISYMAJ,ISYVEC)
               ISYMI  = MULD2H(ISYMAJ,ISYAJI)
C
               KOFF1  = IT2SQ(ISYMAJ,ISYMEM) + 1
               KOFF2  = KPINT + IT2BCT(ISYMI,ISYMEM)
               KOFF3  = KSCR1 + IT2BCD(ISYMAJ,ISYMI)
C
               NTOTAJ = MAX(NT1AM(ISYMAJ),1)
               NTOTI  = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('N','T',NT1AM(ISYMAJ),NRHF(ISYMI),
     &                    NT1AM(ISYMEM),ONE,CTR2(KOFF1),NTOTAJ,
     &                    WORK(KOFF2),NTOTI,ZERO,WORK(KOFF3),NTOTAJ)
C
  120       CONTINUE
C
C---------------------------------------------
C           Resort intermediate result-matrix.
C---------------------------------------------
C
            CALL CCLT_P21I(WORK(KSCR1),ISYAJI,WORK(KEND1),LWRK1,
     &                     IT2BCD,NT2BCD,IT1AM,NT1AM,NVIR)
C
C-----------------------------------------------------------------
C           Final scaling with XLAMDH-element and storage in RHO2.
C-----------------------------------------------------------------
C
            DO 130 B = 1,NVIR(ISYMB)
C
               NDB = ILMVIR(ISYMB) + NBAS(ISYMD)*(B - 1) + IDEL
C
               FACT = -HALF*XLAMDH(NDB)
C
               DO 140 ISYMJ = 1,NSYM
C
                  ISYMAI = MULD2H(ISYMJ,ISYAJI)
                  ISYMBJ = MULD2H(ISYMAI,ISYRES)
C
                  DO 150 J = 1,NRHF(ISYMJ)
C
                     NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
                     NJ  = KSCR1 + IT2BCD(ISYMAI,ISYMJ)
     &                   + NT1AM(ISYMAI)*(J - 1)
C
                     IF (ISYMAI .EQ. ISYMBJ) THEN
C
                        DO 160 NAI = 1,NT1AM(ISYMAI)
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                           + INDEX(NAI,NBJ)
C
                           RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                 + FACT*WORK(NJ + NAI - 1)
C
                           IF (NAI .EQ. NBJ) THEN
C
                              RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                    + FACT*WORK(NJ + NAI - 1)
C
                           ENDIF
C
  160                   CONTINUE
C
                     ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                        + NT1AM(ISYMAI)*(NBJ - 1) + 1
C
                        CALL DAXPY(NT1AM(ISYMAI),FACT,WORK(NJ),1,
     &                             RHO2(NAIBJ),1)
C
                     ELSE IF (ISYMAI .GT. ISYMBJ) THEN
C
                        NAIBJ = IT2AM(ISYMBJ,ISYMAI) + NBJ
C
                        CALL DAXPY(NT1AM(ISYMAI),FACT,WORK(NJ),1,
     &                             RHO2(NAIBJ),NT1AM(ISYMBJ))
C
                     ENDIF
C
  150             CONTINUE
  140          CONTINUE
  130       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_22d */
      SUBROUTINE CC_22D(RHO2,CTR2,ISYVEC,XLAMDH,WORK,LWORK,
     *                  ISYDIM,LUD,DFIL,IV,IOPT)
C
C     Written by Asger Halkier 27/11 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the 22D contribution to RHO2!
C
C     Ove Christiansen 17-9-1996: 
C                 Modified to account for general
C                 non. total symmetric vectors (ISYVEC) and
C                 intermediates (ISYDIM). LUD and DFIL is
C                 used to control from which file the
C                 intermediate is obtained.
C
C                 if iopt = 1 the D intermediate is assumed
C                 to be as in energy calc.
C
C                 if iopt ne. 1 we use the intermediate
C                 on lud with address given according to
C                 transformed vector nr iv (iv is not 1
C                 if several vectors are transformed
C                 simultaneously.)
C
C                 in Lambda vector calc: iv=1,iopt=1
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO2(*), CTR2(*), XLAMDH(*), WORK(LWORK)
      CHARACTER DFIL*(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <maxorb.h>
#include <ccsdio.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYRES = MULD2H(ISYVEC,ISYDIM)
C
      DO 100 ISYMD = 1,NSYM
C
         ISYMB  = ISYMD
         ISYJEM = MULD2H(ISYMD,ISYDIM)
         ISYAIJ = MULD2H(ISYMB,ISYRES)
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         KPINT = 1
         KSCRN = KPINT + NT2BCD(ISYJEM)
         KSCRT = KSCRN + NT2BCD(ISYAIJ)
         KEND1 = KSCRT + NT2BCD(ISYAIJ)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_22D')
         ENDIF
C
         DO 110 IDEL = 1,NBAS(ISYMD)
C
C---------------------------------------------------
C           Read P-intermediate submatrix from disc.
C---------------------------------------------------
C
            ID = IDEL + IBAS(ISYMD)
C
            IF (IOPT .EQ. 1) THEN
               IOFF = IT2DEL(ID) + 1
            ELSE
               IOFF = IT2DLR(ID,IV) + 1
            ENDIF
C
            LEN = NT2BCD(ISYJEM)
C
            IF (LEN .GT. 0) THEN
               CALL GETWA2(LUD,DFIL,WORK(KPINT),IOFF,LEN)
C
            ENDIF
C
C--------------------------------------------------------
C           Contraction of intermediate and trial vector.
C--------------------------------------------------------
C
            DO 120 ISYMAI = 1,NSYM
C
               ISYMJ  = MULD2H(ISYMAI,ISYAIJ)
               ISYMEM = MULD2H(ISYMAI,ISYVEC)
C
               KOFF1  = IT2SQ(ISYMAI,ISYMEM) + 1
               KOFF2  = KPINT + IT2BCT(ISYMJ,ISYMEM)
               KOFF3  = KSCRN + IT2BCD(ISYMAI,ISYMJ)
C
               NTOTAI = MAX(NT1AM(ISYMAI),1)
               NTOTJ  = MAX(NRHF(ISYMJ),1)
C
               CALL DGEMM('N','T',NT1AM(ISYMAI),NRHF(ISYMJ),
     &                    NT1AM(ISYMEM),ONE,CTR2(KOFF1),NTOTAI,
     &                    WORK(KOFF2),NTOTJ,ZERO,WORK(KOFF3),NTOTAI)
C
  120       CONTINUE
C
C----------------------------------------
C           Construct resorted submatrix.
C----------------------------------------
C
            CALL DCOPY(NT2BCD(ISYAIJ),WORK(KSCRN),1,WORK(KSCRT),1)
C
            CALL CCLT_P21I(WORK(KSCRT),ISYAIJ,WORK(KEND1),LWRK1,
     &                     IT2BCD,NT2BCD,IT1AM,NT1AM,NVIR)
C
C-----------------------------------------------------------------
C           Final scaling with XLAMDH-element and storage in RHO2.
C-----------------------------------------------------------------
C
            DO 130 B = 1,NVIR(ISYMB)
C
               NDB = ILMVIR(ISYMB) + NBAS(ISYMD)*(B - 1) + IDEL
C
               FACTN = XLAMDH(NDB)
               FACTT = -HALF*XLAMDH(NDB)
C
               DO 140 ISYMJ = 1,NSYM
C
                  ISYMAI = MULD2H(ISYMJ,ISYAIJ)
                  ISYMBJ = MULD2H(ISYMAI,ISYRES)
C
                  DO 150 J = 1,NRHF(ISYMJ)
C
                     NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
                     NJN = KSCRN + IT2BCD(ISYMAI,ISYMJ)
     &                   + NT1AM(ISYMAI)*(J - 1)
                     NJT = KSCRT + IT2BCD(ISYMAI,ISYMJ)
     &                   + NT1AM(ISYMAI)*(J - 1)
C
                     IF (ISYMAI .EQ. ISYMBJ) THEN
C
                        DO 160 NAI = 1,NT1AM(ISYMAI)
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                           + INDEX(NAI,NBJ)
C
                           RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                 + FACTN*WORK(NJN + NAI - 1)
     &                                 + FACTT*WORK(NJT + NAI - 1)
C
                           IF (NAI .EQ. NBJ) THEN
C
                              RHO2(NAIBJ) = RHO2(NAIBJ)
     &                                    + FACTN*WORK(NJN + NAI - 1)
     &                                    + FACTT*WORK(NJT + NAI - 1)
C
                           ENDIF
C
  160                   CONTINUE
C
                     ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                        + NT1AM(ISYMAI)*(NBJ - 1) + 1
C
                        CALL DAXPY(NT1AM(ISYMAI),FACTN,WORK(NJN),1,
     &                             RHO2(NAIBJ),1)
                        CALL DAXPY(NT1AM(ISYMAI),FACTT,WORK(NJT),1,
     &                             RHO2(NAIBJ),1)
C
                     ELSE IF (ISYMAI .GT. ISYMBJ) THEN
C
                        NAIBJ = IT2AM(ISYMBJ,ISYMAI) + NBJ
C
                        CALL DAXPY(NT1AM(ISYMAI),FACTN,WORK(NJN),1,
     &                             RHO2(NAIBJ),NT1AM(ISYMBJ))
                        CALL DAXPY(NT1AM(ISYMAI),FACTT,WORK(NJT),1,
     &                             RHO2(NAIBJ),NT1AM(ISYMBJ))
C
                     ENDIF
C
  150             CONTINUE
  140          CONTINUE
  130       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_mi */
      SUBROUTINE CC_MI(XMINT,CTR2,IC2SYM,T2AM,IT2SYM,WORK,LWORK)
C
C     Written by Asger Halkier 28/10 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the M-intermediat entering the G-term
C              of both the 2.1-block.
C
C     It is assumed that CTR2 is squared up, and that T2AM is not!
C
C     Ove Christiansen 18-9-1996:
C              Non-total symmetric CTR2 and T2AM amplitudes.
C              Sym:               IC2SYM   IT2SYM
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XMINT(*),CTR2(*),T2AM(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYINT = MULD2H(IC2SYM,IT2SYM)
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      CALL DZERO(XMINT,N3ORHF(ISYINT))
C
      DO 100 ISYME = 1,NSYM
C
         DO 110 E = 1,NVIR(ISYME)
C
            DO 120 ISYML = 1,NSYM
C
               ISYMEL = MULD2H(ISYML,ISYME)
               ISYMDJ = MULD2H(ISYMEL,IT2SYM)
C
               DO 130 L = 1,NRHF(ISYML)
C
                  NEL = IT1AM(ISYME,ISYML) + NVIR(ISYME)*(L - 1) + E
C
                  DO 140 ISYMI = 1,NSYM
C
                     ISYMEI = MULD2H(ISYMI,ISYME)
                     ISYMDK = MULD2H(ISYMEI,IC2SYM)
                     ISYJKL = MULD2H(ISYMI,ISYINT)
                     ISYMJK = MULD2H(ISYJKL,ISYML)
C
                     DO 150 ISYMJ = 1,NSYM
C
                        ISYMD = MULD2H(ISYMJ,ISYMDJ)
                        ISYMK = MULD2H(ISYMJ,ISYMJK)
C
                        LET2SM = NRHF(ISYMJ)*NVIR(ISYMD)
C
                        IF (LWORK .LT. LET2SM) THEN
                           CALL QUIT('Insufficient work space in CC_MI')
                        ENDIF
C
C-----------------------------------------------------
C                       Copy T1-amplitude out of T2AM.
C-----------------------------------------------------
C
                        IF (ISYMDJ .EQ. ISYMEL) THEN
C
                           DO 160 D = 1,NVIR(ISYMD)
C
                              DO 170 J = 1,NRHF(ISYMJ)
C
                                 NDJ   = IT1AM(ISYMD,ISYMJ)
     &                                 + NVIR(ISYMD)*(J - 1) + D
                                 NDJEL = IT2AM(ISYMDJ,ISYMEL)
     &                                 + INDEX(NDJ,NEL)
                                 NJD   = NRHF(ISYMJ)*(D - 1) + J
C
                                 WORK(NJD) = T2AM(NDJEL)
C
  170                         CONTINUE
  160                      CONTINUE
C
                        ELSE IF (ISYMDJ .LT. ISYMEL) THEN
C
                           DO 180 J = 1,NRHF(ISYMJ)
C
                              NDJ   = IT1AM(ISYMD,ISYMJ)
     &                              + NVIR(ISYMD)*(J - 1) + 1
                              NDJEL = IT2AM(ISYMDJ,ISYMEL)
     &                              + NT1AM(ISYMDJ)*(NEL - 1) + NDJ
C
                              CALL DCOPY(NVIR(ISYMD),T2AM(NDJEL),1,
     &                                   WORK(J),NRHF(ISYMJ))
C
  180                      CONTINUE
C
                        ELSE IF (ISYMDJ .GT. ISYMEL) THEN
C
                           DO 190 D = 1,NVIR(ISYMD)
C
                              DO 200 J = 1,NRHF(ISYMJ)
C
                                 NDJ   = IT1AM(ISYMD,ISYMJ)
     &                                 + NVIR(ISYMD)*(J - 1) + D
                                 NELDJ = IT2AM(ISYMEL,ISYMDJ)
     &                                 + NT1AM(ISYMEL)*(NDJ - 1) + NEL
                                 NJD   = NRHF(ISYMJ)*(D - 1) + J
C
                                 WORK(NJD) = T2AM(NELDJ)
C
  200                         CONTINUE
  190                      CONTINUE
C
                        ENDIF
C
C----------------------------------------------------------------------
C                       Contraction of T2AM & CTR2 & storage in result.
C----------------------------------------------------------------------
C
                        DO 210 I = 1,NRHF(ISYMI)
C
                           NEI   = IT1AM(ISYME,ISYMI)
     &                           + NVIR(ISYME)*(I - 1) + E
                           KOFF2 = IT2SQ(ISYMDK,ISYMEI)
     &                           + NT1AM(ISYMDK)*(NEI - 1)
     &                           + IT1AM(ISYMD,ISYMK) + 1
                           KOFF3 = I3ORHF(ISYJKL,ISYMI)
     &                           + NMAIJK(ISYJKL)*(I - 1)
     &                           + IMAIJK(ISYMJK,ISYML)
     &                           + NMATIJ(ISYMJK)*(L - 1)
     &                           + IMATIJ(ISYMJ,ISYMK) + 1
C
                           NTOTJ = MAX(NRHF(ISYMJ),1)
                           NTOTD = MAX(NVIR(ISYMD),1)
C
                           CALL DGEMM('N','N',NRHF(ISYMJ),NRHF(ISYMK),
     &                                NVIR(ISYMD),ONE,WORK,NTOTJ,
     &                                CTR2(KOFF2),NTOTD,ONE,
     &                                XMINT(KOFF3),NTOTJ)
C
  210                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_22am */
      SUBROUTINE CC_22AM(RHO2,XIAJB,XMINT,ISYMIM,WORK,LWORK)
C
C     Written by Asger Halkier 16/11 - 1995
C
C     Version: 1.0
C
C     Purpose: To calculate the 22A-prime contribution to RHO2!
C
C     Ove Christiansen 18-9-1996:
C          ISYMIM is symmetry of intermediate.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO2(*), XIAJB(*), XMINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYRES = ISYMIM 
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMIN = MULD2H(ISYMJ,ISYRES)
C
         IF (NRHF(ISYMJ) .EQ. 0) GOTO 100
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYRES)
C
               IF (ISYMAI .GT. ISYMBJ) GOTO 120
               IF (NVIR(ISYMB) .EQ. 0) GOTO 120
C
C----------------------------------------
C              Work space allocation one.
C----------------------------------------
C
               KSCR1 = 1
               KEND1 = KSCR1 + NT1AM(ISYMAI)
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient work space in CC_22AM')
               ENDIF
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KSCR1),NT1AM(ISYMAI))
C
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  DO 140 ISYMN = 1,NSYM
C
                     ISYMBN = MULD2H(ISYMN,ISYMB)
                     ISYMAM = MULD2H(ISYMBN,ISYMOP)
                     ISYMMI = MULD2H(ISYMN,ISYMIN)
C
C----------------------------------------------
C                    Work space allocation two.
C----------------------------------------------
C
                     KSCR2 = KEND1
                     KEND2 = KSCR2 + NT1AM(ISYMAM)
                     LWRK2 = LWORK - KEND2
C
                     IF (LWRK2 .LT. 0) THEN
                        CALL QUIT('Insufficient work space in CC_22AM')
                     ENDIF
C
                     DO 150 N = 1,NRHF(ISYMN)
C
                        NBN = IT1AM(ISYMB,ISYMN)
     &                      + NVIR(ISYMB)*(N - 1) + B
C
C-------------------------------------------------------
C                       Copy submatrix out of integrals.
C-------------------------------------------------------
C
                        DO 160 NAM = 1,NT1AM(ISYMAM)
C
                           NAMBN = IT2AM(ISYMAM,ISYMBN)
     &                           + INDEX(NAM,NBN)
C
                           WORK(KSCR2 + NAM - 1) = XIAJB(NAMBN)
C
  160                   CONTINUE
C
C--------------------------------------------------------------------
C                       Contraction of integrals with M-intermediate.
C--------------------------------------------------------------------
C
                        DO 170 ISYMA = 1,NSYM
C
                           ISYMM = MULD2H(ISYMA,ISYMAM)
                           ISYMI = MULD2H(ISYMA,ISYMAI)
C
                           KOFF1 = KSCR2 + IT1AM(ISYMA,ISYMM)
                           KOFF2 = I3ORHF(ISYMIN,ISYMJ)
     &                           + NMAIJK(ISYMIN)*(J - 1)
     &                           + IMAIJK(ISYMMI,ISYMN)
     &                           + NMATIJ(ISYMMI)*(N - 1)
     &                           + IMATIJ(ISYMM,ISYMI) + 1
                           KOFF3 = KSCR1 + IT1AM(ISYMA,ISYMI)
C
                           NTOTA = MAX(NVIR(ISYMA),1)
                           NTOTM = MAX(NRHF(ISYMM),1)
C
                           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     &                                NRHF(ISYMM),ONE,WORK(KOFF1),
     &                                NTOTA,XMINT(KOFF2),NTOTM,ONE,
     &                                WORK(KOFF3),NTOTA)
C
  170                   CONTINUE
C
  150                CONTINUE
  140             CONTINUE
C
C-----------------------------------------------
C                 Storage in RHO2 result vector.
C-----------------------------------------------
C
                  IF (ISYMAI .EQ. ISYMBJ) THEN
C
                     DO 180 NAI = 1,NBJ
C
                        NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                        + INDEX(NAI,NBJ)
C
                        RHO2(NAIBJ) = RHO2(NAIBJ)
     &                              + WORK(KSCR1 + NAI - 1)
C
  180                CONTINUE
C
                  ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
                     KOFF4 = IT2AM(ISYMAI,ISYMBJ)
     &                     + NT1AM(ISYMAI)*(NBJ - 1) + 1
C
                     CALL DAXPY(NT1AM(ISYMAI),ONE,WORK(KSCR1),1,
     &                          RHO2(KOFF4),1)
C
                  ENDIF
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
