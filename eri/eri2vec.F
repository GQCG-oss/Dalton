C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id: eri2vec.F,v 1.1.1.1 2001-02-08 13:33:13 hjj Exp $
C
#include <single.h>
C  /* Deck ericor */
      SUBROUTINE ERICOR(COORAO,COORTR,IPNTCR,CORBCH,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      INTEGER PQ, R, S, T, RST, X
      DIMENSION COORTR(NPQBCS,NPRFPQ,3,4),
     &          COORAO(MLTPX,NPQBCS,NPRFPQ,3,4),
     &          XB(8**3), YB(8**3), ZB(8**3),
     &          XC(8**3), YC(8**3), ZC(8**3),
     &          XD(8**3), YD(8**3), ZD(8**3),
     &          IPNTCR(MAXBCH,4), CORBCH(NPQBCS,3)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <symmet.h>
#include <ibtfun.h>
      XPAR(I) = PT(IBTAND(ISYMAX(1,1),I))
      YPAR(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZPAR(I) = PT(IBTAND(ISYMAX(3,1),I))
C
C     Coordinates for identity operations R, S, and T
C     ===============================================
C
      DO 100 IATOM = 1,4
#if defined (VAR_CRY)
         CALL GATHER(NPQBCS,CORBCH(1,1),CORXBT,IPNTCR(1,IATOM))
         CALL GATHER(NPQBCS,CORBCH(1,2),CORYBT,IPNTCR(1,IATOM))
         CALL GATHER(NPQBCS,CORBCH(1,3),CORZBT,IPNTCR(1,IATOM))
#else
         DO 110 PQ = 1, NPQBCS
            CORBCH(PQ,1) = CORXBT(IPNTCR(PQ,IATOM))
            CORBCH(PQ,2) = CORYBT(IPNTCR(PQ,IATOM))
            CORBCH(PQ,3) = CORZBT(IPNTCR(PQ,IATOM))
  110    CONTINUE
#endif
         IF (MLTPX .EQ. 1) THEN
            DO 120 I  = 1, NPRFPQ
            DO 120 PQ = 1, NPQBCS
               COORAO(1,PQ,I,1,IATOM) = CORBCH(PQ,1)
               COORAO(1,PQ,I,2,IATOM) = CORBCH(PQ,2)
               COORAO(1,PQ,I,3,IATOM) = CORBCH(PQ,3)
  120       CONTINUE
         ELSE
            DO 130 I  = 1, NPRFPQ
            DO 130 PQ = 1, NPQBCS
               COORTR(PQ,I,1,IATOM) = CORBCH(PQ,1)
               COORTR(PQ,I,2,IATOM) = CORBCH(PQ,2)
               COORTR(PQ,I,3,IATOM) = CORBCH(PQ,3)
  130       CONTINUE
         END IF
  100 CONTINUE
C
C     Coordinates for remaining operations R, S, and T
C     ================================================
C
      IF (MLTPX.GT.1) THEN
         RST = 0
         DO 300 T = 0, MAXOPR
         IF (IBTAND(T,ISTBLT) .EQ. 0) THEN
            XC0 = XPAR(T)
            YC0 = YPAR(T)
            ZC0 = ZPAR(T)
            DO 310 S = 0, MAXOPR
            IF (IBTAND(S,ISTBLS) .EQ. 0) THEN
               XD0 = XPAR(IBTXOR(S,T))
               YD0 = YPAR(IBTXOR(S,T))
               ZD0 = ZPAR(IBTXOR(S,T))
               DO 320 R = 0, MAXOPR
               IF (IBTAND(R,ISTBLR) .EQ. 0) THEN
                  RST = RST + 1
                  XB(RST) = XPAR(R)
                  YB(RST) = YPAR(R)
                  ZB(RST) = ZPAR(R)
                  XC(RST) = XC0
                  YC(RST) = YC0
                  ZC(RST) = ZC0
                  XD(RST) = XD0
                  YD(RST) = YD0
                  ZD(RST) = ZD0
               END IF
  320          CONTINUE
            END IF
  310       CONTINUE
         END IF
  300    CONTINUE
C
         DO 400 IJ     = 1, NPRFPQ
         DO 400 PQ     = 1, NPQBCS
         DO 400 RST    = 1, MLTPX
            COORAO(RST,PQ,IJ,1,1) = COORTR(PQ,IJ,1,1)
            COORAO(RST,PQ,IJ,2,1) = COORTR(PQ,IJ,2,1)
            COORAO(RST,PQ,IJ,3,1) = COORTR(PQ,IJ,3,1)
C
            COORAO(RST,PQ,IJ,1,2) = COORTR(PQ,IJ,1,2)*XB(RST)
            COORAO(RST,PQ,IJ,2,2) = COORTR(PQ,IJ,2,2)*YB(RST)
            COORAO(RST,PQ,IJ,3,2) = COORTR(PQ,IJ,3,2)*ZB(RST)
C
            COORAO(RST,PQ,IJ,1,3) = COORTR(PQ,IJ,1,3)*XC(RST)
            COORAO(RST,PQ,IJ,2,3) = COORTR(PQ,IJ,2,3)*YC(RST)
            COORAO(RST,PQ,IJ,3,3) = COORTR(PQ,IJ,3,3)*ZC(RST)
C
            COORAO(RST,PQ,IJ,1,4) = COORTR(PQ,IJ,1,4)*XD(RST)
            COORAO(RST,PQ,IJ,2,4) = COORTR(PQ,IJ,2,4)*YD(RST)
            COORAO(RST,PQ,IJ,3,4) = COORTR(PQ,IJ,3,4)*ZD(RST)
  400    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL TITLER('Output from ERICOR','*',103)
         CALL HEADER('Coordinates A in ERICOR',-1)
         CALL OUTPUT(COORAO(1,1,1,1,1),1,NPPX,1,3,NPPX,3,1,LUPRI)
         CALL HEADER('Coordinates B in ERICOR',-1)
         CALL OUTPUT(COORAO(1,1,1,1,2),1,NPPX,1,3,NPPX,3,1,LUPRI)
         CALL HEADER('Coordinates C in ERICOR',-1)
         CALL OUTPUT(COORAO(1,1,1,1,3),1,NPPX,1,3,NPPX,3,1,LUPRI)
         CALL HEADER('Coordinates D in ERICOR',-1)
         CALL OUTPUT(COORAO(1,1,1,1,4),1,NPPX,1,3,NPPX,3,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck eriexp */
      SUBROUTINE ERIEXP(FACINT,EXPAO,EXPTR,RODPE1,RODPE2,RODPF1,RODPF2,
     &                  IPNTPP,ROD1,ROD2,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      INTEGER AB, CD, P, Q, R, S, T, RST, PQ
      DIMENSION
     &          EXPTR (NPQBCS,NPRFPQ,4),
     &          EXPAO (MLTPX,NPQBCS,NPRFPQ,4),
     &          FACINT(MLTPR,MLTPS,MLTPT,NPQBCS,NPRFPQ),
     &          RODPE1(NODPP1,  2), RODPE2(NODPP2,  2),
     &          RODPF1(MAXOPR+1,*), RODPF2(MAXOPR+1,*),
     &          ROD1(NPQBCS,NPRFAB,3), ROD2(NPQBCS,NPRFCD,3),
     &          IPNTPP(MAXBCH,2)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <symmet.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERIEXP','*',103)
      END IF
C
C     Exponents
C     =========
C
      IF (MLTPX .EQ. 1) THEN
         IJ = 0
         DO 400 I = 1, NPRFAB
         DO 400 J = 1, NPRFCD
            IJ = IJ + 1
            DO 410 PQ = 1, NPQBCS
               AB = IPNTPP(PQ,1) + I - 1
               CD = IPNTPP(PQ,2) + J - 1
               EXPAO(1,PQ,IJ,1) = RODPE1(AB,1)
               EXPAO(1,PQ,IJ,2) = RODPE1(AB,2)
               EXPAO(1,PQ,IJ,3) = RODPE2(CD,1)
               EXPAO(1,PQ,IJ,4) = RODPE2(CD,2)
  410       CONTINUE
  400    CONTINUE
      ELSE
         IJ = 0
         DO 500 I = 1, NPRFAB
         DO 500 J = 1, NPRFCD
            IJ = IJ + 1
            DO 510 PQ = 1, NPQBCS
               AB = IPNTPP(PQ,1) + I - 1
               CD = IPNTPP(PQ,2) + J - 1
               EXPTR(PQ,IJ,1) = RODPE1(AB,1)
               EXPTR(PQ,IJ,2) = RODPE1(AB,2)
               EXPTR(PQ,IJ,3) = RODPE2(CD,1)
               EXPTR(PQ,IJ,4) = RODPE2(CD,2)
  510       CONTINUE
  500    CONTINUE
C
         DO 520 IATOM  = 1, 4
         DO 520 IJ     = 1, NPRFPQ
         DO 520 PQ     = 1, NPQBCS
         DO 520 RST    = 1, MLTPX
            EXPAO(RST,PQ,IJ,IATOM) = EXPTR(PQ,IJ,IATOM)
  520    CONTINUE
      END IF
C
C     Preexponential and other overall factors
C     ========================================
C
      SFAC = FMULT(IBTAND(IBTAND(ISTBLA,ISTBLB),IBTAND(ISTBLC,ISTBLD)))
      IF (MLTPR*MLTPS .GT. 1) THEN
         IJ = 0
         DO 600 I = 1, NPRFAB
         DO 600 J = 1, NPRFCD
            IJ = IJ + 1
            DO 700 PQ = 1, NPQBCS
               AB = IPNTPP(PQ,1) + I - 1
               CD = IPNTPP(PQ,2) + J - 1
               DO 800 S = 1, MLTPS
               DO 800 R = 1, MLTPR
                  FACINT(R,S,1,PQ,IJ)=SFAC*RODPF1(R,AB)*RODPF2(S,CD)
  800          CONTINUE
  700       CONTINUE
  600    CONTINUE
      ELSE
         IJ = 0
         DO 610 I = 1, NPRFAB
         DO 610 J = 1, NPRFCD
            IJ = IJ + 1
            DO 710 PQ = 1, NPQBCS
               AB = IPNTPP(PQ,1) + I - 1
               CD = IPNTPP(PQ,2) + J - 1
               FACINT(1,1,1,PQ,IJ)=SFAC*RODPF1(1,AB)*RODPF2(1,CD)
  710       CONTINUE
  610    CONTINUE
      END IF
C
      IF      (MLTPT .EQ. 2) THEN
         DO 902 IJ = 1, NPRFPQ
         DO 902 PQ = 1, NPQBCS
         DO 902 S  = 1, MLTPS
         DO 902 R  = 1, MLTPR
            FACINT(R,S,2,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
  902    CONTINUE
      ELSE IF (MLTPT .EQ. 4) THEN
         DO 904 IJ = 1, NPRFPQ
         DO 904 PQ = 1, NPQBCS
         DO 904 S  = 1, MLTPS
         DO 904 R  = 1, MLTPR
            FACINT(R,S,2,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,3,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,4,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
  904    CONTINUE
      ELSE IF (MLTPT .EQ. 8) THEN
         DO 908 IJ = 1, NPRFPQ
         DO 908 PQ = 1, NPQBCS
         DO 908 S  = 1, MLTPS
         DO 908 R  = 1, MLTPR
            FACINT(R,S,2,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,3,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,4,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,5,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,6,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,7,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
            FACINT(R,S,8,PQ,IJ) = FACINT(R,S,1,PQ,IJ)
  908    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('Exponents in ERIEXP',-1)
         CALL OUTPUT(EXPAO,1,NPPX,1,4,NPPX,4,1,LUPRI)
         CALL HEADER('FACINT in ERIEXP',-1)
         CALL OUTPUT(FACINT,1,NPQBCS*MLTPX,1,NPRFPQ,
     &               NPQBCS*MLTPX,NPRFPQ,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck ericox */
      SUBROUTINE ERICOX(COORAO,IPNTCR,CORBCH,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      INTEGER PQ, R, S, T, RST, X
      DIMENSION COORAO(NPPS,MLTPX,3,4), IPNTCR(MAXBCH,4),
     &          CORBCH(NPQBCS,3)
#include <cbieri.h>
#include <odclss.h>
#include <ericom.h>
#include <hertop.h>
#include <aobtch.h>
#include <symmet.h>
#include <ibtfun.h>
      XPAR(I) = PT(IBTAND(ISYMAX(1,1),I))
      YPAR(I) = PT(IBTAND(ISYMAX(2,1),I))
      ZPAR(I) = PT(IBTAND(ISYMAX(3,1),I))
C
      DO 100 X = 1,4
#if defined (VAR_CRY)
         CALL GATHER(NPQBCS,CORBCH(1,1),CORXBT,IPNTCR(1,X))
         CALL GATHER(NPQBCS,CORBCH(1,2),CORYBT,IPNTCR(1,X))
         CALL GATHER(NPQBCS,CORBCH(1,3),CORZBT,IPNTCR(1,X))
#else
         DO 110 PQ = 1, NPQBCS
            CORBCH(PQ,1) = CORXBT(IPNTCR(PQ,X))
            CORBCH(PQ,2) = CORYBT(IPNTCR(PQ,X))
            CORBCH(PQ,3) = CORZBT(IPNTCR(PQ,X))
  110    CONTINUE
#endif
         IJKL = 0
         DO 120 I = 1, NPRFPQ
         DO 120 PQ = 1,NPQBCS
            IJKL = IJKL + 1
            COORAO(IJKL,1,1,X) = CORBCH(PQ,1)
            COORAO(IJKL,1,2,X) = CORBCH(PQ,2)
            COORAO(IJKL,1,3,X) = CORBCH(PQ,3)
  120    CONTINUE
  100 CONTINUE
C
      IF (SYMTRA) THEN
C
C        Coordinates A
C        =============
C
         DO 200 X   = 1, 3
         DO 200 RST = 2, MLTPX
         DO 200 I   = 1, NPPS
            COORAO(I,RST,X,1) = COORAO(I,1,X,1)
  200    CONTINUE
C
C        Coordinates A, B and C
C        ======================
C
         RST = 0
         DO 300 T = 0, MAXOPR
         IF (IBTAND(T,ISTBLT) .EQ. 0) THEN
            XC = XPAR(T)
            YC = YPAR(T)
            ZC = ZPAR(T)
            DO 310 S = 0, MAXOPR
            IF (IBTAND(S,ISTBLS) .EQ. 0) THEN
               XD = XPAR(IBTXOR(S,T))
               YD = YPAR(IBTXOR(S,T))
               ZD = ZPAR(IBTXOR(S,T))
               DO 320 R = 0, MAXOPR
               IF (IBTAND(R,ISTBLR) .EQ. 0) THEN
                  RST = RST + 1
                  IF (RST .GT. 1) THEN
                     XB = XPAR(R)
                     YB = YPAR(R)
                     ZB = ZPAR(R)
                     DO 400 I = 1, NPPS
                         COORAO(I,RST,1,2) = XB*COORAO(I,1,1,2)
                         COORAO(I,RST,2,2) = YB*COORAO(I,1,2,2)
                         COORAO(I,RST,3,2) = ZB*COORAO(I,1,3,2)
C
                         COORAO(I,RST,1,3) = XC*COORAO(I,1,1,3)
                         COORAO(I,RST,2,3) = YC*COORAO(I,1,2,3)
                         COORAO(I,RST,3,3) = ZC*COORAO(I,1,3,3)
C
                         COORAO(I,RST,1,4) = XD*COORAO(I,1,1,4)
                         COORAO(I,RST,2,4) = YD*COORAO(I,1,2,4)
                         COORAO(I,RST,3,4) = ZD*COORAO(I,1,3,4)
  400                CONTINUE
C
                  END IF
               END IF
  320          CONTINUE
            END IF
  310       CONTINUE
         END IF
  300    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL TITLER('Output from ERICOX','*',103)
         CALL HEADER('Coordinates A in ERICOX',-1)
         CALL OUTPUT(COORAO(1,1,1,1),1,NPPX,1,3,NPPX,3,1,LUPRI)
         CALL HEADER('Coordinates B in ERICOX',-1)
         CALL OUTPUT(COORAO(1,1,1,2),1,NPPX,1,3,NPPX,3,1,LUPRI)
         CALL HEADER('Coordinates C in ERICOX',-1)
         CALL OUTPUT(COORAO(1,1,1,3),1,NPPX,1,3,NPPX,3,1,LUPRI)
         CALL HEADER('Coordinates D in ERICOX',-1)
         CALL OUTPUT(COORAO(1,1,1,4),1,NPPX,1,3,NPPX,3,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck eriexx */
      SUBROUTINE ERIEXX(RODPP1,RODPP2,IPNTPP,FACINT,EXPAO,ROD1,ROD2,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
#include <ericom.h>
      INTEGER AB, CD, P, Q, R, S, T, RST
      DIMENSION RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          ROD1(NPQBCS,NPRFAB,3), ROD2(NPQBCS,NPRFCD,3),
     &          IPNTPP(MAXBCH,2),
     &          FACINT(NPPS,MLTPR,MLTPS,MLTPT),
     &          EXPAO(NPPS,MLTPX,4)
#include <cbieri.h>
#include <odclss.h>
#include <hertop.h>
#include <aobtch.h>
C
      IF (IPRINT .GT. 10) THEN
         CALL TITLER('Output from ERIEXX','*',103)
      END IF
C
#if defined (VAR_CRY)
      DO 100 I = 1,NPRFAB
         CALL GATHER(NPQBCS,ROD1(1,I,1),RODPP1(I,1),IPNTPP(1,1))
         CALL GATHER(NPQBCS,ROD1(1,I,2),RODPP1(I,2),IPNTPP(1,1))
         CALL GATHER(NPQBCS,ROD1(1,I,3),RODPP1(I,3),IPNTPP(1,1))
  100 CONTINUE
      DO 110 I = 1,NPRFCD
         CALL GATHER(NPQBCS,ROD2(1,I,1),RODPP1(I,1),IPNTPP(1,2))
         CALL GATHER(NPQBCS,ROD2(1,I,2),RODPP1(I,2),IPNTPP(1,2))
         CALL GATHER(NPQBCS,ROD2(1,I,3),RODPP1(I,3),IPNTPP(1,2))
  110 CONTINUE
      IJKL = 0
      DO 200 I = 1, NPRFAB
         DO 210 J = 1, NPRFCD
            DO 220 K = 1, NPQBCS
               IJKL = IJKL + 1
               EXPAO (IJKL,1,1)   = ROD1(K,I,1)
               EXPAO (IJKL,1,2)   = ROD1(K,I,2)
               FACINT(IJKL,1,1,1) = ROD1(K,I,3)*ROD2(K,J,3)
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
      DO 300 I = 1, NPRFAB
         KOFF = NPQBCS*NPRFCD*(I - 1) + 1
         CALL DCOPY(NPRFCD*NPQBCS,ROD2(1,1,1),1,EXPAO(KOFF,1,3),1)
         CALL DCOPY(NPRFCD*NPQBCS,ROD2(1,1,2),1,EXPAO(KOFF,1,4),1)
  300 CONTINUE
#else
      IJKL = 0
      DO 400 I = 1, NPRFAB
      DO 400 J = 1, NPRFCD
      DO 400 K = 1, NPQBCS
         AB = IPNTPP(K,1) + I - 1
         CD = IPNTPP(K,2) + J - 1
         IJKL = IJKL + 1
         EXPAO (IJKL,1,1)   = RODPP1(AB,1)
         EXPAO (IJKL,1,2)   = RODPP1(AB,2)
         EXPAO (IJKL,1,3)   = RODPP2(CD,1)
         EXPAO (IJKL,1,4)   = RODPP2(CD,2)
         FACINT(IJKL,1,1,1) = RODPP1(AB,3)*RODPP2(CD,3)
  400 CONTINUE
      IF (SYMTRA) THEN
         DO 500 R = 1, MLTPR
         DO 500 S = 1, MLTPS
            IJKL = 0
            DO 510 I = 1, NPRFAB
            DO 510 J = 1, NPRFCD
            DO 510 K = 1, NPQBCS
               AB = IPNTPP(K,1) + I - 1
               CD = IPNTPP(K,2) + J - 1
               IJKL = IJKL + 1
               FACINT(IJKL,R,S,1) = RODPP1(AB,R+2)*RODPP2(CD,S+2)
  510       CONTINUE
            DO 520 T = 2, MLTPT
            DO 520 I = 1, NPPS
               FACINT(I,R,S,T) = FACINT(I,R,S,1)
  520       CONTINUE
  500    CONTINUE
      END IF
#endif
C
      IF (SYMTRA) THEN
         DO 600 IATOM = 1, 4
         DO 600 RST   = 2, MLTPX
         DO 600 I     = 1, NPPS
            EXPAO(I,RST,IATOM) = EXPAO(I,1,IATOM)
  600    CONTINUE
      END IF
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('FACINT in ERIEXX',-1)
         CALL OUTPUT(FACINT,1,NPPS,1,MLTPX,NPPS,MLTPX,1,LUPRI)
         CALL HEADER('Exponents in ERIEXX',-1)
         CALL OUTPUT(EXPAO,1,NPPX,1,4,NPPX,4,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck erivcs */
      SUBROUTINE ERIVCS(EXPAO,EXPPQ,PQINV,COORAO,COORPQ,PQXYZ,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      INTEGER X
      DIMENSION EXPAO(NPPX,4), COORAO(NPPX,3,4),
     &          EXPPQ(NPPX,2), PQINV(NPPX,2), COORPQ(NPPX,3,4),
     &          PQXYZ(NPPX,3)
#include <ericom.h>
#include <erithr.h>
#include <hertop.h>
#include <ibtfun.h>
C
C     Exponents and coordinates
C     =========================
C
      DO 100 I = 1, NPPX
         EXPPQ(I,1) = EXPAO(I,1) + EXPAO(I,2)
         EXPPQ(I,2) = EXPAO(I,3) + EXPAO(I,4)
         PQINV(I,1) = D1/EXPPQ(I,1)
         PQINV(I,2) = D1/EXPPQ(I,2)
  100 CONTINUE
C
      DO 200 X = 1, 3
      DO 200 I = 1, NPPX
         COORPQ(I,X,1) = (EXPAO(I,1)*COORAO(I,X,1)
     &                  + EXPAO(I,2)*COORAO(I,X,2))*PQINV(I,1)
         COORPQ(I,X,3) = (EXPAO(I,3)*COORAO(I,X,3)
     &                  + EXPAO(I,4)*COORAO(I,X,4))*PQINV(I,2)
  200 CONTINUE
C
      DO 300 I = 1, NPPX
         PQXYZ(I,1) = COORPQ(I,1,1) - COORPQ(I,1,3)
         PQXYZ(I,2) = COORPQ(I,2,1) - COORPQ(I,2,3)
         PQXYZ(I,3) = COORPQ(I,3,1) - COORPQ(I,3,3)
  300 CONTINUE
C
C     A - B
C
      DO 400 X = 1, 3
      DO 400 I = 1, NPPX
         COORPQ(I,X,2) = COORAO(I,X,1) - COORAO(I,X,2)
  400 CONTINUE
C
C     C - D
C
      DO 500 X = 1, 3
      DO 500 I = 1, NPPX
         COORPQ(I,X,4) = COORAO(I,X,3) - COORAO(I,X,4)
  500 CONTINUE
C
C     Local symmetries
C     ================
C
      CALL LOCSYM(PQXYZ,        NPPX,IPQ0,IPQXYZ)
      CALL LOCSYM(COORPQ(1,1,2),NPPX,IAB0,IABXYZ)
      CALL LOCSYM(COORPQ(1,1,4),NPPX,ICD0,ICDXYZ)
C
      IHHXYZ = IPQXYZ
      IHCXYZ = IBTAND(IPQXYZ,IABXYZ)
      ICHXYZ = IBTAND(IPQXYZ,ICDXYZ)
      ICCXYZ = IBTAND(IPQXYZ,IBTAND(IABXYZ,ICDXYZ))
C
      IF (IPRINT .GT. 25) THEN
         CALL HEADER('EXPPQ in ERIVCS',1)
         CALL OUTPUT(EXPPQ,1,NPPX,1,2,NPPX,2,1,LUPRI)
         CALL HEADER('PQINV in ERIVCS',1)
         CALL OUTPUT(PQINV,1,NPPX,1,2,NPPX,2,1,LUPRI)
C
         CALL HEADER('COORPQ(*,*,1) (P) in ERIVCS',1)
         CALL OUTPUT(COORPQ(1,1,1),1,NPPX,1,3,NPPX,3,1,LUPRI)
         IF (JMAXAB .GT. 0) THEN
            CALL HEADER('COORPQ(*,*,2) (A-B) in ERIVCS',1)
            CALL OUTPUT(COORPQ(1,1,2),1,NPPX,1,3,NPPX,3,1,LUPRI)
         END IF
C

         CALL OUTPUT(COORPQ(1,1,3),1,NPPX,1,3,NPPX,3,1,LUPRI)
         IF (JMAXCD .GT. 0) THEN
            CALL HEADER('COORPQ(*,*,4) (C-D) in ERIVCS',1)
            CALL OUTPUT(COORPQ(1,1,4),1,NPPX,1,3,NPPX,3,1,LUPRI)
         END IF
C
         CALL HEADER('PQXYZ in ERIVCS',1)
         CALL OUTPUT(PQXYZ,1,NPPX,1,3,NPPX,3,1,LUPRI)
C
         WRITE (LUPRI,'(/2X,A,E20.10)') ' THRSH', THRSH
         WRITE (LUPRI,'(/,2X,A,I5,5X,3I5)')
     &      'IPQXYZ, IPQ0(3):', IPQXYZ, IPQ0(1), IPQ0(2), IPQ0(3)
         WRITE (LUPRI,'(/,2X,A,I5,5X,3I5)')
     &      'IABXYZ, IAB0(3):', IABXYZ, IAB0(1), IAB0(2), IAB0(3)
         WRITE (LUPRI,'(/,2X,A,I5,5X,3I5)')
     &      'ICDXYZ, ICD0(3):', ICDXYZ, ICD0(1), ICD0(2), ICD0(3)
         WRITE (LUPRI,'(/,2X,A,4I5)')
     &      'IHHXYZ,IHCXYZ,ICHXYZ,ICCXYZ',
     &       IHHXYZ,IHCXYZ,ICHXYZ,ICCXYZ
C
      END IF
      RETURN
      END
C  /* Deck locsym */
      SUBROUTINE LOCSYM(VECTOR,NINNER,I0,IXYZ)
#include <implicit.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION VECTOR(NINNER,3), I0(3)
#include <erithr.h>
C
      I0(1) = 1
      I0(2) = 1
      I0(3) = 1
      DX = D0
      DY = D0
      DZ = D0
      DO 100 I = 1, NINNER
         DX = MAX(DX,ABS(VECTOR(I,1)))
         DY = MAX(DY,ABS(VECTOR(I,2)))
         DZ = MAX(DZ,ABS(VECTOR(I,3)))
  100 CONTINUE
      IF (DX .GT. THRSH) I0(1) = 0
      IF (DY .GT. THRSH) I0(2) = 0
      IF (DZ .GT. THRSH) I0(3) = 0
      IXYZ = I0(1) + 2*I0(2) + 4*I0(3)
      RETURN
      END
