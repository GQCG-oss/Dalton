C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
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
#include <single.h>
C
C  /* Deck makerv*/ 
      SUBROUTINE MAKERV(V2AM,R2AM,U2AM,S2AM,T2AM,
     *                  EV,QQ11,V11,U11,B11,W11,Q11,R11,
     *                  VS11,US11,BS11,WS11,
     *                  QS11,RS11,VT11,UT11,BT11,WT11,QT11,RT11,
     *                  NICDIM,NOCDIM,NSPAIR,NTPAIR,ES,ET,FS,FT,EVS,
     *                  CNINV,QQ2,QQ4,QQ6)            
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      REAL*8  SU, US, SV, VU, VR, RR, UR, RU, UT, TU, QQIJKL, QQIJLK,
     *        EIJAB, EKLAB, EIJ, EMN, EKL, EIJKL, EAB,
     *        SAIBJ, RAKBL, UAKBL, SAKBL, RAIBJ, VAIBJ, UAIBJ,
     *        TAIBJ, TAKBL,
     *        V11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        W11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        B11(NOCDIM,NOCDIM,NOCDIM,NOCDIM,NSPAIR),
     *        U11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        Q11(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *        R11(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
      DIMENSION QQ2(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *          QQ4(NOCDIM,NOCDIM,NOCDIM,NOCDIM),
     *          QQ6(NOCDIM,NOCDIM,NOCDIM,NOCDIM) 
      DIMENSION R2AM(*), V2AM(*), U2AM(*), S2AM(*), T2AM(*), EV(*)
      DIMENSION VS11(NSPAIR,NSPAIR)
      DIMENSION US11(NSPAIR,NSPAIR)
      DIMENSION WS11(NSPAIR,NSPAIR)
      DIMENSION BS11(NSPAIR,NSPAIR,NSPAIR)
      DIMENSION WT11(NSPAIR,NSPAIR)
      DIMENSION BT11(NSPAIR,NSPAIR,NSPAIR)
      DIMENSION QS11(NSPAIR,NSPAIR)
      DIMENSION RS11(NSPAIR,NSPAIR)
      DIMENSION VT11(NSPAIR,NSPAIR),UT11(NSPAIR,NSPAIR)
      DIMENSION QT11(NSPAIR,NSPAIR),RT11(NSPAIR,NSPAIR)
      DIMENSION ES(NSPAIR),ET(NSPAIR),FS(NSPAIR),FT(NSPAIR)
      DIMENSION EVS(NSPAIR),CNINV(NSPAIR,8)
      DIMENSION QQ11(NICDIM,NICDIM,NICDIM,NICDIM)
      LOGICAL DOA, DOB, DOAP, DOBP, ASTRIX, DONINV, LDUM
      INTEGER IDUM
#include <r12int.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      FF = D1 / SQRT(D2)
C
      DO 999 IPDD = 1, 14
C
      DOA = IPDD .LE. 6 
      DOB = .NOT. DOA               
      DONINV = IPDD .EQ. 1 .OR. IPDD .EQ. 4 .OR.
     *         IPDD .EQ. 7 .OR. IPDD .EQ. 9 .OR.
     *         IPDD .EQ. 11 .OR. IPDD .EQ. 13
      DOAP = MOD(IPDD, 3) .EQ. 0 .AND. DOA
      DOBP = IPDD .EQ. 9 .OR. IPDD .EQ. 10 .OR. IPDD .EQ. 13
     *       .OR. IPDD .EQ. 14
      ASTRIX = (IPDD .GE. 4 .AND. IPDD .LE. 6) .OR. 
     *         (IPDD .GE. 11 .AND. IPDD .LE. 14)
      IF (DOA .AND. R12NOA) GOTO 999
      IF (DOB .AND. R12NOB) GOTO 999
      IF (DOAP .AND. R12NOP) GOTO 999
      IF (DOBP .AND. NORXR) GOTO 999
      IF (IPDD .EQ. 2 .OR. IPDD .EQ. 5) THEN
         XPDD = D0
      ELSE
         XPDD = DP5
      END IF

      LU33 = 33
      CALL GPOPEN(LU33,'AUXQ12','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ(LU33,'(4E30.20)') QQ11
      CALL GPCLOSE(LU33,'KEEP')
C
      DO 205 I = 1, NOCDIM**4
         R11(I,1,1,1) = D0
         Q11(I,1,1,1) = D0
         V11(I,1,1,1) = D0
         W11(I,1,1,1) = D0
         U11(I,1,1,1) = D0
  205 CONTINUE
      DO 206 I = 1, NSPAIR*NOCDIM**4
         B11(I,1,1,1,1) = D0
  206 CONTINUE
      CALL DZERO(ES,NSPAIR)
      CALL DZERO(ET,NSPAIR)
C
C     CONSTRUCT PRODUCTS (IA|JB)*(KA|LB)
C
      E2 = D0
      E2S = D0
      E2T = D0
      DO 101 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         DO 102 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 103 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
C
C              COMPUTE MP2 ENERGY
C
               DO 201 I = 1, NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  DO 202 J = 1, NRHF(ISYMJ)
                     KOFFJ = IRHF(ISYMJ) + J           
                     DO 203 A = NRHFS(ISYMA) + 1, NORB1(ISYMA)
                        ISYMJA = MULD2H(ISYMJ,ISYMA)
                        KOFFA = IVIR(ISYMA) + A 
                        DO 204 B = NRHFS(ISYMB) + 1, NORB1(ISYMB)
                           ISYMIB = MULD2H(ISYMI,ISYMB)
                           KOFFB = IVIR(ISYMB) + B
                           NAI = IT1AM(ISYMA,ISYMI) +
     *                           NVIR(ISYMA)*(I-1) + A
                           NBJ = IT1AM(ISYMB,ISYMJ) +
     *                           NVIR(ISYMB)*(J-1) + B
                           NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                             INDEX(NAI,NBJ)
                           NAJ = IT1AM(ISYMA,ISYMJ) +
     *                           NVIR(ISYMA)*(J-1) + A
                           NBI = IT1AM(ISYMB,ISYMI) +
     *                           NVIR(ISYMB)*(I-1) + B
                           NAJBI = IT2AM(ISYMJA,ISYMIB) +
     *                             INDEX(NAJ,NBI)
                           VAIBJ = V2AM(NAIBJ)
                           VAJBI = V2AM(NAJBI)
                           VV = VAIBJ * (D2 * VAIBJ - VAJBI)
                           DENOM = D1 / (EV(KOFFI) + EV(KOFFJ) -
     *                                   EV(KOFFA) - EV(KOFFB))
                           E2 = E2 + VV * DENOM                     
                           IJ = INDEX(KOFFI,KOFFJ)
                           VS = (VAIBJ + VAJBI)**2
                           VT = (VAIBJ - VAJBI)**2
                           ES(IJ) = ES(IJ) + VS * DENOM
                           ET(IJ) = ET(IJ) + VT * DENOM
  204                   CONTINUE
  203                CONTINUE
  202             CONTINUE
  201          CONTINUE
C
               DO 104 ISYMKA = 1, NSYM
                  ISYMLB = ISYMKA 
                  ISYMK = MULD2H(ISYMA,ISYMKA)
                  ISYML = MULD2H(ISYMB,ISYMLB)
                  DO 105 I = 1, NRHF(ISYMI)
                     KOFFI = IRHF(ISYMI) + I
                     DO 106 K = 1, NRHF(ISYMK)
                        KOFFK = IRHF(ISYMK) + K
                        DO 107 A = 1, NVIR(ISYMA)
                           KOFFA = IVIR(ISYMA) + A
                           NAI = IT1AM(ISYMA,ISYMI) +
     *                           NVIR(ISYMA)*(I-1) + A
                           NAK = IT1AM(ISYMA,ISYMK) +
     *                           NVIR(ISYMA)*(K-1) + A
                           DO 108 J = 1, NRHF(ISYMJ)
                              KOFFJ = IRHF(ISYMJ) + J
                              DO 109 L = 1, NRHF(ISYML)
                                 KOFFL = IRHF(ISYML) + L
                                 EIJKL = EV(KOFFI) + EV(KOFFJ)
     *                                 - EV(KOFFK) - EV(KOFFL)
                                 DO 110 B = 1, NVIR(ISYMB)
                                    KOFFB = IVIR(ISYMB) + B
                                    NBJ = IT1AM(ISYMB,ISYMJ) +
     *                                    NVIR(ISYMB)*(J-1) + B
                                    NBL = IT1AM(ISYMB,ISYML) +
     *                                    NVIR(ISYMB)*(L-1) + B
                                    NAIBJ = IT2AM(ISYMIA,ISYMJB) +
     *                                      INDEX(NAI,NBJ)
                                    NAKBL = IT2AM(ISYMKA,ISYMLB) +
     *                                      INDEX(NAK,NBL)
C
                                    EIJAB = EV(KOFFI) + EV(KOFFJ)
     *                                    - EV(KOFFA) - EV(KOFFB)
                                    EKLAB = EV(KOFFK) + EV(KOFFL)
     *                                    - EV(KOFFA) - EV(KOFFB)
C
                                    RAIBJ = R2AM(NAIBJ) 
                                    VAIBJ = V2AM(NAIBJ)
                                    UAIBJ = U2AM(NAIBJ) 
                                    SAIBJ = S2AM(NAIBJ) 
                                    TAIBJ = T2AM(NAIBJ) 
C
                                    RAKBL = R2AM(NAKBL)
                                    UAKBL = U2AM(NAKBL) 
                                    SAKBL = S2AM(NAKBL) 
                                    TAKBL = T2AM(NAKBL) 
C
                                    IF (ASTRIX) THEN
C   
C                                      SU := (Ek+El-Ea-Eb) * <ab|r12|kl>
C                                      US := (Ei+Ej-Ea-Eb) * <ab|r12|ij>
C
                                       SU = EKLAB * RAKBL
                                       US = EIJAB * RAIBJ 
                                    ELSE
C   
C                                      SU := - <ab|[T1+T2,r12]|kl> + <ab|[K1+K2,r12]|kl>
C                                      US := - <ab|[T1+T2,r12]|ij> + <ab|[K1+K2,r12]|kl>
C                                    
                                       SU = UAKBL - SAKBL
                                       US = UAIBJ - SAIBJ
                                    ENDIF
C
C                                   SV := - <ab|(F1+F2-Ei-Ej)r12|kl> =
C                                         - <ab|[F1+F2,r12]|kl> 
C                                         + (Ei+Ej-Ek-El) * <ab|r12|kl>
C
                                    SV = SU + EIJKL * RAKBL
C
                                    VU = VAIBJ * SV
                                    VR = VAIBJ * RAKBL
                                    RR = RAIBJ * RAKBL
                                    UR = UAIBJ * RAKBL + RAIBJ * UAKBL
C      
                                    IF (A .LE. NRHFS(ISYMA) .AND.
     *                                  B .LE. NRHFS(ISYMB)) THEN
C 
C                                      Sum |ij><ij|
C
                                       V11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                                 V11(KOFFI,KOFFJ,KOFFK,KOFFL) - VR
                                       U11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                                 U11(KOFFI,KOFFJ,KOFFK,KOFFL) - UR
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) - RR
C
                                    ELSE IF (A .LE. NRHFS(ISYMA) .AND.  
     *                                       B .GT. NORB1(ISYMB)) THEN
C
C                                      Sum |iq'><iq'|
C
                                       V11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 V11(KOFFI,KOFFJ,KOFFK,KOFFL) + VR
                                       U11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 U11(KOFFI,KOFFJ,KOFFK,KOFFL) + UR
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
C
                                    ELSE IF (A .GT. NORB1(ISYMA) .AND.
     *                                       B .LE. NRHFS(ISYMB)) THEN
C
C                                      Sum |p'j><p'j|
C
                                       V11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 V11(KOFFI,KOFFJ,KOFFK,KOFFL) + VR
                                       U11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 U11(KOFFI,KOFFJ,KOFFK,KOFFL) + UR
                                       Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                                 Q11(KOFFI,KOFFJ,KOFFK,KOFFL) + RR
C
                                    ELSE IF (A. GT. NRHFS(ISYMA) .AND.
     *                                       A .LE. NORB1(ISYMA) .AND.
     *                                       B .GT. NRHFS(ISYMB) .AND. 
     *                                       B .LE. NORB1(ISYMB)) THEN
C
C                                      Sum |ab><ab|
C                                                      
                                       EIJ = D1 / 
     *                                      (EV(KOFFI) + EV(KOFFJ) -
     *                                       EV(KOFFA) - EV(KOFFB)) 
                                       W11(KOFFI,KOFFJ,KOFFK,KOFFL) = 
     *                                 W11(KOFFI,KOFFJ,KOFFK,KOFFL) -
     *                                 VU * EIJ
C
                                       MN = 0
                                       DO 111 M = 1, NOCDIM 
                                          DO 112 N = 1, M
                                             MN = MN + 1
                                             EMN = D1 / (EV(M) + EV(N) -
     *                                       EV(KOFFA) - EV(KOFFB))
                                             EIJ = EV(M) + EV(N) 
     *                                           - EV(KOFFI) - EV(KOFFJ)
                                             EKL = EV(M) + EV(N) 
     *                                           - EV(KOFFK) - EV(KOFFL)
C
C                                            RU := (Em+En-Ea-Eb) * <ab|r12|kl>
C                                            UR := (Em+En-Ea-Eb) * <ab|r12|ij>
C
                                             RU = SU + EKL * RAKBL
                                             UR = US + EIJ * RAIBJ
                                             IF (DOA) THEN
                                                TU = UAKBL
                                                UT = UAIBJ
                                             ELSE IF (DOBP) THEN
                                                TU = UAKBL - SAKBL
                                                UT = UAIBJ - SAIBJ
                                             ELSE
                                                TU = UAKBL - TAKBL
                                                UT = UAIBJ - TAIBJ
                                             END IF
                                             IF (DOAP .OR. DOB) THEN
                                                TU = TU + EKL * RAKBL
                                                UT = UT + EIJ * RAIBJ
                                             END IF
                                             UU = RU * UT + TU * UR 
                                             B11(KOFFI,KOFFJ,    
     *                                       KOFFK,KOFFL,MN) =
     *                                       B11(KOFFI,KOFFJ,    
     *                                       KOFFK,KOFFL,MN) - 
     *                                       UU * EMN
  112                                     CONTINUE
  111                                  CONTINUE
C
                                    END IF
C
  110                            CONTINUE
  109                         CONTINUE
  108                      CONTINUE
  107                   CONTINUE
  106                CONTINUE
  105             CONTINUE
  104          CONTINUE
  103       CONTINUE
  102    CONTINUE
  101 CONTINUE
C
      DO 230 I=1,NSPAIR
        ES(I) = ES(I) * DP25
        ET(I) = ET(I) * DP75
        E2S = E2S + ES(I)
        E2T = E2T + ET(I)
  230 CONTINUE
C
      IJ = 0
      DO 300 I = 1, NOCDIM
         DO 301 J = 1, I
            IJ = IJ + 1
            KL = 0
            DO 302 K = 1, NOCDIM
               DO 303 L = 1, K
                  KL = KL + 1
C
                  IF (R12EOR) THEN
                     RR = -D2 *(QQ6(I,J,K,L) + QQ6(I,J,L,K))
                     RR = RR - (U11(I,J,K,L) + U11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (U11(I,J,K,L) + U11(I,J,L,K))
                     DX = D2
                  END IF 
                  US11(KL,IJ) = RR    
                  IF (I .EQ. J) US11(KL,IJ) = FF * US11(KL,IJ)
                  IF (K .EQ. L) US11(KL,IJ) = FF * US11(KL,IJ)
                  IF (IJ .EQ. KL) US11(KL,IJ) = US11(KL,IJ) - DX
                  IF (R12EOR) THEN
                     RR = -D2 *(QQ6(I,J,K,L) - QQ6(I,J,L,K))
                     RR = RR - (U11(I,J,K,L) - U11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (U11(I,J,K,L) - U11(I,J,L,K))
                     DX = D2
                  END IF 
                  UT11(KL,IJ) = RR                                 
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            UT11(KL,IJ) = UT11(KL,IJ) - DX
C
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) + QQ2(I,J,L,K))
                     RR = RR - (V11(I,J,K,L) + V11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (V11(I,J,K,L) + V11(I,J,L,K))
                     DX = D1
                  END IF 
                  VS11(KL,IJ) = RR
                  IF (I .EQ. J) VS11(KL,IJ) = FF * VS11(KL,IJ)
                  IF (K .EQ. L) VS11(KL,IJ) = FF * VS11(KL,IJ)
                  IF (IJ .EQ. KL) VS11(KL,IJ) = VS11(KL,IJ) + DX
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) - QQ2(I,J,L,K))
                     RR = RR - (V11(I,J,K,L) - V11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (V11(I,J,K,L) - V11(I,J,L,K))
                     DX = D1
                  END IF 
                  VT11(KL,IJ) = RR
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            VT11(KL,IJ) = VT11(KL,IJ) + DX
C
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) + QQ2(I,J,L,K))
                     RR = RR - (V11(I,J,K,L) + V11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (V11(I,J,K,L) + V11(I,J,L,K))
                     DX = D1
                  END IF 
                  RR = RR + (W11(I,J,K,L) + W11(I,J,L,K))
                  WS11(KL,IJ) = RR
                  IF (I .EQ. J) WS11(KL,IJ) = FF * WS11(KL,IJ)
                  IF (K .EQ. L) WS11(KL,IJ) = FF * WS11(KL,IJ) 
                  IF (IJ .EQ. KL) WS11(KL,IJ) = WS11(KL,IJ) + DX
                  IF (R12EOR) THEN
                     RR =      (QQ2(I,J,K,L) - QQ2(I,J,L,K))
                     RR = RR - (V11(I,J,K,L) - V11(I,J,L,K))
                     DX = D0
                  ELSE
                     RR = - (V11(I,J,K,L) - V11(I,J,L,K))
                     DX = D1
                  END IF
                  RR = RR + (W11(I,J,K,L) - W11(I,J,L,K))   
                  WT11(KL,IJ) = RR
                  IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                            WT11(KL,IJ) = WT11(KL,IJ) + DX
C
                  DO 210 MN = 1, NSPAIR
                     IF (R12EOR) THEN
                        RR = -D2 *(QQ6(I,J,K,L) + QQ6(I,J,L,K))
                        RR = RR - (U11(I,J,K,L) + U11(I,J,L,K))
                        DX = D0
                     ELSE
                        RR = - (U11(I,J,K,L) + U11(I,J,L,K))
                        DX = D2
                     END IF 
                     RR = RR + (B11(I,J,K,L,MN) + B11(I,J,L,K,MN))
                     BS11(KL,IJ,MN) = RR
                     IF (I .EQ. J) BS11(KL,IJ,MN) = FF * BS11(KL,IJ,MN)
                     IF (K .EQ. L) BS11(KL,IJ,MN) = FF * BS11(KL,IJ,MN) 
                     IF (IJ .EQ. KL) 
     *                             BS11(KL,IJ,MN) = BS11(KL,IJ,MN) - DX 
                     IF (R12EOR) THEN
                        RR = -D2 *(QQ6(I,J,K,L) - QQ6(I,J,L,K))
                        RR = RR - (U11(I,J,K,L) - U11(I,J,L,K))
                        DX = D0
                     ELSE
                        RR = - (U11(I,J,K,L) - U11(I,J,L,K))
                        DX = D2
                     END IF
                     RR = RR + (B11(I,J,K,L,MN) - B11(I,J,L,K,MN))
                     BT11(KL,IJ,MN) = RR
                     IF (IJ .EQ. KL .AND. I .NE. J .AND. K. NE. L) 
     *                             BT11(KL,IJ,MN) = BT11(KL,IJ,MN) - DX
  210             CONTINUE
C
                  IF (R12EOR) THEN 
                     QQIJKL = QQ4(I,J,K,L)
                     QQIJLK = QQ4(I,J,L,K)
                  ELSE
                     QQIJKL = QQ11(IQ(I),IQ(J),IQ(K),IQ(L))
                     QQIJLK = QQ11(IQ(I),IQ(J),IQ(L),IQ(K))
                  END IF
                  RR =   (QQIJKL + QQIJLK) 
     *                 - (Q11(I,J,K,L) + Q11(I,J,L,K))
                  QS11(KL,IJ) = RR
                  IF (I .EQ. J) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  IF (K .EQ. L) QS11(KL,IJ) = FF * QS11(KL,IJ)
                  RR =   (QQIJKL - QQIJLK)
     *                 - (Q11(I,J,K,L) - Q11(I,J,L,K))
                  QT11(KL,IJ) = RR
C
                  RR =   (QQIJKL + QQIJLK)
     *                 - (Q11(I,J,K,L) + Q11(I,J,L,K))
     *                 + (R11(I,J,K,L) + R11(I,J,L,K))
                  RS11(KL,IJ) = RR                    
                  IF (I .EQ. J) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  IF (K .EQ. L) RS11(KL,IJ) = FF * RS11(KL,IJ)
                  RR =   (QQIJKL - QQIJLK)
     *                 - (Q11(I,J,K,L) - Q11(I,J,L,K))
     *                 + (R11(I,J,K,L) - R11(I,J,L,K))
                  RT11(KL,IJ) = RR
C 
                  US11(KL,IJ) = US11(KL,IJ) * DP5 * DP25
                  VS11(KL,IJ) = VS11(KL,IJ) * DP5
                  WS11(KL,IJ) = WS11(KL,IJ) * DP5 
                  QS11(KL,IJ) = QS11(KL,IJ) * DP25
                  RS11(KL,IJ) = RS11(KL,IJ) * DP25
C
                  UT11(KL,IJ) = UT11(KL,IJ) * D1P5 * DP25 
                  VT11(KL,IJ) = VT11(KL,IJ) * D1P5
                  WT11(KL,IJ) = WT11(KL,IJ) * D1P5  
                  QT11(KL,IJ) = QT11(KL,IJ) * DP75
                  RT11(KL,IJ) = RT11(KL,IJ) * DP75
C
                  DO 211 MN = 1, NSPAIR
                     BS11(KL,IJ,MN) = BS11(KL,IJ,MN) * DP5 * DP25
                     BT11(KL,IJ,MN) = BT11(KL,IJ,MN) * D1P5 * DP25 
  211             CONTINUE
c
  303          CONTINUE
  302       CONTINUE
  301    CONTINUE
  300 CONTINUE
C
      IF (IPRINT .LE. 3 .OR. IPDD .GT. 1) GOTO 998
      CALL AROUND('Singlet <V@> matrix') 
      CALL OUTPUT(VS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Singlet <B@> matrix') 
      CALL OUTPUT(US11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Singlet <S@> matrix') 
      CALL OUTPUT(QS11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Triplet <V@> matrix') 
      CALL OUTPUT(VT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Triplet <B@> matrix') 
      CALL OUTPUT(UT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      CALL AROUND('Triplet <S@> matrix') 
      CALL OUTPUT(QT11,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
C
  998 IJ = 0
      DO 490 I = 1, NOCDIM
         DO 491 J = 1, I
            IJ = IJ + 1
            EVS(IJ) = EV(I) + EV(J)
  491    CONTINUE
  490 CONTINUE
C
      IF (DOB) THEN
         LU43 = 43
         CALL GPOPEN(LU43,'QMATSP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
         READ(LU43,'(4E30.20)') US11
         CALL GPCLOSE(LU43,'KEEP')
         CALL GPOPEN(LU43,'QMATTP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
         READ(LU43,'(4E30.20)') UT11
         CALL GPCLOSE(LU43,'KEEP')
         DO 701 MN = 1, NSPAIR
            DO 702 J = 1, NSPAIR
               DO 703 I = 1, NSPAIR
                  BS11(I,J,MN) = BS11(I,J,MN) - US11(I,J)
                  BT11(I,J,MN) = BT11(I,J,MN) - UT11(I,J)
  703         CONTINUE
  702      CONTINUE
  701   CONTINUE   
        IF (DOBP) THEN         
          LU43 = 43
          CALL GPOPEN(LU43,'XMATSP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') US11  
          CALL GPCLOSE(LU43,'KEEP')
          CALL GPOPEN(LU43,'XMATTP','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
          READ(LU43,'(4E30.20)') UT11   
          CALL GPCLOSE(LU43,'KEEP') 
          DO 704 MN = 1, NSPAIR
             DO 705 J = 1, NSPAIR
                DO 706 I = 1, NSPAIR
                   BS11(I,J,MN) = BS11(I,J,MN) + US11(I,J)
                   BT11(I,J,MN) = BT11(I,J,MN) + UT11(I,J)
  706           CONTINUE
  705        CONTINUE
  704     CONTINUE   
        END IF
      END IF   
C
      IF (IPDD .EQ. 1) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/A'//
     *                   ' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 2) THEN 
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/A PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 3) THEN 
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/A'' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 4) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/A*'//
     *                   ' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 5) THEN 
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/A* PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 6) THEN 
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/A*'' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 7) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B'//
     *                   ' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 8) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 9) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B'''//
     *                   ' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 10) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B'' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 11) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B*'//
     *                   ' PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 12) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B* PAIR ENERGIES:'
      ELSE IF (IPDD .EQ. 13) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER NONINVARIANT R12/B*'''//
     *                   ' PAIR ENERGIES:'  
      ELSE IF (IPDD .EQ. 14) THEN
        WRITE(LUPRI,'(/A/)') ' SECOND-ORDER R12/B*'' PAIR ENERGIES:'  
      END IF
      CALL DZERO(FS,NSPAIR)
      CALL DZERO(FT,NSPAIR)
      F2S = D0
      F2T = D0
      IJ = 0
      DO 500 I = 1, NOCDIM
         DO 501 J = 1, I
            IJ = IJ + 1
            CNINV(IJ,1) = WS11(IJ,IJ)
            CNINV(IJ,2) = BS11(IJ,IJ,IJ)
            IF (ABS(CNINV(IJ,2)) .GT. VCLTHR) THEN
               CNINV(IJ,3) = CNINV(IJ,1) / CNINV(IJ,2)
            ELSE
               CNINV(IJ,3) = D0
            END IF  
            CNINV(IJ,4) = CNINV(IJ,1) * CNINV(IJ,3)
            IF (I .NE. J) THEN
               CNINV(IJ,5) = WT11(IJ,IJ)
               CNINV(IJ,6) = BT11(IJ,IJ,IJ)
               IF (ABS(CNINV(IJ,6)) .GT. VCLTHR) THEN
                  CNINV(IJ,7) = CNINV(IJ,5) / CNINV(IJ,6)
               ELSE
                  CNINV(IJ,7) = D0
               END IF  
               CNINV(IJ,8) = CNINV(IJ,5) * CNINV(IJ,7) 
            ELSE
               CNINV(IJ,5) = D0
               CNINV(IJ,6) = D0
               CNINV(IJ,7) = D0
               CNINV(IJ,8) = D0
            END IF 
            FS(IJ) = CNINV(IJ,4) 
            FT(IJ) = CNINV(IJ,8) 
            IF (.NOT. DONINV) THEN
            CALL R12MP2(WS11(1,IJ),WT11(1,IJ),RS11,RT11,EVS(IJ),
     *      NOCDIM,NSPAIR,FS(IJ),FT(IJ),VS11,VT11,QS11,QT11,
     *      BS11(1,1,IJ),BT11(1,1,IJ),EVS,XPDD,QQ11,IJ,WSMIN,WTMIN)
            WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                       ET(IJ),ET(IJ)+FT(IJ),
     *                                       WSMIN,WTMIN 
            ELSE
            WRITE(LUPRI,'(I4,4F15.9,2G16.6)') IJ,ES(IJ),ES(IJ)+FS(IJ),
     *                                       ET(IJ),ET(IJ)+FT(IJ)
            ENDIF
            F2S = F2S + FS(IJ)
            F2T = F2T + FT(IJ) 
  501    CONTINUE
  500 CONTINUE 
      WRITE(LUPRI,'(/4X,6F15.9)') E2S,E2S+F2S,E2T,E2T+F2T
      IF (IPDD .EQ. 1) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)') 
     * ' Noninvariant MP2-R12/A   correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 600 I = 1, NOCDIM
            DO 601 J = 1, I       
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  601       CONTINUE
  600    CONTINUE
      ELSE IF (IPDD .EQ. 2) THEN
         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/A   correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 3) THEN
         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/A''  correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 4) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)')
     * ' Noninvariant MP2-R12/A*  correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 602 I = 1, NOCDIM
            DO 603 J = 1, I
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  603       CONTINUE
  602    CONTINUE
      ELSE IF (IPDD .EQ. 5) THEN
         WRITE(LUPRI,'(/A,F15.9 )')
     * '              MP2-R12/A*  correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 6) THEN
         WRITE(LUPRI,'(/A,F15.9 )')
     * '              MP2-R12/A*'' correlation energy =',
     *   E2S+E2T+F2S+F2T  
      ELSE IF (IPDD .EQ. 7) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)') 
     * ' Noninvariant MP2-R12/B   correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 604 I = 1, NOCDIM
            DO 605 J = 1, I       
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  605       CONTINUE
  604    CONTINUE
      ELSE IF (IPDD .EQ. 8) THEN
         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/B   correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 9) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)')
     * ' Noninvariant MP2-R12/B''  correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 610 I = 1, NOCDIM
            DO 611 J = 1, I
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  611       CONTINUE
  610    CONTINUE
      ELSE IF (IPDD .EQ. 10) THEN
         WRITE(LUPRI,'(/A,F15.9 )') 
     * '              MP2-R12/B''  correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 11) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)')
     * ' Noninvariant MP2-R12/B*  correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 606 I = 1, NOCDIM
            DO 607 J = 1, I
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  607       CONTINUE
  606    CONTINUE
      ELSE IF (IPDD .EQ. 12) THEN
         WRITE(LUPRI,'(/A,F15.9 )')
     * '              MP2-R12/B*  correlation energy =',
     *   E2S+E2T+F2S+F2T
      ELSE IF (IPDD .EQ. 13) THEN
         WRITE(LUPRI,'(/A,F15.9//A/)')
     * ' Noninvariant MP2-R12/B*'' correlation energy =',
     *   E2S+E2T+F2S+F2T,
     * '  IJ    V(IJ)       U(IJ)       C(IJ)   '//
     *     '    V(IJ)       U(IJ)       C(IJ)   '
         IJ = 0
         DO 608 I = 1, NOCDIM
            DO 609 J = 1, I
            IJ = IJ + 1
            WRITE(LUPRI,'(I4,6E12.4)') IJ,(CNINV(IJ,K),K=1,3),
     *                                (CNINV(IJ,K),K=5,7)
  609       CONTINUE
  608    CONTINUE
      ELSE IF (IPDD .EQ. 14) THEN
         WRITE(LUPRI,'(/A,F15.9 )')
     * '              MP2-R12/B*'' correlation energy =',
     *   E2S+E2T+F2S+F2T  
      END IF
  999 CONTINUE
      RETURN
      END
C  /* Deck ryr*/
      SUBROUTINE RYR(R2AM,SF,TF,Q11,NOCDIM,NSPAIR) 
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      DIMENSION R2AM(*), SF(*), TF(*),
     *          ISB(8), Q11(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
      INTEGER IDUM
      LOGICAL LDUM
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      CALL DZERO(Q11,NOCDIM**4)
      ISB(1) = 0
      DO 100 ISYM = 2, NSYM
         NBASF = NORB1(ISYM-1) + NORB2(ISYM-1)
         NNBASF = NBASF * (NBASF + 1) / 2
         ISB(ISYM) = ISB(ISYM-1) + NNBASF
  100 CONTINUE
      NBASF = NORB1(NSYM) + NORB2(NSYM)
      NNBASF = ISB(NSYM) + NBASF * (NBASF + 1) / 2 
      LUMULB = 34
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN',' ','FORMATTED',IDUM,LDUM) 
      READ (LUMULB,'(4E30.20)') (SF(I), I = 1, NNBASF)
      CALL GPCLOSE(LUMULB,'KEEP')
      CALL GPOPEN(LUMULB,'AUXOVL','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ (LUMULB,'(4E30.20)') (TF(I), I = 1, NNBASF)
      CALL GPCLOSE(LUMULB,'KEEP')
      DO 200 ISYM = 1, NSYM 
         IJ = ISB(ISYM)
         DO 210 I = 1, NORB1(ISYM) + NORB2(ISYM)
            DO 220 J = 1, I
            IJ = IJ + 1
            IF ((I .LE. NORB1(ISYM) .AND. J .GT. NORB1(ISYM)) .OR.
     *          (J .LE. NORB1(ISYM) .AND. I .GT. NORB1(ISYM))) THEN
               SF(IJ) = - SF(IJ)
               TF(IJ) = - TF(IJ)
            END IF
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
      DO 301 ISYMIA = 1, NSYM
        ISYMJB = ISYMIA
        DO 302 ISYMI = 1, NSYM
          ISYMA = MULD2H(ISYMI,ISYMIA)
          ISYMC = ISYMA
          DO 303 ISYMJ = 1, NSYM
            ISYMB = MULD2H(ISYMJ,ISYMJB)
            ISYMD = ISYMB 
            DO 304 ISYMKC = 1, NSYM 
              ISYMLD = ISYMKC
              ISYMK = MULD2H(ISYMC,ISYMKC)
              ISYML = MULD2H(ISYMD,ISYMLD)
              DO 305 I = 1, NRHF(ISYMI)   
                KOFFI = IRHF(ISYMI) + I 
                DO 306 A = 1, NVIR(ISYMA)
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I-1) + A
                  DO 307 J = 1, NRHF(ISYMJ)
                    KOFFJ = IRHF(ISYMJ) + J 
                    DO 308 B = 1, NVIR(ISYMB)
                      NBJ = IT1AM(ISYMB,ISYMJ)
     *                    + NVIR(ISYMB)*(J-1) + B 
                      NAIBJ = IT2AM(ISYMIA,ISYMJB)
     *                      + INDEX(NAI,NBJ)
                      DO 309 K = 1, NRHF(ISYMK)
                        KOFFK = IRHF(ISYMK) + K
                        DO 310 L = 1, NRHF(ISYML)
                          KOFFL = IRHF(ISYML) + L
C
                          IF (A .LE. NORB1(ISYMA) .AND.
     *                        B .LE. NORB1(ISYMB)) THEN
C
                            DO 311 C = 1, NVIR(ISYMC)
                              NCK = IT1AM(ISYMC,ISYMK)
     *                            + NVIR(ISYMC)*(K-1) + C
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              SAC = TF(ISB(ISYMA)+INDEX(A,C))
                              DO 312 D = 1, NVIR(ISYMD) 
                                NDL = IT1AM(ISYMD,ISYML)
     *                              + NVIR(ISYMD)*(L-1) + D
                                NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                + INDEX(NCK,NDL)
                                SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                XSX = XAC * SBD + SAC * XBD
                                Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                          Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                          R2AM(NAIBJ) * XSX * R2AM(NCKDL)
  312                         CONTINUE
  311                       CONTINUE  
C
                          ELSE IF (A .LE. NORB1(ISYMA) .AND.
     *                             B .GT. NORB1(ISYMB)) THEN
C
                            DO 313 C = 1, NVIR(ISYMC)
                              NCK = IT1AM(ISYMC,ISYMK)
     *                            + NVIR(ISYMC)*(K-1) + C
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              SAC = TF(ISB(ISYMA)+INDEX(A,C))
                              D = B
                              NDL = IT1AM(ISYMD,ISYML)
     *                            + NVIR(ISYMD)*(L-1) + D  
                              NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                              + INDEX(NCK,NDL)
                              SBD = TF(ISB(ISYMB)+INDEX(B,D))
                              XBD = SF(ISB(ISYMB)+INDEX(B,D))
                              XSX = XAC * SBD + SAC * XBD
                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                        Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                        R2AM(NAIBJ) * XSX * R2AM(NCKDL) 
                              DO 314 D = 1, NORB1(ISYMD)
                                NDL = IT1AM(ISYMD,ISYML)
     *                              + NVIR(ISYMD)*(L-1) + D
                                NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                + INDEX(NCK,NDL)
                                SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                XSX = XAC * SBD + SAC * XBD
                                Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                          Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                          R2AM(NAIBJ) * XSX * R2AM(NCKDL)
  314                         CONTINUE
  313                       CONTINUE  
C
                          ELSE IF (A .GT. NORB1(ISYMA) .AND.
     *                             B .LE. NORB1(ISYMB)) THEN
C
                            DO 315 D = 1, NVIR(ISYMD)
                              NDL = IT1AM(ISYMD,ISYML)
     *                            + NVIR(ISYMD)*(L-1) + D
                              SBD = TF(ISB(ISYMB)+INDEX(B,D))
                              XBD = SF(ISB(ISYMB)+INDEX(B,D))
                              C = A
                              NCK = IT1AM(ISYMC,ISYMK)
     *                            + NVIR(ISYMC)*(K-1) + C
                              NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                              + INDEX(NCK,NDL)
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              SAC = TF(ISB(ISYMA)+INDEX(A,C))
                              XSX = XAC * SBD + SAC * XBD
                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                        Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                        R2AM(NAIBJ) * XSX * R2AM(NCKDL)
                              DO 316 C = 1, NORB1(ISYMC)
                                NCK = IT1AM(ISYMC,ISYMK)
     *                              + NVIR(ISYMC)*(K-1) + C    
                                NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                + INDEX(NCK,NDL)
                                XAC = SF(ISB(ISYMA)+INDEX(A,C))
                                SAC = TF(ISB(ISYMA)+INDEX(A,C))
                                XSX = XAC * SBD + SAC * XBD
                                Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                          Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                          R2AM(NAIBJ) * XSX * R2AM(NCKDL)
  316                         CONTINUE
  315                       CONTINUE
C
                          ELSE IF (A .GT. NORB1(ISYMA) .AND.
     *                             B .GT. NORB1(ISYMB)) THEN
C
                            DO 317 C = 1, NORB1(ISYMC)
                              NCK = IT1AM(ISYMC,ISYMK)
     *                            + NVIR(ISYMC)*(K-1) + C
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              SAC = TF(ISB(ISYMA)+INDEX(A,C))
                              DO 318 D = 1, NORB1(ISYMD)
                                NDL = IT1AM(ISYMD,ISYML)
     *                              + NVIR(ISYMD)*(L-1) + D
                                NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                + INDEX(NCK,NDL)
                                SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                XSX = XAC * SBD + SAC * XBD
                                Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                          Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                          R2AM(NAIBJ) * XSX * R2AM(NCKDL)
  318                         CONTINUE
  317                       CONTINUE
                            C = A
                            NCK = IT1AM(ISYMC,ISYMK)
     *                          + NVIR(ISYMC)*(K-1) + C  
                            XAC = SF(ISB(ISYMA)+INDEX(A,C))
                            SAC = TF(ISB(ISYMA)+INDEX(A,C)) 
                            DO 319 D = 1, NORB1(ISYMD)
                              NDL = IT1AM(ISYMD,ISYML)
     *                            + NVIR(ISYMD)*(L-1) + D
                              NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                              + INDEX(NCK,NDL)
                              SBD = TF(ISB(ISYMB)+INDEX(B,D))
                              XBD = SF(ISB(ISYMB)+INDEX(B,D))
                              XSX = XAC * SBD + SAC * XBD
                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                        Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                        R2AM(NAIBJ) * XSX * R2AM(NCKDL)
  319                       CONTINUE 
                            D = B
                            NDL = IT1AM(ISYMD,ISYML)
     *                          + NVIR(ISYMD)*(L-1) + D   
                            SBD = TF(ISB(ISYMB)+INDEX(B,D))
                            XBD = SF(ISB(ISYMB)+INDEX(B,D))
                            DO 320 C = 1, NORB1(ISYMC)
                              NCK = IT1AM(ISYMC,ISYMK)
     *                            + NVIR(ISYMC)*(K-1) + C
                              NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                              + INDEX(NCK,NDL)
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              SAC = TF(ISB(ISYMA)+INDEX(A,C))
                              XSX = XAC * SBD + SAC * XBD
                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                        Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                        R2AM(NAIBJ) * XSX * R2AM(NCKDL)
  320                       CONTINUE
                            C = A
                            D = B
                            NCK = IT1AM(ISYMC,ISYMK)
     *                          + NVIR(ISYMC)*(K-1) + C
                            NDL = IT1AM(ISYMD,ISYML)
     *                          + NVIR(ISYMD)*(L-1) + D
                            NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                            + INDEX(NCK,NDL)
                            XAC = SF(ISB(ISYMA)+INDEX(A,C))
                            SAC = TF(ISB(ISYMA)+INDEX(A,C))
                            SBD = TF(ISB(ISYMB)+INDEX(B,D))
                            XBD = SF(ISB(ISYMB)+INDEX(B,D))
                            XSX = XAC * SBD + SAC * XBD
                            Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                      Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                      R2AM(NAIBJ) * XSX * R2AM(NCKDL)
C
                          END IF
C
  310                   CONTINUE
  309                 CONTINUE
  308               CONTINUE
  307             CONTINUE
  306           CONTINUE
  305         CONTINUE
  304       CONTINUE
  303     CONTINUE
  302   CONTINUE
  301 CONTINUE
      IJKL = 0
      FF = D1 / SQRT(D2)
      DO 600 K = 1, NOCDIM
         DO 601 L = 1, K
            DO 602 I = 1, NOCDIM
               DO 603 J = 1, I
                  IJKL = IJKL + 1
                  SF(IJKL) = Q11(I,J,K,L) + Q11(I,J,L,K) 
                  TF(IJKL) = Q11(I,J,K,L) - Q11(I,J,L,K) 
                  IF (I .EQ. J) SF(IJKL) = FF * SF(IJKL)
                  IF (K .EQ. L) SF(IJKL) = FF * SF(IJKL)
                  SF(IJKL) = SF(IJKL) * DP25
                  TF(IJKL) = TF(IJKL) * DP75
  603          CONTINUE
  602       CONTINUE
  601    CONTINUE
  600 CONTINUE
      IF (IPRINT .GT.  3) THEN
         CALL AROUND('Singlet <RXR> matrix') 
         CALL OUTPUT(SF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <RXR> matrix') 
         CALL OUTPUT(TF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck rzr*/
      SUBROUTINE RZR(R2AM,V2AM,SF,TF,Q11,NOCDIM,NSPAIR) 
C
C     This subroutine is not used!
C
C     It evaluates the sum
C
C     Q(ijkl) = Sum_abcd (ia|r12|jb) * [X(ac)*S(bd) + S(ac)*X(db)] * (kc|r12|ld)
C
C     The sum over a,b,c,d runs over both the primary and the secondary basis.
C
C     On input, the array R2AM contains the integrals (ia|r12|jb).
C     NOCDIM is the number of (nonfrozen) occupied orbitals and NSPAIR
C     is the number of pairs of (nonfrozen) occupied orbitals.
C
C     On output, the arrays SF and TF contain the matrices Q(ijkl) for
C     singlet and triplet coupled pairs, respectively.
C
C     V2AM (of length NT2AMX) and Q11 (of length NOCDIM**4) are used
C     for scratch space.
C
C     The one-particle matrices X (AUXFCK) and S (AUXOVL) are read from disk.
C     These are the primary+secondary exchange and overlap matrices in the
C     the orthonormal bases, respectively. It is assumed that these matrices
C     are diagonal in the secondary-secondary block. This is tested for.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0, 
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      PARAMETER (THRDIA = 1D-9)
      DIMENSION R2AM(*), V2AM(*), SF(*), TF(*),
     *          ISB(8), Q11(NOCDIM,NOCDIM,NOCDIM,NOCDIM)
      LOGICAL LDUM
      INTEGER IDUM
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      CALL DZERO(Q11,NOCDIM**4)
      ISB(1) = 0
      DO 100 ISYM = 2, NSYM
         NBASF = NORB1(ISYM-1) + NORB2(ISYM-1)
         NNBASF = NBASF * (NBASF + 1) / 2
         ISB(ISYM) = ISB(ISYM-1) + NNBASF
  100 CONTINUE
      NBASF = NORB1(NSYM) + NORB2(NSYM)
      NNBASF = ISB(NSYM) + NBASF * (NBASF + 1) / 2 
      LUMULB = 34
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ (LUMULB,'(4E30.20)') (SF(I), I = 1, NNBASF)
      CALL GPCLOSE(LUMULB,'KEEP')
      CALL GPOPEN(LUMULB,'AUXOVL','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ (LUMULB,'(4E30.20)') (TF(I), I = 1, NNBASF)
      CALL GPCLOSE(LUMULB,'KEEP')
      DO 200 ISYM = 1, NSYM 
         IJ = ISB(ISYM)
         DO 210 I = 1, NORB1(ISYM) + NORB2(ISYM)
            DO 220 J = 1, I
            IJ = IJ + 1
            IF ((I .LE. NORB1(ISYM) .AND. J .GT. NORB1(ISYM)) .OR.
     *          (J .LE. NORB1(ISYM) .AND. I .GT. NORB1(ISYM))) THEN
               SF(IJ) = - SF(IJ)
               TF(IJ) = - TF(IJ)
            ELSE IF ((I .GT. NORB1(ISYM) .AND. J .GT. NORB1(ISYM))) THEN
               IF (I .NE. J) THEN
                  IF (ABS(SF(IJ)) .GT. THRDIA) THEN
                     WRITE(LUPRI,'(/A/A,3I5,E20.10/)') 
     *               '@ WARNING : Exchange matrix not diagonal',
     *               '            Nondiagonal element is :',
     *                            ISYM,I,J,SF(IJ)
                     IF (ABS(SF(IJ)) .GT. SQRT(THRDIA))
     *                  CALL QUIT('Exchange matrix not diagonal')
                  END IF
                  IF (ABS(TF(IJ)) .GT. THRDIA) THEN
                     WRITE(LUPRI,'(/A/A,3I5,E20.10/)') 
     *               '@ WARNING : Overlap matrix not diagonal',
     *               '            Nondiagonal element is :',
     *                            ISYM,I,J,TF(IJ)
                     IF (ABS(TF(IJ)) .GT. SQRT(THRDIA))
     *                  CALL QUIT('Overlap matrix not diagonal')
                  END IF
               END IF
            END IF
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      DO 301 ISYMKA = 1, NSYM
         ISYMLB = ISYMKA
         DO 302 ISYMK = 1, NSYM
            ISYMA = MULD2H(ISYMK,ISYMKA)
            ISYMC = ISYMA
            ISYMKC = ISYMKA
            DO 303 ISYML = 1, NSYM
               ISYMB = MULD2H(ISYML,ISYMLB)
               ISYMD = ISYMB 
               ISYMLD = ISYMLB
               DO 304 K = 1, NRHF(ISYMK)   
                  DO 305 L = 1, NRHF(ISYML) 
                     DO 306 A = 1, NVIR(ISYMA)
                        NAK = IT1AM(ISYMA,ISYMK)
     *                      + NVIR(ISYMA)*(K-1) + A
                        DO 307 B = 1, NVIR(ISYMB)
                           NBL = IT1AM(ISYMB,ISYML)
     *                         + NVIR(ISYMB)*(L-1) + B 
                           NAKBL = IT2AM(ISYMKA,ISYMLB)
     *                           + INDEX(NAK,NBL)
C
                           V2AM(NAKBL) = D0
C
                           IF (A .LE. NORB1(ISYMA) .AND.
     *                         B .LE. NORB1(ISYMB)) THEN
C
                              DO 308 C = 1, NVIR(ISYMC)
                                 NCK = IT1AM(ISYMC,ISYMK)
     *                               + NVIR(ISYMC)*(K-1) + C
                                 XAC = SF(ISB(ISYMA)+INDEX(A,C))
                                 SAC = TF(ISB(ISYMA)+INDEX(A,C))
                                 DO 309 D = 1, NVIR(ISYMD) 
                                    NDL = IT1AM(ISYMD,ISYML)
     *                                  + NVIR(ISYMD)*(L-1) + D
                                    NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                    + INDEX(NCK,NDL)
                                    SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                    XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                    XSX = XAC * SBD + SAC * XBD
                                    V2AM(NAKBL) = V2AM(NAKBL)
     *                                          + XSX * R2AM(NCKDL)
  309                            CONTINUE
  308                         CONTINUE  
C
                           ELSE IF (A .LE. NORB1(ISYMA) .AND.
     *                              B .GT. NORB1(ISYMB)) THEN
C
                              DO 310 C = 1, NVIR(ISYMC)
                                 NCK = IT1AM(ISYMC,ISYMK)
     *                               + NVIR(ISYMC)*(K-1) + C
                                 XAC = SF(ISB(ISYMA)+INDEX(A,C))
                                 SAC = TF(ISB(ISYMA)+INDEX(A,C))
                                 D = B
                                 NDL = NBL               
                                 NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                 + INDEX(NCK,NDL)
                                 SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                 XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                 XSX = XAC * SBD + SAC * XBD
                                 V2AM(NAKBL) = V2AM(NAKBL)
     *                                       + XSX * R2AM(NCKDL)
                                 DO 311 D = 1, NORB1(ISYMD)
                                    NDL = IT1AM(ISYMD,ISYML)
     *                                  + NVIR(ISYMD)*(L-1) + D
                                    NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                    + INDEX(NCK,NDL)
                                    SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                    XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                    XSX = XAC * SBD + SAC * XBD
                                    V2AM(NAKBL) = V2AM(NAKBL)
     *                                          + XSX * R2AM(NCKDL)
  311                            CONTINUE
  310                         CONTINUE  
C
                           ELSE IF (A .GT. NORB1(ISYMA) .AND.
     *                              B .LE. NORB1(ISYMB)) THEN
C
                              DO 312 D = 1, NVIR(ISYMD)
                                 NDL = IT1AM(ISYMD,ISYML)
     *                               + NVIR(ISYMD)*(L-1) + D
                                 SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                 XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                 C = A
                                 NCK = NAK               
                                 NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                 + INDEX(NCK,NDL)
                                 XAC = SF(ISB(ISYMA)+INDEX(A,C))
                                 SAC = TF(ISB(ISYMA)+INDEX(A,C))
                                 XSX = XAC * SBD + SAC * XBD
                                 V2AM(NAKBL) = V2AM(NAKBL)
     *                                       + XSX * R2AM(NCKDL)
                                 DO 313 C = 1, NORB1(ISYMC)
                                    NCK = IT1AM(ISYMC,ISYMK)
     *                                  + NVIR(ISYMC)*(K-1) + C    
                                    NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                    + INDEX(NCK,NDL)
                                    XAC = SF(ISB(ISYMA)+INDEX(A,C))
                                    SAC = TF(ISB(ISYMA)+INDEX(A,C))
                                    XSX = XAC * SBD + SAC * XBD
                                    V2AM(NAKBL) = V2AM(NAKBL)
     *                                          + XSX * R2AM(NCKDL)
  313                            CONTINUE
  312                         CONTINUE
C
                           ELSE IF (A .GT. NORB1(ISYMA) .AND.
     *                              B .GT. NORB1(ISYMB)) THEN
C
                              DO 314 C = 1, NORB1(ISYMC)
                                 NCK = IT1AM(ISYMC,ISYMK)
     *                               + NVIR(ISYMC)*(K-1) + C
                                 XAC = SF(ISB(ISYMA)+INDEX(A,C))
                                 SAC = TF(ISB(ISYMA)+INDEX(A,C))
                                 DO 315 D = 1, NORB1(ISYMD)
                                    NDL = IT1AM(ISYMD,ISYML)
     *                                  + NVIR(ISYMD)*(L-1) + D
                                    NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                    + INDEX(NCK,NDL)
                                    SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                    XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                    XSX = XAC * SBD + SAC * XBD
                                    V2AM(NAKBL) = V2AM(NAKBL)
     *                                          + XSX * R2AM(NCKDL)
  315                            CONTINUE
  314                         CONTINUE
                              C = A
                              NCK = NAK               
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              SAC = TF(ISB(ISYMA)+INDEX(A,C)) 
                              DO 316 D = 1, NORB1(ISYMD)
                                 NDL = IT1AM(ISYMD,ISYML)
     *                               + NVIR(ISYMD)*(L-1) + D
                                 NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                 + INDEX(NCK,NDL)
                                 SBD = TF(ISB(ISYMB)+INDEX(B,D))
                                 XBD = SF(ISB(ISYMB)+INDEX(B,D))
                                 XSX = XAC * SBD + SAC * XBD
                                 V2AM(NAKBL) = V2AM(NAKBL)
     *                                       + XSX * R2AM(NCKDL)
  316                         CONTINUE 
                              D = B
                              NDL = NBL 
                              SBD = TF(ISB(ISYMB)+INDEX(B,D))
                              XBD = SF(ISB(ISYMB)+INDEX(B,D))
                              DO 317 C = 1, NORB1(ISYMC)
                                 NCK = IT1AM(ISYMC,ISYMK)
     *                               + NVIR(ISYMC)*(K-1) + C
                                 NCKDL = IT2AM(ISYMKC,ISYMLD)
     *                                 + INDEX(NCK,NDL)
                                 XAC = SF(ISB(ISYMA)+INDEX(A,C))
                                 SAC = TF(ISB(ISYMA)+INDEX(A,C))
                                 XSX = XAC * SBD + SAC * XBD
                                 V2AM(NAKBL) = V2AM(NAKBL)
     *                                       + XSX * R2AM(NCKDL)
  317                         CONTINUE
                              C = A
                              D = B
                              NCK = NAK
                              NDL = NBL
                              NCKDL = NAKBL                     
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              SAC = TF(ISB(ISYMA)+INDEX(A,C))
                              SBD = TF(ISB(ISYMB)+INDEX(B,D))
                              XBD = SF(ISB(ISYMB)+INDEX(B,D))
                              XSX = XAC * SBD + SAC * XBD
                              V2AM(NAKBL) = V2AM(NAKBL)
     *                                    + XSX * R2AM(NCKDL)
C
                           END IF
C
  307                   CONTINUE
  306                CONTINUE
  305             CONTINUE
  304          CONTINUE
  303       CONTINUE
  302    CONTINUE
  301 CONTINUE
C
      DO 401 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         DO 402 ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMIA)
            DO 403 ISYMJ = 1, NSYM
               ISYMB = MULD2H(ISYMJ,ISYMJB)
               DO 404 I = 1, NRHF(ISYMI)
                  KOFFI = IRHF(ISYMI) + I
                  DO 405 J = 1, NRHF(ISYMJ)
                     KOFFJ = IRHF(ISYMJ) + J
                     DO 406 A = 1, NVIR(ISYMA)                 
                        ISYMJA = MULD2H(ISYMJ,ISYMA)
                        DO 407 B = 1, NVIR(ISYMB)
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I-1) + A
                           NBJ = IT1AM(ISYMB,ISYMJ)
     *                         + NVIR(ISYMB)*(J-1) + B
                           NAIBJ = IT2AM(ISYMIA,ISYMJB)
     *                           + INDEX(NAI,NBJ)
                           DO 408 ISYMKA = 1, NSYM
                              ISYMLB = ISYMKA
                              ISYMK = MULD2H(ISYMA,ISYMKA)
                              ISYML = MULD2H(ISYMB,ISYMLB)
                              DO 409 K = 1, NRHF(ISYMK)
                                 KOFFK = IRHF(ISYMK) + K
                                 DO 410 L = 1, NRHF(ISYML)
                                    KOFFL = IRHF(ISYML) + L  
                                    NAK = IT1AM(ISYMA,ISYMK)
     *                                  + NVIR(ISYMA)*(K-1) + A
                                    NBL = IT1AM(ISYMB,ISYML)
     *                                  + NVIR(ISYMB)*(L-1) + B
                                    NAKBL = IT2AM(ISYMKA,ISYMLB)
     *                                    + INDEX(NAK,NBL)
                                    Q11(KOFFI,KOFFJ,KOFFK,KOFFL) =
     *                              Q11(KOFFI,KOFFJ,KOFFK,KOFFL) +
     *                              R2AM(NAIBJ) * V2AM(NAKBL)
  410                            CONTINUE
  409                         CONTINUE
  408                      CONTINUE
  407                   CONTINUE
  406                CONTINUE
  405             CONTINUE
  404          CONTINUE   
  403       CONTINUE
  402    CONTINUE
  401 CONTINUE   
      IJKL = 0
      FF = D1 / SQRT(D2)
      DO 600 K = 1, NOCDIM
         DO 601 L = 1, K
            DO 602 I = 1, NOCDIM
               DO 603 J = 1, I
                  IJKL = IJKL + 1
                  SF(IJKL) = Q11(I,J,K,L) + Q11(I,J,L,K) 
                  TF(IJKL) = Q11(I,J,K,L) - Q11(I,J,L,K) 
                  IF (I .EQ. J) SF(IJKL) = FF * SF(IJKL)
                  IF (K .EQ. L) SF(IJKL) = FF * SF(IJKL)
                  SF(IJKL) = SF(IJKL) * DP25
                  TF(IJKL) = TF(IJKL) * DP75
  603          CONTINUE
  602       CONTINUE
  601    CONTINUE
  600 CONTINUE
C
      IF (IPRINT .GT.  3) THEN
         CALL AROUND('Singlet <RXR> matrix') 
         CALL OUTPUT(SF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         CALL AROUND('Triplet <RXR> matrix') 
         CALL OUTPUT(TF,1,NSPAIR,1,NSPAIR,NSPAIR,NSPAIR,1,LUPRI)
         WRITE(LUPRI,*)
      END IF
      RETURN
      END
C  /* Deck comkr*/
      SUBROUTINE COMKR(R2AM,V2AM,S2AM,T2AM,SF,LSF)                       
C
C     This subroutine computes the integrals <pq|[r12,K1+K2]|ij>,
C     utilizing the secondary basis RI approximation.
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
#include <priunit.h>  
      PARAMETER (D0 = 0D0, D1 = 1D0, D2 = 2D0, D3 = 3D0, D1P5 = 1.5D0,
     *           DP5 = 0.5D0, DP25 = 0.25D0, DP75 = 0.75D0)
      PARAMETER (THRDIA = 1D-9)
      DIMENSION R2AM(*), T2AM(*), V2AM(*), S2AM(*), SF(*), ISB(8)
      REAL*8  RR, VV, XBD, XAC, DSM
      LOGICAL LDUM
      INTEGER IDUM
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <r12int.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      ISB(1) = 0
      DO 100 ISYM = 2, NSYM
         NBASF = NORB1(ISYM-1) + NORB2(ISYM-1)
         NNBASF = NBASF * (NBASF + 1) / 2
         ISB(ISYM) = ISB(ISYM-1) + NNBASF
  100 CONTINUE
      NBASF = NORB1(NSYM) + NORB2(NSYM)
      NNBASF = ISB(NSYM) + NBASF * (NBASF + 1) / 2
      IF (NNBASF .GT. LSF) 
     * CALL QUIT('NOT ENOUGH SPACE IN COMKR')
      LUMULB = 34
      CALL GPOPEN(LUMULB,'AUXFCK','UNKNOWN',' ','FORMATTED',IDUM,LDUM)
      READ (LUMULB,'(4E30.20)') (SF(I), I = 1, NNBASF)
      CALL GPCLOSE(LUMULB,'KEEP')
      DO 200 ISYM = 1, NSYM
         DO 210 I = NORB1(ISYM) + 2,  NVIR(ISYM)
            DO 220 J = NORB1(ISYM) + 1 , I - 1
               IJ = ISB(ISYM) + INDEX(I,J)
               IF (ABS(SF(IJ)) .GT. THRDIA) THEN
                  WRITE(LUPRI,'(/A/A,3I5,E20.10/)')
     *            '@ WARNING : Exchange matrix not diagonal', 
     *            '            Nondiagonal element is :',
     *                         ISYM,I,J,SF(IJ)
                  IF (ABS(SF(IJ)) .GT. SQRT(THRDIA))
     *               CALL QUIT('Exchange matrix not diagonal')
               END IF
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      DO 401 ISYMIA = 1, NSYM
         ISYMJB = ISYMIA
         NAIBJ = IT2AM(ISYMIA,ISYMJB)
         KOFFU = IU2AM(ISYMIA,ISYMJB)
         NTOTAI = NT1AM(ISYMIA)
         DO 402 IA = 1, NTOTAI
            DO 403 JB = 1, IA
               NAIBJ = NAIBJ + 1
               T2AM(NAIBJ) = V2AM(KOFFU + (IA-1)*NTOTAI + JB) +
     &                       V2AM(KOFFU + (JB-1)*NTOTAI + IA)
  403        CONTINUE
  402    CONTINUE
  401 CONTINUE
C  
      DO 301 ISYMKA = 1, NSYM
         ISYMLB = ISYMKA
         DO 302 ISYMK = 1, NSYM
            ISYMA = MULD2H(ISYMK,ISYMKA)
            ISYMC = ISYMA
            ISYMKC = ISYMKA
            DO 303 ISYML = 1, NSYM
               ISYMB = MULD2H(ISYML,ISYMLB)
               ISYMD = ISYMB
               ISYMLD = ISYMLB
               DO 304 K = 1, NRHF(ISYMK)
                  DO 305 L = 1, NRHF(ISYML)
                     DO 306 A = 1,  NORB1(ISYMA)
                        NSTRC = NORB1(ISYMC) + 1
                        NENDC = NVIR(ISYMC)
                        NAK = IT1AM(ISYMA,ISYMK)
     *                      + NVIR(ISYMA)*(K-1) + A
                        DO 307 B = 1, NORB1(ISYMB) 
                           NSTRD = NORB1(ISYMD) + 1
                           NENDD = NVIR(ISYMD)
                           NBL = IT1AM(ISYMB,ISYML)
     *                         + NVIR(ISYMB)*(L-1) + B
C
                           IF (NBL .GT. NAK) GOTO 307
C
                           NAKBL = IT2AM(ISYMKA,ISYMLB)
     *                           + INDEX(NAK,NBL)
                           DSM = D0
                           DO 308 C = NSTRC, NENDC       
                              NCK = IT1AM(ISYMC,ISYMK)
     *                            + NVIR(ISYMC)*(K-1) + C 
                              NCKBL = IT2AM(ISYMKC,ISYMLB)
     *                              + INDEX(NCK,NBL) 
                              RR = R2AM(NCKBL)
                              XAC = SF(ISB(ISYMA)+INDEX(A,C))
                              DSM = DSM + XAC * RR
  308                      CONTINUE
                           S2AM(NAKBL) = T2AM(NAKBL) - DSM
                           DSM = D0     
                           DO 309 D = NSTRD, NENDD
                              NDL = IT1AM(ISYMD,ISYML)
     *                            + NVIR(ISYMD)*(L-1) + D
                              NAKDL = IT2AM(ISYMKA,ISYMLD)
     *                              + INDEX(NAK,NDL)
                              RR = R2AM(NAKDL)
                              XBD = SF(ISB(ISYMB)+INDEX(B,D))
                              DSM = DSM + XBD * RR
  309                      CONTINUE
                           S2AM(NAKBL) = S2AM(NAKBL) - DSM         
C
  307                   CONTINUE
  306                CONTINUE
  305             CONTINUE
  304          CONTINUE
  303       CONTINUE
  302    CONTINUE
  301 CONTINUE
C   
      RETURN
      END 
