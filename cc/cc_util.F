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
C  /* Deck ccrhs_cio */
      SUBROUTINE CCRHS_CIO2(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                      ISYVEC,ISYCIM,LUC,CFIL,IV,IOPT,
     *                      IOPTB,LUCB,CBFIL,IOFFB,FACB,
     *                      IOPTE,EMAT1)
C
C
C      vectors        ISYVEC
C      intermediates  ISYCIM
C
C      LUC and CFIL is used to control from which file the
C                   intermediate is obtained.
C
C      if iopt = 1 the C intermediate is assumed
C         to be as in energy calc.
C
C      if iopt ne. 1 we use the intermediate
C         on luc with address given according to
C         transformed vector nr iv (iv is not 1
C         if several vectors are transformed
C         simultaneously.)
C
C      if iopte = 1 add contribution of E1 intermediate to diagonal
C
C      extra contribution in B matrix calculation:
C       if IOPTB=1 contribution form CBAR intermediate is included. 
C       the intermediate is read from file CBFIL (unit LUCB) with
C       the offset IOFFB
C
C      in energy calc: iv=1,iopt=1
C
C
C     PURPOSE:
C             Calculate the C-term 
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*),EMAT1(*)
      DIMENSION WORK(LWORK)
      CHARACTER CFIL*(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      XNORM = 0.0d0
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*) 'I/O in C-term not implemented for '//
     &        'square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_CIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYCIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYCIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
C
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for allocation in '//
     &           'CCRHS_CIO-1')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KSCR3 = KSCR2 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR3 + NT1AM(ISYMCK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for allocation '//
     &              'in CCRHS_CIO-2')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUC,CFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C-----------------------------------------
C              Transform delta index to c.
C-----------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE

C           XNORM = XNORM + 
C    &         DDOT(NDISAI*NT1AM(ISYMDK),WORK(KSCR1),1,WORK(KSCR1),1)
C           WRITE (LUPRI,*) 'CCRHS_CIO2> XNORM,YNORM:',XNORM,YNORM

C
C------------------------------------------
C         add E1 to contrib. to R diagonal:
C------------------------------------------
C
          IF (IOPTE .EQ. 1) THEN

             DO ISYMC = 1, NSYM

               ISYMK  = MULD2H(ISYMC,ISYMCK)
               ISYMA  = MULD2H(ISYMK,ISYMAI)

               DO A = 1, NVIR(ISYMA)
               DO K = 1, NRHF(ISYMK)

                 NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K-1) + A
                 KAK = NAK - IFSTAI + 1

                 IF (NAK.GE.IFSTAI .AND. NAK.LE.ILSTAI) THEN

                   DO C = 1, NVIR(ISYMC)

                      NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) + C

                      KOFF1 = IMATAB(ISYMA,ISYMC)+NVIR(ISYMA)*(C-1)+A

                      KOFF6 = KSCR1-1 + NT1AM(ISYMCK)*(KAK - 1) + NCK

                      WORK(KOFF6) = WORK(KOFF6) - 0.5D0 * EMAT1(KOFF1) 

                   END DO

                 END IF

               END DO
               END DO

             END DO

          END IF
C
C---------------------------------------------
C         add FACB * CBAR to P(ck,#ai)
C---------------------------------------------
C
          IF (IOPTB .EQ. 1) THEN
            DO NAI = IFSTAI, ILSTAI

              IADR = IOFFB+IT2SQ(ISYMCK,ISYMAI)+NT1AM(ISYMCK)*(NAI-1)+1

              CALL GETWA2(LUCB,CBFIL,WORK(KSCR3),IADR,NT1AM(ISYMCK))

              KAI   = NAI - IFSTAI + 1
              KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI-1)

C             CALL DAXPY(NT1AM(ISYMCK),ONE,WORK(KSCR3),1,WORK(KOFF6),1)
             CALL DAXPY(NT1AM(ISYMCK),FACB,WORK(KSCR3),1,WORK(KOFF6),1)

              XNORM = XNORM + DDOT(NT1AM(ISYMCK),WORK(KSCR3),1,
     &                                           WORK(KSCR3),1)

            END DO
          END IF

C           WRITE (LUPRI,*) 'CCRHS_CIO2> norm of CBAR:',XNORM

C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 190 I = 1,NRHF(ISYMI)
C
                  DO 200 A = 1,NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     IF ((NAI .LT. IFSTAI) .OR. (NAI .GT. ILSTAI))
     *                  GOTO 200
C
                     DO 210 ISYMJ = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)
C
                        DO 220 J = 1,NRHF(ISYMJ)
C
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J-1) + A
C
                           CALL CC_PUTC(WORK(KSCR2),OMEGA2,ISYMAI,
     *                                  ISYMAJ,ISYMBI,ISYMBJ,ISYMB,
     *                                  ISYMI,ISYMJ,NAI,NAJ,I,J,
     *                                  IFSTAI)
C
  220                   CONTINUE
  210                CONTINUE
  200             CONTINUE
  190          CONTINUE
  180       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_dio */
      SUBROUTINE CCRHS_DIO2(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                      ISYVEC,ISYDIM,LUD,DFIL,LUC,CFIL,IV,IOPT,
     *                      IOPTB,LUDB,DBFIL,IOFFB,FACB,
     *                      IOPTE,EMAT1)
C
C     asm 20-aug-1994
C
C      vectors        ISYVEC
C      intermediates  ISYDIM 
C      LUD and DFIL are used to control from which file the
C      intermediate is obtained.
C
C      if iopt = 1 the D intermediate is assumed
C         to be as in energy calc.
C
C      if iopt ne. 1 we use the intermediate
C         on luc with address given according to
C         transformed vector nr iv (iv is not 1
C         if several vectors are transformed
C         simultaneously.)
C
C      it iopte = 1 add contribution of E1 intermediate to diagonal
C
C      extra contribution in B matrix calculation:
C       if IOPTB=1 contribution from DBAR 
C       intermediate in included. the intermediate
C       is read from file DBFIL (unit LUDB) with
C       the offset IOFFB
C
C      in energy calc: iv=1,iopt=1
C
C     PURPOSE:
C             Calculate the D-term 
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*),EMAT1(*)
      DIMENSION WORK(LWORK)
      CHARACTER DFIL*(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      XNORM = 0.0d0
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*) 'I/O in D-term not implemented for '//
     &        'square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_DIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYDIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYDIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI 
         IF (LENMIN .EQ. 0) GOTO 100
         LENMIN = LENMIN + NRHFT
C
         NDISAI = ( LWORK - NT1AM(ISYMCK) ) / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_DIO2')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KSCR3 = KSCR2 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR3 + NT1AM(ISYMCK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for allocation '//
     &              'in CCRHS_DIO')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  KSCR4 = KEND
                  IF (LWRK1 .LT. LEN) THEN
                     CALL QUIT('Insufficient memory in CCRHS_DIO2')
                  END IF
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUD,DFIL,WORK(KOFF1),IOFF,LEN)
                     CALL GETWA2(LUC,CFIL,WORK(KSCR4),IOFF,LEN)
C
                     CALL DSCAL(LEN,TWO,WORK(KOFF1),1)
                     CALL DAXPY(LEN,-ONE,WORK(KSCR4),1,WORK(KOFF1),1)
C
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C--------------------------------------
C           Transform delta index to c.
C--------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE
C
C--------------------------------------------------
C           contract to contrib. to R intermediate:
C--------------------------------------------------
C
          IF (IOPTE .EQ. 1) THEN

             DO ISYMC = 1, NSYM

               ISYMK  = MULD2H(ISYMC,ISYMCK)
               ISYMA  = MULD2H(ISYMK,ISYMAI)

               DO A = 1, NVIR(ISYMA)
               DO K = 1, NRHF(ISYMK)

                 NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K-1) + A
                 KAK = NAK - IFSTAI + 1

                 IF (NAK.GE.IFSTAI .AND. NAK.LE.ILSTAI) THEN

                   DO C = 1, NVIR(ISYMC)

                      NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) + C

                      KOFF1 = IMATAB(ISYMA,ISYMC)+NVIR(ISYMA)*(C-1)+A

                      KOFF6 = KSCR1-1 + NT1AM(ISYMCK)*(KAK - 1) + NCK

                      WORK(KOFF6) = WORK(KOFF6) + 0.5D0 * EMAT1(KOFF1)

                   END DO

                 END IF

               END DO
               END DO

             END DO

          END IF
C
C---------------------------------------------
C         add DBAR intermediate to P(ck,#ai)
C---------------------------------------------
C
          IF (IOPTB .EQ. 1) THEN
            DO NAI = IFSTAI, ILSTAI

              IADR = IOFFB+IT2SQ(ISYMCK,ISYMAI)+NT1AM(ISYMCK)*(NAI-1)+1

              CALL GETWA2(LUDB,DBFIL,WORK(KSCR3),IADR,NT1AM(ISYMCK))

              KAI   = NAI - IFSTAI + 1
              KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI-1)

C             CALL DAXPY(NT1AM(ISYMCK),ONE,WORK(KSCR3),1,WORK(KOFF6),1)
             CALL DAXPY(NT1AM(ISYMCK),FACB,WORK(KSCR3),1,WORK(KOFF6),1)

              XNORM = XNORM + DDOT(NT1AM(ISYMCK),WORK(KSCR3),1,
     &                                           WORK(KSCR3),1)

            END DO
          END IF

C           WRITE (LUPRI,*) 'CCRHS_CIO2> norm of DBAR:',XNORM
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 NAI = IFSTAI,ILSTAI
C
               CALL CC_PUTD(WORK(KSCR2),OMEGA2,ISYMAI,ISYMBJ,NAI,IFSTAI)
C
  180       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
*=====================================================================*
      SUBROUTINE CC_EIM(EMAT1,EMAT2,
     *                  RAIM,RBIM,GAIM,GBIM,FOCK,ONEHAM,
     *                  XLAMDHA,XLAMDPA,ISYLMA,XLAMDHB,XLAMDPB,ISYLMB,
     *                  FCKCON,RTRAN,GTRAN,LRELAX,IOPT,ISYMEI)
*---------------------------------------------------------------------*
*
*     Transforms delta index of R intermediates to virtual and
*     and the alpha index of the G intermediates to occupied and
*     calculate the E intermediates in the MO basis
*
*     IOPT = 1:     EMAT1 = FOCK   - XLAMDHA * RAIM
*                   EMAT2 = ONEHAM + XLAMDPA * GAIM
*                   
*                   RBIM,GBIM,XLAMDHB,XLAMDPB,ISYMB is dummy input
*
*
*     IOPT = 2:     EMAT1 = FOCK   - XLAMDHA * RBIM + XLAMDHB * RAIM
*                   EMAT2 = ONEHAM + XLAMDPA * GBIM + XLAMDPB * GAIM
*
*     RTRAN  = FALSE : skip R intermediate contributions
*     GTRAN  = FALSE : skip G intermediate contributions
*     FCKCON = FALSE : skip FOCK & ONEHAM matrix contribution
*     LRELAX = FALSE : skip contributions form XLAMDHB/XLAMDPB
*
*     Symmetries:    ISYMEI  --  EMAT1, EMAT2, FOCK, ONEHAM
*                    ISYLMA  --  XLAMDHA, XLAMDPA   
*                    ISYLMB  --  XLAMDHB, XLAMDPB   
*
*     Christof Haettig 20-6-1998
*
*---------------------------------------------------------------------*
#include <implicit.h>
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DIMENSION EMAT1(*),EMAT2(*),FOCK(*),ONEHAM(*)
      DIMENSION RAIM(*),RBIM(*),GAIM(*),GBIM(*)
      DIMENSION XLAMDHA(*),XLAMDHB(*), XLAMDPA(*), XLAMDPB(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      LOGICAL FCKCON, RTRAN, GTRAN, LRELAX
      INTEGER IOPT
C
      ISYAIM = MULD2H(ISYLMA,ISYMEI)
      ISYBIM = 0
      IF (IOPT.EQ.2) THEN
        ISYAIM = MULD2H(ISYLMB,ISYMEI)
        ISYBIM = MULD2H(ISYLMA,ISYMEI)
      END IF
C
C---------------------------------------------------------
C     Transform the delta index of R intermediate(s) to c.
C     store result in EMAT1
C---------------------------------------------------------
C
      IF ( RTRAN ) THEN
C
         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_EIM> norm^2 of RAIM:',
     &       DDOT(NEMAT1(ISYAIM),RAIM,1,RAIM,1)
           IF (NSYM.EQ.1) THEN
             WRITE (LUPRI,*) 'CC_EIM> RAIM:'
             CALL OUTPUT(RAIM,1,NVIR,1,NBAS,NVIR,NBAS,1,LUPRI)
           END IF
           IF (IOPT.EQ.2) THEN
             WRITE (LUPRI,*) 'CC_EIM> norm^2 of RBIM:',
     &         DDOT(NEMAT1(ISYBIM),RBIM,1,RBIM,1)
             IF (NSYM.EQ.1) THEN
               WRITE (LUPRI,*) 'CC_EIM> RBIM:'
               CALL OUTPUT(RBIM,1,NVIR,1,NBAS,NVIR,NBAS,1,LUPRI)
             END IF
           END IF
         END IF
C
         IF (IOPT.NE.1) THEN
           CALL DZERO(EMAT1,NMATAB(ISYMEI))
         END IF
C
         DO ISYMD = 1,NSYM
C
            ISYMC = MULD2H(ISYMD,ISYLMA)
            ISYMB = MULD2H(ISYMC,ISYMEI)
C
            NVIRB = MAX(NVIR(ISYMB),1)
            NBASD = MAX(NBAS(ISYMD),1)
C
            KOFF1 = IEMAT1(ISYMB,ISYMD) + 1
            KOFF2 = IGLMVI(ISYMD,ISYMC) + 1
            KOFF3 = IMATAB(ISYMB,ISYMC) + 1
C
            IF ( IOPT .EQ. 1) THEN

              CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMD),
     *                   -ONE,RAIM(KOFF1),NVIRB,XLAMDHA(KOFF2),NBASD,
     *                   ZERO,EMAT1(KOFF3),NVIRB)

            ELSE

              CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMD),
     *                   -ONE,RBIM(KOFF1),NVIRB,XLAMDHA(KOFF2),NBASD,
     *                    ONE,EMAT1(KOFF3),NVIRB)
C

              IF (LRELAX) THEN

                ISYMC = MULD2H(ISYMD,ISYLMB)
                ISYMB = MULD2H(ISYMC,ISYMEI)

                NVIRB = MAX(NVIR(ISYMB),1)
                NBASD = MAX(NBAS(ISYMD),1)
C
                KOFF1 = IEMAT1(ISYMB,ISYMD) + 1
                KOFF2 = IGLMVI(ISYMD,ISYMC) + 1
                KOFF3 = IMATAB(ISYMB,ISYMC) + 1
C
                CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMD),
     *                     -ONE,RAIM(KOFF1),NVIRB,XLAMDHB(KOFF2),NBASD,
     *                      ONE,EMAT1(KOFF3),NVIRB)
              END IF

            END IF
C
         END DO
C
      ELSE
C
         CALL DZERO(EMAT1,NMATAB(ISYMEI))
C
      END IF 
C
C---------------------------------------------------------
C     Transform the alpha index of G intermediate(s) to k.
C     store result in EMAT2
C---------------------------------------------------------
C
      IF ( GTRAN ) THEN
C
         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_EIM> norm^2 of GAIM:',
     &       DDOT(NT1AO(ISYAIM),GAIM,1,GAIM,1)
           IF (NSYM.EQ.1) THEN
             WRITE (LUPRI,*) 'CC_EIM> GAIM:'
             CALL OUTPUT(GAIM,1,NBAS,1,NRHF,NBAS,NRHF,1,LUPRI)
           END IF
           IF (IOPT.EQ.2) THEN
             WRITE (LUPRI,*) 'CC_EIM> norm^2 of GBIM:',
     &         DDOT(NT1AO(ISYBIM),GBIM,1,GBIM,1)
             IF (NSYM.EQ.1) THEN
               WRITE (LUPRI,*) 'CC_EIM> GBIM:'
               CALL OUTPUT(GBIM,1,NBAS,1,NRHF,NBAS,NRHF,1,LUPRI)
             END IF
           END IF
         END IF
C
         IF (IOPT.NE.1) THEN
           CALL DZERO(EMAT2,NMATIJ(ISYMEI))
         END IF
C
         DO ISYMA = 1,NSYM
C
            ISYMK = MULD2H(ISYMA,ISYLMA)
            ISYMJ = MULD2H(ISYMK,ISYMEI)
C
            NRHFK = MAX(NRHF(ISYMK),1)
            NBASA = MAX(NBAS(ISYMA),1)
C
            KOFF1 = IT1AO(ISYMA,ISYMJ)  + 1
            KOFF2 = IGLMRH(ISYMA,ISYMK) + 1
            KOFF3 = IMATIJ(ISYMK,ISYMJ) + 1
C
            IF ( IOPT .EQ. 1) THEN

              CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),NBAS(ISYMA),
     *                   ONE,XLAMDPA(KOFF2),NBASA,GAIM(KOFF1),NBASA,
     *                   ZERO,EMAT2(KOFF3),NRHFK)

            ELSE

              CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),NBAS(ISYMA),
     *                   ONE,XLAMDPA(KOFF2),NBASA,GBIM(KOFF1),NBASA,
     *                   ONE,EMAT2(KOFF3),NRHFK)
C
              IF (LRELAX) THEN
C
                ISYMK = MULD2H(ISYMA,ISYLMB)
                ISYMJ = MULD2H(ISYMK,ISYMEI)
C
                NRHFK = MAX(NRHF(ISYMK),1)
                NBASA = MAX(NBAS(ISYMA),1)
C
                KOFF1 = IT1AO(ISYMA,ISYMJ)  + 1
                KOFF2 = IGLMRH(ISYMA,ISYMK) + 1
                KOFF3 = IMATIJ(ISYMK,ISYMJ) + 1
C
                CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),NBAS(ISYMA),
     *                     ONE,XLAMDPB(KOFF2),NBASA,GAIM(KOFF1),NBASA,
     *                     ONE,EMAT2(KOFF3),NRHFK)
              END IF

            END IF
C
         END DO
C
      ELSE
C
         CALL DZERO(EMAT2,NMATIJ(ISYMEI))
C
      ENDIF
C
C---------------------------------------------------------
C     Add the Fock matrix / 1-el Hamiltonian contribution:
C---------------------------------------------------------
C
      IF (FCKCON) THEN
 
         DO ISYMC = 1,NSYM
 
            ISYMB = MULD2H(ISYMC,ISYMEI)
 
            DO C = 1,NVIR(ISYMC)
 
               KOFF1 = IFCVIR(ISYMB,ISYMC) + NORB(ISYMB)*(C - 1)
     *                                     + NRHF(ISYMB) + 1
               KOFF2 = IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C - 1) + 1
 
               CALL DAXPY(NVIR(ISYMB),ONE,FOCK(KOFF1),1,EMAT1(KOFF2),1)
 
            END DO
         END DO
 
         DO ISYMJ = 1,NSYM
 
            ISYMK = MULD2H(ISYMJ,ISYMEI)
 
            DO J = 1,NRHF(ISYMJ)
 
                KOFF1 = IFCRHF(ISYMK,ISYMJ) + NORB(ISYMK)*(J - 1) + 1
                KOFF2 = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(J - 1) + 1
 
                CALL DAXPY(NRHF(ISYMK),ONE,ONEHAM(KOFF1),1,
     *                                     EMAT2(KOFF2),1)
 
            END DO
         END DO

      ENDIF
 
      RETURN
      END
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_GAMMA2(GAMMA,EMAT2,ISYGAM)
*---------------------------------------------------------------------*
*
* Purpose: add EMAT2 to diagonals of the GAMMA intermediate
*
* C. Haettig August 1998
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
      PARAMETER(DMONE=-1.0D0)
      DIMENSION GAMMA(*),EMAT2(*)
C
      DO ISYMKI = 1, NSYM
        ISYMLJ = MULD2H(ISYGAM,ISYMKI)

        IF (ISYMKI.EQ.ISYMLJ .AND. ISYMKI.EQ.1) THEN

            DO ISYML = 1, NSYM
              DO L = 1, NRHF(ISYML)
                NLL   = IMATIJ(ISYML,ISYML) + NRHF(ISYML)*(L-1) + L
                N11LL = IGAMMA(ISYMKI,ISYMLJ) + NLL*(NLL-1)/2 + 1
                CALL DAXPY(NLL,DMONE,EMAT2,1,GAMMA(N11LL),1)
                DO NKI = NLL, NMATIJ(ISYMKI)
                  NLLKI = IGAMMA(ISYMLJ,ISYMKI) + NKI*(NKI-1)/2 + NLL
                  GAMMA(NLLKI) = GAMMA(NLLKI) + DMONE*EMAT2(NKI)
                END DO
              END DO
            END DO
             
        ELSE IF (ISYMKI.LT.ISYMLJ .AND. ISYMKI.EQ.1) THEN

            DO ISYMK = 1, NSYM
              DO K = 1, NRHF(ISYMK)
                NKK   = IMATIJ(ISYMK,ISYMK) + NRHF(ISYMK)*(K-1) + K
                NKK11 = IGAMMA(ISYMKI,ISYMLJ) + NKK
                CALL DAXPY(NMATIJ(ISYMLJ),DMONE,EMAT2,1,
     &                     GAMMA(NKK11),NMATIJ(ISYMKI))
              END DO
            END DO
             
        END IF

      END DO
 
      RETURN
      END
*=====================================================================*
*                   END OF SUBROUTINE CC_GAMMA2                       *
*=====================================================================*
c /* deck cc_iajb2 */
*=====================================================================*
      SUBROUTINE CC_IAJB2( XIAJB, ISYIAJB, IOPT, LRELAX, LOCC,
     &                     LU0IAJB, FN0IAJB, IT2DEL0, XLAMD0, ISYM0,
     &                     LUBIAJB, FNBIAJB, IT2DELB, XLAMDB, ISYMB,
     &                     WORK, LWORK                               )
*---------------------------------------------------------------------*
*
*   Purpose: transformation of (ia|j del) to (ia|jb) integrals
*            last partial transformation: (ia|j del) --> (ia|jb)
*            in a loop over delta index, (ia|j del) read from file
*
*            (ia|bj) += sum_del (ia|j del) * XLAMD(del b)
*
*            IOPT=2 : transform derivative integrals on LUBIAJB
*
*            if LRELAX, add relaxation contributions
*            if LOCC, transform to occupied index instead of a virtual:
*
*               (ia|kj) += sum_del (ia|j del) * XLAMD(del k)
*
*            for IOPT=1 LUBIAJB, FNBIAJB, ..., ISYMB is dummy input
*
*            also used to transform the integrals
*            (ij|del a) to (ij|ba) and (ia|del j) to (ia|bj)
*
*            NOTE: output arrays must be initialized from the outside!!
*
*    Written by Christof Haettig, June 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <mxorb_dec.h>
#  include <ccisao_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <mxorb.h>
#include <ccisao.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LRELAX, LOCC
      CHARACTER*(*) FN0IAJB, FNBIAJB
      INTEGER ISYIAJB, ISYM0, ISYMB, IOPT, LWORK, LU0IAJB, LUBIAJB
      INTEGER IT2DEL0(*), IT2DELB(*)

#if defined (SYS_CRAY)
      REAL ONE, ZERO, DDOT, XNORM
      REAL XIAJB(*), XLAMD0(*), XLAMDB(*), WORK(LWORK)
#else
      DOUBLE PRECISION ONE, ZERO, DDOT, XNORM
      DOUBLE PRECISION XIAJB(*), XLAMD0(*), XLAMDB(*), WORK(LWORK)
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
      
      INTEGER ISYMB0, ISYMBB, ISYDEL, ISY0IAJ, ISYBIAJ
      INTEGER LEN0,IADR0,LENB,IADRB,IDEL,KEND1,LWRK1,KX0IAJ,KXBIAJ

*---------------------------------------------------------------------*
* do the transformation in a loop over AO delta index:
*---------------------------------------------------------------------*
      DO ISYDEL = 1, NSYM

        ISYMB0  = MULD2H(ISYDEL,ISYM0)
C VB mar 2001:
C     If IOPT=1 ISYMB is a dummy argument at the statement below may
C     cause a segmentation fault. Moved inside IF-test instead
C
C        ISYMBB  = MULD2H(ISYDEL,ISYMB)

        IF (IOPT.EQ.1) THEN
          ISY0IAJ = MULD2H(ISYIAJB,ISYMB0)
        ELSE IF (IOPT.EQ.2) THEN
          ISYMBB  = MULD2H(ISYDEL,ISYMB)
          ISY0IAJ = MULD2H(ISYIAJB,ISYMBB)
          ISYBIAJ = MULD2H(ISYIAJB,ISYMB0)
        ELSE
          CALL QUIT('Illegal option in CC_IAJB2.')
        END IF

        KX0IAJ = 1
        KEND1  = KX0IAJ + NT2BCD(ISY0IAJ)
        IF (IOPT.EQ.2) THEN 
          KXBIAJ = KX0IAJ + NT2BCD(ISY0IAJ)
          KEND1  = KXBIAJ + NT2BCD(ISYBIAJ)
        END IF

        LWRK1  = LWORK - KEND1

        IF (LWRK1 .LT. 0) THEN
          WRITE (LUPRI,*) 'LWORK:',LWORK
          WRITE (LUPRI,*) 'LWRK1:',LWRK1
          CALL QUIT('Insufficient work space in CC_IAJB2.')
        END IF

        DO D = 1, NBAS(ISYDEL)

          IDEL  = D + IBAS(ISYDEL)


          IF (IOPT.EQ.1) THEN

C           ----------------------------------------------
C           (ia|jb) += (ia|j del) * XLAMD0(del b)
C           ----------------------------------------------
            LEN0  = NT2BCD(ISY0IAJ)
            IADR0 = IT2DEL0(IDEL)

            CALL GETWA2(LU0IAJB, FN0IAJB, WORK(KX0IAJ), IADR0, LEN0)
      
            CALL CC_IAJB3(WORK(KX0IAJ),ISY0IAJ,IDEL,LOCC,
     &                    XLAMD0,ISYM0,XIAJB)

          END IF

          IF (IOPT.EQ.2) THEN

C           ----------------------------------------------
C           (ia|jb)-bar += (ia|j del)-bar * XLAMD0(del b)
C           ----------------------------------------------
            LENB  = NT2BCD(ISYBIAJ)
            IADRB = IT2DELB(IDEL)

            CALL GETWA2(LUBIAJB, FNBIAJB, WORK(KXBIAJ), IADRB, LENB)

            CALL CC_IAJB3(WORK(KXBIAJ),ISYBIAJ,IDEL,LOCC,
     &                    XLAMD0,ISYM0,XIAJB)


            IF (LRELAX) THEN

C             ----------------------------------------------
C             (ia|jb)-bar += (ia|j del) * XLAMDB(del b-bar)
C             ----------------------------------------------
              LEN0  = NT2BCD(ISY0IAJ)
              IADR0 = IT2DEL0(IDEL)

              CALL GETWA2(LU0IAJB, FN0IAJB, WORK(KX0IAJ), IADR0, LEN0)

              CALL CC_IAJB3(WORK(KX0IAJ),ISY0IAJ,IDEL,LOCC,
     &                      XLAMDB,ISYMB,XIAJB)

            END IF
          END IF

        END DO

      END DO

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IAJB2                          *
*=====================================================================*
c /* deck cc_iajb3 */
*=====================================================================*
      SUBROUTINE CC_IAJB3(XIAJ,ISYAIJ,IDEL,LOCC,XLAMDA,ISYLAM,XIAJB)
*---------------------------------------------------------------------*
*
*   Purpose: transformation to (ia|jb) integrals
*            last partial transformation: (ia|j del) --> (ia|jb)
*
*            (ia|bj) += sum_del (ia|j del) * XLAMDA(del b)
*
*            if LOCC.EQ.TRUE transform to occupied instead:
*               (ia|kj) += sum_del (ia|j del) * XLAMDA(del k)
*
*            also used to transform the integrals
*            (ij|del a) to (ij|ba) and (ia|del j) to (ia|bj)
*
*    Written by Christof Haettig, May 1998.
*    LOCC option added in September 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <mxorb_dec.h>
#  include <ccisao_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <mxorb.h>
#include <ccisao.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

#if defined (SYS_CRAY)
      REAL ONE, ZERO
#else
      DOUBLE PRECISION ONE, ZERO
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      INTEGER IDEL, ISYAIJ, ISYLAM
      
#if defined (SYS_CRAY)
      REAL XLAMDA(*), XIAJ(*), XIAJB(*)
#else
      DOUBLE PRECISION XLAMDA(*), XIAJ(*), XIAJB(*)
#endif
      
      LOGICAL LOCC
      INTEGER ITAIKJ(8,8), NTAIKJ(8), ISYM, ICOUNT
      INTEGER ISYDEL, ISYMB, ISYMJ, ISYMAI, ISYMBJ, ISYMK, ISYMAIK
      INTEGER KOFF0, KOFF1, KLAMD, NBJ

*---------------------------------------------------------------------*
*     set D index and its symmetry and symmetry of B index:
*---------------------------------------------------------------------*
      ISYDEL = ISAO(IDEL)
      ISYMB  = MULD2H(ISYDEL,ISYLAM)
      ISYMK  = MULD2H(ISYDEL,ISYLAM)

      DO  ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO 
        NTAIKJ(ISYM) = ICOUNT
      END DO

      D = IDEL - IBAS(ISYDEL)

*---------------------------------------------------------------------*
*     (ia|jb) = sum_{del} (ia|del j) * XLAMDA(del b)
* or  (ia|jk) = sum_{del} (ia|del j) * XLAMDA(del k)
*---------------------------------------------------------------------*
      DO ISYMJ = 1, NSYM

        ISYMAI  = MULD2H(ISYAIJ,ISYMJ)
        ISYMBJ  = MULD2H(ISYMB,ISYMJ)
        ISYMAIK = MULD2H(ISYMAI,ISYMK)

        DO J = 1, NRHF(ISYMJ)

          KOFF0 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1) + 1

          IF (LOCC) THEN

             DO K = 1, NRHF(ISYMK)
          
               KLAMD = IGLMRH(ISYDEL,ISYMK) + NBAS(ISYDEL)*(K-1) + D
               KOFF1 = ITAIKJ(ISYMAIK,ISYMJ) + NT2BCD(ISYMAIK)*(J-1)
     &               + IT2BCD(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1) + 1

               CALL DAXPY(NT1AM(ISYMAI),XLAMDA(KLAMD),XIAJ(KOFF0),1,
     &                                         XIAJB(KOFF1),1)

             END DO

          ELSE

             DO B = 1, NVIR(ISYMB)
          
               NBJ   = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
               KLAMD = IGLMVI(ISYDEL,ISYMB) + NBAS(ISYDEL)*(B-1) + D
               KOFF1 = IT2SQ(ISYMAI,ISYMBJ) + NT1AM(ISYMAI)*(NBJ-1) + 1

               CALL DAXPY(NT1AM(ISYMAI),XLAMDA(KLAMD),XIAJ(KOFF0),1,
     &                                         XIAJB(KOFF1),1)

             END DO

          END IF

        END DO

      END DO

      RETURN
      END
*=====================================================================*
*                 END OF SUBROUTINE CC_IAJB3                          *
*=====================================================================*
C  /* Deck cc_t2pk */
      SUBROUTINE CC_T2PK(T2PK,T2SQ,ISYAMP,IOPT)
C
C     Pack squared matrix of T2 amplitudes to triangular
C
C     if IOPT.EQ.1 apply P^{ab}_{ij}: 
C            T2PK(ia,jb) = T2SQ(ia,jb) + T2SQ(jb,ia)
C
C     Ch. Haettig, August 1998
C
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0)
      DIMENSION T2PK(*), T2SQ(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      IF (ISYAMP.EQ.1) THEN
C
        DO ISYMAI = 1, NSYM
           DO I = 1, NT1AM(ISYMAI)
              KOFFSQ = IT2SQ(ISYMAI,ISYMAI)+NT1AM(ISYMAI)*(I-1)+1
              KOFFPK = IT2AM(ISYMAI,ISYMAI)+I*(I-1)/2 + 1
              CALL DCOPY(I,T2SQ(KOFFSQ),1,T2PK(KOFFPK),1)
           END DO
        END DO
C
        IF (IOPT.EQ.1) THEN
          DO ISYMAI = 1, NSYM
            DO I = 1, NT1AM(ISYMAI)
              KOFFSQ = IT2SQ(ISYMAI,ISYMAI)+I
              KOFFPK = IT2AM(ISYMAI,ISYMAI)+I*(I-1)/2 + 1
              CALL DAXPY(I,ONE,T2SQ(KOFFSQ),NT1AM(ISYMAI),
     &                         T2PK(KOFFPK),1)
            END DO
          END DO
        END IF
C
      ELSE
C
        DO ISYMAI = 1, NSYM
           ISYMBJ = MULD2H(ISYMAI,ISYAMP)
           IF (ISYMBJ.GT.ISYMAI) THEN
             KOFFSQ = IT2SQ(ISYMAI,ISYMBJ) + 1
             KOFFPK = IT2AM(ISYMAI,ISYMBJ) + 1
             NAMP   = NT1AM(ISYMAI) * NT1AM(ISYMBJ)
             CALL DCOPY(NAMP,T2SQ(KOFFSQ),1,T2PK(KOFFPK),1)
           END IF
        END DO
C
        IF (IOPT.EQ.1) THEN
          DO ISYMAI = 1, NSYM
            ISYMBJ = MULD2H(ISYMAI,ISYAMP)
            IF (ISYMBJ.GT.ISYMAI) THEN
              DO I = 1, NT1AM(ISYMBJ)
                KOFFSQ = IT2SQ(ISYMBJ,ISYMAI)+I
                KOFFPK = IT2AM(ISYMAI,ISYMBJ)+NT1AM(ISYMAI)*(I-1)+1
                CALL DAXPY(NT1AM(ISYMAI),ONE,T2SQ(KOFFSQ),NT1AM(ISYMBJ),
     &                                       T2PK(KOFFPK),1)
              END DO
            END IF
          END DO
        END IF
C
      END IF
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ccrhs_cio2 */
      SUBROUTINE CCRHS_CIO3(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                      ISYVEC,ISYCIM,LUC,CFIL,IV,IOPT,
     *                      IOPTB,LUCB,CBFIL,IOFFB,FACB,
     *                      IOPTE,EMAT1)
C
C
C      vectors        ISYVEC
C      intermediates  ISYCIM
C
C      LUC and CFIL is used to control from which file the
C                   intermediate is obtained.
C
C      if iopt = 1 the C intermediate is assumed
C         to be as in energy calc.
C
C      if iopt ne. 1 we use the intermediate
C         on luc with address given according to
C         transformed vector nr iv (iv is not 1
C         if several vectors are transformed
C         simultaneously.)
C
C      if iopt = 3 we return the C intermediate in OMEGA2
C
C      if iopte = 1 add contribution of E1 intermediate to diagonal
C
C      extra contribution in B matrix calculation:
C       if IOPTB=1 contribution form CBAR intermediate is included. 
C       the intermediate is read from file CBFIL (unit LUCB) with
C       the offset IOFFB
C
C      in energy calc: iv=1,iopt=1
C
C
C     PURPOSE:
C             Calculate the C-term 
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*),EMAT1(*)
      DIMENSION WORK(LWORK)
      CHARACTER CFIL*(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      XNORM = 0.0d0
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*) 'I/O in C-term not implemented for '//
     &        'square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_CIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYCIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYCIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
C
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for allocation '//
     &           'in CCRHS_CIO-1')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KSCR3 = KSCR2 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR3 + NT1AM(ISYMCK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for allocation '//
     &              'in CCRHS_CIO-2')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUC,CFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C-----------------------------------------
C           Transform delta index to c.
C-----------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE
C
C------------------------------------------
C         add E1 to contrib. to R diagonal:
C------------------------------------------
C
          IF (IOPTE .EQ. 1) THEN

             DO ISYMC = 1, NSYM

               ISYMK  = MULD2H(ISYMC,ISYMCK)
               ISYMA  = MULD2H(ISYMK,ISYMAI)

               DO A = 1, NVIR(ISYMA)
               DO K = 1, NRHF(ISYMK)

                 NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K-1) + A
                 KAK = NAK - IFSTAI + 1

                 IF (NAK.GE.IFSTAI .AND. NAK.LE.ILSTAI) THEN

                   DO C = 1, NVIR(ISYMC)

                      NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) + C

                      KOFF1 = IMATAB(ISYMA,ISYMC)+NVIR(ISYMA)*(C-1)+A

                      KOFF6 = KSCR1-1 + NT1AM(ISYMCK)*(KAK - 1) + NCK

                      WORK(KOFF6) = WORK(KOFF6) - 0.5D0 * EMAT1(KOFF1) 

                   END DO

                 END IF

               END DO
               END DO

             END DO

          END IF
C
C---------------------------------------------
C         add FACB * CBAR to P(ck,#ai)
C---------------------------------------------
C
          IF (IOPTB .EQ. 1) THEN
            DO NAI = IFSTAI, ILSTAI

              IADR = IOFFB+IT2SQ(ISYMCK,ISYMAI)+NT1AM(ISYMCK)*(NAI-1)+1

              CALL GETWA2(LUCB,CBFIL,WORK(KSCR3),IADR,NT1AM(ISYMCK))

              KAI   = NAI - IFSTAI + 1
              KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI-1)

             CALL DAXPY(NT1AM(ISYMCK),FACB,WORK(KSCR3),1,WORK(KOFF6),1)

            END DO
          END IF
C
C----------------------------------------------------
C         if IOPT = 3 store C intermediate in OMEGA2:
C----------------------------------------------------
C
          IF (IOPT.EQ.3) THEN

            DO NAI = IFSTAI, ILSTAI

              KAI   = NAI - IFSTAI + 1
              KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI-1)
              KOFF1 = IT2SQ(ISYMCK,ISYMAI)+NT1AM(ISYMCK)*(NAI-1)+1

              CALL DCOPY(NT1AM(ISYMCK),WORK(KOFF6),1,OMEGA2(KOFF1),1) 

            END DO

          ELSE
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 190 I = 1,NRHF(ISYMI)
C
                  DO 200 A = 1,NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     IF ((NAI .LT. IFSTAI) .OR. (NAI .GT. ILSTAI))
     *                  GOTO 200
C
                     DO 210 ISYMJ = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)
C
                        DO 220 J = 1,NRHF(ISYMJ)
C
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J-1) + A
C
                           CALL CC_PUTC(WORK(KSCR2),OMEGA2,ISYMAI,
     *                                  ISYMAJ,ISYMBI,ISYMBJ,ISYMB,
     *                                  ISYMI,ISYMJ,NAI,NAJ,I,J,
     *                                  IFSTAI)
C
  220                   CONTINUE
  210                CONTINUE
  200             CONTINUE
  190          CONTINUE
  180       CONTINUE
C
          END IF
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_dio2 */
      SUBROUTINE CCRHS_DIO3(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                      ISYVEC,ISYDIM,LUD,DFIL,LUC,CFIL,IV,IOPT,
     *                      IOPTB,LUDB,DBFIL,IOFFB,FACB,
     *                      IOPTE,EMAT1)
C
C     asm 20-aug-1994
C
C      vectors        ISYVEC
C      intermediates  ISYDIM 
C      LUD and DFIL are used to control from which file the
C      intermediate is obtained.
C
C      if iopt = 1 the D intermediate is assumed
C         to be as in energy calc.
C
C      if iopt ne. 1 we use the intermediate
C         on luc with address given according to
C         transformed vector nr iv (iv is not 1
C         if several vectors are transformed
C         simultaneously.)
C
C      if iopt = 3 return D intermediate in OMEGA2
C
C      it iopte = 1 add contribution of E1 intermediate to diagonal
C
C      extra contribution in B matrix calculation:
C       if IOPTB=1 contribution from DBAR 
C       intermediate in included. the intermediate
C       is read from file DBFIL (unit LUDB) with
C       the offset IOFFB
C
C      in energy calc: iv=1,iopt=1
C
C     PURPOSE:
C             Calculate the D-term 
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*),EMAT1(*)
      DIMENSION WORK(LWORK)
      CHARACTER DFIL*(*), DBFIL*(*), CFIL*(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccsdio.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      XNORM = 0.0d0
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*) 'I/O in D-term not implemented for '//
     &        'square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_DIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYDIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYDIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI 
         IF (LENMIN .EQ. 0) GOTO 100
         LENMIN = LENMIN + NRHFT
C
         NDISAI = ( LWORK - NT1AM(ISYMCK) ) / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_DIO2')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KSCR3 = KSCR2 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR3 + NT1AM(ISYMCK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for allocation '//
     &              'in CCRHS_DIO')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  KSCR4 = KEND
                  IF (LWRK1 .LT. LEN) THEN
                     CALL QUIT('Insufficient memory in CCRHS_DIO2')
                  END IF
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUD,DFIL,WORK(KOFF1),IOFF,LEN)
                     CALL GETWA2(LUC,CFIL,WORK(KSCR4),IOFF,LEN)
C
                     CALL DSCAL(LEN,TWO,WORK(KOFF1),1)
                     CALL DAXPY(LEN,-ONE,WORK(KSCR4),1,WORK(KOFF1),1)
C
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C--------------------------------------
C           Transform delta index to c.
C--------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE
C
C--------------------------------------------------
C           contract to contrib. to R intermediate:
C--------------------------------------------------
C
          IF (IOPTE .EQ. 1) THEN
             
             DO ISYMC = 1, NSYM

               ISYMK  = MULD2H(ISYMC,ISYMCK)
               ISYMA  = MULD2H(ISYMK,ISYMAI)

               DO A = 1, NVIR(ISYMA)
               DO K = 1, NRHF(ISYMK)

                 NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K-1) + A
                 KAK = NAK - IFSTAI + 1

                 IF (NAK.GE.IFSTAI .AND. NAK.LE.ILSTAI) THEN

                   DO C = 1, NVIR(ISYMC)

                      NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) + C

                      KOFF1 = IMATAB(ISYMA,ISYMC)+NVIR(ISYMA)*(C-1)+A

                      KOFF6 = KSCR1-1 + NT1AM(ISYMCK)*(KAK - 1) + NCK

                      WORK(KOFF6) = WORK(KOFF6) + 0.5D0 * EMAT1(KOFF1)

                   END DO

                 END IF

               END DO
               END DO

             END DO

          END IF
C
C---------------------------------------------
C         add DBAR intermediate to P(ck,#ai)
C---------------------------------------------
C
          IF (IOPTB .EQ. 1) THEN
            DO NAI = IFSTAI, ILSTAI

              IADR = IOFFB+IT2SQ(ISYMCK,ISYMAI)+NT1AM(ISYMCK)*(NAI-1)+1

              CALL GETWA2(LUDB,DBFIL,WORK(KSCR3),IADR,NT1AM(ISYMCK))

              KAI   = NAI - IFSTAI + 1
              KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI-1)

             CALL DAXPY(NT1AM(ISYMCK),FACB,WORK(KSCR3),1,WORK(KOFF6),1)
C
            END DO
          END IF
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
          IF (IOPT.EQ.3) THEN
         
            DO NAI = IFSTAI, ILSTAI

              KAI   = NAI - IFSTAI + 1
              KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI-1)
              KOFF1 = IT2SQ(ISYMCK,ISYMAI)+NT1AM(ISYMCK)*(NAI-1)+1
C
              CALL DCOPY(NT1AM(ISYMCK),WORK(KOFF6),1,OMEGA2(KOFF1),1) 

            END DO

          ELSE
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 NAI = IFSTAI,ILSTAI
C
               CALL CC_PUTD(WORK(KSCR2),OMEGA2,ISYMAI,ISYMBJ,NAI,IFSTAI)
C
  180       CONTINUE
C
          END IF
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
*---------------------------------------------------------------------*
c/* deck cc_cd */
*=====================================================================*
      SUBROUTINE CC_CD(TYPE, ISIDE, OMEGA2, ISYRES, T2AMP, ISYTAM,
     &                 CDBAR, ISYCDBAR, WORK, LWORK        )
*---------------------------------------------------------------------*
*    Purpose: calculate C/D term from C/D in core
*
*    TYPE='C' : calculate C intermediate
*    TYPE='D' : calculate D intermediate
*
*    ISIDE = +1 : contract as needed in Jacobian right transformation
*    ISIDE = -1 : contract as needed in Jacobian left transformation
*
*    N.B. for TYPE='D' and ISIDE=+1 the amplitudes T2AMP assumed to be 
*                                   transformed to 2*t(ai,bj)-t(aj,bi)
*         for TYPE='C' and ISIDE=-1 the amplitudes T2AMP are over-
*                                   written by 2t(aj,bi)+t(ai,bj)
*      
*    Written by Christof Haettig, June 1998.
*    modified for left side stuff in August 1998.
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <ccsdsym_dec.h>
#  include <ccorb_dec.h>
#include <priunit_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccorb.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*(*) TYPE
      INTEGER LWORK, ISYTAM, ISYCDBAR, ISYRES

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL T2AMP(*)      ! dimension (NT2AM(ISYTAM))
      REAL OMEGA2(*)     ! dimension (NT2AM(ISYRES))
      REAL CDBAR(*)      ! dimension (NT2SQ(ISYCDBAR))
      REAL HALF, ONE, TWO, XNORM, DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION T2AMP(*)      ! dimension (NT2AM(ISYTAM))
      DOUBLE PRECISION OMEGA2(*)     ! dimension (NT2AM(ISYRES))
      DOUBLE PRECISION CDBAR(*)      ! dimension (NT2SQ(ISYCDBAR))
      DOUBLE PRECISION HALF, ONE, TWO, XNORM, DDOT
#endif
      PARAMETER (ONE = 1.0d0, HALF = 0.5D0, TWO = 2.0D0)

      INTEGER IOPTZWVI, ISYTIN, ISYCIN, ISIDE
      INTEGER ISYMA, ISYOME, ISYMI, ISYMBJ, ISYMAI
      INTEGER NAI, NI, NBJ, NBJAI, NAIAI, NTNI
      INTEGER KEND1, LWRK1, KTINT, KOMEG, KSCR, INDEX

* statement functions:
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

*---------------------------------------------------------------------*
* check symmetries:
*---------------------------------------------------------------------*
      IF (MULD2H(ISYTAM,ISYCDBAR) .NE. ISYRES) THEN
        WRITE (LUPRI,*) 'ERROR> SYMMETRY MISMATCH IN CC_CD.'
        CALL QUIT('SYMMETRY MISMATCH IN CC_CD.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Entered CDBAR: TYPE, LWORK:',TYPE,LWORK
        CALL FLSHFO(LUPRI)
        XNORM =  0.0d0
      END IF

*---------------------------------------------------------------------*
* set option for CC_ZWVJ routine:
*---------------------------------------------------------------------*
      IF      (TYPE(1:1).EQ.'C' .AND. ISIDE.EQ.+1) THEN

        IOPTZWVI = 2

      ELSE IF (TYPE(1:1).EQ.'D' .AND. ISIDE.EQ.+1) THEN

        IOPTZWVI = 1

      ELSE IF (TYPE(1:1).EQ.'C' .AND. ISIDE.EQ.-1) THEN

        IOPTZWVI = 1

        CALL CCSD_Z2TR(T2AMP,ISYTAM)
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) ' (2EXC + COU) CTR2:'
          CALL CC_PRP(WORK,T2AMP,ISYTAM,0,1)
        END IF

      ELSE IF (TYPE(1:1).EQ.'D' .AND. ISIDE.EQ.-1) THEN

        IOPTZWVI = 1

      ELSE
        CALL QUIT('ILLEGAL OPTION IN CC_CD.')
      END IF

*---------------------------------------------------------------------*
* start loop over virtual orbital index a:
*---------------------------------------------------------------------*
      DO ISYMA = 1, NSYM
        ISYTIN = MULD2H(ISYTAM,ISYMA)
        ISYOME = MULD2H(ISYTIN,ISYCDBAR)

        KTINT = 1
        KOMEG = KTINT + NT2BCD(ISYTIN)
        KSCR  = KOMEG + NT2BCD(ISYOME)
        KEND1 = KSCR  + NT2BCD(ISYOME)

        LWRK1 = LWORK - KEND1
        IF (LWRK1 .LE. 0) THEN
          CALL QUIT('Insufficient work space in CCB_CDBAR.')
        END IF

      DO A = 1, NVIR(ISYMA)
        
C       ----------------------------------------------
C       get t^{a}(bj,i) submatrix of the t amplitudes:
C       ----------------------------------------------
        CALL CCG_TI(WORK(KTINT),ISYTIN,T2AMP,ISYTAM,A,ISYMA)

C       ----------------------------------------
C       calculate OMEGA2^{a}(ck,i) intermediate:
C       ----------------------------------------
        IF (ISIDE.EQ.+1) THEN
          Call CC_ZWVJ(WORK(KOMEG), CDBAR, ISYCDBAR, WORK(KTINT),
     &                 ISYTIN, WORK(KEND1), LWRK1, IOPTZWVI)
        ELSE IF (ISIDE.EQ.-1) THEN
          Call CC_ZWVI(WORK(KOMEG), CDBAR, ISYCDBAR, WORK(KTINT),
     &                 ISYTIN, WORK(KEND1), LWRK1, IOPTZWVI)
        ELSE
          CALL QUIT('ILLEGAL OPTION IN CC_CD.')
        END IF


C       ---------------------------------------------------
C       for D term in right transformation scale with +1/2:
C       ---------------------------------------------------
        IF ( TYPE(1:1).EQ.'D' .AND. ISIDE.EQ.+1 ) THEN
           CALL DSCAL(NT2BCD(ISYOME),HALF,WORK(KOMEG),1)
        END IF

C       ------------------------------------------------------
C       for C term in left transformation scale with -1/2 Pij:
C       ------------------------------------------------------
        IF ( TYPE(1:1).EQ.'C' .AND. ISIDE.EQ.-1 ) THEN
           CALL DSCAL(NT2BCD(ISYOME),-HALF,WORK(KOMEG),1)
           CALL CCLT_P21I(WORK(KOMEG),ISYOME,WORK(KEND1),LWRK1)
        END IF

C       --------------------------------------------------------------
C       for C term in right transformation apply -1/2 * (1 + 2 * Pij):
C       --------------------------------------------------------------
        IF ( TYPE(1:1).EQ.'C' .AND. ISIDE.EQ.+1) THEN
           CALL DCOPY(NT2BCD(ISYOME),WORK(KOMEG),1,WORK(KSCR),1)
           CALL CCLT_P21I(WORK(KSCR),ISYOME,WORK(KEND1),LWRK1)
           CALL DAXPY(NT2BCD(ISYOME),TWO,WORK(KSCR),1,WORK(KOMEG),1)
           CALL DSCAL(NT2BCD(ISYOME),-HALF,WORK(KOMEG),1)
        END IF

C       --------------------------------------------------------
C       for D term in left transformation apply (1 - 1/2 * Pij):
C       --------------------------------------------------------
        IF ( TYPE(1:1).EQ.'D' .AND. ISIDE.EQ.-1) THEN
           CALL DCOPY(NT2BCD(ISYOME),WORK(KOMEG),1,WORK(KSCR),1)
           CALL CCLT_P21I(WORK(KSCR),ISYOME,WORK(KEND1),LWRK1)
           CALL DAXPY(NT2BCD(ISYOME),-HALF,WORK(KSCR),1,WORK(KOMEG),1)
        END IF

C       -----------------------
C       store in result vector:
C       -----------------------
        DO ISYMI = 1, NSYM
          ISYMBJ = MULD2H(ISYOME,ISYMI)
          ISYMAI = MULD2H(ISYMA,ISYMI)

          IF ( MULD2H(ISYMAI,ISYMBJ) .NE. ISYRES ) THEN
             CALL QUIT('Symmetry mismatch in CC_CD.')
          END IF

          DO I = 1, NRHF(ISYMI)

            NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)+A
            NI  = KOMEG + IT2BCD(ISYMBJ,ISYMI)+NT1AM(ISYMBJ)*(I-1)

c           WRITE (LUPRI,*)TYPE,'-term: ISYMA,A,ISYMI,I=',ISYMA,A,ISYMI,I
c           NTNI = KTINT + IT2BCD(ISYMBJ,ISYMI)+NT1AM(ISYMBJ)*(I-1)
c           WRITE (LUPRI,*) 'T amplitudes:'
c           CALL CC_PRP(WORK(NTNI),WORK,ISYMBJ,1,0)
c           WRITE (LUPRI,*) 'contribution:'
c           CALL CC_PRP(WORK(NI),WORK,ISYMBJ,1,0)

            IF      ( ISYMBJ .EQ. ISYMAI ) THEN

               DO NBJ = 1, NT1AM(ISYMBJ)
                 NBJAI = IT2AM(ISYMBJ,ISYMAI) + INDEX(NBJ,NAI)
                 OMEGA2(NBJAI) = OMEGA2(NBJAI) + WORK(NI-1+NBJ)
               END DO
               NAIAI = IT2AM(ISYMAI,ISYMAI) + INDEX(NAI,NAI)
               OMEGA2(NAIAI) = OMEGA2(NAIAI) + WORK(NI-1+NAI)

            ELSE IF ( ISYMBJ .LT. ISYMAI ) THEN

               NBJAI = IT2AM(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)*(NAI-1)+1
               CALL DAXPY(NT1AM(ISYMBJ),ONE,WORK(NI),1,
     &                                     OMEGA2(NBJAI),1)

            ELSE IF ( ISYMBJ .GT. ISYMAI ) THEN
               
               NBJAI = IT2AM(ISYMAI,ISYMBJ) + NAI
               CALL DAXPY(NT1AM(ISYMBJ),ONE,WORK(NI),1,
     &                                     OMEGA2(NBJAI),NT1AM(ISYMAI))

            END IF
     

          END DO
        END DO


      END DO
      END DO


      RETURN
      END
*=====================================================================*
*                      END OF SUBROUTINE CC_CD
*=====================================================================*
*======================================================================
      SUBROUTINE CC_PQIMO(CTR2,ISYCTR,T2AM,ISYTAMP,YINT,IOPTY,
     &                    FILNAM,LUPQIM,IADRPQ,IADR,IV,WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: Calculate the P and Q intermediates in MO basis from the 
*              Lagrangian multipiers CTR2 and the amplitude vector T2AM
*              and write them to a file direct access file FILNAM
*              
*   P{aik,c} = sum_{dl} (2 t(dl,ck) - t(dk,l;del)) Zeta(dl;ai)
*   Q{aik,c} = sum_{dl} t(dl,ck) (2Zeta(di;al) + Zeta(dl;ai))
*
*   for IOPTY = 1 the Y contrib. are substracted from i=k diagonal:
*
*           P{aii,c} = P{aii,c} - 1 * Y{ca}
*           Q{aii,c} = P{aii,c} - 3 * Y{ca}
*
*   IADR contains on INPUT the start address of the P & Q intermed.
*   on file; on OUTPUT it will contain the start address for the
*   next set of P & Q intermediates.
*
*   Christof Haettig, September/October 1998
*
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <ccorb_dec.h>
#  include <maxorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <iratdef_dec.h>
#include <priunit_dec.h>
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <iratdef.h>
#include <priunit.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK, LUPQIM, IOPTY
      CHARACTER*(*) FILNAM

#if defined (SYS_CRAY)
      REAL ONE, TWO, THREE, HALF, XNORM, DDOT, ZERO, YCA
      REAL WORK(*), CTR2(*), T2AM(*), YINT(*)
      REAL SECOND, TIMET, TIMIO, TIMZWV, TIMTAM, DTIME, TIMTI
      REAL TIMZET,TIMSCL
#else
      DOUBLE PRECISION ONE, TWO, THREE, HALF, XNORM, DDOT, ZERO, YCA
      DOUBLE PRECISION WORK(*), CTR2(*), T2AM(*), YINT(*)
      DOUBLE PRECISION SECOND,TIMET,TIMIO,TIMTI,TIMZWV,TIMTAM,DTIME
      DOUBLE PRECISION TIMZET,TIMSCL
#endif
      PARAMETER(ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,HALF=0.5D0,ZERO=0.0D0)

      INTEGER IV, IADR, IADRW
      INTEGER IADRPQ(MAXORB,IV)
      INTEGER ISYTAMP,ISYCTR,ISYMC,ISYTIN,ISYAIK, ISYMAI,ISYMA,ISYMI
      INTEGER KEND1, LWRK1, KOFF, KOFFY, KYIM
      INTEGER IERR, IOPT, IVIRC
      INTEGER KTINT, KTJNT, KPINT, KQINT
     

*----------------------------------------------------------------------
* initialize timings & work space:
*----------------------------------------------------------------------
      TIMET  = SECOND()
      TIMIO  = ZERO
      TIMTI  = ZERO
      TIMZWV = ZERO
      TIMTAM = ZERO
      TIMZET = ZERO
      TIMSCL = ZERO

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'Norm of T2AMP:',
     *                   DDOT(NT2AM(ISYTAMP),T2AM,1,T2AM,1)
c        WRITE (LUPRI,*) 'T2AMP:'
c        CALL CC_PRP(WORK,T2AM,ISYTAMP,0,1)
         WRITE (LUPRI,*) 'Norm of CTR2 :',
     *                   DDOT(NT2SQ(ISYCTR),CTR2,1,CTR2,1)
c        WRITE (LUPRI,*) 'CTR2:'
c        CALL CC_PRSQ(WORK,CTR2,ISYCTR,0,1)
      END IF
*----------------------------------------------------------------------
* calculate the P intermediate in a loop over MO index c:
*----------------------------------------------------------------------
      DO ISYMC = 1, NSYM
      DO C = 1, NVIR(ISYMC)
        IVIRC = IVIR(ISYMC) + C
    
        ISYTIN = MULD2H(ISYMC,ISYTAMP)
        ISYAIK = MULD2H(ISYTIN,ISYCTR)
    
        KTINT = 1
        KTJNT = KTINT + NT2BCD(ISYTIN)
        KPINT = KTJNT + NT2BCD(ISYTIN)
        KQINT = KPINT + NT2BCD(ISYAIK)
        KYIM  = KQINT + NT2BCD(ISYAIK)
        KEND1 = KYIM  + NVIRT
        LWRK1 = LWORK - KEND1

        IF  (LWRK1 .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient memory in CC_PQIM.'
           CALL QUIT('Insufficient memory in CC_PQIM.')
        END IF
C
C       get c batch of T amplitudes:
C
        DTIME = SECOND()
        CALL CCG_TI(WORK(KTINT),ISYTIN,T2AM,ISYTAMP,C,ISYMC)
        TIMTI = TIMTI + SECOND() - DTIME
C
C       calculate 2 x Coul - Exc combination of backtransformed T:
C
        DTIME = SECOND()
        CALL DCOPY(NT2BCD(ISYTIN),WORK(KTINT),1,WORK(KTJNT),1)

        CALL DSCAL(NT2BCD(ISYTIN),-ONE,WORK(KTJNT),1)

        CALL CCLT_P21I(WORK(KTJNT),ISYTIN,WORK(KEND1),LWRK1)
        
        CALL DAXPY(NT2BCD(ISYTIN),TWO,WORK(KTINT),1,WORK(KTJNT),1)
        TIMTAM = TIMTAM + SECOND() - DTIME
C
C       calculate the P intermediate:
C
        IOPT = 1
        DTIME = SECOND()
        CALL CC_ZWVI(WORK(KPINT),CTR2,ISYCTR,WORK(KTJNT),
     &               ISYTIN,WORK(KEND1),LWRK1,IOPT)
        TIMZWV = TIMZWV + SECOND() - DTIME
C
        IF (IOPTY.EQ.1) THEN
C
C       substract contributions from Y intermediate to i=j diagonal:
C
          ISYMA = ISYAIK
          KOFFY = IMATAB(ISYMC,ISYMA) + C
          CALL DCOPY(NVIR(ISYMA),YINT(KOFFY),NVIR(ISYMC),WORK(KYIM),1)

          DO ISYMI = 1, NSYM
           DO I = 1, NRHF(ISYMI)

            ISYMAI = MULD2H(ISYMA,ISYMI)
            KOFF   = KPINT + IT2BCD(ISYMAI,ISYMI) + IT1AM(ISYMA,ISYMI) 
     &             + NT1AM(ISYMAI)*(I-1) + NVIR(ISYMA)*(I-1) 

            CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KYIM),1,WORK(KOFF),1)

           END DO
          END DO
        END IF
C
C       scale with 1/2 for compatibility with Asgers routines
C
        DTIME = SECOND()
        CALL DSCAL(NT2BCD(ISYAIK),HALF,WORK(KPINT),1)
        TIMSCL = TIMSCL + SECOND() - DTIME
C
C       write the intermediate to file:
C
        IADRPQ(IVIRC,IV) = IADR

        DTIME = SECOND()
        CALL  PUTWA2(LUPQIM,FILNAM,WORK(KPINT),IADR,NT2BCD(ISYAIK))
        TIMIO = TIMIO + SECOND() - DTIME

        IADR = IADR + 2*NT2BCD(ISYAIK)

      END DO
      END DO

*----------------------------------------------------------------------
*     Note that IADR contains now the start address for the next
*     set of P & Q intermediates. It must not be changed in the 
*     rest of this subroutine!
*----------------------------------------------------------------------

*----------------------------------------------------------------------
*     calculate (2 x Exc + Coul)/3 combination of ZETA:
*     (we interrupt here the loop over AO to calculate the modified
*      ZETA and accept that we recalculate the backtransformed 
*      amplitude since the transformation of the amplitudes is less
*      expansive than the transformation and reconstruction of ZETA 
*      for each delta. both the transformation of TAMP and the 
*      transformation/reconstruction of ZETA scale with N V^2 O^2. )
*----------------------------------------------------------------------
      DTIME = SECOND()
      CALL CCRHS_T2BT(CTR2,WORK,LWORK,ISYCTR)
      CALL CCSD_T2TP(CTR2,WORK,LWORK,ISYCTR)
      TIMZET = TIMZET + SECOND() - DTIME

*----------------------------------------------------------------------
* calculate the Q intermediate in a loop over MO index c:
*----------------------------------------------------------------------
      DO ISYMC = 1, NSYM
      DO C = 1, NVIR(ISYMC)
        IVIRC = IVIR(ISYMC) + C

        ISYTIN = MULD2H(ISYMC,ISYTAMP)
        ISYAIK = MULD2H(ISYTIN,ISYCTR)
    
        KTINT = 1
        KTJNT = KTINT + NT2BCD(ISYTIN)
        KPINT = KTJNT + NT2BCD(ISYTIN)
        KQINT = KPINT + NT2BCD(ISYAIK)
        KYIM  = KQINT + NT2BCD(ISYAIK)
        KEND1 = KYIM  + NVIRT
        LWRK1 = LWORK - KEND1

        IF  (LWRK1 .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient memory in CC_PQIM.'
           CALL QUIT('Insufficient memory in CC_PQIM.')
        END IF
C
C       calculate delta batch of backtransformed amplitudes:
C
        DTIME = SECOND()
        CALL CCG_TI(WORK(KTINT),ISYTIN,T2AM,ISYTAMP,C,ISYMC)
        TIMTI = TIMTI + SECOND() - DTIME
C
C       calculate Q intermediate:
C
        IOPT = 2
        DTIME = SECOND()
        CALL CC_ZWVI(WORK(KQINT),CTR2,ISYCTR,WORK(KTINT),
     &               ISYTIN,WORK(KEND1),LWRK1,IOPT)
        TIMZWV = TIMZWV + SECOND() - DTIME
C
        IF ( IOPTY.EQ.1 ) THEN
C
C         substract contributions from Y intermediate to i=j diagonal:
C
          ISYMA = ISYAIK
          KOFFY = IMATAB(ISYMC,ISYMA) + C
          CALL DCOPY(NVIR(ISYMA),YINT(KOFFY),NVIR(ISYMC),WORK(KYIM),1)

          DO ISYMI = 1, NSYM
           DO I = 1, NRHF(ISYMI)

            ISYMAI = MULD2H(ISYMA,ISYMI)
            KOFF   = KQINT + IT2BCD(ISYMAI,ISYMI) + IT1AM(ISYMA,ISYMI) 
     &             + NT1AM(ISYMAI)*(I-1) + NVIR(ISYMA)*(I-1) 

            CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KYIM),1,WORK(KOFF),1)

           END DO
          END DO
        END IF
C
C       scale here with 3 to balance factoer 1/3 in CCRHS_T2BT
C       and with 1/2 for compatibility with Asgers routines
C
        DTIME = SECOND()
        CALL DSCAL(NT2BCD(ISYAIK),THREE*HALF,WORK(KQINT),1)
        TIMSCL = TIMSCL + SECOND() - DTIME
C
C       write the intermediate to file:
C
        IADRW = IADRPQ(IVIRC,IV) + NT2BCD(ISYAIK)

        DTIME = SECOND()
        CALL  PUTWA2(LUPQIM,FILNAM,WORK(KQINT),IADRW,NT2BCD(ISYAIK))
        TIMIO = TIMIO + SECOND() - DTIME
        
      END DO
      END DO

*----------------------------------------------------------------------
* recover Zeta vector:
*----------------------------------------------------------------------
      DTIME = SECOND()
      CALL CCSD_T2TP(CTR2,WORK,LWORK,ISYCTR)
      CALL CCRHS_T2TR(CTR2,WORK,LWORK,ISYCTR)
      TIMZET = TIMZET + SECOND() - DTIME

*----------------------------------------------------------------------
* return:
*----------------------------------------------------------------------
      IF (IPRINT.GE.10) THEN
         TIMET = SECOND() - TIMET
         WRITE(LUPRI,*) '  Timings of CC_PQIMO:'
         WRITE(LUPRI,1) 'I/O             ', TIMIO
         WRITE(LUPRI,1) 'CC_TI           ', TIMTI
         WRITE(LUPRI,1) '2C-E of T2      ', TIMTAM
         WRITE(LUPRI,1) '2E+C of L2      ', TIMZET
         WRITE(LUPRI,1) 'CC_ZWV          ', TIMZWV
         WRITE(LUPRI,1) 'scaling etc.    ', TIMSCL
         WRITE(LUPRI,1) 'total CC_PQIMO: ', TIMET  
         CALL FLSHFO(LUPRI)
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'IADRPQ for vector:',IV
        DO ISYMC = 1, NSYM
        DO C = 1, NVIR(ISYMC)
          IVIRC = IVIR(ISYMC) + C
          WRITE (LUPRI,*) ISYMC,C,IVIRC, IADRPQ(IVIRC,IV)
        END DO
        END DO
      END IF

   1  FORMAT(1x,'Time used for',2x,A18,2x,': ',f10.2,' seconds')
 
      RETURN
      END
*======================================================================
*======================================================================
      SUBROUTINE CC_PQIAO(FILMO,LUMO,IADRMO,ISYPQMO,
     &                    FILAO,LUAO,IADRAO,IADRPQ,ITRAN,
     &                    XLAMD,ISYLAM,WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: read MO P and Q intermediates from file FILMO,
*              transform virtual index f back to AO basis,
*              and write the intermediates to file FILAO
*
*   IADRPQ contains on INPUT the start address of the P & Q intermed.
*   on file; on OUTPUT it will contain the start address for the
*   next set of P & Q intermediates.
*
*   Christof Haettig, September 1998
*
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <ccorb_dec.h>
#  include <maxorb_dec.h>
#  include <ccsdsym_dec.h>
#include <priunit_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>


      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*(*) FILAO, FILMO
      INTEGER ISYPQMO, ISYLAM, IADRPQ, ITRAN, LWORK, LUAO, LUMO
      INTEGER IADRAO(MAXORB,*), IADRMO(MAXORB,*)

#if defined (SYS_CRAY)
      REAL ZERO, ONE
      REAL XLAMD(*), WORK(*)
#else
      DOUBLE PRECISION ZERO, ONE
      DOUBLE PRECISION XLAMD(*), WORK(*)
#endif
      PARAMETER( ONE=1.0D0, ZERO=0.0D0 )

      INTEGER ISYMF,ISYEIJ,IVIRF,ISYDEL,NSCRM,INTERM,ISYRES,IDEL,IADR
      INTEGER KEIJF,KPQINT,KSCRM,KXLAM,KOFF,KEND1,LWRK1,NTOTEIJ
      CHARACTER*(1) ISTR(2)
      DATA ISTR /'P','Q'/
     
C     IF (LOCDBG) THEN
C       WRITE (LUPRI,*) 'Entered CC_PQIAO...'
C       CALL FLSHFO(LUPRI)
C       WRITE (LUPRI,*) 'work(1):',WORK(1)
C       CALL FLSHFO(LUPRI)
C       WRITE (LUPRI,*) 'work(lwork):',WORK(LWORK)
C       CALL FLSHFO(LUPRI)
C       WRITE (LUPRI,*) 'XLAMD(1):',XLAMD(1)
C       CALL FLSHFO(LUPRI)
C       WRITE (LUPRI,*) 'XLAMD(NLAMDT):',XLAMD(NLAMDT)
C       CALL FLSHFO(LUPRI)
C     END IF

*----------------------------------------------------------------------
* allocate work space for:
*       --  MO P or Q intermediate
*       --  one (eij) batch of a backtransformed intermediate
*       --  one vector of the transformation matrix
*----------------------------------------------------------------------
      NSCRM = 0
      DO ISYEIJ = 1, NSYM
        NSCRM = MAX(NSCRM,NT2BCD(ISYEIJ))
      END DO

      KPQINT = 1
      KSCRM  = KPQINT + NSCRM*NVIRT
      KXLAM  = KSCRM  + NSCRM
      KEND1  = KXLAM  + NVIRT
      LWRK1  = LWORK  - KEND1

      IF  (LWRK1 .LT. 0) THEN
         WRITE (LUPRI,*) 'Insufficient memory in CC_PQIAO.'
         CALL QUIT('Insufficient memory in CC_PQIAO.')
      END IF
     
*----------------------------------------------------------------------
* set output addresses for FILAO:
*----------------------------------------------------------------------
      ISYRES = MULD2H(ISYPQMO,ISYLAM)
      DO ISYDEL = 1, NSYM
         DO D = 1, NBAS(ISYDEL)
            IDEL = IBAS(ISYDEL) + D
            IADRAO(IDEL,ITRAN) = IADRPQ
            IADRPQ = IADRPQ + 2*NT2BCD(MULD2H(ISYRES,ISYDEL))
         END DO
      END DO

*----------------------------------------------------------------------
* loop over the two intermediates: INTERM = 1 -> P; INTERM = 2 -> Q
*----------------------------------------------------------------------
      DO INTERM = 1, 2
C
C        ---------------------------------
C        read in complete MO intermediate:
C        ---------------------------------
C
         KEIJF = KPQINT
         DO ISYMF = 1, NSYM
            ISYEIJ = MULD2H(ISYPQMO,ISYMF)
            DO F = 1, NVIR(ISYMF)
               IVIRF = IVIR(ISYMF) + F
               IADR  = IADRMO(IVIRF,ITRAN) + NT2BCD(ISYEIJ)*(INTERM-1)
               CALL GETWA2(LUMO,FILMO,WORK(KEIJF),IADR,NT2BCD(ISYEIJ))

               IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'CC_PQIAO> ',ISTR(INTERM),
     &                  ' interm. in MO for IVIRF=',IVIRF
                WRITE (LUPRI,'(5(F12.8))')
     &               (WORK(KEIJF+I),I=0,NT2BCD(ISYEIJ)-1)
               END IF

               KEIJF = KEIJF + NT2BCD(ISYEIJ)
            END DO
         END DO
C
C        --------------------------------------------------------------
C        transform to AO in a loop over delta index and write to FILAO:
C        --------------------------------------------------------------
C
         KEIJF = KPQINT
         DO ISYMF = 1, NSYM
            ISYDEL  = MULD2H(ISYMF,ISYLAM)
            ISYEIJ  = MULD2H(ISYPQMO,ISYMF)
            NTOTEIJ = MAX(NT2BCD(ISYEIJ),1)

            DO D = 1, NBAS(ISYDEL)

               IF ( NVIR(ISYMF).GT.0 ) THEN
                  KOFF = IGLMVI(ISYDEL,ISYMF) + D
                  CALL DCOPY(NVIR(ISYMF),XLAMD(KOFF),NBAS(ISYDEL),
     &                                   WORK(KXLAM),1)

                  CALL DGEMV('N',NT2BCD(ISYEIJ),NVIR(ISYMF),
     &                       ONE,WORK(KEIJF),NTOTEIJ,WORK(KXLAM),1,
     &                       ZERO,WORK(KSCRM),1)
               ELSE
                  CALL DZERO(WORK(KSCRM),NT2BCD(ISYEIJ))
               END IF

               IDEL = IBAS(ISYDEL) + D
               IADR = IADRAO(IDEL,ITRAN) + NT2BCD(ISYEIJ)*(INTERM-1)
               CALL PUTWA2(LUAO,FILAO,WORK(KSCRM),IADR,NT2BCD(ISYEIJ))

               IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'CC_PQIAO> ',ISTR(INTERM),
     &                  ' interm. in AO for IDEL=',IDEL,NVIR(ISYMF)
                WRITE (LUPRI,'(5(F12.8))')
     &               (WORK(KSCRM+I),I=0,NT2BCD(ISYEIJ)-1)
               END IF

            END DO

            KEIJF = KEIJF + NT2BCD(ISYEIJ)*NVIR(ISYMF)

         END DO

      END DO

*----------------------------------------------------------------------
* That's it. Return
*----------------------------------------------------------------------
 
      RETURN
      END
*======================================================================
C  /* Deck ccsd_z2tr */
      SUBROUTINE CCSD_Z2TR(Z2AM,ISYOPE)
C
C     Calculate in place two exchange plus coulomb of Z2 amplitudes.
C
C     Ch. Haettig, August 1998
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION Z2AM(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYOPE)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
            IF (ISYMI .GT. ISYMJ) GOTO 110
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMA = MULD2H(ISYMB,ISYMAB)
C
               IF (ISYMA .GT. ISYMB) GOTO 120
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) THEN
                     NRHFI =  J
                  ELSE
                     NRHFI = NRHF(ISYMI)
                  ENDIF
C
                  DO 140 I = 1,NRHFI
C
                     DO 150 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 160 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                         IF (ISYMAI.EQ.ISYMBJ) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + INDEX(NAI,NBJ)
                         ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI)*(NBJ-1)+NAI
                         ELSE IF (ISYMBJ.LT.ISYMAI) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                         END IF
C
                         IF (ISYMAJ.EQ.ISYMBI) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + INDEX(NAJ,NBI)
                         ELSE IF (ISYMAJ.LT.ISYMBI) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                         ELSE IF (ISYMBI.LT.ISYMAJ) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMBI)*(NAJ-1)+NBI
                         END IF
C
                           XAIBJ = Z2AM(NAIBJ) + TWO*Z2AM(NAJBI)
                           XAJBI = Z2AM(NAJBI) + TWO*Z2AM(NAIBJ)
C
                           Z2AM(NAIBJ) = XAIBJ
                           Z2AM(NAJBI) = XAJBI
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_zwvj */
      SUBROUTINE CC_ZWVJ(ZINT,CTR2,ISYMC2 ,TINT,ISYTIN,
     *                   WORK,LWORK,IOPT)
C
C     like Asgers CC_ZWVI routine, but CTR2 is used transposed
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ZINT(*), CTR2(*), TINT(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      ISYZIN = MULD2H(ISYMC2,ISYTIN)
C
      CALL DZERO(ZINT,NT2BCD(ISYZIN))
C
C-----------------------------------------------------
C     Transpose occupied indices of TINT if requested.
C-----------------------------------------------------
C
      IF (IOPT .EQ. 2) THEN
C
         CALL CCLT_P21I(TINT,ISYTIN,WORK,LWORK)
C
      ENDIF
C
C------------------------
C     Do the calculation.
C------------------------
C
      DO 100 ISYMDK = 1,NSYM
C
         ISYMEI = MULD2H(ISYMDK,ISYMC2)
         ISYMJ  = MULD2H(ISYMDK,ISYTIN)
C
         KOFF1  = IT2SQ(ISYMDK,ISYMEI) + 1
         KOFF2  = IT2BCD(ISYMDK,ISYMJ) + 1
         KOFF3  = IT2BCD(ISYMEI,ISYMJ) + 1
C
         NTOTEI = MAX(NT1AM(ISYMEI),1)
         NTOTDK = MAX(NT1AM(ISYMDK),1)
C
         CALL DGEMM('T','N',NT1AM(ISYMEI),NRHF(ISYMJ),NT1AM(ISYMDK),
     *              ONE,CTR2(KOFF1),NTOTDK,TINT(KOFF2),NTOTDK,ZERO,
     *              ZINT(KOFF3),NTOTEI)
C
  100 CONTINUE
C
C-------------------------------------------
C     Restore TINT intermediate if necessary
C-------------------------------------------
C
      IF (IOPT .EQ. 2) THEN
C
         CALL CCLT_P21I(TINT,ISYTIN,WORK,LWORK)
C
      ENDIF
C
      RETURN
      END
*======================================================================
      SUBROUTINE CC_PQI2(CTR2,ISYCTR,T2AM,ISYTAMP,YINT,IOPTY,
     &                   FILNAM,LUPQIM,IADRPQ,IADR,IV,WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: Calculate the P and Q intermediates from the 
*              Lagrangian multipiers CTR2
*              and the amplitude vector T2AM
*              and write them to a file direct access file FILNAM
*              
*   P{aik,del} = sum_{dl} (2 t(dl,k;del) - t(dk,l;del)) Zeta(dl;ai)
*   Q{aik,del} = sum_{dl} t(dl,k;del) (2Zeta(di;al) + Zeta(dl;ai))
*
*   for IOPTY = 1, the Y intermediate is substracted from the i=k
*   diagonal of the P and Q intermediates.
*
*   IADR contains on INPUT the start address of the P & Q intermed.
*   on file; on OUTPUT it will contain the start address for the
*   next set of P & Q intermediates.
*
*   Christof Haettig & Asger Halkier, August 1998
*
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <ccorb_dec.h>
#  include <maxorb_dec.h>
#  include <ccsdsym_dec.h>
#  include <iratdef_dec.h>
#include <priunit_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <iratdef.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER LWORK, LUPQIM, IOPTY
      CHARACTER*(*) FILNAM

#if defined (SYS_CRAY)
      REAL ONE, TWO, THREE, HALF, XNORM, DDOT, ZERO
      REAL WORK(*), CTR2(*), T2AM(*), YINT(*)
      REAL SECOND, TIMET, TIMIO, TIMZWV, TIMTAM, DTIME, TIMTI
      REAL TIMZET,TIMSCL
#else
      DOUBLE PRECISION ONE, TWO, THREE, HALF, XNORM, DDOT, ZERO
      DOUBLE PRECISION WORK(*), CTR2(*), T2AM(*), YINT(*)
      DOUBLE PRECISION SECOND,TIMET,TIMIO,TIMTI,TIMZWV,TIMTAM,DTIME
      DOUBLE PRECISION TIMZET,TIMSCL
#endif
      PARAMETER(ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,HALF=0.5D0,ZERO=0.0D0)

      INTEGER IV, IADR, IADRW
      INTEGER IADRPQ(MAXORB,IV)
      INTEGER ISYTAMP, ISYCTR, ISYMD, ISYTIN, ISYAIK, ISYMA, ISYMAI
      INTEGER KEND1, LWRK1, KEND2, LWRK2, NVIRD, NBASD
      INTEGER IERR, IOPT, IDEL, ILLL, KOFF1, KOFF2, KOFF, ISYMI
      INTEGER KTINT, KTJNT, KPINT, KQINT, KLAMDH, KLAMDP, KT1AM, KYIM
     

      CALL QENTER('CC_PQI2')
*----------------------------------------------------------------------
* set symmetries and allocate work space:
*----------------------------------------------------------------------
      KLAMDH = 1
      KLAMDP = KLAMDH + NLAMDT
      KT1AM  = KLAMDP + NLAMDT
      KEND1  = KT1AM  + NT1AMX
      LWRK1  = LWORK - KEND1

      IF ( LWRK1 .LT. 0 ) THEN
         WRITE (LUPRI,*) 'Insufficient memory in CC_PQI.'
         CALL QUIT('Insufficient memory in CC_PQI.')
      END IF

      TIMET  = SECOND()
      TIMIO  = ZERO
      TIMTI  = ZERO
      TIMZWV = ZERO
      TIMTAM = ZERO
      TIMZET = ZERO
      TIMSCL = ZERO

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'Norm of T2AMP:',
     &        DDOT(NT2AM(ISYTAMP),T2AM,1,T2AM,1)
         WRITE (LUPRI,*) 'T2AMP:'
         CALL CC_PRP(WORK,T2AM,ISYTAMP,0,1)
         WRITE (LUPRI,*) 'Norm of CTR2 :',
     &        DDOT(NT2SQ(ISYCTR),CTR2,1,CTR2,1)
         WRITE (LUPRI,*) 'CTR2:'
         CALL CC_PRSQ(WORK,CTR2,ISYCTR,0,1)
      END IF
*----------------------------------------------------------------------
* get XLAMD matrices from zero T1 amplitudes:
*----------------------------------------------------------------------
      CALL DZERO(WORK(KT1AM),NT1AMX)
       
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

*----------------------------------------------------------------------
* calculate the P intermediate in a loop over AO index delta:
*----------------------------------------------------------------------
      DO ISYMD = 1, NSYM
      DO ILLL = 1, NBAS(ISYMD)
        IDEL = IBAS(ISYMD) + ILLL
    
        ISYTIN = MULD2H(ISYMD,ISYTAMP)
        ISYAIK = MULD2H(ISYTIN,ISYCTR)
        ISYMA  = ISYAIK
    
        KTINT = KEND1
        KTJNT = KTINT + NT2BCD(ISYTIN)
        KPINT = KTJNT + NT2BCD(ISYTIN)
        KQINT = KPINT + NT2BCD(ISYAIK)
        KYIM  = KQINT + NT2BCD(ISYAIK)
        KEND2 = KYIM  + NVIR(ISYMA)
        LWRK2 = LWORK - KEND2

        IF  (LWRK2 .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient memory in CC_PQIM.'
           CALL QUIT('Insufficient memory in CC_PQIM.')
        END IF
C
C       calculate delta batch of backtransformed amplitudes:
C
        DTIME = SECOND()
        CALL CC_TI(WORK(KTINT),ISYTIN,T2AM,ISYTAMP,WORK(KLAMDH),1,
     &             WORK(KEND2),LWRK2,IDEL,ISYMD)
        TIMTI = TIMTI + SECOND() - DTIME
C
C       calculate 2 x Coul - Exc combination of backtransformed T:
C
        DTIME = SECOND()
        CALL DCOPY(NT2BCD(ISYTIN),WORK(KTINT),1,WORK(KTJNT),1)

        CALL DSCAL(NT2BCD(ISYTIN),-ONE,WORK(KTJNT),1)

        CALL CCLT_P21I(WORK(KTJNT),ISYTIN,WORK(KEND2),LWRK2)
        
        CALL DAXPY(NT2BCD(ISYTIN),TWO,WORK(KTINT),1,WORK(KTJNT),1)
        TIMTAM = TIMTAM + SECOND() - DTIME
C
C       calculate the P intermediate:
C
        IOPT = 1
        DTIME = SECOND()
        CALL CC_ZWVI(WORK(KPINT),CTR2,ISYCTR,WORK(KTJNT),
     &               ISYTIN,WORK(KEND2),LWRK2,IOPT)
        TIMZWV = TIMZWV + SECOND() - DTIME
C        
C       substract  contributions from Y intermediate from i=j diagonal:
C     
        IF (IOPTY.EQ.1) THEN
          KOFF1 = IMATAB(ISYMD,ISYMA) + 1
          KOFF2 = KLAMDH - 1 + IGLMVI(ISYMD,ISYMD) + IDEL - IBAS(ISYMD)
          NBASD = MAX(NBAS(ISYMD),1)
          NVIRD = MAX(NVIR(ISYMD),1)
        
          CALL DZERO(WORK(KYIM),NVIR(ISYMA))

          CALL DGEMV('T',NVIR(ISYMD),NVIR(ISYMA),ONE,YINT(KOFF1),
     &               NVIRD,WORK(KOFF2),NBASD,ZERO,WORK(KYIM),1)

          DO ISYMI = 1, NSYM
           DO I = 1, NRHF(ISYMI)

             ISYMAI = MULD2H(ISYMA,ISYMI)
             KOFF   = KPINT + IT2BCD(ISYMAI,ISYMI) + IT1AM(ISYMA,ISYMI)
     &              + NT1AM(ISYMAI)*(I-1) + NVIR(ISYMA)*(I-1)

             CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KYIM),1,WORK(KOFF),1)

           END DO
          END DO
        END IF
C        
C       scale with 1/2 for compatibility with Asgers routines
C
        DTIME = SECOND()
        CALL DSCAL(NT2BCD(ISYAIK),HALF,WORK(KPINT),1)
        TIMSCL = TIMSCL + SECOND() - DTIME
C
C       write the intermediate to file:
C
        IADRPQ(IDEL,IV) = IADR

        DTIME = SECOND()
        CALL  PUTWA2(LUPQIM,FILNAM,WORK(KPINT),IADR,NT2BCD(ISYAIK))
        TIMIO = TIMIO + SECOND() - DTIME

        IADR = IADR + 2*NT2BCD(ISYAIK)

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_PQI> P interm. in AO for IDEL = ',IDEL
           WRITE (LUPRI,'(5(F12.8))') 
     &          (WORK(KPINT+I),I=0,NT2BCD(ISYAIK)-1)
        END IF
        
      END DO
      END DO

*----------------------------------------------------------------------
*     Note that IADR contains now the start address for the next
*     set of P & Q intermediates. It must not be changed in the 
*     rest of this subroutine!
*----------------------------------------------------------------------

*----------------------------------------------------------------------
*     calculate (2 x Exc + Coul)/3 combination of ZETA:
*     (we interrupt here the loop over AO to calculate the modified
*      ZETA and accept that we recalculate the backtransformed 
*      amplitude since it the transformation of the amplitudes is
*      less expansive than the transformation and restruction of ZETA 
*      for each delta. both the transformation of TAMP and the 
*      transformation/restruction of ZETA scale with N V^2 O^2. )
*----------------------------------------------------------------------
      DTIME = SECOND()
      CALL CCRHS_T2BT(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      CALL CCSD_T2TP(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      TIMZET = TIMZET + SECOND() - DTIME

*----------------------------------------------------------------------
* calculate the Q intermediate in a loop over AO index delta:
*----------------------------------------------------------------------
      DO ISYMD = 1, NSYM
      DO ILLL = 1, NBAS(ISYMD)
        IDEL = IBAS(ISYMD) + ILLL
    
        ISYTIN = MULD2H(ISYMD,ISYTAMP)
        ISYAIK = MULD2H(ISYTIN,ISYCTR)
        ISYMA  = ISYAIK
    
        KTINT = KEND1
        KTJNT = KTINT + NT2BCD(ISYTIN)
        KPINT = KTJNT + NT2BCD(ISYTIN)
        KQINT = KPINT + NT2BCD(ISYAIK)
        KYIM  = KQINT + NT2BCD(ISYAIK)
        KEND2 = KYIM  + NVIR(ISYMA)
        LWRK2 = LWORK - KEND2

        IF  (LWRK2 .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient memory in CC_PQIM.'
           CALL QUIT('Insufficient memory in CC_PQIM.')
        END IF
C
C       calculate delta batch of backtransformed amplitudes:
C
        DTIME = SECOND()
        CALL CC_TI(WORK(KTINT),ISYTIN,T2AM,ISYTAMP,WORK(KLAMDH),1,
     &             WORK(KEND2),LWRK2,IDEL,ISYMD)
        TIMTI = TIMTI + SECOND() - DTIME
C
C       calculate Q intermediate:
C
        IOPT = 2
        DTIME = SECOND()
        CALL CC_ZWVI(WORK(KQINT),CTR2,ISYCTR,WORK(KTINT),
     &               ISYTIN,WORK(KEND2),LWRK2,IOPT)
        TIMZWV = TIMZWV + SECOND() - DTIME
C        
C       substract  contributions from Y intermediate from i=j diagonal:
C     
        IF (IOPTY.EQ.1) THEN
          KOFF1 = IMATAB(ISYMD,ISYMA) + 1
          KOFF2 = KLAMDH - 1 + IGLMVI(ISYMD,ISYMD) + IDEL - IBAS(ISYMD)
          NBASD = MAX(NBAS(ISYMD),1)
          NVIRD = MAX(NVIR(ISYMD),1)
        
          CALL DZERO(WORK(KYIM),NVIR(ISYMA))

          CALL DGEMV('T',NVIR(ISYMD),NVIR(ISYMA),ONE,YINT(KOFF1),
     &               NVIRD,WORK(KOFF2),NBASD,ZERO,WORK(KYIM),1)

          DO ISYMI = 1, NSYM
           DO I = 1, NRHF(ISYMI)

             ISYMAI = MULD2H(ISYMA,ISYMI)
             KOFF   = KQINT + IT2BCD(ISYMAI,ISYMI) + IT1AM(ISYMA,ISYMI)
     &              + NT1AM(ISYMAI)*(I-1) + NVIR(ISYMA)*(I-1)

             CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KYIM),1,WORK(KOFF),1)

           END DO
          END DO
        END IF
C        
C       scale here with 3 to balance factor 1/3 in CCRHS_T2BT
C       and with 1/2 for compatibility with Asgers routines
C
        DTIME = SECOND()
        CALL DSCAL(NT2BCD(ISYAIK),THREE*HALF,WORK(KQINT),1)
        TIMSCL = TIMSCL + SECOND() - DTIME
C
C       write the intermediate to file:
C
        IADRW = IADRPQ(IDEL,IV) + NT2BCD(ISYAIK)

        DTIME = SECOND()
        CALL  PUTWA2(LUPQIM,FILNAM,WORK(KQINT),IADRW,NT2BCD(ISYAIK))
        TIMIO = TIMIO + SECOND() - DTIME

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_PQI> Q interm. in AO for IDEL = ',IDEL
           WRITE (LUPRI,'(5(F12.8))')
     &          (WORK(KQINT+I),I=0,NT2BCD(ISYAIK)-1)
        END IF
        
      END DO
      END DO

*----------------------------------------------------------------------
* recover Zeta vector:
*----------------------------------------------------------------------
      DTIME = SECOND()
      CALL CCSD_T2TP(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      CALL CCRHS_T2TR(CTR2,WORK(KEND1),LWRK1,ISYCTR)
      TIMZET = TIMZET + SECOND() - DTIME

*----------------------------------------------------------------------
* return:
*----------------------------------------------------------------------
      IF (IPRINT.GE.10) THEN
         TIMET = SECOND() - TIMET
         WRITE(LUPRI,*) '  Timings of CC_PQIM:'
         WRITE(LUPRI,1) 'I/O             ', TIMIO
         WRITE(LUPRI,1) 'CC_TI           ', TIMTI
         WRITE(LUPRI,1) '2C-E of T2      ', TIMTAM
         WRITE(LUPRI,1) '2E+C of L2      ', TIMZET
         WRITE(LUPRI,1) 'CC_ZWV          ', TIMZWV
         WRITE(LUPRI,1) 'scaling etc.    ', TIMSCL
         WRITE(LUPRI,1) 'total CC_PQIM:  ', TIMET  
      END IF

   1  FORMAT(1x,'Time used for',2x,A18,2x,': ',f10.2,' seconds')
 
      CALL QEXIT('CC_PQI2')
      RETURN
      END
*======================================================================
C  /* Deck lamda2 */
      SUBROUTINE LAMDA2(XLAMDP,XLAMDH,ISYMX,T1AM,ISYTAM,
     &                  CMOP,CMOH,ISYCMO)
C
C     Calculate the lambda matrices.             asm 05-08-94
C     Extended for generalized symmetry.         ch  15-01-99
C     Extended for CMOP different from CMOH.     ch  01-11-99
C
C
#include <implicit.h>
      PARAMETER (ONE = 1.0D0)
      DIMENSION XLAMDH(*), XLAMDP(*), T1AM(*), CMOH(*), CMOP(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('LAMDA2')
C
      IF (ISYMX.NE.MULD2H(ISYCMO,ISYTAM)) THEN
         CALL QUIT('Symmetry mismatch in LAMDA2.')
      END IF
C
      DO ISYALP = 1,NSYM
C
         ISYMA = MULD2H(ISYCMO,ISYALP)
         ISYMI = MULD2H(ISYTAM,ISYMA)
C
         NBASA = MAX(NBAS(ISYALP),1)
         NVIRA = MAX(NVIR(ISYMA),1)
C
         KOFF1 = IGLMVI(ISYALP,ISYMA) + 1
         KOFF2 = IT1AM(ISYMA,ISYMI)   + 1
         KOFF3 = IGLMRH(ISYALP,ISYMI) + 1
C
         CALL DGEMM('N','N',NBAS(ISYALP),NRHF(ISYMI),NVIR(ISYMA),
     *               ONE,CMOH(KOFF1),NBASA,T1AM(KOFF2),NVIRA,
     *               ONE,XLAMDH(KOFF3),NBASA)
      END DO
C
      DO ISYALP = 1,NSYM
C
         ISYMI = MULD2H(ISYCMO,ISYALP)
         ISYMA = MULD2H(ISYTAM,ISYMI)
C
         NBASA = MAX(NBAS(ISYALP),1)
         NVIRA = MAX(NVIR(ISYMA),1)
C
         KOFF4 = IGLMRH(ISYALP,ISYMI) + 1
         KOFF5 = IT1AM(ISYMA,ISYMI)  + 1
         KOFF6 = IGLMVI(ISYALP,ISYMA) + 1
C
         CALL DGEMM('N','T',NBAS(ISYALP),NVIR(ISYMA),NRHF(ISYMI),
     *              -ONE,CMOP(KOFF4),NBASA,T1AM(KOFF5),NVIRA,
     *               ONE,XLAMDP(KOFF6),NBASA)
C
      END DO
C
      CALL QEXIT('LAMDA2')
      RETURN
      END
C  /* Deck cmo_reorder2 */
      SUBROUTINE CMO_REORDER2(CMO,XLAMDA,ISYCMO)
C
C     Reorder the symmetry ordering of the MO coefficient matrix.
C     First occupied orbitals in different symmetries and then
C     the virtuals in different symmetries.
C
C     Henrik Koch and Alfredo Sanchez.       30-Jun-1994
C
C     revised and extended for general symmetry of CMO 
C     by Christof Haettig, Jan '99
C
#include <implicit.h>
      DIMENSION CMO(*), XLAMDA(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER NCMO(8), ICMO(8,8)
C
      DO ISYM = 1, NSYM
         ICOUNT = 0
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            ICMO(ISYM1,ISYM2) = ICOUNT
            ICOUNT = ICOUNT + NBAS(ISYM1)*NORBS(ISYM2)
         END DO
         NCMO(ISYM) = ICOUNT
      END DO
C
C----------------------
C     Reorder orbitals: 
C----------------------
C
      DO ISYMO = 1,NSYM
C
        ISYAO = MULD2H(ISYMO,ISYCMO)
C
        KOFF1  = ICMO(ISYAO,ISYMO) + NBAS(ISYAO)*NRHFFR(ISYMO) + 1
        KOFF2  = IGLMRH(ISYAO,ISYMO) + 1
        CALL DCOPY(NBAS(ISYAO)*NRHF(ISYMO),CMO(KOFF1),1,XLAMDA(KOFF2),1)
C
        KOFF1  = KOFF1 + NBAS(ISYAO)*NRHFS(ISYMO)
        KOFF2  = IGLMVI(ISYAO,ISYMO) + 1
        CALL DCOPY(NBAS(ISYAO)*NVIR(ISYMO),CMO(KOFF1),1,XLAMDA(KOFF2),1)
C
      END DO
C
C----------------------
C     Print if desired:
C----------------------
C
      IF (LOCDBG) THEN
         CALL AROUND('MO-coefficient matrix in CMO_REORDER2')
         DO ISYMO = 1,NSYM
            ISYAO = MULD2H(ISYMO,ISYCMO)
            WRITE(LUPRI,'(/,/,7X,"Symmetry numbers :",2I5)') ISYAO,ISYMO
            IF (NORBS(ISYMO) .NE. 0) THEN
               WRITE(LUPRI,'(/,/,7X,"Sirius ordering")')
               WRITE(LUPRI,'(7X,"---------------")')
               KOFF1 = ICMO(ISYAO,ISYMO) + 1
               CALL OUTPUT(CMO(KOFF1),1,NBAS(ISYAO),1,NORBS(ISYMO),
     *                     NBAS(ISYAO),NORBS(ISYMO),1,LUPRI)
               WRITE(LUPRI,'(/,/,7X,"CCSD ordering occupied part")')
               WRITE(LUPRI,'(7X,"---------------------------")')
               KOFF2 = IGLMRH(ISYAO,ISYMO) + 1
               CALL OUTPUT(XLAMDA(KOFF2),1,NBAS(ISYAO),1,NRHF(ISYMO),
     *                     NBAS(ISYAO),NRHF(ISYMO),1,LUPRI)
               WRITE(LUPRI,'(/,/,7X,"CCSD ordering virtual part")')
               WRITE(LUPRI,'(7X,"--------------------------")')
               KOFF3 = IGLMVI(ISYAO,ISYMO) + 1
               CALL OUTPUT(XLAMDA(KOFF3),1,NBAS(ISYAO),1,NVIR(ISYMO),
     *                     NBAS(ISYAO),NVIR(ISYMO),1,LUPRI)
            ELSE
               WRITE(LUPRI,'(/,/,7X,"This symmetry is empty")')
            ENDIF
         END DO   
      END IF
C
      RETURN
C
C
      END
C  /* Deck cc_aodens2 */
      SUBROUTINE CC_AODENS2(XLAMDP0,XLAMDH0,ISYLM0,
     &                      XLAMDP1,XLAMDH1,ISYLM1,
     &                      DENS,ICORE,IOPT,WORK,LWORK)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Calculate the AO-density matrix used in contructing
C     the AO Fock matrix.
C
C     IOPT = 0 : vector function: XLAMDP0 x XLAMDH0
C     IOPT = 1 : response:        XLAMDP0 x XLAMDH1
C     IOPT = 2 : derivatives:     XLAMDP0 x XLAMDH1 + XLAMDP1 x XLAMDH0
C
C     Henrik Koch and Alfredo Sanchez.       11-Jul-1994
C     Christof Haettig                       22-Jan-1999
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <dummy.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XLAMDP0(*), XLAMDH0(*), XLAMDP1(*), XLAMDH1(*)
      DIMENSION DENS(*), WORK(LWORK)
#include <inftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      IF (IOPT.GE.1) THEN
         ISYDEN = MULD2H(ISYLM0,ISYLM1)
         CALL DZERO(DENS,N2BST(ISYDEN))
      END IF
C
      DO ISYMK = 1,NSYM
C
         IF (IOPT.EQ.0) THEN
C
            ISYMA = MULD2H(ISYLM0,ISYMK)      ! XLAMDP0
            ISYMB = MULD2H(ISYLM0,ISYMK)      ! XLAMDH0
C
            KOFF1 = 1 + IGLMRH(ISYMA,ISYMK)   ! XLAMDP0
            KOFF2 = 1 + IGLMRH(ISYMB,ISYMK)   ! XLAMDH0
            KOFF3 = 1 + IAODIS(ISYMA,ISYMB)
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
C
            CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYMK),ONE,
     *                 XLAMDP0(KOFF1),NBASA,XLAMDH0(KOFF2),NBASB,ZERO,
     *                 DENS(KOFF3),NBASA)
C
         END IF
C
         IF (IOPT.GE.1) THEN
C
            ISYMA = MULD2H(ISYLM0,ISYMK)      ! XLAMDP0
            ISYMB = MULD2H(ISYLM1,ISYMK)      ! XLAMDH1
C
            KOFF1 = 1 + IGLMRH(ISYMA,ISYMK)   ! XLAMDP0
            KOFF2 = 1 + IGLMRH(ISYMB,ISYMK)   ! XLAMDH1
            KOFF3 = 1 + IAODIS(ISYMA,ISYMB)
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
C
            CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYMK),ONE,
     *                 XLAMDP0(KOFF1),NBASA,XLAMDH1(KOFF2),NBASB,ONE,
     *                 DENS(KOFF3),NBASA)
C
         END IF
C
         IF (IOPT.EQ.2) THEN
C
            ISYMA = MULD2H(ISYLM1,ISYMK)         ! XLAMDP1
            ISYMB = MULD2H(ISYLM0,ISYMK)         ! XLAMDH0
C
            KOFF1 = 1 + IGLMRH(ISYMA,ISYMK)      ! XLAMDP1
            KOFF2 = 1 + IGLMRH(ISYMB,ISYMK)      ! XLAMDH0
            KOFF3 = 1 + IAODIS(ISYMA,ISYMB)
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
C
            CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYMK),ONE,
     *                 XLAMDP1(KOFF1),NBASA,XLAMDH0(KOFF2),NBASB,ONE,
     *                 DENS(KOFF3),NBASA)
C
         END IF
C
      END DO
C
C
C-----------------------------
C     Include frozen orbitals.
C-----------------------------
C
      IF ( (FROIMP.OR.FROEXP) .AND. (ICORE .EQ. 1) ) THEN
C
         IF (IOPT.NE.0) THEN
           WRITE (LUPRI,*) 'CC_AODENS2: ICORE=1 not yet available '//
     &           'for IOPT>0.'
           CALL QUIT('CC_AODENS2: ICORE=1 not yet '//
     &          'available for IOPT>0.')
         END IF
C
         IF (LWORK .LT. NLAMDS) THEN
            CALL QUIT('Insufficient space in CCSD_AODENS')
         ENDIF
C
C-------------------------------------------------
C        Read MO-coefficients from interface file.
C-------------------------------------------------
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
C
         CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
         READ (LUSIFC)
C
         READ (LUSIFC)
         READ (LUSIFC) (WORK(I), I=1,NLAMDS)
C
         CALL GPCLOSE(LUSIFC,'KEEP')
C
C-------------------------------------------------------
C        Add contribution from frozen occupied orbitals.
C-------------------------------------------------------
C
         KOFF1 = 0
         KOFF2 = 0
         DO ISYMK = 1,NSYM
C
            ISYMA = MULD2H(ISYLM0,ISYMK)
            ISYMB = MULD2H(ISYLM0,ISYMK)
C
            DO II = 1,NRHFFR(ISYMK)
C
               K = KFRRHF(II,ISYMK)
C
               DO B = 1,NBAS(ISYMB)
                  DO A = 1,NBAS(ISYMA)
C
                     NAK = KOFF1 + NBAS(ISYMA)*(K - 1) + A
                     NBK = KOFF1 + NBAS(ISYMB)*(K - 1) + B
                     NAB = KOFF2 + NBAS(ISYMA)*(B - 1) + A
C
                     DENS(NAB) = DENS(NAB) + WORK(NAK)*WORK(NBK)
C
                  END DO
               END DO
C
            END DO
C
            KOFF1 = KOFF1 + NBAS(ISYMK)*NORBS(ISYMK)
            KOFF2 = KOFF2 + NBAS(ISYMA)*NBAS(ISYMB)
C
         END DO
C
      ENDIF
C
      END
*=====================================================================*
      SUBROUTINE CC_CDBAR2(CDBAR,EMAT1,FAC,ISYCDBAR)
*---------------------------------------------------------------------*
*
* Purpose: add EMAT1 (scaled by FAC) to the k=i diagonal of the 
*          CDBAR intermediate
*
* C. Haettig August 1998
*---------------------------------------------------------------------*
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION CDBAR(*),EMAT1(*)
C
      DO ISYMAI = 1, NSYM
        ISYMCI = MULD2H(ISYCDBAR,ISYMAI)

        DO ISYMI = 1, NSYM
          ISYMA = MULD2H(ISYMAI,ISYMI)
          ISYMC = MULD2H(ISYMCI,ISYMI)

          DO I = 1, NRHF(ISYMI)
            DO A = 1, NVIR(ISYMA)

              NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
              KOFF = IT2SQ(ISYMCI,ISYMAI) + NT1AM(ISYMCI)*(NAI-1)

              DO C = 1, NVIR(ISYMC)
                NCI = IT1AM(ISYMC,ISYMI)  + NVIR(ISYMC)*(I-1) + C
                NAC = IMATAB(ISYMA,ISYMC) + NVIR(ISYMA)*(C-1) + A
                CDBAR(KOFF+NCI) = CDBAR(KOFF+NCI) + FAC*EMAT1(NAC)
              END DO

            END DO
          END DO
             
        END DO
      END DO
 
      RETURN
      END
*=====================================================================*
*                   END OF SUBROUTINE CC_CDBAR2                       *
*=====================================================================*
C  /* Deck cc_t2mo3 */
      SUBROUTINE CC_T2MO3(RHO1,CTR2,ISYMC2,OMEGA2,
     *                    RHO2,GAMMA,GIMB,GIMC,
     *                    XLAMPB,ISYMPB, XLAMPC,ISYMPC,
     *                    WORK,LWORK,ISYMBF,
     *                    ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
C
C     Henrik Koch and Alfredo Sanchez.       15-July-1994
C
C     Transform the Omega2 vector from the AO basis to the MO
C     basis.
C
C     Ove Christiansen 4-8-1995, Christof Haettig 22-1-1997
C     G-intermediates introduced, Christof Haettig 19-6-1998
C     relaxation contribution to GAMMA intermediate, CH 20-6-1998
C
C     N.B.: this version thus NOT initialize RHO2 & GAMMA !!
C
C     Generalizations for CC response:
C
C      1.ISYMBF is the symmetry of the BF (ali,bej) vector.
C
C      2.Transformation to virtual space (BF intermediate)
C        with one or two non-total symmetric lambda matrices
C        (XLAMPC and XLAMPB(formerly XLAMDP)).
C
C      3.Transformation to occupied space (Gamma intermediate)
C        with total symmetric lambda matrix (XLAMPB) only!
C
C      note that if LGAMMA is true gamma is the gamma vector on return
C      with the same symmetry as the input BF. 
C      if LGAMMA is false the gamma intermediate is not returned.
C      LGAMMA=.TRUE. requires XLAMPB total symmetric.
C
C      ICON is 2 for response to calculate a-tild,ibj and ai,b-tilde,j
C           for the transformation to the virtual space. 
C
C      ICON=4 as for ICON=2 but include extra relaxation contribution
C             for GAMMA intermediate
C
C      NB these changes are only carried through completely and
C      tested for omegor
C
C     Asger Halkier 2/11-1995:
C
C      For ICON equal to 3 the contraction of the (ali,bej) vector with
C      the trialvector CTR2 (i.e the LT21BF-term) is calculated and
C      stored in RHO1!
C      NOTE: ICON=3 requires a total symmetric XLAMPB, generalization
C            for non-total symmetric XLAMPB is not carried through
C
C     Ove Christiansen 4-10-1996:
C 
C      For use in F-matrix generalize ICON .EQ. 3 section
C
C     Christof Haettig 19-6-1998: IOPTG & LO3BF flags
C
C      if IOPTG = 1, the BF-intermediate is contracted with
C      XLAMPB to the G-intermediate: 
C         GB_{alp i} = sum_{bet j} XLAMPB_{bet j}  BF_{alp i bet j} 
C      for ICON=4 calculate also:
C         GC_{alp i} = sum_{bet j} XLAMPC_{bet j}  BF_{alp i bet j} 
C
C      if IOPTG = 2, the BF-intermediate is contracted with
C      XLAMPB to the G-intermediate: 
C         GB_{alpha i} = sum_{beta j} XLAMPB_{beta j} 
C                 * ( 2 BF_{alpha i beta j} - BF_{beta i alpha j} )
C      for ICON=4 calculate also:
C         GC_{alpha i} = sum_{beta j} XLAMPC_{beta j} 
C                 * ( 2 BF_{alpha i beta j} - BF_{beta i alpha j} )
C
C      if LO3BF = .TRUE., it is assumed that the beta index of the
C         BF-intermediate is already transformed to occupied
C
C     Christof Haettig 19-10-1998: 
C      LT21BF section adapted for general symmetry for XLAMPB  
C         RHO1(ck) = sum_dij CTR2(dicj) BF(dikj)
C      with
C         BF(dikj) = XLAMPC(alpha,d) XLAMDPB(beta,k) BF(alpha i beta j)
C      for ICON=6 include extra relaxation contribution in BF(dikj):
C         BF(dikj) = XLAMPC(alpha,d) XLAMDPB(beta,k) BF(alpha i beta j)
C                  + XLAMPB(alpha,d) XLAMDPC(beta,k) BF(alpha i beta j)
C
C     NOTE: Linear response options only valid and debugged for OMEGOR!
C
#include <implicit.h>
#include <maxorb.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*), CTR2(*), OMEGA2(*), RHO2(*), GAMMA(*),
     *          XLAMPC(*), XLAMPB(*), WORK(*), GIMB(*), GIMC(*)
      LOGICAL LGAMMA, LO3BF, LBFZETA
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
#include <cclr.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMO2 = MULD2H(ISYMBF,MULD2H(ISYMPC,ISYMPB))
      ISYMO1 = MULD2H(ISYMO2,ISYMC2)

      DO 100 ISYMJ = 1,NSYM
         DO 110 ISYMI = 1,NSYM
C
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
            ISALBE = MULD2H(ISYMIJ,ISYMBF)
C
            DO 120 ISYBE = 1,NSYM
C
               ISYAL  = MULD2H(ISYBE,ISALBE)
               ISYALI = MULD2H(ISYAL,ISYMI)
               ISYBEJ = MULD2H(ISYBE,ISYMJ)

C
C-----------------------------------------------
C              Dynamic allocation of work space.
C-----------------------------------------------
C
               NVA = MAX(NVIR(MULD2H(ISYAL,ISYMPC)),
     &                   NVIR(MULD2H(ISYAL,ISYMPB)) )
               NRA = MAX(NRHF(MULD2H(ISYAL,ISYMPC)),
     &                   NRHF(MULD2H(ISYAL,ISYMPB)) )
               NVB = MAX(NVIR(MULD2H(ISYBE,ISYMPC)),
     &                   NVIR(MULD2H(ISYBE,ISYMPB)) )
               NRB = MAX(NRHF(MULD2H(ISYBE,ISYMPC)),
     &                   NRHF(MULD2H(ISYBE,ISYMPB)) )
               NVA = MAX(NVA,NRA)
               NVB = MAX(NVB,NRB)
C
               KSCR1 = 1
               KSCR2 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)
               KSCR3 = KSCR2 + NBAS(ISYAL)*NVB
               IF (LGAMMA) THEN
                  KSCR4 = KSCR3 + NVA*NVB
                  KSCR5 = KSCR4 + NBAS(ISYAL)*NRB
                  KEND1 = KSCR5 + NRA*NRB
               ELSE
                  KEND1 = KSCR3 + NVA*NVB
               END IF
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space in CC_T2MO3')
               END IF
C
               DO 130 J = 1,NRHF(ISYMJ)
                  DO 140 I = 1,NRHF(ISYMI)
C
C------------------------------------------
C                    Squareup the AB block.
C------------------------------------------
C
                     IF (((.NOT. OMEGSQ) .AND. (.NOT. OMEGOR)) 
     *                    .OR. LBFZETA ) THEN
C
                     DO 150 B = 1,NBAS(ISYBE)
                        NBJ   = IT1AO(ISYBE,ISYMJ)
     *                        + NBAS(ISYBE)*(J-1) + B
                        DO 155 A = 1,NBAS(ISYAL)
C
                           NAI   = IT1AO(ISYAL,ISYMI)
     *                           + NBAS(ISYAL)*(I-1) + A
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           IF (ISYALI .EQ. ISYBEJ) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYALI .LT. ISYBEJ) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + NT1AO(ISYALI)*(NBJ - 1) + NAI
                           ELSEIF (ISYALI .GT. ISYBEJ) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + NT1AO(ISYBEJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           WORK(NAB) = OMEGA2(NAIBJ)
C
  155                   CONTINUE
  150                CONTINUE
C
                     ENDIF
C
                     IF (OMEGSQ .AND. (.NOT.LBFZETA) ) THEN
C
                     DO 160 B = 1,NBAS(ISYBE)
                        NBJ   = IT1AO(ISYBE,ISYMJ)
     *                        + NBAS(ISYBE)*(J-1) + B
                        DO 165 A = 1,NBAS(ISYAL)
C
                           NAI   = IT1AO(ISYAL,ISYMI)
     *                           + NBAS(ISYAL)*(I-1) + A
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           NAIBJ = IT2AOS(ISYALI,ISYBEJ)
     *                           + NT1AO(ISYALI)*(NBJ - 1) + NAI
                           NBJAI = IT2AOS(ISYBEJ,ISYALI)
     *                           + NT1AO(ISYBEJ)*(NAI - 1) + NBJ
C
                           WORK(NAB) = OMEGA2(NAIBJ) + OMEGA2(NBJAI)
C
  165                   CONTINUE
  160                CONTINUE
C
                     ENDIF
C
                     IF (OMEGOR .AND. (.NOT.LBFZETA) ) THEN
C
                     IF (LO3BF) THEN
C
C                      ISYML   = MULD2H(ISYBE,ISYMPB)
C
                       ISYAIJ  = MULD2H(ISYAL,ISYMIJ)
                       ISYML   = MULD2H(ISYAIJ,ISYMBF)
                       ISYLAL  = MULD2H(ISYAL,ISYML)

                       NIJ = IMATIJ(ISYMI,ISYMJ) 
     *                     + NRHF(ISYMI)*(J-1) + I
C
                       NALIJP = IT2AOIJ(ISYLAL,ISYMIJ)
     *                        + NT1AO(ISYLAL)*(NIJ - 1) 
     *                        + IT1AO(ISYAL,ISYML) + 1
C
                       CALL DCOPY(NRHF(ISYML)*NBAS(ISYAL),
     *                            OMEGA2(NALIJP),1,WORK(KSCR1),1)
C
                     ELSE
C
                       IF (ISYMI .EQ. ISYMJ) THEN
                          NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                          FAC1 = ONE
                          IF (I .GT. J) FAC1 = -ONE
                       ELSE IF (ISYMI .LT. ISYMJ) THEN
                          NIJ = IMIJP(ISYMI,ISYMJ)
     *                        + NRHF(ISYMI)*(J - 1) + I
                          FAC1 = ONE
                       ELSE
                          NIJ = IMIJP(ISYMI,ISYMJ)
     *                        + NRHF(ISYMJ)*(I - 1) + J
                          FAC1 = -ONE
                       ENDIF
C
                       DO 166 B = 1,NBAS(ISYBE)
                          DO 167 A = 1,NBAS(ISYAL)
C
                            IF (ISYAL .EQ. ISYBE) THEN
                               NABP = IAODPK(ISYAL,ISYBE)
     *                              + INDEX(A,B)
                               FAC2 = ONE
                               IF (A .GT. B) FAC2 = -ONE
                            ELSE IF (ISYAL .LT. ISYBE) THEN
                               NABP = IAODPK(ISYAL,ISYBE)
     *                              + NBAS(ISYAL)*(B - 1) + A
                               FAC2 = ONE
                            ELSE
                               NABP = IAODPK(ISYAL,ISYBE)
     *                              + NBAS(ISYBE)*(A - 1) + B
                               FAC2 = -ONE
                            ENDIF
C
                            NABIJP = IT2ORT(ISALBE,ISYMIJ)
     *                             + NNBST(ISALBE)*(NIJ - 1) + NABP
C
                            NABIJM = NT2ORT(ISYMBF)
     *                             + IT2ORT(ISALBE,ISYMIJ)
     *                             + NNBST(ISALBE)*(NIJ - 1) + NABP
C
                            NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                            FAC = FAC1*FAC2
C
                            WORK(NAB) =
     *                        HALF*(OMEGA2(NABIJP) + FAC*OMEGA2(NABIJM))
C
  167                     CONTINUE
  166                  CONTINUE
C
                     ENDIF
C
                     ENDIF
C
C------------------------------------------------------------
C                    Transform the AB block to virtual space.
C------------------------------------------------------------
C
                  IF (ICON.NE.3 .AND. ICON.NE.6 .AND. (.NOT.LO3BF)) THEN
C
                     ISYMA  = MULD2H(ISYAL,ISYMPC)
                     ISYMB  = MULD2H(ISYBE,ISYMPB)
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NVIRA = MAX(NVIR(ISYMA),1)
C
                     KOFF1 = IGLMVI(ISYBE,ISYMB) + 1
                     CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMPB(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPC(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
C--------------------------------------------
C                    Store in omega2 vector.
C--------------------------------------------
C
                     DO 170 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 180 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
                           NAB   = KSCR3 + NVIR(ISYMA)*(B - 1) + A - 1
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
C
                              IF (NAI .GT. NBJ) GOTO 180
C
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              GOTO 180
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
  180                   CONTINUE
  170                CONTINUE

C
                  ENDIF
C
                  IF (ICON.EQ.2 .AND. LO3BF) THEN
C
C------------------------------------------------------------
C                    Transform the AL block to virtual space
C                    and store in Omega2 vector
C------------------------------------------------------------
C
                     ISYMA  = MULD2H(ISYAL,ISYMPB)
                     ISYMB  = MULD2H(ISYML,ISYMPC)
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NVIRB = MAX(NVIR(ISYMB),1)
                     NVIRA = MAX(NVIR(ISYMA),1)
C
                     KOFF1 = IT1AM(ISYMB,ISYML) + 1
                     CALL DGEMM('N','T',NBAS(ISYAL),NVIR(ISYMB),
     *                          NRHF(ISYML),-ONE,WORK(KSCR1),NBASA,
     *                          XLAMPC(KOFF1),NVIRB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPB(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
                     DO B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
C
                        IF (ISYMAI.EQ.ISYMBJ .AND. ISYMI.EQ.ISYMJ
     *                                       .AND. I.EQ.J) THEN
                           NBB = KSCR3 + NVIR(ISYMA)*(B - 1) + B - 1
                           WORK(NBB) = TWO * WORK(NBB)
                        END IF
C
                        DO A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
                           NAB   = KSCR3 + NVIR(ISYMA)*(B - 1) + A - 1
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
                        END DO
                     END DO
C
                  ENDIF
C
C--------------------------------------
C                 CCLR contribution.
C--------------------------------------
C
                  IF ((ICON.EQ.2 .OR. ICON.EQ.4).AND.(.NOT.LO3BF)) THEN
C
                     CALL DZERO(WORK(KSCR2),NVA*NVB)
                     ISYMA = MULD2H(ISYAL,ISYMPB)
                     ISYMB = MULD2H(ISYBE,ISYMPC)
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NVIRA = MAX(NVIR(ISYMA),1)
C
                     KOFF1 = IGLMVI(ISYBE,ISYMB) + 1
                     CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMPC(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPB(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
C--------------------------------------------
C                    Store in omega2 vector.
C--------------------------------------------
C
                     DO 181 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 182 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              IF (NAI .GT. NBJ ) GOTO 182
C                             NAIBJ = IT2AM(ISYALI,ISYBEJ)
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              GOTO 182
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           NAB  = KSCR3+ NVIR(ISYMA)*(B - 1) + A - 1
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
  182                   CONTINUE
  181                CONTINUE
C

                  ENDIF

C
C============================================================
C                 calculating G-intermediate(s)
C============================================================
C
                  IF (IOPTG.GT.0) THEN

                    FACG = ONE
                    IF (IOPTG.EQ.2) FACG = TWO

C                   ---------------------------------------------------
C                   G_{alp i} += F sum_bet LAM_{bet j} BF_{alp i bet j}
C                   G_{alp i} += F BF_{alp i j j}
C                   ---------------------------------------------------
                    IF (MULD2H(ISYBE,ISYMJ).EQ.ISYMPB 
     *                                .AND. .NOT.LO3BF) THEN


                      KLAMD = IGLMRH(ISYBE,ISYMJ) + NBAS(ISYBE)*(J-1)+1
                      KGIM  = IT1AO(ISYAL,ISYMI)  + NBAS(ISYAL)*(I-1)+1
                      NBASA = MAX(NBAS(ISYAL),1)

                      CALL DGEMV('N',NBAS(ISYAL),NBAS(ISYBE),
     *                          FACG,WORK(KSCR1),NBASA,XLAMPB(KLAMD),1,
     *                          ONE,GIMB(KGIM),1)
                        
                    END IF
                      
                    IF (ISYML.EQ.ISYMJ .AND. LO3BF) THEN

                      KGIM  = IT1AO(ISYAL,ISYMI) + NBAS(ISYAL)*(I-1)+1
                      KOFF  = KSCR1 + NBAS(ISYAL)*(J-1)

                      CALL DAXPY(NBAS(ISYAL),FACG,WORK(KOFF),1,
     *                                            GIMB(KGIM),1)

                    END IF

                  END IF

                  IF (IOPTG.EQ.2) THEN

C                   ---------------------------------------------------
C                   G_{bet i} += - sum_alp LAM_{alp j} BF_{alp i bet j}
C                   G_{alp j} += - BF_{alp i i j}
C                   ---------------------------------------------------
c... unused version from F matrix before merge into CC_XIETA version
c                   IF (MULD2H(ISYAL,ISYMJ).EQ.ISYMPB
c    *                                .AND. .NOT.LO3BF) THEN
c
c                     KLAMD = IGLMRH(ISYAL,ISYMJ) + NBAS(ISYAL)*(J-1)+1
c                     KGIM  = IT1AO(ISYBE,ISYMI)  + NBAS(ISYBE)*(I-1)+1
c                     NBASA = MAX(NBAS(ISYAL),1)
c
c                     CALL DGEMV('N',NBAS(ISYAL),NBAS(ISYBE),
c    *                           -ONE,WORK(KSCR1),NBASA,XLAMPB(KLAMD),1,
c    *                           ONE,GIMB(KGIM),1)
c                   END IF

C                   ---------------------------------------------------
C                   G_{alp j} += - sum_bet LAM_{bet i} BF_{alp i bet j}
C                   ---------------------------------------------------
                    IF (MULD2H(ISYBE,ISYMI).EQ.ISYMPB
     *                                .AND. .NOT.LO3BF) THEN

                      KLAMD = IGLMRH(ISYBE,ISYMI) + NBAS(ISYBE)*(I-1)+1
                      KGIM  = IT1AO(ISYAL,ISYMJ)  + NBAS(ISYAL)*(J-1)+1
                      NBASA = MAX(NBAS(ISYAL),1)

                      CALL DGEMV('N',NBAS(ISYAL),NBAS(ISYBE),
     *                          -ONE,WORK(KSCR1),NBASA,XLAMPB(KLAMD),1,
     *                           ONE,GIMB(KGIM),1)
                    END IF

                    IF (ISYML.EQ.ISYMI .AND. LO3BF) THEN
                      
                      KGIM = IT1AO(ISYAL,ISYMJ) + NBAS(ISYAL)*(J-1)+1
                      KOFF = KSCR1 + NBAS(ISYAL)*(I-1)

                      CALL DAXPY(NBAS(ISYAL),-ONE,WORK(KOFF),1,
     *                                            GIMB(KGIM),1)

                    END IF

                  END IF
C
                  IF (ICON.EQ.4 .AND. IOPTG.GT.0) THEN

                    IF (LO3BF) CALL QUIT('Illegal option in CC_T2MO3.')

                    FACG = ONE
                    IF (IOPTG.EQ.2) FACG = TWO

C                   ---------------------------------------------------
C                   G_{alp i} += F sum_bet LAM_{bet j} BF_{alp i bet j}
C                   ---------------------------------------------------
                    IF (MULD2H(ISYBE,ISYMJ).EQ.ISYMPC) THEN

                      KLAMD = IGLMRH(ISYBE,ISYMJ) + NBAS(ISYBE)*(J-1)+1
                      KGIM  = IT1AO(ISYAL,ISYMI)  + NBAS(ISYAL)*(I-1)+1
                      NBASA = MAX(NBAS(ISYAL),1)

                      CALL DGEMV('N',NBAS(ISYAL),NBAS(ISYBE),
     *                          FACG,WORK(KSCR1),NBASA,XLAMPC(KLAMD),1,
     *                          ONE,GIMC(KGIM),1)

                    END IF

                  END IF

                  IF (ICON.EQ.4 .AND. IOPTG.EQ.2) THEN

                    IF (LO3BF) CALL QUIT('Illegal option in CC_T2MO3.')

C                   ---------------------------------------------------
C                   G_{bet i} += - sum_alp LAM_{alp j} BF_{alp i bet j}
C                   ---------------------------------------------------
c... unused version from F matrix before merge into CC_XIETA version
c                   IF (MULD2H(ISYAL,ISYMJ).EQ.ISYMPC) THEN
c
c                     KLAMD = IGLMRH(ISYAL,ISYMJ) + NBAS(ISYAL)*(J-1)+1
c                     KGIM  = IT1AO(ISYBE,ISYMI)  + NBAS(ISYBE)*(I-1)+1
c                     NBASA = MAX(NBAS(ISYAL),1)
c
c                     CALL DGEMV('N',NBAS(ISYAL),NBAS(ISYBE),
c    *                           -ONE,WORK(KSCR1),NBASA,XLAMPC(KLAMD),1,
c    *                           ONE,GIMC(KGIM),1)
c
c                   END IF

C                   ---------------------------------------------------
C                   G_{alp j} += - sum_bet LAM_{bet i} BF_{alp i bet j}
C                   ---------------------------------------------------
                    IF (MULD2H(ISYBE,ISYMI).EQ.ISYMPC) THEN

                      KLAMD = IGLMRH(ISYBE,ISYMI) + NBAS(ISYBE)*(I-1)+1
                      KGIM  = IT1AO(ISYAL,ISYMJ)  + NBAS(ISYAL)*(J-1)+1
                      NBASA = MAX(NBAS(ISYAL),1)

                      CALL DGEMV('N',NBAS(ISYAL),NBAS(ISYBE),
     *                          -ONE,WORK(KSCR1),NBASA,XLAMPC(KLAMD),1,
     *                           ONE,GIMC(KGIM),1)

                    END IF 

                  END IF
C
C============================================================
C                 Section for calculating the LT21BF-term.
C============================================================
C
                     IF (ICON.EQ.3 .OR. ICON.EQ.6) THEN
C
                        IF (.NOT. LO3BF) THEN
                          ISYMK = MULD2H(ISYBE,ISYMPB)
                          ISYMD = MULD2H(ISYAL,ISYMPC)
                          ISYMC = MULD2H(ISYMK,ISYMO1)
                          ISYDI = MULD2H(ISYMD,ISYMI)
                          ISYCJ = MULD2H(ISYMC,ISYMJ)
                        ELSE
                          ISYAIJ  = MULD2H(ISYAL,ISYMIJ)
                          ISYMK   = MULD2H(ISYAIJ,ISYMBF)
                          ISYMD   = MULD2H(ISYAL,ISYMPC)
                          ISYMO1  = MULD2H(ISYMC2,MULD2H(ISYMBF,ISYMPC))
                          ISYMC   = MULD2H(ISYMK,ISYMO1)
                          ISYDI   = MULD2H(ISYMD,ISYMI)
                          ISYCJ   = MULD2H(ISYMC,ISYMJ)
                        END IF
C
                        LENGTH = NBAS(ISYAL)*NRHF(ISYMK)
C
                        CALL DZERO(WORK(KSCR2),LENGTH)
C
C----------------------------------------------------------
C                       Transform the AO-block to MO-basis.
C----------------------------------------------------------
C
                        KOFF1  = IGLMRH(ISYBE,ISYMK) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTBE = MAX(NBAS(ISYBE),1)
C
                        IF (.NOT.LO3BF) THEN
                           CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYMK),
     *                              NBAS(ISYBE),ONE,WORK(KSCR1),NTOTAL,
     *                              XLAMPB(KOFF1),NTOTBE,ZERO,
     *                              WORK(KSCR2),NTOTAL)
                           KOFF0 = KSCR2
                        ELSE
                           KOFF0 = KSCR1
                        END IF
C
                        KOFF2  = IGLMVI(ISYAL,ISYMD) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTK  = MAX(NRHF(ISYMK),1)
C
                        CALL DGEMM('T','N',NRHF(ISYMK),NVIR(ISYMD),
     *                             NBAS(ISYAL),ONE,WORK(KOFF0),NTOTAL,
     *                             XLAMPC(KOFF2),NTOTAL,ZERO,
     *                             WORK(KSCR3),NTOTK)
C
C-----------------------------------------------------------------
C                       Contraction with CTR2 & storage in result.
C-----------------------------------------------------------------
C
                        DO C = 1,NVIR(ISYMC)
C
                           NCJ   = IT1AM(ISYMC,ISYMJ)
     *                           + NVIR(ISYMC)*(J - 1) + C
                           NDICJ = IT2SQ(ISYDI,ISYCJ)
     *                           + NT1AM(ISYDI)*(NCJ - 1)
     *                           + IT1AM(ISYMD,ISYMI)
     *                           + NVIR(ISYMD)*(I - 1) + 1
                           NCK   = IT1AM(ISYMC,ISYMK) + C
C
                           CALL DGEMV('N',NRHF(ISYMK),NVIR(ISYMD),
     *                                -ONE,WORK(KSCR3),NTOTK,
     *                                CTR2(NDICJ),1,ONE,RHO1(NCK),
     *                                NVIR(ISYMC))
C
                        END DO
C
C-----------------------------------------------------------------
C                     for ICON=6 add extra relaxation contribution
C                     where XLAMDB and XLAMDC are exchanged
C-----------------------------------------------------------------
C
                      IF (ICON.EQ.6) THEN
C
                        ISYMK = MULD2H(ISYBE,ISYMPC)
                        ISYMD = MULD2H(ISYAL,ISYMPB)
                        ISYMC = MULD2H(ISYMK,ISYMO1)
                        ISYDI = MULD2H(ISYMD,ISYMI)
                        ISYCJ = MULD2H(ISYMC,ISYMJ)
C
                        LENGTH = NBAS(ISYAL)*NRHF(ISYMK)
C
                        CALL DZERO(WORK(KSCR2),LENGTH)
C
C----------------------------------------------------------
C                       Transform the AO-block to MO-basis.
C----------------------------------------------------------
C
                        KOFF1  = IGLMRH(ISYBE,ISYMK) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTBE = MAX(NBAS(ISYBE),1)
C
                        IF (.NOT.LO3BF) THEN
                           CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYMK),
     *                              NBAS(ISYBE),ONE,WORK(KSCR1),NTOTAL,
     *                              XLAMPC(KOFF1),NTOTBE,ZERO,
     *                              WORK(KSCR2),NTOTAL)
                           KOFF0 = KSCR2
                        ELSE
                           KOFF0 = KSCR1
                        END IF
C
                        KOFF2  = IGLMVI(ISYAL,ISYMD) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTK  = MAX(NRHF(ISYMK),1)
C
                        CALL DGEMM('T','N',NRHF(ISYMK),NVIR(ISYMD),
     *                             NBAS(ISYAL),ONE,WORK(KOFF0),NTOTAL,
     *                             XLAMPB(KOFF2),NTOTAL,ZERO,
     *                             WORK(KSCR3),NTOTK)
C
C-----------------------------------------------------------------
C                       Contraction with CTR2 & storage in result.
C-----------------------------------------------------------------
C
                        DO C = 1,NVIR(ISYMC)
C
                           NCJ   = IT1AM(ISYMC,ISYMJ)
     *                           + NVIR(ISYMC)*(J - 1) + C
                           NDICJ = IT2SQ(ISYDI,ISYCJ)
     *                           + NT1AM(ISYDI)*(NCJ - 1)
     *                           + IT1AM(ISYMD,ISYMI)
     *                           + NVIR(ISYMD)*(I - 1) + 1
                           NCK   = IT1AM(ISYMC,ISYMK) + C
C
                           CALL DGEMV('N',NRHF(ISYMK),NVIR(ISYMD),
     *                                -ONE,WORK(KSCR3),NTOTK,
     *                                CTR2(NDICJ),1,ONE,RHO1(NCK),
     *                                NVIR(ISYMC))
C
                        END DO
C
                      END IF
C
                     ENDIF
C
C-------------------------------------------------------------
C                    Transform the AB block to occupied space.
C-------------------------------------------------------------
C
                     IF (.NOT. LGAMMA) GOTO 999
C
                     IF (ICON.NE.4) THEN
C
C                    IF (ISYMPB.NE.1) THEN
C                      WRITE (LUPRI,*) 
C    &                         'ERROR> CCT2MO not implemented for',
C    &                         ' LGAMMA = .TRUE. and ISYMPB.NE.1 !'
C                      CALL QUIT('ERROR> Illegal parameters in CCT2MO.')
C                    END IF
C
                     IF (.NOT.LO3BF) THEN
                       ISYMK = MULD2H(ISYAL,ISYMPB)
                       ISYML = MULD2H(ISYBE,ISYMPB)
                     ELSE
                       ISYAIJ  = MULD2H(ISYAL,ISYMIJ)
                       ISYML   = MULD2H(ISYAIJ,ISYMBF)
                       ISYMK   = MULD2H(ISYAL,ISYMPB)
                     END IF
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NRHFK = MAX(NRHF(ISYMK),1)
C
                     IF (.NOT.LO3BF) THEN
                       KOFF1 = IGLMRH(ISYBE,ISYML) + 1
                       CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYML),
     *                            NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                            XLAMPB(KOFF1),NBASB,ZERO,WORK(KSCR4),
     *                            NBASA)
                       KOFF0 = KSCR4
                     ELSE
                       KOFF0 = KSCR1
                     END IF
C
                     KOFF2 = IGLMRH(ISYAL,ISYMK) + 1
                     CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYML),
     *                          NBAS(ISYAL),ONE,XLAMPB(KOFF2),NBASA,
     *                          WORK(KOFF0),NBASA,ZERO,WORK(KSCR5),
     *                          NRHFK)
C
C-------------------------------------------
C                    Store the gamma matrix.
C-------------------------------------------
C
                     ISYMKI = MULD2H(ISYMK,ISYMI)
                     ISYMLJ = MULD2H(ISYML,ISYMJ)
C
                     DO L = 1,NRHF(ISYML)
C
                        NLJ = IMATIJ(ISYML,ISYMJ)
     *                      + NRHF(ISYML)*(J - 1) + L
C
                        DO K = 1,NRHF(ISYMK)
C
                           NKL = KSCR5 + NRHF(ISYMK)*(L - 1) + K - 1
C
                           NKI = IMATIJ(ISYMK,ISYMI)
     *                         + NRHF(ISYMK)*(I - 1) + K
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              IF (NKI.LE.NLJ) THEN
                                NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                                + INDEX(NKI,NLJ)
                                GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(NKL)
                              END IF
                           ELSE IF (ISYMKI .LT. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(NKL)
                           ENDIF

 
                        END DO
                     END DO

                     END IF
C
C----------------------------------------------------------
C                    for ICON=4 calculate relaxed GAMMA
C                    consisting of two contributions
C----------------------------------------------------------
C
                     IF (ICON.EQ.4) THEN
C
C                    ---------------------------------
C                    calculate the first contribution:
C                    ---------------------------------
C
                     ISYML = MULD2H(ISYBE,ISYMPC)
                     ISYMK = MULD2H(ISYAL,ISYMPB)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NRHFK = MAX(NRHF(ISYMK),1)
C
                     KOFF1 = IGLMRH(ISYBE,ISYML) + 1
                     CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYML),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMPC(KOFF1),NBASB,ZERO,WORK(KSCR4),
     *                          NBASA)
C
                     KOFF2 = IGLMRH(ISYAL,ISYMK) + 1
                     CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYML),
     *                          NBAS(ISYAL),ONE,XLAMPB(KOFF2),NBASA,
     *                          WORK(KSCR4),NBASA,ZERO,WORK(KSCR5),
     *                          NRHFK)
C
C
C                    ------------------------
C                    Store first contribution
C                    ------------------------
C
                     ISYMKI = MULD2H(ISYMK,ISYMI)
                     ISYMLJ = MULD2H(ISYML,ISYMJ)
C
                     DO L = 1,NRHF(ISYML)
C
                        NLJ = IMATIJ(ISYML,ISYMJ)
     *                      + NRHF(ISYML)*(J - 1) + L
C
                        DO K = 1,NRHF(ISYMK)
C
                           NKL = KSCR5 + NRHF(ISYMK)*(L - 1) + K - 1
C
                           NKI = IMATIJ(ISYMK,ISYMI)
     *                         + NRHF(ISYMK)*(I - 1) + K
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              IF (NKI.LE.NLJ) THEN
                                NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                                + INDEX(NKI,NLJ)
                                GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(NKL)
                              END IF
                           ELSE IF (ISYMKI .LT. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(NKL)
                           ENDIF
C
                        END DO
                     END DO
C
C                    -----------------------------
C                    calculate second contribution
C                    -----------------------------
C
                     ISYML = MULD2H(ISYBE,ISYMPB)
                     ISYMK = MULD2H(ISYAL,ISYMPC)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NRHFK = MAX(NRHF(ISYMK),1)
C
                     KOFF1 = IGLMRH(ISYBE,ISYML) + 1
                     CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYML),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMPB(KOFF1),NBASB,ZERO,WORK(KSCR4),
     *                          NBASA)
C
                     KOFF2 = IGLMRH(ISYAL,ISYMK) + 1
                     CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYML),
     *                          NBAS(ISYAL),ONE,XLAMPC(KOFF2),NBASA,
     *                          WORK(KSCR4),NBASA,ZERO,WORK(KSCR5),
     *                          NRHFK)

C
C                    -------------------------------------------
C                    Add second contribution to the gamma matrix
C                    -------------------------------------------
C
                     ISYMKI = MULD2H(ISYMK,ISYMI)
                     ISYMLJ = MULD2H(ISYML,ISYMJ)
C
                     DO L = 1,NRHF(ISYML)
C
                        NLJ = IMATIJ(ISYML,ISYMJ)
     *                      + NRHF(ISYML)*(J - 1) + L
C
                        DO K = 1,NRHF(ISYMK)
C
                           NKL = KSCR5 + NRHF(ISYMK)*(L - 1) + K - 1
C
                           NKI = IMATIJ(ISYMK,ISYMI)
     *                         + NRHF(ISYMK)*(I - 1) + K
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              IF (NKI.LE.NLJ) THEN
                                NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                                + INDEX(NKI,NLJ)
                                GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(NKL)
                              END IF
                           ELSE IF (ISYMKI .LT. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(NKL)
                           ENDIF
C
                        END DO
                     END DO
C
                     END IF
C
  999                CONTINUE
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C
