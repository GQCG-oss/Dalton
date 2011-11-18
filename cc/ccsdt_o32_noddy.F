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
      SUBROUTINE CCSDT_O32_NODDY(LISTX,IDLSTX,LISTY,IDLSTY,
     *                           LISTXY,IDLSTXY,
     *                           XLAMDP0,XLAMDH0,FOCK0,
     *                           WORK,LWORK)
*---------------------------------------------------------------------*
*     Purpose: compute triples contributions to the rhs vectors for   *
*              the second-order amplitude response equations (O2)     *
*                                                                     *
*              B t^x t^y + A{x} t^y + A{y} t^x                        *
*                                                                     *
*     Note: preliminary version, only the contributions of            *
*           A{x} t^y + A{y} t^x to the triples rhs vector             *
*           are implemented, contributions from B t^x t^y terms,      *
*           correction to singles/doubles terms are not implemented   *
*           partitioning step is not done...                          *
*                                                                     *
*     Written by Christof Haettig, Mai 2003, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <ccr1rsp.h>

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL LOCDBG, HAVET31
      PARAMETER (LOCDBG=.FALSE., HAVET31=.FALSE.)


      CHARACTER*3 LISTXY, LISTX, LISTY
      INTEGER LUDELD,LUCKJD,LUDKBC
      INTEGER LWORK, IDLSTX, IDLSTY, IDLSTXY

      CHARACTER FNDELD*6, FNCKJD*6, FNDKBC*4
      PARAMETER (FNDELD  = 'CKDELD'  , FNCKJD  = 'CKJDEL')
      PARAMETER (FNDKBC  = 'DKBC')

#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      REAL WORK(LWORK)
      REAL HALF, ONE, ZERO, TWO, DDOT, XNORM
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION HALF, ONE, TWO, ZERO, DDOT, XNORM
#endif
      PARAMETER( ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TWO = 2.0D0 )

      CHARACTER MODEL*10, LABELX*8, LABELY*8
      LOGICAL HAVEX, HAVEY, HAVEXY
      INTEGER KEND1, LWRK1, ISYMX, ISYMY, KXPROP, KYPROP, KXYMAT,
     &        KTAMPX3, KTAMPY3, KWYINT, KTHETAY, KWXINT, KTHETAX,
     &        KTAMP03, NCK, NBJ, NAI, KAIBJCK, KCKAIBJ, KBJCKAI,
     &        KXPROP_AO, KYPROP_AO, KFCKBUF, KEND2, KLAMPX, KLAMPY,
     &        KLAMHX, KLAMHY, KTAMX1, KTAMY1, IOPT, LWRK2, IRREP, IERR,
     &        KT3AMPXY, KTAMP02, KTAMPX2, KTAMPY2, KTHEOCC, KWXYV,
     &        KWXYO, KTHVIRT

      INTEGER ILSTSYM

*---------------------------------------------------------------------*
*     Begin:
*---------------------------------------------------------------------*
      IF (NSYM.GT.1) CALL QUIT('No symmetry in CCSDT_O32_NODDY!')

      ISYMX = ILSTSYM(LISTX,IDLSTX)
      ISYMY = ILSTSYM(LISTY,IDLSTY)

      KEND1 = 1

      KXPROP = KEND1
      KYPROP = KXPROP + NBAST*NBAST
      KXYMAT = KYPROP + NBAST*NBAST
      KEND1  = KXYMAT + NBAST*NBAST

      IF (HAVET31) THEN
        KTAMPX3 = KEND1
        KTAMPY3 = KTAMPX3 + NT1AMX*NT1AMX*NT1AMX
        KEND1   = KTAMPY3 + NT1AMX*NT1AMX*NT1AMX
      END IF

      KWYINT  = KEND1
      KTHETAY = KWYINT  + NT1AMX*NT1AMX*NT1AMX
      KWXINT  = KTHETAY + NT1AMX*NT1AMX*NT1AMX
      KTHETAX = KWXINT  + NT1AMX*NT1AMX*NT1AMX
      KTAMP03 = KTHETAX + NT1AMX*NT1AMX*NT1AMX
      KEND1   = KTAMP03 + NT1AMX*NT1AMX*NT1AMX
      LWRK1   = LWORK   - KEND1

      IF (LWRK1.LT.0)CALL QUIT('Out of memory in CCSDT_O32_NODDY (2)')

*---------------------------------------------------------------------*
*     Open some files needed:
*---------------------------------------------------------------------*
      LUDELD  = -1
      CALL WOPEN2(LUDELD,FNDELD,64,0)

      LUCKJD  = -1
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)

      LUDKBC  = -1
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)

*---------------------------------------------------------------------*
*     Precalculate a number of intermediates needed for the
*     contributions from A{x} t^y + A{y} t^x, i.e. the w and theta
*     intemediates and some one-electron matrices:
*---------------------------------------------------------------------*
      KLAMPX = KEND1
      KLAMHX = KLAMPX + NLAMDT
      KLAMPY = KLAMHX + NLAMDT
      KLAMHY = KLAMPY + NLAMDT
      KEND2  = KLAMHY + NLAMDT
      LWRK2  = LWORK  - KEND2
      IF (LWRK2.LT.0)CALL QUIT('Out of memory in CCSDT_O32_NODDY (2a)')

      CALL CCSDT_T32_AAPREP_NODDY(
     *           LISTX,IDLSTX,ISYMX,WORK(KWXINT),WORK(KTHETAX),
     *           WORK(KXPROP),WORK(KLAMPX),WORK(KLAMHX),WORK(KTAMPX3),
     *           LISTY,IDLSTY,ISYMY,WORK(KWYINT),WORK(KTHETAY),
     *           WORK(KYPROP),WORK(KLAMPY),WORK(KLAMHY),WORK(KTAMPY3),
     *           WORK(KTAMP03),WORK(KXYMAT),XLAMDP0,XLAMDH0,FOCK0,
     *           LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     *           LOCDBG,HAVET31,WORK(KEND2),LWRK2)

*---------------------------------------------------------------------*
*     Compute contributions:
*---------------------------------------------------------------------*
      HAVEX  = (LISTX(1:3).EQ.'R1 ')
      HAVEY  = (LISTY(1:3).EQ.'R1 ')
      HAVEXY = HAVEX .OR. HAVEY

      KT3AMPXY = KEND1
      KEND1    = KT3AMPXY + NT1AMX*NT1AMX*NT1AMX

      KTAMP02 = KEND1
      KTAMPX2 = KTAMP02 + NT1AMX*NT1AMX
      KTAMPY2 = KTAMPX2 + NT1AMX*NT1AMX
      KTHEOCC = KTAMPY2 + NT1AMX*NT1AMX
      KWXYO   = KTHEOCC + NT1AMX*NRHFT*NRHFT
      KWXYV   = KWXYO   + NT1AMX*NRHFT*NRHFT
      KEND2   = KWXYV   + NT1AMX*NRHFT*NRHFT

      LWRK2    = LWORK  - KEND2
      IF (LWRK2.LT.NT2AMX)
     &    CALL QUIT('Out of memory in CCSDT_O32_NODDY (3)')

      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMOP,IOPT,MODEL,DUMMY,WORK(KEND2))
      CALL CC_T2SQ(WORK(KEND2),WORK(KTAMP02),ISYM0)

      IOPT = 2
      CALL CC_RDRSP(LISTX,IDLSTX,ISYMX,IOPT,MODEL,DUMMY,WORK(KEND2))
      CALL CCLR_DIASCL(WORK(KEND2),TWO,ISYMX)
      CALL CC_T2SQ(WORK(KEND2),WORK(KTAMPX2),ISYMX)

      IOPT = 2
      CALL CC_RDRSP(LISTY,IDLSTY,ISYMY,IOPT,MODEL,DUMMY,WORK(KEND2))
      CALL CCLR_DIASCL(WORK(KEND2),TWO,ISYMY)
      CALL CC_T2SQ(WORK(KEND2),WORK(KTAMPY2),ISYMY)


      CALL DZERO(WORK(KT3AMPXY),NT1AMX*NT1AMX*NT1AMX)

      CALL CCSDT_O32VIR_NODDY(WORK(KWXINT),WORK(KXPROP),
     &                        WORK(KTAMPX2),HAVEX,
     &                        WORK(KWYINT),WORK(KYPROP),
     &                        WORK(KTAMPY2),HAVEY,
     &                        WORK(KTAMP03),WORK(KTAMP02),
     &                        WORK(KXYMAT),HAVEXY,
     &                        WORK(KTHEOCC),WORK(KWXYV),WORK(KWXYO),
     &                        WORK(KT3AMPXY),.TRUE.,
     &                        DUMMY,.FALSE.,DUMMY,.FALSE. )

      KTHVIRT = KTAMPY2 + NT1AMX*NT1AMX
      KWXYV   = KTHVIRT + NT1AMX*NVIRT*NVIRT
      KWXYO   = KWXYV   + NT1AMX*NVIRT*NVIRT
      KEND2   = KWXYO   + NT1AMX*NVIRT*NVIRT

      LWRK2    = LWORK  - KEND2
      IF (LWRK2.LT.0)
     &    CALL QUIT('Out of memory in CCSDT_O32_NODDY (4)')

      CALL CCSDT_O32OCC_NODDY(WORK(KTHETAX),WORK(KXPROP),HAVEX,
     &                        WORK(KTHETAY),WORK(KYPROP),HAVEY,
     &                        WORK(KT3AMPXY),WORK(KTHVIRT),
     &                        WORK(KWXYV),WORK(KWXYO),.TRUE.,.TRUE.,
     &                        DUMMY,.FALSE.,DUMMY,.FALSE.)

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_O32_NODDY> (incomplete) 2.ord. rhs O^XY:'
C       CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,HALF,WORK(KT3AMPXY),1)
        CALL CCSDT_CLEAN_T3(WORK(KT3AMPXY),NT1AMX,NVIRT,NRHFT)
        CALL PRINT_PT3_NODDY(WORK(KT3AMPXY))
      END IF

*---------------------------------------------------------------------*
* print debug output, close files and return:
*---------------------------------------------------------------------*
      CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_O32_NODDY                       *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_O32VIR_NODDY(WXINT,XPROP,TAMPX2,HAVEX,
     &                              WYINT,YPROP,TAMPY2,HAVEY,
     &                              TAMP03,TAMP02,XYMAT,HAVEXY,
     &                              THEOCC,WXYV,WXYO,
     &                              T3AMPXY,LT3AMPXY,
     &                              WXYOINT,LWXYOINT,
     &                              WXYVINT,LWXYVINT)
*---------------------------------------------------------------------*
*     Purpose: set up loop over virtual, sort Wbar, W, theta as they  *
*              would be sorted in the real code and compute the       *
*              contributions to second-order rhs vectors              *
*                                                                     *
*     Written by Christof Haettig, May 2003, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>  

      LOGICAL HAVEX, HAVEY, HAVEXY, LT3AMPXY, LWXYOINT, LWXYVINT

#if defined (SYS_CRAY)
      REAL XPROP(NORBT,NORBT), YPROP(NORBT,NORBT)
      REAL XYMAT(NORBT,NORBT)
      REAL WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WXINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TAMP03(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WXYOINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WXYVINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL T3AMPXY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THEOCC(NVIRT,NRHFT,NRHFT,NRHFT)
      REAL WXYV(NVIRT,NRHFT,NRHFT,NRHFT)
      REAL WXYO(NVIRT,NRHFT,NRHFT,NRHFT)
      REAL TAMP02(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TAMPX2(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TAMPY2(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL CONTRIB
#else
      DOUBLE PRECISION XPROP(NORBT,NORBT), YPROP(NORBT,NORBT)
      DOUBLE PRECISION XYMAT(NORBT,NORBT)
      DOUBLE PRECISION WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WXINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TAMP03(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WXYOINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WXYVINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION T3AMPXY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THEOCC(NVIRT,NRHFT,NRHFT,NRHFT)
      DOUBLE PRECISION WXYV(NVIRT,NRHFT,NRHFT,NRHFT)
      DOUBLE PRECISION WXYO(NVIRT,NRHFT,NRHFT,NRHFT)
      DOUBLE PRECISION TAMP02(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TAMPX2(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TAMPY2(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION CONTRIB
#endif


      DO D = 1, NVIRT
        DO B = 1, NVIRT

          CALL DZERO(WXYV,NVIRT*NRHFT*NRHFT*NRHFT)
          CALL DZERO(WXYO,NVIRT*NRHFT*NRHFT*NRHFT)

          IF (HAVEX) THEN
c          -----------------------------------------------------------
c          build THEOCC^DB(ai,lm) = w(bm,ai,dl)+w(dl,bm,ai)+w(ai,dl,bm)
c          with the w intermediate for Y
c          (should be available from A density routines !)
c          -----------------------------------------------------------
           DO A = 1, NVIRT
           DO I = 1, NRHFT
             DO L = 1, NRHFT
             DO M = 1, NRHFT
               THEOCC(A,I,L,M) = WYINT(B,M,A,I,D,L)+
     &                           WYINT(D,L,B,M,A,I)+WYINT(A,I,D,L,B,M)
             END DO
             END DO
           END DO
           END DO 
          END IF

          IF (HAVEX) THEN
c          -----------------------------------------------------------
c          Do one-index transformation:  
c              WXYV(ai,lm) -= THEOCC^DB(ci,lm) X_ac 
c              WXYO(ai,lm) += THEOCC^DB(ak,lm) X_ki
c          (should be available by modification of routines for W^X !)
c          -----------------------------------------------------------
           DO A = 1, NVIRT
           DO I = 1, NRHFT
             DO L = 1, NRHFT
             DO M = 1, NRHFT
               DO C = 1, NVIRT
                 WXYV(A,I,L,M) = WXYV(A,I,L,M) - 
     &                THEOCC(C,I,L,M) * XPROP(NRHFT+A,NRHFT+C) 
               END DO
               DO K = 1, NRHFT
                 WXYO(A,I,L,M) = WXYO(A,I,L,M) +
     &                THEOCC(A,K,L,M) * XPROP(K,I) 
               END DO
             END DO
             END DO
           END DO
           END DO 
          END IF

          IF (HAVEY) THEN
c          -----------------------------------------------------------
c          build THEOCC^DB(ai,lm) = w(bm,ai,dl)+w(dl,bm,ai)+w(ai,dl,bm)
c          with the w intermediate for X
c          (should be available from A density routines !)
c          -----------------------------------------------------------
           DO A = 1, NVIRT
           DO I = 1, NRHFT
             DO L = 1, NRHFT
             DO M = 1, NRHFT
               THEOCC(A,I,L,M) = WXINT(B,M,A,I,D,L)+
     &                           WXINT(D,L,B,M,A,I)+WXINT(A,I,D,L,B,M)
             END DO
             END DO
           END DO
           END DO 
          END IF

          IF (HAVEY) THEN
c          -----------------------------------------------------------
c          Do one-index transformation:  
c            WXYV(ai,lm) -= THEOCC^DB(ci,lm) Y_ac 
c            WXYO(ai,lm) += THEOCC^DB(ak,lm) Y_ki
c          (should be available by modification of routines for W^X !)
c          -----------------------------------------------------------
           DO A = 1, NVIRT
           DO I = 1, NRHFT
             DO L = 1, NRHFT
             DO M = 1, NRHFT
               DO C = 1, NVIRT
                 WXYV(A,I,L,M) = WXYV(A,I,L,M) - 
     &                THEOCC(C,I,L,M) * YPROP(NRHFT+A,NRHFT+C) 
               END DO
               DO K = 1, NRHFT
                 WXYO(A,I,L,M) = WXYO(A,I,L,M) +
     &                THEOCC(A,K,L,M) * YPROP(K,I) 
               END DO
             END DO
             END DO
           END DO
           END DO 
          END IF

          IF (HAVEX) THEN
c          -----------------------------------------------------------
c          Add the comutator [[X,T^y_2],T^0_2]
c          (should be (almost) be available from routines for W^X !)
c          -----------------------------------------------------------
           DO A = 1, NVIRT
           DO I = 1, NRHFT
             DO L = 1, NRHFT
             DO M = 1, NRHFT
               DO C = 1, NVIRT
               DO K = 1, NRHFT
                WXYO(A,I,L,M) = WXYO(A,I,L,M) - (
     &           TAMP02(D,L,A,K) * TAMPY2(B,M,C,I) * XPROP(K,NRHFT+C)+
     &           TAMPY2(D,L,A,K) * TAMP02(B,M,C,I) * XPROP(K,NRHFT+C)+
     &           TAMP02(B,M,A,K) * TAMPY2(D,L,C,I) * XPROP(K,NRHFT+C)+
     &           TAMPY2(B,M,A,K) * TAMP02(D,L,C,I) * XPROP(K,NRHFT+C))
               END DO
               END DO
             END DO
             END DO
           END DO
           END DO 
          END IF

          IF (HAVEY) THEN
c          -----------------------------------------------------------
c          Add the comutator [[Y,T^x],T^0_2]
c          (should be (almost) available from routines for W^X !)
c          -----------------------------------------------------------
           DO A = 1, NVIRT
           DO I = 1, NRHFT
             DO L = 1, NRHFT
             DO M = 1, NRHFT
               DO C = 1, NVIRT
               DO K = 1, NRHFT
                WXYO(A,I,L,M) = WXYO(A,I,L,M) - (
     &           TAMP02(D,L,A,K) * TAMPX2(B,M,C,I) * YPROP(K,NRHFT+C)+
     &           TAMPX2(D,L,A,K) * TAMP02(B,M,C,I) * YPROP(K,NRHFT+C)+
     &           TAMP02(B,M,A,K) * TAMPX2(D,L,C,I) * YPROP(K,NRHFT+C)+
     &           TAMPX2(B,M,A,K) * TAMP02(D,L,C,I) * YPROP(K,NRHFT+C))
               END DO
               END DO
             END DO
             END DO
           END DO
           END DO 
          END IF

          IF (HAVEXY) THEN
c          -----------------------------------------------------------
c          Add comutator [(X^Y+Y^X),T^0]:
c          (should be available from routines for W^X !)
c          -----------------------------------------------------------
           DO A = 1, NVIRT
           DO I = 1, NRHFT
             DO L = 1, NRHFT
             DO M = 1, NRHFT
               DO C = 1, NVIRT
                 WXYV(A,I,L,M) = WXYV(A,I,L,M) - 
     &                TAMP03(C,I,D,L,B,M) * XYMAT(NRHFT+A,NRHFT+C) 
               END DO
               DO K = 1, NRHFT
                 WXYO(A,I,L,M) = WXYO(A,I,L,M) +
     &                TAMP03(A,K,D,L,B,M) * XYMAT(K,I) 
               END DO
             END DO
             END DO
           END DO
           END DO 
          END IF

c          -----------------------------------------------------------
c          Assuming that the partioning problem with W intermediates
c          has been solved in the real code, the noddy code uses now
c          a short cut and and applies C^abd_iml while storing in an 
c          N^6 array 
c          -----------------------------------------------------------
           IF (LT3AMPXY) THEN
             DO A = 1, NVIRT
             DO I = 1, NRHFT
               DO L = 1, NRHFT
               DO M = 1, NRHFT
                 CONTRIB = WXYO(A,I,L,M) + WXYV(A,I,L,M)
                 T3AMPXY(A,I,B,M,D,L) = T3AMPXY(A,I,B,M,D,L) + CONTRIB
                 T3AMPXY(B,M,D,L,A,I) = T3AMPXY(B,M,D,L,A,I) + CONTRIB
                 T3AMPXY(D,L,A,I,B,M) = T3AMPXY(D,L,A,I,B,M) + CONTRIB
               END DO
               END DO
             END DO
             END DO 
            END IF

           IF (LWXYOINT) THEN
             DO A = 1, NVIRT
             DO I = 1, NRHFT
              DO L = 1, NRHFT
              DO M = 1, NRHFT
               WXYOINT(A,I,B,M,D,L) = WXYO(A,I,L,M) 
              END DO
              END DO
             END DO
             END DO 
           END IF

           IF (LWXYVINT) THEN
             DO A = 1, NVIRT
             DO I = 1, NRHFT
              DO L = 1, NRHFT
              DO M = 1, NRHFT
               WXYVINT(A,I,B,M,D,L) = WXYV(A,I,L,M) 
              END DO
              END DO
             END DO
             END DO 
           END IF

        END DO
      END DO

C     IF (LWXYOINT) THEN
C       WRITE(LUPRI,*) 'WXYOINT in CCSDT_O32VIR_NODDY:'
C       CALL PRINT_PT3_NODDY(WXYOINT)
C     END IF

C     IF (LWXYVINT) THEN
C       WRITE(LUPRI,*) 'WXYVINT in CCSDT_O32VIR_NODDY:'
C       CALL PRINT_PT3_NODDY(WXYVINT)
C     END IF

      RETURN
      END 
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_O32VIR_NODDY                    *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_O32OCC_NODDY(THETAX,XPROP,HAVEX,
     &                              THETAY,YPROP,HAVEY,
     &                              T3AMPXY,THVIRT,WXYV,WXYO,
     &                              LOTRN,LVTRN,
     &                              TXYOINT,LTXYOINT,
     &                              TXYVINT,LTXYVINT)
*---------------------------------------------------------------------*
*     Purpose: set up loop over occupied, sort Wbar, W, theta as they *
*              would be sorted in the real code and compute the       *
*              contributions to second-order rhs vectors              *
*                                                                     *
*     Written by Christof Haettig, May 2003, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>  

      LOGICAL HAVEX, HAVEY, HAVEXY, LOTRN, LVTRN, LTXYOINT, LTXYVINT

#if defined (SYS_CRAY)
      REAL XPROP(NORBT,NORBT), YPROP(NORBT,NORBT)
      REAL THETAY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TXYOINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TXYVINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL T3AMPXY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THVIRT(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL   WXYO(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL   WXYV(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL CONTRIB
#else
      DOUBLE PRECISION XPROP(NORBT,NORBT), YPROP(NORBT,NORBT)
      DOUBLE PRECISION THETAY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TXYOINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TXYVINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION T3AMPXY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THVIRT(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION   WXYO(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION   WXYV(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION CONTRIB
#endif


      DO L = 1, NRHFT
        DO M = 1, NRHFT

          CALL DZERO(WXYO,NT1AMX*NVIRT*NVIRT)
          CALL DZERO(WXYV,NT1AMX*NVIRT*NVIRT)

          IF (HAVEX) THEN
c          -----------------------------------------------------------
c          build THVIRT^DB(ai,lm) = theta(b-bar m,ai,dl)+
c                   theta(d-bar l,bm,ai) + theta(a-bar i,dl,bm)
c          with the theta intermediate for Y
c          (should be available from A density routines !)
c          -----------------------------------------------------------
           DO I = 1, NRHFT
            DO A = 1, NVIRT
             DO B = 1, NVIRT
              DO D = 1, NVIRT
               THVIRT(B,D,A,I) = THETAY(B,M,A,I,D,L)+
     &                           THETAY(D,L,B,M,A,I)+THETAY(A,I,D,L,B,M)
              END DO
             END DO
            END DO
           END DO 
          END IF

          IF (HAVEX) THEN
c          -----------------------------------------------------------
c          Do one-index transformation:  WXY(bdai) = 
c             + THVIRT^LM(bdci) X_ac - THVIRT^LM(bdak) X_ki
c          -----------------------------------------------------------
           DO I = 1, NRHFT
            DO A = 1, NVIRT
             DO B = 1, NVIRT
             DO D = 1, NVIRT
              IF (LVTRN) THEN
               DO C = 1, NVIRT
                 WXYV(B,D,A,I) = WXYV(B,D,A,I) - 
     &                THVIRT(B,D,C,I) * XPROP(NRHFT+A,NRHFT+C) 
               END DO
              END IF
              IF (LOTRN) THEN
               DO K = 1, NRHFT
                 WXYO(B,D,A,I) = WXYO(B,D,A,I) +
     &                THVIRT(B,D,A,K) * XPROP(K,I) 
               END DO
              END IF
             END DO
             END DO
           END DO
           END DO 
          END IF

          IF (HAVEY) THEN
c          -----------------------------------------------------------
c          build THVIRT^DB(ai,lm) = theta(b-bar m,ai,dl)+
c                   theta(d-bar l,bm,ai) + theta(a-bar i,dl,bm)
c          with the theta intermediate for X
c          (should be available from A density routines !)
c          -----------------------------------------------------------
           DO I = 1, NRHFT
            DO A = 1, NVIRT
             DO B = 1, NVIRT
              DO D = 1, NVIRT
               THVIRT(B,D,A,I) = THETAX(B,M,A,I,D,L)+
     &                           THETAX(D,L,B,M,A,I)+THETAX(A,I,D,L,B,M)
              END DO
             END DO
            END DO
           END DO 
          END IF

          IF (HAVEX) THEN
c          -----------------------------------------------------------
c          Do one-index transformation:  WXY(bdai) = 
c             + THVIRT^LM(bdci) X_ac - THVIRT^LM(bdak) X_ki
c          -----------------------------------------------------------
           DO I = 1, NRHFT
            DO A = 1, NVIRT
             DO B = 1, NVIRT
             DO D = 1, NVIRT
              IF (LVTRN) THEN
               DO C = 1, NVIRT
                 WXYV(B,D,A,I) = WXYV(B,D,A,I) - 
     &                THVIRT(B,D,C,I) * YPROP(NRHFT+A,NRHFT+C) 
               END DO
              END IF
              IF (LOTRN) THEN
               DO K = 1, NRHFT
                 WXYO(B,D,A,I) = WXYO(B,D,A,I) +
     &                THVIRT(B,D,A,K) * YPROP(K,I) 
               END DO
              END IF
             END DO
             END DO
           END DO
           END DO 
          END IF

c          -----------------------------------------------------------
c          Assuming that the partioning problem with W intermediates
c          has been solved in the real code, the noddy code uses now
c          a short cut and and applies C^abd_iml while storing in an 
c          N^6 array 
c          -----------------------------------------------------------
           DO I = 1, NRHFT
           DO A = 1, NVIRT
            DO B = 1, NVIRT
            DO D = 1, NVIRT
             CONTRIB = WXYO(B,D,A,I) + WXYV(B,D,A,I) 
             T3AMPXY(A,I,B,M,D,L) = T3AMPXY(A,I,B,M,D,L) + CONTRIB
             T3AMPXY(B,M,D,L,A,I) = T3AMPXY(B,M,D,L,A,I) + CONTRIB
             T3AMPXY(D,L,A,I,B,M) = T3AMPXY(D,L,A,I,B,M) + CONTRIB
            END DO
            END DO
           END DO
           END DO 

           IF (LTXYVINT) THEN
            DO I = 1, NRHFT
            DO A = 1, NVIRT
             DO B = 1, NVIRT
             DO D = 1, NVIRT
             TXYVINT(A,I,B,M,D,L)=TXYVINT(A,I,B,M,D,L) + WXYV(B,D,A,I)
             TXYVINT(B,M,D,L,A,I)=TXYVINT(B,M,D,L,A,I) + WXYV(B,D,A,I)
             TXYVINT(D,L,A,I,B,M)=TXYVINT(D,L,A,I,B,M) + WXYV(B,D,A,I)
             END DO
             END DO
            END DO
            END DO 
           ENDIF

           IF (LTXYOINT) THEN
            DO I = 1, NRHFT
            DO A = 1, NVIRT
             DO B = 1, NVIRT
             DO D = 1, NVIRT
             TXYOINT(A,I,B,M,D,L)=TXYOINT(A,I,B,M,D,L) + WXYO(B,D,A,I)
             TXYOINT(B,M,D,L,A,I)=TXYOINT(B,M,D,L,A,I) + WXYO(B,D,A,I)
             TXYOINT(D,L,A,I,B,M)=TXYOINT(D,L,A,I,B,M) + WXYO(B,D,A,I)
             END DO
             END DO
            END DO
            END DO 
           ENDIF

        END DO
      END DO

C     IF (LTXYVINT) THEN
C       WRITE(LUPRI,*) 'TXYINT in CCSDT_O32OCC_NODDY:'
C       CALL PRINT_PT3_NODDY(TXYVINT)
C     END IF

      RETURN
      END 
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_O32OCC_NODDY                    *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_T32_AAPREP_NODDY(
     *             LISTX,IDLSTX,ISYMX,WXINT,THETAX,XPROP,
     *             XLAMDPX,XLAMDHX,TAMPX3,
     *             LISTY,IDLSTY,ISYMY,WYINT,THETAY,YPROP,
     *             XLAMDPY,XLAMDHY,TAMPY3,
     *             T3AMP0,XYMAT,XLAMDP0,XLAMDH0,FOCK0,
     *             LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     *             LOCDBG,HAVET31,WORK,LWORK)
*---------------------------------------------------------------------*
* Purpose: prepare intermediates needed for the A{x} t^y + A{y} t^x   *
*          contribution to the triples part of the second-order       *
*          amplitude response                                         *
*                                                                     *
* Written by Christof Haettig, Jun 2003, Friedrichstal                *
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <ccr1rsp.h>
#include <ccnoddy.h>

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER*(3) LISTX, LISTY
      CHARACTER*(*) FNDELD, FNDKBC, FNCKJD
      INTEGER IDLSTX, IDLSTY, LWORK, LUDELD, LUDKBC, LUCKJD,
     &        ISYMX, ISYMY
      LOGICAL LOCDBG, HAVET31

#if defined (SYS_CRAY)
      REAL WXINT(*), THETAX(*), XPROP(*), TAMPX3(*)
      REAL WYINT(*), THETAY(*), YPROP(*), TAMPY3(*)
      REAL T3AMP0(*), XYMAT(*)
      REAL XLAMDP0(*), XLAMDH0(*), FOCK0(*), WORK(*)
      REAL XLAMDPX(*), XLAMDHX(*), XLAMDPY(*), XLAMDHY(*)
      REAL ONE, THREE, FREQX, FREQY
#else
      DOUBLE PRECISION WXINT(*), THETAX(*), XPROP(*), TAMPX3(*)
      DOUBLE PRECISION WYINT(*), THETAY(*), YPROP(*), TAMPY3(*)
      DOUBLE PRECISION T3AMP0(*), XYMAT(*)
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*), WORK(*)
      DOUBLE PRECISION XLAMDPX(*), XLAMDHX(*), XLAMDPY(*), XLAMDHY(*)
      DOUBLE PRECISION ONE, THREE, FREQX, FREQY
#endif
      PARAMETER(ONE = 1.0D0, THREE = 3.0D0)

      CHARACTER*10 MODEL
      CHARACTER*8 LABELX, LABELY
      INTEGER KXPROP_AO, KTAMX1, KYPROP_AO, KTAMY1,
     &        KFCKBUF, KEND2, LWRK2, IRREP, IERR, IOPT,
     &        NCK, NBJ, NAI, KAIBJCK, KCKAIBJ, KBJCKAI,
     &        KFOCKD, KINT1S0, KINT2S0, KLAMPB, KLAMHB, KFOCKB,
     &        KINT1SB, KINT2SB, LUTEMP, LUSIFC

*---------------------------------------------------------------------*
*     Get w^y and theta^y intermediate into memory:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        write(lupri,*) 'entered CCSDT_T32_AAPREP_NODDY...'
        write(lupri,*) 'listx,idlstx:',listx,idlstx
        write(lupri,*) 'listy,idlsty:',listy,idlsty
      END IF

*---------------------------------------------------------------------*
*     Get w^y and theta^y intermediate into memory:
*---------------------------------------------------------------------*
      IF (LISTY(1:3).EQ.'R1 ') THEN

        CALL CCSDT_GWTIC(LISTY,IDLSTY,WYINT,THETAY,
     &                   T3AMP0,HAVET31,TAMPY3,LOCDBG,
     &                   XLAMDP0,XLAMDH0,FOCK0,
     &                   LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     &                   WORK,LWORK)

      ELSE IF (LISTY(1:3).EQ.'RE ' .OR. LISTY(1:3).EQ.'RC ') THEN

        KEND2   = 1

        KFOCKD = KEND2
        KEND2  = KFOCKD + NORBT

        KINT1S0 = KEND2
        KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
        KEND2   = KINT2S0 + NRHFT*NRHFT*NT1AMX
   
        KLAMPB  = KEND2
        KLAMHB  = KLAMPB + NLAMDT
        KFOCKB  = KLAMHB + NLAMDT
        KEND2   = KFOCKB + N2BST(1)

        KINT1SB = KEND2
        KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
        KEND2   = KINT2SB + NRHFT*NRHFT*NT1AMX

        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LT. 0)
     &      CALL QUIT('Out of memory in CCSDT_T32_AAPREP_NODDY...')

C       ---------------------------------------
C       Read precalculated integrals from file:
C           XINT1S0 =  (CK|BD)
C           XINT2S0 =  (CK|LJ)
C       ---------------------------------------
        CALL CCSDT_READ_NODDY(.FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        .FALSE.,DUMMY,DUMMY,
     &                        .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                        .FALSE.,DUMMY,DUMMY,
     &                        .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        REWIND LUSIFC
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ (LUSIFC)
        READ (LUSIFC) (WORK(KFOCKD-1+I), I=1,NORBT)
        CALL GPCLOSE(LUSIFC,'KEEP')


        CALL DZERO(WYINT,NT1AMX*NT1AMX*NT1AMX)

        CALL CCSDT_T31_NODDY(WYINT,LISTY,IDLSTY,FREQY,.FALSE.,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,DUMMY,DUMMY,
     &                       .FALSE.,DUMMY,DUMMY,
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       XLAMDP0,XLAMDH0,FOCK0,DUMMY,WORK(KFOCKD),
     &                       WORK(KEND2),LWRK2) 

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WYINT,1)

        IF (HAVET31) THEN
          CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,WYINT,1,TAMPY3,1)
        END IF

        ! Theta-Y is set to zero for this case
        CALL DZERO(THETAY,NT1AMX*NT1AMX*NT1AMX)

        ! account for the fact that in the real code a combination
        ! of three W intermediates give the complete Tbar3
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,ONE/THREE,WYINT,1)

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (T3AMP0(I),I=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP') 

      ELSE
        CALL QUIT('Unknown or illegal LISTY="'//LISTY(1:3)//
     &            '" in CCSDT_T32_AAPREP_NODDY')
      END IF

C     IF (LOCDBG) THEN
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> W^Y:'
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> LISTY:',LISTY,IDLSTY
C      CALL PRINT_PT3_NODDY(WYINT)
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> Theta^Y:'
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> LISTY:',LISTY,IDLSTY
C      CALL PRINT_PT3_NODDY(THETAY)
C     END IF

*---------------------------------------------------------------------*
*     Get w^x and theta^x intermediate into memory:
*---------------------------------------------------------------------*
      IF (LISTX(1:3).EQ.'R1 ') THEN

        CALL CCSDT_GWTIC(LISTX,IDLSTX,WXINT,THETAX,
     &                   T3AMP0,HAVET31,TAMPX3,LOCDBG,
     &                   XLAMDP0,XLAMDH0,FOCK0,
     &                   LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     &                   WORK,LWORK)

      ELSE IF (LISTX(1:3).EQ.'RE ' .OR. LISTX(1:3).EQ.'RC ') THEN

        KEND2   = 1

        KFOCKD = KEND2
        KEND2  = KFOCKD + NORBT

        KINT1S0 = KEND2
        KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
        KEND2   = KINT2S0 + NRHFT*NRHFT*NT1AMX
   
        KLAMPB  = KEND2
        KLAMHB  = KLAMPB + NLAMDT
        KFOCKB  = KLAMHB + NLAMDT
        KEND2   = KFOCKB + N2BST(1)

        KINT1SB = KEND2
        KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
        KEND2   = KINT2SB + NRHFT*NRHFT*NT1AMX

        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LT. 0)
     &      CALL QUIT('Out of memory in CCSDT_T32_AAPREP_NODDY...')

C       ---------------------------------------
C       Read precalculated integrals from file:
C           XINT1S0 =  (CK|BD)
C           XINT2S0 =  (CK|LJ)
C       ---------------------------------------
        CALL CCSDT_READ_NODDY(.FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        .FALSE.,DUMMY,DUMMY,
     &                        .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                        .FALSE.,DUMMY,DUMMY,
     &                        .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        REWIND LUSIFC
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ (LUSIFC)
        READ (LUSIFC) (WORK(KFOCKD-1+I), I=1,NORBT)
        CALL GPCLOSE(LUSIFC,'KEEP')


        CALL DZERO(WXINT,NT1AMX*NT1AMX*NT1AMX)

        CALL CCSDT_T31_NODDY(WXINT,LISTX,IDLSTX,FREQX,.FALSE.,
     &                       .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                       .FALSE.,DUMMY,DUMMY,
     &                       .FALSE.,DUMMY,DUMMY,
     &                               WORK(KINT1SB),WORK(KINT2SB),
     &                       WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                       XLAMDP0,XLAMDH0,FOCK0,DUMMY,WORK(KFOCKD),
     &                       WORK(KEND2),LWRK2) 

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WXINT,1)

        IF (HAVET31) THEN
          CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,WXINT,1,TAMPX3,1)
        END IF

        ! Theta-Y is set to zero for this case
        CALL DZERO(THETAX,NT1AMX*NT1AMX*NT1AMX)

        ! account for the fact that in the real code a combination
        ! of three W intermediates give the complete Tbar3
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,ONE/THREE,WXINT,1)

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (T3AMP0(I),I=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP') 

      ELSE
        CALL QUIT('Unknown or illegal LISTX="'//LISTY(1:3)//
     &            '" in CCSDT_T32_AAPREP_NODDY')
      END IF

C     IF (LOCDBG) THEN
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> W^X:'
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> LISTX:',LISTX,IDLSTX
C      CALL PRINT_PT3_NODDY(WXINT)
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> Theta^X:'
C      WRITE(LUPRI,*)'CCSDT_T32_AAPREP_NODDY> LISTX:',LISTX,IDLSTX
C      CALL PRINT_PT3_NODDY(THETAX)
C     END IF

*---------------------------------------------------------------------*
*     Allocate some intermediates needed for the next sections:
*---------------------------------------------------------------------*
      KXPROP_AO = 1
      KYPROP_AO = KXPROP_AO + NBAST*NBAST
      KFCKBUF   = KYPROP_AO + NBAST*NBAST
      KEND2     = KFCKBUF   + NBAST*NBAST

      KTAMX1   = KEND2
      KTAMY1   = KTAMX1 + NT1AMX
      KEND2    = KTAMY1 + NT1AMX

      LWRK2    = LWORK  - KEND2
      IF (LWRK2.LT.0) CALL 
     &   QUIT('Out of memory in CCSDT_T32_AAPREP_NODDY')

*---------------------------------------------------------------------*
*     Get property matrices in Lambda-MO basis:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KXPROP_AO),NORBT*NORBT)
      CALL DZERO(WORK(KYPROP_AO),NORBT*NORBT)

      CALL DZERO(XPROP,NORBT*NORBT)
      CALL DZERO(YPROP,NORBT*NORBT)
 
      IF (LISTX(1:3).EQ.'R1 ') THEN
        LABELX = LRTLBL(IDLSTX)

        ! read property integrals from file:
        CALL CCPRPAO(LABELX,.TRUE.,WORK(KXPROP_AO),IRREP,ISYMX,IERR,
     &               WORK(KEND2),LWRK2)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMX)) THEN
          CALL QUIT('CCSDT_T32_NODDY: error reading operator '//LABELX)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(WORK(KXPROP_AO),N2BST(ISYMX))
        END IF
        CALL DCOPY(NORBT*NORBT,WORK(KXPROP_AO),1,XPROP,1)
 
        ! transform property integrals to Lambda-MO basis
        CALL CC_FCKMO(XPROP,XLAMDP0,XLAMDH0,
     &                WORK(KEND2),LWRK2,ISYMX,1,1)
      END IF
 
      IF (LISTY(1:3).EQ.'R1 ') THEN
        LABELY = LRTLBL(IDLSTY)
        ! read property integrals from file:
        CALL CCPRPAO(LABELY,.TRUE.,WORK(KYPROP_AO),IRREP,ISYMY,IERR,
     &               WORK(KEND2),LWRK2)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMY)) THEN
          CALL QUIT('CCSDT_T32_NODDY: error reading operator '//LABELY)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(WORK(KYPROP_AO),N2BST(ISYMY))
        END IF
        CALL DCOPY(NORBT*NORBT,WORK(KYPROP_AO),1,YPROP,1)
 
        ! transform property integrals to Lambda-MO basis
        CALL CC_FCKMO(YPROP,XLAMDP0,XLAMDH0,
     &                WORK(KEND2),LWRK2,ISYMY,1,1)
      END IF

*---------------------------------------------------------------------*
*     Construct sum of one-index transformed property matrices:
*     XYMAT = [X,T^Y_1] + [Y,T^X_1]
*---------------------------------------------------------------------*
      IOPT = 1
      CALL CC_RDRSP(LISTX,IDLSTX,ISYMX,IOPT,MODEL,WORK(KTAMX1),DUMMY)

      CALL CCLR_LAMTRA(XLAMDP0,XLAMDPX,XLAMDH0,XLAMDHX,
     &                 WORK(KTAMX1),ISYMX)

      IOPT = 1
      CALL CC_RDRSP(LISTY,IDLSTY,ISYMY,IOPT,MODEL,WORK(KTAMY1),DUMMY)

      CALL CCLR_LAMTRA(XLAMDP0,XLAMDPY,XLAMDH0,XLAMDHY,
     &                 WORK(KTAMY1),ISYMY)


      CALL DZERO(XYMAT,NORBT*NORBT)
 
      IF (LISTX(1:3).EQ.'R1 ') THEN
        ! add [X,T^Y_1]
        CALL DCOPY(NORBT*NORBT,WORK(KXPROP_AO),1,WORK(KFCKBUF),1)
        CALL CC_FCKMO(WORK(KFCKBUF),XLAMDPY,XLAMDH0,
     &                WORK(KEND2),LWRK2,ISYMX,ISYMY,ISYM0)
        CALL DAXPY(NORBT*NORBT,ONE,WORK(KFCKBUF),1,XYMAT,1)
 
        CALL DCOPY(NORBT*NORBT,WORK(KXPROP_AO),1,WORK(KFCKBUF),1)
        CALL CC_FCKMO(WORK(KFCKBUF),XLAMDP0,XLAMDHY,
     &                WORK(KEND2),LWRK2,ISYMX,ISYM0,ISYMY)
        CALL DAXPY(NORBT*NORBT,ONE,WORK(KFCKBUF),1,XYMAT,1)
      END IF
 
      IF (LISTY(1:3).EQ.'R1 ') THEN
        ! add [Y,T^X_1]
        CALL DCOPY(NORBT*NORBT,WORK(KYPROP_AO),1,WORK(KFCKBUF),1)
        CALL CC_FCKMO(WORK(KFCKBUF),XLAMDPX,XLAMDH0,
     &                WORK(KEND2),LWRK2,ISYMY,ISYMX,ISYM0)
        CALL DAXPY(NORBT*NORBT,ONE,WORK(KFCKBUF),1,XYMAT,1)
 
        CALL DCOPY(NORBT*NORBT,WORK(KYPROP_AO),1,WORK(KFCKBUF),1)
        CALL CC_FCKMO(WORK(KFCKBUF),XLAMDP0,XLAMDHX,
     &                WORK(KEND2),LWRK2,ISYMY,ISYM0,ISYMX)
        CALL DAXPY(NORBT*NORBT,ONE,WORK(KFCKBUF),1,XYMAT,1)
      END IF 

      RETURN
      END
*=====================================================================*
