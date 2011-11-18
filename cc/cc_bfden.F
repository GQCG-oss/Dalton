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
C  /* Deck cc_bfden */
      SUBROUTINE CC_BFDEN( T2AM,   ISYTAM, MINT,   ISYMIN, 
     *                     XLAMD1, ISYML1, XLAMD2, ISYML2,
     *                     XLAMD3, ISYML3, XLAMD4, ISYML4,
     *                     FILNAM, LUBFDEN, IADRBF, IADR, 
     *                     IVEC, IOPT, LO3BF, WORK, LWORK )
*---------------------------------------------------------------------*
*
*     Purpose: Calculate the effective density used for the BF terms
*              write it to the direct access file FILNAM/LUBFDEN
*              start addresses are stored in IADRBF
*
*
*     IOPT .EQ. 0:  two-index back transformation of T --> evaluate:
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2)
*                      = sum_ab T^{a,b}_{i,j} XLAMD1(d,b) XLAMD2(g,a)
*
*        XLAMD1,XLAMD2 : (ordinary) lambda hole matrices
*        XLAMD3,XLAMD4,MINT unused
*
*
*     IOPT .EQ. 1:  vector function --> evaluate:
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2)
*                             + XLAMD1(g,i)*XLAMD3(d,j)
*
*        XLAMD1,XLAMD2,XLAMD3 : (ordinary) lambda hole matrices
*        XLAMD4,MINT unused
*
*
*     IOPT .EQ. 2:   Jacobian left transf. --> evaluate:
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2) + Mint^{g,d}_{i,j}(1,2)
*              + 2 XLAMD1(g,i)*XLAMD3(d,j) - XLAMD3(d,i)*XLAMD1(g,j)
*              + 2 XLAMD3(g,i)*XLAMD1(d,j) - XLAMD1(d,i)*XLAMD3(g,j)
*
*        XLAMD1: zeroth-order lambda particle matrix
*        XLAMD2: zeroth-order lambda particle matrix
*        XLAMD3: zeroth-order chi intermediate
*        XLAMD4: unused
*        MINT  : M intermediate in left transf.
*
*
*
*     IOPT .EQ. 3:   Jacobian right transformation and
*                    B matrix transformation  --> evaluate:
*     !!! not yet tested for this version !!!
*
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2)
*                             + XLAMD1(gam,i)*XLAMD3(del,j)
*                             + XLAMD3(gam,i)*XLAMD1(del,j)
*
*        XLAMD1,XLAMD2,XLAMD3 : lambda hole matrices
*        XLAMD4,MINT unused
*
*
*
*
*     IOPT .EQ. 4: F-matrix transformation --> evaluate
*     !!! not yet tested for this version !!!
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2) + T^{g,d}_{i,j}(2,1)
*              + 2 XLAMD1(g,i)*XLAMD3(d,j) - XLAMD3(d,i)*XLAMD1(g,j)
*              + 2 XLAMD3(g,i)*XLAMD1(d,j) - XLAMD1(d,i)*XLAMD3(g,j)
*
*           XLAMD1: zeroth-order lambda particle matrix
*           XLAMD2: first-order  lambda particle matrix
*           XLAMD3: first-order  one-index backtransf. Zeta1 amplitudes
*           XLAMD4: unused
*
*
*
*     IOPT .EQ. 5: relaxation contrib. to Jacobian right trans. -->
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2) + T^{g,d}_{i,j}(2,1)
*                        + XLAMD1(g,i)*XLAMD3(d,j) 
*                        + XLAMD2(g,i)*XLAMD4(d,j)
*
*           XLAMD1: first-order  lambda hole matrix
*           XLAMD2: zeroth-order lambda hole matrix
*           XLAMD3: zeroth-order lambda hole matrix
*           XLAMD4: first-order  lambda hole matrix
*
*
*     IOPT .EQ. 6: relaxation contrib. to Jacobian left trans. -->
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2) + T^{g,d}_{i,j}(2,1)
*              + 2 XLAMD1(g,i)*XLAMD3(d,j) - XLAMD3(d,i)*XLAMD1(g,j)
*              + 2 XLAMD3(g,i)*XLAMD1(d,j) - XLAMD1(d,i)*XLAMD3(g,j)
*              + 2 XLAMD2(g,i)*XLAMD4(d,j) - XLAMD4(d,i)*XLAMD2(g,j)
*              + 2 XLAMD4(g,i)*XLAMD2(d,j) - XLAMD2(d,i)*XLAMD4(g,j)
*
*           XLAMD1: first-order  lambda hole matrix
*           XLAMD2: zeroth-order lambda hole matrix
*           XLAMD3: zeroth-order chi intermediate
*           XLAMD4: first-order  chi intermediate
*
*     IOPT .EQ. 7: relax. + resp. contribs to eff. dens. for BF(A,QA)
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2) + T^{g,d}_{i,j}(2,1)
*              + XLAMD1(g,i)*XLAMD3(d,j) + XLAMD3(g,i)*XLAMD1(d,j)
*              + XLAMD2(g,i)*XLAMD4(d,j) + XLAMD4(g,i)*XLAMD2(d,j)
*
*           XLAMD1: first-order  lambda hole matrix
*           XLAMD2: zeroth-order lambda hole matrix
*           XLAMD3: T^A tranformed zero-order lambda hole matrix
*           XLAMD4: T^A tranformed first-order lambda hole matrix
*
*
*     IOPT .EQ. 8 used for F matrix transformation
*
*     IOPT .EQ. 9:  symmetrized two-index back transformation of T 
*
*        M^{g,d}_{i,j} = T^{g,d}_{i,j}(1,2) + T^{g,d}_{i,j}(2,1)
*
*        XLAMD1,XLAMD2 : lambda hole matrices
*        XLAMD3,XLAMD4,MINT unused
*
*
*
*     LO3BF - flag, for special treatment in F matrix
*             (unused at present?)
*
*     SYMMETRIES:
* 
*     ISYML1,ISYML2,ISYML3,ISYML4  --  XLAMD1, XLAMD2, XLAMD3, XLAMD4
*     ISYTAM                       --  T2AM
*     ISYMIN                       --  MINT
*
*     Christof Haettig, September 1998
*     IOPT = 7 added, Sonia Coriani, June 1999
*     IOPT = 0 added, Christof Haettig, Februar 2000
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <second.h>
C
      LOGICAL LO3BF
      CHARACTER*(*) FILNAM
      INTEGER LUBFDEN, IVEC, IADR, IADRBF(MXCORB,*)
      INTEGER ISYML1,ISYML2,ISYML3,ISYML4,ISYTAM,ISYMIN,IOPT,LWORK
      INTEGER KEND1,ISYMT1,ISYMT2,ISYDEN,ISYMD,ISYMM1,ISYMM2,NSCRM,NMGD
      INTEGER KSCRM,KMGD,LWRK1,ISYMGD,IOPT2AO,IDEL,IOPBFAO 
C
#if defined (SYS_CRAY)
      REAL XLAMD1(*),XLAMD2(*),XLAMD3(*),XLAMD4(*)
      REAL T2AM(*), MINT(*), WORK(*)
      REAL TIMIO, DTIME, ZERO, ONE, DDOT, XNORM
#else
      DOUBLE PRECISION XLAMD1(*),XLAMD2(*),XLAMD3(*),XLAMD4(*)  
      DOUBLE PRECISION T2AM(*), MINT(*), WORK(*)
      DOUBLE PRECISION TIMIO, DTIME, ZERO, ONE, DDOT, XNORM
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
C-----------------------
C     symmetry analysis:
C-----------------------
C
      ISYMT1 = MULD2H(ISYTAM,ISYML1)
      ISYMT2 = MULD2H(ISYTAM,ISYML2)
      ISYDEN = MULD2H(ISYMT1,ISYML2)
C
      TIMIO = ZERO
C
      IF ( IOPT.NE.0 .AND. IOPT.NE.9) THEN
        IF ( ISYML3.NE.ISYMT2 ) THEN
          WRITE (LUPRI,*) 'ISYML3:',ISYML3
          WRITE (LUPRI,*) 'ISYMT2:',ISYMT2
          CALL QUIT('CC_BFDEN: Symmetry mismatch')
        EN DIF
      END IF
C      
      IF ( (IOPT.EQ.5) .OR. (IOPT.EQ.6) .OR. (IOPT.EQ.7)) THEN
         IF ( MULD2H(ISYML4,ISYML2) .NE. ISYDEN) THEN
            WRITE (LUPRI,*) 'IOPT:',IOPT
            WRITE (LUPRI,*) 'ISYML4:',ISYML4
            WRITE (LUPRI,*) 'ISYML2:',ISYML2
            WRITE (LUPRI,*) 'ISYDEN:',ISYDEN
            CALL QUIT('CC_BFDEN: Symmetry mismatch')
         ENDIF
      END IF
C
C---------------------------------
C     Start loop over Delta index:
C---------------------------------
C
      DO ISYMD = 1, NSYM

         ISYMM1 = MULD2H(ISYMT1,ISYMD)
         ISYMM2 = MULD2H(ISYMT2,ISYMD)
         ISYMGD = MULD2H(ISYDEN,ISYMD)
     
         NSCRM = MAX(NT2BCD(ISYMM1),NT2BCD(ISYMM2))
         NMGD  = NT2BGD(ISYMGD)
         IF (LO3BF .AND. IOPT.EQ.4) NMGD = NMAIJK(ISYMGD)

         KSCRM = 1
         KMGD  = KSCRM + NSCRM
         KEND1 = KMGD  + NMGD
         LWRK1 = LWORK - KEND1

         IF (LWRK1 .LT. 0) THEN
            WRITE (LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
            WRITE (LUPRI,*) 'Insufficient space in CC_BFDEN.'
            CALL QUIT('Insufficient space in CC_BFDEN.')
         ENDIF
C
         DO D = 1, NBAS(ISYMD)

            IDEL = D + IBAS(ISYMD)
C
C           -----------------------------------------------------------
C           calculate one-index AO back-transformed T2 amplitudes
C           back-transforme Delta index with XLAMD1,
C           -----------------------------------------------------------
C
            IOPT2AO = 1
            CALL CC_T2AO(T2AM,XLAMD1,ISYML1,WORK(KSCRM),WORK(KEND1),
     *                   LWRK1,IDEL,ISYMD,ISYTAM,IOPT2AO)
C
C           -----------------------------------------------------------
C           back-transforme Gamma index with XLAMD2,
C           add "F-term" from XLAMD1 x XLAMD3 (MGD is initialized here)
C           -----------------------------------------------------------
C
            IF (IOPT .EQ. 7) THEN
              IOPBFAO = 3
            ELSE IF (IOPT.EQ. 9) THEN
              IOPBFAO = 0
            ELSE
              IOPBFAO = IOPT
            ENDIF

            CALL CC_BFAO(WORK(KMGD),WORK(KSCRM),ISYMM1,XLAMD2,ISYML2,
     *                D,ISYMD,XLAMD3,ISYML3,XLAMD1,ISYML1,ZERO,IOPBFAO)
C
C           -----------------------------------------------------------
C           for IOPT=4,5,6,7,9 add extra response (or relax.) contribs.
C              back-transforme Delta index with XLAMD2,
C              back-transforme Gamma index with XLAMD1,
C              add "F-term" from XLAMD2 x XLAMD4 
C           -----------------------------------------------------------
C
            IF (IOPT.EQ.4 .OR. IOPT.EQ.5 .OR. IOPT.EQ.6 .OR. 
     *          IOPT.EQ.7 .OR. IOPT.EQ.9                     ) THEN

              IOPT2AO = 1
              CALL CC_T2AO(T2AM,XLAMD2,ISYML2,WORK(KSCRM),WORK(KEND1),
     *                     LWRK1,IDEL,ISYMD,ISYTAM,IOPT2AO)
 
              IF (IOPT .EQ. 7) THEN
                 IOPBFAO = 3
              ELSE IF (IOPT.EQ. 9) THEN
                 IOPBFAO = 0
              ELSE
                 IOPBFAO = IOPT
              ENDIF

              CALL CC_BFAO(WORK(KMGD),WORK(KSCRM),ISYMM2,XLAMD1,ISYML1,
     *                  D,ISYMD,XLAMD4,ISYML4,XLAMD2,ISYML2,ONE,IOPBFAO)

            ENDIF
C
C           -----------------------------------------------------------
C           for IOPT=2,6 include also two-index backtransformed
C           M intermediate (transformed with XLAMD1, XLAMD2)
C           -----------------------------------------------------------
C
            IF ( (IOPT.EQ.2) .OR. (IOPT.EQ.6) .OR. (IOPT.EQ.8) ) THEN

              CALL CC_BFMIAO(WORK(KMGD),MINT,ISYMIN,XLAMD1,ISYML1,
     &                       XLAMD2,ISYML2,D,ISYMD,WORK(KEND1),LWRK1)

              IF (IOPT.EQ.6) THEN
                CALL CC_BFMIAO(WORK(KMGD),MINT,ISYMIN,XLAMD2,ISYML2,
     &                         XLAMD1,ISYML1,D,ISYMD,WORK(KEND1),LWRK1)
              END IF

            END IF
C
C           -----------------------------------------------------------
C           save effective density on file:
C           -----------------------------------------------------------
C
            IADRBF(IDEL,IVEC) = IADR

            DTIME = SECOND()
            CALL PUTWA2(LUBFDEN,FILNAM,WORK(KMGD),IADR,NMGD)
            TIMIO = TIMIO + SECOND() - DTIME

            IADR = IADR + NMGD

         END DO
      END DO
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck cc_bfdenf */
      SUBROUTINE CC_BFDENF( ZETA2,  ISYCTR, MINT, ISYMIN,
     *                      XLAMDP, ISYXLP, CHI,  ISYCHI,
     *                      T1AM,   ISYTAM, MGD,  WORK, LWORK  )
*---------------------------------------------------------------------*
*
*     Purpose: Calculate the effective density used for the BZeta term
*              in the F matrix transformation
*
*     MGD^{ij}_{alp,k} = Zeta^{k,alp}_{i,j} + M^{k,alp}_{ij}
*               + 2 Chi(k,i) Lambda(alp,j) - Chi(k,j) Lambda(alp,i)
* 
*     Christof Haettig, November 1998
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      INTEGER LWORK, ISYCTR, ISYMIN, ISYXLP, ISYCHI, ISYTAM

#if defined (SYS_CRAY)
      REAL ZETA2(*), MINT(*), CHI(*), T1AM(*), XLAMDP(*)
      REAL MGD(*), WORK(LWORK)
      REAL ZERO, ONE, HALF, TWO, DDOT
#else
      DOUBLE PRECISION ZETA2(*), MINT(*), CHI(*), T1AM(*), XLAMDP(*)
      DOUBLE PRECISION MGD(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, HALF, TWO, DDOT
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)

      INTEGER ISYALP,ISYMM1,ISYMGD,ISYMAK,ISYMI,ISYMK,ISYMJ,ISYMIJ
      INTEGER KSCRM,KMGD,KOFF1,KOFF2,NIJ,KEND1,LWRK1,IOPT2AO,IOPT
      INTEGER IALP, ISYMKI
C
C-----------------------
C     Symmetry analysis:
C-----------------------
C
      IF (      (ISYMIN.NE.MULD2H(ISYCTR,ISYTAM))
     &     .OR. (ISYMIN.NE.ISYCHI)                ) THEN
         WRITE (LUPRI,*) 'ISYMIN,ISYCTR,ISYTAM,ISYCHI:',
     &             ISYMIN,ISYCTR,ISYTAM,ISYCHI
         WRITE (LUPRI,*) 'CC_BFDENF: Symmetry mismatch.'
         CALL QUIT('CC_BFDENF: Symmetry mismatch.')
      END IF
C
C---------------------------------
C     Start loop over Alpha index:
C---------------------------------
C
      DO ISYALP = 1, NSYM

        ISYMM1 = MULD2H(MULD2H(ISYCTR,ISYXLP),ISYALP)
        ISYMGD = MULD2H(ISYMM1,ISYTAM)

        KSCRM = 1
        KMGD  = KSCRM + NT2BCD(ISYMM1)
        KEND1 = KMGD  + NMAIJK(ISYMGD)
        LWRK1 = LWORK - KEND1

        IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient memory in CC_BFDENF.')
        END IF
      

        DO A = 1, NBAS(ISYALP)
C
C         ---------------------------------------------------------
C         calculate one-index AO back-transformed Zeta2 amplitudes:
C         ---------------------------------------------------------
C
          IALP = A + IBAS(ISYALP)

          IOPT2AO = 1
          CALL CC_T2AO(ZETA2,XLAMDP,ISYXLP,WORK(KSCRM),WORK(KEND1),
     *                 LWRK1,IALP,ISYALP,ISYCTR,IOPT2AO)
C
C         ---------------------------------------------------------
C         transform B index with T1AM to occupied and
C         add 2 Chi(ki) Lambda(alp,j) - Chi(kj) Lambda(alp,i)
C         ---------------------------------------------------------
C
          IOPT = 4
          CALL CC_BFMO(WORK(KMGD),WORK(KSCRM),ISYMM1,T1AM,ISYTAM,
     *                 A,ISYALP,XLAMDP,ISYXLP,CHI,ISYCHI,ZERO,IOPT)
C
C         ---------------------------------------------------------
C         add one-index backtransformed M intermediate:
C         ---------------------------------------------------------
C
          CALL CC_BFMIMO(WORK(KMGD),MINT,ISYMIN,XLAMDP,ISYXLP,A,ISYALP)
C
C         -----------------------
C         sort into result array:
C         -----------------------
C
          DO ISYMJ = 1, NSYM
          DO ISYMI = 1, NSYM

             ISYMIJ = MULD2H(ISYMI,ISYMJ)
             ISYMK  = MULD2H(ISYMGD,ISYMIJ)
             ISYMKI = MULD2H(ISYMK,ISYMI)
             ISYMAK = MULD2H(ISYMK,ISYALP)

             DO J = 1, NRHF(ISYMJ)
             DO I = 1, NRHF(ISYMI)
             DO K = 1, NRHF(ISYMK)

                KOFF1 = IMAIJK(ISYMKI,ISYMJ) + NMATIJ(ISYMKI)*(J-1)
     &                + IMATIJ(ISYMK,ISYMI) + NRHF(ISYMK)*(I-1) + K

                NIJ   = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + I

                KOFF2 = IT2AOIJ(ISYMAK,ISYMIJ) + NT1AO(ISYMAK)*(NIJ-1)
     &                + IT1AO(ISYALP,ISYMK) + NBAS(ISYALP)*(K-1) + A

                MGD(KOFF2) = WORK(KMGD-1+KOFF1)

             END DO
             END DO
             END DO

          END DO
          END DO

        END DO
      END DO
 
      RETURN
      END 
*=====================================================================*
*=====================================================================*
C  /* Deck cc_bfao */
      SUBROUTINE CC_BFAO(XMGD,T2AMP,ISYMT,XLAMD0,ISYML0,D,ISYMD,
     *                   XLAMD1,ISYML1,XLAMD2,ISYML2,FACT,IOPT)
*---------------------------------------------------------------------*
*     Purpose:
*
*         Transform in a three-index batch T2AMP(ci,j) 
*         the virtual index 'c' back to AO basis
*
*         FACT = 0.0  --  initialize XMGD with result
*         FACT = 1.0  --  add result to what is in XMGD
*
*         IOPT=1,3,5 : add F-term contribution
*         IOPT=2,4,6 : add 2*Coul - Exch combination
*         (else do not add anything)
*
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
#if defined (SYS_CRAY)
      REAL XMGD(*),T2AMP(*),XLAMD0(*),XLAMD1(*),XLAMD2(*)
      REAL ONE, TWO, FACT
#else
      DOUBLE PRECISION XMGD(*),T2AMP(*),XLAMD0(*),XLAMD1(*),XLAMD2(*)
      DOUBLE PRECISION ONE, TWO, FACT
#endif
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0)
C
      INTEGER ISYMT, ISYML0, ISYMD, ISYML1, ISYML2, IOPT
      INTEGER ISYMJ,ISYMCI,ISYMDJ,ISYMC,ISYMG,ISYMGI,ISYMGJ,ISYMDI
      INTEGER KOFF1, KOFF2, KOFF3, KOFF4, KOFF5, ISYMI, NVIRC, NBASG
      INTEGER KOFF6
C
      DO ISYMJ = 1,NSYM
C
         ISYMCI = MULD2H(ISYMJ,ISYMT)
         ISYMDJ = MULD2H(ISYMD,ISYMJ)
C
         DO ISYMI = 1,NSYM
C
            ISYMC  = MULD2H(ISYMI,ISYMCI)
            ISYMG  = MULD2H(ISYMC,ISYML0)
            ISYMGI = MULD2H(ISYMG,ISYMI)
            ISYMGJ = MULD2H(ISYMG,ISYMJ)
            ISYMDI = MULD2H(ISYMD,ISYMI)
C
c           IF (IOPT.EQ.2) THEN
c            KOFF5 = IGLMRH(ISYMD,ISYMJ) + 1
c            WRITE (LUPRI,*) 'XLAMD1: block:',ISYMD,ISYMJ,KOFF5
c            CALL OUTPUT(XLAMD1(KOFF5),1,NBAS(ISYMD),1,NRHF(ISYMJ),
c    &                                   NBAS(ISYMD),  NRHF(ISYMJ),1,LUPRI)
c            KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
c            WRITE (LUPRI,*) 'XLAMD2: block:',ISYMG,ISYMI,KOFF5
c            CALL OUTPUT(XLAMD2(KOFF5),1,NBAS(ISYMG),1,NRHF(ISYMI),
c    &                                   NBAS(ISYMG),  NRHF(ISYMI),1,LUPRI)
c           END IF
C
            NVIRC = MAX(NVIR(ISYMC),1)
            NBASG = MAX(NBAS(ISYMG),1)
C
            KOFF1 = IGLMVI(ISYMG,ISYMC) + 1
C
            DO J = 1,NRHF(ISYMJ)
C
C              --------------------------------------------------
C              add T2AMP with virtual index backtransformed to AO
C              --------------------------------------------------
C
               KOFF2 = IT2BCD(ISYMCI,ISYMJ) + IT1AM(ISYMC,ISYMI)
     *               + NT1AM(ISYMCI)*(J - 1) + 1
               KOFF3 = IT2BGD(ISYMGI,ISYMJ) + IT1AO(ISYMG,ISYMI)
     *               + NT1AO(ISYMGI)*(J - 1) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC),
     *                    ONE,XLAMD0(KOFF1),NBASG,T2AMP(KOFF2),NVIRC,
     *                    FACT,XMGD(KOFF3),NBASG)
C
C              ----------------------------------------------------
C              for IOPT=1,3,5 add F term contribution to density
C              ----------------------------------------------------
C
               IF      (IOPT.EQ.1 .OR. IOPT.EQ.3 .OR. IOPT.EQ.5) THEN

                 KOFF4 = IGLMRH(ISYMG,ISYMI) + 1
                 KOFF5 = IGLMRH(ISYMD,ISYMJ) + NBAS(ISYMD)*(J-1) + D
 
                 IF (ISYML2.EQ.ISYMGI .AND. ISYML1.EQ.ISYMDJ) THEN
                   CALL DAXPY(NBAS(ISYMG)*NRHF(ISYMI),XLAMD1(KOFF5),
     *                        XLAMD2(KOFF4),1,XMGD(KOFF3),1)
                 END IF

                 IF (IOPT.EQ.3) THEN
                   IF (ISYML1.EQ.ISYMGI .AND. ISYML2.EQ.ISYMDJ) THEN
                     CALL DAXPY(NBAS(ISYMG)*NRHF(ISYMI),XLAMD2(KOFF5),
     *                          XLAMD1(KOFF4),1,XMGD(KOFF3),1)
                   END IF
                 END IF

C
C              -----------------------------------------------
C              for IOPT=2,4,6,8 add 2*Coul - Exch combination:
C              -----------------------------------------------
C
               ELSE IF (IOPT.EQ.2 .OR. IOPT.EQ.4 .OR. 
     *                  IOPT.EQ.6 .OR. IOPT.EQ.8       ) THEN

                 KOFF4 = IGLMRH(ISYMG,ISYMI) + 1
                 KOFF5 = IGLMRH(ISYMD,ISYMJ) + NBAS(ISYMD)*(J-1) + D
 
                 IF (ISYML2.EQ.ISYMGI .AND. ISYML1.EQ.ISYMDJ) THEN
                  CALL DAXPY(NBAS(ISYMG)*NRHF(ISYMI),TWO*XLAMD1(KOFF5),
     *                        XLAMD2(KOFF4),1,XMGD(KOFF3),1)
                 END IF
               
                 IF (ISYML1.EQ.ISYMGI .AND. ISYML2.EQ.ISYMDJ) THEN
                  CALL DAXPY(NBAS(ISYMG)*NRHF(ISYMI),TWO*XLAMD2(KOFF5),
     *                        XLAMD1(KOFF4),1,XMGD(KOFF3),1)
                 END IF
               
                 IF (ISYML2.EQ.ISYMGJ .AND. ISYML1.EQ.ISYMDI) THEN
                   DO I = 1, NRHF(ISYMI)

                     KOFF4 = IGLMRH(ISYMG,ISYMJ)+ NBAS(ISYMG)*(J-1)+ 1
                     KOFF5 = IGLMRH(ISYMD,ISYMI)+ NBAS(ISYMD)*(I-1)+ D
                     KOFF6 = KOFF3 + NBAS(ISYMG)*(I-1) 

                     CALL DAXPY(NBAS(ISYMG),-XLAMD1(KOFF5),
     *                          XLAMD2(KOFF4),1,XMGD(KOFF6),1)

                   END DO
                 END IF
 
                 IF (ISYML1.EQ.ISYMGJ .AND. ISYML2.EQ.ISYMDI) THEN
                   DO I = 1, NRHF(ISYMI)

                     KOFF4 = IGLMRH(ISYMG,ISYMJ)+ NBAS(ISYMG)*(J-1)+ 1
                     KOFF5 = IGLMRH(ISYMD,ISYMI)+ NBAS(ISYMD)*(I-1)+ D
                     KOFF6 = KOFF3 + NBAS(ISYMG)*(I-1) 

                     CALL DAXPY(NBAS(ISYMG),-XLAMD2(KOFF5),
     *                          XLAMD1(KOFF4),1,XMGD(KOFF6),1)

                   END DO
                 END IF
C
               END IF

            END DO
         END DO
      END DO
C
      RETURN
      END 
*=====================================================================*
*              END OF SUBROUTINE CC_BFAO                              *
*=====================================================================*
*=====================================================================*
C  /* Deck cc_bfmiao */
      SUBROUTINE CC_BFMIAO(XMGD,MINT,ISYMINT,XLAMD1,ISYXL1,
     &                     XLAMD2,ISYXL2,D,ISYMD,WORK,LWORK)
*---------------------------------------------------------------------*
*     Purpose: Transform two indeces of the M intermediate back to AO
*              and add it to the effective density XMGD
*
*              fixed D index is obtained with XLAMD1 (ISYXL1)
*              free  G index is obtained with XLAMD2 (ISYXL2)
*
* Christof Haettig, 01.09.1998
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
#if defined (SYS_CRAY)
      REAL XMGD(*), MINT(*), XLAMD1(*), XLAMD2(*), WORK(*)
      REAL ZERO, ONE, XHALF, XNORM, DDOT
#else
      DOUBLE PRECISION XMGD(*), MINT(*), XLAMD1(*), XLAMD2(*), WORK(*)
      DOUBLE PRECISION ZERO, ONE, XHALF, XNORM, DDOT
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, XHALF = 0.5D0)
 
      INTEGER ISYMD,ISYXL1,ISYXL2,ISYMINT, ISYMN, ISYKIJ, ISYMJ, ISYMKI
      INTEGER ISYKIN, ISYMI, ISYMK, KOFF1, KOFF2, NTOTKI, NTOTD, LWORK
      INTEGER KOFF4, KOFF5, KOFF6, NTOTG, NTOTK, ISYMG, ISYMGI
 
      ISYMN  = MULD2H(ISYMD,ISYXL1)
      ISYKIJ = MULD2H(ISYMN,ISYMINT)
 
      KOFF2  = IGLMRH(ISYMD,ISYMN) + D
      NTOTD  = MAX(NBAS(ISYMD),1)

      IF ( NRHF(ISYMN)*NBAS(ISYMD) .EQ. 0) RETURN
 
      DO ISYMJ = 1,NSYM
 
         ISYMKI = MULD2H(ISYKIJ,ISYMJ)
         ISYKIN = MULD2H(ISYMKI,ISYMN)
         NTOTKI = MAX(NMATIJ(ISYMKI),1)

         IF (LWORK.LT.NTOTKI) THEN
            CALL QUIT('Insufficient memory in CC_BFMIAO.')
         END IF
 
         DO J = 1, NRHF(ISYMJ)
 
             KOFF1 = I3ORHF(ISYKIN,ISYMJ) + NMAIJK(ISYKIN)*(J-1)
     *             + IMAIJK(ISYMKI,ISYMN) + 1
 
             CALL DGEMV('N',NMATIJ(ISYMKI),NRHF(ISYMN),ONE,
     *                  MINT(KOFF1),NTOTKI,XLAMD1(KOFF2),NTOTD,
     *                  ZERO,WORK,1)


             DO ISYMK = 1, NSYM
                ISYMI  = MULD2H(ISYMKI,ISYMK)
                ISYMG  = MULD2H(ISYXL2,ISYMK)
                ISYMGI = MULD2H(ISYMG,ISYMI)
 
                KOFF4 = IGLMRH(ISYMG,ISYMK) + 1
                KOFF5 = IMATIJ(ISYMK,ISYMI) + 1
                KOFF6 = IT2BGD(ISYMGI,ISYMJ) + NT1AO(ISYMGI)*(J-1) 
     *                + IT1AO(ISYMG,ISYMI) + 1

                NTOTG = MAX(NBAS(ISYMG),1)
                NTOTK = MAX(NRHF(ISYMK),1)

                CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NRHF(ISYMK),
     *                     ONE,XLAMD2(KOFF4),NTOTG,WORK(KOFF5),NTOTK,
     *                     ONE,XMGD(KOFF6),NTOTG)

             END DO
 
         END DO

      END DO
 
      RETURN
      END 
*=====================================================================*
*              END OF SUBROUTINE CC_BFMIAO                            *
*=====================================================================*
C  /* Deck cc_bfmo */
      SUBROUTINE CC_BFMO(XMGD,T2AMP,ISYMT,R1AMP,ISYMR,D,ISYMD,
     *                   XLAMD,ISYLAM,CHI,ISYCHI,FACT,IOPT)
*---------------------------------------------------------------------*
*     Purpose:
*
*      Transform in a three-index batch (ci,j) of T2 amplitudes
*      the virtual index 'c' to occupied using response R1 amplitudes
*
*      FACT=0.0 : initialize XMGD
*      FACT=1.0 : add contribution to what is already in XMGD
*
*      IOPT=1,3 : add F-term contribution
*      IOPT=2,4 : add 2*Coul - Exch combination
*      (else do not add anything)
*
* Christof Haettig, 04.07.1998
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION XMGD(*), T2AMP(*), R1AMP(*), XLAMD(*), CHI(*)
C
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
C
      DO ISYMJ = 1,NSYM
C
         ISYMCI = MULD2H(ISYMJ,ISYMT)
         ISYMDJ = MULD2H(ISYMJ,ISYMD)
C
         DO ISYMI = 1,NSYM
C
            ISYMC  = MULD2H(ISYMI,ISYMCI)
            ISYMK  = MULD2H(ISYMC,ISYMR)
            ISYMKI = MULD2H(ISYMK,ISYMI)
            ISYMKJ = MULD2H(ISYMK,ISYMJ)
            ISYMDI = MULD2H(ISYMI,ISYMD)
C
            NVIRC = MAX(NVIR(ISYMC),1)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            KOFF1 = IT1AM(ISYMC,ISYMK) + 1
C
            DO J = 1,NRHF(ISYMJ)
C
C              ----------------------------------------------------
C              add T2AMP with virtual index transformed to occupied
C              ----------------------------------------------------
C
               KOFF2 = IT2BCD(ISYMCI,ISYMJ) + NT1AM(ISYMCI)*(J - 1)
     *               + IT1AM(ISYMC,ISYMI) + 1
               KOFF3 = IMAIJK(ISYMKI,ISYMJ) + NMATIJ(ISYMKI)*(J - 1) 
     *               + IMATIJ(ISYMK,ISYMI) + 1
C
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     *                    -ONE,R1AMP(KOFF1),NVIRC,T2AMP(KOFF2),NVIRC,
     *                    FACT,XMGD(KOFF3),NRHFK)
C
C              ----------------------------------------------------
C              for IOPT=1,3 add F term contribution to density
C              ----------------------------------------------------
C
               IF      (IOPT.EQ.1 .OR. IOPT.EQ.3) THEN

                 KOFF4 = IMATIJ(ISYMK,ISYMI) + 1
                 KOFF5 = IGLMRH(ISYMD,ISYMJ) + NBAS(ISYMD)*(J-1) + D
 
                 IF (ISYCHI.EQ.ISYMKI .AND. ISYLAM.EQ.ISYMDJ) THEN
                   CALL DAXPY(NRHF(ISYMK)*NRHF(ISYMI),XLAMD(KOFF5),
     *                        CHI(KOFF4),1,XMGD(KOFF3),1)
                 END IF

C
C              -------------------------------------------
C              for IOPT=2,4 add 2*Coul - Exch combination:
C              -------------------------------------------
C
               ELSE IF (IOPT.EQ.2 .OR. IOPT.EQ.4) THEN

                 KOFF4 = IMATIJ(ISYMK,ISYMI) + 1
                 KOFF5 = IGLMRH(ISYMD,ISYMJ) + NBAS(ISYMD)*(J-1) + D
 
                 IF (ISYCHI.EQ.ISYMKI .AND. ISYLAM.EQ.ISYMDJ) THEN
                   CALL DAXPY(NRHF(ISYMK)*NRHF(ISYMI),TWO*XLAMD(KOFF5),
     *                        CHI(KOFF4),1,XMGD(KOFF3),1)
                 END IF
               
                 IF (ISYCHI.EQ.ISYMKJ .AND. ISYLAM.EQ.ISYMDI) THEN
                   DO I = 1, NRHF(ISYMI)

                     KOFF4 = IMATIJ(ISYMK,ISYMJ)+ NRHF(ISYMK)*(J-1)+ 1
                     KOFF5 = IGLMRH(ISYMD,ISYMI)+ NBAS(ISYMD)*(I-1)+ D
                     KOFF6 = KOFF3 + NRHF(ISYMK)*(I-1) 

                      CALL DAXPY(NRHF(ISYMK),-XLAMD(KOFF5),
     *                           CHI(KOFF4),1,XMGD(KOFF6),1)

                   END DO
                 END IF
C
               END IF

            END DO
         END DO
      END DO
C
      RETURN
      END 
*=====================================================================*
*              END OF SUBROUTINE CC_BFMO                              *
*=====================================================================*
C  /* Deck cc_bfmimo */
      SUBROUTINE CC_BFMIMO(XMGD,MINT,ISYMINT,XLAMD,ISYML,D,ISYMD)
*---------------------------------------------------------------------*
*     Purpose: Transform one index of the M intermediate back to AO
*              and add it to the effective density XMGD
*
* Christof Haettig, 05.07.1998
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
#if defined (SYS_CRAY)
      REAL XMGD(*), MINT(*), XLAMD(*)
      REAL ONE, XHALF, XNORM, DDOT
#else
      DOUBLE PRECISION XMGD(*), MINT(*), XLAMD(*)
      DOUBLE PRECISION ONE, XHALF, XNORM, DDOT
#endif
      PARAMETER (ONE = 1.0D0, XHALF = 0.5D0)
 
      INTEGER ISYMD, ISYML, ISYMINT, ISYMN, ISYKIJ, ISYMJ, ISYMKI
      INTEGER ISYKIN, ISYMI, ISYMK, KOFF1, KOFF2, KOFF3, NTOTKI, NTOTD
 
      ISYMN  = MULD2H(ISYMD,ISYML)
      ISYKIJ = MULD2H(ISYMN,ISYMINT)
 
      KOFF2  = IGLMRH(ISYMD,ISYMN) + D
      NTOTD  = MAX(NBAS(ISYMD),1)
 
      DO ISYMJ = 1,NSYM
 
         ISYMKI = MULD2H(ISYKIJ,ISYMJ)
         ISYKIN = MULD2H(ISYMKI,ISYMN)
 
         NTOTKI = MAX(NMATIJ(ISYMKI),1)
 
         DO J = 1, NRHF(ISYMJ)
 
             KOFF1 = I3ORHF(ISYKIN,ISYMJ) + NMAIJK(ISYKIN)*(J-1)
     *             + IMAIJK(ISYMKI,ISYMN) + 1
 
             KOFF3 = IMAIJK(ISYMKI,ISYMJ) + NMATIJ(ISYMKI)*(J-1) + 1
 
             CALL DGEMV('N',NMATIJ(ISYMKI),NRHF(ISYMN),XHALF,
     *                  MINT(KOFF1),NTOTKI,XLAMD(KOFF2),NTOTD,
     *                  ONE,XMGD(KOFF3),1)
 
         END DO

      END DO
 
      RETURN
      END 
*=====================================================================*
*              END OF SUBROUTINE CC_BFMIMO                            *
*=====================================================================*
*======================================================================
C  /* Deck cc_bfi */
      SUBROUTINE CC_BFI(BFI,XINT,ISYDIS,XMGD,ISYMGD,D,ISYMD,
     &                  SQRINT,FACTOR,WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: contract effective density with AO integrals to 
*              the BF intermediate
*
*     BFI    :  BF intermediate (symmetry MULD2H(ISYDIS,ISYMGD))
*     XINT   :  delta distribution of AO integrals (symmetry ISYDIS)
*     XMGD   :  effective density for BF term (symmetry ISYMGD)
*     D      :  delta index (symmetry of the orbital ISYMD)
*     FACTOR :  factor with which the contribution is multiplied
*               before it is added to what is already in BFI
*
*     Christof Haettig, September 1998
*     based on Henrik Kochs CC_BF1 routine
*     Sonia Coriani, October 1999
*     -- option to handle squared XINT distributions in input 
*        (SQRINT = true)
*
*======================================================================
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
      LOGICAL LO3BF, SQRINT
      INTEGER IOPT, LWORK, ISYDIS, ISYMGD, ISYMD, IOPTSQ
      DIMENSION XINT(*), BFI(*), XMGD(*), WORK(LWORK)
      PARAMETER (LO3BF = .FALSE., IOPT = 0)
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
C
      ISYRES = MULD2H(ISYDIS,ISYMGD)
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYRES)
         ISYMG  = MULD2H(ISYMAB,ISYDIS)
C
         KSCRAB = 1
         KINDV1 = KSCRAB + N2BST(ISYMAB)
         KINDV2 = KINDV1 + (NNBST(ISYMAB) - 1)/IRAT + 1
         KEND1  = KINDV2 + (NNBST(ISYMAB) - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space in CC_BFI')
         ENDIF
C
C        ------------------------
C        Calculate index vectors:
C        ------------------------
C
         CALL CCSD_INDEX(WORK(KINDV1),WORK(KINDV2),ISYMAB)
C
C        -----------------------------------
C        Batching and work space allocation.
C        -----------------------------------
C
         IF (ISYMG .EQ. ISYMD) THEN
            IMAXG = D
         ELSE IF (ISYMG .LT. ISYMD) THEN
            IMAXG = NBAS(ISYMG)
         ELSE
            IMAXG = 0
         ENDIF
C
         NSIZE  = 2*(NNBST(ISYMAB) + NMIJP(ISYMIJ))
C
         IF (IMAXG.EQ.0 .OR. NSIZE.EQ.0) GOTO 100
C
         IF (LWRK1.LT.NSIZE) THEN
            CALL QUIT('Insufficient memory in CC_BFI.')
         END IF
C
         NMAXG  = MIN(IMAXG,LWRK1/NSIZE)
C
         NBATCH = (IMAXG - 1)/NMAXG + 1
C
         DO 110 IBATCH = 1,NBATCH
C
            NUMG = NMAXG
            IF (IBATCH .EQ. NBATCH) THEN
               NUMG = IMAXG - NMAXG*(NBATCH - 1)
            ENDIF
C
            IG1 = NMAXG*(IBATCH - 1) + 1
            IG2 = NMAXG*(IBATCH - 1) + NUMG
C
            KINTP = KEND1
            KINTM = KINTP + NNBST(ISYMAB)*NUMG
            KT2MP = KINTM + NNBST(ISYMAB)*NUMG
            KT2MM = KT2MP + NUMG*NMIJP(ISYMIJ)
            KEND2 = KT2MM + NUMG*NMIJP(ISYMIJ)
C
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient space in CC_BFI')
            ENDIF
C
C           ------------------------
C           Construct T2MP and T2MM:
C           ------------------------
C
            CALL CCRHS_TPM(WORK(KT2MP),WORK(KT2MM),XMGD,ISYMIJ,
     &                     IG1,NUMG,ISYMG,IOPT,LO3BF)
C
C           ------------------------
C           Construct INTP and INTM:
C           ------------------------
C
            IF (SQRINT) THEN
              IOPTSQ = 1
            ELSE
              IOPTSQ = 0
            END IF

            CALL CCRHS_IPM1(XINT,WORK(KINTP),WORK(KINTM),WORK(KSCRAB),
     *                  WORK(KINDV1),WORK(KINDV2),ISYMAB,ISYMG,
     *                  NUMG,IG1,IG2,IOPTSQ)
C
C           --------------------
C           Scale the diagonals:
C           --------------------
C
            IF ((ISYMG .EQ. ISYMD) .AND. (IBATCH .EQ. NBATCH)) THEN
               KOFF = KINTP + NNBST(ISYMAB)*(NUMG - 1)
               CALL DSCAL(NNBST(ISYMAB),HALF,WORK(KOFF),1)
            ENDIF
C
C           -----------------------------
C           Add the B-term contributions:
C           -----------------------------
C
            KOFF1  = IT2ORT(ISYMAB,ISYMIJ) + 1
            KOFF2  = NT2ORT(ISYRES) + IT2ORT(ISYMAB,ISYMIJ) + 1
C
            NUMGM  = MAX(NUMG,1)
            NTOTAB = MAX(NNBST(ISYMAB),1)
C
            CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                 FACTOR,WORK(KINTP),NTOTAB,WORK(KT2MP),NUMGM,
     *                 ONE,BFI(KOFF1),NTOTAB)
C
            CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                 FACTOR,WORK(KINTM),NTOTAB,WORK(KT2MM),NUMGM,
     *                 ONE,BFI(KOFF2),NTOTAB)
C
  110   CONTINUE
C
  100 CONTINUE
C
      RETURN
C
      END
*=====================================================================*
*=====================================================================*
C  /* Deck cc_bfib */
      SUBROUTINE CC_BFIB(BFI,BSRHF,ISYRHF,XMGD,ISYMGD,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: contract effective density with (**|k delta) integrals
*              to the BF intermediate in the B matrix transformation
*              (special version of the CC_BFI routine for B matrix)
*
*     BFI    : BF(alp i,kj) intermediate (symmetry ISYRHF x ISYMGD)
*     BSRHF  : (**|k delta) integral distribution (symmetry ISYRHF)
*     XMGD   : effective density for BF term (symmetry ISYMGD)
*
*     Written by Christof Haettig July/October 1998
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>

      INTEGER LWORK, ISYMGD, ISYRHF

#if defined (SYS_CRAY)
      REAL BSRHF(*), BFI(*), XMGD(*), WORK(LWORK)
      REAL ONE, HALF
#else
      DOUBLE PRECISION BSRHF(*), BFI(*), XMGD(*), WORK(LWORK)
      DOUBLE PRECISION ONE, HALF
#endif
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0)
C
      INTEGER ISYRES, ISYMIJ, ISYMAK, ISYMG
      INTEGER ISYMI, ISYMJ, ISYMGI, ISYM, ISYBET, ICOUNT
      INTEGER KT2M, KEND1, LWRK1, KOFF1, KOFF3
      INTEGER NGIJ, NIJ, NTOTAK, NBASG
      INTEGER NBSRHF(8), IBSRHF(8,8)
C
C     --------------------------------------
C     precalculate symmetry array for BSRHF:
C     --------------------------------------
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO
C
      ISYRES = MULD2H(ISYRHF,ISYMGD)
C
      DO ISYMIJ = 1, NSYM
C
         ISYMAK = MULD2H(ISYRES,ISYMIJ)
         ISYMG  = MULD2H(ISYMAK,ISYRHF)
C
         KT2M  = 1
         KEND1 = KT2M + NMATIJ(ISYMIJ)*NBAS(ISYMG)
         LWRK1 = LWORK - KEND1
          
         IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_BFIB.')
         END IF
C
C        ---------------------------------------------------
C        resort amplitude density:
C        ---------------------------------------------------
C
         DO ISYMJ = 1, NSYM
C
            ISYMI  = MULD2H(ISYMJ,ISYMIJ)
            ISYMGI = MULD2H(ISYMI,ISYMG)
C
            DO J = 1, NRHF(ISYMJ)
            DO I = 1, NRHF(ISYMI)
C
               NGIJ = IT2BGD(ISYMGI,ISYMJ) + NT1AO(ISYMGI)*(J-1) 
     &              + IT1AO(ISYMG,ISYMI) + NBAS(ISYMG)*(I-1) + 1

               NIJ  = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + I

               KOFF1 = KT2M + NBAS(ISYMG)*(NIJ-1)

               CALL DCOPY(NBAS(ISYMG),XMGD(NGIJ),1,WORK(KOFF1),1)

            END DO
            END DO

         END DO
C
C        --------------------------------------
C        contract with the presorted integrals:
C        --------------------------------------
C
         KOFF1  = IT2AOIJ(ISYMAK,ISYMIJ) + 1
         KOFF3  = IBSRHF(ISYMAK,ISYMG)   + 1
         NTOTAK = MAX(NT1AO(ISYMAK),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NT1AO(ISYMAK),NMATIJ(ISYMIJ),NBAS(ISYMG),
     &              ONE,BSRHF(KOFF3),NTOTAK,WORK(KT2M),NBASG,
     &              ONE,BFI(KOFF1),  NTOTAK)
C
      END DO
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck cc_bfbsort */
      SUBROUTINE CC_BFBSORT(DSRHF,BSRHF,ISYRHF)
*---------------------------------------------------------------------*
*
*     Purpose: presort DSRHF integral array for the BF intermediate
*              calculation in the B matrix transformation
*
*     DSRHF  : (alp bet|k delta) integrals for a fixed delta
*     BSRHF  : integrals sorted as I(alp k;bet)^del
*     ISYRHF : symmetry of the integral arrays
*
*     Written by Christof Haettig July/October 1998
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <symsq.h>

      INTEGER ISYRHF, ISYM, ISYMAK, ISYBET, ISYMK, ISYMAB, ISYALP
      INTEGER ICOUNT, NBSRHF(8), IBSRHF(8,8)
      INTEGER NABK, NAKB, NAK, KOFF1, IJSQ

#if defined (SYS_CRAY)
      REAL DSRHF(*), BSRHF(*)
#else
      DOUBLE PRECISION DSRHF(*), BSRHF(*)
#endif
C
C     --------------------------------------
C     precalculate symmetry array for BSRHF:
C     --------------------------------------
C
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO
C
C     -------------------
C     sort the integrals:
C     -------------------
C
      DO ISYMAK = 1, NSYM
      DO ISYMK  = 1, NSYM
C 
         ISYBET = MULD2H(ISYMAK,ISYRHF)
         ISYALP = MULD2H(ISYMK,ISYMAK)
         ISYMAB = MULD2H(ISYALP,ISYBET)
C
C        --------------------------------------------------------
C        get (alp k;bet) blocks out of (alp bet|k del) integrals:
C        --------------------------------------------------------
C
         DO K = 1, NRHF(ISYMK)
C
            KOFF1  = IDSRHF(ISYMAB,ISYMK) + NNBST(ISYMAB)*(K-1) 
C
            DO A = 1, NBAS(ISYALP)
            DO B = 1, NBAS(ISYBET)
C
               IJSQ = IAODIS(ISYALP,ISYBET) + NBAS(ISYALP)*(B-1) + A
               NABK = KOFF1  + IADRPK( I2BST(ISYMAB) + IJSQ )
               NAK  = IT1AO(ISYALP,ISYMK)   + NBAS(ISYALP)*(K-1) + A
               NAKB = IBSRHF(ISYMAK,ISYBET) +NT1AO(ISYMAK)*(B-1) + NAK
C
               BSRHF(NAKB) = DSRHF(NABK)
C
            END DO
            END DO
C
         END DO
C
      END DO
      END DO
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
C  /* Deck ccrhs_tpm */
      SUBROUTINE CCRHS_TPM(XMP,XMM,XMGD,ISYMIJ,
     &                     IG1,NUMG,ISYMG,IOPT,LO3BF)
*---------------------------------------------------------------------*
*
*     Purpose: construct M+ and M- 'densities' for a batch of G 
*
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      PARAMETER(ONE = 1.0D0, THREE = 3.0D0)
C
      LOGICAL LO3BF
      DIMENSION XMP(*), XMM(*), XMGD(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO ISYMJ = 1,NSYM
C
         ISYMI  = MULD2H(ISYMJ,ISYMIJ)
         ISYMGI = MULD2H(ISYMI,ISYMG)
         ISYMGJ = MULD2H(ISYMJ,ISYMG)
C
         IF (ISYMI .LE. ISYMJ) THEN
C
           NTOTI = NRHF(ISYMI)
C
           DO J = 1,NRHF(ISYMJ)
C
              IF (ISYMI .EQ. ISYMJ) NTOTI = J
C
              DO I = 1,NTOTI
C
                 IF (LO3BF .AND. IOPT.EQ.4) THEN
                   NGIJ = IMAIJK(ISYMGI,ISYMJ) + NMATIJ(ISYMGI)*(J - 1)
     *                  + IMATIJ(ISYMG,ISYMI) + NRHF(ISYMG)*(I-1) + IG1
C
                   NGJI = IMAIJK(ISYMGJ,ISYMI) + NMATIJ(ISYMGJ)*(I - 1)
     *                  + IMATIJ(ISYMG,ISYMJ) + NRHF(ISYMG)*(J-1) + IG1
                 ELSE
                   NGIJ = IT2BGD(ISYMGI,ISYMJ) + NT1AO(ISYMGI)*(J - 1)
     *                  + IT1AO(ISYMG,ISYMI) + NBAS(ISYMG)*(I-1) + IG1
C
                   NGJI = IT2BGD(ISYMGJ,ISYMI) + NT1AO(ISYMGJ)*(I - 1)
     *                  + IT1AO(ISYMG,ISYMJ) + NBAS(ISYMG)*(J-1) + IG1
                 END IF
C
                 IF (ISYMI .EQ. ISYMJ) THEN
                    NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                 ELSE
                    NIJ = IMIJP(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                 ENDIF
C
                 NGIJPM = NUMG*(NIJ - 1) + 1
C
                 IF ( CC2 ) THEN
                    CALL DZERO(XMP(NGIJPM),NUMG)
                    CALL DZERO(XMM(NGIJPM),NUMG)
                 ELSE
                    CALL DCOPY(NUMG,XMGD(NGIJ),1,XMP(NGIJPM),1)
                    CALL DCOPY(NUMG,XMGD(NGIJ),1,XMM(NGIJPM),1)
C
                    CALL DAXPY(NUMG, ONE,XMGD(NGJI),1,XMP(NGIJPM),1)
                    CALL DAXPY(NUMG,-ONE,XMGD(NGJI),1,XMM(NGIJPM),1)
                 ENDIF
 
              END DO

           END DO
 
         END IF

      END DO
C
      RETURN
      END 
*=====================================================================*
*              END OF SUBROUTINE CCRHS_TPM                            *
*=====================================================================*
C  /* Deck cc_bfmf */
      SUBROUTINE CC_BFMF(XMP,XMM,NUMG,FACT,IOPT,LO3BF,
     &                   I,ISYMI,J,ISYMJ,IG1,ISYMG,D,ISYMD,
     &                   XLAMD1,ISYML1,XLAMD2,ISYML2)
*---------------------------------------------------------------------*
*
*     Purpose: add F-term 'density' to M+ and M- arrays
*
*---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      LOGICAL LO3BF
      DIMENSION XMP(*), XMM(*), XLAMD1(*), XLAMD2(*)
C
      IF ((ISYMJ .EQ. MULD2H(ISYML2,ISYMD)) .AND.
     &    (ISYMI .EQ. MULD2H(ISYML1,ISYMG))       ) THEN
C
          KOFF1 = IGLMRH(ISYMD,ISYMJ) + NBAS(ISYMD)*(J - 1) + D
C
          IF (LO3BF .AND. IOPT.EQ.4) THEN
             KOFF2 = IMATIJ(ISYMG,ISYMI) + NRHF(ISYMG)*(I - 1) + IG1
          ELSE
             KOFF2 = IGLMRH(ISYMG,ISYMI) + NBAS(ISYMG)*(I - 1) + IG1
          END IF
C
          CALL DAXPY(NUMG,     XLAMD2(KOFF1),XLAMD1(KOFF2),1,XMP,1)
          CALL DAXPY(NUMG,FACT*XLAMD2(KOFF1),XLAMD1(KOFF2),1,XMM,1)
C
      ENDIF
C
      IF ((ISYMI .EQ. MULD2H(ISYML2,ISYMD)) .AND.
     &    (ISYMJ .EQ. MULD2H(ISYML1,ISYMG))       ) THEN
C
          KOFF1 = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1) + D
C 
          IF (LO3BF .AND. IOPT.EQ.4) THEN
             KOFF2 = IMATIJ(ISYMG,ISYMJ) + NRHF(ISYMG)*(J - 1) + IG1
          ELSE
             KOFF2 = IGLMRH(ISYMG,ISYMJ) + NBAS(ISYMG)*(J - 1) + IG1
          END IF
C
          CALL DAXPY(NUMG,      XLAMD2(KOFF1),XLAMD1(KOFF2),1,XMP,1)
          CALL DAXPY(NUMG,-FACT*XLAMD2(KOFF1),XLAMD1(KOFF2),1,XMM,1)
C
      ENDIF
C
      RETURN
      END 
*=====================================================================*
*              END OF SUBROUTINE CC_BFMF                              *
*=====================================================================*
