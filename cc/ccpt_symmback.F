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
C  /* Deck cc_symmback */
      SUBROUTINE SYMMBACK(MODEL,
     &                    LUIAJK,FNDIAJK,
     &                    LUAIJK,FNDAIJK,
     &                    LUABI1,FNDABI1,
     &                    LUABI3,FNDABI3,
     &                    LUPTIAJB,FNDIAJB,
     &                    ISYRES,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C     Purpose: construct/resort the symmetrized (T) densities for CCSD(T)
C              resort the other (T) densities
C              and drive backtranformation of last index to AO
C
C     S. Coriani, December 2001/January 2002
C-----------------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>

      INTEGER LUPTIAJB, LUAIJK, LUIAJK, LUABI1, LUABI3, LWORK
      INTEGER LUIJKAS, LUIAJKS, LUAIJKS, LUABICS, LUAIBCS, LUIABCS
      INTEGER LUIJKDS, LUIAJDS, LUAIJDS, LUABIDS
      INTEGER LUIAJDE, LUAIJDE, LUAIBDE, LUIABDE
      INTEGER ISYRES, KSTART
      INTEGER KOCC1,KOCC2,KOCCR1,KOCCR2,KOCCS1,LWRK1,KEND1
      INTEGER IOFF, ISYMA, ISYMI, ISYMJ, ISYMK
      INTEGER ISYKJI,ISYMKJ
      integer ISYMIJ, ISYIJK, KIJKA, KIJAK, KIJKAS
      integer LENGTH, KAIB, KBIA, ISYM, KABIS
      integer KDEN4, KDEN3
      integer ISYBIA, ISYAIB, KAIB_C, KBIA_C, KABIS_C, ISYMB
      integer ISYBIC, ISYCIB
      integer KBICS_A 
      integer KOFF2, KOFF3, KOFF1, ISYMAI, ISYMC, ISYMBI
      integer ISYMCI
      integer IOPT, IOPTX
      integer kaijk, kiajk, isyaij, isyden, isymab, isyabi
      integer kcmo,kend0,lwrk0,koff
      integer knoddy1,knoddy2,kaibd,kbida
      integer kdaijb, kdaibj, kscrpk, isymbj, nai, naij, nbj
      integer kaibj, kaijb, kiajb,kdiajb
      integer kd1kjia, kd2jkia, kdiajk, kdaijk,isycmo
      integer kbicr4,ioffbic_a,ioffcib_a

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL ZERO, ONE, TWO, FACT
      REAL DDOT,XNORM2,XNORM,XNABI1,XNABI3
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO, FACT
      DOUBLE PRECISION DDOT,XNORM2,XNORM,XNABI1,XNABI3
#endif


      CHARACTER MODEL*10
      CHARACTER MODELPRI*4, MODELPRI2*14
      CHARACTER*(*) FNDIAJK, FNDAIJK, FNDABI1, FNDABI3, FNDIAJB
      LOGICAL LOCDBG
C
C Locally define density file names for easier debug
C
      CHARACTER FNDSIAJK*9, FNDSAIJK*9, FNDSIJKA*9
      CHARACTER FNDSABIC*9, FNDSAIBC*9, FNDSIABC*9
      CHARACTER FNDSIJKD*9, FNDSIAJD*9, FNDSAIJD*9, FNDSABID*9
      CHARACTER FNDSAIBD*9, FNDSIABD*9
      PARAMETER (FNDSIAJK = 'PT_DSIAJK', FNDSAIJK = 'PT_DSAIJK')
      PARAMETER (FNDSIJKA = 'PT_DSIJKA', FNDSABIC = 'PT_DSABIC')
      PARAMETER (FNDSAIBC = 'PT_DSAIBC', FNDSIABC = 'PT_DSIABC')
c
      PARAMETER (FNDSIJKD = 'PT_DSIJKD', FNDSIAJD = 'PT_DSIAJD')
      PARAMETER (FNDSAIJD = 'PT_DSAIJD', FNDSABID = 'PT_DSABID')
      PARAMETER (FNDSAIBD = 'PT_DSAIBD', FNDSIABD = 'PT_DSIABD')
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      PARAMETER (LOCDBG = .FALSE.)
C
C#include <leinf.h>
C
      CALL QENTER('SYMMBACK')
C
C-------------------------------
C     First allocation
C-------------------------------
C
      KCMO   = 1
      KEND0 = KCMO   + NLAMDS
      LWRK0 = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND0
         CALL QUIT(
     *         'Insufficient memory for initial allocation')
      ENDIF
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE (LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND0),LWRK0)
C
C---------------------------------------------------------
C     Open files for backtransformed/symmetrized densities
C---------------------------------------------------------
C
      LUIJKAS  = -1
      LUIJKDS  = -1
      LUABICS  = -1
      LUABIDS  = -1
C
C Temporarily keep the other "S" files with the resorted densities
C (instead of symmetrized)
C 
      LUIAJKS = -1
      LUAIJKS = -1
      LUAIBCS = -1
      LUIABCS = -1
C
      LUIAJDE  = -1
      LUAIJDE  = -1
      LUAIBDE  = -1
      LUIABDE  = -1
C
C Symmetrized and backtransformed
C
C     d^s_{ijka}
      CALL WOPEN2(LUIJKAS,FNDSIJKA,64,0)
C     d^s_{ijkdelta}
      CALL WOPEN2(LUIJKDS,FNDSIJKD,64,0)
C     d^s_{abic}
      CALL WOPEN2(LUABICS,FNDSABIC,64,0)
C     d^s_{abidelta}
      CALL WOPEN2(LUABIDS,FNDSABID,64,0)
C
C Only backtransformed
C
C     d^s_{iajdelta}
      CALL WOPEN2(LUIAJDE,FNDSIAJD,64,0)
C     d^s_{aijdelta}
      CALL WOPEN2(LUAIJDE,FNDSAIJD,64,0)
C     d^s_{aibdelta}
      CALL WOPEN2(LUAIBDE,FNDSAIBD,64,0)
C     d^s_{iabdelta}
      CALL WOPEN2(LUIABDE,FNDSIABD,64,0)
C
C Only symmetrized
C
C     d^s_{aibc}
      CALL WOPEN2(LUAIBCS,FNDSAIBC,64,0)
C     d^s_{iabc}
      CALL WOPEN2(LUIABCS,FNDSIABC,64,0)
C
C---------------------------------------------------------------------
C - Read in iajk, aijk, iajb (the latter read-in inside backtransform)
C - Resort iajk, aijk
C - Backtransform last index and save on file
C---------------------------------------------------------------------

      KSTART = KEND0
      KOCC1  = KSTART
      KOCC2  = KOCC1 + NCKIJ(ISYRES)
      KOCCR1 = KOCC2 + NCKIJ(ISYRES)  !resorted density 1 diajk
      KOCCR2 = KOCCR1 + NCKIJ(ISYRES) !resorted density 2 daijk
      KEND1  = KOCCR2 + NCKIJ(ISYRES)
      LWRK1  = LWORK - KEND1
C
C-----------
C Initialize
C-----------
C
      CALL DZERO(WORK(KOCC1),NCKIJ(ISYRES))
      CALL DZERO(WORK(KOCC2),NCKIJ(ISYRES))
      CALL DZERO(WORK(KOCCR1),NCKIJ(ISYRES))
      CALL DZERO(WORK(KOCCR2),NCKIJ(ISYRES))
C
C--------------------------------------------------------------
C Read in the entire diajk(kjia) and daijk(jkia) distributions, 
C and keep in WORK(KOCC1) and WORK(KOCC2) for later use
C--------------------------------------------------------------
C
      IOFF = 1
      CALL GETWA2(LUIAJK,FNDIAJK,WORK(KOCC1),IOFF,NCKIJ(ISYRES))
C
!OBS: scale -1 
C
      CALL DSCAL(NCKIJ(ISYRES),-ONE,WORK(KOCC1),1)

      CALL GETWA2(LUAIJK,FNDAIJK,WORK(KOCC2),IOFF,NCKIJ(ISYRES))
C
!OBS: scale -1 
C
      CALL DSCAL(NCKIJ(ISYRES),-ONE,WORK(KOCC2),1)

C
C------------------------------------
C Resort diajk(kjia) to diajk(ai,j,k)
C Resort daijk(jkia) to daijk(ai,j,k)
C------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYKJI = MULD2H(ISYRES,ISYMA)
         DO A = 1, NVIR(ISYMA)
            DO ISYMI = 1, NSYM
               ISYMKJ = MULD2H(ISYKJI,ISYMI)
               DO I = 1, NRHF(ISYMI)
                  DO ISYMJ = 1, NSYM
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYAIJ = MULD2H(ISYMAI,ISYMJ)
                     ISYMK = MULD2H(ISYMKJ,ISYMJ)
                     DO J = 1, NRHF(ISYMJ)
                     DO K = 1, NRHF(ISYMK)
                        KD1KJIA = I3OVIR(ISYKJI,ISYMA) 
     &                          + NMAIJK(ISYKJI)*(A-1) 
     &                          + IMAIJK(ISYMKJ,ISYMI)
     &                          + NMATIJ(ISYMKJ)*(I-1)
     &                          + IMATIJ(ISYMK,ISYMJ)
     &                          + NRHF(ISYMK)*(J-1)
     &                          + K 
     &                          + KOCC1 - 1

                        KD2JKIA = I3OVIR(ISYKJI,ISYMA) 
     &                          + NMAIJK(ISYKJI)*(A-1) 
     &                          + IMAIJK(ISYMKJ,ISYMI)
     &                          + NMATIJ(ISYMKJ)*(I-1)
     &                          + IMATIJ(ISYMJ,ISYMK)
     &                          + NRHF(ISYMJ)*(K-1)
     &                          + J 
     &                          + KOCC2 - 1

                        KDAIJK = ICKITR(ISYAIJ,ISYMK)
     &                          + NCKI(ISYAIJ)*(K-1)
     &                          + ICKI(ISYMAI,ISYMJ)
     &                          + NT1AM(ISYMAI)*(J-1)
     &                          + IT1AM(ISYMA,ISYMI)
     &                          + NVIR(ISYMA)*(I-1) + A

                        KDIAJK = KDAIJK               !d_iajk is stored aijk

                        WORK(KDIAJK + KOCCR1 - 1) = WORK(KD1KJIA)
                        WORK(KDAIJK + KOCCR2 - 1) = WORK(KD2JKIA)

                     END DO         !K
                     END DO         !J
                  END DO            !ISYMJ
               END DO               !I
            END DO                  !ISYMI
         END DO                     !A
      END DO                        !ISYMA
C
C--------------------------------------------------------
C     Read-in d_iajb(ai,bj), unpack and reorder to ai,j,b     
C--------------------------------------------------------
C
      KDIAJB = KEND1
      KDAIBJ = KDIAJB + NT2SQ(ISYRES)
      KSCRPK = KDAIBJ + NT2SQ(ISYRES)
      KEND1  = KSCRPK + NT2AM(ISYRES)
      LWRK1  = LWORK - KEND1

      CALL DZERO(WORK(KDIAJB),NT2SQ(ISYRES))
      CALL DZERO(WORK(KDAIBJ),NT2SQ(ISYRES))

      IF (NT2AM(ISYRES) .GT. 0) THEN
         IOFF = 1
         CALL GETWA2(LUPTIAJB,FNDIAJB,WORK(KSCRPK),
     &               IOFF,NT2AM(ISYRES))
         !square up
         CALL CC_T2SQ(WORK(KSCRPK),WORK(KDAIBJ),1)
         !Husk: scale by 2*
         CALL DSCAL(NT2SQ(ISYRES),TWO,WORK(KDAIBJ),1)
      ENDIF
 
      CALL CC3_T2TP(WORK(KDIAJB),WORK(KDAIBJ),ISYRES)
C
C-----------------------------------------------
C Backtransform
C-----------------------------------------------
C
      IOPTX=2
      CALL BTRIAJDEL(IOPTX,WORK(KCMO),WORK(KOCCR1),WORK(KDIAJB),
     &               ISYRES,LUIAJDE,FNDSIAJD,WORK(KEND1),LWRK1)

      IOPTX=1
      CALL BTRIAJDEL(IOPTX,WORK(KCMO),WORK(KOCCR2),WORK(KDIAJB),
     &           ISYRES,LUAIJDE,FNDSAIJD,WORK(KEND1),LWRK1)
C
C-------------------------------------------------------------------
C STEP A.2
C Construct symmetrized occ.occ.occ.vir and backtransform last index
C d^s_ijka (ijka) = d2_ijak(ijka) + d1_ijka(jika)
C-------------------------------------------------------------------
C
      KOCCS1 = KOCCR1
      KEND1  = KOCCS1 + NCKIJ(ISYRES)
      LWRK1  = LWORK - KEND1
      
      CALL DZERO(WORK(KOCCS1),NCKIJ(ISYRES))

      DO ISYMA = 1, NSYM
         ISYIJK = MULD2H(ISYRES,ISYMA)
         DO A = 1, NVIR(ISYMA)
            DO ISYMK = 1, NSYM
               ISYMIJ = MULD2H(ISYIJK,ISYMK)
               DO K = 1, NRHF(ISYMK)
                  DO ISYMJ = 1, NSYM
                     ISYMI = MULD2H(ISYMIJ,ISYMJ)
                     DO J = 1, NRHF(ISYMJ)
                     DO I = 1, NRHF(ISYMI)

                        KIJAK = I3OVIR(ISYIJK,ISYMA) 
     &                        + NMAIJK(ISYIJK)*(A-1) 
     &                        + IMAIJK(ISYMIJ,ISYMK)
     &                        + NMATIJ(ISYMIJ)*(K-1) 
     &                        + IMATIJ(ISYMI,ISYMJ)
     &                        + NRHF(ISYMI)*(J-1) 
     &                        + I
     &                        + KOCC2 - 1


                        KIJKA = I3OVIR(ISYIJK,ISYMA) 
     &                        + NMAIJK(ISYIJK)*(A-1) 
     &                        + IMAIJK(ISYMIJ,ISYMK)
     &                        + NMATIJ(ISYMIJ)*(K-1) 
     &                        + IMATIJ(ISYMJ,ISYMI)
     &                        + NRHF(ISYMJ)*(I-1) 
     &                        + J
     &                        + KOCC1 - 1

                        KIJKAS = I3OVIR(ISYIJK,ISYMA)
     &                         + NMAIJK(ISYIJK)*(A-1) 
     &                         + IMAIJK(ISYMIJ,ISYMK)
     &                         + NMATIJ(ISYMIJ)*(K-1) 
     &                         + IMATIJ(ISYMI,ISYMJ)
     &                         + NRHF(ISYMI)*(J-1) 
     &                         + I
     &                         + KOCCS1 - 1

                        WORK(KIJKAS) = WORK(KIJAK) + WORK(KIJKA)

                     END DO         !I
                     END DO         !J
                  END DO            !ISYMJ
               END DO               !K
            END DO                  !ISYMK
         END DO                     !A
      END DO                        !ISYMA
C
C--------------
C Backtransform
C--------------
C
      if (locdbg) then
         xnorm = ddot(NCKIJ(ISYRES),WORK(KOCC1),1,WORK(KOCC1),1)
         write(lupri,*) 'Symmback: norm of d_ijka', xnorm
         xnorm = ddot(NCKIJ(ISYRES),WORK(KOCC2),1,WORK(KOCC2),1)
         write(lupri,*) 'Symmback: norm of d_ijak', xnorm
         xnorm = ddot(NCKIJ(ISYRES),WORK(KOCCS1),1,WORK(KOCCS1),1)
         write(lupri,*) 'Symmback: norm of ds_ijka', xnorm
      end if

      ISYCMO = 1
      CALL BTRIJKDEL(WORK(KCMO),ISYCMO,WORK(KOCCS1),ISYRES,
     &               LUIJKDS,FNDSIJKD,WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C STEP A.4
C Construct symmetrized vir.vir.occ.vir density and backtransform
C d^s_abic (abic) = d3_abic(aibc) + d4_abci(biac)
C-----------------------------------------------------------------
C
      LENGTH = 0
      DO ISYM = 1, NSYM
        LENGTH = MAX(LENGTH,NCKATR(ISYM))
      END DO
!
      KAIB  = KSTART                         !til den3
      KBIA  = KAIB  + LENGTH                 !til den4
      KABIS = KBIA  + LENGTH 
      KEND1 = KABIS + LENGTH
      LWRK1 = LWORK - KEND1
!
      CALL DZERO(WORK(KAIB),LENGTH)
      CALL DZERO(WORK(KBIA),LENGTH)
      CALL DZERO(WORK(KABIS),LENGTH)

      xnorm = zero
      xnabi1 = zero
      xnabi3 = zero

      DO ISYMC = 1, NSYM
         ISYAIB = MULD2H(ISYRES,ISYMC)
         ISYBIA = ISYAIB
         ISYABI = ISYAIB

         CALL DZERO(WORK(KAIB),NCKATR(ISYAIB))
         CALL DZERO(WORK(KBIA),NCKATR(ISYAIB))
         CALL DZERO(WORK(KABIS),NCKASR(ISYABI))

         DO C = 1, NVIR(ISYMC)

            !determine the offset on file to the given block aib;c
            !and bia;c

            KAIB_C  = ICKBD(ISYAIB,ISYMC)
     &                  + NCKATR(ISYAIB)*(C-1) + 1

            KBIA_C  = ICKBD(ISYBIA,ISYMC)
     &                  + NCKATR(ISYBIA)*(C-1) + 1 

C Read in aib and bia distributions for a given c. Put in Work restarting
C from the beginning for each C (probabilmente poi per quella dopo devo 
C rileggere dal file)
 
            CALL GETWA2(LUABI1,FNDABI1,WORK(KAIB),KAIB_C,
     &                        NCKATR(ISYAIB))
!OBS: scale -1 
            CALL DSCAL(NCKATR(ISYAIB),-ONE,WORK(KAIB),1)

            xnabi1 = xnabi1 + ddot(nckatr(isyaib),work(kaib),
     &                             1,work(kaib),1)
C
            CALL GETWA2(LUABI3,FNDABI3,WORK(KBIA),KBIA_C,
     &                        NCKATR(ISYBIA))

            xnabi3 = xnabi3 + ddot(nckatr(isyaib),work(kbia),
     &                             1,work(kbia),1)

            DO ISYMB = 1, NSYM
               ISYMAI  = MULD2H(ISYAIB,ISYMB)
               DO B = 1, NVIR(ISYMB)
                 DO ISYMI = 1, NSYM
                    ISYMA  = MULD2H(ISYMAI,ISYMI)
                    ISYMBI = MULD2H(ISYMB,ISYMI)
                    ISYMAB = MULD2H(ISYMA,ISYMB)
                    DO I = 1, NRHF(ISYMI)
                    DO A = 1, NVIR(ISYMA)

                       KOFF1 = KAIB - 1 
     &                       + ICKATR(ISYMAI,ISYMB)
     &                       + NT1AM(ISYMAI)*(B-1)
     &                       + IT1AM(ISYMA,ISYMI)
     &                       + NVIR(ISYMA)*(I-1) + A

                       KOFF2 = KBIA - 1
     &                       + ICKATR(ISYMBI,ISYMA)
     &                       + NT1AM(ISYMBI)*(A-1)
     &                       + IT1AM(ISYMB,ISYMI)
     &                       + NVIR(ISYMB)*(I-1) + B

                       KOFF3 = KABIS - 1
     &                       + ICKASR(ISYMAB,ISYMI)
     &                       + NMATAB(ISYMAB)*(I-1)
     &                       + IMATAB(ISYMA,ISYMB)
     &                       + NVIR(ISYMA)*(B-1) + A
C
                       WORK(KOFF3) = WORK(KOFF1) + WORK(KOFF2)

                    END DO !A
                    END DO !I
                 END DO    !ISYMI
               END DO      !B  
            END DO

            KABIS_C = ICDKVI(ISYABI,ISYMC)
     &                  + NCKASR(ISYABI)*(C-1) + 1
            xnorm = xnorm +
     &      ddot(NCKASR(ISYABI),WORK(KABIS),1,WORK(KABIS),1)
 
            CALL PUTWA2(LUABICS,FNDSABIC,WORK(KABIS),KABIS_C,
     &                        NCKASR(ISYABI))
C
         ENDDO                  !C
      ENDDO                     !ISYMC

      if (locdbg) then
        write(lupri,*) '--------------------------'
        write(lupri,*) 'symmback: norm of ds_abic : ', xnorm
        write(lupri,*) 'symmback: norm of d_abic : ', xnabi1
        write(lupri,*) 'symmback: norm of d_abci : ', xnabi3
      end if
C
C-------------------------------
C Backtransform and dump on file
C-------------------------------
C
      CALL BTRABIDEL(LUABIDS,FNDSABID,LUABICS,FNDSABIC,ISYRES,
     &               WORK(KCMO),WORK(KEND1),LWRK1)

C-----------------------------------------------------------------
C Resort
C d3_iabc(bica) = d3_iabc(bica) (nothing to do) ABI1
C d4_aibc(bica) = d4_aibc(ciba)                 ABI3
C-----------------------------------------------------------------
!1) resort d4 ciba->bica and put on file
!2) resort d3 and d4 bica->biac and then biac->aibc
!3) backtransform last index to delta

! iabc

      IOPT = 2
      FACT = -ONE
      CALL RESORT_DAIBC(IOPT,ISYRES,FACT,LUABI1,FNDABI1,
     &                 LUIABCS,FNDSIABC,WORK(KEND1),LWRK1)
C
C-------------------------------
C Backtransform and dump on file
C-------------------------------
C
      CALL BTRAIBDEL(LUIABDE,FNDSIABD,LUIABCS,FNDSIABC,
     &               ISYRES,WORK(KCMO),WORK(KEND1),LWRK1)

!aibc ------------------

      !first resort ciba-->bica
      !
      IOPT = 1
      FACT = ONE
      CALL RESORT_DAIBC(IOPT,ISYRES,FACT,LUABI3,FNDABI3,
     &                 LUABI3,FNDABI3,WORK(KEND1),LWRK1)

      IOPT = 2
      FACT = ONE
      CALL RESORT_DAIBC(IOPT,ISYRES,FACT,LUABI3,FNDABI3,
     &                 LUAIBCS,FNDSAIBC,WORK(KEND1),LWRK1)

C
C-------------------------------
C Backtransform and dump on file
C-------------------------------
C
      CALL BTRAIBDEL(LUAIBDE,FNDSAIBD,LUAIBCS,FNDSAIBC,
     &               ISYRES,WORK(KCMO),WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Close files and delete densities no longer required
C-----------------------------------------------------------------
C
      CALL WCLOSE2(LUIAJDE,FNDSIAJD,'KEEP')
      CALL WCLOSE2(LUAIJDE,FNDSAIJD,'KEEP')
      CALL WCLOSE2(LUIJKDS,FNDSIJKD,'KEEP')
      CALL WCLOSE2(LUABIDS,FNDSABID,'KEEP')
      CALL WCLOSE2(LUAIBDE,FNDSAIBD,'KEEP')
      CALL WCLOSE2(LUIABDE,FNDSIABD,'KEEP')

      CALL WCLOSE2(LUABICS,FNDSABIC,'DELETE')
      CALL WCLOSE2(LUAIBCS,FNDSAIBC,'DELETE')
      CALL WCLOSE2(LUIABCS,FNDSIABC,'DELETE')
      CALL WCLOSE2(LUIJKAS,FNDSIJKA,'DELETE')
C
C-----------------------------------------------------------------
C     Finished, return
C-----------------------------------------------------------------
C
      CALL QEXIT('SYMMBACK')
      RETURN
      END
C------------------------------------------------------------------
C  /* Deck cc_symmback_1 */
      SUBROUTINE SYMMBACK_1(MODEL,
     &                    LUIAJK,FNDIAJK,
     &                    LUAIJK,FNDAIJK,
     &                    LUABI1,FNDABI1,
     &                    LUABI3,FNDABI3,
     &                    LUPTIAJB,FNDIAJB,
     &                    ISYRES,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C     Purpose: construct/resort the symmetrized (T) densities for CCSD(T)
C              resort the other (T) densities
C              and drive backtranformation of last index to AO
C
C     S. Coriani, December 2001/January 2002
C     Special Modified version for Orbit-Orbit corrections. TO be merged back!!!
C-----------------------------------------------------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
Csonia
#include <ccfop.h>

      INTEGER LUPTIAJB, LUAIJK, LUIAJK, LUABI1, LUABI3, LWORK
      INTEGER LUIJKAS, LUIAJKS, LUAIJKS, LUABICS, LUAIBCS, LUIABCS
      INTEGER LUIJKDS, LUIAJDS, LUAIJDS, LUABIDS
      INTEGER LUIAJDE, LUAIJDE, LUAIBDE, LUIABDE
      INTEGER ISYRES, KSTART
      INTEGER KOCC1,KOCC2,KOCCR1,KOCCR2,KOCCS1,LWRK1,KEND1
      INTEGER IOFF, ISYMA, ISYMI, ISYMJ, ISYMK
      INTEGER ISYKJI,ISYMKJ
      integer iii
      integer ISYMIJ, ISYIJK, KIJKA, KIJAK, KIJKAS
      integer LENGTH, KAIB, KBIA, ISYM, KABIS
      integer KDEN4, KDEN3
      integer ISYBIA, ISYAIB, KAIB_C, KBIA_C, KABIS_C, ISYMB
      integer ISYBIC, ISYCIB
      integer KBICS_A 
      integer KOFF2, KOFF3, KOFF1, ISYMAI, ISYMC, ISYMBI
      integer ISYMCI
      integer IOPT,ioptx
      integer kaijk, kiajk, isyaij, isyden, isymab, isyabi
      integer kcmo,kend0,lwrk0,koff
      integer knoddy1,knoddy2,kaibd,kbida
      integer kdaijb, kdaibj, kscrpk, isymbj, nai, naij, nbj
      integer kaibj, kaijb, kiajb,kdiajb
      integer kd1kjia, kd2jkia, kdiajk, kdaijk,isycmo
      integer kbicr4,ioffbic_a,ioffcib_a

      integer LUABICA, LUABIDA,LUIJKDA, KOCCA1, KABIA

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL ZERO, ONE, TWO, FACT
      real ddot,xnorm2,xnorm,xnabi1,xnabi3
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO, FACT
      double precision ddot,xnorm2,xnorm,xnabi1,xnabi3
#endif


      CHARACTER MODEL*10
      CHARACTER MODELPRI*4, MODELPRI2*14
      CHARACTER*(*) FNDIAJK, FNDAIJK, FNDABI1, FNDABI3, FNDIAJB
      LOGICAL LOCDBG
C
C Locally define density file names for easier debug
C
      CHARACTER FNDSIAJK*9, FNDSAIJK*9, FNDSIJKA*9
      CHARACTER FNDSABIC*9, FNDSAIBC*9, FNDSIABC*9
      CHARACTER FNDSIJKD*9, FNDSIAJD*9, FNDSAIJD*9, FNDSABID*9
      CHARACTER FNDSAIBD*9, FNDSIABD*9
c
c     Special extra antisymmetrized densities for BP2EOO
c
      CHARACTER FNDAABIC*9, FNDAIJKD*9, FNDAABID*9
      PARAMETER (FNDAABIC = 'PT_DAABIC')
      PARAMETER (FNDAABID = 'PT_DAABID',FNDAIJKD='PT_DAIJKD')
c
      PARAMETER (FNDSIAJK = 'PT_DSIAJK', FNDSAIJK = 'PT_DSAIJK')
      PARAMETER (FNDSIJKA = 'PT_DSIJKA', FNDSABIC = 'PT_DSABIC')
      PARAMETER (FNDSAIBC = 'PT_DSAIBC', FNDSIABC = 'PT_DSIABC')
c
      PARAMETER (FNDSIJKD = 'PT_DSIJKD', FNDSIAJD = 'PT_DSIAJD')
      PARAMETER (FNDSAIJD = 'PT_DSAIJD', FNDSABID = 'PT_DSABID')
      PARAMETER (FNDSAIBD = 'PT_DSAIBD', FNDSIABD = 'PT_DSIABD')
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      PARAMETER (LOCDBG = .FALSE.)
C
C#include <leinf.h>
C
      CALL QENTER('SYMMBACK_1')
C
C-------------------------------
C     First allocation
C-------------------------------
C
      KCMO   = 1
      KEND0 = KCMO   + NLAMDS
      LWRK0 = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND0
         CALL QUIT(
     *         'Insufficient memory for initial allocation')
      ENDIF
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE (LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND0),LWRK0)
C
C---------------------------------------------------------
C     Open files for backtransformed/symmetrized densities
C---------------------------------------------------------
C
!      LUIJKAS  = -1
      LUABICS  = -1
      LUIJKDS  = -1
      LUABIDS  = -1
C
C Temporarily keep the other "S" files with the resorted densities
C (instead of symmetrized)
C 
!      LUIAJKS = -1
!      LUAIJKS = -1
      LUAIBCS = -1
      LUIABCS = -1
C
      LUIAJDE  = -1
      LUAIJDE  = -1
      LUAIBDE  = -1
      LUIABDE  = -1
C
C Symmetrized and backtransformed
C
C     d^s_{ijka}
!      CALL WOPEN2(LUIJKAS,FNDSIJKA,64,0)
C     d^s_{ijkdelta}
      CALL WOPEN2(LUIJKDS,FNDSIJKD,64,0)
C     d^s_{abic}
      CALL WOPEN2(LUABICS,FNDSABIC,64,0)
C     d^s_{abidelta}
      CALL WOPEN2(LUABIDS,FNDSABID,64,0)
C
C Only backtransformed
C
C     d^s_{iajdelta}
      CALL WOPEN2(LUIAJDE,FNDSIAJD,64,0)
C     d^s_{aijdelta}
      CALL WOPEN2(LUAIJDE,FNDSAIJD,64,0)
C     d^s_{aibdelta}
      CALL WOPEN2(LUAIBDE,FNDSAIBD,64,0)
C     d^s_{iabdelta}
      CALL WOPEN2(LUIABDE,FNDSIABD,64,0)
C
C Only symmetrized
C
C     d^s_{aibc}
      CALL WOPEN2(LUAIBCS,FNDSAIBC,64,0)
C     d^s_{iabc}
      CALL WOPEN2(LUIABCS,FNDSIABC,64,0)

      if (BP2EOO) then

!        LUIJKAA  = -1
        LUABICA  = -1
        LUIJKDA  = -1
        LUABIDA  = -1
C
C Antisymmetrized and backtransformed
C
C       d^a_{ijka}
!        CALL WOPEN2(LUIJKAA,FNDAIJKA,64,0)
C       d^a_{ijkdelta}
        CALL WOPEN2(LUIJKDA,FNDAIJKD,64,0)
C       d^a_{abic}
        CALL WOPEN2(LUABICA,FNDAABIC,64,0)
C       d^a_{abidelta}
        CALL WOPEN2(LUABIDA,FNDAABID,64,0)

      end if
C
C---------------------------------------------------------------------
C - Read in iajk, aijk, iajb (the latter read-in inside backtransform)
C - Resort iajk, aijk
C - Backtransform last index and save on file
C---------------------------------------------------------------------

      KSTART = KEND0
      KOCC1  = KSTART
      KOCC2  = KOCC1 + NCKIJ(ISYRES)
      KOCCR1 = KOCC2 + NCKIJ(ISYRES)  !resorted density 1 diajk
      KOCCR2 = KOCCR1 + NCKIJ(ISYRES) !resorted density 2 daijk
      KEND1  = KOCCR2 + NCKIJ(ISYRES)
      LWRK1  = LWORK - KEND1
C
C-----------
C Initialize
C-----------
C
      CALL DZERO(WORK(KOCC1),NCKIJ(ISYRES))
      CALL DZERO(WORK(KOCC2),NCKIJ(ISYRES))
      CALL DZERO(WORK(KOCCR1),NCKIJ(ISYRES))
      CALL DZERO(WORK(KOCCR2),NCKIJ(ISYRES))
C
C--------------------------------------------------------------
C Read in the entire diajk(kjia) and daijk(jkia) distributions, 
C and keep in WORK(KOCC1) and WORK(KOCC2) for later use
C--------------------------------------------------------------
C
      IOFF = 1
      CALL GETWA2(LUIAJK,FNDIAJK,WORK(KOCC1),IOFF,NCKIJ(ISYRES))
C
!OBS: scale -1 
C
      CALL DSCAL(NCKIJ(ISYRES),-ONE,WORK(KOCC1),1)

      CALL GETWA2(LUAIJK,FNDAIJK,WORK(KOCC2),IOFF,NCKIJ(ISYRES))
C
!OBS: scale -1 
C
      CALL DSCAL(NCKIJ(ISYRES),-ONE,WORK(KOCC2),1)

C
C------------------------------------
C Resort diajk(kjia) to diajk(ai,j,k)
C Resort daijk(jkia) to daijk(ai,j,k)
C------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYKJI = MULD2H(ISYRES,ISYMA)
         DO A = 1, NVIR(ISYMA)
            DO ISYMI = 1, NSYM
               ISYMKJ = MULD2H(ISYKJI,ISYMI)
               DO I = 1, NRHF(ISYMI)
                  DO ISYMJ = 1, NSYM
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYAIJ = MULD2H(ISYMAI,ISYMJ)
                     ISYMK = MULD2H(ISYMKJ,ISYMJ)
                     DO J = 1, NRHF(ISYMJ)
                     DO K = 1, NRHF(ISYMK)
                        KD1KJIA = I3OVIR(ISYKJI,ISYMA) 
     &                          + NMAIJK(ISYKJI)*(A-1) 
     &                          + IMAIJK(ISYMKJ,ISYMI)
     &                          + NMATIJ(ISYMKJ)*(I-1)
     &                          + IMATIJ(ISYMK,ISYMJ)
     &                          + NRHF(ISYMK)*(J-1)
     &                          + K 
     &                          + KOCC1 - 1

                        KD2JKIA = I3OVIR(ISYKJI,ISYMA) 
     &                          + NMAIJK(ISYKJI)*(A-1) 
     &                          + IMAIJK(ISYMKJ,ISYMI)
     &                          + NMATIJ(ISYMKJ)*(I-1)
     &                          + IMATIJ(ISYMJ,ISYMK)
     &                          + NRHF(ISYMJ)*(K-1)
     &                          + J 
     &                          + KOCC2 - 1

                        KDAIJK = ICKITR(ISYAIJ,ISYMK)
     &                          + NCKI(ISYAIJ)*(K-1)
     &                          + ICKI(ISYMAI,ISYMJ)
     &                          + NT1AM(ISYMAI)*(J-1)
     &                          + IT1AM(ISYMA,ISYMI)
     &                          + NVIR(ISYMA)*(I-1) + A

                        KDIAJK = KDAIJK               !d_iajk is stored aijk

                        WORK(KDIAJK + KOCCR1 - 1) = WORK(KD1KJIA)
                        WORK(KDAIJK + KOCCR2 - 1) = WORK(KD2JKIA)

                     END DO         !K
                     END DO         !J
                  END DO            !ISYMJ
               END DO               !I
            END DO                  !ISYMI
         END DO                     !A
      END DO                        !ISYMA

C
C Save on file
C
!      IOFF = 1
!      CALL PUTWA2(LUIAJKS,FNDSIAJK,WORK(KOCCR1),IOFF,NCKIJ(ISYRES))
!      CALL PUTWA2(LUAIJKS,FNDSAIJK,WORK(KOCCR2),IOFF,NCKIJ(ISYRES))
C
C--------------------------------------------------------
C     Read-in d_iajb(ai,bj), unpack and reorder to ai,j,b     
C--------------------------------------------------------

      KDIAJB = KEND1
      KDAIBJ = KDIAJB + NT2SQ(ISYRES)
      KSCRPK = KDAIBJ + NT2SQ(ISYRES)
      KEND1  = KSCRPK + NT2AM(ISYRES)
      LWRK1  = LWORK - KEND1

      CALL DZERO(WORK(KDIAJB),NT2SQ(ISYRES))
      CALL DZERO(WORK(KDAIBJ),NT2SQ(ISYRES))

      IF (NT2AM(ISYRES) .GT. 0) THEN
         IOFF = 1
         CALL GETWA2(LUPTIAJB,FNDIAJB,WORK(KSCRPK),
     &               IOFF,NT2AM(ISYRES))
         !square up
         CALL CC_T2SQ(WORK(KSCRPK),WORK(KDAIBJ),1)
         !Husk: scale by 2*
         CALL DSCAL(NT2SQ(ISYRES),TWO,WORK(KDAIBJ),1)
      ENDIF
 
      CALL CC3_T2TP(WORK(KDIAJB),WORK(KDAIBJ),ISYRES)
C
C-----------------------------------------------
C Backtransform
C-----------------------------------------------
C
      IOPTX=2
      CALL BTRIAJDEL(IOPTX,WORK(KCMO),WORK(KOCCR1),WORK(KDIAJB),
     &               ISYRES,LUIAJDE,FNDSIAJD,WORK(KEND1),LWRK1)
!    --------------------------------------------------------
!      write(lupri,*)'Symmback: call BTRAIJDEL (aijdel) ' 
!      CALL BTRAIJDEL(WORK(KCMO),WORK(KOCCR2),ISYRES,LUAIJDE,
!     &               FNDSAIJD,WORK(KEND1),LWRK1)
!      write(lupri,*) 'Symmback: exit BTRAIJDEL (aijdelta)'
!      call flshfo(lupri)
!    --------------------------------------------------------
      IOPTX=1
      CALL BTRIAJDEL(IOPTX,WORK(KCMO),WORK(KOCCR2),WORK(KDIAJB),
     &           ISYRES,LUAIJDE,FNDSAIJD,WORK(KEND1),LWRK1)
C
C-------------------------------------------------------------------
C STEP A.2
C
C Construct symmetrized occ.occ.occ.vir and backtransform last index
C d^s_ijka (ijka) = d2_ijak(ijka) + d1_ijka(jika)
C
C If (BP2EOO) antisymmetrize as well
C
C-------------------------------------------------------------------
C
      KOCCS1 = KOCCR1
      KEND1  = KOCCS1 + NCKIJ(ISYRES)
      LWRK1  = LWORK - KEND1
      
      CALL DZERO(WORK(KOCCS1),NCKIJ(ISYRES))

      if (BP2EOO) then

        KOCCA1 = KEND1
        KEND1  = KOCCA1 + NCKIJ(ISYRES)
        LWRK1  = LWORK - KEND1

      end if

      DO ISYMA = 1, NSYM
         ISYIJK = MULD2H(ISYRES,ISYMA)
         DO A = 1, NVIR(ISYMA)
            DO ISYMK = 1, NSYM
               ISYMIJ = MULD2H(ISYIJK,ISYMK)
               DO K = 1, NRHF(ISYMK)
                  DO ISYMJ = 1, NSYM
                     ISYMI = MULD2H(ISYMIJ,ISYMJ)
                     DO J = 1, NRHF(ISYMJ)
                     DO I = 1, NRHF(ISYMI)

                        KIJAK = I3OVIR(ISYIJK,ISYMA) 
     &                        + NMAIJK(ISYIJK)*(A-1) 
     &                        + IMAIJK(ISYMIJ,ISYMK)
     &                        + NMATIJ(ISYMIJ)*(K-1) 
     &                        + IMATIJ(ISYMI,ISYMJ)
     &                        + NRHF(ISYMI)*(J-1) 
     &                        + I
     &                        + KOCC2 - 1


                        KIJKA = I3OVIR(ISYIJK,ISYMA) 
     &                        + NMAIJK(ISYIJK)*(A-1) 
     &                        + IMAIJK(ISYMIJ,ISYMK)
     &                        + NMATIJ(ISYMIJ)*(K-1) 
     &                        + IMATIJ(ISYMJ,ISYMI)
     &                        + NRHF(ISYMJ)*(I-1) 
     &                        + J
     &                        + KOCC1 - 1

                        KIJKAS = I3OVIR(ISYIJK,ISYMA)
     &                         + NMAIJK(ISYIJK)*(A-1) 
     &                         + IMAIJK(ISYMIJ,ISYMK)
     &                         + NMATIJ(ISYMIJ)*(K-1) 
     &                         + IMATIJ(ISYMI,ISYMJ)
     &                         + NRHF(ISYMI)*(J-1) 
     &                         + I
     &                         + KOCCS1 - 1


                        WORK(KIJKAS) = WORK(KIJAK) + WORK(KIJKA)

                        if (BP2EOO) then
                           WORK(KIJKAS - KOCCS1 + KOCCA1) = 
     &                               - WORK(KIJAK) + WORK(KIJKA)
                        end if

                     END DO         !I
                     END DO         !J
                  END DO            !ISYMJ
               END DO               !K
            END DO                  !ISYMK
         END DO                     !A
      END DO                        !ISYMA

C
C
!      IOFF = 1
!      CALL PUTWA2(LUIJKAS,FNDSIJKA,WORK(KOCCS1),IOFF,NCKIJ(ISYRES))
!      CALL WCLOSE2(LUIJKAS,FNDSIJKA,'KEEP')
C
C Backtransform

      if (locdbg) then
         xnorm = ddot(NCKIJ(ISYRES),WORK(KOCC1),1,WORK(KOCC1),1)
         write(lupri,*) 'Symmback: norm of d_ijka', xnorm
         xnorm = ddot(NCKIJ(ISYRES),WORK(KOCC2),1,WORK(KOCC2),1)
         write(lupri,*) 'Symmback: norm of d_ijak', xnorm
         xnorm = ddot(NCKIJ(ISYRES),WORK(KOCCS1),1,WORK(KOCCS1),1)
         write(lupri,*) 'Symmback: norm of ds_ijka', xnorm
      end if

      ISYCMO = 1
      CALL BTRIJKDEL(WORK(KCMO),ISYCMO,WORK(KOCCS1),ISYRES,
     &               LUIJKDS,FNDSIJKD,WORK(KEND1),LWRK1)

      if (BP2EOO) then
         ISYCMO = 1
         CALL BTRIJKDEL(WORK(KCMO),ISYCMO,WORK(KOCCA1),ISYRES,
     &               LUIJKDA,FNDAIJKD,WORK(KEND1),LWRK1)
      end if
C
C-----------------------------------------------------------------
C STEP A.4
C Construct symmetrized vir.vir.occ.vir density and backtransform
C d^s_abic (abic) = d3_abic(aibc) + d4_abci(biac)
C-----------------------------------------------------------------
C
      LENGTH = 0
      DO ISYM = 1, NSYM
        LENGTH = MAX(LENGTH,NCKATR(ISYM))
      END DO
!
      KAIB  = KSTART                         !til den3
      KBIA  = KAIB  + LENGTH                 !til den4
      KABIS = KBIA  + LENGTH 
      KEND1 = KABIS + LENGTH
      LWRK1 = LWORK - KEND1
!
      CALL DZERO(WORK(KAIB),LENGTH)
      CALL DZERO(WORK(KBIA),LENGTH)
      CALL DZERO(WORK(KABIS),LENGTH)


      if (BP2EOO) then

         KABIA = KABIS + LENGTH 
         KEND1 = KABIA + LENGTH
         LWRK1 = LWORK - KEND1
         CALL DZERO(WORK(KABIA),LENGTH)

      end if

      xnorm = zero
      xnabi1 = zero
      xnabi3 = zero

      DO ISYMC = 1, NSYM
         ISYAIB = MULD2H(ISYRES,ISYMC)
         ISYBIA = ISYAIB
         ISYABI = ISYAIB

         CALL DZERO(WORK(KAIB),NCKATR(ISYAIB))
         CALL DZERO(WORK(KBIA),NCKATR(ISYAIB))
         CALL DZERO(WORK(KABIS),NCKASR(ISYABI))

         DO C = 1, NVIR(ISYMC)

            !determine the offset on file to the given block aib;c
            !and bia;c

            KAIB_C  = ICKBD(ISYAIB,ISYMC)
     &                  + NCKATR(ISYAIB)*(C-1) + 1

            KBIA_C  = ICKBD(ISYBIA,ISYMC)
     &                  + NCKATR(ISYBIA)*(C-1) + 1 

C Read in aib and bia distributions for a given c. Put in Work restarting
C from the beginning for each C (probabilmente poi per quella dopo devo 
C rileggere dal file)
 
            CALL GETWA2(LUABI1,FNDABI1,WORK(KAIB),KAIB_C,
     &                        NCKATR(ISYAIB))
!OBS: scale -1 
            CALL DSCAL(NCKATR(ISYAIB),-ONE,WORK(KAIB),1)

            xnabi1 = xnabi1 + ddot(nckatr(isyaib),work(kaib),
     &                             1,work(kaib),1)
C
            CALL GETWA2(LUABI3,FNDABI3,WORK(KBIA),KBIA_C,
     &                        NCKATR(ISYBIA))

            xnabi3 = xnabi3 + ddot(nckatr(isyaib),work(kbia),
     &                             1,work(kbia),1)

            DO ISYMB = 1, NSYM
               ISYMAI  = MULD2H(ISYAIB,ISYMB)
               DO B = 1, NVIR(ISYMB)
                 DO ISYMI = 1, NSYM
                    ISYMA  = MULD2H(ISYMAI,ISYMI)
                    ISYMBI = MULD2H(ISYMB,ISYMI)
                    ISYMAB = MULD2H(ISYMA,ISYMB)
                    DO I = 1, NRHF(ISYMI)
                    DO A = 1, NVIR(ISYMA)

                       KOFF1 = KAIB - 1 
     &                       + ICKATR(ISYMAI,ISYMB)
     &                       + NT1AM(ISYMAI)*(B-1)
     &                       + IT1AM(ISYMA,ISYMI)
     &                       + NVIR(ISYMA)*(I-1) + A

                       KOFF2 = KBIA - 1
     &                       + ICKATR(ISYMBI,ISYMA)
     &                       + NT1AM(ISYMBI)*(A-1)
     &                       + IT1AM(ISYMB,ISYMI)
     &                       + NVIR(ISYMB)*(I-1) + B

                       KOFF3 = KABIS - 1
     &                       + ICKASR(ISYMAB,ISYMI)
     &                       + NMATAB(ISYMAB)*(I-1)
     &                       + IMATAB(ISYMA,ISYMB)
     &                       + NVIR(ISYMA)*(B-1) + A
C
                       WORK(KOFF3) = WORK(KOFF1) + WORK(KOFF2)

                       if (BP2EOO) then

                          WORK(KOFF3 - KABIS + KABIA) = 
     &                        WORK(KOFF1) - WORK(KOFF2)

                       end if

                    END DO !A
                    END DO !I
                 END DO    !ISYMI
               END DO      !B  
            END DO

            KABIS_C = ICDKVI(ISYABI,ISYMC)
     &                  + NCKASR(ISYABI)*(C-1) + 1
            xnorm = xnorm +
     &      ddot(NCKASR(ISYABI),WORK(KABIS),1,WORK(KABIS),1)
 
            CALL PUTWA2(LUABICS,FNDSABIC,WORK(KABIS),KABIS_C,
     &                        NCKASR(ISYABI))

            if (BP2EOO) then
               CALL PUTWA2(LUABICA,FNDAABIC,WORK(KABIA),KABIS_C,
     &                        NCKASR(ISYABI))
            end if

C
         ENDDO                  !C
         !CALL DZERO(WORK(KABIS),NCKASR(ISYABI))
      ENDDO                     !ISYMC
      if (locdbg) then
        write(lupri,*) '--------------------------'
        write(lupri,*) 'symmback: norm of ds_abic : ', xnorm
        write(lupri,*) 'symmback: norm of d_abic : ', xnabi1
        write(lupri,*) 'symmback: norm of d_abci : ', xnabi3
      end if
C
C-------------------------------
C Backtransform and dump on file
C-------------------------------
C
      CALL BTRABIDEL(LUABIDS,FNDSABID,LUABICS,FNDSABIC,ISYRES,
     &               WORK(KCMO),WORK(KEND1),LWRK1)
      if (BP2EOO) then
         CALL BTRABIDEL(LUABIDA,FNDAABID,LUABICA,FNDAABIC,ISYRES,
     &               WORK(KCMO),WORK(KEND1),LWRK1)
      end if

C-----------------------------------------------------------------
C Resort
C d3_iabc(bica) = d3_iabc(bica) (nothing to do) ABI1
C d4_aibc(bica) = d4_aibc(ciba)                 ABI3
C-----------------------------------------------------------------
!1) resort d4 ciba->bica and put on file
!2) resort d3 and d4 bica->biac and then biac->aibc
!3) backtransform last index to delta

! iabc

      IOPT = 2
      FACT = -ONE
      CALL RESORT_DAIBC(IOPT,ISYRES,FACT,LUABI1,FNDABI1,
     &                 LUIABCS,FNDSIABC,WORK(KEND1),LWRK1)
C
C-------------------------------
C Backtransform and dump on file
C-------------------------------
C
      CALL BTRAIBDEL(LUIABDE,FNDSIABD,LUIABCS,FNDSIABC,
     &               ISYRES,WORK(KCMO),WORK(KEND1),LWRK1)

!aibc ------------------

      !first resort ciba-->bica
      !
      IOPT = 1
      FACT = ONE
      CALL RESORT_DAIBC(IOPT,ISYRES,FACT,LUABI3,FNDABI3,
     &                 LUABI3,FNDABI3,WORK(KEND1),LWRK1)

      IOPT = 2
      FACT = ONE
      CALL RESORT_DAIBC(IOPT,ISYRES,FACT,LUABI3,FNDABI3,
     &                 LUAIBCS,FNDSAIBC,WORK(KEND1),LWRK1)

C
C-------------------------------
C Backtransform and dump on file
C-------------------------------
C
      CALL BTRAIBDEL(LUAIBDE,FNDSAIBD,LUAIBCS,FNDSAIBC,
     &               ISYRES,WORK(KCMO),WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Close files and delete densities no longer required
C-----------------------------------------------------------------
C
      CALL WCLOSE2(LUIAJDE,FNDSIAJD,'KEEP')
      CALL WCLOSE2(LUAIJDE,FNDSAIJD,'KEEP')
      CALL WCLOSE2(LUIJKDS,FNDSIJKD,'KEEP')
      CALL WCLOSE2(LUABIDS,FNDSABID,'KEEP')
      CALL WCLOSE2(LUAIBDE,FNDSAIBD,'KEEP')
      CALL WCLOSE2(LUIABDE,FNDSIABD,'KEEP')

      IF (BP2EOO) THEN
        CALL WCLOSE2(LUIJKDA,FNDAIJKD,'KEEP')
        CALL WCLOSE2(LUABIDA,FNDAABID,'KEEP')
      END IF

      CALL WCLOSE2(LUABICS,FNDSABIC,'DELETE')
      CALL WCLOSE2(LUAIBCS,FNDSAIBC,'DELETE')
      CALL WCLOSE2(LUIABCS,FNDSIABC,'DELETE')
!      CALL WCLOSE2(LUIJKAS,FNDSIJKA,'DELETE')
C
C-----------------------------------------------------------------
C     Finished, return
C-----------------------------------------------------------------
C
      CALL QEXIT('SYMMBACK_1')
      RETURN
      END
