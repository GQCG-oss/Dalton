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
*=====================================================================*
      SUBROUTINE CCSDT_ADEN_NODDY(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,DIA,
     *                            WORK,LWORK,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
*---------------------------------------------------------------------*
*     Purpose: compute density corresponding to the triples           *
*              contributions to a L A{O} R contraction                *
*                                                                     *
*     Written by Christof Haettig, Jan 2003, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL LOCDBG, HAVETB3T3, PRINTDENS
      PARAMETER (LOCDBG=.FALSE., HAVETB3T3=.FALSE., PRINTDENS=.FALSE.)

      CHARACTER*3 LISTL, LISTR
      CHARACTER*(*) FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3
      CHARACTER*(*) FN3FOPX, FN3FOP2X
      INTEGER LUDELD,LUCKJD,LUDKBC,LUTOC,LU3VI,LUDKBC3,LU3FOPX,LU3FOP2X
      INTEGER LWORK, IDLSTL, IDLSTR

#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      REAL DIJ(*), DAB(*), DIA(*)
      REAL WORK(LWORK)
      REAL HALF, ONE, ZERO, TWO, DDOT, XNORM
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      DOUBLE PRECISION DIJ(*), DAB(*), DIA(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION HALF, ONE, TWO, ZERO, DDOT, XNORM
#endif
      PARAMETER( ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TWO = 2.0D0 )

      CHARACTER MODEL*10
      INTEGER KEND1, LWRK1, KWBAR, KWYINT, KTHETA, KTBAR3, KTAMP3,
     &        NCK, NBJ, NAI, KAIBJCK, KCKAIBJ, KBJCKAI, IOFF, ISYMR,
     &        KWTILDE, KWBREVE, KTHEOCC, KTVIRT1, KTVIRT2, KEND2, 
     &        LWRK2, KTHVIRT, KSWAP, KDABR, KDIJR, KDIAR, IOPT,
     &        KT0AMP, KTYAMP, KT0FMI, KTYFMI, KT0MI, KTYMI, 
     &        KWDE, KWED, KSCR, KFNEMDL, KEMDLFN, K3AMP0,
     &        NAN, NFN, NFI, NEM, NDL, KDLEMFN, KDLAN, KEMFI,
     &        KL2AMP, KDLFIAN, KDLFNAI, KDLFN, ISYML,
     &        KDAB0, KDAB0R, KDIJ0, KDIJ0R, KT1AMY

      INTEGER ILSTSYM

*---------------------------------------------------------------------*
*     Begin:
*---------------------------------------------------------------------*
      KEND1 = 1

      KDAB0 = KEND1
      KDIJ0 = KDAB0  + NVIRT*NVIRT
      KEND1 = KDIJ0  + NRHFT*NRHFT

      IF (HAVETB3T3) THEN
        KTBAR3 = KEND1
        KTAMP3 = KTBAR3 + NT1AMX*NT1AMX*NT1AMX
        KDABR  = KTAMP3 + NT1AMX*NT1AMX*NT1AMX
        KDIJR  = KDABR  + NVIRT*NVIRT
        KDIAR  = KDIJR  + NRHFT*NRHFT
        KDAB0R = KDIAR  + NVIRT*NRHFT
        KDIJ0R = KDAB0R + NVIRT*NVIRT
        KEND1  = KDIJ0R + NRHFT*NRHFT
      END IF

      IF (NSYM.GT.1) CALL QUIT('No symmetry in CCSDT_ADEN_NODDY!')

*---------------------------------------------------------------------*
*     Get WBAR intermediate into memory:
*---------------------------------------------------------------------*
      KWBAR = KEND1
      KEND1 = KWBAR + NT1AMX*NT1AMX*NT1AMX
      LWRK1 = LWORK - KEND1

      IF (LWRK1.LT.0)CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (1)')

      CALL CCSDT_GWBIC(LISTL,IDLSTL,WORK(KWBAR),
     &                 XLAMDP0,XLAMDH0,FOCK0,
     &                 LUTOC,FNTOC,LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     &                 LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     &                 WORK(KEND1),LWRK1)

      IF (LOCDBG) THEN
        XNORM = DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KWBAR),1,WORK(KWBAR),1)
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> Wbar^Y:',XNORM
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTL,IDLSTL:',LISTL,IDLSTL
        CALL PRINT_PT3_NODDY(WORK(KWBAR))

        IF (HAVETB3T3) THEN
         ! accumulate Tbar^Y amplitudes in WORK and print
         DO NCK = 1, NT1AMX
           DO NBJ = 1, NT1AMX
             DO NAI = 1, NT1AMX
               KAIBJCK = ((NCK-1)*NT1AMX+NBJ-1)*NT1AMX+NAI
               KCKAIBJ = ((NBJ-1)*NT1AMX+NAI-1)*NT1AMX+NCK
               KBJCKAI = ((NAI-1)*NT1AMX+NCK-1)*NT1AMX+NBJ
               WORK(KTBAR3-1+KAIBJCK) =  WORK(KWBAR-1+KAIBJCK) +
     &          WORK(KWBAR-1+KCKAIBJ) + WORK(KWBAR-1+KBJCKAI) 
             END DO
           END DO
         END DO
         WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> Tbar^Y:'
         WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTL,IDLSTL:',LISTL,IDLSTL
         CALL PRINT_PT3_NODDY(WORK(KTBAR3))
        END IF
      END IF

    1 FORMAT(1X,A11,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1,I3,A4,4E20.10)
*---------------------------------------------------------------------*
*     Get w^y and theta intermediate into memory:
*---------------------------------------------------------------------*
      KWYINT = KEND1
      KTHETA = KWYINT + NT1AMX*NT1AMX*NT1AMX
      K3AMP0 = KTHETA + NT1AMX*NT1AMX*NT1AMX
      KEND1  = K3AMP0 + NT1AMX*NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1

      IF (LWRK1.LT.0)CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (2)')

      CALL CCSDT_GWTIC(LISTR,IDLSTR,WORK(KWYINT),WORK(KTHETA),
     &                 WORK(K3AMP0 ),XLAMDP0,XLAMDH0,FOCK0,
     &                 LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     &                 WORK(KEND1),LWRK1)

      ! turn sign of t3, w and theta intermediates
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(K3AMP0),1)
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KWYINT),1)
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KTHETA),1)

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> W^Y:'
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTR,IDLSTR:',LISTR,IDLSTR
        CALL PRINT_PT3_NODDY(WORK(KWYINT))
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> Theta^Y:'
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTR,IDLSTR:',LISTR,IDLSTR
        CALL PRINT_PT3_NODDY(WORK(KTHETA))

        IF (HAVETB3T3) THEN
         ! accumulate T^Y amplitudes in WORK and print
         DO NCK = 1, NT1AMX
           DO NBJ = 1, NT1AMX
             DO NAI = 1, NT1AMX
               KAIBJCK = ((NCK-1)*NT1AMX+NBJ-1)*NT1AMX+NAI
               KCKAIBJ = ((NBJ-1)*NT1AMX+NAI-1)*NT1AMX+NCK
               KBJCKAI = ((NAI-1)*NT1AMX+NCK-1)*NT1AMX+NBJ
               WORK(KTAMP3-1+KAIBJCK) =  
     &            WORK(KWYINT-1+KAIBJCK) + WORK(KTHETA-1+KAIBJCK) +
     &            WORK(KWYINT-1+KCKAIBJ) + WORK(KTHETA-1+KCKAIBJ) + 
     &            WORK(KWYINT-1+KBJCKAI) + WORK(KTHETA-1+KBJCKAI) 
             END DO
           END DO
         END DO
         WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> T^Y:'
         WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTR,IDLSTR:',LISTR,IDLSTR
         CALL PRINT_PT3_NODDY(WORK(KTAMP3))
        END IF
      END IF

*---------------------------------------------------------------------*
*     read doubles amplitudes (zeroth-order and response) from file:
*---------------------------------------------------------------------*
      KL2AMP  = KEND1 
      KEND1   = KL2AMP + NT1AMX*NT1AMX

      KT0AMP  = KEND1
      KTYAMP  = KT0AMP + NT1AMX*NT1AMX
      KEND2   = KTYAMP + NT1AMX*NT1AMX
      LWRK2   = LWORK  - KEND2
      
      IF (LWRK2.LT.NT2AMX)
     &   CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (3a)')

      ! read T^0 doubles amplitudes from file and square up
      IOPT   = 2
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND2))
      CALL CC_T2SQ(WORK(KEND2),WORK(KTYAMP),ISYM0)

      ISYMR = ILSTSYM(LISTR,IDLSTR)

      ! read T^X doubles amplitudes from file and square up
      IOPT = 2
      Call CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,DUMMY,WORK(KEND2))
      Call CCLR_DIASCL(WORK(KEND2),TWO,ISYMR)  
      CALL CC_T2SQ(WORK(KEND2),WORK(KT0AMP),ISYMR)

      ISYML = ILSTSYM(LISTL,IDLSTL)

      ! read left doubles amplitudes from file and square up
      IOPT = 2
      Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,DUMMY,WORK(KEND2))
      CALL CC_T2SQ(WORK(KEND2),WORK(KL2AMP),ISYML)

*---------------------------------------------------------------------*
*     if we have Tbar3 and T3 in core evaluate reference density
*     using the simple formulas:
*         D_ab = +1/2 sum_dlemn tbar(dl,em,an) ty(dl,em,bn)
*         D_ij = -1/2 sum_dlemf ty(dl,em,fi) tbar(dl,em,fj)
*         D_ia = -sum_dln sum_efm tbar(dl,em,fn) ty(em,fi) t0(dl,an)
*                -sum_dln sum_efm tbar(dl,em,fn) t0(em,fi) ty(dl,an)
*                +sum_dlem tbar(dl,em) [ty(dl,em,ai)-ty(dl,ei,am)] 
*---------------------------------------------------------------------*
      IF (HAVETB3T3) THEN
        CALL DZERO(WORK(KDABR),NVIRT*NVIRT)
        CALL DZERO(WORK(KDIJR),NRHFT*NRHFT)
        CALL DZERO(WORK(KDIAR),NVIRT*NRHFT)
        CALL DZERO(WORK(KDAB0R),NVIRT*NVIRT)
        CALL DZERO(WORK(KDIJ0R),NRHFT*NRHFT)

        DO N = 1, NRHFT
          IOFF = NT1AMX*NT1AMX*NVIRT*(N-1)
          CALL DGEMM('T','N',NVIRT,NVIRT,     NT1AMX*NT1AMX,
     &                +HALF,WORK(KTBAR3+IOFF),NT1AMX*NT1AMX,
     &                      WORK(KTAMP3+IOFF),NT1AMX*NT1AMX,
     &                  ONE,WORK(KDABR),NVIRT)

          CALL DGEMM('T','N',NVIRT,NVIRT,     NT1AMX*NT1AMX,
     &                +HALF,WORK(KTBAR3+IOFF),NT1AMX*NT1AMX,
     &                      WORK(K3AMP0+IOFF),NT1AMX*NT1AMX,
     &                  ONE,WORK(KDAB0R),NVIRT)
        END DO

        CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NT1AMX*NVIRT,
     &              -HALF,WORK(KTAMP3),NT1AMX*NT1AMX*NVIRT,
     &                    WORK(KTBAR3),NT1AMX*NT1AMX*NVIRT,
     &               ZERO,WORK(KDIJR),NRHFT)

        CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NT1AMX*NVIRT,
     &              -HALF,WORK(K3AMP0),NT1AMX*NT1AMX*NVIRT,
     &                    WORK(KTBAR3),NT1AMX*NT1AMX*NVIRT,
     &               ZERO,WORK(KDIJ0R),NRHFT)

        DO I = 1, NRHFT
         DO A = 1, NVIRT
          DO N = 1, NRHFT
           DO F = 1, NVIRT

            NAN = NVIRT*(N-1) + A
            NAI = NVIRT*(I-1) + A
            NFN = NVIRT*(N-1) + F
            NFI = NVIRT*(I-1) + F

            DO NDL = 1, NT1AMX
              KDLFIAN = ((NAN-1)*NT1AMX + NFI-1) * NT1AMX + NDL    
              KDLFNAI = ((NAI-1)*NT1AMX + NFN-1) * NT1AMX + NDL    
              KDLFN   =  (NFN-1)*NT1AMX + NDL    

              WORK(KDIAR-1+NAI) = WORK(KDIAR-1+NAI) + 
     &         ( WORK(KTAMP3-1+KDLFNAI) - WORK(KTAMP3-1+KDLFIAN) ) * 
     &           WORK(KL2AMP-1+KDLFN) 
            END DO

            DO NEM = 1, NT1AMX
              DO NDL = 1, NT1AMX
                KDLEMFN = ((NFN-1)*NT1AMX + NEM-1) * NT1AMX + NDL    
                KFNEMDL = ((NDL-1)*NT1AMX + NEM-1) * NT1AMX + NFN    
                KEMDLFN = ((NFN-1)*NT1AMX + NDL-1) * NT1AMX + NEM    
                KDLAN   =  (NAN-1)*NT1AMX + NDL
                KEMFI   =  (NFI-1)*NT1AMX + NEM

                WORK(KDIAR-1+NAI) = WORK(KDIAR-1+NAI) - 
     &            WORK(KTBAR3-1+KDLEMFN) *
     &             (  WORK(KT0AMP-1+KEMFI) * WORK(KTYAMP-1+KDLAN) + 
     &                WORK(KTYAMP-1+KEMFI) * WORK(KT0AMP-1+KDLAN)   )

              END DO
            END DO

           END DO
          END DO
         END DO
        END DO

        IF (LOCDBG) THEN
          WRITE(LUPRI,*) 'The virtual-virtual block of the reference:'
          CALL OUTPUT(WORK(KDABR),1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
          WRITE(LUPRI,*)'The occupied-occupied block of the reference:'
          CALL OUTPUT(WORK(KDIJR),1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
          WRITE(LUPRI,*)'The occupied-virtual block of the reference:'
          CALL OUTPUT(WORK(KDIAR),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
        END IF
      END IF

*---------------------------------------------------------------------*
*     compute the densities:
*        first do some contributions in a loop over one virtual, then 
*        do the remaining contributions in a loop over two occupieds
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KDAB0),NVIRT*NVIRT)
      CALL DZERO(WORK(KDIJ0),NRHFT*NRHFT)
      CALL DZERO(DAB,NVIRT*NVIRT)
      CALL DZERO(DIJ,NRHFT*NRHFT)
      CALL DZERO(DIA,NVIRT*NRHFT)
      
      KWTILDE = KEND2
      KWBREVE = KWTILDE + NT1AMX*NT1AMX
      KTHEOCC = KWBREVE + NT1AMX*NT1AMX
      KTVIRT1 = KTHEOCC + NT1AMX*NT1AMX
      KTVIRT2 = KTVIRT1 + NT1AMX*NT1AMX
      KEND2   = KTVIRT2 + NT1AMX*NT1AMX

      KT0FMI  = KEND2
      KTYFMI  = KT0FMI + NT1AMX*NRHFT
      KT0MI   = KTYFMI + NT1AMX*NRHFT
      KTYMI   = KT0MI  + NRHFT*NRHFT
      KWDE    = KTYMI  + NRHFT*NRHFT
      KWED    = KWDE   + NT1AMX*NRHFT
      KSCR    = KWED   + NT1AMX*NRHFT
      KEND2   = KSCR   + NT1AMX*NRHFT

      LWRK2   = LWORK   - KEND2
      IF (LWRK2.LT.0)
     &  CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (3b)')

      CALL CCSDT_ADENVIR_NODDY(DAB,DIJ,DIA,WORK(KDAB0),WORK(KDIJ0),
     &                         WORK(KWBAR),WORK(KWYINT),
     &                         WORK(KTHETA),WORK(K3AMP0),
     &                         WORK(KWTILDE),WORK(KWBREVE),
     &                         WORK(KTHEOCC),
     &                         WORK(KTVIRT1),WORK(KTVIRT2),
     &                         WORK(KT0AMP), WORK(KTYAMP), WORK(KL2AMP),
     &                         WORK(KT0FMI), WORK(KTYFMI),
     &                         WORK(KT0MI),  WORK(KTYMI),
     &                         WORK(KWDE),   WORK(KWED), WORK(KSCR))
 
      IF (LOCDBG) THEN
        write(lupri,*)'DAB in noddy after VIR'
        call output(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,lupri)
        write(lupri,*)'DIJ in noddy after VIR'
        call output(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,lupri)
        write(lupri,*)'DIA in noddy after VIR'
        call output(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,lupri)
      END IF

      KWTILDE = KEND1
      KTHVIRT = KWTILDE + NVIRT*NVIRT*NT1AMX
      KSWAP   = KTHVIRT + NVIRT*NVIRT*NT1AMX
      KEND2   = KSWAP   + NVIRT
      LWRK2   = LWORK   - KEND2
      IF (LWRK2.LT.0)
     &  CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (3c)')

      CALL CCSDT_ADENOCC_NODDY(DAB,DIJ,DIA,WORK(KWBAR),WORK(KTHETA),
     &                         WORK(KL2AMP),WORK(KWTILDE),
     &                         WORK(KTHVIRT),WORK(KSWAP))
 
      IF (LOCDBG) THEN
        write(lupri,*)'DAB in noddy after OCC'
        call output(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,lupri)
        write(lupri,*)'DIJ in noddy after OCC'
        call output(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,lupri)
        write(lupri,*)'DIA in noddy after OCC'
        call output(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,lupri)
      END IF

*---------------------------------------------------------------------*
* add the one-index transformed zeroth-order density to DIA:
*---------------------------------------------------------------------*
      KT1AMY = KEND1
      KEND2  = KT1AMY + NT1AMX
      LWRK2  = LWORK   - KEND2
      IF (LWRK2.LT.0)
     &  CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (4)')

      ! read T^Y singles amplitudes from file 
      IOPT = 1
      Call CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,WORK(KT1AMY),DUMMY)

      CALL DGEMM('T','N',NVIRT,NRHFT,NVIRT,
     &           -ONE,WORK(KDAB0),NVIRT,WORK(KT1AMY),NVIRT,
     &            ONE,DIA,NVIRT)
 
      CALL DGEMM('N','T',NVIRT,NRHFT,NRHFT,
     &            ONE,WORK(KT1AMY),NVIRT,WORK(KDIJ0),NRHFT,
     &            ONE,DIA,NVIRT)  
 

      IF (HAVETB3T3) THEN
        CALL DGEMM('T','N',NVIRT,NRHFT,NVIRT,
     &             -ONE,WORK(KDAB0R),NVIRT,WORK(KT1AMY),NVIRT,
     &              ONE,WORK(KDIAR),NVIRT)
 
        CALL DGEMM('N','T',NVIRT,NRHFT,NRHFT,
     &              ONE,WORK(KT1AMY),NVIRT,WORK(KDIJ0R),NRHFT,
     &              ONE,WORK(KDIAR),NVIRT)  
      END IF
*---------------------------------------------------------------------*
* print debug output and return:
*---------------------------------------------------------------------*

      IF (PRINTDENS) THEN
        WRITE(LUPRI,*) 'The virtual-virtual block DAB:'
        CALL OUTPUT(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*) 'The virtual-virtual block of the reference:'
         CALL OUTPUT(WORK(KDABR),1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
        END IF
        WRITE(LUPRI,*) 'The occupied-occupied block DIJ:'
        CALL OUTPUT(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*)'The occupied-occupied block of the reference:'
         CALL OUTPUT(WORK(KDIJR),1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
        END IF
        WRITE(LUPRI,*) 'The occupied-virtual block DIA:'
        CALL OUTPUT(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*)'The occupied-virtual block of the reference:'
         CALL OUTPUT(WORK(KDIAR),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
        END IF
        WRITE(LUPRI,*) 'The virtual-virtual block DAB0:'
        CALL OUTPUT(WORK(KDAB0),1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*) 'The virtual-virtual block of the reference:'
         CALL OUTPUT(WORK(KDAB0R),1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
        END IF
        WRITE(LUPRI,*) 'The occupied-occupied block DIJ0:'
        CALL OUTPUT(WORK(KDIJ0),1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*)'The occupied-occupied block of the reference:'
         CALL OUTPUT(WORK(KDIJ0R),1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
        END IF
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_ADEN_NODDY                      *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_ADENVIR_NODDY(DAB,DIJ,DIA,DAB0,DIJ0,
     &                               WBAR,WYINT,THETA,T0AMP,
     &                             WTILDE,WBREVE,THEOCC,TVIRT1,TVIRT2,
     &                             TXAMP,TYAMP,L2AMP,TXFMI,TYFMI,
     &                             TXMI,TYMI,WDE,WED,SCR)
*---------------------------------------------------------------------*
*     Purpose: set up loop over virtual, sort Wbar, W, theta as they  *
*              would be sorted in the real code and compute the       *
*              contributions to the densities                         *
*                                                                     *
*     Written by Christof Haettig, Jan 2003, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      INTEGER ISYMWB, ISYMTH
      PARAMETER (ISYMWB = 1, ISYMTH = 1)

#if defined (SYS_CRAY)
      REAL DIJ(*), DAB(*), DIA(NVIRT,NRHFT), DIJ0(*), DAB0(*)
      REAL WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WTILDE(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WBREVE(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THEOCC(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TVIRT1(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TVIRT2(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TXAMP(NT1AMX,NVIRT,NRHFT)
      REAL TYAMP(NT1AMX,NVIRT,NRHFT)
      REAL L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TXFMI(NVIRT,NRHFT,NRHFT)
      REAL TYFMI(NVIRT,NRHFT,NRHFT)
      REAL TXMI(NRHFT,NRHFT), TYMI(NRHFT,NRHFT)
      REAL WDE(NVIRT,NRHFT,NRHFT)
      REAL WED(NVIRT,NRHFT,NRHFT)
      REAL SCR(NT1AMX*NRHFT)
      REAL HALF, ONE, ZERO, CME, DDOT
#else
      DOUBLE PRECISION DIJ(*),DAB(*),DIA(NVIRT,NRHFT),DIJ0(*),DAB0(*)
      DOUBLE PRECISION WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WTILDE(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WBREVE(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THEOCC(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TVIRT1(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TVIRT2(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TXAMP(NT1AMX,NVIRT,NRHFT)
      DOUBLE PRECISION TYAMP(NT1AMX,NVIRT,NRHFT)
      DOUBLE PRECISION L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TXFMI(NVIRT,NRHFT,NRHFT) 
      DOUBLE PRECISION TYFMI(NVIRT,NRHFT,NRHFT)
      DOUBLE PRECISION TXMI(NRHFT,NRHFT), TYMI(NRHFT,NRHFT)
      DOUBLE PRECISION WDE(NVIRT,NRHFT,NRHFT)
      DOUBLE PRECISION WED(NVIRT,NRHFT,NRHFT)
      DOUBLE PRECISION SCR(NT1AMX*NRHFT)
      DOUBLE PRECISION HALF, ONE, ZERO, CME, DDOT
#endif
      PARAMETER( HALF = 0.5D0, ONE = 1.0D0, ZERO = 0.0D0 )

      INTEGER NDM, NDL
 
*---------------------------------------------------------------------*
*     start loop over D and E to compute contributions to D(IA) using
*     more or less the same algorithms as in CC3_XI_DEN_IA 
*---------------------------------------------------------------------*
      DO D = 1, NVIRT
       DO L = 1, NRHFT 
        DO E = 1, NVIRT

c          ----------------------------------------------
c          get WDE(fnm) = W^DE(fnml) and WED = W^ED(fnml)
c          ----------------------------------------------
           DO M = 1, NRHFT
             CALL DCOPY(NT1AMX,WBAR(1,1,E,M,D,L),1,WDE(1,1,M),1)
             CALL DCOPY(NT1AMX,WBAR(1,1,D,M,E,L),1,WED(1,1,M),1)
           END DO

c          ----------------------------------------------
c          get TXFMI(fmi) = t^x(fm,Ei) and TYFMI
c          and TXMI(mi)   = t^x(Dm,Ei) and TYMI
c          ----------------------------------------------
           DO I = 1, NRHFT
             CALL DCOPY(NT1AMX,TXAMP(1,E,I),1,TXFMI(1,1,I),1)
             CALL DCOPY(NT1AMX,TYAMP(1,E,I),1,TYFMI(1,1,I),1)
             DO M = 1, NRHFT
               NDM = NVIRT*(M-1) + D
               TXMI(M,I) = TXAMP(NDM,E,I)
               TYMI(M,I) = TYAMP(NDM,E,I)
             END DO
           END DO

c          ------------------------------------------
c          D(IA) -= W^DE(fmnl) t^x(fm,Ei) t^y(an,DL)
c                      scr(n,i) 
c          ------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                ONE,WDE,NT1AMX,TXFMI,NT1AMX,
     &                ZERO,SCR,NRHFT)
           CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &               -ONE,TYAMP(1,D,L),NVIRT,SCR,NRHFT,
     &                ONE,DIA,NVIRT)

c          ------------------------------------------
c          D(IA) -= W^DE(fmnl) t^y(fm,Ei) t^x(an,DL)
c                      scr(n,i) 
c          ------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                ONE,WDE,NT1AMX,TYFMI,NT1AMX,
     &                ZERO,SCR,NRHFT)
           CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &               -ONE,TXAMP(1,D,L),NVIRT,SCR,NRHFT,
     &                ONE,DIA,NVIRT)

c          --------------------------------------------
c          resort:   WDE(fmi) -->   WDE(fim)
c                  TXFMI(fmi) --> TXFMI(fim)
c                  TYFMI(fmi) --> TYFMI(fim)
c          --------------------------------------------
           DO M = 1, NRHFT
             DO I = M+1, NRHFT
               CALL DCOPY(NVIRT,WDE(1,M,I),1,SCR,1)
               CALL DCOPY(NVIRT,WDE(1,I,M),1,WDE(1,M,I),1)
               CALL DCOPY(NVIRT,SCR,1,WDE(1,I,M),1)

               CALL DCOPY(NVIRT,TXFMI(1,M,I),1,SCR,1)
               CALL DCOPY(NVIRT,TXFMI(1,I,M),1,TXFMI(1,M,I),1)
               CALL DCOPY(NVIRT,SCR,1,TXFMI(1,I,M),1)

               CALL DCOPY(NVIRT,TYFMI(1,M,I),1,SCR,1)
               CALL DCOPY(NVIRT,TYFMI(1,I,M),1,TYFMI(1,M,I),1)
               CALL DCOPY(NVIRT,SCR,1,TYFMI(1,I,M),1)
             END DO
           END DO

c          ------------------------------------------
c          D(IA) -= W^DE(fnml) t^x(fi,Em) t^y(an,Dl)
c                      scr(n,i) 
c          ------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                ONE,WDE,NT1AMX,TXFMI,NT1AMX,
     &                ZERO,SCR,NRHFT)
           CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &               -ONE,TYAMP(1,D,L),NVIRT,SCR,NRHFT,
     &                ONE,DIA,NVIRT)



c          ------------------------------------------
c          D(IA) -= W^DE(fnml) t^y(fi,Em) t^x(an,Dl)
c                      scr(n,i) 
c          ------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                ONE,WDE,NT1AMX,TYFMI,NT1AMX,
     &                ZERO,SCR,NRHFT)
           CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &               -ONE,TXAMP(1,D,L),NVIRT,SCR,NRHFT,
     &                ONE,DIA,NVIRT)



c          ------------------------------------------
c          D(IA) -= W^ED(fnml) t^x(Dm,Ei) t^y(fn,al)
c                        scr(fni)
c          ------------------------------------------
           CALL DGEMM('N','N',NT1AMX,NRHFT,NRHFT,
     &                ONE,WED,NT1AMX,TXMI,NRHFT,
     &                ZERO,SCR,NT1AMX)
           CALL DGEMM('T','N',NVIRT,NRHFT,NT1AMX,
     &               -ONE,TYAMP(1,1,L),NT1AMX,SCR,NT1AMX,
     &                ONE,DIA,NVIRT)


c          ------------------------------------------
c          D(IA) -= W^ED(fnml) t^y(Dm,Ei) t^x(fn,al)
c                        scr(fni)
c          ------------------------------------------
           CALL DGEMM('N','N',NT1AMX,NRHFT,NRHFT,
     &                ONE,WED,NT1AMX,TYMI,NRHFT,
     &                ZERO,SCR,NT1AMX)
           CALL DGEMM('T','N',NVIRT,NRHFT,NT1AMX,
     &               -ONE,TXAMP(1,1,L),NT1AMX,SCR,NT1AMX,
     &                ONE,DIA,NVIRT)


        END DO
       END DO 
      END DO 

*---------------------------------------------------------------------*
*     start loop over D and L to compute contributions to D(AB), D(IJ)
*---------------------------------------------------------------------*
      DO D = 1, NVIRT
        DO L = 1, NRHFT
  
c          ----------------------------------------------------
c          build Wtilde^DL(em,an) = W^Da(emnl) + 1/2 W^De(anml)
c          ----------------------------------------------------
           CALL DCOPY(NT1AMX*NT1AMX,WBAR(1,1,1,1,D,L),1,WTILDE,1)
           CALL CC_T2MOD(WTILDE,ISYMWB,HALF)

c          ----------------------------------------------------
c          resort Wbreve^DL(em,fj) = W^Df(emlj) 
c          ----------------------------------------------------
           DO J = 1, NRHFT
             CALL DCOPY(NT1AMX*NVIRT,WBAR(1,1,1,L,D,J),1,
     &                             WBREVE(1,1,1,J),1)
           END DO

c          ----------------------------------------------------
c          fetch TVIRT1^DL(fi,em) = Theta(f-bar i,em,dl)
c          ----------------------------------------------------
           CALL DCOPY(NT1AMX*NT1AMX,THETA(1,1,1,1,D,L),1,TVIRT1,1)

c          --------------------------------------------------------
c          build TVIRT2^DL(em,fi) = Theta(em,fl,di)+Theta(fl,em,di)
c          --------------------------------------------------------
           DO I = 1, NRHFT
            DO F = 1, NVIRT
             DO M = 1, NRHFT
             DO E = 1, NVIRT
              TVIRT2(E,M,F,I) = THETA(E,M,F,L,D,I) + THETA(F,L,E,M,D,I)
             END DO
             END DO
            END DO
           END DO

c          -----------------------------------------------------------
c          build THEOCC^DL(em,fn) = w(em,fn,dl)+w(dl,em,fn)+w(fn,dl,em)
c          -----------------------------------------------------------
           DO E = 1, NVIRT
           DO F = 1, NVIRT
             DO M = 1, NRHFT
             DO N = 1, NRHFT

               THEOCC(E,M,F,N) = WYINT(E,M,F,N,D,L)+
     &                           WYINT(D,L,E,M,F,N)+WYINT(F,N,D,L,E,M)

             END DO
             END DO
           END DO
           END DO



*---------------------------------------------------------------------*
* now we have the following intermediates in O^2 N^2 arrays:
*           
*          WTILDE^DL(em,an) = W^Da(emnl) + 1/2 W^De(anml)
*          WBREVE^DL(em,fj) = W^Df(emlj) 
*          THEOCC^DL(em,fn) = w(em,fn,dl)+w(dl,em,fn)+w(fn,dl,em)
*          TVIRT1^DL(fi,em) = theta(f-bar i,em,dl)
*          TVIRT2^DL(em,fi) = theta(e-bar m,fl,di)+theta(f-bar l,em,di)
*
* with these we can compute the contributions with 3 calls to DGEMM
*---------------------------------------------------------------------*

c          --------------------------------------------
c          D(ab) += WTILDE^DL(em,an) * THEOCC^DL(em,bn)
c          --------------------------------------------
           DO N = 1, NRHFT
             CALL DGEMM('T','N',NVIRT,NVIRT,NT1AMX,
     &                  ONE,WTILDE(1,1,1,N),NT1AMX,
     &                      THEOCC(1,1,1,N),NT1AMX,
     &                  ONE,DAB,NVIRT)
           END DO

c          --------------------------------------------
c          D0(ab) += WTILDE^DL(em,an) * T0^DL(em,bn)
c          --------------------------------------------
           DO N = 1, NRHFT
             CALL DGEMM('T','N',NVIRT,NVIRT,NT1AMX,
     &                  ONE,WTILDE(1,1,1,N),NT1AMX,
     &                       T0AMP(1,1,1,N,D,L),NT1AMX,
     &                  ONE,DAB0,NVIRT)
           END DO

c          -----------------------------------------------------
c          D(ia) += [THEOCC^DL(em,ai)-THEOCC^DL(ei,am)]*L(em,DL)
c          -----------------------------------------------------
           DO I = 1, NRHFT
             DO A = 1, NVIRT
               DO M = 1, NRHFT
                 DIA(A,I) = DIA(A,I) + 
     &              DDOT(NVIRT,THEOCC(1,M,A,I),1,L2AMP(1,M,D,L),1) -
     &              DDOT(NVIRT,THEOCC(1,I,A,M),1,L2AMP(1,M,D,L),1)
              END DO
             END DO
           END DO

c          ----------------------------------------------------------
c          add TVIRT1 to THEOCC: THEOCC^DL(em,fi) += TVIRT1^DL(fi,em)
c          ----------------------------------------------------------
           DO M = 1, NRHFT
             DO E = 1, NVIRT
               DO I = 1, NRHFT
                 DO F = 1, NVIRT
                   THEOCC(E,M,F,I) = THEOCC(E,M,F,I)+TVIRT1(F,I,E,M)
                 END DO
               END DO
             END DO
           END DO

c          -----------------------------------------------------------
c          D(ij)+=WTILDE^DL(em,fj)*[THEOCC^DL(em,fi)+TVIRT1^DL(fi,em)]
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                          -ONE,THEOCC,NT1AMX*NVIRT,
     &                             WTILDE,NT1AMX*NVIRT,
     &                           ONE,DIJ,NRHFT)

c          -----------------------------------------------------------
c          D(ij) += WBREVE^DL(em,fj) * TVIRT2^DL(em,fi)
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                          -ONE,TVIRT2,NT1AMX*NVIRT,
     &                               WBREVE,NT1AMX*NVIRT,
     &                           ONE,DIJ,NRHFT)

c          -----------------------------------------------------------
c          D0(ij) += WTILDE^DL(em,fj) * T0AMP^DL(em,fi)
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                          -ONE,T0AMP(1,1,1,1,D,L),NT1AMX*NVIRT,
     &                             WTILDE,NT1AMX*NVIRT,
     &                           ONE,DIJ0,NRHFT)


        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_ADENVIR_NODDY                   *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_ADENOCC_NODDY(DAB,DIJ,DIA,WBAR,THETA,L2AMP,
     &                               WTILDE,THVIRT,SWAP)
*---------------------------------------------------------------------*
*     Purpose: set up loop over occupied, sort Wbar, W, theta as they *
*              would be sorted in the real code and compute the       *
*              contributions to the densities                         *
*                                                                     *
*     Written by Christof Haettig, Jan 2003, Friedrichstal            *
*---------------------------------------------------------------------*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      REAL DIJ(*), DAB(*), DIA(NVIRT,NRHFT)
      REAL WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WTILDE(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL THVIRT(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL SWAP(NVIRT)
      REAL HALF, ONE, DDOT
#else
      DOUBLE PRECISION DIJ(*), DAB(*), DIA(NVIRT,NRHFT)
      DOUBLE PRECISION WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WTILDE(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION THVIRT(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION SWAP(NVIRT)
      DOUBLE PRECISION HALF, ONE, DDOT
#endif
      PARAMETER( HALF = 0.5D0, ONE = 1.0D0 )

*---------------------------------------------------------------------*
*     start loop over L and M
*---------------------------------------------------------------------*
      DO L = 1, NRHFT
        DO M = 1, NRHFT

c          ----------------------------------------------------
c          resort WTILDE^LM(de,an) = W^DE(anml) = W^LM(naed)
c          ----------------------------------------------------
           DO N = 1, NRHFT
             DO A = 1, NVIRT
               DO E = 1, NVIRT
                 DO D = 1, NVIRT
                   WTILDE(D,E,A,N) = WBAR(A,N,E,M,D,L)
                 END DO
               END DO
             END DO
           END DO
        
c          ------------------------------------------------
c          resort THVIRT^LM(de,an) = theta(d-bar l,em,an) + 
c                                    theta(dl,e-bar m,an) 
c          ------------------------------------------------
           DO N = 1, NRHFT
            DO A = 1, NVIRT
             DO E = 1, NVIRT
              DO D = 1, NVIRT
               THVIRT(D,E,A,N) = THETA(D,L,A,N,E,M)+ THETA(E,M,D,L,A,N)
              END DO
             END DO
            END DO
           END DO

c          --------------------------------------------
c          D(ij) += WTILDE^LM(de,fj) * THVIRT^LM(de,fi)
c          --------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NVIRT*NVIRT*NVIRT,
     &                       -HALF,THVIRT,NVIRT*NVIRT*NVIRT,
     &                             WTILDE,NVIRT*NVIRT*NVIRT,
     &                         ONE,DIJ,NRHFT)
C
c          --------------------------------------------
c          add THVIRT^LM(de,an) += theta(dl,em,a-bar n)
c          --------------------------------------------
           DO N = 1, NRHFT
            DO A = 1, NVIRT
             DO E = 1, NVIRT
              DO D = 1, NVIRT
               THVIRT(D,E,A,N) = THVIRT(D,E,A,N) + THETA(A,N,E,M,D,L) 
              END DO
             END DO
            END DO
           END DO
        
c          --------------------------------------------
c          D(ab) += WTILDE^LM(de,an) * THVIRT^LM(de,bn)
c          --------------------------------------------
           DO N = 1, NRHFT
             CALL DGEMM('T','N',NVIRT,NVIRT,NVIRT*NVIRT,
     &                 HALF,WTILDE(1,1,1,N),NVIRT*NVIRT,
     &                      THVIRT(1,1,1,N),NVIRT*NVIRT,
     &                  ONE,DAB,NVIRT)
           END DO

c          -----------------------------------------------------
c          D(ia) += [THVIRT^LM(de,ai)-THVIRT^LM(da,ei)]*L(dL,eM)
c          -----------------------------------------------------
           DO I = 1, NRHFT
             DO A = 1, NVIRT
               DO E = 1, NVIRT
                 DIA(A,I) = DIA(A,I) + 
     &              DDOT(NVIRT,THVIRT(1,E,A,I),1,L2AMP(1,L,E,M),1) -
     &              DDOT(NVIRT,THVIRT(1,A,E,I),1,L2AMP(1,L,E,M),1)
              END DO
             END DO
           END DO


c          --------------------------------------------
c          resort:   WTILDE^LM(de,an) --> WTILDE(da,en)
c                    THVIRT^LM(de,an) --> THVIRT(da,en)
c          --------------------------------------------
           DO N = 1, NRHFT
             DO A = 1, NVIRT
               DO E = A+1, NVIRT
                 CALL DCOPY(NVIRT,WTILDE(1,E,A,N),1,SWAP,1)
                 CALL DCOPY(NVIRT,WTILDE(1,A,E,N),1,WTILDE(1,E,A,N),1)
                 CALL DCOPY(NVIRT,SWAP,1,WTILDE(1,A,E,N),1)

                 CALL DCOPY(NVIRT,THVIRT(1,E,A,N),1,SWAP,1)
                 CALL DCOPY(NVIRT,THVIRT(1,A,E,N),1,THVIRT(1,E,A,N),1)
                 CALL DCOPY(NVIRT,SWAP,1,THVIRT(1,A,E,N),1)
               END DO
             END DO
           END DO

c          --------------------------------------------
c          D(ab) += WTILDE^LM(da,en) * THVIRT^LM(db,en)
c          --------------------------------------------
           DO N = 1, NRHFT
             CALL DGEMM('T','N',NVIRT,NVIRT,NVIRT*NVIRT,
     &                  ONE,WTILDE(1,1,1,N),NVIRT*NVIRT,
     &                      THVIRT(1,1,1,N),NVIRT*NVIRT,
     &                  ONE,DAB,NVIRT)
           END DO


        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_ADENOCC_NODDY                   *
*---------------------------------------------------------------------*
*---------------------------------------------------------------------*
C  /* Deck ccsdt_gwbic */
*=====================================================================*
      SUBROUTINE CCSDT_GWBIC(LISTL,IDLSTL,WBAR,
     &                       XLAMDP0,XLAMDH0,FOCK0,
     &                       LUTOC,FNTOC,LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     &                       LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     &                       WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: set up WBAR as N^6 array in core for use in noddy code
*             
*    W3BAR^Y = ( <L2|[Y,tau3]|HF> + <L3|[Y^,tau3]|HF>
*            +   <L2|[H^Y,tau3]|HF> 
*            +   <L1^Y + L2^Y|[H^,tau3]|HF> )/ (w+epsiln(tau3))
*
*
*    Christof Haettig, spring 2003 based CC3_XI_DEN.
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
#include <cc3t3d.h>
#include <ccl1rsp.h>
#include <ccr1rsp.h>

C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
      CHARACTER LISTL*3, LISTR*3, MODEL*10, LABELY*8

      CHARACTER*5 FN3FOP
      CHARACTER*6 FN3FOP2
      CHARACTER*8 FN3VI2

      PARAMETER (FN3VI2 = 'CC3_VI12')
      PARAMETER (FN3FOP = 'PTFOP'  , FN3FOP2= 'PTFOP2')

      CHARACTER*(*) FNTOC,FN3VI,FNDKBC3,FN3FOPX,FN3FOP2X

      LOGICAL LORX
      INTEGER LWORK, IDLSTL, ISYMY, IDLSTR, LU3VI2, LU3FOPX, LU3FOP2X,
     &        LUTOC, LU3VI, LUDKBC3, LU3FOP, LU3FOP2,
     &        KFOCKD,KFCKBA,KL1AM,KL2TP,KEND0,LWRK0,
     &        KL1BY, KL2TPBY, KFOCKY, KEND1, LWRK1,
     &        IOPT, LENGTH, IOFF, IRREP, ISYM, IERR,
     &        ISYMD, ISYML, ISYMDL, ISYMC, ISYMK, KOFF1, KOFF2, KXIAJB,
     &        KTROC01, KTROC21, KTROC03, KTROC23, KT3OC0, KT3OC02,
     &        KLAMPY, KLAMHY, KFCKYCK, KINTOC, KT1Y,
     &        KTROCCG0,KTROCCL0,KTROCCGY,KTROCCLY, KEND2,LWRK2,
     &        ISYCKB, ISYCKBY,
     &        KTRVI4, KTRVI5, KTRVI6, KTRVI7, KEND3, LWRK3,
     &        KTRVI11,KTRVI12,KTRVI13,
     &        KTRVI14,KTRVI15,KTRVI16,KTRVI17,KTRVI18,KTRVI19,KTRVI20,
     &        KVGDKCB0,KVGDKBC0,KVLDKCB0,KVLDKBC0, KINTVI,
     &        KVGDKCBY,KVGDKBCY,KVLDKCBY,KVLDKBCY, KEND4, LWRK4,
     &        ISYALJ, ISYALJ2, ISYALJBY, ISYALJDY, ISYMBD, ISCKIJ,
     &        ISWMAT, ISYCKD, ISAIBJ, ISYMJ, ISYMAI, ISYAIL,
     &        KSMAT2, KUMAT2, KDIAG, KDIAGW, KINDSQ, KINDSQW, KINDEX,
     &        KINDEX2, KINDEXBY, KINDEXDY, KTMAT,KWMAT, KSMAT4,KUMAT4,
     &        KEND5, LWRK5, LENSQ, LENSQW, ISYMB, ISYMBJ

#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      REAL FREQ, HALF
      REAL WORK(LWORK)
      REAL WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      DOUBLE PRECISION FREQ, HALF
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
#endif
      PARAMETER(HALF = 0.5D0)

* external functions:
      INTEGER IR1TAMP

*---------------------------------------------------------------------*
* some initializations:
*---------------------------------------------------------------------*
      CALL QENTER('CCSDT_GWBIC')

      CALL DZERO(WBAR,NT1AMX*NT1AMX*NT1AMX)

      IF (LISTL(1:3).EQ.'L1 ') THEN
         ! get symmetry, frequency and integral label from common blocks
         ! defined in ccl1rsp.h
         ISYMY   = ISYLRZ(IDLSTL)
         FREQ    = FRQLRZ(IDLSTL)
         LABELY  = LRZLBL(IDLSTL)
         LORX    = LORXLRZ(IDLSTL)
         IF (LORX) CALL QUIT('NO ORBITAL RELAX. IN CCSDT_GWBIC')
      ELSE
        CALL QUIT('Illegal left vector type in CCSDT_GWBIC:'//LISTL)
      END IF

      ! find right response vector:
      IF (LISTL(1:3).EQ.'L1 ') THEN
        LISTR  = 'R1 '
        IDLSTR = IR1TAMP(LABELY,LORX,FREQ,ISYMY)
      ELSE
        CALL QUIT('Unkown list in CCSDT_GWBIC')
      END IF

*---------------------------------------------------------------------*
* open files:
*---------------------------------------------------------------------*
      LU3VI2  = 0
      LU3FOP  = 0
      LU3FOP2 = 0

      CALL WOPEN2(LU3VI2,FN3VI2,64,0) 
      CALL WOPEN2(LU3FOP,FN3FOP,64,0) 
      CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)

*---------------------------------------------------------------------*
* initial allocations, orbital energy, fock matrix and L2 :
*---------------------------------------------------------------------*
      KFOCKD  = 1
      KFCKBA  = KFOCKD  + NORBTS
      KL1AM   = KFCKBA  + NT1AMX 
      KL2TP   = KL1AM   + NT1AM(ISYM0)
      KEND0   = KL2TP   + NT2SQ(ISYM0)
      LWRK0   = LWORK   - KEND0
 
      KL1BY   = KEND0
      KL2TPBY = KL1BY   + NT1AM(ISYMY)
      KFOCKY  = KL2TPBY   + NT2SQ(ISYMY)
      KEND0   = KFOCKY    + N2BST(ISYMY)
 
      KEND1   = KEND0
      LWRK1   = LWORK   - KEND1

C---------------------------------------------------------------------*
C     Matrix with property integrals in MO basis:
C---------------------------------------------------------------------*
      ! read property integrals from file:
      CALL CCPRPAO(LABELY,WORK(KFOCKY),IRREP,ISYM,IERR,
     &             WORK(KEND1),LWRK1)
      IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMY)) THEN
          CALL QUIT('CCSDT_GWBIC: error reading oper. '//LABELY)
      ELSE IF (IERR.LT.0) THEN
          CALL DZERO(WORK(KFOCKY),N2BST(ISYMY))
      END IF

      ! transform property integrals to Lambda-MO basis
      CALL CC_FCKMO(WORK(KFOCKY),XLAMDP0,XLAMDH0,
     &              WORK(KEND1),LWRK1,ISYMY,1,1)

C-------------------------------------
C     Read L2 amplitudes 
C-------------------------------------
      IOPT = 3
      CALL CC_RDRSP('L0 ',0,ISYM0,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KL2TP))
      CALL CC_T2SQ(WORK(KL2TP),WORK(KEND1),ISYM0)
      CALL CC3_T2TP(WORK(KL2TP),WORK(KEND1),ISYM0)

C-------------------------------------
C     Read L1BY and L2BY multipliers 
C-------------------------------------
      IOPT  = 3
      CALL CC_RDRSP(LISTL,IDLSTL,ISYMY,IOPT,MODEL,
     &              WORK(KL1BY),WORK(KL2TPBY))
      CALL CC_T2SQ(WORK(KL2TPBY),WORK(KEND1),ISYMY)
      CALL CC3_T2TP(WORK(KL2TPBY),WORK(KEND1),ISYMY)

C-------------------------------------
C     Read canonical orbital energies:
C-------------------------------------
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
      CALL GPCLOSE(LUSIFC,'KEEP')
 
      ! Delete frozen orbitals in Fock diagonal.
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK0)

C-----------------------------------------------------
C     Sort the fock matrix
C-----------------------------------------------------
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
C
C-----------------------------
C     Memory allocation.
C-----------------------------
      KXIAJB   = KEND1
      KEND1   = KXIAJB  + NT2AM(ISYM0)

      KTROC01 = KEND1
      KTROC21 = KTROC01 + NTRAOC(ISYM0)
      KTROC03 = KTROC21 + NTRAOC(ISYM0)
      KTROC23 = KTROC03 + NTRAOC(ISYM0)
      KT3OC0  = KTROC23 + NTRAOC(ISYM0)
      KT3OC02 = KT3OC0  + NTRAOC(ISYM0)
      KLAMPY  = KT3OC02 + NTRAOC(ISYM0)
      KLAMHY  = KLAMPY  + NLAMDT
      KEND1   = KLAMHY  + NLAMDT
      LWRK1   = LWORK   - KEND1

      KFCKYCK    = KEND1
      KTROCCG0   = KFCKYCK    + NT1AM(ISYMY)
      KTROCCL0   = KTROCCG0   + NTRAOC(ISYM0)
      KTROCCGY   = KTROCCL0   + NTRAOC(ISYM0)
      KTROCCLY   = KTROCCGY   + NTRAOC(ISYMY)
      KEND1      = KTROCCLY   + NTRAOC(ISYMY)
      LWRK1      = LWORK      - KEND1

      KINTOC  = KEND1
      KT1Y    = KINTOC + MAX(NTOTOC(ISYM0),NTOTOC(ISYM0))
      KEND2   = KT1Y + NT1AM(ISYMY)
      LWRK2   = LWORK  - KEND2

      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CCSDT_GWBIC')
      END IF
 
C------------------------
C     Construct L(ia,jb).
C------------------------
      LENGTH = IRAT*NT2AM(ISYM0)
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYM0,1)
 
C------------------------------------------------------------
C     Create Lambda-p^Y and Lambda-h^Y:
C------------------------------------------------------------
      IOPT = 1
      CALL CC_RDRSP(LISTR,IDLSTR,ISYMY,IOPT,MODEL,WORK(KT1Y),DUMMY)
 
      CALL CCLR_LAMTRA(XLAMDP0,WORK(KLAMPY),XLAMDH0,
     *                 WORK(KLAMHY),WORK(KT1Y),ISYMY)

C----------------------------------------------------------------------
C     Calculate the F^Y matrix (kc elements evaluated and stored as ck)
C----------------------------------------------------------------------
      CALL CC3LR_MFOCK(WORK(KFCKYCK),WORK(KT1Y),WORK(KXIAJB),ISYMY)
 
C-----------------------------------------------------------
C     Construct 2*C-E of the integrals.
C     Have integral for both (ij,k,a) and (a,k,j,i)
C-----------------------------------------------------------
      IOFF = 1
      IF (NTOTOC(ISYM0) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fntoc (1)'
cch
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF
 
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC01),XLAMDH0,
     *               WORK(KEND2),LWRK2,ISYM0)
 
      CALL CCSDT_TCMEOCC(WORK(KTROC01),WORK(KTROC21),ISYM0)
 
      CALL CCFOP_SORT(WORK(KTROC01),WORK(KTROC03),ISYM0,1)
 
      CALL CCFOP_SORT(WORK(KTROC21),WORK(KTROC23),ISYM0,1)

C     -------------------------------
C     Read in integrals for WMAT etc.
C     -------------------------------
      IOFF = 1
      IF (NTOTOC(ISYMOP) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fntoc (2)'
cch
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
      ENDIF
 
C     -----------------
C     LAMPY TRANSFORMED
C     -----------------
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROCCGY),
     *                WORK(KLAMHY),ISYMY,WORK(KEND2),LWRK2)
      CALL CC3_SRTOCC(WORK(KTROCCGY),WORK(KTROCCLY),ISYMY)
 
C     -----------------
C     LAMP0 TRANSFORMED
C     -----------------
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROCCG0),
     *                XLAMDH0,ISYM0,WORK(KEND2),LWRK2)
      CALL CC3_SRTOCC(WORK(KTROCCG0),WORK(KTROCCL0),ISYM0)
C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKB  = MULD2H(ISYMD,ISYM0)
         ISYCKBY = MULD2H(ISYMD,ISYMY)
 
         DO D = 1,NVIR(ISYMD)
 
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

            KVGDKCB0  = KEND3
            KVGDKBC0  = KVGDKCB0  + NCKATR(ISYCKB)
            KVLDKCB0  = KVGDKBC0  + NCKATR(ISYCKB)
            KVLDKBC0  = KVLDKCB0  + NCKATR(ISYCKB)
            KEND3     = KVLDKBC0  + NCKATR(ISYCKB)
            LWRK3     = LWORK     - KEND3

            KVGDKCBY  = KEND3
            KVGDKBCY  = KVGDKCBY  + NCKATR(ISYCKBY)
            KVLDKCBY  = KVGDKBCY  + NCKATR(ISYCKBY)
            KVLDKBCY  = KVLDKCBY  + NCKATR(ISYCKBY)
            KEND3     = KVLDKBCY  + NCKATR(ISYCKBY)
            LWRK3     = LWORK     - KEND3
           
            KINTVI = KEND3
            KTRVI6 = KINTVI + MAX(NCKA(ISYCKB),NCKA(ISYCKBY))
            KEND4  = KTRVI6 + NCKATR(ISYCKB) 
            LWRK4  = LWORK  - KEND4
           
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CCSDT_GWBIC')
            END IF
C
C           -------------------------------------------
C           Read 2*C-E of integral used for t3-bar
C           -------------------------------------------
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fn3fop2x (3)'
cch
               CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI4),IOFF,
     &                     NCKATR(ISYCKB))
            ENDIF
C
C           -------------------------------------------------
C           Integrals used for t3-bar for cc3
C           -------------------------------------------------
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fndkbc3 (4)'
cch
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
            CALL CCSDT_SRTVIR(WORK(KTRVI4),WORK(KTRVI5),WORK(KEND4),
     *                        LWRK4,ISYMD,ISYM0)
C
C           ----------------------------------------------------
C           Read virtual integrals used in q3am/u3am for t3-bar.
C           ----------------------------------------------------
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fn3vi (4)'
cch
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF

            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI19),XLAMDP0,
     *                       ISYMD,D,ISYM0,WORK(KEND4),LWRK4)

            IF (LWRK4 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     *                   'CCSDT_GWBIC  (CC3 TRVI)')
            END IF

            CALL CCSDT_SRTVIR(WORK(KTRVI19),WORK(KTRVI18),WORK(KEND4)
     *                        ,LWRK4,ISYMD,ISYM0)

            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
cch
C              write(lupri,*) 'gwic> 1: ioff,len:',ioff,nckatr(isyckb)
cch
               CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI6),IOFF,
     *                     NCKATR(ISYCKB))
            ENDIF

            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDT_GWBIC (2)')
            END IF

            CALL DCOPY(NCKATR(ISYCKB),WORK(KTRVI6),1,WORK(KTRVI7),1)
C
C------------------------------------------------------------
C           Read in, transform and sort integrals used in 
C           construction of W intermediate.
C------------------------------------------------------------
 
C           -------------------------------------------
C           Read in g(c1^h k1^p | delta b2^h) integrals
C           -------------------------------------------
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fn3vi (5)'
cch
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF
 
C           -------------------------------------------------
C           LAMPY TRANSFORMED:
C           transform g(c1^h k1^p | delta b2^h) integrals
C              to g(c1^h k1^p | d2^pY b2^h) with lampY
C           and sort as g(d2^pY k1^p | c1^h b2^h)
C           -------------------------------------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKCBY),
     *                      WORK(KLAMPY),ISYMY,
     *                      ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVGDKCBY),WORK(KEND4),ISYMD,D,ISYMY)

C           -------------------------------------------------
C           LAMP0 TRANSFORMED:
C           transform g(c1^h k1^p | delta b2^h) integrals
C              to g(c1^h k1^p | d2^p b2^h) with lamp0
C           and sort as g(d2^p k1^p | c1^h b2^h)
C           -------------------------------------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKCB0),
     *                      XLAMDP0,ISYM0,
     *                      ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVGDKCB0),WORK(KEND4),ISYMD,D,ISYM0)

C           -------------------------------------------------
C           Read in g(b2^h k1^p | delta c1^h) integrals and 
C           transform and sort --- see above
C           -------------------------------------------------
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fn3vi2 (6)'
cch
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
            ENDIF

C
C           ------------------
C           LAMPY TRANSFORMED:
C           ------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKBCY),
     *                      WORK(KLAMPY),ISYMY,
     *                      ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVGDKBCY),WORK(KEND4),ISYMD,D,ISYMY)

C           ------------------
C           LAMP0 TRANSFORMED:
C           ------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVGDKBC0),XLAMDP0,ISYM0,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVGDKBC0),WORK(KEND4),ISYMD,D,ISYM0)

C           --------------------------------------------
C           Read in L(c1^h k1^p | delta b2^h) integrals:
C           --------------------------------------------
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
cch
C           write(lupri,*) 'read from fn3fop (7):',lu3fop,fn3fop
cch
              CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,NCKA(ISYCKB))
            ENDIF
 
C           ---------------------------------------------
C           LAMPY TRANSFORMED:
C           transform L(c1^h k1^p | delta b2^h) integrals
C              to L(c1^h k1^p | d2^pY b2^h) 
C           and sort as L(d2^pY k1^p | c1^h b2^h)
C           ---------------------------------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKCBY),
     *                          WORK(KLAMPY),ISYMY,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVLDKCBY),WORK(KEND4),ISYMD,D,ISYMY)

C           ---------------------------------------------
C           LAMP0 TRANSFORMED
C           transform L(c1^h k1^p | delta b2^h) integrals
C              to L(c1^h k1^p | d2^p b2^h) 
C           and sort as L(d2^p k1^p | c1^h b2^h)
C           ---------------------------------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKCB0),
     *                      XLAMDP0,ISYM0,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVLDKCB0),WORK(KEND4),ISYMD,D,ISYM0)
 
C           -----------------------------------------------
C           Read in L(b2^h k1^p | delta c1^h) integrals and 
C           transform and sort --- see above
C           -----------------------------------------------
cch
C        write(lupri,*) 'read from fn3fop2 (8)'
cch
            CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
 
C           ------------------
C           LAMPY TRANSFORMED:
C           ------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKBCY),
     *                          WORK(KLAMPY),ISYMY,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVLDKBCY),WORK(KEND4),ISYMD,D,ISYMY)
 
C           ------------------
C           LAMP0 TRANSFORMED:
C           ------------------
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KVLDKBC0),
     *                      XLAMDP0,ISYM0,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            CALL CCSDT_SRVIR3(WORK(KVLDKBC0),WORK(KEND4),ISYMD,D,ISYM0)
 

            DO ISYMB = 1,NSYM

               ISYALJ  = MULD2H(ISYMB,ISYM0)
               ISYALJ2 = MULD2H(ISYMD,ISYM0)
               ISYALJBY= MULD2H(ISYMB,ISYMY)
               ISYALJDY= MULD2H(ISYMD,ISYMY)
               ISYMBD  = MULD2H(ISYMD,ISYMB)
               ISCKIJ  = MULD2H(ISYMBD,ISYM0)
               ISWMAT  = MULD2H(ISCKIJ,ISYMY)
               ISYCKD  = MULD2H(ISYM0,ISYMB)

C              Can use kend3 since we do not need the integrals anymore.
               KSMAT2     = KEND3
               KUMAT2     = KSMAT2    + NCKIJ(ISCKIJ)
               KDIAG      = KUMAT2    + NCKIJ(ISCKIJ)
               KDIAGW     = KDIAG     + NCKIJ(ISCKIJ)
               KINDSQ     = KDIAGW    + NCKIJ(ISWMAT)
               KINDSQW    = KINDSQ    + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX     = KINDSQW   + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEX2    = KINDEX    + (NCKI(ISYALJ)  - 1)/IRAT + 1
               KINDEXBY   = KINDEX2   + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KINDEXDY   = KINDEXBY  + (NCKI(ISYALJBY)  - 1)/IRAT + 1
               KTMAT      = KINDEXDY  + (NCKI(ISYALJDY) - 1)/IRAT + 1
               KWMAT      = KTMAT     + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWMAT))
               KEND4      = KWMAT     + NCKIJ(ISWMAT)
               LWRK4      = LWORK     - KEND4

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
                  CALL QUIT('Insufficient space in CCSDT_GWBIC')
               END IF
 
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
               CALL CC3_DIAG(WORK(KDIAG), WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)
 
C              -----------------------
C              Construct index arrays.
C              -----------------------
               LENSQ  = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               LENSQW = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT)
               CALL CC3_INDEX(WORK(KINDEX),  ISYALJ)
               CALL CC3_INDEX(WORK(KINDEX2), ISYALJ2)
               CALL CC3_INDEX(WORK(KINDEXBY),ISYALJBY)
               CALL CC3_INDEX(WORK(KINDEXDY),ISYALJDY)

               DO B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))

C                 ------------------------------------
C                 Read and transform integrals used in 
C                 first S-bar and U-bar
C                 ------------------------------------
                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B - 1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fndkbc3 (9)'
cch
                     CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI16),IOFF,
     *                           NCKATR(ISYCKD))
                  ENDIF
                  CALL CCSDT_SRVIR3(WORK(KTRVI16),WORK(KEND5),
     *                              ISYMB,B,ISYM0)
                  CALL CCSDT_SRTVIR(WORK(KTRVI16),WORK(KTRVI17),
     *                              WORK(KEND5),LWRK5,ISYMB,ISYM0)
C
C                 ------------------------------------
C                 Read and transform integrals used in 
C                 second S-bar and U-bar
C                 ------------------------------------
                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B-1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fn3fop2x (10)'
cch
                     CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI11),
     *                           IOFF,NCKATR(ISYCKD))
                  ENDIF

                  CALL CCSDT_SRTVIR(WORK(KTRVI11),WORK(KTRVI12),
     *                              WORK(KEND5),LWRK5,ISYMB,
     *                              ISYM0)

                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B - 1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
cch
C              write(lupri,*) 'gwic> 2: ioff,len:',ioff,nckatr(isyckd)
cch
                     CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI13),IOFF,
     *                           NCKATR(ISYCKD))
                  ENDIF

                  IOFF = ICKAD(ISYCKD,ISYMB) 
     *                 + NCKA(ISYCKD)*(B - 1) + 1
                  IF (NCKA(ISYCKD) .GT. 0) THEN
cch
C        write(lupri,*) 'read from fn3vi (10)'
cch
                     CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                           NCKA(ISYCKD))
                  ENDIF

                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI20),
     *                             XLAMDP0,ISYMB,B,ISYM0,
     *                             WORK(KEND5),LWRK5)

C                 ----------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for B,D for T3-BAR:
C                 ----------------------------------------------------
                  CALL DZERO(WORK(KSMAT2),NCKIJ(ISCKIJ))
 
                  CALL CCFOP_SMAT(0.0D0,WORK(KL1AM),ISYM0,WORK(KL2TP),
     *                            ISYM0,WORK(KTMAT),
     *                            WORK(KFCKBA),WORK(KXIAJB),ISYM0,
     *                            WORK(KTRVI14),WORK(KTRVI15),
     *                            WORK(KTRVI4),WORK(KTRVI5),
     *                            WORK(KTROC01),WORK(KTROC21),
     *                            ISYM0,WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KSMAT2),WORK(KEND5),LWRK5,
     *                            WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                            ISYMB,B,ISYMD,D)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT2),1)

                  CALL T3_FORBIDDEN(WORK(KSMAT2),ISYM0,ISYMB,B,ISYMD,D)
C
C                 ----------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for D,B for T3-BAR:
C                 ----------------------------------------------------
                  CALL DZERO(WORK(KSMAT4),NCKIJ(ISCKIJ))

                  CALL CCFOP_SMAT(0.0D0,WORK(KL1AM),ISYM0,WORK(KL2TP),
     *                            ISYM0,WORK(KTMAT),WORK(KFCKBA),
     *                            WORK(KXIAJB),ISYM0,
     *                            WORK(KTRVI16),WORK(KTRVI17),
     *                            WORK(KTRVI11),WORK(KTRVI12),
     *                            WORK(KTROC01),WORK(KTROC21),
     *                            ISYM0,WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KSMAT4),WORK(KEND5),LWRK5,
     *                            WORK(KINDEX2),WORK(KINDSQ),
     *                            LENSQ,ISYMD,D,ISYMB,B)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT4),1)
                  CALL T3_FORBIDDEN(WORK(KSMAT4),ISYM0,ISYMD,D,ISYMB,B)
C
C                 ------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d for t3-bar.
C                 ------------------------------------------------
                  CALL DZERO(WORK(KUMAT2),NCKIJ(ISCKIJ))

                  CALL CCFOP_UMAT(0.0D0,WORK(KL1AM),ISYM0,WORK(KL2TP),
     *                            ISYM0,
     *                            WORK(KXIAJB),ISYM0,WORK(KFCKBA),
     *                            WORK(KTRVI19),WORK(KTRVI7),
     *                            WORK(KTROC03),WORK(KTROC23),ISYM0,
     *                            WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KUMAT2),
     *                            WORK(KTMAT),WORK(KEND5),LWRK5,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KUMAT2),1)
                  CALL T3_FORBIDDEN(WORK(KUMAT2),ISYM0,ISYMB,B,ISYMD,D)
C
C                 ------------------------------------------------
C                 Calculate U(ci,jk) for fixed d,b for t3-bar.
C                 ------------------------------------------------
                  CALL DZERO(WORK(KUMAT4),NCKIJ(ISCKIJ))

                  CALL CCFOP_UMAT(0.0D0,WORK(KL1AM),ISYM0,WORK(KL2TP),
     *                            ISYM0,WORK(KXIAJB),ISYM0,
     *                            WORK(KFCKBA),WORK(KTRVI20),
     *                            WORK(KTRVI13),WORK(KTROC03),
     *                            WORK(KTROC23),ISYM0,
     *                            WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KUMAT4),WORK(KTMAT),
     *                            WORK(KEND5),LWRK5,WORK(KINDSQ),
     *                            LENSQ,ISYMD,D,ISYMB,B)

                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KUMAT4),1)
                  CALL T3_FORBIDDEN(WORK(KUMAT4),ISYM0,ISYMD,D,ISYMB,B)
C
C                 -------------------------------------------
C                 Sum up the S-bar and U-bar to get a real T3-bar
C                 -------------------------------------------
                  CALL CC3_T3BD(ISCKIJ,WORK(KSMAT2),WORK(KSMAT4),
     *                                 WORK(KUMAT2),WORK(KUMAT4),
     *                                 WORK(KTMAT),WORK(KINDSQ),LENSQ)
C
C                 ----------------------------
C                 Calculate <L3|[Y^,tau3]|HF>:
C                 ----------------------------
                  CALL WBARBD_V(WORK(KTMAT),ISCKIJ,WORK(KFOCKY),ISYMY,
     *                 WORK(KWMAT),ISWMAT,WORK(KEND5),LWRK5)
 
                  CALL WBARBD_O(WORK(KTMAT),ISCKIJ,WORK(KFOCKY),ISYMY,
     *                 WORK(KWMAT),ISWMAT,WORK(KEND5),LWRK5)
 
C                 ---------------------------
C                 Calculate <L2|[Y,tau3]|HF>:
C                 ---------------------------
                  CALL WBARBD_T2(B,ISYMB,D,ISYMD,WORK(KL2TP),ISYM0,
     *                           WORK(KFOCKY),ISYMY,WORK(KWMAT),ISWMAT)
 
C                 -----------------------------
C                 calculate <L2|[H^Y,tau3]|HF>:
C                 -----------------------------
                  CALL WBARBD_TMAT(WORK(KL2TP),ISYM0,WORK(KWMAT),
     *                             WORK(KTMAT),ISWMAT,WORK(KFCKYCK),
     *                             ISYMY,WORK(KVLDKBCY),
     *                             WORK(KVLDKCBY),
     *                             WORK(KVGDKBCY),WORK(KVGDKCBY),
     *                             WORK(KTROCCLY),WORK(KTROCCGY),
     *                             ISYMY,WORK(KEND5),LWRK5,
     *                             WORK(KINDEX),WORK(KINDEX2),
     *                             WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)

C                 -----------------------------
C                 Calculate <L2Y|[H^,tau3]|HF>:
C                 -----------------------------
                  CALL WBARBD_TMAT(WORK(KL2TPBY),ISYMY,WORK(KWMAT),
     *                             WORK(KTMAT),ISWMAT,WORK(KFCKBA),
     *                             ISYM0,WORK(KVLDKBC0),
     *                             WORK(KVLDKCB0),
     *                             WORK(KVGDKBC0),WORK(KVGDKCB0),
     *                             WORK(KTROCCL0),WORK(KTROCCG0),
     *                             ISYM0,WORK(KEND5),LWRK5,
     *                             WORK(KINDEXBY),WORK(KINDEXDY),
     *                             WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)

C                 -----------------------------
C                 Calculate <L1Y|[H^,tau3]|HF>:
C                 -----------------------------
                  CALL WBARBD_L1(WORK(KL1BY),ISYMY,WORK(KTMAT),
     *                           WORK(KXIAJB),ISYM0,
     *                           WORK(KWMAT),WORK(KEND5),LWRK5,
     *                           WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)
 
C                 -----------------------------------
C                 Divide by the energy difference and
C                 remove the forbidden elements
C                 -----------------------------------
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,-FREQ,
     *                         ISWMAT,WORK(KWMAT),
     *                         WORK(KDIAGW),WORK(KFOCKD))
C
                  CALL T3_FORBIDDEN(WORK(KWMAT),ISYMY,ISYMB,B,ISYMD,D)
C
C                 -----------------------------
C                 Store WMAT as WBAR(ai,bj,dl):
C                 -----------------------------
                  DO ISYML = 1, NSYM
                   ISYMDL = MULD2H(ISYMD,ISYML)
                   ISAIBJ = MULD2H(ISYMY,ISYMDL)
                   DO L = 1, NRHF(ISYML)
                    DO ISYMJ = 1, NSYM
                     ISYMBJ = MULD2H(ISYMJ,ISYMB)
                     ISYMAI = MULD2H(ISAIBJ,ISYMBJ)
                     ISYAIL = MULD2H(ISYMAI,ISYML)
                     DO J = 1, NRHF(ISYMJ)

                      KOFF1 = ISAIKJ(ISYAIL,ISYMJ)+ NCKI(ISYAIL)*(J-1)
     *                      + ISAIK(ISYMAI, ISYML)+NT1AM(ISYMAI)*(L-1)

                      CALL DCOPY(NT1AMX,WORK(KWMAT+KOFF1),1,
     *                                  WBAR(1,1,B,J,D,L),1)
                     END DO
                    END DO
                   END DO
                  END DO

               ENDDO   ! B
            ENDDO      ! ISYMB
         ENDDO       ! D
      ENDDO          ! ISYMD 

*---------------------------------------------------------------------*
*     Close files and return:
*---------------------------------------------------------------------*
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
      CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')

      CALL QEXIT('CCSDT_GWBIC')
 
      RETURN
      END
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCSDT_GWTIC(LISTR,IDLSTR,WINT,THETA,T0AMP,
     &                       XLAMDP0,XLAMDH0,FOCK0,
     &                       LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     &                       WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: collect w^Y and theta^y intermediates in memory
*             for use in noddy code
*             
*    Written by Christof Haettig, Jan 2003, Friedrichstal
*            
*=====================================================================*
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccr1rsp.h>
#include <ccinftap.h>
#include <inftap.h>
C
      LOGICAL LOCDBG, TOT_T3Y
      PARAMETER(LOCDBG = .FALSE., TOT_T3Y = .TRUE.)
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
C
      CHARACTER*11 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR, FNDKBCR4
      PARAMETER(FN3SRTR  = 'CCSDT_ETA_1', FNCKJDR  = 'CCSDT_ETA_2', 
     *          FNDELDR  = 'CCSDT_ETA_3', FNDKBCR  = 'CCSDT_ETA_4',
     *          FNDKBCR4 = 'CCSDT_ETA_5')

      CHARACTER*(*) FNDKBC, FNDELD, FNCKJD
      
      CHARACTER LISTR*3

      CHARACTER MODEL*10, LABELY*8, CDUMMY*1
 
      INTEGER LWORK, IDLSTR, ISTAMY

      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR, LUDKBCR4,
     &        LUCKJD, LUDKBC, LUDELD, LUFOCK,
     &        KEND1, KFOCKY, KFOCKD, KFCKBA, KEND0, KT2TP, KT1B, KT2B,
     &        LWRK0, KTROC0, KTROC02,LWRK1,KINTOC, KEND2, LWRK2,
     &        KTRVI0, KTRVI1, KTRVI2, KEND3, LWRK3, KTRVI0Y, KINTVI,
     &        KSMAT, KSMAT3, KUMAT, KUMAT3, KDIAG, KDIAGW, KINDSQW,
     &        KINDSQ, KINDEX, KINDEX2, KTMAT, KTRVI8, KTRVI9, KTRVI10,
     &        KEND4, LWRK4, KTRVI8Y, KWMAT, KTROC0Y,
     &        IOPT, IRREP, IMAT, IERR, LENSQ, LENSQW, IOFF,
     &        ISYMC, ISYMK, KOFF1, KOFF2, ISYMD, ISYCKB, ISCKBY,
     &        ISYMB, ISYALJ, ISYALJ2, ISYMBD, ISCKIJ, ISYCKD, ISCKDY,
     &        ISWMAT, ISYML, ISYMDL, ISAIBJ, ISYMJ, ISYMBJ, ISYMAI,
     &        ISYAIL, KOFFA, KEND5, LWRK5
 
#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      REAL WORK(LWORK), FREQY, DDOT
      REAL ZERO, ONE, TWO, HALF
      REAL WINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
      REAL THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
      REAL T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      DOUBLE PRECISION WORK(LWORK), FREQY, DDOT
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      DOUBLE PRECISION WINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
      DOUBLE PRECISION THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
      DOUBLE PRECISION T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
#endif
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
      CALL QENTER('CCSDT_GWTIC')
 
      IF (LISTR(1:3).EQ.'R1 ') THEN
        ISTAMY = ISYLRT(IDLSTR)
        FREQY  = FRQLRT(IDLSTR)
        LABELY = LRTLBL(IDLSTR)
        IF (LORXLRT(IDLSTR)) 
     &   CALL QUIT('NO ORBITAL RELAXED PERTURBATION IN CCSDT_GWTIC')
      ELSE
        ! ups, probably higher-order response, not yet implemented here
        CALL QUIT('Unknown type of right vector in CCSDT_GWTIC.')
      END IF
C
C---------------------
C     Open some files: 
C---------------------
      IF (TOT_T3Y) THEN
         LU3SRTR  = 0
         LUCKJDR  = 0
         LUDELDR  = 0
         LUDKBCR  = 0
         LUDKBCR4 = 0
 
         CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
         CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
         CALL WOPEN2(LUDELDR,FNDELDR,64,0)
         CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
         CALL WOPEN2(LUDKBCR4,FNDKBCR4,64,0)
      ENDIF
C
C-----------------------------------------------
C     initial allocations and preparation of T2:
C-----------------------------------------------
      KEND0  = 1

      IF (LISTR(1:3).EQ.'R1 ') THEN
        KFOCKY = KEND0  
        KEND0  = KFOCKY + N2BST(ISTAMY)
      END IF

      KFOCKD = KEND0
      KFCKBA = KFOCKD + NORBTS
      KEND0  = KFCKBA + NT1AMX

      KT2TP  = KEND0
      KEND0  = KT2TP + NT2SQ(ISYM0)

      IF (TOT_T3Y) THEN
         KT1B   = KEND0
         KT2B   = KT1B + NT1AM(ISTAMY)
         KEND0  = KT2B + NT2SQ(ISTAMY)
      END IF

      LWRK0  = LWORK - KEND0
 
      IF (LWRK0.LT.NT2SQ(ISYM0) .OR. LWRK0.LT.NT2SQ(ISTAMY)) THEN
       CALL QUIT('Not enough memory to construct T2TP (CCSDT_GWTIC)')
      ENDIF
 

      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KT2TP))

      CALL CC_T2SQ(WORK(KT2TP),WORK(KEND0),ISYM0)

      CALL CC3_T2TP(WORK(KT2TP),WORK(KEND0),ISYM0)
 
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
     *    DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)

      IF (TOT_T3Y) THEN
         IOPT = 3
         CALL CC_RDRSP(LISTR,IDLSTR,ISTAMY,IOPT,MODEL,
     *                 WORK(KT1B),WORK(KT2B))
         CALL CCLR_DIASCL(WORK(KT2B),TWO,ISTAMY)
         CALL CC_T2SQ(WORK(KT2B),WORK(KEND0),ISTAMY)
         CALL CC3_T2TP(WORK(KT2B),WORK(KEND0),ISTAMY)
      END IF

C----------------------------------------------------
C     Calculate (ck|de)-tilde(Y) and (ck|lm)-tilde(Y)
C----------------------------------------------------
      IF (TOT_T3Y) THEN
         CALL CC3_BARINT(WORK(KT1B),ISTAMY,XLAMDP0,
     *                   XLAMDH0,WORK(KEND0),LWRK0,
     *                   LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
         CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISTAMY,LU3SRTR,FN3SRTR,
     *                  LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                  IDUMMY,CDUMMY)
C
         CALL CC3_SINT(XLAMDH0,WORK(KEND0),LWRK0,ISTAMY,
     *                 LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
         CALL CC3_TCME(XLAMDP0,ISTAMY,WORK(KEND0),LWRK0,
     *                 IDUMMY,CDUMMY,LUDKBCR,FNDKBCR,
     *                 IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                 IDUMMY,CDUMMY,LUDKBCR4,FNDKBCR4,2)
      ENDIF
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

      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
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

      IF (LOCDBG) THEN
         CALL AROUND('In CCSDT_GWTIC: Triples Fock MO matrix (sort)')
         CALL CC_PRFCKMO(WORK(KFCKBA),ISYM0)
      ENDIF
C
C-----------------------------------------------------------
C     Prepare one-electron operators needed to compute the
C     amplitude response vectors:
C-----------------------------------------------------------
C
      IF (LISTR(1:3).EQ.'R1 ') THEN
        ! read property integrals from file:
        CALL CCPRPAO(LABELY,WORK(KFOCKY),IRREP,IMAT,IERR,
     *               WORK(KEND0),LWRK0)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISTAMY)) THEN
          WRITE(LUPRI,*) 'ISTAMY:',ISTAMY
          WRITE(LUPRI,*) 'IRREP :',IRREP
          WRITE(LUPRI,*) 'IERR  :',IERR
          CALL QUIT('CCSDT_GWTIC: error reading operator '//LABELY)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(WORK(KFOCKY),N2BST(ISTAMY))
        END IF
 
        ! transform property integrals to Lambda-MO basis
        CALL CC_FCKMO(WORK(KFOCKY),XLAMDP0,XLAMDH0,
     &                WORK(KEND0),LWRK0,ISTAMY,1,1)
      END IF

C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
      KTROC0  = KEND0
      KTROC02 = KTROC0  + NTRAOC(ISYM0)
      KEND1   = KTROC02 + NTRAOC(ISYM0)
      IF (TOT_T3Y) THEN
         KTROC0Y = KEND1
         KEND1   = KTROC0Y + NTRAOC(ISTAMY)
      ENDIF
      LWRK1   = LWORK   - KEND1

      KINTOC  = KEND1
      KEND2   = KINTOC + MAX(NTOTOC(ISYM0),NTOTOC(ISYM0))
      LWRK2   = LWORK  - KEND2

      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CCSDT_GWTIC')
      END IF
C
C------------------------
C     Occupied integrals.
C------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYM0) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF

      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of OCC-INT ',
     *    DDOT(NTOTOC(ISYM0),WORK(KINTOC),1,WORK(KINTOC),1)
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),XLAMDP0,
     *               WORK(KEND2),LWRK2,ISYM0)

      CALL CCFOP_SORT(WORK(KTROC0),WORK(KTROC02),ISYM0,1)
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of TROC0 ',
     *     DDOT(NTRAOC(ISYM0),WORK(KTROC0),1,WORK(KTROC0),1)
C
C------------------------------------------
C     Y transformed Occupied integrals.
C-----------------------------------------
C
      IF (TOT_T3Y) THEN
         IOFF = 1
         IF (NTOTOC(ISTAMY) .GT. 0) THEN
            CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,
     *                  NTOTOC(ISTAMY))
         ENDIF
C
         IF (LOCDBG)
     *      WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (Y transformed) ',
     *        DDOT(NTOTOC(ISTAMY),WORK(KINTOC),1,WORK(KINTOC),1)
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0Y),XLAMDP0,
     *                  WORK(KEND2),LWRK2,ISTAMY)
      ENDIF
 
C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKB = MULD2H(ISYMD,ISYM0)
         ISCKBY = MULD2H(ISTAMY,ISYMD)
         IF (LOCDBG) WRITE(LUPRI,*) 'In CCSDT_ETA_CON: ISYCKB :',ISYCKB

         DO D = 1,NVIR(ISYMD)
C
C           ------------------
C           Memory allocation.
C           ------------------
            KTRVI0  = KEND1
            KTRVI1  = KTRVI0 + NCKATR(ISYCKB)
            KTRVI2  = KTRVI1 + NCKATR(ISYCKB)
            KEND3   = KTRVI2 + NCKATR(ISYCKB)
            LWRK3   = LWORK  - KEND3
           
            IF (TOT_T3Y) THEN
               KTRVI0Y  = KEND3
               KEND3    = KTRVI0Y  + NCKATR(ISCKBY)
               LWRK3    = LWORK    - KEND3
            ENDIF
           
            KINTVI = KEND3
            KEND4  = KINTVI + MAX(NCKA(ISYMD),NCKA(ISYCKB))
            LWRK4  = LWORK  - KEND4
           
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CCSDT_GWTIC')
            END IF
C
C           ------------------------------------
C           Integrals used in s3am.
C           ------------------------------------
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                     NCKATR(ISYCKB))
            ENDIF
C
C           ------------------------------------------------------------
C           Read B transformed virtual integrals used for W for TOT_T3Y
C           ------------------------------------------------------------
            IF (TOT_T3Y) THEN
               IOFF = ICKBD(ISCKBY,ISYMD) + NCKATR(ISCKBY)*(D - 1) + 1
               IF (NCKATR(ISCKBY) .GT. 0) THEN
                  CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI0Y),IOFF,
     &                        NCKATR(ISCKBY))
               ENDIF
            ENDIF
 
C           ------------------------------------------------
C           Sort the integrals for s3am:
C           ------------------------------------------------
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISYM0)

C           -------------------------------------------
C           Read virtual integrals used in contraction.
C           -------------------------------------------
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMDH0,
     *                       ISYMD,D,ISYM0,WORK(KEND4),LWRK4)
 
 
            DO ISYMB = 1,NSYM

               ISYALJ  = MULD2H(ISYMB,ISYM0)
               ISYALJ2 = MULD2H(ISYMD,ISYM0)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYM0)
               ISYCKD  = MULD2H(ISYM0,ISYMB)
               ISCKDY  = MULD2H(ISTAMY,ISYMB)
               ISWMAT  = MULD2H(ISCKIJ,ISTAMY)

               IF (LOCDBG) THEN
                  WRITE(LUPRI,*) 'In CCSDT_GWTIC: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CCSDT_GWTIC: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CCSDT_GWTIC: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CCSDT_GWTIC: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CCSDT_GWTIC: ISCKIJ:',ISCKIJ
               ENDIF
C
C              Can use kend3 since we do not need the integrals anymore.
               KSMAT   = KEND3
               KSMAT3  = KSMAT   + NCKIJ(ISCKIJ)
               KUMAT   = KSMAT3  + NCKIJ(ISCKIJ)
               KUMAT3  = KUMAT   + NCKIJ(ISCKIJ)
               KDIAG   = KUMAT3  + NCKIJ(ISCKIJ)
               KDIAGW  = KDIAG   + NCKIJ(ISCKIJ)
               KINDSQW = KDIAGW  + NCKIJ(ISWMAT)
               KINDSQ  = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ)  - 1)/IRAT + 1
               KTMAT   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KTRVI8  = KTMAT   + NCKIJ(ISCKIJ)
               KTRVI9  = KTRVI8  + NCKATR(ISYCKD)
               KTRVI10 = KTRVI9  + NCKATR(ISYCKD)
               KEND4   = KTRVI10 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4

               IF (TOT_T3Y) THEN
                  KTRVI8Y = KEND4
                  KEND4   = KTRVI8Y + NCKATR(ISCKDY)
               ENDIF

               KWMAT   = KEND4
               KEND4   = KWMAT   + NCKIJ(ISWMAT)
               LWRK4   = LWORK   - KEND4

               KINTVI  = KEND4
               KEND5   = KINTVI  + MAX(NCKA(ISYMB),NCKA(ISYCKD))
               LWRK5   = LWORK   - KEND5

               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CCSDT_GWTIC')
               END IF
C
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
               CALL CC3_DIAG(WORK(KDIAG), WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)

               IF (LOCDBG) THEN
                 WRITE(LUPRI,*) 'Norm of DIA  ',
     *              DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,WORK(KDIAG),1)
                 WRITE(LUPRI,*) 'Norm of DIA_W',
     *              DDOT(NCKIJ(ISWMAT),WORK(KDIAGW),1,WORK(KDIAGW),1)
               END IF
C
C              -----------------------
C              Construct index arrays.
C              -----------------------
               LENSQ  = NCKIJ(ISCKIJ)
               LENSQW = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               CALL CC3_INDEX(WORK(KINDEX),ISYALJ)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJ2)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT) 

               DO B = 1,NVIR(ISYMB)
                  IF (LOCDBG) write(lupri,*) 'For B, D : ',B,D
 
C                 ---------------------------------------------
C                 Initialize WMAT:
C                 ---------------------------------------------
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))

C                 ---------------------------------------------
C                 Read and transform integrals used in second S
C                 ---------------------------------------------
                  IOFF = ICKBD(ISYCKD,ISYMB) 
     *                 + NCKATR(ISYCKD)*(B - 1) + 1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI8),IOFF,
     *                           NCKATR(ISYCKD))
                  ENDIF

                  CALL CCSDT_SRTVIR(WORK(KTRVI8),WORK(KTRVI9),
     *                              WORK(KEND4),LWRK4,ISYMB,ISYM0)
C
C                 ----------------------------------------------
C                 Read B transformed virtual integrals used in W
C                 ----------------------------------------------
                  IF (TOT_T3Y) THEN
                     IOFF = ICKBD(ISCKDY,ISYMB) + 
     &                      NCKATR(ISCKDY)*(B - 1) + 1
                     IF (NCKATR(ISCKDY) .GT. 0) THEN
                        CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI8Y),IOFF,
     &                              NCKATR(ISCKDY))
                     ENDIF
                  ENDIF
C
C                 -----------------------------------------
C                 Read virtual integrals used in second U
C                 -----------------------------------------
                  IOFF = ICKAD(ISYCKD,ISYMB) 
     *                 + NCKA(ISYCKD)*(B - 1) + 1
                  IF (NCKA(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                           NCKA(ISYCKD))
                  ENDIF

                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI10),
     *                             XLAMDH0,ISYMB,B,ISYM0,
     *                             WORK(KEND5),LWRK5)
C
C                 --------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for B,D.
C                 --------------------------------------------------
                  CALL CC3_SMAT(0.0D0,WORK(KT2TP),ISYM0,WORK(KTMAT),
     *                          WORK(KTRVI0),WORK(KTRVI2),
     *                          WORK(KTROC0),ISYM0,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMAT),WORK(KEND4),LWRK4,
     *                          WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)

                  CALL T3_FORBIDDEN(WORK(KSMAT),ISYM0,ISYMB,B,ISYMD,D)
C
                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of SMAT :',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,WORK(KSMAT),1)
C
C                 --------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for D,B.
C                 --------------------------------------------------
                  CALL CC3_SMAT(0.0D0,WORK(KT2TP),ISYM0,WORK(KTMAT),
     *                          WORK(KTRVI8),WORK(KTRVI9),
     *                          WORK(KTROC0),ISYM0,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMAT3),WORK(KEND4),LWRK4,
     *                          WORK(KINDEX2),WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)

                  CALL T3_FORBIDDEN(WORK(KSMAT3),ISYM0,ISYMD,D,ISYMB,B)
C
                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of SMAT3:',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KSMAT3),1,WORK(KSMAT3),1)
C
C                 ---------------------------------
C                 Calculate U(ci,jk) for fixed b,d.
C                 ---------------------------------
                  CALL CC3_UMAT(0.0D0,WORK(KT2TP),ISYM0,WORK(KTRVI1),
     *                          WORK(KTROC02),ISYM0,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMAT),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)

                  CALL T3_FORBIDDEN(WORK(KUMAT),ISYM0,ISYMB,B,ISYMD,D)
C
                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of UMAT :',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KUMAT),1,WORK(KUMAT),1)
C
C                 ---------------------------------
C                 Calculate U(ci,jk) for fixed d,b.
C                 ---------------------------------
                  CALL CC3_UMAT(0.0D0,WORK(KT2TP),ISYM0,WORK(KTRVI10),
     *                          WORK(KTROC02),ISYM0,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMAT3),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)

                  CALL T3_FORBIDDEN(WORK(KUMAT3),ISYM0,ISYMD,D,ISYMB,B)

                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of UMAT3:',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KUMAT3),1,WORK(KUMAT3),1)
C
C                 -----------------------------------
C                 Sum up the S and U to get a real T3
C                 -----------------------------------
                  CALL CC3_T3BD(ISCKIJ,WORK(KSMAT),WORK(KSMAT3),
     *                                 WORK(KUMAT),WORK(KUMAT3),
     *                                 WORK(KTMAT),WORK(KINDSQ),LENSQ)
C
C                 ---------------------------------
C                 Use the T3 in TMAT to calculate W
C                 ---------------------------------
                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of TMAT:',
     *               DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,WORK(KTMAT),1)

                  CALL WBD_V(WORK(KTMAT),ISCKIJ,WORK(KFOCKY),ISTAMY,
     *                            WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)

                   IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_V):',
     *                DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)
                  
                  IF (.TRUE.) THEN
C                  ---------------------------------------------------
C                  1)Store this contribution of WMAT in THETA(ai,bj,dl):
C                    but first divide by orbital energies and clean up
C                  2)Store zeroth-order T3 amplitudes in T0AMP:
C                  ---------------------------------------------------
                   CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQY, 
     *                          ISWMAT,WORK(KWMAT), 
     *                          WORK(KDIAGW),WORK(KFOCKD)) 
                   CALL T3_FORBIDDEN(WORK(KWMAT),ISTAMY,ISYMB,B,ISYMD,D)
                  
                   DO ISYML = 1, NSYM
                    ISYMDL = MULD2H(ISYMD,ISYML)
                    ISAIBJ = MULD2H(ISTAMY,ISYMDL)
                    DO L = 1, NRHF(ISYML)
                     DO ISYMJ = 1, NSYM
                      ISYMBJ = MULD2H(ISYMJ,ISYMB)
                      ISYMAI = MULD2H(ISAIBJ,ISYMBJ)
                      ISYAIL = MULD2H(ISYMAI,ISYML)
                      DO J = 1, NRHF(ISYMJ)
                  
                       KOFFA = ISAIKJ(ISYAIL,ISYMJ)+ NCKI(ISYAIL)*(J-1)
     *                       + ISAIK(ISYMAI, ISYML)+NT1AM(ISYMAI)*(L-1)
                  
                       CALL DCOPY(NT1AMX,WORK(KWMAT+KOFFA),1,
     &                                   THETA(1,1,B,J,D,L),1)

                       CALL DCOPY(NT1AMX,WORK(KTMAT+KOFFA),1,
     &                                   T0AMP(1,1,B,J,D,L),1)
                  
                      END DO
                     END DO
                    END DO
                   END DO
                   ! Reinitialize WMAT to avoid a double count:
                   CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
                  END IF

                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_V)',
     *               DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)

                  CALL WBD_O(WORK(KTMAT),ISCKIJ,WORK(KFOCKY),ISTAMY,
     *                            WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)

                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_O)',
     *               DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)
                  CALL WBD_T2(B,ISYMB,D,ISYMD,WORK(KT2TP),ISYM0,
     *                        WORK(KFOCKY),ISTAMY,WORK(KINDSQW),LENSQW,
     *                        WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)

                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_T2):',
     *               DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)

                  IF (TOT_T3Y) THEN
                     CALL WBD_GROUND(WORK(KT2B),ISTAMY,WORK(KTMAT),
     *                               WORK(KTRVI0),WORK(KTRVI8),
     *                               WORK(KTROC0),1,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)
c
                     CALL WBD_GROUND(WORK(KT2TP),ISYM0,WORK(KTMAT),
     *                               WORK(KTRVI0Y),WORK(KTRVI8Y),
     *                               WORK(KTROC0Y),ISTAMY,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)
                  ENDIF

                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQY, 
     *                         ISWMAT,WORK(KWMAT), 
     *                         WORK(KDIAGW),WORK(KFOCKD)) 

                  CALL T3_FORBIDDEN(WORK(KWMAT),ISTAMY,ISYMB,B,ISYMD,D) 

                  IF (LOCDBG) WRITE(LUPRI,*) 'Norm of WMAT (WBD_DIA)',
     *               DDOT(NCKIJ(ISWMAT),WORK(KWMAT),1,WORK(KWMAT),1)
 
C                 -----------------------------
C                 Store WMAT in WINT(ai,bj,dl):
C                 -----------------------------
                  DO ISYML = 1, NSYM
                   ISYMDL = MULD2H(ISYMD,ISYML)
                   ISAIBJ = MULD2H(ISTAMY,ISYMDL)
                   DO L = 1, NRHF(ISYML)
                    DO ISYMJ = 1, NSYM
                     ISYMBJ = MULD2H(ISYMJ,ISYMB)
                     ISYMAI = MULD2H(ISAIBJ,ISYMBJ)
                     ISYAIL = MULD2H(ISYMAI,ISYML)
                     DO J = 1, NRHF(ISYMJ)
    
                      KOFFA = ISAIKJ(ISYAIL,ISYMJ)+ NCKI(ISYAIL)*(J-1)
     *                      + ISAIK(ISYMAI, ISYML)+NT1AM(ISYMAI)*(L-1)

                      CALL DCOPY(NT1AMX,WORK(KWMAT+KOFFA),1,
     &                                  WINT(1,1,B,J,D,L),1)

                     END DO
                    END DO
                   END DO
                  END DO

               ENDDO   ! B
            ENDDO      ! ISYMB

         ENDDO       ! D
      ENDDO          ! ISYMD
c
 
C-------------
C     End
C-------------
      IF (TOT_T3Y) THEN
         CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
         CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
         CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
         CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
         CALL WCLOSE2(LUDKBCR4,FNDKBCR4,'DELETE')
      ENDIF

      CALL QEXIT('CCSDT_GWTIC')
 
      RETURN
      END
*=====================================================================*
