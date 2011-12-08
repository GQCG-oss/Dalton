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
      SUBROUTINE CCSDT_ADEN_NODDY(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,DO_DIA,DIA,DO_XMMAT,XMMAT,
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
#include <ccr2rsp.h>
#include <ccr1rsp.h>
#include <ccer1rsp.h>
#include <ccnoddy.h>

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL LOCDBG, HAVETB3T3, PRINTDENS
      PARAMETER (LOCDBG=.FALSE., HAVETB3T3=.TRUE., PRINTDENS=.FALSE.)

      CHARACTER*3 LISTL, LISTR
      CHARACTER*(*) FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3
      CHARACTER*(*) FN3FOPX, FN3FOP2X
      LOGICAL DO_DIA, DO_XMMAT
      INTEGER LUDELD,LUCKJD,LUDKBC,LUTOC,LU3VI,LUDKBC3,LU3FOPX,LU3FOP2X
      INTEGER LWORK, IDLSTL, IDLSTR

#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      REAL DIJ(*), DAB(*), DIA(*), XMMAT(*)
      REAL WORK(LWORK)
      REAL HALF, ONE, ZERO, TWO, THREE, DDOT, XNORM
      REAL FREQLST, FREQY, FREQX, FREQR
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      DOUBLE PRECISION DIJ(*), DAB(*), DIA(*), XMMAT(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION HALF, ONE, TWO, ZERO, THREE, DDOT, XNORM
      DOUBLE PRECISION FREQLST, FREQY, FREQX, FREQR, FREQL
#endif
      PARAMETER( ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TWO = 2.0D0,
     &           THREE=3.0D0 )

      CHARACTER MODEL*10, LISTX*3, LISTY*3, LABELX*8, LABELY*8
      LOGICAL CUBIC, LORXX, LORXY, HAVEX, HAVEY, HAVEXY, 
     &        SPLIT_WBTB, SPLIT_WT
      INTEGER KEND1, LWRK1, KWBAR, KWYINT, KTHETAY, KTBAR3, KTAMP3,
     &        NCK, NBJ, NAI, KAIBJCK, KCKAIBJ, KBJCKAI, IOFF, ISYMR,
     &        KWTILDE, KWBREVE, KTHEOCC, KTVIRT1, KTVIRT2, KEND2, 
     &        LWRK2, KTHVIRT, KSWAP, KDABR, KDIJR, KDIAR, IOPT,
     &        KT0AMP, KTAMP2, KT0FMI, KTYFMI, KT0MI, KTYMI, 
     &        KWDE, KWED, KSCR, KFNEMDL, KEMDLFN, KT3AMP0, KDLEI,
     &        NAN, NFN, NFI, NEM, NDL, KDLEMFN, KDLAN, KEMFI, NEI,
     &        KL2AMP, KDLFIAN, KDLFNAI, KDLFN, ISYML,
     &        KDAB0, KDAB0R, KDIJ0, KDIJ0R, KTAMP1, KXMMAT, KIMFN,
     &        KTHBAR, KWXINT, KTHETAX, KWXYVINT, KWXYOINT, KXPROP,
     &        KYPROP, KXYMAT, KFOCKD, KTHEOV, KTVIRT3, KSCRVVVO,
     &        ISYMX, ISYMY, IDLSTX, IDLSTY, LUSIFC, LUTEMP,
     &        KSCR1, KINT1T0, KINT2T0, KXIAJB, KYIAJB, 
     &        KLAMPB, KLAMHB, KFOCKB, 
     &        KINT1S0, KINT2S0, KINT1SB, KINT2SB, KDUM

      INTEGER ILSTSYM, IR1TAMP

*---------------------------------------------------------------------*
*     Begin:
*---------------------------------------------------------------------*
      CALL QENTER('CCTADNOD')

      IF (NSYM.GT.1) CALL QUIT('No symmetry in CCSDT_ADEN_NODDY!')

      IF   (LISTR(1:3).EQ.'R1 ') THEN
        ! we compute the first-order amplitude response vectors
        CUBIC  = .FALSE.
        LORXX  = LORXLRT(IDLSTR)
        LORXY  = .FALSE.

      ELSE IF (LISTR(1:3).EQ.'RE ') THEN
        ! we compute eigenvectors
        CUBIC  = .FALSE.
        LORXX  = .FALSE.
        LORXY  = .FALSE.
      ELSE IF (LISTR(1:3).EQ.'R2 ') THEN
        ! we compute the second-order amplitude response vectors
        CUBIC  = .TRUE.

        LISTX  = 'R1 '
        LABELX = LBLR2T(IDLSTR,1)
        ISYMX  = ISYR2T(IDLSTR,1)
        FREQX  = FRQR2T(IDLSTR,1)
        LORXX  = LORXR2T(IDLSTR,1)
        IDLSTX = IR1TAMP(LABELX,LORXX,FREQX,ISYMX)

        LISTY  = 'R1 '
        LABELY = LBLR2T(IDLSTR,2)
        ISYMY  = ISYR2T(IDLSTR,2)
        FREQY  = FRQR2T(IDLSTR,2)
        LORXY  = LORXR2T(IDLSTR,2)
        IDLSTY = IR1TAMP(LABELY,LORXY,FREQY,ISYMY)

      ELSE IF (LISTR(1:3).EQ.'ER1') THEN
        ! we compute first-order response of excited states
       IF (.FALSE.) THEN
        ! use quadratic response code
        ! (for internal check of noddy code)
        CUBIC  = .FALSE.
        LORXX  = LORXER1(IDLSTR)
        LORXY  = .FALSE.
       ELSE
        ! use cubic response code
        ! (this will be most similar to the real code)
        CUBIC  = .TRUE.

        LISTX  = 'R1 '
        LABELX = LBLER1(IDLSTR)
        ISYMX  = ISYOER1(IDLSTR)
        FREQX  = FRQER1(IDLSTR)
        LORXX  = LORXER1(IDLSTR)
        IDLSTX = IR1TAMP(LABELX,LORXX,FREQX,ISYMX)

        LISTY  = 'RE '
        ISYMY  = ISYSER1(IDLSTR)
        FREQY  = EIGER1(IDLSTR)
        LORXY  = .FALSE.
        IDLSTY = ISTER1(IDLSTR)
       END IF
      ELSE
        PRINT *,'LISTR:',LISTR
        CALL QUIT('Unknown LISTR in CCSDT_ADEN_NODDY.')
      END IF
 
      IF (LORXX.OR.LORXY) THEN
       CALL QUIT('Orbital relaxation not allowed in CCSDT_ADEN_NODDY.')
      END IF

      HAVEX  = CUBIC .AND. (LISTX(1:3).EQ.'R1 ')
      HAVEY  = CUBIC .AND. (LISTY(1:3).EQ.'R1 ')
      HAVEXY = HAVEX .OR. HAVEY

      SPLIT_WBTB = CUBIC
      SPLIT_WT   = (LISTR(1:3).NE.'RE ')

      ISYML  = ILSTSYM(LISTL,IDLSTL)
      ISYMR  = ILSTSYM(LISTR,IDLSTR)
      FREQR  = FREQLST(LISTR,IDLSTR)

*---------------------------------------------------------------------*
*     Allocate arrays kept until the end of the routine:
*---------------------------------------------------------------------*
      KEND1 = 1

      KDAB0 = KEND1
      KDIJ0 = KDAB0  + NVIRT*NVIRT
      KEND1 = KDIJ0  + NRHFT*NRHFT

      KFOCKD = KEND1
      KEND1  = KFOCKD + NORBT

      KWBAR = KEND1
      KEND1 = KWBAR + NT1AMX*NT1AMX*NT1AMX

      IF (SPLIT_WBTB) THEN
        KTHBAR = KEND1
        KEND1  = KTHBAR + NT1AMX*NT1AMX*NT1AMX
      END IF

      KWYINT  = KEND1
      KEND1   = KWYINT  + NT1AMX*NT1AMX*NT1AMX

      IF (SPLIT_WT) THEN
        KTHETAY = KWYINT  + NT1AMX*NT1AMX*NT1AMX
        KEND1   = KTHETAY + NT1AMX*NT1AMX*NT1AMX
      END IF

      IF (CUBIC) THEN
        KWXINT  = KEND1
        KTHETAX = KWXINT  + NT1AMX*NT1AMX*NT1AMX
        KEND1   = KTHETAX + NT1AMX*NT1AMX*NT1AMX
        
        KWXYVINT = KEND1
        KWXYOINT = KWXYVINT + NT1AMX*NT1AMX*NT1AMX
        KEND1    = KWXYOINT + NT1AMX*NT1AMX*NT1AMX

        KXPROP = KEND1
        KYPROP = KXPROP + NBAST*NBAST
        KXYMAT = KYPROP + NBAST*NBAST
        KEND1  = KXYMAT + NBAST*NBAST
      ELSE
        KWXYVINT = KWYINT
        KWXYOINT = KWYINT
      END IF

      KT3AMP0 = KEND1 
      KEND1   = KT3AMP0 + NT1AMX*NT1AMX*NT1AMX
      
      KL2AMP  = KEND1 
      KEND1   = KL2AMP + NT1AMX*NT1AMX

      KTAMP1 = KEND1
      KEND1  = KTAMP1 + NT1AMX

      KT0AMP  = KEND1
      KTAMP2  = KT0AMP + NT1AMX*NT1AMX
      KEND1   = KTAMP2 + NT1AMX*NT1AMX

      IF (HAVETB3T3) THEN
        KTBAR3 = KEND1
        KTAMP3 = KTBAR3 + NT1AMX*NT1AMX*NT1AMX
        KDABR  = KTAMP3 + NT1AMX*NT1AMX*NT1AMX
        KDIJR  = KDABR  + NVIRT*NVIRT
        KDIAR  = KDIJR  + NRHFT*NRHFT
        KDAB0R = KDIAR  + NVIRT*NRHFT
        KDIJ0R = KDAB0R + NVIRT*NVIRT
        KEND1  = KDIJ0R + NRHFT*NRHFT

        IF (DO_XMMAT) THEN
          KXMMAT = KEND1
          KEND1  = KXMMAT + NRHFT*NRHFT*NT1AMX
        END IF
      END IF

      LWRK1 = LWORK - KEND1
      IF (LWRK1.LT.0)CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (1)')

      IF (NSYM.GT.1) CALL QUIT('No symmetry in CCSDT_ADEN_NODDY!')

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD-1+I), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     read doubles amplitudes (zeroth-order and response) from file:
*---------------------------------------------------------------------*
      IF (LWRK1.LT.NT2AMX)
     &   CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (3a)')

      ! read T^0 doubles amplitudes from file and square up
      IOPT   = 2
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT0AMP),ISYM0)

      ISYMR = ILSTSYM(LISTR,IDLSTR)

      ! read T^X doubles amplitudes from file and square up
      IOPT = 3
      Call CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,
     &              WORK(KTAMP1),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMR)  
      CALL CC_T2SQ(WORK(KEND1),WORK(KTAMP2),ISYMR)

      ISYML = ILSTSYM(LISTL,IDLSTL)

      ! read left doubles amplitudes from file and square up
      IOPT = 2
      Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,DUMMY,WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KL2AMP),ISYML)

*---------------------------------------------------------------------*
*     Get W-bar (and Theta-bar) intermediate into memory:
*---------------------------------------------------------------------*
      IF (LISTL.EQ.'L0 ') THEN

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODL30,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (WORK(KWBAR+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP') 

        IF (HAVETB3T3) THEN
         CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,WORK(KWBAR),1,WORK(KTBAR3),1)
        END IF

        ! cheat the code by putting 1/3 Tbar^(0) into W-bar
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,ONE/THREE,WORK(KWBAR),1)

        IF (SPLIT_WBTB) THEN
          ! Theta-bar is set to zero for this case
          CALL DZERO(WORK(KTHBAR),NT1AMX*NT1AMX*NT1AMX)
        END IF

      ELSE IF (LISTL.EQ.'L1 ') THEN

        CALL CCSDT_GWBIC(LISTL,IDLSTL,WORK(KWBAR),WORK(KTHBAR),
     &                   SPLIT_WBTB,XLAMDP0,XLAMDH0,FOCK0,
     &                   LUTOC,FNTOC,LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     &                   LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     &                   WORK(KEND1),LWRK1)

        ! turn sign of Wbar intermediate
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KWBAR),1)

        IF (SPLIT_WBTB) THEN
          ! turn also sign of Theta-bar intermediate
          CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KTHBAR),1)

          ! add THBAR intermediate to WBAR to obtain the full WBAR
          CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,ONE,WORK(KTHBAR),1,
     &                                      WORK(KWBAR),1)
        END IF

        IF (LOCDBG) THEN
          XNORM = DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KWBAR),1,WORK(KWBAR),1)
          WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> Wbar^Y:',XNORM
          WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTL,IDLSTL:',LISTL,IDLSTL
          CALL PRINT_PT3_NODDY(WORK(KWBAR))
        END IF

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
        END IF

      ELSEIF (LISTL(1:3).EQ.'LE '.OR.LISTL(1:3).EQ.'M1 '.OR.
     &        LISTL(1:3).EQ.'E0 '.OR.LISTL(1:3).EQ.'N2 '    ) THEN

        KSCR1 = KEND1
        KEND2 = KSCR1 + NT1AMX

        KINT1T0 = KEND2
        KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
        KEND2   = KINT2T0 + NRHFT*NRHFT*NT1AMX
   
        KXIAJB  = KEND2
        KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
        KEND2   = KYIAJB  + NT1AMX*NT1AMX

        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LT. 0)
     &      CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (2)')

C       ---------------------------------------
C       Read precalculated integrals from file:
C             XINT1T0 =  (KC|BD)
C             XINT2T0 =  (KC|LJ)
C             XIAJB   = 2(IA|JB) - (IB|JA)
C             YIAJB   =  (IA|JB)
C       ---------------------------------------
        CALL CCSDT_READ_NODDY(.FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        .TRUE.,WORK(KXIAJB),WORK(KYIAJB),
     &                        .FALSE.,DUMMY,DUMMY,
     &                        .TRUE.,WORK(KINT1T0),WORK(KINT2T0),
     &                        .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

        CALL DZERO(WORK(KWBAR),NT1AMX*NT1AMX*NT1AMX)
        CALL CCSDT_TBAR31_NODDY(WORK(KWBAR),FREQL,LISTL,IDLSTL,
     &                          XLAMDP0,XLAMDH0,FOCK0,
     &                          WORK(KFOCKD),WORK(KSCR1),
     &                          WORK(KXIAJB),WORK(KINT1T0),
     &                          WORK(KINT2T0),WORK(KEND2),LWRK2)

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KWBAR),1)

        IF (HAVETB3T3) THEN
         CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,WORK(KWBAR),1,WORK(KTBAR3),1)
        END IF

        IF (SPLIT_WBTB) THEN
          ! Theta-bar is set to zero for this case
          CALL DZERO(WORK(KTHBAR),NT1AMX*NT1AMX*NT1AMX)
        END IF

        ! account for the fact that in the real code a combination
        ! of three W intermediates give the complete Tbar3
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,ONE/THREE,WORK(KWBAR),1)

      ELSE
        CALL QUIT('Unknown or illegal list in CCDOTRSP_NODDY.')
      END IF

      IF (LOCDBG .AND. HAVETB3T3) THEN
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> Tbar^Y:'
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTL,IDLSTL:',LISTL,IDLSTL
        WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> norm^2:',
     &     DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KTBAR3),1,WORK(KTBAR3),1)
        CALL PRINT_PT3_NODDY(WORK(KTBAR3))
      END IF

*---------------------------------------------------------------------*
*     Compute intermediates needed to represent the first-
*     or second-order response amplitudes:
*---------------------------------------------------------------------*
      IF (.NOT. CUBIC) THEN

       IF (LISTR.EQ.'R1 ') THEN
        ! Get w^y and theta^y intermediate into memory:
        ! (returns in addition the zeroth-order amplitudes in T3AMP0,
        !  and for double checking T3^Y on TAMP3)
        CALL CCSDT_GWTIC(LISTR,IDLSTR,WORK(KWYINT),WORK(KTHETAY),
     &                   WORK(KT3AMP0),HAVETB3T3,WORK(KTAMP3),LOCDBG,
     &                   XLAMDP0,XLAMDH0,FOCK0,
     &                   LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     &                   WORK(KEND1),LWRK1)

        ! turn sign of t3, w and theta intermediates
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KT3AMP0),1)
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KWYINT),1)
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KTHETAY),1)
        IF (HAVETB3T3) 
     &      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KTAMP3),1)

       ELSE

        KEND2   = KEND1

        KINT1S0 = KEND2
        KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
        KEND2   = KINT2S0 + NRHFT*NRHFT*NT1AMX
   
        IF (LISTR.EQ.'RE ') THEN
          KLAMPB  = KEND2
          KLAMHB  = KLAMPB + NLAMDT
          KFOCKB  = KLAMHB + NLAMDT
          KEND2   = KFOCKB + N2BST(1)

          KINT1SB = KEND2
          KINT2SB = KINT1SB + NT1AMX*NVIRT*NVIRT
          KEND2   = KINT2SB + NRHFT*NRHFT*NT1AMX
        ELSE IF (LISTR.EQ.'R2 '.OR. LISTR.EQ.'ER1') THEN
          KSCR1 = KEND2
          KEND2 = KSCR1 + NT1AMX
        END IF


        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LT. 0)
     &      CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (3)')

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


        CALL DZERO(WORK(KWYINT),NT1AMX*NT1AMX*NT1AMX)

        IF (LISTR.EQ.'RE '.OR.LISTR.EQ.'RC ') THEN
          KDUM = KEND2
          CALL CCSDT_T31_NODDY(WORK(KWYINT),LISTR,IDLSTR,FREQR,.FALSE.,
     &                         .FALSE.,WORK(KINT1S0),WORK(KINT2S0),
     &                         .FALSE.,DUMMY,DUMMY,
     &                         .FALSE.,DUMMY,DUMMY,
     &                                 WORK(KINT1SB),WORK(KINT2SB),
     &                         WORK(KLAMPB),WORK(KLAMHB),WORK(KFOCKB),
     &                         XLAMDP0,XLAMDH0,FOCK0,
     &                         WORK(KDUM),WORK(KFOCKD),
     &                         WORK(KEND2),LWRK2) 
        ELSE IF (LISTR.EQ.'R2 '.OR. LISTR.EQ.'ER1') THEN
          CALL CCSDT_T32_NODDY(WORK(KWYINT),LISTR,IDLSTR,FREQR,
     &                         WORK(KINT1S0),WORK(KINT2S0),
     &                         XLAMDP0,XLAMDH0,FOCK0,
     &                         WORK(KFOCKD),DUMMY,DUMMY,
     &                         WORK(KSCR1),WORK(KEND2),LWRK2)
        ELSE
          CALL QUIT('Unknown or illegal LISTR in CCSDT_ADEN_NODDY.')
        END IF

        IF (HAVETB3T3) THEN
         CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,WORK(KWYINT),1,WORK(KTAMP3),1)
        END IF

        IF (SPLIT_WT) THEN
          ! Theta-Y is set to zero for this case
          CALL DZERO(WORK(KTHETAY),NT1AMX*NT1AMX*NT1AMX)
        END IF

        ! account for the fact that in the real code a combination
        ! of three W intermediates give the complete Tbar3
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,ONE/THREE,WORK(KWYINT),1)

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (WORK(KT3AMP0+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP') 

        IF (LOCDBG .AND. HAVETB3T3) THEN
         WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> T^Y:'
         WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> LISTR,IDLSTR:',LISTR,IDLSTR
         WRITE(LUPRI,*) 'CCSDT_ADEN_NODDY> norm^2:',
     &     DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KTAMP3),1,WORK(KTAMP3),1)
         CALL PRINT_PT3_NODDY(WORK(KTAMP3))
        END IF

       END IF

      ELSE ! IF (CUBIC) THEN

        ! Get for both perturbations the intermediates:
        !  w^y and theta^y with same meaning as above
        !  YPROP contains matrix elements of y
        ! In addition some second-order intermediates are set up:
        !  XYMAT contains matrix elements of [X,T1^Y] + [Y,T1^X]
        !  WXYVINT contains theta(e-bar m-bar, f n-bar, d l-bar)
        !  WXYOINT contains theta(e m-bar-bar, f n-bar-bar, d l-bar-bar)
        ! (returns also the zeroth-order amplitudes in T3AMP0,
        !  and for double checking the complete T3^XY on TAMP3)
        CALL CCSDT_T32INT_NODDY(LISTR,IDLSTR,ISYMR,FREQR,
     *                  LISTX,IDLSTX,ISYMX,
     *                  WORK(KWXINT),WORK(KTHETAX),HAVEX,WORK(KXPROP),
     *                  LISTY,IDLSTY,ISYMY,
     *                  WORK(KWYINT),WORK(KTHETAY),HAVEY,WORK(KYPROP),
     *                  WORK(KT3AMP0),HAVEXY,WORK(KXYMAT),
     *                  WORK(KTAMP1),WORK(KTAMP2),
     *                  WORK(KTAMP3),WORK(KWXYOINT),WORK(KWXYVINT),
     *                  XLAMDP0,XLAMDH0,FOCK0,WORK(KFOCKD),
     *                  LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     *                  LOCDBG,HAVETB3T3,
     *                  WORK(KEND1),LWRK1)

        ! turn sign of t3 for compatibility with routines below
        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-ONE,WORK(KT3AMP0),1)
      END IF

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
     &                      WORK(KT3AMP0+IOFF),NT1AMX*NT1AMX,
     &                  ONE,WORK(KDAB0R),NVIRT)
        END DO

        CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NT1AMX*NVIRT,
     &              -HALF,WORK(KTAMP3),NT1AMX*NT1AMX*NVIRT,
     &                    WORK(KTBAR3),NT1AMX*NT1AMX*NVIRT,
     &               ZERO,WORK(KDIJR),NRHFT)

        CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NT1AMX*NVIRT,
     &              -HALF,WORK(KT3AMP0),NT1AMX*NT1AMX*NVIRT,
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
     &             (  WORK(KT0AMP-1+KEMFI) * WORK(KTAMP2-1+KDLAN) + 
     &                WORK(KTAMP2-1+KEMFI) * WORK(KT0AMP-1+KDLAN)   )
   
              END DO
            END DO

           END DO
          END DO
         END DO
        END DO

        IF (DO_XMMAT) THEN

          DO NFN = 1, NT1AMX
            DO M = 1, NRHFT
              DO I = 1, NRHFT

               KIMFN = ((NFN-1)*NRHFT + M-1) * NRHFT + I
               WORK(KXMMAT-1+KIMFN) = ZERO

                DO E = 1, NVIRT
                  NEI = NVIRT*(I-1) + E
                  NEM = NVIRT*(M-1) + E

                  DO NDL = 1, NT1AMX
                    KDLEMFN = ((NFN-1)*NT1AMX + NEM-1) * NT1AMX + NDL 
                    KDLEI   = (NEI-1) * NT1AMX + NDL

                    WORK(KXMMAT-1+KIMFN) = WORK(KXMMAT-1+KIMFN) -
     &                WORK(KTBAR3-1+KDLEMFN) * WORK(KTAMP2-1+KDLEI)
                  END DO
                END DO

              END DO
            END DO
          END DO

        END IF

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
      IF (DO_DIA) CALL DZERO(DIA,NVIRT*NRHFT)
      
      KWTILDE = KEND1
      KTHEOCC = KWTILDE + NT1AMX*NT1AMX
      KEND2   = KTHEOCC + NT1AMX*NT1AMX

      IF (SPLIT_WT) THEN
        KWBREVE = KEND2
        KTVIRT1 = KWBREVE + NT1AMX*NT1AMX
        KTVIRT2 = KTVIRT1 + NT1AMX*NT1AMX
        KEND2   = KTVIRT2 + NT1AMX*NT1AMX
      END IF

      IF (CUBIC) THEN
        KTHEOV  = KEND2
        KTVIRT3 = KTHEOV  + NT1AMX*NT1AMX
        KEND2   = KTVIRT3 + NT1AMX*NT1AMX
      END IF

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

      IF (.NOT.CUBIC) KWXYOINT = KWYINT

      CALL CCSDT_ADENVIR_NODDY(DAB,DIJ,DO_DIA,DIA,
     &                         WORK(KDAB0),WORK(KDIJ0),
     &                         SPLIT_WT,WORK(KWBAR),WORK(KWXYOINT),
     &                         WORK(KTHETAY),WORK(KT3AMP0),
     &                         CUBIC,WORK(KTHBAR),WORK(KWXYVINT),
     &                         WORK(KTHETAX),FREQR,
     &                         WORK(KXPROP),WORK(KYPROP),
     &                         WORK(KTHEOV),WORK(KFOCKD),
     &                         WORK(KWTILDE),WORK(KWBREVE),
     &                         WORK(KTHEOCC),WORK(KTVIRT1),
     &                         WORK(KTVIRT2),WORK(KTVIRT3),
     &                         WORK(KT0AMP),WORK(KTAMP2),WORK(KL2AMP),
     &                         WORK(KT0FMI),WORK(KTYFMI),
     &                         WORK(KT0MI), WORK(KTYMI),
     &                         WORK(KWDE),  WORK(KWED), WORK(KSCR))
 
      IF (LOCDBG) THEN
        write(lupri,*)'DAB in noddy after VIR'
        call output(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,lupri)
        write(lupri,*)'DIJ in noddy after VIR'
        call output(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,lupri)
        if (do_dia) then
          write(lupri,*)'DIA in noddy after VIR'
          call output(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,lupri)
        end if
      END IF

      IF (SPLIT_WT .OR. SPLIT_WBTB) THEN

        KWTILDE = KEND1
        KTHVIRT = KWTILDE + NVIRT*NVIRT*NT1AMX
        KSWAP   = KTHVIRT + NVIRT*NVIRT*NT1AMX
        KEND2   = KSWAP   + NVIRT
        
        IF (CUBIC) THEN
          KSCRVVVO = KEND2
          KEND2    = KSCRVVVO + NVIRT*NVIRT*NT1AMX
        END IF
        
        LWRK2   = LWORK   - KEND2
        IF (LWRK2.LT.0)
     &    CALL QUIT('Out of memory in CCSDT_ADEN_NODDY (3c)')

        CALL CCSDT_ADENOCC_NODDY(DAB,DIJ,DO_DIA,DIA,WORK(KWBAR),
     &                           WORK(KTHETAY),WORK(KL2AMP),
     &                           CUBIC,WORK(KTHBAR),FREQR,WORK(KTHETAX),
     &                           WORK(KWYINT),WORK(KWXINT),
     &                           WORK(KFOCKD),WORK(KXPROP),WORK(KYPROP),
     &                           WORK(KT3AMP0),WORK(KXYMAT),
     &                           WORK(KWTILDE),WORK(KTHVIRT),
     &                           WORK(KSWAP),WORK(KSCRVVVO))
 
        IF (LOCDBG) THEN
          write(lupri,*)'DAB in noddy after OCC'
          call output(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,lupri)
          write(lupri,*)'DIJ in noddy after OCC'
          call output(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,lupri)
          if (do_dia) then
            write(lupri,*)'DIA in noddy after OCC'
            call output(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,lupri)
          end if
        END IF

      END IF

*---------------------------------------------------------------------*
* add the one-index transformed zeroth-order density to DIA:
*---------------------------------------------------------------------*
      IF (DO_DIA) THEN

        ! read T^Y singles amplitudes from file 
        IOPT = 1
        Call CC_RDRSP(LISTR,IDLSTR,ISYMR,IOPT,MODEL,WORK(KTAMP1),DUMMY)

        CALL DGEMM('T','N',NVIRT,NRHFT,NVIRT,
     &             -ONE,WORK(KDAB0),NVIRT,WORK(KTAMP1),NVIRT,
     &              ONE,DIA,NVIRT)
 
        CALL DGEMM('N','T',NVIRT,NRHFT,NRHFT,
     &              ONE,WORK(KTAMP1),NVIRT,WORK(KDIJ0),NRHFT,
     &              ONE,DIA,NVIRT)  
 

        IF (LOCDBG) THEN
          write(lupri,*)'DAB(total) in noddy after OCC'
          call output(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,lupri)
          write(lupri,*)'DIJ(total) in noddy after OCC'
          call output(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,lupri)
          write(lupri,*)'DIA(total) in noddy after OCC'
          call output(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,lupri)
        END IF
  
        IF (HAVETB3T3) THEN
          CALL DGEMM('T','N',NVIRT,NRHFT,NVIRT,
     &             -ONE,WORK(KDAB0R),NVIRT,WORK(KTAMP1),NVIRT,
     &              ONE,WORK(KDIAR),NVIRT)
 
          CALL DGEMM('N','T',NVIRT,NRHFT,NRHFT,
     &                ONE,WORK(KTAMP1),NVIRT,WORK(KDIJ0R),NRHFT,
     &                ONE,WORK(KDIAR),NVIRT)  
        END IF

      END IF
*---------------------------------------------------------------------*
* quick hack solution for the M intermediates (use noddy-noddy):
*---------------------------------------------------------------------*
      IF (DO_XMMAT) THEN
        IF (.NOT. HAVETB3T3) CALL QUIT('DO_XMMAT needs HAVETB3T3!')

        CALL DCOPY(NRHFT*NRHFT*NT1AMX,WORK(KXMMAT),1,XMMAT,1)
      END IF
*---------------------------------------------------------------------*
* print debug output and return:
*---------------------------------------------------------------------*

      IF (PRINTDENS) THEN
        WRITE(LUPRI,*) 'A densities for:'
        WRITE(LUPRI,*) 'LISTL,IDLSTL:',LISTL,IDLSTL
        WRITE(LUPRI,*) 'LISTR,IDLSTR:',LISTR,IDLSTR

        WRITE(LUPRI,*) 'The virtual-virtual block DAB:'
        CALL OUTPUT(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*) 'The virtual-virtual block of the reference:'
         CALL OUTPUT(WORK(KDABR),1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
         CALL DAXPY(NVIRT*NVIRT,-ONE,DAB,1,WORK(KDABR),1)
         WRITE(LUPRI,*)'DAB norm^2(difference):',
     &     DDOT(NVIRT*NVIRT,WORK(KDABR),1,WORK(KDABR),1)
        END IF

        WRITE(LUPRI,*) 'The occupied-occupied block DIJ:'
        CALL OUTPUT(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*)'The occupied-occupied block of the reference:'
         CALL OUTPUT(WORK(KDIJR),1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
         CALL DAXPY(NRHFT*NRHFT,-ONE,DIJ,1,WORK(KDIJR),1)
         WRITE(LUPRI,*)'DIJ norm^2(difference):',
     &     DDOT(NRHFT*NRHFT,WORK(KDIJR),1,WORK(KDIJR),1)
        END IF

       IF (DO_DIA) THEN
        WRITE(LUPRI,*) 'The occupied-virtual block DIA:'
        CALL OUTPUT(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*)'The occupied-virtual block of the reference:'
         CALL OUTPUT(WORK(KDIAR),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,LUPRI)
         CALL DAXPY(NVIRT*NRHFT,-ONE,DIA,1,WORK(KDIAR),1)
         WRITE(LUPRI,*)'DIA norm^2(difference):',
     &     DDOT(NVIRT*NRHFT,WORK(KDIAR),1,WORK(KDIAR),1)
        END IF

        WRITE(LUPRI,*) 'The virtual-virtual block DAB0:'
        CALL OUTPUT(WORK(KDAB0),1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*) 'The virtual-virtual block of the reference:'
         CALL OUTPUT(WORK(KDAB0R),1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,LUPRI)
         CALL DAXPY(NVIRT*NVIRT,-ONE,WORK(KDAB0),1,WORK(KDAB0R),1)
         WRITE(LUPRI,*)'DAB0 norm^2(difference):',
     &     DDOT(NVIRT*NVIRT,WORK(KDAB0R),1,WORK(KDAB0R),1)
        END IF

        WRITE(LUPRI,*) 'The occupied-occupied block DIJ0:'
        CALL OUTPUT(WORK(KDIJ0),1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
        IF (HAVETB3T3) THEN
         WRITE(LUPRI,*)'The occupied-occupied block of the reference:'
         CALL OUTPUT(WORK(KDIJ0R),1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,LUPRI)
         CALL DAXPY(NRHFT*NRHFT,-ONE,WORK(KDIJ0),1,WORK(KDIJ0R),1)
         WRITE(LUPRI,*)'DIJ0 norm^2(difference):',
     &     DDOT(NRHFT*NRHFT,WORK(KDIJ0R),1,WORK(KDIJ0R),1)
        END IF
       END IF
      END IF

      CALL QEXIT('CCTADNOD')
      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_ADEN_NODDY                      *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_ADENVIR_NODDY(DAB,DIJ,DO_DIA,DIA,DAB0,DIJ0,
     &                               SPLIT_WT,WBAR,WYINT,THETAY,T0AMP,
     &                               CUBIC,THBAR,WXYV,THETAX,FREQXY,
     &                               XPROP,YPROP,THEOV,FOCKD,
     &                               WTILDE,WBREVE,THEOCC,
     &                               TVIRT1,TVIRT2,TVIRT3,
     &                               TAMP02,TAMPR2,L2AMP,TXFMI,TYFMI,
     &                               TXMI,TYMI,WDE,WED,SCR)
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

      LOGICAL CUBIC, DO_DIA, SPLIT_WT

#if defined (SYS_CRAY)
      REAL DIJ(NRHFT,NRHFT), DAB(NVIRT,NVIRT),DIA(NVIRT,NRHFT),
     &     DIJ0(*), DAB0(*)
      REAL WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WXYV(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETAY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WTILDE(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WBREVE(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THEOCC(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THEOV(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TVIRT1(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TVIRT2(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TVIRT3(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TAMP02(NT1AMX,NVIRT,NRHFT)
      REAL TAMPR2(NT1AMX,NVIRT,NRHFT)
      REAL L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TXFMI(NVIRT,NRHFT,NRHFT)
      REAL TYFMI(NVIRT,NRHFT,NRHFT)
      REAL TXMI(NRHFT,NRHFT), TYMI(NRHFT,NRHFT)
      REAL WDE(NVIRT,NRHFT,NRHFT)
      REAL WED(NVIRT,NRHFT,NRHFT)
      REAL SCR(NT1AMX*NRHFT), FOCKD(*)
      REAL XPROP(NORBT,NORBT),YPROP(NORBT,NORBT)
      REAL HALF, ONE, ZERO, DDOT, FREQXY
#else
      DOUBLE PRECISION DIJ(NRHFT,NRHFT),DAB(NVIRT,NVIRT),
     &                 DIA(NVIRT,NRHFT),DIJ0(*),DAB0(*)
      DOUBLE PRECISION WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WXYV(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETAY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WTILDE(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WBREVE(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THEOCC(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THEOV(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TVIRT1(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TVIRT2(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TVIRT3(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TAMP02(NT1AMX,NVIRT,NRHFT)
      DOUBLE PRECISION TAMPR2(NT1AMX,NVIRT,NRHFT)
      DOUBLE PRECISION L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TXFMI(NVIRT,NRHFT,NRHFT) 
      DOUBLE PRECISION TYFMI(NVIRT,NRHFT,NRHFT)
      DOUBLE PRECISION TXMI(NRHFT,NRHFT), TYMI(NRHFT,NRHFT)
      DOUBLE PRECISION WDE(NVIRT,NRHFT,NRHFT)
      DOUBLE PRECISION WED(NVIRT,NRHFT,NRHFT)
      DOUBLE PRECISION SCR(NT1AMX*NRHFT), FOCKD(*)
      DOUBLE PRECISION XPROP(NORBT,NORBT),YPROP(NORBT,NORBT)
      DOUBLE PRECISION HALF, ONE, ZERO, DDOT, FREQXY
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
             CALL DCOPY(NT1AMX,TAMP02(1,E,I),1,TXFMI(1,1,I),1)
             CALL DCOPY(NT1AMX,TAMPR2(1,E,I),1,TYFMI(1,1,I),1)
             DO M = 1, NRHFT
               NDM = NVIRT*(M-1) + D
               TXMI(M,I) = TAMP02(NDM,E,I)
               TYMI(M,I) = TAMPR2(NDM,E,I)
             END DO
           END DO

c          ------------------------------------------
c          D(IA) -= W^DE(fmnl) t^x(fm,Ei) t^y(an,DL)
c                      scr(n,i) 
c          ------------------------------------------
           IF (DO_DIA) THEN
             CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                  ONE,WDE,NT1AMX,TXFMI,NT1AMX,
     &                  ZERO,SCR,NRHFT)
             CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &                 -ONE,TAMPR2(1,D,L),NVIRT,SCR,NRHFT,
     &                  ONE,DIA,NVIRT)
           END IF

c          ------------------------------------------
c          D(IA) -= W^DE(fmnl) t^y(fm,Ei) t^x(an,DL)
c                      scr(n,i) 
c          ------------------------------------------
           IF (DO_DIA) THEN
            CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                 ONE,WDE,NT1AMX,TYFMI,NT1AMX,
     &                 ZERO,SCR,NRHFT)
            CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &                -ONE,TAMP02(1,D,L),NVIRT,SCR,NRHFT,
     &                 ONE,DIA,NVIRT)
           END IF

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
           IF (DO_DIA) THEN
            CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                 ONE,WDE,NT1AMX,TXFMI,NT1AMX,
     &                 ZERO,SCR,NRHFT)
            CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &                -ONE,TAMPR2(1,D,L),NVIRT,SCR,NRHFT,
     &                 ONE,DIA,NVIRT)
           END IF


c          ------------------------------------------
c          D(IA) -= W^DE(fnml) t^y(fi,Em) t^x(an,Dl)
c                      scr(n,i) 
c          ------------------------------------------
           IF (DO_DIA) THEN
            CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX,
     &                 ONE,WDE,NT1AMX,TYFMI,NT1AMX,
     &                 ZERO,SCR,NRHFT)
            CALL DGEMM('N','N',NVIRT,NRHFT,NRHFT,
     &                -ONE,TAMP02(1,D,L),NVIRT,SCR,NRHFT,
     &                 ONE,DIA,NVIRT)
           END IF


c          ------------------------------------------
c          D(IA) -= W^ED(fnml) t^x(Dm,Ei) t^y(fn,al)
c                        scr(fni)
c          ------------------------------------------
           IF (DO_DIA) THEN
            CALL DGEMM('N','N',NT1AMX,NRHFT,NRHFT,
     &                 ONE,WED,NT1AMX,TXMI,NRHFT,
     &                 ZERO,SCR,NT1AMX)
            CALL DGEMM('T','N',NVIRT,NRHFT,NT1AMX,
     &                -ONE,TAMPR2(1,1,L),NT1AMX,SCR,NT1AMX,
     &                 ONE,DIA,NVIRT)
           END IF

c          ------------------------------------------
c          D(IA) -= W^ED(fnml) t^y(Dm,Ei) t^x(fn,al)
c                        scr(fni)
c          ------------------------------------------
           IF (DO_DIA) THEN
            CALL DGEMM('N','N',NT1AMX,NRHFT,NRHFT,
     &                 ONE,WED,NT1AMX,TYMI,NRHFT,
     &                 ZERO,SCR,NT1AMX)
            CALL DGEMM('T','N',NVIRT,NRHFT,NT1AMX,
     &                -ONE,TAMP02(1,1,L),NT1AMX,SCR,NT1AMX,
     &                 ONE,DIA,NVIRT)
           END IF

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
           IF (SPLIT_WT) THEN
            DO J = 1, NRHFT
             CALL DCOPY(NT1AMX*NVIRT,WBAR(1,1,1,L,D,J),1,
     &                             WBREVE(1,1,1,J),1)
            END DO
           END IF

c          ==========================================================
c          build TVIRT2^DL intermediate:
c          ==========================================================
           IF (SPLIT_WT) THEN
           IF (.NOT. CUBIC) THEN

c           ------------------------------------------------------------
c           TVIRT2^DL(em,fi) = theta(e-bar m,fl,di)+theta(f-bar l,em,di)
c           ------------------------------------------------------------
            DO I = 1, NRHFT
             DO F = 1, NVIRT
              DO M = 1, NRHFT
              DO E = 1, NVIRT
               TVIRT2(E,M,F,I)=THETAY(E,M,F,L,D,I)+THETAY(F,L,E,M,D,I)
              END DO
              END DO
             END DO
            END DO

           ELSE
c           --------------------------------------------------------
c           TVIRT2^DL(em,fi) =  theta(f-bar-bar l, e m, d i)
c                              +theta(f-bar l, e-bar m, d i)
c                              +theta(f l, e-bar-bar m, d i)
c                              +theta(f-bar l-bar, e m-bar, d i-bar)
c           (uses TVIRT1 array as scratch area!)
c           --------------------------------------------------------
            CALL DZERO(TVIRT2,NT1AMX*NT1AMX)

            ! build TVIRT1^DL(em,fi) = thetay(em,fl,di)+thetay(fl,em,di)
            DO I = 1, NRHFT
             DO F = 1, NVIRT
              DO M = 1, NRHFT
              DO E = 1, NVIRT
               TVIRT1(E,M,F,I)=THETAY(E,M,F,L,D,I)+THETAY(F,L,E,M,D,I)
              END DO
              END DO
             END DO
            END DO

            ! one-index transform. e and f with XPROP
            CALL CCSDT_TX1_NODDY(TVIRT2,TVIRT1,XPROP,.FALSE.,.TRUE.,
     &                           NRHFT,NVIRT,NORBT)

            ! build TVIRT1^DL(em,fi) = thetax(em,fl,di)+thetax(fl,em,di)
            DO I = 1, NRHFT
             DO F = 1, NVIRT
              DO M = 1, NRHFT
              DO E = 1, NVIRT
               TVIRT1(E,M,F,I)=THETAX(E,M,F,L,D,I)+THETAX(F,L,E,M,D,I)
              END DO
              END DO
             END DO
            END DO

            ! one-index transform. e and f with XPROP
            CALL CCSDT_TX1_NODDY(TVIRT2,TVIRT1,YPROP,.FALSE.,.TRUE.,
     &                           NRHFT,NVIRT,NORBT)


            ! divide by orbital energies and remove forbid. elements:
            CALL CCSDT_DIVFCL_NODDY(TVIRT2,D,L,FREQXY,FOCKD,SCR,
     &                              NT1AMX,NVIRT,NRHFT,NORBT)
            CALL DSCAL(NT1AMX*NT1AMX,-1.0D0,TVIRT2,1)
            
            DO I = 1, NRHFT
             DO F = 1, NVIRT
              DO M = 1, NRHFT
              DO E = 1, NVIRT
               TVIRT2(E,M,F,I) = TVIRT2(E,M,F,I) + WXYV(F,L,E,M,D,I)
              END DO
              END DO
             END DO
            END DO

           END IF
           END IF

c          ==========================================================
c          build TVIRT1^DL / TVIRT3^DL intermediates:
c          ==========================================================
           IF (SPLIT_WT) THEN
           IF (.NOT. CUBIC) THEN
c            ----------------------------------------------------
c            TVIRT1^DL(fi,em) = theta(f-bar i,em,dl)
c            ----------------------------------------------------
             CALL DCOPY(NT1AMX*NT1AMX,THETAY(1,1,1,1,D,L),1,TVIRT1,1)

           ELSE
c            -----------------------------------------------------------
c            build TVIRT1^DL(em,fi) = theta(e m, f-bar-bar i, d l)
c                                    +theta(e-bar m, f-bar i, d l)
c            and   TVIRT3^DL(em,fi) = theta(e-bar m, f-bar i, d l)
c            -----------------------------------------------------------
             CALL DZERO(TVIRT1,NT1AMX*NT1AMX)
             CALL DZERO(TVIRT3,NT1AMX*NT1AMX)
            
             DO I = 1, NRHFT
              DO F = 1, NVIRT
               DO M = 1, NRHFT
               DO E = 1, NVIRT
                DO C = 1, NVIRT
                 TVIRT3(E,M,F,I) = TVIRT3(E,M,F,I) 
     &              - THETAY(E,M,C,I,D,L) * XPROP(NRHFT+F,NRHFT+C) 
     &              - THETAY(F,I,C,M,D,L) * XPROP(NRHFT+E,NRHFT+C) 
     &              - THETAX(E,M,C,I,D,L) * YPROP(NRHFT+F,NRHFT+C) 
     &              - THETAX(F,I,C,M,D,L) * YPROP(NRHFT+E,NRHFT+C) 

                 TVIRT1(E,M,F,I) = TVIRT1(E,M,F,I) 
     &              - THETAY(C,I,E,M,D,L) * XPROP(NRHFT+F,NRHFT+C) 
     &              - THETAX(C,I,E,M,D,L) * YPROP(NRHFT+F,NRHFT+C) 
                END DO
               END DO
               END DO
              END DO
             END DO

             CALL DAXPY(NT1AMX*NT1AMX,ONE,TVIRT3,1,TVIRT1,1)
            
             ! divide by orbital energies and remove forbid. elements:
             CALL CCSDT_DIVFCL_NODDY(TVIRT1,D,L,FREQXY,FOCKD,SCR,
     &                               NT1AMX,NVIRT,NRHFT,NORBT)
             CALL DSCAL(NT1AMX*NT1AMX,-1.0D0,TVIRT1,1)

             ! divide by orbital energies and remove forbid. elements:
             CALL CCSDT_DIVFCL_NODDY(TVIRT3,D,L,FREQXY,FOCKD,SCR,
     &                               NT1AMX,NVIRT,NRHFT,NORBT)
             CALL DSCAL(NT1AMX*NT1AMX,-1.0D0,TVIRT3,1)

           END IF
           END IF

c          ============================================================
c          build THEOCC^DL(em,fn) = w(em,fn,dl)+w(dl,em,fn)+w(fn,dl,em)
c          ============================================================
           DO E = 1, NVIRT
           DO F = 1, NVIRT

             DO M = 1, NRHFT
             DO N = 1, NRHFT
 
               THEOCC(E,M,F,N) = WYINT(E,M,F,N,D,L)+
     &                           WYINT(D,L,E,M,F,N)+WYINT(F,N,D,L,E,M)
 
             END DO
             END DO

             IF (CUBIC) THEN
               ! include contributions from 
               ! WXYV(fn,em,dl) = theta(D l-bar, e m-bar, f-bar n-bar) 
               ! and from
               ! WXYV(em,fn,dl) = theta(D l-bar, f n-bar, e-bar m-bar)
               DO M = 1, NRHFT
               DO N = 1, NRHFT
                 THEOCC(E,M,F,N) = THEOCC(E,M,F,N)+
     &                   WXYV(F,N,E,M,D,L) + WXYV(E,M,F,N,D,L)
               END DO
               END DO
             END IF

           END DO
           END DO
           
           IF (CUBIC) THEN
c           -----------------------------------------------------------
c           fetch THEOV^DL(em,fn) = theta(D l-bar,f n-bar,e-bar m-bar)
c           -----------------------------------------------------------
            DO E = 1, NVIRT
            DO F = 1, NVIRT
              DO M = 1, NRHFT
              DO N = 1, NRHFT
          
                THEOV(E,M,F,N) = WXYV(E,M,F,N,D,L)
          
              END DO
              END DO
            END DO
            END DO
           END IF
 
*---------------------------------------------------------------------*
* now we have the following intermediates in O^2 N^2 arrays:
*           
*          WTILDE^DL(em,an) = W^Da(emnl) + 1/2 W^De(anml)
*          WBREVE^DL(em,fj) = W^Df(emlj) 
*
*      for quadratic response:
*          TVIRT2^DL(em,fi) = theta(e-bar m,fl,di)+theta(f-bar l,em,di)
*          TVIRT1^DL(fi,em) = theta(f-bar i,em,dl)
*          THEOCC^DL(em,fn) = w(em,fn,dl)+w(dl,em,fn)+w(fn,dl,em)
*
*      for cubic response:
*          TVIRT2^DL(em,fi) =   theta(e-bar-bar m, f l,di)
*                             + theta(e-bar m, f-bar l,di)
*                             + theta(e m, f-bar-bar l,di)
*                             + theta(e m-bar, f-bar l-bar, d i-bar)
*          TVIRT1^DL(em,fi) =   theta(e m, f-bar-bar i,dl)
*                             + theta(e-bar m, f-bar i,dl)
*          TVIRT3^DL(em,fi) =   theta(e-bar m, f-bar i,dl)
*          THEOCC^DL(em,fn) = wxy(em,fn,dl)+wxy(dl,em,fn)+wxy(fn,dl,em)
*                             + theta(f-bar n-bar, e m-bar, D l-bar)
*                             + theta(e-bar m-bar, f n-bar, D l-bar)
*          THEOV^DL(em,fn)  = theta(e-bar m-bar, f n-bar, D l-bar)
*
* with these we can compute the contributions to D(ij) and D(ab) 
*        for quadratic response with 3 calls to DGEMM, and
*        for cubic response with 7 calls to DGEMM,
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

c          ---------------------------------------------------------
c          compute the contribution from <L_2|[E_ia,T_3^xy]|HF>
c          to D(ia) which can be done in the loop over virtuals:
c
c             D(ia) += [THEOCC^DL(em,ai)-THEOCC^DL(ei,am)] * L(em,DL)
c
c           for cubic response include in addition:
c             D(DL) += THEOV^DL(ai,em) * L(ai,em)
c             D(Di) -= THEOV^DL(am,ei) * L(am,eL)
c          -----------------------------------------------------------
           IF (DO_DIA) THEN
            DO I = 1, NRHFT
              DO A = 1, NVIRT
                DO M = 1, NRHFT
                  DIA(A,I) = DIA(A,I) + 
     &               DDOT(NVIRT,THEOCC(1,M,A,I),1,L2AMP(1,M,D,L),1) -
     &               DDOT(NVIRT,THEOCC(1,I,A,M),1,L2AMP(1,M,D,L),1)
               END DO
              END DO
            END DO
            
            IF (CUBIC) THEN
              DIA(D,L) = DIA(D,L) + DDOT(NT1AMX*NT1AMX,THEOV,1,L2AMP,1)
              
              CALL DGEMV('T',NT1AMX*NVIRT,NRHFT,
     *                   -ONE,THEOV,NT1AMX*NVIRT,L2AMP(1,1,1,L),1,
     *                    ONE,DIA(D,1),NVIRT)
            END IF
           END IF

c          ----------------------------------------------------------
c          add TVIRT1 to THEOCC: THEOCC^DL(em,fi) += TVIRT1^DL(fi,em)
c          (skipped for cubic response)
c          ----------------------------------------------------------
           IF (SPLIT_WT .AND. (.NOT. CUBIC)) THEN
             DO M = 1, NRHFT
               DO E = 1, NVIRT
                 DO I = 1, NRHFT
                   DO F = 1, NVIRT
                     THEOCC(E,M,F,I) = THEOCC(E,M,F,I)+TVIRT1(F,I,E,M)
                   END DO
                 END DO
               END DO
             END DO
           END IF

c          -----------------------------------------------------------
c          D(ij)+=WTILDE^DL(em,fj)*[THEOCC^DL(em,fi)+TVIRT1^DL(fi,em)]
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                          -ONE,THEOCC,NT1AMX*NVIRT,
     &                             WTILDE,NT1AMX*NVIRT,
     &                           ONE,DIJ,NRHFT)

c          -----------------------------------------------------------
c          D(ij) -= WBREVE^DL(em,fj) * TVIRT2^DL(em,fi)
c          -----------------------------------------------------------
           IF (SPLIT_WT) THEN
             CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                            -ONE,TVIRT2,NT1AMX*NVIRT,
     &                                 WBREVE,NT1AMX*NVIRT,
     &                             ONE,DIJ,NRHFT)
           END IF

c          -----------------------------------------------------------
c          D0(ij) -= WTILDE^DL(em,fj) * T0AMP^DL(em,fi)
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                          -ONE,T0AMP(1,1,1,1,D,L),NT1AMX*NVIRT,
     &                             WTILDE,NT1AMX*NVIRT,
     &                           ONE,DIJ0,NRHFT)

*---------------------------------------------------------------------*
* add contributions without counterpart in quadratic response:
*---------------------------------------------------------------------*
         IF (CUBIC) THEN

c          -----------------------------------------------------------
c          D(ij) -= W^DL(em,fj) * TVIRT1^DL(em,fi)
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                        -ONE,TVIRT1,NT1AMX*NVIRT,
     &                  WBAR(1,1,1,1,D,L),NT1AMX*NVIRT,
     &                         ONE,DIJ,NRHFT)

c          ---------------------------------------------------------
c          resort WBREVE^DL(em,fj) = Wbar^De(fj,ml) = WBAR(fj,em,DL)
c          (i.e. transpose the pair indeces em, fj)
c          ---------------------------------------------------------
           DO F = 1, NVIRT
            DO J = 1, NRHFT
              CALL DCOPY(NT1AMX,WBAR(F,J,1,1,D,L),NT1AMX,
     &                          WBREVE(1,1,F,J),1)
            END DO
           END DO

c          -----------------------------------------------------------
c          D(ij) -= 1/2 WBREVE^DL(em,fj) * [THEOV^DL(em,fi) + ...]
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                       -HALF,THEOV,NT1AMX*NVIRT,
     &                             WBREVE,NT1AMX*NVIRT,
     &                         ONE,DIJ,NRHFT)
         
c          -----------------------------------------------------------
c          get in WTILDE^DL(em,fj) = wbar^D(ef;l-bar m j)
c          -----------------------------------------------------------
           DO J = 1, NRHFT
           DO F = 1, NVIRT
             DO M = 1, NRHFT
             DO E = 1, NVIRT
               WTILDE(E,M,F,J) = WBAR(D,L,E,M,F,J) - THBAR(D,L,E,M,F,J)
             END DO
             END DO
           END DO
           END DO

c          -----------------------------------------------------------
c          D(ij) -= W^ef(DL,jm) * Theta^DL(e-bar m, f-bar i)
c          -----------------------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NT1AMX*NVIRT,
     &                        -ONE,TVIRT3,NT1AMX*NVIRT,
     &                             WTILDE,NT1AMX*NVIRT,
     &                         ONE,DIJ,NRHFT)

c          -----------------------------------------------------------
c          D(ab) +=  [ 1/2 W^De(an,ml) + wbar^D(ea; l-bar m n) ] 
c                              * THEOV^DL(em,bn)
c          -----------------------------------------------------------
           CALL DAXPY(NT1AMX*NT1AMX,HALF,WBREVE,1,WTILDE,1)
           DO N = 1, NRHFT
             CALL DGEMM('T','N',NVIRT,NVIRT,NT1AMX,
     &                   ONE,WTILDE(1,1,1,N),NT1AMX,
     &                       THEOV(1,1,1,N),NT1AMX,
     &                   ONE,DAB,NVIRT)
           END DO

          END IF
        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_ADENVIR_NODDY                   *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_TX1_NODDY(TX,TAMP,XMAT,LOCC,LVIRT,
     &                           NRHFT,NVIRT,NORBT)
*---------------------------------------------------------------------*
c     TX(em,fn) +=  sum_k T(ek,fn) X_km - sum_c T(cm,fn) X_ec 
c                  +sum_k T(em,fk) X_kn - sum_c T(em,cn) X_fc 
c     LOCC  = .TRUE.  :  include transformation of occupied index
c     LVIRT = .TRUE.  :  include transformation of virtual index
*---------------------------------------------------------------------*
      IMPLICIT NONE

      LOGICAL LOCC, LVIRT
      INTEGER NVIRT,NRHFT,NORBT,E,F,M,N,C,K

#if defined (SYS_CRAY)
      REAL TAMP(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL TX(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL XMAT(NORBT,NORBT)
#else
      DOUBLE PRECISION TAMP(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION TX(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION XMAT(NORBT,NORBT)
#endif

      DO F = 1, NVIRT
      DO N = 1, NRHFT

       IF (LVIRT) THEN
        DO E = 1, NVIRT
        DO M = 1, NRHFT
         DO C = 1, NVIRT
          TX(E,M,F,N) = TX(E,M,F,N) - 
     &       ( TAMP(C,M,F,N)*XMAT(NRHFT+E,NRHFT+C) +
     &         TAMP(E,M,C,N)*XMAT(NRHFT+F,NRHFT+C)  )
         END DO
        END DO
        END DO
       END IF

       IF (LOCC) THEN
        DO E = 1, NVIRT
        DO M = 1, NRHFT
         DO K = 1, NRHFT
          TX(E,M,F,N) = TX(E,M,F,N) + 
     &       ( TAMP(E,K,F,N)*XMAT(K,M) + TAMP(E,M,F,K)*XMAT(K,N) )
         END DO
        END DO
        END DO
       END IF

      END DO
      END DO 

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_DIVFCL_NODDY(TAMP,NC,NK,FREQ,FOCKD,SCR1,
     &                              NT1AMX,NVIRT,NRHFT,NORBT)
*---------------------------------------------------------------------*
c     1) T^CK(ai,bj) =  T^CK(ai,bj) / (eps_aibjck - freq)
c     2) remove forbidden elements
*---------------------------------------------------------------------*
      IMPLICIT NONE

      INTEGER NVIRT,NRHFT,NORBT,NT1AMX,NAI,NBJ,NCK,NK,NC
      INTEGER NA,NB,NAK,NBK,NI,NJ,NCI,NCJ

#if defined (SYS_CRAY)
      REAL TAMP(NT1AMX,NT1AMX), FOCKD(NORBT), SCR1(NT1AMX)
      REAL FREQ
#else
      DOUBLE PRECISION TAMP(NT1AMX,NT1AMX), FOCKD(NORBT), SCR1(NT1AMX)
      DOUBLE PRECISION FREQ
#endif

      DO NI = 1,NRHFT
         DO NA = 1,NVIRT
            NAI = NVIRT*(NI-1) + NA
            SCR1(NAI) = FOCKD(NRHFT+NA) - FOCKD(NI)
         END DO
      END DO

      NCK = NVIRT*(NK-1) + NC
 
      DO NBJ = 1, NT1AMX
        DO NAI = 1, NT1AMX
          TAMP(NAI,NBJ) = TAMP(NAI,NBJ) /
     &       (SCR1(NAI) + SCR1(NBJ) + SCR1(NCK) - FREQ)
        END DO
      END DO

      DO NB = 1, NVIRT
        NBK = NVIRT*(NK-1) + NB
        DO NA = 1, NVIRT
           NAK = NVIRT*(NK-1) + NA
           TAMP(NAK,NBK) = 0.0d0
         ENDDO
      ENDDO
      DO NJ = 1, NRHFT
         NCJ = NVIRT*(NJ-1) + NC
         DO NI = 1, NRHFT
           NCI = NVIRT*(NI-1) + NC
           TAMP(NCI,NCJ) = 0.0d0
         ENDDO
       ENDDO

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_ADENOCC_NODDY(DAB,DIJ,DO_DIA,DIA,
     &                               WBAR,THETA,L2AMP,
     &                               CUBIC,THBAR,FREQR2,THETAX,
     &                               WYINT,WXINT,FOCKD,XPROP,YPROP,
     &                               T0AMP,XYMAT,
     &                               WTILDE,THVIRT,SWAP,SCR)
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

      LOGICAL CUBIC, DO_DIA

#if defined (SYS_CRAY)
      REAL DIJ(*), DAB(NVIRT,NVIRT), DIA(NVIRT,NRHFT)
      REAL T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WXINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL WTILDE(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL THVIRT(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL SCR(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      REAL FREQR2, XPROP(NORBT,NORBT), YPROP(NORBT,NORBT), FOCKD(*)
      REAL SWAP(NVIRT), XYMAT(NORBT,NORBT)
      REAL HALF, ONE, DDOT
#else
      DOUBLE PRECISION DIJ(*), DAB(NVIRT,NVIRT), DIA(NVIRT,NRHFT)
      DOUBLE PRECISION T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WYINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WXINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION WTILDE(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION THVIRT(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION SCR(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION L2AMP(NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION FREQR2, XPROP(NORBT,NORBT), YPROP(NORBT,NORBT)
      DOUBLE PRECISION SWAP(NVIRT), FOCKD(*), XYMAT(NORBT,NORBT)
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
        
           IF (.NOT. CUBIC) THEN
c            ------------------------------------------------
c            resort THVIRT^LM(de,an) = theta(d-bar l,em,an) + 
c                                      theta(dl,e-bar m,an) 
c            ------------------------------------------------
             DO N = 1, NRHFT
              DO A = 1, NVIRT
               DO E = 1, NVIRT
                DO D = 1, NVIRT
                 THVIRT(D,E,A,N)=THETA(D,L,E,M,A,N)+THETA(E,M,D,L,A,N)
                END DO
               END DO
              END DO
             END DO

           ELSE 
c            ----------------------------------------------------
c            for cubic response transform THVIRT further and add
c            add contributions to obtain in THVIRT^LM the 
c            theta intermediate with double bar on all virtuals
c            ----------------------------------------------------
             CALL DZERO(THVIRT,NVIRT*NVIRT*NT1AMX)
  
             ! get SCR(de,an) = theta(d-bar^y l,em,an) + 
             !    theta(dl,e-bar^y m,an) + theta(dl,em,a-bar^y n)
             DO N = 1, NRHFT
              DO A = 1, NVIRT
               DO E = 1, NVIRT
                DO D = 1, NVIRT
                  SCR(D,E,A,N) = THETA(D,L,E,M,A,N) + 
     &                 THETA(E,M,A,N,D,L) + THETA(A,N,D,L,E,M)
                END DO
               END DO
              END DO
             END DO

             ! transform virtual indeces with X
             call CCSDT_TXVIR_NODDY(THVIRT,SCR,XPROP,0,
     &                              NVIRT,NRHFT,NORBT)
      
             ! get SCR(de,an) = theta(d-bar^x l,em,an) + 
             !    theta(dl,e-bar^x m,an) + theta(dl,em,a-bar^x n)
             DO N = 1, NRHFT
              DO A = 1, NVIRT
               DO E = 1, NVIRT
                DO D = 1, NVIRT
                  SCR(D,E,A,N) = THETAX(D,L,E,M,A,N) + 
     &                 THETAX(E,M,A,N,D,L) + THETAX(A,N,D,L,E,M)
                END DO
               END DO
              END DO
             END DO

             ! transform virtual indeces with Y
             call CCSDT_TXVIR_NODDY(THVIRT,SCR,YPROP,0,
     &                              NVIRT,NRHFT,NORBT)
      
             ! divide by orbital energies and remove forbidden entries
             CALL CCSDT_DIVVIR_NODDY(THVIRT,L,M,FREQR2,
     &                               FOCKD,FOCKD(NRHFT+1),NVIRT,NRHFT)
      
             CALL DSCAL(NVIRT*NVIRT*NT1AMX,-1.0D0,THVIRT,1)
           END IF

c          --------------------------------------------
c          D(ij) += WTILDE^LM(de,fj) * THVIRT^LM(de,fi)
c          --------------------------------------------
           CALL DGEMM('T','N',NRHFT,NRHFT,NVIRT*NVIRT*NVIRT,
     &                       -HALF,THVIRT,NVIRT*NVIRT*NVIRT,
     &                             WTILDE,NVIRT*NVIRT*NVIRT,
     &                         ONE,DIJ,NRHFT)

           IF (.NOT. CUBIC) THEN
c           --------------------------------------------
c           add THVIRT^LM(de,an) += theta(dl,em,a-bar n)
c           --------------------------------------------
            DO N = 1, NRHFT
             DO A = 1, NVIRT
              DO E = 1, NVIRT
               DO D = 1, NVIRT
                THVIRT(D,E,A,N) = THVIRT(D,E,A,N) + THETA(A,N,E,M,D,L) 
               END DO
              END DO
             END DO
            END DO
           END IF
        
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
           IF (DO_DIA) THEN
            DO I = 1, NRHFT
              DO A = 1, NVIRT
                DO E = 1, NVIRT
                  DIA(A,I) = DIA(A,I) + 
     &               DDOT(NVIRT,THVIRT(1,E,A,I),1,L2AMP(1,L,E,M),1) -
     &               DDOT(NVIRT,THVIRT(1,A,E,I),1,L2AMP(1,L,E,M),1)
               END DO
              END DO
            END DO
           END IF

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

           IF (CUBIC) THEN
c            ------------------------------------------------------
c            resort WTILDE^LM(de,an) = theta-bar(d-bar L, e M, a n)
c            ------------------------------------------------------
             DO N = 1, NRHFT
               DO A = 1, NVIRT
                 DO E = 1, NVIRT
                   DO D = 1, NVIRT
                     WTILDE(D,E,A,N) = THBAR(D,L,E,M,A,N)
                   END DO
                 END DO
               END DO
             END DO

c            ------------------------------------------------------
c            set up THVIRT^LM(de,an) = theta(e-bar M, a-bar n, d L)
c            ------------------------------------------------------
             CALL DZERO(THVIRT,NVIRT*NVIRT*NT1AMX)

             DO N = 1, NRHFT
             DO A = 1, NVIRT
               DO E = 1, NVIRT
               DO D = 1, NVIRT
                 DO C = 1, NVIRT 
                   THVIRT(D,E,A,N) = THVIRT(D,E,A,N) -
     &              ( THETA(A,N,C,M,D,L) * XPROP(NRHFT+E,NRHFT+C) +
     &                THETA(E,M,C,N,D,L) * XPROP(NRHFT+A,NRHFT+C) +
     &               THETAX(A,N,C,M,D,L) * YPROP(NRHFT+E,NRHFT+C) +
     &               THETAX(E,M,C,N,D,L) * YPROP(NRHFT+A,NRHFT+C)  )
                 END DO
               END DO
               END DO
             END DO
             END DO

             ! divide by orbital energies and remove forbidden entries
             CALL CCSDT_DIVVIR_NODDY(THVIRT,L,M,FREQR2,
     &                               FOCKD,FOCKD(NRHFT+1),NVIRT,NRHFT)
      
             CALL DSCAL(NVIRT*NVIRT*NT1AMX,-1.0D0,THVIRT,1)


c            ----------------------------------------------------
c            D(ij) += theta-bar(d-bar L,eM,aj) * THVIRT^LM(de,fi)
c            ----------------------------------------------------
             CALL DGEMM('T','N',NRHFT,NRHFT,NVIRT*NVIRT*NVIRT,
     &                          -ONE,THVIRT,NVIRT*NVIRT*NVIRT,
     &                               WTILDE,NVIRT*NVIRT*NVIRT,
     &                           ONE,DIJ,NRHFT)


c            ------------------------------------------------------
c            set up THVIRT^LM(da,en) = theta(d L, a M, e-bar n-bar)
c            ------------------------------------------------------
             CALL DZERO(THVIRT,NVIRT*NVIRT*NT1AMX)

             DO N = 1, NRHFT
             DO A = 1, NVIRT
               DO E = 1, NVIRT
               DO D = 1, NVIRT
                 DO C = 1, NVIRT 
                   THVIRT(D,A,E,N) = THVIRT(D,A,E,N) 
     &              - WYINT(C,N,A,M,D,L) * XPROP(NRHFT+E,NRHFT+C) 
     &              - WXINT(C,N,A,M,D,L) * YPROP(NRHFT+E,NRHFT+C) 
     &              + T0AMP(C,N,A,M,D,L) * XYMAT(NRHFT+E,NRHFT+C)
                 END DO
                 DO K = 1, NRHFT 
                   THVIRT(D,A,E,N) = THVIRT(D,A,E,N) 
     &              +  THETA(E,K,A,M,D,L) * XPROP(K,N)
     &              + THETAX(E,K,A,M,D,L) * YPROP(K,N)
                 END DO
               END DO
               END DO
             END DO
             END DO

             ! divide by orbital energies and remove forbidden entries
             CALL CCSDT_DIVVIR_NODDY(THVIRT,L,M,FREQR2,
     &                               FOCKD,FOCKD(NRHFT+1),NVIRT,NRHFT)
      
             CALL DSCAL(NVIRT*NVIRT*NT1AMX,-1.0D0,THVIRT,1)

c            --------------------------------------------------
c            D(ab) += theta-bar(d-bar L, a M, e n) * 
c                                theta(d L, b M, e-bar n-bar)
c            --------------------------------------------------
             DO N = 1, NRHFT
               DO E = 1, NVIRT
                 CALL DGEMM('T','N',NVIRT,NVIRT,NVIRT,
     &                      ONE,WTILDE(1,1,E,N),NVIRT,
     &                          THVIRT(1,1,E,N),NVIRT,
     &                      ONE,DAB,NVIRT)
               END DO
             END DO


c            --------------------------------------------------------
c            set up THVIRT^LM(de,an) = theta(d L, e-bar M, a n-bar)
c            --------------------------------------------------------
             CALL DZERO(THVIRT,NVIRT*NVIRT*NT1AMX)

             DO N = 1, NRHFT
             DO A = 1, NVIRT
               DO E = 1, NVIRT
               DO D = 1, NVIRT
                 DO C = 1, NVIRT 
                   THVIRT(D,E,A,N) = THVIRT(D,E,A,N) -
     &               WYINT(A,N,C,M,D,L) * XPROP(NRHFT+E,NRHFT+C) -
     &               WXINT(A,N,C,M,D,L) * YPROP(NRHFT+E,NRHFT+C)
                 END DO
                 DO K = 1, NRHFT 
                   THVIRT(D,E,A,N) = THVIRT(D,E,A,N) + 
     &               THETA(E,M,A,K,D,L) * XPROP(K,N) +
     &               THETAX(E,M,A,K,D,L) * YPROP(K,N)
                 END DO
               END DO
               END DO
             END DO
             END DO

             ! divide by orbital energies and remove forbidden entries
             CALL CCSDT_DIVVIR_NODDY(THVIRT,L,M,FREQR2,
     &                               FOCKD,FOCKD(NRHFT+1),NVIRT,NRHFT)
      
             CALL DSCAL(NVIRT*NVIRT*NT1AMX,-1.0D0,THVIRT,1)

c            --------------------------------------------------
c            D(ab) += theta-bar(d-bar L, e M, a n) * 
c                                theta(d L, e-bar M, b n-bar)
c            --------------------------------------------------
             DO N = 1, NRHFT
               CALL DGEMM('T','N',NVIRT,NVIRT,NVIRT*NVIRT,
     &                    ONE,WTILDE(1,1,1,N),NVIRT*NVIRT,
     &                        THVIRT(1,1,1,N),NVIRT*NVIRT,
     &                    ONE,DAB,NVIRT)
             END DO

c            ------------------------------------------------------
c            resort WTILDE^LM(ae,dn) = theta-bar(d-bar n, e M, a L)
c            ------------------------------------------------------
             DO N = 1, NRHFT
               DO A = 1, NVIRT
                 DO E = 1, NVIRT
                   DO D = 1, NVIRT
                     WTILDE(A,E,D,N) = THBAR(D,N,E,M,A,L)
                   END DO
                 END DO
               END DO
             END DO

c            --------------------------------------------------
c            D(ab) += theta-bar(a L, e M, d-bar n) * 
c                                theta(b L, e-bar M, d n-bar)
c            --------------------------------------------------
             CALL DGEMM('N','T',NVIRT,NVIRT,NVIRT*NVIRT*NRHFT,
     &                           ONE,WTILDE,NVIRT,
     &                               THVIRT,NVIRT,
     &                              ONE,DAB,NVIRT)

           END IF

        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*             END OF SUBROUTINE CCSDT_ADENOCC_NODDY                   *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_TXVIR_NODDY(TX,TAMP,XPROP,IOPT,
     &                             NVIRT,NRHFT,NORBT)
*---------------------------------------------------------------------*
c
c    IOPT = 1 : TX(de,fn) -= sum_c X_dc*T(ce,fn)   (transf. 1 index)
c
c    IOPT = 2 : TX(de,fn) -= sum_c X_ec*T(dc,fn)   (transf. 2 index)
c
c    IOPT = 3 : TX(de,fn) -= sum_c X_fc*T(de,cn)   (transf. 3 index)
c
c    IOPT = 0 : transform all three indeces
c
*---------------------------------------------------------------------*
      IMPLICIT NONE

      INTEGER NVIRT, NRHFT, NORBT, IOPT, D, E, F, N, C

#if defined (SYS_CRAY)
      REAL TAMP(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL TX(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL XPROP(NORBT,NORBT)
#else
      DOUBLE PRECISION TAMP(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION TX(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION XPROP(NORBT,NORBT)
#endif

      DO N = 1, NRHFT
       DO F = 1, NVIRT

        IF (IOPT.EQ.0 .OR. IOPT.EQ.1) THEN
          DO D = 1, NVIRT
           DO E = 1, NVIRT
            DO C = 1, NVIRT
              TX(D,E,F,N) = TX(D,E,F,N) - 
     &               XPROP(NRHFT+D,NRHFT+C) * TAMP(C,E,F,N)
            END DO
           END DO
          END DO
        END IF
      
        IF (IOPT.EQ.0 .OR. IOPT.EQ.2) THEN
          DO D = 1, NVIRT
           DO E = 1, NVIRT
            DO C = 1, NVIRT
              TX(D,E,F,N) = TX(D,E,F,N) - 
     &               XPROP(NRHFT+E,NRHFT+C) * TAMP(D,C,F,N)
            END DO
           END DO
          END DO
        END IF
      
        IF (IOPT.EQ.0 .OR. IOPT.EQ.3) THEN
          DO D = 1, NVIRT
           DO E = 1, NVIRT
            DO C = 1, NVIRT
              TX(D,E,F,N) = TX(D,E,F,N) - 
     &               XPROP(NRHFT+F,NRHFT+C) * TAMP(D,E,C,N)
            END DO
           END DO
          END DO
        END IF

       END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_DIVVIR_NODDY(TAMP,L,M,FREQ,FCKDOCC,FCKDVIR,
     &                              NVIRT,NRHFT)
*---------------------------------------------------------------------*
c     1) T^LM(defn) =  T^LM(defn) / (eps_dl + eps_em + eps_fn - freq)
c     2) remove forbidden elements
*---------------------------------------------------------------------*
      IMPLICIT NONE

      INTEGER NVIRT,NRHFT,L,M,N,D,E,F

#if defined (SYS_CRAY)
      REAL TAMP(NVIRT,NVIRT,NVIRT,NRHFT)
      REAL FCKDOCC(NRHFT), FCKDVIR(NVIRT), FREQ, FLMN
#else
      DOUBLE PRECISION TAMP(NVIRT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION FCKDOCC(NRHFT), FCKDVIR(NVIRT), FREQ, FLMN
#endif

      DO N = 1, NRHFT
        FLMN = FCKDOCC(L) + FCKDOCC(M) + FCKDOCC(N) + FREQ
        DO F = 1, NVIRT
          DO E = 1, NVIRT
            DO D = 1, NVIRT
               TAMP(D,E,F,N) = TAMP(D,E,F,N) /
     &           (FCKDVIR(D) + FCKDVIR(E) + FCKDVIR(F) - FLMN)
            END DO
          END DO
        END DO
      END DO

      IF (L.EQ.M) THEN
        DO F = 1, NVIRT
          DO E = 1, NVIRT
            DO D = 1, NVIRT
             TAMP(D,E,F,L) = 0.0d0
            END DO
          END DO
        END DO
      END IF

      DO N = 1, NRHFT
         DO D = 1, NVIRT
           TAMP(D,D,D,N) = 0.0d0
         END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_WXYVINT_NODDY(WXYV,THETAY,THETAX,YPROP,XPROP,
     &                               NVIRT,NRHFT,NORBT)
*---------------------------------------------------------------------*
c     build wxy(em,fn,dl) = wxy(em,fn,dl) +
c        thetay(ek,fn,dl) * xprop(km) + thetax(ek,fn,dl) * yprop(km) +
c        thetay(em,fk,dl) * xprop(kn) + thetax(em,fk,dl) * yprop(kn) +
c        thetay(em,fn,dk) * xprop(kl) + thetax(em,fn,dk) * yprop(kl) 
*---------------------------------------------------------------------*
      IMPLICIT NONE

      INTEGER NRHFT,NVIRT,NORBT, E, M, F, N, D, L, K

#if defined (SYS_CRAY)
      REAL   WXYV(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETAY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      REAL XPROP(NORBT,NORBT), YPROP(NORBT,NORBT)
#else
      DOUBLE PRECISION WXYV(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETAY(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETAX(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION XPROP(NORBT,NORBT), YPROP(NORBT,NORBT)
#endif

      DO D = 1, NVIRT
      DO L = 1, NRHFT
        DO F = 1, NVIRT
        DO N = 1, NRHFT
          DO E = 1, NVIRT
          DO M = 1, NRHFT
            
            DO K = 1, NRHFT
               WXYV(E,M,F,N,D,L) = WXYV(E,M,F,N,D,L) +
     &           THETAY(E,K,F,N,D,L) * XPROP(K,M) + 
     &           THETAX(E,K,F,N,D,L) * YPROP(K,M) + 
     &           THETAY(E,M,F,K,D,L) * XPROP(K,N) + 
     &           THETAX(E,M,F,K,D,L) * YPROP(K,N) + 
     &           THETAY(E,M,F,N,D,K) * XPROP(K,L) + 
     &           THETAX(E,M,F,N,D,K) * YPROP(K,L) 
            END DO
      
          END DO
          END DO
        END DO
        END DO
      END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_T32INT_NODDY(LISTR,IDLSTR,ISYMR,FREQR,
     *                  LISTX,IDLSTX,ISYMX,
     *                  WXINT,THETAX,HAVEX,XPROP,
     *                  LISTY,IDLSTY,ISYMY,
     *                  WYINT,THETAY,HAVEY,YPROP,
     *                  T3AMP0,HAVEXY,XYMAT,
     *                  T1AMXY,T2AMXY,TAMXY3,WXYOINT,WXYVINT,
     *                  XLAMDP0,XLAMDH0,FOCK0,FOCKD,
     *                  LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     *                  LOCDBG,HAVETB3T3,
     *                  WORK,LWORK)
*---------------------------------------------------------------------*
*     calculate a number of intermediates needed in CCSDT_ADEN_NODDY
*     to represent the second-order amplitudes T3^XY
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL LOCDBG, HAVETB3T3, HAVEX, HAVEY, HAVEXY
      CHARACTER LISTX*3, LISTY*3, LISTR*3,
     &          FNDELD*(*), FNDKBC*(*), FNCKJD*(*)
      INTEGER IDLSTX, ISYMX, IDLSTY, ISYMY, IDLSTR, ISYMR, LWORK
      INTEGER LUDELD, LUDKBC, LUCKJD

#if defined (SYS_CRAY)
      REAL WXINT(*), THETAX(*), XPROP(*)
      REAL WYINT(*), THETAY(*), YPROP(*)
      REAL T3AMP0(*), XYMAT(*), WORK(*)
      REAL TAMXY3(*), WXYOINT(*), WXYVINT(*)
      REAL T1AMXY(*), T2AMXY(*)
      REAL XLAMDP0(*), XLAMDH0(*), FOCK0(*), FOCKD(*)
      REAL FREQR, DDOT
#else
      DOUBLE PRECISION WXINT(*), THETAX(*), XPROP(*)
      DOUBLE PRECISION WYINT(*), THETAY(*), YPROP(*)
      DOUBLE PRECISION T3AMP0(*), XYMAT(*), WORK(*)
      DOUBLE PRECISION TAMXY3(*), WXYOINT(*), WXYVINT(*)
      DOUBLE PRECISION T1AMXY(*), T2AMXY(*)
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*), FOCKD(*)
      DOUBLE PRECISION FREQR, DDOT
#endif

      CHARACTER*10 MODEL
      INTEGER KEND1, KLAMPX, KLAMHX, KLAMPY, KLAMHY, LWRK1, KEND2,
     &        KINT1SXY, KINT2SXY, KINT1SX, KINT2SX, KINT1SY, KINT2SY,
     &        KT2AMP0, KINT1S0, KINT2S0, KSCR1, KTAMPX3, KTAMPY3,
     &        LWRK2, IOPT, LUSIFC, KTHEOCC, KWXYO, KWXYV,
     &        KT2AMPX, KT2AMPY, KTHVIRT, KTXY3A, KTXY3B

      CALL QENTER('CCT32INT')
*---------------------------------------------------------------------*
*     allocate memory for local arrays:
*---------------------------------------------------------------------*
      KSCR1    = 1
      KEND1    = KSCR1 + NT1AMX

      KLAMPX   = KEND1
      KLAMHX   = KLAMPX + NLAMDT
      KEND1    = KLAMHX + NLAMDT

      KLAMPY   = KEND1
      KLAMHY   = KLAMPY + NLAMDT
      KEND1    = KLAMHY + NLAMDT

      KT2AMP0  = KEND1
      KEND1    = KT2AMP0 + NT1AMX*NT1AMX

      KT2AMPY  = KEND1
      KT2AMPX  = KT2AMPY + NT1AMX*NT1AMX
      KEND1    = KT2AMPX + NT1AMX*NT1AMX

      KTXY3A  = KEND1
      KTXY3B  = KTXY3A  + NT1AMX*NT1AMX*NT1AMX
      KEND1   = KTXY3B  + NT1AMX*NT1AMX*NT1AMX

      IF (HAVETB3T3) THEN
        KTAMPX3  = KEND1
        KTAMPY3  = KTAMPX3  + NT1AMX*NT1AMX*NT1AMX
        KEND1    = KTAMPY3  + NT1AMX*NT1AMX*NT1AMX
      END IF

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) CALL QUIT('Out of memory in CCSDT_T32INT_NODDY')

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (FOCKD(I), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     Get all the intermediates needed to compute the contributions
*     from the A{x} t^y and A{y} t^x to T^xy, i.e.
*     the w^x, theta^x, w^y, theta^y intermediates and the property
*     integrals and the one-index transformed property integrals
*---------------------------------------------------------------------*
      CALL CCSDT_T32_AAPREP_NODDY(
     *            LISTX,IDLSTX,ISYMX,WXINT,THETAX,XPROP,
     *            WORK(KLAMPX),WORK(KLAMHX),WORK(KTAMPX3),
     *            LISTY,IDLSTY,ISYMY,WYINT,THETAY,YPROP,
     *            WORK(KLAMPY),WORK(KLAMHY),WORK(KTAMPY3),
     *            T3AMP0,XYMAT,XLAMDP0,XLAMDH0,FOCK0,
     *            LUDELD,FNDELD,LUDKBC,FNDKBC,LUCKJD,FNCKJD,
     *            LOCDBG,HAVETB3T3,WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Compute the "easy" contributions to T^XY_3 from the B matrix, 
*     which also in the real code will be accesible explicitly
*     (or as W-like intermediate)
*       <mu_3| [H^AB,T^0_2] + [H^A,T^B_2] + [H^B,T^A_2] |HF>
*---------------------------------------------------------------------*
      KINT1SXY = KEND1
      KINT2SXY = KINT1SXY + NT1AMX*NVIRT*NVIRT
      KEND2    = KINT2SXY + NRHFT*NRHFT*NT1AMX

      KINT1SX = KEND2
      KINT2SX = KINT1SX + NT1AMX*NVIRT*NVIRT
      KEND2   = KINT2SX + NRHFT*NRHFT*NT1AMX
 
      KINT1SY = KEND2
      KINT2SY = KINT1SY + NT1AMX*NVIRT*NVIRT
      KEND2   = KINT2SY + NRHFT*NRHFT*NT1AMX

      LWRK2    = LWORK  - KEND2
      IF (LWRK2 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_T32INT_NODDY')
      ENDIF

      ! read T^0 doubles amplitudes from file and square up
      IOPT   = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND2))
      CALL CC_T2SQ(WORK(KEND2),WORK(KT2AMP0),ISYM0)


      CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SX),WORK(KINT2SX),
     &                       .FALSE.,DUMMY,DUMMY,
     &                       XLAMDP0,XLAMDH0,
     &                       WORK(KLAMPX),WORK(KLAMHX),
     &                       WORK(KEND2),LWRK2)

      CALL CCSDT_INTS1_NODDY(.TRUE.,WORK(KINT1SY),WORK(KINT2SY),
     &                       .FALSE.,DUMMY,DUMMY,
     &                       XLAMDP0,XLAMDH0,
     &                       WORK(KLAMPY),WORK(KLAMHY),
     &                       WORK(KEND2),LWRK2)

      CALL CCSDT_INTS2_NODDY(WORK(KINT1SXY),WORK(KINT2SXY),
     &                       XLAMDP0,XLAMDH0,
     &                       WORK(KLAMPX),WORK(KLAMHX),
     &                       WORK(KLAMPY),WORK(KLAMHY),
     &                       WORK(KEND2),LWRK2) 

      CALL DZERO(WORK(KTXY3B),NT1AMX*NT1AMX*NT1AMX) 

      ! note: this routine will overwrite WORK(KT2AMP0)
      CALL CCSDT_B3AM(WORK(KTXY3B),
     &                WORK(KINT1SXY),WORK(KINT2SXY),FOCKD,
     &                LISTX,IDLSTX,WORK(KINT1SX),WORK(KINT2SX),
     &                LISTY,IDLSTY,WORK(KINT1SY),WORK(KINT2SY),
     &                WORK(KSCR1),WORK(KT2AMP0),WORK(KEND2),LWRK2) 


      ! devide by orbital energy denominators and fix the sign
      ! such that it corresponds to that of T^XY in CCSDT_T32_NODDY
      CALL CCSDT_3AM(WORK(KTXY3B),FREQR,WORK(KSCR1),FOCKD,
     &               .FALSE.,DUMMY,.FALSE.,WORK(KEND1))
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KTXY3B),1)

      CALL CCSDT_CLEAN_T3(WORK(KTXY3B),NT1AMX,NVIRT,NRHFT)


*---------------------------------------------------------------------*
*     Precompute in a similar way the "easy" contributions to T^XY_3 
*     comming from the Jacobian matrix:
*        <mu_3| [[H,T^AB_1],T^0_2] + [H,T^AB_2] |HF>
*---------------------------------------------------------------------*
      KINT1S0  = KEND1
      KINT2S0  = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND2    = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KINT1SXY = KEND2
      KINT2SXY = KINT1SXY + NT1AMX*NVIRT*NVIRT
      KEND2    = KINT2SXY + NRHFT*NRHFT*NT1AMX

      LWRK2    = LWORK  - KEND2
      IF (LWRK2 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_T32INT_NODDY')
      ENDIF

      ! read T^0 doubles amplitudes from file and square up
      IOPT   = 2
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND2))
      CALL CC_T2SQ(WORK(KEND2),WORK(KT2AMP0),ISYM0)

      CALL CCSDT_READ_NODDY(.FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                      .FALSE.,DUMMY,DUMMY,
     &                      .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                      .FALSE.,DUMMY,DUMMY,
     &                      .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                      NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)

      CALL DZERO(WORK(KTXY3A),NT1AMX*NT1AMX*NT1AMX) 

      CALL CCSDT_A3AM(WORK(KTXY3A),T1AMXY,T2AMXY,
     &                DUMMY,ISYMR,WORK(KT2AMP0),T3AMP0,
     &                WORK(KINT1S0), WORK(KINT2S0),
     &                WORK(KINT1SXY),WORK(KINT2SXY),
     &                FOCKD,XLAMDP0,XLAMDH0,
     &                DUMMY,DUMMY,WORK(KSCR1),WORK(KEND2),LWRK2)


      ! devide by orbital energy denominators and fix the sign
      ! such that it corresponds to that of T^XY in CCSDT_T32_NODDY
      CALL CCSDT_3AM(WORK(KTXY3A),FREQR,WORK(KSCR1),FOCKD,
     &               .FALSE.,DUMMY,.FALSE.,WORK(KEND1))
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KTXY3A),1)

      CALL CCSDT_CLEAN_T3(WORK(KTXY3A),NT1AMX,NVIRT,NRHFT)
     

*---------------------------------------------------------------------*
*     precalculate modified ^{XY}W intermediates 
*---------------------------------------------------------------------*
      KTHEOCC = KEND1
      KWXYO   = KTHEOCC + NT1AMX*NRHFT*NRHFT
      KWXYV   = KWXYO   + NT1AMX*NRHFT*NRHFT
      KEND2   = KWXYV   + NT1AMX*NRHFT*NRHFT

      LWRK2    = LWORK  - KEND2
      IF (LWRK2.LT.NT2AMX)
     &    CALL QUIT('Out of memory in CCSDT_T32INT_NODDY')

      ! read T^X doubles amplitudes from file and square up
      IOPT = 2
      Call CC_RDRSP(LISTX,IDLSTX,ISYMX,IOPT,MODEL,DUMMY,WORK(KEND2))
      Call CCLR_DIASCL(WORK(KEND2),2.0D0,ISYMX)  
      CALL CC_T2SQ(WORK(KEND2),WORK(KT2AMPX),ISYMX)

      ! read T^X doubles amplitudes from file and square up
      IOPT = 2
      Call CC_RDRSP(LISTY,IDLSTY,ISYMY,IOPT,MODEL,DUMMY,WORK(KEND2))
      Call CCLR_DIASCL(WORK(KEND2),2.0D0,ISYMY)  
      CALL CC_T2SQ(WORK(KEND2),WORK(KT2AMPY),ISYMY)

      IF (HAVETB3T3) THEN
         CALL DZERO(TAMXY3,NT1AMX*NT1AMX*NT1AMX)
      END IF

      CALL CCSDT_O32VIR_NODDY(WXINT,XPROP,WORK(KT2AMPX),HAVEX,
     &                        WYINT,YPROP,WORK(KT2AMPY),HAVEY,
     &                        T3AMP0,WORK(KT2AMP0),
     &                        XYMAT,HAVEXY,
     &                        WORK(KTHEOCC),WORK(KWXYV),WORK(KWXYO),
     &                        TAMXY3,HAVETB3T3,
     &                        WXYOINT,.TRUE.,WXYVINT,.TRUE.)



      CALL CCSDT_CLEAN_T3(WXYOINT,NT1AMX,NVIRT,NRHFT)
      CALL CCSDT_CLEAN_T3(WXYVINT,NT1AMX,NVIRT,NRHFT)

c     -------------------------------------------------------------
c     add some more contributions to WXYVINT to obtain full
c            theta(e-bar m-bar, f n-bar, d l-bar)
c     -------------------------------------------------------------
      CALL CCSDT_WXYVINT_NODDY(WXYVINT,THETAY,THETAX,YPROP,XPROP,
     &                         NVIRT,NRHFT,NORBT)
      CALL CCSDT_CLEAN_T3(WXYVINT,NT1AMX,NVIRT,NRHFT)


c     -------------------------------------------------------------
c     devide by orbital energy denominators and fix the sign
c     such that it corresponds to that of T^XY in CCSDT_T32_NODDY
c     -------------------------------------------------------------
      CALL CCSDT_3AM(WXYOINT,FREQR,WORK(KSCR1),FOCKD,
     &               .FALSE.,DUMMY,.FALSE.,WORK(KEND1))
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WXYOINT,1)

      CALL CCSDT_3AM(WXYVINT,FREQR,WORK(KSCR1),FOCKD,
     &               .FALSE.,DUMMY,.FALSE.,WORK(KEND1))
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WXYVINT,1)


c     ----------------------------------------------------------------
c     add to WXYOINT 1/3 of the contribution from the B T^X T^Y and 
c     A T^XY transformation such, that we can later generate from 
c     this intermediate the terms with double-bars on occupied indeces
c     ----------------------------------------------------------------
      CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,1.0D0/3.0D0,
     &           WORK(KTXY3B),1,WXYOINT,1)
      CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,1.0D0/3.0D0,
     &           WORK(KTXY3A),1,WXYOINT,1)


*---------------------------------------------------------------------*
*     if the HAVETB3T3 test option is set, complete the T^XY_3 vector:
*---------------------------------------------------------------------*
      IF (HAVETB3T3) THEN
         KTHVIRT = KEND1
         KWXYV   = KTHVIRT + NT1AMX*NVIRT*NVIRT
         KWXYO   = KWXYV   + NT1AMX*NVIRT*NVIRT
         KEND2   = KWXYO   + NT1AMX*NVIRT*NVIRT

         LWRK2    = LWORK  - KEND2
         IF (LWRK2.LT.0)
     &       CALL QUIT('Out of memory in CCSDT_T32INT_NODDY')

         CALL CCSDT_O32OCC_NODDY(THETAX,XPROP,HAVEX,
     &                           THETAY,YPROP,HAVEY,
     &                           TAMXY3,WORK(KTHVIRT),
     &                           WORK(KWXYV),WORK(KWXYO),.TRUE.,.TRUE.,
     &                           DUMMY,.FALSE.,DUMMY,.FALSE.)
         CALL CCSDT_CLEAN_T3(TAMXY3,NT1AMX,NVIRT,NRHFT)


         ! devide by orbital energy denominators and fix the sign
         ! such that it corresponds to that of T^XY in CCSDT_T32_NODDY
         CALL CCSDT_3AM(TAMXY3,FREQR,WORK(KSCR1),FOCKD,
     &                  .FALSE.,DUMMY,.FALSE.,WORK(KEND1))
         CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,TAMXY3,1)


         ! add contributions from B T^X T^Y transformation
         CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,1.0D0,WORK(KTXY3B),1,TAMXY3,1)


         ! add contributions from A T^XY transformation
         CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,1.0D0,WORK(KTXY3A),1,TAMXY3,1)


C        WRITE(LUPRI,*) 'CCSDT_T32INT_NODDY> T^XY:',LISTR,IDLSTR,FREQR
C        CALL CCSDT_CLEAN_T3(TAMXY3,NT1AMX,NVIRT,NRHFT)
C        CALL PRINT_PT3_NODDY(TAMXY3)
      END IF

      CALL QEXIT('CCT32INT')
      RETURN
      END
*=====================================================================*
*---------------------------------------------------------------------*
C  /* Deck ccsdt_gwbic */
*=====================================================================*
      SUBROUTINE CCSDT_GWBIC(LISTL,IDLSTL,WBAR,THETA,SPLIT_WT,
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

      LOGICAL LORX, SPLIT_WT
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
      REAL THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      DOUBLE PRECISION FREQ, HALF
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION WBAR(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
      DOUBLE PRECISION THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT)
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
c
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
        CALL QUIT('Unknown list in CCSDT_GWBIC')
      END IF

*---------------------------------------------------------------------*
* open files:
*---------------------------------------------------------------------*
      LU3VI2  = -1
      LU3FOP  = -1
      LU3FOP2 = -1

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
      CALL CCPRPAO(LABELY,.TRUE.,WORK(KFOCKY),IRREP,ISYM,IERR,
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
      LUSIFC = -1
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

                  IF (SPLIT_WT) THEN
C                  ---------------------------------------------------
C                  1)Store this contribution of WMAT in THETA(ai,bj,dl):
C                    but first divide by orbital energies and clean up
C                  2)Store zeroth-order T3 amplitudes in T0AMP:
C                  ---------------------------------------------------
                   CALL WBD_DIA(B,ISYMB,D,ISYMD,-FREQ,
     *                          ISWMAT,WORK(KWMAT),
     *                          WORK(KDIAGW),WORK(KFOCKD))
                   CALL T3_FORBIDDEN(WORK(KWMAT),ISYMY,ISYMB,B,ISYMD,D)

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
     *                       + ISAIK(ISYMAI, ISYML)+NT1AM(ISYMAI)*(L-1)
                  
                       CALL DCOPY(NT1AMX,WORK(KWMAT+KOFF1),1,
     *                                   THETA(1,1,B,J,D,L),1)
                      END DO
                     END DO
                    END DO
                   END DO
                   ! Reinitialize WMAT to avoid a double count:
                   CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
                  END IF
 
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
     &                       HAVET3Y,TAMP3Y,PRINT_T3,
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
      LOGICAL LOCDBG, TOT_T3Y, HAVET3Y, PRINT_T3
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
      REAL TAMP3Y(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), FOCK0(*)
      DOUBLE PRECISION WORK(LWORK), FREQY, DDOT
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      DOUBLE PRECISION WINT(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
      DOUBLE PRECISION THETA(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
      DOUBLE PRECISION T0AMP(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
      DOUBLE PRECISION TAMP3Y(NVIRT,NRHFT,NVIRT,NRHFT,NVIRT,NRHFT) 
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
c
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
         LU3SRTR  = -1
         LUCKJDR  = -1
         LUDELDR  = -1
         LUDKBCR  = -1
         LUDKBCR4 = -1
 
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
      LUSIFC = -1
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
        CALL CCPRPAO(LABELY,.TRUE.,WORK(KFOCKY),IRREP,IMAT,IERR,
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
C
                  CALL WBD_T2(.FALSE.,B,ISYMB,D,ISYMD,
     *                        WORK(KT2TP),ISYM0,WORK(KT2TP),ISYM0,
     *                        WORK(KFOCKY),ISTAMY,WORK(KINDSQW),LENSQW,
     *                        WORK(KWMAT),ISWMAT,WORK(KEND4),LWRK4)
C
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

      IF (PRINT_T3) THEN
        WRITE(LUPRI,*) 'CCSDT_GWTIC> W^Y:'
        WRITE(LUPRI,*) 'CCSDT_GWTIC> LISTR,IDLSTR:',LISTR,IDLSTR
        CALL PRINT_PT3_NODDY(WINT)
        WRITE(LUPRI,*) 'CCSDT_GWTIC> Theta^Y:'
        WRITE(LUPRI,*) 'CCSDT_GWTIC> LISTR,IDLSTR:',LISTR,IDLSTR
        CALL PRINT_PT3_NODDY(THETA)
      END IF

      IF (HAVET3Y) THEN
       ! accumulate T^Y amplitudes in TAMP3Y and print
       DO K = 1, NRHFT
       DO C = 1, NVIRT
         DO J = 1, NRHFT
         DO B = 1, NVIRT
           DO I = 1, NRHFT
           DO A = 1, NVIRT
             TAMP3Y(A,I,B,J,C,K) =  
     &          WINT(A,I,B,J,C,K) + THETA(A,I,B,J,C,K) +
     &          WINT(C,K,A,I,B,J) + THETA(C,K,A,I,B,J) + 
     &          WINT(B,J,C,K,A,I) + THETA(B,J,C,K,A,I) 
           END DO
           END DO
         END DO
         END DO
       END DO
       END DO
       IF (PRINT_T3) THEN
         WRITE(LUPRI,*) 'CCSDT_GWTIC> T3^Y:'
         WRITE(LUPRI,*) 'CCSDT_GWTIC> LISTR,IDLSTR:',LISTR,IDLSTR
         CALL PRINT_PT3_NODDY(TAMP3Y)
       END IF
      END IF

      CALL QEXIT('CCSDT_GWTIC')
 
      RETURN
      END
*=====================================================================*
