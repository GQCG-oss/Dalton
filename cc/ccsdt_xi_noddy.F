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
*=====================================================================*
      SUBROUTINE CCSDT_XI_NODDY(IOPTRES,XI1,XI2,XI1EFF,XI2EFF,
     &                          FOCKB,FREQ,
     &                          XLAMDP,XLAMDH,CMO,FOCK0,
     &                          IDOTS,DOTPROD,LISTDP,ITRAN,
     &                          NXTRAN,MXVEC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to XKSI vector
*
*             XKSI^eff_1,2 = XSKI_1,2(CCSD) + XKSI_1,2(T3)
*                            - A_1,2;3 (w_3 - w)^1 XKSI_3
*
*        
*     Written by Christof Haettig, April 2002, based on CCSDT_TRIPLE.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.true.)

      CHARACTER LISTDP*3
      INTEGER LWORK, IOPTRES, ITRAN, MXVEC, NXTRAN
      INTEGER IDOTS(MXVEC,NXTRAN)

#if defined (SYS_CRAY)
      REAL DOTPROD(MXVEC,NXTRAN)
      REAL WORK(LWORK), SIXTH, FREQC
      REAL XI1(NT1AMX), XI2(NT2AMX), FREQ, DUMMY, TCON
      REAL XI1EFF(NT1AMX), XI2EFF(NT2AMX), DDOT
      REAL FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NXTRAN)
      DOUBLE PRECISION WORK(LWORK), SIXTH, FREQC
      DOUBLE PRECISION XI1(NT1AMX), XI2(NT2AMX), FREQ, DUMMY
      DOUBLE PRECISION XI1EFF(NT1AMX), XI2EFF(NT2AMX), DDOT, TCON
      DOUBLE PRECISION FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
#endif
      PARAMETER( SIXTH = 1.0D0/6.0D0 )

      CHARACTER*10 MODEL
      INTEGER KEND1, LWRK1, KFOCKD, KIAJB, KYIAJB, KINT1T, KINT2T, 
     &        KXKSI1, KXKSI2, KXKSI3, KINT1S, KINT2S, KLC3AM, KSCR1,
     &        IVEC, IDLSTC, ISYMC, KEND2, LWRK2
      INTEGER IJ, NIJ, LUSIFC, INDEX, IDUMMY, ILSTSYM

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_XI_NODDY')

      IF (DIRECT) CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_XI_NODDY')

      IF (LOCDBG) THEN
        write(lupri,*)  'ccsdt_xi_noddy> frequency:',freq
        write(lupri,*)  'ccsdt_xi_noddy> xksi vector on input:'
        ! call cc_prp(xi1,xi2,1,1,1)
        write(lupri,*)  'norm^2 = ',ddot(nt1amx+nt2amx,xi1,1,xi1,1)
        write(lupri,*)  'ccsdt_xi_noddy> norm^2(FOCKB):',
     &      ddot(norbt*norbt,fockb,1,fockb,1)
        write(lupri,*)  'ccsdt_xi_noddy> norm^2(FOCK0):',
     &      ddot(norbt*norbt,fock0,1,fock0,1)
        WRITE (LUPRI,*) 'FOCKB:'
        CALL OUTPUT(FOCKB,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KFOCKD = 1
      KINT1T = KFOCKD + NORBT
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KIAJB  = KINT2T + NRHFT*NRHFT*NT1AMX
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KINT1S = KYIAJB + NT1AMX*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KXKSI1 = KINT2S + NRHFT*NRHFT*NT1AMX
      KXKSI2 = KXKSI1 + NT1AMX
      KXKSI3 = KXKSI2 + NT1AMX*NT1AMX
      KEND1  = KXKSI3 + NT1AMX*NT1AMX*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_XI_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     Compute xksi_3 and add xksi_1,2(T3) contribution:
*---------------------------------------------------------------------*
      CALL CCSDT_XI3_NODDY(WORK(KXKSI3),FOCKB,FREQ,.TRUE.,XI2,
     &                     WORK(KINT1S),WORK(KINT2S),
     &                    .TRUE.,WORK(KINT1T),WORK(KINT2T),
     &                    .TRUE.,WORK(KIAJB), WORK(KYIAJB),
     &                    XLAMDP,XLAMDH,CMO,FOCK0,WORK(KFOCKD),
     &                    WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Now we split:
*       for IOPTRES < 5 we compute the effective Xi vector
*       for IOPTRES = 5 we compute the contractions Tbar^C Xi
*---------------------------------------------------------------------*
      IF (IOPTRES.GE.1 .AND. IOPTRES.LE.4) THEN

C        -----------------------------------------------------------
C        Calculate -X_nu_3 = Xi_nu_3 / (eps_nu_3 - w)
C        -----------------------------------------------------------
         IF (LWRK1.LT.NT1AMX) CALL QUIT('OUT OF MEMORY - CCSDT_XI3')
         CALL CCSDT_3AM(WORK(KXKSI3),FREQ,WORK(KEND1),WORK(KFOCKD))

C        --------------------------------------------------------
C        Add contribution of xksi_3 to the effective xksi vector:
C        --------------------------------------------------------
         CALL DZERO(WORK(KXKSI2),NT1AMX*NT1AMX)
         CALL DZERO(WORK(KXKSI1),NT1AMX)

         CALL CCSDT_OMEGA1(WORK(KXKSI1),WORK(KIAJB),WORK(KXKSI3))
 
         CALL CCSDT_OMEGA2(WORK(KXKSI2),WORK(KINT1T),WORK(KINT2T),
     *                     WORK(KXKSI3),FOCK0)
 
         DO I = 1,NT1AMX
            XI1EFF(I) = XI1(I) + WORK(KXKSI1+I-1)
         END DO
 
         DO I = 1,NT1AMX
            DO J = 1,I
               IJ = NT1AMX*(I-1) + J
               NIJ = INDEX(I,J)
               XI2EFF(NIJ) = XI2(NIJ) + WORK(KXKSI2+IJ-1)
            END DO
         END DO

         IF (LOCDBG) THEN
           write(lupri,*)  'ccsdt_xi_noddy> xksi vector on output:'
           write(lupri,*)  'norm^2 = ',ddot(nt1amx+nt2amx,xi1,1,xi1,1)
           ! call cc_prp(xi1,xi2,1,1,1)
           write(lupri,*)  'ccsdt_xi_noddy> effective xksi vector:'
           write(lupri,*)  'norm^2 = ',
     &                           ddot(nt1amx+nt2amx,xi1eff,1,xi1eff,1)
           ! call cc_prp(xi1eff,xi2eff,1,1,1)
         END IF

      ELSE IF (IOPTRES.EQ.5) THEN
*---------------------------------------------------------------------*
*       for IOPTRES = 5 we compute the contractions L^C Xi :
*       here in the triples code we compute only the contribution 
*       from the triples result vector: 
*       (Note, that we have -Xksi_3 in memory!!)
*---------------------------------------------------------------------*
        KLC3AM = KEND1 
        KSCR1  = KLC3AM + NT1AMX*NT1AMX*NT1AMX
        KEND2  = KSCR1  + NT1AMX

        LWRK2  = LWORK - KEND2
        IF (LWRK2 .LT. NT2AMX) THEN
           CALL QUIT('Insufficient space in CCSDT_XI_NODDY')
        ENDIF

        IVEC = 1
        DO WHILE (IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

          IDLSTC = IDOTS(IVEC,ITRAN)
          ISYMC  = ILSTSYM(LISTDP,IDLSTC)

          write(lupri,*) 'norm^2(xksi3) - before TBAR31:',
     &     ddot(nt1amx*nt1amx*nt1amx,work(kxksi3),1,work(kxksi3),1)

          IF (LISTDP(1:3).EQ.'L1 ' .OR. LISTDP(1:3).EQ.'LE ') THEN

             CALL CCSDT_TBAR31_NODDY(WORK(KLC3AM),FREQC,LISTDP,IDLSTC,
     &                               XLAMDP,XLAMDH,
     &                               FOCK0,WORK(KFOCKD),
     &                               WORK(KSCR1),
     &                               WORK(KIAJB),WORK(KINT1T),
     &                               WORK(KINT2T),WORK(KEND2),LWRK2)
          ELSE
            CALL QUIT('Unknown or illegal list in CCSDT_XI_NODDY.')
          END IF

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'CCSDT_XI_NODDY> IVEC,ITRAN:',IVEC,ITRAN
            WRITE(LUPRI,*) 'NORM^2(XKSI3):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KXKSI3),1,WORK(KXKSI3),1)
            WRITE(LUPRI,*) 'NORM^2(TC3AM):',
     &        DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KLC3AM),1,WORK(KLC3AM),1)
          END IF

          TCON = - SIXTH *
     &      DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KLC3AM),1,WORK(KXKSI3),1)

          DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN) + TCON

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'CCSDT_XI_NODDY:'
            WRITE(LUPRI,*) 'L_3  x Xi:3:',TCON
          END IF


          IVEC = IVEC + 1

        END DO


      ELSE
        CALL QUIT('Illegal value for IOPTRES IN CCSDT_XI_NODDY')
      END IF

      CALL QEXIT('CCSDT_XI_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XI_NODDY                       *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XI3_NODDY(XI3,FOCKB,FREQ,ADD_DOUBLE,XI2,
     &                           XINT1S,XINT2S,
     &                           DO_INT1,XINT1T,XINT2T,
     &                           DO_INT2,XIAJB,YIAJB,
     &                           XLAMDP,XLAMDH,CMO,FOCK0,FOCKD,
     &                           WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples part of XKSI vector and the integrals
*             XINT1S = (ck|bd) and XINT2S = (ck|lj)
*
*             if (do_int1) compute also
*                 XINT1T = (kc|bd) and 
*                 XINT2T = (kc|lj)
*
*             if (do_int2) compute also
*                 XIAJB  = 2(kc|ld) - (kd|lc)
*                 YIAJB  =  (kc|ld)
*
*             if (add_double) compute also
*                XKSI_1,2(CC3) = XSKI_1,2(CCSD) + XKSI_1,2(T3)
*        
*     Written by Christof Haettig, April 2002, based on CCSDT_TRIPLE.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.true.)

      LOGICAL ADD_DOUBLE, DO_INT1, DO_INT2
      INTEGER LWORK
#if defined (SYS_CRAY)
      REAL WORK(LWORK), FOCKD(*)
      REAL XI2(NT2AMX), FREQ, DUMMY, DDOT
      REAL XI3(NT1AMX,NT1AMX,NT1AMX)
      REAL XINT1T(NT1AMX*NVIRT*NVIRT)
      REAL XINT2T(NRHFT*NRHFT*NT1AMX)
      REAL XINT1S(NT1AMX*NVIRT*NVIRT)
      REAL XINT2S(NRHFT*NRHFT*NT1AMX)
      REAL XIAJB(NT1AMX*NT1AMX)
      REAL YIAJB(NT1AMX*NT1AMX)
      REAL FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
#else
      DOUBLE PRECISION WORK(LWORK), FOCKD(*)
      DOUBLE PRECISION XI2(NT2AMX), FREQ, DUMMY, DDOT
      DOUBLE PRECISION XI3(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1T(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2T(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XINT1S(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2S(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION YIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
#endif

      CHARACTER*10 MODEL
      INTEGER KXINT, KEND1, LWRK1, KEND2, LWRK2, KT3AM,
     &        KSCR1, KRECNR, KT2AM, KXKSI2, IOPT
      INTEGER ISYMD, ILLL, IDEL, ISYDIS, IDUMMY, IJ, NIJ, LUSIFC, INDEX

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      CALL QENTER('CCSDT_XI3_NODDY')

      IF (DIRECT) CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_XI_NODDY')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1  = 1
      KT3AM  = KSCR1  + NT1AMX
      KT2AM  = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KEND1  = KT2AM  + NT1AMX*NT1AMX

      IF (ADD_DOUBLE) THEN
        KXKSI2 = KEND1
        KEND1  = KXKSI2 + NT1AMX*NT1AMX
      END IF

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_XI3_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Loop over distributions of integrals:
*---------------------------------------------------------------------*
      CALL DZERO(XINT1S,NT1AMX*NVIRT*NVIRT)
      CALL DZERO(XINT2S,NT1AMX*NRHFT*NRHFT)

      IF (DO_INT1) THEN
         CALL DZERO(XINT1T,NT1AMX*NVIRT*NVIRT)
         CALL DZERO(XINT2T,NT1AMX*NRHFT*NRHFT)
      END IF
      IF (DO_INT2) THEN
         CALL DZERO(XIAJB,NT1AMX*NT1AMX)
         CALL DZERO(YIAJB,NT1AMX*NT1AMX)
      END IF

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_XI3_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  WORK(KRECNR),DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            IF (DO_INT1) THEN
               CALL CCSDT_TRAN1(XINT1T,XINT2T,XLAMDP,XLAMDH,
     *                          WORK(KXINT),IDEL)
            END IF

            IF (DO_INT2) THEN
               CALL CC3_TRAN2(XIAJB,YIAJB,XLAMDP,XLAMDH,
     *                        WORK(KXINT),IDEL)
            END IF
 
            CALL CCSDT_TRAN3(XINT1S,XINT2S,XLAMDP,XLAMDH,
     *                       WORK(KXINT),IDEL)
 
         END DO   
      END DO  
*---------------------------------------------------------------------*
*     End Loop over distributions of integrals.
*---------------------------------------------------------------------*
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMOP,IOPT,MODEL,DUMMY,WORK(KT3AM))
      CALL CC_T2SQ(WORK(KT3AM),WORK(KT2AM),ISYMOP)

C     --------------------------------------------
C     Calculate the zero-order triples amplitudes:
C     --------------------------------------------
      IF ((NONHF) .AND. (NFIELD .GT. 0)) THEN
        CALL QUIT('FINITE FIELD NOT IMPLEMENTED IN CCSDT_XI3_NODDY')
      END IF
 
      CALL DZERO(WORK(KT3AM),NT1AMX*NT1AMX*NT1AMX) 
      CALL CCSDT_T3AM(WORK(KT3AM),XINT1S,XINT2S,WORK(KT2AM),
     *                WORK(KSCR1),FOCKD,DUMMY,DUMMY)

C     ------------------------------------------------------
C     Calculate triples contribution to doubles xksi vector:
C      xksi_mu_2 = "xksi_mu_2(CCSD)" + <mu_2|[V,T^0_3]|HF>
C     ------------------------------------------------------
      IF (ADD_DOUBLE) THEN
        CALL DZERO(WORK(KXKSI2),NT1AMX*NT1AMX)
        CALL CCSDT_XKSI2_2(WORK(KXKSI2),FOCKB,WORK(KT3AM))
        write(lupri,*)  'ccsdt_xi_noddy> norm of first contrib:',
     &   ddot(nt1amx*nt1amx,work(kxksi2),1,work(kxksi2),1)

        write(lupri,*)  'ccsdt_xi_noddy> norm of xi2 before:',
     &   ddot(nt2amx,xi2,1,xi2,1)
        DO I = 1,NT1AMX
          DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            XI2(NIJ) = XI2(NIJ) + WORK(KXKSI2+IJ-1)
          END DO
        END DO 
C       write(lupri,*), 'ccsdt_xi_noddy> norm of xi2 after contrib:',
C    &   ddot(nt2amx,xi2,1,xi2,1)
      END IF

C     ---------------------------------------------------------
C     Calculate the triples part of the xksi vector (divided by
C     the orbital energy differences and the frequency) and its
C     contribution to the effective xksi vector:
C     ---------------------------------------------------------
      CALL DZERO(XI3,NT1AMX*NT1AMX*NT1AMX)

      CALL CCSDT_XKSI3(XI3,FOCKB,WORK(KT3AM),WORK(KT2AM))
 
      CALL QEXIT('CCSDT_XI3_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XI_NODDY                       *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_T31_NODDY(TB3AM,LISTB,IDLSTB,FREQB,
     &                           XINT1S,XINT2S,
     &                           DO_INT1,XINT1T,XINT2T,
     &                           DO_INT2,XIAJB,YIAJB,
     &                           XLAMPB,XLAMHB,FOCKB,
     &                           XLAMDP,XLAMDH,FOCK0,
     &                           CMO,FOCKD,
     &                           WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples part of first-order response amplitudes
*             and the property integrals as well as the integrals 
*              XINT1S = (ck|bd)-bar and XINT2S = (ck|lj)-bar
*
*             if (do_int1) compute also
*                 XINT1T = (kc|bd) and 
*                 XINT2T = (kc|lj)
*
*             if (do_int2) compute also
*                 XIAJB  = 2(kc|ld) - (kd|lc)
*                 YIAJB  =  (kc|ld)
*
*   Written by Christof Haettig, April 2002, based on CCSDT_XI3_NODDY.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccr1rsp.h>

      LOGICAL LOCDBG, PRINT_T3
      PARAMETER (LOCDBG=.true.,PRINT_T3=.false.)

      LOGICAL ADD_DOUBLE, DO_INT1, DO_INT2

      CHARACTER LISTB*3
      INTEGER LWORK, IDLSTB
#if defined (SYS_CRAY)
      REAL WORK(LWORK), FOCKD(*)
      REAL FREQB, DUMMY, DDOT
      REAL TB3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL XINT1T(NT1AMX*NVIRT*NVIRT)
      REAL XINT2T(NRHFT*NRHFT*NT1AMX)
      REAL XINT1S(NT1AMX*NVIRT*NVIRT)
      REAL XINT2S(NRHFT*NRHFT*NT1AMX)
      REAL XIAJB(NT1AMX*NT1AMX)
      REAL YIAJB(NT1AMX*NT1AMX)
      REAL FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
      REAL XLAMPB(*), XLAMHB(*), TWO
#else
      DOUBLE PRECISION WORK(LWORK), FOCKD(*)
      DOUBLE PRECISION FREQB, DUMMY, DDOT
      DOUBLE PRECISION TB3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1T(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2T(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XINT1S(NT1AMX*NVIRT*NVIRT)
      DOUBLE PRECISION XINT2S(NRHFT*NRHFT*NT1AMX)
      DOUBLE PRECISION XIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION YIAJB(NT1AMX*NT1AMX)
      DOUBLE PRECISION FOCKB(*), CMO(*), XLAMDP(*), XLAMDH(*), FOCK0(*)
      DOUBLE PRECISION XLAMPB(*), XLAMHB(*), TWO
#endif
      PARAMETER (TWO = 2.0D0)

      CHARACTER*10 MODEL, LABELB*8
      INTEGER KTB1AM, KEND0, KSCR1, KTB2AM, KEND1, NAI, NBJ, NCK,
     &        LWRK1, IOPT, KXINT, ISYMD, ILLL, IDEL, KEND2, LWRK2,
     &        KT2AM, ISYDIS, ISYMB, KRECNR, IERR, IRREP, NAIBJCK

*---------------------------------------------------------------------*
*     Do some initial tests and get property integrals
*---------------------------------------------------------------------*
      if (NSYM.NE.1) CALL QUIT('NO SYMMETRY FOR CCSDT_T31_NODDY')

      IF (LISTB(1:3).NE.'R1 ') 
     &   CALL QUIT('LIST MUST BE R1 FOR CCSDT_T31_NODDY')

      ! get symmetry, frequency and integral label from common blocks
      ! defined in ccr1rsp.h
      ISYMB  = ISYLRT(IDLSTB)
      FREQB  = FRQLRT(IDLSTB)
      LABELB = LRTLBL(IDLSTB)

      IF (LORXLRT(IDLSTB)) 
     &  CALL QUIT('NO ORBITAL RELAXED PERTURBATION IN CCSDT_T31_NODDY')

      ! read property integrals from file:
      CALL CCPRPAO(LABELB,FOCKB,IRREP,ISYMB,IERR,WORK,LWORK)
      IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMB)) THEN
        CALL QUIT('CCSDT_T31_NODDY: error reading operator '//LABELB)
      ELSE IF (IERR.LT.0) THEN
        CALL DZERO(FOCKB,N2BST(ISYMB))
      END IF
 
      ! transform property integrals to Lambda-MO basis
      CALL CC_FCKMO(FOCKB,XLAMDP,XLAMDH,WORK,LWORK,ISYMB,1,1)

*---------------------------------------------------------------------*
*     Compute xksi_3/eps_3 contribution to T^B_3 and some integrals:
*---------------------------------------------------------------------*
      CALL CCSDT_XI3_NODDY(TB3AM,FOCKB,FREQB,.FALSE.,DUMMY,
     &                           XINT1S,XINT2S,
     &                    DO_INT1,XINT1T,XINT2T,
     &                    DO_INT2,XIAJB, YIAJB,
     &                    XLAMDP,XLAMDH,CMO,FOCK0,FOCKD,WORK,LWORK)

      IF (LWORK.LT.NT1AMX) CALL QUIT('OUT OF MEMORY - CCSDT_TI3')

      CALL CCSDT_3AM(TB3AM,FREQB,WORK,FOCKD)

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'CCSDT_T31_NODDY:  FREQB =',FREQB
        WRITE(LUPRI,*) 'NORM^2(TB3AM/XKI3)',
     &   DDOT(NT1AMX*NT1AMX*NT1AMX,TB3AM,1,TB3AM,1)
        WRITE(LUPRI,*) 'NORM^2(FOCK0):',
     &   DDOT(NORBT*NORBT,FOCK0,1,FOCK0,1)
        WRITE(LUPRI,*) 'NORM^2(FOCKB):',
     &   DDOT(NORBT*NORBT,FOCKB,1,FOCKB,1)
        WRITE(LUPRI,*) 'NORM^2(XINT1S):',
     &   DDOT(NT1AMX*NVIRT*NVIRT,XINT1S,1,XINT1S,1)
        WRITE(LUPRI,*) 'NORM^2(XINT2S):',
     &   DDOT(NT1AMX*NRHFT*NRHFT,XINT2S,1,XINT2S,1)
        IF (DO_INT1) WRITE(LUPRI,*) 'NORM^2(XINT2T):',
     &   DDOT(NT1AMX*NVIRT*NVIRT,XINT1T,1,XINT1T,1)
        IF (DO_INT1) WRITE(LUPRI,*) 'NORM^2(XINT2T):',
     &   DDOT(NT1AMX*NRHFT*NRHFT,XINT2T,1,XINT2T,1)
        IF (DO_INT2) WRITE(LUPRI,*) 'NORM^2(XIAJB):',
     &   DDOT(NT1AMX*NT1AMX,XIAJB,1,XIAJB,1)
        IF (DO_INT2) WRITE(LUPRI,*) 'NORM^2(YIAJB):',
     &   DDOT(NT1AMX*NT1AMX,YIAJB,1,YIAJB,1)
        WRITE(LUPRI,*) 'NORM^2(XLAMDP):',
     &   DDOT(NLAMDT,XLAMDP,1,XLAMDP,1)
        WRITE(LUPRI,*) 'NORM^2(XLAMDH):',
     &   DDOT(NLAMDT,XLAMDH,1,XLAMDH,1)

C       WRITE(LUPRI,*) 'TB3AM after CCSDT_XI3_NODDY:'
C       CALL OUTPUT(TB3AM,1,NT1AMX*NT1AMX,1,NT1AMX,
C    &                      NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
      END IF

*---------------------------------------------------------------------*
*     Some allocations:
*---------------------------------------------------------------------*
      KTB1AM  = 1
      KEND0   = KTB1AM  + NT1AMX

*---------------------------------------------------------------------*
*     Compute contribution from <mu_3|[H,T^B_2]|HF>/eps_3:
*---------------------------------------------------------------------*
      KSCR1  = KEND0
      KTB2AM = KSCR1  + NT1AMX
      KEND1  = KTB2AM + NT1AMX*NT1AMX
      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_T31_NODDY')
      ENDIF

      IOPT = 3 ! read singles and doubles amplitudes
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KTB1AM),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMB)
      CALL CC_T2SQ(WORK(KEND1),WORK(KTB2AM),ISYMB)

      CALL CCSDT_T3AM_R(TB3AM,FREQB,XINT1S,XINT2S,WORK(KTB2AM),
     &                  WORK(KSCR1),FOCKD)

C       WRITE(LUPRI,*) 'TB3AM after CCSDT_T3AM_R-1:'
C       CALL OUTPUT(TB3AM,1,NT1AMX*NT1AMX,1,NT1AMX,
C    &                      NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
*---------------------------------------------------------------------*
*     Calculate response Lambda matrices:
*---------------------------------------------------------------------*
      CALL CCLR_LAMTRA(XLAMDP,XLAMPB,XLAMDH,XLAMHB,WORK(KTB1AM),ISYMB)

*---------------------------------------------------------------------*
*     Loop over distributions of integrals and compute first-order
*     response versions XINT1S=(ck|db)-bar and XINT2S=(ck|lj)-bar
*---------------------------------------------------------------------*
      CALL DZERO(XINT1S,NT1AMX*NVIRT*NVIRT)
      CALL DZERO(XINT2S,NT1AMX*NRHFT*NRHFT)

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND0
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_T31_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  WORK(KRECNR),DIRECT)
 
C           ----------------------------------------------------------
C           Calculate integrals needed in for remaining contributions:
C           ----------------------------------------------------------

            ! XINT1S = XINT1S + (C-bar K|B D)
            ! XINT2S = XINT2S + (C-bar K|L J)
            CALL CCSDT_TRAN3_R(XINT1S,XINT2S,XLAMDP,XLAMDH,
     &                         XLAMPB,XLAMDH,XLAMDP,XLAMDH,
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C K-bar|B D)
            ! XINT2S = XINT2S + (C K-bar|L J)
            CALL CCSDT_TRAN3_R(XINT1S,XINT2S,XLAMDP,XLAMDH,
     &                         XLAMDP,XLAMHB,XLAMDP,XLAMDH,
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C K|B-bar D)
            ! XINT2S = XINT2S + (C K|L J-bar)
            CALL CCSDT_TRAN3_R(XINT1S,XINT2S,XLAMDP,XLAMDH,
     &                         XLAMDP,XLAMDH,XLAMPB,XLAMHB,
     &                         WORK(KXINT),IDEL)
 
         END DO   
      END DO  

*---------------------------------------------------------------------*
*     Compute contribution from <mu_3|[H^B,T^0_2]|HF>/eps_3:
*---------------------------------------------------------------------*
      KSCR1  = 1
      KT2AM  = KSCR1 + NT1AMX
      KEND1  = KT2AM + NT1AMX*NT1AMX
      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_T31_NODDY')
      ENDIF

      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMOP,IOPT,MODEL,DUMMY,WORK(KEND1))
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMOP)

      CALL CCSDT_T3AM_R(TB3AM,FREQB,XINT1S,XINT2S,WORK(KT2AM),
     &                  WORK(KSCR1),FOCKD)

C       WRITE(LUPRI,*) 'TB3AM after CCSDT_T3AM_R-2:'
C       CALL OUTPUT(TB3AM,1,NT1AMX*NT1AMX,1,NT1AMX,
C    &                      NT1AMX*NT1AMX,NT1AMX,1,LUPRI)
*---------------------------------------------------------------------*
*     Now we need to turn the sign, since
*      T_3 = - { Xksi_3 + <mu_3|[H^B,T^0_2]|HF>
*                       + <mu_3|[H,T^B_2]|HF>   } / eps_3
*---------------------------------------------------------------------*
      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,TB3AM,1)

      CALL CCSDT_CLEAN_T3(TB3AM,NT1AMX,NVIRT,NRHFT)

      IF (PRINT_T3) THEN
        WRITE(LUPRI,*)'CCSDT_T31_AM> first-order T3 vector:'
        DO NAI = 1, NT1AMX
        DO NBJ = 1, NT1AMX
        DO NCK = 1, NT1AMX
          NAIBJCK = ((NAI-1)*NT1AMX+NBJ-1)*NT1AMX + NCK
          WRITE(LUPRI,*) NAIBJCK,TB3AM(NAI,NBJ,NCK)
        END DO
        END DO
        END DO
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_T31_NODDY                      *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XKSI3(XKSI3,FOCKB,T3AM,T2AM)
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
 
#if defined (SYS_CRAY)
      REAL HALF
#else
      DOUBLE PRECISION HALF
#endif
      PARAMETER (HALF = 0.5D0)
 
#if defined (SYS_CRAY)
      REAL AIBJCK
      REAL XKSI3(NT1AMX,NT1AMX,NT1AMX)
      REAL T2AM(NT1AMX,NT1AMX), FOCKB(NORBT,NORBT)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION AIBJCK
      DOUBLE PRECISION XKSI3(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX), FOCKB(NORBT,NORBT)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
    
      INTEGER NCK, NBJ, NAI, NBI, NCI, NAJ, NAK, NDJ, NBL
C
      DO NCK = 1,NT1AMX
         DO J = 1,NRHFT
            DO B = 1,NVIRT
               NBJ = NVIRT*(J-1) + B
               DO NAI = 1,NT1AMX
                  AIBJCK = 0.0D0
 
                  DO D = 1, NVIRT
                     NDJ = NVIRT*(J-1) + D
                     DO L = 1, NRHFT
                        NBL = NVIRT*(L-1) + B
                        AIBJCK = AIBJCK
     *                         - T2AM(NAI,NBL)
     *                          *T2AM(NCK,NDJ)
     *                          *FOCKB(L,NRHFT+D)
                     END DO
                  END DO

                  ! remember: T3 has turned sign
                  DO D = 1, NVIRT
                     NDJ = NVIRT*(J-1) + D
                     AIBJCK = AIBJCK
     *                      - HALF*T3AM(NAI,NDJ,NCK)
     *                            *FOCKB(NRHFT+B,NRHFT+D)
                  END DO
 
                  ! remember: T3 has turned sign
                  DO L = 1, NRHFT
                     NBL = NVIRT*(L-1) + B
                     AIBJCK = AIBJCK
     *                      + HALF*T3AM(NAI,NBL,NCK)
     *                         *FOCKB(L,J)
                  END DO
 
                  XKSI3(NAI,NBJ,NCK) = XKSI3(NAI,NBJ,NCK) + AIBJCK
                  XKSI3(NAI,NCK,NBJ) = XKSI3(NAI,NCK,NBJ) + AIBJCK
                  XKSI3(NBJ,NAI,NCK) = XKSI3(NBJ,NAI,NCK) + AIBJCK
                  XKSI3(NCK,NAI,NBJ) = XKSI3(NCK,NAI,NBJ) + AIBJCK
                  XKSI3(NBJ,NCK,NAI) = XKSI3(NBJ,NCK,NAI) + AIBJCK
                  XKSI3(NCK,NBJ,NAI) = XKSI3(NCK,NBJ,NAI) + AIBJCK
 
               END DO
            END DO
         END DO
      END DO

C----------------------------------------
C     Remove the forbidden elements.
C----------------------------------------

      CALL CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI3                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XKSI2(XKSI2,FOCKB,T3AM)
*---------------------------------------------------------------------*
*     Purpose: compute triples contribution to double excitation 
*              part of the Xski vector
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      REAL XKSI2(NT1AMX,NT1AMX)
      REAL FOCKB(NORBT,NORBT), CONTRIB
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION XKSI2(NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCKB(NORBT,NORBT), CONTRIB
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
    
      INTEGER NAI, NBJ, NBK, NCK, NCJ, INDEX

      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J 

      DO NAI = 1, NT1AMX
        DO J = 1,NRHFT
          DO B = 1,NVIRT
            NBJ = NVIRT*(J-1) + B
            DO K = 1,NRHFT
              NBK = NVIRT*(K-1) + B
              DO C = 1,NVIRT
                NCK = NVIRT*(K-1) + C
                NCJ = NVIRT*(J-1) + C
                CONTRIB = FOCKB(K,NRHFT+C) *
     &                        (T3AM(NAI,NBJ,NCK) - T3AM(NAI,NBK,NCJ))
                XKSI2(NBJ,NAI) = XKSI2(NBJ,NAI) - CONTRIB
                XKSI2(NAI,NBJ) = XKSI2(NAI,NBJ) - CONTRIB
 
              END DO
            END DO
          END DO
        END DO
      END DO
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI2                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_XKSI2_2(OMEGA2,FOCK,T3AM)
*---------------------------------------------------------------------*
*     Purpose: compute triples contribution to double excitation 
*              part of the Xski vector
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      REAL OMEGA2(NT1AMX,NT1AMX)
      REAL FOCK(NORBT,NORBT), XAIBJ
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION OMEGA2(NT1AMX,NT1AMX)
      DOUBLE PRECISION FOCK(NORBT,NORBT), XAIBJ
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif
      
      INTEGER NAI, NBJ, NBK, NAK, NCK, NCJ, NCI
    
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 120 J = 1,NRHFT
               DO 130 B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
C
                  DO 140 K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO 150 C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                        OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) -
     *      (T3AM(NAI,NBJ,NCK) - T3AM(NAI,NBK,NCJ))*FOCK(K,NRHFT+C)
C
C
  150                CONTINUE
  140             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      DO 200 NAI = 1,NT1AMX
         DO 210 NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
  210    CONTINUE
  200 CONTINUE
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_XKSI2                          *
*---------------------------------------------------------------------*
C  /* Deck ccsdt_t3am_r */
      SUBROUTINE CCSDT_T3AM_R(T3AM,FREQ,XINT1,XINT2,T2AM,SCR1,FOCKD)
C
C     response version of T3AM: T3AM needs to be initialized from
C                               the outside and a frequency is included
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(NORBT)
      DIMENSION T2AM(NT1AMX,NT1AMX)
      DIMENSION T3AM2(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (HALF = 0.5D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 50 I = 1,NRHFT
         DO 60 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
            SCR1(NAI) = FOCKD(NRHFT+A) - FOCKD(I)
   60    CONTINUE
   50 CONTINUE
C
      DO 100 NCK = 1,NT1AMX
C
         DO 110 J = 1,NRHFT
            DO 120 B = 1,NVIRT
C
               NBJ = NVIRT*(J-1) + B
C
               DO 130 NAI = 1,NT1AMX
C
                  AIBJCK = 0.0D0
                  DO 140 D = 1,NVIRT
C
                     NDJ = NVIRT*(J-1) + D
C
                     AIBJCK = AIBJCK + XINT1(NCK,B,D)*T2AM(NDJ,NAI)
C
  140             CONTINUE
C
                  DO 150 L = 1,NRHFT
C
                     NBL = NVIRT*(L-1) + B
C
                     AIBJCK = AIBJCK - XINT2(NCK,L,J)*T2AM(NBL,NAI)
C
  150             CONTINUE
C
                  AIBJCK = AIBJCK/(SCR1(NAI)+SCR1(NBJ)+SCR1(NCK)-FREQ)
C
                  T3AM(NAI,NBJ,NCK) = T3AM(NAI,NBJ,NCK) + AIBJCK
                  T3AM(NAI,NCK,NBJ) = T3AM(NAI,NCK,NBJ) + AIBJCK
                  T3AM(NBJ,NAI,NCK) = T3AM(NBJ,NAI,NCK) + AIBJCK
                  T3AM(NCK,NAI,NBJ) = T3AM(NCK,NAI,NBJ) + AIBJCK
                  T3AM(NBJ,NCK,NAI) = T3AM(NBJ,NCK,NAI) + AIBJCK
                  T3AM(NCK,NBJ,NAI) = T3AM(NCK,NBJ,NAI) + AIBJCK
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE

C----------------------------------------
C     Remove the forbidden elements.
C----------------------------------------
  
      CALL CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_T3AM_R                         *
*---------------------------------------------------------------------*
      SUBROUTINE CCSDT_FCK_R(FOCK,XLIAJB,T1AM)
C
C     Compute a response Fock matrix:
C        FOCK(kc) = FOCK(kc) + sum_ia L_kcia t_ia
C
#include <implicit.h>
      DIMENSION XLIAJB(NT1AMX,NT1AMX)
      DIMENSION T1AM(NT1AMX)
      DIMENSION FOCK(NORBT,NORBT)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>

      DO NK = 1, NRHFT
        DO NC = 1, NVIRT
          NCK = (NK-1)*NVIRT + NC
          DO NAI = 1, NT1AMX
            FOCK(NK,NRHFT+NC) = FOCK(NK,NRHFT+NC) + 
     &         XLIAJB(NCK,NAI) * T1AM(NAI)
          END DO
        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_FCK_R                          *
*---------------------------------------------------------------------*
*=====================================================================*
      SUBROUTINE CCSDT_CLEAN_T3(T3AM,NT1AMX,NVIRT,NRHFT)
*---------------------------------------------------------------------*
*     Purpose: remove forbidden elements in a triples vector
*---------------------------------------------------------------------*
      IMPLICIT NONE  

      INTEGER NT1AMX, NVIRT, NRHFT

#if defined (SYS_CRAY)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
#else
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
#endif

      INTEGER NA, NI, NAI, NB, NBI, NC, NCI, NJ, NAJ, NK, NAK
    
      DO NA = 1, NVIRT
       DO NI = 1, NRHFT
         NAI = NVIRT*(NI-1) + NA
         DO NB = 1, NVIRT
            NBI = NVIRT*(NI-1) + NB
            DO NC = 1, NVIRT
               NCI = NVIRT*(NI-1) + NC
               T3AM(NAI,NBI,NCI) = 0.0d0
               T3AM(NAI,NCI,NBI) = 0.0d0
               T3AM(NBI,NAI,NCI) = 0.0d0
               T3AM(NBI,NCI,NAI) = 0.0d0
               T3AM(NCI,NAI,NBI) = 0.0d0
               T3AM(NCI,NBI,NAI) = 0.0d0
            ENDDO
         ENDDO
         DO NJ = 1, NRHFT
            NAJ = NVIRT*(NJ-1) + NA
            DO NK = 1, NRHFT
               NAK = NVIRT*(NK-1) + NA
               T3AM(NAI,NAJ,NAK) = 0.0d0
               T3AM(NAI,NAK,NAJ) = 0.0d0
               T3AM(NAJ,NAI,NAK) = 0.0d0
               T3AM(NAJ,NAK,NAI) = 0.0d0
               T3AM(NAK,NAI,NAJ) = 0.0d0
               T3AM(NAK,NAJ,NAI) = 0.0d0
            ENDDO
         ENDDO
       ENDDO
      ENDDO
 
      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_CLEAN_T3                       *
*---------------------------------------------------------------------*
