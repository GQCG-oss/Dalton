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
C  /* Deck ccsdt_tran1_r */
      SUBROUTINE CCSDT_TRAN1_R(XINT1,XINT2,XLAMDP0,XLAMDH0,
     &                         XLAMDP1,XLAMDH1,AOINT,IDEL)
C
C     XINT1 = XINT1 + (K^p0 C^h0|B^p1 D^h0)
C     XINT2 = XINT2 + (K^p0 C^h0|L^p0 J^h1)
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP0(NBAST,NORBT), XLAMDH0(NBAST,NORBT)
      DIMENSION XLAMDP1(NBAST,NORBT), XLAMDH1(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 IGAM = 1,NBAST
         DO 110 IBET = 1,NBAST
            DO 120 IALP = 1,NBAST
               NAB = INDEX(IALP,IBET)
C
               if (aoint(nab,IGAM) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D)=XINT1(NCK,B,D)+AOINT(NAB,IGAM)
     &                     *XLAMDP0(IBET,K)      *XLAMDH0(IALP,NRHFT+C)
     &                     *XLAMDP1(IGAM,NRHFT+B)*XLAMDH0(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
               DO 300 J = 1,NRHFT
                  DO 310 L = 1,NRHFT
                     DO 320 K = 1,NRHFT
                        DO 330 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT2(NCK,L,J)=XINT2(NCK,L,J)+AOINT(NAB,IGAM)
     &                      * XLAMDP0(IBET,K) * XLAMDH0(IALP,NRHFT+C)
     &                      * XLAMDP0(IGAM,L) * XLAMDH1(IDEL,J)
C
  330                   CONTINUE
  320                CONTINUE
  310             CONTINUE
  300          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_tran3_r */
      SUBROUTINE CCSDT_TRAN3_R(XINT1,XINT2,XLAMDP0,XLAMDH0,
     &                         XLAMDP1,XLAMDH1,XLAMDP2,XLAMDH2,
     &                         AOINT,IDEL)
C
C     XINT1 = XINT1 + (C^p1 K^h1|B^p2 D^h0)
C     XINT2 = XINT2 + (C^p1 K^h1|L^p0 J^h2)
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP0(NBAST,NORBT), XLAMDH0(NBAST,NORBT)
      DIMENSION XLAMDP1(NBAST,NORBT), XLAMDH1(NBAST,NORBT)
      DIMENSION XLAMDP2(NBAST,NORBT), XLAMDH2(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)+AOINT(NAB,G)
     &                      * XLAMDP1(A,NRHFT+C) * XLAMDH1(IB,K)
     &                      * XLAMDP2(G,NRHFT+B) * XLAMDH0(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
               DO 300 J = 1,NRHFT
                  DO 310 L = 1,NRHFT
                     DO 320 K = 1,NRHFT
                        DO 330 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)+AOINT(NAB,G)
     &                      * XLAMDP1(A,NRHFT+C) * XLAMDH1(IB,K)
     &                      * XLAMDP0(G,L)       * XLAMDH2(IDEL,J)
C
  330                   CONTINUE
  320                CONTINUE
  310             CONTINUE
  300          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccfop_tran1 */
      SUBROUTINE CCFOP_TRAN1_R(XINT1,XINT2,XINT3,XINT4,
     &                         XLAMDP0,XLAMDH0,
     &                         XLAMDP1,XLAMDH1,
     &                         XLAMDP2,XLAMDH2,
     &                         AOINT,IDEL)
C
C     XINT1 = (K^p0 C^h0|D^p1 L^h2)    "(O-0 V-0|V-1 O-2)"
C     XINT2 = (K^p0 L^h1|C^p2 D^h0)    "(O-0 O-1|V-2 V-0)"
C     XINT3 = (K^p0 L^h1|M^p0 N^h2)    "(O-0 O-1|O-0 O-2)"
C     XINT4 = (C^p1 D^h0|E^p2 F^h0)    "(V-1 V-0|V-2 V-0)"
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NRHFT,NVIRT,NVIRT,NRHFT)
      DIMENSION XINT2(NRHFT,NRHFT,NVIRT,NVIRT)
      DIMENSION XINT3(NRHFT,NRHFT,NRHFT,NRHFT)
      DIMENSION XINT4(NVIRT,NVIRT,NVIRT,NVIRT)
      DIMENSION XLAMDP0(NBAST,NORBT), XLAMDH0(NBAST,NORBT)
      DIMENSION XLAMDP1(NBAST,NORBT), XLAMDH1(NBAST,NORBT)
      DIMENSION XLAMDP2(NBAST,NORBT), XLAMDH2(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
C
      LOGICAL LDEBUG
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      LDEBUG = .TRUE.
C
C----------------------------------------
C     Calculate integrals :
C----------------------------------------
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
C
               DO NC = 1,NVIRT
                  DO ND = 1,NVIRT
                     DO NK = 1,NRHFT
                        DO NL = 1,NRHFT
C
                           XINT1(NK,NC,ND,NL) = XINT1(NK,NC,ND,NL)
     &                       + AOINT(NAB,G) *
     &                          XLAMDP0(A,NK) * 
     &                          XLAMDH0(IB,NRHFT+NC) *
     &                          XLAMDP1(G,NRHFT+ND) * 
     &                          XLAMDH2(IDEL,NL)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               DO NC = 1,NVIRT
                  DO ND = 1,NVIRT
                     DO NK = 1,NRHFT
                        DO NL = 1,NRHFT
C
                           XINT2(NK,NL,NC,ND) = XINT2(NK,NL,NC,ND)
     &                       + AOINT(NAB,G) *
     &                          XLAMDP0(A,NK) * 
     &                          XLAMDH1(IB,NL) *
     &                          XLAMDP2(G,NRHFT+NC) *
     &                          XLAMDH0(IDEL,NRHFT+ND)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               DO NK = 1,NRHFT
                  DO NL = 1,NRHFT
                     DO NM = 1,NRHFT
                        DO NN = 1,NRHFT
C
                           XINT3(NK,NL,NM,NN) = XINT3(NK,NL,NM,NN)
     &                       + AOINT(NAB,G) *
     &                          XLAMDP0(A,NK) *
     &                          XLAMDH1(IB,NL) *
     &                          XLAMDP0(G,NM) *
     &                          XLAMDH2(IDEL,NN)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
               DO NC = 1,NVIRT
                  DO ND = 1,NVIRT
                     DO NE = 1,NVIRT
                        DO NF = 1,NVIRT
C
                           XINT4(NC,ND,NE,NF) = XINT4(NC,ND,NE,NF)
     &                       + AOINT(NAB,G) * 
     &                          XLAMDP1(A,NRHFT+NC) *
     &                          XLAMDH0(IB,NRHFT+ND) *
     &                          XLAMDP2(G,NRHFT+NE) *
     &                          XLAMDH0(IDEL,NRHFT+NF)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
*=====================================================================*
      SUBROUTINE CCSDT_INTS0_NODDY(DO_INTXY,XIAJB,YIAJB,
     &                             DO_INTS0,XINT1S,XINT2S,
     &                             DO_INTT0,XINT1T,XINT2T,
     &                             XLAMDP,XLAMDH,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute some standard integrals for CC3
*
*             if (do_intxy) compute 
*                 XIAJB  = 2(kc|ld) - (kd|lc)
*                 YIAJB  =  (kc|ld)
*
*             if (do_ints0) compute 
*                 XINT1S = (ck|bd) and 
*                 XINT2S = (ck|lj)
*
*             if (do_intt0) compute 
*                 XINT1T = (kc|bd) and 
*                 XINT2T = (kc|lj)
*
* Written by Christof Haettig, November 2002, based on CCSDT_XI3_NODDY.
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
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0=1)

      LOGICAL DO_INTXY, DO_INTS0, DO_INTT0
      INTEGER LWORK
      REAL*8  WORK(LWORK)
      REAL*8  XINT1T(NT1AMX*NVIRT*NVIRT)
      REAL*8  XINT2T(NRHFT*NRHFT*NT1AMX)
      REAL*8  XINT1S(NT1AMX*NVIRT*NVIRT)
      REAL*8  XINT2S(NRHFT*NRHFT*NT1AMX)
      REAL*8  XIAJB(NT1AMX*NT1AMX)
      REAL*8  YIAJB(NT1AMX*NT1AMX)
      REAL*8  XLAMDP(*), XLAMDH(*)

      INTEGER ISYMD, ILLL, IDEL, ISYDIS, KXINT, KEND1, LWRK1

      CALL QENTER('CCSDT_INTS0_NODDY')

      IF (DIRECT) 
     &  CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_INTS0_NODDY')

*---------------------------------------------------------------------*
*     Loop over distributions of integrals:
*---------------------------------------------------------------------*
      IF (DO_INTS0) THEN
        CALL DZERO(XINT1S,NT1AMX*NVIRT*NVIRT)
        CALL DZERO(XINT2S,NT1AMX*NRHFT*NRHFT)
      END IF
      IF (DO_INTT0) THEN
         CALL DZERO(XINT1T,NT1AMX*NVIRT*NVIRT)
         CALL DZERO(XINT2T,NT1AMX*NRHFT*NRHFT)
      END IF
      IF (DO_INTXY) THEN
         CALL DZERO(XIAJB,NT1AMX*NT1AMX)
         CALL DZERO(YIAJB,NT1AMX*NT1AMX)
      END IF

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYM0)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = 1
            KEND1  = KXINT + NDISAO(ISYDIS)
            LWRK1  = LWORK - KEND1
            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_INTS0_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND1),LWRK1,
     *                  IDUMMY,DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            IF (DO_INTT0) THEN
               CALL CCSDT_TRAN1(XINT1T,XINT2T,XLAMDP,XLAMDH,
     *                          WORK(KXINT),IDEL)
            END IF

            IF (DO_INTXY) THEN
               CALL CC3_TRAN2(XIAJB,YIAJB,XLAMDP,XLAMDH,
     *                        WORK(KXINT),IDEL)
            END IF
 
            IF (DO_INTS0) THEN
              CALL CCSDT_TRAN3(XINT1S,XINT2S,XLAMDP,XLAMDH,
     *                         WORK(KXINT),IDEL)
            END IF
 
         END DO   
      END DO  
*---------------------------------------------------------------------*
*     End Loop over distributions of integrals.
*---------------------------------------------------------------------*
 
      CALL QEXIT('CCSDT_INTS0_NODDY')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_INTS0_NODDY                    *
*---------------------------------------------------------------------*
      SUBROUTINE CCSDT_INIT_NODDY(WORK,LWORK,DO_L03AM)
*---------------------------------------------------------------------*
*
* Purpose: precompute some intermediates for CC3
*
*          for DO_L03AM=.true. the triples part of the zeroth-order
*          langrangian multipliers will be computed and save on file.
*          for this case it will be assumed that all integrals
*          are already present of file
*  
* Written by Christof Haettig, Mai 2003
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <ccnoddy.h>
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0=1)

      INTEGER LWORK
      REAL*8  WORK(LWORK), FF, DDOT

      LOGICAL DO_L03AM

      CHARACTER*10 MODEL
      INTEGER KSCR1, KFOCKD, KEND1, LWRK1, KFOCK0, KFIELD, KFIELDAO,
     &        KLAMP0, KLAMH0, KINT1T0, KINT2T0, KINT1S0, KINT2S0,
     &        KXIAJB, KYIAJB, K0IOVVO, K0IOOVV, K0IOOOO, K0IVVVV,
     &        KT01AM, KT02AM, KT03AM, KT3SCR, LUSIFC, IOPT, LUFOCK,
     &        LUTEMP, ILLL, IDEL, ISYMD, ISYDIS, KEND2, LWRK2,
     &        KXINT

      CALL QENTER('CCTINI')

      IF (DIRECT) 
     &  CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_INIT_NODDY')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KEND1   = KFOCKD + NORBT

      KFOCK0  = KEND1
      KEND1   = KFOCK0  + NORBT*NORBT
     
      IF (NONHF) THEN
        KFIELD   = KEND1
        KFIELDAO = KFIELD   + NORBT*NORBT
        KEND1    = KFIELDAO + NORBT*NORBT
      END IF
      
      KLAMP0  = KEND1
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KINT1S0 = KEND1
      KINT2S0 = KINT1S0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2S0 + NRHFT*NRHFT*NT1AMX

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      K0IOVVO = KEND1
      K0IOOVV = K0IOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      K0IOOOO = K0IOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      K0IVVVV = K0IOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1   = K0IVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KT01AM  = KEND1
      KT02AM  = KT01AM + NT1AMX
      KT03AM  = KT02AM + NT1AMX*NT1AMX
      KT3SCR  = KT03AM + NT1AMX*NT1AMX*NT1AMX
      KEND1   = KT3SCR + NT1AMX*NT1AMX*NT1AMX
              
      LWRK1   = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_INIT_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT01AM),DUMMY)

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT01AM),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*     If needed get external field:
*---------------------------------------------------------------------*
      IF (NONHF) THEN
        CALL DZERO(WORK(KFIELDAO),NORBT*NORBT)
        DO I = 1, NFIELD
          FF = EFIELD(I)
          CALL CC_ONEP(WORK(KFIELDAO),WORK(KEND1),LWRK1,FF,1,LFIELD(I))
        ENDDO
        CALL DCOPY(NORBT*NORBT,WORK(KFIELDAO),1,WORK(KFIELD),1)

        ! calculate external field in zero-order lambda basis
        CALL CC_FCKMO(WORK(KFIELD),WORK(KLAMP0),WORK(KLAMH0),
     *                WORK(KEND1),LWRK1,1,1,1)

        IF (LOCDBG) WRITE(LUPRI,*) 'NORM^2(FIELD):',
     &     DDOT(NORBT*NORBT,WORK(KFIELD),1,WORK(KFIELD),1)
      ENDIF

*---------------------------------------------------------------------*
*     Compute some integrals:
*           XINT1S0 =  (CK|BD)
*           XINT2S0 =  (CK|LJ)
*           XINT1T0 =  (KC|BD)
*           XINT2T0 =  (KC|LJ)
*           XIAJB   = 2(IA|JB) - (IB|JA)
*           YIAJB   =  (IA|JB)
*---------------------------------------------------------------------*
      IF (DO_L03AM) THEN
        CALL CCSDT_READ_NODDY(.FALSE.,WORK(KFOCKD),WORK(KFOCK0),
     &                                WORK(KFIELD),WORK(KFIELDAO),
     &                        .TRUE., WORK(KXIAJB),WORK(KYIAJB),
     &                        .TRUE., WORK(KINT1S0),WORK(KINT2S0),
     &                        .TRUE., WORK(KINT1T0),WORK(KINT2T0),
     &                        .FALSE.,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX) 

      ELSE
        CALL CCSDT_INTS0_NODDY(.TRUE.,WORK(KXIAJB), WORK(KYIAJB),
     &                         .TRUE.,WORK(KINT1S0),WORK(KINT2S0),
     &                         .TRUE.,WORK(KINT1T0),WORK(KINT2T0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KEND1),LWRK1)
      END IF

*---------------------------------------------------------------------*
*     compute and save zeroth-order triples amplitudes:
*---------------------------------------------------------------------*
      IF (.NOT. DO_L03AM) THEN
        ! read T^0 doubles amplitudes from file and square up 
        IOPT   = 2
        Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KT03AM))
        CALL CC_T2SQ(WORK(KT03AM),WORK(KT02AM),ISYM0)

        ! compute triples amplitudes
        CALL CCSDT_T03AM(WORK(KT03AM),WORK(KINT1S0),WORK(KINT2S0),
     *                   WORK(KT02AM),WORK(KSCR1),WORK(KFOCKD),
     *                   WORK(KFIELD),WORK(KT3SCR))

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT03AM),1)

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODT30,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        WRITE(LUTEMP) (WORK(KT03AM+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF

*---------------------------------------------------------------------*
*     compute and save zeroth-order triples multipliers:
*     (note that this section overwrite the cluster amplitudes on
*      KT01AM, KT02AM and KT03AM...)
*---------------------------------------------------------------------*
      IF (DO_L03AM) THEN
        ! read L^0 singles and doubles multipliers from file 
        IOPT   = 3
        Call CC_RDRSP('L0',0,ISYM0,IOPT,MODEL,
     *                 WORK(KT01AM),WORK(KT03AM))
        CALL CC_T2SQ(WORK(KT03AM),WORK(KT02AM),ISYM0)

        ! compute triples multipliers
        CALL CCSDT_L03AM(WORK(KT03AM),WORK(KINT1T0),WORK(KINT2T0),
     *                   WORK(KXIAJB),WORK(KFOCK0),WORK(KT01AM),
     *                   WORK(KT02AM),WORK(KSCR1),WORK(KFOCKD),
     *                   WORK(KFIELD),WORK(KT3SCR))

        CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KT03AM),1)

        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODL30,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        WRITE(LUTEMP) (WORK(KT03AM+I-1), I=1,NT1AMX*NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')

        CALL QEXIT('CCTINI')
        RETURN
      END IF

*---------------------------------------------------------------------*
*     Compute integrals needed for the following contributions:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(K0IOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(K0IOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(K0IOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(K0IVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

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
               CALL QUIT('Insufficient space in CCSDT_INIT_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  IDUMMY,DIRECT) ! NB! DIRECT always .false. here
 
            CALL CCFOP_TRAN1_R(WORK(K0IOVVO),WORK(K0IOOVV),
     &                         WORK(K0IOOOO),WORK(K0IVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

*---------------------------------------------------------------------*
*     Save fock-like intermediates and Lambda matrices
*---------------------------------------------------------------------*
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODFOCK,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE(LUTEMP) (WORK(KFOCKD+I-1),  I=1,NORBT)
      WRITE(LUTEMP) (WORK(KFOCK0+I-1),  I=1,NORBT*NORBT)
      IF (NONHF) THEN
        WRITE(LUTEMP) (WORK(KFIELD+I-1),   I=1,NORBT*NORBT)
        WRITE(LUTEMP) (WORK(KFIELDAO+I-1), I=1,NORBT*NORBT)
      END IF
      CALL GPCLOSE(LUTEMP,'KEEP')

*---------------------------------------------------------------------*
*     Save integrals XINT1T0 and XINT2T0
*---------------------------------------------------------------------*
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODINTT,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE(LUTEMP) (WORK(KINT1T0+I-1), I=1,NT1AMX*NVIRT*NVIRT)
      WRITE(LUTEMP) (WORK(KINT2T0+I-1), I=1,NT1AMX*NRHFT*NRHFT)
      CALL GPCLOSE(LUTEMP,'KEEP')

*---------------------------------------------------------------------*
*     Save integrals XINT1S0 and XINT2S0
*---------------------------------------------------------------------*
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODINTS,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE(LUTEMP) (WORK(KINT1S0+I-1), I=1,NT1AMX*NVIRT*NVIRT)
      WRITE(LUTEMP) (WORK(KINT2S0+I-1), I=1,NT1AMX*NRHFT*NRHFT)
      CALL GPCLOSE(LUTEMP,'KEEP')
 
*---------------------------------------------------------------------*
*     Save integrals XIAJB and YIAJB
*---------------------------------------------------------------------*
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODINTX,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE(LUTEMP) (WORK(KXIAJB+I-1),  I=1,NT1AMX*NT1AMX)
      WRITE(LUTEMP) (WORK(KYIAJB+I-1),  I=1,NT1AMX*NT1AMX)
      CALL GPCLOSE(LUTEMP,'KEEP')
 
*---------------------------------------------------------------------*
*     Save integrals XOOOO, XOVVO, XOOVV and XVVVV
*---------------------------------------------------------------------*
      LUTEMP = -1
      CALL GPOPEN(LUTEMP,FILNODINTV,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE(LUTEMP) (WORK(K0IOVVO+I-1), I=1,NRHFT*NVIRT*NVIRT*NRHFT)
      WRITE(LUTEMP) (WORK(K0IOOVV+I-1), I=1,NRHFT*NVIRT*NVIRT*NRHFT)
      WRITE(LUTEMP) (WORK(K0IOOOO+I-1), I=1,NRHFT*NRHFT*NRHFT*NRHFT)
      WRITE(LUTEMP) (WORK(K0IVVVV+I-1), I=1,NVIRT*NVIRT*NVIRT*NVIRT)
      CALL GPCLOSE(LUTEMP,'KEEP')

      CALL QEXIT('CCTINI')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_INIT_NODDY                     *
*---------------------------------------------------------------------*
      SUBROUTINE CCSDT_READ_NODDY(RD_FOCK0,FOCKD,FOCK0,FIELD,FIELDAO,
     &                            RD_INTXY,XIAJB,YIAJB,
     &                            RD_INTS0,XINT1S0,XINT2S0,
     &                            RD_INTT0,XINT1T0,XINT2T0,
     &                            RD_INTV0,XOVVO,XOOVV,XOOOO,XVVVV,
     &                            NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX)
*---------------------------------------------------------------------*
*
* Purpose: read precomputed intermediates for CC3
*
* Written by Christof Haettig, Mai 2003
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <ccfield.h>
#include <ccnoddy.h>
#include <ccsdinp.h>
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)

      INTEGER NORBT,NLAMDT,NRHFT,NVIRT,NT1AMX, LUTEMP, I

      LOGICAL RD_FOCK0, RD_INTXY, RD_INTS0, RD_INTT0, RD_INTV0, RD_XLAM0

      REAL*8  FOCKD(*), FOCK0(*)
      REAL*8  FIELD(*), FIELDAO(*), XIAJB(*), YIAJB(*)
      REAL*8  XINT1S0(*), XINT2S0(*), XINT1T0(*), XINT2T0(*)
      REAL*8  XOVVO(*), XOOVV(*), XOOOO(*), XVVVV(*)

      CALL QENTER('CCTRDI')

      IF (DIRECT) 
     &  CALL QUIT('DIRECT NOT IMPLEMENTED IN CCSDT_READ_NODDY')

*---------------------------------------------------------------------*
*     Read integrals and triples cluster amplitudes on file
*---------------------------------------------------------------------*
      IF (RD_FOCK0) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODFOCK,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (FOCKD(I),I=1,NORBT)
        READ(LUTEMP) (FOCK0(I),I=1,NORBT*NORBT)
        IF (NONHF) THEN
          READ(LUTEMP) (FIELD(I),  I=1,NORBT*NORBT)
          READ(LUTEMP) (FIELDAO(I),I=1,NORBT*NORBT)
        END IF
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF

      IF (RD_INTT0) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODINTT,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (XINT1T0(I),I=1,NT1AMX*NVIRT*NVIRT)
        READ(LUTEMP) (XINT2T0(I),I=1,NT1AMX*NRHFT*NRHFT)
        CALL GPCLOSE(LUTEMP,'KEEP')
      ENDIF

      IF (RD_INTS0) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODINTS,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (XINT1S0(I),I=1,NT1AMX*NVIRT*NVIRT)
        READ(LUTEMP) (XINT2S0(I),I=1,NT1AMX*NRHFT*NRHFT)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF
 
      IF (RD_INTXY) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODINTX,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (XIAJB(I),I=1,NT1AMX*NT1AMX)
        READ(LUTEMP) (YIAJB(I),I=1,NT1AMX*NT1AMX)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF
 
      IF (RD_INTV0) THEN
        LUTEMP = -1
        CALL GPOPEN(LUTEMP,FILNODINTV,'UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        READ(LUTEMP) (XOVVO(I),I=1,NRHFT*NVIRT*NVIRT*NRHFT)
        READ(LUTEMP) (XOOVV(I),I=1,NRHFT*NVIRT*NVIRT*NRHFT)
        READ(LUTEMP) (XOOOO(I),I=1,NRHFT*NRHFT*NRHFT*NRHFT)
        READ(LUTEMP) (XVVVV(I),I=1,NVIRT*NVIRT*NVIRT*NVIRT)
        CALL GPCLOSE(LUTEMP,'KEEP')
      END IF

      CALL QEXIT('CCTRDI')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_READ_NODDY                     *
*---------------------------------------------------------------------*
