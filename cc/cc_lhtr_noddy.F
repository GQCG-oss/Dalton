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
C  /* Deck ccsd_triple */
      SUBROUTINE CC_LHTR_NODDY(OMEGA1,OMEGA2,T1AM,T2AM,C1AM,C2AM,
     *                         XLAMDP,XLAMDH,WORK,LWORK,IOPT)
C
C     Written by Henrik Koch 19-Sep-1994
C
C     Version 1.0
C
C     Purpose:
C
C     Calculate the iterative triples corrections to the CCSD
C     equations.
C
C     NB! The T2 amplitudes are assumed to be a full square.
C
C
C     NB! It is assumed that the vectors are allocated in the following
C     order:
C           T1AM(*), OMEGA1(*), OMEGA2(*), T2AM(*), SCR(*), WRK(*).
C
C     IOPT = 0 No storage of T3, L3
C     IOPT = 1 Store T3 and L3 in 'T3AMPL' and 'L3AMPL'
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION OMEGA1(*),OMEGA2(*),T1AM(*),T2AM(*), C1AM(*), C2AM(*)
      DIMENSION XLAMDP(*),XLAMDH(*),WORK(LWORK)
C
#include <inforb.h>
#include <infind.h>
#include <blocks.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
COMMENT COMMENT
COMMENT COMMENT
#include <mxcent.h>
#include <eribuf.h>
COMMENT COMMENT
COMMENT COMMENT
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (NSYM .NE. 1) THEN
         CALL QUIT('No symmetry in this part yet')
      ENDIF
C
C------------------------
C     Dynamic Allocation.
C------------------------
C
      KCMO   = 1
      KSCR1  = KCMO   + NLAMDT
      KFOCK  = KSCR1  + NT1AMX
      KFOCKD = KFOCK  + N2BST(ISYMOP)
      KINT1  = KFOCKD + NORBT
      KINT2  = KINT1  + NT1AMX*NVIRT*NVIRT
      KINT1T = KINT2  + NRHFT*NRHFT*NT1AMX
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KINT1S = KINT2T + NRHFT*NRHFT*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KINT3S = KINT2S + NRHFT*NRHFT*NT1AMX
      KINT4S = KINT3S + NT1AMX*NVIRT*NVIRT
      KT3AM  = KINT4S + NRHFT*NRHFT*NT1AMX
      KL3AM  = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KOME1  = KL3AM  + NT1AMX*NT1AMX*NT1AMX
      KOME2  = KOME1  + NT1AMX
      KOME22 = KOME2  + NT1AMX*NT1AMX
      KIAJB  = KOME22 + NT1AMX*NT1AMX
      KSIG2  = KIAJB  + NT1AMX*NT1AMX
      KYIAJB = KSIG2  + NT1AMX*NT1AMX
      KEND1  = KYIAJB + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
C
C     New Integrals.
      KIOVVO = KEND1
      KIOOVV = KIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KIOOOO = KIOOVV + NRHFT*NRHFT*NVIRT*NVIRT
      KIVVVV = KIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1  = KIVVVV + NVIRT*NVIRT*NVIRT*NVIRT
      LWRK1  = LWORK  - KEND1
C
COMMENT COMMENT
COMMENT COMMENT
      ktmat = kend1
      kend1 = ktmat + nvirt*nrhft*nrhft*nrhft
      lwrk1 = lwork - kend1
      call dzero(work(ktmat),nvirt*nrhft*nrhft*nrhft)
COMMENT COMMENT
COMMENT COMMENT
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSD_TRIPLE')
      ENDIF
C
C--------------------------------
C     Initialize integral arrays.
C--------------------------------
C
      CALL DZERO(WORK(KINT1),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT3S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT4S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KT3AM),NT1AMX*NT1AMX*NT1AMX)
      CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME22),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KSIG2),NT1AMX*NT1AMX)
C
      CALL DZERO(WORK(KIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
C
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDT)
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C--------------------------------------------------
C     Read in AO fock matrix and transform to MO
C--------------------------------------------------
C
      LUFCK = 0
C     The T1 transformed density is used to calculate this AO fock matrix.
      CALL GPOPEN(LUFCK,'CC_FCKH','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
C     The CMO transformed density is used to calculate this AO fock matrix.
C      CALL GPOPEN(LUFCK,'CC_FCKREF','OLD',' ','UNFORMATTED',
C     *            IDUMMY,.FALSE.)
C
      REWIND(LUFCK)
      READ(LUFCK)(WORK(KFOCK + I-1),I = 1,N2BST(ISYMOP))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'NODDY : Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFOCK),ISYMOP)
      ENDIF
C
      CALL CC_FCKMO(WORK(KFOCK),XLAMDP,XLAMDH,
     *              WORK(KEND1),LWRK1,1,1,1)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('OUTPUT OF FOCK')
         CALL OUTPUT(WORK(KFOCK),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      ENDIF
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      IF (DIRECT) THEN
         NTOSYM = 1
         DTIME  = SECOND()
c        CALL HERDI1(WORK(KEND1),LWRK1,IPRINT)
         CALL QUIT('Direct not implemented')
         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            NTOT = MAXSHL
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed t2am.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
               DTIME  = SECOND()
c              CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,IPRINT)
               CALL QUIT('Direct not implemented')
               DTIME  = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
COMMENT COMMENT
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_FOP3')
               END IF
COMMENT COMMENT
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCSD_TRIPLE')
               ENDIF
C
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
C
C----------------------------------------------------
C              Calculate integrals needed in CCSD[T].
C----------------------------------------------------
C
c              CALL CCSDT_TRAN1(WORK(KINT1),WORK(KINT2),WORK(KCMO),
c    *                          WORK(KCMO),WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
C               CALL CCSDT_TRAN2(WORK(KIAJB),WORK(KYIAJB),WORK(KCMO),
C     *                          WORK(KXINT),IDEL)
C
               CALL CC3_TRAN2(WORK(KIAJB),WORK(KYIAJB),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN3(WORK(KINT1S),WORK(KINT2S),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CC3_TRAN3(WORK(KINT3S),WORK(KINT4S),XLAMDP,
     *                        XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCFOP_TRAN1(WORK(KIOVVO),WORK(KIOOVV),
     *                          WORK(KIOOOO),WORK(KIVVVV),
     *                          XLAMDP,XLAMDH,WORK(KXINT),IDEL)
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
COMMENT COMMENT
COMMENT COMMENT
C      CALL KHTCME1S(work(kint1s),work(kint3s),work(kint2s),work(kint4s))
C      CALL DCOPY(nt1amx*nvirt*nvirt,work(kint1s),1,work(kint3s),1)
C      CALL DCOPY(nt1amx*nrhft*nrhft,work(kint2s),1,work(kint4s),1)
CC
CC
CC
CC
CC
C      CALL KHTCME1T(work(kint1t),work(kint3s),work(kint2t),work(kint4s))
CC
C      CALL DCOPY(nt1amx*nvirt*nvirt,work(kint1s),1,work(kint3s),1)
C      CALL DCOPY(nt1amx*nrhft*nrhft,work(kint2s),1,work(kint4s),1)
CC
CC
CC
CC
      if (.false.) then
         call new_int_output(work(kiovvo),work(kioovv),work(kioooo),
     *                        work(kivvvv))
C
         xnorm = ddot(NRHFT*NVIRT*NVIRT*NRHFT,work(kiovvo),1,
     *                work(kiovvo),1)
         write(lupri,*) 'Norm of xovvo in noddy : ',xnorm
         xnorm = ddot(NRHFT*NVIRT*NVIRT*NRHFT,work(kioovv),1,
     *                work(kioovv),1)
         write(lupri,*) 'Norm of xoovv in noddy : ',xnorm
      endif
COMMENT COMMENT
COMMENT COMMENT
C-------------------------------------
C     Calculate the T3 amplitudes.
C-------------------------------------
C
      if (.true.) then
      CALL CCSDT_T3AM(WORK(KT3AM),WORK(KINT1S),WORK(KINT2S),T2AM,
     *                WORK(KSCR1),WORK(KFOCKD),WORK(KL3AM),DUMMY)
      endif

      CALL DSCAL(nt1amx*nt1amx*nt1amx,-1.0D0,work(kt3am),1)
C
      IF (IOPT .EQ. 1) THEN
         LUSIFC = 0
         CALL GPOPEN(LUSIFC,'T3AMPL','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
         WRITE (LUSIFC) (WORK(KT3AM-1+I), I=1,NT1AMX*NT1AMX*NT1AMX)
         CALL GPCLOSE(LUSIFC,'KEEP')
      ENDIF
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
      if (.false.) then
         call temppr(work(KT3AM),1)
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C------------------------------------
C     Calculate the CC3 corrections.
C------------------------------------
C
      if (.true.) then
      CALL T3_LEFT1(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,WORK(KT3AM))
      endif
C
      if (.true.) then
      CALL T3_LEFT2(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,WORK(KT3AM))
      endif
C
      if (.true.) then
C
C     L instead of g
C
      CALL T3_LEFT3(WORK(KOME1),WORK(KIAJB),C2AM,WORK(KT3AM))
      endif
C
      DO I = 1, NT1AMX
         IF (abs(work(kome1-1+i)) .gt. 1.0D-12) THEN
            write(lupri,*) 'OMEGA1_NODDY_T3(',i,')',work(kome1-1+i)
         ENDIF
      ENDDO
C
      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      ENDDO
C
C-------------------------------------
C     Calculate the L3 amplitudes.
C-------------------------------------
C
      if (.true.) then
COMMENT      CALL CCSDT_L3AM(WORK(KL3AM),WORK(KINT3S),WORK(KINT4S),
COMMENT     *                WORK(KIAJB),WORK(KFOCK),C1AM,C2AM,
COMMENT     *                WORK(KSCR1),WORK(KFOCKD))
      CALL CCSDT_L3AM(WORK(KL3AM),WORK(KINT1T),WORK(KINT2T),
     *                WORK(KIAJB),WORK(KFOCK),C1AM,C2AM,
     *                WORK(KSCR1),WORK(KFOCKD))
      endif
C
      CALL DSCAL(nt1amx*nt1amx*nt1amx,-1.0D0,work(kl3am),1)
C
      IF (IOPT .EQ. 1) THEN
         LUSIFC = 0
         CALL GPOPEN(LUSIFC,'L3AMPL','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
         WRITE (LUSIFC) (WORK(KL3AM-1+I), I=1,NT1AMX*NT1AMX*NT1AMX)
         CALL GPCLOSE(LUSIFC,'KEEP')
      ENDIF
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
      if (.false.) then
         call temppr(work(KL3AM),3)
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C----------------------------------------------------------
C         Calculate the contributions to Omega2 from L3
C         The result is written out inside the routine.
C----------------------------------------------------------
C
          CALL DZERO(WORK(KOME1),NT1AMX)
          CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
          CALL DZERO(WORK(KOME22),NT1AMX*NT1AMX)
C
       write(lupri,*) 'Before cc3_l3_omega'
          CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     *                             WORK(KIOOOO),WORK(KIOVVO),
     *                             WORK(KIOOVV),WORK(KIVVVV),T2AM)
C
          CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),WORK(KL3AM),
     *                             WORK(KINT1S),WORK(KINT2S))
C
      DO 300 I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
  300 CONTINUE
C
      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         ENDDO
      ENDDO
C
C      do i = 1, nt1amx
C         if (abs(omega1(i)) .gt. 1.0D-12) then
C            write(lupri,*) 'Omega1_total_noddy(',i,') = ',omega1(i)
C         endif
C      enddo
C      do i = 1, nt2am(1)
C         if (abs(omega2(i)) .gt. 1.0D-12) then
C            write(lupri,*) 'Omega2_total_noddy(',i,') = ',omega2(i)
C         endif
C      enddo
C
      XNORM = DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
      write(lupri,*) 'Norm Omega1 at the end of noddy :',xnorm
      XNORM = DDOT(NT2AM(1),OMEGA2,1,OMEGA2,1)
      write(lupri,*) 'Norm Omega2 at the end of noddy :',xnorm
C
      RETURN
      END
C  /* Deck temppr */
      SUBROUTINE TEMPPR(T3AM,IOPT)
C
#include <implicit.h>
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL AROUND('ENTERED TEMPPR')
      WRITE(LUPRI,*) 'IOPT  = ',IOPT
C
      DO NA = 1,NVIRT
      DO NB = 1,NVIRT
      DO NC = 1,NVIRT
         DO NI = 1,NRHFT
         NAI = NVIRT*(NI-1) + NA
            DO NJ = 1,NRHFT
               NBJ = NVIRT*(NJ-1) + NB
               DO NK = 1,NRHFT
                  NCK = NVIRT*(NK-1) + NC
C
               if (abs(t3am(nai,nbj,nck)) .gt. 1.0D-12) then
                 if (iopt .eq. 1) then
           write(lupri,1) 't3am (^{',na,',',nb,',',nc,'}_{',ni,',',
     *                                nj,',',nk,'}) -> (',
     *                      NAI,',',NBJ,',',NCK,') = ',t3am(nai,nbj,nck)
                 else if (iopt .eq. 2) then
           write(lupri,1) 't3bar(^{',na,',',nb,',',nc,'}_{',ni,',',
     *                                nj,',',nk,'}) -> (',
     *                      NAI,',',NBJ,',',NCK,') = ',t3am(nai,nbj,nck)
                 else if (iopt .eq. 3) then
           write(lupri,1) 'L3am (^{',na,',',nb,',',nc,'}_{',ni,',',
     *                                nj,',',nk,'}) -> (',
     *                      NAI,',',NBJ,',',NCK,') = ',t3am(nai,nbj,nck)
                 else
           call quit('Wrong IOPT in temppr')
                 endif
               endif
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
C
   1  FORMAT(1x,A8,I3,A1,I3,A1,I3,A3,I3,A1,I3,A1,I3,A7,I3,A1,I3,A1,I3,
     *       A4,E20.10)
C
      RETURN
      END
C  /* Deck t3_left1 */
      SUBROUTINE T3_LEFT1(OMEGA1,YIAJB,XIAJB,C2AM,T3AM)
C
C
C
#include <implicit.h>
      DIMENSION OMEGA1(NT1AMX), YIAJB(NT1AMX,NT1AMX)
      DIMENSION XIAJB(NT1AMX,NT1AMX)
      DIMENSION C2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO A = 1, NVIRT
      DO I = 1, NRHFT
        NAI = NVIRT*(I-1) + A
        DO D = 1, NVIRT
        DO L = 1, NRHFT
          NDI = NVIRT*(I-1) + D
          NDL = NVIRT*(L-1) + D
          NAL = NVIRT*(L-1) + A
          DO E = 1, NVIRT
          DO M = 1, NRHFT
            NEM = NVIRT*(M-1) + E
            DO F = 1, NVIRT
            DO N = 1, NRHFT
              NEN = NVIRT*(N-1) + E
              NFN = NVIRT*(N-1) + F
              NFM = NVIRT*(M-1) + F
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *                    + T3AM(NDL,NEM,NFN)
     *                    * C2AM(NFM,NDI)*YIAJB(NEN,NAL)
     *                    - T3AM(NDL,NEN,NFM)
     *                    * C2AM(NFM,NDI)*XIAJB(NEN,NAL)
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck t3_left2 */
      SUBROUTINE T3_LEFT2(OMEGA1,YIAJB,XIAJB,C2AM,T3AM)
C
C
C
#include <implicit.h>
      DIMENSION OMEGA1(NT1AMX), YIAJB(NT1AMX,NT1AMX)
      DIMENSION XIAJB(NT1AMX,NT1AMX)
      DIMENSION C2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO A = 1, NVIRT
      DO I = 1, NRHFT
        NAI = NVIRT*(I-1) + A
        DO D = 1, NVIRT
        DO L = 1, NRHFT
          NDL = NVIRT*(L-1) + D
          DO E = 1, NVIRT
          DO M = 1, NRHFT
            NAM = NVIRT*(M-1) + A
            NEI = NVIRT*(I-1) + E
            NEM = NVIRT*(M-1) + E
            DO F = 1, NVIRT
            DO N = 1, NRHFT
              NDN = NVIRT*(N-1) + D
              NFL = NVIRT*(L-1) + F
              NFN = NVIRT*(N-1) + F
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *                    + T3AM(NDL,NEM,NFN)
     *                    * C2AM(NAM,NDN)*YIAJB(NEI,NFL)
     *                    - T3AM(NDN,NEM,NFL)
     *                    * C2AM(NAM,NDN)*XIAJB(NEI,NFL)
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck t3_left3 */
      SUBROUTINE T3_LEFT3(OMEGA1,XIAJB,C2AM,T3AM)
C
C     Note : XIAJB is L and not g.
C
C
#include <implicit.h>
      DIMENSION OMEGA1(NT1AMX), XIAJB(NT1AMX,NT1AMX)
      DIMENSION C2AM(NT1AMX,NT1AMX), T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO A = 1, NVIRT
      DO I = 1, NRHFT
        NAI = NVIRT*(I-1) + A
        DO D = 1, NVIRT
        DO L = 1, NRHFT
          NDL = NVIRT*(L-1) + D
          DO E = 1, NVIRT
          DO M = 1, NRHFT
            NEM = NVIRT*(M-1) + E
            DO F = 1, NVIRT
            DO N = 1, NRHFT
              NEN = NVIRT*(N-1) + E
              NFM = NVIRT*(M-1) + F
              NFN = NVIRT*(N-1) + F
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *                    + (T3AM(NDL,NEM,NFN)-T3AM(NDL,NEN,NFM))
     *                    * C2AM(NEM,NDL)*XIAJB(NFN,NAI)
COMMENT COMMENT
COMMENT COMMENT
C      if (abs((T3AM(NDL,NEM,NFN)-T3AM(NDL,NEN,NFM))
C     *        *C2AM(NEM,NDL)*XIAJB(NFN,NAI)) .gt. 1.0D-12) then
C          write(lupri,'(A,8i3)') 'Contri from A,I,D,L,E,M,F,N : ',
C     *                                        a,i,d,l,e,m,f,n
C          write(lupri,*) 'T3AM, T3AM2, C2AM, INT : ',
C     *            T3AM(NDL,NEM,NFN),T3AM(NDL,NEN,NFM),
C     *            C2AM(NEM,NDL),XIAJB(NFN,NAI)
C      endif
COMMENT COMMENT
COMMENT COMMENT
            ENDDO
            ENDDO
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C  /* Deck cc3_tran3 */
      SUBROUTINE CC3_TRAN3(XINT1,XINT2,XLAMDP,XLAMDH,AOINT,IDEL)
C
C
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
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
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)
     *               + AOINT(NAB,G)*XLAMDP(A,NRHFT+C)*XLAMDH(IB,K)
     *                      *XLAMDP(G,NRHFT+B)*XLAMDH(IDEL,NRHFT+D)
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
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)
     *                  + AOINT(NAB,G)*XLAMDP(A,NRHFT+C)*XLAMDH(IB,K)
     *                                *XLAMDP(G,L)*XLAMDH(IDEL,J)
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
C  /* Deck ccsdt_l3am */
      SUBROUTINE CCSDT_L3AM(T3BAR,XINT1,XINT2,XIAJB,FOCK,T1AM,
     *                      T2AM,SCR1,FOCKD)
C
C
C
#include <implicit.h>
C
      PARAMETER (TWO = 2.0D0)
C
      DIMENSION XINT1(NT1AMX, NVIRT,NVIRT), XINT2(NT1AMX, NRHFT,NRHFT)
      DIMENSION XIAJB(NT1AMX, NT1AMX), FOCK(NORBT,NORBT)
      DIMENSION T3BAR(NT1AMX, NT1AMX,NT1AMX), SCR1(NT1AMX), FOCKD(*)
      DIMENSION T1AM(NT1AMX), T2AM(NT1AMX,NT1AMX)
      LOGICAL LDBG1, LDBG2
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 50 NI = 1,NRHFT
         DO 60 NA = 1,NVIRT
            NAI = NVIRT*(NI-1) + NA
            SCR1(NAI) = FOCKD(NRHFT+NA) - FOCKD(NI)
   60    CONTINUE
   50 CONTINUE
C
      DO 100 NK = 1,NRHFT
      DO 105 NC = 1,NVIRT
C
         NCK = NVIRT*(NK-1) + NC
C
         DO 110 NJ = 1,NRHFT
            NCJ = NVIRT*(NJ-1) + NC
            DO 120 NB = 1,NVIRT
C
               NBJ = NVIRT*(NJ-1) + NB
               NBK = NVIRT*(NK-1) + NB
C
               DO 130 NI = 1,NRHFT
                  NBI = NVIRT*(NI-1) + NB
                  NCI = NVIRT*(NI-1) + NC
               DO 135 NA = 1,NVIRT
C
                  NAI = NVIRT*(NI-1) + NA
                  NAJ = NVIRT*(NJ-1) + NA
                  NAK = NVIRT*(NK-1) + NA
C
                  AIBJCK = 0.0D0
C
          if (.true.) then
C     T1* = TWO T1 => Factor of two
                  AIBJCK = AIBJCK + XIAJB(NBJ,NCK)*T1AM(NAI)
COMMENT COMMENT
COMMENT COMMENT
C            if (abs(XIAJB(NBJ,NCK)*T1AM(NAI)) .gt. 1.0D-10) then
C               write(lupri,*) 'Contribution from 1.st term : ',
C     *                             XIAJB(NBJ,NCK)*T1AM(NAI)
C               write(lupri,*) 'With NAI, NBJ, NCK : ',NAI,NBJ,NCK
C            endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
C     T1* = TWO T1 => Factor of two
                  AIBJCK = AIBJCK - XIAJB(NBJ,NCI)*T1AM(NAK)
COMMENT COMMENT
COMMENT COMMENT
C             if (abs(XIAJB(NBJ,NCK)*T1AM(NAI)) .gt. 1.0D-10) then
C                write(lupri,*) 'Contribution from 2.nd term : ',
C     *                              XIAJB(NBJ,NCK)*T1AM(NAI)
C                write(lupri,*) 'With NAI, NBJ, NCK : ',NAI,NBJ,NCK
C           endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
                  AIBJCK = AIBJCK + FOCK(NK,NRHFT+NC)*
     *                           T2AM(NAI,NBJ)
          endif
C
          if (.true.) then
                  AIBJCK = AIBJCK - FOCK(NJ,NRHFT+NC)*
     *                           T2AM(NAI,NBK)
          endif
COMMENT COMMENT
COMMENT COMMENT
C       if (abs(FOCK(NK,NRHFT+NC)*T2AM(NAI,NBJ)) .gt. 1.0D-12) then
C           write(lupri,*) 'Fock-1 : Fock, L2 : ',
C     *              FOCK(NK,NRHFT+NC),T2AM(NAI,NBJ)
C       endif
C       if (abs(FOCK(NJ,NRHFT+NC)*T2AM(NAI,NBK)) .gt. 1.0D-12) then
C           write(lupri,*) 'Fock-2 : Fock, L2 : ',
C     *              FOCK(NJ,NRHFT+NC),T2AM(NAI,NBK)
C       endif
COMMENT COMMENT
COMMENT COMMENT
C
                  DO 140 ND = 1,NVIRT
C
                     NDI = NVIRT*(NI-1) + ND
                     NDJ = NVIRT*(NJ-1) + ND
                     NDK = NVIRT*(NK-1) + ND
C
          if (.true.) then
                     AIBJCK = AIBJCK + 
     *                        (TWO*XINT1(NCK,ND,NB)-XINT1(NBK,ND,NC))*
     *                        T2AM(NDJ,NAI)
COMMENT     *                        (TWO*XINT1(NCK,NB,ND)-XINT1(NBK,NC,ND))*
COMMENT     *                        T2AM(NDJ,NAI)
C
COMMENT COMMENT
COMMENT COMMENT
C              if (abs((TWO*XINT1(NCK,NB,ND)-XINT1(NBK,NC,ND))*
C     *                        T2AM(NDJ,NAI) .gt. 1.0D-10)) then
C                   write(lupri,*) 'Contribution from 5.th term :',
C     *                    TWO*(TWO*XINT1(NCK,NB,ND)-XINT1(NBK,NC,ND))*
C     *                        T2AM(NDJ,NAI)
C                   write(lupri,*) '2*XINT1, XINT1, T2AM ',
C     *                    TWO*XINT1(NCK,NB,ND),XINT1(NBK,NC,ND),
C     *                        T2AM(NDJ,NAI)
C              endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
                     AIBJCK = AIBJCK - 
     *                        XINT1(NBI,ND,NC)*
     *                        T2AM(NDK,NAJ)
COMMENT     *                        XINT1(NBI,NC,ND)*
COMMENT     *                        T2AM(NDK,NAJ)
COMMENT COMMENT
COMMENT COMMENT
C              if (abs(XINT1(NBI,NC,ND)*
C     *                        T2AM(NDK,NAJ))
C     *                          .gt. 1.0D-10) then
C                   write(lupri,*) 'Contribution from 6.th term :',
C     *                        TWO*XINT1(NBI,NC,ND)*
C     *                        T2AM(NDK,NAJ)
C              endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
  140             CONTINUE
C
                  DO 150 NL = 1,NRHFT
C
                     NAL = NVIRT*(NL-1) + NA
                     NBL = NVIRT*(NL-1) + NB
C
          if (.true.) then
                     AIBJCK = AIBJCK - 
     *                        (TWO*XINT2(NCK,NJ,NL)-XINT2(NCJ,NK,NL))*
     *                        T2AM(NBL,NAI)
COMMENT     *                        (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))*
COMMENT     *                        T2AM(NBL,NAI)
COMMENT COMMENT
COMMENT COMMENT
C         if (abs(TWO*
C     *           (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))*
C     *           T2AM(NBL,NAI)) .gt. 1.0D-10) then
CC
C            write(lupri,*) 'Contribution from L_occ term :',TWO*
C     *                 (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))*
C     *                 T2AM(NBL,NAI)
C            write(lupri,*) 'T2AM, L : ',
C     *                 T2AM(NBL,NAI),
C     *                 (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))
C         endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
                     AIBJCK = AIBJCK + 
     *                        XINT2(NCJ,NI,NL)*
     *                        T2AM(NBK,NAL)
COMMENT     *                        XINT2(NCJ,NL,NI)*
COMMENT     *                        T2AM(NBK,NAL)
COMMENT COMMENT
COMMENT COMMENT
C         if (abs(XINT2(NCJ,NL,NI)*T2AM(NBK,NAL)) .gt. 1.0D-12) then
CC
C            write(lupri,*) 'Contribution from g_occ term :',
C     *                 XINT2(NCJ,NL,NI)*T2AM(NBK,NAL)
C            write(lupri,*) 'T2AM, g : ',
C     *                 T2AM(NBL,NAI),XINT2(NCJ,NL,NI)
C         endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
  150             CONTINUE
C
C
                  AIBJCK = AIBJCK/(SCR1(NAI) + SCR1(NBJ) + SCR1(NCK))
C
COMMENT COMMENT
COMMENT COMMENT
C       if (abs(aibjck).gt.(1.0D-10/(scr1(nai)+scr1(nbj)+scr1(nck)))) 
C     *                                                              then
C          write(lupri,'(1X,A,I3,I3,I3,I3,I3,I3)') 
C     *                  'Contribution for A,B,C,I,J,K = ',
C     *                                 NA,NB,NC,NI,NJ,NK
C          write(lupri,'(A,3I3,e20.10)') 'NAI,NBJ,NCK,AIBJCK = ',
C     *                                   NAI,NBJ,NCK,AIBJCK
CC          write(lupri,*) 'SCR1(NAI),SCR1(NBJ),SCR1(NCK), SUM = ',
CC     *                    SCR1(NAI),SCR1(NBJ),SCR1(NCK),
CC     *                    SCR1(NAI)+SCR1(NBJ)+SCR1(NCK)
CC          write(lupri,*) '       '
C       endif
COMMENT COMMENT
COMMENT COMMENT
                  T3BAR(NAI,NBJ,NCK) = T3BAR(NAI,NBJ,NCK) + AIBJCK
                  T3BAR(NAI,NCK,NBJ) = T3BAR(NAI,NCK,NBJ) + AIBJCK
                  T3BAR(NBJ,NAI,NCK) = T3BAR(NBJ,NAI,NCK) + AIBJCK
                  T3BAR(NCK,NAI,NBJ) = T3BAR(NCK,NAI,NBJ) + AIBJCK
                  T3BAR(NBJ,NCK,NAI) = T3BAR(NBJ,NCK,NAI) + AIBJCK
                  T3BAR(NCK,NBJ,NAI) = T3BAR(NCK,NBJ,NAI) + AIBJCK
C
  135          CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  105 CONTINUE
  100 CONTINUE
C
C------------------------------------------------------
C     Get rid of amplitudes that are not allowed.
C------------------------------------------------------
C
      do a = 1, nvirt
      do i = 1, nrhft
         nai = nvirt*(i-1) + a
         do b = 1, nvirt
            nbi = nvirt*(i-1) + b
            do c = 1, nvirt
               nci = nvirt*(i-1) + c
C
               t3bar(nai,nbi,nci) = 0.0D0
               t3bar(nai,nci,nbi) = 0.0D0
               t3bar(nbi,nai,nci) = 0.0D0
               t3bar(nbi,nci,nai) = 0.0D0
               t3bar(nci,nai,nbi) = 0.0D0
               t3bar(nci,nbi,nai) = 0.0D0
            enddo
         enddo
C
         do j = 1, nrhft
            naj = nvirt*(j-1) + a
            do k = 1, nrhft
               nak = nvirt*(k-1) + a
C
               t3bar(nai,naj,nak) = 0.0D0
               t3bar(nai,nak,naj) = 0.0D0
               t3bar(naj,nai,nak) = 0.0D0
               t3bar(naj,nak,nai) = 0.0D0
               t3bar(nak,nai,naj) = 0.0D0
               t3bar(nak,naj,nai) = 0.0D0
            enddo
         enddo
C
      enddo
      enddo
C
      RETURN
      END
C  /* Deck khtcme */
      SUBROUTINE KHTCME1S(XINT1,XINT3,XINT2,XINT4)
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT3(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XINT4(NT1AMX,NRHFT,NRHFT)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (TWO = 2.0D0)
C
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nb = 1, nvirt
               nbk = nvirt*(nk-1) + nb
               do nd = 1, nvirt
C
               if (abs(xint1(nck,nb,nd)) .gt. 1.0D-12) then
                  write(lupri,*) 'g-S(',nck,nb,nd,') = ',
     *                     xint1(nck,nb,nd)
               endif
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '           '
C
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nb = 1, nvirt
               nbk = nvirt*(nk-1) + nb
               do nd = 1, nvirt
C
               xint3(nck,nb,nd) = two*xint1(nck,nb,nd)
     *                        -     xint1(nbk,nc,nd)
C
               enddo
            enddo
         enddo
      enddo
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nb = 1, nvirt
               nbk = nvirt*(nk-1) + nb
               do nd = 1, nvirt
C
               if (abs(xint3(nck,nb,nd)) .gt. 1.0D-12) then
                  write(lupri,*) 'L-S(',nck,nb,nd,') = ',
     *                     xint3(nck,nb,nd)
               endif
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '           '
C
C     XINT2, XINT4
C
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nl = 1, nrhft
               do nm = 1, nrhft
C
               if (abs(xint2(nck,nl,nm)) .gt. 1.0D-12) then
                  write(lupri,*) 'g-2S(',nck,nl,nm,') = ',
     *                     xint2(nck,nl,nm)
               endif
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '           '
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nl = 1, nrhft
               do nm = 1, nrhft
                  ncm = nvirt*(nm-1) + nc
C
                  xint4(nck,nl,nm) = two*xint2(nck,nl,nm)
     *                             -     xint2(ncm,nl,nk)
               enddo
            enddo
         enddo
      enddo
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nl = 1, nrhft
               do nm = 1, nrhft
C
               if (abs(xint4(nck,nl,nm)) .gt. 1.0D-12) then
                  write(lupri,*) 'L-2S(',nck,nl,nm,') = ',
     *                     xint4(nck,nl,nm)
               endif
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '           '
C
      RETURN
      END
C  /* Deck khtcme */
      SUBROUTINE KHTCME1T(XINT1,XINT3,XINT2,XINT4)
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT3(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XINT4(NT1AMX,NRHFT,NRHFT)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (TWO = 2.0D0)
C
C
C     XINT1, XINT3
C
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nb = 1, nvirt
               nbk = nvirt*(nk-1) + nb
               do nd = 1, nvirt
C
               if (abs(xint1(nck,nb,nd)) .gt. 1.0D-12) then
                  write(lupri,*) 'g-1T(',nck,nb,nd,') = ',
     *                     xint1(nck,nb,nd)
               endif
C
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '           '
C
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nb = 1, nvirt
               nbk = nvirt*(nk-1) + nb
               do nd = 1, nvirt
               ndk = nvirt*(nk-1) + nd
C
               xint3(nck,nb,nd) = two*xint1(nck,nb,nd)
     *                        -     xint1(ndk,nb,nc)
C
               enddo
            enddo
         enddo
      enddo
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nb = 1, nvirt
               nbk = nvirt*(nk-1) + nb
               do nd = 1, nvirt
C
               if (abs(xint3(nck,nb,nd)) .gt. 1.0D-12) then
                  write(lupri,*) 'L-1T(',nck,nb,nd,') = ',
     *                     xint3(nck,nb,nd)
               endif
C
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '           '
C
C     XINT2, XINT4
C
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nl = 1, nrhft
               ncl = nvirt*(nl-1) + nc
               do nm = 1, nrhft
                  ncm = nvirt*(nm-1) + nc
C
               if (abs(xint2(nck,nl,nm)) .gt. 1.0D-12) then
                  write(lupri,*) 'g-2T(',nck,nl,nm,') = ',
     *                     xint2(nck,nl,nm)
               endif
C
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '              '
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nl = 1, nrhft
               ncl = nvirt*(nl-1) + nc
               do nm = 1, nrhft
                  ncm = nvirt*(nm-1) + nc
C
                  xint4(nck,nl,nm) = two*xint2(nck,nl,nm)
     *                             -     xint2(ncl,nk,nm)
C
               enddo
            enddo
         enddo
      enddo
      do nc = 1, nvirt
         do nk = 1, nrhft
            nck = nvirt*(nk-1) + nc
            do nl = 1, nrhft
               do nm = 1, nrhft
C
               if (abs(xint4(nck,nl,nm)) .gt. 1.0D-12) then
                  write(lupri,*) 'L-2T(',nck,nl,nm,') = ',
     *                     xint4(nck,nl,nm)
               endif
C
               enddo
            enddo
         enddo
      enddo
      write(lupri,*) '              '
C
      RETURN
      END
C  /* Deck cc3__l3_omega */
      SUBROUTINE CC3_L3_OMEGA2_NODDY(OMEGA2,L3AM,XINT1,XINT2)
C
C     Kasper Hald, Spring 2002.
C
C     Calculate the doubles of the L3 part of CC3 left hand side doubles
C
#include <implicit.h>
C
      INTEGER INDEX, NAI, NBJ, NBK, NAK, NCK, NCJ, NCI, NDJ, NDK, NDI
      INTEGER NBL, NCL, NAL, NAJ, NBI, AI
C
#if defined (SYS_CRAY)
      REAL OMEGA2(NT1AMX,NT1AMX)
      REAL L3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL XINT1(NT1AMX,NVIRT,NVIRT)
      REAL XINT2(NT1AMX,NRHFT,NRHFT)
      REAL ONE, TWO, XAIBJ
#else
      DOUBLE PRECISION OMEGA2(NT1AMX,NT1AMX)
      DOUBLE PRECISION L3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION ONE, TWO, XAIBJ
#endif
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
C
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C===================================================
C     Calculate the L3 contributions to omega2
C===================================================
C
      CALL DZERO(OMEGA2,NT1AMX*NT1AMX)
C
      DO NI = 1,NRHFT
         DO NA = 1,NVIRT
            NAI = NVIRT*(NI-1) + NA
C
            DO NJ = 1,NRHFT
               DO NB = 1,NVIRT
                  NBJ = NVIRT*(NJ-1) + NB
C
                  DO NK = 1,NRHFT
                     NBK = NVIRT*(NK-1) + NB
                     NAK = NVIRT*(NK-1) + NA
                     DO NC = 1,NVIRT
C
                        NCK = NVIRT*(NK-1) + NC
                        NCJ = NVIRT*(NJ-1) + NC
                        NCI = NVIRT*(NI-1) + NC
C
                        DO ND = 1,NVIRT
C
                           NDJ = NVIRT*(NJ-1) + ND
                           NDK = NVIRT*(NK-1) + ND
                           NDI = NVIRT*(NI-1) + ND
C
                if (.true.) then
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *              L3AM(NBJ,NCI,NDK)*XINT1(NDK,NC,NA)
C
                endif
C
                        ENDDO
C
                        DO NL = 1,NRHFT
C
                           NBL = NVIRT*(NL-1) + NB
                           NCL = NVIRT*(NL-1) + NC
                           NAL = NVIRT*(NL-1) + NA
C
                 if (.true.) then
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI)
     *                     - L3AM(NBJ,NAK,NCL)*XINT2(NCL,NI,NK)
                 endif
C
                        ENDDO
C
                     ENDDO
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------------------------------------
C     Make P^{ab}_{ij} for the T3BAR contributions.
C----------------------------------------------------
C
      DO NAI = 1,NT1AMX
         DO NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
         ENDDO
      ENDDO
C
C=====================================
C     Write out all results
C=====================================
C
C      write(lupri,*) '        '
C      do nai = 1, nt1amx
C         do nbj = 1, nt1amx
C            if (abs(omega2(nai,nbj)) .gt. 1.0D-12) then
C               write(lupri,*) 'Omega2_noddy(',nai,nbj,') = ',
C     *                         omega2(nai,nbj)
C            endif
C         enddo
C      enddo
C      write(lupri,*) '        '
C
      RETURN
      END
C  /* Deck cc3__l3_omega1_noddy */
      SUBROUTINE CC3_L3_OMEGA1_NODDY(OMEGA1,L3AM,
     *                                XOOOO,XOVVO,XOOVV,XVVVV,T2AM)
C
C     Kasper Hald, Spring 2002.
C
C     Calculate the doubles of the L3 part of CC3 left hand side singles
C
#include <implicit.h>
C
      INTEGER INDEX, NAI, NBJ, NBK, NAK, NCK, NCJ, NCI, NDJ, NDK, NDI
      INTEGER NBL, NCL, NAL, NAJ, NBI, AI
C
#if defined (SYS_CRAY)
      REAL OMEGA1(NT1AMX)
      REAL L3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL XOOOO(NRHFT,NRHFT,NRHFT,NRHFT)
      REAL XOVVO(NRHFT,NVIRT,NVIRT,NRHFT)
      REAL XOOVV(NRHFT,NRHFT,NVIRT,NVIRT)
      REAL XVVVV(NVIRT,NVIRT,NVIRT,NVIRT)
      REAL T2AM(NT1AMX,NT1AMX)
      REAL ONE, TWO, XAIBJ
#else
      DOUBLE PRECISION OMEGA1(NT1AMX)
      DOUBLE PRECISION L3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XOOOO(NRHFT,NRHFT,NRHFT,NRHFT)
      DOUBLE PRECISION XOVVO(NRHFT,NVIRT,NVIRT,NRHFT)
      DOUBLE PRECISION XOOVV(NRHFT,NRHFT,NVIRT,NVIRT)
      DOUBLE PRECISION XVVVV(NVIRT,NVIRT,NVIRT,NVIRT)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION ONE, TWO, XAIBJ
#endif
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
C
#include <priunit.h>
#include <inforb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C===================================================
C     Calculate the L3BAR contributions to omega1
C===================================================
C
C-------------------------------------------
C     4 occupied integrals in integrals
C-------------------------------------------
C
      if (.true.) then
      DO NA = 1, NVIRT
       DO NB = 1, NVIRT
        DO NC = 1, NVIRT
         DO NI = 1, NRHFT
          NAI = NVIRT*(NI-1) + NA
          DO NJ = 1, NRHFT
           NBJ = NVIRT*(NJ-1) + NB
           DO NK = 1, NRHFT
            NCK = NVIRT*(NK-1) + NC
            DO NL = 1, NRHFT
             NAL = NVIRT*(NL-1) + NA
             DO NM = 1, NRHFT
                NCM = NVIRT*(NM-1) + NC
                OMEGA1(NAI) = OMEGA1(NAI)
     *           + L3AM(NBJ,NAL,NCM)*T2AM(NBJ,NCK)*XOOOO(NI,NL,NK,NM)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      endif
C
C-------------------------------------------
C     4 virtual integrals in integrals
C-------------------------------------------
C
      if (.true.) then
      DO NA = 1, NVIRT
       DO NB = 1, NVIRT
        DO NC = 1, NVIRT
         DO ND = 1, NVIRT
          DO NE = 1, NVIRT
           DO NI = 1, NRHFT
            NAI = NVIRT*(NI-1) + NA
            NDI = NVIRT*(NI-1) + ND
            DO NJ = 1, NRHFT
             NBJ = NVIRT*(NJ-1) + NB
             DO NK = 1, NRHFT
              NEK = NVIRT*(NK-1) + NE
              NCK = NVIRT*(NK-1) + NC
C
              OMEGA1(NAI) = OMEGA1(NAI)
     *          + L3AM(NBJ,NDI,NEK)*T2AM(NBJ,NCK)*XVVVV(ND,NA,NE,NC)
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      endif
C
C--------------------------------------------------------
C     2 terms with g_{ovvo} and 2 terms with g_{oovv}
C--------------------------------------------------------
C
      if (.true.) then
      DO NA = 1, NVIRT
       DO NB = 1, NVIRT
        DO NC = 1, NVIRT
         DO ND = 1, NVIRT
          DO NI = 1, NRHFT
           NAI = NVIRT*(NI-1) + NA
           NCI = NVIRT*(NI-1) + NC
           NDI = NVIRT*(NI-1) + ND
           DO NJ = 1, NRHFT
            NBJ = NVIRT*(NJ-1) + NB
            DO NK = 1, NRHFT
             NAK = NVIRT*(NK-1) + NA
             NCK = NVIRT*(NK-1) + NC
             NDK = NVIRT*(NK-1) + ND
             DO NL = 1, NRHFT
              NAL = NVIRT*(NL-1) + NA
              NCL = NVIRT*(NL-1) + NC
              NDL = NVIRT*(NL-1) + ND
C
              if (.true.) then
C - L3^{daf}_{lmn} * t^{de}_{lm} * g_{iefn}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NAK,NDL)*T2AM(NBJ,NCK)*XOVVO(NI,NC,ND,NL)
              endif
C
              if (.true.) then
C - L3^{daf}_{lnm} * t^{de}_{lm} * g_{infe}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NAL,NDK)*T2AM(NBJ,NCK)*XOOVV(NI,NL,ND,NC)
              endif
C
              if (.true.) then
C - L3^{def}_{lin} * t^{de}_{lm} * g_{mafn}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NCI,NDL)*T2AM(NBJ,NCK)*XOVVO(NK,NA,ND,NL)
              endif
C
              if (.true.) then
C - L3^{def}_{lni} * t^{de}_{lm} * g_{mnfa}
              OMEGA1(NAI) = OMEGA1(NAI)
     *       - L3AM(NBJ,NCL,NDI)*T2AM(NBJ,NCK)*XOOVV(NK,NL,ND,NA)
              endif
C
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      endif
C
C=====================================
C     Write out all results
C=====================================
C
      write(lupri,*) '        '
      do nai = 1, nt1amx
         if (abs(omega1(nai)) .gt. 1.0D-12) then
            write(lupri,*) 'Omega1_noddy(',nai,')             = ',
     *                      omega1(nai)
         endif
      enddo
      write(lupri,*) '        '
C
      RETURN
      END
C  /* Deck ccfop_tran1 */
      SUBROUTINE CCFOP_TRAN1(XINT1,XINT2,XINT3,XINT4,
     *                       XLAMDP,XLAMDH,AOINT,IDEL)
C
C
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
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
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
COMMENT COMMENT
COMMENT COMMENT
      write(lupri,*) 'Noddy integral code : NBAST, NRHFT, NVIRT : ',
     *                                      nbast, nrhft, nvirt
COMMENT COMMENT
COMMENT COMMENT
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
     *                  + AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NRHFT+NC)
     *                         *XLAMDP(G,NRHFT+ND)*XLAMDH(IDEL,NL)
COMMENT
COMMENT
C       if (abs(AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NRHFT+NC)*
C     *         XLAMDP(G,NRHFT+ND)*XLAMDH(IDEL,NL)) .gt. 1.0D-12) then
C              write(lupri,'(A,4I3,A,e20.10)') 
C     *         'Contri to xint1(',nk,nc,nd,nl,') = ',xint1(nk,nc,nd,nl)
C       endif
COMMENT
COMMENT
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
     *                  + AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NL)
     *                         *XLAMDP(G,NRHFT+NC)*XLAMDH(IDEL,NRHFT+ND)
COMMENT
COMMENT
C       if (abs(AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NL)*
C     *      XLAMDP(G,NRHFT+NC)*XLAMDH(IDEL,NRHFT+ND)) .gt. 1.0D-12) then
C              write(lupri,'(A,4I3,A,e20.10)') 
C     *         'Contri to xint2(',nk,nl,nc,nd,') -> ',xint2(nk,nl,nc,nd)
C       endif
COMMENT
COMMENT
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
     *                  + AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NL)
     *                         *XLAMDP(G,NM)*XLAMDH(IDEL,NN)
C
COMMENT
COMMENT
C       if (abs(AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NL)*
C     *      XLAMDP(G,NM)*XLAMDH(IDEL,NN)) .gt. 1.0D-12) then
C              write(lupri,'(A,8I4)') 'A, IB, G, IDEL : ',
C     *                                a, ib, g, idel
C              write(lupri,'(A,8I4)') 'K, L, M, N     : ',
C     *                               nk,nl,nm,nn
C              write(lupri,'(A,3e20.10)') '3. AOINT, LAM1, LAM2 : ',
C     *              AOINT(NAB,G),XLAMDP(A,NK),XLAMDH(IB,NL)
C              write(lupri,'(A,3e20.10)') '3. Cont, LAM3, LAM4  : ',
C     *                   AOINT(NAB,G)*XLAMDP(A,NK)*XLAMDH(IB,NL)
C     *                         *XLAMDP(G,NM)*XLAMDH(IDEL,NN),
C     *                          XLAMDP(G,NM),XLAMDH(IDEL,NN)
C              write(lupri,*) '     '
C       endif
COMMENT
COMMENT
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
     *           + AOINT(NAB,G)*XLAMDP(A,NRHFT+NC)*XLAMDH(IB,NRHFT+ND)
     *                         *XLAMDP(G,NRHFT+NE)*XLAMDH(IDEL,NRHFT+NF)
C
COMMENT
COMMENT
C
C       if (abs(AOINT(NAB,G)*XLAMDP(A,NRHFT+NC)*XLAMDH(IB,NRHFT+ND)*
C     *      XLAMDP(G,NRHFT+NE)*XLAMDH(IDEL,NRHFT+NF)) .gt. 1.0D-12) then
C              write(lupri,*) 'NAB, G, A, IB, IDEL',NAB,G,A,IB,IDEL
C              write(lupri,'(A,4I3,A,3e18.10)') 
C     *         'Contri to xint4(',nc,nd,ne,nf,') : ',
C     *               AOINT(NAB,G),XLAMDP(A,NRHFT+NC),XLAMDH(IB,NRHFT+ND)
C              write(lupri,'(A,3e18.10)') 
C     *         'Contri to xint4 continued     : ',
C     *               XLAMDP(G,NRHFT+NE),XLAMDH(IDEL,NRHFT+NF),
C     *               AOINT(NAB,G)*XLAMDP(A,NRHFT+NC)*XLAMDH(IB,NRHFT+ND)
C     *              *XLAMDP(G,NRHFT+NE)*XLAMDH(IDEL,NRHFT+NF)
C              write(lupri,*) '       '
C       endif
COMMENT
COMMENT
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
C  /* Deck new_int_output */
      SUBROUTINE NEW_INT_OUTPUT(XINT1,XINT2,XINT3,XINT4)
C
#include <implicit.h>
      DIMENSION XINT1(NRHFT,NVIRT,NVIRT,NRHFT)
      DIMENSION XINT2(NRHFT,NRHFT,NVIRT,NVIRT)
      DIMENSION XINT3(NRHFT,NRHFT,NRHFT,NRHFT)
      DIMENSION XINT4(NVIRT,NVIRT,NVIRT,NVIRT)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL AROUND('ENTERED NEW_INT_OUTPUT')
C
      DO NK = 1,NRHFT
         DO NC = 1,NVIRT
            DO ND = 1,NVIRT
               DO NL = 1,NRHFT
                  if (abs(xint1(nk,nc,nd,nl)) .gt. 1.0D-12) then
               write(lupri,1) 'g-ovvo(',nk,',',nc,',',nd,',',nl,') = ',
     *                 xint1(nk,nc,nd,nl)
                  endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      write(lupri,*) '       '
C
      DO NK = 1,NRHFT
         DO NL = 1,NRHFT
            DO NC = 1,NVIRT
               DO ND = 1,NVIRT
                  if (abs(xint2(nk,nl,nc,nd)) .gt. 1.0D-12) then
               write(lupri,1) 'g-oovv(',nk,',',nl,',',nc,',',nd,') = ',
     *                 xint2(nk,nl,nc,nd)
                  endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      write(lupri,*) '       '
C
      DO NK = 1,NRHFT
         DO NL = 1,NRHFT
            DO NM = 1,NRHFT
               DO NN = 1,NRHFT
                  if (abs(xint3(nk,nl,nm,nn)) .gt. 1.0D-12) then
               write(lupri,1) 'g-oooo(',nk,',',nl,',',nm,',',nn,') = ',
     *                 xint3(nk,nl,nm,nn)
                  endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      write(lupri,*) '       '
C
      DO NC = 1,NVIRT
         DO ND = 1,NVIRT
            DO NE = 1,NVIRT
               DO NF = 1,NVIRT
                  if (abs(xint4(nc,nd,ne,nf)) .gt. 1.0D-12) then
               write(lupri,1) 'g-vvvv(',nc,',',nd,',',ne,',',nf,') = ',
     *                 xint4(nc,nd,ne,nf)
                  endif
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      write(lupri,*) '       '
C
   1  FORMAT(1x,A,I3,A,I3,A,I3,A,I3,A,E20.10)
C
      RETURN
      END
C  /* Deck cc_lhtr_noddy2 */
      SUBROUTINE CC_LHTR_NODDY2(OMEGA1,OMEGA2,OMEGA3,T1AM,T2AM,T3AM,
     *                          C1AM,C2AM,C3AM,XLAMDP,XLAMDH,FOCK,
     *                          WORK,LWORK)
C
C     Written by Henrik Koch 19-Sep-1994
C
C     Version 1.0
C
C     Purpose:
C
C     Calculate the iterative triples corrections to the CCSD
C     equations.
C
C     NB! The T2 amplitudes are assumed to be a full square.
C
C
C     NB! It is assumed that the vectors are allocated in the following
C     order:
C           T1AM(*), OMEGA1(*), OMEGA2(*), T2AM(*), SCR(*), WRK(*).
C
C     IOPT = 0 No storage of T3, L3
C     IOPT = 1 Store T3 and L3 in 'T3AMPL' and 'L3AMPL'
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION OMEGA1(*),OMEGA2(*), OMEGA3(*)
      DIMENSION T1AM(*),T2AM(*), T3AM(*)
      DIMENSION C1AM(*), C2AM(*), C3AM(*)
      DIMENSION XLAMDP(*),XLAMDH(*), FOCK(*), WORK(LWORK)
C
#include <inforb.h>
#include <infind.h>
#include <blocks.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
COMMENT COMMENT
COMMENT COMMENT
#include <mxcent.h>
#include <eribuf.h>
COMMENT COMMENT
COMMENT COMMENT
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (NSYM .NE. 1) THEN
         CALL QUIT('No symmetry in this part yet')
      ENDIF
C
C------------------------
C     Dynamic Allocation.
C------------------------
C
      KCMO   = 1
      KFOCKD = KCMO   + NLAMDT
      KINT1  = KFOCKD + NORBT
      KINT2  = KINT1  + NT1AMX*NVIRT*NVIRT
      KINT1T = KINT2  + NRHFT*NRHFT*NT1AMX
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KINT1S = KINT2T + NRHFT*NRHFT*NT1AMX
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KINT3S = KINT2S + NRHFT*NRHFT*NT1AMX
      KINT4S = KINT3S + NT1AMX*NVIRT*NVIRT
      KOME1  = KINT4S + NRHFT*NRHFT*NT1AMX
      KOME2  = KOME1  + NT1AMX
      KOME22 = KOME2  + NT1AMX*NT1AMX
      KIAJB  = KOME22 + NT1AMX*NT1AMX
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KEND1  = KYIAJB + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
C
C     New Integrals.
      KIOVVO = KEND1
      KIOOVV = KIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KIOOOO = KIOOVV + NRHFT*NRHFT*NVIRT*NVIRT
      KIVVVV = KIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1  = KIVVVV + NVIRT*NVIRT*NVIRT*NVIRT
      LWRK1  = LWORK  - KEND1
C
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSD_TRIPLE')
      ENDIF
C
C--------------------------------
C     Initialize integral arrays.
C--------------------------------
C
      CALL DZERO(WORK(KINT1),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT3S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT4S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME22),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)
C
      CALL DZERO(WORK(KIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      IF (DIRECT) THEN
         NTOSYM = 1
         DTIME  = SECOND()
c        CALL HERDI1(WORK(KEND1),LWRK1,IPRINT)
         CALL QUIT('Direct not implemented')
         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            NTOT = MAXSHL
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed t2am.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
               DTIME  = SECOND()
c              CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,IPRINT)
               CALL QUIT('Direct not implemented')
               DTIME  = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
COMMENT COMMENT
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_FOP3')
               END IF
COMMENT COMMENT
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCSD_TRIPLE')
               ENDIF
C
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
C
C----------------------------------------------------
C              Calculate integrals needed in CCSD[T].
C----------------------------------------------------
C
               CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CC3_TRAN2(WORK(KIAJB),WORK(KYIAJB),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN3(WORK(KINT1S),WORK(KINT2S),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCFOP_TRAN1(WORK(KIOVVO),WORK(KIOOVV),
     *                          WORK(KIOOOO),WORK(KIVVVV),
     *                          XLAMDP,XLAMDH,WORK(KXINT),IDEL)
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C------------------------------------
C     Calculate the corrections from T3.
C------------------------------------
C
      CALL T3_LEFT1(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,T3AM)
C
      CALL T3_LEFT2(WORK(KOME1),WORK(KYIAJB),WORK(KIAJB),
     *              C2AM,T3AM)
C
      CALL T3_LEFT3(WORK(KOME1),WORK(KIAJB),C2AM,T3AM)
C
      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      ENDDO
C
C----------------------------------------------------------
C         Calculate the contributions to Omega2 from C3
C         The result is written out inside the routine.
C----------------------------------------------------------
C
          CALL DZERO(WORK(KOME1),NT1AMX)
          CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
          CALL DZERO(WORK(KOME22),NT1AMX*NT1AMX)
 
          CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),C3AM,
     *                             WORK(KIOOOO),WORK(KIOVVO),
     *                             WORK(KIOOVV),WORK(KIVVVV),T2AM)
 
          CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),C3AM,
     *                             WORK(KINT1S),WORK(KINT2S))
 
      DO 300 I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
  300 CONTINUE
C
C     DO I = 1,NT1AMX
C        DO J = 1,I
C           IJ = NT1AMX*(I-1) + J
C           NIJ = INDEX(I,J)
C           OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
C        ENDDO
C     ENDDO

      CALL DAXPY(NT1AMX*NT1AMX,1.0D0,WORK(KOME2),1,OMEGA2,1)
C
C---------------------------------------------------------
C     Calculate the contributions to omega3 from C1, C2
C---------------------------------------------------------
      do i = 1, nrhft
         work(kfockd-1+i) = 0.0D0
      enddo
      do i = 1, nvirt
         work(kfockd-1+nrhft+i) = 1.0D0/3.0D0
      enddo
C
      CALL CCSDT_L3AM(OMEGA3,WORK(KINT1T),WORK(KINT2T),
     *                WORK(KIAJB),FOCK,C1AM,C2AM,
     *                WORK(KSCR1),WORK(KFOCKD))
C
C----------------------------------------
C     Print norm of result vectors.
C----------------------------------------
C
      XNORM = DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
      write(lupri,*) 'Norm Omega1 at the end of noddy :',xnorm
      XNORM = DDOT(NT2AM(1),OMEGA2,1,OMEGA2,1)
      write(lupri,*) 'Norm Omega2 at the end of noddy :',xnorm
C
      RETURN
      END
