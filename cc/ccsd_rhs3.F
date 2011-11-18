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
C  /* Deck ccrhsn3 */
      SUBROUTINE CCRHSN3(WORK,LWORK)
!
!     Written by Kasper Hald late march 1999
!
!     Based on CCRHSN
!     Written by Henrik Koch 25-Sep-1993
!
!     Purpose:
!
!     Calculation of the triplet Coupled Cluster global int. using
!     AO-integrals directly from disk. No intermediate files are
!     used in this calculation. The intermediates are written to disk
!     even if RSPIM is false.
!
!
!     N.B. The routine is only tested for OMEGOR
!
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <second.h>
#include <r12int.h>
!
      INTEGER LWORK, ISYMTR, LUCSIM, LUDTIM, IERRCS, IERRDT
      INTEGER KOMG1, KOMG2, KLAMDP, KLAMIP, KLAMDH, KDENSI
      INTEGER KFOCK, KEMAT1, KEMAT2, KEND1, LWRK1, JSYM, ISYMH, IC
      INTEGER IF, KENDS2, LWRKS2, NTOSYM, KCCFB1, KINDXB, KT2AMT
      INTEGER KODCL1, KODCL2, KODBC1, KODBC2, KRDBC1, KRDBC2
      INTEGER KODPP1, KODPP2, KRDPP1, KRDPP2, KFREE, LFREE, KENDSV
      INTEGER LWRKSV, ICDEL1, ISYMD1, NTOT, ILLL, NUMDIS, KRECNR
      INTEGER IDEL, IDEL2, ISYMD, ISYDIS, KXINT, KEND2, LWRK2
      INTEGER ISYDEN, ICON, ISYMLH, IOPT, ISYM0
      INTEGER KDSRHF, KEND4, LWRK4, ISYMLP, IV, ISYM, LUBF
      INTEGER KGAMMA, ISYMBF, LUGAM, LUFCK, ISYMEI, LUE1, LUE2
      INTEGER KOFF, KT1AM, KT2AM
!
      INTEGER INDEXA(MXCORB)
!
      INTEGER IDUM
!
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL DTIME, DDOT, FACTC, FACTD, FAKE, FF, HALF, ONE
      REAL PHONEY, RHO1N, RHO2N, FACTE1, FACTE2
      REAL TIMALL, TIMBF, TIMC, TIMDM
      REAL TIMEI, TIMF, TIMFCK, TIMGAM, TIMHER1
      REAL TIMHER2, TIMLAM, TIMRDAO, TIMT2AO, TIMT2BT, TIMT2TR, TIMTRBT
      REAL TWO, XMHALF, XMONE, ZERO, XNORM
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION DTIME, DDOT, FACTC, FACTD
      DOUBLE PRECISION FAKE, FF, HALF, ONE, PHONEY
      DOUBLE PRECISION RHO1N, RHO2N, FACTE1, FACTE2
      DOUBLE PRECISION TIMALL, TIMBF, TIMC, TIMDM
      DOUBLE PRECISION TIMEI, TIMFCK
      DOUBLE PRECISION TIMGAM, TIMHER1, TIMHER2, TIMLAM
      DOUBLE PRECISION TIMRDAO, TIMT2AO, TIMT2BT, TIMT2TR, TIMTRBT
      DOUBLE PRECISION TWO, XMHALF, XMONE, ZERO, XNORM
#endif
!
      PARAMETER (ZERO= 0.0D00, HALF= 0.5D00, ONE= 1.0D00, TWO= 2.0D00)
      PARAMETER (XMHALF = -0.5D00, XMONE= -1.0D00 )
!
      LOGICAL FCKCON, ETRAN, ANTISYM
!
      CHARACTER  CSFIL*6, DTFIL*7
      CHARACTER*10 MODEL
!
!-----------------------------------------------------------
!     For energy calculation trial vector is totalsymmetric.
!-----------------------------------------------------------
!
      ISYMTR = 1
!
!----------------
!     Open files.
!----------------
!
      LUDTIM  = -1
      LUCSIM  = -1
      CSFIL   = 'PMAT_C'
      DTFIL   = 'PMAT_DT'
!
      IF (DUMPCD) THEN
         CALL WOPEN2(LUDTIM,DTFIL,64,0)
         CALL WOPEN2(LUCSIM,CSFIL,64,0)
!
      END IF
!
!----------------------------------
!     Initialize timing parameters.
!----------------------------------
!
      TIMALL  = SECOND()
      TIMBF   = 0.0D00
      TIMDM   = 0.0D00
      TIMEI   = 0.0D00
      TIMFCK  = 0.0D00
      TIMHER1 = 0.0D00
      TIMHER2 = 0.0D00
      TIMLAM  = 0.0D00
      TIMRDAO = 0.0D00
      TIMT2AO = 0.0D00
      TIMT2BT = 0.0D00
      TIMT2TR = 0.0D00
      TIMTRBT = 0.0D00
!
!---------------------------
!     Check inconsistencies.
!---------------------------
!
      IF (NEWGAM) THEN
         IF ((.NOT. DUMPCD) .OR. (.NOT. OMEGOR)) THEN
            WRITE(LUPRI,*) 'NEWGAM requires both DUMPCD and OMEGOR'
            CALL QUIT('ERROR: NEWGAM inconsistency')
         END IF
      END IF
!
!---------------------------------
!     Work space allocation no. 1.
!---------------------------------
!
      KT1AM  = 1
      KT2AM  = KT1AM  + NT1AMX
      KLAMDP = KT2AM  + NT2SQ(1)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
      LWRK1  = LWORK  - KEND1
!
      IF (LWRK1 .LT. NT2AMX) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1+NT2AMX,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHSN3')
      ENDIF
!
!-------------------------------
!     Prepare the t2-amplitudes.
!-------------------------------
!
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KEND1))
!
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMTR)
!
!----------------------------------
!     Calculate the lamda matrices.
!----------------------------------
!
      TIMLAM  = SECOND()
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)
      TIMLAM  = SECOND() - TIMLAM
!
!====================================================
!     Start the loop over distributions of integrals.
!====================================================
!
      KENDS2 = KEND1
      LWRKS2 = LWRK1
!
      IF (DIRECT) THEN
         DTIME  = SECOND()

         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                 WORK(KEND1),LWRK1,IPRERI)
           KEND1 = KFREE
           LWRK1 = LFREE
         END IF

         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
!
      KENDSV = KEND1
      LWRKSV = LWRK1
!
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
!
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
!
         DO 110 ILLL = 1,NTOT
!
!-----------------------------------------------------------------
!           If direct calculate the integrals and transposed t2am.
!-----------------------------------------------------------------
!
            IF (DIRECT) THEN
!
               KEND1 = KENDSV
               LWRK1 = LWRKSV
!
               DTIME  = SECOND()
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                       WORK(KODCL1),WORK(KODCL2),
     &                       WORK(KODBC1),WORK(KODBC2),
     &                       WORK(KRDBC1),WORK(KRDBC2),
     &                       WORK(KODPP1),WORK(KODPP2),
     &                       WORK(KRDPP1),WORK(KRDPP2),
     &                       WORK(KCCFB1),WORK(KINDXB),
     &                       WORK(KEND1), LWRK1,IPRERI)
               END IF
               DTIME   = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
!
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CCRHSN3')
               END IF
!
            ELSE
               NUMDIS = 1
            ENDIF
!
!-----------------------------------------------------
!           Loop over number of distributions in disk.
!-----------------------------------------------------
!
            DO 120 IDEL2 = 1,NUMDIS
!
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
!
               ISYDIS = MULD2H(ISYMD,ISYMOP)
!
!------------------------------------------
!              Adresses for C/D
!------------------------------------------
!
               IT2DEL(IDEL) = ICDEL1
               ICDEL1 = ICDEL1 + NT2BCD(ISYDIS)
!
!------------------------------------------
!              Work space allocation no. 2.
!------------------------------------------
!
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
!
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHSN3')
               ENDIF
!
!
!-----------------------------------------
!              Read in batch of integrals.
!-----------------------------------------
!
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
!
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
!
!------------------------------------------
!              Work space allocation no. 4.
!------------------------------------------
!
               KDSRHF = KEND2
               KEND4  = KDSRHF + NDSRHF(ISYMD)
               LWRK4  = LWORK  - KEND4
!
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND4,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHSN')
               ENDIF
!
!--------------------------------------------------------
!              Transform one index in the integral batch.
!--------------------------------------------------------
!
               DTIME   = SECOND()
               ISYMLP  = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                     ISYMLP,WORK(KEND4),LWRK4,ISYDIS)
               DTIME   = SECOND() - DTIME
               TIMTRBT = TIMTRBT + DTIME
!
!-----------------------------------------------------
!              Calculate (3)D
!-----------------------------------------------------
!
               IF (.NOT. CC2) THEN
!
               DTIME   = SECOND()
!
               FACTD = ONE
               ICON   = 5
               IV = 1
!
                  CALL CCRHS_D3(WORK(KXINT),WORK(KDSRHF),DUMMY,
     *                          WORK(KT2AM),ISYMTR,WORK(KLAMDP),DUMMY,
     *                          WORK(KLAMDH),WORK(KLAMDP),ISYMTR,
     *                          WORK(KLAMDH),ISYMTR,
     *                          DUMMY,WORK(KEND4),LWRK4,IDEL,
     *                          ISYMD,FACTD,ICON,LUDTIM,DTFIL,IV)
!
               DTIME   = SECOND() - DTIME
               TIMDM   = TIMDM  + DTIME
!
               ENDIF
!
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
!
!--------------------
!     Times
!--------------------
!
      TIMALL  = SECOND() - TIMALL
      IF ( IPRINT .GT. 2) THEN
         WRITE(LUPRI,9999) 'RHSN3 - TOT.', TIMALL
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,9999) 'CCRHS_D3   ', TIMDM
         WRITE(LUPRI,9999) 'HERDIS1    ', TIMHER1
         WRITE(LUPRI,9999) 'HERDIS2    ', TIMHER2
         WRITE(LUPRI,9999) 'CCRHS_LAM  ', TIMLAM
         WRITE(LUPRI,9999) 'CCRHS_RDAO ', TIMRDAO
      ENDIF
9999  FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')
!
!-----------------
!     Close files.
!-----------------
!
         CALL WCLOSE2(LUDTIM,DTFIL,'KEEP')
         CALL WCLOSE2(LUCSIM,CSFIL,'KEEP')
!
      RETURN
      END
