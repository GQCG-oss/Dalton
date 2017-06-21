C  /* Deck so_2ndDEcoef */
      SUBROUTINE SO_2NDDECOEF(OMEGA2,T2AM,T2AMSQ,WORK,
     *                           LWORK)
C      SUBROUTINE SO_2NDDECOEF(OMEGA1,OMEGA2,T1AM,T2AM,WORK,LWORK)
C
C     Orignial routine, CCRHSN, is written by Henrik Koch 25-Sep-1993
C
C     Purpose: calculate the 2nd order doubles coeffecients.
C
C     Calculation of the Coupled Cluster vector function using
C     AO-integrals directly from disk. 
C
C     It is assumed that the arrays are stored in the follwoing way:
C           T2AM, OMEGA2, T2AMSQ
C
C     List of Varialbels:
C       OMEGA2: contains the 2nd order doubles amplitudes.
C       T2AM:   Contains the 1st order doubles amplitudes.
C       T2AMSQ: Contains the squared 1st order doubles amplitudes.
C
C
C     Keyword: T2TCOR is assumed to be true, but have been commentet
C               out in case, that it is need, because of lack of mem.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccorb.h"
#include "ccisao.h"
#include "blocks.h"
#include "ccfield.h"
#include "ccsections.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "eribuf.h"
#include "ccnoddy.h"
#include "cbirea.h"
C#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (XMHALF = -0.5D0, XMONE= -1.0D0 )
      PARAMETER (ISYM0 = 1)
C
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION OMEGA2(NT2SQ(1)),T2AM(NT2AMX),T2AMSQ(NT2SQ(1))
      DIMENSION WORK(LWORK)
C
      CHARACTER CFIL*6,DFIL*6
      CHARACTER MODEL*10
C
      INTEGER IGLMRHS(8,8),IGLMVIS(8,8),NGLMDS(8),ICMO(8,8),NCMO(8),
     &        IMAIJM(8,8),NMAIJM(8),
     &        IMATIJM(8,8),NMATIJM(8),NGAMSM(8),IGAMSM(8,8),
     &        IRGIJS(8,8),NRGIJS(8),IR1BASM(8,8),NR1BASM(8),
     &        IR2BASM(8,8),NR2BASM,IR1XBASM(8,8),NR1XBASM(8),
     &        IR2XBASM(8,8),IMATF(8,8),NMATF(8),IMAKLM(8,8),NMAKLM(8)
      INTEGER NADP(8),IADP(8,8),NLAMDX(8),ILAMDX(8,8)
C
      CALL QENTER('SO_2NDDECOEF')
C
C-----------------------------------------------------------
C     For energy calculation trial vector is totalsymmetric.
C-----------------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C Is ISYMOP defines in SOPPA?
C
C----------------
C     Open files.
C----------------
C
      LUC = -1
      LUD = -1
      CFIL = 'PMAT_C'
      DFIL = 'PMAT_D'
C
      CALL WOPEN2(LUC,CFIL,64,0)
      CALL WOPEN2(LUD,DFIL,64,0)
C
C
C----------------------------------
C     Initialize timing parameters.
C----------------------------------
      WRITE(LUPRI, '(A)') 'TIMING PARAMETERS'
C
      TIMALL  = SECOND()
      TIMA    = 0.0D00
      TIMB    = 0.0D00
      TIMC    = 0.0D00
      TIMD    = 0.0D00
      TIMGAM  = 0.0D00
      TIMCMO  = 0.0D00
      TIMRDAO = 0.0D00
      TIMHER1 = 0.0D00
      TIMHER2 = 0.0D00
      TIMT2AO = 0.0D00
      TIMFCK  = 0.0D00
      TIMDM   = 0.0D00
      TIMFCKMO= 0.0D00
      TIMT2TR = 0.0D00
      TIMT2BT = 0.0D00
      TIMTRBT = 0.0D00
C
C-------------------------------------------------------------
C     Work space allocation no. 1.
C--------------------------------------------------------------
C
      KCMO = 1
      KGAMMA = KCMO + NLAMDT
      KEND1 = KGAMMA + NGAMMA(ISYMOP) !længde af kgamma mangler da? 
c
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in SO_2NDDECOEF')
      ENDIF
C
C----------------------------------
C     Calculate the MO coeffecients
C--------------------------------
C
      TIMCMO  = SECOND()
      CALL SO_GETMO(WORK(KCMO),NLAMDT,WORK(KEND1),LWRK1)
      TIMCMO  = SECOND() - TIMCMO
C
C-------------------------------
C     Square-up the t2-amplitudes.
C-------------------------------
C
      CALL CC_T2SQ(T2AM,T2AMSQ,ISYMTR)
C      CALL OUTPUT(T2AM,1,NT2AMX,1,1,NT2AMX,1,1,LUPRI)
C      CALL OUTPUT(T2AMSQ,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C
C
C-----------------------------------------
C     Construct the transposed amplitudes.
C-----------------------------------------
C
C      IF ((.NOT. DIRECT) .AND. T2TCOR) THEN
      IF (.NOT. DIRECT) THEN
C
         KT2AMT = KEND1
         KEND1  = KT2AMT + NT2SQ(1)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient core in CCRHSN')
         END IF
C
         JSYM = 1
         CALL DCOPY(NT2SQ(1),T2AMSQ,1,WORK(KT2AMT),1)
C               Copies T2AMSQ to WORK(KT2AMT)
         CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND1),LWRK1,JSYM)
C         CALL OUTPUT(T2AMSQ,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C         CCSD_T2TP transposes the ij index of the squaed T2-amplitudes.
C           They are put in WORK(KT2AMT)
      END IF
C
C-------------------------------
C     Initialize  OMEGA2
C-------------------------------
C
C      CALL DZERO(OMEGA2,NT2ORT(ISYMOP))
      CALL DZERO(OMEGA2,NT2SQ(1))
C      CALL DZERO(OMEGA2,2*NT2ORT(ISYMOP))
C
C      WRITE(LUPRI,*) 'OMEGA2'
C      CALL OUTPUT(OMEGA2,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C      FLUSH(LUPRI)
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
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
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         WRITE(LUPRI,*) 'ISYMD1', ISYMD1
         FLUSH(LUPRI)
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
            WRITE(LUPRI,*) 'ILLL:', ILLL
            FLUSH(LUPRI)
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed t2am.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRERI)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                        WORK(KODCL1),WORK(KODCL2),
     &                        WORK(KODBC1),WORK(KODBC2),
     &                        WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),
     &                        WORK(KRDPP1),WORK(KRDPP2),
     &                        WORK(KCCFB1),WORK(KINDXB),
     &                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CCRHSN')
               END IF
C
C               IF (T2TCOR) THEN
               KT2AMT = KEND1
               KEND1  = KT2AMT + NT2SQ(1)
               LWRK1  = LWORK  - KEND1
                  IF (LWRK1 .LT. 0) THEN
                     CALL QUIT('Insufficient core in CCRHSN')
                  END IF
C
               JSYM = 1
               CALL DCOPY(NT2SQ(1),T2AMSQ,1,WORK(KT2AMT),1)
               CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND1),LWRK1,JSYM)
C               CALL OUTPUT(T2AMSQ,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C               END IF
C
            ELSE
               NUMDIS = 1
               KRECNR = KENDSV
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
c
               WRITE(LUPRI,*) 'IDEL2',IDEL2
               FLUSH(LUPRI)
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
               IT2DEL(IDEL) = ICDEL1
               ICDEL1 = ICDEL1 + NT2BCD(ISYDIS)
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
                  CALL QUIT('Insufficient space in SO_2NDDECOEF')
               ENDIF
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
C                Reads distribution of AO integrals.
C
C------------------------------------------
C              Work space allocation no. 3.
C------------------------------------------
C
               KSCRM = KEND2
               KEND3 = KSCRM + NT2BCD(ISYMD)
               LWRK3 = LWORK - KEND3
C
               LOUT =NBAS(ISYMD)*NBAS(ISYMD)*NBAS(ISYMD)
C               KSCRM = KEND2
C               KXOUT = KSCRM + NT2BCD(ISYMD)
C               KEND3 = KXOUT + LOUT
C               LWRK3 = LWORK -KEND3
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient space in SO_2NDDECOEF')
               ENDIF
C
C----------------------------------------------------------------
C              Construct T_(ij)^(c delta)
C----------------------------------------------------------------
C
               DTIME   = SECOND()
C               ICON = 1
C               WRITE(LUPRI,*) 'T2AMSQ - før SO_T2AO'
C               CALL OUTPUT(T2AMSQ,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
               ISYMLH = 1
               CALL SO_T2AO(T2AMSQ,WORK(KCMO),ISYMLH,WORK(KSCRM),
     *                         WORK(KEND3),LWRK3,IDEL,ISYMD,
     *                         ISYMTR)
C              WORK(KSCRM) contains Sum_d  T_((ij)^(cd)*Lamda_(delta d)
               DTIME   = SECOND() - DTIME
               TIMT2AO = TIMT2AO + DTIME
C
C               WRITE(LUPRI,*) 'T2AMSQ - før SO_AGBD'
C               CALL OUTPUT(T2AMSQ,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C
C               WRITE(LUPRI,*) 'SCRM - før SO_AGBD'
C               CALL OUTPUT(WORK(KSCRM),1,NT2BCD(ISYMD),1,1,NT2BCD(ISYMD)
C     *                          ,1,1,LUPRI)
C
C               WRITE(LUPRI,*) 'OMEGA2'
C               CALL OUTPUT(OMEGA2,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C
C-----------------------------------
C              Calculate the sum_cd T_ij^cd (alpha gamme|beta delta).
C-----------------------------------
C
               DTIME   = SECOND()
               CALL SO_AGBD(WORK(KXINT),OMEGA2,WORK(KCMO),
     *                         WORK(KCMO),WORK(KSCRM),WORK(KEND3),
     *                         LWRK3,IDEL,ISYMD)
               DTIME   = SECOND() - DTIME
               TIMB    = TIMB     + DTIME
C
C               WRITE(LUPRI,*) 'OMEGA2'
C               CALL OUTPUT(OMEGA2,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C
C               WRITE(LUPRI,*) 'SCRM - efter SO_AGBD'
C               CALL OUTPUT(WORK(KSCRM),1,NT2BCD(ISYMD),1,1,NT2BCD(ISYMD)
C     *                          ,1,1,LUPRI)
C
C               WRITE(LUPRI,*) 'T2AMSQ - efter SO_AGBD'
C               CALL OUTPUT(T2AMSQ,1,NT2SQ(1),1,1,NT2SQ(1),1,1,LUPRI)
C               WRITE(LUPRI,*) 'NDISAO',NDISAO(ISYDIS)
C
C------------------------------------------
C              Work space allocation no. 4.
C------------------------------------------
C
               KDSRHF = KEND3
               KEND4  = KDSRHF + NDSRHF(ISYMD)
               LWRK4  = LWORK  - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND4,'Available : ',LWORK
                  CALL QUIT('Insufficient space in SO_2NDDEWF')
               ENDIF
C
C------------------------------------------------------------------------
C      Transform one index in the integral batch (alpha beta| j delta).
C------------------------------------------------------------------------
C
               DTIME   = SECOND()
               ISYMLP  = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     *                     ISYMLP,WORK(KEND4),LWRK4,ISYDIS)
               DTIME   = SECOND() - DTIME
               TIMTRBT = TIMTRBT + DTIME
C
C-------------------------------------------------------------
C              Calculate the gamma matrix entering the A-term.
C-------------------------------------------------------------
C
C               DTIME   = SECOND()
C               CALL SO_GAMMA(WORK(KDSRHF),WORK(KGAMMA),WORK(KCMO),
C     *                           WORK(KCMO),WORK(KSCRM),WORK(KEND4),
C     *                           LWRK4,IDEL,ISYMD)
C               DTIME   = SECOND() - DTIME
C               TIMGAM  = TIMGAM   + DTIME
C
C------------------------------------------------------------------------------
C              Calculate (ki|a delta).
C------------------------------------------------------------------------------
C
               DTIME   = SECOND()

               CALL SO_KIADELTA(WORK(KXINT),ISYMTR,
     *                         WORK(KCMO),WORK(KCMO),WORK(KCMO),
     *                         ISYMTR, WORK(KEND4), LWRK4,
     *                         IDEL,ISYMD,ZERO,LUC,CFIL)
C
               DTIME   = SECOND() - DTIME
               TIMC    = TIMC     + DTIME
C
C------------------------------------------------------------------------------
C              Calculate (ai|k delta) .
C------------------------------------------------------------------------------
C
               DTIME   = SECOND()
C
C               FACTD = HALF
C
C
               CALL SO_AIKDELTA(WORK(KDSRHF), ISYMTR, WORK(KCMO),
     *                          WORK(KCMO),ISYMTR,ISYMTR, WORK(KEND4),
     *                          LWRK4,IDEL,ISYMD,ZERO,LUD,DFIL)

               DTIME   = SECOND() - DTIME
               TIMD    = TIMD     + DTIME

C------------------------------------------------------------------------------
C              CHECK ROUTINE
               WRITE(LUPRI,'(A)') 'Runing CHECK_ROUTINES'
C------------------------------------------------------------------------------
C 
               CALL CHECK_AGBD(WORK(KXINT),
     *                          NDISAO(ISYMD),NBAS(ISYMD),NRHF,NVIR
     *                         ,WORK(KCMO),NLAMDT,T2AMSQ,NT2SQ(1),
     *                        WORK(KEND3),LWRK3,ILLL,ISYMD,LUPRI)
C
               CALL CHECK_KILJ(WORK(KXINT),
     *                          NDISAO(ISYMD),NBAS(ISYMD),NRHF,NVIR
     *                         ,WORK(KCMO),NLAMDT,T2AMSQ,NT2SQ(1),
     *                        WORK(KEND3),LWRK3,ILLL,ISYMD,LUPRI)
C
               CALL CHECK_KIADELTA(WORK(KXINT),
     *                          NDISAO(ISYMD),NBAS(ISYMD),NRHF,NVIR
     *                         ,WORK(KCMO),NLAMDT,T2AMSQ,NT2SQ(1),
     *                        WORK(KEND3),LWRK3,ILLL,ISYMD,LUPRI)

               CALL CHECK_AIKDELTA(WORK(KXINT),
     *                          NDISAO(ISYMD),NBAS(ISYMD),NRHF,NVIR
     *                         ,WORK(KCMO),NLAMDT,T2AMSQ,NT2SQ(1),
     *                        WORK(KEND3),LWRK3,ILLL,ISYMD,LUPRI)

  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C
C------------------------
C     Recover work space.
C------------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C-------------------------------------------------
C     Transform the Omega2 vector to the MO basis.
C-------------------------------------------------
C
      IF (NT2AM(ISYMOP) .GT. 2*NT2AMX) THEN
         WRITE(LUPRI,*)
     &        'Length of T2AM is smaller than OMEGA2 in MO basis'
         CALL QUIT('Insufficient space in SO_T2MO')
      ENDIF
C
      IF ( .NOT. CC2 ) THEN
C
C---------------------------------------
C        Save the CC amplitudes on disc - WHAT ARE WE SAVEING? CLilli
C---------------------------------------
C
C         WRITE (LURHS1) (T2AM(I), I = 1,NT2AM(ISYMOP))
C
C--------------------------------------------
C        Allocate space for the gamma matrix.
C--------------------------------------------
C
         IF (NEWGAM) THEN
C
            KGAMMA = KEND1
            KEND1  = KGAMMA + NGAMMA(ISYMOP)
            LWRK1  = LWORK  - KEND1
C
            IF (LWRK1 .LT. 0) CALL QUIT('Insufficient memory in GAMMA')
C
         END IF
C
C----------------------------------------------------
C        Transform the Omega2 vector to the MO basis.
C----------------------------------------------------
C
         IF (NT2AM(ISYMOP) .GT. 2*NT2AMX) THEN
            WRITE(LUPRI,*)
     *        'Length of T2AM is smaller than OMEGA2 in AO basis'
            CALL QUIT('Insufficient space in SO_T2MO')
         ENDIF
C
         TIMOME2 = SECOND()
         ISYMBF = ISYMOP
         ICON = 1
CLilli
         ISYMTR = 1
C         CALL SO_T2MO(OMEGA2,T2AM,
C     *                WORK(KCMO),WORK(KCMO),ISYMTR,
C     *                WORK(KEND1),LWRK1)
         CALL SO_T2MO(OMEGA2,OMEGA2,WORK(KCMO),WORK(KCMO),
     *                  ISYMTR,WORK(KEND1),LWRK1)
         CALL DCOPY(NT2AM(ISYMTR),T2AM,1,OMEGA2,1)
         TIMOME2 = SECOND() - TIMOME2
C
C-------------------------------
C        Print the Gamma matrix.
C-------------------------------
C
C         IF (IPRINT .GT. 120) THEN
            CALL AROUND('The Gamma matrix')
            DO 200 ISYM = 1,NSYM
               KOFF = KGAMMA + IGAMMA(ISYM,ISYM)
C               CALL OUTPAK(WORK(KOFF),NMATIJ(ISYM),1,LUPRI)
               CALL OUTPUT(WORK(KOFF),1,NMATIJ(ISYM),1,1,NMATIJ(ISYM),
     *                    1,1,LUPRI)
  200       CONTINUE
C
            WRITE(LUPRI,*) 'Norm of gamma matrix: ',
     *              DDOT(NGAMMA(ISYMOP),WORK(KGAMMA),1,WORK(KGAMMA),1)
C         END IF
C
      ENDIF
C
      WRITE(LUPRI,*) 'Norm of gamma matrix: ',
     *              DDOT(NGAMMA(ISYMOP),WORK(KGAMMA),1,WORK(KGAMMA),1)
      FLUSH(LUPRI)
C
C-------------------------------------------------------
C     Calculate Sum_(kl) t_(kl)^(cd)(ki|lj).
C-------------------------------------------------------
C
      TIMA     = SECOND()
      CALL SO_KILJ(OMEGA2,T2AM,WORK(KGAMMA),WORK(KEND2),LWRK2,
     *               ISYMTR,ISYMTR)
C
      TIMA     = SECOND() - TIMA
C
C----------------------------------------------------------
C      Calculate Sum_(kc) t_(jk)^(cb)*(ki|ac)
C----------------------------------------------------------
C
      ISYVEC = 1
      ISYCIM = 1
      IOPT   = 1
      IVECNR = 1
C
      TIMCIO = SECOND()
C      IF (T2TCOR) THEN
      CALL SO_KIADELTA_IO(OMEGA2,T2AMSQ,WORK(KCMO),
     *                     WORK(KEND2),LWRK2,ISYVEC,ISYCIM,
     *                     LUC,CFIL,IVECNR,IOPT)
C      ELSE
C            ISYM = 1
C            CALL CCSD_T2TP(T2AM,WORK(KEND2),LWRK2,ISYM)
C            CALL CCRHS_CIO(OMEGA2,T2AM,WORK(KLAMDH),WORK(KEND2),
c     *                     LWRK2,ISYVEC,ISYCIM,LUC,CFIL,IVECNR,IOPT)
C            CALL CCSD_T2TP(T2AM,WORK(KEND2),LWRK2,ISYM)
C         ENDIF
C
      TIMCIO  = SECOND() - TIMCIO
C
C------------------------------
C     Transform T2 to 2T2 - T2.
C------------------------------
C
      DTIME    = SECOND()
C      IF (T2TCOR) THEN
      CALL DSCAL(NT2SQ(1),TWO,T2AM,1)
      CALL DAXPY(NT2SQ(1),-ONE,WORK(KT2AMT),1,T2AM,1)
C      ELSE
C         ISYM = 1
C         CALL CCRHS_T2TR(T2AM,WORK(KEND2),LWRK2,ISYM)
C      END IF
      DTIME    = SECOND() - DTIME
      TIMT2TR  = TIMT2TR + DTIME
C
C--------------------------------------------------
C     Calculate Sum_(kc) T_(jk)^(cb)*(ai|kc)
C------------------------------------------------------
C
      ISYDIM = 1
      ISYVEC = 1
      IOPT = 1
      IVECNR = 1
C
      TIMDIO = SECOND()
      CALL SO_AIKDELTA_IO(OMEGA2,T2AM,WORK(KCMO),WORK(KEND2),LWRK2,
     *                  ISYVEC,ISYDIM,LUD,DFIL,IVECNR,IOPT)
      TIMDIO  = SECOND() - TIMDIO
C      END IF
      WRITE(LUPRI,'(A)') 'SO_AIKDELTA_IO - DONE '
C      FLUSH(LUPRI)
C
C------------------------
C     Scale final result.
C------------------------
CC
C     CALL DSCAL(NT2IND,TWO,OMEGA2,1)
CC
C      IF (IPRINT .GT. 25) THEN
C         WRITE(LUPRI,*)
C         CALL AROUND('END OF 2ndDEcoef:OMEGA 2')
C         DO 310 ISYM = 1,NSYM
C            WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
C            KOFF = IT2AM(ISYM,ISYM) + 1
CCLilli            CALL OUTPAK(OMEGA2(KOFF),NT1AM(ISYM),1,LUPRI)
C            CALL OUTPAK(OMEGA2(KOFF),NT2AM(ISYM),1,LUPRI)
C  310    CONTINUE
CC      ENDIF
      TIMALL  = SECOND() - TIMALL
C      IF ( IPRINT .GT. 2) THEN
C         WRITE(LUPRI,9999) 'RHS - TOTAL', TIMALL
C      ENDIF
C      IF (IPRINT .GT. 9) THEN
C         WRITE(LUPRI,9999) 'SO_KILJ    ', TIMA
C         WRITE(LUPRI,9999) 'SO_ACBD    ', TIMB
C         WRITE(LUPRI,9999) 'CCRHS_C    ', TIMC
C         WRITE(LUPRI,9999) 'CCRHS_CIO  ', TIMCIO
C         WRITE(LUPRI,9999) 'CCRHS_C-tot', TIMCIO + TIMC
C         WRITE(LUPRI,9999) 'CCRHS_D    ', TIMD
C         WRITE(LUPRI,9999) 'CCRHS_DIO  ', TIMDIO
C         WRITE(LUPRI,9999) 'CCRHS_D-tot', TIMDIO + TIMD
C         WRITE(LUPRI,9999) 'CCRHS_GAM  ', TIMGAM
C         WRITE(LUPRI,9999) 'SO_GETMO   ', TIMCMO
C         WRITE(LUPRI,9999) 'CCRHS_RDAO ', TIMRDAO
C         WRITE(LUPRI,9999) 'HERDIS1    ', TIMHER1
C         WRITE(LUPRI,9999) 'HERDIS2    ', TIMHER2
C         WRITE(LUPRI,9999) 'CC_T2AO    ', TIMT2AO
C         WRITE(LUPRI,9999) 'CC_FCKMO   ', TIMFCKMO
C         WRITE(LUPRI,9999) 'CCRHS_FCK  ', TIMFCK
C         WRITE(LUPRI,9999) 'CC_AODENS  ', TIMDM
C         WRITE(LUPRI,9999) 'CCRHS_TRBT ', TIMTRBT
C         WRITE(LUPRI,9999) 'CCRHS_T2TR ', TIMT2TR
C         WRITE(LUPRI,9999) 'CCRHS_T2BT ', TIMT2BT
C      ENDIF
C9999  FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')

      LUC = -1
      LUD = -1
      CFIL = 'PMAT_C'
      DFIL = 'PMAT_D'
      CALL WCLOSE2(LUC,CFIL,'KEEP')
      CALL WCLOSE2(LUD,DFIL,'KEEP')

      CALL QEXIT('SO_2NDDECOEF')

      RETURN
      END
