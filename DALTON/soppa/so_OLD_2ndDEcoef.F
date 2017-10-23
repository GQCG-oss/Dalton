C  /* Deck so_2ndDEcoef */
      SUBROUTINE SO_2NDDECOEF(OMEGA2,T2AM,T2AMSQ,WORK,LWORK)
C      SUBROUTINE SO_2NDDECOEF(OMEGA1,OMEGA2,T1AM,T2AM,WORK,LWORK)
C
C     Orignial routine, CCRHSN, is written by Henrik Koch 25-Sep-1993
C
C     Purpose: calculate the 2nd order doubles coeffecients.
C
C     Calculation of the Coupled Cluster vector function using
C     AO-integrals directly from disk. 
C
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
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (XMHALF = -0.5D0, XMONE= -1.0D0 )
      PARAMETER (ISYM0 = 1)
C
C      LOGICAL FCKCON,CC1BSA,ETRAN,LV,LVAJKL,LRES,
C     &        DEBUGV,LVIJKL,LVABKL
C      PARAMETER (DEBUGV = .FALSE.)
C
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION OMEGA2(*),T2AM(*),WORK(LWORK),T2AMSQ(NT2SQ(1))
C
      CHARACTER CFIL*6,DFIL*6
C      CHARACTER CFIL*6,DFIL*6, FN3SRT*8, FNDELD*6, CDUMMY*8
C      CHARACTER FNCKJD*6, FNDKBC*4, FNTOC*8, FN3VI*6, FN3VI2*8
C      CHARACTER FNIADJ*8, FNIJDA*8, CPFIL*8, DPFIL*8
      CHARACTER MODEL*10
C
C      PARAMETER (FNIADJ = 'CCXIADJ0', FNIJDA = 'CCXIJDA0')
C      PARAMETER (CPFIL  = 'CC_CPR12', DPFIL  = 'CC_DPR12')
C

      INTEGER IGLMRHS(8,8),IGLMVIS(8,8),NGLMDS(8),ICMO(8,8),NCMO(8),
     &        IMAIJM(8,8),NMAIJM(8),
     &        IMATIJM(8,8),NMATIJM(8),NGAMSM(8),IGAMSM(8,8),
     &        IRGIJS(8,8),NRGIJS(8),IR1BASM(8,8),NR1BASM(8),
     &        IR2BASM(8,8),NR2BASM,IR1XBASM(8,8),NR1XBASM(8),
     &        IR2XBASM(8,8),IMATF(8,8),NMATF(8),IMAKLM(8,8),NMAKLM(8)
      INTEGER NADP(8),IADP(8,8),NLAMDX(8),ILAMDX(8,8)
C
      LOGICAL MLCC3_RESPONSE
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
C      IF (DEBUG) WRITE(LUPRI,*) 'DUMPCD = ',DUMPCD
      CALL WOPEN2(LUC,CFIL,64,0)
      CALL WOPEN2(LUD,DFIL,64,0)
C
C      END IF
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
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      KLAMDP = 1
      KLAMIP = KLAMDP + NLAMDT
      KLAMDH = KLAMIP + 1
C      KLAMDH = KLAMIP + NLAMDT
      KDENSI = KLAMDH + NLAMDT
      KFOCK  = KDENSI + N2BAST
C      KEMAT1 = KFOCK  + N2BST(ISYMOP)
C      KEMAT2 = KEMAT1 + NEMAT1(ISYMOP)
C      KGAMMA = KEMAT2 + NMATIJ(ISYMOP)
      KGAMMA = KFOCK + NMATIJ(ISYMOP)
      KEND1 = KGAMMA
C
C      IF (CC2 .AND. NONHF) THEN
C        KFCKHF = KEND1
C        KEND1  = KFCKHF + N2BAST
C      END IF
c
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in NAME')
      ENDIF
C
      WRITE(LUPRI, '(A)') 'ALLOCATION DONE'
C------------------------------------
C     Save the CC amplitudes on disk
C------------------------------------
CDon't need this when we don't overwrite T2AM with T2AMSQ
C      LURHS1 = -1
C      CALL GPOPEN(LURHS1,'CCRHS1','UNKNOWN',' ','UNFORMATTED',IDUMMY,
C     &            .FALSE.)
C      REWIND(LURHS1)
C      WRITE (LURHS1) (T1AM(I), I = 1,NT1AMX)
C      WRITE (LURHS1) (T2AM(I), I = 1,NT2AMX)
C
C      WRITE(LUPRI, '(A)') 'CC amplitudes saved'
C
C----------------------------------
C     Calculate the lamda matrices.
C----------------------------------
C
C
      TIMCMO  = SECOND()
      CALL SO_GETMO(WORK(KCMO),NLAMDT,WORK(KEND1),LWRK1)
C      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND1),LWRK1)
      TIMCMO  = SECOND() - TIMCMO
C
      WRITE(LUPRI, '(A)') 'I got the Mo-coeffecients'
      FLUSH(LUPRI)
C-------------------------------
C     Prepare the t2-amplitudes.
C-------------------------------
C
C      CALL DCOPY(NT2AMX,T2AM,1,OMEGA2,1)
C           Copies T2AM to OMEGA2
C      CALL CC_T2SQ(OMEGA2,T2AM,ISYMTR)
C           Squareup the t2-amplitudes distribution.
      CALL CC_T2SQ(T2AM,T2AMSQ,ISYMTR)
C           T2AMSQ is en array fittet to contain the squared up T2AM.
C
      WRITE(LUPRI, '(A)') 'T2-amplitudes prepared'
      FLUSH(LUPRI)
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
         WRITE(LUPRI, '(A)') ' Direct: transposed amplitudes(1)'
         FLUSH(LUPRI)
C
         JSYM = 1
         CALL DCOPY(NT2SQ(1),T2AMSQ,1,WORK(KT2AMT),1)
CLilli         CALL DCOPY(NT2SQ(1),T2AM,1,WORK(KT2AMT),1)
C               Copies T2AM to WORK(KT2AMT)
         CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND1),LWRK1,JSYM)
C         CCSD_T2TP transposes the ij index of the squaed T2-amplitudes.
C           They are put in WORK(KT2AMT)
      END IF
      WRITE(LUPRI, '(A)') 'Transposed amplitudes'
      FLUSH(LUPRI)
C
C-------------------------------
C     Initialize OMEGA1 & OMEGA2
C-------------------------------
C
C      CALL DZERO(OMEGA1,NT1AM(ISYMOP))
C      IF (.NOT. OMEGSQ) THEN
C      IF (OMEGOR) THEN
      CALL DZERO(OMEGA2,2*NT2ORT(ISYMOP))
C           NT2ORT(ISYM)=NNBST(ISYMI)*NMIJP(ISYMJ)
C           NNBST(ISYM)=NBAS(ISYMA)*NBAS(ISYMB)
C           NMIJP(ISYM)=NRHF(ISYMI)*NRHF(ISYMJ)
C   OMEGA2 and OMEGA1 are now 0.
C-------------------------------------
C     Initialize GAMMA, EMAT1 & EMAT2.
C-------------------------------------
C
C      CALL DZERO(WORK(KEMAT1),NEMAT1(ISYMOP))
C      CALL DZERO(WORK(KEMAT2),NMATIJ(ISYMOP))
C
C----------------------------------------
C     Calculate the density matrix.
C     includes core contribution, ic = 1.
C----------------------------------------
C
      TIMDM  = SECOND()
      ISYMH = 1
      IC    = 1
      CALL CC_AODENS(WORK(KLAMDP),WORK(KLAMDH),WORK(KDENSI),ISYMH,
     *               IC,WORK(KEND1),LWRK1)
      TIMDM  = SECOND() - TIMDM
C
      WRITE(LUPRI, '(A)') 'The debsity matrix has been calculated'
      FLUSH(LUPRI)
C------------------------------------------------
C     Read one-electron integrals in Fock-matrix - this i might need to
C                                                   chance?!?
C------------------------------------------------
C
      TIMFCK = SECOND()
      CALL CCRHS_ONEAO(WORK(KFOCK),WORK(KEND1),LWRK1)
      TIMFCK = SECOND() - TIMFCK
C
C------------------------------------------------
C     Read one-electron integrals into Fock-matrix for
C     finite field.
C------------------------------------------------
C
      DO 13 IF = 1, NFIELD
         DTIME  = SECOND()
         FF = EFIELD(IF)
         CALL CC_ONEP(WORK(KFOCK),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
         DTIME  = SECOND() - DTIME
         TIMFCK = TIMFCK + DTIME
 13   CONTINUE
C
      WRITE(LUPRI, '(A)') '1-elec. intrals done. Now for the loop'
      FLUSH(LUPRI)
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
               CALL DCOPY(NT2SQ(1),T2AM,1,WORK(KT2AMT),1)
               CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND1),LWRK1,JSYM)
C               END IF
C
            ELSE
               NUMDIS = 1
               KRECNR = KENDSV
            ENDIF
C
            WRITE(LUPRI, '(A)') 'After direct calculation of t2'
            FLUSH(LUPRI)
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
C                  IF (NOAUXB) THEN
C                     IDUM = 1
C                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
C                  END IF
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
               WRITE(LUPRI, '(A)') 'Allocation no. 2'
               FLUSH(LUPRI)
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHSN')
               ENDIF
C
               WRITE(LUPRI, '(A)') 'Allocation no. 2 done'
               FLUSH(LUPRI)
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
               WRITE(LUPRI, '(A)') 'Batch of integral has been read'
               FLUSH(LUPRI)
C-------------------------------------------
C              Calculate the AO-Fock matrix.
C-------------------------------------------
C
               DTIME   = SECOND()
C
               ISYDEN = 1
               CALL CC_AOFOCK(WORK(KXINT),WORK(KDENSI),WORK(KFOCK),
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,.FALSE.,
     *                        DUMMY,ISYDEN)
               DTIME  = SECOND() - DTIME
               TIMFCK = TIMFCK + DTIME
C
C------------------------------------------
C              Work space allocation no. 3.
C------------------------------------------
C
               KSCRM = KEND2
               KEND3 = KSCRM + NT2BCD(ISYMD)
               LWRK3 = LWORK - KEND3
C
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
               ICON = 1
               ISYMLH = 1
               CALL SO_T2AO(T2AM,WORK(KCMO),ISYMLH,WORK(KSCRM),
     *                         WORK(KEND3),LWRK3,IDEL,ISYMD,
     *                         ISYMTR,ICON)
C              WORK(KSCRM) contains Sum_d  T_((ij)^(cd)*Lamda_(delta d)
               DTIME   = SECOND() - DTIME
               TIMT2AO = TIMT2AO + DTIME
C
               WRITE(LUPRI, '(A)') 'SCRM has calculated'
               FLUSH(LUPRI)
C-----------------------------------
C              Calculate the (ac|bd).
C-----------------------------------
C
               DTIME   = SECOND()
               CALL SO_ACBD(WORK(KXINT),OMEGA2,WORK(KLAMDP),
     *                         WORK(KLAMDH),WORK(KSCRM),WORK(KEND3),
     *                         LWRK3,IDEL,ISYMD)
C               ENDIF
               DTIME   = SECOND() - DTIME
               TIMB    = TIMB     + DTIME
C
               WRITE(LUPRI, '(A)') 'SO_ACBD has been calculated'
               FLUSH(LUPRI)
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
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                     ISYMLP,WORK(KEND4),LWRK4,ISYDIS)
               DTIME   = SECOND() - DTIME
               TIMTRBT = TIMTRBT + DTIME
C
               WRITE(LUPRI, '(A)') 'CCTRBT has run'
               FLUSH(LUPRI)
C-------------------------------------------------------------
C              Calculate the gamma matrix entering the A-term.
C-------------------------------------------------------------
C
               DTIME   = SECOND()
               CALL SO_GAMMA(WORK(KDSRHF),WORK(KGAMMA),WORK(KLAMDP),
     *                           WORK(KLAMDH),WORK(KSCRM),WORK(KEND4),
     *                           LWRK4,IDEL,ISYMD)
               DTIME   = SECOND() - DTIME
               TIMGAM  = TIMGAM   + DTIME
C
               WRITE(LUPRI, '(A)') 'SO_GAMMA - DONE '
               FLUSH(LUPRI)
C------------------------------------------------------------------------------
C              Calculate (ki|a delta).
C------------------------------------------------------------------------------
C
               DTIME   = SECOND()
               FACTC = XMHALF
C
               CALL SO_KIADELTA(WORK(KXINT),ISYMTR,
     *                         WORK(KLAMDP),WORK(KLAMDH),WORK(KLAMDP),
     *                         ISYMTR, WORK(KEND4), LWRK4,
     *                         IDEL,ISYMD,FACTC,LUC,CFIL)
C
               DTIME   = SECOND() - DTIME
               TIMC    = TIMC     + DTIME
C
               WRITE(LUPRI, '(A)') 'SO_KIADELTA - DONE '
               FLUSH(LUPRI)
C------------------------------------------------------------------------------
C              Transform T2 to 2T2 - T2.
C------------------------------------------------------------------------------
C
C               DTIME   = SECOND()
c              IF (T2TCOR) THEN
C               CALL DSCAL(NT2SQ(1),TWO,T2AM,1)
C               CALL DAXPY(NT2SQ(1),-ONE,WORK(KT2AMT),1,T2AM,1)
C               ELSE
C                  ISYM = 1
C                  CALL CCRHS_T2TR(T2AM,WORK(KEND4),LWRK4,ISYM)
C               END IF
C               DTIME   = SECOND() - DTIME
C               TIMT2TR = TIMT2TR  + DTIME
C
C------------------------------------------------------------------------------
CLilli DONT think we need this              Transform the cluster amplitudes.
C------------------------------------------------------------------------------
C
C               CALL CC_MTCME(WORK(KSCRM),WORK(KEND4),LWRK4,
C     *                       ISYMD,ISYMTR)
C
C------------------------------------------------------------------------------
C              Calculate (ai|k delta) .
C------------------------------------------------------------------------------
C
               DTIME   = SECOND()
C
               FACTD = HALF
C
C
               CALL SO_AIKDELTA(WORK(KDSRHF), ISYMTR, WORK(KLAMDH),
     *                          WORK(KLAMDP),ISYMTR,ISYMTR, WORK(KEND4),
     *                          LWRK4,IDEL,ISYMD,FACTD,LUD,DFIL)
C
               DTIME   = SECOND() - DTIME
               TIMD    = TIMD     + DTIME
C
               WRITE(LUPRI, '(A)') 'SO_AIKDELTA - DONE '
               FLUSH(LUPRI)
C---------------------------------------------
C              BackTransform T2 from 2T2 - T2.
C---------------------------------------------
C
C               DTIME   = SECOND()
C               IF (T2TCOR) THEN
C               CALL DAXPY(NT2SQ(1),ONE,WORK(KT2AMT),1,T2AM,1)
C               CALL DSCAL(NT2SQ(1),HALF,T2AM,1)
C               ELSE
C                  ISYM = 1
C                  CALL CCRHS_T2BT(T2AM,WORK(KEND4),LWRK4,ISYM)
C               END IF
C               DTIME   = SECOND() - DTIME
C               TIMT2BT = TIMT2BT  + DTIME
C
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
C      IF (IPRINT .GT. 120) THEN
C         CALL AROUND('After  Delta Loop: Omega1')
C         CALL CC_PRP(OMEGA1,OMEGA2,1,1,0)
C      ENDIF
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
C
         CALL SO_T2MO(OMEGA2,T2AM,WORK(KGAMMA),
     *                WORK(KLAMDP),WORK(KLAMDP),ISYMTR,
     *                WORK(KEND1),LWRK1)
C     Transform the Omega2 vector from the AO basis to the MO
C     basis and store it in WORK(KGAMMA)
         WRITE(LUPRI, '(A)') 'SO_T2MO - DONE'
         FLUSH(LUPRI)
         CALL DCOPY(NT2AM(ISYMTR),T2AM,1,OMEGA2,1)
         TIMOME2 = SECOND() - TIMOME2
C
         IF (IPRINT .GT. 51) THEN
C            RHO1N = DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
            RHO2N = DDOT(NT2AMX,OMEGA2,1,OMEGA2,1)
C            WRITE(LUPRI,*) 'Norm of OMEGA1 -after CC_T2MO: ',RHO1N
            WRITE(LUPRI,*) 'Norm of OMEGA2 -after CC_T2MO: ',RHO2N
         ENDIF
C
C         IF (IPRINT .GT. 120) THEN
C            CALL AROUND('After  T2MO: BF ')
C            CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
C         ENDIF
C
C-------------------------------
C        Print the Gamma matrix.
C-------------------------------
         WRITE(LUPRI,'(A)') 'Priting the Gamma Matrix'
         FLUSH(LUPRI)
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
C--------------------------------------------
C        Restore the CC amplitudes from disk.
C--------------------------------------------
C         WRITE(LUPRI,'(A)') 'Restore cc amplitudes'
C         FLUSH(LUPRI)
C
C         REWIND (LURHS1)
C         READ (LURHS1)
C         READ (LURHS1)
C         READ (LURHS1) (T2AM(I), I = 1,NT2AM(ISYMOP))
C
      ENDIF
C
      WRITE(LUPRI,*) 'Norm of gamma matrix: ',
     *              DDOT(NGAMMA(ISYMOP),WORK(KGAMMA),1,WORK(KGAMMA),1)
      FLUSH(LUPRI)
C------------------------------------------
C     Transform AO Fock matrix to MO basis.
C------------------------------------------
      WRITE(LUPRI,'(A)') 'Transform Fock Matrix to MO basis'
      FLUSH(LUPRI)
C
      TIMFCKMO = SECOND()
      CALL CC_FCKMO(WORK(KFOCK),WORK(KCMO),WORK(KCMO),
     *                 WORK(KEND1),LWRK1,1,1,1)
      TIMFCKMO = SECOND() - TIMFCKMO
C
      WRITE(LUPRI,'(A)') 'CC_FCKMO - done'
      FLUSH(LUPRI)
C---------------------
C     Reallocate T2TP. We already have an array for T2AMSQ.
C---------------------
C
C      IF (DIRECT .AND. T2TCOR) THEN
C      IF (DIRECT) THEN
C
C         KT2AMT = KEND1
C         KEND2  = KT2AMT + NT2SQ(1)
C         LWRK2  = LWORK  - KEND2
C
C         IF (LWRK2. LT. 0) THEN
C            CALL QUIT('Insufficient memory in CCSD_RHS')
C         END IF
C
C      ELSE
C
C         KEND2 = KEND1
C         LWRK2 = LWRK1
C
C      END IF

C----------------------------------------
C        Reconstruct full square of T2AM.  Do I need this with T2AMSQ?
C----------------------------------------
C
C      IF (LWRK1 .LT. NT2AMX) THEN
C      CALL QUIT('Insufficient core in SO_2NDDECOEF')
C      ENDIF
C
C      REWIND (LURHS1)
c      READ (LURHS1)
C      READ (LURHS1) (WORK(KEND1+I-1), I = 1,NT2AMX)
C
C      CALL CC_T2SQ(WORK(KEND1),T2AM,1)
C
C      WRITE(LUPRI,'(A)') 'CC_T2SQ - done'
C      FLUSH(LUPRI)
C
C      IF (IPRINT .GT. 51) THEN
C         RHO1N = DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
C         RHO2N = DDOT(NT2AMX,OMEGA2,1,OMEGA2,1)
C         WRITE(LUPRI,*) 'Norm of OMEGA1 -after cc3_omeg: ',RHO1N
C         WRITE(LUPRI,*) 'Norm of OMEGA2 -after cc3_omeg: ',RHO2N
C      ENDIF
C
C      IF (IPRINT .GT. 120) THEN
C         CALL AROUND('After  CC3_OMEG Omega is ')
C         CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
C      ENDIF
C
C---------------------
C     Reallocate T2TP.
C---------------------
C      WRITE(LUPRI,'(A)') 'Reallocate T2TP'
C      FLUSH(LUPRI)
C
C      IF ((DIRECT .AND. T2TCOR) .OR. (CCSDT .AND. T2TCOR)) THEN
C      IF (DIRECT) THEN
C
C         KT2AMT = KEND1
C         KEND2  = KT2AMT + NT2SQ(1)
C         LWRK2  = LWORK  - KEND2
C
C         IF (LWRK2. LT. 0) THEN
C            CALL QUIT('Insufficient memory in CCSD_RHS')
C         END IF
C
C      ELSE
C
C         KEND2 = KEND1
C         LWRK2 = LWRK1
C
C      END IF
C
C----------------------
C     Recalculate T2TP.
C----------------------
C
C      IF (T2TCOR) THEN
C
C      JSYM = 1
CLilli The full square of T2AM is already in T2AMSQ
C      CALL DCOPY(NT2SQ(1),T2AM,1,WORK(KT2AMT),1)
c      CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND2),LWRK2,JSYM)
C
C      END IF
C      WRITE(LUPRI,'(A)') 'T2TP recalculated'
C      FLUSH(LUPRI)
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
      WRITE(LUPRI,'(A)') 'SO_KILJ-DONE '
      FLUSH(LUPRI)
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
C      CALL SO_KIADELTA_IO(OMEGA2,WORK(KT2AMT),WORK(KLAMDH),
     *                     WORK(KEND2),LWRK2,ISYVEC,ISYCIM,
     *                     LUC,CFIL,IVECNR,IOPT)
C      ELSE
      WRITE(LUPRI,'(A)') 'SO_KIADELTA_IO - DONE '
      FLUSH(LUPRI)
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
      CALL SO_AIKDELTA_IO(OMEGA2,T2AM,WORK(KLAMDH),WORK(KEND2),LWRK2,
     *                  ISYVEC,ISYDIM,LUD,DFIL,IVECNR,IOPT)
      TIMDIO  = SECOND() - TIMDIO
C      END IF
      WRITE(LUPRI,'(A)') 'SO_AIKDELTA_IO - DONE '
      FLUSH(LUPRI)
C
C------------------------
C     Scale final result.
C------------------------
C
C     CALL DSCAL(NT1AM,TWO,OMEGA1,1)
C     CALL DSCAL(NT2IND,TWO,OMEGA2,1)
C
C      IF (IPRINT .GT. 25) THEN
         WRITE(LUPRI,*)
         CALL AROUND('END OF 2ndDEcoef:OMEGA 2')
         DO 310 ISYM = 1,NSYM
            WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
            KOFF = IT2AM(ISYM,ISYM) + 1
CLilli            CALL OUTPAK(OMEGA2(KOFF),NT1AM(ISYM),1,LUPRI)
            CALL OUTPAK(OMEGA2(KOFF),NT2AM(ISYM),1,LUPRI)
  310    CONTINUE
C      ENDIF
      TIMALL  = SECOND() - TIMALL
      IF ( IPRINT .GT. 2) THEN
         WRITE(LUPRI,9999) 'RHS - TOTAL', TIMALL
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,9999) 'SO_KILJ    ', TIMA
         WRITE(LUPRI,9999) 'SO_ACBD    ', TIMB
         WRITE(LUPRI,9999) 'CCRHS_C    ', TIMC
         WRITE(LUPRI,9999) 'CCRHS_CIO  ', TIMCIO
         WRITE(LUPRI,9999) 'CCRHS_C-tot', TIMCIO + TIMC
         WRITE(LUPRI,9999) 'CCRHS_D    ', TIMD
         WRITE(LUPRI,9999) 'CCRHS_DIO  ', TIMDIO
         WRITE(LUPRI,9999) 'CCRHS_D-tot', TIMDIO + TIMD
         WRITE(LUPRI,9999) 'CCRHS_GAM  ', TIMGAM
         WRITE(LUPRI,9999) 'CCRHS_LAM  ', TIMCMO
         WRITE(LUPRI,9999) 'CCRHS_RDAO ', TIMRDAO
         WRITE(LUPRI,9999) 'HERDIS1    ', TIMHER1
         WRITE(LUPRI,9999) 'HERDIS2    ', TIMHER2
         WRITE(LUPRI,9999) 'CC_T2AO    ', TIMT2AO
         WRITE(LUPRI,9999) 'CC_FCKMO   ', TIMFCKMO
         WRITE(LUPRI,9999) 'CCRHS_FCK  ', TIMFCK
         WRITE(LUPRI,9999) 'CCRHS_DM   ', TIMDM
         WRITE(LUPRI,9999) 'CCRHS_TRBT ', TIMTRBT
         WRITE(LUPRI,9999) 'CCRHS_T2TR ', TIMT2TR
         WRITE(LUPRI,9999) 'CCRHS_T2BT ', TIMT2BT
      ENDIF
9999  FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')
C
C-----------------------------------------
C     Restore the CC amplitudes from disk.
C-----------------------------------------
C
C      REWIND (LURHS1)
C      READ(LURHS1) (T1AM(I), I = 1,NT1AMX)
C      READ(LURHS1) (T2AM(I), I = 1,NT2AMX)
C      CALL GPCLOSE(LURHS1,'DELETE')
C
C-----------------
C     Close files.
C-----------------
C
C      IF (DUMPCD) THEN is always true
      CALL WCLOSE2(LUC,CFIL,'KEEP')
      CALL WCLOSE2(LUD,DFIL,'KEEP')
C      END IF
C
C      IF (CCSDT) THEN
C         CALL WCLOSE2(LU3SRT,FN3SRT,'KEEP')
C         CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
C         CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
C         CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
C         CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
C         CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
C         CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
C      ENDIF
C
      CALL QEXIT('SO_2NDDECOEF')
C
      RETURN
      END
