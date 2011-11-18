C  /* Deck so_res_cp */
      SUBROUTINE SO_RES_CP(NNEWTR, DENSIJ, LDENSIJ, DENSAB, LDENSAB,
     &                     T2MP,   LT2MP,  FOCKD,   LFOCKD,
     &                     ISYMTR, NEXCI,  WORK,    LWORK,
     &                     TR1E,   TR1D,   TR2E,    TR2D,
     &                     RES1E,  RES1D,  RES2E,   RES2D,
     &                     RESO1E, RESO1D,
     &                     LTR1E,  LTR1D,  LTR2E,   LTR2D,
     &                     LRES1E, LRES1D, LRES2E,  LRES2D,
     &                     LRESO1E,LRESO1D)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Driver routine for making a linear transformation of
C              a trialvector with the SOPPA hessian and overlap
C              matrices E[2] and S[2]. The trial vector konsists of
C              four parts TR1E, TR1D, TR2E, and TR2D. E refers to
C              excitations and D to de-excitations. 1 refer to the one-
C              particle part and 2 to the two-particle part. The linear
C              transformed trialvector is refered to as the resultvector
C              and is kept in four corresponding arrays. For the linear
C              transformation with E[2] the result vector is in RES1E,
C              RES1D, RES2E, and RES2D and for the linear transformation
C              with S[2] in RESO1E and RESO1D. The linear transformation
C              with the two-particle part of S[2] which equals the
C              identity matrix is not carried out.
C              The linear transformation is driven over atomic orbitals,
C              and E[2] and S[2] are not constructed explicitly.
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
CSPAS:10/11-2003: is merged with maxorb.h
C#include <mxorb.h>
CKeinSPASmehr
#include <aovec.h>
#include <iratdef.h>
#include <eribuf.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), T2MP(LT2MP)
      DIMENSION TR1E(LTR1E),   TR1D(LTR1D),  TR2E(LTR2E),  TR2D(LTR2D)
      DIMENSION RES1E(LTR1E),  RES1D(LTR1D), RES2E(LTR2E), RES2D(LTR2D)
      DIMENSION RESO1E(LTR1E), RESO1D(LTR1D)
      DIMENSION FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
      dimension rint(10,10,10,10)
      logical   sonod
C
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_RES_CP')
C
cKeld      sonod  = .true.
      sonod  = .false.
      if (sonod) then
         call dzero(rint,10000)
      end if
C
C------------------------------------------------------------------
C     Determine the symmetri of the result vector from the symmetry
C     of the trial vector ISYMTR, and the opperator symmtry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      LCMO  = NLAMDT
C
      KCMO    = 1
      KEND1   = KCMO   + LCMO
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_RES_CP.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_RES_CP.1',' ',KEND1,LWORK)
C
C--------------------------------------------------------------------
C     Calculate the lamda matrices which contain the MO coefficients.
C--------------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL LAMMAT(WORK(KCMO),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     &            LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(1)  = SOTIME(1) + DTIME
C
C---------------------------------
C     Work space allocation no. 2.
C---------------------------------
C
      LFOCK   = N2BST(ISYRES)
      LDENS   = N2BST(ISYMTR)
      LBTR1E  = NT1AO(ISYMTR)
      LBTR1D  = NT1AO(ISYMTR)
      LBTJ1E  = NMATAV(ISYMTR)
      LBTJ1D  = NMATAV(ISYMTR)
      LSIGAI1 = NT1AO(ISYRES)
      LSIGAI2 = NT1AO(ISYRES)
      LSIGDA1 = NMATAV(ISYRES)
      LSIGDA2 = NMATAV(ISYRES)
      LAIJ    = NRHFT*NRHFT
      LAAB    = NVIRT*NVIRT
C
      KFOCK   = KEND1
      KDENS   = KFOCK   + LFOCK
      KBTR1E  = KDENS   + LDENS
      KBTR1D  = KBTR1E  + LBTR1E
      KBTJ1E  = KBTR1D  + LBTR1D
      KBTJ1D  = KBTJ1E  + LBTJ1E
      KSIGAI1 = KBTJ1D  + LBTJ1D
      KSIGAI2 = KSIGAI1 + LSIGAI1
      KSIGDA1 = KSIGAI2 + LSIGAI2
      KSIGDA2 = KSIGDA1 + LSIGDA1
      KAIJ    = KSIGDA2 + LSIGDA2
      KAAB    = KAIJ    + LAIJ
      KEND2   = KAAB    + LAAB
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('SO_RES_CP.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_RES_CP.2',' ',KEND2,LWORK)
C
c      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
c      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
c      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
c      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
      CALL SO_OPEN(LUDENS,FNDENS,LDENS)
      CALL SO_OPEN(LUBT1E,FNBT1E,LBTR1E)
      CALL SO_OPEN(LUBT1D,FNBT1D,LBTR1D)
      CALL SO_OPEN(LUBJ1E,FNBJ1E,LBTJ1E)
      CALL SO_OPEN(LUBJ1D,FNBJ1D,LBTJ1D)
C
C---------------------------------------
C     Write new trial vectors to output.
C---------------------------------------
C
c      DO 50 INEWTR = 1,NNEWTR
C
c         CALL AROUND('New trial vector in SO_RES_CP')
C
c         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,NEXCI+INEWTR)
C
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTR1E+I-1),I=1,LTR1E)
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTR2E+I-1),I=1,LTR2E)
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTR1D+I-1),I=1,LTR1D)
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KTR2D+I-1),I=1,LTR2D)
C
C   50 CONTINUE
C
C------------------------------------------------
C     Loop over number of excitations considered.
C------------------------------------------------
C
       inewtr = 1
c      DO 100 INEWTR = 1,NNEWTR
C
C---------------------------------------------------
C        Calculate RPA-density matrices in AO basis.
C---------------------------------------------------
C
c         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,NEXCI+INEWTR)
C
         DTIME     = SECOND()
         CALL SO_AODENS(WORK(KDENS),LDENS,WORK(KCMO),LCMO,
     &                  TR1E,LTR1E,TR1D,LTR1D,ISYMTR,
     &                  WORK(KEND2),LWORK2)
         DTIME     = SECOND()  - DTIME
         SOTIME(6) = SOTIME(6) + DTIME
C
         CALL SO_WRITE(WORK(KDENS),LDENS,LUDENS,FNDENS,INEWTR)
C
C--------------------------------------------
C        Backtransformation of trial vectors.
C--------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_BCKTR(TR1E,LTR1E,TR1D,LTR1D,WORK(KBTR1E),
     &                 LBTR1E,WORK(KBTR1D),LBTR1D,WORK(KBTJ1E),LBTJ1E,
     &                 WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,ISYMTR)
         DTIME     = SECOND()  - DTIME
         SOTIME(7) = SOTIME(7) + DTIME
C
         CALL SO_WRITE(WORK(KBTR1E),LBTR1E,LUBT1E,FNBT1E,INEWTR)
         CALL SO_WRITE(WORK(KBTR1D),LBTR1D,LUBT1D,FNBT1D,INEWTR)
         CALL SO_WRITE(WORK(KBTJ1E),LBTJ1E,LUBJ1E,FNBJ1E,INEWTR)
         CALL SO_WRITE(WORK(KBTJ1D),LBTJ1D,LUBJ1D,FNBJ1D,INEWTR)
C
c  100 CONTINUE
C
c      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
c      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
c      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
c      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
      CALL SO_CLOSE(LUDENS,FNDENS,'KEEP')
      CALL SO_CLOSE(LUBT1E,FNBT1E,'KEEP')
      CALL SO_CLOSE(LUBT1D,FNBT1D,'KEEP')
      CALL SO_CLOSE(LUBJ1E,FNBJ1E,'KEEP')
      CALL SO_CLOSE(LUBJ1D,FNBJ1D,'KEEP')
C
C-----------------------------------------------------------
C     Initialize RES1E, RES1D, RES2E, RES2D, SIGAI1, SIGAI2,
C                SIGDA1, SIGDA2 and FOCK
C-----------------------------------------------------------
C
c      I1 = NEXCI + 1
c      I2 = NEXCI + NNEWTR
C
c      CALL SO_INITIAL(I1,I2,WORK(KRES1E),LRES1E,LURS1E,FNRS1E)
c      CALL SO_INITIAL(I1,I2,WORK(KRES1D),LRES1D,LURS1D,FNRS1D)
c      CALL SO_INITIAL(I1,I2,WORK(KRES2E),LRES2E,LURS2E,FNRS2E)
c      CALL SO_INITIAL(I1,I2,WORK(KRES2D),LRES2D,LURS2D,FNRS2D)
C
      I1 = 1
      I2 = NNEWTR
C
      CALL SO_INITIAL(I1,I2,WORK(KSIGAI1),LSIGAI1,LUSAI1,FNSAI1)
      CALL SO_INITIAL(I1,I2,WORK(KSIGAI2),LSIGAI2,LUSAI2,FNSAI2)
      CALL SO_INITIAL(I1,I2,WORK(KSIGDA1),LSIGDA1,LUSDA1,FNSDA1)
      CALL SO_INITIAL(I1,I2,WORK(KSIGDA2),LSIGDA2,LUSDA2,FNSDA2)
      CALL SO_INITIAL(I1,I2,WORK(KFOCK),LFOCK,LUFOCK,FNFOCK)
C
C----------------------------
C     Initialize AIJ and AAB.
C----------------------------
C
      CALL DZERO(WORK(KAIJ),LAIJ)
      CALL DZERO(WORK(KAAB),LAAB)
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      IF (DIRECT) THEN
         NTOSYM = 1
         DTIME     = SECOND()
C        CALL HERDI1(WORK(KEND2),LWORK2,IPRINT)
         CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &               KODPP1,KODPP2,KRDPP1,KRDPP2,
     &               KFREE,LFREE,KEND2,
     &               WORK(KEND2),LWORK2,IPRERI)
         KEND2  = KFREE
         LWORK2 = LFREE
         DTIME     = SECOND()  - DTIME
         SOTIME(8) = SOTIME(8) + DTIME
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV  = KEND2
      LWORKSV = LWORK2
C
      ICDEL1 = 0
C
      DO 200 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            NTOT = MXCALL
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 210 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND2  = KENDSV
               LWORK2 = LWORKSV
C
               DTIME     = SECOND()
C              CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
C    &                          IPRINT)
               CALL ERIDI2(ILLL,INDEXA,NUMDIS,
     &                     WORK(KODCL1),WORK(KODCL2),
     &                     WORK(KODBC1),WORK(KODBC2),
     &                     WORK(KRDBC1),WORK(KRDBC2),
     &                     WORK(KODPP1),WORK(KODPP2),
     &                     WORK(KRDPP1),WORK(KRDPP2),
     &                     WORK(KEND2), LWORK2,IPRERI)
               DTIME     = SECOND()  - DTIME
               SOTIME(9) = SOTIME(9) + DTIME
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 220 IDEL2 = 1,NUMDIS
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
C              Work space allocation no. 3.
C------------------------------------------
C
               LXINT  = NDISAO(ISYDIS)
               LRECNR = ( (NBUF -1) / IRAT ) + 1
C
               KXINT   = KEND2
               KRECNR  = KXINT  + LXINT
               KEND3   = KRECNR + LRECNR
               LWORK3  = LWORK  - KEND3
C
               CALL SO_MEMMAX ('SO_RES_CP.3',LWORK3)
               IF (LWORK3 .LT. 0) 
     &             CALL STOPIT('SO_RES_CP.3',' ',KEND3,LWORK)
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               DTIME      = SECOND()
CSPAS:24/11-2003:I am not sure whether it should be IDEL2 here or not,
Cbecause it is a new argument
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                     WORK(KRECNR),DIRECT)
               DTIME      = SECOND()   - DTIME
               SOTIME(10) = SOTIME(10) + DTIME
C
C------------------------------------------
C              Work space allocation no. 4.
C------------------------------------------
C
               ISAIJ = MULD2H(ISYMD,1)
C
               IF (NVIR(ISYMD) .GT. 0) THEN
                  LT2M1 = NT2BCD(ISAIJ)
               ELSE
                  LT2M1 = 0
               END IF
C
               KT2M1   = KEND3
	       KEND4   = KT2M1  + LT2M1
               LWORK4  = LWORK  - KEND4
C
               CALL SO_MEMMAX ('SO_RES_CP.4',LWORK4)
               IF (LWORK4 .LT. 0) 
     &             CALL STOPIT('SO_RES_CP.4',' ',KEND4,LWORK)
C
C---------------------------------------------------------
C              Construct the partially back-transformed T2
C              MP-amplitudes.
C---------------------------------------------------------
C
               DTIME      = SECOND()
               CALL SO_T2M1(WORK(KT2M1),LT2M1,T2MP,LT2MP,WORK(KCMO),
     &                      LCMO,IDEL,ISYMD,ISYDIS,WORK(KEND4),LWORK4)
               DTIME      = SECOND()   - DTIME
               SOTIME(12) = SOTIME(12) + DTIME
C
C------------------------------------------
C              Work space allocation no. 5.
C------------------------------------------
C
               LDSRHF = NDSRHF(ISYMD)
C
               KDSRHF  = KEND4
               KEND5   = KDSRHF + LDSRHF
               LWORK5  = LWORK  - KEND5
C
               CALL SO_MEMMAX ('SO_RES_CP.5',LWORK5)
               IF (LWORK5 .LT. 0) 
     &             CALL STOPIT('SO_RES_CP.5',' ',KEND5,LWORK)
C
C----------------------------------------------------------------------
C              Transform one index in the integral batch to an occupied
C              index.
C----------------------------------------------------------------------
C
               DTIME      = SECOND()
               ISYMLP = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     &                     ISYMLP,WORK(KEND5),LWORK5,ISYDIS)
               DTIME      = SECOND()   - DTIME
               SOTIME(13) = SOTIME(13) + DTIME
C
C---------------------------------------------------
C              Noddy code to calculate MO-integrals.
C---------------------------------------------------
C
ckeld          if (sonod) then
ckeld             call so_nod(work(kxint),lxint,work(klamdp),llamdp,
ckeld&                        idel,rint)
ckeld          end if
C
C-------------------------------------------------------
C              Open files to store intermediate results.
C-------------------------------------------------------
C
               CALL SO_OPEN(LUDENS,FNDENS,LDENS)
               CALL SO_OPEN(LUFOCK,FNFOCK,LFOCK)
c               CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
c               CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
c               CALL SO_OPEN(LUTR1E,FNTR1E,FNTR1E,LTR1E)
c               CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
c               CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
c               CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
c               CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
c               CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
               CALL SO_OPEN(LUBT1E,FNBT1E,LBTR1E)
               CALL SO_OPEN(LUBT1D,FNBT1D,LBTR1D)
               CALL SO_OPEN(LUBJ1E,FNBJ1E,LBTJ1E)
               CALL SO_OPEN(LUBJ1D,FNBJ1D,LBTJ1D)
               CALL SO_OPEN(LUSAI1,FNSAI1,LSIGAI1)
               CALL SO_OPEN(LUSAI2,FNSAI2,LSIGAI2)
               CALL SO_OPEN(LUSDA1,FNSDA1,LSIGDA1)
               CALL SO_OPEN(LUSDA2,FNSDA2,LSIGDA2)
C
C---------------------------------------------------------
C              Loop over number of excitations considered.
C---------------------------------------------------------
C
               inewtr = 1
c               DO 230 INEWTR = 1,NNEWTR
C
C----------------------------------------------
C                 Calculate the AO-Fock matrix.
C----------------------------------------------
C
                  CALL SO_READ(WORK(KDENS),LDENS,LUDENS,FNDENS,INEWTR)
                  CALL SO_READ(WORK(KFOCK),LFOCK,LUFOCK,FNFOCK,INEWTR)
C
                  DTIME      = SECOND()
                  CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),WORK(KFOCK),
     &                           WORK(KEND5),LWORK5,IDEL,ISYMD,.FALSE.,
     &                           DUMMY,ISYMTR)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(11) = SOTIME(11) + DTIME
C
                  CALL SO_WRITE(WORK(KFOCK),LFOCK,LUFOCK,FNFOCK,INEWTR)
C
C----------------------------------------------------------------------
C                 Calculate part of the result vectors RES1E and RES1D,
C                 specifically the first and the second term in eqs.
C                 (34,35). Also calculate Aij and Aab in eqs. (43,44).
C----------------------------------------------------------------------
C
c                  CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,
c    &                          NEXCI+INEWTR)
c                  CALL SO_READ(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,
c    &                          NEXCI+INEWTR)
c                  CALL SO_READ(WORK(KTR1E), LTR1E, LUTR1E,FNTR1E,
c    &                          NEXCI+INEWTR)
c                  CALL SO_READ(WORK(KTR1D), LTR1D, LUTR1D,FNTR1D,
c    &                          NEXCI+INEWTR)
C
                  DTIME      = SECOND()
                  CALL SO_RES_A(RES1E,LRES1E,RES1D,LRES1D,
     &                          TR1E,LTR1E,TR1D,LTR1D,
     &                          WORK(KDSRHF),LDSRHF,WORK(KCMO),LCMO,
     &                          WORK(KT2M1),LT2M1,WORK(KAIJ),LAIJ,
     &                          WORK(KAAB),LAAB,INEWTR,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(14) = SOTIME(14) + DTIME
C
C-------------------------------------------------------------------
C                 Calculate the part of the result vectors RES1E and
C                 RES1D which orriginate from the C matrices. See 
C                 eqs. (72) and (73).
C-------------------------------------------------------------------
C
c                  CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,
c    &                          NEXCI+INEWTR)
c                  CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,
c    &                          NEXCI+INEWTR)
C
                  DTIME      = SECOND()
                  CALL SO_RES_TCB(RES1E,LRES1E,RES1D,
     &                            LRES1D,TR2E,LTR2E,
     &                            TR2D,LTR2D,WORK(KDSRHF),LDSRHF,
     &                            WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                            ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(29) = SOTIME(29) + DTIME
C
c                  CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,
c    &                           NEXCI+INEWTR)
c                  CALL SO_WRITE(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,
c    &                           NEXCI+INEWTR)
C
C----------------------------------------------------------------------
C                 Construct C-contribution to 2p2h result vectors RES2E
C                 and RES2D.
C----------------------------------------------------------------------
C
c                  CALL SO_READ(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,
c    &                          NEXCI+INEWTR)
c                  CALL SO_READ(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,
c    &                          NEXCI+INEWTR)
                  CALL SO_READ(WORK(KBTR1E),LBTR1E,LUBT1E,FNBT1E,INEWTR)
                  CALL SO_READ(WORK(KBTR1D),LBTR1D,LUBT1D,FNBT1D,INEWTR)
                  CALL SO_READ(WORK(KBTJ1E),LBTJ1E,LUBJ1E,FNBJ1E,INEWTR)
                  CALL SO_READ(WORK(KBTJ1D),LBTJ1D,LUBJ1D,FNBJ1D,INEWTR)
C
                  DTIME      = SECOND()
                  CALL SO_RES_CB(RES2E,LRES2E,RES2D,
     &                          LRES2D,
     &                          WORK(KDSRHF),LDSRHF,WORK(KBTR1E),LBTR1E,
     &                          WORK(KBTR1D),LBTR1D,WORK(KBTJ1E),LBTJ1E,
     &                          WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,
     &                          IDEL,ISYMD,ISYDIS,ISYMTR,WORK(KEND5),
     &                          LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(15) = SOTIME(15) + DTIME
C
c                  CALL SO_WRITE(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,
c     &                          NEXCI+INEWTR)
c                  CALL SO_WRITE(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,
c     &                          NEXCI+INEWTR)
C
C--------------------------------------------------------------------
C                 Construct SIGMAI1(ALFA,I) and SIGMAI2(ALFA,I) which
C                 are used in SO_RES_B.
C--------------------------------------------------------------------
C
                  CALL SO_READ(WORK(KSIGAI1),LSIGAI1,LUSAI1,FNSAI1,
     &                         INEWTR)
                  CALL SO_READ(WORK(KSIGAI2),LSIGAI2,LUSAI2,FNSAI2,
     &                         INEWTR)
C
                  DTIME      = SECOND()
                  CALL SO_SIGAI(WORK(KSIGAI1),LSIGAI1,WORK(KSIGAI2),
     &                          LSIGAI2,WORK(KT2M1),LT2M1,WORK(KXINT),
     &                          LXINT,WORK(KBTR1E),LBTR1E,WORK(KBTR1D),
     &                          LBTR1D,WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),
     &                          LBTJ1D,WORK(KCMO),LCMO,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(16) = SOTIME(16) + DTIME
C
                  CALL SO_WRITE(WORK(KSIGAI1),LSIGAI1,LUSAI1,FNSAI1,
     &                          INEWTR)
                  CALL SO_WRITE(WORK(KSIGAI2),LSIGAI2,LUSAI2,FNSAI2,
     &                          INEWTR)
C
C--------------------------------------------------------------------
C                 Construct SIGDA1(DELTA,A) and SIGDA2(DELTA,A) which
C                 are used SO_RES_C.
C--------------------------------------------------------------------
C
                  CALL SO_READ(WORK(KSIGDA1),LSIGDA1,LUSDA1,FNSDA1,

     &                         INEWTR)
                  CALL SO_READ(WORK(KSIGDA2),LSIGDA2,LUSDA2,FNSDA2,
     &                         INEWTR)
C
                  DTIME      = SECOND()
                  CALL SO_SIGDA(WORK(KSIGDA1),LSIGDA1,WORK(KSIGDA2),
     &                          LSIGDA2,T2MP,LT2MP,WORK(KDSRHF),LDSRHF,
     &                          WORK(KBTR1E),LBTR1E,WORK(KBTR1D),LBTR1D,
     &                          WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),LBTJ1D,
     &                          WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                          ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(17) = SOTIME(17) + DTIME
C
                  CALL SO_WRITE(WORK(KSIGDA1),LSIGDA1,LUSDA1,FNSDA1,
     &                          INEWTR)
                  CALL SO_WRITE(WORK(KSIGDA2),LSIGDA2,LUSDA2,FNSDA2,
     &                          INEWTR)
C
c  230          CONTINUE
C
C--------------------------
C              Close files.
C--------------------------
C
               CALL SO_CLOSE(LUDENS,FNDENS,'KEEP')
               CALL SO_CLOSE(LUFOCK,FNFOCK,'KEEP')
c               CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
c               CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
c               CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
c               CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
c               CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
c               CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
c               CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
c               CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
               CALL SO_CLOSE(LUBT1E,FNBT1E,'KEEP')
               CALL SO_CLOSE(LUBT1D,FNBT1D,'KEEP')
               CALL SO_CLOSE(LUBJ1E,FNBJ1E,'KEEP')
               CALL SO_CLOSE(LUBJ1D,FNBJ1D,'KEEP')
               CALL SO_CLOSE(LUSAI1,FNSAI1,'KEEP')
               CALL SO_CLOSE(LUSAI2,FNSAI2,'KEEP')
               CALL SO_CLOSE(LUSDA1,FNSDA1,'KEEP')
               CALL SO_CLOSE(LUSDA2,FNSDA2,'KEEP')
C
  220       CONTINUE
C
  210    CONTINUE
C
  200 CONTINUE
C
C
C=================================================
C     End of loop over distributions of integrals.
C=================================================
C
C
C----------------
C     Open files.
C----------------
C
      CALL SO_OPEN(LUFOCK,FNFOCK,LFOCK)
      CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
c      CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
c      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
c      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
c      CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
c      CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
c      CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
c      CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
      CALL SO_OPEN(LUSAI1,FNSAI1,LSIGAI1)
      CALL SO_OPEN(LUSAI2,FNSAI2,LSIGAI2)
      CALL SO_OPEN(LUSDA1,FNSDA1,LSIGDA1)
      CALL SO_OPEN(LUSDA2,FNSDA2,LSIGDA2)
c      CALL SO_OPEN(LURO1E,FNRO1E,LRESO1E)
c      CALL SO_OPEN(LURO1D,FNRO1D,LRESO1D)
C
C------------------------------------------------
C     Loop over number of excitations considered.
C------------------------------------------------
C
       inewtr = 1
c      DO 300 INEWTR = 1,NNEWTR
C
C---------------------------------------------
C        Transform AO Fock matrix to MO basis.
C---------------------------------------------
C
         CALL SO_READ(WORK(KFOCK),LFOCK,LUFOCK,FNFOCK,INEWTR)
C
         DTIME      = SECOND()
         CALL CC_FCKMO(WORK(KFOCK),WORK(KCMO),WORK(KCMO),
     &                    WORK(KEND2),LWORK2,ISYRES,1,1)
         DTIME      = SECOND()   - DTIME
         SOTIME(24) = SOTIME(24) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA two-particle parts to the result 
C        vectors.
C------------------------------------------------------------------
C
c         CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,NEXCI+INEWTR)
C
         DTIME      = SECOND()
         CALL SO_TWOFOCK(RES1E,LRES1E,RES1D,LRES1D,
     &                   WORK(KFOCK),LFOCK,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(25) = SOTIME(25) + DTIME
C
C----------------------------------------------------------------------
C        Add contribution from sigma1(alfa,i) to RES1E and from 
C        sigma2(alfa,i) to RES1D. I.e. the third terms in eqs. (34) and 
C        (35).
C----------------------------------------------------------------------
C
         CALL SO_READ(WORK(KSIGAI1),LSIGAI1,LUSAI1,FNSAI1,INEWTR)
         CALL SO_READ(WORK(KSIGAI2),LSIGAI2,LUSAI2,FNSAI2,INEWTR)
C
         DTIME      = SECOND()
         CALL SO_RES_B(RES1E,LRES1E,RES1D,LRES1D,
     &                 WORK(KSIGAI1),LSIGAI1,WORK(KSIGAI2),LSIGAI2,
     &                 WORK(KCMO),LCMO,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(18) = SOTIME(18) + DTIME
C
C--------------------------------------------------------------------
C        Add contribution from sigda1(delta,a) to RES1E and from
C        sigda2(delta,a) to RES1D. I.e. the fourth terms in eqs. (34)
C        and (35).
C--------------------------------------------------------------------
C
         CALL SO_READ(WORK(KSIGDA1),LSIGDA1,LUSDA1,FNSDA1,INEWTR)
         CALL SO_READ(WORK(KSIGDA2),LSIGDA2,LUSDA2,FNSDA2,INEWTR)
C
         DTIME      = SECOND()
         CALL SO_RES_C(RES1E,LRES1E,RES1D,LRES1D,
     &                 WORK(KSIGDA1),LSIGDA1,WORK(KSIGDA2),LSIGDA2,
     &                 WORK(KCMO),LCMO,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(19) = SOTIME(19) + DTIME
C
C--------------------------------------------------------------
C        Calculate and add the symmetry correcting term to A in
C        eq. (44).
C--------------------------------------------------------------
C
c         CALL SO_READ(WORK(KTR1E), LTR1E, LUTR1E,FNTR1E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KTR1D), LTR1D, LUTR1D,FNTR1D,NEXCI+INEWTR)
C
         DTIME      = SECOND()
         CALL SO_RES_SYM(RES1E,LRES1E,RES1D,LRES1D,
     &                   WORK(KAIJ),LAIJ,WORK(KAAB),LAAB,TR1E,
     &                   LTR1E,TR1D,LTR1D,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(20) = SOTIME(20) + DTIME
C
C---------------------------------------------------------
C        Calculate and add the Fock-term to A in eq. (40).
C---------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_RES_FCK(RES1E,LRES1E,RES1D,LRES1D,
     &                   TR1E,LTR1E,TR1D,
     &                   LTR1D,FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,
     &                   LDENSAB,ISYRES,ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(21) = SOTIME(21) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA one-particle parts to the result
C        vectors.
C------------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_ONEFOCK(RES1E,LRES1E,RES1D,LRES1D,FOCKD,
     &                   LFOCKD,TR1E,LTR1E,TR1D,LTR1D,
     &                   ISYRES,ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(26) = SOTIME(26) + DTIME
C
c         CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,NEXCI+INEWTR)
c         CALL SO_WRITE(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,NEXCI+INEWTR)
C
C-------------------------------------------------
C        Calculate the overlap matrix in eq. (67).
C-------------------------------------------------
C
c         CALL DZERO(WORK(KRESO1E),LRESO1E)
c         CALL DZERO(WORK(KRESO1D),LRESO1D)
C
         DTIME      = SECOND()
         CALL SO_RES_OVL(RESO1E,LRESO1E,RESO1D,LRESO1D,
     &                   TR1E,LTR1E,TR1D,LTR1D,DENSIJ,
     &                   LDENSIJ,DENSAB,LDENSAB,ISYRES,ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(22) = SOTIME(22) + DTIME
C
C-----------------------------------------
C        Calculate the RPA overlap matrix.
C-----------------------------------------
C
         DTIME      = SECOND()
         CALL SO_RES_OVLR(RESO1E,LRESO1E,RESO1D,LRESO1D,
     &                    TR1E,LTR1E,TR1D,LTR1D)
         DTIME      = SECOND()   - DTIME
         SOTIME(23) = SOTIME(23) + DTIME
C
c        CALL SO_WRITE(WORK(KRESO1E),LRESO1E,LURO1E,FNRO1E,NEXCI+INEWTR)
c        CALL SO_WRITE(WORK(KRESO1D),LRESO1D,LURO1D,FNRO1D,NEXCI+INEWTR)
C
C-----------------------------------------------------------------
C        Construct D-contribution to 2p2h result vectors RES2E and
C        RES2D.
C-----------------------------------------------------------------
C
c         CALL SO_READ(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KTR2E), LTR2E, LUTR2E,FNTR2E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KTR2D), LTR2D, LUTR2D,FNTR2D,NEXCI+INEWTR)
C
         DTIME      = SECOND()
         CALL SO_RES_CD(RES2E,LRES2E,RES2D,LRES2D,
     &                  TR2E,LTR2E,TR2D,LTR2D,
     &                  FOCKD,LFOCKD,ISYRES,WORK(KEND2),LWORK2)
         DTIME      = SECOND()   - DTIME
         SOTIME(30) = SOTIME(30) + DTIME
C
c         CALL SO_WRITE(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,NEXCI+INEWTR)
c         CALL SO_WRITE(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,NEXCI+INEWTR)
C
c  300 CONTINUE
C
C---------------------------------------
C     Write new resultvectors to output.
C---------------------------------------
C
c      DO 400 INEWTR = 1,NNEWTR
C
c         CALL AROUND('New linear transformed trial vector')
C
c         CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,NEXCI+INEWTR)
c         CALL SO_READ(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,NEXCI+INEWTR)
C
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KRES1E+I-1),I=1,LRES1E)
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KRES2E+I-1),I=1,LRES2E)
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KRES1D+I-1),I=1,LRES1D)
c         WRITE(LUPRI,'(I8,1X,F14.8)') (I,WORK(KRES2D+I-1),I=1,LRES2D)
C
c  400 CONTINUE
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUFOCK,FNFOCK,'DELETE')
c      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
c      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
c      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
c      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
c      CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
c      CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
c      CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
c      CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
      CALL SO_CLOSE(LUSAI1,FNSAI1,'DELETE')
      CALL SO_CLOSE(LUSAI2,FNSAI2,'DELETE')
      CALL SO_CLOSE(LUSDA1,FNSDA1,'DELETE')
      CALL SO_CLOSE(LUSDA2,FNSDA2,'DELETE')
c      CALL SO_CLOSE(LURO1E,FNRO1E,'KEEP')
c      CALL SO_CLOSE(LURO1D,FNRO1D,'KEEP')
C
ckeld if (sonod) then
ckeld    call so_nod2(work(ktr1e),ltr1e,work(ktr1d),ltr1d,t2mp,lt2mp,
ckeld&                rint,work(kend2),lwork2)
ckeld    call so_nod3(rint,work(kend4),lwork4)
ckeld    if (work(ktr1e) .gt. 0.5d0) then
ckeld       call so_nod4(rint,work(kend4),lwork4)
ckeld    end if
ckeld end if
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_RES_CP')
C
      RETURN
      END
