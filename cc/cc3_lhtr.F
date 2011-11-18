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
C  /* Deck cc3_lhtr_t3 */
      SUBROUTINE CC3_T3_LHTR(ECURR,OMEGA1,T1AM,ISYMT1,T2TP,ISYMT2,C2TP,
     *                       ISYMC2,XLAMDP,XLAMDH,WORK,LWORK,
     *                       LU3SRT,FN3SRT,LUCKJD,FNCKJD)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate the T3 dependent terms that contribute to
C     the left hand side eigenvalue equation.
C
C     ISYMT1 is symmetry of T1AM
C     ISYMT2 is symmetry of T2TP
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <second.h>
C
      INTEGER ISYMT1, ISYMT2, ISYMC2, LWORK
      INTEGER ISYMTR, ISYRES, ISINT2, ISYMIM, KFOCKD, KEND0, LWRK0, KCMO
      INTEGER KTROC0, KXIAJB, KXIAJB2, KEND1, LWRK1, KEND2, LWRK2
      INTEGER LENGTH, IOFF, ISYMD, ISAIJ1, ISYCKB, ISCKB2, KINTOC
      INTEGER KTRVI0, KTRVI2, KTRVI3, KEND3, LWRK3, KEND4, LWRK4
      INTEGER KINTVI, ISYMB, ISYALJ, ISAIJ2, ISYMBD, ISCKIJ
      INTEGER KSMAT, KQMAT, KDIAG, KINDSQ, KINDEX, KTMAT
      INTEGER LENSQ, ISYOPE, IOPTTCME
      INTEGER LU3SRT, LUCKJD, LUDELD, LUDKBC
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), T1AM(*), T2TP(*), C2TP(*)
      REAL XLAMDH(*), XLAMDP(*), WORK(LWORK)
      REAL TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      REAL DTIME, XNORM, DDOT, ECURR
#else
      DOUBLE PRECISION OMEGA1(*), T1AM(*), T2TP(*), C2TP(*)
      DOUBLE PRECISION XLAMDH(*), XLAMDP(*), WORK(LWORK)
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      DOUBLE PRECISION DTIME, XNORM, DDOT, ECURR
#endif
C
      CHARACTER*(*) FN3SRT, FNCKJD
      CHARACTER*1 CDUMMY
      CHARACTER*11 FNDELD, FNDKBC
C
      CALL QENTER('CC3_T3_LHTR')
C
C-----------------------------
C     Open files.
C-----------------------------
C
      CDUMMY = ' '
C
      LUDELD = -1
      LUDKBC = -1
      FNDELD = 'CC3_T3_TMP1'
      FNDKBC = 'CC3_T3_TMP2'
C
      CALL WOPEN2(LUDELD,FNDELD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
c
*     write(lupri,*)'t1am in cc3_lhtr'
*     call print_matai(t1am,1)
C
C-------------------------------------------------------------
C     Set symmetry flags.
C
C     isymres is symmetry of result(omega)
C     isint2 is symmetry of integrals in the triples equation.(int2)
C     isymim is symmetry of S and Q intermediates.(t2*int2)
C      (sym is for all index of S and Q (cbd,klj)
C       thus cklj=b*d*isymim)
C-------------------------------------------------------------
C
      IPRCC = IPRINT
      ISYMTR = ISYMC2
      ISYRES = MULD2H(ISYMTR,ISYMOP)
      ISINT2 = ISYMOP
      ISYMIM = ISYMOP
C
C--------------------
C     Time variables.
C--------------------
C
      TITRAN = 0.0D0
      TISORT = 0.0D0
      TISMAT = 0.0D0
      TIQMAT = 0.0D0
      TIOME1 = 0.0D0
C
C---------------------------------------------------------
C     Transform and sort qmat integrals to smat integrals.
C---------------------------------------------------------
C
      CALL CC3_SORT1(WORK,LWORK,2,ISINT2,LU3SRT,FN3SRT,
     *               LUDELD,FNDELD,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
      CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT2,LUDELD,FNDELD,
     *              LUDKBC,FNDKBC)
C
C--------------------------------------
C     Reorder the t2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP in CC3_LHTR')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMT2),T2TP,1,WORK,1)
      CALL CC3_T2TP(T2TP,WORK,ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NT2SQ(ISYMT2),T2TP,1,T2TP,1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XNORM
      ENDIF
C
C--------------------------------------
C     Reorder the C2-amplitudes i C2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMC2)) THEN
         CALL QUIT('Not enough memory to construct C2TP in CC3_LHTR')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMC2),C2TP,1,WORK,1)
      CALL CC3_T2TP(C2TP,WORK,ISYMC2)
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NT2SQ(ISYMC2),C2TP,1,C2TP,1)
         WRITE(LUPRI,*) 'Norm of C2TP ',XNORM
      ENDIF
C
C---------------------------------------------------------
C     Read canonical orbital energies
C---------------------------------------------------------
C
      KFOCKD = 1
      KCMO   = KFOCKD + NORBTS
      KEND0  = KCMO   + NLAMDS
      LWRK0  = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    > ',KEND0
         CALL QUIT('Insufficient space in CCSDT_OMEG')
      END IF
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND0),LWRK0)
C
C---------------------------------------------
C     Delete frozen orbitals in Fock diagonal.
C---------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC0  = KEND0
      KXIAJB  = KTROC0  + NTRAOC(ISINT2)
      KXIAJB2 = KXIAJB  + NT2AM(ISYMOP)
      KEND1   = KXIAJB2 + NT2AM(ISYMOP)
      LWRK1   = LWORK   - KEND1
C
      KINTOC = KEND1
      KEND2  = KINTOC + MAX(NTOTOC(ISYMOP),NTOTOC(ISINT2))
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    > ',KEND2
         CALL QUIT('Insufficient space in CCSDT_OMEG')
      END IF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      LENGTH = IRAT*NT2AM(ISYMOP)
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
C
      CALL DCOPY(NT2AM(ISYMOP),WORK(KXIAJB),1,WORK(KXIAJB2),1)
      ISYOPE = ISYMOP
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KXIAJB2),1.0D0,ISYOPE,IOPTTCME)
C
      IF ( IPRINT .GT. 55) THEN
         XNORM = DDOT(NT2AM(ISYMOP),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of g-IAJB ',XNORM
         XNORM = DDOT(NT2AM(ISYMOP),WORK(KXIAJB2),1,
     *                WORK(KXIAJB2),1)
         WRITE(LUPRI,*) 'Norm of L-IAJB ',XNORM
      ENDIF
C
C-----------------------
C     Read in integrals.
C-----------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C     Can use WORK(KCMO) here!!!
C----------------------------------------------------------------------
C
      DTIME  = SECOND()
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KCMO),
     *                 WORK(KEND2),LWRK2,ISINT2)
C
      DTIME  = SECOND() - DTIME
      TITRAN = TITRAN   + DTIME
C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 CC3_T3_LHTR : ',XNORM
      ENDIF
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM
C
         ISAIJ1 = MULD2H(ISYMD,ISYRES)
         ISYCKB = MULD2H(ISYMD,ISYMOP)
         ISCKB2 = MULD2H(ISINT2,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISAIJ1:',ISAIJ1
            WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYCKB:',ISYCKB
            WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISCKB2:',ISCKB2
C
         ENDIF
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KTRVI2 = KEND1
         KEND2  = KTRVI2 + NCKATR(ISCKB2)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0 = KEND2
         KTRVI3 = KTRVI0 + NCKATR(ISCKB2)
         KEND3  = KTRVI3 + NCKATR(ISCKB2)
         LWRK3  = LWORK  - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISYCKB),NCKA(ISYMD),NCKA(ISCKB2))
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_LHTR_L3')
         END IF
C
         DO D = 1,NVIR(ISYMD)
C
C-------------------------------------------------------
C           Read virtual integrals used in s3am.
C-------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                     NCKATR(ISCKB2))
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XNORM= DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 CC3_T3_LHTR : ',XNORM
            ENDIF
C
C---------------------------------------
C           Sort the integrals for s3am.
C---------------------------------------
C
            DTIME = SECOND()
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55) THEN
               XNORM= DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XNORM
            ENDIF
C
C-----------------------------------------------
C           Read virtual integrals used in q3am.
C-----------------------------------------------
C
            IOFF = ICKAD(ISCKB2,ISYMD) + NCKA(ISCKB2)*(D - 1) + 1
            IF (NCKA(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     &                     NCKA(ISCKB2))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),XLAMDH,
     *                       ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_LHTR_L3')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND3),ISYMD,D,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55) THEN
               XNORM= DDOT(NCKATR(ISCKB2),WORK(KTRVI3),1,
     *                      WORK(KTRVI3),1)
               WRITE(LUPRI,*) 'Norm of TRVI3 CC3_T3_LHTR : ',XNORM
            ENDIF
C
C---------------------
C           Calculate.
C---------------------
C
            DO ISYMB = 1,NSYM
C
               ISYALJ = MULD2H(ISYMB,ISYMT2)
               ISAIJ2 = MULD2H(ISYMB,ISYRES)
               ISYMBD = MULD2H(ISYMB,ISYMD)
               ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
               IF (IPRINT .GT. 55) THEN
C
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISCKIJ:',ISCKIJ
C
               ENDIF
C
               KSMAT  = KEND3
               KQMAT  = KSMAT  + NCKIJ(ISCKIJ)
               KDIAG  = KQMAT  + NCKIJ(ISCKIJ)
               KINDSQ = KDIAG  + NCKIJ(ISCKIJ)
               KINDEX = KINDSQ + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KTMAT  = KINDEX + (NCKI(ISYALJ) - 1)/IRAT + 1
               KEND4  = KTMAT  + NCKIJ(ISCKIJ)
               LWRK4  = LWORK  - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_LHTR_L3')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
C
               IF (IPRINT .GT. 55) THEN
                  XNORM  = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                    WORK(KDIAG),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XNORM
               ENDIF

C
C-------------------------------------
C              Construct index arrays.
C-------------------------------------
C
               LENSQ = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               CALL CC3_INDEX(WORK(KINDEX),ISYALJ)
C
               DO B = 1,NVIR(ISYMB)
C
C-------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3.
C-------------------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_SMAT(ECURR,T2TP,ISYMT2,WORK(KTMAT),
     *                          WORK(KTRVI0),
     *                          WORK(KTRVI2),WORK(KTROC0),ISINT2,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMAT),WORK(KEND4),LWRK4,
     *                          WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT  ',XNORM
                  ENDIF
C
C--------------------------------------------------
C                 Calculate Q(ci,jk) for fixed b,d.
C--------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_QMAT(ECURR,T2TP,ISYMT2,WORK(KTRVI3),
     *                          WORK(KTROC0),ISINT2,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KQMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT  ',XNORM
                  ENDIF
C
C------------------------------------
C                 Calculate Omega1.
C------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DSCAL(NCKIJ(ISCKIJ),-1.0D0,WORK(KSMAT),1)
                  CALL DSCAL(NCKIJ(ISCKIJ),-1.0D0,WORK(KQMAT),1)
C
                  CALL T3_ONEL1(OMEGA1,WORK(KSMAT),WORK(KQMAT),
     *                          WORK(KTMAT),ISYMIM,WORK(KXIAJB),
     *                          WORK(KXIAJB2),ISINT2,
     *                          C2TP,ISYMC2,WORK(KINDSQ),LENSQ,
     *                          WORK(KEND4),LWRK4,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORM = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                    WRITE(LUPRI,*) 'Norm of Rho1 -after T3_ONEL1',XNORM
                  ENDIF
C
C
                  CALL T3_ONEL2(OMEGA1,WORK(KSMAT),WORK(KQMAT),
     *                          WORK(KTMAT),ISYMIM,WORK(KXIAJB),
     *                          WORK(KXIAJB2),ISINT2,
     *                          C2TP,ISYMC2,WORK(KINDSQ),LENSQ,
     *                          WORK(KEND4),LWRK4,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORM = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                    WRITE(LUPRI,*) 'Norm of Rho1 -after T3_ONEL2',XNORM
                  ENDIF
C
C
                  CALL T3_ONEL3(OMEGA1,WORK(KSMAT),WORK(KQMAT),
     *                          WORK(KTMAT),ISYMIM,WORK(KXIAJB2),ISINT2,
     *                          C2TP,ISYMC2,WORK(KINDSQ),LENSQ,
     *                          WORK(KEND4),LWRK4,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORM = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                    WRITE(LUPRI,*) 'Norm of Rho1 -after T3_ONEL3',XNORM
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TIOME1 = TIOME1   + DTIME
C
C---------------------------------------------------------
C              End for B
C---------------------------------------------------------
C
               ENDDO        ! B
            ENDDO           ! ISYMB
C
C---------------------------------------------------
C           End for D
C---------------------------------------------------
C
         ENDDO     ! D
      ENDDO        ! ISYMD
C
C-----------------------------------------
C     Close and delete files
C-----------------------------------------
C
      CALL WCLOSE2(LUDELD,FNDELD,'DELETE')
      CALL WCLOSE2(LUDKBC,FNDKBC,'DELETE')
C
C-------------------
C     Print timings.
C-------------------
C
      IF (IPRINT .GT. 9) THEN
COMMENT COMMENT
COMMENT COMMENT  Have a look at the timings in this routine.
COMMENT COMMENT
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,1) 'CC3_TRAN  : ',TITRAN
         WRITE(LUPRI,1) 'CC3_SORT  : ',TISORT
         WRITE(LUPRI,1) 'CC3_SMAT  : ',TISMAT
         WRITE(LUPRI,1) 'CC3_QMAT  : ',TIQMAT
         WRITE(LUPRI,1) 'CC3_OME1  : ',TIOME1
         WRITE(LUPRI,*)
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_T3_LHTR')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck t3_onel1 */
      SUBROUTINE T3_ONEL1(OMEGA1,SMAT,QMAT,TMAT,ISYMIM,XIAJB,YIAJB,
     *                    ISYINT,C2TP,ISYMC2,INDSQ,LENSQ,WORK,LWORK,
     *                    ISYMB,B,ISYMD,D)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate the term t^{def}_{lmn} L^{fd}_{mi} g_{nela} 
C                      - t^{def}_{lnm} L^{fd}_{mi} L_{nela}
C
C     XIAJB contains g and YIAJB contains L
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYINT, ISYMC2, LENSQ, LWORK, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYRE1, ISYRES, ISYMBD, ISFLMN, ISYANL, LENGTH
      INTEGER ISYFIM, KTMAT, KC2TEMP, KINT, KEND1, LWRK1
      INTEGER ISYMM, ISYMFI, ISYMF, KOFF1, KOFF2
      INTEGER ISYML, ISYMAN, ISYMA, ISYMN, ISYMLN, NBN, NAN, NAL
      INTEGER ISYMI, ISYMAB, ISYMFM, KOFF3, NUMBFM, NUMBLN, NUMBA
      INTEGER ISYMBN, ISYMAL, KINT2
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), SMAT(*), QMAT(*), TMAT(*), XIAJB(*)
      REAL YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION OMEGA1(*), SMAT(*), QMAT(*), TMAT(*), XIAJB(*)
      DOUBLE PRECISION YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('T3_ONEL1')
C
      ISYRE1 = MULD2H(ISYMIM,ISYMC2)
      ISYRES = MULD2H(ISYRE1,ISYINT)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISFLMN = MULD2H(ISYMIM,ISYMBD)
      ISYANL = MULD2H(ISYMB,ISYINT)
C
      LENGTH = NCKIJ(ISFLMN)
C
C-----------------------------
C     Sort C2
C-----------------------------
C
      ISYFIM = MULD2H(ISYMC2,ISYMD)
C
      KTMAT   = 1
      KC2TEMP = KTMAT   + NCKIJ(ISFLMN)
      KINT    = KC2TEMP + NMAIJA(ISYFIM)
      KINT2   = KINT    + NCKI(ISYANL)
      KEND1   = KINT2   + NCKI(ISYANL)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in T3_ONEL1 (sort)')
      ENDIF
C
      DO ISYMM = 1, NSYM
         ISYMFI = MULD2H(ISYFIM,ISYMM)
         DO ISYMF = 1, NSYM
            ISYMI = MULD2H(ISYMFI,ISYMF)
            ISYMFM = MULD2H(ISYMF,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO I = 1, NRHF(ISYMI)
C
                  KOFF1 = IT2SP(ISYFIM,ISYMD)
     *                  + NCKI(ISYFIM)*(D-1)
     *                  + ICKI(ISYMFI,ISYMM)
     *                  + NT1AM(ISYMFI)*(M-1)
     *                  + IT1AM(ISYMF,ISYMI)
     *                  + NVIR(ISYMF)*(I-1) 
     *                  + 1
C
                  KOFF2 = KC2TEMP
     *                  + ICKI(ISYMFM,ISYMI)
     *                  + NT1AM(ISYMFM)*(I-1)
     *                  + IT1AM(ISYMF,ISYMM)
     *                  + NVIR(ISYMF)*(M-1) 

C
                  CALL DCOPY(NVIR(ISYMF),C2TP(KOFF1),1,WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------------
C     Sort integrals.
C---------------------------
C
      DO ISYML = 1, NSYM
         ISYMAN = MULD2H(ISYANL,ISYML)
         DO ISYMA = 1, NSYM
            ISYMN  = MULD2H(ISYMAN,ISYMA)
            ISYMLN = MULD2H(ISYMN,ISYML)
            ISYMBN = MULD2H(ISYMB,ISYMN)
            ISYMAL = MULD2H(ISYMA,ISYML)
C
            DO N = 1, NRHF(ISYMN)
               NBN = IT1AM(ISYMB,ISYMN) + NVIR(ISYMB)*(N-1) + B
               DO A = 1, NVIR(ISYMA)
                  NAN = IT1AM(ISYMA,ISYMN) + NVIR(ISYMA)*(N-1) + A
                  DO L = 1, NRHF(ISYML)
                     NAL = IT1AM(ISYMA,ISYML) + NVIR(ISYMA)*(L-1) + A
C
                     KOFF1 = IT2AM(ISYMBN,ISYMAL) + INDEX(NBN,NAL)
                     KOFF2 = KINT - 1
     *                     + IMAIJA(ISYMLN,ISYMA)
     *                     + NMATIJ(ISYMLN)*(A-1)
     *                     + IMATIJ(ISYML,ISYMN)
     *                     + NRHF(ISYML)*(N-1)
     *                     + L
                     KOFF3 = KINT2 - 1
     *                     + IMAIJA(ISYMLN,ISYMA)
     *                     + NMATIJ(ISYMLN)*(A-1)
     *                     + IMATIJ(ISYML,ISYMN)
     *                     + NRHF(ISYML)*(N-1)
     *                     + L
C
                     WORK(KOFF2) = XIAJB(KOFF1)
                     WORK(KOFF3) = YIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT for the g term
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,2))
     *           + QMAT(INDSQ(I,1))
C
      WORK(KTMAT-1+I) = SMAT(INDSQ(I,5))
     *                + QMAT(INDSQ(I,4))
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMAB = MULD2H(ISYMA,ISYMB)
         ISYMLN = MULD2H(ISYINT,ISYMAB)
         ISYMFM = MULD2H(ISYFIM,ISYMI)
C
         CALL DZERO(WORK(KEND1),NMATIJ(ISYMLN)*NRHF(ISYMI))
C
         KOFF1 = ISAIKL(ISYMFM,ISYMLN) + 1
         KOFF2 = KC2TEMP
     *         + ICKI(ISYMFM,ISYMI)
         KOFF3 = KEND1
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,TMAT(KOFF1),NUMBFM,WORK(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMFM,ISYMLN)
         KOFF2 = IT2SP(ISYFIM,ISYMD)
     *         + NCKI(ISYFIM)*(D-1)
     *         + ICKI(ISYMFM,ISYMI) + 1
         KOFF3 = KEND1
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,WORK(KOFF1),NUMBFM,C2TP(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KINT
     *         + IMAIJA(ISYMLN,ISYMA)
         KOFF2 = KEND1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMLN),
     *              ONE,WORK(KOFF1),NUMBLN,WORK(KOFF2),NUMBLN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
C
      ENDDO
C
C----------------------------------
C     Construct TMAT for L term
C----------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = - SMAT(INDSQ(I,1))
     *             - QMAT(INDSQ(I,2))
C
      WORK(KTMAT-1+I) = - SMAT(I)
     *                  - QMAT(INDSQ(I,3))
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMAB = MULD2H(ISYMA,ISYMB)
         ISYMLN = MULD2H(ISYINT,ISYMAB)
         ISYMFM = MULD2H(ISYFIM,ISYMI)
C
         CALL DZERO(WORK(KEND1),NMATIJ(ISYMLN)*NRHF(ISYMI))
C
         KOFF1 = ISAIKL(ISYMFM,ISYMLN) + 1
         KOFF2 = KC2TEMP
     *         + ICKI(ISYMFM,ISYMI)
         KOFF3 = KEND1
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,TMAT(KOFF1),NUMBFM,WORK(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMFM,ISYMLN)
         KOFF2 = IT2SP(ISYFIM,ISYMD)
     *         + NCKI(ISYFIM)*(D-1)
     *         + ICKI(ISYMFM,ISYMI) + 1
         KOFF3 = KEND1
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,WORK(KOFF1),NUMBFM,C2TP(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KINT2
     *         + IMAIJA(ISYMLN,ISYMA)
         KOFF2 = KEND1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMLN),
     *              ONE,WORK(KOFF1),NUMBLN,WORK(KOFF2),NUMBLN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
C
      ENDDO
C
C----------------------------
C     End.
C----------------------------
C
      CALL QEXIT('T3_ONEL1')
C
      RETURN
      END
C  /* Deck t3_onel2 */
      SUBROUTINE T3_ONEL2(OMEGA1,SMAT,QMAT,TMAT,ISYMIM,XIAJB,YIAJB,
     *                    ISYINT,C2TP,ISYMC2,INDSQ,LENSQ,WORK,LWORK,
     *                    ISYMB,B,ISYMD,D)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate the term t^{def}_{lmn} L^{ad}_{mn} g_{ielf}
C                       -t^{def}_{nml} L^{ad}_{mn} L_{ielf}
C
C     XIAJB contains g and YIAJB contains L
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYINT, ISYMC2, LENSQ, LWORK, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYRE1, ISYRES, ISYMBD, ISELMN, ISYAMN, ISYELI
      INTEGER LENGTH, KTMAT, KINT1, KINT2, KC2TEMP, KEND1, LWRK1
      INTEGER ISYMN, ISYMAM, ISYMA, ISYMM, ISYMMN, KOFF1, KOFF2, KOFF3
      INTEGER ISYML, ISYMEI, ISYMDL, ISYME, ISYMI, ISYMEL, NDL
      INTEGER NEI, NUMBEL, NUMBMN, NUMBA
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), SMAT(*), QMAT(*), TMAT(*), XIAJB(*)
      REAL YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION OMEGA1(*), SMAT(*), QMAT(*), TMAT(*), XIAJB(*)
      DOUBLE PRECISION YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('T3_ONEL2')
C
      ISYRE1 = MULD2H(ISYMIM,ISYINT)
      ISYRES = MULD2H(ISYRE1,ISYMC2)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISELMN = MULD2H(ISYMIM,ISYMBD)
C
      ISYAMN = MULD2H(ISYMB,ISYMC2)
      ISYELI = MULD2H(ISYMD,ISYINT)
C
      LENGTH = NCKIJ(ISELMN)
C
      KTMAT   = 1
      KINT1   = KTMAT   + NCKIJ(ISELMN)
      KINT2   = KINT1   + NCKI(ISYELI)
      KC2TEMP = KINT2   + NCKI(ISYELI)
      KEND1   = KC2TEMP + NMAIJA(ISYAMN)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in T3_ONEL2 (sort)')
      ENDIF
C
C-----------------------------
C     Sort C2
C-----------------------------
C
      DO ISYMN = 1, NSYM
         ISYMAM = MULD2H(ISYAMN,ISYMN)
         DO ISYMA = 1, NSYM
            ISYMM  = MULD2H(ISYMAM,ISYMA)
            ISYMMN = MULD2H(ISYMM,ISYMN)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = IT2SP(ISYAMN,ISYMB)
     *                  + NCKI(ISYAMN)*(B-1)
     *                  + ICKI(ISYMAM,ISYMN)
     *                  + NT1AM(ISYMAM)*(N-1)
     *                  + IT1AM(ISYMA,ISYMM)
     *                  + NVIR(ISYMA)*(M-1) 
     *                  + 1
C
                  KOFF2 = KC2TEMP - 1
     *                  + IMAIJA(ISYMMN,ISYMA)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M

C
                  CALL DCOPY(NVIR(ISYMA),C2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------------
C     Sort g integrals.
C---------------------------
C
      DO ISYML = 1, NSYM
         ISYMEI = MULD2H(ISYELI,ISYML)
         ISYMDL = MULD2H(ISYML,ISYMD)
         DO ISYME = 1, NSYM
            ISYMI  = MULD2H(ISYMEI,ISYME)
            ISYMEL = MULD2H(ISYME,ISYML)
C
            DO L = 1, NRHF(ISYML)
               NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
               DO E = 1, NVIR(ISYME)
                  DO I = 1, NRHF(ISYMI)
                     NEI = IT1AM(ISYME,ISYMI) + NVIR(ISYME)*(I-1) + E
C
                     KOFF1 = IT2AM(ISYMDL,ISYMEI) + INDEX(NDL,NEI)
                     KOFF2 = KINT1 - 1
     *                     + ICKI(ISYMEL,ISYMI)
     *                     + NT1AM(ISYMEL)*(I-1)
     *                     + IT1AM(ISYME,ISYML)
     *                     + NVIR(ISYME)*(L-1)
     *                     + E
                     KOFF3 = KINT2 - 1
     *                     + ICKI(ISYMEI,ISYML)
     *                     + NT1AM(ISYMEI)*(L-1)
     *                     + IT1AM(ISYME,ISYMI)
     *                     + NVIR(ISYME)*(I-1)
     *                     + E
C
                     WORK(KOFF2) = XIAJB(KOFF1)
                     WORK(KOFF3) = XIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,2))
     *           + QMAT(INDSQ(I,1))
C
      WORK(KTMAT-1+I) = SMAT(INDSQ(I,5))
     *                + QMAT(INDSQ(I,4))
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI  = MULD2H(ISYRES,ISYMA)
         ISYMEL = MULD2H(ISYELI,ISYMI)
         ISYMMN = MULD2H(ISYMEL,ISELMN)
C
         KOFF1 = ISAIKL(ISYMEL,ISYMMN) + 1
         KOFF2 = KINT1
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,TMAT(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ZERO,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMEL,ISYMMN)
         KOFF2 = KINT2
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,WORK(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ONE,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KC2TEMP
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KEND1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              ONE,WORK(KOFF1),NUMBMN,WORK(KOFF2),NUMBMN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
C
      ENDDO
C
C---------------------------
C     Sort L integrals.
C---------------------------
C
      DO ISYML = 1, NSYM
         ISYMEI = MULD2H(ISYELI,ISYML)
         ISYMDL = MULD2H(ISYML,ISYMD)
         DO ISYME = 1, NSYM
            ISYMI  = MULD2H(ISYMEI,ISYME)
            ISYMEL = MULD2H(ISYME,ISYML)
C
            DO L = 1, NRHF(ISYML)
               NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
               DO E = 1, NVIR(ISYME)
                  DO I = 1, NRHF(ISYMI)
                     NEI = IT1AM(ISYME,ISYMI) + NVIR(ISYME)*(I-1) + E
C
                     KOFF1 = IT2AM(ISYMDL,ISYMEI) + INDEX(NDL,NEI)
                     KOFF2 = KINT1 - 1
     *                     + ICKI(ISYMEL,ISYMI)
     *                     + NT1AM(ISYMEL)*(I-1)
     *                     + IT1AM(ISYME,ISYML)
     *                     + NVIR(ISYME)*(L-1)
     *                     + E
                     KOFF3 = KINT2 - 1
     *                     + ICKI(ISYMEI,ISYML)
     *                     + NT1AM(ISYMEI)*(L-1)
     *                     + IT1AM(ISYME,ISYMI)
     *                     + NVIR(ISYME)*(I-1)
     *                     + E
C
                     WORK(KOFF2) = YIAJB(KOFF1)
                     WORK(KOFF3) = YIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = - SMAT(INDSQ(I,1))
     *             - QMAT(INDSQ(I,2))
C
      WORK(KTMAT-1+I) = - SMAT(I)
     *                  - QMAT(INDSQ(I,3))
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI  = MULD2H(ISYRES,ISYMA)
         ISYMEL = MULD2H(ISYELI,ISYMI)
         ISYMMN = MULD2H(ISYMEL,ISELMN)
C
         KOFF1 = ISAIKL(ISYMEL,ISYMMN) + 1
         KOFF2 = KINT1
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,TMAT(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ZERO,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMEL,ISYMMN)
         KOFF2 = KINT2
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,WORK(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ONE,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KC2TEMP
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KEND1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              ONE,WORK(KOFF1),NUMBMN,WORK(KOFF2),NUMBMN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
C
      ENDDO
C
C----------------------------
C     End.
C----------------------------
C
      CALL QEXIT('T3_ONEL2')
C
      RETURN
      END
C  /* Deck t3_onel3 */
      SUBROUTINE T3_ONEL3(OMEGA1,SMAT,QMAT,TMAT,ISYMIM,XIAJB,ISYINT,
     *                    C2TP,ISYMC2,INDSQ,LENSQ,WORK,LWORK,
     *                    ISYMB,B,ISYMD,D)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate the term (t^{def}_{lmn} - t^{def}_{lnm}) L^{de}_{lm} L_{ianf}
C
C     Note : XIAJB is coming in as L and not g.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYINT, ISYMC2, LENSQ, LWORK, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYRE1, ISYRES, ISYMBD, ISFLMN, ISYAIN, ISYFLM, LENGTH
      INTEGER KTMAT, KC2TEMP, KINT, KEND1, LWRK1, ISYMM, ISYMFL
      INTEGER ISYMF, ISYML, ISYMLM, ISYMFM, KOFF1, KOFF2, KOFF3
      INTEGER ISYMI, ISYMAN, ISYMA, ISYMN, ISYMAI, ISYMBN, NBN
      INTEGER NAI, NUMFLM, NUMBAI
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), SMAT(*), QMAT(*), TMAT(*), XIAJB(*)
      REAL C2TP(*), WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION OMEGA1(*), SMAT(*), QMAT(*), TMAT(*), XIAJB(*)
      DOUBLE PRECISION C2TP(*), WORK(LWORK), ZERO, ONE
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('T3_ONEL3')
C
      ISYRE1 = MULD2H(ISYMIM,ISYMC2)
      ISYRES = MULD2H(ISYRE1,ISYINT)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISFLMN = MULD2H(ISYMIM,ISYMBD)
      ISYAIN = MULD2H(ISYMB,ISYINT)
      ISYFLM = MULD2H(ISYMC2,ISYMD)
C
      LENGTH = NCKIJ(ISFLMN)
C
      KTMAT   = 1
      KC2TEMP = KTMAT   + NCKIJ(ISFLMN)
      KINT    = KC2TEMP + NCKI(ISYFLM)
      KEND1   = KINT    + NCKI(ISYAIN)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in T3_ONEL3 (sort)')
      ENDIF
C
C-----------------------------
C     Sort C2
C-----------------------------
C
      DO ISYMM = 1, NSYM
         ISYMFL = MULD2H(ISYFLM,ISYMM)
         DO ISYMF = 1, NSYM
            ISYML = MULD2H(ISYMFL,ISYMF)
            ISYMLM = MULD2H(ISYMM,ISYML)
            ISYMFM = MULD2H(ISYMF,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO L = 1, NRHF(ISYML)
C
                  KOFF1 = IT2SP(ISYFLM,ISYMD)
     *                  + NCKI(ISYFLM)*(D-1)
     *                  + ICKI(ISYMFM,ISYML)
     *                  + NT1AM(ISYMFM)*(L-1)
     *                  + IT1AM(ISYMF,ISYMM)
     *                  + NVIR(ISYMF)*(M-1) 
     *                  + 1
C
                  KOFF2 = KC2TEMP - 1
     *                  + ICKI(ISYMFL,ISYMM)
     *                  + NT1AM(ISYMFL)*(M-1)
     *                  + IT1AM(ISYMF,ISYML)
     *                  + NVIR(ISYMF)*(L-1) 
     *                  + 1

C
                  CALL DCOPY(NVIR(ISYMF),C2TP(KOFF1),1,WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------------
C     Sort integrals.
C---------------------------
C
      DO ISYMI = 1, NSYM
         ISYMAN = MULD2H(ISYAIN,ISYMI)
         DO ISYMA = 1, NSYM
            ISYMN  = MULD2H(ISYMAN,ISYMA)
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMBN = MULD2H(ISYMB,ISYMN)
C
            DO N = 1, NRHF(ISYMN)
               NBN = IT1AM(ISYMB,ISYMN) + NVIR(ISYMB)*(N-1) + B
               DO A = 1, NVIR(ISYMA)
                  DO I = 1, NRHF(ISYMI)
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
                     KOFF1 = IT2AM(ISYMBN,ISYMAI) + INDEX(NBN,NAI)
                     KOFF2 = KINT - 1
     *                     + ICKI(ISYMAI,ISYMN)
     *                     + NT1AM(ISYMAI)*(N-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     WORK(KOFF2) = XIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(I)
     *           - SMAT(INDSQ(I,3))
     *           + QMAT(INDSQ(I,3))
     *           - QMAT(I)
C
      WORK(KTMAT-1+I) = SMAT(INDSQ(I,1))
     *                - SMAT(INDSQ(I,4))
     *                + QMAT(INDSQ(I,2))
     *                - QMAT(INDSQ(I,5))
      ENDDO
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      ISYMN = MULD2H(ISYRES,ISYAIN)
C
      KOFF1 = ISAIKJ(ISYFLM,ISYMN) + 1
      KOFF2 = IT2SP(ISYFLM,ISYMD)
     *      + NCKI(ISYFLM)*(D-1)
     *      + 1
      KOFF3 = KEND1
C
      CALL DZERO(WORK(KOFF3),NRHF(ISYMN))
C
      NUMFLM = MAX(1,NCKI(ISYFLM))
C
      CALL DGEMV('T',NCKI(ISYFLM),NRHF(ISYMN),ONE,
     *           TMAT(KOFF1),NUMFLM,C2TP(KOFF2),1,
     *           ONE,WORK(KOFF3),1)
C
      KOFF1 = KTMAT
     *      + ISAIKJ(ISYFLM,ISYMN)
      KOFF2 = KC2TEMP
      KOFF3 = KEND1
C
      NUMFLM = MAX(1,NCKI(ISYFLM))
C
      CALL DGEMV('T',NCKI(ISYFLM),NRHF(ISYMN),ONE,
     *           WORK(KOFF1),NUMFLM,WORK(KOFF2),1,
     *           ONE,WORK(KOFF3),1)
C
      KOFF1 = KINT
     *      + ICKI(ISYRES,ISYMN)
      KOFF2 = KEND1
      KOFF3 = 1
C
      NUMBAI = MAX(1,NT1AM(ISYRES))
C
      CALL DGEMV('N',NT1AM(ISYRES),NRHF(ISYMN),ONE,
     *           WORK(KOFF1),NUMBAI,WORK(KOFF2),1,
     *           ONE,OMEGA1(KOFF3),1)
C
C----------------------------
C     End.
C----------------------------
C
      CALL QEXIT('T3_ONEL3')
C
      RETURN
      END
C  /* Deck cc3_lhtr_l3 */
      SUBROUTINE CC3_L3_LHTR(ECURR,L1AM,ISYML1,L2TP,ISYML2,T2TP,ISYMT2,
     *                       OMEGA1,OMEGA2,ISYRES,XINT4O,XOVVO,XOOVV,
     *                       XLAMDP,XLAMDH,WORK,LWORK,
     *                       LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUTOC,FNTOC,
     *                       LU3VI,FN3VI,LU4V,FN4V,LUDKBC3,FNDKBC3,
     *                       LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate the L3 dependent terms that contribute to
C     the left hand side eigenvalue equation.
C
C     ISYML1 is symmetry of L1AM
C     ISYML2 is symmetry of L2TP
C     Isyres is symmetry of Omega{1,2}
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <iratdef.h>
#include <ccinftap.h>
#include <second.h>
C
      CHARACTER*6 FNGEI,FNFEI
      CHARACTER*5 FNN1
      PARAMETER( FNGEI = 'N1_GEI' , FNFEI = 'N1_FEI' , FNN1 = 'N1MAT' )
      INTEGER LUGEI,LUFEI,LUN1
C
      INTEGER ISYML1, ISYML2, ISYMT2, ISYRES, LWORK
      INTEGER ISYMTR, ISINT1, ISINT2, ISYMIM, KFOCKD, KCMO, KFCKBA
      INTEGER KEND0, LWRK0, KEND1, LWRK1, KEND2, LWRK2, ISYMC, ISYMK
      INTEGER KOFF1, KOFF2, KXIAJB, KINTOC, LUFCK
      INTEGER LENGTH, ISYOPE, IOPTTCME, IOFF, ISYMD, ISAIJ1, ISYCKB
      INTEGER ISCKB2
      INTEGER KRMAT1, KEND3, LWRK3, KEND4, LWRK4
      INTEGER ISYMB, ISYALJ, ISAIJ2, ISYMBD, ISCKIJ
      INTEGER KDIAG, LENSQ, KINDSQ,  KTMAT
      INTEGER ISCKB1, KTRVI, KTRVI1, KTROC, KTROC1, KVVVV
      INTEGER LUDKBC4, LUCKJD
      INTEGER LUDKBC, LUTOC, LU3VI, LU4V, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER LU3VI2,LU3FOP,LU3FOP2
C
      INTEGER KRBJIA,KW3BXOG1,KW3BXOL1,KW3BXVDG1,KW3BXVDG2
      INTEGER KW3BXVDL1,KW3BXVDL2,ISYALJBL1,ISYALJDL1,KWMAT,KINDEXBL1
      INTEGER KINDEXDL1
C
      INTEGER ISYMN1,ISYMN2,KN2MAT,KINDSQN,LENSQN
C
      INTEGER ISGEI,ISFEI,KGEI,KFEI,IADR
      INTEGER IOPT
c
      integer kx3am
c

c
      integer isymi,isymfge,isymanm
C
#if defined (SYS_CRAY)
      REAL L1AM(*), L2TP(*), T2TP(*), OMEGA1(*), OMEGA2(*), XINT4O(*)
      REAL XOVVO(*), XOOVV(*), XLAMDP(*), XLAMDH(*), WORK(LWORK)
      REAL TITRAN, TISORT,  TICONT, TIOME1
      REAL DTIME, XL2TP, XIAJB, XINT, XTROC0, XTRVI
      REAL RHO1N, RHO2N, XDIA, DDOT, HALF, ONE, ECURR
#else
      DOUBLE PRECISION L1AM(*), L2TP(*), T2TP(*)
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), XINT4O(*)
      DOUBLE PRECISION XOVVO(*),XOOVV(*),XLAMDP(*),XLAMDH(*),WORK(LWORK)
      DOUBLE PRECISION TITRAN, TISORT,  TICONT, TIOME1
      DOUBLE PRECISION DTIME, XL2TP, XIAJB, XINT, XTROC0, XTRVI
      DOUBLE PRECISION RHO1N, RHO2N, XDIA, DDOT, HALF, ONE, ECURR
#endif
C
      CHARACTER*(*) FNCKJD, FNDKBC, FNTOC, FN3VI, FN4V
      CHARACTER*(*) FNDKBC3, FN3FOPX, FN3FOP2X
C
      CHARACTER*5 FN3FOP
      CHARACTER*6 FN3FOP2
      CHARACTER*8 FN3VI2
C
      PARAMETER ( FN3FOP = 'PTFOP', FN3FOP2 = 'PTFOP2', 
     *            FN3VI2 = 'CC3_VI12'                 )
C
      CHARACTER*11 FNDKBC4
      CHARACTER*1 CDUMMY
C
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('CC3_L3_LHTR')
C

C
C----------------------------------------------------
C     Initialise character strings and open files
C----------------------------------------------------
C
      CDUMMY = ' '
      LUDKBC4 = -1
      FNDKBC4 = 'CC3_L3_TMP1'
C
      LU3FOP  = -1
      LU3FOP2 = -1
      LU3VI2  = -1
C
      CALL WOPEN2(LU3FOP,FN3FOP,64,0)
      CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
C
      CALL WOPEN2(LUDKBC4,FNDKBC4,64,0)
C
      IF (.NOT.LVVVV) THEN
         !Open files for N1MAT intermediates
         LUGEI = -1
         LUFEI = -1
         LUN1  = -1
         CALL WOPEN2(LUGEI,FNGEI,64,0)
         CALL WOPEN2(LUFEI,FNFEI,64,0)
         CALL WOPEN2(LUN1,FNN1,64,0)
      END IF
C
C-------------------------------------------------------------
C     Set symmetry flags.
C
C     omega = int1*T2*int2
C     isymres is symmetry of result(omega)
C     isint1 is symmetry of integrals in contraction.(int1)
C     isint2 is symmetry of integrals in the triples equation.(int2)
C     isymim is symmetry of S and Q intermediates.(t2*int2)
C      (sym is for all index of S and Q (cbd,klj)
C       thus cklj=b*d*isymim)
C-------------------------------------------------------------
C
      IPRCC = IPRINT
      ISYMTR = ISYML1
      ISINT1 = ISYMOP
      ISINT2 = ISYMOP
      ISYMIM = MULD2H(ISYMTR,ISYMOP)
C
C--------------------
C     Time variables.
C--------------------
C
      TITRAN = 0.0D0
      TISORT = 0.0D0
      TICONT = 0.0D0
      TIOME1 = 0.0D0
C
C-----------------------------------------------------------
C     Calculate 2*C-E and store 
C     FNDKBC3, FN3FOPX and FN3FOP2X for f.o.p. later. 
C-----------------------------------------------------------
C
      CALL CC3_TCME(XLAMDP,ISINT1,WORK,LWORK,LU3VI,FN3VI,
     *              LUDKBC,FNDKBC,LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *              LUDKBC3,FNDKBC3,LUDKBC4,FNDKBC4,1)
C
C--------------------------------------
C     Reorder the l2-amplitudes i L2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYML2)) THEN
        CALL QUIT('Not enough memory to construct L2TP in CC3_LHTR_L3')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYML2),L2TP,1,WORK,1)
      CALL CC3_T2TP(L2TP,WORK,ISYML2)
C
      IF (IPRINT .GT. 55) THEN
         XL2TP = DDOT(NT2SQ(ISYML2),L2TP,1,L2TP,1)
         WRITE(LUPRI,*) 'Norm of L2TP ',XL2TP
      ENDIF
C
C--------------------------------------
C     Reorder the T2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
        CALL QUIT('Not enough memory to construct T2TP in CC3_LHTR_L3')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMT2),T2TP,1,WORK,1)
      CALL CC3_T2TP(T2TP,WORK,ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XL2TP = DDOT(NT2SQ(ISYMT2),T2TP,1,T2TP,1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XL2TP
      ENDIF
C
C---------------------------------------------------------
C     Read canonical orbital energies and MO coefficients.
C---------------------------------------------------------
C
      IF (.NOT.LVVVV) THEN
         !Symmetries for N1 and N2 intermediates
         ISYMN1 = MULD2H(ISYMIM,ISYMT2)
         ISYMN2 = MULD2H(ISYMIM,ISYMT2)
      END IF
C
      IF (LVVVV) THEN
         KRBJIA = 1
      ELSE 
         KN2MAT = 1
         KRBJIA = KN2MAT + NCKIJ(ISYMN2)
      END IF
      KFOCKD = KRBJIA + NT2SQ(ISYRES)
      KCMO   = KFOCKD + NORBTS
      KFCKBA = KCMO   + NLAMDS
      KEND0  = KFCKBA + N2BST(ISYMOP)
      LWRK0  = LWORK  - KEND0
    
      IF (.NOT.LVVVV) THEN
         KINDSQN = KEND0
         KEND0  = KINDSQN + (6*NCKIJ(ISYMN2) - 1)/IRAT + 1
         LWRK0  = LWORK - KEND0 
      END IF
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_LHTR_L3')
      END IF
C
      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
C
      IF (.NOT.LVVVV) THEN
         CALL DZERO(WORK(KN2MAT),NCKIJ(ISYMN2))
C
         !index array for N2
         LENSQN = NCKIJ(ISYMN2)
         CALL CC3_INDSQ(WORK(KINDSQN),LENSQN,ISYMN2)
      END IF
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND0),LWRK0)
C
C---------------------------------------------
C     Delete frozen orbitals in Fock diagonal.
C---------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
C-----------------------------------------------------
C     Construct the transformed Fock matrix
C-----------------------------------------------------
C
      LUFCK = -1
C     This AO Fock matrix is constructed from the T1 transformed density
      CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
C     This AO Fock matrix is constructed from the CMO transformed density
C      CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
C     *            IDUMMY,.FALSE.)
      REWIND(LUFCK)
      READ(LUFCK)(WORK(KFCKBA + I-1),I = 1,N2BST(ISYMOP))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT. 140) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
      CALL CC_FCKMO(WORK(KFCKBA),XLAMDP,XLAMDH,
     *              WORK(KEND0),LWRK0,1,1,1)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC3_L3_LHTR: Triples Fock MO matrix' )
         CALL CC_PRFCKMO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
C     Sort the fock matrix
C
C
      CALL DCOPY(N2BST(ISINT1),WORK(KFCKBA),1,WORK(KEND0),1)
C
      DO ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISINT1)
C
         DO K = 1,NRHF(ISYMK)
C
            DO C = 1,NVIR(ISYMC)
C
               KOFF1 = KEND0 + IFCVIR(ISYMK,ISYMC) + 
     *                 NORB(ISYMK)*(C - 1) + K - 1
               KOFF2 = KFCKBA + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C - 1
C
               WORK(KOFF2) = WORK(KOFF1)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('In CC3_L3_LHTR: Triples Fock MO matrix (sort)')
         CALL CC_PRFCKMO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
C----------------------------------------
C     If we want to sum the T3 amplitudes
C----------------------------------------
C
      if (.false.) then
         kx3am  = kend0
         kend0 = kx3am + nrhft*nrhft*nrhft*nvirt*nvirt*nvirt
         call dzero(work(kx3am),nrhft*nrhft*nrhft*nvirt*nvirt*nvirt)
         lwrk0 = lwork - kend0
         if (lwrk0 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend0
            call quit('Insufficient space (T3) in CC3_LHTR')
         END IF
      endif
C
C      write(lupri,*) 'WBMAT after dzero'
C      call print_pt3(work(kx3am),ISYML1,4)

C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC  = KEND0
      KTROC1 = KTROC  + NTRAOC(ISINT2)
      KXIAJB = KTROC1 + NTRAOC(ISINT2)
      KEND1  = KXIAJB + NT2AM(ISYMOP)
      LWRK1  = LWORK  - KEND1
C
      KW3BXOG1   = KEND1
      KW3BXOL1   = KW3BXOG1 + NTRAOC(1)
      KEND1   = KW3BXOL1    + NTRAOC(1)
      LWRK1      = LWORK    - KEND1
C
      KINTOC = KEND1
      KEND2  = KINTOC + MAX(NTOTOC(ISYMOP),NTOTOC(ISINT2))
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_LHTR_L3')
      END IF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      LENGTH = IRAT*NT2AM(ISYMOP)
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
C
      ISYOPE = ISYMOP
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYOPE,IOPTTCME)
C
      IF ( IPRINT .GT. 55) THEN
         XIAJB = DDOT(NT2AM(ISYMOP),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XIAJB
      ENDIF
C
C------------------------------------------------------------
C     Read in integrals used in contractions and transform.
C------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC),XLAMDP,
     *                    WORK(KEND2),LWRK2,ISINT2)
C
      CALL CCFOP_SORT(WORK(KTROC),WORK(KTROC1),ISINT2,1)
C
      CALL CC3_LSORT1(WORK(KTROC),ISINT2,WORK(KEND2),LWRK2,5)
C
C     -----------------------------------------------------------------
C     Occupied integrals needed for t3bar0 constructed in terms of WMAT
C     -----------------------------------------------------------------
C
      !kend1 can be used again since KINTOC is not needed any more
      CALL INTOCC_T3BARX(.TRUE.,
     *                   LUTOC,FNTOC,ISYMOP,XLAMDH,1,
     *                   ISINT1,
     *                   DUMMY,IDUMMY,IDUMMY,
     *                   WORK(KW3BXOG1),
     *                   WORK(KW3BXOL1),DUMMY,
     *                   DUMMY,
     *                   WORK(KEND1),LWRK1)

C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM
C
         ISAIJ1 = MULD2H(ISYMD,ISYRES)
         ISYCKB = MULD2H(ISYMD,ISYMOP)
         ISCKB1 = MULD2H(ISINT1,ISYMD)
         ISCKB2 = MULD2H(ISINT2,ISYMD)
C
         IF (.NOT.LVVVV) THEN
            !Symmetry of arrays needed to construct N1MAT
            ISGEI  = MULD2H(ISYMN1,ISYMD)
            ISFEI  = MULD2H(ISYMN1,ISYMD)
         END IF
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISAIJ1 :',ISAIJ1
            WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYCKB :',ISYCKB
            WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISCKB2 :',ISCKB2
C
         ENDIF
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KTRVI  = KEND1
         KTRVI1 = KTRVI  + NCKATR(ISCKB1)
         KRMAT1 = KTRVI1 + NCKATR(ISCKB1)
         KEND2  = KRMAT1 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         IF (LVVVV) THEN
            KVVVV   = KEND2
            KEND3   = KVVVV   + NMAABC(ISCKB2)
            LWRK3   = LWORK  - KEND3
         ELSE
           KEND3    = KEND2
           LWRK3    = LWORK  - KEND3
         END IF
C
         KW3BXVDG1  = KEND3
         KW3BXVDG2  = KW3BXVDG1  + NCKATR(ISCKB2)
         KW3BXVDL1  = KW3BXVDG2  + NCKATR(ISCKB2)
         KW3BXVDL2  = KW3BXVDL1  + NCKATR(ISCKB2)
         KEND3     = KW3BXVDL2  + NCKATR(ISCKB2)
         LWRK3     = LWORK     - KEND3
C
         IF (.NOT.LVVVV) THEN
            !Arrays needed to construct N1MAT
            KGEI  = KEND3
            KFEI  = KGEI  + NCKATR(ISGEI)
            KEND3 = KFEI  + NCKATR(ISFEI)
            LWRK3 = LWORK - KEND3
         END IF
C
         KEND4  = KEND3
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_LHTR_L3')
         END IF
C
C---------------------
C        Sum over D
C---------------------
C
         DO D = 1,NVIR(ISYMD)
C
C------------------------------------
C           Initialize the R1 matrix.
C------------------------------------
C
            CALL DZERO(WORK(KRMAT1),NCKI(ISAIJ1))
            IF (.NOT.LVVVV) THEN
               CALL DZERO(WORK(KGEI),NCKATR(ISGEI))
               CALL DZERO(WORK(KFEI),NCKATR(ISFEI))
            END IF
C
            IF (LVVVV) THEN
C
C              ---------------------------------
C              Read in g_{vvvv} for a given D
C              ---------------------------------
C
               IF (NMAABC(ISCKB2) .GT. 0) THEN
                  IOFF = I3VVIR(ISCKB2,ISYMD)
     *                 + NMAABC(ISCKB2)*(D-1)
     *                 + 1
                  CALL GETWA2(LU4V,FN4V,WORK(KVVVV),IOFF,NMAABC(ISCKB2))
               ENDIF
C
            END IF
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction.
C------------------------------------------------------------
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
               IOFF = ICKBD(ISCKB1,ISYMD) + NCKATR(ISCKB1)*(D - 1) + 1
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI),IOFF,
     &                     NCKATR(ISCKB1))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKB1)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_L3 (TRVI)')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI),WORK(KEND4),ISYMD,D,ISINT1)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
               IOFF = ICKBD(ISCKB1,ISYMD) + NCKATR(ISCKB1)*(D - 1) + 1
               CALL GETWA2(LUDKBC4,FNDKBC4,WORK(KTRVI1),IOFF,
     &                     NCKATR(ISCKB1))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKB1)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_L3 (TRVI1)')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI1),WORK(KEND4),ISYMD,D,ISINT1)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
C           ------------------------------------------------------
C           Integrals needed to construct t3bar0 in terms of WMAT
C           ------------------------------------------------------
C
            CALL INTVIR_T3BARX_D(.TRUE.,
     *                           ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                           LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                           DUMMY,WORK(KW3BXVDG1),
     *                           DUMMY,WORK(KW3BXVDG2),
     *                           DUMMY,WORK(KW3BXVDL1),
     *                           DUMMY,WORK(KW3BXVDL2),
     *                           DUMMY,IDUMMY,XLAMDP,
     *                           1,ISYMD,D,WORK(KEND4),LWRK4)
C
C---------------------
C           Calculate.
C---------------------
C
            DO ISYMB = 1,NSYM
C
               ISYALJ  = MULD2H(ISYMB,ISYML2)
               ISAIJ2  = MULD2H(ISYMB,ISYRES)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYMIM)
C
               ISYALJBL1  = MULD2H(ISYMB,ISYML1)
               ISYALJDL1 = MULD2H(ISYMD,ISYML1)
C
               IF (IPRINT .GT. 55) THEN
C
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISCKIJ:',ISCKIJ
C
               ENDIF
C
C              Can use kend3 since we do not need the integrals anymore.
               KDIAG   = KEND3
               KINDSQ  = KDIAG   + NCKIJ(ISCKIJ)
               KTMAT   = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KEND4   = KTMAT   + NCKIJ(ISCKIJ)
               LWRK4   = LWORK   - KEND4
C
               KWMAT   = KEND4
               KEND4   = KWMAT   + NCKIJ(ISCKIJ)
               LWRK4   = LWORK   - KEND4
C
               KINDEXBL1   = KEND4
               KINDEXDL1  = KINDEXBL1 + (NCKI(ISYALJBL1)-1)/IRAT + 1
               KEND4      = KINDEXDL1  + (NCKI(ISYALJDL1)-1)/IRAT + 1
               LWRK4   = LWORK   - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_LHTR_L3 (inner)')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
C
               IF (IPRINT .GT. 55) THEN
                  XDIA  = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                    WORK(KDIAG),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XDIA
               ENDIF

C
C-------------------------------------
C              Construct index arrays.
C-------------------------------------
C
               LENSQ = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
C
               !Get indeces for t3bar0 constructed in terms of WMAT
               CALL CC3_INDEX(WORK(KINDEXBL1),ISYALJBL1)
               CALL CC3_INDEX(WORK(KINDEXDL1),ISYALJDL1)
C
               DO B = 1,NVIR(ISYMB)
C
C                 --------------------------------------------------
C                 Calculate t3bar0 in terms of W intermediate
C                 --------------------------------------------------
C
                  DTIME = SECOND()
 
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISCKIJ))

                  !<L2Y|[H^,tau3]|HF>
                  CALL WBARBD_TMAT(L2TP,ISYML2,WORK(KWMAT),WORK(KTMAT),
     *                             ISCKIJ,WORK(KFCKBA),ISYMOP,
     *                             WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                             WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                             WORK(KW3BXOL1),WORK(KW3BXOG1),ISINT2,
     *                             WORK(KEND4),LWRK4,WORK(KINDEXBL1),
     *                             WORK(KINDEXDL1),WORK(KINDSQ),LENSQ,
     *                             ISYMB,B,ISYMD,D)
 
                  !<L1Y|[H^,tau3]|HF>
                  CALL WBARBD_L1(L1AM,ISYML1,WORK(KTMAT),WORK(KXIAJB),
     *                           ISYMOP,WORK(KWMAT),WORK(KEND4),LWRK4,
     *                           WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                  !Divide by the energy difference and
                  !remove the forbidden elements
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,ECURR,ISCKIJ,WORK(KWMAT),
     *                         WORK(KDIAG),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMAT),ISYMIM,ISYMB,B,ISYMD,D)

*       call sum_pt3(work(KWMAT),isymb,b,isymd,d,
*    *             1,work(kx3am),4)
C
C-----------------------------------------------------------------------
C                 Calculate the contributions to omega2
C-----------------------------------------------------------------------
C
                  CALL CC3_W3_CY2V(OMEGA2,ISYRES,WORK(KRBJIA),
     *                             WORK(KWMAT),ISCKIJ,
     *                             WORK(KTMAT),WORK(KTRVI),WORK(KTRVI1),
     *                             ISINT1,WORK(KEND4),LWRK4,
     *                             WORK(KINDSQ),LENSQ,
     *                             ISYMB,B,ISYMD,D,.TRUE.)
C
                  CALL CC3_W3_CY2O(OMEGA2,ISYRES,WORK(KWMAT),ISCKIJ,
     *                             WORK(KTMAT),WORK(KTROC),WORK(KTROC1),
     *                             ISINT1,WORK(KEND4),LWRK4,
     *                             WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D,
     *                             .TRUE.)
C
C------------------------------------------------------------------------
C                 Calculate the L3 contribution to omega1
C------------------------------------------------------------------------
C
                  IF (LVVVV) THEN
                     CALL CC3_W3_OMEGA1(OMEGA1,ISYRES,WORK(KWMAT),
     *                                  WORK(KTMAT),ISYMIM,
     *                                  XINT4O,XOVVO,
     *                                  XOOVV,WORK(KVVVV),1,
     *                                  T2TP,ISYMT2,
     *                                  WORK(KEND4),LWRK4,
     *                                  LENSQ,WORK(KINDSQ),
     *                                  ISYMB,B,ISYMD,D,.TRUE.)
                  ELSE
                     CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KWMAT),1)

                     !Construct N1 and N2 intermediates
                     CALL WT2_N1N2(WORK(KWMAT),ISYMIM,
     *                       T2TP,ISYMT2,
     *                       WORK(KGEI),WORK(KFEI),
     *                       ISYMN1,
     *                       WORK(KN2MAT),ISYMN2,
     *                       B,ISYMB,D,ISYMD,
     *                       WORK(KINDSQ),LENSQ,
     *                       WORK(KINDSQN),LENSQN,
     *                       WORK(KEND4),LWRK4,
     *                       .TRUE.)

C
                  END IF
C
                  IF (IPRINT .GT. 55) THEN
                     RHO2N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm (Rho1) after CC3_L3_OMEGA1',
     *                               RHO2N
                  ENDIF
C
C              ----------
C              End B loop
C              ----------
C
               ENDDO   ! B
            ENDDO      ! ISYMB

            IF (.NOT.LVVVV) THEN
C
C              ----------------------------------------------------------
C              Put KGEI(ge,i)^F and KFEI(fe,i)^G (which are intermediates
C              for N1MAT(fge,i) ) to files (for fixed F=D and G=D).
C              ----------------------------------------------------------

               !Put KGEI to file as (gei,F)   (fixed F corresponds to D)
               IADR = ICKBD(ISGEI,ISYMD) + NCKATR(ISGEI)*(D-1) + 1
               CALL PUTWA2(LUGEI,FNGEI,WORK(KGEI),IADR,NCKATR(ISGEI))
C
               !Put KFEI to file as (fei,G)   (fixed G corresponds to D)
               IADR = ICKBD(ISFEI,ISYMD) + NCKATR(ISFEI)*(D-1) + 1
               CALL PUTWA2(LUFEI,FNFEI,WORK(KFEI),IADR,NCKATR(ISFEI))
C
            END IF
C
C          -----------
C           End D loop
C          -----------
C
         ENDDO       ! D
      ENDDO          ! ISYMD

C
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Vccupied  cont ) 
C     in XI2EFF 
C------------------------------------------------------
C
      CALL CC3_RBJIA(OMEGA2,ISYRES,WORK(KRBJIA))
C
      IF (IPRINT .GT. 55) THEN
         RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
         WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC-2',RHO2N
      ENDIF
C
      IF (IPRINT .GT. 220) THEN
         CALL AROUND('After CC3_RACC-2: ')
         CALL CC_PRP(DUMMY,OMEGA2,ISYRES,0,1)
      ENDIF
C
C
      IF (.NOT.LVVVV) THEN
C
         !Read (gei,F) and (fei,G) intermediates from files
         !add them and put the result to a file as (fge,I)
         CALL N1_RESORT(ISYMN1,LUN1,FNN1,LUGEI,FNGEI,LUFEI,FNFEI,
     *                  WORK(KEND0),LWRK0,.FALSE.)
C
         !Calculate <T3|[[H,T2],tau_ai]|HF> except VVVV contribution
         CALL N1N2_G(LUN1,FNN1,
     *                     ISYMN1,
     *                     WORK(KN2MAT),ISYMN2,
     *                     XOVVO,XOOVV,XINT4O,1,
     *                     OMEGA1,ISYRES,
     *                     WORK(KINDSQN),LENSQN,
     *                     WORK(KEND0),LWRK0)
C
         !Calculate VVVV contribution to <T3|[[H,T2],tau_ai]|HF>
         IOPT = 0 !normal Lambda matrices used in backtransformation
         CALL  N1_GV4(IOPT,
     *                LUN1,FNN1,
     *                ISYMN1,
     *                XLAMDP,1, 
     *                XLAMDP,1, 
     *                XLAMDH,1, 
     *                XLAMDH,1, 
     *                OMEGA1,ISYRES,
     *                WORK(KEND0),LWRK0)
C
      END IF
C
C-------------------------------
C     Close and delete files
C-------------------------------
C
*     write(lupri,*)'Omega1 (final) isyres', isyres
*     call PRINT_MATAI(OMEGA1,ISYRES)

*      write(lupri,*) 't3barx  in CC3_LHTR'
*      call print_pt3(work(kx3am),1,4)



      CALL WCLOSE2(LUDKBC4,FNDKBC4,'DELETE')
C
      CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
      CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
C
      IF (.NOT.LVVVV) THEN
         !Close files for N1MAT intermediates
         CALL WCLOSE2(LUGEI,FNGEI,'DELETE')
         CALL WCLOSE2(LUFEI,FNFEI,'DELETE')
         CALL WCLOSE2(LUN1,FNN1,'DELETE')
      END IF
C
C-------------------
C     Print timings.
C-------------------
C
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,1) 'CC3_TRAN  : ',TITRAN
         WRITE(LUPRI,1) 'CC3_SORT  : ',TISORT
         WRITE(LUPRI,1) 'CC3_OME1  : ',TIOME1
         WRITE(LUPRI,*)
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_L3_LHTR')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck cc3_tcme */
      SUBROUTINE CC3_TCME(XLAMDP,ISYINT,WORK,LWORK,LUL3VI,FNL3VI,
     *                    LULDKBC,FNLDKBC,LU3FOPX,FN3FOPX,
     *                    LU3FOP2X,FN3FOP2X,LULDKBC3,FNLDKBC3,
     *                    LULDKBC4,FNLDKBC4,IOPT)
C
C     K. Hald, Spring 2002.
C     Purpose : Calculate 2*C-E of integrals with
C               3 virtual indices and 1 occupied.
C
      IMPLICIT NONE
C
      INTEGER ISYINT, LWORK, ISYMD, ISYMB, KTRVI, KINTVI, KEND1, LWRK1
      INTEGER IOFF, KOFF1, ISYMCK, ISYCKB, ISYCKD, KTRVI1, IOPT
      INTEGER LUL3VI, LULDKBC, LU3FOPX, LU3FOP2X, LULDKBC3, LULDKBC4
      INTEGER LUTMP1, LUTMP2, LUTMP3, LUTMP4
C
#if defined (SYS_CRAY)
      REAL XLAMDP(*), WORK(LWORK), ONE, TWO, HALF
#else
      DOUBLE PRECISION XLAMDP(*), WORK(LWORK), ONE, TWO, HALF
#endif
C
#include <priunit.h>
#include <ccinftap.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      CHARACTER*(*) FNL3VI, FNLDKBC, FN3FOPX, FN3FOP2X
      CHARACTER*(*) FNLDKBC3, FNLDKBC4
      CHARACTER*10 FNTMP1, FNTMP2, FNTMP3, FNTMP4
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CC3_TCME')
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2)) THEN
         CALL QUIT('Wrong IOPT in CC3_TCME')
      ENDIF
C
C-------------------------------
C     For IOPT=1 calculate
C-------------------------------
C
      IF (IOPT .EQ. 1) THEN
C
C------------------------------------
C     Open temporary files
C------------------------------------
C
         LUTMP1 = -1
         LUTMP2 = -1
         LUTMP3 = -1
         LUTMP4 = -1
         FNTMP1 = 'CC3_TCME_1'
         FNTMP2 = 'CC3_TCME_2'
         FNTMP3 = 'CC3_TCME_3'
         FNTMP4 = 'CC3_TCME_4'
C
         CALL WOPEN2(LUTMP1,FNTMP1,64,0)
         CALL WOPEN2(LUTMP2,FNTMP2,64,0)
         CALL WOPEN2(LUTMP3,FNTMP3,64,0)
         CALL WOPEN2(LUTMP4,FNTMP4,64,0)
C
C--------------------------------------------------------
C     Transform the integrals and sort (ckbd) to (ckdb)
C--------------------------------------------------------
C
         DO ISYMD = 1, NSYM
C
            ISYCKB = MULD2H(ISYINT,ISYMD)
C
            KTRVI  = 1
            KINTVI = KTRVI  + NCKATR(ISYCKB)
            KEND1  = KINTVI + NCKA(ISYCKB)
            LWRK1  = LWORK - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Out of memory in CC3_TCME')
            ENDIF
C
            DO D = 1, NVIR(ISYMD)
C
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LUL3VI,FNL3VI,WORK(KINTVI),IOFF,
     *                        NCKA(ISYCKB))
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),XLAMDP,
     *                             ISYMD,D,ISYINT,WORK(KEND1),LWRK1)
C
                  DO ISYMB = 1, NSYM
C
                     ISYMCK = MULD2H(ISYCKB,ISYMB)
                     ISYCKD = MULD2H(ISYINT,ISYMB)
C
                     DO B = 1, NVIR(ISYMB)
C
                        KOFF1 = KTRVI
     *                        + ICKATR(ISYMCK,ISYMB)
     *                        + NT1AM(ISYMCK)*(B-1)
C
                        IOFF  = ICKBD(ISYCKD,ISYMB) 
     *                        + NCKATR(ISYCKD)*(B - 1) 
     *                        + ICKATR(ISYMCK,ISYMD)
     *                        + NT1AM(ISYMCK)*(D-1)
     *                        + 1
C
                        CALL PUTWA2(LUTMP1,FNTMP1,WORK(KOFF1),IOFF,
     *                              NT1AM(ISYMCK))
C
                     ENDDO
                  ENDDO
C
               ENDIF
C
            ENDDO
         ENDDO
C
C
C-----------------------
C     Calculate 2*C-E
C-----------------------
C
         DO ISYMD = 1, NSYM
            ISYCKB = MULD2H(ISYINT,ISYMD)
C
            KTRVI  = 1
            KTRVI1 = KTRVI  + NCKATR(ISYCKB)
            KEND1  = KTRVI1 + NCKATR(ISYCKB)
            LWRK1  = LWORK  - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Out of memory in CC3_TCME')
            ENDIF
C
            DO D = 1, NVIR(ISYMD)
C
               IOFF  = ICKBD(ISYCKB,ISYMD) 
     *               + NCKATR(ISYCKB)*(D - 1) 
     *               + 1
C
               IF (NCKATR(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LUTMP1,FNTMP1,WORK(KTRVI),IOFF,
     *                        NCKATR(ISYCKB))
C
                  CALL CCSDT_SRTVIR(WORK(KTRVI),WORK(KTRVI1),
     *                              WORK(KEND1),LWRK1,ISYMD,ISYINT)
C
                  CALL DSCAL(NCKATR(ISYCKB),TWO,WORK(KTRVI),1)
C
                  CALL DAXPY(NCKATR(ISYCKB),-ONE,WORK(KTRVI1),1,
     *                       WORK(KTRVI),1)
C
                  CALL PUTWA2(LUTMP2,FNTMP2,WORK(KTRVI),IOFF,
     *                        NCKATR(ISYCKB))
                  CALL PUTWA2(LUTMP3,FNTMP3,WORK(KTRVI1),IOFF,
     *                        NCKATR(ISYCKB))
C
                  CALL DAXPY(NCKATR(ISYCKB),ONE,WORK(KTRVI1),1,
     *                       WORK(KTRVI),1)
                  CALL DSCAL(NCKATR(ISYCKB),HALF,WORK(KTRVI),1)
                  CALL DSCAL(NCKATR(ISYCKB),TWO,WORK(KTRVI1),1)
                  CALL DAXPY(NCKATR(ISYCKB),-ONE,WORK(KTRVI),1,
     *                       WORK(KTRVI1),1)
C
                  CALL PUTWA2(LUTMP4,FNTMP4,WORK(KTRVI1),IOFF,
     *                        NCKATR(ISYCKB))
C
               ENDIF
C
            ENDDO
         ENDDO
C
C-----------------------
C     Final sort for L
C-----------------------
C
         DO ISYMD = 1, NSYM
            ISYCKB = MULD2H(ISYINT,ISYMD)
C
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               DO D = 1, NVIR(ISYMD)
C
                  DO ISYMB = 1, NSYM
C
                     ISYMCK = MULD2H(ISYCKB,ISYMB)
                     ISYCKD = MULD2H(ISYMCK,ISYMD)
C
                     KTRVI = 1
                     KEND1 = KTRVI + NT1AM(ISYMCK)
                     LWRK1 = LWORK - KEND1
C
                     DO B = 1, NVIR(ISYMB)
C
                        IOFF  = ICKBD(ISYCKB,ISYMD) 
     *                        + NCKATR(ISYCKB)*(D - 1) 
     *                        + ICKATR(ISYMCK,ISYMB)
     *                        + NT1AM(ISYMCK)*(B-1)
     *                        + 1
C
                        CALL GETWA2(LUTMP2,FNTMP2,WORK(KTRVI),IOFF,
     *                              NT1AM(ISYMCK))
C
                        IOFF  = ICKBD(ISYCKD,ISYMB) 
     *                        + NCKATR(ISYCKD)*(B - 1) 
     *                        + ICKATR(ISYMCK,ISYMD)
     *                        + NT1AM(ISYMCK)*(D-1)
     *                        + 1
C
                        CALL PUTWA2(LU3FOPX,FN3FOPX,WORK(KTRVI),
     *                              IOFF,NT1AM(ISYMCK))
C
                     ENDDO
                  ENDDO
C
               ENDDO
            ENDIF
         ENDDO
C
C----------------------------
C     Final sort for new G
C----------------------------
C
         DO ISYMD = 1, NSYM
            ISYCKB = MULD2H(ISYINT,ISYMD)
C
            DO D = 1, NVIR(ISYMD)
C
               DO ISYMB = 1, NSYM
C
                  ISYMCK = MULD2H(ISYCKB,ISYMB)
                  ISYCKD = MULD2H(ISYMCK,ISYMD)
C
                  KTRVI = 1
                  KEND1 = KTRVI + NT1AM(ISYMCK)
                  LWRK1 = LWORK - KEND1
C
                  DO B = 1, NVIR(ISYMB)
C
                     IOFF  = ICKBD(ISYCKB,ISYMD) 
     *                     + NCKATR(ISYCKB)*(D - 1) 
     *                     + ICKATR(ISYMCK,ISYMB)
     *                     + NT1AM(ISYMCK)*(B-1)
     *                     + 1
C
                     CALL GETWA2(LUTMP3,FNTMP3,WORK(KTRVI),IOFF,
     *                           NT1AM(ISYMCK))
C
                     IOFF  = ICKBD(ISYCKD,ISYMB) 
     *                     + NCKATR(ISYCKD)*(B - 1) 
     *                     + ICKATR(ISYMCK,ISYMD)
     *                     + NT1AM(ISYMCK)*(D-1)
     *                     + 1
C
                     CALL PUTWA2(LULDKBC3,FNLDKBC3,WORK(KTRVI),
     *                           IOFF,NT1AM(ISYMCK))
C
                  ENDDO
               ENDDO
C
            ENDDO
         ENDDO
C
C----------------------------
C     Final sort for new L
C----------------------------
C
         DO ISYMD = 1, NSYM
            ISYCKB = MULD2H(ISYINT,ISYMD)
C
            DO D = 1, NVIR(ISYMD)
C
               DO ISYMB = 1, NSYM
C
                  ISYMCK = MULD2H(ISYCKB,ISYMB)
                  ISYCKD = MULD2H(ISYMCK,ISYMD)
C
                  KTRVI = 1
                  KEND1 = KTRVI + NT1AM(ISYMCK)
                  LWRK1 = LWORK - KEND1
C
                  DO B = 1, NVIR(ISYMB)
C
                     IOFF  = ICKBD(ISYCKB,ISYMD) 
     *                     + NCKATR(ISYCKB)*(D - 1) 
     *                     + ICKATR(ISYMCK,ISYMB)
     *                     + NT1AM(ISYMCK)*(B-1)
     *                     + 1
C
                     CALL GETWA2(LUTMP4,FNTMP4,WORK(KTRVI),IOFF,
     *                           NT1AM(ISYMCK))
C
                     IOFF  = ICKBD(ISYCKD,ISYMB) 
     *                     + NCKATR(ISYCKD)*(B - 1) 
     *                     + ICKATR(ISYMCK,ISYMD)
     *                     + NT1AM(ISYMCK)*(D-1)
     *                     + 1
C
                     CALL PUTWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI),
     *                           IOFF,NT1AM(ISYMCK))
C
                  ENDDO
               ENDDO
C
            ENDDO
         ENDDO
C
C-------------------------------------------
C        Close and delete temporary files.
C-------------------------------------------
C
         CALL WCLOSE2(LUTMP1,FNTMP1,'DELETE')
         CALL WCLOSE2(LUTMP2,FNTMP2,'DELETE')
         CALL WCLOSE2(LUTMP3,FNTMP3,'DELETE')
         CALL WCLOSE2(LUTMP4,FNTMP4,'DELETE')
C
      ENDIF    ! END OF IOPT = 1
C
C----------------------------------------------
C     Sort integrals used in the contraction
C     for both IOPT=1 and IOPT=2
C----------------------------------------------
C
      DO ISYMD = 1, NSYM
         ISYCKB = MULD2H(ISYINT,ISYMD)
C
         DO ISYMB = 1, NSYM
C
            ISYMCK = MULD2H(ISYCKB,ISYMB)
            ISYCKD = MULD2H(ISYINT,ISYMB)
C
            IF (LWORK .LT. NT1AM(ISYMCK)) THEN
               CALL QUIT('Out of memory in CC3_TCME (2)')
            ENDIF
C
            KTRVI = 1
C
            DO D = 1, NVIR(ISYMD)
C
               DO B = 1, NVIR(ISYMB)
C
                  IF (NT1AM(ISYMCK) .GT. 0) THEN
C
                     IOFF  = ICKBD(ISYCKB,ISYMD) 
     *                     + NCKATR(ISYCKB)*(D - 1) 
     *                     + ICKATR(ISYMCK,ISYMB)
     *                     + NT1AM(ISYMCK)*(B-1)
     *                     + 1
C
                     CALL GETWA2(LULDKBC,FNLDKBC,WORK(KTRVI),IOFF,
     *                           NT1AM(ISYMCK))
C
                     IOFF  = ICKBD(ISYCKD,ISYMB) 
     *                     + NCKATR(ISYCKD)*(B - 1) 
     *                     + ICKATR(ISYMCK,ISYMD)
     *                     + NT1AM(ISYMCK)*(D-1)
     *                     + 1
C
                     CALL PUTWA2(LULDKBC4,FNLDKBC4,WORK(KTRVI),IOFF,
     *                           NT1AM(ISYMCK))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_TCME')
C
      RETURN
C
      END
C  /* Deck cc3_LSORT */
      SUBROUTINE CC3_LSORT1(INT1,ISYINT,WORK,LWORK,IOPT)
C
C     Written by K. Hald, Spring 2002.
C
C     Sort INT1 :
C     IOPT = 1  : Sort FROM ljka TO lkja
C     IOPT = 2  : Sort FROM ljka TO jlka
C     IOPT = 3  : Sort FROM ljka TO klja
C     IOPT = 4  : Sort FROM ljka TO jkla
C     IOPT = 5  : Sort FROM ljka TO kjla
C
      IMPLICIT NONE
C
      INTEGER ISYINT, LWORK, IOPT, KOFF1, KOFF2
      INTEGER ISYMA, ISYLJK, KTROC, KEND1, LWRK1, ISYMK, ISYMLJ
      INTEGER ISYMJ, ISYML, ISYMKL, ISYMKJ
C
#if defined (SYS_CRAY)
      REAL INT1(*), WORK(LWORK)
#else
      DOUBLE PRECISION INT1(*), WORK(LWORK)
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      CALL QENTER('CC3_LSORT1')
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2) .AND. (IOPT .NE. 3)
     *                  .AND. (IOPT .NE. 4) .AND. (IOPT .NE. 5)) THEN
         CALL QUIT('Wrong IOPT in CC3_LSORT1')
      ENDIF
C
C--------------------------
C     Sort.
C--------------------------
C
      DO ISYMA = 1, NSYM
         ISYLJK = MULD2H(ISYINT,ISYMA)
C
         KTROC = 1
         KEND1 = KTROC + NMAJIK(ISYLJK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_LSORT')
         ENDIF
C
         IF (NMAJIK(ISYLJK) .GT. 0) THEN
         DO A = 1, NVIR(ISYMA)
            DO ISYMK = 1, NSYM
               ISYMLJ = MULD2H(ISYLJK,ISYMK)
               DO ISYMJ = 1, NSYM
                  ISYML  = MULD2H(ISYMLJ,ISYMJ)
                  ISYMKL = MULD2H(ISYML,ISYMK)
                  ISYMKJ = MULD2H(ISYMJ,ISYMK)
C
                  DO K = 1, NRHF(ISYMK)
                  DO J = 1, NRHF(ISYMJ)
C
                     KOFF1 = ISJIKA(ISYLJK,ISYMA)
     *                     + NMAJIK(ISYLJK)*(A-1)
     *                     + ISJIK(ISYMLJ,ISYMK)
     *                     + NMATIJ(ISYMLJ)*(K - 1)
     *                     + IMATIJ(ISYML,ISYMJ)
     *                     + NRHF(ISYML)*(J - 1) 
     *                     + 1
C
                     IF (IOPT .EQ. 1) THEN
                        KOFF2 = KTROC - 1
     *                        + ISJIK(ISYMKL,ISYMJ)
     *                        + NMATIJ(ISYMKL)*(J - 1)
     *                        + IMATIJ(ISYML,ISYMK)
     *                        + NRHF(ISYML)*(K-1)
     *                        + 1
C
                           CALL DCOPY(NRHF(ISYML),INT1(KOFF1),1,
     *                                WORK(KOFF2),1)
C
                     ELSE IF (IOPT .EQ. 2) THEN
                        KOFF2 = KTROC - 1
     *                        + ISJIK(ISYMLJ,ISYMK)
     *                        + NMATIJ(ISYMLJ)*(K - 1)
     *                        + IMATIJ(ISYMJ,ISYML)
     *                        + J
C
                           CALL DCOPY(NRHF(ISYML),INT1(KOFF1),1,
     *                                WORK(KOFF2),NRHF(ISYMJ))
C
                     ELSE IF (IOPT .EQ. 3) THEN
                        KOFF2 = KTROC - 1
     *                        + ISJIK(ISYMKL,ISYMJ)
     *                        + NMATIJ(ISYMKL)*(J - 1)
     *                        + IMATIJ(ISYMK,ISYML)
     *                        + K
C
                           CALL DCOPY(NRHF(ISYML),INT1(KOFF1),1,
     *                                WORK(KOFF2),NRHF(ISYMK))
C
                     ELSE IF (IOPT .EQ. 4) THEN
                        KOFF2 = KTROC - 1
     *                        + ISJIK(ISYMKJ,ISYML)
     *                        + IMATIJ(ISYMJ,ISYMK)
     *                        + NRHF(ISYMJ)*(K-1)
     *                        + J
C
                           CALL DCOPY(NRHF(ISYML),INT1(KOFF1),1,
     *                                WORK(KOFF2),NMATIJ(ISYMKJ))
C
                     ELSE IF (IOPT .EQ. 5) THEN
                        KOFF2 = KTROC - 1
     *                        + ISJIK(ISYMKJ,ISYML)
     *                        + IMATIJ(ISYMK,ISYMJ)
     *                        + NRHF(ISYMK)*(J-1)
     *                        + K
C
                           CALL DCOPY(NRHF(ISYML),INT1(KOFF1),1,
     *                                WORK(KOFF2),NMATIJ(ISYMKJ))
                     ENDIF
C
                  ENDDO   ! J
                  ENDDO   ! K
               ENDDO      ! ISYMJ
            ENDDO         ! ISYMK
C
            KOFF1 = ISJIKA(ISYLJK,ISYMA)
     *            + NMAJIK(ISYLJK)*(A-1)
     *            + 1
C
            CALL DCOPY(NMAJIK(ISYLJK),WORK(KTROC),1,INT1(KOFF1),1)
C
         ENDDO            ! A
         ENDIF
      ENDDO               ! ISYMA
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CC3_LSORT1')
C
      RETURN
      END
C  /* Deck cc3_lsort2 */
      SUBROUTINE CC3_LSORT2(INT1,ISYINT,WORK,LWORK,IOPT)
C
C     Written by K. Hald, Spring 2002.
C
C     Sort INT1 :
C     IOPT = 1  : Sort FROM akjl TO ajkl
C     IOPT = 2  : Sort FROM akjl TO aklj
C     IOPT = 3  : Sort FROM akjl TO ajlk
C     IOPT = 4  : Sort FROM akjl TO aljk
C     IOPT = 5  : Sort FROM akjl TO alkj
C
      IMPLICIT NONE
C
      INTEGER ISYINT, LWORK, IOPT
      INTEGER ISYMA, ISYLJK, ISYMK, ISYMLJ, ISYMAK, ISYMJ, ISYML
      INTEGER ISYMAJ, ISYMAL, ISYAKJ, ISYAKL, ISYALJ, KOFF1, KOFF2
      INTEGER KTROC, KEND1, LWRK1
C
#if defined (SYS_CRAY)
      REAL INT1(*), WORK(LWORK)
#else
      DOUBLE PRECISION INT1(*), WORK(LWORK)
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      CALL QENTER('CC3_LSORT2')
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2) .AND. (IOPT .NE. 3)
     *                  .AND. (IOPT .NE. 4) .AND. (IOPT .NE. 5)) THEN
         CALL QUIT('Wrong IOPT in CC3_LSORT2')
      ENDIF
C
C--------------------------
C     Sort.
C--------------------------
C
      KTROC = 1
      KEND1 = KTROC + NTRAOC(ISYINT)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_LSORT2')
      ENDIF
C
      DO ISYMA = 1, NSYM
         ISYLJK = MULD2H(ISYINT,ISYMA)
         DO ISYMK = 1, NSYM
            ISYMLJ = MULD2H(ISYLJK,ISYMK)
            ISYMAK = MULD2H(ISYMA,ISYMK)
            DO ISYMJ = 1, NSYM
               ISYML  = MULD2H(ISYMLJ,ISYMJ)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
               ISYMAL = MULD2H(ISYMA,ISYML)
               ISYAKJ = MULD2H(ISYMAK,ISYMJ)
               ISYAKL = MULD2H(ISYMAK,ISYML)
               ISYALJ = MULD2H(ISYMAL,ISYMJ)
C
               DO A = 1, NVIR(ISYMA)
               DO K = 1, NRHF(ISYMK)
               DO J = 1, NRHF(ISYMJ)
C
                  KOFF1 = ISAIKJ(ISYAKJ,ISYML)
     *                  + ICKI(ISYMAK,ISYMJ)
     *                  + NT1AM(ISYMAK)*(J-1)
     *                  + IT1AM(ISYMA,ISYMK)
     *                  + NVIR(ISYMA)*(K-1)
     *                  + A
C
                  IF (IOPT .EQ. 1) THEN
                     KOFF2 = KTROC - 1
     *                     + ISAIKJ(ISYAKJ,ISYML)
     *                     + ICKI(ISYMAJ,ISYMK)
     *                     + NT1AM(ISYMAJ)*(K-1)
     *                     + IT1AM(ISYMA,ISYMJ)
     *                     + NVIR(ISYMA)*(J-1)
     *                     + A
C
                     CALL DCOPY(NRHF(ISYML),INT1(KOFF1),NCKI(ISYAKJ),
     *                          WORK(KOFF2),NCKI(ISYAKJ))
C
                  ELSE IF (IOPT .EQ. 2) THEN
                     KOFF2 = KTROC - 1
     *                     + ISAIKJ(ISYAKL,ISYMJ)
     *                     + NCKI(ISYAKL)*(J-1)
     *                     + ICKI(ISYMAK,ISYML)
     *                     + IT1AM(ISYMA,ISYMK)
     *                     + NVIR(ISYMA)*(K-1)
     *                     + A
C
                     CALL DCOPY(NRHF(ISYML),INT1(KOFF1),NCKI(ISYAKJ),
     *                          WORK(KOFF2),NT1AM(ISYMAK))
C
                  ELSE IF (IOPT .EQ. 3) THEN
                     KOFF2 = KTROC - 1
     *                     + ISAIKJ(ISYALJ,ISYMK)
     *                     + NCKI(ISYALJ)*(K-1)
     *                     + ICKI(ISYMAJ,ISYML)
     *                     + IT1AM(ISYMA,ISYMJ)
     *                     + NVIR(ISYMA)*(J-1)
     *                     + A
C
                     CALL DCOPY(NRHF(ISYML),INT1(KOFF1),NCKI(ISYAKJ),
     *                          WORK(KOFF2),NT1AM(ISYMAK))
C
                  ELSE IF (IOPT .EQ. 4) THEN
                     KOFF2 = KTROC - 1
     *                     + ISAIKJ(ISYALJ,ISYMK)
     *                     + NCKI(ISYALJ)*(K-1)
     *                     + ICKI(ISYMAL,ISYMJ)
     *                     + NT1AM(ISYMAL)*(J-1)
     *                     + IT1AM(ISYMA,ISYML)
     *                     + A
C
                     CALL DCOPY(NRHF(ISYML),INT1(KOFF1),NCKI(ISYAKJ),
     *                          WORK(KOFF2),NVIR(ISYMA))
C
                  ELSE IF (IOPT .EQ. 5) THEN
                     KOFF2 = KTROC - 1
     *                     + ISAIKJ(ISYAKL,ISYMJ)
     *                     + NCKI(ISYAKL)*(J-1)
     *                     + ICKI(ISYMAL,ISYMK)
     *                     + NT1AM(ISYMAL)*(K-1)
     *                     + IT1AM(ISYMA,ISYML)
     *                     + A
C
                     CALL DCOPY(NRHF(ISYML),INT1(KOFF1),NCKI(ISYAKJ),
     *                          WORK(KOFF2),NVIR(ISYMA))
C
                  ENDIF
C
               ENDDO
               ENDDO
               ENDDO
C
            ENDDO
         ENDDO
      ENDDO
C
      CALL DCOPY(NTRAOC(ISYINT),WORK(KTROC),1,INT1,1)
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CC3_LSORT2')
C
      RETURN
      END
C  /* Deck cc3_intstore */
      SUBROUTINE CC3_INTSTORE(LUO3,FNO3,XINT4O,ISYM4O,XLAMDH1,ISYMLH1,
     *                        XLAMDH2,ISYMLH2,LU3V,FN3V,LU4V,FN4V,
     *                        ISYM4V,WORK,LWORK,IOPT)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate g_{OOOO} (O=occ) integrals that are needed for the 
C     CC3 left hand side and return them in XINT4O.
C     Integrals g_{OOO,delta} are read from disc.
C
C     Calculate g_{VVVV} (V=vir) integrals that are needed for the
C     CC3 left hand side and store them on disc.
C     Integrals g_{VVV,delta} are read from disc.
C
C     IOPT = 1 -> DO g-oooo
C     IOPT = 2 -> DO g-vvvv
C     IOPT = 3 -> Do both
C
      IMPLICIT NONE
C
      INTEGER LUO3, ISYM4O, ISYMLH1, ISYMLH2, LU3V, LU4V, ISYM4V, LWORK
      INTEGER IOPT, ISYMD, ISYIJK, ISYML, NTOT, KAOINT, KEND1, LWRK1
      INTEGER IOFF, KOFF1, KOFF2, KOFF3, NTOIJK, NTODEL
      INTEGER ISYABC, ISYDEL, KRES1, IDEL, ISYMC, ISYMAB, ISYABD
C
#if defined (SYS_CRAY)
      REAL XINT4O(*), XLAMDH1(*), XLAMDH2(*)
      REAL WORK(LWORK), ZERO, ONE
      REAL DDOT, XNORM
#else
      DOUBLE PRECISION XINT4O(*), XLAMDH1(*), XLAMDH2(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      CHARACTER FNO3*(*), FN3V*(*), FN4V*(*)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('CC3_INTSTORE')
C
C-----------------------------
C     Sanity check
C-----------------------------
C
      IF ((IOPT .NE. 1) .AND. (IOPT.NE.2) .AND. (IOPT.NE.3)) THEN
         CALL QUIT('Wrong IOPT in CC3_INTSTORE')
      ENDIF
C
C-------------------------------------------------------------
C     Calculate integrals with 4 occupied indices
C-------------------------------------------------------------
C
      IF ((IOPT .EQ.1) .OR. (IOPT.EQ.3)) THEN
         DO ISYMD = 1, NSYM
            IF (NBAS(ISYMD) .GT. 0) THEN
C
               ISYML  = MULD2H(ISYMD,ISYMLH1)
               ISYIJK = MULD2H(ISYM4O,ISYML)
C
               NTOT = NMAIJK(ISYIJK)*NBAS(ISYMD)
C
               KAOINT = 1
               KEND1  = KAOINT + NTOT
               LWRK1  = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Out of memory in CC3_INTSTORE (g-OOOO)')
               ENDIF
C
C------------------------------------
C        Read integrals from disc.
C------------------------------------
C
               IOFF = I3ODEL(ISYIJK,ISYMD) + 1
C
               CALL GETWA2(LUO3,FNO3,WORK(KAOINT),IOFF,NTOT)
C
C-----------------------------------------
C        Transform to four MO index.
C-----------------------------------------
C
               KOFF2 = IGLMRH(ISYMD,ISYML) + 1
               KOFF3 = I3ORHF(ISYIJK,ISYML)
     *               + 1
C
               NTOIJK = MAX(1,NMAIJK(ISYIJK))
               NTODEL = MAX(1,NBAS(ISYMD))
C
               CALL DGEMM('N','N',NMAIJK(ISYIJK),NRHF(ISYML),
     *                    NBAS(ISYMD),ONE,WORK(KAOINT),NTOIJK,
     *                    XLAMDH1(KOFF2),NTODEL,ONE,
     *                    XINT4O(KOFF3),NTOIJK)
C
            ENDIF
C
         ENDDO
C
      ENDIF
C
C----------------------------------------------------
C     Calculate integrals with 4 virtual indices
C     and store on disc
C----------------------------------------------------
C
      IF ((IOPT.EQ.2) .OR. (IOPT.EQ.3)) THEN
         DO ISYMD = 1, NSYM
            IF (NVIR(ISYMD) .GT. 0) THEN
C
               ISYABC = MULD2H(ISYM4V,ISYMD)
               ISYDEL = MULD2H(ISYMD,ISYMLH2)
C
               DO D = 1, NVIR(ISYMD)
C
                  KAOINT = 1
                  KRES1  = KAOINT + NMAABC(ISYABC)
                  KEND1  = KRES1  + NMAABC(ISYABC)
                  LWRK1  = LWORK - KEND1
C
                  IF (LWRK1 .LT. 0) THEN
                     CALL QUIT('Out of memory in CC3_INTSTORE (g-VVVV)')
                  ENDIF
C
                  CALL DZERO(WORK(KRES1),NMAABC(ISYABC))
C
                  DO IDEL = 1, NBAS(ISYDEL)
C
                     IOFF = I3VDEL(ISYABC,ISYDEL) 
     *                    + NMAABC(ISYABC)*(IDEL-1) 
     *                    + 1
C
                     CALL GETWA2(LU3V,FN3V,WORK(KAOINT),IOFF,
     *                           NMAABC(ISYABC))
C
C-----------------------------------------
C        Transform to four MO index.
C-----------------------------------------
C
                     KOFF1 = IGLMVI(ISYDEL,ISYMD)
     *                     + NBAS(ISYDEL)*(D-1)
     *                     + IDEL
C
                     CALL DAXPY(NMAABC(ISYABC),XLAMDH2(KOFF1),
     *                          WORK(KAOINT),1,WORK(KRES1),1)
C
                  ENDDO
C
                  DO ISYMC = 1, NSYM
                     ISYMAB = MULD2H(ISYABC,ISYMC)
                     ISYABD = MULD2H(ISYMAB,ISYMD)
                     IF (NMATAB(ISYMAB) .GT. 0) THEN
                        DO C = 1, NVIR(ISYMC)
C
                           IOFF  = I3VVIR(ISYABD,ISYMC)
     *                           + NMAABC(ISYABD)*(C-1)
     *                           + IMAABC(ISYMAB,ISYMD)
     *                           + NMATAB(ISYMAB)*(D-1)
     *                           + 1
                           KOFF1 = KRES1
     *                           + IMAABC(ISYMAB,ISYMC)
     *                           + NMATAB(ISYMAB)*(C-1)
C
                           CALL PUTWA2(LU4V,FN4V,WORK(KOFF1),IOFF,
     *                                 NMATAB(ISYMAB))
C
                        ENDDO
                     ENDIF
                  ENDDO
C
C
               ENDDO
C
            ENDIF
C
         ENDDO
C
      ENDIF
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CC3_INTSTORE')
C
      RETURN
      END
C  /* Deck cc3_intdel */
      SUBROUTINE CC3_INTDEL(AOINT,ISYMAO,LUINT,FNINT,XLAMDP,ISYMLP,
     *                      XLAMDH,ISYMLH,ISYINT,WORK,LWORK,IDEL,ISYMD)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate integrals that are needed for the CC3 left hand side,
C     and store on file.
C
C     VVV,delta (V=vir.) are needed.
C
      IMPLICIT NONE
C
      INTEGER ISYMAO, LUINT, ISYMLP, ISYMLH, ISYINT, LWORK, IDEL, ISYMD
      INTEGER ISYABG, ISYTMP, ISYABC, KVVVV, KEND1, KEND2, LWRK1, LWRK2
      INTEGER ISYMG, ISYMC, ISALBE, ISYMAB, KINT, KSCR1, KSCR2
      INTEGER KOFF1, KOFF2, KOFF3, ISYMB, ISYMBE, ISYMAL, ISYMA
      INTEGER NBASAL, NBASBE, NVIRA, NAB, NBASG, IOFF
C
#if defined (SYS_CRAY)
      REAL AOINT(*), XLAMDP(*), XLAMDH(*)
      REAL WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION AOINT(*), XLAMDP(*), XLAMDH(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      CHARACTER FNINT*(*)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('CC3_INTDEL')
C
C-------------------------------------------
C     Work space allocation.
C-------------------------------------------
C
      ISYABG = MULD2H(ISYMAO,ISYMD)
C
      ISYTMP = MULD2H(ISYINT,ISYMD)
      ISYABC = MULD2H(ISYTMP,ISYMLH)
C
      KVVVV  = 1
      KEND1  = KVVVV  + NMAABC(ISYABC)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_INTDEL')
      ENDIF
C
      CALL DZERO(WORK(KVVVV),NMAABC(ISYABC))
C
C---------------------------------------------
C     Transform AO-integrals to g_{vvv,delta}
C---------------------------------------------
C
      DO ISYMG = 1, NSYM
         ISYMC  = MULD2H(ISYMG,ISYMLP)
         ISALBE = MULD2H(ISYABG,ISYMG)
         ISYMAB = MULD2H(ISYABC,ISYMC)
         ISYTMP = MULD2H(ISYMAB,ISYMLH)
C
         KINT   = KEND1
         KSCR1  = KINT  + NMATAB(ISYMAB)*NBAS(ISYMG)
         KSCR2  = KSCR1 + N2BST(ISALBE)
         KEND2  = KSCR2 + NEMAT1(ISYTMP)
         LWRK2  = LWORK - KEND2
COMMENT
COMMENT  allocate to much space for kscr2 at the moment
COMMENT
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_INTDEL (2)')
         ENDIF
C
         DO G = 1, NBAS(ISYMG)
C
            KOFF1 = IDSAOG(ISYMG,ISYMD) + NNBST(ISALBE)*(G-1) + 1
            CALL CCSD_SYMSQ(AOINT(KOFF1),ISALBE,WORK(KSCR1))
C
            DO ISYMB = 1,NSYM
C
               ISYMBE = MULD2H(ISYMB,ISYMLH)
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMA  = MULD2H(ISYMAL,ISYMLP)
C
               KOFF1 = KSCR1 
     *               + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMVI(ISYMBE,ISYMB) + 1
               KOFF3 = KSCR2
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMB),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF1),NBASAL,XLAMDH(KOFF2),NBASBE,
     *                    ZERO,WORK(KOFF3),NBASAL)
C
               KOFF1 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF2 = KSCR2
               KOFF3 = KINT 
     *               + NMATAB(ISYMAB)*(G - 1)
     *               + IMATAB(ISYMA,ISYMB)
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMAL),
     *                    ONE,XLAMDP(KOFF1),NBASAL,WORK(KOFF2),NBASAL,
     *                    ZERO,WORK(KOFF3),NVIRA)
C
            ENDDO
C
         ENDDO
C
         KOFF2 = IGLMVI(ISYMG,ISYMC)  + 1
         KOFF3 = KVVVV 
     *         + IMAABC(ISYMAB,ISYMC)
C
         NAB    = MAX(NMATAB(ISYMAB),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NMATAB(ISYMAB),NVIR(ISYMC),NBAS(ISYMG),
     *              ONE,WORK(KINT),NAB,XLAMDP(KOFF2),NBASG,
     *              ONE,WORK(KOFF3),NAB)
C
      ENDDO
C
C----------------------------------------
C     Save the g_{vvv,delta} to disc.
C----------------------------------------
C
      IF (NMAABC(ISYABC) .GT. 0) THEN
         KOFF1 = IDEL - IBAS(ISYMD)
         IOFF = I3VDEL(ISYABC,ISYMD) + NMAABC(ISYABC)*(KOFF1-1) + 1
         CALL PUTWA2(LUINT,FNINT,WORK(KVVVV),IOFF,NMAABC(ISYABC))
      ENDIF
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CC3_INTDEL')
C
      RETURN
      END
C  /* Deck cc3_2o2v */
      SUBROUTINE CC3_2O2V(AOINT,ISYMAO,DSRHF,ISYMDS,XOVVO,XOOVV,
     *                    XLAMP0,ISYMLP0,XLAMH0,ISYMLH0,
     *                    XLAMP1,ISYMLP1,XLAMH1,ISYMLH1,
     *                    ISYINT,WORK,LWORK,IDEL,ISYMD)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate integrals that are needed for the CC3 left hand side.
C
C     VOOV, VVOO (O=occ. V=vir.) are needed.
C
C     (k^p0 l^h1 | c^p1 d^h0)
C     (c^p1 k^h1 | l^p0 d^h0)
C                 ^^^
C            Transformed outside (DSRHF ... cannot be barred)
C
      IMPLICIT NONE
C
      INTEGER ISYMAO, ISYMDS, ISYMLP0, ISYMLH0, ISYMLP1, ISYMLH1
      INTEGER ISYINT, LWORK, IDEL, ISYMD
      INTEGER ISYABJ, ISYTMP, ISYAIJ, KVOO, KEND0, LWRK0, KEND1, LWRK1
      INTEGER KEND2, LWRK2, ISYMJ, ISALBE, ISYMAI, KSCR1, KSCR2
      INTEGER ISYMI, ISYMAL, ISYMBE, ISYMA, KOFF1, KOFF2, KOFF3
      INTEGER NTOTAL, NTOTA, NTOTB, ISYMB
      INTEGER ISYABG, KOOV, ISYMG, ISYMIJ, KINT, NBASAL, NBASBE
      INTEGER NRHFI, NIJ, NBASG, ISAIJD
C
#if defined (SYS_CRAY)
      REAL AOINT(*), DSRHF(*), XOVVO(*), XOOVV(*)
      REAL XLAMP0(*), XLAMH0(*)
      REAL XLAMP1(*), XLAMH1(*)
      REAL WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION AOINT(*), DSRHF(*), XOVVO(*), XOOVV(*)
      DOUBLE PRECISION XLAMP0(*), XLAMH0(*)
      DOUBLE PRECISION XLAMP1(*), XLAMH1(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('CC3_2O2V')
C
C=========================================
C     Calculate the integrals g_{voov}
C=========================================
C
      ISYABJ = MULD2H(ISYMDS,ISYMD)
      ISYTMP = MULD2H(ISYABJ,ISYMLP1)
      ISYAIJ = MULD2H(ISYTMP,ISYMLH1)
C
      KVOO  = 1
      KEND0 = KVOO  + NCKI(ISYAIJ)
      LWRK0 = LWORK - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('0-Insufficient work space area in CC3_2O2V')
      ENDIF
C
      CALL DZERO(WORK(KVOO),NCKI(ISYAIJ))
C
      DO ISYMJ = 1,NSYM
C
         ISALBE = MULD2H(ISYABJ,ISYMJ)
         ISYMAI = MULD2H(ISYAIJ,ISYMJ)
C
         DO J = 1,NRHF(ISYMJ)
C
C------------------------------------------------------------
C           Work space allocation 1 * unpacking of integrals.
C------------------------------------------------------------
C
            KSCR1  = KEND0
            KEND1  = KSCR1  + N2BST(ISALBE)
            LWRK1  = LWORK  - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('1-Insufficient work space area in CC3_2O2V')
            ENDIF
C
            KOFF1 = IDSRHF(ISALBE,ISYMJ) + NNBST(ISALBE)*(J - 1) + 1
C
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISALBE,WORK(KSCR1))
C
            DO ISYMI = 1,NSYM
C
C-----------------------------------------------------------------------
C              Transform remaining AO-indices of integrals to occ. space
C-----------------------------------------------------------------------
C
               ISYMBE = MULD2H(ISYMI,ISYMLH1)
               ISYMAL = MULD2H(ISALBE,ISYMBE)
               ISYMA  = MULD2H(ISYMAL,ISYMLP1)
C
               KSCR2  = KEND1
               KEND2  = KSCR2  + NBAS(ISYMAL)*NRHF(ISYMI)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('2-Insufficient work space area (CC3_2O2V)')
               ENDIF
C
               KOFF1 = KSCR1  + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMRH(ISYMBE,ISYMI) + 1
               KOFF3 = KSCR2 
C
               NTOTA = MAX(NBAS(ISYMAL),1)
               NTOTB = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF1),NTOTA,XLAMH1(KOFF2),NTOTB,
     *                    ZERO,WORK(KOFF3),NTOTA)
C
               KOFF1 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF2 = KSCR2
               KOFF3 = KVOO
     *               + ICKI(ISYMAI,ISYMJ)
     *               + NT1AM(ISYMAI)*(J-1)
     *               + IT1AM(ISYMA,ISYMI)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMAL),
     *                    ONE,XLAMP1(KOFF1),NTOTAL,WORK(KOFF2),NTOTAL,
     *                    ZERO,WORK(KOFF3),NTOTA)
C
            ENDDO
C
         ENDDO
C
      ENDDO
C
C--------------------------------------------------------
C     Transform the current delta index to virtual mo
C--------------------------------------------------------
C
      ISYMB = MULD2H(ISYMD,ISYMLH0)
C
      DO B = 1, NVIR(ISYMB)
C
         KOFF1 = IGLMVI(ISYMD,ISYMB)
     *         + NBAS(ISYMD)*(B-1)
     *         + (IDEL - IBAS(ISYMD))
         KOFF2 = KVOO
         KOFF3 = IT2SP(ISYAIJ,ISYMB)
     *         + NCKI(ISYAIJ)*(B-1)
     *         + 1
C
         CALL DAXPY(NCKI(ISYAIJ),XLAMH0(KOFF1),
     *              WORK(KOFF2),1,XOVVO(KOFF3),1)
C
      ENDDO
C
C=========================================
C     Calculate the integrals g_{oovv}
C=========================================
C
      ISYABG = MULD2H(ISYMAO,ISYMD)
C
      ISAIJD = MULD2H(ISYINT,ISYMLH0)
      ISYAIJ = MULD2H(ISAIJD,ISYMD)
C
      KOOV  = 1
      KEND1  = KOOV  + NCKI(ISYAIJ)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_2O2V (g_{oovv})')
      ENDIF
C
      CALL DZERO(WORK(KOOV),NCKI(ISYAIJ))
C
      DO ISYMG = 1, NSYM
         ISYMA  = MULD2H(ISYMG,ISYMLP1)
         ISALBE = MULD2H(ISYABG,ISYMG)
         ISYMIJ = MULD2H(ISYAIJ,ISYMA)
         ISYTMP = MULD2H(ISYMIJ,ISYMLH1)
C
         KINT   = KEND1
         KSCR1  = KINT  + NMATIJ(ISYMIJ)*NBAS(ISYMG)
         KSCR2  = KSCR1 + N2BST(ISALBE)
         KEND2  = KSCR2 + NT1AO(ISYTMP)
         LWRK2  = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_2O2V (2)')
         ENDIF
C
         DO G = 1, NBAS(ISYMG)
C
            KOFF1 = IDSAOG(ISYMG,ISYMD) + NNBST(ISALBE)*(G-1) + 1
            CALL CCSD_SYMSQ(AOINT(KOFF1),ISALBE,WORK(KSCR1))
C
            DO ISYMJ = 1,NSYM
C
               ISYMBE = MULD2H(ISYMJ,ISYMLH1)
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMI  = MULD2H(ISYMAL,ISYMLP0)
C
               KOFF1 = KSCR1 
     *               + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMRH(ISYMBE,ISYMJ) + 1
               KOFF3 = KSCR2
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF1),NBASAL,XLAMH1(KOFF2),NBASBE,
     *                    ZERO,WORK(KOFF3),NBASAL)
C
               KOFF1 = IGLMRH(ISYMAL,ISYMI) + 1
               KOFF2 = KSCR2
               KOFF3 = KINT 
     *               + NMATIJ(ISYMIJ)*(G - 1)
     *               + IMATIJ(ISYMI,ISYMJ)
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NRHFI  = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYMAL),
     *                    ONE,XLAMP0(KOFF1),NBASAL,WORK(KOFF2),NBASAL,
     *                    ZERO,WORK(KOFF3),NRHFI)
C
            ENDDO
C
         ENDDO
C
         KOFF2 = IGLMVI(ISYMG,ISYMA)  + 1
         KOFF3 = KOOV 
     *         + IMAIJA(ISYMIJ,ISYMA)
C
         NIJ    = MAX(NMATIJ(ISYMIJ),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NMATIJ(ISYMIJ),NVIR(ISYMA),NBAS(ISYMG),
     *              ONE,WORK(KINT),NIJ,XLAMP1(KOFF2),NBASG,
     *              ONE,WORK(KOFF3),NIJ)
C
      ENDDO
C
C--------------------------------------------------------------
C     Do the final contraction of delta and store in XOOVV
C--------------------------------------------------------------
C
      ISYMB = MULD2H(ISYMD,ISYMLH0)
C
      DO B = 1, NVIR(ISYMB)
C
         KOFF1 = IGLMVI(ISYMD,ISYMB)
     *         + NBAS(ISYMD)*(B-1)
     *         + (IDEL - ibas(isymd))
C
         DO ISYMA = 1, NSYM
            ISYMIJ = MULD2H(ISYAIJ,ISYMA)
            DO ISYMI = 1, NSYM
               ISYMJ  = MULD2H(ISYMIJ,ISYMI)
               ISYMAI = MULD2H(ISYMA,ISYMI)
               DO A = 1, NVIR(ISYMA)
                  DO I = 1, NRHF(ISYMI)
C
                     KOFF2 = KOOV - 1
     *                     + IMAIJA(ISYMIJ,ISYMA)
     *                     + NMATIJ(ISYMIJ)*(A-1)
     *                     + IMATIJ(ISYMI,ISYMJ)
     *                     + I
                     KOFF3 = IT2SP(ISYAIJ,ISYMB)
     *                     + NCKI(ISYAIJ)*(B-1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     CALL DAXPY(NRHF(ISYMJ),XLAMH0(KOFF1),
     *                          WORK(KOFF2),NRHF(ISYMI),
     *                          XOOVV(KOFF3),NT1AM(ISYMAI))
                  ENDDO
               ENDDO
C
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CC3_2O2V')
C
      RETURN
      END
C  /* Deck cc3_l3_omega1 */
      SUBROUTINE CC3_L3_OMEGA1(OMEGA1,ISYRES,SMAT,QMAT,TMAT,ISYMIM,
     *                         XOOOO,XOVVO,XOOVV,XVVVV,ISYINT,
     *                         T2TP,ISYMT2,WORK,LWORK,LENSQ,INDSQ,
     *                         ISYMIB,IB,ISYMID,ID)
C
C     Written by K. Hald, Spring 2002.
C
C     Calculate the L3 contributions to omega1.
C
      IMPLICIT NONE
C
      INTEGER ISYMIM, ISYINT, ISYMT2, LWORK, ISYMIB, IB, ISYMID, ID
      INTEGER LENSQ, INDSQ(LENSQ,6)
      INTEGER ISYMBD, ISCKIJ, ISYCKM, LENGTH, ISYTMP, KSCR1
      INTEGER KEND1, LWRK1, KEND2, LWRK2, ISYMCK, ISYMIJ, ISYMDM
      INTEGER ISYMM, KOFF1, KOFF2, KOFF3, NTOTIJ, NTOTCK
      INTEGER KT2TMP, ISYMCM, ISYMK, ISYMC
      INTEGER ISYRES, ISYMI, ISYOOO, NTOIJK, NTOTB, NTOTI, NBI
      INTEGER ISYVVV, ISYEIJ, ISYMKM, ISYME, KSCR2, ISYMEK
      INTEGER ISYMCE, ISYMAC, ISYMA, NTOTCE, NTOTA, ISYENI, ISYMEN
      INTEGER ISYDLM, ISYMN, NTODLM, NTOTE, ISYMDN, ISYDNI
      INTEGER NTOTEN, ISYENF, ISYELM, ISYMLM, ISYML, ISYMFN, ISYFNI
      INTEGER ISYMEL, ISYLMI, NTOTFN, NTOTLM, ISYMEI, KSCR3, ISYMF
      INTEGER ISYMFI, ISYMDL, ISYMIN, NTOTDL, NTOTIN, ISYMD, ISYTMP2
      INTEGER ISYMMN, ISYAMN, ISYDMN, NTOTMN, ISYBMN, ISYMBN, ISYMDI
C
#if defined (SYS_CRAY)
      REAL ZERO, ONE
      REAL OMEGA1(*), SMAT(*), QMAT(*), TMAT(*)
      REAL XOOOO(*), XOVVO(*), XOOVV(*), XVVVV(*)
      REAL T2TP(*), WORK(LWORK)
      REAL ZERO, ONE, DDOT, XNORM
#else
      DOUBLE PRECISION OMEGA1(*), SMAT(*), QMAT(*), TMAT(*)
      DOUBLE PRECISION XOOOO(*), XOVVO(*), XOOVV(*), XVVVV(*)
      DOUBLE PRECISION T2TP(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, DDOT, XNORM
#endif
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO= 0.0D0, ONE= 1.0D0)
C
      CALL QENTER('CC3_L3_OMEGA1')
C
      ISYTMP  = MULD2H(ISYMIM,ISYMT2)
      ISYTMP2 = MULD2H(ISYTMP,ISYINT)
      IF (ISYRES .NE. ISYTMP2) THEN
         CALL QUIT('Symmetry mimatch in CC3_L3_OMEGA1')
      ENDIF
C
      ISYMBD = MULD2H(ISYMIB,ISYMID)
      ISCKIJ = MULD2H(ISYMIM,ISYMBD)
C
      LENGTH = NCKIJ(ISCKIJ)
C
C================================================
C     Calculate contribution from g_{oooo}
C================================================
C
      ISYCKM = MULD2H(ISYMT2,ISYMID)
      ISYTMP = MULD2H(ISCKIJ,ISYCKM) ! Symmetry of intermediate
C
      KSCR1 = 1
      KEND1 = KSCR1 + NMAIJK(ISYTMP)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NMAIJK(ISYTMP))
C
C--------------------------------------------
C     First contribution to intermediate
C--------------------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(I)
     *           + QMAT(INDSQ(I,3))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-1)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      DO ISYMCK = 1, NSYM
C
         ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
         ISYMM  = MULD2H(ISYCKM,ISYMCK)
C
         KOFF1  = ISAIKL(ISYMCK,ISYMIJ)
     *          + 1
         KOFF2  = IT2SP(ISYCKM,ISYMID)
     *          + NCKI(ISYCKM)*(ID-1)
     *          + ICKI(ISYMCK,ISYMM)
     *          + 1
         KOFF3  = KSCR1
     *          + IMAIJK(ISYMIJ,ISYMM)
C
         NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
         NTOTCK = MAX(1,NT1AM(ISYMCK))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIJ),NRHF(ISYMM),NT1AM(ISYMCK),
     *              ONE,TMAT(KOFF1),NTOTCK,T2TP(KOFF2),NTOTCK,
     *              ONE,WORK(KOFF3),NTOTIJ)
C
      ENDDO
C
C--------------------------------------------
C     Second contribution to intermediate
C--------------------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,1))
     *           + QMAT(INDSQ(I,2))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-2)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      KT2TMP = KEND1
      KEND2  = KT2TMP + NCKI(ISYCKM)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Memory exceeded in CC3_L3_OMEGA1')
      ENDIF
C
      DO ISYMK = 1, NSYM
         ISYMCM = MULD2H(ISYCKM,ISYMK)
         DO ISYMC = 1, NSYM
            ISYMM  = MULD2H(ISYMCM,ISYMC)
            ISYMCK = MULD2H(ISYMC,ISYMK)
C
            DO K = 1, NRHF(ISYMK)
               DO M = 1, NRHF(ISYMM)
C
                  KOFF1 = IT2SP(ISYCKM,ISYMID)
     *                  + NCKI(ISYCKM)*(ID-1)
     *                  + ICKI(ISYMCK,ISYMM)
     *                  + NT1AM(ISYMCK)*(M-1)
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K-1)
     *                  + 1
                  KOFF2 = KT2TMP
     *                  + ICKI(ISYMCM,ISYMK)
     *                  + NT1AM(ISYMCM)*(K-1)
     *                  + IT1AM(ISYMC,ISYMM)
     *                  + NVIR(ISYMC)*(M-1)
C
                  CALL DCOPY(NVIR(ISYMC),T2TP(KOFF1),1,WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYMCK = 1, NSYM
C
         ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
         ISYMM  = MULD2H(ISYCKM,ISYMCK)
C
         KOFF1  = ISAIKL(ISYMCK,ISYMIJ)
     *          + 1
         KOFF2  = KT2TMP
     *          + ICKI(ISYMCK,ISYMM)
         KOFF3  = KSCR1
     *          + IMAIJK(ISYMIJ,ISYMM)
C
         NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
         NTOTCK = MAX(1,NT1AM(ISYMCK))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIJ),NRHF(ISYMM),NT1AM(ISYMCK),
     *              ONE,TMAT(KOFF1),NTOTCK,WORK(KOFF2),NTOTCK,
     *              ONE,WORK(KOFF3),NTOTIJ)
C
      ENDDO
C
C------------------------------------------------
C     Contract the intermediate with g_{oooo}
C------------------------------------------------
C
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYOOO = MULD2H(ISYINT,ISYMI)
C
      DO I = 1, NRHF(ISYMI)
         NBI = IT1AM(ISYMIB,ISYMI) + NVIR(ISYMIB)*(I-1) + IB
         KOFF1 = I3ORHF(ISYOOO,ISYMI)
     *         + NMAIJK(ISYOOO)*(I-1)
     *         + 1
         OMEGA1(NBI) = OMEGA1(NBI) 
     *               + DDOT(NMAIJK(ISYOOO),XOOOO(KOFF1),1,WORK(KSCR1),1)
      ENDDO
C
C=============================================
C     Calculate contribution from g_{vvvv}
C=============================================
C
      ISYCKM = MULD2H(ISYMT2,ISYMIB)
      ISYEIJ = ISYCKM
      ISYTMP = MULD2H(ISCKIJ,ISYEIJ)
C
      KSCR1  = 1
      KT2TMP = KSCR1  + NCKATR(ISYTMP)
      KEND1  = KT2TMP + NCKI(ISYEIJ)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (VVVV T2-sort)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NCKATR(ISYTMP))
C
C----------------
C     Sort T2
C----------------
C
      DO ISYMK = 1, NSYM
         ISYMCM = MULD2H(ISYCKM,ISYMK)
         DO ISYMC = 1, NSYM
            ISYMM = MULD2H(ISYMCM,ISYMC)
            ISYMKM = MULD2H(ISYMK,ISYMM)
            ISYMCK = MULD2H(ISYMK,ISYMC)
C
            DO K = 1, NRHF(ISYMK)
               DO M = 1, NRHF(ISYMM)
                  KOFF1 = IT2SP(ISYCKM,ISYMIB)
     *                  + NCKI(ISYCKM)*(IB-1)
     *                  + ICKI(ISYMCK,ISYMM)
     *                  + NT1AM(ISYMCK)*(M-1)
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K-1)
     *                  + 1
                  KOFF2 = KT2TMP - 1
     *                  + IMAIJA(ISYMKM,ISYMC)
     *                  + IMATIJ(ISYMK,ISYMM)
     *                  + NRHF(ISYMK)*(M-1)
     *                  + K
C
                  CALL DCOPY(NVIR(ISYMC),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMKM))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------------
C     First intermediate
C--------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(I)
     *           + QMAT(INDSQ(I,3))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-3)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      DO ISYMCK = 1, NSYM
         ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
         ISYME  = MULD2H(ISYEIJ,ISYMIJ)
C
         KOFF1 = ISAIKL(ISYMCK,ISYMIJ)
     *         + 1
         KOFF2 = KT2TMP
     *         + IMAIJA(ISYMIJ,ISYME)
         KOFF3 = KSCR1
     *         + ICKATR(ISYMCK,ISYME)
C
         NTOTCK = MAX(1,NT1AM(ISYMCK))
         NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
C
         CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYME),
     *              NMATIJ(ISYMIJ),ONE,TMAT(KOFF1),NTOTCK,
     *              WORK(KOFF2),NTOTIJ,ONE,WORK(KOFF3),
     *              NTOTCK)
      ENDDO
C
C---------------------
C     Sort result.
C---------------------
C
      KSCR2  = KEND1
      KEND2  = KSCR2  + NCKATR(ISYTMP)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (Sorting-1)')
      ENDIF
C
      DO ISYMC = 1, NSYM
         ISYMEK = MULD2H(ISYMC,ISYTMP)
         DO ISYMK = 1, NSYM
            ISYME  = MULD2H(ISYMK,ISYMEK)
            ISYMCE = MULD2H(ISYMC,ISYME)
            ISYMCK = MULD2H(ISYMC,ISYMK)
C
            DO K = 1, NRHF(ISYMK)
               DO E = 1, NVIR(ISYME)
C
                  KOFF1 = KSCR1
     *                  + ICKATR(ISYMCK,ISYME)
     *                  + NT1AM(ISYMCK)*(E-1)
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K-1)
                  KOFF2 = KSCR2
     *                  + ICKASR(ISYMCE,ISYMK)
     *                  + NMATAB(ISYMCE)*(K-1)
     *                  + IMATAB(ISYMC,ISYME)
     *                  + NVIR(ISYMC)*(E-1)
C
                  CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),1,
     *                       WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL DCOPY(NCKATR(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
C----------------------------------------
C     Sort and contract with integral.
C----------------------------------------
C
      ISYVVV = MULD2H(ISYINT,ISYMID)
C
      KSCR2 = KEND1
      KEND2 = KSCR2 + NMAABC(ISYVVV)
      LWRK2 = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (Sort/Contract)')
      ENDIF
C
      DO ISYME = 1, NSYM
         ISYMAC = MULD2H(ISYVVV,ISYME)
         DO ISYMC = 1, NSYM
            ISYMA  = MULD2H(ISYMAC,ISYMC)
            ISYMCE = MULD2H(ISYMC,ISYME)
C
            DO A = 1, NVIR(ISYMA)
               DO E = 1, NVIR(ISYME)
C
                  KOFF1 = IMAABC(ISYMAC,ISYME)
     *                  + NMATAB(ISYMAC)*(E-1)
     *                  + IMATAB(ISYMC,ISYMA)
     *                  + NVIR(ISYMC)*(A-1)
     *                  + 1
                  KOFF2 = KSCR2
     *                  + IMAABC(ISYMCE,ISYMA)
     *                  + NMATAB(ISYMCE)*(A-1)
     *                  + IMATAB(ISYMC,ISYME)
     *                  + NVIR(ISYMC)*(E-1)
C
                  CALL DCOPY(NVIR(ISYMC),XVVVV(KOFF1),1,
     *                       WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYMA,ISYRES)
         ISYMCE = MULD2H(ISYMA,ISYVVV)
C
         KOFF1 = KSCR2
     *         + IMAABC(ISYMCE,ISYMA)
         KOFF2 = KSCR1
     *         + ICKASR(ISYMCE,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NTOTCE = MAX(1,NMATAB(ISYMCE))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATAB(ISYMCE),
     *              ONE,WORK(KOFF1),NTOTCE,WORK(KOFF2),NTOTCE,
     *              ONE,OMEGA1(KOFF3),NTOTA)
      ENDDO
C
C---------------------------------------
C     Second contribution.
C---------------------------------------
C
      ISYCKM = MULD2H(ISYMT2,ISYMIB)
      ISYEIJ = ISYCKM
      ISYTMP = MULD2H(ISCKIJ,ISYEIJ)
C
      KSCR1  = 1
      KT2TMP = KSCR1  + NCKATR(ISYTMP)
      KEND1  = KT2TMP + NCKI(ISYCKM)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (T2-sort-2)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NCKATR(ISYTMP))
C
C----------------
C     Sort T2
C----------------
C
      DO ISYMK = 1, NSYM
         ISYMCM = MULD2H(ISYCKM,ISYMK)
         DO ISYMC = 1, NSYM
            ISYMM = MULD2H(ISYMCM,ISYMC)
            ISYMKM = MULD2H(ISYMK,ISYMM)
            ISYMCK = MULD2H(ISYMK,ISYMC)
C
            DO K = 1, NRHF(ISYMK)
               DO M = 1, NRHF(ISYMM)
                  KOFF1 = IT2SP(ISYCKM,ISYMIB)
     *                  + NCKI(ISYCKM)*(IB-1)
     *                  + ICKI(ISYMCK,ISYMM)
     *                  + NT1AM(ISYMCK)*(M-1)
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K-1)
     *                  + 1
                  KOFF2 = KT2TMP - 1
     *                  + IMAIJA(ISYMKM,ISYMC)
     *                  + IMATIJ(ISYMK,ISYMM)
     *                  + NRHF(ISYMK)*(M-1)
     *                  + K
C
                  CALL DCOPY(NVIR(ISYMC),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMKM))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------------
C     Second intermediate
C--------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,1))
     *           + QMAT(INDSQ(I,2))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-3)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      DO ISYMCK = 1, NSYM
         ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
         ISYME  = MULD2H(ISYEIJ,ISYMIJ)
C
         KOFF1 = ISAIKL(ISYMCK,ISYMIJ)
     *         + 1
         KOFF2 = KT2TMP
     *         + IMAIJA(ISYMIJ,ISYME)
         KOFF3 = KSCR1
     *         + ICKATR(ISYMCK,ISYME)
C
         NTOTCK = MAX(1,NT1AM(ISYMCK))
         NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
C
         CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYME),
     *              NMATIJ(ISYMIJ),ONE,TMAT(KOFF1),NTOTCK,
     *              WORK(KOFF2),NTOTIJ,ONE,WORK(KOFF3),
     *              NTOTCK)
      ENDDO
C
C---------------------
C     Sort result.
C---------------------
C
      KSCR2  = KEND1
      KEND2  = KSCR2  + NCKATR(ISYTMP)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (Sorting-1)')
      ENDIF
C
      DO ISYMC = 1, NSYM
         ISYMEK = MULD2H(ISYMC,ISYTMP)
         DO ISYMK = 1, NSYM
            ISYME  = MULD2H(ISYMK,ISYMEK)
            ISYMCE = MULD2H(ISYMC,ISYME)
            ISYMCK = MULD2H(ISYMC,ISYMK)
C
            DO K = 1, NRHF(ISYMK)
               DO E = 1, NVIR(ISYME)
C
                  KOFF1 = KSCR1 - 1
     *                  + ICKATR(ISYMCK,ISYME)
     *                  + NT1AM(ISYMCK)*(E-1)
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K-1)
     *                  + 1
                  KOFF2 = KSCR2 - 1
     *                  + ICKASR(ISYMCE,ISYMK)
     *                  + NMATAB(ISYMCE)*(K-1)
     *                  + IMATAB(ISYMC,ISYME)
     *                  + NVIR(ISYMC)*(E-1)
     *                  + 1
C
                  CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),1,
     *                       WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL DCOPY(NCKATR(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
C----------------------------------------
C     Contract with integral.
C----------------------------------------
C
      ISYVVV = MULD2H(ISYINT,ISYMID)
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYMA,ISYRES)
         ISYMCE = MULD2H(ISYMA,ISYVVV)
C
         KOFF1 = IMAABC(ISYMCE,ISYMA)
     *         + 1
         KOFF2 = KSCR1
     *         + ICKASR(ISYMCE,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NTOTCE = MAX(1,NMATAB(ISYMCE))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATAB(ISYMCE),
     *              ONE,XVVVV(KOFF1),NTOTCE,WORK(KOFF2),NTOTCE,
     *              ONE,OMEGA1(KOFF3),NTOTA)
      ENDDO
C
C================================================
C     - L^{daf}_{lmn} t^{de}_{lm} g_{iefn}
C================================================
C
C--------------------------------
C     First contribution
C--------------------------------
C
      ISYMEN = MULD2H(ISCKIJ,ISYMT2)
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYENI = MULD2H(ISYMEN,ISYMI)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NT1AM(ISYMEN)
      KEND1 = KSCR2 + NCKI(ISYENI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (OVVO-1)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NT1AM(ISYMEN))
      CALL DZERO(WORK(KSCR2),NCKI(ISYENI))
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,3))
     *           + QMAT(I)
      ENDDO
C
      DO ISYME = 1, NSYM
         ISYDLM = MULD2H(ISYME,ISYMT2)
         ISYMN  = MULD2H(ISYMEN,ISYME)
C
         KOFF1 = IT2SP(ISYDLM,ISYME)
     *         + 1
         KOFF2 = ISAIKJ(ISYDLM,ISYMN)
     *         + 1
         KOFF3 = KSCR1
     *         + IT1AM(ISYME,ISYMN)
C
         NTODLM = MAX(1,NCKI(ISYDLM))
         NTOTE  = MAX(1,NVIR(ISYME))
C
         CALL DGEMM('T','N',NVIR(ISYME),NRHF(ISYMN),NCKI(ISYDLM),
     *              -ONE,T2TP(KOFF1),NTODLM,TMAT(KOFF2),NTODLM,
     *              ONE,WORK(KOFF3),NTOTE)
C
      ENDDO
C
      DO ISYME = 1, NSYM
         ISYMN  = MULD2H(ISYMEN,ISYME)
         ISYMDN = MULD2H(ISYMN,ISYMID)
         ISYDNI = MULD2H(ISYMDN,ISYMI)
C
         DO E = 1, NVIR(ISYME)
            DO N = 1, NRHF(ISYMN)
C
               KOFF1 = IT2SP(ISYDNI,ISYME)
     *               + NCKI(ISYDNI)*(E-1)
     *               + ICKI(ISYMDN,ISYMI)
     *               + IT1AM(ISYMID,ISYMN)
     *               + NVIR(ISYMID)*(N-1)
     *               + ID
C
               KOFF2 = KSCR2 - 1
     *               + ICKI(ISYMEN,ISYMI)
     *               + IT1AM(ISYME,ISYMN)
     *               + NVIR(ISYME)*(N-1)
     *               + E
C
               CALL DCOPY(NRHF(ISYMI),XOVVO(KOFF1),NT1AM(ISYMDN),
     *                    WORK(KOFF2),NT1AM(ISYMEN))
C
            ENDDO
         ENDDO
      ENDDO
C
      KOFF1  = KSCR2
     *       + ICKI(ISYMEN,ISYMI)
      KOFF3  = IT1AM(ISYMIB,ISYMI)
     *       + IB
C
      NTOTEN = MAX(1,NT1AM(ISYMEN))
      NTOTB  = MAX(1,NVIR(ISYMIB))
C
      CALL DGEMV('T',NT1AM(ISYMEN),NRHF(ISYMI),ONE,WORK(KOFF1),
     *           NTOTEN,WORK(KSCR1),1,ONE,OMEGA1(KOFF3),NTOTB)
C
C
C--------------------------------
C     Second contribution
C--------------------------------
C
      ISYTMP = MULD2H(ISCKIJ,ISYINT)
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYENF = MULD2H(ISYINT,ISYMI)
      ISYELM = MULD2H(ISYMT2,ISYMID)
C
      KSCR1 = 1
      KEND1 = KSCR1 + NMAIJA(ISYELM)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (OVVO-2)')
      ENDIF
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(I)
     *           + QMAT(INDSQ(I,3))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWORK .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-4)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      DO ISYME = 1, NSYM
         ISYMLM = MULD2H(ISYELM,ISYME)
         DO ISYML = 1, NSYM
            ISYMM  = MULD2H(ISYMLM,ISYML)
            ISYMEL = MULD2H(ISYME,ISYML)
C
            DO L = 1, NRHF(ISYML)
               DO M = 1, NRHF(ISYMM)
C
                  KOFF1 = IT2SP(ISYELM,ISYMID)
     *                  + NCKI(ISYELM)*(ID-1)
     *                  + ICKI(ISYMEL,ISYMM)
     *                  + NT1AM(ISYMEL)*(M-1)
     *                  + IT1AM(ISYME,ISYML)
     *                  + NVIR(ISYME)*(L-1)
     *                  + 1
C
                  KOFF2 = KSCR1 - 1
     *                  + IMAIJA(ISYMLM,ISYME)
     *                  + IMATIJ(ISYMM,ISYML)
     *                  + NRHF(ISYMM)*(L-1)
     *                  + M
C
                  CALL DCOPY(NVIR(ISYME),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMLM))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYME = 1, NSYM
C
         ISYMFN = MULD2H(ISYME,ISYENF)
         ISYMLM = MULD2H(ISCKIJ,ISYMFN)
         ISYFNI = MULD2H(ISYMI,ISYMFN)
         ISYLMI = MULD2H(ISYMI,ISYMLM)
C
         KSCR2 = KEND1
         KEND2 = KSCR2 + NMAIJK(ISYLMI)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (OVVO-3)')
         ENDIF
C
         DO E = 1, NVIR(ISYME)
C
            KOFF1 = ISAIKL(ISYMFN,ISYMLM)
     *            + 1
            KOFF2 = IT2SP(ISYFNI,ISYME)
     *            + NCKI(ISYFNI)*(E-1)
     *            + ICKI(ISYMFN,ISYMI)
     *            + 1
            KOFF3 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
C
            NTOTFN = MAX(1,NT1AM(ISYMFN))
            NTOTLM = MAX(1,NMATIJ(ISYMLM))
C
            CALL DGEMM('T','N',NMATIJ(ISYMLM),NRHF(ISYMI),NT1AM(ISYMFN),
     *                 -ONE,TMAT(KOFF1),NTOTFN,XOVVO(KOFF2),NTOTFN,
     *                 ZERO,WORK(KOFF3),NTOTLM)
C
            KOFF1 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
            KOFF2 = KSCR1
     *            + IMAIJA(ISYMLM,ISYME)
     *            + NMATIJ(ISYMLM)*(E-1)
            KOFF3 = IT1AM(ISYMIB,ISYMI)
     *            + IB
C
            NTOTB  = MAX(1,NVIR(ISYMIB))
            NTOTIJ = MAX(1,NMATIJ(ISYMLM))
C
            CALL DGEMV('T',NMATIJ(ISYMLM),NRHF(ISYMI),ONE,WORK(KOFF1),
     *                 NTOTIJ,WORK(KOFF2),1,ONE,OMEGA1(KOFF3),NTOTB)
C
         ENDDO
      ENDDO
C
C================================================
C     - L^{daf}_{lnm} t^{de}_{lm} g_{infe}
C================================================
C
C--------------------------------
C     First contribution
C--------------------------------
C
      ISYMEN = MULD2H(ISCKIJ,ISYMT2)
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYENI = MULD2H(ISYMEN,ISYMI)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NT1AM(ISYMEN)
      KEND1 = KSCR2 + NCKI(ISYENI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (OOVV-1)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NT1AM(ISYMEN))
      CALL DZERO(WORK(KSCR2),NCKI(ISYENI))
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(I)
     *           + QMAT(INDSQ(I,3))
      ENDDO
C
      DO ISYME = 1, NSYM
         ISYDLM = MULD2H(ISYME,ISYMT2)
         ISYMN  = MULD2H(ISYMEN,ISYME)
C
         KOFF1 = IT2SP(ISYDLM,ISYME)
     *         + 1
         KOFF2 = ISAIKJ(ISYDLM,ISYMN)
     *         + 1
         KOFF3 = KSCR1
     *         + IT1AM(ISYME,ISYMN)
C
         NTODLM = MAX(1,NCKI(ISYDLM))
         NTOTE  = MAX(1,NVIR(ISYME))
C
         CALL DGEMM('T','N',NVIR(ISYME),NRHF(ISYMN),NCKI(ISYDLM),
     *              -ONE,T2TP(KOFF1),NTODLM,TMAT(KOFF2),NTODLM,
     *              ONE,WORK(KOFF3),NTOTE)
C
      ENDDO
C
      DO ISYME = 1, NSYM
         ISYMN  = MULD2H(ISYMEN,ISYME)
         ISYMEI = MULD2H(ISYME,ISYMI)
         ISYMDN = MULD2H(ISYMN,ISYMID)
         ISYMDI = MULD2H(ISYMI,ISYMID)
         ISYDNI = MULD2H(ISYMDN,ISYMI)
C
         DO E = 1, NVIR(ISYME)
            DO N = 1, NRHF(ISYMN)
C
               KOFF1 = IT2SP(ISYDNI,ISYME)
     *               + NCKI(ISYDNI)*(E-1)
     *               + ICKI(ISYMDI,ISYMN)
     *               + NT1AM(ISYMDI)*(N-1)
     *               + IT1AM(ISYMID,ISYMI)
     *               + ID
C
               KOFF2 = KSCR2 - 1
     *               + ICKI(ISYMEN,ISYMI)
     *               + IT1AM(ISYME,ISYMN)
     *               + NVIR(ISYME)*(N-1)
     *               + E
CC
               CALL DCOPY(NRHF(ISYMI),XOOVV(KOFF1),NVIR(ISYMID),
     *                    WORK(KOFF2),NT1AM(ISYMEN))
C
            ENDDO
         ENDDO
      ENDDO
C
      KOFF1  = KSCR2
     *       + ICKI(ISYMEN,ISYMI)
      KOFF3  = IT1AM(ISYMIB,ISYMI)
     *       + IB
C
      NTOTEN = MAX(1,NT1AM(ISYMEN))
      NTOTB  = MAX(1,NVIR(ISYMIB))
C
      CALL DGEMV('T',NT1AM(ISYMEN),NRHF(ISYMI),ONE,WORK(KOFF1),
     *           NTOTEN,WORK(KSCR1),1,ONE,OMEGA1(KOFF3),NTOTB)
C
C
C--------------------------------
C     Second contribution
C--------------------------------
C
      ISYTMP = MULD2H(ISCKIJ,ISYINT)
      ISYMI  = MULD2H(ISYRES,ISYMIB)
      ISYENF = MULD2H(ISYINT,ISYMI)
      ISYELM = MULD2H(ISYMT2,ISYMID)
C
      KSCR1 = 1
      KEND1 = KSCR1 + NMAIJA(ISYELM)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_L3_OMEGA1 (OOVV-2)')
      ENDIF
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,5))
     *           + QMAT(INDSQ(I,4))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWORK .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-5)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      DO ISYME = 1, NSYM
         ISYMLM = MULD2H(ISYELM,ISYME)
         DO ISYML = 1, NSYM
            ISYMM  = MULD2H(ISYMLM,ISYML)
            ISYMEL = MULD2H(ISYME,ISYML)
C
            DO L = 1, NRHF(ISYML)
               DO M = 1, NRHF(ISYMM)
C
                  KOFF1 = IT2SP(ISYELM,ISYMID)
     *                  + NCKI(ISYELM)*(ID-1)
     *                  + ICKI(ISYMEL,ISYMM)
     *                  + NT1AM(ISYMEL)*(M-1)
     *                  + IT1AM(ISYME,ISYML)
     *                  + NVIR(ISYME)*(L-1)
     *                  + 1
C
                  KOFF2 = KSCR1 - 1
     *                  + IMAIJA(ISYMLM,ISYME)
     *                  + IMATIJ(ISYMM,ISYML)
     *                  + NRHF(ISYMM)*(L-1)
     *                  + M
C
                  CALL DCOPY(NVIR(ISYME),T2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMLM))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYME = 1, NSYM
C
         ISYMFN = MULD2H(ISYME,ISYENF)
         ISYMLM = MULD2H(ISCKIJ,ISYMFN)
         ISYFNI = MULD2H(ISYMI,ISYMFN)
         ISYLMI = MULD2H(ISYMI,ISYMLM)
C
         KSCR2 = KEND1
         KSCR3 = KSCR2 + NMAIJK(ISYLMI)
         KEND2 = KSCR3 + NCKI(ISYFNI)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (OOVV-3)')
         ENDIF
C
         DO E = 1, NVIR(ISYME)
C
            DO ISYMF = 1, NSYM
               ISYMN  = MULD2H(ISYMFN,ISYMF)
               ISYMFI = MULD2H(ISYMI,ISYMF)
               DO F = 1, NVIR(ISYMF)
                  DO N = 1, NRHF(ISYMN)
C
                     KOFF1 = IT2SP(ISYFNI,ISYME)
     *                     + NCKI(ISYFNI)*(E-1)
     *                     + ICKI(ISYMFI,ISYMN)
     *                     + NT1AM(ISYMFI)*(N-1)
     *                     + IT1AM(ISYMF,ISYMI)
     *                     + F
C
                     KOFF2 = KSCR3 - 1
     *                     + ICKI(ISYMFN,ISYMI)
     *                     + IT1AM(ISYMF,ISYMN)
     *                     + NVIR(ISYMF)*(N-1)
     *                     + F
C
                     CALL DCOPY(NRHF(ISYMI),XOOVV(KOFF1),NVIR(ISYMF),
     *                          WORK(KOFF2),NT1AM(ISYMFN))
                  ENDDO
               ENDDO
            ENDDO
C
            KOFF1 = ISAIKL(ISYMFN,ISYMLM)
     *            + 1
            KOFF2 = KSCR3
     *            + ICKI(ISYMFN,ISYMI)
            KOFF3 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
C
            NTOTFN = MAX(1,NT1AM(ISYMFN))
            NTOTLM = MAX(1,NMATIJ(ISYMLM))
C
            CALL DGEMM('T','N',NMATIJ(ISYMLM),NRHF(ISYMI),NT1AM(ISYMFN),
     *                 -ONE,TMAT(KOFF1),NTOTFN,WORK(KOFF2),NTOTFN,
     *                 ZERO,WORK(KOFF3),NTOTLM)
C
            KOFF1 = KSCR2
     *            + IMAIJK(ISYMLM,ISYMI)
            KOFF2 = KSCR1
     *            + IMAIJA(ISYMLM,ISYME)
     *            + NMATIJ(ISYMLM)*(E-1)
            KOFF3 = IT1AM(ISYMIB,ISYMI)
     *            + IB
C
            NTOTB  = MAX(1,NVIR(ISYMIB))
            NTOTIJ = MAX(1,NMATIJ(ISYMLM))
C
            CALL DGEMV('T',NMATIJ(ISYMLM),NRHF(ISYMI),ONE,WORK(KOFF1),
     *                 NTOTIJ,WORK(KOFF2),1,ONE,OMEGA1(KOFF3),NTOTB)
C
         ENDDO
      ENDDO
C
C================================================
C     - L^{def}_{lin} t^{de}_{lm} g_{mafn}
C================================================
C
      ISYDLM = MULD2H(ISYMT2,ISYMID)
      ISYTMP = MULD2H(ISYDLM,ISCKIJ)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NMAIJK(ISYTMP)
      KEND1 = KSCR2 + NCKI(ISYDLM)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_L3_OMEGA1 (OVVO-4)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NMAIJK(ISYTMP))
C
C-----------------------------------------------
C     First contribution to intermediate
C-----------------------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(I)
     *           + QMAT(INDSQ(I,3))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-6)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
C
      DO ISYMDL = 1, NSYM
         ISYMM  = MULD2H(ISYDLM,ISYMDL)
         ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
         KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
         KOFF2 = IT2SP(ISYDLM,ISYMID)
     *         + NCKI(ISYDLM)*(ID-1)
     *         + ICKI(ISYMDL,ISYMM)
     *         + 1
         KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
         NTOTDL = MAX(1,NT1AM(ISYMDL))
         NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,T2TP(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
      ENDDO
C
C-----------------------------------------------
C     Second contribution to intermediate
C-----------------------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,1))
     *           + QMAT(INDSQ(I,2))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-6)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      DO ISYMM = 1, NSYM
         ISYMDL = MULD2H(ISYDLM,ISYMM)
         DO ISYMD = 1, NSYM
            ISYML  = MULD2H(ISYMDL,ISYMD)
            ISYMDM = MULD2H(ISYMD,ISYMM)
            DO M = 1, NRHF(ISYMM)
               DO L = 1, NRHF(ISYML)
C
                  KOFF1 = IT2SP(ISYDLM,ISYMID)
     *                  + NCKI(ISYDLM)*(ID-1)
     *                  + ICKI(ISYMDL,ISYMM)
     *                  + NT1AM(ISYMDL)*(M-1)
     *                  + IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L-1)
     *                  + 1
C
                  KOFF2 = KSCR2
     *                  + ICKI(ISYMDM,ISYML)
     *                  + NT1AM(ISYMDM)*(L-1)
     *                  + IT1AM(ISYMD,ISYMM)
     *                  + NVIR(ISYMD)*(M-1)
C
                  CALL DCOPY(NVIR(ISYMD),T2TP(KOFF1),1,
     *                       WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYMDL = 1, NSYM
         ISYMM  = MULD2H(ISYDLM,ISYMDL)
         ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
         KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
         KOFF2 = KSCR2
     *         + ICKI(ISYMDL,ISYMM)
         KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
         NTOTDL = MAX(1,NT1AM(ISYMDL))
         NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,WORK(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
      ENDDO
C
C-------------------------------------------------------
C     Sort intermediate and integrals and contract
C-------------------------------------------------------
C
      KEND1 = KSCR2 + NMAIJK(ISYTMP)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_L3_OMEGA1 (OVVO-5)')
      ENDIF
C
      DO ISYMM = 1, NSYM
         ISYMIN = MULD2H(ISYTMP,ISYMM)
         DO ISYMN = 1, NSYM
            ISYMI  = MULD2H(ISYMIN,ISYMN)
            ISYMMN = MULD2H(ISYMM,ISYMN)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = KSCR1
     *                  + IMAIJK(ISYMIN,ISYMM)
     *                  + NMATIJ(ISYMIN)*(M-1)
     *                  + IMATIJ(ISYMI,ISYMN)
     *                  + NRHF(ISYMI)*(N-1)
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJK(ISYMMN,ISYMI)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M
C
                  CALL DCOPY(NRHF(ISYMI),WORK(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL DCOPY(NMAIJK(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
      ISYAMN = MULD2H(ISYINT,ISYMIB)
C
      DO ISYMA = 1, NSYM
         ISYMMN = MULD2H(ISYMA,ISYAMN)
         ISYBMN = MULD2H(ISYINT,ISYMA)
         DO ISYMM = 1, NSYM
            ISYMN  = MULD2H(ISYMMN,ISYMM)
            ISYMBN = MULD2H(ISYBMN,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = IT2SP(ISYBMN,ISYMA)
     *                  + ICKI(ISYMBN,ISYMM)
     *                  + NT1AM(ISYMBN)*(M-1)
     *                  + IT1AM(ISYMIB,ISYMN)
     *                  + NVIR(ISYMIB)*(N-1)
     *                  + IB
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJA(ISYMMN,ISYMA)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M
C
                  CALL DCOPY(NVIR(ISYMA),XOVVO(KOFF1),NCKI(ISYBMN),
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMMN = MULD2H(ISYAMN,ISYMA)
C
         KOFF1 = KSCR2
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KSCR1
     *         + IMAIJK(ISYMMN,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NTOTMN = MAX(1,NMATIJ(ISYMMN))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              ONE,WORK(KOFF1),NTOTMN,WORK(KOFF2),NTOTMN,
     *              ONE,OMEGA1(KOFF3),NTOTA)
      ENDDO
C
C================================================
C     - L^{def}_{lni} t^{de}_{lm} g_{mnfa}
C================================================
C
      ISYDLM = MULD2H(ISYMT2,ISYMID)
      ISYTMP = MULD2H(ISYDLM,ISCKIJ)
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NMAIJK(ISYTMP)
      KEND1 = KSCR2 + NCKI(ISYDLM)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_L3_OMEGA1 (OVVO-4)')
      ENDIF
C
      CALL DZERO(WORK(KSCR1),NMAIJK(ISYTMP))
C
C-----------------------------------------------
C     First contribution to intermediate
C-----------------------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,3))
     *           + QMAT(I)
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-6)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
C
      DO ISYMDL = 1, NSYM
         ISYMM  = MULD2H(ISYDLM,ISYMDL)
         ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
         KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
         KOFF2 = IT2SP(ISYDLM,ISYMID)
     *         + NCKI(ISYDLM)*(ID-1)
     *         + ICKI(ISYMDL,ISYMM)
     *         + 1
         KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
         NTOTDL = MAX(1,NT1AM(ISYMDL))
         NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,T2TP(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
      ENDDO
C
C-----------------------------------------------
C     Second contribution to intermediate
C-----------------------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = SMAT(INDSQ(I,4))
     *           + QMAT(INDSQ(I,5))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CC3_L3_OMEGA1 (CC_GATHER-6)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      DO ISYMM = 1, NSYM
         ISYMDL = MULD2H(ISYDLM,ISYMM)
         DO ISYMD = 1, NSYM
            ISYML  = MULD2H(ISYMDL,ISYMD)
            ISYMDM = MULD2H(ISYMD,ISYMM)
            DO M = 1, NRHF(ISYMM)
               DO L = 1, NRHF(ISYML)
C
                  KOFF1 = IT2SP(ISYDLM,ISYMID)
     *                  + NCKI(ISYDLM)*(ID-1)
     *                  + ICKI(ISYMDL,ISYMM)
     *                  + NT1AM(ISYMDL)*(M-1)
     *                  + IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L-1)
     *                  + 1
C
                  KOFF2 = KSCR2
     *                  + ICKI(ISYMDM,ISYML)
     *                  + NT1AM(ISYMDM)*(L-1)
     *                  + IT1AM(ISYMD,ISYMM)
     *                  + NVIR(ISYMD)*(M-1)
C
                  CALL DCOPY(NVIR(ISYMD),T2TP(KOFF1),1,
     *                       WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYMDL = 1, NSYM
         ISYMM  = MULD2H(ISYDLM,ISYMDL)
         ISYMIN = MULD2H(ISCKIJ,ISYMDL)
C
         KOFF1 = ISAIKL(ISYMDL,ISYMIN)
     *         + 1
         KOFF2 = KSCR2
     *         + ICKI(ISYMDL,ISYMM)
         KOFF3 = KSCR1
     *         + IMAIJK(ISYMIN,ISYMM)
C
         NTOTDL = MAX(1,NT1AM(ISYMDL))
         NTOTIN = MAX(1,NMATIJ(ISYMIN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMIN),NRHF(ISYMM),NT1AM(ISYMDL),
     *              -ONE,TMAT(KOFF1),NTOTDL,WORK(KOFF2),NTOTDL,
     *              ONE,WORK(KOFF3),NTOTIN)
      ENDDO
C
C-------------------------------------------------------
C     Sort intermediate and integrals and contract
C-------------------------------------------------------
C
      KEND1 = KSCR2 + NMAIJK(ISYTMP)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory CC3_L3_OMEGA1 (OVVO-5)')
      ENDIF
C
      DO ISYMM = 1, NSYM
         ISYMIN = MULD2H(ISYTMP,ISYMM)
         DO ISYMN = 1, NSYM
            ISYMI  = MULD2H(ISYMIN,ISYMN)
            ISYMMN = MULD2H(ISYMM,ISYMN)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = KSCR1
     *                  + IMAIJK(ISYMIN,ISYMM)
     *                  + NMATIJ(ISYMIN)*(M-1)
     *                  + IMATIJ(ISYMI,ISYMN)
     *                  + NRHF(ISYMI)*(N-1)
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJK(ISYMMN,ISYMI)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M
C
                  CALL DCOPY(NRHF(ISYMI),WORK(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL DCOPY(NMAIJK(ISYTMP),WORK(KSCR2),1,WORK(KSCR1),1)
C
      ISYAMN = MULD2H(ISYINT,ISYMIB)
C
      DO ISYMA = 1, NSYM
         ISYMMN = MULD2H(ISYMA,ISYAMN)
         ISYBMN = MULD2H(ISYINT,ISYMA)
         DO ISYMM = 1, NSYM
            ISYMN  = MULD2H(ISYMMN,ISYMM)
            ISYMBN = MULD2H(ISYBMN,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = IT2SP(ISYBMN,ISYMA)
     *                  + ICKI(ISYMBN,ISYMM)
     *                  + NT1AM(ISYMBN)*(M-1)
     *                  + IT1AM(ISYMIB,ISYMN)
     *                  + NVIR(ISYMIB)*(N-1)
     *                  + IB
C
                  KOFF2 = KSCR2 - 1
     *                  + IMAIJA(ISYMMN,ISYMA)
     *                  + IMATIJ(ISYMN,ISYMM)
     *                  + NRHF(ISYMN)*(M-1)
     *                  + N
C
                  CALL DCOPY(NVIR(ISYMA),XOOVV(KOFF1),NCKI(ISYBMN),
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMMN = MULD2H(ISYAMN,ISYMA)
C
         KOFF1 = KSCR2
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KSCR1
     *         + IMAIJK(ISYMMN,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NTOTMN = MAX(1,NMATIJ(ISYMMN))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              ONE,WORK(KOFF1),NTOTMN,WORK(KOFF2),NTOTMN,
     *              ONE,OMEGA1(KOFF3),NTOTA)
      ENDDO
C
C-------------------
C     End.
C-------------------
C
      CALL QEXIT('CC3_L3_OMEGA1')
C
      RETURN
      END
C  /* Deck cc3_sort4o */
      SUBROUTINE CC3_SORT4O(XOOOO,ISYINT,WORK,LWORK)
C
C     Written by K. Hald, Spring 2002.
C
C     Sort (i,j,k,l) to (i,j,l,k)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYINT, LWORK
      INTEGER ISYML, ISYIJK, ISYMK, ISYMIJ, ISYIJL, KOFF1, KOFF2
      INTEGER ISYMJ, ISYMI, ISYMLJ, ISYMIL, KSCR1, KEND1, LWRK1
      INTEGER NJUMP
C
#if defined (SYS_CRAY)
      REAL XOOOO(*), WORK(LWORK)
#else
      DOUBLE PRECISION XOOOO(*), WORK(LWORK)
#endif
C
C
      CALL QENTER('CC3_SORT4O')
C
      KSCR1 = 1
      KEND1 = KSCR1 + N3ORHF(ISYINT)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_SORT4O')
      ENDIF
C
C--------------
C     Sort.
C--------------
C
      DO ISYML = 1, NSYM
         ISYIJK = MULD2H(ISYINT,ISYML)
         DO ISYMK = 1, NSYM
            ISYMIJ = MULD2H(ISYIJK,ISYMK)
            ISYIJL = MULD2H(ISYMIJ,ISYML)
            DO ISYMJ = 1, NSYM
               ISYMI  = MULD2H(ISYMIJ,ISYMJ)
               ISYMIL = MULD2H(ISYMI,ISYML)
               ISYMLJ = MULD2H(ISYML,ISYMJ)
               DO L = 1, NRHF(ISYML)
                  DO K = 1, NRHF(ISYMK)
                     DO J = 1, NRHF(ISYMJ)
C
                        KOFF1 = I3ORHF(ISYIJK,ISYML)
     *                        + NMAIJK(ISYIJK)*(L-1)
     *                        + IMAIJK(ISYMIJ,ISYMK)
     *                        + NMATIJ(ISYMIJ)*(K-1)
     *                        + IMATIJ(ISYMI,ISYMJ)
     *                        + NRHF(ISYMI)*(J-1)
     *                        + 1
C
C                        NJUMP = NMATIJ(ISYMLJ)
C                        KOFF2 = KSCR1 - 1
C     *                        + I3ORHF(ISYIJL,ISYMK)
C     *                        + NMAIJK(ISYIJL)*(K-1)
C     *                        + IMAIJK(ISYMLJ,ISYMI)
C     *                        + IMATIJ(ISYML,ISYMJ)
C     *                        + NRHF(ISYML)*(J-1)
C     *                        + L
C
                        NJUMP = NMATIJ(ISYMLJ)
                        KOFF2 = KSCR1 - 1
     *                        + I3ORHF(ISYIJL,ISYMK)
     *                        + NMAIJK(ISYIJL)*(K-1)
     *                        + IMAIJK(ISYMLJ,ISYMI)
     *                        + IMATIJ(ISYMJ,ISYML)
     *                        + NRHF(ISYMJ)*(L-1)
     *                        + J
C
C                        NJUMP = 1
C                        KOFF2 = KSCR1 - 1
C     *                        + I3ORHF(ISYIJL,ISYMK)
C     *                        + NMAIJK(ISYIJL)*(K-1)
C     *                        + IMAIJK(ISYMIL,ISYMJ)
C     *                        + NMATIJ(ISYMIL)*(J-1)
C     *                        + IMATIJ(ISYMI,ISYML)
C     *                        + NRHF(ISYMI)*(L-1)
C     *                        + 1
C
C                        NJUMP = 1
C                        KOFF2 = KSCR1 - 1
C     *                        + I3ORHF(ISYIJL,ISYMK)
C     *                        + NMAIJK(ISYIJL)*(K-1)
C     *                        + IMAIJK(ISYMIJ,ISYML)
C     *                        + NMATIJ(ISYMIJ)*(L-1)
C     *                        + IMATIJ(ISYMI,ISYMJ)
C     *                        + NRHF(ISYMI)*(J-1)
C     *                        + 1
C
C                        NJUMP = NRHF(ISYMJ)
C                        KOFF2 = KSCR1 - 1
C     *                        + I3ORHF(ISYIJL,ISYMK)
C     *                        + NMAIJK(ISYIJL)*(K-1)
C     *                        + IMAIJK(ISYMIJ,ISYML)
C     *                        + NMATIJ(ISYMIJ)*(L-1)
C     *                        + IMATIJ(ISYMJ,ISYMI)
C     *                        + J
C
C                        NJUMP = NRHF(ISYMJ)
C                        KOFF2 = KSCR1 - 1
C     *                        + I3ORHF(ISYIJL,ISYMK)
C     *                        + NMAIJK(ISYIJL)*(K-1)
C     *                        + IMAIJK(ISYMIL,ISYMJ)
C     *                        + NMATIJ(ISYMIL)*(J-1)
C     *                        + IMATIJ(ISYML,ISYMI)
C     *                        + L
C
C
                        CALL DCOPY(NRHF(ISYMI),XOOOO(KOFF1),1,
     *                             WORK(KOFF2),NJUMP)
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C----------------------------------------------
C     Copy sorted integrals back into XOOOO
C----------------------------------------------
C
      CALL DCOPY(N3ORHF(ISYINT),WORK(KSCR1),1,XOOOO,1)
C
C-------------------
C     End.
C-------------------
C
      CALL QEXIT('CC3_SORT4O')
C
      RETURN
      END
C  /* Deck n1_gv4 */
      SUBROUTINE N1_GV4(IOPT,!0 - normal transformation; 1 - one-index trans.
     *                  LUN1,FNN1,      !  file containing (fge,i)
     *                  ISYMN1,         !  symmetry of N1(fge,i)
     *                  XLAMP1,ISLAMP1, !  LAM^{p}_{al f} 
     *                  XLAMP2,ISLAMP2, !  LAM^{p}_{gam g}
     *                  XLAMH1,ISLAMH1, !  LAM^{h}_{del e}
     *                  XLAMH2,ISLAMH2, !  LAM^{h}_{be a}
     *                  OMEGA1,ISOMEGA1,!  result vector
     *                  WORK,LWORK)
*
**********************************************************************
*                                                                    *
* Calculate that part of <T3|[[H,T2],tau_ai]|HF> contribution        *
* to the left transformation which involves VVVV integrals:          *
*                                                                    *
*    omega_ai = omega_ai + sum_efg g_fage N_fige                     *
*                                                                    *
*                                                                    *
* To avoid an expensive construction of g_fage integrals the term    *
* is calculated AO-direct (al=alpha, be=beta, ga=gamma, de=delta):   *
*                                                                    *
* sum_efg g_fage N_fige =                                            *
* sum_be LAM^{h}_{be a} sum_{al ga de} g_{al be ga de} N_{al i ga de}*
*                                                                    *
* where N_{al i ga de} is obtained by backtransformation of N_fige:  *
*                                                                    *
*  N_{al i ga de}                                                    *
*     = sum_efg LAM^{p}_{al f} LAM^{p}_{ga g} LAM^{h}_{de e} N_fige  *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* STRUCTURE OF THE ROUTINE:                                          *
* =========================                                          *
*                                                                    *
* 1) The general structure involves the loop over AO integrals       *
*    distributions (AO-direct approach) and is taken from            *
*    CC3_INT routine.                                                *
*    Thus "delta" is a general loop index in this routine.           *
*                                                                    *
* 2) Once a distribution of AO integrals is available for            *
*    a fixed "delta", we backtransform N_fige ("delta" still fixed): *
*                                                                    *
*         N_fige -->  N_{al i ga de} (see the last equation above)   *
*                                                                    *
*    This is done by the call to N1_2AORHF routine.                  * 
*                                                                    *
* 3) Next, in N12AO_G3AO routine, a distribution of integrals        *
*    and N_{al i ga de} are read in for a fixed "delta"; they        *
*    are contracted together:                                        *
*                                                                    *
*         g_{al be ga de}  * N_{al i ga de} --> I_{beta,i}           *
*                                                                    *
*    and finally "beta" index is transformed to "a" MO index,        *
*    such that the result vector, omega_ai, is obtained.             *
*                                                                    *
*====================================================================*
* NB !!!                                                             *
* ======                                                             *
*                                                                    *
*  N_fige intermediate is actually sitting on file as (fge,i)        *
*  and is read in for each "I" as (fg,e).                            *
*                                                                    *
*====================================================================*
*                                                                    *
* IOPT:                                                              *
* =====                                                              *
*                                                                    *
*  IOPT = 0 : Use normal MO Lambda matrices to backtransform, i.e.:  *
*                                                                    *
*         N_fige -->  N_{al i ga de}                                 *
*                                                                    *
*         N1_2AORHF routine is thus called once                      *
*         and XLAMP1 and XLAMP2 matrices are the same.               *
*                                                                    *
*                                                                    *
*  IOPT = 1 : Use normal T1-Lambda matrices to backtransform, i.e.:  *
*                                                                    *
*         N_fige -->  N_{al- i ga de} + N_{al i ga- de}              *
*                                                                    *
*         ("-" means one-index backtransformed AO index)             *
*                                                                    *
*         N1_2AORHF routine is thus called twice, the second time    *
*         with XLAMP1 and XLAMP2 matrices interchanged               *
*         (XLAMP1 and XLAMP2 matrices are NOT the same anymore).     *
*                                                                    *
*  NOTE that with IOPT = 1 also double one-index backtransformation: *
*  ****                                                              *
*                                                                    *
*         N_fige -->  N_{al- i ga~ de} + N_{al~ i ga- de}            *
*                                                                    *
*       may be carried out ("~" is in general a transformation       *
*       different from "-").                                         *
*       You decide about that by appropriate choice of XLAMP1 and    *
*       XLAMP2, which are passed to this routine.                    *
*                                                                    *
**********************************************************************
*  F. Pawlowski, 17-Mar-2004, Aarhus.                                *
**********************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h> 
#include <ccisao.h>
#include <r12int.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
C
      CHARACTER*(*) FNN1
      INTEGER       LUN1
C
      INTEGER ISYMN1,ISLAMP1,ISLAMP2,ISLAMH1,ISLAMH2,ISOMEGA1,LWORK
      INTEGER KCCFB1,KINDXB,KEND1,LWRK1
      INTEGER KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2
      INTEGER KODPP1,KODPP2,KRDPP1,KRDPP2,KFREE,LFREE
      INTEGER NTOSYM,KENDSV,LWRKSV,ISYMD1,NTOT,ILLL
      INTEGER KRECNR,NUMDIS,IDEL2,IDEL,ISYMD
      INTEGER ISTMP,ISN12AORHF,KXINT,KN1MAT2AORHF,KEND2,LWRK2
      INTEGER INDEXA(MXCORB),ISYDIS
C
      INTEGER IOPT
c
      integer isymi,isymalbe
      INTEGER IDUM
C
#if defined (SYS_CRAY)
      REAL XLAMP1(*),XLAMP2(*),XLAMH1(*),XLAMH2(*)
      REAL OMEGA1(*),WORK(LWORK)
      real ddot
#else
      DOUBLE PRECISION XLAMP1(*),XLAMP2(*),XLAMH1(*),XLAMH2(*)
      DOUBLE PRECISION OMEGA1(*),WORK(LWORK)
      double precision ddot
#endif
C
      CALL QENTER('NGV4')
C
      !Initial check of IOPT
      IF (IOPT .GT. 1) THEN
         WRITE(LUPRI,*) 'IOPT = ',IOPT
         CALL QUIT('Illegal value of IOPT in N1_GV4')
      END IF
C
      !Because of direct/non-direct logical switches we need to
      !initialise kend1 
      KEND1 = 1
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
            LWRK1 = LWORK - KEND1
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = 1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available: ', LWORK
               WRITE(LUPRI,*) 'Memory needed   : ', KEND1
               CALL QUIT('Insufficient memory in N1_GV4 (1)')
            END IF
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available: ', KEND1 + LWRK1
               WRITE(LUPRI,*) 'Memory needed   : ', KEND1
               CALL QUIT('Insufficient memory in N1_GV4 (2)')
            END IF
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      DO ISYMD1 = 1,NTOSYM
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
         DO ILLL = 1,NTOT
C
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed CTR2.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRINT)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                        WORK(KODCL1),WORK(KODCL2),
     *                        WORK(KODBC1),WORK(KODBC2),
     *                        WORK(KRDBC1),WORK(KRDBC2),
     *                        WORK(KODPP1),WORK(KODPP2),
     *                        WORK(KRDPP1),WORK(KRDPP2),
     *                        WORK(KCCFB1),WORK(KINDXB),
     *                        WORK(KEND1),LWRK1,IPRERI)
               ENDIF
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available: ', LWORK
                  WRITE(LUPRI,*) 'Memory needed   : ', KEND1
                  CALL QUIT('Insufficient memory in N1_GV4 (3)')
               END IF

C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO IDEL2 = 1,NUMDIS
C
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
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)

               !determine symmetry of N_{alpha i gamma delta} (delta is fixed!)
               ISTMP = MULD2H(ISYMN1,ISLAMP1)
               ISTMP = MULD2H(ISTMP,ISLAMP2)
               ISTMP = MULD2H(ISTMP,ISLAMH1)
               ISN12AORHF = MULD2H(ISTMP,ISYMD)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT         = KEND1
               KN1MAT2AORHF  = KXINT + NDISAO(ISYDIS)
               KEND2         = KN1MAT2AORHF + NDSRHFSQ(ISN12AORHF)
               LWRK2         = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available: ', LWORK
                  WRITE(LUPRI,*) 'Memory needed   : ', KEND2
                  CALL QUIT('Insufficient memory in N1_GV4 (4)')
               END IF
C
               CALL DZERO(WORK(KN1MAT2AORHF),NDSRHFSQ(ISN12AORHF))
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C-------------------------------------------------------------------------------C              Backtransform N_fige to N_{alpha i gamma delta} (delta is fixed!)C-------------------------------------------------------------------------------C
               CALL N1_2AORHF(ISYMN1,
     *                        LUN1,FNN1,
     *                        XLAMP1,ISLAMP1,
     *                        XLAMP2,ISLAMP2,
     *                        XLAMH1,ISLAMH1,
     *                        WORK(KN1MAT2AORHF),ISN12AORHF,
     *                        ISYMD,IDEL,   
     *                        WORK(KEND2),LWRK2,
     *                        IOPT)
C
C------------------------------------------------------------------
C              Contract N^{del}_{al,i,gam} with g^{del}_{al,be,gam}
C              and transform "beta" index to "a" MO index.
C------------------------------------------------------------------
C

               CALL N12AO_G3AO(WORK(KN1MAT2AORHF),ISN12AORHF,
     *                         WORK(KXINT),ISYMOP,
     *                         XLAMH2,ISLAMH2,OMEGA1,ISOMEGA1,
     *                         IDEL,ISYMD,WORK(KEND2),LWRK2)
C              
C---------------------------------------
C     END ALL LOOPS
C---------------------------------------
C
            ENDDO    ! IDEL2
         ENDDO       ! ILLL
      ENDDO          ! ISYMD1
C
C----------
C     End.
C----------
C
      CALL QEXIT('NGV4')
C
      RETURN
      END
C  /* Deck n1_2aorhf */
      SUBROUTINE N1_2AORHF(ISYMN1,        !  symmetry of N
     *                     LUN1,FNN1,     !  file containing N(fge,i)
     *                     XLAMP1,ISLAMP1,!  LAM^{p}_{al f} 
     *                     XLAMP2,ISLAMP2,!  LAM^{p}_{gam g}
     *                     XLAMH1,ISLAMH1,!  LAM^{h}_{del e}
     *                     N1MAT2AORHF,ISN12AORHF, !  N_{alpha,i,gamma}
     *                     ISYDEL,IDEL,   !  delta: fixed AO index
     *                     WORK,LWORK,
     *                     IOPT) !0 - normal transformation; 1 - one-index trans
*
************************************************************************
*     Backtransform N_fige to N_{alpha i gamma delta} (delta is fixed!)
*
*     This is just a little driver; the real stuff is done in 
*     N1_2AORHF_1 routine, which is called from here.
*
*     IOPT = 0 : do normal back transformation with MO Lambda matrices
*     IOPT = 1 : do (double) one-index back transformation with T1 matrices
*
*     FP, 28-Apr-2004, Aarhus.
************************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CHARACTER*(*) FNN1
      INTEGER       LUN1
C
      INTEGER ISYMN1,ISLAMP1,ISLAMP2,ISLAMH1,ISN12AORHF,ISYDEL,IDEL
      INTEGER IOPT,LWORK
      INTEGER ISYMI,ISFGE,KFGE,KEND1,LWRK1,IADR
C
#if defined (SYS_CRAY)
      REAL XLAMP1(*),XLAMP2(*),XLAMH1(*)
      REAL N1MAT2AORHF(*),WORK(LWORK)
#else
      DOUBLE PRECISION XLAMP1(*),XLAMP2(*),XLAMH1(*)
      DOUBLE PRECISION N1MAT2AORHF(*),WORK(LWORK)
#endif
C
      CALL QENTER('N1AO')
C
      !Initial check of IOPT
      IF ((IOPT.GT.1) .OR. (IOPT.LT.0)) THEN
         WRITE(LUPRI,*) 'IOPT = ',IOPT
         CALL QUIT('Illegal value of IOPT in N1_2AORHF')
      END IF
C
      DO ISYMI = 1,NSYM
C
          ISFGE = MULD2H(ISYMN1,ISYMI)
C
          KFGE   = 1 !an array used to read in N1MAT(fge,I) from file
          KEND1  = KFGE + NMAABC(ISFGE)
          LWRK1  = LWORK  - KEND1
C
          IF (LWRK1 .LT. 0) THEN
             WRITE(LUPRI,*) 'Memory available: ', LWORK
             WRITE(LUPRI,*) 'Memory needed   : ', KEND1
             CALL QUIT('Insufficient memory in N1_2AORHF (xx)')
          END IF
C
          DO I = 1,NRHF(ISYMI)
C
            CALL DZERO(WORK(KFGE),NMAABC(ISFGE))
C
            IF (NMAABC(ISFGE).GT.0) THEN
               !Read in KFGE(fg,e) array from file
               IADR = IMAABCI(ISFGE,ISYMI) + NMAABC(ISFGE)*(I-1) + 1 
               CALL GETWA2(LUN1,FNN1,WORK(KFGE),IADR,NMAABC(ISFGE))
            END IF
C     
            !Backtransform N_fige to N_{alpha i gamma delta}
            CALL N1_2AORHF_1(WORK(KFGE),ISYMN1,XLAMP1,ISLAMP1,
     *                     XLAMP2,ISLAMP2,XLAMH1,ISLAMH1,
     *                     N1MAT2AORHF,ISN12AORHF,
     *                     ISYDEL,IDEL,ISYMI,I,
     *                     WORK(KEND1),LWRK1)
C
            IF (IOPT.EQ.1) THEN !do one-index backtransformation 
                                !with T1.  Note that XLAMP1 and XLAMP2 
                                !are interchanged
C
               !Backtransform N_fige to N_{alpha i gamma delta}
               CALL N1_2AORHF_1(WORK(KFGE),ISYMN1,XLAMP2,ISLAMP2,
     *                        XLAMP1,ISLAMP1,XLAMH1,ISLAMH1,
     *                        N1MAT2AORHF,ISN12AORHF,
     *                        ISYDEL,IDEL,ISYMI,I,
     *                        WORK(KEND1),LWRK1)
            END IF
C
         END DO !I
      END DO !ISYMI
C
C----------
C     End.
C----------
C
      CALL QEXIT('N1AO')
C
      RETURN
      END
C  /* Deck n1_2aorhf_1 */
      SUBROUTINE N1_2AORHF_1(N1FGE,         !an array containing (fge,I)
     *                       ISYMN1,        !  symmetry of N1(fge,i)
     *                       XLAMP1,ISLAMP1,!  LAM^{p}_{al f} 
     *                       XLAMP2,ISLAMP2,!  LAM^{p}_{gam g}
     *                       XLAMH1,ISLAMH1,!  LAM^{h}_{del e}
     *                       N1MAT2AORHF,ISN12AORHF, !  N_{alpha,i,gamma}
     *                       ISYDEL,IDEL,   !  delta: fixed AO index
     *                       ISYMI,I,       !fixed occupied "I" index   
     *                       WORK,LWORK)
*
**********************************************************************
*                                                                    *
* Backtransform N_fige intermmediate (which has been constructed in  *
* WT2_N1N2 routine) to N_{alpha i gamma delta}, where greek indeces  *
* denote AO indeces:                                                 *
*                                                                    *
* N_{al i gam del}                                                   *
*   =  sum_fge LAM^{p}_{al f} LAM^{p}_{gam g} LAM^{h}_{del e} N_fige *
*                                                                    *
* Calculations are carried out for the FIXED DELTA AO index.         *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* N_{al i gam del} is stored as:                                     *
* ==============================                                     *
*                                                                    *
* N1MAT2AORHF(al gam,i)   !del is fixed                              *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* NOTE: Since "i" MO index is not backtransformed, "i" is our GENERAL*
* ====  LOOP index in this routine and it comes from OUTSIDE.        *
*                                                                    *
*       N_fige is actually read in from LUN1 file (outside), where it*
*       is sitting as (fge,i) (remember that "i" is fixed from outside)*
*                                                                    *
**********************************************************************
*  F. Pawlowski, 16-Mar-2004, Aarhus.                                *
**********************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMN1,ISLAMP1,ISLAMP2,ISLAMH1,ISN12AORHF,IDEL,LWORK
      INTEGER ISTEMP
      INTEGER ISYME,ISFGI,ISYMI,ISFG,ISFGE,ISFGDEL,KFGDEL,KEND1,LWRK1
      INTEGER KOFF1,KOFF2,KOFF3,NFG,NDEL
      INTEGER ISYMG,ISYMF,ISYGAM,ISYAL,ISFGAM,ISALGAM,KFGAM,KEND2,LWRK2
      INTEGER NF,NGAM,NAL,ISYDEL
C
      INTEGER IADR
C
#if defined (SYS_CRAY)
      REAL XLAMP1(*),XLAMP2(*),XLAMH1(*)
      REAL N1MAT2AORHF(*),WORK(LWORK)
      REAL N1FGE(*)
      REAL ONE
      real ddot
#else
      DOUBLE PRECISION XLAMP1(*),XLAMP2(*),XLAMH1(*)
      DOUBLE PRECISION N1MAT2AORHF(*),WORK(LWORK)
      DOUBLE PRECISION N1FGE(*)
      DOUBLE PRECISION ONE
      double precision ddot
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('N1AO_1')
C
      !Initial symmetry check
      ISTEMP = MULD2H(ISLAMP1,ISLAMP2)
      ISTEMP = MULD2H(ISLAMH1,ISTEMP)
      ISTEMP = MULD2H(ISYMN1,ISTEMP)
      ISTEMP = MULD2H(ISYDEL,ISTEMP)
      IF (ISTEMP .NE. ISN12AORHF) THEN
         WRITE(LUPRI,*) 'The following 5 symmetries: '
         WRITE(LUPRI,*) 'ISLAMP1 = ', ISLAMP1
         WRITE(LUPRI,*) 'ISLAMP2 = ', ISLAMP2
         WRITE(LUPRI,*) 'ISLAMH1 = ', ISLAMH1
         WRITE(LUPRI,*) 'ISYMN1  = ', ISYMN1
         WRITE(LUPRI,*) 'ISYDEL  = ', ISYDEL
         WRITE(LUPRI,*) 'when multiplied together should give: '
         WRITE(LUPRI,*) 'ISN12AORHF = ', ISN12AORHF
         CALL QUIT('Symmetry mismatch in N1_2AORHF')
      END IF
C
C     -------------------------------------------------------
C      Calculate the contraction LAM^{h}_{del e} N_fige:     
C
C       sum_fge LAM^{h}_{del e} N_fige --> N_{fig del}
C
C      calculated as:
C
C       N1MAT(fge,i)  *   XLAMH1(del e) !N1MAT(fge,i) sitting on file
C            |           
C            |loop over i
C            V           
C       N1FGE(fg,e)    *   XLAMH1(del e) --> KFGDEL(fg,del)
C
C     -------------------------------------------------------
C
      ISYME = MULD2H(ISLAMH1,ISYDEL) !because: XLAMH1(del e)
      ISFGI = MULD2H(ISYMN1,ISYME)   !because: N_fige
C
      ISFG    = MULD2H(ISFGI,ISYMI)
      ISFGE   = MULD2H(ISFG,ISYME)
      ISFGDEL = MULD2H(ISFG,ISYDEL) 
C
      KFGDEL = 1
      KEND1  = KFGDEL + NMATAB(ISFG)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in N1_2AORHF (1)')
      END IF
C
      CALL DZERO(WORK(KFGDEL),NMATAB(ISFG))
C
      !multiply N1FGE(fg,e)   *   XLAMH1(del e) --> KFGDEL(fg,del)
      KOFF1 = IMAABC(ISFG,ISYME)
     *      + 1
      KOFF2 = IGLMVI(ISYDEL,ISYME)
     *      + IDEL-IBAS(ISYDEL)
      KOFF3 = KFGDEL
C
      NFG  = MAX(NMATAB(ISFG),1)
*     NDEL = MAX(NBAS(ISYDEL),1)
      NDEL = NBAS(ISYDEL)
C
      CALL DGEMV('N',NMATAB(ISFG),NVIR(ISYME),ONE,
     *           N1FGE(KOFF1),NFG,XLAMH1(KOFF2),NDEL,
     *           ONE,WORK(KOFF3),1)
C
C     -------------------------------------------------------
C      So now we have our temporary result KFGDEL(fg,del).
C      We will now backtransform "g" index using XLAMP2(gam g):
C
C      KFGDEL(fg,del)    *    XLAMP2(gam g)
C          |
C          |loop over del
C          V
C      KFGDEL(f,g)       *    XLAMP2(gam g)  -->  KFGAM(f,gam)
C
C     -------------------------------------------------------
C
      DO ISYMG = 1,NSYM !"g" is summation index for this multiplication
                        !but we will stay within this loop till the
                        !end of the routine (otherwise I would have to
                        !loop over "gamma" in the very next contraction)
C
         !the symmetries necessary for this and the NEXT contraction
         ISYMF     = MULD2H(ISFG,ISYMG)
         ISYGAM    = MULD2H(ISLAMP2,ISYMG)!because: XLAMP2(gam g)
         ISYAL     = MULD2H(ISLAMP1,ISYMF)!because: XLAMP1(al  f)
         ISFGAM    = MULD2H(ISYMF,ISYGAM)
         ISALGAM   = MULD2H(ISYAL,ISYGAM)
C
         KFGAM  = KEND1 !we still need KFGDEL, that's why KEND1
         KEND2  = KFGAM + NEMAT1(ISFGAM)
C
         LWRK2 = LWORK - KEND2!for fixed "i" we loop over ISYMG and 
                              !we don't want memory to pile up; 
                              !that's why KEND2
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available: ', LWORK
            WRITE(LUPRI,*) 'Memory needed   : ', KEND2
            CALL QUIT('Insufficient memory in N1_2AORHF (2)')
         END IF
C
         CALL DZERO(WORK(KFGAM),NEMAT1(ISFGAM))
C
         !multiply KFGDEL(f,g) * XLAMP2(gam g) --> KFGAM(f,gam)
         KOFF1 = KFGDEL
     *         + IMATAB(ISYMF,ISYMG)
         KOFF2 = IGLMVI(ISYGAM,ISYMG)
     *         + 1
         KOFF3 = KFGAM
C
         NF   = MAX(NVIR(ISYMF),1)
         NGAM = MAX(NBAS(ISYGAM),1)
C
         CALL DGEMM('N','T',NVIR(ISYMF),NBAS(ISYGAM),
     *              NVIR(ISYMG),
     *              ONE,WORK(KOFF1),NF,XLAMP2(KOFF2),NGAM,
     *              ONE,WORK(KOFF3),NF)
C
C        -------------------------------------------------------
C         So now we have our temporary result KFGAM(f,gam).
C         What remains is to backtransform "f" index using 
C         XLAMP1(al  f) and put the result to the final storage
C         N1MAT2AORHF(al gam,i):
C
C         XLAMP1(al f) * KFGAM(f,gam)  -->  N1MAT2AORHF(al gam,i)
C
C         We are still inside "i" loop so addressing the
C         final storage N1MAT2AORHF(al gam,i) is straightforward.
C        -------------------------------------------------------
C
         !All symmetries have been determined after entering ISYMG loop
C
         !XLAMP1(al f) * KFGAM(f,gam) --> N1MAT2AORHF(al gam,i)
         KOFF1 = IGLMVI(ISYAL,ISYMF)
     *         + 1
         KOFF2 = KFGAM
         KOFF3 = IDSRHFSQ(ISALGAM,ISYMI)
     *         + N2BST(ISALGAM)*(I-1)
     *         + IAODIS(ISYAL,ISYGAM)
     *         + 1
C
         NAL = MAX(NBAS(ISYAL),1)
         NF  = MAX(NVIR(ISYMF),1)
C
         CALL DGEMM('N','N',NBAS(ISYAL),NBAS(ISYGAM),
     *              NVIR(ISYMF),
     *              ONE,XLAMP1(KOFF1),NAL,WORK(KOFF2),NF,
     *              ONE,N1MAT2AORHF(KOFF3),NAL)
C
      END DO !ISYMG
C
C----------
C     End.
C----------
C
      CALL QEXIT('N1AO_1')
C
      RETURN
      END
C
C  /* Deck n12ao_g3ao */
      SUBROUTINE N12AO_G3AO(N1MAT2AORHF,ISN12AORHF, ! N_{alpha,i,gamma}
     *                      AOINT,ISYMAO,           ! g_{al,be,gam,del}
     *                      XLAMH2,ISLAMH2,         ! LAM^{h}_{be a}
     *                      OMEGA1,ISOMEGA1,        !result vector
     *                      IDEL,ISYMD,             ! delta: fixed AO index
     *                      WORK,LWORK)
*
**********************************************************************
*                                                                    *
* Contract N^{del}_{al,i,gam} (sitting as N1MAT2AORHF(al gam,i) )    *
* with  g^{del}_{al,be,gam} integrals and then transform the "beta"  *
* index using LAM^{h}_{be a} to obtain finally OMEGA1(a,i).          *
*                                                                    *
* NOTE  "del" is a fixed AO index (coming from outside).             *
* ====                                                               *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* The routine has the following structure:                           *
* ========================================                           *
*                                                                    *
* 1) Read in g^{del}_{al,be,gam} and sort them to KALGABE(al gam,be) *
*                                                                    *
* 2) Multiply:                                                       *
*                                                                    *
*       KALGABE(al gam,be) * N1MAT2AORHF(al gam,i) --> KBEI(be,i)    *
*                                                                    *
* 3) Transform:                                                      *
*                                                                    *
*       LAM^{h}_{be a} * KBEI(be,i)  -->  OMEGA1(a,i)                *
*                                                                    *
*                                                                    *
*                                                                    *
**********************************************************************
*  F. Pawlowski, 17-Mar-2004, Aarhus.                                *
**********************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISN12AORHF,ISYMAO,ISLAMH2,ISOMEGA1,IDEL,ISYMD,LWORK
      INTEGER ISYABG,KALGABE,KEND1,LWRK1,ISYMG,ISALBE,KSCR1,KEND2,LWRK2
      INTEGER KOFF1,ISYMBE,ISYMAL,ISALGAM,BE,AL,KOFF2
      INTEGER ISBEI,KBEI,ISYALGA,ISYBE,ISYMI,KOFF3,NALGA,NBE
      INTEGER ISYMA,NA
C
#if defined (SYS_CRAY)
      REAL N1MAT2AORHF(*),AOINT(*),XLAMH2(*),OMEGA1(*)
      REAL WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION N1MAT2AORHF(*),AOINT(*),XLAMH2(*),OMEGA1(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('NGAO')
C
C     =======
C     STEP 1:
C     =======--------------------------------------------------------
C     Read in g^{del}_{al,be,gam} and sort them to KALGABE(al gam,be)
C     ---------------------------------------------------------------
C
      ISYABG = MULD2H(ISYMAO,ISYMD) !sym of g^{del}_{al,be,gam}
C
      KALGABE = 1
      KEND1   = KALGABE + NDISAOSQ(ISYABG) 
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in N12AO_G3AO (1)')
      END IF
C
      DO ISYMG = 1, NSYM !loop over gamma
         ISALBE = MULD2H(ISYABG,ISYMG)
C
         KSCR1 = KEND1
         KEND2 = KSCR1 + N2BST(ISALBE) !temp storage of (al,be) integrals
         LWRK2  = LWORK - KEND2        !KEND2 to avoid piling up memory
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available: ', LWORK
            WRITE(LUPRI,*) 'Memory needed   : ', KEND2
            CALL QUIT('Insufficient memory in N12AO_G3AO (2)')
         END IF
C
         DO G = 1, NBAS(ISYMG)
C
            !Read in g^{del}_{al,be,gam} distribution as KSCR1(al,be) 
            !for fixed del (loop from outside) and gam (internal loop)
            KOFF1 = IDSAOG(ISYMG,ISYMD) + NNBST(ISALBE)*(G-1) + 1
            CALL CCSD_SYMSQ(AOINT(KOFF1),ISALBE,WORK(KSCR1))
C
            !Put KSCR1(al,be) to KALGABE(al gam,be)
            DO ISYMBE = 1,NSYM
               ISYMAL  = MULD2H(ISALBE,ISYMBE)
               ISALGAM = MULD2H(ISYMAL,ISYMG)
C
               DO BE = 1,NBAS(ISYMBE)
                  DO AL = 1,NBAS(ISYMAL)
C
                     KOFF1 = KSCR1
     *                     + IAODIS(ISYMAL,ISYMBE)
     *                     + NBAS(ISYMAL)*(BE-1)
     *                     + AL-1
                     KOFF2 = KALGABE
     *                     + I3AO(ISALGAM,ISYMBE)
     *                     + N2BST(ISALGAM)*(BE-1)
     *                     + IAODIS(ISYMAL,ISYMG)
     *                     + NBAS(ISYMAL)*(G-1)
     *                     + AL-1
C
                     WORK(KOFF2) = WORK(KOFF1)
C
                  END DO !AL
               END DO !BE
            END DO !ISYMBE
         END DO !G
      END DO !ISYMG
C
C     =======
C     STEP 2:
C     =======------------------------------------------------------------
C     Multiply: KALGABE(al gam,be) * N1MAT2AORHF(al gam,i) --> KBEI(be,i)
C     -------------------------------------------------------------------
C
      ISBEI = MULD2H(ISYABG,ISN12AORHF)
C
      KBEI  = KEND1 !we still need KALGABE, so KEND1 is kept
      KEND1 = KBEI + NGLMRH(ISBEI)
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in N12AO_G3AO (3)')
      END IF
C
      CALL DZERO(WORK(KBEI),NGLMRH(ISBEI))
C
      DO ISYALGA = 1,NSYM
         ISYBE = MULD2H(ISYABG,ISYALGA)
         ISYMI = MULD2H(ISN12AORHF,ISYALGA)
C
         KOFF1 = KALGABE
     *         + I3AO(ISYALGA,ISYBE)
         KOFF2 = IDSRHFSQ(ISYALGA,ISYMI)
     *         + 1
         KOFF3 = KBEI
     *         + IGLMRH(ISYBE,ISYMI)
C
         NALGA = MAX(N2BST(ISYALGA),1)
         NBE   = MAX(NBAS(ISYBE),1)
C
         CALL DGEMM('T','N',NBAS(ISYBE),NRHF(ISYMI),N2BST(ISYALGA),
     *              ONE,WORK(KOFF1),NALGA,N1MAT2AORHF(KOFF2),NALGA,
     *              ONE,WORK(KOFF3),NBE)
C
      END DO !ISYALGA
C
C     =======
C     STEP 3:
C     =======------------------------------------------------------------
C     Transform:   LAM^{h}_{be a} * KBEI(be,i)  -->  OMEGA1(a,i)
C     -------------------------------------------------------------------
C
      DO ISYBE = 1,NSYM
         ISYMA = MULD2H(ISLAMH2,ISYBE)
         ISYMI = MULD2H(ISOMEGA1,ISYMA)
C
         KOFF1 = IGLMVI(ISYBE,ISYMA)
     *         + 1
         KOFF2 = KBEI
     *         + IGLMRH(ISYBE,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NBE = MAX(NBAS(ISYBE),1)
         NA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYBE),
     *              ONE,XLAMH2(KOFF1),NBE,WORK(KOFF2),NBE,
     *              ONE,OMEGA1(KOFF3),NA)
C
      END DO !ISYBE
C
C----------
C     End.
C----------
C
      CALL QEXIT('NGAO')
C
      RETURN
      END
C
C  /* Deck n1n2_g */
      SUBROUTINE N1N2_G(LUN1,FNN1,   ! file containing (fge,i)
     *                  ISYMN1,      ! symmetry of N1(fge,i)
     *                  N2MAT,ISYMN2,! --> N_anmo
     *                  XOVVO,XOOVV,XOOOO,ISINT, !integrals
     *                  OMEGA1,ISOMEGA1, !result array
     *                  INDSQ,LENSQ,! index associated with N2MAT
     *                  WORK,LWORK)
*
**********************************************************************
*                                                                    *
* Calculate <T3|[[H,T2],tau_ai]|HF> contribution to the left         *
* transformation based on N_fige and N_anmo intermediates            *
* (which have been constructed in WT2_N1N2 routine).                 *
*                                                                    *
* omega_ai = <T3|[[H,T2],tau_ai]|HF>                                 *
*          = - sum_efn g_iefn N_fnae  - sum_efn g_infe N_anfe        *
*            - sum_fmn g_mafn N_fnmi  - sum_fmn g_mnfa N_fimn        *
*            + sum_mno g_inmo N_anmo  + sum_efg g_fage N_fige        *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* NOTE !!!                                                           *
* ========                                                           *
* The very last term (sum_efg g_fage N_fige) involves V4 integrals   *
* and thus requires a special treatment. This term is not calculated *
* here, but in N1N2_GV4 routine.                                     *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* The terms involving N1MAT (N_fige):                                *
* ===================================                                *
*                                                                    *
* omega_ai <--  - sum_efn g_iefn N_fnae - sum_efn g_infe N_anfe      *
*                                                                    *
*      N1(fge,i) IS SITTING ON LUN1 FILE !!!!                        *
*      ==================================                            *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* The terms involving N2MAT (N_anmo):                                *
* ===================================                                *
*                                                                    *
*  omega_ai <--  - sum_fmn g_mafn N_fnmi  - sum_fmn g_mnfa N_fimn    *
*                + sum_mno g_inmo N_anmo                             *
*                                                                    *
*                                                                    *
*--------------------------------------------------------------------*
* omega_ai is stored as:                                             *
* ======================                                             *
*                                                                    *
*   OMEGA1(a,i)                                                      *
*                                                                    *
**********************************************************************
*  F. Pawlowski, 08-Mar-2004, Aarhus.                                *
**********************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CHARACTER*(*) FNN1
      INTEGER       LUN1
C     
      INTEGER ISYMN1,ISYMN2,ISINT,ISOMEGA1,LENSQ,INDSQ(LENSQ,6),LWORK
      INTEGER ISFNM,ISYMA,ISYMI,KOFF1,KOFF2,KOFF3,NFNM,NA
      INTEGER ISFMNI,KFMNI,KEND1,LWRK1
      INTEGER ISYMO,ISINM,ISYMM,ISIN,ISYMN,ISNM
      INTEGER ISFAEN,ISAFEN,ISFNIE,ISFENI,KAFEN,KFENI
      INTEGER IOPT
      INTEGER ISFEN,NFEN
      INTEGER ISFINE
      INTEGER ISFMN,NFMN,IO
C
      INTEGER ISAONM,ISIONM,KAONM,ISONM,NONM
C
      INTEGER ISFAE,ISFEA,ISFEI,KFAE,KFEA,KFEI
      INTEGER ISYME,ISFA,ISYMF,ISFE,ISFNI,ISFN,NFE
C
      INTEGER ISAFE,KAFE,ISAF,ISFIN,ISFI
      INTEGER IADR
#if defined (SYS_CRAY)
      REAL N2MAT(*),XOVVO(*),XOOVV(*),XOOOO(*),OMEGA1(*)
      REAL WORK(LWORK)
      REAL ONE
#else
      
      DOUBLE PRECISION N2MAT(*),XOVVO(*),XOOVV(*),XOOOO(*)
      DOUBLE PRECISION OMEGA1(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('NG')
C
      !Initial symmetry check
      IF (     (ISOMEGA1 .NE. MULD2H(ISYMN1,ISINT))
     *    .OR. (ISOMEGA1 .NE. MULD2H(ISYMN2,ISINT))  ) THEN
         WRITE(LUPRI,*) 'ISOMEGA1 = ', ISOMEGA1
         WRITE(LUPRI,*) 'SHOULD BE EQUAL TO'
         WRITE(LUPRI,*) 'ISYMN1 = ', ISYMN1
         WRITE(LUPRI,*) '( or ISYMN2 = ', ISYMN2, ')'
         WRITE(LUPRI,*) 'TIMES'
         WRITE(LUPRI,*) 'ISINT = ', ISINT
         CALL QUIT('Symmetry mismatch in N1N2_G')
      END IF
C
C     -------------------------------------------------------
C      Calculate the first contribution from N2MAT (N_anmo):
C
C       - sum_fmn g_mafn N_fnmi --> omega_ai
C
C      calculated as:
C
C        - XOVVO(fnm,a)   *   N2MAT(fnm,i)  --> OMEGA1(a,i)
C     -------------------------------------------------------
C
      DO ISFNM = 1,NSYM
         ISYMA = MULD2H(ISINT,ISFNM)
         ISYMI = MULD2H(ISYMN2,ISFNM)
C
         KOFF1 = IT2SP(ISFNM,ISYMA)
     *         + 1
         KOFF2 = ISAIKJ(ISFNM,ISYMI)
     *         + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NFNM = MAX(NCKI(ISFNM),1)
         NA   = MAX(NVIR(ISYMA),1)
C
Caddomega 5
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISFNM),
     *              -ONE,XOVVO(KOFF1),NFNM,N2MAT(KOFF2),NFNM,
     *              ONE,OMEGA1(KOFF3),NA) 
C
      END DO !ISFNM
C
C     -------------------------------------------------------
C      Calculate the second contribution from N2MAT (N_anmo):
C
C       - sum_fmn g_mnfa N_fimn --> omega_ai
C
C      calculated as:
C
C        - XOOVV(fmn,a)   *   N2MAT(fim,n)  
C                                 |
C                                 | indsq(4)
C                                 |
C                                 V
C        - XOOVV(fmn,a)   *   KFMNI(fmn,i)  -->  OMEGA1(a,i)
C     -------------------------------------------------------
C
      ISFMNI = ISYMN2
C
      KFMNI = 1
      KEND1 = KFMNI + NCKIJ(ISFMNI)
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in N1N2_G (1)')
      END IF
C
      !Sort N2MAT(fim,n) to KFMNI(fmn,i)
      DO I = 1,NCKIJ(ISFMNI)
         WORK(KFMNI+I-1) = N2MAT(INDSQ(I,4))
      END DO
C
      !Multiply:    - XOOVV(fmn,a)   *   KFMNI(fmn,i)  -->  OMEGA1(a,i)
      DO ISFMN = 1,NSYM
         ISYMA = MULD2H(ISINT,ISFMN)
         ISYMI = MULD2H(ISYMN2,ISFMN)
C
         KOFF1 = IT2SP(ISFMN,ISYMA)
     *         + 1
         KOFF2 = KFMNI
     *         + ISAIKJ(ISFMN,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NFMN = MAX(NCKI(ISFMN),1)
         NA   = MAX(NVIR(ISYMA),1)
C
C
Caddomega 6
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISFMN),
     *              -ONE,XOOVV(KOFF1),NFMN,WORK(KOFF2),NFMN,
     *              ONE,OMEGA1(KOFF3),NA)
C
      END DO !ISFMN
C
C     -------------------------------------------------------
C      Calculate the third contribution from N2MAT (N_anmo):
C
C         sum_mno N_anmo g_inmo --> omega_ai
C
C      calculated as:
C
C          N2MAT(anm,o)              *      XOOOO(onm,i)  
C           |                                    
C           | <-- indsq(2)
C           |                                    
C           V                                    
C         KAONM(aon,m)
C           |                                    
C           | <-- CALL CC3_SRTVOOO                 
C           |                                    
C           V                                    
C          KAONM(a,onm)              *      XOOOO(onm,i)  -->  OMEGA1(a,i)
C     -------------------------------------------------------
C
      ISAONM = ISYMN2
      ISIONM = ISINT
C
      KAONM = 1
      KEND1 = KAONM + N3VOOO(ISAONM)
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in N1N2_G (2)')
      END IF
C
      DO I = 1,NCKIJ(ISAONM)
         WORK(KAONM+I-1) = N2MAT(INDSQ(I,2))
      END DO


      IF (NSYM .GT. 1) THEN
         IF (LWRK1 .LT. N3VOOO(ISAONM)) THEN
            WRITE(LUPRI,*)'Memory available: ', LWRK1
            WRITE(LUPRI,*)'Memory needed   : ', N3VOOO(ISAONM)
            CALL QUIT('Insufficient space in n1n2_g (2a)')
         END IF
         !Sort KAONM(aon,m) to KAONM(a,onm)
         CALL CC3_SRTVOOO(WORK(KEND1),WORK(KAONM),ISAONM)
         CALL DCOPY(N3VOOO(ISAONM),WORK(KEND1),1,WORK(KAONM),1)
      END IF
C
      !Multiply  KANMO(a,nmo)  *   XOOOO(nmo,i)  -->  OMEGA1(a,i)
      DO ISONM = 1,NSYM
         ISYMA = MULD2H(ISAONM,ISONM)
         ISYMI = MULD2H(ISIONM,ISONM)
C
         KOFF1 = KAONM
     *         + I3VOOO(ISYMA,ISONM)
         KOFF2 = I3ORHF(ISONM,ISYMI)
     *         + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NA   = MAX(NVIR(ISYMA),1)
         NONM   = MAX(NMAIJK(ISONM),1)
C
C
Caddomega 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NMAIJK(ISONM),
     *              ONE,WORK(KOFF1),NA,XOOOO(KOFF2),NONM,
     *              ONE,OMEGA1(KOFF3),NA)
C
      END DO !ISNMO
C
C     -------------------------------------------------------
C      Calculate the first contribution from N1MAT (N_fige):
C
C       - sum_efn g_iefn N_fnae --> omega_ai
C
C      calculated as (for fixed N):
C
C        -  N1MAT(fae,N)    *       XOVVO(fNi,e)  
C               |                       |         
C               |                       |         
C               |                       |         
C               V                       |                      
C        -  KFAE(fa,e)                  |
C               |                       |
C               |                       |
C               |                       |
C               V                       V
C        - KFEA(fe,a)        *       KFEI(fe,i)   -->  OMEGA1(a,i)
C
C     -------------------------------------------------------
C
      ISFAEN = ISYMN1
      ISFNIE = ISINT
C
      DO ISYMN = 1,NSYM
         ISFAE = MULD2H(ISFAEN,ISYMN)
         ISFEA = ISFAE
         ISFEI = MULD2H(ISFNIE,ISYMN)
C
         !Allocations below are used to calculate this and the next 
         !contribution
C
         KFAE  = 1
         KEND1 = KFAE + NMAABC(ISFAE)
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available: ', LWORK
            WRITE(LUPRI,*) 'Memory needed   : ', KEND1
            CALL QUIT('Insufficient memory in N1N2_G (3)')
         END IF
C
         KFEA  = KEND1
         KEND1 = KFEA + NMAABC(ISFEA)
C
         KFEI  = KEND1
         KEND1 = KFEI + NCKATR(ISFEI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available: ', LWORK
            WRITE(LUPRI,*) 'Memory needed   : ', KEND1
            CALL QUIT('Insufficient memory in N1N2_G (4)')
         END IF
C
         DO N = 1,NRHF(ISYMN)
C
            CALL DZERO(WORK(KFEA),NMAABC(ISFEA))
            CALL DZERO(WORK(KFAE),NMAABC(ISFAE))
C
            !Get KFAE(fae)^N from file
            IF (NMAABC(ISFAE).GT.0) THEN
              IADR = IMAABCI(ISFAE,ISYMN) + NMAABC(ISFAE)*(N-1) + 1
              CALL GETWA2(LUN1,FNN1,WORK(KFAE),IADR,NMAABC(ISFAE))
            END IF
C
            !Sort KFAE(fa,e) to KFEA(fe,a)
            DO ISYME = 1,NSYM
               ISFA = MULD2H(ISFAE,ISYME)
               DO ISYMA = 1,NSYM
                  ISYMF = MULD2H(ISFA,ISYMA)
                  ISFE  = MULD2H(ISYMF,ISYME)
                  DO E = 1,NVIR(ISYME)
                     DO A = 1,NVIR(ISYMA)
                        DO F = 1,NVIR(ISYMF)
C
                           KOFF1 = KFAE 
     *                           + IMAABC(ISFA,ISYME)
     *                           + NMATAB(ISFA)*(E-1)
     *                           + IMATAB(ISYMF,ISYMA)
     *                           + NVIR(ISYMF)*(A-1)
     *                           + F-1
                           KOFF2 = KFEA 
     *                           + IMAABC(ISFE,ISYMA)
     *                           + NMATAB(ISFE)*(A-1)
     *                           + IMATAB(ISYMF,ISYME)
     *                           + NVIR(ISYMF)*(E-1)
     *                           + F-1
C
                           WORK(KOFF2) = WORK(KOFF1)
C
                        END DO!F
                     END DO!A
                  END DO!E
               END DO!ISYMA
            END DO!ISYME
C
            !Sort XOVVO(fNi,e) to KFEI(fe,i)
            DO ISYMI = 1,NSYM
               ISFE = MULD2H(ISFEI,ISYMI)
               DO ISYME = 1,NSYM
                  ISYMF = MULD2H(ISFE,ISYME)
                  ISFNI = MULD2H(ISFNIE,ISYME)
                  ISFN  = MULD2H(ISFNI,ISYMI)
                  DO I = 1,NRHF(ISYMI)
                     DO E = 1,NVIR(ISYME)
                        DO F = 1,NVIR(ISYMF)
C
                           KOFF1 = IT2SP(ISFNI,ISYME)
     *                           + NCKI(ISFNI)*(E-1)
     *                           + ICKI(ISFN,ISYMI)
     *                           + NT1AM(ISFN)*(I-1)
     *                           + IT1AM(ISYMF,ISYMN)
     *                           + NVIR(ISYMF)*(N-1)
     *                           + F
C
                           KOFF2 = KFEI 
     *                           + IMAABI(ISFE,ISYMI)
     *                           + NMATAB(ISFE)*(I-1)
     *                           + IMATAB(ISYMF,ISYME)
     *                           + NVIR(ISYMF)*(E-1)
     *                           + F-1
C
                           WORK(KOFF2) = XOVVO(KOFF1)
C
                        END DO!F
                     END DO!E
                  END DO!I
               END DO!ISYME
            END DO!ISYMI
C
            !Multiply:   - KFEA(fe,a)   *   KFEI(fe,i)   -->   OMEGA1(a,i)
            DO ISFE = 1,NSYM
               ISYMA = MULD2H(ISFEA,ISFE)
               ISYMI = MULD2H(ISFEI,ISFE)
C
               KOFF1 = KFEA
     *               + IMAABC(ISFE,ISYMA)
               KOFF2 = KFEI
     *               + IMAABI(ISFE,ISYMI)
               KOFF3 = IT1AM(ISYMA,ISYMI)
     *               + 1
C
               NFE   = MAX(NMATAB(ISFE),1)
               NA    = MAX(NVIR(ISYMA),1)
C
Caddomega3
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATAB(ISFE),
     *                    -ONE,WORK(KOFF1),NFE,WORK(KOFF2),NFE,
     *                    ONE,OMEGA1(KOFF3),NA)
C
            END DO!ISFE

C
C           -------------------------------------------------------
C            Calculate the second contribution from N1MAT (N_fige):
C
C             - sum_efn g_infe N_anfe --> omega_ai
C
C            calculated as:
C
C              -  KAFE(af,e)    *       XOOVV(fiN,e)  
C                     |                       |         
C                     |                       |
C                     |                       |                    
C                     |                       |
C                     V                       V
C              - KFEA(fe,a)     *       KFEI(fe,i)  -->  OMEGA1(a,i)
C
C           -------------------------------------------------------
C
            !Just change the name of KFAE array and ISFAE symmetry...
            ISAFE = ISFAE
            KAFE  = KFAE
C
            !sort KAFE(af,e) to KFEA(fe,a)
            DO ISYME = 1,NSYM
               ISAF = MULD2H(ISAFE,ISYME)
               DO ISYMF = 1,NSYM
                  ISYMA = MULD2H(ISAF,ISYMF)
                  ISFE  = MULD2H(ISYMF,ISYME)
                  DO E = 1,NVIR(ISYME)
                     DO F = 1,NVIR(ISYMF)
                        DO A = 1,NVIR(ISYMA)
C
                           KOFF1 = KAFE
     *                           + IMAABC(ISAF,ISYME)
     *                           + NMATAB(ISAF)*(E-1)
     *                           + IMATAB(ISYMA,ISYMF)
     *                           + NVIR(ISYMA)*(F-1)
     *                           + A-1
                           KOFF2 = KFEA
     *                           + IMAABC(ISFE,ISYMA)
     *                           + NMATAB(ISFE)*(A-1)
     *                           + IMATAB(ISYMF,ISYME)
     *                           + NVIR(ISYMF)*(E-1)
     *                           + F-1
C
                           WORK(KOFF2) = WORK(KOFF1)
C
                        END DO!F
                     END DO!A
                  END DO!E
               END DO!ISYMA
            END DO!ISYME
C
            !Sort XOOVV(fiN,e) to KFEI(fe,i)
            DO ISYMI = 1,NSYM
               ISFE = MULD2H(ISFEI,ISYMI)
               DO ISYME = 1,NSYM
                  ISYMF = MULD2H(ISFE,ISYME)
                  ISFIN = MULD2H(ISFNIE,ISYME)
                  ISFI  = MULD2H(ISFIN,ISYMN)
                  DO I = 1,NRHF(ISYMI)
                     DO E = 1,NVIR(ISYME)
                        DO F = 1,NVIR(ISYMF)
C
                           KOFF1 = IT2SP(ISFIN,ISYME)
     *                           + NCKI(ISFIN)*(E-1)
     *                           + ICKI(ISFI,ISYMN)
     *                           + NT1AM(ISFI)*(N-1)
     *                           + IT1AM(ISYMF,ISYMI)
     *                           + NVIR(ISYMF)*(I-1)
     *                           + F
C
                           KOFF2 = KFEI
     *                           + IMAABI(ISFE,ISYMI)
     *                           + NMATAB(ISFE)*(I-1)
     *                           + IMATAB(ISYMF,ISYME)
     *                           + NVIR(ISYMF)*(E-1)
     *                           + F-1
C
                           WORK(KOFF2) = XOOVV(KOFF1)
C
                        END DO!F
                     END DO!E
                  END DO!I
               END DO!ISYME
            END DO!ISYMI
C
            !Multiply:   - KFEA(fe,a)   *   KFEI(fe,i)   -->   OMEGA1(a,i)
            DO ISFE = 1,NSYM
               ISYMA = MULD2H(ISFEA,ISFE)
               ISYMI = MULD2H(ISFEI,ISFE)
C
               KOFF1 = KFEA
     *               + IMAABC(ISFE,ISYMA)
               KOFF2 = KFEI
     *               + IMAABI(ISFE,ISYMI)
               KOFF3 = IT1AM(ISYMA,ISYMI)
     *               + 1
C
               NFE   = MAX(NMATAB(ISFE),1)
               NA    = MAX(NVIR(ISYMA),1)
C
Caddomega4
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATAB(ISFE),
     *                    -ONE,WORK(KOFF1),NFE,WORK(KOFF2),NFE,
     *                    ONE,OMEGA1(KOFF3),NA)
C
            END DO!ISFE
C
         END DO !N
      END DO !ISYMN
C
C----------
C     End.
C----------
C
      CALL QEXIT('NG')
C
      RETURN
      END
C
C  /* Deck sort_t2_abij */
      SUBROUTINE SORT_T2_ABIJ(T2ABIJ,T2TP,ISYMT2)
C
C-------------------------------
C     Sort T2TP(aij,b) as T2(abi,j)
C-------------------------------
C
C     F. Pawlowski, Aarhus, Winter 2004
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMT2,ISYMB,ISAIJ,ISYMJ,ISAI,ISABI,ISYMI,ISYMA,ISAB
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL T2ABIJ(*),T2TP(*)
#else
      DOUBLE PRECISION T2ABIJ(*),T2TP(*)
#endif
C
      CALL QENTER('SRTABIJ')
C
C     Sort T2TP(aij,b) as T2(abi,j)
C
      DO ISYMB = 1,NSYM
         ISAIJ = MULD2H(ISYMT2,ISYMB)
         DO ISYMJ = 1,NSYM
            ISAI  = MULD2H(ISAIJ,ISYMJ)
            ISABI = MULD2H(ISAI,ISYMB)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISAI,ISYMI)
               ISAB  = MULD2H(ISYMA,ISYMB)
               DO B = 1,NVIR(ISYMB)
                  DO J  =1,NRHF(ISYMJ)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)
C
                           KOFF1 = IT2SP(ISAIJ,ISYMB)
     *                           + NCKI(ISAIJ)*(B-1)
     *                           + ICKI(ISAI,ISYMJ)
     *                           + NT1AM(ISAI)*(J-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
                           KOFF2 = IMAJBAIT(ISABI,ISYMJ)
     *                           + NMAIAB(ISABI)*(J-1)
     *                           + IMAABI(ISAB,ISYMI)
     *                           + NMATAB(ISAB)*(I-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A 
C
                           T2ABIJ(KOFF2) = T2TP(KOFF1)
C
                        END DO !A
                     END DO !I
                  END DO !J
               END DO !B
            END DO !ISYMI
         END DO !ISYMJ
      END DO !ISYMB
C
      CALL QEXIT('SRTABIJ')
C
      RETURN
      END
C  /* Deck sort_t2_abji */
      SUBROUTINE SORT_T2_ABJI(T2ABJI,T2TP,ISYMT2)
C
C-------------------------------
C     Sort T2TP(aij,b) as T2(abj,i)
C-------------------------------
C
C     F. Pawlowski, Aarhus, Winter 2004
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMT2,ISYMB,ISAIJ,ISYMJ,ISAI,ISYMI,ISYMA,ISAB,ISABJ
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL T2ABJI(*),T2TP(*)
#else
      DOUBLE PRECISION T2ABJI(*),T2TP(*)
#endif

C
      CALL QENTER('SRTABJI')
C
C     Sort T2TP(aij,b) as T2(abj,i)
C
      DO ISYMB = 1,NSYM
         ISAIJ = MULD2H(ISYMT2,ISYMB)
         DO ISYMJ = 1,NSYM
            ISAI  = MULD2H(ISAIJ,ISYMJ)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISAI,ISYMI)
               ISAB  = MULD2H(ISYMA,ISYMB)
               ISABJ = MULD2H(ISAB,ISYMJ)
               DO B = 1,NVIR(ISYMB)
                  DO J  =1,NRHF(ISYMJ)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)
C
                           KOFF1 = IT2SP(ISAIJ,ISYMB)
     *                           + NCKI(ISAIJ)*(B-1)
     *                           + ICKI(ISAI,ISYMJ)
     *                           + NT1AM(ISAI)*(J-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
                           KOFF2 = IMAJBAIT(ISABJ,ISYMI)
     *                           + NMAIAB(ISABJ)*(I-1)
     *                           + IMAABI(ISAB,ISYMJ)
     *                           + NMATAB(ISAB)*(J-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
C
                           T2ABJI(KOFF2) = T2TP(KOFF1)
               
C
                        END DO !A
                     END DO !I
                  END DO !J
               END DO !B
            END DO !ISYMI
         END DO !ISYMJ
      END DO !ISYMB
C
      CALL QEXIT('SRTABJI')
C
      RETURN
      END
C  /* Deck wt2_n1n2 */
      SUBROUTINE WT2_N1N2(WMAT,ISWMAT,!ISWMAT: total sysmetry of T3 (6 indeces)
     *                    T2TP,ISYMT2,
     *                    N1GEI,N1FEI, ! OUTPUT arrays: (ge,i)^F and (fe,i)^G
     *                    ISYMN1,      ! symmetry of  N_fige
     *                    N2MAT,ISYMN2,! --> N_anmo
     *                    IB,ISYMIB,ID,ISYMID,
     *                    INDSQ,LENSQ,! index associated with WMAT
     *                    INDSQN,LENSQN,! index associated with N2MAT
     *                    WORK,LWORK,
     *                    W3X)! .false.=>WMAT contains "complete" TMAT,
                              ! .true.=>WMAT contains W intermediate.
*
**********************************************************************
*                                                                    *
* Calculate N_fige and N_anmo intermediates needed for the           *
* calculation of <T3|[[H,T2],tau_ai]|HF> vector.                     *
*                                                                    *
* N_fige = sum_dlm tbar^{dfg}_{lim} t2^{de}_{lm}.                    *
* N_anmo = sum_dle tbar^{dae}_{lno} t2^{de}_{lm}.                    *
*                                                                    *
* tbar maybe either zero- or first-order, BUT:                       *
*   IT IS ALWAYS ASSUMED TO BE SITTING AS W INTERMEDIATE !!!         *
*                                                                    *
*   THEREFORE we introduced W3X flag:                                *
*                                                                    *
*      - W3X = .TRUE. : W intermediate is sitting in WMAT array      *
*                       and ALL the contributions to N_fige          *
*                       and N_anmo in the expressions below          *
*                       are calculated;                              *
*                                                                    *
*      - W3X = .FALSE. : "complete" T (obtained from the call to     *
*                        GET_T3BAR0_BD) is sitting in WMAT array     *
*                        and ONLY ONE contribution to N_fige         *
*                        and ONLY ONE contribution to N_anmo         *
*                        is calculated; these contributions are      *
*                        marked by (*) in the expressions below.     *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* N_fige is calculated as:                                           *
* ======                                                             *
*                                                                    *
* N_fige = sum_dlm (W^df(gmil) + W^gf(dlim) + W^dg(fiml)(*) ) t^{de}_{lm}*
*                                                                    *
* It is stored in two separate intermediates:                        *
*                                                                    *
*    N1GEI(ge,i)^F = sum_dlm (W^df(gmil) + W^gf(dlim)) t^{de}_{lm}   *
* and                                                                *
*    N1FEI(fe,i)^G = sum_dlm W^dg(fiml)(*)  ) t^{de}_{lm}            *
*                                                                    *
* This is done to avoid storage of whole N_fige array (which has     *
* size VVVO) in memory.                                              *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* N_anmo is calculated as:                                           *
* ======                                                             *
*                                                                    *
* N_anmo = sum_dle (W^ae(dlon)(*)  + W^da(eonl) + W^de(anol)) t^{de}_{lm}*
*                                                                    *
* It is stored as:                                                   *
*                                                                    *
*    N2MAT(anm,o)                                                    *
*                                                                    *
**********************************************************************
*  F. Pawlowski, 25-Feb-2004, Aarhus.
**********************************************************************
*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      LOGICAL W3X
C
      INTEGER ISWMAT,ISYMT2,ISYMN1,ISYMN2,IB,ISYMIB,ID,ISYMID,LENSQ
      INTEGER INDSQ(LENSQ,6),LWORK
      INTEGER LENSQN,INDSQN(LENSQN,6)
      INTEGER ISYMF,ISYMG,ISYMFG,ISDLMI,ISYMEI
      INTEGER KTMP,KEND1,LWRK1
      INTEGER ISDLM,ISYME,ISYMI,KOFF1,KOFF2,KOFF3,NDLM,NE,ISFGE
      INTEGER ISYMD,ISYMDF,ISGMIL,ISGILM,ISLME,ISGIE,KLME
      INTEGER KGILM,KEND2,LWRK2
      INTEGER ISLM,ISYMM,ISYMDL,ISYML
      INTEGER ISGI,NGI,NLM,ISFGEI,ISFG
      INTEGER ISYMDG,ISFIML,ISFILM,ISFIE,KFILM
      INTEGER ISFI,NFI
C
      INTEGER ISYMA,ISYMAE,ISDLON,ISONM,KDLON
      INTEGER ISDL,ISON,NDL,NON
      INTEGER ISANMO,ISYMO,ISANM,ISAN,ISYMN
      INTEGER ISYMDA,ISEONL,ISELON,ISELM,KELON,KELM,ISEL,NEL
      INTEGER ISYMDE,ISANOL,ISANOM,KLM,ISANO,NANO,NL
      INTEGER IO
C
      INTEGER KEI,KGIE,KDLMI
      INTEGER ISGEI,ISGE,IADR
C
      INTEGER ISFEI,KFIE,ISFE
c
      integer isymano
C
#if defined (SYS_CRAY)
      REAL WMAT(*),T2TP(*),N2MAT(*)
      REAL WORK(LWORK)
      REAL ONE
      REAL N1GEI(*),N1FEI(*)
#else
      
      DOUBLE PRECISION WMAT(*),T2TP(*),N2MAT(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE
      DOUBLE PRECISION N1GEI(*),N1FEI(*)
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('WTN')
C
      !Initial symmetry check
      IF (     (ISYMN1 .NE. MULD2H(ISWMAT,ISYMT2))
     *    .OR. (ISYMN2 .NE. MULD2H(ISWMAT,ISYMT2)) ) THEN
         WRITE(LUPRI,*) 'ISYMN1 = ', ISYMN1
         WRITE(LUPRI,*) 'AND'
         WRITE(LUPRI,*) 'ISYMN2 = ', ISYMN2
         WRITE(LUPRI,*) 'SHOULD BE EQUAL TO'
         WRITE(LUPRI,*) 'ISWMAT = ', ISWMAT
         WRITE(LUPRI,*) 'TIMES'
         WRITE(LUPRI,*) 'ISYMT2 = ', ISYMT2
         CALL QUIT('Symmetry mismatch in WT2_N1N2')
      END IF
C
C     -------------------------------------------------------
C     Sort T2TP(Dlm,e) to KLME(lm,e) and keep it in memory.
C     It will be needed for both first and third contribution
C     to N_fige.
C     -------------------------------------------------------
C
      !Nomenclature for D is the same as in the first and third
      !contribution to N_fige:
      ISYMD  = ISYMIB 
      ISLME  = MULD2H(ISYMT2,ISYMD)
C
      D = IB
C
      KLME  = 1
      KEND1 = KLME  + NCKI(ISLME) !we want to keep that array...
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in WT2_N1N2 (0)')
      END IF
C
      !Sort T2TP(Dlm,e) to KLME(lm,e)
      DO ISYME = 1,NSYM
         ISDLM = MULD2H(ISYMT2,ISYME)
         ISLM  = MULD2H(ISDLM,ISYMD)
         DO ISYMM = 1,NSYM
            ISYMDL = MULD2H(ISDLM,ISYMM)
            ISYML  = MULD2H(ISLM,ISYMM)
            DO E = 1,NVIR(ISYME)
               DO M = 1,NRHF(ISYMM)
                  DO L = 1,NRHF(ISYML)
C
                  KOFF1 = IT2SP(ISDLM,ISYME)
     *                  + NCKI(ISDLM)*(E-1)
     *                  + ICKI(ISYMDL,ISYMM)
     *                  + NT1AM(ISYMDL)*(M-1)
     *                  + IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L-1)
     *                  + D
                  KOFF2 = KLME
     *                  + IMAIJA(ISLM,ISYME)
     *                  + NMATIJ(ISLM)*(E-1)
     *                  + IMATIJ(ISYML,ISYMM)
     *                  + NRHF(ISYML)*(M-1)
     *                  + L-1
C
                  WORK(KOFF2) = T2TP(KOFF1)
C
                  END DO !L
               END DO !M
            END DO !E
         END DO !ISYMM
      END DO !ISYME
C

C
C     ----------------------------------------------------
      IF (W3X) THEN !Calculate ALL contributions to N_fige
C     ----------------------------------------------------
C

C
C        ----------------------------------------------
C         Calculate the first contribution to N_fige:
C         
C         sum_dlm W^df(gmil) * t^{de}_{lm}  --> N_fige
C
C         calculated as:
C
C         WMAT^DF(gmi,l)       T2TP(Dlm,e)
C              |                    |
C              | indsq(4)           |
C              V                    |
C          KGILM(gil,m)             |(has been sorted above)
C              |                    |
C              | indsq(6)           |
C              V                    V
C          KGILM(gi,lm)  *   KLME(lm,e) --> KGIE(gi,e)
C
C          KGIE(gi,e) + KEI(e,i) --> N1GEI(ge,i)^F
C        ----------------------------------------------
C
         D = IB
         F = ID
C
         ISYMD  = ISYMIB
         ISYMF  = ISYMID
         ISYMDF = MULD2H(ISYMD,ISYMF)
         ISGMIL = MULD2H(ISWMAT,ISYMDF)
         ISGILM = ISGMIL
         ISLME  = MULD2H(ISYMT2,ISYMD)
         ISGIE  = MULD2H(ISGILM,ISLME)
         ISGEI  = ISGIE
C
         KGILM = KEND1 !we want to keep KLME; that's why we start with KEND1
         KGIE  = KGILM + NCKIJ(ISGILM)
         KEND2 = KGIE  + NCKATR(ISGIE)
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available: ', LWORK
            WRITE(LUPRI,*) 'Memory needed   : ', KEND2
            CALL QUIT('Insufficient memory in WT2_N1N2 (1)')
         END IF
C
         CALL DZERO(WORK(KGIE),NCKATR(ISGIE))
C
         !Sort WMAT^DF(gmi,l) to KGILM(gil,m)
         DO I = 1,NCKIJ(ISGMIL)
            WORK(KGILM+I-1) = WMAT(INDSQ(I,4))
         END DO
C
         !If symmetry, sort KGILM(gil,m) to KGILM(gi,lm)
         IF (NSYM .GT. 1) THEN
           IF (LWRK2 .LT. NCKIJ(ISGILM)) THEN
              WRITE(LUPRI,*) 'Memory available: ', LWRK2
              WRITE(LUPRI,*) 'Memory needed   : ', NCKIJ(ISGILM)
              CALL QUIT('Insufficient memory in WT2_N1N2 (2)')
           END IF
           CALL CC_GATHER(NCKIJ(ISGILM),WORK(KEND2),WORK(KGILM),
     *                    INDSQ(1,6))
           CALL DCOPY(NCKIJ(ISGILM),WORK(KEND2),1,WORK(KGILM),1)
         END IF
C
         !Multiply KGILM(gi,lm)  *   KLME(lm,e) --> KGIE(gi,e)
         DO ISLM = 1,NSYM
            ISGI  = MULD2H(ISGILM,ISLM)
            ISYME = MULD2H(ISLME,ISLM)
C
            KOFF1 = KGILM
     *            + ISAIKL(ISGI,ISLM)
            KOFF2 = KLME
     *            + IMAIJA(ISLM,ISYME)
            KOFF3 = KGIE
     *            + ICKATR(ISGI,ISYME)
C
            NGI   = MAX(NT1AM(ISGI),1)
            NLM   = MAX(NMATIJ(ISLM),1)
C
            CALL DGEMM('N','N',NT1AM(ISGI),NVIR(ISYME),NMATIJ(ISLM),
     *                 ONE,WORK(KOFF1),NGI,WORK(KOFF2),NLM,
     *                 ONE,WORK(KOFF3),NGI)
C
         END DO ! ISLM
C
C        ----------------------------------------------
C         Calculate the second contribution to N_fige:
C         
C         sum_dlm W^gf(dlim) * t^{de}_{lm}  --> N_fige
C                     |
C                     | indsq(3)
C                     V
C                   KDLMI(dlm,i)
C
C         T2TP(dlm,e) * KDLMI(dlm,i) --> KEI(e,i)
C
C        ----------------------------------------------
C
         F      = ID
         G      = IB
C
         ISYMF  = ISYMID
         ISYMG  = ISYMIB
         ISYMFG = MULD2H(ISYMF,ISYMG)
         ISDLMI = MULD2H(ISWMAT,ISYMFG)
         ISYMEI = MULD2H(ISYMT2,ISDLMI)
C
         KEI  = KEND2
         KEND2 = KEI  + NT1AM(ISYMEI) 
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available: ', LWORK
            WRITE(LUPRI,*)'Memory needed:    ', KEND2
            CALL QUIT('Insufficient memory in WT2_N1N2 (3)')
         END IF
C
         KDLMI = KEND2
         KEND2 = KDLMI + NCKIJ(ISDLMI)!temporary storage
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available: ', LWORK
            WRITE(LUPRI,*)'Memory needed:    ', KEND2
            CALL QUIT('Insufficient memory in WT2_N1N2 (3x)')
         END IF
C
         CALL DZERO(WORK(KEI),NT1AM(ISYMEI))
C
         !Sort WMAT^GF(dli,m) to KDLMI(dlm,i)
         DO I = 1,NCKIJ(ISDLMI)
            WORK(KDLMI+I-1) = WMAT(INDSQ(I,3))
         END DO
C
         !Multiply T2TP(dlm,e) * KDLMI(dlm,i) --> KEI(e,i)
         DO ISDLM = 1,NSYM
            ISYME = MULD2H(ISYMT2,ISDLM)
            ISYMI = MULD2H(ISDLMI,ISDLM)
C
            KOFF1 = IT2SP(ISDLM,ISYME)
     *            + 1
            KOFF2 = ISAIKJ(ISDLM,ISYMI)
     *            + KDLMI
            KOFF3 = KEI + IT1AM(ISYME,ISYMI)
C
            NDLM = MAX(NMAIJA(ISDLM),1)
            NE   = MAX(NVIR(ISYME),1)
C
            CALL DGEMM('T','N',NVIR(ISYME),NRHF(ISYMI),NMAIJA(ISDLM),
     *                 ONE,T2TP(KOFF1),NDLM,WORK(KOFF2),NDLM,
     *                 ONE,WORK(KOFF3),NE)
C
         END DO !ISDLM
C
         !Put  KGIE(gi,e) and KEI(e,i) together to N1GEI(ge,i)^F
         ISFGEI = MULD2H(ISGIE,ISYMF)
         DO ISYMI = 1,NSYM
            ISFGE = MULD2H(ISFGEI,ISYMI)
            DO ISYME = 1,NSYM
               ISFG = MULD2H(ISFGE,ISYME)
               ISYMG = MULD2H(ISFG,ISYMF)
               ISGI  = MULD2H(ISYMG,ISYMI)
               ISGE  = MULD2H(ISYMG,ISYME)
               DO I = 1,NRHF(ISYMI)
                  DO E = 1,NVIR(ISYME)
                     DO G = 1,NVIR(ISYMG)
C
                        KOFF1 = KGIE
     *                        + ICKATR(ISGI,ISYME)
     *                        + NT1AM(ISGI)*(E-1)
     *                        + IT1AM(ISYMG,ISYMI)
     *                        + NVIR(ISYMG)*(I-1)
     *                        + G-1
                        KOFF2 = KEI
     *                        + IT1AM(ISYME,ISYMI)
     *                        + NVIR(ISYME)*(I-1)
     *                        + E-1
                        KOFF3 = IMAABI(ISGE,ISYMI)
     *                        + NMATAB(ISGE)*(I-1)
     *                        + IMATAB(ISYMG,ISYME)
     *                        + NVIR(ISYMG)*(E-1)
     *                        + G
C
                        IF ((ISYMG.EQ.ISYMIB) .AND. (G.EQ.IB)) THEN
C
                          N1GEI(KOFF3) = N1GEI(KOFF3) + WORK(KOFF1)
     *                                                + WORK(KOFF2)
                        ELSE
C
                          N1GEI(KOFF3) = N1GEI(KOFF3)  + WORK(KOFF1)
                        END IF
C
                     END DO !G
                  END DO !E
               END DO !I
            END DO !ISYME
         END DO !ISYMI
C
C     ----------------------------
      END IF !W3X
C     ----------------------------
C

C
C     ----------------------------------------------
C      Calculate the third contribution to N_fige:
C      
C      sum_dlm W^dg(fiml) * t^{de}_{lm}  --> N_fige
C
C      calculated as:
C
C      WMAT^DG(fim,l)       T2TP(Dlm,e)
C           |                    |
C           | indsq(3)           |
C           V                    |
C       KFILM(fil,m)             |  <-- use already existing KLME array
C           |                    |
C           | indsq(6)           |
C           V                    V
C       KFILM(fi,lm)  *   KLME(lm,e) --> KFIE(fi,e)
C
C       KFIE(fi,e) --> N1FEI(fe,i)^G
C     ----------------------------------------------
C

      D = IB
      G = ID
C
      ISYMD  = ISYMIB
      ISYMG  = ISYMID
      ISYMDG = MULD2H(ISYMD,ISYMG)
      ISFIML = MULD2H(ISWMAT,ISYMDG)
      ISFILM = ISFIML
      ISLME  = MULD2H(ISYMT2,ISYMD)
      ISFIE  = MULD2H(ISFILM,ISLME)
      ISFEI  = ISFIE
C
      KFILM = KEND1 !want to keep KLME; the rest is useless
      KFIE  = KFILM + NCKIJ(ISFIML)
      KEND2 = KFIE  + NCKATR(ISFIE)
      LWRK2 = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND2
         CALL QUIT('Insufficient memory in WT2_N1N2 (4)')
      END IF
C
      CALL DZERO(WORK(KFIE),NCKATR(ISFIE))
C
      !Sort WMAT^DG(fim,l) to KFILM(fil,m)
      DO I = 1,NCKIJ(ISFIML)
         WORK(KFILM+I-1) = WMAT(INDSQ(I,3))
      END DO
C
      !If symmetry sort KFILM(fil,m)  to KFILM(fi,lm)
      IF (NSYM .GT. 1) THEN
        IF (LWRK2 .LT. NCKIJ(ISFILM)) THEN
           WRITE(LUPRI,*) 'Memory available: ', LWRK2
           WRITE(LUPRI,*) 'Memory needed   : ', NCKIJ(ISGILM)
           CALL QUIT('Insufficient memory in WT2_N1N2 (5)')
        END IF 
        CALL CC_GATHER(NCKIJ(ISFILM),WORK(KEND2),WORK(KFILM),INDSQ(1,6))
        CALL DCOPY(NCKIJ(ISFILM),WORK(KEND2),1,WORK(KFILM),1)
      END IF
C
      !Multiply KFILM(fi,lm)   *   KLME(lm,e) --> KFIE(fi,e)
      DO ISLM = 1,NSYM
         ISFI  = MULD2H(ISFILM,ISLM)
         ISYME = MULD2H(ISLME,ISLM)
C
         KOFF1 = KFILM
     *         + ISAIKL(ISFI,ISLM)
         KOFF2 = KLME
     *         + IMAIJA(ISLM,ISYME)
         KOFF3 = KFIE
     *         + ICKATR(ISFI,ISYME)
C
         NFI   = MAX(NT1AM(ISFI),1)
         NLM   = MAX(NMATIJ(ISLM),1)
C
         CALL DGEMM('N','N',NT1AM(ISFI),NVIR(ISYME),NMATIJ(ISLM),
     *              ONE,WORK(KOFF1),NFI,WORK(KOFF2),NLM,
     *              ONE,WORK(KOFF3),NFI)
C
      END DO !  ISLM
C
      !Put KFIE(fi,e) to N1FEI(fe,i)^G
      ISFGEI = MULD2H(ISFIE,ISYMG)
      DO ISYMI = 1,NSYM
         ISFGE = MULD2H(ISFGEI,ISYMI)
         ISFE  = MULD2H(ISFGE,ISYMG)
         DO ISYME = 1,NSYM
            ISFG  = MULD2H(ISFGE,ISYME)
            ISYMF = MULD2H(ISFG,ISYMG)
            ISFI  = MULD2H(ISYMF,ISYMI)
            DO I = 1,NRHF(ISYMI)
               DO E = 1,NVIR(ISYME)
                  DO F = 1,NVIR(ISYMF)
C
                     KOFF1 = KFIE
     *                     + ICKATR(ISFI,ISYME)
     *                     + NT1AM(ISFI)*(E-1)
     *                     + IT1AM(ISYMF,ISYMI)
     *                     + NVIR(ISYMF)*(I-1)
     *                     + F-1
C
                     KOFF2 = IMAABI(ISFE,ISYMI)
     *                     + NMATAB(ISFE)*(I-1)
     *                     + IMATAB(ISYMF,ISYME)
     *                     + NVIR(ISYMF)*(E-1)
     *                     + F
C
                     N1FEI(KOFF2) = N1FEI(KOFF2) + WORK(KOFF1)
C
                  END DO !F
               END DO !E
            END DO !I
         END DO !ISYME
      END DO !ISYMI
C


*--------------------------------------------------------------------*
* Now we go on to calculate the contributions to N_anmo:             *
*                                                ======              *
*                                                                    *
* N_anmo = sum_dle (W^ae(dlon) + W^da(eonl) + W^de(anol)) t^{de}_{lm}*
*                                                                    *
* Stored as:                                                         *
*                                                                    *
*    N(anm,o)                                                        *
*                                                                    *
*--------------------------------------------------------------------*

C
C     ----------------------------------------------
C      Calculate the first contribution to N_anmo:
C      
C      sum_dle W^ae(dlon) * t^{de}_{lm}  --> N_anmo
C
C      calculated as:
C
C      WMAT^AE(dlo,n)       T2TP(dlm,E)
C           |                    |
C           | indsq(6)           | <-- resolve in loops
C           V                    V
C       KDLON(dl,on)  *   T2TP(dl,m) --> KTMP(on,m)
C
C       KTMP(on,m) --> N2MAT(Anm,o)
C     ----------------------------------------------

      A = IB
      E = ID
C
      ISYMA  = ISYMIB
      ISYME  = ISYMID
      ISYMAE = MULD2H(ISYMA,ISYME)
      ISDLON = MULD2H(ISWMAT,ISYMAE)
      ISDLM  = MULD2H(ISYMT2,ISYME)
      ISONM  = MULD2H(ISDLON,ISDLM)
C
      KDLON = 1
      KTMP  = KDLON + NCKIJ(ISDLON)
      KEND1 = KTMP  + NMAIJK(ISONM)
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available: ', LWORK
         WRITE(LUPRI,*) 'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in WT2_N1N2 (6)')
      END IF
C
      CALL DZERO(WORK(KTMP),NMAIJK(ISONM))
C
      !If symmetry, sort WMAT^AE(dlo,n) to KDLON(dl,on)
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(NCKIJ(ISDLON),WORK(KDLON),WMAT,INDSQ(1,6))
      ELSE
         CALL DCOPY(NCKIJ(ISDLON),WMAT,1,WORK(KDLON),1)
      END IF
C
      !Multiply KDLON(dl,on)   *   T2TP(dl,m) --> KTMP(on,m)
      DO ISDL = 1,NSYM
         ISON  = MULD2H(ISDLON,ISDL)
         ISYMM = MULD2H(ISDLM,ISDL)
C
         KOFF1 = KDLON
     *         + ISAIKL(ISDL,ISON)
         KOFF2 = IT2SP(ISDLM,ISYME)
     *         + NCKI(ISDLM)*(E-1)
     *         + ICKI(ISDL,ISYMM)
     *         + 1
         KOFF3 = KTMP
     *         + IMAIJK(ISON,ISYMM)
C
         NDL   = MAX(NT1AM(ISDL),1)
         NON   = MAX(NMATIJ(ISON),1)
C
         CALL DGEMM('T','N',NMATIJ(ISON),NRHF(ISYMM),NT1AM(ISDL),
     *              ONE,WORK(KOFF1),NDL,T2TP(KOFF2),NDL,
     *              ONE,WORK(KOFF3),NON)
C
      END DO ! ISDL
C
      !Put  KTMP(on,m) to N2MAT(Anm,o)
      ISANMO = MULD2H(ISONM,ISYMA)
      DO ISYMO = 1,NSYM
         ISANM = MULD2H(ISANMO,ISYMO)
         DO ISYMM = 1,NSYM
            ISAN  = MULD2H(ISANM,ISYMM)
            ISYMN = MULD2H(ISAN,ISYMA)
            ISON  = MULD2H(ISYMO,ISYMN)
            DO IO = 1,NRHF(ISYMO)
               DO M = 1,NRHF(ISYMM)
                  DO N = 1,NRHF(ISYMN)
C
                     KOFF1 = KTMP
     *                     + IMAIJK(ISON,ISYMM)
     *                     + NMATIJ(ISON)*(M-1)
     *                     + IMATIJ(ISYMO,ISYMN)
     *                     + NRHF(ISYMO)*(N-1)
     *                     + IO-1
                     KOFF2 = ISAIKJ(ISANM,ISYMO)
     *                     + NCKI(ISANM)*(IO-1)
     *                     + ISAIK(ISAN,ISYMM)
     *                     + NT1AM(ISAN)*(M-1)
     *                     + IT1AM(ISYMA,ISYMN)
     *                     + NVIR(ISYMA)*(N-1)
     *                     + A
C
                     N2MAT(KOFF2) = N2MAT(KOFF2) + WORK(KOFF1)
C
                  END DO !N
               END DO !M 
            END DO !IO 
         END DO !ISYMM 
      END DO !ISYMO 
C

C
C     ----------------------------------------------------
      IF (W3X) THEN !Calculate ALL contributions to N_anmo
C     ----------------------------------------------------
C

C
C        ----------------------------------------------
C         Calculate the second contribution to N_anmo:
C         
C         sum_dle W^da(eonl) * t^{ed}_{ml}  --> N_anmo
C
C         calculated as:
C
C         WMAT^DA(eon,l)       T2TP(eml,D)
C              |                    |
C              | indsq(2)           |                     
C              V                    |
C          KELON(elo,n)             | <-- CALL SORT_T2_AJI                 
C              |                    |
C              | indsq(6)           |                     
C              V                    V
C          KELON(el,on)  *   KELM(el,m) --> KTMP(on,m)
C
C          KTMP(on,m) --> N2MAT(Anm,o)
C        ----------------------------------------------

         D = IB
         A = ID
C
         ISYMD  = ISYMIB
         ISYMA  = ISYMID
         ISYMDA = MULD2H(ISYMD,ISYMA)
         ISEONL = MULD2H(ISWMAT,ISYMDA)
         ISELON = ISEONL
         ISELM  = MULD2H(ISYMT2,ISYMD)
         ISONM  = MULD2H(ISELON,ISELM)
C
         KELON = 1
         KELM  = KELON + NCKIJ(ISELON)
         KTMP  = KELM  + NCKI(ISELM)
         KEND1 = KTMP  + NMAIJK(ISONM)
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available: ', LWORK
            WRITE(LUPRI,*) 'Memory needed   : ', KEND1
            CALL QUIT('Insufficient memory in WT2_N1N2 (7)')
         END IF
C
         CALL DZERO(WORK(KTMP),NMAIJK(ISONM))
C
         !Sort WMAT^DA(eon,l) to KELON(elo,n)
         DO I = 1,NCKIJ(ISEONL)
            WORK(KELON+I-1) = WMAT(INDSQ(I,2))
         END DO
C
         !If symmetry sort KELON(elo,n)  to  KELON(el,on)
C
         IF (NSYM .GT. 1) THEN
           IF (LWRK1 .LT. NCKIJ(ISELON)) THEN
              WRITE(LUPRI,*) 'Memory available: ', LWRK1
              WRITE(LUPRI,*) 'Memory needed   : ', NCKIJ(ISELON)
              CALL QUIT('Insufficient memory in WT2_N1N2 (8)')
           END IF
           CALL CC_GATHER(NCKIJ(ISELON),WORK(KEND1),WORK(KELON),
     *                    INDSQ(1,6))
           CALL DCOPY(NCKIJ(ISELON),WORK(KEND1),1,WORK(KELON),1)
         END IF 
C
         !Sort T2TP(eml,D) to KELM(el,m)
         CALL SORT_T2_AJI(WORK(KELM),ISYMD,D,T2TP,ISYMT2)
C
         !Multiply KELON(el,on)   *   KELM(el,m) --> KTMP(on,m)
C
         DO ISEL = 1,NSYM
            ISON  = MULD2H(ISELON,ISEL)
            ISYMM = MULD2H(ISELM,ISEL)
C
            KOFF1 = KELON
     *            + ISAIKL(ISEL,ISON)
            KOFF2 = KELM
     *            + ICKI(ISEL,ISYMM)
            KOFF3 = KTMP
     *            + IMAIJK(ISON,ISYMM)
C
            NEL   = MAX(NT1AM(ISEL),1)
            NON   = MAX(NMATIJ(ISON),1)
C
            CALL DGEMM('T','N',NMATIJ(ISON),NRHF(ISYMM),NT1AM(ISEL),
     *                 ONE,WORK(KOFF1),NEL,WORK(KOFF2),NEL,
     *                 ONE,WORK(KOFF3),NON)
C
         END DO ! ISEL
C
         !Put  KTMP(on,m) to N2MAT(Anm,o)
         ISANMO = MULD2H(ISONM,ISYMA)
         DO ISYMO = 1,NSYM
            ISANM = MULD2H(ISANMO,ISYMO)
            DO ISYMM = 1,NSYM
               ISAN  = MULD2H(ISANM,ISYMM)
               ISYMN = MULD2H(ISAN,ISYMA)
               ISON  = MULD2H(ISYMO,ISYMN)
               DO IO = 1,NRHF(ISYMO)
                  DO M = 1,NRHF(ISYMM)
                     DO N = 1,NRHF(ISYMN)
C
                        KOFF1 = KTMP
     *                        + IMAIJK(ISON,ISYMM)
     *                        + NMATIJ(ISON)*(M-1)
     *                        + IMATIJ(ISYMO,ISYMN)
     *                        + NRHF(ISYMO)*(N-1)
     *                        + IO-1
                        KOFF2 = ISAIKJ(ISANM,ISYMO)
     *                        + NCKI(ISANM)*(IO-1)
     *                        + ISAIK(ISAN,ISYMM)
     *                        + NT1AM(ISAN)*(M-1)
     *                        + IT1AM(ISYMA,ISYMN)
     *                        + NVIR(ISYMA)*(N-1)
     *                        + A
C
                        N2MAT(KOFF2) = N2MAT(KOFF2) + WORK(KOFF1)
C
                     END DO !N
                  END DO !M 
               END DO !IO 
            END DO !ISYMM 
         END DO !ISYMO 
C
C        ----------------------------------------------
C         Calculate the third contribution to N_anmo:
C         
C         sum_dle W^de(anol) * t^{de}_{lm}  --> N_anmo
C
C         calculated as:
C
C         WMAT^DE(ano,l)       T2TP(Dlm,E)
C                                   |
C                                   | <-- CALL SORT_T2_IJ
C                                   V
C         WMAT^DE(ano,l)   *   KLM(l,m) --> KTMP(ano,m)
C
C
C                       indsq(3)
C          KTMP(ano,m) ----------> N2MAT(anm,o)
C        ----------------------------------------------

         D = IB
         E = ID
C
         ISYMD  = ISYMIB
         ISYME  = ISYMID
         ISYMDE = MULD2H(ISYMD,ISYME)
         ISANOL = MULD2H(ISWMAT,ISYMDE)
         ISLM   = MULD2H(ISYMT2,ISYMDE)
         ISANOM = MULD2H(ISANOL,ISLM)
         ISANMO = ISANOM
C
         KTMP  = 1
         KLM   = KTMP + NCKIJ(ISANOM)
         KEND1 = KLM  + NMATIJ(ISLM)
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available: ', LWORK
            WRITE(LUPRI,*) 'Memory needed   : ', KEND1
            CALL QUIT('Insufficient memory in WT2_N1N2 (9)')
         END IF
C
         CALL DZERO(WORK(KTMP),NCKIJ(ISANOM))
C
         !Sort T2TP(Dlm,E) to KLM(l,m)
         CALL SORT_T2_IJ(WORK(KLM),ISYMD,D,ISYME,E,T2TP,ISYMT2)
C
         !Multiply WMAT^DE(ano,l)   *   KLM(l,m) --> KTMP(ano,m)
C
         DO ISYML = 1,NSYM
            ISANO = MULD2H(ISANOL,ISYML)
            ISYMM = MULD2H(ISLM,ISYML)
C
            KOFF1 = ISAIKJ(ISANO,ISYML)
     *            + 1
            KOFF2 = KLM
     *            + IMATIJ(ISYML,ISYMM)
            KOFF3 = KTMP
     *            + ISAIKJ(ISANO,ISYMM)
C
            NANO = MAX(NCKI(ISANO),1)
            NL   = MAX(NRHF(ISYML),1)
C
            CALL DGEMM('N','N',NCKI(ISANO),NRHF(ISYMM),NRHF(ISYML),
     *                 ONE,WMAT(KOFF1),NANO,WORK(KOFF2),NL,
     *                 ONE,WORK(KOFF3),NANO)
C
         END DO ! ISYML
C
         !Put KTMP(ano,m) to N2MAT(anm,o)
         DO I = 1,NCKIJ(ISANOM)
            N2MAT(I) = N2MAT(I) + WORK(INDSQN(I,3))
         END DO
C
C     ----------------------------
      END IF !W3X
C     ----------------------------
C

C
C----------
C     End.
C----------
C
      CALL QEXIT('WTN')
C
      RETURN
      END
C
C  /* Deck n1_resort */
      SUBROUTINE N1_RESORT(ISYMN1,LUN1,FNN1,LUGEI,FNGEI,LUFEI,FNFEI,
     *                     WORK,LWORK,SKIPGEI)
*
**********************************************************************
*                                                                    *
*     Read in (gei,F) array  from file LUGEI and (fei,G) array  from *
*     file LUFEI, sum them up and put to LUN1 file as (fge,i) array  *
*                                                                    *
*     If SKIPGEI = .true. then skip the contribution from LUGEI file.*
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
* F. Pawlowski, 20-Apr-2004, Aarhus.                                 *
**********************************************************************
*
      IMPLICIT NONE
C
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
C
      LOGICAL SKIPGEI
C
      CHARACTER*(*) FNGEI,FNFEI
      CHARACTER*(*) FNN1
      INTEGER   LUGEI,LUFEI,LUN1
C
      INTEGER ISYMN1,LWORK
      INTEGER ISYMI,ISFGE,KFGE,KEND1,LWRK1
      INTEGER ISYMF,ISGE,KGE,KEND2,LWRK2
      INTEGER IADR,ISGEI,ISYME,ISYMG,ISFG,KOFF1,KOFF2
C
      INTEGER ISFE,ISFEI,KFE
C
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('N1RSRT')
C
      DO ISYMI = 1,NSYM !"I" is general loop index in this routine
         ISFGE = MULD2H(ISYMN1,ISYMI)
C
         KFGE  = 1
         KEND1 = KFGE + NMAABC(ISFGE)
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available: ',LWORK
            WRITE(LUPRI,*)'Memory needed:    ',KEND1
            CALL QUIT('Insufficient memory in N1_RESORT (1)')
         END IF
C
         DO I = 1,NRHF(ISYMI)
C
            CALL DZERO(WORK(KFGE),NMAABC(ISFGE))
C
C           ----------------------
            IF (.NOT.SKIPGEI) THEN
C           ----------------------
C

C
C              ----------------------------------------------------
C              Read in (gei,F) from file LUGEI and put to KFGE(fge) 
C              array for fixed I
C              ----------------------------------------------------
C
               DO ISYMF = 1,NSYM
                  ISGE  = MULD2H(ISFGE,ISYMF)
                  ISGEI = MULD2H(ISGE,ISYMI)
C
                  KGE   = KEND1
                  KEND2 = KGE + NMATAB(ISGE)
                  LWRK2 = LWORK - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     WRITE(LUPRI,*)'Memory available: ',LWORK
                     WRITE(LUPRI,*)'Memory needed:    ',KEND2
                     CALL QUIT('Insufficient memory in N1_RESORT (2)')
                  END IF
C
                  DO F = 1,NVIR(ISYMF)
C
                     CALL DZERO(WORK(KGE),NMATAB(ISGE))
C
                     !Read in (ge)^FI
                     IADR = ICKBD(ISGEI,ISYMF) + NCKATR(ISGEI)*(F-1) 
     *                    + IMAABI(ISGE,ISYMI) + NMATAB(ISGE)*(I-1) + 1
                     CALL GETWA2(LUGEI,FNGEI,WORK(KGE),IADR,
     *                           NMATAB(ISGE))
C
                     !Sort (ge)^FI to (fge)^I
                     DO ISYME = 1,NSYM
                        ISYMG = MULD2H(ISGE,ISYME)
                        ISFG  = MULD2H(ISYMF,ISYMG)
C
                        DO E = 1,NVIR(ISYME)
                           DO G = 1,NVIR(ISYMG)
                              KOFF1 = KGE
     *                              + IMATAB(ISYMG,ISYME)
     *                              + NVIR(ISYMG)*(E-1)
     *                              + G-1
                              KOFF2 = KFGE
     *                              + IMAABC(ISFG,ISYME)
     *                              + NMATAB(ISFG)*(E-1)
     *                              + IMATAB(ISYMF,ISYMG)
     *                              + NVIR(ISYMF)*(G-1)
     *                              + F-1
C
                              WORK(KOFF2) = WORK(KOFF2) + WORK(KOFF1)
C
                           END DO !G
                        END DO !E
                     END DO !ISYME
C
                  END DO !F
               END DO !ISYMF
C
C           ----------------------
            END IF !.NOT.SKIPGEI
C           ----------------------
C

C
C           ----------------------------------------------------
C           Read in (fei,G) from file LUFEI and put to KFGE(fge) 
C           array for fixed I
C           ----------------------------------------------------
C
            DO ISYMG = 1,NSYM
               ISFE  = MULD2H(ISFGE,ISYMG)
               ISFEI = MULD2H(ISFE,ISYMI)
C
               KFE   = KEND1
               KEND2 = KFE + NMATAB(ISFE)
               LWRK2 = LWORK - KEND2
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*)'Memory available: ',LWORK
                  WRITE(LUPRI,*)'Memory needed:    ',KEND2
                  CALL QUIT('Insufficient memory in N1_RESORT (3)')
               END IF
C
               DO G = 1,NVIR(ISYMG)
C
                  CALL DZERO(WORK(KFE),NMATAB(ISFE))
C
                  !Read in (fe)^GI
                  IADR = ICKBD(ISFEI,ISYMG) + NCKATR(ISFEI)*(G-1) 
     *                 + IMAABI(ISFE,ISYMI) + NMATAB(ISFE)*(I-1) + 1
                  CALL GETWA2(LUFEI,FNFEI,WORK(KFE),IADR,NMATAB(ISFE))
C
                  !Sort (fe)^GI to (fge)^I
                  DO ISYME = 1,NSYM
                     ISYMF = MULD2H(ISFE,ISYME)
                     ISFG  = MULD2H(ISYMF,ISYMG)
C
                     DO E = 1,NVIR(ISYME)
                        DO F = 1,NVIR(ISYMF)
                           KOFF1 = KFE
     *                           + IMATAB(ISYMF,ISYME)
     *                           + NVIR(ISYMF)*(E-1)
     *                           + F-1
                           KOFF2 = KFGE
     *                           + IMAABC(ISFG,ISYME)
     *                           + NMATAB(ISFG)*(E-1)
     *                           + IMATAB(ISYMF,ISYMG)
     *                           + NVIR(ISYMF)*(G-1)
     *                           + F-1
C
                           WORK(KOFF2) = WORK(KOFF2) + WORK(KOFF1)
C
                        END DO !F
                     END DO !E
                  END DO !ISYME
C
               END DO !G
            END DO !ISYMG
C
                  
C
C           -----------------------------------
C           Put (fge)^I to LUN1 file as (fge,i)
C           -----------------------------------
C
            IADR = IMAABCI(ISFGE,ISYMI) + NMAABC(ISFGE)*(I-1) + 1
            CALL PUTWA2(LUN1,FNN1,WORK(KFGE),IADR,NMAABC(ISFGE))
C
C        -----------------
C        End general loop.
C        -----------------
C
         END DO !I
      END DO !ISYMI
C
C----------
C     End.
C----------
C
      CALL QEXIT('N1RSRT')
C
      RETURN
      END

