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
      REAL SECOND
#else
      DOUBLE PRECISION OMEGA1(*), T1AM(*), T2TP(*), C2TP(*)
      DOUBLE PRECISION XLAMDH(*), XLAMDP(*), WORK(LWORK)
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      DOUBLE PRECISION DTIME, XNORM, DDOT, ECURR
#if !defined (VAR_G77)
      DOUBLE PRECISION SECOND
#endif
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
      LUDELD = 0
      LUDKBC = 0
      FNDELD = 'CC3_T3_TMP1'
      FNDKBC = 'CC3_T3_TMP2'
C
      CALL WOPEN2(LUDELD,FNDELD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
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
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
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
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
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
      DTIME  = SECOND() - DTIME
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
      NUMFLM = MAX(1,NCKI(ISYFLM))
C
      CALL DGEMV('T',NCKI(ISYFLM),NRHF(ISYMN),ONE,
     *           TMAT(KOFF1),NUMFLM,C2TP(KOFF2),1,
     *           ZERO,WORK(KOFF3),1)
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
      INTEGER ISYML1, ISYML2, ISYMT2, ISYRES, LWORK
      INTEGER ISYMTR, ISINT1, ISINT2, ISYMIM, KFOCKD, KCMO, KFCKBA
      INTEGER KEND0, LWRK0, KEND1, LWRK1, KEND2, LWRK2, ISYMC, ISYMK
      INTEGER KOFF1, KOFF2, KTROC0, KTROC2, KXIAJB, KINTOC, LUFCK
      INTEGER LENGTH, ISYOPE, IOPTTCME, IOFF, ISYMD, ISAIJ1, ISYCKB
      INTEGER ISCKB2, KTRVI0, KTRVI2, KTRVI3, KTRVI4, KTRVI5, KTRVI6
      INTEGER KRMAT1, KEND3, LWRK3, KEND4, LWRK4, KINTVI
      INTEGER ISYMB, ISYALJ, ISAIJ2, ISYMBD, ISCKIJ, KSMAT, KQMAT
      INTEGER KDIAG, LENSQ, KINDSQ, KINDEX, KTMAT, KRMAT2
      INTEGER ISCKB1, KTRVI, KTRVI1, KTROC, KTROC1, KVVVV
      INTEGER LUDKBC4, LUCKJD
      INTEGER LUDKBC, LUTOC, LU3VI, LU4V, LUDKBC3, LU3FOPX, LU3FOP2X
C
#if defined (SYS_CRAY)
      REAL L1AM(*), L2TP(*), T2TP(*), OMEGA1(*), OMEGA2(*), XINT4O(*)
      REAL XOVVO(*), XOOVV(*), XLAMDP(*), XLAMDH(*), WORK(LWORK)
      REAL TITRAN, TISORT, TISMAT, TIQMAT, TICONT, TIOME1
      REAL DTIME, XL2TP, XIAJB, XINT, XTROC0, XTRVI, XSMAT
      REAL RHO1N, RHO2N, XDIA, DDOT, HALF, ONE, ECURR
      REAL SECOND
#else
      DOUBLE PRECISION L1AM(*), L2TP(*), T2TP(*)
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), XINT4O(*)
      DOUBLE PRECISION XOVVO(*),XOOVV(*),XLAMDP(*),XLAMDH(*),WORK(LWORK)
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TIQMAT, TICONT, TIOME1
      DOUBLE PRECISION DTIME, XL2TP, XIAJB, XINT, XTROC0, XTRVI, XSMAT
      DOUBLE PRECISION RHO1N, RHO2N, XDIA, DDOT, HALF, ONE, ECURR
#if !defined (VAR_G77)
      DOUBLE PRECISION SECOND
#endif
#endif
C
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <iratdef.h>
#include <ccinftap.h>
C
      CHARACTER*(*) FNCKJD, FNDKBC, FNTOC, FN3VI, FN4V
      CHARACTER*(*) FNDKBC3, FN3FOPX, FN3FOP2X
      CHARACTER*11 FNDKBC4
      CHARACTER*1 CDUMMY
C
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('CC3_L3_LHTR')
C
C----------------------------------------------------
C     Initialise character strings and open files
C----------------------------------------------------
C
      CDUMMY = ' '
      LUDKBC4 = 0
      FNDKBC4 = 'CC3_L3_TMP1'
C
      CALL WOPEN2(LUDKBC4,FNDKBC4,64,0)
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
      TISMAT = 0.0D0
      TIQMAT = 0.0D0
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
      KFOCKD = 1
      KCMO   = KFOCKD + NORBTS
      KFCKBA = KCMO   + NLAMDS
      KEND0  = KFCKBA + N2BST(ISYMOP)
      LWRK0  = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_LHTR_L3')
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
      LUFCK = 0
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
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC  = KEND0
      KTROC0 = KTROC  + NTRAOC(ISINT2)
      KTROC1 = KTROC0 + NTRAOC(ISINT2)
      KTROC2 = KTROC1 + NTRAOC(ISINT2)
      KXIAJB = KTROC2 + NTRAOC(ISINT2)
      KEND1  = KXIAJB + NT2AM(ISYMOP)
      LWRK1  = LWORK  - KEND1
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
C------------------------------------------------------------
C     Read in integrals used in L3 and transform.
C------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),XLAMDH,
     *                 WORK(KEND2),LWRK2,ISINT2)
C
C-----------------------------------------------------------
C     Construct 2*C-E of the integrals.
C-----------------------------------------------------------
C
      CALL CCSDT_TCMEOCC(WORK(KTROC0),WORK(KTROC2),ISINT2)
C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XINT
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 CC3_L3_LHTR : ',XTROC0
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC2),1,
     *                WORK(KTROC2),1)
         WRITE(LUPRI,*) 'Norm of TROC2 CC3_L3_LHTR : ',XTROC0
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
         ISCKB1 = MULD2H(ISINT1,ISYMD)
         ISCKB2 = MULD2H(ISINT2,ISYMD)
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
         KTRVI2 = KTRVI1 + NCKATR(ISCKB1)
         KRMAT1 = KTRVI2 + NCKATR(ISCKB2)
         KEND2  = KRMAT1 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0  = KEND2
         KTRVI3  = KTRVI0  + NCKATR(ISCKB2)
         KTRVI4  = KTRVI3  + NCKATR(ISCKB2)
         KTRVI5  = KTRVI4  + NCKATR(ISCKB2)
         KTRVI6  = KTRVI5  + NCKATR(ISCKB2)
         KVVVV   = KTRVI6  + NCKATR(ISCKB2)
         KEND3   = KVVVV   + NMAABC(ISCKB2)
         LWRK3   = LWORK  - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISYMD),NCKA(ISCKB2))
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
C
C--------------------------------------------
C           Read in g_{vvvv} for a given D
C--------------------------------------------
C
            IF (NMAABC(ISCKB2) .GT. 0) THEN
               IOFF = I3VVIR(ISCKB2,ISYMD)
     *              + NMAABC(ISCKB2)*(D-1)
     *              + 1
               CALL GETWA2(LU4V,FN4V,WORK(KVVVV),IOFF,NMAABC(ISCKB2))
            ENDIF
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
C-----------------------------------------------
C           Integrals used in L3am.
C-----------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI0),IOFF,
     &                     NCKATR(ISCKB2))
            ENDIF
C
            CALL CCSDT_SRVIR3(WORK(KTRVI0),WORK(KEND4),ISYMD,D,ISINT2)
C
C------------------------------------------------------
C           Read 2*C-E of integral used for L3
C------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI4),IOFF,
     &                     NCKATR(ISCKB2))
            ENDIF
C
C-----------------------------------------------------------
C           Sort the integrals for s3am and for t3-bar
C-----------------------------------------------------------
C
            DTIME = SECOND()
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            CALL CCSDT_SRTVIR(WORK(KTRVI4),WORK(KTRVI5),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55) THEN
               XTRVI = DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 ',XTRVI
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI = DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XTRVI
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI = DDOT(NCKATR(ISCKB2),WORK(KTRVI4),1,
     *                      WORK(KTRVI4),1)
               WRITE(LUPRI,*) 'Norm of TRVI4 ',XTRVI
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI = DDOT(NCKATR(ISCKB2),WORK(KTRVI5),1,
     *                      WORK(KTRVI5),1)
               WRITE(LUPRI,*) 'Norm of TRVI5 ',XTRVI
            ENDIF
C
C------------------------------------------------------
C           Calculate virtual integrals used in q3am.
C------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKB))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),XLAMDP,
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_L3_LHTR  (TRVI3)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND3),ISYMD,D,ISINT2)
C
C---------------------------------------------------------------
C           Read virtual integrals used in q3am/u3am for t3-bar.
C---------------------------------------------------------------
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI6),IOFF,
     &                     NCKATR(ISYCKB))
            ENDIF
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_L3_LHTR (TRVI6)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI6),WORK(KEND4),ISYMD,D,ISINT2)
C
            IF (IPRINT .GT. 55) THEN
               XTRVI = DDOT(NCKATR(ISCKB2),WORK(KTRVI6),1,
     *                      WORK(KTRVI6),1)
               WRITE(LUPRI,*) 'Norm of TRVI6 ',XTRVI
            ENDIF
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
               KSMAT   = KEND3
               KQMAT   = KSMAT   + NCKIJ(ISCKIJ)
               KDIAG   = KQMAT   + NCKIJ(ISCKIJ)
               KINDSQ  = KDIAG   + NCKIJ(ISCKIJ)
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KTMAT   = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KRMAT2  = KTMAT   + NCKIJ(ISCKIJ)
               KEND4   = KRMAT2  + NCKI(ISAIJ2)
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
               CALL CC3_INDEX(WORK(KINDEX),ISYALJ)
C
               DO B = 1,NVIR(ISYMB)
C
C-----------------------------------------
C                 Initialize the R2 matrix.
C-----------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
C
C---------------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for for B,D for L3
C                 Scale the smat with a half since have a factor
C                 of two in the routine.
C---------------------------------------------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DZERO(WORK(KSMAT),NCKIJ(ISCKIJ))
C
                  CALL CCFOP_SMAT(ECURR,L1AM,ISYML1,L2TP,ISYML2,
     *                            WORK(KTMAT),
     *                            WORK(KFCKBA),WORK(KXIAJB),ISINT1,
     *                            WORK(KTRVI0),WORK(KTRVI2),
     *                            WORK(KTRVI4),WORK(KTRVI5),
     *                            WORK(KTROC0),WORK(KTROC2),
     *                            ISINT2,WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KSMAT),WORK(KEND4),LWRK4,
     *                            WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                            ISYMB,B,ISYMD,D)
C
                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT),1)
C
                  CALL T3_FORBIDDEN(WORK(KSMAT),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT-L3    ',XSMAT
                  ENDIF
C
C-------------------------------------------------------------------
C                 Calculate Q(ci,jk) for fixed b,d for t3-bar.
C-------------------------------------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DZERO(WORK(KQMAT),NCKIJ(ISCKIJ))
C
                  CALL CCFOP_QMAT(ECURR,L1AM,ISYML1,L2TP,ISYML2,
     *                            WORK(KTMAT),WORK(KFCKBA),
     *                            WORK(KXIAJB),ISINT1,WORK(KTRVI3),
     *                            WORK(KTRVI6),WORK(KTROC0),
     *                            WORK(KTROC2),ISINT2,WORK(KFOCKD),
     *                            WORK(KDIAG),WORK(KQMAT),
     *                            WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                            LENSQ,ISYMB,B,ISYMD,D)
C
                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KQMAT),1)
C
                  CALL T3_FORBIDDEN(WORK(KQMAT),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT-L3  ',XSMAT
                  ENDIF
C
C-----------------------------------------------------------------------
C                 Calculate the contributions to omega2
C-----------------------------------------------------------------------
C
                  CALL CCFOP_CONVIR(WORK(KRMAT2),WORK(KSMAT),
     *                              WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                              WORK(KTRVI),WORK(KTRVI1),ISINT1,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
C
                  CALL CCFOP_CONOCC(OMEGA2,WORK(KRMAT1),
     *                              WORK(KRMAT2),WORK(KSMAT),
     *                              WORK(KTMAT),ISYMIM,
     *                              WORK(KTROC),WORK(KTROC1),ISINT1,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
C
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KSMAT),1)
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KQMAT),1)
C
C------------------------------------------------------------------------
C                 Calculate the L3 contribution to omega1
C------------------------------------------------------------------------
C
                  CALL CC3_L3_OMEGA1(OMEGA1,ISYRES,WORK(KSMAT),
     *                               WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                               XINT4O,XOVVO,XOOVV,WORK(KVVVV),1,
     *                               T2TP,ISYMT2,WORK(KEND4),LWRK4,
     *                               LENSQ,WORK(KINDSQ),ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO2N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm (Rho1) after CC3_L3_OMEGA1',
     *                               RHO2N
                  ENDIF
C
C-------------------------------------------
C                 Accumulate R2 in Omega2
C-------------------------------------------
C
                  CALL CC3_RACC(OMEGA2,WORK(KRMAT2),ISYMB,B,
     *                          ISYRES)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
                     WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC',RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_RACC: ')
                     CALL CC_PRP(DUMMY,OMEGA2,ISYRES,0,1)
                  ENDIF
C
               ENDDO   ! B
            ENDDO      ! ISYMB
C
C---------------------------------------
C           Accumulate R1 in Omega2.
C---------------------------------------
C
            CALL CC3_RACC(OMEGA2,WORK(KRMAT1),ISYMD,D,ISYRES)
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
         ENDDO       ! D
      ENDDO          ! ISYMD
C
C-------------------------------
C     Close and delete files
C-------------------------------
C
      CALL WCLOSE2(LUDKBC4,FNDKBC4,'DELETE')
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
         LUTMP1 = 0
         LUTMP2 = 0
         LUTMP3 = 0
         LUTMP4 = 0
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
