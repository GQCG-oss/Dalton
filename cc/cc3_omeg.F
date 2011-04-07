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
C  /* Deck cc3_omeg */
      SUBROUTINE CC3_OMEG(ECURR,OMEGA1,OMEGA2,T1AM,ISYMT1,T2TP,ISYMT2,
     *                    FOCK,XLAMDP,XLAMDH,WORK,LWORK,LU3SRT,FN3SRT,
     *                    LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,FNDKBC,
     *                    LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Add the triples contribution to the omega vectors
C
C     Ove Christiansen  Jan. 1996:
C
C         Generalisation to non-total sym. case for CC3 rsp.
C
C         ISYMT2 is symmetry of T2TP
C         ISYMT1 is symmetry of T1AM input when CC3LR
C         else ISYMT1 is symmetry of the C1 that the
C         triples integrals have been transformed with.
C
C         Isyres = isymt1*isymt2*isymop
C
#include <implicit.h>
C
      DIMENSION OMEGA1(*),OMEGA2(*),T1AM(*),T2TP(*),
     *          FOCK(*),XLAMDP(*),XLAMDH(*),WORK(LWORK)
C
#include <iratdef.h>
#include <dummy.h>
#include <inftap.h>
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CHARACTER*(*) FN3SRT, FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FN3VI2
      CHARACTER*1 CDUMMY
C
      LOGICAL  FIRST,SECON
      SAVE     FIRST,SECON
      DATA     FIRST,SECON /.TRUE.,.FALSE./
C
      CALL QENTER('CC3_OMEG')
C
      CDUMMY = ' '
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
      ISYMTR = MULD2H(ISYMT1,ISYMT2)
      ISYRES = MULD2H(ISYMTR,ISYMOP)
      ISINT1 = ISYMOP
      ISINT2 = MULD2H(ISYMT1,ISYMOP)
      ISYMIM = MULD2H(ISYMTR,ISYMOP)
      IF (CC3LR) THEN
         ISINT1 = MULD2H(ISYMT1,ISYMOP)
         ISINT2 = ISYMOP
         ISYMIM = ISYMOP
      ENDIF
C
      IF (IPRINT .GT. 20 ) THEN
         WRITE(LUPRI,*) ' In CC3_OMEG: CC1A  = ',CC1A
         WRITE(LUPRI,*) ' In CC3_OMEG: CC1B  = ',CC1B
         WRITE(LUPRI,*) ' In CC3_OMEG: CC3   = ',CC3
         WRITE(LUPRI,*) ' In CC3_OMEG: CC3LR = ',CC3LR
         WRITE(LUPRI,*) ' In CC3_OMEG: ISYMT1, ISYMT2:',ISYMT1,ISYMT2
         WRITE(LUPRI,*) ' In CC3_OMEG: ISYRES, ISYMOP:',ISYRES,ISYMOP
         WRITE(LUPRI,*) ' In CC3_OMEG: ISINT1, ISINT2:',ISINT1,ISINT2
      ENDIF
C
C--------------------
C     Time variables.
C--------------------
C
      TITRAN = 0.0D0
      TISORT = 0.0D0
      TISMAT = 0.0D0
      TIQMAT = 0.0D0
      TICONO = 0.0D0
      TICONV = 0.0D0
      TIOME1 = 0.0D0
C
C---------------------------------------------------------
C     Transform and sort qmat integrals to smat integrals.
C---------------------------------------------------------
C
      IF (.NOT. CC1B) THEN
         CALL CC3_SORT1(WORK,LWORK,2,ISINT2,LU3SRT,FN3SRT,
     *                  LUDELD,FNDELD,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                  IDUMMY,CDUMMY)
         CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT2,LUDELD,FNDELD,
     *                 LUDKBC,FNDKBC)
      ENDIF
C
C--------------------------------------
C     Reorder the t2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP in CC3_OMEG')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMT2),T2TP,1,WORK,1)
      CALL CC3_T2TP(T2TP,WORK,ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XT2TP = DDOT(NT2SQ(ISYMT2),T2TP,1,T2TP,1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XT2TP
      ENDIF
C
C---------------------------------------------------------
C     Read canonical orbital energies and MO coefficients.
C---------------------------------------------------------
C
      KFOCKD = 1
      KCMO   = KFOCKD + NORBTS
      KFCKAK = KCMO   + NLAMDS
      KEND0  = KFCKAK + NT1AM(ISINT1)
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
C------------------------------
C     Resort Fock matrix F(kc).
C------------------------------
C
      IF (.NOT. CC3LR) THEN
C
         DO 50 ISYMC = 1,NSYM
C
            ISYMK = MULD2H(ISYMC,ISYMOP)
C
            DO 60 K = 1,NRHF(ISYMK)
C
               DO 70 C = 1,NVIR(ISYMC)
C
                  KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
                  KOFF2 = KFCKAK + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K - 1) + C - 1
C
                  WORK(KOFF2) = FOCK(KOFF1)
C
   70          CONTINUE
   60       CONTINUE
   50    CONTINUE
C
         IF (CC1A) THEN
            CALL DZERO(WORK(KFCKAK),NT1AM(ISYMOP))
         ENDIF
C
      ENDIF
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC  = KEND0
      KTROC1 = KTROC  + NTRAOC(ISINT1)
      KTROC0 = KTROC1 + NTRAOC(ISINT1)
      KXIAJB = KTROC0 + NTRAOC(ISINT2)
      KEND1  = KXIAJB + NT2AM(ISYMOP)
      LWRK1  = LWORK  - KEND1
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
C-----------------------------------------
C     Calculate Fock matrix used in CC3LR.
C-----------------------------------------
C
      IF (CC3LR) THEN
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
         CALL CC3LR_MFOCK(WORK(KFCKAK),T1AM,WORK(KXIAJB),ISYMT1)
C
         IF (IPRINT .GT. 55) THEN
            XFD   = DDOT(NT1AM(ISINT1),WORK(KFCKAK),1,
     *              WORK(KFCKAK),1)
            WRITE(LUPRI,*) 'Norm of MFOCK',XFD
         ENDIF
C
      ENDIF
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
      IF (.NOT. CC3LR) THEN
         ISYOPE = ISYMOP
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYOPE,IOPTTCME)
      ENDIF
C
      IF ( IPRINT .GT. 55) THEN
         XIAJB = DDOT(NT2AM(ISYMOP),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XIAJB
      ENDIF
C
C------------------------
C     Occupied integrals.
C------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYMOP) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
      ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XINT
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      DTIME = SECOND()
      IF (.NOT. CC3LR) THEN
         IF (CC1A) THEN
            CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC),WORK(KCMO),
     *                       WORK(KEND2),LWRK2)
         ELSE
            CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC),XLAMDH,
     *                       WORK(KEND2),LWRK2)
         END IF
      ELSE
C
         CALL DZERO(WORK(KTROC),NTRAOC(ISINT1))
         CALL CC3LR_QINT(WORK(KTROC),T1AM,WORK(KXIAJB),WORK(KEND2),
     *                   LWRK2,ISYMT1)
C
         IF (IPRINT .GT. 55) THEN
            XTROC = DDOT(NTRAOC(ISINT1),WORK(KTROC),1,
     *                   WORK(KTROC),1)
            WRITE(LUPRI,*) 'Norm of TROC after QINT',XTROC
         ENDIF
C
         IF (LWRK2 .LT. NTRAOC(ISINT1)) THEN
            CALL QUIT('Insufficient space in CC3_OMEG')
         END IF
C
         CALL CCSDT_SRTOCC(WORK(KTROC),WORK(KEND2),ISINT1)
         CALL DCOPY(NTRAOC(ISINT1),WORK(KEND2),1,WORK(KTROC),1)
      ENDIF
C
C-----------------------
C     Read in integrals.
C-----------------------
C
      IF (.NOT. CC1B) THEN
         IOFF = 1
         IF (NTOTOC(ISINT2) .GT. 0) THEN
            CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
         ENDIF
C
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      IF (.NOT. CC1B) THEN
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KCMO),
     *                    WORK(KEND2),LWRK2,ISINT2)
      ELSE
         CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KCMO),
     *                    WORK(KEND2),LWRK2)
      ENDIF
C     
      DTIME  = SECOND() - DTIME
      TITRAN = TITRAN   + DTIME

C
      DTIME = SECOND()
C
      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1,
     *                  WORK(KEND2),LWRK2)
C
      DTIME  = SECOND() - DTIME
      TISORT = TISORT   + DTIME
C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTRAOC(ISINT1),WORK(KTROC),1,
     *                WORK(KTROC),1)
         WRITE(LUPRI,*) 'Norm of TROC ',XTROC
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XINT
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC1 = DDOT(NTRAOC(ISINT1),WORK(KTROC1),1,
     *                WORK(KTROC1),1)
         WRITE(LUPRI,*) 'Norm of TROC1 ',XTROC1
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 ',XTROC0
      ENDIF
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO 100 ISYMD = 1,NSYM
C
         ISAIJ1 = MULD2H(ISYMD,ISYRES)
         ISYCKB = MULD2H(ISYMD,ISYMOP)
         ISCKB1 = MULD2H(ISINT1,ISYMD)
         ISCKB2 = MULD2H(ISINT2,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_OMEG: ISAIJ1:',ISAIJ1
            WRITE(LUPRI,*) 'In CC3_OMEG: ISYCKB:',ISYCKB
            WRITE(LUPRI,*) 'In CC3_OMEG: ISCKB1:',ISCKB1
            WRITE(LUPRI,*) 'In CC3_OMEG: ISCKB2:',ISCKB2
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
         KTRVI0 = KEND2
         KTRVI3 = KTRVI0 + NCKATR(ISCKB2)
         KEND3  = KTRVI3 + NCKATR(ISCKB2)
         LWRK3  = LWORK  - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISYMD),NCKA(ISCKB2))
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CCSDT_OMEG')
         END IF
C
         DO 110 D = 1,NVIR(ISYMD)
C
            IF (IPRINT .GT. 5) THEN
C
               IF (SECON) THEN
C
                  SECON = .FALSE.
C
                  XTIME = SECOND() - TIMED
                  XTOT  = DFLOAT(NVIRT)*XTIME
C
                  CALL AROUND('Time estimate for Triples')
C
                  WRITE(LUPRI,'(/21X,A,F10.2)')
     *                  'Time for the first d : ',XTIME
                  WRITE(LUPRI,'(21X,A,F10.2)')
     *                  'Total time estimate  : ',XTOT
                  WRITE(LUPRI,'(21X,A,I10)')
     *                  'Total memory usage(w): ',INFWRK
C
                  CALL FLSHFO(LUPRI)
C
               ENDIF
C
               IF (FIRST) THEN
                  FIRST = .FALSE.
                  SECON = .TRUE.
                  TIMED = SECOND()
               ENDIF
C
            ENDIF
C
C------------------------------------
C           Initialize the R1 matrix.
C------------------------------------
C
            CALL DZERO(WORK(KRMAT1),NCKI(ISAIJ1))
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction.
C------------------------------------------------------------
C
            IF (.NOT. CC3LR) THEN
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
               ENDIF
C
               IF (CC1A) THEN
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),WORK(KCMO),
     *                             ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
               ELSE
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),XLAMDP,
     *                             ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
               ENDIF
            ENDIF
C
C-----------------------------------------------
C           Read virtual integrals used in s3am.
C-----------------------------------------------
C
            IF (.NOT. CC1B) THEN
C
               IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
               IF (NCKATR(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                        NCKATR(ISCKB2))
               ENDIF
C
            ELSE
C
               IF (CC3LR) THEN
                  IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
                  IF (NCKA(ISYCKB) .GT. 0) THEN
                     CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                           NCKA(ISYCKB))
                  ENDIF
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI0),WORK(KCMO),
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
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
               XTRVI0= DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 ',XTRVI0
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI2= DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XTRVI2
            ENDIF
C
C------------------------------------------------------
C           Read virtual integrals used in contraction.
C------------------------------------------------------
C
            IF (.NOT. CC3LR) THEN
C
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
               ENDIF
               IF (CC1A) THEN
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KCMO),
     *                             ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
               ELSE
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMDP,
     *                             ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
               ENDIF
C
            ENDIF
C
C
C-----------------------------------------------
C           Read virtual integrals used in q3am.
C-----------------------------------------------
C
            IF (.NOT. CC1B) THEN
C
               IOFF = ICKAD(ISCKB2,ISYMD) + NCKA(ISCKB2)*(D - 1) + 1
               IF (NCKA(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     &                        NCKA(ISCKB2))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),XLAMDH,
     *                          ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
            ELSE
C
               IF (CC3LR) THEN
                  IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
                  IF (NCKA(ISYCKB) .GT. 0) THEN
                     CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                           NCKA(ISYCKB))
                  ENDIF
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),WORK(KCMO),
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            ENDIF
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDT_OMEG')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND3),ISYMD,D,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55) THEN
               XTRVI3= DDOT(NCKATR(ISCKB2),WORK(KTRVI3),1,
     *                      WORK(KTRVI3),1)
               WRITE(LUPRI,*) 'Norm of TRVI3 ',XTRVI3
            ENDIF
C
C---------------------------------------------
C           Construct integrals used in CC3LR.
C---------------------------------------------
C
            IF (CC3LR) THEN
               CALL CC3LR_PINT(WORK(KTRVI),WORK(KTRVI1),T1AM,
     *                         WORK(KXIAJB),WORK(KEND4),LWRK4,
     *                         ISYMD,D,ISYMT1)
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI= DDOT(NCKATR(ISCKB1),WORK(KTRVI),1,
     *                      WORK(KTRVI),1)
               WRITE(LUPRI,*) 'Norm of TRVI ',XTRVI
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI1= DDOT(NCKATR(ISCKB1),WORK(KTRVI1),1,
     *                      WORK(KTRVI1),1)
               WRITE(LUPRI,*) 'Norm of TRVI1 ',XTRVI1
            ENDIF
C
C---------------------
C           Calculate.
C---------------------
C
            DO 120 ISYMB = 1,NSYM
C
               ISYALJ = MULD2H(ISYMB,ISYMT2)
               ISAIJ2 = MULD2H(ISYMB,ISYRES)
               ISYMBD = MULD2H(ISYMB,ISYMD)
               ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
               IF (IPRINT .GT. 55) THEN
C
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISCKIJ:',ISCKIJ
C
               ENDIF
C
               KSMAT  = KEND3
               KQMAT  = KSMAT  + NCKIJ(ISCKIJ)
               KDIAG  = KQMAT  + NCKIJ(ISCKIJ)
               KINDSQ = KDIAG  + NCKIJ(ISCKIJ)
               KINDEX = KINDSQ + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KTMAT  = KINDEX + (NCKI(ISYALJ) - 1)/IRAT + 1
               KRMAT2 = KTMAT  + NCKIJ(ISCKIJ)
               KEND4  = KRMAT2 + NCKI(ISAIJ2)
               LWRK4  = LWORK  - KEND4
C
               INFWRK = KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    > ',KEND4
                  CALL QUIT('Insufficient space in CCSDT_OMEG')
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
               DO 130 B = 1,NVIR(ISYMB)
C
C-----------------------------------------
C                 Initialize the R2 matrix.
C-----------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
C
C--------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix.
C--------------------------------------------------
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
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT ',XSMAT
                  ENDIF
C
                  IF (IPRINT .GT. 55) THEN
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KTMAT),1)
                     WRITE(LUPRI,*) 'Norm of TMAT ',XTMAT
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
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT ',XQMAT
                  ENDIF
C
C-----------------------------------------
C                 Contract with integrals.
C-----------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_CONVIR(WORK(KRMAT2),WORK(KSMAT),
     *                            WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                            WORK(KTRVI),WORK(KTRVI1),ISINT1,
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 -after CONVIR ',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 -after CONVIR',XRMAT
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_CONVIR: ')
                     CALL CC_PRP(OMEGA1,OMEGA2,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TICONV = TICONV   + DTIME
C
                  DTIME = SECOND()
                  CALL CC3_CONOCC(OMEGA2,WORK(KRMAT1),WORK(KRMAT2),
     *                            WORK(KSMAT),WORK(KTMAT),ISYMIM,
     *                            WORK(KTROC),WORK(KTROC1),ISINT1,
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 -after CONOCC ',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 -after CONOCC',XRMAT
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
                   WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_CONOCC',RHO1N
                   WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_CONOCC',RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_CONOCC: ')
                     CALL CC_PRP(OMEGA1,OMEGA2,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TICONO = TICONO   + DTIME
C
C----------------------------------
C                 Calculate Omega1.
C----------------------------------
C
                  DTIME = SECOND()
C
                  CALL CC3_ONEL(OMEGA1,OMEGA2,WORK(KRMAT1),
     *                          WORK(KRMAT2),WORK(KFCKAK),WORK(KSMAT),
     *                          WORK(KTMAT),ISYMIM,WORK(KXIAJB),ISINT1,
     *                          WORK(KINDSQ),LENSQ,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_ONEL',RHO1N
                     WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_ONEL',RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_ONEL: ')
                     CALL CC_PRP(OMEGA1,OMEGA2,ISYRES,1,1)
                  ENDIF
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 -after ONEL',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 -after ONEL',XRMAT
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT -after ONEL',XSMAT
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KTMAT),1)
                     WRITE(LUPRI,*) 'Norm of TMAT -after ONEL',XTMAT
                  ENDIF
C
C
                  DTIME  = SECOND() - DTIME
                  TIOME1 = TIOME1   + DTIME
C
C----------------------------------------------------
C                 Accumulate the R2 matrix in Omega2.
C----------------------------------------------------
C
                  CALL CC3_RACC(OMEGA2,WORK(KRMAT2),ISYMB,B,ISYRES)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_RACC',RHO1N
                     WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_RACC',RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_RACC: ')
                     CALL CC_PRP(OMEGA1,OMEGA2,ISYRES,1,1)
                  ENDIF
C
  130          CONTINUE
  120       CONTINUE
C
C----------------------------------------------
C           Accumulate the R1 matrix in Omega2.
C----------------------------------------------
C
            CALL CC3_RACC(OMEGA2,WORK(KRMAT1),ISYMD,D,ISYRES)
C
            IF (IPRINT .GT. 55) THEN
               RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
               RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
               WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_RACC-2',RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_RACC-2',RHO2N
            ENDIF
C
            IF (IPRINT .GT. 220) THEN
               CALL AROUND('After CC3_RACC-2: ')
               CALL CC_PRP(OMEGA1,OMEGA2,ISYRES,1,1)
            ENDIF
C
  110    CONTINUE
  100 CONTINUE
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
         WRITE(LUPRI,1) 'CC3_CONV  : ',TICONV
         WRITE(LUPRI,1) 'CC3_CONO  : ',TICONO
         WRITE(LUPRI,1) 'CC3_OME1  : ',TIOME1
         WRITE(LUPRI,*)
      END IF
C
      CALL QEXIT('CC3_OMEG')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck ccsdt_trocc */
      SUBROUTINE CCSDT_TROCC(XINT,TRINT,XLAMDH,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C
C     Kasper Hald, Summer 2002 -> generalisation to nonsymmetric lambda
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*),TRINT(*),XLAMDH(*),WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSDT_TROCC')
C
      DO 100 ISYMD = 1,NSYM
C
         ISYMK  = ISYMD
         ISYAIJ = MULD2H(ISYMD,ISYMOP)
C
         LENGTH = NCKI(ISYAIJ) * NRHF(ISYMK)
         IF (LENGTH .GT. LWORK) THEN
            CALL QUIT('Not enough space in CCSDT_TROCC')
         END IF
C
         NTOIAJ = MAX(NCKI(ISYAIJ),1)
         NBASD  = MAX(NBAS(ISYMD),1)
C
         KOFF1 = ICKID(ISYAIJ,ISYMD)  + 1
         KOFF2 = ILMRHF(ISYMD) + 1
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),NBAS(ISYMD),
     *              ONE,XINT(KOFF1),NTOIAJ,XLAMDH(KOFF2),NBASD,
     *              ZERO,WORK,NTOIAJ)
C
C        Sort
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMA  = MULD2H(ISYMAI,ISYMI)
               ISYMJI = MULD2H(ISYMJ,ISYMI)
               ISYJIK = MULD2H(ISYMJI,ISYMK)
C
               DO 130 K = 1,NRHF(ISYMK)
C
                  DO 140 J = 1,NRHF(ISYMJ)
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        NTOJIK = NMAJIK(ISYJIK)
C
                        KOFF1 = NCKI(ISYAIJ)*(K - 1)
     *                        + ISAIK(ISYMAI,ISYMJ)
     *                        + NT1AM(ISYMAI)*(J - 1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I - 1) + 1
                        KOFF2 = ISJIKA(ISYJIK,ISYMA)
     *                        + ISJIK(ISYMJI,ISYMK)
     *                        + NMATIJ(ISYMJI)*(K - 1)
     *                        + IMATIJ(ISYMJ,ISYMI)
     *                        + NRHF(ISYMJ)*(I - 1) + J
C
                        CALL DCOPY(NVIR(ISYMA),WORK(KOFF1),1,
     *                             TRINT(KOFF2),NTOJIK)
C
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
C
  100 CONTINUE
C
      CALL QEXIT('CCSDT_TROCC')
C
      RETURN
      END
C  /* Deck ccsdt_trvir */
      SUBROUTINE CCSDT_TRVIR(XINT,TRINT,XLAMDP,ISYMD,D,ISYINT,
     *                       WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Transform (ck|alpha d) integrals to (ck|a d)
C
C     Ove Christiansen 11-1-1996: ISYINT is symmetry of (ck|alpha d)
C
C     Kasper Hald, 01062002 : General symmetry of Lambda
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*),TRINT(*),XLAMDP(*),WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSDT_TRVIR')
C
      ISYDIS = MULD2H(ISYMD,ISYINT)
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMCK = MULD2H(ISYMA,ISYDIS)
C
         NTOTCK = MAX(NT1AM(ISYMCK),1)
         NBASA  = MAX(NBAS(ISYMA),1)
C
         KOFF1  = ICKALP(ISYMCK,ISYMA) + 1
         KOFF2  = ILMVIR(ISYMA) + 1
         KOFF3  = ICKATR(ISYMCK,ISYMA) + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMA),NBAS(ISYMA),
     *              ONE,XINT(KOFF1),NTOTCK,XLAMDP(KOFF2),NBASA,
     *              ZERO,TRINT(KOFF3),NTOTCK)
C
  100 CONTINUE
C
      CALL QEXIT('CCSDT_TRVIR')
C
      RETURN
      END
C  /* Deck cc3_smat */
      SUBROUTINE CC3_SMAT(ECURR,T2TP,ISYMT2,TMAT,TRVIR,TRVIR2,TROCC,
     *                    ISYINT,FOCKD,DIAG,SMAT,WORK,LWORK,INDAJL,
     *                    INDSQ,LENSQ,ISYMB,B,ISYMC,C)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Calculate the S matrix for CCSDT.
C
C     S(ai,bj,ck) * D =
C           = t(ai,dj) (kc|bd) + t(di,bj) (kc|ad)
C           - t(ai,bl) (kc|lj) - t(al,bj) (li|kc)
C
C     S is stored as S(ai,k,j) for fixed b and c
C     (kc|bd) is stored as I(dk,b,c)
C
C     Generalised for symmetry: Ove Christiansen Nov. 1995
C                                                Jan. 1996
C
C     General symmetry: ISYINT is symmetr of integrals and ISYMT2 is
C                       symmetry of T2TP.
C
#include <implicit.h>
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION TRVIR(*),TRVIR2(*),TROCC(*),SMAT(*),FOCKD(*)
      DIMENSION INDAJL(*),DIAG(*),T2TP(*),INDSQ(LENSQ,6),TMAT(*)
      DIMENSION WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_SMAT')
C
C--------------------------------
C     First virtual contribution.
C--------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISYINT,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
      ISYMDK = MULD2H(ISYMBC,ISYINT)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient core in CCSDT_SMAT')
      ENDIF
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMD  = MULD2H(ISYMK,ISYMDK)
         ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIJ,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDK,ISYMB) + NT1AM(ISYMDK)*(B - 1)
     *         + IT1AM(ISYMD,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
         NTOAIJ = MAX(1,NCKI(ISYAIJ))
         NVIRD  = MAX(NVIR(ISYMD),1)
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *              NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIJ,
     *              TRVIR(KOFF2),NVIRD,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
  100 CONTINUE
C
      CALL CC_GATHER(LENGTH,SMAT,WORK,INDSQ(1,3))
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 1. Norm of SMAT ',XSMAT
      ENDIF
C
C---------------------------------
C     Second virtual contribution.
C---------------------------------
C
      ISYAKD = MULD2H(ISYMC,ISYINT)
      ISYDIJ = MULD2H(ISYMB,ISYMT2)
C
      DO 200 ISYMJ = 1,NSYM
C
         ISYMDI = MULD2H(ISYMJ,ISYDIJ)
C
         DO 210 J = 1,NRHF(ISYMJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMD  = MULD2H(ISYMDI,ISYMI)
               ISYMAK = MULD2H(ISYMD,ISYAKD)
               ISYAKI = MULD2H(ISYMAK,ISYMI)
C
               KOFF1 = ICKATR(ISYMAK,ISYMD) + 1
C
               KOFF2 = IT2SP(ISYDIJ,ISYMB)
     *               + NCKI(ISYDIJ)*(B - 1)
     *               + ISAIK(ISYMDI,ISYMJ)
     *               + NT1AM(ISYMDI)*(J - 1)
     *               + IT1AM(ISYMD,ISYMI) + 1
C
               KOFF3 = ISAIKJ(ISYAKI,ISYMJ)
     *               + NCKI(ISYAKI)*(J - 1)
     *               + ISAIK(ISYMAK,ISYMI) + 1
C
               NVIRD  = MAX(NVIR(ISYMD),1)
               NTOTAK = MAX(NT1AM(ISYMAK),1)
C
               CALL DGEMM('N','N',NT1AM(ISYMAK),NRHF(ISYMI),
     *                    NVIR(ISYMD),ONE,TRVIR2(KOFF1),NTOTAK,
     *                    T2TP(KOFF2),NVIRD,ZERO,TMAT(KOFF3),
     *                    NTOTAK)
C
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
c     CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,1))
c     CALL DAXPY(LENGTH,ONE,WORK,1,SMAT,1)
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) + TMAT(INDSQ(I,1))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 2. Norm of SMAT ',XSMAT
      ENDIF
C
C---------------------------------
C     First occupied contribution.
C---------------------------------
C
      ISYAIL = MULD2H(ISYMB,ISYMT2)
      ISYLKJ = MULD2H(ISYMC,ISYINT)
C
      DO 300 ISYMJ = 1,NSYM
C
         ISYMLK = MULD2H(ISYMJ,ISYLKJ)
C
         DO 310 J = 1,NRHF(ISYMJ)
C
            DO 320 ISYMK = 1,NSYM
C
               ISYML  = MULD2H(ISYMK,ISYMLK)
               ISYMAI = MULD2H(ISYAIL,ISYML)
               ISYAIK = MULD2H(ISYMAI,ISYMK)
C
               KOFF1 = IT2SP(ISYAIL,ISYMB)
     *               + NCKI(ISYAIL)*(B - 1)
     *               + ICKI(ISYMAI,ISYML) + 1
               KOFF2 = ISJIKA(ISYLKJ,ISYMC)
     *               + NMAJIK(ISYLKJ)*(C - 1)
     *               + ISJIK(ISYMLK,ISYMJ)
     *               + NMATIJ(ISYMLK)*(J - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
               KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *               + NCKI(ISYAIK)*(J - 1)
     *               + ICKI(ISYMAI,ISYMK) + 1
C
               NTOTAI = MAX(1,NT1AM(ISYMAI))
               NRHFL  = MAX(1,NRHF(ISYML))
C
               CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                    NRHF(ISYML),-ONE,T2TP(KOFF1),NTOTAI,
     *                    TROCC(KOFF2),NRHFL,ONE,SMAT(KOFF3),
     *                    NTOTAI)
C
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 3. Norm of SMAT ',XSMAT
      ENDIF
C
C----------------------------------
C     Second occupied contribution.
C----------------------------------
C
      ISYAJL = MULD2H(ISYMB,ISYMT2)
      ISYLKI = MULD2H(ISYMC,ISYINT)
C
      IF (LWORK .LT. NCKI(ISYAJL)) THEN
         CALL QUIT('Not enough space in CCSDT_SMAT')
      END IF
C
      KOFF = IT2SP(ISYAJL,ISYMB) + NCKI(ISYAJL)*(B - 1) + 1
      CALL CC_GATHER(NCKI(ISYAJL),WORK,T2TP(KOFF),INDAJL)
C
      DO 400 ISYMI = 1,NSYM
C
         ISYMLK = MULD2H(ISYMI,ISYLKI)
C
         DO 410 I = 1,NRHF(ISYMI)
C
            DO 420 ISYMK = 1,NSYM
C
               ISYML  = MULD2H(ISYMK,ISYMLK)
               ISYMAJ = MULD2H(ISYAJL,ISYML)
               ISYAJK = MULD2H(ISYMAJ,ISYMK)
C
               KOFF1 = ICKI(ISYMAJ,ISYML) + 1
C
               KOFF2 = ISJIKA(ISYLKI,ISYMC)
     *               + NMAJIK(ISYLKI)*(C - 1)
     *               + ISJIK(ISYMLK,ISYMI)
     *               + NMATIJ(ISYMLK)*(I - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
C
               KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *               + NCKI(ISYAJK)*(I - 1)
     *               + ICKI(ISYMAJ,ISYMK) + 1
C
               NTOTAJ = MAX(1,NT1AM(ISYMAJ))
               NRHFL  = MAX(1,NRHF(ISYML))
C
               CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMK),
     *                    NRHF(ISYML),ONE,WORK(KOFF1),NTOTAJ,
     *                    TROCC(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAJ)
C
  420       CONTINUE
  410    CONTINUE
  400 CONTINUE
C
c     CALL CC_GATHER(NCKIJ(JSAIKJ),WORK,TMAT,INDSQ(1,5))
c     CALL DAXPY(NCKIJ(JSAIKJ),-ONE,WORK,1,SMAT,1)
C
      DO I = 1,NCKIJ(JSAIKJ)
         SMAT(I) = SMAT(I) - TMAT(INDSQ(I,5))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 4. Norm of SMAT ',XSMAT
      ENDIF
C
C-----------------------------------------
C     Divide by the Fock matrix diagonals.
C-----------------------------------------
C
      NB = IORB(ISYMB) + NRHF(ISYMB) + B
      NC = IORB(ISYMC) + NRHF(ISYMC) + C
C
      EPSIBC = FOCKD(NB) + FOCKD(NC) - ECURR
C
      DO 500 L = 1,NCKIJ(JSAIKJ)
C
         SMAT(L) = SMAT(L)/(DIAG(L) + EPSIBC)
C
  500 CONTINUE
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 5. Norm of SMAT ',XSMAT
      ENDIF
C
C----------------------
C     Print if desired.
C----------------------
C
      IF (IPRCC .GT. 75) THEN
C
         CALL AROUND('The S matrix')
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(2X,A,I4)')  'JSAIKJ ', JSAIKJ
         WRITE(LUPRI,'(2X,A,4I4)') 'isymb,b,isymc,c',ISYMB,B,ISYMC,C
         WRITE(LUPRI,*)
C
         DO 600 ISYMJ = 1,NSYM
C
            ISYAIK = MULD2H(JSAIKJ,ISYMJ)
C
            DO 610 J = 1,NRHF(ISYMJ)
C
               WRITE(LUPRI,'(5X,A,2I4)') 'isymj,j',ISYMJ,J
               WRITE(LUPRI,*)
C
               DO 620 ISYMK = 1,NSYM
C
                  ISYMAI = MULD2H(ISYAIK,ISYMK)
C
                  DO 630 K = 1,NRHF(ISYMK)
C
                     WRITE(LUPRI,'(8X,A,2I4)') 'isymk,k',ISYMK,K
                     WRITE(LUPRI,*)
C
                     DO 640 ISYMI = 1,NSYM
C
                        ISYMA = MULD2H(ISYMAI,ISYMI)
C
                        KOFF = ISAIKJ(ISYAIK,ISYMJ)
     *                       + NCKI(ISYAIK)*(J - 1)
     *                       + ICKI(ISYMAI,ISYMK)
     *                       + NT1AM(ISYMAI)*(K - 1)
     *                       + IT1AM(ISYMA,ISYMI) + 1
C
                        CALL OUTPUT(SMAT(KOFF),1,NVIR(ISYMA),1,
     *                              NRHF(ISYMI),NVIR(ISYMA),
     *                              NRHF(ISYMI),1,LUPRI)
C
  640                CONTINUE
  630             CONTINUE
  620          CONTINUE
  610       CONTINUE
  600    CONTINUE
C
      END IF
C
      CALL QEXIT('CC3_SMAT')
C
      RETURN
      END
C  /* Deck ccsdt_srtvir */
      SUBROUTINE CCSDT_SRTVIR(TRVIR,TRVIR2,WORK,LWORK,ISYMC,ISYINT)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Modified Feb 1995.
C
C     Resort virtual integrals
C
C     (kc|ad) for fixed c
C     initially stored I(dk,a,c) to finally I'(a,k,d,c)
C
C     Ove Christiansen 11-1-1996: ISYINT general symmetry of integral I(dk,a,c)
C
#include <implicit.h>
C
      DIMENSION TRVIR(*),TRVIR2(*)
      DIMENSION WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSDT_SRTVIR')
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMAC = MULD2H(ISYMA,ISYMC)
         ISYMDK = MULD2H(ISYMAC,ISYINT)
C
         DO 110 ISYMK = 1,NSYM
C
            ISYMD  = MULD2H(ISYMDK,ISYMK)
            ISYMAK = MULD2H(ISYMA,ISYMK)
C
            DO 120 A = 1,NVIR(ISYMA)
C
               DO 130 K = 1,NRHF(ISYMK)
C
                  DO 140 D = 1,NVIR(ISYMD)
C
                     NDKA  = ICKATR(ISYMDK,ISYMA)
     *                     + NT1AM(ISYMDK)*(A - 1)
     *                     + IT1AM(ISYMD,ISYMK)
     *                     + NVIR(ISYMD)*(K - 1) + D
C
                     NAKD  = ICKATR(ISYMAK,ISYMD)
     *                     + NT1AM(ISYMAK)*(D - 1)
     *                     + IT1AM(ISYMA,ISYMK)
     *                     + NVIR(ISYMA)*(K - 1) + A
C
                     TRVIR2(NAKD) = TRVIR(NDKA)
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCSDT_SRTVIR')
C
      RETURN
      END
C  /* Deck ccsdt_srtocc */
      SUBROUTINE CCSDT_SRTOCC(TROCC,TROCC2,ISYINT)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Resort occupied integrals
C
C     (ia|jk)
C     initially stored I(ai,j,k) to finally I'(j,i,k,a)
C
C     Ove Christiansen 11-1-1996: ISYINT symmetry of I(ai,j,k)
C
#include <implicit.h>
C
      DIMENSION TROCC(*),TROCC2(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSDT_SRTOCC')
C
      DO 100 ISYMK = 1,NSYM
C
         ISYAIJ = MULD2H(ISYMK,ISYINT)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMA  = MULD2H(ISYMAI,ISYMI)
               ISYMJI = MULD2H(ISYMJ,ISYMI)
               ISYJIK = MULD2H(ISYMJI,ISYMK)
C
               DO 130 K = 1,NRHF(ISYMK)
C
                  DO 140 J = 1,NRHF(ISYMJ)
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        NTOJIK = NMAJIK(ISYJIK)
C
                        KOFF1 = ISAIKJ(ISYAIJ,ISYMK)
     *                        + NCKI(ISYAIJ)*(K - 1)
     *                        + ISAIK(ISYMAI,ISYMJ)
     *                        + NT1AM(ISYMAI)*(J - 1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I - 1) + 1
                        KOFF2 = ISJIKA(ISYJIK,ISYMA)
     *                        + ISJIK(ISYMJI,ISYMK)
     *                        + NMATIJ(ISYMJI)*(K - 1)
     *                        + IMATIJ(ISYMJ,ISYMI)
     *                        + NRHF(ISYMJ)*(I - 1) + J
C
                        CALL DCOPY(NVIR(ISYMA),TROCC(KOFF1),1,
     *                             TROCC2(KOFF2),NTOJIK)
C
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCSDT_SRTOCC')
C
      RETURN
      END
C  /* Deck cc3_index */
      SUBROUTINE CC3_INDEX(INDAJL,ISYALJ)
C
C     Henrik Koch.          Feb 1994
C
#include <implicit.h>
C
      DIMENSION INDAJL(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_INDEX')
C
C----------------------------------
C     Calculation of index vectors.
C----------------------------------
C
      DO 100 ISYML = 1,NSYM
C
         ISYMAJ = MULD2H(ISYML,ISYALJ)
C
         DO 110 L = 1,NRHF(ISYML)
C
            DO 120 ISYMJ = 1,NSYM
C
               ISYMA  = MULD2H(ISYMJ,ISYMAJ)
               ISYMAL = MULD2H(ISYMA,ISYML)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  DO 140 A = 1,NVIR(ISYMA)
C
                     NAJL = ICKI(ISYMAJ,ISYML)
     *                    + NT1AM(ISYMAJ)*(L - 1)
     *                    + IT1AM(ISYMA,ISYMJ)
     *                    + NVIR(ISYMA)*(J - 1) + A
C
                     INDAJL(NAJL) = ICKI(ISYMAL,ISYMJ)
     *                            + NT1AM(ISYMAL)*(J - 1)
     *                            + IT1AM(ISYMA,ISYML)
     *                            + NVIR(ISYMA)*(L - 1) + A
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_INDEX')
C
      RETURN
      END
C  /* Deck cc3_qmat */
      SUBROUTINE CC3_QMAT(ECURR,T2TP,ISYMT2,TRVIR,TROCC,ISYINT,FOCKD,
     *                    DIAG,
     *                    QMAT,WORK,LWORK,INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Calculate Q(ci,k,j) = t(ci,d'j) (kb|dd') - t(ci,dl) (kb|lj)
C
C               (kb|dd') = I(d'k,b,d) = I(fk,b,d)
C
C     Ove Christiansen 9-1-1996.
C
C     General symmetry: ISYINT is symmetr of integrals and ISYMT2 is
C                       symmetry of T2TP.
C
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION T2TP(*),TRVIR(*),TROCC(*),FOCKD(*),QMAT(*)
      DIMENSION WORK(LWORK),DIAG(*),INDSQ(LENSQ,6)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_QMAT')
C
      ISYRES = MULD2H(ISYMT2,ISYINT)
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYMFK = MULD2H(ISYMBD,ISYINT)
      ISCIKJ = MULD2H(ISYMBD,ISYRES)
C
C--------------------------
C     Virtual contribution.
C--------------------------
C
      LENGTH = NCKIJ(ISCIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient core in CC3_QMAT')
      ENDIF
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMF  = MULD2H(ISYMK,ISYMFK)
         ISYCIJ = MULD2H(ISYMK,ISCIKJ)
C
         KOFF1 = IT2SP(ISYCIJ,ISYMF)  + 1
         KOFF2 = ICKATR(ISYMFK,ISYMB) + NT1AM(ISYMFK)*(B - 1)
     *         + IT1AM(ISYMF,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYCIJ,ISYMK) + 1
C
         NTOCIJ = MAX(1,NCKI(ISYCIJ))
         NVIRF  = MAX(NVIR(ISYMF),1)
C
         CALL DGEMM('N','N',NCKI(ISYCIJ),NRHF(ISYMK),
     *              NVIR(ISYMF),ONE,T2TP(KOFF1),NTOCIJ,
     *              TRVIR(KOFF2),NVIRF,ZERO,
     *              WORK(KOFF3),NTOCIJ)
C
  100 CONTINUE
C
      CALL CC_GATHER(LENGTH,QMAT,WORK,INDSQ(1,3))
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(ISCIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CC3_QMAT: 1. Norm of QMAT ',XQMAT
      ENDIF
C
C---------------------------
C     Occupied contribution.
C---------------------------
C
      ISYCIL = MULD2H(ISYMD,ISYMT2)
      ISYLKJ = MULD2H(ISYMB,ISYINT)
C
      DO 220 ISYMJ = 1,NSYM
C
         ISYMLK = MULD2H(ISYLKJ,ISYMJ)
         ISYCIK = MULD2H(ISCIKJ,ISYMJ)
C
         DO 230 J = 1,NRHF(ISYMJ)
C
            DO 240 ISYMK = 1,NSYM
C
               ISYMCI = MULD2H(ISYCIK,ISYMK)
               ISYML  = MULD2H(ISYMLK,ISYMK)
C
               NTOTCI = MAX(NT1AM(ISYMCI),1)
               NRHFL  = MAX(NRHF(ISYML),1)
C
               KOFF1 = IT2SP(ISYCIL,ISYMD)
     *               + NCKI(ISYCIL)*(D - 1)
     *               + ISAIK(ISYMCI,ISYML)  + 1
               KOFF2 = ISJIKA(ISYLKJ,ISYMB) + NMAJIK(ISYLKJ)*(B - 1)
     *               + ISJIK(ISYMLK,ISYMJ)  + NMATIJ(ISYMLK)*(J - 1)
     *               + IMATIJ(ISYML,ISYMK)  + 1
               KOFF3 = ISAIKJ(ISYCIK,ISYMJ) + NCKI(ISYCIK)*(J - 1)
     *               + ISAIK(ISYMCI,ISYMK)  + 1
C
               CALL DGEMM('N','N',NT1AM(ISYMCI),NRHF(ISYMK),NRHF(ISYML),
     *                    -ONE,T2TP(KOFF1),NTOTCI,TROCC(KOFF2),NRHFL,
     *                    ONE,QMAT(KOFF3),NTOTCI)
C
  240       CONTINUE
  230    CONTINUE
  220 CONTINUE
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(ISCIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CC3_QMAT: 2. Norm of QMAT ',XQMAT
      ENDIF
C
C-----------------------------------------
C     Divide by the Fock matrix diagonals.
C-----------------------------------------
C
      NB = IORB(ISYMB) + NRHF(ISYMB) + B
      ND = IORB(ISYMD) + NRHF(ISYMD) + D
C
      EPSIBD = FOCKD(NB) + FOCKD(ND) - ECURR
C
      DO 300 L = 1,NCKIJ(ISCIKJ)
C
         QMAT(L) = QMAT(L)/(DIAG(L) + EPSIBD)
C
  300 CONTINUE
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(ISCIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CC3_QMAT: 3. Norm of QMAT ',XQMAT
      ENDIF
C
C----------------------
C     Print if desired.
C----------------------
C
      IF (IPRCC .GT. 75) THEN
C
         CALL AROUND('The Q matrix')
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(2X,A,4I4)') 'isymb,b,isymd,d',ISYMB,B,ISYMD,D
         WRITE(LUPRI,*)
C
         DO 400 ISYMJ = 1,NSYM
C
            ISYCIK = MULD2H(ISCIKJ,ISYMJ)
C
            DO 410 J = 1,NRHF(ISYMJ)
C
               WRITE(LUPRI,'(5X,A,2I4)') 'isymj,j',ISYMJ,J
               WRITE(LUPRI,*)
C
               DO 420 ISYMK = 1,NSYM
C
                  ISYMCI = MULD2H(ISYCIK,ISYMK)
C
                  DO 430 K = 1,NRHF(ISYMK)
C
                     WRITE(LUPRI,'(8X,A,2I4)') 'isymk,k',ISYMK,K
                     WRITE(LUPRI,*)
C
                     DO 440 ISYMI = 1,NSYM
C
                        ISYMC = MULD2H(ISYMCI,ISYMI)
C
                        KOFF = ISAIKJ(ISYCIK,ISYMJ)
     *                       + NCKI(ISYCIK)*(J - 1)
     *                       + ICKI(ISYMCI,ISYMK)
     *                       + NT1AM(ISYMCI)*(K - 1)
     *                       + IT1AM(ISYMC,ISYMI) + 1
C
                        CALL OUTPUT(QMAT(KOFF),1,NVIR(ISYMC),1,
     *                              NRHF(ISYMI),NVIR(ISYMC),
     *                              NRHF(ISYMI),1,LUPRI)
C
  440                CONTINUE
  430             CONTINUE
  420          CONTINUE
  410       CONTINUE
  400    CONTINUE
C
      END IF
C
      CALL QEXIT('CC3_QMAT')
C
      RETURN
      END
C  /* Deck ccsdt_srvir3 */
      SUBROUTINE CCSDT_SRVIR3(TRVIR,SCR,ISYMD,D,ISYINT)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Sort (kb|dd') = I(bk,d',d) as I'(d'k,b,d) =  I'(fk,b,d)
C
C     Ove Christiansen 11-1-1996: ISYINT = symmetry of I(bk,d',d)
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION TRVIR(*),SCR(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSDT_SRVIR3')
C
      ISYBKF = MULD2H(ISYMD,ISYINT)
C
      DO 100 ISYMF = 1,NSYM
C
         ISYMBK = MULD2H(ISYBKF,ISYMF)
C
         DO 110 F = 1,NVIR(ISYMF)
C
            DO 120 ISYMK = 1,NSYM
C
               ISYMB  = MULD2H(ISYMBK,ISYMK)
               ISYMFK = MULD2H(ISYMF,ISYMK)
C
               DO 130 K = 1,NRHF(ISYMK)
C
                  KOFF1 = ICKATR(ISYMBK,ISYMF) + NT1AM(ISYMBK)*(F - 1)
     *                  + IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1) + 1
                  KOFF2 = ICKATR(ISYMFK,ISYMB)
     *                  + IT1AM(ISYMF,ISYMK) + NVIR(ISYMF)*(K - 1) + F
C
                  CALL DCOPY(NVIR(ISYMB),TRVIR(KOFF1),1,
     *                       SCR(KOFF2),NT1AM(ISYMFK))
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL DCOPY(NCKATR(ISYBKF),SCR,1,TRVIR,1)
C
      CALL QEXIT('CCSDT_SRVIR3')
C
      RETURN
      END
C  /* Deck cc3_convir */
      SUBROUTINE CC3_CONVIR(RMAT,SMAT,QMAT,TMAT,ISYMIM,TRVIR,
     *                      TRVIR1,ISYINT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMB,B,ISYMD,D)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Set up combinations of Q's and S's and contract with integrals.
C
C     Ove Christiansen 9-1-1996:
C
C     General symmetry: ISYMIM is symmetry of SMAT and TMAT intermdiates.
C                       ISYINT is symmetry of FOCKAK and XIAJB
C                       ISYRES = ISYMIM*ISYINT
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RMAT(*),SMAT(*),QMAT(*)
      DIMENSION TMAT(*),TRVIR(*),TRVIR1(*)
      DIMENSION WORK(LWORK),INDSQ(LENSQ,6)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CONVIR')
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
      LENGTH = NCKIJ(ISCKIJ)
C
C------------------------
C     First virtual term.
C------------------------
C
      IF (LWORK .LT. NCKIJ(ISCKIJ)) THEN
         CALL QUIT('Insufficient core in CC3_CONVIR')
      ENDIF
C
c     CALL DCOPY(LENGTH,SMAT,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,1))
c     CALL DAXPY(LENGTH,-TWO,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,QMAT,INDSQ(1,1))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,2))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,QMAT,INDSQ(1,2))
c     CALL DAXPY(LENGTH,-TWO,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,QMAT,INDSQ(1,3))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
      DO I = 1,LENGTH
C
         TMAT(I) =     SMAT(I) 
     *           - TWO*SMAT(INDSQ(I,1))
     *           +     QMAT(INDSQ(I,1))
     *           +     SMAT(INDSQ(I,2))
     *           - TWO*QMAT(INDSQ(I,2))
     *           +     QMAT(INDSQ(I,3))
C
      ENDDO
C
C---------------------------
C     Contract with (ac|kd).
C---------------------------
C
      DO 200 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYCKI = MULD2H(ISCKIJ,ISYMJ)
C
         KSCR1  = 1
         KEND1  = KSCR1 + NT1AM(ISYMAI)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC3_CONVIR')
         ENDIF
C
         DO 210 J = 1,NRHF(ISYMJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYMAI,ISYMI)
C
               NTOTCK = MAX(NT1AM(ISYMCK),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               KOFF1 = ICKATR(ISYMCK,ISYMA) + 1
               KOFF2 = ISAIKJ(ISYCKI,ISYMJ)
     *               + NCKI(ISYCKI)*(J - 1)
     *               + ISAIK(ISYMCK,ISYMI)  + 1
               KOFF3 = ISAIK(ISYMAI,ISYMJ)
     *               + NT1AM(ISYMAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NT1AM(ISYMCK),
     *                    ONE,TRVIR(KOFF1),NTOTCK,TMAT(KOFF2),NTOTCK,
     *                    ONE,RMAT(KOFF3),NVIRA)
C
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
C-------------------------
C     Second virtual term.
C-------------------------
C
c     CALL DCOPY(LENGTH,SMAT,1,TMAT,1)
c     CALL DSCAL(LENGTH,-TWO,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,1))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,QMAT,INDSQ(1,2))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,QMAT,INDSQ(1,3))
c     CALL DAXPY(LENGTH,-TWO,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,QMAT,INDSQ(1,4))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,5))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
      DO I = 1,LENGTH
C
         TMAT(I) = - TWO*SMAT(I) 
     *             +     SMAT(INDSQ(I,1))
     *             +     QMAT(INDSQ(I,2))
     *             - TWO*QMAT(INDSQ(I,3))
     *             +     QMAT(INDSQ(I,4))
     *             +     SMAT(INDSQ(I,5))
C
      ENDDO
C
C---------------------------
C     Contract with (ad|kc).
C---------------------------
C
      DO 300 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYCKI = MULD2H(ISCKIJ,ISYMJ)
C
         KSCR1  = 1
         KEND1  = KSCR1 + NT1AM(ISYMAI)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC3_CONVIR')
         ENDIF
C
         DO 310 J = 1,NRHF(ISYMJ)
C
            DO 320 ISYMI = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYMAI,ISYMI)
C
               NTOTCK = MAX(NT1AM(ISYMCK),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               KOFF1 = ICKATR(ISYMCK,ISYMA) + 1
               KOFF2 = ISAIKJ(ISYCKI,ISYMJ)
     *               + NCKI(ISYCKI)*(J - 1)
     *               + ISAIK(ISYMCK,ISYMI)  + 1
               KOFF3 = ISAIK(ISYMAI,ISYMJ)
     *               + NT1AM(ISYMAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NT1AM(ISYMCK),
     *                    ONE,TRVIR1(KOFF1),NTOTCK,TMAT(KOFF2),NTOTCK,
     *                    ONE,RMAT(KOFF3),NVIRA)
C
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
C
      CALL QEXIT('CC3_CONVIR')
C
      RETURN
      END
C  /* Deck cc3_conocc */
      SUBROUTINE CC3_CONOCC(OMEGA2,RMAT1,RMAT2,SMAT,TMAT,ISYMIM,TROCC,
     *                      TROCC1,ISYINT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMIB,IB,ISYMID,ID)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Set up combinations of Q's and S's and contract with integrals.
C
C     Ove Christiansen 9-1-1996:
C
C     General symmetry: ISYMIM is symmetry of SMAT and TMAT intermediates.
C                       (including isymib*isymid)
C                       ISYINT is symmetry of integrals in TROCC and TROCC1.
C                       ISYRES = ISYMIM*ISYINT
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION OMEGA2(*),RMAT1(*),RMAT2(*),SMAT(*),TMAT(*)
      DIMENSION TROCC(*),TROCC1(*),WORK(LWORK),INDSQ(LENSQ,6)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CONOCC')
C
      IF (LWORK .LT. LENSQ) THEN
         CALL QUIT('Insufficient core in CONOCC')
      ENDIF
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
C-------------------------
C     First occupied term.
C-------------------------
C
      C = ID
      B = IB
C
      ISYMC = ISYMID
      ISYMB = ISYMIB
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKL = MULD2H(ISYMBC,ISYMIM)
C
      LENGTH = NCKIJ(JSAIKL)
C
C----------------------------------
C     Setup combinations of smat's.
C----------------------------------
C
c     CALL DCOPY(LENGTH,SMAT,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,3))
c     CALL DAXPY(LENGTH,-TWO,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,4))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
      DO I = 1,LENGTH
C
         TMAT(I) =       SMAT(I) 
     *             - TWO*SMAT(INDSQ(I,3))
     *             +     SMAT(INDSQ(I,4))
C
      ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: 1. Norm of TMAT = ',XTMAT
      ENDIF
C
C-----------------------
C     First contraction.
C-----------------------
C
      DO 200 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYMKL = MULD2H(JSAIKL,ISYMAI)
         ISYKLJ = MULD2H(ISYMKL,ISYMJ)
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         KOFF1  = ISAIKL(ISYMAI,ISYMKL) + 1
         KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *          + NMAJIK(ISYKLJ)*(C - 1)
     *          + ISJIK(ISYMKL,ISYMJ) + 1
         KOFF3  = ISAIK(ISYMAI,ISYMJ) + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *              -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *              ONE,RMAT2(KOFF3),NTOTAI)
C
  200 CONTINUE
C
      IF (IPRINT .GT. 55) THEN
         XRMAT = DDOT(NCKI(ISYRES),RMAT2,1,RMAT2,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: Norm of RMAT2 =  ',XRMAT
      ENDIF
C
C--------------------------
C     Second occupied term.
C--------------------------
C
      B = ID
      C = IB
C
      ISYMB = ISYMID
      ISYMC = ISYMIB
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKL = MULD2H(ISYMBC,ISYMIM)
C
      LENGTH = NCKIJ(JSAIKL)
C
C----------------------------------
C     Setup combinations of smat's.
C----------------------------------
C
c     CALL DCOPY(LENGTH,SMAT,1,TMAT,1)
c     CALL DSCAL(LENGTH,-TWO,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,3))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,5))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
      DO I = 1,LENGTH
C
         TMAT(I) = - TWO*SMAT(I) 
     *             +     SMAT(INDSQ(I,3))
     *             +     SMAT(INDSQ(I,5))
C
      ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: 2. Norm of TMAT = ',XTMAT
      ENDIF
C
C------------------------
C     Second contraction.
C------------------------
C
      DO 400 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYMKL = MULD2H(JSAIKL,ISYMAI)
         ISYKLJ = MULD2H(ISYMKL,ISYMJ)
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         KOFF1  = ISAIKL(ISYMAI,ISYMKL) + 1
         KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *          + NMAJIK(ISYKLJ)*(C - 1)
     *          + ISJIK(ISYMKL,ISYMJ) + 1
         KOFF3  = ISAIK(ISYMAI,ISYMJ) + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *              -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *              ONE,RMAT1(KOFF3),NTOTAI)
C
  400 CONTINUE
C
      IF (IPRINT .GT. 55) THEN
         XRMAT = DDOT(NCKI(ISYRES),RMAT1,1,RMAT1,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: Norm of RMAT1 =  ',XRMAT
      ENDIF
C
C-------------------------
C     Third occupied term.
C-------------------------
C
      A = ID
      B = IB
C
      ISYMA = ISYMID
      ISYMB = ISYMIB
C
      ISYMAB = MULD2H(ISYMA,ISYMB)
      JSCKLI = MULD2H(ISYMAB,ISYMIM)
C
      LENGTH = NCKIJ(JSCKLI)
C
C----------------------------------
C     Setup combinations of smat's.
C----------------------------------
C
c     CALL CC_GATHER(LENGTH,TMAT,SMAT,INDSQ(1,5))
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,2))
c     CALL DAXPY(LENGTH,-TWO,WORK,1,TMAT,1)
C
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,3))
c     CALL DAXPY(LENGTH,ONE,WORK,1,TMAT,1)
C
      DO I = 1,LENGTH
C
         TMAT(I) =       SMAT(INDSQ(I,5)) 
     *             - TWO*SMAT(INDSQ(I,2))
     *             +     SMAT(INDSQ(I,3))
C
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: 3. Norm of TMAT = ',XTMAT
      ENDIF
C
C-----------------------
C     Third contraction.
C-----------------------
C
      DO 600 ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYMI  = MULD2H(ISYMAI,ISYMA)
         ISYCKL = MULD2H(ISYMI,JSCKLI)
C
         IF (LWORK .LT. NRHF(ISYMI)*NRHF(ISYMJ)) THEN
            CALL QUIT('Insufficient memory in CC3_CONOCC')
         END IF
C
         NTOCKL = MAX(NCKI(ISYCKL),1)
         NRHFI  = MAX(NRHF(ISYMI),1)
C
         KOFF1  = ISAIKJ(ISYCKL,ISYMI) + 1
         KOFF2  = ISAIKJ(ISYCKL,ISYMJ) + 1
         KOFF3  = 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NCKI(ISYCKL),
     *              ONE,TMAT(KOFF1),NTOCKL,TROCC1(KOFF2),NTOCKL,
     *              ZERO,WORK(KOFF3),NRHFI)
C
         DO 610 J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            IF (ISYMAI.EQ.ISYMBJ) THEN
C
               DO 620 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                 + INDEX(NAI,NBJ)
C
                  KOFF6 = NRHF(ISYMI)*(J - 1) + I
C
                  IF (NAI .EQ. NBJ) WORK(KOFF6) = TWO*WORK(KOFF6)
C
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) - WORK(KOFF6)
C
  620          CONTINUE
C
            ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO 630 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ-1) + NAI
C
                  KOFF6 = NRHF(ISYMI)*(J - 1) + I
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) - WORK(KOFF6)
C
  630          CONTINUE
C
            ELSE IF (ISYMBJ .LT. ISYMAI) THEN
C
               DO 640 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMBJ)*(NAI-1) + NBJ
C
                  KOFF6 = NRHF(ISYMI)*(J - 1) + I
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) - WORK(KOFF6)
C
  640          CONTINUE
C
            ENDIF
C
  610    CONTINUE
C
  600 CONTINUE
C
      CALL QEXIT('CC3_CONOCC')
C
      RETURN
      END
C  /* Deck ccsdt_srtoc2 */
      SUBROUTINE CCSDT_SRTOC2(TROCC,TROCC1,ISYINT,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Sort occupied integrals I(kl,j,c) as I'(ck,l,j)
C
C     Ove Christiansen 16-1-1996: ISYINT sym of I(kl,j,c)
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION TROCC(*),TROCC1(*)
      DIMENSION WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSDT_SRTOC2')
C
      KOFF2 = 0
      DO 100 ISYMC = 1,NSYM
C
         ISYKLJ = MULD2H(ISYMC,ISYINT)
C
         DO 110 C = 1,NVIR(ISYMC)
C
            DO 120 ISYMJ = 1,NSYM
C
               ISYMKL = MULD2H(ISYKLJ,ISYMJ)
               ISYCKL = MULD2H(ISYMJ,ISYINT)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  DO 140 ISYML = 1,NSYM
C
                     ISYMK  = MULD2H(ISYMKL,ISYML)
                     ISYMCK = MULD2H(ISYCKL,ISYML)
C
                     DO 150 L = 1,NRHF(ISYML)
C
                        DO 160 K = 1,NRHF(ISYMK)
C
                           KOFF1 = ISAIKJ(ISYCKL,ISYMJ)
     *                           + NCKI(ISYCKL)*(J - 1)
     *                           + ISAIK(ISYMCK,ISYML)
     *                           + NT1AM(ISYMCK)*(L - 1)
     *                           + IT1AM(ISYMC,ISYMK)
     *                           + NVIR(ISYMC)*(K - 1) + C
C
                           KOFF2 = KOFF2 + 1
C
                           TROCC1(KOFF1) = TROCC(KOFF2)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCSDT_SRTOC2')
C
      RETURN
      END
C  /* Deck ccsd_tcmepk */
      SUBROUTINE CCSD_TCMEPK(T2AM,SCAL,ISYOPE,IOPT)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Kasper Hald Jan. 2001.
C     Added IOPT and IOPT = 2.
C
C     IOPT = 1
C     Calculate in place two coulomb minus exchange of t2 amplitudes.
C     IOPT = 2
C     Calculate in place minus exchange of t2 amplitudes.
C
C     N.B.: quick hack to fix symmetry bug, 
C           but now crappy and propably inefficient...
C           Ch. Haettig
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION T2AM(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_TCMEPK')
C
C----------------------------
C     Sanity check for IOPT.
C----------------------------
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2)) THEN
         CALL QUIT('Wrong IOPT in CCSD_TCMEPK')
      ENDIF
C
C     IF (ISYOPE.NE.1) THEN
C       WRITE (LUPRI,*) 'He were are in CCSD_TCMEPK...'
C     END IF
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYOPE)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
            IF (ISYMI .GT. ISYMJ) GOTO 110
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMA = MULD2H(ISYMB,ISYMAB)
C
               IF (ISYMA .GT. ISYMB) GOTO 120
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
C              IF (ISYOPE.NE.1) THEN
C                WRITE (LUPRI,*) '... and we entered also the loop...'
C              END IF
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) THEN
                     NRHFI =  J
                  ELSE
                     NRHFI = NRHF(ISYMI)
                  ENDIF
C
                  DO 140 I = 1,NRHFI
C
                     DO 150 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 160 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                         IF (ISYMAI.EQ.ISYMBJ) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + INDEX(NAI,NBJ)
                         ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI)*(NBJ-1)+NAI
                         ELSE IF (ISYMBJ.LT.ISYMAI) THEN
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMBJ)*(NAI-1)+NBJ
                         END IF
C
                         IF (ISYMAJ.EQ.ISYMBI) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + INDEX(NAJ,NBI)
                         ELSE IF (ISYMAJ.LT.ISYMBI) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ)*(NBI-1)+NAJ
                         ELSE IF (ISYMBI.LT.ISYMAJ) THEN
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMBI)*(NAJ-1)+NBI
                         END IF
C
                         IF (IOPT .EQ. 1) THEN
                            XAIBJ = TWO*T2AM(NAIBJ) - T2AM(NAJBI)
                            XAJBI = TWO*T2AM(NAJBI) - T2AM(NAIBJ)
                         ELSE IF (IOPT .EQ. 2) THEN
                            XAIBJ = - T2AM(NAJBI)
                            XAJBI = - T2AM(NAIBJ)
                         ENDIF
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------------------------
C     Scale diagonal elements of result.
C---------------------------------------
C
      IF (ISYOPE .NE. 1) THEN
         CALL QEXIT('CCSD_TCMEPK')
         RETURN
      ENDIF
C
      DO 200 ISYMAI = 1,NSYM
         DO 210 NAI = 1,NT1AM(ISYMAI)
            NAIAI = IT2AM(ISYMAI,ISYMAI) + INDEX(NAI,NAI)
            T2AM(NAIAI) = SCAL*T2AM(NAIAI)
  210    CONTINUE
  200 CONTINUE
C
      CALL QEXIT('CCSD_TCMEPK')
C
      RETURN
      END
C  /* Deck cc3_onel */
      SUBROUTINE CC3_ONEL(OMEGA1,OMEGA2,RMAT1,RMAT2,FOCKAK,SMAT,TMAT,
     *                    ISYMIM,XIAJB,ISYINT,INDSQ,LENSQ,WORK,LWORK,
     *                    ISYMIB,IB,ISYMID,ID)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Calculate Omega1 and Fock contibution to Omega2.
C
C     Ove Christiansen 9-1-1996:
C
C     General symmetry: ISYMIM is symmetry of SMAT and TMAT intermdiates.(incl isymd,isymb)
C                       ISYINT is symmetry of FOCKAK and XIAJB
C                       ISYRES = ISYMIM*ISYINT
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION OMEGA1(*),OMEGA2(*),RMAT1(*),RMAT2(*),FOCKAK(*),SMAT(*)
      DIMENSION TMAT(*),XIAJB(*),WORK(LWORK),INDSQ(LENSQ,6)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_ONEL')
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      B = IB
      C = ID
C
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
C----------------------------------
C     First contribution to Omega2.
C----------------------------------
C
      ISYMK  = MULD2H(ISYMC,ISYINT)
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYMBC,ISYMIM)
      ISYAIJ = MULD2H(ISYMK,JSAIKJ)
      ISYMCK = MULD2H(ISYMC,ISYMK)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CC3_ONEL')
      END IF
C
c     CALL CC_GATHER(LENGTH,TMAT,SMAT,INDSQ(1,4))
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,3))
c     CALL DAXPY(LENGTH,-ONE,WORK,1,TMAT,1)
C
      DO I = 1,LENGTH
C
         TMAT(I) =   SMAT(INDSQ(I,4)) 
     *             - SMAT(INDSQ(I,3))
C
      ENDDO
C
      NCK = IT1AM(ISYMC,ISYMK) + C
C
      KOFF1 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
      NTOAIJ = MAX(NCKI(ISYAIJ),1)
      NTOTC  = MAX(NVIR(ISYMC),1)
C
      CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *           FOCKAK(NCK),NTOTC,ONE,RMAT2,1)
C
C----------------------------------
C     First contribution to Omega1.
C----------------------------------
C
      ISYMI  = MULD2H(ISYMC,ISYRES)
      ISYAKJ = MULD2H(ISYMB,ISYINT)
C
      IF ((.NOT. CC3LR) .AND. (NRHF(ISYMI) .NE. 0)) THEN
C
         IF (LWORK .LT. NCKI(ISYAKJ)) THEN
            CALL QUIT('Not enough core in CC3_ONEL')
         END IF
C
C        Construct M(ak,j) = L(ak,bj)
C        ---------------------------
C
         DO 100 ISYMJ = 1,NSYM
C
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAK = MULD2H(ISYMJ,ISYAKJ)
C
            DO 110 J = 1,NRHF(ISYMJ)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO 120 NAK = 1,NT1AM(ISYMAK)
C
                  NAKBJ = IT2AM(ISYMAK,ISYMBJ) + INDEX(NAK,NBJ)
                  NAKJ  = ICKI(ISYMAK,ISYMJ)
     *                  + NT1AM(ISYMAK)*(J - 1) + NAK
C
                  WORK(NAKJ) = XIAJB(NAKBJ)
C
  120          CONTINUE
  110       CONTINUE
  100    CONTINUE
C
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOAKJ = MAX(NCKI(ISYAKJ),1)
C
         KOFF1 = ISAIKJ(ISYAKJ,ISYMI) + 1
         KOFF2 = IT1AM(ISYMC,ISYMI) + C
C
         CALL DGEMV('T',NCKI(ISYAKJ),NRHF(ISYMI),ONE,TMAT(KOFF1),
     *              NTOAKJ,WORK,1,ONE,OMEGA1(KOFF2),NTOTC)
C
      ENDIF
C
C-----------------------------------
C     Second contribution to Omega2.
C-----------------------------------
C
      ISYMK  = MULD2H(ISYMB,ISYINT)
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYMBC,ISYMIM)
      ISYAIJ = MULD2H(ISYMK,JSAIKJ)
      ISYMBK = MULD2H(ISYMB,ISYMK)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CC3_ONEL')
      END IF
C
c     CALL CC_GATHER(LENGTH,TMAT,SMAT,INDSQ(1,5))
c     CALL DAXPY(LENGTH,-ONE,SMAT,1,TMAT,1)
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(INDSQ(I,5)) 
     *             - SMAT(I)
      ENDDO
C
      NBK = IT1AM(ISYMB,ISYMK) + B
C
      KOFF1 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
      NTOAIJ = MAX(NCKI(ISYAIJ),1)
      NTOTB  = MAX(NVIR(ISYMB),1)
C
      CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *           FOCKAK(NBK),NTOTB,ONE,RMAT1,1)
C
C-----------------------------------
C     Second contribution to Omega1.
C-----------------------------------
C
      IF (.NOT. CC3LR) THEN
C
C        Symmetry sorting if symmetry
C        ----------------------------
C
         IF (NSYM .GT. 1) THEN
            CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
            CALL DCOPY(LENGTH,WORK,1,TMAT,1)
         ENDIF
C
         ISYMKJ = MULD2H(ISYMBC,ISYINT)
         ISYMAI = ISYRES
C
C        Construct M(k,j) = L(ck,bj)
C        ---------------------------
C
         DO 200 ISYMJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMJ,ISYMKJ)
            ISYMCK = MULD2H(ISYMC,ISYMK)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO 210 J = 1,NRHF(ISYMJ)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO 220 K = 1,NRHF(ISYMK)
C
                  NKJ = IMATIJ(ISYMK,ISYMJ)+ NRHF(ISYMK)*(J - 1) + K
                  NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
                  NCKBJ = IT2AM(ISYMCK,ISYMBJ) + INDEX(NCK,NBJ)
C
                  WORK(NKJ) = XIAJB(NCKBJ)
C
  220          CONTINUE
  210       CONTINUE
  200    CONTINUE
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
C
         KOFF1 = ISAIKL(ISYMAI,ISYMKJ) + 1
C
         CALL DGEMV('N',NT1AM(ISYMAI),NMATIJ(ISYMKJ),ONE,TMAT(KOFF1),
     *              NTOTAI,WORK,1,ONE,OMEGA1,1)
C
      ENDIF
C
C----------------------------------
C     Third contribution to Omega2.
C----------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAKIJ = MULD2H(ISYMBC,ISYMIM)
      ISYMIJ = MULD2H(ISYMBC,ISYRES)
      ISYMAK = MULD2H(JSAKIJ,ISYMIJ)
C
      LENGTH = NCKIJ(JSAKIJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CC3_ONEL')
      END IF
C
c     CALL CC_GATHER(LENGTH,TMAT,SMAT,INDSQ(1,1))
c     CALL DAXPY(LENGTH,-ONE,SMAT,1,TMAT,1)
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(INDSQ(I,1)) 
     *             - SMAT(I)
      ENDDO
C
C     Symmetry sorting if symmetry
C     ----------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      NTOTAK = MAX(NT1AM(ISYMAK),1)
      NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
C
      KOFF1 = ISAIKL(ISYMAK,ISYMIJ) + 1
C
      CALL DGEMV('T',NT1AM(ISYMAK),NMATIJ(ISYMIJ),ONE,TMAT(KOFF1),
     *           NTOTAK,FOCKAK,1,ZERO,WORK,1)
C
      DO 300 ISYMJ = 1,NSYM
C
         ISYMI  = MULD2H(ISYMIJ,ISYMJ)
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMCI = MULD2H(ISYMC,ISYMI)
C
         DO 310 J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            IF (ISYMCI .EQ. ISYMBJ) THEN
C
               DO 320 I = 1,NRHF(ISYMI)
C
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NCI = IT1AM(ISYMC,ISYMI)  + NVIR(ISYMC)*(I - 1) + C
C
                  IF (NCI .EQ. NBJ) WORK(NIJ) = TWO*WORK(NIJ)
C
                  NCIBJ = IT2AM(ISYMCI,ISYMBJ) + INDEX(NCI,NBJ)
C
                  OMEGA2(NCIBJ) = OMEGA2(NCIBJ) + WORK(NIJ)
C
  320          CONTINUE
C
            ELSE IF (ISYMCI .LT. ISYMBJ) THEN
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NCI = IT1AM(ISYMC,ISYMI)  + NVIR(ISYMC)*(I - 1) + C
C
                  NCIBJ = IT2AM(ISYMCI,ISYMBJ)
     *                  + NT1AM(ISYMCI)*(NBJ-1) + NCI
C
                  OMEGA2(NCIBJ) = OMEGA2(NCIBJ) + WORK(NIJ)
C
  330          CONTINUE
C
            ELSE IF (ISYMBJ .LT. ISYMCI) THEN
C
               DO 340 I = 1,NRHF(ISYMI)
C
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NCI = IT1AM(ISYMC,ISYMI)  + NVIR(ISYMC)*(I - 1) + C
C
                  NCIBJ = IT2AM(ISYMBJ,ISYMCI)
     *                  + NT1AM(ISYMBJ)*(NCI-1) + NBJ
C
                  OMEGA2(NCIBJ) = OMEGA2(NCIBJ) + WORK(NIJ)
C
  340          CONTINUE
C
            ENDIF
C
  310    CONTINUE
C
  300 CONTINUE
C
C----------------------------------
C     Third contribution to Omega1.
C----------------------------------
C
  333 CONTINUE
C
      IF (.NOT. CC3LR) THEN
C
         ISYMKJ = MULD2H(ISYMBC,ISYINT)
         ISYMAI = ISYRES
C
C        Construct M(k,j) = L(ck,bj)
C        ---------------------------
C
         DO 400 ISYMJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMJ,ISYMKJ)
            ISYMCK = MULD2H(ISYMC,ISYMK)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO 410 J = 1,NRHF(ISYMJ)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO 420 K = 1,NRHF(ISYMK)
C
                  NKJ = IMATIJ(ISYMK,ISYMJ)+ NRHF(ISYMK)*(J - 1) + K
                  NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
                  NCKBJ = IT2AM(ISYMCK,ISYMBJ) + INDEX(NCK,NBJ)
C
                  WORK(NKJ) = XIAJB(NCKBJ)
C
  420          CONTINUE
  410       CONTINUE
  400    CONTINUE
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
C
         KOFF1 = ISAIKL(ISYMAI,ISYMKJ) + 1
C
         CALL DGEMV('N',NT1AM(ISYMAI),NMATIJ(ISYMKJ),ONE,TMAT(KOFF1),
     *              NTOTAI,WORK,1,ONE,OMEGA1,1)
C
      ENDIF
C
      CALL QEXIT('CC3_ONEL')
C
C
      RETURN
      END
C  /* Deck cc3lr_mfock */
      SUBROUTINE CC3LR_MFOCK(FOCK,T1AM,XIAJB,ISYMTR)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Calculate modified Fock matrix iused in CC3LR.
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION FOCK(*),T1AM(*),XIAJB(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3LR_MFOCK')
C
      ISYMBJ = ISYMTR
      ISYMAI = MULD2H(ISYMBJ,ISYMOP)
C
      IF (ISYMAI .EQ. ISYMBJ) THEN
C
         CALL DZERO(FOCK,NT1AM(ISYMAI))
C
         DO 100 NBJ = 1,NT1AM(ISYMBJ)
            DO 110 NAI = 1,NT1AM(ISYMAI)
C
               NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
               FOCK(NAI) = FOCK(NAI) + XIAJB(NAIBJ)*T1AM(NBJ)
C
  110       CONTINUE
  100    CONTINUE
C
      ENDIF
C
      IF (ISYMAI .GT. ISYMBJ) THEN
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
C
         KOFF1 = IT2AM(ISYMAI,ISYMBJ) + 1
C
         CALL DGEMV('N',NT1AM(ISYMAI),NT1AM(ISYMBJ),ONE,XIAJB(KOFF1),
     *              NTOTAI,T1AM,1,ZERO,FOCK,1)
C
      ENDIF
C
      IF (ISYMAI .LT. ISYMBJ) THEN
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
C
         KOFF1 = IT2AM(ISYMAI,ISYMBJ) + 1
C
         CALL DGEMV('T',NT1AM(ISYMAI),NT1AM(ISYMBJ),ONE,XIAJB(KOFF1),
     *              NTOTAI,T1AM,1,ZERO,FOCK,1)
C
      ENDIF
C
      CALL QEXIT('CC3LR_MFOCK')
C
      RETURN
      END
C  /* Deck cc3lr_pint */
      SUBROUTINE CC3LR_PINT(P2INT,P1INT,T1AM,XIAJB,WORK,LWORK,
     *                       ISYMD,D,ISYMTR)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Calculate modified Fock matrix iused in CC3LR.
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION P1INT(*),P2INT(*),T1AM(*),XIAJB(*),WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3LR_PINT')
C
      ISCKAD = MULD2H(ISYMOP,ISYMTR)
      ISYCKA = MULD2H(ISYMD,ISCKAD)
C
      DO 100 ISYMDL = 1,NSYM
C
         ISYML  = MULD2H(ISYMDL,ISYMD)
         ISYMA  = MULD2H(ISYML,ISYMTR)
         ISYMCK = MULD2H(ISYMDL,ISYMOP)
C
         IF (LWORK .LT. NT1AM(ISYMCK)*NRHF(ISYML)) THEN
            CALL QUIT('Insuffiecient work space in CC3LR_PINT')
         ENDIF
C
         IF (ISYMCK .NE. ISYMDL)
     *      CALL QUIT('Symmetry inconsistent in CC3LR_PINT')
C
C-------------------------------------------------------
C        Unpack part of integrals before multiplication.
C-------------------------------------------------------
C
         DO 110 L = 1,NRHF(ISYML)
C
            NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L - 1) + D
C
            DO 120 NCK = 1,NT1AM(ISYMCK)
C
               NCKDL = IT2AM(ISYMCK,ISYMDL) + INDEX(NCK,NDL)
               NCKL  = NT1AM(ISYMCK)*(L - 1) + NCK
C
               WORK(NCKL) = XIAJB(NCKDL)
C
  120       CONTINUE
  110    CONTINUE
C
C-------------------------------------------
C        Calculate transformed P1 integrals.
C-------------------------------------------
C
         NVIRA  = MAX(NVIR(ISYMA),1)
         NTOTCK = MAX(NT1AM(ISYMCK),1)
C
         KOFF1 = IT1AM(ISYMA,ISYML) + 1
         KOFF2 = ICKATR(ISYMCK,ISYMA) + 1
C
         CALL DGEMM('N','T',NT1AM(ISYMCK),NVIR(ISYMA),NRHF(ISYML),
     *              ONE,WORK,NTOTCK,T1AM(KOFF1),NVIRA,ZERO,P1INT(KOFF2),
     *              NTOTCK)
C
  100 CONTINUE
C
C
      DO 200 ISYMDL = 1,NSYM
C
         ISYMCK = MULD2H(ISYMDL,ISYMOP)
         ISYML  = MULD2H(ISYMD,ISYMDL)
         ISYMA  = MULD2H(ISYML,ISYMTR)
C
         IF (LWORK .LT. NT1AM(ISYMCK)*NRHF(ISYML)) THEN
            CALL QUIT('Insuffiecient work space in CC3LR_PINT')
         ENDIF
C
         IF (ISYMCK .NE. ISYMDL)
     *      CALL QUIT('Symmetry inconsistent in CC3LR_PINT')
C
C-------------------------------------------------------
C        Unpack part of integrals before multiplication.
C-------------------------------------------------------
C
         DO 210 L = 1,NRHF(ISYML)
C
            DO 220 ISYMK = 1,NSYM
C
               ISYMC  = MULD2H(ISYMK,ISYMCK)
               ISYMDK = MULD2H(ISYMD,ISYMK)
               ISYMCL = MULD2H(ISYMC,ISYML)
C
               DO 230 K = 1,NRHF(ISYMK)
C
                 NDK = IT1AM(ISYMD,ISYMK) + NVIR(ISYMD)*(K - 1) + D
C
                  DO 240 C = 1,NVIR(ISYMC)
C
                     NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
                     NCL = IT1AM(ISYMC,ISYML) + NVIR(ISYMC)*(L - 1) + C
C
                     NDKCL = IT2AM(ISYMDK,ISYMCL) + INDEX(NDK,NCL)
                     NCKL  = NT1AM(ISYMCK)*(L - 1) + NCK
C
                     WORK(NCKL) = XIAJB(NDKCL)
C
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
C
C-------------------------------------------
C        Calculate transformed P2 integrals.
C-------------------------------------------
C
         NVIRA  = MAX(NVIR(ISYMA),1)
         NTOTCK = MAX(NT1AM(ISYMCK),1)
C
         KOFF1 = IT1AM(ISYMA,ISYML) + 1
         KOFF2 = ICKATR(ISYMCK,ISYMA) + 1
C
         CALL DGEMM('N','T',NT1AM(ISYMCK),NVIR(ISYMA),NRHF(ISYML),
     *              ONE,WORK,NTOTCK,T1AM(KOFF1),NVIRA,ZERO,P2INT(KOFF2),
     *              NTOTCK)
C
  200 CONTINUE

C
C---------------------------
C     Scale integrals by -1.
C---------------------------
C
      CALL DSCAL(NCKATR(ISYCKA),-ONE,P1INT,1)
      CALL DSCAL(NCKATR(ISYCKA),-ONE,P2INT,1)
C
      CALL QEXIT('CC3LR_PINT')
C
      RETURN
      END
C  /* Deck cc3lr_qint */
      SUBROUTINE CC3LR_QINT(QINT,T1AM,XIAJB,WORK,LWORK,ISYMTR)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     int(clki) = int(cidj)*tdi
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION QINT(*),T1AM(*),XIAJB(*),WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3LR_QINT')
C
      DO 100 ISYMK = 1,NSYM
C
         DO 110 ISYMD = 1,NSYM
C
            ISYMDK = MULD2H(ISYMD,ISYMK)
            ISYMCL = MULD2H(ISYMDK,ISYMOP)
            ISYMI  = MULD2H(ISYMD,ISYMTR)
C
            ISYCLK = MULD2H(ISYMCL,ISYMK)
C
            NTOTCL = MAX(NT1AM(ISYMCL),1)
C
            IF (LWORK .LT. NT1AM(ISYMCL)*NVIR(ISYMD)) THEN
               CALL QUIT('Insufficient work space in CC3LR_QINT')
            ENDIF
C
            DO 120 K = 1,NRHF(ISYMK)
C
C----------------------------------------
C              Copy integrals from XIAJB.
C----------------------------------------
C
               DO 130 D = 1,NVIR(ISYMD)
C
                  NDK = IT1AM(ISYMD,ISYMK) + NVIR(ISYMD)*(K - 1) + D
C
                  DO 140 NCL = 1,NT1AM(ISYMCL)
C
                     NCLDK = IT2AM(ISYMCL,ISYMDK) + INDEX(NCL,NDK)
                     NCLD  = NT1AM(ISYMCL)*(D - 1) + NCL
C
                     WORK(NCLD) = XIAJB(NCLDK)
C
  140             CONTINUE
  130          CONTINUE
C
C------------------------------------------------
C              Contract integrals and amplitudes.
C------------------------------------------------
C
               DO 150 I = 1,NRHF(ISYMI)
C
                  KOFF1 = IT1AM(ISYMD,ISYMI) + NVIR(ISYMD)*(I - 1) + 1
                  KOFF2 = ISAIKJ(ISYCLK,ISYMI)
     *                  + NCKI(ISYCLK)*(I - 1)
     *                  + ISAIK(ISYMCL,ISYMK)
     *                  + NT1AM(ISYMCL)*(K - 1) + 1
C
                  CALL DGEMV('N',NT1AM(ISYMCL),NVIR(ISYMD),ONE,WORK,
     *                       NTOTCL,T1AM(KOFF1),1,ZERO,QINT(KOFF2),1)
C
  150          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3LR_QINT')
C
      RETURN
      END
C  /* Deck cc3_sortocc */
      SUBROUTINE CC3_SORTOCC(XOCC,WORK,LWORK,ISYINT)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Resort occupied integrals
C
C     initially stored I(ai,j,k) to finally I'(ai,k,j)
C
#include <implicit.h>
C
      DIMENSION XOCC(*)
      DIMENSION WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_SORTOCC')
C
      DO 100 ISYMK = 1,NSYM
C
         ISYAIJ = MULD2H(ISYMK,ISYINT)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
            ISYAIK = MULD2H(ISYMAI,ISYMK)
C
            DO 120 K = 1,NRHF(ISYMK)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  KOFF1 = ISAIKJ(ISYAIJ,ISYMK)
     *                  + NCKI(ISYAIJ)*(K - 1)
     *                  + ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + 1
C
                  KOFF2 = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ISAIK(ISYMAI,ISYMK)
     *                  + NT1AM(ISYMAI)*(K - 1) + 1
C
                  CALL DCOPY(NT1AM(ISYMAI),XOCC(KOFF1),1,
     *                       WORK(KOFF2),1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL DCOPY(NCKIJ(ISYINT),WORK,1,XOCC,1)
C
      CALL QEXIT('CC3_SORTOCC')
C
      RETURN
      END
C  /* Deck cc3_sint */
      SUBROUTINE CC3_SINT(XLAMD,WORK,LWORK,ISYINT,LUDELD,FNDELD,
     *                    LUDKBC,FNDKBC)
C
C     Henrik Koch and Poul Joergensen Jan 1995
C
C     Construct s-matrix integrals from q-matrix integrals.
C
C     Ove Christiansen 11-1-1996: ISYINT; symmetry of integrals.
C                                 (C1 transformed integrals)
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION XLAMD(*),WORK(LWORK)
C
      CHARACTER*(*) FNDELD, FNDKBC
C
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_SINT')
C
C----------------------------------
C     Read and transform integrals.
C----------------------------------
C
      DO 100 ISYMB = 1,NSYM
C
         IF (NVIR(ISYMB) .EQ. 0) GOTO 100
C
         ISYCKD = MULD2H(ISYINT,ISYMB)
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KSORT = 1
         KEND1 = KSORT + NCKATR(ISYCKD)
         LWRK1 = LWORK - KEND1
C
         MINWRK = NCKA(ISYCKD) + NCKATR(ISYCKD)
C
         IF (LWRK1 .LT. MINWRK) THEN
            CALL QUIT('Insufficient space in CC3_SINT')
         END IF
C
C------------------------
C        Batch structure.
C------------------------
C
         NDISB = MIN(LWRK1/MINWRK,NVIR(ISYMB))
         NUMB  = NDISB
COMMENT COMMENT
COMMENT COMMENT Old version gave wrong results
COMMENT COMMENT if not all batches were in memory
COMMENT COMMENT
C         NBATB =(NDISB - 1)/NVIR(ISYMB) + 1
         NBATB =(NVIR(ISYMB) - 1)/NDISB + 1
COMMENT COMMENT
COMMENT COMMENT
COMMENT COMMENT
C
         ILASTB = 0
         DO 110 IBATB = 1,NBATB
C
            IFRSTB = ILASTB + 1
            ILASTB = ILASTB + NDISB
            IF (ILASTB .GT. NVIR(ISYMB)) THEN
               ILASTB = NVIR(ISYMB)
               NUMB   = ILASTB - IFRSTB + 1
            END IF
C
C------------------------------
C           Dynamic allocation.
C------------------------------
C
            KTRANS = KEND1
            KEND2  = KTRANS + NUMB*NCKATR(ISYCKD)
C
            KUNTRA = KEND2
            KEND3  = KUNTRA + NUMB*NCKA(ISYCKD)
            LWRK3  = LWORK  - KEND3
C
            IF (LWRK3 .LT.0) THEN
               CALL QUIT('Insufficient space in CC3_SINT')
            END IF
C
C---------------------------------------------------
C           Read integrals I(ck,del,b) = (ck|b del).
C---------------------------------------------------
C
            IOFF = ICKAD(ISYCKD,ISYMB)
     *           + NCKA(ISYCKD)*(IFRSTB - 1) + 1
            LENG = NCKA(ISYCKD)*NUMB
C
            IF (NCKA(ISYCKD) .GT. 0) THEN
               CALL GETWA2(LUDELD,FNDELD,WORK(KUNTRA),IOFF,LENG)
            ENDIF
C
C---------------------------------------
C           Loop over orbitals in batch.
C---------------------------------------
C
            DO 120 B = IFRSTB,ILASTB
C
               IB = B - IFRSTB + 1
C
C--------------------------------
C              Integral adresses.
C--------------------------------
C
               KSCR1 = KUNTRA + NCKA(ISYCKD)*(IB - 1)
               KSCR2 = KTRANS + NCKATR(ISYCKD)*(IB - 1)
C
C---------------------------------------------
C              Transform integrals with XLAMD.
C---------------------------------------------
C
               DO 130 ISYMD = 1,NSYM
C
                  ISYMCK = MULD2H(ISYMD,ISYCKD)
C
                  NTOTCK = MAX(NT1AM(ISYMCK),1)
                  NBASD  = MAX(NBAS(ISYMD),1)
C
                  KOFF1  = KSCR1 + ICKA(ISYMCK,ISYMD)
                  KOFF2  = ILMVIR(ISYMD) + 1
                  KOFF3  = KSCR2 + ICKATR(ISYMCK,ISYMD)
C
                  CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMD),
     *                       NBAS(ISYMD),ONE,WORK(KOFF1),NTOTCK,
     *                       XLAMD(KOFF2),NBASD,ZERO,
     *                       WORK(KOFF3),NTOTCK)
C
  130          CONTINUE
C
C--------------------------------------------------------------
C              Sort integrals I(ck,d,b) to I(dk,c,b) for fix b.
C--------------------------------------------------------------
C
               CALL CC3_CKDSOR(WORK(KSCR2),WORK(KSORT),
     *                         NCKATR(ISYCKD),ISYCKD)
C
  120       CONTINUE
C
C---------------------------------------------------------------
C           Write to disk MO integrals I(dk,c,#b) as I(dk,#b,c).
C---------------------------------------------------------------
C
            KSRT  = KEND2
            KEND3 = KSRT  + NUMB*NCKATR(ISYCKD)
            LWRK3 = LWORK - KEND3
C
            IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient space in CC3_SINT')
            END IF
C
            KOFF2 = KSRT
            KOFF3 = KSRT
            DO 140 ISYMC = 1,NSYM
C
               ISYMDK = MULD2H(ISYCKD,ISYMC)
               ISYDKB = MULD2H(ISYMDK,ISYMB)
C
               LENGT1 = NT1AM(ISYMDK)
               LENGT2 = NT1AM(ISYMDK)*NUMB
C
               DO 150 C = 1,NVIR(ISYMC)
C
C----------------------
C                 Sort.
C----------------------
C
                  DO 160 IB = 1,NUMB
C
                     KOFF1 = KTRANS
     *                     + NCKATR(ISYCKD)*(IB - 1)
     *                     + ICKATR(ISYMDK,ISYMC)
     *                     + NT1AM(ISYMDK)*(C - 1)
C
                     CALL DCOPY(LENGT1,WORK(KOFF1),1,WORK(KOFF2),1)
C
                     KOFF2 = KOFF2 + LENGT1
C
  160             CONTINUE
C
C----------------------
C                 Dump.
C----------------------
C
                  IOFF  = ICKBD(ISYDKB,ISYMC)
     *                  + NCKATR(ISYDKB)*(C - 1)
     *                  + ICKATR(ISYMDK,ISYMB)
     *                  + NT1AM(ISYMDK)*(IFRSTB - 1) + 1
C
                  IF (LENGT2 .GT. 0) THEN
                     CALL PUTWA2(LUDKBC,FNDKBC,WORK(KOFF3),IOFF,
     *                           LENGT2)
                     KOFF3 = KOFF3 + LENGT2
                  ENDIF
C
  150          CONTINUE
  140       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_SINT')
C
      RETURN
      END
C  /* Deck cc3_diag */
      SUBROUTINE CC3_DIAG(DIAG,FOCKD,JSAIKJ)
C
C     Henrik Koch Feb 1995
C
C     Construct diagonal matrix D(aikj) = F(a) - F(k) - F(j) - F(i)
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION DIAG(*),FOCKD(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_DIAG')
C
      KOFF = 0
      DO 100 ISYMJ = 1,NSYM
C
         ISYAIK = MULD2H(JSAIKJ,ISYMJ)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            NJ = IORB(ISYMJ) + J
C
            DO 120 ISYMK = 1,NSYM
C
               ISYMAI = MULD2H(ISYAIK,ISYMK)
C
               DO 130 K = 1,NRHF(ISYMK)
C
                  NK = IORB(ISYMK) + K
C
                  DO 140 ISYMI = 1,NSYM
C
                     ISYMA  = MULD2H(ISYMAI,ISYMI)
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        NI = IORB(ISYMI) + I
C
                        DO 160 A = 1,NVIR(ISYMA)
C
                           NA = IORB(ISYMA) + NRHF(ISYMA) + A
C
                           KOFF = KOFF + 1
C
                           DIAG(KOFF) = FOCKD(NA) - FOCKD(NI)
     *                                - FOCKD(NJ) - FOCKD(NK)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_DIAG')
C
      RETURN
      END
C  /* Deck cc3_indsq */
      SUBROUTINE CC3_INDSQ(INDSQ,LENGTH,ISCKIJ)
C
C     Henrik Koch Feb 1995
C
C     Construct index arrays used triples code.
C
C     index(ckij)  = (cikj) (1)
C                    (cijk) (2)
C                    (ckji) (3)
C                    (cjki) (4)
C                    (cjik) (5)
C
C     index(ck,ij) = (ckij) (6)
C
#include <implicit.h>
C
      DIMENSION INDSQ(LENGTH,6)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_INDSQ')
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYCKI = MULD2H(ISYMJ,ISCKIJ)
         ISYCIK = MULD2H(ISYMJ,ISCKIJ)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYCKJ = MULD2H(ISYMI,ISCKIJ)
            ISYCJK = MULD2H(ISYMI,ISCKIJ)
            ISYMCK = MULD2H(ISYMI,ISYCKI)
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
            DO 120 ISYMK = 1,NSYM
C
               ISYCIJ = MULD2H(ISYMK,ISCKIJ)
               ISYCJI = MULD2H(ISYMK,ISCKIJ)
               ISYMC  = MULD2H(ISYMK,ISYMCK)
               ISYMCI = MULD2H(ISYMC,ISYMI)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                   DO 140 I = 1,NRHF(ISYMI)
C
                      NIJ = IMATIJ(ISYMI,ISYMJ)
     *                    + NRHF(ISYMI)*(J - 1) + I
C
                      DO 150 K = 1,NRHF(ISYMK)
C
                         KOFFA = ISAIKJ(ISYCKI,ISYMJ)
     *                         + NCKI(ISYCKI)*(J - 1)
     *                         + ISAIK(ISYMCK,ISYMI)
     *                         + NT1AM(ISYMCK)*(I - 1)
     *                         + IT1AM(ISYMC,ISYMK)
     *                         + NVIR(ISYMC)*(K - 1)
C
                         KOFFB = ISAIKL(ISYMCK,ISYMIJ)
     *                         + NT1AM(ISYMCK)*(NIJ - 1)
     *                         + IT1AM(ISYMC,ISYMK)
     *                         + NVIR(ISYMC)*(K - 1)
C
                         KOFF1 = ISAIKJ(ISYCIK,ISYMJ)
     *                         + NCKI(ISYCIK)*(J - 1)
     *                         + ISAIK(ISYMCI,ISYMK)
     *                         + NT1AM(ISYMCI)*(K - 1)
     *                         + IT1AM(ISYMC,ISYMI)
     *                         + NVIR(ISYMC)*(I - 1)
C
                         KOFF2 = ISAIKJ(ISYCIJ,ISYMK)
     *                         + NCKI(ISYCIJ)*(K - 1)
     *                         + ISAIK(ISYMCI,ISYMJ)
     *                         + NT1AM(ISYMCI)*(J - 1)
     *                         + IT1AM(ISYMC,ISYMI)
     *                         + NVIR(ISYMC)*(I - 1)
C
                         KOFF3 = ISAIKJ(ISYCKJ,ISYMI)
     *                         + NCKI(ISYCKJ)*(I - 1)
     *                         + ISAIK(ISYMCK,ISYMJ)
     *                         + NT1AM(ISYMCK)*(J - 1)
     *                         + IT1AM(ISYMC,ISYMK)
     *                         + NVIR(ISYMC)*(K - 1)
C
                         KOFF4 = ISAIKJ(ISYCJK,ISYMI)
     *                         + NCKI(ISYCJK)*(I - 1)
     *                         + ISAIK(ISYMCJ,ISYMK)
     *                         + NT1AM(ISYMCJ)*(K - 1)
     *                         + IT1AM(ISYMC,ISYMJ)
     *                         + NVIR(ISYMC)*(J - 1)
C
                         KOFF5 = ISAIKJ(ISYCJI,ISYMK)
     *                         + NCKI(ISYCJI)*(K - 1)
     *                         + ISAIK(ISYMCJ,ISYMI)
     *                         + NT1AM(ISYMCJ)*(I - 1)
     *                         + IT1AM(ISYMC,ISYMJ)
     *                         + NVIR(ISYMC)*(J - 1)
C
                         DO 160 C = 1,NVIR(ISYMC)
C
                            INDSQ(KOFFA + C,1) = KOFF1 + C
                            INDSQ(KOFFA + C,2) = KOFF2 + C
                            INDSQ(KOFFA + C,3) = KOFF3 + C
                            INDSQ(KOFFA + C,4) = KOFF4 + C
                            INDSQ(KOFFA + C,5) = KOFF5 + C
C
                            INDSQ(KOFFB + C,6) = KOFFA + C
C
  160                    CONTINUE
  150                 CONTINUE
  140              CONTINUE
  130           CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_INDSQ')
C
      RETURN
      END
C  /* Deck cc3_t2tp */
      SUBROUTINE CC3_T2TP(T2TP,T2AM,ISYMT2)
C
C     Henrik Koch Feb 1995
C
C     Reorder t2 amplitudes as:
C
C     t2tp(aijb) = t2am(ai,bj)
C
C
#include <implicit.h>
C
      DIMENSION T2AM(*),T2TP(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_T2TP')
C
      DO 100 ISYMB = 1,NSYM
C
         ISYAIJ = MULD2H(ISYMB,ISYMT2)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAI = MULD2H(ISYMBJ,ISYMT2)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ   = IT1AM(ISYMB,ISYMJ)
     *                  + NVIR(ISYMB)*(J - 1) + B
C
                  KOFF1 = IT2SQ(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ - 1) + 1
C
                  KOFF2 = IT2SP(ISYAIJ,ISYMB)
     *                  + NCKI(ISYAIJ)*(B - 1)
     *                  + ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + 1
C
                  CALL DCOPY(NT1AM(ISYMAI),T2AM(KOFF1),1,T2TP(KOFF2),1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_T2TP')
C
      RETURN
      END
C  /* Deck cc_gather */
      SUBROUTINE CC_GATHER(LENGTH,VEC1,VEC2,INDEX)
C
C     Henrik Koch Feb 1995
C
C     Fortran version of Cray's gather routine.
C
C
#include <implicit.h>
C
      DIMENSION VEC1(*),VEC2(*),INDEX(*)
C
      CALL QENTER('CC_GATHER')
C
#if defined (SYS_CRAY)
      CALL GATHER(LENGTH,VEC1,VEC2,INDEX)
#else
      DO 100 I = 1,LENGTH
C
         VEC1(I) = VEC2(INDEX(I))
C
  100 CONTINUE
#endif
C
      CALL QEXIT('CC_GATHER')
C
      RETURN
      END
C  /* Deck cc3_racc */
      SUBROUTINE CC3_RACC(OMEGA2,RMAT,ISYMB,B,ISYRES)
C
C     Written by HK and ASM Maj 1995.
C
C     Purpose :  Accumulate the R matrix into the Omega vector.
C
C     OC: 9-1-1996 general symmetry, ISYRES is symmetry of resulting vector.
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0)
C
      DIMENSION OMEGA2(*),RMAT(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_RACC')
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMJ,ISYAIJ)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
         IF (NT1AM(ISYMAI) .LE. 0 ) GO TO 100
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            IF (ISYMAI .EQ. ISYMBJ) THEN
C
               NBJJ = ISAIK(ISYMBJ,ISYMJ)
     *              + NT1AM(ISYMBJ)*(J - 1) + NBJ
C
               RMAT(NBJJ) = TWO * RMAT(NBJJ)
C
               DO 230 NAI = 1,NT1AM(ISYMAI)
C
                  KOFF1 = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  KOFF2 = ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + NAI
C
                  OMEGA2(KOFF1) = OMEGA2(KOFF1) + RMAT(KOFF2)
C
  230          CONTINUE
C
            ENDIF
C
            IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO 240 NAI = 1,NT1AM(ISYMAI)
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ-1) + NAI
C
                  KOFF6 = ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) +  NAI
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) + RMAT(KOFF6)
C
  240          CONTINUE
C
            ENDIF
C
            IF (ISYMAI .GT. ISYMBJ) THEN
C
               DO 250 NAI = 1,NT1AM(ISYMAI)
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMBJ)*(NAI-1) + NBJ
C
                  KOFF6 = ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + NAI
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) + RMAT(KOFF6)
C
  250          CONTINUE
C
            ENDIF
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_RACC')
C
      RETURN
      END
C  /* Deck cc3_trocc */
      SUBROUTINE CC3_TROCC(XINT,TRINT,XLAMDH,WORK,LWORK,ISYINT)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Transform (ai|j delta) integrals to (ai|j k) and sort as (ij,k,a)
C
C     Debugged for ISYINT .NE. 1 Ove Christiansen 11-1-1996
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*),TRINT(*),XLAMDH(*),WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      CALL QENTER('CC3_TROCC')
C
      IF (LWORK .LT. NTRAOC(ISYINT)) THEN
         CALL QUIT('Insufficient space in CC3_TROCC')
      END IF
C
C     write out integral norm
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTOTOC(ISYINT),XINT,1,XINT,1)
         WRITE(LUPRI,*) 'In CC3_TROCC: Norm of INT = ',XTROC
      ENDIF
C
C     Transform
C
      DO 100 ISYMD = 1,NSYM
C
         ISYMK  = ISYMD
         ISYAIJ = MULD2H(ISYMD,ISYINT)
C
         NTOIAJ = MAX(NCKI(ISYAIJ),1)
         NBASD  = MAX(NBAS(ISYMD),1)
C
         KOFF1 = ICKID(ISYAIJ,ISYMD)  + 1
         KOFF2 = ILMRHF(ISYMD) + 1
         KOFF3 = ICKITR(ISYAIJ,ISYMK)  + 1
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),NBAS(ISYMD),
     *              ONE,XINT(KOFF1),NTOIAJ,XLAMDH(KOFF2),NBASD,
     *              ZERO,TRINT(KOFF3),NTOIAJ)
C
  100 CONTINUE
C
C     write out integral norm
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTRAOC(ISYINT),TRINT,1,TRINT,1)
         WRITE(LUPRI,*) 'In CC3_TROCC: Norm of transformed INT = ',XTROC
      ENDIF
C
C
C     Intechange j and k
C
      DO 200 ISYMK = 1,NSYM
C
         ISYAIJ = MULD2H(ISYMK,ISYINT)
C
         DO 210 ISYMJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
            ISYAIK = MULD2H(ISYMAI,ISYMK)
C
            DO 220 K = 1,NRHF(ISYMK)
C
               DO 230 J = 1,NRHF(ISYMJ)
C
                  KOFF1 = ISAIKJ(ISYAIJ,ISYMK)
     *                  + NCKI(ISYAIJ)*(K - 1)
     *                  + ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + 1
C
                  KOFF2 = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ISAIK(ISYMAI,ISYMK)
     *                  + NT1AM(ISYMAI)*(K - 1) + 1
C
                  CALL DCOPY(NT1AM(ISYMAI),TRINT(KOFF1),1,
     *                       WORK(KOFF2),1)
C
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
C     write out integral norm
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTRAOC(ISYINT),TRINT,1,TRINT,1)
         WRITE(LUPRI,*) 'In CC3_TROCC: Norm of INT interchaged = ',XTROC
      ENDIF
C
C     Resort
C
      DO 300 ISYMK = 1,NSYM
C
         ISYAIJ = MULD2H(ISYMK,ISYINT)
C
         DO 310 ISYMJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
C
            DO 320 ISYMI = 1,NSYM
C
               ISYMA  = MULD2H(ISYMAI,ISYMI)
               ISYMJI = MULD2H(ISYMJ,ISYMI)
               ISYJIK = MULD2H(ISYMJI,ISYMK)
C
               DO 330 K = 1,NRHF(ISYMK)
C
                  DO 340 J = 1,NRHF(ISYMJ)
C
                     DO 350 I = 1,NRHF(ISYMI)
C
                        NTOJIK = NMAJIK(ISYJIK)
C
                        KOFF1 = ISAIKJ(ISYAIJ,ISYMK)
     *                        + NCKI(ISYAIJ)*(K - 1)
     *                        + ISAIK(ISYMAI,ISYMJ)
     *                        + NT1AM(ISYMAI)*(J - 1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I - 1) + 1
                        KOFF2 = ISJIKA(ISYJIK,ISYMA)
     *                        + ISJIK(ISYMJI,ISYMK)
     *                        + NMATIJ(ISYMJI)*(K - 1)
     *                        + IMATIJ(ISYMJ,ISYMI)
     *                        + NRHF(ISYMJ)*(I - 1) + J
C
                        CALL DCOPY(NVIR(ISYMA),WORK(KOFF1),1,
     *                             TRINT(KOFF2),NTOJIK)
C
  350                CONTINUE
  340             CONTINUE
  330          CONTINUE
C
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
C
C     write out integral norm
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTRAOC(ISYINT),TRINT,1,TRINT,1)
         WRITE(LUPRI,*) 'In CC3_TROCC: Norm of INT Resorted  = ',XTROC
      ENDIF
C
C
      CALL QEXIT('CC3_TROCC')
C
      RETURN
      END
C  /* Deck cc3_filop */
      SUBROUTINE CC3_FILOP()
C
C     Ove Christiansen 30-11-1995
C
C     Open files for storage of integrals in triple routines
C
#include <implicit.h>
C
#include <priunit.h>
#include <iratdef.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_FILOP')
C
      IF (IPRINT .GT. 15) CALL AROUND('Open  files for triples ')
C
C      CALL WOPEN2(LUTOC,FNTOC,64,0)
C      CALL WOPEN2(LU3VI,FN3VI,64,0)
C      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
C      CALL WOPEN2(LU3VI3,FN3VI3,64,0)
C      CALL WOPEN2(LU3VI4,FN3VI4,64,0)
C      CALL WOPEN2(LU3VI5,FN3VI5,64,0)
C      CALL WOPEN2(LU3SRT,FN3SRT,64,0)
C      CALL WOPEN2(LUDELD,FNDELD,64,0)
C      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
C      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
C      CALL WOPEN2(LU3SRT2,FN3SRT2,64,0)
C      CALL WOPEN2(LU3SRT3,FN3SRT3,64,0)
C      CALL WOPEN2(LUCKJD2,FNCKJD2,64,0)
C      CALL WOPEN2(LUCKJD3,FNCKJD3,64,0)
C      CALL WOPEN2(LUDELD2,FNDELD2,64,0)
C      CALL WOPEN2(LUDELD3,FNDELD3,64,0)
C      CALL WOPEN2(LUDELD4,FNDELD4,64,0)
C      CALL WOPEN2(LUDELD5,FNDELD5,64,0)
C      CALL WOPEN2(LUDKBC2,FNDKBC2,64,0)
C      CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
C      CALL WOPEN2(LUDKBC4,FNDKBC4,64,0)
C      CALL WOPEN2(LUDKBC5,FNDKBC5,64,0)
C      CALL WOPEN2(LU3FOP,FN3FOP,64,0)
C      CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
C
      CALL QEXIT('CC3_FILOP')
C
      RETURN
      END
C  /* Deck cc3_filcl */
      SUBROUTINE CC3_FILCL()
C
C     Ove Christiansen 30-11-1995
C
C     Close files for storage of integrals in triple routines
C
#include <implicit.h>
C
#include <priunit.h>
#include <iratdef.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC3_FILCL')
C
      IF (IPRINT .GT. 15) CALL AROUND('Close files for triples ')
C
C      CALL WCLOSE2(LUTOC,FNTOC,'DELETE')
C      CALL WCLOSE2(LU3VI,FN3VI,'DELETE')
C      CALL WCLOSE2(LU3VI2,FN3VI2,'DELETE')
C      CALL WCLOSE2(LU3VI3,FN3VI3,'DELETE')
C      CALL WCLOSE2(LU3VI4,FN3VI4,'DELETE')
C      CALL WCLOSE2(LU3VI5,FN3VI5,'DELETE')
C      CALL WCLOSE2(LU3SRT,FN3SRT,'DELETE')
C      CALL WCLOSE2(LUDELD,FNDELD,'DELETE')
C      CALL WCLOSE2(LUDKBC,FNDKBC,'DELETE')
C      CALL WCLOSE2(LUCKJD,FNCKJD,'DELETE')
C      CALL WCLOSE2(LU3SRT2,FN3SRT2,'DELETE')
C      CALL WCLOSE2(LU3SRT3,FN3SRT3,'DELETE')
C      CALL WCLOSE2(LUCKJD2,FNCKJD2,'DELETE')
C      CALL WCLOSE2(LUCKJD3,FNCKJD3,'DELETE')
C      CALL WCLOSE2(LUDELD2,FNDELD2,'DELETE')
C      CALL WCLOSE2(LUDELD3,FNDELD3,'DELETE')
C      CALL WCLOSE2(LUDELD4,FNDELD4,'DELETE')
C      CALL WCLOSE2(LUDELD5,FNDELD5,'DELETE')
C      CALL WCLOSE2(LUDKBC2,FNDKBC2,'DELETE')
C      CALL WCLOSE2(LUDKBC3,FNDKBC3,'DELETE')
C      CALL WCLOSE2(LUDKBC4,FNDKBC4,'DELETE')
C      CALL WCLOSE2(LUDKBC5,FNDKBC5,'DELETE')
C      CALL WCLOSE2(LU3FOP,FN3FOP,'DELETE')
C      CALL WCLOSE2(LU3FOP2,FN3FOP2,'DELETE')
C
      CALL QEXIT('CC3_FILCL')
C
      END
