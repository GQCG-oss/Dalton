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
C  /* Deck cc3_omeg3 */
      SUBROUTINE CC3_OMEG3(ECURR,OMEGA1,OMEGA2P,OMEGA2M,T1AM,ISYMT1,
     *                     T2TP,ISYMT2,FOCK,XLAMDP,XLAMDH,WORK,LWORK,
     *                     LU3SRT,FN3SRT,LUTOC,FNTOC,LUCKJD,FNCKJD)
C
C     Routine to calculate triplet excitation
C     energies for triple methods.
C     Kasper Hald Jan. 2001.
C
C     Based on cc3_omeg by:
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C     Ove Christiansen  Jan. 1996:
C
C
      IMPLICIT NONE
C
#include <iratdef.h>
#include <dummy.h>
#include <inftap.h>
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <second.h>
C
      INTEGER LWORK, ISYMTR, ISYMT1, ISYMT2, ISYRES, ISINT1, ISINT2
      INTEGER ISYMIM, KEND0, LWRK0, KEND1, LWRK1, KEND2, LWRK2
      INTEGER KEND3, LWRK3, KEND4, LWRK4, KOFF1, KOFF2
      INTEGER KFOCKD, KCMO, KFCKAK, KTROC, KTROC1, KTROC0
      INTEGER KSMAT, KQMAT, KDIAG, KINDSQ, KINDEX, KTMAT, KRMAT1, KRMAT2
      INTEGER KRMAT3, KRMAT4, KXIAJB, KINTVI, KINTOC, KTRVI, KTRVI0
      INTEGER KTRVI1, KTRVI2, KTRVI3, KTRVI4, KTRVI5, KTRVI6, KTRVI7
      INTEGER ISYMC, ISYMK, IOPT, IOFF, ISYOPE, LENGTH, LENSQ
      INTEGER ISYMB, ISYALJ, ISAIJ2, ISYMBD, ISCKIJ, ISAIJ1
      INTEGER ISYMD, ISYCKB, ISCKB1, ISCKB2, LU3SRT, LUTOC, LUCKJD
      INTEGER LUDELD, LUDKBC
!
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2P(*),OMEGA2M(*),T1AM(*),T2TP(*)
      REAL FOCK(*),XLAMDP(*),XLAMDH(*),WORK(LWORK)
      REAL XT2TP, XIAJB, XINT, XTROC, XTROC0, XTROC1
      REAL XTRVI0, XTRVI, XTRVI1, XTRVI2, XTRVI3
      REAL XFD, XDIA, XSMAT, XQMAT, XTMAT, XRMAT, RHO1N
      REAL RHO2N, ECURR
      REAL DTIME, TITRAN, TISORT, TISMAT, TIQMAT
      REAL TIOME1, TICONV, TICONO, DDOT, XMONE, ONE, MFACTOR
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2P(*),OMEGA2M(*),T1AM(*),T2TP(*)
      DOUBLE PRECISION FOCK(*),XLAMDP(*),XLAMDH(*),WORK(LWORK)
      DOUBLE PRECISION XT2TP, XIAJB, XINT, XTROC, XTROC0, XTROC1
      DOUBLE PRECISION XTRVI0, XTRVI, XTRVI1, XTRVI2, XTRVI3
      DOUBLE PRECISION XFD, XDIA, XSMAT, XQMAT, XTMAT, XRMAT, RHO1N
      DOUBLE PRECISION RHO2N, ECURR
      DOUBLE PRECISION DTIME, TITRAN, TISORT, TISMAT, TIQMAT
      DOUBLE PRECISION TIOME1, TICONV, TICONO, DDOT, XMONE, ONE, MFACTOR
#endif
      PARAMETER (XMONE = -1.0D0, ONE = 1.0D0)
C
      CHARACTER*(*) FN3SRT, FNTOC, FNCKJD
      CHARACTER*11 FNDELD, FNDKBC
      CHARACTER CDUMMY*1
C
      CALL QENTER('CC3_OMEG3')
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
      CDUMMY = ' '
C
      IPRCC = IPRINT
      ISYMTR = MULD2H(ISYMT1,ISYMT2)
      ISYRES = MULD2H(ISYMTR,ISYMOP)
      ISINT1 = MULD2H(ISYMT1,ISYMOP)
      ISINT2 = ISYMOP
      ISYMIM = ISYMOP
C
      IF (IPRINT .GT. 20 ) THEN
         WRITE(LUPRI,*) ' In CC3_OMEG3: CC1A  = ',CC1A
         WRITE(LUPRI,*) ' In CC3_OMEG3: CC1B  = ',CC1B
         WRITE(LUPRI,*) ' In CC3_OMEG3: CC3   = ',CC3
         WRITE(LUPRI,*) ' In CC3_OMEG3: CC3LR = ',CC3LR
         WRITE(LUPRI,*) ' In CC3_OMEG3: ISYMT1, ISYMT2:',ISYMT1,ISYMT2
         WRITE(LUPRI,*) ' In CC3_OMEG3: ISYRES, ISYMOP:',ISYRES,ISYMOP
         WRITE(LUPRI,*) ' In CC3_OMEG3: ISINT1, ISINT2:',ISINT1,ISINT2
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
      TICONV = 0.0D0
      TICONO = 0.0D0
      TIOME1 = 0.0D0
C
C---------------------------
C     Open files
C---------------------------
C
      LUDELD = -1
      LUDKBC = -1
      FNDELD = 'CC3_OMEG3_1'
      FNDKBC = 'CC3_OMEG3_2'
      CALL WOPEN2(LUDELD,FNDELD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
C
C---------------------------------------------------------
C     Transform and sort qmat integrals to smat integrals.
C---------------------------------------------------------
C
      CALL CC3_SORT1(WORK,LWORK,2,ISINT2,LU3SRT,FN3SRT,LUDELD,FNDELD,
     *               IDUMMY,CDUMMY,IDUMMY,CDUMMY,IDUMMY,CDUMMY)
      CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT2,LUDELD,FNDELD,
     *              LUDKBC,FNDKBC)
C
C--------------------------------------
C     Reorder the t2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP in CC3_OMEG3')
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
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_OMEG3')
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
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_OMEG3')
      END IF
C
C-----------------------------------------
C     Calculate Fock matrix used in CC3LR.
C-----------------------------------------
C
      LENGTH = IRAT*NT2AM(ISYMOP)
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
C
      ISYOPE = ISYMOP
      IOPT = 2
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYOPE,IOPT)
C
      CALL CC3LR_MFOCK(WORK(KFCKAK),T1AM,WORK(KXIAJB),ISYMT1)
C
      IF (IPRINT .GT. 55) THEN
         XFD   = DDOT(NT1AM(ISINT1),WORK(KFCKAK),1,
     *           WORK(KFCKAK),1)
         WRITE(LUPRI,*) 'Norm of MFOCK',XFD
      ENDIF
C
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
C
      CALL DZERO(WORK(KTROC),NTRAOC(ISINT1))
      CALL CC3LR_QINT(WORK(KTROC),T1AM,WORK(KXIAJB),WORK(KEND2),
     *                LWRK2,ISYMT1)
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTRAOC(ISINT1),WORK(KTROC),1,
     *                WORK(KTROC),1)
         WRITE(LUPRI,*) 'Norm of TROC after QINT',XTROC
      ENDIF
C
      IF (LWRK2 .LT. NTRAOC(ISINT1)) THEN
         CALL QUIT('Insufficient space in CC3_OMEG3')
      END IF
C
      CALL CCSDT_SRTOCC(WORK(KTROC),WORK(KEND2),ISINT1)
      CALL DCOPY(NTRAOC(ISINT1),WORK(KEND2),1,WORK(KTROC),1)
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
C----------------------------------------------------------------------
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KCMO),
     *                 WORK(KEND2),LWRK2,ISINT2)
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
            WRITE(LUPRI,*) 'In CC3_OMEG3: ISAIJ1:',ISAIJ1
            WRITE(LUPRI,*) 'In CC3_OMEG3: ISYCKB:',ISYCKB
            WRITE(LUPRI,*) 'In CC3_OMEG3: ISCKB1:',ISCKB1
            WRITE(LUPRI,*) 'In CC3_OMEG3: ISCKB2:',ISCKB2
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
         KRMAT4 = KRMAT1 + NCKI(ISAIJ1)
         KEND2  = KRMAT4 + NCKI(ISAIJ1)
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
            CALL QUIT('Insufficient space in CC3_OMEG3')
         END IF
C
         DO 110 D = 1,NVIR(ISYMD)
C
C----------------------------------------
C           Initialize the R1/R4 matrix.
C----------------------------------------
C
            CALL DZERO(WORK(KRMAT1),NCKI(ISAIJ1))
            CALL DZERO(WORK(KRMAT4),NCKI(ISAIJ1))
C
C-----------------------------------------------
C           Read virtual integrals used in s3am.
C-----------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                     NCKATR(ISCKB2))
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
C-----------------------------------------------
C           Read virtual integrals used in q3am.
C-----------------------------------------------
C
            IOFF = ICKAD(ISCKB2,ISYMD) + NCKA(ISCKB2)*(D - 1) + 1
            IF (NCKA(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                     NCKA(ISCKB2))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),XLAMDH,
     *                       ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_OMEG3')
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
     *                       WORK(KTRVI3),1)
                WRITE(LUPRI,*) 'Norm of TRVI3 ',XTRVI3
            ENDIF
C
C---------------------------------------------
C           Construct integrals used in CC3LR.
C---------------------------------------------
C
            CALL CC3LR_PINT(WORK(KTRVI),WORK(KTRVI1),T1AM,
     *                      WORK(KXIAJB),WORK(KEND4),LWRK4,
     *                      ISYMD,D,ISYMT1)
C
           IF (IPRINT .GT. 55) THEN
               XTRVI= DDOT(NCKATR(ISCKB1),WORK(KTRVI),1,
     *                      WORK(KTRVI),1)
               WRITE(LUPRI,*) 'Norm of TRVI ',XTRVI
C
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
                  WRITE(LUPRI,*) 'In CC3_OMEG3: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_OMEG3: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_OMEG3: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_OMEG3: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CC3_OMEG3: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_OMEG3: ISCKIJ:',ISCKIJ
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
               KRMAT3 = KRMAT2 + NCKI(ISAIJ2)
               KEND4  = KRMAT3 + NCKI(ISAIJ2)
               LWRK4  = LWORK  - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_OMEG3')
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
C-----------------------------------------------------
C                 Initialize the R2/R3 matrices.
C-----------------------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
                  CALL DZERO(WORK(KRMAT3),NCKI(ISAIJ2))
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
                  CALL CC3_CONVIR(WORK(KRMAT3),WORK(KSMAT),
     *                            WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                            WORK(KTRVI),WORK(KTRVI1),ISINT1,
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 after CONVIR ',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 after CONVIR',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT3),1,
     *                       WORK(KRMAT3),1)
                     WRITE(LUPRI,*) 'Norm of RMAT3 after CONVIR',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT4),1,
     *                       WORK(KRMAT4),1)
                     WRITE(LUPRI,*) 'Norm of RMAT4 after CONVIR ',XRMAT
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_CONVIR: Rho2+')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,0,1)
                     CALL AROUND('After CC3_CONVIR: Rho2-')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TICONV = TICONV   + DTIME
C
                  DTIME = SECOND()
                  CALL CC3_CONOCC3(OMEGA2P,OMEGA2M,WORK(KRMAT1),
     *                            WORK(KRMAT2),WORK(KSMAT),WORK(KTMAT),
     *                            ISYMIM,WORK(KTROC),WORK(KTROC1),
     *                            ISINT1,WORK(KEND4),LWRK4,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 after CONOCC3',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 after CONOCC3',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT3),1,
     *                       WORK(KRMAT3),1)
                     WRITE(LUPRI,*) 'Norm of RMAT3 after CONOCC3',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT4),1,
     *                       WORK(KRMAT4),1)
                     WRITE(LUPRI,*) 'Norm of RMAT4 after CONOCC3',XRMAT
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_CONOCC',
     *                               RHO1N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
                     WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_CONOCC3',
     *                               RHO2N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
                     WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_CONOCC3',
     *                               RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_CONOCC: Rho2+')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,0,1)
                     CALL AROUND('After CC3_CONOCC: Rho2-')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
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
                  CALL CC3_ONEL3(OMEGA1,OMEGA2P,OMEGA2M,WORK(KRMAT4),
     *                          WORK(KRMAT3),WORK(KFCKAK),WORK(KSMAT),
     *                          WORK(KTMAT),ISYMIM,WORK(KXIAJB),ISINT1,
     *                          WORK(KINDSQ),LENSQ,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 after CC3_ONEL3',
     *                               RHO1N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
                     WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_ONEL3',
     *                               RHO2N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
                     WRITE(LUPRI,*) 'Norm of Rho2- after CC3_ONEL3',
     *                               RHO2N
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 after ONEL',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 after ONEL',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT3),1,
     *                       WORK(KRMAT3),1)
                     WRITE(LUPRI,*) 'Norm of RMAT3 after ONEL',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT4),1,
     *                       WORK(KRMAT4),1)
                     WRITE(LUPRI,*) 'Norm of RMAT4 after ONEL',XRMAT
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT after ONEL',XSMAT
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KTMAT),1)
                     WRITE(LUPRI,*) 'Norm of TMAT after ONEL',XTMAT
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_ONEL3: Rho1 & Rho2+ ')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,1,1)
                     CALL AROUND('After CC3_ONEL3: Rho2- ')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TIOME1 = TIOME1   + DTIME
C
C-------------------------------------------------------------------
C                 Accumulate the R2 and R3 matrix in Omega2+ & Omega2-.
C-------------------------------------------------------------------
C
                  MFACTOR = ONE
                  IOPT    = 1
                  CALL CC3_RACC3(OMEGA2P,OMEGA2M,WORK(KRMAT2),ISYMB,B,
     *                           ISYRES,MFACTOR,IOPT)
C
                  MFACTOR = XMONE
                  IOPT    = 1
                  CALL CC3_RACC3(OMEGA2P,OMEGA2M,WORK(KRMAT3),ISYMB,B,
     *                           ISYRES,MFACTOR,IOPT)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 after CC3_RACC3',
     *                               RHO1N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
                     WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_RACC3',
     *                               RHO2N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
                     WRITE(LUPRI,*) 'Norm of Rho2- after CC3_RACC3',
     *                               RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_RACC3: Rho1 & Rho2+')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,1,1)
                     CALL AROUND('After CC3_RACC3: Rho2-')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
                  ENDIF
C
  130          CONTINUE
  120       CONTINUE
C
C---------------------------------------------------------------------
C           Accumulate the R1 & R4 matrices in Omega2+ & Omega2-.
C---------------------------------------------------------------------
C
            MFACTOR = ONE
            IOPT    = 1
            CALL CC3_RACC3(OMEGA2P,OMEGA2M,WORK(KRMAT1),ISYMD,D,ISYRES,
     *                     MFACTOR,IOPT)
            MFACTOR = XMONE
            IOPT    = 1
            CALL CC3_RACC3(OMEGA2P,OMEGA2M,WORK(KRMAT4),ISYMD,D,ISYRES,
     *                     MFACTOR,IOPT)
C
            IF (IPRINT .GT. 55) THEN
               RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
               WRITE(LUPRI,*) 'Norm of Rho1 after CC3_RACC3-2',RHO1N
               RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
               WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_RACC3-2',RHO2N
               RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
               WRITE(LUPRI,*) 'Norm of Rho2- after CC3_RACC3-2',RHO2N
            ENDIF
C
            IF (IPRINT .GT. 220) THEN
               CALL AROUND('After CC3_RACC3-2: Rho1 & Rho2+')
               CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,1,1)
               CALL AROUND('After CC3_RACC3-2: Rho2-')
               CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
            ENDIF
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------
C     Close files
C-----------------------
C
      CALL WCLOSE2(LUDELD,FNDELD,'DELETE')
      CALL WCLOSE2(LUDKBC,FNDKBC,'DELETE')
C
C-------------------
C     Print timings.
C-------------------
C
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) '** Timings in CC3_OMEG3 **'
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
      CALL QEXIT('CC3_OMEG3')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck cc3_onel3 */
      SUBROUTINE CC3_ONEL3(OMEGA1,OMEGA2P,OMEGA2M,RMAT1,RMAT2,FOCKAK,
     *                     SMAT,TMAT,ISYMIM,XIAJB,ISYINT,INDSQ,LENSQ,
     *                     WORK,LWORK,ISYMIB,IB,ISYMID,ID)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Calculate Omega1 and Fock contibution to Omega2.
C
C     Ove Christiansen 9-1-1996:
C
C     General symmetry: ISYMIM is symmetry of SMAT and TMAT 
C                       intermdiates.(incl isymd,isymb)
C                       ISYINT is symmetry of FOCKAK and XIAJB
C                       ISYRES = ISYMIM*ISYINT
C
C     K. Hald, Jan. 2001 : Adapted to triplet.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMIM, ISYINT, LENSQ, LWORK, ISYMIB, IB, ISYMID, ID
      INTEGER ISYRES, ISYMB, ISYMC, ISYMK, ISYMBC, ISYAIJ
      INTEGER ISYMCK, LENGTH, NCK, NTOAIJ, NTOTC, ISYMI
      INTEGER ISYAKJ, ISYMJ, ISYMBJ, ISYMAK, NBJ, NAK, NAKBJ, NAKJ
      INTEGER NTOAKJ, NBK, ISYMBK, NTOTB, ISYMKJ, ISYMAI, NCI, NIJ
      INTEGER NCIBJ, NTOTIJ, NTOTAK, ISYMIJ, JSAIKJ, KOFF1, KOFF2
      INTEGER NKJ, NCKBJ, NTOTAI, JSAKIJ, ISYMCI
      INTEGER INDEX, INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2P(*), OMEGA2M(*), RMAT1(*)
      REAL RMAT2(*), FOCKAK(*),SMAT(*), TMAT(*), XIAJB(*)
      REAL WORK(LWORK), ZERO, ONE, TWO, HALF
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2P(*), OMEGA2M(*), RMAT1(*)
      DOUBLE PRECISION RMAT2(*), FOCKAK(*),SMAT(*), TMAT(*), XIAJB(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE, TWO, HALF
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_ONEL3')
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
         CALL QUIT('Not enough core in CCSDT_ONEL')
      END IF
C
c     CALL CC_GATHER(LENGTH,TMAT,SMAT,INDSQ(1,4))
c     CALL CC_GATHER(LENGTH,WORK,SMAT,INDSQ(1,3))
c     CALL DAXPY(LENGTH,-ONE,WORK,1,TMAT,1)
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(INDSQ(I,4)) 
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
            CALL QUIT('Not enough core in CCSDT_ONEL')
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
         CALL QUIT('Not enough core in CCSDT_ONEL')
      END IF
C
c     CALL CC_GATHER(LENGTH,TMAT,SMAT,INDSQ(1,5))
c     CALL DAXPY(LENGTH,-ONE,SMAT,1,TMAT,1)
C
      DO I = 1,LENGTH
         TMAT(I) =  SMAT(INDSQ(I,5))
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
         CALL QUIT('Not enough core in CCSDT_ONEL')
      END IF
C
c     CALL CC_GATHER(LENGTH,TMAT,SMAT,INDSQ(1,1))
c     CALL DAXPY(LENGTH,-ONE,SMAT,1,TMAT,1)
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(INDSQ(I,1)) 
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
                  NCIBJ = IT2AM(ISYMCI,ISYMBJ) + INDEX(NCI,NBJ)
C
                  IF (NCI .EQ. NBJ) THEN
C
                     OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + TWO*WORK(NIJ)
                     OMEGA2M(NCIBJ) = ZERO
C
                  ELSE IF (NCI .GT. NBJ) THEN
C
                     OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
                     OMEGA2M(NCIBJ) = OMEGA2M(NCIBJ) - HALF*WORK(NIJ)
C
                  ELSE IF (NCI .LT. NBJ) THEN
C
                     OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
                     OMEGA2M(NCIBJ) = OMEGA2M(NCIBJ) + HALF*WORK(NIJ)
C
                  ENDIF
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
                  OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
                  OMEGA2M(NCIBJ) = OMEGA2M(NCIBJ) + HALF*WORK(NIJ)
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
                  OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
                  OMEGA2M(NCIBJ) = OMEGA2M(NCIBJ) - HALF*WORK(NIJ)
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
      CALL QEXIT('CC3_ONEL3')
C
      RETURN
      END
C  /* Deck cc3_racc3 */
      SUBROUTINE CC3_RACC3(OMEGA2P,OMEGA2M,RMAT,ISYMB,B,ISYRES,MFACTOR,
     *                     IOPT)
C
C     Written by HK and ASM Maj 1995.
C
C     Purpose :  Accumulate the R matrix into the Omega vector.
C
C     OC: 9-1-1996 general symmetry, ISYRES is symmetry of resulting vector.
C
C     K. Hald, Jan 2001 : Adapted to triplet.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMB, ISYRES, ISYAIJ, ISYMJ, ISYMAI, ISYMBJ
      INTEGER NBJ, NBJJ, NAI, KOFF1, KOFF2, KOFF5, KOFF6, IOPT
      INTEGER INDEX
C
#if defined (SYS_CRAY)
      REAL OMEGA2P(*), OMEGA2M(*), RMAT(*), MFACTOR
      REAL ZERO, HALF, TWO, FACT
#else
      DOUBLE PRECISION OMEGA2P(*), OMEGA2M(*), RMAT(*), MFACTOR
      DOUBLE PRECISION ZERO, HALF, TWO, FACT
#endif
      PARAMETER (TWO = 2.0D0, HALF = 0.5D0, ZERO = 0.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_RACC3')
C
C---------------------------
C     Sanity check.
C---------------------------
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2))
     *    CALL QUIT('Wrong IOPT in CC3_RACC3')
C
C---------------------------
C     Calculate.
C---------------------------
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
C
      FACT = HALF * MFACTOR
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
                  IF (NAI .EQ. NBJ) THEN
                     IF (IOPT .EQ. 1) THEN
                        OMEGA2P(KOFF1) = OMEGA2P(KOFF1) + RMAT(KOFF2)
                     ENDIF
                     OMEGA2M(KOFF1) = ZERO
                  ELSE IF (NAI .GT. NBJ) THEN
                     IF (IOPT .EQ. 1) THEN
                        OMEGA2P(KOFF1) = OMEGA2P(KOFF1) + RMAT(KOFF2)
                     ENDIF
                     OMEGA2M(KOFF1) = OMEGA2M(KOFF1) + FACT*RMAT(KOFF2)
                  ELSE
                     IF (IOPT .EQ. 1) THEN
                        OMEGA2P(KOFF1) = OMEGA2P(KOFF1) + RMAT(KOFF2)
                     ENDIF
                     OMEGA2M(KOFF1) = OMEGA2M(KOFF1) - FACT*RMAT(KOFF2)
                  ENDIF
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
C
                  IF (IOPT .EQ. 1) THEN
                     OMEGA2P(KOFF5) = OMEGA2P(KOFF5) + RMAT(KOFF6)
                  ENDIF
                  OMEGA2M(KOFF5) = OMEGA2M(KOFF5) - FACT*RMAT(KOFF6)
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
C
                  IF (IOPT .EQ. 1) THEN
                     OMEGA2P(KOFF5) = OMEGA2P(KOFF5) + RMAT(KOFF6)
                  ENDIF
                  OMEGA2M(KOFF5) = OMEGA2M(KOFF5) + FACT*RMAT(KOFF6)
C
  250          CONTINUE
C
            ENDIF
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_RACC3')
C
      RETURN
      END
C  /* Deck cc3_conocc3 */
      SUBROUTINE CC3_CONOCC3(OMEGA2P,OMEGA2M,RMAT1,RMAT2,SMAT,TMAT,
     *                       ISYMIM,TROCC,TROCC1,ISYINT,WORK,LWORK,
     *                       INDSQ,LENSQ,ISYMIB,IB,ISYMID,ID)
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
C     K. Hald, Jan 2001, Adapted to triplet.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMIM, ISYINT, LWORK, LENSQ, ISYMIB, IB, ISYMID, ID
      INTEGER ISYRES, ISYMB, ISYMJ, ISYMBJ, ISYMC, ISYMBC, ISYCKL
      INTEGER ISYMA, ISYMI, ISYMAI, ISYMAB, ISYMKL, ISYKLJ, JSAIKL
      INTEGER KOFF1, KOFF2, KOFF3, KOFF5, KOFF6, LENGTH, NTOTAI, NTOTKL
      INTEGER NAI, NBJ, NRHFI, NTOCKL, JSCKLI
      INTEGER INDEX, INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL OMEGA2P(*), OMEGA2M(*), RMAT1(*), RMAT2(*)
      REAL SMAT(*),TMAT(*), TROCC(*),TROCC1(*),WORK(LWORK)
      REAL XRMAT, XTMAT, DDOT, ZERO, ONE, TWO, HALF
#else
      DOUBLE PRECISION OMEGA2P(*), OMEGA2M(*), RMAT1(*), RMAT2(*)
      DOUBLE PRECISION SMAT(*),TMAT(*), TROCC(*),TROCC1(*),WORK(LWORK)
      DOUBLE PRECISION XRMAT, XTMAT, DDOT, ZERO, ONE, TWO, HALF
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CONOCC3')
C
      IF (LWORK .LT. LENSQ) THEN
         CALL QUIT('Insufficient core in CONOCC3')
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
         WRITE(LUPRI,*) 'In CC3_CONOCC3: 1. Norm of TMAT = ',XTMAT
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
         WRITE(LUPRI,*) 'In CC3_CONOCC3: Norm of RMAT2 =  ',XRMAT
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
         WRITE(LUPRI,*) 'In CC3_CONOCC3: 2. Norm of TMAT = ',XTMAT
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
         WRITE(LUPRI,*) 'In CC3_CONOCC3: Norm of RMAT1 =  ',XRMAT
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
         WRITE(LUPRI,*) 'In CC3_CONOCC3: 3. Norm of TMAT = ',XTMAT
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
            CALL QUIT('Insufficient memory in CC3_CONOCC3')
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
                  IF (NAI .EQ. NBJ) THEN
                     OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - TWO*WORK(KOFF6)
                     OMEGA2M(KOFF5) = ZERO
                  ELSE IF (NAI .LT. NBJ) THEN
                     OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - WORK(KOFF6)
                     OMEGA2M(KOFF5) = OMEGA2M(KOFF5) + HALF*WORK(KOFF6)
                  ELSE IF (NAI .GT. NBJ) THEN
                     OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - WORK(KOFF6)
                     OMEGA2M(KOFF5) = OMEGA2M(KOFF5) - HALF*WORK(KOFF6)
                  ENDIF
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
                  OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - WORK(KOFF6)
                  OMEGA2M(KOFF5) = OMEGA2M(KOFF5) + HALF*WORK(KOFF6)
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
                  OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - WORK(KOFF6)
                  OMEGA2M(KOFF5) = OMEGA2M(KOFF5) - HALF*WORK(KOFF6)
C
  640          CONTINUE
C
            ENDIF
C
  610    CONTINUE
C
  600 CONTINUE
C
      CALL QEXIT('CC3_CONOCC3')
C
      RETURN
      END
C  /* Deck cc3_convir3 */
      SUBROUTINE CC3_CONVIR3(RMAT,SMAT,QMAT,TMAT,ISYMIM,TRVIR,
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
C     K. Hald, Jan 2001, Adapted to triplet.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMIM, ISYINT, LWORK, LENSQ, ISYMB, ISYMD
      INTEGER KOFF1, KOFF2, KOFF3, LENGTH
      INTEGER KSCR1, KEND1, LWRK1
      INTEGER ISYRES, ISYMBD, ISCKIJ, ISYMJ, ISYMBJ, ISYMAI
      INTEGER ISYCKI, ISYMI, ISYMA, ISYMCK
      INTEGER NTOTCK, NVIRA
      INTEGER INDEX, INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL RMAT(*),SMAT(*),QMAT(*), TMAT(*),TRVIR(*)
      REAL TRVIR1(*), WORK(LWORK), ZERO, ONE, TWO
#else
      DOUBLE PRECISION RMAT(*),SMAT(*),QMAT(*), TMAT(*),TRVIR(*)
      DOUBLE PRECISION TRVIR1(*), WORK(LWORK), ZERO, ONE, TWO
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CONVIR3')
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
         CALL QUIT('Insufficient core in CCSDT_CONVIR')
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
            CALL QUIT('Insufficient work space in CCSDT_CONVIR')
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
            CALL QUIT('Insufficient work space in CCSDT_CONVIR')
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
      CALL QEXIT('CC3_CONVIR3')
C
      RETURN
      END
C  /* Deck cc3_omeg33 */
      SUBROUTINE CC3_OMEG33(ECURR,OMEGA1,OMEGA2P,OMEGA2M,
     *                      T1AM,ISYMT1,T2TP,
     *                      ISYMT2,C2AMP,C2AMM,ISYMC2,ISYMT12,
     *                      FOCK,XLAMDP,XLAMDH,WORK,LWORK,LU3SRT,FN3SRT,
     *                      LU3VI,FN3VI,LU3VI2,FN3VI2,LUCKJD,FNCKJD,
     *                      LUCKJD2,FNCKJD2,LUCKJD3,FNCKJD3,
     *                      LU3SRT2,FN3SRT2,LU3SRT3,FN3SRT3,LUTOC,FNTOC)
C
C     K. Hald, Spring 2001.
C
C     Based on CC3_OMEG
C
C     Add the triples contribution to the result vectors
C
C
C         ISYMT2 is symmetry of T2TP
C         ISYMT1 is symmetry of T1AM input when CC3LR
C         else ISYMT1 is symmetry of the C1 that the
C         triples integrals have been transformed with.
C         ISYMC2 is symmetry of C2+ and C2-.
C
C     FN3VI3 and FN3VI4 are used to store various integrals,
C     but contains name correct integrals after virtint
C
      IMPLICIT NONE
C
#include <iratdef.h>
#include <dummy.h>
#include <inftap.h>
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <second.h>
C
      INTEGER ISYMT1, ISYMT2, ISYMC2, ISYMT12, LWORK, ISYMTR, ISYMTR2
      INTEGER KEND0, LWRK0, KEND1, LWRK1, KEND2, LWRK2, KEND3, LWRK3
      INTEGER KEND4, LWRK4, BSTART, BEND, IOPT
      INTEGER KOFF1, KOFF2, ISYMIM, ISYRES, ISYOPE, LENSQ, LENGTH
      INTEGER KSMAT, KSMAT2, KQMAT, KDIAG, KINDSQ, KINDEX, KINDEX2
      INTEGER KINDEX3, KINDEX4, KTMAT, KRMAT1, KRMAT2, KRMAT3, KRMAT4
      INTEGER KCMO, KFOCKD, ISINT1, ISINT2, ISINT22, ISYMB, ISYMD
      INTEGER ISYMBD, ISYALJ, ISYMC, ISYMK, IOFF, ISAIJ1, ISAIJ2
      INTEGER ISYCKB, ISCKB1, ISCKIJ, ISCKB2, KTRVI, KTRVI0, KTRVI1
      INTEGER KTRVI2, KTRVI3, KTRVI4, KTRVI5, KTRVI6, KTRVI7
      INTEGER KTRVIB0, KTRVIB2, KTRVIB4, KTRVIB5, KTRVIB6, KTRVIB7
      INTEGER KINTVI, KINTOC, KTROC, KTROC0, KTROC1, KTROC2, KTROC3
      INTEGER KXIAJB, KFCKAK, ISYALJ2, ISYALJ3, ISYALJ4
      INTEGER KTROC4, KTROC5, ISYCKD, ISCKD1, ISCKD2, ISCKD3, ISCKB3
      INTEGER KTRVI8, KTRVI9, KXIAJB2
      INTEGER KRES2P, KRES2M
      INTEGER KR3
      INTEGER LUDKBC2, LUDKBC3, LUDKBC4, LUDKBC5, LU3VI3, LU3VI4
      INTEGER LUDELD, LUDKBC, LU3SRT, LU3VI, LU3VI2, LUCKJD
      INTEGER LUCKJD2, LUCKJD3, LU3SRT2, LU3SRT3, LUTOC
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), OMEGA2P(*), OMEGA2M(*), T1AM(*)
      REAL T2TP(*), C2AMP(*), C2AMM(*),FOCK(*),XLAMDP(*)
      REAL XLAMDH(*),WORK(LWORK), RHO1N, RHO2N
      REAL XIAJB, XINT, XTROC, XTROC0, XTROC1, XT2TP
      REAL XTRVI, XTRVI0, XTRVI1, XTRVI2, XTRVI3
      REAL XSMAT, XTMAT, XQMAT, XRMAT, XDIA, XFD
      REAL TITRAN, TISORT, TISMAT, TICONV, TICONO, TIIO, TIVINT, TIOME1
      REAL TIME1P, TIME2P,TIME3P, TIME1M, TIME2M, TSORTP, TSORTM
      REAL XTIME, DTIME, DDOT, XMONE, FACTOR, ONE, ECURR
#else
      DOUBLE PRECISION OMEGA1(*), OMEGA2P(*), OMEGA2M(*), T1AM(*)
      DOUBLE PRECISION T2TP(*), C2AMP(*), C2AMM(*),FOCK(*),XLAMDP(*)
      DOUBLE PRECISION XLAMDH(*),WORK(LWORK), RHO1N, RHO2N
      DOUBLE PRECISION XIAJB, XINT, XTROC, XTROC0, XTROC1, XT2TP
      DOUBLE PRECISION XTRVI, XTRVI0, XTRVI1, XTRVI2, XTRVI3
      DOUBLE PRECISION XSMAT, XTMAT, XQMAT, XRMAT, XDIA, XFD
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TICONV, TICONO, TIIO
      DOUBLE PRECISION TIVINT, TIOME1
      DOUBLE PRECISION TIME1P, TIME2P, TIME3P, TIME1M, TIME2M
      DOUBLE PRECISION TSORTP, TSORTM
      DOUBLE PRECISION XTIME, DTIME, DDOT, XMONE, FACTOR, ONE, ECURR
#endif
C
      PARAMETER (XMONE = -1.0D0, ONE = 1.0D0)
C
      LOGICAL LDEBUG
C
      CHARACTER*12 FNDKBC2, FNDKBC3, FNDKBC4, FNDKBC5, FN3VI3, FN3VI4
      CHARACTER*12 FNDELD, FNDKBC
      CHARACTER*1 CDUMMY
      CHARACTER*(*) FN3SRT, FN3VI, FN3VI2, FNCKJD, FNCKJD2, FNCKJD3
      CHARACTER*(*) FN3SRT2, FN3SRT3, FNTOC
C
      CALL QENTER('CC3_OMEG33')
C
      LDEBUG = .FALSE.
C
      CDUMMY = ' '
C
C-------------------------------------------------------------
C     Set symmetry flags.
C
C     omega = int1*T2*int2
C     isymres is symmetry of result(omega)
C     isint1 is symmetry of integrals in contraction.
C     isint2 is symmetry of integrals in the triples equation
C     used with T2.
C     isintt12 is symmetry of integrals in contraction
C     isint22 is symmetry of integrals in the triples equation
C     used with C2AMP and C2AMM.
C     isymim is symmetry of S and Q intermediates.(t2*int2)
C      (sym is for all index of S and Q (cbd,klj)
C       thus cklj=b*d*isymim)
C-------------------------------------------------------------
C
      IPRCC   = IPRINT
      ISYMTR  = MULD2H(ISYMT1,ISYMT2)
      ISYMTR2 = MULD2H(ISYMT12,ISYMC2)
      IF (ISYMTR .NE. ISYMTR2) THEN
          CALL QUIT('Symmetry error in CC3_OMEG33')
      ENDIF
      ISYRES  = MULD2H(ISYMTR,ISYMOP)
      ISINT1  = ISYMOP
      ISINT2  = MULD2H(ISYMT1,ISYMOP)
      ISINT22 = MULD2H(ISYMT12,ISYMOP)
      ISYMIM  = MULD2H(ISYMTR,ISYMOP)
      IF (CC3LR) THEN
         ISINT1 = MULD2H(ISYMT1,ISYMOP)
         ISINT2 = ISYMOP
         ISYMIM = ISYMOP
      ENDIF
C
      IF (IPRINT .GT. 20 .OR. LDEBUG) THEN
        WRITE(LUPRI,*) ' In CC3_OMEG33: CC1A  = ',CC1A
        WRITE(LUPRI,*) ' In CC3_OMEG33: CC1B  = ',CC1B
        WRITE(LUPRI,*) ' In CC3_OMEG33: CC3   = ',CC3
        WRITE(LUPRI,*) ' In CC3_OMEG33: CC3LR = ',CC3LR
        WRITE(LUPRI,*) ' In CC3_OMEG33: ISYMT1 , ISYMT2:',ISYMT1,ISYMT2
        WRITE(LUPRI,*) ' In CC3_OMEG33: ISYMT12, ISYMC2:',ISYMT12,ISYMT2
        WRITE(LUPRI,*) ' In CC3_OMEG33: ISYRES , ISYMOP:',ISYRES,ISYMOP
        WRITE(LUPRI,*) ' In CC3_OMEG33: ISINT1 , ISINT2:',ISINT1,ISINT2
        WRITE(LUPRI,*) ' In CC3_OMEG33: ISINT22 ',ISINT22
      ENDIF
C
C--------------------
C     Time variables.
C--------------------
C
      TITRAN = 0.0D0
      TISORT = 0.0D0
      TISMAT = 0.0D0
      TICONO = 0.0D0
      TIVINT = 0.0D0
      TIIO   = 0.0D0
      TICONV = 0.0D0
      TIME1P = 0.0D0
      TIME2P = 0.0D0
      TIME3P = 0.0D0
      TIME1M = 0.0D0
      TIME2M = 0.0D0
      TSORTP = 0.0D0
      TSORTM = 0.0D0
      TIOME1 = 0.0D0
C
C-------------------------------
C     Open temp. files.
C-------------------------------
C
      LUDKBC2 = -1
      LUDKBC3 = -1
      LUDKBC4 = -1
      LUDKBC5 = -1
      LU3VI3  = -1
      LU3VI4  = -1
      LUDELD  = -1
      LUDKBC  = -1
C
      FNDKBC2 = 'CC3_OMEG33_1'
      FNDKBC3 = 'CC3_OMEG33_2'
      FNDKBC4 = 'CC3_OMEG33_3'
      FNDKBC5 = 'CC3_OMEG33_4'
      FN3VI3  = 'CC3_OMEG33_5'
      FN3VI4  = 'CC3_OMEG33_6'
      FNDELD  = 'CC3_OMEG33_7'
      FNDKBC  = 'CC3_OMEG33_8'
C
      CALL WOPEN2(LUDKBC2,FNDKBC2,64,0)
      CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
      CALL WOPEN2(LUDKBC4,FNDKBC4,64,0)
      CALL WOPEN2(LUDKBC5,FNDKBC5,64,0)
      CALL WOPEN2(LU3VI3,FN3VI3,64,0)
      CALL WOPEN2(LU3VI4,FN3VI4,64,0)
      CALL WOPEN2(LUDELD,FNDELD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
C
C--------------------------------------------------
C     Allocate and initialise two result vectors
C--------------------------------------------------
C
      KRES2P  = 1
      KRES2M  = KRES2P + NT2SQ(ISYRES)
      KEND0   = KRES2M + NT2SQ(ISYRES)
C
      CALL DZERO(WORK(KRES2P),NT2SQ(ISYRES))
      CALL DZERO(WORK(KRES2M),NT2SQ(ISYRES))
C
C---------------------------------------------------------
C     Transform and sort qmat integrals to smat integrals.
C
C     CHECK IF THE IOPT = 2 (SORT) and IOPT = 1 (SINT)
C     IS NECESSARY (MAYBE CALCULATED IN T3 PART)
C
C---------------------------------------------------------
C
      CALL CC3_SORT1(WORK,LWORK,2,ISINT22,LU3SRT,FN3SRT,
     *               LUDELD,FNDELD,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
      CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT22,LUDELD,FNDELD,
     *              LUDKBC,FNDKBC)
      CALL CC3_SORT3(WORK,LWORK,1,ISINT2,LU3SRT2,FN3SRT2,
     *               LU3SRT3,FN3SRT3,LU3VI4,FN3VI4,LU3VI3,FN3VI3)
      CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT2,LU3VI4,FN3VI4,
     *              LUDKBC2,FNDKBC2)
      CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT2,LU3VI3,FN3VI3,
     *              LUDKBC3,FNDKBC3)
      CALL CC3_SORT3(WORK,LWORK,2,ISINT2,LU3SRT2,FN3SRT2,
     *               LU3SRT3,FN3SRT3,LU3VI4,FN3VI4,IDUMMY,CDUMMY)
      CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT2,LU3VI4,FN3VI4,
     *              LUDKBC4,FNDKBC4)
      CALL CC3_SORT2(WORK,LWORK,ISINT22,LU3SRT,FN3SRT,LU3VI4,FN3VI4)
      CALL CC3_SINT(XLAMDH,WORK,LWORK,ISINT22,LU3VI4,FN3VI4,
     *              LUDKBC5,FNDKBC5)
C
C---------------------------------------------
C     Reorder the t2 amplitudes in T2TP.
C---------------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP in CC3_OMEG33')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMT2),T2TP,1,WORK,1)
      CALL CC3_T2TP(T2TP,WORK,ISYMT2)
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XT2TP = DDOT(NT2SQ(ISYMT2),T2TP,1,T2TP,1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XT2TP
      ENDIF
C
C-----------------------------------------------------
C     Reorder the R2+ and R2- amplitudes in C2AMP.
C-----------------------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMC2)) THEN
         CALL QUIT('Not enough memory to construct C2AMP in CC3_OMEG33')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMC2),C2AMP,1,WORK,1)
      CALL CC3_T2TP(C2AMP,WORK,ISYMC2)
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XT2TP = DDOT(NT2SQ(ISYMC2),C2AMP,1,C2AMP,1)
         WRITE(LUPRI,*) 'Norm of C2AMM ',XT2TP
      ENDIF
C
C
      IF (LWORK .LT. NT2SQ(ISYMC2)) THEN
         CALL QUIT('Not enough memory to construct C2AMM in CC3_OMEG33')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMC2),C2AMM,1,WORK,1)
      CALL CC3_T2TP(C2AMM,WORK,ISYMC2)
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XT2TP = DDOT(NT2SQ(ISYMC2),C2AMM,1,C2AMM,1)
         WRITE(LUPRI,*) 'Norm of C2AMM ',XT2TP
      ENDIF
C
C---------------------------------------------------------
C     Read canonical orbital energies and MO coefficients.
C---------------------------------------------------------
C
      KFOCKD = KEND0
      KCMO   = KFOCKD + NORBTS
      KFCKAK = KCMO   + NLAMDS
      KEND0  = KFCKAK + NT1AM(ISINT1)
      LWRK0  = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_OMEG33')
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
     *               + NVIR(ISYMC)*(K - 1) + C - 1
C
               WORK(KOFF2) = FOCK(KOFF1)
C
   70       CONTINUE
   60    CONTINUE
   50 CONTINUE
C
C---------------------------------------------------------
C     Transform the "virtual" integrals and store them
C     on file.
C---------------------------------------------------------
C
      DTIME = SECOND()
      CALL CC3_VIRTINT(XLAMDP,WORK(KEND0),LWRK0,ISINT1,ISINT2,
     *                 LU3VI,FN3VI,LU3VI2,FN3VI2,LU3VI4,FN3VI4,
     *                 LU3VI3,FN3VI3)
      DTIME  = SECOND() - DTIME
      TIVINT = TIVINT   + DTIME
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC   = KEND0
      KTROC1  = KTROC  + NTRAOC(ISINT1)
      KTROC0  = KTROC1 + NTRAOC(ISINT1)
      KTROC2  = KTROC0 + NTRAOC(ISINT22)
      KTROC3  = KTROC2 + NTRAOC(ISINT2)
      KXIAJB  = KTROC3 + NTRAOC(ISINT2)
      KXIAJB2 = KXIAJB + NT2AM(ISYMOP)
      KEND1   = KXIAJB2 + NT2AM(ISYMOP)
C
      IF (LDEBUG) THEN
         KR3    = KEND1
         KEND1  = KR3 + NVIRT*NVIRT*NVIRT*NRHFT*NRHFT*NRHFT
         CALL DZERO(WORK(KR3),NVIRT*NVIRT*NVIRT*NRHFT*NRHFT*NRHFT)
      ENDIF
C
      KTROC4 = KEND1
      KTROC5 = KTROC4 + MAX(NTRAOC(ISINT2),NTRAOC(ISINT22))
      KEND1  = KTROC5 + NTRAOC(ISINT2)
      LWRK1  = LWORK  - KEND1
C
      KINTOC = KEND1
      KEND2  = KINTOC + MAX(NTOTOC(ISYMOP),NTOTOC(ISINT2))
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_OMEG33')
      END IF
C
C---------------------------------
C     Construct -Exchange(ia,jb)
C---------------------------------
C
      LENGTH = IRAT*NT2AM(ISYMOP)
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
C
      ISYOPE = ISYMOP
      IOPT = 2
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYOPE,IOPT)
C
C------------------------
C     Occupied integrals.
C------------------------
C
      DTIME = SECOND()
      IOFF = 1
      IF (NTOTOC(ISYMOP) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
      ENDIF
      DTIME  = SECOND() - DTIME
      TIIO = TIIO   + DTIME
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XINT  = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CC3_OMEG33_OC-INT ',XINT
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      DTIME = SECOND()
         CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC),XLAMDH,
     *                    WORK(KEND2),LWRK2)
      DTIME  = SECOND() - DTIME
      TITRAN = TITRAN   + DTIME
C
C------------------------------------
C     Read in integrals part 2.
C------------------------------------
C
      DTIME = SECOND()
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUCKJD2,FNCKJD2,WORK(KINTOC),IOFF,
     *               NTOTOC(ISINT2))
      ENDIF
      DTIME  = SECOND() - DTIME
      TIIO = TIIO   + DTIME
C
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ik,j,a)
C----------------------------------------------------------------------
C
      DTIME = SECOND()
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC2),WORK(KCMO),
     *                 WORK(KEND2),LWRK2,ISINT2)
      DTIME  = SECOND() - DTIME
      TIIO = TIIO   + DTIME
C
C------------------------------------
C     Read in integrals part 3.
C------------------------------------
C
      DTIME = SECOND()
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUCKJD3,FNCKJD3,WORK(KINTOC),IOFF,
     *               NTOTOC(ISINT2))
      ENDIF
      DTIME  = SECOND() - DTIME
      TIIO = TIIO   + DTIME
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ik,j,a)
C----------------------------------------------------------------------
C
      DTIME = SECOND()
      CALL CC3_TROCC3(WORK(KINTOC),WORK(KTROC3),WORK(KTROC4),
     *                WORK(KCMO),WORK(KEND2),LWRK2,ISINT2)
C
      DTIME  = SECOND() - DTIME
      TITRAN = TITRAN   + DTIME
C
C---------------------------------------------------
C     Construct the special occupied integrals.
C     (g^3 and g^4)
C---------------------------------------------------
C
      DTIME = SECOND()
      CALL DAXPY(NTRAOC(ISINT2),XMONE,WORK(KTROC2),1,WORK(KTROC3),1)
      CALL DAXPY(NTRAOC(ISINT2),XMONE,WORK(KTROC4),1,WORK(KTROC2),1)
C
      CALL CCSDT_SRTOC3(WORK(KTROC3),WORK(KTROC5),ISINT2,
     *                  WORK(KEND2),LWRK2)
C
      DTIME  = SECOND() - DTIME
      TISORT = TISORT   + DTIME
C
C-----------------------
C     Read in integrals.
C-----------------------
C
      DTIME = SECOND()
      IOFF = 1
      IF (NTOTOC(ISINT22) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT22))
      ENDIF
      DTIME  = SECOND() - DTIME
      TIIO = TIIO   + DTIME
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ik,j,a)
C----------------------------------------------------------------------
C
      DTIME = SECOND()
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KCMO),
     *                 WORK(KEND2),LWRK2,ISINT22)
C
      DTIME  = SECOND() - DTIME
      TITRAN = TITRAN   + DTIME

C
      DTIME = SECOND()
C
      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1,
     *                  WORK(KEND2),LWRK2)
C
      CALL DZERO(WORK(KTROC4),MAX(NTRAOC(ISINT2),NTRAOC(ISINT22)))
      CALL CCSDT_SRTOC3(WORK(KTROC0),WORK(KTROC4),ISINT22,
     *                  WORK(KEND2),LWRK2)
C
      DTIME  = SECOND() - DTIME
      TISORT = TISORT   + DTIME
C
C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
C
         XTROC = DDOT(NTRAOC(ISINT1),WORK(KTROC),1,
     *                WORK(KTROC),1)
         WRITE(LUPRI,*) 'Norm of TROC ',XTROC
C
         XTROC0 = DDOT(NTRAOC(ISINT22),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 ',XTROC0
C
         XTROC1 = DDOT(NTRAOC(ISINT1),WORK(KTROC1),1,
     *                WORK(KTROC1),1)
         WRITE(LUPRI,*) 'Norm of TROC1 ',XTROC1
C
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC2),1,
     *                WORK(KTROC2),1)
         WRITE(LUPRI,*) 'Norm of TROC2 ',XTROC0
C
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC3),1,
     *                WORK(KTROC3),1)
         WRITE(LUPRI,*) 'Norm of TROC3 ',XTROC0
C
         XTROC0 = DDOT(NTRAOC(ISINT22),WORK(KTROC4),1,
     *                WORK(KTROC4),1)
         WRITE(LUPRI,*) 'Norm of TROC4 ',XTROC0
C
         XINT  = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XINT
C
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
         ISCKB3 = MULD2H(ISINT22,ISYMD)
         ISYALJ3= MULD2H(ISYMD,ISYMT2)
         ISYALJ4= MULD2H(ISYMD,ISYMC2)
C
         IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
C
            WRITE(LUPRI,*) 'In CC3_OMEG33: ISAIJ1:',ISAIJ1
            WRITE(LUPRI,*) 'In CC3_OMEG33: ISYCKB:',ISYCKB
            WRITE(LUPRI,*) 'In CC3_OMEG33: ISCKB1:',ISCKB1
            WRITE(LUPRI,*) 'In CC3_OMEG33: ISCKB2:',ISCKB2
            WRITE(LUPRI,*) 'In CC3_OMEG33: ISCKB3:',ISCKB3
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
         KRMAT1 = KTRVI2 + NCKATR(ISCKB3)
         KRMAT3 = KRMAT1 + NCKI(ISAIJ1)
         KEND2  = KRMAT3 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0 = KEND2
         KTRVI3 = KTRVI0 + NCKATR(ISCKB3)
         KTRVI4 = KTRVI3 + NCKATR(ISCKB2)
         KTRVI5 = KTRVI4 + NCKATR(ISCKB2)
         KTRVI6 = KTRVI5 + NCKATR(ISCKB2)
         KTRVI7 = KTRVI6 + MAX(NCKATR(ISCKB2),NCKATR(ISCKB3))
         KINDEX3= KTRVI7 + NCKATR(ISCKB3)
         KINDEX4= KINDEX3 + (NCKI(ISYALJ3) - 1)/IRAT + 1
         KEND3  = KINDEX4 + (NCKI(ISYALJ4) - 1)/IRAT + 1
         LWRK3  = LWORK  - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISYMD),NCKA(ISCKB2))
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_OMEG33')
         END IF
C
C------------------------------------------
C        Construct the first index arrays.
C------------------------------------------
C
         CALL CC3_INDEX(WORK(KINDEX3),ISYALJ3)
         CALL CC3_INDEX(WORK(KINDEX4),ISYALJ4)
C
C-----------------------------------------------------------
C        Initialise the result vectors for cc3_convir33.
C-----------------------------------------------------------
C
         CALL DZERO(WORK(KRES2P),NT2SQ(ISYRES))
         CALL DZERO(WORK(KRES2M),NT2SQ(ISYRES))
C
C--------------------------
C        Sum over D.
C--------------------------
C
         DO 110 D = 1,NVIR(ISYMD)
C
            CALL DZERO(WORK(KRMAT1),NCKI(ISAIJ1))
            CALL DZERO(WORK(KRMAT3),NCKI(ISAIJ1))
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction.
C           The integrals are stored transformed!
C------------------------------------------------------------
C
            DTIME = SECOND()
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI3,FN3VI3,WORK(KTRVI),IOFF,
     &                     NCKATR(ISYCKB))
            ENDIF
C
C-----------------------------------------------------------
C           Read virtual integrals used in s3am part 1.
C-----------------------------------------------------------
C
            IOFF = ICKBD(ISCKB3,ISYMD) + NCKATR(ISCKB3)*(D - 1) + 1
            IF (NCKATR(ISCKB3) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                     NCKATR(ISCKB3))
C
               IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                  XTRVI0= DDOT(NCKATR(ISCKB3),WORK(KTRVI0),1,
     *                         WORK(KTRVI0),1)
                  WRITE(LUPRI,*) 'Norm of TRVI0 readin',XTRVI0
               ENDIF
            ENDIF
C
            DTIME  = SECOND() - DTIME
            TIIO = TIIO   + DTIME
C
C---------------------------------------
C           Sort the integrals for s3am.
C---------------------------------------
C
            DTIME = SECOND()
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT22)
C
            IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
               XTRVI0= DDOT(NCKATR(ISCKB3),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 ',XTRVI0
C
               XTRVI2= DDOT(NCKATR(ISCKB3),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XTRVI2
            ENDIF
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
C-----------------------------------------------------------
C           Read virtual integrals used in s3am part 2.
C-----------------------------------------------------------
C
            DTIME = SECOND()
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC2,FNDKBC2,WORK(KTRVI4),IOFF,
     &                     NCKATR(ISCKB2))
C
               IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                  XTRVI0= DDOT(NCKATR(ISCKB2),WORK(KTRVI4),1,
     *                         WORK(KTRVI4),1)
                  WRITE(LUPRI,*) 'Norm of TRVI4 readin',XTRVI0
               ENDIF
            ENDIF
            DTIME  = SECOND() - DTIME
            TIIO = TIIO   + DTIME
C
C
C-------------------------------------------------
C           Sort the integrals for s3am part 2.
C-------------------------------------------------
C
            DTIME = SECOND()
            CALL CCSDT_SRTVIR(WORK(KTRVI4),WORK(KTRVI5),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
               XTRVI2= DDOT(NCKATR(ISCKB2),WORK(KTRVI5),1,
     *                      WORK(KTRVI5),1)
               WRITE(LUPRI,*) 'Norm of TRVI5 ',XTRVI2
            ENDIF
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
C-----------------------------------------------------------
C           Read virtual integrals used in s3am part 2b.
C-----------------------------------------------------------
C
            DTIME = SECOND()
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC4,FNDKBC4,WORK(KTRVI6),IOFF,
     &                     NCKATR(ISCKB2))
            ENDIF
C
            DTIME  = SECOND() - DTIME
            TIIO = TIIO   + DTIME
C
C----------------------------------------------------
C           Sort the integrals for s3am part 2b.
C----------------------------------------------------
C
            DTIME = SECOND()
            CALL CCSDT_SRTVIR(WORK(KTRVI6),WORK(KTRVI4),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
               XTRVI0= DDOT(NCKATR(ISCKB2),WORK(KTRVI4),1,
     *                      WORK(KTRVI4),1)
               WRITE(LUPRI,*) 'Norm of TRVI4 ',XTRVI0
C
               XTRVI2= DDOT(NCKATR(ISCKB2),WORK(KTRVI5),1,
     *                      WORK(KTRVI5),1)
               WRITE(LUPRI,*) 'Norm of TRVI5 ',XTRVI2
            ENDIF
C
C-------------------------------------------------------------------
C           Read virtual integrals used in s3am part 3a.
C           g integrals. Need only the transformed.
C-------------------------------------------------------------------
C
            DTIME = SECOND()
C
            CALL DZERO(WORK(KTRVI6),MAX(NCKATR(ISCKB2),NCKATR(ISCKB3)))
            IOFF = ICKBD(ISCKB3,ISYMD) + NCKATR(ISCKB3)*(D - 1) + 1
            IF (NCKATR(ISCKB3) .GT. 0) THEN
               CALL GETWA2(LUDKBC5,FNDKBC5,WORK(KTRVI6),IOFF,
     &                     NCKATR(ISCKB3))
            ENDIF
C
            DTIME  = SECOND() - DTIME
            TIIO = TIIO   + DTIME
C
C---------------------------------------
C           Sort the integrals for s3am.
C---------------------------------------
C
            DTIME = SECOND()
            CALL CCSDT_SRTVIR(WORK(KTRVI6),WORK(KTRVI7),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT22)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
               XTRVI2= DDOT(NCKATR(ISCKB3),WORK(KTRVI7),1,
     *                      WORK(KTRVI7),1)
               WRITE(LUPRI,*) 'Norm of TRVI7 ',XTRVI2
            ENDIF
C
C-------------------------------------------------------------------
C           Read virtual integrals used in s3am part 3b.
C           This is the g^1 integral. Do not need the transformed.
C-------------------------------------------------------------------
C
            DTIME = SECOND()
C
            CALL DZERO(WORK(KTRVI6),MAX(NCKATR(ISCKB2),NCKATR(ISCKB3)))
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI6),IOFF,
     &                     NCKATR(ISCKB2))
            ENDIF
C
C-----------------------------------------------------------
C           Read virtual integrals used in contraction.
C-----------------------------------------------------------
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI4,FN3VI4,WORK(KTRVI1),IOFF,
     &                     NCKATR(ISYCKB))
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
            DTIME  = SECOND() - DTIME
            TIIO = TIIO   + DTIME
C
            DTIME = SECOND()
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),XLAMDH,
     *                       ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_OMEG33')
            END IF
C
            DTIME  = SECOND() - DTIME
            TITRAN = TITRAN   + DTIME
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND3),ISYMD,D,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
               XTRVI3= DDOT(NCKATR(ISCKB2),WORK(KTRVI3),1,
     *                      WORK(KTRVI3),1)
               WRITE(LUPRI,*) 'Norm of TRVI3 ',XTRVI3
            ENDIF
C
C---------------------
C           Calculate.
C---------------------
C
            DO 120 ISYMB = 1,NSYM
C
               ISYALJ  = MULD2H(ISYMB,ISYMT2)
               ISYALJ2 = MULD2H(ISYMB,ISYMC2)
               ISAIJ2  = MULD2H(ISYMB,ISYRES)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYMIM)
C
               ISYCKD = MULD2H(ISYMB,ISYMOP)
               ISCKD1 = MULD2H(ISINT1,ISYMB)
               ISCKD2 = MULD2H(ISINT2,ISYMB)
               ISCKD3 = MULD2H(ISINT22,ISYMB)
C
               IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
C
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISCKIJ:',ISCKIJ
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISYCKD:',ISYCKD
                  WRITE(LUPRI,*) 'In CC3_OMEG33: ISCKD1:',ISCKD1
C
               ENDIF
C
               KSMAT   = KEND3
               KSMAT2  = KSMAT   + NCKIJ(ISCKIJ)
               KQMAT   = KSMAT2  + NCKIJ(ISCKIJ)
               KDIAG   = KQMAT   + NCKIJ(ISCKIJ)
               KINDSQ  = KDIAG   + NCKIJ(ISCKIJ)
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KTMAT   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KRMAT2  = KTMAT   + NCKIJ(ISCKIJ)
               KRMAT4  = KRMAT2  + NCKI(ISAIJ2)
               KEND4   = KRMAT4  + NCKI(ISAIJ2)
               LWRK4   = LWORK   - KEND4
               KTRVI8  = KEND4
               KTRVI9  = KTRVI8  + NCKATR(ISCKD1)
               KTRVIB0 = KTRVI9  + NCKATR(ISCKD1)
               KTRVIB2 = KTRVIB0 + NCKATR(ISCKD3)
               KTRVIB4 = KTRVIB2 + NCKATR(ISCKD3)
               KTRVIB5 = KTRVIB4 + NCKATR(ISCKD2)
               KTRVIB6 = KTRVIB5 + NCKATR(ISCKD2)
               KTRVIB7 = KTRVIB6 + MAX(NCKATR(ISCKD2),NCKATR(ISCKD3))
               KEND4   = KTRVIB7 + NCKATR(ISCKD3)
               LWRK4   = LWORK   - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_OMEG33')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
C
               IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
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
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJ2)
C
C--------------------------------------
C              Sum over B
C--------------------------------------
C
               IF (ISYMB .GT. ISYMD) THEN
                    BSTART = 1
                    BEND   = NVIR(ISYMB)
               ELSE IF (ISYMB .EQ. ISYMD) THEN
                    BSTART = D + 1
                    BEND   = NVIR(ISYMB)
               ELSE
                     BSTART = 1
                     BEND   = 0
               ENDIF
C
C
               DO 130 B = BSTART,BEND
C
C----------------------------------------------------------
C                 Initialize the R2/R4/SMAT/TMAT matrices.
C----------------------------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
                  CALL DZERO(WORK(KRMAT4),NCKI(ISAIJ2))
                  CALL DZERO(WORK(KSMAT),NCKIJ(ISCKIJ))
                  CALL DZERO(WORK(KSMAT2),NCKIJ(ISCKIJ))
                  CALL DZERO(WORK(KTMAT),NCKIJ(ISCKIJ))
C
C---------------------------------------------------------
C                 Read the integrals for given B that is
C                 needed to calculate the SMAT.
C---------------------------------------------------------
C
                  DTIME = SECOND()
                  IOFF = ICKBD(ISCKD3,ISYMB) + NCKATR(ISCKD3)*(B-1) + 1
                  IF (NCKATR(ISCKD3) .GT. 0) THEN
                     CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVIB0),IOFF,
     &                           NCKATR(ISCKD3))
C
                     IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                        XTRVI0= DDOT(NCKATR(ISCKD3),WORK(KTRVIB0),1,
     *                               WORK(KTRVIB0),1)
                        WRITE(LUPRI,*) 'Norm of TRVIB0 readin',XTRVI0
                     ENDIF
                  ENDIF
                  DTIME  = SECOND() - DTIME
                  TIIO = TIIO   + DTIME
C
C---------------------------
C                 Sort.
C---------------------------
C
                  DTIME = SECOND()
                  CALL CCSDT_SRTVIR(WORK(KTRVIB0),WORK(KTRVIB2),
     *                              WORK(KEND4),LWRK4,ISYMB,ISINT22)
C
                  DTIME  = SECOND() - DTIME
                  TISORT = TISORT   + DTIME
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XTRVI0= DDOT(NCKATR(ISCKD3),WORK(KTRVIB0),1,
     *                            WORK(KTRVIB0),1)
                     WRITE(LUPRI,*) 'Norm of TRVIB0 ',XTRVI0
C
                     XTRVI2= DDOT(NCKATR(ISCKD3),WORK(KTRVIB2),1,
     *                            WORK(KTRVIB2),1)
                     WRITE(LUPRI,*) 'Norm of TRVIB2 ',XTRVI2
                  ENDIF
C
C-----------------------------------------
C                 Readin part2.
C-----------------------------------------
C
                  DTIME = SECOND()
                  IOFF = ICKBD(ISCKD2,ISYMB) + NCKATR(ISCKD2)*(B-1) + 1
                  IF (NCKATR(ISCKD2) .GT. 0) THEN
                     CALL GETWA2(LUDKBC2,FNDKBC2,WORK(KTRVIB4),IOFF,
     &                           NCKATR(ISCKD2))
C
                     IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                        XTRVI0= DDOT(NCKATR(ISCKD2),WORK(KTRVIB4),1,
     *                               WORK(KTRVIB4),1)
                        WRITE(LUPRI,*) 'Norm of TRVIB4 readin',XTRVI0
                     ENDIF
                  ENDIF
                  DTIME  = SECOND() - DTIME
                  TIIO = TIIO   + DTIME
C
C-------------------------------------------------
C                 Sort the integrals part 2.
C-------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CCSDT_SRTVIR(WORK(KTRVIB4),WORK(KTRVIB5),
     *                              WORK(KEND4),LWRK4,ISYMB,ISINT2)
C
                  DTIME  = SECOND() - DTIME
                  TISORT = TISORT   + DTIME
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XTRVI2= DDOT(NCKATR(ISCKD2),WORK(KTRVIB5),1,
     *                            WORK(KTRVIB5),1)
                     WRITE(LUPRI,*) 'Norm of TRVI5 ',XTRVI2
                  ENDIF
C
C-----------------------------------------------------------
C                 Read virtual integrals part 2b.
C-----------------------------------------------------------
C
                  DTIME = SECOND()
                  IOFF = ICKBD(ISCKD2,ISYMB) + NCKATR(ISCKD2)*(B-1) + 1
                  IF (NCKATR(ISCKD2) .GT. 0) THEN
                     CALL GETWA2(LUDKBC4,FNDKBC4,WORK(KTRVIB6),IOFF,
     &                           NCKATR(ISCKD2))
                  ENDIF
                  DTIME  = SECOND() - DTIME
                  TIIO = TIIO   + DTIME
C
C----------------------------------------------------
C                 Sort the integrals for s3am part 2b.
C----------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CCSDT_SRTVIR(WORK(KTRVIB6),WORK(KTRVIB4),
     &                              WORK(KEND4),LWRK4,ISYMB,ISINT2)
C
                  DTIME  = SECOND() - DTIME
                  TISORT = TISORT   + DTIME
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XTRVI0= DDOT(NCKATR(ISCKD2),WORK(KTRVIB4),1,
     *                            WORK(KTRVIB4),1)
                     WRITE(LUPRI,*) 'Norm of TRVI4 ',XTRVI0
C
                     XTRVI2= DDOT(NCKATR(ISCKD2),WORK(KTRVIB5),1,
     *                            WORK(KTRVIB5),1)
                     WRITE(LUPRI,*) 'Norm of TRVI5 ',XTRVI2
                  ENDIF
C
C----------------------------------------------------
C                 Read some more integrals.
C----------------------------------------------------
C
                  DTIME = SECOND()
                  CALL DZERO(WORK(KTRVIB6),
     &                       MAX(NCKATR(ISCKD2),NCKATR(ISCKD3)))
                  IOFF = ICKBD(ISCKD3,ISYMB) + NCKATR(ISCKD3)*(B-1) + 1
                  IF (NCKATR(ISCKD3) .GT. 0) THEN
                     CALL GETWA2(LUDKBC5,FNDKBC5,WORK(KTRVIB6),IOFF,
     &                           NCKATR(ISCKD3))
                  ENDIF
                  DTIME  = SECOND() - DTIME
                  TIIO = TIIO   + DTIME
C
C---------------------------------------------------
C                 Sort again.
C---------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CCSDT_SRTVIR(WORK(KTRVIB6),WORK(KTRVIB7),
     *                              WORK(KEND4),LWRK4,ISYMB,ISINT22)
C
                  DTIME  = SECOND() - DTIME
                  TISORT = TISORT   + DTIME
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XTRVI2= DDOT(NCKATR(ISCKD3),WORK(KTRVIB7),1,
     *                            WORK(KTRVIB7),1)
                     WRITE(LUPRI,*) 'Norm of TRVI7 ',XTRVI2
                  ENDIF
C
C---------------------------------------------
C                 Read for the last time.
C---------------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DZERO(WORK(KTRVIB6),
     &                       MAX(NCKATR(ISCKD2),NCKATR(ISCKD3)))
                  IOFF = ICKBD(ISCKD2,ISYMB) + NCKATR(ISCKD2)*(B-1) + 1
                  IF (NCKATR(ISCKD2) .GT. 0) THEN
                     CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVIB6),IOFF,
     &                           NCKATR(ISCKD2))
                  ENDIF
C
C---------------------------------------------------------------
C                 Read the contraction integrals for given B
C---------------------------------------------------------------
C
                  IOFF = ICKBD(ISYCKD,ISYMB)+NCKATR(ISYCKD)*(B - 1)+1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LU3VI3,FN3VI3,WORK(KTRVI8),IOFF,
     &                           NCKATR(ISYCKD))
                  ENDIF
                  IOFF = ICKBD(ISYCKD,ISYMB)+NCKATR(ISYCKD)*(B - 1)+1
                  IF (NCKATR(ISYCKD) .GT. 0) THEN
                     CALL GETWA2(LU3VI4,FN3VI4,WORK(KTRVI9),IOFF,
     &                           NCKATR(ISYCKD))
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TIIO = TIIO   + DTIME
C
C-------------------------------------------------------
C                 Calculate the part of the R3 matrix
C-------------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_SMAT3(ECURR,T2TP,ISYMT2,C2AMP,C2AMM,ISYMC2,
     *                           WORK(KTMAT),WORK(KTRVI0),
     *                           WORK(KTRVI2),
     *                           WORK(KTRVI4),WORK(KTRVI5),
     *                           WORK(KTRVI6),WORK(KTRVI7),
     *                           WORK(KTROC0),WORK(KTROC2),
     *                           WORK(KTROC3),WORK(KTROC4),
     *                           WORK(KTROC5),ISINT2,
     *                           ISINT22,WORK(KFOCKD),WORK(KDIAG),
     *                           WORK(KSMAT),WORK(KEND4),LWRK4,
     *                           WORK(KINDEX),WORK(KINDEX2),
     *                           WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                  CALL CC3_SMAT3(ECURR,T2TP,ISYMT2,C2AMP,C2AMM,ISYMC2,
     *                           WORK(KTMAT),WORK(KTRVIB0),
     *                           WORK(KTRVIB2),
     *                           WORK(KTRVIB4),WORK(KTRVIB5),
     *                           WORK(KTRVIB6),WORK(KTRVIB7),
     *                           WORK(KTROC0),WORK(KTROC2),
     *                           WORK(KTROC3),WORK(KTROC4),
     *                           WORK(KTROC5),ISINT2,
     *                           ISINT22,WORK(KFOCKD),WORK(KDIAG),
     *                           WORK(KSMAT2),WORK(KEND4),LWRK4,
     *                           WORK(KINDEX3),WORK(KINDEX4),
     *                           WORK(KINDSQ),LENSQ,ISYMD,D,ISYMB,B)
C
                  CALL DAXPY(NCKIJ(ISCKIJ),XMONE,WORK(KSMAT2),1,
     *                       WORK(KSMAT),1)
C
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (LDEBUG) THEN
C
CKH   DO NOT REMOVE THIS CALL. USED FOR DEBUGGING.
CKH   USES A STATIC ALL. ARRAY OF N^6 THUS THE SUBROUTINE
CKH   (AND THEREFORE ALL CALLS) SHOULD BE OUTCOMMENTED
C     
C                     CALL SUM_R3(WORK(KSMAT),ISYMD,D,ISYMB,B,
C     *                           NCKIJ(ISCKIJ),WORK(KR3))
C
                  ENDIF
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT ',XSMAT
                  ENDIF
C
C-----------------------------------------------
C                 Contract with integrals.
C-----------------------------------------------
C
C
                  DTIME = SECOND()
                  IOPT = 1
                  CALL CC3_CONVIR33(OMEGA2P,WORK(KRES2P),OMEGA2M,
     *                              WORK(KRES2M),WORK(KRMAT1),
     *                              WORK(KRMAT2),WORK(KRMAT3),
     *                              WORK(KRMAT4),WORK(KSMAT),
     *                              WORK(KTMAT),ISYMIM,WORK(KTRVI),
     *                              WORK(KTRVI1),WORK(KTRVI8),
     *                              WORK(KTRVI9),ISINT1,WORK(KEND4),
     *                              LWRK4,WORK(KINDSQ),LENSQ,
     *                              ISYMB,B,ISYMD,D,IOPT,
     *                              TIME1P,TIME2P,TIME3P,TIME1M,
     *                              TIME2M,TSORTP,TSORTM)
C
                  IOPT = 2
                  CALL CC3_CONVIR33(OMEGA2P,WORK(KRES2P),OMEGA2M,
     *                              WORK(KRES2M),WORK(KRMAT1),
     *                              WORK(KRMAT2),WORK(KRMAT3),
     *                              WORK(KRMAT4),WORK(KSMAT),
     *                              WORK(KTMAT),ISYMIM,WORK(KTRVI),
     *                              WORK(KTRVI1),WORK(KTRVI8),
     *                              WORK(KTRVI9),ISINT1,WORK(KEND4),
     *                              LWRK4,WORK(KINDSQ),LENSQ,
     *                              ISYMB,B,ISYMD,D,IOPT,
     *                              TIME1P,TIME2P,TIME3P,TIME1M,
     *                              TIME2M,TSORTP,TSORTM)
C
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'RMAT1 norm -after CONVIR33',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'RMAT2 norm -after CONVIR33',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT3),1,
     *                       WORK(KRMAT3),1)
                     WRITE(LUPRI,*) 'RMAT3 norm -after CONVIR33',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT4),1,
     *                       WORK(KRMAT4),1)
                     WRITE(LUPRI,*) 'RMAT4 norm -after CONVIR33',XRMAT
                  ENDIF
C
                  IF (IPRINT .GT. 220 .OR. LDEBUG) THEN
                     CALL AROUND('After CC3_CONVIR: Rho+')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,0,1)
                     CALL AROUND('After CC3_CONVIR: Rho-')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TICONV = TICONV   + DTIME
C
                  DTIME = SECOND()
C
                  IOPT = 1
                  CALL CC3_CONOCC33(OMEGA2P,OMEGA2M,WORK(KRMAT1),
     *                              WORK(KRMAT2),WORK(KRMAT3),
     *                              WORK(KRMAT4),WORK(KSMAT),
     *                              WORK(KTMAT),ISYMIM,WORK(KTROC),
     *                              WORK(KTROC1),ISINT1,WORK(KEND4),
     *                              LWRK4,WORK(KINDSQ),LENSQ,ISYMB,B,
     *                              ISYMD,D,IOPT)
                  IOPT = 2
                  CALL CC3_CONOCC33(OMEGA2P,OMEGA2M,WORK(KRMAT1),
     *                              WORK(KRMAT2),WORK(KRMAT3),
     *                              WORK(KRMAT4),WORK(KSMAT),
     *                              WORK(KTMAT),ISYMIM,WORK(KTROC),
     *                              WORK(KTROC1),ISINT1,WORK(KEND4),
     *                              LWRK4,WORK(KINDSQ),LENSQ,ISYMB,B,
     *                              ISYMD,D,IOPT)
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'RMAT1 norm -after CONOCC33',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'RMAT2 norm -after CONOCC33',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT3),1,
     *                       WORK(KRMAT3),1)
                     WRITE(LUPRI,*) 'RMAT3 norm -after CONOCC33',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT4),1,
     *                       WORK(KRMAT4),1)
                     WRITE(LUPRI,*) 'RMAT4 norm -after CONOCC33',XRMAT
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Rho1 norm after CC3_CONOCC33',
     *                                                       RHO1N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
                     WRITE(LUPRI,*) 'Rho2+ norm after CC3_CONOCC33',
     *                                                       RHO2N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
                     WRITE(LUPRI,*) 'Rho2- norm after CC3_CONOCC33',
     *                                                       RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220 .OR. LDEBUG) THEN
                     CALL AROUND('After CC3_CONOCC33: Rho2+ ')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,0,1)
                     CALL AROUND('After CC3_CONOCC33: Rho2- ')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TICONO = TICONO   + DTIME
C
C--------------------------------------------------------------------
C                 Calculate Omega1 and the Fock cont. to omega2P/M
C--------------------------------------------------------------------
C
                  DTIME = SECOND()
C
                  IOPT = 1
                  CALL CC3_ONEL33(OMEGA1,OMEGA2P,OMEGA2M,WORK(KRMAT1),
     *                            WORK(KRMAT2),WORK(KRMAT3),WORK(KRMAT4)
     *                            ,WORK(KFCKAK),WORK(KSMAT),
     *                            WORK(KTMAT),ISYMIM,WORK(KXIAJB),
     *                            WORK(KXIAJB2),
     *                            ISINT1,WORK(KINDSQ),LENSQ,WORK(KEND4),
     *                            LWRK4,ISYMB,B,ISYMD,D,IOPT)
                  IOPT = 2
                  CALL CC3_ONEL33(OMEGA1,OMEGA2P,OMEGA2M,WORK(KRMAT1),
     *                            WORK(KRMAT2),WORK(KRMAT3),WORK(KRMAT4)
     *                            ,WORK(KFCKAK),WORK(KSMAT),
     *                            WORK(KTMAT),ISYMIM,WORK(KXIAJB),
     *                            WORK(KXIAJB2),
     *                            ISINT1,WORK(KINDSQ),LENSQ,WORK(KEND4),
     *                            LWRK4,ISYMB,B,ISYMD,D,IOPT)
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 -after ONEL33',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 -after ONEL33',XRMAT
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT -after ONEL33',XSMAT
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KTMAT),1)
                     WRITE(LUPRI,*) 'Norm of TMAT -after ONEL33',XTMAT
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_ONEL33',
     *                                                    RHO1N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
                     WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_ONEL33',
     *                                                    RHO2N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
                     WRITE(LUPRI,*) 'Norm of Rho2- after CC3_ONEL33',
     *                                                    RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220 .OR. LDEBUG) THEN
                     CALL AROUND('After CC3_ONEL33: Rho1 & Rho2+ ')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,1,1)
                     CALL AROUND('After CC3_ONEL:33 Rho2- ')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TIOME1 = TIOME1   + DTIME
C
C----------------------------------------------------
C                 Accumulate the R2 matrix in Omega2.
C----------------------------------------------------
C
                  CALL CC3_RACC(OMEGA2P,WORK(KRMAT2),ISYMB,B,ISYRES)
                  IOPT = 2
                  FACTOR = ONE
                  CALL CC3_RACC3(DUMMY,OMEGA2M,WORK(KRMAT4),ISYMB,B,
     *                          ISYRES,FACTOR,IOPT)
C
                  IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_RACC',RHO1N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
                     WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_RACC',RHO2N
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
                     WRITE(LUPRI,*) 'Norm of Rho2- after CC3_RACC',RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220 .OR. LDEBUG) THEN
                     CALL AROUND('After CC3_RACC: Rho1 & Rho2+')
                     CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,1,1)
                     CALL AROUND('After CC3_RACC: Rho2-')
                     CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
                  ENDIF
C
  130          CONTINUE
C
  120       CONTINUE
C
C----------------------------------------------
C           Accumulate the R1 matrix in Omega2.
C----------------------------------------------
C
            CALL CC3_RACC(OMEGA2P,WORK(KRMAT1),ISYMD,D,ISYRES)
            IOPT = 2
            FACTOR = ONE
            CALL CC3_RACC3(DUMMY,OMEGA2M,WORK(KRMAT3),ISYMD,D,ISYRES,
     *                     FACTOR,IOPT)
C
            IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
               RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
               WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_RACC-2',RHO1N
               RHO2N = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
               WRITE(LUPRI,*) 'Norm of Rho2+ after CC3_RACC-2',RHO2N
               RHO2N = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
               WRITE(LUPRI,*) 'Norm of Rho2- after CC3_RACC-2',RHO2N
            ENDIF
C
            IF (IPRINT .GT. 220 .OR. LDEBUG) THEN
               CALL AROUND('After CC3_RACC-2: Rho1 & Rho2+ ')
               CALL CC_PRP(OMEGA1,OMEGA2P,ISYRES,1,1)
               CALL AROUND('After CC3_RACC-2: Rho2- ')
               CALL CC_PRP(OMEGA1,OMEGA2M,ISYRES,0,1)
            ENDIF
C
  110    CONTINUE
C
C-----------------------------------------------------------
C        Accumulate the result vectors from cc3_convir33
C        and last loop stop.
C-----------------------------------------------------------
C
         DTIME = SECOND()
         CALL CC3_SORTPLUS(OMEGA2P,WORK(KRES2P),ISYRES)
         CALL CC3_SORTMINUS(OMEGA2M,WORK(KRES2M),ISYRES)
         DTIME = SECOND() - DTIME
         TISORT = TISORT + DTIME
C
  100 CONTINUE
C
C-----------------------------------------
C     Close and delete files.
C-----------------------------------------
C
      CALL WCLOSE2(LUDKBC2,FNDKBC2,'DELETE')
      CALL WCLOSE2(LUDKBC3,FNDKBC3,'DELETE')
      CALL WCLOSE2(LUDKBC4,FNDKBC4,'DELETE')
      CALL WCLOSE2(LUDKBC5,FNDKBC5,'DELETE')
      CALL WCLOSE2(LU3VI3,FN3VI3,'DELETE')
      CALL WCLOSE2(LU3VI4,FN3VI4,'DELETE')
      CALL WCLOSE2(LUDELD,FNDELD,'DELETE')
      CALL WCLOSE2(LUDKBC,FNDKBC,'DELETE')
C
C----------------------------------------
C     Print timings and R3 if LDEBUG.
C----------------------------------------
C
      IF (LDEBUG) THEN
C
CKH   DO NOT REMOVE THIS CALL. USED FOR DEBUGGING.
CKH   USES A STATIC ALL. ARRAY OF N^6 THUS THE SUBROUTINE
CKH   (AND THEREFORE ALL CALLS) SHOULD BE OUTCOMMENTED
C
C         WRITE(LUPRI,*) 'THE R3 AMPLITUDES : '
C         CALL PRINT_R3(WORK(KR3),ISYMIM)
C
      ENDIF
C
      IF (IPRINT .GT. 9 .OR. LDEBUG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,1) 'CC3_TRAN  : ',TITRAN
         WRITE(LUPRI,1) 'CC3_SORT  : ',TISORT
         WRITE(LUPRI,1) 'CC3_SMAT  : ',TISMAT
         WRITE(LUPRI,1) 'CC3_CONV  : ',TICONV
         WRITE(LUPRI,2) '1.ST +    : ',TIME1P
         WRITE(LUPRI,2) '2.ND +    : ',TIME2P
         WRITE(LUPRI,2) '3.RD +    : ',TIME3P
         WRITE(LUPRI,2) '1.ST -    : ',TIME1M
         WRITE(LUPRI,2) '2.ND -    : ',TIME2M
         WRITE(LUPRI,1) 'CC3_CONO  : ',TICONO
         WRITE(LUPRI,1) 'CC3_OME1  : ',TIOME1
         WRITE(LUPRI,1) 'CC3_IO    : ',TIIO
         WRITE(LUPRI,1) 'CC3_IVINT : ',TIVINT
         WRITE(LUPRI,*)
      END IF
C
      CALL QEXIT('CC3_OMEG33')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
    2 FORMAT(7X,'CC3_CONV consists of',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck cc3_smat3 */
      SUBROUTINE CC3_SMAT3(ECURR,T2TP,ISYMT2,C2AMP,C2AMM,ISYMC2,
     *                     TMAT,TRVIR,
     *                     TRVIR2,TRVIR4,TRVIR5,TRVIR6,TRVIR7,TROCC,
     *                     TROCC2,TROCC3,TROCC4,TROCC5,ISYINT,ISYINT2,
     *                     FOCKD,DIAG,SMAT,WORK,LWORK,INDAJL,INDAJL2,
     *                     INDSQ,LENSQ,ISYMB,B,ISYMC,C)
C
C     Written by Kasper Hald, Spring 2001.
C     Based on CC3_SMAT by
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C     Generalised for symmetry: Ove Christiansen Nov. 1995
C                                                Jan. 1996
C
C     Calculate the S matrix for approximate triples.
C
C     S(ck,bj,ai) = sum(d)[
C                 + t(ck,di)g(1)(ajbd) + t(bj,di)g(2)(ckad)
C                 + t(bi,dk)g(2)(ajcd) + r+(di,bj) (ad|ck)
C                 + r-(ai,dk) (bj|cd)  + r-(ck,di) (aj|bd)  ]
C                 - sum(l)[
C                 + t(ck,al)g(3)(bilj) + t(bj,al)g(4)(ckli)
C                 + t(bi,cl)g(4)(ajlk) + r+(al,bj) (ck|li)
C                 - r-(cl,ai) (bj|lk)  - r-(bl,ck) (aj|li) ]
C
C     S is stored as S(ai,k,j) for fixed b and c
C
C
C     General symmetry: ISYINT is symmetry of integrals for T2, 
C                       ISYMT2 is the symmetry of T2TP,
C                       ISYINT2 is the symmetry of integrals for C2,
C                   and ISYMC2 is the symmetry of C2AMP and C2AMM.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMT2, ISYMC2, ISYINT, ISYINT2, LENSQ, INDSQ(LENSQ,6)
      INTEGER ISYMB, ISYMC, LENGTH, ISYMA, ISYMAJ, ISYAIK, ISYMBC
      INTEGER KOFF, KOFF1, KOFF2, KOFF3
      INTEGER NTOAIJ, NVIRD, ISYAKD, ISYDIJ, ISYMJ, ISYMDI, ISYMI
      INTEGER ISYMAK, ISYAKI, NTOTAK, ISYAIL, ISYLKJ, ISYMLK
      INTEGER NTOTAI, NRHFL, ISYAJL, ISYLKI, NB, NC, NTOTAJ, ISYAJK
      INTEGER ISYMAI, ISYML, ISYAIJ, ISYMD, ISYMK, ISYMDK, ISYMDK2
      INTEGER JSAIKJ, INDAJL(*), INDAJL2(*), LWORK, ISYRES, ISYRES2
      INTEGER ISYBJL, ISYCJL, ISYMBJ, ISYMCJ, ISYMJL
      INTEGER NTOAIK, NTOTL
C
#if defined (SYS_CRAY)
      REAL TRVIR(*),TRVIR2(*),TROCC(*),SMAT(*),FOCKD(*)
      REAL C2AMP(*), C2AMM(*), XSMAT, DDOT, EPSIBC
      REAL DIAG(*),T2TP(*),TMAT(*), TROCC2(*), TROCC3(*)
      REAL TROCC4(*), TROCC5(*)
      REAL TRVIR4(*),TRVIR5(*),TRVIR6(*),TRVIR7(*)
      REAL WORK(LWORK), ZERO, ONE, XMONE, TWO, HALF
      REAL FACTOR, XMHALF, XDOT, ECURR
#else
      DOUBLE PRECISION TRVIR(*),TRVIR2(*),TROCC(*),SMAT(*),FOCKD(*)
      DOUBLE PRECISION C2AMP(*), C2AMM(*), XSMAT, DDOT, EPSIBC
      DOUBLE PRECISION DIAG(*),T2TP(*),TMAT(*), TROCC2(*), TROCC3(*)
      DOUBLE PRECISION TROCC4(*), TROCC5(*)
      DOUBLE PRECISION TRVIR4(*),TRVIR5(*),TRVIR6(*),TRVIR7(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE, XMONE, TWO, HALF
      DOUBLE PRECISION FACTOR, XMHALF, XDOT, ECURR
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, XMONE = -1.0D0)
      PARAMETER(TWO = 2.0D0, HALF = 0.5D0, XMHALF = -0.5D0)
C
      LOGICAL LDEBUG
C
      CALL QENTER('CC3_SMAT3')
C
      LDEBUG = .FALSE.
C
C--------------------------
C     Initialitation.
C--------------------------
C
      ISYMBC  = MULD2H(ISYMB,ISYMC)
      ISYRES  = MULD2H(ISYINT,ISYMT2)
      ISYRES2 = MULD2H(ISYINT2,ISYMC2)
      JSAIKJ  = MULD2H(ISYMBC,ISYRES)
      ISYMDK  = MULD2H(ISYMBC,ISYINT)
      ISYMDK2 = MULD2H(ISYMBC,ISYINT2)
C
      IF (ISYRES .NE. ISYRES2) THEN
         CALL QUIT('Symmetry error in CC3_SMAT3')
      ENDIF
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient core in CC3_SMAT3')
      ENDIF
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 'Entered CC3_SMAT3 with B = ',B,'   &   C = ',C
         WRITE(LUPRI,*) 'ISYMB = ',ISYMB,' and ISYMC = ',ISYMC
         WRITE(LUPRI,*) 'ISYMT2 = ',ISYMT2,' and ISYMC2 = ',ISYMC2
         WRITE(LUPRI,*) 'ISYINT = ',ISYINT,' and ISYINT2 = ',ISYINT2
         WRITE(LUPRI,*) 'JSAIKJ = ',JSAIKJ,' and ISYRES=ISYRES2=',ISYRES
         WRITE(LUPRI,*) 'ISYMDK = ',ISYMDK,' and ISYMDK2 = ',ISYMDK2
C
         XDOT = DDOT(NT2SQ(ISYMT2),T2TP,1,T2TP,1)
         WRITE(LUPRI,*) 'Norm of T2  Amplitudes in SMAT3 = ',XDOT
         XDOT = DDOT(NT2SQ(ISYMC2),C2AMP,1,C2AMP,1)
         WRITE(LUPRI,*) 'Norm of R2+ Amplitudes in SMAT3 = ',XDOT
         XDOT = DDOT(NT2SQ(ISYMC2),C2AMM,1,C2AMM,1)
         WRITE(LUPRI,*) 'Norm of R2- Amplitudes in SMAT3 = ',XDOT
      ENDIF
C
C---------------------------------------------
C     First virtual contribution for T2.
C---------------------------------------------
C
      FACTOR = HALF
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
     *              NVIR(ISYMD),factor,T2TP(KOFF1),NTOAIJ,
     *              TRVIR6(KOFF2),NVIRD,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
  100 CONTINUE
C
      CALL DZERO(TMAT,LENGTH)
      CALL CC_GATHER(LENGTH,TMAT,WORK,INDSQ(1,3))
      CALL DAXPY(LENGTH,ONE,TMAT,1,SMAT,1)
C
      CALL DZERO(TMAT,LENGTH)
      CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      CALL DAXPY(LENGTH,XMONE,TMAT,1,SMAT,1)
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 1. Norm of SMAT ',XSMAT
      ENDIF
C
C--------------------------------------------------
C     First virtual contribution from C2+.
C--------------------------------------------------
C
      FACTOR = ONE
      DO ISYMK = 1,NSYM
C
         ISYMD  = MULD2H(ISYMK,ISYMDK2)
         ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIJ,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDK2,ISYMB) + NT1AM(ISYMDK2)*(B - 1)
     *         + IT1AM(ISYMD,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
         NTOAIJ = MAX(1,NCKI(ISYAIJ))
         NVIRD  = MAX(NVIR(ISYMD),1)
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *              NVIR(ISYMD),factor,C2AMM(KOFF1),NTOAIJ,
     *              TRVIR(KOFF2),NVIRD,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
      ENDDO
C
      CALL DZERO(TMAT,LENGTH)
      CALL CC_GATHER(LENGTH,TMAT,WORK,INDSQ(1,3))
      CALL DAXPY(LENGTH,ONE,TMAT,1,SMAT,1)
C
      CALL DZERO(TMAT,LENGTH)
      CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      CALL DAXPY(LENGTH,XMONE,TMAT,1,SMAT,1)
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 2. Norm of SMAT ',XSMAT
      ENDIF
C
C---------------------------------------------
C     Second virtual contribution for T2.
C---------------------------------------------
C
      ISYAKD = MULD2H(ISYMC,ISYINT)
      ISYDIJ = MULD2H(ISYMB,ISYMT2)
C
      FACTOR = HALF
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
     *                    NVIR(ISYMD),factor,TRVIR5(KOFF1),NTOTAK,
     *                    T2TP(KOFF2),NVIRD,ZERO,TMAT(KOFF3),
     *                    NTOTAK)
C
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) + TMAT(INDSQ(I,1)) - TMAT(INDSQ(I,4))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 3. Norm of SMAT ',XSMAT
      ENDIF
C
C---------------------------------------------
C     Second virtual contribution for R2-
C---------------------------------------------
C
      ISYAKD = MULD2H(ISYMC,ISYINT2)
      ISYDIJ = MULD2H(ISYMB,ISYMC2)
C
      FACTOR = ONE
C
      DO ISYMJ = 1,NSYM
C
         ISYMDI = MULD2H(ISYMJ,ISYDIJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
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
     *                    NVIR(ISYMD),factor,TRVIR2(KOFF1),NTOTAK,
     *                    C2AMM(KOFF2),NVIRD,ZERO,TMAT(KOFF3),
     *                    NTOTAK)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) + TMAT(INDSQ(I,1)) - TMAT(INDSQ(I,4))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 4. Norm of SMAT ',XSMAT
      ENDIF
C
C---------------------------------------------
C     Third virtual contribution for T2.
C---------------------------------------------
C
      ISYAKD = MULD2H(ISYMC,ISYINT)
      ISYDIJ = MULD2H(ISYMB,ISYMT2)
C
      FACTOR = XMHALF
C
      DO ISYMJ = 1,NSYM
C
         ISYMDI = MULD2H(ISYMJ,ISYDIJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
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
     *                    NVIR(ISYMD),factor,TRVIR4(KOFF1),NTOTAK,
     *                    T2TP(KOFF2),NVIRD,ZERO,TMAT(KOFF3),
     *                    NTOTAK)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) + TMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 5. Norm of SMAT ',XSMAT
      ENDIF
C
C-----------------------------------
C     Virtual contribution for R2+
C-----------------------------------
C
      ISYAKD = MULD2H(ISYMC,ISYINT2)
      ISYDIJ = MULD2H(ISYMB,ISYMC2)
C
      FACTOR = XMHALF
C
      DO ISYMJ = 1,NSYM
C
         ISYMDI = MULD2H(ISYMJ,ISYDIJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
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
     *                    NVIR(ISYMD),factor,TRVIR7(KOFF1),NTOTAK,
     *                    C2AMP(KOFF2),NVIRD,ZERO,TMAT(KOFF3),
     *                    NTOTAK)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) + TMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 6. Norm of SMAT ',XSMAT
      ENDIF
C
C--------------------------------------------
C     First occupied contribution from T2.
C--------------------------------------------
C
      ISYAIL = MULD2H(ISYMB,ISYMT2)
      ISYLKJ = MULD2H(ISYMC,ISYINT)
C
      CALL DZERO(TMAT,LENGTH)
      FACTOR = HALF
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
     *                    NRHF(ISYML),factor,T2TP(KOFF1),NTOTAI,
     *                    TROCC2(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                    NTOTAI)
C
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
C
      DO I = 1, LENGTH
         SMAT(I) = SMAT(I) + TMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 7. Norm of SMAT ',XSMAT
      ENDIF
C
C--------------------------------------------
C     First occupied contribution from R2-.
C--------------------------------------------
C
      ISYAIL = MULD2H(ISYMB,ISYMC2)
      ISYLKJ = MULD2H(ISYMC,ISYINT2)
C
      CALL DZERO(TMAT,LENGTH)
      FACTOR = XMONE
C
      DO ISYMJ = 1,NSYM
C
         ISYMLK = MULD2H(ISYMJ,ISYLKJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMK = 1,NSYM
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
     *                    NRHF(ISYML),factor,C2AMM(KOFF1),NTOTAI,
     *                    TROCC(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO ! ISYMK
         ENDDO ! J
      ENDDO !ISYMJ
C
      DO I = 1, LENGTH
         SMAT(I) = SMAT(I) + TMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 8. Norm of SMAT ',XSMAT
      ENDIF
C
C---------------------------------------------
C     Second occupied contribution for T2.
C---------------------------------------------
C
      ISYAJL = MULD2H(ISYMB,ISYMT2)
      ISYLKI = MULD2H(ISYMC,ISYINT)
C
      IF (LWORK .LT. NCKI(ISYAJL)) THEN
         CALL QUIT('Not enough space in CC3_SMAT3')
      END IF
C
      FACTOR = HALF
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
     *                    NRHF(ISYML),factor,WORK(KOFF1),NTOTAJ,
     *                    TROCC3(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAJ)
C
  420       CONTINUE
  410    CONTINUE
  400 CONTINUE
C
      DO I = 1,NCKIJ(JSAIKJ)
         SMAT(I) = SMAT(I) + TMAT(INDSQ(I,2)) - TMAT(INDSQ(I,5))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 9. Norm of SMAT ',XSMAT
      ENDIF
C
C---------------------------------------------
C     Second occupied contribution for R2-.
C---------------------------------------------
C
      ISYAJL = MULD2H(ISYMB,ISYMC2)
      ISYLKI = MULD2H(ISYMC,ISYINT2)
C
      IF (LWORK .LT. NCKI(ISYAJL)) THEN
         CALL QUIT('Not enough space in CC3_SMAT3')
      END IF
C
      FACTOR = ONE
C
      KOFF = IT2SP(ISYAJL,ISYMB) + NCKI(ISYAJL)*(B - 1) + 1
      CALL CC_GATHER(NCKI(ISYAJL),WORK,C2AMM(KOFF),INDAJL2)
C
      DO ISYMI = 1,NSYM
C
         ISYMLK = MULD2H(ISYMI,ISYLKI)
C
         DO I = 1,NRHF(ISYMI)
C
            DO ISYMK = 1,NSYM
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
     *                    NRHF(ISYML),factor,WORK(KOFF1),NTOTAJ,
     *                    TROCC(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAJ)
C
            ENDDO ! ISYMK
         ENDDO    ! I
      ENDDO       ! ISYMI
C
      DO I = 1,NCKIJ(JSAIKJ)
         SMAT(I) = SMAT(I) + TMAT(INDSQ(I,2)) - TMAT(INDSQ(I,5))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 10. Norm of SMAT ',XSMAT
      ENDIF
C
C--------------------------------------------
C     First occupied contribution from R2+.
C--------------------------------------------
C
      ISYBJL = MULD2H(ISYMC,ISYMC2)
      ISYCJL = MULD2H(ISYMB,ISYMC2)
      ISYMJL  = MULD2H(ISYMC,ISYCJL)
      IF (LWORK .LE. NMATIJ(ISYMJL))
     *    CALL QUIT('OUT OF WORKSPACE IN CC33_SMAT3 (Sort)')
C
       CALL CC33_T2SORT(C2AMP,ISYMC2,WORK,ISYMB,B,ISYMC,C)
C
      CALL DZERO(TMAT,LENGTH)
      FACTOR = HALF
C
      DO ISYMJ = 1,NSYM
C
         ISYML  = MULD2H(ISYMJL,ISYMJ)
         ISYAIK = MULD2H(ISYINT2,ISYML)
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMCJ = MULD2H(ISYMC,ISYMJ)
C
         KOFF1 = ISAIKJ(ISYAIK,ISYML)
     *         + 1
         KOFF2 = IMATIJ(ISYML,ISYMJ)
     *         + 1
         KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *         + 1
C
         NTOAIK = MAX(1,NCKI(ISYAIK))
         NTOTL  = MAX(1,NRHF(ISYML))
C
         CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *              NRHF(ISYML),factor,TROCC4(KOFF1),NTOAIK,
     *              WORK(KOFF2),NTOTL,ONE,TMAT(KOFF3),
     *              NTOAIK)
C
      ENDDO
C
      DO I = 1, LENGTH
         SMAT(I) = SMAT(I) + TMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 11. Norm of SMAT ',XSMAT
      ENDIF
C
C--------------------------------------------
C     Third occupied contribution from T2.
C--------------------------------------------
C
      ISYAIL = MULD2H(ISYMB,ISYMT2)
      ISYLKJ = MULD2H(ISYMC,ISYINT)
C
      ISYBJL = MULD2H(ISYMC,ISYMT2)
      ISYCJL = MULD2H(ISYMB,ISYMT2)
      ISYMJL  = MULD2H(ISYMC,ISYCJL)
       IF (LWORK .LE. NMATIJ(ISYMJL))
     *     CALL QUIT('OUT OF WORKSPACE IN CC33_smat3 (SORT)')
C
       CALL CC33_T2SORT(T2TP,ISYMT2,WORK,ISYMB,B,ISYMC,C)
C
      CALL DZERO(TMAT,LENGTH)
      FACTOR = HALF
C
      DO ISYMJ = 1,NSYM
C
         ISYMLK = MULD2H(ISYMJ,ISYLKJ)
         ISYML  = MULD2H(ISYMJL,ISYMJ)
         ISYAIK = MULD2H(ISYINT,ISYML)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMCJ = MULD2H(ISYMC,ISYMJ)
C
         KOFF1 = ISAIKJ(ISYAIK,ISYML)
     *         + 1
         KOFF2 = IMATIJ(ISYML,ISYMJ)
     *         + 1
         KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *         + 1
         NTOAIK = MAX(1,NCKI(ISYAIK))
         NTOTL = MAX(1,NRHF(ISYML))
C
         CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *              NRHF(ISYML),factor,TROCC5(KOFF1),NTOAIK,
     *              WORK(KOFF2),NTOTL,ONE,TMAT(KOFF3),NTOAIK)
C
      ENDDO
C
      DO I = 1, LENGTH
         SMAT(I) = SMAT(I) + TMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT3: 12. Norm of SMAT ',XSMAT
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
         SMAT(L) = -SMAT(L)/(DIAG(L) + EPSIBC)
C
  500 CONTINUE
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: Norm of SMAT (end) ',XSMAT
      ENDIF
C
C
      CALL QEXIT('CC3_SMAT3')
C
      RETURN
      END
C  /* Deck cc3_sort3 */
      SUBROUTINE CC3_SORT3(WORK,LWORK,IOPT,ISYINT,LU3SRT2,FN3SRT2,
     *                     LU3SRT3,FN3SRT3,LUDELD2,FNDELD2,
     *                     LUDELD3,FNDELD3)
C
C     Written by K. Hald, Feb. 2001.
C     Based on cc3_sort1 written by
C     Henrik Koch and Alfredo Sanchez.       28-May-1995
C
C     Sort virtual integrals for perturbative triples.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccinftap.h>
#include <ccsdsym.h>
!
      INTEGER LWORK, ISYINT, MAXCK, ISYMCK, ISYMD, ISYCKB, LENMIN
      INTEGER NDISTR, NBATCH, KSCR1, KSCR2, KSCR3, KSCR4, KSCR5, KSCR6
      INTEGER KEND1, NUMD, IBATCH, ID1, LENGTH, IOFF, ISYMB, ISYCKD
      INTEGER ID, ISYMK, ISYMC, ISYMBK, NTOTBK, KOFF1, KOFF2
      INTEGER IOPT, LU3SRT2, LU3SRT3, LUDELD2, LUDELD3
C      
!
#if defined (SYS_CRAY)
      REAL WORK(LWORK), ONE, XMONE
      REAL DDOT, XNORM
#else
      DOUBLE PRECISION WORK(LWORK), ONE, XMONE
      DOUBLE PRECISION DDOT, XNORM
#endif
!
      PARAMETER (ONE = 1.0D0, XMONE = -1.0D0)
!
      CHARACTER*(*) FN3SRT2, FN3SRT3, FNDELD2, FNDELD3
!
      CALL QENTER('CC3_SORT3')
!
C---------------------------------
C     Sanity check of iopt.
C---------------------------------
C
      IF (IOPT .NE. 1 .AND. IOPT .NE. 2) THEN
         CALL QUIT('Wrong iopt in cc3_sort3')
      ENDIF
C
C-----------------------------------------
C     Start loop over symmetries of delta.
C-----------------------------------------
C
      MAXCK = 0
      DO ISYMCK = 1,NSYM
         IF (NT1AM(ISYMCK) .GT. MAXCK) MAXCK = NT1AM(ISYMCK)
      ENDDO
C
      DO ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .NE. 0) THEN
C
C--------------------------
C        Memory allocation.
C--------------------------
C
           ISYCKB = MULD2H(ISYMD,ISYINT)
C
           LENMIN = NCKATR(ISYCKB) + MAXCK
           NDISTR = MIN(LWORK/LENMIN,NBAS(ISYMD))
C
           IF (NDISTR .EQ. 0) THEN
              CALL QUIT('Insufficient work space in CC3_SORT3')
           ENDIF
C
           NBATCH = (NBAS(ISYMD) - 1)/NDISTR + 1
C
           KSCR1 = 1
           KSCR2 = KSCR1 + NCKATR(ISYCKB)*NDISTR
           KSCR3 = KSCR2 + MAXCK*NDISTR
           KSCR4 = KSCR3 + NCKATR(ISYCKB)*NDISTR
           KSCR5 = KSCR4 + MAXCK*NDISTR
           KSCR6 = KSCR5 + NCKATR(ISYCKB)*NDISTR
           KEND1 = KSCR6 + MAXCK*NDISTR
C
           DO IBATCH = 1,NBATCH
C
              NUMD = NDISTR
              IF (IBATCH .EQ. NBATCH) THEN
                 NUMD = NBAS(ISYMD) - NDISTR*(NBATCH - 1)
              ENDIF
C
              ID1 = NDISTR*(IBATCH - 1) + 1
C
C--------------------------
C           Read integrals.
C--------------------------
C
              LENGTH = NCKATR(ISYCKB)*NUMD
C
              IOFF = ICKDAO(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(ID1 - 1) + 1
C
              IF (LENGTH .GT. 0) THEN
                 CALL GETWA2(LU3SRT2,FN3SRT2,WORK(KSCR1),IOFF,LENGTH)
                 CALL GETWA2(LU3SRT3,FN3SRT3,WORK(KSCR3),IOFF,LENGTH)
              ENDIF
C
C----------------------------------------------------
C            Copy g(c,k,b-bar,del) to scr5
C----------------------------------------------------
C
              CALL DCOPY(LENGTH,WORK(KSCR3),1,WORK(KSCR5),1)
C
C--------------------------------------------------------------------
C        For iopt = 1.
C            Sort g(c-bar,k,b,del) to g(b,k,del,c-bar)
C            Sort g(c-bar,k,b,del) + g(c,k-bar,b,del) to g(c,k,del,b)
C            Sort g(c,k,b-bar,del) to g(c,k,del,b)
C        For iopt = 2.
C            Sort g(c-bar,k,b,del) + g(c,k-bar,b,del) to g(b,k,del,c)
C            Sort g(c,k,b-bar,del) to g(b,k,del,c)
C--------------------------------------------------------------------
C
              DO ISYMB = 1,NSYM
C
                 ISYMCK = MULD2H(ISYCKB,ISYMB)
                 ISYCKD = MULD2H(ISYMCK,ISYMD)
C
                 DO B = 1,NVIR(ISYMB)
C
                    DO I = 1,NUMD
C
                       ID = ID1 + I - 1
C
                       DO ISYMK = 1,NSYM
C
                          ISYMC  = MULD2H(ISYMCK,ISYMK)
                          ISYMBK = MULD2H(ISYMB,ISYMK)
C
                          NTOTBK = MAX(NT1AM(ISYMBK),1)
C
                          IF (IOPT .EQ. 1) THEN
                             DO K = 1,NRHF(ISYMK)
C
                                KOFF1 = KSCR5
     *                                + NCKATR(ISYCKB)*(I - 1)
     *                                + ICKATR(ISYMBK,ISYMC)
     *                                + IT1AM(ISYMB,ISYMK)
     *                                + NVIR(ISYMB)*(K - 1) + B - 1
C
                                KOFF2 = KSCR6
     *                                + NT1AM(ISYMCK)*(I - 1)
     *                                + IT1AM(ISYMC,ISYMK)
     *                                + NVIR(ISYMC)*(K - 1)
C
                                CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),
     *                                     NTOTBK,WORK(KOFF2),1)
C
                             ENDDO   ! K loop
C
                          ELSE IF (IOPT .EQ. 2) THEN
C
                             DO K = 1,NRHF(ISYMK)
C
                                KOFF1 = KSCR1
     *                                + NCKATR(ISYCKB)*(I - 1)
     *                                + ICKATR(ISYMBK,ISYMC)
     *                                + IT1AM(ISYMB,ISYMK)
     *                                + NVIR(ISYMB)*(K - 1) + B - 1
C
                                KOFF2 = KSCR2
     *                                + NT1AM(ISYMCK)*(I - 1)
     *                                + IT1AM(ISYMC,ISYMK)
     *                                + NVIR(ISYMC)*(K - 1)
C
                                CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),
     *                                     NTOTBK,WORK(KOFF2),1)
C
                                KOFF1 = KSCR3
     *                                + NCKATR(ISYCKB)*(I - 1)
     *                                + ICKATR(ISYMBK,ISYMC)
     *                                + IT1AM(ISYMB,ISYMK)
     *                                + NVIR(ISYMB)*(K - 1) + B - 1
C
                                KOFF2 = KSCR4
     *                                + NT1AM(ISYMCK)*(I - 1)
     *                                + IT1AM(ISYMC,ISYMK)
     *                                + NVIR(ISYMC)*(K - 1)
C
                                CALL DCOPY(NVIR(ISYMC),WORK(KOFF1),
     *                                     NTOTBK,WORK(KOFF2),1)
C
                             ENDDO
C
                          ENDIF
C
                       ENDDO  ! ISYMK loop
                       
C
                       IF (IOPT .EQ. 1) THEN
                          KOFF1 = KSCR1
     *                          + NCKATR(ISYCKB)*(I - 1)
     *                          + ICKATR(ISYMCK,ISYMB)
     *                          + NT1AM(ISYMCK)*(B - 1)
                          KOFF2 = KSCR2
     *                          + NT1AM(ISYMCK)*(I - 1)
C
                          CALL DCOPY(NT1AM(ISYMCK),WORK(KOFF1),1,
     *                            WORK(KOFF2),1)
C
                          KOFF1 = KSCR3
     *                          + NCKATR(ISYCKB)*(I - 1)
     *                          + ICKATR(ISYMCK,ISYMB)
     *                          + NT1AM(ISYMCK)*(B - 1)
                          KOFF2 = KSCR4
     *                          + NT1AM(ISYMCK)*(I - 1)
C
                          CALL DCOPY(NT1AM(ISYMCK),WORK(KOFF1),1,
     *                               WORK(KOFF2),1)
                       ENDIF
C
                    ENDDO  ! I loop
C
C
C----------------------------------------
C       Add the contributions :
C       g^(1) = g(c-bar,k,b,del) + g(c,k-bar,b,del) - g(b,k,c-bar,del)
C       g^(2) = g(c-bar,k,b,del) + g(c,k-bar,b,del) - g(c,k,b-bar,del)
C----------------------------------------
C
                    CALL DAXPY(NT1AM(ISYMCK)*NUMD,XMONE,WORK(KSCR2),1,
     *                         WORK(KSCR4),1)
C
                    IF (IOPT .EQ. 1) THEN
                       CALL DAXPY(NT1AM(ISYMCK)*NUMD,XMONE,WORK(KSCR2),
     *                            1,WORK(KSCR6),1)
                    endif
C
C----------------------------------------
C                 Write sorted integrals.
C----------------------------------------
C
                    LENGTH = NT1AM(ISYMCK)*NUMD
C
                    IF (LENGTH .GT. 0) THEN
C
                       IOFF = ICKAD(ISYCKD,ISYMB)
     *                      + NCKA(ISYCKD)*(B - 1)
     *                      + ICKA(ISYMCK,ISYMD)
     *                      + NT1AM(ISYMCK)*(ID1 - 1) + 1
C
                          IF (IOPT .EQ. 1) THEN
C                     THIS IS THE OUTPUT FOR G1
                             CALL PUTWA2(LUDELD3,FNDELD3,WORK(KSCR6),
     *                                   IOFF,LENGTH)
C                     THIS IS THE OUTPUT FOR G2
                             CALL PUTWA2(LUDELD2,FNDELD2,WORK(KSCR4),
     *                                   IOFF,LENGTH)
                          ELSE IF (IOPT .EQ. 2) THEN
C                     THIS IS THE OUTPUT FOR THE SPECIAL G2
                             CALL PUTWA2(LUDELD2,FNDELD2,WORK(KSCR4),
     *                                   IOFF,LENGTH)
                          ENDIF
                    ENDIF
C
                 ENDDO  ! B loop
              ENDDO  ! ISYMB loop
C
           ENDDO  !IBATCH loop
         ENDIF
      ENDDO  ! ISYMD loop
C
      CALL QEXIT('CC3_SORT3')
C
      RETURN
      END
C  /* Deck cc3_onel33 */
      SUBROUTINE CC3_ONEL33(OMEGA1,OMEGA2P,OMEGA2M,RMAT1,RMAT2,RMAT3,
     *                      RMAT4,FOCKAK,SMAT,TMAT,ISYMIM,XIAJB,XIAJB2,
     *                      ISYINT,INDSQ,LENSQ,WORK,LWORK,ISYMIB,IB,
     *                      ISYMID,ID,IOPT)
C
C     K. Hald, Spring 2001
C
C     Based on cc3_onel by :
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C     Ove Christiansen 9-1-1996
C
C     Calculate Omega1 and Fock contibution to Omega2.
C
C     General symmetry: ISYMIM is symmetry of SMAT and TMAT 
C                       intermdiates.(incl isymd,isymb)
C                       ISYINT is symmetry of FOCKAK and XIAJB
C                       ISYRES = ISYMIM*ISYINT
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMIM, ISYINT, LENSQ, LWORK, ISYMIB, IB, ISYMID, ID
      INTEGER ISYRES, ISYMB, ISYMC, ISYMK, ISYMBC, ISYAIJ
      INTEGER ISYMCK, LENGTH, NCK, NTOAIJ, NTOTC, ISYMI
      INTEGER ISYAKJ, ISYMJ, ISYMBJ, ISYMAK, NBJ, NAK, NAKBJ, NAKJ
      INTEGER NTOAKJ, NBK, ISYMBK, NTOTB, ISYMKJ, ISYMAI, NCI, NIJ
      INTEGER NCIBJ, NTOTIJ, NTOTAK, ISYMIJ, JSAIKJ, KOFF1, KOFF2
      INTEGER NKJ, NCKBJ, NTOTAI, JSAKIJ, ISYMCI
      INTEGER ISYMLJ, ISYML, ISYMCL, NCL, NCLBJ, NLJ
      INTEGER INDEX, INDSQ(LENSQ,6), IOPT
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2P(*), OMEGA2M(*), RMAT1(*)
      REAL RMAT2(*), RMAT3(*), RMAT4(*), FOCKAK(*), SMAT(*)
      REAL TMAT(*), XIAJB(*), XIAJB2(*)
      REAL WORK(LWORK), ZERO, ONE, TWO, HALF
      REAL XMTWO, XMONE, FACT, XDOT, DDOT
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2P(*), OMEGA2M(*), RMAT1(*)
      DOUBLE PRECISION RMAT2(*), RMAT3(*), RMAT4(*), FOCKAK(*),SMAT(*)
      DOUBLE PRECISION TMAT(*), XIAJB(*), XIAJB2(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE, TWO, HALF
      DOUBLE PRECISION XMTWO, XMONE, FACT, XDOT, DDOT
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
      PARAMETER (XMTWO = -2.0D0, XMONE = -1.0D0)
C
      LOGICAL LDEBUG
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_ONEL33')
C
      LDEBUG = .FALSE.
C
C----------------------
C     General setup.
C----------------------
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      IF (IOPT. EQ. 2) THEN
        B = IB
        C = ID
        ISYMB = ISYMIB
        ISYMC = ISYMID
        FACT = ONE
      ELSE
        C = IB
        B = ID
        ISYMC = ISYMIB
        ISYMB = ISYMID
        FACT = XMONE
      ENDIF
C
C----------------------------------------
C     First contribution to Omega2P.
C----------------------------------------
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
         CALL QUIT('Not enough core in CCSDT_ONEL33')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) = FACT*(TWO*SMAT(I) + SMAT(INDSQ(I,4)))
      ENDDO
C
      NCK = IT1AM(ISYMC,ISYMK) + C
C
      KOFF1 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
      NTOAIJ = MAX(NCKI(ISYAIJ),1)
      NTOTC  = MAX(NVIR(ISYMC),1)
C
      IF (IOPT .EQ. 1) THEN
         CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *              FOCKAK(NCK),NTOTC,ONE,RMAT1,1)
      ELSE
         CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *              FOCKAK(NCK),NTOTC,ONE,RMAT2,1)
      ENDIF
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         IF (IOPT .EQ. 1) THEN
           XDOT = DDOT(NCKI(ISYAIJ),RMAT1,1,RMAT1,1)
           WRITE(LUPRI,*) 'CC3_ONEL33 : NORM RMAT1 (1) = ',XDOT
         ELSE
           XDOT = DDOT(NCKI(ISYAIJ),RMAT2,1,RMAT2,1)
           WRITE(LUPRI,*) 'CC3_ONEL33 : NORM RMAT2 (1) = ',XDOT
         ENDIF
      ENDIF
C
C----------------------------------------
C     First contribution to Omega2M.
C----------------------------------------
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
         CALL QUIT('Not enough core in CCSDT_ONEL')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) =  TWO*FACT*SMAT(I)
      ENDDO
C
      NCK = IT1AM(ISYMC,ISYMK) + C
C
      KOFF1 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
      NTOAIJ = MAX(NCKI(ISYAIJ),1)
      NTOTC  = MAX(NVIR(ISYMC),1)
C
      IF (IOPT .EQ. 1) THEN
         CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *              FOCKAK(NCK),NTOTC,ONE,RMAT3,1)
      ELSE
         CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *              FOCKAK(NCK),NTOTC,ONE,RMAT4,1)
      ENDIF
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         IF (IOPT .EQ. 1) THEN
           XDOT = DDOT(NCKI(ISYAIJ),RMAT3,1,RMAT3,1)
           WRITE(LUPRI,*) 'CC3_ONEL33 : NORM RMAT3 (1) = ',XDOT
         ELSE
           XDOT = DDOT(NCKI(ISYAIJ),RMAT4,1,RMAT4,1)
           WRITE(LUPRI,*) 'CC3_ONEL33 : NORM RMAT4 (1) = ',XDOT
         ENDIF
      ENDIF
C
C----------------------------------
C     First contribution to Omega1.
C----------------------------------
C
      ISYMI  = MULD2H(ISYMC,ISYRES)
      ISYAKJ = MULD2H(ISYMB,ISYINT)
C
C-----------------------------------------
C      Start by calculating the correct
C      linear combination in TMAT.
C------------------------------------------
C
      DO I = 1,LENGTH
         TMAT(I) =  SMAT(INDSQ(I,3)) - TWO*SMAT(INDSQ(I,2))
     *           -  SMAT(INDSQ(I,4))
      ENDDO
C
      IF ((.NOT. CC3LR) .AND. (NRHF(ISYMI) .NE. 0)) THEN
C
         IF (LWORK .LT. NCKI(ISYAKJ)) THEN
            CALL QUIT('Not enough core in CCSDT_ONEL33')
         END IF
C
C        Construct M(ak,j) = g(ak,bj)
C        ---------------------------
C        Interchanged the a and b index here!!!!
C
         DO ISYMJ = 1,NSYM
C
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAK = MULD2H(ISYMJ,ISYAKJ)
C
            DO J = 1,NRHF(ISYMJ)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO NAK = 1,NT1AM(ISYMAK)
C
                  NAKBJ = IT2AM(ISYMAK,ISYMBJ) + INDEX(NAK,NBJ)
                  NAKJ  = ICKI(ISYMAK,ISYMJ)
     *                  + NT1AM(ISYMAK)*(J - 1) + NAK
C
                  WORK(NAKJ) = FACT*XMTWO*XIAJB(NAKBJ)
C
               ENDDO
            ENDDO
         ENDDO
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
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XDOT = DDOT(NT1AM(ISYMCK),OMEGA1,1,OMEGA1,1)
         WRITE(LUPRI,*) 'CC3_ONEL33 :OMEGA1 NORM (1) = ',XDOT
      ENDIF
C
C-------------------------------------
C     Second contribution to Omega1.
C-------------------------------------
C
      ISYMAI = MULD2H(ISYMIM,ISYINT)
      ISYAKJ = MULD2H(ISYMB,ISYINT)
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYMLJ = MULD2H(ISYINT,ISYMBC)
C
C-----------------------------------------
C      Start by calculating the correct
C      linear combination in TMAT.
C------------------------------------------
C
      DO I = 1,LENGTH
         TMAT(I) =  SMAT(INDSQ(I,4))
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
C     Construct M(l,j) = g(cl,bj)
C     ---------------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYML  = MULD2H(ISYMLJ,ISYMJ)
         ISYMCL = MULD2H(ISYMC,ISYML)
C
         DO J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            DO L = 1,NRHF(ISYML)
C
               NCL = IT1AM(ISYMC,ISYML) + NVIR(ISYMC)*(L - 1) + C
C
               NCLBJ = IT2AM(ISYMCL,ISYMBJ) + INDEX(NCL,NBJ)
               NLJ  = IMATIJ(ISYML,ISYMJ)
     *              + NRHF(ISYML)*(J-1) + L
C
               WORK(NLJ) = FACT*XMTWO*XIAJB(NCLBJ)
C
            ENDDO
         ENDDO
      ENDDO
C
      NTOTC  = MAX(NVIR(ISYMC),1)
      NTOAKJ = MAX(NCKI(ISYAKJ),1)
      NTOTAI = MAX(1,NT1AM(ISYMAI))
C
      KOFF1 = ISAIKL(ISYMAI,ISYMLJ) + 1
      KOFF2 = 1
C
      CALL DGEMV('N',NT1AM(ISYMAI),NMATIJ(ISYMLJ),ONE,TMAT(KOFF1),
     *           NTOTAI,WORK,1,ONE,OMEGA1(KOFF2),1)
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XDOT = DDOT(NT1AM(ISYMAI),OMEGA1,1,OMEGA1,1)
         WRITE(LUPRI,*) 'CC3_ONEL33 : OMEGA1 NORM (2) = ',XDOT
      ENDIF
C
C---------------------------------------
C     Second contribution to Omega2P.
C---------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAKIJ = MULD2H(ISYMBC,ISYMIM)
      ISYMIJ = MULD2H(ISYMBC,ISYRES)
      ISYMAK = MULD2H(JSAKIJ,ISYMIJ)
C
      LENGTH = NCKIJ(JSAKIJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CC3_ONEL33')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) = -FACT*(SMAT(I)+SMAT(INDSQ(I,4)))
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
                  NCIBJ = IT2AM(ISYMCI,ISYMBJ) + INDEX(NCI,NBJ)
C
                  IF (NCI .EQ. NBJ) THEN
C
                     OMEGA2P(NCIBJ) = ZERO
C
                  ELSE IF (NCI .GT. NBJ) THEN
C
                     OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
C
                  ELSE IF (NCI .LT. NBJ) THEN
C
                     OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
C
                  ENDIF
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
                  OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
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
                  OMEGA2P(NCIBJ) = OMEGA2P(NCIBJ) + WORK(NIJ)
C
  340          CONTINUE
C
            ENDIF
C
  310    CONTINUE
C
  300 CONTINUE
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XDOT = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
         WRITE(LUPRI,*) 'OMEGA2P NORM (END) = ',XDOT
      ENDIF
C
      CALL QEXIT('CC3_ONEL33')
C
      RETURN
      END
CC  /* DECK SUM_R3 */
CC      SUBROUTINE SUM_R3(SMAT,ISYMB,B,ISYMD,D,NCKIJ,R3SUM)
CC
CC     Sum up the R3 amplitudes.
CC     This routine should only be commented in when debugging
CC     the CC3 triplet triple amplitudes, since we statically allocate
CC     the r3sum array.
CC
CC     HOWEVER : PLEASE DO NOT REMOVE THIS ROUTINE
CC
CC     K. Hald, Spring 2001.
CC
C      IMPLICIT NONE
CC
C#include <priunit.h>
C#include <ccorb.h>
CC
C      INTEGER ISYMB, ISYMD, B, D, nckij, i
C      INTEGER KOFF1, KOFF2
CC
C#if defined (SYS_CRAY)
C      REAL SMAT(*), r3sum(nvirt,nvirt,nvirt*nrhft*nrhft*nrhft)
C      REAL one, xmone, fact
C#else
C      DOUBLE PRECISION SMAT(*), one, xmone, fact
C      DOUBLE PRECISION r3sum(nvirt,nvirt,nvirt*nrhft*nrhft*nrhft)
C#endif
CC
C      PARAMETER (one = 1.0D0, xmone = -1.0D0)
C      LOGICAL LDEBUG
CC
C      CALL QENTER('SUM_R3')
CC
C      LDEBUG = .FALSE.
CC
C      IF ((ISYMB .EQ. ISYMD) .AND. (B .EQ. D)) THEN
C          CALL QUIT('SUM_R3 CALLED WITH B = D')
C      ENDIF
CC
C      KOFF1 = 0
C      KOFF2 = 0
CC
C      DO I = 1, ISYMB-1
C        KOFF1 = KOFF1 + NVIR(I)
C      ENDDO
C      DO I = 1, ISYMD-1
C        KOFF2 = KOFF2 + NVIR(I)
C      ENDDO
CC
C      DO I = 1, NCKIJ
CC
C         IF (LDEBUG) THEN
C            IF (ABS(SMAT(I)) .GT. 1.0d-9) THEN
C              WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',I,' WITH SMAT = ',
C     *                           SMAT(I)
C              WRITE(LUPRI,*) 'KOFF1 = ',KOFF1,' KOFF2 = ',KOFF2
C            endif
C         ENDIF
CC
C         R3SUM(B+KOFF1,D+KOFF2,I) = R3SUM(B+KOFF1,D+KOFF2,I) + SMAT(I)
C         R3SUM(D+KOFF2,B+KOFF1,I) = R3SUM(D+KOFF2,B+KOFF1,I) - SMAT(I)
CC
C      ENDDO
CC
C      CALL QEXIT('SUM_R3')
CC
C      RETURN
C      END
CC  /* Deck print_r3 */
C      SUBROUTINE PRINT_R3(R3SUM,ISYAMP)
CC
CC     Print the R3 amplitudes.
CC     This routine should only be commented in when debugging
CC     the CC3 triplet triple amplitudes, since we statically allocate
CC     the r3sum array.
CC
CC     HOWEVER : PLEASE DO NOT REMOVE THIS ROUTINE
CC
CC     K. Hald, Spring 2001.
C#include <priunit.h>
C#include <ccorb.h>
C#include <ccsdsym.h>
CC
C      integer isyamp, koff1, isckij, koff2, isyma, isymb, isymi, isymab
CC
C#if defined (SYS_CRAY)
C      REAL r3sum(nvirt,nvirt,nvirt*nrhft*nrhft*nrhft)
C#else
C      DOUBLE PRECISION r3sum(nvirt,nvirt,nvirt*nrhft*nrhft*nrhft)
C#endif
CC
C      CALL QENTER('PRINT_R3')
CC
C      do isyma = 1, nsym
CC
C        koff1 = 0
C        do isymi = 1, isyma-1
C          koff1 = koff1 + nvir(isymi)
C        enddo
CC
C        do isymb = 1, nsym
CC
C          koff2 = 0
C          do isymi = 1, isymb-1
C            koff2 = koff2 + nvir(isymi)
C          enddo
CC
C          isymab = muld2h(isyma,isymb)
C          isckij = muld2h(isymab,isyamp)
CC
C          do a = 1, nvir(isyma)
C            do b = 1, nvir(isymb)
C              do i = 1, nckij(isckij)
C                 if (abs(r3sum(a+koff1,b+koff2,i)) .gt. 1.0D-9) then
C              write(lupri,*) 'R3(',a+koff1,',',b+koff2,',',i,') = ',
C     *                        r3sum(a+koff1,b+koff2,i)
C                 endif
C              enddo
C            enddo
C          enddo
C        enddo
C      enddo
CC
C      CALL QEXIT('PRINT_R3')
CC
C      RETURN
C      END
C  /* Deck cc3_sort2 */
      SUBROUTINE CC3_SORT2(WORK,LWORK,ISYINT,LU3SRT,FN3SRT,
     *                     LUDELD5,FNDELD5)
C
C     KH April 2001 based on cc3_sort1 by
C     Henrik Koch and Alfredo Sanchez.       28-May-1995
C
C     Sort virtual integrals for perturbative triples.
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccinftap.h>
#include <ccsdsym.h>
      DIMENSION WORK(LWORK)
C
      CALL QENTER('CC3_SORT2')
C
C-----------------------------------------
C     Start loop over symmetries of delta.
C-----------------------------------------
C
      MAXCK = 0
      DO 50 ISYMCK = 1,NSYM
         IF (NT1AM(ISYMCK) .GT. MAXCK) MAXCK = NT1AM(ISYMCK)
   50 CONTINUE
C
      DO 100 ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 100
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         ISYCKB = MULD2H(ISYMD,ISYINT)
C
         LENMIN = NCKATR(ISYCKB) + MAXCK
         NDISTR = MIN(LWORK/LENMIN,NBAS(ISYMD))
C
         IF (NDISTR .EQ. 0) THEN
            CALL QUIT('Insufficient work space in CC3_SORT2')
         ENDIF
C
         NBATCH = (NBAS(ISYMD) - 1)/NDISTR + 1
C
         KSCR1 = 1
         KSCR2 = KSCR1 + NCKATR(ISYCKB)*NDISTR
         KEND1 = KSCR2 + MAXCK*NDISTR
C
         DO 110 IBATCH = 1,NBATCH
C
            NUMD = NDISTR
            IF (IBATCH .EQ. NBATCH) THEN
               NUMD = NBAS(ISYMD) - NDISTR*(NBATCH - 1)
            ENDIF
C
            ID1 = NDISTR*(IBATCH - 1) + 1
C
C--------------------------
C           Read integrals.
C--------------------------
C
            LENGTH = NCKATR(ISYCKB)*NUMD
C
            IOFF = ICKDAO(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(ID1 - 1) + 1
C
            IF (LENGTH .GT. 0) THEN
               CALL GETWA2(LU3SRT,FN3SRT,WORK(KSCR1),IOFF,LENGTH)
            ENDIF
C
C-----------------------------------------------------
C           Sort integrals (bk,del,c) from (ck,b,del).
C-----------------------------------------------------
C
            DO 150 ISYMC = 1,NSYM
C
               ISYMBK = MULD2H(ISYCKB,ISYMC)
               ISYBKD = MULD2H(ISYMBK,ISYMD)
C
               DO 160 C = 1,NVIR(ISYMC)
C
                  DO 170 I = 1,NUMD
C
                     ID = ID1 + I - 1
C
                     DO 180 ISYMK = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMBK,ISYMK)
                        ISYMCK = MULD2H(ISYMC,ISYMK)
C
                        NTOTCK = MAX(NT1AM(ISYMCK),1)
C
                        DO 190 K = 1,NRHF(ISYMK)

C
                           KOFF1 = KSCR1
     *                           + NCKATR(ISYCKB)*(I - 1)
     *                           + ICKATR(ISYMCK,ISYMB)
     *                           + IT1AM(ISYMC,ISYMK)
     *                           + NVIR(ISYMC)*(K - 1) + C - 1
C
                           KOFF2 = KSCR2
     *                           + NT1AM(ISYMBK)*(I - 1)
     *                           + IT1AM(ISYMB,ISYMK)
     *                           + NVIR(ISYMB)*(K - 1)
C
                           CALL DCOPY(NVIR(ISYMB),WORK(KOFF1),NTOTCK,
     *                                WORK(KOFF2),1)
C
  190                   CONTINUE
  180                CONTINUE
  170             CONTINUE
C
C----------------------------------------
C                 Write sorted integrals.
C----------------------------------------
C
                  LENGTH = NT1AM(ISYMBK)*NUMD
C
                  IF (LENGTH .GT. 0) THEN
C
                     IOFF = ICKAD(ISYBKD,ISYMC)
     *                    + NCKA(ISYBKD)*(C - 1)
     *                    + ICKA(ISYMBK,ISYMD)
     *                    + NT1AM(ISYMBK)*(ID1 - 1) + 1
C
                     CALL PUTWA2(LUDELD5,FNDELD5,WORK(KSCR2),
     *                           IOFF,LENGTH)
                  ENDIF
C
  160          CONTINUE
  150       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_SORT2')
C
      RETURN
      END
C  /* Deck cc3_trocc3 */
      SUBROUTINE CC3_TROCC3(XINT,TRINT,TRINT2,XLAMDH,WORK,LWORK,ISYINT)
C
C     Adapted to triplet by Kasper Hald, April 2001.
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C
C     Debugged for ISYINT .NE. 1 Ove Christiansen 11-1-1996
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      INTEGER LWORK, ISYINT, KOFF1, KOFF2, KOFF3, ISYMD, ISYMK
      INTEGER ISYAIJ, NTOIAJ, NBASD, ISYMJ, ISYMAI, ISYAIK, ISYMI
      INTEGER ISYMA, ISYMJI, ISYJIK, NTOJIK, ISYMJK
C
#if defined (SYS_CRAY)
      REAL XINT(*),TRINT(*), TRINT2(*), XLAMDH(*),WORK(LWORK)
      REAL XTROC, DDOT, ZERO, ONE
#else
      DOUBLE PRECISION XINT(*),TRINT(*), TRINT2(*), XLAMDH(*)
      DOUBLE PRECISION WORK(LWORK), XTROC, DDOT, ZERO, ONE
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('CC3_TROCC3')
C
      IF (LWORK .LT. NTRAOC(ISYINT)) THEN
         CALL QUIT('Insufficient space in CC3_TROCC3')
      END IF
C
C     write out integral norm
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTOTOC(ISYINT),XINT,1,XINT,1)
         WRITE(LUPRI,*) 'In CC3_TROCC3: Norm of INT = ',XTROC
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
         WRITE(LUPRI,*) 'CC3_TROCC3: Norm of transformed INT = ',XTROC
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
            ISYMJK = MULD2H(ISYMJ,ISYMK)
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
                        KOFF3 = ISJIKA(ISYJIK,ISYMA)
     *                        + ISJIK(ISYMJK,ISYMI)
     *                        + NMATIJ(ISYMJK)*(I - 1)
     *                        + IMATIJ(ISYMJ,ISYMK)
     *                        + NRHF(ISYMJ)*(K - 1) + J
C
                        CALL DCOPY(NVIR(ISYMA),WORK(KOFF1),1,
     *                             TRINT(KOFF2),NTOJIK)
                        CALL DCOPY(NVIR(ISYMA),WORK(KOFF1),1,
     *                             TRINT2(KOFF3),NTOJIK)
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
         WRITE(LUPRI,*) 'In CC3_TROCC3: Norm of INT Resorted  = ',XTROC
      ENDIF
C
      CALL QEXIT('CC3_TROCC3')
C
      RETURN
      END
C  /* Deck ccsdt_srtoc3 */
      SUBROUTINE CCSDT_SRTOC3(TROCC,TROCC1,ISYINT,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Sort occupied integrals I(kl,j,c) as I'(cj,l,k)
C
C     Ove Christiansen 16-1-1996: ISYINT sym of I(kl,j,c)
C     KH April 2001: Adapted to triplet.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      INTEGER ISYINT, LWORK, KOFF1, KOFF2, ISYMC, ISYKLJ, ISYMJ
      INTEGER ISYCKL, ISYCJL, ISYMKL
      INTEGER ISYMCL, ISYML, ISYMK, ISYMCK
C
#if defined (SYS_CRAY)
      REAL TROCC(*),TROCC1(*), WORK(LWORK), ZERO, ONE, TWO
#else
      DOUBLE PRECISION TROCC(*),TROCC1(*), WORK(LWORK), ZERO, ONE, TWO
#endif

C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
C
      CALL QENTER('CCSDT_SRTOC3')
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
                     ISYCJL = MULD2H(ISYINT,ISYMK)
                     ISYMCL = MULD2H(ISYMC,ISYML)
C
                     DO 150 L = 1,NRHF(ISYML)
C
                        DO 160 K = 1,NRHF(ISYMK)
C
                           KOFF1 = ISAIKJ(ISYCJL,ISYMK)
     *                           + NCKI(ISYCJL)*(K - 1)
     *                           + ISAIK(ISYMCL,ISYMJ)
     *                           + NT1AM(ISYMCL)*(J - 1)
     *                           + IT1AM(ISYMC,ISYML)
     *                           + NVIR(ISYMC)*(L - 1) + C
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
      CALL QEXIT('CCSDT_SRTOC3')
C
      RETURN
      END
C  /* Deck cc33_t2sort */
      SUBROUTINE CC33_T2SORT(C2,ISYAMP,C2TRANS,ISYMB,B,ISYMC,C)
C
C     Written April 2001 by KH.
C
C     Sort the C2 amplitudes from C2(cl,bj) stored as (cljb)
C     to a sorting of (ljcb) for a given b and c.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYAMP, LWORK, ISYMB, ISYMC
      INTEGER ISYCLJ, ISYMLJ, ISYML, ISYMJ, ISYMCL, ISYMBJ
      integer KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL C2(*), C2TRANS(*)
#else
      DOUBLE PRECISION C2(*), C2TRANS(*)
#endif
C
      CALL QENTER('CC33_T2SORT')
C
C----------------------------------
C     Symmetry initialitation
C----------------------------------
C
      ISYCLJ  = MULD2H(ISYAMP,ISYMB)
      ISYMLJ  = MULD2H(ISYCLJ,ISYMC)
C
      DO ISYML = 1, NSYM
C
         ISYMJ  = MULD2H(ISYMLJ,ISYML)
         ISYMCL = MULD2H(ISYMC,ISYML)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
C----------------------------------
C     Sort.
C----------------------------------
C
         DO L = 1, NRHF(ISYML)
C
            DO J = 1, NRHF(ISYMJ)
C
               KOFF1 = IT2SP(ISYCLJ,ISYMB)
     *               + NCKI(ISYCLJ)*(B-1)
     *               + ISAIK(ISYMCL,ISYMJ)
     *               + NT1AM(ISYMCL)*(J-1)
     *               + IT1AM(ISYMC,ISYML)
     *               + NVIR(ISYMC)*(L-1) + C
C
               KOFF2 = IMATIJ(ISYML,ISYMJ)
     *               + NRHF(ISYML)*(J-1) + L
C
               C2TRANS(KOFF2) = C2(KOFF1)
C
            ENDDO
         ENDDO
      ENDDO
C
      CALL QEXIT('CC33_T2SORT')
C
      RETURN
      END
C  /* Deck cc3_conocc33 */
      SUBROUTINE CC3_CONOCC33(OMEGA2P,OMEGA2M,RMAT1,RMAT2,RMAT3,RMAT4,
     *                        SMAT,TMAT,ISYMIM,TROCC,TROCC1,ISYINT,
     *                        WORK,LWORK,INDSQ,LENSQ,ISYMIB,IB,ISYMID,
     *                        ID,IOPT)
C
C     Henrik Koch and Alfredo Sanchez.         Dec   1994
C
C     Set up combinations of S's and contract with integrals.
C
C     Ove Christiansen 9-1-1996:
C
C     General symmetry: ISYMIM is symmetry of SMAT.
C                       (including isymib*isymid)
C                       ISYINT is symmetry of integrals in 
C                       TROCC and TROCC1.
C                       ISYRES = ISYMIM*ISYINT
C
C     Changed to triplet. Kasper Hald          April 2001
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER LWORK, ISYMIM, ISYINT, LENSQ, ISYMIB, IB, ISYMID, ID
      INTEGER IOPT, KOFF1, KOFF2, KOFF3, KOFF5, KOFF6, LENGTH
      INTEGER ISYRES, INDEX
      INTEGER NAI, NBJ, NRHFI, NTOCKL, ISYCKL, ISYMI, JSCKLI, ISYMAB
      INTEGER ISYMA, NTOTKL, NTOTAI, ISYKLJ, ISYMKL, ISYMAI, ISYMBJ
      INTEGER ISYMJ, JSAIKL, ISYMBC, ISYMB, ISYMC, ISYAIJ
      INTEGER INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL OMEGA2P(*), OMEGA2M(*), RMAT1(*), RMAT2(*)
      REAL RMAT3(*), RMAT4(*)
      REAL SMAT(*), TMAT(*), TROCC(*), TROCC1(*)
      REAL WORK(LWORK), FACT, ZERO, ONE, TWO, XMONE
      REAL DDOT, XRMAT, XTMAT, XDOT
#else
      DOUBLE PRECISION OMEGA2P(*), OMEGA2M(*), RMAT1(*), RMAT2(*)
      DOUBLE PRECISION RMAT3(*), RMAT4(*)
      DOUBLE PRECISION SMAT(*), TMAT(*), TROCC(*), TROCC1(*)
      DOUBLE PRECISION WORK(LWORK), FACT, ZERO, ONE, TWO, XMONE
      DOUBLE PRECISION DDOT, XRMAT, XTMAT, XDOT
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, XMONE = -1.0D0)
C
      LOGICAL LDEBUG
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CONOCC33')
C
      LDEBUG = .FALSE.
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
C------------------------
C     Sanity check.
C------------------------
C
      IF (LWORK .LT. LENSQ) THEN
         CALL QUIT('Insufficient core in CONOCC33')
      ENDIF
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2))
     *      CALL QUIT('WRONG IOPT IN CC3_CONOCC33')
C
C---------------------------------
C     First occupied R2+ term.
C---------------------------------
C
      IF (IOPT .EQ. 1) THEN
         C = ID
         B = IB
         ISYMC = ISYMID
         ISYMB = ISYMIB
         FACT = ONE
      ELSE
         C = IB
         B = ID
         ISYMC = ISYMIB
         ISYMB = ISYMID
         FACT = XMONE
      ENDIF
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
      DO I = 1,LENGTH
         TMAT(I) = -FACT*TWO*(SMAT(INDSQ(I,1))-SMAT(I)+SMAT(INDSQ(I,5)))
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
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
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
         IF (IOPT .EQ. 1) THEN
            CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *                 -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *                 ONE,RMAT2(KOFF3),NTOTAI)
         ELSE
            CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *                 -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *                 ONE,RMAT1(KOFF3),NTOTAI)
         ENDIF
C
  200 CONTINUE
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         ISYAIJ = MULD2H(ISYRES,ISYMB)
         IF (IOPT .EQ. 1) THEN
            XRMAT = DDOT(NCKI(ISYAIJ),RMAT2,1,RMAT2,1)
            WRITE(LUPRI,*) 'CC3_CONOCC33: RMAT2 norm (1) =  ',XRMAT
         ELSE
            XRMAT = DDOT(NCKI(ISYAIJ),RMAT1,1,RMAT1,1)
            WRITE(LUPRI,*) 'CC3_CONOCC33: RMAT1 norm (1) =  ',XRMAT
         ENDIF
      ENDIF
C
C---------------------------------
C     First occupied R2- term.
C---------------------------------
C
      IF (IOPT .EQ. 1) THEN
         C = ID
         B = IB
         ISYMC = ISYMID
         ISYMB = ISYMIB
         FACT = ONE
      ELSE
         C = IB
         B = ID
         ISYMC = ISYMIB
         ISYMB = ISYMID
         FACT = XMONE
      ENDIF
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
      DO I = 1,LENGTH
         TMAT(I) = FACT*TWO*SMAT(I)
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
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: 1. Norm of TMAT = ',XTMAT
      ENDIF
C
C-----------------------
C     First contraction.
C-----------------------
C
      DO ISYMJ = 1,NSYM
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
         IF (IOPT .EQ. 1) THEN
            CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *                 -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *                 ONE,RMAT4(KOFF3),NTOTAI)
         ELSE
            CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *                 -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *                 ONE,RMAT3(KOFF3),NTOTAI)
         ENDIF
C
      ENDDO ! ISYMJ
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         isyaij = muld2h(isyres,isymb)
         IF (IOPT .EQ. 1) THEN
            XRMAT = DDOT(NCKI(ISYAIJ),RMAT4,1,RMAT4,1)
            WRITE(LUPRI,*) 'CC3_CONOCC33: Norm of RMAT4 (1)= ',XRMAT
         ELSE
            XRMAT = DDOT(NCKI(ISYAIJ),RMAT3,1,RMAT3,1)
            WRITE(LUPRI,*) 'CC3_CONOCC33: Norm of RMAT3 (1)= ',XRMAT
         ENDIF
      ENDIF
C
C---------------------------------
C     Second occupied R2- term.
C---------------------------------
C
      IF (IOPT .EQ. 1) THEN
         B     = ID
         C     = IB
         ISYMB = ISYMID
         ISYMC = ISYMIB
         FACT  = ONE
      ELSE
         B     = IB
         C     = ID
         ISYMB = ISYMIB
         ISYMC = ISYMID
         FACT  = XMONE
      ENDIF
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
      DO I = 1,LENGTH
          TMAT(I) = FACT*TWO*SMAT(INDSQ(I,2))
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
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC33: 2. Norm of TMAT = ',XTMAT
      ENDIF
C
C------------------------
C     Second contraction.
C------------------------
C
      DO ISYMJ = 1,NSYM
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
     *              ZERO,WORK(KOFF3),NTOTAI)
C
      ENDDO
C
C----------------------------------------
C     Sort result into result vector.
C----------------------------------------
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
      IF (IOPT .EQ. 1) THEN
         CALL CC3_SORTIJ(WORK,RMAT3,ISYAIJ)
      ELSE
         CALL CC3_SORTIJ(WORK,RMAT4,ISYAIJ)
      ENDIF
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         IF (IOPT .EQ. 1) THEN
            XRMAT = DDOT(NCKI(ISYAIJ),RMAT3,1,RMAT3,1)
            WRITE(LUPRI,*) 'CC3_CONOCC33: Norm of RMAT3 (1)= ',XRMAT
         ELSE
            XRMAT = DDOT(NCKI(ISYAIJ),RMAT4,1,RMAT4,1)
            WRITE(LUPRI,*) 'CC3_CONOCC33: Norm of RMAT4 (1)= ',XRMAT
         ENDIF
      ENDIF
C
C--------------------------------
C     Second occupied R2+ term.
C--------------------------------
C
      IF (IOPT .EQ. 1) THEN
         A     = ID
         B     = IB
         ISYMA = ISYMID
         ISYMB = ISYMIB
         FACT  = ONE
      ELSE
         A     = IB
         B     = ID
         ISYMA = ISYMIB
         ISYMB = ISYMID
         FACT  = XMONE
      ENDIF
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
      DO I = 1,LENGTH
         TMAT(I) = FACT*(TWO*SMAT(INDSQ(I,1))-SMAT(I)-SMAT(INDSQ(I,4)))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC33: 3. Norm of TMAT = ',XTMAT
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
            CALL QUIT('Insufficient memory in CC3_CONOCC33')
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
                  IF (NAI .NE. NBJ) THEN
                     OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - WORK(KOFF6)
                  ENDIF
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
                  OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - WORK(KOFF6)
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
                  OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - WORK(KOFF6)
C
  640          CONTINUE
C
            ENDIF
C
  610    CONTINUE
C
  600 CONTINUE
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
           XDOT = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
           WRITE(LUPRI,*) 'CC3_CONOCC33 : Norm OMEGA2P (1) = ',XDOT
      ENDIF
C
      CALL QEXIT('CC3_CONOCC33')
C
      RETURN
      END
C  /* Deck cc3_sortij */
      SUBROUTINE CC3_SORTIJ(TEMPRMAT,RESRMAT,ISYAIJ)
C
C     Written by Kasper Hald, April 2001.
C
C     Sort the cc3 result vector ai,j to aj,i
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYAIJ, ISYMJ, ISYMAI, ISYMA, ISYMI, ISYMAJ, NAI
      INTEGER NAJ, KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL TEMPRMAT(*), RESRMAT(*)
#else
      DOUBLE PRECISION TEMPRMAT(*), RESRMAT(*)
#endif
C
      CALL QENTER('CC3_SORTIJ')
C
      DO ISYMJ = 1, NSYM
C
         ISYMAI = MULD2H(ISYAIJ,ISYMJ)
C
         DO ISYMA = 1, NSYM
C
            ISYMI  = MULD2H(ISYMAI,ISYMA)
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
            DO A = 1, NVIR(ISYMA)
C
               DO J = 1, NRHF(ISYMJ)
C
                  NAJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J-1) + A
C
                  DO I = 1, NRHF(ISYMI)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
                     KOFF1 = ISAIK(ISYMAI,ISYMJ)
     *                     + NT1AM(ISYMAI)*(J-1) + NAI
                     KOFF2 = ISAIK(ISYMAJ,ISYMI)
     *                     + NT1AM(ISYMAJ)*(I-1) + NAJ
C
                     RESRMAT(KOFF2) = RESRMAT(KOFF2) + TEMPRMAT(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL QEXIT('CC3_SORTIJ')
C
      RETURN
      END
C  /* Deck cc3_virtint */
      SUBROUTINE CC3_VIRTINT(XLAMDP,WORK,LWRK,ISINT1,ISINT2,
     *                       LU3VI,FN3VI,LU3VI2,FN3VI2,LU3VI4,FN3VI4,
     *                       LU3VI3,FN3VI3)
C
C     Written by Kasper Hald, April 2001.
C
C     Transform the integrals from (ai|del c) to (ai|bc)
C     and store them on file.
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccinftap.h>
C
      INTEGER LWRK, ISINT1, ISINT2, ISYMD, ISYCKB, ISCKB1
      INTEGER ISCKB2, KENDLO, LWRKLO, IOFF, KTRVI, KINTVI
      INTEGER LU3VI, LU3VI2, LU3VI4, LU3VI3
C
#if defined (SYS_CRAY)
      REAL XLAMDP(*), WORK(LWRK)
#else
      DOUBLE PRECISION XLAMDP(*), WORK(LWRK)
#endif
C
      CHARACTER*(*) FN3VI, FN3VI2, FN3VI4, FN3VI3
C
      CALL QENTER('CC3_VIRTINT')
C
      DO ISYMD = 1, NSYM
C
         ISYCKB = MULD2H(ISYMD,ISYMOP)
         ISCKB1 = MULD2H(ISYMD,ISINT1)
         ISCKB2 = MULD2H(ISYMD,ISINT2)
C
         KTRVI = 1
         KINTVI = KTRVI  + NCKATR(ISCKB1)
         KENDLO = KINTVI + MAX(NCKA(ISYMD),NCKA(ISCKB2))
         LWRKLO = LWRK - KENDLO
C
         IF (LWRKLO .LE. 0) THEN
            CALL QUIT('OUT OF MEMORY IN CC3_VIRTINT')
         ENDIF
C
         DO D = 1, NVIR(ISYMD)
C
C-----------------------------------------------
C     Read and transform the first integral.
C     Afterwards write to file.
C-----------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D-1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),XLAMDP,
     *                       ISYMD,D,ISYMOP,WORK(KENDLO),LWRKLO)
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
              IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D-1) + 1
              CALL PUTWA2(LU3VI3,FN3VI3,WORK(KTRVI),IOFF,NCKATR(ISCKB1))
            ENDIF
C
C
C-----------------------------------------------
C     Read and transform the second integral.
C     Afterwards write to file.
C-----------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),XLAMDP,
     *                       ISYMD,D,ISYMOP,WORK(KENDLO),LWRKLO)
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
              IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D-1) + 1
              CALL PUTWA2(LU3VI4,FN3VI4,WORK(KTRVI),IOFF,NCKATR(ISCKB1))
            ENDIF
C
         ENDDO
      ENDDO
C
      CALL QEXIT('CC3_VIRTINT')
C
      RETURN
      END
C  /* Deck cc3_convir33 */
      SUBROUTINE CC3_CONVIR33(OMEGA2P,RES2P,OMEGA2M,RES2M,RMAT1,RMAT2,
     *                        RMAT3,RMAT4,SMAT,TMAT,ISYMIM,TRVIR,TRVIR1,
     *                        TRVIR2,TRVIR3,ISYINT,WORK,LWORK,INDSQ,
     *                        LENSQ,ISYMB1,B1,ISYMD1,D1,IOPT,TIME1P,
     *                        TIME2P,TIME3P,TIME1M,TIME2M,TSORTP,TSORTM)
C
C     Written by Kasper Hald, April 2001.
C
C     Based on CC3_CONVIR by
C     Henrik Koch and Alfredo Sanchez. Dec 1994  &
C     Ove Christiansen 9-1-1996:
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <second.h>
C
      INTEGER LWORK, ISYMIM, ISYINT, LENSQ, ISYMB1, B1
      INTEGER ISYMD1, D1, IOPT, ISYRES, ISYMBD, ISCKIJ, LENGTH, ISYMJ 
      INTEGER KSCR1, KEND1, LWRK1, ISYMI, ISYMCK, ISYMA, NTOTCK, KOFF1
      INTEGER KOFF2, KOFF3, ISYMBJ, ISYMAI, ISYCKI, ISYMD, ISYMB, NVIRA
      INTEGER INDSQ(LENSQ,6), INDEX, KOFF4, KOFF5, NAI, NBJ, ISYAIB
      INTEGER NTOBIJ, NTOTL, ISYBIJ, ISYML, ISYMAL, ISYAEL, KRESP
      INTEGER ISBIJL, ISYMED, ISYME, ISYMAB, NTOTA, NMAXAI, KINT1
      INTEGER ISYCKA, ISYDOT, ISYMEL
C
#if defined (SYS_CRAY)
      REAL RES2P(*), RES2M(*)
      REAL OMEGA2P(*), OMEGA2M(*), RMAT1(*), RMAT2(*)
      REAL RMAT3(*), RMAT4(*), SMAT(*), TMAT(*), TRVIR(*), TRVIR1(*)
      REAL TRVIR2(*), TRVIR3(*), WORK(LWORK)
      REAL FACT, FACT2, ZERO, ONE, TWO, XMONE, XDOT, DDOT
      REAL TIME1P, TIME2P, TIME3P, TIME1M, TIME2M, DTIME
      REAL TSORTP, TSORTM
#else
      DOUBLE PRECISION RES2P(*), RES2M(*)
      DOUBLE PRECISION OMEGA2P(*), OMEGA2M(*), RMAT1(*), RMAT2(*)
      DOUBLE PRECISION RMAT3(*), RMAT4(*), SMAT(*), TMAT(*), TRVIR(*)
      DOUBLE PRECISION TRVIR1(*), TRVIR2(*), TRVIR3(*), WORK(LWORK)
      DOUBLE PRECISION FACT, FACT2, ZERO, ONE, TWO, XMONE, XDOT, DDOT
      DOUBLE PRECISION TIME1P, TIME2P, TIME3P, TIME1M, TIME2M, DTIME
      DOUBLE PRECISION TSORTP, TSORTM
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, XMONE = -1.0D0)
C
      LOGICAL LDEBUG
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CONVIR33')
C
      LDEBUG = .FALSE.
C
C------------------------
C     Sanity check.
C------------------------
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2)) THEN
         CALL QUIT('Wrong IOPT in CC3_CONVIR3')
      ENDIF
C
C------------------------------------------
C     First virtual contribution to (+).
C------------------------------------------
C
      DTIME = SECOND()
C
      IF (IOPT .EQ. 1) THEN
        ISYMB = ISYMB1
        B     = B1
        ISYMD = ISYMD1
        D     = D1
        FACT  = ONE
      ELSE
        ISYMB = ISYMD1
        B     = D1
        ISYMD = ISYMB1
        D     = B1
        FACT  = XMONE
      ENDIF
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
      LENGTH = NCKIJ(ISCKIJ)
C
      IF (LWORK .LT. NCKIJ(ISCKIJ)) THEN
         CALL QUIT('Insufficient core in CCSDT_CONVIR33 (1)')
      ENDIF
C
      DO I = 1,LENGTH
          TMAT(I) = -FACT*TWO*(SMAT(I)+SMAT(INDSQ(I,4)))
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
            CALL QUIT('Insufficient work space in CC3_CONVIR33 (1)')
         ENDIF
C
         DO 210 J = 1,NRHF(ISYMJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYMAI,ISYMI)
               ISYCKA = MULD2H(ISYMCK,ISYMA)
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
               IF (IOPT .EQ. 1) THEN
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NT1AM(ISYMCK),ONE,TRVIR1(KOFF1),NTOTCK,
     *                       TMAT(KOFF2),NTOTCK,ONE,RMAT2(KOFF3),NVIRA)
               ELSE
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NT1AM(ISYMCK),ONE,TRVIR3(KOFF1),NTOTCK,
     *                       TMAT(KOFF2),NTOTCK,ONE,RMAT1(KOFF3),NVIRA)
               ENDIF
C
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         ISYDOT = MULD2H(ISYMB,ISYRES)
         IF (IOPT .EQ. 1) THEN
           XDOT = DDOT(NCKI(ISYDOT),RMAT2,1,RMAT2,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm RMAT2 (1) = ',XDOT
         ELSE
           XDOT = DDOT(NCKI(ISYDOT),RMAT1,1,RMAT1,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm RMAT1 (1) = ',XDOT
         ENDIF
      ENDIF
C
      DTIME  = SECOND() - DTIME
      TIME1P = TIME1P + DTIME
C
C--------------------------------------------
C     Second virtual contribution to (+).
C--------------------------------------------
C
      DTIME = SECOND()
C
      IF (IOPT .EQ. 1) THEN
        ISYMB = ISYMB1
        B     = B1
        ISYMD = ISYMD1
        D     = D1
        FACT  = ONE
      ELSE
        ISYMB = ISYMD1
        B     = D1
        ISYMD = ISYMB1
        D     = B1
        FACT  = XMONE
      ENDIF
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
      LENGTH = NCKIJ(ISCKIJ)
C
      IF (LWORK .LT. NCKIJ(ISCKIJ)) THEN
         CALL QUIT('Insufficient core in CC3_CONVIR33 (2)')
      ENDIF
C
      DO I = 1,LENGTH
          TMAT(I) = FACT*(SMAT(I)+TWO*SMAT(INDSQ(I,4)))
      ENDDO
C
C---------------------------
C     Contract with (ac|kd).
C---------------------------
C
      DO ISYMJ = 1,NSYM
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
            CALL QUIT('Insufficient work space in CC3_CONVIR33 (2)')
         ENDIF
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYMAI,ISYMI)
               ISYCKA = MULD2H(ISYMCK,ISYMA)
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
               IF (IOPT .EQ. 1) THEN
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NT1AM(ISYMCK),ONE,TRVIR(KOFF1),NTOTCK,
     *                       TMAT(KOFF2),NTOTCK,ONE,RMAT2(KOFF3),NVIRA)
               ELSE
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NT1AM(ISYMCK),ONE,TRVIR2(KOFF1),NTOTCK,
     *                       TMAT(KOFF2),NTOTCK,ONE,RMAT1(KOFF3),NVIRA)
               ENDIF
C
            ENDDO  ! ISYMI
         ENDDO     ! J
      ENDDO        ! ISYMJ
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         ISYDOT = MULD2H(ISYMB,ISYRES)
         IF (IOPT .EQ. 1) THEN
           XDOT = DDOT(NCKI(ISYDOT),RMAT2,1,RMAT2,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm RMAT2 (2) = ',XDOT
         ELSE
           XDOT = DDOT(NCKI(ISYDOT),RMAT1,1,RMAT1,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm RMAT1 (2) = ',XDOT
         ENDIF
      ENDIF
C
      DTIME  = SECOND() - DTIME
      TIME2P = TIME2P + DTIME
C
C--------------------------------------------
C     First virtual contribution to (-).
C--------------------------------------------
C
      DTIME = SECOND()
C
      IF (IOPT .EQ. 1) THEN
        ISYMB = ISYMB1
        B     = B1
        ISYMD = ISYMD1
        D     = D1
        FACT  = ONE
      ELSE
        ISYMB = ISYMD1
        B     = D1
        ISYMD = ISYMB1
        D     = B1
        FACT  = XMONE
      ENDIF
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
      LENGTH = NCKIJ(ISCKIJ)
C
      IF (LWORK .LT. NCKIJ(ISCKIJ)) THEN
         CALL QUIT('Insufficient core in CC3_CONVIR33 (2)')
      ENDIF
C
      DO I = 1,LENGTH
          TMAT(I) = -FACT*TWO*SMAT(INDSQ(I,1))
      ENDDO
C
C---------------------------
C     Contract with (ac|kd).
C---------------------------
C
      DO ISYMJ = 1,NSYM
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
            CALL QUIT('Insufficient work space in CC3_CONVIR33 (2)')
         ENDIF
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYMAI,ISYMI)
               ISYCKA = MULD2H(ISYMCK,ISYMA)
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
               IF (IOPT .EQ. 1) THEN
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NT1AM(ISYMCK),ONE,TRVIR(KOFF1),NTOTCK,
     *                       TMAT(KOFF2),NTOTCK,ONE,RMAT4(KOFF3),NVIRA)
               ELSE
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NT1AM(ISYMCK),ONE,TRVIR2(KOFF1),NTOTCK,
     *                       TMAT(KOFF2),NTOTCK,ONE,RMAT3(KOFF3),NVIRA)
               ENDIF
C
            ENDDO  ! ISYMI
         ENDDO     ! J
      ENDDO        ! ISYMJ
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
         ISYDOT = MULD2H(ISYMB,ISYRES)
         IF (IOPT .EQ. 1) THEN
           XDOT = DDOT(NCKI(ISYDOT),RMAT4,1,RMAT4,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm RMAT4 (1) = ',XDOT
         ELSE
           XDOT = DDOT(NCKI(ISYDOT),RMAT3,1,RMAT3,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm RMAT3 (1) = ',XDOT
         ENDIF
      ENDIF
C
      DTIME  = SECOND() - DTIME
      TIME1M = TIME1M + DTIME
C
C---------------------------------------
C     Third virtual contribution to (+).
C---------------------------------------
C
      DTIME = SECOND()
C
      IF (IOPT .EQ. 1) THEN
        ISYME = ISYMB1
        E     = B1
        ISYMD = ISYMD1
        D     = D1
        FACT  = ONE
      ELSE
        ISYME = ISYMD1
        E     = D1
        ISYMD = ISYMB1
        D     = B1
        FACT  = XMONE
      ENDIF
C
      NMAXAI = 0
      DO ISYMI = 1, NSYM
         IF (NT1AM(ISYMI) .GT. NMAXAI) NMAXAI = NT1AM(ISYMI)
      ENDDO
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
      ISYMED = MULD2H(ISYME,ISYMD)
      ISYAEL = MULD2H(ISYINT,ISYMD)
      ISYMAL = MULD2H(ISYAEL,ISYME)
      ISBIJL = MULD2H(ISYMIM,ISYMED)
C
      KINT1 = 1
      KEND1 = KINT1 + NMAXAI
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LE. 0) THEN
         CALL QUIT('Insufficient work space in CC3_CONVIR33 (3)')
      ENDIF
C
      LENGTH = NCKIJ(ISBIJL)
      DO I = 1,LENGTH
         TMAT(I) = -FACT*(TWO*SMAT(I)+SMAT(INDSQ(I,4)))
      ENDDO
C
C--------------------------------
C     Calculate
C--------------------------------
C
      DO ISYML = 1, NSYM
C
         ISYMA  = MULD2H(ISYMAL,ISYML)
         ISYBIJ = MULD2H(ISBIJL,ISYML)
         ISYMEL = MULD2H(ISYME,ISYML)
C
C--------------------------------
C      Sort the integrals.
C--------------------------------
C
         DO L = 1, NRHF(ISYML)
            DO A = 1, NVIR(ISYMA)
               KOFF1 = KINT1 - 1
     *               + NVIR(ISYMA)*(L-1) + A
               KOFF2 = ICKATR(ISYMEL,ISYMA)
     *               + NT1AM(ISYMEL)*(A-1)
     *               + IT1AM(ISYME,ISYML)
     *               + NVIR(ISYME)*(L-1) + E
C
               IF (IOPT .EQ. 1) THEN
                  WORK(KOFF1) = TRVIR(KOFF2)
               ELSE
                  WORK(KOFF1) = TRVIR2(KOFF2)
               ENDIF
            ENDDO
         ENDDO
C
C-------------------------------------
C        Contract.
C-------------------------------------
C
C
         NTOTL  = MAX(1,NRHF(ISYML))
         NTOBIJ = MAX(1,NCKI(ISYBIJ))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         KOFF1  = ISAIKJ(ISYBIJ,ISYML)
     *          + 1
         KOFF2  = KINT1
         KOFF3 = 1
     *         + IT2SP(ISYBIJ,ISYMA)
C
         IF (IOPT .EQ. 1) THEN
            CALL DGEMM('N','T',NCKI(ISYBIJ),NVIR(ISYMA),NRHF(ISYML),
     *                 ONE,TMAT(KOFF1),NTOBIJ,WORK(KOFF2),NTOTA,ONE,
     *                 RES2P(KOFF3),NTOBIJ)
         ELSE
            CALL DGEMM('N','T',NCKI(ISYBIJ),NVIR(ISYMA),NRHF(ISYML),
     *                 ONE,TMAT(KOFF1),NTOBIJ,WORK(KOFF2),NTOTA,ONE,
     *                 RES2P(KOFF3),NTOBIJ)
         ENDIF
C
      ENDDO  !ISYML
C
      DTIME  = SECOND() - DTIME
      TIME3P = TIME3P + DTIME
C
C---------------------------------------
C     Sort result into result vector.
C---------------------------------------
C
      DTIME = SECOND()
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
           XDOT = DDOT(NT2AM(ISYRES),OMEGA2P,1,OMEGA2P,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm OMEGA2P (1) = ',XDOT
      ENDIF
C
      DTIME  = SECOND() - DTIME
      TSORTP = TSORTP + DTIME
C
C------------------------------------------
C     Second virtual contribution to (-).
C------------------------------------------
C
      DTIME = SECOND()
C
      IF (IOPT .EQ. 1) THEN
        ISYME = ISYMB1
        E     = B1
        ISYMD = ISYMD1
        D     = D1
        FACT  = ONE
      ELSE
        ISYME = ISYMD1
        E     = D1
        ISYMD = ISYMB1
        D     = B1
        FACT  = XMONE
      ENDIF
C
      NMAXAI = 0
      DO ISYMI = 1, NSYM
         IF (NT1AM(ISYMI) .GT. NMAXAI) NMAXAI = NT1AM(ISYMI)
      ENDDO
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
      ISYMED = MULD2H(ISYME,ISYMD)
      ISYAEL = MULD2H(ISYINT,ISYMD)
      ISYMAL = MULD2H(ISYAEL,ISYME)
      ISBIJL = MULD2H(ISYMIM,ISYMED)
C
      KINT1 = 1
      KEND1 = KINT1 + NMAXAI
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LE. 0) THEN
         CALL QUIT('Insufficient work space in CC3_CONVIR33 (3)')
      ENDIF
C
      CALL DZERO(WORK(KINT1),NMAXAI)
C
      LENGTH = NCKIJ(ISBIJL)
      DO I = 1,LENGTH
         TMAT(I) = -FACT*SMAT(INDSQ(I,3))
      ENDDO
C
C--------------------------------
C     Calculate
C--------------------------------
C
      DO ISYML = 1, NSYM
C
         ISYMA  = MULD2H(ISYMAL,ISYML)
         ISYBIJ = MULD2H(ISBIJL,ISYML)
         ISYMEL = MULD2H(ISYME,ISYML)
C
C--------------------------------
C      Sort the integrals.
C--------------------------------
C
         DO L = 1, NRHF(ISYML)
            DO A = 1, NVIR(ISYMA)
               KOFF1 = KINT1 - 1
     *               + NVIR(ISYMA)*(L-1) + A
               KOFF2 = ICKATR(ISYMEL,ISYMA)
     *               + NT1AM(ISYMEL)*(A-1)
     *               + IT1AM(ISYME,ISYML)
     *               + NVIR(ISYME)*(L-1) + E
C
               IF (IOPT .EQ. 1) THEN
                  WORK(KOFF1) = TRVIR(KOFF2)
               ELSE
                  WORK(KOFF1) = TRVIR2(KOFF2)
               ENDIF
            ENDDO
         ENDDO
C
C-------------------------------------
C        Contract.
C-------------------------------------
C
C
         NTOTL  = MAX(1,NRHF(ISYML))
         NTOBIJ = MAX(1,NCKI(ISYBIJ))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         KOFF1  = ISAIKJ(ISYBIJ,ISYML)
     *          + 1
         KOFF2  = KINT1
         KOFF3  = 1
     *          + IT2SP(ISYBIJ,ISYMA)
C
         IF (IOPT .EQ. 1) THEN
            CALL DGEMM('N','T',NCKI(ISYBIJ),NVIR(ISYMA),NRHF(ISYML),
     *                 ONE,TMAT(KOFF1),NTOBIJ,WORK(KOFF2),NTOTA,ONE,
     *                 RES2M(KOFF3),NTOBIJ)
         ELSE
            CALL DGEMM('N','T',NCKI(ISYBIJ),NVIR(ISYMA),NRHF(ISYML),
     *                 ONE,TMAT(KOFF1),NTOBIJ,WORK(KOFF2),NTOTA,ONE,
     *                 RES2M(KOFF3),NTOBIJ)
         ENDIF
C
      ENDDO  !ISYML
C
      DTIME  = SECOND() - DTIME
      TIME2M = TIME2M + DTIME
C
C---------------------------------------
C     Sort result into result vector.
C---------------------------------------
C
      DTIME = SECOND()
C
      IF (IPRINT .GT. 55 .OR. LDEBUG) THEN
           XDOT = DDOT(NT2AM(ISYRES),OMEGA2M,1,OMEGA2M,1)
           WRITE(LUPRI,*) 'CC3_CONVIR33 : Norm OMEGA2M (1) = ',XDOT
      ENDIF
C
      DTIME  = SECOND() - DTIME
      TSORTM = TSORTM + DTIME
C
      CALL QEXIT('CC3_CONVIR33')
C
      RETURN
      END
C  /* Deck cc3_sortminus */
      SUBROUTINE CC3_SORTMINUS(OMEGA2M,RES2M,ISYRES)
C
C     K. Hald, April 2001.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMA, ISYMB, ISYMAB, ISYMI, ISYMAI, ISYAIB, ISYMJ
      INTEGER ISYMBJ, ISYBIJ, NAI, NBJ, KOFF4, KOFF5, INDEX, ISYRES
C
#if defined (SYS_CRAY)
      REAL OMEGA2M(*), RES2M(*), ZERO, ONE, XMONE, FACT2
#else
      DOUBLE PRECISION OMEGA2M(*), RES2M(*), ZERO, ONE, XMONE, FACT2
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, XMONE = -1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_SORTMINUS')
C
C---------------------------------
C     Sort the (-) part.
C---------------------------------
C
      DO ISYMA = 1, NSYM
C
         DO ISYMB = 1, NSYM
C
            ISYMAB = MULD2H(ISYMA,ISYMB)
C
            DO ISYMI = 1, NSYM
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYAIB = MULD2H(ISYMAI,ISYMB)
               ISYMJ  = MULD2H(ISYAIB,ISYRES)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYBIJ = MULD2H(ISYMBJ,ISYMI)
C
               DO A = 1, NVIR(ISYMA)
               DO I = 1, NRHF(ISYMI)
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  DO B = 1, NVIR(ISYMB)
                  DO J = 1, NRHF(ISYMJ)
C
                     NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                     KOFF4 =  IT2SP(ISYBIJ,ISYMA)
     *                     + NCKI(ISYBIJ)*(A-1)
     *                     + ISAIK(ISYMBJ,ISYMI)
     *                     + NT1AM(ISYMBJ)*(I-1) + NBJ
C
                     IF (ISYMAI.EQ.ISYMBJ) THEN
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                       + INDEX(NAI,NBJ)
                        IF (NAI .EQ. NBJ) THEN
                           FACT2 = ZERO
                        ELSE IF (NAI .GT. NBJ) THEN
                           FACT2 = ONE
                        ELSE
                           FACT2 = XMONE
                        ENDIF
                     ELSE IF (ISYMAI .LT. ISYMBJ) THEN
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1) + NAI
                        FACT2 = XMONE
                     ELSE IF (ISYMBJ .LT. ISYMAI) THEN
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                        FACT2 = ONE
                     ENDIF
C
                     OMEGA2M(KOFF5) = OMEGA2M(KOFF5) -FACT2*RES2M(KOFF4)
C
                  ENDDO  ! J
                  ENDDO  ! B
               ENDDO     ! I
               ENDDO     ! A
            ENDDO        ! ISYMI
         ENDDO           ! ISYMB
      ENDDO              ! ISYMA
C
      CALL QEXIT('CC3_SORTMINUS')
C
      RETURN
      END
C  /* Deck cc3_sortplus */
      SUBROUTINE CC3_SORTPLUS(OMEGA2P,RES2P,ISYRES)
C
C     K. Hald, April 2001.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMA, ISYMB, ISYMAB, ISYMI, ISYMAI, ISYAIB, ISYMJ
      INTEGER ISYMBJ, ISYBIJ, NAI, NBJ, KOFF4, KOFF5, INDEX, ISYRES
C
#if defined (SYS_CRAY)
      REAL OMEGA2P(*), RES2P(*), ZERO
#else
      DOUBLE PRECISION OMEGA2P(*), RES2P(*), ZERO
#endif
C
      PARAMETER(ZERO = 0.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_SORTPLUS')
C
C--------------------------------------
C     Sort the (+) result.
C--------------------------------------
C
      DO ISYMA = 1, NSYM
C
         DO ISYMB = 1, NSYM
C
            ISYMAB = MULD2H(ISYMA,ISYMB)
C
            DO ISYMI = 1, NSYM
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYAIB = MULD2H(ISYMAI,ISYMB)
               ISYMJ  = MULD2H(ISYAIB,ISYRES)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYBIJ = MULD2H(ISYMBJ,ISYMI)
C
               DO A = 1, NVIR(ISYMA)
               DO I = 1, NRHF(ISYMI)
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  DO B = 1, NVIR(ISYMB)
                  DO J = 1, NRHF(ISYMJ)
C
                     NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                     KOFF4 = IT2SP(ISYBIJ,ISYMA)
     *                     + NCKI(ISYBIJ)*(A-1)
     *                     + ISAIK(ISYMBJ,ISYMI)
     *                     + NT1AM(ISYMBJ)*(I-1) + NBJ
C
                     IF (ISYMAI.EQ.ISYMBJ) THEN
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                       + INDEX(NAI,NBJ)
                        IF (NAI .EQ. NBJ) THEN
                           RES2P(KOFF4) = ZERO
                        ENDIF
                     ELSE IF (ISYMAI .LT. ISYMBJ) THEN
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1) + NAI
                     ELSE IF (ISYMBJ .LT. ISYMAI) THEN
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                     ENDIF
C
                     OMEGA2P(KOFF5) = OMEGA2P(KOFF5) - RES2P(KOFF4)
C
                  ENDDO  ! J
                  ENDDO  ! B
               ENDDO     ! I
               ENDDO     ! A
            ENDDO        ! ISYMI
         ENDDO           ! ISYMB
      ENDDO              ! ISYMA
C
      CALL QEXIT('CC3_SORTPLUS')
C
      RETURN
      END
