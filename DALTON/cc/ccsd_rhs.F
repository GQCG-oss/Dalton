!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
C  /* Deck ccrhsn */
      SUBROUTINE CCRHSN(OMEGA1,OMEGA2,T1AM,T2AM,WORK,LWORK,APROXR12)
C
C     Written by Henrik Koch 25-Sep-1993
C
C     Version 3.0
C
C     Purpose:
C
C     Calculation of the Coupled Cluster vector function using
C     AO-integrals directly from disk. 
C
C
C     NB! It is assumed that the vectors are allocated in the following
C     order:
C           T1AM(*), OMEGA1(*), OMEGA2(*), T2AM(*),  WORK(*).
C
C     some changes for CC2 with non-Hatree-Fock fields (NONHF=.true.)
C     to allow for finite difference also w.r.t. orbital coefficients
C     (i.e. the CMO vector), spring 2000, Ch. Haettig
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
#include "r12int.h"
#include "ccr12int.h"
#include "qm3.h"
!#include "qmmm.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (XMHALF = -0.5D0, XMONE= -1.0D0 )
      PARAMETER (ISYM0 = 1)
C
      LOGICAL FCKCON,CC1BSA,ETRAN,CC2R12,CCSDR12,LV,LVAJKL,LRES,
     &        DEBUGV,LVIJKL,LVABKL
      PARAMETER (DEBUGV = .FALSE.)
C
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION OMEGA1(*),OMEGA2(*),T1AM(*),T2AM(*),WORK(LWORK)
C
      CHARACTER CFIL*6,DFIL*6, FN3SRT*8, FNDELD*6, CDUMMY*8
      CHARACTER FNCKJD*6, FNDKBC*4, FNTOC*8, FN3VI*6, FN3VI2*8
      CHARACTER FNIADJ*8, FNIJDA*8, CPFIL*8, DPFIL*8
      CHARACTER*(*) APROXR12
      CHARACTER MODEL*10
C
      PARAMETER (FNIADJ = 'CCXIADJ0', FNIJDA = 'CCXIJDA0')
      PARAMETER (CPFIL  = 'CC_CPR12', DPFIL  = 'CC_DPR12')
C

      INTEGER IGLMRHS(8,8),IGLMVIS(8,8),NGLMDS(8),ICMO(8,8),NCMO(8),
     &        IMAIJM(8,8),NMAIJM(8),
     &        IMATIJM(8,8),NMATIJM(8),NGAMSM(8),IGAMSM(8,8),
     &        IRGIJS(8,8),NRGIJS(8),IR1BASM(8,8),NR1BASM(8),
     &        IR2BASM(8,8),NR2BASM,IR1XBASM(8,8),NR1XBASM(8),
     &        IR2XBASM(8,8),IMATF(8,8),NMATF(8),IMAKLM(8,8),NMAKLM(8)
      INTEGER NADP(8),IADP(8,8),NLAMDX(8),ILAMDX(8,8)
C
      CALL QENTER('CCRHSN')
C
      CC2R12  = CC2 .AND. LMULBS
      IF (LMULBS.AND. .NOT.(CC2R12 .OR. CCS .OR. CIS)) THEN
        CCSDR12 = .TRUE.
        IF (IANR12.EQ.2) CALL QUIT('CCSD(R12) only implemented for '//
     &                             'Ansaetze 1 and 3')
      ELSE
        CCSDR12 = .FALSE.
      END IF
      IF (LMULBS.AND.NONHF.AND.IANR12.NE.1) THEN
        CALL QUIT('CC-R12 with finite fields only implemented for '//
     &             'Ansatz 1')
      END IF
CTesT
C     CCSDR12 = .TRUE.
C     DUMPCD = .TRUE.
CTesT
C
C-----------------------------------------------------------
C     For energy calculation trial vector is totalsymmetric.
C-----------------------------------------------------------
C
      ISYMTR = 1
C
C-----------------------------------------
C     Save CC1B flag and if CC1A set true.
C-----------------------------------------
C
      CC1BSA = CC1B
      IF ( CC1A ) CC1B = .TRUE.
C
      IF ( IPRINT .GT. 10 ) THEN
C
         WRITE(LUPRI,*) ' In ccsd_rhs : '
         WRITE(LUPRI,*) ' CCSD, CC2: ',CCSD,CC2
         WRITE(LUPRI,*) ' CC1A, CC1B, CC3: ', CC1A, CC1B, CC3
C
      ENDIF
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
      IF (DEBUG) WRITE(LUPRI,*) 'DUMPCD = ',DUMPCD
      IF (DUMPCD) THEN
         CALL WOPEN2(LUC,CFIL,64,0)
         CALL WOPEN2(LUD,DFIL,64,0)
C
      END IF
C
      IF (CCSDT) THEN
C
         LU3SRT = -1
         LUCKJD = -1
         LUDELD = -1
         LUDKBC = -1
         LUTOC  = -1
         LU3VI  = -1
         LU3VI2 = -1
         FN3SRT = 'CC3_SORT'
         FNCKJD = 'CKJDEL'
         FNDELD = 'CKDELD'
         FNDKBC = 'DKBC'
         FNTOC  = 'CCSDT_OC'
         FN3VI  = 'CC3_VI'
         FN3VI2 = 'CC3_VI12'
C
         CALL WOPEN2(LU3SRT,FN3SRT,64,0)
         CALL WOPEN2(LUCKJD,FNCKJD,64,0)
         CALL WOPEN2(LUDELD,FNDELD,64,0)
         CALL WOPEN2(LUDKBC,FNDKBC,64,0)
         CALL WOPEN2(LUTOC,FNTOC,64,0)
         CALL WOPEN2(LU3VI,FN3VI,64,0)
         CALL WOPEN2(LU3VI2,FN3VI2,64,0)
C
      ENDIF

      IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
         LUIADJ = -1
         LUIJDA = -1
         CALL WOPEN2(LUIADJ,FNIADJ,64,0)
         CALL WOPEN2(LUIJDA,FNIJDA,64,0)
         LUCP = -1
         LUDP = -1
         CALL WOPEN2(LUCP,CPFIL,64,0)
         CALL WOPEN2(LUDP,DPFIL,64,0)
      END IF
C
C----------------------------------
C     Initialize timing parameters.
C----------------------------------
C
      TIMALL  = SECOND()
      TIMA    = 0.0D00
      TIMB    = 0.0D00
      TIMBF   = 0.0D00
      TIMC    = 0.0D00
      TIMD    = 0.0D00
      TIME    = 0.0D00
      TIMEP   = 0.0D00
      TIMF    = 0.0D00
      TIMFP   = 0.0D00
      TIMG    = 0.0D00
      TIMGP   = 0.0D00
      TIMH    = 0.0D00
      TIMI    = 0.0D00
      TIMJ    = 0.0D00
      TIMGAM  = 0.0D00
      TIMEI   = 0.0D00
      TIMLAM  = 0.0D00
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
      TIMRDAOR12 = 0.0D00
      TIMINTR12  = 0.0D00
C
C---------------------------
C     Check inconsistencies.
C---------------------------
C
      IF (NEWGAM) THEN
         IF ((.NOT. DUMPCD) .OR. (.NOT. OMEGOR)) THEN
            WRITE(LUPRI,*) 'NEWGAM requires both DUMPCD and OMEGOR'
            CALL QUIT('ERROR: NEWGAM inconsistency')
         END IF
      END IF
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      KLAMDP = 1
      KLAMIP = KLAMDP + NLAMDT
      IF (.NOT. DUMPCD) THEN
          KLAMDH = KLAMIP + NLAMDT
      ELSE
          KLAMDH = KLAMIP + 1
      END IF
      KDENSI = KLAMDH + NLAMDT
      KFOCK  = KDENSI + N2BAST
      KEMAT1 = KFOCK  + N2BST(ISYMOP)
      KEMAT2 = KEMAT1 + NEMAT1(ISYMOP)
      KGAMMA = KEMAT2 + NMATIJ(ISYMOP)
      IF (NEWGAM) THEN
         KEND1 = KGAMMA
      ELSE
         KEND1 = KGAMMA + NGAMMA(ISYMOP)
      END IF
      IF (CC2 .AND. NONHF) THEN
        KFCKHF = KEND1
        KEND1  = KFCKHF + N2BAST
      END IF
c
      IF (CCR12) THEN
         KVIJKL  = KEND1
         KEND1   = KVIJKL + NTR12SQ(1)
      END IF

      IF (CCR12) THEN
         CALL CC_R12OFFS23(IGLMRHS,IGLMVIS,NGLMDS,ICMO,NCMO,
     &                     IMAIJM,NMAIJM,IMAKLM,NMAKLM,
     &                     IMATIJM,NMATIJM,
     &                     IGAMSM,NGAMSM,IRGIJS,NRGIJS,
     &                     IR1BASM,NR1BASM,IR2BASM,NR2BASM,IR1XBASM,
     &                     NR1XBASM,IR2XBASM,IMATF,NMATF)
         KLAMDHS = KEND1 
         KLAMDPS = KLAMDHS + NGLMDS(1)
         KEND1   = KLAMDPS + NGLMDS(1)
 
         CALL LAMMATS(WORK(KLAMDPS),WORK(KLAMDHS),T1AM,
     &                1,.FALSE.,.FALSE.,
     &                NGLMDS,IGLMRHS,IGLMVIS,ICMO,WORK(KEND1),LWRK1)
      END IF

      IF (CCR12.AND..NOT.USEVABKL) THEN
         KVAJKL = KEND1
         KEND1  = KVAJKL + NVAJKL(1)
      END IF
C
      IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
         DO ISYM = 1, NSYM
           NLAMDX(ISYM) = 0
           NADP(ISYM) = 0
           DO ISYM2 = 1, NSYM
             ISYM1 = MULD2H(ISYM,ISYM2)
             ILAMDX(ISYM1,ISYM2) = NLAMDX(ISYM)
             NLAMDX(ISYM) = NLAMDX(ISYM) +
     &           (MBAS1(ISYM1)+MBAS2(ISYM1))*(NORB1(ISYM2)+NORB2(ISYM2))
             IADP(ISYM1,ISYM2) = NADP(ISYM)
             NADP(ISYM) = NADP(ISYM) + 
     &                    NVIR(ISYM1)*(MBAS1(ISYM2)+MBAS2(ISYM2))
           END DO
         END DO
 
         KFCKVAO = KEND1 
         KEND1   = KFCKVAO + NEMAT1(1)

         KE1PIM = KEND1
         KEND1  = KE1PIM + NADP(1)
      ELSE
         KE1PIM = KEND1
      END IF
C
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHSN')
      ENDIF
C
C------------------------------------
C     Save the CC amplitudes on disk.
C------------------------------------
C
      LURHS1 = -1
      CALL GPOPEN(LURHS1,'CCRHS1','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LURHS1)
      WRITE (LURHS1) (T1AM(I), I = 1,NT1AMX)
      WRITE (LURHS1) (T2AM(I), I = 1,NT2AMX)
C
C----------------------------------
C     Calculate the lamda matrices.
C----------------------------------
C
      TIMLAM  = SECOND()
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND1),LWRK1)
      TIMLAM  = SECOND() - TIMLAM
C
C-----------------------------------------
C     Calculate the inverse xlamdp matrix.
C-----------------------------------------
C
      IF (.NOT. DUMPCD)
     *   CALL CCSD_INVLDP(WORK(KLAMDP),WORK(KLAMIP),WORK(KEND1),LWRK1)
C
C-----------------------------------
C     initialize R12 vector function
C-----------------------------------
      IF (CCR12) CALL DZERO(WORK(KVIJKL),NTR12SQ(1))
      IF (CCR12.AND..NOT.USEVABKL) THEN
        IOPT = 2
        CALL CC_R12MKVAMKL0(WORK(KVIJKL),NTR12SQ(1),IOPT,WORK(KLAMDH),
     &                      1,WORK(KEND1),LWRK1)
        IF (RSPIM) THEN
          IOPT = 1
          CALL CC_R12MKVAMKL0(WORK(KVAJKL),NVAJKL(1),IOPT,WORK(KLAMDH),
     &                        1,WORK(KEND1),LWRK1)
        END IF
      END IF
      IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
        ! initialize the Fhat(a,del) matrix
        CALL DZERO(WORK(KFCKVAO),NEMAT1(1))
      END IF
C
C-------------------------------
C     Prepare the t2-amplitudes.
C-------------------------------
C
      CALL DCOPY(NT2AMX,T2AM,1,OMEGA2,1)
      CALL CC_T2SQ(OMEGA2,T2AM,ISYMTR)
C
C-----------------------------------------
C     Construct the transposed amplitudes.
C-----------------------------------------
C
      IF (CCSDT .OR. CCSDR12) THEN
         KEND1T = KEND1
         LWRK1T = LWRK1
      ENDIF
C
      IF ((.NOT. DIRECT) .AND. T2TCOR) THEN
C
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
C
      END IF
C
C-------------------------------
C     Initialize OMEGA1 & OMEGA2
C-------------------------------
C
      CALL DZERO(OMEGA1,NT1AM(ISYMOP))
      IF (.NOT. OMEGSQ) THEN
         IF (OMEGOR) THEN
            CALL DZERO(OMEGA2,2*NT2ORT(ISYMOP))
         ELSE
            CALL DZERO(OMEGA2,NT2AO(ISYMOP))
         ENDIF
      ELSE
         CALL DZERO(OMEGA2,NT2AOS(ISYMOP))
      ENDIF
C
C-------------------------------------
C     Initialize GAMMA, EMAT1 & EMAT2.
C-------------------------------------
C
      IF (.NOT. NEWGAM) CALL DZERO(WORK(KGAMMA),NGAMMA(ISYMOP))
      CALL DZERO(WORK(KEMAT1),NEMAT1(ISYMOP))
      CALL DZERO(WORK(KEMAT2),NMATIJ(ISYMOP))
      IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
        CALL DZERO(WORK(KE1PIM),NADP(1))
      END IF
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
C------------------------------------------------
C     Read one-electron integrals in Fock-matrix.
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
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C SLV98,OC
C-------------------------------------
C
      IF (CCSLV .AND. (.NOT. CCMM )) THEN
         CALL CCSL_RHSTG(WORK(KFOCK),WORK(KEND1),LWRK1)
      ENDIF
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C CCMM02,JA+AO
C-------------------------------------
C
      IF (CCMM) THEN
         IF (.NOT. NYQMMM) THEN
            CALL CCMM_RHSTG(WORK(KFOCK),WORK(KEND1),LWRK1)
         ELSE IF (NYQMMM) THEN
            IF (HFFLD ) THEN
              CALL CCMM_ADDGHF(WORK(KFOCK),WORK(KEND1),LWRK1)
            ELSE
              CALL CCMM_ADDG(WORK(KFOCK),WORK(KEND1),LWRK1)
            END IF
         END IF
      ENDIF
C
      IF (IPRINT .GT.15) THEN
         CALL AROUND( 'Fock AO matrix after ff/slv/mm contribution' )
         CALL CC_PRFCKAO(WORK(KFOCK),1)
      ENDIF
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
               IF (T2TCOR) THEN
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
               END IF
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
                  CALL QUIT('Insufficient space in CCRHSN')
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
C-----------------------------------------------------------
C              Calculate transformed integrals used in t3am.
C-----------------------------------------------------------
C
               IF (CCSDT .AND. ((.NOT. CC1B) .OR. (.NOT. CC1A))) THEN
C
                  CALL CC3_T3INT(WORK(KXINT),WORK(KLAMDP),WORK(KLAMDH),
     *                           T1AM,1,WORK(KEND2),LWRK2,IDEL,ISYMD,1,
     *                           LU3SRT,FN3SRT,LUCKJD,FNCKJD)
C
               ENDIF
C
C-------------------------------------------------------------------
C              Calculate additional integrals needed for CCSD(R12)/2
C-------------------------------------------------------------------
C
               IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
                  LENIAJ = NT2BCD(ISYDIS)

                  KXIADJ = KEND2
                  KXIJDA = KXIADJ + LENIAJ
                  KEND3  = KXIJDA + LENIAJ
                  LWRK3 = LWORK - KEND3
                  IF (LWRK3 .LT. 0) THEN
                    WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                    CALL QUIT('Insufficient space in CCRHSN')
                  ENDIF

                  CALL DZERO(WORK(KXIADJ),LENIAJ)
                  CALL DZERO(WORK(KXIJDA),LENIAJ)

                  DO ISYGAM = 1, NSYM
                     ISYALBE = MULD2H(ISYDIS,ISYGAM)
                  DO G = 1, NBAS(ISYGAM)
                     IGAM = G + IBAS(ISYGAM)  
                 
                     KOFFG = KXINT + IDSAOG(ISYGAM,ISYDIS)
     &                       + NNBST(ISYALBE)*(G-1)

                     CALL CC_IAJB( WORK(KOFFG), ISYALBE, DUMMY, ISYM0,
     &                             IDEL, IGAM, .FALSE., IDUMMY,
     &                             DUMMY, WORK(KXIADJ), WORK(KXIJDA),
     &                             DUMMY, DUMMY, DUMMY, 
     &                             WORK(KLAMDP), WORK(KLAMDH), ISYM0,
     &                             DUMMY, DUMMY, ISYM0,
     &                             WORK(KLAMDP), WORK(KLAMDH), ISYM0,
     &                             DUMMY, DUMMY, ISYM0,
     &                             WORK(KEND3), LWRK3,   3,    
     &                             .FALSE., .FALSE.,  .TRUE.,   
     &                             .FALSE., .FALSE.,  0      )
                  END DO
                  END DO

c                 ------------------------------------
c                 update Fhat_{del a}:
c                 ------------------------------------
                  D = IDEL - IBAS(ISYMD)
                  CALL CC_FCKDELA(D,ISYMD,WORK(KFCKVAO),ISYM0,
     &                            WORK(KXIJDA),WORK(KXIADJ),IEMAT1)

C                 ------------------------------------
C                 transform (ia|del j) to L(ia|del j):
C                 ------------------------------------
                  CALL DSCAL(LENIAJ, TWO,WORK(KXIADJ),1)
                  CALL DAXPY(LENIAJ,-ONE,WORK(KXIJDA),1,
     *                                           WORK(KXIADJ),1)
                  
C                 --------------------------------------------
C                 write 3-index transformed integrals to disk:
C                 --------------------------------------------
                  IADR = IT2DEL(IDEL) + 1
                  CALL PUTWA2(LUIADJ,FNIADJ,WORK(KXIADJ),IADR,LENIAJ)
                  CALL PUTWA2(LUIJDA,FNIJDA,WORK(KXIJDA),IADR,LENIAJ)

               END IF
C
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
                  CALL QUIT('Insufficient space in CCRHSN')
               ENDIF
C
C----------------------------------------------------------------
C              Construct the partially transformed T2-amplitudes.
C----------------------------------------------------------------
C
               DTIME   = SECOND()
               ICON = 1
               ISYMLH = 1
               CALL CC_T2AO(T2AM,WORK(KLAMDH),ISYMLH,WORK(KSCRM),
     *                         WORK(KEND3),LWRK3,IDEL,ISYMD,
     *                         ISYMTR,ICON)
               DTIME   = SECOND() - DTIME
               TIMT2AO = TIMT2AO + DTIME
C
C-----------------------------------
C              Calculate the F-term.
C-----------------------------------
C
               DTIME   = SECOND()
               IF (.NOT. OMEGOR) THEN
                  CALL CCRHS_F(WORK(KXINT),OMEGA2,WORK(KLAMDH),
     *                         WORK(KEND3),LWRK3,IDEL,ISYMD)
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMF    = TIMF     + DTIME
C
C-------------------------------------------------------
C              Calculate the F-term in MO basis for CC2.
C-------------------------------------------------------
C
               IF ( CC2 ) THEN
                  DTIME = SECOND() - TIMFP
                  IOPT = 1
                  CALL GETTIM(T0,W0)
                  LVIJKL = .NOT.USEVABKL .AND. CC2R12
                  LVAJKL = LVIJKL .AND. RSPIM
                  CALL CC_MOFCON(WORK(KXINT),OMEGA2,
     *                           WORK(KLAMDP),WORK(KLAMDH),
     *                           WORK(KLAMDP),WORK(KLAMDH),
     *                           WORK(KEND3),LWRK3,IDEL,
     *                           ISYMD,ISYMTR,IOPT,
     *                           WORK(KVIJKL),LVIJKL,IANR12,
     *                           WORK(KVAJKL),LVAJKL,TIMFP)
                  CALL GETTIM(T1,W1)
                  TIMMOFCPU = T1-T0
                  TIMMOFWAL = W1-W0
                  DTIME   = (SECOND() - TIMFP) - DTIME
                  TIMF    = TIMF     + DTIME
               ENDIF
C
C-----------------------------------
C              Calculate the B-term.
C-----------------------------------
C
               DTIME   = SECOND()
               IF ((.NOT. OMEGOR) .AND. (.NOT. CC2)) THEN
                  CALL CCRHS_B(WORK(KXINT),OMEGA2,WORK(KLAMDP),
     *                         WORK(KLAMDH),WORK(KSCRM),WORK(KEND3),
     *                         LWRK3,IDEL,ISYMD)
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMB    = TIMB     + DTIME
C
C------------------------------------------
C              Calculate the B and F terms.
C------------------------------------------
C
               DTIME   = SECOND()
               IF (OMEGOR .AND. ( .NOT. CC2) ) THEN
                  IOPT = 1
                  CALL CC_BF(WORK(KXINT),OMEGA2,WORK(KLAMDH),1,
     *                       WORK(KLAMDH),1,WORK(KLAMDH),1,
     *                       WORK(KSCRM),ISYMD,DUMMY,ISYMD,
     *                       WORK(KEND3),LWRK3,IDEL,ISYMD,IOPT)
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMBF   = TIMBF    + DTIME
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
                  CALL QUIT('Insufficient space in CCRHSN')
               ENDIF
C
C--------------------------------------------------------
C              Transform one index in the integral batch.
C--------------------------------------------------------
C
               DTIME   = SECOND()
               ISYMLP  = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                     ISYMLP,WORK(KEND4),LWRK4,ISYDIS)
               DTIME   = SECOND() - DTIME
               TIMTRBT = TIMTRBT + DTIME
C
C-------------------------------------------------------------
C              Calculate the gamma matrix entering the A-term.
C-------------------------------------------------------------
C
               DTIME   = SECOND()
               IF ((.NOT. CC2) .AND. (.NOT. NEWGAM)) THEN
                 CALL CCRHS_GAM(WORK(KDSRHF),WORK(KGAMMA),WORK(KLAMDP),
     *                           WORK(KLAMDH),WORK(KSCRM),WORK(KEND4),
     *                           LWRK4,IDEL,ISYMD)
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMGAM  = TIMGAM   + DTIME
C
C-----------------------------------
C              Calculate the C-term.
C-----------------------------------
C
               DTIME   = SECOND()
C
               IF ( RSPIM ) THEN
                  FACTC = XMONE
               ELSE
                  FACTC = XMHALF
               ENDIF
C
               ICON = 2
               IV = 1
C
               IF (CCSDR12 .AND. (IANR12.EQ.2 .OR.IANR12.EQ.3)) THEN
                 IOPTR12 = 1
                 IOPTE = 1
               ELSE
                 IOPTR12 = 0
                 IOPTE = 0
               END IF
C
               IF (.NOT. T2TCOR) THEN
                  CALL CCRHS_C(WORK(KXINT),WORK(KDSRHF),OMEGA2,
     *                         T2AM,ISYMOP,WORK(KLAMDP),WORK(KLAMIP),
     *                         WORK(KLAMDH),WORK(KLAMDP),ISYMTR,
     *                         WORK(KLAMDP),ISYMTR,
     *                         WORK(KSCRM),WORK(KE1PIM),WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,FACTC,ICON,IOPTR12,
     *                         IOPTE,LUC,CFIL,LUCP,CPFIL,IV)
               ELSE
                  CALL CCRHS_C(WORK(KXINT),WORK(KDSRHF),OMEGA2,
     *                         WORK(KT2AMT),ISYMOP,
     *                         WORK(KLAMDP),WORK(KLAMIP),
     *                         WORK(KLAMDH),WORK(KLAMDP),ISYMTR,
     *                         WORK(KLAMDP),ISYMTR,
     *                         WORK(KSCRM),WORK(KE1PIM),WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,FACTC,ICON,IOPTR12,
     *                         IOPTE,LUC,CFIL,LUCP,CPFIL,IV)
               END IF
CTesT
C              WRITE(LUPRI,*) 'E1PIM after CCRHS_C:'
C              WRITE(LUPRI,*) 'Norm^2: ',
C    &                       DDOT(NADP(1),WORK(KE1PIM),1,WORK(KE1PIM),1)
C              DO ISYM = 1,NSYM
C                CALL OUTPUT(WORK(KE1PIM+IADP(ISYM,ISYM)),
C    &                       1,NVIR(ISYM),1,MBAS1(ISYM)+MBAS2(ISYM),
C    &                       NVIR(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
C    &                       1, LUPRI)
C              END DO
C              CALL FLSHFO(LUPRI)
CTesT
C
               DTIME   = SECOND() - DTIME
               TIMC    = TIMC     + DTIME
C
C---------------------------------------
C              Transform T2 to 2T2 - T2.
C---------------------------------------
C
               DTIME   = SECOND()
               IF (T2TCOR) THEN
                  CALL DSCAL(NT2SQ(1),TWO,T2AM,1)
                  CALL DAXPY(NT2SQ(1),-ONE,WORK(KT2AMT),1,T2AM,1)
               ELSE
                  ISYM = 1
                  CALL CCRHS_T2TR(T2AM,WORK(KEND4),LWRK4,ISYM)
               END IF
               DTIME   = SECOND() - DTIME
               TIMT2TR = TIMT2TR  + DTIME
C
C-----------------------------------------------
C              Transform the cluster amplitudes.
C-----------------------------------------------
C
               CALL CC_MTCME(WORK(KSCRM),WORK(KEND4),LWRK4,
     *                       ISYMD,ISYMTR)
C
C-----------------------------------
C              Calculate the D-term.
C-----------------------------------
C
               DTIME   = SECOND()
C
               IF ( RSPIM ) THEN
                  FACTD = ONE
               ELSE
                  FACTD = HALF
               ENDIF
C
               ICON   = 2
               IV = 1
C
               IF (CCSDR12 .AND. (IANR12.EQ.2 .OR. IANR12.EQ.3)) THEN
                 IOPTR12 = 1
                 IOPTE = 1
               ELSE
                 IOPTR12 = 0
                 IOPTE = 0
               END IF
C
               IF ( .NOT. CC2) THEN
                  CALL CCRHS_D(WORK(KXINT),WORK(KDSRHF),OMEGA2,T2AM,
     *                         ISYMTR,WORK(KLAMDP),WORK(KLAMIP),
     *                         WORK(KLAMDH),WORK(KLAMDP),ISYMTR,
     *                         WORK(KLAMDH),ISYMTR,
     *                         WORK(KSCRM),WORK(KE1PIM),WORK(KEND4),
     *                         LWRK4,IDEL,ISYMD,FACTD,ICON,IOPTR12,
     *                         IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
               ENDIF
CTesT
C              WRITE(LUPRI,*) 'E1PIM after CCRHS_D:'
C                            WRITE(LUPRI,*) 'Norm^2: ',
C    &                       DDOT(NADP(1),WORK(KE1PIM),1,WORK(KE1PIM),1)
C              DO ISYM = 1,NSYM
C                CALL OUTPUT(WORK(KE1PIM+IADP(ISYM,ISYM)),
C    &                       1,NVIR(ISYM),1,MBAS1(ISYM)+MBAS2(ISYM),
C    &                       NVIR(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
C    &                       1, LUPRI)
C              END DO
C              CALL FLSHFO(LUPRI)
CTesT
C
               DTIME   = SECOND() - DTIME
               TIMD    = TIMD     + DTIME
C
C----------------------------------------
C              Calculate E-intermediates.
C----------------------------------------
C
               DTIME   = SECOND()
               IF ((.NOT. CC2) .OR. RSPIM) THEN
                  CALL CCRHS_EI(WORK(KDSRHF),WORK(KEMAT1),WORK(KEMAT2),
     *                          T2AM,WORK(KSCRM),WORK(KLAMDP),
     *                          WORK(KLAMDH),WORK(KEND4),LWRK4,
     *                          IDEL,ISYMD,ISYDIS,ISYMTR)
               ENDIF
               DTIME   = SECOND() - DTIME
               TIMEI   = TIMEI    + DTIME
C
C-----------------------------------
C              Calculate the G-term.
C-----------------------------------
C
               DTIME   = SECOND()
               ISYMP1 = 1
               ISYMH1 = 1
               CALL CCRHS_G(WORK(KDSRHF),OMEGA1,WORK(KLAMDP),ISYMP1,
     *                      WORK(KLAMDH),ISYMH1,WORK(KSCRM),WORK(KEND4),
     *                      LWRK4,ISYDIS,ISYMD,ISYMTR)
               DTIME   = SECOND() - DTIME
               TIMG    = TIMG     + DTIME
C
C-----------------------------------
C              Calculate the H-term.
C-----------------------------------
C
               DTIME   = SECOND()
               CALL CCRHS_H(WORK(KDSRHF),OMEGA1,WORK(KLAMDP),
     *                      WORK(KLAMDH),WORK(KSCRM),WORK(KEND4),
     *                      LWRK4,ISYDIS,ISYMD,ISYMTR)
               DTIME   = SECOND() - DTIME
               TIMH    = TIMH     + DTIME
C
C---------------------------------------------
C              BackTransform T2 from 2T2 - T2.
C---------------------------------------------
C
               DTIME   = SECOND()
               IF (T2TCOR) THEN
                  CALL DAXPY(NT2SQ(1),ONE,WORK(KT2AMT),1,T2AM,1)
                  CALL DSCAL(NT2SQ(1),HALF,T2AM,1)
               ELSE
                  ISYM = 1
                  CALL CCRHS_T2BT(T2AM,WORK(KEND4),LWRK4,ISYM)
               END IF
               DTIME   = SECOND() - DTIME
               TIMT2BT = TIMT2BT  + DTIME
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C     ------------------------------------------------------------------
C     save the special fock matrix computed for CCSD(R12) ansaetze 2/3
C     ------------------------------------------------------------------
      IF (CCSDR12 .AND. (IANR12.EQ.2 .OR. IANR12.EQ.3)) THEN
        LUFHATAD = -1
        CALL GPOPEN(LUFHATAD,'CCFHATADEL','UNKNOWN',' ','UNFORMATTED',
     &              IDUMMY,.FALSE.)
        REWIND(LUFHATAD)
        WRITE(LUFHATAD) (WORK(KFCKVAO-1+I),I=1,NEMAT1(1))
        CALL GPCLOSE(LUFHATAD,'KEEP')
      END IF
C
C     ------------------------------------------------------------------
C     for CCSD(R12) ansaetze 2/3 do here the C, D, and E terms requiring 
C     the calculation of integrals with delta from the  auxiliary basis
C     ------------------------------------------------------------------
C
      IF (CCSDR12 .AND. (IANR12.EQ.2 .OR. IANR12.EQ.3)) THEN
         CALL CCSDR12AO(CCSDR12,
     &                  T2AM,WORK(KLAMDP),WORK(KLAMDH),
     &                  FNIADJ,LUIADJ,FNIJDA,LUIJDA,
     &                  CPFIL,LUCP,DPFIL,LUDP,WORK(KE1PIM),
     &                  TIMINTR12,TIMRDAOR12,TIMTRBT,
     &                  TIMC,TIMD,TIMT2TR,TIMT2BT,
     &                  WORK(KEND1T),LWRK1T)
      END IF
CTesT
C        WRITE(LUPRI,*) 'E1PIM after CCSDR12AO:'
C        WRITE(LUPRI,*) 'Norm^2: ',
C    &                 DDOT(NADP(1),WORK(KE1PIM),1,WORK(KE1PIM),1)
C        DO ISYM = 1,NSYM
C          CALL OUTPUT(WORK(KE1PIM+IADP(ISYM,ISYM)),
C    &                 1,NVIR(ISYM),1,MBAS1(ISYM)+MBAS2(ISYM),
C    &                 NVIR(ISYM),MBAS1(ISYM)+MBAS2(ISYM),
C    &                 1, LUPRI)
C        END DO
C        CALL FLSHFO(LUPRI)
CTesT
C
C------------------------
C     Recover work space.
C------------------------
C
      IF (CCSDT) THEN
         KEND1 = KEND1T
         LWRK1 = LWRK1T
      ELSE
         KEND1 = KENDS2
         LWRK1 = LWRKS2
      ENDIF
C
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('After  Delta Loop: Omega1')
         CALL CC_PRP(OMEGA1,OMEGA2,1,1,0)
      ENDIF
C     
C     ----------------------------------------------------------------
C     for CC2-R12 ansatz 3 add (ai|bj)-hat x (ka|r12|lb) to V intermediate
C     (Note: this requires that omega2 contains the integral 
C            (ai|bj)-hat in packed triangular storage)
C     ----------------------------------------------------------------
C
      IF (CC2 .AND. CCR12 .AND. IANR12.EQ.3) THEN
         ! get R12 integrals
         lunit = -1
         call gpopen(lunit,fr12r12,'unknown',' ','unformatted',
     &              idum,.false.)
         read(lunit)(t2am(i),i=1,nt2r12(1))  
         call gpclose(lunit,'KEEP')
         
         CALL CC_R12MI2(WORK(KVIJKL),T2AM,OMEGA2,1,1,-1.0d0,
     &                  WORK(KEND1),LWRK1)

         ! restore amplitudes stored as full square matrix
         IF (LWRK1.LT.NT2AMX) CALL QUIT('Out of memory in CCRHSN')
         REWIND (LURHS1)
         READ (LURHS1)
         READ (LURHS1) (WORK(KEND1+I-1), I = 1,NT2AMX)
         CALL CC_T2SQ(WORK(KEND1),T2AM,1)
      END IF
C
C---------------------------------------------------------------------
C     for CC2 and NONHF=.true. calculate Fock matrix entering E-terms:
C       the SCF Fock matrix is in principle given by the SCF orbital
C       energies, but in recomputing it here from the SCF AO-Fock
C       matrix computed in CCSD_IAJB allows to do finite difference
C       on the vector function with respect to the CMO vector
C       (see CC_FDXI & CC_FDETA). Note the SCF AO-Fock matrix read
C       from file includes the `relaxed' external fields, so we
C       only have to add the unrelaxed fields.
C---------------------------------------------------------------------
C
      DTIME = SECOND()
      IF ((CC2 .OR. CCR12) .AND. NONHF) THEN
        KFIELD = KEND1
        KEND2  = KFIELD + N2BAST
        IF (CC2) THEN
          KCMO   = KEND2
          KEND2  = KCMO + MAX(NLAMDT,NLAMDS)
        END IF 
        IF (CCR12) THEN
          if (isymop.ne.1) call quit('Symmetry problem in CCSD_RHS')
          kvxintsq   = kend2
          kxint    = kvxintsq + nr12r12sq(isymop)
          kxintsq  = kxint + nr12r12p(1)
          ktr12    = kxintsq + nr12r12sq(1)
          ktr12sq  = ktr12 + ntr12am(1)
          kxir12   = ktr12sq + ntr12sq(1)
          kend2    = kxir12 + ntr12sq(1)
        END IF
        LWRK2  = LWORK  - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insufficient memory in CCRHSN.')
        END IF

        CALL DZERO(WORK(KFIELD),N2BAST)
        IF (CCR12) THEN
          CALL DZERO(WORK(KVXINTSQ),NR12R12SQ(1))
        END IF
        DO  IF = 1, NFIELD
          IF ( NHFFIELD(IF) ) THEN
            DTIME = SECOND()
            CALL CC_ONEP(WORK(KFIELD),WORK(KEND2),LWRK2,EFIELD(IF),1,
     *                   LFIELD(IF))
            TIMFCKMO = TIMFCKMO + SECOND() - DTIME
            IF (CCR12) THEN
              DTIME = SECOND()
              CALL CC_R12RDVXINT(WORK(KVXINTSQ),WORK(KEND2),LWRK2,
     &                         EFIELD(IF),1, LFIELD(IF))
              TIMEP = TIMEP + SECOND() - DTIME
            END IF
          ELSE IF (.NOT. NHFFIELD(IF) .AND. CCR12) THEN
            CALL QUIT('CCR12 response can only handle unrelaxed '//
     &                'orbitals (w.r.t. the perturbation)')
          END IF
        END DO

        IF (CCR12) THEN
            DTIME = SECOND()
            ! read R12 amplitudes and reorder to full square
            iopt=32
            call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(ktr12))
            iopt = 1
            call ccr12unpck2(work(ktr12),1,work(ktr12sq),'N',iopt)
 
            ! read R12 overlap matrix and reorder to full square
            lunit = -1
            call gpopen(lunit,fccr12x,'old',' ','unformatted',idummy,
     &                  .false.)
            rewind(lunit)
 8888       read(lunit) ian
            read(lunit) (work(kxint-1+i), i=1, nr12r12p(1))
            if (ian.ne.ianr12) goto 8888 
            call gpclose(lunit,'KEEP')
            iopt = 2
            call ccr12unpck2(work(kxint),1,work(kxintsq),'N',iopt)

            ! calculate R12 response contribution to Omega_{kilj}: 
            CALL CC_R12XI(work(kxir12),1,'T',work(ktr12sq),1,
     &                    work(kxintsq),work(kvxintsq),1,work(kfield),
     &                    work(klamdp),work(klamdh),'N',work(kend2),
     &                    lwrk2)

C           ! transpose Xi: in Xi the r12-pair index (kl) is leading,
C           ! in Vijkl the occ. index pair (ij) is leading!!! 
C           call cclr_trsqr12(work(kxir12),1)
            
            ! add it to Omega_FP term = VIJKL 
            call daxpy(ntr12sq(1),one,work(kxir12),1,work(kvijkl),1)

            TIMEP = TIMEP + SECOND() - DTIME
        END IF

        IF (CC2) THEN
          DTIME = SECOND()
          LUSIFC = -1
          CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *                IDUMMY,.FALSE.)
          REWIND(LUSIFC)
          CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
          READ(LUSIFC)
          READ(LUSIFC)
          READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
          CALL GPCLOSE(LUSIFC,'KEEP')
          
          CALL CMO_REORDER(WORK(KCMO),WORK(KEND2),LWRK2)   

          LUFCK = -1
          CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *                IDUMMY,.FALSE.)
          REWIND(LUFCK)
          READ(LUFCK)(WORK(KFCKHF + I-1),I = 1,N2BST(ISYMOP))
          CALL GPCLOSE(LUFCK,'KEEP' )
          
          ! SCF Fock matrix in transformed using CMO vector
          CALL CC_FCKMO(WORK(KFCKHF),WORK(KCMO),WORK(KCMO),
     *                  WORK(KEND2),LWRK2,1,1,1)
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C SLV98,OC
C-------------------------------------
C
        IF (CCSLV .AND. (.NOT. CCMM )) THEN
           CALL CCSL_RHSTG(WORK(KFIELD),WORK(KEND2),LWRK2)
        ENDIF
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C CCMM02,JA+AO
C-------------------------------------
C
        IF (CCMM) THEN
           IF (.NOT. NYQMMM) THEN
              CALL CCMM_RHSTG(WORK(KFIELD),WORK(KEND2),LWRK2)
           ELSE IF (NYQMMM) THEN
              IF ( HFFLD ) THEN
                CALL CCMM_ADDGHF(WORK(KFIELD),WORK(KEND2),LWRK2)
              ELSE
                CALL CCMM_ADDG(WORK(KFIELD),WORK(KEND2),LWRK2)
              END IF
           END IF
        ENDIF
C
C----------------------------------
C
          ! unrelaxed fields are transformed using the Lambda matrices
          CALL CC_FCKMO(WORK(KFIELD),WORK(KLAMDP),WORK(KLAMDH),
     *                  WORK(KEND2),LWRK2,1,1,1)

          CALL DAXPY(N2BAST,ONE,WORK(KFIELD),1,WORK(KFCKHF),1)
          TIMFCKMO = TIMFCKMO + SECOND() - DTIME
        END IF

      END IF

C
C------------------------------------------------------------------
C     for CCSD(R12) add the R12 contribution to the BF intermediate
C     which at this place is (still) stored in OMEGA2:
C------------------------------------------------------------------
      IF (CCSDR12) THEN
        TIMR12CPU = 0.0d0
        TIMR12WAL = 0.0d0
        CALL GETTIM(T0,W0)
c
        IOPT = 0
        IAMP = 0
        CALL CCRHS_BP(OMEGA2,1,IOPT,IAMP,DUMMY,IDUMMY,IDUMMY,DUMMY,
     &                IDUMMY,DUMMY,WORK(KEND1),LWRK1)
c
        CALL GETTIM(T1,W1)
        IF (IPRINT .GT. 9) THEN
          WRITE(LUPRI,*)'Time used for CCRHS_BP cpu:', T1-T0 
          WRITE(LUPRI,*)'Time used for CCRHS_BP wall:', W1-W0 
        END IF
        TIMR12CPU = TIMR12CPU + (T1-T0)
        TIMR12WAL = TIMR12WAL + (W1-W0)
      END IF
C
C-------------------------------------------------
C     for CC-R12:
C-------------------------------------------------
C
      IF (CCR12) THEN
        TIMR12CPU = 0.0d0
        TIMR12WAL = 0.0d0
        CALL GETTIM(T0,W0)
        IF (.NOT.USEVABKL) THEN
          LVIJKL = .TRUE.
          LVAJKL = RSPIM 
          LVABKL = .FALSE.
          IOPTBAS = 1
          IF (R12CBS .AND. (IANR12.NE.1)) IOPTBAS = 2
          FACTERM23 = TWO
          CALL CC_MOFCONR12(WORK(KLAMDH),1,WORK(KLAMDHS),
     &                      WORK(KLAMDPS),WORK(KLAMDHS),ISYMTR,
     &                      WORK(KVIJKL),FACTERM23,WORK(KVAJKL),IDUMMY,
     &                      LVIJKL,LVAJKL,LVABKL,IOPTBAS,
     &                      TIMRDAOR12,TIMFP,TIMINTR12,
     &                      IGLMRHS,NGLMDS,IMAIJM,NMAIJM,
     &                      IMAKLM,NMAKLM,WORK(KEND1),LWRK1)
C
C         write V(alpha jtilde,kl) to disk
C
          IF (RSPIM) THEN
            IF (IANR12.EQ.2.OR.IANR12.EQ.3) THEN
C             calculate contributions for ansatz 2
              ISYMH = ISYMTR
              ISYMV = 1
              CALL CC_R12MKVAJ2(WORK(KVAJKL),ISYMV,WORK(KLAMDH),ISYMH,
     &                          WORK(KLAMDHS),ISYMH,WORK(KEND1),LWRK1)
            END IF
C           WRITE(LUPRI,*)'write Vajtkl on disk'
            LUVAJTKL = -1
            CALL GPOPEN(LUVAJTKl,FVAJTKL,'UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND(LUVAJTKL)
            WRITE(LUVAJTKL) (WORK(KVAJKL+I-1), I = 1,NVAJKL(1))
            CALL GPCLOSE(LUVAJTKL,'KEEP')
          END IF
        ELSE
          KVABKL = KEND1
          KVAJKL = KVABKL + NVABKL(1)
          KEND2  = KVAJKL + NVAJKL(1)
          LWRK2  = LWORK - KEND2
          IF (LWRK2.LT.0) THEN
            CALL QUIT('Insufficient work space in ccrhsn')
          END IF
          ISYMC = 1
          LV = .TRUE.
          LVIJKL = .TRUE.
          LVAJKL = RSPIM
c
          CALL CC_R12MKVTF(WORK(KVABKL),WORK(KVAJKL),WORK(KVIJKL),
     &                     WORK(KLAMDH),ISYMC,
     &                     LV,LVIJKL,LVAJKL,FVAJTKL,WORK(KEND2),LWRK2)
c
        END IF
        CALL GETTIM(T1,W1)
        IF (IPRINT .GT. 9) THEN
          WRITE(LUPRI,*)'Time used for F''-term cpu:', T1-T0
          WRITE(LUPRI,*)'Time used for F''-term wall:',W1-W0
        END IF
        TIMR12CPU = TIMR12CPU + (T1-T0)
        TIMR12WAL = TIMR12WAL + (W1-W0)
        CALL GETTIM(T1,W1)
        TIMMOFR12CPU = T1-T0
        TIMMOFR12WAL = W1-W0
C
C-------------------------------------------------
C     for CC2-R12:
C-------------------------------------------------
C
        IF (CC2) THEN
          IF (IANR12.EQ.2 .OR. IANR12.EQ.3) THEN
            ISYMV = 1 
            ISYMH = 1
            CALL GETTIM(T0,W0)
            CALL CC_R12INTF2(WORK(KVIJKL),WORK(KLAMDH),ISYMH,
     &                       WORK(KLAMDHS),ISYMV,WORK(KLAMDHS),ISYMH,
     &                       WORK(KEND1),LWRK1)
            CALL GETTIM(T1,W1)
            TIMINTF2CPU = T1-T0
            TIMINTF2WAL = W1-W0
            TIMR12CPU = TIMR12CPU + (T1-T0)
            TIMR12WAL = TIMR12WAL + (W1-W0)
          END IF

          IF ((IANR12.EQ.2.OR.IANR12.EQ.3).AND.DEBUGV) THEN
c           symmetrize Vijkl
            ISYMV = 1
            KVSYM = KEND1
            KEND1 = KVSYM + NTR12SQ(1)
            CALL SYMV(WORK(KVIJKL),ISYMV,WORK(KVSYM),
     &                NRHF,IMATIJ,ITR12SQT,NMATIJ,WORK(KEND1),LWRK1)

c           write V^ij_kl on file to calculate later numerically V bar
            LUVIJKL = -1
            CALL GPOPEN(LUVIJKL,FVIJKL,'UNKNOWN',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            WRITE(LUVIJKL)(WORK(KVSYM-1+I),I=1,NTR12SQ(1))
            CALL GPCLOSE(LUVIJKL,'KEEP')
            WRITE(LUPRI,*)'VIJKL WRITTEN ON FILE'

            DO ISYMIJ = 1, NSYM
               ISYMKL = MULD2H(ISYMIJ,ISYMTR)
               WRITE(LUPRI,*) 'ISYMIJ,ISYMKL:',ISYMIJ,ISYMKL
               CALL OUTPUT(WORK(KVSYM+ITR12SQT(ISYMIJ,ISYMKL)),1,
     &              NMATIJ(ISYMIJ),1,NMATKL(ISYMKL),NMATIJ(ISYMIJ),
     &              NMATKL(ISYMKL),1,LUPRI)
            END DO
          END IF
        END IF
C
        IF (IANR12.EQ.2 .OR. IANR12.EQ.3) THEN
          LRES = .FALSE.
          CALL GETTIM(T0,W0)
          CALL CCRHS_EPP(WORK(KVIJKL),T2AM,1,WORK(KEND1),LWRK1,
     &                   APROXR12,LRES,IDUMMY,CDUMMY,IDUMMY)
          CALL GETTIM(T1,W1)
          TIMEPPCPU = T1-T0
          TIMEPPWAL = W1-W0
          TIMR12CPU = TIMR12CPU + (T1-T0)
          TIMR12WAL = TIMR12WAL + (W1-W0)
c
          CALL GETTIM(T0,W0)
          IOPTE = 0
          CALL CCRHS_HP(OMEGA1,WORK(KLAMDH),ISYMH,WORK(KLAMDH),ISYMH,
     &                  WORK(KEND1),LWRK1,0,1,CDUMMY,IDUMMY,IDUMMY,
     &                  IOPTE)  
          CALL GETTIM(T1,W1)
          TIMHPCPU = T1-T0
          TIMHPWAL = W1-W0
          TIMR12CPU = TIMR12CPU + (T1-T0)
          TIMR12WAL = TIMR12WAL + (W1-W0)
c
          CALL GETTIM(T0,W0)
          CALL CCRHS_IP(OMEGA1,T1AM,1,WORK(KLAMDH),ISYMH,0,1,
     &                  CDUMMY,IDUMMY,IDUMMY,WORK(KEND1),LWRK1)
c
          CALL GETTIM(T1,W1)
          TIMIPCPU = T1-T0
          TIMIPWAL = W1-W0
          TIMR12CPU = TIMR12CPU + (T1-T0)
          TIMR12WAL = TIMR12WAL + (W1-W0)
        END IF
C
        IF (CCSDR12) THEN
          CALL GETTIM(T0,W0)
c
          CALL CCRHS_BPP(WORK(KVIJKL),T2AM,1,.FALSE.,
     &                   FVCDKL,1,WORK(KEND1),LWRK1)
c
          CALL GETTIM(T1,W1)
          IF (IPRINT .GT. 9) THEN
            WRITE(LUPRI,*)'Time used for CCRHS_BPP cpu:', T1-T0 
            WRITE(LUPRI,*)'Time used for CCRHS_BPP wall:',W1-W0 
          END IF
          TIMR12CPU = TIMR12CPU + (T1-T0)
          TIMR12WAL = TIMR12WAL + (W1-W0)
        END IF
c
        ISYMV = 1 
        CALL GETTIM(T0,W0)
c
        CALL CCRHS_EP(WORK(KVIJKL),ISYMV,.FALSE.,DUMMY,
     &                WORK(KEND1),LWRK1,0,
     &                CDUMMY,IDUMMY,CDUMMY,IDUMMY,IDUMMY,APROXR12,
     &                BRASCL,KETSCL)
        CALL GETTIM(T1,W1)
        TIMEPCPU = T1-T0
        TIMEPWAL = W1-W0
        TIMR12CPU = TIMR12CPU + (T1-T0)
        TIMR12WAL = TIMR12WAL + (W1-W0)
c
        IF ((IANR12.EQ.2.OR.IANR12.EQ.3).AND.DEBUGV) THEN
c         write V^ij_kl on file to calculate later numerically RHOR12 
          LUVIJKL = -1
          CALL GPOPEN(LUVIJKL,FVIJKL,'UNKNOWN',' ','UNFORMATTED',
     &                IDUMMY,.FALSE.)
          WRITE(LUVIJKL)(WORK(KVIJKL-1+I),I=1,NTR12SQ(1))
          CALL GPCLOSE(LUVIJKL,'KEEP')
          WRITE(LUPRI,*)'VIJKL WRITTEN ON FILE'
        END IF
c
        CALL GETTIM(T0,W0)
        CALL CCRHS_GP(OMEGA1,WORK(KLAMDP),
     &       WORK(KEND1),LWRK1,0,1,CDUMMY,IDUMMY,IDUMMY)
        CALL GETTIM(T1,W1)
        TIMGPCPU = T1-T0
        TIMGPWAL = W1-W0
        TIMR12CPU = TIMR12CPU + (T1-T0)
        TIMR12WAL = TIMR12WAL + (W1-W0)
c       TIMGP = TIMGP + ( SECOND() - DTIME ) 

      END IF !CCR12
C
C-------------------------------------------------
C     Transform the Omega2 vector to the MO basis.
C-------------------------------------------------
C
      IF (NT2AM(ISYMOP) .GT. 2*NT2AMX) THEN
         WRITE(LUPRI,*)
     &        'Length of T2AM is smaller than OMEGA2 in MO basis'
         CALL QUIT('Insufficient space in CC_T2MO')
      ENDIF
C
      IF ( .NOT. CC2 ) THEN
C
C---------------------------------------
C        Save the CC amplitudes on disk.
C---------------------------------------
C
         WRITE (LURHS1) (T2AM(I), I = 1,NT2AM(ISYMOP))
C
C----------------------------------------------------------------------
C        Write Omega2 vector to disk if needed in response calculation.
C----------------------------------------------------------------------
C
         IF ( RSPIM ) THEN
C
            LUBF = -1
            CALL GPOPEN(LUBF,'CC_BFIM','UNKNOWN',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
            REWIND(LUBF)
            WRITE(LUBF) (OMEGA2(I),I = 1,2*NT2ORT(1))
            CALL GPCLOSE(LUBF,'KEEP')
C
         ENDIF
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
            CALL QUIT('Insufficient space in CC_T2MO')
         ENDIF
C
         TIMOME2 = SECOND()
         ISYMBF = ISYMOP
         ICON = 1
         
         CALL CC_T2MO(FAKE,PHONEY,ISYMOP,OMEGA2,T2AM,WORK(KGAMMA),
     *                WORK(KLAMDP),WORK(KLAMDP),ISYMTR,
     *                WORK(KEND1),LWRK1,ISYMBF,ICON)
         CALL DCOPY(NT2AM(ISYMTR),T2AM,1,OMEGA2,1)
         TIMOME2 = SECOND() - TIMOME2
C
         IF (IPRINT .GT. 51) THEN
            RHO1N = DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
            RHO2N = DDOT(NT2AMX,OMEGA2,1,OMEGA2,1)
            WRITE(LUPRI,*) 'Norm of OMEGA1 -after CC_T2MO: ',RHO1N
            WRITE(LUPRI,*) 'Norm of OMEGA2 -after CC_T2MO: ',RHO2N
         ENDIF
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND('After  T2MO: BF ')
            CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
         ENDIF
C
C---------------------------------------------------------------------
C        Write Gamma vector to disk if needed in response calculation.
C---------------------------------------------------------------------
C
         IF ( RSPIM ) THEN
C
            LUGAM = -1
            CALL GPOPEN(LUGAM,'CC_GAMIM','UNKNOWN',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
            REWIND(LUGAM)
            WRITE(LUGAM)(WORK(KGAMMA+I-1),I = 1,NGAMMA(ISYMOP))
            CALL GPCLOSE(LUGAM,'KEEP')
C
         ENDIF
C
C-------------------------------
C        Print the Gamma matrix.
C-------------------------------
C
         IF (IPRINT .GT. 120) THEN
            CALL AROUND('The Gamma matrix')
            DO 200 ISYM = 1,NSYM
               KOFF = KGAMMA + IGAMMA(ISYM,ISYM)
               CALL OUTPAK(WORK(KOFF),NMATIJ(ISYM),1,LUPRI)
  200       CONTINUE
C
            WRITE(LUPRI,*) 'Norm of gamma matrix: ',
     *              DDOT(NGAMMA(ISYMOP),WORK(KGAMMA),1,WORK(KGAMMA),1)
         END IF
C
C--------------------------------------------
C        Restore the CC amplitudes from disk.
C--------------------------------------------
C
         REWIND (LURHS1)
         READ (LURHS1)
         READ (LURHS1)
         READ (LURHS1) (T2AM(I), I = 1,NT2AM(ISYMOP))
C
      ENDIF
C
C---------------------------------------
C     Write out AO fock as intermediate.
C---------------------------------------
C
      IF ( RSPIM ) THEN
C
         LUFCK = -1
         CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         REWIND(LUFCK)
         WRITE(LUFCK)(WORK(KFOCK + I-1),I = 1,N2BST(ISYMOP))
         CALL GPCLOSE(LUFCK,'KEEP' )
C
         IF (IPRINT .GT.150) THEN
            CALL AROUND( 'Fock AO matrix written to disk' )
            CALL CC_PRFCKAO(WORK(KFOCK),1)
         ENDIF
C
      ENDIF
C
C------------------------------------------
C     Transform AO Fock matrix to MO basis.
C------------------------------------------
C
      TIMFCKMO = SECOND()
      CALL CC_FCKMO(WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND1),LWRK1,1,1,1)
      TIMFCKMO = SECOND() - TIMFCKMO
C
C---------------------
C     Reallocate T2TP.
C---------------------
C
      IF (DIRECT .AND. T2TCOR) THEN
C
         KT2AMT = KEND1
         KEND2  = KT2AMT + NT2SQ(1)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2. LT. 0) THEN
            CALL QUIT('Insufficient memory in CCSD_RHS')
         END IF
C
      ELSE
C
         KEND2 = KEND1
         LWRK2 = LWRK1
C
      END IF
C
C--------------------------------------------------------------
C     Add connected triples corrections to the vector function.
C--------------------------------------------------------------
C
      IF (CCSDT) THEN

         IF (NODDY_OMEGA) THEN
C          Unrelaxed noddy
C          CALL CC_FOPTRIPLES(OMEGA1,DUMMY,DUMMY,T1AM,T2AM,
C    *                        WORK(KLAMDP),WORK(KLAMDH),
C    *                        WORK(KEND2),LWRK2)
C          Original noddy part
C          Used to calculate Finite difference CC3
           CALL CCSD_TRIPLE(OMEGA1,OMEGA2,T1AM,T2AM,WORK(KFOCK),
     *                      WORK(KLAMDP),WORK(KLAMDH),WORK(KEND2),
     *                      LWRK2)
         ELSE
 
C          Normal triples part
           CALL CC3_OMEG(0.0D0,OMEGA1,OMEGA2,T1AM,ISYMTR,T2AM,ISYMTR,
     *                   WORK(KFOCK),WORK(KLAMDP),WORK(KLAMDH),
     *                   WORK(KEND2),LWRK2,LU3SRT,FN3SRT,LUDELD,
     *                   FNDELD,LUCKJD,FNCKJD,LUDKBC,FNDKBC,
     *                   LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2)
        
         END IF
C 
C
C----------------------------------------
C        Reconstruct full square of T2AM.
C----------------------------------------
C
         IF (LWRK1 .LT. NT2AMX) THEN
            CALL QUIT('Insufficient core in CCRHSN')
         ENDIF
C
         REWIND (LURHS1)
         READ (LURHS1)
         READ (LURHS1) (WORK(KEND1+I-1), I = 1,NT2AMX)
C
         CALL CC_T2SQ(WORK(KEND1),T2AM,1)
C
      ENDIF
C
      IF (IPRINT .GT. 51) THEN
         RHO1N = DDOT(NT1AMX,OMEGA1,1,OMEGA1,1)
         RHO2N = DDOT(NT2AMX,OMEGA2,1,OMEGA2,1)
         WRITE(LUPRI,*) 'Norm of OMEGA1 -after cc3_omeg: ',RHO1N
         WRITE(LUPRI,*) 'Norm of OMEGA2 -after cc3_omeg: ',RHO2N
      ENDIF
C
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('After  CC3_OMEG Omega is ')
         CALL CC_PRP(OMEGA1,OMEGA2,1,1,1)
      ENDIF
C
C---------------------
C     Reallocate T2TP.
C---------------------
C
      IF ((DIRECT .AND. T2TCOR) .OR. (CCSDT .AND. T2TCOR)) THEN
C
         KT2AMT = KEND1
         KEND2  = KT2AMT + NT2SQ(1)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2. LT. 0) THEN
            CALL QUIT('Insufficient memory in CCSD_RHS')
         END IF
C
      ELSE
C
         KEND2 = KEND1
         LWRK2 = LWRK1
C
      END IF
C
C----------------------
C     Recalculate T2TP.
C----------------------
C
      IF (T2TCOR) THEN
C
         JSYM = 1
         CALL DCOPY(NT2SQ(1),T2AM,1,WORK(KT2AMT),1)
         CALL CCSD_T2TP(WORK(KT2AMT),WORK(KEND2),LWRK2,JSYM)
C
      END IF
C----------------------
C     Calculate J-term.
C----------------------
C
      TIMJ     = SECOND()
      CALL CCRHS_J(OMEGA1,1,WORK(KFOCK))
      TIMJ     = SECOND() - TIMJ
C
C----------------------
C     Calculate A-term.
C----------------------
C
      IOPT = 1
      TIMA     = SECOND()
      IF (.NOT. CC2) THEN
         CALL CCRHS_A(OMEGA2,T2AM,WORK(KGAMMA),WORK(KEND2),LWRK2,
     *                ISYMTR,ISYMTR,IOPT)
      ENDIF
      TIMA     = SECOND() - TIMA
C
C------------------------------------------------------------------
C     Calculate E-term.
C     Write out the matrices if response calculation is to be done.
C------------------------------------------------------------------
C
      TIME     = SECOND()
      IF (CC2 .AND. (.NOT.RSPIM)) THEN
         IF (.NOT. NONHF) THEN
           ISIDE = 1
           CALL CC2_FCK(OMEGA2,T2AM,WORK(KEND2),LWRK2,ISYMTR,
     *                  WORK(KLAMDP),WORK(KLAMDH),ISIDE)
         ELSE
           ETRAN  = .FALSE.
           FCKCON = .TRUE.
           ISYMEI = ISYMTR
           CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),WORK(KLAMDH),
     *                     WORK(KFCKHF),WORK(KEND2),LWRK2,FCKCON,
     *                     ETRAN,ISYMEI)
           CALL CCRHS_E(OMEGA2,T2AM,WORK(KEMAT1),WORK(KEMAT2),
     *                  WORK(KEND2),LWRK2,ISYMTR,ISYMOP)
         END IF
      ENDIF
C
      IF (CCR12 .AND. (IANR12.EQ.2 .OR. IANR12.EQ.3)) THEN
          LRES = .FALSE.
          CALL GETTIM(T0,W0)
          CALL CCRHS_EPPP(OMEGA2,WORK(KEND2),LWRK2,APROXR12,LRES,
     &                    IDUMMY,CDUMMY,IDUMMY,ISYMTR)
          CALL GETTIM(T1,W1)
          TIMEPPPCPU = T1-T0
          TIMEPPPWAL = W1-W0
          TIMR12CPU = TIMR12CPU + (T1-T0)
          TIMR12WAL = TIMR12WAL + (W1-W0)
      END IF
c
      IF ((.NOT.CC2) .OR. RSPIM) THEN
C
         ETRAN  = .TRUE.
         FCKCON = .TRUE.
         ISYMEI = ISYMTR
         CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),WORK(KLAMDH),
     *                   WORK(KFOCK),WORK(KEND2),LWRK2,FCKCON,
     *                   ETRAN,ISYMEI)
C
         IF (CCR12) THEN
           KEMAT2P = KEND2
           KEND3   = KEMAT2P + NMATIJ(ISYMOP)
           LWRK3   = LWORK - KEND3
           IF (LWRK3.LT.0) CALL QUIT('Insufficient memory in CCSD_RHS')

           CALL GETTIM(T0,W0)
           CALL DZERO(WORK(KEMAT2P),NMATIJ(ISYMOP))
           CALL CCRHS_EINTP(WORK(KEMAT2P),WORK(KLAMDP),
     &                      WORK(KEND3),LWRK3,0,1,CDUMMY,IDUMMY,
     &                      IDUMMY,CDUMMY,IDUMMY) 
           CALL GETTIM(T1,W1)
           IF (IPRINT .GT. 9) THEN
             WRITE(LUPRI,*)'Time used for CCRHS_EINTP cpu:', T1-T0
             WRITE(LUPRI,*)'Time used for CCRHS_EINTP wall:',W1-W0
           END IF
           TIMR12CPU = TIMR12CPU + (T1-T0)
           TIMR12WAL = TIMR12WAL + (W1-W0)
C         
           ! add R12 contribution to usual E_ij intermediate
           CALL DAXPY(NMATIJ(ISYMOP),ONE,WORK(KEMAT2P),1,WORK(KEMAT2),1)
C
           IF (IANR12.NE.1) THEN
             ! add R12 contribution to usual E_ab intermediate
             IOPTE = 1
             CALL GETTIM(T0,W0)
             CALL CCRHS_HP(WORK(KEMAT1),WORK(KLAMDH),1,WORK(KLAMDH),1,
     &                     WORK(KEND2),LWRK2,0,1,CDUMMY,IDUMMY,IDUMMY,
     &                     IOPTE)
             CALL GETTIM(T1,W1)
             TIMR12CPU = TIMR12CPU + (T1-T0)
             TIMR12WAL = TIMR12WAL + (W1-W0)
           END IF
C
           IF (IPRINT .GT. 9) THEN
             WRITE(LUPRI,*)'Time used for R12 part in CCSD_RHS cpu:', 
     &                      TIMR12CPU
             WRITE(LUPRI,*)'Time used for R12 part in CCSD_RHS wall:', 
     &                      TIMR12WAL
           END IF
         END IF
C
         IF ( RSPIM ) THEN
C
            LUE1 = -1
            CALL GPOPEN(LUE1,'CC_E1IM','UNKNOWN',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
            REWIND(LUE1)
            WRITE(LUE1)(WORK(KEMAT1+ I-1),I = 1,NMATAB(ISYMOP))
            CALL GPCLOSE(LUE1,'KEEP' )
C
            LUE2 = -1
            CALL GPOPEN(LUE2,'CC_E2IM','UNKNOWN',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
            REWIND(LUE2)
            WRITE(LUE2)(WORK(KEMAT2+ I-1),I = 1,NMATIJ(ISYMOP))
            CALL GPCLOSE(LUE2,'KEEP' )
C
            IF (CCR12) THEN
             LUE2P = -1
             CALL GPOPEN(LUE2P,'CC_E2PIM','UNKNOWN',' ','UNFORMATTED',
     *                   IDUMMY,.FALSE.)
             REWIND(LUE2P)
             WRITE(LUE2P)(WORK(KEMAT2P+ I-1),I = 1,NMATIJ(ISYMOP))
             CALL GPCLOSE(LUE2P,'KEEP' )
            END IF
C
            IF (IPRINT.GT.40) THEN
               CALL AROUND( 'E-intermediates written to disk ')
               CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYMOP,1)
            ENDIF
            IF (DEBUG) THEN
               XNORM1 = DDOT(NMATAB(1),WORK(KEMAT1),1,WORK(KEMAT1),1)
               XNORM2 = DDOT(NMATIJ(1),WORK(KEMAT2),1,WORK(KEMAT2),1)
               WRITE(LUPRI,*) 'Norm of E1 intermediate:',XNORM1
               WRITE(LUPRI,*) 'Norm of E2 intermediate:',XNORM2
               IF (CCR12) THEN
                 XNORM2=DDOT(NMATIJ(1),WORK(KEMAT2P),1,WORK(KEMAT2P),1)
                 WRITE(LUPRI,*) 'Norm of E2P intermediate:',XNORM2
               END IF
            END IF
C
         ENDIF
C
         IF (.NOT.CC2) THEN
C
            CALL CCRHS_E(OMEGA2,T2AM,WORK(KEMAT1),WORK(KEMAT2),
     *                   WORK(KEND2),LWRK2,ISYMTR,ISYMOP)
C
         ENDIF
C
      ENDIF
      TIME     = SECOND() - TIME
C
C--------------------------------------
C     If (DUMPCD) calculate the C-term.
C--------------------------------------
C
      IF (DUMPCD .AND. (.NOT. CC2)) THEN
C
         ISYVEC = 1
         ISYCIM = 1
         IOPT   = 1
         IVECNR = 1
C
         TIMCIO = SECOND()
         IF (T2TCOR) THEN
            CALL CCRHS_CIO(OMEGA2,WORK(KT2AMT),WORK(KLAMDH),
     *                     WORK(KEND2),LWRK2,ISYVEC,ISYCIM,
     *                     LUC,CFIL,IVECNR,IOPT)
         ELSE
            ISYM = 1
            CALL CCSD_T2TP(T2AM,WORK(KEND2),LWRK2,ISYM)
            CALL CCRHS_CIO(OMEGA2,T2AM,WORK(KLAMDH),WORK(KEND2),
     *                     LWRK2,ISYVEC,ISYCIM,LUC,CFIL,IVECNR,IOPT)
            CALL CCSD_T2TP(T2AM,WORK(KEND2),LWRK2,ISYM)
         ENDIF
C
         TIMCIO  = SECOND() - TIMCIO
C
      ENDIF
C
C------------------------------
C     Transform T2 to 2T2 - T2.
C------------------------------
C
      DTIME    = SECOND()
      IF (T2TCOR) THEN
         CALL DSCAL(NT2SQ(1),TWO,T2AM,1)
         CALL DAXPY(NT2SQ(1),-ONE,WORK(KT2AMT),1,T2AM,1)
      ELSE
         ISYM = 1
         CALL CCRHS_T2TR(T2AM,WORK(KEND2),LWRK2,ISYM)
      END IF
      DTIME    = SECOND() - DTIME
      TIMT2TR  = TIMT2TR + DTIME
C
C--------------------------------------
C     If (DUMPCD) calculate the D-term.
C--------------------------------------
C
      IF (DUMPCD .AND. (.NOT. CC2)) THEN
C
         ISYDIM = 1
         ISYVEC = 1
         IOPT = 1
         IVECNR = 1
C
         TIMDIO = SECOND()
         CALL CCRHS_DIO(OMEGA2,T2AM,WORK(KLAMDH),WORK(KEND2),LWRK2,
     *                  ISYVEC,ISYDIM,LUD,DFIL,IVECNR,IOPT)
         TIMDIO  = SECOND() - TIMDIO
      END IF
C
C----------------------
C     Calculate I-term.
C----------------------
C
      TIMI     = SECOND()
      CALL CCRHS_I(OMEGA1,T2AM,WORK(KFOCK),WORK(KEND2),LWRK2,ISYMTR,1)
      TIMI     = SECOND() - TIMI

C-----------------------------------------------------------------
C     Add the remaining CCSDR12 C-, D- and E-contributions 
C     to Omega2
C-----------------------------------------------------------------
C
      IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
C
         !calculate t(bj,p'k) amplitudes:
         CALL CC_R12MKTBJPK(T2AM,WORK(KEND2),LWRK2)
C
         !Read CMO-Matrix incl. aux.-orbitals:
         KCMOX = KEND2
         KEND3 = KCMOX + NLAMDX(1)
         LWRK3 = LWORK - KEND3
         IF (LWRK3.LT.0) CALL QUIT('Insuff. memory for CCRHS_CIO')
         CALL CC_R12CMO(WORK(KCMOX),WORK(KEND3),LWRK3)
C
         !calculate E1P_(ap') intermediate:
         CALL CCRHS_E1PIM(WORK(KE1PIM),WORK(KCMOX),ILAMDX,WORK(KLAMDH),
     &                    WORK(KEND3),LWRK3)
C
CTesT
C        WRITE(LUPRI,*) 'E1PIM after transformation to MO:'
C        KOFF = 0
C        DO ISYM = 1,NSYM
C          CALL OUTPUT(WORK(KE1PIM+KOFF),
C    &                 1,NVIR(ISYM),1,NORB2(ISYM),
C    &                 NVIR(ISYM),NORB2(ISYM),1, LUPRI)
C          KOFF = KOFF + NVIR(ISYM)*NORB2(ISYM)
C        END DO
C        WRITE(LUPRI,*) 'Norm^2: ',
C    &                  DDOT(KOFF,WORK(KE1PIM),1,WORK(KE1PIM),1)
C        CALL FLSHFO(LUPRI)
C
C        WRITE(LUPRI,*) "OMEGA2 before C', D', E' contr.:"
C        WRITE(LUPRI,*) "Norm^2: ", DDOT(NT2AMX,OMEGA2,1,OMEGA2,1)
C        DO ISYM = 1,NSYM
C           WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
C           KOFF = IT2AM(ISYM,ISYM) + 1
C           CALL OUTPAK(OMEGA2(KOFF),NT1AM(ISYM),1,LUPRI)
C        END DO
CTesT
C
         ISYVEC = 1
         ISYCIM = 1
         ISYDIM = 1
         IVECNR = 1
         IOPT = 1
         IOPTB = 0
         IOPTE = 1
         CALL CCRHS_CIO2(OMEGA2,T2AM,WORK(KCMOX),
     *                   WORK(KEND3),LWRK3,ISYVEC,ISYCIM,
     *                   LUCP,CPFIL,IVECNR,IOPT,IOPTB,IDUMMY,
     *                   DUMMY,IDUMMY,DUMMY,IOPTE,WORK(KE1PIM),.TRUE.)
 
         CALL CCRHS_DIO2(OMEGA2,T2AM,WORK(KCMOX),
     *                   WORK(KEND3),LWRK3,ISYVEC,ISYDIM,
     *                   LUDP,DPFIL,IDUMMY,DUMMY,IVECNR,IOPT,
     *                   IOPTB,IDUMMY,DUMMY,IDUMMY,DUMMY,
     *                   IOPTE,WORK(KE1PIM),.TRUE.)
C
CTesT
C        WRITE(LUPRI,*) "OMEGA2 after C', D', E' contr.:"
C        WRITE(LUPRI,*) "Norm^2: ", DDOT(NT2AMX,OMEGA2,1,OMEGA2,1)
C        DO ISYM = 1,NSYM
C           WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
C           KOFF = IT2AM(ISYM,ISYM) + 1
C           CALL OUTPAK(OMEGA2(KOFF),NT1AM(ISYM),1,LUPRI)
C        END DO
C        STOP
CTesT         
      END IF
C
C-----------------------------------------------------------------
C     Calculate the C and D contributions to the R12 result vector
C-----------------------------------------------------------------
C
      IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
         ! save vector function on disk, since memory is needed  
         LUOME1 = -1
         CALL GPOPEN(LUOME1,'CCOME1','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
cch
         write(lupri,*) 'Norm of OMEGA2:',ddot(nt2amx,OMEGA2,1,OMEGA2,1)
cch
         REWIND(LUOME1)
         WRITE(LUOME1) (OMEGA2(I), I = 1,NT2AMX)
         CALL GPCLOSE(LUOME1,'KEEP')

         ! read cluster amlitudes in packed form into memory
         REWIND(LURHS1)
         READ(LURHS1)
         READ(LURHS1) (OMEGA2(I), I = 1,NT2AMX)

         !Read CMO-Matrix incl. aux.-orbitals:
         KCMOX = KEND2
         KEND3 = KCMOX + NLAMDX(1)
         LWRK3 = LWORK - KEND3
         IF (LWRK3.LT.0) CALL QUIT('Insuff. memory for CCRHS_CIO')
         CALL CC_R12CMO(WORK(KCMOX),WORK(KEND3),LWRK3)

         CALL CCSDR12CD(CCSDR12,
     &                  T2AM,1,OMEGA2,1,1,
     &                  FNIADJ,LUIADJ,FNIJDA,LUIJDA,IT2DEL,
     &                  WORK(KLAMDH),1,
     &                  WORK(KCMOX),ILAMDX,
     &                  WORK(KEND3),LWRK3)

         ! restore vector function
         LUOME1 = -1
         CALL GPOPEN(LUOME1,'CCOME1','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
         REWIND(LUOME1)
         READ(LUOME1) (OMEGA2(I), I = 1,NT2AMX)
         CALL GPCLOSE(LUOME1,'DELETE')
cch
         write(lupri,*) 'Norm of OMEGA2:',ddot(nt2amx,OMEGA2,1,OMEGA2,1)
cch       
      END IF
C
C------------------------
C     Scale final result.
C------------------------
C
C     CALL DSCAL(NT1AM,TWO,OMEGA1,1)
C     CALL DSCAL(NT2IND,TWO,OMEGA2,1)
C
      IF (IPRINT .GT. 25) THEN
         CALL AROUND('END OF CCRHS:OMEGA 1')
         DO 300 ISYM = 1,NSYM
            WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
            KOFF = IT1AM(ISYM,ISYM) + 1
            CALL OUTPUT(OMEGA1(KOFF),1,NVIR(ISYM),1,NRHF(ISYM),
     *                  NVIR(ISYM),NRHF(ISYM),1,LUPRI)
  300    CONTINUE
         WRITE(LUPRI,*)
         CALL AROUND('END OF CCRHS:OMEGA 2')
         DO 310 ISYM = 1,NSYM
            WRITE(LUPRI,*) 'Symmetry block number : ',ISYM
            KOFF = IT2AM(ISYM,ISYM) + 1
            CALL OUTPAK(OMEGA2(KOFF),NT1AM(ISYM),1,LUPRI)
  310    CONTINUE
      ENDIF
      TIMALL  = SECOND() - TIMALL
      IF ( IPRINT .GT. 2) THEN
         WRITE(LUPRI,9999) 'RHS - TOTAL', TIMALL
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,9999) 'CCRHS_A    ', TIMA
         WRITE(LUPRI,9999) 'CCRHS_B    ', TIMB
         WRITE(LUPRI,9999) 'CCRHS_BF   ', TIMBF
         WRITE(LUPRI,9999) 'CCRHS_C    ', TIMC
         WRITE(LUPRI,9999) 'CCRHS_CIO  ', TIMCIO
         WRITE(LUPRI,9999) 'CCRHS_C-tot', TIMCIO + TIMC
         WRITE(LUPRI,9999) 'CCRHS_D    ', TIMD
         WRITE(LUPRI,9999) 'CCRHS_DIO  ', TIMDIO
         WRITE(LUPRI,9999) 'CCRHS_D-tot', TIMDIO + TIMD
         WRITE(LUPRI,9999) 'CCRHS_E    ', TIME
         WRITE(LUPRI,9999) 'CCRHS_EI   ', TIMEI
         WRITE(LUPRI,9999) 'CCRHS_E-tot', TIMEI + TIME
         WRITE(LUPRI,9999) 'CCRHS_F    ', TIMF
         WRITE(LUPRI,9999) 'CCRHS_G    ', TIMG
         WRITE(LUPRI,9999) 'CCRHS_H    ', TIMH
         WRITE(LUPRI,9999) 'CCRHS_I    ', TIMI
         WRITE(LUPRI,9999) 'CCRHS_J    ', TIMJ
         WRITE(LUPRI,9999) 'CCRHS_GAM  ', TIMGAM
         WRITE(LUPRI,9999) 'CCRHS_LAM  ', TIMLAM
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
         IF (CCR12.AND.(IANR12.EQ.1)) THEN
           WRITE(LUPRI,9999) 'CCRHS_FP   ', TIMFP
           WRITE(LUPRI,9999)'CC_MOFCON cpu:', TIMMOFCPU 
           WRITE(LUPRI,9999)'CC_MOFCON wall:', TIMMOFWAL 
           WRITE(LUPRI,9999)'CC_MOFCONR12 cpu:', TIMMOFR12CPU 
           WRITE(LUPRI,9999)'CC_MOFCONR12 wall:', TIMMOFR12WAL 
           WRITE(LUPRI,9999)'CCRHS_EP cpu:', TIMEPCPU
           WRITE(LUPRI,9999)'CCRHS_EP wall:',TIMEPWAL
           WRITE(LUPRI,9999)'CCRHS_GP cpu:', TIMGPCPU
           WRITE(LUPRI,9999)'CCRHS_GP wall:', TIMGPWAL
           WRITE(LUPRI,9999) 'INTEG. R12 ', TIMINTR12
           WRITE(LUPRI,9999) 'RDAO   R12 ', TIMRDAOR12
         ELSE IF (CCR12.AND.(IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
           WRITE(LUPRI,9999)'CC_MOFCON cpu:', TIMMOFCPU 
           WRITE(LUPRI,9999)'CC_MOFCON wall:', TIMMOFWAL 
           WRITE(LUPRI,9999)'CC_R12INTF2 cpu:', TIMINTF2CPU 
           WRITE(LUPRI,9999)'CC_R12INTF2 wall:', TIMINTF2WAL 
           WRITE(LUPRI,9999)'CC_MOFCONR12 cpu:', TIMMOFR12CPU 
           WRITE(LUPRI,9999)'CC_MOFCONR12 wall:', TIMMOFR12WAL 
           WRITE(LUPRI,9999)'CCRHS_EPP cpu:', TIMEPPCPU 
           WRITE(LUPRI,9999)'CCRHS_EPP wall:', TIMEPPWAL 
           WRITE(LUPRI,9999)'CCRHS_EPPP cpu:', TIMEPPPCPU 
           WRITE(LUPRI,9999)'CCRHS_EPPP wall:', TIMEPPPWAL 
           WRITE(LUPRI,9999)'CCRHS_HP cpu:', TIMHPCPU 
           WRITE(LUPRI,9999)'CCRHS_HP wall:', TIMHPWAL 
           WRITE(LUPRI,9999)'CCRHS_IP cpu:', TIMIPCPU 
           WRITE(LUPRI,9999)'CCRHS_IP wall:', TIMIPWAL 
           WRITE(LUPRI,9999)'CCRHS_EP cpu:', TIMEPCPU
           WRITE(LUPRI,9999)'CCRHS_EP wall:',TIMEPWAL
           WRITE(LUPRI,9999)'CCRHS_GP cpu:', TIMGPCPU
           WRITE(LUPRI,9999)'CCRHS_GP wall:', TIMGPWAL
           WRITE(LUPRI,9999)'R12 cpu:', TIMR12CPU
           WRITE(LUPRI,9999)'R12 wall:', TIMR12WAL
           WRITE(LUPRI,9999) 'INTEG. R12 ', TIMINTR12
           WRITE(LUPRI,9999) 'RDAO   R12 ', TIMRDAOR12
         END IF
      ENDIF
9999  FORMAT(7x,'Time used in',2x,A12,2x,': ',f10.2,' seconds')
C
C-----------------------------------------
C     Restore the CC amplitudes from disk.
C-----------------------------------------
C
      REWIND (LURHS1)
      READ(LURHS1) (T1AM(I), I = 1,NT1AMX)
      READ(LURHS1) (T2AM(I), I = 1,NT2AMX)
      CALL GPCLOSE(LURHS1,'DELETE')
C
C-----------------
C     Close files.
C-----------------
C
      IF (DUMPCD) THEN
         CALL WCLOSE2(LUC,CFIL,'KEEP')
         CALL WCLOSE2(LUD,DFIL,'KEEP')
      END IF
C
      IF (CCSDT) THEN
         CALL WCLOSE2(LU3SRT,FN3SRT,'KEEP')
         CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
         CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
         CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
         CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
         CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
         CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      ENDIF

      IF (CCSDR12 .AND. (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
         CALL WCLOSE2(LUIADJ,FNIADJ,'KEEP')
         CALL WCLOSE2(LUIJDA,FNIJDA,'KEEP')
         CALL WCLOSE2(LUCP,CPFIL,'KEEP')
         CALL WCLOSE2(LUDP,DPFIL,'KEEP')
      END IF
C
C-----------------------
C     Restore CC1B flag.
C-----------------------
C
      CC1B = CC1BSA
C
      CALL QEXIT('CCRHSN')
C
      RETURN
      END
C  /* Deck ccrhs_e */
      SUBROUTINE CCRHS_E(OMEGA2,T2AM,EMAT1,EMAT2,WORK,LWORK,
     *                  ISYMTR,ISYMIM)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 20-Jan-1994
C     Symmetry 3-aug
C     Contraction of EI intermediates with double excitaion amplitudes.
C     It is assumed that the fock matrix is included. OC 13-1-1995
C
C     Purpose: Calculate E-terms
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DIMENSION EMAT1(*),EMAT2(*)
      DIMENSION T2AM(*),OMEGA2(*)
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C--------------------------------------------------------------
C     Contract and accumulate the first intermediate in OMEGA2.
C--------------------------------------------------------------
C
      ISYAIBJ = MULD2H(ISYMTR,ISYMIM)
C
      DO 300 ISYMAI = 1,NSYM
C
         ISYMCJ = MULD2H(ISYMAI,ISYMTR)
         ISYMBJ = MULD2H(ISYMAI,ISYAIBJ)
C
         IF (LWORK .LT. NT1AM(ISYMBJ)) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_E1')
         END IF
C
         DO 310 NAI = 1,NT1AM(ISYMAI)
C
            CALL DZERO(WORK,NT1AM(ISYMBJ))
C
            DO 320 ISYMJ = 1,NSYM
C
               ISYMC  = MULD2H(ISYMJ,ISYMCJ)
               ISYMB  = MULD2H(ISYMJ,ISYMBJ)
C
               NVIRB = MAX(NVIR(ISYMB),1)
               NVIRC = MAX(NVIR(ISYMC),1)
C
               KOFF1 = IMATAB(ISYMB,ISYMC) + 1
               KOFF2 = IT2SQ(ISYMCJ,ISYMAI) + NT1AM(ISYMCJ)*(NAI - 1)
     *                  + IT1AM(ISYMC,ISYMJ) + 1
               KOFF3 = IT1AM(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMJ),
     *                    NVIR(ISYMC),ONE,EMAT1(KOFF1),NVIRB,
     *                    T2AM(KOFF2),NVIRC,
     *                    ONE,WORK(KOFF3),NVIRB)
  320       CONTINUE
C
            IF (ISYMAI .EQ. ISYMBJ ) THEN
C
               WORK(NAI) = TWO*WORK(NAI)
               DO 330 NBJ = 1,NT1AM(ISYMBJ)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(NBJ)
  330          CONTINUE
C
            ENDIF
C
            IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO 340 NBJ = 1,NT1AM(ISYMBJ)
                 NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(NBJ)
  340          CONTINUE
C
            ENDIF
C
            IF (ISYMBJ .LT. ISYMAI) THEN
C
               DO 350 NBJ = 1,NT1AM(ISYMBJ)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(NBJ)
  350          CONTINUE
C
            ENDIF
C
  310    CONTINUE
  300 CONTINUE
C
C-----------------------------------------------------
C     Contract and accumulate the second intermediate.
C-----------------------------------------------------
C
C
      DO 400 ISYMAI = 1,NSYM
C
         ISYMBK = MULD2H(ISYMAI,ISYMTR)
         ISYMBJ = MULD2H(ISYMAI,ISYAIBJ)
C
         IF (LWORK .LT. NT1AM(ISYMBJ)) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_E1')
         END IF
C
         DO 410 NAI = 1,NT1AM(ISYMAI)
C
            CALL DZERO(WORK,NT1AM(ISYMBJ))
C
            DO 420 ISYMB = 1,NSYM
C
               ISYMJ  = MULD2H(ISYMB,ISYMBJ)
               ISYMK  = MULD2H(ISYMJ,ISYMIM)
C
               NVIRB = MAX(NVIR(ISYMB),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KOFF1 = IT2SQ(ISYMBK,ISYMAI) + NT1AM(ISYMBK)*(NAI - 1)
     *               + IT1AM(ISYMB,ISYMK) + 1
               KOFF2 = IMATIJ(ISYMK,ISYMJ) + 1
               KOFF3 = IT1AM(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMJ),
     *                    NRHF(ISYMK),ONE,T2AM(KOFF1),NVIRB,
     *                    EMAT2(KOFF2),NRHFK,
     *                    ONE,WORK(KOFF3),NVIRB)
  420       CONTINUE
C
C
            IF (ISYMAI .EQ. ISYMBJ ) THEN
C
               WORK(NAI) = TWO*WORK(NAI)
C
               DO 430 NBJ = 1,NT1AM(ISYMBJ)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(NBJ)
  430          CONTINUE
C
            ENDIF
C
            IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO 440 NBJ = 1,NT1AM(ISYMBJ)
                 NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(NBJ)
  440          CONTINUE
C
            ENDIF
C
            IF (ISYMBJ .LT. ISYMAI) THEN
C
               DO 450 NBJ = 1,NT1AM(ISYMBJ)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) - WORK(NBJ)
  450          CONTINUE
C
            ENDIF
C
  410    CONTINUE
  400 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_i */
      SUBROUTINE CCRHS_I(OMEGA1,T2AM,FOCK,WORK,LWORK,ISYMT2,ISYMCK)
C
C     Written by Henrik Koch & Ove Christiansen 20-Jan-1994
C
C     Purpose: Calculate I-term.
C
#include "implicit.h"
      PARAMETER(ONE=1.0D0)
      DIMENSION OMEGA1(*),WORK(*)
      DIMENSION T2AM(*),FOCK(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      ISYMAI = MULD2H(ISYMT2,ISYMCK)
C
      KSCR1 = 1
      KEND1 = KSCR1 + NT1AM(ISYMCK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in CCRHS_I')
      END IF
C
      DO 110 ISYMK = 1,NSYM
C
         ISYMC = MULD2H(ISYMK,ISYMCK)
C
         NVIRC = MAX(NVIR(ISYMC),1)
C
         DO 120 K = 1,NRHF(ISYMK)
C
            KOFF1 = IFCVIR(ISYMK,ISYMC) + K
            KOFF2 = KSCR1 + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
C
            CALL DCOPY(NVIR(ISYMC),FOCK(KOFF1),NORB(ISYMK),
     *                 WORK(KOFF2),1)
C
  120    CONTINUE
C
  110 CONTINUE
C
      NTOTAI = MAX(NT1AM(ISYMAI),1)
C
      KOFF3 = IT2SQ(ISYMAI,ISYMCK) + 1
C
      CALL DGEMV('N',NT1AM(ISYMAI),NT1AM(ISYMCK),ONE,T2AM(KOFF3),
     *           NTOTAI,WORK(KSCR1),1,ONE,OMEGA1,1)
C
      RETURN
      END
      SUBROUTINE CCRHS_A(OMEGA2,T2AM,GAMMA,WORK,LWORK,ISYGAM,ISYVEC,
     *                   IOPT)
C
C     Written by Henrik Koch & Ove Christiansen 20-Jan-1994
C
C     Generalised to non. total sym gamma (isygam) og non. tot. sym
C     double excitation vector (isyvec) Ove Christiansen 29-7-1995
C
C     Generalised to handle left hand side contribution (IOPT 2) as
C     well as usual contributions (IOPT 1) by Asger Halkier 22/11-95.
C
C     Purpose: Calculate A-term.
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0, ONE=1.0D0)
      DIMENSION OMEGA2(*),GAMMA(*),T2AM(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      ISAIBJ = MULD2H(ISYGAM,ISYVEC)
C
      DO 100 ISYMLJ = 1,NSYM
C
         ISYMKI = MULD2H(ISYMLJ,ISYGAM)
C
         KSCR1 = 1
         KEND1 = KSCR1 + NMATIJ(ISYMKI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_A')
         END IF
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYML = MULD2H(ISYMJ,ISYMLJ)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 L = 1,NRHF(ISYML)
C
                  IF (IOPT .EQ. 1) THEN
C
                     NLJ = IMATIJ(ISYML,ISYMJ)
     *                   + NRHF(ISYML)*(J - 1) + L
C
                  ELSE IF (IOPT .EQ. 2) THEN
C
                     NLJ = IMATIJ(ISYMJ,ISYML)
     *                   + NRHF(ISYMJ)*(L - 1) + J
C
                  ENDIF
C
                  DO 140 ISYMK = 1,NSYM
C
                     ISYMI = MULD2H(ISYMK,ISYMKI)
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        DO 160 K = 1,NRHF(ISYMK)
C
                           IF (IOPT .EQ. 1) THEN
C
                              NKI = IMATIJ(ISYMK,ISYMI)
     *                            + NRHF(ISYMK)*(I - 1) + K
C
                           ELSE IF (IOPT .EQ. 2) THEN
C
                              NKI = IMATIJ(ISYMI,ISYMK)
     *                            + NRHF(ISYMI)*(K - 1) + I
C
                           ENDIF
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + INDEX(NKI,NLJ)
                           ELSE
                              IF (ISYMKI .LT. ISYMLJ) THEN
                                 NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                                 + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              ELSE
                                 NKILJ = IGAMMA(ISYMLJ,ISYMKI)
     *                                 + NMATIJ(ISYMLJ)*(NKI - 1) + NLJ
                              ENDIF
                           ENDIF
C
                           NSTO = IMATIJ(ISYMK,ISYMI)
     *                          + NRHF(ISYMK)*(I - 1) + K
C
                           WORK(KSCR1 + NSTO - 1) = GAMMA(NKILJ)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
C
                  DO 170 ISYMB = 1,NSYM
C
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
                     ISYMAI = MULD2H(ISYMBJ,ISAIBJ)
                     ISYMBL = MULD2H(ISYMB,ISYML)
                     ISYMAK = MULD2H(ISYVEC,ISYMBL)
C
                     KSCR2 = KEND1
                     KEND2 = KSCR2 + NT1AM(ISYMAI)
                     LWRK2 = LWORK - KEND2
C
                     IF (LWRK2 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCRHS_A')
                     END IF
C
                     IF (ISYMAI .GT. ISYMBJ) GOTO 170
C
                     DO 180 B = 1,NVIR(ISYMB)
C
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
                        NBL = IT1AM(ISYMB,ISYML)
     *                      + NVIR(ISYMB)*(L - 1) + B
C
                        CALL DZERO(WORK(KSCR2),NT1AM(ISYMAI))
C
                        DO 190 ISYMI = 1,NSYM
C
                           ISYMK = MULD2H(ISYMI,ISYMKI)
                           ISYMA = MULD2H(ISYMK,ISYMAK)
C
                           NVIRA = MAX(NVIR(ISYMA),1)
                           NRHFK = MAX(NRHF(ISYMK),1)
C
                           KOFF1 = IT2SQ(ISYMAK,ISYMBL)
     *                           + NT1AM(ISYMAK)*(NBL - 1)
     *                           + IT1AM(ISYMA,ISYMK) + 1
                           KOFF2 = KSCR1 + IMATIJ(ISYMK,ISYMI)
                           KOFF3 = KSCR2 + IT1AM(ISYMA,ISYMI)
C
                           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                                NRHF(ISYMK),ONE,T2AM(KOFF1),
     *                                NVIRA,WORK(KOFF2),NRHFK,ZERO,
     *                                WORK(KOFF3),NVIRA)
C
  190                   CONTINUE
C
                        IF (ISYMAI .EQ. ISYMBJ) THEN
                           NTOT = NBJ
                        ELSE
                           NTOT = NT1AM(ISYMAI)
                        ENDIF
C
                        DO 200 NAI = 1,NTOT
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ENDIF
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                                   + WORK(KSCR2 + NAI - 1)
C
  200                   CONTINUE
C
  180                CONTINUE
  170             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_j */
      SUBROUTINE CCRHS_J(OMEGA1,ISYM,FOCK)
C
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C
C     Purpose: Calculate J-term.
C
#include "implicit.h"
      DIMENSION FOCK(*),OMEGA1(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C--------------------
C     Calculate term.
C--------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYM)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            DO 120 A = 1,NVIR(ISYMA)
C
               KOFF1 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
               KOFF2 = IFCRHF(ISYMA,ISYMI) + NORB(ISYMA)*(I - 1)
     *               + NRHF(ISYMA) + A
C
               OMEGA1(KOFF1) = OMEGA1(KOFF1) + FOCK(KOFF2)
C
  120       CONTINUE
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_fckmo */
      SUBROUTINE CC_FCKMO(FOCK,XLAMDP,XLAMDH,WORK,LWORK,ISYFAO,
     *                    ISYMPA,ISYMHO)
C
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 11-July-1994
C
C     Ove Christiansen 14-7-1994 generalized to 1. non-tot.symmetric FockAO
C                                               2. non-tot symmetric LAM.
C
C     Filip Pawlowski  03-Jan-2007: introduced zeroing of fock array as
C     a symmetry bug fix
C
C              ISYFAO is the symmtry of the AO fock matrix.
C              isymp(isymh) is the symmetry of lamp(lamh)
C
C     Purpose: Calculate MO Fock Matrix.
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, ZERO = 0.0D0)
      DIMENSION FOCK(*),XLAMDP(*),XLAMDH(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      ISYML  = MULD2H(ISYMPA,ISYMHO)
      ISYFMO = MULD2H(ISYML,ISYFAO)
C
      KFOCK = 1
      KEND1 = KFOCK + N2BST(ISYFMO)
      LWRK1 = LWORK - KEND1
C
      KOFF1 = 1
      KOFF2 = KFOCK
C
      DO 100 ISYMQ = 1,NSYM
C
         ISYMB = MULD2H(ISYMQ,ISYMHO)
         ISYMA = MULD2H(ISYMB,ISYFAO)
         ISYMP = MULD2H(ISYFMO,ISYMQ)
C
C-----------------------------------------
C        Dynamic allocation of work space.
C-----------------------------------------
C
         NBQ = NBAS(ISYMB)*NORB(ISYMQ)
         NAP = NBAS(ISYMA)*NORB(ISYMP)
         KLAMDA = KEND1
         KSCR1  = KLAMDA + MAX(NAP,NBQ)
         KEND2  = KSCR1  + NBAS(ISYMA)*NORB(ISYMQ)
         LWRK2  = LWORK  - KEND2
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient space in CC_FCKMO')
         ENDIF
C
C-----------------------------------------
C        Copy transformation coefficients.
C-----------------------------------------
C
         NTOTR = NBAS(ISYMB)*NRHF(ISYMQ)
         KOFF  = IGLMRH(ISYMB,ISYMQ) + 1
         CALL DCOPY(NTOTR,XLAMDH(KOFF),1,WORK(KLAMDA),1)
C
         NTOTV = NBAS(ISYMB)*NVIR(ISYMQ)
         KOFF  = IGLMVI(ISYMB,ISYMQ) + 1
         CALL DCOPY(NTOTV,XLAMDH(KOFF),1,WORK(KLAMDA+NTOTR),1)
C
C----------------------------------------
C        Do first partial transformation.
C----------------------------------------
C
         NBASB = MAX(NBAS(ISYMB),1)
         NBASA = MAX(NBAS(ISYMA),1)
C
         KOFF1 = IAODIS(ISYMA,ISYMB) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMA),NORB(ISYMQ),NBAS(ISYMB),
     *              ONE,FOCK(KOFF1),NBASA,WORK(KLAMDA),NBASB,
     *              ZERO,WORK(KSCR1),NBASA)
C
C-----------------------------------------
C        Copy transformation coefficients.
C-----------------------------------------
C
         NTOTR = NBAS(ISYMA)*NRHF(ISYMP)
         KOFF  = IGLMRH(ISYMA,ISYMP) + 1
         CALL DCOPY(NTOTR,XLAMDP(KOFF),1,WORK(KLAMDA),1)
C
         NTOTV = NBAS(ISYMA)*NVIR(ISYMP)
         KOFF  = IGLMVI(ISYMA,ISYMP) + 1
         CALL DCOPY(NTOTV,XLAMDP(KOFF),1,WORK(KLAMDA+NTOTR),1)
C
C-----------------------------------------
C        Do second partial transformation.
C-----------------------------------------
C
         NBASA = MAX(NBAS(ISYMA),1)
         NORBP = MAX(NORB(ISYMP),1)
C
         CALL DGEMM('T','N',NORB(ISYMP),NORB(ISYMQ),NBAS(ISYMA),ONE,
     *              WORK(KLAMDA),NBASA,WORK(KSCR1),NBASA,ZERO,
     *              WORK(KOFF2),NORBP)
C
         KOFF2 = KOFF2 + NORB(ISYMP)*NORB(ISYMQ)
C
  100 CONTINUE
C
C-----------------------------------------------------
C     Reorder the Fock matrix in occupied and virtual.
C-----------------------------------------------------
C
      KOFF1 = KFOCK
      KOFF2 = 1
      KOFF3 = NLRHFR(ISYFMO)  + 1

      CALL DZERO(FOCK,N2BST(ISYFMO))

      DO 110 ISYMQ = 1,NSYM
C
         ISYMP = MULD2H(ISYMQ,ISYFMO)
C
         NTOTR = NORB(ISYMP)*NRHF(ISYMQ)
         CALL DCOPY(NTOTR,WORK(KOFF1),1,FOCK(KOFF2),1)
C
         NTOTV = NORB(ISYMP)*NVIR(ISYMQ)
         CALL DCOPY(NTOTV,WORK(KOFF1+NTOTR),1,FOCK(KOFF3),1)
C
         KOFF1 = KOFF1 + NORB(ISYMP)*NORB(ISYMQ)
         KOFF2 = KOFF2 + NORB(ISYMP)*NRHF(ISYMQ)
         KOFF3 = KOFF3 + NORB(ISYMP)*NVIR(ISYMQ)
C
  110 CONTINUE
C
      END
C  /* Deck ccrhs_h */
      SUBROUTINE CCRHS_H(DSRHF,OMEGA1,XLAMDP,XLAMDH,SCRM,
     *                   WORK,LWORK,ISYDIS,ISYDEL,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to do linear transformation parts by
C     OC 30-1-1995
C
C     Purpose: Calculate H-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DSRHF(*),OMEGA1(*),XLAMDH(*),WORK(LWORK)
      DIMENSION XLAMDP(*),SCRM(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_H1(DSRHF,OMEGA1,SCRM,XLAMDP,XLAMDH,WORK,LWORK,
     *              ISYDIS,ISYDEL,ISYMTR)
C
C
      RETURN
      END
      SUBROUTINE CCRHS_H1(DSRHF,OMEGA1,SCRM,XLAMDP,XLAMDH,WORK,LWORK,
     *                    ISYDIS,ISYDEL,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to do linear transformation parts by
C     OC 30-1-1995
C
C     Purpose: Calculate H-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      DIMENSION DSRHF(*),OMEGA1(*),SCRM(*)
      DIMENSION XLAMDP(*),XLAMDH(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C--------------------------------------
C     Calculate contribution.
C--------------------------------------
C
      ISYAKL = MULD2H(ISYMTR,ISYDEL)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMGB = MULD2H(ISYML,ISYDIS)
         ISYMAK = MULD2H(ISYML,ISYAKL)
         ISYMKI = ISYMGB
C
         KSCR1 = 1
         KEND1 = KSCR1 + N2BST(ISYMGB)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
            CALL QUIT('Insufficient space in CCRHS_H1')
         ENDIF
C
         DO 110 L = 1,NRHF(ISYML)
C
            KOFF1 = IDSRHF(ISYMGB,ISYML) + NNBST(ISYMGB)*(L - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMGB,WORK(KSCR1))
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMB = ISYMI
               ISYMG = MULD2H(ISYMB,ISYMGB)
               ISYMK = ISYMG
               ISYMA = MULD2H(ISYMK,ISYMAK)
C
               KSCR2 = KEND1
               KSCR3 = KSCR2 + NBAS(ISYMG)*NRHF(ISYMI)
               KEND2 = KSCR3 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2 = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHS_H1')
               ENDIF
C
               NBASG = MAX(NBAS(ISYMG),1)
               NBASB = MAX(NBAS(ISYMB),1)
               NRHFK = MAX(NRHF(ISYMK),1)
               NVIRA = MAX(NVIR(ISYMA),1)
C
               KOFF2 = KSCR1 + IAODIS(ISYMG,ISYMB)
               KOFF3 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NBAS(ISYMB),
     *                    ONE,WORK(KOFF2),NBASG,XLAMDH(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR2),NBASG)
C
               KOFF4 = ILMRHF(ISYMK) + 1
C
               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,XLAMDP(KOFF4),NBASG,WORK(KSCR2),NBASG,
     *                    ZERO,WORK(KSCR3),NRHFK)
C
               KOFF5 = IT2BCD(ISYMAK,ISYML) + NT1AM(ISYMAK)*(L - 1)
     *               + IT1AM(ISYMA,ISYMK) + 1
               KOFF6 = IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *                    -ONE,SCRM(KOFF5),NVIRA,WORK(KSCR3),NRHFK,
     *                    ONE,OMEGA1(KOFF6),NVIRA)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_g */
      SUBROUTINE CCRHS_G(DSRHF,OMEGA1,XLAMP1,ISYMP1,XLAMH1,ISYMH1,SCRM,
     *                   WORK,LWORK,ISYDIS,ISYDEL,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to calculated term of linear transformation
C     and handle different transformations on integral indices by OC 30-1-1995
C
C     G(a,i) = sum(cdk)[t(ci,dk)*Lackd]
C     G(a,i)for fixed del = sum(ck)[M(ci,k)*L(alfa gamma k]
C
C     XLAMP1 is the transformation matrix for a ; XLAMP or a oneindex transformed.
C     XLAMH1 is the transformation matrix for c ; XLAMH or a oneindex transformed.
C     DSRHF is the (alfa gamma | k) array for a given delta.
C
C     not implemented yet with DSRHF and SCRM index transformed.
C
C     tested for energy with symmetry: ordinary XLAM matrices  - OC 10-2-1995
C     tested for linear transformation without symmetry.       - OC spring 95.
C
C     Purpose: Calculate G-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DSRHF(*),OMEGA1(*),XLAMP1(*),WORK(LWORK)
      DIMENSION XLAMH1(*),SCRM(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      ISYINT = MULD2H(ISYMH1,ISYMOP)
      ISYALI = MULD2H(ISYINT,ISYMTR)
C
      KSCR1  = 1
      KEND1  = KSCR1  + NT1AO(ISYALI)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_G')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_G1(DSRHF,OMEGA1,SCRM,XLAMP1,ISYMP1,XLAMH1,ISYMH1,
     *              WORK(KSCR1),WORK(KEND1),LWRK1,ISYDIS,ISYDEL,ISYMTR)
C
C
      RETURN
      END
      SUBROUTINE CCRHS_G1(DSRHF,OMEGA1,SCRM,XLAMP1,ISYMP1,XLAMH1,ISYMH1,
     *                    SCR1,WORK,LWORK,ISYDIS,ISYDEL,ISYMTR)
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch & Ove Christiansen 19-Jan-1994
C     Generalized to calculated term of linear transformation
C     by OC 30-1-1995
C
C     Purpose: Calculate G-term.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      PARAMETER(TWO=2.0D0)
      DIMENSION DSRHF(*),OMEGA1(*),SCRM(*),SCR1(*)
      DIMENSION XLAMP1(*),XLAMH1(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYINT = MULD2H(ISYMH1,ISYMOP)
      ISYALI = MULD2H(ISYINT,ISYMTR)
      ISYMAI = MULD2H(ISYALI,ISYMP1)
      ISYCIK = MULD2H(ISYMTR,ISYDEL)
C
      CALL DZERO(SCR1,NT1AO(ISYMAI))
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMAG = MULD2H(ISYMK,ISYDIS)
         ISYMCI = MULD2H(ISYMK,ISYCIK)
         ISYMGI = MULD2H(ISYALI,ISYMAG)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMAG)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
            CALL QUIT('Insufficient space in CCRHS_G1')
         ENDIF
C
         DO 110 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYMAG,ISYMK) + NNBST(ISYMAG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR10))
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMG = MULD2H(ISYMI,ISYMGI)
               ISYMA = MULD2H(ISYMG,ISYMAG)
               ISYMC = ISYMG
C
               NBASG = MAX(NBAS(ISYMG),1)
               NBASA = MAX(NBAS(ISYMA),1)
               NVIRC = MAX(NVIR(ISYMC),1)
C
               KSCR11 = KEND1
               KEND2  = KSCR11 + NBAS(ISYMG)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCRHS_G1')
               ENDIF
C
               KOFF2 = IGLMVI(ISYMG,ISYMC) + 1
               KOFF3 = IT2BCD(ISYMCI,ISYMK) + NT1AM(ISYMCI)*(K - 1)
     *               + IT1AM(ISYMC,ISYMI) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC),
     *                    ONE,XLAMH1(KOFF2),NBASG,SCRM(KOFF3),NVIRC,
     *                    ZERO,WORK(KSCR11),NBASG)
C
               KOFF4 = KSCR10 + IAODIS(ISYMA,ISYMG)
               KOFF6 = IT1AO(ISYMA,ISYMI) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF4),NBASA,WORK(KSCR11),NBASG,
     *                    ONE,SCR1(KOFF6),NBASA)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C----------------------------------------------
C     Accumulation into OMEGA1 in the MO basis.
C----------------------------------------------
C
      DO 200 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMAI)
         ISYMAL= MULD2H(ISYMI,ISYALI)
C
         NBASA = MAX(NBAS(ISYMA),1)
         NVIRA = MAX(NVIR(ISYMA),1)
C
         KOFF1 = IGLMVI(ISYMAL,ISYMA) + 1
         KOFF2 = IT1AO(ISYMA,ISYMI) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMA),ONE,
     *              XLAMP1(KOFF1),NBASA,SCR1(KOFF2),NBASA,ONE,
     *              OMEGA1(KOFF3),NVIRA)
C
  200 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_ei */
      SUBROUTINE CCRHS_EI(DSRHF,EMAT1,EMAT2,T2AM,SCRM,XLAMDP,XLAMDH,
     *                   WORK,LWORK,IDEL,ISYMD,ISYDIS,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch 12-Jan-1994
C     Symmetry 2-aug
C     Modified slightly by Ove Christiansen 31-1-95 for construction of
C     linear transformation intermediates. ISYMTR = SYM OF T2-VEC
C
C     Purpose: Calculate E-intermediates.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION EMAT1(*), EMAT2(*)
      DIMENSION DSRHF(*),WORK(LWORK)
      DIMENSION T2AM(*),SCRM(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KSCR1  = 1
      KSCR2  = KSCR1  + NT2BCD(ISYDIS)
      KSCR3  = KSCR2  + NT2BGD(ISYDIS)
      KEND1  = KSCR3  + NT2BGD(ISYDIS)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_EI')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_EI1(DSRHF,EMAT1,EMAT2,T2AM,SCRM,
     *              WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),
     *              XLAMDP,XLAMDH,WORK(KEND1),LWRK1,IDEL,
     *              ISYMD,ISYDIS,ISYMTR)
C
      RETURN
      END
      SUBROUTINE CCRHS_EI1(DSRHF,EMAT1,EMAT2,T2AM,SCRM,SCR1,SCR2,
     *                    SCR3,XLAMDP,XLAMDH,WORK,LWORK,IDEL,
     *                    ISYMD,ISYDIS,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch 12-Jan-1994
C     Symmetry 2-aug
C     Modified slightly by Ove Christiansen 31-1-95 for construction of
C     linear transformation intermediates. ISYMTR = SYM OF T2-VEC
C
C     Purpose: Calculate E-intermediates.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      DIMENSION DSRHF(*)
      DIMENSION EMAT1(*),EMAT2(*)
      DIMENSION T2AM(*),SCRM(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C
C===================================
C     First intermediate I(b,delta).
C===================================
C
C-------------------------------------------------------
C     Construct the integrals I(dl,m) = (l d | m delta).
C-------------------------------------------------------
C
      DO 100 ISYMM = 1,NSYM
C
         ISYMAG = MULD2H(ISYMM,ISYDIS)
         ISYMDL = ISYMAG
         ISYMGL = ISYMAG
C
         DO 110 M = 1,NRHF(ISYMM)
C
            KSCR1 = 1
            KEND1 = KSCR1 + N2BST(ISYMAG)
            LWRK1 = LWORK - KEND1
            IF (LWRK1. LT. 0) THEN
               CALL QUIT('Insufficient core in CCRHS_EI1')
            END IF
C
            KOFF1 = IDSRHF(ISYMAG,ISYMM) + NNBST(ISYMAG)*(M - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR1))
C
            DO 120 ISYML = 1,NSYM
C
               ISYMD1 = MULD2H(ISYML,ISYMDL)
               ISYMA  = ISYML
               ISYMG  = ISYMD1
C
               NBASA = MAX(NBAS(ISYMA),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NVIRD = MAX(NVIR(ISYMD1),1)
C
               KOFF2 = KSCR1 + IAODIS(ISYMA,ISYMG)
               KOFF3 = ILMRHF(ISYML) + 1
               KOFF4 = IT2BGD(ISYMGL,ISYMM) + NT1AO(ISYMGL)*(M - 1)
     *               + IT1AO(ISYMG,ISYML) + 1
C
               CALL DGEMM('T','N',NBAS(ISYMG),NRHF(ISYML),NBAS(ISYMA),
     *                    ONE,WORK(KOFF2),NBASA,XLAMDP(KOFF3),NBASA,
     *                    ZERO,SCR2(KOFF4),NBASG)
C
               KOFF5 = ILMVIR(ISYMD1) + 1
               KOFF6 = IT2BCD(ISYMDL,ISYMM) + NT1AM(ISYMDL)*(M - 1)
     *               + IT1AM(ISYMD1,ISYML) + 1
C
               CALL DGEMM('T','N',NVIR(ISYMD1),NRHF(ISYML),
     *                    NBAS(ISYMG),ONE,XLAMDH(KOFF5),NBASG,
     *                    SCR2(KOFF4),NBASG,ZERO,SCR1(KOFF6),NVIRD)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C-------------------------------------------------------
C     Contract the integrals in SCR1 with t2 amplitudes.
C-------------------------------------------------------
C
      DO 200 ISYMM = 1,NSYM
C
         ISYMDL = MULD2H(ISYMM,ISYDIS)
         ISYMBM = MULD2H(ISYMDL,ISYMTR)
         ISYMB  = MULD2H(ISYMBM,ISYMM)
C
         DO 210 M = 1,NRHF(ISYMM)
C
            NT1DL = MAX(NT1AM(ISYMDL),1)
C
            KBM   = IT1AM(ISYMB,ISYMM) + NVIR(ISYMB)*(M - 1) + 1
            KOFF1 = IT2SQ(ISYMDL,ISYMBM)
     *            + NT1AM(ISYMDL)*(KBM - 1) + 1
            KOFF2 = IT2BCD(ISYMDL,ISYMM)
     *            + NT1AM(ISYMDL)*(M - 1) + 1
            KOFF3 = IEMAT1(ISYMB,ISYMD)
     *            + (IDEL - IBAS(ISYMD) - 1)*NVIR(ISYMB) + 1
C
            CALL DGEMV('T',NT1AM(ISYMDL),NVIR(ISYMB),ONE,T2AM(KOFF1),
     *                 NT1DL,SCR1(KOFF2),1,ONE,EMAT1(KOFF3),1)
C
  210    CONTINUE
C
  200 CONTINUE
C
C================================
C     Second intermediate I(k,j).
C================================
C
C-------------------------------------------
C     Transform the SCRM amplitudes to SCR3.
C-------------------------------------------
C
      DO 300 ISYMJ = 1,NSYM
C
         ISYMDJ = MULD2H(ISYMD,ISYMJ)
         ISYMEM = MULD2H(ISYMDJ,ISYMTR)
         ISYMGM = ISYMEM
C
         DO 310 J = 1,NRHF(ISYMJ)
C
            DO 320 ISYMM = 1,NSYM
C
               ISYME = MULD2H(ISYMM,ISYMEM)
               ISYMG = ISYME
C
               NVIRE = MAX(NVIR(ISYME),1)
               NBASG = MAX(NBAS(ISYMG),1)
C
               KOFF1 = ILMVIR(ISYME) + 1
               KOFF2 = IT2BCD(ISYMEM,ISYMJ) + NT1AM(ISYMEM)*(J - 1)
     *               + IT1AM(ISYME,ISYMM) + 1
               KOFF3 = IT2BGD(ISYMGM,ISYMJ) + NT1AO(ISYMGM)*(J - 1)
     *               + IT1AO(ISYMG,ISYMM) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMM),NVIR(ISYME),
     *                    ONE,XLAMDH(KOFF1),NBASG,SCRM(KOFF2),NVIRE,
     *                    ZERO,SCR3(KOFF3),NBASG)
C
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
C
C----------------------------------------------------------------
C     Contract the integrals in SCR2 with the amplitudes in SCR3.
C----------------------------------------------------------------
C
      DO 400 ISYMJ = 1,NSYM
C
         ISYMDJ = MULD2H(ISYMD,ISYMJ)
         ISYMEM = MULD2H(ISYMDJ,ISYMTR)
         ISYMGM = ISYMEM
         ISYMK  = MULD2H(ISYMGM,ISYDIS)
C
         NT1GM = MAX(NT1AO(ISYMGM),1)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         KOFF1 = IT2BGD(ISYMGM,ISYMK) + 1
         KOFF2 = IT2BGD(ISYMGM,ISYMJ) + 1
         KOFF3 = IMATIJ(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),NT1AO(ISYMGM),
     *              ONE,SCR2(KOFF1),NT1GM,SCR3(KOFF2),NT1GM,
     *              ONE,EMAT2(KOFF3),NRHFK)
C
  400 CONTINUE
C
      RETURN
      END
C  /* Deck cc_aofock */
      SUBROUTINE CC_AOFOCK(XINT,DENSIT,FOCK,WORK,LWORK,IDEL,
     *                      ISYMD,LAUXD,IBASX,ISYDEN)
C
C     Written by Asger Halkier and Henrik Koch 27-4-95.
C
C     Debugged Ove Christiansen august 1995
C
C     Purpose: Calculate the two electron contribution to the
C              AO-fock matrix using matrix vector routines.
C
C     Obs: It can be done as F(g>=d) = G(a>=b) I(a>=b,g,d) where
C          G(a>=b) = D(a,b) + D(b,a), the diagonal properly scaled
C
C     Adapted for R12: LAUXD=.TRUE.: Delta runs also over aux-functions
C     Christian Neiss, spring 2006 
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "r12int.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),DENSIT(*)
      DIMENSION FOCK(*),WORK(LWORK)
      LOGICAL   LAUXD
      INTEGER   IBASX(8),NBAS2(8),NGDP(8),IGDP(8,8)
C
      IF (LAUXD) THEN
        DO ISYM = 1, NSYM
          NBAS2(ISYM) = MBAS1(ISYM)+MBAS2(ISYM)
        END DO
        DO ISYM = 1, NSYM
          NGDP(ISYM) = 0
          DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            IGDP(ISYM1,ISYM2) = NGDP(ISYM)
            NGDP(ISYM) = NGDP(ISYM) + MBAS1(ISYM1)*NBAS2(ISYM2)
          END DO
        END DO
      END IF
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      DO 100 ISYMG = 1,NSYM
C
         IF (NBAS(ISYMG) .EQ. 0) GOTO 100
C
         ISYMAB = MULD2H(ISYMG,ISYDIS)
C
         NDISTG = MIN(LWORK/MAX(N2BST(ISYMAB),1),NBAS(ISYMG))
C
         IF (NDISTG .LT. 1) THEN
            CALL QUIT('Insufficient work space in CC_AOFOCK1')
         ENDIF
C
         NBATCH = (NBAS(ISYMG) - 1)/NDISTG + 1
C
C-------------------------------------
C        Start the loops over batches.
C-------------------------------------
C
         DO 110 IBATCH = 1,NBATCH
C
            NUMG = NDISTG
            IF (IBATCH .EQ. NBATCH) THEN
               NUMG = NBAS(ISYMG) - NDISTG*(NBATCH - 1)
            ENDIF
C
            IG1 = NDISTG*(IBATCH - 1) + 1
            IG2 = NDISTG*(IBATCH - 1) + NUMG
C
            KOFF2 = 1
            DO 120 IG = IG1,IG2
C
               KOFF1 = IDSAOG(ISYMG,ISYDIS)
     *               + (IG - 1)*NNBST(ISYMAB) + 1
C
               CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,WORK(KOFF2))
C
               KOFF2 = KOFF2 + N2BST(ISYMAB)
C
  120       CONTINUE
C
            IF (ISYMAB .EQ. ISYDEN) THEN
C
               IF (LAUXD) THEN
                 KGD = IGDP(ISYMG,ISYMD)
     *               + (IDEL-IBAS(ISYMD)-IBASX(ISYMD) - 1)*NBAS(ISYMG) 
     *               + IG1
               ELSE
                 KGD = IAODIS(ISYMG,ISYMD)
     *               + (IDEL-IBAS(ISYMD) - 1)*NBAS(ISYMG) + IG1
               END IF
C
               NTOBST = MAX(N2BST(ISYMAB),1)
C
               CALL DGEMV('T',N2BST(ISYMAB),NUMG,TWO,WORK,NTOBST,
     *                    DENSIT,1,ONE,FOCK(KGD),1)
C
            ENDIF
C
            ISYMA = MULD2H(ISYMD,ISYDEN)
            ISYMB = MULD2H(ISYMA,ISYMAB)
C
            IF (LAUXD) THEN
              KAD = IGDP(ISYMA,ISYMD)
     *            + NBAS(ISYMA)*(IDEL-IBAS(ISYMD)-IBASX(ISYMD) - 1) + 1
            ELSE
              KAD = IAODIS(ISYMA,ISYMD)
     *            + NBAS(ISYMA)*(IDEL - IBAS(ISYMD) - 1) + 1
            END IF
C
            DO 130 IG = IG1,IG2
C
               KOFF1 = (IG - IG1)*N2BST(ISYMAB)
     *               + IAODIS(ISYMA,ISYMB) + 1
               KGB   = IAODIS(ISYMG,ISYMB) + IG
C
               NTOTA = MAX(NBAS(ISYMA),1)
               NTOTG = MAX(NBAS(ISYMG),1)
C
C              CALL DGEMV('T',NBAS(ISYMA),NBAS(ISYMB),-ONE,WORK(KOFF1),
C    *                    NTOTA,DENSIT(KAD),1,ONE,FOCK(KGB),NTOTG)
C
               CALL DGEMV('N',NBAS(ISYMA),NBAS(ISYMB),-ONE,WORK(KOFF1),
     *                    NTOTA,DENSIT(KGB),NTOTG,ONE,FOCK(KAD),1)
C
  130       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_d */
      SUBROUTINE CCRHS_D(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                   XLAMDP,XLAMIP,XLAMDH,
     *                   XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                   SCRM,E1PIM,WORK,LWORK,IDEL,ISYMD,FACTD,ICON,
     *                   IOPTR12,IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
C
C     Written by Henrik Koch 9-Jan-1994
C
C     Generalisation for CCLR by Ove Christiansen august-september 1995
C     (right transformation) and september 1996 (F-matrix).
C
C     adapted for CCSDR12, C. Neiss, spring 2006
C     IOPTR12 = 1 Calculate both conv. D and r12 D' intermediates
C                 T2-dependent contr. to D' interm. is added with a prefactor
C                 of 2*FACTD
C     IOPTE   = 1 Calculate the T-dependent part of the
C                 E_{a delta')^1' intermediate (on E1PIM).
C 
C     Purpose: Calculate D-term.
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),DSRHF(*),OMEGA2(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*),SCRM(*)
      DIMENSION XLAMPC(*),XLAMHC(*)
      DIMENSION T2AM(*),E1PIM(*)
      CHARACTER DFIL*(*),DPFIL
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      IF (ISYMT2 .NE. ISYMPC) CALL QUIT('Symmetry Mismatch in CCRHS_D' )
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KSCR1  = 1
      KSCR2  = KSCR1  + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3  = KSCR2  + NT2BGD(ISYDIS)
      IF (ICON .EQ. 2) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF
      IF (IOPTR12.EQ.1) THEN
         KSCR4  = KEND1
         KEND1  = KSCR4  + NT2BCD(ISYAIK)
      END IF

      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_D')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_D1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *              SCRM,E1PIM,WORK(KSCR1),
     *              WORK(KSCR2),WORK(KSCR3),WORK(KSCR4),XLAMDP,XLAMIP,
     *              XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *              WORK(KEND1),LWRK1,ISYDIS,IDEL,
     *              ISYMD,FACTD,ICON,IOPTR12,IOPTE,
     *              LUD,DFIL,LUDP,DPFIL,IV)
C
      RETURN
      END
      SUBROUTINE CCRHS_D1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                    SCRM,E1PIM,SCR1,SCR2,SCR3,SCR4,
     *                    XLAMDP,XLAMIP,XLAMDH,XLAMPC,ISYMPC,XLAMHC,
     *                    ISYMHC,WORK,LWORK,ISYDIS,IDEL,ISYDEL,FACTD,
     *                    ICON,IOPTR12,IOPTE,LUD,DFIL,LUDP,DPFIL,IV)
C
C     Written by Henrik Koch 3-Jan-1994
C
C     Modification by Ove Christiansen 25-7-1995 to allow for a
C     general factor (FACTD). NB: Assumes DUMCD. 
C     - calculate intermediates for CCLR.
C
C     29-9-1995 (17-9-1996 F-matrix.) Ove Christiansen:  
C     
C                 1. If Icon = 2 both contributions are calculated,
C                    for total sym. case. Otherwise 
C                    a.ICON = 1 only the integral Laikc(del)
C                               = La-bar,i,k,c + La,i-bar,k,c
C                      for Jacobian right transformation
C                    b.ICON = 3 
C                          La-bar,i,k,c + La,i-bar,k,c + Tx*Int
C                      for F-matrix times vector.
C                              
C                 2. Allow for general transformation matrix for
C                    alpha to a(XLAMPC) and for i (XLAMHC).
C                    (the extra i transformation introduces new
C                     blocks which is only entered when icon = 1 or 3)
C
C                 3. If icon diff. from 2 (we have linear response)
C                    The D intermediate is stored according to
C                    the number of simultaneous trial vector
C                    given by IV. This is ensured using IT2DLR.
C
C
C     This to calculate terms in 2C1 right transformation in CCLR.
C
C     adapted for CCSDR12, C. Neiss spring 2006
C
C     Purpose: Calculate D-term.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccsdinp.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,XMHALF=-0.5D0)
      PARAMETER(TWO=2.0D0)
      DIMENSION XINT(*),OMEGA2(*),T2AM(*),DSRHF(*)
      DIMENSION SCRM(*),E1PIM(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*)
      DIMENSION XLAMPC(*),XLAMHC(*)
      DIMENSION WORK(LWORK)
      INTEGER   NADP(8),IADP(8,8),IBASX(8)
      CHARACTER DFIL*(*),DPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "r12int.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (IOPTE.EQ.1) THEN
        IF (.NOT.CCR12) CALL QUIT('IOPTE only implemented for CC-R12')
        IBASX(1) = 0
        DO ISYM = 2, NSYM
          IBASX(ISYM) = IBASX(ISYM-1) + MBAS2(ISYM-1)
        END DO
        DO ISYM = 1, NSYM
          NADP(ISYM) = 0
          DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            IADP(ISYM1,ISYM2) = NADP(ISYM)
            NADP(ISYM) = NADP(ISYM) +
     &                   NVIR(ISYM1)*(MBAS1(ISYM2)+MBAS2(ISYM2))
          END DO
        END DO
      END IF
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C-------------------------------------------------------
C     Calculate the integrals K(k,dl) = (k d | l delta).
C-------------------------------------------------------
C
      IF (ICON .GE. 2) THEN
C
         DO 100 ISYMK = 1,NSYM
C
            ISYMAG = MULD2H(ISYMK,ISYDIS)
C
            DO 110 K = 1,NRHF(ISYMK)
C
               ISYMDL = MULD2H(ISYMK,ISYDIS)
C
               KSCR10 = 1
               KEND1  = KSCR10 + N2BST(ISYMAG)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_D1')
               END IF
C
               KOFF1 = IDSRHF(ISYMAG,ISYMK) + NNBST(ISYMAG)*(K-1) + 1
               CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR10))
C
               DO 120 ISYML = 1,NSYM
C
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  ISYMA = ISYML
                  ISYMG = ISYMD
C
                  NBASA = MAX(NBAS(ISYMA),1)
                  NBASG = MAX(NBAS(ISYMG),1)
                  NVIRD = MAX(NVIR(ISYMD),1)
C
                  KSCR11 = KEND1
                  KEND2  = KSCR11 + NBAS(ISYMG)*NRHF(ISYML)
                  LWRK2  = LWORK  - KEND2
C
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D1')
                  END IF
C
                  KOFF2 = KSCR10 + IAODIS(ISYMA,ISYMG)
                  KOFF3 = ILMRHF(ISYML) + 1
C
                  CALL DGEMM('T','N',NBAS(ISYMG),NRHF(ISYML),
     *                       NBAS(ISYMA),ONE,WORK(KOFF2),NBASA,
     *                       XLAMDP(KOFF3),NBASA,
     *                       ZERO,WORK(KSCR11),NBASG)
C
                  KOFF5 = ILMVIR(ISYMD) + 1
                  KOFF6 = IT2BCD(ISYMDL,ISYMK) + NT1AM(ISYMDL)*(K - 1)
     *                  + IT1AM(ISYMD,ISYML) + 1
C
                  CALL DGEMM('T','N',NVIR(ISYMD),NRHF(ISYML),
     *                       NBAS(ISYMG),ONE,XLAMDH(KOFF5),NBASG,
     *                       WORK(KSCR11),NBASG,
     *                       ZERO,SCR1(KOFF6),NVIRD)
C
  120          CONTINUE
C
  110       CONTINUE
C
  100    CONTINUE
C
C---------------------------------
C        Transpose integral array.
C---------------------------------
C
         CALL CC_MTCME(SCR1,WORK,LWORK,ISYDIS,1)
C
         IF (LWORK .LT. NT2BCD(ISYDIS)) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_D1')
         END IF
C
         DO 130 ISYMK = 1,NSYM
C
            ISYMDL = MULD2H(ISYMK,ISYDIS)
C
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 140 K = 1,NRHF(ISYMK)
C
               KOFF1 = IT2BCD(ISYMDL,ISYMK) + NT1AM(ISYMDL)*(K - 1) + 1
               KOFF2 = IT2BCT(ISYMK,ISYMDL) + K
C
               CALL DCOPY(NT1AM(ISYMDL),SCR1(KOFF1),1,WORK(KOFF2),NRHFK)
C
  140       CONTINUE
C
  130    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYDIS),WORK,1,SCR1,1)
C
C-----------------------------------------
C        Calculate the first contribution.
C        sum(2*t(ai,dl)-t(di,al))*L(ldkc)
C-----------------------------------------
C
         IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
            CALL QUIT('Insufficient work space in CCRHS_D1')
         ENDIF
C
         DO 200 ISYMK = 1,NSYM
C
            ISYMDL = MULD2H(ISYMK,ISYDIS)
            ISYMAI = MULD2H(ISYAIK,ISYMK)
C
            NRHFK  = MAX(NRHF(ISYMK),1)
            NTOTDL = MAX(NT1AM(ISYMDL),1)
C
            KOFF1  = IT2BCT(ISYMK,ISYMDL) + 1
            KOFF2  = IT2SQ(ISYMDL,ISYMAI) + 1
            KOFF3  = IT2BCT(ISYMK,ISYMAI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                 ONE,SCR1(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                 WORK(KOFF3),NRHFK)
C
  200    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR1,1)
C
         !save a copy of first contribution:
         IF (IOPTR12.EQ.1) THEN
           CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR4,1)
         END IF
C
      ENDIF
C
C----------------------------------------------------------
C     Calculate the integrals K(k,ai) = (k i | alfa delta).
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_D1')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_D1')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KSCR11),NRHFK,XLAMDH(KOFF5),NBASG,
     *                    ZERO,WORK(KSCR12),NRHFK)
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C
  330          CONTINUE
C
C-------------------------------------------------------
C              In 2C1 linear transformation extra  cont.
C-------------------------------------------------------
C
               IF ((ICON .EQ. 1) .OR. (ICON.EQ.3)) THEN
C
                  ISYMI  = MULD2H(ISYMG,ISYMHC)
                  ISYMAI = MULD2H(ISYMA,ISYMI)
C
                  KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D1')
                  END IF
C
                  KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMHC(KOFF5),NBASG,
     *                       ZERO,WORK(KSCR12),NRHFK)
C
                  DO 331 I = 1,NRHF(ISYMI)
C
                     NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                     KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                     + NRHF(ISYMK)*(NAI - 1) + 1
                     KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR3(KOFF8),1)
C
  331             CONTINUE
C
               ENDIF
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
C
      CALL DSCAL(NT2BGD(ISYDIS),-ONE,SCR2,1)
C
      ISALIK = MULD2H(ISYDIS,ISYMHC)
C
      CALL DSCAL(NT2BGD(ISALIK),-ONE,SCR3,1)
C
      DO 340 ISYMK = 1,NSYM
C
         ISYALG = MULD2H(ISYMK,ISYDIS)
         ISYALI = MULD2H(ISYMHC,ISYALG)
         NT1AOM = MAX(NT1AO(ISYALG),NT1AO(ISYALI))
C
         KSCR10 = 1
         KSCR11 = KSCR10 + N2BST(ISYALG)
         KEND1  = KSCR11 + NT1AOM
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_D1')
         END IF
C
         DO 350 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYALG,ISYMK) + NNBST(ISYALG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYALG,WORK(KSCR10))
C
            ISYALI = ISYALG
            CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
C------------------------------
C           Usual contribution.
C------------------------------
C
            DO 360 ISYMI = 1,NSYM
C
               ISYMAL = MULD2H(ISYMI,ISYALI)
               ISYMG  = ISYMI
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASG = MAX(NBAS(ISYMG),1)
C
               KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
               KOFF3 = ILMRHF(ISYMI) + 1
               KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF2),NBASAL,XLAMDH(KOFF3),NBASG,
     *                    ZERO,WORK(KOFF4),NBASAL)
C
  360       CONTINUE
C
            NRHFK = MAX(NRHF(ISYMK),1)
            KOFF5 = IT2BGT(ISYMK,ISYALI) + K
            CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,SCR2(KOFF5),
     *                 NRHFK)
C
C----------------------------------------------------
C           In 2C1 linear tronsformation extra  cont.
C----------------------------------------------------
C
            IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
               ISYALI = MULD2H(ISYALG,ISYMHC)
C
               CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
               DO 361 ISYMI = 1,NSYM
C
                  ISYMAL = MULD2H(ISYMI,ISYALI)
                  ISYMG  = MULD2H(ISYMI,ISYMHC)
C
                  NBASAL = MAX(NBAS(ISYMAL),1)
                  NBASG  = MAX(NBAS(ISYMG),1)
C
                  KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
                  KOFF3 = IGLMRH(ISYMG,ISYMI) + 1
                  KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KOFF2),NBASAL,
     *                       XLAMHC(KOFF3),NBASG,
     *                       ZERO,WORK(KOFF4),NBASAL)
C
  361          CONTINUE
C
               NRHFK = MAX(NRHF(ISYMK),1)
               KOFF5 = IT2BGT(ISYMK,ISYALI) + K
               CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,
     *                    SCR3(KOFF5),NRHFK)
C
            ENDIF
C
  350    CONTINUE
C
  340 CONTINUE
C
      IF (DUMPCD) GOTO 700
C
      IF (CCR12) CALL QUIT('CCSDR12 requires DUMPCD=.TRUE.')
C
C-----------------------------------------
C     Back transformation to the AO basis.
C-----------------------------------------
C
      DO 400 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYDIS)
C
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 410 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            NBASA = MAX(NBAS(ISYMA),1)
C
            DO 420 I = 1,NRHF(ISYMI)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I - 1) + 1
C
               KOFF1 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
               KOFF2 = ILMVIR(ISYMA) + 1
               KOFF3 = IT2BGT(ISYMK,ISYMAI) + NRHF(ISYMK)*(MAI - 1) + 1
C
               CALL DGEMM('N','T',NRHF(ISYMK),NBAS(ISYMA),NVIR(ISYMA),
     *                    HALF,SCR1(KOFF1),NRHFK,XLAMIP(KOFF2),NBASA,
     *                    ONE,SCR2(KOFF3),NRHFK)
C
  420       CONTINUE
C
  410    CONTINUE
C
  400 CONTINUE
C
C
      DO 500 ISYMK = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         DO 510 K = 1,NRHF(ISYMK)
C
            DO 520 ISYMJ = 1,NSYM
C
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NVIRB = MAX(NVIR(ISYMB),1)
C
               KOFF1 = ILMVIR(ISYMB) + 1
               KOFF2 = IT2BCD(ISYMBJ,ISYMK) + NT1AM(ISYMBJ)*(K - 1)
     *               + IT1AM(ISYMB,ISYMJ) + 1
               KOFF3 = IT2BGD(ISYMBJ,ISYMK) + NT1AO(ISYMBJ)*(K - 1)
     *               + IT1AO(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMB),NRHF(ISYMJ),NVIR(ISYMB),
     *                    ONE,XLAMIP(KOFF1),NBASB,SCRM(KOFF2),NVIRB,
     *                    ZERO,SCR3(KOFF3),NBASB)
C
  520       CONTINUE
C
  510    CONTINUE
C
  500 CONTINUE
C
C---------------------------------------
C     Calculate the second contribution.
C---------------------------------------
C
      DO 600 ISYMAI = 1,NSYM
C
         ISYMK  = MULD2H(ISYMAI,ISYDIS)
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         IF (NRHF(ISYMK) .EQ. 0) GOTO 600
C
         IF (LWORK .LT. NT1AO(ISYMBJ)) THEN
            CALL QUIT('Insufficient work space in CCRHS_D1')
         ENDIF
C
         NTOTBJ = MAX(NT1AO(ISYMBJ),1)
         NRHFK  = MAX(NRHF(ISYMK),1)
C
         IF (.NOT. OMEGSQ) THEN
C
            DO 610 NAI = 1,NT1AO(ISYMAI)
C
               KOFF1 = IT2BGD(ISYMBJ,ISYMK) + 1
               KOFF2 = IT2BGT(ISYMK,ISYMAI)
     *               + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMV('N',NT1AO(ISYMBJ),NRHF(ISYMK),ONE,
     *                    SCR3(KOFF1),NTOTBJ,SCR2(KOFF2),1,
     *                    ZERO,WORK,1)
C
               IF (ISYMAI .EQ. ISYMBJ) THEN
                  WORK(NAI) = TWO*WORK(NAI)
               ENDIF
C
               DO 620 NBJ = 1,NT1AO(ISYMBJ)
                  NAIBJ = IT2AO(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*WORK(NBJ)
  620          CONTINUE
C
  610       CONTINUE
C
         ELSE
C
            KOFF1 = IT2BGD(ISYMBJ,ISYMK)  + 1
            KOFF2 = IT2BGT(ISYMK,ISYMAI)  + 1
            KOFF3 = IT2AOS(ISYMBJ,ISYMAI) + 1
C
            CALL DGEMM('N','N',NT1AO(ISYMBJ),NT1AO(ISYMAI),NRHF(ISYMK),
     *                 HALF,SCR3(KOFF1),NTOTBJ,SCR2(KOFF2),NRHFK,
     *                 ONE,OMEGA2(KOFF3),NT1AO(ISYMBJ))
C
         ENDIF
C
  600 CONTINUE
C
      GOTO 999
C
C-------------------
C     I/O algorithm.
C-------------------
C
  700 CONTINUE
C
C--------------------------------------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     for 2C1 transformation this means lamdpc is a C1 transformed lambda.
C--------------------------------------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      DO 710 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 720 ISYMI = 1,NSYM
C
            ISYMA  = MULD2H(ISYMI,ISYMAI)
            ISYMAL = MULD2H(ISYMPC,ISYMA)
            ISYALI = MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 730 I = 1,NRHF(ISYMI)
C
               NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1) + 1
               MALI  = IT1AO(ISYMAL,ISYMI)  + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI - 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTD ,SCR1(KOFF3),NRHFK)
C
               IF (IOPTE.EQ.1) THEN
                 IF (ISYMI.EQ.ISYMK) THEN
                   KOFF3 = IT2BCT(ISYMK,ISYMAI) +
     &                     NRHF(ISYMK)*(NAI - 1) + I
                   IF (IDEL.GT.NBAST) THEN
                     D = IDEL-IBASX(ISYDEL)-NBAST+MBAS1(ISYDEL)
                   ELSE
                     D = IDEL-IBAS(ISYDEL)
                   END IF
                   KOFFE = IADP(ISYMA,ISYDEL) +
     &                     NVIR(ISYMA)*(D-1) + 1
                   CALL DAXPY(NVIR(ISYMA),-0.5D0,SCR1(KOFF3),
     &                        NRHF(ISYMK),E1PIM(KOFFE),1)
                 END IF
               END IF
C
  730       CONTINUE
  720    CONTINUE
  710 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     I is C1 transformed.
C-----------------------------------------------
C
      IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
         ISYAIK = MULD2H(ISYDIS,ISYMHC)
C
         DO 750 ISYMAI = 1,NSYM
C
            ISYMK = MULD2H(ISYMAI,ISYAIK)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 760 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAL= ISYMA
               ISYALI= MULD2H(ISYMAL,ISYMI)
               NBASAL = MAX(NBAS(ISYMAL),1)
C
               DO 770 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I - 1) + 1
                  MALI = IT1AO(ISYMAL,ISYMI)
     *                 + NBAS(ISYMAL)*(I - 1) + 1
C
                  KOFF1 = IT2BGT(ISYMK,ISYALI)
     *                  + NRHF(ISYMK)*(MALI - 1) + 1
                  KOFF2 = ILMVIR(ISYMA) + 1
                  KOFF3 = IT2BCT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),
     *                       NBAS(ISYMAL),ONE,SCR3(KOFF1),NRHFK,
     *                       XLAMDP(KOFF2),NBASAL,
     *                       ONE,SCR1(KOFF3),NRHFK)
C
  770          CONTINUE
  760       CONTINUE
  750    CONTINUE
C
      ENDIF
C
C---------------------------------------
C     Dump to disk the new contribution.
C---------------------------------------
C
C
      IF ( ICON .EQ. 2 ) THEN
         IOFF = IT2DEL(IDEL) + 1
      ELSE
         IOFF = IT2DLR(IDEL,IV) + 1
      ENDIF
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUD,DFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
      IF (IOPTR12.EQ.1) THEN
        CALL DAXPY(NT2BCD(ISYAIK),FACTD,SCR4,1,SCR1,1)
        IF (NT2BCD(ISYAIK) .GT. 0) THEN
          CALL PUTWA2(LUDP,DPFIL,SCR1,IOFF,NT2BCD(ISYAIK))
        END IF
      END IF
C
  999 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_c */
      SUBROUTINE CCRHS_C(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,
     *                   XLAMDP,XLAMIP,XLAMDH,
     *                   XLAMPC,ISYMPC,XLAMHC,ISYMHC,SCRM,E1PIM,
     *                   WORK,LWORK,IDEL,ISYMD,FACTC,ICON,IOPTR12,
     *                   IOPTE,LUC,CFIL,LUCP,CPFIL,IV)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 27-July-1994
C     Generalisation for CCLR by Ove Christiansen august-september 1995
C     (right transformation) and september 1996 (F-matrix).
C
C     Extended for CCSDR12, C. Neiss spring 2006
C     IOPTR12 = 1 Calculate both conv. C and r12 C' intermediates;
C                 T2-dependent contr. to C' interm. is added with a prefactor
C                 of 2*FACTC
C     IOPTE   = 1 Calculate the T-dependent part of the 
C                 E_{a delta')^1' intermediate (on E1PIM).
C
C     Purpose: Calculate C-term.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION XINT(*),DSRHF(*),OMEGA2(*),XLAMDH(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMIP(*),SCRM(*),XLAMPC(*),XLAMHC(*)
      DIMENSION T2AM(*),E1PIM(*)
      CHARACTER CFIL*(*),CPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C--------------------------------------
C     Dynamic allocation of work space.
C--------------------------------------
C
      KSCR1 = 1
      KSCR2 = KSCR1 + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3 = KSCR2 + NT2BGD(ISYDIS)
      IF (ICON .EQ. 2) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF
      IF (IOPTR12.EQ.1) THEN
         KSCR4  = KEND1
         KEND1  = KSCR4  + NT2BCD(ISYAIK)
      END IF
      
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in CCRHS_C')
      END IF
C
C--------------------------------------
C     Transpose the cluster amplitudes.
C--------------------------------------
C
      IF (ICON .GE. 2) THEN
         IF (.NOT. T2TCOR) THEN
            CALL CCSD_T2TP(T2AM,WORK(KEND1),LWRK1,ISYMT2)
         ENDIF
         IF (.NOT. DUMPCD) CALL CCSD_T2MTP(SCRM,WORK(KEND1),LWRK1,ISYMD)
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      IF (.NOT. CC2) THEN
         CALL CCRHS_C1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,SCRM,E1PIM,
     *                 WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),WORK(KSCR4),
     *                 XLAMDP,XLAMIP,XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,
     *                 WORK(KEND1),LWRK1,
     *                 ISYDIS,IDEL,ISYMD,FACTC,ICON,IOPTR12,IOPTE,
     *                 LUC,CFIL,LUCP,CPFIL,IV)
      ENDIF
C
C--------------------------------------
C     Transpose the cluster amplitudes.
C--------------------------------------
C
      IF (ICON .GE. 2) THEN
         IF (.NOT. T2TCOR) THEN
            CALL CCSD_T2TP(T2AM,WORK(KEND1),LWRK1,ISYMT2)
         ENDIF
         IF (.NOT. DUMPCD) CALL CCSD_T2MTP(SCRM,WORK(KEND1),LWRK1,ISYMD)
      ENDIF
C
      RETURN
      END
      SUBROUTINE CCRHS_C1(XINT,DSRHF,OMEGA2,T2AM,ISYMT2,SCRM,E1PIM,
     *                    SCR1,SCR2,SCR3,SCR4,XLAMDP,XLAMIP,
     *                    XLAMDH,XLAMPC,ISYMPC,XLAMHC,ISYMHC,WORK,
     *                    LWORK,ISYDIS,IDEL,ISYDEL,FACTC,ICON,IOPTR12,
     *                    IOPTE,LUC,CFIL,LUCP,CPFIL,IV)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 27-July-1994
C
C     modification by Ove Christiansen 25-7-1995 to allow for a
C     general factor (FACTC) ( assumes DUMCD )
C     and - calculate intermediates for CCLR.
C
C     modification by Ove Christiansen 17-9-1996 for calculating
C     local C-intermediate for F-matrix transformation.
C
C     Thus:
C
C     Modification to calculate terms in 2C1 right transformation in CCLR:
C
C                         1. if icon = 2 both contributions are calculated,
C                            otherwise if ICON = 
C                            1:only the integral (ki | ac)
C                              = (k i-bar | a c) + (k i | a-bar c)
C
C                         3: (k i-bar | a c) + (k i | a-bar c)
C                              + FACTC*Sum(xT*int)
C                                where xT may be non total symmetric.
C
C                         2. Allow for general transformation matrix for
C                            alpha to a(XLAMPC) and for i (XLAMHC).
C                            (the extra i transformation introduces new
C                             blocks which is only entered when 
C                             icon =1 or 3)
C
C                         3. If icon diff. from 2 (we have linear response)
C                            The C intermediate is stored according to
C                            the number of simultaneous trial vector
C                            given by IV. This is ensured using IT2DLR.
C
C     Thus in energy calc: icon = 2,fact = 1/2
C     For right transformation: 
C         icon=1,fact=anything, iv = current vector being transformed
C     For F-matrix transformation:
C         icon=3,fact=1.0, NB - not implemented several vectors yet.
C
C     extended for CCSDR12, C. Neiss spring 2006
C
C     Purpose: Calculate C-term intermediate.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccsdinp.h"
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,XMHALF=-0.5D0)
      PARAMETER (TWO=2.0D0)
      DIMENSION XINT(*),OMEGA2(*),T2AM(*),DSRHF(*)
      DIMENSION SCRM(*),E1PIM(*)
      DIMENSION SCR1(*), SCR2(*), SCR3(*), SCR4(*)
      DIMENSION XLAMDP(*),XLAMIP(*),XLAMDH(*),XLAMPC(*),XLAMHC(ISYMHC)
      DIMENSION WORK(LWORK)
      INTEGER   NADP(8),IADP(8,8),IBASX(8)
      CHARACTER CFIL*(*),CPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "r12int.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (IOPTE.EQ.1) THEN
        IF (.NOT.CCR12) CALL QUIT('IOPTE only implemented for CC-R12')
        IBASX(1) = 0
        DO ISYM = 2, NSYM
          IBASX(ISYM) = IBASX(ISYM-1) + MBAS2(ISYM-1)
        END DO
        DO ISYM = 1, NSYM
          NADP(ISYM) = 0
          DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            IADP(ISYM1,ISYM2) = NADP(ISYM)
            NADP(ISYM) = NADP(ISYM) + 
     &                   NVIR(ISYM1)*(MBAS1(ISYM2)+MBAS2(ISYM2))
          END DO
        END DO
      END IF  
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      ISAIK2 = MULD2H(ISYDIS,ISYMT2)
      IF (ISYAIK .NE. ISAIK2) CALL QUIT('Symmetry mismatch in CCRHS_C')
C
C-------------------------------------------------------
C     Calculate the integrals K(k,dl) = (k d | l delta).
C-------------------------------------------------------
C
      IF (ICON .GE. 2) THEN
C
         DO 100 ISYML = 1,NSYM
C
            ISYMAG = MULD2H(ISYML,ISYDIS)
C
            DO 110 L = 1,NRHF(ISYML)
C
               KSCR10 = 1
               KEND1  = KSCR10 + N2BST(ISYMAG)
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_C1')
               END IF
C
               KOFF1 = IDSRHF(ISYMAG,ISYML) + NNBST(ISYMAG)*(L-1) + 1
               CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK(KSCR10))
C
               DO 120 ISYMDL = 1,NSYM
C
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  ISYMK = MULD2H(ISYMDL,ISYDIS)
                  ISYMA = ISYMK
                  ISYMG = ISYMD
C
                  NBASA = MAX(NBAS(ISYMA),1)
                  NBASG = MAX(NBAS(ISYMG),1)
                  NRHFK = MAX(NRHF(ISYMK),1)
C
                  KSCR11 = KEND1
                  KEND2  = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_C1')
                  END IF
C
                  KOFF2 = ILMRHF(ISYMK) + 1
                  KOFF3 = IAODIS(ISYMA,ISYMG) + 1
C
                  CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),
     *                       NBAS(ISYMA),ONE,XLAMDP(KOFF2),NBASA,
     *                       WORK(KOFF3),NBASA,
     *                       ZERO,WORK(KSCR11),NRHFK)
C
                  NDL   = IT1AM(ISYMD,ISYML)
     *                  + NVIR(ISYMD)*(L - 1) + 1
                  KOFF5 = ILMVIR(ISYMD) + 1
                  KOFF6 = IT2BCT(ISYMK,ISYMDL)
     *                  + NRHF(ISYMK)*(NDL - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMD),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMDH(KOFF5),NBASG,
     *                       ZERO,SCR1(KOFF6),NRHFK)
C
  120          CONTINUE
C
  110       CONTINUE
C
  100    CONTINUE
C
C-----------------------------------------
C        Calculate the first contribution.
C        Sum(dl)T(al,di)*I(lckd)
C-----------------------------------------
C
         IF (LWORK .LT. NT2BCD(ISYAIK)) THEN
            CALL QUIT('Insufficient work space in CCRHS_C1')
         ENDIF
C
         DO 200 ISYMK  = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIK,ISYMK)
            ISYMDL = MULD2H(ISYDIS,ISYMK)
C
            NRHFK  = MAX(NRHF(ISYMK),1)
            NTOTDL = MAX(NT1AM(ISYMDL),1)
C
            KOFF1 = IT2BCT(ISYMK,ISYMDL) + 1
            KOFF2 = IT2SQ(ISYMDL,ISYMAI) + 1
            KOFF3 = IT2BCT(ISYMK,ISYMAI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NT1AM(ISYMAI),NT1AM(ISYMDL),
     *                 ONE,SCR1(KOFF1),NRHFK,T2AM(KOFF2),NTOTDL,ZERO,
     *                 WORK(KOFF3),NRHFK)
C
  200    CONTINUE
C
         CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR1,1)
C
         !save a copy for first contribution:
         IF (IOPTR12.EQ.1) THEN
           CALL DCOPY(NT2BCD(ISYAIK),WORK,1,SCR4,1)
         END IF
C
      ENDIF
C
C----------------------------------------------------------
C     Calculate the integrals K(k,ai) = (k i | alfa delta).
C----------------------------------------------------------
C
      DO 300 ISYMA = 1,NSYM
C
         ISYMBG = MULD2H(ISYMA,ISYDIS)
C
         KSCR10 = 1
         KEND1  = KSCR10 + N2BST(ISYMBG)
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space for allocation in CCRHS_C1')
         END IF
C
         DO 310 A = 1,NBAS(ISYMA)
C
            KOFF1 = IDSAOG(ISYMA,ISYDIS) + NNBST(ISYMBG)*(A - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMBG,WORK(KSCR10))
C
            DO 320 ISYMG = 1,NSYM
C
               ISYMI  = ISYMG
               ISYMB  = MULD2H(ISYMG,ISYMBG)
               ISYMK  = ISYMB
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NBASG = MAX(NBAS(ISYMG),1)
               NRHFK = MAX(NRHF(ISYMK),1)
C
               KSCR11 = KEND1
               KSCR12 = KSCR11 + NRHF(ISYMK)*NBAS(ISYMG)
               KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Not enough space for '//
     &                 'allocation in CCRHS_C1')
               END IF
C
               KOFF2 = ILMRHF(ISYMK) + 1
               KOFF3 = KSCR10 + IAODIS(ISYMB,ISYMG)
C
               CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMB),
     *                    ONE,XLAMDP(KOFF2),NBASB,WORK(KOFF3),NBASB,
     *                    ZERO,WORK(KSCR11),NRHFK)
C
               KOFF5 = ILMRHF(ISYMI) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KSCR11),NRHFK,XLAMDH(KOFF5),NBASG,
     *                    ZERO,WORK(KSCR12),NRHFK)
C
C
               DO 330 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
                  KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                  CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR2(KOFF8),1)
C
  330          CONTINUE
C
C
C-------------------------------------------------------
C              In 2C1 linear transformation extra  cont.
C-------------------------------------------------------
C
               IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
                  ISYMI  = MULD2H(ISYMG,ISYMHC)
                  ISYMAI = MULD2H(ISYMA,ISYMI)
C
                  KEND2  = KSCR12 + NRHF(ISYMK)*NRHF(ISYMI)
                  LWRK2  = LWORK  - KEND2
                  IF (LWRK2 .LT. 0) THEN
                     CALL QUIT('Not enough space for '//
     &                    'allocation in CCRHS_D1')
                  END IF
C
                  KOFF5 = IGLMRH(ISYMG,ISYMI) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NBAS(ISYMG),ONE,WORK(KSCR11),NRHFK,
     *                       XLAMHC(KOFF5),NBASG,
     *                       ZERO,WORK(KSCR12),NRHFK)
C
                  DO 331 I = 1,NRHF(ISYMI)
C
                     NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                     KOFF8 = IT2BGT(ISYMK,ISYMAI)
     *                     + NRHF(ISYMK)*(NAI - 1) + 1
                     KOFF7 = KSCR12 + NRHF(ISYMK)*(I - 1)
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF7),1,SCR3(KOFF8),1)
C
  331             CONTINUE
C
               ENDIF
C
  320       CONTINUE
C
  310    CONTINUE
C
  300 CONTINUE
C
      IF (DUMPCD) GOTO 800
C
      IF (CCR12) CALL QUIT('CCSDR12 requires DUMPCD=.TRUE.')
C
C-----------------------------------------
C     Back transformation to the AO basis.
C-----------------------------------------
C
      DO 400 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYDIS)
C
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 410 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            NBASA = MAX(NBAS(ISYMA),1)
C
            DO 420 I = 1,NRHF(ISYMI)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I - 1) + 1
C
               KOFF1 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
               KOFF2 = ILMVIR(ISYMA) + 1
               KOFF3 = IT2BGT(ISYMK,ISYMAI) + NRHF(ISYMK)*(MAI - 1) + 1
C
               CALL DGEMM('N','T',NRHF(ISYMK),NBAS(ISYMA),NVIR(ISYMA),
     *                    XMHALF,SCR1(KOFF1),NRHFK,XLAMIP(KOFF2),NBASA,
     *                    ONE,SCR2(KOFF3),NRHFK)
C
  420       CONTINUE
C
  410    CONTINUE
C
  400 CONTINUE
C
C
      DO 500 ISYMK = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         DO 510 K = 1,NRHF(ISYMK)
C
            DO 520 ISYMJ = 1,NSYM
C
               ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
               NBASB = MAX(NBAS(ISYMB),1)
               NVIRB = MAX(NVIR(ISYMB),1)
C
               KOFF1 = ILMVIR(ISYMB) + 1
               KOFF2 = IT2BCD(ISYMBJ,ISYMK) + NT1AM(ISYMBJ)*(K - 1)
     *               + IT1AM(ISYMB,ISYMJ) + 1
               KOFF3 = IT2BGD(ISYMBJ,ISYMK) + NT1AO(ISYMBJ)*(K - 1)
     *               + IT1AO(ISYMB,ISYMJ) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMB),NRHF(ISYMJ),NVIR(ISYMB),
     *                    ONE,XLAMIP(KOFF1),NBASB,SCRM(KOFF2),NVIRB,
     *                    ZERO,SCR3(KOFF3),NBASB)
C
  520       CONTINUE
C
  510    CONTINUE
C
  500 CONTINUE
C
C---------------------------------------
C     Calculate the second contribution.
C
C     Alfredo will introduce the batching over ai before the
C     end of august 1994.
C---------------------------------------
C
      DO 600 ISYMAI = 1,NSYM
C
         ISYMK  = MULD2H(ISYMAI,ISYDIS)
         ISYMBJ = MULD2H(ISYMK,ISYDEL)
C
         IF (NRHF(ISYMK) .EQ. 0) GOTO 600
C
         IF (LWORK .LT. NT1AO(ISYMBJ)) THEN
            CALL QUIT('Insufficient work space in CCRHS_C1')
         ENDIF
C
         NTOTBJ = MAX(NT1AO(ISYMBJ),1)
C
         DO 610 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            DO 620 I = 1,NRHF(ISYMI)
C
               DO 630 A = 1,NBAS(ISYMA)
C
                  NAI = IT1AO(ISYMA,ISYMI) + NBAS(ISYMA)*(I-1) + A
C
                  KOFF1 = IT2BGD(ISYMBJ,ISYMK) + 1
                  KOFF2 = IT2BGT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMV('N',NT1AO(ISYMBJ),NRHF(ISYMK),ONE,
     *                       SCR3(KOFF1),NTOTBJ,SCR2(KOFF2),1,
     *                       ZERO,WORK,1)
C
                  IF (.NOT. OMEGSQ) THEN
C
C
                  IF (ISYMAI .EQ. ISYMBJ) THEN
                     WORK(NAI) = TWO*WORK(NAI)
                  ENDIF
C
                  DO 640 ISYMJ = 1,NSYM
C
                     ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     DO 650 J = 1,NRHF(ISYMJ)
C
                        NAJ = IT1AO(ISYMA,ISYMJ) + NBAS(ISYMA)*(J-1) + A
C
                        DO 660 B = 1,NBAS(ISYMB)
C
                           NBI = IT1AO(ISYMB,ISYMI)
     *                         + NBAS(ISYMB)*(I-1) + B
                           NBJ = IT1AO(ISYMB,ISYMJ)
     *                         + NBAS(ISYMB)*(J-1) + B
C
                           NAIBJ = IT2AO(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                           NAJBI = IT2AO(ISYMAJ,ISYMBI) + INDEX(NAJ,NBI)
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)-HALF*WORK(NBJ)
                           OMEGA2(NAJBI) = OMEGA2(NAJBI)-WORK(NBJ)
C
  660                   CONTINUE
  650                CONTINUE
  640             CONTINUE
C
C
                  ELSE
C
C
                  KOFF = IT2AOS(ISYMBJ,ISYMAI)
     *                 + NT1AO(ISYMBJ)*(NAI - 1) + 1
                  CALL DAXPY(NT1AO(ISYMBJ),-HALF,WORK,1,OMEGA2(KOFF),1)
C
                  DO 740 ISYMJ = 1,NSYM
C
                     ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     NBI = IT1AO(ISYMB,ISYMI) + NBAS(ISYMB)*(I-1) + 1
      
C
                     DO 750 J = 1,NRHF(ISYMJ)
C
                        NAJ = IT1AO(ISYMA,ISYMJ) + NBAS(ISYMA)*(J-1) + A
                        NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J-1) + 1
C
                        NBIAJ = IT2AOS(ISYMBI,ISYMAJ)
     *                        + NT1AO(ISYMBI)*(NAJ - 1) + NBI
C
                        CALL DAXPY(NBAS(ISYMB),-ONE,WORK(NBJ),1,
     *                             OMEGA2(NBIAJ),1)
C
  750                CONTINUE
  740             CONTINUE
C
                  ENDIF
C
  630          CONTINUE
  620       CONTINUE
C
  610    CONTINUE
  600 CONTINUE
C
      GOTO 999
C
C-------------------
C     I/O algorithm.
C-------------------
C
  800 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C-----------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      IF ( ICON .EQ. 1 ) CALL DZERO(SCR1,NT2BCD(ISYAIK))
C
      DO 810 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 820 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMAI)
            ISYMAL= MULD2H(ISYMPC,ISYMA)
            ISYALI= MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 830 I = 1,NRHF(ISYMI)
C
               NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
               MALI = IT1AO(ISYMAL,ISYMI) + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI- 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTC,SCR1(KOFF3),NRHFK)
C
               IF (IOPTE.EQ.1) THEN
                 IF (ISYMI.EQ.ISYMK) THEN
                   KOFF3 = IT2BCT(ISYMK,ISYMAI) + 
     &                     NRHF(ISYMK)*(NAI - 1) + I
                   IF (IDEL.GT.NBAST) THEN
                     D = IDEL-IBASX(ISYDEL)-NBAST+MBAS1(ISYDEL)
                   ELSE
                     D = IDEL-IBAS(ISYDEL)
                   END IF
C                  WRITE(LUPRI,*)'ISYDEL, IDEL, D:',ISYDEL, IDEL, D
                   KOFFE = IADP(ISYMA,ISYDEL) + 
     &                     NVIR(ISYMA)*(D-1) + 1
                   CALL DAXPY(NVIR(ISYMA),1.5D0,SCR1(KOFF3),NRHF(ISYMK),
     &                        E1PIM(KOFFE),1)
                 END IF
               END IF
C
  830       CONTINUE
  820    CONTINUE
  810 CONTINUE
C
C-----------------------------------------------
C     Transform the alpha index of K(k,ai) to a.
C     I is C1 transformed.
C-----------------------------------------------
C
      IF ((ICON .EQ. 3) .OR. (ICON .EQ. 1)) THEN
C
         ISYAIK = MULD2H(ISYDIS,ISYMHC)
C
         DO 850 ISYMAI = 1,NSYM
C
            ISYMK = MULD2H(ISYMAI,ISYAIK)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            DO 860 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAL= ISYMA
               ISYALI= MULD2H(ISYMAL,ISYMI)
               NBASAL = MAX(NBAS(ISYMAL),1)
C
               DO 870 I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I - 1) + 1
                  MALI = IT1AO(ISYMAL,ISYMI)
     *                 + NBAS(ISYMAL)*(I - 1) + 1
C
                  KOFF1 = IT2BGT(ISYMK,ISYALI)
     *                  + NRHF(ISYMK)*(MALI - 1) + 1
                  KOFF2 = ILMVIR(ISYMA) + 1
                  KOFF3 = IT2BCT(ISYMK,ISYMAI)
     *                  + NRHF(ISYMK)*(NAI - 1) + 1
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),
     *                       NBAS(ISYMAL),ONE,SCR3(KOFF1),NRHFK,
     *                       XLAMDP(KOFF2),NBASAL,
     *                       ONE,SCR1(KOFF3),NRHFK)
C
  870          CONTINUE
  860       CONTINUE
  850    CONTINUE
C
      ENDIF
C---------------------------------------------------------
C     Dump to disk the new contribution.
C     energy calc icon = 2
C     rsp calc. write to position given by it2dlr(idel,iv)
C---------------------------------------------------------
C
      IF ( ICON .EQ. 2 ) THEN
C
         IOFF = IT2DEL(IDEL) + 1
C
      ELSE
C
         IOFF = IT2DLR(IDEL,IV) + 1
C
      ENDIF
C
      IF (NT2BCD(ISYAIK) .GT. 0) THEN
         CALL PUTWA2(LUC,CFIL,SCR1,IOFF,NT2BCD(ISYAIK))
      ENDIF
C
      IF (IOPTR12.EQ.1) THEN
        CALL DAXPY(NT2BCD(ISYAIK),FACTC,SCR4,1,SCR1,1)
        IF (NT2BCD(ISYAIK) .GT. 0) THEN
          CALL PUTWA2(LUCP,CPFIL,SCR1,IOFF,NT2BCD(ISYAIK))
        END IF
      END IF
C
  999 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_gam */
      SUBROUTINE CCRHS_GAM(DSRHF,GAMMA,XLAMDP,XLAMDH,SCRM,
     *                     WORK,LWORK,IDEL,ISYMD)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 21-July-1994
C
C     Purpose: Calculate the gamma intermediate.
C
#include "implicit.h"
      DIMENSION DSRHF(*),GAMMA(*),SCRM(*)
      DIMENSION WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KLAMDA = 1
      KEND1  = KLAMDA + NRHF(ISYMD)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_GAM')
      ENDIF
C
C---------------------------------------
C     Copy XLAMDH vector for given IDEL.
C---------------------------------------
C
      KOFF1 = ILMRHF(ISYMD) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NRHF(ISYMD),XLAMDH(KOFF1),NBAS(ISYMD),WORK(KLAMDA),1)
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMAG = MULD2H(ISYML,ISYDIS)
C
C---------------------------
C        Dynamic allocation.
C---------------------------
C
         KSCR1  = KEND1
         KSCR2  = KSCR1  + N2BST(ISYMAG)
         KSCR3  = KSCR2  + NT1AO(ISYMAG)
         KSCR4  = KSCR3  + NT1AM(ISYMAG)
         KEND2  = KSCR4  + NMATIJ(ISYMAG)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
            CALL QUIT('Insufficient space in CCRHS_GAM')
         ENDIF
C
         CALL CCRHS_GAM1(DSRHF,GAMMA,SCRM,WORK(KLAMDA),XLAMDP,XLAMDH,
     *                   WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),
     *                   WORK(KSCR4),WORK(KEND2),LWRK2,ISYMD,ISYML,
     *                   ISYMAG)
C
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CCRHS_GAM1(DSRHF,GAMMA,SCRM,XLAM,
     *              XLAMDP,XLAMDH,SCR1,SCR2,SCR3,SCR4,WORK,
     *              LWORK,ISYMD,ISYML,ISYMAG)
C
C     Written by Henrik Koch 3-Jan-1994
C
C     Purpose: Calculate the gamma intermediate.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DSRHF(*),GAMMA(*),SCRM(*),XLAM(*)
      DIMENSION SCR1(*),SCR2(*),SCR3(*),SCR4(*),WORK(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMKC = ISYMAG
C
      DO 100 L = 1,NRHF(ISYML)
C
         KOFF1 = IDSRHF(ISYMAG,ISYML) + NNBST(ISYMAG)*(L - 1) + 1
C
         CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,SCR1)
C
         DO 110 ISYMG = 1,NSYM
C
            ISYMA = MULD2H(ISYMG,ISYMAG)
            ISYMK = ISYMA
            ISYMC = ISYMG
            ISYMI = ISYMG
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASG = MAX(NBAS(ISYMG),1)
            NRHFK = MAX(NRHF(ISYMK),1)
C
            KOFF2 = ILMRHF(ISYMK) + 1
            KOFF3 = IAODIS(ISYMA,ISYMG) + 1
            KOFF4 = IT1AOT(ISYMK,ISYMG) + 1
C
            CALL DGEMM('T','N',NRHF(ISYMK),NBAS(ISYMG),NBAS(ISYMA),
     *                 ONE,XLAMDP(KOFF2),NBASA,SCR1(KOFF3),NBASA,
     *                 ZERO,SCR2(KOFF4),NRHFK)
C
            KOFF5 = ILMVIR(ISYMC) + 1
            KOFF6 = IT1AMT(ISYMK,ISYMC) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMC),NBAS(ISYMG),
     *                 ONE,SCR2(KOFF4),NRHFK,XLAMDH(KOFF5),NBASG,
     *                 ZERO,SCR3(KOFF6),NRHFK)
C
            KOFF7 = ILMRHF(ISYMI) + 1
            KOFF8 = IMATIJ(ISYMK,ISYMI) + 1
C
            CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NBAS(ISYMG),
     *                 ONE,SCR2(KOFF4),NRHFK,XLAMDH(KOFF7),NBASG,
     *                 ZERO,SCR4(KOFF8),NRHFK)
C
  110    CONTINUE
C
         DO 120 ISYMJ = 1,NSYM
C
            ISYMLJ = MULD2H(ISYML,ISYMJ)
            ISYMKI = MULD2H(ISYMLJ,ISYMOP)
            ISYMCI = MULD2H(ISYMJ,ISYMD)
C
            KSCR5 = 1
            KEND1 = KSCR5 + NMATIJ(ISYMKI)
C
            IF (ISYMKI .GT. ISYMLJ) GOTO 120
C
            DO 130 J = 1,NRHF(ISYMJ)
C
               DO 140 ISYMI = 1,NSYM
C
                  ISYMC = MULD2H(ISYMI,ISYMCI)
                  ISYMK = MULD2H(ISYMI,ISYMKI)
C
                  NVIRC = MAX(NVIR(ISYMC),1)
                  NRHFK = MAX(NRHF(ISYMK),1)
C
                  KOFF2 = IT1AMT(ISYMK,ISYMC) + 1
                  KOFF3 = IT2BCD(ISYMCI,ISYMJ)
     *                  + NT1AM(ISYMCI)*(J - 1)
     *                  + IT1AM(ISYMC,ISYMI) + 1
                  KOFF4 = KSCR5 + IMATIJ(ISYMK,ISYMI)
C
                  CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),
     *                       NVIR(ISYMC),ONE,SCR3(KOFF2),NRHFK,
     *                       SCRM(KOFF3),NVIRC,ZERO,WORK(KOFF4),NRHFK)
C
  140          CONTINUE
C
               IF (ISYMJ .EQ. ISYMD) THEN
                  CALL DAXPY(NMATIJ(ISYMKI),XLAM(J),SCR4,1,
     *                       WORK(KSCR5),1)
               ENDIF
C
               NLJ = IMATIJ(ISYML,ISYMJ) + NRHF(ISYML)*(J - 1) + L
C
               IF (ISYMOP .EQ. 1) THEN
                  KKILJ = IGAMMA(ISYMKI,ISYMLJ) + NLJ*(NLJ-1)/2
                  DO 150 NKI = 1,NLJ
C
                     KOFF = KSCR5 + NKI - 1
                     NKILJ = KKILJ + NKI
                     GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C
  150             CONTINUE
               ELSE
                  KKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                  + NMATIJ(ISYMKI)*(NLJ - 1)
                  DO 160 NKI = 1,NMATIJ(ISYMKI)
C
                     KOFF = KSCR5 + NKI - 1
                     NKILJ = KKILJ + NKI
                     GAMMA(NKILJ) = GAMMA(NKILJ) + WORK(KOFF)
C
  160             CONTINUE
               END IF
C
  130       CONTINUE
  120    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_b */
      SUBROUTINE CCRHS_B(XINT,OMEGA2,XLAMDP,XLAMDH,SCRM,
     *                   WORK,LWORK,IDEL,ISYMD)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 18-July-1994
C
C     Purpose: Calculate B-term.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),OMEGA2(*),XLAMDH(*),WORK(LWORK)
      DIMENSION XLAMDP(*),SCRM(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KMGD   = 1
      KEND1  = KMGD   + NT2BGD(ISYMD)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_B')
      ENDIF
C
C-----------------------------
C     Prepare the data arrays.
C-----------------------------
C
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMCI = MULD2H(ISYMJ,ISYMD)
         ISYMGI = ISYMCI
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMC = MULD2H(ISYMI,ISYMCI)
            ISYMG = ISYMC
C
            NVIRC = MAX(NVIR(ISYMC),1)
            NBASG = MAX(NBAS(ISYMG),1)
C
            KOFF1 = ILMVIR(ISYMC) + 1
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               KOFF2 = IT2BCD(ISYMCI,ISYMJ) + IT1AM(ISYMC,ISYMI)
     *               + NT1AM(ISYMCI)*(J - 1) + 1
               KOFF3 = IT2BGD(ISYMGI,ISYMJ) + IT1AO(ISYMG,ISYMI)
     *               + NT1AO(ISYMGI)*(J - 1) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC),
     *                    ONE,XLAMDH(KOFF1),NBASG,SCRM(KOFF2),NVIRC,
     *                    ZERO,WORK(KOFF3),NBASG)
C
               IF (ISYMG .EQ. ISYMD) THEN
                  KOFF4 = KOFF3 + IDEL - IBAS(ISYMD) - 1
                  CALL DSCAL(NRHF(ISYMI),HALF,WORK(KOFF4),NBAS(ISYMG))
               ENDIF
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_B1(XINT,OMEGA2,WORK(KMGD),WORK(KEND1),LWRK1,IDEL,ISYMD)
C
      RETURN
      END
      SUBROUTINE CCRHS_B1(XINT,OMEGA2,XMGD,WORK,LWORK,IDEL,ISYMD)
C
C     Written by Henrik Koch 3-Jan-1994
C
C     Purpose: Calculate B-term.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),OMEGA2(*),XMGD(*)
      DIMENSION WORK(LWORK)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYDIS = MULD2H(ISYMOP,ISYMD)
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      IF (OMEGSQ) GOTO 200
C
      DO 100 ISYMB = 1,NSYM
C
         ISYMAG = MULD2H(ISYMB,ISYDIS)
C
         KSCR1 = 1
         KEND1 = KSCR1 + N2BST(ISYMAG)
         LWRK1 = LWORK - KEND1
C
         DO 110 B = 1,NBAS(ISYMB)
C
            KOFF1 = IDSAOG(ISYMB,ISYDIS) + NNBST(ISYMAG)*(B - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAG,WORK(KSCR1))
C
            DO 120 ISYMJ = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYMOP)
               ISYMGI = MULD2H(ISYMJ,ISYMD)
C
               KSCR2 = KEND1
               KEND2 = KSCR2 + NT1AO(ISYMAI)
               LWRK2 = LWORK - KEND2
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  CALL DZERO(WORK(KSCR2),NT1AO(ISYMAI))
C
                  DO 140 ISYMI = 1,NSYM
C
                     ISYMG  = MULD2H(ISYMI,ISYMGI)
C
                     IF (ISYMG .GT. ISYMD) GOTO 140
C
                     ISYMA  = MULD2H(ISYMG,ISYMAG)
C
                     KOFF2 = KSCR1 + IAODIS(ISYMA,ISYMG)
                     KOFF3 = IT2BGD(ISYMGI,ISYMJ)
     *                     + NT1AO(ISYMGI)*(J - 1) + 1
     *                     + IT1AO(ISYMG,ISYMI)
                     KOFF4 = KSCR2 + IT1AO(ISYMA,ISYMI)
C
                     NBASA = MAX(NBAS(ISYMA),1)
                     NBASG = MAX(NBAS(ISYMG),1)
C
                     IF (ISYMG .LT. ISYMD) THEN
                        NTOTG = NBAS(ISYMG)
                     ELSE
                        NTOTG = IDEL - IBAS(ISYMD)
                     ENDIF
C
                     CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),
     *                          NTOTG,ONE,WORK(KOFF2),NBASA,
     *                          XMGD(KOFF3),NBASG,ZERO,WORK(KOFF4),
     *                          NBASA)
C
  140             CONTINUE
C
C---------------------------------------
C                 Accumulate the result.
C---------------------------------------
C
                  NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J - 1) + B
C
                  IF (ISYMAI .EQ. ISYMBJ) THEN
                     WORK(KSCR2+NBJ-1) = TWO*WORK(KSCR2+NBJ-1)
                  ENDIF
C
                  DO 150 NAI = 1,NT1AO(ISYMAI)
                     KOFF5 = KSCR2 + NAI - 1
                     NAIBJ = IT2AO(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KOFF5)
  150             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
C
  200 CONTINUE
C
      DO 300 ISYMB = 1,NSYM
C
         ISYMAG = MULD2H(ISYMB,ISYDIS)
C
         KSCR1 = 1
         KEND1 = KSCR1 + N2BST(ISYMAG)
         LWRK1 = LWORK - KEND1
C
         DO 310 B = 1,NBAS(ISYMB)
C
            KOFF1 = IDSAOG(ISYMB,ISYDIS) + NNBST(ISYMAG)*(B - 1) + 1
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAG,WORK(KSCR1))
C
            DO 320 ISYMJ = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYMOP)
               ISYMGI = MULD2H(ISYMJ,ISYMD)
C
               KSCR2 = KEND1
               KEND2 = KSCR2 + NT1AO(ISYMAI)
               LWRK2 = LWORK - KEND2
C
               DO 330 J = 1,NRHF(ISYMJ)
C
                  NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J - 1) + B
C
                  DO 340 ISYMI = 1,NSYM
C
                     ISYMG  = MULD2H(ISYMI,ISYMGI)
C
                     IF (ISYMG .GT. ISYMD) GOTO 340
C
                     ISYMA  = MULD2H(ISYMG,ISYMAG)
C
                     KOFF2 = KSCR1 + IAODIS(ISYMA,ISYMG)
                     KOFF3 = IT2BGD(ISYMGI,ISYMJ)
     *                     + NT1AO(ISYMGI)*(J - 1) + 1
     *                     + IT1AO(ISYMG,ISYMI)
C
                     KOFF4 = IT2AOS(ISYMAI,ISYMBJ)
     *                     + NT1AO(ISYMAI)*(NBJ - 1)
     *                     + IT1AO(ISYMA,ISYMI) + 1
C
                     NBASA = MAX(NBAS(ISYMA),1)
                     NBASG = MAX(NBAS(ISYMG),1)
C
                     IF (ISYMG .LT. ISYMD) THEN
                        NTOTG = NBAS(ISYMG)
                     ELSE
                        NTOTG = IDEL - IBAS(ISYMD)
                     ENDIF
C
                     CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),
     *                          NTOTG,ONE,WORK(KOFF2),NBASA,
     *                          XMGD(KOFF3),NBASG,ONE,OMEGA2(KOFF4),
     *                          NBASA)
C
  340             CONTINUE
C
  330          CONTINUE
  320       CONTINUE
C
  310    CONTINUE
  300 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_f */
      SUBROUTINE CCRHS_F(XINT,OMEGA2,XLAMDH,WORK,LWORK,IDEL,ISYMD)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 13-July-1994
C
C     Purpose: Calculate F-term.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (HALF = 0.5D0)
      DIMENSION XINT(*),OMEGA2(*)
      DIMENSION XLAMDH(*),WORK(LWORK)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
C
      ISYMJ  = ISYMD
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KLAMDA = 1
      KEND1  = KLAMDA + NRHF(ISYMJ)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_F')
      ENDIF
C
C---------------------------------------
C     Copy XLAMDH vector for given IDEL.
C---------------------------------------
C
      KOFF = ILMRHF(ISYMJ) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NRHF(ISYMD),XLAMDH(KOFF),NBAS(ISYMD),WORK(KLAMDA),1)
C
      IF (OMEGSQ) THEN
         CALL DSCAL(NRHF(ISYMD),HALF,WORK(KLAMDA),1)
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 ISYMB = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYMOP)
C
         IF (ISYMAI .GT. ISYMBJ) GOTO 100
C
         KOFF1 = IDSAOG(ISYMB,ISYDIS) + 1
C
         IF (.NOT. OMEGSQ) THEN
            KOFF2 = IT2AO(ISYMAI,ISYMBJ) + 1
         ELSE
            KOFF2 = IT2AOS(ISYMAI,ISYMBJ) + 1
         ENDIF
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
         KSCR1 = KEND1
         KSCR2 = KSCR1 + N2BST(ISYMAI)
         KEND2 = KSCR2 + NT1AO(ISYMAI)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
            CALL QUIT('Insufficient space in CCRHS_F')
         ENDIF
C
         CALL CCRHS_F1(XINT(KOFF1),OMEGA2(KOFF2),WORK(KLAMDA),
     *                 WORK(KSCR1),WORK(KSCR2), XLAMDH,ISYMJ,
     *                 ISYMB,ISYMAI)
C
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CCRHS_F1(XINT,OMEGA2,XLAM,SCR1,SCR2,XLAMDH,ISYMJ,
     *                    ISYMB,ISYMAI)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 13-July-1994
C
C     Purpose: Calculate F-term.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*),OMEGA2(*),XLAM(*)
      DIMENSION SCR1(*),SCR2(*)
      DIMENSION XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      DO 100 B = 1,NBAS(ISYMB)
C
         KOFF1 = NNBST(ISYMAI)*(B-1) + 1
         CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAI,SCR1)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMG = ISYMI
            ISYMA = MULD2H(ISYMI,ISYMAI)
C
            KOFF2 = IAODIS(ISYMA,ISYMG) + 1
            KOFF3 = ILMRHF(ISYMI) + 1
            KOFF4 = IT1AO(ISYMA,ISYMI) + 1
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASG = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     *                 ONE,SCR1(KOFF2),NBASA,XLAMDH(KOFF3),NBASG,
     *                 ZERO,SCR2(KOFF4),NBASA)
C
  110    CONTINUE
C
         IF (.NOT. OMEGSQ) THEN
            DO 120 J = 1,NRHF(ISYMJ)
C
               NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J-1) + B
C
               IF (ISYMOP .EQ. 1) THEN
                  NTOTAI = NBJ
                  KOFF5  = NBJ*(NBJ - 1)/2 + 1
               ELSE
                  NTOTAI = NT1AO(ISYMAI)
                  KOFF5  = NT1AO(ISYMAI)*(NBJ - 1) + 1
               ENDIF
C
               IF (XLAM(J) .NE. ZERO) THEN
                  CALL DAXPY(NTOTAI,XLAM(J),SCR2,1,OMEGA2(KOFF5),1)
               ENDIF
C
  120       CONTINUE
         ELSE
            DO 130 J = 1,NRHF(ISYMJ)
C
               NBJ = IT1AO(ISYMB,ISYMJ) + NBAS(ISYMB)*(J-1) + B
C
               KOFF5  = NT1AO(ISYMAI)*(NBJ - 1) + 1
C
               CALL DAXPY(NT1AO(ISYMAI),XLAM(J),SCR2,1,OMEGA2(KOFF5),1)
C
  130       CONTINUE
         ENDIF
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cctrbt */
      SUBROUTINE CCTRBT(XINT,DSRHF,XLAMDP,ISYMLP,WORK,LWORK,ISYDIS)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 12-July-1994
C
C     Ove Christiansen 14-6-1996: General sym. lambda matrix ISYMLP
C
C     Purpose: Transform integral batch.
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION XINT(*),DSRHF(*),XLAMDP(*),WORK(LWORK)
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      KOFF1 = 1
      KOFF2 = 1
      KOFF3 = 1
      DO 100 ISYMG = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMLP,ISYMG)
         ISYMAB = MULD2H(ISYMG,ISYDIS)
C
         NNBSAB = MAX(NNBST(ISYMAB),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         KOFF2  = 1 + IGLMRH(ISYMG,ISYMJ)
         KOFF3  = 1 + IDSRHF(ISYMAB,ISYMJ)
C
         CALL DGEMM('N','N',NNBST(ISYMAB),NRHF(ISYMJ),NBAS(ISYMG),
     *              ONE,XINT(KOFF1),NNBSAB,XLAMDP(KOFF2),NBASG,
     *              ZERO,DSRHF(KOFF3),NNBSAB)
C
         KOFF1 = KOFF1 + NNBST(ISYMAB)*NBAS(ISYMG)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrdao */
      SUBROUTINE CCRDAO(XINT,IDELTA,IDEL2,WORK,LWORK,IRECNR,DIRECT)
C
C     Written by Henrik Koch 25-Sep-1993
C
C     Purpose: Read distribution of AO integrals.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxash.h"
#include "iratdef.h"
C
      LOGICAL FIRST, DIRECT
      DIMENSION XINT(*),WORK(LWORK)
      DIMENSION IRECNR(*)
C
      CHARACTER*8 NAME(8)
C
#include "ccorb.h"
C
C   #include "infind.h" replaced by: #include <ccisao.h>
C   (WK/UniKA/28-04-2003).
C
#include "ccisao.h"
C
#include "ccsdsym.h"
#include "cbieri.h"
#include "eribuf.h"
#include "ccpack.h"
#include "r12int.h"
C
      SAVE FIRST
      DATA FIRST /.TRUE./
C
      DATA NAME  /'CCAOIN_1','CCAOIN_2','CCAOIN_3','CCAOIN_4',
     *            'CCAOIN_5','CCAOIN_6','CCAOIN_7','CCAOIN_8'/
      COMMON/SORTIO/LUAOIN(8)
C
      CALL QENTER('CCRDAO')
C
      ISYMD  = ISAO(IDELTA)
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      IF (.NOT. DIRECT) THEN
C
         NFILE = LUAOIN(ISYMD)
         IF (NFILE.LE.0) THEN
           NFILE = 0
           CALL WOPEN2(NFILE,NAME(ISYMD),64,0)
           LUAOIN(ISYMD) = NFILE
         END IF
C
         LENGTH  = NDISAO(ISYDIS)
         NBYTE   = NPCKINT(IDELTA)
         IOFF    = IOFFINT(IDELTA)
         NDWORDS = LENGTH
C
         IF (LPACKINT) NDWORDS = (NBYTE+7)/8
C
         CALL GETWA2(NFILE,NAME(ISYMD),XINT,IOFF,NDWORDS)
C
         IF (LPACKINT) THEN
            DTIME = SECOND()
            CALL UNPCKR8(XINT,LENGTH,XINT,NBYTE,
     &                   IPCKTABINT,LPACKINT)
            PCKTIME = PCKTIME + SECOND() - DTIME
         END IF
C
         GOTO 999
      ENDIF
C
C----------------------------
C     Construct index arrays.
C----------------------------
C
      KADR1 = 1
      KADR2 = KADR1 + (NBAST + 1)/IRAT + 1
      KEND1 = KADR2 + (NBAST*NBAST + 1)/IRAT + 1
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in CCRDAO')
      END IF
C
      CALL CCRD_INIT(WORK(KADR1),WORK(KADR2),ISYDIS)
C
C--------------------
C     Construct XINT.
C--------------------
C
      IF (U21INT) THEN
        CALL DZERO(XINT,2*NDISAO(ISYDIS))
      ELSE
        CALL DZERO(XINT,NDISAO(ISYDIS))
      END IF
C
C     Buffer allocation
C
      KIBUF = KEND1
      KRBUF = KIBUF + (NIBUF*LBUF-1)/IRAT + 1
      KEND2 = KRBUF + LBUF
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CCRDAO')
      ENDIF
C
      CALL CCRDA1(XINT,WORK(KIBUF),WORK(KRBUF),IDELTA,IDEL2,
     *            WORK(KADR1),WORK(KADR2),IRECNR)
C
 999  CONTINUE
      CALL QEXIT('CCRDAO')
      RETURN
      END
C  /* Deck ccrda1 */
      SUBROUTINE CCRDA1(XINT,IBUF4,RBUF,IDELTA,IDEL2,KADR1,KADR2,
     *                  IRECNR)
C
C     Written by Henrik Koch 25-Sep-1993
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "ibtpar.h"
#include "ccorb.h"
#include "mxcent.h"
#include "eribuf.h"
      DIMENSION XINT(*)
      DIMENSION KADR1(NBAST),KADR2(NBAST,NBAST)
      DIMENSION RBUF(LBUF)
      INTEGER*4 IBUF4(LBUF*NIBUF), LENGTH4
      INTEGER   INDX4(4,LBUF)
      DIMENSION IRECNR(*)
      CHARACTER*8 FAODER
      LOGICAL OLDDX
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
#include "ccinftap.h"
#include "nuclei.h"
#include "inftap.h"
#include "eritap.h"
#include "chrnos.h"
#include "r12int.h"
C

C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      IF (NEWDIS) THEN
C
         NEWDIS = .FALSE.
C
         IF (LUINTR .LE. 0) THEN
            CALL GPOPEN(LUINTR,'AOTWODIS','UNKNOWN',' ',
     &        'UNFORMATTED',IDUMMY,.FALSE.)
         END IF
         REWIND (LUINTR)
         DO 50 I = 1,NBUFX(0)
            READ(LUINTR) IRECNR(I)
   50    CONTINUE
C
      ENDIF

      IF (LUAORC(0) .LE. 0) THEN
            LBFINP = LBUF
C     
#if defined (SYS_NEC)
            LRECL =   LBFINP + NIBUF*LBFINP/2 + 1    ! in integer*8 units
#else
            LRECL = 2*LBFINP + NIBUF*LBFINP   + 1    ! in integer*4 units
#endif
            FAODER = 'AO2DIS'//CHRNOS(0)//CHRNOS(0)
            CALL GPOPEN(LUAORC(0),FAODER,'UNKNOWN','DIRECT',
     &           'UNFORMATTED',LRECL,OLDDX)
            IF (U21INT)
     &      CALL GPOPEN(LU21INT,'AOTDIS00','UNKNOWN','DIRECT',
     &           'UNFORMATTED',LRECL,OLDDX)
      END IF
C
      ICOUNT = 0
      IDUM = 1
C
         DO 100 J = 1,NBUFX(0)
C
            IRECJ = IRECNR(J)
            IF (NOAUXB.AND..NOT.LOOPDP) THEN
               IDUM = 1
               CALL IJKAUX(IRECJ,IDUM,IDUM,IDUM)
            END IF
            IF (IRECJ .EQ. IDELTA) THEN
               ICOUNT = ICOUNT + 1
               NREC   = J
               READ(LUAORC(0),ERR=2000,REC=NREC) RBUF,IBUF4,LENGTH4
               LENGTH = LENGTH4
               CALL AOLAB4_cc(IBUF4,NIBUF,NBITS,INDX4,LENGTH)
               DO 110 I = 1,LENGTH
                  IP = INDX4(4,I)
                  IQ = INDX4(3,I)
                  IR = INDX4(2,I)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IP,IQ,IR,IDUM)
                  END IF
                  IADR = KADR1(IR) + KADR2(IP,IQ) + 1
                  XINT(IADR) = RBUF(I)
              
  110          CONTINUE
               IF (U21INT) THEN
               READ(LU21INT,ERR=2000,REC=NREC) RBUF,IBUF4,LENGTH4
               LENGTH = LENGTH4
               CALL AOLAB4_cc(IBUF4,NIBUF,NBITS,INDX4,LENGTH)
               DO 115 I = 1,LENGTH
                  IP = INDX4(4,I)
                  IQ = INDX4(3,I)
                  IR = INDX4(2,I)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IP,IQ,IR,IDUM)
                  END IF 
                  IADR = KADR1(IR) + KADR2(IP,IQ) + 1
                  XINT(IADR + IOFFU21) = RBUF(I)
  115          CONTINUE
               ENDIF
            ENDIF
C
  100    CONTINUE
C
C
      CALL GPCLOSE(LUAORC(0),'KEEP')
      IF (U21INT) CALL GPCLOSE(LU21INT,'KEEP')
C
      RETURN
 2000 CALL QUIT('Error in CCRDA1 reading Integral distribution')
      END
C  /* Deck lammat */
      SUBROUTINE LAMMAT(XLAMDP,XLAMDH,T1AM,WORK,LWORK)
C
C     Written by Henrik Koch 19-oct-1990.
C
C     PURPOSE:
C             Calculate transformation matrices used in ccsd
C             calculations.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      DIMENSION XLAMDH(*),XLAMDP(*),WORK(LWORK),T1AM(*)
#include "inftap.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "r12int.h"
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      CALL QENTER('LAMMAT')
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KCMO  = 1
      KEND  = KCMO + NLAMDS
      LWRK1 = LWORK - KEND
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient spaces in LAMMAT')
      ENDIF
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
C     LABEL is used instead of 'TRCCINT' (WK/UniKA/04-11-2002).
      CALL MOLLAB(LABEL,LUSIFC,LUPRI)
      READ (LUSIFC)
C
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
CHF
c      WRITE(LUPRI,*)'CMO out in lammat'
c      CALL OUTPUT(WORK(KCMO),1,NLAMDS,1,NLAMDS,NLAMDS,NLAMDS,1,LUPRI)
CHF
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------------------
C     Reorder the MO-coefficient matrix.
C---------------------------------------
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND),LWRK1)
C
C-------------------------------------------
C     Calculate the transformation matrices.
C-------------------------------------------
C
      CALL DCOPY(NLAMDT,WORK(KCMO),1,XLAMDH,1)
      CALL DCOPY(NLAMDT,WORK(KCMO),1,XLAMDP,1)
C
      CALL LAMDA1(XLAMDP,XLAMDH,T1AM,WORK(KCMO))
C
      IF (IPRINT .GT. 200 .OR. LOCDBG) THEN
C
         CALL AROUND('Lambda Particle matrix in LAMMAT')
         KOFF1 = 1
         KOFF2 = NLMRHF + 1
         DO 200 ISYM = 1,NSYM
            WRITE(LUPRI,1) ISYM
            WRITE(LUPRI,2)
            WRITE(LUPRI,3)
            IF (NRHF(ISYM) .EQ. 0) THEN
               WRITE(LUPRI,4)
               GOTO 210
            ENDIF
            CALL OUTPUT(XLAMDP(KOFF1),1,NBAS(ISYM),1,NRHF(ISYM),
     *                  NBAS(ISYM),NRHF(ISYM),1,LUPRI)
  210       WRITE(LUPRI,5)
            WRITE(LUPRI,6)
            IF (NVIR(ISYM) .EQ. 0) THEN
               WRITE(LUPRI,4)
               GOTO 220
            ENDIF
            CALL OUTPUT(XLAMDP(KOFF2),1,NBAS(ISYM),1,NVIR(ISYM),
     *                  NBAS(ISYM),NVIR(ISYM),1,LUPRI)
C
  220       CONTINUE
            KOFF1 = KOFF1 + NBAS(ISYM)*NRHF(ISYM)
            KOFF2 = KOFF2 + NBAS(ISYM)*NVIR(ISYM)
  200    CONTINUE
C
         CALL AROUND('Lambda Hole matrix in LAMMAT')
         KOFF1 = 1
         KOFF2 = NLMRHF + 1
         DO 300 ISYM = 1,NSYM
            WRITE(LUPRI,1) ISYM
            WRITE(LUPRI,7)
            WRITE(LUPRI,8)
            IF (NRHF(ISYM) .EQ. 0) THEN
               WRITE(LUPRI,4)
               GOTO 310
            ENDIF
            CALL OUTPUT(XLAMDH(KOFF1),1,NBAS(ISYM),1,NRHF(ISYM),
     *                  NBAS(ISYM),NRHF(ISYM),1,LUPRI)
  310       WRITE(LUPRI,9)
            WRITE(LUPRI,10)
            IF (NVIR(ISYM) .EQ. 0) THEN
               WRITE(LUPRI,4)
               GOTO 320
            ENDIF
            CALL OUTPUT(XLAMDH(KOFF2),1,NBAS(ISYM),1,NVIR(ISYM),
     *                  NBAS(ISYM),NVIR(ISYM),1,LUPRI)
C
  320       CONTINUE
            KOFF1 = KOFF1 + NBAS(ISYM)*NRHF(ISYM)
            KOFF2 = KOFF2 + NBAS(ISYM)*NVIR(ISYM)
  300    CONTINUE
C
      END IF
C
      CALL QEXIT('LAMMAT')
      RETURN
C
    1 FORMAT(/,/,7X,'Symmetry number :',I5)
    2 FORMAT(/,/,7X,'Lambda particle occupied part')
    3 FORMAT(7X,'-----------------------------')
    4 FORMAT(/,/,7X,'This symmetry is empty')
    5 FORMAT(/,/,7X,'Lambda particle virtual part')
    6 FORMAT(7X,'----------------------------')
    7 FORMAT(/,/,7X,'Lambda hole occupied part')
    8 FORMAT(7X,'-------------------------')
    9 FORMAT(/,/,7X,'Lambda hole virtual part')
   10 FORMAT(7X,'------------------------')
C
      END
C  /* Deck lamda1 */
      SUBROUTINE LAMDA1(XLAMDP,XLAMDH,T1AM,CMO)
C
C     Calculate the lambda matrices.             asm 05-08-94
C
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0)
      DIMENSION XLAMDH(*),XLAMDP(*)
      DIMENSION T1AM(*),CMO(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      DO 100 ISYMP = 1,NSYM
C
         ISYMI = ISYMP
         ISYMB = ISYMI
         ISYMA = ISYMP
         ISYMJ = ISYMA
C
         NBASP = MAX(NBAS(ISYMP),1)
         NVIRB = MAX(NVIR(ISYMB),1)
         NVIRA = MAX(NVIR(ISYMA),1)
C
         KOFF1 = ILMVIR(ISYMB) + 1
         KOFF2 = IT1AM(ISYMB,ISYMI) + 1
         KOFF3 = ILMRHF(ISYMI) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMP),NRHF(ISYMI),NVIR(ISYMB),
     *              ONE,CMO(KOFF1),NBASP,T1AM(KOFF2),NVIRB,
     *              ONE,XLAMDH(KOFF3),NBASP)
C
         KOFF4 = ILMRHF(ISYMJ) + 1
         KOFF5 = IT1AM(ISYMA,ISYMJ) + 1
         KOFF6 = ILMVIR(ISYMJ) + 1
C
         CALL DGEMM('N','T',NBAS(ISYMP),NVIR(ISYMA),NRHF(ISYMJ),
     *              -ONE,CMO(KOFF4),NBASP,T1AM(KOFF5),NVIRA,
     *              ONE,XLAMDP(KOFF6),NBASP)
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck sqmatr */
      SUBROUTINE SQMATR(NDIM,PKMAT,SQMAT)
C
C     Written by Henrik Koch 19-oct-1990.
C
C     PURPOSE:
C             Square up packed matrix.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PKMAT(*),SQMAT(NDIM,NDIM)
C
      DO 100 I = 1,NDIM
         DO 110 J = 1,I
C
            IJ = I*(I-1)/2 + J
            SQMAT(I,J) = PKMAT(IJ)
            SQMAT(J,I) = PKMAT(IJ)
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_t2ao */
      SUBROUTINE CC_T2AO(T2AM,XLAMDH,ISYMLH,SCRM,WORK,LWORK,
     *                   IDEL,ISYMD,ISYMTR,IOPT)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Written by Henrik Koch 22-dec-1993.
C     Symmetry due to Alfredo Sanchez and Henrik Koch 11-July 1994
C     Nontotal symmetric amplitudes Ove Christiansen 14-2-1995.
C     LAMDH is still assumed tot. sym.
C     Asger Halkier 13/2-1996: Generalised to handle "non-direct"
C     AO-index gamma in lampda matrix (IOPT = 2), as well as the
C     usual "direct" delta AO-index (IOPT = 1).
C     Ove Christiansen 16-6-1996:
C     Generalised to non-total symmetric Lamdba matrices.
C     PURPOSE:
C             Tdjci -> Tci,j (delta)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION T2AM(*),XLAMDH(*)
      DIMENSION SCRM(*),WORK(LWORK)
#include "ccorb.h"
#include "ccsdsym.h"
C
C-----------------------------------------------------
C     Calculate the transformed t2-amplitude and save.
C-----------------------------------------------------
C
      ISYDVI = MULD2H(ISYMLH,ISYMD)
      ISYMM = MULD2H(ISYMTR,ISYDVI)
      CALL DZERO(SCRM,NT2BCD(ISYMM))
C
      IF ( LWORK .LT. NVIR(ISYDVI)) THEN
         CALL QUIT('Insufficient core in CC_T2AO')
      ENDIF
C
      CALL DZERO(WORK,NVIR(ISYDVI))
C
      IF (IOPT .EQ. 1) THEN
         KOFF1 = IGLMVI(ISYMD,ISYDVI) + IDEL - IBAS(ISYMD)
      ELSE IF (IOPT .EQ. 2) THEN
         KOFF1 = IGLMVI(ISYMD,ISYDVI) + IDEL
      ENDIF
      CALL DCOPY(NVIR(ISYDVI),XLAMDH(KOFF1),NBAS(ISYMD),WORK,1)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMDJ = MULD2H(ISYMJ,ISYDVI)
         ISYMCI = MULD2H(ISYMTR,ISYMDJ)
C
         NTOTCI = MAX(NT1AM(ISYMCI),1)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            KDJ   = IT1AM(ISYDVI,ISYMJ) + NVIR(ISYDVI)*(J-1) + 1
            KOFF2 = IT2SQ(ISYMCI,ISYMDJ)
     *            + NT1AM(ISYMCI)*(KDJ - 1) + 1
            KOFF3 = IT2BCD(ISYMCI,ISYMJ)
     *            + NT1AM(ISYMCI)*(J-1) + 1
C
            CALL DGEMV('N',NT1AM(ISYMCI),NVIR(ISYDVI),ONE,
     *                 T2AM(KOFF2),NTOTCI,WORK,1,ZERO,
     *                 SCRM(KOFF3),1)
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck trsrec */
      SUBROUTINE TRSREC(NDIM1,NDIM2,XMAT1,XMAT2)
C
C     Written by Henrik Koch 19-oct-1990.
C
C     PURPOSE:
C            Transpose rectangular matrix.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XMAT1(NDIM1,NDIM2),XMAT2(NDIM2,NDIM1)
C
      DO 100 I = 1,NDIM1
         DO 110 J = 1,NDIM2
C
            XMAT2(J,I) = XMAT1(I,J)
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccrhs_oneao */
      SUBROUTINE CCRHS_ONEAO(FOCK,WORK,LWRK)
C
C     Written by Henrik Koch & Ove Christiansen 24-jan-1994.
C     Symmetry due to Alfredo Sanchez and Henrik Koch 11-July 1994
C
C     PURPOSE:
C             Read one electron integrals into matrix.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      DIMENSION FOCK(*),WORK(*)
#include "inftap.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      LOGICAL EX
C
      IF (LWRK .LT. NNBST(ISYMOP))
     *           CALL QUIT('Insufficient space in CCRHS_ONEAO')
C
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK,NNBST(ISYMOP))
      CALL CCSD_SYMSQ(WORK,ISYMOP,FOCK)
C
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('One electron AO-integrals in fock matrix')
         KOFF1 = 1
         DO 110 ISYMB = 1,NSYM
            WRITE(LUPRI,*) 'Symmetry number : ',ISYMB
            NBASB = NBAS(ISYMB)
            CALL OUTPUT(FOCK(KOFF1),1,NBASB,1,NBASB,NBASB,NBASB,1,LUPRI)
            KOFF1 = KOFF1 + NBAS(ISYMB)*NBAS(ISYMB)
  110    CONTINUE
C
      ENDIF
      RETURN
      END
C  /* Deck cc_t2sq */
      SUBROUTINE CC_T2SQ(T2AM,T2SQ,ISYM)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Henrik Koch and Alfredo Sanchez.       11-July-1994
C     Modified by Ove Christiansen 24-1-1995 to handle
C     a general non total symmetric vector.
C     Squareup the t2-amplitudes distribution.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      DIMENSION T2AM(*),T2SQ(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      IF (ISYM.EQ.1) THEN
         KOFF1 = 1
         KOFF2 = 1
         DO 100 ISYMBJ = 1,NSYM
            CALL SQMATR(NT1AM(ISYMBJ),T2AM(KOFF1),T2SQ(KOFF2))
            KOFF1 = KOFF1 + NT1AM(ISYMBJ)*(NT1AM(ISYMBJ)+1)/2
            KOFF2 = KOFF2 + NT1AM(ISYMBJ)*NT1AM(ISYMBJ)
  100    CONTINUE
C
      ELSE
C
         KOFF = 1
         DO 200 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYM,ISYMBJ)
C
            IF (ISYMBJ.GT.ISYMAI) THEN
C
               NAMP = NT1AM(ISYMAI)*NT1AM(ISYMBJ)
               KOFF1 = IT2SQ(ISYMAI,ISYMBJ) + 1
               CALL DCOPY(NAMP,T2AM(KOFF),1,T2SQ(KOFF1),1)
               NAI = MAX(NT1AM(ISYMAI),1)
               NBJ = MAX(NT1AM(ISYMBJ),1)
               KOFF2 = IT2SQ(ISYMBJ,ISYMAI) + 1
               CALL TRM(T2AM(KOFF),NAI,NT1AM(ISYMAI),NT1AM(ISYMBJ),
     *                     T2SQ(KOFF2),NBJ)
               KOFF = KOFF + NAMP
C
            ENDIF
C
  200    CONTINUE
C
      ENDIF
C
      RETURN
      END
C  /* Deck trm */
      SUBROUTINE TRM(A,LDA,M,N,B,LDB)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Transpose matrix A dim m,n in array with logical dim. lda.
C     and put result into B with logical dim. ldb.
C     Use dcopy for vectorization.
C
C     Ove Christiansen 14-2-1995
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
C
      DIMENSION A(LDA,*),B(LDB,*)
C
      DO 100 I = 1, N
C
         CALL DCOPY(M,A(1,I),1,B(I,1),LDB)
C
 100  CONTINUE
C
      RETURN
      END
C  /* Deck cc_aodens */
      SUBROUTINE CC_AODENS(XLAMDP,XLAMDH,DENS,ISYMH,IC,WORK,LWORK)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Henrik Koch and Alfredo Sanchez.       11-July-1994
C
C     Calculate the AO-density matrix used in constructing
C     the AO Fock matrix.
C
C
C     Ove Christiansen 13-7-1995
C         generalise to non-totalsymmetric lambda matrices
C         for C1 transformation.
C         ISYMH is the symmetry of the transformed LAMBDAH
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XLAMDP(*), XLAMDH(*), DENS(*), WORK(LWORK)
#include "inftap.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "r12int.h"
C
      KOFF1 = 1
      KOFF2 = 1
      KOFF3 = 1
C
      DO 100 ISYMB = 1,NSYM
C
         ISYMA = MULD2H(ISYMH,ISYMB)
         ISYMK = ISYMA
         NBASA = MAX(NBAS(ISYMA),1)
         NBASB = MAX(NBAS(ISYMB),1)
C
         KOFF1 = 1 + IGLMRH(ISYMA,ISYMK)
         KOFF2 = 1 + IGLMRH(ISYMB,ISYMK)
C
         CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYMK),ONE,
     *              XLAMDP(KOFF1),NBASA,XLAMDH(KOFF2),NBASB,ZERO,
     *              DENS(KOFF3),NBASA)
C
         KOFF3 = KOFF3 + NBAS(ISYMA)*NBAS(ISYMB)
C
  100 CONTINUE
C
C
C-----------------------------
C     Include frozen orbitals.
C-----------------------------
C
      IF ((FROIMP .OR. FROEXP).AND.(IC .EQ. 1)) THEN
C
         IF (LWORK .LT. NLAMDS) THEN
            CALL QUIT('Insufficient space in CCSD_AODENS')
         ENDIF
C
C-------------------------------------------------
C        Read MO-coefficients from interface file.
C-------------------------------------------------
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
C
C        Use LABEL instead of 'TRCCINT ' (WK/UniKA/04-11-2002).
         CALL MOLLAB(LABEL,LUSIFC,LUPRI)
         READ (LUSIFC)
C
         READ (LUSIFC)
         READ (LUSIFC) (WORK(I), I=1,NLAMDS)
C
         CALL GPCLOSE(LUSIFC,'KEEP')
C
C-------------------------------------------------------
C        Add contribution from frozen occupied orbitals.
C-------------------------------------------------------
C
         KOFF1 = 0
         KOFF2 = 0
         DO 200 ISYMK = 1,NSYM
C
            ISYMA = ISYMK
            ISYMB = ISYMK
C
            DO 210 II = 1,NRHFFR(ISYMK)
C
               K = KFRRHF(II,ISYMK)
C
               DO 220 B = 1,NBAS(ISYMB)
                  DO 230 A = 1,NBAS(ISYMA)
C
                     NAK = KOFF1 + NBAS(ISYMA)*(K - 1) + A
                     NBK = KOFF1 + NBAS(ISYMB)*(K - 1) + B
                     NAB = KOFF2 + NBAS(ISYMA)*(B - 1) + A
C
                     DENS(NAB) = DENS(NAB) + WORK(NAK)*WORK(NBK)
C
  230             CONTINUE
  220          CONTINUE
C
  210       CONTINUE
C
            KOFF1 = KOFF1 + NBAS(ISYMK)*NORBS(ISYMK)
            KOFF2 = KOFF2 + NBAS(ISYMA)*NBAS(ISYMB)
C
  200    CONTINUE
C
      ENDIF
C
      END
C  /* Deck cc_t2mo */
      SUBROUTINE CC_T2MO(RHO1,CTR2,ISYMC2,OMEGA2,RHO2,GAMMA,XLAMDP,
     *                   XLAMPC,ISYMPC,WORK,LWORK,ISYMBF,ICON)
C
C     Henrik Koch and Alfredo Sanchez.       15-July-1994
C
C     Transform the Omega2 vector from the AO basis to the MO
C     basis.
C
C     Ove Christiansen 4-8-1995:
C
C     Generalizations for CC response.
C
C        1.ISYMBF is the symmetry of the BF (ali,bej) vector.
C        2.Transform with a non total symmetric lambda matrix.
C          (one with sym 1 and one with sym isympc)
C
C        note that if newgam is true gamma is the gamma vector on return
C        with the same symmetry as the input BF. (transformed with xlamdp)
C
C        if newgam is false the gamma intermediate is not returned.
C
C        ICON is 2 for response to calculat a-tild,ibj and ai,b-tilde,j
C
C        NB these changes are only carried through completely and
C        tested for omegor
C
C     Asger Halkier 2/11-1995:
C
C        For ICON equal to 3 the contraction of the (ali,bej) vector with
C        the trialvector CTR2 (i.e the LT21BF-term) is calculated and
C        stored in RHO1!
C
C     Ove Christiansen 4-10-1996:
C 
C        For use in F-matrix generalize ICON .EQ. 3 section
C
C     NOTE: Linear response options only valid and debugged for OMEGOR!
C
C     Christian Neiss 09/11/2005:
C        ICON .EQ. 4: transform only beta index to occupied using XLAMDP
C        (--> only total-symmetric transf. allowed); result is
C        added on GAMMA; RHO2 will not be used
C        Dimension of GAMMA = NT2AOIJ(ISYMO2)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION RHO1(*), CTR2(*), OMEGA2(*), RHO2(*), GAMMA(*),
     *          XLAMDP(*), WORK(*), XLAMPC(*)
#include "ccorb.h"
#include "ccsdsym.h"
#include "symsq.h"
#include "cclr.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYMO2 = MULD2H(ISYMBF,ISYMPC)
      ISYMO1 = MULD2H(ISYMO2,ISYMC2)
C
      IF ((ICON.EQ.1).OR.(ICON.EQ.2)) THEN
         CALL DZERO(RHO2,NT2AM(ISYMO2))
      ENDIF
C
      DO 100 ISYMJ = 1,NSYM
         DO 110 ISYMI = 1,NSYM
C
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
            ISALBE = MULD2H(ISYMIJ,ISYMBF)
            ISYMAB = MULD2H(ISYMIJ,ISYMO2)
C
            DO 120 ISYBE = 1,NSYM
C
               ISYAL  = MULD2H(ISYBE,ISALBE)
               ISYALI = MULD2H(ISYAL,ISYMI)
               ISYBEJ = MULD2H(ISYBE,ISYMJ)
C
C-----------------------------------------------
C              Dynamic allocation of work space.
C-----------------------------------------------
C
               ISYMA = MULD2H(ISYAL,ISYMPC)
               NVA = MAX(NVIR(ISYMA),NVIR(ISYAL))
               NRA = MAX(NRHF(ISYMA),NRHF(ISYAL))
               ISYMB = MULD2H(ISYBE,ISYMPC)
               NVB = MAX(NVIR(ISYMB),NVIR(ISYBE),NRHF(ISYBE))
               NRB = MAX(NRHF(ISYMB),NRHF(ISYBE))
C
               KSCR1 = 1
               IF (ICON.NE.4) THEN
                  KSCR2 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)
                  KSCR3 = KSCR2 + NBAS(ISYAL)*NVB
                  IF (NEWGAM) THEN
                     KSCR4 = KSCR3 + NVA*NVB
                     KSCR5 = KSCR4 + NBAS(ISYAL)*NRB
                     KEND1 = KSCR5 + NRA*NRB
                  ELSE
                     KEND1 = KSCR3 + NVA*NVB
                  END IF
               ELSE
                  KSCR4 = KSCR1 + NBAS(ISYAL)*NBAS(ISYBE)
                  KEND1 = KSCR4 + NBAS(ISYAL)*NRB
               END IF
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Not enough space in CC_T2MO')
               END IF
C
               DO 130 J = 1,NRHF(ISYMJ)
                  DO 140 I = 1,NRHF(ISYMI)
C
C------------------------------------------
C                    Squareup the AB block.
C------------------------------------------
C
                     IF ((.NOT. OMEGSQ) .AND. (.NOT. OMEGOR)) THEN
C
                     DO 150 B = 1,NBAS(ISYBE)
                        NBJ   = IT1AO(ISYBE,ISYMJ)
     *                        + NBAS(ISYBE)*(J-1) + B
                        DO 155 A = 1,NBAS(ISYAL)
C
                           NAI   = IT1AO(ISYAL,ISYMI)
     *                           + NBAS(ISYAL)*(I-1) + A
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           IF (ISYMO2 .EQ. 1) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYALI .LT. ISYBEJ) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + NT1AO(ISYALI)*(NBJ - 1) + NAI
                           ELSEIF (ISYALI .GT. ISYBEJ) THEN
                              NAIBJ = IT2AO(ISYALI,ISYBEJ)
     *                              + NT1AO(ISYBEJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           WORK(NAB) = OMEGA2(NAIBJ)
C
  155                   CONTINUE
  150                CONTINUE
C
                     ENDIF
C
                     IF (OMEGSQ) THEN
C
                     DO 160 B = 1,NBAS(ISYBE)
                        NBJ   = IT1AO(ISYBE,ISYMJ)
     *                        + NBAS(ISYBE)*(J-1) + B
                        DO 165 A = 1,NBAS(ISYAL)
C
                           NAI   = IT1AO(ISYAL,ISYMI)
     *                           + NBAS(ISYAL)*(I-1) + A
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           NAIBJ = IT2AOS(ISYALI,ISYBEJ)
     *                           + NT1AO(ISYALI)*(NBJ - 1) + NAI
                           NBJAI = IT2AOS(ISYBEJ,ISYALI)
     *                           + NT1AO(ISYBEJ)*(NAI - 1) + NBJ
C
                           WORK(NAB) = OMEGA2(NAIBJ) + OMEGA2(NBJAI)
C
  165                   CONTINUE
  160                CONTINUE
C
                     ENDIF
C
                     IF (OMEGOR) THEN
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                        FAC1 = ONE
                        IF (I .GT. J) FAC1 = -ONE
                     ELSE IF (ISYMI .LT. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMI)*(J - 1) + I
                        FAC1 = ONE
                     ELSE
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMJ)*(I - 1) + J
                        FAC1 = -ONE
                     ENDIF
C
                     DO 166 B = 1,NBAS(ISYBE)
                        DO 167 A = 1,NBAS(ISYAL)
C
                           IF (ISYAL .EQ. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + INDEX(A,B)
                              FAC2 = ONE
                              IF (A .GT. B) FAC2 = -ONE
                           ELSE IF (ISYAL .LT. ISYBE) THEN
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYAL)*(B - 1) + A
                              FAC2 = ONE
                           ELSE
                              NABP = IAODPK(ISYAL,ISYBE)
     *                             + NBAS(ISYBE)*(A - 1) + B
                              FAC2 = -ONE
                           ENDIF
C
                           NABIJP = IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
C
                           NABIJM = NT2ORT(ISYMBF)
     *                            + IT2ORT(ISALBE,ISYMIJ)
     *                            + NNBST(ISALBE)*(NIJ - 1) + NABP
C
                           NAB   = KSCR1 + NBAS(ISYAL)*(B - 1) + A - 1
C
                           FAC = FAC1*FAC2
C
                           WORK(NAB) =
     *                       HALF*(OMEGA2(NABIJP) + FAC*OMEGA2(NABIJM))
C
  167                   CONTINUE
  166                CONTINUE
C
                     ENDIF
C
C------------------------------------------------------------
C                    Transform the AB block to virtual space.
C------------------------------------------------------------
C
                     IF ((ICON.EQ.1).OR.(ICON.EQ.2)) THEN
C
                     ISYMA = MULD2H(ISYAL,ISYMPC)
                     ISYMB = ISYBE
                     ISYMAI = MULD2H(ISYMA,ISYMI)
                     ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NVIRA = MAX(NVIR(ISYMA),1)
C
                     KOFF1 = ILMVIR(ISYBE) + 1
C
                     CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                          NBASA)
C
                     KOFF2 = IGLMVI(ISYAL,ISYMA) + 1
C
                     CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                          NBAS(ISYAL),ONE,XLAMPC(KOFF2),NBASA,
     *                          WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                          NVIRA)
C
C--------------------------------------------
C                    Store the omega2 vector.
C--------------------------------------------
C
                     DO 170 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 180 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
                           NAB   = KSCR3 + NVIR(ISYMA)*(B - 1) + A - 1
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
C
                              IF (NAI .GT. NBJ) GOTO 180
C
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              GOTO 180
chjaaj: next two lines are commented because it is dead code
c                             NAIBJ = IT2AM(ISYMAI,ISYMBJ)
c    *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
  180                   CONTINUE
  170                CONTINUE
C
                     ENDIF
C
C--------------------------------------
C                    CCLR contribution.
C--------------------------------------
C
                     IF (ICON .EQ. 2 ) THEN
C
                        CALL DZERO(WORK(KSCR2),NVA*NVB)
                        ISYMA = ISYAL
                        ISYMB = MULD2H(ISYBE,ISYMPC)
                        ISYMAI = MULD2H(ISYMA,ISYMI)
                        ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
                        NBASA = MAX(NBAS(ISYAL),1)
                        NBASB = MAX(NBAS(ISYBE),1)
                        NVIRA = MAX(NVIR(ISYMA),1)
C
                        KOFF1 = IGLMVI(ISYBE,ISYMB) + 1
C
                        CALL DGEMM('N','N',NBAS(ISYAL),NVIR(ISYMB),
     *                             NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                             XLAMPC(KOFF1),NBASB,ZERO,WORK(KSCR2),
     *                             NBASA)
C
                        KOFF2 = ILMVIR(ISYAL) + 1
C
                        CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                             NBAS(ISYAL),ONE,XLAMDP(KOFF2),NBASA,
     *                             WORK(KSCR2),NBASA,ZERO,WORK(KSCR3),
     *                             NVIRA)
C
C--------------------------------------------
C                    Store the omega2 vector.
C--------------------------------------------
C
                     DO 181 B = 1,NVIR(ISYMB)
                        NBJ   = IT1AM(ISYMB,ISYMJ)
     *                        + NVIR(ISYMB)*(J-1) + B
                        DO 182 A = 1,NVIR(ISYMA)
C
                           NAI   = IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1) + A
C
                           IF (ISYMAI .EQ. ISYMBJ) THEN
                              IF (NAI .GT. NBJ ) GOTO 182
                              NAIBJ = IT2AM(ISYALI,ISYBEJ)
     *                              + INDEX(NAI,NBJ)
                           ELSEIF (ISYMAI .LT. ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                              + NT1AM(ISYMAI)*(NBJ - 1) + NAI
                           ELSEIF (ISYMAI .GT. ISYMBJ) THEN
                              GOTO 182
chjaaj: next two lines are commented because it is dead code
c                             NAIBJ = IT2AM(ISYMAI,ISYMBJ)
c    *                              + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
                           ENDIF
C
                           NAB  = KSCR3+ NVIR(ISYMA)*(B - 1) + A - 1
                           RHO2(NAIBJ) = RHO2(NAIBJ) + WORK(NAB)
C
  182                   CONTINUE
  181                CONTINUE
C
                     ENDIF
C
C============================================================
C                    Section for calculating the LT21BF-term.
C============================================================
C
                     IF (ICON .EQ. 3) THEN
C
                        ISYMK = ISYBE
                        ISYMD = MULD2H(ISYAL,ISYMPC)
                        ISYMC = MULD2H(ISYMK,ISYMO1)
                        ISYDI = MULD2H(ISYMD,ISYMI)
                        ISYCJ = MULD2H(ISYMC,ISYMJ)
C
                        LENGTH = NBAS(ISYAL)*NRHF(ISYMK)
C
                        CALL DZERO(WORK(KSCR2),LENGTH)
C
C----------------------------------------------------------
C                       Transform the AO-block to MO-basis.
C----------------------------------------------------------
C
                        KOFF1  = ILMRHF(ISYMK) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTBE = MAX(NBAS(ISYBE),1)
C
                        CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYMK),
     *                             NBAS(ISYBE),ONE,WORK(KSCR1),NTOTAL,
     *                             XLAMDP(KOFF1),NTOTBE,ZERO,
     *                             WORK(KSCR2),NTOTAL)
C
                        KOFF2  = IGLMVI(ISYAL,ISYMD) + 1
C
                        NTOTAL = MAX(NBAS(ISYAL),1)
                        NTOTK  = MAX(NRHF(ISYMK),1)
C
                        CALL DGEMM('T','N',NRHF(ISYMK),NVIR(ISYMD),
     *                             NBAS(ISYAL),ONE,WORK(KSCR2),NTOTAL,
     *                             XLAMPC(KOFF2),NTOTAL,ZERO,
     *                             WORK(KSCR3),NTOTK)
C
C-----------------------------------------------------------------
C                       Contraction with CTR2 & storage in result.
C-----------------------------------------------------------------
C
                        DO 47 C = 1,NVIR(ISYMC)
C
                           NCJ   = IT1AM(ISYMC,ISYMJ)
     *                           + NVIR(ISYMC)*(J - 1) + C
                           NDICJ = IT2SQ(ISYDI,ISYCJ)
     *                           + NT1AM(ISYDI)*(NCJ - 1)
     *                           + IT1AM(ISYMD,ISYMI)
     *                           + NVIR(ISYMD)*(I - 1) + 1
                           NCK   = IT1AM(ISYMC,ISYMK) + C
C
                           CALL DGEMV('N',NRHF(ISYMK),NVIR(ISYMD),
     *                                -ONE,WORK(KSCR3),NTOTK,
     *                                CTR2(NDICJ),1,ONE,RHO1(NCK),
     *                                NVIR(ISYMC))
C
  47                    CONTINUE
C
                     ENDIF
C
C-------------------------------------------------------------
C                    Transform the AB block to occupied space.
C-------------------------------------------------------------
C
                     IF (.NOT.(NEWGAM.OR.(ICON.EQ.4))) GOTO 999
C
                     NBASA = MAX(NBAS(ISYAL),1)
                     NBASB = MAX(NBAS(ISYBE),1)
                     NRHFA1 = MAX(NRHF(ISYAL),1)
C
                     KOFF1 = ILMRHF(ISYBE) + 1
C
                     CALL DGEMM('N','N',NBAS(ISYAL),NRHF(ISYBE),
     *                          NBAS(ISYBE),ONE,WORK(KSCR1),NBASA,
     *                          XLAMDP(KOFF1),NBASB,ZERO,WORK(KSCR4),
     *                          NBASA)
C
                     IF (ICON.NE.4) THEN
C
                     KOFF2 = ILMRHF(ISYAL) + 1
C
                     CALL DGEMM('T','N',NRHF(ISYAL),NRHF(ISYBE),
     *                          NBAS(ISYAL),ONE,XLAMDP(KOFF2),NBASA,
     *                          WORK(KSCR4),NBASA,ZERO,WORK(KSCR5),
     *                          NRHFA1)
C
C-------------------------------------------
C                    Store the gamma matrix.
C-------------------------------------------
C
                     ISYMK = ISYAL
                     ISYML = ISYBE
C
                     ISYMKI = MULD2H(ISYMK,ISYMI)
                     ISYMLJ = MULD2H(ISYML,ISYMJ)
C
                     DO 190 L = 1,NRHF(ISYML)
C
                        NLJ = IMATIJ(ISYML,ISYMJ)
     *                      + NRHF(ISYML)*(J - 1) + L
C
                        DO 200 K = 1,NRHF(ISYMK)
C
                           NKL = KSCR5 + NRHF(ISYMK)*(L - 1) + K - 1
C
                           NKI = IMATIJ(ISYMK,ISYMI)
     *                         + NRHF(ISYMK)*(I - 1) + K
C
                           IF (ISYMKI .EQ. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + INDEX(NKI,NLJ)
                              GAMMA(NKILJ) = WORK(NKL)
                           ELSE IF (ISYMKI .LT. ISYMLJ) THEN
                              NKILJ = IGAMMA(ISYMKI,ISYMLJ)
     *                              + NMATIJ(ISYMKI)*(NLJ - 1) + NKI
                              GAMMA(NKILJ) = WORK(NKL)
                           ENDIF
C
  200                   CONTINUE
  190                CONTINUE
C
                     ELSE
C
C------------------------------------------------------------------
C                    Store "half-transformed" GAMMA for ICON .EQ. 4
C------------------------------------------------------------------
C
                     ISYML = ISYBE
                     ISYLAL = MULD2H(ISYAL,ISYML)
C
                     NIJ = IMATIJ(ISYMI,ISYMJ) + 
     *                     NRHF(ISYMI)*(J-1) + I
C
                     NALIJ = IT2AOIJ(ISYLAL,ISYMIJ) +
     *                       NT1AO(ISYLAL)*(NIJ-1) + 
     *                       IT1AO(ISYAL,ISYML) + 1
C
                     CALL DAXPY(NBAS(ISYAL)*NRHF(ISYML),ONE,
     *                          WORK(KSCR4),1,GAMMA(NALIJ),1)
C
                     END IF
C
  999                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsd_t2mtp */
      SUBROUTINE CCSD_T2MTP(SCRM,WORK,LWORK,ISYMD)
C
C     Alfredo Sanchez and Henrik Koch 26-July 1994
C
C     PURPOSE:
C             Transpose ij index of the T2M-amplitudes.
C
#include "implicit.h"
      DIMENSION SCRM(*)
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C-------------------------------------------
C     Calculate the transposed t2-amplitude.
C-------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMCI = MULD2H(ISYMJ,ISYMD)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMI = 1,ISYMJ
C
               ISYMC  = MULD2H(ISYMI,ISYMCI)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
C
               IF (LWORK .LT. NVIR(ISYMC)) THEN
                  CALL QUIT('Insufficient space in CCSD_T2MTP')
               ENDIF
C
               IF (ISYMI .EQ. ISYMJ) THEN
                  NRHFI = J - 1
               ELSE
                  NRHFI = NRHF(ISYMI)
               END IF
C
               DO 130 I = 1,NRHFI
C
                  NCIJ = IT2BCD(ISYMCI,ISYMJ) + NT1AM(ISYMCI)*(J-1)
     *                 + IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I-1) + 1
                  NCJI = IT2BCD(ISYMCJ,ISYMI) + NT1AM(ISYMCJ)*(I-1)
     *                 + IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J-1) + 1
C
                  CALL DCOPY(NVIR(ISYMC),SCRM(NCIJ),1,WORK,1)
                  CALL DCOPY(NVIR(ISYMC),SCRM(NCJI),1,SCRM(NCIJ),1)
                  CALL DCOPY(NVIR(ISYMC),WORK,1,SCRM(NCJI),1)
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsd_t2tp */
      SUBROUTINE CCSD_T2TP(T2AM,WORK,LWORK,ISYM)
C
C     Alfredo Sanchez and Henrik Koch 26-July 1994
C
C     PURPOSE:
C             Transpose ij index of the T2-amplitudes.
C             The amplitudes are assumed to be a square matrix.
C
#include "implicit.h"
      DIMENSION T2AM(*)
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C-------------------------------------------
C     Calculate the transposed t2-amplitude.
C-------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYM)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
C
                  DO 140 ISYMI = 1,ISYMJ
C
                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     IF (LWORK .LT. NVIR(ISYMA)) THEN
                        CALL QUIT('Insufficient space in CCSD_T2TP')
                     ENDIF
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NRHFI = J - 1
                     ELSE
                        NRHFI = NRHF(ISYMI)
                     END IF
C
                     DO 150 I = 1,NRHFI
C
                        NBI = IT1AM(ISYMB,ISYMI)+NVIR(ISYMB)*(I-1)+B
C
                        NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1)
     *                        + IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+1
C
                        NAJBI = IT2SQ(ISYMAJ,ISYMBI)
     *                        + NT1AM(ISYMAJ)*(NBI-1)
     *                        + IT1AM(ISYMA,ISYMJ)+NVIR(ISYMA)*(J-1)+1
C
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAIBJ),1,WORK,1)
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAJBI),1,
     *                             T2AM(NAIBJ),1)
                        CALL DCOPY(NVIR(ISYMA),WORK,1,T2AM(NAJBI),1)
C
  150               CONTINUE
C
  140             CONTINUE
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsd_invldp */
      SUBROUTINE CCSD_INVLDP(XLAMDP,XLAMIP,WORK,LWORK)
C
C     Alfredo Sanchez and Henrik Koch 26-July 1994
C
C     PURPOSE:
C             Invert the lambda particle matrix.
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
#include "iratdef.h"
      DIMENSION XLAMDP(*), XLAMIP(*)
      DIMENSION WORK(LWORK)
      DIMENSION DET(2)
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      DO 100 ISYMA = 1,NSYM
C
         KSCR  = 1
         KEND1 = KSCR + NBAS(ISYMA)*NORB(ISYMA)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for '//
     &           'allocation in CCSD_INVLDP')
         END IF
C
         NTOTR = NBAS(ISYMA)*NRHF(ISYMA)
C
         KOFF1 = ILMRHF(ISYMA) + 1
C
         CALL DCOPY(NTOTR,XLAMDP(KOFF1),1,WORK(KSCR),1)
C
         NTOTV = NBAS(ISYMA)*NVIR(ISYMA)
         KOFF2 = ILMVIR(ISYMA) + 1
         KOFF3 = KSCR + NTOTR
C
         CALL DCOPY(NTOTV,XLAMDP(KOFF2),1,WORK(KOFF3),1)
C
C
         NBASA = MAX(NBAS(ISYMA),1)
C
#if defined (SYS_xxx)
         NAUX = INT(32.5D0*DFLOAT(NBAS(ISYMA))) + 1
         IF (LWRK1. LT. NAUX) THEN
            CALL QUIT('Not enough space for DGEICD in CCSD_INVLDP')
         END IF

         CALL DGEICD(WORK(KSCR),NBASA,NBAS(ISYMA),0,RCOND,DET,
     *               WORK(KEND1),LWRK1)
#else
         NBASA2 = MAX(NBAS(ISYMA),1)
         NBASA1 = NBAS(ISYMA)
C
         KIPVT = KEND1
         KEND2 = KIPVT + NBAS(ISYMA)/IRAT + 1
         LWRK2 = LWORK - KEND2
         IF (LWRK2. LT. NBASA1) THEN
            CALL QUIT('Not enough space for DGEDI in CCSD_INVLDP')
         END IF
C
#if !defined (SYS_CRAY)
         IF (NBAS(ISYMA) .GT. 1) THEN
            CALL DGEFA(WORK(KSCR),NBAS(ISYMA),NBAS(ISYMA),
     *                 WORK(KIPVT),IERR)
         END IF
C
         CALL DGEDI(WORK(KSCR),NBASA2,NBASA1,WORK(KIPVT),DET,
     *              WORK(KEND2),1)
#else
         IF (NBAS(ISYMA) .GT. 1) THEN
            CALL SGEFA(WORK(KSCR),NBAS(ISYMA),NBAS(ISYMA),
     *                 WORK(KIPVT),IERR)
         END IF
C
         CALL SGEDI(WORK(KSCR),NBASA2,NBASA1,WORK(KIPVT),DET,
     *              WORK(KEND2),1)
#endif
#endif
C
         DO 110 I = 1,NRHF(ISYMA)
C
            KOFF1 = KSCR + I - 1
            KOFF2 = ILMRHF(ISYMA) + NBAS(ISYMA)*(I-1) + 1
C
            CALL DCOPY(NBAS(ISYMA),WORK(KOFF1),NBAS(ISYMA),
     *                 XLAMIP(KOFF2),1)
C
  110    CONTINUE
C
         DO 120 A = 1,NVIR(ISYMA)
C
            KOFF1 = KSCR + NRHF(ISYMA) + A - 1
            KOFF2 = ILMVIR(ISYMA) + NBAS(ISYMA)*(A-1) + 1
C
            CALL DCOPY(NBAS(ISYMA),WORK(KOFF1),NBAS(ISYMA),
     *                 XLAMIP(KOFF2),1)
C
  120    CONTINUE
C
100   CONTINUE
C
C------------------
C     Test section.
C------------------
C
      IF (IPRINT .GT. 120) THEN
C
      CALL AROUND('The inverse lambda matrix. Occupied part')
      DO 199 ISYMI = 1,NSYM
C
         WRITE(LUPRI,*) 'The symmetry of the block :',ISYMI
C
         KOFF1 = ILMRHF(ISYMI) + 1
C
         CALL OUTPUT(XLAMIP(KOFF1),1,NBAS(ISYMI),1,NRHF(ISYMI),
     *               NBAS(ISYMI),NRHF(ISYMI),1,LUPRI)
C
  199 CONTINUE
C
      CALL AROUND('Test of the occupied part of inverse xlamdp')
      DO 200 ISYMI = 1,NSYM
C
         WRITE(LUPRI,*) 'The symmetry of the block :',ISYMI
C
         NBASI = MAX(NBAS(ISYMI),1)
         NRHFI = MAX(NRHF(ISYMI),1)
C
         KOFF1 = ILMRHF(ISYMI) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMI),NBAS(ISYMI),ONE,
     *              XLAMDP(KOFF1),NBASI,XLAMIP(KOFF1),NBASI,ZERO,
     *              WORK,NRHFI)
C
         CALL OUTPUT(WORK,1,NRHF(ISYMI),1,NRHF(ISYMI),NRHF(ISYMI),
     *               NRHF(ISYMI),1,LUPRI)
C
  200 CONTINUE
C
      CALL AROUND('The inverse lambda matrix. Virtual part')
      DO 209 ISYMA = 1,NSYM
C
         WRITE(LUPRI,*) 'The symmetry of the block :',ISYMI
C
         KOFF1 = ILMVIR(ISYMA) + 1
C
         CALL OUTPUT(XLAMIP(KOFF1),1,NBAS(ISYMA),1,NVIR(ISYMA),
     *               NBAS(ISYMA),NVIR(ISYMA),1,LUPRI)
C
  209 CONTINUE
C
      CALL AROUND('Test of the virtual part of inverse xlamdp')
      DO 210 ISYMA = 1,NSYM
C
         WRITE(LUPRI,*) 'The symmetry of the block :',ISYMA
C
         NBASA = MAX(NBAS(ISYMA),1)
         NVIRA = MAX(NVIR(ISYMA),1)
C
         KOFF1 = ILMVIR(ISYMA) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMA),NBAS(ISYMA),ONE,
     *              XLAMDP(KOFF1),NBASA,XLAMIP(KOFF1),NBASA,ZERO,
     *              WORK,NVIRA)
C
         CALL OUTPUT(WORK,1,NVIR(ISYMA),1,NVIR(ISYMA),NVIR(ISYMA),
     *               NVIR(ISYMA),1,LUPRI)
C
  210 CONTINUE
C
      ENDIF
C
      RETURN
      END
C  /* Deck ccrhs_t2tr */
      SUBROUTINE CCRHS_T2TR(T2AM,WORK,LWORK,ISYM)
C
C     Alfredo Sanchez and Henrik Koch 30-July 1994
C
C     PURPOSE:
C             Calculate two coulomb minus exchange of t2 amplitudes.
C             The amplitudes are assumed to be a square matrix.
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION T2AM(*)
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
C----------------------------------------------------------
C     Calculate two coulomb minus exchange of t2-amplitude.
C----------------------------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYM)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
C
                  DO 140 ISYMI = 1,ISYMJ
C
                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     KSCR1 = 1
                     KSCR2 = KSCR1 + NVIR(ISYMA)
                     KEND1 = KSCR2 + NVIR(ISYMA)
                     LWRK1 = LWORK - KEND1
                     IF (LWRK1 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCRHS_T2TR')
                     ENDIF
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NRHFI = J - 1
                     ELSE
                        NRHFI = NRHF(ISYMI)
                     END IF
C
                     DO 150 I = 1,NRHFI
C
                        NBI = IT1AM(ISYMB,ISYMI)+NVIR(ISYMB)*(I-1)+B
C
                        NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1)
     *                        + IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+1
C
                        NAJBI = IT2SQ(ISYMAJ,ISYMBI)
     *                        + NT1AM(ISYMAJ)*(NBI-1)
     *                        + IT1AM(ISYMA,ISYMJ)+NVIR(ISYMA)*(J-1)+1
C
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAIBJ),1,
     *                             WORK(KSCR1),1)
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAJBI),1,
     *                             WORK(KSCR2),1)
C
                        CALL DSCAL(NVIR(ISYMA),TWO,T2AM(NAIBJ),1)
                        CALL DSCAL(NVIR(ISYMA),TWO,T2AM(NAJBI),1)
C
                        CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KSCR2),1,
     *                             T2AM(NAIBJ),1)
                        CALL DAXPY(NVIR(ISYMA),-ONE,WORK(KSCR1),1,
     *                             T2AM(NAJBI),1)
C
  150               CONTINUE
C
  140             CONTINUE
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('Two coulomb minus exchamge of t2am')
         DO 200 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYM)
            KOFF = IT2SQ(ISYMAI,ISYMBJ) + 1
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Symmetry block:',ISYMBJ
            CALL OUTPUT(T2AM(KOFF),1,NT1AM(ISYMAI),1,NT1AM(ISYMBJ),
     *                  NT1AM(ISYMAI),NT1AM(ISYMBJ),1,LUPRI)
  200    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck ccrhs_t2bt */
      SUBROUTINE CCRHS_T2BT(T2AM,WORK,LWORK,ISYM)
C
C     Alfredo Sanchez and Henrik Koch 30-July 1994
C
C     PURPOSE:
C             Back transform t2 amplitudes.
C             The amplitudes are assumed to be a square matrix.
C
#include "implicit.h"
      PARAMETER(ONETHD = 1.0D0/3.0D0,TWOTHD = 2.0D0/3.0D0)
      DIMENSION T2AM(*)
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
C----------------------------------
C     Back transform t2-amplitudes.
C----------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMAI = MULD2H(ISYMBJ,ISYM)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
C
                  DO 140 ISYMI = 1,ISYMJ
C
                     ISYMA  = MULD2H(ISYMI,ISYMAI)
                     ISYMAJ = MULD2H(ISYMA,ISYMJ)
                     ISYMBI = MULD2H(ISYMB,ISYMI)
C
                     KSCR1 = 1
                     KSCR2 = KSCR1 + NVIR(ISYMA)
                     KEND1 = KSCR2 + NVIR(ISYMA)
                     LWRK1 = LWORK - KEND1
                     IF (LWRK1 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCRHS_T2TR')
                     ENDIF
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NRHFI = J - 1
                     ELSE
                        NRHFI = NRHF(ISYMI)
                     END IF
C
                     DO 150 I = 1,NRHFI
C
                        NBI = IT1AM(ISYMB,ISYMI)+NVIR(ISYMB)*(I-1)+B
C
                        NAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1)
     *                        + IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+1
C
                        NAJBI = IT2SQ(ISYMAJ,ISYMBI)
     *                        + NT1AM(ISYMAJ)*(NBI-1)
     *                        + IT1AM(ISYMA,ISYMJ)+NVIR(ISYMA)*(J-1)+1
C
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAIBJ),1,
     *                             WORK(KSCR1),1)
                        CALL DCOPY(NVIR(ISYMA),T2AM(NAJBI),1,
     *                             WORK(KSCR2),1)
C
                        CALL DSCAL(NVIR(ISYMA),TWOTHD,T2AM(NAIBJ),1)
                        CALL DSCAL(NVIR(ISYMA),TWOTHD,T2AM(NAJBI),1)
C
                        CALL DAXPY(NVIR(ISYMA),ONETHD,WORK(KSCR2),1,
     *                             T2AM(NAIBJ),1)
                        CALL DAXPY(NVIR(ISYMA),ONETHD,WORK(KSCR1),1,
     *                             T2AM(NAJBI),1)
C
  150               CONTINUE
C
  140             CONTINUE
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      IF (IPRINT .GT. 120) THEN
         CALL AROUND('Back-transformed t2am')
         DO 200 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYM)
            KOFF = IT2SQ(ISYMAI,ISYMBJ) + 1
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Symmetry block:',ISYMBJ
            CALL OUTPUT(T2AM(KOFF),1,NT1AM(ISYMAI),1,NT1AM(ISYMBJ),
     *                  NT1AM(ISYMAI),NT1AM(ISYMBJ),1,LUPRI)
  200    CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck cc_mtcme */
      SUBROUTINE CC_MTCME(SCRM,WORK,LWORK,ISYMD,ISYMTR)
C
C     Alfredo Sanchez and Henrik Koch 26-July 1994
C     General non. total sym. Ampl. Ove Christiansen 15-2-1994.
C
C     PURPOSE:
C             Calculate 2 coublomb minus exchange of the T2M-amplitudes.
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION SCRM(*)
      DIMENSION WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      ISYMM = MULD2H(ISYMD,ISYMTR)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMCI = MULD2H(ISYMJ,ISYMM)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMI = 1,ISYMJ
C
               ISYMC  = MULD2H(ISYMI,ISYMCI)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
C
               KSCR1 = 1
               KSCR2 = KSCR1 + NVIR(ISYMC)
               KEND1 = KSCR2 + NVIR(ISYMC)
               LWRK1 = LWORK - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient space in CCSD_T2MTP')
               ENDIF
C
               IF (ISYMI .EQ. ISYMJ) THEN
                  NRHFI = J - 1
               ELSE
                  NRHFI = NRHF(ISYMI)
               END IF
C
               DO 130 I = 1,NRHFI
C
                  NCIJ = IT2BCD(ISYMCI,ISYMJ) + NT1AM(ISYMCI)*(J-1)
     *                 + IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I-1) + 1
                  NCJI = IT2BCD(ISYMCJ,ISYMI) + NT1AM(ISYMCJ)*(I-1)
     *                 + IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J-1) + 1
C
                  CALL DCOPY(NVIR(ISYMC),SCRM(NCIJ),1,WORK(KSCR1),1)
                  CALL DCOPY(NVIR(ISYMC),SCRM(NCJI),1,WORK(KSCR2),1)
                  CALL DSCAL(NVIR(ISYMC),TWO,SCRM(NCIJ),1)
                  CALL DSCAL(NVIR(ISYMC),TWO,SCRM(NCJI),1)
                  CALL DAXPY(NVIR(ISYMC),-ONE,WORK(KSCR2),1,
     *                       SCRM(NCIJ),1)
                  CALL DAXPY(NVIR(ISYMC),-ONE,WORK(KSCR1),1,
     *                       SCRM(NCJI),1)
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsd_index */
      SUBROUTINE CCSD_INDEX(INDV1,INDV2,ISYMAB)
C
C     Written by Henrik Koch 17-aug-1994.
C
C
#include "implicit.h"
#include "maxorb.h"
      DIMENSION INDV1(*), INDV2(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "symsq.h"
C
      NAB = 0
      DO 100 ISYMB = 1,NSYM
C
         ISYMA = MULD2H(ISYMB,ISYMAB)
C
         IF (ISYMA .GT. ISYMB) GOTO 100
C
         NTOTA = NBAS(ISYMA)
C
         DO 110 B = 1,NBAS(ISYMB)
C
            IF (ISYMAB .EQ. 1) NTOTA = B
C
            DO 120 A = 1,NTOTA
C
               NAB = NAB + 1
C
               NRAB = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(B - 1) + A
               NRBA = IAODIS(ISYMB,ISYMA) + NBAS(ISYMB)*(A - 1) + B
C
               INDV1(NAB) = NRAB
               INDV2(NAB) = NRBA
C
  120       CONTINUE
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CCRHS_IPM(XINT,XINTP,XINTM,SCRAB,INDV1,INDV2,
     *                     ISYMAB,ISYMG,NUMG,IG1,IG2)
C
C     Written by Henrik Koch 17-aug-1994.
C
C
C     Purpose: Making plus and minus combination of integrals.
C              (a>=g|bd) -> K+ og K- where
C                           K+- = (ag|bd) +- (bg|ad) a<=b,g<=d
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),XINTP(*),XINTM(*),SCRAB(*)
      DIMENSION INDV1(*), INDV2(*)
#include "ccorb.h"
#include "ccsdsym.h"
#include "symsq.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYDIS = MULD2H(ISYMAB,ISYMG)
C
C
      DO 100 G = IG1,IG2
C
         IG = G - IG1 + 1
C
         DO 110 ISYMB = 1,NSYM
C
            ISYMA  = MULD2H(ISYMB,ISYMAB)
            ISYMAG = MULD2H(ISYMA,ISYMG)
C
            NTOTA  = MAX(NBAS(ISYMA),1)
            NTOTAG = MAX(NNBST(ISYMAG),1)
C
            DO 120 A = 1,NBAS(ISYMA)
C
               IF (ISYMA .EQ. ISYMG) THEN
                  KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                  + INDEX(G,A)
               ELSE IF (ISYMA .LT. ISYMG) THEN
                  KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                  + NBAS(ISYMA)*(G - 1) + A
               ELSE
                  KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                  + NBAS(ISYMG)*(A - 1) + G
               ENDIF
C
               KOFF2 = IAODIS(ISYMA,ISYMB) + A
C
               CALL DCOPY(NBAS(ISYMB),XINT(KOFF1),NTOTAG,
     *                    SCRAB(KOFF2),NTOTA)
C
  120       CONTINUE
C
  110    CONTINUE
C
         KOFF = NNBST(ISYMAB)*(IG - 1)
C
#if !defined (SYS_CRAY)
         DO 130 I = 1,NNBST(ISYMAB)
C
            XINTP(KOFF + I) = SCRAB(INDV1(I))
            XINTM(KOFF + I) = SCRAB(INDV2(I))
C
  130    CONTINUE
#else
         CALL GATHER(NNBST(ISYMAB),XINTP(KOFF + 1),SCRAB,INDV1)
         CALL GATHER(NNBST(ISYMAB),XINTM(KOFF + 1),SCRAB,INDV2)
#endif
C
  100 CONTINUE
C
C
      NTOT = NNBST(ISYMAB)*NUMG
C
      CALL DAXPY(NTOT,ONE,XINTM,1,XINTP,1)
      CALL DSCAL(NTOT,-TWO,XINTM,1)
      CALL DAXPY(NTOT,ONE,XINTP,1,XINTM,1)
C
      RETURN
      END
C  /* Deck ccrhs_cio */
      SUBROUTINE CCRHS_CIO(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                     ISYVEC,ISYCIM,LUC,CFIL,IV,IOPT)
C
C     asm 17-aug-1994
C
C     Ove Christiansen 30-7-1995: modified to account for general
C                                 non. total symmetric vectors (ISYVEC) and
C                                 intermediates (ISYCIM). LUC and CFIL is
C                                 used to control from which file the
C                                 intermediate is obtained.
C
C                                 if iopt = 1 the C intermediate is assumed
C                                    to be as in energy calc.
C
C                                 if iopt ne. 1 we use the intermediate
C                                    on luc with address given according to
C                                    transformed vector nr iv (iv is not 1
C                                    if several vectors are transformed
C                                    simultaneously.)
C
C                                 in energy calc: iv=1,iopt=1
C
C     PURPOSE:
C             Calculate the C-term making I/O
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*)
      DIMENSION WORK(LWORK)
      CHARACTER CFIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*)
     &        'I/O in C-term not implemented for square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_CIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYCIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYCIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
C
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for '//
     &           'allocation in CCRHS_CIO-1')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR2 + NDISAI*NT1AO(ISYMDK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for '//
     &              'allocation in CCRHS_CIO-2')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUC,CFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C-----------------------------------------
C              Transform delta index to c.
C-----------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 190 I = 1,NRHF(ISYMI)
C
                  DO 200 A = 1,NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     IF ((NAI .LT. IFSTAI) .OR. (NAI .GT. ILSTAI))
     *                  GOTO 200
C
                     DO 210 ISYMJ = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMJ,ISYMBJ)
                        ISYMAJ = MULD2H(ISYMA,ISYMJ)
                        ISYMBI = MULD2H(ISYMB,ISYMI)
C
                        DO 220 J = 1,NRHF(ISYMJ)
C
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J-1) + A
C
                           CALL CC_PUTC(WORK(KSCR2),OMEGA2,ISYMAI,
     *                                  ISYMAJ,ISYMBI,ISYMBJ,ISYMB,
     *                                  ISYMI,ISYMJ,NAI,NAJ,I,J,
     *                                  IFSTAI)
C
  220                   CONTINUE
  210                CONTINUE
  200             CONTINUE
  190          CONTINUE
  180       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_putc */
      SUBROUTINE CC_PUTC(SCR2,OMEGA2,ISYMAI,ISYMAJ,ISYMBI,ISYMBJ,
     *                   ISYMB,ISYMI,ISYMJ,NAI,NAJ,I,J,IFSTAI)
C
C     Ove Christiansen 30-10-1995: Put in C contribution in omega vector
C                                  avoid troble on cray with optimization.
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION SCR2(*),OMEGA2(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF ( ISYMAI .EQ. ISYMBJ ) THEN
C
         DO 100 B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + INDEX(NAI,NBJ)
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                 - HALF * SCR2(KOFF9)
C
  100    CONTINUE
C
      ENDIF
C
      IF ( ISYMAI .LT. ISYMBJ ) THEN
C
         DO 200 B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMAI)*(NBJ - 1) + NAI
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                 - HALF * SCR2(KOFF9)
C
  200    CONTINUE
C
      ENDIF
C
      IF ( ISYMBJ .LT. ISYMAI ) THEN
C
         DO 300 B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *                 - HALF * SCR2(KOFF9)
C
  300    CONTINUE
C
      ENDIF
C
      IF (ISYMAJ .EQ. ISYMBI) THEN
C
         DO 400 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + INDEX(NAJ,NBI)
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
  400    CONTINUE
C
      ENDIF
C
      IF (ISYMAJ .LT. ISYMBI) THEN
C
         DO 500 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + NT1AM(ISYMAJ)*(NBI - 1)
     *            + NAJ
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
  500    CONTINUE
C
      ENDIF
C
      IF (ISYMBI .LT. ISYMAJ) THEN
C
         DO 600 B = 1,NVIR(ISYMB)
C
            NBI = IT1AM(ISYMB,ISYMI)
     *          + NVIR(ISYMB)*(I-1) + B
            NBJ = IT1AM(ISYMB,ISYMJ)
     *          + NVIR(ISYMB)*(J-1) + B
            KOFF9 = NT1AM(ISYMBJ)*(NAI - IFSTAI) + NBJ
            NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *            + NT1AM(ISYMBI)*(NAJ - 1)
     *            + NBI
            OMEGA2(NAJBI) = OMEGA2(NAJBI) - SCR2(KOFF9)
C
  600    CONTINUE
C
      ENDIF
C
      END
C  /* Deck ccrhs_dio */
      SUBROUTINE CCRHS_DIO(OMEGA2,T2AM,XLAMDH,WORK,LWORK,
     *                     ISYVEC,ISYDIM,LUD,DFIL,IV,IOPT)
C
C     asm 20-aug-1994
C
C     Ove Christiansen 30-7-1995: Modified to account for general
C                                 non. total symmetric vectors (ISYVEC) and
C                                 intermediates (ISYDIM). LUD and DFIL is
C                                 used to control from which file the
C                                 intermediate is obtained.
C
C                                 if iopt = 1 the D intermediate is assumed
C                                    to be as in energy calc.
C
C                                 if iopt ne. 1 we use the intermediate
C                                    on luc with address given according to
C                                    transformed vector nr iv (iv is not 1
C                                    if several vectors are transformed
C                                    simultaneously.)
C
C                                 in energy calc: iv=1,iopt=1
C
C     PURPOSE:
C             Calculate the D-term making I/O
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION OMEGA2(*),T2AM(*),XLAMDH(*)
      DIMENSION WORK(LWORK)
      CHARACTER DFIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
C      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (OMEGSQ) THEN
         WRITE(LUPRI,*)
     &        'I/O in D-term not implemented for square Omega2'
         CALL QUIT('OMEGSQ = .TRUE.  in CCRHS_DIO')
      END IF
C
      ISAIBJ = MULD2H(ISYVEC,ISYDIM)
C
      DO 100 ISYMAI = 1,NSYM
C
         IF (NT1AM(ISYMAI) .EQ. 0) GOTO 100
C
C
         ISYMBJ = MULD2H(ISYMAI,ISAIBJ)
         ISYMCK = MULD2H(ISYVEC,ISYMBJ)
         ISYMDK = MULD2H(ISYDIM,ISYMAI)
C
C------------------------
C        Batch structure.
C------------------------
C
         NT1AI = NT1AM(ISYMAI)
C
         LENAI  = NT1AO(ISYMDK)
         LENMIN = 2*LENAI
         IF (LENMIN .EQ. 0) GOTO 100
C
         NDISAI = LWORK / LENMIN
         IF (NDISAI .LT. 1) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_DIO')
         END IF
         NDISAI = MIN(NDISAI,NT1AI)
C
         NBATAI = (NT1AI - 1) / NDISAI + 1
C
C--------------------------
C        Loop over batches.
C--------------------------
C
         ILSTAI = 0
         DO 110 IBATAI = 1,NBATAI
C
            IFSTAI = ILSTAI + 1
            ILSTAI = ILSTAI + NDISAI
            IF (ILSTAI .GT. NT1AI) THEN
               ILSTAI = NT1AI
               NDISAI = ILSTAI - IFSTAI + 1
            END IF
C
C-----------------------------
C           Memory allocation.
C-----------------------------
C
            KSCR1 = 1
            KSCR2 = KSCR1 + NDISAI*NT1AO(ISYMDK)
            KEND  = KSCR2 + NDISAI*NT1AO(ISYMDK)
            LWRK1 = LWORK - KEND
C
            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient space for '//
     &              'allocation in CCRHS_DIO')
            END IF
C
C----------------------------------
C           Construct P(del k,#ai).
C----------------------------------
C
            KOFF1 = KSCR1
            DO 120 ISYDEL = 1,NSYM
C
               ISYMK = MULD2H(ISYDEL,ISYMDK)
C
               DO 130 IDELTA = 1,NBAS(ISYDEL)
C
                  ID = IDELTA + IBAS(ISYDEL)
C
                  IF (IOPT .EQ. 1 ) THEN
                     IOFF = IT2DEL(ID) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ELSE
                     IOFF = IT2DLR(ID,IV) + IT2BCT(ISYMK,ISYMAI)
     *                    + NRHF(ISYMK)*(IFSTAI - 1) + 1
                  ENDIF
C
                  LEN  = NDISAI*NRHF(ISYMK)
C
                  IF (LEN .GT. 0) THEN
                     CALL GETWA2(LUD,DFIL,WORK(KOFF1),IOFF,LEN)
                  ENDIF
C
                  DO 140 NAI = IFSTAI,ILSTAI
C
                     KAI = NAI - IFSTAI + 1
C
                     KOFF2 = KOFF1 + NRHF(ISYMK)*(KAI - 1)
                     KOFF3 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                     + IT1AO(ISYDEL,ISYMK) + IDELTA - 1
C
                     CALL DCOPY(NRHF(ISYMK),WORK(KOFF2),1,WORK(KOFF3),
     *                          NBAS(ISYDEL))
C
  140             CONTINUE
C
                  KOFF1 = KOFF1 + LEN
C
  130          CONTINUE
  120       CONTINUE
C
C--------------------------------------
C           Transform delta index to c.
C--------------------------------------
C
            DO 150 NAI = IFSTAI,ILSTAI
C
               KAI = NAI - IFSTAI + 1
C
               DO 160 ISYMC = 1,NSYM
C
                  ISYDEL = ISYMC
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  NBASD = MAX(NBAS(ISYDEL),1)
                  NVIRC = MAX(NVIR(ISYMC),1)
C
                  KOFF4 = ILMVIR(ISYDEL) + 1
                  KOFF5 = KSCR2 + NT1AO(ISYMDK)*(KAI - 1)
     *                  + IT1AO(ISYDEL,ISYMK)
                  KOFF6 = KSCR1 + NT1AM(ISYMCK)*(KAI - 1)
     *                  + IT1AM(ISYMC,ISYMK)
C
                  CALL DGEMM('T','N',NVIR(ISYMC),NRHF(ISYMK),
     *                       NBAS(ISYDEL),ONE,XLAMDH(KOFF4),NBASD,
     *                       WORK(KOFF5),NBASD,ZERO,WORK(KOFF6),
     *                       NVIRC)
C
  160          CONTINUE
  150       CONTINUE
C
C--------------------------------------------
C           Contract P(ck,#ai) with T(bj,ck).
C--------------------------------------------
C
            NT1BJ = MAX(NT1AM(ISYMBJ),1)
            NT1CK = MAX(NT1AM(ISYMCK),1)
C
            KOFF7 = IT2SQ(ISYMBJ,ISYMCK) + 1
C
            CALL DGEMM('N','N',NT1AM(ISYMBJ),NDISAI,NT1AM(ISYMCK),
     *                 ONE,T2AM(KOFF7),NT1BJ,WORK(KSCR1),NT1CK,
     *                 ZERO,WORK(KSCR2),NT1BJ)
C
C------------------------------
C           Scale the diagonal.
C------------------------------
C
            IF (ISYMBJ .EQ. ISYMAI) THEN
C
               DO 170 NAI = IFSTAI,ILSTAI
                  KOFF8 = KSCR2 + NT1AM(ISYMBJ)*(NAI-IFSTAI) + NAI - 1
                  WORK(KOFF8) = TWO * WORK(KOFF8)
  170          CONTINUE
C
            END IF
C
C-----------------------------------------------
C           Add the result to the packed omega2.
C-----------------------------------------------
C
            DO 180 NAI = IFSTAI,ILSTAI
C
               CALL CC_PUTD(WORK(KSCR2),OMEGA2,ISYMAI,ISYMBJ,NAI,IFSTAI)
C
  180       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cc_putd */
      SUBROUTINE CC_PUTD(SCR2,OMEGA2,ISYMAI,ISYMBJ,NAI,IFSTAI)
C
C     Ove Christiansen 30-10-1995: Put in D contribution in omega vector
C                                  avoid troble on cray with optimization.
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION SCR2(*),OMEGA2(*)
C
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccsdio.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF ( ISYMAI .EQ. ISYMBJ) THEN
         DO 190 NBJ = 1,NT1AM(ISYMBJ)
C
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + INDEX(NAI,NBJ)
C
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*SCR2(KOFF9)
C
  190    CONTINUE
C
      ENDIF
C
      IF ( ISYMAI .LT. ISYMBJ) THEN
         DO 200 NBJ = 1,NT1AM(ISYMBJ)
C
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMAI)*(NBJ - 1) + NAI
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*SCR2(KOFF9)
C
  200    CONTINUE
C
      ENDIF
C
      IF (ISYMBJ .LT. ISYMAI) THEN
         DO 210 NBJ = 1,NT1AM(ISYMBJ)
C
            KOFF9 = NT1AM(ISYMBJ)*(NAI-IFSTAI)+ NBJ
            NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *            + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
            OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + HALF*SCR2(KOFF9)
C
  210    CONTINUE
C
      ENDIF
C
      END
C  /* Deck ccrd_init */
      SUBROUTINE CCRD_INIT(KADR1,KADR2,ISYDIS)
C
C     asm 22-aug-1994
C
C     Purpose: Construct index arrays for CCRDAO
C
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
C
      DIMENSION KADR1(NBAST),KADR2(NBAST,NBAST)
C
#include "ccsdsym.h"
C
      ICOUN1 = 0
      DO 100 ISYMG = 1,NSYM
C
         ISYMAB = MULD2H(ISYMG,ISYDIS)
C
         DO 110 G = 1,NBAS(ISYMG)
            NG = IBAS(ISYMG) + G
C
            KADR1(NG) = ICOUN1
            ICOUN1 = ICOUN1 + NNBST(ISYMAB)
C
  110    CONTINUE
  100 CONTINUE
C
C
      DO 200 ISYMAB = 1,NSYM
C
         ICOUN2 = 0
         DO 210 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            IF (ISYMB .GT. ISYMA) THEN

               DO 220 B = 1,NBAS(ISYMB)
                  NB = IBAS(ISYMB) + B
C
                  DO 230 A = 1,NBAS(ISYMA)
                     NA = IBAS(ISYMA) + A
C
                     KADR2(NA,NB) = ICOUN2
                     KADR2(NB,NA) = ICOUN2
C
                     ICOUN2 = ICOUN2 + 1
C
  230             CONTINUE
  220          CONTINUE
C
            ELSE IF (ISYMA .EQ. ISYMB) THEN
C
               DO 240 B = 1,NBAS(ISYMB)
                  NB = IBAS(ISYMB) + B
C
                  DO 250 A = 1,B
                     NA = IBAS(ISYMA) + A
C
                     KADR2(NA,NB) = ICOUN2
                     KADR2(NB,NA) = ICOUN2
C
                     ICOUN2 = ICOUN2 + 1
C
  250             CONTINUE
  240          CONTINUE
C
            END IF
C
  210    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck cc2_fck */
      SUBROUTINE CC2_FCK(OMEGA2,T2AM,WORK,LWORK,ISYMTR,
     *                   XLAMDP,XLAMDH,ISIDE)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     hko 5-jan-1995
C     sym debugged 25-1-1995 oc
C     CC2 finite diff. fix - march 1997 oc
C
C     Purpose: Fock contribution in CC2 model.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
      DIMENSION OMEGA2(*),T2AM(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*)
      LOGICAL FCKCON,ETRAN
C
#include "inftap.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "leinf.h"
#include "ccsections.h"
#include "qm3.h"
!#include "qmmm.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      KSCR1 = 1
      KEND  = KSCR1 + NORBTS
      LWRK  = LWORK - KEND
C
      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient space in CC2_FCK')
      ENDIF
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(I), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KSCR1),WORK(KEND),LWRK)
C
      IF (IPRINT .GT. 80) THEN
            CALL AROUND('CC2_FCK - Orbital energies. ')
            write (LUPRI,*) (WORK(I), I=1,NORBT)
            CALL AROUND('CC2_FCK - start - : RHO2 ')
            CALL CC_PRP(DUMMY,OMEGA2,ISYMTR,0,1)
            CALL AROUND('CC2_FCK - start - : T2AM ')
            CALL CC_PRSQ(DUMMY,T2AM,ISYMTR,0,1)
      ENDIF
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 ISYMBJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 140 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 160 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ)).OR.(ISYMBJ.LT.ISYMAI))
     *                          GOTO 160
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           MAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI)*(NBJ - 1) + NAI
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)
     *         + T2AM(MAIBJ)*(WORK(MA) + WORK(MB) - WORK(MI) - WORK(MJ))
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF (((NFIELD.GT.0).OR.CCSLV)
     *     .AND.NONHF.AND.(ISIDE.NE.0)) THEN
C
         KFOCK  = 1
         KEMAT1 = KFOCK  + N2BST(ISYMOP)
         KEMAT2 = KEMAT1 + NEMAT1(ISYMOP)
         KEND1  = KEMAT2 + NMATIJ(ISYMOP)
         LWRK1  = LWORK  - KEND1
C
         CALL DZERO(WORK(KFOCK),N2BST(ISYMOP))
         CALL DZERO(WORK(KEMAT1),NEMAT1(ISYMOP))
         CALL DZERO(WORK(KEMAT2),NMATIJ(ISYMOP))
         DO 13 IF = 1, NFIELD
            FF =  EFIELD(IF)
            CALL CC_ONEP(WORK(KFOCK),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
 13      CONTINUE
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C SLV98,OC
C-------------------------------------
C
         IF (CCSLV .AND. (.NOT. CCMM )) THEN
            CALL CCSL_RHSTG(WORK(KFOCK),WORK(KEND1),LWRK1)
         ENDIF
C
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C CCMM02,JA+AO
C-------------------------------------
C
         IF (CCMM) THEN
            IF (.NOT. NYQMMM) THEN
               CALL CCMM_RHSTG(WORK(KFOCK),WORK(KEND1),LWRK1)
            ELSE IF (NYQMMM .AND. (.NOT. HFFLD)) THEN 
             !WRITE(LUPRI,*) 'About to add difference density contri'
             CALL CCMM_ADDGDIFF(WORK(KFOCK),WORK(KEND1),LWRK1)
            ELSE IF (NYQMMM .AND. HFFLD) THEN
             ! WRITE(LUPRI,*) 'You are using a hffld so no corr. needed'
              CONTINUE
            END IF
         ENDIF
C
C-------------------------------------
C
         CALL CC_FCKMO(WORK(KFOCK),XLAMDP,XLAMDH,WORK(KEND1),
     *                 LWRK1,ISYMOP,1,1)
         ETRAN  = .FALSE.
         FCKCON = .TRUE.
         ISYMEI = ISYMOP
         CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),XLAMDH,
     *                   WORK(KFOCK),WORK(KEND1),LWRK1,
     *                   FCKCON,ETRAN,ISYMEI)
C
         IF (ISIDE .EQ. -1 ) THEN
           CALL CC_EITR(WORK(KEMAT1),WORK(KEMAT2),WORK(KEND1),LWRK1,
     *                  ISYMEI)
         ENDIF
C
         CALL CCRHS_E(OMEGA2,T2AM,WORK(KEMAT1),WORK(KEMAT2),
     *                WORK(KEND1),LWRK1,ISYMTR,ISYMEI)
C
      ENDIF
C
      IF (IPRINT .GT. 80) THEN
            CALL AROUND('CC2_FCK - end - : RHO2 (RHO1=dummy')
            CALL CC_PRP(DUMMY,OMEGA2,ISYMTR,0,1)
            CALL AROUND('CC2_FCK - end - : T2AM (T1AM=dummy')
            CALL CC_PRSQ(DUMMY,T2AM,ISYMTR,0,1)
      ENDIF
C
      RETURN
      END
C  /* Deck ccrhs_efck */
      SUBROUTINE CCRHS_EFCK(EMAT1,EMAT2,XLAMDH,FOCK,WORK,LWORK,
     *                      FCKCON,ETRAN,ISYMEI)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Transforms E-intermediates to molecular basis and add Fock Matrix.
C
C     Written by Henrik Koch & Ove Christiansen 20-Jan-1994
C     Symmetry 3-aug HK, Separated from contraction OC 13-2-1995
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION EMAT1(*), EMAT2(*)
      DIMENSION WORK(LWORK),FOCK(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      LOGICAL FCKCON,ETRAN
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KSCR1  = 1
      KEND1  = KSCR1  + NMATAB(ISYMEI)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_E')
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CCRHS_EFCK1(EMAT1,EMAT2,FOCK,WORK(KSCR1),XLAMDH,
     *                 WORK(KEND1),LWRK1,FCKCON,ETRAN,ISYMEI)
C
      RETURN
      END
      SUBROUTINE CCRHS_EFCK1(EMAT1,EMAT2,FOCK,SCR1,XLAMDH,
     *                       WORK,LWORK,FCKCON,ETRAN,ISYMEI)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Transforms E-intermediates to molecular basis and add Fock Matrix.
C
C     Written by Henrik Koch & Ove Christiansen 20-Jan-1994
C     Symmetry 3-aug HK, Separated from contraction OC 13-2-1995
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DIMENSION EMAT1(*),EMAT2(*),SCR1(*)
      DIMENSION XLAMDH(*),FOCK(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      LOGICAL FCKCON,ETRAN
C
C---------------------------------------------
C     Transform the delta index of EMAT1 to c.
C---------------------------------------------
C
      IF ( ETRAN ) THEN
         DO 100 ISYMD = 1,NSYM
C
            ISYMC = ISYMD
            ISYMB = MULD2H(ISYMD,ISYMEI)
C
            NVIRB = MAX(NVIR(ISYMB),1)
            NBASD = MAX(NBAS(ISYMD),1)
C
            KOFF1 = IEMAT1(ISYMB,ISYMD) + 1
            KOFF2 = ILMVIR(ISYMD) + 1
            KOFF3 = IMATAB(ISYMB,ISYMC) + 1
C
            CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMD),
     *                 ONE,EMAT1(KOFF1),NVIRB,XLAMDH(KOFF2),NBASD,
     *                 ZERO,SCR1(KOFF3),NVIRB)
C
  100    CONTINUE
C
         CALL DSCAL(NMATAB(ISYMEI),-ONE,SCR1,1)
C
      ELSE
C
         CALL DZERO(SCR1,NMATAB(ISYMEI))
         CALL DZERO(EMAT2,NMATIJ(ISYMEI))
C
      ENDIF
C
C--------------------------------
C     Add the Fock contributions.
C--------------------------------
C
      IF (FCKCON) THEN
C
         DO 200 ISYMC = 1,NSYM
C
            ISYMB = MULD2H(ISYMC,ISYMEI)
C
            DO 210 C = 1,NVIR(ISYMC)
C
               KOFF1 = IFCVIR(ISYMB,ISYMC)  + NORB(ISYMB)*(C - 1)
     *               + NRHF(ISYMB) + 1
               KOFF2 = IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C - 1) + 1
C
               CALL DAXPY(NVIR(ISYMB),ONE,FOCK(KOFF1),1,SCR1(KOFF2),1)
C
  210       CONTINUE
  200    CONTINUE
C
         DO 220 ISYMJ = 1,NSYM
C
            ISYMK = MULD2H(ISYMJ,ISYMEI)
C
            DO 230 J = 1,NRHF(ISYMJ)
C
                KOFF1 = IFCRHF(ISYMK,ISYMJ)  + NORB(ISYMK)*(J - 1) + 1
                KOFF2 = IMATIJ(ISYMK,ISYMJ) + NRHF(ISYMK)*(J - 1) + 1
C
                CALL DAXPY(NRHF(ISYMK),ONE,FOCK(KOFF1),1,EMAT2(KOFF2),1)
C
  230       CONTINUE
  220    CONTINUE
C
      ENDIF
C
C-----------------------------------
C     Put E1 transformed back in E1.
C-----------------------------------
C
      CALL DCOPY(NMATAB(ISYMEI),SCR1,1,EMAT1,1)
C
      RETURN
      END
C  /* Deck cc_mofcon */
      SUBROUTINE CC_MOFCON(XINT,OMEGA2,XLAMDP,XLAMDH,XLAMPC,XLAMHC,
     *                     WORK,LWORK,IDEL,ISYMD,ISYMTR,IOPT,
     *                     VIJKL,CC2R12,IANR12,VAJKL,MKVAJKL,TIMR12)
C
C     Written by Asger Halkier and Henrik Koch 3-5-95.
C
C     Debugged By Ove Christiansen 25-7-1995
C
C     Purpose: To calculate the F-term's contribution to the
C              vector function using matrix vector routines.
C
C     N.B. This routine assumes AO-symmetric integrals, and can therefor
C          not be used directly for calculations with London-orbitals!!!
C
#include "implicit.h"
#include "maxorb.h"
#include "priunit.h"
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "dummy.h"
#include "ccr12int.h"
      PARAMETER(ZERO = 0.0D0,ONE = 1.0D0,XMONE=-1.0D0,TWO = 2.0D0)
      LOGICAL CC2R12,MKVAJKL,LRES
      INTEGER IANR12
      DIMENSION XINT(*),OMEGA2(*)
      DIMENSION XLAMPC(*),XLAMHC(*),XLAMDH(*),XLAMDP(*)
      DIMENSION WORK(LWORK),VIJKL(*),VAJKL(*)
      CHARACTER*8 FILBACK
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      call qenter('mofcon')
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      
      KEND0 = 1
             
      IF (CC2R12) THEN 
         KGAIJD = KEND0 
         KEND0 = KGAIJD + ND2IJG(ISYDIS)
      END IF
     
      LWRK0 = LWORK - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Lwrk0 = ',LWRK0
         CALL QUIT('Insufficient work space area in CC_MOFCON')
      ENDIF
C
      DO 100 ISYMG = 1,NSYM
C
         IF (NBAS(ISYMG) .EQ. 0) GOTO 100
C
         ISALBE = MULD2H(ISYMG,ISYDIS)
         ISYMAI = MULD2H(ISALBE,ISYMTR)
         ISYMJ  = ISYMG
C
C-----------------------------------------
C        Dynamic allocation of work space.
C-----------------------------------------
C
         KSCR1 = KEND0
         KSCR2 = KSCR1 + NNBST(ISALBE)*NRHF(ISYMJ)
         KSCR3 = KSCR2 + N2BST(ISALBE)
         KSCR4 = KSCR3 + NT1AM(ISYMAI)
         KEND1 = KSCR4 + NT1AM(ISYMAI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Lwrk1 = ',LWRK1
            CALL QUIT('Insufficient work space area in CC_MOFCON')
         ENDIF
C
C--------------------------------
C        Do first transformation.
C--------------------------------
C
         KOFF1 = IDSAOG(ISYMG,ISYDIS) + 1
         KOFF2 = ILMRHF(ISYMJ) + 1
C
         NTALBE = MAX(NNBST(ISALBE),1)
         NTOTG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NNBST(ISALBE),NRHF(ISYMJ),NBAS(ISYMG),
     *              ONE,XINT(KOFF1),NTALBE,XLAMDH(KOFF2),NTOTG,
     *              ZERO,WORK(KSCR1),NTALBE)

C---------------------------------------------------------
C                 compute contributions to V(alpha j,kl)
C---------------------------------------------------------
         IF (MKVAJKL .AND. IANR12.EQ.1) THEN
          DTIME = SECOND() 
          IF (NBAS(ISYMG).GT.0 .OR. NRHF(ISYMJ).GT.0) THEN
            KGABJD = KSCR1
            KSCR5 = KGABJD + NNBST(ISALBE)*NRHF(ISYMJ)
            KEND2 = KSCR5 + NBAST*NBAST
            LWRK2 = LWORK - KEND2

            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient space in CC_MOFCON')
            END IF

            KOFF1 = 1 + IDSAOG(ISYMG,ISYDIS)
            FILBACK = FNBACK
            CALL R12MKVAMKL(FILBACK,WORK(KGABJD),WORK(KGABJD),VAJKL,
     &           XLAMDH,1,DUMMY,DUMMY,XINT(KOFF1),XINT(KOFF1),
     &           IDEL,ISYMD,ISYMJ,
     &           ISALBE,ISYMG,WORK(KSCR5),IDUMMY,IGLMRHS,
     &           NGLMDS,WORK(KEND2),LWRK2)
          END IF
          TIMR12 = TIMR12 + (SECOND()-DTIME)
         END IF
C
C-----------------------------------
C        Last index transformations.
C-----------------------------------
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            KOFF1 = KSCR1 + NNBST(ISALBE)*(J - 1)
C
            CALL CCSD_SYMSQ(WORK(KOFF1),ISALBE,WORK(KSCR2))
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMBE = ISYMI
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMA  = MULD2H(ISYMAL,ISYMTR)
C
               KSCR5 = KEND1 
               KEND2 = KSCR5 + NBAS(ISYMAL)*NRHF(ISYMI)
               LWRK2 = LWORK - KEND2 
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient space for 2. trf. '//
     &                 'in CC_MOFCON')
               ENDIF
C
               KOFF2 = KSCR2 + IAODIS(ISYMAL,ISYMBE)
               KOFF3 = ILMRHF(ISYMI) + 1
               KOFF4 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF5 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
               NTOTAL = MAX(NBAS(ISYMAL),1)
               NTOTBE = MAX(NBAS(ISYMBE),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF2),NTOTAL,XLAMDH(KOFF3),NTOTBE,
     *                    ZERO,WORK(KSCR5),NTOTAL)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMAL),
     *                    ONE,XLAMPC(KOFF4),NTOTAL,WORK(KSCR5),NTOTAL,
     *                    ZERO,WORK(KOFF5),NTOTA)
C
C              -----------------------------------------
C                 save g(aijd) as three index array 
C              ----------------------------------------
C
               IF (CC2R12.AND.IANR12.EQ.1) THEN
                  DTIME = SECOND()

                  ISYMIJ = MULD2H(ISYMI,ISYMJ)
                  DO I = 1, NRHF(ISYMI)
                     DO A = 1, NBAS(ISYMAL)
                        IDXAI = NBAS(ISYMAL)*(I-1)+A
                        IDXIJ = IMATIJ(ISYMI,ISYMJ)+NRHF(ISYMI)*(J-1)+I
                        IDXAIJ = ID2IJG(ISYMIJ,ISYMAL)+
     &                          NBAS(ISYMAL)*(IDXIJ-1)+A 
                        WORK(KGAIJD-1+IDXAIJ) = WORK(KSCR5-1+IDXAI) 
                     END DO
                  END DO

                  TIMR12 = TIMR12 + (SECOND()-DTIME)
               END IF          

               IF (IOPT .EQ. 2) THEN
C
                  ISYMBE = MULD2H(ISYMI,ISYMTR)
                  ISYMAL = MULD2H(ISYMBE,ISALBE)
                  ISYMA  = ISYMAL
C
                  IF (LWRK1 .LT. NBAS(ISYMAL)*NRHF(ISYMI)) THEN
                     CALL QUIT('Insufficient space for 2. '//
     &                    'trf. in CC_MOFCON')
                  ENDIF
C
                  KOFF2 = KSCR2 + IAODIS(ISYMAL,ISYMBE)
                  KOFF3 = IGLMRH(ISYMBE,ISYMI) + 1
                  KOFF4 = ILMVIR(ISYMA) + 1
                  KOFF5 = KSCR3 + IT1AM(ISYMA,ISYMI)
C
                  NTOTAL = MAX(NBAS(ISYMAL),1)
                  NTOTBE = MAX(NBAS(ISYMBE),1)
                  NTOTA  = MAX(NVIR(ISYMA),1)
C
                  CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     *                       NBAS(ISYMBE),ONE,WORK(KOFF2),NTOTAL,
     *                       XLAMHC(KOFF3),NTOTBE,ZERO,WORK(KEND1),
     *                       NTOTAL)
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *                       NBAS(ISYMAL),ONE,XLAMDP(KOFF4),NTOTAL,
     *                       WORK(KEND1),NTOTAL,ONE,WORK(KOFF5),NTOTA)
C
               ENDIF
C

  120       CONTINUE
C
C--------------------------------------------------
C           Storing the result in the omega2-array.
C--------------------------------------------------
C
            ISYMB  = ISYMD
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO 130 B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
               NDB = ILMVIR(ISYMB) + NBAS(ISYMD)*(B - 1)
     *             + IDEL - IBAS(ISYMD)
C
               CALL DZERO(WORK(KSCR4),NT1AM(ISYMAI))
C
               XLB  = XLAMDP(NDB)
C
               CALL DAXPY(NT1AM(ISYMAI),XLB,WORK(KSCR3),1,WORK(KSCR4),1)
C
               IF (ISYMBJ .EQ. ISYMAI) THEN
C
                  NTOTAI = NBJ
C
                  IF (IOPT .EQ. 2) THEN
                     NTOTAI = NT1AM(ISYMAI)
                     WORK(KSCR4+NBJ-1) = TWO*WORK(KSCR4+NBJ-1)
                  ENDIF
C
                  DO 140 NAI = 1,NTOTAI
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
C
  140             CONTINUE
C
               ENDIF
C
               IF (ISYMAI .LT. ISYMBJ) THEN
C
                  DO 150 NAI = 1,NT1AM(ISYMAI)
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMAI)*(NBJ - 1) + NAI
C
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
C
  150             CONTINUE
C
               ENDIF
C
               IF ((ISYMBJ .LT. ISYMAI) .AND. (IOPT .EQ. 2)) THEN
C
                  DO 160 NAI = 1,NT1AM(ISYMAI)
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                     + NT1AM(ISYMBJ)*(NAI - 1) + NBJ
C
                     OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + WORK(KSCR4+NAI-1)
C
  160             CONTINUE
C
               ENDIF
C
  130       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE

      IF (CC2R12.AND.IANR12.EQ.1) THEN
         DTIME = SECOND()
         FACTERM23 = ONE
         CALL CC_R12MKVKL(WORK(KGAIJD),VIJKL,FACTERM23,XLAMDH,IGLMRH,
     &                    ISYMD,
     &                    ISYMTR,IDEL,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     &                    IDUMMY,.FALSE.,
     &                    WORK(KEND0),LWRK0,FNBACK)

         TIMR12 = TIMR12 + (SECOND()-DTIME)
      END IF      
 
      call qexit('mofcon')
      RETURN
      END
C  /* Deck cc_onep */
      SUBROUTINE CC_ONEP(FOCK,WORK,LWRK,FF,ISYMPT,LABPT)
C
C     Ove Christiansen 22-jan-1996.
C
C     PURPOSE:
C             Read one electron perturbation integrals
C             into FOCK AO-matrix.
C
C             If ISYMPT is input as -1 CC_ONEP returns
C             ISYMPT as correct irrep of operator given
C             by label. This thus assumes that FOCK is allocated
C             as max dim = n2bst(1).
C
C     Asger Halkier 6/2 - 1995: Fieldstrength now passed to
C     routine through the variable FF.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "iratdef.h"
      LOGICAL LOCDBG, LSYUNK
      PARAMETER (LOCDBG = .FALSE.)
      DIMENSION FOCK(*),WORK(*)
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "symsq.h"
C
      CHARACTER LABPT*(*)
C
      LSYUNK = .FALSE.
      IF (ISYMPT .EQ.-1) THEN
         LSYUNK =.TRUE.
         ISYMPT = 1
      ENDIF
C
C
      IF (IPRINT .GT. 20 ) THEN
         DN = DDOT(N2BST(ISYMPT),FOCK,1,FOCK,1)
         WRITE(LUPRI,*) 'IN ONEP: FOCK in norm:',DN
      ENDIF
C
      K2    = 1
      KEND1 = K2    + N2BST(ISYMPT)
      LEND1 = LWRK  - KEND1
C
      IF (LEND1 .LT. 0 )CALL QUIT('Insufficient space in CC_ONEP')
C
      CALL DZERO(WORK(K2),N2BST(ISYMPT))
      IERR = -1
      CALL CCPRPAO(LABPT,.TRUE.,WORK(K2),IRREP,ISYM,IERR,
     &             WORK(KEND1),LEND1)
      IF (IERR.GT.0) THEN
         CALL QUIT('CC_ONEP: I/O error while reading operator '
     &                             //LABPT(1:8))
      ELSE IF (IERR.LT.0) THEN
        CALL DZERO(WORK(K2),N2BST(ISYMPT))
      ELSE IF ((IERR.EQ.0 .AND. IRREP.NE.ISYMPT).AND.(.NOT.LSYUNK)) THEN
         CALL QUIT('CC_ONEP: symmetry mismatch for operator '
     &                             //LABPT(1:8))
      END IF
C
      IF (IPRINT .GT. 50 .OR. LOCDBG) THEN
         CALL AROUND( ' In CC_ONEP: one el. pert. integrals')
         CALL CC_PRFCKAO(WORK(K2),IRREP)
      ENDIF
C
      IF (LSYUNK) ISYMPT = IRREP
C
      CALL DAXPY(N2BST(ISYMPT),FF,WORK(K2),1,FOCK,1)
C
      IF (IPRINT .GT. 50 ) THEN
         CALL AROUND( ' In CC_ONEP: Fock AO matrix with oneel. pert')
         CALL CC_PRFCKAO(FOCK,ISYMPT)
      ENDIF
C
      IF (IPRINT .GT. 20 ) THEN
         DN = DDOT(N2BST(ISYMPT),FOCK,1,FOCK,1)
         WRITE(LUPRI,*) 'IN ONEP: FOCK out norm:',DN
      ENDIF
C
      RETURN
      END
C  /* Deck cc_bf */
      SUBROUTINE CC_BF(XINT,OMEGA2,XLAMD1,ISYML1,XLAMD2,
     *                 ISYML2,XLAMD3,ISYML3,
     *                 SCRM,ISYMM1,SCRM2,ISYMM2,WORK,LWORK,
     *                 IDEL,ISYMD,IOPT)
C
C     Written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 18-July-1994
C     Generalized by Asger Halkier and Henrik Koch 19/9 - 1995
C     to handle left-hand-side transformation contribution as well.
C     Righthand generalizations and debugging Ove Christiansen 23-9-1995
C
C     Ove Christiansen 24-9-1996: Generalization for calculating
C           terms similar to B and F-terms in the transformation
C           of vectors with the F-matrix. 
C
C
C     Purpose: Calculate B-term and F-term in the orthonormal basis.
C
C     IOPT equals one for energy-calculations and two or three for
C     response calculations (2 for left trans. and 3 for right trans.)
C     IOPT eq. 4 for F*vector contributions.
C
C     XLAMD1 is always a true lamda matrix whereas XLAMD2
C     is an AO transformed trialvector in the case af a
C     response calculation.
C
C
C     24-9-1996:
C
C     IF (IOPT .EQ. 1):
C                       scrm is ordinary t2: tci,j(delta)
C                       XLAMD1 and XLAMD2 is ordinary lamda Hole
C                       matrices.
C                       (XLAMD1(gam,i)*XLAMD1(del,j))
C
C     IF (IOPT .EQ. 2/3) 
C                       scrm is left/right vector transformed 
C                       to tci,j(delta): vector general symmetry
C                       lambda particle/hole matrix is tot.sym.
C                       XLAMD1 is ordinary lambda particle/hole matrix.
C                       XLAMD2 is transformed (barred) 
C                       lambda particle/hole matrix.
C                       (XLAMD1(gam,i)*XLAMD2(del,j)
C                       +XLAMD2(gam,i)*XLAMD1(del,j))
C
C     IF (IOPT .EQ. 4) 
C                       scrm is left/right vector transformed 
C                       to tci,j(delta): vector general symmetry
C                       lambda particle matrix is transformed.
C
C                       scrm2 is left/right vector transformed 
C                       to tci,j(delta): vector general symmetry
C                       lambda particle matrix is tot.sym. ordinary
C                       lambda particle matrixes.
C
C                       Total transformed vector to be contracted 
C                       with integrals is therefore
C 
C                       XLAMD1 is an ordinary lambda particle matrix.
C                       XLAMD2 is a double transformed 
C                              lambda particle matrix.
C                              (both R1 and L1)
C                       XLAMD3 is R1-transformed lambda particle matrix. 
C
C                       sum(gam,del)(T(gam-bar,i,j,del)+T(gam,i,j,del-bar)
C                                   +3(XLAMD1(del,j)*XLAMD2(gam,i)
C                                     +XLAMD1(del,j)*XLAMD2(gam,i)))
C
C     The symmetry input to this routine is somewhat redundant but
C     hopefully logical and flexible:
C     Isymm1 is symmetry of SCRM
C     Isymm2 is symmetry of SCRM2
C     Isyml1 is symmetry of XLAMD1
C     Isyml2 is symmetry of XLAMD2
C     Isyml3 is symmetry of XLAMD3
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),OMEGA2(*),XLAMD1(*),XLAMD2(*),XLAMD3(*)  
      DIMENSION SCRM(*),SCRM2(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      ISYMGD = MULD2H(ISYMM1,ISYML1)
C
      KMGD   = 1
      KEND1  = KMGD   + NT2BGD(ISYMGD)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CC_BF')
      ENDIF
C
C-----------------------------
C     Prepare the data arrays.
C-----------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMCI = MULD2H(ISYMJ,ISYMM1)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMC  = MULD2H(ISYMI,ISYMCI)
            ISYMG  = MULD2H(ISYMC,ISYML1)
            ISYMGI = MULD2H(ISYMG,ISYMI)
C
            NVIRC = MAX(NVIR(ISYMC),1)
            NBASG = MAX(NBAS(ISYMG),1)
C
            KOFF1 = IGLMVI(ISYMG,ISYMC) + 1
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               KOFF2 = IT2BCD(ISYMCI,ISYMJ) + IT1AM(ISYMC,ISYMI)
     *               + NT1AM(ISYMCI)*(J - 1) + 1
               KOFF3 = IT2BGD(ISYMGI,ISYMJ) + IT1AO(ISYMG,ISYMI)
     *               + NT1AO(ISYMGI)*(J - 1) + 1
C
               CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC),
     *                    ONE,XLAMD1(KOFF1),NBASG,SCRM(KOFF2),NVIRC,
     *                    ZERO,WORK(KOFF3),NBASG)
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
C---------------------------------------------------------
C     Calculate extra contribution to T2 double AO transf.
C     if F-matrix transformation.
C---------------------------------------------------------
C
      IF (IOPT .EQ. 4) THEN
C
         IF (MULD2H(ISYML3,ISYMM2).NE.ISYMGD) THEN
            CALL QUIT('CC_BF: Symmetry mismatch')
         ENDIF
         DO 200 ISYMJ = 1,NSYM
C
            ISYMCI = MULD2H(ISYMJ,ISYMM2)
C
            DO 210 ISYMI = 1,NSYM
C
               ISYMC  = MULD2H(ISYMI,ISYMCI)
               ISYMG  = MULD2H(ISYMC,ISYML3)
               ISYMGI = MULD2H(ISYMG,ISYMI) 
C
               NVIRC = MAX(NVIR(ISYMC),1)
               NBASG = MAX(NBAS(ISYMG),1)
C
               KOFF1 = IGLMVI(ISYMG,ISYMC) + 1
C
                  DO 220 J = 1,NRHF(ISYMJ)
C
                  KOFF2 = IT2BCD(ISYMCI,ISYMJ) + IT1AM(ISYMC,ISYMI)
     *                  + NT1AM(ISYMCI)*(J - 1) + 1
                  KOFF3 = IT2BGD(ISYMGI,ISYMJ) + IT1AO(ISYMG,ISYMI)
     *                  + NT1AO(ISYMGI)*(J - 1) + 1
C
                  CALL DGEMM('N','N',NBAS(ISYMG),NRHF(ISYMI),NVIR(ISYMC)
     *                      ,ONE,XLAMD3(KOFF1),NBASG,SCRM2(KOFF2),NVIRC,
     *                       ONE,WORK(KOFF3),NBASG)
C
  220          CONTINUE
C
  210       CONTINUE
C
  200    CONTINUE
C
      ENDIF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      CALL CC_BF1(XINT,OMEGA2,WORK(KMGD),ISYMGD,XLAMD1,ISYML1,
     *            XLAMD2,ISYML2,WORK(KEND1),LWRK1,
     *            IDEL,ISYMD,IOPT)
C
      RETURN
      END
C  /* Deck cc_bf1 */
      SUBROUTINE CC_BF1(XINT,OMEGA2,XMGD,ISYMGD,XLAMD1,ISYML1,
     *                  XLAMD2,ISYML2,WORK,LWORK,
     *                  IDEL,ISYMD,IOPT)
C
C     Written by Henrik Koch 3-Jan-1994
C
C     Purpose: Calculate B-term.
C 
C     See CC_BF( for more info.
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
      PARAMETER(FOURTH = 0.25D0, TWO = 2.0D0, THREE = 3.0D0)
      DIMENSION XINT(*),OMEGA2(*),XMGD(*),XLAMD1(*),XLAMD2(*)
      DIMENSION WORK(LWORK)
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYDIS = MULD2H(ISYMOP,ISYMD)
      ISYRES = MULD2H(ISYDIS,ISYMGD)
      ISYCH  = MULD2H(ISYML2,ISYMD)
C
      IF (ISYML1 .NE. 1) CALL QUIT('CC_BF: Symmetry of '//
     &     'XLAMD1 must be 1')
      IF (ISYML2 .NE. MULD2H(ISYMGD,ISYMD)) 
     *            CALL QUIT('Symmetry mismatch in CC_BF1')
C
C================================
C     Calculate the contribution.
C================================
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYRES)
         ISYMG  = MULD2H(ISYMAB,ISYDIS)
         D      = IDEL - IBAS(ISYMD)
C
         KSCRAB = 1
         KINDV1 = KSCRAB + N2BST(ISYMAB)
         KINDV2 = KINDV1 + (NNBST(ISYMAB) - 1)/IRAT + 1
         KEND1  = KINDV2 + (NNBST(ISYMAB) - 1)/IRAT + 1
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space in CC_BF1')
         ENDIF
C
C--------------------------------
C        Calculate index vectors.
C--------------------------------
C
         CALL CCSD_INDEX(WORK(KINDV1),WORK(KINDV2),ISYMAB)
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         NSIZE  = 2*(NNBST(ISYMAB) + NMIJP(ISYMIJ))
C
         IF ((NNBST(ISYMAB) .EQ. 0) .OR.
     *       (NMIJP(ISYMIJ) .EQ. 0)) GOTO 100
C
         IF (ISYMG .EQ. ISYMD) THEN
            IMAXG = D
         ELSE IF (ISYMG .LT. ISYMD) THEN
            IMAXG = NBAS(ISYMG)
         ELSE
            GOTO 100
         ENDIF
C
         IF (IMAXG.EQ.0) GOTO 100
C
         IF (LWRK1.LT.NSIZE) THEN
           CALL QUIT('Insufficient memory in CC_BF1.')
         END IF
C
         NMAXG  = MIN(IMAXG,LWRK1/NSIZE)
         NBATCH = (IMAXG - 1)/NMAXG + 1
C
         DO 110 IBATCH = 1,NBATCH
C
            NUMG = NMAXG
            IF (IBATCH .EQ. NBATCH) THEN
               NUMG = IMAXG - NMAXG*(NBATCH - 1)
            ENDIF
C
            IG1 = NMAXG*(IBATCH - 1) + 1
            IG2 = NMAXG*(IBATCH - 1) + NUMG
C
            KINTP = KEND1
            KINTM = KINTP + NNBST(ISYMAB)*NUMG
            KT2MP = KINTM + NNBST(ISYMAB)*NUMG
            KT2MM = KT2MP + NUMG*NMIJP(ISYMIJ)
            KEND2 = KT2MM + NUMG*NMIJP(ISYMIJ)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient space in CC_BF1')
            ENDIF
C
C-----------------------------------
C           Construct T2MP and T2MM.
C-----------------------------------
C
            DO 200 ISYMJ = 1,NSYM
C
               ISYMI  = MULD2H(ISYMJ,ISYMIJ)
               ISYMGI = MULD2H(ISYMI,ISYMG)
               ISYMGJ = MULD2H(ISYMJ,ISYMG)
C
               IF (ISYMI .GT. ISYMJ) GOTO 200
C
               NTOTI = NRHF(ISYMI)
C
               DO 210 J = 1,NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) NTOTI = J
C
                  DO 220 I = 1,NTOTI
C
                     NGIJ = IT2BGD(ISYMGI,ISYMJ)
     *                    + NT1AO(ISYMGI)*(J - 1)
     *                    + IT1AO(ISYMG,ISYMI)
     *                    + NBAS(ISYMG)*(I - 1) + IG1
C
                     NGJI = IT2BGD(ISYMGJ,ISYMI)
     *                    + NT1AO(ISYMGJ)*(I - 1)
     *                    + IT1AO(ISYMG,ISYMJ)
     *                    + NBAS(ISYMG)*(J - 1) + IG1
C
                     IF (ISYMI .EQ. ISYMJ) THEN
                        NIJ = IMIJP(ISYMI,ISYMJ) + INDEX(I,J)
                     ELSE
                        NIJ = IMIJP(ISYMI,ISYMJ)
     *                      + NRHF(ISYMI)*(J - 1) + I
                     ENDIF
C
                     NGIJPM = NUMG*(NIJ - 1)
C
                     KOFFP = KT2MP + NGIJPM
                     KOFFM = KT2MM + NGIJPM
C
                     IF (CC2) THEN
                        CALL DZERO(WORK(KOFFP),NUMG)
                        CALL DZERO(WORK(KOFFM),NUMG)
                     ELSE
                        CALL DCOPY(NUMG,XMGD(NGIJ),1,WORK(KOFFP),1)
                        CALL DCOPY(NUMG,XMGD(NGIJ),1,WORK(KOFFM),1)
C
                        CALL DAXPY(NUMG,ONE,XMGD(NGJI),1,WORK(KOFFP),1)
                        CALL DAXPY(NUMG,-ONE,XMGD(NGJI),1,WORK(KOFFM),1)
                     ENDIF
C
C-------------------------------------------------
C                    Add the F-term contributions.
C-------------------------------------------------
C
                     FACT = ONE
C
                     IF ((IOPT .EQ. 2) .OR. (IOPT .EQ. 4)) THEN
                        FACT = THREE
                     ENDIF
C
                     IF ((ISYMJ .EQ. ISYCH).AND.(ISYMI .EQ. ISYMG)) THEN
C
                        KOFF1 = IGLMRH(ISYMD,ISYMJ)
     &                        + NBAS(ISYMD)*(J - 1) + D
                        KOFF2 = ILMRHF(ISYMI) + NBAS(ISYMG)*(I - 1) +IG1
C
                        CALL DAXPY(NUMG,XLAMD2(KOFF1),XLAMD1(KOFF2),1,
     *                             WORK(KOFFP),1)
                        CALL DAXPY(NUMG,FACT*XLAMD2(KOFF1),
     *                             XLAMD1(KOFF2),1,WORK(KOFFM),1)
C
                     ENDIF
C
                     IF ((ISYMI .EQ. ISYCH).AND.(ISYMJ .EQ. ISYMG)) THEN
C
                        KOFF1 = IGLMRH(ISYMD,ISYMI)
     &                        + NBAS(ISYMD)*(I - 1) + D
                        KOFF2 = ILMRHF(ISYMJ) + NBAS(ISYMG)*(J - 1) +IG1
C
                        CALL DAXPY(NUMG,XLAMD2(KOFF1),XLAMD1(KOFF2),1,
     *                             WORK(KOFFP),1)
                        CALL DAXPY(NUMG,-FACT*XLAMD2(KOFF1),
     *                             XLAMD1(KOFF2),1,WORK(KOFFM),1)
C
                     ENDIF
C
C---------------------------------------------------------------
C                    For response calculation add permuted terms.
C---------------------------------------------------------------
C
                     IF (IOPT .GE. 2) THEN
C
                        ISHELP = MULD2H(ISYMG,ISYML2)
C
                        IF ((IOPT .EQ. 2) .OR. (IOPT .EQ. 4)) THEN
                           FACT = THREE
                        ENDIF
C
                        IF ((ISYMJ .EQ. ISYMD) .AND.
     &                      (ISYMI .EQ. ISHELP)) THEN
C
                           KOFF1 = ILMRHF(ISYMJ)
     &                           + NBAS(ISYMD)*(J - 1) + D
                           KOFF2 = IGLMRH(ISYMG,ISYMI)
     &                           + NBAS(ISYMG)*(I - 1) +IG1
C
                           CALL DAXPY(NUMG,XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFP),1)
                           CALL DAXPY(NUMG,FACT*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFM),1)
C
                        ENDIF
C
                        IF ((ISYMI .EQ. ISYMD) .AND.
     &                      (ISYMJ .EQ. ISHELP)) THEN
C
                           KOFF1 = ILMRHF(ISYMI)
     &                           + NBAS(ISYMD)*(I - 1) + D
                           KOFF2 = IGLMRH(ISYMG,ISYMJ)
     &                           + NBAS(ISYMG)*(J - 1) + IG1
C
                           CALL DAXPY(NUMG,XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFP),1)
                           CALL DAXPY(NUMG,-FACT*XLAMD1(KOFF1),
     &                                XLAMD2(KOFF2),1,WORK(KOFFM),1)
C
                        ENDIF
C
                     ENDIF
C
  220             CONTINUE
C
  210          CONTINUE
C
  200       CONTINUE
C
C-----------------------------------
C           Construct INTP and INTM.
C-----------------------------------
C
            CALL CCRHS_IPM(XINT,WORK(KINTP),WORK(KINTM),WORK(KSCRAB),
     *                     WORK(KINDV1),WORK(KINDV2),ISYMAB,ISYMG,
     *                     NUMG,IG1,IG2)
C
C-------------------------------
C           Scale the diagonals.
C-------------------------------
C
            IF ((ISYMG .EQ. ISYMD) .AND. (IBATCH .EQ. NBATCH)) THEN
               KOFF = KINTP + NNBST(ISYMAB)*(NUMG - 1)
               CALL DSCAL(NNBST(ISYMAB),HALF,WORK(KOFF),1)
            ENDIF
C
C----------------------------------------
C           Add the B-term contributions.
C----------------------------------------
C
            NUMGM  = MAX(NUMG,1)
            NTOTAB = MAX(NNBST(ISYMAB),1)
C
            KOFF1 = IT2ORT(ISYMAB,ISYMIJ) + 1
C
            CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                 ONE,WORK(KINTP),NTOTAB,WORK(KT2MP),NUMGM,
     *                 ONE,OMEGA2(KOFF1),NTOTAB)
C
            KOFF2 = NT2ORT(ISYRES) + IT2ORT(ISYMAB,ISYMIJ) + 1
C
            CALL DGEMM('N','N',NNBST(ISYMAB),NMIJP(ISYMIJ),NUMG,
     *                 ONE,WORK(KINTM),NTOTAB,WORK(KT2MM),NUMGM,
     *                 ONE,OMEGA2(KOFF2),NTOTAB)
C
  110    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
