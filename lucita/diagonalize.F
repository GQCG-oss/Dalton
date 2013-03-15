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

***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE CIEIG5(MV7,INICI,EROOTL,VEC1,
     &                  VEC2,MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                  LU41,NDIM,NBLK,NROOTSL,MAXVEC,MXCIIT,LUINCI,
     &                  IPRT,PEIGVC,NPRDET,H0,IPNTR,NP1,NP2,NQ,H0SCR,
     &                  EIGSHF,ICISTR,LBLK,IDIAG,VEC3,THRES_E,
     &                  NBATCH,LBATCHB,LBATCHE,LBLOCK,IBLOCK,INIDEG,
     &                  E_THRE,C_THRE,E_CONV,C_CONV,ICLSSEL,IBLK_TO_CLS,
     &                  NCLS,CLS_C,CLS_E,CLS_CT,CLS_ET,CLS_A,ICLS_L,
     &                  BLKS_A,CLS_DEL,CLS_DELT,ISKIPEI,I2BLK,ICLS_A2,
     &                  MXLNG,IROOTHOMING,IBASSPC,EBASC,CBASC,NSPC,
     &                  MULSPC,IPAT,LPAT,ISPC,IBLOCKAR,
     &                  IBLKDSTND,RCCTOS,IPROCLIST,IGROUPLIST
#ifdef VAR_MPI
     &                  ,LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                  LU6LIST,LU7LIST,LUCLIST
#endif
     &                  )

*
* Master routine for CI diagonalization
*
* Modified to handle PQ - preconditioner , May 1990
* PICO,MICDV4 added spring of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KSIOIO, KSBLTP, KSVST, KLXLBT, KLXLEBT, KLXI1BT, KLXIBT
!               for addressing of WORK
C
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LUCLIST(*)
      DIMENSION LU7LIST(*), LU6LIST(*)
#endif
#include "priunit.h"
#include "parluci.h"
#include "files.inc"
#include "mxpdim.inc"
C     Definition of c and sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
C        /ORBINP/ : NACOB used
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "gasstr.inc"
#include "cgas.inc"
C        Used: NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
      LOGICAL CA
C
C
      DIMENSION VEC1(*),VEC2(*)
C     DIMENSION INIDET(100)
      PARAMETER( LLWRK =100000)
      integer, allocatable  :: iscr1(:)
      real(8), allocatable  :: scr1(:)
*. Output from Subspace dagonalization
      DIMENSION H0(*),IPNTR(*),H0SCR(*),PEIGVC(*)
*.
      DIMENSION EROOTL(NROOTSL)
      INTEGER   RCCTOS(*)

!     use dynamic instead of static allocation (+ common block)
      allocate(iscr1(llwrk))
      allocate(scr1(llwrk))
C
C     parameter used in parallel scheme /sk march 2007/
C
      IAMSURE = 1
*
      IF( IPRT.GT. 1 )  WRITE(luwrt,'(/A)')
     &'          *** information from ci diagonalization  ***'
!     NTEST = 9999 ! NTEST=0
      NTEST = 0
C
C               ====================================
C  1 :               INITIAL VARIATIONAL SUBSPACE
C               ====================================
C
      IF( INICI .EQ. 0 ) THEN
        IF(NPRDET .EQ. 0 ) THEN
C         ==================================================
C          Initial guess from lowest elements of CI diagonal
C         ==================================================
C
C         in order treat degeneracies, the lowest 4 * NROOTSL
C         elements are obtained
C
          NFINDM = MIN(NDIM,4*NROOTSL)
C
          IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
              CALL FNDMND_PAR(LUDIA,LBLK,VEC1,NFINDM,NFINDA,
     &                    ISCR1(1+2*NFINDM),SCR1(1+2*NFINDM),ISCR1,
     &                    SCR1,IBLOCKAR,IBLKDSTND,NUM_BLOCKS2,IPRT)
#endif
          ELSE
             CALL FNDMND(LUDIA,LBLK,VEC1,NFINDM,NFINDA,
     &            ISCR1(1+2*NFINDM),SCR1(1+2*NFINDM),ISCR1,SCR1,IPRT)
          END IF
          CALL REWINE(LU1,-1)
          IBASE = 1
          TEST = 1.0D-10
          DO 100 IROOTSI = 1, NROOTSL
*. Number of degenerate elements
            NDEG = 1
            XVAL = SCR1(IBASE)
   90       CONTINUE
            IF(IBASE-1+NDEG+1.LE.NFINDA) THEN
              IF (ABS(SCR1(IBASE-1+NDEG+1)-XVAL).LE.TEST) THEN
                NDEG = NDEG + 1
                GOTO 90
              END IF
            END IF
*
            IF (INIDEG.EQ.0.AND.NDEG.GT.1) THEN
              WRITE(luwrt,'(/A,I5/A)')
     &        ' INFO !  DEGENERATE INITIAL VECTORS FOR CI ROOT',IROOTSI,
     &        ' INFO !  ONLY ONE PARTNER USED, CHECK YOUR RESULTS!'
              NDEG = 1
            END IF

*. Initial guess in compressed form in SCR1
            SCALE = NDEG
            SCALE =1.0D0 / SQRT(SCALE)
            DO 250 II = 1,NDEG
*. Anti symmetric combination
              IF(INIDEG.EQ.-1) THEN
                SCR1(II) = (-1.0D0)**II * SCALE
*. Symmetric combination
              ELSE IF (INIDEG.EQ.1.OR.INIDEG.EQ.0) THEN
                SCR1(II) =  SCALE
              END IF
  250       CONTINUE
            IF(IDIAG.EQ.2) THEN
              JPACK = 1
            ELSE
              JPACK = 0
            END IF
            IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
               CALL WRSVCD_PAR(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,
     &                      IBLOCKAR,IBLKDSTND,NUM_BLOCKS2,JPACK,
     &                      IROOTSI,LU1LIST)
#endif
            ELSE
               CALL WRSVCD(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
     &                  LUDIA,JPACK)
            END IF
            IBASE = IBASE + NDEG
  100     CONTINUE
        ELSE
* =====================================
*. Initial approximations are in PEIGVC
* =====================================
          CALL REWINE(LU1,-1)
          IF(IDIAG.EQ.2) THEN
            JPACK = 1
          ELSE
            JPACK = 0
          END IF
          DO 1984 IROOTSI = 1, NROOTSL
            CALL WRSVCD(LU1,LBLK,VEC1,IPNTR,
     &           PEIGVC((IROOTSI-1)*NPRDET+1),NPRDET,NDIM,LUDIA,JPACK)
 1984     CONTINUE
        END IF
      END IF
*
* ======================================================
* Initial CI vectors are already on file LU1, do nothing
* ======================================================
*
      IF( INICI .LT. 0 ) THEN
*.      vectors assumed already in LU1/ILU1 (parallel)
        IF(IPRT.GT. 1 )
     &  WRITE(luwrt,*)' Initial CI vector assumed in place '
      END IF
*
*                 ========================
* 2 :                  Diagonalization
*                 ========================
*
!     release scratch memory
      deallocate(iscr1)

      CALL QENTER('CIEIG')
* Inverse iteration modified Davidson with 2 vectors in core
      IF(IPRT .GE. 5 ) THEN
         WRITE(luwrt,'(A,I3)')
     &   '  Number of roots to be converged..  ',NROOTSL
         WRITE(luwrt,'(A,I3)')
     &   '  Largest allowed number of vectors..',MAXVEC
         WRITE(luwrt,'(A,I3)')
     &   '  Allowed number of CI iterations  ..',MXCIIT
      END IF

*     allocations for SCR1
      KRNRM  = 1
      KEIG   = KRNRM + MXCIIT*NROOTSL
      KFIN   = KEIG  + MXCIIT*NROOTSL
      KAPROJ = KFIN + NROOTSL
      KAVEC  = KAPROJ + MAXVEC*(MAXVEC+1)/2
      KWORK  = KAVEC + MAXVEC ** 2
      KLFREE = KWORK + MAXVEC*(MAXVEC+1)
      IF( IPRT .GE. 5 ) THEN
         WRITE(luwrt,*) ' KRNRM KEIG KFIN KAPROJ KAVEC KWORK KLFREE '
         WRITE(luwrt,'(6I8)')KRNRM,KEIG,KFIN,KAPROJ,KAVEC,KWORK,KLFREE
      END IF
      IF( KLFREE-1 .GT. LLWRK) THEN
           WRITE(luwrt,'(A,2I5)' )
     &     ' Not enough memory in CIEIG5 : neeeded and available ',
     &     KLFREE-1, LLWRK
           WRITE(luwrt,'(A,2I5)' )
     &     ' Increase parameter LLWRK in CIEIG5 to   ', KLFREE-1
           Call Abend2( ' insufficient memory in cieig5 ' )
       END IF
*
      IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
C
C     partition CI vector with max. batch length using LBLOCK
C
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'CIEIG5')
C
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
C     Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
C
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
C     Arrays giving allowed type combinations
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KSIOIO),0)
C
      KSVST = 1
C     Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
C     Arrays for partitioning of sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLXLBT ,NTTS   ,'ADDL  ' ,1,'CLBT  ')
      CALL MEMMAN(KLXLEBT ,NTTS  ,'ADDL  ' ,1,'CLEBT ')
      CALL MEMMAN(KLXI1BT,NTTS   ,'ADDL  ' ,1,'CI1BT ')
      CALL MEMMAN(KLXIBT ,8*NTTS ,'ADDL  ' ,1,'CIBT  ')
C     partition the CI-vactor in batches
      ITTSS_ORD = 2
C
      CALL PART_CIV_PAR1(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,
     &                   LBLOCK_LUCI,
     &                   WORK(KSIOIO),ISMOST(1,ISSM),
     &                   NBATCH,WORK(KLXLBT),WORK(KLXLEBT),
     &                   WORK(KLXI1BT),WORK(KLXIBT),0,ITTSS_ORD,
     &                   IBLKDSTND,IBLOCKI)
*
       I_USE_MICDV6 = 1
       CALL MICDV6_PAR(VEC1,VEC2,SCR1(KRNRM),SCR1(KEIG),
     &                 EROOTL,MXCIIT,NDIM,NROOTSL,
     &                 MAXVEC,NROOTSL,SCR1(KAPROJ),SCR1(KAVEC),
     &                 SCR1(KWORK),IPRT,NPRDET,H0,IPNTR,NP1,NP2,NQ,
     &                 H0SCR,LBLK,EIGSHF,THRES_E,IROOTHOMING,
     &                 luwrt,IBLOCKAR,IBLKDSTND,
     &                 RCCTOS,
     &                 LU1LIST,LU2LIST,LU3LIST,LU4LIST,LU5LIST,
     &                 LU6LIST,LU7LIST,LUCLIST,NBATCH,
     &                 WORK(KLXLBT),WORK(KLXLEBT),WORK(KLXI1BT),
     &                 WORK(KLXIBT),IPROCLIST,IGROUPLIST)

#endif
      ELSE ! (LUCI_NMPROC .GT. 1) THEN
       I_USE_MICDV6 = 1
       CALL MICDV6(VEC1,VEC2,LU1,LU2,SCR1(KRNRM),SCR1(KEIG),
     &             EROOTL,MXCIIT,NDIM,LU3,LU4,LU5,LU6,LU7,
     &             LUDIA,NROOTSL,
     &             MAXVEC,NROOTSL,SCR1(KAPROJ),SCR1(KAVEC),
     &             SCR1(KWORK) ,IPRT,NPRDET,H0,IPNTR,NP1,NP2,NQ,
     &             H0SCR,-1,EIGSHF,THRES_E,IROOTHOMING,luwrt)
      END IF ! (LUCI_NMPROC .GT. 1) THEN

      ENOT = SCR1(KFIN)
!     release scratch memory
      deallocate(scr1)

      CALL QEXIT('CIEIG')

#ifdef VAR_MPI
      IF (LUCI_NMPROC .GT. 1) THEN
C       eliminate local memory
        IDUM = 0
        CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'CIEIG5')
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE CLASS_TRUNC(NCLS,ICLS_L,CLS_CT,CLS_ET,CLS_C,CLS_E,
     &                       E_CONV,ICLS_A,N_TRN_CLS,E_TRUNC,W_TRUNC,
     &                       IPRNT)
*
* Decide which classes of parameters that can be eliminated
*
* Jeppe Olsen, Jan 97
*              March '97 updated
*
*. Note in current version all energy contributions are
*. positive.
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION ICLS_L(NCLS)
      DIMENSION CLS_CT(NCLS),CLS_ET(NCLS),CLS_C(NCLS),CLS_E(NCLS)
*. Output
      INTEGER ICLS_A(NCLS)
*. Giving  the truncated classes
*  Additional output is
*              N_TRN_CLS : Number of truncated classes
*              E_TRUNC   : Estimated error of eliminating these classes
*                          (compared to expansion with only largest
*                           coefficients included).
*                          E_TRUNC is thus the error that should be
*                          added to the estmated error arising from
*                          eliminating small terms
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRNT)
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' Welcome to TRUNC_CLASS '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*)
     &  ' Required threshold for convergence of energy',E_CONV
      END IF
*
*. First find total energy correction and
*. largest class contribution
*
      E_TOT = 0.0D0
      E_CLS_MAX= 0.0D0
      DO JCLS = 1, NCLS
       E_TOT = E_TOT + CLS_ET(JCLS)
       E_CLS_MAX = MAX( E_CLS_MAX,CLS_ET(JCLS))
      END DO
C?    WRITE(6,*) ' E_CLS_MAX = ', E_CLS_MAX
*
*. The truncation is done in two steps, first an overall
*  threshold for deleting classes is constructed,
*  and all classes with energy contributions less than this
*  threshold are excluded.
*. Then one takes a pass through remaining classes and
*  eliminates untill the quota is filled
*
*. Threshold for including determinants, start with E_CLS_MAX
*. and decrease until quota is filled
*
      ILOOP = 0
      STEP = 1.2
 1000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
          E_THRES =  E_CLS_MAX
        ELSE
          E_THRES = E_THRES/STEP
        END IF
        E_ELI = 0.0D0
        N_ELI = 0
        DO JCLS = 1, NCLS
          IF(ABS(CLS_ET(JCLS)).LE.E_THRES) THEN
            E_ELI = E_ELI + CLS_ET(JCLS)
            N_ELI = N_ELI + 1
          END IF
        END DO
        IF(ILOOP.EQ.10000) THEN
          WRITE(6,*) ' Loop count exceeded 10000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          Call Abend2( 'CLASS_TRUNC: Iloop.eq.10000' )
        END If
      IF(E_ELI.GT. E_CONV) GOTO 1000
      IF(NTEST.GE.10) WRITE(6,*)
     & ' Overall Threshold for eliminating classes',E_THRES
*
      N_PAS_CLS = 0
      L_PAS_CLS = 0
      N_TRN_CLS = 0
      L_TRN_CLS = 0
      N_ACT_CLS = 0
      L_ACT_CLS = 0
*
      E_TRUNC = 0.0D0
      E_TRUNCT = 0.0D0
      W_TRUNCT= 0.0D0
      W_TRUNC = 0.0D0
*
*. Eliminate classes with energy contribution less than E_THRES
*
      IONE = 1
      CALL ISETVC(ICLS_A,IONE,NCLS)
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).NE.0.0D0.AND
     &     .ABS(CLS_ET(JCLS)).LE. E_THRES) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        END IF
      END DO
      IF(NTEST.GE.10) THEN
      WRITE(6,*)
     &  ' Number of classes with contributions less than E_THRES',
     &   N_TRN_CLS
      WRITE(6,*) ' Energy contributions from these classes ',
     &   E_TRUNCT
      END IF
* Eliminate remaining classes until thres hold is obtained
      DO JCLS = 1, NCLS
        IF(CLS_CT(JCLS).EQ.0.0D0) THEN
*. Passive class, no contribution before truncation
          N_PAS_CLS = N_PAS_CLS + 1
          L_PAS_CLS = L_PAS_CLS + ICLS_L(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1
     &  .AND.ABS(E_TRUNCT+CLS_ET(JCLS)).LT.E_CONV) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        ELSE IF(ICLS_A(JCLS).EQ.1) THEN
*. Class is active
          N_ACT_CLS = N_ACT_CLS + 1
          L_ACT_CLS = L_ACT_CLS + ICLS_L(JCLS)
           ICLS_A(JCLS) = 1
        END IF
      END DO
*. Correct for missing minus in first order correction (not here, not now )
      E_TRUNC  =  E_TRUNC
      E_TRUNCT =  E_TRUNCT
*
      IF(NTEST.GE.10) THEN
      WRITE(6,*)
      WRITE(6,*) ' Estimated complete energy contribution of ',
     &            '  eliminated classes ', -E_TRUNCT
      WRITE(6,*) ' Estimated truncated energy contribution of ',
     &            '  eliminated classes ', -E_TRUNC
      WRITE(6,*) ' Truncation error in eliminated classes was ',
     &              -E_TRUNC+E_TRUNCT
      WRITE(6,*) ' Estimated truncated weight contribution of ',
     &            '  eliminated classes ', W_TRUNC
      WRITE(6,*) ' Estimated energy contribution without trunc',
     &             -E_TOT
      WRITE(6,*) ' Energy contribution of of active classes ',
     &             -E_TOT+E_TRUNCT
      END IF

C?    WRITE(6,*)
C?    WRITE(6,*) '  Class      Number     Dimension  '
C?    WRITE(6,*) ' ================================= '
C?    WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Passive  ', N_PAS_CLS,L_PAS_CLS
C?    WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Truncated', N_TRN_CLS,L_TRN_CLS
C?    WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Active   ', N_ACT_CLS,L_ACT_CLS
*
*. Eliminate classes with energy contributions less than E_TOT* FAC2
*  (temporary elimination, these classes may be invoked later in
*  the iterative sequence). Only the active classes are examined.
*
      FAC2 = 0.1
      ILOOP = 0
      STEP = 1.2
      E_TEMP_TRUNC =  ABS(E_TOT)*FAC2
 2000 CONTINUE
        ILOOP = ILOOP+1
        IF(ILOOP.EQ.1) THEN
          E_THRES2 =   E_TEMP_TRUNC
        ELSE
          E_THRES2 = E_THRES2/STEP
        END IF
        E_ELI2 = 0.0D0
        N_ELI2 = 0
        DO JCLS = 1, NCLS
          IF(ICLS_A(JCLS).EQ.1.AND.
     &       ABS(CLS_ET(JCLS)).LE.E_THRES2) THEN
            E_ELI2 = E_ELI2 + ABS(CLS_ET(JCLS))
            N_ELI2 = N_ELI2 + 1
          END IF
        END DO
        IF(ILOOP.EQ.10000) THEN
          WRITE(6,*) ' Loop count exceeded 10000'
          WRITE(6,*) ' I am afraid I am in an infinite loop'
          WRITE(6,*) ' So I will stop '
          Call Abend2( 'CLASS_TRUNC: Iloop.eq.10000' )
        END If
      IF(E_ELI2.GT. E_TEMP_TRUNC) GOTO 2000
      IF(NTEST.GE.10) THEN
      WRITE(6,*)
     & ' Temporary elimination of classes with total contribution'
      WRITE(6,*) ' less than ',E_TEMP_TRUNC
      WRITE(6,*) ' gives threshold for temporary elimination',E_THRES2
      END IF
*
*. Eliminate classes with energy contribution less than E_THRES2
*
      N_TRN_CLS2 = 0
      L_TRN_CLS2 = 0
      E_TRUNC2 = 0
      W_TRUNC2 = 0
      DO JCLS = 1, NCLS
        IF(ICLS_A(JCLS).EQ.1.AND.CLS_CT(JCLS).NE.0.0D0.AND
     &     .ABS(CLS_ET(JCLS)).LE. E_THRES2) THEN
          N_TRN_CLS = N_TRN_CLS + 1
          N_TRN_CLS2= N_TRN_CLS2+ 1
          L_TRN_CLS = L_TRN_CLS + ICLS_L(JCLS)
          L_TRN_CLS2= L_TRN_CLS2+ ICLS_L(JCLS)
          E_TRUNC  = E_TRUNC  + CLS_E (JCLS)
          E_TRUNC2 = E_TRUNC2 + CLS_E (JCLS)
          E_TRUNCT = E_TRUNCT + CLS_ET(JCLS)
          W_TRUNC  = W_TRUNC + CLS_C (JCLS)
          W_TRUNC2 = W_TRUNC2+ CLS_C (JCLS)
          W_TRUNCT = W_TRUNCT+ CLS_CT(JCLS)
          ICLS_A(JCLS) = 0
        END IF
      END DO
      IF(NTEST.GE.1) THEN
      WRITE(6,*)
      WRITE(6,'(A,F25.12)')
     &  ' Energy contributions from eliminated classes ', -E_TRUNC
      WRITE(6,'(A,F25.12)')
     &  ' Norm of eliminated classes                   ',SQRT(W_TRUNC)
      END IF
*
      N_ACT_CLS = N_ACT_CLS - N_TRN_CLS2
      L_ACT_CLS = L_ACT_CLS - L_TRN_CLS2
*
      E_TRUNC = - E_TRUNC
      E_TRUNCT = -  E_TRUNCT
*
      if (NTEST.ge.5) then
      WRITE(6,*)
      WRITE(6,*) '             Number     Dimension  '
      WRITE(6,*) ' ================================= '
      WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Passive  ', N_PAS_CLS,L_PAS_CLS
      WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Truncated', N_TRN_CLS,L_TRN_CLS
      WRITE(6,'(1H ,A,5X,I5,5X,I10)') '(Temp)   ',
     &N_TRN_CLS2,L_TRN_CLS2
      WRITE(6,'(1H ,A,5X,I5,5X,I10)') 'Active   ', N_ACT_CLS,L_ACT_CLS
*
      WRITE(6,*)
      WRITE(6,*) ' Information about classes '
      WRITE(6,*) ' =========================='
      WRITE(6,*)
      WRITE(6,*)
     & ' Class    Dimension      E          E(Trunc)      C',
     & '         C(Trunc)   Active'
      WRITE(6,*)
     & ' =============================================================',
     & '================'
       DO JCLS = 1, NCLS
         IF(ABS(CLS_ET(JCLS)).GT.0.0D0.OR.ABS(CLS_CT(JCLS)).GT.0.0D0)
     &   THEN
           WRITE(6,'(1H ,2X,I4,1X,I10,3X,4(E12.6,1X),1X,I2)')
     &     JCLS,ICLS_L(JCLS),CLS_ET(JCLS),CLS_E(JCLS),
     &     CLS_CT(JCLS),CLS_C(JCLS),ICLS_A(JCLS)
         END IF
       END DO
      end if
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
C        newden.f propint.f newrou.f  reo.f annstr.f\
C        adast.f decode_line.f rsbb2bn.f adaada.f rsbb2an.f\
C        advice_sigma.f traci.f  rsbb2bvn.f
      SUBROUTINE CLS_TO_BASE(CLS_E,EBASC,CLS_C,CBASC,NCLS,NSPC,
     &                       IBASSPC,IPRNT)
*
* Class info => Base space info
* for energy abd wf correction
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBASSPC(NCLS)
      DIMENSION CLS_C(NCLS),CLS_E(NCLS)
*. Output
      DIMENSION EBASC(*),CBASC(*)
*
      ZERO = 0.0D0
      CALL SETVEC(EBASC,NSPC,ZERO)
      CALL SETVEC(CBASC,NSPC,ZERO)
*
      DO ICLS = 1, NCLS
        ISPC = IBASSPC(ICLS)
        EBASC(ISPC) = EBASC(ISPC) + CLS_E(ICLS)
        CBASC(ISPC) = CBASC(ISPC) + CLS_C(ICLS)
      END DO
*
      NTESTL = 000
      NTEST = max(IPRNT,NTESTL)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =============================================='
        WRITE(6,*) ' Contribution  to energy and wf per base space '
        WRITE(6,*) ' =============================================='
        WRITE(6,*)
        WRITE(6,'(A)') '  Class         Energy          wf '
        WRITE(6,'(A)') ' ==========================================='
        DO ISPC = 1, NSPC
          WRITE(6,'(2X,I3,3X,E12.6,2X,E12.6)')
     &          ISPC,EBASC(ISPC),CBASC(ISPC)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIAG_BLKS(A,X,LGAS,NOBPSM,MXPOBS,NSMOB,NGAS,
     &                     SCR1,SCR2)
*
* A matrix A and an orbital partitioning LGAS is given.
* Diagonalize diagonal blocks of A
*
*. A is assumed to contain only active orbitals and is assumed to be in
* packed form and symmetry ordered
*
*
* Input
* =====
* A : Input matrix
* LGAS : Orbital partitioning
* NOBPSM : Number of orbitals per symmetry
* MXPOBS : Max number of orbital symmetries
* NSMOB  : Number of orbital symmetries
* NGAS   : Number of orbital partitionings
*
* Output
* ======
* X : Eigenvector expansion, sorted according to eigenvalues in
*     each subspace
* A : Correspinding eigenvalues
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION A(*),X(*)
      DIMENSION LGAS(MXPOBS,*)
      DIMENSION NOBPSM(*)
*. Scratch : Number of orbitals **2 ( atmost )
      DIMENSION SCR1(*), SCR2(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DIAG_BLKS'
        WRITE(6,*) ' ========='
        WRITE(6,*) ' NSMOB NGAS ', NSMOB,NGAS
      END IF
*
      DO ISMOB = 1, NSMOB
*. offset for symmetryblocks in matrices
        IF(ISMOB.EQ.1) THEN
          IOFFMTP = 1
          IOFFMTC = 1
        ELSE
          IOFFMTP = IOFFMTP + NOBPSM(ISMOB-1)*(NOBPSM(ISMOB-1)+1)/2
          IOFFMTC = IOFFMTC + NOBPSM(ISMOB-1) ** 2
        END IF
*. Zero symmetry block of eigenvector matrix to avoid interactions
*. between different blocks
        LOBPS = NOBPSM(ISMOB)
        ZERO = 0.0D0
        CALL SETVEC(X(IOFFMTC),ZERO,LOBPS**2)
*. Loop over subbloks, extract,  diagonalize, and expand
        DO ITPOB = 1, NGAS
          IF(ITPOB.EQ.1) THEN
            IOFFOB=1
          ELSE
            IOFFOB = IOFFOB + LGAS(ISMOB,ITPOB-1)
          END IF
          LOB = LGAS(ISMOB,ITPOB)
*. Extract
          IJ2 = 0
          DO IOB = IOFFOB,IOFFOB+LOB-1
            JOBMX = IOB
            DO JOB = IOFFOB,JOBMX
              IJ1 = IOFFMTP -1 + IOB*(IOB-1)/2+JOB
              IJ2 = IJ2 + 1
              SCR1(IJ2) = A(IJ1)
            END DO
          END DO
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Extracted block of matrix for  ISMOB,ITPOB = ',
     &      ISMOB,ITPOB
            CALL PRSYM(SCR1,LOB)
          END IF
*, Diagonalize
C         CALL EIGEN_LUCI(WORK(KMAT1-1+IOFFP),WORK(KMAT2-1+IOFFC),LORB,0,1)
          CALL EIGEN_LUCI(SCR1, SCR2, LOB,0,1)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Corresponding eigenvalues and eigenvectors'
            WRITE(6,*) (SCR1(I*(I+1)/2),I=1, LOB)
            WRITE(6,*)
            CALL WRTMT_LU(SCR2,LOB,LOB,LOB,LOB)
          END IF
*. Expand eigenvector to full symmetry block
          IJ2 = 0
          DO JOB = IOFFOB,IOFFOB+LOB-1
            DO IOB = IOFFOB,IOFFOB+LOB-1
              IJ1 = IOFFMTC -1 + (JOB-1)*LOBPS + IOB
              IJ2 = IJ2 + 1
              X(IJ1) = SCR2(IJ2)
            END DO
          END DO
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Complete eigenvector matrix'
        CALL APRBLM2(X,NOBPSM,NOBPSM,NSMOB,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE EKTPERT(F,S,NDIM,NORD,EN,C,
     &                   VEC1,VEC2,VEC3,AMAT1,AMAT2,AMAT3)
*
* Perturbation expansion of generalized eigenvalue problem
* Special version for EKT problem where there are singularities in
* the zero order matrices.
*
* Ordering the matrices so the occupied orbitals come first, and
* then the virtual orbitals, the  zero order matrices have the form
*
*     (x *                     )
*     ( x*               0     )
*     (**x                     )
*     (                        )
*     (                        )
*     (                        )
*     (                        )
*     (                        )
*     (  0              0      )
*     (                        )
*     (                        )
*     (                        )
*
*
*
* The matrices in the  eigenvalue problem FC = E SC
* are separated into orders
*
* F = Sum(k=0,NORD) F(K)
* S = Sum(K=0,NORD) S(K)
*
*. Obtain corrections to energy and wave functions
*
*. The normalization condition used is C(K)T S(0) C(0) = 0
*
* The energy corrections become
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
* and the wave function corrections in the occupied orbital space are
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(I=0,N-K)E(N-K-L)S(L)C(K))
*
* Whereas they read in the virtual space
*
* C(N-2) = (F(2)-E(0)S(2))-1 (-Sum(K=0,N-3)F(N-K)C(K)
*                             +Sum(K=0,N-3)Sum(I=0,N-K)E(N-K-L)S(L)C(K)
*                             + -(F(2)-E0*S(2)) C(N-2)OCC )
* Where only the
*
* The zero order matrices F(0),S(0) are assumed diagonal
*
* Jeppe and Dage, Oct 22 1997

*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION F(NDIM**2,*),S(NDIM**2,*)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(0:NORD)
*. Scratch
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
      DIMENSION AMAT1(NDIM,NDIM),AMAT2(NDIM,NDIM)
      DIMENSION AMAT3(NDIM,NDIM)
*
      NTEST = 1

*. Zero order wavefunction in virtual space
* obtained by solving - in the virtual space, the
* equations
*
* (F(2)-E(0)S(2))(virt,virt) C0(virt) = -(F(2)-E(0)S(2))C0(occ))(virt)
*
*
*. Obtain the number of occupied orbitals by examining S(0)
      NOCC = 0
      DO I = 1, NDIM
       IF(S((I-1)*NDIM+I,1) .NE. 0) NOCC = NOCC + 1
      END DO
      NVIRT = NDIM - NOCC
      WRITE(6,*) ' Number of occupied orbitals ', NOCC
      WRITE(6,*) ' Number of virtual  orbitals ', NVIRT
*
*. Zero order energy
*
C          MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(F,C,VEC1,NDIM,NDIM,0)
      C0FC0 = INPROD(VEC1,C,NDIM)
*
      CALL MATVCB(S,C,VEC1,NDIM,NDIM,0)
      C0SC0 = INPROD(VEC1,C,NDIM)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'C(0)T F(0) C(0) = ', C0FC0
        WRITE(6,*) 'C(0)T S(0) C(0) = ', C0SC0
      END IF
      E0 = C0FC0/C0SC0
      EN(0) = E0
*. Save diagonal of F(0) - E(0)C(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = F((I-1)*NDIM+I,1)-E0*S((I-1)*NDIM+I,1)
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Zero order diagonal '
        CALL WRTMT_LU(VEC3,1,NDIM,1,NDIM)
      END IF
*
*. The virtual-virtual  part of the operator F(2)-E(0)S(2)
*. will be used to solve linear equations. Set up inverse
*
      ONE = 1.0D0
      CALL VECSUM(AMAT1,F(1,3),S(1,3),ONE,-E0,NDIM**2)
      DO IOCC = 1, NOCC
        DO J = 1, NDIM
          AMAT1(IOCC,J) = 0.0D0
          AMAT1(J,IOCC) = 0.0D0
        END DO
        AMAT1(IOCC,IOCC) = 0.001
      END DO
C?    write(6,*) ' Input matrix to INVMAT'
C?    CALL WRTMT_LU(AMAT1,NDIM,NDIM,NDIM,NDIM)
      CALL COPVEC(AMAT1,AMAT3,NDIM**2)
      CALL INVMAT(AMAT1,AMAT2,NDIM,NDIM)
*     ^ Returns inverse matrix in AMAT1
      WRITE(6,*) ' Inverse matrix obtained by INVMAT'
      CALL WRTMT_LU(AMAT1,NDIM,NDIM,NDIM,NDIM)
*. Inverse matrix by diagonalization
      CALL COPVEC(AMAT3,AMAT1,NDIM**2)
C          INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
      CALL INVERT_BY_DIAG(AMAT1,AMAT2,AMAT3,VEC1,NDIM)
*. Zero occ-occ part
      DO IOCC = 1, NOCC
        DO J = 1, NDIM
          AMAT1(IOCC,J) = 0.0D0
          AMAT1(J,IOCC) = 0.0D0
        END DO
      END DO
*. -(F(2)-E(0)S(2)) Will be used in the future, reconstruct in AMAT2
      CALL VECSUM(AMAT2,F(1,3),S(1,3),ONE,-E0,NDIM**2)
      ONEM = -1.0D0
      CALL SCALVE(AMAT2,ONEM,NDIM**2)
*. Obtain virtual part of C0
* (F(2)-E(0)S(2))(virt,virt) C0(virt) = -(F(2)-E(0)S(2))C0(occ))(virt)
      CALL MATVCB(AMAT2,C(1,1),VEC1,NDIM,NDIM,0)
      CALL MATVCB(AMAT1,VEC1,VEC2,NDIM,NDIM,0)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Virtual part of C0 '
        CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
      END IF
      CALL COPVEC(VEC2(1+NOCC),C(1+NOCC,1),NVIRT)
C?    WRITE(6,*) ' complete zero order correction vector '
C?    CALL WRTMT_LU(C(1,1),1,NDIM,1,NDIM)
*. And zero it once again ( COnstructed below )
      ZERO = 0.0D0
      CALL SETVEC(C(1+NOCC,1),ZERO,NVIRT)
*
*. And then start the iterations
      DO IORD = 1, NORD
*
* =================================================
* The (IORD-2) wf corrections in the virtual space
* =================================================
*
        IF(IORD.GE.2) THEN
          ZERO = 0.0D0
          CALL SETVEC(VEC2,ZERO,NDIM)
*. Note : Only the occupied part of C(N-2) is included in RHS of
*         expression. The virtual part was carefully zeroed !
          DO K = 0, IORD -2
            CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
            ONE = 1.0D0
            ONEM = -1.0D0
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
          END DO
*
          DO K = 0, IORD -2
            DO L = 0, IORD -K
              CALL MATVCB(S(1,L+1),C(1,K+1),VEC1,NDIM,NDIM,0)
              CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
            END DO
          END DO
*. Multiply with (E(2)-E0S(2))-1
          CALL MATVCB(AMAT1,VEC2,VEC1,NDIM,NDIM,0)
*. And save
          CALL COPVEC(VEC1(1+NOCC),C(1+NOCC,IORD-2+1),NVIRT)
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' occ+virtual part for order = ',IORD-2
            CALL WRTMT_LU(C(1,IORD-2+1),1,NDIM,1,NDIM)
          END IF
*
        END IF
*       ^ End of construction of virtual part of C(VIRT,IORD-2)
*
*  =================
*. Energy correction
*  =================
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      - Sum(I=0,N-1)Sum(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
        EN(IORD) = 0.0D0
        DO I = 1, IORD
          CALL MATVCB(F(1,I+1),C(1,IORD-I+1),VEC1,NDIM,NDIM,0)
          EN(IORD) = EN(IORD)+INPROD(C,VEC1,NDIM)
        END DO
C?      write(6,*) ' First term to En ', EN(IORD)
        DO I = 0,IORD-1
          DO J = 1, IORD-I
            CALL MATVCB(S(1,J+1),C(1,I+1),VEC1,NDIM,NDIM,0)
            EN(IORD) = EN(IORD)-EN(IORD-I-J)*INPROD(C,VEC1,NDIM)
C?          write(6,*) ' second term to EN: I J EN ',I,J,EN(IORD)
          END DO
        END DO
        EN(IORD) = EN(IORD)/C0SC0
        WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ===========================================
*. Wave function corrections, occupied of IORD
*  ===========================================
*
* The occupied part obtained from
*
* C(N) = (F(0)-E(0)S(0))-1 (-Sum(K=0,N-1)F(N-K)C(K)
*                           +Sum(K=0,N-1)Sum(L=0,N-K)E(N-K-L)S(L)C(K))
        ZERO = 0.0D0
        CALL SETVEC(VEC2,ZERO,NDIM)
        DO K = 0, IORD -1
          CALL MATVCB(F(1,IORD-K+1),C(1,K+1),VEC1,NDIM,NDIM,0)
          ONE = 1.0D0
          ONEM = -1.0D0
          CALL VECSUM(VEC2,VEC2,VEC1,ONE,ONEM,NDIM)
        END DO
C?      write(6,*) ' first term to rhs '
C?      CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
*
        DO K = 0, IORD -1
          DO L = 0, IORD -K
            CALL MATVCB(S(1,L+1),C(1,K+1),VEC1,NDIM,NDIM,0)
            CALL VECSUM(VEC2,VEC2,VEC1,ONE,EN(IORD-K-L),NDIM)
C?          write(6,*) ' second term to rhs, K,L = ', K,L
C?          CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
          END DO
        END DO
*. Check overlap with S(0) times zero order state ( should be zero )
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
        CALL MATVCB(S(1,1),C(1,1),VEC1,NDIM,NDIM,0)
        C0SSC0 = INPROD(VEC1,VEC1,NDIM)
        OVLAP = INPROD(VEC1,VEC2,NDIM)
        WRITE(6,*) ' OVLAP = ', OVLAP
        FACTOR = -OVLAP/C0SSC0
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Vector before ortho'
          CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
        END IF
        CALL VECSUM(VEC2,VEC2,VEC1  ,ONE,FACTOR,NDIM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Vector before DIAVC2'
          CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
        END IF
*. Multiply with (F(0)-E(0)S(0))-1
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
*
        JEPZAP = 1
        IF(IORD.EQ.1.AND.JEPZAP.EQ.1) THEN
          WRITE(6,*) ' First order correction zapped '
          ZERO = 0.0D0
          CALL SETVEC(VEC1,ZERO,NOCC)
        END IF
*
        CALL COPVEC(VEC1(1),C(1,IORD+1),NOCC)
*. The virtual part is still not known
        ZERO = 0.0D0
        CALL SETVEC(C(1+NOCC,IORD+1),ZERO,NVIRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Eigenfunction correction ', IORD
          CALL WRTMT_LU(C(1,IORD+1),1,NDIM,1,NDIM)
        END IF
      END DO
*
      WRITE(6,*) ' Energy corrections : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) '   Order        Correction '
      WRITE(6,*) ' ==============================='
      DO IORD = 1, NORD
        WRITE(6,'(1H ,3X,I3,E25.13)')IORD,EN(IORD)
      END DO
*
      ETOT = E0
      DO IORD = 1, NORD
        ETOT = ETOT + EN(IORD)
      END DO
      WRITE(6,'(A,E25.13)') ' Zero-order energy ', E0
      WRITE(6,'(A,E25.13)') ' Sum(K=0,NORD) E(K) ', ETOT
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0DIAG(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,NROOT,
     &           EIGVAL,EIGVEC,SCR,NTESTG,ECORE)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           approximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Obtain the lowest NROOT eigenvectors
*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* NROOT : Number of roots to be obtained
*
* ======
* Output
* ======
* EIGVAL(*) : Converged eigen values
* EIGVEC(IROOT,*) : Complete eigenvector IROOT
*
* Note : The NROOT eigenpairs to be obtained are assumed
*        to be ' slightly ' perturbed eigensolutions
*        of PHP.
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL HPQTVM
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*)
* Output
      DIMENSION EIGVAL(*),EIGVEC(NP1DM+NP2DM+NQDM,*)
* Scratch
      DIMENSION SCR(*)
*.SCR Should atleast be dimensioned ??????
      LOGICAL CONVER
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
* There are two routes :2 Iterative diagonalization of complete matrix
*                       1 complete diagonalizations of partitioned
*                         matrices.
*.
* ========
* Route 1 :
* ========
* The Q-space can be partitioned into the P -space
* to give the effective eigenvalue equation
*
* (PHP - PHQ  (QHQ-E)**-1 QHP ) VP = E VP
*
* This leads to a simple iterative scheme

*
* Note : Only NROOT = 1 tested
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
      CALL QENTER('H0DIA')
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0DIAG speaking '
        WRITE(6,*) ' =============== '
      END IF
*
C?    write(6,*) ' QHQ as delivered '
C?    call WRTMT_LU(QHQ,1,NQDM,1,NQDM)
C?    write(6,*) ' PHP as delivered '
C?    call PRSYM(PHP,NP1DM)
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
*. A bit of memory
*
      I12 = 2
      IF(I12.EQ.1. OR. I12. EQ.3 ) THEN
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. A PQ matrix
      KLPQ = KLFREE
      KLFREE = KLFREE + NPDM * NQDM
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
*
*. Initial eigenvalues
      CALL COPVEC(PHP,SCR(KLPP1),NPDM*(NPDM+1)/2)
      CALL EIGEN_LUCI(SCR(KLPP1),SCR(KLPP2),NPDM,0,1)
*. Extract eigenvalues
      CALL XTRCDI(SCR(KLPP1),EIGVAL,NROOT,1)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Initial set of eigenvalues '
        CALL WRTMT_LU(EIGVAL,1,NROOT,1,NROOT)
      END IF
*. Largest allowed number of iterations
      MAXIT = 5
      DO 1000 IROOT = 1, NROOT
        CONVER = .FALSE.
        EINI = EIGVAL(IROOT)
        DO 900 ITER = 1, MAXIT
          IF(NTEST.GE.2) WRITE(6,*) ' Info from iteration ', ITER
*. Current eigenvalue and eigenvector
          E = EIGVAL(IROOT)
          CALL COPVEC(SCR(KLPP2+(IROOT-1)*NPDM),SCR(KLV1),NPDM)
* ==============================
* HPP - PHQ (QHQ - E) **-1 * QHP
* ==============================
*. QHP in KLPQ
C         TRPMAT(XIN,NROW,NCOL,XOUT)
          CALL TRPMAT(PHQ,NP1DM,NQDM,SCR(KLPQ))
*.Multiply with (QHQ-E)**-1
          DO 30 IP1 = 1, NP1DM
C           DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
            IOFF = KLPQ + (IP1-1)*NQDM
            CALL DIAVC3_LUCI(SCR(IOFF),SCR(IOFF),QHQ,-E,NQDM,XDUMMY)
   30     CONTINUE
*. Multiply with PHQ
          CALL MATML4(SCR(KLPP1),PHQ,SCR(KLPQ),NP1DM,NP1DM,
     &                NP1DM,NQDM,NQDM,NP1DM,0)
*.
C?      write(6,*) ' PHQ (QHQ-E)-1 QHP matrix '
C?      CALL WRTMT_LU(SCR(KLPP1),NP1DM,NP1DM,NP1DM,NP1DM)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          SIGNTP = 1.0
          CALL TRIPAK_LUCI(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM,SIGNTP)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL PRSYM(SCR(KLPP1),NPDM)
          END IF
*.Diagonalize
          CALL EIGEN_LUCI(SCR(KLPP1),SCR(KLPP2),NPDM,0,1)
*. Extract eigenvalues
           EIGVAL(IROOT) = SCR(KLPP1-1+IROOT*(IROOT+1)/2)
           IF(NTEST.GE.2)
     &     WRITE(6,*) ' Eigenvalue ', EIGVAL(IROOT)
           IF(NTEST.GE.10) THEN
             WRITE(6,*) ' P-space eigenvector '
             CALL WRTMT_LU(SCR(KLPP2+(IROOT-1)*NPDM),
     &            1,NPDM,1,NPDM)
           END IF
*. Check for convergence
           EVALDF = ABS(E-EIGVAL(IROOT))
           EVECOV = SQRT(INPROD(SCR(KLV1),SCR(KLPP2+(IROOT-1)*NPDM),
     &                          NPDM) )
           IF(EVALDF.LT.1.0D-7.AND.EVECOV.GT.0.999D0) CONVER = .TRUE.
           IF(CONVER) GOTO 901
  900    CONTINUE
  901    CONTINUE
*. P-part of eigenvector
         CALL COPVEC(SCR(KLPP2+(IROOT-1)*NPDM),EIGVEC(1,IROOT),
     &               NPDM)
*. Obtain Q part of eigenvector
*.    -(QHQ-E)**-1 HQP XP

         CALL MATML4(SCR(KLV1),PHQ,EIGVEC(1,IROOT),
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
*
C             DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3_LUCI(EIGVEC(NPDM+1,IROOT),SCR(KLV1),QHQ,
     &               -EIGVAL(IROOT),NQDM,XDUMMY)
         CALL SCALVE(EIGVEC(NPDM+1,IROOT),-1.0D0,NQDM)
*. Normalize
         XNORM = INPROD(EIGVEC(1,IROOT),EIGVEC(1,IROOT),NPQDM)
         SCALE = 1.0D0/SQRT(XNORM)
         CALL SCALVE(EIGVEC(1,IROOT),SCALE,NPQDM)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Initial and final eigenvalue ',
     &     EINI,EIGVAL(IROOT)
           WRITE(6,*)
     &   ' Part of eigenvector in Q space',
     &     SQRT(ABS(XNORM-1.0D0)/XNORM)
C?         WRITE(6,*) ' Eigenvector in PQ space '
C?         CALL WRTMT_LU(EIGVEC(1,IROOT),1,NPQDM,1,NPQDM)
         END IF
*
 1000 CONTINUE
       END IF
       IF( I12. EQ.2 .OR. I12 .EQ.3 ) THEN
*
*. Iterative scheme
*
*
*. Initial eigenvalues and eigenvectors
*
       KLPP1 = 1
       KLFREE = KLPP1 + NP1DM*(NP1DM+1)/2
       KLPP2 = KLFREE
       KLFREE = KLPP2 + NP1DM*NP1DM
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LU4 = 37
        LU5 = 38
        LUDIA = 39
        CALL COPVEC(PHP,SCR(KLPP1),NP1DM*(NP1DM+1)/2)
        CALL EIGEN_LUCI(SCR(KLPP1),SCR(KLPP2),NP1DM,0,1)
*. Extract eigenvalues and eigenvectors on LU1
        CALL XTRCDI(SCR(KLPP1),EIGVAL,NROOT,1)
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' Initial set of eigenvalues '
          CALL WRTMT_LU(EIGVAL,1,NROOT,1,NROOT)
        END IF
*. Eigenvectors on LU1
        CALL SETVEC(EIGVEC(1,1),0.0D0,NPQDM)
        CALL REWINE(LU1,-1)
        DO 510 IROOT = 1, NROOT
          CALL COPVEC(SCR(KLPP2+(IROOT-1)*NP1DM),EIGVEC(1,1),NP1DM)
          CALL TODSC_LUCI(EIGVEC(1,1),NPQDM,-1,LU1)
  510   CONTINUE
*
*. Iterations
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LUDIA)
*. Davidson CI diagonalization
        MINST = 1
        NBLK = 1
        INICI = -1
        MAXCIT = 15
        IPRTCI = MAX(NTEST-2,0)
        MXVCCI = 3 * NROOT
        SHIFT = 0.0D0
        IDIAG = 1
        ICISTR = 1
        IDIAG = 1
        THRES= 1.0D-8
* FIXME !!!
C  Only 33 arguments in this call, 70 are required.
C
c       CALL CIEIG5(HPQTVM,INICI,EIGVAL,SCR(KLV1),SCR(KLV2),
c    &            MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,NPQDM ,
c    &            NBLK,NROOT,MXVCCI,MAXCIT,LU1,IPRTCI,
c    &            DUMMY,0,DUMMY,IDUMMY,
c    &            0,0,0,DUMMY,ECORE,ICISTR,NPQDM,IDIAG,DUMMY,THRES)
C       CALL CIEIG5(HPQTVM,INICI,ENOT,SCR(KLV1),SCR(KLV2),
C    &            MINST,LUDIA,LU1,LU2,LU3,LU4,LU5,NPQDM ,
C    &            NBLK,NROOT,MXVCCI,MAXCIT,LU1,IPRTCI,
C    &            DUMMY,0,DUMMY,IDUMMY,
C    &            0,0,0,DUMMY,0.0D0,1,NPQDM,1,DUMMY)
        CALL REWINE(LU1,-1)
        DO 1286 JROOT = 1, NROOT
         CALL FRMDSC_LUCI(EIGVEC(1,JROOT),NPQDM,-1,LU1,IMZERO,IAMPACK)
 1286   CONTINUE
      END IF
*
      CALL QEXIT('H0DIA')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0INTSPC(IH0SPC,NPTSPC,IOCPTSPC,NOCTPA,NOCTPB,
     &                    IOCA,IOCB,NGAS,MXPNGAS,INTH0SPC,NELFTP)
*
* Set up INTH0SPC : Division of CI space, so only determinants
*                   belonging to the same space  have nonvanishing
*                   matrix elements of H0
*
* =====
* Input
* =====
*
* IH0SPC : ne. 0 : Interacting subspaces have been defined
*          .eq.0 : Interacting subspaces not defined, let
*                  evrything interact
* NPTSPC : Number of subspaces defined
* IOCPTSPC : Allowed occumulated occupation of each subspace
* NOCTPA :  Number of alpha occupation types
* NOCTPB : Number of beta occupation types
* IOCA : Occupation  of alpha string
* IOCB : Occupation  of beta string
*
* Jeppe Olsen, January 1996
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      DIMENSION IOCPTSPC(2,MXPNGAS,*)
      DIMENSION IOCA(MXPNGAS,*),IOCB(MXPNGAS,*)
      DIMENSION NELFTP(*)
*. Output
      DIMENSION INTH0SPC(NOCTPA,NOCTPB)
*
      IF(IH0SPC.EQ.0) THEN
*. All interactions allowed
        IONE = 1
        CALL ISETVC(INTH0SPC,IONE,NOCTPA*NOCTPB)
      ELSE
*. Explicit construction of matrix giving partitioning of
*  subspaces
        IZERO = 0
        CALL ISETVC(INTH0SPC,IZERO,NOCTPA*NOCTPB)
*
        DO ISPC = 1, NPTSPC
          DO IATP = 1, NOCTPA
            DO IBTP = 1, NOCTPB
              IAMOKAY = 1
              IEL = 0
C?            WRITE(6,*) ' ISPC IATP IBTP ', ISPC,IATP,IBTP
              DO IGAS = 1, NGAS
               IEL = IEL
     &             + NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
C?             WRITE(6,*) ' IGAS IEL ', IGAS,IEL
C?             WRITE(6,*)
C?   &          ' Limits :',IOCPTSPC(1,IGAS,ISPC),IOCPTSPC(2,IGAS,ISPC)
               IF(IEL.LT.IOCPTSPC(1,IGAS,ISPC).OR.
     &            IEL.GT.IOCPTSPC(2,IGAS,ISPC)    ) IAMOKAY = 0
              END DO
C?            WRITE(6,*) ' IAMOKAY = ', IAMOKAY
*. Allowed
              IF(IAMOKAY.EQ.1.AND.INTH0SPC(IATP,IBTP).EQ.0) THEN
                INTH0SPC(IATP,IBTP) = ISPC
              END IF
            END DO
          END DO
        END DO
      END IF
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Output from  H0INTSPC '
        WRITE(6,*) ' ======================'
        WRITE(6,*)
        WRITE(6,*) ' Output from H0INTSPC '
        WRITE(6,*)
        CALL IWRTMA(INTH0SPC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
     &           X,RHS,S,SCR,NTESTG)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Solve the set of equations
*
*     ( H0+S ) X = RHS

*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* Modified to allow solution by conjugate gradient, March 1993
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* RHS   : Right hand side of equations
*
* ======
* Output
* ======
* X : solution to linear equations
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONVER
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*),RHS(*)
* Output
      DIMENSION X(*)
* Scratch
      DIMENSION SCR(*), ERROR(20+1)
*.SCR Should atleast be dimensioned 2 *(NP1DM+NP2DM)** 2 + 2 NPQDM
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
      EXTERNAL HPQTVM
*.
* The Q-space can be partitioned into the P -space
* to give the effective linear equation
*
* (PHP+S - PHQ  (QHQ+S)**-1 QHP ) XP = RHSP - HPQ(QHQ+S)-1 RHSQ
*
* This leads to a simple iterative scheme
*
      CALL QENTER('H0LNS')
      NTESTL =  00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0LNSL speaking '
        WRITE(6,*) ' =============== '
      END IF
*
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
      IROUTE = 2
*
      IF( IROUTE.EQ.1. OR. IROUTE. EQ.3 ) THEN
*. Solve by partitioning theory
*. A bit of memory
*
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
* =========================
*  RHSP - HPQ(QHQ+S)-1 RHSQ
* =========================
*          DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
      CALL DIAVC3_LUCI(SCR(KLV1),RHS(1+NPDM),QHQ,S,NQDM,XDUMMY)
      CALL MATML4(SCR(KLV2),PHQ,SCR(KLV1),NP1DM,1,NP1DM,NQDM,NQDM,1,0)
      CALL VECSUM(SCR(KLV1),RHS,SCR(KLV2),1.0D0,-1.0D0,NP1DM)
      CALL COPVEC(RHS(1+NP1DM),SCR(KLV1+NP1DM),NP2DM)
* ===============================
* (PHP+S - PHQ  (QHQ+S)**-1 QHP )
* ===============================
C          XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
      CALL XDIXT2(SCR(KLPP1),PHQ,QHQ,NP1DM,NQDM,S,SCR(KLV2))
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGNTP)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          SIGNTP = 1.0
          CALL TRIPAK_LUCI(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM,SIGNTP)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
C                ADDDIA(A,FACTOR,NDIM,IPACK)
          CALL ADDDIA(SCR(KLPP1),S,NPDM,1)
*. Pack to full matrix
          CALL TRIPAK_LUCI(SCR(KLPP2),SCR(KLPP1),2,NPDM,NPDM,SIGNTP)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL WRTMT_LU(SCR(KLPP2),NPDM,NPDM,NPDM,NPDM)
          END IF
*.Solve p equations by inverting and multiplying
C               INVMAT(A,B,MATDIM,NDIM)
           CALL INVMAT(SCR(KLPP2),SCR(KLPP1),NPDM,NPDM)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
           CALL MATVCB(SCR(KLPP2),SCR(KLV1),X,NPDM,NPDM,0)
*. q part of solution
* ==================================
* XQ = (QHQ+SHIFT)**-1 (RHS Q - QHP XP)
* ==================================

         CALL MATML4(SCR(KLV1),PHQ,X,
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
         CALL VECSUM(SCR(KLV2),RHS(NPDM+1),SCR(KLV1),1.0D0,
     &               -1.0D0,NQDM)
*
C             DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3_LUCI(X(NPDM+1),SCR(KLV2),QHQ,
     &               S,NQDM,XDUMMY)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(X,1,NPQDM,1,NPQDM)
         END IF
      END IF
*
      IF (IROUTE. EQ. 2 .OR. IROUTE .EQ. 3 ) THEN
*. Use preconditioned conjugate gradient
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LUDIA = 37
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LUDIA)
*. Initial Guess
        CALL REWINE(LU1,-1)
        CALL SETVEC(SCR(KLV1),0.0D0,NPQDM)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LU1)
*. Right hand side
        CALL REWINE(LU2,-1)
        CALL TODSC_LUCI(RHS,NPQDM,-1,LU2)
*
        MAXIT = 20
        CONVER = .FALSE.
        TEST = 1.0D-9 * SQRT(INPROD(RHS,RHS,NPQDM))
        SHIFT = S
        ILNPRT = MAX(NTEST-10,0)
        CALL MINGCG(HPQTVM,LU1,LU2,LU3,LUDIA,SCR(KLV1),SCR(KLV2),
     &              MAXIT,CONVER,TEST,S,ERROR,NPQDM,0,ILNPRT)
        CALL REWINE(LU1,-1)
        CALL FRMDSC_LUCI(SCR(KLV1),NPQDM,-1,LU1,IMZERO,IAMPACK)
        CALL COPVEC(SCR(KLV1),X,NPQDM)
*
         IF(NTEST.GE.50) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(X,1,NPQDM,1,NPQDM)
         END IF
*
      END IF
*
      CALL QEXIT('H0LNS')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
* Disc version
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  LUOUT       LUDIA        LUIN
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal on file LUDIA
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements on LUDIA corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , September 1993
*
*
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File Containing diagonal of H0
* LUIN  : File Containing input vector
* LBLK : Defines format of files
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
COLD  REAL * 8  INPROD
*
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
!     work is h0scr in calling routine
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.10)
     & write(6,*) ' H0M1TV , NPQDM = ', NPQDM
*
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
*. Obtain subspace components of input vector
C            GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,NTESTG)
        IZERO = 0
        CALL GATVCD(LUIN,LBLK,NPQDM,IPNTR,WORK(KLV1),VEC1,
     &              IZERO)
*. Solve linear equations in subspace
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
      END IF
*
*. Calculate inverse diagonal and scatter results from subspace,
*. Write to file LUOUT
C     DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
C    &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
      CALL DMTVDS(VEC1,VEC2,LUDIA,LUIN,LUOUT,SHIFT,1,1,
     &            IPNTR,WORK(KLV2),NPQDM,LBLK,XH0PSX)
*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TD '
        WRITE(6,*) ' ========================= '
*. Note : works only if result vector is first file on LUOUT
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
        WRITE(6,*) ' Overlap between input and output vector',
     &               XH0PSX
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0M1TV(DIAG,VECIN,VECUT,NVAR,NPQDM,IPNTR,
     &                  H0,SHIFT,WORK,XH0PSX,
     &                  NP1,NP2,NQ,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal Diag
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements in DIAG corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , May 1990

*
*
* =====
* Input
* =====
* DIAG : Diagonal of matrix
* VECIN : Input vector
* NVAR : Dimension of full matrix
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
* ======
* Output
* ======
* VECUT : Output vector (you guessed ?? ), can occupy same space
*         as VECIN or DIAG
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      REAL * 8  INPROD
*
      DIMENSION DIAG(*),VECIN(*),VECUT(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
!               work is h0scr in calling routine
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
C?    write(6,*) ' H0M1TV , NPQDM = ', NPQDM
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
        CALL GATVEC(WORK(KLV1),VECIN,IPNTR,NPQDM)
* X(T)(DIAG+SHIFT)-1 X in subspace, for later subtraction
        CALL GATVEC(WORK(KLV2),DIAG,IPNTR,NPQDM)
        CALL DIAVC3_LUCI(WORK(KLV2),WORK(KLV1),
     &       WORK(KLV2),SHIFT,NPQDM,X1)
       ELSE
         X1 = 0.0D0
       END IF
*
      CALL DIAVC3_LUCI(VECUT,VECIN,DIAG,SHIFT,NVAR,X2)
*
      IF(NPQDM .NE. 0 ) THEN
C                H0LNSL(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
C    &           X,RHS,S,SCR)
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
C?     write(6,*) ' KLPHP KLPHQ KLQHQ ',KLPHP,KLPHQ,KLQHQ
*
         CALL H0LNSL(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &               WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &               NTEST )
         X3 = INPROD(WORK(KLV1),WORK(KLV2),NPQDM)
         CALL SCAVEC(VECUT,WORK(KLV2),IPNTR,NPQDM)
      ELSE
         X3 = 0.0D0
      END IF
      XH0PSX  = X2 - X1 + X3
C?    write(6,*) ' XH0PSX x1 x2 x3 ', XH0PSX,X1,X2,X3


*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TV '
        WRITE(6,*) ' ========================= '
        CALL WRTMT_LU(VECUT,1,NVAR,1,NVAR)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0MAT(INTSPC,NROOT,ONEBOD,H0,SBEVC,SBEVL,ISBDET,
     &                 ISBIA,ISBIB,ISBCNF,
     &                 LUHDIA,LBLK,
     &                 MXP1,MXP2,MXQ,
     &                 MP1CSF,MP2CSF,MQCSF,NOCOB,
     &                 NPRCIV,NOCSF,IREFSM,IPRT,IPROCC,
     &                 VEC1,VEC2,H0SCR,IDC,PSSIGN,ECORE)
* Obtain preconditioner space corresponding to internalt space INTSPC
* Obtain Hamiltonian matrices correponding to this subspacw
* Obtain the first Nroot eigensolutions.
*
* Construct Preconditioner blocks of Hamilton matrix
*
* ======
*.Output
* ======
*
* CSF : NP1CSF,NP2CSF,NQCSF : Number of CSF's in the 3 primary subspaces
*
* NPRCIV : Number of parameters in preconditioner space
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*8 KIOCOC
!               for addressing of WORK
#include "mxpdim.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "stinf.inc"
#include "spinfo.inc"
#include "csfbas.inc"
#include "wrkspc.inc"
      DIMENSION ONEBOD(*)
*
      NTEST = 0
      NTEST = MAX(IPRT,NTEST)
*
      MXSBDT = MXP1+MXP2+MXQ
      IF(NTEST.NE.0.AND.MXSBDT.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ========================================== '
        WRITE(6,*) '    Information about CI preconditioner '
        WRITE(6,*) ' ========================================== '
        WRITE(6,*)
       END IF
       WRITE(6,*) ' H0MAT : ecore ', ECORE
* Info on actual internal subspace
      IATP = IASTFI(INTSPC)
      IBTP = IBSTFI(INTSPC)
      MNR1 = MNR1IC(INTSPC)
      MXR1 = MXR1IC(INTSPC)
      MNR3 = MNR3IC(INTSPC)
      MXR3 = MXR3IC(INTSPC)
      NAEL = NAELCI(INTSPC)
      NBEL = NBELCI(INTSPC)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Allowed combination of alpha and beta strings
      CALL MEMMAN(KIOCOC,NOCTPA*NOCTPB,'ADDL  ',2,'IOCOC ')
      CALL IAIBCM_GAS(MNR1,MXR3,NOCTPA,NOCTPB,WORK(KEL1(IATP)),
     &            WORK(KEL3(IATP)),WORK(KEL1(IBTP)),WORK(KEL3(IBTP)),
     &            WORK(KIOCOC),NTEST)
*
      IF(IDC.EQ.1) THEN
        ICOMBI = 0
        PSIGN = 0.0D0
      ELSE
        PSIGN = PSSIGN
        ICOMBI = 1
      END IF
*
      IF( NOCSF .NE. 0) THEN
*.Combinations expansion, PQ preconditioner
*
        IHAMSM = 1
        IWAY = 1
* strings are unsigned
        ISTSGN = 0
        CALL H0SD(LUHDIA,LBLK,VEC1,IWAY,NSBDET,NAEL,NBEL,
     &            ISMOST(1,IREFSM),WORK(KIOCOC),
     &            IHAMSM,H0,ONEBOD,NOCOB,0,
     &            ECORE,ICOMBI,PSIGN,NPRCIV,SBEVC,
     &            SBEVL,1,NCIVAR,ISBDET,ISBIA,ISBIB,NROOT,
     &            MXP1,MXP2,MXQ,
     &            MP1CSF,MP2CSF,MQCSF,
     &            WORK(KOCSTR(IATP)),WORK(KOCSTR(IBTP)),
     &            ISTSGN,IDUMMY,IDUMMY,
     &            INTSPC,IPRT,IPROCC)
      END IF
*
      IF(IPRT.NE.0 .AND. NPRCIV.LT.NSBDET) THEN
        NSBDET = NPRCIV
        WRITE(6,*)
     &  '  Number of elements in primary space REDUCED to . ',NSBDET
      END IF
*
      IF(NTEST.NE.0.AND.MXSBDT.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================================== '
        WRITE(6,*) '    End of information about CI preconditioner '
        WRITE(6,*) ' ============================================== '
        WRITE(6,*)
       END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0SD(LUDIA,LBLK,CIDIA,IWAY,MXPRDT,NAEL,NBEL,
     &                ISMOST,IOCOC,
     &                IHSYM,HAMIL,ONEBOD,NORB,
     &                NINOB,ECORE,ICOMBI,PSIGN,NPRDET,EIGVEC,
     &                EIGVAL,IDODGN,NDET,IDET,IA,IB,NROOT,
     &                MXP1DM,MXP2DM,MXQDM,NP1DM,NP2DM,NQDM,
     &                IASTR,IBSTR,ISTSGN,IASGN,IBSGN,
     &                INTSPC,IPRT,IPROCC)
*
* Construct an P1 P2 Q preconditioner matrix in the SD basis
* and/or find and print the lowest IPROCC elements
*
*.. Subspace :
* ============
*
*     IWAY = 1 : choose the lowest values of
*                the CI diagonal. The number of combinations is
*                chosen so that no degenerate levels are splitted.
*                The number of combinations used, NPRDET, can
*                therefore be lower than mxsbdt.
*     IWAY = 2 : Choose the first NPRDET combinations,
*                stupid, but convenient for testing
*
* IDET contains adresses of elements chosen
*
* IDODGN .GT. 0 : DIAGONALIZE CONSTRUCTED HAMILTON MATRIX.
*                 EIGVL CONTAINS EIGENVALUES ON RETURN
*                 EIGVC CONTAINS EIGENVECTORS(COLUMNS) ON RETURN
*
* Obtained from EXPHAM, January 1993
*
* Put into LUCIA form June 1993
* Combinations enabled September 1993
*
*. Core energy : not added to matrix, added to resulting eigenvalues
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*8 KL1, KL2, KL3, KSCR, KLIAST, KLIBST, KFREE
!               for addressing of WORK
*. Input
      DIMENSION ISMOST(*)
      DIMENSION ONEBOD(*),CIDIA(*)
      INTEGER IASTR(NAEL,*),IBSTR(NBEL,*)
*. Output
      DIMENSION HAMIL(*)
      DIMENSION IDET(MXPRDT), EIGVAL(NROOT), EIGVEC(NROOT * MXPRDT )
      DIMENSION IA(MXPRDT),IB(MXPRDT)
*
#include "mxpdim.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*
      CALL QENTER('EXPHAM')
      NTEST = 0
      NTEST = MAX(NTEST,IPRT)
      MXDM1 = MXP1DM + MXP2DM + MXQDM
      MXDM = MAX(MXDM1,IPROCC)
      IF (NTEST .GE. 5 ) THEN
        WRITE(6,*) ' MXDM ', MXDM
        write(6,*) ' MXP1DM MXP2DM MXQDM ', MXP1DM, MXP2DM,MXQDM
        write(6,*) ' LUDIA and LBLK ', LUDIA,LBLK
        write(6,*) ' INTSPC and IPRT ', INTSPC,IPRT
        write(6,*) ' NROOT ', NROOT
        write(6,*) ' IPROCC ', IPROCC
      END IF
*
* ====================
* 0 : Select subspace
* ====================
*
      IF( IWAY .EQ. 1) THEN
*.      Find number of combinations less or equal to MXDM
*       that does not separate degenerate pairs .
*. Used Scratch space : 6 * MXDM ( a bit extravagant )
        CALL MEMMAN(KL1,MAX(NAEL,3*MXDM+1),'ADDL  ',1,'KL1   ')
        CALL MEMMAN(KL2,MAX(NBEL,2*MXDM+1),'ADDL  ',2,'KL2   ')
        CALL MEMMAN(KL3,2*MXDM+1,'ADDL  ',2,'KL3   ')
        CALL FNDMND(LUDIA,LBLK,CIDIA,MXDM,NPRDET,WORK(KL1),
     &              WORK(KL2),IDET,WORK(KL3),NTEST)
      ELSE IF ( IWAY .EQ. 2 ) THEN
CTF     Only for debugging purposes. Fixed argument list.
        IZERO = 0
        CALL ISETVC(IDET,IZERO,NPRDET)
      END IF
*. Check for degeneracies on the boundaries between P1, P2 and Q space
*. P1 - P2
      IF( MXP1DM .GT. 0 ) THEN
        IF(MXP1DM.GE.NPRDET ) THEN
          NP1DM = NPRDET
        ELSE
          IIDET = MXP1DM
 101      CONTINUE
          IF( ABS(WORK(KL3-1+IIDET+1)-WORK(KL3-1+IIDET))
     &       .LE. 0.000001) THEN
             IIDET = IIDET - 1
             GOTO 101
          END IF
          NP1DM = IIDET
        END IF
      ELSE
        NP1DM = 0
      END IF
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 )
     & WRITE(6,*) ' Actual dimension of P1 Space ', NP1DM
*. P2 - Q space
      IF(MXP2DM.GT.0) THEN
        IF(MXP1DM+MXP2DM.GE.NPRDET) THEN
          NP2DM = NPRDET - NP1DM
        ELSE
          IIDET = MXP1DM + MXP2DM
 102      CONTINUE
          IF( ABS(WORK(KL3-1+IIDET+1)-WORK(KL3-1+IIDET))
     &       .LE. 0.0000001) THEN
             IIDET = IIDET - 1
             GOTO 102
          END IF
          NP2DM = IIDET - NP1DM
        END IF
      ELSE
        NP2DM = 0
      END IF
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 )
     & WRITE(6,*) ' Actual dimension of P2 Space ', NP2DM
*. Q space
      IF(MXQDM.NE.0) THEN
        NQDM = NPRDET - NP1DM - NP2DM
      ELSE
        NQDM = 0
      END IF
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 )
     & WRITE(6,*) ' Actual dimension of Q Space ', NQDM
*
      NPDM = NP1DM + NP2DM
      NPRDET = NP1DM + NP2DM + NQDM
*. Copy over to MX numbers ( not nice but .. )
      MXP1DM = NP1DM
      MXP2DM = NP2DM
      MXQDM = NQDM
*
      IF( NTEST .GE. 2 .AND. MXDM1 .NE. 0 )
     & WRITE(6,*)
     & ' Total number of combinations in subspace : ', NPRDET
      IF(NTEST .GE. 10 ) THEN
        WRITE(6,*) ' IDET IN EXPHAM '
        CALL IWRTMA(IDET,1,NPRDET,1,NPRDET)

      END IF
*
* =============================================================
* Alpha and beta strings corresponding to selected combinations
* =============================================================
*
*
*. Convert determinant numbers to string numbers
*
*. P dets
      CALL STRFDT_LUCI(INTSPC,ISMOST,IOCOC,NPDM,IDET,IA,IB,ICOMBI)
*. Q dets
      IF(NQDM.NE.0)
     &CALL STRFDT_LUCI(INTSPC,ISMOST,IOCOC,NQDM,
     &            IDET(1+NPDM),IA(1+NPDM),IB(1+NPDM),ICOMBI)
*. and remaining to be printed
      IF(NPDM+NQDM.LT.IPROCC) THEN
        LREST = IPROCC - NPDM - NQDM
        CALL STRFDT_LUCI(INTSPC,ISMOST,IOCOC,LREST,
     &       IDET(1+NPDM+NQDM),IA(1+NPDM+NQDM),
     &       IB(1+NPDM+NQDM),ICOMBI)
       END IF
*
      IF(IPROCC.NE.0) THEN
*. Print occupation of lowest IPROCC lowest SD 's
        WRITE(6,*)
        WRITE(6,'(A)')
     &  ' ========================================================='
        WRITE(6,'(A,I4,A)')
     &  ' Occupation and energy of lowest ', IPROCC, ' combinations'
        WRITE(6,'(A)')
     &  ' ========================================================='
        WRITE(6,*)
        DO I = 1, IPROCC
*
          WRITE(6,'(A,I8,A,F18.10)')
     &    '  Energy  of combination ',IDET(I),' is ',
     &       WORK(KL3-1+I)+ECORE
          WRITE(6,'(A)')
     &    '  Corresponding alpha - and beta string '
            WRITE(6,'(4X,10I4)')
     &      (IASTR(IEL,IA(I)),IEL = 1, NAEL )
            WRITE(6,'(4X,10I4)')
     &      (IBSTR(IEL,IB(I)),IEL = 1, NBEL )
C         ELSE
C           WRITE(6,'(4X,10(1X,A6))')
C    &      (IOBLAB(IASTR(IEL,IA)),IEL = 1, NAEL )
C           WRITE(6,'(4X,10(1X,A6))')
C    &      (IOBLAB(IBSTR(IEL,IB)),IEL = 1, NBEL )
C         END IF
        END DO
      END IF
      MXDM = MXDM1
      IF(MXDM.EQ.0) GOTO 9999
*
* ==================================
* Obtain Hamiltonian matrix elements
* ==================================
*
*. Pointers
      KLPHP = 1
      KLPHQ = KLPHP + NPDM*(NPDM+1)/2
      KLQHQ = KLPHQ + NP1DM * NQDM
*. PHP Hamiltonian
      CALL QENTER('DIHDJ')
*
* Scratch space for DIHDJ
      LSCR = 4 * NORB
      CALL MEMMAN(KSCR,LSCR,'ADDL  ',1,' H0SCR')
      CALL MEMMAN(KLIAST,NAEL,'ADDL  ',1,'LIA   ')
      CALL MEMMAN(KLIBST,NBEL,'ADDL  ',1,'LIB   ')
*
      ISTSGN = 0
      ECOREP = 0.0D0
      CALL DIHDJ_LUCI(IA,IB,NPDM,IA,IB,
     &           NPDM,NAEL,NBEL,WORK(KSCR),LSCR,
     &           NORB,ONEBOD,HAMIL(KLPHP),
     &           1,NINOB,ECOREP,ICOMBI,PSIGN,IASTR,IBSTR,IASTR,
     &           IBSTR,ISTSGN,IASGN,IBSGN,IASGN,IBSGN,WORK(KLIAST),
     &           WORK(KLIBST),NDIF0,NDIF1,NDIF2,NTEST)
*. PHQ Hamiltonian
      IF(NQDM.NE.0)
     &CALL DIHDJ_LUCI(IA,IB,NP1DM,IA(1+NPDM),IB(1+NPDM),NQDM,
     &           NAEL,NBEL,WORK(KSCR),LSCR,
     &           NORB,ONEBOD,HAMIL(KLPHQ),
     &           0,NINOB,ECOREP,ICOMBI,PSIGN,IASTR,IBSTR,IASTR,
     &           IBSTR,ISTSGN,IASGN,IBSGN,IASGN,IBSGN,
     &           WORK(KLIAST),
     &           WORK(KLIBST),NDIF0,NDIF1,NDIF2,NTEST)
      CALL QEXIT('DIHDJ')
*. QHQ Hamiltonian
      IF(LUDIA.LE.0) THEN
        DO 0607 IIDET = NPDM + 1, NPDM+NQDM
          HAMIL(KLQHQ + IIDET - NPDM-1) =
     &    CIDIA(IDET(IIDET))
 0607   CONTINUE
      ELSE IF (LUDIA.GT. 0 ) THEN
C            GATVCD(LU   ,LBLK,NGAT,IGAT,XGAT,SEGMNT,IPRT)
        IZERO = 0
        CALL GATVCD(LUDIA,LBLK,NQDM,IDET(NPDM+1),
     &              HAMIL(KLQHQ),CIDIA,IZERO)
      END IF
*
      IF(NTEST .GE. 20 ) THEN
        IF(NQDM .NE. 0 ) THEN
          WRITE(6,*) ' PHP, PHQ and QHQ parts of H0 '
        END IF
        WRITE(6,*) ' PHP '
        WRITE(6,*) ' === '
        CALL PRSYM(HAMIL(KLPHP),NPDM)
        IF(NQDM.NE.0) THEN
          WRITE(6,*) ' PHQ '
          WRITE(6,*) ' === '
          CALL WRTMT_LU(HAMIL(KLPHQ),NP1DM,NQDM,NP1DM,NQDM)
          WRITE(6,*) ' QHQ '
          WRITE(6,*) ' === '
          CALL WRTMT_LU(HAMIL(KLQHQ),1,NQDM,1,NQDM)
        END IF
      END IF
*
* =================================
* Diagonalize to obtain lowest roots
* =================================
*
      IF(IDODGN .GT. 0 ) THEN
        CALL MEMMAN(KFREE,IDUMMY,'FREE  ',IDUMMY,'CDUMMY')
        CALL H0DIAG(HAMIL(KLPHP),HAMIL(KLPHQ),HAMIL(KLQHQ),
     &       NP1DM,NP2DM,NQDM,NROOT,EIGVAL,EIGVEC,WORK(KFREE),
     &       NTEST,ECORE )
        IF( NTEST .GE. 2 ) THEN
          WRITE(6,*) ' Eigenvalues of subspace Hamiltonian '
C         CALL WRTMT_LU(EIGVAL,1,NROOT,1,NROOT)
          WRITE(6,'(5F18.10)') (EIGVAL(IROOT),IROOT=1,NROOT)
        END IF
        IF( NTEST .GE. 50 ) THEN
          WRITE(6,*) ' Eigenvectors of subspace Hamiltonian '
          CALL WRTMT_LU(EIGVEC,NPRDET,NROOT,NPRDET,NROOT)
*
          WRITE(6,*) ' Subspace determinants '
          CALL IWRTMA(IDET,1,NPRDET,1,NPRDET)
        END IF
      END IF
*
* ===============================
* Analyze the Root approximations
* ===============================
*
C     IF( IPRT .GT. 0 ) THEN
C     DO 1869 IROOT = 1, NROOT
C       IOFF = (IROOT-1)*NPRDET + 1
C       CALL SETVEC(CIDIA,0.0D0,NDET)
C       CALL SCAVEC(CIDIA,EIGVEC(IOFF),IDET,NPRDET)
C       WRITE(6,*)
C       WRITE(6,'(A,I3)') '  Information about root ... ',IROOT
C       WRITE(6,'(A)')    '  ******************************'
C       WRITE(6,*)
C       WRITE(6,'(A,E15.8)') '   Energy .... ',EIGVAL(IROOT)
C       CUTOFF = 0.1D0
CTOBE   CALL ANACI(CIDIA,ICSYM,CUTOFF,100)
C1869 CONTINUE
C     END IF
*
 9999 CONTINUE
*
      CALL QEXIT('EXPHAM')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE HINTV(LURHS,LUX,SHIFT,SHIFT_DIA,VECIN,VECOUT,
     &                LBLK,LUPROJ,LUPROJ2)
*
* Solve  (H+Shift)X = RHS
*
* Where H is matrix rep of operator defined by /OPER/ in space defined by
* /CANDS/
*
* If ICISTR.EQ.1 VECIN contains RHS, else RHS is assumed  on LURHS
* Output : solution is on LUX
*
* Jeppe Olsen, Winter of 1996
*
* Jan. 98 : SHIFT_DIA added
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLCLBT, KLCLEBT, KLCI1BT, KLCIBT, KLCIOIO, KLCBLTP
!               for addressing of WORK
      REAL*8 INPROD , INPRDD
      LOGICAL CONVER
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
#include "mxpdim.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
#include "wrkspc.inc"
#include "clunit.inc"
* SCRATCH files used : LUSC3,LUSC34,LUSC35,LUSC37
*. These are not used in SIMPRT !
#include "oper.inc"
#include "crun.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "stinf.inc"
#include "csm.inc"
*. Max number of iterations is picked from MAXCIT in crun
      EXTERNAL H0TVM
      DIMENSION ERROR(100)

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'HINTV ')
*
      NTEST = 0
*
*
* 1 : Construct diagonal  on LUSC3
* ================================
*
*. Use type of H0 as type of zero order operator  in all spaces
       IF(IAPR.NE.0) THEN
*. Not exact Hamiltonian
         IF(IH0INSPC(1).EQ.1.OR.IH0INSPC(1).EQ.3) THEN
*. Mp operator
           I12 = 1
           CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
           CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
         ELSE IF (IH0INSPC(1).EQ.2.OR.IH0INSPC(1).EQ.4) THEN
*. EN diagonal
           I12 = 2
         END IF
       END IF
       ECOREX = SHIFT_DIA
*
* Partitioning and blockstructure of CI vector
*
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      NTTS = MXNTTS
      CALL MEMMAN(KLCLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLCLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLCI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLCIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM(ISSPC,WORK(KLCIOIO))
*
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTP),0)
*
*. Batches  of C vector
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KLCBLTP),WORK(KNSTSO(IATP)),
     &              WORK(KNSTSO(IBTP)),
     &              NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLCIOIO),
     &              ISMOST(1,ISSM),NBATCH,WORK(KLCLBT),
     &              WORK(KLCLEBT),WORK(KLCI1BT),WORK(KLCIBT),0,
     &              ITTSS_ORD)
*. Number of BLOCKS
        NBLOCK = IFRMR(WORK(KLCI1BT),1,NBATCH)
     &         + IFRMR(WORK(KLCLBT),1,NBATCH) - 1
C?      WRITE(6,*) ' HINTV : Number of blocks ', NBLOCK
        CALL GASDIAT(VECIN,LUSC3,ECOREX,ICISTR,I12,
     &               WORK(KLCBLTP),NBLOCK,WORK(KLCIBT),IDUMMY)
C      CALL GASDIAT(VECIN,ISSM,ISSPC,LUSC3,ECOREX,ICISTR,I12)
*. Clean up time
       IF(IH0INSPC(1).EQ.1.OR.IH0INSPC(1).EQ.3) THEN
*. MP operator
         CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
         CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
       END IF
*
* 2 : Solve linear set of equations
* ==================================
*
      ZERO = 0.0D0
      IF(LBLK.GT.0 ) THEN
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        WRITE(6,*) ' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1'
        Call Abend2(' PRESENT SCHEME DOES NOT WORK FOR ICISTR = 1')
*. Two vectors can be in core
*. Initial Guess on LUX
        NDIM = LBLK
        CALL SETVEC(VECOUT,ZERO,NDIM)
        CALL REWINE(LUX,-1)
        CALL TODSC_LUCI(VECOUT,NDIM,-1,LUX)
*. Right hand side on LUSC34
        CALL REWINE(LUSC34,-1)
        CALL TODSC_LUCI(VECIN,NDIM,-1,LUSC34)
*. Max number of its and convergence thresholds are picked up from
* corresponding eigenvalue info
        CONVER = .FALSE.
        TEST = SQRT(THRES_E) * SQRT(INPROD(VECIN,VECIN,NDIM))
        ILNPRT = MAX(NTEST-10,0)
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXITLE = ',MXITLE
        CALL MINGCG(MV8,LUX,LUSC34,LUSC35,LUSC3,VECIN,VECOUT,
     &              MXIT_LOC,CONVER,TEST,SHIFT,ERROR,NDIM,
     &              LUPROJ,LUPROJ2,ILNPRT)
        CALL REWINE(LUX,-1)
        CALL FRMDSC_LUCI(VECOUT,NDIM,-1,LUX,IMZERO,IAMPACK)
*
         IF(NTEST.GE.5) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(VECOUT,1,NDIM,1,NDIM)
         END IF
*
      ELSE IF(LBLK.LE.0)   THEN
*
*. Use path allowing us to work with segments of vectors
*
*. Initial guess on LUX
        CALL SETVCD(LUSC3,LUX,VECOUT,ZERO,1,LBLK)
*. (Right hand side vector is assumed in place)
*. convergence threshold is picked up from
* corresponding eigenvalue info
        TEST =
     &  SQRT(THRES_E) * SQRT(INPRDD(VECIN,VECOUT,LURHS,LURHS,1,-1))
        ILNPRT = NTEST
        SHIFT2 = 0.0D0
        CONVER = .FALSE.
        MXIT_LOC = MXITLE
C?      WRITE(6,*) ' HINTV : MXIT_LOC ',MXIT_LOC
        CALL MICGCG(H0TVM,LUX,LURHS,LUSC34,LUSC35,LUSC37,LUSC3,
     &              VECIN,VECOUT,MXIT_LOC,
     &              CONVER,TEST,SHIFT2,ERROR,NDIM,LUPROJ,LUPROJ2,
     &              ILNPRT)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Solution to linear set of Equations '
          CALL WRTVCD(VECIN,LUX,1,LBLK)
        END IF
*
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'HINTV ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                   IABSPC,JABSPC,IABOP,JABOP,IIF,JDOH2,
     &                   IDOH2,IMZERO,IDIAG)
*
* Decide upon the treatment of matrix element
*
* <IASM IATP IBSM IBTP | H(apr) | JASM, JATP JBSM, JBTP>
*
* and do preparations (IIF = 1 )
* or counteract preparations (IIF = 2)
*
* Jeppe Olsen, The last day of January 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      IF((IABSPC.NE.JABSPC).OR.
     &   (IABSPC.EQ.JABSPC.AND.(IABOP.EQ.1.OR.IABOP.EQ.2).AND.
     &    (JASM.NE.IASM.OR.IATP.NE.JATP.OR.IBTP.NE.JBTP))) THEN
*. Zero
        IMZERO = 1
      ELSE
*. Not zero !
        IMZERO = 0
      END IF
*. Diagonal approximation?
        IF(IABOP.EQ.1.OR.IABOP.EQ.2) THEN
          IDIAG = 1
        ELSE
          IDIAG = 0
        END IF
*. Moller Plesset or normal operator ?
        IF(IABOP.EQ.1.OR.IABOP.EQ.3.OR.IABOP.EQ.5)THEN
          IMP = 1
        ELSE
          IMP = 0
        END IF
C     END IF
*. Two - or one- electron operator
      IF(IABOP.EQ.1.OR.IABOP.EQ.3) THEN
        IDOH2 = 0
      ELSE
C       IDOH2 = JDOH2
        IDOH2 = 1
      END IF
*. Put MP integrals in place ( Or put good old one-electron integrals
*. back where they belong
      IF(IMP.EQ.1) THEN
        CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
        CALL SWAPVE(WORK(KINT1O),WORK(KFIO),NINT1)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
         WRITE(6,*) ' HMATAPR speaking '
         WRITE(6,*) ' ================='
         WRITE(6,*) ' | IASM IATP IBSM IBTP > :',
     &              '|',IASM,IATP,IBSM,IBTP,'>'
         WRITE(6,*) ' | JASM JATP JBSM JBTP > :',
     &              '|',JASM,JATP,JBSM,JBTP,'>'
         WRITE(6,*) ' IABSPC,JABSPC :', IABSPC,JABSPC
         WRITE(6,*) ' IABOP ', IABOP
         WRITE(6,*) ' Results : IMP IDIAG IDOH2 IMZERO : ',
     &                IMP,IDIAG,IDOH2, IMZERO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE HTV(VECIN,VECOUT,LUIN,LUOUT)
*
* Full operator times vector
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "oper.inc"
*
* complete operator in action
      IPERTOP = 0
      I12 = 2
*
      CALL MV7(VECIN,VECOUT,LUIN,LUOUT)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE INVERT_BY_DIAG(A,B,SCR,VEC,NDIM)
*
* Invert symmetric  - hopefully nonsingular - matrix A
* by diagonalization
*
* Jeppe Olsen, Oct 97 to check INVMAT
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output matrix
      DIMENSION A(*)
*. Scratch matrices and vector
      DIMENSION B(*),SCR(*),VEC(*)
*
      NTEST = 01
*. Reform a to symmetric packed form
C          TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGNTP)
      SIGNTP = 1.0
      CALL TRIPAK_LUCI(A,SCR,1,NDIM,NDIM,SIGNTP)
*. Diagonalize
      CALL EIGEN_LUCI(SCR,B,NDIM,0,1)
      CALL COPDIA(SCR,VEC,NDIM,1)
      IF( NTEST .GE. 1 ) THEN
        WRITE(6,*) ' Eigenvalues of matrix : '
        CALL WRTMT_LU(VEC,NDIM,1,NDIM,1)
      END IF
*. Invert diagonal elements - without safety at the moment
      DO I = 1, NDIM
       IF(ABS(VEC(I)).GT.1.0D-15) THEN
         VEC(I) = 1.0D0/VEC(I)
       ELSE
         VEC(I) = 0.0D0
         WRITE(6,*) ' Singular mode inactivated '
       END IF
      END DO
*. and obtain inverse matrix by transformation
C     XDIAXT(XDX,X,DIA,NDIM,SCR)
      CALL XDIAXT(A,B,VEC,NDIM,SCR)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Inverse matrix from INVERSE_BY_DIAG'
        CALL WRTMT_LU(A,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MATPERT(H0,V,NDIM,NORD,EN,C,VEC1,VEC2,VEC3,ECORE)
*
* Perturbation expansion of simple  eigenvalue problem
*
* Explicit matrix version
*
*
*. Obtain corrections to energy and wawe functions
*
*. The normalization condition used is C(K)T  C(0) = 0
*
* The energy corrections are
*
* E(n) = Sum(I=1,N) C(0)TF(I)C(N-I)
*      _ SUM(I=0,N-1)SUM(J=1,N-I)E(N-I-J)C(0)T S(J) C(I)
*
*
*
* Jeppe Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. Input
      DIMENSION H0(NDIM**2),V(NDIM**2)
*. Input and output (C(0) is supposed to be delivered here
      DIMENSION C(NDIM,*)
*. Output
      DIMENSION EN(0:NORD)
*. Scratch
      DIMENSION VEC1(NDIM),VEC2(NDIM),VEC3(NDIM)
*
*. Zero order energy
C  MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
      CALL MATVCB(H0,C,VEC1,NDIM,NDIM,0)
      E0   = INPROD(VEC1,C,NDIM)
*
      WRITE(6,*) 'E0  = ', E0
      EN(0) = E0
*. Save diagonal of H0 - E(0) in VEC3
      DO I = 1, NDIM
        VEC3(I) = H0((I-1)*NDIM+I)-E0
      END DO
C?    WRITE(6,*) ' Zero order diagonal '
C?    CALL WRTMT_LU(VEC3,1,NDIM,1,NDIM)
*. And then start the iterations
      DO IORD = 1, NORD
*
*  =================
*. Energy correction
*  =================
*
* E(n) =  C(0)T V C(N-1)
        CALL MATVCB(V,C(1,IORD+1-1),VEC1,NDIM,NDIM,0)
        EN(IORD) = INPROD(C,VEC1,NDIM)
C?      WRITE(6,*) ' Energy correction I,E(I) ',IORD,EN(IORD)
*
*  ==========================
*. Wave function corrections
*  ==========================
*
* C(N) = (H(0)-E(0))-1 (-VC(N-1)
*                           +Sum(K=1,N)E(K)C(N-K))
        CALL MATVCB(V,C(1,IORD+1-1),VEC2,NDIM,NDIM,0)
        ONEM = -1.0D0
        CALL SCALVE(VEC2,ONEM,NDIM)
C?      write(6,*) ' first term to rhs '
C?      CALL WRTMT_LU(VEC2,1,NDIM,1,NDIM)
*
        ONE = 1.0D0
        DO K = 1, IORD
          CALL VECSUM(VEC2,VEC2,C(1,IORD+1-K),ONE,EN(K),NDIM)
        END DO
*. Check overlap with zero order state ( should be zero )
        OVLAP = INPROD(C(1,1),VEC2,NDIM)
        FACTOR = -OVLAP
C?      WRITE(6,*) ' OVLAP = ',OVLAP
        CALL VECSUM(VEC2,VEC2,C(1,1),ONE,FACTOR,NDIM)
*. Multiply with (H0(0)-E(0))-1
C            DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
        ZERO = 0.0D0
        CALL DIAVC2(VEC1,VEC2,VEC3,ZERO,NDIM)
*
        CALL COPVEC(VEC1,C(1,IORD+1),NDIM)
*
C?      WRITE(6,*) ' Eigenfunction correction ', IORD
C?      CALL WRTMT_LU(C(1,IORD+1),1,NDIM,1,NDIM)
      END DO
*
      WRITE(6,*) ' Energy corrections : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      WRITE(6,*) '   Order             Correction '
      WRITE(6,*) ' ===================================='
      DO IORD = 1, NORD
        WRITE(6,'(1H ,3X,I3,E20.8)')IORD,EN(IORD)
      END DO
*
      ETOT = E0 + ECORE
      DO IORD = 1, NORD
        ETOT = ETOT + EN(IORD)
      END DO
      WRITE(6,*) ' Zero-order energy ', E0 + ECORE
      WRITE(6,*) ' Sum(K=0,NORD) E(K) ', ETOT
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MICDV6(VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LU4,LU5,LU6,LU7,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRTXX,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,
     &                  E_CONV,IROOTHOMING,LUWRTOUT)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version
*
* From MICDV4, Jeppe Olsen, April 1997
*              Modified : Oct 97 : root homing added
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scratch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
* On input LU1 is supposed to hold initial guesses to eigenvectors
*
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "parluci.h"
       DIMENSION VEC1(*),VEC2(*)
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
*      IROOTHOMING : Do roothoming, i.e. select the
*      eigenvectors in iteration n+1 as the approximations
*      with largest overlap with the previous space
*
      DIMENSION FINEIG(nroot)
      LOGICAL CONVER,LUCEND
      REAL*8 INPRDD, INPROD
      DOUBLE PRECISION tottime,endtime,starttime
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
*     New characters used for timing in send of VEC from LU3
      CHARACTER CPULUSTEP*12, WALLLUSTEP*12
*     New characters used for timing in cpback of sigma vec and calc of
*     aproj(ij)
      CHARACTER CPUFSTEP*12, WALLFSTEP*12
*     New characters used for timing in step3
      CHARACTER CPUTSTEP*12, WALLTSTEP*12
*     XJEP is used for ROOTHOMING
      real(8), allocatable ::  xjep(:)
      integer, allocatable :: ixjep(:)
      character (len=13)   :: convergence
      logical, allocatable :: rtcnv(:)

      allocate(rtcnv(nroot))
      rtcnv  = .false.
      lucend = .FALSE.

      MYNEW_ID   = 0
      MYNEW_COMM = 0
C
      WRITE(LUWRTOUT,'(/10X,A)')
     &' **************************************************'
      WRITE(LUWRTOUT,'(10X,A)')
     &'     entering MICDV6 (sequential solver routine)   '
      WRITE(LUWRTOUT,'(10X,A )')
     &' **************************************************'
*
       IPICO = 0
       IF(IPICO.NE.0) THEN
         MAXVEC = MIN(MAXVEC,2)
       END IF
*
!      IPRT = 9999 ! IPRT = 0000  DEBUG
       IPRT = 0000
       lbl  = 0
c      IPRT = IPRTXX
#ifdef LUCI_DEBUG
         WRITE(LUWRTOUT,*) ' Initial set of eigenvectors '
         CALL REWINE(LU1,-1)
         DO IROOT = 1, NROOT
           WRITE(LUWRTOUT,1066) IROOT, LBLK, LBL
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
#endif
*
       IOLSTM = 1
       IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     & WRITE(LUWRTOUT,*) ' Inverse iteration modified Davidson '
       IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     & WRITE(LUWRTOUT,*) ' Normal Davidson method '
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(LUWRTOUT,*)' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
         WRITE(LUWRTOUT,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
         WRITE(LUWRTOUT,*) ' Raise MXCIV to be at least 2 * Nroot '
         WRITE(LUWRTOUT,*) ' Enforced stop on MICDV6 '
         Call Abend1( 20 )
       END IF
*
      if(iroothoming.eq.1)then
        write(luwrtout,*) ' Root homing performed '
        allocate(xjep( maxvec*3))
        allocate(ixjep(maxvec*3))
      end if
*
      KAPROJ = 1
      KFREE  = KAPROJ+ MAXVEC*(MAXVEC+1)/2
      TEST   = 1.0D-4
      CONVER = .FALSE.
*
* ===================
*.Initial iteration
* ===================
      ITER = 1
      CALL REWINE(LU1,-1)
      CALL REWINE(LU2,-1)
*
      DO IVEC = 1,NINVEC
        CALL REWINE(LU5,-1)
        CALL REWINE(LU6,-1)
        CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
*       Check timings for sigma vector generation
        CALL GETTIM(CPUSIG1,WALLSIG1)

        CALL MV7(VEC1,VEC2,LU5,LU6)

        CALL GETTIM(CPUSIG2,WALLSIG2)
        WALLTID = SECTID(WALLSIG2-WALLSIG1)
        if(ivec .eq. 1)then
          write(luwrtout,9777) WALLTID
        else
          write(luwrtout,9778) WALLTID
        end if
*. Move sigma to LU2, LU2 is positioned at end of vector IVEC - 1
        CALL REWINE(LU6,-1)
        CALL COPVCD(LU6,LU2,VEC1,0,LBLK)
*. Projected matrix
        CALL REWINE(LU2,-1)
          DO JVEC = 1, IVEC
             CALL REWINE(LU5,-1)
             IJ = IVEC*(IVEC-1)/2 + JVEC
             APROJ(IJ) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
          END DO ! loop over IVEC
       END DO ! loop over NINVEC
*
       IF( IPRT .GE.5 ) THEN
         WRITE(LUWRTOUT,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*. Diagonalize initial projected matrix
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
       END DO
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(LUWRTOUT,*) 
     &   ' Eigenvalues of initial iteration (with shift)'
         WRITE(LUWRTOUT,'(5F18.13)')
     &   ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
       END IF
       IF( IPRT  .GE. 5 ) THEN
        WRITE(luwrt,*) ' Eigenvalues of initial iteration (no shift)'
        WRITE(luwrt,'(5F18.13)')
     &  ( EIG(1,IROOT),IROOT=1,NROOT)
       END IF
*. Transform vectors on LU1 so they become the actual
*. eigenvector approximations
C     TRAVCD(VEC1,VEC2,NVECIN,NVECOUT,LUIN,LUOUT,
C    &                  ICOPY,LBLK,LUSCR1,LUSCR2)
       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU1,LU3,1,LBLK,LU4,LU5)
*. And the sigma vectors
       CALL REWINE(LU3,-1)
       CALL TRAVCD(VEC1,VEC2,AVEC,NINVEC,NROOT,LU2,LU3,1,LBLK,LU4,LU5)
*
       IF(IPRT.GE.600) THEN
         WRITE(LUWRTOUT,*) ' Initial set of eigenvectors '
         CALL REWINE(LU1,-1)
         DO IROOT = 1, NROOT
           WRITE(LUWRTOUT,1066) IROOT, LBLK, LBL
 1066 FORMAT(/'lets see IROOT =', 3I4/)
           CALL WRTVCD(VEC1,LU1,0,LBLK)
         END DO
*
         WRITE(LUWRTOUT,*) ' Initial set of sigma vectors '
         CALL REWINE(LU2,-1)
         DO IROOT = 1, NROOT
           WRITE(LUWRTOUT,1067) IROOT
 1067 FORMAT(/'lets see IROOT again=', I4/)
           CALL WRTVCD(VEC1,LU2,0,LBLK)
         END DO
       END IF
*. And the corresponding Hamiltonian matrix, no problems
*. with numerical stabilities, so just use eigenvalues
       CALL dzero(APROJ,NROOT*(NROOT+1)/2)
       DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
       END DO
*
       NVEC = NROOT
CJAN20 IF (MAXIT .EQ. 1 ) GOTO  901
 12    CONTINUE 
       IF (MAXIT .EQ. 1 ) GOTO  1001
*
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
        CALL GETTIM(CPUITR1,WALLITR1)
        write(LUWRTOUT,*)
        write(LUWRTOUT,'(A21,3X,I3)') ' Info from iteration ',ITER
        write(LUWRTOUT,*) '_______________________'
        ITER = ITER + 1
*
* ===============================
*.1 New directions to be included
* ===============================
*
* 1.1 : R = H*X - EIGAPR*X
*
       IADD = 0
       CONVER = .TRUE.
*
       CALL REWINE(LU1,-1)
       CALL REWINE(LU2,-1)
*
       DO 100 IROOT = 1, NROOT
*. Save current eigenvector IROOT on LU7
         CALL SKPVCD(LU1,IROOT-1,VEC1,1,LBLK)
         CALL REWINE(LU7,-1)
         CALL COPVCD(LU1,LU7,VEC1,0,LBLK)
*. Calculate (HX - EX ) and store on LU5
*. Current eigenvector is  on LU7, corresponding sigma vector
*. on LU2
         EIGAPR = EIG(ITER-1,IROOT)
         ONE = 1.0D0
*
         CALL REWINE(LU7,-1)
         CALL REWINE(LU5,-1)
         FACTOR = - EIGAPR
         CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU2,LU7,LU5,0,LBLK)
*
         IF ( IPRT  .GE. 10 ) THEN
           WRITE(LUWRTOUT,*) '  ( HX - EX ) '
           CALL WRTVCD(VEC1,LU5,1,LBLK)
         END IF
*  Strange place to put convergence but ....
         RNORM = SQRT( INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK) )
         RNRM(ITER-1,IROOT) = RNORM
         WRITE(LUWRTOUT,'(A19,I8,1p,1E15.5,0p,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
*
         CALL FLSHFO(LUWRTOUT)
*
         IF(RNORM.LT. TEST .OR.
     &      (ITER.GT.2.AND.
     &      ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.E_CONV)) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
* =====================================================================
*. 1.2 : Multiply with inverse Hessian approximation to get new directio
* =====================================================================
*. (H0-E) -1 *(HX-EX) on LU6
CSK         print*,'NP1,NP2,NQ ',NP1,NP2,NQ
         IF( .NOT. RTCNV(IROOT) ) THEN
           IF(IPRT.GE.3) THEN
             WRITE(LUWRTOUT,*) ' Correction vector added for root',IROOT
           END IF
           IADD = IADD + 1
           CALL REWINE(LUDIA,-1)
           CALL REWINE(LU5,-1)
           CALL REWINE(LU6,-1)
           CALL H0M1TD(LU6,LUDIA,LU5,LBLK,NP1+NP2+NQ,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XH0IX,
     &                 NP1,NP2,NQ,VEC1,VEC2,IPRT)
C     H0M1TD(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
C    &                  H0,SHIFT,WORK,XH0PSX,
C    &                  NP1,NP2,NQ,VEC1,VEC2,NTESTG)
           IF ( IPRT  .GE. 600) THEN
             WRITE(LUWRTOUT,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU6,1,LBLK)
           END IF
CSK        print*,'IOLSTM ',IOLSTM
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
* (H0 - E )-1  * X on LU5
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
             CALL REWINE(LUDIA,-1)
*
             CALL H0M1TD(LU5,LUDIA,LU7,LBLK,Np1+Np2+NQ,
     &                   IPNTR,H0,-EIGAPR,H0SCR,XH0IX,
     &                   NP1,NP2,NQ,VEC1,VEC2,IPRT)
*
* Gamma = X(T) * (H0 - E) ** -1 * X
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
             GAMMA = INPRDD(VEC1,VEC2,LU5,LU7,0,LBLK)
CSK             WRITE(LUWRTOUT,*) ' Gamma = ', Gamma
             IF(IPRT.GE.1000) WRITE(LUWRTOUT,*) ' Gamma = ', Gamma
* is X an eigen vector for (H0 - 1 ) - 1
             CALL REWINE(LU5,-1)
             CALL REWINE(LU7,-1)
              VNORM =
     &        SQRT(VCSMDN(VEC1,VEC2,-GAMMA,1.0D0,LU7,LU5,0,LBLK))
CSK              write(LUWRTOUT,*) ' VNORM = ', VNORM
              IF(IPRT.GE.1000) write(LUWRTOUT,*) ' VNORM = ', VNORM
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
             IF(IPRT.GE.5) WRITE(LUWRTOUT,*) ' Olsen Correction active '
CSK                WRITE(LUWRTOUT,*) ' Olsen Correction active '
                DELTA = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                FACTOR = -(DELTA/GAMMA)
                IF(IPRT.GE.5) WRITE(LUWRTOUT,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(LUWRTOUT,*)   DELTA,GAMMA,FACTOR
CSK                WRITE(LUWRTOUT,*) ' DELTA,GAMMA,FACTOR'
CSK                WRITE(LUWRTOUT,*)   DELTA,GAMMA,FACTOR
                CALL VECSMD(VEC1,VEC2,1.0D0,FACTOR,LU6,LU5,LU7,1,LBLK)
                CALL COPVCD(LU7,LU6,VEC1,1,LBLK)
*
                IF(IPRT.GE.600) THEN
                  WRITE(LUWRTOUT,*) ' Modified trial vector '
                  CALL WRTVCD(VEC1,LU6,1,LBLK)
                END IF
*
              END IF
            END IF
*. 1.3 Orthogonalize to all previous vectors
*.. Vectors on LU1
*
           CALL REWINE( LU1 ,LBLK)
           DO IVEC = 1, NROOT
               CALL REWINE(LU6,-1)
               WORK(IVEC) = -INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
           END DO
           IF(IPRT.GE.1000) THEN
             Write(LUWRTOUT,*) ' Overlap with vectors on LU1'
             call WRTMT_LU(work,1,nroot,1,nroot)
           END IF
           ONE = 1.0D0
           CALL REWINE(LU1,-1)
           CALL MVCSMD(LU1,WORK,LU7,LU5,VEC1,VEC2,NROOT,1,LBLK)
           CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU6,LU5,1,LBLK)
           IF(IPRT.GE.1000) THEN
             write(LUWRTOUT,*) ' orthogonalized to vectors on LU1'
             CALL WRTVCD(VEC1,LU5,1,LBLK)
           END IF
*.. Vectors on LU3
           IF(NVEC+IADD-1-NROOT.GT.0) THEN
             CALL REWINE( LU3 ,LBLK)
             DO IVEC = 1, NVEC+IADD-1-NROOT
                 CALL REWINE(LU6,-1)
                 WORK(IVEC) = -INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
             END DO
             ONE = 1.0D0
             CALL REWINE(LU3,-1)
             CALL MVCSMD(LU3,WORK,LU7,LU6,VEC1,VEC2,NVEC+IADD-1-NROOT,
     &            1,LBLK)
             CALL VECSMD(VEC1,VEC2,ONE,ONE,LU7,LU5,LU6,1,LBLK)
           ELSE
             CALL REWINE(LU3,-1)
             CALL COPVCD(LU5,LU6,VEC1,1,LBLK)
           END IF
*. 1.4 Normalize vector
           SCALE = INPRDD(VEC1,VEC1,LU6,LU6,1,LBLK)
           FACTOR = 1.0D0/SQRT(SCALE)
           CALL REWINE(LU6,LBLK)
           CALL SCLVCD(LU6,LU3,FACTOR,VEC1,0,LBLK)
           IF(IPRT.GE.600) THEN
             CALL SCLVCD(LU6,LU7,FACTOR,VEC1,1,LBLK)
             WRITE(LUWRTOUT,*) '   normalized     (D-E)-1 *( HX - EX ) '
             CALL WRTVCD(VEC1,LU7,1,LBLK)
           END IF
*
         END IF
  100 CONTINUE
  101 CONTINUE
*
      IF( CONVER ) THEN
         GOTO  1001
      END IF
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
*
* ====================================================
*  2 : Optimal combination of new and old directions
* ====================================================
*
*  2.1: Multiply new directions with matrix
*
      IF(IPRT.GE.1000) THEN
       IF (luci_MYPROC.EQ.luci_MASTER) THEN
        WRITE(LUWRTOUT,*) ' Vectors on LU3'
        WRITE(LUWRTOUT,*) 'NVEC-NROOT+IADD =',NVEC-NROOT+IADD
        CALL REWINE(LU3,-1)
        DO IVEC = 1, NVEC-NROOT+IADD
          CALL WRTVCD(VEC1,LU3,0,LBLK)
        END DO
       END IF
       CALL REWINE(LU3,-1)
      END IF
*
      CALL SKPVCD(LU3,NVEC-NROOT,VEC1,1,LBLK)
      CALL SKPVCD(LU4,NVEC-NROOT,VEC1,1,LBLK)
*
*     
 150  DO IVEC = 1, IADD
        CALL REWINE(LU5,LBLK)
        CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
*Check timings for sigma vector generation
        CALL GETTIM(CPUSIG1,WALLSIG1)
        CALL MV7(VEC1,VEC2,LU5,LU6)
        CALL GETTIM(CPUSIG2,WALLSIG2)
        WALLTID = SECTID(WALLSIG2-WALLSIG1)
        if(ivec .eq. 1)then
          WRITE(LUWRTOUT,9400) WALLTID
        else
          WRITE(LUWRTOUT,9401) WALLTID
        end if
        CALL REWINE(LU6,LBLK)
        CALL COPVCD(LU6,LU4,VEC1,0,LBLK)
*2.2 Augment projected matrix
        CALL REWINE( LU1,LBLK)
        DO JVEC = 1, NROOT
          CALL REWINE(LU6,LBLK)
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU1,LU6,0,LBLK)
        END DO

        CALL REWINE(LU3,LBLK)
C       DO JVEC = NROOT+1, NVEC+IADD
        DO JVEC = NROOT+1, NVEC+IVEC
         CALL REWINE(LU6,LBLK)
         IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
         APROJ(IJ) = INPRDD(VEC1,VEC2,LU3,LU6,0,LBLK)
        END DO
      END DO
*     /\ End do over new trial vectors
*
*. 2.3 Diagonalize projected matrix
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NVEC,0,1)
*. Test : transform projected matrix
C TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
C     CALL TRAN_SYM_BLOC_MAT(APROJ,AVEC,1,NVEC,XJEP(1000),XJEP(1))
C     WRITE(LUWRTOUT,*) ' Explicitly transformed matrix'
C     CALL PRSYM(XJEP(1000),NVEC)

      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
        C1NRM = SQRT(C0VAR**2+C1VAR**2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -(C1NRM/SQRT(1.0D0+C1NRM**2))
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(LUWRTOUT,*) 
     &  ' Var and Pert solution, energy and coefficients'
        WRITE(LUWRTOUT,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(LUWRTOUT,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
*
      IF(IROOTHOMING.EQ.1) THEN
*
*. Reorder roots so the NROOT with the largest overlap with
*  the original roots become the first
*
*. Norm of wavefunction in previous space
       DO IVEC = 1, NVEC
         XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                AVEC(1+(IVEC-1)*NROOT),NROOT)
       END DO
       WRITE(LUWRTOUT,*)
     & ' Norm of projections to previous vector space '
       CALL WRTMT_LU(XJEP,1,NVEC,1,NVEC)
*. My sorter arranges in increasing order, multiply with minus 1
*  so the eigenvectors with largest overlap comes out first
       ONEM = -1.0D0
       CALL SCALVE(XJEP,ONEM,NVEC)
       CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
       IF(NSORT.LT.NVEC) THEN
         WRITE(LUWRTOUT,*) ' Warning : Some elements lost in sorting '
         WRITE(LUWRTOUT,*) ' NVEC,NSORT = ', NSORT,NVEC
       END IF
       IF(IPRT.GE.3) THEN
         WRITE(LUWRTOUT,*) ' New roots choosen as vectors '
         CALL IWRTMA(IXJEP,1,NROOT,1,NROOT)
       END IF
*. Reorder
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         CALL COPVEC(AVEC(1+(IOLD-1)*NVEC),XJEP(1+(INEW-1)*NVEC),NVEC)
       END DO
       CALL COPVEC(XJEP,AVEC,NROOT*NVEC)
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
       END DO
       DO INEW = 1, NVEC
         WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
       END DO
*
       IF(IPRT.GE.3) THEN
         WRITE(LUWRTOUT,*) ' Reordered WORK and AVEC arrays '
         CALL PRSYM(WORK,NVEC)
         CALL WRTMT_LU(AVEC,NVEC,NVEC,NVEC,NVEC)
       END IF
*
      END IF
*     ^ End of root homing procedure
      DO IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
      END DO
*
      IF(IPRT .GE. 3 ) THEN
        WRITE(LUWRTOUT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
        WRITE(LUWRTOUT,'(5F18.13)')
     &  ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
        WRITE(LUWRTOUT,'(A)') ' Norm of Residuals (Previous it) '
        WRITE(LUWRTOUT,'(5F18.13)')
     &  ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
      END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRTOUT,*) ' Projected matrix'
        CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRTOUT,*) ' eigen pairs '
        WRITE(LUWRTOUT,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMT_LU(AVEC,NVEC,NROOT,NVEC,NROOT)
      END IF
*
**  perhaps reset or assemble converged eigenvectors
*
  901 CONTINUE
*
*. Reset
*
*
* case 1 : Only NROOT vectors can be stored
*          save current eigenvector approximations
* Case 2 : Atleast 2*NROOT eigenvectors can be saved
*          Current eigenvactor approximations+
*          vectors allowing generation of previous approxs.
*
*
C     IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
*
*. 3.1 : Current wave function approximations, collect on LU7
*
*  check timing of part3 in loop iteration
       CALL GETTIM(CPUSTEP3,WALLSTEP3)
*
        IF(IPRT.GE.1000) THEN
        write(LUWRTOUT,*) ' I am going to reset '
        write(LUWRTOUT,*) ' nroot, nvec ', nroot,nvec
        END IF
        IF(IPRT.GE.1000) THEN
          WRITE(LUWRTOUT,*) ' Initial vectors on LU1'
          CALL REWINE(LU1,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
          WRITE(LUWRTOUT,*) ' Initial vectors on LU3'
          CALL REWINE(LU3,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU3,0,LBLK)
          END DO
        END IF
*
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU1 to LU5
          CALL REWINE(LU1,-1)
          CALL MVCSMD(LU1,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of c reset, part 1'
*. and add LU3 to LU5
          CALL REWINE(LU3,-1)
          ONE = 1.0D0
          CALL MVCSMD2(LU3,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of c reset, part 2'
C              MVCSMD2(LUIN,FAC,FACLUOUT,LUOUT,LUSCR,
C    &         VEC1,VEC2,NVEC,IREW,LBLK)
          XNORM = INPRDD(VEC1,VEC1,LU5,LU5,1,LBLK)
          CALL REWINE(LU5,LBLK)
          SCALE  = 1.0D0/SQRT(XNORM)
CSK          WRITE(LUWRTOUT,*) ' SCALE = ', SCALE
          WORK(IROOT) = SCALE
*. scale LU5 => LU7
          CALL REWINE(LU5,-1)
          CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of c reset, part 3'
        END DO
*. Transfer C vectors from LU7 to LU1
        CALL REWINE( LU7,LBLK)
        CALL REWINE( LU1,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCD(LU7,LU1,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(LUWRTOUT,*) ' end of c reset, part 4'
        IF(IPRT.GE.1000) THEN
          WRITE(LUWRTOUT,*) ' Reset C-vectors on LU1 '
          CALL REWINE(LU1,-1)
          DO IVEC = 1, NROOT
             CALL WRTVCD(VEC1,LU1,0,LBLK)
          END DO
        END IF
*. 3.2 : corresponding sigma vectors
        CALL REWINE( LU7,LBLK)
        DO IROOT = 1, NROOT
*. From LU2 to LU5
          CALL REWINE(LU2,-1)
          CALL MVCSMD(LU2,AVEC((IROOT-1)*NVEC+1),
     &    LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of s reset, part 1'
*. and add LU4 to LU5
          CALL REWINE(LU4,-1)
          CALL MVCSMD2(LU4,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &    LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of s reset, part 2'
          SCALE  = WORK(IROOT)
*. scale LU5 => LU7
          CALL REWINE(LU5,-1)
          CALL SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &    write(LUWRTOUT,*) ' end of s reset, part 3'
        END DO
*. Transfer sigma  vectors from LU7 to LU2
        CALL REWINE( LU7,LBLK)
        CALL REWINE( LU2,LBLK)
        DO IVEC = 1,NROOT
          CALL COPVCD(LU7,LU2,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(LUWRTOUT,*) ' end of s reset, part 4'
        NNVEC = NROOT
*. Overlap between first vectors on LU1 and LU2
        E11 = INPRDD(VEC1,VEC2,LU1,LU2,1,LBLK)
CSK        WRITE(6,*) ' <Lu1!Lu2> = ', E11
        IF(3*NROOT.LE.MAXVEC) THEN
*
*. Orthogonalize the
*. last set of correction vectors to the current
*. eigenvectors on LU1, and save on LU2
*. Overlap with root approximations
*. Start of last set of trial vectors
          ISTART = NVEC-NROOT-IADD+1
          CALL SKPVCD(LU3,ISTART-1,VEC1,1,LBLK)
*
          CALL REWINE(LU5,-1)
          DO JVEC = 1, IADD
*. Orthogonalize to vectors on LU1
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU3,LU7,VEC1,0,LBLK)
            IF(IPRT.GE.1000) THEN
            write(LUWRTOUT,*) ' Initial vector on LU7'
            CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
            CALL REWINE(LU1,-1)
            DO IROOT = 1, NROOT
              CALL REWINE(LU7,-1)
              WORK(IROOT+(JVEC-1)*2*NROOT) =
     &        -INPRDD(VEC1,VEC2,LU7,LU1,0,LBLK)
            END DO
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' c, part1 finito'
*. And to trial vectors on LU5
            CALL REWINE(LU5,-1)
            DO KVEC = 1, JVEC-1
              CALL REWINE(LU7,-1)
              WORK(NROOT+KVEC+(JVEC-1)*2*NROOT) =
     &        -INPRDD(VEC1,VEC2,LU7,LU5,0,LBLK)
            END DO
            WORK(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' c, part2 finito'
*
            ONE = 1.0D0
            CALL MVCSMD2(LU1,WORK(1+(JVEC-1)*2*NROOT),ONE ,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            ONE = 1.0D0
            CALL MVCSMD2(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000) THEN
              write(LUWRTOUT,*) ' c, part4 finito'
              write(LUWRTOUT,*) ' Vector after sec ort'
              CALL WRTVCD(VEC1,LU7,1,LBLK)
            END IF
*
            FACTOR = INPRDD(VEC1,VEC1,LU7,LU7,1,LBLK)
            SCALE = 1.0D0/SQRT(FACTOR)
            CALL SCALVE(WORK((JVEC-1)*2*NROOT+1),SCALE,
     &           NROOT+JVEC)
            CALL REWINE(LU7,-1)
            CALL SCLVCD(LU7,LU5,SCALE,VEC1,0,LBLK)
C                SCLVCD(LU5,LU7,SCALE,VEC1,0,LBLK)
           IF(IPRT.GE.1000)
     &       write(LUWRTOUT,*) ' c, part6 finito'
         END DO
*        /\ End of loop over orthogonalized directions
*. Transfer modified directions to LU3
         CALL REWINE(LU3,-1)
         CALL REWINE(LU5,-1)
         DO JVEC =1, IADD
           CALL COPVCD(LU5,LU3,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000) THEN
           write(LUWRTOUT,*) ' c, part7 finito'
           WRITE(LUWRTOUT,*) ' Additional trial vectors on LU3'
           CALL REWINE(LU3,-1)
           DO JVEC = 1, IADD
            CALL WRTVCD(VEC1,LU3,0,LBLK)
           END DO
         END IF
* Sigma vectors corresponding to orthogonalized directions
         CALL SKPVCD(LU4,ISTART-1,VEC1,1,LBLK)
         CALL REWINE(LU5,-1)
         DO JVEC = 1, IADD
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU4,LU7,VEC1,0,LBLK)
*
            FACT = WORK(NROOT+JVEC+(JVEC-1)*2*NROOT)
            CALL MVCSMD2(LU2,WORK(1+(JVEC-1)*2*NROOT),FACT,
     &      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
*
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' s, part 1 finito '
            ONE = 1.0D0
            CALL MVCSMD2(LU5,WORK(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &           LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' s, part 2 finito '
            CALL REWINE(LU7,-1)
            CALL COPVCD(LU7,LU5,VEC1,0,LBLK)
            IF(IPRT.GE.1000)
     &      write(LUWRTOUT,*) ' s, part 3 finito '
         END DO
*        /\ End of loop over orthogonalized directions
*. Copy back to LU4
         CALL REWINE(LU4,-1)
         CALL REWINE(LU5,-1)
         DO JVEC = 1, IADD
            CALL COPVCD(LU5,LU4,VEC1,0,LBLK)
         END DO
         IF(IPRT.GE.1000)
     &   write(LUWRTOUT,*) ' s, part 4 finito '
C         WRITE(LUWRTOUT,*) 'I was in 3.1 setting NNVEC=NROOT+IADD', 
C     &luci_myproc,(NROOT + IADD)
         NNVEC = NROOT + IADD
       END IF
       IF(IPRT.GE.1000) THEN
         WRITE(LUWRTOUT,*) ' Additional sigma vectors on LU4'
         CALL REWINE(LU4,-1)
         DO JVEC = 1, IADD
           CALL WRTVCD(VEC1,LU4,0,LBLK)
         END DO
       END IF
*
*      /\ End if more than NROOT vectors to be reset
       NVEC = NNVEC
*       write(LUWRTOUT,*)'NVEC(1)',NVEC NVEC = IADD (iadd can change!!!)+ NROOT
C      END IF
*.     ^ End if reset
*. New subspace
*. Calculate subspace Hamiltonian from actual vectors
*. on disc
       IF(IPRT.GE.1000) write(LUWRTOUT,*) ' Subspace hamiltonian'
       CALL REWINE(LU1,-1)
       CALL REWINE(LU3,-1)
CSK       ididaproj = 0
CSK       write(LUWRTOUT,*)'NVEC = ', NVEC
       DO IVEC = 1, NVEC
*
         CALL REWINE(LU5,-1)
         IF(IVEC.LE.NROOT) THEN
           CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
         ELSE
           CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
         END IF
*
         CALL REWINE(LU2,-1)
         DO JVEC = 1, MIN(IVEC,NROOT)
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
         END DO
         CALL REWINE(LU4,-1)
         DO JVEC = NROOT+1,IVEC
           CALL REWINE(LU5,-1)
           IJ = IVEC*(IVEC-1)/2+JVEC
           APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
         END DO
CSK        ididaproj = ididaproj + 1
       END DO
       if (IPRT.ge.2) then
         write(LUWRTOUT,*) ' Reset hamiltonian'
         call prsym(aproj,nvec)
       end if
CSK       write(LUWRTOUT,*)'actual dimension of aproj',MAXVEC*(MAXVEC+1)/2
CSK       write(LUWRTOUT,*)'calculated elements of aproj',ididaproj
*  check timing of part3 in loop iteration
       CALL GETTIM(CPUSTEP4,WALLSTEP4)
       WALLTSTEP = SECTID(WALLSTEP4-WALLSTEP3)
       WRITE(LUWRTOUT,9600) WALLTSTEP
*
 998   CONTINUE
*  Timing of this iteration
       CALL GETTIM(CPUITR2,WALLITR2)
       WALLTID = SECTID(WALLITR2-WALLITR1)
       WRITE(LUWRTOUT,9300) WALLTID
*. End of resetting business
 999  CONTINUE
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
* ( End of loop over iterations )
C
C
      IF( .NOT. CONVER ) THEN
*       CONVERGENCE WAS NOT OBTAINED
        IF(IPRT .GE. 2 )
     &  WRITE(LUWRTOUT,1170) MAXIT
 1170   FORMAT('   Convergence was not obtained in ',I3,' iterations')
      ELSE
*       CONVERGENCE WAS OBTAINED
        ITER = ITER - 1
        IF (IPRT .GE. 2 )
     &  WRITE(LUWRTOUT,1180) ITER
 1180   FORMAT(1X,' Convergence was obtained in ',I3,' iterations')
      END IF
*
      do iroot = 1, nroot
        write(LUWRTOUT,'(/1x,a)')    '-------------------'
        write(LUWRTOUT,'( 1x,a,i6)') 'root number =',IROOT
        write(LUWRTOUT,'(1x,a/)')    '-------------------'
        do i = 1, iter
          write(LUWRTOUT,'( 1x,a9,i4,f20.10,4x,1p,e10.3)')
     & 'iteration',i,eig(i,iroot)+eigshf,rnrm(i,iroot)
        end do
      end do

      write(LUWRTOUT,'(/1x,a)') '**********************************'//
     &           '******************************'
      write(LUWRTOUT,*) '  iter   root         energy          '//
     &           'residual     convergence'
      write(LUWRTOUT,'(1x,a )') '**********************************'//
     &           '******************************'
      DO 1601 IROOT = 1, NROOT
         FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
         if(rtcnv(iroot)) then
           convergence = ' converged'
         else
           convergence = 'NOT converged'
         end if
         WRITE(LUWRTOUT,'(i7,i7,f20.10,4x,1p,E10.3,4x,a)')
     &         ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &         convergence
 1601 CONTINUE
*
      if(iroothoming.eq.1)then
        deallocate(ixjep)
        deallocate( xjep)
      end if
      deallocate(rtcnv)

      CALL FLSHFO(LUWRTOUT)
*
      RETURN
 9300 FORMAT(' >>>  WALL TIME FOR CURRENT ITERATION            : ',A)
 9400 FORMAT(/' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9401 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9600 FORMAT(/' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION   : ',A)
 9777 FORMAT(/' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
 9778 FORMAT(' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MICGCG(MV8,LU1,LU2,LU3,LU4,LU5,LUDIA,VEC1,VEC2,
     &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
     &                  LUPROJ,LUPROJ2,IPRT)
*
* Solve set of linear equations
*
*             AX = B
*
* with preconditioned conjugate gradient method for
* case where two complete vectors can be stored in core
*
* Initial appriximation to solution must reside on LU1
* LU2 must contain B.All files are  overwritten
*
*
* Final solution vector is stored in LU1
* A scalar w can be added to the diagonal of the preconditioner
*
* If LUPROJ .NE. 0 , the optimization subspace is restricted to be orthogonal
* to the first vector in LUPROJ.
* The vector used to orthogonalize is saved on LUPROJ2
*
* Version using blocks of vectors
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*),ERROR(MAXIT+1)
      REAL*8 INPRDD
      LOGICAL CONVER
*
      EXTERNAL MV8
*
      CALL QENTER('MICGC')
      NTEST = 005
      NTEST = MAX(NTEST,IPRT)
      IF(NTEST.GE.5) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' Welcome to MICGCG '
        WRITE(6,*) ' =================='
        WRITE(6,*)
*
C?    WRITE(6,*) ' NTEST ,LU1,LU2,LU3 = ', NTEST,LU1,LU2,LU3
      END IF
      CONVER = .FALSE.
      ITER = 1
*
      LBLK = -1
*
      ONE = 1.0D0
      ONEM = -1.0D0
      ZERO = 0.0D0
*. Overlap between LUPROJ and LUPROJ2
      IF(LUPROJ.GT.0) THEN
        X12 = INPRDD(VEC1,VEC2,LUPROJ,LUPROJ2,1,LBLK)
      ELSE
        X12 = 0.0D0
      END IF
C?    WRITE(6,*) ' MICGCG : X12 = ', X12
*
* =============
* Initial point
* =============
*
*.R = B - (A)*X on LU2
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector on LU1 '
        CALL WRTVCD(VEC1,LU1,1,LBLK)
        WRITE(6,*) ' Vector on LU2 '
        CALL WRTVCD(VEC1,LU2,1,LBLK)
      END IF
      CALL MV8(VEC1,VEC2,LU1,LU3)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector on LU3 '
        CALL WRTVCD(VEC1,LU3,1,LBLK)
      END IF
*
C          VECSMD(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
      CALL VECSMD(VEC1,VEC2,ONE,ONEM,LU2,LU3,LU4,1,LBLK)
      CALL COPVCD(LU4,LU2,VEC1,1,LBLK)
*
*
      RNORM = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
      ERROR(1) = SQRT(RNORM)
      IF(ERROR(1).LE.TEST) THEN
*. Convergence in one shot- you are lucky or have
* supplied a vaninshing RHS
        NITER = 0
        CONVER = .TRUE.
        GOTO 1001
      END IF
*
*. Preconditioner H times initial residual, H * R on LU4
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Diagonal and input to diagonal '
       CALL WRTVCD(VEC1,LUDIA,1,LBLK)
       CALL WRTVCD(VEC1,LU2  ,1,LBLK)
       WRITE(6,*) ' SHIFT = ', W
      END IF
      CALL DMTVCD(VEC1,VEC2,LUDIA,LU2,LU4,W,1,1,LBLK)
      IF(LUPROJ.NE.0) THEN
        OVLAP = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
        FACTOR = -OVLAP/X12
        CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU4,LUPROJ2,LU3,1,LBLK)
        CALL COPVCD(LU3,LU4,VEC1,1,LBLK)
        OVLAP2 = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
        WRITE(6,*) ' Updated overlap of trial vector ', OVLAP2
      END IF
*. GAMMA = <R!H!R>
      GAMMA = INPRDD(VEC1,VEC2,LU2,LU4,1,LBLK)
*. P = RHO * H*R on LU3
      RHO = 1.0D0
      CALL SCLVCD(LU4,LU3,RHO,VEC1,1,LBLK)
*.S = AP on LU4
      CALL MV8(VEC1,VEC2,LU3,LU4)
*
* ====================
* Loop over iterations
* ====================
*
      NITER = 0
      DO 1000 ITER = 1, MAXIT
*
* Vectors on files :
*     X on LU1
*     R on LU2
*     P on LU3
*  S=AP on LU4
*     H on LUDIA

        NITER = NITER + 1
       IF ( NTEST .GE. 10 )
     & WRITE(6,*) ' INFORMATION FROM ITERATION... ',ITER
*.    D = <P!S>
        D = INPRDD(VEC1,VEC2,LU3,LU4,1,LBLK)
        C = RHO * GAMMA
        A = C/D
*.    R = R - A * S on LU2
        CALL VECSMD(VEC1,VEC2,ONE,-A,LU2,LU4,LU5,1,LBLK)
        CALL COPVCD(LU5,LU2,VEC1,1,LBLK)
*
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Residual on LU2 '
          CALL WRTVCD(VEC1,LU2,1,LBLK)
        END IF
*.    new residual has been obtained , check for convergence
        RNORM = INPRDD(VEC1,VEC2,LU2,LU2,1,LBLK)
        RNORME = MAX(RNORM,0.0D0)
        ERROR(ITER+1) = SQRT(RNORME)
*.    X = X + A * P
C?      WRITE(6,*) ' MICGCG : A = ', A
        CALL VECSMD(VEC1,VEC2,ONE,A,LU1,LU3,LU5,1,LBLK)
        CALL COPVCD(LU5,LU1,VEC1,1,LBLK)
*
        IF( SQRT(RNORME) .LT. TEST ) THEN
           CONVER = .TRUE.
           GOTO 1001
        ELSE
           CONVER = .FALSE.
*
* ============================
*. Prepare for next iteration
* ============================
*
*.H * R on LU4
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Diagonal and input to diagonal '
             CALL WRTVCD(VEC1,LUDIA,1,LBLK)
             CALL WRTVCD(VEC1,LU2  ,1,LBLK)
             WRITE(6,*) ' SHIFT = ', W
           END IF
*
           CALL DMTVCD(VEC1,VEC2,LUDIA,LU2,LU4,W,1,1,LBLK)
           IF(NTEST.GE.100) THEN
             WRITE(6,*) ' Preconditioner times residual '
             CALL WRTVCD(VEC1,LU4,1,LBLK)
           END IF
           IF(LUPROJ.NE.0) THEN
             OVLAP = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
             FACTOR = -OVLAP/X12
             CALL VECSMD(VEC1,VEC2,ONE,FACTOR,LU4,LUPROJ2,LU5,1,LBLK)
             CALL COPVCD(LU5,LU4,VEC1,1,LBLK)
             OVLAP2 = INPRDD(VEC1,VEC2,LUPROJ,LU4,1,LBLK)
C?           WRITE(6,*) ' Updated overlap of trial vector ', OVLAP2
*. Overlap between X and LUPROJ
             OVLAP3 = INPRDD(VEC1,VEC2,LUPROJ,LU1,1,LBLK)
             WRITE(6,*) ' Overlap between LU1 and LUPROJ ', OVLAP3
           END IF
*. GAMMA = <R!H!R>
           GAMMA = INPRDD(VEC1,VEC2,LU2,LU4,1,LBLK)
           B = GAMMA/C
*. P = RHO*(H*R + B*P) on LU3
           RHO = 1.0D0
           CALL VECSMD(VEC1,VEC2,ONE,B,LU4,LU3,LU5,1,LBLK)
           CALL COPVCD(LU5,LU3,VEC1,1,LBLK)
*.S = AP on LU4
           CALL MV8(VEC1,VEC2,LU3,LU4)
*.End of prepations for next iteration
        END IF
*
 1000 CONTINUE
 1001 CONTINUE
      IF(NTEST .GT. 0 ) THEN
*
      IF(CONVER) THEN
       WRITE(6,1010) NITER  ,ERROR(NITER+1)
 1010  FORMAT(1H0,'  convergence was obtained in...',I3,' iterations',/,
     +        1H ,'  norm of residual..............',E13.8)
      ELSE
       WRITE(6,1020) MAXIT ,ERROR(MAXIT+1)
 1020  FORMAT(1H0,' convergence was not obtained in',I3,'iterations',/,
     +        1H ,' norm of residual...............',E13.8)
      END IF
*
      END IF
*
      IF(NTEST.GE. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(1H0,' solution to set of linear equations')
       CALL WRTVCD(VEC1,LU1,1,LBLK)
C?     write(6,*) ' Matrix times solutiom through another cal to MV 8'
C?     CALL MV8(VEC1,VEC2,0,0)
C?     call WRTMT_LU(vec2,1,nvar,1,nvar)
      END IF
C
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(1H0,10X,'iteration point     norm of residual')
      DO 350 I=1,NITER+1
       II=I-1
       WRITE(6,1050)II,ERROR(I)
 1050  FORMAT(1H ,12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
*
      CALL QEXIT('MICGC')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MINDV4(MV7,
     &                  VEC1,VEC2,LU1,LU2,RNRM,EIG,FINEIG,MAXIT,NVAR,
     &                  LU3,LUDIA,NROOT,MAXVEC,NINVEC,
     &                  APROJ,AVEC,WORK,IPRT,
     &                  NPRDIM,H0,IPNTR,NP1,NP2,NQ,H0SCR,EIGSHF,
     &                  IOLSEN,IPICO)
*
* Davidson algorithm , requires two vectors in core
* Multi root version
*
* Allows updating of preconditioning matrix so this is
* the current eigenvector approximation
* is an eigenvector for the preconditioner
*
* Jeppe Olsen Sept 89
*             Jan  92 : MV7 entry
*
* Input :
* =======
*        MV7 : Name of routine performing matrix*vector calculation
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    complete vector
*        LU2,LU3   : Scatch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
* On input LU1 is supposed to hold initial guess to eigenvectors
*
* IOLSEN : Use inverse iteration modified Davidson
* IPICO  : Use perturbation estimate of new vector instead of
*          variational method
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION VEC1(*),VEC2(*)
       REAL * 8   INPROD
       DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
       DIMENSION APROJ(*),AVEC(*),WORK(*)
       DIMENSION H0(*),IPNTR(1)
       DIMENSION H0SCR(*)
*
* Dimensioning required of local vectors
*      APROJ  : MAXVEC*(MAXVEC+1)/2
*      AVEC   : MAXVEC ** 2
*      WORK   : MAXVEC*(MAXVEC+1)/2
*      H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
       DIMENSION FINEIG(1)
       LOGICAL CONVER,RTCNV(10)
*
       EXTERNAL MV7
*
       WRITE(6,*) ':::::::::::::::::::'
       WRITE(6,*) '  Entering MINDV4  '
       WRITE(6,*) ':::::::::::::::::::'
       IOLSTM = IOLSEN
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Variational'
       IF(IPRT.GT.1.AND.(IOLSEN.NE.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Inverse iteration modified Davidson, Perturbational'
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.EQ.0))
     & WRITE(6,*) ' Normal Davidson, Variational '
       IF(IPRT.GT.1.AND.(IOLSEN.EQ.0.AND.IPICO.NE.0))
     & WRITE(6,*) ' Normal Davidson, Perturbational'
       IF( MAXVEC .LT. 2 * NROOT ) THEN
         WRITE(6,*) ' SORRY MINDV2 WOUNDED , MAXVEC .LT. 2*NROOT '
         Call Abend2( ' ENFORCED STOP IN MINDV2' )
       END IF
*
       IF(IPICO.NE.0) THEN
         MAXVEC = 2*NROOT
       END IF
*
       KAPROJ = 1
       KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
       TEST = 1.0D-7
       CONVER = .FALSE.
       DO 1234 MACRO = 1,1
*
*.   INITAL ITERATION
       ITER = 1
       CALL REWINE( LU1 ,-1)
       CALL REWINE( LU2 ,-1)
       DO 10 IVEC = 1,NINVEC
         CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL MV7(VEC1,VEC2,0,0)
         CALL TODSC_LUCI(VEC2,NVAR,-1  ,LU2)
*        PROJECTED MATRIX
         CALL REWINE( LU1,-1)
         DO 8 JVEC = 1, IVEC
           CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
    8    CONTINUE
   10  CONTINUE
*
       IF( IPRT .GE.10 ) THEN
         WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
         CALL PRSYM(APROJ,NINVEC)
       END IF
*  DIAGONALIZE INITIAL PROJECTED MATRIX
       CALL COPVEC(APROJ,WORK(KAPROJ),NINVEC*(NINVEC+1)/2)
       CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NINVEC,0,1)
       DO 20 IROOT = 1, NROOT
         EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
   20  CONTINUE
*
       IF( IPRT  .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Initial set of eigenvalues '
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF),IROOT=1,NROOT)
       END IF
       NVEC = NINVEC
       IF (MAXIT .EQ. 1 ) GOTO  901
*
** LOOP OVER ITERATIONS
*
 1000 CONTINUE
      IF(IPRT  .GE. 10 ) THEN
       WRITE(6,*) ' INFO FORM ITERATION .... ', ITER
      END IF


        ITER = ITER + 1
*
** 1          NEW DIRECTION TO BE INCLUDED
*
*   1.1 : R = H*X - EIGAPR*X
       IADD = 0
       CONVER = .TRUE.
       DO 100 IROOT = 1, NROOT
         CALL SETVEC(VEC1,0.0D0,NVAR)
*
         CALL REWINE( LU2,-1)
         DO 60 IVEC = 1, NVEC
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           FACTOR = AVEC((IROOT-1)*NVEC+IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   60    CONTINUE
*
         EIGAPR = EIG(ITER-1,IROOT)
         CALL REWINE( LU1,-1)
         DO 50 IVEC = 1, NVEC
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
           FACTOR = -EIGAPR*AVEC((IROOT-1)*NVEC+ IVEC)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   50    CONTINUE
           IF ( IPRT  .GE.600 ) THEN
             WRITE(6,*) '  ( HX - EX ) '
             CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
           END IF
*  STRANGE PLACE TO TEST CONVERGENCE , BUT ....
         RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
         RNRM(ITER-1,IROOT) = RNORM
         IF(RNORM.LT. TEST ) THEN
            RTCNV(IROOT) = .TRUE.
         ELSE
            RTCNV(IROOT) = .FALSE.
            CONVER = .FALSE.
         END IF
         IF( ITER .GT. MAXIT) GOTO 100
*.  1.2 : MULTIPLY WITH INVERSE HESSIAN APROXIMATION TO GET NEW DIRECTIO
         IF( .NOT. RTCNV(IROOT) ) THEN
           IADD = IADD + 1
           CALL REWINE( LUDIA,-1)
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL H0M1TV(VEC2,VEC1,VEC1,NVAR,NPRDIM,IPNTR,
     &                 H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
           IF ( IPRT  .GE. 600) THEN
             WRITE(6,*) '  (D-E)-1 *( HX - EX ) '
             CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
           END IF
*
           IF(IOLSTM .NE. 0 ) THEN
* add Olsen correction if neccessary
              CALL REWINE(LU3,-1)
              CALL TODSC_LUCI(VEC1,NVAR,-1,LU3)
* Current eigen vector
              CALL REWINE( LU1,-1)
              CALL SETVEC(VEC1,0.0D0,NVAR)
              DO 59 IVEC = 1, NVEC
                CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
                FACTOR = AVEC((IROOT-1)*NVEC+ IVEC)
                CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
   59         CONTINUE
              IF ( IPRT  .GE. 600 ) THEN
                WRITE(6,*) ' And X  '
                CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
              END IF
              CALL TODSC_LUCI(VEC1,NVAR,-1,LU3)
* (H0 - E )-1  * X
              CALL REWINE( LUDIA,-1)
              CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
              CALL H0M1TV(VEC2,VEC1,VEC2,NVAR,NPRDIM,IPNTR,
     &                   H0,-EIGAPR,H0SCR,XDUMMY,NP1,NP2,NQ,
     &                 IPRT)
              CALL TODSC_LUCI(VEC2,NVAR,-1,LU3)
* Gamma = X(T) * (H0 - E) ** -1 * X
              GAMMA = INPROD(VEC2,VEC1,NVAR)
* is X an eigen vector for (H0 - 1 ) - 1
              CALL VECSUM(VEC2,VEC1,VEC2,GAMMA,-1.0D0,NVAR)
              VNORM = SQRT(MAX(0.0D0,INPROD(VEC2,VEC2,NVAR)))
              IF(VNORM .GT. 1.0D-7 ) THEN
                IOLSAC = 1
              ELSE
                IOLSAC = 0
              END IF
              IF(IOLSAC .EQ. 1 ) THEN
                IF(IPRT.GE.5) WRITE(6,*) ' Olsen Correction active '
                CALL REWINE(LU3,-1)
                CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU3,IMZERO,IAMPACK)
                DELTA = INPROD(VEC1,VEC2,NVAR)
                CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
                FACTOR = -DELTA/GAMMA
                IF(IPRT.GE.5) WRITE(6,*) ' DELTA,GAMMA,FACTOR'
                IF(IPRT.GE.5) WRITE(6,*)   DELTA,GAMMA,FACTOR
                CALL VECSUM(VEC1,VEC1,VEC2,FACTOR,1.0D0,NVAR)
                IF(IPRT.GE.600) THEN
                  WRITE(6,*) '  Modified new trial vector '
                  CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
                END IF
              ELSE
                IF(IPRT.GT.0) WRITE(6,*)
     &          ' Inverse correction switched of'
                CALL REWINE(LU3,-1)
                CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
              END IF
            END IF
*. 1.3 ORTHOGONALIZE TO ALL PREVIOUS VECTORS
           XNRMI =    INPROD(VEC1,VEC1,NVAR)
           CALL REWINE( LU1 ,-1)

           DO 80 IVEC = 1,NVEC+IADD-1
             CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
   80      CONTINUE
*. 1.4 Normalize vector and check for linear dependency
           SCALE = INPROD(VEC1,VEC1,NVAR)
           IF(ABS(SCALE)/XNRMI .LT. 1.0D-10) THEN
*. Linear dependency
             IADD = IADD - 1
             IF ( IPRT  .GE. 10 ) THEN
               WRITE(6,*) '  Trial vector linear dependent so OUT !!! '
             END IF
           ELSE
             C1NRM = SQRT(SCALE)
             FACTOR = 1.0D0/SQRT(SCALE)
             CALL SCALVE(VEC1,FACTOR,NVAR)
*
             CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU1)
             IF ( IPRT  .GE.600 ) THEN
               WRITE(6,*) 'ORTHONORMALIZED (D-E)-1 *( HX - EX ) '
               CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
             END IF
           END IF
*
         END IF
  100 CONTINUE
      IF( CONVER ) GOTO  901
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
**  2 : OPTIMAL COMBINATION OF NEW AND OLD DIRECTION
*
*  2.1: MULTIPLY NEW DIRECTION WITH MATRIX
       CALL REWINE( LU1,-1)
       CALL REWINE( LU2,-1)
       DO 110 IVEC = 1, NVEC
         CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU1,IMZERO,IAMPACK)
         CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
  110  CONTINUE
*
      DO 150 IVEC = 1, IADD
        CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        CALL MV7(VEC1,VEC2,0,0)
        CALL TODSC_LUCI(VEC2,NVAR,-1  ,LU2)
*   AUGMENT PROJECTED MATRIX
        CALL REWINE( LU1,-1)
        DO 140 JVEC = 1, NVEC+IVEC
          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
          CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
  140   CONTINUE
  150 CONTINUE
*  DIAGONALIZE PROJECTED MATRIX
      NVEC = NVEC + IADD
      CALL COPVEC(APROJ,WORK(KAPROJ),NVEC*(NVEC+1)/2)
      CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NVEC,0,1)
      IF(IPICO.NE.0) THEN
        E0VAR = WORK(KAPROJ)
        C0VAR = AVEC(1)
        C1VAR = AVEC(2)
*. overwrite with pert solution
        AVEC(1) = 1.0D0/SQRT(1.0D0+C1NRM**2)
        AVEC(2) = -C1NRM/SQRT(1.0D0+C1NRM**2)
        E0PERT = AVEC(1)**2*APROJ(1)
     &         + 2.0D0*AVEC(1)*AVEC(2)*APROJ(2)
     &         + AVEC(2)**2*APROJ(3)
        WORK(KAPROJ) = E0PERT
        WRITE(6,*) ' Var and Pert solution, energy and coefficients'
        WRITE(6,'(4X,3E15.7)') E0VAR,C0VAR,C1VAR
        WRITE(6,'(4X,3E15.7)') E0PERT,AVEC(1),AVEC(2)
      END IF
      DO 160 IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
 160  CONTINUE
*
       IF(IPRT .GE. 3 ) THEN
         WRITE(6,'(A,I4)') ' Eigenvalues of iteration ..', ITER
         WRITE(6,'(5F18.13)')
     &   ( (EIG(ITER,IROOT)+EIGSHF) ,IROOT=1,NROOT)
       END IF
*
      IF( IPRT  .GE. 5 ) THEN
        WRITE(6,*) ' PROJECTED MATRIX AND EIGEN PAIRS '
        CALL PRSYM(APROJ,NVEC)
        WRITE(6,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMT_LU(AVEC,NVEC,NROOT,MAXVEC,NROOT)
      END IF
*
**  PERHAPS RESET OR ASSEMBLE CONVERGED EIGENVECTORS
*
  901 CONTINUE
*
      IPULAY = 1
      IF(IPULAY.EQ.1 .AND. MAXVEC.EQ.3 .AND.NVEC.GE.2.
     &   .AND. .NOT.CONVER) THEN
* Save trial vectors : 1 -- current trial vector
*                      2 -- previous trial vector orthogonalized
        CALL REWINE( LU3,-1)
        CALL REWINE( LU1,-1)
*. Current trial vector
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2200 IVEC = 1, NVEC
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
         CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2200   CONTINUE
        SCALE = INPROD(VEC1,VEC1,NVAR)
        SCALE  = 1.0D0/SQRT(SCALE)
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU3)
* Previous trial vector orthonormalized
        CALL REWINE(LU1,-1)
        CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU1,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        SCALE2 = INPROD(VEC2,VEC2,NVAR)
        SCALE2 = 1.0D0/SQRT(SCALE2)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC_LUCI(VEC2,NVAR,-1,LU3)
*
        CALL REWINE( LU1,-1)
        CALL REWINE( LU3,-1)
        DO 2411 IVEC = 1,2
          CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC1,NVAR,-1,  LU1)
 2411   CONTINUE
*. Corresponding sigma vectors
        CALL REWINE ( LU3,-1)
        CALL REWINE( LU2,-1)
        CALL SETVEC(VEC1,0.0D0,NVAR)
        DO 2250 IVEC = 1, NVEC
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
          FACTOR =  AVEC(IVEC)
          CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
 2250   CONTINUE
*
        CALL SCALVE(VEC1,SCALE,NVAR)
        CALL TODSC_LUCI(VEC1,NVAR,-1,  LU3)
* Sigma vector corresponding to second vector on LU1
        CALL REWINE(LU2,-1)
        CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL SCALVE(VEC2,SCALE2,NVAR)
        CALL TODSC_LUCI(VEC2,NVAR,-1,LU3)
*
        CALL REWINE( LU2,-1)
        CALL REWINE( LU3,-1)
        DO 2400 IVEC = 1,2
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC2,NVAR,-1  ,LU2)
 2400   CONTINUE
        NVEC = 2
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 2410 IROOT = 1,NVEC
          AVEC((IROOT-1)*NVEC+IROOT) = 1.0D0
 2410   CONTINUE
*.Projected hamiltonian
       CALL REWINE( LU1 ,-1)
       DO 2010 IVEC = 1,NVEC
         CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
         CALL REWINE( LU2,-1)
         DO 2008 JVEC = 1, IVEC
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           IJ = IVEC*(IVEC-1)/2 + JVEC
           APROJ(IJ) = INPROD(VEC1,VEC2,NVAR)
 2008    CONTINUE
 2010  CONTINUE
      END IF
      IF(NVEC+NROOT.GT.MAXVEC .OR. CONVER .OR. MAXIT .EQ.ITER)THEN
        CALL REWINE( LU3,-1)
        DO 320 IROOT = 1, NROOT
          CALL REWINE( LU1,-1)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 200 IVEC = 1, NVEC
            CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  200     CONTINUE
*
          SCALE = INPROD(VEC1,VEC1,NVAR)
          SCALE  = 1.0D0/SQRT(SCALE)
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU3)
  320   CONTINUE
        CALL REWINE( LU1,-1)
        CALL REWINE( LU3,-1)
        DO 411 IVEC = 1,NROOT
          CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC1,NVAR,-1,  LU1)
  411   CONTINUE
* CORRESPONDING SIGMA VECTOR
        CALL REWINE ( LU3,-1)
        DO 329 IROOT = 1, NROOT
          CALL REWINE( LU2,-1)
          CALL SETVEC(VEC1,0.0D0,NVAR)
          DO 250 IVEC = 1, NVEC
            CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
            FACTOR =  AVEC((IROOT-1)*NVEC+IVEC)
            CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FACTOR,NVAR)
  250     CONTINUE
*
          CALL SCALVE(VEC1,SCALE,NVAR)
          CALL TODSC_LUCI(VEC1,NVAR,-1,  LU3)
  329   CONTINUE
* PLACE C IN LU1 AND HC IN LU2
        CALL REWINE( LU2,-1)
        CALL REWINE( LU3,-1)
        DO 400 IVEC = 1,NROOT
          CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
          CALL todsc_luci(VEC2,NVAR,-1  ,LU2)
  400   CONTINUE
        NVEC = NROOT
*
        CALL SETVEC(AVEC,0.0D0,NVEC**2)
        DO 410 IROOT = 1,NROOT
          AVEC((IROOT-1)*NROOT+IROOT) = 1.0D0
  410   CONTINUE
*
        CALL SETVEC(APROJ,0.0D0,NVEC*(NVEC+1)/2)
        DO 420 IROOT = 1, NROOT
          APROJ(IROOT*(IROOT+1)/2 ) = EIG(ITER,IROOT)
  420   CONTINUE
C
      END IF
C
C     IF( ITER .LT. MAXIT .AND. .NOT. CONVER) GOTO 1000
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
 1001 CONTINUE
*. Place first eigenvector in vec1
      CALL REWINE(LU1,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)



* ( End of loop over iterations )
*
*
*
      IF( .NOT. CONVER ) THEN
*        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(6,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
*        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(6,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
*. Final eigenvalues
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
 1601   CONTINUE
*
      IF ( IPRT .GT. 1 ) THEN
        DO 1600 IROOT = 1, NROOT
          WRITE(6,*)
          WRITE(6,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(6,*)
     &    '============================================'
          WRITE(6,*)
          WRITE(6,1190) FINEIG(IROOT)
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(6,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL REWINE( LU1,-1)
            CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
            CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
          END IF
          WRITE(6,1300)
 1300     FORMAT(1H0,' Summary of iterations ',/,1H
     +          ,' ----------------------')
          WRITE(6,1310)
 1310     FORMAT
     &    (1H0,' Iteration point        Eigenvalue         Residual ')
          DO 1330 I=1,ITER
 1330     WRITE(6,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      END IF
*
      IF(IPRT .EQ. 1 ) THEN
        DO 1607 IROOT = 1, NROOT
          WRITE(6,'(A,2I3,F14.7,2E10.3)')
     &    ' >>> CI-OPT Iter Root E g-norm g-red',
     &                 ITER,IROOT,FINEIG(IROOT),
     &                 RNRM(ITER,IROOT),
     &                 RNRM(1,IROOT)/RNRM(ITER,IROOT)
 1607   CONTINUE
      END IF
 1234 CONTINUE
C
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MINGCG(MV8,LU1,LU2,LU3,LUDIA,VEC1,VEC2,
     &                  MAXIT,CONVER,TEST,W,ERROR,NVAR,
     &                  LUPROJ,IPRT)
*
* Solve set of linear equations
*
*             AX = B
*
* with preconditioned conjugate gradient method for
* case where two complete vectors can be stored in core
*
* Initial appriximation to solution must reside on LU1
* LU2 must contain B.All files are  overwritten
*
*
* Final solution vector is stored in LU1
* A scalar w can be added to the diagonal of the preconditioner
*
* If LUPROJ .NE. 0 , the optimization subspace is restricted to be orthogonal
* to the first vector in LUPROJ.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*),ERROR(MAXIT+1)
      REAL*8 INPROD
      LOGICAL CONVER
*
      EXTERNAL MV8
*
      CONVER = .FALSE.
      ITER = 1
      NTEST = 0
      NTEST = MAX(NTEST,IPRT)
*
* =============
* Initial point
* =============
*
*.R = B - (A)*X
      CALL REWINE(LU1,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINE(LU2,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
      CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-1.0D0,NVAR)
*
      RNORM = SQRT( INPROD(VEC1,VEC1,NVAR) )
      ERROR(1) = RNORM
      CALL REWINE(LU2,-1)
      CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU2)
*. Preconditioner H times initial vector , H * R
*.H * R
      CALL REWINE(LUDIA,-1)
      CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
      CALL DIAVC2(VEC2,VEC1,VEC2,W,NVAR)
      IF(LUPROJ.NE.0) THEN
        CALL REWINE(LUPROJ,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
        OVLAP = INPROD(VEC1,VEC2,NVAR)
        CALL VECSUM(VEC2,VEC2,VEC1,1.0D0,-OVLAP,NVAR)
        CALL REWINE(LU2,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU2,IMZERO,IAMPACK)
      END IF
*. GAMMA = <R!H!R>
      GAMMA = INPROD(VEC1,VEC2,NVAR)
*. P = RHO * H*R
      RHO = 1.0D0
      CALL SCALVE(VEC2,RHO,NVAR)
      CALL REWINE(LU3,-1)
      CALL TODSC_LUCI(VEC2,NVAR,-1  ,LU3)
      CALL COPVEC(VEC2,VEC1,NVAR)
*.S = AP
      CALL MV8(VEC1,VEC2,0,0)
      CALL REWINE (LU3,-1)
      CALL FRMDSC_LUCI(VEC1,NVAR,-1,LU3,IMZERO,IAMPACK)
*
* ====================
* Loop over iterations
* ====================
*
      NITER = 0
      DO 1000 ITER = 1, MAXIT
*.    P is assumed in VEC1 and S = A*P in VEC2

        NITER = NITER + 1
       IF ( NTEST .GE. 10 )
     & WRITE(6,*) ' INFORMATION FROM ITERATION... ',ITER
*.    D = <P!S>
        D = INPROD(VEC1,VEC2,NVAR)
        C = RHO * GAMMA
        A = C/D
*.    R = R - A * S
        CALL REWINE(LU2,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU2,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-A,NVAR)
        CALL REWINE(LU2,-1)
        CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU2)
*.    new residual has been obtained , check for convergence
        RNORM = INPROD(VEC1,VEC1,NVAR)
        ERROR(ITER+1) = SQRT(RNORM)
*.    X = X + A * P
        CALL REWINE(LU1,-1)
        CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU1,IMZERO,IAMPACK)
        CALL REWINE(LU3,-1)
        CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
        CALL VECSUM(VEC1,VEC2,VEC1,1.0D0,A,NVAR)
        CALL REWINE(LU1,-1)
        CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU1)
*
        IF( SQRT(RNORM) .LT. TEST ) THEN
           CONVER = .TRUE.
           GOTO 1001
        ELSE
           CONVER = .FALSE.
*
* ============================
*. Prepare for next iteration
* ============================
*
*.       H * R
           CALL REWINE(LU2,-1)
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU2,IMZERO,IAMPACK)
           CALL REWINE(LUDIA,-1)
           CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LUDIA,IMZERO,IAMPACK)
           CALL DIAVC2(VEC1,VEC2,VEC1 ,W,NVAR)
           IF(LUPROJ.NE.0) THEN
             CALL REWINE(LUPROJ,-1)
             CALL FRMDSC_LUCI(VEC2,NVAR,-1,LUPROJ,IMZERO,IAMPACK)
             OVLAP = INPROD(VEC1,VEC2,NVAR)
             CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,-OVLAP,NVAR)
             CALL REWINE(LU2,-1)
             CALL FRMDSC_LUCI(VEC2,NVAR,-1,LU2,IMZERO,IAMPACK)
           END IF
           GAMMA = INPROD(VEC1,VEC2,NVAR)
           B = GAMMA/C
*.       P = RHO*(H*R + B*P)
           CALL REWINE(LU3,-1)
           CALL FRMDSC_LUCI(VEC2,NVAR,-1  ,LU3,IMZERO,IAMPACK)
           CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,B,NVAR)
*.       Define next RHO
           RHO = 1.0D0
           CALL SCALVE(VEC1,RHO,NVAR)
           CALL REWINE(LU3,-1)
           CALL TODSC_LUCI(VEC1,NVAR,-1  ,LU3)
*.       S = MATRIX * P
           CALL MV8(VEC1,VEC2,0,0)
           CALL REWINE(LU3,-1)
           CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU3,IMZERO,IAMPACK)
*.End of prepations for next iteration
        END IF
*
*
 1000 CONTINUE
 1001 CONTINUE
      IF(NTEST .GT. 0 ) THEN
      IF(CONVER) THEN
       WRITE(6,1010) NITER  ,ERROR(NITER+1)
 1010  FORMAT(1H0,'  convergence was obtained in...',I3,' iterations',/,
     +        1H ,'  norm of residual..............',F13.8)
      ELSE
       WRITE(6,1020) MAXIT ,ERROR(MAXIT +1 )
 1020  FORMAT(1H0,' convergence was not obtained in',I3,'iterations',/,
     +        1H ,' norm of residual...............',F13.8)
      END IF
      END IF
C
      IF(NTEST.GT. 50 ) THEN
       WRITE(6,1025)
 1025  FORMAT(1H0,' solution to set of linear equations')
       CALL REWINE(LU1,-1)
       CALL FRMDSC_LUCI(VEC1,NVAR,-1  ,LU1,IMZERO,IAMPACK)
       CALL WRTMT_LU(VEC1,1,NVAR,1,NVAR)
C?     write(6,*) ' Matrix times solutiom through another cal to MV 8'
C?     CALL MV8(VEC1,VEC2,0,0)
C?     call WRTMT_LU(vec2,1,nvar,1,nvar)
      END IF
C
      IF(NTEST.GT.0) THEN
      WRITE(6,1040)
 1040 FORMAT(1H0,10X,'iteration point     norm of residual')
      DO 350 I=1,NITER+1
       II=I-1
       WRITE(6,1050)II,ERROR(I)
 1050  FORMAT(1H ,12X,I5,13X,E15.8)
  350 CONTINUE
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MIXHONE(H1,H2,NREPTP,IREPTP,NOBTP,NSMOB)
*
* Replace selected type blocks of H1 with the corresponding blocks
* in H2
*
*. H1 and H2 are assumed to be in symmetry order !
*. -and total symmetric
*
*     Jeppe Olsen, March 14 1996 ( Still snowing in Lund )
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "orbinp.inc"
*. Specific input
      DIMENSION IREPTP(*)
      DIMENSION H2(*)
*. Input and output
      DIMENSION H1(*)
*
      DO ISMOB = 1, NSMOB
        IF (ISMOB.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISMOB-1)*(NTOOBS(ISMOB-1)+1)/2
        END IF
*. Loop over types for given symmetry
        DO ITP = 1, NOBTP
          IF(ITP.EQ.1) THEN
           IOBOFF = 1
          ELSE
            IOBOFF = IOBOFF + NOBPTS(ITP-1,ISMOB)
          END IF
          DO JTP = 1, ITP
            IF(JTP.EQ.1) THEN
             JOBOFF = 1
            ELSE
              JOBOFF = JOBOFF + NOBPTS(JTP-1,ISMOB)
            END IF
*. Number of elements in this type-type block
            LIOB = NOBPTS(ITP,ISMOB)
            LJOB = NOBPTS(JTP,ISMOB)
*
*. Should this block of H1 be replaced by corresponding block of H2
            IF(ITP.EQ.JTP) THEN
              IMOVE = 0
              DO KTP = 1, NREPTP
                IF(IREPTP(KTP).EQ.ITP) IMOVE = 1
              END DO
*
              IF(IMOVE.EQ.1) THEN
C?              WRITE(6,*) ' Block transfer ISMOB ITP JTP ',
C?   &          ISMOB,ITP,JTP
                DO IOB = IOBOFF,IOBOFF+LIOB-1
                  DO JOB = JOBOFF, IOB
                    H1(IOFF-1+IOB*(IOB-1)/2+JOB)
     &            = H2(IOFF-1+IOB*(IOB-1)/2+JOB)
                  END DO
                END DO
              END IF
*
            END IF
          END DO
        END DO
      END DO
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' MIXHONE in action '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' NSMOB NOBTP ', NSMOB,NOBTP
        WRITE(6,*) ' Types to be changed '
        CALL IWRTMA(IREPTP,1,NREPTP,1,NREPTP)
        WRITE(6,*) ' output H1 and H2 '
C       APRBLM2(A,LROW,LCOL,NBLK,ISYM)
        CALL APRBLM2(H1,NTOOBS,NTOOBS,NSMOB,1)
        CALL APRBLM2(H2,NTOOBS,NTOOBS,NSMOB,1)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ONEEL_MAT_DISC(H,IHSM,NSM,NRPSM,NCPSM,LUH,IFT)
*
* Transfer one-electron matrix H between memory and disc file in
* LUCIA format
*
* IFT = 1 => From disc ( read)
* IFT = 2 => To   disc (write)
*
* Note : File LUH is supposed to be at start of correct integral block
*
* Jeppe Olsen, Feb. 98
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER NRPSM(NSM),NCPSM(*)
      DIMENSION H(*)
*
*. Order of integrals are
*
*     Loop over Symmetry of row index => Symmetry of column  index
*      Loop over columns in symmetry block
*        Loop over rows in symmetry block
*        End of loop over rows in symmetry block
*      End of Loop over columns in symmetry block
*     End of loop over symmetry of row index
*
* Each symmetry block is thus given in complete form
* Note all integrals are in a single record
*
*. Length of list
C              NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,IPACK)
      LENGTH = NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,0)
*. and read/write
      WRITE(6,*) ' ONEEL, LUH = ', LUH
      IF(IFT.EQ.1) THEN
        DO IJ = 1, LENGTH
          READ(LUH,'(E22.15)') H(IJ)
        END DO
      ELSE IF (IFT.EQ.2) THEN
        DO IJ = 1, LENGTH
          WRITE(LUH,'(E22.15)') H(IJ)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE PERTCTL(ISM,ISPC,EREF)
*
* Master routine for perturbation calculations
* (Largest allowed order : 1000)
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLSMAT, KLVMAT, KLH0MAT, KLEN, KLSCR, KVEC1, KVEC2, 
     &          KLCIOIO
!               for addressing of WORK
      EXTERNAL MV7
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "gasstr.inc"

*. And defining perturbation operator
#include "oper.inc"
*
#include "csfbas.inc"
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      COMMON/SPINFO/MULTSP,MS2P,
     &              MINOP,MAXOP,NTYP,NDPCNT(MXPCTP),NCPCNT(MXPCTP),
     &              NCNATS(MXPCTP,MXPCSM),NDTASM(MXPCSM),NCSASM(MXPCSM),
     &              NCNASM(MXPCSM)
      COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H,ECORE_HEX
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      DIMENSION EN(1000)
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PERTCT')


*
      WRITE(6,*) '**************************************'
      WRITE(6,*) '*                                    *'
      WRITE(6,*) '*   Perturbation calculation         *'
      WRITE(6,*) '*                                    *'
      WRITE(6,*) '**************************************'
*
      WRITE(6,*)
      WRITE(6,*) '  Largest order of correction vector ',NPERT
      WRITE(6,*)
C     IF(IPART.EQ.1) THEN
C       WRITE(6,*) ' Moller-Plesset Partitioning'
C     ELSE IF(IPART.EQ.2) THEN
C       WRITE(6,*) ' Epstein-Nesbet Partitioning'
C     ELSE IF (IPART.EQ.3) THEN
C       WRITE(6,*) ' One-body Hamiltonian read in '
C     END IF
      IF(IE0AVEX.EQ.1) THEN
          WRITE(6,*)
     &  ' expectation value of H0 used as zero order energy '
      ELSE IF( IE0AVEX.EQ.2) THEN
          WRITE(6,*)
     &  ' Exact energy of reference used as zero order energy'
      ELSE IF( IE0AVEX.EQ.3) THEN
          WRITE(6,*)
     &  ' Read in energy is used as zero order energy'
      END IF
      WRITE(6,*) ' Root choosen as zero order state ', IRFROOT
*
*
*. 0 : Initialization
*
      IF(NOCSF.EQ.0) THEN
        WRITE(6,*) ' Please turn off csf''s '
        STOP'NO CSF''s in PERTCTL !! '
      END IF
*
      NTEST = 10
      IPRNT = NTEST
      NTEST = MAX(NTEST,IPRNT)
      NDET = XISPSM(ISM,ISPC)
      NEL = NELCI(ISPC)
      WRITE(6,*) ' ISM ISPC ', ISM,ISPC
      WRITE(6,*) ' Number of determinants in internal space ',NDET
*.Transfer to CANDS
      ICSM = ISM
      ISSM = ISM
      ICSPC = ISPC
      ISSPC = ISPC
      WRITE(6,*) ' PERTCTL : ICSPC ISSSPC : ', ICSPC,ISSPC
      NVAR = NDET
      IF(IPRNT.GE.5)
     &WRITE(6,*) '  NVAR in REFCI ', NVAR
*. Arrays for S, V, H0 over correction vectors
      LENNY = (NPERT+1)*(NPERT+3)/2
      WRITE(6,*) ' LENNY ', LENNY
      CALL MEMMAN(KLSMAT ,LENNY,'ADDL  ',2,'LSMAT ')
      CALL MEMMAN(KLVMAT ,LENNY,'ADDL  ',2,'VSMAT ')
      CALL MEMMAN(KLH0MAT,LENNY,'ADDL  ',2,'H0SMAT')
*. Energy correction and scratch vector
      CALL MEMMAN(KLEN ,2*NPERT+2,'ADDL  ',2,'EN    ' )
      CALL MEMMAN(KLSCR,NPERT+1,'ADDL  ',2,'SCR   ' )


      CALL MEMMAN(KVEC1,LBLOCK,'ADDL  ',2,'VEC1  ')
*. Vec2 will also be used as scratch in explicit hamiltonian generation
*. for CSF's
      IF(NOCSF.NE.0) THEN
        LBLOC2 = LBLOCK
      ELSE
        CALL LCNHCN(LSCR)
        LBLOC2 = MAX(LBLOCK,LSCR)
      END IF
      CALL MEMMAN(KVEC2,LBLOC2,'ADDL  ',2,'VEC2  ')
*. Sblock is used in general nowadays so, allocate an extra block
      I_USE_SBLOCK=1
      IF(I_USE_SBLOCK.EQ.1) THEN
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*
        NOCTPA = NOCTYP(IATP)
        NOCTPB = NOCTYP(IBTP)
*
        CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL IAIBCM(ISPC,WORK(KLCIOIO))
        CALL MXRESC(WORK(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &              NSMST,NSTFSMSPGP,MXPNSMST,
     &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &              NELFSPGP,
     &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
        IF(IPRCIX.GE.2) THEN
          WRITE(6,*) 'PERTCT : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                         MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
           WRITE(6,*) 'PERTCT : MXADKBLK ', MXADKBLK
        END IF
        LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
        IF(IPRCIX.GE.2)
     &  WRITE(6,*) ' Space for resolution matrices ',LSCR2
        LSCR12 = MAX(LBLOCK,2*LSCR2)
        CALL MEMMAN(KVEC3,LSCR12,'ADDL  ',2,'KC2   ')
      END IF
*
*. 1 : Construct zero order operator : FI + FA
*
*. Copy root defining zero order operator to first vectors
      CALL REWINE(LUHC,-1)
      CALL REWINE(LUC,-1)
      WRITE(6,*) ' Root used to define Zero order op ', IH0ROOT
      DO JROOT = 1, IH0ROOT
        CALL REWINE(LUSC36,-1)
        CALL COPVCD(LUC,LUSC36,WORK(KVEC1),0,-1)
      END DO
      CALL COPVCD(LUSC36,LUHC,WORK(KVEC1),1,-1)
*. Construct corresponding one-body density matrix
COLD  CALL DENSI1(WORK(KRHO1),WORK(KVEC1),WORK(KVEC2),LUHC,LUSC36,1)
C     KRHO2 = 1
      CALL DENSI2(1,WORK(KRHO1),WORK(KRHO2),WORK(KVEC1),WORK(KVEC2),
     &     LUHC,LUSC36,EXPS2)
*
*. Initialize with proper zero order root
*
      WRITE(6,*) ' After DENSI2, LUC LUSC36',LUC,LUSC36
      CALL REWINE(LUC,-1)
      DO JROOT = 1, IRFROOT
        CALL REWINE(LUSC36,-1)
        CALL COPVCD(LUC,LUSC36,WORK(KVEC1),0,-1)
      END DO
      CALL COPVCD(LUSC36,LUC,WORK(KVEC1),1,-1)
*
      LU18 = 18
*. Will MP operator be invoked
      IUSEMP = 0
      DO ISPC = 1, NPTSPC
        IF(IH0INSPC(ISPC).EQ.1.OR.IH0INSPC(ISPC).EQ.3
     &     .OR.IH0INSPC(ISPC).EQ.5) IUSEMP = 1
      END DO
*
      IF(IUSEMP.EQ.1) THEN
        WRITE(6,*) ' Moller-Plesset operator will be used '
      ELSE
        WRITE(6,*) ' Moller-Plesset operator will not be used '
      END IF
*. Will operator be read in ( not clean : Use IPART = 3 in PERTU
* keyword
      IF(IUSEMP.EQ.1) THEN
        CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
        CALL FIFAM(WORK(KFI))
        CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
        ECORE_H = 0.0D0
        IF(IUSE_PH.EQ.1) THEN
         CALL FI(WORK(KFI),ECORE_H,0)
        END IF
*. Should a part of original one electron operator be
*  copied ( For mix exact Hamiltonian/Fock arroaches )
        IF(NH0EXSPC.NE.0) THEN
C             MIXHONE(H1,H2,NSMOB,NREPTP,IREPTP,NSMOB,NOBTP)
         CALL MIXHONE(WORK(KFI),WORK(KINT1),NH0EXSPC,IH0EXSPC,NGAS,
     &                NSMOB)
        END IF
      ELSE
        CALL COPVEC(WORK(KINT1),WORK(KFI),NINT1)
        CALL COPVEC(WORK(KINT1O),WORK(KFIO),NINT1)
      END IF
      IF (IPART.EQ.3) THEN
*. Read in from file 18
        REWIND (LU18)
        CALL FRMDSC_LUCI(WORK(KFI),NINT1,-1,LU18,IMZERO,IAMPACK)
        CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
        ECORE_H = 0
        IF(IUSE_PH.EQ.1) THEN
         CALL FI(WORK(KFI),ECORE_H,0)
        END IF
        WRITE(6,*) ' H0 read in from LU18 '
        CALL APRBLM2(WORK(KFI),NTOOBS,NTOOBS,NSMOB,ISM)
*. Continue as mormal MP a piece of dirty code can never harm
        IPART = 1
        MPORENP = 1
      END IF
*. Save H0 for future generations
      REWIND  LU18
      CALL TODSC_LUCI(WORK(KFI),NINT1,-1,LU18)
      REWIND LU18

*. No explicit construction of diagonal
      IDIDIA = 1
      IF(IDIDIA.EQ.0) THEN
*
*. 2 : Diagonal with FI + FA
*
*. swap H and FI + FA
        IF(IPART.EQ.1) THEN
          CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
          CALL SWAPVE(WORK(KFIO),WORK(KINT1O),NINT1)
        END IF
        IF(ICISTR.GE.2) CALL REWINE(LUDIA,-1)
*. Transfer to COPER
        IPERTOP = 1
        IF(IPART.EQ.1) THEN
          I12 = 1
        ELSE
          I12 = 2
        END IF
        ECOREP = ECORE_H
* FIXME !!!
        call quit('Wrong invocation of GASDIAT (8 arguments needed)')
C       CALL GASDIAT(WORK(KVEC1),LUDIA,ECOREP,ICISTR,I12)
* FIXME !!!
        IF(NOCSF.EQ.1.AND.ICISTR.EQ.1) THEN
          CALL REWINE(LUDIA,-1)
          CALL TODSC_LUCI(WORK(KVEC1),NVAR,-1,LUDIA)
C       ELSE IF(ICISTR.EQ.1.AND.NOCSF.EQ.0) THEN
C         CALL CSDIAG(WORK(KVEC2),WORK(KVEC1),NCNATS(1,ISM),NTYP,
C    &                WORK(KICTS(1)),NDPCNT,NCPCNT,0,
C    &                0,IDUM,IPRNT)
C         CALL REWINE(LUDIA,-1)
C         CALL TODSC_LUCI(WORK(KVEC2),NVAR,-1,LUDIA)
*. For transfer to H0CSF
C         CALL COPVEC(WORK(KVEC2),WORK(KVEC1),NVAR)
        END IF
*. swap H and FI + FA to get things in right place !
        IF(IPART.EQ.1)  THEN
          CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
          CALL SWAPVE(WORK(KFIO),WORK(KINT1O),NINT1)
        END IF
      END IF
*
* Transfer control to perturbation iterater
*
*. IS there a pert of Hamiltonian that is no diagonal
*. (requires solution of linear equations )
      IH0DIA = 1
      DO ISPC = 1, NPTSPC
        IF(IH0INSPC(ISPC).EQ.3.OR.IH0INSPC(ISPC).EQ.4.OR.
     &     IH0INSPC(ISPC).EQ.5) IH0DIA=0
      END DO
*
      IF(IH0DIA.EQ.0) THEN
        WRITE(6,*) ' Nondiagonal Approximate Hamiltonian '
      ELSE
        WRITE(6,*) ' Diagonal approximate Hamiltonian '
      END IF

*
* Nondiagonal form of perturbations : Currently indicated by
* operator type 3 and 4
      IF(ICISTR.EQ.1) THEN
        LBLK = NVAR
      ELSE
        LBLK = - 1
      END IF
*
*. Transfer to COPER
*.  Perturbation matrix
       IPERTOP = 1
       IF(IPART.EQ.1) THEN
        I12 = 1
       ELSE
        I12 = 2
       END IF
*
      IF(IE0AVEX.EQ.3) THEN
        EREF = E0READ-ECORE
        WRITE(6,*) ' Zero order energy read in - ECORE ',EREF
      END IF
      CALL SIMPRT(LUC,LUSC36,LUHC,WORK(KLEN),WORK(KLSCR),
     &            NPERT,WORK(KVEC1),WORK(KVEC2),
     &            LUSC1,LUSC2,LBLK,IH0DIA,LUDIA,WORK(KLSMAT),
     &            WORK(KLVMAT),WORK(KLH0MAT),ECORE,ECORE_H,
     &            ECORE_HEX,EREF,IE0AVEX,LUSC39)
*. Analyze space spanned by zero order state and correction vectors
      CALL PERT_SUBSPACE(NPERT,WORK(KLH0MAT),
     &     WORK(KLVMAT),WORK(KLSMAT),ECORE)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PERTCT')
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SIMPRT(LURF,LUN,LUVN,
     &           EN,SCR,MAXORD,VEC1,VEC2,LU1,LU2,
     &           LBLK,IH0DIA,LUH0,S,V,H0,ECORE,ECORE_H,ECORE_HEX,
     &           EREF,IE0AVEX,LUHI0)
*
* Solve the perturbation equations
*
* E(n) = <0!V!0(n-1)>
*
* !0(n)> = (H-E)-1( E(n)!0> -V!0(n-1)>
*                  + sum(l=1,n-1) E(l)!0(n-l)> )
*
* Simplified version for total symmetric perturbation
*
* Alternative expressions for the energy corrections
* are invoked using the 2n+1 rule
*
* E(n+1) =  <0!v!n>
*
*        = <k!v!n-k> - sum(j=0,k-1)sum(m=k-j,n-j)<j+1!n-m-j>E(m)
*
* A note on the perturbation :
*
* The perturbation is in general of the form
*
* H0 = QH(apr)Q + E0P, P = |0><0>, Q = 1-P
*
* Where Q is some approcimation to the hamiltonian.
*
* In order to calculate ((H0-alpha)** -1 |x>, where <x|0>=0 , one must
* distinguish between two cases
*
* =====================================
* 1 |0> is an eigenfunction for H(apr)
* =====================================
*
* In this case (H0-alpha)** -1 |x> = (H(apr)-alpha)** -1 |x>
*
* =========================================
* 2 |0> is  not an eigenfunction for H(apr)
* =========================================
*
*               (H0-alpha)** -1 |x> = (H(apr)-alpha)** -1 |x>
*             - (H(apr)-alpha)** -1 |0> <x|(H(apr)-alpha)** -1|0>
*                                       -------------------------
*                                       <X|(H(apr)-alpha)** -1|0>
*
*
* Input
* =====
* LURF : file containing reference vector
* LUN : file number for file to contain perturbation vectors
* LUN : file number for file to contain perturbation vectors
* MAXORD : Order through which the equations should be solved
* VEC1,VEC2 : Scratch vectors ,omplete or blocks of vectorS
* LU1, LU2, : scratch files
* S V : vectors of size MAXORD*(MAXORD-1)/2
* Ecore : Core energy
* Eref  : Exact energy of reference state
* IE0AVEX : choice of zero order energy : 1 => E0 = <0|H0|0>
*                                         2 => E0 = EREF ( as supplied )
*
*
* Output
* ======
* LUN : contains the MAXORD correction vectors
* LUVN : contains the perturbation times the last correction vector
* EN : Contains the energy corrections through order 2*MAXORD+1
*
* Internal links
* ===============
*
* Solutions of linear eqs : HINTV
* Hamiltonian times vector: MV7
* H0 times vector         : H0TVM
*
*. Please do not go beyond perturbation level 100
* Jeppe Olsen ,  Summer of 94
*                Winter 96 : Nondiagonal H(apr), general H(apr),
*                option for diagonal on disc eliminated
*                (only direct calculation allowed now)
*
*                Winter of 99 : (H0-1)** -1 |0> on LUHI0 added
      IMPLICIT REAL*8(A-H,O-Z)
      REAL * 8 INPRDD
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION EN(*)
      DIMENSION S(*),V(*),H0(*)
*. For communicating with H0TVM
      COMMON/CENOT/E0
      COMMON/CSHIFT/SHIFT,IPROJ
*. For communicating with MV7
#include "oper.inc"
*. A bit of  scratch
      DIMENSION SCR(*)
*.
      DIMENSION XTEST(10000)
*
      NTEST = 5
*. Use direct diagonal routines
      IDIDIA = 1

      ONE = 1.0D0
      ONEM = -1.0D0
      ZERO = 0.0D0
      CALL SETVEC(V,ZERO,(MAXORD+1)*(MAXORD+1+1)/2)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to SIMPRT '
        WRITE(6,*) ' LBLK = ', LBLK
        WRITE(6,*) ' ECORE_H = ', ECORE_H
        WRITE(6,*) ' ECORE_HEX = ', ECORE_HEX
      END IF
*
      IF(NTEST.GE.600) THEN
        WRITE(6,*) ' initial reference '
        CALL WRTVCD(VEC1,LURF,1,LBLK)
      END IF
*
* ===============================================================
* 1 :                   Initialization
* ===============================================================
*
*
        IF(IE0AVEX.GE.2) THEN
          E0RF = EREF
        ELSE
*  ===============
*. E0RF = <0!H(apr)!0>
*  ===============
*
          IF(IH0DIA.NE.0) THEN
*. Diagonal H0, simple
            CALL REWINE(LU1,-1)
            CALL REWINE(LURF,-1)
            CALL DIA0TRM_GAS(1,LURF,LU1,VEC1,VEC2,0.0D0)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' LU1 according to DIATRM '
              CALL WRTVCD(VEC1,LU1,1,LBLK)
            END IF
          ELSE
*. multiply with H(apr)
            E0 = 0.0D0
            IPROJ = 0
            SHIFT = ECORE_H
            IPERTOP = 1
CJAN25      CALL MV7(VEC1,VEC2,LURF,LU1)
            CALL H0TVM(VEC1,VEC2,LURF,LU1)
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' H(apr) times vector '
              CALL WRTVCD(VEC1,LU1,1,LBLK)
            END IF
          END IF
          E0RF = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
        END IF
*
        WRITE(6,*) ' E0RF = ', E0RF
        ENERGY = E0RF
*. Check of |0> is an eigenfunction for H(apr) ( not H0 ! )
* Calculate H(apr)|0> - <0!H(apr)!0> |0>
        IPERTOP = 1
        CALL MV7(VEC1,VEC2,LURF,LU1)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H(apr) times vector, part II '
          CALL WRTVCD(VEC1,LU1,1,LBLK)
        END IF
        HAPR00 = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
        CALL VECSMD(VEC1,VEC2,ONE,-HAPR00,LU1,LURF,LUVN,1,LBLK)
        XNORM = INPRDD(VEC1,VEC2,LUVN,LUVN,1,LBLK)
*
        IF(ABS(XNORM/HAPR00) .LE. 1.0D-12) THEN
         IHAPREIG = 1
        ELSE
         IHAPREIG = 0
        END IF
*
        IF(IHAPREIG.EQ.0) THEN
*         HAPRM100 =  <0!(H(apr)-E0)**-1 |0>
          CALL DIA0TRM_GAS(2,LURF,LU1,VEC1,VEC2,-E0RF)
          HAPRM100  = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
          WRITE(6,*) ' HAPRM100', HAPRM100
*. Obtain (H0-E0) ** (-1) |0> ( diagonal approx )
          CALL DIA0TRM_GAS(2,LURF,LUHI0,VEC1,VEC2,-E0RF)
C         SHIFT = -(E0RF-ECORE_H)
C         SHIFT_DIA = -E0RF
C         E0 = E0RF
C         IAPR = 1
C         IPERTOP = 1
C         IPROJ = 0
C         WRITE(6,*) ' SHIFT before call to HINTV ', SHIFT
C         CALL COPVCD(LURF,LU2,VEC1,1,LBLK)
C         CALL HINTV(LU2,LUHI0,SHIFT,SHIFT_DIA,VEC1,VEC2,LBLK,0,0 )
        ENDIF


      WRITE(6,*) '  HAPR00,  XNORM, IHAPREIG, HAPRM100 : ',
     &              HAPR00,  XNORM, IHAPREIG, HAPRM100

*
*
*. V times initial vector  on  LUVN
* ==================================
*
*. H0+V !0(0)> on LU1 (ECORE_HEX missing)
        IAPR = 0
        IPERTOP = 0
        CALL REWINE(LURF,-1)
        CALL REWINE(LU1,-1)
        CALL HTV(VEC1,VEC2,LURF,LU1)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Output from HTV '
          CALL WRTVCD(VEC1,LU1,1,LBLK)
        END IF
*
* V|0> = (H - H0) !0(0)> = (H' + ECORE_HEX - E0) !0(0)> on LUVN
        E0RFM =  - (E0RF-ECORE_HEX)
        CALL VECSMD(VEC1,VEC2,ONE,E0RFM,LU1,LURF,LUVN,1,LBLK)
*
*.S(11),V(11),EN(1)
*=================
        S(1) = INPRDD(VEC1,VEC2,LURF,LURF,1,LBLK)
        V(1) = INPRDD(VEC1,VEC2,LUVN,LURF,1,LBLK)
        H0(1) = E0RF
        EN(1) = V(1)
        IF(NTEST.GE.1) WRITE(6,*)
     &  ' Energy correction , n and E(n) ', 1,EN(1)
*
*
* =======================================================================
*.2               Loop over orders of correction vectors
* =======================================================================
*
      DO 1000 IORD = 1, MAXORD
*
* On entrance :  correction vectors 1 - IORD-1 on LUN
*                V!0(IORD-1)>                  on LUVN
*
*
* !0(n)>
* ======
*
*  E(n)!0(0)> -V!0(n-1)> on LU1
        CALL VECSMD(VEC1,VEC2,EN(IORD),ONEM,LURF,LUVN,LU1,1,LBLK)
*.  sum( l = 1,  IORD -1) (E( n- l) !0(l)> on LU2
        IF(IORD.GT.1) THEN
          DO II = 1, IORD -1
            SCR(II) = EN(IORD-II)
          END DO
          CALL MVCSMD(LUN,SCR,LU2,LUVN,VEC1,VEC2,IORD-1,1,LBLK)
*. add on LUVN
          CALL VECSMD(VEC1,VEC2,ONE,ONE,LU2,LU1,LUVN,1,LBLK)
        ELSE
          CALL COPVCD(LU1,LUVN,VEC1,1,LBLK)
          CALL REWINE(LUN,-1)
        END IF
*. project !0> component out, SAVE on LU2
        OVLAP = INPRDD(VEC1,VEC2,LURF,LUVN,1,LBLK)
        IF(NTEST.GE.2)  write(6,*) ' ovlap1  ', OVLAP
        CALL  VECSMD(VEC1,VEC2,ONE,-OVLAP,LUVN,LURF,LU2,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' RHS of lin.eq ,order ', IORD
          CALL WRTVCD(VEC1,LU2,1,LBLK)
        END IF
*. Multiply with (H0-E0)-1, save result on LUN
        IF(IH0DIA.NE.0) THEN
*. Multiply with inverted Diagonal
          CALL REWINE(LU2,-1)
          CALL REWINE(LU1,-1)
          CALL DIA0TRM_GAS(2,LU2,LU1,VEC1,VEC2,-E0RF)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' new correction vector of order ', IORD
            CALL WRTVCD(VEC1,LU1,1,LBLK)
          END IF
          IF(IHAPREIG.EQ.0) THEN
*. Orthogonalize with (H(apr)-E0)**-1|0>
            OVLAP = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
            IF(NTEST.GE.2)  write(6,*) ' ovlap2  ', OVLAP
*. Set (H(apr)-E0) ** 1 |0>
            CALL REWINE(LURF,-1)
            CALL REWINE(LUVN,-1)
            CALL DIA0TRM_GAS(2,LURF,LUVN,VEC1,VEC2,-E0RF)
            FACTOR = - OVLAP/HAPRM100
            CALL  VECSMD(VEC1,VEC2,ONE,FACTOR,LU1,LUVN,LU2,1,LBLK)
*. Save on LU1
            CALL COPVCD(LU2,LU1,VEC1,1,LBLK)
          END IF
        ELSE
*. Solve set of linear equations
          SHIFT = -(E0RF-ECORE_H)
          SHIFT_DIA = -E0RF
          E0 = E0RF
          IAPR = 1
          IPERTOP = 1
          IPROJ = 1
          WRITE(6,*) ' SHIFT before call to HINTV ', SHIFT
          CALL HINTV(LU2,LU1,SHIFT,SHIFT_DIA,VEC1,VEC2,LBLK,LURF,LUHI0)
          IF(NTEST.GE.1000) THEN
            WRITE(6,*) ' new correction vector '
            CALL WRTVCD(VEC1,LU1,1,LBLK)
          END IF
        END If
        X0N = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
        WRITE(6,*) ' Overlap <0!N> ', X0N
*Save on LUN
        CALL REWINE(LU1,-1)
        CALL COPVCD(LU1,LUN,VEC1,0,LBLK)
*. V!0(n)> on LUVN = (H - H0 )!0(n)> = (H'(holeform) - (H0 -ECORE_HEX))!0(n)>
* ================
*
*. H0+V !0(n)> on LU2 ( except ECORE_HEX ) |0(n)>
        CALL REWINE(LU1,-1)
        CALL REWINE(LU2,-1)
        CALL HTV(VEC1,VEC2,LU1,LU2)
*
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H times correction vector '
          CALL WRTVCD(VEC1,LU2,1,LBLK)
        END IF
*. Test : E(N+1) = <0! V !n> = <0! H !n>
        ENN = INPRDD(VEC1,VEC2,LURF,LU2,1,LBLK)
        IF(NTEST.GE.1) WRITE(6,*) ' TEST : ENN = ', ENN
*
*
* H0 |0(n)> = Q H apr |0(n)> on LUVN ( and include -ECORE_HEX missi
*
        E0 = E0RF
        SHIFT = ECORE_H-ECORE_HEX
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Input to H0TVM '
          CALL WRTVCD(VEC1,LU1,1,LBLK)
        END IF
        IPERTOP = 1
        IPROJ = 1
        CALL H0TVM(VEC1,VEC2,LU1,LUVN)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H0 times correction vector '
          CALL WRTVCD(VEC1,LUVN,1,LBLK)
        END IF
*. Project |0> out
        OVLAP = INPRDD(VEC1,VEC2,LURF,LUVN,1,LBLK)
        IF(NTEST.GE.2)  write(6,*) ' ovlap3  ', OVLAP
        OVLAPM = -OVLAP
        CALL  VECSMD(VEC1,VEC2,ONE,OVLAPM,LUVN,LURF,LU1,1,LBLK)
        CALL COPVCD(LU1,LUVN,VEC1,1,LBLK)
        OVLAP = INPRDD(VEC1,VEC2,LU1,LURF,1,LBLK)
*
* (H - H0) !0(n)> on LUVN
C       CALL VECSMD(VEC1,VEC2,ONE,ONEM,LU2,LUVN,LU1,1,LBLK)
        CALL VECSMD(VEC1,VEC2,ONE,ONEM,LU2,LU1,LUVN,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' V ! 0(n)> '
          CALL WRTVCD(VEC1,LUVN,1,LBLK)
        END IF
*. E(n+1) = <0(n)!V!0>
* ====================
        EN(IORD+1) = INPRDD(VEC1,VEC2,LUVN,LURF,1,LBLK)
        ENERGY = ENERGY + EN(IORD)
        IF(NTEST.GE.1) WRITE(6,*)
     &  ' Energy correction , n and E(n) ', IORD+1,EN(IORD+1)
*
*. Augment matrices H0, S and V
* ============================
*
*
*  H0(ij) = <0(i-1)!H0!0(j-1)>
*
      CALL REWINE(LUN,-1)
      DO JORD = 0, IORD
        IJ = (IORD+1)*(IORD+1-1)/2 + JORD+1
        IF(JORD.NE.0) THEN
          CALL REWINE(LU1,-1)
          H0(IJ) = INPRDD(VEC1,VEC2,LUN ,LU1,0,LBLK)
        ELSE
          H0(IJ) = INPRDD(VEC1,VEC2,LURF,LU1,1,LBLK)
        END IF
      END DO
*
*  s(ij) = <0(i-1)!0(j-1)>
*
*. Place correction vector !0(n)> on LU1
      CALL SKPVCD(LUN,IORD-1,VEC1,1,LBLK)
      CALL REWINE(LU1,-1)
      CALL COPVCD(LUN,LU1,VEC1,0,LBLK)
*
      CALL REWINE(LUN,-1)
      DO JORD = 0, IORD
        IJ = (IORD+1)*(IORD+1-1)/2 + JORD+1
        IF(JORD.NE.0) THEN
          CALL REWINE(LU1,-1)
          S(IJ) = INPRDD(VEC1,VEC2,LU1,LUN,0,LBLK)
        ELSE
          S(IJ) = 0.0D0
        END IF
      END DO
*
*  v(ij) = <0(i-1)!v!0(j-1)>
*
      CALL REWINE(LUN,-1)
      DO JORD = 0, IORD
        IJ = (IORD+1)*(IORD+1-1)/2 + JORD+1
        IF(JORD.NE.0) THEN
          CALL REWINE(LUVN,-1)
          V(IJ) = INPRDD(VEC1,VEC2,LUN,LUVN,0,LBLK)
        ELSE
          V(IJ) = INPRDD(VEC1,VEC2,LURF,LUVN,1,LBLK)
        END IF
      END DO
C!    IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Updated S matrix '
        CALL PRSYM(S,IORD+1)
        WRITE(6,*) ' Updated V matrix '
        CALL PRSYM(V,IORD+1)
        WRITE(6,*) ' Updated H0 matrix '
        CALL PRSYM(H0,IORD+1)
C!    END IF
*
*. Obtain additional energy expressions by 2n+1 rule.
* ===================================================
*
        DO N = IORD+1,2*IORD+1
*. E(N) = <K!V!N-K-1> - Sum(j=0,k-1)sum(m=k-j,n-1-j)<j+1!n-1-m-j>E(m)
*. Use K = IORD
           X = V((IORD+1)*(IORD+1-1)/2+N-IORD-1+1)
           DO J = 0,IORD-1
             DO M = IORD-J,N-1-J
               II = MAX(J+1,N-1-M-J)
               JJ = MIN(J+1,N-1-M-J)
               IJ = (II+1)*(II+1-1)/2 + JJ + 1
               X = X - S(IJ)*EN(M)
             END DO
           END DO
           EN(N) = X
        END DO
*
        IF(NTEST.GE.2 ) THEN
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*) ' =========================================== '
          WRITE(6,*) ' Energy corrections obtained in iteration ', IORD
          WRITE(6,*) ' =========================================== '
          WRITE(6,*)
          WRITE(6,*)
     &    '   Order       Energy correction      Total Energy '
          WRITE(6,*)
     &   ' ========================================================='
          ENERGY = E0RF+ECORE-ECORE_HEX
          DO JORD = 1, 2*IORD+1
            ENERGY = ENERGY + EN(JORD)
            WRITE(6,'(4X,I2,8X,1E18.10,6X,1E18.10)')
     &      JORD,EN(JORD),ENERGY
          END DO
        END IF
*. I can't wait to see the output, so XFLUSH
        LUOUT = 6
C       CALL  XFLUSH(LUOUT)
*
 1000 CONTINUE
*
      WRITE(6,*)
      WRITE(6,*) ' Zero order energy : ', E0RF+ECORE-ECORE_HEX
      WRITE(6,*)
      WRITE(6,*) ' =========================================== '
      WRITE(6,*) ' Energy corrections obtained as <0!V!0(n-1)> '
      WRITE(6,*) ' =========================================== '
      WRITE(6,*)
      WRITE(6,*)
     &'   Order         Energy correction          Total Energy '
      WRITE(6,*)
     &' ========================================================='
      ENERGY = E0RF+ECORE-ECORE_HEX
      DO IORD = 1, MAXORD
        ENERGY = ENERGY + EN(IORD)
        WRITE(6,'(4X,I2,8X,1E20.12,6X,1E22.14)')
     &  IORD,EN(IORD),ENERGY
      END DO
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) ' =========================================== '
      WRITE(6,*) ' Energy corrections obtained from 2n+1 rule  '
      WRITE(6,*) ' =========================================== '
      WRITE(6,*)
*
      WRITE(6,*)
     &  '   Order         Energy correction        Total Energy '
      WRITE(6,*)
     & ' ========================================================='
      DO JORD =MAXORD+1, 2*MAXORD+1
        ENERGY = ENERGY + EN(JORD)
        WRITE(6,'(4X,I2,8X,1E20.12,6X,1E22.14)')
     &  JORD,EN(JORD),ENERGY
      END DO
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Final S matrix '
        CALL PRSYM(S,MAXORD+1)
        WRITE(6,*) ' Final V matrix '
        CALL PRSYM(V,MAXORD+1)
        WRITE(6,*) ' Final H0 matrix '
        CALL PRSYM(H0,MAXORD+1)
      END IF

*
      RETURN
      END
