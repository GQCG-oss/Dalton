C  /* Deck cc_lanczos_newdrv */
C--------------------------------------------------------------------
      SUBROUTINE CC_LANCZOS_NEWDRV(WORK,LWORK,APROXR12)
C--------------------------------------------------------------------
C     Updated version of asymmetric lanczos drive
C     Purpose: build the Lanczos recursive chain and diagonalize
C     Handling of RHS vectors is for the time being done
C     like in CC_XETST
C     FIXME: RESTRUCTURE THE DUMP-ON-FILE FOR CCSDR3
C     SYMMETRY inserted in noddy way.
C
C     Sonia Coriani & Kristian Sneskov, August 2010
C     Sonia Coriani & Janusz Cukras, August 2013
C     Sonia Coriani, restart option introduced January 2015
C     Sonia Coriani, alternative start seeds
C--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
#include "cclrlancz.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "maxorb.h"
#include "ccroper.h"
#include "ccr1rsp.h"
#include "cco1rsp.h"
#include "ccx1rsp.h"
#include "ccqlrlcz.h"
#include "ccfro.h"
#include "cclists.h"
#include "dummy.h"
#include "codata.h"
#include "ccexcinf.h"
#include "ccexci.h"
!try intro sym...
#include "cclrinf.h"
!solvent info
#include "ccslvinf.h"
#include "qm3.h"
!mlcc multi_level
#include "ccpert.h"
!cvs: core-valence separation
#include "cclanczcvs.h"

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_LANCZOS_NEWDRV> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
      INTEGER MXXETRAN, MXVEC, MXFQTRAN
      PARAMETER (MXXETRAN = 20, MXVEC = 10, MXFQTRAN=50)
      integer ifile, kRfull, kLfull, KoffR,KoffL,istart,NrExc,koffre
      INTEGER LWORK, LENGTH,Jkaede, kPTcsi,ketaQ,kinverse
      INTEGER KFTRANS, kftransIm
      INTEGER kLe1,lensym, length_check,msym,NTRIP,IEX,krod,lentot
      INTEGER JOLDTST,iexidx
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FREQ, ECURR, rxs, test
      REAL DDOT
      REAL ZERO, ONE, TWO, FOUR, FIVE, Q1norm, factor
      REAL ovlpLR, xnorm1, ovlpRR, ovlpLL
      real sum11,csinorm, etanorm
#else 
      !sonia: enlarge alpha by 1 to check residuals
      DOUBLE PRECISION WORK(LWORK)
      !DOUBLE PRECISION ALPHA(JCHAIN+1), BETA(jChain), xgamma(jchain)
      DOUBLE PRECISION rxs, test
      DOUBLE PRECISION FREQ, ECURR, PT1, xresid, Tif, Tfi
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, FIVE, Q1norm, factor
      DOUBLE PRECISION thrshold
      DOUBLE PRECISION ovlpLR, xnorm1
      DOUBLE PRECISION ovlpRR, ovlpLL, sum11, csinorm,etanorm, excien
      
#endif 
      PARAMETER (ZERO = 0.0D0, ONE  = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0, FIVE = 5.0D0)
      PARAMETER (thrshold = 1.0d-12)

      CHARACTER*(3) FILXI, FILETA, LISTL, APROXR12
      CHARACTER*(8) LABEL
      CHARACTER*(10) MODEL, CROUT
      LOGICAL LORX, LTWO
      INTEGER IOPTRES, N2VEC, NSIMTR, ISIDE, NSIM, IST, IERR
      INTEGER IXETRAN(MXDIM_XEVEC,MXXETRAN), NXETRAN, IDXR1HF
      !INTEGER IFQTRAN(MXDIM_FTRAN,JCHAIN), NFQTRAN, ICHAIN
      INTEGER ISYAMP,ISYM,ksym, ICHAIN
      INTEGER ISYM0, ISYHOP, IOPT, IOPER, ITRAN, IVEC, IDUM, ISYETA
      INTEGER IOPTRW,III
      INTEGER KEND0, ITEST, KEND3, LWRK3, IORDER, ILSTETA, ISIGN
      INTEGER KEND1, KDUM, KEND2, LWRK1, LWRK2
      INTEGER KTTRI, KWI, KEIVAL, KEIVER, KEIVEL
      INTEGER IDLSTL,ISYCTR,IDXR1,INUM,IREAL
      integer KEND6, LWRK6, KEND61, LWRK61, KFQQ, KFQQR
      INTEGER KETA1, KETA2, KRHO1, KRHO2, kp1,kq1
      INTEGER IRELAX, KCTR1, KCTR2
      INTEGER KETA12, kpnewa1, ks1, kr1, kpnewa2, kqnew1, kqnew2
      INTEGER KPNEW1, KPNEW2, KAQNEW1, KAQNEW2, kqold1, kqold2, kpold1
      INTEGER kpold2, KEND4, KEND5, LWRK4, LWRK5, ICPLX,IDXVEC,IOFF,JOFF
      integer iopera,ioperb,isyma,isymb,ioptseed
      integer kalpha,kbeta,kgamma
      CHARACTER*(8) LABELA,LABELB
      logical lorxa,lorxb,lcomplex
*
*     Files on which we dump the q and p(T) matrices
*
      integer LUQMAT, LUPMAT, LUTTRMAT
      CHARACTER FQMAT*7
      PARAMETER (FQMAT ='CCLAN_Q' )
      CHARACTER FPMAT*8
      PARAMETER (FPMAT ='CCLAN_PT')
      !CHARACTER  FTTRMAT*9
      !PARAMETER (FTTRMAT ='CCLAN_TTR') 
      CHARACTER*7 FNTmat
      PARAMETER (FNTmat='CCLAN_T') 
      INTEGER LUTmat


* external functions (lists):
      INTEGER ILSTSYM
      INTEGER IROPER
      INTEGER IR1TAMP
      INTEGER IRHSR1
      INTEGER IETA1
      INTEGER IQLLIST

      CALL QENTER('CC_LANCZOS_NEWDRV')
      CALL AROUND('CC_LANCZOS_NEWDRV (version 2)')
      write(lupri,*) 'CORE-VALENCE SEPARATION REQUESTED? ', LCVSLCZ
      write(lupri,*) 'IONISATION REQUESTED? ', LIONIZLCZ

      !CALL IZERO(NCCEXCI,3*8)
      IERR = 0

      !Open file to save the Triadiagonal Matrix
      !LUTTRMAT = -1
      !CALL WOPEN2(LUTTRMAT,FTTRMAT,64,0)

*---------------------------------------------------------------------*
* set up information for rhs vectors we want
*---------------------------------------------------------------------*
*  number of simultaneous transformations:

      NXETRAN = 1
*.......................................................................
*  compute the xi and eta vector for a non-relaxed one-electron 
*          perturbation (does only require that the integrals for 
*          the operator are available on file)
* This info will have to be passed from input routine!!!!!!!!!!!!
* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
      ITEST = 2
      IORDER = 1

      !LABEL  = LABELO
      !LORX   = .FALSE.
      !LTWO   = .FALSE.
      !FREQ   = ZERO
      !ISYHOP = 1  !change according labelo symmetry

      !this loop should be outside everyhing but we are hacking
      !NLROP SHOULD BE 1

      DO IOPER = 1, NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        LORXA  = LALORX(IOPER)
        LORXB  = LBLORX(IOPER)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        !LPDBSA = LPDBSOP(IOPERA)
        !LPDBSB = LPDBSOP(IOPERB)
       END DO

       FREQ   = ZERO
*---------------------------------------------------------------------*
      ! allow extensions in the vector lists for the next few lines
!      if (isyma.eq.isymb) then
!
!
!           CALL CC_SETXE('Eta',IXETRAN,IEDOTS,MXTRAN,MXVEC,
!     &                   0,IOPERA,IKAPA,0,0,0,IR1VECB,ITRAN,IVEC)
!
!           CALL CC_SETXE('Eta',IXETRAN,IEDOTS,MXTRAN,MXVEC,
!     &                      0,IOPERB,IKAPB,0,0,0,IR1VECA,ITRAN,IVEC)
!
!           CALL CC_SETXE('Xi ',IXETRAN,IXDOTS,MXTRAN,MXVEC,
!     &                      0,IOPERA,IKAPA,0,0,0,IL1VECB,ITRAN,IVEC)
!
!           CALL CC_SETXE('Xi ',IXETRAN,IXDOTS,MXTRAN,MXVEC,
!     &                      0,IOPERB,IKAPB,0,0,0,IL1VECA,ITRAN,IVEC)
!

      ! allow extensions in the vector lists for the next few lines
      if (isyma.eq.isymb) then
       isyhop  = isymb
       label   = labela
       ltwo    = .false.
       lorx    = .false.
       write(lupri,*)'isyhop',isyhop,' label=', label
       LOPROPN = .TRUE.
       LO1OPN  = .TRUE.
       LX1OPN  = .TRUE.
       !LQLOPN  = .TRUE.
 
       IRELAX = 0
       IF (LORXA) THEN
         IRELAX = IR1TAMP(LABELA,LORXA,FREQ,ISYMA)
       END IF
       IF (LORXB) THEN
         IRELAX = IR1TAMP(LABELB,LORXB,FREQ,ISYMB)
       END IF

       IF (LORX) THEN
         IRELAX = IR1TAMP(LABEL,LORX,FREQ,ISYHOP)
       END IF
       LPDBSOP(IROPER(LABEL,ISYHOP)) = LTWO

       ! set the IXETRAN array for one (XI,ETA) pair
       IXETRAN(1,1) = IROPER(LABELA,ISYHOP)
       IXETRAN(2,1) = 0
       IXETRAN(3,1) = IRHSR1(LABELA,LORX,FREQ,ISYHOP)
       IXETRAN(4,1) = IETA1(LABELA,LORX,FREQ,ISYHOP)
       IXETRAN(5,1) = IRELAX

       IXETRAN(1,2) = IROPER(LABELB,ISYHOP)
       IXETRAN(2,2) = 0
       IXETRAN(3,2) = IRHSR1(LABELB,LORX,FREQ,ISYHOP)
       IXETRAN(4,2) = IETA1(LABELB,LORX,FREQ,ISYHOP)
       IXETRAN(5,2) = IRELAX
C
C Set dimensions of full space and read/write options for
C vectors on files
C
      IF (CCS) THEN
         LENGTH = NT1AM(ISYHOP)
         IOPTRW = 1
         N2VEC  = 0
         LENTOT = NT1AMX
      ELSE 
         LENTOT = NT1AMX+NT2AMX
         IOPTRW = 3
         N2VEC  = 1
         LENGTH = NT1AM(ISYHOP) + NT2AM(ISYHOP)
      END IF
C
C Check dimension of Lanczos space
C
      if (lchainrestart) then
        write(lupri,*)'This is a restart from a previous Lanczos calc'
        write(lupri,*)'Previous LCHAIN =', JCHAINOLD
        write(lupri,*)'New LCHAIN length (total)=', JCHAIN
      end if
      WRITE(LUPRI,*) 'LANCZOS: Chain length by input, JCHAIN=', JCHAIN

      Jkaede = Jchain

      IF (JCHAIN.GT.LENGTH) THEN
         WRITE(LUPRI,*) 'WARNING: JCHAIN EXCEED FULL SPACE DIMENSION'
         Jkaede = length
         WRITE(LUPRI,*) 'WE RESET IT TO dim{AMAT} =   ', jkaede
      END IF

      ! disallow again extension in the vector lists
      LOPROPN = .FALSE.
      LO1OPN  = .FALSE.
      LX1OPN  = .FALSE.

      if (locdbg) then
        CALL AROUND('CC_LANCZOS: test of call to CC_XIETA setup')
        WRITE (LUPRI,*) 'ITEST  =',ITEST
        WRITE (LUPRI,*) 'ISYHOP =',ISYHOP
        WRITE (LUPRI,*) 'LABELA  =',LABELA
        WRITE (LUPRI,*) 'LABELB  =',LABELB
        WRITE (LUPRI,*) 'LTWO   =',LTWO
        WRITE (LUPRI,*) 'LORX   =',LORX
        WRITE (LUPRI,*) 'FREQ   =',FREQ
        WRITE (LUPRI,*) 'IROPER =',IXETRAN(1,1),' and ', IXETRAN(1,2)
        WRITE (LUPRI,*) 'ILEFT  =',IXETRAN(2,1),' and ', IXETRAN(2,2)
        WRITE (LUPRI,*) 'IRHSR1 =',IXETRAN(3,1),' and ', IXETRAN(3,2)
        WRITE (LUPRI,*) 'IETA1  =',IXETRAN(4,1),' and ', IXETRAN(4,2)
        WRITE (LUPRI,*) 'IRELAX =',IXETRAN(5,1),' and ', IXETRAN(5,2)
        WRITE (LUPRI,*) 'JCHAIN =', JCHAIN
        WRITE (LUPRI,*) 'JKAEDE =', JKAEDE

      end if

*---------------------------------------------------------------------*
* call CC_XIETA to calculate the Xi and Eta vectors:
* they are dumped directly on file within XIETA
*---------------------------------------------------------------------*

      KEND0 = 1
      KEND1 = KEND0 + 1
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND1, 'Available:', LWORK
         CALL QUIT('Insuff work space in CC_LANCZOS_DRV (1)')
      END IF

      LISTL  = 'L0 '
      FILXI  = 'O1 '
      FILETA = 'X1 '
      IOPTRES = 3
      CALL AROUND('CC_LANCZOS_DRV: call CC_XIETA (both Xi and Eta)')
      CALL CC_XIETA(IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL, 
     &              FILXI,  idummy, dummy,
     &              FILETA, idummy, dummy,
     &              .FALSE.,0,  WORK(KEND1), LWRK1 )
!
!     open files for Q and P matrices
!
      LUQMAT = -1
      LUPMAT = -1
      CALL WOPEN2(LUQMAT,FQMAT,64,0)
      CALL WOPEN2(LUPMAT,FPMAT,64,0)

*---------------------------------------------------------------------*
*     Allocate for RHS vector to be read in after having been
*     computed by the XIETA routine
*---------------------------------------------------------------------*

      KQnew1 = 1
      if (ccs) then
        Kend1  = KQnew1 + NT1AM(ISYHOP)
        KQnew2 = Kend1
      else
        KQnew2 = KQnew1 + NT1AM(ISYHOP)
        KEND1  = KQnew2 + NT2AM(ISYHOP)
      end if
      LWRK1  = LWORK   - KEND1
!
      KAQnew1 = KEND1
      if (ccs) then
        kend2   = KAQnew1 + NT1AM(ISYHOP)
        KAQnew2 = kend2
      else
        KAQnew2 = KAQnew1 + NT1AM(ISYHOP)
        KEND2   = KAQnew2 + NT2AM(ISYHOP)
      end if
      LWRK2   = LWORK - KEND2
!
      IF (LWRK1.LT.0 .OR. LWRK2.LT.0) THEN
         WRITE(LUPRI,*) 'Need:', KEND1,',',KEND2, 'Available:', LWORK
         CALL QUIT('Insuff memory in CC_LANCZOS_NEWDRV')
      END IF
C     -------------------------------------------------
C     read in the Csi^Z vector from CC_XIETA and print
C     we do this also in the restart case to recover the
C     norm of the vector for the LRF calculation
C     -------------------------------------------------
      NSIMTR = NXETRAN
      ECURR  = ZERO


      !To Janusz: NOW Need to read 2 vectors, ITRAN=1,2...
      if (Sn_seeds) then
         ITRAN = 1
         IVEC = IXETRAN(3,ITRAN)
         CALL CC_RDRSP('O1 ',IVEC,ISYHOP,IOPTRW,MODEL,
     &                WORK(KQnew1),WORK(KQnew2))
         !
         !I read O1 as dirty trick to get MODEL
         !
         write(lupri,*)'LANCZOS: USE SN Q1 SEEDS'
         call flshfo(lupri)
         ioptseed = 1   !right seed
         CALL DZERO(WORK(KQNew1),length)
         call cc_lanczos_Sn_seeds(ioptseed,work(KQnew1),
     &        labela,isyhop,excien,work(kend2),lwrk2)
      q1norm = dsqrt(abs(ddot(length,WORK(KQnew1),1,WORK(KQnew1),1)))
      write(lupri,*)'q1norm AxRe: ', q1norm, ' exci en:', excien
      else
         write(lupri,*)'LANCZOS: USE CSIA SEEDS'
         !ITRAN = 2
         ITRAN = 1
         IVEC = IXETRAN(3,ITRAN)
         CALL CC_RDRSP('O1 ',IVEC,ISYHOP,IOPTRW,MODEL,
     &                WORK(KQnew1),WORK(KQnew2))
      q1norm = dsqrt(abs(ddot(length,WORK(KQnew1),1,WORK(KQnew1),1)))
      write(lupri,*)'q1norm CSIX: ', q1norm
      end if

      if (locdbg) then
        CALL AROUND('The Q mat seed vector Q1:')
        !CALL CC_PRP(WORK(KQnew1),WORK(KQnew2),ISYHOP,1,N2VEC)
      end if

      q1norm = dsqrt(abs(ddot(length,WORK(KQnew1),1,WORK(KQnew1),1)))
      write(lupri,*)'q1norm csinorm: ', q1norm
      !write(lupri,*)'SONIA: TRDRV Analyze before prj'
      !CALL CC_PRAM0(WORK(KQNew1),PT1,ISYHOP,.FALSE.)
      !write(lupri,*)'SONIA: Alternative Analyze before prj'
      !CALL CCSD_TCMEPK_skod(WORK(KQnew2),1.0d0,ISYHOP,1)
!
      IF (LCVSLCZ) THEN
         write(lupri,*)'SONIA: TRDRV Project out ai, i=/icore'
         !CALL CC_FREEZE_VAL(WORK(KQnew1),ISYHOP,.false.)
         !CALL CCSD_FREEZE_EXCI2(WORK(KQnew1),ISYHOP,
         !write(lupri,*) MAXCORE,MAXION
         !write(lupri,*) NRHFCORE
         !write(lupri,*) IRHFCORE
         !write(lupri,*) NVIRION
         !write(lupri,*) IVIRION
         CALL CC_FREEZE_EXCI(WORK(KQnew1),ISYHOP,
     &        MAXCORE,MAXION,
     &        NRHFCORE,IRHFCORE,
     &        NVIRION,IVIRION,LBOTHLCZ)
      END IF

      csinorm = q1norm
      factor = ONE/q1norm
      CALL DSCAL(length,factor,WORK(KQnew1),1)
      !first q vector is normalized against itself
      if (locdbg) then
        CALL AROUND('NORMALIZED Q1 vector (= Xi^Z):')
        !CALL CC_PRP(WORK(KQnew1),WORK(KQnew2),ISYHOP,1,N2VEC)
      end if
      q1norm = sqrt(abs(ddot(length,WORK(KQnew1),1,WORK(KQnew1),1)))
      write(lupri,*)'q1norm (2): ', q1norm
!
      !dump q1 onto the QMAT file
      !if (locdbg) then
        !CALL AROUND('analytical Q1 vector:')
        !CALL CC_PRP(WORK(KQnew1),WORK(KQnew2),ISYHOP,1,N2VEC)
      !end if

      if (.not.lchainrestart) then
         write(lupri,*)'Dumping Q1 vector on file (Qmat)'
         CALL CC_WVEC(LUQMAT,FQMAT,length,length,1,WORK(KQnew1))
      end if
!
!      Read in Xi^Z as starting guess for p1 vector
!
      ISYCTR = ILSTSYM(LISTL,IDLSTL)
      ISYETA = MULD2H(ISYCTR,ISYHOP)

      KALPHA  = KEND2
      KBETA   = KALPHA + JCHAIN + 1
      KGAMMA  = KBETA  + JCHAIN
      KPnew1  = KGAMMA + JCHAIN
      !KPnew1  = KEND2
      KPnew2  = KPnew1   + NT1AM(ISYETA)
      KPnewA1 = KPnew2   + NT2AM(ISYETA)
      KPnewA2 = KPnewA1  + NT1AM(ISYETA)
      KEND3   = KPnewA2  + NT2AM(ISYETA)
      LWRK3   = LWORK   - KEND3
!
      IF (LWRK3.LT.0) THEN
         WRITE(LUPRI,*) 'Need:', KEND3, 'Available:', LWORK
         CALL QUIT('Insufficient memory in CC_LANCZOS.')
      END IF
C     -------------------------------------------------
C     read the eta vector from CC_XIETA 
C     we do this also in the restart case to recover the
C     norm of the vector for the LRF calculation
C     -------------------------------------------------
      if (Sn_seeds) then
!         ITRAN = 1
!         IVEC = IXETRAN(4,ITRAN) 
!         WRITE (LUPRI,*) 'CC_LANCZOS reading ETA^X : IVEC = ',IVEC
!         CALL CC_RDRSP('X1 ',IVEC,ISYETA,IOPTRW,MODEL,
!     &                WORK(KPnew1),WORK(KPnew2))

         write(lupri,*)'LANCZOS: USE SN P1 SEEDS'
         call flshfo(lupri)
         ioptseed = 2   !left seed
         CALL DZERO(WORK(KPNew1),length)
         call cc_lanczos_Sn_seeds(ioptseed,work(KPnew1),
     &        labela,isyhop,excien,work(kend3),lwrk3)
         q1norm = ddot(length,WORK(KPnew1),1,WORK(KPnew1),1)
         write(lupri,*)'LE*Ax norm from ETAC: ', q1norm
      else
         write(lupri,*)'LANCZOS: USE ETA^X SEEDS'
         ITRAN = 1
         IVEC = IXETRAN(4,ITRAN) 
         WRITE (LUPRI,*) 'CC_LANCZOS reading ETA^X : IVEC = ',IVEC
         CALL CC_RDRSP('X1 ',IVEC,ISYETA,IOPTRW,MODEL,
     &                WORK(KPnew1),WORK(KPnew2))
         q1norm = ddot(length,WORK(KPnew1),1,WORK(KPnew1),1)
         write(lupri,*)'etaC  norm from XIETA: ', q1norm

      end if

      ! q1norm = ddot(length,WORK(KPnew1),1,WORK(KQnew1),1)
      ! write(lupri,*)'p1*q1 norm from XIETA: ', q1norm

      !save a copy for later use (resp amplitudes)
      !call dcopy(length,work(KPnew1),1,work(keta),1)
      !binormalize vs q1
!----
      IF (LCVSLCZ) THEN
         write(lupri,*)'SONIA: P1 Project out ai, i=/icore'
         !CALL CC_FREEZE_VAL(WORK(KPnew1),ISYHOP,.false.)
         !CALL CCSD_FREEZE_EXCI2(WORK(KPnew1),ISYHOP,
         CALL CC_FREEZE_EXCI(WORK(KPnew1),ISYHOP,
     &        MAXCORE,MAXION,
     &        NRHFCORE,IRHFCORE,
     &        NVIRION,IVIRION,LBOTHLCZ)
      END IF

      q1norm = ddot(length,WORK(KPnew1),1,WORK(KQnew1),1)
      write(lupri,*)'p1*q1 norm: ', q1norm
      etanorm = q1norm
      factor = ONE/q1norm
      CALL DSCAL(length,factor,WORK(KPnew1),1)
      !normalize vs q vector, not against itself
      if (locdbg) then
         CALL AROUND('BINORMALIZED P1 vector (= eta^Z):')
         !CALL CC_PRP(WORK(KPnew1),WORK(KPnew2),ISYETA,1,N2VEC)
      end if
      q1norm = ddot(length,WORK(KPnew1),1,WORK(KQnew1),1)
      write(lupri,*)'p1*q1 norm (2): ', q1norm
      !save a copy for later use (resp amplitudes)
      !call dcopy(length,work(KPnew1),1,work(keta),1)
!----
      !dump p1 onto the PMAT file
      if (lchainrestart) then

         !apro file formattato per salvare i coefficienti
         write(lupri,*)'LANCZOS RESTART OPEN LU_TMAT'
         call flshfo(lupri)
         LUTmat=-1
         CALL GPOPEN(LUTmat,FNTmat,'OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

         !use Aq for q0 and q for q1
         !use pA for p0 and p for p1
         CALL DZERO(WORK(Kalpha),Jchainold)
         CALL DZERO(WORK(Kbeta), Jchainold)
         CALL DZERO(WORK(Kgamma),Jchainold)
         write(lupri,*)'Loading restart T matrix'
         !CALL CC_RVEC(LUTTRMAT,FTTRMAT,1,1,1,JOLDTST)
         rewind LuTmat
         read(LuTmat,*) JOLDTST

         if (joldtst==jchainold) then
             write(lupri,*) 'T^old matrix long', joldtst
             write(lupri,*) 'and I restart from',JCHAINOLD
         else
             write(lupri,*) 'T^old matrix long:',joldtst
             write(lupri,*) 'I will only read:',JCHAINOLD
         end if
         read(LuTmat,*) (WORK(KALPHA+j-1),j=1,JCHAINOLD) 
         read(LuTmat,*) (WORK(KBETA+j-1),j=1,JCHAINOLD) 
         read(LuTmat,*) (WORK(Kgamma+j-1),j=1,JCHAINOLD) 
!
!SUBROUTINE CC_WVEC(LU,FIL,LLEN,LEN,NR,VEC)
!IOFF = 1 + LLEN*(NR-1)
!
!         CALL CC_RVEC3(LUTTRMAT,FTTRMAT,1,
!     &                JCHAINOLD,1,1,WORK(KALPHA))
!         CALL CC_RVEC3(LUTTRMAT,FTTRMAT,Joldtst,
!     &                JCHAINOLD,2,1,WORK(KBETA))
!         CALL CC_RVEC3(LUTTRMAT,FTTRMAT,Joldtst,
!     &                JCHAINOLD,3,1,WORK(KGAMMA))

         !this step is not strictly necessary
         !we recompute the norms 
!         write(lupri,*)'Csinorm built:', csinorm
!         CALL CC_RVEC(LUTTRMAT,FTTRMAT,JCHAINOLD,
!     *               1,4,CSINORM)
!         write(lupri,*)'Csinorm read :', csinorm
!         write(lupri,*)'Etanorm built:', etanorm
!         CALL CC_RVEC(LUTTRMAT,FTTRMAT,JCHAINOLD,
!     *               1,5,ETANORM)
!         write(lupri,*)'Etanorm read :', etanorm

         if (locdbg) then
            DO i=1,jchainold
               write(lupri,*)'a:',WORK(KALPHA+i-1)
            END DO
            DO i=1,jchainold
               write(lupri,*)'b:',WORK(KBETA+i-1)
            END DO
            DO i=1,jchainold
               write(lupri,*)'g:',WORK(KGAMMA+i-1)
            END DO
         end if

         K = jchainold-1
         CALL CC_RVEC(LUQMAT,FQMAT,LENGTH,LENGTH,
     &                K,WORK(KAQNew1))
!       write(lupri,*)'q0 2',ddot(length,WORK(KAqnew1),1,WORK(KAqnew1),1)

         K = jchainold
         CALL CC_RVEC(LUQMAT,FQMAT,LENGTH,LENGTH,
     &                K,WORK(KQNEW1))
!       write(lupri,*)'q1',ddot(length,WORK(KQnew1),1,WORK(KQnew1),1)

         K = jchainold-1
         CALL CC_RVEC(LUPMAT,FPMAT,LENGTH,LENGTH,
     &                K,WORK(KPnewA1))
!       write(lupri,*)'p0',ddot(length,WORK(KPnewA1),1,WORK(KPnewA1),1)

         K = jchainold
         CALL CC_RVEC(LUPMAT,FPMAT,LENGTH,LENGTH,
     &                K,WORK(KPNEW1))
!       write(lupri,*)'p1',ddot(length,WORK(KPnew1),1,WORK(KPnew1),1)

         write(lupri,*)'Calling MAKECHAIN_rst to build Chain'
         call cc_lanczos_makechain_rst(WORK(KAqnew1),WORK(KPnewA1),
     &           WORK(KQnew1),WORK(KPnew1),Length,
     &           WORK(KALPHA),WORK(KBETA),WORK(KGAMMA),
     &           ISYHOP,WORK(KEND3),LWRK3,APROXR12,JCHAINOLD,JKAEDE,
     &           LUQMAT,FQMAT,LUPMAT,FPMAT,LUTMAT,FNTMAT)
!         
         write(lupri,*)'LANCZOS RESTART CLOSE LU_TMAT'
         call flshfo(lupri)
         call gpclose(luTmat,'KEEP')
      else
        write(lupri,*)'NO restart, regular makechain routine'
        write(lupri,*)'Dump P1 vector on file (PTmat)'
        IVEC  = 1
        CALL CC_WVEC(LUPMAT,FPMAT,length,length,ivec,WORK(KPnew1))

        write(lupri,*)'LANCZOS: p1*q1 norm', 
     &             ddot(length,WORK(KPNEW1),1,WORK(KQNEW1),1)
!
!   make chain
!            
        write(lupri,*)'Calling MAKECHAIN routine to build Chain'
        call cc_lanczos_makechain(WORK(KQnew1),WORK(KPnew1),Length,
     &           WORK(KALPHA),WORK(KBETA),WORK(KGAMMA),
     &           ISYHOP,WORK(KEND3),LWRK3,APROXR12,JKAEDE,
     &           LUQMAT,FQMAT,LUPMAT,FPMAT)

       end if


!-----------------------------------------
! ERROR ESTIMATE....
!-----------------------------------------
!     compute q_jchain+1
!----------------------------------------
!BUILD AND DIAGONALIZE TRIDIAGONAL MATRIX
!----------------------------------------

      !I believe I can re-initialize memory here!!!!
      !KTTRI = 1
      !take care, now we save q_jchain+1
!
!     Allocate for the eigenvalues and eigenvectors of the T matrix
!
      KEIVAL = KEND3                   !real eigenvalues
      KWI    = KEIVAL + JKAEDE         !imaginary eigenvalues
      KEIVER = KWI    + jkaede         !R eigenvectors
      KEIVEL = KEIVER + jkaede*jkaede  !L eigenvectors
      KEND6  = KEIVEL + jkaede*jkaede
      LWRK6  = LWORK - KEND6
      if (LWRK6.lt.0) then
        write(lupri,*)'Need:', kend6, ' Available:', lwork
        call quit('Not enough memory in Lanczos')
      end if
      CALL DZERO(WORK(KEIVAL),JKAEDE)
      CALL DZERO(WORK(KWI),JKAEDE)
      CALL DZERO(WORK(KEIVER),JKAEDE*JKAEDE)
      CALL DZERO(WORK(KEIVEL),JKAEDE*JKAEDE)

      lcomplex = .false.
!      call cc_lanczos_Tmatrix(work(kalpha),work(kbeta),work(kgamma),
!     &     Jkaede,work(keival),work(kwi),work(keiveR),work(keiveL),
!     &     csinorm,etanorm,luttrmat,fttrmat,work(kend6),lwrk6,lcomplex,
!     &     IPRLRLCZ)
      call cc_lanczos_Tmatrix(work(kalpha),work(kbeta),work(kgamma),
     &     Jkaede,work(keival),work(kwi),work(keiveR),work(keiveL),
     &     csinorm,etanorm,work(kend6),lwrk6,
     &     lcomplex,IPRLRLCZ)
!
!biorthonormalize lanczos R and L eigenvectors
!
      call CC_BIORTOCOMP(Jkaede,work(keival),work(kwi),
     &                   work(keiver),work(keivel),
     &                   WORK(KEND6),LWRK6)
!
! Test biorthonormality
!
!      ioff = 0
!      do i=1,jkaede
!         joff = 0
!         do j = 1, jkaede
!           ovlpLR = DDOT(jkaede,WORK(KEIVEL+IOFF),1,WORK(KEIVER+JOFF),1)
!           joff = joff+jkaede
!            write(lupri,*)'<L_',i,'|R_',j,'> = ', ovlpLR
!         end do
!         ioff = ioff+jkaede
!      end do
!---------------------------------------------------
!check residual according to Lanczos iteration check 
!---------------------------------------------------
        !fix this as I messed up things a bit...
        !call cc_resid_lanczos(beta(jchain),length,jkaede,

        call cc_resid_lanczos(work(kbeta+jkaede-1),length,jkaede,
     &       work(kqnew1),work(keiver),work(keival),work(kwi))
C
C------------------------------------------------------------------
C     Setup for build of QQ^F matrix 
C-------------------------------------------------------------------
C
!
!   make q vectors lists and prepare Fq lists
!
      !allow extension in the vector lists
      LQLOPN  = .true.
      do i=1,jkaede
         ICHAIN = i
         ivec = IQLLIST(LABELB,LORX,ICHAIN,FREQ,ISYHOP)
      end do
      ! disallow again extension in the vector lists
      LQLOPN  = .FALSE.

      if (.not.Sn_Seeds) then
      write(lupri,*)'Call cc_lanczos_qlist'
      call cc_lanczos_qlist(IQLLIST,LABELB,LUQMAT,FQMAT,
     &                      jkaede,isyhop,model,
     &                      WORK(KEND6),LWRK6)
      end if
!
      if (lcomplex) then
         kftransIm = kend6
         kftrans   = kftransIm + jkaede*jkaede
      else
         kftransIm = kend6
         kftrans   = kend6
      end if
      kfqq    = kftrans + jkaede*jkaede
      kfqqr   = kfqq    + jkaede*jkaede
      kend61  = kfqqr   + jkaede*jkaede
      lwrk61  = lwork - kend61
      IF (LWRK61 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need:', KEND6, 'Available:', LWORK
         CALL QUIT('Insufficient memory in CC_LANCZOS FQQ')
      END IF

       if (Sn_seeds) then
          WRITE(LUPRI,*) 'FMATRIX is skipped for the SEED option'
          call dzero(work(kfqq),jkaede*jkaede)
       else

       IF (OLDFMAT) THEN
          WRITE(LUPRI,*) 'Calling OLD LCZ_FMATRIX to get QQF matrix'
          CALL CC_LANCZOS_FMATRIX_OLD(JKAEDE,ISYHOP,IOPTRW,
     &                 WORK(kfqq),WORK(KEND61),LWRK61)

       ELSE
         WRITE(LUPRI,*) 'Calling NEW LCZ_FMATRIX to get QQF matrix'
         CALL CC_LANCZOS_FMATRIX(JKAEDE,ISYHOP,LABELB,IOPTRW,
     &                 WORK(kfqq),WORK(KEND61),LWRK61)
       END IF
       end if
!
       if (locdbg) then
        CALL AROUND('>>> The QQF matrix (in driver) <<<')
        CALL OUTPUT(WORK(KFQQ),1,JKAEDE,1,JKAEDE,JKAEDE,JKAEDE,1,LUPRI)
       end if
C
C------------------------------------------------------------------
C     Generate trans^F
C-------------------------------------------------------------------
C
       CALL DZERO(work(kFtrans),jkaede*jkaede)
       if (lcomplex) CALL DZERO(work(kFtransIm),jkaede*jkaede)
       if (Sn_seeds) then
          WRITE(LUPRI,*) 'Ftrans skipped for SEED'
       else
       WRITE(LUPRI,*) 'Calling Lanczos_Ftrans to make transF (calF)'
       call flshfo(lupri)
       CALL CC_Lanczos_FQQRR(IOPTRW,JKAEDE,WORK(KFQQ),
     &               WORK(KEIVER),
     &               WORK(KFtrans),WORK(KFtransIm),WORK(KFQQR),
     &               WORK(KEND61),LWRK61,
     &               WORK(KEIVAL),WORK(KWI),isyhop,lcomplex)
       end if
!       end if

      if (locdbg) then
         !CALL AROUND('>>> The trans^F matrix <<<')
         !CALL OUTPUT(Ftrans,1,jkaede,1,jkaede,jkaede,jkaede,1,LUPRI)
      end if
C
C-------------------------------------------------------------
C     Compute the response functions at required freqs and for
C     required gammas from input
C     Results are dumped on a separate file!!!!
C-------------------------------------------------------------
C
      if (sum_rules) then
        CALL CC_LANCZOS_SUMRULES(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
     &                  WORK(KEIVER),
     &                  WORK(KEIVEL),WORK(KFTRANS),ETANORM,CSINORM)
      else if (Sn_seeds) then
        write(lupri,*) ' Call CC_LANCZOS_T_if'
        CALL CC_LANCZOS_T_if(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
     &                 WORK(KEIVER),WORK(KEIVEL),ETANORM,CSINORM,
     &                 excien)

      else
        CALL CC_LANCZOS_LR(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
     &                  WORK(KEIVER),
     &                  WORK(KEIVEL),WORK(KFTRANS),ETANORM,CSINORM)
!        CALL CC_LANCZOS_LR(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
!     *                  WORK(KEIVER),WORK(KEIVEL),
!     *                  WORK(KFTRANS),WORK(KFTRANSIm),ETANORM,CSINORM)
        !compute oscillator strengths as residue of Lanczos LR
        CALL CC_LANCZOS_OSCSTR(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
     &                 WORK(KEIVER),
     &                 WORK(KEIVEL),WORK(KFTRANS),ETANORM,CSINORM)
      end if
C
C--------------------------------------------------------------
C ANALYZE VECTORS AND DUMP ON FILE IF REQUIRED
C--------------------------------------------------------------
C
      !For comparison with Standard code we need to renormalize

!      CALL WCLOSE2(LUQMAT,FQMAT,'KEEP')
!      CALL WCLOSE2(LUPMAT,FPMAT,'KEEP')
!      LUQMAT = -1
!      LUPMAT = -1
!      CALL WOPEN2(LUQMAT,FQMAT,64,0)
!      CALL WOPEN2(LUPMAT,FPMAT,64,0)
!
!
      !disable test of moments by expansion in full space
      if (.false.) then
         nrexc = 4  !HACK HACK
         do i=1,nrexc
            KRfull = kftrans
            kLfull = KRfull + length
            kend6  = KLfull + length
            kq1    = kend6 
            kp1    = kq1 + length
            kend6  = kp1 + length
            lwrk6  = lwork  - kend6
            IF (LWRK6 .LT. 0) THEN
              WRITE(LUPRI,*) 'Need:', KEND6, 'Available:', LWORK
              CALL QUIT('NEWDRV: Insuff space in full RE/LE build ')
            END IF
            CALL DZERO(WORK(KRfull),length)
            CALL DZERO(WORK(KLfull),length)
        WRITE(LUPRI,*) 'CC_build_Rfull3 2 build R vector nr', i
              !Lanczos eigevector is "jkaede" long, NOT "length"
              !IDENTIFY THE LANCZOS EIGENVECTORS TO BE TRANFORMED
              KoffR = keiver + (i-1)*jkaede
              KoffL = keivel + (i-1)*jkaede
!
! Change this to use the Q files instead of FQMAT...
!
              call CC_build_Rfull3(LUQMAT,FQMAT,Jkaede,length,
     &                    WORK(KRfull),work(koffR),
     &                    work(keival+i-1),
     &                    work(kwi+i-1),
     &                    WORK(KEND6),LWRK6)
!
              call CC_build_Rfull3(LUPMAT,FPMAT,Jkaede,length,
     &                    WORK(KLfull),work(koffL),
     &                    work(keival+i-1),
     &                    work(kwi+i-1),
     &                    WORK(KEND6),LWRK6)
!!
!!            normalize chosen RE vector 
!!
              iopt=1
              WRITE(LUPRI,*) 'cc_normalize selected (full) RE vector'
              call cc_normalize(length,1,work(kRfull),work(krfull),
     &                    WORK(KEIVAL+i-1),
     &                    WORK(KWI+i-1),
     &                    isyhop,iopt)
!!
!!            binormalize chosen LE vector versus corresponding RE one.
!!
              iopt=2
              WRITE(LUPRI,*) 'cc_normalize: binormalize LE vector'
              call cc_normalize(length,1,work(kRfull),work(kLfull),
     &                    WORK(KEIVAL+i-1),
     &                    WORK(KWI+i-1),
     &                    isyhop,iopt)
              ioff=0
              IOPTRW = 3
              write(lupri,*)'I=',I, ' isyofe(isyhop)=', isyofe(isyhop)
              koffre = krfull
              WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     &        '- XR EIGENVALUE NO.',I,' - EIGENVALUE',
     &                        WORK(KEIVAL-1+I),
     &        '- XR EIGENVECTOR to be dumped on file:'
              !CALL CC_PRAM(WORK(KoffRe),PT1,ISYHOP,.FALSE.)
              call flshfo(lupri)
!
! upload p1 from file
!
             CALL CC_RVEC(LUPMAT,FPMAT,length,length,1,WORK(KP1))
             Tif=ddot(length,WORK(Krfull),1,WORK(KP1),1)
             WRITE (LUPRI,*)'< P1| Ri (full)>', i, Tif*etanorm,
     &       Tif
!! 
!! dump on file 
!!              call cc_dump_onfile2(length,i,WORK(KoffRe),isyhop,'RE ',
!!     &                          work(kend6), lwrk6,
!!     &                          model,ioptrw)
!!
!! ANALYSE LEFT VECTOR
!!
              WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     &        '- LPT EIGENVALUE NO.',I,' - EIGENVALUE',
     &                        WORK(KEIVAL-1+I),
     &        '- LPT EIGENVECTOR to be dumped on file:'
!!              CALL CC_PRAM(WORK(kLfull),PT1,ISYHOP,.FALSE.)
!! dump on file 
!!              call cc_dump_onfile2(length,i,WORK(KLfull),isyhop,'LE ',
!!     &                          work(kend6), lwrk6,
!!     &                          model,ioptrw)
              !
              !dot with Li*Q1 to check T_1f
              !
              IVEC = 1
              CALL CC_RDRSP('QL ',ivec,ISYHOP,IOPTRW,MODEL,
     &                WORK(KQ1),WORK(KQ1+NT1AM(ISYHOP)))
              xresid=ddot(length,WORK(KQ1),1,WORK(KP1),1)
              WRITE (LUPRI,*)'< P1|Q1 > test', i, xresid

              Tfi=ddot(length,WORK(KLfull),1,WORK(KQ1),1)
              WRITE (LUPRI,*)'< Li(full)| Q1 >', i, Tfi*csinorm,
     &        TFi
              write(lupri,'(50A)'), '--'
              write(lupri,*)'Tif*Tfi (f=',i,':) = ', 
     &                       Tif*Tfi*etanorm*csinorm, ' wf=',
     &        WORK(KEIVAL-1+I), ' a.u.'
              write(lupri,*)'------------------------------------'
          end do
      end if
!!-------------
!
      if ((locdbg).or.ABSANALYZE) then

!---------------------------------------------------------------------
! Initialization of stuff needed to be able to dump on RE and LE files
!---------------------------------------------------------------------
!
!        NCCEXCI(ISYHOP,1) = NCCEXCI(ISYHOP,1)+JKAEDE
!        MSYM set equal to NSYM (coming from outside, ccorb.h)
!
        msym=NSYM

      IF (IPRINT.GT.2) THEN
        WRITE(LUPRI,*) 'IN LANCZOS BEFORE DUMPING RE,LE VECTS'
        write(lupri,*)'MSYM =', MSYM, ' inside Lanczos (NSYM)'
        WRITE (LUPRI,*) 'NCCEXCI for singlet in LANCZOS:',
     &                  (NCCEXCI(ISYM,1),ISYM=1,msym)
        WRITE(LUPRI,*) 'NEXCI: ',NEXCI
        WRITE(LUPRI,*) 'Singlet: ',(NCCEXCI(J,1),J=1,MSYM)
        WRITE(LUPRI,*) 'Triplet: ',(NCCEXCI(J,3),J=1,MSYM)
        WRITE(LUPRI,*) 'ISYOFE:',(ISYOFE(J), J=1,MSYM)
        WRITE(LUPRI,*) 'ITROFE:',(ISYOFE(J), J=1,MSYM)
        WRITE(LUPRI,*) 'ISYEXC:',(ISYEXC(J), J=1,NEXCI)
        WRITE(LUPRI,*) 'IMULTE:',(IMULTE(J), J=1,NEXCI)
        WRITE(LUPRI,*) 'EIGVAL:',(EIGVAL(J), J=1,NEXCI)
      ENDIF
        call flshfo(lupri)
C     ----------------------------
C     set up symmetry information: DO I REALLY NEED TO DO THIS????
C     ----------------------------
!        NEXCI  = 0
!        NTRIP  = 0
!        DO ISYM = 1,MSYM
!          ISYOFE(ISYM) = NEXCI
!          ITROFE(ISYM) = ISYOFE(ISYM) + NCCEXCI(ISYM,1)
!          NEXCI        = ITROFE(ISYM) + NCCEXCI(ISYM,3)
!          NTRIP        = NTRIP        + NCCEXCI(ISYM,3)
!          DO IEX = ISYOFE(ISYM)+1, NEXCI
!             ISYEXC(IEX) = ISYM
!          END DO
!          DO IEX = ISYOFE(ISYM)+1, ITROFE(ISYM)
!             IMULTE(IEX) = 1
!          END DO
!          DO IEX = ITROFE(ISYM)+1, NEXCI
!             IMULTE(IEX) = 3
!          END DO
!        END DO
!-------------------------------------------------------------------
!       write(lupri,*) '>> Check residuals before normalizing R <<'
!-------------------------------------------------------------------
!
!      iside=+1
!      !second input entry is the number of excitations we are passing
!      WRITE(LUPRI,*) 'Check residuals for unnormalized pseudoeigenvcs'
!      call cc_check_resid(ECURR,iside,Jkaede,length,isyhop,
!     &                    WORK(KRfull),WORK(KEIVAL),
!     &                    WORK(KWI),WORK(KEND6),LWRK6,APROXR12,N2VEC,
!     &                    .true.)
!
!-------------------------------------------------------------------
!       write(lupri,*) '>> Normalize approx right eigenvecs of A <<'
!-------------------------------------------------------------------
!
!        iopt=1
!        call cc_normalize(length,jkaede,work(krfull),work(kend6),
!     &                    WORK(KEIVAL),WORK(KWI),isyhop,iopt)
!
!        IOFF = 0
!        DO I = 1,jkaede !loop on nr of available pseudoeigenvectors
!          xnorm1 = DDOT(length,WORK(KRfull+IOFF),1,WORK(KRfull+IOFF),1)
!          write(lupri,*)'unorm <R|R>', xnorm1,'& norm',sqrt(abs(xnorm1))
!          CALL DSCAL(length,(one/sqrt(abs(xnorm1))),WORK(KRfull+IOFF),1)
!          xnorm1 = DDOT(length,WORK(KRfull+IOFF),1,WORK(KRfull+IOFF),1)
!          write(lupri,*) 'Norm or R (full) after scaling', sqrt(xnorm1)
!          WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
!     *  ' - XR EIGENVALUE NO.',I, ' - EIGENVALUE',WORK(KEIVAL+I-1),
!     *  ' - XR EIGENVECTOR :'
!            CALL CC_PRP(WORK(KRfull+IOFF),
!     &                  WORK(KRfull+NT1AM(isyhop)+IOFF),isyhop,1,1)
!            CALL CC_PRAM_LANCZOS(WORK(KRfull+IOFF),PT1,ISYHOP,.TRUE.,
!     &                           work(kend6),lwrk6)
!            CALL CC_PRAM(WORK(KRfull+ioff),PT1,ISYHOP,.FALSE.)
!
!           IOFF = IOFF + length
!        END DO
C----------------------------------
!      iside=+1
!      !second input entry is the number of excitations we are passing
!      WRITE(LUPRI,*) 'Checking residuals for all pseudoeigenvectors'
!      call cc_check_resid(ECURR,iside,Jkaede,length,isyhop,
!     &                    WORK(KRfull),WORK(KEIVAL),
!     &                    WORK(KWI),WORK(KEND6),LWRK6,APROXR12,N2VEC,
!     &                    .true.)
C---------------------------------------------------------------------
        if (dump_eigfil) then
           !determine index of first eigenvalue above threshold
           istart=0
           !do i=1,nccexci(isyhop,1)
           do i=1,jkaede
              istart = i
              if (WORK(KEIVAL+I-1).ge.EIG_RANGE(1)) exit
           end do
           write(lupri,*)'Index of FIRST FREQ >= EIG_RANGE(1)',istart
           !determine how many in chosen interval
           nrexc = 0
           !do i=1,nccexci(isyhop,1)
           do i=1,jkaede
!           write(lupri,*)'WORK(KEIVAL+I-1)',WORK(KEIVAL+I-1),
!     &      ' EIG_RANGE(1):', EIG_RANGE(1), ' EIG_RANGE(2):', 
!     &        EIG_RANGE(2)
              IF ((WORK(KEIVAL+I-1).ge.EIG_RANGE(1)).and.
     &            (WORK(KEIVAL+I-1).le.EIG_RANGE(2))) then
                 nrexc = nrexc + 1
              end if
           end do
           write(lupri,*)'NR OF FREQS IN SELECTED EIG_RANGE',NrExc
!From now on Ftrans is not needed to I overwrite it
!Tranform one vector at a time
!
           do i=1,nrexc
              KRfull = kftrans
              kLfull = KRfull + length
              kend6  = KLfull + length
              lwrk6  = lwork  - kend6
              IF (LWRK6 .LT. 0) THEN
              WRITE(LUPRI,*) 'Need:', KEND6, 'Available:', LWORK
              CALL QUIT('NEWDRV: Insuff space in full RE/LE build ')
              END IF
              CALL DZERO(WORK(KRfull),length)
              CALL DZERO(WORK(KLfull),length)
        WRITE(LUPRI,*) 'CC_build_Rfull3 2 build R vector nr', i+istart-1
              !Lanczos eigevector is "jkaede" long, NOT "length"
              !IDENTIFY THE LANCZOS EIGENVECTORS TO BE TRANFORMED
              KoffR = keiver + (istart+i-1-1)*jkaede
              KoffL = keivel + (istart+i-1-1)*jkaede
!
! Change this to use the Q files instead of FQMAT...
!
              call CC_build_Rfull3(LUQMAT,FQMAT,Jkaede,length,
     &                    WORK(KRfull),work(koffR),
     &                    work(keival+istart+i-1-1),
     &                    work(kwi+istart+i-1-1),
     &                    WORK(KEND6),LWRK6)

              call CC_build_Rfull3(LUPMAT,FPMAT,Jkaede,length,
     &                    WORK(KLfull),work(koffL),
     &                    work(keival+istart+i-1-1),
     &                    work(kwi+istart+i-1-1),
     &                    WORK(KEND6),LWRK6)
!
!            normalize chosen RE vector 
!
              iopt=1
              WRITE(LUPRI,*) 'cc_normalize selected (full) RE vector'
              call cc_normalize(length,1,work(kRfull),work(krfull),
     &                    WORK(KEIVAL+istart+i-1-1),
     &                    WORK(KWI+istart+i-1-1),
     &                    isyhop,iopt)
!
!            binormalize chosen LE vector versus corresponding RE one.
!
              iopt=2
              WRITE(LUPRI,*) 'cc_normalize: binormalize LE vector'
              call cc_normalize(length,1,work(kRfull),work(kLfull),
     &                    WORK(KEIVAL+istart+i-1-1),
     &                    WORK(KWI+istart+i-1-1),
     &                    isyhop,iopt)
              ioff=0
              IOPTRW = 3
              write(lupri,*)'I=',I, ' isyofe(isyhop)=', isyofe(isyhop)
              koffre = krfull
              WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     &        '- XR EIGENVALUE NO.',Istart+I-1,' - EIGENVALUE',
     &                        WORK(KEIVAL+IStart-1+I-1),
     &        '- XR EIGENVECTOR to be dumped on file:'
              CALL CC_PRAM(WORK(KoffRe),PT1,ISYHOP,.FALSE.)
              call flshfo(lupri)
! dump on file 
              write(lupri,*)'WORK(KoffRe):', WORK(KoffRe)
              call flshfo(lupri)

              !WRITE(LUPRI,*) 'IN LANCZOS'
              !WRITE(LUPRI,*) 'NEXCI: ',NEXCI
              !WRITE(LUPRI,*) 'Singlet: ',(NCCEXCI(J,1),J=1,MSYM)
              !WRITE(LUPRI,*) 'Triplet: ',(NCCEXCI(J,3),J=1,MSYM)
              !WRITE(LUPRI,*) 'ISYOFE:',(ISYOFE(J), J=1,MSYM)
              !WRITE(LUPRI,*) 'ITROFE:',(ISYOFE(J), J=1,MSYM)
              !WRITE(LUPRI,*) 'ISYEXC:',(ISYEXC(J), J=1,NEXCI)
              !WRITE(LUPRI,*) 'IMULTE:',(IMULTE(J), J=1,NEXCI)

              NEXCI=NEXCI+1
              NCCEXCI(ISYHOP,1) = NCCEXCI(ISYHOP,1) + 1

              !WRITE(LUPRI,*) 'Singlet: ',(NCCEXCI(J,1),J=1,MSYM)
              !WRITE(LUPRI,*) 'Triplet: ',(NCCEXCI(J,3),J=1,MSYM)
              !WRITE(LUPRI,*) 'ISYOFE:',(ISYOFE(J), J=1,MSYM)
              !WRITE(LUPRI,*) 'ITROFE:',(ISYOFE(J), J=1,MSYM)
              ISYEXC(NEXCI) = isyhop
              IMULTE(NEXCI) = 1
              EIGVAL(NEXCI) = WORK(KEIVAL+IStart-1+I-1)
              iexidx = nexci

              call cc_dump_onfile2(length,iexidx,WORK(KoffRe),isyhop,
     &                          'RE ',work(kend6), lwrk6,
     &                           model,ioptrw)
!
! ANALYSE LEFT VECTOR
!
              WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     &        '- LPT EIGENVALUE NO.',Istart+I-1,' - EIGENVALUE',
     &                        WORK(KEIVAL+IStart-1+I-1),
     &        '- LPT EIGENVECTOR to be dumped on file:'
              CALL CC_PRAM(WORK(kLfull),PT1,ISYHOP,.FALSE.)
! dump on file 
              call cc_dump_onfile2(length,iexidx,WORK(KLfull),isyhop,
     &                          'LE ',
     &                          work(kend6), lwrk6,
     &                          model,ioptrw)
           end do

        if (.false.) then  

!I assume from now on that KRFUL is not available.
!determine if we have space enough for the requested RE vectors in full space
!From now on Ftrans is not needed to I overwrite it
!
           !KRfull = kfqq
           KRfull = kftrans
           kend6  = KRfull + length*NrExc
           lwrk6  = lwork  - kend6
           IF (LWRK6 .LT. 0) THEN
             WRITE(LUPRI,*) 'Need:', KEND6, 'Available:', LWORK
             CALL QUIT('Insuff space in full RE build: reduce NrExc')
           END IF
           CALL DZERO(WORK(KRfull),length*NrExc)
           WRITE(LUPRI,*) 'Call CC_build_Rfull2 2 build NrExc R vectors'
           !pas paa: each Lanczos eigevectors is "jkaede" long, NOT "length"
           KoffR = keiver + (istart - 1)*jkaede
!           write(lupri,*)'At build2, keivel=', keivel,' koffl=',koffl
!           write(lupri,*)'At build2, keival=', keival
!           write(lupri,*)'At build2, keival+istart-1=',keival+istart-1
!           write(lupri,*)'At build2, work(keival+istart-1) ',
!     &                    work(keival+istart-1)
!
! Change this to use the Q files instead of FQMAT...
!
           call CC_build_Rfull2(LUQMAT,FQMAT,NrExc,Jkaede,length,
     &                    WORK(KRfull),work(koffR),
     &                    work(keival+istart-1),
     &                    work(kwi+istart-1),
     &                    WORK(KEND6),LWRK6)
!
!          normalize chosen NrExc RE vectors 
!
           iopt=1
           !koffR = kRfull + length*(istart-1) 
           WRITE(LUPRI,*) 'cc_normalize: normalize selected RE vectors'
           call cc_normalize(length,nrexc,work(kRfull),work(krfull),
     &                    WORK(KEIVAL+istart-1),WORK(KWI+istart-1),
     &                    isyhop,iopt)

!           !check residual dump given pseudo RE vectors on file
!           koffR=KRFULL+length*(istart-1)
!           write(lupri,*)'R vecs RFULL', KRFULL, ' koffR = ',koffR
!           write(lupri,*)'Start RE vector is nr=',istart
!           write(lupri,*)'Checking selected residuals'
!           call cc_check_resid(ECURR,iside,nrexc,length,isyhop,
!     &                    WORK(koffR),
!     &                    WORK(KEIVAL+istart-1),
!     &                    WORK(KWI+istart-1),WORK(KEND6),LWRK6,
!     &                    APROXR12,N2VEC,
!     &                    .true.)
!           IOPTRW = 3
!           call cc_dump_onfile(length,nrexc,WORK(koffR),isyhop,'RE ',
!     &                          model,ioptrw)
           ioff=0
           do i=1,nrexc
              IOPTRW = 3
              write(lupri,*)'I=',I, ' isyofe(isyhop)=', isyofe(isyhop)
              koffre = krfull + (i-1)*length
              WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     &        '- XR EIGENVALUE NO.',I,' - EIGENVALUE',
     &                        WORK(KEIVAL+IStart-1+I-1),
     &        '- XR EIGENVECTOR to be dumped on file:'
!!
!!              !check residual 
!!              call dcopy(length,WORK(KoffR+IOFF),1,work(kend6),1)
!!              iside=1
!!              call CC_ATRR(ECURR,ISYHOP,ISIDE,WORK(KEND6),LWRK6,.false.,
!!     &                     dummy,.false.)
!!              krod = keival + istart - 1 + i - 1
!!              call daxpy(length,-work(krod),WORK(KoffR+IOFF),
!!     &                   1,work(kend6),1)
!!              xresid=sqrt(ddot(length,WORK(Kend6),1,WORK(Kend6),1))
!!              write(lupri,*)'Norm resid AR-wR for w=', 
!!     &                       work(krod),' =',xresid
!              CALL CC_WRRSP('RE ',I,ISYHOP,IOPTRW,MODEL,WORK(KEND6),
!     *                   WORK(KOFFRE),
!     *                   WORK(KOFFRE+NT1AM(ISYHOP)),
!     *                   WORK(KEND6),LWRK6)
!!              CALL CC_WRRSP('RE ',I,ISYHOP,IOPTRW,MODEL,WORK(KEND6),
!!     *                   WORK(KoffR+IOFF),
!!     *                   WORK(KoffR+IOFF+NT1AM(ISYHOP)),
!!     *                   WORK(KEND6),LWRK6)
!              !write(lupri,*)'Below RE_i=', i, ' of address= ',koffR+ioff
              CALL CC_PRAM(WORK(KoffRe),PT1,ISYHOP,.FALSE.)
              !IOFF = IOFF + length
           end do
!
           end if
           !determine if we have space enough for the corresponding LE vectors

           if (.false.) then
           KLfull = kftrans
           kend6  = KLfull + length*NrExc
           lwrk6  = lwork  - kend6
           IF (LWRK6 .LT. 0) THEN
              WRITE(LUPRI,*) 'Need:', KEND6, 'Available:', LWORK
              CALL QUIT('Insuff work space in LE build: reduce NrExc')
           END IF
           CALL DZERO(WORK(KLfull),length*NrExc)
           WRITE(LUPRI,*) 'Call CC_build_Rfull2 2 build NrExc L vectors'
           !pas paa: each Lanczos eigevectors is "jkaede" long, NOT "length"
           KoffL = keivel + (istart - 1)*jkaede
!           write(lupri,*)'At build2, keivel=', keivel,' koffl=',koffl
!           write(lupri,*)'At build2, keival=', keival
!           write(lupri,*)'At build2, keival+istart-1=',keival+istart-1
!           write(lupri,*)'At build2, work(keival+istart-1) ',
!     &                    work(keival+istart-1)
           call CC_build_Rfull2(LUPMAT,FPMAT,NrExc,Jkaede,length,
     &                    WORK(KLfull),work(koffL),
     &                    work(keival+istart-1),
     &                    work(kwi+istart-1),
     &                    WORK(KEND6),LWRK6)
!
!           binormalize chosen NrExc LE vectors versus corresponding RE ones.
!
           iopt=2
           !koffR = kRfull + length*(istart-1) 
           WRITE(LUPRI,*) 'cc_normalize: binormalize LE vectors'
           call cc_normalize(length,nrexc,work(kRfull),work(klfull),
     &                    WORK(KEIVAL+istart-1),WORK(KWI+istart-1),
     &                    isyhop,iopt)
           do i=1,NrExc
              koffl = klfull + (i-1)*length
!!              koffL = kLfull + ioff
!!              koffR = kRfull + length*(istart-1) + ioff
!!              xnorm1=DDOT(length,WORK(KOFFL),1,WORK(KOFFR),1)
!!              write(lupri,*) '<LP|QR> Norm before scaling', xnorm1
!!              CALL DSCAL(length,(one/(xnorm1)),WORK(KOFFL),1)
!!              xnorm1=DDOT(length,WORK(KoffL),1,WORK(KOFFR),1)
!!              write(lupri,*) '<LP|QR> Norm after scaling', xnorm1

              WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     &        '- LPT EIGENVALUE NO.',I,' - EIGENVALUE',
     &                        WORK(KEIVAL+IStart-1+I-1),
     &        '- LPT EIGENVECTOR to be dumped on file:'
              CALL CC_PRAM(WORK(KoffL),PT1,ISYHOP,.FALSE.)
!!              CALL CC_PRAM_LANCZOS(WORK(KOFFL),PT1,1,.TRUE.,
!!     &                             work(kend6),lwrk6)
!!              write(lupri,*) 'Writing on file LE vector nr ', I
              CALL CC_WRRSP('LE ',I,ISYHOP,IOPTRW,MODEL,WORK(KEND6),
     &                       WORK(Koffl),
     &                       WORK(Koffl+NT1AM(ISYHOP)),
     &                       WORK(KEND6),LWRK6)
              call flshfo(lupri)
              !IOFF = IOFF + length
          end do
!
!          iside=-1
!          write(lupri,*) 'Checking left residual'
!          call cc_check_resid(ECURR,iside,nrexc,length,isyhop,
!     &                    WORK(kLfull),
!     &                    WORK(KEIVAL+istart-1),
!     &                    WORK(KWI+istart-1),WORK(KEND6),LWRK6,
!     &                    APROXR12,N2VEC,
!     &                    .true.)
         ! write(lupri,*) 'the NrExc Lfull MATRIX'
         !CALL OUTPUT(WORK(KLfull),1,length,1,NrExc,length,NrExc,1,LUPRI)
         end if
        end if  !dump_eigfil
!!
       end if   !analyze

C-------------------------------------------------------------
C     finished
C-------------------------------------------------------------

      end if !isyma=isymb

      !CALL WCLOSE2(LUTTRMAT,FTTRMAT,'KEEP')
      CALL WCLOSE2(LUQMAT,FQMAT,'KEEP')
      CALL WCLOSE2(LUPMAT,FPMAT,'KEEP')

      CALL QEXIT('CC_LANCZOS_NEWDRV')
      RETURN
      END
