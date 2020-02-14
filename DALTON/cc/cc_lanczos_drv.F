C  /* Deck cc_lanczos_drv */
C--------------------------------------------------------------------
      SUBROUTINE CC_LANCZOS_drv(WORK,LWORK,APROXR12)
C--------------------------------------------------------------------
C     Purpose: Noddy code to build the Lanczos recursive chain
C     Handling of RHS vectors is for the time being done
C     like in CC_XETST
C     Sonia Coriani & Kristian Sneskov, August 2010
C     FIXME: RESTRUCTURE THE DUMP-ON-FILE FOR CCSDR3
C     SYMMETRY inserted in noddy way.
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

* local parameters:
      CHARACTER MSGDBG*(18)
      PARAMETER (MSGDBG='[debug] CC_LANCZOS_DRV> ')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .true.)
      INTEGER MXXETRAN, MXVEC, IPRLEE, MXFQTRAN
      PARAMETER (MXXETRAN = 20, MXVEC = 10, IPRLEE=1, MXFQTRAN=50)
      integer ifile, kRfull, kLfull, KoffR,KoffL,istart,NrExc,koffre
      INTEGER LWORK, LENGTH,Jkaede, kFtrans
      INTEGER kLe1,lensym, length_check,msym,NTRIP,IEX,krod,lentot
#if defined (SYS_CRAY) 
      REAL WORK(LWORK), ALPHA(JCHAIN), BETA(jChain), xgamma(jchain)
      REAL FREQ, ECURR, rxs, test
      REAL DDOT
      REAL ZERO, ONE, TWO, FOUR, FIVE, Q1norm, factor
      REAL ovlpLR, xnorm1, ovlpRR, ovlpLL
      real sum11,csinorm, etanorm
#else 
      !sonia: enlarge alpha by 1 to check residuals
      DOUBLE PRECISION WORK(LWORK), ALPHA(JCHAIN+1), BETA(jChain)
      DOUBLE PRECISION xgamma(jchain), rxs, test
      DOUBLE PRECISION FREQ, ECURR, PT1, xresid
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, FIVE, Q1norm, factor
      DOUBLE PRECISION thrshold
      DOUBLE PRECISION ovlpLR, xnorm1
      DOUBLE PRECISION ovlpRR, ovlpLL, sum11, csinorm,etanorm
      
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
      INTEGER IFQTRAN(MXDIM_FTRAN,JCHAIN), NFQTRAN, ICHAIN
      INTEGER ISYAMP,ISYM,ksym
      INTEGER ISYM0, ISYHOP, IOPT, IOPER, ITRAN, IVEC, IDUM, ISYETA
      INTEGER IOPTRW,III
      INTEGER ITEST, IORDER, ILSTETA, ISIGN
      INTEGER KDUM, LWRK1, LWRK2, LWRK3,LWRK4,LWRK5,LWRK6,LWRK61
      INTEGER KEND0,KEND1,KEND2,KEND3,KEND4,KEND5,KEND6,KEND61
      INTEGER KTTRI, KWI, KEIVAL, KEIVER, KEIVEL
      INTEGER IDLSTL,ISYCTR,IDXR1,INUM,IREAL
      INTEGER KETA1, KETA2, KRHO1, KRHO2
      INTEGER IRELAX, KCTR1, KCTR2
      INTEGER KETA12, kpnewa1, ks1, kr1, kpnewa2, kqnew1, kqnew2
      INTEGER KPNEW1, KPNEW2, KAQNEW1, KAQNEW2, kqold1, kqold2, kpold1
      INTEGER kpold2, ICPLX,IDXVEC,IOFF,JOFF
      integer iopera,ioperb,isyma,isymb
      CHARACTER*(8) LABELA,LABELB
      logical lorxa,lorxb
*
* files expected to exist in RHS & LHS transformations
*
      INTEGER LUQ11,LUQ12, LUAQ11, LUAQ12, LUDUM
      CHARACTER*8 FRHS1,FRHS2,FRHO1,FRHO2,CDUM
      CHARACTER FILEX*(10), LISTI*(4)
      PARAMETER (FRHS1 ='CCXS_Q11' ,FRHS2 ='CCXS_Q12')
      PARAMETER (FRHO1 ='CCXSAQ11',FRHO2 ='CCXSAQ12')
      CHARACTER*8 FP11,FP12,FP1A1,FP1A2
      INTEGER LUP11,LUP12, LUP1A1, LUP1A2
      PARAMETER (FP11 ='CCXS_P11' ,FP12 ='CCXS_P12')
      PARAMETER (FP1A1 ='CCXSP1A1',FP1A2 ='CCXSP1A2')
*
*     Files on which we dump the q and p(T) matrices
*
      integer LUQMAT, LUPMAT
      CHARACTER*8 FQMAT,FPMAT
      PARAMETER (FQMAT ='CCLAN_Q' ,FPMAT ='CCLAN_PT')

* external functions (lists):
      INTEGER ILSTSYM
      INTEGER IROPER
      INTEGER IR1TAMP
      INTEGER IRHSR1
      INTEGER IETA1
      INTEGER IQLLIST

      CALL QENTER('CC_LANCZOS_DRV')

      CALL IZERO(NCCEXCI,3*8)

*---------------------------------------------------------------------*
* set up information for rhs vectors we want
*---------------------------------------------------------------------*
*  number of simultaneous transformations:

      NXETRAN = 1
      NFQTRAN = 0  !transformations of QL vectors with Fmatrix
                   !it will be equal to Jkaede at the end

      !initialize the list of F*q transformations
      DO I=1,JCHAIN
         DO J = 1, MXDIM_FTRAN
          IFQTRAN(J,I) = 0
         END DO
      END DO
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
      if (isyma.eq.isymb) then
      isyhop  = isymb
      label   = labela
      ltwo    = .false.
      lorx    = .false.
      write(lupri,*)'isyhop',isyhop,' label=', label
      LOPROPN = .TRUE.
      LO1OPN  = .TRUE.
      LX1OPN  = .TRUE.
      LQLOPN  = .TRUE.
      !ISYHOP  = IROPER(LABELO,ksym)  !symmetric operator (dipole)
      !write(lupri,*)'ksym=',ksym,' isyhop=', isyhop
 
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
      IXETRAN(1,1) = IROPER(LABEL,ISYHOP)
      IXETRAN(2,1) = 0
      IXETRAN(3,1) = IRHSR1(LABEL,LORX,FREQ,ISYHOP)
      IXETRAN(4,1) = IETA1(LABEL,LORX,FREQ,ISYHOP)
      IXETRAN(5,1) = IRELAX
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
      WRITE(LUPRI,*) 'LANCZOS: Chain length by input, JCHAIN=', JCHAIN
      Jkaede = Jchain
      IF (JCHAIN.GT.LENGTH) THEN
         WRITE(LUPRI,*) 'WARNING: JCHAIN EXCEED FULL SPACE DIMENSION'
         Jkaede = length
         WRITE(LUPRI,*) 'WE RESET IT TO dim{AMAT} =   ', jkaede
      END IF

      ! Generate index of vector Q on file
      ! set the IFQTRAN array for the first Qkaede vector
      !
      do i=1,jkaede
         ICHAIN = i
         !ifile = IQLLIST(LABEL,LORX,ICHAIN,FREQ,ISYHOP)
         IFQTRAN(1,i) = 0
         IFQTRAN(2,i) = IQLLIST(LABEL,LORX,ICHAIN,FREQ,ISYHOP)
         IFQTRAN(3,i) = IFQTRAN(2,i)
      end do

      ! disallow again extension in the vector lists
      LOPROPN = .FALSE.
      LO1OPN  = .FALSE.
      LX1OPN  = .FALSE.
      LQLOPN  = .FALSE.

      if (locdbg) then
        CALL AROUND('CC_LANCZOS: test of call to CC_XIETA setup')
        WRITE (LUPRI,*) 'ITEST  =',ITEST
        WRITE (LUPRI,*) 'ISYHOP =',ISYHOP
        WRITE (LUPRI,*) 'LABEL  =',LABEL
        WRITE (LUPRI,*) 'LTWO   =',LTWO
        WRITE (LUPRI,*) 'LORX   =',LORX
        WRITE (LUPRI,*) 'FREQ   =',FREQ
        WRITE (LUPRI,*) 'IROPER =',IXETRAN(1,1)
        WRITE (LUPRI,*) 'ILEFT  =',IXETRAN(2,1)
        WRITE (LUPRI,*) 'IRHSR1 =',IXETRAN(3,1)
        WRITE (LUPRI,*) 'IETA1  =',IXETRAN(4,1)
        WRITE (LUPRI,*) 'IRELAX =',IXETRAN(5,1)
        WRITE (LUPRI,*) 'JCHAIN =', JCHAIN
        WRITE (LUPRI,*) 'JKAEDE =', JKAEDE

        CALL AROUND('CC_LANCZOS_DRV: test of call to CC_FMATRIX setup')
        !do i=1,jkaede
        !WRITE (LUPRI,*) 'IL0    =',IFQTRAN(1,i)
        !WRITE (LUPRI,*) 'IQL    =',IFQTRAN(2,i)
        !WRITE (LUPRI,*) 'IFQL   =',IFQTRAN(3,i)
        !write (lupri,*)'---'
        !end do
      end if

*---------------------------------------------------------------------*
* call CC_XIETA to calculate the Xi and Eta vectors:
* they are dumped directly on file withing XIETA
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
     &              .FALSE.,0,WORK(KEND1), LWRK1 )
!
!     open files for Q and P matrices
!     FIXME: replace all this with call to ATTR
!
      LUQMAT = -1
      LUPMAT = -1
      CALL WOPEN2(LUQMAT,FQMAT,64,0)
      CALL WOPEN2(LUPMAT,FPMAT,64,0)

      LUQ11=-1
      LUQ12=-1 
      LUAQ11=-1 
      LUAQ12=-1
      LUP11=-1
      LUP12=-1 
      LUP1A1=-1 
      LUP1A2=-1

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
      LWRK1  = LWORK - KEND1
      if (LWRK1.LT.0) then
         WRITE(LUPRI,*) 'Need:', KEND1, 'Available:', LWORK
         CALL QUIT('Insufficient memory in CC_LANCZOS_DRV (1)')
      end if
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
      IF (LWRK2.LT.0) THEN
         WRITE(LUPRI,*) 'Need:',KEND2, 'Available:', LWORK
         CALL QUIT('Insufficient memory in CC_LANCZOS_DRV (2)')
      END IF
C     -------------------------------------------------
C     read in the Csi^Z vector from CC_XIETA and print
C     -------------------------------------------------
      NSIMTR = NXETRAN
      ECURR  = ZERO
      ITRAN = 1
      IVEC = IXETRAN(3,ITRAN)
      !WRITE (LUPRI,*) 'CC_LANCZOS: IVEC = ',IVEC
      CALL CC_RDRSP('O1 ',IVEC,ISYHOP,IOPTRW,MODEL,
     &                WORK(KQnew1),WORK(KQnew2))

      if (locdbg) then
        CALL AROUND('analytical Xi vector:')
        !CALL CC_PRP(WORK(KQnew1),WORK(KQnew2),ISYHOP,1,N2VEC)
      end if
      q1norm = sqrt(abs(ddot(length,WORK(KQnew1),1,WORK(KQnew1),1)))
      write(lupri,*)'q1norm: ', q1norm
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
      !save a copy for later use (resp amplitudes)
      !call dcopy(length,work(KQnew1),1,work(kcsi),1)
!
C     --------------------------------
C     Compute the A*q1 transformation
C     --------------------------------
!     Define file names on which to write the transformed vectors
!     It looks like singles and doubles part are kept on different files...

      IVEC  = 1
      !open files
      CALL CC_FILOP(FRHO1,LUAQ11,FRHO2,LUAQ12,CDUM,LUDUM,
     *              FRHS1,LUQ11,FRHS2,LUQ12,CDUM,LUDUM,
     *              CDUM,LUDUM,CDUM,LUDUM)

      CALL CC_WVEC(LUQ11,FRHS1,NT1AM(ISYHOP),
     *             NT1AM(ISYHOP),1,WORK(KQnew1))
      if (.not.CCS) then
         CALL CC_WVEC(LUQ12,FRHS2,NT2AM(ISYHOP),
     *                   NT2AM(ISYHOP),1,WORK(KQnew2))
      end if

      !dump q1 onto the QMAT file
      !if (locdbg) then
        !CALL AROUND('analytical Q1 vector:')
        !CALL CC_PRP(WORK(KQnew1),WORK(KQnew2),ISYHOP,1,N2VEC)
      !end if

      CALL CC_WVEC(LUQMAT,FQMAT,length,length,1,WORK(KQnew1))
!
!try to dump into the Qvec_1 file as well for later use in the
!Ftran list
!
      IFILE = IQLLIST(LABEL,LORX,1,FREQ,ISYHOP)

      write(lupri,*)'BEFORE CC_WRRSP: IFILE=', 
     &               IFILE, ' & ifqt=', IFqtran(2,1)
      WRITE (LUPRI,*) 'LABEL  =', LABEL
      WRITE (LUPRI,*) 'LORX   =',LORX
      WRITE (LUPRI,*) 'FREQ   =',FREQ
      WRITE (LUPRI,*) 'ISYHOP =',ISYHOP

!      call cc_rwpre('QL ',IFILE,ISYHOP,LISTI,IDXVEC,FILEX)
!      write(lupri,*)' FILEX: ', filex
!      write(lupri,*)' listi: ', listi
!      write(lupri,*)' idxvec: ', idxvec
!      call flshfo(lupri)

      IFILE = IQLLIST(LABEL,LORX,1,FREQ,ISYHOP)
      CALL CC_WRRSP('QL ',IFILE,ISYHOP,IOPTRW,MODEL,WORK(KEND2),
     &           WORK(KQnew1),WORK(KQnew2),WORK(KEND2),LWRK2)

      call dcopy(length,work(KQnew1),1,work(KAQnew1),1)
      iside=1
      call CC_ATRR(ECURR,ISYHOP,ISIDE,WORK(KAQnew1),LWRK2,.false.,Cdum,
     &                   APROXR12,.false.)
!
!      CALL DZERO(WORK(KAQnew1),NT1AM(isyhop))
!      if (.not.ccs) CALL DZERO(WORK(KAQnew2),NT2AM(isyhop))
!
! Compute the transformed vectors A*q1=A*csi^X
!
!      WRITE (LUPRI,*) 'just before cc_trdrv'
!
!      iside=+1
!      nsim = 1
!      ist =  1  !A*q
!      CALL CC_TRDRV(ECURR,FRHO1,LUAQ11,FRHO2,LUAQ12,CDUM,LUDUM,
!     *                    FRHS1,LUQ11,FRHS2,LUQ12,CDUM,LUDUM,
!     *                    .false.,.false.,ITRAN,FREQ,cdum,ludum,
!     *                    cdum,ludum,
!     *                    ISIDE,IST,NSIM,WORK(KEND2),LWRK2,APROXR12)
!      WRITE (LUPRI,*) 'just after cc_trdrv'
!
!      CALL CC_RVEC(LUAQ11,FRHO1,NT1AM(isyhop),NT1AM(isyhop),1,
!     &             WORK(KAQnew1))
!      WRITE (LUPRI,*) 'just after cc_rvec'
!      if (.not.ccs) 
!     &   CALL CC_RVEC(LUAQ12,FRHO2,NT2AM(isyhop),NT2AM(isyhop),1,
!     &        WORK(KAQnew2))

      if (locdbg) then
        Write(lupri,*)'---AQ vector reread from file---'
        !CALL CC_PRP(WORK(KAQnew1),WORK(KAQnew2),1,1,N2VEC)
      end if
!
! NOW ETA^X VECTOR IN PRINCIPLE ANOTHER SYMMETRY
!

      !Read in Xi^Z as starting guess for p1 vector
      ISYCTR = ILSTSYM(LISTL,IDLSTL)
      ISYETA = MULD2H(ISYCTR,ISYHOP)

      KPnew1  = KEND2
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
C     -------------------------------------------------
      ITRAN = 1
      IVEC = IXETRAN(4,ITRAN) 
      WRITE (LUPRI,*) 'CC_LANCZOS: IVEC = ',IVEC
      CALL CC_RDRSP('X1 ',IVEC,ISYETA,IOPTRW,MODEL,
     &                WORK(KPnew1),WORK(KPnew2))

      if (locdbg) then
        CALL AROUND('P1 vector (eta^X) before biorthonorm:')
        !CALL CC_PRP(WORK(KPnew1),WORK(KPnew2),ISYETA,1,N2VEC)
        !recompute using different drive to make sure we got it right
!        call DZERO(NT1AM(ISYETA)+NT2AM(ISYETA),WORK(KPnew1))
!        call CC_ETAC(ISYETA,LABEL,WORK(KPnew1),'L0',1,0,
!     *                   DUMMY,WORK(KEND3),LWRK3)
!        CALL AROUND('P1 vector (eta^X) from cc_etac:')
!        CALL CC_PRP(WORK(KPnew1),WORK(KPnew2),ISYETA,1,N2VEC)
      end if
      !save a copy for later use (resp amplitudes)
      !call dcopy(length,work(KPnew1),1,work(keta),1)
      !binormalize vs q1
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
      !dump p1 onto the PMAT file
      IVEC  = 1
      CALL CC_WVEC(LUPMAT,FPMAT,length,length,ivec,WORK(KPnew1))

      write(lupri,*)'LANCZOS: p1*q1 norm', 
     &             ddot(length,WORK(KPNEW1),1,WORK(KQNEW1),1)
!
! Compute first alpha(1)=<p1|A*q>
!
      alpha(1) = ddot(length,WORK(KPNEW1),1,WORK(KAQNEW1),1)
      write(lupri,*)'LANCZOS: first alpha in chain', alpha(1)

      IVEC  = 1
      !open files for p^T*A transformation
      CALL CC_FILOP(FP1A1,LUP1A1,FP1A2,LUP1A2,CDUM,LUDUM,
     *              FP11,LUP11,FP12,LUP12,CDUM,LUDUM,
     *              CDUM,LUDUM,CDUM,LUDUM)

      CALL CC_WVEC(LUP11,FP11,NT1AM(ISYETA),
     *             NT1AM(ISYETA),ivec,WORK(KPnew1))
      CALL DZERO(WORK(KPnewA1),NT1AM(ISYETA))
      if (.not.CCS) then
        CALL CC_WVEC(LUP12,FP12,NT2AM(ISYETA),
     *                   NT2AM(ISYETA),ivec,WORK(KPnew2))
        CALL DZERO(WORK(KPnewA2),NT2AM(ISYETA))
      end if


      iside=-1
      nsim = 1
      ist =  1  !p^T A
      CALL CC_TRDRV(ECURR,FP1A1,LUP1A1,FP1A2,LUP1A2,CDUM,LUDUM,
     *                    FP11,LUP11,FP12,LUP12,CDUM,LUDUM,
     *                    .false.,.false.,ITRAN,FREQ,cdum,ludum,
     *                    cdum,ludum,
     *                    ISIDE,IST,NSIM,WORK(KEND3),LWRK3,APROXR12)

      CALL CC_RVEC(LUP1A1,FP1A1,NT1AM(ISYETA),NT1AM(ISYETA),1,
     &             WORK(KPnewA1))
      if (.not.CCS)
     &   CALL CC_RVEC(LUP1A2,FP1A2,NT2AM(isyeta),NT2AM(isyeta),1,
     &                WORK(KPnewA2))
      !Write(lupri,*)'---p*A vector reread from file---'
      !CALL CC_PRP(WORK(KPnewA1),WORK(KPnewA2),1,1,N2VEC)

      KR1 = KEND3
      KS1 = Kr1 + length
      KEND4 = Ks1 + length
      LWRK4 = LWORK - KEND4
      IF (LWRK4.LT.0) THEN
           WRITE(LUPRI,*) 'Need:', KEND4, 'Available:', LWORK
           CALL QUIT('Insufficient memory in CC_LANCZOS.')
      END IF

      CALL DZERO(WORK(KR1),length)
      call dcopy(length,work(KAQnew1),1,work(KR1),1)
      CALL DAXPY(length,-alpha(1),WORK(KQnew1),1,WORK(KR1),1)

      CALL DZERO(WORK(KS1),length)
      call dcopy(length,work(KPNEWA1),1,work(KS1),1)
      CALL DAXPY(length,-alpha(1),WORK(KPNEW1),1,WORK(KS1),1)
      rxs = ddot(length,WORK(KR1),1,WORK(KS1),1)
      beta(1)   = sqrt(abs(rxs))
      xgamma(1) = sign(beta(1),rxs)
      write(lupri,'(A,3f16.10)') 'LANCZOS 1 : rxs, beta, gamma:', 
     &                  rxs, beta(1),xgamma(1)
!      write(lupri,*) 'LANCZOS: sxr:', 
!     &       ddot(length,WORK(Ks1),1,WORK(Kr1),1)

      if (abs(rxs).lt.thrshold) call quit('rxs at start = 0')

C Prepare to make chain. Allocate workspace - at all times we only
C need 2 (of each kind) vectors in memory. The rest are stored on disk.

      KQold1 = KEND4
      KQold2 = KQold1 + NT1AM(ISYHOP)
      KPold1 = KQold2 + NT2AM(ISYHOP)
      KPold2 = KPOLD1 + NT1AM(ISYeta)
      KEND5  = KPold2 + NT2AM(ISYeta)
      LWRK5  = LWORK  - KEND5

      IF (LWRK5.LT.0) THEN
           WRITE(LUPRI,*) 'Need:', KEND5, 'Available:', LWORK
           CALL QUIT('Insufficient memory in CC_LANCZOS, 5.')
      END IF
      
      length_check = jkaede

      do i=2,jkaede+1
        
        CALL DCOPY(length,work(KQnew1),1,work(KQold1),1)
        CALL DCOPY(length,work(KPnew1),1,work(KPold1),1)
        CALL DCOPY(length,work(KR1),1,work(KQnew1),1)
        factor = ONE/beta(i-1)
        CALL DSCAL(length,factor,WORK(KQnew1),1)
        CALL DCOPY(length,work(KS1),1,work(KPnew1),1)
        factor = ONE/xgamma(i-1)
        CALL DSCAL(length,factor,WORK(KPnew1),1)
        !In principle should get x_i now
        !dump qnew onto the QMAT file, idem for p
        CALL CCBIORT(LUQMAT,FQMAT,
     &                LUPMAT,FPMAT,i,i-1,length,
     *                work(KQNEW1),work(KPNEW1),
     &                isyhop,1.0d0-10,WORK(KEND5),LWRK5)
      if (i.ne.(jchain+1)) then
        CALL CC_WVEC(LUQMAT,FQMAT,length,
     *               length,i,WORK(KQnew1))

        IFILE = IQLLIST(LABEL,LORX,i,FREQ,ISYHOP)
        CALL CC_WRRSP('QL ',IFILE,ISYHOP,IOPTRW,MODEL,WORK(KEND5),
     &          WORK(KQnew1),WORK(KQnew2),WORK(KEND5),LWRK5)

        CALL CC_WVEC(LUPMAT,FPMAT,length,
     *               length,i,WORK(KPnew1))

C     --------------------------------
C     Compute the A*Qnew transformation
C     --------------------------------

        CALL CC_WVEC(LUQ11,FRHS1,NT1AM(ISYHOP),
     *             NT1AM(ISYHOP),1,WORK(KQnew1))
        CALL DZERO(WORK(KAQnew1),NT1AM(ISYHOP))
        if (.not.CCS) then
          CALL CC_WVEC(LUQ12,FRHS2,NT2AM(ISYHOP),
     *                   NT2AM(ISYHOP),1,WORK(KQnew2))
          CALL DZERO(WORK(KAQnew2),NT2AM(ISYHOP))
        end if
C
        iside=+1
        nsim = 1
        ist =  1
        CALL CC_TRDRV(ECURR,FRHO1,LUAQ11,FRHO2,LUAQ12,CDUM,LUDUM,
     *                    FRHS1,LUQ11,FRHS2,LUQ12,CDUM,LUDUM,
     *                    .false.,.false.,ITRAN,FREQ,cdum,ludum,
     *                    cdum,ludum,
     *                    ISIDE,IST,NSIM,WORK(KEND5),LWRK5,APROXR12)
C
        CALL CC_RVEC(LUAQ11,FRHO1,NT1AM(isyhop),NT1AM(isyhop),1,
     &               WORK(KAQnew1))
        if (.not.CCS) 
     &   CALL CC_RVEC(LUAQ12,FRHO2,NT2AM(isyhop),NT2AM(isyhop),1,
     &                WORK(KAQnew2))
        !Write(lupri,*)'---AQnew vector reread from file, ichain:', i
        !CALL CC_PRP(WORK(KAQnew1),WORK(KAQnew2),1,1,N2VEC)

        CALL CC_WVEC(LUP11,FP11,NT1AM(ISYeta),
     *               NT1AM(ISYeta),1,WORK(KPnew1))
        CALL DZERO(WORK(KPnewA1),NT1AM(isyeta))
        if (.not.CCS) then
          CALL CC_WVEC(LUP12,FP12,NT2AM(ISYeta),
     *               NT2AM(ISYeta),1,WORK(KPnew2))
          CALL DZERO(WORK(KPnewA2),NT2AM(isyeta))
        end if

        iside=-1
        nsim = 1
        ist =  1
        CALL CC_TRDRV(ECURR,FP1A1,LUP1A1,FP1A2,LUP1A2,CDUM,LUDUM,
     *                    FP11,LUP11,FP12,LUP12,CDUM,LUDUM,
     *                    .false.,.false.,ITRAN,FREQ,cdum,ludum,
     *                    cdum,ludum,
     *                    ISIDE,IST,NSIM,WORK(KEND5),LWRK5,APROXR12)

        CALL CC_RVEC(LUP1A1,FP1A1,NT1AM(isyeta),NT1AM(isyeta),1,
     &               WORK(KPnewA1))
        if (.not.CCS) 
     &    CALL CC_RVEC(LUP1A2,FP1A2,NT2AM(isyeta),NT2AM(isyeta),1,
     &         WORK(KPnewA2))
        !Write(lupri,*)'---pnew*A vector reread from file---'
        !CALL CC_PRP(WORK(KPnewA1),WORK(KPnewA2),1,1,N2VEC)

        alpha(i) = ddot(length,WORK(KPNEW1),1,WORK(KAQNEW1),1)
        write(lupri,*)'At ichain:',i,' alpha(ichain)=', alpha(i)
        !
        CALL DZERO(WORK(KR1),length)
        call dcopy(length,work(kAQnew1),1,work(kr1),1)
        CALL DAXPY(length,-alpha(i),WORK(KQnew1),1,WORK(KR1),1)
        CALL DAXPY(length,-xgamma(i-1),WORK(KQold1),1,WORK(KR1),1)

        CALL DZERO(WORK(KS1),length)
        call dcopy(length,work(kPnewA1),1,work(ks1),1)
        CALL DAXPY(length,-alpha(i),WORK(KPnew1),1,WORK(KS1),1)
        CALL DAXPY(length,-beta(i-1),WORK(KPold1),1,WORK(KS1),1)

        rxs = ddot(length,WORK(Kr1),1,WORK(KS1),1)
        !WRITE(LUPRI,*) 'rxs=',rxs,' for ichain=',i
        beta(i) = sqrt(abs(rxs))
        xgamma(i) = sign(beta(i),rxs)
        write(lupri,'(A,i5,3f16.10)') 'ichain,rxs, beta, gamma:', 
     &                  i, rxs,beta(i),xgamma(i)
        if (abs(rxs).lt.thrshold) then
           write(lupri,*)'|rxs|< thrshold at ichain=', i
           write(lupri,*)'Chain ends here, jchain reset to ', i
           length_check = i
           go to 111
        end if   
       end if
      END DO

 111  jkaede = length_check

      !reset F mat list????
      !DO I=1,JCHAIN
      !   DO J = 1, MXDIM_FTRAN
      !    IFQTRAN(J,I) = 0
      !   END DO
      !END DO
      !do i=1,jkaede
      !   ICHAIN = i
      !   IFQTRAN(1,i) = 0
      !   IFQTRAN(2,i) = IQLLIST(LABEL,LORX,ICHAIN,FREQ,ISYHOP)
      !   IFQTRAN(3,i) = IFQTRAN(2,i)
      !end do

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

      write(lupri,*)' AT T build jkaede =', jkaede

      KTTRI  = KEND5                   !The T matrix
      KEIVAL = KTTRI  + JKAEDE*JKAEDE  !real eigenvalues
      KWI    = KEIVAL + JKAEDE         !imaginary eigenvalues
      KEIVER = KWI    + jkaede         !R eigenvectors
      KEIVEL = KEIVER + jkaede*jkaede   !L eigenvectors
      KEND6  = KEIVEL + jkaede*jkaede
      LWRK6  = LWORK - KEND6
      if (LWRK6.lt.0) then
        write(lupri,*)'Need:', kend6, ' Available:', lwork
        call quit('Not enough memory in Lanczos')
      end if
      CALL DZERO(WORK(KTTRI),JKAEDE*JKAEDE)
      CALL DZERO(WORK(KEIVAL),JKAEDE)
      CALL DZERO(WORK(KEIVER),JKAEDE*JKAEDE)
      CALL DZERO(WORK(KEIVEL),JKAEDE*JKAEDE)
      CALL DZERO(WORK(KWI),JKAEDE)
C
C Buid tridiagonal T matrix
C
      do i=1,jkaede
        WORK(KTTRI+i-1+jkaede*(i-1))=alpha(i)
        IF (i.lt.jkaede) WORK(KTTRI+i-1+jkaede*i)     = xgamma(i)
        IF (i.gt.1)      WORK(KTTRI+i-1+jkaede*(i-2)) = beta(i-1)
      end do
      CALL AROUND('--- The T tridiagonal matrix --- ') 
      CALL OUTPUT(WORK(KTTRI),1,jkaede,1,jkaede,jkaede,jkaede,1,LUPRI)
C
C Diagonalize T matrix and find both R and L
C
      !husk dgesv
      IERR=0
      CALL DGEEV('V','V',jkaede,WORK(KTTRI),jkaede,
     &           WORK(KEIVAL),WORK(KWI),WORK(KEIVEL), jkaede, 
     &           WORK(KEIVER),jkaede,
     &           WORK(KEND6), LWRK6, IERR)
      WRITE(LUPRI,*) 'Printing INFO after call to DGEEV: ', IERR

      IF (LOCDBG .OR. IERR .NE. 0) THEN
          WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES real part:'
          WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
          !CALL OUTPUT(WORK(KEIVAL),1,jkaede,1,1,jkaede,jkaede,1,LUPRI)
          WRITE (LUPRI,'(/A)')
     *    ' REDUCED EIGENVALUES imaginary part:'
          WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
          !CALL OUTPUT(WORK(KWI),1,jkaede,1,1,jkaede,jkaede,1,LUPRI)
      END IF
      IF ( IERR.NE.0 ) THEN
         WRITE(LUPRI,'(/A,I5)')
     *   ' EIGENVALUE PROBLEM NOT CONVERGED IN RG/DGEEV, IERR =',IERR
         CALL QUIT(' LANCZOS: EIGENVALUE EQUATION NOT CONVERGED ')
      END IF
C
      !WRITE (LUPRI,*) ' WARNING: I SKIP REORDER FTB SONIA'
      !WRITE (LUPRI,*) ' REORDER REINTRODUCED'
      !reorder and Normalize both KEIVER and KEIVAL
      CALL RGORD2(jkaede,jkaede,WORK(KEIVAL),WORK(KWI),
     &            WORK(KEIVER),WORK(KEIVEL),.FALSE.)

      ICPLX = 0
      DO I=1,jkaede
         IF (WORK(KWI+I-1) .NE. ZERO) THEN
               ICPLX = ICPLX + 1
               WRITE(LUPRI,'(I10,1P,2D15.8,A/)') I,WORK(KEIVAL+I-1),
     *               WORK(KWI+I-1),
     *            ' *** CCRED  WARNING **** COMPLEX VALUE.'
         END IF
      END DO
C
      IF (IPRLEE .GE.10.OR. ICPLX .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' REDUCED right EIGENVALUES real part:'
         CALL OUTPUT(WORK(KEIVAL),1,jkaede,1,1,jkaede,1,1,LUPRI)
      ENDIF
C
      IF (IPRLEE .GE.10.OR. ICPLX .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' REDUCED Right EIGENVALUES imag part:'
         CALL OUTPUT(WORK(KWI),1,jkaede,1,1,jkaede,1,1,LUPRI)
      END IF
C
      IF (IPRLEE .GE. 10) THEN
         WRITE (LUPRI,'(/A)') ' REDUCED right EIGENVECTORS :'
         CALL OUTPUT(WORK(KEIVER),1,jkaede,1,
     &                            jkaede,jkaede,jkaede,1,LUPRI)
      END IF
C
      IF (ICPLX .GT. 0)
     *   WRITE(LUPRI,*) '**WARNING CCRED: COMPLEX EIGENVALUES.'
C
      IF (IPRLEE .GE. 10) THEN
         WRITE (LUPRI,'(/A)') 'reordered REDUCED left EIGENVECTORS :'
        CALL OUTPUT(WORK(KEIVEL),1,jkaede,1,
     &                           jkaede,jkaede,jkaede,1,LUPRI)
      END IF

      IF (IPRLEE .GE. 1) THEN
            WRITE (LUPRI,'(//A)')
     *      ' Final Eigenvectors in Lanczos basis for :'
         IOFF = 0
         DO I = 1,jkaede
!            WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
!     *      ' - EIGENVALUE NO.',I, ' - EIGENVALUE',WORK(KEIVAL+I-1),
!     *      ' - EIGENVECTOR :'
            WRITE (LUPRI,'(A,I4,A,F13.7,A,F13.7)')
     *   ' - EIGENVALUE NO.',I, 
     &   ' - EIGENVALUE (au) ',WORK(KEIVAL+I-1), 
     &   ' - Im part ', WORK(KWI+I-1)
!            write(lupri,*)' Right Eigenvector '
!            WRITE (LUPRI,'(5F15.8)') (WORK(KEIVER+IOFF+J-1),J=1,jkaede)
!            write(lupri,*)' Left Eigenvector '
!            WRITE (LUPRI,'(5F15.8)') (WORK(KEIVEL+IOFF+J-1),J=1,jkaede)
!            write(lupri,*)'                  '
            IOFF = IOFF + jkaede
         END DO
      END IF
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
        call cc_resid_lanczos(beta(jchain),length,jkaede,
     &       work(kqnew1),work(keiver),work(keival),work(kwi))
!
!------------------------------------------------------------------------
! BUILD R=Q*R_lanczos, if renormalized in full space they give standard 
! Right eigenvectors. We use them to build Ftrans and, if ANALYZE, to check
! the nature of the excitation
!------------------------------------------------------------------------
!
      !length of eigenvectors is now?
      !isymRL = muld2h(isyhop,isyeta)
      !lengthR = nt1am(isymRL)
      !if (.not.CCS) lengthR=lengthR + nt2am(isymRL)
      !pseudoeigenvectors have symmetry of operators. Lanczos eigenvectors
      !are totalsymmetric (or isyeta*isyhop)

      KRfull = kend6
      kend6  = KRfull + length*jkaede
      lwrk6  = lwork - kend6
      IF (LWRK6 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need:', KEND6, 'Available:', LWORK
         CALL QUIT('Insufficient work space in CC_LANCZOS R=Q*R_l')
      END IF

      CALL DZERO(WORK(KRfull),length*JKAEDE)
      WRITE(LUPRI,*) 'Call CC_build_Rfull to build R vectors'
      call CC_build_Rfull(LUQMAT,FQMAT,
     &                    Jkaede,length,WORK(KRfull),work(keiver),
     &                    work(keival),work(kwi),
     &                    WORK(KEND6),LWRK6)
      WRITE(LUPRI,*) 'Out of CC_build_Rfull'
      call flshfo(lupri)

C
C------------------------------------------------------------------
C     Perform F*Q matrix transformations and save on separate files:
C-------------------------------------------------------------------
C
       IOPT = 3
       WRITE(LUPRI,*) 'Calling CC_FMATRIX to perform FQ tranforms'
       call flshfo(lupri)
       CALL CC_FMATRIX(IFQTRAN, NQLLBL, 'L0 ', 'QL ', IOPT, 'FQ ',
     &                 IDUMMY, DUMMY, 0, WORK(KEND6), LWRK6)
       WRITE(LUPRI,*) 'out of CC_FMATRIX'
       call flshfo(lupri)
C
C------------------------------------------------------------------
C     Generate trans^F
C-------------------------------------------------------------------
C
       kftrans = kend6
       kend61  = kftrans + jkaede*jkaede
       lwrk61  = lwork - kend61
       IF (LWRK61.LT. 0) THEN
         WRITE(LUPRI,*) 'Need:', KEND61, 'Available:', LWORK
         CALL QUIT('Insufficient work space in LANCZOS Ftrans')
       END IF
       
       CALL DZERO(work(kFtrans),jkaede*jkaede)
       WRITE(LUPRI,*) 'Calling Lanczos_Ftrans to make transF (calF)'
       call flshfo(lupri)
       CALL CC_Lanczos_Ftrans(IOPTRW,JKAEDE,LENGTH,IFQTRAN,WORK(KRfull),
     &               WORK(KEIVER),work(kFtrans),WORK(KEND61),LWRK61,
     &               WORK(KEIVAL),WORK(KWI),isyhop)

      if (locdbg) then
         !CALL AROUND('--- The trans^F matrix ---')
         !CALL OUTPUT(Ftrans,1,jkaede,1,jkaede,jkaede,jkaede,1,LUPRI)
      end if
C
C-------------------------------------------------------------
C     Compute the response functions at requires freqs and for
C     required gammas from input
C     Results are damped on a separate file!!!!
C-------------------------------------------------------------
C
      if (sum_rules) then
        CALL CC_LANCZOS_SUMRULES(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
     *                  WORK(KEIVER),WORK(KEIVEL),
     &                  work(kFTRANS),ETANORM,CSINORM)
      else
        CALL CC_LANCZOS_LR(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
     *                  WORK(KEIVER),WORK(KEIVEL),
     &                  work(kFTRANS),ETANORM,CSINORM)
        !compute oscillator strengths as residue of Lanczos LR
        CALL CC_LANCZOS_OSCSTR(LABEL,JKAEDE,WORK(KEIVAL),WORK(KWI),
     *                 WORK(KEIVER),WORK(KEIVEL),
     &                 work(kFTRANS),ETANORM,CSINORM)
      end if
C
C--------------------------------------------------------------
C ANALYZE VECTORS AND DUMP ON FILE IF REQUIRED
C--------------------------------------------------------------
C
      !For comparison with Standard code we need to renormalize
      if ((locdbg).or.ABSANALYZE) then

!---------------------------------------------------------------------
! Initialization of stuff needed to be able to dump on RE and LE files
!---------------------------------------------------------------------
        CALL IZERO(NCCEXCI,3*8)
!
!
        NCCEXCI(ISYHOP,1) = JKAEDE
        !FIXME: MSYM SHOULD COME FROM OUTSIDE !msym=8 or 1?????
        msym=NSYM
        write(lupri,*)'MSYM =', MSYM, ' inside Lanczos (NSYM)'
        WRITE (LUPRI,*) 'NCCEXCI for singlet in LANCZOS:',
     &                  (NCCEXCI(ISYM,1),ISYM=1,msym)
C     ----------------------------
C     set up symmetry information: DO I REALLY NEED TO DO THIS????
C     ----------------------------
        NEXCI  = 0
        NTRIP  = 0
        DO ISYM = 1,MSYM
          ISYOFE(ISYM) = NEXCI
          ITROFE(ISYM) = ISYOFE(ISYM) + NCCEXCI(ISYM,1)
          NEXCI        = ITROFE(ISYM) + NCCEXCI(ISYM,3)
          NTRIP        = NTRIP        + NCCEXCI(ISYM,3)
          DO IEX = ISYOFE(ISYM)+1, NEXCI
             ISYEXC(IEX) = ISYM
          END DO
          DO IEX = ISYOFE(ISYM)+1, ITROFE(ISYM)
             IMULTE(IEX) = 1
          END DO
          DO IEX = ITROFE(ISYM)+1, NEXCI
             IMULTE(IEX) = 3
          END DO
        END DO
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

        iopt=1
        call cc_normalize(length,jkaede,work(krfull),work(kend6),
     &                    WORK(KEIVAL),WORK(KWI),isyhop,iopt)

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
           do i=1,nccexci(isyhop,1)
              istart = i
              if (WORK(KEIVAL+I-1).ge.EIG_RANGE(1)) exit
           end do
           write(lupri,*)'Index of FIRST FREQ >= EIG_RANGE(1)',istart
           !determine how many in chosen interval
           nrexc = 0
           do i=1,nccexci(isyhop,1)
              IF ((WORK(KEIVAL+I-1).ge.EIG_RANGE(1)).and.
     &           (WORK(KEIVAL+I-1).le.EIG_RANGE(2))) then
                 nrexc = nrexc + 1
              end if
           end do
           write(lupri,*)'NR OF FREQS IN SELECTED EIG_RANGE',NrExc
           !check residual dump given pseudo RE vectors on file
           koffR=KRFULL+length*(istart-1)
           write(lupri,*)'R vecs RFULL', KRFULL, ' koffR = ',koffR
           write(lupri,*)'Start RE vector is nr=',istart
           write(lupri,*)'Checking selected residuals'
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
              koffre = koffr + (i-1)*length
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
              CALL CC_WRRSP('RE ',I,ISYHOP,IOPTRW,MODEL,WORK(KEND6),
     *                   WORK(KOFFRE),
     *                   WORK(KOFFRE+NT1AM(ISYHOP)),
     *                   WORK(KEND6),LWRK6)
!!              CALL CC_WRRSP('RE ',I,ISYHOP,IOPTRW,MODEL,WORK(KEND6),
!!     *                   WORK(KoffR+IOFF),
!!     *                   WORK(KoffR+IOFF+NT1AM(ISYHOP)),
!!     *                   WORK(KEND6),LWRK6)
!              !write(lupri,*)'Below RE_i=', i, ' of address= ',koffR+ioff
              CALL CC_PRAM(WORK(KoffRe),PT1,ISYHOP,.FALSE.)
              !IOFF = IOFF + length
           end do
!
           !determine if we have space enough for the corresponding LE vectors

           KLfull = kend6
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
           koffR = kRfull + length*(istart-1) 
           WRITE(LUPRI,*) 'cc_normalize: binormalize LE vectors'
           call cc_normalize(length,nrexc,work(koffR),work(klfull),
     &                    WORK(KEIVAL+istart-1),WORK(KWI+istart-1),
     &                    isyhop,iopt)
           !ioff = 0
           !koffl = kLfull
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
!!
       end if

C-------------------------------------------------------------
C     finished
C-------------------------------------------------------------

      !Close files not required any more (A lh & rh transform)
      !Could be closed before this
      CALL CC_FILCL(FRHO1,LUAQ11,FRHO2,LUAQ12,CDUM,LUDUM,
     *              FRHS1,LUQ11,FRHS2,LUQ12,CDUM,LUDUM,
     *              CDUM,LUDUM,CDUM,LUDUM)
      CALL CC_FILCL(FP1A1,LUP1A1,FP1A2,LUP1A2,CDUM,LUDUM,
     *              FP11,LUP11,FP12,LUP12,CDUM,LUDUM,
     *              CDUM,LUDUM,CDUM,LUDUM)

      end if !isyma=isymb

      CALL QEXIT('CC_LANCZOS_DRV')
      RETURN
      END
