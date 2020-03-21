#define SRDFT_DEBUG -1
!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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


C*****************************************************************************
      SUBROUTINE SRDFT1JT(ND_SIM,EXCMAT,DMAT,DOLND,CORX,CORY,CORZ,
     &                  WEIGHT,NBUF,GSO,NCNT,DMATGAOSR,DOGGA,DOERG,
     &                  DO_MOLGRAD,DOATR,TRIPLET,do_metagga,RESULTS,
     &                  WORK,LWORK,IPRINT_in)
C*****************************************************************************
C     Calculate DFT energy and potential contributions for
C     SR-DFT hybrid methods
C
C     Based on the routine SRDFT1 of J. K. Pedersen
C
C    Input: ND_SIM: number of input density matrices and output potential matrices
C                  DOERG   true: energy and Vxc potential matrix
C                  DO_MOLGRAD   true: molecular gradient
C                  DOATR   true: electronic Hessian, non-linear Exc second derivative terms
C                  TRIPLET true: triplet operators
C                  DMAT(1) = charge density matrix
C                  DMAT(2) = 1-index transformed charge density matrix (DOATR)
C                       or = 1-index transformed spin density matrix (TRIPLET and DOATR)
C                  DMAT(3) = spin density matrix (DFT_SPINDNS)
C                  DMAT(4) = 1-index transformed spin density matrix (DFT_SPINDNS and DOATR)
C            (if SRHYBR then different ordering of DMAT matrices)
C            (Note that the order of DMATs also fits for gradient calculations, DOATR false:
C             for energy and gradien calculations DMAT(2) will be the DVAO matrix needed
C             for Coulomb and exchange Fock matrices and not needed for SRDFT.)
C
C                  EXCMAT(1) = charge density potential matrix
C                  EXCMAT(2) = spin density potential matrix (DFT_SPINDNS)
C
C     Created : 15-02-05, J. Toulouse
C     Modified: 17-08-09, J. Toulouse, add spin density
C     Modified: Feb. 2016 E. Hedegaard add spin gradient
C*****************************************************************************
      implicit none
#include "priunit.h"
#include "mxcent.h"
#include "dummy.h"
! infpar.h : MYNUM, NODTOT
#include "maxorb.h"
#include "gnrinf.h"
#include "infpar.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "inftap.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"
#include "chrnos.h"
      ! Need intention on many of the inputs for this subroutine!
      real*8, parameter :: D0 = 0.0D0
      real*8, parameter :: DP5 = 0.5D0
      real*8, parameter :: D1 = 1.0D0
      real*8, parameter :: D2 = 2.0D0
      real*8, parameter :: THRELCTOT = 1.0D-7
      real*8, parameter :: THRESRHO  = 1.0D-8

      logical, intent(in) :: DOLND, DOGGA, DOERG, DO_MOLGRAD, DOATR,
     &                       TRIPLET, do_metagga
      integer, intent(in) :: ND_SIM, NBUF, IPRINT_in

      logical :: FROMVX
      real*8 :: EXCMAT(NBAST,NBAST,*), DMAT(N2BASX,*),
     &          CORX(NBUF), CORY(NBUF), CORZ(NBUF), WEIGHT(NBUF),
     &          GSO(NBAST*NTYPSO), NCNT(NBAST),
     &          DMATGAOSR(NBAST,2),
     &          COR(3), RESULTS(11), WORK(LWORK)

      character(len=17) :: QUADNAME

      integer, save :: NTOT_DFTGRID_save = -1

      integer :: LWORK, i, IPNT, IPNT_END,
     &           IPNT_ST, IPRINT, IPRINT_x, JDXSAO, JDXTAO, LDUMMY,
     &           LUQUAD, my_first_point, my_last_point, N_RS_GT_RT,
     &           NORDER, NPNT1, NPNTS, NPNTS_PER_NODE, NPOINT, NSKIP1,
     &           KSO1_gx, KSO1_gy, KSO1_gz
      real*8 :: DELCTRN, EcInt, ELCTRN, ELCTRSKIP, ExInt,
     &          SPINDEN, THRELCTRN, THRINT, WGHT, WGHTC, WGHTX,
     &          TAU_X, TAU_Y, TAU_Z

      real*8 :: RHO(4), RHG(3,4), GRD2(6), tau(4), laplace(4)
C     ... RHO(1) = rhoc = rho_a + rho_b
C     ... RHO(2) = rhos = rho_a - rho_b
C     ... RHO(3) = rhoa = rho_a
C     ... RHO(4) = rhob = rho_b
C     ... RHG(1:3,1) = RHGC(1:3) = grad ( rho_a + rho_b ) = grad ( rho_c )
C     ... RHG(1:3,2) = RHGS(1:3) = grad ( rho_a - rho_b ) = grad ( rho_s )
C     ... RHG(1:3,3) = RHGA(1:3) = grad ( rho_a )
C     ... RHG(1:3,4) = RHGB(1:3) = grad ( rho_b )
C     ... GRD2(1) = grdcc = RHG(:,1) . RHG(:,1)
C     ... GRD2(2) = grdss = RHG(:,1) . RHG(:,2)
C     ... GRD2(3) = grdcs = RHG(:,2) . RHG(:,2)
C     ... GRD2(4) = grdaa = RHG(:,3) . RHG(:,3)
C     ... GRD2(5) = grdbb = RHG(:,4) . RHG(:,4)
C     ... GRD2(6) = grdab = RHG(:,3) . RHG(:,4)

!     these are not initialized to zero!
      real*8 :: Ex,  d1Ex(9),  d2Ex(45)
      ! d1Ex(1:5) =   d1Exdrhoc, d1Exdrhos, d1Exdgrdcc, d1Exdgrdss, d1Exdgrdcs

      real*8 :: Ec,  d1Ec(9),  d2Ec(45)
      ! d1Ec(1:5) =   d1Ecdrhoc, d1Ecdrhos, d1Ecdgrdcc, d1Ecdgrdss, d1Ecdgrdcs

      real*8 :: Exc, d1Exc(9), d2Exc(45)
      ! d1Exc(1:5) =  d1Excdrhoc, d1Excdrhos, d1Excdgrdcc, d1Excdgrdss, d1Excdgrdcs

#if SRDFT_DEBUG > 10
      integer :: j, k
      real*8 :: print_thr
      print_thr = 0.0d0
#endif

      IPRINT = MAX(IPRINT_in,SRDFT_DEBUG)
!     Debug
      IF (IPRINT .GE. 5) THEN
      WRITE(LUPRI,*) 'SRDFT1JT: entering, ND_SIM = ',ND_SIM
      WRITE(LUPRI,*) 'SRDFT1JT: MYNUM, NODTOT = ',MYNUM,NODTOT
      WRITE(LUPRI,*) 'SRDFT1JT: LUPRI,IPRINT  =',LUPRI, IPRINT
      WRITE(LUPRI,*) 'SRDFT1JT: CHIVAL =',CHIVAL
      WRITE(LUPRI,*) 'SRDFT1JT: HFXFAC =',HFXFAC
      WRITE(LUPRI,*) 'SRDFT1JT: DOERG,DOATR,DO_MOLGRAD,DOGGA,TRIPLET  ='
     &   ,DOERG,DOATR,DO_MOLGRAD,DOGGA,TRIPLET
      WRITE(LUPRI,*) 'SRDFT1JT: ERFEXP(0:2)=',ERFEXP(0:2)
      WRITE(LUPRI,*) 'SRDFT1JT: DOSRX_PBERI=',DOSRX_PBERI
      WRITE(LUPRI,*) 'SRDFT1JT: DFT_SPINDNS=',DFT_SPINDNS
      flush(lupri)
      IF (IPRINT .GE. 12) THEN
      WRITE(LUPRI,*) 'D_TOT_AO matrix:'
      call output(DMAT(1,1),1,nbast,1,nbast,nbast,nbast,-1,lupri)
      IF (DOATR) THEN
         WRITE(LUPRI,*) 'DX_TOT_AO matrix for DOATR true:'
         call output(DMAT(1,2),1,nbast,1,nbast,nbast,nbast,-1,lupri)
      END IF
      IF (DFT_SPINDNS) THEN
         WRITE(LUPRI,*) 'D_SPIN_AO matrix:'
         call output(DMAT(1,3),1,nbast,1,nbast,nbast,nbast,-1,lupri)
         IF (DOATR) THEN
            WRITE(LUPRI,*) 'DX_SPIN_AO matrix for DOATR:'
            call output(DMAT(1,4),1,nbast,1,nbast,nbast,nbast,-1,lupri)
         END IF
      END IF
      END IF ! iprint .ge. 12
      flush(lupri)
      END IF ! iprint .ge. 5

!     Is the requested calculation implemented ?

      IF (ND_SIM .NE. 1) THEN
         WRITE(LUPRI,*) 'SRDFT1: illegal ND_SIM ',ND_SIM
         CALL QUIT('SRDFT1jt: illegal ND_SIM value')
      END IF

      IF (DO_MOLGRAD) THEN
         CALL QUIT(
     &  'SRDFT1JT: molecular gradient (DO_MOLGRAD) not implemented yet')
      END IF

      IF (DOLND) THEN
         CALL QUIT(
     &   'SRDFT1JT: London orbitals (DOLND) not implemented yet')
      END IF

      IF (DFT_SPINDNS .AND. TRIPLET) THEN
         CALL QUIT('SRDFT1JT: triplet only for closed shell')
      END IF
      IF (DOERG .AND. TRIPLET) THEN
         CALL QUIT('SRDFT1JT: triplet only for response')
      END IF

C     Initializations -------------------------------------------------------

      RHG(:,:)   = 0.0d0 ! so initialized for LDA
      GRD2(:)    = 0.0d0
      tau(:)     = 0.0d0
      laplace(:) = 0.0d0

      ! edh the arrays are also zeroed in SRDFTEXCJT
      d1Ex(:)  = 0.0d0
      d1Ec(:)  = 0.0d0
      !---------------

!     (Integrated) exchange and correlation energies
      ExInt  = D0
      EcInt  = D0

C     DFT grid (NTOT_DFTGRID is the total number of grid points)
      IF (.NOT.DFTGRID_DONE_OLD .OR. NTOT_DFTGRID_save .LE. 0) THEN
         CALL MAKE_DFTGRID(WORK,LWORK,NTOT_DFTGRID,1,.FALSE.)
         CALL CONDFT
         NTOT_DFTGRID_save = NTOT_DFTGRID
         DFTGRID_DONE_OLD = .TRUE.
      ELSE
         NTOT_DFTGRID = NTOT_DFTGRID_save
         ! because wrong NTOT_DFTGRID has been transferred from master in SRDFT_PAR_NODE
      END IF

C     Electron number
      ELCTRN = D0

C     Integrated spin density
      SPINDEN = D0

C     Thresholds for screening
C     THRELCTRN makes sure that error in no. of electrons
C     from grid integration less than DFTHR
      !real*8, parameter :: THRELCTOT = 1.0D-7
      !real*8, parameter :: THRESRHO  = 1.0D-8
      DFTHRI    = D0 ! real*8 from dftcom.h
      DFTHRL    = D0 ! real*8 from dftcom.h
      THRELCTRN = D0

C     Skipped grid points and nb of electrons
C     NSKIP1 is not used anywhere? -E. R. K.
      NSKIP1 = 0
      ELCTRSKIP = D0
C
C     Counter for how many times RHO(2) > RHO(1)
      N_RS_GT_RT = 0

C     For quadrature integration
C
C     Make quadname: Can take 9999 procs
C
      IF (MYNUM .EQ. 0) THEN
         QUADNAME = 'DALTON.QUAD'
      ELSE
         QUADNAME = 'DALTON.QUAD.n'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
      END IF
C
      LUQUAD = -1
      CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)
C
      IF (MYNUM .EQ. 0) THEN
         NPNTS_PER_NODE = NTOT_DFTGRID
         my_first_point = 1
         my_last_point  = NTOT_DFTGRID
      ELSE
         NPNTS_PER_NODE = (NTOT_DFTGRID-1)/NODTOT + 1
         my_first_point = NPNTS_PER_NODE*(MYNUM-1) + 1
         my_last_point  = MIN(NPNTS_PER_NODE*MYNUM, NTOT_DFTGRID)
      END IF
      if (iprint .ge. 5) then
         write(lupri,*)
     &      'MYNUM,NPNTS_PER_NODE,my_first_point,my_last_point',
     &       MYNUM,NPNTS_PER_NODE,my_first_point,my_last_point
         flush(lupri)
      end if
      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
      if (iprint .ge. 5) then
         write(lupri,*) 'MYNUM, NPOINT from LUQUAD',MYNUM,NPOINT
      end if

      IF (NPOINT.GT.0) THEN
         NPNT1 = NPNTS + 1       ! first integration point from this record
         NPNTS = NPNTS + NPOINT  ! last  integration point from this record
         IPNT_ST  = max(NPNT1,my_first_point) + 1 - NPNT1
         IPNT_END = min(NPNTS,my_last_point)  + 1 - NPNT1
         if (iprint .ge. 5) then
            write(lupri,*)
     &      'MYNUM,IPNT_ST,IPNT_END', MYNUM,IPNT_ST,IPNT_END
            flush(lupri)
         end if
         IF (IPNT_ST .GT. IPNT_END) THEN
            if (iprint .ge. 5) write(lupri,*) 'Skipping this record'
            READ (LUQUAD) ! skip this record
            GO TO 200
         END IF
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)

C   Loop over grid points -----------------------------------------------

      DO 300 IPNT = IPNT_ST, IPNT_END

C     Quadrature weight
      WGHT  = WEIGHT(IPNT)
      WGHTX = (1.0D0-HFXFAC)*WGHT
      WGHTC = WGHT

C     Get AOs
C     DFTHRI is always zero, thus THRINT is always zero
C       Is this a bug or feature?
      THRINT = DFTHRI/WGHT
      COR(1) = CORX(IPNT)
      COR(2) = CORY(IPNT)
      COR(3) = CORZ(IPNT)

!     print for every 1000th grid point if iprint .gt. 10
      IPRINT_x = 0
      IF (IPRINT .GT. 10) THEN
      IF( MOD(IPNT,1000) .EQ. 1) THEN
!     IF( IPNT .gt. 18000 .and. IPNT .lt. 19000
!    &   .and. mod(IPNT,100) .eq. 1) then
         IPRINT_x = iprint
      END IF
      END IF

      CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,NBAST,DOLND,DOGGA,
     &            THRINT,IPRINT_x)

C     Density
      CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHO(1),
     &       DMATGAOSR(1,1),THRINT,IPRINT_x)

C     Contribution to electron number (N_alpha + N_beta)
      DELCTRN = WGHT*RHO(1)

      IF (iprint_x .gt. 10) then
C        Print grid information
         WRITE (LUPRI,'(A,I8,3F15.8,1P,D15.8)') 'Grid point',
     &      IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
         WRITE (LUPRI,*) 'RHOC, DELCTRN, ELCTRN',
     &      RHO(1),DELCTRN,ELCTRN
      END IF

C     Screening
      IF (ABS(DELCTRN).LE.THRELCTRN .OR.
     &          RHO(1).LE.THRESRHO) THEN

         IF (iprint_x .gt. 10) then
            write(lupri,*) 'Skipping this point'
         END IF
         NSKIP1 = NSKIP1 + 1
         ELCTRSKIP = ELCTRSKIP + DELCTRN
         THRELCTRN = (THRELCTOT-ELCTRSKIP)/100
         GO TO 300
C        ... skip this point, go to next grid point
      END IF

C     Spin density rho_s ( in RHO(2) )
      IF (DFT_SPINDNS) THEN
C        DMAT(1,3) is the spin density matrix
         CALL GETRHO_OLD(DMAT(1,3),GSO(KSO0),RHO(2),
     &                   DMATGAOSR(1,2),THRINT,IPRINT)

         IF (DFT_LOCALSPIN) THEN
C        In our local spin model we can get unphysical densities
C        because we allow DS > DV in some regions.
            IF (RHO(2) .GT. RHO(1)) THEN
                RHO(2) = RHO(1)
                N_RS_GT_RT = N_RS_GT_RT + 1
            END IF
         END IF
C        Contribution to spin density (N_alpha - N_beta)
         SPINDEN = SPINDEN + WGHT*RHO(2)
      ELSE
         RHO(2) = 0.0D0
      END IF
      RHO(3) = 0.5d0 * ( RHO(1) + RHO(2) ) ! rho_a
      RHO(4) = 0.5d0 * ( RHO(1) - RHO(2) ) ! rho_b
      ! For small values of rho to avoid NaN
      do i = 1,4
         if (abs(rho(i)) .lt. 1.d-30) rho(i) = 1.d-30
      end do
      if (rho(1) .gt. 1.d-10) then
         if ((rho(1) - abs(rho(2))) .le. 0.0d0) then
            call quit('rho_c is less or equal to |rho_s|')
         end if
         if ((rho(1) - abs(rho(2)) - 1.d-4*rho(1)) .le. 0.0d0 ) then
            call quit('rho_c and rho_s are too close for finite diff.')
         end if
      end if

C     Gradients of densities
      IF (DOGGA) THEN
C       GSO(KSO1)*DMATGAOSR
        CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMATGAOSR(1,1),
     &     1,D0,RHG(1,1),1)
C
        IF (DFT_SPINDNS) THEN
           CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMATGAOSR(1,2),
     &        1,D0,RHG(1,2),1)
        ELSE
           RHG(1:3,2) = 0.0D0 ! RHGS is zero
        ENDIF
        DO I = 1,3
           RHG(I,3) = 0.5d0 * ( RHG(I,1) + RHG(I,2) ) ! RHGA
           RHG(I,4) = 0.5d0 * ( RHG(I,1) - RHG(I,2) ) ! RHGB
        END DO
!       Calculate the six variables related to four forms of grad(rho)
!       (grdcc, grdss, grdcs, grdaa, grdbb, grdab)
        grd2(1) = rhg(1,1)**2 + rhg(2,1)**2 + rhg(3,1)**2  ! grdcc
        grd2(2) = rhg(1,2)**2 + rhg(2,2)**2 + rhg(3,2)**2  ! grdss
        grd2(3) = rhg(1,1)*rhg(1,2) + rhg(2,1)*rhg(2,2)    ! grdcs = grad(rhoc) . grad(rhos)
     &          + rhg(3,1)*rhg(3,2)
        grd2(4) = rhg(1,3)**2 + rhg(2,3)**2 + rhg(3,3)**2  ! grdaa
        grd2(5) = rhg(1,4)**2 + rhg(2,4)**2 + rhg(3,4)**2  ! grdbb
        grd2(6) = rhg(1,3)*rhg(1,4) + rhg(2,3)*rhg(2,4)    ! grdab = grad(rhoa) . grad(rhob)
     &          + rhg(3,3)*rhg(3,4)
      END IF

      if (do_metagga) then
         KSO1_gx = KSO1
         KSO1_gy = KSO1 + NBAST
         KSO1_gz = KSO1 + 2*NBAST
         CALL GETRHO_OLD(DMAT(1,1),GSO(KSO1_gx),TAU_X,
     &       DMATGAOSR(1,1),THRINT,IPRINT_x)
         CALL GETRHO_OLD(DMAT(1,1),GSO(KSO1_gy),TAU_Y,
     &       DMATGAOSR(1,1),THRINT,IPRINT_x)
         CALL GETRHO_OLD(DMAT(1,1),GSO(KSO1_gz),TAU_Z,
     &       DMATGAOSR(1,1),THRINT,IPRINT_x)
         TAU(1) = 0.5D0*(TAU_X + TAU_Y + TAU_Z)
         IF (DFT_SPINDNS) THEN
C           DMAT(1,3) is the spin density matrix
            CALL GETRHO_OLD(DMAT(1,3),GSO(KSO1_gx),TAU_X,
     &          DMATGAOSR(1,1),THRINT,IPRINT_x)
            CALL GETRHO_OLD(DMAT(1,3),GSO(KSO1_gy),TAU_Y,
     &          DMATGAOSR(1,1),THRINT,IPRINT_x)
            CALL GETRHO_OLD(DMAT(1,3),GSO(KSO1_gz),TAU_Z,
     &          DMATGAOSR(1,1),THRINT,IPRINT_x)
            TAU(2) = 0.5D0*(TAU_X + TAU_Y + TAU_Z)
         ELSE
            TAU(2) = 0.0D0
         END IF
         TAU(3) = 0.5d0 * ( TAU(1) + TAU(2) ) ! tau_a
         TAU(4) = 0.5d0 * ( TAU(1) - TAU(2) ) ! tau_b
         ! For small values of tau to avoid NaN
         !do i = 1,4
            ! this is not a stable condition. TPSSc and TPSSx
            ! have reversed tau dependency.
            !if (abs(tau(i)) .lt. 1.d-30) tau(i) = 1.d-30
         !end do
      ENDIF
C
C     Number of electrons
      ELCTRN = ELCTRN + DELCTRN

      IF (DOERG) THEN
C
C      Exchange and correlation energies and derivatives
       NORDER = 1  ! first-order derivatives
C
       CALL SRDFTEXCJT(RHO,GRD2,tau,laplace,CHIVAL,NORDER,ERFEXP,
     &                 .false.,Ex,d1Ex,d2Ex,Ec,d1Ec,d2Ec,VLAMBDA)
C      Summation with quadrature weights
       ExInt     = ExInt + WGHTX*Ex
       EcInt     = EcInt + WGHTC*Ec

C      Sum of exchange and correlation
       d1Exc(1:9) = WGHTX*d1Ex(1:9) + WGHTC*d1Ec(1:9)

       IF (iprint_x .gt. 10) then
          write(lupri,*)'rhoc, rhcg ',rho(1),rhg(1:3,1)
          write(lupri,*)'rhos, rhsg ',rho(2),rhg(1:3,2)
          write(lupri,*)'rhoa, rhag ',rho(3),rhg(1:3,3)
          write(lupri,*)'rhob, rhbg ',rho(4),rhg(1:3,4)
          write(lupri,*)'grd2       ',grd2(1:6)
          write(lupri,*)'tau        ',tau(1:4)
          write(lupri,*)'Ex, ExInt',Ex,ExInt
          write(lupri,*)'Ec, EcInt',Ec,EcInt
          write(lupri,*)'WGHTX, WGHTC  ', WGHTX, WGHTC
          write(lupri,*)'d1Exdrhoc  ...', d1Ex(1:9)
          write(lupri,*)'d1Ecdrhoc  ...', d1Ec(1:9)
          write(lupri,*)'d1Excdrhoc ...', d1Exc(1:9)
       end if

C      Exchange-correlation contribution to Kohn-Sham matrix (charge and spin components)
       FROMVX = .FALSE.
       if (domcsrdft .and. docisrdft) then  !  only true for start-CI  (cf. routine OPST)
       ! Ignore spin in start CI ...  FIXME TODO
           CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG(1,1),
     &                 d1Exc(1)  ,d1Exc(3)    ,DOGGA,FROMVX,DFTHRL)
!                      d1Excdrhoc,d1Excdgrdcc,DOGGA,FROMVX,DFTHRL)
       else if (domcsrdft .and. .not. docisrdft) then
          if (do_metagga) then
             if (dft_spindns) then
               CALL DFTKSM_MGGA(EXCMAT, GSO(KSO0), GSO(KSO1), 0.0d0,
     &                             RHG, d1Exc, DFTHRL, .TRUE.)
             else
               CALL DFTKSM_MGGA(EXCMAT, GSO(KSO0), GSO(KSO1), 0.0d0,
     &                             RHG, d1Exc, DFTHRL, .FALSE.)
             end if
          else if (dogga .and. dft_spindns) then ! special treatment if gga + spin
              CALL DFTKSMGGASPIN(EXCMAT,GSO(KSO0),GSO(KSO1),RHG,
     &             d1Exc,DFTHRL)
!                  d1Exc = (d1Excdrhoc,d1Excdrhos, d1Excdgrdcc,d1Excdgrdss,d1Excdgrdcs)
          else                              ! GGA without spin density or LDA
             CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG(1,1),
     &                d1Exc(1)  ,d1Exc(3)    ,DOGGA,FROMVX,DFTHRL)
!                     d1Excdrhoc,d1Excdgrdcc,DOGGA,FROMVX,DFTHRL)
             if (dft_spindns .and. .not. dogga) then ! i.e. LDA with spin density
                 CALL DFTKSM(EXCMAT(1,1,2),GSO(KSO0),GSO(KSO1),RHG(1,2),
     &                 d1Exc(2)  ,0.0d0,DOGGA,FROMVX,DFTHRL)
!    &                 d1Excdrhos,0.0d0,DOGGA,FROMVX,DFTHRL)
             end if
          end if
       else ! if hfsrdft
          IF (DFT_SPINDNS) CALL QUIT('no HFsrDFT with open shells')
          if (do_metagga) then
             CALL DFTKSM_MGGA(EXCMAT, GSO(KSO0), GSO(KSO1), 0.0d0,
     &                            RHG, d1Exc, DFTHRL, .FALSE.)
          else
             CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG(1,1),
     &                d1Exc(1)  ,d1Exc(3)    ,DOGGA,FROMVX,DFTHRL)
!                     d1Excdrhoc,d1Excdgrdcc,DOGGA,FROMVX,DFTHRL)
          end if
       end if

      ENDIF  !DOERG

C     Hessian transformation for second-order terms
      IF (DOATR) THEN
C        Notation: "DT" for charge density, "DS" for spin density,
C                  "DX*" for 1-index transformed or transition density
C                  "AO" for in AO basis
C        Three cases:
         if (dft_spindns) then  ! open-shell non-singlet reference
         ! four input matrices: DTAO,  DXTAO, DSAO,  DXSAO
         ! two output matrices: EXCMAT(:,:,1) = W2CAO, EXCMAT(:,:,2) = W2SAO
            JDXTAO = 2
            JDXSAO = 4
         else if (triplet) then !  singlet reference, triplet perturbation
         ! two input matrices: DTAO,  DXSAO
         ! one output matrix : EXCMAT(:,:,1) = W2SAO
         ! Note that DMAT(1,2) = DXSAO when TRIPLET = .TRUE.
            JDXTAO = 2
            JDXSAO = 2
         else                   !  singlet reference, singlet perturbation
         ! two input matrices: DTAO,  DXTAO
         ! one output matrix : EXCMAT(:,:,1) = W2CAO
            JDXTAO = 2
            JDXSAO = 2
         end if

#if SRDFT_DEBUG > 41
       print_thr = 1000.0d0
       do k = 1,2
         do i = 1, nbast
           do j = 1, i - 1
             if (abs(excmat(j,i,k)) > print_thr) then
                write(lupri,*)'edh 0 k, j, i', k, i, j, excmat(j,i,k)
             endif
           end do
         end do
       end do
#endif

! SRDFTLTRJT: SR short-range DFT density functional theory LTR linear transformation JT Julien Toulouse
        CALL SRDFTLTRJT(JWOPSY,DMAT(1,JDXTAO),DMAT(1,JDXSAO),EXCMAT,
     &                 WGHTX,WGHTC,GSO(KSO0),GSO(KSO1),CHIVAL,VLAMBDA,
     &                 ERFEXP,RHO,RHG,GRD2,tau,laplace,DMATGAOSR,
     &                 DOGGA,do_metagga,TRIPLET,IPRINT_X)

      ENDIF ! DOATR

#if SRDFT_DEBUG > 41
       print_thr = 50000.0d0
       do k = 1,2
         do i = 1, nbast
           do j = 1, i - 1
             if (abs(excmat(j,i,k)) > print_thr) then
                write(lupri,*)'edh 1 k, j, i', k, i, j, excmat(j,i,k)
             endif
           end do
         end do
       end do
#endif

  300 CONTINUE
C End of loop over grid points ------------------------------------------------
C
         GO TO 200
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
      END IF  !NPOINT.GT.0

C Termination --------------------------------------------------------------

#if SRDFT_DEBUG > 21
       print_thr = 10.0d0
       write(lupri,*) 'print_thr', print_thr
       do k = 1,2
         do i = 1, nbast
           do j = 1, i - 1
             if (abs(excmat(j,i,k)) > print_thr) then
                write(lupri,*)'edh 2 k, j, i', k, i, j, excmat(j,i,k)
             endif
           end do
         end do
       end do
#endif
      CALL GPCLOSE(LUQUAD,'KEEP')

      RESULTS(1) = ELCTRN
      RESULTS(2) = ELCTRSKIP
      RESULTS(3) = SPINDEN
      RESULTS(4) = ExInt
      RESULTS(5) = EcInt
      RESULTS(6) = N_RS_GT_RT
      ! LSRHYBR
      RESULTS(7:11) = 0.0D0

      RETURN
      END


C*****************************************************************************
      SUBROUTINE SRDFTEXCJT(RHO,GRD2,tau,laplace,CHIVAL,NORDER,ERFEXP,
     &           triplet,Ex,d1Ex,d2Ex,Ec,d1Ec,d2Ec,VLAMBDA)
c*****************************************************************************
C    Driver for short-range DFT exchange and correlation
C    energies and derivatives
C
C    Based on the routine SRDFTEXC of J. K. Pedersen
C
C    Created : 15-02-05, J. Toulouse
C    Modified: 17-08-09, J. Toulouse, add spin density
C    Modified: 15-08-16, E. Hedegaard, add spin density gradient
C*****************************************************************************
      implicit none
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "dftinf.h"
#include "dftcom.h"

      real*8, intent(in) :: rho(4), grd2(6), tau(4), laplace(4),
     &                      CHIVAL, VLAMBDA
!     grd2(1:6) = (grdcc, grdss, grdcs, grdaa, grdbb, grdab)
      integer, intent(in) :: NORDER
      logical, intent(in) :: triplet

      real*8, intent(out) :: Ex, d1Ex(9), d2Ex(45)
      real*8, intent(out) :: Ec, d1Ec(9), d2Ec(45)

!     Maybe ERFEXP should be set correctly before calling the 
!     subroutine, such that no input that is not an output is
!     modified
!       -- E.R. Kjellgren
      logical :: ERFEXP(0:2)
     
!     Integer to control request for special cases of derivatives of the
!     functionals.
!     special_case = 1 : Requests a Hessian for triplet response based
!                        on a singlet reference.
!                  = 2 : No more special cases yet.
      integer :: special_case
      real*8 :: mu
      integer, save :: icall = 0

C     exchange functionals
      external ESRX_PBEHSEERF, ESRX_PBETCSERF, ESRX_PBERIERF
      external ESRX_PBEGWSERF

C     correlation functionals
      external ESRC_PBETCSJERF, ESRC_PBETCSERF, ESRC_PBERIERF
      external ESRC_LYPRIERF
      external ESRC_PBEWIERF, ESRC_PBEGWSERF, ESRC_PBELOERF

C     spin functionals
      external ESRX_LDAERFSPIN
      external ESRC_LDAERFSPIN
      external ESRC_MULOC_GGA, ESRC_MULOD_GGA, ESRC_MULOE_GGA

C     sr mdc func
      external ESRC_MD_LDAERF

      if (triplet) then
         special_case = 1
      else
         special_case = 0
      end if

c  Exchange  ------------------------------------------------------------

      Ex      = 0.0d0
      d1Ex(:) = 0.0d0
      d2Ex(:) = 0.0d0

      icall = icall + 1
      ! Doesnt seem like this do anything? Should be removed?
      IF (ERFEXP(0)) THEN
         write (lupri,*)
     &   'WARNING srdftexcjt: erfexp true in call no. ',icall
         erfexp(1:2) = .false.
!        CALL QUIT('erfgau not implemented in SRDFTEXCJT')
      END IF


!     LDA exchange functional
      IF (DOSRX_LDA) THEN
         call LDAX_ERF(rho, norder, chival, Ex, d1Ex, d2Ex)

!     short-range exchange PBE of HSE
      ELSE IF (DOSRX_PBEHSE) THEN
        call DESRX_SPINGGA(ESRX_PBEHSEERF,rho,grd2,CHIVAL,
     &                     norder,Ex,d1Ex,d2Ex)

!     short-range exchange PBE of TCS
      ELSE IF (DOSRX_PBETCS) THEN
        call DESRX_SPINGGA(ESRX_PBETCSERF,rho,grd2,CHIVAL,
     &                     norder,Ex,d1Ex,d2Ex)

!     short-range exchange PBE of GWS ! JT 11-08-09
      ELSE IF (DOSRX_PBEGWS) THEN
         call PBEX_GWS_ERF(rho, grd2, norder, chival,
     &                      Ex, d1Ex, d2Ex)

!     short-range exchange RI
      ELSE IF (DOSRX_PBERI) THEN ! added Manu 01-02-2006
        call DESRX_SPINGGA(ESRX_PBERIERF,rho,grd2,CHIVAL,
     &                     norder,Ex,d1Ex,d2Ex)

!     TPSS exchange functional
      ELSE IF (DOSRX_TPSS) THEN
         call TPSSX_GWS_ERF(rho, grd2, tau, norder, chival,
     &                       Ex, d1Ex, d2Ex)

!     wPBE exchange functional
      ELSE IF (DOSRX_wPBE) THEN
         call wPBEX(rho, grd2, norder, chival, Ex, d1Ex, d2Ex)

C linear complement spin-dependent exchange LDA functional, Kamal SHARKAS 09-05-11
      ELSE IF (DOLAX_LDAS) THEN
         call ELX_LDALAMSPIN (rho(1),rho(2),VLAMBDA,Ex)
         call DELSPIN(ESRX_LDAERFSPIN,rho(1),rho(2),VLAMBDA,DFT_SPINDNS,
     &                NORDER,d1Ex,d2Ex)

C linear complement exchange PBE functional, Kamal SHARKAS 19-05-11
      ELSE IF (DOLAX_PBEGWS) THEN
         call ELAX_PBEGWS(rho,grd2,VLAMBDA,Ex)
         call DELA(ESRX_PBEGWSERF,rho,grd2,VLAMBDA,NORDER,d1Ex,d2Ex)

      ENDIF


C  Correlation  ----------------------------------------------------------

      Ec      = 0.0d0
      d1Ec(:) = 0.0d0
      d2Ec(:) = 0.0d0

!     short-range spin-dependent correlation LDA ! JT 18-08-09
      IF (DOSRC_LDA_PW92) THEN
         call PW92C_ERF(rho, norder, special_case, chival, 
     &                  Ec, d1Ec, d2Ec)

!     short-range correlation PBE of TCS
      ELSE IF (DOSRC_PBETCS) THEN
        call DESR(ESRC_PBETCSERF,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

!     short-range correlation PBE of TCSJ
      ELSE IF (DOSRC_PBETCSJ) THEN
        call DESR(ESRC_PBETCSJERF,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

!     Modified short-range correlation PBE of GWS
      ELSE IF (DOSRC_PBEGWS) THEN
        call PBEC_GWS_ERF(rho, grd2, norder, special_case, chival,
     &                    Ec, d1Ec, d2Ec)

!     rational interpolation
      ELSE IF (DOSRC_PBERI) THEN
        call DESR(ESRC_PBERIERF,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

!     rational interpolation LYP
      ELSE IF (DOSRC_LYPRI) THEN
        call DESR(ESRC_LYPRIERF,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

!     PBE local energy formula
      ELSE IF (DOSRC_PBELO) THEN
        call DESR(ESRC_PBELOERF,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

!     weighted interpolation
      ELSE IF (DOSRC_PBEWI) THEN
        call DESR(ESRC_PBEWIERF,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

!     TPSS correlation spin-unpolarized functional
      ELSE IF (DOSRC_TPSS) THEN
         call TPSSC_GWS_ERF(rho, grd2, tau, norder, special_case,
     &                      chival, Ec, d1Ec, d2Ec)

!     PBE correlation spin-unpolarized functional
      ELSE IF (DOC_PBE_nomu) THEN
         call PBEC_nomu(rho, grd2, norder, special_case, Ec, d1Ec, d2Ec)

!     VWN5 correlation spin-unpolarized functional, with no range-separation
      ELSE IF (DOC_VWN5_nomu) THEN
         call VWN5C_nomu(rho, norder, special_case, Ec, d1Ec, d2Ec)

!     VWN5 correlation spin-unpolarized functional
      ELSE IF(DOSRC_VWN5 .OR. DOSRC_LDA) THEN
         call VWN5C_ERF(rho, norder, special_case, chival, 
     &                  Ec, d1Ec, d2Ec)

C linear Complement spin-dependent (NON-Scaled) correlation LDA functional
      ELSE IF (DOLANSC_LDAS) THEN
         call ELC_LDALAMSPIN (rho(1),rho(2),VLAMBDA,DFT_SPINDNS,Ec)
         call DELSPINnSC(ESRC_LDAERFSPIN,rho(1),rho(2),VLAMBDA,
     &                DFT_SPINDNS,NORDER,d1Ec,d2Ec)

C Complement (NON-Scaled) correlation PBE functional
      ELSE IF (DOLANSC_PBEGWS) THEN
        call ELAC_PBEGWS(rho,grd2,VLAMBDA,Ec)
        call DELANSC(ESRC_PBEGWSERF,rho,grd2,VLAMBDA,NORDER,d1Ec,d2Ec)

C Complement Scaled correlation PBE functional
      ELSE IF (DOLASC_PBEGWS) THEN
        call ELASC_PBEGWS(rho,grd2,VLAMBDA,Ec)
        call DELASC(ESRC_PBEGWSERF,rho,grd2,VLAMBDA,NORDER,d1Ec,d2Ec)

!     short-range spin-dependent correlation LDA with multideterminant reference (for OEP) ! JT 26-08-11
!     ... or for RSDHf calculations (Manu 06-12-2012)
      ELSE IF (DOSRC_MD_LDA) THEN
        call DESRSPIN(ESRC_MD_LDAERF,rho(1),rho(2),
     &             CHIVAL,DFT_SPINDNS,NORDER,Ec,d1Ec,d2Ec)

!     Short-range LDA functional with local CHIVAL parameter ! MNP 13-09-11
      ELSE IF (DOSRC_MULOC_GGA) THEN
C        IF (SRCMULOFAC) THEN
C           xfac = XMULFAC_READIN
C        ELSE
C           xfac = 1.0d0/4.0d0
C        END IF
C        mu = sqrt(grd2(1))/rho(1)
C        mu = xfac*mu
C        mu_new = max(mu,CHIVAL)
        call ESRC_MULOC_GGA(rho,grd2,CHIVAL,Ec)
        call DESR(ESRC_MULOC_GGA,rho,grd2,CHIVAL,NORDER,
     &             Ec, d1Ec, d2Ec)

!     Short-range LDA functional with local CHIVAL parameter ! MNP 13-09-11
      ELSE IF (DOSRC_MULOD_GGA) THEN
        call DESR(ESRC_MULOD_GGA,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

!     Short-range LDA functional with local CHIVAL parameter ! MNP 13-09-11
      ELSE IF(DOSRC_MULOE_GGA) THEN
        call DESR(ESRC_MULOE_GGA,rho,grd2,CHIVAL,NORDER,Ec,d1Ec,d2Ec)

      ENDIF

      RETURN
      END ! subroutine SRDFTEXCJT


C*****************************************************************************
      SUBROUTINE SRDFTLTRJT(KSYMOP,DTXMAT,DSXMAT,EXCMAT,
     &                      WGHTX,WGHTC,GAO,GAO1,CHIVAL,VLAMBDA,
     &                      ERFEXP,RHO,RHG,GRD2,tau,laplace,DGAO,
     &                      DOGGA,do_metagga,TRIPLET,IPRINT_X)
C*****************************************************************************
C     Generate the 1-index transformed DFT Hessian
C     needed in second-order optimization for short-range DFT hybrids.
C
C     Based on the routine SRDFTLTR of J. K. Pedersen
C
C     Created : 15-02-05, J. Toulouse
C     Modified: 17-08-09, J. Toulouse,  add spin density
C     Modified: 16-02-16, E. Hedegaard, add spin gradient
C     Modified: April-18, H. J. Aa. Jensen, new simpler structure for ESRX_ *and ESRC_*
C*****************************************************************************
         implicit none
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"

      integer, intent(in) :: KSYMOP, IPRINT_X
      logical, intent(in) :: DOGGA,ERFEXP(0:2),TRIPLET, do_metagga
      real*8, intent(in) :: rho(4), RHG(3,4), GRD2(6), tau(4),
     &                      laplace(4)
      real*8, intent(in) :: CHIVAL,VLAMBDA
      real*8, intent(in) :: DTXMAT(NBAST,NBAST),DSXMAT(NBAST,NBAST)
      real*8, intent(in) :: GAO(NBAST),GAO1(NBAST,3), DGAO(NBAST,2)
      real*8, intent(in) :: WGHTX,WGHTC
      
      real*8, intent(inout) :: EXCMAT(NBAST,NBAST,2)

      integer, save :: n_call = 0

!     local
      integer :: I, J, ISYM, JSYM, JEND, JSTR, ISTR, IEND
      integer :: NORDER, p, q

      real*8 :: ddot

      real*8 :: BR(4), B0(2), B3(3,2), BMAX, B4(2), B5(2)
      real*8 :: AR(2), AB(2), AX, AY, AZ, RX(2), RY(2), RZ(2)

      real*8 :: FAC0R, FAC0S, FACR, FACS, FAC_AR, FAC_RB, FAC_SB
      real*8 :: FAC_RC, FAC_SC
      real*8 :: G0, GX, GY, GZ, Xi

!     X, C, and XC energies, and derivatives
      real*8 :: Ex, Ec, d1Ex(9), d1Ec(9), d2Ex(45), d2Ec(45)
      real*8 :: Exc, d1Exc(9), d2Exc(45)
!     
!     r = rho, g=gamma, t=tau, e=eta
!     gamma = grad_rho . grad_rho
!     tau   = kinetic energy density
!     eta   = laplace_rho
!
!           dr_c   dr_s   dg_cc   dg_ss   dg_cs   dt_c   dt_s   de_c   de_s
!  d1E?(idx)  1      2      3       4       5       6      7       8      9
!  
!  d2E?(idx)
!           dr_c   dr_s   dg_cc   dg_ss   dg_cs   dt_c   dt_s   de_c   de_s
!     dr_c    1
!     dr_s    2      3
!     dg_cc   4      5      6
!     dg_ss   7      8      9      10
!     dg_cs  11     12     13      14      15
!     dt_c   16     17     18      19      20      21
!     dt_s   22     23     24      25      26      27     28
!     de_c   29     30     31      32      33      34     35      36
!     de_s   37     38     39      40      41      42     43      44     45

      n_call = n_call + 1

      IF (DFT_SPINDNS .AND. TRIPLET) THEN
         CALL QUIT('SRDFTJT.F: Triplet with open-shell not implemented')
      END IF

      ! DGAO(1,1) = DTXMAT*GAO;
      !  if (.not.triplet) DTXMAT is 1-index transformed or transition total (charge) density matrix
      !  if (     triplet) DTXMAT is in fact 1-index transformed or transition spin (triplet) density matrix = DSXMAT
      if (triplet .and. .not.dft_spindns) then
         j = 2 ! triplet (spin)   terms in B0(2), B2(:,2)
      else
         j = 1 ! singlet (charge) terms in B0(1), B2(:,1)
      end if
      CALL DGEMV('N',NBAST,NBAST,1.0d0,
     &           DTXMAT,NBAST,GAO,1,0.0d0,DGAO(1,1),1)

      B0(j) = DDOT(NBAST,DGAO(1,1),1,GAO,1)
      BMAX = ABS(B0(j))
      IF (DOGGA) THEN
C        B3(1:3,j) = GAO1'*DGAO(1,1)
         CALL DGEMV('T',NBAST,3,1.0d0,GAO1,NBAST,
     &         DGAO(1,1),1,0.0d0,B3(1,j),1)

C        DGAO(1,1)= DTXMAT'*GAO !!! HJAAJ: is this necessary ? We have D*XMAT symmetrized!
         CALL DGEMV('T',NBAST,NBAST,1.0d0,DTXMAT,NBAST,
     &         GAO,1,0.0d0,DGAO(1,1),1)
C        B3(1:3,j) = B3(1:3,j) + GAO1(1:NBAST,1:3)'*DGAO((1:NBAST, 1)
         CALL DGEMV('T',NBAST,3,1.0d0,GAO1,NBAST,
     &         DGAO(1,1),1,1.0d0,B3(1,j),1)
         BMAX = MAX(BMAX, ABS(B3(1,j)),ABS(B3(2,j)),ABS(B3(3,j)))
      END IF
      if (do_metagga) then
         B4(j) = 0.0d0
         do p = 1, nbast
            GX = GAO1(p,1)
            GY = GAO1(p,2)
            GZ = GAO1(p,3)
            do q = 1, nbast
               B4(j) = B4(j) + ( GX*GAO1(q,1) + GY*GAO1(q,2) +
     &                           GZ*GAO1(q,3) ) * DTXMAT(p,q)
            end do
         end do
         B4(j) = B4(j)*0.5d0
      end if

      IF (DFT_SPINDNS) THEN
         ! DGAO(1,2) = DSXMAT*GAO; DSXMAT is 1-index transformed or transition spin density matrix
         CALL DGEMV('N',NBAST,NBAST,1.0d0,DSXMAT,NBAST,
     &              GAO,1,0.0d0,DGAO(1,2),1)
         B0(2) = DDOT(NBAST,DGAO(1,2),1,GAO,1)
         BMAX = MAX( BMAX,ABS(B0(2)) )
         IF (DOGGA) THEN
C           B3(2) = GAO1'*DGAO(1,2)
            CALL DGEMV('T',NBAST,3,1.0d0,GAO1,NBAST,
     &         DGAO(1,2),1,0.0d0,B3(1,2),1)

C           DGAO(1,2) = DSXMAT'*GAO !!! HJAAJ: is this necessary ? We have D*XMAT symmetrized!
            CALL DGEMV('T',NBAST,NBAST,1.0d0,DSXMAT,NBAST,
     &         GAO,1,0.0d0,DGAO(1,2),1)
C           B3(1:3,2) = B3(1:3,2) + GAO1(1:NBAST,1:3)T*DGAO((1:NBAST, 2)
            CALL DGEMV('T',NBAST,3,1.0d0,GAO1,NBAST,
     &         DGAO(1,2),1,1.0d0,B3(1,2),1)
            BMAX = MAX(BMAX, ABS(B3(1,2)),ABS(B3(2,2)),ABS(B3(3,2)))
         END IF
         if (do_metagga) then
            B4(2) = 0.0d0
            do p = 1, nbast
               GX = GAO1(p,1)
               GY = GAO1(p,2)
               GZ = GAO1(p,3)
               do q = 1, nbast
                  B4(2) = B4(2) + ( GX*GAO1(q,1) + GY*GAO1(q,2) +
     &                              GZ*GAO1(q,3) ) * DSXMAT(p,q)
               end do
            end do
            B4(2) = B4(2)*0.5d0
         end if
      END IF
C
      IF (BMAX.GT.DFTHRL) THEN
C

C Exchange and correlation derivatives ------------------------------------

      NORDER = 2  ! zero-order, first-order and second-order derivatives
      CALL SRDFTEXCJT(RHO,GRD2,tau,laplace,CHIVAL,NORDER,ERFEXP,triplet,
     &                Ex,d1Ex,d2Ex,Ec,d1Ec,d2Ec,VLAMBDA)
      Exc         = WGHTX*  Ex       + WGHTC*  Ec
      d1Exc(1:9)  = WGHTX*d1Ex(1:9)  + WGHTC*d1Ec(1:9)
      d2Exc(1:45) = WGHTX*d2Ex(1:45) + WGHTC*d2Ec(1:45)

C Linear transformation -----------------------------------------------------

            ! It is used that d2/dxdy = d2/dydx etc. with x and y being the rho_c, gamma_ss, etc.
            ! "_b" in e.g. rhoc_b is for rhoc calculated from DTXMAT (from trial vector B)

            ! B0(1) = rhoc_b
            ! B0(2) = rhos_b
            ! B3(:,1) = grad(rhoc_b)
            ! B3(:,2) = grad(rhos_b)
            ! B4(1) = tauc_b
            ! B4(2) = taus_b
            ! B5(1) = laplace(rhoc_b)
            ! B5(2) = laplace(rhos_b)
            ! BR(1) = grad(rhoc) . grad(rhoc_b)
            ! BR(2) = grad(rhos) . grad(rhos_b)
            ! BR(3) = grad(rhoc) . grad(rhos_b)
            ! BR(4) = grad(rhos) . grad(rhoc_b)
            ! AX, AY, AZ : grad(chi_i chi_j)
            ! BX, BY, BZ : grad(phi_p) . grad(phi_q) 
            ! CX, AC, CZ : laplace(chi_i chi_j)
            ! AR(1) = AX*RHG(1,1)   + AY*RHG(2,1)   + AZ*RHG(3,1)   ! grad(chi_i chi_j) . grad(rhoc)
            ! AR(2) = AX*RHG(1,2)   + AY*RHG(2,2)   + AZ*RHG(3,2)   ! grad(chi_i chi_j) . grad(rhos)
            ! AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1) ! grad(chi_i chi_j) . grad(rhoc_b)
            ! AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2) ! grad(chi_i chi_j) . grad(rhos_b)
            
            ! FAC0R  = All terms in EXCMAT(1) multiplied with G0*GAO(J)
            ! FACR   = All terms in EXCMAT(1) multiplied with AR(1)
            ! FAC0S  = All terms in EXCMAT(2) multiplied with G0*GA0(J)
            ! FACS   = All terms in EXCMAT(2) multiplied with AR(2)
            ! FAC_AR = All terms multiplied with AR(2) for EXCMAT(1),
            !                                and AR(1) for EXCMAT(2)
            ! FAC_RB = All terms multiplied with G1*GAO1(J) in EXCMAT(1)
            ! FAC_SB = All terms multiplied with G1*GAO1(J) in EXCMAT(2)
            ! FAC_RC = All terms multiplied with G2*GAO(J) + G0*GAO2(J)
            !                                    + 2*G1*GAO1(J) in EXCMAT(1)
            ! FAC_SC = All terms multiplied with G2*GAO(J) + G0*GAO2(J)
            !                                    + 2*G1*GAO1(J) in EXCMAT(2)
            !              FAC_RC and FAC_SC, will be implemented later.
            !              all terms related to Laplacian rho, will be 
            !              implemented later. And is therefore commented
            !              out in the below code.

            ! important to check for meta-gga before gga, becuase they
            ! are both true in the case of meta-ggas
            if (dogga) then
               ! grd[rhoc^(X,Y,Z)]
               RX(1) = RHG(1,1)
               RY(1) = RHG(2,1)
               RZ(1) = RHG(3,1)
               if (dft_spindns .or. triplet) then
                  ! grd[rhos^(X,Y,Z)]
                  RX(2) = RHG(1,2)
                  RY(2) = RHG(2,2)
                  RZ(2) = RHG(3,2)
               end if
            end if
            if (do_metagga) then
               if (dft_spindns) then
                  BR(1) = B3(1,1)*RX(1) + B3(2,1)*RY(1) + B3(3,1)*RZ(1)
                  BR(2) = B3(1,2)*RX(2) + B3(2,2)*RY(2) + B3(3,2)*RZ(2)
                  BR(3) = B3(1,2)*RX(1) + B3(2,2)*RY(1) + B3(3,2)*RZ(1)
                  BR(4) = B3(1,1)*RX(2) + B3(2,1)*RY(2) + B3(3,1)*RZ(2)
                  FAC0R = d2Exc(1)*B0(1) + d2Exc(2)*B0(2)
     &                    + 2.0d0*d2Exc(4)*BR(1) + 
     &                    2.0d0*d2Exc(7)*BR(2) + d2Exc(11)*BR(4) +
     &                    d2Exc(11)*BR(3) + d2Exc(16)*B4(1) +
     &                    d2Exc(22)*B4(2)
      !&                    + d2Exc(29)*B5(1) + d2Exc(37)*B5(2)
                  FAC0S = d2Exc(2)*B0(1) + d2Exc(3)*B0(2) 
     &                    + 2.0d0*d2Exc(5)*BR(1) +
     &                    2.0d0*d2Exc(8)*BR(2) + d2Exc(12)*BR(4) +
     &                    d2Exc(12)*BR(3) + d2Exc(17)*B4(1) +
     &                    d2Exc(23)*B4(2)
      !&                    + d2Exc(30)*B5(1) + d2Exc(38)*B5(2)
                  FACR  = 2.0d0*d2Exc(4)*B0(1) + 2.0d0*d2Exc(5)*B0(2) 
     &                    + 4.0d0*d2Exc(6)*BR(1) +
     &                    4.0d0*d2Exc(9)*BR(2) + 2.0d0*d2Exc(13)*BR(4) +
     &                    2.0d0*d2Exc(13)*BR(3) + 2.0d0*d2Exc(18)*B4(1)
     &                    + 2.0d0*d2Exc(24)*B4(2)
      !&                    + 2.0d0*d2Exc(31)*B5(1) + 2.0d0*d2Exc(39)*B5(2)
                  FACS  = 2.0d0*d2Exc(7)*B0(1) + 2.0d0*d2Exc(8)*B0(2) 
     &                    + 4.0d0*d2Exc(9)*BR(1) +
     &                    4.0d0*d2Exc(10)*BR(2) + 2.0d0*d2Exc(14)*BR(4)
     &                   + 2.0d0*d2Exc(14)*BR(3) + 2.0d0*d2Exc(19)*B4(1)
     &                    + d2Exc(25)*B4(2)
      !&                    + 2.0d0*d2Exc(32)*B5(1) + 2.0d0*d2Exc(40)*B5(2)
                  FAC_AR = d2Exc(11)*B0(1) + d2Exc(12)*B0(2) +
     &                     2.0d0*d2Exc(13)*BR(1) + 2.0d0*d2Exc(14)*BR(2)
     &                     + d2Exc(15)*BR(4) + d2Exc(15)*BR(3) +
     &                     d2Exc(20)*B4(1) + d2Exc(26)*B4(2)
      !&                     + d2Exc(33)*B5(1) + d2Exc(41)*B5(2)
                  FAC_RB = d2Exc(16)*B0(1) + d2Exc(17)*B0(2) +
     &                     2.0d0*d2Exc(18)*BR(1) + 2.0d0*d2Exc(19)*BR(2)
     &                     + d2Exc(20)*BR(4) + d2Exc(20)*BR(3) +
     &                     d2Exc(21)*B4(1) + d2Exc(27)*B4(2)
      !&                     + d2Exc(34)*B5(1) + d2Exc(42)*B5(2)
                  FAC_SB = d2Exc(22)*B0(1) + d2Exc(23)*B0(2) +
     &                     2.0d0*d2Exc(24)*BR(1) + d2Exc(25)*BR(2) + 
     &                     d2Exc(26)*BR(4) + d2Exc(26)*BR(3) +
     &                     d2Exc(27)*B4(1) + d2Exc(28)*B4(2)
      !&                     + d2Exc(35)*B5(1) + d2Exc(43)*B5(2)
                  !FAC_RC = d2Exc(29)*B0(1) + d2Exc(30)*B0(2) + 
      !&                     2.0d0*d2Exc(31)*BR(1) + 2.0d0*d2Exc(32)*BR(2) + 
      !&                     d2Exc(33)*BR(4) + d2Exc(33)*BR(3) + 
      !&                     d2Exc(36)*B5(1) + d2Exc(44)*B5(2) +
      !&                     d2Exc(34)*B4(1) + d2Exc(35)*B4(2)
                  !FAC_SC = d2Exc(37)*B0(1) + d2Exc(38)*B0(2) +
      !&                     2.0d0*d2Exc(39)*BR(1) + 2.0d0*d2Exc(40)*BR(2) +
      !&                     d2Exc(41)*BR(4) + d2Exc(41)*BR(3) +
      !&                     d2Exc(44)*B5(1) + d2Exc(45)*B5(2) +
      !&                     d2Exc(42)*B4(1) + d2Exc(43)*B4(2)
               else if (triplet) then
                  BR(3) = B3(1,2)*RX(1) + B3(2,2)*RY(1) + B3(3,2)*RZ(1)
                  FAC0S = d2Exc(3)*B0(2) + d2Exc(12)*BR(3) 
     &                    + d2Exc(23)*B4(2)
      !&                    + d2Exc(38)*B5(2)
                  FAC_AR = d2Exc(12)*B0(2) + d2Exc(15)*BR(3) +
     &                     d2Exc(26)*B4(2)
      !&                     + d2Exc(41)*B5(2)
                  FAC_SB = d2Exc(23)*B0(2) + d2Exc(26)*BR(3) +
     &                     d2Exc(28)*B4(2)
      !&                     + d2Exc(43)*B5(2)
                  !FAC_SC = d2Exc(38)*B0(2) + d2Exc(41)*BR(3) +
      !&                     d2Exc(45)*B5(2) + d2Exc(43)*B4(2)
               else ! closed-shell metaGGA
                  BR(1) = B3(1,1)*RX(1) + B3(2,1)*RY(1) + B3(3,1)*RZ(1)
                  FAC0R = d2Exc(1)*B0(1) + 2.0d0*d2Exc(4)*BR(1) 
     &                    + d2Exc(16)*B4(1)
      !&                    + d2Exc(29)*B5(1)
                  FACR  = 2.0d0*d2Exc(4)*B0(1) + 4.0d0*d2Exc(6)*BR(1) 
     &                   + 2.0d0*d2Exc(18)*B4(1)
      !&                    + 2.0d0*d2Exc(31)*B5(1)
                  FAC_RB = d2Exc(16)*B0(1) + 2.0d0*d2Exc(18)*BR(1) +
     &                     d2Exc(21)*B4(1)
      !&                     + d2Exc(34)*B5(1)
                  !FAC_RC = d2Exc(29)*B0(1) + 2.0d0*d2Exc(31)*BR(1) +
      !&                     d2Exc(36)*B5(1) + d2Exc(34)*B4(1)
               end if
            else if (dogga) then
               if (dft_spindns) then
                  BR(1) = B3(1,1)*RX(1) + B3(2,1)*RY(1) + B3(3,1)*RZ(1)
                  BR(2) = B3(1,2)*RX(2) + B3(2,2)*RY(2) + B3(3,2)*RZ(2)
                  BR(3) = B3(1,2)*RX(1) + B3(2,2)*RY(1) + B3(3,2)*RZ(1)
                  BR(4) = B3(1,1)*RX(2) + B3(2,1)*RY(2) + B3(3,1)*RZ(2)
                  FAC0R = d2Exc(1)*B0(1) + 2.0d0*d2Exc(4)*BR(1) 
     &                    + d2Exc(2)*B0(2) +2.0d0*d2Exc(7)*BR(2)
     &                    + d2Exc(11)*BR(4) + d2Exc(11)*BR(3)
                  FAC0S = d2Exc(2)*B0(1) + d2Exc(3)*B0(2) 
     &                    + 2.0d0*d2Exc(5)*BR(1) +2.0d0*d2Exc(8)*BR(2)
     &                    + d2Exc(12)*BR(4) + d2Exc(12)*BR(3)
                  FACR  = 2.0d0*d2Exc(4)*B0(1) + 4.0d0*d2Exc(6)*BR(1) 
     &                    + 4.0d0*d2Exc(9)*BR(2) +2.0d0*d2Exc(5)*B0(2)
     &                   + 2.0d0*d2Exc(13)*BR(4) + 2.0d0*d2Exc(13)*BR(3)
                  FACS  = 2.0d0*d2Exc(8)*B0(2) + 2.0d0*d2Exc(7)*B0(1) +
     &                    4.0d0*d2Exc(9)*BR(1) + 4.0d0*d2Exc(10)*BR(2)
     &                   + 2.0d0*d2Exc(14)*BR(4) + 2.0d0*d2Exc(14)*BR(3)
                  FAC_AR = d2Exc(11)*B0(1) + 2.0d0*d2Exc(13)*BR(1) +
     &                      d2Exc(12)*B0(2) + 2.0d0*d2Exc(14)*BR(2) +
     &                      d2Exc(15)*BR(3) + d2EXC(15)*BR(4)
               else if (triplet) then
                  BR(3) = B3(1,2)*RX(1) + B3(2,2)*RY(1) + B3(3,2)*RZ(1)
                  FAC0S = d2Exc(3)*B0(2) + d2Exc(12)*BR(3)
                  FAC_AR = d2Exc(12)*B0(2) + d2Exc(15)*BR(3)
               else ! normal GGA
                  BR(1) = B3(1,1)*RX(1) + B3(2,1)*RY(1) + B3(3,1)*RZ(1)
                  FAC0R = d2Exc(1)*B0(1) + 2.0d0*d2Exc(4)*BR(1)
                  FACR  = 2.0d0*d2Exc(4)*B0(1) + 4.0d0*d2Exc(6)*BR(1)
               end if
            else !LDA
               if (dft_spindns) then
                  FAC0R = d2Exc(1)*B0(1) + d2Exc(2)*B0(2)
                  FAC0S = d2Exc(2)*B0(1) + d2Exc(3)*B0(2)
               else if (triplet) then
                  FAC0S = d2Exc(3)*B0(2)
               else ! normal LDA
                  FAC0R = d2Exc(1)*B0(1)
               end if
            end if

#if SRDFT_DEBUG > 11
           if (modulo(n_call, 1000) .eq. 0) then
              write(lupri,*) 'hjaaj: rho ',rho(1:4)
              write(lupri,*) 'hjaaj: rhgc',rhg(1:3,1)
              write(lupri,*) 'hjaaj: grd2',grd2(1:6)
              write(lupri,*) 'debug: tau',tau(1:4)
              write(lupri,*) 'hjaaj: Ex, Ec, gga, spindns, triplet',
     &          Ex, Ec, dogga, dft_spindns, triplet
              write(lupri,*) 'debug: do_metagga',do_metagga
              write(lupri,*) 'hjaaj: d1Ex  ',d1Ex(:)
              write(lupri,*) 'hjaaj: d1Ec  ',d1Ec(:)
              write(lupri,*) 'hjaaj: d2Ex(:)',d2Ex(:)
              write(lupri,*) 'hjaaj: d2Ec(:)',d2Ec(:)

              write(lupri,*) 'hjaaj: d2Exc(1)',d2Exc(1)
              write(lupri,*) 'hjaaj:d2Exc(4),d2Exc(6)',d2Exc(4),d2Exc(6)
              write(lupri,*) 'hjaaj:d2Exc(3),d2Exc(2)',d2Exc(3),d2Exc(2)
              write(lupri,*) 'hjaaj:d1Exc(3),d1Exc(4)',d1Exc(3),d1Exc(4)
            write(lupri,*) 'hjaaj:d2Exc(8),d2Exc(7),d2Exc(9),d2Exc(10)',
     &                             d2Exc(8),d2Exc(7),d2Exc(9),d2Exc(10)
              write(lupri,*) 'debug: d1Exc(1:9)', d1Exc(1:9)
              write(lupri,*) 'debug: d2Exc(1:45)', d2Exc(1:45)
              write(lupri,*) 'hjaaj:RX,RY,RZ',RHG(1,1),RHG(2,1),RHG(3,1)
              write(lupri,*) 'hjaaj: B0, BMAX',B0(1:2), BMAX
              write(lupri,*) 'hjaaj: B4',B4(1:2)
              write(lupri,*) 'hjaaj: B3 matrix',B3(1:3,1:2)
              call output(b3,1,3,1,2,3,2,1,lupri)
              write(lupri,*) 'hjaaj: BR',BR(1:2)
              write(lupri,*) 'hjaaj: FAC0R, FAC0S, FACR, FACS, FAC_AR',
     &                               FAC0R, FAC0S, FACR, FACS, FAC_AR
              write(lupri,*) 'debug: FAC_RB, FAC_SB',FAC_RB, FAC_SB
              write(lupri,*) 'hjaaj: current EXCMAT left corner; before'
              call output(excmat,1,5,1,5,nbast,nbast,1,lupri)

              write(lupri,*) 'debug: d2Exc(23)*B4(2)', d2Exc(23)*B4(2)
              write(lupri,*) 'debug: d2Exc(26)*B4(2)', d2Exc(26)*B4(2)
              write(lupri,*) 'debug: d2Exc(23)*B0(2)', d2Exc(23)*B0(2)
              write(lupri,*) 'debug: d2Exc(26)*BR(3)', d2Exc(26)*BR(3)
              write(lupri,*) 'debug: d2Exc(28)*B4(2)', d2Exc(28)*B4(2)
              write(lupri,*) 'debug: d2Exc(16)*B4(1)', d2Exc(16)*B4(1)
              write(lupri,*) 'debug: d2Exc(18)*B4(1)', d2Exc(18)*B4(1)
              write(lupri,*) 'debug: d2Exc(16)*B0(1)', d2Exc(16)*B0(1)
              write(lupri,*) 'debug: d2Exc(18)*BR(1)', d2Exc(18)*BR(1)
              write(lupri,*) 'debug: d2Exc(21)*B4(1)', d2Exc(21)*B4(1)
           end if
#endif
            IF (NSYM.EQ.1) THEN
               if (do_metagga) then
                  if (triplet) then
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        GX = GAO1(I,1)
                        GY = GAO1(I,2)
                        GZ = GAO1(I,3)
                        DO J = 1, I
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           Xi = GX*GAO1(J,1)+GY*GAO1(J,2)+GZ*GAO1(J,3)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1) +
     &                         FAC0S*G0*GAO(J) + FAC_AR*AR(1) +
     &                         2.0d0*d1Exc(4)*AB(2) + FAC_SB*Xi*0.5d0
      !&                         + FAC_SC*??
                        END DO
                     END DO
                  else if (dft_spindns) then
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        GX = GAO1(I,1)
                        GY = GAO1(I,2)
                        GZ = GAO1(I,3)
                        DO J = 1, I
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AR(2) = AX*RX(2) + AY*RY(2) + AZ*RZ(2)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           Xi = GX*GAO1(J,1)+GY*GAO1(J,2)+GZ*GAO1(J,3)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1) +
     &                         FAC0R*G0*GAO(J) + FACR*AR(1) +
     &                         FAC_RB*Xi*0.5d0 + FAC_AR*AR(2) + 
     &                         2.0d0*d1Exc(3)*AB(1) + d1Exc(5)*AB(2)
      !&                         + FAC_RC*??
                           EXCMAT(J,I,2) = EXCMAT(J,I,2) + 
     &                         FAC0S*G0*GAO(J) + FACS*AR(2) +
     &                         FAC_SB*Xi*0.5d0 + FAC_AR*AR(1) +
     &                         2.0d0*d1Exc(4)*AB(2) + d1Exc(5)*AB(1)
      !&                         + FAC_SC*??
                        END DO
                     END DO
                  else ! Closed-shell meta-GGA
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        GX = GAO1(I,1)
                        GY = GAO1(I,2)
                        GZ = GAO1(I,3)
                        DO J = 1, I
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           Xi = GX*GAO1(J,1)+GY*GAO1(J,2)+GZ*GAO1(J,3)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1) +
     &                         FAC0R*G0*GAO(J) + FACR*AR(1) 
     &                         + 2.0d0*d1Exc(3)*AB(1)
     &                         + FAC_RB*Xi*0.5d0
      !&                         + FAC_RC*??
                        END DO
                     END DO
                  end if
               ELSE IF (DOGGA) THEN
                  IF (TRIPLET) THEN ! closed-shell triplet GGA response
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        GX = GAO1(I,1)
                        GY = GAO1(I,2)
                        GZ = GAO1(I,3)
                        DO J = 1, I
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           ! spin matrix W[b]2s in EXCMAT no. 1
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0S*G0*GAO(J)
     &                                     + 2.0d0*d1Exc(4)*AB(2) 
     &                                     + FAC_AR*AR(1)
                        END DO
                     END DO
                  ELSE IF (DFT_SPINDNS) THEN ! open-shell GGA response
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        GX = GAO1(I,1)
                        GY = GAO1(I,2)
                        GZ = GAO1(I,3)
                        DO J = 1, I
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AR(2) = AX*RX(2) + AY*RY(2) + AZ*RZ(2)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           ! charge matrix W[b]2c in EXCMAT no. 1
                           ! spin   matrix W[b]2s in EXCMAT no. 2
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
     &                       + d1Exc(5)*AB(2) + 2.0d0*d1Exc(3)*AB(1) 
     &                       + FACR*AR(1) + FAC_AR*AR(2)
                           EXCMAT(J,I,2) = EXCMAT(J,I,2)+FAC0S*G0*GAO(J)
     &                       + d1Exc(5)*AB(1) + 2.0d0*d1Exc(4)*AB(2) 
     &                       + FACS*AR(2) + FAC_AR*AR(1)
                        END DO
                     END DO
                  ELSE ! closed-shell singlet GGA response
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        GX = GAO1(I,1)
                        GY = GAO1(I,2)
                        GZ = GAO1(I,3)
                        DO J = 1, I
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           ! charge matrix W[b]2c in EXCMAT no. 1
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
     &                                   + FACR*AR(1) 
     &                                   + 2.0d0*d1Exc(3)*AB(1)
                        END DO
                     END DO
                  END IF
               ELSE ! Doing LDA calculations from here
                  IF (TRIPLET) THEN
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        DO J = 1, I
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0S*G0*GAO(J)
                        END DO
                     END DO
                  ELSE IF (DFT_SPINDNS) THEN
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        DO J = 1, I
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
                           EXCMAT(J,I,2) = EXCMAT(J,I,2)+FAC0S*G0*GAO(J)
                        END DO
                     END DO
                  ELSE ! closed-shell LDA
                     DO I = 1, NBAST
                        G0 = GAO(I)
                        DO J = 1, I
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
                        END DO
                     END DO
                  END IF
               END IF ! GGA or LDA check
            ELSE ! NSYM .gt. 1
               ! Only two indentations is used in symmstry code insted of three.
               !   This is to not exceed 72 chars all the time. - E.R.K
               DO ISYM = 1, NSYM
                 ISTR = IBAS(ISYM) + 1
                 IEND = IBAS(ISYM) + NBAS(ISYM)
                 JSYM = MULD2H(ISYM,KSYMOP)
                 IF (ISYM.GE.JSYM) THEN
                   JSTR = IBAS(JSYM) + 1
                   JEND = IBAS(JSYM) + NBAS(JSYM)
                   IF (do_metagga) THEN
                     IF (triplet) THEN
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         GX = GAO1(I,1)
                         GY = GAO1(I,2)
                         GZ = GAO1(I,3)
                         DO J = JSTR, min(I, JEND)
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           Xi = GX*GAO1(J,1)+GY*GAO1(J,2)+GZ*GAO1(J,3)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1) +
     &                         FAC0S*G0*GAO(J) + FAC_AR*AR(1) +
     &                         2.0d0*d1Exc(4)*AB(2) + FAC_SB*Xi*0.5d0
      !&                         + FAC_SC*??
                         END DO
                       END DO
                     ELSE IF (dft_spindns) THEN
                       do i = istr, iend
                         G0 = GAO(I)
                         GX = GAO1(I,1)
                         GY = GAO1(I,2)
                         GZ = GAO1(I,3)
                         DO J = JSTR, min(I, JEND)
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AR(2) = AX*RX(2) + AY*RY(2) + AZ*RZ(2)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           Xi = GX*GAO1(J,1)+GY*GAO1(J,2)+GZ*GAO1(J,3)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1) +
     &                         FAC0R*G0*GAO(J) + FACR*AR(1) +
     &                         FAC_RB*Xi*0.5d0 + FAC_AR*AR(2) + 
     &                         2.0d0*d1Exc(3)*AB(1) + d1Exc(5)*AB(2)
      !&                         + FAC_RC*??
                           EXCMAT(J,I,2) = EXCMAT(J,I,2) + 
     &                         FAC0S*G0*GAO(J) + FACS*AR(2) +
     &                         FAC_SB*Xi*0.5d0 + FAC_AR*AR(1) +
     &                         2.0d0*d1Exc(4)*AB(2) + d1Exc(5)*AB(1)
      !&                         + FAC_SC*??
                         END DO
                       END DO
                     ELSE ! Closed-shell MGGA singlet response
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         GX = GAO1(I,1)
                         GY = GAO1(I,2)
                         GZ = GAO1(I,3)
                         do J = JSTR, min(I, JEND)
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           Xi = GX*GAO1(J,1)+GY*GAO1(J,2)+GZ*GAO1(J,3)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1) +
     &                         FAC0R*G0*GAO(J) + FACR*AR(1) 
     &                         + 2.0d0*d1Exc(3)*AB(1)
     &                         + FAC_RB*Xi*0.5d0
      !&                         + FAC_RC*??
                         END DO
                       END DO
                     END IF
                   ELSE IF (DOGGA) THEN
                     IF (TRIPLET) THEN
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         GX = GAO1(I,1)
                         GY = GAO1(I,2)
                         GZ = GAO1(I,3)
                         DO J = JSTR, MIN(I,JEND)
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           ! spin matrix W[b]2s in EXCMAT no. 1
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0S*G0*GAO(J)
     &                                     + 2.0d0*d1Exc(4)*AB(2) 
     &                                     + FAC_AR*AR(1)
                         END DO
                       END DO
                     ELSE IF (DFT_SPINDNS) THEN
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         GX = GAO1(I,1)
                         GY = GAO1(I,2)
                         GZ = GAO1(I,3)
                         DO J = JSTR, MIN(I,JEND)
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AR(2) = AX*RX(2) + AY*RY(2) + AZ*RZ(2)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
                           ! charge matrix W[b]2c in EXCMAT no. 1
                           ! spin   matrix W[b]2s in EXCMAT no. 2
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
     &                       + d1Exc(5)*AB(2) + 2.0d0*d1Exc(3)*AB(1) 
     &                       + FACR*AR(1) + FAC_AR*AR(2)
                           EXCMAT(J,I,2) = EXCMAT(J,I,2)+FAC0S*G0*GAO(J)
     &                       + d1Exc(5)*AB(1) + 2.0d0*d1Exc(4)*AB(2) 
     &                       + FACS*AR(2) + FAC_AR*AR(1)
                         END DO
                       END DO
                     ELSE ! Closed-shell GGA singlet response
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         GX = GAO1(I,1)
                         GY = GAO1(I,2)
                         GZ = GAO1(I,3)
                         DO J = JSTR, MIN(I,JEND)
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
                           ! charge matrix W[b]2c in EXCMAT no. 1
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
     &                                     + FACR*AR(1) 
     &                                     + 2.0d0*d1Exc(3)*AB(1)
                         END DO
                       END DO
                     END IF
                   ELSE ! Doing LDA calculations from here
                     IF (TRIPLET) THEN
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         DO J = JSTR, MIN(I,JEND)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0S*G0*GAO(J)
                         END DO
                       END DO
                     ELSE IF (DFT_SPINDNS) THEN
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         DO J = JSTR, MIN(I,JEND)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
                           EXCMAT(J,I,2) = EXCMAT(J,I,2)+FAC0S*G0*GAO(J)
                         END DO
                       END DO
                     ELSE
                       DO I = ISTR, IEND
                         G0 = GAO(I)
                         DO J = JSTR, MIN(I,JEND)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
                         END DO
                       END DO
                     END IF
                   END IF
                  END IF
               END DO
            END IF
      END IF   ! IF (BMAX.GT.DFTHRL) THEN
#if SRDFT_DEBUG > 11
           if (modulo(n_call, 1000) .eq. 0) then
         write(lupri,*) ':: current EXCMAT left corner - after'
         call output(excmat,1,5,1,5,nbast,nbast,1,lupri)
           end if
#endif
      RETURN
      END
! end of srdftjt.F

