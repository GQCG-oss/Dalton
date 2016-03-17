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

#define noSRDFT_DEBUG

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRDFT1JT(ND_SIM,EXCMAT,EDFTY,DMAT,DOLND,CORX,CORY,CORZ,
     &                  WEIGHT,NBUF,GSO,NCNT,DMAGAOSR,DOGGA,
     &                  DOERG,DOGRD,DOATR,WORK,LWORK,IPRINT_inp)
C*****************************************************************************
C     Calculate DFT energy and potential contributions for
C     SR-DFT hybrid methods
C
C     Based on the routine SRDFT1 of J. K. Pedersen
C
C     Created : 15-02-05, J. Toulouse
C     Modified: 17-08-09, J. Toulouse, add spin density
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
#include "dummy.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "inftap.h"
#include "infpar.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"
#include "gnrinf.h"
#include "chrnos.h"

      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)

      LOGICAL   DOLND,DOGGA,DOERG,DOGRD,DOATR,FROMVX

#ifdef SRDFT_DEBUG
      !edh debug
      real*8    print_thr
#endif

      DIMENSION EXCMAT(NBAST,NBAST,*), DMAT(N2BASX,*), 
     &          CORX(NBUF),CORY(NBUF),CORZ(NBUF),WEIGHT(NBUF),
     &          GSO(NBAST*NTYPSO),NCNT(NBAST),
     &          DMAGAOSR(NBAST,2),
     &          COR(3),WORK(LWORK)
C
      DIMENSION RHO(4), RHG(3,4)
C     ... RHO(1) = RHOT  = rho_a + rho_b
C     ... RHO(2) = RHOS  = rho_a - rho_b
C     ... RHO(3) = RHO2A = 2 rho_a
C     ... RHO(4) = RHO2B = 2 rho_b
C     ... RHG(1:3,1) = RHGT(1:3)  = grad ( rho_a + rho_b ) = grad ( rho_t )
C     ... RHG(1:3,2) = RHGS(1:3)  = grad ( rho_a - rho_b ) = grad ( rho_s )
C     ... RHG(1:3,3) = RHG2A(1:3) = 2 grad ( rho_a )
C     ... RHG(1:3,4) = RHG2B(1:3) = 2 grad ( rho_b )

!     these are not initialized to zero!
      double precision d1Ex(6), d1Ec(6), d2Ex(10), d2Ec(10)
      double precision d1Exdrhot, d1Exdrhos, d1Exdgrdt, d1Exdgrds
      double precision d1Ecdrhot, d1Ecdrhos, d1Ecdgrdt, d1Ecdgrds
      double precision d1Exdngrd2, d1Exdngrds2
      double precision d1Ecdngrd2, d1Ecdngrds2

      double precision d1Excdrhot, d1Excdrhos, d1Excdgrdt, d1Excdgrds
      double precision d1Excdngrd2, d1Excdngrds2


      CHARACTER*16 QUADNAME

      IPRINT = IPRINT_inp
!#ifdef SRDFT_DEBUG
      IPRINT = MAX(IPRINT,5)
!     Debug
      WRITE(LUPRI,*) 'SRDFT1JT: entering, ND_SIM = ',ND_SIM
      WRITE(LUPRI,*) 'SRDFT1JT: LUPRI  =',LUPRI, IPRINT
      WRITE(LUPRI,*) 'SRDFT1JT: CHIVAL =',CHIVAL
      WRITE(LUPRI,*) 'SRDFT1JT: HFXFAC =',HFXFAC
      WRITE(LUPRI,*) 'SRDFT1JT: DOERG,DOATR,DOGRD,DOGGA  =',
     &   DOERG,DOATR,DOGRD,DOGGA
      WRITE(LUPRI,*) 'SRDFT1JT: ERFEXP(0:2) =',ERFEXP(0:2)
      WRITE(LUPRI,*) 'SRDFT1JT: DOSRX_PBERI=',DOSRX_PBERI
      WRITE(LUPRI,*) 'DFT_SPINDNS=',DFT_SPINDNS
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
C
!#endif
C     Initializations -------------------------------------------------------

      ! edh debug
#ifdef SRDFT_DEBUG
      print_thr = 0.0d0
#endif
      ! edh the arrays are also zeroed in SRDFTEXCJT
      d1Ex(1:6)  = 0.0d0
      d1Exdrhot  = 0.0d0
      d1Exdrhos  = 0.0d0
      d1Exdgrdt  = 0.0d0
      d1Exdgrds  = 0.0d0
      d1Exdngrd2 = 0.0d0
      d1Exdngrds2 = 0.0d0
      !----------------
      d1Ec(1:6)  = 0.0d0
      d1Ecdrhot  = 0.0d0
      d1Ecdrhos  = 0.0d0
      d1Ecdgrdt  = 0.0d0
      d1Ecdgrds  = 0.0d0
      d1Ecdngrd2 = 0.0d0
      d1Ecdngrds2 = 0.0d0
      !---------------
      d1Excdrhot  = 0.0d0
      d1Excdrhos  = 0.0d0
      d1Excdgrdt  = 0.0d0
      d1Excdgrds  = 0.0d0
      d1Excdngrd2 = 0.0d0
      d1Excdngrds2 = 0.0d0


      IF (ND_SIM .NE. 1) THEN
         WRITE(LUPRI,*) 'SRDFT1: illegal ND_SIM ',ND_SIM
         CALL QUIT('SRDFT1: illegal ND_SIM value')
      END IF

      IF (DOGRD) THEN
         CALL QUIT('SRDFTJT: DOGRD not implemented yet')
      END IF

!     (Integrated) exchange and correlation energies
      ExInt  = D0
      EcInt  = D0
      IF (DOERG) EDFTY  = D0

C     DFT grid (NTOTAL is the total number of grid points)
      IF (.NOT.GRDONE_OLD) THEN
         CALL DFTGRD(WORK,LWORK,NTOTAL,1,.FALSE.)
         CALL CONDFT
         GRDONE_OLD = .TRUE.
      END IF

C     Electron number
      ELCTRN = D0

C     Integrated spin density
      SPINDEN = D0

C     Thresholds for screening
C     THRELCTRN makes sure that error in no. of electrons
C     from grid integration less than DFTHR
      THRELCTOT = 1.0D-7
      THRELCTRN = D0
      THRESRHO  = 1.0D-8
      THRDRHO   = D0
      DFTHR0    = D0
      DFTHRI    = D0
      DFTHRL    = D0

C     Skipped grid points and nb of electrons
      NSKIP1 = 0     
      ELCTRSKIP = D0 
C
C     Counter for how many times RHO(2) > RHO(1)
      NRSGTR = 0

C     For quadrature integration
C
C     Make quadname: Can take 9999 procs
C
      QUADNAME = 'DALTON.QUAD.'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
C
      LUQUAD = -1
      CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)
C
      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
      kpoint = npoint/10

      IF (NPOINT.GT.0) THEN
C      jpnt = 0  !debug
      NPNTS = NPNTS + NPOINT
      CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)

C   Loop over grid points -----------------------------------------------
      DO 300 IPNT = 1, NPOINT

C     Print grid information
      IF (IPRINT .GT. 100) THEN
        WRITE (LUPRI,'(2X,I6,4F12.6)') 
     &        IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
      END IF

C     Quadrature weight
      WGHT  = WEIGHT(IPNT)
      WGHTX = (1.0D0-HFXFAC)*WGHT
      WGHTC = WGHT

C     Get AOs
      THRINT = DFTHRI/WGHT
      COR(1) = CORX(IPNT)
      COR(2) = CORY(IPNT)
      COR(3) = CORZ(IPNT)
      CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,NBAST,DOLND,DOGGA,
     &            THRINT,IPRINT)

C     Density
      CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHO(1),RHOT13,
     &       DMAGAOSR(1,1),THRINT,IPRINT)

C     Contribution to electron number (N_alpha + N_beta)
      DELCTRN = WGHT*RHO(1)

C     Screening
!      print *,'DELCTRN=',DELCTRN
!      print *,'THRELCTRN=',THRELCTRN
      IF (ABS(DELCTRN).LE.THRELCTRN .OR. 
     &          RHO(1).LE.THRESRHO) THEN
         NSKIP1 = NSKIP1 + 1
         ELCTRSKIP = ELCTRSKIP + DELCTRN
         THRELCTRN = (THRELCTOT-ELCTRSKIP)/100
         GO TO 300
C        ... skip this point, go to next grid point
      END IF

C     Spin density
      IF (DFT_SPINDNS) THEN
C        DMAT(1,3) is the spin density matrix
         CALL GETRHO_OLD(DMAT(1,3),GSO(KSO0),RHO(2),RHOS13,
     &                   DMAGAOSR(1,2),THRINT,IPRINT)
C        In our local spin model we can get unphysical densities
C        because we allow DS > DV in some regions.          
         IF (DFT_LOCALSPIN) THEN
            IF (RHO(2) .GT. RHO(1)) THEN
                RHO(2) = RHO(1)
                NRSGTR = NRSGTR + 1 
            END IF
         END IF
C        Contribution to spin density (N_alpha - N_beta)
         SPINDEN = SPINDEN + WGHT*RHO(2)
      ELSE
        RHO(2) = 0.0D0
      END IF
      RHO(3) = RHO(1) + RHO(2)
      RHO(4) = RHO(1) - RHO(2)

C     Gradient of density
      IF (DOGGA) THEN
C       GSO(KSO1)*DMAGAOSR
        CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAOSR(1,1),
     &     1,D0,RHG(1,1),1)
C
        IF (DFT_SPINDNS) THEN
           CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAOSR(1,2),
     &        1,D0,RHG(1,2),1)
        ELSE
           RHG(1,2) = 0.0D0 ! RHGS is zero
           RHG(2,2) = 0.0D0
           RHG(3,2) = 0.0D0
        ENDIF
        DO I = 1,3
           RHG(I,3) = RHG(I,1) + RHG(I,2) ! RHG2A
           RHG(I,4) = RHG(I,1) - RHG(I,2) ! RHG2B
        END DO
      ENDIF
C
C     Number of electrons
      ELCTRN = ELCTRN + DELCTRN

      IF (DOERG) THEN
C
C     Exchange and correlation energies and derivatives
      NORDER = 1  ! first-order derivatives
C
      CALL SRDFTEXCJT(RHO,RHG,CHIVAL,NORDER,ERFEXP,
     >                Ex,d1Ex,d2Ex,Ec,d1Ec,d2Ec,VLAMBDA,ipnt)
C      Summation with quadrature weigths
       ExInt     = ExInt + WGHTX*Ex
       EcInt     = EcInt + WGHTC*Ec

C      Potentials with quadrature weigths
!
       if (DOGGA .and. DFT_SPINDNS) then
           if (ipnt .eq. 1) write(lupri,*)'case 1 active'
           d1Exdrhot   = d1Ex(1)
           d1Exdgrdt   = d1Ex(2)
           d1Exdngrd2  = d1Ex(3)
           d1Exdrhos   = d1Ex(4)
           d1Exdgrds   = d1Ex(5)
           d1Exdngrds2 = d1Ex(6)
           !-----------------
           d1Ecdrhot   = d1Ec(1)
           d1Ecdgrdt   = d1Ec(2)
           d1Ecdngrd2  = d1Ec(3)
           d1Ecdrhos   = d1Ec(4)
           d1Ecdgrds   = d1Ec(5)
           d1Ecdngrds2 = d1Ec(6)
       else if (DFT_SPINDNS .and. .not. DOGGA) then
           if (ipnt .eq. 1) write(lupri,*)'case 2 active'
           d1Exdrhot  = d1Ex(1)
           d1Exdrhos  = d1Ex(2)
           !-----------------
           d1Ecdrhot  = d1Ec(1)
           d1Ecdrhos  = d1Ec(2)
       else if (DOGGA .and. .not. DFT_SPINDNS) then
           if (ipnt .eq. 1) write(lupri,*)'case 3 active'
           d1Exdrhot  = d1Ex(1)
           d1Exdgrdt  = d1Ex(2)
           d1Exdngrd2 = d1Ex(3)
           !-----------------
           d1Ecdrhot  = d1Ec(1) 
           d1Ecdgrdt  = d1Ec(2)
           d1Ecdngrd2 = d1Ec(3)
       else
           if (ipnt .eq. 1) write(lupri,*)'case 4 active'
           d1Exdrhot  = d1Ex(1)
           d1Ecdrhot  = d1Ec(1)
       end if

       d1Exdrhot  = WGHTX*d1Exdrhot
       d1Exdrhos  = WGHTX*d1Exdrhos
       d1Exdgrdt  = WGHTX*d1Exdgrdt
       d1Exdngrd2 = WGHTX*d1Exdngrd2
       d1Exdgrds  = WGHTX*d1Exdgrds
       d1Exdngrds2 = WGHTX*d1Exdngrd2

       d1Ecdrhot  = WGHTC*d1Ecdrhot
       d1Ecdrhos  = WGHTC*d1Ecdrhos
       d1Ecdgrdt  = WGHTC*d1Ecdgrdt
       d1Ecdngrd2 = WGHTC*d1Ecdngrd2
       d1Ecdgrds  = WGHTC*d1Ecdgrds
       d1Ecdngrds2 = WGHTC*d1Ecdngrd2

C      Sum of exchange and correlation
       d1Excdrhot  = d1Exdrhot  + d1Ecdrhot
       d1Excdrhos  = d1Exdrhos  + d1Ecdrhos
       d1Excdgrdt  = d1Exdgrdt  + d1Ecdgrdt
       d1Excdngrd2 = d1Exdngrd2 + d1Ecdngrd2
       d1Excdgrds  = d1Exdgrds  + d1Ecdgrds
       d1Excdngrds2 = d1Exdngrds2 + d1Ecdngrds2

       if (ipnt .eq. 1) then
          write(lupri,*)'d1Excdrhot', d1Excdrhot
          write(lupri,*)'d1Exdrhot', d1Exdrhot
          write(lupri,*)'d1Ecdrhot', d1Ecdrhot

          write(lupri,*)'d1Excdgrdt', d1Excdgrdt
          write(lupri,*)'d1Exdgrdt', d1Exdgrdt
          write(lupri,*)'d1Ecdgrdt', d1Ecdgrdt

          write(lupri,*)'d1Excdngrd2', d1Excdngrd2
          write(lupri,*)'d1Exdngrd2', d1Exdngrd2
          write(lupri,*)'d1Ecdngrd2', d1Ecdngrd2
       end if

C      Exchange-correlation contribution to Kohn-Sham matrix (charge and spin components)
       ! Ignore spin in start CI ...
       if (ipnt .eq. 1) then
          write(lupri,*)'DOERG,DOATR',DOERG,DOATR
          write(lupri,*)'excmat before DFTKSM'
         !CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
       end if

       FROMVX = .FALSE.
       if (domcsrdft .and. docisrdft) then        ! (DOCISRDFT .and. DOMCSRDFT) is only true for start-CI  (cf. routine OPST)
           CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG(1,1),
     &                 d1Excdrhot,d1Excdngrd2,DOGGA,FROMVX,DFTHRL,IPNT)
           if (ipnt .eq. 1) then
              write(lupri,*)'excmat after DFTKSM'
              !CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
           endif
       else if (domcsrdft .and. .not. docisrdft) then
          if (dogga. and. dft_spindns) then       ! special treatment if gga + spin
              CALL DFTKSMGGASPIN(EXCMAT,GSO(KSO0),GSO(KSO1),RHG,
     &             d1Excdrhot,d1Excdngrd2,d1Excdngrd2,d1Excdngrds2,
     &             DFTHRL,IPNT)
          else                                    ! pure GGA or pure spin
             CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG(1,1),
     &                d1Excdrhot,d1Excdngrd2,DOGGA,FROMVX,DFTHRL, IPNT)
             if (dft_spindns .and. .not. dogga) then
                 CALL DFTKSM(EXCMAT(1,1,2),GSO(KSO0),GSO(KSO1),RHG(1,2),
     &                       d1Excdrhos,0.0d0,DOGGA,FROMVX,DFTHRL,IPNT)
             end if
          end if
       end if
C
C
       ENDIF  !DOERG

C      Hessian transformation for second-order terms
       IF (DOATR) THEN
C        matrix 1 is DTAO
C        matrix 2 is DXTAO
C        matrix 3 is DSAO
C        matrix 4 is DXSAO
         JDXTAO = 2
         JDXSAO = 4

#ifdef SRDFT_DEBUG
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
        if (ipnt .eq. 1) then
            write(lupri,*)'edh before SRDFTLTRJT' 
            WRITE(LUPRI,*) ' <<< Kohn-Sham matrix >>>'
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
            WRITE(LUPRI,*)' <<< Kohn-Sham matrix from spin-density>>>'
            CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,
     &                  NBAST,NBAST,NBAST,-1,LUPRI)
        end if

        CALL SRDFTLTRJT(JWOPSY,DMAT(1,JDXTAO),DMAT(1,JDXSAO),EXCMAT,
     &                 WGHTX,WGHTC, GSO(KSO0),GSO(KSO1),
     &                 CHIVAL,VLAMBDA,ERFEXP,RHO,RHG,DMAGAOSR,DOGGA,
     &                 IPNT)

        if (ipnt .eq. 1) then
            write(lupri,*)'edh after SRDFTLTRJT'
            WRITE(LUPRI,*)'<<< Kohn-Sham matrix >>>'
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
            WRITE(LUPRI,*)' <<< Kohn-Sham matrix from spin-density>>>'
            CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,
     &                  NBAST,NBAST,NBAST,-1,LUPRI)
        end if
    
       ENDIF

#ifdef SRDFT_DEBUG
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

C
  300  CONTINUE
C End of loop over points ------------------------------------------------
C
       GO TO 200
       ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
       END IF  !NPOINT.GT.0

C Termination --------------------------------------------------------------

#ifdef SRDFT_DEBUG
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

      IF (DOERG) THEN
         EDFTY = ExInt + EcInt
         FAC = DP5
      ELSE IF (DOATR) THEN
         FAC = D1
      ELSE
         CALL QUIT('both DOERG and DOATR are false')
      END IF

      IF (DFT_SPINDNS) THEN
        NXCMAT = 2
      ELSE
        NXCMAT = 1
      ENDIF

      DO K = 1, NXCMAT
        DO I = 1, NBAST
          DO J = 1, I - 1
            AVERAG = FAC*(EXCMAT(I,J,K) + EXCMAT(J,I,K))
            EXCMAT(I,J,K) = AVERAG
            EXCMAT(J,I,K) = AVERAG
          END DO
        END DO
      END DO

C     Print section
!     IF (IPRINT .GE. 2) THEN
      IF (IPRINT .GE. 1) THEN
      IF (DOERG) THEN
         ELCTOT = ELCTRN + ELCTRSKIP
         WRITE (LUPRI,'(4(/1X,A,F18.10))')
     &'No. of non-skipped electrons from num. int.  :',ELCTRN,
     &'No. of     skipped electrons   in num. int.  :',ELCTRSKIP,
     &'Sum of non-skipped and skipped electrons     :',ELCTOT,
     &'Integrated spin density (N_alpha-N_beta)     :',SPINDEN

         WRITE (LUPRI,'(/1X,A,F18.10)')
     &'Exchange energy                                     :',ExInt
         WRITE (LUPRI,'(1X,A,F18.10)')
     &'Correlation energy                                  :',EcInt
         WRITE (LUPRI,'(1X,A,F18.10)')
     &'DFT exchange-correlation energy                     :',EDFTY
         WRITE (LUPRI,*)
     &'Number of times RHO(2) > RHO(1)                     :',NRSGTR
      END IF 

      IF (DOATR) THEN
         WRITE(LUPRI,'(/A)') ' srDFT linear transformation finished'
      END IF

!      IF (IPRINT.GE.5) THEN
            WRITE(LUPRI,'(//A/)') ' <<< Kohn-Sham matrix >>>'
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
            IF (DFT_SPINDNS) THEN
               WRITE(LUPRI,'(//A/)')
     &         ' <<< Kohn-Sham matrix from spin-density>>>'
               CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,NBAST,NBAST,NBAST,
     &         -1,LUPRI)
            END IF
!      END IF

      END IF  ! IPRINT .GE. 2

      RETURN 
      END

C****************************************************************************
      SUBROUTINE SRDFTEXCJT(RHO,RHG, CHIVAL,NORDER,ERFEXP,
     >           Ex,d1Ex,d2Ex,Ec,d1Ec,d2Ec,VLAMBDA, ipnt)
C*****************************************************************************
C    Driver for short-range DFT exchange and correlation
C    energies and derivatives
C
C    Based on the routine SRDFTEXC of J. K. Pedersen
C
C    Created : 15-02-05, J. Toulouse
C    Modified: 17-08-09, J. Toulouse, add spin density
C    Modified: 15-08-16, E. Hedegård, add spin density gradient
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "dftinf.h"
#include "dftcom.h"

      double precision rho(4),rhg(3,4),CHIVAL,VLAMBDA ! input
      double precision rhogrd(4) ! local variable for norm of RHG(1:3,*)
      integer NORDER
      logical ERFEXP(0:2)

      double precision Ex, d1Ex(6), d2Ex(10)
      double precision Ec, d1Ec(6), d2Ec(10)

C The new terms
      double precision mu, mu_new

      save icall
      data icall/0/
C     functionals
      external ESRX_PBEHSEERF,  ESRX_PBETCSERF, ESRX_PBERIERF
      external ESRX_PBEGWSERF
      external ESRC_PBETCSJERF, ESRC_PBETCSERF, ESRC_PBERIERF
      external ESRC_LYPRIERF
      external ESRC_PBEWIERF,   ESRC_PBEGWSERF, ESRC_PBELOERF     
C     spin functionals
      external ESRX_LDAERFSPIN, ESRC_LDAERFSPIN
      external ESRX_PBEGWSERFSPIN, ESRC_PBEGWSERFSPIN
      external ESRC_MULOC_GGA, ESRC_MULOD_GGA, ESRC_MULOE_GGA   
C     sr mdc func 
      external ESRC_MD_LDAERF
       
      ! for debug prints (counter for grid points)
      integer inpt 

c  Exchange  ------------------------------------------------------------

      Ex         = 0.0d0
      d1Ex(1:6)  = 0.0d0
      d2Ex(1:10) = 0.0d0

      icall = icall + 1
      IF (ERFEXP(0)) THEN
         write (lupri,*)
     &   'WARNING srdftexcjt: erfexp true in call no. ',icall
         erfexp(1:2) = .false.
!        CALL QUIT('erfgau not implemented in SRDFTEXCJT')
      END IF

      IF (DFT_SPINDNS) THEN
         do i = 1,4
            rhogrd(i) = sqrt(rhg(1,i)**2 + rhg(2,i)**2 + rhg(3,i)**2)
         end do
      ELSE
         rhogrd(1) = sqrt(rhg(1,1)**2 + rhg(2,1)**2 + rhg(3,1)**2)
         rhogrd(2) = 0.0d0
         rhogrd(3) = 0.0d0
         rhogrd(4) = 0.0d0
      END IF

!     short-range exchange PBE of HSE
      IF(DOSRX_PBEHSE) THEN
        IF(.not. DFT_SPINDNS) THEN     
           call DESR(ESRX_PBEHSEERF,rho,RHOGRD,CHIVAL,NORDER,
     >        Ex,d1Ex,d2Ex,ipnt)
        ELSE
           call quit(' srx_pbehse not implemented for DFT_SPINDNS')
        ENDIF
      ENDIF !  (DOSRX_PBEHSE)

!     short-range exchange PBE of TCS
      IF(DOSRX_PBETCS) THEN
        call DESR(ESRX_PBETCSERF,rho,RHOGRD,CHIVAL,NORDER,Ex,d1Ex,d2Ex,
     &            ipnt) 
      ENDIF

!     short-range exchange PBE of GWS ! JT 11-08-09
      IF(DOSRX_PBEGWS) THEN
        call DESR(ESRX_PBEGWSERF,rho,RHOGRD,CHIVAL,NORDER,Ex,d1Ex,d2Ex,
     &            ipnt)
      ENDIF

C  Manu 01-02-2006 begin
!     short-range exchange RI
      IF(DOSRX_PBERI) THEN
        call DESR(ESRX_PBERIERF,rho,RHOGRD,CHIVAL,NORDER,Ex,d1Ex,d2Ex,
     &            ipnt)
      ENDIF
C  Manu 01-02-2006 end

!     short-range spin-dependent exchange LDA ! JT 11-08-09
      IF(DOSRX_LDA_S) THEN
        call DESRSPIN(ESRX_LDAERFSPIN,rho(1),rho(2),
     >             CHIVAL,DFT_SPINDNS,NORDER,Ex,d1Ex,d2Ex,ipnt)
      ENDIF
      IF(DOSRX_PBEGWS_S) THEN
        call DESRSPINGGA(ESRX_PBEGWSERFSPIN,rho,RHOGRD,CHIVAL,
     >                   norder,Ex,d1Ex,d2Ex,ipnt)
       ENDIF 
C linear complement spin-dependent exchange LDA functional, Kamal SHARKAS 09-05-11
      IF(DOLAX_LDAS)  THEN
        call ELX_LDALAMSPIN (rho(1),rho(2),VLAMBDA,DFT_SPINDNS,Ex)
        ! d1Ex(2), d1Ex(3) here ...
        call DELSPIN(ESRX_LDAERFSPIN,rho(1),rho(2),VLAMBDA,DFT_SPINDNS,
     >               NORDER,d1Ex,d2Ex)
      ENDIF

C linear complement exchange PBE functional, Kamal SHARKAS 19-05-11
      IF(DOLAX_PBEGWS) THEN

        call ELAX_PBEGWS(rho,RHOGRD,VLAMBDA,Ex)
        ! d1Ex(2), d1Ex(3) here ...
        call DELA(ESRX_PBEGWSERF,rho,RHOGRD,VLAMBDA,NORDER,
     >            d1Ex,d2Ex)
      ENDIF


C  Correlation  ----------------------------------------------------------

      Ec         = 0.0d0
      d1Ec(1:6)  = 0.0d0
      d2Ec(1:10) = 0.0d0

!     short-range correlation PBE of TCS
      IF(DOSRC_PBETCS) THEN
!        call VSRCPBEERF(ec,d1ecdrho,d1ecdgrd,d1ecdgrda2,d2ecdrhodrho,
!     >                  d2ecdgrddgrd,d2ecdrhodgrd,rho,rhogrd,chival,
!     >                  norder)
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_PBETCSERF,rho,RHOGRD,CHIVAL,NORDER,
     >            Ec,d1Ec,d2Ec,ipnt)
      ENDIF

!     short-range correlation PBE of TCSJ
      IF(DOSRC_PBETCSJ) THEN
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_PBETCSJERF,rho,RHOGRD,CHIVAL,NORDER,
     >            Ec,d1Ec,d2Ec,ipnt)
      ENDIF

!     short-range correlation PBE of GWS
      IF(DOSRC_PBEGWS) THEN
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_PBEGWSERF,rho,RHOGRD,CHIVAL,NORDER,
     >            Ec,d1Ec,d2Ec,ipnt)
      ENDIF

!     rational interpolation
      IF(DOSRC_PBERI) THEN
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_PBERIERF,rho,RHOGRD,CHIVAL,NORDER,
     >             Ec,d1Ec,d2Ec,ipnt)
      ENDIF

!     rational interpolation LYP
      IF(DOSRC_LYPRI) THEN
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_LYPRIERF,rho,RHOGRD,CHIVAL,NORDER,Ec,d1Ec,d2Ec,
     &            ipnt)
      ENDIF


!     PBE local energy formula
      IF(DOSRC_PBELO) THEN
        call DESR(ESRC_PBELOERF,rho,RHOGRD,CHIVAL,NORDER,Ec,d1Ec,d2Ec,
     &            ipnt)
      ENDIF
!     weighted interpolation
      IF(DOSRC_PBEWI) THEN
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_PBEWIERF,rho,RHOGRD,CHIVAL,NORDER,Ec,d1Ec,d2Ec,
     &            ipnt)
      ENDIF

!     short-range spin-dependent correlation LDA ! JT 18-08-09
      IF(DOSRC_LDA_S) THEN
        ! d1Ec(2), d1Ec(3) here ...
        call DESRSPIN(ESRC_LDAERFSPIN,rho(1),rho(2),
     >             CHIVAL,DFT_SPINDNS,NORDER,Ec,d1Ec,d2Ec,ipnt)
      ENDIF
      IF(DOSRC_PBEGWS_S) THEN
        call DESRSPINGGA(ESRC_PBEGWSERFSPIN,rho,RHOGRD,CHIVAL,
     >                   norder,Ec,d1Ec,d2Ec,ipnt)
      ENDIF
C Kamal SHARKAS 20-05-11 beg
C linear Complement spin-dependent (NON-Scaled) correlation LDA functional
       IF(DOLANSC_LDAS) THEN
        call ELC_LDALAMSPIN (rho(1),rho(2),VLAMBDA,DFT_SPINDNS,Ec)
        ! d1Ec(2), d1Ec(3) here ...
        call DELSPINnSC(ESRC_LDAERFSPIN,rho(1),rho(2),VLAMBDA,
     >                DFT_SPINDNS,NORDER,d1Ec,d2Ec)
       ENDIF

C Complement (NON-Scaled) correlation PBE functional
      IF(DOLANSC_PBEGWS) THEN
        call ELAC_PBEGWS(rho,RHOGRD,VLAMBDA,Ec)
        ! d1Ec(3), d1Ec(3) here ...
        call DELANSC(ESRC_PBEGWSERF,rho,RHOGRD,VLAMBDA,NORDER,d1Ec,d2Ec)
      ENDIF

C Complement Scaled correlation PBE functional
      IF(DOLASC_PBEGWS) THEN

        call ELASC_PBEGWS(rho,RHOGRD,VLAMBDA,Ec)
        ! d1Ec(3), d1Ec(3) here ...
        call DELASC(ESRC_PBEGWSERF,rho,RHOGRD,VLAMBDA,NORDER,d1Ec,d2Ec)

      ENDIF
C Kamal SHARKAS 20-05-11 end

!     short-range spin-dependent correlation LDA with multideterminant reference (for OEP) ! JT 26-08-11
!     ... or for RSDHf calculations (Manu 06-12-2012)  
      IF(DOSRC_MD_LDA) THEN
!CASS        Write(LUPRI,*)'La nouvelle fonctionnelle'
        ! d1Ec(3), d1Ec(3) here ...
        call DESRSPIN(ESRC_MD_LDAERF,rho(1),rho(2),
     >             CHIVAL,DFT_SPINDNS,NORDER,Ec,d1Ec,d2Ec,ipnt)
      ENDIF
!CASS      Write(LUPRI,*)'EC=',Ec

!     Short-range LDA functional with local CHIVAL parameter ! MNP 13-09-11 
      IF(DOSRC_MULOC_GGA) THEN
C      IF (SRCMULOFAC) THEN
C         xfac = XMULFAC_READIN
C      ELSE
C         xfac = 1.0d0/4.0d0 
C      END IF
C      mu = RHOGRD(1)/rho(1)
C      mu = xfac*mu
C      mu_new = max(mu,CHIVAL)
!      print *,'xfac,SRCMULOFAC,mu_new ==>', xfac, SRCMULOFAC,mu_new
        call ESRC_MULOC_GGA(rho,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_MULOC_GGA,rho,RHOGRD,CHIVAL,NORDER,
     >             Ec, d1Ec, d2Ec,ipnt)
      ENDIF

!     Short-range LDA functional with local CHIVAL parameter ! MNP 13-09-11 
      IF(DOSRC_MULOD_GGA) THEN
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_MULOD_GGA,rho,RHOGRD,CHIVAL,NORDER,Ec,d1Ec,d2Ec,
     &            ipnt)
      ENDIF
!     Short-range LDA functional with local CHIVAL parameter ! MNP 13-09-11 
      IF(DOSRC_MULOE_GGA) THEN
        ! d1Ec(3), d1Ec(3) here ...
        call DESR(ESRC_MULOE_GGA,rho,RHOGRD,CHIVAL,NORDER,Ec,d1Ec,d2Ec,
     &            ipnt)
      ENDIF

      RETURN
      END

C**************************************************************************
      SUBROUTINE SRDFTLTRJT(KSYMOP,DTXMAT,DSXMAT,EXCMAT,
     &                      WGHTX,WGHTC,GAO,GAO1,
     &                      CHIVAL,VLAMBDA,ERFEXP,RHO,RHG, DGAO,DOGGA,
     &                      IPNT)
C**************************************************************************
C     Generate the 1-index transformed DFT Hessian
C     needed in second-order optimization for short-range DFT hybrids.
C
C     Based on the routine SRDFTLTR of J. K. Pedersen
C
C     Created : 15-02-05, J. Toulouse
C     Modified: 17-08-09, J. Toulouse, add spin density
C     Modified: 16-02-16, E. Hedegård, add spin-gradient
C*************************************************************************** 
         implicit none
#include "priunit.h"
#include "mxcent.h"
#include "pi.h"
#include "maxorb.h"
#include "inforb.h"
#include "infinp.h"
#include "nuclei.h"
#include "dftinf.h"
#include "dftcom.h"

!     input
      integer KSYMOP
      logical DOGGA,ERFEXP(0:2)
      double precision RHO(4), RHG(3,4)
      double precision CHIVAL,VLAMBDA 
      double precision rhogrd(4) ! local variable for norm of RHG(1:3,*)
      double precision DTXMAT(NBAST,NBAST),DSXMAT(NBAST,NBAST)
      double precision GAO(NBAST),GAO1(NBAST,3), DGAO(NBAST,2)
      double precision WGHTX,WGHTC
!     output
      double precision EXCMAT(NBAST,NBAST,2)

!     local
      integer I, J, ISYM, JSYM, JEND, JSTR, ISTR, IEND
      integer I_MAT, N_MAT, NORDER

      double precision ddot

      double precision ZNV(2), BR(2), B0(2), B3(3,2), BMAX
      double precision AR(2), AB(2), AX, AY, AZ, RX(2), RY(2), RZ(2)

      double precision FAC0R, FAC0S, FACR, FACS
      double precision FRR, FRZ, FZZ, FZ0
      double precision FRS, FSS, FSZ, FSZ_S, FRZ_S, FZZ_S, FZ0_S, FZZ_SS
      double precision G0, GX, GY, GZ

!     XC energies, and derivatives
      double precision Ex, Ec, d1Ex(6), d1Ec(6), d2Ex(10), d2Ec(10)
      double precision d1Exdrhot, d1Exdrhos, d1Exdgrdt, d1Exdgrds
      double precision d1Exdngrd2, d1Exdngrds2 !d1Ex(6)
      double precision d1Ecdrhot, d1Ecdrhos, d1Ecdgrdt, d1Ecdgrds !d1Ex(6)
      double precision d1Ecdngrd2, d1Ecdngrds2

      double precision d2Exdrhotdrhot, d2Exdrhotdrhos, d2Exdrhosdrhos
      double precision d2Exdrhotdgrdt, d2Exdrhosdgrdt, d2Exdgrdtdgrdt
      double precision d2Exdrhotdgrds, d2Exdrhosdgrds, d2Exdgrdtdgrds
      double precision d2Exdgrdsdgrds ! d2Ex(10)

      double precision d2Ecdrhotdrhot, d2Ecdrhotdrhos, d2Ecdrhosdrhos
      double precision d2Ecdrhotdgrdt, d2Ecdrhosdgrdt, d2Ecdgrdtdgrdt
      double precision d2Ecdrhotdgrds, d2Ecdrhosdgrds, d2Ecdgrdtdgrds
      double precision d2Ecdgrdsdgrds ! d2Ec(10)

      ! edh for debug
      integer ipnt

      if (ipnt .eq. 1) write(lupri,*)'entering SRDFTLTRJT'
      ! For small values of rho to avoid NaN
      do i = 1,4
         if (abs(rho(i)) .lt. 1.d-30) rho(i) = 1.d-30
      end do

      if (ipnt .eq. 1) then 
         write(lupri,*)'DGAO(1,1) and DGAI(1,2)'
         call output(DGAO(1,1),1,nbast,1,nbast,nbast,nbast,-1,lupri)
         call output(DGAO(1,2),1,nbast,1,nbast,nbast,nbast,-1,lupri)
      end if 


      ! DGAO(1,1) = DTXMAT*GAO  
      CALL DGEMV('N',NBAST,NBAST,1.0d0,
     &            DTXMAT,NBAST,GAO,1,0.0d0,DGAO(1,1),1)
      IF (DFT_SPINDNS) THEN
        N_MAT = 2
        ! DGAO(1,2) = DSXMAT*GAO
        CALL DGEMV('N',NBAST,NBAST,1.0d0,DSXMAT,
     &              NBAST,GAO,1,0.0d0,DGAO(1,2),1)
      ELSE
        N_MAT = 1
      END IF
C
      BMAX = 0.0D0
      DO I_MAT = 1, N_MAT ! 1->total, 2->spin
         ! rho and rho_s
         B0(I_MAT) = DDOT(NBAST,DGAO(1,I_MAT),1,GAO,1)
         BMAX = MAX( BMAX,ABS(B0(I_MAT)) )
         IF (DOGGA) THEN
            rhogrd(I_MAT) =
     &       sqrt( rhg(1,I_MAT)**2 + rhg(2,I_MAT)**2 + rhg(3,I_MAT)**2 )
C           B3(1) = GAO1'*DGAO(1,1); B3(2) = GAO1'*DGAO(1,2)
            CALL DGEMV('T',NBAST,3,1.0d0,GAO1,NBAST,
     &         DGAO(1,I_MAT),1,0.0d0,B3(1,I_MAT),1)
C           DGAO(1,1)= DTXMAT'*GAO; DGAO(1,2) = DSXMAT'*GAO !!! HJAAJ: is this necessary ? We have D*XMAT symmetrized!
            IF (I_MAT .EQ. 1) THEN  
               CALL DGEMV('T',NBAST,NBAST,1.0d0,DTXMAT,NBAST,
     &                    GAO,1,0.0d0,DGAO(1,I_MAT),1)
            ELSE IF (I_MAT .EQ. 2) THEN 
               CALL DGEMV('T',NBAST,NBAST,1.0d0,DSXMAT,NBAST,
     &                    GAO,1,0.0d0,DGAO(1,I_MAT),1)
            END IF 
C           B3(1,I_MAT) = B3(1,I_MAT) + GAO1'*DGAO(1, I_MAT)
            CALL DGEMV('T',NBAST,3,1.0d0,GAO1,NBAST,
     &         DGAO(1,I_MAT),1,1.0d0,B3(1,I_MAT),1)
            BMAX = MAX(BMAX,
     &       ABS(B3(1,I_MAT)),ABS(B3(2,I_MAT)),ABS(B3(3,I_MAT)))
         END IF
      END DO
      !edh debug
      !write(lupri,*)'test B3(1,1); B3(1,2)', B3(1,1), B3(1,2)
C
      IF (BMAX.GT.DFTHRL) THEN 
C

C Exchange and correlation derivatives ------------------------------------

      NORDER = 2  ! first-order and second-order derivatives
      CALL SRDFTEXCJT(RHO,RHG,CHIVAL,NORDER,ERFEXP,
     >                Ex,d1Ex,d2Ex,Ec,d1Ec,d2Ec,VLAMBDA,ipnt)
       if (ipnt .eq. 1) write(lupri,*)'d1Ex', d1Ex
       if (ipnt .eq. 1) write(lupri,*)'d2Ex', d2Ex
       if (DOGGA .and. DFT_SPINDNS) then
           if (ipnt .eq. 1) write(lupri,*)'case 1 active'
           d1Exdrhot      = d1Ex(1)
           d1Exdgrdt      = d1Ex(2)
           d1Exdngrd2     = d1Ex(3)
           d1Exdrhos      = d1Ex(4)
           d1Exdgrds      = d1Ex(5)
           d1Exdngrds2    = d1Ex(6)
           !-----------------------
           d1Ecdrhot      = d1Ec(1)
           d1Ecdgrdt      = d1Ec(2)
           d1Ecdngrd2     = d1Ec(3)
           d1Ecdrhos      = d1Ec(4)
           d1Ecdgrds      = d1Ec(5)
           d1Ecdngrds2    = d1Ec(6)
           !-----------------------
           d2Exdrhotdrhot = d2Ex(1)
           d2Exdrhotdrhos = d2Ex(2)
           d2Exdrhosdrhos = d2Ex(3)
           d2Exdrhotdgrdt = d2Ex(4)
           d2Exdrhosdgrdt = d2Ex(5)
           d2Exdgrdtdgrdt = d2Ex(6)
           d2Exdrhotdgrds = d2Ex(7)
           d2Exdrhosdgrds = d2Ex(8)
           d2Exdgrdtdgrds = d2Ex(9)
           d2Exdgrdsdgrds = d2Ex(10)
           !-----------------------
           d2Ecdrhotdrhot = d2Ec(1)
           d2Ecdrhotdrhos = d2Ec(2)
           d2Ecdrhosdrhos = d2Ec(3)
           d2Ecdrhotdgrdt = d2Ec(4)
           d2Ecdrhosdgrdt = d2Ec(5)
           d2Ecdgrdtdgrdt = d2Ec(6)
           d2Ecdrhotdgrds = d2Ec(7)
           d2Ecdrhosdgrds = d2Ec(8)
           d2Ecdgrdtdgrds = d2Ec(9)
           d2Ecdgrdsdgrds = d2Ec(10)
       else if (DFT_SPINDNS .and. .not. DOGGA) then
           if (ipnt .eq. 1) write(lupri,*)'case 2 active'
           d1Exdrhot      = d1Ex(1)
           d1Exdrhos      = d1Ex(2)
           !----------------------
           d1Ecdrhot      = d1Ec(1)
           d1Ecdrhos      = d1Ec(2)
           !----------------------
           d2Exdrhotdrhot = d2Ex(1)
           d2Exdrhotdrhos = d2Ex(2)
           d2Exdrhosdrhos = d2Ex(3)
           !----------------------
           d2Ecdrhotdrhot = d2Ec(1)
           d2Ecdrhotdrhos = d2Ec(2)
           d2Ecdrhosdrhos = d2Ec(3)
       else if (DOGGA .and. .not. DFT_SPINDNS) then
           if (ipnt .eq. 1) write(lupri,*)'case 3 active'
           d1Exdrhot      = d1Ex(1)
           d1Exdgrdt      = d1Ex(2)
           d1Exdngrd2     = d1Ex(3)
           !----------------------
           d1Ecdrhot      = d1Ec(1) 
           d1Ecdgrdt      = d1Ec(2)
           d1Ecdngrd2     = d1Ec(3)
           !----------------------
           d2Exdrhotdrhot = d2Ex(1) 
           d2Exdgrdtdgrdt = d2Ex(2)
           d2Exdrhotdgrdt = d2Ex(3)
           !----------------------
           d2Ecdrhotdrhot = d2Ec(1)
           d2Ecdgrdtdgrdt = d2Ec(2)
           d2Ecdrhotdgrdt = d2Ec(3)
       else ! (.not. DOGGA .and. .not. DFT_SPINDNS)
           if (ipnt .eq. 1) write(lupri,*)'case 4 active'
           d1Exdrhot      = d1Ex(1)
           d1Ecdrhot      = d1Ec(1)
           d2Exdrhotdrhot = d2Ex(1) 
           d2Ecdrhotdrhot = d2Ec(1)
       end if 



C Total contribution -------------------------------------------------------

      FRR = WGHTX*d2Exdrhotdrhot + WGHTC*d2Ecdrhotdrhot
      IF (DOGGA .AND. DFT_SPINDNS) THEN
         ZNV(1) = 1.0d0/RHOGRD(1)
         IF (RHOGRD(2) > 1.0D-30) THEN
             ZNV(2) = 1.0d0/RHOGRD(2)
         ELSE 
           ZNV(2) = 0.0D0
         END IF 
         FZ0   = ZNV(1)*(WGHTX*d1Exdgrdt + WGHTC*d1Ecdgrdt)
         ! gradient
         FRZ   = WGHTX*d2Exdrhotdgrdt + WGHTC*d2Ecdrhotdgrdt
         FZZ   = WGHTX*d2Exdgrdtdgrdt + WGHTC*d2Ecdgrdtdgrdt
         ! spin
         FZ0_S = ZNV(2)*(WGHTX*d1Exdgrds + WGHTC*d1Ecdgrds)
         FSZ   = WGHTX*d2Exdrhosdgrdt + WGHTC*d2Ecdrhosdgrdt
         ! spin-gradient
         FSZ_S = WGHTX*d2Exdrhosdgrds + WGHTC*d2Ecdrhosdgrds
         FRZ_S = WGHTX*d2Exdrhotdgrds + WGHTC*d2Ecdrhotdgrds
         FZZ_S = WGHTX*d2Exdgrdtdgrds + WGHTC*d2Ecdgrdtdgrds
         FZZ_SS= WGHTX*d2Exdgrdsdgrds + WGHTC*d2Ecdgrdsdgrds

         ! normalized grd[rhot^(X,Y,Z)]
         RX(1) = ZNV(1)*RHG(1,1)
         RY(1) = ZNV(1)*RHG(2,1)
         RZ(1) = ZNV(1)*RHG(3,1)
         ! normalized grd[rhos^(X,Y,Z)]
         RX(2) = ZNV(2)*RHG(1,2)
         RY(2) = ZNV(2)*RHG(2,2)
         RZ(2) = ZNV(2)*RHG(3,2)
         ! pure spin --- SPIN + GGA NOT FINISHED YET!!!
         FSS = WGHTX*d2Exdrhosdrhos + WGHTC*d2Ecdrhosdrhos
         FRS = WGHTX*d2Exdrhotdrhos + WGHTC*d2Ecdrhotdrhos
      ELSE IF (DFT_SPINDNS .AND. .NOT. DOGGA) THEN
         FSS = WGHTX*d2Exdrhosdrhos + WGHTC*d2Ecdrhosdrhos
         FRS = WGHTX*d2Exdrhotdrhos + WGHTC*d2Ecdrhotdrhos
      ELSE IF (DOGGA .AND. .NOT. DFT_SPINDNS) THEN
         ZNV(1) = 1.0d0/RHOGRD(1)
         FZ0 = ZNV(1)*(WGHTX*d1Exdgrdt + WGHTC*d1Ecdgrdt)
         FRZ = WGHTX*d2Exdrhotdgrdt + WGHTC*d2Ecdrhotdgrdt
         FZZ = WGHTX*d2Exdgrdtdgrdt + WGHTC*d2Ecdgrdtdgrdt
         ! normalized grd[rhot^(X,Y,Z)]
         RX(1) = ZNV(1)*RHG(1,1)
         RY(1) = ZNV(1)*RHG(2,1)
         RZ(1) = ZNV(1)*RHG(3,1)
      ELSE 
         call quit('SRDFTLTRJT: we should not end here')
      END IF

!#ifdef SRDFT_DEBUG
          if (ipnt .eq. 1) then
           !if (abs(FSS) .gt. 0.d0 .or. abs(FRS) .gt. 0.d0) then
              write(lupri,*) 'hjaaj: rho',rho(1:4)
              write(lupri,*) 'hjaaj: ZNV',ZNV(1:2)
              write(lupri,*) 'hjaaj: B0',B0(1:2)
              write(lupri,*) 'hjaaj: Ex, Ec',Ex, Ec
              write(lupri,*) 'hjaaj: d1Exdrhot  ',d1Exdrhot
              write(lupri,*) 'hjaaj: d1Exdgrdt  ',d1Exdgrdt
              write(lupri,*) 'hjaaj: d1Exdngrd2 ',d1Exdngrd2
              write(lupri,*) 'hjaaj: d1Exdrhos  ',d1Exdrhos
              write(lupri,*) 'hjaaj: d1Exdgrds  ',d1Exdgrds
              write(lupri,*) 'hjaaj: d1Exdngrds2',d1Exdngrds2

              write(lupri,*) 'hjaaj: d1Ecdrhot  ',d1Ecdrhot
              write(lupri,*) 'hjaaj: d1Ecdgrdt  ',d1Ecdgrdt
              write(lupri,*) 'hjaaj: d1Ecdngrd2 ',d1Ecdngrd2
              write(lupri,*) 'hjaaj: d1Ecdrhos  ',d1Ecdrhos
              write(lupri,*) 'hjaaj: d1Ecdgrds  ',d1Ecdgrds
              write(lupri,*) 'hjaaj: d1Ecdngrds2',d1Ecdngrds2

              write(lupri,*) 'hjaaj: d2Exdrhotdrhot',d2Exdrhotdrhot
              write(lupri,*) 'hjaaj: d2Exdrhotdrhos',d2Exdrhotdrhos
              write(lupri,*) 'hjaaj: d2Exdrhosdrhos',d2Exdrhosdrhos
              write(lupri,*) 'hjaaj: d2Exdrhotdgrdt',d2Exdrhotdgrdt
              write(lupri,*) 'hjaaj: d2Exdrhosdgrdt',d2Exdrhosdgrdt
              write(lupri,*) 'hjaaj: d2Exdgrdtdgrdt',d2Exdgrdtdgrdt
              write(lupri,*) 'hjaaj: d2Exdrhotdgrds',d2Exdrhotdgrds
              write(lupri,*) 'hjaaj: d2Exdrhosdgrds',d2Exdrhosdgrds
              write(lupri,*) 'hjaaj: d2Exdgrdtdgrds',d2Exdgrdtdgrds
              write(lupri,*) 'hjaaj: d2Exdgrdsdgrds',d2Exdgrdsdgrds

              write(lupri,*) 'hjaaj: d2Ecdrhotdrhot',d2Ecdrhotdrhot
              write(lupri,*) 'hjaaj: d2Ecdrhotdrhos',d2Ecdrhotdrhos
              write(lupri,*) 'hjaaj: d2Ecdrhosdrhos',d2Ecdrhosdrhos
              write(lupri,*) 'hjaaj: d2Ecdrhotdgrdt',d2Ecdrhotdgrdt
              write(lupri,*) 'hjaaj: d2Ecdrhosdgrdt',d2Ecdrhosdgrdt
              write(lupri,*) 'hjaaj: d2Ecdgrdtdgrdt',d2Ecdgrdtdgrdt
              write(lupri,*) 'hjaaj: d2Ecdrhotdgrds',d2Ecdrhotdgrds
              write(lupri,*) 'hjaaj: d2Ecdrhosdgrds',d2Ecdrhosdgrds
              write(lupri,*) 'hjaaj: d2Ecdgrdtdgrds',d2Ecdgrdtdgrds
              write(lupri,*) 'hjaaj: d2Ecdgrdsdgrds',d2Ecdgrdsdgrds

              write(lupri,*) 'FZ0', FZ0
              write(lupri,*) 'hjaaj: FRR, FRZ, FZZ,',FRR, FRZ, FZZ
              write(lupri,*) 'hjaaj: FZ0_S',FZ0_S
              write(lupri,*) 'hjaaj: FSZ_S,FRZ_S,FZZ_S,FZZ_SS',
     &                               FSZ_S,FRZ_S,FZZ_S,FZZ_SS
           !end if
          end if 
!#endif 


C Linear transformation -----------------------------------------------------

!        IF (DOHFSRDFT) THEN   !hjaaj July 2010: this was a
!        misunderstanding of the code in SRDFTLTR which had  IF
!        (.NOT.SRHYBR .OR. DHFSRDFT) ; for normal MC-srDFT SRHYBR would
!        have been false, but in this code this part was skipped until
!        now if DOMCSRDFT

            ! It is used that FSR = FRS, FSZ = FZS, FSZ_S = FZS_S
            if (dogga .and. dft_spindns) then
               BR(1) = B3(1,1)*RX(1) + B3(2,1)*RY(1) + B3(3,1)*RZ(1)
               BR(2) = B3(1,2)*RX(2) + B3(2,2)*RY(2) + B3(3,2)*RZ(2)
               FAC0R = FRR*B0(1) + FRZ*BR(1) + FRS*B0(2) + FRZ_S*BR(2)
               FAC0S = FRS*B0(1) + FSS*B0(2) + FSZ*BR(1) + FSZ_S*BR(2)
               FACR  = FRZ*B0(1) + FZZ*BR(1) + FZZ_S*BR(2) +  FSZ*B0(2)
               FACS  = FSZ_S*B0(2) + FRZ_S*B0(1) +
     &                 FZZ_S*BR(1) + FZZ_SS*BR(2)
            else if (dogga .and. .not. dft_spindns) then
               BR(1) = B3(1,1)*RX(1) + B3(2,1)*RY(1) + B3(3,1)*RZ(1)
               FAC0R = FRR*B0(1) + FRZ*BR(1)
               FACR = FRZ*B0(1) + FZZ*BR(1)
            else if (dft_spindns .and. .not. dogga) then
               FAC0R = FRR*B0(1) + FRS*B0(2)
               FAC0S = FRS*B0(1) + FSS*B0(2) ! as FSR = FRS
            else
               FAC0R = FRR*B0(1)
            end if
            if (ipnt .eq. 1) then
               write(lupri,*) 'hjaaj: BR', BR(1:2)
               write(lupri,*) 'hjaaj: FAC0R,FAC0S',FAC0R,FAC0S
            end if
            IF (NSYM.EQ.1) THEN
               DO I = 1, NBAST
                  G0 = GAO(I)
                  IF (DOGGA) THEN
                     GX = GAO1(I,1)
                     GY = GAO1(I,2)
                     GZ = GAO1(I,3)
                  END IF
                  DO J = 1, I
                     EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
                     if (dogga .and. dft_spindns) then
                        EXCMAT(J,I,2) = EXCMAT(J,I,2)+FAC0S*G0*GAO(J)

                        AX = GX*GAO(J) + G0*GAO1(J,1)
                        AY = GY*GAO(J) + G0*GAO1(J,2)
                        AZ = GZ*GAO(J) + G0*GAO1(J,3)
                        AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                        AR(2) = AX*RX(2) + AY*RY(2) + AZ*RZ(2)
                        AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
     &                          - AR(1)*BR(1)
                        AB(2) = AX*B3(1,2) + AY*B3(2,2) + AZ*B3(3,2)
     &                          - AR(2)*BR(2)
                        ! grad + spin-grad
                        EXCMAT(J,I,1) = EXCMAT(J,I,1)
     &                                + FACR*AR(1) + FZ0*AB(1)
                        EXCMAT(J,I,2) = EXCMAT(J,I,2)
     &                                + FACS*AR(2) + FZ0_S*AB(2)
                     else if (dft_spindns .and. .not. dogga) then
                        EXCMAT(J,I,2) = EXCMAT(J,I,2)+FAC0S*G0*GAO(J)
                     else if (dogga .and. .not. dft_spindns) then
                           AX = GX*GAO(J) + G0*GAO1(J,1)
                           AY = GY*GAO(J) + G0*GAO1(J,2)
                           AZ = GZ*GAO(J) + G0*GAO1(J,3)
                           AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                           AB(1) = AX*B3(1,1) + AY*B3(2,1) + AZ*B3(3,1)
     &                        - AR(1)*BR(1)
                            EXCMAT(J,I,1) = EXCMAT(J,I,1)
     &                                    + FACR*AR(1)+FZ0*AB(1)
                     end if
                  END DO
               END DO
            ELSE
               DO ISYM = 1, NSYM
                  ISTR = IBAS(ISYM) + 1
                  IEND = IBAS(ISYM) + NBAS(ISYM)
                  JSYM = MULD2H(ISYM,KSYMOP)
                  IF (ISYM.GE.JSYM) THEN
                     JSTR = IBAS(JSYM) + 1
                     JEND = IBAS(JSYM) + NBAS(JSYM)
                     DO I = ISTR, IEND
                        G0 = GAO(I)
                        IF (DOGGA) THEN
                           GX = GAO1(I,1)
                           GY = GAO1(I,2)
                           GZ = GAO1(I,3)
                        END IF
                        DO J = JSTR, MIN(I,JEND)
                           EXCMAT(J,I,1) = EXCMAT(J,I,1)+FAC0R*G0*GAO(J)
!                          IF (DFT_SPINDNS) THEN
!                             EXCMAT(J,I,2) = EXCMAT(J,I,2)
!    &                                      + FAC0S*G0*GAO(J)
!                          END IF
!                          IF (DOGGA) THEN
!                             AX = GX*GAO(J) + G0*GAO1(J,1)
!                             AY = GY*GAO(J) + G0*GAO1(J,2)
!                             AZ = GZ*GAO(J) + G0*GAO1(J,3)
!                             AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
!                             AB(1) = AX*B3(1,1) +
!    &                        AY*B3(2,1)+ AZ*B3(3,1) - AR(1)*BR(1)
!                             EXCMAT(J,I,1) = EXCMAT(J,I,1)
!    &                        + FACR*AR(1) + FZ0*AB(1)
!                          END IF
                           if (dogga .and. dft_spindns) then
                              ! spin
                              EXCMAT(J,I,2) = EXCMAT(J,I,2) +
     &                                        FAC0S*G0*GAO(J)
                              ! grad + spin-grad
                              AX = GX*GAO(J) + G0*GAO1(J,1)
                              AY = GY*GAO(J) + G0*GAO1(J,2)
                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
                              AR(1) = AX*RX(1)  + AY*RY(1) + AZ*RZ(1)
                              AR(2) = AX*RX(2)  + AY*RY(2) + AZ*RZ(2)
                              AB(1) = AX*B3(1,1) + AY*B3(2,1)
     &                              + AZ*B3(3,1) - AR(1)*BR(1)
                              AB(2) = AX*B3(1,2) + AY*B3(2,2)
     &                              + AZ*B3(3,2) - AR(2)*BR(2)
                              EXCMAT(J,I,1) = EXCMAT(J,I,1)
     &                                      + FACR*AR(1)+FZ0*AB(1)
                              EXCMAT(J,I,2) = EXCMAT(J,I,2)
     &                                      + FACS*AR(2) + FZ0_S*AB(2)
                           else if (dft_spindns .and. .not. dogga) then
                              EXCMAT(J,I,2) = EXCMAT(J,I,2)
     &                                      + FAC0S*G0*GAO(J)
                           else if (dogga) then
                              AX = GX*GAO(J) + G0*GAO1(J,1)
                              AY = GY*GAO(J) + G0*GAO1(J,2)
                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
                              AR(1) = AX*RX(1) + AY*RY(1) + AZ*RZ(1)
                              AB(1) = AX*B3(1,1) + AY*B3(2,1)
     &                              + AZ*B3(3,1) - AR(1)*BR(1)
                              EXCMAT(J,I,1) = EXCMAT(J,I,1)
     &                                      + FACR*AR(1)+FZ0*AB(1)
                            end if
                        END DO
                     END DO
                  END IF
               END DO
            END IF
!        END IF   ! IF (DOHFSRDFT) THEN
      END IF   ! IF (BMAX.GT.DFTHRL) THEN
c
      RETURN
      END
