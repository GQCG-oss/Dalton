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
      SUBROUTINE GETRHO_OLD(DMAT,GSO,RHO,DMATGAO,DFTHRI,IPRINT)
C*****************************************************************************
C
C     T. Helgaker feb 01
C     (RHO13 removed Aug 18)
C
C Output:
C    RHO
C    DMATGAO(i) = sum(j) dmat(i,j) * gso(j)
C
C*****************************************************************************
      implicit none
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "symmet.h"
#include "inforb.h"
      real*8, parameter :: D0 = 0.0d0
      real*8, parameter :: D1 = 1.0d0
      real*8, parameter :: D2 = 2.0d0
      real*8, parameter :: DP5 = 0.5d0
      real*8, parameter :: DP3 = D1/3.0d0
C
      integer, intent(in) :: IPRINT
      real*8, intent(in) :: DMAT(NBAST,NBAST), GSO(NBAST), DFTHRI
      real*8, intent(out) :: RHO, DMATGAO(NBAST)
C
      integer :: ISTR, ISYM, NBASI
      real*8 :: DDOT
C
      IF (IPRINT .GT. 100) THEN
         write(lupri,*) 'GETRHO_OLD NBAST, NSYM',NBAST,NSYM
         write(lupri,*) 'DFTHRI, DMAT(1,1)',DFTHRI,DMAT(1,1),GSO(1)
         write(lupri,*) 'IBAS',ibas(1:nsym)
         write(lupri,*) 'NBAS',nbas(1:nsym)
      END IF
      IF (NSYM.EQ.1) THEN
         CALL DSYMV('U',NBAST,D1,DMAT,NBAST,GSO,1,D0,DMATGAO,1)
      ELSE
         CALL DZERO(DMATGAO,NBAST)
         DO ISYM = 1, NSYM
            ISTR = IBAS(ISYM) + 1
            NBASI= NBAS(ISYM)
            CALL DSYMV('U',NBASI,D1,DMAT(ISTR,ISTR),NBAST,GSO(ISTR),1,
     &                 D0,DMATGAO(ISTR),1)
         END DO
      END IF
      RHO = DDOT(NBAST,GSO,1,DMATGAO,1)
      RETURN
      END    


C*****************************************************************************
      SUBROUTINE DFTKSMGGASPIN(EXCMAT, GAO, GAO1, RHG, d1E, DFTHRL)
C*****************************************************************************
C
C     Erik Hedegaard (feb. 2016) based on DFTKSM by T. Helgaker
C     Revised Apr. 2018 by Hans Joergen Aa. Jensen
C
C     RHG(1:3,1:4): grad(rhoc), grad(rhos), grad(rhoa), grad(rhob)
C     VXC  = d1E(1) = d(e_xc) / d(rhoc)
C     VSC  = d1E(2) = d(e_xc) / d(rhos)
C     VXB  = d1E(3) = d(e_xc) / d(grdcc)
C     VSB  = d1E(4) = d(e_xc) / d(grdss)
C     VXSB = d1E(5) = d(e_xc) / d(grdcs)
C
C*****************************************************************************
      implicit none
! inforb.h : NBAST, NSYM
#include "inforb.h"
      real*8, intent(in) :: GAO(NBAST), GAO1(NBAST,3), RHG(3,2)
      real*8, intent(in) :: d1E(5), DFTHRL
      real*8, intent(inout) :: EXCMAT(NBAST,NBAST,2)

      ! local variables:
      real*8 :: VXB, VXC, VSB, VSC, VXSB
      real*8 :: FX, FY, FZ, SX, SY, SZ
      real*8 :: FC, FCS
      integer :: I, J, ISYM, ISTR, IEND
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      VXC  = d1E(1)
      VSC  = d1E(2)
      VXB  = d1E(3) * 4.0d0
      VSB  = d1E(4) * 4.0d0
      VXSB = d1E(5) * 2.0d0
      FX = VXB*RHG(1,1) + VXSB*RHG(1,2)
      FY = VXB*RHG(2,1) + VXSB*RHG(2,2)
      FZ = VXB*RHG(3,1) + VXSB*RHG(3,2)
      SX = VSB*RHG(1,2) + VXSB*RHG(1,1)
      SY = VSB*RHG(2,2) + VXSB*RHG(2,1)
      SZ = VSB*RHG(3,2) + VXSB*RHG(3,1)
      IF (NSYM.EQ.1) THEN
         DO J = 1, NBAST
            FC = VXC*GAO(J)+FX*GAO1(J,1)+FY*GAO1(J,2)+FZ*GAO1(J,3)
            IF (abs(FC).GT.DFTHRL) THEN
               CALL DAXPY(NBAST,FC,GAO,1,EXCMAT(1,J,1),1)
            END IF
            FCS = VSC*GAO(J)+SX*GAO1(J,1)+SY*GAO1(J,2)+SZ*GAO1(J,3)
            IF (abs(FCS).GT.DFTHRL) THEN
               CALL DAXPY(NBAST,FCS,GAO,1,EXCMAT(1,J,2),1)
            END IF
         END DO
      ELSE
         DO ISYM = 1, NSYM
            ISTR = IBAS(ISYM) + 1
            IEND = IBAS(ISYM) + NBAS(ISYM)
            DO J = ISTR, IEND
               FC = VXC*GAO(J)+FX*GAO1(J,1)+FY*GAO1(J,2)+FZ*GAO1(J,3)
               IF (abs(FC).GT.DFTHRL) THEN
                  DO I = ISTR, IEND
                     EXCMAT(I,J,1) = EXCMAT(I,J,1) + FC*GAO(I)
                  END DO
               END IF
               FCS = VSC*GAO(J)+SX*GAO1(J,1)+SY*GAO1(J,2)+SZ*GAO1(J,3)
               IF (abs(FCS).GT.DFTHRL) THEN
                  DO I = ISTR, IEND
                     EXCMAT(I,J,2) = EXCMAT(I,J,2) + FCS*GAO(I)
                  END DO
               END IF
            END DO
         END DO
      END IF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTKSM_MGGA(EXCMAT, GAO, GAO1, GAO2, RHG, d1E,
     &                           DFTHRL, DO_SPIN)
C*****************************************************************************
C
C     Erik Kjellgren (oct. 2018) based on DFTKSMGGASPIN by E. Hedegaard
C
C     RHG(1:3,1:4): grad(rhoc), grad(rhos), grad(rhoa), grad(rhob)
C     VXC  = d1E(1) = d(e_xc) / d(rhoc)
C     VSC  = d1E(2) = d(e_xc) / d(rhos)
C     VXB  = d1E(3) = d(e_xc) / d(grdcc)
C     VSB  = d1E(4) = d(e_xc) / d(grdss)
C     VXSB = d1E(5) = d(e_xc) / d(grdcs)
C     VXT  = d1E(6) = d(e_xc) / d(tauc)
C     VST  = d1E(7) = d(e_xc) / d(taus)
C     VXL  = d1E(8) = d(e_xc) / d(laplace(rhoc))
C     VSL  = d1E(9) = d(e_xc) / d(laplace(rhos))
C
C     GAO2 is just a dummy variable for now, since it is not calculated
C     in the code. It is put in now to make it easier to implement MGGA
C     functional in the future that depends on laplace(rho). GAO2 is
C     supposed to be defined as:
C     GAO2 = laplace(Omega) = lapalce(GAOp)*GAOq + GAOp*laplace(GAOq)
C                              + 2*GRAD(GAOp)*GRAD(GAOq)
C     When calling this subroutine just set GAO2 = 0.0d0
C     Later remember to change it to a vector when actually need in the
C     declaration.
C
C*****************************************************************************
      implicit none
! inforb.h : NBAST, NSYM
#include "inforb.h"
      logical, intent(in) :: DO_SPIN
      real*8, intent(in) :: GAO(NBAST), GAO1(NBAST,3), RHG(3,2)
      real*8, intent(in) :: GAO2, d1E(9), DFTHRL
      real*8, intent(inout) :: EXCMAT(NBAST,NBAST,2)

      ! local variables:
      real*8 :: VXB, VXC, VSB, VSC, VXSB, VXT, VST
      !real*8 :: VXL, VSL
      real*8 :: FX, FY, FZ, SX, SY, SZ
      real*8 :: FC, FCS
      integer :: I, J, ISYM, ISTR, IEND
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      VXC  = d1E(1)
      VSC  = d1E(2)
      VXB  = d1E(3) * 4.0d0
      VSB  = d1E(4) * 4.0d0
      VXSB = d1E(5) * 2.0d0
      VXT  = d1E(6) * 0.5d0
      VST  = d1E(7) * 0.5d0
      !VXL  = d1E(8), no laplace(rho) implementation yet
      !VSL  = d1E(9), no laplace(rho) implementation yet
      FX = VXB*RHG(1,1) + VXSB*RHG(1,2)
      FY = VXB*RHG(2,1) + VXSB*RHG(2,2)
      FZ = VXB*RHG(3,1) + VXSB*RHG(3,2)
      SX = VSB*RHG(1,2) + VXSB*RHG(1,1)
      SY = VSB*RHG(2,2) + VXSB*RHG(2,1)
      SZ = VSB*RHG(3,2) + VXSB*RHG(3,1)
      IF (NSYM.EQ.1) THEN
         DO J = 1, NBAST
            FC = VXC*GAO(J)+FX*GAO1(J,1)+FY*GAO1(J,2)+FZ*GAO1(J,3)
            IF (abs(FC).GT.DFTHRL) THEN
               CALL DAXPY(NBAST,FC,GAO,1,EXCMAT(1,J,1),1)
            END IF
            ! Meta-GGA tau part of FC
            DO I = 1, NBAST
               EXCMAT(I,J,1) = EXCMAT(I,J,1) 
     &                         + VXT * (GAO1(J,1)*GAO1(I,1)
     &                                 + GAO1(J,2)*GAO1(I,2)
     &                                 + GAO1(J,3)*GAO1(I,3))
C     &                        + VXL * (GAO2(J)*GAO(I)
C     &                               + GAO(J)*GAO2(I)
C     &                               + 2.0d0*GAO1(J,1)*GAO1(I,1)
C     &                               + 2.0d0*GAO1(J,2)*GAO1(I,2)
C     &                               + 2.0d0*GAO1(J,3)*GAO1(I,3)
            END DO
            IF (DO_SPIN) THEN
               FCS = VSC*GAO(J)+SX*GAO1(J,1)+SY*GAO1(J,2)+SZ*GAO1(J,3)
               IF (abs(FCS).GT.DFTHRL) THEN
                  CALL DAXPY(NBAST,FCS,GAO,1,EXCMAT(1,J,2),1)
               END IF
               ! Meta-GGA tau part of FCS
               DO I = 1, NBAST
                  EXCMAT(I,J,2) = EXCMAT(I,J,2) 
     &                           + VST * (GAO1(J,1)*GAO1(I,1)
     &                                   + GAO1(J,2)*GAO1(I,2)
     &                                   + GAO1(J,3)*GAO1(I,3))
C     &                           + VSL * (GAO2(J)*GAO(I)
C     &                                  + GAO(J)*GAO2(I)
C     &                                  + 2.0d0*GAO1(J,1)*GAO1(I,1)
C     &                                  + 2.0d0*GAO1(J,2)*GAO1(I,2)
C     &                                  + 2.0d0*GAO1(J,3)*GAO1(I,3)
              END DO
            END IF
         END DO
      ELSE
         DO ISYM = 1, NSYM
            ISTR = IBAS(ISYM) + 1
            IEND = IBAS(ISYM) + NBAS(ISYM)
            DO J = ISTR, IEND
               FC = VXC*GAO(J)+FX*GAO1(J,1)+FY*GAO1(J,2)+FZ*GAO1(J,3)
               IF (abs(FC).GT.DFTHRL) THEN
                  DO I = ISTR, IEND
                     EXCMAT(I,J,1) = EXCMAT(I,J,1) + FC*GAO(I)
                  END DO
               END IF
               DO I = ISTR, IEND
                  EXCMAT(I,J,1) = EXCMAT(I,J,1)
     &                            + VXT * (GAO1(J,1)*GAO1(I,1)
     &                                    + GAO1(J,2)*GAO1(I,2)
     &                                    + GAO1(J,3)*GAO1(I,3))
C     &                           + VXL * (GAO2(J)*GAO(I)
C     &                                  + GAO(J)*GAO2(I)
C     &                                  + 2.0d0*GAO1(J,1)*GAO1(I,1)
C     &                                  + 2.0d0*GAO1(J,2)*GAO1(I,2)
C     &                                  + 2.0d0*GAO1(J,3)*GAO1(I,3)
               END DO
               IF (DO_SPIN) THEN
                 FCS = VSC*GAO(J)+SX*GAO1(J,1)+SY*GAO1(J,2)+SZ*GAO1(J,3)
                 IF (abs(FCS).GT.DFTHRL) THEN
                    DO I = ISTR, IEND
                       EXCMAT(I,J,2) = EXCMAT(I,J,2) + FCS*GAO(I)
                    END DO
                 END IF
                 DO I = ISTR, IEND
                   EXCMAT(I,J,2) = EXCMAT(I,J,2)
     &                             + VST * (GAO1(J,1)*GAO1(I,1)
     &                                     + GAO1(J,2)*GAO1(I,2)
     &                                     + GAO1(J,3)*GAO1(I,3))
C     &                            + VSL * (GAO2(J)*GAO(I)
C     &                                   + GAO(J)*GAO2(I)
C     &                                   + 2.0d0*GAO1(J,1)*GAO1(I,1)
C     &                                   + 2.0d0*GAO1(J,2)*GAO1(I,2)
C     &                                   + 2.0d0*GAO1(J,3)*GAO1(I,3)
                 END DO
               END IF
            END DO
         END DO
      END IF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTKSM(EXCMAT,GAO,GAO1,RHG,VXC,VXB,DOGGA,FROMVX,DFTHRL)
C*****************************************************************************
C
C     T. Helgaker oct 2000
C
C     VXC = d(e_xc) / d(rhoc)
C     VXB = d(e_xc) / d(grdcc)
C
C*****************************************************************************
      implicit none
#include "inforb.h"
      real*8, parameter :: D2 = 2.0D0
C
      logical, intent(in) :: DOGGA, FROMVX
      real*8, intent(in) :: GAO(NBAST), GAO1(NBAST,3), RHG(3),
     &                      VXC, VXB, DFTHRL
      real*8, intent(inout) :: EXCMAT(NBAST,NBAST)
C
      integer :: I, IEND, ISTR, ISYM, J
      real*8 :: FC, FX, FY, FZ, GJ, GVXC
C
C     Exchange-correlation contribution to Kohn-Sham matrix
C
      IF (DOGGA .AND. .NOT.FROMVX) THEN
         FX = 4.0d0*VXB*RHG(1)
         FY = 4.0d0*VXB*RHG(2)
         FZ = 4.0d0*VXB*RHG(3)
         IF (NSYM.EQ.1) THEN
            DO J = 1, NBAST
               FC = VXC*GAO(J)+FX*GAO1(J,1)+FY*GAO1(J,2)+FZ*GAO1(J,3)
               IF (abs(FC).GT.DFTHRL) THEN
                  CALL DAXPY(NBAST,FC,GAO,1,EXCMAT(1,J),1)
               END IF
            END DO
         ELSE
            DO ISYM = 1, NSYM
               ISTR = IBAS(ISYM) + 1
               IEND = IBAS(ISYM) + NBAS(ISYM)
               DO J = ISTR, IEND 
                  FC = VXC*GAO(J)+FX*GAO1(J,1)+FY*GAO1(J,2)+FZ*GAO1(J,3)
                  IF (abs(FC).GT.DFTHRL) THEN
                     DO I = ISTR, IEND 
                        EXCMAT(I,J) = EXCMAT(I,J) + FC*GAO(I)
                     END DO
                  END IF
               END DO
            END DO
         END IF
      ELSE
         DO ISYM = 1, NSYM
            ISTR = IBAS(ISYM) + 1
            IEND = IBAS(ISYM) + NBAS(ISYM)
            DO J = ISTR, IEND 
               GJ   = GAO(J)
               GVXC = D2*VXC*GJ
               IF (abs(GVXC).GT.DFTHRL) THEN
                  DO I = ISTR, J - 1
                     EXCMAT(I,J) = EXCMAT(I,J) + GVXC*GAO(I)
                  END DO
                  EXCMAT(J,J) = EXCMAT(J,J) + VXC*GJ*GJ 
               END IF
            END DO
         END DO
      END IF
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTFRC(DMAT,GAO,GAO1,GAO2,VXC,VXB,RHX,RHY,RHZ,DOGGA)
C*****************************************************************************
C
C     Exchange-correlation contribution to molecular gradient
C
C     T. Helgaker sep 99/oct 00/feb 01
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
      LOGICAL DOGGA
      REAL*8  DMAT(NBAST,NBAST), 
     &        GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6)
C
#include "inforb.h"
#include "energy.h"
#include "symmet.h"
#include "shells.h"
C
      V2 = D2*VXC
      DO IX = 1, 3
         IF (IX.EQ.1) THEN
            K1 = 1
            K2 = 2
            K3 = 3
         ELSE IF (IX.EQ.2) THEN
            K1 = 2
            K2 = 4
            K3 = 5
         ELSE
            K1 = 3
            K2 = 5
            K3 = 6
         END IF
         DO IREPA = 0, MAXREP
            ISTR = IBAS(IREPA+1) + 1
            IEND = IBAS(IREPA+1) + NBAS(IREPA+1)
            IRPAX = IEOR(IREPA,ISYMAX(IX,1))
            IORBA = 0
            DO ISHELA = 1, KMAX
               ISCOOR = IPTCNT(3*(NCENT(ISHELA) - 1) + IX,0,1)
               DO ICOMPA = 1, KHKT(ISHELA)
                  IORBA = IORBA + 1
                  IA = IPTSYM(IORBA,IREPA)
                  KA = IPTSYM(IORBA,IRPAX)
                  IF (KA.GT.0) THEN
                     IF (DOGGA) THEN
                        GA  = GAO1(KA,IX)
                        GAX = RHX*GA
                        GAY = RHY*GA
                        GAZ = RHZ*GA
                        GA2 = RHX*GAO2(KA,K1) + RHY*GAO2(KA,K2)
     &                                        + RHZ*GAO2(KA,K3)
                        GD = D0 
                        GF = D0
                        DO IB = ISTR, IEND 
                           GD = GD + DMAT(IB,IA)*GAO(IB)
                           GF = GF + DMAT(IB,IA)*(GAO(IB)*GA2
     &                                          + GAO1(IB,1)*GAX
     &                                          + GAO1(IB,2)*GAY
     &                                          + GAO1(IB,3)*GAZ)
                        END DO
                        FRC = V2*GD*GA + VXB*GF
                     ELSE
                        GD = D0 
                        DO IB = ISTR, IEND 
                           GD = GD + GAO(IB)*DMAT(IB,IA)
                        END DO
                        FRC = V2*GD*GAO1(KA,IX)
                     END IF
                     GRADFT(ISCOOR) = GRADFT(ISCOOR) - FRC 
                  END IF
               END DO
            END DO
         END DO
      END DO
      RETURN
      END


C*****************************************************************************
      SUBROUTINE DFTLTR(KSYMOP,DTRMAT,EXCMAT,GAO,GAO1,C0,C1,C2,HES,
     &                  RHO,RHO13,RHOGRD,RHG,WDRC,WVWN,WBCK,WLYP,DODRC,
     &                  DOVWN,DOBCK,DOLYP,DOGGA,TRPLET,DOHES,DTGAO)
C*****************************************************************************
C
C     T. Helgaker sep 99/oct 00
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
C
#include "inforb.h"
#include "nuclei.h"
#include "dftcom.h"
C
      INTEGER A, B
      LOGICAL DODRC, DOVWN, DOBCK, DOLYP, DOGGA, TRPLET, DOHES
      DIMENSION DTRMAT(NBAST,NBAST),
     &          GAO(NBAST), GAO1(NBAST,3), 
     &          EXCMAT(NBAST,NBAST), 
     &          C0(NORBT), C1(NORBT,3), C2(NORBT),
     &          HES(NOCCT,NVIRT,NOCCT,NVIRT),RHG(3),DTGAO(NBAST)
      DIMENSION B3(3)
C
      IF (DOHES .AND. NSYM.NE.1) THEN
         WRITE (LUPRI,'(2X,A,/A)') 
     &      ' Symmetry not implemented for explicit Hessian in DFTLTR.',
     &      ' Program aborted.'
         CALL QUIT('Symmetry not implementd for DOHES in DFTLTR')
      END IF
      IF (DOHES) THEN
         B0 = D1
      ELSE
         CALL DGEMV('N',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
         B0 = DDOT(NBAST,DTGAO,1,GAO,1)
      END IF
C
C     ***************
C     ***** LDA *****
C     ***************
C
      IF (.NOT.DOGGA) THEN
         IF (abs(B0).GT.DFTHRL) THEN 
C
C           Calculated VT
C
            IF (DODRC) THEN
               CALL V1DRC(VDRC0,VDRC1,RHO,RHO13)
            ELSE
               VDRC0 = 0D0
               VDRC1 = 0D0
            ENDIF
            IF (DOVWN) THEN 
               IF (.NOT.TRPLET) THEN
                  CALL V1VWN(VVWN0,FRRVWN,RHO,RHO13)
               ELSE
                  CALL VTVWN(FRRVWN,RHO,RHO13)
                  FRRVWN = DP5*FRRVWN 
               END IF
            ELSE
               VVWN0 = 0D0
               FRRVWN = 0D0
            END IF
            VT = (WDRC*VDRC1 + WVWN*FRRVWN)*B0
C
C           Linear transformation
C
            IF (.NOT.DOHES) THEN
               IF (NSYM.EQ.1) THEN
                  DO I = 1, NBAST
                     GVI = VT*GAO(I)
                     DO J = 1, I
                        EXCMAT(J,I) = EXCMAT(J,I) + GVI*GAO(J)
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
                           GVI = VT*GAO(I)
                           DO J = JSTR, MIN(I,JEND) 
                              EXCMAT(J,I) = EXCMAT(J,I) + GVI*GAO(J)
                           END DO
                        END DO
                     END IF
                  END DO
               END IF
C
C           Explicit Hessian
C
            ELSE
               CALL DGEMM('T','N',NORBT,1,NBAST,1.D0,
     &                    DTRMAT,NBAST,
     &                    GAO,NBAST,0.D0,
     &                    C0,NORBT)
               DO 300 B = NOCCT + 1, NORBT
               DO 300 J = 1, NOCCT
               DO 300 A = NOCCT + 1, NORBT
               DO 300 I = 1, NOCCT
                  IA = A - NOCCT
                  IB = B - NOCCT
                  HES(I,IA,J,IB) = HES(I,IA,J,IB) 
     &                           + VT*C0(A)*C0(I)*C0(B)*C0(J)
  300          CONTINUE
            END IF
         END IF
      ELSE
C
C        ***************
C        ***** GGA *****
C        ***************
C
C        B0, BX=B3(1), BY=B3(2), BZ=B3(3), BMAX
C
         IF (DOHES) THEN
            BMAX = D1
         ELSE
C           B3 = GAO1'*DTGAO
            CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D0,B3,1)
C           DTGAO= DTRMAT'*GAO
            CALL DGEMV('T',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
C           B3 = B3 + GAO1'*DTGAO
            CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D1,B3,1)
            BMAX = MAX(abs(B0),abs(B3(1)),abs(B3(2)),abs(B3(3)))
         END IF
C
         IF (BMAX.GT.DFTHRL) THEN
C
C           ZNV, FZ0, FRR, FRZ, FZZ
C
            IF (DODRC) CALL V1DRC(VDRC0,FRRDRC,RHO,RHO13)
            IF (DOBCK) CALL V1BCK(FR0BCK,FZ0BCK,FRRBCK,FRZBCK,
     &                            FZZBCK,RHO,RHOGRD)
            IF (DOVWN) THEN
               IF (.NOT.TRPLET) THEN
                  CALL V1VWN(VVWN0,FRRVWN,RHO,RHO13)
               ELSE
                  CALL VTVWN(FRRVWN,RHO,RHO13)
                  FRRVWN = DP5*FRRVWN 
               END IF
            END IF 
            IF (DOLYP) THEN
               RHOA = DP5*RHO
               RHGA = (DP5*RHOGRD)**2
C              CALL V1LYP(FR0LYP,FZ0LYP,FRRLYP,FRZLYP,
C    &                    FZZLYP,RHO,RHOGRD)
               CALL GLYPCO(DF1000,DF0100,DF0010,DF0001,
     &                     DF00001,RHO,RHO13,RHOGRD,.TRUE.)
               CALL VTLYP (DF2000,DF0200,DF1100,DF1010,
     &                     DF0101,DF1001,DF0110,DF10001,
     &                     DF01001,RHOA,RHOA,RHGA,RHGA,RHGA)
               IF (.NOT.TRPLET) THEN
                  FZ0LYP = DP5*(DF0010 + DF00001)*RHOGRD
                  FRRLYP = DP5*(DF2000 + DF1100)
                  FRZLYP = DP5*(DF1010 + DF1001+DF10001)*RHOGRD
                  FZZLYP = FZ0LYP/RHOGRD 
               ELSE
                  FZ0LYP = DP5*(DF0010 - DF00001)*RHOGRD
                  FRRLYP = DP5*(DF2000 - DF1100)
                  FRZLYP = DP5*(DF1010 - DF1001)*RHOGRD
                  FZZLYP = FZ0LYP/RHOGRD 
               END IF
            ELSE 
               FZ0LYP = 0D0
               FRRLYP = 0D0
               FRZLYP = 0D0
               FZZLYP = 0D0
            END IF
C
            ZNV = D1/RHOGRD
            FZ0 = ZNV*(WBCK*FZ0BCK + WLYP*FZ0LYP)
            FRR = WDRC*FRRDRC + WVWN*FRRVWN 
     &          + WBCK*FRRBCK + WLYP*FRRLYP 
            FRZ = WBCK*FRZBCK + WLYP*FRZLYP 
            FZZ = WBCK*FZZBCK + WLYP*FZZLYP 
C 
            RX = ZNV*RHG(1) 
            RY = ZNV*RHG(2)
            RZ = ZNV*RHG(3)
C
C           Linear transformation
C
            IF (.NOT.DOHES) THEN
               BR = B3(1)*RX + B3(2)*RY + B3(3)*RZ
               FAC0 = FRR*B0 + FRZ*BR
               FACR = FRZ*B0 + FZZ*BR
               IF (NSYM.EQ.1) THEN
                  DO I = 1, NBAST
                     G0 = GAO(I)
                     GX = GAO1(I,1)
                     GY = GAO1(I,2)
                     GZ = GAO1(I,3)
                     DO J = 1, I 
                        A0 = G0*GAO(J)
                        AX = GX*GAO(J) + G0*GAO1(J,1)
                        AY = GY*GAO(J) + G0*GAO1(J,2)
                        AZ = GZ*GAO(J) + G0*GAO1(J,3)
                        AR = AX*RX + AY*RY + AZ*RZ
                        AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) - AR*BR
                        EXCMAT(J,I) = EXCMAT(J,I)+FAC0*A0+FACR*AR+FZ0*AB
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
                           GX = GAO1(I,1)
                           GY = GAO1(I,2)
                           GZ = GAO1(I,3)
                           DO J = JSTR, MIN(I,JEND) 
                              A0 = G0*GAO(J)
                              AX = GX*GAO(J) + G0*GAO1(J,1)
                              AY = GY*GAO(J) + G0*GAO1(J,2)
                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
                              AR = AX*RX + AY*RY + AZ*RZ
                              AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) -AR*BR
                              EXCMAT(J,I) = EXCMAT(J,I) + FAC0*A0 
     &                                    + FACR*AR + FZ0*AB
                           END DO
                        END DO
                     END IF
                  END DO
               END IF
C
C           Explicit Hessian
C
            ELSE
               CALL DGEMM('T','N',NORBT,1,NBAST,1.D0,
     &                    DTRMAT,NBAST,
     &                    GAO,NBAST,0.D0,
     &                    C0,NORBT)
               CALL DGEMM('T','N',NORBT,3,NBAST,1.D0,
     &                    DTRMAT,NBAST,
     &                    GAO1,NBAST,0.D0,
     &                    C1,NORBT)
C
               DO I = 1, NORBT
                 C2(I) = RX*C1(I,1)+RY*C1(I,2)+RZ*C1(I,3)
               END DO
C
               DO B = NOCCT + 1, NORBT 
               DO J = 1, NOCCT 
                 IB  = B - NOCCT
                 GBJ = C0(B)*C0(J)
                 PBJ = C2(B)*C0(J)   + C2(J)*C0(B)
                 CBX = C0(B)*C1(J,1) + C1(B,1)*C0(J)
                 CBY = C0(B)*C1(J,2) + C1(B,2)*C0(J) 
                 CBZ = C0(B)*C1(J,3) + C1(B,3)*C0(J)
                 DO A = NOCCT + 1, B
                 DO I = 1, NOCCT
                   IA  = A - NOCCT
                   GAI = C0(A)*C0(I)
                   PAI = C2(A)*C0(I) + C2(I)*C0(A)
                   CAB = CBX*(C0(A)*C1(I,1) + C1(A,1)*C0(I))
     &                 + CBY*(C0(A)*C1(I,2) + C1(A,2)*C0(I))
     &                 + CBZ*(C0(A)*C1(I,3) + C1(A,3)*C0(I))
                   HES(I,IA,J,IB) = HES(I,IA,J,IB) 
     &                            + FRR*GAI*GBJ
     &                            + FRZ*(PAI*GBJ + GAI*PBJ)
     &                            + FZZ*PAI*PBJ
     &                            + FZ0*(CAB - PAI*PBJ)
                  END DO
                  END DO
               END DO
               END DO
            END IF
         END IF
      END IF
      RETURN
      END
