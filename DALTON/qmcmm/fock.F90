module qmcmm_fock

   implicit none

   public qmnpmm_fock
   public scf_qmnpmm_out

   private

contains

      SUBROUTINE QMNPMM_FOCK(DCAO,DVAO,RCPMAT,RMMMAT,FMAT,EQMNP,WORK,   &
     &                       LWORK)
!
! Purpose:
!     Computes QM/NP/MM contribution to Fock matrix.
!
! Input:
!   DCAO   - Inactive density matrix.
!   DVAO   - Active density matrix.
!   RCPMAT - Inverted Relay matrix or initial quess for MQ vector (NP)
!   RMMMAT - Inverter Relay matrix or initial quess for M vector (MM)
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array
! Output:
!   FMAT   - QM/NP/MM contribution to Fock matrix
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
      use qmcmm, only: getdim_relmat

#include "qmnpmm.h"
#include "priunit.h"
#include "inforb.h"
!
      integer :: lwork
      real(8) :: DCAO(*), DVAO(*), RCPMAT(*), RMMMAT(*), FMAT(*),        &
     &          WORK(LWORK)

      real(8), allocatable :: mqvec(:)
      real(8), allocatable :: fvvec(:)
      real(8), allocatable :: fao(:)

      real(8) :: eqmnp
      integer :: i
      integer :: idimension
!
      idimension = getdim_relmat(.false.)
!
      EQMNP = 0.0d0
!
      IF (.NOT.MQITER) THEN
         IF (.NOT.(DOMMSUB.AND.DOMMPOL)) THEN
            allocate(mqvec(idimension))
            allocate(fvvec(idimension))
            allocate(fao(nnbasx))
!           Determine electric field/potential vector
            CALL GET_FVVEC(DCAO,DVAO,FVVEC,idimension,WORK,     &
     &                     lwork)
!           Determine induced momemnts/charges
            ! rcpmat becomes complex -> complex fock operator
            CALL DGEMV('N',idimension,idimension,1.0d0,RCPMAT,idimension,FVVEC,1,0.0d0,  &
     &                 MQVEC,1)
            deallocate(fvvec)
          if (iprtlvl > 14) then
             write(lupri, '(/,2x,a)') '*** Computed MQ vector start ***'
             do i = 1, idimension
                write(lupri, '(i8, f18.8)') i, mqvec(i)
             end do
             write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
          end if
!           Compute induced dipoles & charges contribution to
!           Fock/Kohn-Sham matrix
            CALL GET_INDMQ_FOCK(DCAO,DVAO,MQVEC,idimension,FAO, &
     &                          WORK,lwork)
            deallocate(mqvec)
!           Add energy contributions
            IF (DONPPOL) EQMNP = EQMNP+EESOLMNP+ENSOLMNP
            IF (DONPCAP) EQMNP = EQMNP+EESOLQNP+ENSOLQNP
!           Computer permanent charges in MM region contribution to
!           Fock/Kohn-Sham matrix and add energy contributions
            IF (DOMMSUB) THEN
               CALL GET_PERMQ_FOCK(DCAO,DVAO,FAO,WORK,    &
     &                             lwork)
               EQMNP = EQMNP+EESOLQMM+ENSOLQMM
            END IF
         ELSE
             stop 'not implemented'
         END IF
!        Pack matrix and add it to Fock matrix
         CALL PKSYM1(FAO,FMAT,NBAS,NSYM,+1)
         if (allocated(fao)) deallocate(fao)
      END IF

      end subroutine
      SUBROUTINE GET_FVVEC(DCAO,DVAO,FVVEC,idimension,WORK,LWORK)
!
! Purpose:
!     Computes electric field/potential vector.
!
! Input:
!   DCAO   - Inactive density matrix.
!   DVAO   - Active density matrix.
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array.
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   idimension   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "priunit.h"
#include "qmnpmm.h"
!
      integer :: idimension, lwork
      real(8) :: DCAO(*), DVAO(*), FVVEC(idimension), WORK(LWORK)
!
      integer :: i
!
      fvvec = 0.0d0
!     determine QM region contributions to FV vector
      CALL GET_QMNUCFV(FVVEC,idimension)
      CALL GET_QMELEFV(FVVEC,idimension,DCAO,DVAO,WORK,lwork)
!
      CALL GET_QLAGRAN(FVVEC,idimension)
!     Print final FV vector
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
         write(lupri, '(/,2x,a)') '*** Computed FV vector start ***'
         do i = 1, idimension
            write(lupri, '(i8, f18.8)') i, fvvec(i)
         end do
         write(lupri, '(/,2x,a)') '*** Computed FV vector end ***'
      END IF
!     Add MM region contribution to FV vector
      IF (DOMMSUB.AND.(.NOT.DOMMPOL)) THEN
         CALL GET_MMFV(FVVEC,idimension)
!         Print final FV vector with MM contribution
          IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed FV+MM vector ***'
           CALL OUTPUT(FVVEC,1,idimension,1,1,idimension,1,1,LUPRI)
          END IF
      END IF
!
      end subroutine

   pure subroutine get_qmnucfv(fvvec, idimension)
!
! Purpose:
!     Computes contribution to electric field/potential vector from
!     nuclei in QM region.
!
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   idimension   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "priunit.h"
#include "qmnpmm.h"
#include "mxcent.h"
#include "nuclei.h"
!
      integer, intent(in)    :: idimension
      real(8), intent(inout) :: fvvec(idimension)

      real(8) :: rij(3)
      real(8) :: fact, rad, rad3
      integer :: i, j
      integer :: ioff
      integer :: istart

!     Electric field due to nuclei in QM region
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            IOFF = 3*(I-1)
            DO J=1,NUCIND
              RIJ(1) = NPCORD(1,I)-CORD(1,J)
              RIJ(2) = NPCORD(2,I)-CORD(2,J)
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD
              FACT = CHARGE(J)/RAD3
              FVVEC(IOFF+1) = FVVEC(IOFF+1)+FACT*RIJ(1)
              FVVEC(IOFF+2) = FVVEC(IOFF+2)+FACT*RIJ(2)
              FVVEC(IOFF+3) = FVVEC(IOFF+3)+FACT*RIJ(3)
            END DO
         END DO
      END IF
!     Determine potential if capacitancy is used
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            DO J=1,NUCIND
              RIJ(1) = NPCORD(1,I)-CORD(1,J)
              RIJ(2) = NPCORD(2,I)-CORD(2,J)
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              FVVEC(IOFF) = FVVEC(IOFF)+CHARGE(J)/RAD
            END DO
!           Fix me: add MM region polarizable points
         END DO
      END IF
!
      end subroutine
      SUBROUTINE GET_QMELEFV(FVVEC,idimension,DCAO,DVAO,WORK,LWORK)
!
! Purpose:
!     Computes contributio to electric field/potential vector
!     for electrons in QM regions.
!
! Input:
!   DCAO   - Inactive density matrix.
!   DVAO   - Active density matrix.
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array.
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   idimension   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"

!
      integer :: idimension, lwork
      real(8) :: DCAO(*), DVAO(*), FVVEC(idimension), WORK(LWORK)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      integer :: INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      real(8) :: RSAVORG(3)

      integer :: i, j, ioff, istart
      integer :: kpatom, nocomp
      real(8) :: relfld, rvpot
      real(8), external :: ddot
      real(8), allocatable :: intao(:)
!
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!     Electric field due to electrons in QM region
      IF (DONPPOL) THEN
         allocate(intao(3*nnbasx))
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 3
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
!
         DO I=1,TNPATM
             IOFF = (I-1)*3
             DIPORG(1) = NPCORD(1,I)
             DIPORG(2) = NPCORD(2,I)
             DIPORG(3) = NPCORD(3,I)
             intao = 0.0d0
             CALL GET1IN(INTAO,'NEFIELD',NOCOMP,WORK,     &
     &                   lwork,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             IF (NISHT.GT.0) THEN
                RELFLD = DDOT(NNBASX,DCAO,1,INTAO,1)
                FVVEC(IOFF+1) = FVVEC(IOFF+1)+RELFLD
                RELFLD = DDOT(NNBASX,DCAO,1,INTAO(NNBASX + 1),1)
                FVVEC(IOFF+2) = FVVEC(IOFF+2)+RELFLD
                RELFLD = DDOT(NNBASX,DCAO,1,INTAO(2*NNBASX + 1),1)
                FVVEC(IOFF+3) = FVVEC(IOFF+3)+RELFLD
             END IF
             IF (NASHT.GT.0) THEN
                RELFLD = DDOT(NNBASX,DVAO,1,INTAO,1)
                FVVEC(IOFF+1) = FVVEC(IOFF+1)+RELFLD
                RELFLD = DDOT(NNBASX,DVAO,1,INTAO(NNBASX + 1),1)
                FVVEC(IOFF+2) = FVVEC(IOFF+2)+RELFLD
                RELFLD = DDOT(NNBASX,DVAO,1,INTAO(2*NNBASX + 1),1)
                FVVEC(IOFF+3) = FVVEC(IOFF+3)+RELFLD
             END IF
         END DO
         RUNQM3 = .FALSE.
         deallocate(intao)
      END IF
!     Potential due to electrons in QM region
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
         allocate(intao(nnbasx))
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
!
         DO I=1,TNPATM
             IOFF = ISTART+I
             DIPORG(1) = NPCORD(1,I)
             DIPORG(2) = NPCORD(2,I)
             DIPORG(3) = NPCORD(3,I)
             intao = 0.0d0
             CALL GET1IN(INTAO,'NPETES ',NOCOMP,WORK,     &
     &                   lwork,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             IF (NISHT.GT.0) THEN
                RVPOT = DDOT(NNBASX,DCAO,1,INTAO,1)
                FVVEC(IOFF) = FVVEC(IOFF)+RVPOT
             END IF
             IF (NASHT.GT.0) THEN
                RVPOT = DDOT(NNBASX,DVAO,1,INTAO,1)
                FVVEC(IOFF) = FVVEC(IOFF)+RVPOT
             END IF
         END DO
         RUNQM3 = .FALSE.
         deallocate(intao)
      END IF
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!
      end subroutine

   pure subroutine get_qlagran(fvvec, idimension)
!
! Purpose:
!     Determines charge contrain fro electric/field potential vector.
!
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   idimension   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "qmnpmm.h"
!
      integer, intent(in)    :: idimension
      real(8), intent(inout) :: FVVEC(idimension)

      integer :: i
!
      IF (DONPCAP) THEN
         DO I=1,TNPBLK
            FVVEC(idimension) = FVVEC(idimension)+NPCHRG(I)
         END DO
      END IF
!
   end subroutine


   pure subroutine get_mmfv(fvvec, idimension)
!
! Purpose:
!     Computes contribution to electric field/potential vector from
!     permanent point charges in MM region.
!
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   idimension   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "priunit.h"
#include "qmnpmm.h"
!
      integer, intent(in)    :: idimension
      real(8), intent(inout) :: fvvec(idimension)

      real(8) :: rij(3)
      real(8) :: fact, rad, rad3
      integer :: i, j
      integer :: ioff
      integer :: istart

!     Electric field due to point charges in MM region
      IF (DONPPOL) THEN
         DO I=1,TNPATM
            IOFF = 3*(I-1)
            DO J=1,TMMATM
              RIJ(1) = NPCORD(1,I)-mm_cord(1,J)
              RIJ(2) = NPCORD(2,I)-mm_cord(2,J)
              RIJ(3) = NPCORD(3,I)-mm_cord(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD
              FACT = MMFM0(MMFTYP(J))/RAD3
              FVVEC(IOFF+1) = FVVEC(IOFF+1)+FACT*RIJ(1)
              FVVEC(IOFF+2) = FVVEC(IOFF+2)+FACT*RIJ(2)
              FVVEC(IOFF+3) = FVVEC(IOFF+3)+FACT*RIJ(3)
            END DO
         END DO
      END IF
!     Determine potential if capacitancy is used
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
         DO I=1,TNPATM
            IOFF = ISTART+I
            DO J=1,TMMATM
              RIJ(1) = NPCORD(1,I)-mm_cord(1,J)
              RIJ(2) = NPCORD(2,I)-mm_cord(2,J)
              RIJ(3) = NPCORD(3,I)-mm_cord(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              FVVEC(IOFF) = FVVEC(IOFF)+MMFM0(MMFTYP(J))/RAD
            END DO
!           Fix me: add MM region polarizable points
         END DO
      END IF
!
   end subroutine


      SUBROUTINE GET_INDMQ_FOCK(DCAO,DVAO,MQVEC,idimension,FCAO,WORK,LWORK)
!
! Purpose:
!     Computes electric field/potential vector.
!
! Input:
!   DCAO   - Inactive density matrix.
!   DVAO   - Active density matrix.
!   MQVEC  - Induced dipoles & charges vector.
!   idimension   - Size of induced dipoles & charges vector.
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array.
! Output:
!   FCAO   - Computed Fock matrix (AO basis).
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
      use qmcmm, only: comp_dampvmat

#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "pi.h"
!
      integer :: lwork, idimension
      real(8) :: DCAO(*), DVAO(*), MQVEC(idimension), FCAO(*), WORK(LWORK)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      integer :: INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      real(8) :: RSAVORG(3), RIJ(3)

      integer :: kpatom
      integer :: i, j, ioff, istart, nocomp
      real(8) :: rad, rad3, rfact, rdval, fact, ripol, rqval, rdim
      real(8) :: enucfx, enucfy, enucfz
      real(8) :: elex, eley, elez
      real(8), external :: derf
      real(8), external :: ddot
      real(8), allocatable :: intao(:)
!
      CALL DZERO(FCAO,NNBASX)
!
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!
      ENSOLQNP = 0.0d0
      EESOLQNP = 0.0d0
      ENSOLMNP = 0.0d0
      EESOLMNP = 0.0d0
!     Induced dipole moment in NP region interaction with QM region
      IF (DONPPOL.AND.NOVDAMP) THEN
!        Electronic interaction part
         allocate(intao(3*nnbasx))
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 3
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
         DO I=1,TNPATM
             IOFF = 3*(I-1)
             DIPORG(1) = NPCORD(1,I)
             DIPORG(2) = NPCORD(2,I)
             DIPORG(3) = NPCORD(3,I)
             intao = 0.0d0
             CALL GET1IN(INTAO,'NEFIELD',NOCOMP,WORK,     &
     &                   lwork,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             CALL DSCAL(NNBASX,MQVEC(IOFF+1),INTAO,1)
             CALL DSCAL(NNBASX,MQVEC(IOFF+2),INTAO(NNBASX + 1),1)
             CALL DSCAL(NNBASX,MQVEC(IOFF+3),INTAO(2*NNBASX + 1),1)
             CALL DAXPY(NNBASX,-1.0d0,INTAO,1,FCAO,1)
             CALL DAXPY(NNBASX,-1.0d0,INTAO(NNBASX + 1),1,FCAO,1)
             CALL DAXPY(NNBASX,-1.0d0,INTAO(2*NNBASX + 1),1,FCAO,1)
             IF (NISHT.GT.0) THEN
               ELEX = DDOT(NNBASX,DCAO,1,INTAO,1)
               ELEY = DDOT(NNBASX,DCAO,1,INTAO(NNBASX + 1),1)
               ELEZ = DDOT(NNBASX,DCAO,1,INTAO(2*NNBASX + 1),1)
               EESOLMNP = EESOLMNP+ELEX+ELEY+ELEZ
             END IF
             IF (NASHT.GT.0) THEN
               ELEX = DDOT(NNBASX,DVAO,1,INTAO,1)
               ELEY = DDOT(NNBASX,DVAO,1,INTAO(NNBASX + 1),1)
               ELEZ = DDOT(NNBASX,DVAO,1,INTAO(2*NNBASX + 1),1)
               EESOLMNP = EESOLMNP+ELEX+ELEY+ELEZ
             END IF
         END DO
         RUNQM3 = .FALSE.
         deallocate(intao)
         EESOLMNP = -0.5d0*EESOLMNP
!        Nuclear interaction part
         DO I=1,TNPATM
             IOFF = 3*(I-1)
             ENUCFX = 0.0d0
             ENUCFY = 0.0d0
             ENUCFZ = 0.0d0
             DO J=1,NUCIND
              RIJ(1) = NPCORD(1,I)-CORD(1,J)
              RIJ(2) = NPCORD(2,I)-CORD(2,J)
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD
              FACT = CHARGE(J)/RAD3
              ENUCFX = ENUCFX + FACT*RIJ(1)
              ENUCFY = ENUCFY + FACT*RIJ(2)
              ENUCFZ = ENUCFZ + FACT*RIJ(3)
            END DO
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+1)*ENUCFX
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+2)*ENUCFY
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+3)*ENUCFZ
         END DO
         ENSOLMNP = -0.5d0*ENSOLMNP
      END IF
!     Induced charges in NP region interaction with QM region
      IF (DONPCAP.AND.NOVDAMP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
!        Electronic interaction part
         allocate(intao(nnbasx))
!        Set integrals evaluation flags
         KPATOM = 0
         NOCOMP = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         RUNQM3 = .TRUE.
!
         DO I=1,TNPATM
             IOFF = ISTART+I
             DIPORG(1) = NPCORD(1,I)
             DIPORG(2) = NPCORD(2,I)
             DIPORG(3) = NPCORD(3,I)
             intao = 0.0d0
             CALL GET1IN(INTAO,'NPETES ',NOCOMP,WORK,     &
     &                   lwork,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)

             WRITE(LUPRI,'(/,2X,A)') '*** jaimes AO matrix ***'
             CALL OUTPAK(INTAO,NBAST,1,LUPRI)
             CALL DSCAL(NNBASX,MQVEC(IOFF),INTAO,1)
             CALL DAXPY(NNBASX,1.0d0,INTAO,1,FCAO,1)
             IF (NISHT.GT.0) THEN
               EESOLQNP = EESOLQNP+DDOT(NNBASX,DCAO,1,INTAO,1)
             END IF
             IF (NASHT.GT.0) THEN
               EESOLQNP = EESOLQNP+DDOT(NNBASX,DVAO,1,INTAO,1)
             END IF
         END DO
         RUNQM3 = .FALSE.
         deallocate(intao)
!
         EESOLQNP = 0.5d0*EESOLQNP
!        Nuclear interaction part
         DO I=1,TNPATM
            IOFF = ISTART+I
            DO J=1,NUCIND
              RIJ(1) = NPCORD(1,I)-CORD(1,J)
              RIJ(2) = NPCORD(2,I)-CORD(2,J)
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              ENSOLQNP = ENSOLQNP+MQVEC(IOFF)*CHARGE(J)/RAD
            END DO
         END DO
         ENSOLQNP = 0.5d0*ENSOLQNP
      END IF
!
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!
!     DAMPED INTERACTION CASE
!
!     Induced dipole moment in NP region interaction (damped) with QM region
      IF (DONPPOL.AND.(.NOT.NOVDAMP)) THEN
         allocate(intao(nnbasx))
         CALL COMP_DAMPVMAT(INTAO, MQVEC)
         IF (NISHT.GT.0) THEN
            EESOLQNP = EESOLQNP+DDOT(NNBASX,DCAO,1,INTAO,1)
         END IF
         IF (NASHT.GT.0) THEN
            EESOLQNP = EESOLQNP+DDOT(NNBASX,DVAO,1,INTAO,1)
         END IF
         CALL DAXPY(NNBASX,1.0d0,INTAO,1,FCAO,1)
         deallocate(intao)
         EESOLQNP = 0.5d0*EESOLQNP
!        Nuclear interaction part: induced dipoles and induced charges
         RDIM = DSQRT(2.0d0)/SQRTPI
!        Nuclear interaction part
         DO I=1,TNPATM
             IOFF = 3*(I-1)
             ENUCFX = 0.0d0
             ENUCFY = 0.0d0
             ENUCFZ = 0.0d0
             RIPOL = NPFPOL(NPFTYP(I))/3.0d0
             RDVAL = (RDIM*RIPOL)**(1.0D0/3.0D0)
             DO J=1,NUCIND
              RIJ(1) = NPCORD(1,I)-CORD(1,J)
              RIJ(2) = NPCORD(2,I)-CORD(2,J)
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD
              FACT = CHARGE(J)/RAD3
              RFACT = 2.0d0*RAD*DEXP(-(RAD*RAD)/(RDVAL*RDVAL))
              RFACT = -RFACT/(SQRTPI*RDVAL)
              RFACT = RFACT + DERF(RAD/RDVAL)
              ENUCFX = ENUCFX + FACT*RIJ(1)*RFACT
              ENUCFY = ENUCFY + FACT*RIJ(2)*RFACT
              ENUCFZ = ENUCFZ + FACT*RIJ(3)*RFACT
            END DO
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+1)*ENUCFX
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+2)*ENUCFY
            ENSOLMNP = ENSOLMNP + MQVEC(IOFF+3)*ENUCFZ
         END DO
         ENSOLMNP = -0.5d0*ENSOLMNP
!
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
         DO I=1,TNPATM
            IOFF = ISTART+I
            RQVAL = RDIM*NPFCAP(NPFTYP(I))
            DO J=1,NUCIND
              RIJ(1) = NPCORD(1,I)-CORD(1,J)
              RIJ(2) = NPCORD(2,I)-CORD(2,J)
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RFACT = DERF(RAD/RQVAL)
              ENSOLQNP = ENSOLQNP+MQVEC(IOFF)*CHARGE(J)*RFACT/RAD
            END DO
         END DO
         ENSOLQNP = 0.5d0*ENSOLQNP
      END IF
      end subroutine
      SUBROUTINE GET_PERMQ_FOCK(DCAO,DVAO,FCAO,WORK,LWORK)
!
! Purpose:
!     Computes permanent charges in MM region contribution to Fock matrix.
!
! Input:
!   DCAO   - Inactive density matrix.
!   DVAO   - Active density matrix.
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array.
! Output:
!   FCAO   - Computed Fock matrix (AO basis).
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "inforb.h"
#include "mxcent.h"
#include "iratdef.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
!
      integer :: lwork
      real(8) :: DCAO(*), DVAO(*), FCAO(*), WORK(LWORK)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      integer :: INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      real(8) :: RSAVORG(3), RIJ(3), rad
      integer :: i, j, kpatom
      integer :: nocomp
      real(8), external :: ddot
      real(8), allocatable :: intao(:)
!
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!
      ENSOLQMM = 0.0d0
      EESOLQMM = 0.0d0
!     Induced charges in NP region interaction with QM region
!     Fix me MM region shift
!     Electronic interaction part
      allocate(intao(nnbasx))
!     Set integrals evaluation flags
      KPATOM = 0
      NOCOMP = 1
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      RUNQM3 = .TRUE.
!
      DO I=1,TMMATM
         DIPORG(1) = mm_cord(1,I)
         DIPORG(2) = mm_cord(2,I)
         DIPORG(3) = mm_cord(3,I)
         intao = 0.0d0
         CALL GET1IN(INTAO,'NPETES ',NOCOMP,WORK,         &
     &               lwork,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,        &
     &               TRIMAT,DUMMY,EXP1VL,DUMMY,0)
          CALL DSCAL(NNBASX,MMFM0(MMFTYP(I)),INTAO,1)
          CALL DAXPY(NNBASX,1.0d0,INTAO,1,FCAO,1)
          IF (NISHT.GT.0) THEN
              EESOLQMM = EESOLQMM+DDOT(NNBASX,DCAO,1,INTAO,1)
          END IF
          IF (NASHT.GT.0) THEN
              EESOLQMM = EESOLQMM+DDOT(NNBASX,DVAO,1,INTAO,1)
          END IF
      END DO
      RUNQM3 = .FALSE.
      deallocate(intao)
!     Nuclear interaction part
      DO I=1,TMMATM
         DO J=1,NUCIND
            RIJ(1) = mm_cord(1,I)-CORD(1,J)
            RIJ(2) = mm_cord(2,I)-CORD(2,J)
            RIJ(3) = mm_cord(3,I)-CORD(3,J)
            RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
            ENSOLQMM = ENSOLQMM+MMFM0(MMFTYP(I))*CHARGE(J)/RAD
         END DO
      END DO
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!
      end subroutine


   subroutine scf_qmnpmm_out()
! Purpose:
!     Prints QM/NP/MM contribution to SCF energy and it's decomposition into
!     components.
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "priunit.h"
#include "qmnpmm.h"
#include "infopt.h"

      IF (IPRTLVL.GE.15) THEN
        WRITE(LUPRI,'(/5X,A/)')                                         &
     &        'QM/NP/MM calculation converged     :'
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'NP charges nuclear contribution    :', ENSOLQNP
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'NP charges electronic contribution :', EESOLQNP
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'NP dipoles nuclear contribution    :', ENSOLMNP
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'NP dipoles electronic contribution :', EESOLMNP
        IF (DOMMSUB) THEN
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'MM charges nuclear contribution    :', ENSOLQMM
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'MM charges electronic contribution :', EESOLQMM
        END IF
        WRITE(LUPRI,'(5X,A,F20.12/)')                                   &
     &        'Total QM/NP/MM energy              :', ESOLT
      ELSE
        WRITE(LUPRI,'(/5X,A/)')                                         &
     &        'QM/NP/MM calculation converged :'
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'NP charges contribution        :',                       &
     &         ENSOLQNP+EESOLQNP
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'NP dipoles contribution        :',                       &
     &        ENSOLMNP+EESOLMNP
        IF (DOMMSUB) THEN
        WRITE(LUPRI,'(5X,A,F20.12)')                                    &
     &        'MM charges contribution        :',                       &
     &        ENSOLQMM+EESOLQMM
        END IF
        WRITE(LUPRI,'(5X,A,F20.12/)')                                   &
     &        'Total QM/NP/MM energy          :', ESOLT
      END IF
!
   end subroutine

end module
