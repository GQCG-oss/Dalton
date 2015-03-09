module qmcmm_fock

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

#include "implicit.h"
#include "dummy.h"
#include "qmnpmm.h"
#include "priunit.h"
#include "inforb.h"
!
      DIMENSION DCAO(*), DVAO(*), RCPMAT(*), RMMMAT(*), FMAT(*),        &
     &          WORK(LWORK)
!
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)

      real(8), allocatable :: mqvec(:)
      real(8), allocatable :: fvvec(:)
      real(8), allocatable :: fao(:)
!
      KFREE = 1
      LFREE = LWORK
!
      CALL  GETDIM_RELMAT(IDIM,.FALSE.)
!
      EQMNP = D0
!
      IF (.NOT.MQITER) THEN
         IF (.NOT.(DOMMSUB.AND.DOMMPOL)) THEN
            allocate(mqvec(idim))
            allocate(fvvec(idim))
            allocate(fao(nnbasx))
!           Determine electric field/potential vector
            CALL GET_FVVEC(DCAO,DVAO,FVVEC,IDIM,WORK(KFREE),     &
     &                     LFREE)
!           Determine induced momemnts/charges
            ! rcpmat becomes complex -> complex fock operator
            CALL DGEMV('N',IDIM,IDIM,D1,RCPMAT,IDIM,FVVEC,1,D0,  &
     &                 MQVEC,1)
            deallocate(fvvec)
          if (iprtlvl > 14) then
             write(lupri, '(/,2x,a)') '*** Computed MQ vector start ***'
             do i = 1, idim
                write(lupri, '(i8, f18.8)') i, mqvec(i)
             end do
             write(lupri, '(/,2x,a)') '*** Computed MQ vector end ***'
          end if
!           Compute induced dipoles & charges contribution to
!           Fock/Kohn-Sham matrix
            CALL GET_INDMQ_FOCK(DCAO,DVAO,MQVEC,IDIM,FAO, &
     &                          WORK(KFREE),LFREE)
            deallocate(mqvec)
!           Add energy contributions
            IF (DONPPOL) EQMNP = EQMNP+EESOLMNP+ENSOLMNP
            IF (DONPCAP) EQMNP = EQMNP+EESOLQNP+ENSOLQNP
!           Computer permanent charges in MM region contribution to
!           Fock/Kohn-Sham matrix and add energy contributions
            IF (DOMMSUB) THEN
               CALL GET_PERMQ_FOCK(DCAO,DVAO,FAO,WORK(KFREE),    &
     &                             LFREE)
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
!  /* Deck get_fvvec */
      SUBROUTINE GET_FVVEC(DCAO,DVAO,FVVEC,IDIM,WORK,LWORK)
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
!   IDIM   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      DIMENSION DCAO(*), DVAO(*), FVVEC(IDIM), WORK(LWORK)
!
      KFREE = 1
      LFREE = LWORK
!
      CALL DZERO(FVVEC,IDIM)
!     determine QM region contributions to FV vector
      CALL GET_QMNUCFV(FVVEC,IDIM)
      CALL GET_QMELEFV(FVVEC,IDIM,DCAO,DVAO,WORK(KFREE),LFREE)
!
      CALL GET_QLAGRAN(FVVEC,IDIM)
!     Print final FV vector
      IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
         write(lupri, '(/,2x,a)') '*** Computed FV vector start ***'
         do i = 1, idim
            write(lupri, '(i8, f18.8)') i, fvvec(i)
         end do
         write(lupri, '(/,2x,a)') '*** Computed FV vector end ***'
      END IF
!     Add MM region contribution to FV vector
      IF (DOMMSUB.AND.(.NOT.DOMMPOL)) THEN
         CALL GET_MMFV(FVVEC,IDIM)
!         Print final FV vector with MM contribution
          IF ((IPRTLVL.GE.15).AND.(.NOT.MQITER)) THEN
            WRITE(LUPRI,'(/,2X,A)') '*** Computed FV+MM vector ***'
           CALL OUTPUT(FVVEC,1,IDIM,1,1,IDIM,1,1,LUPRI)
          END IF
      END IF
!
      RETURN
      end subroutine
!  /* Deck get_qmnucfv */
      SUBROUTINE GET_QMNUCFV(FVVEC,IDIM)
!
! Purpose:
!     Computes contribution to electric field/potential vector from
!     nuclei in QM region.
!
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   IDIM   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "mxcent.h"
#include "nuclei.h"
!
      DIMENSION FVVEC(IDIM)
!
      DIMENSION RIJ(3)
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
      RETURN
      end subroutine
!  /* Deck get_qmelefv */
      SUBROUTINE GET_QMELEFV(FVVEC,IDIM,DCAO,DVAO,WORK,LWORK)
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
!   IDIM   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
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
      DIMENSION DCAO(*), DVAO(*), FVVEC(IDIM), WORK(LWORK)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3)
!
      KFREE = 1
      LFREE = LWORK
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!     Electric field due to electrons in QM region
      IF (DONPPOL) THEN
         CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
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
             CALL DZERO(WORK(KINTAO),3*NNBASX)
             CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             IF (NISHT.GT.0) THEN
                RELFLD = DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)
                FVVEC(IOFF+1) = FVVEC(IOFF+1)+RELFLD
                RELFLD = DDOT(NNBASX,DCAO,1,WORK(KINTAO+NNBASX),1)
                FVVEC(IOFF+2) = FVVEC(IOFF+2)+RELFLD
                RELFLD = DDOT(NNBASX,DCAO,1,WORK(KINTAO+2*NNBASX),1)
                FVVEC(IOFF+3) = FVVEC(IOFF+3)+RELFLD
             END IF
             IF (NASHT.GT.0) THEN
                RELFLD = DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)
                FVVEC(IOFF+1) = FVVEC(IOFF+1)+RELFLD
                RELFLD = DDOT(NNBASX,DVAO,1,WORK(KINTAO+NNBASX),1)
                FVVEC(IOFF+2) = FVVEC(IOFF+2)+RELFLD
                RELFLD = DDOT(NNBASX,DVAO,1,WORK(KINTAO+2*NNBASX),1)
                FVVEC(IOFF+3) = FVVEC(IOFF+3)+RELFLD
             END IF
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_QMELEFV',WORK,1,1,KFREE,LFREE)
      END IF
!     Potential due to electrons in QM region
      IF (DONPCAP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
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
             CALL DZERO(WORK(KINTAO),NNBASX)
             CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             IF (NISHT.GT.0) THEN
                RVPOT = DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)
                FVVEC(IOFF) = FVVEC(IOFF)+RVPOT
             END IF
             IF (NASHT.GT.0) THEN
                RVPOT = DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)
                FVVEC(IOFF) = FVVEC(IOFF)+RVPOT
             END IF
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_QMELEFV',WORK,1,1,KFREE,LFREE)
      END IF
!     Restore origin coordinates
      DIPORG(1) = RSAVORG(1)
      DIPORG(2) = RSAVORG(2)
      DIPORG(3) = RSAVORG(3)
!
      RETURN
      end subroutine
!  /* Deck get_qlag */
      SUBROUTINE GET_QLAGRAN(FVVEC,IDIM)
!
! Purpose:
!     Determines charge contrain fro electric/field potential vector.
!
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   IDIM   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "qmnpmm.h"
!
      DIMENSION FVVEC(IDIM)
!
      IF (DONPCAP) THEN
         DO I=1,TNPBLK
            FVVEC(IDIM) = FVVEC(IDIM)+NPCHRG(I)
         END DO
      END IF
!
      RETURN
      end subroutine
!  /* Deck get_mmfv */
      SUBROUTINE GET_MMFV(FVVEC,IDIM)
!
! Purpose:
!     Computes contribution to electric field/potential vector from
!     permanent point charges in MM region.
!
! Output:
!   FVVEC  - Electric field/potential vector at NP/MM centers.
!   IDIM   - Size of electric field/potential vector
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
!
      DIMENSION FVVEC(IDIM)
!
      DIMENSION RIJ(3)
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
      RETURN
      end subroutine
!  /* Deck get_fvvec */
      SUBROUTINE GET_INDMQ_FOCK(DCAO,DVAO,MQVEC,IDIM,FCAO,WORK,LWORK)
!
! Purpose:
!     Computes electric field/potential vector.
!
! Input:
!   DCAO   - Inactive density matrix.
!   DVAO   - Active density matrix.
!   MQVEC  - Induced dipoles & charges vector.
!   IDIM   - Size of induced dipoles & charges vector.
!   WORK   - Temporary memory array
!   LWORK  - Size of temporary memory array.
! Output:
!   FCAO   - Computed Fock matrix (AO basis).
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
      use qmcmm, only: comp_dampvmat

#include "implicit.h"
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
      DOUBLE PRECISION MQVEC
!
      DIMENSION DCAO(*), DVAO(*), MQVEC(IDIM), FCAO(*), WORK(LWORK)
!
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, DMP5 = -0.5D0 )
      PARAMETER (DP5 = 0.5D0, D2 = 2.0D0, D3 = 3.0D0)
      PARAMETER (D13 = 1.0D0/3.0D0)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3), RIJ(3)
!
      KFREE = 1
      LFREE = LWORK
!
      CALL DZERO(FCAO,NNBASX)
!
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!
      ENSOLQNP = D0
      EESOLQNP = D0
      ENSOLMNP = D0
      EESOLMNP = D0
!     Induced dipole moment in NP region interaction with QM region
      IF (DONPPOL.AND.NOVDAMP) THEN
!        Electronic interaction part
         CALL MEMGET('REAL',KINTAO,3*NNBASX,WORK,KFREE,LFREE)
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
             CALL DZERO(WORK(KINTAO),3*NNBASX)
             CALL GET1IN(WORK(KINTAO),'NEFIELD',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)
             CALL DSCAL(NNBASX,MQVEC(IOFF+1),WORK(KINTAO),1)
             CALL DSCAL(NNBASX,MQVEC(IOFF+2),WORK(KINTAO+NNBASX),1)
             CALL DSCAL(NNBASX,MQVEC(IOFF+3),WORK(KINTAO+2*NNBASX),1)
             CALL DAXPY(NNBASX,DM1,WORK(KINTAO),1,FCAO,1)
             CALL DAXPY(NNBASX,DM1,WORK(KINTAO+NNBASX),1,FCAO,1)
             CALL DAXPY(NNBASX,DM1,WORK(KINTAO+2*NNBASX),1,FCAO,1)
             IF (NISHT.GT.0) THEN
               ELEX = DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)
               ELEY = DDOT(NNBASX,DCAO,1,WORK(KINTAO+NNBASX),1)
               ELEZ = DDOT(NNBASX,DCAO,1,WORK(KINTAO+2*NNBASX),1)
               EESOLMNP = EESOLMNP+ELEX+ELEY+ELEZ
             END IF
             IF (NASHT.GT.0) THEN
               ELEX = DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)
               ELEY = DDOT(NNBASX,DVAO,1,WORK(KINTAO+NNBASX),1)
               ELEZ = DDOT(NNBASX,DVAO,1,WORK(KINTAO+2*NNBASX),1)
               EESOLMNP = EESOLMNP+ELEX+ELEY+ELEZ
             END IF
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_INDMQ_FOCK',WORK,1,1,KFREE,LFREE)
         EESOLMNP = DMP5*EESOLMNP
!        Nuclear interaction part
         DO I=1,TNPATM
             IOFF = 3*(I-1)
             ENUCFX = D0
             ENUCFY = D0
             ENUCFZ = D0
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
         ENSOLMNP = DMP5*ENSOLMNP
      END IF
!     Induced charges in NP region interaction with QM region
      IF (DONPCAP.AND.NOVDAMP) THEN
         ISTART = 0
         IF (DONPPOL) ISTART = 3*TNPATM
!        Fix me MM region shift
!        Electronic interaction part
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
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
             CALL DZERO(WORK(KINTAO),NNBASX)
             CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),     &
     &                   LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,    &
     &                   TRIMAT,DUMMY,EXP1VL,DUMMY,0)

             WRITE(LUPRI,'(/,2X,A)') '*** jaimes AO matrix ***'
             CALL OUTPAK(WORK(KINTAO),NBAST,1,LUPRI)
             CALL DSCAL(NNBASX,MQVEC(IOFF),WORK(KINTAO),1)
             CALL DAXPY(NNBASX,D1,WORK(KINTAO),1,FCAO,1)
             IF (NISHT.GT.0) THEN
               EESOLQNP = EESOLQNP+DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)
             END IF
             IF (NASHT.GT.0) THEN
               EESOLQNP = EESOLQNP+DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)
             END IF
         END DO
         RUNQM3 = .FALSE.
         CALL MEMREL('GET_INDMQ_FOCK',WORK,1,1,KFREE,LFREE)
!
         EESOLQNP = DP5*EESOLQNP
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
         ENSOLQNP = DP5*ENSOLQNP
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
         CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
         CALL COMP_DAMPVMAT(WORK(KINTAO), MQVEC)
         IF (NISHT.GT.0) THEN
            EESOLQNP = EESOLQNP+DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)
         END IF
         IF (NASHT.GT.0) THEN
            EESOLQNP = EESOLQNP+DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)
         END IF
         CALL DAXPY(NNBASX,D1,WORK(KINTAO),1,FCAO,1)
         CALL MEMREL('GET_INDMQ_FOCK',WORK,1,1,KFREE,LFREE)
         EESOLQNP = DP5*EESOLQNP
!        Nuclear interaction part: induced dipoles and induced charges
         RDIM = DSQRT(D2)/SQRTPI
!        Nuclear interaction part
         DO I=1,TNPATM
             IOFF = 3*(I-1)
             ENUCFX = D0
             ENUCFY = D0
             ENUCFZ = D0
             RIPOL = NPFPOL(NPFTYP(I))/D3
             RDVAL = (RDIM*RIPOL)**D13
             DO J=1,NUCIND
              RIJ(1) = NPCORD(1,I)-CORD(1,J)
              RIJ(2) = NPCORD(2,I)-CORD(2,J)
              RIJ(3) = NPCORD(3,I)-CORD(3,J)
              RAD = DSQRT(RIJ(1)*RIJ(1)+RIJ(2)*RIJ(2)+RIJ(3)*RIJ(3))
              RAD3 = RAD*RAD*RAD
              FACT = CHARGE(J)/RAD3
              RFACT = D2*RAD*DEXP(-(RAD*RAD)/(RDVAL*RDVAL))
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
         ENSOLMNP = DMP5*ENSOLMNP
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
         ENSOLQNP = DP5*ENSOLQNP
      END IF
      RETURN
      end subroutine
!  /* Deck get_permq_fock */
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
#include "implicit.h"
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
      DIMENSION DCAO(*), DVAO(*), FCAO(*), WORK(LWORK)
!
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, DMP5 = -0.5D0 )
      PARAMETER (DP5 = 0.5D0)
!
      LOGICAL TOFILE,TRIMAT,EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      DIMENSION RSAVORG(3), RIJ(3)
!
      KFREE = 1
      LFREE = LWORK
!     Save origin coordinates
      RSAVORG(1) = DIPORG(1)
      RSAVORG(2) = DIPORG(2)
      RSAVORG(3) = DIPORG(3)
!
      ENSOLQMM = D0
      EESOLQMM = D0
!     Induced charges in NP region interaction with QM region
!     Fix me MM region shift
!     Electronic interaction part
      CALL MEMGET('REAL',KINTAO,NNBASX,WORK,KFREE,LFREE)
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
         CALL DZERO(WORK(KINTAO),NNBASX)
         CALL GET1IN(WORK(KINTAO),'NPETES ',NOCOMP,WORK(KFREE),         &
     &               LFREE,LABINT,INTREP,INTADR,I,TOFILE,KPATOM,        &
     &               TRIMAT,DUMMY,EXP1VL,DUMMY,0)
          CALL DSCAL(NNBASX,MMFM0(MMFTYP(I)),WORK(KINTAO),1)
          CALL DAXPY(NNBASX,D1,WORK(KINTAO),1,FCAO,1)
          IF (NISHT.GT.0) THEN
              EESOLQMM = EESOLQMM+DDOT(NNBASX,DCAO,1,WORK(KINTAO),1)
          END IF
          IF (NASHT.GT.0) THEN
              EESOLQMM = EESOLQMM+DDOT(NNBASX,DVAO,1,WORK(KINTAO),1)
          END IF
      END DO
      RUNQM3 = .FALSE.
      CALL MEMREL('GET_PERMQ_FOCK',WORK,1,1,KFREE,LFREE)
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
      RETURN
      end subroutine
!     /* Deck scf_qmnpmm_out */
      SUBROUTINE SCF_QMNPMM_OUT()
! Purpose:
!     Prints QM/NP/MM contribution to SCF energy and it's decomposition into
!     components.
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "qmnpmm.h"
#include "infopt.h"
!
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
      RETURN
      end subroutine

end module
