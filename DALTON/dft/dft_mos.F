!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
      SUBROUTINE DFTMOMO(CMO2,WORK,LFREE,IPRINT)
C
C     June 2008; PBHT MO overlap diagnostic
C     Reference:
C       Excitation energies in density functional theory: An evaluation
C       and a diagnostic test
C       Michael J. G. Peach, Peter Benfield, Trygve Helgaker, and David J. Tozer,
C       J.  Chem. Phys. 128, 044118 (2008), DOI:10.1063/1.2831900
C
#include "implicit.h"
#include "dummy.h"
#include "inforb.h"
#include "priunit.h"
      PARAMETER (D1 = 1.0D0)
      EXTERNAL DFTMO2
      DIMENSION WORK(LFREE), CMO2(NORBT,NORBT,2)

      CALL QENTER('DFTMOMO')
      KDMAT = 1
      KCMOS = KDMAT + N2BASX
      KCMO2 = KCMOS + N2BASX
      KLAST = KCMO2 + 2*N2ORBX 
      LWRK  = LFREE - KLAST +1
      IF(KLAST.GT.LFREE) CALL QUIT('Insufficient memory in DFTMOMO')
      CALL DZERO(WORK(KCMO2),2*N2ORBX)
      CALL DFTDNS(WORK(KDMAT),WORK(KLAST),LWRK,0)
      CALL DFTGMO(WORK(KCMOS),IPRINT)
Ckr
C     Kickstart for parallel calculations
C
      CALL KICK_SLAVES_MOMO(NBAST,NORBT,WORK(KDMAT),WORK(KCMOS),IPRINT)
      CALL DFTINT(WORK(KDMAT),1,0,.FALSE.,WORK(KLAST),LWRK,
     &            DFTMO2,WORK(KCMOS),ELE)
      CALL DFT_MOMO_COLLECT(WORK(KCMO2),N2ORBX,WORK(KLAST),LWRK)
      IF(IPRINT.GE.0) WRITE(LUPRI,'(A,F20.14)')
     &     ' Electrons in DFTMOMO:', ELE
      CALL DCOPY(2*N2ORBX,WORK(KCMO2),1,CMO2,1)
C
      DO J = 1, NORBT 
      DO I = 1, NORBT
         IF(I.NE.J)CMO2(I,J,2)=CMO2(I,J,2)/SQRT(CMO2(I,I,2)*CMO2(J,J,2))
      END DO
      END DO
      DO I = 1 ,NORBT
         CMO2(I,I,2) = D1
      END DO 
C
      IF (IPRINT.GE.15) THEN
         CALL HEADER('Matrix < abs(phi_p) abs(phi_q) >',-1)
         CALL OUTPUT(CMO2,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         CALL HEADER('Matrix < phi_p^2 phi_q^2 >',-1)
         CALL OUTPUT(CMO2(1,1,2),1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
      END IF
C
      CALL QEXIT('DFTMOMO')
      RETURN
      END
C /* Deck dftgmo */
      SUBROUTINE DFTGMO(CMO,IPRINT)
C
C     T. Helgaker 
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "inforb.h"
#include "inftap.h"
C 
      DIMENSION CMO(NCMOT)
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) 
      CALL READT(LUSIFC,NCMOT,CMO)
      IF (IPRINT.GT.100) THEN
         CALL HEADER('MOS in DFTGMO ',-1)
         CALL OUTPUT(CMO,1,NBAST,1,NORBT,NBAST,NORBT,-1,LUPRI)
      END IF
      RETURN
      END
C /* Deck dftmo2 */
      SUBROUTINE DFTMO2(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                  RHOA,GRADA,DST,VFA,XCPOT,COORD,WGHT,WORK)
#include "implicit.h"
#include "inforb.h"
      DIMENSION NBLCNT(*), NBLOCKS(*), GAO(*), WGHT(*), WORK(*)
      DIMENSION DST(*), VFA(*), XCPOT(*)
      KCMOS = 1
      KCMO2 = KCMOS + N2BASX ! must be same off-set as in DFTMOMO
      CALL DFTMOX(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,WGHT,
     &   WORK(KCMOS),WORK(KCMO2))
      RETURN
      END
C /* Deck dftmoX */
      SUBROUTINE DFTMOX(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,WGHT,CMO,CMO2)
#include "implicit.h"
#include "inforb.h"
      INTEGER P, Q
      DIMENSION GAO(NBLEN,NBAST), WGHT(NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          CMOAO(NBLEN,NORBT),
     &          CMO(NCMOT),
     &          CMO2(NORBT,NORBT,2)
C
      CALL DZERO(CMOAO,NBLEN*NORBT)
      DO ISYM = 1, NSYM
         IENDI = ICMO(ISYM) - IBAS(ISYM)
         DO P = IORB(ISYM) + 1, IORB(ISYM) + NORB(ISYM)
            DO JBL = 1, NBLCNT(ISYM)
            DO J = NBLOCKS(1,JBL,ISYM), NBLOCKS(2,JBL,ISYM)
               DO I = 1, NBLEN
                CMOAO(I,P) = CMOAO(I,P) + CMO(IENDI + J)*GAO(I,J)
               END DO
            END DO
            END DO
            IENDI = IENDI + NBAS(ISYM)
         END DO
      END DO
C
      DO Q = 1, NORBT
      DO P = 1, NORBT
         DO I = 1, NBLEN 
            CMO2(P,Q,1) = CMO2(P,Q,1) 
     &                  + WGHT(I)*ABS(CMOAO(I,P))*ABS(CMOAO(I,Q))
#ifdef DFTMOS_DEBUG
! The following two lines produce the MO overlap matrix, that is a unit matrix,
! in CMO2(:,:,2), if the code is correct.
            CMO2(P,Q,2) = CMO2(P,Q,2) 
     &                  + WGHT(I)*CMOAO(I,P)*CMOAO(I,Q)
#else
            CMO2(P,Q,2) = CMO2(P,Q,2) 
     &                  + WGHT(I)*(CMOAO(I,P)**2)*(CMOAO(I,Q)**2)
#endif
         END DO
      END DO
      END DO 
      RETURN
      END
C
      SUBROUTINE KICK_SLAVES_MOMO(NBAST,NORBT,DMAT,CMOS,IPRINT)
#if defined (VAR_MPI)
#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
C defined parallel calculation types  
#include "iprtyp.h"
C     
      DIMENSION DMAT(NBAST,NBAST), CMOS(NBAST,NORBT)
C
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = DFT_MOMO_WORK
         CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL DFTINTBCAST
         CALL MPI_BCAST(DMAT,NBAST*NBAST,MPI_DOUBLE_PRECISION,MASTER,
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(CMOS,NBAST*NORBT,MPI_DOUBLE_PRECISION,MASTER,
     &                  MPI_COMM_WORLD,IERR)
      END IF
      RETURN
#endif
      END
C
#if defined (VAR_MPI)
      SUBROUTINE DFT_MOMO_SLAVE(WORK,LWORK,IPRINT)
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "inforb.h"
C
      DIMENSION WORK(LWORK)
#include "dftcom.h"
#include "mpif.h"
      EXTERNAL DFTMO2
      LOGICAL DFT_ISGGA
      EXTERNAL DFT_ISGGA
C
      KDMAT = 1
      KCMOS = KDMAT + N2BASX
      KCMO2 = KCMOS + N2BASX
      KLAST = KCMO2 + 2*N2ORBX
      LWRK  = LWORK - KLAST +1
      IF(KLAST.GT.LWORK) CALL QUIT('NOMEM IN DFT_MOMO_SLAVE')
      CALL DFTINTBCAST
      CALL MPI_BCAST(WORK(KDMAT),NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WORK(KCMOS),NBAST*NORBT,MPI_DOUBLE_PRECISION,
     &               MASTER,MPI_COMM_WORLD,IERR)
      CALL DZERO(WORK(KCMO2),2*N2ORBX)
      CALL DFTINT(WORK(KDMAT),1,0,.FALSE.,WORK(KLAST),LWRK,
     &            DFTMO2,WORK(KCMOS),ELE)
      CALL DFT_MOMO_COLLECT(WORK(KCMO2),N2ORBX,WORK(KLAST),LWRK)
      RETURN
      END
#endif      
C
      SUBROUTINE DFT_MOMO_COLLECT(CMO2,IDIM,WORK,LWORK)
#if defined (VAR_MPI)
#include "implicit.h"
#include "mxcent.h"
#include "mpif.h"
      DIMENSION CMO2(IDIM,2), WORK(LWORK)
      CALL DCOPY(2*IDIM,CMO2,1,WORK,1)
      CALL MPI_REDUCE(WORK,CMO2,2*IDIM,MPI_DOUBLE_PRECISION,
     &                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      RETURN
#endif
      END
      SUBROUTINE KRDUMMY
      RETURN
      END
