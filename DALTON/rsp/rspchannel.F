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
      SUBROUTINE CHANNEL_VAR()
      IMPLICIT NONE
#include "infvar.h"
#include "channel.h"
#include "wrkrsp.h"
      INTEGER IOP, CIOP, J, CORBJ
      INTEGER, allocatable :: CHANNEL_JWOP(:,:)
C
C Core RPA: redefine the excitation manifold,
C restrict to excitations involving the few chosen
C orbitals in CHANNEL_ORB
C
      allocate( CHANNEL_JWOP(2,MAXWOP) )
      CIOP = 0
      DO J = 1, CHANNEL_NORB
         CORBJ = CHANNEL_ORB(J)
         DO IOP = 1,KZWOPT
            IF (JWOP(1,IOP).EQ.CORBJ) THEN
               CIOP = CIOP + 1
               CHANNEL_JWOP(1,CIOP) = JWOP(1,IOP)
               CHANNEL_JWOP(2,CIOP) = JWOP(2,IOP)
            END IF
         END DO
      END DO
      CALL ICOPY(2*CIOP,CHANNEL_JWOP,1,JWOP,1)
      KZWOPT = CIOP
      KZYWOP = 2*CIOP
      NWOPT = CIOP
C
      deallocate( CHANNEL_JWOP )
      RETURN
      END
      SUBROUTINE CHANNEL_VIR()
      IMPLICIT NONE
#include "priunit.h"
#include "infvar.h"
#include "channel.h"
#include "wrkrsp.h"
#include "inforb.h"
      INTEGER IOP, CIOP, J, CORBJ,ISYM,IMAX
      INTEGER MAXVIRT(8)
      INTEGER, allocatable :: CHANNEL_JWOP(:,:)
C
C Core RPA: redefine the excitaiton manifold
C to delete high virtual orbitals
C
      allocate( CHANNEL_JWOP(2,MAXWOP) )
      DO ISYM=1,NSYM
         IF (NASH(ISYM).NE.0) THEN
            WRITE (LUPRI,'(/A)')
     &           'WARNING: virtual channel restriction not tested'
     &           //' for open shells.'
         ENDIF
         IF (CHANNEL_VIRT(ISYM).EQ.-1) THEN
            MAXVIRT(ISYM) = IORB(ISYM) + NISH(ISYM) + NASH(ISYM) +
     &           NSSH(ISYM)
         ELSE
            MAXVIRT(ISYM) = IORB(ISYM) + NISH(ISYM) + NASH(ISYM) +
     &           CHANNEL_VIRT(ISYM)
         ENDIF
         WRITE (LUPRI,*) 'MAXVIRT(',ISYM,')=',MAXVIRT(ISYM)
      ENDDO
      CIOP = 0
      DO IOP = 1,KZWOPT
         DO ISYM=1,NSYM
            IMAX= IORB(ISYM) + NISH(ISYM) + NASH(ISYM)
            IF (JWOP(2,IOP).GE.IMAX.AND.
     &          JWOP(2,IOP).LE.MAXVIRT(ISYM)) THEN
               CIOP = CIOP + 1
               CHANNEL_JWOP(1,CIOP) = JWOP(1,IOP)
               CHANNEL_JWOP(2,CIOP) = JWOP(2,IOP)
c               WRITE(LUPRI,*) 'Keeping',JWOP(1,IOP),JWOP(2,IOP)
               GOTO 10
            ELSE
c               WRITE(LUPRI,*) 'Skipping',JWOP(1,IOP),JWOP(2,IOP)
            END IF
         ENDDO
 10      CONTINUE
      END DO
      CALL ICOPY(2*CIOP,CHANNEL_JWOP,1,JWOP,1)
      KZWOPT = CIOP
      KZYWOP = 2*CIOP
      NWOPT = CIOP
C
      deallocate( CHANNEL_JWOP )
      RETURN
      END
