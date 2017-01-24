!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
