C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
      SUBROUTINE DFTMOMO(CMO2,WORK,LFREE,IPRINT)
#include <implicit.h>
#include <dummy.h>
#include <inforb.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0)
      EXTERNAL DFTMO2
      DIMENSION WORK(LFREE), CMO2(NORBT,NORBT,2)
      KDMAT = 1
      KCMOS = KDMAT + N2BASX
      KCMO2 = KCMOS + N2BASX
      KLAST = KCMO2 + 2*N2ORBX 
      LWRK  = LFREE - KLAST +1
      IF(KLAST.GT.LFREE) CALL QUIT('NOMEM IN DFTMOMO')
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
      DO I = 1, NORBT
      DO J = 1, NORBT 
         IF(I.NE.J)CMO2(I,J,2)=CMO2(I,J,2)/SQRT(CMO2(I,I,2)*CMO2(J,J,2))
      END DO
      END DO
      DO I = 1 ,NORBT
         CMO2(I,I,2) = D1
      END DO 
C
      IF (IPRINT.GT.3) THEN
         CALL HEADER('Matrix < abs(phi_p) abs(phi_q) >',-1)
         CALL OUTPUT(CMO2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         CALL HEADER('Matrix < phi_p^2 phi_q^2 >',-1)
         CALL OUTPUT(CMO2(1,1,2),1, NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
C /* Deck dftgmo */
      SUBROUTINE DFTGMO(CMO,IPRINT)
C
C     T. Helgaker 
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <iratdef.h>
#include <inforb.h>
C 
      DIMENSION CMO(NBAST,NORBT)
C
#include <inftap.h>
#include <infvar.h>
#include <nuclei.h>
#include <dftcom.h>
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) 
      CALL READI(LUSIFC,IRAT*NCMOT,CMO)
      IF (IPRINT.GT.100) THEN
         CALL HEADER('MOS in DFTGMO ',-1)
         CALL OUTPUT(CMO,1,NBAST,1,NORBT,NBAST,NORBT,1,LUPRI)
      END IF
      RETURN
      END
C /* Deck dftmo2 */
      SUBROUTINE DFTMO2(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,
     &                  RHOA,GRADA,COORD,WGHT,WORK)
#include <implicit.h>
#include <inforb.h>
      DIMENSION WORK(*)
      K1 = 1
      K2 = K1 + NBAST*NORBT
      CALL DFTMOX(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,WGHT,WORK(K1),WORK(K2))
      RETURN
      END
C /* Deck dftmoX */
      SUBROUTINE DFTMOX(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO,WGHT,CMO,CMO2)
#include <implicit.h>
#include <inforb.h>
      INTEGER P, Q
      DIMENSION GAO(NBLEN,NBAST), WGHT(NBLEN),
     &          NBLCNT(8),NBLOCKS(2,LDAIB,8),
     &          CMOAO(NBLEN,NORBT),
     &          CMO(NCMOT),
     &          CMO2(NORBT,NORBT,2)
C     
      CALL DZERO(CMOAO,NBLEN*NORBT)
      IEND = 0
      DO ISYM = 1, NSYM
         IENDI = IEND
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
         IEND = IEND + NORB(ISYM)*NBAS(ISYM)
      END DO
C
      DO Q = 1, NORBT
      DO P = 1, NORBT
         DO I = 1, NBLEN 
            CMO2(P,Q,1) = CMO2(P,Q,1) 
     &                  + WGHT(I)*ABS(CMOAO(I,P))*ABS(CMOAO(I,Q))
            CMO2(P,Q,2) = CMO2(P,Q,2) 
     &                  + WGHT(I)*(CMOAO(I,P)**2)*(CMOAO(I,Q)**2)
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
C     
      DIMENSION DMAT(NBAST,NBAST), CMOS(NBAST,NORBT)
C
      IF (MYNUM .EQ. MASTER) THEN
         IPRTYP = 15
         CALL MPI_BCAST(IPRTYP,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(IPRINT,1,MPI_INTEGER,MASTER,MPI_COMM_WORLD,IERR)
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
