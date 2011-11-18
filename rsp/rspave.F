!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
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
#include <single.h>
C
#ifdef UNDEF
===========================================================================
/* Comdeck rev_log */
Revision 1.2  2001/01/17 10:41:29  vebjornb
Calls to *MPA*B* in arhpack.F have been replaced with DGEMM calls

Revision 1.2  2000/05/01 12:13:01  hjj
Removed obsolete KAVER test.

Written july / august 1992 Hinne Hettema

Purpose of routines :
=====================
1. perform averaging of start vectors and linear transformations in response
   calculations (RSPAVE)
2. symmetrize the E[2] matrix (RFANTI ... RSPFCA)

List of updates:
================
===========================================================================
#endif
C  /* Deck rspave */
      SUBROUTINE RSPAVE (VECS,NVDIM,NVSIM)
C
C 23-Jul-1992 Hinne Hettema
C
#include <implicit.h>
      DIMENSION VECS(NVDIM,NVSIM)
C
C Used from common blocks:
C  INFINP : SUPSYM
C  INFVAR : NWOPT,JWOPSY
C  INFRSP : RSPSUP
C  WRKRSP : KSYMOP, KZWOPT
C  INFPRI : IPRAVE
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <infpri.h>
C
#include <infrsp.h>
#include <wrkrsp.h>
C
      CALL QENTER('RSPAVE')
C
      IF ( IPRAVE .GT. 20 ) THEN
         WRITE(LUPRI,'(A)')    ' Test output from RSPAVE : '
         WRITE(LUPRI,'(A)')    ' ========================= '
         WRITE(LUPRI,'(A,I8)') ' NVDIM     =  ', NVDIM
         WRITE(LUPRI,'(A,I8)') ' NVSIM     =  ', NVSIM
         WRITE(LUPRI,'(A,I8)') ' NWOPT     =  ', NWOPT
         WRITE(LUPRI,'(A,I8)') ' JWOPSY    =  ', JWOPSY
         WRITE(LUPRI,'(A,I8)') ' KSYMOP    =  ', KSYMOP
         WRITE(LUPRI,'(A,I8)') ' KZWOPT    =  ', KZWOPT
         WRITE(LUPRI,'(A,L8)') ' SUPSYM    =  ', SUPSYM
         WRITE(LUPRI,'(A,L8)') ' RSPSUP    =  ', RSPSUP
         WRITE(LUPRI,'(A)')    ' ========================= '
      END IF
      IF (.NOT. SUPSYM) GO TO 9999
      IF (NWOPT  .EQ. 0) GO TO 9999
      IF (JWOPSY .NE. 1) GO TO 910
      IF (NVDIM  .LT. NWOPT) GO TO 920
      IF (JWOPSY .NE. KSYMOP) GO TO 930
      IF (NWOPT  .NE. KZWOPT) GO TO 940
      IF (SUPSYM) THEN
         IF ( IPRAVE .GT. 80 ) THEN
            WRITE(LUPRI,'(A)')    ' Vector before averaging : '
            WRITE(LUPRI,'(A)')    ' ========================= '
            CALL OUTPUT(VECS,1,NVDIM,1,NVSIM,NVDIM,NVSIM,1,LUPRI)
         END IF
         CALL AVERSS(VECS,NVDIM,NVSIM)
         IF ( IPRAVE .GT. 80 ) THEN
            WRITE(LUPRI,'(A)')    ' Vector after averaging : '
            WRITE(LUPRI,'(A)')    ' ======================== '
            CALL OUTPUT(VECS,1,NVDIM,1,NVSIM,NVDIM,NVSIM,1,LUPRI)
         END IF
      END IF
C
 9999 CONTINUE
      IF (IPRAVE .GT. 80) THEN
         WRITE(LUPRI,'(A)')    ' End of RSPAVE'
      END IF
      CALL QEXIT('RSPAVE')
      RETURN
C
C *** Error sections
C
  910 CONTINUE
      WRITE (LUERR,9010) JWOPSY
      CALL QTRACE(LUERR)
      CALL QUIT('RSPAVE ERROR, operator symmetry JWOPSY .ne. 1')
 9010 FORMAT(/' ERROR-RSPAVE, operator symmetry .ne. 1; JWOPSY =',I8)
C
  920 CONTINUE
      WRITE (LUERR,9020) NVDIM,NWOPT
      CALL QTRACE(LUERR)
      CALL QUIT('RSPAVE ERROR, vector length lt NWOPT')
 9020 FORMAT(/' ERROR-RSPAVE, vector length',I8,' is less than NWOPT',
     *   I8)
C
  930 CONTINUE
      WRITE (LUERR,9020) JWOPSY, KSYMOP
      CALL QTRACE(LUERR)
      CALL QUIT('RSPAVE ERROR, JWOPSY .NE. KSYMOP')
 9030 FORMAT(/' ERROR-RSPAVE, JWOPSY = ',I8,' KSYMOP = ', I8)
C
  940 CONTINUE
      WRITE (LUERR,9020) NWOPT, KZWOPT
      CALL QTRACE(LUERR)
      CALL QUIT('RSPAVE ERROR, NWOPT .NE. KZWOPT')
 9040 FORMAT(/' ERROR-RSPAVE, NWOPT = ',I8,' KZWOPT = ', I8)
C
      END
C  /* Deck rfanti */
      SUBROUTINE RFANTI(NSIM,EVECS,ZYMAT,WRK,LWRK)
C
C Subroutine to obtain antisymmetric part of the
C MCSCF gradient from total MCSCF Fock matrix
C
#include <implicit.h>
C
      DIMENSION EVECS(*),ZYMAT(NORBT,NORBT,NSIM),WRK(LWRK)
C
C Used from common blocks
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inftap.h>
#include <infopt.h>
#include <infpri.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      PARAMETER ( DHALF = 0.50D0 )
C
      CALL QENTER('RFANTI')
C
C     Layout work space
C
      KGRA   = 1
      KFOCK  = KGRA    + NSIM * N2ORBX
      KANTI  = KFOCK   + N2ORBT
      KWRK   = KANTI   + N2ORBX
      LWRK1  = LWRK    - KWRK
C
      IF (LWRK1 .LT. 0 ) THEN
         CALL QTRACE(LUERR)
         CALL QUIT(' Not enough work space in RFANTI ')
      END IF
C
C  -- Get the Fock matrix from the interface file:
C     Comment : we know that RSPCI is False, so no test.
C
      REWIND(LUSIFC)
      CALL MOLLAB(LBSIFC,LUSIFC,LUERR)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READT(LUSIFC,N2ORBT,WRK(KFOCK))
      IF (IPRAVE .GT. 80 ) THEN
         WRITE(LUPRI,'(//A)')  ' Fock matrix for use in RFANTI'
         WRITE(LUPRI,'(A)')    ' ============================='
         DO 10 ISYM = 1, NSYM
            I2ORBI = I2ORB(ISYM)
            NORBI  = NORB(ISYM)
            WRITE(LUPRI,'(/A,I4)')  ' Symmetry representation: ', ISYM
            WRITE(LUPRI,'(A)')      ' =============================='
            CALL OUTPUT(WRK(KFOCK+I2ORBI),1,NORBI,1,NORBI,NORBI,NORBI,
     &                  1,LUPRI)
 10      CONTINUE
      END IF
C
C  -- Get the antisymmetric component of the Fock matrix
C
      CALL DZERO(WRK(KANTI), NSIM*N2ORBX)
      DO 100 ISYM = 1, NSYM
         NORBI  = NORB(ISYM)
         IORBI  = IORB(ISYM)
         I2ORBI = I2ORB(ISYM)
         CALL RSPANT(WRK(KFOCK+I2ORBI),NORBI,IORBI,NORBT,WRK(KANTI))
  100 CONTINUE
C
      IF (IPRAVE .GT.  80) THEN
         WRITE(LUPRI,'(//A)') ' Antisymmetry component'
         WRITE(LUPRI,'(A)')   ' ======================'
         CALL OUTPUT(WRK(KANTI),1,NORBT,1,NORBT,NORBT,NORBT,
     &               1,LUPRI)
         DO 90 ISIM = 1, NSIM
            WRITE(LUPRI,'(//A)') ' ZYMAT matrix used :'
            WRITE(LUPRI,'(A)')   ' ==================='
            CALL OUTPUT(ZYMAT(1,1,ISIM),1,NORBT,1,NORBT,NORBT,NORBT,
     &                  1,LUPRI)
  90     CONTINUE
      END IF
C
C  -- Do matrix multiplication with ZYMAT to get correction
C     0.5 sum(s) ( k(q,s) g(p,s) - g(s,q) k(s,p) )
C
      DO 200 ISIM = 1, NSIM
         CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.D0,
     &              ZYMAT(1,1,ISIM),NORBT,
     &              WRK(KANTI),NORBT,0.D0,
     &              WRK(KGRA+N2ORBX*(ISIM-1)),NORBT)
         CALL DGEMM('T','N',NORBT,NORBT,NORBT,-1.D0,
     &              WRK(KANTI),NORBT,
     &              ZYMAT(1,1,ISIM),NORBT,1.D0,
     &              WRK(KGRA+N2ORBX*(ISIM-1)),NORBT)
         IF (IPRAVE .GT.  80) THEN
            WRITE(LUPRI,'(//A,I4)') ' ISIM = ', ISIM
            WRITE(LUPRI,'(A)') ' Multiplication on ZYMAT'
            WRITE(LUPRI,'(A)') ' ======================='
            CALL OUTPUT(WRK(KGRA+N2ORBX*(ISIM-1)),1,NORBT,1,
     &                  NORBT,NORBT,NORBT,1,LUPRI)
         END IF
  200 CONTINUE
      CALL DSCAL(NSIM*N2ORBX,DHALF,WRK(KGRA),1)
C
C     and correct the gradient
C
      CALL RSPCFA(NSIM,WRK(KGRA),EVECS)
C
      CALL QEXIT('RFANTI')
C
      RETURN
      END
C  /* Deck rspant */
      SUBROUTINE RSPANT(FOCK,NORB,IOFF,NORBT,GMAT)
C
C Subroutine to compute antisymmetric part of Fock matrix
C
#include <implicit.h>
      DIMENSION FOCK(NORB,NORB), GMAT(NORBT,NORBT)
C
      CALL QENTER('RSPANT')
C
      DO 100 I = 1, NORB
         DO 200 J = 1, NORB
            GMAT(IOFF+I,IOFF+J) = FOCK(I,J) - FOCK(J,I)
            GMAT(IOFF+J,IOFF+I) = FOCK(J,I) - FOCK(I,J)
200      CONTINUE
100   CONTINUE
C
      CALL QEXIT('RSPANT')
C
      RETURN
      END
C  /* Deck rspcfa */
      SUBROUTINE RSPCFA(NSIM,GRA,EVECS)
#include <implicit.h>
C
      DIMENSION GRA(NORBT,NORBT,*), EVECS(KZYVAR,*)
C
#include <maxorb.h>
#include <maxash.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      CALL QENTER('RSPCFA')
C
      KYCONF = KZCONF + KZVAR
C
C Distribute gradient matrices into sigma vector in EVECS
C
      DO 100 ISIM = 1, NSIM
         DO 200 IG = 1,KZWOPT
            K     = JWOP(1,IG)
            L     = JWOP(2,IG)
            EVECS(KZCONF+IG,ISIM) = EVECS(KZCONF+IG,ISIM)
     *            + GRA(L,K,ISIM)
            EVECS(KYCONF+IG,ISIM) = EVECS(KYCONF+IG,ISIM)
     *            + GRA(K,L,ISIM)
 200    CONTINUE
 100  CONTINUE
C
      CALL QEXIT('RSPCFA')
C
      RETURN
      END
