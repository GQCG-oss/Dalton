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
C  /* Deck rlmset */
      SUBROUTINE RLMSET_CBISOL(LMAX)
C
C  17-Oct-2007 hjaaj
C
C  Set some /cbisol/ variables which are used in some calls,
C  also when DORLM is false.
C
#include <implicit.h>
C
C  CBISOL : LCAVMX,LMTOT,LMNTOT
C
#include <cbisol.h>
C
C     transfer information to CBISOL (for SOLNUC)
      LCAVMX = LMAX
      LMTOT  = (LCAVMX+1) ** 2
      LMNTOT = (LCAVMX+1)*(LCAVMX+2)*(LCAVMX+3) / 6
      RETURN
      END
C  /* Deck ireplm */
      FUNCTION IREPLM(L,M)
C
C     Symmetry of RLM integrals
C
C     tuh 01.02.90
C
#include "implicit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "symmet.h"
#include "ibtfun.h"
C
      IREPLM = 0
      IF (MOD(L + M,2) .EQ. 1)  IREPLM = ISYMAX(3,1)
      IF (MOD(ABS(M),2) .EQ .1) IREPLM = IBTXOR(IREPLM,ISYMAX(1,1))
      IF (M.LT.0) IREPLM=IBTXOR(IREPLM,IBTXOR(ISYMAX(1,1),ISYMAX(2,1)))
      RETURN
      END
C  /* Deck wrtsol */
      SUBROUTINE WRTSOL(RLMINT,L,NBAST,NNBASX,NOPTYP,INTREP,
     &                  SYMTOT,IPRINT)
C
C     22-Jan-1991 hjaaj
C     Write Solvent integrals to LUSOL
C
C  SYMTOT true: only totally-symmetric integrals written to LUSOL.
C
#include "implicit.h"
      DIMENSION RLMINT(NNBASX,NOPTYP), INTREP(NOPTYP)
      LOGICAL SYMTOT
#include "priunit.h"
C
C Used from common blocks:
C  INFTAP : LUSOL
C
#include <inftap.h>
C
      IF (NOPTYP .NE. 2*L + 1) THEN
         WRITE (LUPRI,*) 'WRTSOL ERROR: NOPTYP .ne. 2*L+1'
         CALL QUIT('WRTSOL ERROR')
      END IF
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(/A/A,L8/)') ' --- Output from WRTSOL ---',
     &      ' SYMTOT =',SYMTOT
      END IF
      DO 100 M = -L,L
         INTREP(L+1+M) = IREPLM(L,M) + 1
         IF (SYMTOT .AND. INTREP(L+1+M) .GT. 1) THEN
            INTREP(L+1+M) = -INTREP(L+1+M)
         END IF
         IF (IPRINT .GT. 2) WRITE (LUPRI,'(A,3I4)')
     &         ' l, m, symmetry of R(l,m)',L,M,INTREP(L+1+M)
  100 CONTINUE
      NOPTP8 = MAX(8,NOPTYP)
      WRITE (LUSOL) L,(INTREP(I),I=1,NOPTP8)
      DO 200 M = -L, L
      IF (INTREP(L+1+M) .GT. 0) THEN
C        Find R(l,m), they come in the order 0,1,-1,2,-2,...
         IF (M .GT. 0) THEN
            I = 2*M
         ELSE
            I = 2*(-M) + 1
         END IF
         IF (IPRINT .GT. 4) THEN
            WRITE (LUPRI,'(//A,2I4)') ' Integrals of operator R(l,m)'/
     &         /'; l,m =',L,M
            WRITE (LUPRI,'(A,I2)') ' Symmetry of operator:',
     &         INTREP(L+1+M)
            CALL OUTPAK(RLMINT(1,I),NBAST,1,LUPRI)
         END IF
         WRITE (LUSOL) (RLMINT(J,I), J = 1,NNBASX)
      END IF
  200 CONTINUE
      RETURN
      END
C  /* Deck rlmnuc */
      SUBROUTINE RLMNUC(WORK,LWORK,IPRINT)
C
C  9-Feb-1993 hjaaj
C
C interface to SOLNUC for Hermit
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
      DIMENSION WORK(LWORK)
C
C Used from common blocks:
C  SYMMET : MAXREP,NAOS()
C  CBISOL : LCAVMX,LMTOT,LMNTOT
C  ORGCOM : CAVORG
C  INFTAP : LUSOL
C
#include <symmet.h>
#include <cbisol.h>
#include <orgcom.h>
#include <inftap.h>
C
      CHARACTER*8 LABSOL(4)
      DIMENSION   FILL4(4)
      DATA LABSOL/'********',' HERMIT ',' HERMIT ','SOLVRLM '/,
     *     FILL4 /4*1.0D20/
C     Memory allocation
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
      CALL MEMGET('REAL',KTLMN,LMNTOT,WORK,KFREE,LFREE)
C     Calculate TLMN = RNUC in Cartesian basis
      CALL SOLNUC(.FALSE.,0,.FALSE.,WORK(KTLMN),DUMMY,DUMMY,WORK,KFREE,
     &            LFREE,IPRINT)
C     Transform TLMN to RNUC(l,m) in spherical basis and
C     write nuclear contribution, RNUC(l,m) on LUSOL
      LMTOT4 = MAX(4,LMTOT)
      CALL MEMGET('REAL',KRNUC,LMTOT4,WORK,KFREE,LFREE)
      CALL TLMTRA(1,WORK(KTLMN),WORK(KRNUC),WORK(KFREE),LFREE,IPRINT)
      IF (IPRINT .GT. 2) THEN
         CALL HEADER('Nuclear solvent integrals Tn(l,m)',-1)
         WRITE (LUPRI,'(1X,1P,4D18.8)') (WORK(I),I=KRNUC,KRNUC+LMTOT-1)
      END IF
C
C     ***** Label and info record on LUSOL *****
C     ***** Nuclear contributions          *****
C
      NBAST = 0
      DO 100 IREP = 1, MAXREP+1
         NBAST = NBAST + NAOS(IREP)
  100 CONTINUE
C
      CALL GETDAT(LABSOL(2),LABSOL(3))
      REWIND LUSOL
      WRITE (LUSOL) LABSOL
      WRITE (LUSOL) LCAVMX, (LCAVMX+1)*(LCAVMX+1), NBAST, CAVORG
      CALL WRITT(LUSOL,LMTOT4,WORK(KRNUC))
      CALL MEMREL('RLMNUC',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck solnuc */
      SUBROUTINE SOLNUC(PROPTY,MAXDIF,DIFDIP,TLMND,HSOLNN,FCM,
     &                  WORK,KFRSAV,LFRSAV,IPRINT)
C
C 12-Jan-1993 hjaaj+kvm
C
#include "implicit.h"
      DIMENSION TLMND(*),FCM(*), WORK(*), HSOLNN(*)
      LOGICAL   PROPTY, DIFDIP
C
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
C
C Used from common blocks:
C  CBISOL : LCAVMX,LMNTOT
C  NUCLEI : NUCDEP
C
#include "nuclei.h"
#include "cbisol.h"
C
      IF (IPRINT .GT. 4) THEN
         CALL HEADER('Output from SOLNUC',-1)
         WRITE (LUPRI,*) 'PROPTY ',PROPTY
         WRITE (LUPRI,*) 'MAXDIF ',MAXDIF
         WRITE (LUPRI,*) 'LFRSAV ',LFRSAV
         IF (PROPTY .AND. MAXDIF .GT. 0 .AND. IPRINT .GT. 8) THEN
            WRITE (LUPRI,*) 'FCM factors'
            IOFF = 0
            DO 100 LCAV = 0,LCAVMX
               J = (LCAV+1)*(LCAV+2)/2
               WRITE (LUPRI,*) '   L_cav =',LCAV
               WRITE (LUPRI,'(5F15.8)') (FCM(IOFF+I),I=1,J)
               IOFF = IOFF + J
  100       CONTINUE
         END IF
      END IF
C
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      CALL MEMGET('REAL',KDX0,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDY0,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDZ0,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLO,LMNTOT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KMO,LMNTOT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNO,LMNTOT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSCTRA,MXCOOR*MXCOOR,WORK,KFREE,LFREE)
      LWRK  = LFRSAV - KFREE + 1
      CALL SOLNC0(TLMND,FCM,WORK(KLO),WORK(KMO),WORK(KNO),
     &            WORK(KDX0),WORK(KDY0),WORK(KDZ0),
     &            HSOLNN,WORK(KCSTRA),WORK(KSCTRA),WORK(KFREE),LWRK,
     &            PROPTY,MAXDIF,DIFDIP,IPRINT)
C
      CALL MEMREL('SOLNUC',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck solnc0 */
      SUBROUTINE SOLNC0(TLMND,FCM,LO,MO,NO,DX0,DY0,DZ0,HSOLNN,
     &                  CSTRA,SCTRA,WORK,LWORK,PROPTY,MAXDIF,DIFDIP,
     &                  IPRINT)
C
C     12. Jan. 1993 HJAaJ+KVM
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
      DIMENSION TLMND(LMNTOT,*), FCM(LMNTOT), WORK(LWORK)
      DIMENSION LO(*), MO(*), NO(*), HSOLNN(MXCOOR,MXCOOR)
      DIMENSION DX0(*), DY0(*), DZ0(*), CSTRA(*), SCTRA(*)
      LOGICAL   PROPTY, DIFDIP
C
      PARAMETER (D0 = 0.0 D00, D1 = 1.0 D00, D2 = 2.0 D00 )
      PARAMETER (DP5 = 0.5 D00, DP25 = 0.25 D00 )
C
C Used from common blocks:
C  CBISOL : LCAVMX,NCNTCV
C  TAYSOL : GSOLNN()
C  ORGCOM : CAVORG(3)
C
#include "cbisol.h"
#include "nuclei.h"
#include "symmet.h"
#include "taysol.h"
#include "orgcom.h"
C
C     Local arrays
C
C
C     Statement functions
C
#include "ibtfun.h"
C
C     Cartesian orders
C     ----------------
C
      CALL RLMSET(LCAVMX,LO,MO,NO)
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Test output from SOLNUC',-1)
         WRITE (LUPRI,*) 'LMN index, LO,MO,NO'
         DO 10 I = 1,LMNTOT
   10       WRITE (LUPRI,'(5I10)') I,LO(I),MO(I),NO(I)
      END IF
C
C     Find cavity coordinates and indices
C     ===================================
C
      PXCAV = CAVORG(1)
      PYCAV = CAVORG(2)
      PZCAV = CAVORG(3)
      IF (.NOT.PROPTY) THEN
C        ... Hermit
         CALL DZERO(TLMND,LMNTOT)
      ELSE
         NERR = 0
         IF (NCNTCV .NE. NUCIND) NERR = NERR + 1
         IF (CORD(1,NUCIND) .NE. CAVORG(1)) NERR = NERR + 1
         IF (CORD(2,NUCIND) .NE. CAVORG(2)) NERR = NERR + 1
         IF (CORD(3,NUCIND) .NE. CAVORG(3)) NERR = NERR + 1
         IF (NERR .GT. 0) THEN
            WRITE (LUPRI,*) 'SOLNUC errors for ABACUS'
            WRITE (LUPRI,*) 'NCNTCV,NUCIND',NCNTCV,NUCIND
            WRITE (LUPRI,*) 'CAVORG(1:3)     ',CAVORG
            WRITE (LUPRI,*) 'CORD(1:3,NUCIND)',
     &         CORD(1,NUCIND),CORD(2,NUCIND),CORD(3,NUCIND)
            CALL QUIT('SOLNUC errors for ABACUS')
         END IF
         NXCAV = 3*NUCDEP - 2
         NYCAV = 3*NUCDEP - 1
         NZCAV = 3*NUCDEP
         IF (MAXDIF .GE. 2 .OR. DIFDIP) 
     &        CALL DZERO(TLMND,3*NUCDEP*LMNTOT)
      END IF
C
C     loop over atoms
C     ===============
C
      IATOM = 0
      DO 720 ICENT = 1, NUCIND
         IF (PROPTY .AND. ICENT .EQ. NUCIND) GO TO 720
C        ... skip cavity center
         MULCNT = ISTBNU(ICENT)
      DO 700 ISYMOP = 0, MAXOPR
      IF (IBTAND(ISYMOP,MULCNT) .NE. 0) GO TO 700
         IATOM = IATOM + 1
      IF (CHARGE(ICENT) .EQ. D0) GO TO 700
         PX    = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
         PY    = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
         PZ    = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
         PXDIF = PX - PXCAV
         PYDIF = PY - PYCAV
         PZDIF = PZ - PZCAV
         NXA   = 3*IATOM - 2
         NYA   = 3*IATOM - 1
         NZA   = 3*IATOM
C
         IF (IPRINT .GE. 10) THEN
            WRITE(LUPRI,'(/A,I3,/A,3F15.8)')
     *         ' >>> SOLNUC test print for atom',IATOM,
     *         '     with coordinates :',PX,PY,PZ
         END IF
         CALL SOLNC1(PROPTY,MAXDIF,DIFDIP,IATOM,CHARGE(ICENT),
     &               PXDIF,PYDIF,PZDIF,TLMND,
     &               DX0,DY0,DZ0, LO,MO,NO, FCM, HSOLNN,IPRINT)
  700 CONTINUE
  720 CONTINUE
C
      IF (PROPTY) THEN
         CALL TRACOR(CSTRA,SCTRA,1,MXCOOR,0)
         IF (IATOM .NE. NUCDEP-1) THEN
            WRITE (LUPRI,*) 'SOLNUC error; IATOM .ne. NUCDEP-1'
            WRITE (LUPRI,*) 'IATOM, NUCDEP',IATOM,NUCDEP
            CALL QUIT('SOLNUC error: counted # of atoms .ne. NUCDEP')
         END IF
         IF (MAXDIF .GE. 1) THEN
            GSOLNN(NXCAV) = - DSUM(IATOM,GSOLNN(1),3)
            GSOLNN(NYCAV) = - DSUM(IATOM,GSOLNN(2),3)
            GSOLNN(NZCAV) = - DSUM(IATOM,GSOLNN(3),3)
            CALL TRACTS(GSOLNN,3*NUCDEP,CSTRA)
C           CALL TRACTS(VEC,NCOOR,CSTRA)
         END IF
         IF (MAXDIF .GE. 2) THEN
            DO 800 I = 1, 3*IATOM
               HSOLNN(NXCAV,I) = - DSUM(IATOM,HSOLNN(1,I),3)
               HSOLNN(NYCAV,I) = - DSUM(IATOM,HSOLNN(2,I),3)
               HSOLNN(NZCAV,I) = - DSUM(IATOM,HSOLNN(3,I),3)
               HSOLNN(I,NXCAV) = HSOLNN(NXCAV,I)
               HSOLNN(I,NYCAV) = HSOLNN(NYCAV,I)
               HSOLNN(I,NZCAV) = HSOLNN(NZCAV,I)
  800       CONTINUE
            HSOLNN(NXCAV,NXCAV) = - DSUM(IATOM,HSOLNN(1,NXCAV),3)
            HSOLNN(NXCAV,NYCAV) = - DSUM(IATOM,HSOLNN(1,NYCAV),3)
            HSOLNN(NXCAV,NZCAV) = - DSUM(IATOM,HSOLNN(1,NZCAV),3)
            HSOLNN(NYCAV,NYCAV) = - DSUM(IATOM,HSOLNN(2,NYCAV),3)
            HSOLNN(NYCAV,NZCAV) = - DSUM(IATOM,HSOLNN(2,NZCAV),3)
            HSOLNN(NZCAV,NZCAV) = - DSUM(IATOM,HSOLNN(3,NZCAV),3)
            HSOLNN(NYCAV,NXCAV) = HSOLNN(NXCAV,NYCAV)
            HSOLNN(NZCAV,NXCAV) = HSOLNN(NXCAV,NZCAV)
            HSOLNN(NZCAV,NYCAV) = HSOLNN(NYCAV,NZCAV)
         END IF
         IF (MAXDIF .GE. 2 .OR. DIFDIP) THEN
            DO 810 I = 1,LMNTOT
               TLMND(I,NXCAV) = - DSUM(IATOM,TLMND(I,1),LMNTOT*3)
               TLMND(I,NYCAV) = - DSUM(IATOM,TLMND(I,2),LMNTOT*3)
               TLMND(I,NZCAV) = - DSUM(IATOM,TLMND(I,3),LMNTOT*3)
  810       CONTINUE
            IF (MAXDIF .GE. 2) THEN
               CALL TRHCTS(HSOLNN,3*NUCDEP,MAXDIF,TLMND,LMNTOT,
     &                     CSTRA,WORK,LWORK)
            ELSE
               CALL TRHCTS(DUMMY,3*NUCDEP,MAXDIF,TLMND,LMNTOT,
     &                     CSTRA,WORK,LWORK)
            END IF
         END IF
C
         IF (IPRINT .GE. 2 .AND. MAXDIF .GE.1) THEN
            CALL HEADER('Nuclear contributions to solvent gradient',-1)
            CALL PRIGRD(GSOLNN,CSTRA,SCTRA)
         END IF
         IF (IPRINT .GE. 2 .AND. MAXDIF .GE. 2) THEN
            CALL HEADER('Nuclear contributions to solvent Hessian',-1)
            CALL PRIHES(HSOLNN,'CENTERS',CSTRA,SCTRA)
         END IF
      END IF
C
C     End of SOLNC0.
C
      RETURN
      END
C  /* Deck tracts */
      SUBROUTINE TRACTS(VEC,NCOOR,CSTRA)
C
C 12. Jan. 1993 hjaaj+kvm (based on TRAGRD)
C
C Transform Cartesian vector (e.g. gradient) to
C           symmetry  vector
C
#include "implicit.h"
#include "mxcent.h"
      DIMENSION VEC(MXCOOR), CSTRA(MXCOOR,MXCOOR)
      DIMENSION SVEC(MXCOOR)
C
      CALL DGEMM('N','N',NCOOR,1,NCOOR,1.D0,
     &           CSTRA(1,1),MXCOOR,
     &           VEC,MXCOOR,0.D0,
     &           SVEC,MXCOOR)
      CALL DCOPY(NCOOR,SVEC,1,VEC,1)
      RETURN
      END
C  /* Deck trhcts */
      SUBROUTINE TRHCTS(HMAT,NCOOR,MAXDIF,TLMND,LMNTOT,CSTRA,WORK,LWORK)
C
C tuh - based on TRACTS
C
C Transform Cartesian matrix (e.g. Hessian) to
C           symmetry  vector
C
#include "implicit.h"
#include "mxcent.h"
      DIMENSION HMAT(MXCOOR,MXCOOR), TLMND(LMNTOT,NCOOR), 
     &          CSTRA(MXCOOR,MXCOOR), WORK(LWORK)
      NDIM = MAX(MXCOOR*MXCOOR,LMNTOT*NCOOR)
      IF (LWORK .LT. NDIM) CALL STOPIT('TRHCTS',' ',LWORK,NDIM)
C
C     Transform Hessian
C
      IF (MAXDIF .GE. 2) THEN
         CALL DGEMM('N','N',NCOOR,NCOOR,NCOOR,1.D0,
     &              CSTRA(1,1),MXCOOR,
     &              HMAT,MXCOOR,0.D0,
     &              WORK,MXCOOR)
         CALL DGEMM('N','T',NCOOR,NCOOR,NCOOR,1.D0,
     &              WORK,MXCOOR,
     &              CSTRA(1,1),MXCOOR,0.D0,
     &              HMAT,MXCOOR)
      END IF
C
C     Transform TLMND
C
      CALL DGEMM('N','T',LMNTOT,NCOOR,NCOOR,1.D0,
     &           TLMND,LMNTOT,
     &           CSTRA(1,1),MXCOOR,0.D0,
     &           WORK,LMNTOT)
      CALL DCOPY(LMNTOT*NCOOR,WORK,1,TLMND,1)
      RETURN
      END
C  /* Deck solnc1 */
      SUBROUTINE SOLNC1(PROPTY,MAXDIF,DIFDIP,IATOM,CHRG,PXDIF,PYDIF,
     &                  PZDIF,TLMND,DX0,DY0,DZ0,LO,MO,NO,
     &                  FCM,HSOLNN,IPRINT)
C
C              FCM(k) = -2 sum(lm) g(l) T(lm) C_k(lm)
C                 where C_k(lm) transforms from Cartesian moments
C                 to spherical moments.
C
C from common:
C              LMNTOT  = (LCAVMX+1)*(LCAVMX+2)*(LCAVMX+3)/6
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DMP5 = -0.5D0)
C
      DIMENSION DX0(-2:LCAVMX), DY0(-2:LCAVMX), DZ0(-2:LCAVMX)
      DIMENSION LO(*), MO(*), NO(*), HSOLNN(MXCOOR,MXCOOR)
      DIMENSION TLMND(LMNTOT,*), FCM(LMNTOT), HSOL(6)
      LOGICAL   PROPTY, DIFDIP
C
C Used from common blocks:
C  CBISOL : LCAVMX,LMNTOT,?
C  LMNS   : [LMN]VALU[AB](MXAQN)
C  TAYSOL : GSOLNN(), HSOLNN()
C
#include "cbisol.h"
#include "lmns.h"
#include "taysol.h"
C
C     Statement functions
C     -------------------
C
      DX1(I) = I * DX0(I-1)
      DY1(I) = I * DY0(I-1)
      DZ1(I) = I * DZ0(I-1)
      DX2(I) = I * (I-1) * DX0(I-2)
      DY2(I) = I * (I-1) * DY0(I-2)
      DZ2(I) = I * (I-1) * DZ0(I-2)
C
C
C     Cartesian integrals
C     -------------------
C
C     Note: DX0,DY0,DZ0 starts at -2
C           in order to avoid "if" tests in construction of
C           three-dimensional integrals from one-dimensional integrals.
C           These values must thus be initialized to zero.
C
C        One-dimensional integrals
C
         CALL DZERO(DX0(-2),LCAVMX+3)
         CALL DZERO(DY0(-2),LCAVMX+3)
         CALL DZERO(DZ0(-2),LCAVMX+3)
         DX0(0) = D1
         DY0(0) = D1
         DZ0(0) = D1
CDIR$    NORECURRENCE
         DO 400 L = 1, LCAVMX
            DX0(L) = PXDIF ** L
            DY0(L) = PYDIF ** L
            DZ0(L) = PZDIF ** L
  400    CONTINUE
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,*) 'Test output from SOLNC1'
            WRITE (LUPRI,*) 'IATOM, CHRG',IATOM,CHRG
            WRITE (LUPRI,*) 'P?DIF ',PXDIF,PYDIF,PZDIF
            WRITE (LUPRI,*) 'L, DX0, DY0, DZ0'
            WRITE (LUPRI,'(I10,3F20.10)')
     &         (I,DX0(I),DY0(I),DZ0(I),I=0,LCAVMX)
         END IF
C
C        Three-dimensional integrals
C
         IF (.NOT. PROPTY) THEN
C           solvent TN(l,m,n) for Hermit
            DO 500 I = 1, LMNTOT
               TLMND(I,1) = TLMND(I,1)
     &                    + CHRG*DX0(LO(I))*DY0(MO(I))*DZ0(NO(I))
  500       CONTINUE
            GO TO 9999
         END IF
C
         NXA = 3*IATOM - 2
         NYA = 3*IATOM - 1
         NZA = 3*IATOM
         ETMPNN = D0
         DO 600 I = 1, LMNTOT
C           *  undifferentiated integrals
            ETMPNN = ETMPNN + FCM(I)*DX0(LO(I))*DY0(MO(I))*DZ0(NO(I))
            IF (MAXDIF .GE. 1) THEN
C           *  0X00,0Y00,0Z00 terms
               D0X00 = CHRG*DX1(LO(I))*DY0(MO(I))*DZ0(NO(I))
               D0Y00 = CHRG*DX0(LO(I))*DY1(MO(I))*DZ0(NO(I))
               D0Z00 = CHRG*DX0(LO(I))*DY0(MO(I))*DZ1(NO(I))
               GSOLNN(NXA) = GSOLNN(NXA) - FCM(I)*D0X00
               GSOLNN(NYA) = GSOLNN(NYA) - FCM(I)*D0Y00
               GSOLNN(NZA) = GSOLNN(NZA) - FCM(I)*D0Z00
               IF (MAXDIF .GE. 2 .OR. DIFDIP) THEN
                  TLMND(I,NXA) = TLMND(I,NXA) + D0X00
                  TLMND(I,NYA) = TLMND(I,NYA) + D0Y00
                  TLMND(I,NZA) = TLMND(I,NZA) + D0Z00
               END IF
            END IF
  600    CONTINUE
         ESOLNN = ESOLNN + DMP5 * CHRG * ETMPNN
         IF (MAXDIF .GE. 2) THEN
            CALL DZERO(HSOL,6)
            DO 620 I = 1, LMNTOT
               HSOL(1) = HSOL(1)-FCM(I)*DX2(LO(I))*DY0(MO(I))*DZ0(NO(I))
               HSOL(2) = HSOL(2)-FCM(I)*DX1(LO(I))*DY1(MO(I))*DZ0(NO(I))
               HSOL(3) = HSOL(3)-FCM(I)*DX1(LO(I))*DY0(MO(I))*DZ1(NO(I))
               HSOL(4) = HSOL(4)-FCM(I)*DX0(LO(I))*DY2(MO(I))*DZ0(NO(I))
               HSOL(5) = HSOL(5)-FCM(I)*DX0(LO(I))*DY1(MO(I))*DZ1(NO(I))
               HSOL(6) = HSOL(6)-FCM(I)*DX0(LO(I))*DY0(MO(I))*DZ2(NO(I))
  620       CONTINUE
            HSOLNN(NXA,NXA) = HSOLNN(NXA,NXA) + CHRG*HSOL(1)
            HSOLNN(NXA,NYA) = HSOLNN(NXA,NYA) + CHRG*HSOL(2)
            HSOLNN(NXA,NZA) = HSOLNN(NXA,NZA) + CHRG*HSOL(3)
            HSOLNN(NYA,NXA) = HSOLNN(NYA,NXA) + CHRG*HSOL(2)
            HSOLNN(NYA,NYA) = HSOLNN(NYA,NYA) + CHRG*HSOL(4)
            HSOLNN(NYA,NZA) = HSOLNN(NYA,NZA) + CHRG*HSOL(5)
            HSOLNN(NZA,NXA) = HSOLNN(NZA,NXA) + CHRG*HSOL(3)
            HSOLNN(NZA,NYA) = HSOLNN(NZA,NYA) + CHRG*HSOL(5)
            HSOLNN(NZA,NZA) = HSOLNN(NZA,NZA) + CHRG*HSOL(6)
         END IF
C
 9999 RETURN
      END
#ifdef UNDEF
/* Comdeck to_do */
921013-kvm+hjaaj
  flyt RLMSET udenfor (?)
#endif
#ifdef UNDEF
/* Comdeck log */
Last revision: 27-May-1992 kvm+hjaaj
               10-Dec-1992 th+kvm+hjaaj
#endif
C  /* Deck solint */
      SUBROUTINE SOLINT(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,
     &                  RLMINT,RLMTAB,CORPX,CORPY,CORPZ,MAXDIF,
     &                  EXPPI,FCM,WORK,LWORK,IPRINT)
C
C 27-May-1992 kvm+hjaaj
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      DIMENSION RLMINT(*), RLMTAB(*), FCM(*), WORK(LWORK)
C
C Used from common blocks:
C  ONECOM : NHKT[AB]
C  CBISOL : LCAVMX,LMNTOT
C
#include "onecom.h"
#include "cbisol.h"
C
      ITMAX = NHKTA-1 + NHKTB-1 + MAXDIF
      ITMAX = MIN(LCAVMX,ITMAX)
      IF (IPRINT .GT. 4) THEN
         CALL HEADER('Output from SOLINT',-1)
         WRITE (LUPRI,*) 'MAXDIF',MAXDIF
         WRITE (LUPRI,'(A,3F20.12)') 'CORP? ',CORPX,CORPY,CORPZ
         WRITE (LUPRI,*) 'LWORK ',LWORK
         WRITE (LUPRI,*) 'EXPPI ',EXPPI
         WRITE (LUPRI,*) 'JMAX? ',JMAXA,JMAXB,JMAXT,JMAXD,JMAXM
         IF (MAXDIF .GT. 0 .AND. IPRINT .GT. 8) THEN
            WRITE (LUPRI,*) 'FCM factors'
            IOFF = 0
            DO 100 LCAV = 0,LCAVMX
               J = (LCAV+1)*(LCAV+2)/2
               WRITE (LUPRI,*) '   L_cav =',LCAV
               WRITE (LUPRI,'(5F15.8)') (FCM(IOFF+I),I=1,J)
               IOFF = IOFF + J
  100       CONTINUE
         END IF
      END IF
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
      CALL MEMGET('REAL',KHMU,(LCAVMX+1)*(LCAVMX+1)*3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDX0,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDY0,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDZ0,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDX1,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDY1,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDZ1,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDX2,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDY2,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDZ2,LCAVMX+3,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLO ,LMNTOT  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KMO ,LMNTOT  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNO ,LMNTOT  ,WORK,KFREE,LFREE)
C     IF (KFREE .GT. LWORK) CALL STOPIT('SOLINT',' ',LWORK,KFREE)
      CALL SOLIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,MAXDIF,
     &            RLMINT,RLMTAB,CORPX,CORPY,CORPZ,EXPPI,
     &            WORK(KDX0),WORK(KDY0),WORK(KDZ0),
     &            WORK(KDX1),WORK(KDY1),WORK(KDZ1),
     &            WORK(KDX2),WORK(KDY2),WORK(KDZ2),
     &            WORK(KLO), WORK(KMO), WORK(KNO),
     &            WORK(KHMU),ITMAX,FCM, IPRINT)
      CALL MEMREL('SOLINT',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck solin1 */
      SUBROUTINE SOLIN1(ODC,JMAXA,JMAXB,JMAXT,JMAXD,JMAXM,MAXDIF,
     &                  RLMINT,RLMTAB,CORPX,CORPY,CORPZ,EXPPI,
     &                  DX0,DY0,DZ0,DX1,DY1,DZ1,DX2,DY2,DZ2,LO,MO,NO,
     &                  HMULT,ITMAX,FCM,IPRINT)
C
Chjaaj-920527: RLMINT(KCKTAB,LMNTOT,7)
C                RLMINT(*,*,1)   : undifferentiated integrals
C                RLMINT(*,*,2:7) : first derivative integrals
C              RLMTAB(KCKTAB,21)  : second derivative expval contr.
C              LMNTOT  = (LCAVMX+1)*(LCAVMX+2)*(LCAVMX+3)/6
C              FCM(k) = -2 sum(lm) g(l) T(lm) C_k(lm)
C                 where C_k(lm) transforms from Cartesian moments
C                 to spherical moments.
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.0D0)
C
      DIMENSION ODC(0:JMAXA,0:JMAXB,0:JMAXT,0:JMAXD,0:JMAXM,3)
      DIMENSION DX0(-2:LCAVMX), DY0(-2:LCAVMX), DZ0(-2:LCAVMX),
     &          DX1(-2:LCAVMX), DY1(-2:LCAVMX), DZ1(-2:LCAVMX),
     &          DX2(-2:LCAVMX), DY2(-2:LCAVMX), DZ2(-2:LCAVMX)
      DIMENSION LO(*), MO(*), NO(*)
      DIMENSION HMULT(0:LCAVMX,0:LCAVMX,3)
      DIMENSION RLMINT(KCKTAB,LMNTOT,7), RLMTAB(KCKTAB,21), FCM(LMNTOT)
C
C Used from common blocks:
C  ONECOM : KCKTA, KCKTB, KCKTAB,
C           (KCKTAB = KCKTA*KCKTB)
C  LMNS   : [LMN]VALU[AB](MXAQN)
C
#include "cbisol.h"
#include "onecom.h"
#include "lmns.h"
#include "orgcom.h"
C
      INTEGER T,U,V
C
C     Cartesian orders
C     ----------------
C
      CALL RLMSET(LCAVMX,LO,MO,NO)
C
C     Hermitian integrals
C     -------------------
C
      CALL MOMHER(HMULT,ITMAX,LCAVMX,CAVORG,CORPX,CORPY,CORPZ,EXPPI)
C
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Test output from SOLIN1',-1)
         WRITE (LUPRI,*) 'KCKT? ',KCKTA,KCKTB
         WRITE (LUPRI,*) 'CAVORG(1:3)',(CAVORG(I),I=1,3)
         WRITE (LUPRI,*) 'LMN index, LO,MO,NO'
         DO 10 I = 1,LMNTOT
   10       WRITE (LUPRI,'(5I10)') I,LO(I),MO(I),NO(I)
         WRITE (LUPRI,*) 'HMULT(0:ITMAX,0:LCAVMX,1) matrix (X)'
         CALL OUTPUT(HMULT(0,0,1),1,ITMAX+1,1,LCAVMX+1,
     &               LCAVMX+1,LCAVMX+1,1,LUPRI)
         WRITE (LUPRI,*) 'HMULT(0:ITMAX,0:LCAVMX,2) matrix (Y)'
         CALL OUTPUT(HMULT(0,0,2),1,ITMAX+1,1,LCAVMX+1,
     &               LCAVMX+1,LCAVMX+1,1,LUPRI)
         WRITE (LUPRI,*) 'HMULT(0:ITMAX,0:LCAVMX,3) matrix (Z)'
         CALL OUTPUT(HMULT(0,0,3),1,ITMAX+1,1,LCAVMX+1,
     &               LCAVMX+1,LCAVMX+1,1,LUPRI)
      END IF
C
C     Cartesian integrals
C     -------------------
C
C     Note: DX0,DY0,DZ0 starts at -2 and DZ1,DY1,DZ1 starts at -1
C           in order to avoid "if" tests in construction of
C           three-dimensional integrals from one-dimensional integrals.
C           These values must thus be initialized to zero.
C
      IRADR   = 0
      DO 300 ICOMPA = 1,KCKTA
         LVALA = LVALUA(ICOMPA)
         MVALA = MVALUA(ICOMPA)
         NVALA = NVALUA(ICOMPA)
      DO 300 ICOMPB = 1, KCKTB
         LVALB = LVALUB(ICOMPB)
         MVALB = MVALUB(ICOMPB)
         NVALB = NVALUB(ICOMPB)
C
C        One-dimensional integrals
C
         CALL DZERO(DX0(-2),LCAVMX+3)
         CALL DZERO(DY0(-2),LCAVMX+3)
         CALL DZERO(DZ0(-2),LCAVMX+3)
         IF (MAXDIF .GE. 1) THEN
            CALL DZERO(DX1(-2),LCAVMX+3)
            CALL DZERO(DY1(-2),LCAVMX+3)
            CALL DZERO(DZ1(-2),LCAVMX+3)
         END IF
         IF (MAXDIF .GE. 2) THEN
           CALL DZERO(DX2(-2),LCAVMX+3)
           CALL DZERO(DY2(-2),LCAVMX+3)
           CALL DZERO(DZ2(-2),LCAVMX+3)
         END IF
         DO 400 L = 0, LCAVMX
            DO 500 T = 0, MIN(LVALA+LVALB,L)
              DX0(L) = DX0(L) + ODC(LVALA,LVALB,T,0,0,1)*HMULT(T,L,1)
  500       CONTINUE
            DO 510 U = 0, MIN(MVALA+MVALB,L)
              DY0(L) = DY0(L) + ODC(MVALA,MVALB,U,0,0,2)*HMULT(U,L,2)
  510       CONTINUE
            DO 520 V = 0, MIN(NVALA+NVALB,L)
              DZ0(L) = DZ0(L) + ODC(NVALA,NVALB,V,0,0,3)*HMULT(V,L,3)
  520       CONTINUE
            IF (MAXDIF .GE. 1) THEN
               DO 501 T = 0, MIN(LVALA+LVALB+1,L)
                 DX1(L) = DX1(L) + ODC(LVALA,LVALB,T,1,0,1)*HMULT(T,L,1)
  501          CONTINUE
               DO 511 U = 0, MIN(MVALA+MVALB+1,L)
                 DY1(L) = DY1(L) + ODC(MVALA,MVALB,U,1,0,2)*HMULT(U,L,2)
  511          CONTINUE
               DO 521 V = 0, MIN(NVALA+NVALB+1,L)
                 DZ1(L) = DZ1(L) + ODC(NVALA,NVALB,V,1,0,3)*HMULT(V,L,3)
  521          CONTINUE
            END IF
            IF (MAXDIF .GE. 2) THEN
               DO 502 T = 0, MIN(LVALA+LVALB+2,L)
                 DX2(L) = DX2(L) + ODC(LVALA,LVALB,T,2,0,1)*HMULT(T,L,1)
  502          CONTINUE
               DO 512 U = 0, MIN(MVALA+MVALB+2,L)
                 DY2(L) = DY2(L) + ODC(MVALA,MVALB,U,2,0,2)*HMULT(U,L,2)
  512          CONTINUE
               DO 522 V = 0, MIN(NVALA+NVALB+2,L)
                 DZ2(L) = DZ2(L) + ODC(NVALA,NVALB,V,2,0,3)*HMULT(V,L,3)
  522          CONTINUE
            END IF
  400    CONTINUE
         IF (IPRINT .GT. 10) THEN
            WRITE (LUPRI,*) 'ICOMP? ',ICOMPA,ICOMPB
            WRITE (LUPRI,*) '[LMN]VALA',LVALA,MVALA,NVALA
            WRITE (LUPRI,*) '[LMN]VALB',LVALB,MVALB,NVALB
            WRITE (LUPRI,*) 'L, DX0, DX1, DX2'
            WRITE (LUPRI,'(I10,3F20.10)')
     &         (I,DX0(I),DX1(I),DX2(I),I=0,LCAVMX)
            WRITE (LUPRI,*) 'M, DY0, DY1, DY2'
            WRITE (LUPRI,'(I10,3F20.10)')
     &         (I,DY0(I),DY1(I),DY2(I),I=0,LCAVMX)
            WRITE (LUPRI,*) 'N, DZ0, DZ1, DZ2'
            WRITE (LUPRI,'(I10,3F20.10)')
     &         (I,DZ0(I),DZ1(I),DZ2(I),I=0,LCAVMX)
         END IF
C
C        Three-dimensional integrals
C
         IRADR    = IRADR    + 1
         DO 600 I = 1, LMNTOT
C           *  undifferentiated integrals
            RLMINT(IRADR,I,1) = RLMINT(IRADR,I,1)
     +                       + DX0(LO(I))*DY0(MO(I))*DZ0(NO(I))
            IF (MAXDIF .GE. 1) THEN
C           *  0X00,0Y00,0Z00 terms
               RLMINT(IRADR,I,2) = RLMINT(IRADR,I,2)
     &                        + DX1(LO(I))*DY0(MO(I))*DZ0(NO(I))
               RLMINT(IRADR,I,3) = RLMINT(IRADR,I,3)
     &                        + DX0(LO(I))*DY1(MO(I))*DZ0(NO(I))
               RLMINT(IRADR,I,4) = RLMINT(IRADR,I,4)
     &                        + DX0(LO(I))*DY0(MO(I))*DZ1(NO(I))
C           *  000X,000Y,000Z terms
               RLMINT(IRADR,I,5) = RLMINT(IRADR,I,5)
     &                        - LO(I)*DX0(LO(I)-1)*DY0(MO(I))*DZ0(NO(I))
               RLMINT(IRADR,I,6) = RLMINT(IRADR,I,6)
     &                        - MO(I)*DX0(LO(I))*DY0(MO(I)-1)*DZ0(NO(I))
               RLMINT(IRADR,I,7) = RLMINT(IRADR,I,7)
     &                        - NO(I)*DX0(LO(I))*DY0(MO(I))*DZ0(NO(I)-1)
            END IF
  600    CONTINUE
         IF (MAXDIF .GE. 2) THEN
            DO 620 I = 1, LMNTOT
C           *  XX00,XY00,XZ00,YY00,YZ00,ZZ00 terms
               RLMTAB(IRADR,1) = RLMTAB(IRADR,1)
     &        +FCM(I)*DX2(LO(I))*DY0(MO(I))*DZ0(NO(I))
               RLMTAB(IRADR,2) = RLMTAB(IRADR,2)
     &        +FCM(I)*DX1(LO(I))*DY1(MO(I))*DZ0(NO(I))
               RLMTAB(IRADR,3) = RLMTAB(IRADR,3)
     &        +FCM(I)*DX1(LO(I))*DY0(MO(I))*DZ1(NO(I))
               RLMTAB(IRADR,4) = RLMTAB(IRADR,4)
     &        +FCM(I)*DX0(LO(I))*DY2(MO(I))*DZ0(NO(I))
               RLMTAB(IRADR,5) = RLMTAB(IRADR,5)
     &        +FCM(I)*DX0(LO(I))*DY1(MO(I))*DZ1(NO(I))
               RLMTAB(IRADR,6) = RLMTAB(IRADR,6)
     &        +FCM(I)*DX0(LO(I))*DY0(MO(I))*DZ2(NO(I))
C           *  0X0X,0X0Y,0X0Z terms
               RLMTAB(IRADR,7) = RLMTAB(IRADR,7)
     &        -FCM(I)*LO(I)*DX1(LO(I)-1)*DY0(MO(I))*DZ0(NO(I))
               RLMTAB(IRADR,8) = RLMTAB(IRADR,8)
     &        -FCM(I)*MO(I)*DX1(LO(I))*DY0(MO(I)-1)*DZ0(NO(I))
               RLMTAB(IRADR,9) = RLMTAB(IRADR,9)
     &        -FCM(I)*NO(I)*DX1(LO(I))*DY0(MO(I))*DZ0(NO(I)-1)
C           *  0Y0X,0Y0Y,0Y0Z terms
               RLMTAB(IRADR,10) = RLMTAB(IRADR,10)
     &        -FCM(I)*LO(I)*DX0(LO(I)-1)*DY1(MO(I))*DZ0(NO(I))
               RLMTAB(IRADR,11) = RLMTAB(IRADR,11)
     &        -FCM(I)*MO(I)*DX0(LO(I))*DY1(MO(I)-1)*DZ0(NO(I))
               RLMTAB(IRADR,12) = RLMTAB(IRADR,12)
     &        -FCM(I)*NO(I)*DX0(LO(I))*DY1(MO(I))*DZ0(NO(I)-1)
C           *  0Z0X,0Z0Y,0Z0Z terms
               RLMTAB(IRADR,13) = RLMTAB(IRADR,13)
     &        -FCM(I)*LO(I)*DX0(LO(I)-1)*DY0(MO(I))*DZ1(NO(I))
               RLMTAB(IRADR,14) = RLMTAB(IRADR,14)
     &        -FCM(I)*MO(I)*DX0(LO(I))*DY0(MO(I)-1)*DZ1(NO(I))
               RLMTAB(IRADR,15) = RLMTAB(IRADR,15)
     &        -FCM(I)*NO(I)*DX0(LO(I))*DY0(MO(I))*DZ1(NO(I)-1)
C           *  00XX,00XY,00XZ,00YY,00YZ,00ZZ terms
               RLMTAB(IRADR,16) = RLMTAB(IRADR,16)
     &        +FCM(I)*LO(I)*(LO(I)-1)*DX0(LO(I)-2)*DY0(MO(I))*DZ0(NO(I))
               RLMTAB(IRADR,17) = RLMTAB(IRADR,17)
     &        +FCM(I)*LO(I)*MO(I)*DX0(LO(I)-1)*DY0(MO(I)-1)*DZ0(NO(I))
               RLMTAB(IRADR,18) = RLMTAB(IRADR,18)
     &        +FCM(I)*LO(I)*NO(I)*DX0(LO(I)-1)*DY0(MO(I))*DZ0(NO(I)-1)
               RLMTAB(IRADR,19) = RLMTAB(IRADR,19)
     &        +FCM(I)*MO(I)*(MO(I)-1)*DX0(LO(I))*DY0(MO(I)-2)*DZ0(NO(I))
               RLMTAB(IRADR,20) = RLMTAB(IRADR,20)
     &        +FCM(I)*MO(I)*NO(I)*DX0(LO(I))*DY0(MO(I)-1)*DZ0(NO(I)-1)
               RLMTAB(IRADR,21) = RLMTAB(IRADR,21)
     &        +FCM(I)*NO(I)*(NO(I)-1)*DX0(LO(I))*DY0(MO(I))*DZ0(NO(I)-2)
  620       CONTINUE
         END IF
  300 CONTINUE
      RETURN
      END
C  /* Deck momher */
      SUBROUTINE MOMHER(HMULT,ITMAX,LCAVMX,
     &                  ORIGIN,CORPX,CORPY,CORPZ,EXPPI)
C
C     910531-kvm+hjaaj: extracted from MOMINT
C     Revised 920828-hjaaj
C
C     This subroutine calculates Hermitian integrals
C     used in SOLIN1 and MOMINT
C
#include "implicit.h"
#include "pi.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
      DIMENSION HMULT(0:LCAVMX,0:LCAVMX,3), ORIGIN(3)
C
C
C     Overlap integral
C
      SHGTF = SQRT(PI*EXPPI)
C
      EXPPIH = DP5*EXPPI
      CORPXA = CORPX - ORIGIN(1)
      CORPYA = CORPY - ORIGIN(2)
      CORPZA = CORPZ - ORIGIN(3)
C
      HMULT(0,0,1) = SHGTF
      HMULT(0,0,2) = SHGTF
      HMULT(0,0,3) = SHGTF
      DO 100 IO = 1, LCAVMX
         DO 200 IT = 0, MIN(IO,ITMAX+LCAVMX-IO)
            HX = D0
            HY = D0
            HZ = D0
            IF (IT .GT. 0) THEN
               HX = HX + IT*HMULT(IT-1,IO-1,1)
               HY = HY + IT*HMULT(IT-1,IO-1,2)
               HZ = HZ + IT*HMULT(IT-1,IO-1,3)
            END IF
            IF (IT .LE. IO-1) THEN
               HX = HX + CORPXA*HMULT(IT,IO-1,1)
               HY = HY + CORPYA*HMULT(IT,IO-1,2)
               HZ = HZ + CORPZA*HMULT(IT,IO-1,3)
            END IF
            IF (IT .LE. IO-2) THEN
               HX = HX + EXPPIH*HMULT(IT+1,IO-1,1)
               HY = HY + EXPPIH*HMULT(IT+1,IO-1,2)
               HZ = HZ + EXPPIH*HMULT(IT+1,IO-1,3)
            END IF
            HMULT(IT,IO,1) = HX
            HMULT(IT,IO,2) = HY
            HMULT(IT,IO,3) = HZ
  200    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck rlmset */
      SUBROUTINE RLMSET(LCAVMX,LO,MO,NO)
C
C 27-May-1992 kvm+hjaaj
C
#include "implicit.h"
C
      DIMENSION LO(*), MO(*), NO(*)
C
C     Set up list of L,M,N for all orders for Cartesian integrals
C
      DO 100 L = 0, LCAVMX
         LSTART = 1 + L*(L+1)*(L+2)/6
         CALL LMNVAL(L+1,(L+1)*(L+2)/2,LO(LSTART),MO(LSTART),NO(LSTART))
  100 CONTINUE
      RETURN
C
      END
#ifdef UNDEF
/* Comdeck log */
960417- kr
 TLMD is also needed for dipole gradient, thus TLMD for MAXDIF .ge. 2 or DIFDIP
921209- get cavity center from CBISOL; only calc. TLMD for MAXDIF .ge. 2;
        only calc. GSOLT for MAXDIF .ge. 1
921013-kvm+hjaaj
 removed FAC, center C always different from center A
#endif
C  /* Deck avesol */
      SUBROUTINE AVESOL(RLMINT,RLMTAB,FCM,MAXDIF,DIFDIP,ISYMOP,MAXCMP,
     &                  DSHELL,TLMD,LMNO,IPRINT,HSOLT2)
C
C
C     RLMINT(pq,LMN,a) = R(a)(LMN)pq
C     RLMTAB(pq,ab) = t(ab)pq = sum(LMN) fcm(LMN) * R(ab)(LMN)pq
C     where a,b nuclear coordinates; p,q ao indicies, LMN cartesian moment
C
C     Feb. 90, kvm+tuh (then called RLMAVX; based on AVENA2)
C     Expectation terms from solvent part;
C           <u(A)/Rlm(C)/u(B)>
C     Consider Rlm(C) derivatives explicit;
C           SUM(pq) Dpq(tpq)c =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)c):
C     MAY 92
C     Corrections(counters are being stepped)kvm+hjaaj
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.00 D00, D1 = 1.00 D00, D2 = 2.00D00)
      PARAMETER (DP5 = 0.5 D00)
C
C Used from common blocks:
C  ONECOM : ONECEN, ?
C  CBISOL : NCNTCV
C
#include "onecom.h"
#include "symmet.h"
#include "cbisol.h"
#include "taysol.h"
C
      LOGICAL DOAX, DOAY, DOAZ, DOBX, DOBY, DOBZ,
     *        DOCX, DOCY, DOCZ, DIFDIP
      DIMENSION DSHELL(MAXCMP), FCM(LMNTOT), TLMD(LMNTOT,*)
      DIMENSION RLMINT(KCKTAB,LMNTOT,7), RLMTAB(KCKTAB,21),
     &          LMNO(LMNTOT,3), HSOLT2(MXCOOR,MXCOOR)
#include "ibtfun.h"
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
      NCX = 3*NCNTCV - 2
      NCY = 3*NCNTCV - 1
      NCZ = 3*NCNTCV
      FAB = D1
      IF (NCENTA .EQ. NCENTB) FAB = D2
C
      IF (IPRINT .GT. 4) THEN
         CALL HEADER('Output from AVESOL',-1)
         WRITE (LUPRI,*) 'IPRINT',IPRINT
         WRITE (LUPRI,*) 'MAXDIF,ISYMOP',MAXDIF,ISYMOP
         WRITE (LUPRI,*) 'MAXCMP,KCKTAB,LMNTOT',MAXCMP,KCKTAB,LMNTOT
         WRITE (LUPRI,*) 'ONECEN', ONECEN
         IF (MAXDIF .GE. 1 .AND. IPRINT .GT. 6) THEN
            WRITE (LUPRI,*) 'DSHELL(1:MAXCMP)'
            WRITE (LUPRI,'(5F15.10)') (DSHELL(I),I=1,MAXCMP)
            WRITE (LUPRI,*) 'Primitive undiff. RLMINT(KCKTAB,LMNTOT)'
            CALL OUTPUT(RLMINT(1,1,1),1,KCKTAB,1,LMNTOT,
     &                  KCKTAB,LMNTOT,1,LUPRI)
            WRITE (LUPRI,*) 'primitive RLMINT for X0 first derivative'
            CALL OUTPUT(RLMINT(1,1,2),1,KCKTAB,1,LMNTOT,
     &                  KCKTAB,LMNTOT,1,LUPRI)
            WRITE (LUPRI,*) 'primitive RLMINT for Y0 first derivative'
            CALL OUTPUT(RLMINT(1,1,3),1,KCKTAB,1,LMNTOT,
     &                  KCKTAB,LMNTOT,1,LUPRI)
            WRITE (LUPRI,*) 'primitive RLMINT for Z0 first derivative'
            CALL OUTPUT(RLMINT(1,1,4),1,KCKTAB,1,LMNTOT,
     &                  KCKTAB,LMNTOT,1,LUPRI)
            WRITE (LUPRI,*) 'primitive RLMINT for 0X first derivative'
            CALL OUTPUT(RLMINT(1,1,5),1,KCKTAB,1,LMNTOT,
     &                  KCKTAB,LMNTOT,1,LUPRI)
            WRITE (LUPRI,*) 'primitive RLMINT for 0Y first derivative'
            CALL OUTPUT(RLMINT(1,1,6),1,KCKTAB,1,LMNTOT,
     &                  KCKTAB,LMNTOT,1,LUPRI)
            WRITE (LUPRI,*) 'primitive RLMINT for 0Z first derivative'
            CALL OUTPUT(RLMINT(1,1,7),1,KCKTAB,1,LMNTOT,
     &                  KCKTAB,LMNTOT,1,LUPRI)
         END IF
      END IF
C
C        Multiply densities and integrals
C
      DA0000 = D0
      DA0X00 = D0
      DA0Y00 = D0
      DA0Z00 = D0
      DA000X = D0
      DA000Y = D0
      DA000Z = D0
      IF (MAXDIF .GE. 2) THEN
         DAXX00 = D0
         DAXY00 = D0
         DAXZ00 = D0
         DAYY00 = D0
         DAYZ00 = D0
         DAZZ00 = D0
         DA00XX = D0
         DA00XY = D0
         DA00XZ = D0
         DA00YY = D0
         DA00YZ = D0
         DA00ZZ = D0
         DA0X0X = D0
         DA0X0Y = D0
         DA0X0Z = D0
         DA0Y0X = D0
         DA0Y0Y = D0
         DA0Y0Z = D0
         DA0Z0X = D0
         DA0Z0Y = D0
         DA0Z0Z = D0
      END IF
      DO 300 ICOMP = 1, MAXCMP
         DT0000 = D0
         DT0X00 = D0
         DT0Y00 = D0
         DT0Z00 = D0
         DT000X = D0
         DT000Y = D0
         DT000Z = D0
         DO 200 I = 1, LMNTOT
            DT0000 = DT0000 + RLMINT(ICOMP,I,1)*FCM(I)
            DT0X00 = DT0X00 + RLMINT(ICOMP,I,2)*FCM(I)
            DT0Y00 = DT0Y00 + RLMINT(ICOMP,I,3)*FCM(I)
            DT0Z00 = DT0Z00 + RLMINT(ICOMP,I,4)*FCM(I)
            DT000X = DT000X + RLMINT(ICOMP,I,5)*FCM(I)
            DT000Y = DT000Y + RLMINT(ICOMP,I,6)*FCM(I)
            DT000Z = DT000Z + RLMINT(ICOMP,I,7)*FCM(I)
C
  200    CONTINUE
C
         DENS = DSHELL(ICOMP)
C
         DA0000 = DA0000 + DT0000 * DENS
         DA0X00 = DA0X00 + DT0X00 * DENS
         DA0Y00 = DA0Y00 + DT0Y00 * DENS
         DA0Z00 = DA0Z00 + DT0Z00 * DENS
         DA000X = DA000X + DT000X * DENS
         DA000Y = DA000Y + DT000Y * DENS
         DA000Z = DA000Z + DT000Z * DENS
         IF (MAXDIF .GE. 2) THEN
            DAXX00 = DAXX00 + DENS*RLMTAB(ICOMP,1)
            DAXY00 = DAXY00 + DENS*RLMTAB(ICOMP,2)
            DAXZ00 = DAXZ00 + DENS*RLMTAB(ICOMP,3)
            DAYY00 = DAYY00 + DENS*RLMTAB(ICOMP,4)
            DAYZ00 = DAYZ00 + DENS*RLMTAB(ICOMP,5)
            DAZZ00 = DAZZ00 + DENS*RLMTAB(ICOMP,6)
            DA0X0X = DA0X0X + DENS*RLMTAB(ICOMP,7)
            DA0X0Y = DA0X0Y + DENS*RLMTAB(ICOMP,8)
            DA0X0Z = DA0X0Z + DENS*RLMTAB(ICOMP,9)
            DA0Y0X = DA0Y0X + DENS*RLMTAB(ICOMP,10)
            DA0Y0Y = DA0Y0Y + DENS*RLMTAB(ICOMP,11)
            DA0Y0Z = DA0Y0Z + DENS*RLMTAB(ICOMP,12)
            DA0Z0X = DA0Z0X + DENS*RLMTAB(ICOMP,13)
            DA0Z0Y = DA0Z0Y + DENS*RLMTAB(ICOMP,14)
            DA0Z0Z = DA0Z0Z + DENS*RLMTAB(ICOMP,15)
            DA00XX = DA00XX + DENS*RLMTAB(ICOMP,16)
            DA00XY = DA00XY + DENS*RLMTAB(ICOMP,17)
            DA00XZ = DA00XZ + DENS*RLMTAB(ICOMP,18)
            DA00YY = DA00YY + DENS*RLMTAB(ICOMP,19)
            DA00YZ = DA00YZ + DENS*RLMTAB(ICOMP,20)
            DA00ZZ = DA00ZZ + DENS*RLMTAB(ICOMP,21)
         END IF
  300 CONTINUE
C
C        ***** three-center case *****
C
      IAX  = IPTCNT(NAX,0,1)
      IAY  = IPTCNT(NAY,0,1)
      IAZ  = IPTCNT(NAZ,0,1)
      IBX  = IPTCNT(NBX,0,1)
      IBY  = IPTCNT(NBY,0,1)
      IBZ  = IPTCNT(NBZ,0,1)
      ICX  = IPTCNT(NCX,0,1)
      ICY  = IPTCNT(NCY,0,1)
      ICZ  = IPTCNT(NCZ,0,1)
      DOAX = IAX .NE. 0
      DOAY = IAY .NE. 0
      DOAZ = IAZ .NE. 0
      DOBX = IBX .NE. 0
      DOBY = IBY .NE. 0
      DOBZ = IBZ .NE. 0
      DOCX = ICX .NE. 0
      DOCY = ICY .NE. 0
      DOCZ = ICZ .NE. 0
C
C           Undifferentiated expectation value
C           of SUM(pq) Dpq*((-2 Sum(lm) gl Tlm) Rlm) * (0.5)
C                      DENS*sum(lm) (FCM(lm) * R(lm)) * DP5
C
      IF (IPRINT .GT. 4) THEN
         WRITE (LUPRI,*) 'old ESOLTT  ',ESOLTT
         WRITE (LUPRI,*) 'contribution',DP5*DA0000
         WRITE (LUPRI,*) 'new ESOLTT  ',ESOLTT + DP5*DA0000
      END IF
      ESOLTT = ESOLTT + DA0000 * DP5
C
      IF (MAXDIF .GE. 1) THEN
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider u(A) derivatives;
C           SUM(pq) Dpq(tpq)a =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)a):
C
         IF (DOAX) GSOLTT(IAX) = GSOLTT(IAX) + DA0X00
         IF (DOAY) GSOLTT(IAY) = GSOLTT(IAY) + DA0Y00
         IF (DOAZ) GSOLTT(IAZ) = GSOLTT(IAZ) + DA0Z00
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider u(B) derivatives;
C           SUM(pq) Dpq(tpq)b =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)b):
C
         IF (DOBX) GSOLTT(IBX) = GSOLTT(IBX)
     *                      - SIGNBX*(DA0X00 + DA000X)
         IF (DOBY) GSOLTT(IBY) = GSOLTT(IBY)
     *                      - SIGNBY*(DA0Y00 + DA000Y)
         IF (DOBZ) GSOLTT(IBZ) = GSOLTT(IBZ)
     *                      - SIGNBZ*(DA0Z00 + DA000Z)
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider Rlm(C) derivatives;
C           SUM(pq) Dpq(tpq)c =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)c):
C
         IF (DOCX) GSOLTT(ICX) = GSOLTT(ICX) + DA000X
         IF (DOCY) GSOLTT(ICY) = GSOLTT(ICY) + DA000Y
         IF (DOCZ) GSOLTT(ICZ) = GSOLTT(ICZ) + DA000Z
      END IF
C
      IF (MAXDIF .GE. 2 .OR. DIFDIP) THEN
C
         IF (MAXDIF .GE. 2) THEN
            DO 600 IREP = 0, MAXREP
               CHIB = PT(IBTAND(ISYMOP,IREP))
               CSBX = CHIB*SIGNBX
               CSBY = CHIB*SIGNBY
               CSBZ = CHIB*SIGNBZ
               IAX  = IPTCNT(NAX,IREP,1)
               IAY  = IPTCNT(NAY,IREP,1)
               IAZ  = IPTCNT(NAZ,IREP,1)
               IBX  = IPTCNT(NBX,IREP,1)
               IBY  = IPTCNT(NBY,IREP,1)
               IBZ  = IPTCNT(NBZ,IREP,1)
               ICX  = IPTCNT(NCX,IREP,1)
               ICY  = IPTCNT(NCY,IREP,1)
               ICZ  = IPTCNT(NCZ,IREP,1)
               DOAX = IAX .NE. 0
               DOAY = IAY .NE. 0
               DOAZ = IAZ .NE. 0
               DOBX = IBX .NE. 0
               DOBY = IBY .NE. 0
               DOBZ = IBZ .NE. 0
               DOCX = ICX .NE. 0
               DOCY = ICY .NE. 0
               DOCZ = ICZ .NE. 0
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider u(A), u(A) second derivatives;
C           SUM(pq) Dpq(tpq)aa =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)aa):
C
               IF (DOAX)
     *              HSOLT2(IAX,IAX) = HSOLT2(IAX,IAX) + DAXX00
               IF (DOAX.AND.DOAY)
     *              HSOLT2(IAX,IAY) = HSOLT2(IAX,IAY) + DAXY00
               IF (DOAX.AND.DOAZ)
     *              HSOLT2(IAX,IAZ) = HSOLT2(IAX,IAZ) + DAXZ00
               IF (DOAY)
     *              HSOLT2(IAY,IAY) = HSOLT2(IAY,IAY) + DAYY00
               IF (DOAY.AND.DOAZ)
     *              HSOLT2(IAY,IAZ) = HSOLT2(IAY,IAZ) + DAYZ00
               IF (DOAZ)
     *              HSOLT2(IAZ,IAZ) = HSOLT2(IAZ,IAZ) + DAZZ00
C     
C           Consider u(A), u(B) second derivatives;
C           SUM(pq) Dpq(tpq)ab =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)ab):
C
               IF (DOAX.AND.DOBX)
     *              HSOLT2(IAX,IBX) = HSOLT2(IAX,IBX)
     *                              - FAB*CSBX*(DAXX00 + DA0X0X)
               IF (DOAX.AND.DOBY)
     *              HSOLT2(IAX,IBY) = HSOLT2(IAX,IBY)
     *                              - CSBY*(DAXY00 + DA0X0Y)
               IF (DOAX.AND.DOBZ)
     *              HSOLT2(IAX,IBZ) = HSOLT2(IAX,IBZ)
     *                              - CSBZ*(DAXZ00 + DA0X0Z)
               IF (DOAY.AND.DOBX)
     *              HSOLT2(IAY,IBX) = HSOLT2(IAY,IBX)
     *                              - CSBX*(DAXY00 + DA0Y0X)
               IF (DOAY.AND.DOBY)
     *              HSOLT2(IAY,IBY) = HSOLT2(IAY,IBY)
     *                              - FAB*CSBY*(DAYY00 + DA0Y0Y)
               IF (DOAY.AND.DOBZ)
     *              HSOLT2(IAY,IBZ) = HSOLT2(IAY,IBZ)
     *                              - CSBZ*(DAYZ00 + DA0Y0Z)
               IF (DOAZ.AND.DOBX)
     *              HSOLT2(IAZ,IBX) = HSOLT2(IAZ,IBX)
     *                              - CSBX*(DAXZ00 + DA0Z0X)
               IF (DOAZ.AND.DOBY)
     *              HSOLT2(IAZ,IBY) = HSOLT2(IAZ,IBY)
     *                             - CSBY*(DAYZ00 + DA0Z0Y)
               IF (DOAZ.AND.DOBZ)
     *              HSOLT2(IAZ,IBZ) = HSOLT2(IAZ,IBZ)
     *                              - FAB*CSBZ*(DAZZ00 + DA0Z0Z)
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider u(A), Rlm(C) second derivatives;
C           SUM(pq) Dpq(tpq)ac =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)ac):
C
               IF (DOAX.AND.DOCX)
     *              HSOLT2(IAX,ICX) = HSOLT2(IAX,ICX) + DA0X0X
               IF (DOAX.AND.DOCY)
     *              HSOLT2(IAX,ICY) = HSOLT2(IAX,ICY) + DA0X0Y
               IF (DOAX.AND.DOCZ)
     *              HSOLT2(IAX,ICZ) = HSOLT2(IAX,ICZ) + DA0X0Z
               IF (DOAY.AND.DOCX)
     *              HSOLT2(IAY,ICX) = HSOLT2(IAY,ICX) + DA0Y0X
               IF (DOAY.AND.DOCY)
     *              HSOLT2(IAY,ICY) = HSOLT2(IAY,ICY) + DA0Y0Y
               IF (DOAY.AND.DOCZ)
     *              HSOLT2(IAY,ICZ) = HSOLT2(IAY,ICZ) + DA0Y0Z
               IF (DOAZ.AND.DOCX)
     *              HSOLT2(IAZ,ICX) = HSOLT2(IAZ,ICX) + DA0Z0X
               IF (DOAZ.AND.DOCY)
     *              HSOLT2(IAZ,ICY) = HSOLT2(IAZ,ICY) + DA0Z0Y
               IF (DOAZ.AND.DOCZ)
     *              HSOLT2(IAZ,ICZ) = HSOLT2(IAZ,ICZ) + DA0Z0Z
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider u(B), u(B) second derivatives;
C           SUM(pq) Dpq(tpq)bb =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)bb)
C
               IF (DOBX)
     *              HSOLT2(IBX,IBX) = HSOLT2(IBX,IBX) +
     *              (DAXX00 + DA00XX + DA0X0X + DA0X0X)
               IF (DOBX.AND.DOBY)
     *              HSOLT2(IBX,IBY) = HSOLT2(IBX,IBY) + SIGNBX*SIGNBY*
     *              (DAXY00 + DA00XY + DA0X0Y + DA0Y0X)
               IF (DOBX.AND.DOBZ)
     *              HSOLT2(IBX,IBZ) = HSOLT2(IBX,IBZ) + SIGNBX*SIGNBZ*
     *              (DAXZ00 + DA00XZ + DA0X0Z + DA0Z0X)
               IF (DOBY)
     *              HSOLT2(IBY,IBY) = HSOLT2(IBY,IBY) +
     *              (DAYY00 + DA00YY + DA0Y0Y + DA0Y0Y)
               IF (DOBY.AND.DOBZ)
     *              HSOLT2(IBY,IBZ) = HSOLT2(IBY,IBZ) + SIGNBY*SIGNBZ*
     *              (DAYZ00 + DA00YZ + DA0Y0Z + DA0Z0Y)
               IF (DOBZ)
     *              HSOLT2(IBZ,IBZ) = HSOLT2(IBZ,IBZ) +
     *              (DAZZ00 + DA00ZZ + DA0Z0Z + DA0Z0Z)
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider Rlm(C), u(B) second derivatives;
C           SUM(pq) Dpq(tpq)bc =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)bc):
C
               IF (DOBX.AND.DOCX) HSOLT2(IBX,ICX) = HSOLT2(IBX,ICX)
     *              - CSBX*(DA0X0X + DA00XX)
               IF (DOBX.AND.DOCY) HSOLT2(IBX,ICY) = HSOLT2(IBX,ICY)
     *              - CSBX*(DA0X0Y + DA00XY)
               IF (DOBX.AND.DOCZ) HSOLT2(IBX,ICZ) = HSOLT2(IBX,ICZ)
     *              - CSBX*(DA0X0Z + DA00XZ)
               IF (DOBY.AND.DOCX) HSOLT2(IBY,ICX) = HSOLT2(IBY,ICX)
     *              - CSBY*(DA0Y0X + DA00XY)
               IF (DOBY.AND.DOCY) HSOLT2(IBY,ICY) = HSOLT2(IBY,ICY)
     *              - CSBY*(DA0Y0Y + DA00YY)
               IF (DOBY.AND.DOCZ) HSOLT2(IBY,ICZ) = HSOLT2(IBY,ICZ)
     *              - CSBY*(DA0Y0Z + DA00YZ)
               IF (DOBZ.AND.DOCX) HSOLT2(IBZ,ICX) = HSOLT2(IBZ,ICX)
     *              - CSBZ*(DA0Z0X + DA00XZ)
               IF (DOBZ.AND.DOCY) HSOLT2(IBZ,ICY) = HSOLT2(IBZ,ICY)
     *              - CSBZ*(DA0Z0Y + DA00YZ)
               IF (DOBZ.AND.DOCZ) HSOLT2(IBZ,ICZ) = HSOLT2(IBZ,ICZ)
     *              - CSBZ*(DA0Z0Z + DA00ZZ)
C
C           Consider RLM(C), RLM(C) second derivatives;
C           SUM(pq) Dpq(tpq)cc =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)cc):
C
               IF (DOCX)
     *              HSOLT2(ICX,ICX) = HSOLT2(ICX,ICX) + DA00XX
               IF (DOCX.AND.DOCY)
     *              HSOLT2(ICX,ICY) = HSOLT2(ICX,ICY) + DA00XY
               IF (DOCX.AND.DOCZ)
     *              HSOLT2(ICX,ICZ) = HSOLT2(ICX,ICZ) + DA00XZ
               IF (DOCY)
     *              HSOLT2(ICY,ICY) = HSOLT2(ICY,ICY) + DA00YY
               IF (DOCY.AND.DOCZ)
     *              HSOLT2(ICY,ICZ) = HSOLT2(ICY,ICZ) + DA00YZ
               IF (DOCZ)
     *              HSOLT2(ICZ,ICZ) = HSOLT2(ICZ,ICZ) + DA00ZZ
 600        CONTINUE
         END IF
C
C
         CALL RLMSET(LCAVMX,LMNO(1,1),LMNO(1,2),LMNO(1,3))
         DO 700 I = 1, LMNTOT
            ISX = MOD(LMNO(I,1),2)*ISYMAX(1,1)
            ISY = MOD(LMNO(I,2),2)*ISYMAX(2,1)
            ISZ = MOD(LMNO(I,3),2)*ISYMAX(3,1)
            IREP = IBTXOR(ISX,IBTXOR(ISY,ISZ))
            CHIB = PT(IBTAND(ISYMOP,IREP))
            CSBX = CHIB*SIGNBX
            CSBY = CHIB*SIGNBY
            CSBZ = CHIB*SIGNBZ
            IAX  = IPTCNT(NAX,IREP,1)
            IAY  = IPTCNT(NAY,IREP,1)
            IAZ  = IPTCNT(NAZ,IREP,1)
            IBX  = IPTCNT(NBX,IREP,1)
            IBY  = IPTCNT(NBY,IREP,1)
            IBZ  = IPTCNT(NBZ,IREP,1)
            ICX  = IPTCNT(NCX,IREP,1)
            ICY  = IPTCNT(NCY,IREP,1)
            ICZ  = IPTCNT(NCZ,IREP,1)
            TA0X00 = D0
            TA0Y00 = D0
            TA0Z00 = D0
            TA000X = D0
            TA000Y = D0
            TA000Z = D0
            DO 730 ICOMP = 1,MAXCMP
               DENS = DSHELL(ICOMP)
               TA0X00 = TA0X00 + RLMINT(ICOMP,I,2)*DENS
               TA0Y00 = TA0Y00 + RLMINT(ICOMP,I,3)*DENS
               TA0Z00 = TA0Z00 + RLMINT(ICOMP,I,4)*DENS
               TA000X = TA000X + RLMINT(ICOMP,I,5)*DENS
               TA000Y = TA000Y + RLMINT(ICOMP,I,6)*DENS
               TA000Z = TA000Z + RLMINT(ICOMP,I,7)*DENS
  730       CONTINUE
            TBX = - TA0X00 - TA000X
            TBY = - TA0Y00 - TA000Y
            TBZ = - TA0Z00 - TA000Z
            IF (IAX .NE. 0) TLMD(I,IAX) = TLMD(I,IAX) + TA0X00
            IF (IAY .NE. 0) TLMD(I,IAY) = TLMD(I,IAY) + TA0Y00
            IF (IAZ .NE. 0) TLMD(I,IAZ) = TLMD(I,IAZ) + TA0Z00
            IF (IBX .NE. 0) TLMD(I,IBX) = TLMD(I,IBX) + CSBX*TBX
            IF (IBY .NE. 0) TLMD(I,IBY) = TLMD(I,IBY) + CSBY*TBY
            IF (IBZ .NE. 0) TLMD(I,IBZ) = TLMD(I,IBZ) + CSBZ*TBZ
            IF (ICX .NE. 0) TLMD(I,ICX) = TLMD(I,ICX) + TA000X
            IF (ICY .NE. 0) TLMD(I,ICY) = TLMD(I,ICY) + TA000Y
            IF (ICZ .NE. 0) TLMD(I,ICZ) = TLMD(I,ICZ) + TA000Z
  700    CONTINUE
      END IF
      RETURN
      END
C  /* Deck fcmfac */
      SUBROUTINE FCMFAC(LCAVMX,LMNTOT,LMTOT,GLM,TLM,FCM,
     *                  WORK,LWORK,IPRINT)
C
C 13-Oct-1992 hjaaj
C
C FCM(LMN) = -2 * g(l) * sum(m) [T(lm) * C(LMN,m)]
C  where l = L+M+N
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION GLM(LMTOT), TLM(LMTOT), FCM(LMNTOT), WORK(LWORK)
C
C
      MXLM  = 2*LCAVMX + 1
      MXXYZ = (LCAVMX+1)*(LCAVMX+2)/2
      KTRAMA = 1
      KWRK1  = KTRAMA + MXLM*MXXYZ
      IF (KWRK1 .GT. LWORK) CALL STOPIT('FCMFAC',' ',LWORK,KWRK1)
      LWRK1  = LWORK + 1 - KWRK1
      CALL FCMFA1(LCAVMX,LMNTOT,LMTOT,GLM,TLM,FCM,
     *            WORK(KTRAMA),MXXYZ,MXLM,WORK(KWRK1),LWRK1,IPRINT)
      IF (IPRINT .GE. 3) THEN
         WRITE (LUPRI,*) 'FCMFAC: Solvent test information'
         WRITE (LUPRI,*) 'LCAVMX',LCAVMX
         WRITE (LUPRI,*) 'LMNTOT',LMNTOT
         WRITE (LUPRI,*) 'LMTOT ',LMTOT
         WRITE (LUPRI,*) 'IPRINT',IPRINT
         WRITE (LUPRI,*) 'LWORK ',LWORK
         WRITE (LUPRI,'(/A/A)')
     &   '       LM           GLM(LM)             TLM(LM)',
     &   '     =====     ===============     ==============='
         WRITE (LUPRI,'(I10,2F20.12)') (I,GLM(I),TLM(I),I=1,LMTOT)
         WRITE (LUPRI,'(/A/A)')
     &   '      LMN          FCM(LMN)',
     &   '     =====     ==============='
         WRITE (LUPRI,'(I10,F20.12)') (I,FCM(I),I=1,LMNTOT)
      END IF
      RETURN
      END
C  /* Deck fcmfa1 */
      SUBROUTINE FCMFA1(LCAVMX,LMNTOT,LMTOT,GLM,TLM,FCM,
     *                  TRAMAT,MXXYZ,MXLM,WORK,LWORK,IPRINT)
C
C 13-Oct-1992 hjaaj
C
#include "implicit.h"
      PARAMETER (DM2 = -2.0D0, D0 = 0.0D0)
      PARAMETER (MORDER = 1, MINTEG = 0)
C
      DIMENSION GLM(LMTOT), TLM(LMTOT), FCM(LMNTOT)
      DIMENSION TRAMAT(MXXYZ,MXLM), WORK(LWORK)
C
      JLM  = 0
      JXYZ = 0
      DO 800 IO = 0,LCAVMX
         NLM  = 2*IO + 1
         NXYZ = (IO+1)*(IO+2)/2
         CALL SPHCOM(IO,TRAMAT,MXLM,MXXYZ,MORDER,MINTEG,WORK,LWORK,
     &               IPRINT)
         DO 400 I = 1,NXYZ
            FCMI = D0
            DO 200 J = 1,NLM
               FCMI = FCMI + GLM(JLM+J) * TLM(JLM+J) * TRAMAT(I,J)
  200       CONTINUE
            FCM(JXYZ+I) = DM2 * FCMI
  400    CONTINUE
         JLM  = JLM  + NLM
         JXYZ = JXYZ + NXYZ
  800 CONTINUE
      IF (JLM .NE. LMTOT) THEN
         CALL QUIT('Final JLM .ne. LMTOT in FCMFAC')
      END IF
      IF (JXYZ.NE. LMNTOT) THEN
         CALL QUIT('Final JXYZ .ne. LMNTOT in FCMFAC')
      END IF
C
      RETURN
      END
#ifndef PRG_DIRAC
C  /* Deck magsol */
      SUBROUTINE MAGSOL(DERMAT,NCOMP,IPRINT,INTPRI,WORK,LWORK)
#include "implicit.h"
#include "iratdef.h"
#include "priunit.h"
#include "mxcent.h"
C
      DIMENSION DERMAT(*), WORK(LWORK)
#include "cbisol.h"
#include "inforb.h"
C
C
      CALL QENTER('MAGSOL')
      KGLM   = 1
      KTLM   = KGLM + LMTOT
      KAOMAT = KTLM + LMTOT
      KTEMPM = KAOMAT + (LCAVMX + 1)*(LCAVMX + 2)*NCOMP*NNBASX/2
      KTRAMA = KTEMPM + NCOMP*NNBASX*(2*LCAVMX + 1)
      KLAST  = KTRAMA + (LCAVMX + 1)*(LCAVMX + 2)*(2*LCAVMX + 1)/2
      LWRK   = LWORK - KLAST + 1
      IF (2*KLAST .GT. LWORK) THEN
         KTRAMA = KTLM + LMTOT
         KIX    = KTRAMA + (LCAVMX + 1)*(LCAVMX + 2)*(2*LCAVMX + 1)/2
         KIY    = KIX + (9*MXCENT + 1)/IRAT
         KIZ    = KIY + (9*MXCENT + 1)/IRAT
         KINTRP = KIZ + (9*MXCENT + 1)/IRAT
         KINTAD = KINTRP + (9*MXCENT + 1)/IRAT
         KLAST  = KINTAD + (9*MXCENT + 1)/IRAT
         LWRK   = LWORK - KLAST + 1
         WRITE (LUPRI,'(A)') ' Because of memory restrictions, diff.'//
     &        'solvent integrals written to disk'
         IF ((KLAST + (2*LCAVMX + 1)*NNBASX) .GT. LWORK)
     &        CALL STOPIT('MAGSOL','TOFILE',KLAST,LWORK)
         CALL MAGSO2(DERMAT,WORK(KGLM),WORK(KTLM),WORK(KIX),
     &               WORK(KIY),WORK(KIZ),WORK(KINTRP),WORK(KINTAD),
     &               WORK(KTRAMA),NCOMP,IPRINT,INTPRI,WORK(KLAST),LWRK)
      ELSE
         CALL MAGSO1(DERMAT,WORK(KGLM),WORK(KTLM),
     &               WORK(KAOMAT),WORK(KTRAMA),
     &               WORK(KTEMPM),NCOMP,IPRINT,INTPRI,WORK(KLAST),LWRK)
      END IF
      CALL QEXIT('MAGSOL')
      RETURN
      END
C  /* Deck magso1 */
      SUBROUTINE MAGSO1(DERMAT,GLM,TLM,AOMAT1,
     &                  TRAMAT,TEMPMA,NCOMP,IPRINT,INTPRI,WORK,LWORK)
C
C     kr and kvm, Dec.-94
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "iratdef.h"
      PARAMETER (MORDER = 1, MINTEG=0, D2 = 2.0D0)
      DIMENSION GLM(LMTOT), TLM(LMTOT),
     &          AOMAT1(NNBASX,NCOMP,(LCAVMX + 1)*(LCAVMX + 2)/2),
     &          DERMAT(*), TEMPMA(NNBASX,NCOMP,-LCAVMX:LCAVMX),
     &          TRAMAT((LCAVMX + 1)*(LCAVMX + 2)*(2*LCAVMX + 1)/2),
     &          WORK(LWORK)
#include <symmet.h>
#include <cbisol.h>
#include <orgcom.h>
#include <mgsolt.h>
#include <inforb.h>
#include <infinp.h>
#include <inftap.h>
C
      DIMENSION CAVCNT(3)
      LOGICAL   FOUND
C
C     Construct GLM, TLM, and density matrix
C
      CALL SOLFL(GLM,EPDIEL,RCAV,LCAVMX)
      CALL RD_SIRIFC('SOLTLM',FOUND,TLM,WORK,LWORK)
      IF (.NOT. FOUND)
     &    CALL QUIT('MAGSO1 ERROR: "SOLTLM" not found on SIRIFC.')
C
      CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUPRI)
      READ (LUSOL) LMAXSS, LMTOTX, NNNBAS, CAVCNT
      IF (LMAXSS .LT. LSOLMX) THEN
         WRITE (LUPRI,'(//2A,2(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' insufficient number of intgrals on LUSOL',
     *   ' l max from SIRIUS input :',LSOLMX,
     *   ' l max from LUSOL  file  :',LMAXSS
         CALL QUIT('MAGSO1: lmax on LUSOL is too small')
      END IF
      IF ((LMAXSS+1)**2 .NE. LMTOTX) THEN
         WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' LUSOL file info inconsistent',
     *   ' l_max               :',LMAXSS,
     *   ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *   ' LMTOTX              :',LMTOTX
         CALL QUIT('MAGSO1: LUSOL info not internally consistent')
      END IF
      IF (NNNBAS .NE. NBAST) THEN
         WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' LUSOL file info inconsistent with SIRIUS input',
     *   ' NBAST - LUSOL       :',NNNBAS,
     *   ' NBAST - SIRIUS      :',NBAST
         CALL QUIT('MAGSO1: LUSOL info inconsistent with SIRIUS input.')
      END IF
c      IF (CAVCNT(1) .NE. CAVORG(1) .OR.
c     &    CAVCNT(2) .NE. CAVORG(2) .OR.
c     &    CAVCNT(3) .NE. CAVORG(3) ) THEN
c         WRITE (LUPRI,'(//2A,3(/A,3F20.15))') ' >>> SOLRH1 ERROR,',
c     &   ' LUSOL center of cavity not consistent with ABACUS value',
c     &   ' CAVORG(1:3) from LUSOL        :',CAVCNT,
c     &   ' CAVORG(1:3) from common block :',CAVORG
c         CALL QUIT('MAGSO1: LUSOL cavity center .ne. ABACUS value.')
c      END IF
C
      READ (LUSOL)
      LM = 0
      ICOUNT = 0
      DOALL_MGSOLT  = .TRUE.
      CALL DZERO(DERMAT,NCOMP*N2BASX)
C
      DO 100 L = 0, LSOLMX
         MCOMP = (L + 1)*(L + 2)*NCOMP/2
C
C     Construct differentiated integrals in spherical harmonic basis
C
         CALL GET1PR(AOMAT1,'LONSOL1',MCOMP,'TRIANG',.TRUE.,WORK,
     &               LWORK,L,INTPRI)
         CALL SPHCOM(L,TRAMAT,2*L + 1,(L + 1)*(L + 2)/2,MORDER,MINTEG,
     &               WORK,LWORK,IPRINT)
         CALL DZERO(TEMPMA(1,1,-LCAVMX),NCOMP*(2*LCAVMX + 1)*NNBASX)
         IF (IPRINT .GE. 10) THEN
            CALL AROUND('SOLINT  integrals in cartesian basis')
            DO 103 ICOMP = 1,(L + 1)*(L + 2)/2
               DO 104 ICOOR = 1, NCOMP
                  WRITE (LUPRI,'(2X,A,I3,5X,A,I3)') 'Component: ',
     &                  ICOMP, 'Coordinate direction: ',ICOOR
                  CALL OUTPAK(AOMAT1(1,ICOOR,ICOMP),NBAST,1,LUPRI)
 104           CONTINUE
 103        CONTINUE
         END IF
         DO 430 ICOOR = 1, NCOMP
            DO 400 ICAR = 1, (L + 1)*(L + 2)/2
               IL = 0
               DO 410 ISPH = -L, L
                  IL = IL + 1
                  IELMNT = (L + 1)*(L + 2)*(IL - 1)/2 + ICAR
                  COEF = TRAMAT(IELMNT)
                  DO 420 IBAS1 = 1, NNBASX
                     TEMPMA(IBAS1,ICOOR,ISPH) = TEMPMA(IBAS1,ICOOR,ISPH)
     &                     + AOMAT1(IBAS1,ICOOR,ICAR)*COEF
 420              CONTINUE
 410           CONTINUE
 400        CONTINUE
 430     CONTINUE
         IF (IPRINT .GE. 10) THEN
            CALL AROUND('SOLINT1 integrals in spherical harmonic basis')
            DO 101 ICOMP = -L, L
               DO 102 ICOOR = 1, NCOMP
                  WRITE (LUPRI,'(2X,A,I3,5X,A,I3)') 'Component: ',
     &                  ICOMP, 'Coordinate direction: ',ICOOR
                  CALL OUTPAK(TEMPMA(1,ICOOR,ICOMP),NBAST,1,LUPRI)
 102           CONTINUE
 101        CONTINUE
         END IF
C
C     Construct and -2g(l)sum(m) T(lm) t(lm)a_mn
C
         DO 200 M = -L, L
            ICOUNT = ICOUNT + 1
            DO 190 ICOOR = 1, NCOMP
               FAC2 = -D2*GLM(ICOUNT)*TLM(ICOUNT)
               IOFF = 1 + (ICOOR - 1)*NNBASX
               CALL DAXPY(NNBASX,FAC2,TEMPMA(1,ICOOR,M),1,
     &                    DERMAT(IOFF),1)
 190        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
C        We need to square things away, in particular certain matrices
C
      CALL DCOPY(NCOMP*NNBASX,DERMAT,1,TEMPMA(1,1,-LCAVMX),1)
      DO 300 ICOMP = 1, NCOMP
         IOFFS = 1 + (ICOMP - 1)*N2BASX
         CALL DAPTGE(NBAST,TEMPMA(1,ICOMP,-LCAVMX),DERMAT(IOFFS))
 300  CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL AROUND('Solvent contributions to gradient in MAGSO1')
         IOFF = 1
         WRITE (LUPRI,'(2X,A)') 'X coordinate'
         CALL OUTPUT(DERMAT(IOFF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IOFF = IOFF + N2BASX
         WRITE (LUPRI,'(2X,A)') 'Y coordinate'
         CALL OUTPUT(DERMAT(IOFF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IOFF = IOFF + N2BASX
         WRITE (LUPRI,'(2X,A)') 'Z coordinate'
         CALL OUTPUT(DERMAT(IOFF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     We ought to be more or less finished now
C
      CALL GPCLOSE(LUSOL,'KEEP')
      RETURN
      END
C  /* Deck magso2 */
      SUBROUTINE MAGSO2(DERMAT,GLM,TLM,IX,IY,IZ,INTREP,INTADR,
     &                  TRAMAT,NCOMP,IPRINT,INTPRI,WORK,LWORK)
C
C     kr , Nov.-96
C     Writes integrals to file in order to save memory
C     Rewritten to save even more memory, May 4-97
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "iratdef.h"
#include "priunit.h"
      PARAMETER (MORDER = 1, MINTEG=0, D2 = 2.0D0)
      CHARACTER*2 LBL3(3), LBL6(6), LABXY
      CHARACTER*8 LABINT(9*MXCENT), LABEL, RTNLBL(2)
      DIMENSION GLM(LMTOT), TLM(LMTOT), INTREP(*), INTADR(*),
     &          DERMAT(*), 
     &          TRAMAT((LCAVMX + 1)*(LCAVMX + 2)*(2*LCAVMX + 1)/2),
     &          IX(9*MXCENT), IY(9*MXCENT), IZ(9*MXCENT),
     &          WORK(LWORK)
      DIMENSION CAVCNT(3)
      LOGICAL   FOUND,FNDLAB
C
#include <symmet.h>
#include <cbisol.h>
#include <mgsolt.h>
#include <orgcom.h>
#include <inforb.h>
#include <infinp.h>
#include <inftap.h>
#include <chrnos.h>
C
      DATA (LBL3(I), I = 1, 3) /'X ','Y ', 'Z '/
      DATA (LBL6(I), I = 1, 6) /'XX','XY','XZ','YY','YZ','ZZ'/
C
C     Construct GLM, TLM, and density matrix
C
      CALL SOLFL(GLM,EPDIEL,RCAV,LCAVMX)
      CALL RD_SIRIFC('SOLTLM',FOUND,TLM,WORK,LWORK)
      IF (.NOT. FOUND)
     &    CALL QUIT('MAGSO2 ERROR: "SOLTLM" not found on SIRIFC.')
C
      CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUPRI)
      READ (LUSOL) LMAXSS, LMTOTX, NNNBAS, CAVCNT
      IF (LMAXSS .LT. LSOLMX) THEN
         WRITE (LUPRI,'(//2A,2(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' insufficient number of intgrals on LUSOL',
     *   ' l max from SIRIUS input :',LSOLMX,
     *   ' l max from LUSOL  file  :',LMAXSS
         CALL QUIT('MAGSO1: lmax on LUSOL is too small')
      END IF
      IF ((LMAXSS+1)**2 .NE. LMTOTX) THEN
         WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' LUSOL file info inconsistent',
     *   ' l_max               :',LMAXSS,
     *   ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *   ' LMTOTX              :',LMTOTX
         CALL QUIT('MAGSO1: LUSOL info not internally consistent')
      END IF
      IF (NNNBAS .NE. NBAST) THEN
         WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' LUSOL file info inconsistent with SIRIUS input',
     *   ' NBAST - LUSOL       :',NNNBAS,
     *   ' NBAST - SIRIUS      :',NBAST
         CALL QUIT('MAGSO1: LUSOL info inconsistent with SIRIUS input.')
      END IF
c      IF (CAVCNT(1) .NE. CAVORG(1) .OR.
c     &    CAVCNT(2) .NE. CAVORG(2) .OR.
c     &    CAVCNT(3) .NE. CAVORG(3) ) THEN
c         WRITE (LUPRI,'(//2A,3(/A,3F20.15))') ' >>> SOLRH1 ERROR,',
c     &   ' LUSOL center of cavity not consistent with ABACUS value',
c     &   ' CAVORG(1:3) from LUSOL        :',CAVCNT,
c     &   ' CAVORG(1:3) from common block :',CAVORG
c         CALL QUIT('MAGSO1: LUSOL cavity center .ne. ABACUS value.')
c      END IF
C
      READ (LUSOL)
      LM = 0
      DOALL_MGSOLT = .FALSE.
      CALL DZERO(DERMAT,NCOMP*N2BASX)
C
C     Calculate all integrals, order by order, component by component
C
      DO 100 L = 0, LSOLMX
         DO 110 ICOOR = 1, NCOMP
            ICOMP = ICOOR
            MCOMP = (L + 1)*(L + 2)/2
C
C     Construct differentiated integrals in spherical harmonic basis
C
            KPATOM = 0
            CALL GET1IN(DUMMY,'LONSOL1',0,WORK,LWORK,LABINT,INTREP,
     &                  INTADR,L,.TRUE.,KPATOM,.TRUE.,DUMMY,.FALSE.,
     &                  DUMMY,INTPRI)
 110     CONTINUE 
 100  CONTINUE 
C
C     All integrals have been calculated. Now we read them from disk,
C     transforms to spherical basis, and adds them to the appropriate
C     differentiated matrices
C
      CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      KTMP   = 1
      KSTART = KTMP + (2*LCAVMX + 1)*NNBASX
      KEND   = KSTART + (LCAVMX + 1)*(LCAVMX + 2)*NNBASX/2
      IF (KEND .GT. LWORK) CALL STOPIT('MAGSO2',' ',KEND,LWORK)
      DO 140 ICOOR = 1, NCOMP
         LABXY = LBL3(ICOOR)
         ICOUNT = 0
         DO 120 L = 0, LSOLMX
            CALL DZERO(WORK(KTMP),(2*LCAVMX + 1)*NNBASX)
            CALL SPHCOM(L,TRAMAT,2*L + 1,(L + 1)*(L + 2)/2,MORDER,
     &                  MINTEG,WORK(KEND),LWORK,IPRINT)
            CALL LMNVAL(L + 1,(L + 1)*(L + 2)/2,IX,IY,IZ)
            DO 130 LCAR = 1, (L + 1)*(L + 2)/2
               NX = IX(LCAR)
               NY = IY(LCAR)
               NZ = IZ(LCAR)
               LABEL = LABXY//CHRNOS(NX/10)//CHRNOS(MOD(NX,10))
     &                      //CHRNOS(NY/10)//CHRNOS(MOD(NY,10))
     &                      //CHRNOS(NZ/10)//CHRNOS(MOD(NZ,10))
               IADR = (LCAR - 1)*NNBASX + KSTART
               REWIND (LUPROP)
               CALL MOLLAB(LABEL,LUPROP,LUPRI)
               CALL READT(LUPROP,NNBASX,WORK(IADR))
               IF (IPRINT .GE. 10) THEN
                  CALL AROUND('SOLINT  integrals in cartesian basis')
                  WRITE (LUPRI,'(2X,A,I3,5X,A,I3)') 'Component: ',
     &                 ICOMP, 'Coordinate direction: ',ICOOR
                  CALL OUTPAK(WORK(IADR),NBAST,1,LUPRI)
               END IF
 130        CONTINUE 
C
C     We now have all Cartesian components of a given direction of the
C     property solvent integrals. Need to transform them to spherical basis
C
            DO 150 ICAR = 1, (L + 1)*(L + 2)/2
               IADR = (ICAR - 1)*NNBASX + KSTART
               IL = 0
               DO 410 ISPH = -L, L
                  IL = IL + 1
                  IELMNT = (L + 1)*(L + 2)*(IL - 1)/2 + ICAR
                  COEF = TRAMAT(IELMNT)
                  IPOINT = KTMP + (IL - 1)*NNBASX
                  DO 420 IBAS1 = 0, NNBASX - 1
                     WORK(IPOINT + IBAS1) = WORK(IPOINT + IBAS1)
     &                                    + WORK(IBAS1 + IADR)*COEF
 420              CONTINUE
 410           CONTINUE
 150        CONTINUE
            IF (IPRINT .GE. 10) THEN
               CALL AROUND('SOLINT1 integrals in spherical '//
     &              'harmonic basis')
               IL = 0
               DO 101 ICOMP = -L, L
                  IL = IL + 1
                  IPOINT = KTMP + (IL - 1)*NNBASX
                  WRITE (LUPRI,'(2X,A,I3,5X,A,I3)') 'Component: ',
     &                 ICOMP, 'Coordinate direction: ',ICOOR
                  CALL OUTPAK(WORK(IPOINT),NBAST,1,LUPRI)
 101           CONTINUE
            END IF
C
C     Construct and -2g(l)sum(m) T(lm) t(lm)a_mn
C
            IL = 0
            DO 200 M = -L, L
               ICOUNT = ICOUNT + 1
               IL = IL + 1
               IPOINT = KTMP + (IL - 1)*NNBASX
               FAC2 = -D2*GLM(ICOUNT)*TLM(ICOUNT)
               IOFFD = 1 + (ICOOR - 1)*NNBASX
               CALL DAXPY(NNBASX,FAC2,WORK(IPOINT),1,
     &                    DERMAT(IOFFD),1)
 200        CONTINUE
 120     CONTINUE 
 140  CONTINUE 
C
C        We need to square things away, in particular certain matrices
C
      CALL DCOPY(NCOMP*NNBASX,DERMAT,1,WORK(1),1)
      DO 300 ICOMP = 1, NCOMP
         IOFFD = 1 + (ICOMP - 1)*NNBASX
         IOFFS = 1 + (ICOMP - 1)*N2BASX
         CALL DAPTGE(NBAST,WORK(IOFFD),DERMAT(IOFFS))
 300  CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL AROUND('Solvent contributions to gradient in MAGSO2')
         IOFF = 1
         WRITE (LUPRI,'(2X,A)') 'X coordinate'
         CALL OUTPUT(DERMAT(IOFF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IOFF = IOFF + N2BASX
         WRITE (LUPRI,'(2X,A)') 'Y coordinate'
         CALL OUTPUT(DERMAT(IOFF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IOFF = IOFF + N2BASX
         WRITE (LUPRI,'(2X,A)') 'Z coordinate'
         CALL OUTPUT(DERMAT(IOFF),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     We ought to be more or less finished now. Change EOF label in order
C     to overwrite the integrals
C
      CALL GPCLOSE(LUSOL,'KEEP')
      LABEL = 'X 000000'
      REWIND (LUPROP)
      IF(.NOT.FNDLAB(LABEL,LUPROP)) THEN
         WRITE (LUPRI,'(/A)') ' Integral label '//LABEL//
     &        ' not found on file AOPROPER.'
         CALL QUIT('Integral label not found in hersol')
      END IF
      BACKSPACE LUPROP
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      CALL NEWLB2('EOFLABEL',RTNLBL,LUPROP,LUPRI)
      CALL GPCLOSE(LUPROP,'KEEP')
      RETURN
      END
C  /* Deck magsolexp */
      SUBROUTINE MAGSOLEXP(EXPVAL,EXPSPH,EXPCAR,INTREP,INTADR,LABINT,
     &                     GLM,TLM,TRAMAT,INTPRI,DENMAT,WORK,LWORK)
C
C     K.Ruud, July-05
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <iratdef.h>
      PARAMETER (MORDER = 1, MINTEG=0, D2 = 2.0D0)
      LOGICAL EXP1VL
      DIMENSION GLM(LMTOT), TLM(LMTOT),
     &          EXPVAL(6), EXPSPH(6,-LCAVMX:LCAVMX),
     &          EXPCAR(6,(LCAVMX + 1)*(LCAVMX + 2)/2),
     &          INTREP(*), INTADR(*),
     &          TRAMAT((LCAVMX + 1)*(LCAVMX + 2)*(2*LCAVMX + 1)/2),
     &          DENMAT(*), WORK(LWORK)
      CHARACTER*8 LABINT(*)
#include <symmet.h>
#include <cbisol.h>
#include <orgcom.h>
#include <mgsolt.h>
#include <inforb.h>
#include <infinp.h>
#include <inftap.h>
C
      DIMENSION CAVCNT(3)
      LOGICAL   FOUND
C
C     Construct GLM, TLM, and density matrix
C
      CALL SOLFL(GLM,EPDIEL,RCAV,LCAVMX)
      CALL RD_SIRIFC('SOLTLM',FOUND,TLM,WORK,LWORK)
      IF (.NOT. FOUND)
     &    CALL QUIT('MAGSO1 ERROR: "SOLTLM" not found on SIRIFC.')
C
      CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUPRI)
      READ (LUSOL) LMAXSS, LMTOTX, NNNBAS, CAVCNT
      IF (LMAXSS .LT. LSOLMX) THEN
         WRITE (LUPRI,'(//2A,2(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' insufficient number of intgrals on LUSOL',
     *   ' l max from SIRIUS input :',LSOLMX,
     *   ' l max from LUSOL  file  :',LMAXSS
         CALL QUIT('MAGSO1: lmax on LUSOL is too small')
      END IF
      IF ((LMAXSS+1)**2 .NE. LMTOTX) THEN
         WRITE (LUPRI,'(//2A,3(/A,I5))') ' >>> MAGSO1 ERROR,',
     *   ' LUSOL file info inconsistent',
     *   ' l_max               :',LMAXSS,
     *   ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *   ' LMTOTX              :',LMTOTX
         CALL QUIT('MAGSO1: LUSOL info not internally consistent')
      END IF
C
      READ (LUSOL)
      LM = 0
      ICOUNT = 0
      DOALL_MGSOLT  = .TRUE.
      EXP1VL = .TRUE.
      CALL DZERO(EXPVAL,6)
C
      NCOMP = 6
      DO 100 L = 0, LSOLMX
         CALL DZERO(EXPSPH,(2*LCAVMX + 1)*6)
         CALL DZERO(EXPCAR,(LCAVMX + 1)*(LCAVMX + 2)*3)
         MCOMP = (L + 1)*(L + 2)*NCOMP/2
C
C     Construct differentiated integrals in spherical harmonic basis
C
         NPATOM = 0
         JCOMP  = 0
         CALL GET1IN(DUMMY,'LONSOL2',JCOMP,WORK,LWORK,LABINT,INTREP,
     &               INTADR,L,.FALSE.,NPATOM,.TRUE.,EXPCAR,EXP1VL,
     &               DENMAT,INTPRI)
         CALL SPHCOM(L,TRAMAT,2*L + 1,(L + 1)*(L + 2)/2,MORDER,MINTEG,
     &               WORK,LWORK,INTPRI)
         DO 430 ICOOR = 1, NCOMP
            DO 400 ICAR = 1, (L + 1)*(L + 2)/2
               IL = 0
               DO 410 ISPH = -L, L
                  IL = IL + 1
                  IELMNT = (L + 1)*(L + 2)*(IL - 1)/2 + ICAR
                  COEF = TRAMAT(IELMNT)
                  EXPSPH(ICOOR,ISPH) = EXPSPH(ICOOR,ISPH)
     &                               + EXPCAR(ICOOR,ICAR)*COEF
 410           CONTINUE
 400        CONTINUE
 430     CONTINUE
C
C     Construct and -2g(l)sum(m) T(lm) t(lm)a_mn
C
         DO 200 M = -L, L
            ICOUNT = ICOUNT + 1
            DO 190 ICOOR = 1, NCOMP
               FAC2 = -D2*GLM(ICOUNT)*TLM(ICOUNT)
               EXPVAL(ICOOR) = EXPVAL(ICOOR) + EXPSPH(ICOOR,M)*FAC2
 190        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
C     We ought to be more or less finished now
C
      CALL GPCLOSE(LUSOL,'KEEP')
      RETURN
      END
#endif /* ifndef PRG_DIRAC */
! -- end of hersol.F --
