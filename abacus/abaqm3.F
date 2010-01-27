C/* deck qm3first */
      SUBROUTINE QM3FIRST(VEC2,NCOMP,IPRINT,INTPRI,WORK,LWORK)
C
c This routine adds the MM solvent contribution 
c to the 1st order (NCOMP=3) magnetic field perturbation due to the use of
c London orbitals
c Copenhagen, January 2006, K.Ruud
c
c VEC2 (OUTPUT)  : contribution to be added to perturbed Fock matrix
c NCOMP (INPUT)  : number of independent perturbation (3)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include "nuclei.h"
#include "qm3.h"
#include <inforb.h>
C#include <infpri.h>
C#include <infvar.h>
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION WORK(*),VEC2(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C
C     Construct charges for all MM centers
C
      L = 0
      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)
      CALL DZERO(VEC2,3*N2BASX)
      DO I = 1, ISYTP
         IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO J = NSYSBG(I), NSYSED(I)
               DO K = 1,NSISY(I)
                  L = L + 1
                  KMAT = 1
                  KLAST = KMAT + 3*N2BASX
                  LWRK = LWORK - KLAST + 1
                  IATNOW = NUCIND + L
C
                  KPATOM = 0
                  NOSIM = 3
                  TOFILE = .FALSE.
                  TRIMAT = .FALSE.
                  EXP1VL = .FALSE.
                  DIPORG(1) = CORD(1,IATNOW)
                  DIPORG(2) = CORD(2,IATNOW)
                  DIPORG(3) = CORD(3,IATNOW)
                  CALL GET1IN(WORK(KMAT),'PCMBSOL',NOSIM,WORK(KLAST),
     &                        LWRK,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                        KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)

                  IF (IQM3PR .GE. 12) THEN
                     WRITE (LUPRI,'(/A,I3,A)')
     *                    ' N(',L,')_ao matrix (Bx): '
                     CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A,I3,A)')
     *                    ' N(',L,')_ao matrix (By): '
                     CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A,I3,A)')
     *                    ' N(',L,')_ao matrix (Bz): '
                     CALL OUTPUT(WORK(KMAT+2*N2BASX),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                  END IF
C  
                  FAC1 = -1.0D0*CHARGE(IATNOW)
                  CALL DAXPY(3*N2BASX,FAC1,WORK(KMAT),1,VEC2,1)
               END DO
            END DO
         END IF
      END DO
C
C     Calculate MM corrections for all distributed polarizabilities
C
C     Remember that the first coordinate for the polarizability is the COM for the QM system.
C     Therefor we increase L with 1. If however we have distributed polarizability for the QM 
C     system we need to skip the first NUALIS(0) points

      L = L + NUALIS(0)

      LM = 0
      KRAX   = 1
      KRAY   = KRAX   + NCOMS
      KRAZ   = KRAY   + NCOMS
      KENSAX = KRAZ   + NCOMS
      KENSAY = KENSAX + NCOMS
      KENSAZ = KENSAY + NCOMS
      KLST   = KENSAZ + NCOMS
      CALL DZERO(WORK(KRAX),6*NCOMS)
C
      CALL CC_GET31('CC_RA',NCOMS,WORK(KRAX),
     &              WORK(KRAY),WORK(KRAZ))
      CALL CC_GET31('ENSAFILE',NCOMS,WORK(KENSAX),
     &              WORK(KENSAY),WORK(KENSAZ))
C
      DO I = 1, ISYTP
         IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO J = NSYSBG(I), NSYSED(I)
               DO K = 1, NUALIS(I)
                  LM = LM + 1
                  KMAT = KLST
                  KLAST = KMAT + 9*N2BASX
                  LWRK = LWORK - KLAST + 1
                  IATNOW = NUCIND + L + LM
C
                  KPATOM = 0
                  NOSIM = 9
                  TOFILE = .FALSE.
                  TRIMAT = .FALSE.
                  EXP1VL = .FALSE.
                  DIPORG(1) = CORD(1,IATNOW)
                  DIPORG(2) = CORD(2,IATNOW)
                  DIPORG(3) = CORD(3,IATNOW)

                  CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                        LWRK,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                        KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)
                  IF (IQM3PR .GE. 15) THEN
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ex Bx matrix:'
                     CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ex By matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ex Bz matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ey Bx matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ey By matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ey Bz matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ez Bx matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ez By matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra_ao Ez Bz matrix:'
                     CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                           NBAST,NBAST,1,LUPRI)
                  END IF

                  IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
                     FACx = -ALPIMM(I,K)*(WORK(KRAX + LM - 1)
     &                    + 0.5D0 * WORK(KENSAX + LM - 1))
                     FACy = -ALPIMM(I,K)*(WORK(KRAY + LM - 1)
     &                    + 0.5D0 * WORK(KENSAY + LM - 1))
                     FACz = -ALPIMM(I,K)*(WORK(KRAZ + LM - 1)
     &                    + 0.5D0 * WORK(KENSAZ + LM - 1))
                     CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,VEC2,1)
                     CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &                          VEC2(1+N2BASX),1)
                     CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &                          VEC2(1+2*N2BASX),1)
                     CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &                          VEC2,1)
                     CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &                          VEC2(1+N2BASX),1)
                     CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &                          VEC2(1+2*N2BASX),1)
                     CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &                          VEC2,1)
                     CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &                          VEC2(1+N2BASX),1)
                     CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &                          VEC2(1+2*N2BASX),1)
                  END IF
               END DO
            END DO
         END IF
      END DO
C
C     **print out section**
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND(
     &        'First order solvent contributions to gradient in MAGQM3')
         WRITE (LUPRI,'(2X,A)') 'X coordinate'
         CALL OUTPUT(VEC2(1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Y coordinate'
         CALL OUTPUT(VEC2(1+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &        LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Z coordinate'
         CALL OUTPUT(VEC2(1+2*N2BASX),1,NBAST,1,NBAST,
     &        NBAST,NBAST,1,LUPRI)
      END IF
C
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
      RETURN
      END

C/* deck qm3first_p */
      SUBROUTINE QM3FIRST_P(VEC2,NCOMP,IPRINT,INTPRI,WORK,LWORK)
C
C The parallel version of QM3FIRST, Tromso Oct 08, Arnfinn
C
c This routine adds the MM solvent contribution 
c to the 1st order (NCOMP=3) magnetic field perturbation due to the use of
c London orbitals
c Copenhagen, January 2006, K.Ruud
c
c VEC2 (OUTPUT)  : contribution to be added to perturbed Fock matrix
c NCOMP (INPUT)  : number of independent perturbation (3)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include "nuclei.h"
#include "qm3.h"
#include <inforb.h>
Cahs#include <infpar.h>
Cahs#include <mtags.h>
Cahs#if defined(VAR_MPI)
Cahs#include <mpif.h>
Cahs#endif
C#include <infpri.h>
C#include <infvar.h>
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION WORK(*),VEC2(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C
      CALL QENTER('QM3FIRST_P')
C
C     Construct charges for all MM centers
C
      L = 0
      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)
      CALL DZERO(VEC2,3*N2BASX)
      DO I = 1, ISYTP
         IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            CALL QM3FIRST_M1(VEC2,IATNOW,NOSIM,I,L,
     &                       WORK(1),LWORK)
         ENDIF
      ENDDO
C
C     Calculate MM corrections for all distributed polarizabilities
C
C     Remember that the first coordinate for the polarizability is the COM for the QM system.
C     Therefor we increase L with 1. If however we have distributed polarizability for the QM 
C     system we need to skip the first NUALIS(0) points
C
      L = L + NUALIS(0)
      LM = 0
      KRAX   = 1
      KRAY   = KRAX   + NCOMS
      KRAZ   = KRAY   + NCOMS
      KENSAX = KRAZ   + NCOMS
      KENSAY = KENSAX + NCOMS
      KENSAZ = KENSAY + NCOMS
      KLST   = KENSAZ + NCOMS
      LWORK1 = LWORK - KLST +1
      CALL DZERO(WORK(KRAX),6*NCOMS)
C
      CALL CC_GET31('CC_RA',NCOMS,WORK(KRAX),
     &              WORK(KRAY),WORK(KRAZ))
      CALL CC_GET31('ENSAFILE',NCOMS,WORK(KENSAX),
     &              WORK(KENSAY),WORK(KENSAZ))
C
      DO I = 1, ISYTP
         IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            CALL QM3FIRST_M2(WORK(KRAX),WORK(KRAY),WORK(KRAZ),
     &                       WORK(KENSAX),WORK(KENSAY),WORK(KENSAZ),
     &                       LM,VEC2,IATNOW,NOSIM,I,L,WORK(KLST),LWORK1)
         ENDIF
      ENDDO
C
C     **print out section**
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND(
     &        'First order solvent contributions to gradient in MAGQM3')
         WRITE (LUPRI,'(2X,A)') 'X coordinate'
         CALL OUTPUT(VEC2(1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Y coordinate'
         CALL OUTPUT(VEC2(1+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &        LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Z coordinate'
         CALL OUTPUT(VEC2(1+2*N2BASX),1,NBAST,1,NBAST,
     &        NBAST,NBAST,1,LUPRI)
      END IF
C
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
      CALL QEXIT('QM3FIRST_P')
C
      RETURN
      END
C
C-----------------------------------------------------
C
      SUBROUTINE QM3FIRST_M1(VEC2,IATNOW,NOSIM,INUM,LNUM,WORK,LWORK)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include "nuclei.h"
#include "qm3.h"
#include <inforb.h>
#include <magone.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      INTEGER IATNOW, NOSIM
      DIMENSION VEC2(3*N2BASX)
      DIMENSION WORK(LWORK)
      DIMENSION ISX(MXCORB),IPATOM(MXCENT)
C
      CALL QENTER('QM3FIRST_M1')
C
      IF (TOFILE) CALL QUIT('Parallel calculations do not allow '//
     &                     'for storing integrals on disk')

      KVEC2 = 1
      JVEC2 = KVEC2 + 3*N2BASX
      KLAST = JVEC2 + 3*N2BASX
      IPRTYP = 19
C
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)
C
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)

      DO J = NSYSBG(INUM), NSYSED(INUM)
         DO K = 1, NSISY(INUM)
            LNUM = LNUM + 1
            IWHO = -1
            CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
            CALL MPIXSEND(LNUM, 1, 'INTEGER', NWHO, MPTAG2)
         END DO
      END DO

C     Send end message to all slaves

      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves

      CALL DZERO(WORK(KVEC2),3*N2BASX)
      CALL MPI_REDUCE(WORK(KVEC2),WORK(JVEC2),3*N2BASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      FAC = 1.0D0
      CALL DAXPY(3*N2BASX,FAC,WORK(JVEC2),1,VEC2,1)

      CALL QEXIT('QM3FIRST_M1')

      RETURN
      END


C-----------------------------------------------------
      SUBROUTINE QM3FIRST_S1(WRK,LWRK,IPRTMP)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include "nuclei.h"
#include "qm3.h"
#include <inforb.h>
#include <magone.h>

#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)
      INTEGER IATNOW, NOSIM
      DIMENSION ISX(MXCORB),IPATOM(MXCENT)

      CALL QENTER('QM3FIRST_S1')

      KVEC2  = 1
      KLAST = KVEC2 + 3*N2BASX

      CALL DZERO(WRK(KVEC2),3*N2BASX)

      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)

 187  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(L, 1, 'INTEGER', MASTER, MPTAG2)

      IF (L .GT. 0) THEN

         KMAT = KLAST
         KLAST = KMAT + 3*N2BASX
         LWRK1 = LWRK - KLAST + 1
         IATNOW = NUCIND + L
C
         KPATOM = 0
         NOSIM = 3
         TOFILE = .FALSE.
         TRIMAT = .FALSE.
         EXP1VL = .FALSE.
         DIPORG(1) = CORD(1,IATNOW)
         DIPORG(2) = CORD(2,IATNOW)
         DIPORG(3) = CORD(3,IATNOW)
         CALL GET1IN(WRK(KMAT),'PCMBSOL',NOSIM,WRK(KLAST),
     &        LWRK1,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &        KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)

         IF (IQM3PR .GE. 12) THEN
            WRITE (LUPRI,'(/A,I3,A)')
     *           ' N(',L,')_ao matrix (Bx): '
            CALL OUTPUT(WRK(KMAT),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A,I3,A)')
     *           ' N(',L,')_ao matrix (By): '
            CALL OUTPUT(WRK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A,I3,A)')
     *           ' N(',L,')_ao matrix (Bz): '
            CALL OUTPUT(WRK(KMAT+2*N2BASX),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
         END IF
C  
         FAC1 = -1.0D0*CHARGE(IATNOW)
         CALL DAXPY(3*N2BASX,FAC1,WRK(KMAT),1,WRK(KVEC2),1)
         GO TO 187
      END IF


C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KVEC2),WRK(KVEC2),3*N2BASX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QM3FIRST_S1')

      RETURN
      END

C-----------------------------------------------------
      SUBROUTINE QM3FIRST_M2(RAX,RAY,RAZ,ENSAX,ENSAY,ENSAZ,LM,
     &                       VEC2,IATNOW,NOSIM,INUM,LNUM,WORK,LWORK)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include "nuclei.h"
#include "qm3.h"
#include <inforb.h>
#include <magone.h>

#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      INTEGER IATNOW, NOSIM

      DIMENSION VEC2(*)
      DIMENSION RAX(*),RAY(*),RAZ(*),ENSAX(*),ENSAY(*),ENSAZ(*)
      DIMENSION WORK(LWORK)
      DIMENSION ISX(MXCORB),IPATOM(MXCENT)

      CALL QENTER('QM3FIRST_M2')

      IF (TOFILE) CALL QUIT('Parallel calculations do not allow '//
     &                     'for storing integrals on disk')

      KVEC2 = 1
      JVEC2 = KVEC2 + 3*N2BASX
      KLAST = JVEC2 + 3*N2BASX

      IPRTYP = 20

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NCOMS,1,'INTEGER',MASTER)
      CALL MPIXBCAST(RAX,NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(RAY,NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(RAZ,NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(ENSAX,NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(ENSAY,NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(ENSAZ,NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LNUM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MDLWRD(INUM),7,'CHARACTER',MASTER)

      DO J = NSYSBG(INUM), NSYSED(INUM)
         DO K = 1, NUALIS(INUM)
            LM = LM + 1
            FAC = -ALPIMM(INUM,K)
            IWHO = -1
            CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
            CALL MPIXSEND(LM, 1, 'INTEGER', NWHO, MPTAG2)
            CALL MPIXSEND(FAC, 1, 'DOUBLE', NWHO, MPTAG2)
         END DO
      END DO

C     Send end message to all slaves

      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND, 1, 'INTEGER', NWHO, MPTAG2)
         CALL MPIXSEND(FAC, 1, 'DOUBLE', NWHO, MPTAG2)
      END DO

C     Collect data from all slaves
      CALL DZERO(WORK(KVEC2),3*N2BASX)
      CALL MPI_REDUCE(WORK(KVEC2),WORK(JVEC2),3*N2BASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      FAC = 1.0D0
      CALL DAXPY(3*N2BASX,FAC,WORK(JVEC2),1,VEC2,1)

      CALL QEXIT('QM3FIRST_M2')

      RETURN
      END

C-----------------------------------------------------
      SUBROUTINE QM3FIRST_S2(WRK,LWRK,IPRTMP)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include "nuclei.h"
#include "qm3.h"
#include <inforb.h>
#include <magone.h>

#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)
      INTEGER IATNOW, NOSIM
      DIMENSION ISX(MXCORB),IPATOM(MXCENT)
      CHARACTER MDLWORD*7

      CALL QENTER('QM3FIRST_S2')

      CALL MPIXBCAST(NCOMS,1,'INTEGER',MASTER)

      KVEC2 = 1
      KRAX   = KVEC2 + 3*N2BASX
      KRAY   = KRAX + NCOMS
      KRAZ   = KRAY + NCOMS
      KENSAX = KRAZ + NCOMS
      KENSAY = KENSAX + NCOMS
      KENSAZ = KENSAY + NCOMS
      KLAST  = KENSAZ + NCOMS

      CALL DZERO(WRK(KVEC2),3*N2BASX)

      CALL MPIXBCAST(WRK(KRAX),NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KRAY),NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KRAZ),NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KENSAX),NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KENSAY),NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KENSAZ),NCOMS,'DOUBLE',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(L,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(MDLWORD,7,'CHARACTER',MASTER)

 187  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(LM, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(FAC, 1, 'DOUBLE', MASTER, MPTAG2)

      IF (LM .GT. 0) THEN
         KMAT  = KLAST
         KLAST = KMAT + 9*N2BASX
         LWRK1 = LWRK - KLAST + 1
         IATNOW = NUCIND + L + LM
C
         KPATOM = 0
         NOSIM = 9
         TOFILE = .FALSE.
         TRIMAT = .FALSE.
         EXP1VL = .FALSE.
         DIPORG(1) = CORD(1,IATNOW)
         DIPORG(2) = CORD(2,IATNOW)
         DIPORG(3) = CORD(3,IATNOW)
         CALL GET1IN(WRK(KMAT),'EFIELB1',NOSIM,WRK(KLAST),
     &        LWRK1,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &        KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)
         IF (IQM3PR .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' Rra_ao Ex Bx matrix:'
            CALL OUTPUT(WRK(KMAT),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ex By matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ex Bz matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ey Bx matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ey By matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ey Bz matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ez Bx matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ez By matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rra_ao Ez Bz matrix:'
            CALL OUTPUT(WRK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &           NBAST,NBAST,1,LUPRI)
         END IF
         IF (MDLWORD .EQ. 'SPC_E01') THEN
            FACx = FAC*(WRK(KRAX + LM - 1)
     &           + 0.5D0 * WRK(KENSAX + LM - 1))
            FACy = FAC*(WRK(KRAY + LM - 1)
     &           + 0.5D0 * WRK(KENSAY + LM - 1))
            FACz = FAC*(WRK(KRAZ + LM - 1)
     &           + 0.5D0 * WRK(KENSAZ + LM - 1))
            CALL DAXPY(N2BASX,FACX,WRK(KMAT),1,WRK(KVEC2),1)
            CALL DAXPY(N2BASX,FACX,WRK(KMAT+N2BASX),1,
     &           WRK(KVEC2+N2BASX),1)
            CALL DAXPY(N2BASX,FACX,WRK(KMAT+2*N2BASX),1,
     &           WRK(KVEC2+2*N2BASX),1)
            CALL DAXPY(N2BASX,FACY,WRK(KMAT+3*N2BASX),1,
     &           WRK(KVEC2),1)
            CALL DAXPY(N2BASX,FACY,WRK(KMAT+4*N2BASX),1,
     &           WRK(KVEC2+N2BASX),1)
            CALL DAXPY(N2BASX,FACY,WRK(KMAT+5*N2BASX),1,
     &           WRK(KVEC2+2*N2BASX),1)
            CALL DAXPY(N2BASX,FACZ,WRK(KMAT+6*N2BASX),1,
     &           WRK(KVEC2),1)
            CALL DAXPY(N2BASX,FACZ,WRK(KMAT+7*N2BASX),1,
     &           WRK(KVEC2+N2BASX),1)
            CALL DAXPY(N2BASX,FACZ,WRK(KMAT+8*N2BASX),1,
     &           WRK(KVEC2+2*N2BASX),1)
         END IF
         GO TO 187
      END IF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KVEC2),WRK(KVEC2),3*N2BASX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL QEXIT('QM3FIRST_S2')

      RETURN
      END
C/* deck qm3b2 */
      SUBROUTINE QM3B2(EXPVAL,DENMAT,WORK,LWORK)
C
c This routine adds the MM solvent contribution 
c to the 2nd-order magnetic field perturbation due to the use of
c London orbitals
c Copenhagen, January 2006, K.Ruud
c
c EXPVAL (OUTPUT): contribution to be added to diamagnetic magnetizability
c
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include "nuclei.h"
#include "qm3.h"
#include <inforb.h>
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION WORK(*), EXPVAL(6), DENMAT(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C
C     Construct charges for all MM centers
C
      L = 0
      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)
      CALL DZERO(EXPVAL,6)
      DO I = 1, ISYTP
         IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO J = NSYSBG(I), NSYSED(I)
               DO K = 1,NSISY(I)
                  L = L + 1
                  KTMP = 1
                  KLAST = KTMP + 6
                  LWRK = LWORK - KLAST + 1
                  IATNOW = NUCIND + L
C
                  KPATOM = 0
                  NOSIM = 0
                  TOFILE = .FALSE.
                  TRIMAT = .FALSE.
                  EXP1VL = .TRUE.
                  DIPORG(1) = CORD(1,IATNOW)
                  DIPORG(2) = CORD(2,IATNOW)
                  DIPORG(3) = CORD(3,IATNOW)
                  CALL GET1IN(DUMMY,'PCMB2SL',NOSIM,WORK(KLAST),
     &                        LWRK,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                        KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                        IQM3PR)
                  IF (IQM3PR .GE. 4) THEN
                     WRITE (LUPRI,'(/A,I3,A)')
     *                    ' N(',L,') correction matrix '
                     CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
                  END IF
C  
                  FAC1 = -1.0D0*CHARGE(IATNOW)
                  CALL DAXPY(6,FAC1,WORK(KTMP),1,EXPVAL,1)
               END DO
            END DO
         END IF
      END DO
C
C     Calculate MM corrections for all distributed polarizabilities
C
C     Remember that the first coordinate for the polarizability is the COM for the QM system.
C     Therefor we increase L with 1. If however we have distributed polarizability for the QM
C     system we need to skip the first NUALIS(0) points

      L = L + NUALIS(0)

      LM = 0
      KRAX   = 1
      KRAY   = KRAX   + NCOMS
      KRAZ   = KRAY   + NCOMS
      KENSAX = KRAZ   + NCOMS
      KENSAY = KENSAX + NCOMS
      KENSAZ = KENSAY + NCOMS
      KLST   = KENSAZ + NCOMS
      CALL DZERO(WORK(KRAX),6*NCOMS)
C
      CALL CC_GET31('CC_RA',NCOMS,WORK(KRAX),
     &              WORK(KRAY),WORK(KRAZ))
      CALL CC_GET31('ENSAFILE',NCOMS,WORK(KENSAX),
     &              WORK(KENSAY),WORK(KENSAZ))
C
      DO I = 1, ISYTP
         IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO J = NSYSBG(I), NSYSED(I)
               DO K = 1, NUALIS(I)
                  LM = LM + 1
                  KTMP = KLST
                  KLAST = KTMP + 18
                  LWRK = LWORK - KLAST + 1
                  IATNOW = NUCIND + L + LM
C
                  KPATOM = 0
                  NOSIM = 0
                  TOFILE = .FALSE.
                  TRIMAT = .FALSE.
                  EXP1VL = .TRUE.
                  DIPORG(1) = CORD(1,IATNOW)
                  DIPORG(2) = CORD(2,IATNOW)
                  DIPORG(3) = CORD(3,IATNOW)
                  CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                        LWRK,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                        KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                        IQM3PR)
                  IF (IQM3PR .GE. 4) THEN
                     WRITE (LUPRI,'(/A)') ' Rra Ex 2.order matrix:'
                     CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra Ey 2.order matrix:'
                     CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
                     WRITE (LUPRI,'(/A)') ' Rra Ez 2.order matrix:'
                     CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
                  END IF
C
                  IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
                     FACx = -ALPIMM(I,K)*(WORK(KRAX + LM - 1)
     &                    + 0.5D0 * WORK(KENSAX + LM - 1))
                     FACy = -ALPIMM(I,K)*(WORK(KRAY + LM - 1)
     &                    + 0.5D0 * WORK(KENSAY + LM - 1))
                     FACz = -ALPIMM(I,K)*(WORK(KRAZ + LM - 1)
     &                    + 0.5D0 * WORK(KENSAZ + LM - 1))
                     DO I2 = 1, 6
                        EXPVAL(I2) = EXPVAL(I2) + WORK(KTMP+I2-1)*FACX +
     &                      WORK(KTMP+I2+5)*FACY + WORK(KTMP+I2+11)*FACZ
                     END DO
                  END IF
               END DO
            END DO
         END IF
      END DO
C
C     **print out section**
C
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
      RETURN
      END
C-------------------------------------------------------------------------------
C*******************************************************************
C/* deck qmmmfirst */
      SUBROUTINE QMMMFIRST(VEC2,IPRINT,INTPRI,WORK,LWORK)
C
c This routine adds the MM environmental contribution 
c to the 1st order (NCOMP=3) magnetic field perturbation due to the use of
c London orbitals
c Lund, Nov. 2008, J. Kongsted
c Based on QM3FIRST by K. Ruud
c
c VEC2 (OUTPUT)  : contribution to be added to perturbed Fock matrix
c NCOMP (INPUT)  : number of independent perturbation (3)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include <nuclei.h>
#include <qm3.h>
#include <qmmm.h>
#include <inforb.h>
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL,LOCDEB
      LOGICAL LCOMB
      DIMENSION WORK(*),VEC2(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      LOCDEB = .FALSE.

C     If we have BOTH permanent dipoles AND induced dipoles then these
C     can be added together and we can treat the effect of those
C     simultanously. Then no seperate loop for dipoles and
C     polarizabilities but only one combined loop.

      LCOMB = .FALSE.
      IF ( (IPOLTP .GT. 0) .AND. (NMULT .GE. 1) ) LCOMB = .TRUE.

C     SPLNMR overrules these settings   
      IF (SPLNMR) LCOMB = .FALSE.

C     We don't have these integrals yet !     
      IF (NMULT .GE. 2) THEN
        CALL QUIT('Quadrupole contr. not implemented in QMMMFIRST')     
      ENDIF

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(VEC2,3*N2BASX)

C     1. Construct the charge contribution for all MM centers

      IF (NMULT .GE. 0) THEN

        KMAT = 1
        KLAST = KMAT + 3*N2BASX
        LWRK = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),3*NNBASX)

        DO 100 I = 1,MMCENT

           IF (ABS(MUL0MM(I)) .LE. THRMM) GOTO 100

           KPATOM = 0
           NOSIM = 3
           TOFILE = .FALSE.
           TRIMAT = .FALSE.
           EXP1VL = .FALSE.
           DIPORG(1) = MMCORD(1,I)
           DIPORG(2) = MMCORD(2,I)
           DIPORG(3) = MMCORD(3,I)

           CALL DZERO(WORK(KMAT),3*NNBASX)

           CALL GET1IN(WORK(KMAT),'PCMBSOL',NOSIM,WORK(KLAST),
     &                 LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                 KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)

           IF ( (INTPRI .GE. 12) .OR. (LOCDEB) ) THEN
              WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (Bx): '
              CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
              WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (By): '
              CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
              WRITE (LUPRI,'(/A,I3,A)')
     *             ' N(',I,')_ao matrix (Bz): '
              CALL OUTPUT(WORK(KMAT+2*N2BASX),1,NBAST,1,NBAST,
     &                    NBAST,NBAST,1,LUPRI)
           END IF

           FAC1 = -1.0D0*MUL0MM(I)
           CALL DAXPY(3*N2BASX,FAC1,WORK(KMAT),1,VEC2,1)

  100   CONTINUE

      ENDIF !( NMULT .GE. 0)

C     2. Construct the permanent dipole contribution for all MM centers

      IF ( (NMULT .GE. 1) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from permanent dipoles'
        ENDIF

        KMAT = 1
        KLAST = KMAT + 9*N2BASX
        LWRK = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),9*NNBASX)

        DO 200 I = 1,MMCENT

          DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2
          DNORM = SQRT(DNORM2)

          IF (ABS(DNORM) .LE. THRMM) GOTO 200

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)
          IF ( (INTPRI .GE. 15) .OR. (LOCDEB) ) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*MUL1MM(1,I)
          FACY = -1.0D0*MUL1MM(2,I)
          FACZ = -1.0D0*MUL1MM(3,I)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,VEC2,1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)

  200   CONTINUE

      ENDIF !( NMULT .GE. 1)


C     3. Construct the induced dipole contribution for all MM centers
C        NOTE: This could be done more effectively by adding the
C        permanent dipoles and induced dipoles and then performing the
C        London orbital contribution. See the LCOMB option below.

      IF ( (IPOLTP .GT. 0) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from induced dipoles'
        ENDIF

        KMAT    = 1
        KINDMOM = KMAT    + 9*N2BASX ! INDMOM contains the ind. dips.
        KLAST   = KINDMOM + 3*NNZAL
        LWRK    = LWORK   - KLAST + 1 

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),9*NNBASX)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMFIRST'
           WRITE(LUPRI,*)
        ENDIF 

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

        IL = 0 ! Counter for induced dipole moments

        DO 300 I = 1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 300

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)
          IF ( (INTPRI .GE. 15) .OR. (LOCDEB) ) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KINDMOM+IL+0)
          FACY = -1.0D0*WORK(KINDMOM+IL+1)
          FACZ = -1.0D0*WORK(KINDMOM+IL+2)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,VEC2,1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)

          IL = IL+3

  300   CONTINUE

      ENDIF !(IPOLTP .GT. 0)

C     3a. If both permanent dipoles AND polarizabilities construct TOTAL 
C         dipole contribution for all MM centers

      IF (LCOMB)  THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from total dipoles'
        ENDIF

        KMAT    = 1
        KINDMOM = KMAT    + 9*N2BASX ! INDMOM contains the ind. dips.
        KTOTMOM = KINDMOM + 3*NNZAL  ! TOTMOM contains the TOTAL dips.
        KLAST   = KTOTMOM + 3*MMCENT
        LWRK    = LWORK   - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMFIRST',-KLAST,LWORK)

        CALL DZERO(WORK(KMAT),9*NNBASX)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)
        CALL DZERO(WORK(KTOTMOM),3*MMCENT)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMFIRST'
           WRITE(LUPRI,*)
        ENDIF 

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

        IL = 0
        IJ = 0
        DO 400 I = 1,MMCENT
          IF (ZEROAL(I) .EQ. -1) THEN
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I)
          ELSE IF (ZEROAL(I) .EQ. 1) THEN 
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I) + WORK(KINDMOM+IJ+0)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I) + WORK(KINDMOM+IJ+1)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I) + WORK(KINDMOM+IJ+2)
            IJ=IJ+3
          ENDIF
          IL=IL+3
  400   CONTINUE

        DO 500 I = 1,MMCENT
           IL = 3*(I - 1)
C         See if the total dipole moment at this site is zero
          DNORM2 = WORK(KTOTMOM+IL+0)*WORK(KTOTMOM+IL+0)
     &           + WORK(KTOTMOM+IL+1)*WORK(KTOTMOM+IL+1)
     &           + WORK(KTOTMOM+IL+2)*WORK(KTOTMOM+IL+2)
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) GOTO 500

          CALL DZERO(WORK(KMAT),9*NNBASX)

          KPATOM = 0
          NOSIM = 9
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .FALSE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(WORK(KMAT),'EFIELB1',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,INTPRI)
          IF ( (INTPRI .GE. 15) .OR. (LOCDEB) ) THEN
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bx matrix:'
             CALL OUTPUT(WORK(KMAT),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ex Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*2),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*3),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*4),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ey Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*5),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bx matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*6),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez By matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*7),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield_ao Ez Bz matrix:'
             CALL OUTPUT(WORK(KMAT+N2BASX*8),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KTOTMOM+IL+0)
          FACY = -1.0D0*WORK(KTOTMOM+IL+1)
          FACZ = -1.0D0*WORK(KTOTMOM+IL+2)

          CALL DAXPY(N2BASX,FACX,WORK(KMAT),1,VEC2,1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACX,WORK(KMAT+2*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+3*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+4*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACY,WORK(KMAT+5*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+6*N2BASX),1,
     &               VEC2,1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+7*N2BASX),1,
     &               VEC2(1+N2BASX),1)
          CALL DAXPY(N2BASX,FACZ,WORK(KMAT+8*N2BASX),1,
     &               VEC2(1+2*N2BASX),1)

  500   CONTINUE

      ENDIF !(LCOMB)

C     **print out section**
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND(
     &        'First order qmmm contr. to gradient in QMMMFIRST')
         WRITE (LUPRI,'(2X,A)') 'X coordinate'
         CALL OUTPUT(VEC2(1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Y coordinate'
         CALL OUTPUT(VEC2(1+N2BASX),1,NBAST,1,NBAST,NBAST,NBAST,1,
     &        LUPRI)
         WRITE (LUPRI,'(2X,A)') 'Z coordinate'
         CALL OUTPUT(VEC2(1+2*N2BASX),1,NBAST,1,NBAST,
     &        NBAST,NBAST,1,LUPRI)
      END IF
C
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
      RETURN
      END
C*****************************************************************************
C/* deck qmmmb2 */
      SUBROUTINE QMMMB2(EXPVAL,DENMAT,WORK,LWORK,IPRINT)
C
c This routine adds the /MM contribution 
c to the 2nd-order magnetic field perturbation due to the use of
c London orbitals. Works in combination with new QMMM code.
c Lund, Nov. 2008, J. Kongsted, based on QM3B3 by K. Ruud.
c
c EXPVAL (OUTPUT): contribution to be added to diamagnetic magnetizability
c
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <orgcom.h>
#include <nuclei.h>
#include <qm3.h>
#include <qmmm.h>
#include <inforb.h>
C
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      LOGICAL LOCDEB, LCOMB
      DIMENSION WORK(*), EXPVAL(6), DENMAT(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C
      LOCDEB = .FALSE.

C     If we have BOTH permanent dipoles AND induced dipoles then these
C     can be added together and we can treat the effect of those
C     simultanously. Then no seperate loop for dipoles and
C     polarizabilities but only one combined loop.

      LCOMB = .FALSE.
      IF ( (IPOLTP .GT. 0) .AND. (NMULT .GE. 1) ) LCOMB = .TRUE.

C     SPLNMR overrules these settings   
      IF (SPLNMR) LCOMB = .FALSE.

C     We don't have these integrals yet !     
      IF (NMULT .GE. 2) THEN
        CALL QUIT('Quadrupole contr. not implemented in QMMMB2')
      ENDIF

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL DZERO(EXPVAL,6)

C     1. Construct the charge contribution for all MM centers

      IF (NMULT .GE. 0) THEN

        KTMP  = 1
        KLAST = KTMP  + 6
        LWRK  = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),6)

        DO 100 I = 1,MMCENT

          IF (ABS(MUL0MM(I)) .LE. THRMM) GOTO 100

          CALL DZERO(WORK(KTMP),6)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)

          CALL GET1IN(DUMMY,'PCMB2SL',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A,I3,A)')
     *       ' N(',L,') QMMM correction matrix '
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
          END IF

          FAC1 = -1.0D0*MUL0MM(I)
          CALL DAXPY(6,FAC1,WORK(KTMP),1,EXPVAL,1)

  100   CONTINUE

      ENDIF !( NMULT .GE. 0)

C     2. Construct the permanent dipole contribution for all MM centers

      IF ( (NMULT .GE. 1) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from permanent dipoles'
        ENDIF

        KTMP  = 1
        KLAST = KTMP  + 18
        LWRK  = LWORK - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),18)

        DO 200 I = 1,MMCENT

          DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2
          DNORM = SQRT(DNORM2)

          IF (ABS(DNORM) .LE. THRMM) GOTO 200

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*MUL1MM(1,I)
          FACY = -1.0D0*MUL1MM(2,I)
          FACZ = -1.0D0*MUL1MM(3,I)

          DO I2 = 1,6
             EXPVAL(I2) = EXPVAL(I2) + WORK(KTMP+I2-1)*FACX +
     &                    WORK(KTMP+I2+5)*FACY + WORK(KTMP+I2+11)*FACZ
          END DO

  200   CONTINUE

      ENDIF !( NMULT .GE. 1)

C     3. Construct the induced dipole contribution for all MM centers
C        NOTE: This could be done more effectively by adding the
C        permanent dipoles and induced dipoles and then performing the
C        London orbital contribution. See the LCOMB option below.

      IF ( (IPOLTP .GT. 0) .AND. (.NOT. LCOMB) ) THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from induced dipoles'
        ENDIF

        KTMP    = 1
        KINDMOM = KTMP    + 18        !INDMOM contains the ind. dips.
        KLAST   = KINDMOM + 3*NNZAL
        LWRK    = LWORK   - KLAST + 1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),18)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMFB2'
           WRITE(LUPRI,*)
        ENDIF

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

        IL = 0 ! Counter for induced dipole moments

        DO 300 I = 1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 300

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KINDMOM+IL+0)
          FACY = -1.0D0*WORK(KINDMOM+IL+1)
          FACZ = -1.0D0*WORK(KINDMOM+IL+2)

          DO I2 = 1,6
             EXPVAL(I2) = EXPVAL(I2) + WORK(KTMP+I2-1)*FACX +
     &                    WORK(KTMP+I2+5)*FACY + WORK(KTMP+I2+11)*FACZ
          END DO

          IL = IL+3

  300   CONTINUE

      ENDIF !(IPOLTP .GT. 0)

C     3a. If both permanent dipoles AND polarizabilities construct TOTAL 
C         dipole contribution for all MM centers

      IF (LCOMB)  THEN

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'London gradient from total dipoles'
        ENDIF

        KTMP    = 1
        KINDMOM = KTMP    + 18       ! INDMOM contains the ind. dips.
        KTOTMOM = KINDMOM + 3*NNZAL  ! TOTMOM contains the TOTAL dips.
        KLAST   = KTOTMOM + 3*MMCENT
        LWRK    = LWORK   - KLAST +1

        IF (LWRK .LT. 0) CALL ERRWRK('QMMMB2',-KLAST,LWORK)

        CALL DZERO(WORK(KTMP),18)
        CALL DZERO(WORK(KINDMOM),3*NNZAL)
        CALL DZERO(WORK(KTOTMOM),3*MMCENT)

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in QMMMB2'
           WRITE(LUPRI,*)
        ENDIF

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDMOM))

        IL = 0
        IJ = 0
        DO 400 I = 1,MMCENT
          IF (ZEROAL(I) .EQ. -1) THEN
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I)
          ELSE IF (ZEROAL(I) .EQ. 1) THEN
            WORK(KTOTMOM+IL+0) = MUL1MM(1,I) + WORK(KINDMOM+IJ+0)
            WORK(KTOTMOM+IL+1) = MUL1MM(2,I) + WORK(KINDMOM+IJ+1)
            WORK(KTOTMOM+IL+2) = MUL1MM(3,I) + WORK(KINDMOM+IJ+2)
            IJ=IJ+3
          ENDIF
          IL=IL+3
  400   CONTINUE

        DO 500 I = 1,MMCENT
          IL = 3*(I - 1)
C         See if the total dipole moment at this site is zero
          DNORM2 = WORK(KTOTMOM+IL+0)*WORK(KTOTMOM+IL+0)
     &           + WORK(KTOTMOM+IL+1)*WORK(KTOTMOM+IL+1)
     &           + WORK(KTOTMOM+IL+2)*WORK(KTOTMOM+IL+2)
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) GOTO 500

          CALL DZERO(WORK(KTMP),18)

          KPATOM = 0
          NOSIM = 0
          TOFILE = .FALSE.
          TRIMAT = .FALSE.
          EXP1VL = .TRUE.
          DIPORG(1) = MMCORD(1,I)
          DIPORG(2) = MMCORD(2,I)
          DIPORG(3) = MMCORD(3,I)
          CALL GET1IN(DUMMY,'EFIELB2',NOSIM,WORK(KLAST),
     &                LWRK,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,WORK(KTMP),EXP1VL,DENMAT,
     &                IPRINT)

          IF (IPRINT .GE. 4) THEN
             WRITE (LUPRI,'(/A)') ' Efield Ex 2.order matrix:'
             CALL OUTPAK(WORK(KTMP),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ey 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+6),3,1,LUPRI)
             WRITE (LUPRI,'(/A)') ' Efield Ez 2.order matrix:'
             CALL OUTPAK(WORK(KTMP+12),3,1,LUPRI)
          END IF

          FACX = -1.0D0*WORK(KTOTMOM+IL+0)
          FACY = -1.0D0*WORK(KTOTMOM+IL+1)
          FACZ = -1.0D0*WORK(KTOTMOM+IL+2)

          DO I2 = 1,6
             EXPVAL(I2) = EXPVAL(I2) + WORK(KTMP+I2-1)*FACX +
     &                    WORK(KTMP+I2+5)*FACY + WORK(KTMP+I2+11)*FACZ
          END DO

  500   CONTINUE

      ENDIF !(LCOMB)

C     Put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ
C
      RETURN
      END
C****************************************************************************

