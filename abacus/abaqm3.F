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

                  IF (IQM3PR .GE. 4) THEN
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
      LUCCRA = -1
      LUCCEA = -1
      CALL CC_GET31(LUCCRA,'CC_RA',NCOMS,WORK(KRAX),
     &              WORK(KRAY),WORK(KRAZ))
      CALL CC_GET31(LUCCEA,'ENSAFILE',NCOMS,WORK(KENSAX),
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
      LUCCRA = -1
      LUCCEA = -1
      CALL CC_GET31(LUCCRA,'CC_RA',NCOMS,WORK(KRAX),
     &              WORK(KRAY),WORK(KRAZ))
      CALL CC_GET31(LUCCEA,'ENSAFILE',NCOMS,WORK(KENSAX),
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

