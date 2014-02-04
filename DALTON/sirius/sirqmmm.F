!
!...   Copyright (c) 2013 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2013 (2013), see http://daltonprogram.org"
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
C FILE: sirius/sirqmmm.F
C Oct. 2009: JMO and AHS
C   Moved all routines relevant for the new QMMM code to this file and
C   added parallel QMMM routines.
*******************************************************************************
C  /* Deck qmmmfck */
      SUBROUTINE QMMMFCK(DCAO,DVAO,FSOL,EQMMM,WRK,LWRK,IPRINT)

#include "implicit.h"
#include "dummy.h"
#include "inftap.h"
#include "priunit.h"
#include "mxcent.h"
#include "qmmm.h"
#include "thrzer.h"
#include "iratdef.h"
#include "codata.h"
#include "maxash.h"
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"

      DIMENSION DCAO(*), DVAO(*)
      DIMENSION FSOL(*), WRK(LWRK)

      CALL QENTER('QMMMFCK')
      KDTAO = 1
      KTAO = KDTAO + NNBASX
      KEND = KTAO + NNBASX
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMMFCK',-KEND,LWRK)


C     Get total density
      IF (NASHT .EQ. 0) THEN
            CALL PKSYM1(WRK(KDTAO),DCAO,NBAS,NSYM,-1)
      ELSE
            DO I = 1,NNBAST
               IF (HSFOCK) THEN
                  WRK(KTAO-1+I) = DCAO(I)
               ELSE
                  WRK(KTAO-1+I) = DCAO(I) + DVAO(I)
               END IF
            END DO
            CALL PKSYM1(WRK(KDTAO),WRK(KTAO),NBAS,NSYM,-1)
      END IF

C     Modify the fock operator. Modification returned in FSOL. 
C     QMMM contribution to the energy returned in EQMMM.
      CALL QMMM_FCK_AO(FSOL,WRK(KDTAO),EQMMM,WRK(KEND),LWRK1,IPRINT)

      CALL QEXIT('QMMMFCK')
      RETURN
      END

C******************************************************************************
C  /* Deck qmmm_fck_ao */
      SUBROUTINE QMMM_FCK_AO(FSOL,DCAO,ESOLT,WRK,LWRK,IPRINT)   

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "mmtimes.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "infpar.h"

      DIMENSION WRK(LWRK)
      DIMENSION FSOL(*),DCAO(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      LOGICAL EXCENT,LOCDEB
      INTEGER NZERAL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      LOGICAL FIRST
      SAVE    FIRST
      DATA    FIRST /.TRUE./

      CALL QENTER('QMMM_FCK_AO')

      LOCDEB = .FALSE.

      KTAO   = 1
      KWRK1  = KTAO   + NNBASX
      LWRK1  = LWRK   - KWRK1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMM_FCK_AO',-KWRK1,LWRK)

      CALL DZERO(WRK(KTAO),NNBASX)

C     The different static energy contributions
      ECHART = 0.0D0
      EDIPT  = 0.0D0
      EQUADT = 0.0D0

C     Backup diporg. We use diporg to transfer coordinates to int. program.

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      IF (MMTIME) DTIME = SECOND()
#if defined(VAR_MPI)
      IF (NODTOT.GE.1) THEN
C     All the corrections to the Fock/KS operator due to the static
C     multipoles when the calculation is done in parallel
        CALL PARQMMM_M(DCAO,WRK(KTAO),ESOLT,LOCDEB,
     &                  WRK(KWRK1),LWRK1,IPRINT)
      ELSE
#endif
C     1) The charge correction to the Fock/KS operator
        IF (NMULT .GE. 0) CALL QMMM_CHARGE(DCAO,ESOLT,WRK(KTAO),
     &                       LOCDEB,FIRST,WRK(KWRK1),LWRK1,IPRINT)
C     2) The dipole correction to the Fock/KS operator
        IF (NMULT .GE. 1) CALL QMMM_DIPOLE(DCAO,ESOLT,WRK(KTAO),
     &                       LOCDEB,FIRST,WRK(KWRK1),LWRK1,IPRINT)

C     3) The quadrupole correction to the Fock/KS operator
        IF (NMULT .GE. 2) CALL QMMM_QUADPOLE(DCAO,ESOLT,WRK(KTAO),
     &                       LOCDEB,FIRST,WRK(KWRK1),LWRK1,IPRINT)


#if defined(VAR_MPI)
      ENDIF ! nodtot .ge. 1
#endif
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMMULPOL = TMMMULPOL + DTIME
      ENDIF

      IF ( (IPRINT.GT.1) .OR. (LOCDEB) ) THEN
        write(lupri,*)
        write(lupri,*) 'MM-charge QM density interaction energy:',ECHART
        write(lupri,*) 'MM-dipole QM density interaction energy:',EDIPT
        write(lupri,*) 'MM-quadr. QM density interaction energy:',EQUADT
      ENDIF 

C     5) The polarization correction to the Fock/KS operator
      
      IF (MMTIME) DTIME = SECOND()
      IF (IPOLTP .GT. 0) CALL QMMM_POLARI(DCAO,ESOLT,WRK(KTAO),
     &                       LOCDEB,WRK(KWRK1),LWRK1,IPRINT)
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMPOL = TMMPOL + DTIME
      ENDIF

C     Finally, put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL PKSYM1(WRK(KTAO),FSOL,NBAS,NSYM,+1)
      CALL QEXIT('QMMM_FCK_AO')

      IF (FIRST) THEN
         FIRST = .FALSE.
      END IF

      RETURN
      END
C******************************************************************************
C  /* Deck qmmm_charge */
      SUBROUTINE QMMM_CHARGE(DCAO,ESOLT,TAO,LOCDEB,FIRST,
     &                       WRK,LWRK,IPRINT)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "infpar.h"

      DIMENSION WRK(LWRK), TAO(NNBASX), DCAO(*)
      LOGICAL LOCDEB, FIRST

      CALL QENTER('QMMM_CHARGE')

      KTAO   = 1
      KNSEL  = KTAO   + NNBASX
      KNSNUC = KNSEL  + MMCENT
      KLAST  = KNSNUC + MMCENT
      LWRK2  = LWRK   - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_CHARGE 1',-KLAST,LWRK)

      FAC1   = 1.0D0      
      EXPNST = 0.0D0
      ECHCH  = 0.0D0

      CALL DZERO(WRK(KTAO),NNBASX)

      DO 100 I = 1,MMCENT 

         DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *           (MMCORD(2,I)-QMCOM(2))**2 +
     *           (MMCORD(3,I)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
           WRK(KNSEL + I - 1)  = 0.0D0
           WRK(KNSNUC + I - 1) = 0.0D0
           IF (LOCDEB) THEN
             WRITE(LUPRI,*) 'Skipping charge ', I
           ENDIF 
           GOTO 100
         ENDIF

         CALL CHARGE_ITER(I,DCAO,WRK(KNSEL+I-1),WRK(KNSNUC+I-1),LOCDEB,
     *                    WRK(KTAO),WRK(KLAST),LWRK2,IPRINT)
         EXPNST = EXPNST + WRK(KNSEL+I-1)
         ECHCH  = ECHCH  + WRK(KNSNUC+I-1)

 100  CONTINUE
C Transfering the QM nuclei - MM multipole energy contribution
C to the CC part of the code. We start with the charge contribution
      ENUMUL = 0.0D0
      ENUMUL = ECHCH

      IF (FIRST) THEN
C     Write integrals to file        
         LUQMMM = -1
         IF (LUQMMM .LT. 0) THEN
            CALL GPOPEN(LUQMMM,'MU0INT','UNKNOWN','','',
     &                  IDUMMY,.FALSE.)
         ENDIF
         REWIND(LUQMMM)
         WRITE (LUQMMM) (WRK(KTAO+J-1), J=1,NNBASX)
         CALL GPCLOSE(LUQMMM,'KEEP')
      ENDIF

      CALL DAXPY(NNBASX,1.0D0,WRK(KTAO),1,TAO,1)

      ECHART = EXPNST + ECHCH
      ESOLT  = ECHART

      IF ( (IPRINT.GT.5) .OR. (LOCDEB) ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Charge-electronic Charge-nuclear Total'
        DO 102 I = 1,MMCENT
          ELTEMP = WRK(KNSEL + I - 1) + WRK(KNSNUC + I - 1)
          WRITE(LUPRI,*) I,WRK(KNSEL + I - 1),WRK(KNSNUC + I - 1),ELTEMP
  102   CONTINUE

        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EXPNST, ECHCH, EXPNST+ECHCH
        WRITE(LUPRI,*)
      ENDIF

      CALL QEXIT('QMMM_CHARGE')

      RETURN
      END
C******************************************************************************
C  /* Deck qmmm_dipole */
      SUBROUTINE QMMM_DIPOLE(DCAO,ESOLT,TAO,LOCDEB,FIRST,
     &                       WRK,LWRK,IPRINT)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "infpar.h"

      DIMENSION WRK(LWRK), TAO(NNBASX), DCAO(*)
      LOGICAL LOCDEB, FIRST

      CALL QENTER('QMMM_DIPOLE')

      KTAO   = 1
      KNSEL  = KTAO + NNBASX
      KNSNUC = KNSEL + MMCENT
      KLAST  = KNSNUC + MMCENT
      LWRK2  = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_DIPOLE 1',-KLAST,LWRK)

      FAC1   = 1.0D0      
      FACM1  = -1.0D0
      EMUL1T = 0.0D0
      ELOCT  = 0.0D0

      CALL DZERO(WRK(KTAO),NNBASX)

      DO 200 I = 1,MMCENT 

         DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
            WRK(KNSEL + I - 1)  = 0.0D0
            WRK(KNSNUC + I - 1) = 0.0D0
            IF (LOCDEB) THEN
               WRITE(LUPRI,*) 'Skipping dipole ', I
            ENDIF
            GOTO 200
         ENDIF

         CALL DIPOLE_ITER(I,DCAO,WRK(KNSEL+I-1),WRK(KNSNUC+I-1),LOCDEB,
     *                   WRK(KTAO),WRK(KLAST),LWRK2,IPRINT)
         EMUL1T = EMUL1T + WRK(KNSEL+I-1)
         ELOCT  = ELOCT  + WRK(KNSNUC+I-1)

 200  CONTINUE

C Add up QM nuclei - multipole energy contributions to be used in CC
      ENUMUL = ENUMUL + ELOCT
      IF (FIRST) THEN
C     Write integrals to file
         LUQMMM = -1
         IF (LUQMMM .LT. 0) THEN
            CALL GPOPEN(LUQMMM,'MU1INT','UNKNOWN','','',
     &                  IDUMMY,.FALSE.)
         ENDIF
         REWIND(LUQMMM)
         WRITE (LUQMMM) (WRK(KTAO+J-1), J=1,NNBASX)
         CALL GPCLOSE(LUQMMM,'KEEP')
      ENDIF

      CALL DAXPY(NNBASX,1.0D0,WRK(KTAO),1,TAO,1)

      EDIPT = EMUL1T + ELOCT
      ESOLT = ESOLT + EDIPT

      IF ( (IPRINT.GT.5) .OR. (LOCDEB) ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Dipole-electronic Dipole-nuclear Total'
        DO 202 I = 1,MMCENT
          ETEMP = WRK(KNSEL + I - 1) + WRK(KNSNUC + I - 1)
          WRITE(LUPRI,*) I,WRK(KNSEL + I - 1),WRK(KNSNUC + I - 1),ETEMP
  202   CONTINUE
  
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EMUL1T, ELOCT, EMUL1T+ELOCT
        WRITE(LUPRI,*) 
      ENDIF

      CALL QEXIT('QMMM_DIPOLE')

      RETURN
      END
C******************************************************************************
C  /* Deck qmmm_quadpole */
      SUBROUTINE QMMM_QUADPOLE(DCAO,ESOLT,TAO,LOCDEB,FIRST,
     &                       WRK,LWRK,IPRINT)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "infpar.h"

      DIMENSION WRK(LWRK), TAO(NNBASX), DCAO(*)
      LOGICAL LOCDEB, FIRST

      CALL QENTER('QMMM_QUADPOLE')

      KTAO   = 1
      KNSEL  = KTAO   + NNBASX
      KNSNUC = KNSEL  + MMCENT
      KLAST  = KNSNUC + MMCENT
      LWRK2  = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_QUADPOLE 1',-KLAST,LWRK)

      FAC1   = 1.0D0
      EMUL2T = 0.0D0
      ELOCT  = 0.0D0

      CALL DZERO(WRK(KTAO),NNBASX)

      DO 300 I = 1,MMCENT

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
           WRK(KNSEL + I - 1)  = 0.0D0
           WRK(KNSNUC + I - 1) = 0.0D0
           IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping quadrupole ', I
           ENDIF
           GOTO 300
        ENDIF

        CALL QUADPOLE_ITER(I,DCAO,WRK(KNSEL+I-1),WRK(KNSNUC+I-1),LOCDEB,
     *                    WRK(KTAO),WRK(KLAST),LWRK2,IPRINT)
        EMUL2T = EMUL2T + WRK(KNSEL+I-1)
        ELOCT  = ELOCT  + WRK(KNSNUC+I-1)

 300  CONTINUE

C Add up QM nuclei - multipole energy contributions to be used in CC
      ENUMUL = ENUMUL + ELOCT
      IF (FIRST) THEN
C     Write integrals to file
         LUQMMM = -1
         IF (LUQMMM .LT. 0) THEN
            CALL GPOPEN(LUQMMM,'MU2INT','UNKNOWN','','',
     &                  IDUMMY,.FALSE.)
         ENDIF
         REWIND(LUQMMM)
         WRITE (LUQMMM) (WRK(KTAO+J-1), J=1,NNBASX)
         CALL GPCLOSE(LUQMMM,'KEEP')
      ENDIF

      CALL DAXPY(NNBASX,1.0D0,WRK(KTAO),1,TAO,1)

      EQUADT = EMUL2T + ELOCT
      ESOLT  = ESOLT  + EQUADT

      IF ( (IPRINT.GT.5) .OR. (LOCDEB) ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Quadr.-electronic Quadr.-nuclear Total'
        DO 302 I = 1,MMCENT
          ETEMP = WRK(KNSEL + I - 1) + WRK(KNSNUC + I - 1)
          WRITE(LUPRI,*) I,WRK(KNSEL + I - 1),WRK(KNSNUC + I - 1),ETEMP
  302   CONTINUE
  
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EMUL2T, ELOCT, EMUL2T+ELOCT
        WRITE(LUPRI,*) 
      ENDIF

      CALL QEXIT('QMMM_QUADPOLE')

      RETURN
      END

C******************************************************************************
C  /* Deck qmmm_polari */
      SUBROUTINE QMMM_POLARI(DCAO,ESOLT,TAO,LOCDEB,
     &                       WRK,LWRK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "qmmm.h"
#include "mmtimes.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "infpar.h"

      DIMENSION WRK(LWRK), TAO(NNBASX), DCAO(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, EXCENT, LOCDEB, LSKIP
      INTEGER NZERAL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_POLARI')

C     Zero out a list of centers having zero polarizability. We don't
C     know yet the number of zero polarizabilities so we take the
C     worst case, i.e. MXMMCT, for the length of this list

      DO 443 I=1,MXMMCT
        ZEROAL(I) = 0
 443  CONTINUE
c
c     Check if the polarizability is equal to zero; if so put -1 on
c     the list for this center. If not equal to zero put +1 on the
c     list for this center and if not touched upon leave zero

      LIZA = 1   ! Counts centers having polarizability equal to zero

      DO 400 I=1,MMCENT

        IF (IPOLTP .EQ. 1) THEN
          ANORM2 = 3*(POLIMM(I)**2)
          ANORM  = SQRT(ANORM2)
          IF (ANORM .LE. THRMM) THEN
            ZEROAL(I) = -1
            LIZA = LIZA + 1
          ELSE 
            ZEROAL(I) = 1
          ENDIF
        ENDIF

        IF (IPOLTP .EQ. 2) THEN
          ANORM2 = POLMM(1,I)**2 + 2*(POLMM(2,I)**2) +
     &                             2*(POLMM(3,I)**2) +
     &             POLMM(4,I)**2 + 2*(POLMM(5,I)**2) +
     &             POLMM(6,I)**2
          ANORM  = SQRT(ANORM2)
          IF (ANORM .LE. THRMM) THEN
            ZEROAL(I) = -1
            LIZA = LIZA + 1
          ELSE
            ZEROAL(I) = 1
          ENDIF
        ENDIF

 400  CONTINUE

      NZERAL = LIZA - 1
      NNZAL  = MMCENT - NZERAL  ! Number of MM centers with ALPHA .NE. 0

      IF ( (IPRINT.GT.1) .OR. (LOCDEB) ) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Number of polarizable sites: ', NNZAL
        WRITE(LUPRI,*)
      ENDIF

      IF (MMMAT) THEN

        KINVMAT   = 1
        KINDMOM   = KINVMAT + 3*NNZAL*(3*NNZAL+1)/2 ! for packed response matrix
        KMAT      = KINDMOM + 3*NNZAL ! List for induced dipoles
        KIPVT     = KMAT    + 3*NNBASX ! For Rr_a integrals
        KWRKV     = KIPVT   + 3*NNZAL ! For matrix inv.
        KTAO      = KWRKV   + 3*NNZAL ! For matrix inv.
        KWRK2     = KTAO    + NNBASX
        LWRK2     = LWRK    - KWRK2 + 1

        IF (LWRK2 .LT. 0) THEN
          CALL ERRWRK('QMMM_POLARI 1',-KWRK2,LWRK)
        ENDIF

        CALL DZERO(WRK(KINVMAT), 3*NNZAL*(3*NNZAL+1)/2)
        CALL DZERO(WRK(KINDMOM), 3*NNZAL)
        CALL DZERO(WRK(KIPVT), 3*NNZAL)
        CALL DZERO(WRK(KWRKV), 3*NNZAL)
        CALL DZERO(WRK(KMAT), 3*NNBASX)

C       FIXDIP assumes induced dipoles are calculated in a previous run.
C       Mainly due to debugging. Assumes identical molecules and order
C       of atoms in previous and current run.

        IF (.NOT. FIXDIP) THEN
          CALL GET_IND_DIPOLES_1(DCAO,NNZAL,WRK(KINVMAT),WRK(KINDMOM),
     &                           WRK(KWRK2),WRK(KIPVT),WRK(KWRKV),
     &                           LWRK2,IPRINT)
        ELSE
          WRITE(LUPRI,*) 'Ind. dips. from a prev. calc. read from file'
          CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WRK(KINDMOM))
        ENDIF
      ELSE IF (MMITER) THEN

        KINDMOM = 1                 
        KMAT    = KINDMOM + 3*NNZAL      ! List for induced dipoles
        KTAO    = KMAT    + 3*NNBASX
        KWRK2   = KTAO    + NNBASX
        LWRK2   = LWRK    - KWRK2 + 1

        IF (LWRK2 .LT. 0) THEN
          CALL ERRWRK('QMMM_POLARI 2',-KWRK2,LWRK)
        ENDIF

        CALL DZERO(WRK(KINDMOM),(3*NNZAL))

C       FIXDIP assumes induced dipoles are calculated in a previous run.
C       Mainly due to debugging. Assumes identical molecules and order
C       of atoms in previous and current run.

        IF (.NOT. FIXDIP) THEN
          IF (MMTIME) DTIME = SECOND()
          CALL GET_IND_DIPOLES_2(DCAO,NNZAL,WRK(KINDMOM),
     &                           WRK(KWRK2),LWRK2,IPRINT)
          IF (MMTIME) THEN
            DTIME = SECOND() - DTIME
            TMMGID2 = TMMGID2 + DTIME
          ENDIF
        ELSE
          WRITE(LUPRI,*) 'Ind. dips. from a prev. calc. read from file'
          CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WRK(KINDMOM))
        ENDIF

      ENDIF

C     Compute polarization contributions to the Fock/KS matrix and
C     total solvation energy

      FACM1 = -1.0D0
      IINIM = 0   ! important should be zero due to the indexing used !

      EDELD  = 0.0D0            ! For interaction with electronic density
      EDNUC  = 0.0D0            ! For interaction with QM nuclei
      ED0MOM = 0.0D0            ! For interaction with point-charges
      ED1MOM = 0.0D0            ! For interaction with permanent dipoles
      ED2MOM = 0.0D0            ! For interaction with quadrupoles
      EDMULT = 0.0D0            ! For interaction with permanent multipoles
      EPOLT  = 0.0D0            ! Total polarization energy

      CALL DZERO(WRK(KTAO),NNBASX)

#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
        CALL MM_POLAR_CONTR_M(DCAO(1),WRK(KTAO),WRK(KINDMOM),
     &                      WRK(KWRK2),LWRK2,IPRINT)
      ELSE
#endif
        KEDALL = KWRK2
        KWRK3  = KEDALL + 6
        LWRK3  = LWRK - KWRK3 + 1

        DO 500 I=1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 500

          CALL GET_POL_CONTR(I,WRK(KINDMOM+IINIM),WRK(KEDALL),
     &                        DCAO,WRK(KTAO),WRK(KWRK3),LWRK3)

          EDELD  = EDELD  + WRK(KEDALL)
          EDNUC  = EDNUC  + WRK(KEDALL+1)
          ED0MOM = ED0MOM + WRK(KEDALL+2)
          ED1MOM = ED1MOM + WRK(KEDALL+3)
          ED2MOM = ED2MOM + WRK(KEDALL+4)

          IINIM = IINIM + 3

 500    CONTINUE

        EDMULT = ED0MOM + ED1MOM + ED2MOM

#if defined(VAR_MPI)
      ENDIF                     ! IF (NODTOT .GE. 1) ... ELSE
#endif
      CALL DAXPY(NNBASX,1.0D0,WRK(KTAO),1,TAO,1)

      EPOLT  = EDELD + EDNUC + EDMULT

      ESOLT = ESOLT + EPOLT

      IF (IPRINT .GT. 1) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,5001) 
        WRITE(LUPRI,*)
        WRITE(LUPRI,5002) EDELD
        WRITE(LUPRI,5003) EDNUC
        WRITE(LUPRI,5004) EDMULT
        WRITE(LUPRI,*)
        WRITE(LUPRI,5005) EPOLT
        WRITE(LUPRI,*)
      ENDIF

C      IF (MMPROP) CALL MM_PROPS(WRK(KWRK2),LWRK2,IPRINT)

 5001 FORMAT(' Polarization energy: ')
 5002 FORMAT('      Electronic contribution:   ',F15.9)
 5003 FORMAT('      Nuclear contribution:      ',F15.9)
 5004 FORMAT('      Multipole contribution:    ',F15.9)
 5005 FORMAT('      Total:                     ',F15.9)



      CALL QEXIT('QMMM_POLARI')

      RETURN
      END
C
C******************************************************************************
C  /* Deck qmmmfckmo */
      SUBROUTINE QMMMFCKMO(CMO,FSOL,WRK,LWRK,IPRINT)
C
C     Construct the QMMM contribution to the Fock-matrix in MO basis
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "inforb.h"
#include "infopt.h"
C
      DIMENSION CMO(*), FSOL(*), WRK(LWRK)
C
      CALL QENTER('QMMMFCKMO')
C
      KDV     = 1
      KDENS   = KDV     + N2BASX
      KDVS    = KDENS   + NNBASX
      KFSOLAO = KDVS    + NNBASX
      KUCMO   = KFSOLAO + NNBASX
      KZERO   = KUCMO   + NORBT*NBAST
      KWRK    = KZERO   + NNBASX
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMMFCKMO',-KWRK,LWRK)

      CALL DZERO(WRK(KZERO),NNBASX)

C     Construct the density matrix
      CALL FCKDEN((NISHT.GT.0),.FALSE.,WRK(KDV),
     *            DUMMY,CMO,DUMMY,WRK(KWRK),LWRK1)

      CALL DGEFSP(NBAST,WRK(KDV),WRK(KDVS))
      CALL PKSYM1(WRK(KDVS),WRK(KDENS),NBAS,NSYM,1)

C     Construct the QMMM contribution to the Fock-matrix in AO
C     For the openshell density we Put in zero as this is now included
C     in 
C     KDENS already.
      CALL QMMMFCK(WRK(KDENS),WRK(KZERO),WRK(KFSOLAO),ESOLT,
     *             WRK(KWRK),LWRK1,IPRINT)

C     Transform to mo
      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL UTHU(WRK(KFSOLAO),FSOL,WRK(KUCMO),WRK(KWRK),
     &             NBAST,NORBT)
C
      CALL QEXIT('QMMMFCKMO')
      RETURN
      END
C******************************************************************************
C  /* Deck GET_IND_DIPOLES_1 */
      SUBROUTINE GET_IND_DIPOLES_1(DCAO,POLDIM,INVMAT,INDMOM,WRK,IPVT,
     &                             WRKV,LWRK,IPRINT)
C
C A subroutine to calculate induced dipole moments
C
C Input:
C
C   DCAO    - density matrix in AO basis
C   POLDIM  - number of polarizable MM centers. Actually in common as
C             NNZAL.
C
C Output:
C
C   INVMAT  - the classical response matrix, i.e. [ALPHA^(-1) - T]^(-1)
C   INDMOM  - a vector containing induced dipole moments
C
C From Common   
C
C   ZEROAL  - a vector containing +1 for polarizable MM centers and -1
C             for non-polarizable
C
C Oct. 2009: JMO 
C   Changed the routines used to construct the classical response matrix
C   to more efficient ones that use the fact that it is symmetric.
C Sep 2010: JMO & KS
C   Starting sharing of DFT/MM and CC/MM field routines
C Oct 2010: AHS 
C   Sharing of parallel and serial routines
C Jan 2011: JMO
C   Construct the classical response matrix using packed storage.
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "infpar.h"

      CHARACTER LLAB
      LOGICAL EXCENT,FNDLAB, LSKIP
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB
      INTEGER POLDIM, IPVT
      DOUBLE PRECISION WRK, INVMAT, INDMOM, WRKV
      DIMENSION INVMAT(3*POLDIM*(3*POLDIM+1)/2)
      DIMENSION INDMOM(3*POLDIM)
      DIMENSION IPVT(3*POLDIM)
      DIMENSION WRKV(3*POLDIM)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)
      DIMENSION DCAO(*)

      CHARACTER*8 LABINT(9*MXCENT)

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('GET_IND_DIPOLES_1')

      LOCDEB = .FALSE.

      IF (POLDIM .NE. NNZAL) THEN
        WRITE(LUPRI,*) 'ERROR in no. of polarizabilities'
        CALL QUIT('ERROR in GET_IND_DIPOLES_1')
      ENDIF

C     Allocate memory for electric field integrals and electric fields
C     (the order KELF KELFEL KELFNU has to be kept because of QMMM_POLARI_M1! AHS)
      KMAT    = 1
      KELF    = KMAT + 3*NNBASX      ! For electric field integrals
      KELFEL  = KELF + 3*POLDIM      ! For total OR (if SPLDIP) multipole electric field
      IF (SPLDIP) THEN
        KELFNU = KELFEL + 3*POLDIM   ! For electronic electric field
        KIMMUL = KELFNU + 3*POLDIM   ! For nuclear electric field
        KIMNUC = KIMMUL + 3*POLDIM   ! For induced moments due to permanent multipoles
        KIMELD = KIMNUC + 3*POLDIM   ! For induced moments due to QM nuclei
        KEND   = KIMELD + 3*POLDIM   ! For induced moments due to electronic density
      ELSE 
        KEND   = KELFEL
      ENDIF
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('GET_IND_DIPOLES_1',-KEND,LWRK)

      CALL DZERO(WRK(KMAT),3*NNBASX)
      CALL DZERO(WRK(KELF),3*POLDIM)
      IF (SPLDIP) THEN
        CALL DZERO(WRK(KELFEL),3*POLDIM)
        CALL DZERO(WRK(KELFNU),3*POLDIM)
        CALL DZERO(WRK(KIMMUL),3*POLDIM)
        CALL DZERO(WRK(KIMNUC),3*POLDIM)
        CALL DZERO(WRK(KIMELD),3*POLDIM)
      ENDIF

C     Form F vector due to permanent MM moments
#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
        CALL MM_FIELD_M1(DCAO(1),WRK(KELF),POLDIM,
     &                      WRK(KEND),LWRK1,IPRINT)
      ELSE
#endif
        LRI = 1                 ! Row index in the large matrix

        DO 200 I=1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 200

          CALL GET_FIELD(I,LRI,WRK(KELF),WRK(KELFEL),WRK(KELFNU),
     &                          DCAO,LOCDEB,WRK(KEND),LWRK1)
          LRI = LRI + 3
 200    CONTINUE

#if defined(VAR_MPI)
      ENDIF
#endif
      NDIM = 3*POLDIM

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'Done generating the F-Vector'
        WRITE(LUPRI,*) 'Done generating the interaction matrix'
        WRITE(LUPRI,*) 'F-Vector'
        DO 777 KK=1,NDIM
          WRITE(LUPRI,*) WRK(KELF+KK-1)
 777    CONTINUE
      ENDIF

C     If needed, construct the [ALPHA^(-1) - T]^(-1) matrix and write it to 
C     file. ELSE: read matrix from the file. CONMAT = CONstruct MATrix

      IF (CONMAT) THEN

        CALL MAKE_QMMM_INVERSE_RESPONSE_MATRIX(INVMAT,POLDIM) ! Construct inverse response matrix

        IF (IPRINT .GT. 1) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) ' The classical response matrix is'//
     &                   ' explicitly constructed. '
          WRITE(LUPRI,*) ' Dimension is: ',NDIM
          WRITE(LUPRI,*)
        ENDIF

        IF ((IPRINT.GT.15) .OR. (LOCDEB)) THEN
          WRITE(LUPRI,*)'Matrix to be inverted: '
          DO I = 1, NDIM*(NDIM+1)/2
            WRITE(LUPRI,*) INVMAT(I)
          END DO
        END IF

        IF (IPRINT.GT.1) CALL TIMER('START ',TIMSTR,TIMEND)

C       Construct the classical response matrix
        CALL DSPTRF('L', NDIM, INVMAT, IPVT, INFO)
        IF (INFO .NE. 0) THEN
          CALL QUIT('ERROR: construction of the classical'//
     &              ' response matrix failed!')
        END IF
        CALL DSPTRI('L', NDIM, INVMAT, IPVT, WRKV, INFO)
        IF (INFO .NE. 0) THEN
          CALL QUIT('ERROR: construction of the classical response'// 
     &              ' matrix failed!')
        END IF
        IF(IPRINT.GT.1) CALL TIMER('MATINV',TIMSTR,TIMEND)

        IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
          WRITE(LUPRI,*)'Classical response matrix: '
          DO I = 1, NDIM*(NDIM+1)/2
            WRITE(LUPRI,*) INVMAT(I)
          END DO
        END IF

C       We write the classical response matrix to file

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     $               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF

        REWIND(LUQMMM)
        CALL WRTIEF(INVMAT,NDIM*(NDIM+1)/2,'QQMMMMAT',LUQMMM)
        CALL GPCLOSE(LUQMMM,'KEEP')

        IF (RELMAT) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'The classical response matrix saved in QMMMIM.'
         WRITE(LUPRI,*)
         STOP
        ENDIF

        CONMAT = .FALSE.

      ELSE  ! read the inverted matrix from the file

        IF (IPRINT .GT. 5) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) ' The classical response matrix is'//
     &                   ' read from the file. '
          WRITE(LUPRI,*)
        ENDIF

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
        REWIND(LUQMMM)

        CALL DZERO(INVMAT, NDIM*(NDIM+1)/2)
        
        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,NDIM*(NDIM+1)/2,INVMAT)
        ELSE
          CALL QUIT('Problem reading the classical response matrix'//
     &              ' from QMMMIM file')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
          WRITE(LUPRI,*) ' The classical response matrix is'//
     &                   ' read from the QMMMIM file: '
          DO I = 1, NDIM*(NDIM+1)/2
            WRITE(LUPRI,*) INVMAT(I)
          END DO
        ENDIF

      ENDIF

      IF (IPRINT .GT. 1) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,1051)
        WRITE(LUPRI,1050)
        WRITE(LUPRI,1051)
        WRITE(LUPRI,*)
      ENDIF

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'F-Vector'
        DO 899 I=1,NDIM
        WRITE(LUPRI,*) WRK(KELF+I-1)
 899    CONTINUE       

      ENDIF

      IF (SPLDIP) THEN
        CALL DSPMV('L', NDIM, D1, INVMAT, WRK(KELF), 1, D0,
     &             WRK(KIMMUL), 1)
        CALL DSPMV('L', NDIM, D1, INVMAT, WRK(KELFNU), 1, D0,
     &             WRK(KIMNUC), 1)
        CALL DSPMV('L', NDIM, D1, INVMAT, WRK(KELFEL), 1, D0,
     &             WRK(KIMELD), 1)
      ELSE
        CALL DSPMV('L', NDIM, D1, INVMAT, WRK(KELF), 1, D0, INDMOM, 1)
      ENDIF

C     Write the nonzero induced dipoles to files. Only if not fixdip.
      IF ( (.NOT. FIXDIP) .AND. (SPLDIP) ) THEN
        CALL PUT_TO_FILE_1('INDUCED_DIPOLES_MUL',POLDIM,WRK(KIMMUL))
        CALL PUT_TO_FILE_1('INDUCED_DIPOLES_NUC',POLDIM,WRK(KIMNUC))
        CALL PUT_TO_FILE_1('INDUCED_DIPOLES_ELE',POLDIM,WRK(KIMELD))
      ENDIF


      IF (SPLDIP) THEN

        DO 420 I=1,NDIM
          INDMOM(I) = WRK(KIMMUL+I-1) + WRK(KIMNUC+I-1) +
     &                WRK(KIMELD+I-1)
 420    CONTINUE

        IIMIEL = 1
        IIMINU = 1
        IIMIMU = 1

        WRITE(LUPRI,*)
        WRITE(LUPRI,1040)
        WRITE(LUPRI,*)
        WRITE(LUPRI,1000)
        WRITE(LUPRI,1010)
        WRITE(LUPRI,1000)
        DO 421 I=1,MMCENT
          IF (ZEROAL(I) .EQ. -1) THEN
            DIPX = 0.0D0
            DIPY = 0.0D0
            DIPZ = 0.0D0
          ELSE
            DIPX = WRK(KIMELD+IIMIEL-1+0)
            DIPY = WRK(KIMELD+IIMIEL-1+1)
            DIPZ = WRK(KIMELD+IIMIEL-1+2)
            IIMIEL = IIMIEL + 3
          ENDIF
          WRITE(LUPRI,1020) I,DIPX,DIPY,DIPZ
 421    CONTINUE
        WRITE(LUPRI,1000)
        WRITE(LUPRI,*)

        WRITE(LUPRI,*)
        WRITE(LUPRI,1041)
        WRITE(LUPRI,*)
        WRITE(LUPRI,1000)
        WRITE(LUPRI,1010)
        WRITE(LUPRI,1000)
        DO 422 I=1,MMCENT
          IF (ZEROAL(I) .EQ. -1) THEN
            DIPX = 0.0D0
            DIPY = 0.0D0
            DIPZ = 0.0D0
          ELSE
            DIPX = WRK(KIMNUC+IIMINU-1+0)
            DIPY = WRK(KIMNUC+IIMINU-1+1)
            DIPZ = WRK(KIMNUC+IIMINU-1+2)
            IIMINU = IIMINU + 3
          ENDIF
          WRITE(LUPRI,1020) I,DIPX,DIPY,DIPZ
 422    CONTINUE
        WRITE(LUPRI,1000)
        WRITE(LUPRI,*)

        WRITE(LUPRI,*)
        WRITE(LUPRI,1042)
        WRITE(LUPRI,*)
        WRITE(LUPRI,1000)
        WRITE(LUPRI,1010)
        WRITE(LUPRI,1000)
        DO 423 I=1,MMCENT
          IF (ZEROAL(I) .EQ. -1) THEN
            DIPX = 0.0D0
            DIPY = 0.0D0
            DIPZ = 0.0D0
          ELSE
            DIPX = WRK(KIMMUL+IIMIMU-1+0)
            DIPY = WRK(KIMMUL+IIMIMU-1+1)
            DIPZ = WRK(KIMMUL+IIMIMU-1+2)
            IIMIMU = IIMIMU + 3
          ENDIF
          WRITE(LUPRI,1020) I,DIPX,DIPY,DIPZ
 423    CONTINUE
        WRITE(LUPRI,1000)
        WRITE(LUPRI,*)

      ENDIF

      IF (IPRINT .GT. 1) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,1030)
        WRITE(LUPRI,*)
        WRITE(LUPRI,1000)
        WRITE(LUPRI,1010)
        WRITE(LUPRI,1000)
      ENDIF 

      IINIM = 1

      DO 500 I=1,MMCENT
        IF (ZEROAL(I) .EQ. -1) THEN
          DIPX = 0.0D0
          DIPY = 0.0D0
          DIPZ = 0.0D0
        ELSE
          DIPX = INDMOM(IINIM+0)
          DIPY = INDMOM(IINIM+1)
          DIPZ = INDMOM(IINIM+2)
          IINIM = IINIM + 3
        ENDIF
        IF (IPRINT .GT. 1) WRITE(LUPRI,1020) I,DIPX,DIPY,DIPZ
 500  CONTINUE

      IF (IPRINT .GT. 1) THEN
        WRITE(LUPRI,1000)
        WRITE(LUPRI,*)
      ENDIF

C     Finally, write the nonzero induced dipoles to file
      IF (.NOT. FIXDIP) THEN
        CALL PUT_TO_FILE_1('INDUCED_DIPOLES',POLDIM,INDMOM)
      ENDIF


 1040 FORMAT(' Due to electronic density: ')
 1041 FORMAT(' Due to nuclei: ')
 1042 FORMAT(' Due to permanent multipoles: ')
 1050 FORMAT('   Induced dipole moments   ')
 1051 FORMAT(2X,'=',22('-'),'=',2X)
 1030 FORMAT(' Total induced dipole moments: ')
 1000 FORMAT(1X,51('='))
 1010 FORMAT(' | Site  |      X      |      Y      |      Z      |')
 1020 FORMAT(1X,I6,3(4X,F10.6))

      CALL QEXIT('GET_IND_DIPOLES_1')
      RETURN
      END
C******************************************************************************
C  /* Deck GET_CHARGE_ELFLD */
      SUBROUTINE GET_CHARGE_ELFLD(Q,XORI,YORI,ZORI,
     &                            XTAR,YTAR,ZTAR,
     &                            EFX,EFY,EFZ)
C
C     Calculates the electric field strength due to electric point
C     charge.
C
C     INPUT:
C
C       Q              - the magnitude of the point charge
C       XORI,YORI,ZORI - position of the point charge
C       XTAR,YTAR,ZTAR - position of the point where electric field is to be calculated
C       
C     OUTPUT:
C
C       EFX,EFY,EFZ    - components of the electric field strength vector
C
C KA, 2008 Oct. 22
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C

      DOUBLE PRECISION Q,XORI,YORI,ZORI
      DOUBLE PRECISION XTAR,YTAR,ZTAR
      DOUBLE PRECISION EFX,EFY,EFZ

      CALL QENTER('GET_CHARGE_ELFLD')

      EFX = 0.0D0
      EFY = 0.0D0
      EFZ = 0.0D0

      DIST2 = 0.0D0
      DIST2 = DIST2 + (XTAR - XORI)**2
      DIST2 = DIST2 + (YTAR - YORI)**2
      DIST2 = DIST2 + (ZTAR - ZORI)**2
      DIST  = SQRT(DIST2)
      DIST3 = DIST**3

      EFX = Q*(XTAR - XORI)/DIST3
      EFY = Q*(YTAR - YORI)/DIST3
      EFZ = Q*(ZTAR - ZORI)/DIST3

      CALL QEXIT('GET_CHARGE_ELFLD')

      RETURN
      END
C******************************************************************************
C  /* Deck GET_DIPOLE_ELFLD */
      SUBROUTINE GET_DIPOLE_ELFLD(MJUX,MJUY,MJUZ,
     &                            XORI,YORI,ZORI,
     &                            XTAR,YTAR,ZTAR,
     &                            EFX,EFY,EFZ)
C
C     Calculates the electric field strength due to electric dipole
C     moment.
C
C     INPUT:
C
C       MJUX,MJUY,MJUZ - the components of the dipole moment
C       XORI,YORI,ZORI - position of the dipole moment
C       XTAR,YTAR,ZTAR - position of the point where electric field is
C                        to be calculated
C       
C     OUTPUT:
C
C       EFX,EFY,EFZ    - components of the electric field strength
C                        vector
C
C KA, 2008 Oct. 22
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C
      DOUBLE PRECISION MJUX,MJUY,MJUZ
      DOUBLE PRECISION XORI,YORI,ZORI
      DOUBLE PRECISION XTAR,YTAR,ZTAR
      DOUBLE PRECISION EFX,EFY,EFZ

      CALL QENTER('GET_DIPOLE_ELFLD')

      EFX = 0.0D0
      EFY = 0.0D0
      EFZ = 0.0D0

      DIST2 = 0.0D0
      DIST2 = DIST2 + (XTAR - XORI)**2
      DIST2 = DIST2 + (YTAR - YORI)**2
      DIST2 = DIST2 + (ZTAR - ZORI)**2
      DIST  = SQRT(DIST2)
      DIST3 = DIST**3
      DIST5 = DIST**5

      EFX = EFX + MJUX*((3*(XTAR - XORI)*(XTAR - XORI))/DIST5 -
     &      (1.0/DIST3))
      EFX = EFX + MJUY* (3*(XTAR - XORI)*(YTAR - YORI))/DIST5
      EFX = EFX + MJUZ* (3*(XTAR - XORI)*(ZTAR - ZORI))/DIST5

      EFY = EFY + MJUX* (3*(YTAR - YORI)*(XTAR - XORI))/DIST5
      EFY = EFY + MJUY*((3*(YTAR - YORI)*(YTAR - YORI))/DIST5 -
     &      (1.0/DIST3))
      EFY = EFY + MJUZ* (3*(YTAR - YORI)*(ZTAR - ZORI))/DIST5

      EFZ = EFZ + MJUX* (3*(ZTAR - ZORI)*(XTAR - XORI))/DIST5
      EFZ = EFZ + MJUY* (3*(ZTAR - ZORI)*(YTAR - YORI))/DIST5
      EFZ = EFZ + MJUZ*((3*(ZTAR - ZORI)*(ZTAR - ZORI))/DIST5 -
     &      (1.0/DIST3))

      CALL QEXIT('GET_DIPOLE_ELFLD')

      RETURN
      END
C******************************************************************************
C  /* Deck GET_QUADRUPOLE_ELFLD */
      SUBROUTINE GET_QUADRUPOLE_ELFLD(QXX,QXY,QXZ,
     &                                QYY,QYZ,QZZ,
     &                                XORI,YORI,ZORI,
     &                                XTAR,YTAR,ZTAR,
     &                                EFX,EFY,EFZ)
C
C     Calculates the electric field strength due to electric quadrupole
C     moment.
C
C     INPUT:
C
C       QXX,QXY,QXZ,QYY,QYZ,QZZ - the components of the symmetric
C                                 quadrupole moment tensor
C       XORI,YORI,ZORI          - position of the quadrupole moment
C       XTAR,YTAR,ZTAR          - position of the point where electric field is
C                                 to be calculated
C       
C     OUTPUT:
C
C       EFX,EFY,EFZ             - components of the electric field strength
C                                 vector
C
C KA, 2008 Oct. 22
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C
      DOUBLE PRECISION QXX,QXY,QXZ,QYY,QYZ,QZZ
      DOUBLE PRECISION XORI,YORI,ZORI
      DOUBLE PRECISION XTAR,YTAR,ZTAR
      DOUBLE PRECISION EFX,EFY,EFZ

      DOUBLE PRECISION QTENS,ELFVEC,CORDO,CORDT
      DIMENSION QTENS(3,3),ELFVEC(3),CORDO(3),CORDT(3)

      CALL QENTER('GET_QUADRUPOLE_ELFLD')

      EFX = 0.0D0
      EFY = 0.0D0
      EFZ = 0.0D0

      DIST2 = 0.0D0
      DIST2 = DIST2 + (XTAR - XORI)**2
      DIST2 = DIST2 + (YTAR - YORI)**2
      DIST2 = DIST2 + (ZTAR - ZORI)**2
      DIST  = SQRT(DIST2)
      DIST5 = DIST**5
      DIST7 = DIST**7

      QTENS(1,1) = QXX
      QTENS(1,2) = QXY
      QTENS(1,3) = QXZ
      QTENS(2,1) = QXY
      QTENS(2,2) = QYY
      QTENS(2,3) = QYZ
      QTENS(3,1) = QXZ
      QTENS(3,2) = QYZ
      QTENS(3,3) = QZZ

      CORDO(1) = XORI
      CORDO(2) = YORI
      CORDO(3) = ZORI

      CORDT(1) = XTAR
      CORDT(2) = YTAR
      CORDT(3) = ZTAR

      ELFVEC(1) = 0.0D0
      ELFVEC(2) = 0.0D0
      ELFVEC(3) = 0.0D0

      DO 100 I=1,3
        DO 110 J=1,3
          DO 120 K=1,3

            ELEM = 0.0D0
            ELEM = (15*(CORDT(K) - CORDO(K))*
     &                 (CORDT(J) - CORDO(J))*
     &                 (CORDT(I) - CORDO(I)))/
     &                  DIST7
            IF (K .EQ. J) THEN
              ELEM = ELEM - (3*(CORDT(I) - CORDO(I))/DIST5)
            ENDIF
            IF (I .EQ. K) THEN
              ELEM = ELEM - (3*(CORDT(J) - CORDO(J))/DIST5)
            ENDIF
            IF (I .EQ. J) THEN
              ELEM = ELEM - (3*(CORDT(K) - CORDO(K))/DIST5)
            ENDIF
            ELEM = ELEM*QTENS(K,J)

            ELFVEC(I) = ELFVEC(I) + ELEM

 120      CONTINUE
 110    CONTINUE
        ELFVEC(I) = ELFVEC(I)/2.0
 100  CONTINUE

      EFX = ELFVEC(1)
      EFY = ELFVEC(2)
      EFZ = ELFVEC(3)

      CALL QEXIT('GET_QUADRUPOLE_ELFLD')

      RETURN
      END
C******************************************************************************
C  /* Deck Put_To_File_1 */
      SUBROUTINE PUT_TO_FILE_1(FLNAME,NULOOP,DDATA)
C
#include "implicit.h"
#include "dummy.h"
C
      CHARACTER*(*) FLNAME
      INTEGER   NMBU,NULOOP
      DIMENSION DDATA(*)
C
      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'UNKNOWN',' ','FORMATTED',IDUMMY,.FALSE.)
C
      REWIND (NMBU)
      LM = 1
      DO 820 L = 1,NULOOP
        WRITE(NMBU,'(I5,3E25.15)') L,DDATA(LM),DDATA(LM+1),DDATA(LM+2)
        LM = LM + 3
  820 CONTINUE
C
      CALL GPCLOSE(NMBU,'KEEP')
C
      END
C
C******************************************************************************
C**************************************************************
C  /* Deck Get_From_File_1 */
      SUBROUTINE GET_FROM_FILE_1(FLNAME,NULOOP,DDATA)
C**************************************************************
C
#include "implicit.h"
#include "dummy.h"
C
      CHARACTER*(*) FLNAME
      INTEGER   NMBU,NULOOP
      DIMENSION DDATA(*)
C
      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'UNKNOWN',' ','FORMATTED',IDUMMY,.FALSE.)
C
      REWIND (NMBU)
      LM = 1
      DO 820 L = 1,NULOOP
        READ(NMBU,'(I5,3E25.15)') LK,DDATA(LM),DDATA(LM+1),DDATA(LM+2)
        LM = LM + 3
  820 CONTINUE
C
      IF (LK.NE.NULOOP) THEN
        CALL QUIT('Problem in dimension in GET_FROM_FILE_1')
      ENDIF

      CALL GPCLOSE(NMBU,'KEEP')
C
      END
C
C******************************************************************************
C  /* Deck MM_PROPS */
      SUBROUTINE MM_PROPS(WRK,LWRK,IPRINT)
C
C  Calculates properties of the MM region.
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "qm3.h"
#include "qmmm.h"
#include "infpri.h"

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0) 

      DIMENSION WRK(LWRK)

      LOGICAL LOCDEB,FNDLAB

      CALL QENTER('MM_PROPS')

      LOCDEB = .FALSE. 

      WRITE(LUPRI,*) ' -------------------------------------- '
      WRITE(LUPRI,*) '     Output from MM property module     '
      WRITE(LUPRI,*) ' ---------------------------------------'
      WRITE(LUPRI,*)

      KINVMAT    = 1
      KFULLMAT   = KINVMAT  + 3*NNZAL*(3*NNZAL+1)/2
      KBMATS     = KFULLMAT  + (3*NNZAL)*(3*NNZAL)
      KPOLMAT    = KBMATS   + (3*NNZAL)*3
      KPOLCORD   = KPOLMAT  + 3*3
      KEND       = KPOLCORD + 3*NNZAL
      LWRK1      = LWRK     - KEND

      IF (LWRK1 .LT. 0) CALL ERRWRK('MM_PROPS',-KEND,LWRK)

      CALL DZERO(WRK(KINVMAT),3*NNZAL*(3*NNZAL+1)/2)
      CALL DZERO(WRK(KFULLMAT),(3*NNZAL)*(3*NNZAL))
      CALL DZERO(WRK(KBMATS),(3*NNZAL)*3)
      CALL DZERO(WRK(KPOLMAT),(3*3))

      CALL MM_DIPANDCHARGE(WRK(KEND),LWRK1,IPRINT)

      IF (MMMAT) THEN
        CALL MM_POLARIZABILITY(WRK(KINVMAT),WRK(KFULLMAT),WRK(KBMATS),
     &                         WRK(KPOLMAT),IPRINT)

        CALL MM_OPTROT(WRK(KINVMAT),WRK(KFULLMAT),WRK(KPOLCORD),IPRINT)
      ELSE 
        WRITE(LUPRI,*) 'MM properties skipped since MMITER' 
      ENDIF

      WRITE(LUPRI,*) ' ---------------------------------------' 
      WRITE(LUPRI,*)

      CALL QEXIT('MM_PROPS')
      RETURN
      END
C******************************************************************************
C  /* Deck MM_POLARIZABILITY */
      SUBROUTINE MM_POLARIZABILITY(INVMAT,FULLMAT,BMATS,POLMAT,IPRINT)
C
C  Contracts the Relay matrix to the group and molecular
C  polarizabilities
C
#include "implicit.h"
#include "priunit.h"
#include "infpri.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"

      LOGICAL FNDLAB,LOCDEB
      DOUBLE PRECISION INVMAT,FULLMAT,BMATS,POLMAT
      DIMENSION INVMAT(3*NNZAL*(3*NNZAL+1)/2)
      DIMENSION FULLMAT(3*NNZAL,3*NNZAL)
      DIMENSION BMATS(3*NNZAL,3)
      DIMENSION POLMAT(3,3)

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0) 

      CALL QENTER('MM_POLARIZABILITY')

      LOCDEB = .FALSE. 

C     Read the relay matrix from file

      LUQMMM = -1
      IF (LUQMMM .LT. 0) THEN
        CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &             'UNFORMATTED',IDUMMY,.FALSE.)
      ENDIF
      REWIND(LUQMMM)

      N = 3*NNZAL
      IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
        CALL READT(LUQMMM,N*(N+1)/2,INVMAT)
      ELSE
        CALL QUIT('Problem reading the matrix from the QMMMIM file.')
      ENDIF

      CALL GPCLOSE(LUQMMM,'KEEP')

      L = 1
      DO J = 1, N
        K = J*(J-1)/2
        DO I = J, N
          FULLMAT(I,J) = INVMAT(L)
          L = L + 1
        END DO
        M = J*N-K
        L = 1 + M
      END DO

      DO I = 2, N
        DO J = 1, I-1
          FULLMAT(J,I) = FULLMAT(I,J)
        ENDDO
      ENDDO

      IF ( (IPRINT .GE. 15) .OR. (LOCDEB) ) THEN
        WRITE(LUPRI,*) 'Relay mat. is read from file MM_POLARIZABILITY'
        CALL OUTPUT(FULLMAT,1,N,1,N,N,N,1,LUPRI)
      ENDIF

C     Contract the Relay matrix

      K1=1
      DO 100 K = 1,NNZAL
        J1 = 1
        DO 101 J = 1,NNZAL
          BMATS(K1,1)    = BMATS(K1,1) + FULLMAT(K1,J1)
          BMATS(K1,2)    = BMATS(K1,2) + FULLMAT(K1,J1+1)
          BMATS(K1,3)    = BMATS(K1,3) + FULLMAT(K1,J1+2)
          BMATS(K1+1,1)  = BMATS(K1+1,1) + FULLMAT(K1+1,J1)
          BMATS(K1+1,2)  = BMATS(K1+1,2) + FULLMAT(K1+1,J1+1)
          BMATS(K1+1,3)  = BMATS(K1+1,3) + FULLMAT(K1+1,J1+2)
          BMATS(K1+2,1)  = BMATS(K1+2,1) + FULLMAT(K1+2,J1)
          BMATS(K1+2,2)  = BMATS(K1+2,2) + FULLMAT(K1+2,J1+1)
          BMATS(K1+2,3)  = BMATS(K1+2,3) + FULLMAT(K1+2,J1+2)
          J1 = J1 + 3
  101   CONTINUE

        IF (LOCDEB) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) 'Polarizability for group ', K
          WRITE(LUPRI,*) BMATS(K1,1),BMATS(K1,2),BMATS(K1,3)
          WRITE(LUPRI,*) BMATS(K1+1,1), BMATS(K1+1,2), BMATS(K1+1,3)
          WRITE(LUPRI,*) BMATS(K1+2,1), BMATS(K1+2,2), BMATS(K1+2,3) 
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) 'Isotropic polarizability '
          TEMP = BMATS(K1,1)+BMATS(K1+1,2)+BMATS(K1+2,3)
          WRITE(LUPRI,*) 1.0D0/3.0D0*TEMP
          WRITE(LUPRI,*)
        ENDIF

        K1 = K1 +3

  100 CONTINUE

C     Contract to molecular polarizability

      K1=1
      DO 102 J = 1,NNZAL
        POLMAT(1,1)  = POLMAT(1,1) + BMATS(K1,1)
        POLMAT(1,2)  = POLMAT(1,2) + BMATS(K1,2)
        POLMAT(1,3)  = POLMAT(1,3) + BMATS(K1,3)
        POLMAT(2,1)  = POLMAT(2,1) + BMATS(K1+1,1)
        POLMAT(2,2)  = POLMAT(2,2) + BMATS(K1+1,2)
        POLMAT(2,3)  = POLMAT(2,3) + BMATS(K1+1,3)
        POLMAT(3,1)  = POLMAT(3,1) + BMATS(K1+2,1)
        POLMAT(3,2)  = POLMAT(3,2) + BMATS(K1+2,2)
        POLMAT(3,3)  = POLMAT(3,3) + BMATS(K1+2,3)
        K1 = K1 + 3
  102 CONTINUE

      N=3
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'Molecular polarizability of the MM region'
      CALL OUTPUT(POLMAT,1,N,1,N,N,N,1,LUPRI)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'Isotropic polarizability '
      TEMP = POLMAT(1,1)+POLMAT(2,2)+POLMAT(3,3)
      WRITE(LUPRI,*) 1.0D0/3.0D0*TEMP
      WRITE(LUPRI,*)

      XI = FLOAT(NNZAL)
      XXI = DBLE(XI)
      TEMP = 1.0D0/3.0D0*TEMP/XXI
      WRITE(LUPRI,*) 'Isotropic polarizability pr. pol. site'
      WRITE(LUPRI,*) TEMP
      WRITE(LUPRI,*)

      CALL QEXIT('MM_POLARIZABILITY')
      RETURN
      END
C******************************************************************************
C  /* Deck MM_DIPANDCHARGE */
      SUBROUTINE MM_DIPANDCHARGE(WRK,LWRK,IPRINT)
C
C     Calculates the MM total charge and dipole moment
C
#include "implicit.h"
#include "priunit.h"
#include "infpri.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)
      DIMENSION WRK(LWRK)
      LOGICAL LOCDEB

      CALL QENTER('MM_DIPANDCHARGE')

      LOCDEB = .FALSE.

      KINDMOM = 1                   
      KLAST   =  KINDMOM + 3*NNZAL 
      LWRK1    = LWRK - KLAST

      IF (LWRK1 .LT. 0) CALL ERRWRK('MM_DIPANDCHARGE',-KLAST,LWRK)

      CALL DZERO(WRK(KINDMOM),3*NNZAL)

      XDIPIND = 0.0D0
      YDIPIND = 0.0D0
      ZDIPIND = 0.0D0

      IF (IPOLTP .GT. 0) THEN

        IF (LOCDEB) THEN
           WRITE(LUPRI,*)
           WRITE(LUPRI,*) 'Ind. dips read from file in MM_DIPANDCHARGE'
           WRITE(LUPRI,*)
        ENDIF

        CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WRK(KINDMOM))

C       Add induced dipoles  

        IJ = 0
        DO 100 I=1,NNZAL
          XDIPIND = XDIPIND + WRK(KINDMOM+IJ+0)
          YDIPIND = YDIPIND + WRK(KINDMOM+IJ+1)
          ZDIPIND = ZDIPIND + WRK(KINDMOM+IJ+2)
          IJ = IJ +3
  100   CONTINUE

      ENDIF

C     Add permanent dipoles

      XDIPP = 0.0D0
      YDIPP = 0.0D0
      ZDIPP = 0.0D0

      IF (NMULT .GE. 1) THEN

        DO 101 I=1,MMCENT
          XDIPP = XDIPP + MUL1MM(1,I)
          YDIPP = YDIPP + MUL1MM(2,I)
          ZDIPP = ZDIPP + MUL1MM(3,I)
  101   CONTINUE

      ENDIF

C     Add charges

      QMMT = 0.0D0
      XQ   = 0.0D0
      YQ   = 0.0D0
      ZQ   = 0.0D0

      IF (NMULT .GE. 0) THEN
        
        DO 102 I=1,MMCENT
          QMMT = QMMT + MUL0MM(I)
          XQ = XQ + MMCORD(1,I)*MUL0MM(I)
          YQ = YQ + MMCORD(2,I)*MUL0MM(I)
          ZQ = ZQ + MMCORD(3,I)*MUL0MM(I)
  102   CONTINUE

      ENDIF

      IF (NMULT .GE. 0) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' MM total charge: ', QMMT
        IF (ABS(QMMT) .GT. THRMM) THEN
          WRITE(LUPRI,*) ' The MM region is charged '
        ENDIF
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' MM total charge dipole moment (x,y,z): '
        WRITE(LUPRI,*)   XQ,YQ,ZQ
        WRITE(LUPRI,*)
      ENDIF

      IF (NMULT .GE. 1) THEN
        WRITE(LUPRI,*) ' MM total permanent dipole moment (x,y,z): '
        WRITE(LUPRI,*)   XDIPP,YDIPP,ZDIPP
        WRITE(LUPRI,*)
      ENDIF

      IF (IPOLTP .GT. 0) THEN 
        WRITE(LUPRI,*) ' MM total induced dipole moment (x,y,z): '
        WRITE(LUPRI,*)   XDIPIND,YDIPIND,ZDIPIND
        WRITE(LUPRI,*)
      ENDIF

C     Add all contributions to the dipule moment

      XDIP = XQ+XDIPP+XDIPIND
      YDIP = YQ+YDIPP+YDIPIND
      ZDIP = ZQ+ZDIPP+ZDIPIND

      IF ( (NMULT .GE. 0) .OR. (IPOLTP .GT. 0) ) THEN
        WRITE(LUPRI,*) ' MM total dipole moment (x,y,z): '
        WRITE(LUPRI,*)   XDIP,YDIP,ZDIP
        WRITE(LUPRI,*)
      ENDIF

      CALL QEXIT('MM_DIPANDCHARGE')
      RETURN
      END
C******************************************************************************
C  /* Deck MM_OPTROT */
      SUBROUTINE MM_OPTROT(INVMAT,FULLMAT,POLCORD,IPRINT)
C
C  Contracts the Relay matrix to the molecular optical rotation (beta)
C
#include "implicit.h"
#include "priunit.h"
#include "infpri.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"

      LOGICAL FNDLAB,LOCDEB
      DOUBLE PRECISION INVMAT,FULLMAT,BMAT,TEMP
      DIMENSION FULLMAT(3*NNZAL,3*NNZAL)
      DIMENSION POLCORD(3,NNZAL),BMAT(3,3)
      DIMENSION INVMAT(3*NNZAL*(3*NNZAL+1)/2)

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0) 

      CALL QENTER('MM_OPTROT')

      LOCDEB = .FALSE. 

C     Read the relay matrix from file

      LUQMMM = -1
      IF (LUQMMM .LT. 0) THEN
        CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &             'UNFORMATTED',IDUMMY,.FALSE.)
      ENDIF
      REWIND(LUQMMM)

      N = 3*NNZAL
      IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
        CALL READT(LUQMMM,N*(N+1)/2,INVMAT)
      ELSE
        CALL QUIT('Problem reading the matrix from the QMMMIM file.')
      ENDIF

      CALL GPCLOSE(LUQMMM,'KEEP')

      L = 1
      DO J = 1, N
        K = J*(J-1)/2
        DO I = J, N
          FULLMAT(I,J) = INVMAT(L)
          L = L + 1
        END DO
        M = J*N-K
        L = 1 + M
      END DO

      DO I = 2, N
        DO J = 1, I-1
          FULLMAT(J,I) = FULLMAT(I,J)
        ENDDO
      ENDDO

      IF ( (IPRINT .GE. 15) .OR. (LOCDEB) )THEN
        WRITE(LUPRI,*) 'Response mat. is read from file MM_OPTROT'
        CALL OUTPUT(FULLMAT,1,N,1,N,N,N,1,LUPRI)
      ENDIF

C     Construct an array of coordinates having polarizabilities

      IL = 1
      DO 100 I=1,MMCENT

        IF (ZEROAL(I) .EQ. -1) GOTO 100

        POLCORD(1,IL) = MMCORD(1,I)
        POLCORD(2,IL) = MMCORD(2,I)
        POLCORD(3,IL) = MMCORD(3,I)

        IL = IL + 1

 100  CONTINUE

      IF ( (IL-1) .NE. NNZAL) THEN
        CALL QUIT('Problem in coordinate dimension in MM_OPTROT.')
      ENDIF

      BETA = 0.0D0
      DO 101 I=1,NNZAL-1
        DO 102 J=I+1,NNZAL
    
          K=(I-1)*3+1
          L=(J-1)*3+1
          BMAT(1,1) = FULLMAT(K,L)
          BMAT(1,2) = FULLMAT(K,L+1)
          BMAT(1,3) = FULLMAT(K,L+2)
          BMAT(2,1) = FULLMAT(K+1,L)
          BMAT(2,2) = FULLMAT(K+1,L+1)
          BMAT(2,3) = FULLMAT(K+1,L+2)
          BMAT(3,1) = FULLMAT(K+2,L)
          BMAT(3,2) = FULLMAT(K+2,L+1)
          BMAT(3,3) = FULLMAT(K+2,L+2)
          XDIST = POLCORD(1,J) - POLCORD(1,I)
          YDIST = POLCORD(2,J) - POLCORD(2,I)
          ZDIST = POLCORD(3,J) - POLCORD(3,I)

          BETA = BETA + XDIST*(BMAT(3,2)-BMAT(2,3))
     *                + YDIST*(BMAT(1,3)-BMAT(3,1))
     *                + ZDIST*(BMAT(2,1)-BMAT(1,2))

 102    CONTINUE
 101  CONTINUE

      BETA = D6I*BETA

      WRITE(LUPRI,*) 'Isotropic OPTROT (beta)'
      WRITE(LUPRI,*)  BETA
      WRITE(LUPRI,*)
c
      CALL QEXIT('MM_OPTROT')
      RETURN
      END
C******************************************************************************
C  /* Deck GET_IND_DIPOLES_2 */
      SUBROUTINE GET_IND_DIPOLES_2(DCAO,POLDIM,INDMOM,
     &                             WRK,LWRK,IPRINT)
C
C A subroutine to calculate induced dipole moments by simple Jacobi iteration
C
C Input:
C
C   DCAO    - density matrix in AO basis
C   POLDIM  - the number of polarizable MM centers. 
C             (Actually in common as NNZAL....)
C
C Output:
C
C   INDMOM  - a vector containing induced dipole moments
C
C From Common   
C
C   ZEROAL  - a vector containing +1 for polarizable MM centers and -1
C             for non-polarizable
C
C Sep 2010 - JMO & KS:
C Started sharing of DFT/MM and CC/MM field routines
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "mmtimes.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
#include "infpar.h"

      LOGICAL EXCENT,LOCDEB,DIPCON, LSKIP
      LOGICAL TOFILE,TRIMAT,EXP1VL
      INTEGER POLDIM
      DOUBLE PRECISION INDMOM
      DIMENSION INDMOM(3*POLDIM),WRK(LWRK), DCAO(*)

      DOUBLE PRECISION EVEC,TTENS,ATMAT,DIP
      DIMENSION EVEC(3),TTENS(3,3)
      DIMENSION ATMAT(3,3),DIP(3)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CHARACTER*8 LABINT(9*MXCENT)

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('GET_IND_DIPOLES_2')

      LOCDEB = .FALSE.

      IF (POLDIM .NE. NNZAL) THEN
        WRITE(LUPRI,*) 'ERROR in no. of polarizabilities'
        CALL QUIT('ERROR in GET_IND_DIPOLES_2')
      ENDIF

      IF (SPLDIP) THEN
        WRITE(LUPRI,*) 'Split not implemented for iterative QMMM'
      ENDIF

C     Allocate memory for electric field integrals and electric fields
      KMAT    = 1                      ! For electric field integrals
      KELF    = KMAT   + 3*NNBASX      ! For total electric field     
      KEND    = KELF   + 3*POLDIM     
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('GET_IND_DIPOLES_2',-KEND,LWRK)

      CALL DZERO(WRK(KMAT),3*NNBASX)
      CALL DZERO(WRK(KELF),3*POLDIM)

C     1. Form F vector due to permanent MM moments

      IF (MMTIME) DTIME = SECOND()
#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
        CALL MM_FIELD_M2(DCAO(1),WRK(KELF),POLDIM,
     &                       WRK(KEND),LWRK1,IPRINT)
      ELSE
#endif
        LRI = 1   

        DO 200 I=1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 200

          CALL GET_FIELD(I,LRI,WRK(KELF),WRK(KEND),WRK(KEND),
     *                   DCAO,LOCDEB,WRK(KEND),LWRK1)
          LRI = LRI + 3

 200    CONTINUE

#if defined(VAR_MPI)
      ENDIF
#endif
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMPOL2 = TMMPOL2 + DTIME
      ENDIF

      NDIM = 3*POLDIM

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'F-Vector'
        DO 899 I=1,NDIM
        WRITE(LUPRI,*) WRK(KELF+I-1)
 899    CONTINUE       
      ENDIF

C     Convert the F-vector into induced dipole moments 

      IOPT = 1 ! read file with ind. momens from previous SCF iteration.
      IF (MMTIME) DTIME = SECOND()
      CALL F2QMMM(WRK(KELF),POLDIM,INDMOM,WRK(KEND),LWRK1,
     *            IOPT,IPRINT)

      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMF2 = TMMF2 + DTIME
      ENDIF
      IF (IPRINT .GT. 1) THEN
C       Write induced moments at each MM site to the DAL.OUT file
        WRITE(LUPRI,*)
        WRITE(LUPRI,1030)
        WRITE(LUPRI,*)
        WRITE(LUPRI,1000)
        WRITE(LUPRI,1010)
        WRITE(LUPRI,1000)
      ENDIF

      IINIM = 1

      DO 500 I=1,MMCENT
        IF (ZEROAL(I) .EQ. -1) THEN
          DIPX = 0.0D0
          DIPY = 0.0D0
          DIPZ = 0.0D0
        ELSE
          DIPX = INDMOM(IINIM+0)
          DIPY = INDMOM(IINIM+1)
          DIPZ = INDMOM(IINIM+2)
          IINIM = IINIM + 3
        ENDIF
        IF (IPRINT .GT. 1) WRITE(LUPRI,1020) I,DIPX,DIPY,DIPZ
 500  CONTINUE

      IF (IPRINT .GT. 1) THEN
        WRITE(LUPRI,1000)
        WRITE(LUPRI,*)
      ENDIF

C     Write the nonzero induced dipoles to file
      IF (.NOT. FIXDIP) THEN
        CALL PUT_TO_FILE_1('INDUCED_DIPOLES',POLDIM,INDMOM)
      ENDIF

 1050 FORMAT('   Induced dipole moments   ')
 1051 FORMAT(2X,'=',22('-'),'=',2X)
 1030 FORMAT(' Total induced dipole moments: ')
 1000 FORMAT(1X,51('='))
 1010 FORMAT(' | Site  |      X      |      Y      |      Z      |')
 1020 FORMAT(1X,I6,3(4X,F10.6))

      CALL QEXIT('GET_IND_DIPOLES_2')
      RETURN
      END
C******************************************************************************
C  /* Deck F2QMMM */
      SUBROUTINE F2QMMM(ELF,POLDIM,INDMOM,WRK,LWRK,IOPT,IPRINT)
C
C Converts a field vector into induced dipoles using iterative procedures.
C
C     Input: ELF
C     Output: INDMOM
C
C     INDMOM is the induced dipole moments
C     INDDIA is the diagonal part of the induced dipole moments,
C            i.e. the part corresponding directly to the F ELF vector.
C JK

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "mmtimes.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "infpar.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"


      LOGICAL EXCENT,LOCDEB,DIPCON
      INTEGER POLDIM
      DOUBLE PRECISION INDMOM,ELF
      DIMENSION INDMOM(3*POLDIM),ELF(3*POLDIM)
      DIMENSION WRK(LWRK)

      DOUBLE PRECISION AMAT,EVEC,MY0,TTENS,ATMAT,DIP
      DOUBLE PRECISION MY
      DIMENSION AMAT(3,3),EVEC(3),MY0(3),TTENS(3,3)
      DIMENSION ATMAT(3,3),DIP(3),MY(3)

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      LOGICAL FIRST
      SAVE    FIRST
      DATA    FIRST /.TRUE./

      CALL QENTER('F2QMMM')

      BTIME = SECOND()

      LOCDEB = .FALSE.

c      IF (IOPT .EQ. 1) THRESL = THMMIT
c      IF (IOPT .EQ. 2) THRESL = SQRT(THMMIT)/10.0D0

      THRESL = THMMIT

      IF (FIRST) NMMAC = 0

      NDIM = 3*POLDIM

      KINDP   = 1                   ! For the previos induced dipole (super) vector
      KINDDIA = KINDP   + 3*POLDIM  ! For the diagonal part of the induced moments
      KEND    = KINDDIA + 3*POLDIM
      LWRK1   = LWRK   - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('F2QMMM 1',-KEND,LWRK)

      CALL DZERO(WRK(KINDP),3*POLDIM)
      CALL DZERO(WRK(KINDDIA),3*POLDIM)

      KVEC  = KEND
      IF (MMDIIS) THEN
        KEND = KVEC + (MXMMIT+1)*NDIM
        LWRK1 = LWRK - KEND + 1
        IF (LWRK1 .LT. 0) CALL ERRWRK('F2QMMM 2',-KEND,LWRK)

        CALL DZERO(WRK(KVEC),(MXMMIT+1)*NDIM)
      ENDIF

C     Convert the F-vector into induced dipole moments 
C     by neglecting the off diagonal elements (the T tensor)
C     These moments are used as the initial guess.

      LRI = 1

      DO 400 I=1,MMCENT

        IF (ZEROAL(I) .EQ. -1) GOTO 400

C       Get the polarizability tensor for this site
        DO 401 K=1,3
          DO 402 J=1,3
            AMAT(K,J)  = 0.0D0
 402      CONTINUE
 401    CONTINUE

        IF (IPOLTP .EQ. 1)  THEN
          DO 403 J=1,3
            AMAT(J,J) = POLIMM(I)
 403      CONTINUE
        ELSE IF (IPOLTP .EQ. 2)  THEN
          AMAT(1,1) = POLMM(1,I)
          AMAT(1,2) = POLMM(2,I)
          AMAT(1,3) = POLMM(3,I)
          AMAT(2,1) = POLMM(2,I)
          AMAT(2,2) = POLMM(4,I)
          AMAT(2,3) = POLMM(5,I)
          AMAT(3,1) = POLMM(3,I)
          AMAT(3,2) = POLMM(5,I)
          AMAT(3,3) = POLMM(6,I)
        ENDIF

C       Now get the F-vector for this site
        EVEC(1) = ELF(LRI+0)
        EVEC(2) = ELF(LRI+1)
        EVEC(3) = ELF(LRI+2)

C       Calculate the induced dipole moment
        NLDIM = 3
        NTOTI = MAX(NLDIM,1)
        CALL DGEMV('N',NLDIM,NLDIM,D1,AMAT,NTOTI,EVEC,1,D0,MY0,1)

        WRK(KINDDIA-1+LRI+0) = MY0(1)
        WRK(KINDDIA-1+LRI+1) = MY0(2)
        WRK(KINDDIA-1+LRI+2) = MY0(3)

        LRI = LRI + 3
 400  CONTINUE

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'My-Vector: Diagonal contribution'
        DO 404 I=1,NDIM
        WRITE(LUPRI,*) WRK(KINDDIA+I-1)
 404    CONTINUE
      ENDIF

      CALL DCOPY(NDIM,WRK(KINDDIA),1,INDMOM,1)

      IF (IOPT .EQ. 1) THEN
        IF (.NOT. FIRST) THEN
          CALL GET_FROM_FILE_1('INDUCED_DIPOLES',POLDIM,WRK(KINDP))
        ELSE 
          CALL DCOPY(NDIM,WRK(KINDDIA),1,WRK(KINDP),1)
        ENDIF
      ENDIF 

      IF (IOPT .EQ. 2) CALL DCOPY(NDIM,WRK(KINDDIA),1,WRK(KINDP),1)

      IF (MMDIIS) THEN
        CALL DCOPY(NDIM,WRK(KINDP),1,WRK(KVEC),1)
      ENDIF

      IF (LOCDEB) WRITE(LUPRI,*) 'Done generating the F-Vector'

C     Now iterate...  !

      IF (MMTIME) DTIME = SECOND()

      LM = 0
      DIPCON = .FALSE.
#if defined(VAR_MPI)
      IF (NODTOT .GE. 1) THEN
        CALL MMITER_INDDIP_M(POLDIM,WRK(KINDP),INDMOM,WRK(KVEC),
     *                WRK(KINDDIA),WRK(KEND),LWRK1,LOCDEB,DIPCON,LM)
      ELSE
#endif
        DO 999 ITER = 1, MXMMIT
          LM = LM + 1

          LRI = 1
          DO 405 I=1,MMCENT
            IF (ZEROAL(I) .EQ. -1) GOTO 405

            LCI = 1
            DO 409 J=1,MMCENT

              IF (ZEROAL(J) .EQ. -1) GOTO 409 

              CALL GET_MY(I,J,WRK(KINDP+LCI-1),MY)
              INDMOM(LRI+0) = INDMOM(LRI+0) + MY(1)
              INDMOM(LRI+1) = INDMOM(LRI+1) + MY(2)
              INDMOM(LRI+2) = INDMOM(LRI+2) + MY(3)
              LCI = LCI + 3
 409        CONTINUE

            LRI = LRI + 3
 405      CONTINUE

          TERROR=0.0D0
          DO 414 I=1,NDIM
            TERROR = TERROR + (INDMOM(I)-WRK(KINDP+I-1))*
     &                    (INDMOM(I)-WRK(KINDP+I-1)) 
 414      CONTINUE

          IF ( (LOCDEB) .OR. (IPRINT .GE. 15) ) THEN
            LMAX = 0
            TMAX = 0.0D0
            DO 413 I=1,NDIM
              TDIFF = ABS(INDMOM(I) - WRK(KINDP-1+I))
              IF (TDIFF .GT. TMAX) THEN
                TMAX = TDIFF
                LMAX = I 
              ENDIF
 413        CONTINUE
            IF (LMAX .NE. 0) THEN
              WRITE(LUPRI,*) 'Maximum deviation (element) is ',
     *                        TMAX, LMAX
            ENDIF 
          ENDIF

          IF (ABS(TERROR) .LT. THRESL) THEN
            DIPCON = .TRUE.
            GOTO 9000
          ELSE
            DIPCON = .FALSE. 
            IF (LOCDEB )WRITE(LUPRI,*) 'TERROR ',TERROR
            IF (MMDIIS) THEN
              CALL DCOPY(NDIM,INDMOM,1,WRK(KVEC+ITER*NDIM),1)
              CALL MM_DIIS_EXTRAPOLATION(WRK(KVEC),ITER,NDIM,WRK(KINDP),
     *                               WRK(KEND),LWRK1,IPRINT)
            ELSE 
              CALL DCOPY(NDIM,INDMOM,1,WRK(KINDP),1)
            ENDIF
C       If no convergence in last iteration keep the values for the
C       induced dipoles, i.e. not only the diagonal part
            IF (ITER .NE. MXMMIT) CALL DCOPY(NDIM,WRK(KINDDIA),1,
     *                                       INDMOM,1)
          ENDIF

 999    CONTINUE

 9000   CONTINUE !Done

#if defined(VAR_MPI)
      ENDIF !parallel mmiter
#endif
      IF (MMTIME) THEN
        DTIME = SECOND() - DTIME
        TMMITER = TMMITER + DTIME
      ENDIF

      LM = LM - 1
      IF (DIPCON) THEN
        IF (IPRINT .GT. 1) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Done with induced dipoles in ',LM,' iterations'
         WRITE(LUPRI,*)
        ENDIF 
      ELSE
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) 'WARNING: Induced dipoles NOT converged'
        WRITE(LUPRI,*)
      ENDIF

      NMMAC = NMMAC + LM
      IF (IPRINT .GT. 1) THEN
        WRITE(LUPRI,*) 'Acc. iterations:', NMMAC
      ENDIF

      IF (FIRST) FIRST = .FALSE.

      BTIME = SECOND() - BTIME
      TF2QMMM = TF2QMMM + BTIME

      CALL QEXIT('F2QMMM')
      RETURN
      END
C******************************************************************************
C  /* Deck MM_DIIS_EXTRAPOLATION */
      SUBROUTINE MM_DIIS_EXTRAPOLATION(VEC,ITER,NDIM,RESVEC,WRK,LWRK,
     *                                 IPRINT)
C
C     Find the optimal DIIS vector of previously iterated induced dipoles. 
C
C     Input: VEC, ITER, NDIM   
C     Output: RESVEC
C
C     VEC is the collection of previos induced dipole vectors
C     RESVEC is the result vector 
C     NDIM is 3*(the number of polarizable sites)
C     ITER is the iteration number
C JK

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"

      LOGICAL LOCDEB

      INTEGER NDIM,ITER

      DOUBLE PRECISION VEC,RESVEC

      DIMENSION VEC(NDIM,(MXMMIT+1))
      DIMENSION RESVEC(NDIM)
      DIMENSION WRK(LWRK)

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )

      CALL QENTER('MM_DIIS_EXTRAPOLATION')

      LOCDEB = .FALSE.

      IF (ITER .LE. MXMMDI) THEN
        NDIIS= ITER+1
        IOFF = 0
      ELSE
        NDIIS = MXMMDI+1
        IOFF  = ITER-MXMMDI
      ENDIF

      KDIIS   = 1                   
      KVECA   = KDIIS   + NDIIS*NDIIS  
      KPVT    = KVECA   + NDIIS
      KEND    = KPVT    + NDIIS
      LWRK1   = LWRK   - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('MM_DIIS_EXTRAPOLATION',-KEND,LWRK)

      CALL DZERO(WRK(KDIIS),NDIIS*NDIIS)
      CALL DZERO(WRK(KVECA),NDIIS)
      CALL DZERO(WRK(KPVT),NDIIS)

      WRK(KDIIS) = D0
      WRK(KVECA) = -1.0D0
   
      DO 100 I=1,NDIIS-1
        WRK(KDIIS+I) = -1.0D0
        WRK(KVECA+I) = D0
 100  CONTINUE

      DO 101 I=2,NDIIS
        DO 102 J=1,NDIIS
          IF (J .EQ. 1) THEN 
            WRK(KDIIS+(I-1)*NDIIS+(J-1)) = -1.0D0
          ELSE
            TEMP=0.0D0
            DO 103 K=1,NDIM
              TEMP = TEMP + (VEC(K,I+IOFF)-VEC(K,I-1+IOFF))*
     *                      (VEC(K,J+IOFF)-VEC(K,J-1+IOFF))
  103       CONTINUE 
            WRK(KDIIS+(I-1)*NDIIS+(J-1)) = TEMP
          ENDIF
 102    CONTINUE 
 101  CONTINUE

      IF (LOCDEB) THEN
        N=NDIIS
        WRITE(LUPRI,*) 'DIIS matrix in iteration ',ITER
        CALL OUTPUT(WRK(KDIIS),1,N,1,N,N,N,1,LUPRI)

        WRITE(LUPRI,*) 'B-DIIS Vector',ITER
        DO 104 I=1,NDIIS
          WRITE(LUPRI,*) WRK(KVECA+I-1)
 104    CONTINUE
      ENDIF

      CALL DGESV(NDIIS,1,WRK(KDIIS),NDIIS,WRK(KPVT),WRK(KVECA),
     *           NDIIS,INFO)
      IF (INFO .NE. 0) THEN
         CALL QUIT('Error in MM_DIIS_EXTRAPOLATION')
      END IF

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'A-DIIS Vector',ITER
        DO 105 I=1,NDIIS
          WRITE(LUPRI,*) WRK(KVECA+I-1)
 105    CONTINUE
      ENDIF

      TEMP = D0
      DO 106 I=2,NDIIS
        TEMP = TEMP + WRK(KVECA+I-1)
 106  CONTINUE

      IF (ABS(TEMP-D1) .GT. 1.0D-08) THEN
        WRITE(LUPRI,*) 'WARNING: Sum of lambdas in MM_DIIS is ', TEMP
      ENDIF

      CALL DZERO(RESVEC,NDIM)

      DO 107 I=2,NDIIS
       CALL DAXPY(NDIM,WRK(KVECA+I-1),VEC(1,I+IOFF),1,RESVEC,1)
 107  CONTINUE

      IF (LOCDEB .OR. (IPRINT .GE. 15)) THEN
        WRITE(LUPRI,*) 'Guess induced dipole vector from MM_DIIS',ITER
        DO 108 I=1,NDIM
          WRITE(LUPRI,*) RESVEC(I)
 108    CONTINUE
      ENDIF

      CALL DCOPY(NDIM,RESVEC,1,VEC(1,NDIIS+IOFF),1)

      IF (.FALSE.) THEN ! Damp procedure
        IF (ITER .GE. 2) THEN
          TEMP1 = 0.0D0
          TEMP2 = 0.0D0
          DO 200 I=1,NDIM
            TEMP1 = TEMP1 + (VEC(I,NDIIS)-VEC(I,NDIIS-1))**2
            TEMP2 = TEMP2 + (VEC(I,NDIIS)-VEC(I,NDIIS-2))**2
 200      CONTINUE
          TLAM1 = 1.0D0/TEMP1
          TLAM2 = 1.0D0/TEMP2
          TLAM  = TLAM1/(TLAM1+TLAM2) 
          TLAMM = 1.0D0-TLAM
          CALL DAXPY(NDIM,TLAM,VEC(1,NDIIS),1,RESVEC,1)
          CALL DAXPY(NDIM,TLAMM,VEC(1,NDIIS-1),1,RESVEC,1)
          CALL DCOPY(NDIM,RESVEC,1,VEC(1,NDIIS),1)
        ELSE
          CALL DAXPY(NDIM,1.0D0,VEC(1,NDIIS),1,RESVEC,1)
        ENDIF
      ENDIF

      CALL QEXIT('MM_DIIS_EXTRAPOLATION')
      RETURN
      END
C******************************************************************************
C  /* Deck MAKE_QMMM_INVERSE_RESPONSE_MATRIX */
      SUBROUTINE MAKE_QMMM_INVERSE_RESPONSE_MATRIX(INVMAT,POLDIM) ! Construct inverse response matrix
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C
      DOUBLE PRECISION AMATS
      DIMENSION AMATS(6)
      LOGICAL EXCENT
      INTEGER POLDIM, IPVT
      DOUBLE PRECISION INVMAT, WRKV
      DIMENSION INVMAT(3*POLDIM*(3*POLDIM+1)/2)
      DIMENSION IPVT(3)
      DIMENSION WRKV(3)

      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('MAKE_QMMM_INV_RESP_MATRIX')

      DO I=1,6
          AMATS(I)  = 0.0D0
      END DO

C     Construct packed inverse response matrix

      M = 0

      DO I = 1, MMCENT

        IF (ZEROAL(I) .EQ. -1) CYCLE

C       Isotropic polarizability is easy to invert
        IF ((IPOLTP .EQ. 1) .AND. (CONMAT)) THEN
          AMATS(1) = 1.0D0/POLIMM(I)
          AMATS(4) = 1.0D0/POLIMM(I)
          AMATS(6) = 1.0D0/POLIMM(I)
        ENDIF

C       Anisotropic polarizability inversion
        IF ((IPOLTP .EQ. 2) .AND. (CONMAT)) THEN
          AMATS(1) = POLMM(1,I)
          AMATS(2) = POLMM(2,I)
          AMATS(3) = POLMM(3,I)
          AMATS(4) = POLMM(4,I)
          AMATS(5) = POLMM(5,I)
          AMATS(6) = POLMM(6,I)

C         Factorization
          CALL DSPTRF('L', 3, AMATS, IPVT, INFO)
          IF (INFO .NE. 0) THEN
            DIST1 = 1.0D50
            DO K = 1, NUCIND
              DIST2 = SQRT((CORD(1,K)-MMCORD(1,I))**2 +
     *                     (CORD(2,K)-MMCORD(2,I))**2 +
     *                     (CORD(3,K)-MMCORD(3,I))**2)
              IF (DIST2 .LT. DIST1) THEN
                CLDIST = DIST2
                DIST1 = DIST2
              END IF
            END DO
            CLDIST = CLDIST*0.5291772108
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,*) 'WARNING: problems with
     & polarizability at site:', I
            WRITE(LUPRI,*) 'Distance to closest QM nucleus is (Å):',
     & CLDIST
            WRITE(LUPRI,*) 'Polarizability (xx, xy, xz, yy, yz, zz):'
            DO K = 1, 6
              WRITE(LUPRI,*) POLMM(K,I)
            END DO
            CALL QUIT('Error during factorization of polarizability!')
          END IF

C         Inversion
          CALL DSPTRI('L', 3, AMATS, IPVT, WRKV, INFO)
          IF (INFO .NE. 0) THEN
            CALL QUIT('Error during inversion of local polarizability')
          END IF
        ENDIF

        DO L = 3, 1, -1
          DO J = I, MMCENT
            IF (ZEROAL(J) .EQ. -1) CYCLE
            IF (J .EQ. I) THEN
              IF (L .EQ. 3) THEN
                DO K = 1, L
                  INVMAT(M+K) = AMATS(K)
                END DO
              ELSE IF (L .EQ. 2) THEN
                DO K = 1, L
                  INVMAT(M+K) = AMATS(3+K)
                END DO
              ELSE IF (L .EQ. 1) THEN
                  INVMAT(M+1) = AMATS(5+1)
              END IF
              M = M + L
            ELSE
              IF (NOMB) THEN
                DO K = 1, 3
                  INVMAT(M+K) = 0.0D0
                END DO
                M = M + 3
                CYCLE
              END IF

              R = 0.0D0; R2 = 0.0D0
              R2 = (MMCORD(1,I)-MMCORD(1,J))**2 +
     &             (MMCORD(2,I)-MMCORD(2,J))**2 +
     &             (MMCORD(3,I)-MMCORD(3,J))**2
              R = SQRT(R2)
              R3 = R**3
              R5 = R**5

              IF (R .GT. RCUTMM) THEN
                M = M + 3
                CYCLE
              ENDIF

              EXCENT = .FALSE.
              IF (NEWEXC) THEN
                DO N = 1, NEXLST
                  IF (EXLIST(1,I) .EQ. EXLIST(N,J)) EXCENT = .TRUE.
                ENDDO
              ELSE
                DO N = 1, NEXLST
                  IF (EXLIST(N,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
                END DO
              ENDIF

              IF (EXCENT) THEN
                DO K = 1, 3
                  INVMAT(M+K) = 0.0D0
                END DO
                M = M + 3
              ELSE

C               Include damping in the exponential form
C               JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
                IF (MMDAMP) THEN
                  IF (IPOLTP .EQ. 1) THEN
                    TEMPI = POLIMM(I)
                    TEMPJ = POLIMM(J)
                  ELSE IF (IPOLTP .EQ. 2) THEN
                    TEMPI =  (POLMM(1,I)+POLMM(4,I)
     &                       +POLMM(6,I))*D3I
                    TEMPJ =  (POLMM(1,J)+POLMM(4,J)
     &                       +POLMM(6,J))*D3I
                  ENDIF
                  TEMP = (TEMPI*TEMPJ)**D6I
                  SCREEN = 2.1304*R/TEMP
                  FE = 1.0D0-(1.0D0 + SCREEN + 0.5D0*SCREEN**2)
     &                       *EXP(-SCREEN)
                  FT = FE - (D6I*SCREEN**3)*EXP(-SCREEN)
                ELSE
                  FE = 1.0D0
                  FT = 1.0D0
                ENDIF

                IF (L .EQ. 3) THEN

                  DO K = 1, 3
                    T = FT*3.0D0*(MMCORD(1,I) - MMCORD(1,J))*
     &                        (MMCORD(K,I) - MMCORD(K,J))
                    T = T/R5
                    IF (K .EQ. 1) T = T - FE*1.0D0/R3
                    INVMAT(M+K) = -1.0D0*T
                  END DO

                ELSE IF (L .EQ. 2) THEN
                  DO K = 1, 3
                    T = FT*3.0D0*(MMCORD(2,I) - MMCORD(2,J))*
     &                        (MMCORD(K,I) - MMCORD(K,J))
                    T = T/R5
                    IF (K .EQ. 2) T = T - FE*1.0D0/R3
                    INVMAT(M+K) = -1.0D0*T
                  END DO
                ELSE IF (L .EQ. 1) THEN
                  DO K = 1, 3
                    T = FT*3.0D0*(MMCORD(3,I) - MMCORD(3,J))*
     &                        (MMCORD(K,I) - MMCORD(K,J))
                    T = T/R5
                    IF (K .EQ. 3) T = T - FE*1.0D0/R3
                    INVMAT(M+K) = -1.0D0*T
                  END DO
                END IF
                M = M + 3
              END IF
            END IF
          END DO
        END DO
      END DO

      CALL QEXIT('MAKE_QMMM_INV_RESP_MATRIX')

      RETURN
      END
C******************************************************************************
C 'Inside loops' routines (can be used both by parallel and sequential code) 
C Arnfinn Oct. 2010
C******************************************************************************
C  /* Deck charge_iter */
      SUBROUTINE CHARGE_ITER(I,DCAO,ENSEL,ENSNUC,LOCDEB,
     &                       TAO,WRK,LWRK,IPRTMP)
C
C     Calculate the energy contribution due to the charge on a MM cite
C
C     Input:
C       I      - MM cite I
C       DCAO   - density matrix
C       LOCDEB - local debugging
C       
C     Output:
C       ENSEL  - Energy due to QM electrons 
C       ENSNUC - Energy due to QM nuclear 
C       TAO    - Integrals
C
#include "implicit.h"
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "gnrinf.h"
#include "orgcom.h"
#include "priunit.h"

      DIMENSION WRK(LWRK), DCAO(NNBASX), TAO(NNBASX)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('CHARGE_ITER')

      IF (ABS(MUL0MM(I)) .LE. THRMM) THEN
         ENSEL  = 0.0D0
         ENSNUC = 0.0D0
         CALL QEXIT('CHARGE_ITER')
         RETURN
      ENDIF

      FAC1   =  1.0D0

      KMAT = 1
      KLAST = KMAT + NNBASX
      LWRK2 = LWRK - KLAST + 1
      IF (LWRK2 .LT. 0) CALL ERRWRK('CHARGE_ITER',-KLAST,LWRK)

      CALL DZERO(WRK(KMAT),NNBASX)

      KPATOM = 0
      NOSIM  = 1
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,I)
      DIPORG(2) = MMCORD(2,I)
      DIPORG(3) = MMCORD(3,I)

      IF (LOCDEB) THEN
C     Test for numerical int.
         CORZSAVE  = DIPORG(3)
         KMAT1  = KLAST
         KMAT2  = KMAT1 + NNBASX
         KMAT3  = KMAT2 + NNBASX
         KLAST1 = KMAT3 + NNBASX
         LWRK3  = LWRK - KLAST1 + 1

         IF (LWRK3 .LT. 0) CALL ERRWRK('CHARGE_ITER 2',-KLAST1,LWRK)

         CALL DZERO(WRK(KMAT1),3*NNBASX)
 
         DIPORG(3) = DIPORG(3) + 0.01
         RUNQM3=.TRUE.
         CALL GET1IN(WRK(KMAT1),'NPETES ',NOSIM,WRK(KLAST1),
     &                   LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
         DIPORG(3) = DIPORG(3) - 0.02
         CALL GET1IN(WRK(KMAT2),'NPETES ',NOSIM,WRK(KLAST1),
     &                   LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
         DIPORG(3) = CORZSAVE
         CALL GET1IN(WRK(KMAT3),'NPETES ',NOSIM,WRK(KLAST1),
     &                   LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
         RUNQM3=.FALSE.
C     Gradient        
         FM1 = -1.0D0
         FSCAL = 1.0D0/0.02
         CALL DAXPY(NNBASX,FM1,WRK(KMAT2),1,WRK(KMAT1),1)
         CALL DSCAL(NNBASX,FSCAL,WRK(KMAT1),1)
         FSCAL = -1.0D0
         CALL DSCAL(NNBASX,FSCAL,WRK(KMAT1),1)
         WRITE (LUPRI,'(/A)') 'E_z num matrix in QMMM_FCK_AO'
         CALL OUTPAK(WRK(KMAT1),NBAST,1,LUPRI)
         DIPORG(3) = CORZSAVE
      ENDIF

      RUNQM3=.TRUE.
      CALL GET1IN(WRK(KMAT),'NPETES ',NOSIM,WRK(KLAST),
     &               LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3=.FALSE.

      IF ( (IPRTMP.GT.15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A)') 'Pot. energy matrix in QMMM_CHARGE'
         CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
      ENDIF

      CALL DSCAL(NNBASX,MUL0MM(I),WRK(KMAT),1)
      EXPNS=DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
      ENSEL = EXPNS

      CALL DAXPY(NNBASX,FAC1,WRK(KMAT),1,TAO,1)

C     Now the QM nuclear contribution

      ECHCHL  = 0.0D0
      DO 101 J = 1,NUCIND
         XDIS   = CORD(1,J) - MMCORD(1,I) 
         YDIS   = CORD(2,J) - MMCORD(2,I) 
         ZDIS   = CORD(3,J) - MMCORD(3,I) 
         DIST2  = XDIS**2+YDIS**2+ZDIS**2
         DIST   = SQRT(DIST2) 
         ECHCHL = ECHCHL + CHARGE(J)*MUL0MM(I)/DIST
 101  CONTINUE

      ENSNUC = ECHCHL

      CALL QEXIT('CHARGE_ITER')

      RETURN
      END
C
C******************************************************************************
C  /* Deck dipole_iter */
      SUBROUTINE DIPOLE_ITER(I,DCAO,ENSEL,ENSNUC,LOCDEB,
     *                       TAO,WRK,LWRK,IPRTMP)

#include "implicit.h"
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "gnrinf.h"
#include "orgcom.h"
#include "priunit.h"

      DIMENSION WRK(LWRK), DCAO(NNBASX), TAO(NNBASX)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('DIPOLE_ITER')

C     See if the dipole moment at this site is zero
      DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2     
      DNORM = SQRT(DNORM2)
      IF (ABS(DNORM) .LE. THRMM) THEN
         ENSEL = 0.0D0
         ENSNUC = 0.0D0
         CALL QEXIT('DIPOLE_ITER')
         RETURN
      ENDIF 

      FAC1   =  1.0D0
      FACM1  = -1.0D0

      KMAT = 1
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1
      IF (LWRK2 .LT. 0) CALL ERRWRK('DIPOLE_ITER',-KLAST,LWRK)

      CALL DZERO(WRK(KMAT),3*NNBASX)

      KPATOM = 0
      NOSIM  = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,I)
      DIPORG(2) = MMCORD(2,I)
      DIPORG(3) = MMCORD(3,I)

      IF (LOCDEB) THEN
C     TEST for numerical int.
         CORZSAVE  = DIPORG(3)

         KMAT1  = KLAST
         KMAT2  = KMAT1 + 3*NNBASX
         KMAT3  = KMAT2 + 3*NNBASX
         KLAST1 = KMAT3 + 3*NNBASX
         LWRK3  = LWRK - KLAST1 + 1

         IF (LWRK3 .LT. 0) CALL ERRWRK('QMMM_DIPOLE 2',-KLAST1,LWRK)

         CALL DZERO(WRK(KMAT1),9*NNBASX)
          
         DIPORG(3) = DIPORG(3) + 0.01
         RUNQM3=.TRUE.
         CALL GET1IN(WRK(KMAT1),'NEFIELD',NOSIM,WRK(KLAST1),
     &                  LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                  KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
         DIPORG(3) = DIPORG(3) - 0.02
         CALL GET1IN(WRK(KMAT2),'NEFIELD',NOSIM,WRK(KLAST1),
     &                  LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                  KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
         DIPORG(3) = CORZSAVE
         CALL GET1IN(WRK(KMAT3),'NEFIELD',NOSIM,WRK(KLAST1),
     &                  LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                  KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
         RUNQM3=.FALSE.
C      Gradient
         FM1 = -1.0D0
         FSCAL = 1.0D0/0.02
         CALL DAXPY(3*NNBASX,FM1,WRK(KMAT2),1,WRK(KMAT1),1)
         CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT1),1)
         FSCAL = -1.0D0
         CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT1),1)
         WRITE (LUPRI,'(/A)') 'E_xz num matrix in QMMM_FCK_AO'
         CALL OUTPAK(WRK(KMAT1),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') 'E_zz num matrix in QMMM_FCK_AO'
         CALL OUTPAK(WRK(KMAT1+2*NNBASX),NBAST,1,LUPRI)
         
         DIPORG(3) = CORZSAVE
      ENDIF

      RUNQM3=.TRUE.
      CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3=.FALSE.

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF 
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 123 M=1,NUCIND
               DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                    (DIPORG(2)-CORD(2,M))**2 +
     &                    (DIPORG(3)-CORD(3,M))**2
               IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
               ENDIF
 123        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                (DIPORG(2)-QMCOM(2))**2 +
     &                (DIPORG(3)-QMCOM(3))**2
         ENDIF
         DIST = SQRT(DIST)

         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,I)+POLMM(4,I)+POLMM(6,I))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(I)
            ENDIF 
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIST/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIST))**3
         ENDIF
         CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
      ENDIF

      IF ( (IPRTMP.GT.15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A)') ' E_x_matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
         
         WRITE (LUPRI,'(/A)') ' E_y matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

         WRITE (LUPRI,'(/A)') ' E_z matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
      END IF

      CALL DSCAL(NNBASX,MUL1MM(1,I),WRK(KMAT),1)
      CALL DSCAL(NNBASX,MUL1MM(2,I),WRK(KMAT+NNBASX),1)
      CALL DSCAL(NNBASX,MUL1MM(3,I),WRK(KMAT+2*NNBASX),1)

      CALL DAXPY(NNBASX,FACM1,WRK(KMAT),1,TAO,1)
      CALL DAXPY(NNBASX,FACM1,WRK(KMAT+NNBASX),1,TAO,1)
      CALL DAXPY(NNBASX,FACM1,WRK(KMAT+2*NNBASX),1,TAO,1)

      EXCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
      EYCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
      EZCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

      ENSEL  = EXCOMP + EYCOMP + EZCOMP

C     Now the QM nuclear contribution. Note that we switch the sign here
C     by writing CORD(1,J) - MMCORD(1,I)

      ELOC     = 0.0D0
      DO 201 J = 1,NUCIND
         XDIS   = CORD(1,J) - MMCORD(1,I) 
         YDIS   = CORD(2,J) - MMCORD(2,I) 
         ZDIS   = CORD(3,J) - MMCORD(3,I) 
         DIST2  = XDIS**2+YDIS**2+ZDIS**2
         DIST   = SQRT(DIST2) 
         DIST3  = DIST2*DIST
         ELOC   = ELOC
     *           + CHARGE(J)*MUL1MM(1,I)*XDIS/DIST3
     *           + CHARGE(J)*MUL1MM(2,I)*YDIS/DIST3 
     *           + CHARGE(J)*MUL1MM(3,I)*ZDIS/DIST3
 201  CONTINUE
      ENSNUC = ELOC

      CALL QEXIT('DIPOLE_ITER')

      RETURN
      END
C******************************************************************************
C  /* Deck quadpole_iter */
      SUBROUTINE QUADPOLE_ITER(I,DCAO,ENSEL,ENSNUC,LOCDEB,
     &                       TAO,WRK,LWRK,IPRTMP)

#include "implicit.h"
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "gnrinf.h"
#include "orgcom.h"
#include "priunit.h"

      DIMENSION WRK(LWRK), DCAO(NNBASX), TAO(NNBASX)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0 )

      CALL QENTER('QUADPOLE_ITER')

      FAC1   =  1.0D0
      FACM1  = -1.0D0

      KMAT = 1
      KLAST = KMAT + 6*NNBASX
      LWRK2 = LWRK - KLAST + 1
      IF (LWRK2 .LT. 0) CALL ERRWRK('QUADPOLE_ITER',-KLAST,LWRK)

C       See if the quadrupole moment at this site is zero
      DNORM2 = MUL2MM(1,I)**2+MUL2MM(2,I)**2+MUL2MM(3,I)**2
     *         + MUL2MM(4,I)**2+MUL2MM(5,I)**2+MUL2MM(6,I)**2
      DNORM = SQRT(DNORM2)
      IF (ABS(DNORM) .LE. THRMM) THEN
         ENSEL = 0.0D0
         ENSNUC = 0.0D0
         CALL QEXIT('QUADPOLE_ITER')
         RETURN
      ENDIF

      CALL DZERO(WRK(KMAT),6*NNBASX)

      KPATOM = 0
      NOSIM = 6
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.
      DIPORG(1) = MMCORD(1,I)
      DIPORG(2) = MMCORD(2,I)
      DIPORG(3) = MMCORD(3,I)

      RUNQM3=.TRUE.
      CALL GET1IN(WRK(KMAT),'ELFGRDC',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
      RUNQM3=.FALSE.

      IF ( (IPRTMP.GT.15) .OR. (LOCDEB) ) THEN
         WRITE (LUPRI,'(/A)') ' E_xx_matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
            
         WRITE (LUPRI,'(/A)') ' E_xy matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

         WRITE (LUPRI,'(/A)') ' E_xz matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)

         WRITE (LUPRI,'(/A)') ' E_yy_matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT+3*NNBASX),NBAST,1,LUPRI)

         WRITE (LUPRI,'(/A)') ' E_yz_matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT+4*NNBASX),NBAST,1,LUPRI)

         WRITE (LUPRI,'(/A)') ' E_zz_matrix in QMMM_FCK:'
         CALL OUTPAK(WRK(KMAT+5*NNBASX),NBAST,1,LUPRI)
      END IF

      CALL DSCAL(NNBASX,MUL2MM(1,I),WRK(KMAT),1)
      CALL DSCAL(NNBASX,D2*MUL2MM(2,I),WRK(KMAT+NNBASX),1)
      CALL DSCAL(NNBASX,D2*MUL2MM(3,I),WRK(KMAT+2*NNBASX),1)
      CALL DSCAL(NNBASX,MUL2MM(4,I),WRK(KMAT+3*NNBASX),1)
      CALL DSCAL(NNBASX,D2*MUL2MM(5,I),WRK(KMAT+4*NNBASX),1)
      CALL DSCAL(NNBASX,MUL2MM(6,I),WRK(KMAT+5*NNBASX),1)

      FACS = 0.5D0
      CALL DSCAL(6*NNBASX,FACS,WRK(KMAT),1)
C
C     The integrals contains a factor of -1. Therefore daxpy with fac1
      CALL DAXPY(NNBASX,FAC1,WRK(KMAT),1,TAO(1),1)
      CALL DAXPY(NNBASX,FAC1,WRK(KMAT+NNBASX),1,TAO(1),1)
      CALL DAXPY(NNBASX,FAC1,WRK(KMAT+2*NNBASX),1,TAO(1),1)
      CALL DAXPY(NNBASX,FAC1,WRK(KMAT+3*NNBASX),1,TAO(1),1)
      CALL DAXPY(NNBASX,FAC1,WRK(KMAT+4*NNBASX),1,TAO(1),1)
      CALL DAXPY(NNBASX,FAC1,WRK(KMAT+5*NNBASX),1,TAO(1),1)

C     Contract with the density to get the expectation values.  The
C     factor of 1/2 in the Taylor expansion has been included.  Also,
C     the off-diagonal elements have been scaled by 2 in order to
C     include all contributions (the off -diagonal parts are related by
C     symmetry)

C     Since the integrals contains a factor of -1 no -DDOT here.

      EMU2XX=DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
      EMU2XY=DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
      EMU2XZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)
      EMU2YY=DDOT(NNBASX,DCAO,1,WRK(KMAT+3*NNBASX),1)
      EMU2YZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+4*NNBASX),1)
      EMU2ZZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+5*NNBASX),1)

      EQTOT = EMU2XX + EMU2XY + EMU2XZ + EMU2YY + EMU2YZ + EMU2ZZ
      ENSEL = EQTOT

C     Now the QM nuclear contribution

      ELOC     = 0.0D0
      DO 301 J = 1,NUCIND
         XDIS   = CORD(1,J) - MMCORD(1,I)  
         YDIS   = CORD(2,J) - MMCORD(2,I) 
         ZDIS   = CORD(3,J) - MMCORD(3,I)  
         DIST2  = XDIS**2+YDIS**2+ZDIS**2
         DIST   = SQRT(DIST2)
         DIST3  = DIST2*DIST
         DIST5  = DIST3*DIST2
C
         TXX    = (3.0D0*XDIS*XDIS - DIST2)/DIST5
         TXY    =  3.0D0*XDIS*YDIS/DIST5
         TXZ    =  3.0D0*XDIS*ZDIS/DIST5
         TYY    = (3.0D0*YDIS*YDIS - DIST2)/DIST5
         TYZ    =  3.0D0*YDIS*ZDIS/DIST5
         TZZ    = (3.0D0*ZDIS*ZDIS - DIST2)/DIST5

         ELOC   =   ELOC
     *           +   CHARGE(J)*MUL2MM(1,I)*TXX
     *           + 2*CHARGE(J)*MUL2MM(2,I)*TXY
     *           + 2*CHARGE(J)*MUL2MM(3,I)*TXZ
     *           +   CHARGE(J)*MUL2MM(4,I)*TYY
     *           + 2*CHARGE(J)*MUL2MM(5,I)*TYZ
     *           +   CHARGE(J)*MUL2MM(6,I)*TZZ
 301  CONTINUE
        
C     Remember the factor of 1/2 from the Taylor expansion
      ELOC   = 0.5D0*ELOC

      ENSNUC = ELOC

      CALL QEXIT('QUADPOLE_ITER')

      RETURN
      END

C******************************************************************************
C  /* Deck get_field */
      SUBROUTINE GET_FIELD(I,LRI,ELF,ELFEL,ELFNU,DCAO,
     &                         LOCDEB,WRK,LWRK)

      IMPLICIT NONE
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "gnrinf.h"
#include "orgcom.h"
#include "priunit.h"

      INTEGER LWRK, I, LRI
      DOUBLE PRECISION WRK, DCAO, ELF, ELFEL, ELFNU
      DIMENSION WRK(LWRK), DCAO(NNBASX)
      DIMENSION ELF(*), ELFEL(*), ELFNU(*)
      LOGICAL LSKIP, LOCDEB

      INTEGER KMAT, KEND, LWRK1
      DOUBLE PRECISION EXELCO, EYELCO, EZELCO, DDOT

      CALL QENTER('GET_FIELD')

      KMAT  = 1
      KEND  = KMAT + 3*NNBASX
      LWRK1 = LWRK - KEND + 1

C     Calculate field due to MM multipoles
      CALL CCMM_FMUL(ELF,LRI,I)

C     Add QM region contribution to the F vector

C     A) electronic contribution

      CALL DZERO(WRK(KMAT),3*NNBASX)

      LSKIP = .FALSE.

      CALL CCMM_EPSAO(WRK(KMAT),I,LSKIP,WRK(KEND),LWRK1)

      IF (LSKIP) THEN
         CALL QEXIT('GET_FIELD')
         RETURN
      END IF

      EXELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
      EYELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
      EZELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

      IF (SPLDIP) THEN
         ELFEL(LRI+0) = EXELCO
         ELFEL(LRI+1) = EYELCO
         ELFEL(LRI+2) = EZELCO
      ELSE
         ELF(LRI+0) = ELF(LRI+0) + EXELCO
         ELF(LRI+1) = ELF(LRI+1) + EYELCO
         ELF(LRI+2) = ELF(LRI+2) + EZELCO
      ENDIF

      IF (LOCDEB) THEN
         WRITE(LUPRI,*) 'electronic field:',EXELCO,EYELCO,EZELCO
      ENDIF
        
C     B) nuclear contribution
      IF (SPLDIP) THEN
         CALL CCMM_FNUC(ELFNU,LRI,I)
      ELSE  
         CALL CCMM_FNUC(ELF,LRI,I)
      END IF

      CALL QEXIT('GET_FIELD')

      RETURN
      END
C******************************************************************************
C  /* Deck get_pol_contr */
      SUBROUTINE GET_POL_CONTR(I,DINDMOM,EDALL,DCAO,TAO,
     &                         WRK,LWRK)

#include "implicit.h"
#include "mxcent.h"
#include "inforb.h"
#include "nuclei.h"
#include "qm3.h"
#include "qmmm.h"
#include "gnrinf.h"
#include "orgcom.h"
#include "priunit.h"

      DIMENSION WRK(LWRK), DCAO(NNBASX), TAO(NNBASX)
      DIMENSION DINDMOM(3),EDALL(6)
      LOGICAL LSKIP, EXCENT
      PARAMETER ( DMINV2 = -0.50D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0 )

      CALL QENTER('GET_POL_CONTR')

      FACM1 = -1.0D0

      KMAT = 1
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      LSKIP = .FALSE. 

      CALL CCMM_EPSAO(WRK(KMAT),I,LSKIP,WRK(KLAST),LWRK2)

      IF (LSKIP) THEN
         CALL QEXIT('GET_POL_CONTR')
         RETURN
      ENDIF

      CALL DZERO(EDALL,6)
      CALL DSCAL(NNBASX,DINDMOM(1),WRK(KMAT),1)
      CALL DSCAL(NNBASX,DINDMOM(2),WRK(KMAT+NNBASX),1)
      CALL DSCAL(NNBASX,DINDMOM(3),WRK(KMAT+2*NNBASX),1)

      CALL DAXPY(NNBASX,FACM1,WRK(KMAT),1,TAO,1)
      CALL DAXPY(NNBASX,FACM1,WRK(KMAT+NNBASX),1,TAO,1)
      CALL DAXPY(NNBASX,FACM1,WRK(KMAT+2*NNBASX),1,TAO,1)

C     Polarization contribution to the total energy 

C     A) Electronic contribution

      EXCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
      EYCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
      EZCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

      ET = 0.0D0
      ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
      EDALL(1) = ET

C     B) Nuclear contribution

      EFNUCX = 0.0D0
      EFNUCY = 0.0D0
      EFNUCZ = 0.0D0

      DO 510 J=1,NUCIND
         CALL GET_CHARGE_ELFLD(CHARGE(J),
     &                         CORD(1,J),CORD(2,J),CORD(3,J),
     &                         MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                         ELFLDX,ELFLDY,ELFLDZ)
         EFNUCX = EFNUCX + ELFLDX
         EFNUCY = EFNUCY + ELFLDY
         EFNUCZ = EFNUCZ + ELFLDZ
 510  CONTINUE

      IF (QMDAMP) THEN
         IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
         ENDIF
         IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIQM = 9.99D+99
            MHIT = 0
            DO 125 M=1,NUCIND
               DIQMC = (MMCORD(1,I)-CORD(1,M))**2 +
     &                 (MMCORD(2,I)-CORD(2,M))**2 +
     &                 (MMCORD(3,I)-CORD(3,M))**2
               IF (DIQMC .LE. DIQM) THEN
                  DIQM = DIQMC
                  MHIT = M
               ENDIF
 125        CONTINUE
         ELSE IF (IDAMP .EQ. 2) THEN
            DIQM = (MMCORD(1,I)-QMCOM(1))**2 +
     &             (MMCORD(2,I)-QMCOM(2))**2 +
     &             (MMCORD(3,I)-QMCOM(3))**2
         ENDIF
         DIQM = SQRT(DIQM)

         IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
               TEMPI =  D3I*(POLMM(1,I)+POLMM(4,I)+POLMM(6,I))
            ELSE IF (IPOLTP .EQ. 1) THEN
               IF (IPOLTP .EQ. 1) TEMPI = POLIMM(I)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304*DIQM/TEMP
            DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
         ELSE
            DFACT = (1-exp(-ADAMP*DIQM))**3
         ENDIF

         EFNUCX = EFNUCX*DFACT
         EFNUCY = EFNUCY*DFACT
         EFNUCZ = EFNUCZ*DFACT
      END IF

      EXCOMP = DINDMOM(1)*EFNUCX
      EYCOMP = DINDMOM(2)*EFNUCY
      EZCOMP = DINDMOM(3)*EFNUCZ

      ET = 0.0D0
      ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
      EDALL(2) = ET

C     C) Multipole contribution

      EF0MX = 0.0D0
      EF0MY = 0.0D0
      EF0MZ = 0.0D0
      EF1MX = 0.0D0
      EF1MY = 0.0D0
      EF1MZ = 0.0D0
      EF2MX = 0.0D0
      EF2MY = 0.0D0
      EF2MZ = 0.0D0
      EF3MX = 0.0D0
      EF3MY = 0.0D0
      EF3MZ = 0.0D0

C     Get electric fields due to permanent moments

      DO 520 J=1,MMCENT

         IF (J .EQ. I) GOTO 520

         EXCENT = .FALSE.
         IF (NEWEXC) THEN
           DO L = 1, NEXLST
             IF (EXLIST(1,I) .EQ. EXLIST(L,J)) EXCENT = .TRUE.
           ENDDO
         ELSE
           DO L = 1, NEXLST
             IF (EXLIST(L,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
           END DO
         ENDIF

         IF (.NOT. EXCENT) THEN

C     C.1  Point-charge contribution

            IF ( (NMULT .GE. 0) .AND. 
     &           (ABS(MUL0MM(J)) .GT. THRMM) ) THEN

               CALL GET_CHARGE_ELFLD(MUL0MM(J),
     &                        MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                        MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                        ELFLDX,ELFLDY,ELFLDZ)

               EF0MX = EF0MX + ELFLDX
               EF0MY = EF0MY + ELFLDY
               EF0MZ = EF0MZ + ELFLDZ
            ENDIF

C     C.2  Dipole contribution
            
            IF (NMULT .GE. 1) THEN

               CALL GET_DIPOLE_ELFLD(MUL1MM(1,J),MUL1MM(2,J),
     &                        MUL1MM(3,J),
     &                        MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                        MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                        ELFLDX,ELFLDY,ELFLDZ)

               EF1MX = EF1MX + ELFLDX
               EF1MY = EF1MY + ELFLDY
               EF1MZ = EF1MZ + ELFLDZ

            ENDIF

C     C.3  Quadrupole contribution

            IF (NMULT .GE. 2) THEN

               CALL GET_QUADRUPOLE_ELFLD(
     &                        MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     &                        MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     &                        MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                        MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                        ELFLDX,ELFLDY,ELFLDZ)

               EF2MX = EF2MX + ELFLDX
               EF2MY = EF2MY + ELFLDY
               EF2MZ = EF2MZ + ELFLDZ

            ENDIF

         ENDIF

 520  CONTINUE

C     Point-charge contribution

      IF (NMULT .GE. 0) THEN
            
         EXCOMP = DINDMOM(1)*EF0MX
         EYCOMP = DINDMOM(2)*EF0MY
         EZCOMP = DINDMOM(3)*EF0MZ

         ET = 0.0D0
         ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
         EDALL(3) = ET

      ENDIF

C     Dipole contribution

      IF (NMULT .GE. 1) THEN

         EXCOMP = DINDMOM(1)*EF1MX
         EYCOMP = DINDMOM(2)*EF1MY
         EZCOMP = DINDMOM(3)*EF1MZ

         ET = 0.0D0
         ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
         EDALL(4) = ET

      ENDIF

C     Quadrupole contribution

      IF (NMULT .GE. 2) THEN

         EXCOMP = DINDMOM(1)*EF2MX
         EYCOMP = DINDMOM(2)*EF2MY
         EZCOMP = DINDMOM(3)*EF2MZ

         ET = 0.0D0
         ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
         EDALL(5) = ET

      ENDIF

      CALL QEXIT('GET_POL_CONTR')

      RETURN
      END
C******************************************************************************
C  /* Deck get_my */
      SUBROUTINE GET_MY(I,J,DIP,MY)

C Input: I,J,DIP
C Output: MY
C Get the polarizability tensor MY at site I due to polarizability
C tensor DIP at site J.
      IMPLICIT NONE
#include "mxcent.h"
#include "qmmm.h"

      LOGICAL EXCENT
      INTEGER I, J, K, L, NLDIM, NTOTI
      DOUBLE PRECISION AMAT,TTENS,ATMAT,DIP,MY
      DOUBLE PRECISION TEMPJ, TEMP, SCREEN, TEMPI,DIST,DIST2,DIST3,DIST5
      DOUBLE PRECISION ELEM,FEIJ,FTIJ, D0, D1, D3I, D6I
      DIMENSION AMAT(3,3),TTENS(3,3),ATMAT(3,3),DIP(3),MY(3)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('GET_MY')

      CALL DZERO(MY,3)
      IF (J .NE. I) THEN
        EXCENT = .FALSE.
        IF (NEWEXC) THEN
          DO L = 1, NEXLST
            IF (EXLIST(1,I) .EQ. EXLIST(L,J)) EXCENT = .TRUE.
          ENDDO
        ELSE
          DO L = 1, NEXLST
            IF (EXLIST(L,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
          END DO
        ENDIF

        IF (.NOT. EXCENT) THEN
C     Get the polarizability tensor for this site
          DO K=1,3
            DO L=1,3
              AMAT(K,L)  = 0.0D0
            ENDDO
          ENDDO

          IF (IPOLTP .EQ. 1)  THEN
            DO L=1,3
              AMAT(L,L) = POLIMM(I)
            ENDDO
          ELSE IF (IPOLTP .EQ. 2)  THEN
            AMAT(1,1) = POLMM(1,I)
            AMAT(1,2) = POLMM(2,I)
            AMAT(1,3) = POLMM(3,I)
            AMAT(2,1) = POLMM(2,I)
            AMAT(2,2) = POLMM(4,I)
            AMAT(2,3) = POLMM(5,I)
            AMAT(3,1) = POLMM(3,I)
            AMAT(3,2) = POLMM(5,I)
            AMAT(3,3) = POLMM(6,I)
          ENDIF

C     Now calculate the T tensor for these sites
          DIST2 = 0.0D0
          DO K=1,3
            DIST2 = DIST2 + (MMCORD(K,I) - MMCORD(K,J))**2
          ENDDO
          DIST  = SQRT(DIST2)
          DIST3 = DIST**3
          DIST5 = DIST**5

          DO K=1,3
            DO L=1,3

C     Include damping in the exponential form
C     JPC A 102 (1998) 2399
              IF (MMDAMP) THEN
                IF (IPOLTP .EQ. 1) THEN
                  TEMPI = POLIMM(I)
                  TEMPJ = POLIMM(J)
                ELSE IF (IPOLTP .EQ. 2) THEN
                  TEMPI =  D3I*(POLMM(1,I)+POLMM(4,I)+POLMM(6,I))
                  TEMPJ =  D3I*(POLMM(1,J)+POLMM(4,J)+POLMM(6,J))
                ENDIF
                TEMP = (TEMPI*TEMPJ)**(D6I)
                SCREEN = 2.1304*DIST/TEMP
                FEIJ = 1.0D0-(1.0D0+SCREEN+0.5D0*SCREEN**2)
     &                     *EXP(-SCREEN)
                FTIJ = FEIJ - (1.0D0/6.0D0*SCREEN**3)
     &                     *EXP(-SCREEN)
              ELSE
                FEIJ = D1
                FTIJ = D1
              ENDIF

              ELEM = FTIJ*3*(MMCORD(K,I) - MMCORD(K,J))*
     &                          (MMCORD(L,I) - MMCORD(L,J))
              ELEM = ELEM/DIST5
              IF (K .EQ. L) ELEM = ELEM - (FEIJ*1.0/DIST3)
              TTENS(K,L) = ELEM
            ENDDO
          ENDDO

C     calculate alpha*T
          CALL DGEMM('N','N',3,3,3,1.D0,AMAT,3,
     &                   TTENS,3,0.D0,ATMAT,3)

          NLDIM = 3
          NTOTI = MAX(NLDIM,1)
          CALL DGEMV('N',NLDIM,NLDIM,D1,ATMAT,NTOTI,DIP,1,D0,MY,1)
        ENDIF
      ENDIF

      CALL QEXIT('GET_MY')

      RETURN
      END
C******************************************************************************
C  /* Deck qmmmtimes */
      SUBROUTINE QMMMTIMES(WORD)

      IMPLICIT NONE

#include "priunit.h"
#include "mxcent.h"
#include "qmmm.h"
#include "mmtimes.h"
      CHARACTER*(*) WORD
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO = 0.0D0)
      CALL QENTER('QMMMTIMES')

      IF (.NOT. MMTIME) THEN
        CALL QEXIT('QMMMTIMES')
        RETURN
      ENDIF

      WRITE(LUPRI,*) '  - QM/MM times:'
      IF (WORD .EQ. 'SIRIUS') THEN
        WRITE(LUPRI,1) 'QMMMFCK      ',TMMFCK
        WRITE(LUPRI,1) 'QMMM MULPOLES',TMMMULPOL
        WRITE(LUPRI,1) 'QMMM_POLARI  ',TMMPOL
        IF (MMITER) THEN
          WRITE(LUPRI,*) '    - MMITER times:'
          WRITE(LUPRI,2) 'GET_IND_DIPOLES_2',TMMGID2
          WRITE(LUPRI,2) 'GET_FIELD        ',TMMPOL2
          WRITE(LUPRI,2) 'F2QMMM           ',TF2QMMM
          WRITE(LUPRI,2) 'the iteration    ',TMMITER
          TMMGID2 = ZERO
          TMMPOL2 = ZERO
          TF2QMMM = ZERO
          TMMITER = ZERO
        ENDIF
      ELSEIF (WORD .EQ. 'RESPONSE') THEN
        WRITE(LUPRI,1) 'QMMMRSP',TMMRSP
        WRITE(LUPRI,2) 'QMMMLNO      ',TMMLNO
        WRITE(LUPRI,3) 'QMMMLNO0     ',TMMLNO0
        WRITE(LUPRI,3) 'QMMMLNO1     ',TMMLNO1
        WRITE(LUPRI,3) 'QMMMLNO2     ',TMMLNO2
        WRITE(LUPRI,3) 'QMMMLNO3     ',TMMLNO3
        WRITE(LUPRI,3) 'QMMMLNO4     ',TMMLNO4
        WRITE(LUPRI,2) 'QMMMQRO      ',TMMQRO
        WRITE(LUPRI,3) 'QMMMQRO0     ',TMMQRO0
        WRITE(LUPRI,3) 'QMMMQRO1     ',TMMQRO1
        WRITE(LUPRI,3) 'QMMMQRO2     ',TMMQRO2
        WRITE(LUPRI,3) 'QMMMQRO3     ',TMMQRO3
        WRITE(LUPRI,3) 'QMMMQRO4     ',TMMQRO4
        WRITE(LUPRI,2) 'QMMMCRO      ',TMMCRO
        IF (MMITER) THEN
          WRITE(LUPRI,2) 'F2QMMM       ',TF2QMMM
          WRITE(LUPRI,2) 'the iteration',TMMITER
        ENDIF
      ELSEIF (WORD .EQ. 'ABACUS') THEN
        WRITE(LUPRI,1) 'QMMMFIRST',TMMFIRST
        WRITE(LUPRI,1) 'QMMMB2   ',TMMB2
      ENDIF

 1    FORMAT( '     - total time used in ',A,': ',F10.2,' seconds')
 2    FORMAT( '       - total time used in ',A,': ',F10.2,' seconds')
 3    FORMAT( '         - total time used in ',A,': ',F10.2,' seconds')

      CALL QEXIT('QMMMTIMES')

      RETURN
      END
C
#if defined(VAR_MPI)
C*****************************************************
C Parallel routines for QM/MM SIRIUS
C Arnfinn, Odense/Tromso Oct. 2009 - Oct. 2010
C As little as possible of calculations are done here;
C Mostly calling routines shared by the serial code.
C*****************************************************
C  /* Deck parqmmm_m */
      SUBROUTINE PARQMMM_M(DCAO,TAO,ESOLT,LOCDEB,WRK,LWRK,
     &                         IPRINT)

#include "implicit.h"
!  mxcoor in nuclei.h
#include "mxcent.h"
!  nnbasx, icmo, nbast, 
#include "inforb.h"
!  nctot, cord, charge, nucind, nucdep
#include "nuclei.h"
! luprop
#include "inftap.h"
! npatom, ipatom
#include "cbiher.h"
! qmcom, isytp, qmdamp
#include "qm3.h"
! mmcent, mul0mm, mul1mm etc, rcutmm, delfld, nmult, nexlst, exlist
! nnzal (updates), spldip, zeroal (updates?), idamp, ipoltp, 
! From potfile: mmcent, nmult, ipoltp, nexlst, neleme
!               exlist
#include "qmmm.h"
C ----
#include "maxorb.h"
! MXSHEL
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
C#include "cbiher.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK), TAO(NNBASX), DCAO(NNBASX)
      LOGICAL LOCDEB

      CALL QENTER('PARQMMM_M')

      KNSNUC   = 1
      KNSNUC2  = KNSNUC  + MMCENT
      KNSEL    = KNSNUC2 + MMCENT
      KNSEL2   = KNSEL   + MMCENT
      KTAO     = KNSEL2  + MMCENT
      KTAO2    = KTAO    + NNBASX
      KLAST    = KTAO2   + NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('PARQMMM_M',-KLAST,LWRK)

      CALL DZERO(WRK(KNSNUC),MMCENT)
      CALL DZERO(WRK(KNSNUC2),MMCENT)
      CALL DZERO(WRK(KNSEL),MMCENT)
      CALL DZERO(WRK(KNSEL2),MMCENT)
      CALL DZERO(WRK(KTAO),NNBASX)
      CALL DZERO(WRK(KTAO2),NNBASX)
      ECHTMP = 0.0D0
      EDITMP = 0.0D0
      EQUTMP = 0.0D0

C     Wake up slaves (Rock and roll all nite)

      IPRTYP = PARQMMM__WORK
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

C     Send data to slaves (Lick it up)

      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)
      IF (NMULT .GE. 0) CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      IF (NMULT .GE. 1) CALL MPIXBCAST(MUL1MM,3*MMCENT,'DOUBLE',MASTER)
      IF (NMULT .GE. 2) CALL MPIXBCAST(MUL2MM,6*MMCENT,'DOUBLE',MASTER)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(DCAO,NNBASX,'DOUBLE',MASTER)

      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMITER,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMPROP,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMDIIS,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(LOCDEB,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMDAMP,1,'LOGICAL',MASTER)

C     The loop (Shock me)
      LNUM = 0
      DO L = 1,MMCENT
        LNUM = LNUM + 1
        IWHO = -1
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(LNUM,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Send end message to all slaves (Rock bottom)
      LEND = -1
      DO ISLAVE = 1,NODTOT
        IWHO = -1
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves (Great expectations) 

      CALL MPI_REDUCE(WRK(KNSNUC2),WRK(KNSNUC),MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KNSEL2),WRK(KNSEL),MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KTAO2),WRK(KTAO),NNBASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      IF (NMULT .GE. 0) CALL MPI_REDUCE(ECHTMP,ECHART,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      IF (NMULT .GE. 1) CALL MPI_REDUCE(EDITMP,EDIPT,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      IF (NMULT .GE. 2) CALL MPI_REDUCE(EQUTMP,EQUADT,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      CALL DAXPY(NNBASX,1.0D0,WRK(KTAO),1,TAO(1),1)
      
      ECHCH  = 0.0D0
      EXPNST = 0.0D0

      DO I = 1, MMCENT
         ECHCH  = ECHCH  + WRK(KNSNUC + I - 1)
         EXPNST = EXPNST + WRK(KNSEL  + I - 1)
      END DO

      ENUMUL = ECHCH
      ESOLT = ECHART + EDIPT + EQUADT

      CALL QEXIT('PARQMMM_M')

      RETURN
      END
C******************************************************************************
C  /* Deck parqmmm_s */
      SUBROUTINE PARQMMM_S(WRK,LWRK,IPRTMP)

#include "implicit.h"
!  mxcoor in nuclei.h
#include "mxcent.h"
!  nnbasx, icmo, nbast, 
#include "inforb.h"
!  nctot, cord, charge, nucind, nucdep
#include "nuclei.h"
! luprop
#include "inftap.h"
! npatom, ipatom
#include "cbiher.h"
! qmcom, isytp, qmdamp
#include "qm3.h"
! mmcent, mul0mm, mul1mm etc, rcutmm, delfld, nmult, nexlst, exlist
! nnzal (updates), spldip, zeroal (updates?), idamp, ipoltp, 
! From potfile: mmcent, nmult, ipoltp, nexlst, neleme
!               exlist
#include "qmmm.h"
#include "maxorb.h"
! MXSHEL
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
! qmmm
#include "gnrinf.h"
! diporg
#include "orgcom.h"
#include "priunit.h"

      DIMENSION WRK(LWRK)
      LOGICAL LOCDEB

      CALL QENTER('PARQMMM_S')

      QMMM = .TRUE.

C     Receiving data from master (I was made for lovin' you)
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)
      IF (NMULT .GE. 0) CALL MPIXBCAST(MUL0MM,MMCENT,'DOUBLE',MASTER)
      IF (NMULT .GE. 1) CALL MPIXBCAST(MUL1MM,3*MMCENT,'DOUBLE',MASTER)
      IF (NMULT .GE. 2) CALL MPIXBCAST(MUL2MM,6*MMCENT,'DOUBLE',MASTER)

      KNSNUC = 1
      KNSEL  = KNSNUC + MMCENT
      KTAO   = KNSEL  + MMCENT
      KDCAO  = KTAO   + NNBASX
      KLAST  = KDCAO  + NNBASX
      LWRK2  = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('PARQMMM_S',-KLAST,LWRK)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)
      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(ICMO,8,'INTEGER',MASTER)

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(WRK(KDCAO),NNBASX,'DOUBLE',MASTER)

      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(MMITER,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMPROP,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMDIIS,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(LOCDEB,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(MMDAMP,1,'LOGICAL',MASTER)

C     Do the work (I love it load)

      CALL DZERO(WRK(KNSEL),MMCENT)
      CALL DZERO(WRK(KNSNUC),MMCENT)
      CALL DZERO(WRK(KTAO),NNBASX)

      ECHART = 0.0D0
      EDIPT  = 0.0D0
      EQUADT = 0.0D0

 1    CONTINUE

      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(I,1,'INTEGER',MASTER,MPTAG2)

      IF (I.GT.0) THEN

         DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *           (MMCORD(2,I)-QMCOM(2))**2 +
     *           (MMCORD(3,I)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
           GOTO 1
         ENDIF

C------------------------------------------------- 
C        Charge contributions:
C------------------------------------------------- 
         CALL CHARGE_ITER(I,WRK(KDCAO),ENSEL,ENSNUC,LOCDEB,
     *                 WRK(KTAO),WRK(KLAST),LWRK2,IPRTMP)
         WRK(KNSEL+I-1)  = WRK(KNSEL+I-1)  + ENSEL
         WRK(KNSNUC+I-1) = WRK(KNSNUC+I-1) + ENSNUC
         ECHART = ECHART + ENSEL + ENSNUC
         IF (NMULT .LT. 1) GOTO 1

C------------------------------------------------- 
C        Dipole contributions:
C------------------------------------------------- 
         CALL DIPOLE_ITER(I,WRK(KDCAO),ENSEL,ENSNUC,LOCDEB,
     *                    WRK(KTAO),WRK(KLAST),LWRK2,IPRTMP)
         WRK(KNSEL+I-1)  = WRK(KNSEL+I-1)  + ENSEL
         WRK(KNSNUC+I-1) = WRK(KNSNUC+I-1) + ENSNUC
         EDIPT = EDIPT + ENSEL + ENSNUC
         IF (NMULT .LT. 2) GOTO 1

C------------------------------------------------- 
C        Quadrupole contributions:
C------------------------------------------------- 
         CALL QUADPOLE_ITER(I,WRK(KDCAO),ENSEL,ENSNUC,LOCDEB,
     *                    WRK(KTAO),WRK(KLAST),LWRK2,IPRTMP)
         WRK(KNSEL+I-1)  = WRK(KNSEL+I-1)  + ENSEL
         WRK(KNSNUC+I-1) = WRK(KNSNUC+I-1) + ENSNUC
         EQUADT = EQUADT + ENSEL + ENSNUC
         GOTO 1
      ENDIF

C     Send data to master (Do you love me?)
      CALL MPI_REDUCE(WRK(KNSNUC),MPI_IN_PLACE,MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KNSEL),MPI_IN_PLACE,MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KTAO),MPI_IN_PLACE,NNBASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      IF (NMULT .GE. 0) CALL MPI_REDUCE(ECHART,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      IF (NMULT .GE. 1) CALL MPI_REDUCE(EDIPT,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      IF (NMULT .GE. 2) CALL MPI_REDUCE(EQUADT,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      DIPORG(1) = OBKPX 
      DIPORG(2) = OBKPY 
      DIPORG(3) = OBKPZ 

      CALL QEXIT('PARQMMM_S')
      RETURN
      END
C******************************************************************************
C  /* Deck mm_field_m1 */
      SUBROUTINE MM_FIELD_M1(DCAO,ELF,POLDIM,WRK,LWRK,IPRINT)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
C defined parallel calculation types  
#include "iprtyp.h"

      INTEGER POLDIM
      DIMENSION WRK(LWRK), ELF(*)

      CALL QENTER('MM_FIELD_M1')

      KELF  = 1
      KEND  = KELF + 3*POLDIM
      IF (SPLDIP) THEN
         KELFEL  = KEND
         KELFNU  = KELFEL + 3*POLDIM
         KEND    = KELFNU + 3*POLDIM
      ENDIF
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('MM_FIELD_M1',-KEND,LWRK)

C     Beginning of parallel section

      IPRTYP = MM_FIELD_1_WORK

C     Wake up slaves

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

C     Send data to slaves

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(POLDIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(SPLDIP,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(CONMAT,1,'LOGICAL',MASTER)

      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NEXLST,1,'INTEGER',MASTER)
      DO N=1, NEXLST
         CALL MPIXBCAST(EXLIST(N,1:MMCENT),MMCENT,'INTEGER',MASTER)
      ENDDO
C      CALL MPIXBCAST(EXLIST,NEXLST*MMCENT,'INTEGER',MASTER)

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,MXCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF

C      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(DELFLD,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(DCAO,NNBASX,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZEROAL,MMCENT,'INTEGER',MASTER)

C     Start parallelized loop
      LRI = 1
      DO 100 L = 1,MMCENT
        IWHO = -1
        IF (ZEROAL(L) .EQ. -1) GOTO 100
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(L,1,'INTEGER',NWHO,MPTAG2)
        CALL MPIXSEND(LRI,1,'INTEGER',NWHO,MPTAG2)
        LRI = LRI + 3
 100  CONTINUE

C     Send end message to all slaves

      LEND = -1
      DO ISLAVE = 1,NODTOT
        IWHO = -1
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
        CALL MPIXSEND(LRI,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves

      CALL DZERO(WRK(KELF),3*POLDIM)
      CALL MPI_REDUCE(WRK(KELF),ELF,3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      IF (SPLDIP) THEN 
         CALL DZERO(WRK(KELFEL),3*POLDIM)
         CALL MPI_REDUCE(WRK(KELFEL),ELF(3*POLDIM+1),3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

         CALL DZERO(WRK(KELFNU),3*POLDIM)
         CALL MPI_REDUCE(WRK(KELFNU),ELF(6*POLDIM+1),3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      ENDIF

      CALL QEXIT('MM_FIELD_M1')

      RETURN
      END
C******************************************************************************
C  /* Deck mm_field_s1 */
      SUBROUTINE MM_FIELD_S1(WRK,LWRK,IPRTMP)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"

      DIMENSION WRK(LWRK)
      INTEGER POLDIM

      CALL QENTER('MM_FIELD_S1')

      QMMM = .TRUE.
      IPQMMM = IPRTMP

C     Receiving data from master

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(POLDIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(SPLDIP,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(CONMAT,1,'LOGICAL',MASTER)

      KELF    = 1
      IF (SPLDIP) THEN
         KELFEL = KELF   + 3*POLDIM
         KELFNU = KELFEL + 3*POLDIM
         KDCAO  = KELFNU + 3*POLDIM
      ELSE
         KDCAO  = KELF   + 3*POLDIM
      ENDIF
      KMAT    = KDCAO   + NNBASX
      KLAST   = KMAT    + 3*NNBASX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MM_FIELD_S1',-KLAST,LWRK)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)
C
      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NEXLST,1,'INTEGER',MASTER)
      DO N=1, NEXLST
         CALL MPIXBCAST(EXLIST(N,1:MMCENT),MMCENT,'INTEGER',MASTER)
      ENDDO
C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,MXCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(DELFLD,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KDCAO),NNBASX,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZEROAL,MMCENT,'INTEGER',MASTER)

C     Do the work

      CALL DZERO(WRK(KELF),3*POLDIM)
      IF (SPLDIP) THEN
         CALL DZERO(WRK(KELFEL),3*POLDIM)
         CALL DZERO(WRK(KELFNU),3*POLDIM)
      ENDIF

 200  CONTINUE

      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(I,1,'INTEGER',MASTER,MPTAG2)
      CALL MPIXRECV(LRI,1,'INTEGER',MASTER,MPTAG2)

      IF (I.GT.0) THEN
         CALL GET_FIELD(I,LRI,WRK(KELF),WRK(KELFEL),WRK(KELFNU),
     &                       WRK(KDCAO),.FALSE.,WRK(KLAST),LWRK2)
         GO TO 200
      ENDIF

      CALL MPI_REDUCE(WRK(KELF),MPI_IN_PLACE,3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      IF (SPLDIP) THEN
         CALL MPI_REDUCE(WRK(KELFEL),MPI_IN_PLACE,3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
         CALL MPI_REDUCE(WRK(KELFNU),MPI_IN_PLACE,3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      ENDIF

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MM_FIELD_S1')

      RETURN
      END
C
C******************************************************************************
C  /* Deck mm_field_m2 */
      SUBROUTINE MM_FIELD_M2(DCAO,ELF,POLDIM,WRK,LWRK,IPRINT)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
C defined parallel calculation types  
#include "iprtyp.h"

      INTEGER POLDIM
      DIMENSION WRK(LWRK), ELF(*)

      CALL QENTER('MM_FIELD_M2')

      KELF  = 1
      KEND  = KELF + 3*POLDIM
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('MM_FIELD_M2',-KEND,LWRK)

C     Beginning of parallel section

      IPRTYP = MM_FIELD_2_WORK

C     Wake up slaves

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

C     Send data to slaves

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(POLDIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NEXLST,1,'INTEGER',MASTER)
      DO N=1, NEXLST
         CALL MPIXBCAST(EXLIST(N,1:MMCENT),MMCENT,'INTEGER',MASTER)
      ENDDO

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,MXCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF
C     <-

      CALL MPIXBCAST(DELFLD,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(DCAO,NNBASX,'DOUBLE',MASTER)

      LRI = 1   ! important should be one due to the indexing used !

C     Start parallelized loop
      DO 100 L = 1,MMCENT
        IWHO = -1
        IF (ZEROAL(L) .EQ. -1) GOTO 100
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(L,1,'INTEGER',NWHO,MPTAG2)
        CALL MPIXSEND(LRI,1,'INTEGER',NWHO,MPTAG2)
        LRI = LRI + 3
 100  CONTINUE

C     Send end message to all slaves

      LEND = -1
      DO ISLAVE = 1,NODTOT
        IWHO = -1
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
        CALL MPIXSEND(LRI,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves

      CALL DZERO(WRK(KELF),3*POLDIM)
      CALL MPI_REDUCE(WRK(KELF),ELF,3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      CALL QEXIT('MM_FIELD_M2')

      RETURN
      END
C******************************************************************************
C  /* Deck mm_field_s2 */
      SUBROUTINE MM_FIELD_S2(WRK,LWRK,IPRTMP)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"

      DIMENSION WRK(LWRK)
      INTEGER POLDIM

      CALL QENTER('MM_FIELD_S2')

      QMMM = .TRUE.
      SPLDIP = .FALSE.          ! Not implemented for iterative QMMM
      IPQMMM = IPRTMP

C     Receiving data from master

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(POLDIM,1,'INTEGER',MASTER)

      KELF    = 1
      KDCAO   = KELF    + 3*POLDIM
      KMAT    = KDCAO   + NNBASX
      KLAST   = KMAT    + 3*NNBASX
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MM_FIELD_S2',-KLAST,LWRK)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)
C
      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NEXLST,1,'INTEGER',MASTER)
      DO N=1, NEXLST
         CALL MPIXBCAST(EXLIST(N,1:MMCENT),MMCENT,'INTEGER',MASTER)
      ENDDO

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,MXCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF
C     <-

      CALL MPIXBCAST(DELFLD,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KDCAO),NNBASX,'DOUBLE',MASTER)

C     Do the work

      CALL DZERO(WRK(KELF),3*POLDIM)

 200  CONTINUE

      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(I,1,'INTEGER',MASTER,MPTAG2)
      CALL MPIXRECV(LRI,1,'INTEGER',MASTER,MPTAG2)

      IF (I.GT.0) THEN
         CALL GET_FIELD(I,LRI,WRK(KELF),WRK(KLAST),WRK(KLAST),
     *                  WRK(KDCAO),.FALSE.,WRK(KLAST),LWRK2)
         GOTO 200
      ENDIF

      CALL MPI_REDUCE(WRK(KELF),MPI_IN_PLACE,3*POLDIM,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MM_FIELD_S2')

      RETURN
      END
C******************************************************************************
C  /* Deck mm_polar_contr_m */
      SUBROUTINE MM_POLAR_CONTR_M(DCAO,TAO,CINDMOM,WRK,LWRK,IPRINT)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"
C defined parallel calculation types  
#include "iprtyp.h"

      DIMENSION WRK(LWRK), TAO(NNBASX), CINDMOM(*)

      CALL QENTER('MM_POLAR_CONTR_M')

      KTAO    = 1
      KTAO2   = KTAO   + NNBASX
      KREC    = KTAO2  + NNBASX
      KWRK2   = KREC   + 6
      LWRK2   = LWRK    - KWRK2 + 1

      IF (LWRK2 .LT. 0) THEN
         CALL ERRWRK('MM_POLAR_CONTR_M',-KWRK2,LWRK)
      ENDIF

      EDELD  = 0.0D0            ! For interaction with electronic density
      EDNUC  = 0.0D0            ! For interaction with QM nuclei
      ED0MOM = 0.0D0            ! For interaction with point-charges
      ED1MOM = 0.0D0            ! For interaction with permanent dipoles
      ED2MOM = 0.0D0            ! For interaction with quadrupoles
      EDMULT = 0.0D0            ! For interaction with permanent multipoles

C     Beginning of parallel section

      IPRTYP = MM_POLAR_CONTR_WORK

C     Wake up slaves

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

C     Send data to slaves

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)

      CALL MPIXBCAST(CINDMOM,3*NNZAL,'DOUBLE',MASTER)

      CALL MPIXBCAST(NEXLST,1,'INTEGER',MASTER)
      DO II = 1,NEXLST
         CALL MPIXBCAST(EXLIST(II,1:MMCENT),MMCENT,'INTEGER',MASTER)
      ENDDO

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,MXCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(DELFLD,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(DCAO,NNBASX,'DOUBLE',MASTER)

      IINIM = 0   ! important should be zero due to the indexing used !

C     Start parallelized loop
      DO 100 L = 1,MMCENT
        IWHO = -1
        IF (ZEROAL(L) .EQ. -1) GOTO 100
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(L,1,'INTEGER',NWHO,MPTAG2)
        CALL MPIXSEND(IINIM,1,'INTEGER',NWHO,MPTAG2)
        IINIM = IINIM + 3
 100  CONTINUE

C     Send end message to all slaves

      LEND = -1
      DO ISLAVE = 1,NODTOT
        IWHO = -1
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
        CALL MPIXSEND(IINIM,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves

      CALL DZERO(WRK(KTAO),NNBASX)
      CALL DZERO(WRK(KTAO2),NNBASX)
      CALL DZERO(WRK(KREC),6)
      CALL MPI_REDUCE(WRK(KTAO2),WRK(KTAO),NNBASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      CALL DAXPY(NNBASX,1.0D0,WRK(KTAO),1,TAO(1),1)

      CALL MPI_REDUCE(WRK(KREC+0),EDELD,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KREC+1),EDNUC,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KREC+2),ED0MOM,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KREC+3),ED1MOM,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KREC+4),ED2MOM,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      EDMULT = ED0MOM + ED1MOM + ED2MOM

      CALL QEXIT('MM_POLAR_CONTR_M')

      RETURN
      END
C******************************************************************************
C  /* Deck mm_polar_contr_s */
      SUBROUTINE MM_POLAR_CONTR_S(WRK,LWRK,IPRTMP)

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "iratdef.h"
#include "maxash.h"
#include "maxorb.h"

#include "qmmm.h"
#include "qm3.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "scbrhf.h"
#include "maxaqn.h"
#include "symmet.h"
#include "orgcom.h"
#include "infinp.h"
#include "nuclei.h"
#include "codata.h"
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
#include "cbiher.h"
#include "gnrinf.h"

      DIMENSION WRK(LWRK)

      CALL QENTER('MM_POLAR_CONTR_S')

      QMMM = .TRUE.
      IPQMMM = IPRTMP

C     Receiving data from master

      CALL MPIXBCAST(NNZAL,1,'INTEGER',MASTER)

      KDCAO   = 1
      KTAO    = KDCAO   + NNBASX
      KMAT    = KTAO    + NNBASX
      KINDMOM = KMAT    + 3*NNBASX
      KEDALL  = KINDMOM + 3*NNZAL
      KLAST   = KEDALL + 6
      LWRK2   = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('MM_POLAR_CONTR_S',-KLAST,LWRK)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      CALL MPIXBCAST(NMULT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KINDMOM),3*NNZAL,'DOUBLE',MASTER)

      CALL MPIXBCAST(NEXLST,1,'INTEGER',MASTER)
      DO II = 1,NEXLST
         CALL MPIXBCAST(EXLIST(II,1:MMCENT),MMCENT,'INTEGER',MASTER)
      ENDDO

C     Damping
      CALL MPIXBCAST(QMDAMP,1,'LOGICAL',MASTER)
      IF (QMDAMP) THEN
         CALL MPIXBCAST(IDAMP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NQMNUC,1,'INTEGER',MASTER)
         CALL MPIXBCAST(QMPOL,MXCENT,'DOUBLE',MASTER)
         CALL MPIXBCAST(ADAMP,1,'DOUBLE',MASTER)
         CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF

      CALL MPIXBCAST(DELFLD,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KDCAO),NNBASX,'DOUBLE',MASTER)

C     Compute polarization contributions to the Fock/KS matrix and
C     total solvation energy

      CALL DZERO(WRK(KTAO),NNBASX)

C     Compute polarization contributions to the Fock/KS matrix and
C     total solvation energy

      EDELD  = 0.0D0            ! For interaction with electronic density
      EDNUC  = 0.0D0            ! For interaction with QM nuclei
      ED0MOM = 0.0D0            ! For interaction with point-charges
      ED1MOM = 0.0D0            ! For interaction with permanent dipoles
      ED2MOM = 0.0D0            ! For interaction with quadrupoles

 20   CONTINUE

      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(I,1,'INTEGER',MASTER,MPTAG2)
      CALL MPIXRECV(IINIM,1,'INTEGER',MASTER,MPTAG2)

      IF (I.GT.0) THEN
         CALL GET_POL_CONTR(I,WRK(KINDMOM+IINIM),WRK(KEDALL),
     &                      WRK(KDCAO),WRK(KTAO),WRK(KLAST),LWRK2)
         EDELD  = EDELD  + WRK(KEDALL)
         EDNUC  = EDNUC  + WRK(KEDALL + 1)
         ED0MOM = ED0MOM + WRK(KEDALL + 2)
         ED1MOM = ED1MOM + WRK(KEDALL + 3)
         ED2MOM = ED2MOM + WRK(KEDALL + 4)
         GOTO 20
      ENDIF


      CALL MPI_REDUCE(WRK(KTAO),MPI_IN_PLACE,NNBASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      CALL MPI_REDUCE(EDELD,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(EDNUC,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(ED0MOM,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(ED1MOM,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(ED2MOM,MPI_IN_PLACE,1, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL QEXIT('MM_POLAR_CONTR_S')

      RETURN
      END
C******************************************************************************
C  /* Deck mmiter_inddip_m */
      SUBROUTINE MMITER_INDDIP_M(POLDIM,INDP1,INDMOM,VEC,INDDIA,
     *                    WRK,LWRK,LOCDEB,DIPCON,LM)

#include "implicit.h"
C      IMPLICIT NONE

#include "priunit.h"
#include "mxcent.h"
#include "qmmm.h"
#include "maxorb.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif
C defined parallel calculation types  
#include "iprtyp.h"

      INTEGER POLDIM, POLARRAY
      DIMENSION POLARRAY(POLDIM)

      LOGICAL LOCDEB,DIPCON

      DOUBLE PRECISION INDMOM,INDDIA,INDP1
      DIMENSION INDMOM(3*POLDIM),VEC(MXMMIT+3*POLDIM), INDDIA(3*POLDIM)
      DIMENSION WRK(LWRK),INDP1(3*POLDIM)

      DOUBLE PRECISION TERROR,TDIFF,TMAX
      DOUBLE PRECISION DIP,MY
      DIMENSION DIP(3),MY(3)

      CALL QENTER('MMITER_INDDIP_M')

      DIPCON = .FALSE.

      THRESL = THMMIT
      NDIM = 3*POLDIM

C     Make a vector of pol sites
      L = 0
      DO 1 I=1,MMCENT
        IF (ZEROAL(I) .EQ. -1) GOTO 1
        L = L + 1
        POLARRAY(L) = I
 1    CONTINUE

C     Beginning of parallel section 

      IPRTYP = MMITER_INDDIP_WORK

C     Wake up slaves

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

      CALL MPIXBCAST(POLDIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NODTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(POLARRAY,POLDIM,'INTEGER',MASTER)
      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      IF (IPOLTP .EQ. 1) THEN
        CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ELSE IF (IPOLTP .EQ. 2) THEN
        CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF

      KINDP1  = 1
      KINDP2  = KINDP1  + NDIM
      KLAST   = KINDP2  + NDIM
      LWRK2   = LWRK - KLAST + 1

      DO 100 ITER = 1, MXMMIT
        LM = LM + 1
        DO ISLAVE = 1, NODTOT
          IWHO = -1
          NRUN = 1
          CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
          CALL MPIXSEND(NRUN, 1, 'INTEGER', NWHO, MPTAG2)
        ENDDO

        CALL MPIXBCAST(INDP1,NDIM,'DOUBLE',MASTER)
        CALL DZERO(WRK(KINDP1),NDIM)
        CALL DZERO(WRK(KINDP2),NDIM)

        CALL MPI_REDUCE(WRK(KINDP1),WRK(KINDP2),NDIM, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

        CALL DAXPY(NDIM,1.0D0,WRK(KINDP2),1,INDMOM,1)

        TERROR=0.0D0
        DO I=1,NDIM
          TERROR = TERROR + (INDMOM(I)-INDP1(I))*
     &                      (INDMOM(I)-INDP1(I))
        ENDDO

        IF ( (LOCDEB) .OR. (IPRINT .GE. 15) ) THEN
          LMAX = 0
          TMAX = 0.0D0
          DO I=1,NDIM
            TDIFF = ABS(INDMOM(I)-INDP1(I))
            IF (TDIFF .GT. TMAX) THEN
              TMAX = TDIFF
              LMAX = I 
            ENDIF
          ENDDO
          IF (LMAX .NE. 0) THEN
            WRITE(LUPRI,*) 'Maximum deviation (element) is ',TMAX, LMAX
          ENDIF 
        ENDIF


        IF (ABS(TERROR) .LT. THRESL) THEN
          DIPCON = .TRUE.
          GOTO 200
        ELSE
          DIPCON = .FALSE. 
          IF (LOCDEB )WRITE(LUPRI,*) 'TERROR ',TERROR
          IF ( MMDIIS ) THEN
            CALL DCOPY(NDIM,INDMOM,1,VEC(ITER*NDIM+1),1)
            CALL MM_DIIS_EXTRAPOLATION(VEC,ITER,NDIM,INDP1,
     *                                 WRK(KLAST),LWRK2,IPRINT)
          ELSE 
            CALL DCOPY(NDIM,INDMOM,1,INDP1,1)
          ENDIF
C     If no convergence in last iteration keep the values for the
C     induced dipoles, i.e. not only the diagonal part
          IF (ITER .NE. MXMMIT) CALL DCOPY(NDIM,INDDIA,1,
     *                                       INDMOM,1)
        ENDIF

 100  CONTINUE

 200  CONTINUE                  !Done

C     End message to slaves
      NRUN = -1
      DO ISLAVE = 1, NODTOT
        IWHO = -1
        CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
        CALL MPIXSEND(NRUN, 1, 'INTEGER', NWHO, MPTAG2)
      ENDDO

      CALL QEXIT('MMITER_INDDIP_M')

      RETURN
      END
C******************************************************************************
C  /* Deck mmiter_inddip_s */
      SUBROUTINE MMITER_INDDIP_S(WRK,LWRK,IPRINT)

#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mxcent.h"
#include "qmmm.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include "mpif.h"
#endif

      INTEGER POLDIM, POLARRAY
      DIMENSION WRK(LWRK), POLARRAY(:)
      ALLOCATABLE POLARRAY
      DOUBLE PRECISION DIP,MY 
      DIMENSION DIP(3),MY(3)
      LOGICAL RUN
      
      CALL QENTER('MMITER_INDDIP_S')

      CALL MPIXBCAST(POLDIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NODTOT,1,'INTEGER',MASTER)

      ALLOCATE(POLARRAY(POLDIM))
      CALL MPIXBCAST(POLARRAY,POLDIM,'INTEGER',MASTER)
      CALL MPIXBCAST(IPOLTP,1,'INTEGER',MASTER)
      IF (IPOLTP .EQ. 1) THEN
        CALL MPIXBCAST(POLIMM,MMCENT,'DOUBLE',MASTER)
      ELSE IF (IPOLTP .EQ. 2) THEN
        CALL MPIXBCAST(POLMM,6*MMCENT,'DOUBLE',MASTER)
      ENDIF

      NSLICE = POLDIM/NODTOT
      ISTART = (MYNUM-1)*NSLICE + 1
      IEND   = ISTART + NSLICE - 1
C     check if there is leftovers
      IF ( (NODTOT*NSLICE) .LT. POLDIM) THEN
        LEFT = POLDIM - NODTOT*NSLICE
        IF (MYNUM .LE. LEFT) THEN
          ISTART = ISTART + MYNUM - 1
          IEND   = IEND + MYNUM
        ELSE
          ISTART = ISTART + LEFT
          IEND   = IEND + LEFT
        ENDIF
      END IF

      NDIM = 3*POLDIM

      KINDP1  = 1
      KINDP2  = KINDP1  + NDIM
      KLAST   = KINDP2  + NDIM

      LWRK2 = LWRK - KLAST + 1
      IF (LWRK2 .LT. 0) CALL ERRWRK('MMITER_INDDIP_S',-KLAST,LWRK)

      CALL DZERO(WRK(KINDP2),NDIM)

 20   CONTINUE

      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(IRUN,1,'INTEGER',MASTER,MPTAG2)

      IF (IRUN .EQ. 1) THEN
        CALL DZERO(WRK(KINDP2),NDIM)
        CALL MPIXBCAST(WRK(KINDP1),NDIM,'DOUBLE',MASTER)
        LRI = 1 + 3*(ISTART-1)
        DO L = ISTART, IEND
          I = POLARRAY(L)
          LCI = 1
          DO K = 1, POLDIM
            J = POLARRAY(K)
            CALL GET_MY(I,J,WRK(KINDP1+LCI-1),MY)
            WRK(KINDP2+LRI-1+0) = WRK(KINDP2+LRI-1+0) + MY(1)
            WRK(KINDP2+LRI-1+1) = WRK(KINDP2+LRI-1+1) + MY(2)
            WRK(KINDP2+LRI-1+2) = WRK(KINDP2+LRI-1+2) + MY(3)
            LCI = LCI + 3
          ENDDO
          LRI = LRI + 3
        ENDDO

        CALL MPI_REDUCE(WRK(KINDP2),MPI_IN_PLACE,NDIM, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

        GOTO 20
      ENDIF

      DEALLOCATE(POLARRAY)

      CALL QEXIT('MMITER_INDDIP_S')

      RETURN
      END

#endif
