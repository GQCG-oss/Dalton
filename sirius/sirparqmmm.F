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
C FILE: sirius/sirparqmmm.F
C Magnus and Arnfinn
C**************************************************************************************
C  /* Deck qmmm_fck_ao_p */
      SUBROUTINE QMMM_FCK_AO_P(FSOL,DCAO,ESOLT,WRK,LWRK,IPRINT)   

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>
#include <infpar.h>


      DIMENSION WRK(LWRK)
      DIMENSION FSOL(*),DCAO(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      LOGICAL EXCENT
      INTEGER NZERAL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_FCK_AO_P')

      KTAO   = 1
      KNSEL  = KTAO   + NNBASX
      KNSNUC = KNSEL  + MMCENT
      KWRK1  = KNSNUC + MMCENT
      LWRK1  = LWRK   - KWRK1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMM_FCK_AO_P 1',-KWRK1,LWRK)

      CALL DZERO(WRK(KTAO),NNBASX)
      CALL DZERO(WRK(KNSEL),MMCENT)
      CALL DZERO(WRK(KNSNUC),MMCENT)

C     The different static energy contributions
      ECHART = 0.0D0
      EDIPT  = 0.0D0
      EQUADT = 0.0D0
      EOCTUP = 0.0D0

C     Backup diporg. We use diporg to transfer coordinates to int. program.

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

C     1) The charge correction to the Fock/KS operator

      IF (NMULT .GE. 0) THEN
C         IF (NODTOT .GE. 1) THEN
         IF (.FALSE.) THEN
            CALL QMMM_CHARGE_M()
         ELSE
            CALL QMMM_CHARGE(DCAO,ESOLT,WRK(KTAO),WRK(KNSEL),
     &                       WRK(KNSNUC),WRK(KWRK1),LWRK1,IPRINT)
         ENDIF
      ENDIF

C     2) The dipole correction to the Fock/KS operator

      IF (NMULT .GE. 1) THEN
C         IF (NODTOT .GE. 1) THEN
         IF (.FALSE.) THEN
            CALL QMMM_DIPOLE_M()
         ELSE
            CALL QMMM_DIPOLE(DCAO,ESOLT,WRK(KTAO),WRK(KNSEL),
     &                       WRK(KNSNUC),WRK(KWRK1),LWRK1,IPRINT)
         ENDIF
      ENDIF

C     3) The quadrupole correction to the Fock/KS operator

      IF (NMULT .GE. 2) THEN
C         IF (NODTOT .GE. 1) THEN
         IF (.FALSE.) THEN
            CALL QMMM_QUADPOLE_M(DCAO,ESOLT,WRK(KTAO),WRK(KNSEL),
     &                       WRK(KNSNUC),WRK(KWRK1),LWRK1,IPRINT)
         ELSE
            CALL QMMM_QUADPOLE(DCAO,ESOLT,WRK(KTAO),WRK(KNSEL),
     &                       WRK(KNSNUC),WRK(KWRK1),LWRK1,IPRINT)
         ENDIF
      ENDIF

C     4) The octupole correction to the Fock/KS operator. We do the integrals numerically here.

      IF (NMULT .GE. 3) THEN
C         IF (NODTOT .GE. 1) THEN
         IF (.FALSE.) THEN
            CALL QMMM_OCTUPOLE_M()
         ELSE
            CALL QMMM_OCTUPOLE(DCAO,ESOLT,WRK(KTAO),WRK(KNSEL),
     &                       WRK(KNSNUC),WRK(KWRK1),LWRK1,IPRINT)
         ENDIF
      ENDIF

      IF (IPRINT.GT.1) THEN
        write(lupri,*)
        write(lupri,*) 'MM-charge QM density interaction energy:',ECHART
        write(lupri,*) 'MM-dipole QM density interaction energy:',EDIPT
        write(lupri,*) 'MM-quadr. QM density interaction energy:',EQUADT
        write(lupri,*) 'MM-octup. QM density interaction energy:',EOCTUP
      ENDIF 

C     5) The polarization correction to the Fock/KS operator

      IF (IPOLTP .GT. 0) THEN
C         IF (NODTOT .GE. 1) THEN
         IF (.FALSE.) THEN
            CALL QMMM_POLARI_M()
         ELSE
            CALL QMMM_POLARI(DCAO,ESOLT,EPOLT,WRK(KTAO),
     &                       WRK(KWRK1),LWRK1,IPRINT)
         ENDIF
      ENDIF

C     Finally, put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL PKSYM1(WRK(KTAO),FSOL,NBAS,NSYM,+1)

      CALL QEXIT('QMMM_FCK_AO_P')

      RETURN
      END
C**************************************************************************************
C  /* Deck qmmm_charge */
      SUBROUTINE QMMM_CHARGE(DCAO,ESOLT,TAO,SEL,SNUC,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>


      DIMENSION WRK(LWRK), TAO(NNBASX), SEL(MMCENT), SNUC(MMCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_CHARGE')

      KMAT  = 1
      KLAST = KMAT + NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO_P 2',-KLAST,LWRK)

      FAC1 = 1.0D0      
      EXPNST = 0.0D0
      ECHCH = 0.0D0


      DO 100 I = 1,MMCENT 

         IF (ABS(MUL0MM(I)) .LE. THRMM) THEN
           SEL(I)  = 0.0D0
           SNUC(I) = 0.0D0
           GOTO 100
         ENDIF 

         DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *           (MMCORD(2,I)-QMCOM(2))**2 +
     *           (MMCORD(3,I)-QMCOM(3))**2
         DIST = SQRT(DIST2)

         IF (DIST .GT. RCUTMM) THEN
           SEL(I)  = 0.0D0
           SNUC(I) = 0.0D0
           GOTO 100
         ENDIF

         CALL DZERO(WRK(KMAT),NNBASX)

         KPATOM = 0
         NOSIM = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         DIPORG(1) = MMCORD(1,I)
         DIPORG(2) = MMCORD(2,I)
         DIPORG(3) = MMCORD(3,I)

       
         RUNQM3=.TRUE.
         CALL GET1IN(WRK(KMAT),'NPETES ',NOSIM,WRK(KLAST),
     &               LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
         RUNQM3=.FALSE.

         IF (IPRINT.GT.15) THEN
           WRITE (LUPRI,'(/A)') 'Pot. energy matrix in QMMM_FCK_AO_P'
           CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
         ENDIF

         CALL DSCAL(NNBASX,MUL0MM(I),WRK(KMAT),1)
         EXPNS=DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
         SEL(I) = EXPNS

         CALL DAXPY(NNBASX,FAC1,WRK(KMAT),1,TAO(1),1)
         EXPNST=EXPNST+EXPNS

C        Now the QM nuclear contribution

         ECHCHL  = 0.0D0
         DO 101 J = 1,NUCIND
           XDIS   = CORD(1,J) - MMCORD(1,I) 
           YDIS   = CORD(2,J) - MMCORD(2,I) 
           ZDIS   = CORD(3,J) - MMCORD(3,I) 
           DIST2  = XDIS**2+YDIS**2+ZDIS**2
           DIST   = SQRT(DIST2) 
           ECHCHL = ECHCHL + CHARGE(J)*MUL0MM(I)/DIST
  101    CONTINUE
         ECHCH = ECHCH + ECHCHL
         SNUC(I) = ECHCHL
  100 CONTINUE

      ECHART = EXPNST + ECHCH
      ESOLT  = ECHART

      IF (IPRINT.GT.5) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Charge-electronic Charge-nuclear Total'
        DO 102 I = 1,MMCENT
          ELTEMP = SEL(I) + SNUC(I)
          WRITE(LUPRI,*) I,SEL(I),SNUC(I),ELTEMP
  102   CONTINUE

        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EXPNST, ECHCH, EXPNST+ECHCH
        WRITE(LUPRI,*)
      ENDIF


      CALL QEXIT('QMMM_CHARGE')

      RETURN
      END
C**************************************************************************************
C  /* Deck qmmm_dipole */
      SUBROUTINE QMMM_DIPOLE(DCAO,ESOLT,TAO,SEL,SNUC,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>

      DIMENSION WRK(LWRK), TAO(NNBASX), SEL(MMCENT), SNUC(MMCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_DIPOLE')



      KMAT  = 1
      KLAST = KMAT + 3*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO_P 4',-KLAST,LWRK)

      FAC1   = 1.0D0      
      FACM1  = -1.0D0
      EMUL1T = 0.0D0
      ELOCT  = 0.0D0

      DO 200 I = 1,MMCENT 

C         See if the dipole moment at this site is zero
        DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2     
        DNORM = SQRT(DNORM2)
        IF (ABS(DNORM) .LE. THRMM) THEN
C         we use KNS to collect the energy terms
          SEL(I)  = 0.0D0
          SNUC(I) = 0.0D0
          GOTO 200
        ENDIF 

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
          SEL(I)  = 0.0D0
          SNUC(I) = 0.0D0
          GOTO 200
        ENDIF

        CALL DZERO(WRK(KMAT),3*NNBASX)

        KPATOM = 0
        NOSIM = 3
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        DIPORG(1) = MMCORD(1,I)
        DIPORG(2) = MMCORD(2,I)
        DIPORG(3) = MMCORD(3,I)


        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
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
     &                (DIPORG(2)-CORD(2,M))**2 +
     &                (DIPORG(3)-CORD(3,M))**2
              IF (DISTC .LE. DIST) THEN
                DIST = DISTC
                MHIT = M
              ENDIF
  123       CONTINUE
          ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &             (DIPORG(2)-QMCOM(2))**2 +
     &             (DIPORG(3)-QMCOM(3))**2
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

c          write(lupri,*) 'TEMPI',TEMPI
c          write(lupri,*) 'QMPOL(MHIT)',QMPOL(MHIT)
c          write(lupri,*) 'MHIT',MHIT
c          write(lupri,*) 'TEMP',TEMP
c          write(lupri,*) 'SIJ',SIJ
c          write(lupri,*) 'DFACT',DFACT

        IF (IPRINT.GT.15) THEN
          WRITE (LUPRI,'(/A)') ' E_x_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_y matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_z matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
        END IF

        CALL DSCAL(NNBASX,MUL1MM(1,I),WRK(KMAT),1)
        CALL DSCAL(NNBASX,MUL1MM(2,I),WRK(KMAT+NNBASX),1)
        CALL DSCAL(NNBASX,MUL1MM(3,I),WRK(KMAT+2*NNBASX),1)

        CALL DAXPY(NNBASX,FACM1,WRK(KMAT),1,TAO(1),1)
        CALL DAXPY(NNBASX,FACM1,WRK(KMAT+NNBASX),1,TAO(1),1)
        CALL DAXPY(NNBASX,FACM1,WRK(KMAT+2*NNBASX),1,TAO(1),1)

        EXCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
        EYCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
        EZCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

        SEL(I) = EXCOMP + EYCOMP + EZCOMP
        EMUL1T = EMUL1T + SEL(I)

C       Now the QM nuclear contribution. Note that we switch the sign here by writing 
C       CORD(1,J) - MMCORD(1,I)

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
  201   CONTINUE
        SNUC(I) = ELOC
        ELOCT = ELOCT + SNUC(I)

  200 CONTINUE

      EDIPT = EMUL1T + ELOCT
      ESOLT = ESOLT + EDIPT

      IF (IPRINT.GT.5) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Dipole-electronic Dipole-nuclear Total'
        DO 202 I = 1,MMCENT
          ETEMP = SEL(I) + SNUC(I)
          WRITE(LUPRI,*) I,SEL(I),SNUC(I),ETEMP
  202   CONTINUE
  
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EMUL1T, ELOCT, EMUL1T+ELOCT
        WRITE(LUPRI,*) 
      ENDIF

      CALL QEXIT('QMMM_DIPOLE')

      RETURN
      END
C**************************************************************************************
C  /* Deck qmmm_quadpole */
      SUBROUTINE QMMM_QUADPOLE(DCAO,ESOLT,TAO,SEL,SNUC,WRK,LWRK,
     &                         IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>

      DIMENSION WRK(LWRK), TAO(NNBASX), SEL(MMCENT), SNUC(MMCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_QUADPOLE')


      KMAT  = 1
      KLAST = KMAT + 6*NNBASX
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO_P 6',-KLAST,LWRK)

      FAC1   = 1.0D0
      FACM1  = -1.0D0
      EMUL2T = 0.0D0
      ELOCT  = 0.0D0

      DO 300 I = 1,MMCENT

C       See if the quadrupole moment at this site is zero
        DNORM2 = MUL2MM(1,I)**2+MUL2MM(2,I)**2+MUL2MM(3,I)**2
     *         + MUL2MM(4,I)**2+MUL2MM(5,I)**2+MUL2MM(6,I)**2
        DNORM = SQRT(DNORM2)
        IF (ABS(DNORM) .LE. THRMM) THEN
          SEL(I)  = 0.0D0
          SNUC(I) = 0.0D0
C         we use KNS to collect the energy terms
          GOTO 300
        ENDIF

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
          SEL(I)  = 0.0D0
          SNUC(I) = 0.0D0
          GOTO 300
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
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.

        IF (IPRINT.GT.15) THEN
          WRITE (LUPRI,'(/A)') ' E_xx_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xy matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xz matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yy_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+3*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yz_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+4*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_zz_matrix in QMMM_FCK_AO_P:'
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
C       The integrals contains a factor of -1. Therefore daxpy with fac1
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT),1,TAO(1),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+NNBASX),1,TAO(1),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+2*NNBASX),1,TAO(1),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+3*NNBASX),1,TAO(1),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+4*NNBASX),1,TAO(1),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+5*NNBASX),1,TAO(1),1)

C       Contract with the density to get the expectation values. 
C       The factor of 1/2 in the Taylor expansion has been included.
C       Also, the off-diagonal elements have been scaled by 2 in order 
C       to include all contributions (the off -diagonal parts are related by symmetry)

C       Since the integrals contains a factor of -1 no -DDOT here.

        EMU2XX=DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
        EMU2XY=DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
        EMU2XZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)
        EMU2YY=DDOT(NNBASX,DCAO,1,WRK(KMAT+3*NNBASX),1)
        EMU2YZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+4*NNBASX),1)
        EMU2ZZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+5*NNBASX),1)

        SEL(I) = EMU2XX+EMU2XY+EMU2XZ
     *                 + EMU2YY+EMU2YZ+EMU2ZZ
        EMUL2T = EMUL2T + SEL(I)

C       Now the QM nuclear contribution

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
  301   CONTINUE

C       Remember the factor of 1/2 from the Taylor expansion
        ELOC   = 0.5D0*ELOC

        SNUC(I) = ELOC
        ELOCT = ELOCT + SNUC(I)

  300 CONTINUE

      EQUADT = EMUL2T + ELOCT
      ESOLT  = ESOLT  + EQUADT

      IF (IPRINT.GT.5) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Quadr.-electronic Quadr.-nuclear Total'
        DO 302 I = 1,MMCENT
          ETEMP = SEL(I) + SNUC(I)
          WRITE(LUPRI,*) I,SEL(I),SNUC(I),ETEMP
  302   CONTINUE
  
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EMUL2T, ELOCT, EMUL2T+ELOCT
        WRITE(LUPRI,*) 
      ENDIF

      CALL QEXIT('QMMM_QUADPOLE')

      RETURN
      END

C**************************************************************************************
C  /* Deck qmmm_octupole */
      SUBROUTINE QMMM_OCTUPOLE(DCAO,ESOLT,TAO,SEL,SNUC,WRK,LWRK,
     &                         IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>

      DIMENSION WRK(LWRK), TAO(NNBASX), SEL(MMCENT), SNUC(MMCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_OCTUPOLE')

      KMAT1  = 1
      KMAT2  = KMAT1 + 6*NNBASX
      KMAT3  = KMAT2 + 6*NNBASX
      KMAT4  = KMAT3 + 6*NNBASX
      KLAST  = KMAT4 + 6*NNBASX
      LWRK2  = LWRK  - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO_P 7',-KLAST,LWRK)

      FAC1   = 1.0D0
      FACM1  = -1.0D0
      EMUL3T = 0.0D0
      ELOCT  = 0.0D0

      DO 410 I = 1,MMCENT

C       See if the octupole moment at this site is zero
        DNORM2 = MUL3MM(1,I)**2+MUL3MM(2,I)**2+MUL3MM(3,I)**2
     *         + MUL3MM(4,I)**2+MUL3MM(5,I)**2+MUL3MM(6,I)**2
     *         + MUL3MM(7,I)**2+MUL3MM(8,I)**2+MUL3MM(9,I)**2
     *         + MUL3MM(10,I)**2
        DNORM = SQRT(DNORM2)
        IF (ABS(DNORM) .LE. THRMM) THEN
          SEL(I)  = 0.0D0
          SNUC(I) = 0.0D0
C         we use KNS to collect the energy terms
          GOTO 410
        ENDIF

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
          SEL(I)  = 0.0D0
          SNUC(I) = 0.0D0
          GOTO 410
        ENDIF

        CALL DZERO(WRK(KMAT1),6*NNBASX)
        CALL DZERO(WRK(KMAT2),6*NNBASX)
        CALL DZERO(WRK(KMAT3),6*NNBASX)
        CALL DZERO(WRK(KMAT4),6*NNBASX)

        KPATOM = 0
        NOSIM = 6
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        DIPORG(1) = MMCORD(1,I)
        DIPORG(2) = MMCORD(2,I)
        DIPORG(3) = MMCORD(3,I)

        XSAVE = DIPORG(1)
        YSAVE = DIPORG(2)
        ZSAVE = DIPORG(3)

C       Consider first the derivative wrt. z 

        DIPORG(3) = DIPORG(3) + DELFLD
        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT1),'ELFGRDC',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.
        DIPORG(3) = ZSAVE

        DIPORG(3) = DIPORG(3) - DELFLD
        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT2),'ELFGRDC',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.
        DIPORG(3) = ZSAVE

C       Calculate gradient of these integrals

        FSCAL = 1.0D0/(2.0D0*DELFLD)
        CALL DAXPY(6*NNBASX,FACM1,WRK(KMAT2),1,WRK(KMAT1),1)
        CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT1),1)
        CALL DZERO(WRK(KMAT2),6*NNBASX)

        IF (IPRINT.GT.15) THEN
          WRITE (LUPRI,'(/A)') ' E_xxz_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT1),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xyz matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT1+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xzz matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT1+2*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yyz_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT1+3*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yzz_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT1+4*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_zzz_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT1+5*NNBASX),NBAST,1,LUPRI)
        END IF

C       Continue with the derivatives wrt. y   

        DIPORG(2) = DIPORG(2) + DELFLD
        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT2),'ELFGRDC',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.
        DIPORG(2) = YSAVE

        DIPORG(2) = DIPORG(2) - DELFLD
        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT3),'ELFGRDC',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.
        DIPORG(2) = YSAVE

C       Calculate gradient of these integrals

        FSCAL = 1.0D0/(2.0D0*DELFLD)
        CALL DAXPY(6*NNBASX,FACM1,WRK(KMAT3),1,WRK(KMAT2),1)
        CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT2),1)
        CALL DZERO(WRK(KMAT3),6*NNBASX)

        IF (IPRINT.GT.15) THEN
          WRITE (LUPRI,'(/A)') ' E_xxy_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT2),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xyy matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT2+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xzy matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT2+2*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yyy_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT2+3*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yzy_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT2+4*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_zzy_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT2+5*NNBASX),NBAST,1,LUPRI)
        END IF

C       Continue with the derivatives wrt. x

        DIPORG(1) = DIPORG(1) + DELFLD
        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT3),'ELFGRDC',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.
        DIPORG(1) = XSAVE

        DIPORG(1) = DIPORG(1) - DELFLD
        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT4),'ELFGRDC',NOSIM,WRK(KLAST),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.
        DIPORG(1) = XSAVE

C       Calculate gradient of these integrals

        FSCAL = 1.0D0/(2.0D0*DELFLD)
        CALL DAXPY(6*NNBASX,FACM1,WRK(KMAT4),1,WRK(KMAT3),1)
        CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT3),1)
        CALL DZERO(WRK(KMAT4),6*NNBASX)

        IF (IPRINT.GT.15) THEN
          WRITE (LUPRI,'(/A)') ' E_xxx_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT3),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xyx matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT3+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xzx matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT3+2*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yyx_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT3+3*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yzx_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT3+4*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_zzx_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT3+5*NNBASX),NBAST,1,LUPRI)
        END IF

        CALL DSCAL(NNBASX,MUL3MM(1,I),WRK(KMAT3),1) ! xxx
        CALL DSCAL(NNBASX,D3*MUL3MM(2,I),WRK(KMAT2),1) ! xxy
        CALL DSCAL(NNBASX,D3*MUL3MM(3,I),WRK(KMAT1),1) ! xxz
        CALL DSCAL(NNBASX,D3*MUL3MM(4,I),WRK(KMAT2+NNBASX),1) ! xyy
        CALL DSCAL(NNBASX,D6*MUL3MM(5,I),WRK(KMAT1+NNBASX),1) ! xyz
        CALL DSCAL(NNBASX,D3*MUL3MM(6,I),WRK(KMAT1+2*NNBASX),1) ! xzz
        CALL DSCAL(NNBASX,MUL3MM(7,I),WRK(KMAT2+3*NNBASX),1) ! yyy
        CALL DSCAL(NNBASX,D3*MUL3MM(8,I),WRK(KMAT1+3*NNBASX),1) ! yyz
        CALL DSCAL(NNBASX,D3*MUL3MM(9,I),WRK(KMAT1+4*NNBASX),1) ! yzz
        CALL DSCAL(NNBASX,MUL3MM(10,I),WRK(KMAT1+5*NNBASX),1) ! zzz

C       Include 1/6 from the Taylor expansion
        CALL DSCAL(6*NNBASX,D6I,WRK(KMAT1),1)
        CALL DSCAL(6*NNBASX,D6I,WRK(KMAT2),1)
        CALL DSCAL(6*NNBASX,D6I,WRK(KMAT3),1)
C
C       The integrals contains a factor of -1 since they are obtained as derivative if 
C       the field gradient and these contains a factor of -1. Therefore no -1 in the daxpy.
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT3),1,TAO(1),1) ! xxx 
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT2),1,TAO(1),1) ! xxy
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT1),1,TAO(1),1) ! xxz
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT2+NNBASX),1,TAO(1),1) ! xyy
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+NNBASX),1,TAO(1),1) ! xyz
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+2*NNBASX),1,TAO(1),1) ! xzz 
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT2+3*NNBASX),1,TAO(1),1) ! yyy
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+3*NNBASX),1,TAO(1),1) ! yyz
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+4*NNBASX),1,TAO(1),1) ! yzz
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+5*NNBASX),1,TAO(1),1) ! zzz

C       Contract with the density to get the expectation values. 
C       The factor of 1/6 in the Taylor expansion has been included.
C       Also, the off-diagonal elements have been scaled by 3 and the xyz element by 6 in order 
C       to include all contributions (the off -diagonal parts are related by symmetry)

C       Since the integrals contains a factor of -1 no -DDOT here.

        E3XXX=DDOT(NNBASX,DCAO,1,WRK(KMAT3),1)
        E3XXY=DDOT(NNBASX,DCAO,1,WRK(KMAT2),1)
        E3XXZ=DDOT(NNBASX,DCAO,1,WRK(KMAT1),1)
        E3XYY=DDOT(NNBASX,DCAO,1,WRK(KMAT2+NNBASX),1)
        E3XYZ=DDOT(NNBASX,DCAO,1,WRK(KMAT1+NNBASX),1)
        E3XZZ=DDOT(NNBASX,DCAO,1,WRK(KMAT1+2*NNBASX),1)
        E3YYY=DDOT(NNBASX,DCAO,1,WRK(KMAT2+3*NNBASX),1)
        E3YYZ=DDOT(NNBASX,DCAO,1,WRK(KMAT1+3*NNBASX),1)
        E3YZZ=DDOT(NNBASX,DCAO,1,WRK(KMAT1+4*NNBASX),1)
        E3ZZZ=DDOT(NNBASX,DCAO,1,WRK(KMAT1+5*NNBASX),1)

        SEL(I) = E3XXX+E3XXY+E3XXZ+E3XYY+E3XYZ
     *         + E3XZZ+E3YYY+E3YYZ+E3YZZ+E3ZZZ
        EMUL3T = EMUL3T + SEL(I)

C       Now the QM nuclear contribution. Note that we switch the sign here by writing
C       CORD(1,J) - MMCORD(1,I)


        ELOC     = 0.0D0
        DO 411 J = 1,NUCIND
          XDIS   = CORD(1,J) - MMCORD(1,I)  
          YDIS   = CORD(2,J) - MMCORD(2,I) 
          ZDIS   = CORD(3,J) - MMCORD(3,I)  
          DIST2  = XDIS**2+YDIS**2+ZDIS**2
          DIST   = SQRT(DIST2)
          DIST3  = DIST2*DIST
          DIST5  = DIST3*DIST2
          DIST7  = DIST5*DIST2
C
          TXXX=(15.0D0*XDIS*XDIS*XDIS-3*DIST2*(XDIS+XDIS+XDIS))/DIST7
          TXXY=(15.0D0*XDIS*XDIS*YDIS-3*DIST2*(          YDIS))/DIST7
          TXXZ=(15.0D0*XDIS*XDIS*ZDIS-3*DIST2*(          ZDIS))/DIST7
          TXYY=(15.0D0*XDIS*YDIS*YDIS-3*DIST2*(XDIS          ))/DIST7
          TXYZ=(15.0D0*XDIS*YDIS*ZDIS)/DIST7
          TXZZ=(15.0D0*XDIS*ZDIS*ZDIS-3*DIST2*(XDIS          ))/DIST7
          TYYY=(15.0D0*YDIS*YDIS*YDIS-3*DIST2*(YDIS+YDIS+YDIS))/DIST7
          TYYZ=(15.0D0*YDIS*YDIS*ZDIS-3*DIST2*(          ZDIS))/DIST7
          TYZZ=(15.0D0*YDIS*ZDIS*ZDIS-3*DIST2*(YDIS          ))/DIST7
          TZZZ=(15.0D0*ZDIS*ZDIS*ZDIS-3*DIST2*(ZDIS+ZDIS+ZDIS))/DIST7

          ELOC   =   ELOC
     *           +   CHARGE(J)*MUL3MM(1,I)*TXXX
     *           +   CHARGE(J)*MUL3MM(2,I)*TXXY*3
     *           +   CHARGE(J)*MUL3MM(3,I)*TXXZ*3
     *           +   CHARGE(J)*MUL3MM(4,I)*TXYY*3
     *           +   CHARGE(J)*MUL3MM(5,I)*TXYZ*6
     *           +   CHARGE(J)*MUL3MM(6,I)*TXZZ*3
     *           +   CHARGE(J)*MUL3MM(7,I)*TYYY
     *           +   CHARGE(J)*MUL3MM(8,I)*TYYZ*3
     *           +   CHARGE(J)*MUL3MM(9,I)*TYZZ*3
     *           +   CHARGE(J)*MUL3MM(10,I)*TZZZ
  411   CONTINUE

C       Remember the factor of 1/6 from the Taylor expansion
        ELOC   = D6I*ELOC

        SNUC(I) = ELOC
        ELOCT = ELOCT + SNUC(I)

  410 CONTINUE

      EOCTUP = EMUL3T + ELOCT
      ESOLT  = ESOLT  + EOCTUP

      IF (IPRINT.GT.5) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Octup.-electronic Octup.-nuclear Total'
        DO 412 I = 1,MMCENT
          ETEMP = SEL(I) + SNUC(I)
          WRITE(LUPRI,*) I,SEL(I),SNUC(I),ETEMP
  412   CONTINUE
  
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EMUL3T, ELOCT, EMUL3T+ELOCT
        WRITE(LUPRI,*) 
      ENDIF


      CALL QEXIT('QMMM_OCTUPOLE')

      RETURN
      END

C**************************************************************************************
C  /* Deck qmmm_polari */
      SUBROUTINE QMMM_POLARI(DCAO,ESOLT,EPOLT,TAO,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>


      DIMENSION WRK(LWRK), TAO(NNBASX)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, EXCENT
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

      IF (IPRINT.GT.1) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Number of polarizable sites: ', NNZAL
        WRITE(LUPRI,*)
      ENDIF

      IF (MMMAT) THEN

C       There is no need to construct ALPHA^(-1) and T matrices
C       explicitly, we can directly form the [ALPHA^(-1) - T] matrix
        KALINMT = 1
        KINVAMT = KALINMT + (3*NNZAL)*(3*NNZAL) ! [ALPHA^(-1) - T] matrix
        KUMAT   = KINVAMT + (3*NNZAL)*(3*NNZAL) ! [ALPHA^(-1) - T]^(-1) matrix
        KINDMOM = KUMAT   + (3*NNZAL)*(3*NNZAL) ! UMAT is just for testing
        KMAT    = KINDMOM + (3*NNZAL) ! List for induced dipoles
        KIPVT   = KMAT    + (3*NNBASX) ! For Rr_a integrals
        KWRKV   = KIPVT   + (3*NNZAL) ! For matrix inv.
        KWRK2   = KWRKV   + (3*NNZAL) ! For matrix inv.
        LWRK2   = LWRK    - KWRK2 + 1

        IF (LWRK2 .LT. 0) THEN
          CALL ERRWRK('QMMM_FCK_AO_P pol. part 1',-KWRK2,LWRK)
        ENDIF

        CALL DZERO(WRK(KALINMT),(3*NNZAL)*(3*NNZAL))
        CALL DZERO(WRK(KINVAMT),(3*NNZAL)*(3*NNZAL))
        CALL DZERO(WRK(KUMAT),(3*NNZAL)*(3*NNZAL))
        CALL DZERO(WRK(KINDMOM),(3*NNZAL))
        CALL DZERO(WRK(KIPVT),(3*NNZAL))
        CALL DZERO(WRK(KWRKV),(3*NNZAL))
        CALL DZERO(WRK(KMAT),(3*NNBASX))

C       FIXDIP assumes induced dipoles are calculated in a previous run.
C       Mainly due to debugging. Assumes identical molecules and order
C       of atoms in previous and current run.

        IF (.NOT. FIXDIP) THEN
          CALL GET_IND_DIPOLES_2(DCAO,NNZAL,
     &                           WRK(KALINMT),WRK(KINVAMT),WRK(KUMAT),
     &                           WRK(KINDMOM),
     &                           WRK(KWRK2),WRK(KIPVT),WRK(KWRKV),
     &                           LWRK2,IPRINT)
        ELSE
          WRITE(LUPRI,*) 'Ind. dips. from a prev. calc. read from file'
          CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WRK(KINDMOM))
        ENDIF

      ELSE IF (MMITER) THEN

        KINDMOM = 1                               ! List for induced dipoles
        KMAT    = KINDMOM + (3*NNZAL)
        KWRK2   = KMAT    + (3*NNBASX)             
        LWRK2   = LWRK    - KWRK2 + 1

        IF (LWRK2 .LT. 0) THEN
          CALL ERRWRK('QMMM_FCK_AO_P pol. part 2',-KWRK2,LWRK)
        ENDIF

        CALL DZERO(WRK(KINDMOM),(3*NNZAL))

C       FIXDIP assumes induced dipoles are calculated in a previous run.
C       Mainly due to debugging. Assumes identical molecules and order
C       of atoms in previous and current run.

        IF (.NOT. FIXDIP) THEN
          CALL GET_IND_DIPOLES_1(DCAO,NNZAL,WRK(KINDMOM),WRK(KWRK2),
     &                             LWRK2,IPRINT)
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
      ED3MOM = 0.0D0            ! For interaction with octupoles
      EDMULT = 0.0D0            ! For interaction with permanent multipoles
      EPOLT  = 0.0D0            ! Total polarization energy

      DO 500 I=1,MMCENT

        IF (ZEROAL(I) .EQ. -1) GOTO 500

        DNORM2 = WRK(KINDMOM+IINIM+0)**2 + 
     *           WRK(KINDMOM+IINIM+1)**2 +
     *           WRK(KINDMOM+IINIM+2)**2 
        DNORM = SQRT(DNORM2)
        IF (ABS(DNORM) .LE. THRMM) THEN
          IINIM = IINIM+3
          GOTO 500
        ENDIF

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
          IINIM = IINIM+3
          GOTO 500
        ENDIF

        CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+0),WRK(KMAT),1)
        CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+1),WRK(KMAT+NNBASX),1)
        CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+2),WRK(KMAT+2*NNBASX),1)
        CALL DZERO(WRK(KMAT),3*NNBASX)

        KPATOM = 0
        NOSIM = 3
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        DIPORG(1) = MMCORD(1,I)
        DIPORG(2) = MMCORD(2,I)
        DIPORG(3) = MMCORD(3,I)

        RUNQM3 = .TRUE.
        CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIM,WRK(KWRK2),
     &              LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3 = .FALSE.

        IF (QMDAMP) THEN
          IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
          ENDIF
          IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) )THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 124 M=1,NUCIND
              DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                (DIPORG(2)-CORD(2,M))**2 +
     &                (DIPORG(3)-CORD(3,M))**2
              IF (DISTC .LE. DIST) THEN
                DIST = DISTC
                MHIT = M
              ENDIF
  124       CONTINUE
          ELSE IF (IDAMP .EQ. 2) THEN
            DIST = (DIPORG(1)-QMCOM(1))**2 +
     &             (DIPORG(2)-QMCOM(2))**2 +
     &             (DIPORG(3)-QMCOM(3))**2
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

        IF (IPRINT.GT.15) THEN
          WRITE (LUPRI,'(/A)') ' E_x_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_y matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_z matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
        END IF

        CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+0),WRK(KMAT),1)
        CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+1),WRK(KMAT+NNBASX),1)
        CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+2),WRK(KMAT+2*NNBASX),1)

        CALL DAXPY(NNBASX,FACM1,WRK(KMAT),1,TAO(1),1)
        CALL DAXPY(NNBASX,FACM1,WRK(KMAT+NNBASX),1,TAO(1),1)
        CALL DAXPY(NNBASX,FACM1,WRK(KMAT+2*NNBASX),1,TAO(1),1)

C       Polarization contribution to the total energy 

C       A) Electronic contribution

        EXCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
        EYCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
        EZCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

        ET = 0.0D0
        ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
        EDELD = EDELD + ET

C       B) Nuclear contribution

        EFNUCX = 0.0D0
        EFNUCY = 0.0D0
        EFNUCZ = 0.0D0

        DO 510 J=1,NUCIND
          CALL GET_CHARGE_ELFLD(CHARGE(J),
     &                          CORD(1,J),CORD(2,J),CORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
          EFNUCX = EFNUCX + ELFLDX
          EFNUCY = EFNUCY + ELFLDY
          EFNUCZ = EFNUCZ + ELFLDZ
 510    CONTINUE

        IF (QMDAMP) THEN
          IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
          ENDIF
          IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIQM = 9.99D+99
            MHIT = 0
            DO 125 M=1,NUCIND
              DIQMC = (MMCORD(1,I)-CORD(1,M))**2 +
     &                (MMCORD(2,I)-CORD(2,M))**2 +
     &                (MMCORD(3,I)-CORD(3,M))**2
              IF (DIQMC .LE. DIQM) THEN
                DIQM = DIQMC
                MHIT = M
              ENDIF
  125       CONTINUE
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

        EXCOMP = WRK(KINDMOM+IINIM+0)*EFNUCX
        EYCOMP = WRK(KINDMOM+IINIM+1)*EFNUCY
        EZCOMP = WRK(KINDMOM+IINIM+2)*EFNUCZ

        ET = 0.0D0
        ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
        EDNUC = EDNUC + ET

C       C) Multipole contribution

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

C       Get electric fields due to permanent moments

        DO 520 J=1,MMCENT

          IF (J .EQ. I) GOTO 520

          EXCENT = .FALSE.

          DO 522 L=1,NEXLST
            IF (EXLIST(L,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
 522      CONTINUE

          IF (.NOT. EXCENT) THEN

C           C.1  Point-charge contribution

            IF ( (NMULT .GE. 0) .AND. 
     &         (ABS(MUL0MM(J)) .GT. THRMM) ) THEN

              CALL GET_CHARGE_ELFLD(MUL0MM(J),
     &                        MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                        MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                        ELFLDX,ELFLDY,ELFLDZ)

              EF0MX = EF0MX + ELFLDX
              EF0MY = EF0MY + ELFLDY
              EF0MZ = EF0MZ + ELFLDZ

            ENDIF

C           C.2  Dipole contribution

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

C           C.3  Quadrupole contribution

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

C           C.3  Octupole contribution

            IF (NMULT .GE. 3) THEN

              CALL GET_OCTUPOLE_ELFLD(
     &                        MUL3MM(1,J),MUL3MM(2,J),MUL3MM(3,J),
     &                        MUL3MM(4,J),MUL3MM(5,J),MUL3MM(6,J),
     &                        MUL3MM(7,J),MUL3MM(8,J),MUL3MM(9,J),
     &                        MUL3MM(10,J),
     &                        MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                        MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                        ELFLDX,ELFLDY,ELFLDZ)

              EF3MX = EF3MX + ELFLDX
              EF3MY = EF3MY + ELFLDY
              EF3MZ = EF3MZ + ELFLDZ

            ENDIF

          ENDIF

 520    CONTINUE

C       Point-charge contribution

        IF (NMULT .GE. 0) THEN
            
          EXCOMP = WRK(KINDMOM+IINIM+0)*EF0MX
          EYCOMP = WRK(KINDMOM+IINIM+1)*EF0MY
          EZCOMP = WRK(KINDMOM+IINIM+2)*EF0MZ

          ET = 0.0D0
          ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
          ED0MOM = ED0MOM + ET

        ENDIF

C         Dipole contribution

        IF (NMULT .GE. 1) THEN

          EXCOMP = WRK(KINDMOM+IINIM+0)*EF1MX
          EYCOMP = WRK(KINDMOM+IINIM+1)*EF1MY
          EZCOMP = WRK(KINDMOM+IINIM+2)*EF1MZ

          ET = 0.0D0
          ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
          ED1MOM = ED1MOM + ET

        ENDIF

C       Quadrupole contribution

        IF (NMULT .GE. 2) THEN

          EXCOMP = WRK(KINDMOM+IINIM+0)*EF2MX
          EYCOMP = WRK(KINDMOM+IINIM+1)*EF2MY
          EZCOMP = WRK(KINDMOM+IINIM+2)*EF2MZ

          ET = 0.0D0
          ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
          ED2MOM = ED2MOM + ET

        ENDIF

C       Octupole contribution

        IF (NMULT .GE. 3) THEN

          EXCOMP = WRK(KINDMOM+IINIM+0)*EF3MX
          EYCOMP = WRK(KINDMOM+IINIM+1)*EF3MY
          EZCOMP = WRK(KINDMOM+IINIM+2)*EF3MZ

          ET = 0.0D0
          ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
          ED3MOM = ED3MOM + ET

        ENDIF


        IINIM = IINIM + 3

 500  CONTINUE

      EDMULT = ED0MOM + ED1MOM + ED2MOM + ED3MOM
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
C**************************************************************************************
C  /* Deck qmmm_charge_m */
      SUBROUTINE QMMM_CHARGE_M()

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>


      CALL QENTER('QMMM_CHARGE_M')


      CALL QEXIT('QMMM_CHARGE_M')

      RETURN
      END

C**************************************************************************************
C  /* Deck qmmm_dipole_m */
      SUBROUTINE QMMM_DIPOLE_M()

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>


      CALL QENTER('QMMM_DIPOLE_M')


      CALL QEXIT('QMMM_DIPOLE_M')

      RETURN
      END
C**************************************************************************************
C  /* Deck qmmm_quadpole_m */
      SUBROUTINE QMMM_QUADPOLE_M(DCAO,ESOLT,TAO,SEL,SNUC,WRK,LWRK,
     &                         IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>


      DIMENSION WRK(LWRK), TAO(NNBASX), SEL(MMCENT), SNUC(MMCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT), DCAO(*)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_QUADPOLE_M')


      KNSNUC = 1
      KTAO  = KNSNUC + MMCENT
      KNSEL = KTAO + NNBASX
      KLAST = KNSEL + MMCENT
      LWRK2 = LWRK - KLAST + 1

      IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO_P 6',-KLAST,LWRK)

      FAC1   = 1.0D0
      FACM1  = -1.0D0

C     Wake up slaves (Rock and roll all nite)

      IPRTYP = 25
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPQMMM,1,'INTEGER',MASTER)

C     Send data to slaves (Lick it up)

C     These will be sent before. Not necessary to resend?
      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
C **************************************************


      CALL MPIXBCAST(TAO(1),NNBASX,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL2MM,6*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)

C     From old routine
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(DCAO,NNBASX,'DOUBLE',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)



      CALL DZERO(TAO(1),NNBASX)

C     Start parallelized loop
      LNUM = 0
      DO L = 1,MMCENT
        LNUM = LNUM + 1
        IWHO = -1
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(LNUM,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Send end message to all slaves

      LEND = -1
      DO ISLAVE = 1,NODTOT
        IWHO = -1
        CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
        CALL MPIXSEND(LEND,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves (I love it load)

      CALL DZERO(WRK(KNSNUC),MMCENT)
      CALL MPI_REDUCE(WRK(KNSNUC),SNUC,MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL DZERO(WRK(KNSEL),MMCENT)
      CALL MPI_REDUCE(WRK(KNSEL),SEL,MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL DZERO(WRK(KTAO),NNBASX)
      CALL MPI_REDUCE(WRK(KTAO),TAO,NNBASX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      EMUL2T = 0.0D0
      ELOCT  = 0.0D0
      DO I = 1, MMCENT
         ELOCT = ELOCT + SNUC(I)
         EMUL2T = EMUL2T + SEL(I)
      END DO

      EQUADT = EMUL2T + ELOCT
      ESOLT  = ESOLT  + EQUADT

      IF (IPRINT.GT.5) THEN
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Center Quadr.-electronic Quadr.-nuclear Total'
        DO 302 I = 1,MMCENT
          ETEMP = SEL(I) + SNUC(I)
          WRITE(LUPRI,*) I,SEL(I),SNUC(I),ETEMP
  302   CONTINUE
  
        WRITE(LUPRI,*)
        WRITE(LUPRI,*) ' Total '
        WRITE(LUPRI,*) EMUL2T, ELOCT, EMUL2T+ELOCT
        WRITE(LUPRI,*) 
      ENDIF

      CALL QEXIT('QMMM_QUADPOLE_M')

      RETURN
      END
C**************************************************************************************
C  /* Deck qmmm_quadpole_s */
      SUBROUTINE QMMM_QUADPOLE_S(WRK,LWRK,IPRTMP)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>
C ----
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>
#include <gnrinf.h>


      DIMENSION WRK(LWRK)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_QUAD_S')

      FAC1   = 1.0D0
      FACM1  = -1.0D0


      QMMM = .TRUE.
C     Receiving data from master (I was made for lovin' you)

      CALL MPIXBCAST(MMCENT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)

      KNSNUC = 1
      KNSEL  = KNSNUC + MMCENT
      KTAO   = KNSEL  + MMCENT
      KDCAO  = KTAO   + NNBASX
      KMAT   = KDCAO  + NNBASX
      KLAST  = KMAT   + 6*NNBASX

      LWRK1  = LWRK - KLAST + 1

      CALL MPIXBCAST(WRK(KTAO),NNBASX,'DOUBLE',MASTER)
      CALL MPIXBCAST(MUL2MM,6*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(MMCORD,3*MMCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(QMCOM,3,'DOUBLE',MASTER)

      CALL MPIXBCAST(RCUTMM,1,'DOUBLE',MASTER)

C     From old routine
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KDCAO),NNBASX,'DOUBLE',MASTER)

      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

C     Do the work

      CALL DZERO(WRK(KNSEL),MMCENT)
      CALL DZERO(WRK(KNSNUC),MMCENT)

 20   CONTINUE

      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(I,1,'INTEGER',MASTER,MPTAG2)

      IF (I.GT.0) THEN
C       See if the quadrupole moment at this site is zero
        DNORM2 = MUL2MM(1,I)**2+MUL2MM(2,I)**2+MUL2MM(3,I)**2
     *         + MUL2MM(4,I)**2+MUL2MM(5,I)**2+MUL2MM(6,I)**2
        DNORM = SQRT(DNORM2)
        IF (ABS(DNORM) .LE. THRMM) THEN
          WRK(KNSEL + I - 1)  = 0.0D0
          WRK(KNSNUC + I - 1) = 0.0D0
C         we use KNS to collect the energy terms
          GOTO 20
        ENDIF

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
          WRK(KNSEL + I - 1)  = 0.0D0
          WRK(KNSNUC + I - 1) = 0.0D0
          GOTO 20
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
     &              LWRK1,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRTMP)
        RUNQM3=.FALSE.

        IF (IPRTMP.GT.15) THEN
          WRITE (LUPRI,'(/A)') ' E_xx_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xy matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_xz matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yy_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+3*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_yz_matrix in QMMM_FCK_AO_P:'
          CALL OUTPAK(WRK(KMAT+4*NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_zz_matrix in QMMM_FCK_AO_P:'
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
C       The integrals contains a factor of -1. Therefore daxpy with fac1
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT),1,WRK(KTAO),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+NNBASX),1,WRK(KTAO),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+2*NNBASX),1,WRK(KTAO),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+3*NNBASX),1,WRK(KTAO),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+4*NNBASX),1,WRK(KTAO),1)
        CALL DAXPY(NNBASX,FAC1,WRK(KMAT+5*NNBASX),1,WRK(KTAO),1)

C       Contract with the density to get the expectation values. 
C       The factor of 1/2 in the Taylor expansion has been included.
C       Also, the off-diagonal elements have been scaled by 2 in order 
C       to include all contributions (the off -diagonal parts are related by symmetry)

C       Since the integrals contains a factor of -1 no -DDOT here.

        EMU2XX=DDOT(NNBASX,WRK(KDCAO),1,WRK(KMAT),1)
        EMU2XY=DDOT(NNBASX,WRK(KDCAO),1,WRK(KMAT+NNBASX),1)
        EMU2XZ=DDOT(NNBASX,WRK(KDCAO),1,WRK(KMAT+2*NNBASX),1)
        EMU2YY=DDOT(NNBASX,WRK(KDCAO),1,WRK(KMAT+3*NNBASX),1)
        EMU2YZ=DDOT(NNBASX,WRK(KDCAO),1,WRK(KMAT+4*NNBASX),1)
        EMU2ZZ=DDOT(NNBASX,WRK(KDCAO),1,WRK(KMAT+5*NNBASX),1)

        WRK(KNSEL + I - 1) = EMU2XX+EMU2XY+EMU2XZ
     *                 + EMU2YY+EMU2YZ+EMU2ZZ
        EMUL2T = EMUL2T + WRK(KNSEL + I - 1)

C       Now the QM nuclear contribution

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
  301   CONTINUE

C       Remember the factor of 1/2 from the Taylor expansion
        ELOC   = 0.5D0*ELOC

        WRK(KNSNUC + I - 1) = ELOC
C        ELOCT = ELOCT + WRK(KNSNUC + I - 1)
        GO TO 20
      ENDIF


C     Send data to master (Do you love me?)
      CALL MPI_REDUCE(WRK(KNSNUC),WRK(KNSNUC),MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KNSEL),WRK(KNSEL),MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      CALL MPI_REDUCE(WRK(KTAO),WRK(KTAO),MMCENT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      CALL QEXIT('QMMM_QUAD_S')

      RETURN
      END
C**************************************************************************************
C  /* Deck qmmm_octupole_m */
      SUBROUTINE QMMM_OCTUPOLE_M()

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>


      CALL QENTER('QMMM_OCTUPOLE_M')


      CALL QEXIT('QMMM_OCTUPOLE_M')

      RETURN
      END
C**************************************************************************************
C  /* Deck qmmm_polari_m */
      SUBROUTINE QMMM_POLARI_M()

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <codata.h>


      CALL QENTER('QMMM_POLARI_M')


      CALL QEXIT('QMMM_POLARI_M')

      RETURN
      END
