C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
C
C      CBN+JK
C*************************************************************************
C  /* Deck qm3fck */
      SUBROUTINE QM3FCK(DCAO,DVAO,FSOL,EQM3T,ENSEL,EPOL,EELEL,
     *                  WRK,LWRK,IPRINT)

#include <implicit.h>
#include <dummy.h>
#include <inftap.h>
#include <priunit.h>
#include <qm3.h>
#include <mxcent.h>
#include <thrzer.h>
#include <iratdef.h>
#include <codata.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infpar.h>

      DIMENSION DCAO(*), DVAO(*)
      DIMENSION FSOL(*), WRK(LWRK)
      LOGICAL HSFOCK
      DIMENSION DFTNS(MXQM3)


      CALL QENTER('QM3FCK')
      KDTAO = 1
      KTAO = KDTAO + NNBASX
      KEND = KTAO + NNBASX
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3FCK',-KEND,LWRK)

      HSFOCK=((DODFT .AND. NASHT .GE. 1) .OR. HSROHF)

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

C     From COMMON
      ENSQM3 = 0.0D0
      EPOQM3 = 0.0D0

C     Calculate RRa, Obar and ENSA and write these to file
      CALL QM3_RRA(WRK(KDTAO),WRK(KEND),LWRK1,IPRINT) 

C     Calculate Ns and keep in memory (DFTNS).
C     (Run QM3_NSP instead of QM3_NS if #nodes > 1, Arnfinn nov. -08)
      IF (NODTOT.GE.1) THEN 
         IF (.NOT. INTDIR) THEN
            CALL QUIT('Include .INTDIR in *QM3 when run in parallel')
         END IF
         CALL QM3_NSP(WRK(KDTAO),DFTNS,WRK(KEND),LWRK1,IPRINT)
      ELSE
         CALL QM3_NS(WRK(KDTAO),DFTNS,WRK(KEND),LWRK1,IPRINT)
      END IF

C     Modify the fock operator. Modification returned in FSOL. 
      CALL QM3_FMO(FSOL,WRK(KEND),LWRK1,IPRINT)

C     Calculate the QM3 contribution to the energy (returned in EQM3T)
      CALL QM3_ENERGY(DFTNS,ENSEL,EPOL,EELEL,EQM3T,WRK(KEND),LWRK1,
     &                IPRINT)

      CALL QEXIT('QM3FCK')
      RETURN
      END
C ***********************************************************************
C  /* Deck qm3_rra */
      SUBROUTINE QM3_RRA(DCAO,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
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
#include "nuclei.h"


      DIMENSION WRK(LWRK)
      DIMENSION EELEC(3,MXQM3)
      DIMENSION FFJ(3)
      LOGICAL   LOINDM
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

C From here. 
C This is NOT good. In fact very bad. Needs to be fixed. 
C Neede for f.ex. spin-spin couplings since nucdep is used here for allocation.
C Can only be done AFTER coordinates for charges and polarizabilities have been defined 
C since nudep defines no. of coordinates !! Only active when requsted by the user.

      IF (REDCNT) NUCDEP=NUCIND 

C To here

      IF (LOSPC) RETURN

      CALL QENTER('QM3_RRA')

      KRRAOx = 1
      KRRAOy = KRRAOx + NNBASX
      KRRAOz = KRRAOy + NNBASX
      KRAx   = KRRAOz + NNBASX
      KRAy   = KRAx   + NCOMS
      KRAz   = KRAy   + NCOMS
      KWRK1  = KRAz   + NCOMS
      LWRK1  = LWRK   - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3_RRA',-KWRK1,LWRK)

      CALL DZERO(WRK(KRRAOx),3*NNBASX)
      CALL DZERO(WRK(KRAx),3*NCOMS)

      DO 879 I = 1, MXQM3
        DO 880 J = 1, 3
          EELEC(J,I) = 0.0D0
  880   CONTINUE
  879 CONTINUE

      IF (.NOT. INTDIR) THEN
        IF (IQM3PR .GT. 15) THEN
          WRITE(LUPRI,*) 'QM3RAINT: Read in integrals'
        ENDIF

        CALL GPOPEN(LUQM3E,'ELFDMM','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND(LUQM3E)
      END IF

      LM = 0

      IF (INTDIR) THEN
        L = NUSITE + NUALIS(0)
        OBKPX = DIPORG(1)
        OBKPY = DIPORG(2)
        OBKPZ = DIPORG(3)
      END IF


      DO 520 I = 1, ISYTP
        IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
          DO 521 J = NSYSBG(I), NSYSED(I)
            DO 522 K = 1, NUALIS(I)
              LM = LM + 1

              IF (INTDIR) THEN
                KMAT = KWRK1
                KLAST = KMAT + 3*NNBASX
                LWRK2 = LWRK - KLAST + 1
                IATNOW = NUCIND + L + LM

                KPATOM = 0
                NOSIMI = 3
                TOFILE = .FALSE.
                TRIMAT = .TRUE.
                EXP1VL = .FALSE.
                DIPORG(1) = CORD(1,IATNOW)
                DIPORG(2) = CORD(2,IATNOW)
                DIPORG(3) = CORD(3,IATNOW)

                RUNQM3 = .TRUE.
                CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIMI,WRK(KLAST),
     &                      LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                      KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)
                RUNQM3 = .FALSE.

                IF (IPRINT.GT.15) THEN
                  WRITE (LUPRI,'(/A)') ' Rra_ao_x matrix:'
                  CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

                  WRITE (LUPRI,'(/A)') ' Rra_ao_y matrix:'
                  CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

                  WRITE (LUPRI,'(/A)') ' Rra_ao_z matrix:'
                  CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
                END IF

                IF (QMDAMP) THEN
                  DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                   (DIPORG(2)-QMCOM(2))**2 +
     &                   (DIPORG(3)-QMCOM(3))**2
                  DIST = SQRT(DIST)
                  DFACT = (1-exp(-ADAMP*DIST))**3
                  CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
                ENDIF

                WRK(KRAx+LM-1) = DDOT(NNBASX,WRK(KMAT),1,DCAO,1)
                WRK(KRAy+LM-1) = DDOT(NNBASX,WRK(KMAT+NNBASX),1,DCAO,1)
                WRK(KRAz+LM-1) = DDOT(NNBASX,WRK(KMAT+2*NNBASX),
     *                                1,DCAO,1)
              ELSE   

                CALL READT(LUQM3E,NNBASX,WRK(KRRAOx))
                WRK(KRAx+LM-1) = DDOT(NNBASX,WRK(KRRAOx),1,DCAO,1)

                IF (IPRINT.GT.15) THEN
                  CALL AROUND('Rra_ao_x matrix:')
                  CALL OUTPAK(WRK(KRRAOx),NBAST,1,LUPRI)
                END IF 

                CALL READT(LUQM3E,NNBASX,WRK(KRRAOy))
                WRK(KRAy+LM-1) = DDOT(NNBASX,WRK(KRRAOy),1,DCAO,1)

                IF (IPRINT.GT.15) THEN
                  CALL AROUND('Rra_ao_y matrix:')
                  CALL OUTPAK(WRK(KRRAOy),NBAST,1,LUPRI)
                END IF

                CALL READT(LUQM3E,NNBASX,WRK(KRRAOz))
                WRK(KRAz+LM-1) = DDOT(NNBASX,WRK(KRRAOz),1,DCAO,1)

                IF (IPRINT.GT.15) THEN
                  CALL AROUND('Rra_ao_z matrix:')
                  CALL OUTPAK(WRK(KRRAOz),NBAST,1,LUPRI)
                END IF

              END IF
  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE

      IF (INTDIR) THEN
        DIPORG(1) = OBKPX
        DIPORG(2) = OBKPY
        DIPORG(3) = OBKPZ
      END IF


      IF (IQM3PR .GE. 5) THEN
        WRITE(LUPRI,'(/A/A/A)')
     *  ' +==========================================================+',
     *  ' | COM| <Rra>_x         | <Rra>_y         | <Rra>_z         |',
     *  ' +==========================================================+'

         NUM = 0

         DO 215 I = 1, ISYTP
           IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
             DO 216 J = NSYSBG(I), NSYSED(I)
               DO 217 K=1,NUALIS(I)

                 NUM = NUM + 1

                 WRITE(LUPRI,'(A,I3,A,3(F16.10,A)/A)')
     *  ' | ', J,'|', WRK(KRAx + NUM - 1),' |',
     *           WRK(KRAy + NUM - 1),' |', WRK(KRAz + NUM - 1),' |',
     *  ' +----------------------------------------------------------+'
  217          CONTINUE
  216        CONTINUE
           END IF
  215    CONTINUE
         WRITE(LUPRI,'(//,A)')
      END IF

      IF (LM .NE. NCOMS) THEN
        CALL QUIT('Error in no. of center of masses in QM3RAINT')
      END IF

      DO 534 LM = 1,NCOMS
        EELEC(1,LM) = WRK(KRAx+LM-1)
        EELEC(2,LM) = WRK(KRAy+LM-1)
        EELEC(3,LM) = WRK(KRAz+LM-1)
  534 CONTINUE

C     If RELMOM is true we want to include the external field(s) in
C     the determination of the induced dipole moments

        FFJ(1) = 0.0D0
        FFJ(2) = 0.0D0
        FFJ(3) = 0.0D0

      IF (RELMOM) THEN
        DO 330 IF =1, NFIELD
          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'ZDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
  330   CONTINUE
      END IF

      IF (FIXMOM) THEN
        WRITE(LUPRI,'(/A)')'FIXMOM: NO ITER. DETERM. OF MYIND'
      ELSE IF (LOSPC) THEN
        WRITE(LUPRI,'(/A)')'All MM models are SPC, INDMOM not called'
      ELSE
        LOINDM = .FALSE.
        CALL INDMOM(EELEC,LOINDM,FFJ)
      END IF

      CALL QM3_OBAR(FFJ)

      CALL CC_PUT31('CC_RA',NCOMS,WRK(KRAx),WRK(KRAy),WRK(KRAz))

      IF (.NOT. INTDIR) CALL GPCLOSE(LUQM3E,'KEEP')

      CALL QEXIT('QM3_RRA')
      RETURN
      END
C*********************************************************************
C  /* Deck qm3_ns */
      SUBROUTINE QM3_NS(DCAO,DFTNS,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
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
#include "nuclei.h"


      DIMENSION WRK(LWRK) 
      DIMENSION DFTNS(*), DCAO(NNBASX)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      
      CALL QENTER('QM3_NS')

C From here. 
C This is NOT good. In fact very bad. Needs to be fixed. 
C Neede for f.ex. spin-spin couplings since nucdep is used here for allocation.
C Can only be done AFTER coordinates for charges and polarizabilities have been defined 
C since nudep defines no. of coordinates !! Only active when requsted by the user.

      IF (REDCNT) NUCDEP=NUCIND

C To here


      KNSAO  = 1
      KWRK1  = KNSAO + NNBASX
      LWRK1  = LWRK  - KWRK1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3_NS',-KWRK1,LWRK)

      CALL DZERO(WRK(KNSAO),NNBASX)

      ENSEL = 0.0D0

      IF (.NOT. INTDIR) THEN
         CALL GPOPEN(LUQM3P,'POTMM','UNKNOWN',' ',
     &              'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQM3P)
      ENDIF

      FAC1 = -1.0D0

      L = 0

      IF (INTDIR) THEN
         OBKPX = DIPORG(1)
         OBKPY = DIPORG(2)
         OBKPZ = DIPORG(3)
      ENDIF

      DO 525 I = 1, ISYTP
         IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO 526 J = NSYSBG(I), NSYSED(I)
               DO 527 K = 1,NSISY(I)

                  L = L +1

                  IF (INTDIR) THEN

                     KMAT = KWRK1
                     KLAST = KMAT + NNBASX
                     LWRK2 = LWRK - KLAST + 1
                     IATNOW = NUCIND + L 

                     KPATOM = 0
                     NOSIMI = 1
                     TOFILE = .FALSE.
                     TRIMAT = .TRUE.
                     EXP1VL = .FALSE.
                     DIPORG(1) = CORD(1,IATNOW)
                     DIPORG(2) = CORD(2,IATNOW)
                     DIPORG(3) = CORD(3,IATNOW)

                     RUNQM3=.TRUE.
                     CALL GET1IN(WRK(KMAT),'NPETES ',NOSIMI,WRK(KLAST),
     &                        LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                        KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)
                     RUNQM3=.FALSE.
                     IF (OLDTG) THEN
                        FACQ = -1.0D0*CHAOLD(IATNOW)
                     ELSE
                        FACQ = -1.0D0*CHARGE(IATNOW)
                     ENDIF

                     CALL DSCAL(NNBASX,FACQ,WRK(KMAT),1)

                     IF (IPRINT.GT.15) THEN
                        WRITE (LUPRI,'(/A)') 'NUCPOT matrix'
                        CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
                     ENDIF 

                     DFTNS(L) = DDOT(NNBASX,DCAO,1,WRK(KMAT),1)

                  ELSE
                     CALL READT(LUQM3P,NNBASX,WRK(KNSAO))
                     IF (IPRINT.GT.15) THEN
                        WRITE (LUPRI,'(/A,I3,A)')
     *                       ' N(',L,')_ao matrix: '
                        CALL OUTPAK(WRK(KNSAO),NBAST,1,LUPRI)
                     ENDIF 

                     DFTNS(L) = DDOT(NNBASX,DCAO,1,WRK(KNSAO),1)
                     ENSEL     = ENSEL - DFTNS(L) 
                  ENDIF
 527           CONTINUE
 526        CONTINUE
         END IF
 525  CONTINUE       

      IF (INTDIR) THEN
         DIPORG(1) = OBKPX
         DIPORG(2) = OBKPY
         DIPORG(3) = OBKPZ
      ENDIF

      IF (.NOT. INTDIR) CALL GPCLOSE(LUQM3P,'KEEP')

C-------------------
C     Print section:
C-------------------
C
      IF (IQM3PR .GT. 3) THEN
         WRITE(LUPRI,'(/A/A/A)')
     *        ' +======================+',
     *        ' |Site| <N_s>           |',
     *        ' +======================+'

         LS = 0

         DO 215 I = 1, ISYTP
            IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
               DO 216 J = NSYSBG(I), NSYSED(I)
                  DO 217 K = 1, NSISY(I)

                     LS = LS + 1

                     WRITE(LUPRI,'(A,I3,A,F16.10,A/A)')
     *                      ' | ', LS,'|', DFTNS(LS),' |',
     *                      ' +----------------------+'
 217              CONTINUE
 216           CONTINUE
            END IF
 215     CONTINUE
         WRITE(LUPRI,'()')
      END IF
      CALL QEXIT('QM3_NS')
      RETURN
      END


C*********************************************************************
C  /* Deck qm3_nsp */
      SUBROUTINE QM3_NSP(DCAO,DFTNS,WRK,LWRK,IPRINT)
C
C     The parallel version of the routine QM3_NS, Arnfinn nov. -08
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
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
#include "nuclei.h"
#include "mtags.h"
#include "infpar.h"

      DIMENSION WRK(LWRK) 
      DIMENSION DFTNS(MXQM3), DCAO(NNBASX)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

#include <ibtfun.h>
      IPRTYP = 16
      
      CALL QENTER('QM3_NSP')

C From here.  TODO TODO TODO
C This is NOT good. In fact very bad. Needs to be fixed. 
C Needed for f.ex. spin-spin couplings since nucdep is used here for allocation.
C Can only be done AFTER coordinates for charges and polarizabilities have been defined 
C since nudep defines no. of coordinates !! Only active when requsted by the user.

      IF (REDCNT) NUCDEP=NUCIND

C To here

CHJAaJ-b KDFTNS not used /June 09
C     KDFTNS = 1
C     KWRK1  = KDFTNS + MXQM3
CHJAaJ-e
      KWRK1  = 1
      LWRK1  = LWRK   - KWRK1 + 1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3_NSP',-KWRK1,LWRK)

      ENSEL = 0.0D0

      FAC1 = -1.0D0

      L = 0

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      DO 525 I = 1, ISYTP
         IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            CALL QM3_NSPM(IPRTYP,DCAO,DFTNS,WRK(KWRK1),
     &                    LWRK1,IPRINT,EXP1VL,.FALSE.,ENSEL,FAC1,L,I)

         END IF
 525  CONTINUE       

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ



C-------------------
C     Print section:
C-------------------
C
      IF (IQM3PR .GT. 3) THEN
         WRITE(LUPRI,'(/A/A/A)')
     *        ' +======================+',
     *        ' |Site| <N_s>           |',
     *        ' +======================+'

         LS = 0
         
         DO 215 I = 1, ISYTP
            IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
               DO 216 J = NSYSBG(I), NSYSED(I)
                  DO 217 K = 1, NSISY(I)
                     
                     LS = LS + 1

                     WRITE(LUPRI,'(A,I3,A,F16.10,A/A)')
     *                      ' | ', LS,'|', DFTNS(LS),' |',
     *                      ' +----------------------+'
 217              CONTINUE
 216           CONTINUE
            END IF
 215     CONTINUE
         WRITE(LUPRI,'()')
      END IF
      CALL QEXIT('QM3_NSP')
      RETURN
      END
C*********************************************************************
      SUBROUTINE QM3_NSPM(IPRTYP,DCAO,DFTNS,WRK,LWRK,IPRINT,
     &                     EXP1VL,TOFILE,ENSEL,FAC1,LNUM,INUM)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
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
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include <mpif.h>
#endif

#include <cbiher.h>

C     The master routine

      LOGICAL TOFILE, EXP1VL
      DIMENSION WRK(LWRK)
      DIMENSION DCAO(NNBASX)
      DIMENSION DFTNS(MXQM3)

      IF (TOFILE) CALL QUIT('Parallel calculations do not allow '//
     &     'for storing NS-integrals on disk')

      KDCAO=1
      KDFTNS=KDCAO+NNBASX
      KLAST=KDFTNS+MXQM3
      LWRK1=LWRK-KLAST+1
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
C
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(DCAO,NNBASX,'DOUBLE',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(ENSEL,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(FAC1,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      CALL DZERO(DFTNS,MXQM3)

C     Loop over all MM nuclei

C      ENSEL = 0.0D0
C      FAC1 = -1.0D0

      DO J = NSYSBG(INUM), NSYSED(INUM)
         DO K = 1, NSISY(INUM)
            LNUM = LNUM + 1
            IWHO = -1
            CALL MPIXRECV(NWHO, 1, 'INTEGER', IWHO, MPTAG1)
            CALL MPIXSEND(LNUM,1,'INTEGER',NWHO,MPTAG2)
         END DO
      END DO

C     Send end message to all slaves

      LEND = -1
      DO ISLAVE = 1, NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MPTAG1)
         CALL MPIXSEND(LEND ,1,'INTEGER',NWHO,MPTAG2)
      END DO

C     Collect data from all slaves      
      CALL DZERO(WRK(KDFTNS),MXQM3)
      CALL MPI_REDUCE(WRK(KDFTNS),DFTNS,MXQM3,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      RETURN
      END

C*********************************************************************

      SUBROUTINE QM3_NSPS(WRK,LWRK,IPRTMP)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <gnrinf.h>
#include <qm3.h>
#include <mxcent.h>
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
#include "nuclei.h"
#include "infpar.h"
#include "mtags.h"
#if defined(VAR_MPI)
#include <mpif.h>
#endif

#include <cbiher.h>

C     The slave routine

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION WRK(LWRK)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      IQM3PR = IPRTMP
      QM3 = .TRUE.

      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      KDCAO  = 1
      KDFTNS = KDCAO  + NNBASX
      KLAST1 = KDFTNS + MXQM3
      LWRK1  = LWRK - KLAST1 + 1
      CALL MPIXBCAST(WRK(KDCAO),NNBASX,'DOUBLE',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(ENSEL,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(FAC1,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)

      KMAT   = KLAST1
      KLAST2 = KMAT   + NNBASX
      LWRK2  = LWRK1  - KLAST2 + KLAST1

      CALL DZERO(WRK(KDFTNS),MXQM3)

C Run loop over MM nuclear charges

 20   CONTINUE

      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MPTAG1)
      CALL MPIXRECV(L,1,'INTEGER',MASTER,MPTAG2)
      CALL DZERO(WRK(KMAT),NNBASX)
      IF (L .GT. 0) THEN
         IATNOW = NUCIND + L
         KPATOM = 0
         NOSIMI = 1
         TOFILE = .FALSE.
         TRIMAT = .TRUE.
         EXP1VL = .FALSE.
         DIPORG(1) = CORD(1,IATNOW)
         DIPORG(2) = CORD(2,IATNOW)
         DIPORG(3) = CORD(3,IATNOW)

         RUNQM3=.TRUE.
         CALL GET1IN(WRK(KMAT),'NPETES ',NOSIMI,WRK(KLAST2),
     &               LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &               KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)

         RUNQM3=.FALSE.
         IF (OLDTG) THEN
            FACQ = -1.0D0*CHAOLD(IATNOW)
         ELSE
            FACQ = -1.0D0*CHARGE(IATNOW)
         ENDIF
         CALL DSCAL(NNBASX,FACQ,WRK(KMAT),1)
         WRK(KDFTNS + L - 1) = DDOT(NNBASX,WRK(KDCAO),1,WRK(KMAT),1)
         ENSEL     = ENSEL - WRK(KDFTNS + L - 1)
         GO TO 20
      END IF

C     No more Ns to calculate

      CALL MPI_REDUCE(WRK(KDFTNS),WRK(KDFTNS),MXQM3,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      RETURN
      END

C*********************************************************************  
C  /* Deck qm3_fmo */
      SUBROUTINE QM3_FMO(FSOL,WRK,LWRK,IPRINT)   

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
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
#include "nuclei.h"


      DIMENSION WRK(LWRK)
      DIMENSION FSOL(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)


      IF (LOSPC .AND. (.NOT. OLDTG)) RETURN

      CALL QENTER('QM3_FMO')

      KRRAOx = 1
      KRRAOy = KRRAOx + NNBASX
      KRRAOz = KRRAOy + NNBASX
      KRAx   = KRRAOz + NNBASX
      KRAy   = KRAx   + NCOMS
      KRAz   = KRAy   + NCOMS
      KENSAx = KRAz   + NCOMS
      KENSAy = KENSAx + NCOMS
      KENSAz = KENSAy + NCOMS
      KTAO   = KENSAz + NCOMS
      KWRK1  = KTAO   + NNBASX
      LWRK1  = LWRK   - KWRK1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3_FMO',-KWRK1,LWRK)

      CALL DZERO(WRK(KRRAOx),3*NNBASX)
      CALL DZERO(WRK(KRAx),6*NCOMS)
      CALL DZERO(WRK(KTAO),NNBASX)

      IF (.NOT. LOSPC) THEN  

        CALL CC_GET31('CC_RA',NCOMS,
     *                 WRK(KRAx),WRK(KRAy),WRK(KRAz))

        CALL CC_GET31('ENSAFILE',NCOMS,
     *                 WRK(KENSAx),WRK(KENSAy),WRK(KENSAz))

        IF (IQM3PR .GT. 5) THEN
          WRITE(LUPRI,'(/A/A/A)')
     *  ' +==========================================================+',
     *  ' | COM| <Rra>_x         | <Rra>_y         | <Rra>_z         |',
     *  ' +==========================================================+'

          NUM = 0

          DO 215 I = 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
              DO 216 J = NSYSBG(I), NSYSED(I)
                DO 217 K=1,NUALIS(I)

                  NUM = NUM + 1

                  WRITE(LUPRI,'(A,I3,A,F16.10,A,F16.10,A,F16.10,A/A)')
     *  ' | ', J,'|', WRK(KRAx + NUM - 1),' |',
     *              WRK(KRAy + NUM - 1),' |', WRK(KRAz + NUM - 1),' |',
     *  '+----------------------------------------------------------+'
  217           CONTINUE
  216         CONTINUE
            END IF
  215     CONTINUE

          WRITE(LUPRI,'(/A/A/A)')
     *  ' +==========================================================+',
     *  ' | COM| ENSA_x          | ENSA_y          | ENSA_z          |',
     *  ' +==========================================================+'

          NUM = 0

          DO 415 I = 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN

              NUM = NUM + 1

              DO 416 J = NSYSBG(I), NSYSED(I)
                DO 417 K=1,NUALIS(I)
                  WRITE(LUPRI,'(A,I3,A,F16.10,A,
     &                          F16.10,A,F16.10,A/A)')
     *  ' | ', J,'|', WRK(KENSAx + NUM - 1),' |',
     *           WRK(KENSAy + NUM - 1),' |', WRK(KENSAz + NUM - 1),' |',
     *  ' +----------------------------------------------------------+'
  417           CONTINUE
  416         CONTINUE
            END IF
  415     CONTINUE
        END IF

        IF (.NOT. INTDIR) THEN
          CALL GPOPEN(LUQM3E,'ELFDMM','OLD',' ',
     &                'UNFORMATTED',IDUMMY,.FALSE.)
          REWIND(LUQM3E)
        ENDIF

        LM = 0

      IF (INTDIR) THEN
        L = NUSITE + NUALIS(0)
        OBKPX = DIPORG(1)
        OBKPY = DIPORG(2)
        OBKPZ = DIPORG(3)
      ENDIF

        DO 520 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 521 J = NSYSBG(I), NSYSED(I)
              DO 522 K = 1, NUALIS(I)
                LM = LM + 1

                IF (INTDIR) THEN
                  KMAT = KWRK1
                  KLAST = KMAT + 3*NNBASX
                  LWRK2 = LWRK - KLAST + 1
                  IATNOW = NUCIND + L + LM

                  KPATOM = 0
                  NOSIMI = 3
                  TOFILE = .FALSE.
                  TRIMAT = .TRUE.
                  EXP1VL = .FALSE.
                  DIPORG(1) = CORD(1,IATNOW)
                  DIPORG(2) = CORD(2,IATNOW)
                  DIPORG(3) = CORD(3,IATNOW)

                  RUNQM3=.TRUE.
                  CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIMI,WRK(KLAST),
     &                       LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                       KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)
                  RUNQM3=.FALSE.

                  IF (QMDAMP) THEN
                    DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                     (DIPORG(2)-QMCOM(2))**2 +
     &                     (DIPORG(3)-QMCOM(3))**2
                    DIST = SQRT(DIST)
                    DFACT = (1-exp(-ADAMP*DIST))**3
                    CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
                  ENDIF

                ELSE

                  CALL READT(LUQM3E,NNBASX,WRK(KRRAOx))
                  CALL READT(LUQM3E,NNBASX,WRK(KRRAOy))
                  CALL READT(LUQM3E,NNBASX,WRK(KRRAOz))

                ENDIF

                IF (MDLWRD(I) .EQ. 'SPC_E01') THEN
 
                  FACx =  - ALPIMM(I,K) * (WRK(KRAx + LM - 1)
     &                  + 0.5D0 * WRK(KENSAx + LM - 1))
                  FACy =  - ALPIMM(I,K) * (WRK(KRAy + LM - 1)
     &                  + 0.5D0 * WRK(KENSAy + LM - 1))
                  FACz =  - ALPIMM(I,K) * (WRK(KRAz + LM - 1)
     &                  + 0.5D0 * WRK(KENSAz + LM - 1))

                  IF (INTDIR) THEN
                    CALL DAXPY(NNBASX,FACx,WRK(KMAT),1,WRK(KTAO),1)
                    CALL DAXPY(NNBASX,FACy,WRK(KMAT+NNBASX),1,
     *                         WRK(KTAO),1)
                    CALL DAXPY(NNBASX,FACz,WRK(KMAT+2*NNBASX),1,
     *                         WRK(KTAO),1)
                  ELSE
                    CALL DAXPY(NNBASX,FACx,WRK(KRRAOx),1,WRK(KTAO),1)
                    CALL DAXPY(NNBASX,FACy,WRK(KRRAOy),1,WRK(KTAO),1)
                    CALL DAXPY(NNBASX,FACz,WRK(KRRAOz),1,WRK(KTAO),1)
                  ENDIF

                ENDIF
  522         CONTINUE
  521       CONTINUE
          END IF
  520   CONTINUE

        IF (.NOT. INTDIR) CALL GPCLOSE(LUQM3E,'KEEP')

        IF (IQM3PR .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' QM/MM contribution:'
            CALL OUTPAK(WRK(KTAO),NBAST,1,LUPRI)
        END IF

      END IF ! LOSPC

      IF (OLDTG) THEN

C       We use RRAOx in this case for the pot.energy integrals
        CALL DZERO(WRK(KRRAOx),NNBASX)

        IF (.NOT. INTDIR) THEN
          CALL GPOPEN(LUQM3P,'POTMM','UNKNOWN',' ',
     &              'UNFORMATTED',IDUMMY,.FALSE.)
          REWIND (LUQM3P)
        ENDIF

        FAC1 = -1.0D0

        L = 0

        IF (INTDIR) THEN
          OBKPX = DIPORG(1)
          OBKPY = DIPORG(2)
          OBKPZ = DIPORG(3)
        ENDIF 

        DO 525 I = 1, ISYTP
          IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
            DO 526 J = NSYSBG(I), NSYSED(I)
              DO 527 K = 1,NSISY(I)

                L = L +1

                IF (INTDIR) THEN

                  KMAT = KWRK1
                  KLAST = KMAT + NNBASX
                  LWRK2 = LWRK - KLAST + 1

                  IATNOW = NUCIND + L

                  CALL DZERO(WRK(KMAT),NNBASX)

                  KPATOM = 0
                  NOSIMI = 1
                  TOFILE = .FALSE.
                  TRIMAT = .TRUE.
                  EXP1VL = .FALSE.
                  DIPORG(1) = CORD(1,IATNOW)
                  DIPORG(2) = CORD(2,IATNOW)
                  DIPORG(3) = CORD(3,IATNOW)

                  RUNQM3=.TRUE.
                  CALL GET1IN(WRK(KMAT),'NPETES ',NOSIMI,WRK(KLAST),
     &                        LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                        KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)
                  RUNQM3=.FALSE.

                  FACQ =  -1.0D0*CHAOLD(IATNOW)

                  CALL DSCAL(NNBASX,FACQ,WRK(KMAT),1)

                  IF (IPRINT.GT.15) THEN
                    WRITE (LUPRI,'(/A)') 'NUCPOT matrix in QM3_FMO'
                    CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
                  ENDIF

                  CALL DAXPY(NNBASX,FAC1,WRK(KMAT),
     *                       1,WRK(KTAO),1)

                ELSE   

                  CALL READT(LUQM3P,NNBASX,WRK(KRRAOx))

                  IF (IPRINT.GT.15) THEN
                    WRITE (LUPRI,'(/A,I3,A)')
     *                 ' N(',L,')_ao matrix: '
                    CALL OUTPAK(WRK(KRRAOx),NBAST,1,LUPRI)
                  ENDIF  

                  CALL DAXPY(NNBASX,FAC1,WRK(KRRAOx),
     *                       1,WRK(KTAO),1)

                ENDIF

  527         CONTINUE
  526       CONTINUE
          END IF
  525   CONTINUE

        IF (INTDIR) THEN
          DIPORG(1) = OBKPX
          DIPORG(2) = OBKPY
          DIPORG(3) = OBKPZ
        ENDIF

        IF (.NOT. INTDIR)  CALL GPCLOSE(LUQM3P,'KEEP')

      END IF

      CALL PKSYM1(WRK(KTAO),FSOL,NBAS,NSYM,+1)

      CALL QEXIT('QM3_FMO')
      RETURN
      END
C*********************************************************************
C  /* Deck qm3_energy */
      SUBROUTINE QM3_ENERGY(DFTNS,EEL,EPOL,EELEL,EQM3T,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
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

      DIMENSION WRK(LWRK)
      DIMENSION FFJ(3), DFTNS(*)

      CALL QENTER('QM3_ENERGY')
      IF ( .NOT. (LOSPC) ) THEN

        KOMMSx = 1
        KOMMSy = KOMMSx + NCOMS
        KOMMSz = KOMMSy + NCOMS
        KRAx   = KOMMSz + NCOMS
        KRAy   = KRAx   + NCOMS
        KRAz   = KRAy   + NCOMS
        KWRK1  = KRAz   + NCOMS
        LWRK1  = LWRK   - KWRK1

        IF (LWRK1.LT.0) CALL QUIT( 'Too little work in QM3_ENERGY')

        CALL DZERO(WRK(KOMMSx),6*NCOMS)

      END IF

C     First we see if any static fields are to be added

        FFJ(1) = 0.0D0
        FFJ(2) = 0.0D0
        FFJ(3) = 0.0D0

      IF (RELMOM) THEN
        DO 330 IF =1, NFIELD
          IF (LFIELD(IF) .EQ. 'XDIPLEN') FFJ(1) = FFJ(1) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'YDIPLEN') FFJ(2) = FFJ(2) + EFIELD(IF)
          IF (LFIELD(IF) .EQ. 'ZDIPLEN') FFJ(3) = FFJ(3) + EFIELD(IF)
  330   CONTINUE
      END IF

C-------------------------------------
C     Add up the energy contributions:
C     1) Epol:
C-------------------------------------

      EQM3T = 0.0D0

      IF (.NOT. LOSPC) THEN
        CALL CC_GET31('OBARFILE',NCOMS,
     *                WRK(KOMMSx),WRK(KOMMSy),WRK(KOMMSz))

        CALL CC_GET31('CC_RA',NCOMS,
     *                WRK(KRAx),WRK(KRAy),WRK(KRAz))

        KS = 0

        DO 110 I = 1, ISYTP
          IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            DO 120 J = NSYSBG(I), NSYSED(I)
              DO 121 K = 1, NUALIS(I)
                KS = KS+1
                RAT = 0.0D0
                RAT =  0.5D0 * ALPIMM(I,K) * WRK(KRAx + KS - 1) *
     &                (WRK(KRAx + KS - 1) + WRK(KOMMSx + KS - 1))
     &               + 0.5D0 * ALPIMM(I,K) * WRK(KRAy + KS - 1) *
     &                (WRK(KRAy + KS - 1) + WRK(KOMMSy + KS - 1))
     &               + 0.5D0 * ALPIMM(I,K) * WRK(KRAz + KS - 1) *
     &                (WRK(KRAz + KS - 1) + WRK(KOMMSz + KS - 1))
                EQM3T = EQM3T - RAT
  121         CONTINUE
  120       CONTINUE
          END IF
  110   CONTINUE

        TEMP1 = EQM3T

        CALL QM3_OTILDE(OTILDE,FFJ)

      ELSE

        OTILDE = 0.0D0
        TEMP1 = 0.0D0

      END IF

      EQM3T = EQM3T + OTILDE
      EPOL = EQM3T

C--------
C 2) Eel:
C--------
C
      ETEMP = 0.0D0
      L = 0

      DO 130 I = 1, ISYTP
        IF (MDLWRD(I)(1:3) .EQ. 'SPC') THEN
          DO 140 J = NSYSBG(I), NSYSED(I)
            DO 150 K = 1, NSISY(I)
              L = L + 1
              ETEMP = ETEMP - DFTNS(L)
  150       CONTINUE
  140     CONTINUE
        END IF
  130 CONTINUE

      EELEL = ETEMP

      EQM3T = EQM3T + ENUQM3 + EELEL

      EEL = ENUQM3 + EELEL

C-------------
C 3) E(QM/MM):
C-------------

      EQM3T = EQM3T + ECLVDW

C-----------------------------------------------
C   Testing the energy contributions one by one:
C-----------------------------------------------

      IF (IQM3PR .GT. 1) THEN
        WRITE(LUPRI,'(//A/A/A)')
     *       ' +====================================='
     *      //'======================================'
     *      //'==================================+',

     *       ' |<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< '
     *      //'The different energy contributions outlined'
     *      //' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>|',

     *       ' +====================================='
     *      //'======================================'
     *      //'==================================+'
        WRITE(LUPRI,'(A)')
     *       ' | Evdw                | E-nuclear           |'
     *      //' -Sum_s <N_s>        |'
     *      //' -alpha*Sum_a[    >  | O-tilde             |',

     *       ' | contribution        | contribution        |'
     *      //'                     |'
     *      //'<Rra>*{<Rra>+OmmS}]  | contribution        |',

     *       ' +-------------------------------------'
     *      //'--------------------------------------'
     *      //'----------------------------------+'
        WRITE(LUPRI,'(A,F20.16,A,F20.16,A,F20.16,A,
     &                F20.16,A,F20.16,A)')
     *       ' |',ECLVDW,' |',ENUQM3,' |',EELEL,' |',TEMP1,' |',
     *      OTILDE,' |'
        WRITE(LUPRI,'(A/)')
     *       ' +====================================='
     *      //'======================================'
     *      //'==================================+'
      END IF

      CALL QEXIT('QM3_ENERGY')
      RETURN
      END
C*************************************************************************
C  /* Deck qm3grad */
      SUBROUTINE QM3GRAD(CREF,CMO,INDXCI,G,EQM3T,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <mxcent.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
#include <inflin.h>

      DIMENSION CREF(*),CMO(*),INDXCI(*)
      DIMENSION G(*),WRK(LWRK)
      LOGICAL FIRST, FNDLAB
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, DCVAL = D2 )
      CHARACTER*8 STAR8, SOLVDI, EODATA
      SAVE        FIRST
      DATA        FIRST/.TRUE./, STAR8/'********'/
      DATA        SOLVDI/'SOLVDIAG'/, EODATA/'EODATA  '/

      CALL QENTER('QM3GRAD')

      KDIASH  = 1
      KGRDLM  = KDIASH  + NVAR
      KUCMO   = KGRDLM  + NVARH
      KFSOLAO = KUCMO   + NORBT*NBAST
      KFSOLMO = KFSOLAO + NNBASX
      KDIALM  = KFSOLMO + NNORBX
      KDV     = KDIALM  + NVAR
      KDENS   = KDV     + N2BASX
      KDVS    = KDENS   + NNBASX
      KWRK    = KDVS    + NNBASX
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3GRAD',-KWRK,LWRK)

C     Construct the ao density matrix from the molecular orbial
C     coefficients.

      CALL FCKDEN((NISHT.GT.0),.FALSE.,WRK(KDV),
     *            DUMMY,CMO,DUMMY,WRK(KWRK),LWRK1)

      CALL DGEFSP(NBAST,WRK(KDV),WRK(KDVS))
      CALL PKSYM1(WRK(KDVS),WRK(KDENS),NBAS,NSYM,1)
    

      CALL DZERO(WRK(KDIASH),NVAR)
      CALL DZERO(WRK(KGRDLM),NVARH)

      EQM3T = 0.0D0
      ENSEL = 0.0D0
      EPOL = 0.0D0
      EELEL = 0.0D0

C     Make effective solvent operator and transform to mo basis
      CALL QM3FCK(WRK(KDENS),WRK(KDENS),WRK(KFSOLAO),EQM3T,ENSEL,
     &            EPOL,EELEL,WRK(KWRK),LWRK1,IPRINT)
      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL UTHU(WRK(KFSOLAO),WRK(KFSOLMO),WRK(KUCMO),WRK(KWRK),
     &             NBAST,NORBT)
 
      IF (.NOT. OLDTG) EQM3T = EQM3T - EELEL 

C     Calculate gradients
      IF (NWOPT .GT. 0) THEN
         CALL SOLGO(D2,WRK(KDENS),WRK(KFSOLMO),WRK(KGRDLM+NCONF))
      ENDIF

      IF (NWOPPT .GT. 0) THEN
         CALL OR1DIA(WRK(KDENS),WRK(KFSOLMO),WRK(KDIALM+NCONST),
     &               IPRINT)
      END IF

      DO 420 I = 0,(NVAR-1)
         WRK(KDIASH+I) = WRK(KDIASH+I)
     &                 + WRK(KDIALM+I)
     &                 + WRK(KGRDLM+I) * WRK(KGRDLM+I)
  420    CONTINUE

      FAC=1.0D0
      CALL DAXPY(NVARH,FAC,WRK(KGRDLM),1,G,1)

      IF (LUIT2 .GT. 0) THEN
         NW4 = MAX(NWOPT,4)
         REWIND LUIT2
         IF (FNDLAB(EODATA,LUIT2)) BACKSPACE LUIT2
         WRITE (LUIT2) STAR8,STAR8,STAR8,SOLVDI
         IF (NWOPT .GT. 0) CALL WRITT(LUIT2,NW4,WRK(KDIASH+NCONF))
         WRITE (LUIT2) STAR8,STAR8,STAR8,EODATA
      END IF

      FIRST = .FALSE.
      CALL QEXIT('QM3GRAD')
      RETURN
      END
C
C*************************************************************************
C  /* Deck qm3lin */
      SUBROUTINE QM3LIN(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &                  SOVECS,ORBLIN,WRK,LWRK)
C
#include <implicit.h>
#include <inflin.h>
#include <inforb.h>
#include <dummy.h>
#include <priunit.h>
#include <infpri.h>

      DIMENSION BOVECS(*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION SOVECS(*),WRK(LWRK)
      LOGICAL   ORBLIN

      CALL QENTER('QM3LIN')

      KDV   = 1
      KDVS  = KDV   + N2BASX
      KDENS = KDVS  + NNBASX
      KWRK  = KDENS + NNBASX
      LWRK1 = LWRK  - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3LIN',-KWRK,LWRK)

C     Construct the ao density matrix from the molecular orbial
C     coefficients.
      CALL FCKDEN((NISHT.GT.0),.FALSE.,WRK(KDV),
     *            DUMMY,CMO,DUMMY,WRK(KWRK),LWRK1)

      CALL DGEFSP(NBAST,WRK(KDV),WRK(KDVS))
      CALL PKSYM1(WRK(KDVS),WRK(KDENS),NBAS,NSYM,1)

C
      IF ( NOSIM .GT.0 ) THEN
         IF ( .NOT.ORBLIN ) THEN
            NSVAR  = NVARPT
         ELSE
            NSVAR  = NWOPPT
         END IF


      CALL QM3HESS(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &             WRK(KDENS),SOVECS,NSVAR,WRK(KWRK),LWRK1)

      END IF
      CALL QEXIT('QM3LIN')
      RETURN
      END
C
C*************************************************************************
C  /* Deck qm3hess */
      SUBROUTINE QM3HESS(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &                   DV,SVEC,NSVEC,WRK,LWRK)
C
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <maxash.h>
#include <infind.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <orgcom.h>
#include <infinp.h>
#include "nuclei.h"
C
      DIMENSION BOVECS(NWOPPT,*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION DV(*),SVEC(NSVEC,*),WRK(LWRK)
      LOGICAL   FULHES
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C
      CALL QENTER('QM3HESS')

      CALL QUIT('QM3HESS needs further debuging. Sorry !')
C
      KRXYO   = 1
      KUCMO   = KRXYO   + NOSIM*NNORBX
      KUBO    = KUCMO   + NORBT*NBAST
      KFSOLAO = KUBO    + NOSIM*N2ORBX
      KFSOLMO = KFSOLAO + NNBASX
      KUFSOL  = KFSOLMO + NNORBX
      KTGX    = KUFSOL  + N2ORBX
      KRRAOx  = KTGX    + NNORBX
      KRRAOy  = KRRAOx  + NNBASX
      KRRAOz  = KRRAOy  + NNBASX
      KTRAO   = KRRAOz  + NNBASX
      KTRMO   = KTRAO   + NNBASX
      KUTRX   = KTRMO   + NNORBX
      KUTR    = KUTRX   + N2ORBX
      KTRX    = KUTR    + N2ORBX
      KWRK    = KTRX    + NNORBX
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3HESS',-KWRK,LWRK)

      IF (JWOPSY .NE. LSYMPT .OR. NWOPPT .NE. NWOPT) THEN
         WRITE (LUPRI,*) 'QM3HESS ERROR: JWOPSY .ne. LSYMPT .or.'//
     *      ' NWOPPT .ne. NWOPT'
         WRITE (LUPRI,*) 'LSYMPT,JWOPSY:',LSYMPT,JWOPSY
         WRITE (LUPRI,*) 'NWOPPT,NWOPT :',NWOPPT,NWOPT
         CALL QUIT('QM3HESS ERROR,lsympt.ne.jwopsy or nwoppt.ne.nwopt')
      END IF

C     Determine if full Hessian or only orbital Hessian

      FULHES = (NSVEC .EQ. NVARPT)
      IF (FULHES) THEN
         JSOVEC = 1 + NCONST
      ELSE
         JSOVEC = 1
      END IF

      IF (IQM3PR .GT. 15) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(orb) on entry'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *               NSVEC,NOSIM,1,LUPRI)
      END IF

C     1. Construct effective operators Tgxo and RRayo

      CALL DZERO(WRK(KRXYO),NOSIM*NNORBX )
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)

C     1.a. Unpack symmetry blocked CMO
      CALL UPKCMO(CMO,WRK(KUCMO))
C     1.b. Calculate unpacked orbital trial vectors in UBO

      IF (NOSIM.GT.0) THEN
         DO 210 IOSIM = 1,NOSIM
            JUBO = KUBO + (IOSIM-1)*N2ORBX
            CALL UPKWOP(NWOPPT,JWOP,BOVECS(1,IOSIM),WRK(JUBO))
  210    CONTINUE
      END IF


      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)

C     1.c. Construct the Tg operator in ao and transform to mo
      CALL QM3FCKMO(CMO,WRK(KFSOLMO),WRK(KWRK),LWRK1,IPRINT)

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

      DO 200 IOSIM = 1,NOSIM
         JRXYO = KRXYO + (IOSIM-1)*NNORBX
         JUBO = KUBO + (IOSIM-1)*N2ORBX

C        1.d. Calculate one-index transformed Tg integrals

         CALL DZERO(WRK(KUTRX),N2ORBX)
         CALL TR1UH1(WRK(JUBO),WRK(KUFSOL),WRK(KUTRX),1)
         CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTGX))

         FAC = 1.0D0
         CALL DAXPY(NNORBX,FAC,WRK(KTGX),1,WRK(JRXYO),1)


         WRITE (LUPRI,'(/A,I3,A,I3)')
     *      ' >>> SOLLNO - Rxo matrix no.',IOSIM,' of',NOSIM
         CALL OUTPAK(WRK(JRXYO),NORBT,1,LUPRI)


         IF (.NOT. INTDIR) THEN
           CALL GPOPEN(LUQM3E,'ELFDMM','OLD',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
           REWIND(LUQM3E)
         ENDIF

C        1.e. Readin Rra integrals and transform to mo

         LM = 0

         IF (INTDIR) THEN
           L = NUSITE + NUALIS(0)
           OBKPX = DIPORG(1)
           OBKPY = DIPORG(2)
           OBKPZ = DIPORG(3)
         ENDIF

         DO 520 I = 1, ISYTP
           IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
             DO 521 J = NSYSBG(I), NSYSED(I)
               DO 522 K = 1, NUALIS(I)
                 LM = LM + 1

                 FAC = -ALPIMM(I,K)

C                Calculate one-index transformed Rra.x integrals

                 IF (INTDIR) THEN
                   KMAT = KWRK
                   KLAST = KMAT + 3*NNBASX
                   LWRK2 = LWRK - KLAST + 1
                   IATNOW = NUCIND + L + LM

                   CALL DZERO(WRK(KMAT),3*NNBASX)

                   KPATOM = 0
                   NOSIMI = 3
                   TOFILE = .FALSE.
                   TRIMAT = .TRUE.
                   EXP1VL = .FALSE.
                   DIPORG(1) = CORD(1,IATNOW)
                   DIPORG(2) = CORD(2,IATNOW)
                   DIPORG(3) = CORD(3,IATNOW)

                   RUNQM3=.TRUE.
                   CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIMI,WRK(KLAST),
     &                         LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &                         KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IQM3PR)
                   RUNQM3=.FALSE.

                   IF (QMDAMP) THEN
                     DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                      (DIPORG(2)-QMCOM(2))**2 +
     &                      (DIPORG(3)-QMCOM(3))**2
                     DIST = SQRT(DIST)
                     DFACT = (1-exp(-ADAMP*DIST))**3
                     CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
                   ENDIF

                   CALL DZERO(WRK(KTRMO),NNORBX)
                   CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &                       NBAST,NORBT)
                 ELSE
                   CALL DZERO(WRK(KRRAOx),NNBASX)
                   CALL READT(LUQM3E,NNBASX,WRK(KRRAOx))
                   CALL DZERO(WRK(KTRMO),NNORBX)
                   CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                 ENDIF

                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)
                 CALL DZERO(WRK(KTRX),NNORBX)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL TR1UH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX),1)
                 CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTRX))

                 FACx = 0.0D0
                 DO 523 IW = 1,NISHT
                   IG = ISX(IW)
                   FACx = FACx + 2.0D0*WRK(KTRX -1 + IROW(IG+1))
  523            CONTINUE
                 FACx = FACx * FAC
  
                 CALL DAXPY(NNORBX,FACx,WRK(KTRMO),1,WRK(JRXYO),1)

C                Calculate one-index transformed Rra.y integrals

                 IF (INTDIR) THEN
                   CALL DZERO(WRK(KTRMO),NNORBX)
                   CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KLAST),NBAST,NORBT)
                 ELSE
                   CALL DZERO(WRK(KRRAOy),NNBASX)
                   CALL READT(LUQM3E,NNBASX,WRK(KRRAOy))
                   CALL DZERO(WRK(KTRMO),NNORBX)
                   CALL UTHU(WRK(KRRAOy),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                 ENDIF

                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)
                 CALL DZERO(WRK(KTRX),NNORBX)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL TR1UH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX),1)
                 CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTRX))

                 FACy = 0.0D0
                 DO 524 IW = 1,NISHT
                   IG = ISX(IW)
                   FACy = FACy + 2.0D0*WRK(KTRX -1 + IROW(IG+1))
  524            CONTINUE
                 FACy = FACy * FAC

                 CALL DAXPY(NNORBX,FACy,WRK(KTRMO),1,WRK(JRXYO),1)

C                Calculate one-index transformed Rra.z integrals

                 IF (INTDIR) THEN
                   CALL DZERO(WRK(KTRMO),NNORBX)
                   CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KLAST),NBAST,NORBT)
                 ELSE
                   CALL DZERO(WRK(KRRAOz),NNBASX)
                   CALL READT(LUQM3E,NNBASX,WRK(KRRAOz))
                   CALL DZERO(WRK(KTRMO),NNORBX)
                   CALL UTHU(WRK(KRRAOz),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                 ENDIF

                 CALL DZERO(WRK(KUTR),N2ORBX)
                 CALL DZERO(WRK(KUTRX),N2ORBX)
                 CALL DZERO(WRK(KTRX),NNORBX)
                 CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                 CALL TR1UH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX),1)
                 CALL DGETSP(NORBT,WRK(KUTRX),WRK(KTRX))

                 FACz = 0.0D0
                 DO 525 IW = 1,NISHT
                   IG = ISX(IW)
                   FACz = FACz + 2.0D0*WRK(KTRX -1 + IROW(IG+1))
  525            CONTINUE
                 FACz = FACz * FAC

                 CALL DAXPY(NNORBX,FACz,WRK(KTRMO),1,WRK(JRXYO),1)

  522          CONTINUE
  521        CONTINUE
           END IF
  520    CONTINUE
C
         IF (.NOT. INTDIR) CALL GPCLOSE(LUQM3E,'KEEP')

         IF (INTDIR) THEN
           DIPORG(1) = OBKPX
           DIPORG(2) = OBKPY
           DIPORG(3) = OBKPZ
         END IF

  200 CONTINUE

      IF (IQM3PR .GT. 15) THEN
         DO 600 IOSIM = 1,NOSIM
            JRXYO = KRXYO + (IOSIM-1)*NNORBX
            WRITE (LUPRI,'(/A,I3,A,I3)')
     *         ' >>> SOLLNO - (Rxo + Ryo) matrix no.',IOSIM,' of',NOSIM
            CALL OUTPAK(WRK(JRXYO),NORBT,1,LUPRI)
  600    CONTINUE
      END IF

C     2. Add effective operators Tgxo and RRayo contribution to SVEC(NSVEC,NOSIM) 
C     Tell SOLGO only to use the NWOPPT first JWOP entries
      DO 800 IOSIM = 1,NOSIM
         JRXYO  = KRXYO  + (IOSIM-1)*NNORBX
         CALL SOLGO(D2,DV,WRK(JRXYO),SVEC(JSOVEC,IOSIM))
  800 CONTINUE
C
      IF (IQM3PR .GT. 15) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(orb) on exit'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *               NSVEC,NOSIM,1,LUPRI)
      END IF

      CALL QEXIT('QM3HESS')
      RETURN
      END
C*************************************************************************
C  /* Deck qm3fck */
      SUBROUTINE QM3FCKMO(CMO,FSOL,WRK,LWRK,IPRINT)
C
C     Construct the QM/MM contribution to the Fock-matrix in MO basis
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <qm3.h>
#include <inforb.h>
#include <infopt.h>
C
      DIMENSION CMO(*), FSOL(*), WRK(LWRK)
C
      CALL QENTER('QM3FCKMO')
C
      KDV     = 1
      KDENS   = KDV     + N2BASX
      KDVS    = KDENS   + NNBASX
      KFSOLAO = KDVS    + NNBASX
      KUCMO   = KFSOLAO + NNBASX
      KWRK    = KUCMO   + NORBT*NBAST
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3FCKMO',-KWRK,LWRK)

C     Construct the density matrix
      CALL FCKDEN((NISHT.GT.0),.FALSE.,WRK(KDV),
     *            DUMMY,CMO,DUMMY,WRK(KWRK),LWRK1)

      CALL DGEFSP(NBAST,WRK(KDV),WRK(KDVS))
      CALL PKSYM1(WRK(KDVS),WRK(KDENS),NBAS,NSYM,1)

C     Construct the QM/MM contribution to the Fock-matrix in AO
      CALL QM3FCK(WRK(KDENS),WRK(KDENS),WRK(KFSOLAO),ESOLT,
     *            ENSQM3,EPOQM3,EELEL,WRK(KWRK),LWRK1,IPRINT)

C     Transform to mo
      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL UTHU(WRK(KFSOLAO),FSOL,WRK(KUCMO),WRK(KWRK),
     &             NBAST,NORBT)
C
      CALL QEXIT('QM3FCKMO')
      RETURN
      END
*******************************************************************************
C  /* Deck qmmmfck */
      SUBROUTINE QMMMFCK(DCAO,DVAO,FSOL,EQMMM,WRK,LWRK,IPRINT)

#include <implicit.h>
#include <dummy.h>
#include <inftap.h>
#include <priunit.h>
#include <mxcent.h>
#include <qmmm.h>
#include <thrzer.h>
#include <iratdef.h>
#include <codata.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
#include <orgcom.h>
#include <infpar.h>


      DIMENSION DCAO(*), DVAO(*)
      DIMENSION FSOL(*), WRK(LWRK)
      LOGICAL HSFOCK


      CALL QENTER('QMMMFCK')
      KDTAO = 1
      KTAO = KDTAO + NNBASX
      KEND = KTAO + NNBASX
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMMFCK',-KEND,LWRK)

      HSFOCK=((DODFT .AND. NASHT .GE. 1) .OR. HSROHF)

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
C      IF (NODTOT .GE. 1) THEN
      IF (.TRUE.) THEN
         PRINT *,'arnfinn before call qmmm_fck_ao_p'
         CALL QMMM_FCK_AO_P(FSOL,WRK(KDTAO),EQMMM,WRK(KEND),LWRK1,
     & IPRINT)
      ELSE
         CALL QMMM_FCK_AO(FSOL,WRK(KDTAO),EQMMM,WRK(KEND),LWRK1,IPRINT)
      ENDIF
      CALL QEXIT('QMMMFCK')
      RETURN
      END

C**************************************************************************************
C  /* Deck qmmm_fck_ao */
      SUBROUTINE QMMM_FCK_AO(FSOL,DCAO,ESOLT,WRK,LWRK,IPRINT)   

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


      DIMENSION WRK(LWRK)
      DIMENSION FSOL(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      LOGICAL EXCENT,LOCDEB,LZERO
      INTEGER NZERAL,RI,CI,CZA
      DOUBLE PRECISION INVAMT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0, D3 = 3.0D0, D6 = 6.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('QMMM_FCK_AO')

      LOCDEB = .FALSE.

      KTAO   = 1
      KNSEL  = KTAO   + NNBASX
      KNSNUC = KNSEL  + MMCENT
      KWRK1  = KNSNUC + MMCENT
      LWRK1  = LWRK   - KWRK1

      IF (LWRK1 .LT. 0) CALL ERRWRK('QMMM_FCK_AO 1',-KWRK1,LWRK)

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

        KMAT  = KWRK1
        KLAST = KMAT + NNBASX
        LWRK2 = LWRK - KLAST + 1

        IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO 2',-KLAST,LWRK)

        FAC1 = 1.0D0      
        EXPNST = 0.0D0
        ECHCH = 0.0D0

        DO 100 I = 1,MMCENT 

           IF (ABS(MUL0MM(I)) .LE. THRMM) THEN
             WRK(KNSEL+I-1)  = 0.0D0
             WRK(KNSNUC+I-1) = 0.0D0
             GOTO 100
           ENDIF 

           DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *             (MMCORD(2,I)-QMCOM(2))**2 +
     *             (MMCORD(3,I)-QMCOM(3))**2
           DIST = SQRT(DIST2)

           IF (DIST .GT. RCUTMM) THEN
             WRK(KNSEL+I-1)  = 0.0D0
             WRK(KNSNUC+I-1) = 0.0D0
             IF (LOCDEB) THEN
               WRITE(LUPRI,*) 'Skipping charge ', I
             ENDIF 
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

           IF (LOCDEB) THEN
C            Test for numerical int.
             CORZSAVE  = DIPORG(3)
             KMAT1  = KLAST
             KMAT2  = KMAT1 + NNBASX
             KMAT3  = KMAT2 + NNBASX
             KLAST1 = KMAT3 + NNBASX
             LWRK3  = LWRK - KLAST1 + 1

             IF (LWRK3 .LT. 0) CALL ERRWRK('QMMM_FCK_AO 3',-KLAST1,LWRK)

             CALL DZERO(WRK(KMAT1),3*NNBASX)
 
             DIPORG(3) = DIPORG(3) + 0.01
             RUNQM3=.TRUE.
             CALL GET1IN(WRK(KMAT1),'NPETES ',NOSIM,WRK(KLAST1),
     &                   LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
             DIPORG(3) = DIPORG(3) - 0.02
             CALL GET1IN(WRK(KMAT2),'NPETES ',NOSIM,WRK(KLAST1),
     &                   LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
             DIPORG(3) = CORZSAVE
             CALL GET1IN(WRK(KMAT3),'NPETES ',NOSIM,WRK(KLAST1),
     &                   LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                   KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
             RUNQM3=.FALSE.
C  Gradient        
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
     &                 LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                 KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
           RUNQM3=.FALSE.

           IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
             WRITE (LUPRI,'(/A)') 'Pot. energy matrix in QMMM_FCK_AO'
             CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
           ENDIF

           CALL DSCAL(NNBASX,MUL0MM(I),WRK(KMAT),1)
           EXPNS=DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
           WRK(KNSEL+I-1) = EXPNS

           CALL DAXPY(NNBASX,FAC1,WRK(KMAT),1,WRK(KTAO),1)
           EXPNST=EXPNST+EXPNS

C          Now the QM nuclear contribution

           ECHCHL  = 0.0D0
           DO 101 J = 1,NUCIND
             XDIS   = CORD(1,J) - MMCORD(1,I) 
             YDIS   = CORD(2,J) - MMCORD(2,I) 
             ZDIS   = CORD(3,J) - MMCORD(3,I) 
             DIST2  = XDIS**2+YDIS**2+ZDIS**2
             DIST   = SQRT(DIST2) 
             ECHCHL = ECHCHL + CHARGE(J)*MUL0MM(I)/DIST
  101      CONTINUE
           ECHCH = ECHCH + ECHCHL
           WRK(KNSNUC+I-1) = ECHCHL

  100   CONTINUE

        ECHART = EXPNST + ECHCH
        ESOLT  = ECHART

        IF ( (IPRINT.GT.5) .OR. (LOCDEB) ) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Center Charge-electronic Charge-nuclear Total'
         DO 102 I = 1,MMCENT
           ELTEMP = WRK(KNSEL+I-1) + WRK(KNSNUC+I-1)
           WRITE(LUPRI,*) I,WRK(KNSEL+I-1),WRK(KNSNUC+I-1),ELTEMP
  102    CONTINUE

         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Total '
         WRITE(LUPRI,*) EXPNST, ECHCH, EXPNST+ECHCH
         WRITE(LUPRI,*)
        ENDIF

      END IF ! NMULT .GE. 0

C     2) The dipole correction to the Fock/KS operator
C
      IF (NMULT .GE. 1) THEN

        PRINT *,'arnfinn bfr qmmmdipole esolt',ESOLT
        KMAT  = KWRK1
        KLAST = KMAT + 3*NNBASX
        LWRK2 = LWRK - KLAST + 1

        IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO 4',-KLAST,LWRK)

        FAC1   = 1.0D0      
        FACM1  = -1.0D0
        EMUL1T = 0.0D0
        ELOCT  = 0.0D0

        DO 200 I = 1,MMCENT 

C         See if the dipole moment at this site is zero
          DNORM2 = MUL1MM(1,I)**2+MUL1MM(2,I)**2+MUL1MM(3,I)**2     
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) THEN
C           we use KNS to collect the energy terms
            WRK(KNSEL+I-1)  = 0.0D0
            WRK(KNSNUC+I-1) = 0.0D0
            GOTO 200
          ENDIF 

          DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *            (MMCORD(2,I)-QMCOM(2))**2 +
     *            (MMCORD(3,I)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
            WRK(KNSEL+I-1)  = 0.0D0
            WRK(KNSNUC+I-1) = 0.0D0
            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping dipole ', I
            ENDIF
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

          IF (LOCDEB) THEN
C           TEST for numerical int.
            CORZSAVE  = DIPORG(3)

            KMAT1  = KLAST
            KMAT2  = KMAT1 + 3*NNBASX
            KMAT3  = KMAT2 + 3*NNBASX
            KLAST1 = KMAT3 + 3*NNBASX
            LWRK3  = LWRK - KLAST1 + 1

            IF (LWRK3 .LT. 0) CALL ERRWRK('QMMM_FCK_AO 5',-KLAST1,LWRK)

            CALL DZERO(WRK(KMAT1),9*NNBASX)
 
            DIPORG(3) = DIPORG(3) + 0.01
            RUNQM3=.TRUE.
            CALL GET1IN(WRK(KMAT1),'NEFIELD',NOSIM,WRK(KLAST1),
     &                  LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                  KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
            DIPORG(3) = DIPORG(3) - 0.02
            CALL GET1IN(WRK(KMAT2),'NEFIELD',NOSIM,WRK(KLAST1),
     &                  LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                  KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
            DIPORG(3) = CORZSAVE
            CALL GET1IN(WRK(KMAT3),'NEFIELD',NOSIM,WRK(KLAST1),
     &                  LWRK3,LABINT,INTREP,INTADR,I,TOFILE,
     &                  KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
            RUNQM3=.FALSE.
C  Gradient
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
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
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
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
                IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
                ENDIF
  123         CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
              DIST = (DIPORG(1)-QMCOM(1))**2 +
     &               (DIPORG(2)-QMCOM(2))**2 +
     &               (DIPORG(3)-QMCOM(3))**2
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

          IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
            WRITE (LUPRI,'(/A)') ' E_x_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_y matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_z matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
          END IF

          CALL DSCAL(NNBASX,MUL1MM(1,I),WRK(KMAT),1)
          CALL DSCAL(NNBASX,MUL1MM(2,I),WRK(KMAT+NNBASX),1)
          CALL DSCAL(NNBASX,MUL1MM(3,I),WRK(KMAT+2*NNBASX),1)

          CALL DAXPY(NNBASX,FACM1,WRK(KMAT),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FACM1,WRK(KMAT+NNBASX),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FACM1,WRK(KMAT+2*NNBASX),1,WRK(KTAO),1)

          EXCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
          EYCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
          EZCOMP = -DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

          WRK(KNSEL+I-1) = EXCOMP + EYCOMP + EZCOMP
          EMUL1T = EMUL1T + WRK(KNSEL+I-1)

C         Now the QM nuclear contribution. Note that we switch the sign here by writing 
C         CORD(1,J) - MMCORD(1,I)

          ELOC     = 0.0D0
          DO 201 J = 1,NUCIND
            XDIS   = CORD(1,J) - MMCORD(1,I) 
            YDIS   = CORD(2,J) - MMCORD(2,I) 
            ZDIS   = CORD(3,J) - MMCORD(3,I) 
            DIST2  = XDIS**2+YDIS**2+ZDIS**2
            DIST   = SQRT(DIST2) 
            DIST3  = DIST2*DIST
            ELOC   = ELOC
     *             + CHARGE(J)*MUL1MM(1,I)*XDIS/DIST3
     *             + CHARGE(J)*MUL1MM(2,I)*YDIS/DIST3 
     *             + CHARGE(J)*MUL1MM(3,I)*ZDIS/DIST3
  201     CONTINUE
          WRK(KNSNUC+I-1) = ELOC
          ELOCT = ELOCT + WRK(KNSNUC+I-1)

  200   CONTINUE

        EDIPT = EMUL1T + ELOCT
        ESOLT = ESOLT + EDIPT

        IF ( (IPRINT.GT.5) .OR. (LOCDEB) ) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Center Dipole-electronic Dipole-nuclear Total'
         DO 202 I = 1,MMCENT
           ETEMP = WRK(KNSEL+I-1) + WRK(KNSNUC+I-1)
           WRITE(LUPRI,*) I,WRK(KNSEL+I-1),WRK(KNSNUC+I-1),ETEMP
  202    CONTINUE
  
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Total '
         WRITE(LUPRI,*) EMUL1T, ELOCT, EMUL1T+ELOCT
         WRITE(LUPRI,*) 
        ENDIF
        PRINT *,'arnfinn after qmmmdipole esolt',ESOLT
        PRINT *, EMUL1T, ELOCT, EMUL1T+ELOCT

      END IF ! (NMULT .GE. 1)

C     3) The quadrupole correction to the Fock/KS operator

      IF (NMULT .GE. 2) THEN

        KMAT  = KWRK1
        KLAST = KMAT + 6*NNBASX
        LWRK2 = LWRK - KLAST + 1

        IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO 6',-KLAST,LWRK)

        FAC1   = 1.0D0
        FACM1  = -1.0D0
        EMUL2T = 0.0D0
        ELOCT  = 0.0D0

        DO 300 I = 1,MMCENT

C         See if the quadrupole moment at this site is zero
          DNORM2 = MUL2MM(1,I)**2+MUL2MM(2,I)**2+MUL2MM(3,I)**2
     *           + MUL2MM(4,I)**2+MUL2MM(5,I)**2+MUL2MM(6,I)**2
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) THEN
            WRK(KNSEL+I-1)  = 0.0D0
            WRK(KNSNUC+I-1) = 0.0D0
C           we use KNS to collect the energy terms
            GOTO 300
          ENDIF

          DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *            (MMCORD(2,I)-QMCOM(2))**2 +
     *            (MMCORD(3,I)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
            WRK(KNSEL+I-1)  = 0.0D0
            WRK(KNSNUC+I-1) = 0.0D0
            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping quadrupole ', I
            ENDIF
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
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
          RUNQM3=.FALSE.

          IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
            WRITE (LUPRI,'(/A)') ' E_xx_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xy matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xz matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yy_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+3*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yz_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+4*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_zz_matrix in QMMM_FCK_AO:'
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
C         The integrals contains a factor of -1. Therefore daxpy with fac1
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT+NNBASX),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT+2*NNBASX),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT+3*NNBASX),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT+4*NNBASX),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT+5*NNBASX),1,WRK(KTAO),1)

C         Contract with the density to get the expectation values. 
C         The factor of 1/2 in the Taylor expansion has been included.
C         Also, the off-diagonal elements have been scaled by 2 in order 
C         to include all contributions (the off -diagonal parts are related by symmetry)

C         Since the integrals contains a factor of -1 no -DDOT here.

          EMU2XX=DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
          EMU2XY=DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
          EMU2XZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)
          EMU2YY=DDOT(NNBASX,DCAO,1,WRK(KMAT+3*NNBASX),1)
          EMU2YZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+4*NNBASX),1)
          EMU2ZZ=DDOT(NNBASX,DCAO,1,WRK(KMAT+5*NNBASX),1)

          WRK(KNSEL+I-1) = EMU2XX+EMU2XY+EMU2XZ
     *                   + EMU2YY+EMU2YZ+EMU2ZZ
          EMUL2T = EMUL2T + WRK(KNSEL+I-1)

C         Now the QM nuclear contribution

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
     *             +   CHARGE(J)*MUL2MM(1,I)*TXX
     *             + 2*CHARGE(J)*MUL2MM(2,I)*TXY
     *             + 2*CHARGE(J)*MUL2MM(3,I)*TXZ
     *             +   CHARGE(J)*MUL2MM(4,I)*TYY
     *             + 2*CHARGE(J)*MUL2MM(5,I)*TYZ
     *             +   CHARGE(J)*MUL2MM(6,I)*TZZ
  301     CONTINUE

C         Remember the factor of 1/2 from the Taylor expansion
          ELOC   = 0.5D0*ELOC

          WRK(KNSNUC+I-1) = ELOC
          ELOCT = ELOCT + WRK(KNSNUC+I-1)

  300   CONTINUE

        EQUADT = EMUL2T + ELOCT
        ESOLT  = ESOLT  + EQUADT

        IF ( (IPRINT.GT.5) .OR. (LOCDEB) ) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Center Quadr.-electronic Quadr.-nuclear Total'
         DO 302 I = 1,MMCENT
           ETEMP = WRK(KNSEL+I-1) + WRK(KNSNUC+I-1)
           WRITE(LUPRI,*) I,WRK(KNSEL+I-1),WRK(KNSNUC+I-1),ETEMP
  302    CONTINUE
  
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Total '
         WRITE(LUPRI,*) EMUL2T, ELOCT, EMUL2T+ELOCT
         WRITE(LUPRI,*) 
        ENDIF

      END IF ! (NMULT .GE. 2)

C     4) The octupole correction to the Fock/KS operator. We do the integrals numerically here.

      IF (NMULT .GE. 3) THEN

        KMAT1  = KWRK1
        KMAT2  = KMAT1 + 6*NNBASX
        KMAT3  = KMAT2 + 6*NNBASX
        KMAT4  = KMAT3 + 6*NNBASX
        KLAST  = KMAT4 + 6*NNBASX
        LWRK2  = LWRK  - KLAST + 1

        IF (LWRK2 .LT. 0) CALL ERRWRK('QMMM_FCK_AO 7',-KLAST,LWRK)

        FAC1   = 1.0D0
        FACM1  = -1.0D0
        EMUL3T = 0.0D0
        ELOCT  = 0.0D0

        DO 410 I = 1,MMCENT

C         See if the octupole moment at this site is zero
          DNORM2 = MUL3MM(1,I)**2+MUL3MM(2,I)**2+MUL3MM(3,I)**2
     *           + MUL3MM(4,I)**2+MUL3MM(5,I)**2+MUL3MM(6,I)**2
     *           + MUL3MM(7,I)**2+MUL3MM(8,I)**2+MUL3MM(9,I)**2
     *           + MUL3MM(10,I)**2
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) THEN
            WRK(KNSEL+I-1)  = 0.0D0
            WRK(KNSNUC+I-1) = 0.0D0
C           we use KNS to collect the energy terms
            GOTO 410
          ENDIF

          DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *            (MMCORD(2,I)-QMCOM(2))**2 +
     *            (MMCORD(3,I)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
            WRK(KNSEL+I-1)  = 0.0D0
            WRK(KNSNUC+I-1) = 0.0D0
            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping octupole ', I
            ENDIF
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

C Consider first the derivative wrt. z 

          DIPORG(3) = DIPORG(3) + DELFLD
          RUNQM3=.TRUE.
          CALL GET1IN(WRK(KMAT1),'ELFGRDC',NOSIM,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
          RUNQM3=.FALSE.
          DIPORG(3) = ZSAVE

          DIPORG(3) = DIPORG(3) - DELFLD
          RUNQM3=.TRUE.
          CALL GET1IN(WRK(KMAT2),'ELFGRDC',NOSIM,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
          RUNQM3=.FALSE.
          DIPORG(3) = ZSAVE

C  Calculate gradient of these integrals

           FSCAL = 1.0D0/(2.0D0*DELFLD)
           CALL DAXPY(6*NNBASX,FACM1,WRK(KMAT2),1,WRK(KMAT1),1)
           CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT1),1)
           CALL DZERO(WRK(KMAT2),6*NNBASX)

          IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
            WRITE (LUPRI,'(/A)') ' E_xxz_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT1),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xyz matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT1+NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xzz matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT1+2*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yyz_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT1+3*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yzz_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT1+4*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_zzz_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT1+5*NNBASX),NBAST,1,LUPRI)
          END IF

C Continue with the derivatives wrt. y   

          DIPORG(2) = DIPORG(2) + DELFLD
          RUNQM3=.TRUE.
          CALL GET1IN(WRK(KMAT2),'ELFGRDC',NOSIM,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
          RUNQM3=.FALSE.
          DIPORG(2) = YSAVE

          DIPORG(2) = DIPORG(2) - DELFLD
          RUNQM3=.TRUE.
          CALL GET1IN(WRK(KMAT3),'ELFGRDC',NOSIM,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
          RUNQM3=.FALSE.
          DIPORG(2) = YSAVE

C  Calculate gradient of these integrals

           FSCAL = 1.0D0/(2.0D0*DELFLD)
           CALL DAXPY(6*NNBASX,FACM1,WRK(KMAT3),1,WRK(KMAT2),1)
           CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT2),1)
           CALL DZERO(WRK(KMAT3),6*NNBASX)

          IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
            WRITE (LUPRI,'(/A)') ' E_xxy_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT2),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xyy matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT2+NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xzy matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT2+2*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yyy_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT2+3*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yzy_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT2+4*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_zzy_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT2+5*NNBASX),NBAST,1,LUPRI)
          END IF

C Continue with the derivatives wrt. x

          DIPORG(1) = DIPORG(1) + DELFLD
          RUNQM3=.TRUE.
          CALL GET1IN(WRK(KMAT3),'ELFGRDC',NOSIM,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
          RUNQM3=.FALSE.
          DIPORG(1) = XSAVE

          DIPORG(1) = DIPORG(1) - DELFLD
          RUNQM3=.TRUE.
          CALL GET1IN(WRK(KMAT4),'ELFGRDC',NOSIM,WRK(KLAST),
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
          RUNQM3=.FALSE.
          DIPORG(1) = XSAVE

C  Calculate gradient of these integrals

           FSCAL = 1.0D0/(2.0D0*DELFLD)
           CALL DAXPY(6*NNBASX,FACM1,WRK(KMAT4),1,WRK(KMAT3),1)
           CALL DSCAL(3*NNBASX,FSCAL,WRK(KMAT3),1)
           CALL DZERO(WRK(KMAT4),6*NNBASX)

          IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
            WRITE (LUPRI,'(/A)') ' E_xxx_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT3),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xyx matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT3+NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_xzx matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT3+2*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yyx_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT3+3*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_yzx_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT3+4*NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_zzx_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT3+5*NNBASX),NBAST,1,LUPRI)
          END IF

          CALL DSCAL(NNBASX,MUL3MM(1,I),WRK(KMAT3),1)              ! xxx
          CALL DSCAL(NNBASX,D3*MUL3MM(2,I),WRK(KMAT2),1)           ! xxy
          CALL DSCAL(NNBASX,D3*MUL3MM(3,I),WRK(KMAT1),1)           ! xxz
          CALL DSCAL(NNBASX,D3*MUL3MM(4,I),WRK(KMAT2+NNBASX),1)    ! xyy
          CALL DSCAL(NNBASX,D6*MUL3MM(5,I),WRK(KMAT1+NNBASX),1)    ! xyz
          CALL DSCAL(NNBASX,D3*MUL3MM(6,I),WRK(KMAT1+2*NNBASX),1)  ! xzz
          CALL DSCAL(NNBASX,MUL3MM(7,I),WRK(KMAT2+3*NNBASX),1)     ! yyy
          CALL DSCAL(NNBASX,D3*MUL3MM(8,I),WRK(KMAT1+3*NNBASX),1)  ! yyz
          CALL DSCAL(NNBASX,D3*MUL3MM(9,I),WRK(KMAT1+4*NNBASX),1)  ! yzz
          CALL DSCAL(NNBASX,MUL3MM(10,I),WRK(KMAT1+5*NNBASX),1)    ! zzz

C         Include 1/6 from the Taylor expansion
          CALL DSCAL(6*NNBASX,D6I,WRK(KMAT1),1)
          CALL DSCAL(6*NNBASX,D6I,WRK(KMAT2),1)
          CALL DSCAL(6*NNBASX,D6I,WRK(KMAT3),1)
C
C         The integrals contains a factor of -1 since they are obtained as derivative if 
C         the field gradient and these contains a factor of -1. Therefore no -1 in the daxpy.
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT3),1,WRK(KTAO),1)           ! xxx 
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT2),1,WRK(KTAO),1)           ! xxy
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT1),1,WRK(KTAO),1)           ! xxz
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT2+NNBASX),1,WRK(KTAO),1)    ! xyy
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+NNBASX),1,WRK(KTAO),1)    ! xyz
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+2*NNBASX),1,WRK(KTAO),1)  ! xzz 
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT2+3*NNBASX),1,WRK(KTAO),1)  ! yyy
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+3*NNBASX),1,WRK(KTAO),1)  ! yyz
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+4*NNBASX),1,WRK(KTAO),1)  ! yzz
          CALL DAXPY(NNBASX,FAC1,WRK(KMAT1+5*NNBASX),1,WRK(KTAO),1)  ! zzz

C         Contract with the density to get the expectation values. 
C         The factor of 1/6 in the Taylor expansion has been included.
C         Also, the off-diagonal elements have been scaled by 3 and the xyz element by 6 in order 
C         to include all contributions (the off -diagonal parts are related by symmetry)

C         Since the integrals contains a factor of -1 no -DDOT here.

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

          WRK(KNSEL+I-1) = E3XXX+E3XXY+E3XXZ+E3XYY+E3XYZ
     *                   + E3XZZ+E3YYY+E3YYZ+E3YZZ+E3ZZZ
          EMUL3T = EMUL3T + WRK(KNSEL+I-1)

C         Now the QM nuclear contribution. Note that we switch the sign here by writing
C         CORD(1,J) - MMCORD(1,I)


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
     *             +   CHARGE(J)*MUL3MM(1,I)*TXXX
     *             +   CHARGE(J)*MUL3MM(2,I)*TXXY*3
     *             +   CHARGE(J)*MUL3MM(3,I)*TXXZ*3
     *             +   CHARGE(J)*MUL3MM(4,I)*TXYY*3
     *             +   CHARGE(J)*MUL3MM(5,I)*TXYZ*6
     *             +   CHARGE(J)*MUL3MM(6,I)*TXZZ*3
     *             +   CHARGE(J)*MUL3MM(7,I)*TYYY
     *             +   CHARGE(J)*MUL3MM(8,I)*TYYZ*3
     *             +   CHARGE(J)*MUL3MM(9,I)*TYZZ*3
     *             +   CHARGE(J)*MUL3MM(10,I)*TZZZ
  411     CONTINUE

C         Remember the factor of 1/6 from the Taylor expansion
          ELOC   = D6I*ELOC

          WRK(KNSNUC+I-1) = ELOC
          ELOCT = ELOCT + WRK(KNSNUC+I-1)

  410   CONTINUE

        EOCTUP = EMUL3T + ELOCT
        ESOLT  = ESOLT  + EOCTUP

        IF ( (IPRINT.GT.5) .OR. (LOCDEB) ) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Center Octup.-electronic Octup.-nuclear Total'
         DO 412 I = 1,MMCENT
           ETEMP = WRK(KNSEL+I-1) + WRK(KNSNUC+I-1)
           WRITE(LUPRI,*) I,WRK(KNSEL+I-1),WRK(KNSNUC+I-1),ETEMP
  412    CONTINUE
  
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' Total '
         WRITE(LUPRI,*) EMUL3T, ELOCT, EMUL3T+ELOCT
         WRITE(LUPRI,*) 
        ENDIF

      END IF ! (NMULT .GE. 3)

      IF ( (IPRINT.GT.1) .OR. (LOCDEB) ) THEN
        write(lupri,*)
        write(lupri,*) 'MM-charge QM density interaction energy:',ECHART
        write(lupri,*) 'MM-dipole QM density interaction energy:',EDIPT
        write(lupri,*) 'MM-quadr. QM density interaction energy:',EQUADT
        write(lupri,*) 'MM-octup. QM density interaction energy:',EOCTUP
      ENDIF 

C     5) The polarization correction to the Fock/KS operator

      IF (IPOLTP .GT. 0) THEN

C       Zero out a list of centers having zero polarizability. We don't
C       know yet the number of zero polarizabilities so we take the
C       worst case, i.e. MXMMCT, for the length of this list

        DO 443 I=1,MXMMCT
          ZEROAL(I) = 0
 443    CONTINUE
c
c       Check if the polarizability is equal to zero; if so put -1 on
c       the list for this center. If not equal to zero put +1 on the
c       list for this center and if not touched upon leave zero

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
     &                                2*(POLMM(3,I)**2) +
     &                POLMM(4,I)**2 + 2*(POLMM(5,I)**2) +
     &                POLMM(6,I)**2
             ANORM  = SQRT(ANORM2)
             IF (ANORM .LE. THRMM) THEN
               ZEROAL(I) = -1
               LIZA = LIZA + 1
             ELSE
               ZEROAL(I) = 1
             ENDIF
           ENDIF

 400    CONTINUE

        NZERAL = LIZA - 1
        NNZAL  = MMCENT - NZERAL  ! Number of MM centers with ALPHA .NE. 0

        IF ( (IPRINT.GT.1) .OR. (LOCDEB) ) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) ' Number of polarizable sites: ', NNZAL
          WRITE(LUPRI,*)
        ENDIF

        IF (MMMAT) THEN

C         There is no need to construct ALPHA^(-1) and T matrices
C         explicitly, we can directly form the [ALPHA^(-1) - T] matrix
          KALINMT = KWRK1
          KINVAMT = KALINMT + (3*NNZAL)*(3*NNZAL)   ! [ALPHA^(-1) - T] matrix
          KUMAT   = KINVAMT + (3*NNZAL)*(3*NNZAL)   ! [ALPHA^(-1) - T]^(-1) matrix
          KINDMOM = KUMAT   + (3*NNZAL)*(3*NNZAL)   ! UMAT is just for testing
          KMAT    = KINDMOM + (3*NNZAL)             ! List for induced dipoles
          KIPVT   = KMAT    + (3*NNBASX)            ! For Rr_a integrals
          KWRKV   = KIPVT   + (3*NNZAL)             ! For matrix inv.
          KWRK2   = KWRKV   + (3*NNZAL)             ! For matrix inv.
          LWRK2   = LWRK    - KWRK2 + 1

          IF (LWRK2 .LT. 0) THEN
            CALL ERRWRK('QMMM_FCK_AO pol. part 1',-KWRK2,LWRK)
          ENDIF

          CALL DZERO(WRK(KALINMT),(3*NNZAL)*(3*NNZAL))
          CALL DZERO(WRK(KINVAMT),(3*NNZAL)*(3*NNZAL))
          CALL DZERO(WRK(KUMAT),(3*NNZAL)*(3*NNZAL))
          CALL DZERO(WRK(KINDMOM),(3*NNZAL))
          CALL DZERO(WRK(KIPVT),(3*NNZAL))
          CALL DZERO(WRK(KWRKV),(3*NNZAL))
          CALL DZERO(WRK(KMAT),(3*NNBASX))

C         FIXDIP assumes induced dipoles are calculated in a previous run.
C         Mainly due to debugging. Assumes identical molecules and order
C         of atoms in previous and current run.

          IF (.NOT. FIXDIP) THEN
            CALL GET_IND_DIPOLES_2(DCAO,NNZAL,
     &                             WRK(KALINMT),WRK(KINVAMT),WRK(KUMAT),
     &                             WRK(KINDMOM),
     &                             WRK(KWRK2),WRK(KIPVT),WRK(KWRKV),
     &                             LWRK2,IPRINT)
          ELSE
           WRITE(LUPRI,*) 'Ind. dips. from a prev. calc. read from file'
            CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WRK(KINDMOM))
          ENDIF

        ELSE IF (MMITER) THEN

          KINDMOM = KWRK1                           ! List for induced dipoles
          KMAT    = KINDMOM + (3*NNZAL)
          KWRK2   = KMAT    + (3*NNBASX)             
          LWRK2   = LWRK    - KWRK2 + 1

          IF (LWRK2 .LT. 0) THEN
            CALL ERRWRK('QMMM_FCK_AO pol. part 2',-KWRK2,LWRK)
          ENDIF

          CALL DZERO(WRK(KINDMOM),(3*NNZAL))

C         FIXDIP assumes induced dipoles are calculated in a previous run.
C         Mainly due to debugging. Assumes identical molecules and order
C         of atoms in previous and current run.

          IF (.NOT. FIXDIP) THEN
           CALL GET_IND_DIPOLES_1(DCAO,NNZAL,WRK(KINDMOM),WRK(KWRK2),
     &                             LWRK2,IPRINT)
          ELSE
           WRITE(LUPRI,*) 'Ind. dips. from a prev. calc. read from file'
            CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WRK(KINDMOM))
          ENDIF

        ENDIF

C       Compute polarization contributions to the Fock/KS matrix and
C       total solvation energy

        FACM1 = -1.0D0
        IINIM = 0 ! important should be zero due to the indexing used !

        EDELD  = 0.0D0   ! For interaction with electronic density
        EDNUC  = 0.0D0   ! For interaction with QM nuclei
        ED0MOM = 0.0D0   ! For interaction with point-charges
        ED1MOM = 0.0D0   ! For interaction with permanent dipoles
        ED2MOM = 0.0D0   ! For interaction with quadrupoles
        ED3MOM = 0.0D0   ! For interaction with octupoles
        EDMULT = 0.0D0   ! For interaction with permanent multipoles
        EPOLT  = 0.0D0   ! Total polarization energy

        DO 500 I=1,MMCENT

          IF (ZEROAL(I) .EQ. -1) GOTO 500

          DNORM2 = WRK(KINDMOM+IINIM+0)**2 + 
     *             WRK(KINDMOM+IINIM+1)**2 +
     *             WRK(KINDMOM+IINIM+2)**2 
          DNORM = SQRT(DNORM2)
          IF (ABS(DNORM) .LE. THRMM) THEN
            IINIM = IINIM+3
            GOTO 500
          ENDIF

          DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *            (MMCORD(2,I)-QMCOM(2))**2 +
     *            (MMCORD(3,I)-QMCOM(3))**2
          DIST = SQRT(DIST2)

          IF (DIST .GT. RCUTMM) THEN
              IINIM = IINIM+3
              IF (LOCDEB) THEN
                WRITE(LUPRI,*) 'Skipping induced dipole ', I
              ENDIF
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
     &                LWRK2,LABINT,INTREP,INTADR,I,TOFILE,
     &                KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
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
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
                IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
                ENDIF
  124         CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
              DIST = (DIPORG(1)-QMCOM(1))**2 +
     &               (DIPORG(2)-QMCOM(2))**2 +
     &               (DIPORG(3)-QMCOM(3))**2
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

          IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
            WRITE (LUPRI,'(/A)') ' E_x_matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_y matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

            WRITE (LUPRI,'(/A)') ' E_z matrix in QMMM_FCK_AO:'
            CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
          END IF

          CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+0),WRK(KMAT),1)
          CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+1),WRK(KMAT+NNBASX),1)
          CALL DSCAL(NNBASX,WRK(KINDMOM+IINIM+2),WRK(KMAT+2*NNBASX),1)

          CALL DAXPY(NNBASX,FACM1,WRK(KMAT),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FACM1,WRK(KMAT+NNBASX),1,WRK(KTAO),1)
          CALL DAXPY(NNBASX,FACM1,WRK(KMAT+2*NNBASX),1,WRK(KTAO),1)

C         Polarization contribution to the total energy 

C         A) Electronic contribution

          EXCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
          EYCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
          EZCOMP = DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

          ET = 0.0D0
          ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
          EDELD = EDELD + ET

C         B) Nuclear contribution

          EFNUCX = 0.0D0
          EFNUCY = 0.0D0
          EFNUCZ = 0.0D0

          DO 510 J=1,NUCIND
            CALL GET_CHARGE_ELFLD(CHARGE(J),
     &                            CORD(1,J),CORD(2,J),CORD(3,J),
     &                            MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                            ELFLDX,ELFLDY,ELFLDZ)
            EFNUCX = EFNUCX + ELFLDX
            EFNUCY = EFNUCY + ELFLDY
            EFNUCZ = EFNUCZ + ELFLDZ
 510      CONTINUE

          IF (QMDAMP) THEN
            IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
              CALL QUIT('ERROR in no. of assigned QM polarizabilities')
            ENDIF
            IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
              DIQM = 9.99D+99
              MHIT = 0
              DO 125 M=1,NUCIND
                DIQMC = (MMCORD(1,I)-CORD(1,M))**2 +
     &                  (MMCORD(2,I)-CORD(2,M))**2 +
     &                  (MMCORD(3,I)-CORD(3,M))**2
                IF (DIQMC .LE. DIQM) THEN
                  DIQM = DIQMC
                  MHIT = M
                ENDIF
  125         CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
              DIQM = (MMCORD(1,I)-QMCOM(1))**2 +
     &               (MMCORD(2,I)-QMCOM(2))**2 +
     &               (MMCORD(3,I)-QMCOM(3))**2
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

C         C) Multipole contribution

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

C         Get electric fields due to permanent moments

          DO 520 J=1,MMCENT

            IF (J .EQ. I) GOTO 520

            EXCENT = .FALSE.

            DO 522 L=1,NEXLST
              IF (EXLIST(L,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
 522        CONTINUE

            IF (.NOT. EXCENT) THEN

C             C.1  Point-charge contribution

              IF ( (NMULT .GE. 0) .AND. 
     &           (ABS(MUL0MM(J)) .GT. THRMM) ) THEN

                CALL GET_CHARGE_ELFLD(MUL0MM(J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)

                EF0MX = EF0MX + ELFLDX
                EF0MY = EF0MY + ELFLDY
                EF0MZ = EF0MZ + ELFLDZ

              ENDIF

C             C.2  Dipole contribution

              IF (NMULT .GE. 1) THEN

                CALL GET_DIPOLE_ELFLD(MUL1MM(1,J),MUL1MM(2,J),
     &                          MUL1MM(3,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)

                EF1MX = EF1MX + ELFLDX
                EF1MY = EF1MY + ELFLDY
                EF1MZ = EF1MZ + ELFLDZ

              ENDIF

C             C.3  Quadrupole contribution

              IF (NMULT .GE. 2) THEN

                CALL GET_QUADRUPOLE_ELFLD(
     &                          MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     &                          MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)

                EF2MX = EF2MX + ELFLDX
                EF2MY = EF2MY + ELFLDY
                EF2MZ = EF2MZ + ELFLDZ

              ENDIF

C             C.3  Octupole contribution

              IF (NMULT .GE. 3) THEN

                CALL GET_OCTUPOLE_ELFLD(
     &                          MUL3MM(1,J),MUL3MM(2,J),MUL3MM(3,J),
     &                          MUL3MM(4,J),MUL3MM(5,J),MUL3MM(6,J),
     &                          MUL3MM(7,J),MUL3MM(8,J),MUL3MM(9,J),
     &                          MUL3MM(10,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)

                EF3MX = EF3MX + ELFLDX
                EF3MY = EF3MY + ELFLDY
                EF3MZ = EF3MZ + ELFLDZ

              ENDIF

            ENDIF

 520      CONTINUE

C         Point-charge contribution

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

C         Quadrupole contribution

          IF (NMULT .GE. 2) THEN

            EXCOMP = WRK(KINDMOM+IINIM+0)*EF2MX
            EYCOMP = WRK(KINDMOM+IINIM+1)*EF2MY
            EZCOMP = WRK(KINDMOM+IINIM+2)*EF2MZ

            ET = 0.0D0
            ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
            ED2MOM = ED2MOM + ET

          ENDIF

C         Octupole contribution

          IF (NMULT .GE. 3) THEN

            EXCOMP = WRK(KINDMOM+IINIM+0)*EF3MX
            EYCOMP = WRK(KINDMOM+IINIM+1)*EF3MY
            EZCOMP = WRK(KINDMOM+IINIM+2)*EF3MZ

            ET = 0.0D0
            ET = ET + DMINV2*(EXCOMP + EYCOMP + EZCOMP)
            ED3MOM = ED3MOM + ET

          ENDIF


          IINIM = IINIM + 3

 500    CONTINUE

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

C        IF (MMPROP) CALL MM_PROPS(WRK(KWRK2),LWRK2,IPRINT)

 5001   FORMAT(' Polarization energy: ')
 5002   FORMAT('      Electronic contribution:   ',F15.9)
 5003   FORMAT('      Nuclear contribution:      ',F15.9)
 5004   FORMAT('      Multipole contribution:    ',F15.9)
 5005   FORMAT('      Total:                     ',F15.9)

      ENDIF   ! Polarization part

C     Finally, put back the dipole origin

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      CALL PKSYM1(WRK(KTAO),FSOL,NBAS,NSYM,+1)
      CALL QEXIT('QMMM_FCK_AO')

      RETURN
      END
C*************************************************************************
C  /* Deck qmmmfckmo */
      SUBROUTINE QMMMFCKMO(CMO,FSOL,WRK,LWRK,IPRINT)
C
C     Construct the QMMM contribution to the Fock-matrix in MO basis
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <qmmm.h>
#include <inforb.h>
#include <infopt.h>
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
C     For the openshell density we Put in zero as this is now included in 
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
C**************************************************************************************
C  /* Deck GET_IND_DIPOLES_2 */
      SUBROUTINE GET_IND_DIPOLES_2(DCAO,POLDIM,ALINMT,INVAMT,UMAT,
     &                             INDMOM,WRK,IPVT,WRKV,LWRK,IPRINT)
C
C A subroutine to calculate induced dipole moments
C
C Input:
C
C   DCAO    - density matrix in AO basis
C   POLDIM  - a number of polarizable MM centers. Actually in common as
C             NNZAL.
C
C Output:
C
C   ALINMT  - [ALPHA^(-1) - T] matrix
C   INVAMT  - [ALPHA^(-1) - T]^(-1) matrix
C   INDMOM  - a vector containing induced dipole moments
C
C From Common   
C
C   ZEROAL  - a vector containing +1 for polarizable MM centers and -1
C             for non-polarizable
C
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

      CHARACTER LLAB
      LOGICAL EXCENT,FNDLAB
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB
      INTEGER POLDIM
      DOUBLE PRECISION ALINMT,INVAMT,INDMOM
      DIMENSION ALINMT((3*POLDIM),(3*POLDIM))
      DIMENSION INVAMT((3*POLDIM),(3*POLDIM))
      DIMENSION INDMOM(3*POLDIM)
      DIMENSION IPVT(3*POLDIM),WRKV(3*POLDIM),WRK(LWRK)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CHARACTER*8 LABINT(9*MXCENT)

      DOUBLE PRECISION AMATS,AMATSI
      DIMENSION AMATS(3,3),AMATSI(3,3)
      DOUBLE PRECISION UMAT
      DIMENSION UMAT((3*POLDIM),(3*POLDIM))

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('GET_IND_DIPOLES_2')

      LOCDEB = .FALSE.

      IF (POLDIM .NE. NNZAL) THEN
        WRITE(LUPRI,*) 'ERROR in no. of polarizabilities'
        CALL QUIT('ERROR in GET_IND_DIPOLES_2')
      ENDIF

C     Allocate memory for electric field integrals and electric fields
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
      IF (LWRK1 .LT. 0) CALL ERRWRK('GET_IND_DIPOLES_2',-KEND,LWRK)

      CALL DZERO(WRK(KMAT),3*NNBASX)
      CALL DZERO(WRK(KELF),3*POLDIM)
      IF (SPLDIP) THEN
        CALL DZERO(WRK(KELFEL),3*POLDIM)
        CALL DZERO(WRK(KELFNU),3*POLDIM)
        CALL DZERO(WRK(KIMMUL),3*POLDIM)
        CALL DZERO(WRK(KIMNUC),3*POLDIM)
        CALL DZERO(WRK(KIMELD),3*POLDIM)
      ENDIF

      DO 100 I=1,3
        DO 110 J=1,3
          AMATS(I,J)  = 0.0D0
          AMATSI(I,J) = 0.0D0
 110    CONTINUE
 100  CONTINUE

      FACM1 = -1.0D0

C     Construct [ALPHA^(-1) - T] matrix
C     AND
C     form F vector due to permanent MM moments

      LRI = 1   ! Row index in the large matrix

      DO 200 I=1,MMCENT

        IF (ZEROAL(I) .EQ. -1) GOTO 200

        IF ((IPOLTP .EQ. 1) .AND. (CONMAT)) THEN
          DO 210 J=1,3
            AMATSI(J,J) = 1.0/POLIMM(I)
 210      CONTINUE
        ENDIF

        IF ((IPOLTP .EQ. 2) .AND. (CONMAT)) THEN
          AMATS(1,1) = POLMM(1,I)
          AMATS(1,2) = POLMM(2,I)
          AMATS(1,3) = POLMM(3,I)
          AMATS(2,1) = POLMM(2,I)
          AMATS(2,2) = POLMM(4,I)
          AMATS(2,3) = POLMM(5,I)
          AMATS(3,1) = POLMM(3,I)
          AMATS(3,2) = POLMM(5,I)
          AMATS(3,3) = POLMM(6,I)
          CALL DF_MATINV(AMATS,AMATSI,3)
        ENDIF

        LCI = 1  ! Column index in the large matrix

        DO 300 J=1,MMCENT

          IF ((J .EQ. I) .AND. (CONMAT)) THEN

            IF (ZEROAL(J) .EQ. -1) GOTO 300

            DO 310 K=1,3
              DO 311 L=1,3
                ALINMT(LRI+K-1,LCI+L-1) = AMATSI(K,L)
 311          CONTINUE
 310        CONTINUE

          ENDIF

          IF (J .NE. I) THEN

            DIST2 = (MMCORD(1,I)-MMCORD(1,J))**2 +
     *              (MMCORD(2,I)-MMCORD(2,J))**2 +
     *              (MMCORD(3,I)-MMCORD(3,J))**2 
            DIST = SQRT(DIST2)

            IF (DIST .GT. RCUTMM) THEN
              LCI = LCI + 3
              IF (LOCDEB) THEN
                WRITE(LUPRI,*) 'Skipping element in T^2 ', I,J
              ENDIF
              GOTO 300
            ENDIF

            EXCENT = .FALSE.

            DO 321 L=1,NEXLST
              IF (EXLIST(L,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
 321        CONTINUE

            IF (.NOT. (EXCENT)) THEN

C             Form F vector here due to permanent MM moments

C             A) point-charge contribution

              IF ( (NMULT .GE. 0) .AND.  
     &             (ABS(MUL0MM(J)) .GT. THRMM) ) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_CHARGE_ELFLD(MUL0MM(J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

C             B) dipole contribution

              IF (NMULT .GE. 1) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_DIPOLE_ELFLD(MUL1MM(1,J),MUL1MM(2,J),
     &                          MUL1MM(3,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

C             C) quadrupole contribution

              IF (NMULT .GE. 2) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_QUADRUPOLE_ELFLD(
     &                          MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     &                          MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

C             D) octupole contribution

              IF (NMULT .GE. 3) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_OCTUPOLE_ELFLD(
     &                          MUL3MM(1,J),MUL3MM(2,J),MUL3MM(3,J),
     &                          MUL3MM(4,J),MUL3MM(5,J),MUL3MM(6,J),
     &                          MUL3MM(7,J),MUL3MM(8,J),MUL3MM(9,J),
     &                          MUL3MM(10,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

            ENDIF

            IF ((ZEROAL(J) .EQ. -1) .OR. (.NOT. (CONMAT))) GOTO 300

            IF (EXCENT) THEN

              DO 330 K=1,3
                DO 331 L=1,3
                  ALINMT(LRI+K-1,LCI+L-1) = 0.0D0
 331            CONTINUE
 330          CONTINUE

            ELSE

              DIST2 = 0.0D0
              DO 340 K=1,3
                DIST2 = DIST2 + (MMCORD(K,I) - MMCORD(K,J))**2
 340          CONTINUE
              DIST  = SQRT(DIST2)
              DIST3 = DIST**3
              DIST5 = DIST**5

              DO 350 K=1,3
                DO 351 L=1,3

C               Include damping in the exponential form
C               JPC A 102 (1998) 2399

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
                  ELEM = -1.0*ELEM

                  ALINMT(LRI+K-1,LCI+L-1) = ELEM

 351            CONTINUE
 350          CONTINUE

            ENDIF

          ENDIF

          LCI = LCI + 3

 300    CONTINUE

C       Add QM region contribution to the F vector

C       A) electronic contribution

        CALL DZERO(WRK(KMAT),3*NNBASX)

        KPATOM = 0
        NOSIM  = 3
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.

        DIPORG(1) = MMCORD(1,I)
        DIPORG(2) = MMCORD(2,I)
        DIPORG(3) = MMCORD(3,I)

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 + 
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2 
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
          LRI = LRI + 3
          IF (LOCDEB) THEN
            WRITE(LUPRI,*) 'Skipping site in F vector', I
          ENDIF
          GOTO 200
        ENDIF

        RUNQM3 = .TRUE.
        CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIM,WRK(KEND),
     &              LWRK1,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)  
        RUNQM3 = .FALSE.

          IF (QMDAMP) THEN
            IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
              CALL QUIT('ERROR in no. of assigned QM polarizabilities')
            ENDIF
            IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
              DIST = 9.99D+99
              MHIT = 0
              DO 126 M=1,NUCIND
                DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                  (DIPORG(2)-CORD(2,M))**2 +
     &                  (DIPORG(3)-CORD(3,M))**2
                IF (DISTC .LE. DIST) THEN
                  DIST = DISTC
                  MHIT = M
                ENDIF
  126         CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
              DIST = (DIPORG(1)-QMCOM(1))**2 +
     &               (DIPORG(2)-QMCOM(2))**2 +
     &               (DIPORG(3)-QMCOM(3))**2
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

          IF (LOCDEB) THEN 
           write(lupri,*) 'TEMPI',TEMPI
           write(lupri,*) 'QMPOL(MHIT)',QMPOL(MHIT)
           write(lupri,*) 'MHIT',MHIT
           write(lupri,*) 'TEMP',TEMP
           write(lupri,*) 'SIJ',SIJ
           write(lupri,*) 'DFACT',DFACT
          ENDIF

        IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
          WRITE (LUPRI,'(/A)') ' E_x_matrix in GET_IND_DIPOLES_2:'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_y matrix in GET_IND_DIPOLES_2:'
          CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_z matrix in GET_IND_DIPOLES_2:'
          CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
        END IF

        EXELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
        EYELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
        EZELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

        IF (SPLDIP) THEN
          WRK(KELFEL+LRI-1+0) = EXELCO
          WRK(KELFEL+LRI-1+1) = EYELCO
          WRK(KELFEL+LRI-1+2) = EZELCO
        ELSE
          WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + EXELCO
          WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + EYELCO
          WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + EZELCO
        ENDIF

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'electronic field:',EXELCO,EYELCO,EZELCO
        ENDIF

C       B) nuclear contribution

        TJKX = 0.0D0
        TJKY = 0.0D0
        TJKZ = 0.0D0

        DO 220 J = 1,NUCIND
          CALL GET_CHARGE_ELFLD(CHARGE(J),
     &                          CORD(1,J),CORD(2,J),CORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)

          IF (QMDAMP) THEN
            IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
              CALL QUIT('ERROR in no. of assigned QM polarizabilities')
            ENDIF
            IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) )THEN
              DIQM = 9.99D+99
              MHIT = 0
              DO 127 M=1,NUCIND
                DIQMC = (MMCORD(1,I)-CORD(1,M))**2 +
     &                  (MMCORD(2,I)-CORD(2,M))**2 +
     &                  (MMCORD(3,I)-CORD(3,M))**2
                IF (DIQMC .LE. DIQM) THEN 
                  DIQM = DIQMC
                  MHIT = M
                ENDIF
  127         CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
              DIQM = (MMCORD(1,I)-QMCOM(1))**2 +
     &               (MMCORD(2,I)-QMCOM(2))**2 +
     &               (MMCORD(3,I)-QMCOM(3))**2
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

            ELFLDX = ELFLDX*DFACT
            ELFLDY = ELFLDY*DFACT
            ELFLDZ = ELFLDZ*DFACT
          ENDIF

          IF (SPLDIP) THEN
            WRK(KELFNU+LRI-1+0) = WRK(KELFNU+LRI-1+0) + ELFLDX
            WRK(KELFNU+LRI-1+1) = WRK(KELFNU+LRI-1+1) + ELFLDY
            WRK(KELFNU+LRI-1+2) = WRK(KELFNU+LRI-1+2) + ELFLDZ
          ELSE 
            WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
            WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
            WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ
          ENDIF
          TJKX = TJKX + ELFLDX
          TJKY = TJKY + ELFLDY
          TJKZ = TJKZ + ELFLDZ
 220    CONTINUE

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'nuclear field:',TJKX,TJKY,TJKZ
        ENDIF

        LRI = LRI + 3

 200  CONTINUE

      NDIM = 3*POLDIM

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'Done generating the F-Vector'
        WRITE(LUPRI,*) 'Done generating the interaction matrix'
        WRITE(LUPRI,*) 'F-Vector'
        DO 777 KK=1,NDIM
          WRITE(LUPRI,*) WRK(KELF+KK-1)
 777    CONTINUE
      ENDIF

C     If needed, invert the [ALPHA^(-1) - T] matrix and write it to 
C     file. ELSE: read matrix from the file. CONMAT = CONstruct MATrix

      IF (CONMAT) THEN

        IF (IPRINT .GT. 1) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) ' The Relay matrix is'//
     &                   ' explicitly constructed. '
          WRITE(LUPRI,*) ' Dimension is: ',NDIM
          WRITE(LUPRI,*)
        ENDIF

        IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
          WRITE(LUPRI,*)'Matrix to be inverted: '
          CALL OUTPUT(ALINMT,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
        END IF

C       Invert matrix
        CALL DGEINV(NDIM,ALINMT,INVAMT,IPVT,WRKV,INFO)

        IF (INFO .NE. 0) THEN
          CALL QUIT('ERROR in matrix inversion from QM/MMpol')
        END IF

        IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
 
          WRITE(LUPRI,*)'Inverted matrix (Relay matrix): '
          CALL OUTPUT(INVAMT,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)

          DO 999 I=1,NDIM
            DO 998 J=1,NDIM
              UMAT(I,J) = 0.0D0
              DO 997 K=1,NDIM
                UMAT(I,J) = UMAT(I,J) + ALINMT(I,K)*INVAMT(K,J)
 997          CONTINUE
 998        CONTINUE
 999      CONTINUE

          WRITE(LUPRI,*)'A * Ainv : '
          CALL OUTPUT(UMAT,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)

        END IF

C       We write the [ALPHA^(-1) - T]^(-1) matrix to file

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     $               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF

        REWIND(LUQMMM)
        CALL WRTIEF(INVAMT,NDIM*NDIM,'QQMMMMAT',LUQMMM)
        CALL GPCLOSE(LUQMMM,'KEEP')

        CONMAT = .FALSE.

      ELSE  ! read the inverted matrix from the file

        IF (IPRINT .GT. 5) THEN
          WRITE(LUPRI,*)
          WRITE(LUPRI,*) ' The Relay matrix is'//
     &                   ' read from the file. '
          WRITE(LUPRI,*)
        ENDIF

        LUQMMM = -1
        IF (LUQMMM .LT. 0) THEN
          CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUMMY,.FALSE.)
        ENDIF
        REWIND(LUQMMM)

        DO I=1,NDIM
          DO J=1,NDIM
            INVAMT(I,J) = 0.0D0
          ENDDO
        ENDDO

        IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
          CALL READT(LUQMMM,NDIM*NDIM,INVAMT)
        ELSE
          CALL QUIT('Problem reading the Relay matrix from QMMMIM file')
        ENDIF

        CALL GPCLOSE(LUQMMM,'KEEP')

        IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
          WRITE(LUPRI,*) ' The Relay matrix is'//
     &                   ' read from the QMMMIM file: '
          CALL OUTPUT(INVAMT,1,NDIM,1,NDIM,NDIM,NDIM,1,LUPRI)
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

      DO 400 I=1,NDIM
        DO 410 J=1,NDIM
          IF (SPLDIP) THEN
            WRK(KIMMUL+I-1) = WRK(KIMMUL+I-1) +
     &                        INVAMT(I,J)*WRK(KELF+J-1)
            WRK(KIMNUC+I-1) = WRK(KIMNUC+I-1) +
     &                        INVAMT(I,J)*WRK(KELFNU+J-1)
            WRK(KIMELD+I-1) = WRK(KIMELD+I-1) +
     &                        INVAMT(I,J)*WRK(KELFEL+J-1)
          ELSE
            INDMOM(I) = INDMOM(I) + INVAMT(I,J)*WRK(KELF+J-1)
          ENDIF
 410    CONTINUE
 400  CONTINUE

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

      CALL QEXIT('GET_IND_DIPOLES_2')
      RETURN
      END

C**************************************************************************************
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
C**************************************************************************************
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
C**************************************************************************************
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
      QTENS(2,3) = QXZ
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
C**************************************************************************************
C  /* Deck Put_To_File_1 */
      SUBROUTINE PUT_TO_FILE_1(FLNAME,NULOOP,DDATA)
C**************************************************************
C
#include <implicit.h>
#include <dummy.h>
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
C**************************************************************
C  /* Deck Get_From_File_1 */
      SUBROUTINE GET_FROM_FILE_1(FLNAME,NULOOP,DDATA)
C**************************************************************
C
#include <implicit.h>
#include <dummy.h>
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
C**************************************************************
C  /* Deck MM_PROPS */
      SUBROUTINE MM_PROPS(WRK,LWRK,IPRINT)
C
C  Calculates properties of the MM region.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <qm3.h>
#include <qmmm.h>
#include <infpri.h>

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

      KINVAMT    = 1
      KBMATS     = KINVAMT  + (3*NNZAL)*(3*NNZAL)
      KPOLMAT    = KBMATS   + (3*NNZAL)*3
      KPOLCORD   = KPOLMAT  + 3*3
      KEND       = KPOLCORD + 3*NNZAL
      LWRK1      = LWRK     - KEND

      IF (LWRK1 .LT. 0) CALL ERRWRK('MM_PROPS',-KEND,LWRK)

      CALL DZERO(WRK(KINVAMT),(3*NNZAL)*(3*NNZAL))
      CALL DZERO(WRK(KBMATS),(3*NNZAL)*3)
      CALL DZERO(WRK(KPOLMAT),(3*3))

      CALL MM_DIPANDCHARGE(WRK(KEND),LWRK1,IPRINT)

      IF (MMMAT) THEN
        CALL MM_POLARIZABILITY(WRK(KINVAMT),WRK(KBMATS),
     &                         WRK(KPOLMAT),IPRINT)

        CALL MM_OPTROT(WRK(KINVAMT),WRK(KPOLCORD),IPRINT)
      ELSE 
        WRITE(LUPRI,*) 'MM properties skipped since MMITER' 
      ENDIF

      WRITE(LUPRI,*) ' ---------------------------------------' 
      WRITE(LUPRI,*)

      CALL QEXIT('MM_PROPS')
      RETURN
      END

C**************************************************************************************
C  /* Deck MM_POLARIZABILITY */
      SUBROUTINE MM_POLARIZABILITY(INVAMT,BMATS,POLMAT,IPRINT)
C
C  Contracts the Relay matrix to the group and molecular
C  polarizabilities
C
#include <implicit.h>
#include <priunit.h>
#include <infpri.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>

      LOGICAL FNDLAB,LOCDEB
      DOUBLE PRECISION INVAMT,BMATS,POLMAT
      DIMENSION INVAMT(3*NNZAL,3*NNZAL)
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
        CALL READT(LUQMMM,N*N,INVAMT)
      ELSE
        CALL QUIT('Problem reading the matrix from the QMMMIM file.')
      ENDIF

      CALL GPCLOSE(LUQMMM,'KEEP')

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) )THEN
        WRITE(LUPRI,*) 'Relay mat. is read from file MM_POLARIZABILITY'
        CALL OUTPUT(INVAMT,1,N,1,N,N,N,1,LUPRI)
      ENDIF

C     Contract the Relay matrix

      K1=1
      DO 100 K = 1,NNZAL
        J1 = 1
        DO 101 J = 1,NNZAL
          BMATS(K1,1)    = BMATS(K1,1) + INVAMT(K1,J1)
          BMATS(K1,2)    = BMATS(K1,2) + INVAMT(K1,J1+1)
          BMATS(K1,3)    = BMATS(K1,3) + INVAMT(K1,J1+2)
          BMATS(K1+1,1)  = BMATS(K1+1,1) + INVAMT(K1+1,J1)
          BMATS(K1+1,2)  = BMATS(K1+1,2) + INVAMT(K1+1,J1+1)
          BMATS(K1+1,3)  = BMATS(K1+1,3) + INVAMT(K1+1,J1+2)
          BMATS(K1+2,1)  = BMATS(K1+2,1) + INVAMT(K1+2,J1)
          BMATS(K1+2,2)  = BMATS(K1+2,2) + INVAMT(K1+2,J1+1)
          BMATS(K1+2,3)  = BMATS(K1+2,3) + INVAMT(K1+2,J1+2)
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

C**************************************************************************************
C  /* Deck MM_DIPANDCHARGE */
      SUBROUTINE MM_DIPANDCHARGE(WRK,LWRK,IPRINT)
C
C     Calculates the MM total charge and dipole moment
C
#include <implicit.h>
#include <priunit.h>
#include <infpri.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>

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

      CALL DZERO(WRK(KINDMOM),(3*NNZAL))

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

C**************************************************************************************
C  /* Deck MM_OPTROT */
      SUBROUTINE MM_OPTROT(INVAMT,POLCORD,IPRINT)
C
C  Contracts the Relay matrix to the molecular optical rotation (beta)
C
#include <implicit.h>
#include <priunit.h>
#include <infpri.h>
#include <mxcent.h>
#include <qmmm.h>
#include <qm3.h>

      LOGICAL FNDLAB,LOCDEB
      DOUBLE PRECISION INVAMT,BMAT
      DIMENSION INVAMT(3*NNZAL,3*NNZAL)
      DIMENSION POLCORD(3,NNZAL),BMAT(3,3)

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
        CALL READT(LUQMMM,N*N,INVAMT)
      ELSE
        CALL QUIT('Problem reading the matrix from the QMMMIM file.')
      ENDIF

      CALL GPCLOSE(LUQMMM,'KEEP')

      IF ( (IPRRSP .GE. 15) .OR. (LOCDEB) )THEN
        WRITE(LUPRI,*) 'Relay mat. is read from file MM_OPTROT'
        CALL OUTPUT(INVAMT,1,N,1,N,N,N,1,LUPRI)
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
          BMAT(1,1) = INVAMT(K,L)
          BMAT(1,2) = INVAMT(K,L+1)
          BMAT(1,3) = INVAMT(K,L+2)
          BMAT(2,1) = INVAMT(K+1,L)
          BMAT(2,2) = INVAMT(K+1,L+1)
          BMAT(2,3) = INVAMT(K+1,L+2)
          BMAT(3,1) = INVAMT(K+2,L)
          BMAT(3,2) = INVAMT(K+2,L+1)
          BMAT(3,3) = INVAMT(K+2,L+2)
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

C**************************************************************************************
C  /* Deck GET_IND_DIPOLES_1 */
      SUBROUTINE GET_IND_DIPOLES_1(DCAO,POLDIM,INDMOM,WRK,LWRK,IPRINT)
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

      LOGICAL EXCENT,LOCDEB,DIPCON
      LOGICAL TOFILE,TRIMAT,EXP1VL
      INTEGER POLDIM
      DOUBLE PRECISION INDMOM
      DIMENSION INDMOM(3*POLDIM),WRK(LWRK)

      DOUBLE PRECISION AMAT,EVEC,MY0,TTENS,ATMAT,DIP
      DOUBLE PRECISION MY
      DIMENSION AMAT(3,3),EVEC(3),MY0(3),TTENS(3,3)
      DIMENSION ATMAT(3,3),DIP(3),MY(3)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CHARACTER*8 LABINT(9*MXCENT)

      PARAMETER ( D2 = 2.0D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('GET_IND_DIPOLES_1')

      LOCDEB = .FALSE.

      IF (POLDIM .NE. NNZAL) THEN
        WRITE(LUPRI,*) 'ERROR in no. of polarizabilities'
        CALL QUIT('ERROR in GET_IND_DIPOLES_1')
      ENDIF

      IF (SPLDIP) THEN
        WRITE(LUPRI,*) 'Split not implemented for iterative QMMM'
      ENDIF

C     Allocate memory for electric field integrals and electric fields
      KMAT    = 1                      ! For electric field integrals
      KELF    = KMAT   + 3*NNBASX      ! For total electric field     
      KEND    = KELF   + 3*POLDIM     
      LWRK1 = LWRK - KEND
      IF (LWRK1 .LT. 0) CALL ERRWRK('GET_IND_DIPOLES_1',-KEND,LWRK)

      CALL DZERO(WRK(KMAT),3*NNBASX)
      CALL DZERO(WRK(KELF),3*POLDIM)

C     1. Form F vector due to permanent MM moments

      LRI = 1   

      DO 200 I=1,MMCENT

        IF (ZEROAL(I) .EQ. -1) GOTO 200

        DO 300 J=1,MMCENT

          IF (J .NE. I) THEN

            DIST2 = (MMCORD(1,I)-MMCORD(1,J))**2 +
     *              (MMCORD(2,I)-MMCORD(2,J))**2 +
     *              (MMCORD(3,I)-MMCORD(3,J))**2 
            DIST = SQRT(DIST2)

            IF (LOCDEB) THEN
              WRITE(LUPRI,*) 'Skipping site in F-vector', I
            ENDIF
            IF (DIST .GT. RCUTMM) GOTO 300

            EXCENT = .FALSE.

            DO 321 L=1,NEXLST
              IF (EXLIST(L,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
 321        CONTINUE

            IF (.NOT. (EXCENT)) THEN

C             Form F vector here due to permanent MM moments

C             A) point-charge contribution

              IF (NMULT .GE. 0) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_CHARGE_ELFLD(MUL0MM(J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

C             B) dipole contribution

              IF (NMULT .GE. 1) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_DIPOLE_ELFLD(MUL1MM(1,J),MUL1MM(2,J),
     &                          MUL1MM(3,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

C             C) quadrupole contribution

              IF (NMULT .GE. 2) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_QUADRUPOLE_ELFLD(
     &                          MUL2MM(1,J),MUL2MM(2,J),MUL2MM(3,J),
     &                          MUL2MM(4,J),MUL2MM(5,J),MUL2MM(6,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

C             D) octupole contribution

              IF (NMULT .GE. 3) THEN

                ELFLDX = 0.0D0
                ELFLDY = 0.0D0
                ELFLDZ = 0.0D0
                CALL GET_OCTUPOLE_ELFLD(
     &                          MUL3MM(1,J),MUL3MM(2,J),MUL3MM(3,J),
     &                          MUL3MM(4,J),MUL3MM(5,J),MUL3MM(6,J),
     &                          MUL3MM(7,J),MUL3MM(8,J),MUL3MM(9,J),
     &                          MUL3MM(10,J),
     &                          MMCORD(1,J),MMCORD(2,J),MMCORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)
                WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
                WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
                WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

              ENDIF

            ENDIF
          ENDIF
 300    CONTINUE

C       Add QM region contribution to the F vector
C       A) electronic contribution

        CALL DZERO(WRK(KMAT),3*NNBASX)

        KPATOM = 0
        NOSIM  = 3
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.

        DIPORG(1) = MMCORD(1,I)
        DIPORG(2) = MMCORD(2,I)
        DIPORG(3) = MMCORD(3,I)

        DIST2 = (MMCORD(1,I)-QMCOM(1))**2 +
     *          (MMCORD(2,I)-QMCOM(2))**2 +
     *          (MMCORD(3,I)-QMCOM(3))**2 
        DIST = SQRT(DIST2)

        IF (DIST .GT. RCUTMM) THEN
          LRI = LRI + 3
          IF (LOCDEB) THEN
           WRITE(LUPRI,*) 'Skipping site in F-vector', I
          ENDIF
          GOTO 200
        ENDIF

        RUNQM3 = .TRUE.
        CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIM,WRK(KEND),
     &              LWRK1,LABINT,INTREP,INTADR,I,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)  
        RUNQM3 = .FALSE.

        IF (QMDAMP) THEN
          IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
            CALL QUIT('ERROR in no. of assigned QM polarizabilities')
          ENDIF
          IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIST = 9.99D+99
            MHIT = 0
            DO 126 M=1,NUCIND
              DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &                (DIPORG(2)-CORD(2,M))**2 +
     &                (DIPORG(3)-CORD(3,M))**2
              IF (DISTC .LE. DIST) THEN
                DIST = DISTC
                MHIT = M
              ENDIF
  126       CONTINUE
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

        IF ( (IPRINT.GT.15) .OR. (LOCDEB) ) THEN
          WRITE (LUPRI,'(/A)') ' E_x_matrix in GET_IND_DIPOLES_1:'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_y matrix in GET_IND_DIPOLES_1:'
          CALL OUTPAK(WRK(KMAT+NNBASX),NBAST,1,LUPRI)

          WRITE (LUPRI,'(/A)') ' E_z matrix in GET_IND_DIPOLES_1:'
          CALL OUTPAK(WRK(KMAT+2*NNBASX),NBAST,1,LUPRI)
        END IF

        EXELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT),1)
        EYELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT+NNBASX),1)
        EZELCO = DDOT(NNBASX,DCAO,1,WRK(KMAT+2*NNBASX),1)

        WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + EXELCO
        WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + EYELCO
        WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + EZELCO

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'electronic field:',EXELCO,EYELCO,EZELCO
        ENDIF

C       B) nuclear contribution

        TJKX = 0.0D0
        TJKY = 0.0D0
        TJKZ = 0.0D0

        DO 220 J = 1,NUCIND
          CALL GET_CHARGE_ELFLD(CHARGE(J),
     &                          CORD(1,J),CORD(2,J),CORD(3,J),
     &                          MMCORD(1,I),MMCORD(2,I),MMCORD(3,I),
     &                          ELFLDX,ELFLDY,ELFLDZ)

          IF (QMDAMP) THEN
            IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
              CALL QUIT('ERROR in no. of assigned QM polarizabilities')
            ENDIF
            IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) )THEN
              DIQM = 9.99D+99
              MHIT = 0
              DO 127 M=1,NUCIND
                DIQMC = (MMCORD(1,I)-CORD(1,M))**2 +
     &                  (MMCORD(2,I)-CORD(2,M))**2 +
     &                  (MMCORD(3,I)-CORD(3,M))**2
                IF (DIQMC .LE. DIQM) THEN 
                  DIQM = DIQMC
                  MHIT = M
                ENDIF
  127         CONTINUE
            ELSE IF (IDAMP .EQ. 2) THEN
              DIQM = (MMCORD(1,I)-QMCOM(1))**2 +
     &               (MMCORD(2,I)-QMCOM(2))**2 +
     &               (MMCORD(3,I)-QMCOM(3))**2
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

            ELFLDX = ELFLDX*DFACT
            ELFLDY = ELFLDY*DFACT
            ELFLDZ = ELFLDZ*DFACT
          ENDIF

          WRK(KELF+LRI-1+0) = WRK(KELF+LRI-1+0) + ELFLDX
          WRK(KELF+LRI-1+1) = WRK(KELF+LRI-1+1) + ELFLDY
          WRK(KELF+LRI-1+2) = WRK(KELF+LRI-1+2) + ELFLDZ

          TJKX = TJKX + ELFLDX
          TJKY = TJKY + ELFLDY
          TJKZ = TJKZ + ELFLDZ
 220    CONTINUE

        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'nuclear field:',TJKX,TJKY,TJKZ
        ENDIF

        LRI = LRI + 3

 200  CONTINUE

      NDIM = 3*POLDIM

      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'F-Vector'
        DO 899 I=1,NDIM
        WRITE(LUPRI,*) WRK(KELF+I-1)
 899    CONTINUE       
      ENDIF

C     Convert the F-vector into induced dipole moments 

      IOPT = 1 ! read file with ind. momens from previous SCF iteration.
      CALL F2QMMM(WRK(KELF),POLDIM,INDMOM,WRK(KEND),LWRK1,
     *            IOPT,IPRINT)

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

      CALL QEXIT('GET_IND_DIPOLES_1')
      RETURN
      END
C********************************************************************************************
C  /* Deck F2QMMM */
      SUBROUTINE F2QMMM(ELF,POLDIM,INDMOM,WRK,LWRK,IOPT,IPRINT)
C
C Converts a field vector into induced dipoles using iterative procedures.
C
C     Input: ELF
C     Output: INDMOM
C
C     INDMOM is the induced dipole moments
C     INDDIA is the diagonal part of the induced dipole moments, i.e. the part corresponding 
C            directly to the F ELF vector. 
C JK

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

      IF (MMDIIS) THEN
        KVEC  = KEND
        KEND2 = KVEC + (MXMMIT+1)*NDIM
        LWRK2 = LWRK - KEND2 + 1
        IF (LWRK2 .LT. 0) CALL ERRWRK('F2QMMM 2',-KEND2,LWRK)

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

      DIPCON = .FALSE.
      LM = 0
      DO 999 ITER = 1, MXMMIT
      LM = LM + 1

      LRI = 1
      DO 405 I=1,MMCENT
        IF (ZEROAL(I) .EQ. -1) GOTO 405

        LCI = 1
        DO 409 J=1,MMCENT

          IF (ZEROAL(J) .EQ. -1) GOTO 409 

            IF (J .NE. I) THEN

            EXCENT = .FALSE.
            DO 921 L=1,NEXLST
              IF (EXLIST(L,I) .EQ. EXLIST(1,J)) EXCENT = .TRUE.
 921        CONTINUE

            IF (.NOT. EXCENT) THEN
C             Get the polarizability tensor for this site
              DO 406 K=1,3
                DO 407 L=1,3
                  AMAT(K,L)  = 0.0D0
 407            CONTINUE
 406          CONTINUE

              IF (IPOLTP .EQ. 1)  THEN
                DO 408 L=1,3
                  AMAT(L,L) = POLIMM(I)
 408            CONTINUE
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

C             Now calculate the T tensor for these sites
              DIST2 = 0.0D0
              DO 410 K=1,3
                DIST2 = DIST2 + (MMCORD(K,I) - MMCORD(K,J))**2
 410          CONTINUE
              DIST  = SQRT(DIST2)
              DIST3 = DIST**3
              DIST5 = DIST**5

              DO 411 K=1,3
                DO 412 L=1,3

C               Include damping in the exponential form
C               JPC A 102 (1998) 2399
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
 412            CONTINUE
 411          CONTINUE

C             calculate alpha*T
              CALL DGEMM('N','N',3,3,3,1.D0,AMAT,3,
     &                   TTENS,3,0.D0,ATMAT,3)

              DIP(1) = WRK(KINDP+LCI-1+0)
              DIP(2) = WRK(KINDP+LCI-1+1)
              DIP(3) = WRK(KINDP+LCI-1+2)

              NLDIM = 3
              NTOTI = MAX(NLDIM,1)
              CALL DGEMV('N',NLDIM,NLDIM,D1,ATMAT,NTOTI,DIP,1,D0,MY,1)

              INDMOM(LRI+0) = INDMOM(LRI+0) + MY(1)
              INDMOM(LRI+1) = INDMOM(LRI+1) + MY(2)
              INDMOM(LRI+2) = INDMOM(LRI+2) + MY(3)

            ENDIF
          ENDIF

          LCI = LCI + 3
 409    CONTINUE

        LRI = LRI + 3
 405  CONTINUE

      TERROR=0.0D0
      DO 414 I=1,NDIM
        TERROR = TERROR + (INDMOM(I)-WRK(KINDP+I-1))*
     &                    (INDMOM(I)-WRK(KINDP+I-1)) 
 414  CONTINUE

      IF ( (LOCDEB) .OR. (IPRINT .GE. 15) ) THEN
        LMAX = 0
        TMAX = 0.0D0
        DO 413 I=1,NDIM
          TDIFF = ABS(INDMOM(I) - WRK(KINDP-1+I))
          IF (TDIFF .GT. TMAX) THEN
            TMAX = TDIFF
            LMAX = I 
          ENDIF
 413    CONTINUE
        IF (LMAX .NE. 0) THEN
          WRITE(LUPRI,*) 'Maximum deviation (element) is ',TMAX, LMAX
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
     *                               WRK(KEND2),LWRK2,IPRINT)
        ELSE 
          CALL DCOPY(NDIM,INDMOM,1,WRK(KINDP),1)
        ENDIF
C       If no convergence in last iteration keep the values for the induced dipoles, 
C       i.e. not only the diagonal part 
        IF (ITER .NE. MXMMIT) CALL DCOPY(NDIM,WRK(KINDDIA),1,INDMOM,1)
      ENDIF

  999 CONTINUE

 9000 CONTINUE !Done

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

      IF (FIRST) THEN
         FIRST = .FALSE.
      END IF

      CALL QEXIT('F2QMMM')
      RETURN
      END
C*********************************************************************************************************
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
C*********************************************************************************************************
C  /* Deck domep */
      SUBROUTINE DOMEP(CMO,WRK,KFRSAV,LFRSAV)

#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
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
#include <gnrinf.h>


      DIMENSION WRK(LFRSAV)
      DIMENSION CMO(*)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LSAVE
      LOGICAL LOCDEB
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*2 UNITS
      PARAMETER ( D1 = 1.0D0 )
      PARAMETER ( D2 = 2.0D0, DMINV2 = -0.50D0 )
      PARAMETER ( D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      CALL QENTER('DOMEP')

C     Only scf for the moment 

      KFREE  = KFRSAV
      LFREE  = LFRSAV

      LOCDEB = .TRUE.

      LUMEP = -1
      CALL GPOPEN(LUMEP,'MEP.INP','OLD',' ',
     &           'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUMEP)
      READ(LUMEP,*) NCENTERS
      READ(LUMEP,'(A2)') UNITS
      CALL UPCASE(UNITS)

      IF (UNITS .EQ. 'AA') THEN
        SCAL = XTANG
      ELSE IF (UNITS .EQ. 'AU') THEN
        SCAL = 1.0D0
      ELSE
         CALL QUIT('Unknown units in POTENTIAL.INP')
      ENDIF
   
      CALL MEMGET('REAL',KDV   ,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENC ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENV ,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDENS ,NNBASX,WRK,KFREE,LFREE)

      CALL DZERO(WRK(KDV),NNASHX)
      CALL DZERO(WRK(KDENC),N2BASX)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DZERO(WRK(KDENS),NNBASX)

      IF (NASHT .EQ. 1) THEN
         WRK(KDV) = D1
      ELSE IF (HSROHF) THEN
         DO I = 1, NASHT
            II = I*(I+1)/2
            WRK(KDV+II-1) = D1
         END DO
      ENDIF

      CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDENC),WRK(KDENV),
     &            CMO,WRK(KDV),WRK(KFREE),LFREE)
      CALL DAXPY(N2BASX,1.0D0,WRK(KDENV),1,WRK(KDENC),1)
      CALL DZERO(WRK(KDENV),N2BASX)
      CALL DGEFSP(NBAST,WRK(KDENC),WRK(KDENV))
      CALL PKSYM1(WRK(KDENV),WRK(KDENS),NBAS,NSYM,1)

      CALL MEMGET('REAL',KXMEP ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KYMEP ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZMEP ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KMEP  ,NCENTERS,WRK,KFREE,LFREE)

      CALL MEMGET('REAL',KMEFX ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KMEFY ,NCENTERS,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KMEFZ ,NCENTERS,WRK,KFREE,LFREE)

      CALL MEMGET('REAL',KMAT  ,NNBASX,WRK,KFREE,LFREE)
      LENGTH = 3*NNBASX
      CALL MEMGET('REAL',KMATEF,LENGTH,WRK,KFREE,LFREE)

      CALL DZERO(WRK(KMAT),NNBASX)
      CALL DZERO(WRK(KMATEF),3*NNBASX)

      DO 101 J = 1,NCENTERS
        READ(LUMEP,*) WRK(KXMEP+J-1),WRK(KYMEP+J-1),WRK(KZMEP+J-1)
        WRK(KXMEP+J-1) = WRK(KXMEP+J-1)/SCAL
        WRK(KYMEP+J-1) = WRK(KYMEP+J-1)/SCAL
        WRK(KZMEP+J-1) = WRK(KZMEP+J-1)/SCAL
  101 CONTINUE

      CLOSE(LUMEP)

      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      DO 102 J = 1,NCENTERS
        KPATOM = 0
        NOSIM = 1
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        IPRINT = 0
        DIPORG(1) = WRK(KXMEP+J-1)
        DIPORG(2) = WRK(KYMEP+J-1)
        DIPORG(3) = WRK(KZMEP+J-1)

        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMAT),'NPETES ',NOSIM,WRK(KFREE),
     &              LFREE,LABINT,INTREP,INTADR,J,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.

        IF (LOCDEB)  THEN
          WRITE (LUPRI,'(/A)') 'Pot. energy matrix in DOMEP'
          CALL OUTPAK(WRK(KMAT),NBAST,1,LUPRI)
        ENDIF

        WRK(KMEP+J-1) = DDOT(NNBASX,WRK(KDENS),1,WRK(KMAT),1)

        ECHCHL = 0.0D0
        DO 103 I = 1,NUCIND
          XDIS   = WRK(KXMEP+J-1) - CORD(1,I) 
          YDIS   = WRK(KYMEP+J-1) - CORD(2,I)
          ZDIS   = WRK(KZMEP+J-1) - CORD(3,I)
          DIST2  = XDIS**2+YDIS**2+ZDIS**2
          DIST   = SQRT(DIST2)
          ECHCHL = ECHCHL + CHARGE(I)/DIST
  103   CONTINUE

        WRK(KMEP+J-1) = WRK(KMEP+J-1) + ECHCHL

  102 CONTINUE

C     Calculate also electric the electric field at these points

      DO 104 J = 1,NCENTERS
        KPATOM = 0
        NOSIM = 3
        TOFILE = .FALSE.
        TRIMAT = .TRUE.
        EXP1VL = .FALSE.
        IPRINT = 0
        DIPORG(1) = WRK(KXMEP+J-1)
        DIPORG(2) = WRK(KYMEP+J-1)
        DIPORG(3) = WRK(KZMEP+J-1)

        QMMM = LSAVE
        IF (.NOT. LSAVE) QMMM = .TRUE.
        RUNQM3=.TRUE.
        CALL GET1IN(WRK(KMATEF),'NEFIELD',NOSIM,WRK(KFREE),
     &              LFREE,LABINT,INTREP,INTADR,J,TOFILE,
     &              KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRINT)
        RUNQM3=.FALSE.
        QMMM = LSAVE

        IF (LOCDEB)  THEN
           WRITE (LUPRI,'(/A)') ' Ex matrix:'
           CALL OUTPAK(WRK(KMATEF),NBAST,1,LUPRI)

           WRITE (LUPRI,'(/A)') ' Ey matrix:'
           CALL OUTPAK(WRK(KMATEF+NNBASX),NBAST,1,LUPRI)

           WRITE (LUPRI,'(/A)') ' Ez matrix:'
           CALL OUTPAK(WRK(KMATEF+2*NNBASX),NBAST,1,LUPRI)
        END IF

        WRK(KMEFX+J-1)=DDOT(NNBASX,WRK(KDENS),1,WRK(KMATEF),1)
        WRK(KMEFY+J-1)=DDOT(NNBASX,WRK(KDENS),1,WRK(KMATEF+NNBASX),1)
        WRK(KMEFZ+J-1)=DDOT(NNBASX,WRK(KDENS),1,WRK(KMATEF+2*NNBASX),1)

        ENUCX = 0.0D0
        ENUCY = 0.0D0
        ENUCZ = 0.0D0
        DO 105 I = 1,NUCIND
          XDIS   = WRK(KXMEP+J-1) - CORD(1,I)
          YDIS   = WRK(KYMEP+J-1) - CORD(2,I)
          ZDIS   = WRK(KZMEP+J-1) - CORD(3,I)
          DIST2  = XDIS**2+YDIS**2+ZDIS**2
          DIST   = SQRT(DIST2)
          DIST3  = DIST**3
          ENUCX  = ENUCX + CHARGE(I)*XDIS/DIST3
          ENUCY  = ENUCY + CHARGE(I)*YDIS/DIST3
          ENUCZ  = ENUCZ + CHARGE(I)*ZDIS/DIST3
  105   CONTINUE

        WRK(KMEFX+J-1)=WRK(KMEFX+J-1)+ENUCX
        WRK(KMEFY+J-1)=WRK(KMEFY+J-1)+ENUCY
        WRK(KMEFZ+J-1)=WRK(KMEFZ+J-1)+ENUCZ 

  104 CONTINUE

      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      LUMEP = -1
      CALL GPOPEN(LUMEP,'MEP.OUT','NEW',' ',
     &           'FORMATTED',IDUMMY,.FALSE.)
      REWIND(LUMEP)
      WRITE(LUMEP,*) NCENTERS
      WRITE(LUMEP,'(A2)') UNITS

      DO 106 J = 1,NCENTERS
        WRITE(LUMEP,1000) WRK(KXMEP+J-1)*SCAL,WRK(KYMEP+J-1)*SCAL,
     &                    WRK(KZMEP+J-1)*SCAL,WRK(KMEP+J-1),
     &                    WRK(KMEFX+J-1),WRK(KMEFY+J-1),WRK(KMEFZ+J-1)
  106 CONTINUE

      CLOSE(LUMEP) 

      CALL MEMREL('DOMEP',WRK,KFRSAV,KFRSAV,KFREE,LFREE)

 1000 FORMAT(7(F15.10,2X))

      CALL QEXIT('DOMEP')
      RETURN
      END
C*********************************************************************************************************
C  /* Deck GET_OCTUPOLE_ELFLD */
      SUBROUTINE GET_OCTUPOLE_ELFLD(QXXX,QXXY,QXXZ,QXYY,QXYZ,
     &                              QXZZ,QYYY,QYYZ,QYZZ,QZZZ,
     &                              XORI,YORI,ZORI,
     &                              XTAR,YTAR,ZTAR,
     &                              EFX,EFY,EFZ)
C
C     Calculates the electric field strength due to electric ocupoles
C     moment.
C
C     INPUT:
C
C       Qij               - the symmetry independent components of the octupole moment
C       XORI,YORI,ZORI    - position of the octupole moment
C       XTAR,YTAR,ZTAR    - position of the point where electric field is
C                           to be calculated
C       
C     OUTPUT:
C
C       EFX,EFY,EFZ       - components of the electric field strength
C                           vector
C
C JK, 2009 July
C
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
C
      DOUBLE PRECISION QXXX,QXXY,QXXZ,QXYY,QXYZ
      DOUBLE PRECISION QXZZ,QYYY,QYYZ,QYZZ,QZZZ
      DOUBLE PRECISION XORI,YORI,ZORI
      DOUBLE PRECISION XTAR,YTAR,ZTAR
      DOUBLE PRECISION EFX,EFY,EFZ

      DOUBLE PRECISION QTENS,ELFVEC,CORDO,CORDT
      DIMENSION QTENS(3,3,3),ELFVEC(3),CORDO(3),CORDT(3)

      CALL QENTER('GET_OCTUPOLE_ELFLD')

      EFX = 0.0D0
      EFY = 0.0D0
      EFZ = 0.0D0

      DIST2  = 0.0D0
      DIST2  = DIST2 + (XTAR - XORI)**2
      DIST2  = DIST2 + (YTAR - YORI)**2
      DIST2  = DIST2 + (ZTAR - ZORI)**2
      DIST   = SQRT(DIST2)
      DIST9  = DIST**9
      DISTI9 = 1.0D0/DIST9
      DIST7  = DIST**7
      DISTI7 = 1.0D0/DIST7
      DIST5  = DIST**5
      DISTI5 = 1.0D0/DIST5

      QTENS(1,1,1) = QXXX
      QTENS(1,1,2) = QXXY
      QTENS(1,1,3) = QXXZ
      QTENS(1,2,1) = QXXY
      QTENS(1,2,2) = QXYY
      QTENS(1,2,3) = QXYZ
      QTENS(1,3,1) = QXXZ
      QTENS(1,3,2) = QXYZ
      QTENS(1,3,3) = QXZZ
      QTENS(2,1,1) = QXXY
      QTENS(2,1,2) = QXYY
      QTENS(2,1,3) = QXYZ
      QTENS(2,2,1) = QXYY
      QTENS(2,2,2) = QYYY
      QTENS(2,2,3) = QYYZ
      QTENS(2,3,1) = QXYZ
      QTENS(2,3,2) = QYYZ
      QTENS(2,3,3) = QYZZ
      QTENS(3,1,1) = QXXZ
      QTENS(3,1,2) = QXYZ
      QTENS(3,1,3) = QXZZ
      QTENS(3,2,1) = QXYZ
      QTENS(3,2,2) = QYYZ
      QTENS(3,2,3) = QYZZ
      QTENS(3,3,1) = QXZZ
      QTENS(3,3,2) = QYZZ
      QTENS(3,3,3) = QZZZ

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
            DO 130 L=1,3

              ELEM = 0.0D0
              ELEM = (105*(CORDT(L) - CORDO(L))*
     &                    (CORDT(K) - CORDO(K))*
     &                    (CORDT(J) - CORDO(J))*
     &                    (CORDT(I) - CORDO(I)))*
     &                     DISTI9
              IF (L .EQ. K) THEN
                ELEM = ELEM - 15.0D0*DISTI7*(CORDT(I) - CORDO(I))*
     *                                      (CORDT(J) - CORDO(J))
              ENDIF
              IF (L .EQ. J) THEN
                ELEM = ELEM - 15.0D0*DISTI7*(CORDT(I) - CORDO(I))*
     *                                      (CORDT(K) - CORDO(K))
              ENDIF
              IF (K .EQ. J) THEN
                ELEM = ELEM - 15.0D0*DISTI7*(CORDT(I) - CORDO(I))*
     *                                      (CORDT(L) - CORDO(L))
              ENDIF
              IF (L .EQ. I) THEN
                ELEM = ELEM - 15.0D0*DISTI7*(CORDT(J) - CORDO(J))*
     *                                      (CORDT(K) - CORDO(K))
              ENDIF
              IF (K .EQ. I) THEN
                ELEM = ELEM - 15.0D0*DISTI7*(CORDT(J) - CORDO(J))*
     *                                      (CORDT(L) - CORDO(L))
              ENDIF
              IF (J .EQ. I) THEN
                ELEM = ELEM - 15.0D0*DISTI7*(CORDT(K) - CORDO(K))*
     *                                      (CORDT(L) - CORDO(L))
              ENDIF
         
              IF ((I .EQ. J) .AND. (K .EQ. L)) THEN
                ELEM = ELEM + 3.0D0*DISTI5
              ENDIF
              IF ((I .EQ. K) .AND. (J .EQ. L)) THEN
                ELEM = ELEM + 3.0D0*DISTI5
              ENDIF
              IF ((I .EQ. L) .AND. (J .EQ. K)) THEN
                ELEM = ELEM + 3.0D0*DISTI5
              ENDIF
              ELEM = ELEM*QTENS(L,K,J)

              ELFVEC(I) = ELFVEC(I) + ELEM

 130        CONTINUE
 120      CONTINUE
 110    CONTINUE
        ELFVEC(I) = ELFVEC(I)/6.0
 100  CONTINUE

      EFX = ELFVEC(1)
      EFY = ELFVEC(2)
      EFZ = ELFVEC(3)

      CALL QEXIT('GET_OCTUPOLE_ELFLD')

      RETURN
      END
C**************************************************************************************

