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
C-----------------------------------------------------------------------
C  /* Deck qm3lno */
      SUBROUTINE QM3LNO(NOSIM,BOVECS,CREF,CMO,XINDX,UDV,DV,UDVTR,DVTR,
     &                  EVECS,WRK,LWRK)
C
C  CBN+JK, Nov. 05
C
C  Purpose:  Calculate SCF/DFT E^2 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector.
C  
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>


      DIMENSION BOVECS(*), CMO(*), XINDX(*), UDV(*), EVECS(KZYVAR,*)
      DIMENSION UDVTR(*), DVTR(*)
      DIMENSION WRK(*), DV(*), CREF(*)
      LOGICAL LOPEN 

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      LOPEN = .FALSE.

      CALL QENTER('QM3LNO')
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3LNO')
         RETURN
      END IF

      KUCMO   = 1
      KUBO    = KUCMO   + NORBT*NBAST
      KFSOLMO = KUBO    + NOSIM*N2ORBX
      KUFSOL  = KFSOLMO + NNORBX
      KUTGX   = KUFSOL  + N2ORBX
      KRXYO   = KUTGX   + N2ORBX
      IF (TRPLET) THEN
        KRXYOT  = KRXYO   + NOSIM*N2ORBX
        KRRAOx  = KRXYOT  + NOSIM*N2ORBX
      ELSE
        KRRAOx  = KRXYO   + NOSIM*N2ORBX
      ENDIF
      KURXAC  = KRRAOx  + NNBASX
      KTRMO   = KURXAC  + N2ASHX
      KUTR    = KTRMO   + NNORBX
      KUTRX   = KUTR    + N2ORBX
      KWRK    = KUTRX   + N2ORBX
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3LNO',-KWRK,LWRK)

      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KUTGX),N2ORBX)
      CALL DZERO(WRK(KRXYO),NOSIM*N2ORBX)
      IF (TRPLET) CALL DZERO(WRK(KRXYOT),NOSIM*N2ORBX)
      CALL DZERO(WRK(KRRAOx),NNBASX)
      CALL DZERO(WRK(KURXAC),N2ASHX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)

      IF (IPRRSP .GE. 140) THEN
          WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM QM3LNO >>>'
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on entry'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF


      CALL UPKCMO(CMO,WRK(KUCMO))

      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WRK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WRK(KUBO),1)
      END IF

C--------------------------------------------------------
C 1. Construct one-index transformed T^g. 
C    First read qm/mm mo fock matrix from interface file.
C    This includes N_s contribution if OLDTG = .TRUE.. 
C    If OLDTG = .FALSE. N_s is included in one-electron 
C    part of vacuum Hamiltonian
C--------------------------------------------------------

       IF (INTDIR) THEN
         OBKPX = DIPORG(1)
         OBKPY = DIPORG(2)
         OBKPZ = DIPORG(3)
       ENDIF

      IF (LUSIFC .LE. 0) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         LOPEN = .TRUE.
      END IF
     
      REWIND(LUSIFC)

      MMORBT = MAX(4,NNORBT)
      IF (MMORBT .NE. NNORBX) THEN
        CALL QUIT('ERROR IN DIMENSION OF QM/MM FOCK CONTR. IN QM3LNO')
      END IF

      CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
      CALL READT (LUSIFC,NNORBX,WRK(KFSOLMO))

      IF (LOPEN) THEN
         CALL GPCLOSE(LUSIFC,'KEEP')
      END IF

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

      DO 200 IOSIM = 1,NOSIM

         JRXYO = KRXYO + (IOSIM-1)*N2ORBX
         IF (TRPLET) JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
         JUBO  = KUBO  + (IOSIM-1)*N2ORBX

         CALL DZERO(WRK(KUTGX),N2ORBX)

         CALL ONEXH1(WRK(JUBO),WRK(KUFSOL),WRK(KUTGX))

         FAC = 1.0D0
         IF (TRPLET) THEN
           CALL DAXPY(N2ORBX,FAC,WRK(KUTGX),1,WRK(JRXYOT),1)
         ELSE
           CALL DAXPY(N2ORBX,FAC,WRK(KUTGX),1,WRK(JRXYO),1)
         ENDIF

C        if hf/dft closed shell no contribution can arrise from
C        the one-index transformed electric field if triplet.
C        Also, if LGSPOL neglect the time-dependece of the expectation value
C        contained in the interaction operator, i.e. keep the ground state polarization only.
C
         IF ( ((NASHT.EQ.0) .AND. (TRPLET)) .OR. (LGSPOL) ) GOTO 300

C-------------------------------------------------------
C 2. Construct -alpha \sum_a <0|Q_a|0>Rra
C-------------------------------------------------------

         IF (.NOT. INTDIR) THEN
            LUCCEF=-1
            CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &                  'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND(LUCCEF)
         ENDIF
C
C        Readin Rra integrals (ao) and transform to mo
C        ---------------------------------------------
C
         LM = 0

         IF (INTDIR) L = NUSITE + NUALIS(0)

         DO 520 I = 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
               DO 521 J = NSYSBG(I), NSYSED(I)
                  DO 522 K = 1, NUALIS(I)
                     LM = LM + 1

                     FAC = -ALPIMM(I,K)
C x-component
                     CALL DZERO(WRK(KRRAOx),NNBASX)
                     CALL DZERO(WRK(KTRMO),NNORBX)
                     CALL DZERO(WRK(KUTR),N2ORBX)
                     CALL DZERO(WRK(KUTRX),N2ORBX)
                     CALL DZERO(WRK(KURXAC),N2ASHX)

                     IF (INTDIR) THEN
                        KMAT = KWRK
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
                        CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIMI,
     &                         WRK(KLAST),LWRK2,LABINT,INTREP,INTADR,
     &                         IATNOW,TOFILE,KPATOM,TRIMAT,DUMMY,EXP1VL,
     &                         DUMMY,IPRRSP)
                        RUNQM3 = .FALSE.

                        IF (QMDAMP) THEN
                          DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                           (DIPORG(2)-QMCOM(2))**2 +
     &                           (DIPORG(3)-QMCOM(3))**2
                          DIST = SQRT(DIST)
                          DFACT = (1-exp(-ADAMP*DIST))**3
                          CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
                        ENDIF

                        CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),
     &                           WRK(KLAST),NBAST,NORBT)
                        CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                        CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                        IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),
     &                                                WRK(KURXAC))
                        IF (TRPLET) THEN
                          FACx = FAC*SLVTLM(UDVTR,WRK(KURXAC),
     &                                      WRK(KUTRX),TAC)
                          CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,
     &                               WRK(JRXYOT),1)
                        ELSE
                          FACx = FAC*SLVQLM(UDV,WRK(KURXAC),
     &                                      WRK(KUTRX),TAC)
                          CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,
     &                               WRK(JRXYO),1)
                        ENDIF
                     ELSE
                        CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                        CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                        CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                        CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                        IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),
     &                                                WRK(KURXAC))
                        IF (TRPLET) THEN
                         FACx = FAC*SLVTLM(UDVTR,WRK(KURXAC),
     &                                     WRK(KUTRX),TAC)
                         CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,
     &                              WRK(JRXYOT),1)
                        ELSE
                         FACx = FAC*SLVQLM(UDV,WRK(KURXAC),
     &                                     WRK(KUTRX),TAC)
                         CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,
     &                              WRK(JRXYO),1)
                        ENDIF
                     ENDIF
C y-component
                     CALL DZERO(WRK(KRRAOx),NNBASX)
                     CALL DZERO(WRK(KTRMO),NNORBX)
                     CALL DZERO(WRK(KUTR),N2ORBX)
                     CALL DZERO(WRK(KUTRX),N2ORBX)
                     CALL DZERO(WRK(KURXAC),N2ASHX)

                     IF (INTDIR) THEN
                        CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),
     &                       WRK(KUCMO),WRK(KLAST),NBAST,NORBT)
                     ELSE
                        CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                        CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                     ENDIF
                     CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                     CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                     IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),
     &                                             WRK(KURXAC))
                     IF (TRPLET) THEN
                      FACy = FAC*SLVTLM(UDVTR,WRK(KURXAC),
     &                                  WRK(KUTRX),TAC)
                      CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYOT),1)
                     ELSE
                      FACy = FAC*SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
                      CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYO),1)
                     ENDIF

C z-component
                     CALL DZERO(WRK(KRRAOx),NNBASX)
                     CALL DZERO(WRK(KTRMO),NNORBX)
                     CALL DZERO(WRK(KUTR),N2ORBX)
                     CALL DZERO(WRK(KUTRX),N2ORBX)
                     CALL DZERO(WRK(KURXAC),N2ASHX)
                     
                     IF (INTDIR) THEN
                        CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),
     &                       WRK(KUCMO),WRK(KLAST),NBAST,NORBT)
                     ELSE
                        CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                        CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),
     &                       WRK(KWRK),NBAST,NORBT)
                     ENDIF
                     CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
                     CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
                     IF (NASHT .GT. 0) CALL GETACQ(WRK(KUTRX),
     &                                             WRK(KURXAC))
                     IF (TRPLET) THEN
                       FACz = FAC*SLVTLM(UDVTR,WRK(KURXAC),
     &                                   WRK(KUTRX),TAC)
                       CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYOT),1)
                     ELSE
                       FACz = FAC*SLVQLM(UDV,WRK(KURXAC),WRK(KUTRX),TAC)
                      CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYO),1)
                     ENDIF

 522              CONTINUE
 521           CONTINUE
            END IF
 520     CONTINUE
C
         IF (.NOT. INTDIR) CALL GPCLOSE(LUCCEF,'KEEP')
 300     CONTINUE 
 200  CONTINUE

      IF (INTDIR) THEN
         DIPORG(1) = OBKPX
         DIPORG(2) = OBKPY
         DIPORG(3) = OBKPZ
      ENDIF

C--------------------------------------------
C 3.  Add effective operators to response
C--------------------------------------------


      IF (TRPLET) THEN
         CALL SLVSOR(.TRUE.,.FALSE.,NOSIM,UDVTR,EVECS(1,1),WRK(KRXYO))
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYOT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYO))
      ENDIF

      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on exit'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF

      CALL QEXIT('QM3LNO')
      RETURN
      END 
C
C-----------------------------------------------------------------------
      SUBROUTINE QM3LNO_P(NOSIM,BOVECS,CREF,CMO,XINDX,UDV,DV,UDVTR,
     &               EVECS,WRK,LWRK)
C                        
C  Arnfinn, Aug. -08: The parallel version of QM3LNO
C
C  Purpose:  Calculate SCF/DFT E^2 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector.
C
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
C#include <infpar.h>

      DIMENSION BOVECS(*), CMO(*), XINDX(*), UDV(*), EVECS(KZYVAR,*)
      DIMENSION UDVTR(*)
      DIMENSION WRK(*), DV(*), CREF(*)
      LOGICAL LOPEN

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      IF (NASHT .GT. 0) THEN
         CALL QUIT('QM3 only implemented for clsoed shell')
      ENDIF

      LOPEN = .FALSE.

      CALL QENTER('QM3LNO_P')

      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3LNO_P')
         RETURN
      END IF

      KUCMO   = 1
      KUBO    = KUCMO   + NORBT*NBAST
      KFSOLMO = KUBO    + NOSIM*N2ORBX
      KUFSOL  = KFSOLMO + NNORBX
      KUTGX   = KUFSOL  + N2ORBX
      KRXYO   = KUTGX   + N2ORBX
      IF (TRPLET) THEN
         KRXYOT  = KRXYO   + NOSIM*N2ORBX
         KRRAOx  = KRXYOT  + NOSIM*N2ORBX
      ELSE
         KRRAOx  = KRXYO   + NOSIM*N2ORBX
      ENDIF
      KTRMO   = KRRAOx  + NNBASX
      KUTR    = KTRMO   + NNORBX
      KUTRX   = KUTR    + N2ORBX
      KWRK    = KUTRX   + N2ORBX
      LWRK1   = LWRK    - KWRK

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3LNO_P',-KWRK,LWRK)

      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KUBO),NOSIM*N2ORBX)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KUTGX),N2ORBX)
      CALL DZERO(WRK(KRXYO),NOSIM*N2ORBX)
      IF (TRPLET) CALL DZERO(WRK(KRXYOT),NOSIM*N2ORBX)
      CALL DZERO(WRK(KRRAOx),NNBASX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)
      CALL DZERO(WRK(KUTRX),N2ORBX)

      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM QM3LNO >>>'
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on entry'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF


      CALL UPKCMO(CMO,WRK(KUCMO))

      IF (NOSIM.GT.0) THEN
         CALL RSPZYM(NOSIM,BOVECS,WRK(KUBO))
         CALL DSCAL(NOSIM*N2ORBX,-1.0D0,WRK(KUBO),1)
      END IF
C--------------------------------------------------------
C 1. Construct one-index transformed T^g.
C    First read qm/mm mo fock matrix from interface file.
C    This includes N_s contribution if OLDTG = .TRUE..
C    If OLDTG = .FALSE. N_s is included in one-electron
C    part of vacuum Hamiltonian
C--------------------------------------------------------

      IF (INTDIR) THEN
         OBKPX = DIPORG(1)
         OBKPY = DIPORG(2)
         OBKPZ = DIPORG(3)
      ENDIF

      IF (LUSIFC .LE. 0) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         LOPEN = .TRUE.
      END IF

      REWIND(LUSIFC)

      MMORBT = MAX(4,NNORBT)
      IF (MMORBT .NE. NNORBX) THEN
         CALL QUIT('ERROR IN DIMENSION OF QM/MM FOCK CONTR. IN QM3LNO')
      END IF

      CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
      CALL READT (LUSIFC,NNORBX,WRK(KFSOLMO))

      IF (LOPEN) THEN
         CALL GPCLOSE(LUSIFC,'KEEP')
      END IF

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

      DO 200 IOSIM = 1,NOSIM

         JRXYO = KRXYO + (IOSIM-1)*N2ORBX
         IF (TRPLET) JRXYOT = KRXYOT + (IOSIM-1)*N2ORBX
         JUBO  = KUBO  + (IOSIM-1)*N2ORBX

         CALL DZERO(WRK(KUTGX),N2ORBX)
         CALL ONEXH1(WRK(JUBO),WRK(KUFSOL),WRK(KUTGX))

         FAC = 1.0D0

         IF (TRPLET) THEN
            CALL DAXPY(N2ORBX,FAC,WRK(KUTGX),1,WRK(JRXYOT),1)
         ELSE
            CALL DAXPY(N2ORBX,FAC,WRK(KUTGX),1,WRK(JRXYO),1)
         ENDIF

C        Since only hf/dft closed shell no contribution can arrise from
C        the one-index transformed electric field if triplet.
C        Also, if LGSPOL neglect the time-dependece of the expectation value
C        contained in the interaction operator, i.e. keep the ground state 
C        polarization only.
C

         IF (TRPLET. OR. LGSPOL) GOTO 300

C-------------------------------------------------------
C 2. Construct -alpha \sum_a <0|Q_a|0>Rra
C-------------------------------------------------------

         IF (.NOT. INTDIR) THEN
            LUCCEF=-1
            CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
            REWIND(LUCCEF)
         ENDIF
C
C        Readin Rra integrals (ao) and transform to mo
C        ---------------------------------------------
C
         LM = 0

         IF (INTDIR) L = NUSITE + NUALIS(0)

         DO 520 I = 1, ISYTP
            IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
               CALL QM3LNO_M(WRK(KRRAOx),WRK(KTRMO),WRK(KUTR),
     &                      WRK(KUTRX),IATNOW,NOSIM,WRK(KWRK),LWRK1,I,
     &                      LM,L,WRK(JUBO),WRK(JRXYO),WRK(KUCMO),
     &                      .FALSE.,IPRINT)
            END IF
 520     CONTINUE
         IF (.NOT. INTDIR) CALL GPCLOSE(LUCCEF,'KEEP')
 300     CONTINUE
 200  CONTINUE

      IF (INTDIR) THEN
         DIPORG(1) = OBKPX
         DIPORG(2) = OBKPY
         DIPORG(3) = OBKPZ
      ENDIF

C--------------------------------------------
C 3.  Add effective operators to response
C--------------------------------------------

      IF (TRPLET) THEN
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYOT))
      ELSE
         CALL SLVSOR(.TRUE.,.TRUE.,NOSIM,UDV,EVECS(1,1),WRK(KRXYO))
      ENDIF

      IF (IPRRSP .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> QM3LNO - svec(orb,1) on exit'
         WRITE (LUPRI,'(/A)') ' Z - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,1,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' Y - PART OR VECTOR '
         CALL OUTPUT(EVECS(KZVAR+KZCONF+1,1),1,KZWOPT,1,1,KZWOPT,
     *               1,1,LUPRI)
      END IF

      CALL QEXIT('QM3LNO_P')
      RETURN
      END
C 
C
C---------------------------------------------------------------------
      SUBROUTINE QM3LNO_M(RRAOx,TRMO,UTR,UTRX,IATNOW,NOSIM,WRK,LWRK,
     &                   INUM,LM,LNUM,UBO,RXYO,UCMO,TOFILE,IPRINT)
C
C  Arnfinn, Aug.-08: The master routine for QM3LNO_P
C
#include <implicit.h>
#include <maxorb.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif
#include <cbiher.h>

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      INTEGER IATNOW, INUM, LM, LNUM, NOSIM
      DIMENSION WRK(LWRK)
      DIMENSION RRAOx(NNBASX), TRMO(NNORBX), UTR(N2ORBX), UTRX(N2ORBX)
      DIMENSION UBO(N2ORBX), RXYO(N2ORBX), UCMO(NORBT*NBAST)

      CALL QENTER('QM3LNO_M')

      IF (TOFILE) CALL QUIT('Parallel calculations do not allow '//
     &                     'for storing integrals on disk')

      IPRTYP = 17
      KRXYO = 1
      JRXYO = KRXYO + N2ORBX
      KLAST = JRXYO + N2ORBX

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)
C
      CALL MPIXBCAST(MXORB,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LNUM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(UBO,N2ORBX,'DOUBLE',MASTER)
      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)


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

      CALL DZERO(WRK(KRXYO),N2ORBX)
      CALL MPI_REDUCE(WRK(KRXYO),WRK(JRXYO),N2ORBX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)
      FAC = 1.0D0
      CALL DAXPY(N2ORBX,FAC,WRK(JRXYO),1,RXYO,1)

      CALL QEXIT('QM3LNO_M')

      RETURN
      END
C---------------------------------------------------------------
      SUBROUTINE QM3LNO_S(WRK,LWRK,IPRTMP)
C
C  Arnfinn, Aug.-08: The slave routine for QM3LNO_P
C
#include <implicit.h>
#include <maxorb.h>
#include <gnrinf.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <maxash.h>
#include <infind.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <dummy.h>
#include <inftap.h>
#include <ccinftap.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <infinp.h>
#include <nuclei.h>
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif

#include <cbiher.h>

      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)
      INTEGER IATNOW, NOSIM

      CALL QENTER('QM3LNO_S')

      IPRRSP = IPRTMP
      QM3 = .TRUE.


      CALL MPIXBCAST(MXORB,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(L,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)


      KRRAOx = 1
      KTRMO  = KRRAOx + NNBASX
      KUTR   = KTRMO  + NNORBX
      KUTRX  = KUTR   + N2ORBX
      JUBO   = KUTRX  + N2ORBX
      JRXYO  = JUBO   + N2ORBX
      KUCMO  = JRXYO  + N2ORBX
      KLAST1 = KUCMO  + NORBT*NBAST

      CALL MPIXBCAST(WRK(JUBO),N2ORBX,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)

      CALL DZERO(WRK(JRXYO),N2ORBX)

 187  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(LM, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(FAC, 1, 'DOUBLE', MASTER, MPTAG2)

      IF (LM .GT. 0) THEN

C x-component

         CALL DZERO(WRK(KRRAOx),NNBASX)
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
         CALL DZERO(WRK(KUTRX),N2ORBX)

         KMAT = KLAST1
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
         CALL GET1IN(WRK(KMAT),'NEFIELD',NOSIMI,
     &              WRK(KLAST),LWRK2,LABINT,INTREP,INTADR,
     &              IATNOW,TOFILE,KPATOM,TRIMAT,DUMMY,EXP1VL,
     &              DUMMY,IPRRSP)
         RUNQM3 = .FALSE.
         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),
     &            WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
         FACx = FAC * QM3QLM(WRK(KUTRX))
         CALL DAXPY(N2ORBX,FACx,WRK(KUTR),1,WRK(JRXYO),1)
C y-component
         CALL DZERO(WRK(KRRAOx),NNBASX)
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
         CALL DZERO(WRK(KUTRX),N2ORBX)

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),
     &            WRK(KUCMO),WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
         FACy = FAC * QM3QLM(WRK(KUTRX))
         CALL DAXPY(N2ORBX,FACy,WRK(KUTR),1,WRK(JRXYO),1)
C z-component
         CALL DZERO(WRK(KRRAOx),NNBASX)
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)
         CALL DZERO(WRK(KUTRX),N2ORBX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),
     &            WRK(KUCMO),WRK(KLAST),NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
         CALL ONEXH1(WRK(JUBO),WRK(KUTR),WRK(KUTRX))
         FACz = FAC * QM3QLM(WRK(KUTRX))
         CALL DAXPY(N2ORBX,FACz,WRK(KUTR),1,WRK(JRXYO),1)
         GO TO 187
      END IF

C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(JRXYO),WRK(JRXYO),N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QM3LNO_S')

      RETURN
      END

C ---------------------------------------------------------------------

C  /* Deck qm3qro */
      SUBROUTINE QM3QRO(VEC1,VEC2,ETRS,XINDX,ZYM1,ZYM2,
     &                 UDV,WRK,LWRK,KZYVR,KZYV1,KZYV2,
     &                 IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP)
C
C  CBN+JK, Dec. 05
C
C  Purpose:  Calculate SCF/DFT E^3 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector.
C
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include <rspprp.h>
#include <infcr.h>
#include <inftap.h>
#include <qrinf.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <ccinftap.h>
#include "nuclei.h"

      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(LWRK),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2)
      DIMENSION MJWOP(2,MAXWOP,8)
      LOGICAL LCON, LORB, LREF, LOPEN
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QM3QRO')
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3QRO')
         RETURN
      END IF

      LOPEN = .FALSE.

      IF (NASHT.GE.1) CALL QUIT('ONLY CLOSED SHELL FOR QR-QM/MM')

      KCREF   = 1
      KTRES   = KCREF   + NCREF
      KUCMO   = KTRES   + N2ORBX
      KFSOLMO = KUCMO   + NORBT*NBAST
      KUFSOL  = KFSOLMO + NNORBX
      KTLMA   = KUFSOL  + N2ORBX
      KTLMB   = KTLMA   + N2ORBX
      KRRAOx  = KTLMB   + N2ORBX
      KTRMO   = KRRAOx  + NNBASX
      KUTR    = KTRMO   + NNORBX
      KWRK    = KUTR    + N2ORBX
      LWRK1   = LWRK    - KWRK 

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3QRO',-KWRK,LWRK1)

      CALL DZERO(WRK(KCREF),NCREF)
      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL DZERO(WRK(KTLMB),N2ORBX)
      CALL DZERO(WRK(KRRAOx),NNBASX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      NSIM = 1

C We assume no symmetry in the DFT/MM calculations although we have kept 
C some symmetry options below ...

      ISYMT = 1

C Get the reference state

      CALL GETREF(WRK(KCREF),MZCONF(1))


C Unpack the response vectors

      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)


C Unpack symmetry blocked CMO

      CALL UPKCMO(CMO,WRK(KUCMO))

C--------------------------------------------------------
C 1. Construct two-index transformed T^g.
C    First read qm/mm mo fock matrix from interface file.
C    This includes N_s contribution if OLDTG = .TRUE..
C    If OLDTG = .FALSE. N_s is included in one-electron
C    part of vacuum Hamiltonian
C--------------------------------------------------------

      IF (LUSIFC .LE. 0) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         LOPEN = .TRUE.
      END IF

      MMORBT = MAX(4,NNORBT)
      IF (MMORBT .NE. NNORBX) THEN
        CALL QUIT('ERROR IN DIMENSION OF QM/MM FOCK CONTR. IN QM3QRO')
      END IF
C
      REWIND(LUSIFC)

      CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
      CALL READT (LUSIFC,NNORBX,WRK(KFSOLMO))

      IF (LOPEN) THEN
         CALL GPCLOSE(LUSIFC,'KEEP')
      END IF

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

C     Create the effective operator:
C
C     TRES = Tg(k1,k2) + A1(k2) + A12

      IF (ISYMT.EQ.ISYMV1 .AND. ISYMT.EQ.ISYMV2
     &    .AND. ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
          CALL DZERO(WRK(KTLMA),N2ORBX)
          CALL DZERO(WRK(KTLMB),N2ORBX)
          CALL OITH1(ISYMV1,ZYM1,WRK(KUFSOL),WRK(KTLMA),ISYMT)
          CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV1))
          CALL DAXPY(N2ORBX,0.5D0,WRK(KTLMB),1,WRK(KTRES),1)
      END IF

C     If LGSPOL neglect the time-dependece of the expectation value 
C     contained in the interaction operator, i.e. keep the ground state polarization only. 

      IF (LGSPOL) GOTO 600


      IF (.NOT. INTDIR) THEN
        CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND(LUCCEF)
      ENDIF

C     Readin Rra integrals (ao) and transform to mo
C     ---------------------------------------------

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

C             A1(k2) = F1*Rra(k2)
C             A12 = F2*Rra
C
C             F1 = 2*<0| Rra(k1) |0>
C             F2 = <0| Rra(k1,k2) |0>

C x-component
              F1=D0
              F2=D0
              CALL DZERO(WRK(KRRAOx),NNBASX)
              CALL DZERO(WRK(KTRMO),NNORBX)
              CALL DZERO(WRK(KUTR),N2ORBX)

              IF (INTDIR) THEN
                KMAT = KWRK
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
     &                      KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
                RUNQM3 = .FALSE.

                IF (QMDAMP) THEN
                  DIST = (DIPORG(1)-QMCOM(1))**2 +
     &                   (DIPORG(2)-QMCOM(2))**2 +
     &                   (DIPORG(3)-QMCOM(3))**2
                  DIST = SQRT(DIST)
                  DFACT = (1-exp(-ADAMP*DIST))**3
                  CALL DSCAL(3*NNBASX,DFACT,WRK(KMAT),1)
                ENDIF

                CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &                    NBAST,NORBT)
              ELSE
                CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                    NBAST,NORBT)
              ENDIF 

              CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

              IF (ISYMT.EQ.ISYMV1) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'QM3QRO')
                 F1 = 1.0D0*F1*FAC

                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)
              END IF

              IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL DZERO(WRK(KTLMB),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
                 CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'C3SOL')
                 F2 = 0.5D0*F2*FAC

                 CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
              END IF

C y-component
              F1=D0
              F2=D0
              CALL DZERO(WRK(KRRAOx),NNBASX)
              CALL DZERO(WRK(KTRMO),NNORBX)
              CALL DZERO(WRK(KUTR),N2ORBX)

              IF (INTDIR) THEN
                CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                    WRK(KLAST),NBAST,NORBT)
              ELSE
                CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                    NBAST,NORBT)
              ENDIF

              CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

              IF (ISYMT.EQ.ISYMV1) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'QM3QRO')
                 F1 = 1.0D0*F1*FAC

                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)

                 CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)
              END IF

              IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL DZERO(WRK(KTLMB),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
                 CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'QM3QRO')
                 F2 = 0.5D0*F2*FAC

                 CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
              END IF

C z-component
              F1=D0
              F2=D0
              CALL DZERO(WRK(KRRAOx),NNBASX)
              CALL DZERO(WRK(KTRMO),NNORBX)
              CALL DZERO(WRK(KUTR),N2ORBX)

              IF (INTDIR) THEN
                CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &                    WRK(KLAST),NBAST,NORBT)
              ELSE
                CALL READT(LUCCEF,NNBASX,WRK(KRRAOx))
                CALL UTHU(WRK(KRRAOx),WRK(KTRMO),WRK(KUCMO),WRK(KWRK),
     &                    NBAST,NORBT)
              ENDIF 

              CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

              IF (ISYMT.EQ.ISYMV1) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'C3SOL')
                 F1 = 1.0D0*F1*FAC

                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KUTR),WRK(KTLMA),ISYMT)

                 CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)   
              END IF

              IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
                 CALL DZERO(WRK(KTLMA),N2ORBX)
                 CALL DZERO(WRK(KTLMB),N2ORBX)

                 CALL OITH1(ISYMV1,ZYM1,WRK(KUTR),WRK(KTLMA),ISYMT)
                 CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
                 CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'C3SOL')
                 F2 = 0.5D0*F2*FAC

                 CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
              END IF

  522       CONTINUE
  521     CONTINUE
        END IF
  520 CONTINUE

      IF (INTDIR) THEN
        DIPORG(1) = OBKPX
        DIPORG(2) = OBKPY
        DIPORG(3) = OBKPZ
      ENDIF

      IF (.NOT. INTDIR) CALL GPCLOSE(LUCCEF,'KEEP')

  600 CONTINUE ! if lgspol

C       Make the gradient
C
C     / <0| [qj ,TRES] |0> \
C     |          0         |
C     | <0| [qj+,TRES] |0> |
C      \         0         /
C
      ISYMDN = 1
      OVLAP  = D1
      JSPIN = 0
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = NCREF
      NZCVEC = NCREF

      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KTRES),
     *            XINDX,MJWOP,WRK(KWRK),LWRK1,LORB,LCON,LREF)
C
      CALL QEXIT('QM3QRO')
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE QM3QRO_P(VEC1,VEC2,ETRS,XINDX,ZYM1,ZYM2,
     &                 UDV,WRK,LWRK,KZYVR,KZYV1,KZYV2,
     &                 IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP)
C  Parallel version of QM3QRO, Arnfinn, Aug. -08
C  CBN+JK, Dec. 05
C
C  Purpose:  Calculate SCF/DFT E^3 contribution from a surrounding
C            polarizable MM medium to an orbital trial vector.
C
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include <rspprp.h>
#include <infcr.h>
#include <inftap.h>
#include <qrinf.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <ccinftap.h>
#include "nuclei.h"

      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(LWRK),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2)
      DIMENSION MJWOP(2,MAXWOP,8)
      LOGICAL LCON, LORB, LREF, LOPEN
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)

      CALL QENTER('QM3QRO_P')
      IF ( (.NOT. (OLDTG)) .AND. (LOSPC) ) THEN
         CALL QEXIT('QM3QRO_P')
         RETURN
      END IF

      LOPEN = .FALSE.
  
      IF (TRPLET) CALL QUIT('TRIPLET NOT IMPLEMENTED FOR QR-QM/MM YET')

      KCREF   = 1
      KTRES   = KCREF   + NCREF
      KUCMO   = KTRES   + N2ORBX
      KFSOLMO = KUCMO   + NORBT*NBAST
      KUFSOL  = KFSOLMO + NNORBX
      KTLMA   = KUFSOL  + N2ORBX
      KTLMB   = KTLMA   + N2ORBX
      KRRAOx  = KTLMB   + N2ORBX
      KTRMO   = KRRAOx  + NNBASX
      KUTR    = KTRMO   + NNORBX
      KWRK    = KUTR    + N2ORBX
      LWRK1   = LWRK    - KWRK 

      IF (LWRK1 .LT. 0) CALL ERRWRK('QM3QRO_P',-KWRK,LWRK1)

      CALL DZERO(WRK(KCREF),NCREF)
      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL DZERO(WRK(KUCMO),NORBT*NBAST)
      CALL DZERO(WRK(KFSOLMO),NNORBX)
      CALL DZERO(WRK(KUFSOL),N2ORBX)
      CALL DZERO(WRK(KTLMA),N2ORBX)
      CALL DZERO(WRK(KTLMB),N2ORBX)
      CALL DZERO(WRK(KRRAOx),NNBASX)
      CALL DZERO(WRK(KTRMO),NNORBX)
      CALL DZERO(WRK(KUTR),N2ORBX)

      NSIM = 1

C We assume no symmetry in the DFT/MM calculations although we have kept 
C some symmetry options below ...

      ISYMT = 1

C Get the reference state

      CALL GETREF(WRK(KCREF),MZCONF(1))


C Unpack the response vectors

      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)


C Unpack symmetry blocked CMO

      CALL UPKCMO(CMO,WRK(KUCMO))

C--------------------------------------------------------
C 1. Construct two-index transformed T^g.
C    First read qm/mm mo fock matrix from interface file.
C    This includes N_s contribution if OLDTG = .TRUE..
C    If OLDTG = .FALSE. N_s is included in one-electron
C    part of vacuum Hamiltonian
C--------------------------------------------------------

      IF (LUSIFC .LE. 0) THEN
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         LOPEN = .TRUE.
      END IF

      MMORBT = MAX(4,NNORBT)
      IF (MMORBT .NE. NNORBX) THEN
        CALL QUIT('ERROR IN DIMENSION OF QM/MM FOCK CONTR. IN QM3QRO_P')
      END IF
C
      REWIND(LUSIFC)

      CALL MOLLAB('SOLVTMAT',LUSIFC,LUPRI)
      CALL READT (LUSIFC,NNORBX,WRK(KFSOLMO))

      IF (LOPEN) THEN
         CALL GPCLOSE(LUSIFC,'KEEP')
      END IF

      CALL DSPTSI(NORBT,WRK(KFSOLMO),WRK(KUFSOL))

C     Create the effective operator:
C
C     TRES = Tg(k1,k2) + A1(k2) + A12

      IF (ISYMT.EQ.ISYMV1 .AND. ISYMT.EQ.ISYMV2
     &    .AND. ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
          CALL DZERO(WRK(KTLMA),N2ORBX)
          CALL DZERO(WRK(KTLMB),N2ORBX)
          CALL OITH1(ISYMV1,ZYM1,WRK(KUFSOL),WRK(KTLMA),ISYMT)
          CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &              MULD2H(ISYMT,ISYMV1))
          CALL DAXPY(N2ORBX,0.5D0,WRK(KTLMB),1,WRK(KTRES),1)
      END IF

C     If LGSPOL neglect the time-dependece of the expectation value 
C     contained in the interaction operator, i.e. keep the ground state polarization only. 

      IF (LGSPOL) GOTO 600


      IF (.NOT. INTDIR) THEN
        CALL GPOPEN(LUCCEF,'ELFDMM','OLD',' ',
     &                 'UNFORMATTED',IDUMMY,.FALSE.)
        REWIND(LUCCEF)
      ENDIF

C     Readin Rra integrals (ao) and transform to mo
C     ---------------------------------------------

      LM = 0

      IF (INTDIR) THEN
        L = NUSITE + NUALIS(0)
        OBKPX = DIPORG(1)
        OBKPY = DIPORG(2)
        OBKPZ = DIPORG(3)
      ENDIF 

      DO 520 I = 1, ISYTP
         IF (MDLWRD(I)(1:5) .EQ. 'SPC_E') THEN
            CALL QM3QRO_M(WRK(KTRES),IATNOW,NOSIM,WRK(KWRK),LWRK1,I,
     &           ISYMT,LM,L,WRK(KUCMO),.FALSE.,ZYM1,ZYM2,ISYMV1,
     &           ISYMV2)
         END IF
 520  CONTINUE

      IF (INTDIR) THEN
         DIPORG(1) = OBKPX
         DIPORG(2) = OBKPY
         DIPORG(3) = OBKPZ
      ENDIF

      IF (.NOT. INTDIR) CALL GPCLOSE(LUCCEF,'KEEP')

  600 CONTINUE ! if lgspol

C       Make the gradient
C
C     / <0| [qj ,TRES] |0> \
C     |          0         |
C     | <0| [qj+,TRES] |0> |
C      \         0         /
C
      ISYMDN = 1
      OVLAP  = D1
      JSPIN = 0
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = NCREF
      NZCVEC = NCREF

      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KTRES),
     *            XINDX,MJWOP,WRK(KWRK),LWRK1,LORB,LCON,LREF)
C
      CALL QEXIT('QM3QRO_P')
      RETURN
      END
C 
C
C---------------------------------------------------------------------
      SUBROUTINE QM3QRO_M(TRES,IATNOW,NOSIM,WRK,LWRK,INUM,ISYMT,
     &                   LM,LNUM,UCMO,TOFILE,ZYM1,ZYM2,ISYMV1,
     &                   ISYMV2)
C
C  Arnfinn Aug. -08
C  The master routine for QM3QRO_P
C
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include <rspprp.h>
#include <infcr.h>
#include <inftap.h>
#include <qrinf.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <ccinftap.h>
#include "nuclei.h"
#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif

      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION ZYM1(*), ZYM2(*)
      INTEGER IATNOW, NOSIM
      DIMENSION WRK(LWRK)
      DIMENSION UCMO(NORBT*NBAST),TRES(N2ORBX)
      DIMENSION ISX(MXCORB),IPATOM(MXCENT)

      CALL QENTER('QM3QRO_M')

      IF (TOFILE) CALL QUIT('Parallel calculations do not allow '//
     &                     'for storing integrals on disk')

      IPRTYP = 18
      KTRES = 1
      JTRES = KTRES + N2ORBX
      KLAST = JTRES + N2ORBX

      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRRSP,1,'INTEGER',MASTER)

      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)

      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LNUM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(UCMO,NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(ZYM1,NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ZYM2,NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)

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

      CALL DZERO(WRK(KTRES),N2ORBX)
      CALL MPI_REDUCE(WRK(KTRES),WRK(JTRES),N2ORBX,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,
     &                   IERR)

      FAC = 1.0D0
      CALL DAXPY(N2ORBX,FAC,WRK(JTRES),1,TRES,1)

      CALL QEXIT('QM3QRO_M')

      RETURN
      END
C---------------------------------------------------------------
      SUBROUTINE QM3QRO_S(WRK,LWRK,IPRTMP)
C
C  Arnfinn Aug. -08
C  The slave routine for QM3QRO
C
#include <implicit.h>
#include <dummy.h>
#include <maxorb.h>
#include <gnrinf.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include <rspprp.h>
#include <infcr.h>
#include <inftap.h>
#include <qrinf.h>
#include <qm3.h>
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <orgcom.h>
#include <ccinftap.h>
#include "nuclei.h"

#include <infpar.h>
#include <mtags.h>
#if defined(VAR_MPI)
#include <mpif.h>
#endif

      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, EXP1VL, TRIMAT
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(LWRK)
      INTEGER IATNOW, NOSIM
      DIMENSION ISX(MXCORB), IPATOM(MXCENT)

      CALL QENTER('QM3QRO_S')
      IPRRSP = IPRTMP

      CALL MPIXBCAST(NNBASX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NNORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(N2ORBX,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NORBT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NBAST,1,'INTEGER',MASTER)

      QM3 = .TRUE.

      KRRAOx = 1
      KTRMO  = KRRAOx + NNBASX
      KUTR   = KTRMO  + NNORBX
      KUTRX  = KUTR   + N2ORBX
      JRXYO  = KUTRX  + N2ORBX
      KUCMO  = JRXYO  + N2ORBX
      KTLMA  = KUCMO  + NORBT*NBAST
      KTLMB  = KTLMA  + N2ORBX
      KTRES  = KTLMB  + N2ORBX
      KZYM1  = KTRES  + N2ORBX
      KZYM2  = KZYM1  + NORBT*NORBT
      KLAST1 = KZYM2  + NORBT*NORBT

      CALL DZERO(WRK(KTRES),N2ORBX)

      CALL MPIXBCAST(ISYTP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NCTOT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IQM3PR,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TOFILE,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(INTDIR,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(FORQM3,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NUCIND,1,'INTEGER',MASTER)
      CALL MPIXBCAST(L,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NOSIM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(WRK(KUCMO),NORBT*NBAST,'DOUBLE',MASTER)
      CALL MPIXBCAST(CORD,3*MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(CHARGE,MXCENT,'DOUBLE',MASTER)
      CALL MPIXBCAST(NUCDEP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(LUPROP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(NPATOM,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPATOM,MXCENT,'INTEGER',MASTER)
      CALL MPIXBCAST(ISX,MXCORB,'INTEGER',MASTER)
      CALL MPIXBCAST(KSYMOP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(TRPLET,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(WRK(KZYM1),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(WRK(KZYM2),NORBT*NORBT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ISYMV1,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMV2,1,'INTEGER',MASTER)
      CALL MPIXBCAST(ISYMT,1,'INTEGER',MASTER)

      CALL DZERO(WRK(JRXYO),N2ORBX)

 187  CONTINUE

      CALL MPIXSEND(MYNUM, 1, 'INTEGER', MASTER, MPTAG1)
      CALL MPIXRECV(LM, 1, 'INTEGER', MASTER, MPTAG2)
      CALL MPIXRECV(FAC, 1, 'DOUBLE', MASTER, MPTAG2)

      IF (LM .GT. 0) THEN

C x-component
         F1=D0
         F2=D0
         CALL DZERO(WRK(KRRAOx),NNBASX)
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         KMAT = KLAST1
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
     &           LWRK2,LABINT,INTREP,INTADR,IATNOW,TOFILE,
     &           KPATOM,TRIMAT,DUMMY,EXP1VL,DUMMY,IPRRSP)
         RUNQM3 = .FALSE.
         CALL UTHU(WRK(KMAT),WRK(KTRMO),WRK(KUCMO),WRK(KLAST),
     &           NBAST,NORBT)

         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))
C
         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'QM3QRO')
            F1 = 1.0D0*F1*FAC

            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)
         END IF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)

            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'C3SOL')
            F2 = 0.5D0*F2*FAC

            CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
         END IF

C y-component
         F1=D0
         F2=D0
         CALL DZERO(WRK(KRRAOx),NNBASX)
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+NNBASX),WRK(KTRMO),WRK(KUCMO),
     &           WRK(KLAST),NBAST,NORBT)

         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)

            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'QM3QRO')
            F1 = 1.0D0*F1*FAC

            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)

            CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)
         END IF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)

            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'QM3QRO')
            F2 = 0.5D0*F2*FAC

            CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
         END IF

C z-component
         F1=D0
         F2=D0
         CALL DZERO(WRK(KRRAOx),NNBASX)
         CALL DZERO(WRK(KTRMO),NNORBX)
         CALL DZERO(WRK(KUTR),N2ORBX)

         CALL UTHU(WRK(KMAT+2*NNBASX),WRK(KTRMO),WRK(KUCMO),
     &           WRK(KLAST),NBAST,NORBT)

         CALL DSPTSI(NORBT,WRK(KTRMO),WRK(KUTR))

         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)

            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F1,200,'C3SOL')
            F1 = 1.0D0*F1*FAC

            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KUTR),WRK(KTLMA),ISYMT)

            CALL DAXPY(N2ORBX,F1,WRK(KTLMA),1,WRK(KTRES),1)   
         END IF

         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)

            CALL OITH1(ISYMV1,WRK(KZYM1),WRK(KUTR),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,WRK(KZYM2),WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,D1,F2,200,'C3SOL')
            F2 = 0.5D0*F2*FAC

            CALL DAXPY(N2ORBX,F2,WRK(KUTR),1,WRK(KTRES),1)
         END IF
         GO TO 187
      END IF


C     No more integrals to calculate
      CALL MPI_REDUCE(WRK(KTRES),WRK(KTRES),N2ORBX,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

      CALL QEXIT('QM3QRO_S')

      RETURN
      END
C-----------------------------------------------------------------------
C  /* Deck qm3qlm */
      FUNCTION QM3QLM(RLM)
C
#include <implicit.h>
      DIMENSION RLM(NORBT,*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
C Used from common blocks:
C  INFORB: NISHT,NASHT
C  INFIND: IROW(*), ISX(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <priunit.h>

      TELM  = D0

      DO 300 IW = 1,NISHT
         I = ISX(IW)
         TELM = TELM + D2*RLM(I,I)
  300 CONTINUE

      QM3QLM = TELM
      RETURN
      END

