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
C  /* Deck orblin */
      SUBROUTINE ORBLIN(NOSIM,BOVECS,CMO,GORB,DV,PV,FC,FV,
     &                  WRK,KFRSAV,LFRSAV)
C
C December 1989 Hans Joergen Aa. Jensen
C
C
#include <implicit.h>
      DIMENSION BOVECS(*),CMO(*),GORB(*),DV(*),PV(*),FC(*),FV(*),
     &          WRK(*)
#include <dummy.h>
C
      CALL QENTER('ORBLIN')
      NCSIM = -1
      CALL LINTRN(NCSIM,NOSIM,DUMMY,BOVECS,
     *            CMO,DUMMY,DUMMY,GORB,DV,PV,
     *            FC,FV,DUMMY,DUMMY,DUMMY,WRK,KFRSAV,LFRSAV)
C     CALL LINTRN(NCSIM,NOSIM,BCVECS,BOVECS,
C    *            CMO,CREF,EACTIV,GORB,DV,PV,
C    *            FC,FV,FCAC,H2AC,INDXCI,WRK,KFRSAV,LFRSAV)
      CALL QEXIT('ORBLIN')
      RETURN
      END
C  /* Deck lintrn */
      SUBROUTINE LINTRN(NCSIM,NOSIM,BCVECS,BOVECS,
     *                  CMO,CREF,EACTIV,GORB,DV,PV,
     *                  FC,FV,FCAC,H2AC,INDXCI,WRK,KFRSAV,LFRSAV)
C
C Rewritten Jan 1990 hjaaj
C Originally written Nov-1983 by Hans Jorgen Aa. Jensen in Uppsala
C Revisions:
C   Feb 90: entry point LINMEM
C
C MOTECC-90: The purpose of this module, LINTRN, and the algorithms used
C            are described in Chapter 8 Section B.5 of MOTECC-90
C            "The Direct Iterative NEO Algorithm"
C
C Purpose:
C  Calculate the simultaneous linear transformation defined by
C  the K matrix of the SIM(ultaneous) B vectors :
C
C  SCVECS(i,j)  =  sum(k=1,NCONF) K(i,k) * BCVECS(k,j)
C                  (i = 1,NVAR; j = 1,NCSIM)
C
C  SOVECS(i,j)  =  sum(k=1,NWOPT) K(i,NCONF+k) * BOVECS(k,j)
C                  (i = 1,NVAR; j = 1,NOSIM)
C
C NOTE: 1) all BCVECS are assumed orthogonal to CREF;
C          the linear transformation of CREF can be calculated by GRAD.
C       2) NCSIM .lt. 0 flags only sigma vectors of orbital Hessian
C          wanted (frozen CI coefficients), used by ORBLIN. Then
C          no SCVECS and SOVECS will be:
C          SOVECS(i,j) = sum(k=1,NWOPT) K(NCONF+i,NCONF+k) * BOVECS(k,j)
C                        (i = 1,NWOPT; j = 1,NOSIM)
C       3) No guarantee for correct CREF part in sigma vectors.
C
C Input:
C     GORB(NWOPH)   total orbital gradient
C
C Output:
C  SCVECS; will start at WRK(1) on output; followed by
C  SOVECS; if NCSIM .ge. 0 then SOVECS(NVAR,NOSIM),
C                               starting at WRK(1+NCSIM*NVAR)
C                          else SOVECS(NWOPT,NOSIM),
C                               starting at WRK(1)
C
#include <implicit.h>
      DIMENSION BCVECS(NCONDI,*), BOVECS(NWOPDI,*)
      DIMENSION CREF(*), GORB(*), CMO(*), DV(*), PV(*)
      DIMENSION FC(*), FV(*), FCAC(*), H2AC(*)
      DIMENSION INDXCI(*),WRK(*)
#include <iratdef.h>
C
C *** local constants
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C   INFINP : FLAG(*),LSYM
C   INFVAR : NCONF,NWOPT,NVAR,NWOPH,JWOP(2,*),JWOPSY
C   INFORB : N2ORBX,...
C   INFLIN : NCONRF,NCONST,NWOPPT,NVARPT,LSYMRF
C   INFDIM : NCONDI,NWOPDI,LACIMX,LBCIMX,?
C   INFPRI : P6FLAG(*), ?
C CBGETDIS : DISTYP,IADINT,IADH2,IADH2X
C   INFTIM : NCALLS,TIMCPU,TIMWAL  ! IDTIM is index for these
C
#include <maxorb.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <inflin.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
#include <cbgetdis.h>
#include <dftcom.h> 
! gnrinf.h: QM3 /CBN+JK 03.01.06
#include <gnrinf.h>
#include <wrkrsp.h>
C ---
      PARAMETER (IDTIM1 = 6, IDTIM2 = 8, D4 = 4.0D0)
#include <inftim.h>
C
      LOGICAL SOLVNT, ORBLIN, MEMLIN, DFTADX
C
C
      CALL QENTER('LINTRN')
      KFREE = KFRSAV
      LFREE = LFRSAV
      MEMLIN = .FALSE.
      GO TO 10
         ENTRY LINMEM(NCSIM,NOSIM,LNEED)
         KFREE  = 1
         LFREE  = 0
         MEMLIN = .TRUE.
   10 CONTINUE
C
C
      SOLVNT = FLAG(16)
C
C     Check for some parameter errors
C
      NUMERR = 0
C     check INFDIM parameters
      IF (NCONF .NE. NCONDI .AND. NCONF .GT. 0) NUMERR = NUMERR + 1
      IF (NWOPT .NE. NWOPDI .AND. NWOPT .GT. 0) NUMERR = NUMERR + 1
C     check INFLIN parameters
      IF (LSYM  .NE. LSYMRF) NUMERR = NUMERR + 1
      IF (JWOPSY.NE. LSYMPT) NUMERR = NUMERR + 1
      IF (NCONF .NE. NCONST) NUMERR = NUMERR + 1
      IF (NWOPT .NE. NWOPPT) NUMERR = NUMERR + 1
      IF (NVAR  .NE. NVARPT) NUMERR = NUMERR + 1
      IF (NUMERR .GT. 0) THEN
         WRITE (LUPRI,'(//A/I8,A/)') ' >>> FATAL ERROR (LINTRN)',
     *      NUMERR,' dimension error(s) in INFDIM and/or INFLIN'
         WRITE (LUPRI,'(A,2I8)')
     &      'NCONF ,NCONDI :',NCONF,NCONDI,
     &      'NWOPT ,NWOPDI :',NWOPT,NWOPDI,
     &      'LSYM  ,LSYMRF :',LSYM,LSYMRF,
     &      'JWOPSY,LSYMPT :',JWOPSY,LSYMPT,
     &      'NCONF ,NCONST :',NCONF,NCONST,
     &      'NWOPT ,NWOPPT :',NWOPT,NWOPPT,
     &      'NVAR  ,NVARPT :',NVAR,NVARPT
         CALL QTRACE(LUPRI)
         CALL QUIT('FATAL ERROR (LINTRN) INFDIM or '//
     &             'INFLIN DIMENSION ERROR')
      END IF
C
C *** Allocate work area for SIRTR1,ORBSIG,CISIG
C
      IF (NCSIM.GE.0) THEN
         ORBLIN = .FALSE.
         IDTIM = IDTIM1
         MCSIM = NCSIM
         MOSIM = NOSIM
      ELSE
         ORBLIN = .TRUE.
         IDTIM = IDTIM2
         MCSIM = 0
         MOSIM = 0
      END IF
      IF (MCSIM + NOSIM .LE. 0) THEN
         WRITE (LUPRI,*) 'LINTRN ERROR: NCSIM + NOSIM .le. 0'
         CALL QTRACE(LUPRI)
         CALL QUIT('ERROR (LINTRN) NCSIM + NOSIM .le. 0')
      END IF
      CALL GETTIM(T0,W0)
C-- ORBSIG, CISIG
C MAERKE flyt SCVEC/SOVEC efter SIRTR1 variable
C MAERKE og CALL DCOPY((KSEND+1-KSCVEC),WRK(KSCVEC),1,WRK,1)
C MAERKE til sidst for at flytte sigma vektorerne forrest;
C MAERKE derved behoeves PTV og SVECS ikke i core samtidigt/891130
      KSCVEC = KFREE
      KSOVEC = KSCVEC + MCSIM*NVAR
      IF (.NOT.ORBLIN) THEN
         KSEND  = KSOVEC + NOSIM*NVAR
      ELSE
         KSEND  = KSOVEC + NOSIM*NWOPT
      END IF
C-- SIRTR1, ..., ORBSIG, SOLLIN
      KDTV   = KSEND
      KWSOL  = KDTV  + MCSIM*NNASHX
      LWSOL  = LFREE - KWSOL
C-- SIRTR1, ..., CISIG
      KEMYX  = KWSOL
C-- SIRTR1, ..., ORBSIG
      KFXC   = KEMYX  + MOSIM
      KFXV   = KFXC   + NOSIM*N2ORBX
      KFTV   = KFXV   + NOSIM*N2ORBX
C-- SIRTR1, ORBSIG, CISIG
      LFTV   = MCSIM*N2ORBX
      IF (DOMCSRDFT) LFTV = 2*LFTV
C     ... Extra allocation for "FTC" in MCSCF-SRDFT
C         Will be constructed in TR1FCY
      KH2XAC = KFTV   + LFTV
      KFXQ   = KH2XAC + MOSIM*NNASHX*NNASHX
      KFTQ   = KFXQ   + NOSIM*NASHT*NORBT
C-- SIRTR1, ORBSIG
      K3     = KFTQ   + MCSIM*NASHT*NORBT
      L3     = LFREE  - K3
      KUBO   = K3
      K3A    = KUBO  + NOSIM*N2ORBX
      L3A    = LFREE - K3A
C-- CISIGO/CISIGC
      IF ( .NOT.ORBLIN ) THEN
         KFXCAC = K3
         KCW    = KFXCAC + MOSIM*NNASHX
         LCW    = LFREE  - KCW
         MONE   = LACIMX + MOSIM*LBCIMX
         MTWO   = LACIMX + MCSIM*LBCIMX
         K4     = KCW    + MAX(MONE,MTWO)
      ELSE
         KFXCAC = K3
         KCW    = K3
         LCW    = 0
         K4     = K3
      END IF
C--
      KEND   = MAX(K3A,K4) - 1
      IF (MEMLIN) THEN
         LNEED = KEND
         RETURN
      END IF
      IF (KEND.GT.LFREE) CALL ERRWRK('LINTRN',KEND,LFREE)
      CALL GETTIM(T1,W1)
C
C *** If CI (i.e. NWOPT = 0) call special routine
C     to calculate sigma vector
C
C     KDTV need to be zeroed out in case of FCI solvent calculation
C     (DTV only set in sirtr1). K.Ruud-July 97
C
      CALL DZERO(WRK(KDTV),MCSIM*NNASHX)
      IF (NWOPT.EQ.0.AND. .NOT.DOMCSRDFT) THEN
         T2 = T1
         T3 = T1
         T4 = T1
         T5 = T1
         W2 = W1
         W3 = W1
         W4 = W1
         W5 = W1
         DISTYP = 1
         IADINT = IADH2
         CALL CISIGC(NCSIM,BCVECS,WRK(KSCVEC),NCONF,FCAC,H2AC,
     *               INDXCI,WRK(KCW),LCW)
C        CALL CISIGC(NSIM,BCVECS,SCVECS,LSCVEC,FCAC,H2AC,INDXCI,
C    *               WRK,LFREE)
         DO 60 ICSIM = 1,NCSIM
            ISC0 = KSCVEC - 1 + (ICSIM-1)*NCONF
            DO 50 I = 1,NCONF
               WRK(ISC0+I) = D2 * (WRK(ISC0+I) - EACTIV*BCVECS(I,ICSIM))
   50       CONTINUE
   60    CONTINUE
         GO TO 900
      END IF
C
C
C *** Calculate various transformed Fock matrices and H2XAC
C
C
C
      DO 200 IOSIM = 1,NOSIM
         JUBO = KUBO + (IOSIM-1)*N2ORBX
         CALL UPKWOP(NWOPT,JWOP,BOVECS(1,IOSIM),WRK(JUBO))
  200 CONTINUE
      DFTADX = DFTADD
      DFTADD = .FALSE.
      CALL SIRTR1(MCSIM,NOSIM,CMO,
     &            BCVECS,CREF,INDXCI,WRK(KDTV),WRK(KFTV),WRK(KFTQ),
     &            WRK(KUBO),DV,PV,FC,FV,
     &            WRK(KFXC),WRK(KFXV),WRK(KFXQ),WRK(KEMYX),WRK(KH2XAC),
     &            WRK,K3A,L3A,ORBLIN,.FALSE.)
      DFTADD = DFTADX
C     CALL SIRTR1(NCSIM,NOSIM,CMO,
C    &            BCVECS,CREF,INDXCI,DTV,FTV,FTQ,
C    &            UBO,DV,PV,FC,FV,FXC,FXV,FXQ,EMYX,H2XAC,
C    &            WRK,KFRSAV,LFRSAV,ORBLIN,ADDTR1)
C
C ***
C
      CALL GETTIM(T2,W2)
C
C --- initialize all sigma vectors to zero
C     (We have now finished with PTV)
C
      CALL DZERO(WRK(KSCVEC),(KSEND-1))
C
C --- calculate orbital part of sigma vector(s)
C
      CALL ORBSIG(NCSIM,NOSIM,WRK(KSCVEC),WRK(KSOVEC),
     *            WRK(KDTV),FC,WRK(KFTV),WRK(KFTQ),
     *            BOVECS,DV,WRK(KFXC),WRK(KFXV),WRK(KFXQ),
     *            GORB,WRK(K3),L3)
C     CALL ORBSIG(NCSIM,NOSIM,SCVECS,SOVECS,
C    *            DTV,FC,FTV,FTQ,BOVECS,DV,FXC,FXV,FXQ,
C    *            GORB,WRK,LWRK)
C
C
C ***
C
      CALL GETTIM(T3,W3)
      IF (.NOT.ORBLIN) THEN
         IF (NCSIM .GT. 0) THEN
            DISTYP = 1
            IADINT = IADH2
            CALL CISIGC(NCSIM,BCVECS,WRK(KSCVEC),NVAR,FCAC,H2AC,
     *                  INDXCI,WRK(KCW),LCW)
C           CALL CISIGC(NSIM,BCVECS,SCVECS,LSCVEC,FCAC,H2AC,
C    *                  INDXCI,WRK,LFREE)
            DO 540 ICSIM = 1,NCSIM
               JSCVEC = (KSCVEC-1) + (ICSIM-1)*NVAR
               DO 520 I = 1,NCONF
                  WRK(JSCVEC+I) =
     *               D2 * ( WRK(JSCVEC+I) - EACTIV * BCVECS(I,ICSIM) )
  520          CONTINUE
  540       CONTINUE
            IF (DOMCSRDFT) THEN
C              ... special MCSRDFT contribution to csf sigma-vectors,
C                  stored in FTV(1,NCSIM+ICSIM) in SIRTR1 after the 
C                  standard FTV matrices.
C                  In paper: V^([2c]xc-SR) matrix.
               ESRLTR   = D0
C              ... expectation value is not zero, but it doesn't matter
C                  because CREF is projected out of SCVEC later
               KSRCVEC  = KCW
               KSRLTRAC = KSRCVEC + NCONF
               KSRW     = KSRLTRAC + NNASHX
               LSRW     = LFREE - KSRW
               DO ICSIM = 1,NCSIM
                  JFTC  = KFTV + (NCSIM+ICSIM-1)*N2ORBX
                  JSCVEC = KSCVEC + (ICSIM-1)*NVAR
                  CALL GETAC1(WRK(JFTC),WRK(KSRW))
                  CALL DGETSP(NASHT,WRK(KSRW),WRK(KSRLTRAC))
                  IF (IPRLIN .GT. 10) THEN
                     WRITE(LUPRI,'(/A,I3)')
     &               'MCSRDFT "FTC"="V^([2c],xc-SR)" no.',ICSIM
             CALL OUTPUT(WRK(JFTC),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                     WRITE(LUPRI,'(/A,I3)') 'MCSRDFT SRLTRAC no.',ICSIM
                     CALL OUTPAK(WRK(KSRLTRAC),NASHT,1,LUPRI)
                  END IF
                  CALL SOLGC(CREF,WRK(KSRLTRAC),ESRLTR,WRK(KSRCVEC),
     &                       INDXCI,WRK(KSRW),LSRW)
                  CALL DAXPY(NCONF,D1,WRK(KSRCVEC),1,WRK(JSCVEC),1)
               ENDDO
            END IF
         END IF
         IF (NOSIM .GT. 0 .AND. NCONF .GT. 0) THEN
         IF (NCONF .GT. 1 .OR. LSYMPT .NE. 1) THEN
C
C
C *** Extract FXCAC from FXC and
C     modify  FXCAC, H2XAC as needed by CI routines
C
            DO 600 IOSIM = 1,NOSIM
               JFXC   = KFXC   + (IOSIM-1)*N2ORBX
               JFXCAC = KFXCAC + (IOSIM-1)*NNASHX
               JH2XAC = KH2XAC + (IOSIM-1)*NNASHX*NNASHX
               CALL GETAC1(WRK(JFXC),WRK(KCW))
               CALL DGETSP(NASHT,WRK(KCW),WRK(JFXCAC))
               IF (P6FLAG(22) .OR. IPRLIN .GE. 22) THEN
                  WRITE (LUPRI,6010) IOSIM,NOSIM
                  CALL OUTPAK(WRK(JFXCAC),NASHT,1,LUPRI)
               END IF
  600       CONTINUE
 6010       FORMAT (/' Modified FXCAC matrix (no.',I3,' of',I3,')')
C
            DISTYP = 1
            IADINT = IADH2X
            CALL CISIGO(NOSIM,WRK(KSOVEC),CREF,WRK(KEMYX),
     *                  WRK(KFXCAC),WRK(KH2XAC),INDXCI,WRK(KCW),LCW)
C           CALL CISIGO(NOSIM,SOVECS,CREF,EMYX,FXCAC,H2XAC,
C    *                  INDXCI,WRK,LFREE)
         ELSE
            DO 800 IOSIM = 1,NOSIM
               WRK(KSOVEC+(IOSIM-1)*NVAR) = D0
  800       CONTINUE
         END IF
         END IF
      END IF
C
C *** end of electronic contribution to LINTRN
C
  900 CONTINUE
      CALL GETTIM(T4,W4)
      DFTADX = DFTADD
      DFTADD = .TRUE.
      IF (DFTRUN) THEN
         IF(NASHT.NE.0) CALL QUIT('QC-SCF not implemented for RO-DFT')
         DO IOSIM = 1,NOSIM
           CALL DZERO(WRK(KFXC),N2ORBX)
           CALL UPKWOP(NWOPT,JWOP,BOVECS(1,IOSIM),WRK(KUBO))
           KSYMOX = KSYMOP
           KSYMOP = 1
           CALL DFT_LIN_RESPF(1,WRK(KFXC),CMO,WRK(KUBO),
     &                        .FALSE.,KSYMOP,WRK(K3A),L3A)
           KSYMOP = KSYMOX
           CALL PKWOP(NWOPT,JWOP,WRK(KUBO),WRK(KFXC))
           IF (IPRLIN.GT.20) THEN
              WRITE (LUPRI,'(/A)') 'Orbital part of linearly'//
     &             ' transformed B vector before DFT contribution'
              IF (IPRLIN.GT.100) THEN
                 PRFAC = 0.0D0
              ELSE
                 PRFAC = 0.1D0
              END IF
              CALL PRKAP (NWOPT,WRK(KSOVEC),PRFAC,LUPRI)
           END IF
           CALL DAXPY(NWOPT,D4,WRK(KUBO),1,WRK(KSOVEC+(IOSIM-1)*NVAR),1)
         END DO
      END IF
      DFTADD = DFTADX 
C
C *** Now the solvent contribution, if any
C
      IF (SOLVNT) THEN
         CALL SOLLIN(NCSIM,NOSIM,BCVECS,BOVECS,CREF,CMO,INDXCI,
     &               DV,WRK(KDTV),WRK(KSCVEC),WRK(KSOVEC),ORBLIN,
     &               WRK(KWSOL),LWSOL)
         CALL GETTIM(T5,W5)
      ELSEIF (PCM) THEN
         MXSIM = MAX(NOSIM,NCSIM)
         KVTEX = KWSOL
         KWPCM = KVTEX + MXSIM*NTS
         LWPCM = LWSOL - MXSIM*NTS
         CALL PCMLIN(NCSIM,NOSIM,BCVECS,BOVECS,CREF,CMO,INDXCI,
     &               DV,WRK(KDTV),WRK(KSCVEC),WRK(KSOVEC),ORBLIN,
     &               WRK(KVTEX),WRK(KWPCM),LWPCM)
         CALL GETTIM(T5,W5)
      ELSE
         T5 = T4
         W5 = W4
      END IF
C
C---------------------------
C CBN 03.01.06
C---------------------------
C
      IF (QM3) THEN
         CALL QM3LIN(NOSIM,BOVECS,CREF,CMO,INDXCI,
     &               WRK(KSOVEC),ORBLIN,WRK(K3A),L3A)
         CALL GETTIM(T5,W5)
      END IF
C---------------------------
C CBN 03.01.06
C---------------------------
C
C
      NCALLS(IDTIM)   = NCALLS(IDTIM)   + 1
      TIMCPU(1,IDTIM) = TIMCPU(1,IDTIM) + T5 - T0
      TIMCPU(2,IDTIM) = TIMCPU(2,IDTIM) + T2 - T1
      TIMCPU(3,IDTIM) = TIMCPU(3,IDTIM) + T3 - T2
      TIMCPU(4,IDTIM) = TIMCPU(4,IDTIM) + T4 - T3
      TIMCPU(5,IDTIM) = TIMCPU(5,IDTIM) + T5 - T4
      TIMWAL(1,IDTIM) = TIMWAL(1,IDTIM) + W5 - W0
      TIMWAL(2,IDTIM) = TIMWAL(2,IDTIM) + W2 - W1
      TIMWAL(3,IDTIM) = TIMWAL(3,IDTIM) + W3 - W2
      TIMWAL(4,IDTIM) = TIMWAL(4,IDTIM) + W4 - W3
      TIMWAL(5,IDTIM) = TIMWAL(5,IDTIM) + W5 - W4
C
C
      CALL QEXIT('LINTRN')
      RETURN
C     end of lintrn.
      END
C  /* Deck orbsig */
      SUBROUTINE ORBSIG(NCSIM,NOSIM,SCVECS,SOVECS,
     *                  DTV,FC,FTV,FTQ,BOVECS,DV,FXC,FXV,FXQ,
     *                  GORB,WRK,LWRK)
C
C 30-Nov-1989 hjaaj
C
C New interface to ORBSIG, old ORBSIG now ORBSI2
C
C MOTECC-90: The algorithms used in this module, ORBSIG, are
C            described in Chapter 8 Appendix 8A of MOTECC-90
C            "Calculation of the Orbital Part of the Gradient
C            Vectors and Sigma Vectors"
C
C
#include <implicit.h>
      DIMENSION SCVECS(*), SOVECS(*)
      DIMENSION DTV(*), FC(*), FTV(*), FTQ(*)
      DIMENSION BOVECS(NWOPDI,*), DV(*), FXC(*), FXV(*), FXQ(*)
      DIMENSION GORB(*), WRK(LWRK)
C
C Used from common blocks:
C   INFORB : N2ORBX,...
C   INFVAR : NVAR, JWOP
C   INFDIM : NWOPDI
C   INFLIN : IPRLIN
C   INFPRI : P6FLAG()
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infvar.h>
#include <infdim.h>
#include <inflin.h>
#include <infpri.h>
#include <infinp.h>
C
      LOGICAL ORBLIN
C
      CALL QENTER('ORBSIG')
      ORBLIN = NCSIM.LT.0
      KUBO   = 1
      IF (NOSIM .GT. 0) THEN
         KUGORB = KUBO   + NOSIM*N2ORBX
         KUDV   = KUGORB + N2ORBX
      ELSE
         KUGORB = KUBO
         KUDV   = KUGORB
      END IF
      KW1    = KUDV   + N2ASHX
      IF (KW1 .GT. LWRK) CALL ERRWRK('ORBSIG',KW1,LWRK)
C
C     first unpack B-orb to antisymmetric matrix UBO  for
C       one-index transformations:
C
      IF (NOSIM.GT.0) THEN
         DO 210 IOSIM = 1,NOSIM
            JUBO = KUBO + (IOSIM-1)*N2ORBX
            CALL UPKWOP(NWOPT,JWOP,BOVECS(1,IOSIM),WRK(JUBO))
#if defined (VAR_DEBUG)
            IF (P6FLAG(19) .OR. IPRLIN .GE. 19) THEN
               WRITE (LUPRI,2110) IOSIM,NOSIM
               CALL OUTPUT(WRK(JUBO),1,NORBT,1,NORBT,NORBT,NORBT,
     &                     1,LUPRI)
            END IF
 2110 FORMAT (/,' ORBSIG: Orbital trial vector unpacked to ',
     *        'matrix form (no.',I3,' of',I3,')')
#endif
  210    CONTINUE
         CALL UPKWOP(NWOPH,JWOP,GORB,WRK(KUGORB))
      END IF
      CALL ORBSI2(NCSIM,NOSIM,SCVECS,SOVECS,
     *            DTV,FC,FTV,FTQ,DV,FXC,FXV,FXQ,
     *            WRK(KUBO),WRK(KUGORB),WRK(KUDV))
C
C     Special MCSRDFT correction, corresponding to the effective operator
C        V^([2c],xc-SR)
C     is saved in FTV(1,NCSIM+ICSIM) in SIRTR1.
C
      IF (DOMCSRDFT .AND. NCSIM .GT. 0) THEN
#ifndef MOD_SRDFT
         call quit('srdft not included in this version')
#else
         CALL DSPTSI(NASHT,DV,WRK(KUDV))
         DO ICSIM = 1,NCSIM
            JSCVEC = (NCONF + 1) + (ICSIM-1)*NVAR
            JFSR   = 1 + (NCSIM+ICSIM-1)*N2ORBX
            CALL SRDFTSO(2.0D0,WRK(KUDV),FTV(JFSR),SCVECS(JSCVEC))
         END DO
#endif
      END IF
C
C     *** Make microcanonical average on SVECS
C         (AVERAG returns immediately if no averageing requested)
C
      IF (NCSIM .GT. 0) THEN
         CALL AVERAG(SCVECS(1+NCONF),NVAR,NCSIM)
      END IF
      IF (NOSIM .GT. 0) THEN
         IF (ORBLIN) THEN
            CALL AVERAG(SOVECS,NWOPT,NOSIM)
         ELSE
            CALL AVERAG(SOVECS(1+NCONF),NVAR,NOSIM)
         END IF
      END IF
      CALL QEXIT('ORBSIG')
      RETURN
      END
C  /* Deck orbsi2 */
      SUBROUTINE ORBSI2(NCSIM,NOSIM,SCVECS,SOVECS,
     *                  DTV,FC,FTV,FTQ,DV,FXC,FXV,FXQ,
     *                  UBO,UGORB,UDV)
C
C Written 29-Nov-1983 by Hans Agren and Hans Jorgen Aa. Jensen
C Revisions:
C  10-Oct-1984 / 2-May-1984 hjaaj
C  29-Dec-1984 hjaaj  (skip CI-orbital coupling contribution
C                      if NCONF has been set .eq. 0;
C                      removed redundant "2. ind. act." sct.)
C   7-Jan-1985 hjaaj  (removed BOVLP vector)
C  13-Jan-1985 hjaaj  (corrected act-act rotations)
C  30 Jan  '85 hjaaj  (ver 2; separated conf. and orb. trial vectors)
C  31-Oct-1989 hjaaj:  UBO(NORBT,NORBT) instead of UTBO(N2ORBT)
C                      UGORB(NORBT,NORBT) instead of UGORB(N2ORBT)
C
C Purpose:
C
C   To construct the orbital part of the sigma vectors;
C
C   If NCSIM.lt.0 (used by ORBLIN),
C   output is the orbital Hessian sigma vectors.
C
C Input:
C   The active one-electron density matrix DV (packed) and
C     transition density matrix DTV
C   The inactive and active transition Fock matrices FC,FTV
C   The inactive and active 1-index transformed Fock matrices FXC,FXV
C   The unpacked gradient and kappa vectors UGORB and UBO.
C   The Fock "Q" matrices in FXQ and FTQ
C
C Earlier version:
C   The overlap between B(*,icsim) and CREF(*) in BOVLP(icsim),
C   now they are assumed to be orthogonal.
C Output:
C   orbital sigma vectors in SCVECS(NCONF+1:NVAR) and
C   SOVECS(NCONF+1:NVAR)
C
C Scratch:
C   UDV; dimension: NASHT*NASHT
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION SCVECS(NVAR,*), SOVECS(*)
      DIMENSION DTV(*), FC(*), FTV(NORBT,NORBT,*),
     &          FTQ(NORBT,NASHDI,*)
      DIMENSION DV(*), FXC(NORBT,NORBT,*), FXV(NORBT,NORBT,*),
     &          FXQ(NORBT,NASHDI,*)
      DIMENSION UBO(NORBT,NORBT,*),UGORB(NORBT,NORBT), UDV(NASHDI,*)
C
      PARAMETER (HALF = 0.5D0)
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0, D8 = 8.0D0)
C
C Used from common blocks:
C   INFVAR : NCONF,NWOPT,NVAR,JWOP(2,*)
C   INFORB : NORBT,NNASHX,NNORBT,?
C   INFIND : IROW(*),ISMO(*)
C   INFDIM : NASHDI
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>

C
C ***** Step 1: calculate sigma vectors of configuration trial vectors
C
      IF (NCSIM.LE.0) GO TO 1000
C
C ***
C
      DO 900 ICSIM = 1,NCSIM
C
C **  Unpack DTV into UDV
C
         JDTV = 1 + (ICSIM-1)*NNASHX
         CALL DSPTSI(NASHT,DTV(JDTV),UDV)
C
         JSYMK = 0
         JSYML = 0
         DO 300 IG = 1,NWOPT
            IG1 = NCONF + IG
            K = JWOP(1,IG)
            L = JWOP(2,IG)
            ISYMK = ISMO(K)
            ISYML = ISMO(L)
            IF (ISYMK.NE.JSYMK) THEN
               JSYMK = ISYMK
               NORBK = NORB(ISYMK)
               IIORBK = IIORB(ISYMK)
               IORBK = IORB(ISYMK)
               NISHK = NISH(ISYMK)
               NASHK = NASH(ISYMK)
               IASHK = IASH(ISYMK)
            END IF
            IF (ISYML.NE.JSYML) THEN
               JSYML = ISYML
               NORBL = NORB(ISYML)
               IIORBL = IIORB(ISYML)
               IORBL = IORB(ISYML)
               NISHL = NISH(ISYML)
               NASHL = NASH(ISYML)
               IASHL = IASH(ISYML)
            END IF
            NK = K - IORBK
            NL = L - IORBL
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF (ITYPK .EQ. JTINAC) THEN
C **        first index inactive:
              SCVECS(IG1,ICSIM) = SCVECS(IG1,ICSIM) + D4*FTV(L,K,ICSIM)
C    *          + D8*BOVLP(ICSIM)*FC(LK)
C             HJ-850107-we assume BOVLP(ICSIM) = <BCSIM|CREF> = 0
            ELSE
C **        first index active:
              NKW = ISW(K) - NISHT
              TEMP = D0
              NX = NISHL
              DO 100 NXW = IASHL+1,IASHL+NASHL
                NX = NX + 1
                IF (NX.LE.NL) THEN
                  IFCXL = IIORBL+IROW(NL)+NX
                ELSE
                  IFCXL = IIORBL+IROW(NX)+NL
                END IF
                TEMP = TEMP + UDV(NXW,NKW)*FC(IFCXL)
  100         CONTINUE
              SCVECS(IG1,ICSIM) = SCVECS(IG1,ICSIM) +
     *          D2*(TEMP + FTQ(L,NKW,ICSIM))
C
            END IF
            IF (ITYPL .EQ. JTACT) THEN
C **        second index active:
               NLW = ISW(L) - NISHT
               TEMP = D0
               NX = NISHK
               DO 200 NXW = IASHK+1,IASHK+NASHK
                 NX = NX + 1
                 IF (NX.LE.NK) THEN
                   IFCXK = IIORBK+IROW(NK)+NX
                 ELSE
                   IFCXK = IIORBK+IROW(NX)+NK
                 END IF
                 TEMP = TEMP + UDV(NXW,NLW)*FC(IFCXK)
  200          CONTINUE
               SCVECS(IG1,ICSIM) = SCVECS(IG1,ICSIM) -
     *           D2*(TEMP + FTQ(K,NLW,ICSIM))
            END IF
C ***    next IG
  300    CONTINUE
C
  900 CONTINUE
C
C ***** Step 2: calculate sigma vectors of orbital trial vectors
C
 1000 IF (NOSIM.LE.0) GO TO 2000
C
C **  Unpack DV into UDV
C
      CALL DSPTSI(NASHT,DV,UDV)
C
      IF (NCSIM.GE.0) THEN
         ISO = NCONF
      ELSE
         ISO = 0
      END IF
      DO 1900 IOSIM = 1,NOSIM
         CALL DGETRN(UBO(1,1,IOSIM),NORBT,NORBT)
C        ... transpose to get UTBO
         JSYMK = 0
         JSYML = 0
         DO 1300 IG = 1,NWOPT
            K = JWOP(1,IG)
            L = JWOP(2,IG)
            ISYMK = ISMO(K)
            ISYML = ISMO(L)
            IF (ISYMK.NE.JSYMK) THEN
               JSYMK = ISYMK
               NORBK = NORB(ISYMK)
               IORBK = IORB(ISYMK)
               NISHK = NISH(ISYMK)
               NASHK = NASH(ISYMK)
               IASHK = IASH(ISYMK)
            END IF
            IF (ISYML.NE.JSYML) THEN
               JSYML = ISYML
               NORBL = NORB(ISYML)
               IORBL = IORB(ISYML)
               NISHL = NISH(ISYML)
               NASHL = NASH(ISYML)
               IASHL = IASH(ISYML)
            END IF
C
C **        Multiply together unpacked kappa (UBO) and gradient (UGORB)
C           matrices (here UBO(r,s) = borb(sr) and GORB(r,s) = GORB(rs))
C
            SOVECS(ISO + IG) = SOVECS(ISO + IG) - HALF*
     *        ( DDOT(NORBL,UBO(IORBL+1,K,IOSIM),1,UGORB(IORBL+1,L),1)
     *        - DDOT(NORBK,UBO(IORBK+1,L,IOSIM),1,UGORB(IORBK+1,K),1) )
C
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF (ITYPK .EQ. JTINAC) THEN
C **        first index inactive:
               SOVECS(ISO + IG) = SOVECS(ISO + IG) +
     *            D4*(FXC(L,K,IOSIM) + FXV(L,K,IOSIM))
            ELSE
C **        first index active:
              NKW = ISW(K) - NISHT
              TEMP = D0
              DO 1100 NX = 1,NASHK
                NXW = IASHK + NX
                IX  = IORBK + NISHK + NX
                TEMP = TEMP + UDV(NXW,NKW)*FXC(IX,L,IOSIM)
 1100         CONTINUE
              SOVECS(ISO + IG) = SOVECS(ISO + IG) +
     *          D2*(TEMP + FXQ(L,NKW,IOSIM))
C
            END IF
            IF (ITYPL .EQ. JTACT) THEN
C **        second index active:
               NLW = ISW(L) - NISHT
               TEMP = D0
               DO 1200 NX = 1,NASHL
                 NXW = IASHL + NX
                 IX  = IORBL + NISHL + NX
                 TEMP = TEMP + UDV(NXW,NLW)*FXC(IX,K,IOSIM)
 1200          CONTINUE
               SOVECS(ISO + IG) = SOVECS(ISO + IG) -
     *           D2*(TEMP + FXQ(K,NLW,IOSIM))
            END IF
C ***    next IG
 1300    CONTINUE
C
         IF (NCSIM.GE.0) THEN
            ISO = ISO + NVAR
         ELSE
            ISO = ISO + NWOPT
         END IF
         CALL DGETRN(UBO(1,1,IOSIM),NORBT,NORBT)
C        ... transpose to get UBO back
 1900 CONTINUE
C
 2000 CONTINUE
      RETURN
C
C ***** End of subroutine ORBSIG
C
      END
