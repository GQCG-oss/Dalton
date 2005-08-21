C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
C$Id: sirsolpcm.F,v 1.4 2005-08-21 11:43:12 klopper Exp $
C
#include <single.h>
C  /* Deck pcmgrd */
      SUBROUTINE PCMGRD(CREF,CMO,INDXCI,DV,G,ESOLT,WRK,LFREE)
C
C     Written by R. Cammi, L. Frediani, B. Mennucci, and K.Ruud,
C                Pisa June 2001, based on SOLGRD by hjj.
C
C   Purpose:  calculate MCSCF energy and gradient contribution
C             from a surrounding medium
C
C   Output:
C    G          MCSCF gradient with solvation contribution added
C    ESOLT      total solvation energy
C
#include <implicit.h>
#include <priunit.h>
#include <pi.h>
C
      DIMENSION CREF(*), CMO(*), INDXCI(*)
      DIMENSION DV(*),   G(*),   WRK(*)
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0,
     &            DCVAL = D2, FPI = 4.0D0 * PI )
#include <thrzer.h>
#include <dummy.h>
#include <iratdef.h>
#include <mxcent.h>
C
C Used from common blocks:
C   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFIND: IROW(*)
C   INFTAP: LUSOL,  LUIT2
C   INFPRI: LUERR,  LUPRI,   IPRSOL
C
#include <maxash.h>
#include <maxorb.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <nuclei.h>
#include <orgcom.h>
#include <pcmnuclei.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infinp.h>
#include <inftap.h>
#include <infpri.h>
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
C
      LOGICAL     FIRST, FNDLAB, EXP1VL, TOFILE, TRIMAT
      CHARACTER*8 STAR8, SOLVDI, EODATA, LABINT(9*MXCENT)
      SAVE        FIRST
      DATA        FIRST/.TRUE./, STAR8/'********'/
      DATA        SOLVDI/'SOLVDIAG'/, EODATA/'EODATA  '/
C
C     Statement functions;
C     define automatic arrays (dynamic core allocation)
C
C
      CALL QENTER('PCMGRD')
C
C
C     Core allocation
C        DIASH  diagonal of solvent contribution to Hessian
C        GRDLM  TELM gradient for current l,m value in the l,m loop
C        UCMO   CMO unpacked (i.e. no symmetry blocking)
C        RLMAC  active-active subblock of RLM
C        RLM    R(l,m) integrals for current l,m value in l,m loop
C
C        We may able to reuse J1 for storing J2
C
      KJENAO = 1
      KJEN   = KJENAO + NNBASX
      KJ1AO  = KJEN   + NNORBX
      KJ1    = KJ1AO  + NNBASX*NSYM
      KJENAC = KJ1    + NNORBX
      KJ1AC  = KJENAC + NNASHX
      KDIASH  = KJ1AC  + NNASHX
      KGRDLM = KDIASH + NVAR
      KUCMO  = KGRDLM + NVARH
      KJ2GRD = KUCMO  + NORBT*NBAST
      KPOT   = KJ2GRD + NVAR
      KDENC  = KPOT   + NTS
      KDENV  = KDENC  + NNBASX
      KW10   = KDENV  + NNBASX
C     1.3 rest of CSF contribution
      KDIALM = KW10
      KW20   = KDIALM + NVAR
C
C     Allocations for non-equlibrium energy solvation
      IF (NONEQ) THEN
         KMPOT  = KW20
         KQSEGR = KMPOT  + NTS * NTS 
         KPOTGR = KQSEGR + NTS
         KW21   = KPOTGR + NTS
      ELSE
         KW21   = KW20
      END IF
C
      LW21   = LFREE - KW21
C
      KWRK1  = KW21
      LWRK1  = LFREE  - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('PCMGRD',-KWRK1,LFREE)
C     
      IF (IPRPCM .GE. 130) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *        ' >>> PCMGRD - gtot (input) - non-zero elements',
     *        '     NCONF, NWOPT =',NCONF,NWOPT
         DO 40 I = 1,NCONF
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' conf #',I,G(I)
 40      CONTINUE
         DO 50 I = NCONF+1,NVAR
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *           ' orb  #',I,G(I)
 50      CONTINUE
      END IF
C
C
C     Read and check dimension information (if first read) and
C     nuclear contributions to ERLM (always).
C
      IF (IPRPCM .GE. 20 .AND. NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' >>> PCMGRD - DV matrix :'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF
C
C     Unpack symmetry blocked CMO
C     Loop over l,m expansion
C
      CALL UPKCMO(CMO,WRK(KUCMO))
      CALL DZERO(WRK(KDIASH),NVAR)
C
C     Read JEN = J(en) + J(ne) in ao basis and transform to mo basis.
C     Extract active-active block in WRK(KJENAC).
C     
      LUPBKP = LUPROP
      IF (LUPROP .LT. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      IF (.NOT. (FNDLAB('NE-PCMIN',LUPROP))) THEN
         CALL PCMJMAT(WRK(KJENAO),NNBASX,WRK(KWRK1),LWRK1)
      END IF
      REWIND (LUPROP)
      CALL REAPCM('NE-PCMIN','PCMGRD  ',LUPROP,WRK(KJENAO),NNBASX)
      CALL UTHU(WRK(KJENAO),WRK(KJEN),WRK(KUCMO),WRK(KWRK1),
     &          NBAST,NORBT)
      IF (NASHT .GT. 0) THEN
         CALL GETAC2(WRK(KJEN),WRK(KJENAC))
      END IF
      IF (IPRPCM .GE. 15) THEN
         WRITE (LUPRI,'(/A)') ' JEN_ao matrix:'
         CALL OUTPAK(WRK(KJENAO),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' JEN_mo matrix:'
         CALL OUTPAK(WRK(KJEN),  NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
            WRITE (LUPRI,'(/A)') ' JEN_ac matrix:'
            CALL OUTPAK(WRK(KJENAC),NASHT,1,LUPRI)
         END IF
      END IF
C     
C     Expextation value of JEN (=PB)
C     
      TJEN     = SOLELM(DV,WRK(KJENAC),WRK(KJEN),TJENAC)
      IF (IPRPCM .GE. 6) THEN
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- JEN expectation value(=PB) :',TJEN
         WRITE (LUPRI,'(A,F17.8)')
     *   ' --- active part of JEN(=PB)    :',TJENAC
      END IF
Cbm   PB=-TJEN
Cbm   WRITE(LUPRI,*)'PB =',PB
C     
C     Read J2 in ao basis and transform to mo basis.
C     Extract active-active block in WRK(KJ2AC).
C     
      CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDENC),WRK(KDENV),
     &            CMO,DV,WRK(KWRK1),LW21)
      CALL DAXPY(NNBASX,1.0D0,WRK(KDENV),1,WRK(KDENC),1)
      EXP1VL = .TRUE.
      TOFILE = .FALSE.
      CALL J1INT(WRK(KPOT),EXP1VL,WRK(KDENC),1,TOFILE,'NPETES ',
     &           1,WRK(KWRK1),LW21)
      CALL V2Q(WRK(KWRK1),LW21,WRK(KPOT),QSE,QET,.FALSE.)
      CALL DSCAL(NTS,-1.0D0,QSE,1)
C     
C     Read J1 (=ELECTROSTATIC POTENTIAL) in ao basis and transform to mo
C     basis.  Extract active-active block in WRK(KJ1AC).
C     
      XI = DIPORG(1)
      YI = DIPORG(2)
      ZI = DIPORG(3)
Ckr
Ckr should be parallelized, but a lot of information to send
Ckr However, in general not used for SCF and DFT optimizations
Ckr
      DO I = 1 , NTSIRR
         L = 1
         NCOMP = NSYM
         DIPORG(1) = XCTS(I)
         DIPORG(2) = YCTS(I)
         DIPORG(3) = ZCTS(I)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WRK(KJ1AO),'NPETES ',NCOMP,WRK(KWRK1),LW21,LABINT,
     &               INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,DUMMY,EXP1VL,
     &               DUMMY,IPRRSP)
         CALL UTHU(WRK(KJ1AO),WRK(KJ1),WRK(KUCMO),WRK(KWRK1),
     &        NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KJ1),WRK(KJ1AC))
         END IF
         IF (NONEQ) THEN
            WRK(KPOT + I - 1) = SOLELM(DV,WRK(KJ1AC),WRK(KJ1),TJ1AC)
         END IF
         IF (IPRPCM .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' J1_ao matrix:'
            CALL OUTPAK(WRK(KJ1AO),NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' J1_mo matrix:'
            CALL OUTPAK(WRK(KJ1),  NORBT,1,LUPRI)
            IF (NASHT .GT. 0) THEN
               WRITE (LUPRI,'(/A)') ' J1_ac matrix:'
               CALL OUTPAK(WRK(KJ1AC),NASHT,1,LUPRI)
            END IF
         END IF
         CALL DAXPY(NNORBX,-QSE(I),WRK(KJ1),1,WRK(KJEN),1)
C
C     Due to the computational cost of building the CI gradient, and the lack
C     of importance on convergence rates, we skip the construction of the 
C     solvent contribution to the diagonal Hessian. K.Ruud, Oct.-01
C
C         CALL DZERO(WRK(KGRDLM),NVARH)
C         IF (NCONF .GT. 1) THEN
C            CALL SOLGC(CREF,WRK(KJ1AC),TJ1AC,WRK(KGRDLM),INDXCI,
C     &                 WRK(KWRK1),LWRK1)
CC           CALL SOLGC(CREF,RLMAC,TELMAC,GLMCI,INDXCI,WRK,LWRK)
CC         END IF
C         IF (NWOPT .GT. 0) THEN
C            CALL SOLGO(DCVAL,DV,WRK(KJ1),WRK(KGRDLM+NCONF))
C         END IF
C         READ(LUGRDQ)WRK(KJ2GRD)
C         DO J =  NCONF, NVAR - 1 
C            WRK(KDIASH+J) = WRK(KDIASH+J) - WRK(KGRDLM+J)*WRK(KJ2GRD+J)
C         ENDDO
      ENDDO
Ckr
Ckr   End of parallelization loop
Ckr
C
C     Expextation value of J + X(0) = PB + PX
C     
      IF (NASHT .GT. 0) THEN
         CALL GETAC2(WRK(KJEN),WRK(KJENAC))
      END IF
      TJEN     = SOLELM(DV,WRK(KJENAC),WRK(KJEN),TJENAC)
      IF (IPRPCM .GE. 6) THEN
         CALL OUTPAK(WRK(KJENAC),NASHT,1,LUPRI)
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- TJEN expectation value :',TJEN
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- active part of TJEN    :',TJENAC
      END IF
C     
      CALL DZERO(WRK(KGRDLM),NVARH)
      IF (NCONF .GT. 1) THEN
         CALL SOLGC(CREF,WRK(KJENAC),TJENAC,WRK(KGRDLM),INDXCI,
     &              WRK(KWRK1),LWRK1)
C        CALL SOLGC(CREF,RLMAC,TELMAC,GLMCI,INDXCI,WRK,LWRK)
      END IF
      IF (NWOPT .GT. 0) THEN
         CALL SOLGO(DCVAL,DV,WRK(KJEN),WRK(KGRDLM+NCONF))
      END IF
C     
C     
C     Obtain DIALM = diagonal TE(l,m) Hessian
C                  = 2 ( <i|R(l,m)|i> - TE(l,m) )
C     Add the DIALM contribution and the GRDLM contribution
C     to solvent Hessian diagonal.
C     
      CALL SOLDIA(TJENAC,WRK(KJENAC),INDXCI,
     *            WRK(KJEN),DV,WRK(KDIALM),WRK(KW21),LW21)
C     CALL SOLDIA(TELM,RLMAC,INDXCI,RLM,DV,DIAG,WRK,LFREE)
C     
      DO 420 I = 0,(NVAR-1)
         WRK(KDIASH+I) = WRK(KDIASH+I)
     *                 - WRK(KDIALM+I)
  420 CONTINUE
C     
C     test orthogonality
C     
      IF (IPRPCM .GE. 120) THEN
         WRITE (LUPRI,'(/A)')' >>> PCMGRD - grdj1, grdj2, dialm, '//
     &                      'diash, cref'
         DO 430 I = 1,NCONF
            WRITE (LUPRI,'(A,I10,3F10.6)') ' conf #',I,
     *            WRK(KDIALM-1+I),
     *            WRK(KDIASH-1+I),CREF(I)
  430    CONTINUE
      END IF
C     
      TEST = DDOT(NCONF,CREF,1,WRK(KGRDLM),1)
      IF (ABS(TEST) .GT. 1.D-8) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A,/A,1P,D12.4)')
     *   ' >>> PCMGRD WARNING, for B',
     *   ' <CREF | GRADB > =',TEST
      END IF
C     
C      TEST = DDOT(NCONF,CREF,1,WRK(KGRDJ1),1)
      IF (ABS(TEST) .GT. 1.D-8) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A,/A,1P,D12.4)')
     *   ' >>> PCMGRD WARNING, for J1 ',
     *   ' <CREF | GRADJ1 > =',TEST
      END IF
C      TEST = DDOT(NCONF,CREF,1,WRK(KGRDJ2),1)
      IF (ABS(TEST) .GT. 1.D-8) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A,/A,1P,D12.4)')
     *   ' >>> PCMGRD WARNING, for J2',
     *   ' <CREF | GRADJ2 > =',TEST
      END IF
C     
C     Add PCM gradient contribution to MCSCF gradient
C     
      CALL DAXPY(NVARH,-D1,WRK(KGRDLM),1,G,1)
      IF (IPRPCM .GE. 140) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *      ' >>> PCMGRD - grdB, gtot (accum) - non-zero grdlm',
     *      '     NCONF, NWOPT =',NCONF,NWOPT
         DO 440 I = 1,NCONF
            IF (WRK(KGRDLM-1+I) .NE. D0)
     *         WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' conf #',I,WRK(KGRDLM-1+I),G(I)
  440    CONTINUE
         DO 450 I = NCONF+1,NVAR
            IF (WRK(KGRDLM-1+I) .NE. D0)
     *         WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' orb  #',I,WRK(KGRDLM-1+I),G(I)
  450    CONTINUE
      END IF
C     
C     
C     
      IF (IPRPCM .GE. 130) THEN
         WRITE (LUPRI,'(/A/A,2I10)')
     *      ' >>> PCMGRD - gtot (output) - non-zero elements',
     *      '     NCONF, NWOPT =',NCONF,NWOPT
         DO 840 I = 1,NCONF
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *      ' conf #',I,G(I)
  840    CONTINUE
         DO 850 I = NCONF+1,NVAR
            IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *      ' orb  #',I,G(I)
  850    CONTINUE
      END IF
C
C     Calculate PCM solvent energy contributions 
C
        PC = D0
	DO ITS = 1, NTSIRR
        DO JATOM = 1, NUCDEP
          R2 = (PCMCORD(1,JATOM)-XCTS(ITS))**2 + (PCMCORD(2,JATOM)-
     *         YCTS(ITS))**2 + (PCMCORD(3,JATOM)-ZCTS(ITS))**2
          R = DSQRT(R2)
          PC = PC + NSYM*QSE(ITS) * PCMCHG(JATOM) / R
        ENDDO
      ENDDO
      ESOLT = - DP5 * TJEN + DP5 * PC + UNZ
      IF (IPRPCM .GE. 6) THEN
         WRITE(LUPRI,*) 'ESOLT, TJEN, PC, UNZ'
         WRITE(LUPRI,'(3f20.15)') TJEN, PC, UNZ
      END IF
C
C     We have used "effective" nuclear charges in non-equilibrium
C     calculations, and we need to modify the solvent energy to account
C     for this (related to "slow" charges).
C
Clf If we are in a nonequilibrium case we correct the solvation energy.
C
      IF (NONEQ) THEN
         CALL NEQENR(ESOLT,WRK(KMPOT),WRK(KQSEGR),WRK(KPOTGR),WRK(KPOT),
     $        WRK(KWRK1),LWRK1)
      END IF
C
C     Write solvent contribution to Hessian diagonal on LUIT2
C
      IF (LUIT2 .GT. 0) THEN
         NC4 = MAX(NCONF,4)
         NW4 = MAX(NWOPT,4)
         REWIND LUIT2
         IF (FNDLAB(EODATA,LUIT2)) BACKSPACE LUIT2
         WRITE (LUIT2) STAR8,STAR8,STAR8,SOLVDI
         IF (NCONF .GT. 1) CALL WRITT(LUIT2,NC4,WRK(KDIASH))
         IF (NWOPT .GT. 0) CALL WRITT(LUIT2,NW4,WRK(KDIASH+NCONF))
         WRITE (LUIT2) STAR8,STAR8,STAR8,EODATA
      END IF
C
C
      FIRST = .FALSE.
      IF (LUPBKP .LT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('PCMGRD')
      RETURN
C     end of pcmgrd.
      END
C  /* Deck pcmfck */
      SUBROUTINE PCMFCK(DC,DV,FSOL,CMO,MOBAS,ESOLT,WRK,LFREE,IPRINT)
C
C     Written by R. Cammi, L. Frediani, B. Mennucci, and K.Ruud,
C                Pisa June 2001, based on SOLGRD by hjj.
C
C   Purpose:  calculate effective solvent Fock operator
C
C   Output:
C    FSOL     MCSCF gradient with solvation contribution added
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <thrzer.h>
#include <iratdef.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <codata.h>
C
C Used from common blocks:
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFPRI: LUERR,  LUPRI
C
#include <pcm.h>
#include <pcmlog.h>
#include <nuclei.h>
#include <pcmnuclei.h>
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
#include <maxaqn.h>
#include <symmet.h>
C
#include <orgcom.h>
C
      PARAMETER (ANTOAU=1.0D+00/XTANG)
C
      LOGICAL FNDLAB, EXP1VL, TOFILE
      DIMENSION DV(NNBAST), DC(NNBAST)
      DIMENSION FSOL(*), CMO(*), WRK(LFREE)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0,  D2 = 2.0D0, DP5 = 0.5D0 )
C
      LOGICAL     MOBAS
      LOGICAL     FIRST
      SAVE        FIRST
      DATA        FIRST/.TRUE./
#include <sdpre.h>
C
C     Statement functions;
C     define automatic arrays (dynamic core allocation)
C
      CALL QENTER('PCMFCK')
C
C
C     Core allocation
C        UCMO   CMO unpacked (i.e. no symmetry blocking)
C        RLM    R(l,m) integrals for current l,m value in l,m loop
C
C
      KJENAO = 1
      KJ1AO  = KJENAO + NNBASX
      KJAO   = KJ1AO  + NNBASX
      KTJ2   = KJAO   + NNBASX
      KPKDEN = KTJ2   + NTS
      KWRK1  = KPKDEN + NNBASX
      LWRK1  = LFREE  - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('PCMFCK',-KWRK1,LFREE)
C
      ESOLT = D0
      QET = D0
C
Clf   I don't think the next call is necessary, anyway it cannot take a long time..
C
      CALL DZERO(WRK(KTJ2),NTS)
C
C        Read JEN in ao basis 
C
      IF (LUPROP .LT. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',
     &     'SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      IF (.NOT. (FNDLAB('NE-PCMIN',LUPROP))) THEN
         CALL PCMJMAT(WRK(KJENAO),NNBASX,WRK(KWRK1),LWRK1)
      END IF
      REWIND (LUPROP)
      CALL REAPCM('NE-PCMIN','PCMFCK  ',LUPROP,WRK(KJENAO),NNBASX)
      IF (IPRINT .GE. 15) THEN
         WRITE (LUPRI,'(/A)') ' JNE_ao matrix:'
         CALL OUTPAK(WRK(KJENAO),NBAST,1,LUPRI)
      END IF
C
      IF (.NOT. MOBAS) THEN
C
C     Loop over all tesserae. We temporarily use FSOL for the unpacked
C     inactive density matrix
C
         CALL PKSYM1(FSOL,DC,NBAS,NSYM,-1)
C
C benedetta: EN contribution
C    PB = interaction solute electrons - polarization nuclear charges
C             PB = sum_i [Vel(i) * QSN(i)]
C
         PB = -DDOT(NNBASX,FSOL,1,WRK(KJENAO),1)
         EXP1VL = .TRUE.
         TOFILE = .FALSE.
         CALL J1INT(WRK(KTJ2),EXP1VL,DC,1,TOFILE,'NPETES ',
     &              1,WRK(KWRK1),LWRK1)
         CALL V2Q(WRK(KWRK1),WRK(KTJ2),QSE,QET,.FALSE.)
      ELSE
         DO ITS = 1, NTSIRR
            QET = QET + QSE(ITS)
         END DO
         QET = QET * DBLE (MAXREP+1)
      END IF
C
C     Normalization of electronic apparent charges
C     For ICOMPCM=1 the correction on charge ITS is done proportionally to
C     the area of the tessera ITS
C     For ICOMPCM=2 the correction is done with a constant factor.
C
C for Ruud: I'm not sure that NRHFEL is the number of electrons, and where it
C           is stored.
      DNE = SDPRE(2*NISHT+NACTEL)
      IF (NONEQ) THEN
         TCH = DNE * (EPSINF - D1) / EPSINF
      ELSE
         TCH = DNE * (EPS - D1) / EPS
      END IF         
C
C     Option 1
C
      QETN = D0
      FE = D1
Clf dirty trick
Clf      ICOMPCM = 0
CLF
      IF(ICOMPCM.EQ.1) THEN
         SUPTOT = STOT * ANTOAU * ANTOAU
         DIFF = TCH - QET
         DO ITS = 1, NTSIRR
            QSE(ITS) = QSE(ITS) + DIFF * AS(ITS) / SUPTOT
            QETN = QETN + QSE(ITS)
         ENDDO
C
C     Option 2
C
      ELSE IF(ICOMPCM.EQ.2) THEN
         FE = TCH / QET
         DO ITS = 1, NTSIRR
            QSE(ITS) = QSE(ITS) * FE
            QETN = QETN + QSE(ITS)
         ENDDO
      ENDIF
      QETN = QETN * DBLE(MAXREP + 1)
C
C     Print total induced charge before and after renormalization
C
      IF(IPRPCM.GE.5)THEN
         IF(ICOMPCM.NE.0) WRITE(LUPRI,1000) QET, TCH, QETN
         IF(ICOMPCM.EQ.0) WRITE(LUPRI,1100) QET, TCH
      END IF
C      
 1000 FORMAT(1X,'ELECTRONIC APPARENT CHARGE',F10.5,' THEORETICAL ',
     * F10.5,' RENORMALIZED',F10.5)
 1100 FORMAT(1X,'ELECTRONIC APPARENT CHARGE',F10.5,' THEORETICAL ',
     * F10.5,' NOT RENORMALIZED')
C
      EXP1VL = .FALSE.
      TOFILE = .FALSE.
      CALL J1INT(QSE,EXP1VL,WRK(KJENAO),1,TOFILE,'NPETES ',
     &           1,WRK(KWRK1),LWRK1)
C
C
C Benedetta: NE and EE contributions
C
C    PC = interaction solute nuclei - polarization electronic charges
C             PC = sum_i [Vnuc(i) * QSE(i)]
C
C    PX = interaction solute electrons - polarization electronic charges
C             PX = sum_i [Vel(i) * QSE(i)]
C
        PC = D0
	DO ITS = 1, NTSIRR
        DO JATOM = 1, NUCDEP
          R2 = (PCMCORD(1,JATOM)-XCTS(ITS))**2 + (PCMCORD(2,JATOM)-
     *         YCTS(ITS))**2 + (PCMCORD(3,JATOM)-ZCTS(ITS))**2
          R = DSQRT(R2)
          PC = PC + NSYM * QSE(ITS) * PCMCHG(JATOM) / R
        ENDDO
      ENDDO
      PX = DDOT(NNBASX,FSOL,1,WRK(KJENAO),1) - PB
C
C     Electrostatic contribution:
C     GES = 1/2 (PB + PC) + 1/2 PX + UNZ
C     UNZ = 1/2 interaction nuclei - polarization nuclear charges
      ESOLT = - DP5*(PX+PB) + DP5*PC + UNZ
      print *, 'px,pb,pc,unz',px,pb,pc,unz,esolt
C     
C     Symmetry pack TAO using WRK(KJAO)
C     Transform to MO basis
C
      CALL PKSYM1(WRK(KJENAO),WRK(KJAO),NBAS,NSYM,+1)
      IF (IPRINT .GE. 15) THEN
         WRITE (LUPRI,'(/A)')
     &      ' PCMFCK effective solvent B_ao matrix (not packed):'
         CALL OUTPAK(WRK(KJENAO),NBAST,1,LUPRI)
      END IF
      CALL DZERO(FSOL,NNORBT)
      IF (MOBAS) THEN 
         DO 700 ISYM = 1,NSYM
            JCMO = ICMO(ISYM) + 1
            J1AO = KJAO + IIBAS(ISYM)
            J1MO = 1 + IIORB(ISYM)
            CALL UTHU(WRK(J1AO),FSOL(J1MO),CMO(JCMO),
     &                WRK(KWRK1),NBAS(ISYM),NORB(ISYM))
  700    CONTINUE
       ELSE 
          CALL DCOPY(NNBASX,WRK(KJAO),1,FSOL(1),1)
       END IF
C
C
       IF (IPRINT .GE. 6 .AND. MOBAS) THEN
         WRITE (LUPRI,'(/A)')
     &      ' PCMFCK: effective solvent one-electron matrix (MO basis):'
         CALL OUTPKB(FSOL,NORB,NSYM,1,LUPRI)
      END IF
      FIRST = .FALSE.
      CALL QEXIT('PCMFCK')
      RETURN
C     end of pcmfck.
      END

C  /* Deck pcmlin */
      SUBROUTINE PCMLIN(NCSIM,NOSIM,BCVECS,BOVECS,CREF,CMO,INDXCI,
     &                  DV,DTV,SCVECS,SOVECS,ORBLIN,VTEX,WRK,LWRK)
C
C Copyright 25-Apr-1990 Hans Joergen Aa. Jensen
C Common driver for PCMLNC and PCMLNO
C
#include <implicit.h>
      DIMENSION BCVECS(1),BOVECS(1),CREF(1),CMO(1),INDXCI(1)
      DIMENSION DV(1),DTV(1),SCVECS(1),SOVECS(1),WRK(LWRK)
      DIMENSION VTEX(*)
      LOGICAL   ORBLIN
C
C Used from common blocks:
C   INFLIN : NWOPPT,NVARPT
C
#include <maxorb.h>
#include <inflin.h>
C
      CALL QENTER('PCMLIN')
C
      IF (NCSIM .GT. 0) THEN
         CALL PCMLNC(NCSIM,BCVECS,CREF,CMO,INDXCI,
     *               DV,DTV,SCVECS,VTEX,WRK,LWRK)
C        CALL PCMLNC(NCSIM,BCVEC,CREF,CMO,INDXCI,
C    *               DV,DTV,SVEC,VTEX,WRK,LFREE)
C
      END IF
      IF ( NOSIM .GT.0 ) THEN
         IF ( .NOT.ORBLIN ) THEN
            NSVAR  = NVARPT
         ELSE
            NSVAR  = NWOPPT
         END IF
         CALL PCMLNO(NOSIM,BOVECS,CREF,CMO,INDXCI,
     *               DV, SOVECS,NSVAR,VTEX,WRK,LWRK)
C        CALL PCMLNO(NOSIM,BOVEC,CREF,CMO,INDXCI,
C    *               DV, SVEC,NSVEC,VTEX,,WRK,LFREE)
      END IF
      CALL QEXIT('PCMLIN')
      RETURN
      END
C  /* Deck pcmlnc */
      SUBROUTINE PCMLNC(NCSIM,BCVEC,CREF,CMO,INDXCI,
     *                  DV,DTV,SVEC,VTEX,WRK,LFREE)
C
C  Copyright 24-Jan-1987 Hans Jorgen Aa. Jensen
C  Revised for detci 17-Jul-1990 hjaaj
C
C  Purpose:  Calculate MCSCF Hessian contribution from a
C            surrounding medium to a csf trial vector.
C            Cavity radius = Rsol and dielectric constant = EPsol.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <dummy.h>
      DIMENSION BCVEC(*),  CREF(*), CMO(*) 
      DIMENSION INDXCI(*), DV(*),   DTV(NNASHX,*)
      DIMENSION SVEC(NVARPT,*),     WRK(*)
C
#include <iratdef.h>
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 )
#include <thrzer.h>
C
C
C  Used from common blocks:
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : NWOPH
C    INFLIN : NCONST, NVARPT, NWOPPT
C    INFPRI : LUERR
C
#include <maxash.h>
#include <maxorb.h>
#include <pcmdef.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inflin.h>
#include <inftap.h>
#include <infpri.h>
#include <pcm.h>
#include <pcmlog.h>
      DIMENSION VTEX(NTS,NCSIM)
C
      CALL QENTER('PCMLNC')
C
C     Core allocation
C
      KUCMO  = 1
      KJPX   = KUCMO  + NORBT*NBAST
      KJPXAC = KJPX   + NNORBX
      KPCMXB = KJPXAC + NNASHX
      KXBAC  = KPCMXB + NCSIM*NNORBX
      KEXBAC = KXBAC  + NCSIM*NNASHX
      KW10   = KEXBAC + NCSIM
C
      KCHRG   = KW10
      KW20    = KCHRG + NNBASX
C     2.1 read rlmao in ao basis and transform to rlm in mo basis
      KDV     = KDV
      KJPCMAO = KDV + NCSIM*NNBASX
      KW30    = KJPCMAO + NNBASX
C
      LW30    = LFREE   + 1 - KW30
C
      KNEED   = KW30
      IF (KNEED .GT. LFREE) CALL ERRWRK('PCMLNC',KNEED,LFREE)
C
      IF (IPRPCM .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> PCMLNC - svec(ci,1) on entry'
         DO I = 1,NCONST
            IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *      ' conf #',I,SVEC(I,1)
         END DO
         WRITE (LUPRI,'(/A)') ' >>> PCMLNC - svec(orb) on entry'
         CALL OUTPUT(SVEC(1+NCONST,1),1,NWOPPT,1,NCSIM,
     *               NSVEC,NCSIM,1,LUPRI)
      END IF                                                                
C
      CALL DZERO(WRK(KJPX),NNORBX)
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
      CALL DCOPY(NTS,QSN,1,WRK(KCHRG),1)
      CALL DAXPY(NTS,1.0D0,QSE,1,WRK(KCHRG),1)
      CALL J1INT(WRK(KCHRG),.FALSE.,WRK(KJPXAO),1,.FALSE.,'NPETES ',
     &           1,WRK(KW30),LW30)
      CALL UTHU(WRK(KJPXAO),WRK(KJPX),WRK(KUCMO),WRK(KW30),NBAST,NORBT)
      IF (NASHT .GT. 0) THEN
         CALL GETAC2(WRK(KJPX),WRK(KJPXAC))
      END IF
      IF (IPRPCM .GE. 15) THEN
         WRITE (LUPRI,'(/A)') ' J+X_mo matrix:'
         CALL OUTPAK(WRK(KJPX),  NORBT,1,LUPRI)
         IF (NASHT .GT. 0) THEN
            WRITE (LUPRI,'(/A)') ' J+X_ac matrix:'
            CALL OUTPAK(WRK(KJPXAC),NASHT,1,LUPRI)
         END IF
      END IF
C     
C     Expectation value of J+X
C     
      EXPJPX = SOLELM(DV,WRK(KJPXAC),WRK(KJPX),ACJPX)
      IF (IPRPCM .GE. 6) THEN
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- J+X expectation value :',EXPJPX
         WRITE (LUPRI,'(A,F17.8)')
     *        ' --- active part of J1X    :',ACJPX
      END IF
C
      DO ICSIM = 1, NCSIM
         CALL FCKDEN(.FALSE.,.TRUE.,DUMMY,WRK(KDV),CMO,DTV(1,ICSIM),
     &               WRK(KW30),LW30)
         CALL J1INT(WRK(KCHRG),.TRUE.,WRK(KDV),1,.FALSE.,'NPETES ',
     &              1,WRK(KW30),LW30)
         CALL DSCAL(NTS,-1.0D0,WRK(KCHRG),1)
         CALL V2Q(WRK(KW30),WRK(KCHRG),VTEX(1,ICSIM),QET,.FALSE.)
      END DO
C
      DO ICSIM = 1, NCSIM
         CALL DZERO(WRK(KJPCMAO),NNBASX)
         CALL J1INT(VTEX(1,ICSIM),.FALSE.,WRK(KJPCMAO),.FALSE.,1,
     &              'NPETES ',1,WRK(KW30),LW30)
         JPCMXB = KPCMXB + (ICSIM - 1)*NNORBX
         CALL UTHU(WRK(KJPCMAO),WRK(JPCMXB),WRK(KUCMO),WRK(KW30),
     &             NBAST,NORBT)
      END DO
C
      DO ICSIM = 1, NCSIM
         JPCMXB = KPCMXB + (ICSIM - 1)*NNORBX
         JXBAC  = KXBAC  + (ICSIM - 1)*NNASHX
         CALL GETAC2(WRK(JPCMXB),WRK(JXBAC))
         TEMP = SOLELM(DV,WRK(JXBAC),WRK(JPCMXB),XBAC)
         WRK(KEXBAC - 1 + ICSIM) = XBAC
      END DO
C
C     Make gradient-like contribution
C

      CALL SOLSC(NCSIM,0,BCVEC,CREF,SVEC,WRK(KXBAC),WRK(KJPXAC),
     *           WRK(KEXBAC),ACJPX,INDXCI,WRK(KW30),LW30)
C     CALL SOLSC(NCSIM,NOSIM,BCVECS,CREF,SVECS,
C    *           RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
C
C
C
C     ... orbital part of sigma vector(s)
C
      IF (IPRPCM .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> PCMLNC - svec(ci,1) on exit 2'
         DO I = 1,NCONST
            IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *      ' conf #',I,SVEC(I,1)
         END DO
         WRITE (LUPRI,'(/A)') ' >>> PCMLNC - svec(orb) on exit 2'
         CALL OUTPUT(SVEC(1+NCONST,1),1,NWOPPT,1,NCSIM,
     *               NSVEC,NCSIM,1,LUPRI)
      END IF                                                                
      IF (NWOPPT .GT. 0) THEN
         MWOPH  = NWOPH
         NWOPH  = NWOPPT
C        ... tell SOLGO only to use the NWOPPT first JWOP entries
         JSVECO = 1 + NCONST
         JPCMXB = KPCMXB
         DO 800 ICSIM = 1,NCSIM
            CALL SOLGO(D2,DV,WRK(JPCMXB),SVEC(JSVECO,ICSIM))
            IF (IPRSOL.GT.101) THEN
               WRITE(LUPRI,*)
     *         ' ORBITAL PART OF LIN TRANS CONF VEC NO:', ICSIM
               WRITE(LUPRI,*)' X(B,0) CONTRIBUTION'
               CALL OUTPUT(SVEC(JSVECO,ICSIM),1,NWOPPT,1,1,
     *                                        NWOPPT,1,1,LUPRI)
            END IF
            CALL SOLGO(D0,DTV(1,ICSIM),WRK(KJPX),SVEC(JSVECO,ICSIM))
            IF (IPRSOL.GT.101) THEN
               WRITE(LUPRI,*)
     *         ' ORBITAL PART OF LIN TRANS CONF VEC NO:', ICSIM
               WRITE(LUPRI,*)' + (j+X) CONTRIBUTION'
               CALL OUTPUT(SVEC(JSVECO,ICSIM),1,NWOPPT,1,1,
     *                                        NWOPPT,1,1,LUPRI)
            END IF
            JPCMXB = JPCMXB + NNORBX
  800    CONTINUE
         NWOPH  = MWOPH
      END IF
C
      IF (IPRPCM .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> PCMLNC - svec(ci,1) on exit'
         DO I = 1,NCONST
            IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *      ' conf #',I,SVEC(I,1)
         END DO
         WRITE (LUPRI,'(/A)') ' >>> PCMLNC - svec(orb) on exit'
         CALL OUTPUT(SVEC(1+NCONST,1),1,NWOPPT,1,NCSIM,
     *               NSVEC,NCSIM,1,LUPRI)
      END IF                                                                
      IF (LUPBKP .LT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('PCMLNC')
      RETURN
C     end of pcmlnc.
      END
C  /* Deck pcmlno */
      SUBROUTINE PCMLNO(NOSIM,BOVECS,CREF,CMO,INDXCI,
     *                  DV,SVEC,NSVEC,VTEX,WRK,LFREE)
C
C  Copyright 24-Jan-1987 Hans Jorgen Aa. Jensen
C
C  Purpose:  Calculate MCSCF Hessian contribution from a
C            surrounding medium to an orbital trial vector.
C            Cavity radius = Rsol and dielectric constant = EPsol.
C
C  NSVEC          may be NVAR or NWOPT, dependent on LINTRN
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      DIMENSION BOVECS(NWOPPT,*), CREF(*), CMO(*)
      DIMENSION INDXCI(*),        DV(*)
      DIMENSION SVEC(NSVEC,*),    WRK(*)
C
#include <iratdef.h>
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0, DP5 = 0.5D0 )
#include <thrzer.h>
C
C
C  Used from common blocks:
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : JWOP
C    INFLIN : NWOPPT, NVARPT, NCONST, NCONRF
C    INFPRI : LUERR,  IPRSOL
C
#include <pcmdef.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
      CHARACTER*9 LABINT(9*MXCENT)
      DIMENSION VTEX(NTS,NOSIM), QTEX(NTS,NOSIM), INTREP(9*MXCENT),
     &          INTADR(9*MXCENT)
#include <infinp.h>
#include <pcm.h>
#include <pcmlog.h>
#include "orgcom.h"
#include <inforb.h>
#include <infvar.h>
#include <inflin.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL FULHES, TOFILE, TRIMAT, EXP1VL
C
C
      CALL QENTER('PCMLNO')
C
C     Determine if full Hessian or only orbital Hessian
C
      FULHES = (NSVEC .EQ. NVARPT)
      IF (FULHES) THEN
         JSOVEC = 1 + NCONST
      ELSE
         JSOVEC = 1
      END IF
C
      IF (IPRSOL .GE. 40) THEN
         WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM PCMLNO >>>'
      END IF
      IF (IPRPCM .GE. 140) THEN
         IF (FULHES) THEN
            WRITE (LUPRI,'(/A)') ' >>> PCMLNO - svec(ci,1) on entry'
            DO 30 I = 1,NCONST
               IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *              ' conf #',I,SVEC(I,1)
 30         CONTINUE
         END IF
         WRITE (LUPRI,'(/A)') ' >>> PCMLNO - svec(orb) on entry'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *        NSVEC,NOSIM,1,LUPRI)
      END IF
C
C     Core allocation
C
      KUBO   = 1
      KW10   = KUBO   + NOSIM*N2ORBX
C
      KUCMO  = KW10   
      KPCMX  = KUCMO  + NORBT*NBAST
      KPCMXA = KPCMX  + NOSIM*NNORBX
      KJPXAC = KPCMXA + NOSIM*NNASHX
      KW20   = KJPXAC + NOSIM
C
      KJ1AO  = KW20
      KJ1SQ  = KJ1AO  + NNBASX*NSYM
      KJ1    = KJ1SQ  + N2ORBX
      KJ1X   = KJ1    + NNORBX
      KJ1XSQ = KJ1X   + NOSIM*NNORBX
      KJ1XAC = KJ1XSQ + N2ORBX
      KW30   = KJ1XAC  + NOSIM * NNASHX
C
C     3.0 SOLSC
      IF (FULHES) THEN
         KOVLP = KW30
         KW50   = KOVLP  + NOSIM
         LW50   = LFREE  + 1 - KW50
      ELSE
         KW50   = 0
         LW50   = LFREE + 1 - KW30
      END IF
C
      KNEED = MAX(KW20,KW30)
      IF (KNEED .GT. LFREE) CALL ERRWRK('PCMLNO',KNEED,LFREE)
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Calculate unpacked orbital trial vectors in UBO
C
C
      IF (NOSIM.GT.0) THEN
         DO 210 IOSIM = 1,NOSIM
            JUBO = KUBO + (IOSIM-1)*N2ORBX
            CALL UPKWOP(NWOPPT,JWOP,BOVECS(1,IOSIM),WRK(JUBO))
            IF (IPRSOL .GE. 55) THEN
               WRITE (LUPRI,2110) IOSIM,NOSIM
               CALL OUTPUT(WRK(JUBO),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                     LUPRI)
            END IF
  210    CONTINUE
      END IF
 2110 FORMAT (/,' Orbital trial vector unpacked to matrix form (no.',
     *        I3,' of',I3,')')
C
C     Contributions from all tesserae are included.
C
      CALL DZERO(WRK(KPCMX),NOSIM*NNORBX)
      XI = DIPORG(1)
      YI = DIPORG(2)
      ZI = DIPORG(3)
Ckr
Ckr   Parallelizable, but not top priority
Ckr
      DO I = 1, NTSIRR
         L = 1
         NCOMP = NSYM
         DIPORG(1) = XCTS(I)
         DIPORG(2) = YCTS(I)
         DIPORG(3) = ZCTS(I)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WRK(KJ1A0),'NPETES ',NCOMP,WRK(KW50),LW50,
     &               LABINT,INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,DUMMY,
     &               EXP1VL,DUMMY,IPRPCM)
         CALL UTHU(WRK(KJ1AO),WRK(KJ1),WRK(KUCMO),WRK(KW30),
     *        NBAST,NORBT)
         CALL DSPTSI(NORBT,WRK(KJ1),WRK(KJ1SQ))
         DO IOSIM = 1, NOSIM
            JUBO = KUBO + (IOSIM - 1) * N2ORBX
            JJ1X = KJ1X + (IOSIM - 1) * NNORBX
            JJ1XAC = KJ1XAC + (IOSIM - 1) * NNASHX
            CALL DZERO(WRK(KJ1XSQ),N2ORBX)
            CALL TR1UH1(WRK(JUBO),WRK(KJ1SQ),WRK(KJ1XSQ),1)
            CALL DGETSP(NORBT,WRK(KJ1XSQ),WRK(JJ1X))
            IF (NASHT .GT. 0) THEN
               CALL GETAC2(WRK(JJ1X),WRK(JJ1XAC))
            END IF
            IF (IPRPCM .GE. 15) THEN
               WRITE (LUPRI,'(/A)') ' J1X_mo matrix:'
               CALL OUTPAK(WRK(JJ1X),  NORBT,1,LUPRI)
               IF (NASHT .GT. 0) THEN
                  WRITE (LUPRI,'(/A)') ' J1X_ac matrix:'
                  CALL OUTPAK(WRK(JJ1XAC),NASHT,1,LUPRI)
               END IF
            END IF
C     
C     Expectation value of transformed potential on tesserae
C     
            VTEX(I,IOSIM)     = SOLELM(DV,WRK(JJ1XAC),WRK(JJ1X),TJ1XAC)
            IF (IPRPCM .GE. 6) THEN
               WRITE (LUPRI,'(A,F17.8)')
     *              ' --- J1X expectation value :',VTEX(I,IOSIM)
               WRITE (LUPRI,'(A,F17.8)')
     *              ' --- active part of J1X    :',TJ1XAC
            END IF
         ENDDO
         CALL DAXPY(NOSIM*NNORBX,QSN(I)+QSE(I),WRK(KJ1X),1,
     &              WRK(KPCMX),1)
      ENDDO
Ckr
Ckr   End parallelization
Ckr
      DIPORG(1) = XI
      DIPORG(2) = YI
      DIPORG(3) = ZI
C
      DO IOSUM = 1, NOSUM
         CALL V2Q(WRK(KW50),VTEX(1,IOSIM),QTEX(1,IOSIM),QTEXS,.FALSE.)
         IF (IPRPCM .GE. 6) THEN
            WRITE (LUPRI,'(A,F17.8)')
     *           ' --- J2X expectation value :',VTEX(I,IOSIM)
         END IF
      ENDDO
C
C     TRANSFORMED CHARGES MULTIPLIED TO POTENTIALS
C
Ckr
Ckr   Parallelizable
Ckr
      DO I = 1 , NTSIRR
         L = 1
         NCOMP = NSYM
         DIPORG(1) = XCTS(I)
         DIPORG(2) = YCTS(I)
         DIPORG(3) = ZCTS(I)
         EXP1VL    = .FALSE.
         TOFILE    = .FALSE.
         KPATOM    = 0
         TRIMAT    = .TRUE.
         CALL GET1IN(WRK(KJ1AO),'NPETES ',NCOMP,WRK(KW50),LW50,LABINT,
     &               INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,DUMMY,
     &               EXP1VL,DUMMY,IPRRSP)
         CALL UTHU(WRK(KJ1AO),WRK(KJ1),WRK(KUCMO),WRK(KW30),
     *        NBAST,NORBT)
         DO IOSIM = 1 , NOSIM
            JPCMX = KPCMX + (IOSIM - 1) * NNORBX 
            CALL DAXPY(NNORBX,-DP5*QTEX(I,IOSIM),WRK(KJ1),1,
     &                WRK(JPCMX),1)
         ENDDO
      ENDDO
Ckr
Ckr   End parallelizable loop
Ckr
C
C     Note: the 1/2 sum(t) [ gsol(ts)b(rt) - gsol(tr)b(st) ]
C           correction to the one-index transformation in Ryo
C           is added in lintrn, where g(rs) = gstd(rs) + gsol(rs).
C
      IF (LSYMRF .EQ. LSYMST) THEN
         NCOLIM = 1
      ELSE
         NCOLIM = 0
      END IF
      IF (FULHES .AND. NCONST .GT. NCOLIM) THEN
C
C        ... CSF part of sigma vectors
C
         DO 700 IOSIM = 1,NOSIM
            JPCMX  = KPCMX  + (IOSIM-1)*NNORBX
            JPCMXA = KPCMXA + (IOSIM-1)*NNASHX
            IF (LSYMRF .EQ. LSYMST) THEN
               WRK(KOVLP-1+IOSIM) = DDOT(NCONRF,CREF,1,SVEC(1,IOSIM),1)
            ELSE
               WRK(KOVLP-1+IOSIM) = D0
            END IF
            CALL GETAC2(WRK(JPCMX),WRK(JPCMXA))
            TJPX   = SOLELM(DV,WRK(JPCMXA),WRK(JPCMX),TJPXAC)
            IF (IPRPCM .GE. 40) WRITE (LUPRI,'(/A,I5,3F15.10)')
     *         ' IOSIM, C_OVLP, TJPXAC, TJPX :',
     *         IOSIM,WRK(KOVLP-1+IOSIM),TJPXAC,TJPX
            WRK(KJPXAC-1+IOSIM) = TJPXAC
  700    CONTINUE
         CALL SOLSC(0,NOSIM,DUMMY,CREF,SVEC,WRK(KPCMXA),DUMMY,
     *              WRK(KJPXAC),DUMMY,INDXCI,WRK(KW50),LW50)
C        CALL SOLSC(NCSIM,NOSIM,BCVECS,CREF,SVECS,
C    *              RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
      END IF
C
      MWOPH  = NWOPH
      NWOPH  = NWOPPT
C    ... tell SOLGO only to use the NWOPPT first JWOP entries
      DO 800 IOSIM = 1,NOSIM
         JPCMX  = KPCMX + (IOSIM-1)*NNORBX
         CALL SOLGO(D2,DV,WRK(JPCMX),SVEC(JSOVEC,IOSIM))
  800 CONTINUE
      NWOPH  = MWOPH
      IF (FULHES .AND. NCONRF.GT.1 .AND. LSYMRF.EQ.LSYMST) THEN
C
C        ... test orthogonality
C
         DO 900 IOSIM = 1,NOSIM
            TEST = DDOT(NCONRF,CREF,1,SVEC(1,IOSIM),1)
     *           - WRK(KOVLP-1+IOSIM)
            IF (ABS(TEST) .GT. 1.D-8) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(/A,I5,/A,1P,D12.4)')
     *            ' >>> SOLLNO WARNING, for IOSIM =',IOSIM,
     *            ' <CREF | SVEC_solvent(iosim) > =',TEST
            END IF
  900    CONTINUE
      END IF
C
C        ... test print
C
      IF (FULHES) THEN
         IF (IPRSOL .GE. 140) THEN
            WRITE (LUPRI,'(/A)') ' >>> PCMLNO - svec(ci,1) on exit'
            DO 930 I = 1,NCONST
               IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *         ' conf #',I,SVEC(I,1)
  930       CONTINUE
         END IF
      END IF
      IF (IPRSOL .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> PCMLNO - svec(orb) on exit'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *               NSVEC,NOSIM,1,LUPRI)
      END IF
C
      IF (LUPBKP .LT. 0) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('PCMLNO')
      RETURN
C     ... end of pcmlno.
      END

C  /* Deck NEQENR */
      SUBROUTINE NEQENR(ESOLT,POTMAT,QSEGR,POTGR,POT,WORK,LWRK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <pi.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0,
     &            DCVAL = D2, FPI = 4.0D0 * PI )
      DIMENSION POTMAT(NTS,*),QSEGR(*),POTGR(*),
     $     POT(*),WORK(*)

      CALL PCMINR(QSEGR,POTGR)
      VTOT0QS = D0
      VTOTQS = D0
      UELOR0 = D0
      UELOR  = D0
      FACFAS = (EPSINF - D1) / (EPS - D1)
      FACSLO = (EPS - EPSINF) / (EPS - D1)
      DO ITS = 1, NTS
         POTMAT(ITS,ITS) = 1.07D0 * DSQRT ( FPI/AS(ITS) )
         XI = XCTS(ITS)
         YI = YCTS(ITS)
         ZI = ZCTS(ITS)
         DO JTS = 1, NTS
            IF (JTS .NE. ITS) THEN
               XJ = XCTS(JTS)
               YJ = YCTS(JTS)
               ZJ = ZCTS(JTS)
               DRIJ = DSQRT((XI-XJ)**2 + (YI-YJ)**2 + (ZI-ZJ)**2)
               POTMAT(ITS,JTS) = D1/DRIJ
            END IF
         END DO
      END DO
      CALL SYMPCM(POTMAT,NTS,WORK,LWRK,NTSIRR)
      DO ITS = 1, NTSIRR
         SSI = D0
         SSN = D0
         DO JTS = 1, NTSIRR
C     bm
C     we need Q(el+N)^f for the excited state
C     calculated as sum of three contributions: Q(el)^ff + Q(N)^ff + Q(el+N)^fs
C     where Q(el)^ff=QSE and Q(N)^ff + Q(el+N)^fs = QSN - QSLO
c     bm
            QSLO  = FACSLO*QSEGR(JTS)
            SSI = SSI + QSEGR(JTS) * POTMAT(ITS,JTS)
            SSN = SSN - (QSE(JTS) + QSN(JTS) + QSLO) * POTMAT(ITS,JTS)
         ENDDO
Clf            IF (IPRPCM .GE. 140) THEN
         QSLO  = FACSLO*QSEGR(iTS)
         WRITE(LUPRI,*)'ITS,QSE(its),QSN(iTS),QSLO'
         WRITE(LUPRI,*) ITS,QSE(iTS),QSN(iTS),QSLO
Clf            END IF
C     
C     In the expression of the interaction of the slow charge with the solute, a nuclear part 
C     should be included (namely that related to the nuclear potential VNUC), however this is
C     here not computed as it cancels out
C     
         VTOT0QS = VTOT0QS - QSEGR(ITS) * POTGR(ITS)
         VTOTQS  = VTOTQS  - QSEGR(ITS) * POT(ITS)
         UELOR0  = UELOR0  + QSEGR(ITS) * FACFAS*SSI
         UELOR   = UELOR   + QSEGR(ITS) * SSN
      END DO
C
C  SYMMETRY NORMALIZATION OF UELOR E UEOLOR0
C
      UELOR  = UELOR  * (MAXREP + 1)
      UELOR0 = UELOR0 * (MAXREP + 1)


C     
      PORT =(VTOT0QS+UELOR0)*DP5*FACSLO
      QSSLW = (UELOR+VTOTQS)*DP5*FACSLO
c     lf         IF (IPRPCM .GE. 6) THEN
      WRITE(LUPRI,*)'GS(EL+N)-QS',VTOT0QS*FACSLO
      WRITE(LUPRI,*)'EXC(EL+N)-QS',VTOTQS*FACSLO
      WRITE(LUPRI,*)'QF(GS)-QS   = UELOR0 ',UELOR0*FACSLO
      WRITE(LUPRI,*)'QF(EXC)-QS = UELOR ',UELOR*FACSLO
      WRITE(LUPRI,*)'PORT =',PORT,' QSSLW = ',QSSLW
Clf      WRITE(LUPRI,*)'-TJEN',-TJEN
      WRITE(LUPRI,*)'PC=', PC, '  UNZ=', UNZ
Clf      WRITE(LUPRI,*)'PX=', -TJEN-PB
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)'ESOLT', ESOLT
      WRITE(LUPRI,*)'CORRECTION TO ESOLT', QSSLW-PORT
c     lf         END IF
C
      ESOLT = ESOLT + QSSLW - PORT
Clf      
      IF (IPRPCM .GE. 6) THEN         
         WRITE (LUPRI,*) 'Contributions to',
     $        'QSSLW = (UELOR+VTOTQS)*DP5*FACSLO'
         WRITE (LUPRI,*) 'UELOR,VTOTQS,FACSLO',
     $        UELOR,VTOTQS,FACSLO
         WRITE(LUPRI,*) 'ESOLT WITH NONEQ CONTRIB'
         WRITE(LUPRI,*) ESOLT, QSSLW, PORT 
      END IF
      RETURN
      END
