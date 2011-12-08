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
C===========================================================================
C
C File   : sirius/sirsol.F
C
C Module for spherical dielectric continuum model.
C
C Original code written in fall 1986 in Minneapolis, USA,
C when we all four visited Jan Almloef there.
C
C Original reference for the MCSCRF solvation model:
C          K. V. Mikkelsen, H. Agren, H. J. Aa. Jensen, and T. Helgaker.
C          J. Chem. Phys. 89 (1988) 3086-3095
C
C Intial polarization added later.
C
C===========================================================================
#ifdef REV_LOG
TO DO :: <  read CAVORG from LUSOL if label is new SOLVLRM
hjaaj: someone (Kurt, Trygve, Kenneth?) had changed THRZER in SOLGRD 1.0D-10
       in Oslo stating
       "Default threshold (in COMDECK THRZER) for zero is a bit too tight here"
       I don't see why that should be true so I have removed the change until
       further.  There must be a bug in SOLNC0 if abs(Tn(l,m)) .gt. 1.0D-14
       when symmetry of (l,m) is not totally symmetric.  /950824-hjaaj

Revision 1.2  2000/05/24 12:53:44  hjj
SOLFCK: corrected sign error in FCSOL introduced with AO code;
        fixed AO code for open shell and for all information
        needed for call WR_SIRIFC in sirdiis.F

980826-hjaaj: revised SOLLNC and SOLLNO for LSYMPT .ne. 1
950829-hjaaj: (ver n07) fixed error in SOLLNC (/INFVAR/ with NWOPH was missing)
950623/950501-hjaaj: *** version n06 ***
  new label on LUSOL (was SOLVINT, new SOLVRLM) [some SOLVINT forgotten 950501]
94-03-05 -kvm+ha: Korrigerat inertial polarisations del
930208-hjaaj: SOLLNO: fixed bug for NOSIM .gt. 1
930111-hjaaj: SOLFCK: zero WRK(KTAO) (bug fix)
921210-hjaaj -- preparation for ellipsoidal cavity: RSOL(3)
SOLFL : now general utility routine (nothing from common);
SOLFCK: new routine for calculating solvent t operator matrix
920714-hjaaj: SOLGRD: use F17.8 instead of F15.10 for T(lm) etc.
920601-hjaaj+ha+kvm -- implemented inertial polarization
SOLFL: EPSLON in parameter list
920522-hjaaj: added SOLREF deck
summer 1990-hjaaj: rewritten for determinant CI and symmetry
fall 1986 Minneapolis: original code
===========================================================================
#endif
#ifdef UNDEF
/* Comdeck hj_todo */
900723: timing of SOLGRD (the following is removed from SIROPT)
      PARAMETER (IDTGRD = 1)
C     ... note : IDTGRD = IDTIM in GRAD
#include <inftim.h>
...
         CALL GETTIM(T1,W1)
         ... call solgrd ...
         CALL GETTIM(T2,W2)
         TIMCPU(1,IDTGRD) = TIMCPU(1,IDTGRD) + T2 - T1
         TIMWAL(1,IDTGRD) = TIMWAL(1,IDTGRD) + W2 - W1
C        MAERKE: code 8 for SOLGRD contribution ???? 900723
900720: ORBLIN option to SOLDIA for no CSF diagonal ?? (for absorp.).
-- SOLGRD must be called by Abacus for solvent diagonal,
   i.e. NCONST,NWOPPT will be necessary. (or augment ABADIA).
NB. ABADIA should be in sirius, not abacus, to ensure correct argument
   list to sirius and detci routines.
-- Also consider -TELM contributions to diagonal Hessian and SOLLIN.
900718: absorption for solvent calculation.
#endif
C  /* Deck solref */
      SUBROUTINE SOLREF(IWUNIT)
C
C     Copyright 920522 Hans Joergen Aa. Jensen
C     Reference for solvent routines in Sirius;
C     called by SIRIUS.PRTINP
C
      WRITE(IWUNIT,100)
  100 FORMAT(/T2,'Original reference for MCSCRF solvation model:'
     *   /T5,'K.V.Mikkelsen, H.Agren, H.J.Aa.Jensen, and T. Helgaker',
     *   /T5,'J. Chem. Phys. 89 (1988) 3086-3095')
      RETURN
      END
C  /* Deck soldia */
      SUBROUTINE SOLDIA(TELM,RLMAC,INDXCI,RLM,DV,DIAG,WRK,LFREE)
C
C Copyright 29-Nov-1986 Hans Agren
C 25-Apr-1990 hjaaj -- detci version
C
C Purpose:
C  Calculate the diagonal of the solvent TE(LM) Hessian matrix for
C  a specific value of LM
C
C Input:
C   RLMAC: one-electron solvent integrals, active indices
C   DV   : one-electron density matrix (active part)
C   RLM  : one-electron solvent integrals
C
C Output:
C   DIAG: diagonal part of Te(LM) Hessian matrix
C
#include <implicit.h>
      DIMENSION RLM(*),RLMAC(*), DV(*), DIAG(*), WRK(*)
      DIMENSION INDXCI(*)
C
      PARAMETER (D2=2.0D0)
#include <dummy.h>
C
C Used from common blocks:
C   INFORB : NASHT
C   INFLIN : LSYMST,NCONST,NWOPPT
C   INFPRI : IPRSOL
C
#include <inforb.h>
#include <inflin.h>
#include <infpri.h>
C
      CALL QENTER('SOLDIA')
C
C     CSF part of diagonal.
C
      IF (NCONST .GT. 0) THEN
         CALL DZERO(DIAG,NCONST)
         IF (NASHT .GT. 1) THEN
            CALL CIDIAG(LSYMST,.TRUE.,RLMAC,DUMMY,INDXCI,DIAG,WRK,LFREE)
C           CALL CIDIAG(ICSYM,NOH2,FCAC,H2AC,XNDXCI,DIAGC,WRK,LFREE)
C
            DO 400 NA = 1,NCONST
               DIAG(NA) = D2 * (DIAG(NA) - TELM)
  400       CONTINUE
         END IF
      END IF
C
C     Now the orbital part of the diagonal:
C
      IF (NWOPPT .GT. 0) THEN
         CALL OR1DIA(DV,RLM,DIAG(1+NCONST),IPRSOL)
      END IF
C
      CALL QEXIT('SOLDIA')
      RETURN
      END
C  /* Deck solelm */
      FUNCTION SOLELM(DV,RLMAC,RLM,TELMAC)
C
C   Copyright 30-Nov-1986 Hans Jorgen Aa. Jensen
C   l.r. 6-May-1990 hjaaj
C
#include <implicit.h>
      DIMENSION DV(*), RLMAC(*), RLM(*)
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
C
      IF (NASHT .GT. 0) THEN
         TELM = D2*DDOT(NNASHX,DV,1,RLMAC,1)
         DO 100 I = 1,NASHT
            TELM = TELM - DV(IROW(I+1))*RLMAC(IROW(I+1))
  100    CONTINUE
      ELSE
         TELM   = D0
      END IF
      TELMAC = TELM
C
      DO 300 IW = 1,NISHT
         I = ISX(IW)
         TELM = TELM + D2*RLM(IROW(I+1))
  300 CONTINUE
C
      SOLELM = TELM
      RETURN
      END
C  /* Deck solfl */
      SUBROUTINE SOLFL(FLVEC,EPSLON,RSOL,LSOLMX)
C
C   Copyright 2-May-1987 Hans Joergen Aa. Jensen
C   Revised 10-Dec-1992 hjaaj (RSOL,LSOLMX in param.list)
C
C   Purpose: Calculate FLVEC, the f(l) factors in the solvent
C            energy expression for a surrounding medium,
C            cavity radius = Rsol and dielectric constant = EPSLON
C
#include <implicit.h>
      DIMENSION FLVEC(*),RSOL(3)
C
C     Parameters: FLFAC = - 1 / (8 pi epsilon_null)   (S.I.)
C                       = - 1 / 2                     (a.u.)
C
      PARAMETER ( FLFAC = -0.5D0 , D1 = 1.0D0 )
C
C
      IF (RSOL(1) .NE. RSOL(2) .OR. RSOL(1) .NE. RSOL(3)) THEN
         CALL QUIT('SOLFL error: only spherical cavity implemented')
      END IF
      RSOLAV = RSOL(1)
      IFL = 0
      DO 140 L = 0,LSOLMX
         RL = L
         FL = FLFAC * RSOLAV**(-(2*L+1))
     *              * (RL + D1) * (EPSLON - D1)
     *              / (RL + EPSLON*(RL + D1))
         DO 120 M = -L,L
            IFL = IFL + 1
            FLVEC(IFL) = FL
  120    CONTINUE
  140 CONTINUE
C
      RETURN
C     ... end of SOLFL.
      END
C  /* Deck solinr */
      SUBROUTINE SOLINR(FLVEC,FLINR,TLMSI)
C
C Copyright 1-Jun-1992 hjaaj+ha+kvm
C
C Calculate inertial f(l) factors;
C Read initial state T(lm) values
C
C Input : FLVEC with optical f(l) values
C Output: FLINR with inertial f(l) values
C         TLMSI with initial state T(lm) values for inertial contr.
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION FLVEC(NLMSOL), FLINR(NLMSOL), TLMSI(NLMSOL)
      PARAMETER (DM1 = -1.0D0)
C
C Used from common blocks:
C   INFINP: NLMSOL, LSOLMX, EPSTAT, RSOL(3)
C   INFORB: NSYM, NBAS(8)
C   INFTAP: LUSIFC, FNSIFC
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
C
      DIMENSION MBAS(8)
      LOGICAL FIRST
      SAVE    FIRST
      DATA    FIRST /.TRUE./
C
      CALL QENTER('SOLINR')
      CALL SOLFL(FLINR,EPSTAT,RSOL,LSOLMX)
      CALL DAXPY(NLMSOL,DM1,FLVEC,1,FLINR,1)
C
C     If (first call) then check if initial state information
C     is compatible with this calculation -- same geometry,
C     basis set, RSOL, EPSOL, LSOLMX, ...
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,FNSIFC,'OLD',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      IF (FIRST) THEN
         REWIND LUSIFC
         CALL MOLLAB('SOLVINFO',LUSIFC,lupri)
C        WRITE (LUSIFC) EPSOL,EPSTAT,RSOL,RSOL,RSOL,LSOLMX
C        WRITE (LUSIFC) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C        CALL WRITT (LUSIFC,NLMSO4,ERLM(1,1)))
C        CALL WRITT (LUSIFC,NLMSO4,ERLM(1,2)))
C        WRITE (LUSIFC) NSYM, NBAS
         READ (LUSIFC) EPSOLI,EPSTAI,RSOLX,RSOLY,RSOLZ,LSOLMI
         READ (LUSIFC) DUMINP,POTNUI
         READ (LUSIFC)
         READ (LUSIFC)
         READ (LUSIFC) NSYMI, MBAS
         NERR = 0
         IF (EPSOLI .NE. EPSTAT)  NERR = NERR + 1
         IF (EPSOLI .NE. EPSTAI)  NERR = NERR + 1
         IF (RSOLX  .NE. RSOL(1)) NERR = NERR + 1
         IF (RSOLY  .NE. RSOL(2)) NERR = NERR + 1
         IF (RSOLZ  .NE. RSOL(3)) NERR = NERR + 1
         IF (LSOLMI .NE. LSOLMX)  NERR = NERR + 1
         IF (POTNUI .NE. POTNUC)  NERR = NERR + 1
         IF (NSYMI  .NE. NSYM  )  NERR = NERR + 1
         DO 100 ISYM = 1,8
            IF (NBAS(ISYM) .NE. MBAS(ISYM)) NERR = NERR + 1
  100    CONTINUE
         IF (NERR .GT. 0) THEN
            WRITE (LUPRI,'(//A/A,2(/T20,A,T50,A))')
     *         ' SIRINR ERROR, ".INERSI" information on SIRIFC',
     *         '       is not compatible with this calculation',
     *         ' Information from SIRIFC','    This calculation',
     *         '========================','========================'
            WRITE (LUPRI,1001) 'EPSTAT' ,EPSOLI,EPSTAT
            WRITE (LUPRI,1001) 'RSOL_x' ,RSOLX ,RSOL(1)
            WRITE (LUPRI,1001) 'RSOL_y' ,RSOLY ,RSOL(2)
            WRITE (LUPRI,1001) 'RSOL_z' ,RSOLZ ,RSOL(3)
            WRITE (LUPRI,1001) 'POTNUC' ,POTNUI,POTNUC
            WRITE (LUPRI,1002) 'LSOLMX' ,LSOLMI,LSOLMX
            WRITE (LUPRI,1002) 'NSYM  ' ,NSYMI ,NSYM
            WRITE (LUPRI,1003) 'NBAS(8)',MBAS  ,NBAS
            IF (EPSOLI .NE. EPSTAI) WRITE (LUPRI,1011)
     &         'EPSOL .ne. EPSTAT on SIRIFC; EPSOL and EPSTAT =',
     &         EPSOLI,EPSTAI
            CALL QUIT('SIRINR: SIRIFC incompat. with .INERSF input')
         END IF
         FIRST = .FALSE.
      END IF
 1001 FORMAT(T4,A,T20,F20.10,T50,F20.10)
 1002 FORMAT(T4,A,T20,I20,T50,I20)
 1003 FORMAT(T4,A,T20,8I3,T50,8I3)
 1011 FORMAT(/T4,A,2F20.10)
C
C     Read TLMSI from LUSIFC
      REWIND LUSIFC
      CALL MOLLAB('SOLVINFO',LUSIFC,lupri)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (TLMSI(I),I=1,NLMSOL)
C
      CALL QEXIT('SOLINR')
      RETURN
      END
C  /* Deck solgc */
      SUBROUTINE SOLGC(CREF,RLMAC,TELMAC,GLMCI,INDXCI,WRK,LWRK)
C
C Copyright 30-Nov-1986 Ha.
C Completely rewritten for detci 25-Apr-1990 hjaaj
C
C Purpose:
C
C  To calculate solvent CI gradient contribution for a value of LM
C  (or several values summed) in RLMAC.
C  ( ER(l,m) configuration gradient contribution. )
C
C
C Input:
C  CREF(NCONF)   CI reference vector
C  RLMAC(NNASHX) solvent integrals
C  TELMAC        CREF expectation value.
C  INDXCI(*)     CI index vector
C
C Output:
C  GLMCI(NCONF)  = CI gradient vector for Rlm =
C                2 ( < i | R(l,m) | cref> - TE(l,m) cref(i) )
C
C Scratch:
C  WRK(LWRK)
C
#include <implicit.h>
      DIMENSION CREF(*),RLMAC(*),GLMCI(*),INDXCI(*),WRK(LWRK)
C
      PARAMETER ( D2 = 2.0D0 )
C
C Used from common blocks:
C   INFORB : N2ASHX
C   INFVAR : NCONF
C
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
C
      CALL QENTER('SOLGC ')
      KURAC = 1
      KWRK1 = KURAC + N2ASHX
      LWRK1 = LWRK + 1 - KWRK1
      CALL DSPTSI(NASHT,RLMAC,WRK(KURAC))
      CALL CIPRP(1,CREF,GLMCI,NCONF,WRK(KURAC),INDXCI,WRK(KWRK1),LWRK1)
C     CALL CIPRP(NSIM,CREF,SCVECS,LSCVEC,PRPAC,XNDXCI,WRK,LFREE)
C
      CALL DSCAL(NCONF,D2,GLMCI,1)
      CALL DAXPY(NCONF,(-D2*TELMAC),CREF,1,GLMCI,1)
C
      CALL QEXIT('SOLGC ')
      RETURN
      END
C  /* Deck solgo */
      SUBROUTINE SOLGO(DCVAL,DV,RLM,GOLM)
C
C  Copyright 29-Nov-1986 Hans Joergen Aa. Jensen
C  Revised for detci and no symmetry 6-May-1990 hjaaj.
C
C  Purpose: Add Te(l,m) orbital gradient contribution
C           to GOLM.
C
#include <implicit.h>
      DIMENSION DV(*), RLM(*), GOLM(*)
C
C  Used from common blocks:
C    INFORB: NISHT,NASHT
C    INFVAR: NWOPH,JWOP(2,*)
C    INFIND: IROW(*),IOBTYP(*),ISX(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
#include <infind.h>
C
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
      CALL QENTER('SOLGO ')
      DCFAC = D2 * DCVAL
C     ... DCVAL is the value of DC = <0| Eii |B>,
C         = 2 <0|0> = 2 normally,
C         = 2 <0|B> = 0 for transition density matrix
      DO 300 IG = 1,NWOPH
         K     = JWOP(1,IG)
         L     = JWOP(2,IG)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IF (K .LT. L) THEN
            KL = IROW(L) + K
         ELSE
            KL = IROW(K) + L
         END IF
         IF (ITYPK .EQ. JTINAC) THEN
C           first index inactive
            GOLM(IG) = GOLM(IG) + DCFAC * RLM(KL)
         ELSE
C           first index active
            NK   = ICH(K)
            TEMP = D0
            DO 100 NX = 1,NASHT
               IF (NX .LE. NK) THEN
                  DVKX  = DV(IROW(NK)+NX)
               ELSE
                  DVKX  = DV(IROW(NX)+NK)
               END IF
               IX = ISX(NISHT + NX)
               IF (IX .LE. L) THEN
                  RLMXL = RLM(IROW(L)+IX)
               ELSE
                  RLMXL = RLM(IROW(IX)+L)
               END IF
               TEMP = TEMP + DVKX*RLMXL
  100       CONTINUE
            GOLM(IG) = GOLM(IG) + D2 * TEMP
         END IF
         IF (ITYPL .EQ. JTACT) THEN
C           second index active
            NL   = ICH(L)
            TEMP = D0
            DO 200 NX = 1,NASHT
               IF (NX .LE. NL) THEN
                  DVLX  = DV(IROW(NL)+NX)
               ELSE
                  DVLX  = DV(IROW(NX)+NL)
               END IF
               IX = ISX(NISHT + NX)
               IF (IX .LE. K) THEN
                  RLMXK = RLM(IROW(K)+IX)
               ELSE
                  RLMXK = RLM(IROW(IX)+K)
               END IF
               TEMP = TEMP + DVLX*RLMXK
  200       CONTINUE
            GOLM(IG) = GOLM(IG) - D2 * TEMP
         END IF
  300 CONTINUE
C
      CALL QEXIT('SOLGO ')
      RETURN
C     end of solgo.
      END
C  /* Deck solgrd */
      SUBROUTINE SOLGRD(CREF,CMO,INDXCI,DV,G,ESOLT,ERLM,WRK,LFREE)
C
C   Copyright 29-Nov-1986 Hans Joergen Aa. Jensen
C
C   Purpose:  calculate MCSCF energy and gradient contribution
C             from a surrounding medium, cavity radius = Rsol
C             and dielectric constant = EPsol.
C
C   Output:
C    G          MCSCF gradient with solvation contribution added
C    ESOLT      total solvation energy
C    ERLM(lm,1) contains Esol(l,m) contribution to ESOLT
C    ERLM(lm,2) contains Tsol(l,m)
C
#include <implicit.h>
C
      DIMENSION CREF(*), CMO(*), INDXCI(*)
      DIMENSION DV(*),   G(*),   ERLM(NLMSOL,2),  WRK(*)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, DCVAL = D2 )
#include <thrzer.h>
#include <iratdef.h>
C
C Used from common blocks:
C   INFINP: NLMSOL, LSOLMX, EPSOL, RSOL(3)
C   INFVAR: NCONF,  NWOPT,  NVAR,   NVARH
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFIND: IROW(*)
C   INFTAP: LUSOL,  LUIT2
C   INFPRI: IPRSOL
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL     FIRST, FNDLAB
      PARAMETER  (MXLMAX = 50)
      DIMENSION   ISYRLM(2*MXLMAX+1)
      CHARACTER*8 STAR8, SOLVDI, EODATA
      SAVE        FIRST
      DATA        FIRST/.TRUE./, STAR8/'********'/
      DATA        SOLVDI/'SOLVDIAG'/, EODATA/'EODATA  '/
C
C     Statement functions;
C     define automatic arrays (dynamic core allocation)
C
      FLVEC(LM) = WRK(LM)
      FLINR(LM) = WRK(KFLINR-1+LM)
      TLMSI(LM) = WRK(KTLMSI-1+LM)
C
      CALL QENTER('SOLGRD')
C
      IF (LSOLMX .GT. MXLMAX) THEN
         WRITE (LUERR,*) 'ERROR SOLGRD, increase MXLMAX parameter'
         WRITE (LUERR,*) ' LSOLMX =',LSOLMX
         WRITE (LUERR,*) ' MXLMAX =',MXLMAX
         CALL QUIT('ERROR SOLGRD, increase MXLMAX parameter')
      END IF
C
C     Core allocation
C        FLVEC  f(l) factors in solvent energy expression
C        DIASH  diagonal of solvent contribution to Hessian
C        GRDLM  TELM gradient for current l,m value in the l,m loop
C        UCMO   CMO unpacked (i.e. no symmetry blocking)
C        RLMAC  active-active subblock of RLM
C        RLM    R(l,m) integrals for current l,m value in l,m loop
C
C     If (INERSF)
C       (i.e. If (inertial polarization contribution to final state))
C        FLINR  f(l) factors for inertial pol. contrib.
C        TLMSI  T(lm) values for initial state
C     end if
C
      KFLVEC = 1
C     ... NOTE: KFLVEC = 1 assumed in FLVEC(LM) definition above.
      IF (INERSF) THEN
         KFLINR = KFLVEC + NLMSOL
         KTLMSI = KFLINR + NLMSOL
         KDIASH = KTLMSI + NLMSOL
      ELSE
         KFLINR = 1
         KTLMSI = 1
         KDIASH = KFLVEC + NLMSOL
      END IF
      KGRDLM = KDIASH + NVAR
      KUCMO  = KGRDLM + NVARH
      KRLMAC = KUCMO  + NORBT*NBAST
      KRLM   = KRLMAC + NNASHX
      KW10   = KRLM   + NNORBX
C     1.1 read rlmao in ao basis and transform to rlm in mo basis
      KRLMAO = KW10
      KW20   = KRLMAO + NNBASX
C     1.2 diagonal contribution for current l,m value in the l,m loop
      KDIALM = KW10
      KW21   = KDIALM + NVAR
      LW21   = LFREE  - KW21
C     1.3 rest of CSF contribution
      KW22   = KW10
C
      KWRK1  = MAX(KW20,KW21,KW22)
      LWRK1  = LFREE  - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('SOLGRD',-KWRK1,LFREE)
C
         IF (IPRSOL .GE. 130) THEN
            WRITE (LUPRI,'(/A/A,2I10)')
     *         ' >>> SOLGRD - gtot (input) - non-zero elements',
     *         '     NCONF, NWOPT =',NCONF,NWOPT
            DO 40 I = 1,NCONF
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' conf #',I,G(I)
   40       CONTINUE
            DO 50 I = NCONF+1,NVAR
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' orb  #',I,G(I)
   50       CONTINUE
         END IF
C
C     Calculate f(l) factors
C     If (INERSF) FLVEC factors describe the optical polarization
C             and FLINR factors describe the inertial polarization
C     else        FLVEC may describe optical or static polarization
C
      CALL SOLFL(WRK(KFLVEC),EPSOL,RSOL,LSOLMX)
      IF (INERSF) THEN
         CALL SOLINR(WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI))
      END IF
      IF ((IPRSOL .GE. 5 .AND. FIRST) .OR. IPRSOL .GE. 15) THEN
         IF (.NOT.INERSF) THEN
            WRITE (LUPRI,'(//A/A)')
     *      ' >>> SOLGRD:  l     f(l) factor',
     *      '             === ================='
         ELSE
            WRITE (LUPRI,'(//A/A)')
     *      ' >>> SOLGRD:  l  optical f(l) factor inertial f(l) factor',
     *      '             === =================== ===================='
         END IF
         DO 140 L = 0,LSOLMX
            LL = (L+1)*(L+1)
            FL = FLVEC(LL)
            IF (INERSF) THEN
               FLI = FLINR(LL)
               WRITE (LUPRI,'(I15,F17.10,F21.10)') L, FL, FLI
            ELSE
               WRITE (LUPRI,'(I15,F16.10)') L, FL
            END IF
  140    CONTINUE
      END IF
C
C     Read and check dimension information (if first read) and
C     nuclear contributions to ERLM (always).
C
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,lupri)
      IF (FIRST) THEN
         READ (LUSOL) LMAXSS, LMTOT, NNNBAS
         IF (LMAXSS .LT. LSOLMX) THEN
            WRITE (LUERR,'(//2A,2(/A,I5))') ' >>> SOLGRD ERROR,',
     *      ' insufficient number of intgrals on LUSOL',
     *      ' l max from SIRIUS input :',LSOLMX,
     *      ' l max from LUSOL  file  :',LMAXSS
            CALL QUIT('SOLGRD: lmax on LUSOL is too small')
         END IF
         IF ((LMAXSS+1)**2 .NE. LMTOT) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> SOLGRD ERROR,',
     *      ' LUSOL file info inconsistent',
     *      ' l_max               :',LMAXSS,
     *      ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *      ' LMTOT               :',LMTOT
            CALL QUIT('SOLGRD: LUSOL info not internally consistent')
         END IF
         IF (NNNBAS .NE. NBAST) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> SOLGRD ERROR,',
     *      ' LUSOL file info inconsistent with SIRIUS input',
     *      ' NBAST - LUSOL       :',NNNBAS,
     *      ' NBAST - SIRIUS      :',NBAST
            CALL QUIT('SOLGRD: LUSOL info not '//
     &                'consistent with SIRIUS input.')
         END IF
      ELSE
         READ (LUSOL)
      END IF
      CALL READT(LUSOL,NLMSOL,ERLM(1,2))
C
      IF (IPRSOL .GE. 20 .AND. NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLGRD - DV matrix :'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF
      IF (IPRSOL .GE. 7) THEN
         WRITE (LUPRI,'(/A/)')
     *      ' l, m, Tn(lm) - the nuclear contributions :'
         LM = 0
         DO 220 L = 0,LSOLMX
            DO 210 M = -L,L
               LM = LM + 1
               WRITE (LUPRI,'(2I5,F15.10)') L,M,ERLM(LM,2)
  210       CONTINUE
            WRITE (LUPRI,'()')
  220    CONTINUE
      END IF
C
C     Unpack symmetry blocked CMO
C     Loop over l,m expansion
C
      CALL UPKCMO(CMO,WRK(KUCMO))
      IF (IPRSOL .GE. 6)
     *   WRITE (LUPRI, '(//A/)') ' >>> SOLGRD: START LOOP OVER LM'
      CALL DZERO(WRK(KDIASH),NVAR)
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
            WRITE (LUERR,*) 'ERROR SOLGRD: L from LUSOL not as expected'
            WRITE (LUERR,*) 'L from 520 loop:',L
            WRITE (LUERR,*) 'L from LUSOL   :',L1
            CALL QUIT('ERROR SOLGRD: L from LUSOL not as expected')
         END IF
      DO 500 M = -L,L
         LM = LM + 1
         IF (IPRSOL .GE. 15) THEN
            WRITE (LUPRI,'(/A,2I5/A)') ' >>> l, m :',L,M,
     *                                ' ===================='
            WRITE (LUPRI,'(A,I2)') ' Symmetry :',ABS(ISYRLM(L+M+1))
         END IF
         IF (ISYRLM(L+M+1) .NE. 1) THEN
Chj         IF (ABS(ERLM(LM,2)) .GT. THRZER) THEN
Chj numerical errors have seen to give 1.0D-12,
Chj so we now use 1.0D-10 for test /Feb 2005
            IF (ABS(ERLM(LM,2)) .GT. 1.0D-10) THEN
               WRITE (LUPRI,*) 'ERROR SOLGRD for l,m',L,M
               WRITE (LUPRI,*) 'Symmetry :',ISYRLM(L+M+1)
               WRITE (LUPRI,*) 'Tn(l,m) .ne. 0, but =',ERLM(LM,2)
               CALL QUIT( 'ERROR SOLGRD: Tn(l,m) not 0 as expected')
            END IF
            ERLM(LM,2) = D0
C           ... to fix round-off errors in Tn(l,m) calculation
            IF (ISYRLM(L+M+1) .GT. 1) READ (LUSOL)
            GO TO 500
         END IF
C
C        Read R(l,m) in ao basis and transform to mo basis.
C        Extract active-active block in RLMAC(1) = WRK(KRLMAC).
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KRLM),WRK(KUCMO),WRK(KWRK1),
     &             NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KRLM),WRK(KRLMAC))
         END IF
         IF (IPRSOL .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' Rlm_ao matrix:'
            CALL OUTPAK(WRK(KRLMAO),NBAST,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' Rlm_mo matrix:'
            CALL OUTPAK(WRK(KRLM),  NORBT,1,LUPRI)
            IF (NASHT .GT. 0) THEN
               WRITE (LUPRI,'(/A)') ' Rlm_ac matrix:'
               CALL OUTPAK(WRK(KRLMAC),NASHT,1,LUPRI)
            END IF
         END IF
C
C        Add electronic contribution TE(l,m) to T(l,m)
C
         TELM     = SOLELM(DV,WRK(KRLMAC),WRK(KRLM),TELMAC)
         IF (IPRSOL .GE. 6) THEN
            WRITE (LUPRI,'(A,2I5,/A,3F17.8)')
     *      ' >>> l, m :',L,M,
     *      '     Te(lm), Tn(lm), T(lm) :',
     *         TELM,ERLM(LM,2),ERLM(LM,2)-TELM
            IF (IPRSOL .GE. 10) WRITE (LUPRI,'(A,F17.8)')
     *      ' --- active part of Te(lm) :',TELMAC
            IF (INERSF) WRITE (LUPRI,'(A,F17.8)')
     *      ' --- inertial T(lm) value  :',TLMSI(LM)
         END IF
         ERLM(LM,2) = ERLM(LM,2) - TELM
      IF (ABS(ERLM(LM,2)) .LE. THRZER) THEN
         ERLM(LM,2) = D0
         GO TO 500
      END IF
C     ... test introduced 880109 hjaaj
C         (the only possible problem is the DO 420 loop,
C          but I think (w.o. having checked) that this
C          contribution to the Hessian diagonal also will be
C          zero if ERLM(LM,2) zero).
C
C        Calculate orbital TE(l,m) gradient contribution
C        and part of csf contribution.
C
         CALL DZERO(WRK(KGRDLM),NVARH)
         IF (NCONF .GT. 1) THEN
            CALL SOLGC(CREF,WRK(KRLMAC),TELMAC,WRK(KGRDLM),INDXCI,
     &                 WRK(KWRK1),LWRK1)
C           CALL SOLGC(CREF,RLMAC,TELMAC,GLMCI,INDXCI,WRK,LWRK)
         END IF
         IF (NWOPT .GT. 0) THEN
            CALL SOLGO(DCVAL,DV,WRK(KRLM),WRK(KGRDLM+NCONF))
         END IF
C
C
C        Obtain DIALM = diagonal TE(l,m) Hessian
C                     = 2 ( <i|R(l,m)|i> - TE(l,m) )
C        Add the DIALM contribution and the GRDLM contribution
C        to solvent Hessian diagonal.
C
         CALL SOLDIA(TELMAC,WRK(KRLMAC),INDXCI,
     *               WRK(KRLM),DV,WRK(KDIALM),WRK(KW21),LW21)
C        CALL SOLDIA(TELM,RLMAC,INDXCI,RLM,DV,DIAG,WRK,LFREE)
C
         FAC1 = - D2 * FLVEC(LM) * ERLM(LM,2)
         IF (INERSF) THEN
            FAC1 = FAC1 - FLINR(LM) * D2 * TLMSI(LM)
         END IF
         FAC2 =   D2 * FLVEC(LM)
         DO 420 I = 0,(NVAR-1)
            WRK(KDIASH+I) = WRK(KDIASH+I)
     *                    + FAC1 * WRK(KDIALM+I)
     *                    + FAC2 * WRK(KGRDLM+I) * WRK(KGRDLM+I)
  420    CONTINUE
C
C        test orthogonality
C
         IF (IPRSOL .GE. 120) THEN
           WRITE (LUPRI,'(/A)')' >>> SOLGRD - grdlm, dialm, diash, cref'
           DO 430 I = 1,NCONF
              WRITE (LUPRI,'(A,I10,4F12.6)') ' conf #',I,
     *        WRK(KGRDLM-1+I),WRK(KDIALM-1+I),WRK(KDIASH-1+I),CREF(I)
  430      CONTINUE
         END IF
         TEST = DDOT(NCONF,CREF,1,WRK(KGRDLM),1)
         IF (ABS(TEST) .GT. 1.D-8) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A,I5,/A,1P,D12.4)')
     *      ' >>> SOLGRD WARNING, for LM =',LM,
     *      ' <CREF | GRADlm > =',TEST
         END IF
C
C        Add TE(l,m) gradient contribution to MCSCF gradient
C        g  =  g  -  2 f(l) * T(l,m) * (dTE(l,m)/d(lambda))
C
         FAC = - D2 * FLVEC(LM) * ERLM(LM,2)
         IF (INERSF) THEN
            FAC = FAC - FLINR(LM) * D2 * TLMSI(LM)
         END IF
         CALL DAXPY(NVARH,FAC,WRK(KGRDLM),1,G,1)
         IF (IPRSOL .GE. 140) THEN
            WRITE (LUPRI,'(/A/A,2I10)')
     *         ' >>> SOLGRD - grdlm, gtot (accum) - non-zero grdlm',
     *         '     NCONF, NWOPT =',NCONF,NWOPT
            DO 440 I = 1,NCONF
               IF (WRK(KGRDLM-1+I) .NE. D0)
     *            WRITE (LUPRI,'(A,I10,3F15.10)')
     *            ' conf #',I,FAC*WRK(KGRDLM-1+I),G(I)
  440       CONTINUE
            DO 450 I = NCONF+1,NVAR
               IF (WRK(KGRDLM-1+I) .NE. D0)
     *            WRITE (LUPRI,'(A,I10,3F15.10)')
     *            ' orb  #',I,FAC*WRK(KGRDLM-1+I),G(I)
  450       CONTINUE
         END IF
C
  500 CONTINUE
  520 CONTINUE
C
C     500 is end of (l,m) loop.
C
C
         IF (IPRSOL .GE. 130) THEN
            WRITE (LUPRI,'(/A/A,2I10)')
     *         ' >>> SOLGRD - gtot (output) - non-zero elements',
     *         '     NCONF, NWOPT =',NCONF,NWOPT
            DO 840 I = 1,NCONF
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' conf #',I,G(I)
  840       CONTINUE
            DO 850 I = NCONF+1,NVAR
               IF (G(I) .NE. D0) WRITE (LUPRI,'(A,I10,3F15.10)')
     *         ' orb  #',I,G(I)
  850       CONTINUE
         END IF
C
C
C     Calculate ER(l,m) energy contributions and add them up
C
      ESOLT = D0
      DO 900 LM = 1,NLMSOL
         ERLM(LM,1) = FLVEC(LM) * ERLM(LM,2) * ERLM(LM,2)
         IF (INERSF) THEN
            ERLM(LM,1) = ERLM(LM,1)
     *                 + FLINR(LM) * ERLM(LM,2) * D2 * TLMSI(LM)
     *                 - FLINR(LM) * TLMSI(LM) * TLMSI(LM)
         END IF
         ESOLT    = ESOLT     + ERLM(LM,1)
  900 CONTINUE
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
      CALL QEXIT('SOLGRD')
      RETURN
C     end of solgrd.
      END
C  /* Deck solfck */
      SUBROUTINE SOLFCK(DCAO,DVAO,FSOL,CMO,TLM,MOBAS,ESOLT,
     &                  WRK,LFREE,IPRINT)
C
C   Copyright 10-Dec-1992 Hans Joergen Aa. Jensen
C   (based on SOLGRD)
C   Modified Dec.-98 to allow the Fock operator to be returned
C   in AO basis, bmkvmkr.
C
C   Purpose:  calculate solvent t operator matrix, where
C             t_pq = -2 sum(lm) g_lm T_lm R_lm_pq
C
C   Output:
C    FSOL     MCSCF gradient with solvation contribution added
C
#include <implicit.h>
C
      DIMENSION FSOL(*),CMO(*),TLM(NLMSOL,2),DCAO(*),DVAO(*),WRK(LFREE)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0)
#include <thrzer.h>
#include <iratdef.h>
C
C Used from common blocks:
C   INFINP: NLMSOL, LSOLMX, EPSOL, RSOL(3)
C   INFORB: NNASHX, NNBASX, NNORBX, etc.
C   INFTAP: LUSOL
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <scbrhf.h>
C
      LOGICAL     FIRST, MOBAS
      PARAMETER  (MXLMAX = 50)
      DIMENSION   ISYRLM(2*MXLMAX+1)
      SAVE        FIRST
      DATA        FIRST/.TRUE./
C
C     Statement functions;
C     define automatic arrays (dynamic core allocation)
C
      FLVEC(LM) = WRK(LM)
      FLINR(LM) = WRK(KFLINR-1+LM)
      TLMSI(LM) = WRK(KTLMSI-1+LM)
C
      CALL QENTER('SOLFCK')
C
      IF (LSOLMX .GT. MXLMAX) THEN
         WRITE (LUERR,*) 'ERROR SOLFCK, increase MXLMAX parameter'
         WRITE (LUERR,*) ' LSOLMX =',LSOLMX
         WRITE (LUERR,*) ' MXLMAX =',MXLMAX
         CALL QUIT('ERROR SOLFCK, increase MXLMAX parameter')
      END IF
C
C     Core allocation
C        FLVEC  f(l) factors in solvent energy expression
C        UCMO   CMO unpacked (i.e. no symmetry blocking)
C        RLM    R(l,m) integrals for current l,m value in l,m loop
C
C     If (INERSF)
C       (i.e. If (inertial polarization contribution to final state))
C        FLINR  f(l) factors for inertial pol. contrib.
C        TLMSI  T(lm) values for initial state
C     end if
C
      KFLVEC = 1
C     ... NOTE: KFLVEC = 1 assumed in FLVEC(LM) definition above.
      IF (INERSF) THEN
         KFLINR = KFLVEC + NLMSOL
         KTLMSI = KFLINR + NLMSOL
         KTAO   = KTLMSI + NLMSOL
      ELSE
         KFLINR = 1
         KTLMSI = 1
         KTAO   = KFLVEC + NLMSOL
      END IF
      KRLMAO = KTAO   + NNBASX
      KWRK1  = KRLMAO + NNBASX
C     1.1 read rlmao in ao basis and add to tao in ao basis
      IF (MOBAS) THEN
         KTNLM  = -999 999 999
         KDTAO  = -999 999 999
      ELSE
         KTNLM  = KWRK1  + NNBASX
         KDTAO  = KTNLM  + NLMSOL
         KWRK1  = KDTAO  + NNBASX
      END IF
      LWRK1  = LFREE  - KWRK1
      IF (LWRK1 .LT. 0) CALL ERRWRK('SOLFCK',-KWRK1,LFREE)
C
C
C     Calculate f(l) factors
C     If (INERSF) FLVEC factors describe the optical polarization
C             and FLINR factors describe the inertial polarization
C     else        FLVEC may describe optical or static polarization
C
      CALL SOLFL(WRK(KFLVEC),EPSOL,RSOL,LSOLMX)
      IF (INERSF) THEN
         CALL SOLINR(WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI))
      END IF
      IF (IPRINT .GE. 5 .AND. FIRST) THEN
         IF (.NOT.INERSF) THEN
            WRITE (LUPRI,'(//A/A)')
     *      ' >>> SOLFCK:  l     f(l) factor',
     *      '             === ================='
         ELSE
            WRITE (LUPRI,'(//A/A)')
     *      ' >>> SOLFCK:  l  optical f(l) factor inertial f(l) factor',
     *      '             === =================== ===================='
         END IF
         DO 140 L = 0,LSOLMX
            LL = (L+1)*(L+1)
            FL = FLVEC(LL)
            IF (INERSF) THEN
               FLI = FLINR(LL)
               WRITE (LUPRI,'(I15,F17.10,F21.10)') L, FL, FLI
            ELSE
               WRITE (LUPRI,'(I15,F16.10)') L, FL
            END IF
  140    CONTINUE
      END IF
C
C     Read and check dimension information (if first read)
C
C     Read nuclear contributions to ERLM (TN(l,m))
C     and position for R(l,m).
C
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,lupri)
      IF (FIRST) THEN
         READ (LUSOL) LMAXSS, LMTOT, NNNBAS
         IF (LMAXSS .LT. LSOLMX) THEN
            WRITE (LUERR,'(//2A,2(/A,I5))') ' >>> SOLFCK ERROR,',
     *      ' insufficient number of intgrals on LUSOL',
     *      ' l max from SIRIUS input :',LSOLMX,
     *      ' l max from LUSOL  file  :',LMAXSS
            CALL QUIT('SOLFCK: lmax on LUSOL is too small')
         END IF
         IF ((LMAXSS+1)**2 .NE. LMTOT) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> SOLFCK ERROR,',
     *      ' LUSOL file info inconsistent',
     *      ' l_max               :',LMAXSS,
     *      ' (l_max + 1) ** 2    :',(LMAXSS+1)**2,
     *      ' LMTOT               :',LMTOT
            CALL QUIT('SOLFCK: LUSOL info not internally consistent')
         END IF
         IF (NNNBAS .NE. NBAST) THEN
            WRITE (LUERR,'(//2A,3(/A,I5))') ' >>> SOLFCK ERROR,',
     *      ' LUSOL file info inconsistent with SIRIUS input',
     *      ' NBAST - LUSOL       :',NNNBAS,
     *      ' NBAST - SIRIUS      :',NBAST
            CALL QUIT('SOLFCK: LUSOL info not consistent '//
     &                'with SIRIUS input.')
         END IF
      ELSE
         READ (LUSOL)
      END IF
      IF (MOBAS) THEN
         READ (LUSOL)
C     ... skip Tn(lm)
      ELSE
         CALL READT(LUSOL,NLMSOL,WRK(KTNLM))
C        get total dens.mat. DTAO = DCAO + DVAO
         IF (NASHT .EQ. 0) THEN
            CALL PKSYM1(WRK(KDTAO),DCAO,NBAS,NSYM,-1)
         ELSE
            DO I = 1,NNBAST
               WRK(KTAO-1+I) = DCAO(I) + DVAO(I)
            END DO
            CALL PKSYM1(WRK(KDTAO),WRK(KTAO),NBAS,NSYM,-1)
         END IF
      END IF
C
C     Loop over l,m expansion
C
      IF (IPRINT .GE. 10)
     *   WRITE (LUPRI, '(//A)') ' >>> SOLFCK: START LOOP OVER LM'
      CALL DZERO(WRK(KTAO),NNBASX)
      LM = 0
      ESOLT = D0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
         IF (L1 .NE. L) THEN
            WRITE (LUPRI,*) 'ERROR SOLFCK: L from LUSOL not as expected'
            WRITE (LUPRI,*) 'L from 520 loop:',L
            WRITE (LUPRI,*) 'L from LUSOL   :',L1
            CALL QUIT('ERROR SOLFCK: L from LUSOL not as expected')
         END IF
      DO 500 M = -L,L
         LM = LM + 1
         IF (MOBAS) THEN
            FAC = - D2 * FLVEC(LM) * TLM(LM,2)
         ELSE
            FAC = D0
         END IF
         IF (INERSF) THEN
            FAC = FAC - FLINR(LM) * D2 * TLMSI(LM)
         END IF
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI,'(/A,2I5/A/A,I2)')
     &         ' >>> l, m :',L,M,
     &         ' ====================',
     &         ' Symmetry   :',ABS(ISYRLM(L+M+1))
            IF (MOBAS) WRITE (LUPRI,'(A,1P,D16.8)')
     &         ' T(l,m)     :',TLM(LM,2),
     &         ' Factor     :',FAC
         END IF
         IF (ISYRLM(L+M+1) .LT. 1) THEN
            WRITE(LUPRI,*) 'SOLFCK FATAL ERROR for l,m:',L,M
            WRITE(LUPRI,*) 'ISYRLM(l,m) =',ISYRLM(L+M+1)
            WRITE(LUPRI,*) 'This RLMAO not available on LUSOL !'
            CALL QUIT('SOLFCK error: a needed Rlmao not on LUSOL')
         END IF
C
C        Read R(l,m) in ao basis and add to tao
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         IF (IPRINT .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' Rlm_ao matrix:'
            CALL OUTPAK(WRK(KRLMAO),NBAST,1,LUPRI)
         END IF
C
         IF (.NOT. MOBAS) THEN
            TLMEXP = WRK(KTNLM + LM -1)
     &             - DDOT(NNBASX,WRK(KDTAO),1,WRK(KRLMAO),1)
            FAC = FAC - D2 * FLVEC(LM) * TLMEXP
            TLM(LM,2) = TLMEXP
            TLM(LM,1) = FLVEC(LM) * TLMEXP * TLMEXP
            IF (INERSF) THEN
               TLM(LM,1) = TLM(LM,1)
     *                   + FLINR(LM) * TLM(LM,2) * D2 * TLMSI(LM)
     *                   - FLINR(LM) * TLMSI(LM) * TLMSI(LM)
            END IF
            ESOLT = ESOLT + TLM(LM,1)
            IF (IPRINT .GE. 10) WRITE (LUPRI,'(A,1P,2D16.8)')
     &         ' T(l,m,1:2) :',TLM(LM,1),TLM(LM,2),
     &         ' Factor     :',FAC
         END IF
         CALL DAXPY(NNBASX,FAC,WRK(KRLMAO),1,WRK(KTAO),1)
C
  500 CONTINUE
  520 CONTINUE
C
C     500 is end of (l,m) loop.
C
C
C     Symmetry pack TAO using WRK(KRLMAO)
C     Transform to MO basis
C
      CALL PKSYM1(WRK(KTAO),WRK(KRLMAO),NBAS,NSYM,+1)
      IF (IPRINT .GE. 15) THEN
         WRITE (LUPRI,'(/A)')
     &      ' SOLFCK effective solvent t_ao matrix (not packed):'
         CALL OUTPAK(WRK(KTAO),NBAST,1,LUPRI)
         WRITE (LUPRI,'(/A)')
     &      ' SOLFCK effective solvent t_ao matrix (packed):'
         CALL OUTPKB(WRK(KRLMAO),NBAS,NSYM,1,LUPRI)
      END IF
      IF (MOBAS) THEN 
         DO 700 ISYM = 1,NSYM
            JCMO = ICMO(ISYM) + 1
            J1AO = KRLMAO + IIBAS(ISYM)
            J1MO = 1 + IIORB(ISYM)
            CALL UTHU(WRK(J1AO),FSOL(J1MO),CMO(JCMO),
     &                WRK(KWRK1),NBAS(ISYM),NORB(ISYM))
 700     CONTINUE
      ELSE 
         CALL DCOPY(NNBAST,WRK(KRLMAO),1,FSOL,1)
      END IF
C
C
      IF (IPRINT .GE. 6 .AND. MOBAS) THEN
         WRITE (LUPRI,'(/A)')
     &      ' SOLFCK: effective solvent one-electron matrix:'
         CALL OUTPKB(FSOL,NORB,NSYM,1,LUPRI)
      END IF
      FIRST = .FALSE.
      CALL QEXIT('SOLFCK')
      RETURN
C     end of solfck.
      END
C  /* Deck sollin */
      SUBROUTINE SOLLIN(NCSIM,NOSIM,BCVECS,BOVECS,CREF,CMO,INDXCI,
     &                  DV,DTV,SCVECS,SOVECS,ORBLIN,WRK,LWRK)
C
C Copyright 25-Apr-1990 Hans Joergen Aa. Jensen
C Common driver for SOLLNC and SOLLNO
C
#include <implicit.h>
      DIMENSION BCVECS(*),BOVECS(*),CREF(*),CMO(*),INDXCI(*)
      DIMENSION DV(*),DTV(*),SCVECS(*),SOVECS(*),WRK(LWRK)
      LOGICAL   ORBLIN
C
C Used from common blocks:
C   INFINP : LSOLMX,NLMSOL, EPSOL, INERSF, RSOL
C   INFLIN : NWOPPT,NVARPT
C
#include <maxorb.h>
#include <infinp.h>
#include <inflin.h>
C
      CALL QENTER('SOLLIN')
      KFRSAV = 1
      LFRSAV = LWRK
      KFREE  = KFRSAV
      LFREE  = LWRK
C
C     Calculate f(l) factors;
C     Rsol and EPsol dependence is located here.
C
      CALL MEMGET('REAL',KFLVEC,NLMSOL,WRK,KFREE,LFREE)
      CALL SOLFL(WRK(KFLVEC),EPSOL,RSOL,LSOLMX)
      IF (INERSF) THEN
         CALL MEMGET('REAL',KFLINR,NLMSOL,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLMSI,NLMSOL,WRK,KFREE,LFREE)
         CALL SOLINR(WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI))
      ELSE
         CALL MEMGET('REAL',KFLINR,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLMSI,0,WRK,KFREE,LFREE)
      END IF
      CALL MEMGET('INTE',KSYRLM,(2*LSOLMX+1),WRK,KFREE,LFREE)
C
      IF (NCSIM .GT. 0) THEN
         CALL SOLLNC(NCSIM,BCVECS,CREF,CMO,INDXCI,
     *               WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI),
     *               DV,DTV,SCVECS,WRK(KSYRLM),WRK(KFREE),LFREE)
C        CALL SOLLNC(NCSIM,BCVEC,CREF,CMO,INDXCI,
C    *               FLVEC,FLINR,TLMSI,
C    *               DV,DTV,SVEC,ISYRLM,WRK,LFREE)
C
      END IF
      IF ( NOSIM .GT.0 ) THEN
         IF ( .NOT.ORBLIN ) THEN
            NSVAR  = NVARPT
         ELSE
            NSVAR  = NWOPPT
         END IF
         CALL SOLLNO(NOSIM,BOVECS,CREF,CMO,INDXCI,
     *               WRK(KFLVEC),WRK(KFLINR),WRK(KTLMSI),
     *               DV, SOVECS,NSVAR,WRK(KSYRLM),WRK(KFREE),LFREE)
C        CALL SOLLNO(NOSIM,BOVEC,CREF,CMO,INDXCI,
C    *               FLVEC,FLINR,TLMSI,
C    *               DV, SVEC,NSVEC,ISYRLM,WRK,LFREE)
      END IF
      CALL MEMREL('SOLLIN',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SOLLIN')
      RETURN
      END
C  /* Deck sollnc */
      SUBROUTINE SOLLNC(NCSIM,BCVEC,CREF,CMO,INDXCI,
     *                  FLVEC,FLINR,TLMSI,
     *                  DV,DTV,SVEC,ISYRLM,WRK,LFREE)
C
C  Copyright 24-Jan-1987 Hans Jorgen Aa. Jensen
C  Revised for detci 17-Jul-1990 hjaaj
C  Revised for LSYMPT .ne. 1 26-Aug-1998 hjaaj
C
C  Purpose:  Calculate MCSCF Hessian contribution from a
C            surrounding medium to a csf trial vector.
C            Cavity radius = Rsol and dielectric constant = EPsol.
C
#include <implicit.h>
      DIMENSION BCVEC(*),  CREF(*), CMO(*),  ISYRLM(*)
      DIMENSION FLVEC(NLMSOL), FLINR(NLMSOL), TLMSI(NLMSOL)
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
C    INFINP : NLMSOL, LSOLMX, INERSF
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : NWOPH
C    INFLIN : NCONST, NVARPT, NWOPPT
C    INFTAP : LUSOL
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inflin.h>
#include <inftap.h>
#include <infpri.h>
C
C     Statement functions
C
      TN(   LM) = WRK((KTNLM -1)+LM)
C
      CALL QENTER('SOLLNC')
C
C     Core allocation
C
      KRXC   = 1
      KRYC   = KRXC   + NCSIM*NNORBX
      KTXCAC = KRYC   +       NNORBX
      KW10   = KTXCAC + NCSIM
C
      KTNLM  = KW10
      KUCMO  = KTNLM  + NLMSOL
      KRLMAC = KUCMO  + NORBT*NBAST
      KRLM   = KRLMAC + NNASHX
      KW20   = KRLM   + NNORBX
C     2.1 read rlmao in ao basis and transform to rlm in mo basis
      KRLMAO = KW20
      KW21   = KRLMAO + NNBASX
C     3.0 SOLSC
      KRXCAC = KW10
      KRYCAC = KRXCAC + NCSIM*NNASHX
      KW30   = KRYCAC +       NNASHX
      LW30   = LFREE  + 1 - KW30
C
      KNEED  = MAX(KW21,KW30)
      IF (KNEED .GT. LFREE) CALL ERRWRK('SOLLNC',KNEED,LFREE)
C
      IF (IPRSOL .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLLNC - svec(ci,1) on entry'
         DO I = 1,NCONST
            IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *      ' conf #',I,SVEC(I,1)
         END DO
         WRITE (LUPRI,'(/A)') ' >>> SOLLNC - svec(orb) on entry'
         CALL OUTPUT(SVEC(1+NCONST,1),1,NWOPPT,1,NCSIM,
     *               NSVEC,NCSIM,1,LUPRI)
      END IF
C
C     Read nuclear contributions to ERLM (TN(l,m))
C     and position for R(l,m).
C
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,lupri)
      READ (LUSOL)
      CALL READT(LUSOL,NLMSOL,WRK(KTNLM))
C
C
C     Construct effective operators Rxc and Ryc
C     =========================================
C
C     Zero Rxc and Ryc storage
C
      CALL DZERO(WRK(KRXC), NCSIM * NNORBX )
      CALL DZERO(WRK(KRYC),         NNORBX )
      TYCAC = D0
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Loop over l,m expansion.
C
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
      DO 500 M = -L,L
         LM = LM + 1
         IF (ISYRLM(L+M+1) .NE. LSYMPT .AND. ISYRLM(L+M+1) .NE. 1) THEN
            IF (ISYRLM(L+M+1) .GT. 0) THEN
               READ (LUSOL)
            ELSE IF (ISYRLM(L+M+1) .EQ. -LSYMPT
     &          .OR. ISYRLM(L+M+1) .EQ. -1) THEN
             WRITE(LUPRI,*)'SOLLNC ERROR, needed RLM not on LUSOL file!'
             WRITE(LUPRI,*)'- You must not specify ".NOT ALL RLM" under'
             WRITE(LUPRI,*)'  "*ONEINT" when running ABACUS or RESPONSE'
             WRITE(LUPRI,*)'  with non-totally symmetric perturbations.'
             WRITE(LUPRI,*)'L,M,ISYRLM(l,m), LSYMPT',
     &                    L,M,ISYRLM(L+M+1),LSYMPT
             CALL QUIT('SOLLNC ERROR, needed RLM not on file!')
            END IF
            GO TO 500
         END IF
C
C        Read R(l,m) in ao basis and transform to mo basis.
C        Extract active-active block in RLMAC = WRK(KRLMAC).
C        Electronic contribution TE(l,m).
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KRLM),WRK(KUCMO),WRK(KW21),
     *             NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KRLM),WRK(KRLMAC))
         END IF
         TLM   = TN(LM) - SOLELM(DV,WRK(KRLMAC),WRK(KRLM),TELMAC)
C
         JRXC  = KRXC
         DO 200 ICSIM = 1,NCSIM
            TTT = SOLELM(DTV(1,ICSIM),WRK(KRLMAC),WRK(KRLM),TAC)
            FACX  = D2 * FLVEC(LM) * TAC
            IF (ABS(FACX) .GT. THRZER)
     *         CALL DAXPY(NNORBX,FACX,WRK(KRLM),1,WRK(JRXC),1)
            JRXC  = JRXC + NNORBX
  200    CONTINUE
C
         FACY  = - D2 * FLVEC(LM) * TLM
         IF (INERSF) THEN
            FACY = FACY - FLINR(LM) * D2 * TLMSI(LM)
         END IF
         TYCAC = TYCAC + FACY*TELMAC
         IF (ABS(FACY) .GT. THRZER)
     *      CALL DAXPY(NNORBX,FACY,WRK(KRLM),1,WRK(KRYC),1)
C
  500 CONTINUE
  520 CONTINUE
C
C
C     Calculate Rxc and Ryc contributions to SCVECS(NVAR,NCSIM)
C     =========================================================
C
C     Rxc contribution is gradient-like,
C     Ryc contribution is "lintrn"-like:
C
C     scvecs(i) = scvecs(i) + [ d(Rxc) / di ]
C                           + sum(j) [ d2(Ryc) / di dj ] * bcvecs(j)
C
C
C     ... CSF part of sigma vectors
C
      JRXC   = KRXC
      JRXCAC = KRXCAC
      DO 600 ICSIM = 1,NCSIM
         CALL GETAC2(WRK(JRXC),WRK(JRXCAC))
         TXC    = SOLELM(DV,WRK(JRXCAC),WRK(JRXC),TXCAC)
         WRK(KTXCAC-1+ICSIM) = TXCAC
         JRXC   = JRXC   + NNORBX
         JRXCAC = JRXCAC + NNASHX
  600 CONTINUE
C
      CALL GETAC2(WRK(KRYC),WRK(KRYCAC))
C
      CALL SOLSC(NCSIM,0,BCVEC,CREF,SVEC,WRK(KRXCAC),WRK(KRYCAC),
     *           WRK(KTXCAC),TYCAC,INDXCI,WRK(KW30),LW30)
C     CALL SOLSC(NCSIM,NOSIM,BCVECS,CREF,SVECS,
C    *           RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
C
C
C     ... orbital part of sigma vector(s)
C
      IF (NWOPPT .GT. 0) THEN
         MWOPH  = NWOPH
         NWOPH  = NWOPPT
C        ... tell SOLGO only to use the NWOPPT first JWOP entries
         JSVECO = 1 + NCONST
         JRXC   = KRXC
         DO 800 ICSIM = 1,NCSIM
            CALL SOLGO(D2,DV,WRK(JRXC),SVEC(JSVECO,ICSIM))
            IF (IPRSOL.GT.101) THEN
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               WRITE(LUPRI,*)' TXC CONTRIBUTION'
               CALL OUTPUT(SVEC(JSVECO,ICSIM),1,NWOPPT,1,1,
     *                                        NWOPPT,1,1,LUPRI)
            END IF
            CALL SOLGO(D0,DTV(1,ICSIM),WRK(KRYC),SVEC(JSVECO,ICSIM))
            IF (IPRSOL.GT.101) THEN
               WRITE(LUPRI,*)
     *         ' orbital part of LINEAR TRANSFORMED CONF VEC No',ICSIM
               WRITE(LUPRI,*)' + TG CONTRIBUTION'
               CALL OUTPUT(SVEC(JSVECO,ICSIM),1,NWOPPT,1,1,
     *                                        NWOPPT,1,1,LUPRI)
            END IF
            JRXC   = JRXC   + NNORBX
  800    CONTINUE
         NWOPH  = MWOPH
      END IF
C
      IF (IPRSOL .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLLNC - svec(ci,1) on exit'
         DO I = 1,NCONST
            IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *      ' conf #',I,SVEC(I,1)
         END DO
         WRITE (LUPRI,'(/A)') ' >>> SOLLNC - svec(orb) on exit'
         CALL OUTPUT(SVEC(1 + NCONST,1),1,NWOPPT,1,NCSIM,
     *               NSVEC,NCSIM,1,LUPRI)
      END IF
C
      CALL QEXIT('SOLLNC')
      RETURN
C     end of sollnc.
      END
C  /* Deck sollno */
      SUBROUTINE SOLLNO(NOSIM,BOVECS,CREF,CMO,INDXCI,
     *                  FLVEC,FLINR,TLMSI,
     *                  DV,SVEC,NSVEC,ISYRLM,WRK,LFREE)
C
C  Copyright 24-Jan-1987 Hans Jorgen Aa. Jensen
C  Revised for LSYMPT .ne. 1 26-Aug-1998 hjaaj
C
C  Purpose:  Calculate MCSCF Hessian contribution from a
C            surrounding medium to an orbital trial vector.
C            Cavity radius = Rsol and dielectric constant = EPsol.
C
C  NSVEC          may be NVAR or NWOPT, dependent on LINTRN
C
#include <implicit.h>
      DIMENSION BOVECS(NWOPPT,*), CREF(*), CMO(*)
      DIMENSION INDXCI(*),        DV(*),   ISYRLM(*)
      DIMENSION FLVEC(NLMSOL), FLINR(NLMSOL), TLMSI(NLMSOL)
      DIMENSION SVEC(NSVEC,*),    WRK(*)
C
#include <iratdef.h>
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 )
#include <thrzer.h>
#include <dummy.h>
C
C
C  Used from common blocks:
C    INFINP : NLMSOL, LSOLMX, INERSF
C    INFORB : NNASHX, NNORBX, NNBASX, etc.
C    INFVAR : JWOP
C    INFLIN : NWOPPT, NVARPT, NCONST, NCONRF
C    INFTAP : LUSOL
C    INFPRI : IPRSOL
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inflin.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL FULHES
C
C     Statement functions
C
      TN(   LM) = WRK((KTNLM -1)+LM)
C
      CALL QENTER('SOLLNO')
      IF (JWOPSY .NE. LSYMPT .OR. NWOPPT .NE. NWOPT) THEN
         WRITE (LUPRI,*) 'SOLLNO ERROR: JWOPSY .ne. LSYMPT .or.'//
     *      ' NWOPPT .ne. NWOPT'
         WRITE (LUPRI,*) 'LSYMPT,JWOPSY:',LSYMPT,JWOPSY
         WRITE (LUPRI,*) 'NWOPPT,NWOPT :',NWOPPT,NWOPT
         CALL QUIT('SOLLNO ERROR,lsympt.ne.jwopsy .or. nwoppt.ne.nwopt')
      END IF
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
         WRITE (LUPRI,'(//A)') ' >>> TEST OUTPUT FROM SOLLNO >>>'
      END IF
      IF (IPRSOL .GE. 140) THEN
         IF (FULHES) THEN
            WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(ci,1) on entry'
            DO 30 I = 1,NCONST
               IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *              ' conf #',I,SVEC(I,1)
 30         CONTINUE
         END IF
         WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(orb) on entry'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *               NSVEC,NOSIM,1,LUPRI)
      END IF
C
C
C     Core allocation
C
      KRXYO  = 1
      KTEXY  = KRXYO  + NOSIM*NNORBX
      KUBO   = KTEXY  + NOSIM
      KW10   = KUBO   + NOSIM*N2ORBX
C
      KTNLM  = KW10
      KUCMO  = KTNLM  + NLMSOL
      KRLMAC = KUCMO  + NORBT*NBAST
      KRLM   = KRLMAC + NNASHX
      KW20   = KRLM   + NNORBX
C     2.1 read rlmao in ao basis and transform to rlm in mo basis
      KRLMAO = KW20
      KURLM  = KRLMAO + NNBASX
      KURXLM = KURLM  + N2ORBX
      KW21   = KURXLM + N2ORBX
C
      KRXLM  = KW20
      KW22   = KRXLM  + NNORBX
C     3.0 SOLSC
      IF (FULHES) THEN
         KRXYOA = KW10
         KTXYOA = KRXYOA + NOSIM*NNASHX
         KOVLP  = KTXYOA + NOSIM
         KW30   = KOVLP  + NOSIM
         LW30   = LFREE  + 1 - KW30
      ELSE
         KW30   = 0
         LW30   = 0
      END IF
C
      KNEED = MAX(KW21,KW22,KW30)
      IF (KNEED .GT. LFREE) CALL ERRWRK('SOLLNO',KNEED,LFREE)
C
C     Read nuclear contributions to ERLM (TN(l,m))
C     and position for R(l,m).
C
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,lupri)
      READ (LUSOL)
      CALL READT(LUSOL,NLMSOL,WRK(KTNLM))
C
C
C
C     Construct effective operators Rxo and Ryo
C     =========================================
C
C     Zero Rxo and Ryo storage
C
      CALL DZERO(WRK(KRXYO), NOSIM * NNORBX )
      CALL DZERO(WRK(KTEXY), NOSIM )
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
              CALL OUTPUT(WRK(JUBO),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
  210    CONTINUE
      END IF
 2110 FORMAT (/,' Orbital trial vector unpacked to matrix form (no.',
     *        I3,' of',I3,')')
C
C
C     Loop over l,m expansion.
C
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
      DO 500 M = -L,L
         LM = LM + 1
         IF (ISYRLM(L+M+1) .NE. LSYMPT .AND. ISYRLM(L+M+1) .NE. 1) THEN
            IF (ISYRLM(L+M+1) .GT. 0) THEN
               READ (LUSOL)
            ELSE IF (ISYRLM(L+M+1) .EQ. -LSYMPT
     &          .OR. ISYRLM(L+M+1) .EQ. -1) THEN
             WRITE(LUPRI,*)'SOLLNO ERROR, needed RLM not on LUSOL file!'
             WRITE(LUPRI,*)'- You must not specify ".NOT ALL RLM" under'
             WRITE(LUPRI,*)'  "*ONEINT" when running ABACUS or RESPONSE'
             WRITE(LUPRI,*)'  with non-totally symmetric perturbations.'
             WRITE(LUPRI,*)'L,M,ISYRLM(l,m), LSYMPT',
     &                    L,M,ISYRLM(L+M+1),LSYMPT
             CALL QUIT('SOLLNO ERROR, needed RLM not on file!')
            END IF
            GO TO 500
         END IF
C
C        Read R(l,m) in ao basis and transform to mo basis.
C        Extract active-active block in RLMAC = WRK(KRLMAC).
C        Electronic contribution TE(l,m).
C
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KRLM),WRK(KUCMO),WRK(KW21),
     *             NBAST,NORBT)
         IF (NASHT .GT. 0) THEN
            CALL GETAC2(WRK(KRLM),WRK(KRLMAC))
         END IF
         TELM  = SOLELM(DV,WRK(KRLMAC),WRK(KRLM),TELMAC)
         TLM   = TN(LM) - TELM
C
         IF (IPRSOL .GE. 40) THEN
            WRITE (LUPRI,'(/A,I5,4(/A,F15.10))')
     *         ' LM       ',LM,
     *         ' TELM     ',TELM,
     *         ' TELMAC   ',TELMAC,
     *         ' TLM      ',TLM,
     *         ' FLVEC(LM)',FLVEC(LM)
         END IF
         IF (IPRSOL .GE. 70) THEN
            WRITE (LUPRI,'(/A)') ' Rlm_mo matrix:'
            CALL OUTPAK(WRK(KRLM),  NORBT,1,LUPRI)
            IF (NASHT .GT. 0) THEN
               WRITE (LUPRI,'(/A)') ' Rlm_ac matrix:'
               CALL OUTPAK(WRK(KRLMAC),NASHT,1,LUPRI)
            END IF
         END IF
C
         CALL DSPTSI(NORBT,WRK(KRLM),WRK(KURLM))
         FACY  = - D2 * FLVEC(LM) * TLM
         IF (INERSF) THEN
            FACY  = FACY - FLINR(LM) * D2 * TLMSI(LM)
         END IF
         DO 200 IOSIM = 1,NOSIM
            JRXYO = KRXYO + (IOSIM-1)*NNORBX
C
C           Calculate one-index transformed integrals
C
            CALL DZERO(WRK(KURXLM),N2ORBX)
            JUBO = KUBO + (IOSIM-1)*N2ORBX
            CALL TR1UH1(WRK(JUBO),WRK(KURLM),WRK(KURXLM),1)
            CALL DGETSP(NORBT,WRK(KURXLM),WRK(KRXLM))
C           MAERKE: what if LSYMPT .ne. 1 ???
C           900720: I think OK, i.e. that only RLM of sym 1
C           such that RXLM is of sym LSYMPT is needed.
C           CALL TR1UH1(UBOVEC,H1,H1X,IH1SYM)
C
            IF (NASHT .GT. 0) CALL GETAC2(WRK(KRXLM),WRK(KRLMAC))
            IF (IPRSOL .GE. 60) THEN
               WRITE (LUPRI,'(/A,I5)') ' Rlm_X_mo matrix, IOSIM =',IOSIM
               CALL OUTPAK(WRK(KRXLM),NORBT,1,LUPRI)
               IF (NASHT .GT. 0) THEN
                  WRITE (LUPRI,'(/A)') ' Rlm_X_ac matrix:'
                  CALL OUTPAK(WRK(KRLMAC),NASHT,1,LUPRI)
               END IF
            END IF
            FACX  = D2 * FLVEC(LM)
     *                 * SOLELM(DV,WRK(KRLMAC),WRK(KRXLM),TAC)
            WRK(KTEXY-1+IOSIM) = WRK(KTEXY-1+IOSIM) + FACX * TELM
            IF (IPRSOL .GE. 40)
     *         WRITE (LUPRI,'(/A,2I5,2F15.10)')
     *         ' LM, IOSIM, FACX, FACY :',LM,IOSIM,FACX,FACY
            IF (ABS(FACX) .GT. THRZER)
     *         CALL DAXPY(NNORBX,FACX,WRK(KRLM), 1,WRK(JRXYO),1)
            IF (ABS(FACY) .GT. THRZER)
     *         CALL DAXPY(NNORBX,FACY,WRK(KRXLM),1,WRK(JRXYO),1)
  200    CONTINUE
C
  500 CONTINUE
  520 CONTINUE
C
C
      IF (IPRSOL .GE. 50) THEN
         DO 600 IOSIM = 1,NOSIM
            JRXYO = KRXYO + (IOSIM-1)*NNORBX
            WRITE (LUPRI,'(/A,I3,A,I3)')
     *         ' >>> SOLLNO - (Rxo + Ryo) matrix no.',IOSIM,' of',NOSIM
            CALL OUTPAK(WRK(JRXYO),NORBT,1,LUPRI)
  600    CONTINUE
      END IF
C
C
C
C     Calculate Rxo and Ryo contributions to SVEC(NSVEC,NOSIM)
C     =========================================================
C
C     Rxo and Ryo contributions are gradient-like:
C
C     svec(i) = svec(i) + [ d(Rxo) / di ] + [ d(Ryo) / di ]
C             = svec(i) + [ d(Rxo + Ryo) / di ]
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
            JRXYO  = KRXYO  + (IOSIM-1)*NNORBX
            JRXYOA = KRXYOA + (IOSIM-1)*NNASHX
            IF (LSYMRF .EQ. LSYMST) THEN
               WRK(KOVLP-1+IOSIM) = DDOT(NCONRF,CREF,1,SVEC(1,IOSIM),1)
            ELSE
               WRK(KOVLP-1+IOSIM) = D0
            END IF
            CALL GETAC2(WRK(JRXYO),WRK(JRXYOA))
            TXYO   = SOLELM(DV,WRK(JRXYOA),WRK(JRXYO),TXYOAC)
            IF (IPRSOL .GE. 40) WRITE (LUPRI,'(/A,I5,3F15.10)')
     *         ' IOSIM, C_OVLP, TXYOAC, TXYO :',
     *         IOSIM,WRK(KOVLP-1+IOSIM),TXYOAC,TXYO
            WRK(KTXYOA-1+IOSIM) = TXYOAC
  700    CONTINUE
         CALL SOLSC(0,NOSIM,DUMMY,CREF,SVEC,WRK(KRXYOA),DUMMY,
     *              WRK(KTXYOA),DUMMY,INDXCI,WRK(KW30),LW30)
C        CALL SOLSC(NCSIM,NOSIM,BCVECS,CREF,SVECS,
C    *              RXAC,RYAC,TRXAC,TRYAC,INDXCI,WRK,LWRK)
      END IF
C
      MWOPH  = NWOPH
      NWOPH  = NWOPPT
C    ... tell SOLGO only to use the NWOPPT first JWOP entries
      DO 800 IOSIM = 1,NOSIM
         JRXYO  = KRXYO  + (IOSIM-1)*NNORBX
         CALL SOLGO(D2,DV,WRK(JRXYO),SVEC(JSOVEC,IOSIM))
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
            WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(ci,1) on exit'
            DO 930 I = 1,NCONST
               IF (SVEC(I,1) .NE. D0) WRITE (LUPRI,'(A,I10,F15.10)')
     *         ' conf #',I,SVEC(I,1)
  930       CONTINUE
         END IF
      END IF
      IF (IPRSOL .GE. 140) THEN
         WRITE (LUPRI,'(/A)') ' >>> SOLLNO - svec(orb) on exit'
         CALL OUTPUT(SVEC(JSOVEC,1),1,NWOPPT,1,NOSIM,
     *               NSVEC,NOSIM,1,LUPRI)
      END IF
C
      CALL QEXIT('SOLLNO')
      RETURN
C     ... end of sollno.
      END
C  /* Deck solsc */
      SUBROUTINE SOLSC(NCSIM,NOSIM,BCVECS,CREF,SVECS,
     *                 RXAC,RYAC,TRXAC,TRYAC,XNDXCI,WRK,LWRK)
C
C   Copyright 5-May-1987 Hans Joergen Aa. Jensen
C   Rewritten for determinant CI 19-Jul-1990 hjaaj
C
C   Purpose:
C      Perform CI transformation of a one-electron operator,
C      as the solvent operator for "cavity in dielectric" model.
C
C
C  Input:   ordinary and 1-index
C           transformed inactive solvent-matrices (RXAC and RYAC),
C           tral and reference CI-vectors (BCVECS and CREF)
C
C  Output:  New CI-part of sigma vector (SVECS)
C
C Scratch:  WRK
C
#include <implicit.h>
      DIMENSION BCVECS(NCONST,*),  SVECS(NVARPT,*)
      DIMENSION CREF(*),  RXAC(NNASHX,*),RYAC(NNASHX)
      DIMENSION TRXAC(*), XNDXCI(*), WRK(LWRK)
C
      PARAMETER ( D2 = 2.0D0 )
C     ... note: factors of 2 (D2) below takes care of the "2" in
C         2 ( <j / R / B> - T bj )
C
C Used from common blocks:
C   INFORB : NNASHX,N2ASHX,NASHT
C   INFLIN : LSYMRF,LSYMST,LSYMPT,NCONRF,NCONST
C
#include <priunit.h>
#include <inforb.h>
#include <inflin.h>
#include <infpri.h>
C
      CALL QENTER('SOLSC ')
C
C     ... input check
C
      IF (NCSIM .GT. 0 .AND. NOSIM .GT. 0) THEN
         WRITE (LUERR,'(//A/A,2I5)')
     *   ' >>> SOLSC ERROR, both NCSIM .gt. 0 and NOSIM .gt. 0',
     *   '     NCSIM, NOSIM =',NCSIM,NOSIM
         CALL QTRACE(LUERR)
         CALL QUIT('SOLSC ERROR, both NCSIM .gt. 0 and NOSIM .gt. 0')
      END IF
C
      IF (NCSIM .GT. 0) THEN
         KURYAC = 1
         KURXAC = KURYAC + N2ASHX
         KW1    = KURXAC + N2ASHX
         LW1    = LWRK + 1 - KW1
         CALL DSPTSI(NASHT,RYAC,WRK(KURYAC))
         CALL DSCAL(N2ASHX,D2,WRK(KURYAC),1)
         DO 18 ICSIM = 1,NCSIM
            CALL DSPTSI(NASHT,RXAC(1,ICSIM),WRK(KURXAC))
            CALL DSCAL(N2ASHX,D2,WRK(KURXAC),1)
C           First SVECS(NA) = SVECS(NA) + < NA | 2*RX | CREF >
            IF (IPRSOL.GT.100) THEN
               WRITE(LUPRI,*)ICSIM,' RXAC MATRIX SCALED WITH 2'
               CALL OUTPUT(WRK(KURXAC),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
            END IF
            CALL CISIGD(LSYMRF,LSYMST,NCONRF,NCONST,CREF,SVECS(1,ICSIM),
     &                  WRK(KURXAC),DUMMY,
     &                  .TRUE.,.FALSE.,XNDXCI,0,0,WRK(KW1),LW1)
            IF (IPRSOL.GT.125) THEN
               WRITE(LUPRI,*)' LINEAR TRANS WITH RXAC ADDED '
               CALL OUTPUT(SVECS(1,ICSIM),1,NCONST,1,1,NCONST,1,1,LUPRI)
            END IF
C           then  SVECS(NA) = SVECS(NA) + < NA | 2*RY | BCVECS >
            IF (IPRSOL.GT.100) THEN
               WRITE(LUPRI,*)' RYAC MATRIX SCALED WITH 2'
               CALL OUTPUT(WRK(KURYAC),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
            END IF
            CALL CISIGD(LSYMST,LSYMST,NCONST,NCONST,
     &                  BCVECS(1,ICSIM),SVECS(1,ICSIM),
     &                  WRK(KURYAC),DUMMY,
     &                  .TRUE.,.FALSE.,XNDXCI,0,0,WRK(KW1),LW1)
            IF (IPRSOL.GT.125) THEN
               WRITE(LUPRI,*)' LINEAR TRANS WITH RYAC ADDED '
               CALL OUTPUT(SVECS(1,ICSIM),1,NCONST,1,1,NCONST,1,1,LUPRI)
            END IF
            IF (LSYMRF .EQ. LSYMST) THEN
               DO 200 NA = 1,NCONST
                  SVECS(NA,ICSIM) = SVECS(NA,ICSIM)
     *                            - D2 * TRXAC(ICSIM) * CREF(NA)
     *                            - D2 * TRYAC        * BCVECS(NA,ICSIM)
  200          CONTINUE
            ELSE
               FAC = - D2 * TRYAC
               CALL DAXPY(NCONST,FAC,BCVECS(1,ICSIM),1,SVECS(1,ICSIM),1)
            END IF
            IF (IPRSOL.GT.125) THEN
               WRITE(LUPRI,*)' -2*TRXAC*CREF -2*TRYAC*BCVECS ADDED '
               CALL OUTPUT(SVECS(1,ICSIM),1,NCONST,1,1,NCONST,1,1,LUPRI)
            END IF
   18    CONTINUE
      END IF
      IF (NOSIM .GT. 0) THEN
         KURXAC = 1
         KW1    = KURXAC + N2ASHX
         LW1    = LWRK + 1 - KW1
         DO 22 IOSIM = 1,NOSIM
            CALL DSPTSI(NASHT,RXAC(1,IOSIM),WRK(KURXAC))
            CALL DSCAL(N2ASHX,D2,WRK(KURXAC),1)
C           now SVECS(NA) = SVECS(NA) + <NA | 2*RXY | CREF >
            IF (IPRSOL.GT.100) THEN
               WRITE(LUPRI,*)' RXAC MATRIX SCALED WITH 2, iosim',IOSIM
               CALL OUTPUT(WRK(KURXAC),1,NASHT,1,NASHT,
     *                                 NASHT,NASHT,1,LUPRI)
            END IF
            CALL CISIGD(LSYMRF,LSYMST,NCONRF,NCONST,CREF,SVECS(1,IOSIM),
     &               WRK(KURXAC),DUMMY,
     &               .TRUE.,.FALSE.,XNDXCI,0,0,WRK(KW1),LW1)
            IF (IPRSOL.GT.125) THEN
               WRITE(LUPRI,*)' LINEAR TRANS WITH RXAC ADDED '
               CALL OUTPUT(SVECS(1,IOSIM),1,NCONST,1,1,NCONST,1,1,LUPRI)
            END IF
C           Remove CREF component of SVECS (inactive contribution
C           not needed because of this projection)
            IF (LSYMRF .EQ. LSYMST) THEN
               DO 260 NA = 1,NCONRF
                  SVECS(NA,IOSIM) = SVECS(NA,IOSIM)
     *                            - D2 * TRXAC(IOSIM) * CREF(NA)
  260          CONTINUE
               IF (IPRSOL.GT.125) THEN
                  WRITE(LUPRI,*)' -2*TRXAC*CREF  ADDED, iosim',IOSIM
                  CALL OUTPUT(SVECS(1,IOSIM),1,NCONRF,1,1,
     *                        NCONRF,1,1,LUPRI)
               END IF
            END IF
   22    CONTINUE
      END IF
C
C     CALL CISIGD(ICSYM,IHCSYM,NCDET,NHCDET, C,HC, FCAC,H2AC,
C    &            NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK,LFREE)
C
      CALL QEXIT('SOLSC ')
      RETURN
      END
C --- end of sirsol.F ---
