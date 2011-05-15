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
#ifdef OLD_REV_LOG
!===========================================================================
!
!Revision 1.5  2000/05/24 12:43:25  hjj
! 1) attempted to get WLKREJ true working again (not tested)
! 2) new parameter lists to SOLFCK and WR_SIRIFC
! 3) write CREF in determinants on SIRIFC in WR_SIRIFC
!
!Revision 1.4  2000/05/01 14:07:33  hjj
! removed obsolete KAVER,KSYM(4) in WRITE (LUSIFC)
!
!revision 1.3 2000/05/01 07:21:48 hjj
! polishing of warnings so that 'grep WARNING' gets more relevant information.
!----------------------------
!revision 1.2 2000/04/14 10:10:52 hjj
!\hjaaj jan 2000:
! also screening for SIRCNO in OPTST
!\hjaaj dec 1999:
! 1) implemented dynamic Fock matrix screening for QCHF and DIRFCK
! 2) polishing
!===========================================================================
#endif
C  /* Deck siropt */
      SUBROUTINE SIROPT(WRK,LFREE,ICONV)
C
C 31-Oct-1983-V1/6-May-1984-V2 -- Hans Jorgen Aa. Jensen
C Revisions
C   6-Aug-1984/16-May-1984 - hjaaj
C  14-May-1985 - hjaaj (update     -- sirupd)
C  10-Nov-1985 - hjaaj (absorption -- sirorb)
C  14-May-1986 - hjaaj (nr eq's    -- sirnr)
C  Last revision 15-Mar-1991 - hjaaj
C
C Purpose: Control MCSCF optimization.
C
C MOTECC-90: The purpose of this module, SIROPT, and the algorithms used
C            are described in Chapter 8 Section C.0 of MOTECC-90
C            "Implementation"
C
C Output:
C  ICONV  .eq.0  MCSCF not converged.
C         .gt.0  MCSCF converged.
C         -----  ICONV can be used for conditional call of
C                subsequent steps requiring a converged w.f.
C
#ifdef VAR_MPI
      use par_mcci_io
#endif

#include <implicit.h>

C
      DIMENSION WRK(LFREE)
C
#include <thrldp.h>
      PARAMETER (MAXBCK = 5)
      PARAMETER (FAKABS = 0.40D0, STPABS = 0.20D0)
      PARAMETER (STPLCL = 0.15D0, SHFLCL = 1.D-3)
      PARAMETER (D0 = 0.D0, D1 = 1.D0, D2 = 2.0D0)
#include <iratdef.h>
#include <dummy.h>
C
C used from common blocks:
C  INFINP: LSYM,FLAG(*),...,ITRLVL,ITRFIN,THRPWF,THRCGR,MAXUIT,
C          INERSI,...
C  INFVAR: NCONF,NWOPT,NVAR,NVARH
C  INFORB: NCMOT,NNASHX,...
C  INFOPT: EMCSCF,...
C  INFDIM: LCINDX,MAXRL,?
C  INFTAP: LUIT1,?
C  CBIHR2: IFTHRS,USRSCR
C  CBIREA: LMULBS
C  R12INT: R12CAL
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <pcmlog.h>
#include <infvar.h>
#include <inforb.h>
#include <infopt.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
#include <itinfo.h>
#include <gnrinf.h>
#include <cbihr2.h>
#include <cbirea.h>
#include <r12int.h>
C------------------------
#include <dfterg.h>
#include <dftcom.h>
C ---
C
C *** local:
      CHARACTER*8 STAR8, RTNLBL(2)
      LOGICAL ABSORP, UPDATE, SOLVNT, GEOWLK, RHFCLC
      LOGICAL CHKSTP, DONR, NRALW,   LOCAL
C *** data:
      DATA STAR8/'********'/
C
      CALL QENTER('SIROPT')
C
      IF (IPRSTAT .GT. 0) CALL TIMOPT('START',LUSTAT,WRK,LFREE)
C
      IF (NASHT .GT. 1) THEN
         WRITE (LUW4,960)
         IF (LUPRI .NE. LUW4) WRITE (LUPRI,960)
      ELSE
         WRITE (LUW4,961)
         IF (LUPRI .NE. LUW4) WRITE (LUPRI,961)
      END IF
  960 FORMAT(//T9,     'SIRIUS MCSCF optimization (SIROPT)'
     *       /' ================================================'/)
  961 FORMAT(//T9,     'SIRIUS QC-HF optimization (SIROPT)'
     *       /' ================================================'/)
C
      TOPT  = SECOND()
      TIMST = TOPT
C
C *** Initializations
C     convergence flag ICONV to "not converged"
C     UPDATE , if use UPDATE for MC optimization
C
      ICONV  = 0
      RHFCLC = (NASHT .LE. 1 .OR. HSROHF)
      SOLVNT = FLAG(16)
      ABSORP = ( FLAG(51) .OR. FLAG(52) .OR. FLAG(53) )
#if defined (VAR_NOEXCABS)
      IF ( ABSORP .AND. ISTATE .GT. 1 ) THEN
C        WRITE (LUW4,'(//A/A/)') ' *** Orbital absorption requested,',
C    *      '     orbital absorption is not currently defined for'//
C    *      ' excited states and request is ignored.'
C        ABSORP   = .FALSE.
C        FLAG(51) = .FALSE.
C        FLAG(52) = .FALSE.
C        FLAG(53) = .FALSE.
         NWARN = NWARN + 1
         WRITE (LUW4,'(//A/)') ' *** WARNING *** ABSORPTION FOR'//
     *      ' AN EXCITED STATE IS EXPERIMENTAL'
      END IF
#endif
      IF (ABSORP .AND. MAXAPM .LE. 0) THEN
         WRITE (LUW4,'(//A/A/)') ' *** Orbital absorption requested,',
     *      '     request is ignored because zero iterations specified.'
         ABSORP   = .FALSE.
         FLAG(51) = .FALSE.
         FLAG(52) = .FALSE.
         FLAG(53) = .FALSE.
      END IF
CHJAaJ Nov 09: KTRLVL=0 does not work in TR1H2M in current Dalton
CHJ   IF (FLAG(51) .OR. FLAG(52)) THEN
CHJ      KTRLVL = 0
CHJ   ELSE
         KTRLVL = ITRLVL
CHJ   END IF
#if defined (VAR_NOITRABS)
C --- 851216 -- EXPLICIT ORBITAL HESSIAN DIAGONAL IS WRONG
C               IF KTRLVL = 0 !
      KTRLVL = ITRLVL
#endif
C     If FLAG(39) always NR (e.g. for core hole calculations /860515-hj)
C     921028-hjaaj: FLAG(39) copied to local variable, because
C     it may be reset to true in local region
      NRALW  = FLAG(39)
      DONR   = NRALW
      UPDATE = FLAG(19)
      ITRLSV = ITRLVL
      LTRLVL = -999
      NWOPS  = NWOPH
      NVARS  = NVARH
C
C     dec 99 hjaaj:
C     1) save initial screening threshold for Direct Fock matrix construction
C     2) set screening threshold to accuracy for converged gradient for first
C        iteration (only used if DIRFCK true)
C
      IFTHSV = IFTHRS
      GNORM(3) = THRGRD
C
C *** Core allocation:
C
C.. 1.0 .. (OPTST, READMO, GETCIX, GRAD and SIRNEO)
      KCINDX = 1
      KCMO   = KCINDX + LCINDX
C.. 1.1.1 .. (GETCIX, TRACTL, READMO)
      KWRK10 = KCMO   + NCMOT
      LWRK10 = LFREE  + 1 - KWRK10
C.. 1.2 .. (GRAD, SIRORB and SIRNEO)
      IF (DFT_SPINDNS) THEN
         NDV = 2
      ELSE
         NDV = 1
      END IF
      KGTOT  = KWRK10
      KDV    = KGTOT  + NVARH
      KPV    = KDV    + NDV*NNASHX
      KFC    = KPV    + NNASHX*NNASHX
      KFV    = KFC    + NNORBT
!HJJ - FRAN KFS added for the SR_DFT Spin Dens. contribution
      KFS    = KFV    + NNORBT
      IF (DFT_SPINDNS) THEN
         KSRHYB = KFS + NNORBT
      ELSE
         KSRHYB = KFS
      END IF
!HJJ - FRAN END
      IF (SRHYBR .AND. MCTYPE.NE.0) THEN
         KFQ = KSRHYB + NNORBT
      ELSE
         KFQ = KSRHYB
      END IF
      KFCAC  = KFQ    + NASHT*NORBT
      KCREF  = KFCAC  + NNASHX
      KH2AC  = KCREF  + NCONF
      KWRK12 = KH2AC  + NNASHX*NNASHX
C.. 1.2.1 .. (GRAD,SOLGRD,WR_SIRIFC)
      KERLM  = KWRK12
      IF (SOLVNT) THEN
         KWRKG  = KERLM + 2*NLMSOL
      ELSE
         KWRKG  = KWRK12
      END IF
      LWRKG  = LFREE  - (KWRKG-1)
      IF (LWRKG.LT.0) CALL ERRWRK('SIROPT.GRAD',-KWRKG,LFREE)
C.. 1.2.2 .. (SIRORB)
      KWOCTL = KFCAC
      LWOCTL = LFREE  - KWOCTL
      IF (ABSORP .AND. LWOCTL .LE. 0)
     *   CALL ERRWRK('SIROPT.SIRORB',-KWOCTL,LFREE)
      IF (.NOT.FLAG(20)) THEN
C.. 1.2.3 .. (SIRSTP, SIRNEO)
         KIBNDX = KWRK12
         KWNEO  = KIBNDX + (MAXRL + MAXRTS*3)/IRAT + 1
         LWRK2  = KWNEO  + 3*NVAR
C       (KWNEO space for vectors as above;
C        ?????? estimate of space for DTV,PTV,FXC,...; MAERKE
C        3*NVAR space for one Bvec, one Svec, and diagonal)
         IF (LWRK2 .GT. LFREE) CALL ERRWRK('SIROPT.SIRNEO',-LWRK2,LFREE)
         LWNEO  = LFREE + 1 - KWNEO
      END IF
C.. 1.3 .. (SIRNEO)
C.. 1.4 .. (SIRUPD)
      KREFU  = KWRK10
      KGRDU  = KREFU + NVAR
      KWU    = KGRDU + NVAR
      LWU    = LFREE - KWU
C..
      NCONF4 = MAX(4,NCONF)
C
C     Get CI information if active orbitals
C
      IF (.NOT.RHFCLC)
     &   CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KWRK10),LWRK10,0)
C
C ****************************
C *** GET STARTING VECTORS ***
C ****************************
C
C   OPTST constructs LROOTS starting vectors for
C   the very first macroiteration.
C   The C-vectors are constructed in OPTST or from the result of
C   another run (another geometry). Later the resulting vectors
C   from the previous macroiteration is read in.
C   ITMAC and EMCOLD are initialized before CALL OPTST
C   so they can be defined by OPTST when restart.
C
C   ITRLVL is temporarily set to KTRLVL, so we don't perform larger
C   integral transformation in OPTST than needed.
C
      ITMAC  = 0
      EMCOLD = D0
      DEPRED = D0
      STPLEN = D0
      JTRLVL = ITRLVL
      ITRLVL = KTRLVL
C
C jan00 hjaaj: set screening threshold for SIRCNO in OPTST for direct
C    SCF (canonical orbitals are not needed with high precision).
C
      IF (DIRFCK .AND. .NOT. USRSCR) THEN
         IFTHRS  = MIN( IFTHSV, 9 )
         IF (IPRI6 .GE. 0) WRITE (LUPRI,'(/A,I5)')
     &      ' Fock matrix screening setting for SIRCNO: IFTHRS =',
     &      IFTHRS
      END IF
      CALL OPTST(WRK(KCINDX),WRK(KWRK10),LWRK10)
C     CALL OPTST(INDXCI,WRK,LFREE)
      ITRLVL = JTRLVL
C
C     set parameters
C       ABSORP may be reset by OPTST, if restart
C       DONR   may be reset by OPTST, if restart
C       CHKSTP false for first macro iteration, unless restart
C
      ABSORP = ( FLAG(51) .OR. FLAG(52) .OR. FLAG(53) )
      IF (.NOT. FLAG(54)) ABSORP = ABSORP .AND. (MAXABS .GT. 0)
      DONR   = FLAG(39)
      GEOWLK = (ICI0 .EQ. 6)
      IF ( GEOWLK ) THEN
         JWSTEP = 1
         CHKSTP = .FALSE.
         ISTEP  = 0
         EMCGEO = EMCOLD
         DEPGEO = DEPRED
         THRGEO = MIN(1.D-2,0.1D0*SQRT( ABS(DEPGEO) ))
      ELSE
         JWSTEP = 0
         IF (DEPRED .NE. D0) THEN
            CHKSTP = .TRUE.
            ISTEP  = -1
         ELSE
            CHKSTP = .FALSE.
            ISTEP  = 0
         END IF
      END IF
C
C-841209-hj: force beta positive
      BETA   = ABS(BETA)
      IF (BETA   .EQ. D0)     BETA   = D1
      IF (BETMIN .LE. D0 .OR.
     *    BETMIN .GT. BETMAX) BETMIN = D1
      IF (BETMAX .LT. D1)     BETMAX = D1
C
C
      TIMST = SECOND() - TIMST
      IF (IPRI6 .GT. 1) WRITE (LUPRI,'(/A,F15.2,A)')
     *   ' >>> Time used in SIROPT startup :',TIMST,' seconds.'
C
C ***********************************************
C *** START LOOP OVER MACRO ITERATIONS **********
C ***********************************************
C
      ITMACN = 0
      ITMICT = 0
      ITBCK  = 0
      ITABS  = 0
C*****NON-STANDARD LOOP (MACRO ITERATIONS)
  100 CONTINUE
      IF (.NOT.ABSORP) ITABS = -1
C
      TIMMAC = SECOND()
      IF (ITBCK .GT. 0) THEN
         WRITE (LUW4,975) ITMAC,ITBCK
         IF (LUPRI.NE.LUW4) WRITE (LUPRI,975) ITMAC,ITBCK
      ELSE IF (ITABS .GT. 0) THEN
         WRITE (LUW4,977) ITMAC
         IF (LUPRI.NE.LUW4) WRITE (LUPRI,977) ITMAC
      ELSE
         ITMAC  = ITMAC  + 1
         ITMACN = ITMACN + 1
         WRITE (LUW4,970) ITMAC
         IF (LUPRI.NE.LUW4) WRITE (LUPRI,970) ITMAC
      END IF
  970 FORMAT(//' <<< MACRO ITERATION',I3,' >>>'
     *        /' --------------------------')
  975 FORMAT(//' <<< MACRO ITERATION',I3,' BACKSTEP NO.',I3,' >>>')
  977 FORMAT(//' <<< MACRO ITERATION',I3,' AFTER ABSORPTION >>>')
C
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
      CALL DZERO(DINFO,LDINFO)
      CALL IZERO(IINFO,LIINFO)
      SHFLVL = D0
C
C     First retrieve orbitals :
C
      CALL READMO (WRK(KCMO),9,WRK(KWRK10),LWRK10)
C
C *** Transform two-electron integrals
C
      TIMITR = SECOND()
      IF (.NOT.FLAG(14)) THEN
         JTRLVL = ITRLVL
         IF (ITABS .EQ. 0) JTRLVL = KTRLVL
         IF (FLAG(20) .OR. ITMACN .GE. MAXMAC) JTRLVL = 0
C        ... stop after gradient: only 1. order transform. needed
         CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KWRK10),LWRK10)
C        CALL TRACTL(ITRLVL,CMO,WRK,LFREE)
         LTRLVL = JTRLVL
      ELSE
         LTRLVL = -999
      END IF
      TIMITR = SECOND() - TIMITR
      IF (NWOPT.EQ.0 .OR. FLAG(34)) THEN
C        ... We do not need the integrals for optimization.
C        ... flag(34) = (int.transf. not needed in optimization)
         FLAG(14) = .TRUE.
      ELSE
         FLAG(14) = .FALSE.
      END IF
C
C *** Retrieve CI coefficients:
C
      CALL RDCREF(WRK(KCREF))
C
C         --- Normalize CI vector (i.e. CREF)
C         6-Dec-1984-hjaaj: SIRSAV changed so CREF normally normalized now
C
      IF (NCONF .GT. 1) THEN
         C0NORM = DNORM2(NCONF,WRK(KCREF),1)
         THREQL = SQRT(NCONF*THRLDP)
         IF (ABS(C0NORM-D1) .GT. THREQL) THEN
            IF (C0NORM.EQ.D0) THEN
               WRITE (LUW4,1110)
               IF (LUPRI.NE.LUW4) WRITE (LUPRI,1110)
               WRITE (LUERR,1110)
               CALL QUIT('ERROR (SIROPT) CI vector has zero norm')
            END IF
            C0NORM = D1/C0NORM
            NWARN = NWARN + 1
            WRITE (LUPRI,1100) C0NORM
            CALL DSCAL(NCONF,C0NORM,WRK(KCREF),1)
         END IF
      END IF
 1100 FORMAT(/' WARNING: Reference CI vector normalized, factor:',
     &       F15.10)
 1110 FORMAT(/' Optimization control: FATAL ERROR',
     *        ' - zero norm of reference CI vector')
C
      IF (P6FLAG(32) .AND. NCONF.GT.1) THEN
C        use PRWF here ?? 890925-hjaaj
         WRITE (LUPRI,2010) ITMAC,THRPWF
 2010    FORMAT(/' Reference CI vector in iteration',I3,
     *          /' -----------------------------------',
     *         //' Cutoff for print:',F7.4,/)
         CALL PRVEC(NCONF,WRK(KCREF),1,THRPWF,200,LUPRI)
         CALL PRMGN(NCONF,WRK(KCREF),1,12,LUPRI)
#if defined (VAR_OLDCODE) || defined (VAR_SECSEC)
         DO 202 I = 1,NCONF
            CREFI = WRK((KCREF-1)+I)
            IF (ABS(CREFI) .GE. THRPWF) WRITE (LUPRI,2020) I,CREFI
  202    CONTINUE
 2010    FORMAT(/' Reference CI vector in iteration',I3,
     *          /' -----------------------------------',
     *         //' Cutoff for print:',F7.4,
     *         //' Configuration no.           value',
     *          /' -----------------           -----')
 2020    FORMAT(I16,F20.10)
#endif
      END IF
C
C *** Calculate gradient for CREF and some matrices/arrays needed later.
C
C Dec99 hjaaj: set screening threshold for direct SCF gradient
C              GNORM(3) = grdnrm from prev.iter. (THRGRD in iter 1)
C              (based on code in sirdiis.F)
C Nov07 hjaaj: we expect quadratic convergence and GNORM(3) is from
C              prev.iter., thus GNORM(3)**2 
C              NORBT factor: statistical error is ca. sqrt(# elements) ~ NORBT
C
      IF (DIRFCK .AND. .NOT. USRSCR) THEN
         IF (UPDATE) THEN
           ! Hessian update gives max a factor 0.1, so 0.01 should be safe
           GRDACC  = MAX(THRGRD,0.01D0*GNORM(3))
         ELSE
           ! quadratic convergence; IFTHRS .ge. 9 takes care of linear region
           GRDACC  = MAX(THRGRD,GNORM(3)**2)
         END IF
         GRDACC  = GRDACC / NORBT
         IFTHRS  = MIN( IFTHSV, NINT(-LOG10(GRDACC)) + 1)
         IFTHRS  = MAX( 9, IFTHRS )
         IF (IPRI6 .GE. 0) WRITE (LUPRI,'(/A,I3,A)')
     &   ' Fock matrix screening setting for this macro iteration: 10^('
     &   ,-IFTHRS,')'
      END IF
      IF (.NOT.UPDATE) THEN
         CALL GRAD (WRK(KCREF),WRK(KCMO),WRK(KCINDX),WRK(KDV),
     *              WRK(KPV),WRK(KFC),WRK(KFV),WRK(KFQ),WRK(KFCAC),
     *              WRK(KH2AC),WRK(KGTOT),EMY,EACTIV,
     &              WRK(KWRKG),LWRKG)
C        CALL GRAD (CREF,CMO,INDXCI,DV,PV,FC,FV,FQ,
C    *              FCAC,H2AC,G,EMY,EACTIV,WRK,LFREE)
      ELSE
         CALL UPDGRD(0,WRK(KCREF),WRK(KCMO),WRK(KDV),WRK(KPV),
     *               WRK(KGTOT),EMY,EACTIV,WRK(KFC),WRK(KFV),WRK(KFQ),
     *               DUMMY,WRK(KCINDX),WRK(KWRKG),LWRKG)
C        CALL UPDGRD(IGRTYP,REF1,CMO1,DV,PV,GRD,EMY,
C    *               EACTIV,FC,FV,DIAOR,INDXCI,WRK,LFREE)
      END IF
C
C     hj-861130 solvent contribution, ESOLT is saved in /INFOPT/
C
      IF (SOLVNT) THEN
         CALL SOLGRD(WRK(KCREF),WRK(KCMO),WRK(KCINDX),WRK(KDV),
     &               WRK(KGTOT),ESOLT,WRK(KERLM),WRK(KWRKG),LWRKG)
C        CALL SOLGRD(CREF,CMO,INDXCI,DV,G,ESOLT,ERLM,WRK,LFREE)
      ELSEIF (PCM) THEN
         CALL PCMGRD(WRK(KCREF),WRK(KCMO),WRK(KCINDX),WRK(KDV),
     %               WRK(KGTOT),ESOLT,WRK(KWRKG),LWRKG)
      ELSE
         ESOLT = D0
      END IF
C----------------
C CBN+JK 03.01.06
C----------------
      IF (QM3) THEN
         CALL QM3GRAD(WRK(KCREF),WRK(KCMO),
     &                WRK(KCINDX),WRK(KGTOT),
     &                EQM3,WRK(KWRKG),LWRKG,IPRSOL)
      ELSE
         EQM3=D0
      END IF
C----------------
C CBN+JK 03.01.06
C----------------
      GNRMSV = GNORM(3)
      CALL GRDINF(GNORM,WRK(KGTOT))
C
C     HJ-860808 no absorption if orb.grd. .lt. fakabs*ci.grd.
C
      IF ( ITABS .EQ. 0 .AND. GNORM(2) .LT. FAKABS*GNORM(1) ) THEN
         IF (IPRI6.GE.5) WRITE (LUPRI,'(/A/A/)')
     *   ' *** no orbital absorption this macro',
     *   '     orbital gradient significantly smaller than CI gradient.'
         ITABS = -1
      END IF
C
C----------------
      EMCSCF = POTNUC + EMY + EACTIV + ESOLT + EQM3
C----------------
C
      IF (ITABS .GT. 0) THEN
         WRITE (LUW4,1005) EMCSCF,ITABS
         IF (LUPRI.NE.LUW4) WRITE (LUPRI,1005) EMCSCF,ITABS
      ELSE
         IF (DOHFSRDFT.OR.DOMCSRDFT) THEN
            WRITE (LUW4,1011) EMCSCF,ITMAC
            IF (LUPRI.NE.LUW4) WRITE (LUPRI,1011) EMCSCF,ITMAC
         ELSE
            WRITE (LUW4,1010) EMCSCF,ITMAC
            IF (LUPRI.NE.LUW4) WRITE (LUPRI,1010) EMCSCF,ITMAC
         END IF
      END IF
      IF (DOHFSRDFT.OR.DOMCSRDFT) THEN
         EMYTMP = EMY - ESRDFTY
      ELSE
         EMYTMP = EMY
      ENDIF
      IF (P4FLAG(2)) THEN
         WRITE (LUW4,1015) POTNUC,EMYTMP,EACTIV
         IF (DOHFSRDFT.OR.DOMCSRDFT) WRITE(LUW4,1018) ESRDFTY
      END IF
chj   IF (P6FLAG(1)) THEN
         IF (LUPRI.NE.LUW4 .OR. .NOT.P4FLAG(2)) THEN
            WRITE (LUPRI,1015) POTNUC,EMYTMP,EACTIV
            IF (DOHFSRDFT.OR.DOMCSRDFT) WRITE(LUPRI,1018) ESRDFTY
         END IF
chj   END IF
      IF (SOLVNT) WRITE (LUW4,1016) ESOLT
      IF (PCM)    WRITE (LUW4,1019) ESOLT
      IF (QM3)    WRITE (LUW4,1017) EQM3
      WRITE (LUW4,1020) GNORM(3),GNORM(1),GNORM(2)
      IF (LUPRI.NE.LUW4) THEN
         IF (SOLVNT.OR.PCM) WRITE (LUPRI,1016) ESOLT
         IF (QM3) WRITE (LUPRI,1017) EQM3
         WRITE (LUPRI,1020) GNORM(3),GNORM(1),GNORM(2)
      END IF
 1005 FORMAT(/' Total MCSCF energy       :',F30.15,T54,
     *        '(after absorp',I4,')' )
 1010 FORMAT(/' Total MCSCF energy       :',F30.15,T60,
     *        '(MACRO ',I4,')' )
 1011 FORMAT(/' Total MC-SRDFT energy    :',F30.15,T60,
     *        '(MACRO ',I4,')' )
 1015 FORMAT(/' - Nuclear repulsion      :',F30.15,
     *       /' - Inactive energy        :',F30.15,
     *       /' - Active energy          :',F30.15)
 1016 FORMAT( ' - Dielec. solvation ener.:',F30.15)
 1017 FORMAT( ' - QM/MM solvation energy :',F30.15)
 1018 FORMAT( ' - SRDFT effective energy :',F30.15)
 1019 FORMAT( ' - PCM solvation energy   :',F30.15)
 1020 FORMAT(/' Norm of total gradient   :',F27.12,
     *       /' -    of CI gradient      :',F27.12,
     *       /' -    of orbital gradient :',F27.12)

      IF (lim_poppri .gt. 0) THEN
         call sirpop('MCITER',WRK(KDV),WRK(KWRKG),LWRKG)
      END IF

C
C
      IF (SOLVNT .AND. IPRI6 .GT. 1) THEN
         WRITE (LUPRI,'(//A/)')
     *   ' Multipole analysis of dielectric solvation energy :'
         JERLM = KERLM
         DO 244 LVAL = 0,LSOLMX
            NMVAL = 2*LVAL + 1
            ERLTOT = DSUM(NMVAL,WRK(JERLM),1)
            WRITE (LUPRI,'(A5,I3,T12,F12.8)')
     *         '  L =',LVAL,ERLTOT
            IF (IPRI6 .GT. 5) THEN
               WRITE(LUPRI,'(A11,5F12.8/,(T12,5F12.8))')
     *            ' Components', (WRK(JERLM+MVAL),MVAL=0,(NMVAL-1))
               WRITE(LUPRI,'()')
            END IF
            JERLM = JERLM + NMVAL
  244    CONTINUE
      END IF
C
      IF (P4FLAG(15) .AND. NCONF .GT. 1) THEN
CHJ: in future version maybe CALL PRWF here (and above for CREF)
         PTEST = MAX(1.D-10,GNORM(1)*THRCGR)
         WRITE (LUW4,3010) ITMAC, PTEST
         CALL PRVEC(NCONF,WRK(KGTOT),1,-PTEST,200,LUW4)
         CALL PRMGN(NCONF,WRK(KGTOT),1,12,LUW4)
#if defined (VAR_OLDCODE) || defined (VAR_SECSEC)
         DO 301 I = 1,NCONF
            GCI = WRK((KGTOT-1)+I)
            IF (ABS(GCI) .GE. PTEST) WRITE (LUW4,3011) I,GCI
  301    CONTINUE
 3011 FORMAT(I16,F20.10)
#endif
      END IF
      IF (P4FLAG(14) .AND. NWOPT.GT.0) THEN
         WRITE (LUW4,3020) ITMAC
         IF (MPRI4.GT.100) THEN
            PRFAC = 0.0D0
         ELSE IF (MPRI4.GT.10) THEN
            PRFAC = 0.1D0
         ELSE
            PRFAC = 0.2D0
         END IF
         CALL PRKAP (NWOPT,WRK(KGTOT+NCONF),PRFAC,LUW4)
      END IF
 3010 FORMAT(//' Configuration gradient iteration no.',I3,
     *        /' ---------------------------------------'
     *       //' Cutoff for print:',1P,D10.2,
     *       //' Configuration no.           value'
     *        /' -----------------           -----')
 3020 FORMAT(//' Orbital gradient iteration no.',I3,
     *        /' ---------------------------------')
C
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
C
C *** Test for convergence and trust region size
C
      IF (GNORM(3) .LE. THRGRD) THEN
         ICONV = 1
         ISTEP = -2
      END IF
C
      IF (FLAG(20)) GO TO 700
C
      IF (ITMACN .GE. MAXMAC) ISTEP = -2
C
      IF ( CHKSTP ) THEN
C        istep = -2 : converged or stop after gradient,
C                     get step info for analysis
C        istep = -1 : restart prediction available
C        istep =  0 : normal second-order step check
         CALL SIRSTP(ISTEP,WRK(KCMO),WRK(KIBNDX),WRK(KCINDX),
     *               WRK(KWNEO),LWNEO)
C        CALL SIRSTP(ISTEP,CMO,IBNDX,INDXCI,WRK,LFREE)
C
         IF ( ICONV .NE. 1) THEN
C        Check if RATIO is OK for local UPDATE/DONR
C        NB! Only disable DONR if not flag(39) ".NR ALWAYS"
          IF ( UPDATE .OR. (DONR .AND. .NOT. NRALW) ) THEN
C         if ( update ) then this is local update,
C           otherwise CHKSTP would be false
            RATCHK = ABS(DINFO(6) - D1)
            IF (RATCHK .GT. 0.02) THEN
               IF (UPDATE) THEN
                   WRITE (LUPRI,'(/A,F10.5/A)')
     &     ' Local UPDATE disabled because 0.02 .lt. abs(ratio-1) =',
     &            RATCHK, ' We must therefore recalculate gradient.'
                  UPDATE = .FALSE.
                  DONR = NRALW
Chj-sep99: next 3 lines work, the rest from go to 100 to 'ELSE' does not, why???
                  CHKSTP = .FALSE.
                  itrlvl = itrlsv
                  go to 100
Chj-sep99: previous 3 lines work, the rest from go to 100 to 'ELSE' does not, why???
#ifdef HJ_CHECK_LATER
                  IF (ITRLSV .GT. ITRLVL) THEN
C                    need higher int. transf. level for NEO/NR
                     ITRLVL = ITRLSV
                     JTRLVL = ITRLVL
                     CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KWRK10),LWRK10)
                     LTRLVL = JTRLVL
                  END IF
                  CALL GRAD (WRK(KCREF),WRK(KCMO),WRK(KCINDX),WRK(KDV),
     *                  WRK(KPV),WRK(KFC),WRK(KFV),WRK(KFQ),WRK(KFCAC),
     *                  WRK(KH2AC),WRK(KGTOT),EMY,EACTIV,
     &                  WRK(KWRKG),LWRKG)
                  IF (SOLVNT) THEN
                     CALL SOLGRD(WRK(KCREF),WRK(KCMO),WRK(KCINDX),
     *                  WRK(KDV),WRK(KGTOT),ESOLT,
     *                  WRK(KERLM),WRK(KWRKG),LWRKG)
                  ELSEIF (PCM) THEN
                     CALL PCMGRD(WRK(KCREF),WRK(KCMO),WRK(KCINDX),
     *                  WRK(KDV),WRK(KGTOT),ESOLT,
     *                  WRK(KWRKG),LWRKG)
                  END IF
#endif
               ELSE
                  WRITE (LUPRI,'(/A,F10.5)')
     &            ' Local NR disabled because 0.02 .lt. abs(ratio-1) =',
     &            RATCHK
                  DONR = .FALSE.
               END IF
            END IF
          END IF
         END IF
      ELSE
         EMCOLD = EMCSCF
         ISTEP  = 0
      END IF
      DINFO(9) = RTRUST
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
C     --- Exit if converged ---
      IF (ICONV .EQ. 1) GO TO 700
C     --- Exit if max macro iter reached ---
      IF (ITMACN .GE. MAXMAC) GO TO 700
C
      IF (ISTEP.GT.0) THEN
         ISTEP  = 0
         UPDATE = .FALSE.
         DONR   = NRALW
         ITRLVL = ITRLSV
         ITBCK  = ITBCK + 1
         IF (ITBCK.LE.MAXBCK) THEN
            WRITE(LUW4,9000)
            IF (LUPRI.NE.LUW4) WRITE(LUPRI,9000)
         ELSE
            WRITE (LUERR,9100) ITBCK-1
            WRITE (LUW4 ,9100) ITBCK-1
            IF (LUPRI.NE.LUW4) WRITE (LUPRI,9100) ITBCK-1
         END IF
         GO TO 700
      ELSE IF (ISTEP.LT.0) THEN
C        close to convergence, no step control performed
         IF (ITMACN .GT. 1 .AND. GNORM(3) .GT. GNRMSV) THEN
C            ... hjaaj: no test if restart (ITMACN .eq. 1)
            WRITE (LUW4,'(/A/A)')
     &         ' ERROR EXIT: gradient norm increasing in local region.',
     &         ' Probable cause: numerical round-off errors'//
     &         ' or too many integrals neglected.'
            IF (LUPRI .NE. LUW4) WRITE (LUPRI,'(/A/A)')
     &         ' ERROR EXIT: gradient norm increasing in local region.',
     &         ' Probable cause: numerical round-off errors'//
     &         ' or too many integrals neglected.'
c           GO TO 750
            write(lupri,*) 'Well, let us try anyway ...'
         END IF
         ITBCK = 0
      ELSE
         ITBCK = 0
      END IF
 9000 FORMAT(//' Optimization control: Step was too large, we backstep')
 9100 FORMAT(//' Optimization control: ERROR EXIT NO CONVERGENCE'
     *       //' Step was too large and maximum number of backsteps,',
     *         I2,', reached'/)
C
C     Check ratio for geometry optimization (ABACUS interface)
C
      IF (GEOWLK .AND. JWSTEP .EQ. 1) THEN
         IF (GNORM(3) .LE. THRGEO) THEN
            DEPSAV = DEPRED
            EMCOLD = EMCGEO
            DEPRED = DEPGEO
            CALL SIRSTP(JWSTEP,DUMMY,DUMMY,DUMMY,DUMMY,1)
            EMCOLD = EMCSCF
            DEPRED = DEPSAV
            IF (WLKREJ) GO TO 700
         END IF
      END IF
C
C **********
C We have not converged (yet),
C
C
C if (ABSORP) then call SIRORB for orbital optimization.
C
      CHKSTP = .TRUE.
      IF (ITABS .EQ. 0) THEN
         CHKSTP = .FALSE.
         CALL SIRORB(MAXAPM,ITABS,WRK(KCMO),WRK(KGTOT+NCONF),WRK(KDV),
     &      WRK(KPV),WRK(KFC),WRK(KFV),WRK(KFQ), WRK(KWOCTL),LWOCTL)
C        CALL SIRORB(MAXITO,NITORB,CMO,GORB,DV,PV,FC,FV,FQ,WRK,LFREE)
         GO TO 700
C        ... change 870607/hjaaj, was GO TO 100
C            use GO TO 700 to get timing statistics.
       ELSE IF (ABSORP) THEN
         ITABS  =  0
      END IF
C
C If (UPDATE) then call SIRUPD to do UPDATE optimization,
C else call SIRNEO to find next step vector (using the NEO
C optimization technique).
C
      IF (UPDATE) THEN
         CHKSTP = .FALSE.
         ITRLVL = 0
         NWOPH  = NWOPT
         NVARH  = NVAR
         CALL DCOPY(NCONF,WRK(KCREF),1,WRK(KREFU),1)
         CALL DCOPY(NVAR,WRK(KGTOT),1,WRK(KGRDU),1)
         CALL SIRUPD(ICONV,MAXUIT,WRK(KREFU),WRK(KGRDU),WRK(KCMO),
     &               WRK(KCINDX),WRK(KWU),LWU)
C        CALL SIRUPD(ICONV,MAXIT,REF1,GRD,CMO1,INDXCI,WRK,LFREE)
         LTRLVL = ITRLVL
         ITRLVL = ITRLSV
         NWOPH  = NWOPS
         NVARH  = NVARS
         NREDL  = 0
         UPDATE = .FALSE.
         IF (ICONV .GT. 0) THEN
            IF (FLAG(25)) THEN
               WRITE (LUPRI,'(//A/)')
     &            ' >>> after UPDATE: recalculating gradient for SIRIFC'
               ICONV    = 0
               UPDATE   = .FALSE.
               IF (.NOT. RHFCLC) FLAG(14) = .TRUE.
               GO TO 700
            END IF
            GO TO 800
         ELSE
            IF (ISTATE .EQ. 1 .AND. NROOTS .EQ. 1) THEN
               GO TO 700
C              ... change 870607/hjaaj, was GO TO 100
C                  use GO TO 700 to get timing statistics.
            ELSE
               WRITE (LUW4,'(/A)')
     &  ' ERROR EXIT: NEO cannot proceed after UPDATE for NROOTS .gt. 1'
               IF (LUPRI .NE. LUW4) WRITE (LUPRI,'(/A)')
     &  ' ERROR EXIT: NEO cannot proceed after UPDATE for NROOTS .gt. 1'
               GO TO 750
C              ... too few start vectors on LUIT1
            END IF
         END IF
      END IF
C
C
C ************************************************************
C     second order step
C ************************************************************
C
C
      IF (DONR) THEN
         CALL SIRNR (WRK(KCREF),WRK(KGTOT),WRK(KCMO),WRK(KCINDX),
     *               WRK(KDV),WRK(KPV),
     *               WRK(KFC),WRK(KFV),WRK(KFCAC),
     *               WRK(KH2AC),EACTIV,WRK(KIBNDX),
     &               WRK(KWNEO),LWNEO)
C        CALL SIRNR (CREF,G,CMO,INDXCI,DV,PV,
C    *               FC,FV,FCAC,H2AC,EACTIV,IBNDX,WRK,LFREE)
C
C        Absorption should not be used after SIRNR
C        (because SIRNR can not check if root flip has occurred
C         for excited states, and anyway SIRNR is normally only
C         used in the local region).
C
         ABSORP   = .FALSE.
         FLAG(51) = .FALSE.
         FLAG(52) = .FALSE.
         FLAG(53) = .FALSE.
      ELSE
         CALL SIRNEO(WRK(KCREF),WRK(KGTOT),WRK(KCMO),WRK(KCINDX),
     *               WRK(KDV),WRK(KPV),WRK(KFC),WRK(KFV),WRK(KFCAC),
     *               WRK(KH2AC),EACTIV,WRK(KIBNDX),
     &               WRK(KWNEO),LWNEO)

C        CALL SIRNEO(CREF,G,CMO,INDXCI,DV,PV,
C    *               FC,FV,FCAC,H2AC,EACTIV,IBNDX,WRK,LFREE)
C
C        Test if SIRNEO has disabled absorption
C        (because of root flipping for excited state optimization)
C
         IF (ABSORP) ABSORP = ( FLAG(51) .OR. FLAG(52) .OR. FLAG(53) )
      END IF
C
C     Check if local region where absorption normally is
C     switched off (unless FLAG(54) has been set in input)
C
      IF (ABSORP .AND. .NOT.FLAG(54)) THEN
         LOCAL = (STPLEN .LT. STPABS*RTRUST) .AND. ITBCK .EQ. 0
         IF (LOCAL) ITABS = -1
C        ... no absorption next iteration
         LOCAL = (LOCAL .AND. STPC .LT. D2*STPO)
C        ... only switch absorption completely off if also
C            STPC .lt. 2*STPO (else we may get a large orbital step
C            again).
      IF (LOCAL .OR. ITMACN .GE. MAXABS) THEN
         WRITE (LUPRI,'(//A)') ' *** Orbital absorption switched off'
         IF (LOCAL) THEN
            WRITE (LUPRI,'(A/)')
     *      '     because local region has been reached.'
         ELSE
            WRITE (LUPRI,'(A/)')
     *      '     because max macro iterations w/ absorption reached.'
         END IF
         ABSORP   = .FALSE.
         FLAG(51) = .FALSE.
         FLAG(52) = .FALSE.
         FLAG(53) = .FALSE.
      END IF
      END IF
C
C     Test if switch to update method or switch to NR linear equations
C     (condition: local region, defined as step .le. STPLCL*RTRUST)
C
      IF (.NOT.FLAG(37) .OR. .NOT.FLAG(38) .OR. .NOT.NRALW) THEN
C     If (not noupdate or not neo_always or not nr_always) then
         LOCAL = (STPLEN .LT. STPLCL*RTRUST .AND. STPC .LT. STPO)
         I = ISTATE + 1 - IINFO(13)
         SHFLVL = DINFO(20+I)
         IF (IPRSTAT .GT. 5) WRITE (LUSTAT,'(/A,F20.10)')
     *      ' Level shift used for UPDATE/SIRNR test:',SHFLVL
         LOCAL = LOCAL .AND. (ABS(SHFLVL) .LT. SHFLCL)
         LOCAL = LOCAL .AND. ITBCK .EQ. 0
         IF ( LOCAL ) THEN
            IF (.NOT. FLAG(37)) THEN
               WRITE (LUPRI,'(/A)') ' Local region, UPDATE activated'
               UPDATE = .TRUE.
               DONR   = .FALSE.
               ITRLVL = 0
            ELSE IF (ISTATE .GT. 1) THEN
C              ... no reason to switch to NR for ground state
               WRITE (LUPRI,'(/A)') ' Local region, DONR activated'
               DONR     = .TRUE.
               FLAG(39) = .TRUE.
            END IF
         END IF
      END IF
C
C save information for final summary print-out
C
  700 CONTINUE
      TIMMAC = SECOND() - TIMMAC
      IF (IPRI6 .GT. 0) WRITE (LUPRI,'(/A,F15.2,A)')
     *   ' >>> Time used in this macro iteration :',TIMMAC,' seconds.'
C
      IINFO(1)  = ITMAC
      DINFO(1)  = EMY
      DINFO(2)  = EACTIV
      DINFO(3)  = EMCSCF
      DINFO(16) = TIMMAC
      DINFO(18) = TIMITR
      WRITE (LUINF) DINFO,IINFO
C
      IF (IPRSTAT .GT. 1) THEN
         WRITE (LUSTAT,'(//A/A/)')
     *   ' Iteration no., total time, integral transformation time,',
     *   ' micro it. time, and lin.trn. time'
         WRITE (LUSTAT,'(I5,4F10.2)')
     *   ITMAC,TIMMAC,TIMITR,DINFO(17),DINFO(19)
      END IF
C     ( DINFO(17) = TIMMIC and DINFO(19) = TIMLIN )
C
C CI case special output
C
      IF (NWOPT.EQ.0 .AND. NCONF.GT.1) THEN
         WRITE (LUW4,7011) DEPRED,(EMCSCF+DEPRED-POTNUC),EMCSCF+DEPRED
         IF (LUPRI.NE.LUW4)
     *   WRITE (LUPRI,7011) DEPRED,(EMCSCF+DEPRED-POTNUC),EMCSCF+DEPRED
      END IF
 7011 FORMAT(/' (SIROPT) CI energy lowering:  ',F30.15,
     *       /T11,'CI electronic energy:',F30.15,
     *       /T11,'CI total energy:     ',F30.15)
C
      IF (ICONV .GT. 0) GO TO 800
      IF (FLAG(20))     GO TO 890
      IF (ITBCK .GT. 0 .AND. ITBCK.LE.MAXBCK) GO TO 100
      IF (ITBCK .EQ. 0 .AND. ITMACN.LT.MAXMAC .AND.
     *    ITMICT.LT.MAXMIC) GO TO 100
C
C ***********************************
C *** END OF MACRO ITERATION LOOP ***
C ***********************************
C
C We did not converge...
C
      NWARN = NWARN + 1
      TOPT = SECOND() - TOPT
      WRITE (LUSTAT,7050) ITMACN,MAXMAC,ITMICT,MAXMIC,ITBCK,MAXBCK,TOPT
      WRITE (LUW4,7050)  ITMACN,MAXMAC,ITMICT,MAXMIC,ITBCK,MAXBCK,TOPT
      IF (LUPRI.NE.LUW4)
     *WRITE (LUPRI,7050)  ITMACN,MAXMAC,ITMICT,MAXMIC,ITBCK,MAXBCK,TOPT
 7050 FORMAT(
     *   /' *** Optimization control WARNING: MCSCF not converged ***',
     *   /5X,  'Maximum number of iterations or backsteps reached:',
     *   /5X,  'Number of macro iterations used',T45,I5,
     *   /5X,  '   Maximum',T45,I5,
     *   /5X,  'Number of micro iterations used',T45,I5,
     *   /5X,  '   Maximum',T45,I5,
     *   /5X,  'Number of back steps this macro',T45,I5,
     *   /5X,  '   Maximum',T45,I5,
     *   /5X,  'Total number of CPU seconds used',F13.2)
      GO TO 890
  750 CONTINUE
      NWARN = NWARN + 1
      TOPT = SECOND() - TOPT
      WRITE (LUSTAT,7051) ITMACN,ITMICT,TOPT
      WRITE (LUW4,7051)  ITMACN,ITMICT,TOPT
      IF (LUPRI.NE.LUW4)
     *WRITE (LUPRI,7051)  ITMACN,ITMICT,TOPT
 7051 FORMAT(
     *   /' *** Optimization control WARNING: MCSCF not converged ***',
     *   /5X,  'Number of macro iterations used',T45,I5,
     *   /5X,  'Number of micro iterations used',T45,I5,
     *   /5X,  'Total number of CPU seconds used',F13.2)
      GO TO 890
C
C ***********************************************************
C
C We have converged...
C
  800 CONTINUE
      TOPT = SECOND() - TOPT
      WRITE (LUW4,8000) ITMAC,ITMICT,TOPT
      IF (LUPRI.NE.LUW4) WRITE (LUPRI,8000) ITMAC,ITMICT,TOPT
 8000 FORMAT(/' *** Optimization control: MCSCF converged ***',
     *       /5X,  'Number of macro iterations used',T45,I5,
     *       /5X,  'Number of micro iterations used',T45,I5,
     *       /5X,  'Total number of CPU seconds used',F13.2)
C
      IF (SOLVNT) THEN
         KFCSOL = KWRKG
         KWRKG  = KFCSOL + NNORBT
         LWRKG  = LFREE  - (KWRKG-1)
         CALL SOLFCK(DUMMY,DUMMY,WRK(KFCSOL),WRK(KCMO),
     &               WRK(KERLM),.TRUE.,DUMSOL,
     &               WRK(KWRKG),LWRKG,IPRSOL)
         CALL DAXPY(NNORBT,D1,WRK(KFC),1,WRK(KFCSOL),1)
      ELSE IF (PCM) THEN
         KFCSOL = KWRKG
         KWRKG  = KFCSOL + NNORBT
         LWRKG  = LFREE  - (KWRKG-1)
         CALL PCMFCK(DUMMY,DUMMY,WRK(KFCSOL),WRK(KCMO),
     &               .TRUE.,DUMSOL, WRK(KWRKG),LWRKG,IPRSOL)
         CALL DSCAL(NNORBT,-1.0D0,WRK(KFCSOL),1)
         CALL DAXPY(NNORBT,D1,WRK(KFC),1,WRK(KFCSOL),1)
C------------------
C CBN+JK 03.01.06
C------------------
      ELSE IF (QM3) THEN
         KFCSOL = KWRKG
         KWRKG  = KFCSOL + NNORBT
         LWRKG  = LFREE  - (KWRKG-1)
         CALL QM3FCKMO(WRK(KCMO),WRK(KFCSOL),WRK(KWRKG),
     &                 LWRKG,IPRSOL)
         CALL DAXPY(NNORBT,D1,WRK(KFC),1,WRK(KFCSOL),1)
      ELSE
         KFCSOL = KFC
      END IF

      IF (NASHT .EQ. 0) THEN
C        hjaaj dec 2002: if extended to cover nasht.eq.1 as in sirdiis.F
C        then remember to add FCSOL and FV to FD for RHFENR call
         WRITE (LUW4,'(//A)') ' *** SCF orbital energy analysis ***'
         IF (SOLVNT .OR. QM3)
     &     WRITE (LUW4,'(A)') '    (incl. solvent contribution)'
         CALL RHFENR(IPRI4,LUW4,WRK(KFCSOL), WRK(KWRKG),LWRKG)
         IF (LUPRI .NE. LUW4) THEN
            WRITE (LUPRI,'(//A)') ' *** SCF orbital energy analysis ***'
            IF (SOLVNT .OR. QM3)
     &        WRITE (LUPRI,'(A)') '    (incl. solvent contribution)'
            CALL RHFENR(IPRI6,LUPRI,WRK(KFCSOL), WRK(KWRKG),LWRKG)
         END IF
C        CALL RHFENR(IPRINT,LUPRI,FD,SCRA,LSCRA)
      END IF
C
C
C *** LOCALIZATION
C
      IF (BOYORB.OR.PIPORB) CALL SIRLOC(WRK(KCMO),WRK(KWRKG),LWRKG)
C
      IF (LMULBS .OR. R12CAL) THEN
C        Construct orthonormal auxiliary basis sets (WK/UniKA/04-11-2002).
         LWRKG  = LFREE  - (KWRKG-1)
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL R12AUX(WRK(KWRKG),LWRKG)
         CALL TIMER('R12AUX',TIMSTR,TIMEND)
      END IF
C
C     If (FLAG(25) .OR. INERSI) write LUSIFC
C     modified hjaaj-aug99: only if this is final MO level
C       FLAG(25) is true for .ABACUS (geometry optimization)
C                         or .RESPONSE (response calculation)
C                         or .INTERFACE (request for LUSIFC)
C       INERSI is true if this is initial state calculation
C              in a solvent calculation with inertial polarization
C
      LSTLVL = 1
      IF (FLAG(21) .AND. DOMC) LSTLVL = 0
C     ... this is RHF which will be followed by MCSCF
      IF (LSTLVL .EQ. 1 .AND. (FLAG(25) .OR. INERSI)) THEN
         IF (GEOWLK) THEN
            IF (.NOT. OPTNEW) WRITE (LUW4,'(//A/)')
     &         ' Check ratio for geometry walk with converged MC :'
            JWSTEP = 1
            EMCOLD = EMCGEO
            DEPRED = DEPGEO
            CALL SIRSTP(JWSTEP,DUMMY,DUMMY,DUMMY,DUMMY,1)
         ELSE
            WLKREJ = .FALSE.
         END IF
C
Chj      if walk is rejected: do not update information on SIRIFC
         IF (.NOT.WLKREJ)
     &   CALL WR_SIRIFC(WRK(KCMO),WRK(KDV),WRK(KPV),
     &               WRK(KFCSOL),WRK(KFV),WRK(KFQ),
     &               WRK(KCREF),WRK(KFCAC),WRK(KH2AC),
     &               WRK(KWRKG),LWRKG,WRK(KGTOT+NCONF),WRK(KERLM),
     &              .TRUE.,WRK(KCINDX))
C        CALL WR_SIRIFC(CMO,DV,PV,FC,FV,FQ,CREF,FCAC,H2AC,WRK,LFREE,
C    &               GORB,ERLM,ORBHES,XINDX)
C        ... note: GORB will in a normal macro iteration be
C                  overwritten in SIRNEO or SIRNR, but none of
C                  these routines have been called after GRAD
C                  when the wave function is converged.
C
      END IF
C
C     Transform integrals for abacus or response or
C     for other programs, if requested.
C     modified hjaaj-aug99: only if this is final wf level
C
      LSTLVL = 1
      IF (FLAG(21) .AND. (DOMP2 .OR. DOCI .OR. DOMC)) LSTLVL = 0
C     ... if (true) this QCHF will be followed by higher wf level
      IF (LSTLVL .EQ. 1 .AND. ABS(ITRFIN) .LE. 10) THEN
         IF (IPRI6 .GE. 5) WRITE (LUPRI,'(/2A,I3,A)')
     *   ' >>> Transforming integrals for Abacus/response with',
     *   ' ITRLVL =',ITRFIN,' <<<'
         IF (ABS(ITRFIN) .EQ. LTRLVL) THEN
            IF (IPRI6 .GE. 5) WRITE (LUPRI,'(/A)')
     *         ' Transformation skipped, was already performed.'
         ELSE
            JTRLVL = ITRFIN
            CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KWRK10),LWRK10)
         END IF
      END IF
C
C ***********************************************************
C
  890 CONTINUE
C     restore control variables
      FLAG(39) = NRALW
      IFTHRS = IFTHSV
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
C
C write timing information for optimization to LUSTAT
C
      IF (IPRSTAT.GT.0) CALL TIMOPT('END',LUSTAT,WRK,LFREE)
C *** end of subroutine SIROPT
      CALL QEXIT('SIROPT')
      RETURN
      END
C  /* Deck rhfenr */
      SUBROUTINE RHFENR(IPRINT,LUPRI,FD,SCRA,LSCRA)
C
C 27-Oct-1990 Hans Joergen Aa. Jensen (revised 17-Dec-2002 hjaaj)
C
C Purpose:
C  Check if orbitals are canonical Hartree-Fock orbitals
C  and print orbital energies
C
C Input:
C  FD; the total inactive + active Fock matrix
C
C Scratch:
C  SCRA; general scratch area
C
#include <implicit.h>
      DIMENSION FD(*),SCRA(*)
C
C
      PARAMETER (D0 = 0.0D0,    EMYCNV = 1.D-4,
     &           DBIG = 1.D+12, GAPMIN = 0.1D0)
#include <dummy.h>
C
C  Used from common blocks:
C     INFORB : NISH(*), NISHT, ...
Chj-sep99: do not use NRHF(*) and NRHFT, because they are not
C          correct in DIIS if AUTOCC and occupations have changed.
C     SCBRHF : IOPRHF,NFRRHF(*)
C     INFIND : IROW(*),ISW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <molde.h>
#include <infinp.h>
#include <inforb.h>
#include <scbrhf.h>
#include <infind.h>
C
      CALL QENTER('RHFENR')
C
      IF (NASHT .GT. 0) THEN
C        hjaaj Dec 2002:
C        allow printing of orbital energies for open-shell systems.
         WRITE (LUPRI,'(/A/A/A/A)')
     &' Orbital energy analysis for an open-shell system.',
     &' Orbital energies are not well-defined for open-shell systems',
     &'   here is used block diagonalization of'//
     &  ' the FD=FC+FV Fock matrix.',
     &" NOTE that Koopmans' theorem is not fulfilled for this case."
C        WRITE (LUPRI,'(/A)')
C    *      ' Orbital energy analysis skipped for open-shell systems'
C        GO TO 9999
C hjaaj TODO:   EKT orbital energies may be defined later.
      END IF
      IF (IPRINT .LT. 10) WRITE (LUPRI,'(/A)') ' Only the five'//
     &   ' lowest virtual orbital energies printed in each symmetry.'
      IF (IPRINT .GE. 2) WRITE (LUPRI,'(/A,I5/A,8I5)')
     *   ' Number of electrons :',2*NISHT,
     *   ' Orbital occupations :',(NISH(I),I=1,NSYM)
C
      KORBEN = 1
      LNEED  = KORBEN + NORBT 
      IF (LNEED .GT. LSCRA) CALL ERRWRK('RHFENR',LNEED,LSCRA)
      CALL DZERO(SCRA,NORBT)
C
C     Print orbital energies.
C     Check if orbitals are canonical orbitals
C
      IF (SUPSYM) THEN
         IF (DODFT) THEN
            WRITE (LUPRI,1841)
         ELSE
            WRITE (LUPRI,1741)
         END IF
      ELSE
         IF (DODFT) THEN
            WRITE (LUPRI,1840)
         ELSE
            WRITE (LUPRI,1740)
         END IF
      END IF
      VALCON =  D0
      EHOMO  = -DBIG
      IHOMO  = 0
      ESOMO  = -DBIG
      ISOMO  = 0
      ELUMO  =  DBIG
      ILUMO  = 0
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         NFRHFI= NFRRHF(ISYM)
         IORBI = IORB(ISYM)
         ISSYM = IIORB(ISYM)
         EHOMOI= EHOMO
         ESOMOI= ESOMO
         ELUMOI= ELUMO
         DO 300 I = 1,NORBI
            IJOFF = ISSYM + IROW(I)
            SCRA(IORBI+I) = FD(IJOFF+I)
            IF (I .GT. NFRHFI) THEN
C           ... not frozen
               IW = ISW(IORBI+I)
               IF (IW.LE.NISHT) THEN 
C              ... closed shell (inactive) orbital
                  EHOMOI = MAX(EHOMOI,SCRA(IORBI+I))
                  DO J = NFRHFI+1,I-1
                     VALCON = MAX(VALCON,ABS(FD(IJOFF+J)))
                  END DO
               ELSE IF (IW.LE.NOCCT) THEN
C              ... open shell (active) orbital
                  ESOMOI = MAX(ESOMOI,SCRA(IORBI+I))
C                 No 'canonical' test for the open shell active orbitals
C                 because the Fock matrix is not zero for these
               ELSE
C              ... virtual (secondary) orbital
                  ELUMOI = MIN(ELUMOI,SCRA(IORBI+I))
                  DO J = NFRHFI+1,I-1
                     VALCON = MAX(VALCON,ABS(FD(IJOFF+J)))
                  END DO
C                 ... skipping the open shell active orbitals
C                     because the Fock matrix is not zero for these
                  DO J = NOCCT+1,I-1
                     VALCON = MAX(VALCON,ABS(FD(IJOFF+J)))
                  END DO
               END IF
            END IF
  300    CONTINUE
         IF (EHOMOI .GT. EHOMO) THEN
            EHOMO = EHOMOI
            IHOMO = ISYM
         END IF
         IF (ESOMOI .GT. ESOMO) THEN
            ESOMO = ESOMOI
            ISOMO = ISYM
         END IF
         IF (ELUMOI .LT. ELUMO) THEN
            ELUMO = ELUMOI
            ILUMO = ISYM
         END IF
         IF (IPRINT .GE. 10) THEN
            NLAST = NORBI
         ELSE
            NLAST = MIN(NORBI,NOCC(ISYM)+5)
         END IF
         IF (NFRHFI .GT. 0) WRITE (LUPRI,'(A,I3,A)')
     &' Note: the first',NFRHFI,' orbitals in next symmetry are frozen'
         IF (SUPSYM) THEN
            WRITE (LUPRI,1751)
     &         ISYM,(SCRA(IORBI+I),ISSMO(IORBI+I),I=1,NLAST)
         ELSE
            WRITE (LUPRI,1750) ISYM,(SCRA(IORBI+I),I=1,NLAST)
         END IF
  200 CONTINUE
      IF (MOLDEN) CALL MOMOS(2,SCRA,DUMMY,DUMMY,DUMMY,DUMMY)
C     take care of 1-electron case and no virtuals case
      IF (EHOMO .EQ. -DBIG) EHOMO = D0
      IF (ELUMO .EQ.  DBIG) ELUMO = D0
      WRITE (LUPRI,'(2(/A,F15.8,A,I2,A)/A/A,F15.8,A)')
     &   '    E(LUMO) :',ELUMO,' au (symmetry',ILUMO,')',
     &   '  - E(HOMO) :',EHOMO,' au (symmetry',IHOMO,')',
     &   '  ------------------------------------------',
     &   '    gap     :',ELUMO-EHOMO,' au'
      IF (NASHT .GT. 0) THEN
         WRITE (LUPRI,'(/A,F15.8,A,I2,A)')
     &   'and E(SOMO) :',ESOMO,' au (symmetry',ISOMO,')'
      END IF
      IF (ELUMO-EHOMO .LT. GAPMIN) THEN
         WRITE(LUPRI,'(/A)')
     &   ' INFO: E(LUMO) - E(HOMO) small or negative.'
      END IF
      IF ( VALCON .GT. EMYCNV ) THEN
         WRITE(LUPRI,'(/A//A,1P,D10.2)') ' NOTE:'//
     *   ' MOLECULAR ORBITALS ARE NOT CANONICAL HARTREE-FOCK ORBITALS',
     *   ' Largest off-diagonal Fock matrix element is',VALCON
         IF (IPRINT .GT. 10) THEN
            WRITE(LUPRI,'(/A)') ' THE TOTAL FOCK MATRIX :'
            CALL OUTPKB(FD,NORB,NSYM,1,LUPRI)
         END IF
      ELSE IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI,'(/A/A,1P,D10.2)')
     *      ' Deviation from canonical HARTREE-FOCK orbitals :',
     *      ' Largest off-diagonal Fock matrix element is',VALCON
      END IF
 9999 CALL QEXIT('RHFENR')
      RETURN
C
C1730 FORMAT(//' Hartree-Fock electronic energy:',F30.15)
C1732 FORMAT( /' Hartree-Fock total      energy:',F30.15)
 1740 FORMAT( /' Sym       Hartree-Fock orbital energies')
 1741 FORMAT( /' Sym       Hartree-Fock orbital energies',
     &         ' (supersymmetry in parenthesis)')
 1840 FORMAT( /' Sym       Kohn-Sham orbital energies')
 1841 FORMAT( /' Sym       Kohn-Sham orbital energies',
     &         ' (supersymmetry in parenthesis)')
 1750 FORMAT(/I3,1X,5F15.8/,(4X,5F15.8) )
 1751 FORMAT(/I3,1X,4(F15.8,' (',I2,')')/,(4X,4(F15.8,' (',I2,')')) ) 
C
C
C *** end of subroutine RHFENR
C
      END
C  /* Deck optst */
      SUBROUTINE OPTST (INDXCI,WRK,LFREE)
C
C  Written december 83 by Hans Joergen Aa. Jensen and Hans Agren
C  Revisions:
C    29-Jul-1984 hjaaj // 19-May-1984 hjaaj/ha
C     7-Jan-1985 hjaaj (ICI0=5 option)
C    22-Apr-1985 hjaaj (ICI0=6 option)
C     6-Jul-1985 hjaaj (CICTL for ICI0=1 option)
C     4-Aug-1986 hjaaj (removed obsolete ICI0=2, =3 options)
C
C MOTECC-90: Some of the algorithms used in OPTST are
C            described in Chapter 8 Section E.6 of MOTECC-90
C            "Initial CI iterations for Frozen Orbitals"
C
C Purpose:
C   Obtain CI vector for the expansion point in the first macro
C   iteration and a start guess for the other vectors in a
C   simultaneous expansion algorithm for solving the eigenvalue/
C   eigenvector problem. To write molecular orbitals on LUIT1.
C   The routine is controlled by the ICIO parameter.
C
C if (ICI0<0) then
C    select -ICIO as start configuration
C else
C
C    ICI0>100 : transform to natural orbitals
C               Should not be used, use FLAG instead.
C
C    mod(ICI0,100) --
C     =1 : Compute start guess of CI-vector from a new H-diagonal.
C          The diagonal is computed in and returned from the HDIAG
C          routine. This call is preceded by a transformation call.
C          The CI-vector is afterwards written on LUIT1.
C          If MAXCIT .gt. 0, or the ISTATE configuration is degenerate,
C          then call CICTL to do MAXCIT iterations (3 if deg.).
C
C     =2 : *OBSOLETE* (860804 hjaaj); now as ICI0 = 4 (880616 hjaaj)
C
C     =3 : *OBSOLETE* (860804 hjaaj); now as ICI0 = 4 (880616 hjaaj)
C
C     =4 : Start CI-vector is already on LUIT1 (which is checked).
C
C     =5 : Restart, as ICI0=4 plus read EMCOLD,BETA,RTRUST from LUIT1.
C          ('NEOSAVE' in SIRSAV)
C
C     =6 : Geometry walk, new geometry ('GEOSAVE' in SIRSAV)
C
C     =7 : Start using MC update information ('UPDSAVE' in SIRSAV)
C end if
C
#include <implicit.h>
C
      DIMENSION INDXCI(*),WRK(LFREE)
C
      PARAMETER ( D0=0.0D0, DP5=0.5D0, D1=1.0D0 )
#include <dummy.h>
C
C Used from common blocks:
C   INFINP: ISTATE,LROOTS,FLAG(*),MAXCIT,THRGRD,ITRLVL,?
C   INFVAR: NCONF (p.t. only NCONF)
C   INFORB: NCMOT,?
C   INFOPT: EMCOLD,BETA,RTRUST (read from LUIT1 when restart, ICI0=5)
C   INFDIM: IVEX4,INTOT
C   INFTAP: LUIT1,?
C   INFPRI: ?
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infopt.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
#include <gnrinf.h>
#include <dftcom.h>
C
C *** Externals:
C
      LOGICAL FNDLAB, FNDLB2
C
C *** Local variables:
C
      LOGICAL     HFCALC,    RSTFLG(4), TRCNO,     LSAV15
      LOGICAL     DFT_SPINDNS_SAVE, DFT_LOCALSPIN_SAVE
      CHARACTER*8 TABLE(8),  LBLINF(2), LBLDAT(2), STAR8
      CHARACTER*6 CNOLBL
C
C *** data:
C
      DATA TABLE/'EODATA  ', 'OLDORB  ', 'STARTVEC', 'CIDIAG2 ',
     *           'RESTART ', 'NEWORB  ', 'GEOWALK ', 'NATOCCAC'/
      DATA STAR8/'********'/
C
      CALL QENTER('OPTST ')
C     Transfer ICI0 to local variable ICSTRT
      ICSTRT = ICI0
C
      CALL GETDAT(LBLDAT(1),LBLDAT(2))
      NC4    = MAX(4,NCONF)
      NCMOT4 = MAX(4,NCMOT)
      TRCNO  = .FALSE.
      HFCALC = (NCONF .EQ. 1)
   10 CONTINUE
C
C *** Go to appropriate section to obtain start guess for CI vectors
C     (as specified with the input parameter ICI0 (here ICSTRT))
C     However, if RHF or NCONF.eq.1, simply write CREF(1) = D1 to LUIT1
C
      IF (HFCALC) THEN
         ICSTRT = MOD(ICSTRT,100)
         IF (ICSTRT .EQ. 1) ICSTRT = -1
         TRCNO  = FLAG(15)
         CNOLBL = 'ONLYFD'
      ELSE IF (FLAG(15)) THEN
         TRCNO  = .TRUE.
         IF (FLAG(48)) THEN
            CNOLBL = 'ONLYFD'
         ELSE IF (FLAG(46)) THEN
            CNOLBL = 'ONLYNO'
         ELSE
            CNOLBL = 'FD+NO '
         END IF
         ICSTRT = MOD(ICSTRT,100)
      ELSE IF (ICSTRT .GE. 100) THEN
         TRCNO  = .TRUE.
         CNOLBL = 'ONLYNO'
         ICSTRT   = MOD(ICSTRT,100)
      ELSE
         TRCNO  = .FALSE.
      END IF
C
C *** ICSTRT < 0 ***
C === Start guess = CSF no. -ICSTRT
C
      IF (ICSTRT.LT.0) THEN
         JCSTRT = -ICSTRT
C        ... 930521-hjaaj: JCSTRT introduced to
C            fix AIX xlf v.2.3 error for 'xlf -O'.
         IF (JCSTRT.GT.NCONF) THEN
            WRITE (LUPRI,5010) JCSTRT,NCONF
            WRITE (LUERR,5010) JCSTRT,NCONF
            CALL QTRACE(LUERR)
            CALL QUIT(
     &         'ERROR (OPTST) non-existent configuration specified')
         END IF
         IF (LROOTS.GT.1) THEN
            WRITE (LUPRI,5020) ICSTRT
            WRITE (LUERR,5020) ICSTRT
            CALL QTRACE(LUERR)
            CALL QUIT('(OPTST) ICSTRT must be pos. when LROOTS .gt. 1')
         END IF
         KCMO   = 1
         KCREF  = KCMO   + NCMOT
         KW11   = KCREF  + NCONF
         IF (KW11 .GT. LFREE) CALL ERRWRK('OPTST.CISTRT 0',KW11,LFREE)
C
C        Read NEWORB and write as OLDORB
C
         REWIND LUIT1
         CALL MOLLB2(TABLE(6),LBLINF,LUIT1,LUERR)
         CALL READT(LUIT1,NCMOT,WRK(KCMO))
         CALL NEWIT1
         WRITE (LUIT1) STAR8,LBLINF,TABLE(2)
         CALL WRITT(LUIT1,NCMOT4,WRK(KCMO))
C
C        Write STARTVEC
C
         CALL DZERO(WRK(KCREF),NCONF)
         WRK(KCREF-1+JCSTRT) = D1
         LBLINF(1) = '   1   1'
         LBLINF(2) = 'CISTRT 0'
         WRITE (LUIT1) STAR8,LBLINF,TABLE(3)
         CALL WRITT(LUIT1,NC4,WRK(KCREF))
C
C        Write NEWORB
C
         WRITE (LUIT1) STAR8,LBLDAT(1),LBLINF(2),TABLE(6)
         CALL WRITT(LUIT1,NCMOT4,WRK(KCMO))
C        single configuration has always natural orbitals
         IF (CNOLBL .EQ. 'ONLYNO') THEN
            TRCNO  = .FALSE.
         ELSE
            CNOLBL = 'ONLYFD'
         END IF
         GO TO 9000
      END IF
 5010 FORMAT(/' OPTST FATAL ERROR,'
     *       /' specified start configuration no.',I10,
     *        ' is non-existent (last is no.',I10,')')
 5020 FORMAT(/' OPTST FATAL ERROR,'
     *       /' negative CISTRT option not allowed',I10,
     *       /' when LROOTS is greater then one.')
C
C
C *** ICSTRT > 0 ***
C
C
      IF (ICSTRT .EQ. 2 .OR. ICSTRT .EQ. 3) ICSTRT = 4
C     ... 880616 hjaaj (CISTRT = 2 has now same meaning as in CIST).
      GO TO (100,50,50,400,400,400,400) ICSTRT
   50 CONTINUE
         WRITE (LUPRI,5030) ICSTRT
         WRITE (LUERR ,5030) ICSTRT
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR, illegal control option in OPTST')
 5030    FORMAT(/' OPTST FATAL ERROR, ICSTRT =',I4,' is not defined')
C
C *** ICSTRT = 1 ***
C *** Compute H-diagonal and then select starting configurations
C     (If MAXCIT .eq. 0 then CICTL selects smallest diag. elements of H;
C      if MAXCIT .gt. 0 then CICTL follows up with MAXCIT CI iterations)
C
  100 CONTINUE
C
         KCMO   = 1
         KW21   = KCMO   + NCMOT
         LW21   = LFREE  - KW21
         KECI   = KW21
         KICROO = KECI   + LROOTS
         KW22   = KICROO + LROOTS
         LW22   = LFREE  - KW22
         IF (LW22 .LT. 0) CALL ERRWRK('OPTST.CICTL',-KW22,LFREE)
C
C      A  single configuration will already have natural orbitals
         IF (MAXCIT .EQ. 0) THEN
            IF (CNOLBL .EQ. 'ONLYNO') THEN
               TRCNO  = .FALSE.
            ELSE
               CNOLBL = 'ONLYFD'
            END IF
         END IF
C
C        Read NEWORB and write as OLDORB
C
         REWIND LUIT1
         CALL MOLLB2(TABLE(6),LBLINF,LUIT1,LUERR)
         CALL READT(LUIT1,NCMOT,WRK(KCMO))
         CALL NEWIT1
         WRITE (LUIT1) STAR8,LBLINF,TABLE(2)
         CALL WRITT(LUIT1,NCMOT4,WRK(KCMO))
C
C ---    Call TRACTL first, if needed
C        (do it out here unless we later will transform to NO,
C         CICTL only does an ITRLVL = 0 transformation)
C
         IF (.NOT.FLAG(14) .AND. .NOT.TRCNO) THEN
            CALL TRACTL(ITRLVL,WRK(KCMO),WRK(KW21),LW21)
C           CALL TRACTL(ITRLVL,CMO,WRK,LFREE)
            FLAG(14) = .TRUE.
         END IF
C
         IF (ICHECK .GT. 0) THEN
            IF (MAXCIT .LE. 0) THEN
               WRITE (LUPRI,'(//A/A)')
     &         ' >>> INFO OPTST, symmetry check of CI vectors'//
     &         ' only works after CI iterations.',
     &         ' ".MAX CI ITERATIONS" reset to 3'
               MAXCIT = 3
            END IF
            IF (LROOTS .EQ. NROOTS) THEN
C           ... if (LROOTS .gt. NROOTS) then LROOTS has been
C               specified in input; use this value
               LROOTS = MAX(LROOTS,2*ISTATE+2)
               LROOTS = MIN(LROOTS,NCONF)
            END IF
         END IF
C
         THRCIX = DP5*THRGRD
         LSAV15   = FLAG(15)
         FLAG(15) = .FALSE.
         ISTASV   = ISTACI
         IF (ICHECK .EQ. 1 .AND. ISTATE .GT. 1) THEN
C        ... we do not know which root will become ISTATE
            ISTACI = 0
         ELSE
            ISTACI = ISTATE
         END IF
C
         IF (DOMCSRDFT) THEN
            ADDSRI    = .FALSE.
            DOCISRDFT = .TRUE.
            IF (DFT_SPINDNS .OR. DFT_LOCALSPIN)THEN
               WRITE(LUPRI,'(/A)')
     &           'INFO : spin density ignored in initial CI iterations.'
            END IF
            DFT_SPINDNS_SAVE = DFT_SPINDNS
            DFT_SPINDNS = .FALSE.
            DFT_LOCALSPIN_SAVE = DFT_LOCALSPIN
            DFT_LOCALSPIN = .FALSE.
         END IF
         CALL CICTL(1,LROOTS,MAXCIT,THRCIX,WRK(KCMO),INDXCI,
     *              WRK(KECI),ICONV,WRK(KICROO),WRK(KW22),LW22)
C        CALL CICTL(ICICTL,NCROOT,MAXITC,THRCIX,CMO,INDXCI,
C    *              ECI,ICONV,ICROOT,WRK,LFREE)
         IF (DOMCSRDFT) THEN
            ADDSRI    = .TRUE.
            DOCISRDFT = .FALSE.
            DFT_SPINDNS   = DFT_SPINDNS_SAVE
            DFT_LOCALSPIN = DFT_LOCALSPIN_SAVE
         END IF
         FLAG(15) = LSAV15
         ISTACI   = ISTASV
C
C ---    Call CICHCK to remove CI vectors that do not have the same
C        symmetry as ISTATE (ICHECK =1) or as the state of the lowest
C        symmetry (ICHECK=2).  (This call can only be made after some
C        CI-iterations have been done.)
C
         IF (ICHECK .GE. 0 .AND. MAXCIT .GT. 0) THEN
            CALL CICHCK(WRK,LFREE,ICHECK)
C           CALL CICHCK(WRK,LWRK,CCHECK)
         END IF
C
         IF (MAXCIT .GT. 0) LROOTS = NROOTS
C
C
      GO TO 9000
C
C *** ICSTRT = 4, 5, 6, or 7 ***
C     (See what codes means at beginning of this subroutine)
C
C === start guess should already be on LUIT1 (SIRIUS.RST)
C     (but we test to be sure!)
C
  400 CONTINUE
         REWIND LUIT1
         IF ( .NOT.FNDLAB(TABLE(2),LUIT1) ) GO TO 490
            READ (LUIT1,END=490,ERR=490)
         IF ( .NOT.FNDLB2(TABLE(3),LBLINF,LUIT1) ) GO TO 490
         IF ( ICSTRT.EQ.5 .AND. LBLINF(2).EQ.'NRSAVE  ')
     &      FLAG(39) = .TRUE.
C        ... if user has specified .RESTART, then respect shift to NR
         IF (FLAG(39)) GO TO 495
C        ... if (ONLYNR) one vector is all that is needed,
C            so no reason to do more checking.
C        We now check if sufficient vectors using label information.
         READ (LBLINF(1),'(2I4)',ERR=405) LROIT1,JSTA
         IF (LROIT1 .LT. LROOTS) GO TO 493
         GO TO 495
C        error reading LBLINF(1) -- use the old check
  405    LROIT1 = 0
         DO 410 I = 1,LROOTS
            READ (LUIT1,END=493,ERR=493) LBLINF
            IF (LBLINF(1) .EQ. STAR8) GO TO 493
            LROIT1 = LROIT1 + 1
  410    CONTINUE
      GO TO 495
C
  490 CONTINUE
         REWIND LUIT1
         IF (ICSTRT .EQ. 4) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//A/A//)')
     *       '@@*** (OPTST) MC start vectors not found on SIRIUS.RST',
     *       '@@*** WARNING lowest diagonal elements used.'
            ICSTRT = 1
            GO TO 10
         ELSE IF (ICSTRT .EQ. 5) THEN
            WRITE (LUPRI,'(//A/A//)')
     *       '@@*** (OPTST) MC start vectors not found on SIRIUS.RST,',
     *       '@@            .RESTART not possible.'
            CALL QTRACE(LUPRI)
            CALL QUIT(
     &         'MC start vectors for .RESTART not found on SIRIUS.RST.')
         ELSE
            WRITE (LUPRI,'(//A/A/A//)')
     *       ' *** (OPTST) MC start vectors not found on SIRIUS.RST,',
     *       '             no recover when CISTRT 6 or 7',
     *       '             optimization abandoned.'
            CALL QTRACE(LUPRI)
            CALL QUIT(
     &     '(OPTST) Start vectors not found on SIRIUS.RST as expected.')
         END IF
  493 CONTINUE
         REWIND LUIT1
         IF ( LROIT1 .GE. NROOTS) THEN
            WRITE (LUPRI,'(//A,I3,A/T14,A,I3,A,I3/)')
     *         ' *** (OPTST) only',LROIT1,' start vectors found ',
     *         'on SIRIUS.RST,', 'LROOTS reset from',LROOTS,' to',LROIT1
            LROOTS = LROIT1
         ELSE
            GO TO 490
         END IF
  495 CONTINUE
C
C ===================
C
C     READ RESTART INFORMATION
C
C ===================
C
      REWIND LUIT1
      IF (ICSTRT.EQ.5) THEN
C        "NEOSAVE" or "NRSAVE" restart;
C        1) get trust radius and beta
C        2) if (absorption) check if absorption switched off
         IF ( .NOT.FNDLB2(TABLE(5),LBLINF,LUIT1) ) THEN
          WRITE (LUPRI,'(//A)')
     &    ' (OPTST) Error, restart information not found on SIRIUS.RST.'
          CALL QUIT(
     &    '(OPTST) Error, restart information not found on SIRIUS.RST.')
         END IF
         READ (LUIT1) EMCOLD,DEPRED,BETA,RDUMMY,RTRUST,ITMAC
         WRITE (LUPRI,5000) EMCOLD,ITMAC,BETA,RTRUST,EMCOLD+DEPRED
      IF (LBLINF(2) .NE. STAR8) THEN
C        this is a new restart file (6-Aug-1986)
         READ (LUIT1) DUM,DUM,DUM,DUM,RSTFLG
         IF ((FLAG(51) .OR. FLAG(52) .OR. FLAG(53)) .AND.
     *      .NOT.FLAG(54)) THEN
            FLAG(51) = RSTFLG(1)
            FLAG(52) = RSTFLG(2)
            FLAG(53) = RSTFLG(3)
            IF (.NOT. (FLAG(51) .OR. FLAG(52) .OR. FLAG(53)) )
     *         WRITE (LUPRI,'(/A/)')
     *          ' *** absorption switched off in a previous iteration.'
         END IF
         IF (.NOT.FLAG(38) .AND. .NOT. FLAG(39)) THEN
C           if (not always neo and not always nr) then
            FLAG(39) = RSTFLG(4)
            IF (FLAG(39)) WRITE (LUPRI,'(/A/)')
     *          ' *** Switched from NEO to NR in a previous iteration.'
         END IF
      END IF
         TRCNO = .FALSE.
C        ... no TRCNO when restart, presumably done in prev. iteration
      ELSE IF (ICSTRT.EQ.6) THEN
C        "GEOSAVE" restart
         IF (ISTATE .GT. 1) THEN
            IF (FLAG(39)) THEN
               WRITE (LUPRI,6008) ISTATE
            ELSE
               WRITE (LUPRI,6010)
               CALL QTRACE(LUPRI)
               CALL QUIT('(OPTST) NEO and ISTATE .gt. 1 for GEO WALK')
            END IF
         END IF
         IF ( FNDLAB(TABLE(7),LUIT1) ) THEN
            READ  (LUIT1)      EMCOLD,DEPRED,REJWMI,REJWMA
            WRITE (LUPRI,6000)EMCOLD,DEPRED,EMCOLD+DEPRED,REJWMI,REJWMA
         ELSE 
C
C     We do not complain about lack of GEOWALK if .OPTIMI, K.Ruud, Nov.29-96
C     OK! But we must anyway reset ICSTRT and ICI0 to 4 hjaaj/961213.
C
            IF (.NOT. OPTNEW) THEN
               WRITE (LUPRI,'(//A/A/)')
     *         '@@ INFO : Label "GEOWALK " not found on SIRIUS.RST,',
     *         '@@ INFO : calculation continues without'//
     *         ' test for rejection of geometry change.'
            END IF
            ICSTRT = 4
            ICI0   = 4
C           ICI0 modified such that SIROPT does not perform
C           walk rejection check
         END IF
      ELSE IF (ICSTRT.EQ.7) THEN
C        "UPDSAVE" restart
         IF (ISTATE .GT. 1) THEN
            WRITE (LUPRI,7010)
            CALL QTRACE(LUPRI)
            CALL QUIT('ERROR (OPTST) ISTATE .gt. 1 for ICSTRT .eq. 7')
         END IF
      END IF
      REWIND LUIT1
      GO TO 9000
 5000 FORMAT(/' (OPTST) RESTART, Old MCSCF energy:',T50,F30.15,
     *       /T18,'for macro iteration no.',T49,I10
     *       /T23,'BETA value:',T49,F10.2,
     *       /T23,'Trust radius:',T44,F15.5,
     *       /T18,'Predicted new MCSCF energy:',T50,F30.15)
 6000 FORMAT(//' (OPTST) Geometry walk',
     *       /T16,'Old MCSCF energy:       ',F30.15,
     *       /T16,'Predicted energy change:',F30.15,
     *       /T16,'Predicted new energy:   ',F30.15,
     *      //T16,'Ratios for rejecting walk,',
     *       /T16,'  Minimum :',F14.5,
     *       /T16,'  Maximum :',F14.5)
 6008 FORMAT(//' (OPTST) INFO: Newton-Raphson used in geometry walk',
     *         ' for state no.',I3/
     *         '         INFO: no absolute guarantee for staying',
     *         ' on same electronic surface'//)
 6010 FORMAT(//' (OPTST) ERROR, ICSTRT = 6 is not',
     *         ' allowed for NEO excited state optimization')
 7010 FORMAT(//' (OPTST) ERROR, ICSTRT = 7 is not',
     *         ' allowed for excited state optimization')
C
C
C
 9000 CONTINUE
      REWIND LUIT1
      IF ( .NOT.FNDLAB(TABLE(6),LUIT1) ) THEN
         CALL READMO(WRK,7,DUMMY,0)
         CALL NEWORB(LBLDAT(2),WRK,.FALSE.)
      END IF
      IF (TRCNO) THEN
C     -- check if already done, e.g. by CISAVE 900720/hjaaj
         REWIND (LUIT1)
         IF (FNDLB2(TABLE(6),LBLINF,LUIT1)) THEN
            IF (LBLINF(2) .EQ. '(CNOORB)') TRCNO = .FALSE.
         END IF
      END IF
      IF (TRCNO) THEN
C
C     -- Transform to natural orbitals --
C
         KCMO   = 1
         KCVECS = KCMO   + NCMOT
         KAOCC  = KCVECS + LROOTS*NCONF
         KW31   = KAOCC  + NASHT
         LW31   = LFREE  - KW31
         IF (LW31 .LE. 0) CALL ERRWRK('OPTST.SIRCNO',-KW31,LFREE)
         CALL READMO(WRK(KCMO),9,DUMMY,0)
         IF (HFCALC) THEN
            WRK(KCVECS) = D1
            ICREF  = 1
            NCVECS = 1
            JSTA   = 1
         ELSE
            REWIND (LUIT1)
C           ... search for "STARTVEC"
            CALL MOLLB2(TABLE(3),LBLINF,LUIT1,LUERR)
            READ (LBLINF(1),'(2I4)',ERR=405) LROIT1,JSTA
            IF (JSTA .LT. 0) THEN
               ICREF  = 1
               NCVECS = 1
               JSTA   = -ISTATE
            ELSE IF (LROOTS .LT. ISTATE) THEN
C              ... to handle .NR ALWAYS calculation following
C                  NEO calculation, e.g. for double core hole.
               ICREF  = 1
               NCVECS = 1
               JSTA   = -ISTATE
               DO 9005 I = 1,ISTATE-1
                  READ (LUIT1)
 9005          CONTINUE
            ELSE
               ICREF  = ISTATE
               NCVECS = LROOTS
               JSTA   = ISTATE
            END IF
            JCVECS = KCVECS
            DO 9010 I = 1,NCVECS
               CALL READT(LUIT1,NCONF,WRK(JCVECS))
               JCVECS = JCVECS + NCONF
 9010       CONTINUE
         END IF
         CALL SIRCNO(CNOLBL,NCVECS,ICREF,WRK(KCVECS),
     *               WRK(KCMO),WRK(KAOCC),INDXCI,WRK(KW31),1,LW31)
C        CALL SIRCNO(KEYCNO,NCVECS,ICREF,CREF,CMO,AOCC,
C    *               INDXCI,WRK,KFRSAV,LFRSAV)
C
         REWIND (LUIT1)
         CALL MOLLAB(TABLE(3),LUIT1,lupri)
C        search for "STARTVEC"
         WRITE (LBLINF(1),'(2I4)') NCVECS,JSTA
         LBLINF(2) = 'CNOSAVE '
         BACKSPACE LUIT1
C        write startvec with new vectors
         WRITE (LUIT1) STAR8,LBLINF,TABLE(3)
         JCVECS = KCVECS
         DO 9020 I = 1,NCVECS
            CALL WRITT(LUIT1,NC4,WRK(JCVECS))
            JCVECS = JCVECS + NCONF
 9020    CONTINUE
C        write the natural orbitals to LUIT1 label NEWORB
         LBLINF(2) = '(CNOORB)'
         WRITE (LUIT1) STAR8,LBLDAT(1),LBLINF(2),TABLE(6)
         CALL WRITT(LUIT1,NCMOT4,WRK(KCMO))
         IF (NASHT .EQ. 1 .OR.
     &      (CNOLBL .NE. 'ONLYFD' .AND. NASHT .GT. 1)) THEN
            WRITE (LUIT1) STAR8,LBLDAT(1),LBLINF(2),TABLE(8)
            NAT4 = MAX(NASHT,4)
            CALL WRITT(LUIT1,NAT4,WRK(KAOCC))
         END IF
         IF (.NOT.FLAG(34)) FLAG(14) = .FALSE.
C        if (int.transf. needed in optimization)
      END IF
C
C *** Append label EODATA, if not already there
C
      REWIND LUIT1
      IF ( .NOT.FNDLAB(TABLE(1),LUIT1) ) THEN
         WRITE (LUIT1) STAR8,LBLDAT,TABLE(1)
      END IF
      REWIND LUIT1
C
C
C *** End of subroutine OPTST.
C
C     Set LROOTS negative to tell SIRNEO to use LROOTS and not NROOTS.
C
      LROOTS = -LROOTS
      CALL QEXIT('OPTST ')
      RETURN
      END
C  /* Deck sirstp */
      SUBROUTINE SIRSTP(ISTEP,CMO,IBNDX,INDXCI,WRK,LFREE)
C
C Written 8-Apr-1984 by Hans Jorgen Aa. Jensen and Hans Agren
C Last revision 14-May-1984 hjaaj
C               21-May-1986 hjaaj (geowlk step control)
C
C Purpose:
C     STEP CONTROL -- check if step was too large, if actual
C   energy deviates too much from the one predicted by last
C   macro iteration.
C     If step is too large, read back reduced L-matrix from
C   last macro iteration and find the level-shift (beta value)
C   giving the step length corresponding to the revised trust
C   radius. Then read the b-vectors from last macro iteration
C   to find the revised eigenvectors and write those on LUIT1
C   (in SIRSAV).
C     Signal to calling program with ISTEP = 0 if step is OK,
C   with ISTEP = 1 if step is to large.
C
C MOTECC-90: The purpose of this module, SIRSTP, and the algorithms
C            used are described in Chapter 8 Section C.5 of MOTECC-90
C            "Step-Control Algorithm"
C
C
C Input:
C  ISTEP, =1    check geometry walk step
C         =-1   restart step (step cannot be rejected,
C               but trust radius will be adjusted)
C         =-2   wave function is converged, only generate info
C         else  normal MC step check
C Output:
C  ISTEP, = 0 if step is OK
C         = 1 if step is too large
C         =-1 if no check because close to convergence
C
C Scratch:
C  WRK
C
#include <implicit.h>
      DIMENSION CMO(*), IBNDX(*), INDXCI(*), WRK(LFREE)
      PARAMETER ( DP1 = 0.1D0, DP5 = 0.5D0,
     *            D0  = 0.0D0, D1  = 1.0D0, D2  = 2.0D0 )
      PARAMETER ( THDE = 1.D-10 )
#include <dummy.h>
C
C Used from common blocks:
C  INFINP: ISTATE,FLAG(*),?
C  INFORB: NCMOT
C  INFVAR: NWOPT
C  INFOPT: NREDL,EMCSCF,EMCOLD,DEPRED,BETA,GAMMA,SHFLVL,
C          STPLEN,STPMAX,STPINC,STPRED,
C          RTRUST,RTTOL,RATMIN,RATGOD,RATREJ,REJWMI,REJWMA, GRDNRM
C  ITINFO: DINFO(*)
C  INFTAP: LUIT1
C  INFPRI: P4FLAG(*),P6FLAG(*)
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <infopt.h>
#include <itinfo.h>
#include <inftap.h>
#include <infpri.h>
#include <gnrinf.h>
C
C -- local:
C
      CHARACTER*8 TABLE(3), RTNLBL(2)
      LOGICAL GEOWLK, NOREJT
C
C
C -- data:
C
      DATA TABLE  /'OLDORB  ','RESTART ','LREDUCED'/
C
C
      CALL QENTER('SIRSTP')
C
C *****************************************************************
C *** Trust region control
C
      GEOWLK = .FALSE.
      NOREJT = .FALSE.
      IF (ISTEP .EQ. 1 .AND. .NOT. OPTNEW) THEN
         GEOWLK = .TRUE.
         WRITE (LUW4,'(//A/)')
     *      ' >>> step control for geometry walk (ABACUS interface)'
      ELSE IF (ISTEP .LT. 0) THEN
         NOREJT = .TRUE.
      END IF
      ISTEP = 0
      DEACT = EMCSCF - EMCOLD
      IF ( ABS(DEPRED/EMCSCF).GT.THDE ) THEN
         RATIO = DEACT / DEPRED
         IF (GEOWLK .OR. P4FLAG(7)) WRITE (LUW4,1100) DEACT,DEPRED,RATIO
      ELSE
         IF (DEPRED .NE. 0.0D0) THEN
           RATIO = DEACT / DEPRED
           WRITE(LUW4,1110) DEACT, DEPRED, RATIO
           WRITE(LUW4,'(A)') ' Close to convergence, ratio set to one.'
         ELSE
           IF (GEOWLK .OR. P4FLAG(7)) WRITE (LUW4,1110) DEACT,DEPRED
         END IF
         RATIO = D1
         ISTEP = -1
      END IF
 1100 FORMAT (/' (SIRSTP) Energy difference;',
     *        /' actual, predicted and ratio:',2F15.10,F12.6)
 1110 FORMAT (/' (SIRSTP) Close to convergence, ratio set to one.',
     *        /' Energy difference; actual and predicted:',1P,2D15.5)
C
C special test if geometry walk
C
      IF (GEOWLK) THEN
         IF (RATIO .LT. REJWMI .OR. RATIO. GT. REJWMA) THEN
            WLKREJ = .TRUE.
            WRITE (LUW4,1210)
            IF (LUPRI .NE. LUW4) WRITE(LUPRI,1210)
            ICI0   = 0
            GOTO 9999
C
C            We don't want to abort the program just because of a long step
C
C            WRITE (LUERR,1210)
C            CALL QTRACE(LUERR)
C            CALL QUIT('***step control:geometry walk step too long***')
         ELSE
            WRITE (LUW4,'(//A/)') ' *** GEOMETRY WALK STEP ACCEPTED.'
         END IF
         EMCOLD = EMCSCF
         ICI0   = 0
         GO TO 9999
      END IF
 1210 FORMAT(//' *** GEOMETRY WALK STEP TOO LONG,',
     *        /'     DECREASE WALK TRUST RADIUS AND TRY AGAIN')
C
C save information for final summary print-out
C
      DINFO(4) = DEPRED
      DINFO(5) = DEACT
      DINFO(6) = RATIO
C
C
C     Normal MCSCF optimization (not first iteration in geometry walk)
C
      RATVGD = D1 - DP1*(D1 - RATGOD)
      RTSAVE = RTRUST
      IF (ISTATE .EQ. 1 .AND. .NOT. FLAG(35)) THEN
C     (case 1: ground state optimization)
C      flag(35): tight step control also for ground state optimization
         IF (RATIO .LT. RATMIN) THEN
            RTRUST = STPRED*RTRUST
            IF (RATIO .LT. RATREJ .AND. .NOT. NOREJT) THEN
               ISTEP  = 1
               RTRUST = MIN(RTRUST, STPRED*STPLEN)
            END IF
         ELSE IF (RATIO.GT.RATGOD) THEN
            RTRUST = MIN(STPMAX,STPINC*RTRUST)
            IF (RATIO.GT.RATVGD) RTRUST = MIN(STPMAX,STPINC*RTRUST)
C           extra increase if ratio is very good
         END IF
      ELSE
C     (case 2: excited state optimization)
         IF (RATIO .LT. RATMIN .OR. RATIO .GT. (D2-RATMIN)) THEN
            RTRUST = STPRED*RTRUST
            IF (.NOT. NOREJT) THEN
            IF (RATIO .LT. RATREJ .OR. RATIO .GT. (D2-RATREJ)) THEN
               ISTEP  = 1
               RTRUST = MIN(RTRUST, STPRED*STPLEN)
            END IF
            END IF
         ELSE IF (RATIO .GT. RATGOD .AND. RATIO .LT. (D2-RATGOD)) THEN
            RTRUST = MIN(STPMAX,STPINC*RTRUST)
            IF (RATIO .GT. RATVGD .AND. RATIO .LT. (D2-RATVGD))
     *         RTRUST = MIN(STPMAX,STPINC*RTRUST)
C           extra increase if ratio is very good
         END IF
      END IF
C
C *****************************************************************
C *** Step is acceptable ( or restart ) :
C
      IF (ISTEP.LE.0) THEN
         EMCOLD = EMCSCF
         GO TO 9000
      END IF
C
C *****************************************************************
C *** Step is too large :
C
      IF (.NOT.P4FLAG(7)) WRITE (LUW4,1100) DEACT,DEPRED,RATIO
      IF (LUPRI.NE.LUW4)   WRITE (LUPRI,1100) DEACT,DEPRED,RATIO
      WRITE (LUW4,7000)
      IF (LUPRI.NE.LUW4 .AND. .NOT.P6FLAG(10)) WRITE (LUPRI,7000)
 7000 FORMAT(/' (SIRSTP) step is too large -- backstep')
C
C     Core allocation:
C
      NNREDL = NREDL*(NREDL+1)/2
      KREDL  = 1
      KEVEC  = KREDL  + NNREDL
      KEVAL  = KEVEC  + NREDL*NREDL
      KWRK11 = KEVAL  + NREDL
      KWRK12 = KWRK11 + NNREDL
      KWRK13 = KWRK12 + NREDL
      LWRK11 = LFREE  + 1 - KWRK11
      KXKAP  = KWRK11
      KWRK21 = KXKAP  + NWOPT
      LWRK21 = LFREE  + 1 - KWRK21
C
C     Recover old orbitals, IBNDX and REDL from previous macro iteration
C
      REWIND LUIT1
      CALL MOLLAB(TABLE(1),LUIT1,lupri)
      CALL READT(LUIT1,NCMOT,CMO)
      CALL MOLLB2 (TABLE(2),RTNLBL,LUIT1,LUERR)
      IF (RTNLBL(2) .NE. 'NEOSAVE ') THEN
         WRITE (LUW4,'(/A/2A/A)')
     *   ' Sorry, backstep only implemented for NEO',
     *   ' Label from last iterations : ',RTNLBL(2),
     *   ' Program cannot continue (hint: try ".NEO ALWAYS" option).'
         CALL QTRACE(LUW4)
         CALL QUIT('Sorry, backstep only implemented for NEO')
      END IF
      READ (LUIT1) DUM,DUM,DUM,DUM,DUM,IDUM,
     *             MREDL,(IBNDX(I),I=1,MREDL)
      IF (MREDL.NE.NREDL) THEN
C        this should never occur ...
         WRITE (LUERR,'(/1X,A,2I5)')
     *      'SIRSTP: NREDL on LUIT1 inconsistent with NREDL in common',
     *      MREDL,NREDL
         CALL QTRACE(LUERR)
         CALL QUIT('ERROR (SIRSTP) unexpected value of NREDL on LUIT1')
      END IF
      CALL MOLLAB(TABLE(3),LUIT1,lupri)
      CALL READT(LUIT1,NNREDL,WRK(KREDL))
C
C
      CALL NEORED (2,0,0,WRK(KWRK11),DUMMY,BETA,DUMMY,IBNDX,
     *             NREDL,WRK(KREDL),WRK(KEVAL),WRK(KEVEC),LWRK11)
C     CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,G,IBNDX,
C    *             NREDL,REDL,EVAL,EVEC, LBVECS)
C
C
      CALL UPDBET(WRK(KREDL),WRK(KEVEC),WRK(KEVAL),WRK(KWRK11),LWRK11)
C     CALL UPDBET(REDL,EVEC,EVAL,WRK,LFREE)
      IF (P6FLAG(10)) THEN
         WRITE (LUPRI,8010) BETA,RTRUST
         CALL OUTPAK(WRK(KREDL),NREDL,1,LUPRI)
         WRITE (LUPRI,8012)
         WRITE (LUPRI,8015) (WRK(KEVAL-1+I),I=1,NREDL)
         LIMP = MIN(((ISTATE+5)/4) * 4, NREDL)
         CALL OUTPUT(WRK(KEVEC),1,NREDL,1,LIMP,NREDL,NREDL,1,LUPRI)
      END IF
 8010 FORMAT(/' (SIRSTP) step was too long;',
     *       /T11,'new beta and trust radius: ',2F15.8,
     *       /' - The updated reduced L matrix:')
 8012 FORMAT(/' - Eigenvalues and eigenvectors of ',
     *        'the updated reduced L matrix:')
 8015 FORMAT(/,(10X,1P,4D15.6))
C
C *** Save the revised (restricted step) eigenvectors (CI part only)
C     and rotate the orbitals corresponding to the revised eigenvector.
C
C
      SHFLVL = WRK(KEVAL-1+ISTATE)
      CALL SIRSAV ('NEOSAVE',CMO,IBNDX,WRK(KREDL),WRK(KEVEC),
     *             WRK(KXKAP),INDXCI,WRK(KWRK21),LWRK21)
C
C     CALL SIRSAV (KEYWRD,CMO,IBNDX,REDL,EVEC,XKAP,INDXCI,WRK,LFREE)
C
C === Revise predicted energy difference
C     between next macro iteration and this one.
C
        DEPRED = ( GAMMA*(D1-GAMMA) + DP5*GAMMA*GAMMA*
     *           (WRK(KEVEC+(ISTATE-1)*NREDL) ** (-2)) ) *
     *           WRK(KEVAL-1+ISTATE) / (BETA*BETA)
C
C
C *****************************************************************
C *** End of subroutine SIRSTP
C
 9000 CONTINUE
      IF (P4FLAG(7)) WRITE (LUW4,1200) RTSAVE,RTRUST
 1200 FORMAT(/' (SIRSTP) Old and new MC trust radius:',2F15.10)
 9999 CALL QEXIT('SIRSTP')
      RETURN
      END
C  /* Deck timopt */
      SUBROUTINE TIMOPT(KEY,LUPRI,WRK,LFREE)
C
C Written 18-Dec-1984 Hans Joergen Aa. Jensen.
C Revised 891201-hjaaj: LUPRI in parameter list.
C
C Purpose:
C  Write timing statistics to LUPRI for central (time-consuming)
C  routines used in MC optimization, that is from SIROPT.
C
#include <implicit.h>
      CHARACTER*(*) KEY
C
      PARAMETER ( NTYP = 7 )
      DIMENSION WRK(7,2,NTYP)
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
C
C
C Used from common blocks:
C   INFTIM : NCALLS,TIMCPU,TIMWAL
C
#include <inftim.h>
C
C Local variables:
C
      INTEGER INAME(NTYP),IEND(NTYP)
      CHARACTER*8 NAME(7,NTYP)
      DATA INAME /1,2,6,7,10,8,3/, IEND /7,4,5,5,5,5,7/
      DATA NAME /'    GRAD','   MAKDM','  FCKMAT',
     *           '   GETH2','  ORDIAG','  ORBGRD','  CIGRAD',
     *           '  CISIGC','    diag','  ONEELC',
     *           '  TWOELC','        ','        ','        ',
     *           '  LINTRN','  SIRTR1','  ORBSIG',
     *           '   CISIG','  SOLLIN','        ','        ',
     *           '  MAKTDM','  "diag"',' "oneel"',
     *           ' "twoel"','"sym PV"','        ','        ',
     *           '  CISIGO','  DIAG1X','  ONEELX',
     *           '  TWOELX','   other','        ','        ',
     *           '  ORBLIN','  SIRTR1','  ORBSIG',
     *           '   CISIG','  SOLLIN','        ','        ',
     *           ' UPDGRAD','   MAKDM','  FCKMAT',
     *           '   GETH2','  ORDIAG','  ORBGRD','  CIGRAD'/
C
C
      IF (KEY(1:3) .EQ. 'STA') THEN
         DO 140 ITYP = 1,20
            NCALLS(ITYP) = 0
            NVECS (ITYP) = 0
            DO 120 IBLCK = 1,7
               TIMCPU(IBLCK,ITYP) = D0
               TIMWAL(IBLCK,ITYP) = D0
  120       CONTINUE
  140    CONTINUE
         RETURN
      END IF
C
      WRITE (LUPRI,10)
      WRITE (LUPRI,15)
C
      DO 400 I = 1,NTYP
         II     = INAME(I)
         JEND   = IEND(I)
         MCALLS = NCALLS(II)
         IF (MCALLS .EQ. 0) THEN
            WRITE (LUPRI,25) NAME(1,I)
            WRITE (LUPRI,15)
            GO TO 400
         END IF
         IF (NVECS(II) .GT. 0) THEN
            FAC = NVECS(II)
         ELSE
            FAC = MCALLS
         END IF
         FAC = D1 / FAC
         DO 300 J = 1,JEND
            WRK(J,1,I) = FAC * TIMCPU(J,II)
            WRK(J,2,I) = FAC * TIMWAL(J,II)
  300    CONTINUE
C
         IF (NVECS(II) .LE. 0) THEN
            WRITE (LUPRI,20) NAME(1,I),MCALLS,TIMCPU(1,II),WRK(1,1,I),
     *                       TIMWAL(1,II),WRK(1,2,I)
         ELSE
            WRITE (LUPRI,21) NAME(1,I),MCALLS,NVECS(II),
     *                       TIMCPU(1,II),WRK(1,1,I),
     *                       TIMWAL(1,II),WRK(1,2,I)
         END IF
         IF (JEND.GT.1)
     *      WRITE (LUPRI,30) (NAME(J,I),TIMCPU(J,II),WRK(J,1,I),
     *                       TIMWAL(J,II),WRK(J,2,I),J=2,JEND)
         WRITE (LUPRI,15)
  400 CONTINUE
C
      RETURN
C
   10 FORMAT(//' (TIMOPT) Timing of MC optimization:',
     *       //'    Name . # calls .    total CPU .   CPU / call .',
     *         '   total wall .  wall / call .')
   15 FORMAT(1X,78('-'))
   20 FORMAT(/1X,A8,I10,4F15.2)
   21 FORMAT(/1X,A8,I10/'   # vecs',I10,4F15.2)
   25 FORMAT(/1X,A8,' has not been called in this calculation.')
   30 FORMAT(/1X,A8,10X,4F15.2)
      END
C  /* Deck wr_sifc */
      SUBROUTINE WR_SIRIFC(CMO,DV,PV,FC,FV,FQ,CREF,FCAC,H2AC,WRK,LFREE,
     *                  GORB,ERLM,ORBHES,XINDX)
C
C MOTECC: this routine writes interface file for post-programs.
C         The structure of the file is described here.
C
C Written 14-Feb-1985 Hans Joergen Aa. Jensen
C Revisions:
C  21-Nov-1985 hjaaj (for RHF)
C   8-Feb-1987 hjaaj (full Fock matrix)
C     Oct-1988 hjaaj (for non-guga mcscf)
C   6-Sep-1989 hjaaj (corrected FCDIA calculation)
C   2-Sep-1999 hjaaj (do not use FV for NASHT.eq.0)
C
C Purpose:
C  Write information to LUSIFC needed for (1) Trygve Helgaker's
C  first and second order geometry derivatives program and/or
C  (2) response calculation.
C
C Suggestions:
C  890906-hjaaj: save parameter telling if orbitals are canonical HF
C                orbitals or natural orbitals or ?
C
C  The following records are written:
C
C    0) label LBSIFC ("SIR IPH ")
C    1) POTNUC,EMY,EACTIV,EMCSCF
C    2) NISHT,NASHT,...
C    3) CMO
C    4) CREF
C    5) DV
C    6) FOCK
C    7) PV
C    8) FC
C    9) FV
C   10) FCAC
C   11) H2AC
C   12) GORB
C   If (GUGA) then
C      13) label "CIDIAG1 "
C      14) CIDIAG1
C   end if
C   15) label "CIDIAG2 "
C   16) CIDIAG2
C   17) label "ORBDIAG "
C   18) ORBDIAG
C
C    *) label "SIRFLAGS"
C    *) FLAG(1:NFLAG)
C    *) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C   If (fields) then
C    *) label "EXTFIELD"
C    *) NFIELD
C    *) (EFIELD(I),I=1,NFIEL4)
C    *) (LFIELD(I),I=1,NFIEL4)
C   end if
C   If (solvent) then
C    *) label "SOLVINFO"
C    *) EPSOL,EPSTAT,RSOL(1:3),LSOLMX
C    *) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C    *) ERLM(LM,1), LM = 1,NLMSOL)
C    *) (TRLM(LM), LM = 1,NLMSOL)    where TRLM(i) = ERLM(i,2)
C    *) NSYM, NBAS
C    *) label "SOLVTMAT"
C    *) TMAT(1:NNORBX)
C   end if
C   If (CSF's) then
C    *) label "CREFDETS"
C    *) CREF in dets
C   end if
C    *) label 'TRCCINT"
C    *) NSYM, NORBT, ...
C    *) FCDIA, ISMO
C    *) CMO
C
#include <implicit.h>
      DIMENSION CMO(*),DV(*),PV(*),FC(*),FV(*),FQ(NORBT,*),XINDX(*)
      DIMENSION ERLM(NLMSOL,2),GORB(*),CREF(*),FCAC(*),H2AC(*),WRK(*)
C
      PARAMETER ( D0=0.0D0, D2=2.0D0 )
#include <dummy.h>
C
C Used from common blocks:
C   INFINP : ABAIPH,ISTATE,ISPIN,NACTEL,FLAG(*),INERSI
C   INFOPT : EMY, EACTIV, EMCSCF, ESOLT
C   INFVAR : NCONF,NWOPT
C   INFORB : NSYM,NISHT,NASHT,NNASHX,NNASHY,NOCCT,NNOCCX,NCMOT,...
C   INFIND : IROW(*)
C   INFTAP : LUSIFC
C   INFDIM : NASHDI,NASHMA,NORBMA
C   INFPRI : ?
C CBGETDIS : IADH2
C   CBIREA : LMULBS
C   R12INT : R12CAL
C
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <priunit.h>
#include <pcm.h>
#include <pcmlog.h>
#include <infinp.h>
#include <infopt.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infdim.h>
#include <infpri.h>
#include <cbgetdis.h>
#include <cbirea.h>
#include <r12int.h>
C----------------
C CBN+JK 03.01.06 gnrinf.h included
#include <gnrinf.h>
C----------------
C
      LOGICAL     FNDLAB, ORBHES, LGLOCAL, EXP1VL, TOFILE
C
      CHARACTER*8 LAB123(3), TABLE(10), LABEOD, LABCC, LABABA
      CHARACTER*8 LABAUX, LABUNI
      DATA LAB123/'********','********','********'/
      DATA TABLE /'CIDIAG1 ','CIDIAG2 ','ORBDIAG ','SIRFLAGS',
     *            'EXTFIELD','SOLVINFO','SOLVTMAT','CREFDETS',
     *            'PCMINFO ','PCMJXMAT'/
      DATA LABEOD, LABCC  /'EODATA  ', 'TRCCINT '/
      DATA LABUNI, LABAUX /'FULLBAS ', 'AUXILBAS'/
C
      CALL QENTER('WR_SIRIFC')
      LGLOCAL = BOYORB.OR.PIPORB
      CALL GETDAT(LAB123(2),LAB123(3))
C     ... place date in LAB123(2) and time in LAB123(3)

C
C *** Now create LUSIFC and write ...
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,FNSIFC,
     &            'UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND (LUSIFC)

      IF (IPRI6 .GE. 5) WRITE (LUPRI,'(/A)')
     *   ' >>> Writing Sirius interface file for'//
     *   ' post-programs in WR_SIRIFC. <<<'
C
C     Calculate Fock matrix
C
      KFCDIA = 1
      KFOCK  = KFCDIA + NORBT
      KWRK1  = KFOCK  + N2ORBT
      KFCI   = KWRK1
      KUDV   = KFCI   + NORBMA*NORBMA
      KFCFV  = KUDV   + N2ASHX
      KEND   = KFCFV  + NNORBT
      KEND   = MAX(KEND,KWRK1+NCONF,KWRK1+NWOPT)
      IF (KEND .GT. LFREE) CALL ERRWRK('WR_SIRIFC',KEND,LFREE)
C
      IF (NASHT .EQ. 0) THEN
         DO I = 1,NNORBT
            WRK((KFCFV-1)+I) = D2*FC(I)
         END DO
      ELSE
         CALL DSPTSI(NASHT,DV,WRK(KUDV))
         DO I = 1,NNORBT
            WRK((KFCFV-1)+I) = D2*(FC(I) + FV(I))
         END DO
      END IF
C
      CALL DZERO(WRK(KFOCK),N2ORBT)
      DO 500 ISYM = 1,NSYM
         NORBI  = NORB(ISYM)
         IORBI  = IORB(ISYM)
C
         JFCDIA = KFCDIA - 1 + IORBI
         DO 170 I = 1,NORBI
            WRK(JFCDIA+I) = FC( IIORB(ISYM) + IROW(I+1) )
  170    CONTINUE
C
         JFOCK  = KFOCK + I2ORB(ISYM)
         NOCCI  = NOCC(ISYM)
      IF (NOCCI .EQ. 0) THEN
         IF (NORBI .GT. 0) CALL DZERO(WRK(JFOCK),N2ORB(ISYM))
      ELSE
         NISHI  = NISH(ISYM)
         NASHI  = NASH(ISYM)
         NSSHI  = NSSH(ISYM)
C        inactive-general + junk in active-general and secondary-general
         IF (NISHI .GT. 0) THEN
            JFCFV  = KFCFV + IIORB(ISYM)
            CALL DSPTSI(NORBI,WRK(JFCFV),WRK(JFOCK))
         END IF
C        active-general
         IF (NASHI .GT. 0) THEN
            CALL DSPTSI(NORBI,FC(1+IIORB(ISYM)),WRK(KFCI))
            IASHI  = IASH(ISYM)
            JUDV   = KUDV + IASHI + IASHI*NASHT
            JFCI   = KFCI + NISHI
            IFOCKJ = JFOCK + NISHI
            CALL DGEMM('N','N',NASHI,NORBI,NASHI,1.D0,
     &                 WRK(JUDV),NASHT,
     &                 WRK(JFCI),NORBI,0.D0,
     &                 WRK(IFOCKJ),NORBI)
            DO 360 J = 1,NORBI
               IFOCKJ = JFOCK - 1 + NISHI + (J-1)*NORBI
               DO 340 NI = 1,NASHI
                  WRK(IFOCKJ + NI) = WRK(IFOCKJ + NI)
     &                             + FQ(IORBI + J, IASHI + NI)
  340          CONTINUE
  360       CONTINUE
         END IF
C        secondary-general is zero
         IF (NSSHI .GT. 0) THEN
            DO 460 J = 1,NORBI
               IFOCKJ = JFOCK - 1 + NOCCI + (J-1)*NORBI
               DO 440 NI = 1,NSSHI
                  WRK(IFOCKJ + NI) = D0
  440          CONTINUE
  460       CONTINUE
         END IF
      END IF
C
  500 CONTINUE
C
      IF (P6FLAG(14)) THEN
         WRITE (LUPRI,'(/A)') ' Full Fock matrix from WR_SIRIFC :'
         CALL OUTPTB(WRK(KFOCK),NORB,NSYM,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' FC diagonal elememts from WR_SIRIFC :'
         WRITE (LUPRI,'(/,(1X,5F15.8))') (WRK(KFCDIA-1+I),I=1,NORBT)
      END IF
C
C
      WRITE (LUSIFC) LAB123,LBSIFC
C     ... the label LBSIFC is initialized in block data sbdtap
      WRITE (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,
     &               ISTATE,ISPIN,NACTEL,LSYM
      WRITE (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,
     &               NCDETS, NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,
     &               NSYM,MULD2H, NRHF,NFRO,
     &               NISH,NASH,NORB,NBAS,
     &               NELMN1, NELMX1, NELMN3, NELMX3, MCTYPE,
     &               NAS1, NAS2, NAS3
C
C 880920-hjaaj: later write label here for orbitals
C
      NC4    = MAX(4,NCONF)
      NW4    = MAX(4,NWOPT)
      NWH4   = MAX(4,NWOPH)
      NCMOT4 = MAX(4,NCMOT)
      MMORBT = MAX(4,NNORBT)
      M2ORBT = MAX(4,N2ORBT)
      MMASHX = MAX(4,NNASHX)
      MMASHY = MAX(4,NNASHY)
      M2ASHY = MAX(4,NNASHX*NNASHX)
C
      CALL WRITT (LUSIFC,NCMOT4,CMO) ! rec no. 3
      CALL WRITT (LUSIFC,NC4,CREF)   ! rec no. 4
      CALL WRITT (LUSIFC,MMASHX,DV)  ! rec no. 5
C Write Fock matrix.
      CALL WRITT (LUSIFC,M2ORBT,WRK(KFOCK)) ! rec no. 6
      CALL WRITT (LUSIFC,M2ASHY,PV)  ! rec no. 7
      CALL WRITT (LUSIFC,MMORBT,FC)  ! rec no. 8
      IF (NASHT .EQ. 0) THEN
         CALL DZERO(WRK(KFOCK),NNORBT)
         CALL WRITT (LUSIFC,MMORBT,WRK(KFOCK))
      ELSE
         CALL WRITT (LUSIFC,MMORBT,FV) ! rec. no. 9
      END IF
C
      CALL WRITT (LUSIFC,MMASHX,FCAC) ! rec. no. 10
      IF (IADH2 .GE. 0) THEN
#if defined (VAR_NEWCODE)
... 27-Jun-1987/hjaaj -- noget lignende : MAERKE
         DO 100 IJ = 1,NNASHX
            CALL READDI(LUDA2,IADH2+IJ,IRAT*NNASHX,H2AC)
            CALL WRITT (LUSIFC,MMASHX,H2AC)
  100    CONTINUE
... men der mangler en label som fortaeller response at H2AC on disk.
#else
         WRITE (LUPRI,*)
     &      'SIRCI.WR_SIRIFC: ".DISKH2" not implemented here yet.'
         CALL QTRACE(LUPRI)
         CALL QUIT('SIRCI.WR_SIRIFC: ".DISKH2" not implemented here.')
#endif
      ELSE
         CALL WRITT (LUSIFC,(MMASHX*MMASHX),H2AC) ! rec no. 11
      END IF
      CALL WRITT (LUSIFC,NWH4,GORB) ! rec no. 12
C
C *** now write diagonals of L-matrix needed for LINTRN
C     and for conjugate gradient "next" algorithm
C
      IF (NCONF .GT. 1) THEN
         REWIND LUIT2
         IF ( FNDLAB(TABLE(1),LUIT2) ) THEN
C        ... CIDIAG1 is only used if CASGUGA
            CALL READT (LUIT2,NCONF,WRK(KWRK1))
            WRITE (LUSIFC) LAB123,TABLE(1)
            CALL WRITT (LUSIFC,NC4,WRK(KWRK1))
         ELSE
            REWIND LUIT2
         END IF
C
         CALL MOLLAB(TABLE(2),LUIT2,lupri)
         CALL READT (LUIT2,NCONF,WRK(KWRK1))
         WRITE (LUSIFC) LAB123,TABLE(2)
         CALL WRITT (LUSIFC,NC4,WRK(KWRK1))
      END IF
C
      IF (NWOPT .GT. 0 .AND. ORBHES) THEN
         REWIND LUIT2
         CALL MOLLAB(TABLE(3),LUIT2,lupri)
         CALL READT (LUIT2,NWOPT,WRK(KWRK1))
         WRITE (LUSIFC) LAB123,TABLE(3)
         CALL WRITT (LUSIFC,NW4,WRK(KWRK1))
      END IF
      WRITE (LUSIFC) LAB123,TABLE(4)
      WRITE (LUSIFC) (FLAG(I),I=1,NFLAG)
      WRITE (LUSIFC) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
C     Write External fields, if any
      IF (NFIELD .GT. 0) THEN
         WRITE (LUSIFC) LAB123,TABLE(5)
         WRITE (LUSIFC) NFIELD,DUMMY,DUMMY,DUMMY,DUMMY
         NFIEL4 = MAX(4,NFIELD)
         WRITE (LUSIFC) (EFIELD(I),I=1,NFIEL4)
         WRITE (LUSIFC) (LFIELD(I),I=1,NFIEL4)
      END IF

C     Write Solvent information, if Solvent calculation.
      IF (FLAG(16)) THEN
         KTMAT = KWRK1
         KWRK2 = KTMAT + NNORBT
         LWRK2 = LFREE - KWRK2 + 1
         IF (KWRK2 .GT. LFREE)
     &      CALL ERRWRK('WR_SIRIFC-SOLVENT',KWRK2,LFREE)
         WRITE (LUSIFC) LAB123,TABLE(6)
         WRITE (LUSIFC) EPSOL,EPSTAT,RSOL,LSOLMX
         WRITE (LUSIFC) GRDNRM,POTNUC,EMY,EACTIV,ESOLT,EMCSCF
         NLMSO4 = MAX(4,NLMSOL)
         CALL WRITT (LUSIFC,NLMSO4,ERLM(1,1))
         CALL WRITT (LUSIFC,NLMSO4,ERLM(1,2))
         WRITE (LUSIFC) NSYM, NBAS
         CALL SOLFCK(DUMMY,DUMMY,WRK(KTMAT),CMO,ERLM(1,1),.TRUE.,
     &               DUMSOL,WRK(KWRK2),LWRK2,IPRSOL)
         WRITE (LUSIFC) LAB123,TABLE(7)
         CALL WRITT (LUSIFC,MMORBT,WRK(KTMAT))
      END IF
C     Write PCM solvent info if PCM calculation
      IF (PCM) THEN
         KTMAT = 1
         KDCAO = KTMAT + NNORBT
         KDVAO = KDCAO + N2BASX
         KPOT  = KDVAO + N2BASX
         KWRK2 = KPOT  + NTS
         LWRK2 = LFREE - KWRK2 + 1
         IF (KWRK2 .GT. LFREE) CALL ERRWRK('WR_SIRIFC-PCM',KWRK2,LFREE)
C
         CALL PCMFCK(DUMMY,DUMMY,WRK(KTMAT),CMO,.TRUE.,
     *               DUMSOL,WRK(KWRK2),LWRK2,IPRSOL)
         WRITE (LUSIFC) LAB123,TABLE(10)
         CALL WRITT (LUSIFC,MMORBT,WRK(KTMAT))
         WRITE (LUSIFC) LAB123,TABLE(9)
         WRITE (LUSIFC) NTS
         WRITE (LUSIFC) (QSE(ITS),ITS = 1, NTS)
         WRITE (LUSIFC) (QSN(ITS),ITS = 1, NTS)
         WRITE (LUSIFC) (-QSE(ITS)-QSN(ITS),ITS = 1, NTS)
C
Clf         WRITE (LUSIFC) (QSE(ITS)-QSN(ITS),ITS = 1, NTS)
Clf
C
C     Calculate expectation value of the electronic potential
C     at the tessera
C
         EXP1VL = .TRUE.
         TOFILE = .FALSE.
         CALL DZERO(WRK(KDCAO),2*N2BASX)
         CALL FCKDEN((NISHT.GT.0),(NASHT.GT.0),WRK(KDCAO),
     &               WRK(KDVAO),CMO,DV,WRK(KWRK2),LWRK2)
         CALL DGEFSP(NBAST,WRK(KDCAO),WRK(KWRK2))
         CALL DCOPY(NNORBX,WRK(KWRK2),1,WRK(KDCAO),1)
ckr         CALL PKSYM1(WRK(KWRK2),WRK(KDCAO),NBAS,NSYM,1)
         IF (NASHT .GT. 0) THEN
            CALL DGEFSP(NBAST,WRK(KDVAO),WRK(KWRK2))
            CALL DCOPY(NNORBX,WRK(KWRK2),1,WRK(KDVAO),1)
ckr            CALL PKSYM1(WRK(KWRK2),WRK(KDVAO),NBAS,NSYM,1)
         END IF
Ckr
Ckr      Check whether DV+DC is good enough for open-shell systems (should
Ckr      be the case).
Ckr
         CALL DAXPY(NNBASX,1.0D0,WRK(KDVAO),1,WRK(KDCAO),1)
         CALL J1INT(WRK(KPOT),EXP1VL,WRK(KDCAO),1,TOFILE,'NPETES ',
     &              1,WRK(KWRK2),LWRK2)
         WRITE (LUSIFC) (WRK(KPOT + I - 1), I = 1, NTS)
         IF (NEQRSP) CALL V2Q(WRK(KWRK2),WRK(KPOT),QSENEQ,QETNEQ,NEQRSP)
      END IF
C--------------
C CBN+JK, 03.01.06
C-------------------
C Write qm3 solvent part of the Fock matrix to interface file

      IF (QM3) THEN
         KTMAT = KWRK1
         KFSOLAO = KTMAT + NNORBT
         KUCMO = KFSOLAO + NNBASX
         KWRK2 = KUCMO + NORBT*NBAST
         LWRK2 = LFREE - KWRK2 + 1
         IF (KWRK2 .GT. LFREE) CALL ERRWRK('WR_SIRIFC-QM3',KWRK2,LFREE)
         CALL QM3_FMO(WRK(KFSOLAO),WRK(KWRK2),LWRK2,IPRSOL)
         CALL UPKCMO(CMO,WRK(KUCMO))
         CALL UTHU(WRK(KFSOLAO),WRK(KTMAT),WRK(KUCMO),WRK(KWRK2),
     &             NBAST,NORBT)
         WRITE (LUSIFC) LAB123,TABLE(7)
         CALL WRITT (LUSIFC,MMORBT,WRK(KTMAT))
      END IF
C--------------
C CBN+JK, 03.01.06
C-------------------
C
C     Write CREF in determinants here, if NCDETS .gt. NCONF
C
      IF (NCDETS .GT. NCONF) THEN
         IF (FLAG(27)) THEN
            WRITE (LUPRI,*)
     &      'WR_SIRIFC ERROR, .DETERMINANTS but NCDETS.gt.NCONF:',
     &      NCDETS, NCONF
            CALL QUIT('WR_SIRIFC: NCDETS.gt. NCONF for .DETERMINANTS')
         END IF
         WRITE (LUSIFC) LAB123,TABLE(8) ! label "CREFDETS"
         KWRK1 = 1 + NCDETS
         LWRK1 = LFREE - KWRK1
         CALL GETDETS(LSYM,NCONF,NCDETS,XINDX,CREF,WRK,WRK(KWRK1),LWRK1)
         NC4 = MAX(4,NCDETS)
         CALL WRITT(LUSIFC,NC4,WRK)
      END IF
C
C     Interface for CC modules
C
      IF (LMULBS) THEN
         IOFF1 = 1
         IOFF2 = 1 + NORBT
         NBAST1 = 0
         NCMOT1 = 0
         DO ISYM = 1, NSYM
            NBAST1 = NBAST1 + MBAS1(ISYM)
            NCMOT1 = NCMOT1 + MBAS1(ISYM) * NORB(ISYM)
            DO I = 1, NORB(ISYM)
               CALL DCOPY(MBAS1(ISYM),CMO(IOFF1),1,WRK(IOFF2),1)
               IOFF1 = IOFF1 + NBAS(ISYM)
               IOFF2 = IOFF2 + MBAS1(ISYM)
            ENDDO
         ENDDO
         NCMOTX = MAX(4,NCMOT1)
         WRITE (LUSIFC) LAB123, LABCC ! label "TRCCINT "
         WRITE (LUSIFC) NSYM,NORBT,NBAST1,NCMOT1,(NOCC(I),I=1,NSYM),
     *        (NORB(I),I=1,NSYM),(MBAS1(I),I=1,NSYM),POTNUC,EMCSCF
         WRITE (LUSIFC) (WRK(KFCDIA-1+I),I=1,NORBT),(ISMO(I),I=1,NORBT),
     *                  DUMMY,DUMMY,DUMMY
         WRITE (LUSIFC) (WRK(I),I=1+NORBT,NCMOTX+NORBT)
      ELSE
         WRITE (LUSIFC) LAB123, LABCC
         WRITE (LUSIFC) NSYM,NORBT,NBAST,NCMOT,(NOCC(I),I=1,NSYM),
     *        (NORB(I),I=1,NSYM),(NBAS(I),I=1,NSYM),POTNUC,EMCSCF
         WRITE (LUSIFC) (WRK(KFCDIA-1+I),I=1,NORBT),(ISMO(I),I=1,NORBT),
     *                  DUMMY,DUMMY,DUMMY
         WRITE (LUSIFC) (CMO(I),I=1,NCMOT4)
      END IF
C
      IF (LMULBS .OR. R12CAL) THEN
C        Interface for auxiliary basis set (WK/UniKA/04-11-2002).
         DO ISYM = 1, NSYM
           IF (NRXR12(ISYM).GT.(NORB1(ISYM)-NOCC(ISYM))) THEN
             WRITE(LUPRI,*) 'You specified more virtual R12 orbitals '//
     &                   'than there are virtual orbitals in this '//
     &                   'symmery:'
             WRITE(LUPRI,*) 'NVIR(',ISYM,') = ',NORB1(ISYM)-NOCC(ISYM)
             WRITE(LUPRI,*) 'NRXR12(',ISYM,') = ',NRXR12(ISYM)
             CALL QUIT('Too many virtual R12 orbitals')
           ELSE
             !redefine NOCC
             NOCC(ISYM) = NOCC(ISYM) + NRXR12(ISYM)
           END IF
         END DO
         LUMULB = 34
         CALL GPOPEN(LUMULB,'AUXBAS','UNKNOWN','SEQUENTIAL',
     &               'FORMATTED',IDUMMY,LDUMMY)
         READ (LUMULB,*) NDUMM
         NAUXT = 0
         MOAUX = 0
         NFULLT = 0
         NMOFLT = 0
         NOCCT = 0
         NORBFT = 0
         DO ISYM = 1, NSYM
            NOCCT  = NOCCT + NOCC(ISYM)
            NORBFT = NORBFT + NORB1(ISYM) + NORB2(ISYM)
         END DO
         NBASC = NBAST + NOCCT
         NORBFT = NORBFT + NOCCT
         IOFFEV = KFCDIA + NBASC
         IOFFMO = IOFFEV + NBASC * 2
         DO ISYM = 1, NSYM
            READ (LUMULB,*) IDUMM,NORBI,NAUXI,NBASI
            NAUXT = NAUXT + NAUXI
            MOAUX = MOAUX + NAUXI * NBASI
            NFULLT = NFULLT + NBASI
            NMOFLT = NMOFLT + NBASI * (NORBI + NAUXI + NOCC(ISYM))
            DO I = 1, NAUXI
               READ (LUMULB,*) NDUMM
               READ (LUMULB,'(4E30.20)') WRK(IOFFEV)
               IOFFEV = IOFFEV + 1
               READ (LUMULB,'(4E30.20)') (WRK(IOFFMO+K), K = 0, NBASI-1)
               IOFFMO = IOFFMO + NBASI
            END DO
         END DO
         CALL GPCLOSE(LUMULB,'KEEP')
         WRITE (LUSIFC) LAB123, LABAUX ! label "AUXILBAS"
         WRITE (LUSIFC) NSYM,NAUXT,NBAST,MOAUX,
     *                  (0,I=1,NSYM),
     *                  (NORB2(I),I=1,NSYM),
     *                  (NBAS(I),I=1,NSYM),POTNUC,EMCSCF
         WRITE (LUSIFC) (WRK(NBASC+I),I=1,MAX(NAUXT,4))
         WRITE (LUSIFC) (WRK(3*NBASC+I),I=1,MAX(MOAUX,4))
C        Interface for primary and secondary basis sets (WK/UniKA/04-11-2002).
         IF (LGLOCAL) THEN
            NOCL = NOCC(1)
            NBASL = NBAS(1)
            NORB1L = NORB1(1)
            NORB2L = NORB2(1)
            NORBFTL = NORB1L + NORB2L + NOCL
            NMOFLTL = NORBFTL * NBASL 
            IEIGL = IOFFMO
            IOCCL = IEIGL + NORBFTL
            ICMOL = IOCCL + NBASL * NOCL
            CALL DZERO(WRK(IEIGL),NORBFTL)
            OPEN(99,FILE='LOCMO',FORM='FORMATTED')
            IJ = 0
            DO IJ = 1, NBASL * NORB1L
              READ(99,'(D30.20)') WRK(ICMOL+IJ-1)
            END DO
            CLOSE(99)
            CALL DCOPY(NBASL*NOCL,WRK(ICMOL),1,WRK(IOCCL),1)
            CALL DCOPY(NORB2L*NBASL,WRK(3*NBASC+1),1,
     *           WRK(IOCCL+NBASL*(NORB1L+NOCL)),1)
            WRITE (LUSIFC) LAB123, LABUNI
            WRITE (LUSIFC) 1,NORBFTL,NBAST,NMOFLTL,NOCL,NORB1L+
     *                     NORB2L+NOCL,NBASL,POTNUC,EMCSCF
            WRITE (LUSIFC) (WRK(IEIGL+I-1),I=1,MAX(NORBFTL,4))
            WRITE (LUSIFC) (WRK(IOCCL+I-1),I=1,MAX(NMOFLTL,4))
         ELSE
            IOFF1 = KFCDIA
            IOFF2 = KFCDIA + NBASC
            IOFFEV = IOFF2 + NBASC
            IOFM1 = 1
            IOFM2 = IOFFEV + NBASC
            IOFMO = IOFM2 + NBASC*NBASC
            NFBAST = 0
            DO ISYM = 1, NSYM
               CALL DCOPY(NOCC(ISYM),WRK(IOFF1),1,WRK(IOFFEV),1)
               IOFFEV = IOFFEV + NOCC(ISYM)
               CALL DCOPY(NORB1(ISYM),WRK(IOFF1),1,WRK(IOFFEV),1)
               IOFFEV = IOFFEV + NORB1(ISYM)
               IOFF1 = IOFF1 + NORB1(ISYM)
               CALL DCOPY(NORB2(ISYM),WRK(IOFF2),1,WRK(IOFFEV),1)
               IOFF2 = IOFF2 + NORB2(ISYM)
               IOFFEV = IOFFEV + NORB2(ISYM)
               NBASI = NBAS(ISYM)
               CALL DCOPY(NOCC(ISYM)*NBASI,CMO(IOFM1),1,WRK(IOFMO),1)
               IOFMO = IOFMO + NOCC(ISYM)*NBASI
               CALL DCOPY(NORB1(ISYM)*NBASI,CMO(IOFM1),1,WRK(IOFMO),1)
               IOFMO = IOFMO + NORB1(ISYM)*NBASI
               IOFM1 = IOFM1 + NORB1(ISYM)*NBASI
               CALL DCOPY(NORB2(ISYM)*NBASI,WRK(IOFM2),1,WRK(IOFMO),1)
               IOFM2 = IOFM2 + NORB2(ISYM)*NBASI
               IOFMO = IOFMO + NORB2(ISYM)*NBASI
               NFBAST = NFBAST +
     *         (NORB1(ISYM) + NORB2(ISYM) + NOCC(ISYM)) * NBASI
            END DO
            WRITE (LUSIFC) LAB123, LABUNI
            WRITE (LUSIFC) NSYM,NORBFT,NBAST,NMOFLT,(NOCC(I),I=1,NSYM),
     *                     (NORB1(I)+NORB2(I)+NOCC(I),I=1,NSYM),
     *                     (NBAS(I),I=1,NSYM),POTNUC,EMCSCF
            WRITE (LUSIFC) (WRK(2*NBASC+I),I=1,MAX(NORBFT,4))
            WRITE (LUSIFC) (WRK((3+NBASC)*NBASC+I),I=1,MAX(NMOFLT,4))
         END IF
      END IF
C
C *** end of subroutine WR_SIRIFC
C
      WRITE (LUSIFC) LAB123,LABEOD
      CALL GPCLOSE(LUSIFC,'KEEP')
      CALL QEXIT('WR_SIRIFC')
      RETURN
      END
