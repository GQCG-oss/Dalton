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
#ifdef REVISION_LOG
===========================================================================
Revision 1.2  2000/05/01 14:05:16  hjj
removed KAVER reference.
===========================================================================
#endif
C  /* Deck sirnr */
      SUBROUTINE SIRNR(CREF,G,CMO,INDXCI,DV,PV,
     *                 FC,FV,FCAC,H2AC,EACTVN,IBNDX,WRK,LFREE)
C
C     13-May-1986 Hans Joergen Aa. Jensen
C
C     Solve Newton-Raphson set of linear equations
C     for local second-order MC step.
C
C MOTECC-90: The purpose of this module, SIRNR, and the algorithms used
C            are described in Chapter 8 Sections B.5 and E.1 of
C            MOTECC-90 "The Direct Iterative NEO Algorithm" and "Split
C            Configuration and Orbital Trial Vectors"
C
C Input:
C  G, gradient (destroyed in SIRNR, see output list)
C
C Output:
C  G, XKAP is returned in G
C
C Common blocks:
C  ITMICT, total number of micro iterations (is updated).
C  MAXJT, maximum number of micro iterations in one macro iteration
C
C Scratch:
C  WRK
C
C Local:
C  MAXSIM, the number of B and Sigma vectors we have space
C          for in core simultaneously.
C
C
C
#include <implicit.h>
      DIMENSION CREF(*),G(*),CMO(*),DV(*),PV(*)
      DIMENSION FC(*),FV(*), FCAC(*),H2AC(*)
      DIMENSION INDXCI(*),IBNDX(*),WRK(*)
C
#include <iratdef.h>
#include <ibndxdef.h>
C
C -- local variables and constants:
C
      PARAMETER (D0=0.0D0, D1=1.0D0)
      PARAMETER (THREPS = 1.D-12)
C
C Used from common blocks:
C  MAXRL, max. dimension of the reduced L matrix
C         (the projected L matrix)
C
C  INFINP: ISTATE,?
C  INFVAR: NCONF,NWOPT,NVAR
C  INFORB: NNASHX,... (to calculate MAXSIM)
C  INFDIM: MAXRL,LACIMX,LBCIMX,LPHPMX,MAXPHP,?
C  INFOPT: RTRUST,RTTOL,ITMAC,NREDL,?
C  ITINFO: DINFO(*),IINFO(*)
C  INFTAP: LUIT3,LUIT5,?
C  INFPRI: P6FLAG,...
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <infopt.h>
#include <itinfo.h>
#include <inftap.h>
#include <infpri.h>
C
      SAVE NCONSV,NWOPSV,MXRLSV,LFRESV, LA,LB,LC,LBC,MAXSIM
      SAVE IGORB,KREDH,KREDG,KREDV,KLDIA,KWRK1,LWRK1
C
      DATA NCONSV,NWOPSV,MXRLSV,LFRESV /-1,-1,-1,-1/
C
      CALL QENTER('SIRNR ')
C
C *** Work space allocation
C
      IF ( NCONF .NE. NCONSV .OR. NWOPT .NE. NWOPSV .OR.
     *     MAXRL .NE. MXRLSV .OR. LFREE .NE. LFRESV) THEN
         NCONSV = NCONF
         NWOPSV = NWOPT
         MXRLSV = MAXRL
         LFRESV = LFREE
         IGORB = 1 + NCONF
         CALL PHPINI(LPHPMX,NCONF,NWOPT,MAXPHP)
C        CALL PHPINI(LPHPT,NCVAR,NOVAR,MAXPHP)
         KREDG = 1
         KREDH = KREDG + MAXRL
         KREDV = KREDH + (MAXRL*(MAXRL+1))/2
         KLDIA = KREDV + MAXRL
         KWRK1 = KLDIA + LPHPMX
C
         LWRK1 = (LFREE+1) - KWRK1
C
C === allocation for LINTRN
C
C        work space needed for LINTRN should be .lt.
C        LA + NCSIM*LB + NOSIM*LC
C
C        CALL LINMEM(NCSIM,NOSIM,LNEED)
         CALL LINMEM(1,0,LC1)
         CALL LINMEM(2,0,LC2)
         CALL LINMEM(0,1,LO1)
         CALL LINMEM(0,2,LO2)
         LB   = LC2 - LC1 + NCONF
         LC   = LO2 - LO1 + NWOPT
         LA   = MAX(LC1-LB,LO1-LC)
         LBC = MAX( LB , LC )
         MAXSIM = (LWRK1 - LA) / (2 + LBC)
         MAXSIM = MIN(MAXSIM,2)
         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,9011) NCONF,NWOPT,MAXRL,
     *      LA,LB,LC,LBC,LFREE,LWRK1,MAXSIM
         IF (MAXSIM.LE.0) THEN
            KWRK2 = KWRK1 + LA + LBC
            CALL ERRWRK('SIRNR.LINTRN',KWRK2,LFREE)
         END IF
      END IF
 9011 FORMAT(/' (SIRNR) LINTRN needs (at most) ',
     *        'LA + NCSIM*LB + NOSIM*LC',
     *       /' NCONF,NWOPT,MAXRL               :',3I10,
     *       /' LA,LB,LC,LBC                    :',4I10,
     *       /' LFREE(total), LFREE(for LINTRN) :',2I10,
     *       /' MAXSIM                          :',I10)
C
C
      NCOMAX = MAX(NCONF,NWOPT)
      NCONF4 = MAX(4,NCONF)
      NWOPT4 = MAX(4,NWOPT)
C
C     space for LINTRN
C
      KBVECS = KWRK1
      KSVECS = KBVECS + MAXSIM*NCOMAX
      KWRK2  = KSVECS + MAXSIM*NVAR
C     space for SNRRED
      KWRK2  = MAX(KWRK2,KBVECS + (MAXRL*(MAXRL+1)/2) + 2*MAXRL)
      IF (KWRK2 .GT. LFREE)
     *   CALL ERRWRK('SIRNR.(LINTRN,SNRRED)',KWRK2,LFREE)
C
      LBVECS = (LFREE+1) - KBVECS
      LSVECS = (LFREE+1) - KSVECS
C
C *******************************************************
C  Initialize various entities we will need,
C     NREDH is dimension of reduced Hessian
C     NREDL is number of b-vectors on LUIT3 (CREF + NREDH trial vectors)
C
      NREDH = 0
      NREDL = 1
      DO 50 I = 1,MAXRL
   50    IBNDX(I) = 0
C
C *** ***
C     Generate and save PHP matrix on LUIT2 (SIRPHP).
C     Construct the trial vectors for the first micro iteration
C
      LBVST  = LFREE - KBVECS
      CALL SIRPHP(WRK(KLDIA),EACTVN,INDXCI,FCAC,H2AC,
     *            WRK(KBVECS),LBVST)
C     CALL SIRPHP(DIAGL,EACTVN,XNDXCI,FCAC,H2AC,WRK,LWRK)
      CALL SNRST (NCSIM,NOSIM,CREF,G,IBNDX,EACTVN,WRK(KLDIA),
     *            WRK(KBVECS),LBVST)
C     CALL SNRST (NCSIM,NOSIM,CREF,G,IBNDX,EACTVN,DIAGL,BVECS,LBVECS)
      NTSIM = NCSIM + NOSIM
C
C *************************************
C *** Start of micro iteration loop ***
C *************************************
C
      TIMMIC = SECOND()
      TIMLIN = D0
      NLIN   = 0
      NCLIN  = 0
      NOLIN  = 0
      ITMIC  = 0
C
  100 ITMIC  = ITMIC + 1
C
         IF (NTSIM .GT. MAXSIM) THEN
            WRITE (LUERR,'(//A,2I5)')
     *         ' SIRNR-FATAL, too many trial vectors; NTSIM,MAXSIM =',
     *         NTSIM,MAXSIM
            CALL QUIT('SIRNR-FATAL: insufficient space for '//
     &                'trial vectors')
         END IF
C
C        test against max. dimension of reduced Hessian
C
         IF (NREDH + NTSIM .GT. MAXRL) GO TO 940
C                                    ------------v
C
C
C        *** Calculate linear transformation induced by K-matrix
C        NOTE: CI BVEC's are assumed orthogonal to CREF, otherwise
C              LINTRN will not return correct result for
C              K(orb,CI) * b-vec (see comments in ORBSIG)
C
         TIMLIN = TIMLIN - SECOND()
         NLIN  = NLIN  + 1
         NCLIN = NCLIN + NCSIM
         NOLIN = NOLIN + NOSIM
         CALL LINTRN(NCSIM,NOSIM,WRK(KBVECS),
     *               WRK(KBVECS+NCSIM*NCONF),CMO,CREF,EACTVN,
     *               G(IGORB),DV,PV,FC,FV,FCAC,H2AC,
     *               INDXCI,WRK(KSVECS),1,LSVECS)
C        CALL LINTRN(NCSIM,NOSIM,BCVECS,BOVECS,
C    *               CMO,CREF,EACTVN,GORB,DV,PV,
C    *               FC,FV,FCAC,H2AC,INDXCI,WRK,KFRSAV,LFRSAV)
C
         TIMLIN = TIMLIN + SECOND()
C
C        To save SVECS for PHP = PLP, that is for beta = 0, we must
C        subtract CREF part of sigma vectors:
C
         ISCVEC = KSVECS
         DO 780 ITSIM = 1,NTSIM
            TT = DDOT(NCONF,CREF,1,WRK(ISCVEC),1)
            IF (ABS(TT) .GT. THREPS) THEN
               IF (P6FLAG(17))
     *            WRITE (LUPRI,7800) ITMIC,(NREDH+ITSIM),TT
               CALL DAXPY(NCONF,(-TT),CREF,1,WRK(ISCVEC),1)
            END IF
            ISCVEC = ISCVEC + NVAR
  780    CONTINUE
 7800    FORMAT(/' (SIRNR) micro iteration',I3,/T11,
     *      'CI reference vector projected out of sigma vector no.',
     *      I3,/T11,'Overlap was:',1P,D20.12)
C
C        *** Save SVECS on LUIT5
C           (SVECS is placed at WRK(KSVECS) by LINTRN)
C
         REWIND LUIT5
         DO 600 I = 1,NREDL
            IF (NCONF.GT.1) READ (LUIT5)
            IF (NWOPT.GT.0) READ (LUIT5)
  600    CONTINUE
         ISC = KSVECS
         ISO = ISC + NCSIM*NVAR
         DO 150 ISIM = 1,NTSIM
            IF (IBNDX(NREDL+ISIM) .EQ. JBCNDX) THEN
               ISVEC = ISC
               ISC   = ISC + NVAR
            ELSE
               ISVEC = ISO
               ISO   = ISO + NVAR
            END IF
            IF (NCONF.GT.1) CALL WRITT (LUIT5,NCONF4,WRK(ISVEC))
            IF (NWOPT.GT.0) THEN
               ISVEC = ISVEC + NCONF
               CALL WRITT (LUIT5,NWOPT4,WRK(ISVEC))
            END IF
  150    CONTINUE
         REWIND LUIT5
C
         NREDL = NREDL + NTSIM
         NREDH = NREDL - 1
         CALL SNRRED (3,IPRI6,NREDH,NTSIM,
     *                WRK(KREDH),WRK(KREDG),WRK(KREDV),G,IBNDX,
     *                WRK(KBVECS),WRK(KBVECS+NCSIM*NCONF),WRK(KSVECS))
C        CALL SNRRED (ICTL,IPRNR,NREDH,N,REDH,REDG,SOLNR,G,
C    *                IBNDX,BCNRVE,BONRVE,SNRVEC)
C
         IF (P6FLAG(10)) THEN
            WRITE (LUPRI,8010) ITMAC,ITMIC
            IF (P6FLAG(34)) THEN
               CALL OUTPAK(WRK(KREDH),NREDH,1,LUPRI)
            END IF
            WRITE (LUPRI,8012) (-WRK(KREDG-1+I),I = 1,NREDH)
            WRITE (LUPRI,8018) (WRK(KREDV-1+I),I = 1,NREDH)
         END IF
 8010    FORMAT(/' (SIRNR) The reduced Hessian iteration no. (',
     *           I3,',',I3,')')
 8012    FORMAT(/' - the reduced gradient',//,(4X,1P,5D15.6))
 8018    FORMAT(/' - and the reduced solution vector',
     *         //,(4X,1P,5D15.6))
C
#if defined (VAR_OLDCODE)
         IF (NREDH .LE. 4) THEN
            JCONV = 0
         ELSE
            JCONV = 1
         END IF
CHJ-860923: SIRNR is only called in local regions where we have
C           a 1-to-1 relation between the residual and the next
C           gradient; hence we want to check convergence as soon
C           as possible.
#endif
         JCONV = 1
         CALL SNRNEX (IPRI4,JCONV,NCSIM,NOSIM,IBNDX,G,WRK(KLDIA),
     *                WRK(KREDV),CMO,CREF,DV,PV,FC,FV,
     *                WRK(KBVECS),LBVECS)
C        CALL SNRNEX (IPRNR,JCONV,NCSIM,NOSIM,IBNDX,G,DIAGL,
C    *                SOLNR,CMO,CREF,DV,PV,FC,FV,WRK,LFREE)
C
         IF (JCONV.EQ.1) GO TO 9000
C                        ------------v
C     *** micro iterations converged if jconv .eq. 1 ***
C
         NTSIM = NCSIM + NOSIM
         IF (JCONV .EQ. -1 .OR. NTSIM .EQ. 0) THEN
C           (linear dependency)
            GO TO 900
C           -----------v
         END IF
C        (not converged
         IF (ITMIC.GE.MAXJT) GO TO 910
C                            -----------v
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         GO TO 100
C     ^-----------
C
C ***********************************
C *** End of micro iteration loop ***
C ***********************************
C
C     Iterations stopped because of linear dependence
C
  900 CONTINUE
C
      WRITE (LUPRI,9020) ITMAC,ITMIC
      WRITE (LUERR,9020) ITMAC,ITMIC
      GO TO 9000
 9020 FORMAT (/' (SIRNR) iteration (',I3,',',I3,')',
     *        /' microiterations stopped because of linear dependence')
C
C     Iterations stopped because maximum number of micro iterations
C     reached
C
  910 CONTINUE
      WRITE (LUPRI,9030) ITMAC,ITMIC
      WRITE (LUERR,9030) ITMAC,ITMIC
      GO TO 9000
 9030 FORMAT (/' (SIRNR) iteration (',I3,',',I3,')',
     *        /' maximum number of microiterations reached')
C
  940 CONTINUE
      WRITE (LUPRI,9400)  ITMAC,ITMIC,MAXRL,NREDH
      WRITE (LUERR,9400) ITMAC,ITMIC,MAXRL,NREDH
      GO TO 9000
 9400 FORMAT (/' (SIRNR) iteration (',I3,',',I3,')',
     *        /' Microiterations stopped because maximum dimension',
     *         ' of reduced L',I4,
     *        /' was about to be passed; dim(red L) =',I4)
C
C *******************************
C *** Save new CI reference vector on LUIT1 and calculate
C     orbital transformation vector "XKAP" in G vector.
C
C
 9000 CONTINUE
      SHFLVL = D0
C
C *** Calculate predicted energy difference
C     between next macro iteration and this one.
C     NOTE: DEPRED is used in SIRSAV
C
      DEPRED = SNRPRD(NREDH,WRK(KREDG),WRK(KREDH),WRK(KREDV))
C     DEPRED = SNRPRD(NREDH,REDG,REDH,SOLEQ)
C
      CALL SIRSAV ('NRSAVE',CMO,IBNDX,WRK(KREDG),WRK(KREDV),
     1             G,INDXCI,WRK(KWRK1),LWRK1)
C     CALL SIRSAV (KEYWRD,CMO,IBNDX,REDL,EVEC,XKAP,INDXCI,WRK,LFREE)
C
      ITMICT = ITMICT + ITMIC
      TIMMIC = SECOND() - TIMMIC
C
C
      IF (IPRI6 .GE. 2) WRITE (LUPRI,9500) NLIN,NCLIN,NOLIN
      IF (IPRI6 .GE. 5) WRITE (LUPRI,9510) NREDH,(IBNDX(I+1),I=1,NREDH)
 9500 FORMAT(/' (SIRNR) LINTRN has been called',I3,' times'
     *       /T11,'with',I3,' CI and',I3,' orbital trial vectors.')
 9510 FORMAT(/T11,'IBNDX index vector for the',I4,' trial vectors:'
     *       /,(T11,5I3,3X,5I3,3X,5I3,3X,5I3))
C
C     save information for final summary print-out
C
      IINFO(2) = ITMIC
      ITMIC    = 0
      IINFO(3) = NREDH
      IINFO(4) = NLIN
      IINFO(5) = NCLIN
      IINFO(6) = NOLIN
C
      DINFO(7)  = D0
C     DINFO(7)  = BETA
      DINFO(8)  = STPLEN
      DINFO(10) = D1
C     DINFO(10) = GAMMA
      DINFO(17) = TIMMIC
      DINFO(19) = TIMLIN
C
      IINFO(13) = ISTATE
      IINFO(14) = ISTATE
      DINFO(21) = D0
      DINFO(31) = D0
      DINFO(41) = WRK(KREDV)
      DINFO(51) = WRK(KREDV+1)
C
C *********************
C *** End of SIRNR ***
C *********************
C
      CALL QEXIT('SIRNR ')
      RETURN
      END
C  /* Deck snrst */
      SUBROUTINE SNRST (NCSIM,NOSIM,CREF,G,IBNDX,
     *                  EACTVN,DIAGL,BVECS,LBVECS)
C
C Written 11-Feb-1985 Hans Joergen Aa. Jensen (was in SIRNEO)
C Revision 9-Jan-1986 hjaaj+pj (+DIAGL, will now also be called
C                               from SIRNR)
C 14-May-1986 PJ
C 14-Jan-1987 hjaaj: solvent contributions to diagonal.
C
C Purpose:
C  NR start. Find the trial vectors for the first micro
C  iteration and define optimization parameters in common /infopt/.
C
#include <implicit.h>
C
      DIMENSION CREF(*),G(*),BVECS(*),IBNDX(*),DIAGL(*)
C
#include <ibndxdef.h>
C
#include <thrldp.h>
      PARAMETER ( D0 = 0.0D0 , DP5 = 0.5D0 , D2 = 2.0D0 )
      PARAMETER ( DTEST = 1.0D-4 )
C
C Used from common blocks:
C  INFINP: FLAG(*),ISTATE,LROOTS,NROOTS,CORHOL
C  INFVAR: NCONF,NWOPT,NVAR,JWOP(2,*)
C  INFOPT: NROOTC,JROOT(*)
C  INFTAP: LUIT2, LUIT3, LUIT5
C
#include <maxorb.h>
#include <priunit.h>
#include <infpri.h>
#include <infinp.h>
#include <infvar.h>
#include <infopt.h>
#include <inftap.h>
C
C
      LOGICAL     GETBVO
      CHARACTER*8 LABIT2(3)
      DATA  LABIT2/'CIDIAG2 ','ORBDIAG ','SOLVDIAG'/
C
      CALL QENTER('SNRST ')
C
C     (space for both CI and orbital trial vectors)
C
      MAXSIM = LBVECS / NVAR
      IF (MAXSIM .LE. 0) CALL ERRWRK('SNRST',NVAR,LBVECS)
C
C
      NCOMAX = MAX(NCONF,NWOPT)
      NCONFX = MAX(4,NCONF)
      NWOPTX = MAX(4,NWOPT)
C
C     Check diagonal for negative elements
C     (so the user is warned)
C
      IF (NCONF .GT. 1) THEN
         NNEG = 0
         DO 100 I = 1,NCONF
            IF (DIAGL(I).LT.D0) NNEG = NNEG + 1
  100    CONTINUE
         IF (NNEG .GT. 0) THEN
            NWARN = NWARN + 1
            WRITE (LUW4,'(/A,I5,A)') ' ***WARNING from SNRST:',NNEG,
     *         ' negative elements of CI diagonal'
         END IF
      END IF
C
      IF (NWOPT .GT. 0) THEN
         NSML = 0
         ISML = 0
         DSML = DTEST
         NNEG = 0
         DO 105 I = NCONF+1,NVAR
            IF (DIAGL(I).LT.DTEST) THEN
               IF (ABS(DIAGL(I)) .LT. DTEST) THEN
                  NSML  = NSML + 1
                  IF (ABS(DIAGL(I)) .LT. ABS(DSML)) THEN
                     ISML = I
                     DSML = DIAGL(I)
                  END IF
C                 WRITE(LUW4,'(/A,I8,A,F16.9,/A,F12.9)')
C    *            ' **WARNING DIAGONAL ORBITAL HESSIAN',I-NCONF,
C    *            ' IS',DIAGL(I),
C    *            ' WHICH IS SMALLER THAN THE TOLERANCE',DTEST
                  DIAGL(I) = SIGN(DTEST,DIAGL(I))
C                 WRITE (LUW4,'(A,F16.9)') ' VALUE USED :',DIAGL(I)
               END IF
               IF (DIAGL(I).LT.D0) NNEG = NNEG + 1
            END IF
  105    CONTINUE
         IF (NSML .GT. 0) THEN
            NWARN = NWARN + 1
            WRITE (LUW4,'(/A,I5,A/A,1P,D12.5)')
     *         ' ***WARNING from SNRST:',NSML,
     *         ' small elements of orbital Hessian diagonal',
     *         ' These elements have been changed to +/-',DTEST
            ISML = ISML - NCONF
            WRITE (LUW4,'(A,1P,D12.3,A,2I5)')
     *         ' Smallest element found was',DSML,
     *         ' for orbital rotation',JWOP(1,ISML),JWOP(2,ISML)
         END IF
         IF (NNEG .GT. 0) THEN
         IF (CORHOL) THEN
            WRITE (LUW4,'(/A,/I5,A)') ' ***core hole info from SNRST:',
     *         NNEG,' negative elements of orbital Hessian diagonal'
         ELSE
            NWARN = NWARN + 1
            WRITE (LUW4,'(/A,I5,A)') ' ***WARNING from SNRST:',NNEG,
     *         ' negative elements of orbital Hessian diagonal'
         END IF
         END IF
      END IF
C
C     First basis vector for L-reduced is CREF.
C     Corresponding sigma vector is the gradient.
C
      REWIND (LUIT3)
      CALL WRITT(LUIT3,NCONFX,CREF)
      REWIND (LUIT5)
      IF (NCONF.GT.1) CALL WRITT(LUIT5,NCONFX,G)
      IF (NWOPT.GT.0) CALL WRITT(LUIT5,NWOPTX,G(1+NCONF))
      IBNDX(1) = JBCNDX
      NBPREV = 1
C
C
      KBC = 1
C
      IF (NCONF.GT.1) THEN
         NCSIM = 1
         IBNDX(NBPREV+NCSIM) = JBCNDX
         IPRPHP = MAX(0,IPRI6-5)
         KWPHP  = KBC + NCONF
         LWPHP  = LBVECS + 1 - KWPHP
         CALL NEXCI(.FALSE.,D0,NCONF,BVECS(KBC),DUMMY,G,DIAGL,
     *              IPRPHP,BVECS(KWPHP),LWPHP)
C        CALL NEXCI(OLSEN,ENER,NCVAR,D,XVEC,
C    &              RES,DIAG,IPRPHP,WRK,LWRK)
         KBO = KBC + NCONF
      ELSE
         NCSIM = 0
         KBO = KBC
      END IF
      GETBVO = (NWOPT.GT.0 .AND. (FLAG(55) .OR. NCONF.EQ.1) )
   25 IF ( GETBVO ) THEN
         NOSIM = 1
         IBNDX(NBPREV+NCSIM+NOSIM) = JBONDX
         DO 30 I = 1,NWOPT
            BVECS((KBO-1)+I) = G(NCONF+I) / DIAGL(NCONF+I)
   30    CONTINUE
      ELSE
         NOSIM = 0
      END IF
C
C     Do microcanonical average on new BOVECS
C     (AVERAG returns immediately if no averageing requested)
C
      IF (NOSIM.GT.0) THEN
         CALL AVERAG(BVECS(KBO),NWOPT,NOSIM)
      END IF
C
C     Orthogonalize and renormalize start vectors.
C     NBPREV is updated in ORTBVC.
C
      KBPREV = KBO + NOSIM*NWOPT
      KTT0   = KBPREV + NCOMAX
      IF (KTT0+NCSIM+NOSIM .GT. LBVECS)
     *   CALL ERRWRK('SNRST.ORTBVC',(KTT0+NCSIM+NOSIM),LBVECS)
      NDMBC  = MAX(1,NCONF)
      NDMBO  = MAX(1,NWOPT)
      THRLDV = NVAR*THRLDP
      CALL ORTBVC (NCSIM,NOSIM,NDMBC,NDMBO,NBPREV,IBNDX,LUIT3,
     *             BVECS(KBC),BVECS(KBO),
     *             THRLDV,BVECS(KTT0),BVECS(KBPREV))
C     CALL ORTBVC(NBC,NBO,NDMBC,NDMBO,NBPREV,IBNDX,LUBVC,
C    *            BCVECS,BOVECS,THRLDP,TT0,BPREV)
C
C
      IF (NCSIM.EQ.0) THEN
         IF (NOSIM.GT.0) THEN
            IF (KBC.NE.KBO) THEN
               DO 120 I = 0,NWOPT-1
                  BVECS(KBC+I) = BVECS(KBO+I)
  120          CONTINUE
            END IF
         ELSE IF (NWOPT.GT.0 .AND. .NOT.GETBVO) THEN
C           if csf vector didn't work out, try an orbital if possible
            GETBVO = .TRUE.
            GO TO 25
C    ^--------------
         ELSE
            WRITE(LUPRI,'(//A)')
     *       ' ** fatal error: NO INITIAL TRIAL VECTORS FORMED IN SNRST'
            CALL QUIT('SNRST-FATAL, no initial trial vectors found.')
         END IF
      END IF
C
      CALL QEXIT('SNRST ')
      RETURN
C
C *** end of subroutine SNRST
C
      END
C  /* Deck snrnex */
      SUBROUTINE SNRNEX (IPRNR,JCONV,NCSIM,NOSIM,IBNDX,G,DIAGL,
     *                   SOLNR,CMO,CREF,DV,PV,FC,FV,WRK,LFREE)
C
C PURPOSE: 1) CONSTRUCT RESIDUAL R(I) = (A-B)*X(I)+F
C             FOR SOLUTION X(I) OF REDUCED N-R EQUATION
C          2) TEST FOR CONVERGENCE OF N-R SOLUTION
C             CONVERGENCE CRITERION:  // R(I) // .LE. THRQN
C          3) USE GENERALIZED CONJUGATE GRADIENT ALGORITHM
C             TO DETERMINE NEXT GUESS OF TRIAL VECTOR
C             IN CSF SPACE AND ORBITAL OPTIMIZATION
C             IN NEXKAP TO DETERMINE TRIAL ORBITAL
C             VECTORS. EITHER THE ORBITAL OR THE CSF PARAMETERS
C             ARE IMPROVED IN ONE ITERATION NEVER BOTH.
C
C MOTECC-90: Some of the algorithms used in this module, SNRNEX,
C            are described in Chapter 8 Section E.3 of MOTECC-90
C            "Convergence of Solution vectors in Direct NEO algorithms"
C
C JCONV  input: if JCONV .lt. 0 do not calculate new trial vectors
C               if JCONV .eq. 0 do not check convergence
C       output: =  1 converged
C               =  0 not converged
C               = -1 not converged, linear dependency among all
C                                   trial vectors.
C
C PJ JAN 1985
C
#include <implicit.h>
      DIMENSION IBNDX(*),G(*),DIAGL(*),SOLNR(*)
      DIMENSION CMO(*),CREF(*),DV(*),PV(*),FC(*),FV(*),WRK(*)
C
C
      PARAMETER ( FKAPC = 0.2D0, FKAPO = 0.01D0)
C
C
#include <thrldp.h>
      PARAMETER ( MAXTST = 30 , QCNRAT = 0.80D0 )
      PARAMETER ( DP1 = 0.1D0, D1=1.0D0 , D0=0.0D0 )
      PARAMETER ( THRSIM = 1.D-6 )
#include <dummy.h>
C
#include <ibndxdef.h>
C
C
C Used from common blocks:
C
C  INFINP: FLAG(*),THRGRD
C  INFVAR: NCONF,NWOPT,NVAR
C  INFOPT: NREDL,GRDNRM
C  INFTAP: LUIT3,LUIT5
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <infopt.h>
#include <inftap.h>
#include <infpri.h>
C
C
C     local variables:
C
      LOGICAL KAPOPT, FORCEO, CTRUNC, SYMTRZ
C
#include <ibtdef.h>
C
C     **********
C     ********** End of declarations, start execution.
C     ********** First, define some entities we will need.
C     **********
C
      CALL QENTER('SNRNEX')
      NBPREV = NREDL
      IF (NWOPT.GT.0) THEN
         KAPOPT = FLAG(41)
         IGORB  = 1 + NCONF
      ELSE
         KAPOPT = .FALSE.
      END IF
      CTRUNC = FLAG(77)
C     ... "truncate trial vectors", i.e. zero small elements.
      SYMTRZ = FLAG(78)
C     ... symmetrize CI vector
      THRLDV = NVAR*THRLDP
      TLINDP = SQRT(THRLDV)
      TLINDP = 10*TLINDP
      IF (THRGRD .LT. TLINDP) THEN
         WRITE(LUPRI,'(//A/A,2(/A,D10.2))')
     *   ' SNRNEX ERROR: Gradient threshold is less than the',
     *   '               numerical linear dependence limit',
     *   ' Convergence (gradient) threshold',THRGRD,
     *   ' Linear dependence limit         ',TLINDP
         CALL QUIT('SNRNEX: convergence threshold below lin.dep. limit')
      END IF
C
C
C Space for SNRNEX:
C
C MAXNEX: MAXIMUM NUMBER OF SIMULTANEOUS VECTORS IN SNRNEX
C
C     MAXNEX = (LFREE-NVAR)/(NVAR+NWOPT)
      MAXNEX = 1
C
C     CALCULATE CONVERGENCE THRESHOLD
C
      IF (JCONV .EQ. 0) THEN
C     IF (not check convergence) only exit if space spanned
         THRQN = TLINDP
      ELSE
         THRQN = CNVTHR(.FALSE.)
C        THRQN = CNVTHR(LINCNV)
      END IF
C
C     SPACE ALLOCATION
C
      KB    = 1
      KS    = KB + MAXNEX*NVAR
      KWKAP = KS + NVAR
C
C     CONSTRUCT RESIDUAL IN BNRVEC
C     RESIDUAL: (A-B)*X(I)+F
C
C     PUT F IN BNRVEC
C
      CALL DCOPY(NVAR,G,1,WRK(KB),1)
      CALL DZERO(WRK(KWKAP),NWOPT)
C
C     CONSTRUCT  (A-B)*X(I) WHERE X(I) IS THE SOLUTION TO THE
C     I'TH SET OF NEWTON-RAPHSON EQUATIONS
C
      REWIND LUIT5
C
C     FIRST VECTOR CORRESPOND TO THE CI VECTOR AT THE EXPANSION
C     POINT.
C     FORCEO is used to ensure we have at least one orbital trial
C     vector.
C
      IF (NCONF.GT.1) READ (LUIT5)
      IF (NWOPT.EQ.0) THEN
         FORCEO = .FALSE.
      ELSE
         READ (LUIT5)
         FORCEO = .TRUE.
      END IF
      DO 900 K=1,NREDL-1
         IF(NCONF.GT.1) THEN
            CALL READT(LUIT5,NCONF,WRK(KS))
         ELSE
            WRK(KS) = D0
         END IF
         IF (NWOPT.GT.0) CALL READT(LUIT5,NWOPT,WRK(KS+NCONF))
         EVAL1 = SOLNR(K)
         IF (IBNDX(K+1).EQ.JBONDX)THEN
            FORCEO = .FALSE.
            CALL DAXPY(NVAR,EVAL1,WRK(KS),1,WRK(KB),1)
         ELSE
            CALL DAXPY(NCONF,EVAL1,WRK(KS),1,WRK(KB),1)
            CALL DAXPY(NWOPT,EVAL1,WRK(KS+NCONF),1,WRK(KWKAP),1)
         END IF
 900  CONTINUE
      CALL DAXPY(NWOPT,D1,WRK(KWKAP),1,WRK(KB+NCONF),1)
      CALL DAXPY(NWOPT,D1,G(1+NCONF),1,WRK(KWKAP),1)
C
C     RESIDUAL IS NOW IN WRK(KB) AND MODIFIED GRADIENT IN WRK(KWKAP)
C
C     TEST FOR CONVERGENCE OF THE N-R SOLUTIONS
C
      QCNORM = DNORM2(NCONF,WRK(KB),1)
      QONORM = DNORM2(NWOPT,WRK(KB+NCONF),1)
      QBTEST = SQRT(QCNORM**2+QONORM**2)
      IF (QBTEST .LE. THRGRD) THEN
C        This should be sufficient for MC optimization to be
C        converged
         THRQNX = THRQN
      ELSE IF (QCNORM .LT. DP1*THRQN) THEN
C        exclude case where next check is not appropriate :
C        if CSF converged (e.g. for Hartree-Fock)
         THRQNX = THRQN
      ELSE IF (QONORM .GT. QCNRAT*QCNORM) THEN
C        no convergence check if orbital residual bigger than QCNRAT of
C        CSF residual, this ensures greater stability in optimization.
         THRQNX = TLINDP
      ELSE
         THRQNX = THRQN
      END IF
      IF (IPRNR .GT. 5) THEN
         WRITE (LUW4,'(/A/I10,1P,4D16.4)')
     *   ' Dimension    CSF residual   orb. residual  total residual'//
     *   '   conv. thresh.',
     *   (NREDL-1),QCNORM,QONORM,QBTEST,THRQNX
      END IF
      NCSIM = 0
      NOSIM = 0
      IF (QBTEST .GT. THRQNX) THEN
         JCONV = 0
         THRKAP = MAX(FKAPC*QCNORM,FKAPO*QONORM)
         IF (FORCEO .OR. QONORM.GT.QCNRAT*QCNORM) THEN
            NOSIM=NOSIM+1
            IBNDX(NREDL+1)=JBONDX
         ELSE
            NCSIM=NCSIM+1
            IBNDX(NREDL+1)=JBCNDX
         END IF
      ELSE
         JCONV = 1
         GO TO 3999
      END IF
C
C     USE CONJUGATE GRAD ALGORITHM TO FORM NEW CSF TRIAL VECTORS
C     OR SOLVE LINEAR EQUATIONS TO GET IMPROVED ORBITAL VECTORS
C
      IF (IBNDX(NREDL+1).EQ.JBONDX) THEN
         IF (KAPOPT) THEN
            CALL DCOPY(NWOPT,WRK(KWKAP),1,WRK(KB),1)
         ELSE
            IOFF=KB+NCONF-1
            DO 1350 I=1,NWOPT
               WRK(KB-1+I)=WRK(IOFF+I) / DIAGL(NCONF+I)
 1350       CONTINUE
         END IF
      ELSE
         IPRPHP = MAX(0,IPRI6-5)
         KWPHP  = KWKAP
         LWPHP  = LFREE + 1 - KWPHP
         CALL NEXCI(.FALSE.,D0,NCONF,WRK(KB),DUMMY,WRK(KB),DIAGL,
     *              IPRPHP,WRK(KWPHP),LWPHP)
C        CALL NEXCI(OLSEN,ENER,NCVAR,D,XVEC,
C    &              RES,DIAG,IPRPHP,WRK,LWRK)
      END IF
C
C
C     ORTHOGONALIZE TRIAL VECTORS AND EXAMINE FOR LINEAR DEPENDENCE
C
      KWRK1 = KB + NVAR
      LWRK1 = LFREE - KWRK1
C
C
      IF (KAPOPT .AND. NOSIM .GT. 0) THEN
         SHIFNR = D0
         DAMP   = D0
         CALL NEXKAP(NOSIM,WRK(KB),NREDL,LUIT3,LUIT5,IBNDX,1,
     *               DUMMY,THRKAP,SHIFNR,DAMP,DIAGL(1+NCONF),CMO,
     *               G(IGORB),DV,PV,FC,FV,WRK(KWRK1),LWRK1)
      END IF
C
C     Do microcanonical average on new BOVECS
C     (AVERAG returns immediately if no averageing requested)
C
      IF (NOSIM .GT. 0) THEN
         CALL AVERAG(WRK(KB),NWOPT,NOSIM)
      END IF
C
      IF (CTRUNC .AND. NCSIM .EQ. 1) THEN
C        ... increase efficiency by zeroing elements in
C            trial vectors abs less than 1.d-3 abs largest element.
C            Do not truncate orbital trial vectors as that
C            will destroy KAPOPT performance.
         FACTRN = 1.0D-3
         IBMAX  = IDAMAX(NCONF,WRK(KB),1)
         BTEST  = FACTRN * ABS(WRK((KB-1)+IBMAX))
         DO 4300 I = KB,(KB-1)+NCONF
            IF (ABS(WRK(I)) .LT. BTEST) WRK(I) = D0
 4300    CONTINUE
      END IF
      IF (SYMTRZ .AND. NCSIM .EQ. 1) THEN
C        ... symmetrize CI vector so that numerical
C            inaccuracy won't break symmetry;
C            WARNING : if symmetry were broken in residual vectors,
C                      and thus in sigma vectors, then this will
C                      lead to constant residual.
C
         CNORM = DNORM2(NCONF,WRK(KB),1)
         TLINDV = SQRT(THRLDV)
         IF (CNORM .GT. TLINDV) THEN
            CNORM = D1 / CNORM
            CALL DSCAL(NCONF,CNORM,WRK(KB),1)
            CALL VSYMTR(NCONF,WRK(KB),THRSIM)
         END IF
      END IF
C
      KWRK2  = KWRK1 + NCSIM + NOSIM
      NDMBC  = MAX(1,NCONF)
      NDMBO  = MAX(1,NWOPT)
      CALL ORTBVC (NCSIM,NOSIM,NDMBC,NDMBO,NBPREV,IBNDX,LUIT3,
     *            WRK(KB),WRK(KB),THRLDV,WRK(KWRK1),WRK(KWRK2))
      NTEST = NCSIM + NOSIM
C
C     LINEAR DEPENDENCE BETWEEN TRIAL VECTORS
C
      IF (NTEST.EQ.0) JCONV = -1
C
 3999 CONTINUE
C
C
C     Output:
C     =======
C
      IF (IPRNR .LE. 5 .AND. IPRNR .GT. 0 .AND. JCONV .NE. 0) THEN
         WRITE (LUW4,'(/A/I10,1P,4D16.4)')
     *   ' Dimension    CSF residual   orb. residual  total residual'//
     *   '   conv. thresh.',
     *   (NREDL-1),QCNORM,QONORM,QBTEST,THRQNX
      END IF
      IF (JCONV .EQ. 1 .AND. IPRNR .GT. 0) THEN
C        NEWTON-RAPHSON LINEAR EQUATIONS HAVE CONVERGED
         WRITE(LUW4,5030)
      ELSE IF (JCONV .EQ. -1) THEN
C        LINEAR DEPENDENCE BETWEEN TRIAL VECTORS
         WRITE(LUW4,5010)
      END IF
 5030 FORMAT(/' *** Newton-Raphson eq.s converged')
 5010 FORMAT(/' *** Microiterations stopped due to linear ',
     *        'dependence between all new trial vectors')
C
C     END OF SNRNEX
C
      CALL QEXIT('SNRNEX')
      RETURN
      END
C  /* Deck snrred */
      SUBROUTINE SNRRED (ICTL,IPRNR,NREDH,N,REDH,REDG,SOLNR,G,
     *                   IBNDX,BCNRVE,BONRVE,SNRVEC)
C
C Written May 86 PJ based on NEORED
C
C INPUT:
C  ICTL, flow control
C         =1 ,extend the reduced PROJECTED HESSIAN and gradients
C         =2 ,find new solutions to reduced problem
C         =3 ,both
C  NREDH, dimension of new reduced PROJECTED HESSIAN matrix
C  N, number of new b-vectors
C  BCNRVE, the new csf b-vectors (is destroyed in NRRED)
C  BONRVE, the new orbital b-vectors (is destroyed in NRRED)
C  SNRVEC, the sigma vectors of the N new b-vectors
C
C  (The reduced PROJECTED HESSIAN matrix is the projection
C   on the basis of b-vectors.)
C
C Output:
C  REDH, the new, extended reduced PROJECTED HESSIAN matrix
C        (dimension: NREDH)
C  REDG, the new, extended reduced minus gradient vector
C  SOLNR, solution to the  set of NEWTON RAPHSON equation
C
C Scratch:
C  BCNRVE are used for B vectors from unit LUIT3 (section 1)
C  SNRVEC is used for reduced PROJECTED HESSIAN matrix (section 2)
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION REDH(*),REDG(*),SOLNR(*),G(*)
      DIMENSION IBNDX(*),BCNRVE(*),BONRVE(*),SNRVEC(NVAR,*)
C
#include <ibndxdef.h>
C
C  INFIND : IROW(*),?
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
C
#include <ibtdef.h>
C
      CALL QENTER('SNRRED')
C
C Section 1:extend reduced PROJECTED HESSIAN matrix with N new b-vectors
C
      IF (IBTAND(ICTL,1).EQ.0) GO TO 1000
      IF (N.LT.1) GO TO 1000
      IF (NREDH .GE. LIROW)
     *   CALL QUIT('SNRRED error: NREDH exceeds LIROW')
C        ... IROW(NREDH+1) used below
      IREDH = NREDH - N
C
C New b-vectors (are in B(o,c)NRVE)
C EXTEND THE REDUCED MINUS GRAD VECTOR, AND
C EXTEND THE REDUCED PROJECTED HESSIAN MATRIX
C
C
      ISTBC = 1
      ISTBO = 1
      DO 30 K = 1,N
         KL = IROW(IREDH+K) + IREDH
         IF(IBNDX(IREDH+K+1).EQ.JBCNDX) THEN
            REDG(IREDH+K) = -DDOT(NCONF,BCNRVE(ISTBC),1,G,1)
            DO 40 L=1,K
               REDH(KL+L) = DDOT(NCONF,BCNRVE(ISTBC),1,SNRVEC(1,L),1)
 40         CONTINUE
            ISTBC=ISTBC+NCONF
         ELSE
            REDG(IREDH+K) = -DDOT(NWOPT,BONRVE(ISTBO),1,G(NCONF+1),1)
            DO 41 L=1,K
               REDH(KL+L) =
     *            DDOT(NWOPT,BONRVE(ISTBO),1,SNRVEC(NCONF+1,L),1)
 41         CONTINUE
            ISTBO=ISTBO+NWOPT
         END IF
 30   CONTINUE
C
C Old b-vectors.
C Rewind LUIT3 and skip CREF, the first b-vector
C (THE PROJECTED HESSIAN MATRIX HAS NO CREF COMPONENT)
C
      REWIND LUIT3
      READ (LUIT3)
      LL = 0
      DO 15 J = 1,IREDH
         IF (IBNDX(J+1).EQ.JBCNDX) THEN
            CALL READT(LUIT3,NCONF,BCNRVE)
            LL=LL+1
            DO 27 K=1,N
               REDH(LL+IROW(IREDH+K)) =
     *            DDOT(NCONF,BCNRVE,1,SNRVEC(1,K),1)
   27        CONTINUE
          ELSE
             CALL READT(LUIT3,NWOPT,BCNRVE)
             LL=LL+1
             DO 26 K=1,N
                REDH(LL+IROW(IREDH+K)) =
     *             DDOT(NWOPT,BCNRVE,1,SNRVEC(NCONF+1,K),1)
   26        CONTINUE
          END IF
   15  CONTINUE
C
C *********************************************************
C Section 2: FIND SOLUTION TO THE NEWTON-RAPHSON EQ.
C
 1000 CONTINUE
      IF (IBTAND(ICTL,2).EQ.0) GO TO 2000
C
C SOLVE  N-R EQ.
C Copy REDH to SNRVEC and REDG to SOLNR
C
      KMA = IROW(NREDH+1)
      CALL DCOPY(KMA,REDH,1,SNRVEC,1)
      CALL DCOPY(NREDH,REDG,1,SOLNR,1)
      LIWRK = 1 + KMA / NVAR
      KIWRK = 1 + KMA - (LIWRK-1)*NVAR
      CALL DSPSOL(NREDH,1,SNRVEC,SOLNR,SNRVEC(KIWRK,LIWRK),INFO)
      IF (IPRNR.GE.12 .OR. INFO.NE.0) THEN
         IF (INFO.NE.0) WRITE(LUPRI,7000) INFO
         WRITE(LUPRI,'(/A)')' *** REDUCED HESSIAN MATRIX ***'
         CALL OUTPAK(REDH,NREDH,1,LUPRI)
         WRITE(LUPRI,'(/A)')' *** REDUCED MINUS GRADIENT VECTOR ***'
         CALL OUTPUT(REDG,1,1,1,NREDH,1,NREDH,1,LUPRI)
         WRITE(LUPRI,'(/A)')' * SOLUTION TO NEWTON RAPHSON EQUATIONS *'
         CALL OUTPUT(SOLNR,1,1,1,NREDH,1,NREDH,1,LUPRI)
      END IF
      IF (INFO.NE.0)THEN
         WRITE(LUERR,7000) INFO
         CALL QTRACE(LUERR)
         CALL QUIT('SNRRED: NO SOLUTION FOUND TO LINEAR EQUATIONS')
      END IF
 7000 FORMAT(/' ***SNRRED*** SOLUTION NOT OBTAINED TO LINEAR EQUATIONS',
     *       /' CHECK IF HESSIAN MATRIX IS SINGULAR, INFO =',I5)
C
C *** End of subroutine SNRRED
C
 2000 CONTINUE
      CALL QEXIT('SNRRED')
      RETURN
      END
C  /* Deck snrprd */
      FUNCTION SNRPRD(NREDH,REDG,REDH,SOLEQ)
C
C     May-1986 hjaaj
C
C     CALCULATE SECOND ORDER CHANGE IN THE TOTAL ENERGY
C
C     REDG  = minus reduced gradient
C     REDH  = reduced Hessian
C     SOLEQ = step vector in reduced basis
C
#include <implicit.h>
      DIMENSION REDG(*),REDH(*),SOLEQ(*)
C
      PARAMETER ( D0 = 0.0D0 , DP5 = 0.5D0 )
C
      SUM = D0
      IJ = 0
      DO 200 I = 1,NREDH
         SUM1 = D0
         DO 100 J = 1,(I-1)
            IJ = IJ + 1
            SUM1 = SUM1 + SOLEQ(J) * REDH(IJ)
 100     CONTINUE
         IJ = IJ + 1
         SUM = SUM
     *       + ( -REDG(I) + DP5*SOLEQ(I)*REDH(IJ) + SUM1 ) * SOLEQ(I)
 200  CONTINUE
C
      SNRPRD = SUM
      RETURN
      END
