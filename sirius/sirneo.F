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
C  /* Deck sirneo */
      SUBROUTINE SIRNEO(CREF,G,CMO,INDXCI,DV,PV,
     *                  FC,FV,FCAC,H2AC,EACTVN,IBNDX,WRK,LFREE)
C
C Written 12-Dec-1983 by Hans Jorgen Aa. Jensen and Hans Agren.
C Revisions:
C  20-May-1984 hjaaj
C  29-Oct-1984 hjaaj (new termination criteria included)
C  12-Nov-1984 hjaaj (modified 841029 change when ISTATE .gt. 1)
C   7-Jan-1985 hjaaj (new LUIT5 SIGMA vector format)
C  11-Jan-1985 hjaaj (IBNDX index array, dimension MAXRL)
C
C Purpose:
C  Find the next step vector in an MCSCF optimization using the
C  NEO algorithm.
C
C MOTECC-90: The purpose of this module, SIRNEO, and the algorithms used
C            are described in Chapter 8 Section C.0 (items 4.0 to 4.5)
C            of MOTECC-90 "Implementation", and in Sections B.5 and E.1
C            "The Direct Iterative NEO Algorithm" and
C            "Split Configuration and Orbital Trial Vectors"
C
C Input:
C  G, total gradient (destroyed in SIRNEO, see output list)
C
C Output:
C  G, XKAP is returned in G
C
C Common blocks:
C  BETA, damping factor (from Fletcher's restricted step algorithm)
C  ITMICT, total number of micro iterations (is updated).
C  MAXJT, maximum number of micro iterations in one macro iteration
C
C Scratch:
C  WRK
C
C Local:
C  MAXSIM, the number of B and Sigma vectors we have space
C          for in core simultaneously.
C  BINMEM, .true. if all NROOTS b-vectors in memory simultaneously
C
C
#include <implicit.h>
      DIMENSION CREF(*),G(*),CMO(*),DV(*),PV(*)
      DIMENSION FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION INDXCI(*),IBNDX(*),WRK(*)
C
#include <iratdef.h>
#include <ibndxdef.h>
C
C -- local variables and constants:
C
      LOGICAL BINMEM, BETUPD, ROTFLP
C
      PARAMETER (DHALF=0.5D0, D0=0.0D0, D1=1.0D0)
      PARAMETER (RKVAD = 0.10D0, THREPS = 1.D-12)
C
C Used from common blocks:
C  MAXRL, max. dimension of the reduced L matrix
C         (the projected L matrix)
C
C   INFINP : NROOTS,?
C   INFVAR : NCONF,NWOPT,NVAR
C   INFORB : NNASHX,... (to calculate MAXSIM)
C   INFDIM : MAXRL,?
C   INFOPT : BETA,BETMIN,RTRUST,RTTOL,RBETA,ITMAC,NREDL,?
C   ITINFO : DINFO(*),IINFO(*)
C   INFTAP : LUIT3,LUIT5,?
C   INFPRI : P6FLAG,...
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
      SAVE IGORB,KREDL,KEVEC,KEVAL,KWRK1,LWRK1,KWRK2
C
      DATA NCONSV,NWOPSV,MXRLSV,LFRESV /-1,-1,-1,-1/
C
      CALL QENTER('SIRNEO')
C
C *** Work space allocation
C
      NCOMAX = MAX(NCONF,NWOPT)
      NCONF4 = MAX(4,NCONF)
      NWOPT4 = MAX(4,NWOPT)
C
C
      IF ( NCONF .NE. NCONSV .OR. NWOPT .NE. NWOPSV .OR.
     *     MAXRL .NE. MXRLSV .OR. LFREE .NE. LFRESV) THEN
         NCONSV = NCONF
         NWOPSV = NWOPT
         MXRLSV = MAXRL
         LFRESV = LFREE
         IGORB  = 1 + NCONF
         KREDL  = 1
         KEVEC  = KREDL + (MAXRL*(MAXRL+1))/2
         KEVAL  = KEVEC + MAXRL*MAXRL
         KWRK1  = KEVAL + MAXRL
C
         LWRK1  = (LFREE+1) - KWRK1
C
C === allocation for LINTRN
C
C work space needed for LINTRN should be .lt. LA + NCSIM*LB + NOSIM*LC
C
C        CALL LINMEM(NCSIM,NOSIM,LNEED)
         CALL LINMEM(1,0,LC1)
         CALL LINMEM(2,0,LC2)
         CALL LINMEM(0,1,LO1)
         CALL LINMEM(0,2,LO2)
         LB   = LC2 - LC1 + NCONF
         LC   = LO2 - LO1 + NWOPT
         LA   = MAX(LC1-LB,LO1-LC)
         LBC  = MAX(LB,LC)
         MAXSIM = (LWRK1 - LA) / (2 + LBC)
         I      = 2*ABS(LROOTS)
         MAXSIM = MIN(MAXSIM,I,NVAR)
         IF (IPRSTAT .GT. 1) WRITE (LUSTAT,9011) NCONF,NWOPT,MAXRL,
     *      LA,LB,LC,LBC,LFREE,LWRK1,MAXSIM
         IF (FLAG(62)) THEN
            IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A)')
     *         ' INFO: MAXSIM reduced to 1 by request (flag(62))'
            MAXSIM = MIN(MAXSIM,1)
         END IF
         IF (MAXSIM.LE.0) THEN
            KWRK2 = KWRK1 + LA + LBC
            CALL ERRWRK('SIRNEO.LINTRN',-KWRK2,LFREE)
         END IF
      END IF
 9011 FORMAT(/' (SIRNEO) LINTRN needs (at most)',
     *        ' LA + NCSIM*LB + NOSIM*LC',
     *       /' NCONF,NWOPT,MAXRL               :',3I10,
     *       /' LA,LB,LC,LBC                    :',4I10,
     *       /' LFREE(total), LFREE(for LINTRN) :',2I10,
     *       /' MAXSIM                          :',I10)
C
C     space for LINTRN
C
      KBVECS = KWRK1
      KSVECS = KBVECS + MAXSIM*NCOMAX
      KWRK2  = KSVECS + MAXSIM*NVAR
C     space for NEORED
      KWRK2  = MAX(KWRK2,KBVECS + (MAXRL*(MAXRL+1)/2) + 2*MAXRL)
      IF (KWRK2 .GT. LFREE)
     *   CALL ERRWRK('SIRNEO.(LINTRN,NEORED)',KWRK2,LFREE)
C
      LBVECS = (LFREE+1) - KBVECS
      LSVECS = (LFREE+1) - KSVECS
#if !defined (VAR_RPATES)
C
C     If flag(61) calculate Hessian matrix explicitly for test
C     ("LINTST")
C
      IF (FLAG(61)) THEN
         CALL LINTST(CMO,CREF,EACTVN,G(IGORB),DV,PV,FC,FV,FCAC,H2AC,
     *               INDXCI,WRK,LFREE)
      END IF
#endif
C
C *******************************************************
C  Fill in first element of REDL (B(1)=CREF, S(1)=betaG)
C  Initialize various entities we will need
C  Initialize step scaling to no scaling
C  Set min no. of trial vectors before convergence check
C  (one extra if not optimal orbital trial vectors, i.e. not flag(41))
C
      WRK(KREDL) = D0
      NREDL      = 1
      BETUPD     = .FALSE.
      DO 50 I = 1,MAXRL
   50    IBNDX(I) = 0
      GAMMA = D1
      IF (NCONF .EQ. 1) THEN
         JCNTST = NROOTS + 1
      ELSE
         JCNTST = NROOTS + 2
C        ... we also want at least one csf trial vector
C            when nconf .gt. 1
      END IF
      IF (.NOT. FLAG(41)) JCNTST = JCNTST + 1
C     ... if (not optimal orbital trial vectors) jcntst = jcntst + 1
C
C *** ***
C     Generate and save PHP matrix on LUIT2 (SIRPHP).
C     Construct the trial vectors for the first micro iteration
C
      CALL PHPINI(LPHPMX,NCONF,NWOPT,MAXPHP)
C     CALL PHPINI(LPHPT,NCVAR,NOVAR,MAXPHP)
      IF (MAXPHP .GT. 0) THEN
         BETMIN = D1
C PHPMAERKE : such that temp. DOPHP true in NEONEX in local region.
      END IF
      KLDIAG = (LFREE+1) - LPHPMX
      LBVST  = KLDIAG - KBVECS
      CALL SIRPHP(WRK(KLDIAG),EACTVN,INDXCI,FCAC,H2AC,
     *            WRK(KBVECS),LBVST)
C     CALL SIRPHP(DIAGL,EACTVN,XNDXCI,FCAC,H2AC,WRK,LWRK)
C     Note that DIAGL may be overwritten in NEOST
      CALL NEOST (NCSIM,NOSIM,CREF,G,WRK(KBVECS),LBVST,
     *            IBNDX,WRK(KLDIAG))
C     CALL NEOST (NCSIM,NOSIM,CREF,G,BVECS,LBVECS,IBNDX,DIAGL)
      NTSIML = NCSIM + NOSIM
      IF (NTSIML.LT.0) THEN
         BINMEM = .FALSE.
         NTSIML = -NTSIML
      ELSE
         BINMEM = (NTSIML.LE.MAXSIM)
      END IF
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
         IF (IPRI6.ge.10) WRITE(LUPRI,'(/A,I5)')
     &      ' ---> Micro iteration no.',ITMIC
C
         ITMICB = 0
  200    CONTINUE
            ITMICB = ITMICB + 1
            NTSIM = MIN(MAXSIM,NTSIML)
            IF (IPRI6.ge.10) WRITE(LUPRI,'(A,2I5)')
     &      ' ---> Batch no., NTSIM',ITMICB,NTSIM
C
C           test against max. dimension of reduced L
C
            IF (NREDL+NTSIM .GT. MAXRL) GO TO 940
C                                       ------------v
C
C
            IF (.NOT.BINMEM) THEN
               REWIND LUIT3
               DO 210 I = 1,NREDL
  210             READ (LUIT3)
               NCSIM = 0
               NOSIM = 0
               DO 320 ISIM = 1,NTSIM
                  IF (IBNDX(NREDL+ISIM) .EQ. JBCNDX) THEN
                     NCSIM = NCSIM + 1
                  ELSE
                     NOSIM = NOSIM + 1
                  END IF
  320          CONTINUE
               ISTBC = KBVECS
               ISTBO = ISTBC + NCSIM*NCONF
               DO 340 ISIM = 1,NTSIM
                  IF (IBNDX(NREDL+ISIM) .EQ. JBCNDX) THEN
                     CALL READT(LUIT3,NCONF,WRK(ISTBC))
                     ISTBC = ISTBC + NCONF
                  ELSE
                     CALL READT(LUIT3,NWOPT,WRK(ISTBO))
                     ISTBO = ISTBO + NWOPT
                  END IF
  340          CONTINUE
            END IF
#if defined (VAR_RPATES)
C
C           *** rpa test for poul *** 860820-hj
C
            IF ( FLAG(61) ) THEN
               WRITE (LUPRI,'(////A/)') ' <<<<< RPATEST START >>>>>'
               CALL RPATES(NCSIM,NOSIM,DV,FC,FV,WRK(KBVECS),
     *                     WRK(KBVECS+NCSIM*NCONF),WRK(KSVECS),LSVECS)
C              CALL RPATES(NCSIM,NOSIM,DV,FC,FV,ZOVEC,ZCVEC,WRK,LFREE)
               WRITE (LUPRI,'(//A/)') ' <<<<< RPATEST END >>>>>'
C              set print flags for wanted test output
               P6FLAG(11) = .TRUE.
               P6FLAG(19) = .TRUE.
C DTV          P6FLAG(20) = .TRUE.
               P6FLAG(21) = .TRUE.
C FXCAC        P6FLAG(22) = .TRUE.
C PTV          P6FLAG(27) = .TRUE.
C H2XAC        P6FLAG(28) = .TRUE.
C FXQ, FTQ     P6FLAG(33) = .TRUE.
            END IF
#endif
C
C           *** Calculate linear transformation induced by K-matrix
C           NOTE: CI BVEC's are assumed orthogonal to CREF, otherwise
C                 LINTRN will not return correct result for
C                 K(orb,CI) * b-vec (see comments in ORBSIG)
C
            IF (IPRI6 .gt. 100) THEN
            IF (NCSIM .gt. 0) THEN
               WRITE (LUPRI,'(/A,I3,A/A,20I3)')
     &           'Printing',NCSIM,' configuration trial vectors:',
     &           'IBNDX: ',(IBNDX(NREDL+ISIM),ISIM=1,NTSIM)
               CALL OUTPUT(WRK(KBVECS),1,NCONF,1,NCSIM,
     &              NCONF,NCSIM,-1,LUPRI)
            END IF
            IF (NOSIM .gt. 0) THEN
               WRITE (LUPRI,'(/A,I3,A/A,20I3)')
     &           'Printing',NOSIM,' orbital trial vectors:',
     &           'IBNDX: ',(IBNDX(NREDL+ISIM),ISIM=1,NTSIM)
               CALL OUTPUT(WRK(KBVECS+NCSIM*NCONF),1,NWOPT,1,NOSIM,
     &              NWOPT,NOSIM,-1,LUPRI)
            END IF
            END IF
            TIMLIN = TIMLIN - SECOND()
            NLIN   = NLIN   + 1
            NCLIN  = NCLIN  + NCSIM
            NOLIN  = NOLIN  + NOSIM
            CALL LINTRN(NCSIM,NOSIM,WRK(KBVECS),
     *                  WRK(KBVECS+NCSIM*NCONF),CMO,CREF,EACTVN,
     *                  G(IGORB),DV,PV,FC,FV,FCAC,H2AC,
     *                  INDXCI,WRK(KSVECS),1,LSVECS)
C           CALL LINTRN(NCSIM,NOSIM,BCVECS,BOVECS,
C    *                  CMO,CREF,EACTVN,GORB,DV,PV,
C    *                  FC,FV,FCAC,H2AC,INDXCI,WRK,KFRSAV,LFRSAV)
C
            TIMLIN = TIMLIN + SECOND()
C
C           Beta dependent correction is done in NEORED, the SVECS
C           saved on LUIT5 correspond to BETA .eq. 0 (apart from
C           s1, the gradient, which corresponds to BETA .eq. 1).
C           With this convention, LUIT5 can also be used by NEXKAP
C           without problems.
C
C           To save SVECS for PLP, that is for beta = 0, we must
C           subtract CREF part of sigmavectors:
C
            ISCVEC = KSVECS
            DO 780 ITSIM = 1,NTSIM
               TT = DDOT(NCONF,CREF,1,WRK(ISCVEC),1)
               IF (ABS(TT) .GT. THREPS) THEN
                  IF (P6FLAG(17))
     *               WRITE (LUPRI,7800) ITMIC,(NREDL+ITSIM),TT
                  CALL DAXPY(NCONF,(-TT),CREF,1,WRK(ISCVEC),1)
               END IF
               ISCVEC = ISCVEC + NVAR
  780       CONTINUE
 7800       FORMAT(/' (SIRNEO) micro iteration',I3,/T11,
     *      'CI reference vector projected out of sigma vector no.',
     *      I3,/T11,'Overlap was:',1P,D20.12)
C
C           *** Save SVECS on LUIT5
C              (SVECS is placed at WRK(KSVECS) by LINTRN)
C
            IF (IPRI6 .gt. 100) THEN
               WRITE (LUPRI,'(/A,I3,A/A,20I3)')
     &           'Printing',NTSIM,' sigma vectors:',
     &           'IBNDX: ',(IBNDX(NREDL+ISIM),ISIM=1,NTSIM)
               CALL OUTPUT(WRK(KSVECS),1,NVAR,1,NTSIM,
     &              NVAR,NTSIM,-1,LUPRI)
            END IF
            REWIND LUIT5
            DO 600 I = 1,NREDL
               IF (NCONF.GT.1) READ (LUIT5)
               IF (NWOPT.GT.0) READ (LUIT5)
  600       CONTINUE
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
  150       CONTINUE
            REWIND LUIT5
C
            NREDL  = NREDL  + NTSIM
            NTSIML = NTSIML - NTSIM
            IF (NTSIML.LE.0) GO TO 800
C        v----------------------------
            CALL NEORED (1,NCSIM,NOSIM,WRK(KBVECS),WRK(KSVECS),BETA,
     *                   G,IBNDX,NREDL,WRK(KREDL),
     *                   WRK(KEVAL),WRK(KEVEC),LBVECS)
C           CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,
C    *                   G,IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
            GO TO 200
C        ^-----------
C
C ===
C
  800    CONTINUE
         CALL NEORED (3,NCSIM,NOSIM,WRK(KBVECS),WRK(KSVECS),BETA,
     *                G,IBNDX,NREDL,WRK(KREDL),WRK(KEVAL),WRK(KEVEC),
     *                LBVECS)
C        CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,
C    *                G,IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
C
         IF (P6FLAG(51)) THEN
C           check symmetry of REDL
            CALL NEORED_CHK(IBNDX,WRK(KWRK1),LWRK1)
         END IF
C
C        *** update BETA in UPDBET as if eigenvector was exact
C           (dynamic update of BETA in micro iterations,
C            dimension of reduced L must be .ge. 3 before this
C            this is reasonable (so at least CREF, GC', and
C            an orbital trial vector from NEONEX))
C
C
         IF (NREDL .GE. JCNTST) THEN
            CALL UPDBET(WRK(KREDL),WRK(KEVEC),WRK(KEVAL),
     *                  WRK(KWRK1),LWRK1)
            JCONV = 1
            BETUPD = .TRUE.
         ELSE
            JCONV = 0
         END IF
         IF (P6FLAG(10)) THEN
            WRITE (LUPRI,8010) ITMAC,ITMIC
            IF (P6FLAG(34)) CALL OUTPAK(WRK(KREDL),NREDL,1,LUPRI)
            WRITE (LUPRI,8012)
            WRITE (LUPRI,8015) BETA,(WRK((KEVAL-1)+I),I=1,NREDL)
            LIMP = MIN(((ISTATE+5)/4) * 4, NREDL)
            WRITE (LUPRI,8018)
            CALL OUTPUT(WRK(KEVEC),1,NREDL,1,LIMP,NREDL,NREDL,1,LUPRI)
         END IF
 8010    FORMAT(/' (SIRNEO) The reduced L matrix iteration no. (',
     *           I3,',',I3,')')
 8012    FORMAT(/' - Beta and eigenvalues')
 8015    FORMAT(/1X,F9.2,1P,4D15.6,/,(10X,1P,4D15.6))
 8018    FORMAT(/' - and the first few eigenvectors')
C
C        If RBETA has been evaluated :
C        1) check for extra termination criterion :
C           EVAL(ISTATE+1) .lt. 0 and not flag(36) : GO TO 920
C
C        2) else if RBETA is outside quadratic region (RKVAD)
C           test only to linear convergence in NEONEX (JCONV = -1)
C
         IF (BETUPD) THEN
            IF (.NOT.FLAG(36) .AND. WRK(KEVAL+ISTATE).LT.D0) GO TO 920
            IF (RBETA .GE. MIN(RKVAD,(RTRUST/RTTOL))) JCONV = -1
         END IF
C
C        hjaaj aug 2006: changed ESHIFT from 0.0 to 0.01 to make sure
C                        that precondition in NEONEX should not
C                        give a factor of more than 100
         ESHIFT = 0.01D0
C
         CALL NEONEX(G,WRK(KEVAL),WRK(KEVEC),EACTVN,ESHIFT,
     *               CREF,CMO,DV,PV,FC,FV,
     *               IBNDX,NCSIM,NOSIM,JCONV,
     &               WRK(KBVECS),LBVECS)
C        CALL NEONEX(G,EVAL,EVEC,EACTVN,ESHIFT,
C    *               CREF,CMO,DV,PV,FC,FV,
C    *               IBNDX,NCSIM,NOSIM,JCONV,WRK,LFREE)
C
         IF (JCONV.EQ.0) GO TO 9000
C                        ------------v
C     *** micro iterations converged if jconv .eq. 0 ***
C
         NTSIML = NCSIM + NOSIM
         IF (NTSIML .EQ. 0) THEN
C           (linear dependency, cannot continue)
            GO TO 900
C           -----------v
         ELSE IF (NTSIML.GT.0) THEN
            BINMEM = (NTSIML.LE.MAXSIM)
         ELSE
            BINMEM = .FALSE.
            NTSIML = -NTSIML
         END IF
         IF (ITMIC.GE.MAXJT) GO TO 910
C                            -----------v
C        (not converged, max no. of microiterations reached)
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         GO TO 100
C     ^-----------
C
C ***********************************
C *** End of micro iteration loop ***
C ***********************************
C
C Iterations stopped because of linear dependence
C
  900 CONTINUE
C
      WRITE (LUPRI,9020) ITMAC,ITMIC
      WRITE (LUERR,9020) ITMAC,ITMIC
      GO TO 9000
 9020 FORMAT (/' (SIRNEO) iteration (',I3,',',I3,')',
     *        /' microiterations stopped because of linear dependence')
C
C     Iterations stopped because maximum number of micro iterations
C     reached
C
  910 CONTINUE
      WRITE (LUPRI,9030) ITMAC,ITMIC
      WRITE (LUERR,9030) ITMAC,ITMIC
      GO TO 9000
 9030 FORMAT (/' (SIRNEO) iteration (',I3,',',I3,')',
     *        /' maximum number of microiterations reached')
C
C     Stopped because of termination criterion  EVAL(ISTATE+1) .lt. 0
C
  920 CONTINUE
      WRITE (LUPRI,9200) ITMAC,ITMIC,WRK(KEVAL+ISTATE)
      GO TO 9000
 9200 FORMAT (/' (SIRNEO) iteration (',I3,',',I3,')',
     *        /' Micro iteration termination criterion met,',
     *         ' EVAL(ISTATE+1) =',1P,D10.2)
C
  940 CONTINUE
      BETUPD = .FALSE.
      WRITE (LUPRI,9400) ITMAC,ITMIC,MAXRL,NREDL
      WRITE (LUERR,9400) ITMAC,ITMIC,MAXRL,NREDL
      GO TO 9000
 9400 FORMAT (/' (SIRNEO) iteration (',I3,',',I3,')',
     *        /' Microiterations stopped because maximum dimension',
     *         ' of reduced L',I4,
     *        /' was about to be passed; dim(red L) =',I4)
C
C *******************************
C *** Save orbitals, CI part of eigenvectors, and REDL on LUIT1.
C     If UPDBET has not been called (i.e. BETUPD false), call it
C     to get correct BETA and GAMMA.
C
C
 9000 CONTINUE
      IF (.NOT. BETUPD)
     &   CALL UPDBET(WRK(KREDL),WRK(KEVEC),WRK(KEVAL),WRK(KWRK1),LWRK1)
      SHFLVL = WRK(KEVAL-1+ISTATE)
C
C *** Calculate predicted energy difference
C     between next macro iteration and this one.
C     NOTE: DEPRED is used in SIRSAV
C
C     NWOPT=0, CI, gets special treatment (to get exact
C     energy change.
C
      IF (NWOPT .GT. 0) THEN
         DEPRED = ( WRK(KEVAL-1+ISTATE) / (BETA*BETA) )
     *          * (GAMMA*(D1-GAMMA) + DHALF*
     *             GAMMA*GAMMA*(WRK(KEVEC+(ISTATE-1)*NREDL)**(-2)))
      ELSE
         DEPRED = DHALF*WRK(KEVAL-1+ISTATE)
      END IF
C
      CALL SIRSAV ('NEOSAVE',CMO,IBNDX,WRK(KREDL),WRK(KEVEC),
     1             G,INDXCI,WRK(KWRK1),LWRK1)
C     CALL SIRSAV (KEYWRD,CMO,IBNDX,REDL,EVEC,XKAP,INDXCI,WRK,LFREE)
C
      ITMICT = ITMICT   + ITMIC
      TIMMIC = SECOND() - TIMMIC
C
C
      IF (IPRI6 .GE. 2) WRITE (LUPRI,9500) NLIN,NCLIN,NOLIN
      IF (IPRI6 .GT. 5 .AND. NCLIN .GT. 0)
     *   WRITE (LUPRI,9510) NREDL-1,(IBNDX(I),I=2,NREDL)
 9500 FORMAT(/' (SIRNEO) LINTRN has been called',I3,' times'
     *       /T11,'with',I3,' CI and',I3,' orbital trial vectors.')
 9510 FORMAT(/T11,'IBNDX index vector for the',I4,' trial vectors:'
     *       /,(T11,5I3,3X,5I3,3X,5I3,3X,5I3))
C
C     If absorption for excited state, test if root flipping
C     has occurred such that absorption must be disabled
C     (860521-hjaaj)
C
      IF (ISTATE.GT.1 .AND. (FLAG(51) .OR. FLAG(52) .OR. FLAG(53))) THEN
         IF (ROTFLP(ISTATE,IBNDX,NREDL,WRK(KREDL),WRK(KWRK1))) THEN
            WRITE (LUPRI,'(///2A//)')
     *         ' *** Orbital absorption switched off because root',
     *         ' flipping has occurred.'
            FLAG(51) = .FALSE.
            FLAG(52) = .FALSE.
            FLAG(53) = .FALSE.
         END IF
      END IF
C
C
C     Save information for final summary print-out
C
      IINFO(2)  = ITMIC
      ITMIC     = 0
      IINFO(3)  = NREDL
      IINFO(4)  = NLIN
      IINFO(5)  = NCLIN
      IINFO(6)  = NOLIN
C
      DINFO(7)  = BETA
      DINFO(8)  = STPLEN
      DINFO(10) = GAMMA
      DINFO(17) = TIMMIC
      DINFO(19) = TIMLIN
C
      IOF  = MAX(0,MIN(ISTATE-5,NROOTS-10))
      IEND = MIN(NROOTS,IOF+10) - IOF
      IINFO(13) = IOF + 1
      IINFO(14) = IEND
      IOEVAL = KEVAL-1+IOF
      I1EVEC = KEVEC+IOF*NREDL
      DO 950 I = 1,IEND
         DINFO(20+I) = WRK(IOEVAL+I)
         DINFO(30+I) = WRK(I1EVEC)
         DINFO(40+I) = WRK(I1EVEC+1)
         DINFO(50+I) = WRK(I1EVEC+2)
         I1EVEC = I1EVEC + NREDL
  950 CONTINUE
C
C *********************
C *** End of SIRNEO ***
C *********************
C
      CALL QEXIT('SIRNEO')
      RETURN
      END
C  /* Deck neonex */
      SUBROUTINE NEONEX(G,EVAL,EVEC,EACTVN,ESHIFT,CREF,CMO,DV,PV,
     *                  FC,FV,IBNDX,NCSIM,NOSIM,JCONV,WRK,LFREE)
C
C Written 1-Mar-1985 by Hans Jorgen Aa. Jensen VERSION 2.0
C This newwritten version replaces earlier versions.
C Revisions:
C   4-May-1987 hjaaj SOLVNT contributions to diagonal
C
C Purpose:
C  Find the next trial vectors for the NEO algorithm.
C  (Second line in parameter list is only used when NEXKAP called.)
C
C MOTECC-90: Some of the algorithms used in this module, NEONEX,
C            are described in Chapter 8 Section E.3 of MOTECC-90
C            "Convergence of Solution vectors in Direct NEO algorithms"
C
C Input:
C  G(NVAR), MC gradient.
C  EVAL(NREDL), eigenvalues of reduced L matrix.
C  EVEC(NREDL,NREDL), corresponding eigenvectors.
C  EACTVN, active energy of CREF.
C  ESHIFT, shift of "Davidson" denominator
C  CREF,CMO,DV,PV,FC,FV matrices passed on to NEXKAP
C
C Input and output:
C  IBNDX(*) input:  type of NREDL previous trial vectors,
C           output: extended with type of NCSIM + NOSIM new do.
C  JCONV input:  .gt. 0; test for quadratic convergence
C                .lt. 0; test for linear convergence
C                .eq. 0; do not test for convergence
C        output: .eq. 0; ISTATE converged
C                .eq. 1; ISTATE not converged
C
C Output:
C  NCSIM abs value is number of CI trial vectors written to LUIT3
C        if positive they are returned in WRK(1)
C  NOSIM abs value is number of orb trial vectors written to LUIT3
C        if positive they are returned in WRK(1+NCSIM*NCONF)
C
C  note: if NCSIM+NOSIM .eq. 0 and JCONF .eq. 0, then not converged
C  ----  and no new trial vectors could be found with the algorithm
C        used (e.g. all removed because of linear dependency).
C
C Scratch:
C  WRK(LFREE)
C
#include <implicit.h>
      DIMENSION CREF(*), CMO(*), DV(*), PV(*), FC(*), FV(*)
      DIMENSION G(*), EVAL(*), EVEC(NREDL,*), IBNDX(*), WRK(*)
C
C     ********** Define parameters, common blocks, and local variables.
C
C     (extensive testing by PJ for reasonable convergence in
C      NEXKAP (850301) has been rationalized into FKAP* factors)
C
#include <thrldp.h>
      PARAMETER (QCNRAT = 0.80D0, QCNFAC = 1.70D0, ABSCNV = 0.70D0)
      PARAMETER (FKAPC  = 0.20D0, FKAPO  = 0.01D0)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DHALF = 0.5D0, DP1 = 0.1D0)
      PARAMETER ( THRSIM = 1.D-6 , DTEST  = 0.01D0 )
#include <ibndxdef.h>
C
C Common blocks:
C   INFOPT: NROOTC, number of residual vectors to construct
C           (number of eigenvalues/-vectors of REDL to use)
C Used from common blocks:
C   INFINP : FLAG(*),THRGRD,JOLSEN
C   INFVAR : NCONF,NWOPT,NVAR
C   INFOPT : BETA,RTRUST,NREDL,NROOTC,JROOT(*),GRDNRM
C   INFDIM : LPHPMX
C   INFTAP : LUIT2,LUIT3,LUIT5
C   INFPRI : LUERR
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <pcmlog.h>
#include <infvar.h>
#include <infopt.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C-----------------------
C CBN+JK 03.01.06
C----------------------
#include <gnrinf.h>
C-----------------------
C CBN+JK 03.01.06
C----------------------
C
C     local variables:
C
      LOGICAL KAPOPT, KAPNEO, CITYPE, TSTCNV, LINCNV
      LOGICAL FORCEO, SOLVNT, CTRUNC, SYMTRZ, OLSEN1, DOPHP
      CHARACTER*8 LABIT2(4)
C
      SAVE QCNOLD
C
#include <ibtdef.h>
C
C     data statements:
C
      DATA QCNOLD/0.0D0/
      DATA LABIT2/'CIDIAG2 ','BETADIAG','ORBDIAG ','SOLVDIAG'/
C
C     **********
C     ********** End of declarations, start execution
C     ********** First, define some entities we will need.
C     **********
C
      CALL QENTER('NEONEX')
      SOLVNT = FLAG(16)
      CTRUNC = FLAG(77)
C     ... "truncate trial vectors", i.e. zero small elements.
      SYMTRZ = FLAG(78)
C     ... symmetrize CI vector
#if defined (VAR_NOPHP)
      DOPHP  = .FALSE.
#else
      DOPHP  = (MAXPHP .GT. 0)
C----------------------
C CBN+JK 03.01.06
C----------------------
      IF (SOLVNT .OR. PCM .OR. QM3) DOPHP = .FALSE.
C-----------------------
C CBN+JK 03.01.06
C----------------------
      IF (BETA .NE. D1) DOPHP = .FALSE.
C     ... only use PHP for beta = 1
C     PHPMAERKE: implement for beta .ne. 1
#endif
      IF (NWOPT.GT.0) THEN
         KAPOPT = FLAG(41)
#if defined (VAR_KAPNEO)
CHJ-860921 KAPNEO disabled, is not working satisfactorily
         KAPNEO = .NOT.FLAG(42)
#endif
      ELSE
         KAPOPT = .FALSE.
      END IF
      THRLDV = NVAR*THRLDP
      TLINDP = SQRT(THRLDV)
      TLINDP = 10*TLINDP
      IF (THRGRD .LT. TLINDP) THEN
         WRITE(LUPRI,'(//A/A,2(/A,D10.2))')
     *   ' NEONEX ERROR: Gradient threshold is less than the',
     *   '               numerical linear dependence limit',
     *   ' Convergence (gradient) threshold',THRGRD,
     *   ' Linear dependence limit         ',TLINDP
         CALL QUIT('NEONEX: convergence threshold below lin.dep. limit')
      END IF
      OLSEN1 = .FALSE.
C     OLSEN1 = (JOLSEN .AND. DOPHP)
C PHPMAERKE : JOLSEN not implemented for NEO yet
C
      IF (JCONV.EQ.0) THEN
         TSTCNV = .FALSE.
      ELSE
         TSTCNV = .TRUE.
      END IF
      IF (JCONV.LT.0) THEN
         LINCNV = .TRUE.
      ELSE
         LINCNV = .FALSE.
      END IF
C
      JCONV  = 1
      NCSIM  = 0
      NOSIM  = 0
      NBVIT3 = NREDL
      NBPREV = NBVIT3
      BETM1  = BETA - D1
      IF (NCONF.GT.1) THEN
         NCONFX = NCONF
         IGX    = 1
      ELSE
         NCONFX = 0
         IGX    = 1 + NCONF
CKR         YCNRM  = D0
CKR         QCNRM  = D0
      END IF
CKR      IF (NWOPT.GT.0) THEN
CKR      ELSE
CKR         YONRM  = D0
CKR         QONRM  = D0
CKR      END IF
C
C     Because all of Y?NRM and Q?NRM multiplied with FAC, we better initialize
C     them to avoid FPE on e.g. Cray-T3E, K.Ruud, May 1-96
C
      YNRM  = D0
      YONRM = D0
      YCNRM = D0
      QNRM  = D0
      QONRM = D0
      QCNRM = D0
C
      NVARX = NCONFX + NWOPT
C
C     **********
C     ********** Find NSIMC: how many we can treat simultaneously.
C     **********
C     MWRK0 is space needed for arrays independent of NSIMC.
C     MWRK1 is space needed per residual vector.
C
      MWRK0 = NVAR  + LPHPMX
      MWRK1 = NVARX + NWOPT
      IF (OLSEN1) MWRK1 = MWRK1  + NCONF
      NSIMC = MIN(NROOTC,((LFREE-MWRK0)/MWRK1))
      IF (FLAG(62)) THEN
         NSIMC = MIN(NSIMC, 1)
      END IF
      IF (IPRSTAT .GT. 7) THEN
         WRITE (LUSTAT,'(/,(A,I10))')
     *      ' (NEONEX) number of non-conv. vectors',NROOTC,
     *      '          number of residuals/batch  ',NSIMC
         WRITE (LUSTAT,'(A,10I5/,(T29,10I5))')
     *      ' Index of non-conv. vectors:',
     *      (JROOT(I),I=1,NROOTC)
      END IF
      IF (NSIMC.EQ.0) CALL ERRWRK('NEONEX',(MWRK0+MWRK1),LFREE)
C
C     **********
C     ********** Loop over batches:
C     **********
C
      JSIMX  = 0
      NSIMX  = NSIMC
      NSIML  = NROOTC
      IBATCH = 0
C     ********** REPEAT UNTIL (ALL BATCHES DONE)
   10 CONTINUE
         IBATCH = IBATCH + 1
C
C        **********
C        ********** Calculate residual vectors.
C        **********
C
C        ********** a) core allocation
C          RVCS for residual vectors
C          WKAP for K(oc) * y(c) (if KAPOPT)
C               and calculation of orb. trial vectors (always).
C
         KRVCS = 1
         KWKAP = KRVCS + NSIMX*NVARX
         KTEMP = KWKAP + NSIMX*NWOPT
         IF (OLSEN1) THEN
            KXVECS = KTEMP
            KSC    = KXVECS + NSIMX*NCONF
         ELSE
            KXVECS = KRVCS
            KSC    = KTEMP
         END IF
         KSO   = KSC   + NCONFX
         CALL DZERO (WRK(KRVCS),(KSC-KRVCS))
C
C        ********** b) calculate RVCS except for K(oc) * y(c)
C                   which is saved in WKAP
C
C        ********** sigma vector part
C                   (first sigma vector is BETA*G)
C
         REWIND LUIT5
         IF (NCONFX.GT.0) READ (LUIT5)
         IF (NWOPT.GT.0)  READ (LUIT5)
         JRVCS = KRVCS
         DO 200 ISIMX = JSIMX+1,JSIMX+NSIMX
            FAC   = BETA*EVEC(1,JROOT(ISIMX))
            CALL DAXPY(NVARX,FAC,G(IGX),1,WRK(JRVCS),1)
            JRVCS = JRVCS + NVARX
  200    CONTINUE
            IF (IPRI6 .gt. 210) THEN
               WRITE (LUPRI,'(/A,I3,A)')
     &           'Printing',NSIMX,' gradient part of residual vectors:'
               WRITE (LUPRI,*) 'FAC, beta, evec(1,1)',
     &                          fac,beta,evec(1,jroot(jsimx+1))
               CALL OUTPUT(WRK(KRVCS),1,NVAR,1,NSIMX,
     &              NVAR,NSIMX,-1,LUPRI)
            END IF
C
         DO 220 IREDL = 2,NREDL
            CITYPE = (IBNDX(IREDL) .EQ. JBCNDX)
            JRVCS  = KRVCS
            IF (NCONFX.GT.0) CALL READT(LUIT5,NCONFX,WRK(KSC))
            IF (NWOPT.GT.0) THEN
               CALL READT(LUIT5,NWOPT,WRK(KSO))
               IF (KAPOPT) JWKAP = KWKAP
            END IF
            DO 210 ISIMX = JSIMX+1,JSIMX+NSIMX
               FAC = EVEC(IREDL,JROOT(ISIMX))
               IF (CITYPE .AND. KAPOPT) THEN
                  CALL DAXPY(NCONF,FAC,WRK(KSC),1,WRK(JRVCS),1)
                  CALL DAXPY(NWOPT,FAC,WRK(KSO),1,WRK(JWKAP),1)
                  JWKAP = JWKAP + NWOPT
               ELSE
                  CALL DAXPY(NVARX,FAC,WRK(KSC),1,WRK(JRVCS),1)
               END IF
               JRVCS = JRVCS + NVARX
  210       CONTINUE
  220    CONTINUE
C
C        ********** add WKAP to RVCS to finish sigma vector part
C
         JCRV  = KRVCS
         JORV  = JCRV + NCONFX
         IF (KAPOPT) JWKAP = KWKAP
         DO 230 ISIMX = JSIMX+1,JSIMX+NSIMX
            IF (KAPOPT) THEN
               CALL DAXPY(NWOPT,D1,WRK(JWKAP),1,WRK(JORV),1)
               JWKAP = JWKAP  + NWOPT
            END IF
            IF (JROOT(ISIMX) .EQ. ISTATE .AND. P4FLAG(6)) THEN
               IF(NCONFX.GT.0) YCNRM=DNORM2(NCONF,WRK(JCRV),1)
               IF(NWOPT.GT.0)  YONRM=DNORM2(NWOPT,WRK(JORV),1)
               YNRM  = SQRT(YCNRM**2+YONRM**2)
            END IF
            JCRV  = JCRV + NVARX
            JORV  = JORV + NVARX
  230    CONTINUE
C
C        ******** subtract E times eigenvector to finish residual.
C                 Skip first trial vector: CREF, CREF contributions
C                  always add up to zero. Furthermore, beware that
C                  CREF has been projected out from all SIGMA vectors.
C                 Check for later if we have any orbital trial vectors
C                  (or NWOPT .EQ. 0), because then we should include
C                  an orbital trial vector.
C
            IF (IPRI6 .gt. 210) THEN
               WRITE (LUPRI,'(/A,I3,A)')
     &           'Printing',NSIMX,' sigma part of residual vectors:'
               CALL OUTPUT(WRK(KRVCS),1,NVAR,1,NSIMX,
     &              NVAR,NSIMX,-1,LUPRI)
            END IF
         KB = KSC
         REWIND LUIT3
         READ (LUIT3)
         IF (NWOPT .EQ. 0) THEN
            FORCEO = .FALSE.
         ELSE
            FORCEO = .TRUE.
         END IF
         DO 330 IREDL = 2,NREDL
            CITYPE = (IBNDX(IREDL) .EQ. JBCNDX)
            IF (CITYPE) THEN
               CALL READT(LUIT3,NCONF,WRK(KB))
               JRV = KRVCS
               LBV = NCONF
            ELSE
               FORCEO = .FALSE.
               CALL READT(LUIT3,NWOPT,WRK(KB))
               JRV = KRVCS + NCONFX
               LBV = NWOPT
            END IF
            DO 310 ISIMX = JSIMX+1,JSIMX+NSIMX
               JROOTX = JROOT(ISIMX)
               FAC = - EVAL(JROOTX) * EVEC(IREDL,JROOTX)
               CALL DAXPY(LBV,FAC,WRK(KB),1,WRK(JRV),1)
               JRV = JRV + NVARX
  310       CONTINUE
  330    CONTINUE
C
         IF (P6FLAG(43)) THEN
            WRITE (LUPRI,'(/A,(10I4))')
     *         ' (NEONEX) residual vectors for roots',
     *         (JROOT(ISIMX),ISIMX=JSIMX+1,JSIMX+NSIMX)
            CALL OUTPUT(WRK(KRVCS),1,NVARX,1,NSIMX,NVARX,NSIMX,-1,LUPRI)
         END IF
         IF (KAPOPT .AND. P6FLAG(44)) THEN
            WRITE (LUPRI,'(/A)')
     *         ' (NEONEX) o-c coupling of residual vectors:'
            CALL OUTPUT(WRK(KWKAP),1,NWOPT,1,NSIMX,NWOPT,NSIMX,-1,LUPRI)
         END IF
C
C        **********
C        ********** Convergence check, contract RVCS
C        **********
C
         IF (TSTCNV) THEN
            THRQN = CNVTHR(LINCNV)
            THRQX = THQMIN
            IF (P4FLAG(6)) THEN
            IF (LINCNV) THEN
               WRITE(LUW4,'(/A,F15.8)')
     *         ' (NEONEX) Linear (global) convergence threshold:',THRQN
            ELSE
               WRITE(LUW4,'(/A,F15.8)')
     *         ' (NEONEX) Quadratic (local) convergence threshold:',
     *         THRQN
            END IF
            END IF
         ELSE
C        if (not test convergence) only exit if space spanned
            THRQN = TLINDP
            THRQX = TLINDP
         END IF
         IF (P4FLAG(6)) WRITE (LUW4,9010) NREDL
 9010    FORMAT(/' Residual norms when dim(red L) =',I4,
     *          /' NEO root     CSF        orbital          total'/)
 9020    FORMAT(I5,3F15.8,2A)
 9022    FORMAT(5X,3F15.8,' sigma part of previous vector')
 9025    FORMAT(/' Residual norm when dim(red L) =',I4,
     *          /' NEO root     CSF        orbital          total',
     *          /I5,3F15.8,A)
C
C        IBNDX(IBNOFF+*) is used temporarily to keep track of trial
C                        vector type: CI type or orbital type, or both,
C                        or none.
C                        bit 1 set: CI type
C                        bit 2 set: orbital type
C
C
         IF (KAPOPT) KTHKAP = KB
         JCRV  = KRVCS
         JORV  = JCRV + NCONFX
         NCVEC = 0
         NOVEC = 0
         IBNOFF = NBVIT3
         DO 1100 ISIMX = JSIMX+1,JSIMX+NSIMX
            JROOTX = JROOT(ISIMX)
            IF (NCONFX.GT.0) QCNRM = DNORM2(NCONF,WRK(JCRV),1)
            IF (NWOPT.GT.0)  QONRM = DNORM2(NWOPT,WRK(JORV),1)
            QNRM  = SQRT(QCNRM**2+QONRM**2)
C
            IF (JROOTX .EQ. ISTATE) THEN
C              renormalize so residual same as if we
C              had solved level-shifted N-R linear equations.
               FAC   = D1 / ABS(BETA*EVEC(1,ISTATE))
               YCNRM = FAC * YCNRM
               YONRM = FAC * YONRM
               YNRM  = FAC * YNRM
               QCNRM = FAC * QCNRM
               QONRM = FAC * QONRM
               QNRM  = FAC * QNRM
               IF (QNRM .LT. THRGRD) THEN
C                 MC optimization may be converged
                  THRQNX = THRQN
               ELSE IF (QCNRM .LT. DP1*THRQN) THEN
C                 exclude case where next check is not appropriate :
C                 if CSF converged (e.g. for Hartree-Fock)
                  THRQNX = THRQN
               ELSE IF (QONRM .GT. QCNRAT*QCNRM) THEN
C                 no convergence check if orbital residual
C                 bigger than QCNRAT of CSF residual, this leads
C                 to MUCH greater stability in the optimization.
                  THRQNX = TLINDP
               ELSE IF (QNRM  .GT. ABSCNV*THRQN   .AND.
     *                  QCNRM .GT. QCNFAC*QCNOLD) THEN
C                 if residual in quadratic "grey" area then
C                 no convergence check if CI residual has increased
C                 by a factor QCNFAC, because this means strong,
C                 unaccounted for coupling between orbital and CSF.
                  THRQNX = TLINDP
               ELSE
                  THRQNX = THRQN
               END IF
               QCNOLD = QCNRM
            ELSE
C              not ISTATE, sharp convergence not necessary
               THRQNX = THRQX
            END IF
C
C           convergence check
C
            IF (.NOT. TSTCNV) THRQNX = TLINDP
            IF (QNRM .LE. THRQNX) THEN
C              Converged, if we should check for convergence,
C              or norm equal zero (space spanned)
               IF (P4FLAG(6)) THEN
                  WRITE (LUW4,9020) JROOTX,QCNRM,QONRM,QNRM,' converged'
               ELSE IF (IPRI4 .GT. 2) THEN
                  WRITE (LUW4,9025) NREDL,JROOTX,QCNRM,QONRM,QNRM,
     *               ' converged'
               END IF
               IF (JROOTX .EQ. ISTATE) THEN
C                 we have finished
                  JCONV  = 0
                  NROOTC = 0
                  WRITE (LUW4,9030)
                  GO TO 8100
               ELSE
                  IF (.NOT.P4FLAG(6) .AND. IPRI4 .GT. 2)
     *               WRITE (LUW4,9040) JROOTX
                  JROOT(ISIMX) = 0
                  IBNDX(IBNOFF+ISIMX) = 0
               END IF
            ELSE
               IF (P4FLAG(6)) THEN
                  IF (JROOTX .EQ. ISTATE) THEN
                     WRITE (LUW4,9020) JROOTX,QCNRM,QONRM,QNRM,' ***'
                     WRITE (LUW4,9022) YCNRM,YONRM,YNRM
                  ELSE
                     WRITE (LUW4,9020) JROOTX,QCNRM,QONRM,QNRM
                  END IF
               END IF
               IF (JROOTX .EQ. ISTATE .AND. FORCEO) THEN
C                 assure first trial vector is an orbital vector
                  NOVEC = NOVEC + 1
                  IBNDX(IBNOFF+ISIMX) = 2
                  IF (KAPOPT) THEN
                     WRK((KTHKAP-1)+NOVEC) =
     *               MAX(FKAPC*QCNRM, FKAPO*QONRM)
                  END IF
                  FACTRN = 1.0 D-3
                  CTEST = FACTRN * QONRM
                  DO I = JORV+1,JORV+NWOPT
                     IF (ABS(WRK(I)) .LT. CTEST) WRK(I) = D0
                  END DO
               ELSE IF (QCNRAT*QCNRM.GT.QONRM) THEN
C                 CI trial vector, prepare for Davidson:
                  NCVEC = NCVEC + 1
                  IBNDX(IBNOFF+ISIMX) = 1
                  FACTRN = 1.0 D-3
                  CTEST = FACTRN * QCNRM
                  DO I = JCRV+1,JCRV+NCONF
                     IF (ABS(WRK(I)) .LT. CTEST) WRK(I) = D0
                  END DO
               ELSE
C                 orbital trial vector
                  NOVEC = NOVEC + 1
                  IBNDX(IBNOFF+ISIMX) = 2
                  IF (KAPOPT) THEN
                     WRK((KTHKAP-1)+NOVEC) =
     *               MAX(FKAPC*QCNRM, FKAPO*QONRM)
                  END IF
                  FACTRN = 1.0 D-3
                  CTEST = FACTRN * QONRM
                  DO I = JORV+1,JORV+NWOPT
                     IF (ABS(WRK(I)) .LT. CTEST) WRK(I) = D0
                  END DO
               END IF
            END IF
            JCRV = JCRV + NVARX
            JORV = JORV + NVARX
 1100    CONTINUE
C
 9030    FORMAT(/' (NEONEX) NEO vector is converged.')
 9040    FORMAT(/' (NEONEX) NEO auxiliary root no.',I3,
     *           ' is converged.')
C
C        Use IBNDX(IBNOFF+NSIMX+*) for keeping track of which root
C        goes with each trial vector
C
         KBCVCS = KRVCS
         KBOVCS = KBCVCS + NCVEC*NCONFX
         JCRV   = KRVCS
         JORV   = JCRV + NCONFX
         IF (KAPOPT) JWKAP  = KWKAP
         JCRVX  = KBCVCS
         JORVX  = KWKAP
         ICSIMX = IBNOFF + NSIMX
         IOSIMX = ICSIMX + NCVEC
         DO 1200 ISIMX = JSIMX+1,JSIMX+NSIMX
            JROOTX = JROOT(ISIMX)
            IBTYPE = IBNDX(IBNOFF+ISIMX)
            IF (IBTAND(IBTYPE,1) .NE. 0) THEN
               ICSIMX = ICSIMX + 1
               IBNDX(ICSIMX) = JROOTX
               IF (JCRV.NE.JCRVX) THEN
                  CALL DCOPY(NCONF,WRK(JCRV),1,WRK(JCRVX),1)
               END IF
               JCRVX = JCRVX + NCONF
            END IF
            IF (IBTAND(IBTYPE,2) .NE. 0) THEN
               IOSIMX = IOSIMX + 1
               IBNDX(IOSIMX) = JROOTX
CHJ-860921: currently KAPOPT only works correctly for ISTATE
C              IF (KAPOPT) THEN
               IF (KAPOPT .AND. JROOTX .EQ. ISTATE) THEN
                  FAC = D1 / (BETA * EVEC(1,JROOTX))
                  DO 1280 K = 1,NWOPT
                     WRK((JORVX-1)+K) = FAC*WRK((JWKAP-1)+K)+G(NCONF+K)
 1280             CONTINUE
               ELSE
                  CALL DCOPY(NWOPT,WRK(JORV),1,WRK(JORVX),1)
               END IF
               JORVX = JORVX + NWOPT
            END IF
            JCRV = JCRV + NVARX
            JORV = JORV + NVARX
            IF (KAPOPT) JWKAP = JWKAP + NWOPT
 1200    CONTINUE
C
C
C
         IF (NWOPT.GT.0) THEN
            CALL DCOPY((NOVEC*NWOPT),WRK(KWKAP),1,WRK(KBOVCS),1)
            KTEMP = KBOVCS + NOVEC*NWOPT
            IF (KAPOPT) THEN
               CALL DCOPY(NOVEC,WRK(KTHKAP),1,WRK(KTEMP),1)
               KTHKAP = KTEMP
               KTEMP  = KTHKAP + NOVEC
            END IF
         ELSE
            KTEMP = KBOVCS
         END IF
C
C        **********
C        ********** Use Davidson's shift to get CI trial vectors
C        **********
C
C        D = CDIAG(I) + DSHIFT = 1/2 * [ 2*(CIDIAG(I)-EACTVN) - EVAL ]
C          = 1/2 * [ LDIAG(I) - EVAL ]
C
C        hjaaj aug 2006: new default of 0.01 for ESHIFT
         IF (ESHIFT.NE.0.01D0) WRITE (LUW4,9120) ESHIFT
 9120    FORMAT(/' (NEONEX) Davidson level shift shifted by ESHIFT =',
     *          F20.15)
C
C
C
         IF (NCVEC.GT.0) THEN
            KCDIAG = KTEMP
            IF (.NOT.DOPHP) THEN
               LCDIAG = NVAR
            ELSE
               LCDIAG = LPHPMX
            END IF
            KWPHP = KCDIAG + LCDIAG
            LWPHP = LFREE + 1 - KWPHP
            CALL SRWPHP(WRK(KCDIAG),LCDIAG,.FALSE.)
C           CALL SRWPHP(PHPVEC,LPHPMX,PHPWRT)
            IF (BETM1.NE.D0) THEN
C
C              Calculate beta-dependent part of diagonal for beta .ne.1
C              correction. (It is scaled by 1/2 compared to the formula
C              because CI L diagonal is scaled by 1/2 -- the true CI
C              block of the L matrix is twice the CI matrix).
C
               DO 1600 I = 1,NCONF
                  WRK((KCDIAG-1)+I) = WRK((KCDIAG-1)+I)
     *               + BETM1 * CREF(I) * G(I)
 1600          CONTINUE
            END IF
C
C           find lowest diagonal value ? maybe thru common from GRAD ?
C           in common CDLOW = min(CDIAG) and ICDLOW ?
C           any "DLOW" as in old NEXT to find extra "ESHIFT"?
            ECSHFT = DHALF*ESHIFT - EACTVN
            ICRV   = KBCVCS - 1
            ICXOFF = IBNOFF + NSIMX
            DO 2400 ICX = ICXOFF + 1, ICXOFF + NCVEC
               JROOTX = IBNDX(ICX)
            IF (.NOT. DOPHP) THEN
               DSHIFT = ECSHFT - DHALF * EVAL(JROOTX)
               DO 2200 ICONF = 1,NCONF
                  D = WRK(KCDIAG-1+ICONF) + DSHIFT
                  IF (D.LT.DTEST) THEN
                     IF (D.GT.(-DTEST)) D = SIGN(DTEST,D)
                  END IF
                  WRK(ICRV+ICONF) = WRK(ICRV+ICONF) / D
 2200          CONTINUE
            ELSE
               DSHIFT = DHALF * EVAL(JROOTX) - ECSHFT
CPHP           IF (OLSEN1) THEN
CPHP              JXVEC = KXVECS + (ICVEC-1)*NCONF
CPHP           ELSE
                  JXVEC = 999 999 999
CPHP           END IF
               IPRPHP = MAX(0,IPRI6-5)
               CALL NEXCI(OLSEN1,DSHIFT,NCONF,WRK(ICRV+1),WRK(JXVEC),
     *                  WRK(ICRV+1),WRK(KCDIAG),IPRPHP,WRK(KWPHP),LWPHP)
C              CALL NEXCI(OLSEN,ENER,NCVAR,D,XVEC,
C    &                    RES,DIAG,IPRPHP,WRK,LWRK)
            END IF
               ICRV = ICRV + NCONF
 2400       CONTINUE
         END IF
C
C        **********
C        **********
C        **********
C
         IF (NOVEC.GT.0) THEN
            KODIAG = KTEMP
            KWO    = KODIAG + NWOPT
            REWIND LUIT2
C           Find label ORBDIAG, for diagonal of orbital Hessian.
            CALL MOLLAB(LABIT2(3),LUIT2,lupri)
            CALL READT(LUIT2,NWOPT,WRK(KODIAG))
C
C-----------------------
C CBN+JK 03.01.06
C----------------------
            IF (SOLVNT .OR. PCM .OR. QM3) THEN
C-----------------------
C CBN+JK 03.01.06
C----------------------
               CALL MOLLAB(LABIT2(4),LUIT2,lupri)
               IF (NCONF .GT. 1) READ (LUIT2)
               CALL READT (LUIT2,NWOPT,WRK(KWO))
               CALL DAXPY (NWOPT,D1,WRK(KWO),1,WRK(KODIAG),1)
            END IF
C
C           (find lowest diagonal value ? maybe thru common from GRAD ?
C            in common ODLOW = min(ODIAG) and IODLOW ?
C            any "DLOW" as in old NEXT to find extra "ESHIFT"?
C            If (KAPOPT) and any extra shift is wanted, then
C            add it to DLAM below. ) 850308/hjaaj-comment
            IOXOFF = IBNOFF + NSIMX + NCVEC
#if defined (VAR_OLDCODE)
            IF (.NOT.KAPOPT) THEN
               NOVECX = NOVEC
            ELSE
               KDLAM = KWO
               KWKAP = KDLAM + NOVEC
               LWKAP = LFREE - KWKAP
#if defined (VAR_KAPNEO)
               IF (KAPNEO .AND. NROOTC.EQ.1 .AND. ISTATE.EQ.1) THEN
                  NRSEQ = 0
                  DAMP  = BETA
                  STPCI = D0
                  DO 3200 IREDL = 2,NREDL
                     IF (IBNDX(IREDL).EQ.JBCNDX) STPCI =
     *                  STPCI + EVEC(IREDL,1)*EVEC(IREDL,1)
 3200             CONTINUE
                  STPCI = STPCI / ABS(EVEC(1,1) * BETA)
                  IF (STPCI .GE. 0.8D0*RTRUST) THEN
C                    ignore KAPNEO this time
                     NRSEQ = 1
                     WRK(KDLAM) = EVAL(1)
                     ROTRST = D0
                  ELSE
                     ROTRST = SQRT(RTRUST*RTRUST - STPCI*STPCI)
                  END IF
                  NOVECX = 0
               ELSE
#endif
                  NRSEQ  = 1
                  ROTRST = D0
                  DO 3400 IOX = IOXOFF + 1, IOXOFF + NOVEC
                     JROOTX = IBNDX(IOX)
CHJ-860921: currently KAPOPT onlyt works correctly for ISTATE
C                    WRK(KDLAM+(IOX-1)) = EVAL(JROOTX)
                     IF (JROOTX .EQ. ISTATE) THEN
                        WRK(KDLAM)  = EVAL(ISTATE)
                        WRK(KTHKAP) = WRK(KTHKAP+(IOX-1))
                        IOXSAV = IOX
                        NOVECX = NOVEC - 1
                        GO TO 3450
                     END IF
 3400             CONTINUE
                  NOVECX = NOVEC
                  GO TO 3500
#if defined (VAR_KAPNEO)
               END IF
#endif
 3450          CONTINUE
C              IF (IPRI6 .GE. 25) WRITE (LUPRI,'(//A,1P,D10.2)')
C    *            ' (NEONEX) calling NEXKAP, THRKAP(ISTATE)=',
C    *            WRK(KTHKAP+(ISTATE-1))
C              CALL NEXKAP(NOVEC,WRK(KBOVCS),NBVIT3,LUIT3,LUIT5,IBNDX,
C    *                     NRSEQ,ROTRST,WRK(KTHKAP),WRK(KDLAM),DAMP,
C    *                     WRK(KODIAG),CMO,G(1+NCONF),DV,PV,FC,FV,
C    *                     WRK(KWKAP),LWKAP)
               IF (IPRI6 .GE. 25) WRITE (LUPRI,'(//A,1P,D10.2)')
     *            ' (NEONEX) calling NEXKAP for ISTATE, THRKAP =',
     *            WRK(KTHKAP)
               CALL NEXKAP(1,WRK(KBOVCS+(IOXSAV-1)*NWOPT),
     *                     NBVIT3,LUIT3,LUIT5,IBNDX,NRSEQ,ROTRST,
     *                     WRK(KTHKAP),WRK(KDLAM),DAMP,
     *                     WRK(KODIAG),CMO,G(1+NCONF),DV,PV,FC,FV,
     *                     WRK(KWKAP),LWKAP)
C              CALL NEXKAP(NOSIM,GDM,NBFIL3,IFIL3,IFIL5,IBNDX,
C    *                     NRSEQ,ROTRST,THRKAP,EVAL,DAMP,DIAOR,
C    *                     CMO,GORB,DV,PV,FC,FV,WRK,LFREE)
            END IF
C
C           make NOVECX Davidson/conj.grad.inverse it.
C           type trial vectors
C
 3500       IF (NOVECX .GT. 0) THEN
               IORV   = KBOVCS - 1
               IODIAG = KODIAG - 1
               DO 3800 IOX = IOXOFF + 1, IOXOFF + NOVEC
                  JROOTX = IBNDX(IOX)
                  IF (KAPOPT .AND. JROOTX .EQ. ISTATE) GO TO 3700
                  DSHIFT = ESHIFT - EVAL(JROOTX)
                     DO 3600 IWOPT = 1,NWOPT
                        D = WRK(IODIAG+IWOPT) + DSHIFT
                        IF (D.LT.DTEST) THEN
                           IF (D.GT.(-DTEST)) D = SIGN(DTEST,D)
                        END IF
                        WRK(IORV+IWOPT) = WRK(IORV+IWOPT) / D
 3600                CONTINUE
 3700             IORV = IORV + NWOPT
 3800          CONTINUE
            END IF
#endif
#if !defined (VAR_OLDCODE)
            IORV   = KBOVCS - 1
            IODIAG = KODIAG - 1
            DO 3800 IOX = IOXOFF + 1, IOXOFF + NOVEC
               JROOTX = IBNDX(IOX)
               IF (KAPOPT .AND. JROOTX .EQ. ISTATE) THEN
C
C              improved orbital trial vector
C
                  LWO    = LFREE - KWO
                  ROTRST = D0
                  DLAM   = EVAL(ISTATE)
                  THRKAP = WRK(KTHKAP+(IOX-(IOXOFF+1)))
                  IF (IPRI6 .GE. 25) WRITE (LUPRI,'(//A,1P,D10.2)')
     *               ' (NEONEX) calling NEXKAP for ISTATE, THRKAP =',
     *               THRKAP
                  CALL NEXKAP(1,WRK(IORV+1),NBVIT3,LUIT3,LUIT5,IBNDX,
     *                        1,ROTRST,THRKAP,DLAM,DAMP,WRK(KODIAG),
     *                        CMO,G(1+NCONF),DV,PV,FC,FV,
     &                        WRK(KWO),LWO)
C                 CALL NEXKAP(NOSIM,GDM,NBFIL3,IFIL3,IFIL5,IBNDX,
C    *                        NRSEQ,ROTRST,THRKAP,EVAL,DAMP,DIAOR,
C    *                        CMO,GORB,DV,PV,FC,FV,WRK,LFREE)
               ELSE
C
C              Davidson/conj.grad.inverse it. type trial vectors
C
                  DSHIFT = ESHIFT - EVAL(JROOTX)
                  DO 3600 IWOPT = 1,NWOPT
                     D = WRK(IODIAG+IWOPT) + DSHIFT
                     IF (D.LT.DTEST) THEN
                        IF (D.GT.(-DTEST)) D = SIGN(DTEST,D)
                     END IF
                     WRK(IORV+IWOPT) = WRK(IORV+IWOPT) / D
 3600             CONTINUE
               END IF
               IORV = IORV + NWOPT
 3800       CONTINUE
#endif
         END IF
C
C        Do microcanonical average on new BOVECS
C        (AVERAG returns immediately if no averageing requested)
C
         IF (NOVEC.GT.0) THEN
            CALL AVERAG(WRK(KBOVCS),NWOPT,NOVEC)
         END IF
C
C        **********
C        ********** Orthogonalize these NCVEC+NOVEC new b-vectors
C        ********** against previous b-vectors and each other.
C        **********
C
C        ORTBVC also defines IBNDX for the new trial vectors
C               and updates NBPREV.
C        LUIT3 is positioned, ready for the new trial vectors.
C
         IF (SYMTRZ .OR. CTRUNC) THEN
C           ... symmetrize CI vector so that numerical
C               inaccuracy won't break symmetry;
C               WARNING : if symmetry were broken in residual vectors,
C                         and thus in sigma vectors, then this will
C                         lead to constant residual.
C
C           ... increase efficiency by zeroing elements in
C               trial vectors abs less than 1.d-3 abs largest element.
C               Do not truncate orbital trial vectors as that
C               will destroy KAPOPT performance.
            FACTRN = 1.0 D-3
            ICRV   = KBCVCS - 1
            DO 4400 ICX = 1,NCVEC
            IF (CTRUNC) THEN
               ICMAX = IDAMAX(NCONF,WRK(ICRV+1),1)
               CTEST = FACTRN * ABS(WRK(ICRV+ICMAX))
               DO 4300 ICONF = ICRV+1,ICRV+NCONF
                  IF (ABS(WRK(ICONF)) .LT. CTEST) WRK(ICONF) = D0
 4300          CONTINUE
               ICRV = ICRV + NCONF
            END IF
            IF (SYMTRZ) THEN
               CNORM = DNORM2(NCONF,WRK(ICRV+1),1)
               IF (CNORM .GT. TLINDP) THEN
                  CNORM = D1 / CNORM
                  CALL DSCAL(NCONF,CNORM,WRK(ICRV+1),1)
                  CALL VSYMTR(NCONF,WRK(ICRV+1),THRSIM)
               END IF
            END IF
 4400       CONTINUE
         END IF
C
         KBPREV = KTEMP
         KTT0   = KBPREV + MAX(NCONF,NWOPT)
C
         NDMBC  = MAX(1,NCONF)
         NDMBO  = MAX(1,NWOPT)
         CALL ORTBVC(NCVEC,NOVEC,NDMBC,NDMBO,NBPREV,IBNDX,LUIT3,
     *               WRK(KBCVCS),WRK(KBOVCS),THRLDV,WRK(KTT0),
     *               WRK(KBPREV))
C        CALL ORTBVC(NBC,NBO,NDMBC,NDMBO,NBPREV,IBNDX,LUBVC,
C    *               BCVECS,BOVECS,THRLDP,TT0,BPREV)
C
C        **********
C        ********** Check if finished all NROOTC residuals,
C        ********** if not go back to 10 for next batch.
C        **********
C
         NCSIM = NCSIM + NCVEC
         NOSIM = NOSIM + NOVEC
         NSIML = NSIML - NSIMX
         JSIMX = JSIMX + NSIMX
         NSIMX = MIN(NSIMC,NSIML)
      IF (NSIMX.GT.0) GO TO 10
C     ^--- END REPEAT UNTIL (ALL BATCHES DONE, i.e. NSIMX = 0)
C
C     **********
C     ********** If "in memory", let BOVECS follow BCVECS w/o gap
C     **********
C
      IF (IBATCH .EQ. 1) THEN
         KBONEW = KBCVCS + NCSIM*NCONF
         IF (KBONEW.LT.KBOVCS .AND. NOSIM.GT.0) THEN
            DO 5000 I = 0,NOSIM*NWOPT-1
               WRK(KBONEW+I) = WRK(KBOVCS+I)
 5000       CONTINUE
            KBOVCS = KBONEW
         END IF
      ELSE
         NCSIM = -NCSIM
         NOSIM = -NOSIM
      END IF
C
C     **********
C     ********** Update NROOTC and JROOT by removing converged
C     ********** vectors
C     **********
C
      N = 0
      DO 8000 ISIM = 1,NROOTC
         IF (JROOT(ISIM).NE.0) THEN
            N = N + 1
            JROOT(N) = JROOT(ISIM)
         END IF
 8000 CONTINUE
      NROOTC = N
C
C     **********
C     ********** End of subroutine NEONEX
C     **********
C
 8100 CONTINUE
      CALL QEXIT('NEONEX')
      RETURN
      END
C  /* Deck cnvthr */
      FUNCTION CNVTHR(LINCNV)
C
C  3-Nov-1986 Hans Joergen Aa. Jensen
C
C  Function for evaluating a reasonable convergence threshold
C  for NEO and Newton-Raphson (NR) equations.  The residual is
C  related to a prediction for the gradient at the next set of MC
C  parameters and this is used to estimate when the quadratic region
C  has been traversed and any further accuracy in the solution vector
C  will be a waste of effort.
C
#include <implicit.h>
      LOGICAL LINCNV
C
      PARAMETER (THQCXX = 0.10D0, RTTHRG = 4.00D0)
      PARAMETER (THGFAC = 0.98D0, THGFC2 = 0.70D0)
C
C Used from common blocks:
C   INFINP : THRGRD
C   INFOPT : GRDNRM,GNORM(*),THQLIN,THQKVA
C
#include <maxorb.h>
#include <infinp.h>
#include <infopt.h>
C
C      CNVTHR  = qa) dynamic threshold from 8-Apr-1984 /hj
C                    idea: need to converge down below GRDNRM squared
C                          to get quadratic convergence, further
C                          convergence is wasted.
C
C                qb) never less than THGFAC = 0.98 times what has
C                    been requested as total convergence of gradient,
C                    THGFAC is a safety factor for numerical
C                    inaccuracy (860921/hj)
C
C                l ) globally only use linear convergence
C
      IF (LINCNV) THEN
         THRQN = THQLIN*GRDNRM
      ELSE
C
CHJ/PJ-860923: After meticulous study of H2O and H2O2 we decided to use
C              these convergence thresholds, and hope they carry over.
C              We noticed the quadratic factor is usually smaller
C              for the CSF gradient than for the orbital gradient,
C              except in the local region. Thus we divide the
C              quadratic factor estimate, THQKVA, by THQCXX for the
C              CSF part.
C              If THRQN ends within RTTHRG times THRGRD, we take
C              the chance an try to converge completely in this
C              macro. We converge to THGFC2 times THRGRD because it
C              is likely the gradient in next iteration will be
C              higher than the residual.
C
         GCNRM = GNORM(1)
         GONRM = GNORM(2)
         THRQN = THQKVA * (GCNRM*GCNRM)
         IF (GCNRM .GE. 5.D-3) THRQN = THQCXX * THRQN
         THRQN = MAX( THRQN, THQKVA*(GONRM*GONRM) )
         THRQN = MIN( THRQN, THQLIN*GRDNRM )
         IF (THRQN .LE. THRGRD) THEN
            THRQN = THGFAC*THRGRD
         ELSE IF (THRQN .LE. RTTHRG*THRGRD) THEN
            THRQN = THGFC2*THRGRD
         END IF
      END IF
      CNVTHR = THRQN
      RETURN
C     end of function CNVTHR.
      END
C  /* Deck neored */
      SUBROUTINE NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,
     *                   G,IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
C
C Written 30-Mar-1984 by Hans Jorgen Aa. Jensen
C
C Revisions:
C   9-Apr-1984 hjaaj
C  31-Oct-1984 hjaaj (use new b-vec. directly from BVECS)
C
C Input:
C  ICTL,  flow control
C         =1, (BIT 1 set) extend the reduced L-matrix
C         =2, (BIT 2 set) find the new eigenvalues and -vectors
C         =3, (BIT 1 and BIT 2 set) both
C  REDL,  the old reduced L-matrix (dimension: NREDL-N)
C  N,     number of new b-vectors
C  BVECS, the N new b-vectors (is destroyed in NEORED)
C  SVECS, the sigma vectors of the N new b-vectors
C  G,     the gradient
C  BETVAL,the beta value to use when solving eigenvalue eq.
C
C  (The reduced L-matrix is the projection of the L-matrix
C   on the basis of b-vectors.)
C
C Output:
C  REDL, the new, extended reduced L-matrix (dimension: NREDL)
C  EVEC, eigenvectors of the new reduced L-matrix
C  EVAL, eigenvalues of the new reduced L-matrix
C
C Scratch:
C  BVECS, is used for B vectors from unit LUIT3 (section 1)
C         and for the reduced L-matrix (section 2)
C         dimension: NVAR and NREDL*(NREDL+1)/2+2*NREDL, rsp.
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION BVECS(*),SVECS(NVAR,*), REDL(*),IBNDX(*),G(*)
      DIMENSION EVEC(NREDL,*),EVAL(*)
#include <ibndxdef.h>
C
C Used from common blocks:
C   INFVAR : NCONF,NWOPT
C   INFIND : IROW(*)
C   INFTAP : LUIT3
C
#include <maxash.h>
#include <maxorb.h>
#include <infind.h>
#include <inftap.h>
C
C local:
C
      LOGICAL EXTEND,SOLVE
C
#include <ibtdef.h>
C
      CALL QENTER('NEORED')
C
C Section 1: extend reduced L-matrix with N new b-vectors
C
      EXTEND = (IBTAND(ICTL,1) .NE. 0)
      SOLVE  = (IBTAND(ICTL,2) .NE. 0)
      NTSIM = NCSIM + NOSIM
      IF (EXTEND .AND. NTSIM.GT.0) THEN
         IF (NREDL .GE. LIROW)
     *      CALL QUIT('NEORED error: NREDL exceeds LIROW')
C        ... IROW(NREDL+1) used below
         IREDL = NREDL - NTSIM
         KBC = 1
         KBO = KBC + NCSIM*NCONF
C
C New b-vectors (are in BVECS)
C
         IBC = 1
         IBO = IBC + NCSIM
         DO 200 K = 1,NTSIM
            KL = IROW(IREDL+K) + IREDL
            IF (IBNDX(IREDL+K) .EQ. JBCNDX) THEN
               IK = IBC
               IBC = IBC + 1
            ELSE
               IK = IBO
               IBO = IBO + 1
            END IF
            ISTBC = KBC
            ISTBO = KBO
            DO 100 L=1,K
               IF (IBNDX(IREDL+L) .EQ. JBCNDX) THEN
                  REDL(KL+L) = DDOT(NCONF,BVECS(ISTBC),1,SVECS(1,IK),1)
                  ISTBC = ISTBC + NCONF
               ELSE
                  REDL(KL+L) =
     *               DDOT(NWOPT,BVECS(ISTBO),1,SVECS(NCONF+1,IK),1)
                  ISTBO = ISTBO + NWOPT
               END IF
  100       CONTINUE
  200    CONTINUE
C
C Old b-vectors.
C Rewind LUIT3
C
         IF (IREDL.GT.0) THEN
            REWIND LUIT3
C
C First overlap with s1 for first column (s1 = gradient G):
C
            READ (LUIT3)
            ISTBC = KBC
            ISTBO = KBO
            DO 300 ISIM = 1,NTSIM
               KL = IROW(IREDL+ISIM) + 1
               IF (IBNDX(IREDL+ISIM) .EQ. JBCNDX) THEN
                  REDL(KL) = DDOT(NCONF,BVECS(ISTBC),1,G,1)
                  ISTBC = ISTBC + NCONF
               ELSE
                  REDL(KL) = DDOT(NWOPT,BVECS(ISTBO),1,G(1+NCONF),1)
                  ISTBO = ISTBO + NWOPT
               END IF
  300       CONTINUE
C
            DO 15 L = 2,IREDL
               IBC = 1
               IBO = IBC + NCSIM
               IF (IBNDX(L) .EQ. JBCNDX) THEN
                  CALL READT(LUIT3,NCONF,BVECS)
                  DO 25 K=1,NTSIM
                     KL = IROW(IREDL+K) + L
                     IF (IBNDX(K) .EQ. JBCNDX) THEN
                        IK = IBC
                        IBC = IBC + 1
                     ELSE
                        IK = IBO
                        IBO = IBO + 1
                     END IF
                     REDL(KL) = DDOT(NCONF,BVECS,1,SVECS(1,IK),1)
   25             CONTINUE
               ELSE
                  CALL READT(LUIT3,NWOPT,BVECS)
                  DO 20 K=1,NTSIM
                     KL = IROW(IREDL+K) + L
                     IF (IBNDX(K) .EQ. JBCNDX) THEN
                        IK = IBC
                        IBC = IBC + 1
                     ELSE
                        IK = IBO
                        IBO = IBO + 1
                     END IF
                     REDL(KL) = DDOT(NWOPT,BVECS,1,SVECS(1+NCONF,IK),1)
   20             CONTINUE
               END IF
   15       CONTINUE
         END IF
      END IF
C
C *********************************************************
C Section 2: find eigenvalues and -vectors of reduced L-matrix
C
      IF (SOLVE) THEN
C
C Diagonalize reduced L.
C Copy REDL to BVECS for diagonalization and multiply first
C column by BETVAL.
C
         KMA = IROW(NREDL+1)
         CALL DCOPY(KMA,REDL(1),1,BVECS,1)
         I1 = 1
         DO 1050 I = 1,NREDL
            BVECS(I1) = BETVAL*BVECS(I1)
            I1 = I1 + I
 1050    CONTINUE
         CALL DUNIT(EVEC,NREDL)
         KBIG  = 1     + KMA
         KJBIG = KBIG  + NREDL
         LNEED = KJBIG + NREDL
         IF (LNEED .GT. LBVECS) CALL ERRWRK('NEORED.JACO',LNEED,LBVECS)
         CALL JACO(BVECS,EVEC,NREDL,NREDL,NREDL,
     *             BVECS(KBIG),BVECS(KJBIG))
         DO 1100 I=1,NREDL
            EVAL(I)=BVECS(IROW(I+1))
 1100    CONTINUE
         CALL ORDER (EVEC,EVAL,NREDL,NREDL)
      END IF
C
C *** End of subroutine NEORED
C
      CALL QEXIT('NEORED')
      RETURN
      END
C  /* Deck neost */
      SUBROUTINE NEOST (NCSIM,NOSIM,CREF,G,BVECS,LBVECS,IBNDX,DIAGL)
C
C Written 11-Feb-1985 Hans Joergen Aa. Jensen (was in SIRNEO)
C Revision 9-Jan-1986 hjaaj+pj (+DIAGL, will now also be called
C                               from SIRNR)
C
C Purpose:
C  NEO Start and NR start. Find the trial vectors for the first micro
C  iteration and define optimization parameters in common /infopt/.
C
C *** Construct the starting vectors for the first micro iteration
C     by orthogonalizing the CI vectors written on LUIT1 to the
C     normalized gradient fragments and to each other.
C
#include <implicit.h>
C
      DIMENSION CREF(*),G(*),BVECS(*),IBNDX(*),DIAGL(*)
C
#include <ibndxdef.h>
C
#include <dummy.h>
#include <thrldp.h>
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 )
      PARAMETER ( GLIM = 0.4D0 , DLIM = 0.01D0 , GTRLIM = 0.8D0 )
C
C Used from common blocks:
C   INFINP : ISTATE,LROOTS,NROOTS,IROOT(*)
C   INFVAR : NCONF,NWOPT,NVAR
C   INFOPT : NROOTC,JROOT(*)
C   INFDIM : MAXPHP,LPHPMX
C   INFTAP : LUIT1,LUIT3
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <infopt.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C
C
      LOGICAL BINMEM, OTRIAL, FIRST, GTRIAL, GETBVC, GETBVO
      CHARACTER*8 LABIT1, LABINF(2)
      DATA LABIT1/'STARTVEC'/
C
      CALL QENTER('NEOST ')
C
      NCOMAX = MAX(NCONF,NWOPT)
      GETBVC = (NCONF .GT. 1 .AND. GNORM(1) .GT. GLIM*GNORM(2))
      OTRIAL = FLAG(80)
C
C     Use gradient as first trial vector if global
C     (GRDNRM .GT. GTRLIM); otherwise assume NEO eigenvalue
C     close to Zero and divide by diagonal /920521-hjaaj
C
      GTRIAL = (GRDNRM .GT. GTRLIM)
C
C     (space for both CI and orbital trial vectors)
C
      MAXSIM = LBVECS / NVAR
      IF (MAXSIM .LE. 0) CALL ERRWRK('NEOST',NVAR,LBVECS)
C
C --- copy root information to JROOT
C     for use in micro iterations
C
      DO 10 I = 1,NROOTS
         JROOT(I) = IROOT(I)
   10 CONTINUE
      NROOTC = NROOTS
C
C
C
      NCONFX = MAX(4,NCONF)
      IF (NWOPT.EQ.0) THEN
         NWOPTX = 0
      ELSE
         NWOPTX = MAX(4,NWOPT)
      END IF
C
C     First basis vector for L-reduced is CREF.
C     Corresponding sigma vector is the gradient.
C
      REWIND (LUIT3)
      CALL WRITT(LUIT3,NCONFX,CREF)
      REWIND (LUIT5)
      IF (NCONF.GT.1)  CALL WRITT(LUIT5,NCONFX,G)
      IF (NWOPTX.GT.0) CALL WRITT(LUIT5,NWOPTX,G(1+NCONF))
      IBNDX(1) = JBCNDX
      NBPREV = 1
C
C
C
      NCSIM = 0
      NOSIM = 0
      NCVEC = 0
      NOVEC = 0
C
C     If LROOTS negative, this is first macro iteration and
C     we should take abs(LROOTS) trial vectors instead of
C     the NROOTS we use for simultaneous expansion.
C
      IF (LROOTS .GE. 0) THEN
         NTSIM  = MIN(MAXSIM,NROOTS)
         NROOTX = NROOTS
      ELSE
         LROOTS = -LROOTS
         NTSIM  = MIN(MAXSIM,LROOTS)
         NROOTX = LROOTS
      END IF
C
C
      KBC = 1
      KBO = KBC + NTSIM*NCONF
C
C
      ISTBC = KBC
      ISTBO = KBO
      IF (GETBVC) THEN
         IF (GTRIAL) THEN
            CALL DCOPY(NCONF,G,1,BVECS(KBC),1)
         ELSE
#if defined (VAR_NOPHP)
            IBC = KBC - 1
            DO 20 I = 1,NCONF
               IF (ABS(DIAGL(I)) .GT. DLIM) THEN
                  BVECS(IBC+I) = G(I) / DIAGL(I)
               ELSE
                  BVECS(IBC+I) = G(I) / DLIM
               END IF
   20       CONTINUE
#else
            IPRPHP = MAX(0,IPRI6-5)
            KWPHP  = KBC + NCONF
            LWPHP  = LBVECS + 1 - KWPHP
            CALL NEXCI(.FALSE.,D0,NCONF,BVECS(KBC),DUMMY,
     *                 G,DIAGL,IPRPHP,BVECS(KWPHP),LWPHP)
C           CALL NEXCI(OLSEN,ENER,NCVAR,D,XVEC,
C       &              RES,DIAG,IPRPHP,WRK,LWRK)
#endif
         END IF
         NCVEC = NCVEC + 1
         ISTBC = ISTBC + NCONF
      END IF
      GETBVO = ( NWOPT.GT.0 .AND.
     *   (OTRIAL .OR. .NOT.GETBVC .OR. FLAG(55) .OR. .NOT.FLAG(41)) )
      if (domcsrdft .and. NCVEC .gt. 0) GETBVO = .false.
C     861107-hj/pj: only include Gorb/Ldiag if Hartree-Fock
C                   or specifically requested (NEXKAP is more
C                   efficient for getting orbital trial vec.s)
C     920221-hjaaj: NCONF.eq.1 test ended with no trial vectors
C                   for NDET=2 case, therefore now NCONF.le.10
C     930323-hjaaj: NCONF.le.1 again; instead getbvo below if
C                   no trial vectors
C     941220-hjaaj: NCONF test is now via GETBVC
      IF (GETBVO) THEN
         IF (GTRIAL) THEN
            CALL DCOPY(NWOPT,G(NCONF+1),1,BVECS(KBO),1)
         ELSE
            IBO = KBO - NCONF - 1
            DO 30 I = NCONF+1,NVAR
               IF (ABS(DIAGL(I)) .GT. DLIM) THEN
                  BVECS(IBO+I) = G(I) / DIAGL(I)
               ELSE
                  BVECS(IBO+I) = G(I) / DLIM
               END IF
   30       CONTINUE
         END IF
         NOVEC = NOVEC + 1
         ISTBO = ISTBO + NWOPT
      END IF
      NROOTX = NROOTX - 1
C
      BINMEM = .TRUE.
      FIRST  = .TRUE.
      JSTATE = 0
   40 CONTINUE
      IF (NROOTX .GT. 0) THEN
         IF (FIRST) THEN
            REWIND LUIT1
            CALL MOLLB2(LABIT1,LABINF,LUIT1,LUERR)
            READ (LABINF(1),'(4X,I4)') JSTA
            IF (JSTA .LT. 0) THEN
               NWARN = NWARN + 1
               WRITE (LUW4,'(/A/A,I5/2A)')
     &         ' WARNING from NEOST; only gradient trial vector',
     &         '     number of trial vectors requested:',NROOTX+1,
     &         '     reason: label on LUITS was ',LABINF(2)
               NROOTX = 0
               GO TO 40
            END IF
            IF (LABINF(2) .NE. 'NEOSAVE ') OTRIAL = .FALSE.
            IF (OTRIAL) THEN
C           ... if orbital trial vectors, test if available
               READ (LUIT1,ERR=44,END=44) (DIAGL(I),I=1,NVAR)
               GO TO 46
   44          CONTINUE
                  WRITE (LUPRI,'(//A/)')
     *            ' >>>NEOST : no orbital trial vectors this iteration',
     *            '            only CI part available on LUIT1.'
                  OTRIAL = .FALSE.
   46          CONTINUE
               REWIND LUIT1
               CALL MOLLAB(LABIT1,LUIT1,lupri)
            END IF
         ELSE
            NCVEC = 0
            NOVEC = 0
            ISTBC = KBC
            ISTBO = KBO
         END IF
         MROOTX = MIN(NROOTX,NTSIM)
   60    CONTINUE
            JSTATE = JSTATE + 1
            IF (JSTATE.EQ.ISTATE) THEN
               READ (LUIT1)
            ELSE
               NROOTX = NROOTX - 1
               MROOTX = MROOTX -  1
               IF (.NOT.OTRIAL) THEN
                  CALL READT(LUIT1,NCONF,BVECS(ISTBC))
                  ISTBC = ISTBC + NCONF
                  NCVEC = NCVEC + 1
               ELSE
C                 860205-hjaaj -- prepared for using
C                 orbital start trial vectors from previous iter.
C                 ... implemented 870804-hjaaj.
                  CALL READT(LUIT1,NVAR,DIAGL)
                  IF (NCONF .GT. 1) THEN
                     CALL DCOPY(NCONF,DIAGL,1,BVECS(ISTBC),1)
                     ISTBC = ISTBC + NCONF
                     NCVEC = NCVEC + 1
                  END IF
                  IF (NWOPT .GT. 0) THEN
                     CALL DCOPY(NWOPT,DIAGL(1+NCONF),1,BVECS(ISTBO),1)
                     ISTBO = ISTBO + NWOPT
                     NOVEC = NOVEC + 1
                  END IF
               END IF
            END IF
         IF (MROOTX .GT. 0) GO TO 60
      END IF
C
C     Do microcanonical average on new BOVECS
C     (AVERAG returns immediately if no averageing requested)
C
   80 IF (NOVEC.GT.0) THEN
         CALL AVERAG(BVECS(KBO),NWOPT,NOVEC)
      END IF
C
C     Orthogonalize and renormalize start vectors.
C     NBPREV is updated in ORTBVC.
C
      KBPREV = MAX(ISTBC,ISTBO)
      KTT0   = KBPREV + NCOMAX
      NDMBC  = MAX(1,NCONF)
      NDMBO  = MAX(1,NWOPT)
      THRLDV = NVAR*THRLDP
      CALL ORTBVC(NCVEC,NOVEC,NDMBC,NDMBO,NBPREV,IBNDX,LUIT3,
     *            BVECS(KBC),BVECS(KBO),
     *            THRLDV,BVECS(KTT0),BVECS(KBPREV))
C     CALL ORTBVC(NBC,NBO,NDMBC,NDMBO,NBPREV,IBNDX,LUBVC,
C    *            BCVECS,BOVECS,THRLDP,TT0,BPREV)
C
      NCSIM = NCSIM + NCVEC
      NOSIM = NOSIM + NOVEC
      IF (NROOTX.GT.0) THEN
         BINMEM = .FALSE.
         FIRST  = .FALSE.
         GO TO 40
      END IF
C     930323-hjaaj: if no trial vectors and getbvo false, it must be
C     because this is really a one-configuration case; do getbvo
      IF (NCSIM + NOSIM .EQ. 0 .AND. .NOT. GETBVO) THEN
         GETBVO = .TRUE.
         NCVEC = 0
         NOVEC = 1
         ISTBC = KBC
         ISTBO = KBO + NWOPT
         CALL DCOPY(NWOPT,G(NCONF+1),1,BVECS(KBO),1)
         GO TO 80
      END IF
C
C
      IF (BINMEM) THEN
         KBONEW = KBC + NCSIM*NCONF
         IF (KBONEW.NE.KBO .AND. NOSIM.GT.0) THEN
            DO 120 I = 0,NOSIM*NWOPT-1
               BVECS(KBONEW+I) = BVECS(KBO+I)
  120       CONTINUE
         END IF
      ELSE
         NCSIM = -NCSIM
         NOSIM = -NOSIM
      END IF
C
      CALL QEXIT('NEOST ')
      RETURN
C
C *** end of subroutine NEOST
C
      END
C  /* Deck ortbvc */
      SUBROUTINE ORTBVC (NBC,NBO,NDMBC,NDMBO,NBPREV,IBNDX,LUBVC,
     *                   BCVECS,BOVECS,THRLDP,TT0,BPREV)
C
C Written 12-Apr-1984 by Hans Jorgen Aa. Jensen
C Revisions
C  11-Feb-1985 hjaaj (new version)
C  17-Nov-1988 pj (for abacus WITH SYMMETRY)
C   4-May-1994 hjaaj (Only eliminate input vectors if norm < THRLDP)
C Purpose:
C  Orthogonalize new b-vectors against all previous b-vectors
C  and among themselves, and renormalize.
C  (Orthogonalization is performed twice if round-off is large,
C   if larger than THRRND).
C
C Input:
C  BCVECS, non-orthogonal configurational b-vectors
C  BOVECS, non-orthogonal orbital b-vectors
C  NBC,    number of BCVECS
C  NBO,    number of BOVECS
C  THRLDP, threshold for linear dependence
C
C Output:
C  TT0(NBC+NBO), norm of b-vectors after orthogonalization but
C               before normalization (!).
C
C Scratch:
C  BPREV(nvar), scratch array for old b-vectors on LUBVC
C
#include <implicit.h>
      DIMENSION BCVECS(NDMBC,*),BOVECS(NDMBO,*)
      DIMENSION IBNDX(*),TT0(*),BPREV(*)
#include <ibndxdef.h>
      PARAMETER (THRRND = 1.D-2, THRTT = 1.D-4)
      PARAMETER (D1 = 1.0D0)
C
#include <priunit.h>
#include <infpri.h>
C
      CALL QENTER('ORTBVC')
C
      TLINDP = SQRT(THRLDP)
      NDMBC4 = MAX(4,NDMBC)
      NDMBO4 = MAX(4,NDMBO)
      IBC    = NBPREV
      IBO    = IBC + NBC
C
C     Initialize IBNDX for the new vectors
C
      DO 50 I = 1,NBC
         IBNDX(IBC+I) = JBCNDX
   50 CONTINUE
      DO 60 I = 1,NBO
         IBNDX(IBO+I) = JBONDX
   60 CONTINUE
C
C     Normalize input BCVECS and BOVECS
C
      DO 100 I = 1,NBC
         TT = DNORM2(NDMBC,BCVECS(1,I),1)
         IF (TT.LE.THRLDP) THEN
            TT0(I)       = TT
            IBNDX(IBC+I) = 0
            WRITE (LUPRI,3240) 'CSF',I,TT
         ELSE
            IF (TT.LT.THRTT) THEN
               CALL DSCAL(NDMBC,(D1/TT),BCVECS(1,I),1)
               TT = DNORM2(NDMBC,BCVECS(1,I),1)
            END IF
            CALL DSCAL(NDMBC,(D1/TT),BCVECS(1,I),1)
         END IF
  100 CONTINUE
C
      DO 200 I = 1,NBO
         TT = DNORM2(NDMBO,BOVECS(1,I),1)
         IF (TT.LE.THRLDP) THEN
            TT0(NBC+I)   = TT
            IBNDX(IBO+I) = 0
            WRITE (LUPRI,3240) 'orbital',I,TT
         ELSE
            IF (TT.LT.THRTT) THEN
               CALL DSCAL(NDMBO,(D1/TT),BOVECS(1,I),1)
               TT = DNORM2(NDMBO,BOVECS(1,I),1)
            END IF
            CALL DSCAL(NDMBO,(D1/TT),BOVECS(1,I),1)
         END IF
  200 CONTINUE
C
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
C
C     Orthogonalize new b-vectors agains previous b-vectors
C
      REWIND LUBVC
      DO 2000 K=1,NBPREV
         IF (IBNDX(K) .EQ. JBCNDX .AND. NBC .GT. 0) THEN
C        if (CI-type) then
            CALL READT(LUBVC,NDMBC,BPREV)
            DO 1910 J = 1,NBC
               IF (IBNDX(IBC+J).EQ.JBCNDX) THEN
                  TT = -DDOT(NDMBC,BPREV,1,BCVECS(1,J),1)
                  CALL DAXPY(NDMBC,TT,BPREV,1,BCVECS(1,J),1)
               END IF
 1910       CONTINUE
         ELSE IF (IBNDX(K) .EQ. JBONDX .AND. NBO .GT. 0) THEN
C        else if (orbital-type) then
            CALL READT(LUBVC,NDMBO,BPREV)
            DO 1920 J = 1,NBO
               IF (IBNDX(IBO+J).EQ.JBONDX) THEN
                  TT = -DDOT(NDMBO,BPREV,1,BOVECS(1,J),1)
                  CALL DAXPY(NDMBO,TT,BPREV,1,BOVECS(1,J),1)
               END IF
 1920       CONTINUE
         ELSE
            IF (IBNDX(K).NE.JBONDX .AND. IBNDX(K).NE.JBCNDX) THEN
               WRITE (LUERR,'(/1X,A,2I5)')
     *           'ORTBVC, illegal IBNDX value; K,IBNDX(K) =',K,IBNDX(K)
               CALL QUIT('SIRIUS ERROR: illegal IBNDX value in ORTBVC')
            END IF
            READ (LUBVC)
         END IF
 2000 CONTINUE
C
C     Normalize and orthogonalize new vectors against each other
C
      DO 2400 I = 1,NBC
      IF (IBNDX(IBC+I).EQ.0) GO TO 2400
C        ... orthogonalize using that prev. vectors are normalized
         DO 2300 J = 1,(I-1)
         IF (IBNDX(IBC+J).EQ.0) GO TO 2300
            TT = -DDOT(NDMBC,BCVECS(1,J),1,BCVECS(1,I),1)
            CALL DAXPY(NDMBC,TT,BCVECS(1,J),1,BCVECS(1,I),1)
 2300    CONTINUE
         TT = DNORM2(NDMBC,BCVECS(1,I),1)
         IF (ITURN.EQ.1) TT0(I)=TT
         IF (TT .LE. TLINDP) THEN
            WRITE (LUPRI,3250) 'CSF',I,TT
            IBNDX(IBC+I) = 0
         ELSE
            IF (TT .LT. THRRND) IROUND = IROUND+1
            IF (TT.LT.THRTT) THEN
               CALL DSCAL(NDMBC,(D1/TT),BCVECS(1,I),1)
               TT = DNORM2(NDMBC,BCVECS(1,I),1)
            END IF
            CALL DSCAL(NDMBC,(D1/TT),BCVECS(1,I),1)
         END IF
 2400 CONTINUE
C
      DO 2600 I = 1,NBO
      IF (IBNDX(IBO+I).EQ.0) GO TO 2600
C        ... orthogonalize using that prev. vectors are normalized
         DO 2500 J = 1,(I-1)
         IF (IBNDX(IBO+J).EQ.0) GO TO 2500
            TT = -DDOT(NDMBO,BOVECS(1,J),1,BOVECS(1,I),1)
            CALL DAXPY(NDMBO,TT,BOVECS(1,J),1,BOVECS(1,I),1)
 2500    CONTINUE
         TT = DNORM2(NDMBO,BOVECS(1,I),1)
         IF (ITURN.EQ.1) TT0(NBC+I) = TT
         IF (TT .LE. TLINDP) THEN
            WRITE (LUPRI,3250) 'orbital',I,TT
            IBNDX(IBO+I) = 0
         ELSE
            IF (TT .LT. THRRND) IROUND = IROUND+1
            IF (TT.LT.THRTT) THEN
               CALL DSCAL(NDMBO,(D1/TT),BOVECS(1,I),1)
               TT = DNORM2(NDMBO,BOVECS(1,I),1)
            END IF
            CALL DSCAL(NDMBO,(D1/TT),BOVECS(1,I),1)
         END IF
 2600 CONTINUE
      IF (IROUND.GT.0 .AND. ITURN.EQ.1) GO TO 1500
C
 3240 FORMAT(/' (ORTBVC) new ',A,' b-vector no.',I3,
     *       /' is removed at input because of too small norm;'
     *       /' norm of input vector is',1P,D12.5)
 3250 FORMAT(/' (ORTBVC) new ',A,' b-vector no.',I3,
     *       /' is removed because of linear dependence;'
     *       /' norm of vector after Gram-Schmidt''s orthogonalization',
     *        1P,D12.5)
C
C Save new vector together with old ones on LUBVC
C
      NLAST = NBPREV + NBC + NBO
      J = 0
      DO 4300 I = 1,NBC
         IF (IBNDX(IBC+I) .NE. 0) THEN
            J = J + 1
            IF (J .LT. I) THEN
               CALL DCOPY(NDMBC,BCVECS(1,I),1,BCVECS(1,J),1)
            END IF
            CALL WRITT(LUBVC,NDMBC4,BCVECS(1,J))
         END IF
 4300 CONTINUE
      IF (J .LT. NBC) THEN
         WRITE (LUPRI,4310) NBC,J
         NBC = J
      END IF
C
      J = 0
      DO 4400 I = 1,NBO
         IF (IBNDX(IBO+I) .NE. 0) THEN
            J = J + 1
            IF (J .LT. I) THEN
               CALL DCOPY (NDMBO,BOVECS(1,I),1,BOVECS(1,J),1)
            END IF
            CALL WRITT(LUBVC,NDMBO4,BOVECS(1,J))
         END IF
 4400 CONTINUE
      IF (J .LT. NBO) THEN
         WRITE (LUPRI,4410) NBO,J
         NBO = J
      END IF
C
C     update NBPREV
C
      DO 5100 I = NBPREV+1,NBPREV+NBC
         IBNDX(I) = JBCNDX
 5100 CONTINUE
      NBPREV = NBPREV + NBC
      DO 5200 I = NBPREV+1,NBPREV+NBO
         IBNDX(I) = JBONDX
 5200 CONTINUE
      NBPREV = NBPREV + NBO
      DO 5300 I = NBPREV+1,NLAST
         IBNDX(I) = 0
 5300 CONTINUE
C
 4310 FORMAT(/' (ORTBVC) NBCSIM reduced from',I3,' to',I3)
 4410 FORMAT(/' (ORTBVC) NBOSIM reduced from',I3,' to',I3)
C
C *** End of subroutine ORTBVC
C
      CALL QEXIT('ORTBVC')
      RETURN
      END
C  /* Deck updbet */
      SUBROUTINE UPDBET (REDL,EVEC,EVAL,WRK,LFREE)
C
C Written 31-Mar-1984 by Hans Joergen Aa. Jensen
C
C Revisions 26-Oct-1984 / 19-May-1984 hjaaj
C  20-Dec-1994 hjaaj (check if EVEC(1,ISTATE) .eq. 0)
C
C  This is the routine for update of beta in the reduced matrix
C  to give the desired step length (trust region control).
C
C MOTECC-90: The purpose of this module, UPDBET, and the algorithms
C            used are described in Chapter 8 Section C.3 of MOTECC-90
C            "The Dynamical Update of the Damping Factor"
C
C
C
#include <implicit.h>
C
      DIMENSION REDL(*),EVEC(NREDL,*),EVAL(*)
      PARAMETER (ITBMAX = 100)
      DIMENSION WRK(ITBMAX,4)
C
      PARAMETER (BETINC = 0.2D0, SECFAC = 0.1D0)
      PARAMETER (DHALF=0.5D0, D0=0.0D0, D1=1.0D0, D2=2.0D0)
      PARAMETER (DP1 = 0.1D0, DP9 = 0.9D0)
#include <priunit.h>
#include <dummy.h>
#include <thrzer.h>
C
C Used from common blocks:
C   INFINP : ISTATE
C   INFVAR : NWOPT
C   INFOPT : BETA,BETMIN,(BETMAX,)RTRUST,RTTOL,RBETA,GRDNRM,NREDL
C   INFPRI : P6FLAG(*)
C
#include <maxorb.h>
#include <infinp.h>
#include <infvar.h>
#include <infopt.h>
#include <infpri.h>
C
C local:
C
      LOGICAL SECANT
C
      CALL QENTER('UPDBET')
C
C     space for NEORED
C
      LWRK1 = LFREE - 3*ITBMAX
C
C Initialize step scaling to no scaling:
C
      GAMMA = D1
      BETL  = ABS(BETA)
C
C --- If CI always select BETA = 1:
C
      IF (NWOPT.EQ.0) THEN
         BETL = D1
         IF (BETL.NE.BETA) THEN
            BETA = D1
            CALL NEORED (2,0,0,WRK(1,4),DUMMY,BETL,DUMMY,
     *                   IDUMMY,NREDL,REDL,EVAL,EVEC,LWRK1)
C           CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,G,
C    *                   IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
         END IF
         GO TO 9999
      END IF
C
C --- This was not a CI, find optimal BETA
C     (branch for ground state or excited state section)
C
      ITB = 0
      IF (ISTATE.GT.1) THEN
         BETL = BETMIN - BETINC
         GO TO 400
      END IF
C
C ***
C *** Ground state section (ISTATE.EQ.1):
C ***
C
      SECANT = .FALSE.
      RTTOL1 = RTTOL
      RTTOL2 = D1 / (DP9 + DP1*RTTOL)
      FBI2   = D0
  100 CONTINUE
      ITB = ITB + 1
      IF (BETL .LT. BETMIN) THEN
         BETL = BETMIN
         CALL NEORED (2,0,0,WRK(1,4),DUMMY,BETL,DUMMY,
     *                IDUMMY,NREDL,REDL,EVAL,EVEC,LWRK1)
C        CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,G,
C    *                IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
      END IF
      IF (ABS(EVEC(1,ISTATE)) .LT. THRZER) THEN
         IF (ISTATE.LT.NREDL) THEN
            NWARN = NWARN + 1
            WRITE(LUPRI,'(/A)') 'UPDBET WARNING: '//
     &      'Symmetry break, will try to switch ISTATE and ISTATE+1 '

            CALL DSWAP(NREDL,EVEC(1,ISTATE),1,EVEC(1,ISTATE+1),1)
            EVALX = EVAL(ISTATE)
            EVAL(ISTATE) = EVAL(ISTATE+1)
            EVAL(ISTATE+1) = EVALX
         END IF
         IF (ABS(EVEC(1,ISTATE)) .LT. THRZER) GO TO 910
      END IF
      RBETL = SQRT(EVEC(1,ISTATE) ** (-2) - D1) / BETL
      WRK(ITB,1) = BETL
      WRK(ITB,2) = RBETL
      WRK(ITB,3) = EVAL(ISTATE)
      FBETL = RBETL - RTRUST
C
      IF (RBETL.GT.RTTOL1*RTRUST) THEN
C     (step too long, increase beta)
         IF (SECANT) GO TO 200
         IF (FBI2*FBETL .LT. D0) THEN
            SECANT = .TRUE.
            GO TO 200
         END IF
         FBI2 = FBETL
         BI2  = BETL
         BETL = BETL * MIN(D2,(RBETL/RTRUST))
         IF (ITB .EQ. 1) BETL = MAX(BETL,D1)
         GO TO 300
      ELSE IF (RBETL.LT.RTTOL2*RTRUST) THEN
C          if test for RBETL.LT.RTRUST then algorithm may not converge
C          as RBETL may approach RTRUST from below.
C     (step shorter then RTRUST, if NR-like
C      -- i.e. beta = betmin -- exit loop)
#if defined (VAR_OLDCODE)
CHJ-860925: not needed any more after change in residual calc.
CHJ-S-850211: revised to set BETA = D1 when local.
         IF (BETL .EQ. BETMIN) THEN
            IF (RBETL .LT. DHALF*RTRUST) BETL = MAX(D1,BETMIN)
            IF (BETL .EQ. BETMIN) THEN
               GO TO 800
            ELSE
               GO TO 700
            END IF
         END IF
CHJ-E-850211
#else
         IF (BETL .EQ. BETMIN) GO TO 800
#endif
         IF (SECANT) GO TO 200
         IF (FBI2*FBETL .LT. D0) THEN
            SECANT = .TRUE.
            GO TO 200
         END IF
         IF (BETL .LE. D1) THEN
            BETL = BETMIN
         ELSE
            BETL = BETL * (RBETL/RTRUST)
         END IF
         GO TO 300
      ELSE
C     (step within tolerance on RTRUST, exit loop)
         GO TO 800
      END IF
C
C Use the modified "Illinois method", a variation of the secant method,
C  to obtain next guess for BETA (Illinois: SECFAC = one half).
C  (the "Pegasus method", SECFAC = FBI2 / (FBI2+FBETL),
C   did not work as well)
C
  200 CONTINUE
      IF (FBI2*FBETL .LT. D0) THEN
         FBI1 = FBI2
         BI1  = BI2
      ELSE
         FBI1 = SECFAC*FBI1
      END IF
      FBI2 = FBETL
      BI2  = BETL
      BETL = (FBI2*BI1 - FBI1*BI2) / (FBI2 - FBI1)
C
C Solve for eigenvalues and -vectors with new beta
C
  300 IF (ITB.EQ.ITBMAX) GO TO 900
      CALL NEORED (2,0,0,WRK(1,4),DUMMY,BETL,DUMMY,
     *             IDUMMY,NREDL,REDL,EVAL,EVEC,LWRK1)
C     CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,G,
C    *             IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
      GO TO 100
C
C ***
C *** Excited state section (ISTATE.GT.1):
C
  400 CONTINUE
         ITB = ITB + 1
         BETL = BETL + BETINC
         CALL NEORED (2,0,0,WRK(1,4),DUMMY,BETL,DUMMY,
     *                IDUMMY,NREDL,REDL,EVAL,EVEC,LWRK1)
C        CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,G,
C    *                IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
         IF (ABS(EVEC(1,ISTATE)) .LT. THRZER) THEN
            IF (ISTATE.LT.NREDL) THEN
               NWARN = NWARN + 1
               WRITE(LUPRI,'(/A)') 'UPDBET WARNING: '//
     &         'Symmetry break, will try to switch ISTATE and ISTATE+1'
               CALL DSWAP(NREDL,EVEC(1,ISTATE),1,EVEC(1,ISTATE+1),1)
               EVALX = EVAL(ISTATE)
               EVAL(ISTATE) = EVAL(ISTATE+1)
               EVAL(ISTATE+1) = EVALX
            END IF
            IF (ABS(EVEC(1,ISTATE)) .LT. THRZER) GO TO 910
         END IF
         RBETL = SQRT(EVEC(1,ISTATE) ** (-2) - D1) / BETL
         WRK(ITB,1) = BETL
         WRK(ITB,2) = RBETL
         WRK(ITB,3) = EVAL(ISTATE)
         IF (RBETL .LT. RTTOL*RTRUST) THEN
            J = ITB
            GAMMA = D1
            GO TO 800
         END IF
      IF (ITB.EQ.1) GO TO 400
         IF (RBETL .GT. WRK(ITB-1,2)) THEN
            J = MAX(ITB-2,1)
            BETL = WRK(J,1)
            GAMMA = RTRUST / WRK(J,2)
            IF (ITB.LT.ITBMAX) GO TO 700
         END IF
      IF (ITB.LT.ITBMAX) GO TO 400
C
      GAMMA = RTRUST / RBETL
      GO TO 800
C
C ***
C *** normal exit
C
  700 CONTINUE
      ITB = ITB + 1
      CALL NEORED (2,0,0,WRK(1,4),DUMMY,BETL,DUMMY,
     *             0,NREDL,REDL,EVAL,EVEC,LWRK1)
C     CALL NEORED (ICTL,NCSIM,NOSIM,BVECS,SVECS,BETVAL,G,
C    *             IBNDX,NREDL,REDL,EVAL,EVEC,LBVECS)
      IF (ABS(EVEC(1,ISTATE)) .LT. THRZER) THEN
         IF (ISTATE.LT.NREDL) THEN
            NWARN = NWARN + 1
            WRITE(LUPRI,'(/A)') 'UPDBET WARNING: '//
     &      'Symmetry break, will try to switch ISTATE and ISTATE+1 '
            CALL DSWAP(NREDL,EVEC(1,ISTATE),1,EVEC(1,ISTATE+1),1)
            EVALX = EVAL(ISTATE)
            EVAL(ISTATE) = EVAL(ISTATE+1)
            EVAL(ISTATE+1) = EVALX
         END IF
         IF (ABS(EVEC(1,ISTATE)) .LT. THRZER) GO TO 910
      END IF
      RBETL  = SQRT(EVEC(1,ISTATE) ** (-2) - D1) / BETL
      WRK(ITB,1) = BETL
      WRK(ITB,2) = RBETL
      WRK(ITB,3) = EVAL(ISTATE)
  800 CONTINUE
      IF (P6FLAG(8) .AND. BETA.NE.BETL) THEN
         WRITE (LUPRI,8010) BETL,ITB,NREDL,GAMMA
         CALL OUTPUT (WRK,1,ITB,1,3,ITBMAX,3,1,LUPRI)
      END IF
 8010 FORMAT(/' (UPDBET) New beta =',F10.5,' found in',I3,' iterations.'
     *       /' NREDL =',I3,', step scaling factor, gamma, is:',F10.5,
     *       /' Reduced L values of beta, step length, and ',
     *        'level shift:')
      BETA  = BETL
      RBETA = RBETL
 9999 CALL QEXIT('UPDBET')
      RETURN
C
C ***
C *** ERROR exits
C
  900 CONTINUE
      WRITE (LUPRI,9010) ITBMAX
      CALL OUTPUT (WRK,1,ITBMAX,1,3,ITBMAX,3,1,LUPRI)
      CALL QTRACE(LUPRI)
      CALL QUIT('ERROR: could not find new beta (UPDBET)')
 9010 FORMAT(/' (UPDBET) ERROR, new beta not found in',I3,
     *        ' iterations.',
     *       /' Reduced L values of beta, step length, and ',
     *        'level shift:')
  910 CONTINUE
      WRITE (LUPRI,9110) ISTATE
      CALL OUTPAK(REDL,NREDL,1,LUPRI)
      WRITE (LUPRI,9112) BETL,(EVAL(I),I=1,NREDL)
      LIMP = MIN(((ISTATE+5)/4) * 4, NREDL)
      WRITE (LUPRI,9118)
      CALL OUTPUT(EVEC,1,NREDL,1,LIMP,NREDL,NREDL,1,LUPRI)
      CALL QTRACE(LUPRI)
      CALL QUIT(
     &   'ERROR (UPDBET): no overlap between NEO eig.vec. and CREF')
 9110 FORMAT(/' (UPDBET) ERROR, no overlap between NEO eigenvector',
     *  /' no. ISTATE =',I3,' and the reference vector.',
     *  /' Probable cause: NEO eigenvector no. ISTATE breaks symmetry.',
     *  /' - The reduced L matrix :')
 9112 FORMAT(/' - Beta and eigenvalues :'
     *       /F10.2,1P,4D15.6,/,(10X,1P,4D15.6))
 9118 FORMAT(/' - and the first few eigenvectors :')
C
C -- end of subroutine UPDBET
C
      END
C  /* Deck neored_chk */
      SUBROUTINE NEORED_CHK(IBNDX,WRK,LFREE)
C
C 8-Nov-1986 Hans Joergen Aa. Jensen
C
C Purpose:
C   Calculate full reduced Hessian for symmetry check;
C   asymmetry would indicate probable bug in LINTRN module.
C
#include <implicit.h>
      DIMENSION IBNDX(*), WRK(LFREE)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DIFTST = 1.0D-8)
#include <ibndxdef.h>
C
C  Used from common blocks:
C   INFVAR : NCONF,NWOPT,NVAR
C   INFOPT : NREDL
C   INFTAP : LUIT3, LUIT5
C
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <infopt.h>
#include <inftap.h>
C
      LOGICAL CITYPE
C
      REDH(I,J) = WRK((J-1)*NREDH+I)
      REDS(I,J) = WRK((KREDS-1)+(J-1)*NREDL+I)
C
      NREDH  = NREDL - 1
      NCOMAX = MAX(NCONF,NWOPT)
      KREDH  = 1
      KREDS  = KREDH  + NREDH*NREDH
      KBVEC  = KREDS  + NREDL*NREDL
      KSVECS = KBVEC  + NCOMAX
      LSVECS = LFREE  - KSVECS
      IF (LSVECS .LT. NREDH*NVAR) THEN
C        insufficient space for the algorithm used
         WRITE (LUPRI,'(//2A,A/)') ' %%% Insufficient space',
     *      ' for symmetry check of reduced Hessian.',
     *      '     Nothing done.'
         GO TO 9999
      END IF
C
      JSVECI = KSVECS
      REWIND LUIT5
      IF (NCONF .GT. 1) READ (LUIT5)
      IF (NWOPT .GT. 0) READ (LUIT5)
      DO 100 I = 1,NREDH
         IF (NCONF .GT. 1) THEN
            CALL READT(LUIT5,NCONF,WRK(JSVECI))
         ELSE
            WRK(JSVECI) = D0
         END IF
         JSVECI = JSVECI + NCONF
         IF (NWOPT .GT. 0) CALL READT(LUIT5,NWOPT,WRK(JSVECI))
         JSVECI = JSVECI + NWOPT
  100 CONTINUE
C
C     REDH(I,*) = B(*,I)t S(*,*)
C
      REWIND LUIT3
      READ (LUIT3)
      DO 200 I = 1,NREDH
         CITYPE = (IBNDX(1+I) .EQ. JBCNDX)
         IF (CITYPE) THEN
            CALL READT(LUIT3,NCONF,WRK(KBVEC))
            CALL DGEMM('T','N',1,NREDH,NCONF,1.D0,
     &                 WRK(KBVEC),NCONF,
     &                 WRK(KSVECS),NVAR,0.D0,
     &                 WRK(I),NREDH)
C           note: WRK(I) is REDH(I,1)
         ELSE
            CALL READT(LUIT3,NWOPT,WRK(KBVEC))
            CALL DGEMM('T','N',1,NREDH,NWOPT,1.D0,
     &                 WRK(KBVEC),NWOPT,
     &                 WRK(KSVECS+NCONF),NVAR,0.D0,
     &                 WRK(I),NREDH)
         END IF
  200 CONTINUE
C
C     Reduced overlap matrix, get b vectors
C
      REWIND LUIT3
      JBVEC = KBVEC
      DO 220 I = 1,NREDL
         CITYPE = (IBNDX(I) .EQ. JBCNDX)
         IF (CITYPE) THEN
            CALL READT(LUIT3,NCONF,WRK(JBVEC))
            JBVEC = JBVEC + NCONF
         ELSE
            CALL READT(LUIT3,NWOPT,WRK(JBVEC))
            JBVEC = JBVEC + NWOPT
         END IF
  220 CONTINUE
C
C     REDS(I,J) = B(*,I)t B(*,J)
C
      CALL DZERO(WRK(KREDS),(NREDL*NREDL))
      JSJI   = KREDS - 1
      JBVECI = KBVEC
      DO 280 I = 1,NREDL
         CITYPE = (IBNDX(I) .EQ. JBCNDX)
         IF (CITYPE) THEN
            LBVEC = NCONF
         ELSE
            LBVEC = NWOPT
         END IF
         JBVECJ = KBVEC
         DO 260 J = 1,NREDL
            JSJI = JSJI + 1
            IF (IBNDX(J) .EQ. IBNDX(I)) THEN
               WRK(JSJI) = DDOT(LBVEC,WRK(JBVECJ),1,WRK(JBVECI),1)
            END IF
            IF (IBNDX(J) .EQ. JBCNDX) THEN
               JBVECJ = JBVECJ + NCONF
            ELSE
               JBVECJ = JBVECJ + NWOPT
            END IF
  260    CONTINUE
         IF (CITYPE) THEN
            JBVECI = JBVECI + NCONF
         ELSE
            JBVECI = JBVECI + NWOPT
         END IF
  280 CONTINUE
C
C     check reduced matrices
C
      WRITE (LUPRI,'(//A/A/A,1P,D10.2,A)')
     *    ' %%% Non-symmetric elements in redH',
     *    '     and wrong elements in redS',
     *    '     ( threshold is',DIFTST,' ) :'
      NERROR = 0
      DO 340 I = 2,NREDH
         DO 320 J = 1,(I-1)
            DIFF = ABS( REDH(I,J) - REDH(J,I) )
            IF (DIFF .GT. DIFTST) THEN
               NERROR = NERROR + 1
               WRITE (LUPRI,'(A,2I5,2F20.14)')
     *            ' i, j, redh(i,j), redh(j,i) :',
     *            I,J,REDH(I,J),REDH(J,I)
            END IF
            IF (ABS(REDS(I+1,J+1)) .GT. DIFTST .OR.
     *          ABS(REDS(J+1,I+1)) .GT. DIFTST) THEN
               NERROR = NERROR + 1
               WRITE (LUPRI,'(A,2I5,2F20.14)')
     *            ' i, j, reds(i,j), reds(j,i) :',
     *            I,J,REDS(I+1,J+1),REDS(J+1,I+1)
            END IF
  320    CONTINUE
         IF (ABS(REDS(I+1,  1)) .GT. DIFTST .OR.
     *       ABS(REDS(  1,I+1)) .GT. DIFTST) THEN
            NERROR = NERROR + 1
            WRITE (LUPRI,'(A,I10,2F20.14)')
     *         ' i, <bi , cref>, <cref, bi> :',
     *         I,REDS(I+1,1),REDS(1,I+1)
         END IF
         IF (ABS(REDS(I+1,I+1) - D1) .GT. DIFTST) THEN
            NERROR = NERROR + 1
            WRITE (LUPRI,'(A,2I5,F20.14)')
     *         ' i, i, reds(i,i)            :',
     *         I,I,REDS(I+1,I+1)
         END IF
  340 CONTINUE
      IF (NERROR .EQ. 0) WRITE (LUPRI,'(A)') ' None, redH is symmetric.'
C
      WRITE (LUPRI,'(/A,(/15I5))') ' IBNDX vector:',
     &   (IBNDX(I), I = 1,NREDL)
      WRITE (LUPRI,'(//A)') ' Reduced Hessian (as bT*s) :'
      CALL OUTPUT(WRK(KREDH),1,NREDH,1,NREDH,NREDH,NREDH,1,LUPRI)
      WRITE (LUPRI,'(//A)') ' Reduced overlap (as bT*b) :'
      CALL OUTPUT(WRK(KREDS),1,NREDL,1,NREDL,NREDL,NREDL,1,LUPRI)
C
      IF (NERROR .GT. 0) CALL QUIT('ERROR detected in NEORED_CHK.')
 9999 RETURN
      END
