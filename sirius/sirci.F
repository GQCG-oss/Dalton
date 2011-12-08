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
C  /* Deck cimain */
      SUBROUTINE CIMAIN(WORK,LWORK,ICONV)
C
C 28-Jun-1985 hjaaj
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION WORK(LWORK)
#include <dummy.h>
C
      PARAMETER (MXCIDFT = 20)
C
C Used from common blocks:
C   INFINP: ICI0,NROOCI,MXCIMA,THRCI,LSYM
C   INFORB: NCMOT,NASHT
C   INFDIM: LCINDX
C   dftcom: DFTADD
C
#include <maxorb.h>
#include <gnrinf.h>
#include <infinp.h>
#include <inforb.h>
#include <infdim.h>
#include <dftcom.h>
C
      LOGICAL ADDSRI_SAVE
C
      CALL QENTER('CIMAIN')
      ICIST  = ICI0
      NCROOT = NROOCI
      THRCIX = THRCI
C
C
      KFREE  = 1
      LFREE  = LWORK
      CALL MEMGET('REAL',KCINDX,LCINDX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KECI  ,NCROOT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KICROO,NCROOT,WORK,KFREE,LFREE)
C
      ADDSRI_SAVE = ADDSRI
      IF (DOCISRDFT) THEN
#ifndef MOD_SRDFT
         call quit('CI-srdft not implemented in this version')
#else
C        ... do not add SR integrals to Fock matrices for CI-DFT
C            (will be true on entry if CI-srDFT preceded by normal RHF)
         ADDSRI = .FALSE.
C        ... Dont risk adding DFT contributions in later call of FCKMAT
C            The SR-DFT terms needed in CI-DFT is handled by SRFMAT.
C            Maybe this should be added to the generel code and not 
C            just the CIDFT code if one wants to run a CI with KS orbitals !
C            (if the CISRDFT was specified together with a presceding
C             HFSRDFT call, we assume the HFSRDFT is done by now!)
         DFTADD    = .FALSE.
         DOHFSRDFT = .FALSE.
#endif
      ENDIF
C
      IF (NASHT .GT. 1 .AND. .NOT. HSROHF)
     *   CALL GETCIX(WORK(KCINDX),LSYM,LSYM,WORK(KFREE),LFREE,0)
      CALL READMO(WORK(KCMO),9,WORK(KFREE),LFREE)
C
C     =================
      IF (.NOT.DOCISRDFT .OR. ISTACI .LE. 0) THEN
C     =================
         CALL CICTL(ICIST,NCROOT,MXCIMA,THRCIX,WORK(KCMO),WORK(KCINDX),
     *              WORK(KECI),ICONV,WORK(KICROO),WORK(KFREE),LFREE)
C     =================
      ELSE
C     =================
C
C     CIDFT Macro-iteration control
C
#ifndef MOD_SRDFT
         call quit('srdft not implemented in this version')
#else
         CALL MEMGET('REAL',KECID,NCROOT,WORK,KFREE,LFREE)
         IMCITR = 0
         DO I = 1,MXCIDFT
            IMCITR = IMCITR + 1
            CALL CICTL(ICIST,NCROOT,MXCIMA,THRCIX,WORK(KCMO),
     *                 WORK(KCINDX),WORK(KECI),ICONV,WORK(KICROO),
     *                 WORK(KFREE),LFREE)
            ICIST  = 5
C           ... restart CI in next macro CIDFT iteration ...
            ECIDFT = WORK(KECI + (ISTACI-1) )
            WORK(KECID + I - 1) = ECIDFT
            IF(IMCITR.NE.1) EDIFF = ECIDFT - WORK(KECID + I - 2) 
            IF (IMCITR .GE. MXCIDFT) THEN
               WRITE(LUPRI,'(///A//A)') 
     &  '  >>>>> CI-srDFT optimization -- Summary <<<<<',
     &  '  ** Maximum number of Macro CI-srDFT iterations reached **'
               GO TO 5
            ELSE IF ((ABS(EDIFF).LT.THRCIDFT).AND.(IMCITR.NE.1)) THEN
               WRITE(LUPRI,'(///A//A,I3,A)') 
     &         '  >>>>> CI-DFT optimization -- Summary <<<<<',
     &         '  ** Convergence reached in',IMCITR,
     &         ' CI-DFT macro iterations.'
               GO TO 5
            ENDIF
         ENDDO
5        CONTINUE
            WRITE(LUPRI,'(3(/2X,A))') 
     &      '---------------------------------------------------------',
     &      '       Itr.           Energy                Delta(E)',
     &      '---------------------------------------------------------'
            WRITE(LUPRI,'(9X,I2,3X,F20.12,A)')
     &      1,WORK(KECID),'             ------'
            DO I = 2,IMCITR
               DE = WORK(KECID+I-2) - WORK(KECID+I-1)
               WRITE(LUPRI,'(9X,I2,2(3X,F20.12))')
     &         I,WORK(KECID + I - 1),DE
            ENDDO
            WRITE(LUPRI,'(2X,A/)') 
     &      '---------------------------------------------------------'
            WRITE(LUPRI,'(/2X,A,F20.12/)')
     &      'Final CI-DFT energy  : ',WORK(KECID + IMCITR -1)
#endif
C
C     =================
      ENDIF
C     =================
      CALL MEMREL('CIMAIN',WORK,1,1,KFREE,LFREE)
      ADDSRI = ADDSRI_SAVE
C
      CALL QEXIT('CIMAIN')
      RETURN
      END
C  /* Deck cictl */
      SUBROUTINE CICTL (ICICTL,NCROOT,MAXITC,THRCIX,CMO,INDXCI,
     *                  ECI,ICONV,ICROOT,WRK,LWRK)
C
C 20-Jun-1985 hjaaj; l.r. 18-Sep-1990 hjaaj
C
C Purpose:
C  Perform CI optimization.
C  If MAXITC .le. 0 then return after CIST
C  (e.g. when called from OPTST to find lowest diagonal elements).
C
C Input:
C  ICICTL     control of methods used
C  NCROOT     number of CI states to converge (from below)
C  MAXITC     maximum number of iterations
C  THRCIX     threshold for convergence
C  CMO(*)     m.o. coefficients
C
C Output:
C  ECI(*)     final CI energies.
C  ICONV      number of CI states converged
C             (NCROOT-ICONV states are not converged)
C  ICROOT(*)  index of CI states not converged
C
#include <implicit.h>
      DIMENSION CMO(*), INDXCI(*), ECI(*), ICROOT(*), WRK(*)
C
#include <iratdef.h>
#include <thrldp.h>
      PARAMETER ( D0 = 0.0D0 , DP5 = 0.5D0, D1 = 1.0D0 )
#include <dummy.h>
      LOGICAL CINMEM
C
C Used from common blocks:
C   CICB01 : NCIRED,NJCR,JCROOT(*),THRCIT
C   INFINP : POTNUC,FLAG(*),LSYM
C   INFVAR : NCONF
C   INFORB : NNASHX
C   INFDIM : LCINDX,LACIMX,LBCIMX,LPHPMX,MAXPHP
C   INFTAP : LUIT1,LUIT3,LUIT5
C
#include <maxorb.h>
#include <priunit.h>
#include <cicb01.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C
C     Local variables
C
      LOGICAL DISKH2, NATORB, MOISNO
C
C
C     *****************************************************************
C     Analyze input control specification
C     *****************************************************************
C
C
      CALL QENTER('CICTL ')
      TIMTOT = SECOND()
      WRITE (LUW4,'(//A///)')
     *   ' <<<<< Output from SIRIUS CI module (CICTL) >>>>>'
      IF (LUPRI .NE. LUW4) WRITE (LUPRI,'(//A///)')
     *   ' <<<<< Output from SIRIUS CI module (CICTL) >>>>>'
      IF (NASHT .LE. 1) THEN
         WRITE(LUPRI,'(//A,I8/)')
     *      ' *** ERROR *** CICTL CALLED, BUT NASHT =',NASHT
         CALL QTRACE(LUPRI)
         CALL QUIT('*** ERROR (CICTL) NASHT .LE. 1')
      END IF
      IF (NCROOT .LE. 0) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//A,I4,A/)')
     *      ' *** WARNING,',NCROOT,' CI roots requested, nothing done.'
         ICONV = 0
         CALL QTRACE(LUPRI)
         GO TO 9999
      END IF
      IF (NCROOT .GT. NCONF) THEN
         WRITE (LUPRI,'(//A,I5/A,I5)')
     *      ' *** number of CI roots to calculate reduced from',NCROOT,
     *      '     to the number of configurations =',NCONF
         NCROOT = NCONF
      END IF
      IF (NCROOT .GT. MXCIRT) THEN
         WRITE (LUPRI,'(//A,I5,A/A,I5,A)')
     *      ' *** ERROR in sirci,',NCROOT,' CI roots requested',
     *      '     which exceeds maximum of',MXCIRT,', increase MXCIRT.'
         CALL QUIT('SIRCI, number of CI roots exceeds predefined '//
     &             'maximum')
      END IF
C
C
C     *****************************************************************
C     Core allocation
C     *****************************************************************
C
C
      KFCAC  = 1
      IF (DOCISRDFT) THEN
C        allocate space for Short Range ACtive effective 1-el. integrals
         KSRAC  = KFCAC  + NNASHX
         KH2AC  = KSRAC  + NNASHX
      ELSE
         KSRAC  = KFCAC
         KH2AC  = KFCAC  + NNASHX
      END IF
C     ... first attempt to keep H2AC in core :
      DISKH2 = .FALSE.
      IF (FLAG(69)) DISKH2 = .TRUE.
   10 CONTINUE
C
      IF (DISKH2) THEN
         KCDIAG = KH2AC  + NNASHX
C        ... allocate space for one distribution (* * / u v) of H2AC
      ELSE
         KCDIAG = KH2AC  + NNASHX*NNASHX
C        ... allocate space for whole H2AC
      END IF
      CALL PHPINI(LPHPMX,NCONF,0,MAXPHP)
C     CALL PHPINI(LPHPT,NCVAR,NOVAR,MAXPHP)
      KW1    = KCDIAG + LPHPMX
      LW1    = LWRK   - KW1
C
C
      KREDCI = KW1
      KEVAL  = KREDCI + (MAXRC*MAXRC+MAXRC) /2
      KEVEC  = KEVAL  + MAXRC
      KRESID = KEVEC  + MAXRC * MAXRC
      KW2    = KRESID + MAXRC
      LW2    = LWRK   - KW2
C
C
C     work space needed for CILIN is .lt. LA + NCSIM*LB
C
      LA = LACIMX
      LB = 2*NCONF + LBCIMX
C
      MAXSIM = MIN(NCROOT, (LW2 - LA) / LB )
      IF (MAXSIM .LE. 0) THEN
         IF (.NOT. DISKH2) THEN
C           ... try put H2AC on disk :
            DISKH2 = .TRUE.
            GO TO 10
         END IF
         LNEED = KW2 + LA + LB
         CALL ERRWRK('CICTL.CILIN',-LNEED,LWRK )
      END IF
      IF (IPRSTAT .GT. 1) WRITE (LUSTAT,'(/A,I8/A,4I12/)')
     &   ' CICTL :  MAXSIM =',MAXSIM,
     &   ' KW1, KW2, LWRK  =',KW1,KW2,LWRK
C
   11 KCVECS = KW2
      KSVECS = KCVECS + MAXSIM*NCONF
      LSVECS = LWRK   - KSVECS
      KW3    = KSVECS + MAXSIM*NCONF
      LW3    = LWRK   - KW3
      MW3    = MAX(KW3 + 2 * MAXRC, KSVECS + MAXRC*MAXRC)
      IF (MW3 .GT. LWRK ) THEN
         MAXSIM = MAXSIM - 1 - (MW3 - LWRK ) / (2*NCONF)
         WRITE (LUSTAT,'(/A,I8/)')
     &      ' CICTL: MAXSIM REDUCED TO',MAXSIM
         IF (MAXSIM .GT. 0) GO TO 11
         CALL ERRWRK('CICTL.CIRED',MW3,LWRK )
      END IF
C
      NCONF4 = MAX(4,NCONF)
C
C
C     *****************************************************************
C     Initialization
C     Set up initial guess for CI vectors
C     *****************************************************************
C
C
      NCIRED = 0
      THRCIT = THRCIX
      NJCR   = NCROOT
      JCONV  = NCROOT
      TIMITR = SECOND()
      IF (.NOT. FLAG(14)) THEN
      IF (DONEVPT) THEN
         CALL TRACTL(4,CMO,WRK,LWRK)
      ELSE
         CALL TRACTL(0,CMO,WRK,LWRK)
      ENDIF
C        CALL TRACTL(ITRLVL,CMO,WRK,LWRK)
         FLAG(14) = .TRUE.
      END IF
      TIMITR = SECOND() - TIMITR
C
C     *** Obtain FCAC and H2AC integral matrices
C
      CALL CIINTS(DISKH2,CMO,
     *            EMY,WRK(KFCAC),WRK(KH2AC),WRK(KW1),LW1)
C     CALL CIINTS(DISKH2,CMO,EMY,FCAC,H2AC,WRK,LFREE)
C
C     *** Calculate the diagonal of the CI matrix.
C
      CALL CIDIAG(LSYM,.FALSE.,WRK(KFCAC),WRK(KH2AC),INDXCI,
     &            WRK(KCDIAG),WRK(KW1),LW1)
C     CALL CIDIAG(ICSYM,NOH2,FCAC,H2AC,XNDXCI,DIAGC,WRK,LFREE)
C
      EMY = EMY + POTNUC
      DO 50 I = 0,NCONF-1
         WRK(KCDIAG+I) = WRK(KCDIAG+I) + EMY
   50 CONTINUE
#if defined (VAR_CIHAMTST)
C
C     --- 890131 -- explicit CI Hamiltonian matrix ---
C
      IF ( FLAG(43) ) THEN
         WRITE (LUPRI,'(//A/)') ' CICTL: test call of CIHAM'
         IWAY   = 1
         MXPRDT = MIN(20,NCONF)
         IPRHAM = 10
C
         KIDET  = KW1
         KCIDIA = KIDET  + MXPRDT/IRAT + 1
         KWHAM  = KCIDIA + NCONF
         LWHAM  = LWRK   - KWHAM
         CALL DCOPY(NCONF,WRK(KCDIAG),1,WRK(KCIDIA),1)
         CALL CIHAM(IWAY,MXPRDT,LSYM,WRK(KIDET),WRK(KCIDIA),NCONF,
     *              EMY,WRK(KFCAC),WRK(KH2AC),INDXCI,
     *              WRK(KWHAM),LWHAM,IPRHAM)
C        CALL CIHAM(IWAY,MXPRDT,ICSYM,IDET,CIDIA,NCONF,
C    *              ECORE,FCAC,H2AC,INDXCI,WRK,LWRK,IPRINT)
      END IF
#endif
C
C    Get PHP CI matrix (with correct CI energies)
C
      IF (MAXPHP .GE. 1) THEN
         IPWAY = 1
         IPRPHP = IPRI6 - 3
         CALL PHPGET(LSYM,NCONF,INDXCI,WRK(KFCAC),WRK(KH2AC),
     *               WRK(KCDIAG),EMY,IPWAY,IPRPHP,WRK(KW1),LW1)
      END IF
C
C
C    ******************************
C    *** GET INITIAL CI VECTORS ***
C    ******************************
C
C    CIST constructs NCROOT start vectors for CI optimization.
C    The C-vectors are constructed in CIST or from the result of
C    another run (another geometry).
C
      NCSIM = NCROOT
      CALL CIST(ICICTL,NCSIM,CMO,WRK(KCDIAG),EMY,WRK(KW1),LW1)
      IF (NCSIM .LE. 0) THEN
         WRITE (LUPRI,'(//A,I5)')
     *      ' *** ERROR after CIST, NCSIM =',NCSIM
         CALL QTRACE(LUPRI)
         CALL QUIT('*** ERROR (CICTL) NCSIM .le. 0 after CIST')
      END IF
C
C     *****************************************************************
C     CI-DFT contributions
C     *****************************************************************
C
      IF (DOCISRDFT) THEN
#ifndef MOD_SRDFT
         call quit('srdft not implemented in this version')
#else
         KDVREF = KW2
         KFSR = KDVREF + NNASHX
         KW2A = KFSR   + NNORBT
         LW2A = LWRK   - KW2A
         IF (ISTACI .EQ. 0) THEN
            CALL DZERO(WRK(KDVREF),NNASHX)
         ELSE
            IST = ABS(ISTACI)
            CALL GETDVREF(WRK(KDVREF),IST,INDXCI,WRK(KW2A),LW2A)
         END IF
C
         CALL SRFMAT(WRK(KFSR),CMO,WRK(KDVREF),EJCSR,EJVSR,EDSR,EDFT,
     &               WRK(KW2A),LW2A,IPRI4)
C        ... Correct EMY and CI diagonal
         EMYDFT = EJCSR + EJVSR + EDSR + EDFT
         EMY = EMY + EMYDFT
C        ... extract and save SR over active indices
         CALL GETAC(WRK(KFSR),WRK(KSRAC))
C        ... add the contributions over the active indices
         CALL DAXPY(NNASHX,D1,WRK(KSRAC),1,WRK(KFCAC),1)
C        ... fold SRAC for 
C         1) Adding it to CI-diagonal to get better pre-conditioning
C         2) Later printing of CIDFT specific energy contributions
         CALL DSPTGE(NASHT,WRK(KSRAC),WRK(KFSR))
         CALL DGEFSP(NASHT,WRK(KFSR),WRK(KSRAC))
         ESRDV = DDOT(NNASHX,WRK(KDVREF),1,WRK(KSRAC),1)
         DO I = 0,NCONF-1
            WRK(KCDIAG+I) = WRK(KCDIAG+I) + EMYDFT + ESRDV
         END DO
#endif
      ENDIF
C
C    *******************************
      IF (MAXITC .LE. 0) GO TO 9999
C    *******************************
      CINMEM = (NCSIM .LE. MAXSIM)
C
      IF (NCSIM .GT. MXCIRT) THEN
         WRITE (LUPRI,'(//A,I5,A,I5)')
     *   ' *** ERROR after CIST, NCSIM =',NCSIM,
     *   ', maximum (MXCIRT) =',MXCIRT
         CALL QTRACE(LUPRI)
         CALL QUIT('*** ERROR (CICTL) MXCIRT parameter .lt. NCSIM')
      END IF
      DO 100 I = 1,NCSIM
         JCROOT(I) = I
         WRK(KRESID + I-1) = D0
  100 CONTINUE
C
      REWIND LUIT1
      CALL MOLLAB('STARTVEC',LUIT1,LUPRI)
      REWIND LUIT3
      MXCSIM = 2*MAXSIM - 1
C     reuse space for orthonormalization,
C     except reserving last NCONF elements for CPREV
      KCPREV = KCVECS + MXCSIM*NCONF
      NCPREV = 0
      DO 490 I = 1,NCSIM,MXCSIM
         NXSIM = MIN(MXCSIM,(NCSIM-I+1))
         JCVECS = KCVECS
         DO 420 J = 1,NXSIM
            CALL READT(LUIT1,NCONF,WRK(JCVECS))
            JCVECS = JCVECS + NCONF
  420    CONTINUE
         THRLDV = NCONF * THRLDP
         CALL CIORT(NXSIM,NCPREV,LUIT3,JCROOT(NCPREV+1),
     *              WRK(KCVECS),THRLDV,WRK(KCPREV))
  490 CONTINUE
      IF (NCPREV .NE. NCSIM) THEN
         WRITE (LUPRI,'(//A,I8,A)')
     *      ' *** ERROR ***',(NCSIM-NCPREV),' OF THE TRIAL CI VECTORS'//
     *      ' WERE LINEAR DEPENDENT (CICTL).'
         CALL QUIT('***ERROR (CICTL) LIN.DEP. AMONG INITIAL CI '//
     &             'TRIAL VEC.S')
      END IF
C
C     *****************************************************************
C     CI iteration loop
C     *****************************************************************
C
C
      TIMLIN = D0
      TIMCIT = SECOND()
      NCLIN  = 0
      ITCI   = 0
C
  500 ITCI = ITCI + 1
         IF (IPRI4 .GT. 5)
     *      WRITE (LUW4,'(//A,I4,/)') ' *** CI ITERATION NO.',ITCI
        IF (LUPRI.NE.LUW4 .AND. IPRI6 .GE. 5) WRITE (LUPRI,'(//A,I4,/)')
     *      ' *** CI ITERATION NO.',ITCI
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         REWIND LUIT5
         DO 220 I = 1,NCIRED
  220       READ (LUIT5)
C
C        REPEAT UNTIL (NCSIM .LE. 0)
C
  240    CONTINUE
            NXSIM = MIN(MAXSIM,NCSIM)
            IF (.NOT. CINMEM) THEN
               REWIND LUIT3
               DO 310 I = 1,NCIRED
  310             READ (LUIT3)
               DO 320 I = 1,NXSIM
                  JCVECS = KCVECS + (I-1)*NCONF
                  CALL READT(LUIT3,NCONF,WRK(JCVECS))
  320          CONTINUE
            END IF
            TIMLIN = TIMLIN - SECOND()
            NCLIN  = NCLIN + NXSIM
            CALL CILIN(NXSIM,WRK(KCVECS),WRK(KFCAC),WRK(KH2AC),
     *                 INDXCI,WRK(KSVECS),LSVECS)
            TIMLIN = TIMLIN + SECOND()
            DO 340 I = 1,NXSIM
               JSVECS = KSVECS + (I-1)*NCONF
               CALL WRITT(LUIT5,NCONF4,WRK(JSVECS))
  340       CONTINUE
            NCIRED = NCIRED + NXSIM
            CALL CIRED(1,NCIRED,WRK(KREDCI),NXSIM,WRK(KCVECS),
     *                 WRK(KSVECS),DUMMY,DUMMY,WRK(KW3),LW3)
            NCSIM = NCSIM - NXSIM
            IF (NCSIM .GT. 0) GO TO 240
C        ^----------------------------- END REPEAT
C
         CALL CIRED(2,NCIRED,WRK(KREDCI),0,WRK(KCVECS),
     *              DUMMY,WRK(KEVAL),WRK(KEVEC),WRK(KSVECS),LSVECS)
C        CALL CIRED(ICTL,NCIRED,REDCI,NCNEW,CVECS,SVECS,EVAL,EVEC,WRK,LWRK)
         IF (P6FLAG(10) .OR. P6FLAG(34)) THEN
            WRITE (LUPRI,'(/A,I3)')
     *         ' The reduced CI matrix, iteration',ITCI
            IF (P6FLAG(34)) CALL OUTPAK(WRK(KREDCI),NCIRED,-1,LUPRI)
            MCROOT = MIN(NCROOT+2,NCIRED)
            WRITE (LUPRI,'(/I4,A//,(10X,1P,6D15.6))')
     *         MCROOT,' lowest eigenvalues of reduced CI matrix :',
     *         (WRK(KEVAL-1+I),I=1,MCROOT)
            IF (P6FLAG(34) .OR. IPRI6 .GT. 10) THEN
               WRITE (LUPRI,'(/A,I3,A)')
     *            ' - and the lowest',NCROOT,' eigenvectors :'
               CALL OUTPUT(WRK(KEVEC),1,NCIRED,1,NCROOT,
     *                     NCIRED,NCIRED,-1,LUPRI)
            END IF
         END IF
         IF (P6FLAG(51)) THEN
C           ... check symmetry of reduced CI matrix and CI overlap
            CALL CIRCHK(WRK(KSVECS),LSVECS)
         END IF
         ESHIFT = D0
         NCSIM  = NCROOT
         CALL CINEX(WRK(KEVAL),WRK(KEVEC),EMY,ESHIFT,WRK(KCDIAG),
     *              NCSIM,JCONV,WRK(KRESID),WRK(KW2),LW2)
      IF (JCONV .GT. 0) THEN
         IF (NCSIM .EQ. 0) THEN
            GO TO 8100
         ELSE IF (NCSIM .LT. 0) THEN
            CINMEM = .FALSE.
            NCSIM  = -NCSIM
         ELSE
            CINMEM = (NCSIM .LE. MAXSIM)
         END IF
         IF (ITCI         .GE. MAXITC) GO TO 8000
         IF (NCIRED+NCSIM .GE. MAXRC ) GO TO 8200
         GO TO 500
      END IF
      WRITE (LUW4,'(//A,I3,A/)')
     *   ' *** CI converged in',ITCI,' iterations.'
      IF (LUPRI .NE. LUW4) WRITE (LUPRI,'(//A,I3,A/)')
     *   ' *** CI converged in',ITCI,' iterations.'
C
C
C     *****************************************************************
C     Save final result
C     *****************************************************************
C
 6000 CONTINUE
      TIMCIT = SECOND() - TIMCIT
      DO 810 I = 1,NCROOT
         ECI(I) = EMY + WRK(KEVAL-1+I)
  810 CONTINUE
      IF (DOCISRDFT) THEN
         WRITE(LUPRI,'(/A,I5)') 'SR 0-el. energy for input .STATE',
     &                            ISTACI
         WRITE(LUPRI,'(A/)') '-------------------------------------'
         WRITE(LUPRI,805) '  SR core    Hartree energy   :',EJCSR
         WRITE(LUPRI,805) '- SR valence Hartree energy   :',EJVSR
         WRITE(LUPRI,805) '+ SR Exchange-correlation     :',EDFT
         WRITE(LUPRI,805) '- SR Exchange-correlation pot.:',EDSR
         WRITE(LUPRI,806) '= Total eff. SR 0-el. energy  :',EMYDFT
C
       DO IST = 1,NCROOT
         KDVREF = KW2
         KW2A   = KDVREF + NNASHX
         LW2A   = LWRK   - KW2A
         CALL GETDVREF(WRK(KDVREF),IST,INDXCI,WRK(KW2A),LW2A)
C        ... DVREF is equal to DV for state abs(ISTI),
C            i.e. 1-el. energy can be calculated
         ESRDV = DDOT(NNASHX,WRK(KDVREF),1,WRK(KSRAC),1)
         WRITE(LUPRI,'(/A,I5)') 'CI-DFT energy for state no.',IST
         WRITE(LUPRI,'(A/)') '-------------------------------------'
         WRITE(LUPRI,805) 'SR eff. 1-el. energy          :',ESRDV
         EJTOT = ESRDV + EDSR + EJVSR + EJCSR
         WRITE(LUPRI,805) 'SR total Hartree energy       :',EJTOT
         ETDFT = EMYDFT + ESRDV
         WRITE(LUPRI,805) 'SR eff. total DFT energy      :',ETDFT
         ELRCI = ECI(IST) - ETDFT
         WRITE(LUPRI,805) 'LR total CI  energy           :',ELRCI
         WRITE(LUPRI,806) 'Total CI-DFT energy           :',ECI(IST)
        END DO
  805   FORMAT( 1X,A,F25.12)
  806   FORMAT(/1X,A,F25.12)
      ENDIF
      J = LUW4
  815 WRITE (J,'(//A/,(I5,0P,F25.15,1P,D15.2))')
     *   ' CI energies and residuals:',
     *   (I,ECI(I),WRK(KRESID+I-1), I = 1,NCROOT)

      IF(NCROOT .EQ. 1) THEN 
C JT 22-03-05 : add pretty output line if only one root
        WRITE (J,'(//A,(F25.15))')' Final CI energy:',ECI(1)
      ENDIF

      IF (J .NE. LUPRI) THEN
         J = LUPRI
         GO TO 815
      END IF
C
C     Save CI vectors and orbitals on LUIT1
C
      CALL CISAVE(NCROOT,NCIRED,ECI,WRK(KRESID),WRK(KEVEC),CMO,INDXCI,
     *            WRK(KW2),LW2)
C     CALL CISAVE(NCROOT,NCIRED,ECI,RESID,EVEC,CMO,INDXCI,WRK,LFREE)
C
      ICONV = 0
      J = 0
      DO 820 I = 1,NCROOT
         ICROOT(I) = 0
         IF (JCROOT(I) .EQ. 0) THEN
            ICONV = ICONV + 1
         ELSE
            J = J + 1
            ICROOT(J) = JCROOT(I)
         END IF
  820 CONTINUE
C
      IF (FLAG(67)) THEN
         IPRNO = MAX(IPRI4,4)
         KCVEC = KW2
         KOCCNO= KCVEC + NCONF
         KUNO  = KOCCNO+ NOCCT
         KW3   = KUNO  + N2ASHX
         LW3   = LWRK  - KW3
C
         REWIND LUIT1
         CALL MOLLAB('STARTVEC',LUIT1,LUPRI)
      DO 900 JSTATE = 1,NCROOT
         CALL READT(LUIT1,NCONF,WRK(KCVEC))
C
         WRITE (LUW4,'(//A)')
     *      ' >>> SIRCI.CICTL natural orbital analysis <<<'
         IF (JSTATE .EQ. ISTACI) THEN
            WRITE (LUW4,2000) JSTATE
            IF (LUPRI .NE. LUW4) WRITE (LUPRI,2000) JSTATE
            NATORB = FLAG(68)
         ELSE
            WRITE (LUW4,2001) JSTATE
            IF (LUPRI .NE. LUW4) WRITE (LUPRI,2001) JSTATE
            NATORB = .FALSE.
         END IF
 2000    FORMAT(/' CI eigenvector no.',I4,' (= the reference state)',
     *          /' ------------------------------------------------')
 2001    FORMAT(/' CI eigenvector no.',I4,
     *          /' ----------------------')
C
         IF (P6FLAG(32) .OR. NATORB) THEN
C           NATORB is true if CINO and reference state
            WRITE (LUPRI,'(/A,F8.5)')
     &         ' Printout of CI-coefficients abs greater than:',THRPWF
            CALL PRWF(WRK(KCVEC),INDXCI,LUPRI,THRPWF)
C           CALL PRWF(C,INDXCI,IOUT,THRPWF)
         END IF
C
C        MAERKE insert calculation of spin-density matrix here,
C        if ISPIN .ne. 1
C
C        insert option for calculation of two-electron density matrix,
C        using MAKDM, for deeper analysis of state vectors.
C
         IF (NATORB) THEN
            KCMO = KW3
            KW4  = KCMO + NCMOT
            CALL DCOPY(NCMOT,CMO,1,WRK(KCMO),1)
         ELSE
            KCMO = KW3
            KW4  = KCMO
         END IF
         KWNO = 1
         LWNO = LWRK - KW4 + 1
         CALL GETNO(WRK(KCVEC),LSYM,WRK(KOCCNO),WRK(KUNO),WRK(KCMO),
     *              .TRUE.,NATORB,.FALSE.,MOISNO,INDXCI,
     *              WRK(KW4),KWNO,LWNO,LUW4,IPRNO)
C        CALL GETNO(CVEC,ICSYM,OCCNO,UNO,CMO,ORDER,NATORB,NOAVER,
C    *              MOISNO,INDXCI,WRK,KFRSAV,LFRSAV,LUPRI,IPRINT)
         IF (NATORB) THEN
            WRITE (LUW4,'(/A/)')
     *         ' The natural orbitals are saved with label "NEWORB  ".'
            CALL NEWORB('CINOSAVE',WRK(KCMO),.TRUE.)
         END IF
  900 CONTINUE
      END IF
C
C     Generate density and dump to file.
C
C      CALL MCRHO(CMO,INDXCI,WRK(KW2),LW2,IPRSIR)
C

      IF (FLAG(4) .AND. FLAG(25)) THEN
C        ... if (cicalc and save for response) then
         IF (DOCISRDFT) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//1X,A/)')
     &         'WARNING, CI RESPONSE is not possible after CI-srDFT'
         ELSE IF (FLAG(67) .AND. FLAG(68)) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//1X,A/)')
     &         'WARNING, CI RESPONSE is not possible after CI-NO'
         ELSE
            IF (.NOT. FLAG(14)) THEN
C              ... orbitals have been changed by natural orbital
C                  transformation in CISAVE
               WRITE (LUPRI,'(//A/A/A)') ' Molecular orbitals'//
     *            ' were modified after CI convergence',
     *            ' integrals are transformed to new orbital basis',
     *            ' before saving information for post-programs'
               CALL TRACTL(0,CMO,WRK,LWRK)
C              CALL TRACTL(ITRLVL,CMO,WRK,LWRK)
               FLAG(14) = .TRUE.
C
C              *** Obtain FCAC and H2AC integral matrices
C
               CALL CIINTS(DISKH2,CMO,
     *                     EMY,WRK(KFCAC),WRK(KH2AC),WRK(KW1),LW1)
C              CALL CIINTS(DISKH2,CMO,EMY,FCAC,H2AC,WRK,LFREE)
C
C              *** Calculate the diagonal of the CI matrix.
C
               CALL CIDIAG(LSYM,.FALSE.,WRK(KFCAC),WRK(KH2AC),INDXCI,
     &                     WRK(KCDIAG),WRK(KW1),LW1)
C              CALL CIDIAG(ICSYM,NOH2,FCAC,H2AC,XNDXCI,DIAGC,WRK,LFREE)
C
               EMY = EMY + POTNUC
            ELSE
C              subtract EMY from diagonal again, as expected in CIRSPS
               DO 950 I = 0,NCONF-1
                  WRK(KCDIAG+I) = WRK(KCDIAG+I) - EMY
  950          CONTINUE
            END IF
C
            CALL CIRSPS(ECI(ISTACI),EMY,CMO,WRK(KFCAC),WRK(KH2AC),
     *                  WRK(KCDIAG),INDXCI,WRK(KW2),LW2)
C           CALL CIRSPS(ECIREF,EMY,CMO,FCAC,H2AC,HDIAG,INDXCI,WRK,LFREE)
         END IF
      END IF
C
      IF (IPRI6 .GT. 0) THEN
         TIMTOT = SECOND() - TIMTOT
         WRITE (LUPRI,'(2(/A,F10.2),2(/A,I3,A,F10.2))')
     *      ' Total time used in CICTL                   :',TIMTOT,
     *      '   Integral transformation                  :',TIMITR,
     *      '   Total for',ITCI,' CI iterations               :',TIMCIT,
     *      '   hereof used in',NCLIN,' linear transformations :',TIMLIN
      END IF
C
 9999 CONTINUE
C
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
      CALL QEXIT('CICTL ')
      RETURN
C
C
C
 8000 CONTINUE
      WRITE (LUPRI,'(//A,I5/I15,A)')
     *   ' *** Reached maximum number of CI iterations:',MAXITC,
     *   JCONV,' CI roots are not converged.'
      GO TO 6000
C
 8100 CONTINUE
      NWARN = NWARN + 1
      WRITE (LUPRI,'(//A/I15,A)')
     *   ' *** WARNING, all CI trial vectors linear dependent',
     *   JCONV,' CI roots are not converged.'
      GO TO 6000
 8200 CONTINUE
      NWARN = NWARN + 1
      WRITE (LUPRI,'(//A,I5/I15,A)')
     *   ' *** WARNING, reached maximum dimension of reduced matrix:',
     *   MAXRC,JCONV,' CI roots are not converged.'
      GO TO 6000
C
C     *****************************************************************
C     End of CICTL
C     *****************************************************************
C
      END
C  /* Deck ciints */
      SUBROUTINE CIINTS(DISKH2,CMO,EMY,FCAC,H2AC,WRK,LFREE)
C
C Written by Hans Joergen Jensen 13-Dec-1984
C Revision for LUCITA Dec 2010 /HJAaJ (moved CALL CIDIAG outside)
C
C Purpose:
C  CIIN*: Calculate integrals needed for CI. EMY, FCAC and H2AC
C         are needed by calling program.
C
#include <implicit.h>

      DIMENSION CMO(*), FCAC(*),H2AC(*),WRK(*)
      LOGICAL   DISKH2

#include <priunit.h>
#include <dummy.h>
C
C Used from common blocks:
C   INFORB : NNASHX,NNBAST,NNORBT
C   INFINP : SRHYBR
C
#include <maxorb.h>
#include <inforb.h>
#include <infinp.h>

      LOGICAL SRHYBR_SAVE
C
      CALL QENTER('CIINTS')
C
C *** Section 0 ***
C
C work space allocation:
C
      NF = 1
      SRHYBR_SAVE = SRHYBR
      IF (SRHYBR) THEN
         SRHYBR = .FALSE.
         WRITE(LUPRI,'(//A)') ' INFO: .SRHYBR ignored in CIINTS'
!        NF = 3
      END IF
      KFC   = 1
      KW0   = KFC   + NF*NNORBT
      LW0   = LFREE - KW0
C
      KW1   = 1
      LW1   = LFREE - KW1
C
      KEND  = MAX(KW0,KW1)
      IF (KEND.GT.LFREE) CALL ERRWRK('CIINTS',KEND,LFREE)
C
C *** Section 1 ***
C
C     Calculate F-matrix FCAC,
C     and integrals H2AC.
C
      CALL FCKMAT(.TRUE.,DUMMY,CMO,EMY,WRK(KFC),DUMMY,
     &             WRK(KW0),LW0)
C     CALL FCKMAT(ONLYFC,DV,CMO,EMCMY,FC,FV,WRK,LFREE)
      CALL GETAC(WRK(KFC),FCAC)
      IF (DISKH2) THEN
         CALL CIIN2B(H2AC,WRK(KW1),LW1)
      ELSE
         CALL CIIN2A(H2AC,WRK(KW1),LW1)
      END IF
C
C *** End of subroutine CIINTS
C
      SRHYBR = SRHYBR_SAVE
      CALL QEXIT('CIINTS')
      RETURN
      END
C  /* Deck ciin2a */
      SUBROUTINE CIIN2A(H2AC,WRK,LWRK)
C
C Written by Hans Agren 27-Jun-1987
C Rewritten 891017 hjaaj using NXTH2M
C
C Revisions:
C  Purpose : To obtain the necessary two-electron
C            integrals to perform a CI calculation.
C
#include <implicit.h>
C
      DIMENSION H2AC(NNASHX,NNASHX),WRK(LWRK)
C
C Used from common blocks:
C   INFORB : NNASHX,NNBAST,?
C   INFIND : ISMO,?
C CBGETDIS : IADH2
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <cbgetdis.h>
C
      DIMENSION NEEDMU(6)
      DATA NEEDMU /0,0,1,0,0,0/
C     ... only active-active distributions needed (type 3)
      CALL QENTER('CIIN2A')
C
C    get array of 2-electron integrals with active indices
C    for ci-use.
C
C
      IF (NASHT .LE. 1) THEN
         CALL QTRACE(LUPRI)
         CALL QUIT('PROGRAM ERROR, CIIN2A called with NASHT .le. 1')
      END IF
C
C *** Initialize H2AC ***
C
      CALL DZERO(H2AC,NNASHX*NNASHX)
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KH2CD,N2ORBX,WRK,KFREE,LFREE)
C
C ****************************************************************
C     Loop over active-active distributions
C
      IDIST = 0
  100 CALL NXTH2M(IC,ID,WRK(KH2CD),NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
C     We know that IC and ID are both active orbitals
C     because we have specified so in NEEDMU.
C
#if defined (VAR_DEBUG)
      IF (IOBTYP(IC).NE.JTACT .OR. IOBTYP(ID).NE.JTACT) THEN
         CALL QUIT('CIIN2A: NXTH2M did not return '//
     &             'active-active distribution')
      END IF
#endif
         IF (IC.LT.ID) THEN
            ISWAP = IC
            IC = ID
            ID = ISWAP
         END IF
         ICSYM = ISMO(IC)
         IDSYM = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         NCW  = ICH(IC)
         NDW  = ICH(ID)
         NCDW = IROW(NCW) + NDW
C
C        Add active-active elements to H2AC
C
         CALL ADH2AC(H2AC(1,NCDW),WRK(KH2CD),ICDSYM)
C
C        Go get next active-active distributions
C
      GO TO 100
C
C     arrive at 800 when finished with all active-active distributions
C
  800 CONTINUE
C
C *** IADH2 .lt. 0 means that H2AC is in core.
C
      IADH2 = -1
C
C
      IF (P6FLAG(30)) THEN
         NCDW = 0
         DO 300 NCW = 1,NASHT
         DO 300 NDW = 1,NCW
            NCDW = NCDW + 1
            WRITE(LUPRI,3100)NCW,NDW
            CALL OUTPAK(H2AC(1,NCDW),NASHT,-1,LUPRI)
  300    CONTINUE
      END IF
 3100 FORMAT(/' CIIN2A: Two-electron integrals H2AC, distribution',
     *        2I4)
C     end of CIIN2A
C
      CALL QEXIT('CIIN2A')
      RETURN
      END
C  /* Deck ciin2b */
      SUBROUTINE CIIN2B(H2ACCD,WRK,LWRK)
C
C 27-6 1987 Hans Agren
C
C Purpose:
C
C    To get array of 2-electron integrals with active indices
C    and to store them on direct access file (LUDA2=24), one
C    CD-distribution per record.
C
C    Special SIRIUS ordering of MO integrals assumed (only
C    one distribution in each record).
C
C
C
C Scratch: H2ACCD (active 2-el integrals)
C          WRK
C
C ****************************************************************
#include <implicit.h>
      DIMENSION H2ACCD(NNASHX), WRK(LWRK)
#include <iratdef.h>
C
C
C Used from common blocks:
C   INFORB : NNASHX,...
C   INFIND : ISMO,ICH
C   INFTAP : LUH2AC
C   INFPRI : P6FLAG
C CBGETDIS : IADH2
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
      LOGICAL OLDDX
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
#include <cbgetdis.h>

      DIMENSION NEEDMU(6)
      DATA NEEDMU /0,0,1,0,0,0/
C     ... only active-active distributions needed (type 3)

      CALL QENTER('CIIN2B')
      IF (NASHT .LE. 1) THEN
         CALL QTRACE(LUPRI)
         CALL QUIT('PROGRAM ERROR, CIIN2B called with NASHT .le. 1')
      END IF
C
C     Define off-set on LUH2AC for H2AC.
C     (IADH2 .lt. 0 means that H2AC is in core.)
C
C     MAERKE 900302: lav mekanisme til at samordne allok. af
C     H2AC, H2XAC, H2DAC.
C
      IADH2 = 0
C
C open direct access file
C
      CALL GPOPEN(LUH2AC,'MOH2AC.DA','NEW','DIRECT',' ',NNASHX*IRAT,
     &            OLDDX)
C
C
C ****************************************************************
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KH2CD,N2ORBX,WRK,KFREE,LFREE)
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDMU(6)
C
      IDIST = 0
  100 CALL NXTH2M(IC,ID,WRK(KH2CD),NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
         IF (IC.LT.ID) THEN
            ISWAP = IC
            IC    = ID
            ID    = ISWAP
         END IF
         ICSYM  = ISMO(IC)
         IDSYM  = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         NCW    = ICH(IC)
         NDW    = ICH(ID)
         NCDW   = IROW(NCW) + NDW
C
C        Collect active-active elements in H2ACCD
C
         CALL DZERO(H2ACCD,NNASHX)
         CALL ADH2AC(H2ACCD,WRK(KH2CD),ICDSYM)
C
C        Write this CD-distribution on disk:
C
         CALL WRITDX(LUH2AC,IADH2+NCDW,IRAT*NNASHX,H2ACCD)
         IF (P6FLAG(30)) THEN
            WRITE (LUPRI,3100) NCW,NDW
            CALL OUTPAK(H2ACCD,NASHT,-1,LUPRI)
         END IF
 3100    FORMAT(/' CIIN2B: Two-electron integrals H2AC, distribution',
     *           2I4)
C
C        Go get next active-active Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all needed Mulliken distributions
C
  800 CONTINUE
C
C *** end of subroutine CIIN2B
C
      CALL QEXIT('CIIN2B')
      RETURN
      END
C  /* Deck cilin */
      SUBROUTINE CILIN(NCSIM,BCVECS,FCAC,H2AC,INDXCI,WRK,LFREE)
C
C Written 30-Apr-1985 by Hans Jorgen Aa. Jensen
C (based on LINTRN)
C
C Purpose:
C  Do direct CI: calculate the simultaneous linear transformation
C  defined by the CI matrix H of the NCSIM(ultaneous) B CI vectors :
C
C  SCVECS(i,j)  =  sum(k=1,NCONF) H(i,k) * BCVECS(k,j)
C                  (i = 1,NCONF; j = 1,NCSIM)
C
C Input:
C
C Output:
C  SCVECS; will start at WRK(1) on output; followed by
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION BCVECS(NCONF,*), FCAC(*), H2AC(*)
      DIMENSION INDXCI(*),WRK(LFREE)
C
C *** local constants
C
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C   INFVAR : NCONF
C CBGETDIS : DISTYP,IADINT,IADH2
C
#include <maxorb.h>
#include <cbgetdis.h>
C
      CALL QENTER('CILIN ')
C
C *** Allocate work area for CISIGC
C
      KSCVEC = 1
      KCW    = KSCVEC + NCSIM*NCONF
      LCW    = LFREE  - KCW
C--
      IF (LCW .LT. 0) CALL ERRWRK('CILIN',-KCW,LFREE)
C
C *** call CISIGC to calculate CI sigma vector
C
      DISTYP = 1
      IADINT = IADH2
      CALL CISIGC(NCSIM,BCVECS,WRK(KSCVEC),NCONF,
     *            FCAC,H2AC,INDXCI,WRK(KCW),LCW)
C     CALL CISIGC(NSIM,BCVECS,SCVECS,LSCVEC,FCAC,H2AC,INDXCI,WRK,LFREE)
C
C *** MAERKE implement solvent contribution.
C
C *** end of SUBROUTINE CILIN
C
      CALL QEXIT('CILIN ')
      RETURN
      END
C  /* Deck cinex */
      SUBROUTINE CINEX(EVAL,EVEC,EMY,ESHIFT,CDIAG,
     *                 NCSIM,JCONV,RESID,WRK,LFREE)
C
C Written 21-Jun-1985 by Hans Jorgen Aa. Jensen
C (based on NEONEX version 2)
C Revisions:
C   1-Jul-1987 hjaaj: CTRUNC and SYMTRZ.
C
C Purpose:
C  Find the next trial vectors for the CI optimization.
C
C Input:
C  EVAL(NCIRED),        eigenvalues of reduced CI matrix.
C  EVEC(NCIRED,NCIRED), corresponding eigenvectors.
C  ESHIFT,              shift of "Davidson" denominator (usually zero)
C  CDIAG(*),            diagonal of CI matrix + any PHP information
C  NCSIM                = NCROOT, number of roots we want to converge to
C
C
C Output:
C  JCONV number of roots not converged
C  NCSIM abs value is number of CI trial vectors written to LUIT3
C        if positive they are returned in WRK(1)
C
C  note: if NCSIM .eq. 0 and JCONV .gt. 0, then not converged
C  ----  and no new trial vectors could be found with the algorithm
C        used (e.g. all removed because of linear dependency).
C
C Scratch:
C  WRK(LFREE)
C
#include <implicit.h>
      DIMENSION EVAL(*), EVEC(NCIRED,*), CDIAG(*), RESID(*), WRK(*)
C
C     ********** Define parameters, common blocks, and local variables.
C
C
      PARAMETER ( THRAUX = 0.1D0 )
      PARAMETER ( THKC   = 0.10D0, THKO   = 0.001D0 )
      PARAMETER ( THLIU1 = 0.01D0, THLIU2 = 0.0D0 )
      PARAMETER ( THRSIM = 1.D-6 , DTEST  = 0.01D0 )
#include <thrldp.h>
      PARAMETER ( D0=0.0D0, D1 = 1.0D0, DHALF=0.5D0 )
C
C Common blocks:
C  NJCR, number of residual vectors to construct
C          (number of eigenvalues/-vectors of REDL to use)
C Used from common blocks:
C   CICB01 : NJCR,NKCR,JCROOT(*),KCROOT(*),NCIRED,THRCIT
C   INFINP : FLAG(*),JOLSEN,RESPHP,ISTACI
C   INFVAR : NCONF
C   INFTAP : LUIT3,LUIT5
C   INFPRI : P4FLAG(*)
C   PRIUNIT: IPRSTAT
C
#include <maxorb.h>
#include <priunit.h>
#include <cicb01.h>
#include <infinp.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
C
C     local variables:
C
      LOGICAL FIRST, ALLRES, WCRESD, CTRUNC, SYMTRZ
C
C
C     data statement:
C
      SAVE FIRST
      DATA FIRST /.TRUE./
C
C     **********
C     ********** End of declarations, start execution.
C     ********** First, define some entities we will need.
C     **********
C
      CALL QENTER('CINEX ')
      IF (FIRST) THEN
         WRITE (LUW4,'(//A,F15.8/)')
     *      ' Convergence threshold for CI optimization :',THRCIT
         IF (LUPRI .NE. LUW4) WRITE (LUPRI,'(//A,F15.8/)')
     *      ' Convergence threshold for CI optimization :',THRCIT
         IF (.NOT.DOCISRDFT .AND. ISTACI.GT.0 .AND. NROOCI .GT. 1) THEN
            WRITE (LUW4,1010) THRAUX
            IF (LUPRI.NE.LUW4) WRITE (LUPRI,1010) THRAUX
         END IF
         TLINDP = 100*NCONF*THRLDP
         TLINDP = SQRT(TLINDP)
         IF (THRCIT .LT. TLINDP) THEN
            WRITE(LUPRI,'(//A/A,2(/A,D10.2))')
     *      ' CINEX ERROR: Convergence threshold is less than the',
     *      '              numerical linear dependence limit',
     *      ' Convergence threshold   ',THRGRD,
     *      ' Linear dependence limit ',TLINDP
            CALL QUIT('CINEX:convergence threshold below lin.dep.limit')
         END IF
         FIRST = .FALSE.
      END IF
 1010 FORMAT(' Convergence threshold for auxiliary roots:',F15.8/)
      NCIT3 = NCIRED
C
      WCRESD = FLAG(76)
C     ... use energy weighted residuals for convergence
C         this gives more uniform accuracy in state vectors
      CTRUNC = FLAG(77)
C     ... "truncate ci trial vectors", i.e. zero small elements.
      SYMTRZ = FLAG(78)
C     ... symmetrize CI vector
C
C     ALLRES = FLAG(59)
      ALLRES = .TRUE.
C     ... if flag(59) check convergence of all roots,,
C         roots may have switched order if CI matrix has symmetry
C         890915-hjaaj: ALLRES true always now.
C
      IF (.NOT.ALLRES) NCSIM = NJCR
C     ... else NCSIM = NCROOT (set by calling routine)
      DO 1 I = 1,NCSIM
         KCROOT(I) = JCROOT(I)
         IF (ALLRES) JCROOT(I) = I
    1 CONTINUE
C
C     **********
C     ********** Find MAXRCV: how many we can treat simultaneously.
C     **********
C     MWRK0 is space needed for arrays independent of MAXRCV.
C     MWRK1 is space needed per residual vector.
C
      MWRK0 = NCONF * 2
      IF (JOLSEN) THEN
         MWRK1 = NCONF * 2
C        solution vector plus residual vector
      ELSE
         MWRK1 = NCONF
C        only residual vector
      END IF
      MAXRCV = MIN(NCSIM, (LFREE-MWRK0) / MWRK1)
      IF (IPRSTAT .GE. 5) THEN
         WRITE (LUSTAT,'(/,(A,I10))')
     *      ' (CINEX) number of residual vectors',NCSIM,
     *      '         number of residuals/batch ',MAXRCV
         WRITE (LUSTAT,'(/A,I4,A/,(10I5))') '         Index of',NJCR,
     *      ' non-conv. vectors from last iteration:',
     *      (KCROOT(I),I=1,NJCR)
      END IF
      IF (MAXRCV.LE.0) CALL ERRWRK('CINEX',(MWRK0+MWRK1),LFREE)
C
      NJCR = NCSIM
C     ... set NJCR to actual number of residual vectors
C         (if ALLRES it is larger than no. non-conv. vectors)
C
C     CHECK IF WE SHOULD INCREASE NCSIM (Liu simultaneous expansion)
C     BECAUSE     EVAL(last root to converge)
C     IS CLOSE TO EVAL(first root not to converge).
C
C     CRIT1 / THLIU1 : for absolute test
C     CRIT2 / THLIU2 : maybe for relative test, we haven't decided yet
C                      what the algorithm should look like. 851213/hj
C
      MJCR   = NJCR
      LSTJCR = 0
      DO 2 I = 1,NJCR
         LSTJCR = MAX(LSTJCR,JCROOT(I))
    2 CONTINUE
      IJCR = LSTJCR
    3 IF (IJCR .LT. NCIRED .AND. IJCR .LT. MXCIRT) THEN
         CRIT1 = EVAL(IJCR+1) - EVAL(IJCR)
C        CRIT2 = ABS(EVAL(IJCR)) / (EVAL(NCIRED) - EVAL(1)) ???
         IF ( CRIT1 .LT. THLIU1 ) THEN
C           here ought to be a symmetry test of IJCR vs. IJCR + 1
            MJCR = MJCR + 1
            IJCR = IJCR + 1
            JCROOT(MJCR) = IJCR
            GO TO 3
         END IF
      END IF
      IF (IPRSTAT .GT. 0 .AND. MJCR .GT. NJCR) THEN
         WRITE (LUSTAT,'(//A,2(/A,I4),//A)')
     *      ' -- CINEX, number of trial vectors increased --',
     *      ' Number of CI roots                   :',NJCR,
     *      ' Number of simultaneous trial vectors :',MJCR,
     *      ' Lowest eigenvalues of reduced CI matrix:'
         I = MIN(NCIRED, MJCR+2)
         WRITE (LUSTAT,'(5X,5F12.8)') (EVAL(J), J = 1,I)
      END IF
C
C     **********
C     ********** Loop over batches
C     **********
C
      ISTCNV = 0
      NCSIM  = 0
      JSIMX  = 0
      NSIML  = MJCR
      IBATCH = 0
   10 CONTINUE
      IBATCH = IBATCH + 1
      NSIMX  = MIN(MAXRCV,NSIML)
C
C     **********
C     ********** Calculate residual vectors.
C     **********
C
C     ********** a) core allocation
C       RVCS for residual vectors
C
      MSIMX = 0
      DO 205 ISIMX = JSIMX+1, JSIMX+NSIMX
         IF (JCROOT(ISIMX) .NE. 0) THEN
            MSIMX = MSIMX + 1
            JCROOT(JSIMX+MSIMX) = JCROOT(ISIMX)
         END IF
  205 CONTINUE
      DO 206 ISIMX = JSIMX+MSIMX+1,JSIMX+NSIMX
         JCROOT(ISIMX) = 0
  206 CONTINUE
C
      KRVECS = 1
      KSCR   = KRVECS + MSIMX*NCONF
      IF (JOLSEN) THEN
         KXVECS = KSCR
         KSCR   = KXVECS + MSIMX*NCONF
      ELSE
         KXVECS = KRVECS
      END IF
      CALL DZERO (WRK(KRVECS),(KSCR-KRVECS))
C
C     ********** b) calculate RVCS
C
C     ********** sigma vector part
C
      REWIND LUIT5
      DO 220 ICIRED = 1,NCIRED
         JRVECS = KRVECS
         CALL READT(LUIT5,NCONF,WRK(KSCR))
         DO 210 ISIMX = JSIMX+1,JSIMX+MSIMX
            FAC = EVEC(ICIRED,JCROOT(ISIMX))
            CALL DAXPY(NCONF,FAC,WRK(KSCR),1,WRK(JRVECS),1)
            JRVECS = JRVECS + NCONF
  210    CONTINUE
  220 CONTINUE
C
C     ********** subtract E times eigenvector to finish residual
C     note that if JOLSEN false, then KXVECS = KRVECS.
C
      REWIND LUIT3
      DO 330 ICIRED = 1,NCIRED
         CALL READT(LUIT3,NCONF,WRK(KSCR))
         JXVECS = KXVECS
         DO 310 ISIMX = JSIMX+1, JSIMX+MSIMX
            JCROTX = JCROOT(ISIMX)
            IF (JOLSEN) THEN
               FAC = EVEC(ICIRED,JCROTX)
            ELSE
               FAC = - EVAL(JCROTX) * EVEC(ICIRED,JCROTX)
            END IF
            CALL DAXPY(NCONF,FAC,WRK(KSCR),1,WRK(JXVECS),1)
            JXVECS = JXVECS + NCONF
  310    CONTINUE
  330 CONTINUE
      IF (JOLSEN) THEN
         JXVECS = KXVECS
         JRVECS = KRVECS
         DO 340 ISIMX = JSIMX+1, JSIMX+MSIMX
            JCROTX = JCROOT(ISIMX)
            FAC = - EVAL(JCROTX)
            CALL DAXPY(NCONF,FAC,WRK(JXVECS),1,WRK(JRVECS),1)
            JXVECS = JXVECS + NCONF
            JRVECS = JRVECS + NCONF
  340    CONTINUE
      END IF
      IF (P6FLAG(43)) THEN
         WRITE (LUPRI,9001) (JCROOT(ISIMX),ISIMX=JSIMX+1,JSIMX+MSIMX)
         CALL OUTPUT(WRK(KRVECS),1,NCONF,1,MSIMX,NCONF,MSIMX,-1,LUPRI)
      END IF
 9001 FORMAT(/' (CINEX) residual vectors for roots',(10I4))
C
C     **********
C     ********** Convergence check, contract RVCS
C     **********
C
C
      IF (P4FLAG(6)) THEN
         IF (WCRESD) THEN
            WRITE (LUW4,9012)
         ELSE
            WRITE (LUW4,9010)
         END IF
      END IF
      JRVECS  = KRVECS
      DO 1100 ISIMX = JSIMX+1,JSIMX+MSIMX
         JCROTX = JCROOT(ISIMX)
         IF (JCROTX .EQ. 0) GO TO 1100
         ECROOT = EMY + EVAL(JCROTX)
         IF (JCROTX .LE. LSTJCR) THEN
C           this is one of the roots we want to converge
            QCNRM = DNORM2(NCONF,WRK(JRVECS),1)
C
            IF (WCRESD) THEN
C
C              Weighted residual gives more uniform vector
C              (wave function) accuracy, however, unweighted
C              is equivalent to CI part of MC gradient in
C              SIRIUS.SIROPT MCSCF optimization
C
               QCFAC = MAX( DTEST, ABS(EVAL(JCROTX)) )
               QCNRM = QCNRM / QCFAC
            END IF
C
C           save residual for final print
C
            RESID(JCROTX) = QCNRM
C
            IF (DOCISRDFT .OR. ISTACI .LE. 0) THEN
               THRQC = THRCIT
            ELSE
               IF (JCROTX .EQ. ISTACI) THEN
                  THRQC = THRCIT
               ELSE
                  THRQC = THRAUX
               END IF
            END IF
            IF (QCNRM .LE. THRQC) THEN
               IF (JCROTX .EQ. LSTJCR) THEN
                  DO 1095 IJCR = NJCR+1,MJCR
                     JCROOT(IJCR) = 0
 1095             CONTINUE
               END IF
               JCROOT(ISIMX) = 0
               IF (P4FLAG(6)) THEN
                  WRITE (LUW4,9030) JCROTX,ECROOT,QCNRM
               ELSE IF (IPRI4 .GT. 5) THEN
                  IF (WCRESD) THEN
                     WRITE (LUW4,9012)
                  ELSE
                     WRITE (LUW4,9010)
                  END IF
                  WRITE (LUW4,9030) JCROTX,ECROOT,QCNRM
               END IF
               IF (JCROTX .EQ. ISTACI) THEN
                  ISTCNV = 1
                  WRITE (LUW4,9014)
               END IF
            ELSE
               IF (P4FLAG(6)) WRITE(LUW4,9020) JCROTX,ECROOT,QCNRM
            END IF
         END IF
         JRVECS = JRVECS + NCONF
 1100 CONTINUE
C
 9010 FORMAT(/' CI root     eigenvalue        residual norm'/)
 9012 FORMAT(/' CI root     eigenvalue     weighted residual norm'/)
 9014 FORMAT(/' The requested root number is now converged.')
 9020 FORMAT(I5,2F20.10)
 9030 FORMAT(I5,2F20.10,' converged')
C
C     set all to converged if ISTACI converged
C
      IF (ISTCNV .NE. 0) THEN
         DO 1120 IJCR = 1,NJCR
            JCROOT(IJCR) = 0
 1120    CONTINUE
      END IF
C
      KCVECS = KRVECS
      JRVECS = KRVECS
      JCVECS = KCVECS
      JXVECS = KXVECS
      JYVECS = KXVECS
      NCVEC  = 0
      DO 1200 ISIMX = JSIMX+1,JSIMX+MSIMX
         JCROTX = JCROOT(ISIMX)
         IF (JCROTX .NE. 0) THEN
            NCVEC = NCVEC + 1
            KCROOT(NCVEC) = JCROTX
            IF (JRVECS.NE.JCVECS) THEN
               CALL DCOPY(NCONF,WRK(JRVECS),1,WRK(JCVECS),1)
               IF (JOLSEN) THEN
                  CALL DCOPY(NCONF,WRK(JXVECS),1,WRK(JYVECS),1)
               END IF
            END IF
            JCVECS = JCVECS + NCONF
            JYVECS = JYVECS + NCONF
         END IF
         JRVECS = JRVECS + NCONF
         JXVECS = JXVECS + NCONF
 1200 CONTINUE
      KCPREV = KCVECS + NCVEC*NCONF
C     ... KCPREV is used in CIORT, it is OK that it overlaps with KXVECS
      KWPHP  = KXVECS + NCVEC*NCONF
      LWPHP  = LFREE + 1 - KWPHP
C
C
C     **********
C     ********** Use Davidson's shift to get CI trial vectors
C     **********
C
C     D = ( CIDIAG(I) - EVAL ) + ESHIFT
C
      IF (ESHIFT.NE.D0 .AND. IPRI4.GT.5) WRITE (LUW4,9120) ESHIFT
 9120 FORMAT(/' (CINEX) Davidson level shift shifted by ESHIFT =',
     *        F20.15)
C
C
      ICRV = KCVECS  - 1
      DO 2400 ICVEC = 1,NCVEC
         KCROTX = KCROOT(ICVEC)
#if defined (VAR_NOPHP)
         DSHIFT = ESHIFT - EVAL(KCROTX) - EMY
         DO 2200 ICONF = 1,NCONF
            D = CDIAG(ICONF) + DSHIFT
            IF (D.LT.DTEST) THEN
               IF (D.GT.(-DTEST)) D = SIGN(DTEST,D)
            END IF
            WRK(ICRV+ICONF) = WRK(ICRV+ICONF) / D
 2200    CONTINUE
#else
         DSHIFT = EMY + EVAL(KCROTX) - ESHIFT
         IF (JOLSEN) THEN
            JXVEC = KXVECS + (ICVEC-1)*NCONF
         ELSE
            JXVEC = 999 999 999
         END IF
         IPRPHP = MAX(0,IPRI6-5)
         CALL NEXCI(JOLSEN,DSHIFT,NCONF,WRK(ICRV+1),WRK(JXVEC),
     *              WRK(ICRV+1),CDIAG,IPRPHP,WRK(KWPHP),LWPHP)
C        CALL NEXCI(OLSEN,ENER,NCVAR,D,XVEC,
C    &              RES,DIAG,IPRPHP,WRK,LWRK)
#endif
         IF (CTRUNC) THEN
C           ... increase efficiency by zeroing elements in
C               trial vector abs less than 1.d-3 abs largest element.
            FACTRN = 1.D-3
            ICMAX  = IDAMAX(NCONF,WRK(ICRV+1),1)
            CTEST  = FACTRN * ABS(WRK(ICRV+ICMAX))
            DO 2300 ICONF = ICRV+1,ICRV+NCONF
               IF (ABS(WRK(ICONF)) .LT. CTEST) WRK(ICONF) = D0
 2300       CONTINUE
         END IF
         IF (SYMTRZ) THEN
C           ... symmetrize CI vector so that numerical
C               inaccuracy won't break symmetry;
C               WARNING : if symmetry were broken in residual vectors,
C                         and thus in sigma vectors, then this will
C                         lead to constant residual.
C
            CNORM = DNORM2(NCONF,WRK(ICRV+1),1)
            TLINDP = SQRT(NCONF*THRLDP)
            IF (CNORM .GT. TLINDP) THEN
               CNORM = D1 / CNORM
               CALL DSCAL(NCONF,CNORM,WRK(ICRV+1),1)
               CALL VSYMTR(NCONF,WRK(ICRV+1),THRSIM)
            END IF
         END IF
         ICRV = ICRV + NCONF
 2400 CONTINUE
C
C     **********
C     **********
C     **********
C
C
C     **********
C     ********** Orthogonalize these NCVEC new b-vectors
C     ********** against previous b-vectors and each other.
C     **********
C
C     CIORT updates NCIT3.
C     LUIT3 is positioned, ready for the new trial vectors.
C
C
      THRLDV = NCONF * THRLDP
      CALL CIORT (NCVEC,NCIT3,LUIT3,KCROOT,WRK(KCVECS),THRLDV,
     *            WRK(KCPREV))
C     CALL CIORT (NCNEW,NCPREV,LUCVEC,KCROOT,CVEC,THRLDP,CPREV)
C
C     **********
C     ********** Check if finished all NJCR residuals,
C     ********** if not go back to 10 for next batch.
C     **********
C
      NCSIM = NCSIM + NCVEC
      NSIML = NSIML - NSIMX
      JSIMX = JSIMX + NSIMX
      IF (NSIML.GT.0) GO TO 10
C     ^-----------------------
C
C
C     **********
C     ********** Test if "in memory"
C     **********
C
C NECgh971126  Fopp has problems with this kind of IF
C NECgh971126      IF (IBATCH .EQ. 1) THEN
C NECgh971126      ELSE
      IF (IBATCH .NE. 1) THEN
         NCSIM = -NCSIM
      END IF
C
C     **********
C     ********** Update NJCR and JCROOT by removing converged vectors.
C     **********
C
      NJ = 0
      DO 8000 ISIMX = 1,NJCR
         IF (JCROOT(ISIMX).GT.0 .AND. JCROOT(ISIMX).LE.LSTJCR) THEN
            NJ = NJ + 1
            JCROOT(NJ) = JCROOT(ISIMX)
         END IF
 8000 CONTINUE
      NJCR  = NJ
      JCONV = NJCR
C
C     **********
C     ********** End of subroutine CINEX
C     **********
C
      CALL QEXIT('CINEX ')
      RETURN
      END
C  /* Deck ciort */
      SUBROUTINE CIORT (NCNEW,NCPREV,LUCVEC,KCROOT,
     *                  CVEC,THRLDP,CPREV)
C
C Written 23-Jun-1985 by Hans Jorgen Aa. Jensen
C
C Purpose:
C  Orthogonalize new CI vectors against all previous CI vectors
C  and among themselves, and renormalize.
C  (Orthogonalization is performed twice if round-off is large,
C   if larger than THRRND).
C
C Input:
C  NCNEW,     number of new CVEC
C  NCPREV,    number of previous, orthonormal vectors on LUCVEC
C  LUCVEC,    file unit with previous vectors
C  CVEC(*,*), non-orthogonal configurational vectors
C  KCROOT(*), index of input CVEC
C             (if KCROOT(i) = 0 then vector no. i is skipped)
C  THRLDP,    threshold for linear dependence
C
C Output:
C  NCNEW,     number of linear indep. vectors in CVEC
C  NCPREV,    number of vectors now on LUCVEC
C  KCROOT(*), index of output CVEC
C  CVEC,      NCNEW orthonormal vectors
C
C Scratch:
C  CPREV(*),    scratch array for old vectors on LUCVEC
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION KCROOT(*), CVEC(NCONF,*), CPREV(*)
C
      PARAMETER ( THRRND=1.D-2, THRTT=1.D-4, D1=1.0D0 )
C
C Used from common blocks:
C   INFVAR : NCONF
C
#include <maxorb.h>
#include <priunit.h>
#include <infpri.h>
C
C
      CALL QENTER('CIORT ')
      IF (NCNEW .EQ. 0) GO TO 9999
C
      NCONF4 = MAX(4,NCONF)
      TLINDP = SQRT(THRLDP)
C
C Normalize input CVEC
C
      DO 100 I = 1,NCNEW
         TT = DNORM2(NCONF,CVEC(1,I),1)
         IF (TT.LE.TLINDP) THEN
            WRITE (LUPRI,3240) I,KCROOT(I),TT
            KCROOT(I) = 0
         ELSE
            IF (TT.LT.THRTT) THEN
               TT = D1 / TT
               CALL DSCAL(NCONF,TT,CVEC(1,I),1)
               TT = DNORM2(NCONF,CVEC(1,I),1)
            END IF
            TT = D1 / TT
            CALL DSCAL(NCONF,TT,CVEC(1,I),1)
         END IF
  100 CONTINUE
C
      ITURN=0
 1500 ITURN=ITURN+1
      IROUND=0
C
C Orthogonalize new vectors agains previous vectors
C
      REWIND LUCVEC
      DO 2000 I=1,NCPREV
         CALL READT(LUCVEC,NCONF,CPREV)
         DO 1910 J = 1,NCNEW
            IF (KCROOT(J) .NE. 0) THEN
               TT = -DDOT(NCONF,CPREV,1,CVEC(1,J),1)
               CALL DAXPY(NCONF,TT,CPREV,1,CVEC(1,J),1)
            END IF
 1910    CONTINUE
 2000 CONTINUE
C
C Orthogonalize new vectors against each other and normalization
C
      DO 3200 I = 1,NCNEW
      IF (KCROOT(I) .EQ. 0) GO TO 3200
C        Orthogonalize against previous new vectors
C        (which have been normalized in the code following 2300)
         DO 2300 J = 1,(I-1)
            IF (KCROOT(J) .NE. 0) THEN
               TT = -DDOT(NCONF,CVEC(1,J),1,CVEC(1,I),1)
               CALL DAXPY(NCONF,TT,CVEC(1,J),1,CVEC(1,I),1)
            END IF
 2300    CONTINUE
C        Normalize
         TT = DNORM2(NCONF,CVEC(1,I),1)
         IF (TT .LE. TLINDP) THEN
            WRITE (LUPRI,3250) I,KCROOT(I),TT
            KCROOT(I) = 0
         ELSE
            IF (TT .LT. THRRND) IROUND=IROUND+1
            IF (TT.LT.THRTT) THEN
               TT = D1 / TT
               CALL DSCAL(NCONF,TT,CVEC(1,I),1)
               TT = DNORM2(NCONF,CVEC(1,I),1)
            END IF
            TT = D1 / TT
            CALL DSCAL(NCONF,TT,CVEC(1,I),1)
         END IF
 3200 CONTINUE
C
C
      IF (IROUND.GT.0 .AND. ITURN.EQ.1) GO TO 1500
      IF (IROUND.GT.0) CALL QUIT('CIORT: round-off errors, see code.')
 3240 FORMAT(/' (CIORT) trial vector no.',I4,' for root',I3,
     *        ' is removed because of linear dependence;',
     *       /' norm of vector on input',1P,D12.5)
 3250 FORMAT(/' (CIORT) trial vector no.',I4,' for root',I3,
     *        ' is removed because of linear dependence;',
     *       /' norm of vector after Gram-Schmidt''s orthogonalization',
     *        1P,D12.5)
C
C Save new vector together with old ones on LUCVEC
C
      J = 0
      DO 4300 I = 1,NCNEW
         IF (KCROOT(I) .NE. 0) THEN
            J = J + 1
            IF (J .LT. I) THEN
               KCROOT(J) = KCROOT(I)
               CALL DCOPY(NCONF,CVEC(1,I),1,CVEC(1,J),1)
            END IF
            CALL WRITT(LUCVEC,NCONF4,CVEC(1,J))
         END IF
 4300 CONTINUE
      IF (J .LT. NCNEW) THEN
         WRITE (LUPRI,4310) NCNEW,J
         NCNEW = J
      END IF
C
C     update NCPREV
C
      NCPREV = NCPREV + NCNEW
C
 4310 FORMAT(/' (CIORT) NCNEW reduced from',I3,' to',I3)
C
C *** End of subroutine CIORT
C
 9999 CALL QEXIT('CIORT ')
      RETURN
      END
C  /* Deck circhk */
      SUBROUTINE CIRCHK(WRK,LFREE)
C
C 1-Aug-1987 Hans Joergen Aa. Jensen, based on REDCHK from SIRNEO.
C
C Purpose:
C   Calculate full reduced CI-matrix for symmetry check in
C   CI calculation;
C   asymmetry would indicate probable bug in CILIN module.
C
#include <implicit.h>
      DIMENSION WRK(LFREE)
      PARAMETER ( D1 = 1.0D0, DIFTST = 1.0D-8 )
C
C  Used from common blocks:
C   INFVAR : NCONF
C   CICB01 : NCIRED
C   INFTAP : LUIT3, LUIT5
C
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <cicb01.h>
#include <inftap.h>
C
C
      REDH(I,J) = WRK((J-1)*NCIRED+I)
      REDS(I,J) = WRK((KREDS-1)+(J-1)*NCIRED+I)
C
      CALL QENTER('CIRCHK')
      KREDH  = 1
      KREDS  = KREDH  + NCIRED*NCIRED
      KSVEC  = KREDS  + NCIRED*NCIRED
      KBVECS = KSVEC  + NCONF
      LBVECS = LFREE  - KBVECS
      IF (LBVECS .LT. NCIRED*NCONF) THEN
C        insufficient space for the algorithm used
         WRITE (LUPRI,'(//2A,/A/)') ' %%% Insufficient space',
     *      ' for symmetry check of reduced CI matrix.',
     *      '     Nothing done.'
         GO TO 9999
      END IF
C
      JBVECI = KBVECS
      REWIND LUIT3
      DO 100 I = 1,NCIRED
         CALL READT(LUIT3,NCONF,WRK(JBVECI))
         JBVECI = JBVECI + NCONF
  100 CONTINUE
C
C     REDH(*,I) = B(*,*)t S(*,I)
C
      I1 = KREDH
      REWIND LUIT5
      DO 200 I = 1,NCIRED
         CALL READT(LUIT5,NCONF,WRK(KSVEC))
         CALL DGEMM('T','N',NCIRED,1,NCONF,1.D0,
     &              WRK(KBVECS),NCONF,
     &              WRK(KSVEC),NCONF,0.D0,
     &              WRK(I1),NCIRED)
         I1 = I1 + NCIRED
  200 CONTINUE
C
C     Reduced overlap matrix, get b vectors
C
C
C     REDS(I,J) = B(*,I)t B(*,J)
C
      CALL DGEMM('T','N',NCIRED,NCIRED,NCONF,1.D0,
     &           WRK(KBVECS),NCONF,
     &           WRK(KBVECS),NCONF,0.D0,
     &           WRK(KREDS),NCIRED)
C
C     check reduced matrices
C
      WRITE (LUPRI,'(//A/A/A,1P,D10.2,A)')
     *    ' %%% Non-symmetric elements in redH, reduced CI matrix,',
     *    '     and wrong elements in redS, reduced CI overlap',
     *    '     ( threshold is',DIFTST,' ) :'
      NERROR = 0
      DO 340 I = 2,NCIRED
         DO 320 J = 1,(I-1)
            DIFF = ABS( REDH(I,J) - REDH(J,I) )
            IF (DIFF .GT. DIFTST) THEN
               NERROR = NERROR + 1
               WRITE (LUPRI,'(A,2I5,2F20.14)')
     *            ' i, j, redh(i,j), redh(j,i) :',
     *            I,J,REDH(I,J),REDH(J,I)
            END IF
            IF (ABS(REDS(I,J)) .GT. DIFTST .OR.
     *          ABS(REDS(J,I)) .GT. DIFTST) THEN
               NERROR = NERROR + 1
               WRITE (LUPRI,'(A,2I5,2F20.14)')
     *            ' i, j, reds(i,j), reds(j,i) :',
     *            I,J,REDS(I,J),REDS(J,I)
            END IF
  320    CONTINUE
         IF (ABS(REDS(I,I) - D1) .GT. DIFTST) THEN
            NERROR = NERROR + 1
            WRITE (LUPRI,'(A,2I5,2F20.14)')
     *         ' i, i, reds(i,i)            :',
     *         I,I,REDS(I,I)
         END IF
  340 CONTINUE
      IF (NERROR .EQ. 0) WRITE (LUPRI,'(A)') ' None, redH is symmetric.'
C
      WRITE (LUPRI,'(//A)') ' Reduced CI matrix (as bT*s) :'
      CALL OUTPUT(WRK(KREDH),1,NCIRED,1,NCIRED,NCIRED,NCIRED,-1,LUPRI)
      WRITE (LUPRI,'(//A)') ' Reduced overlap (as bT*b) :'
      CALL OUTPUT(WRK(KREDS),1,NCIRED,1,NCIRED,NCIRED,NCIRED,-1,LUPRI)
C
      IF (NERROR .GT. 0) THEN
         CALL QTRACE(LUPRI)
         CALL QUIT('ERROR detected in CIRCHK.')
      END IF
C
 9999 CALL QEXIT('CIRCHK')
      RETURN
      END
C  /* Deck cired */
      SUBROUTINE CIRED (ICTL,NREDCI,REDCI,NCNEW,
     *                  CVECS,SVECS,EVAL,EVEC,WRK,LWRK)
C
C Written 25-Jun-1985 by Hans Jorgen Aa. Jensen
C
C Input:
C  ICTL, flow control
C         =1, extend the reduced CI-matrix
C         =2, find the new eigenvalues and -vectors
C         =3, both
C  REDCI,  the old reduced CI-matrix (dimension: NREDCI-NCNEW)
C  NREDCI, dimension of new reduced CI-matrix
C  NCNEW,  number of new CI-vectors
C  CVECS,  the NCNEW new CI-vectors (is destroyed)
C  SVECS,  the sigma vectors of the NCNEW new CI-vectors
C
C  (The reduced CI-matrix is the projection of the CI-matrix
C   on the basis of CI-vectors.)
C
C Output:
C  REDCI, the new, extended reduced CI-matrix (dimension: NREDCI)
C  EVAL,  eigenvalues of the new reduced CI matrix
C  EVEC,  eigenvectors of the new reduced CI matrix
C
C Scratch:
C  CVECS, is used for CI vectors from unit LUIT3 (section 1)
C         and for the reduced CI matrix (section 2)
C         dimension: NCONF and NREDCI*(NREDCI+1)/2, rsp.
C  WRK, real scratch array for BIG and JBIG in JACO, dim: 2*NREDCI
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION CVECS(*), SVECS(NCONF,*)
      DIMENSION REDCI(*), EVEC(NREDCI,*),EVAL(*), WRK(LWRK)
C
C Used from common blocks:
C   INFVAR : NCONF
C   INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
C
C
      CALL QENTER('CIRED ')
C
C Section 1: extend reduced CI matrix with NCNEW new CI-vectors
C
      IF (MOD(ICTL,2) .NE. 1) GO TO 1000
      IF (NCNEW.LT.1) GO TO 1000
      IF (NREDCI .GE. LIROW) 
     &     CALL QUIT('CIRED error: NREDCI exceeds LIROW')
C     IROW(NREDCI+1) used below
      IREDCI = NREDCI - NCNEW
C
C New CI-vectors (are in CVECS)
C
      K1 = 1
      DO 140 K = 1,NCNEW
         KL = IROW(IREDCI+K) + IREDCI
         DO 120 L=1,K
            REDCI(KL+L) = DDOT(NCONF,CVECS(K1),1,SVECS(1,L),1)
  120    CONTINUE
         K1 = K1 + NCONF
  140 CONTINUE
C
C Old CI-vectors.
C Rewind LUIT3
C
      IF (IREDCI .GT. 0) THEN
         REWIND LUIT3
         LL = 0
         DO 240 J = 1,IREDCI
            CALL READT(LUIT3,NCONF,CVECS)
            LL = LL + 1
            DO 220 K=1,NCNEW
               REDCI(IROW(IREDCI+K) + LL ) =
     *            DDOT(NCONF,CVECS,1,SVECS(1,K),1)
  220       CONTINUE
  240    CONTINUE
      END IF
C
C *********************************************************
C Section 2: find eigenvalues and -vectors of reduced CI matrix
C
 1000 CONTINUE
      IF (MOD(ICTL,4) .LE. 1) GO TO 2000
C
C Diagonalize reduced CI matrix.
C Copy REDCI to CVECS for diagonalization.
C
      LMA = IROW(NREDCI+1)
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KREDCI,LMA   ,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF  ,NREDCI,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF ,NREDCI,WRK,KFREE,LFREE)
      CALL DCOPY(LMA,REDCI,1,WRK(KREDCI),1)
      CALL DUNIT(EVEC,NREDCI)
      CALL JACO(WRK(KREDCI),EVEC,NREDCI,NREDCI,NREDCI,
     &          WRK(KBUF),WRK(KIBUF))
      II = KREDCI - 1
      DO 1100 I = 1,NREDCI
         II = II + I
         EVAL(I) = WRK(II)
 1100 CONTINUE
      CALL MEMREL('CIRED.JACO',WRK,1,1,KFREE,LFREE)
      CALL ORDER (EVEC,EVAL,NREDCI,NREDCI)
C
C *** End of subroutine CIRED
C
 2000 CONTINUE
      CALL QEXIT('CIRED ')
      RETURN
      END
C  /* Deck cisave */
      SUBROUTINE CISAVE(NCROOT,NCIRED,ECI,RESID,EVEC,CMO,INDXCI,
     *                  WRK,LFREE)
C
C 27-Sep-1986 Hans Joergen Aa. Jensen
C
C revised 870627-hjaaj, write ci energies to LUIT1 for OCE program
C         890925-hjaaj, write ci residuals to LUIT1 also
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION ECI(*), RESID(*), EVEC(*), CMO(*), INDXCI(*), WRK(*)
#include <dummy.h>
C
C Used from common blocks :
C   INFINP : NROOTS
C   INFOPT : NREDL
C   INFTAP : LUIT1
C
#include <maxorb.h>
#include <infinp.h>
#include <infopt.h>
#include <inftap.h>
#include <infpri.h>
C
      CHARACTER*8 LABENR(4), LABEOD(4)
      DATA LABENR/'********','********','********','ENERGIES'/
      DATA LABEOD/'********','********','********','EODATA  '/
C
      CALL QENTER('CISAVE')
C
      NREDL  = NCIRED
      ISTSAV = ISTATE
      NRTSAV = NROOTS
      ISTATE = ISTACI
      NROOTS = NCROOT
      CALL SIRSAV('CISAVE',CMO,DUMMY,DUMMY,EVEC,DUMMY,INDXCI,WRK,LFREE)
C     CALL SIRSAV(KEYWRD,CMO,IBNDX,REDL,EVEC,XKAP,INDXCI,WRK,LFREE)
      ISTATE = ISTSAV
      NROOTS = NRTSAV
C
C     save ECI on LUIT1
C
      REWIND LUIT1
      CALL MOLLAB('NEWORB  ',LUIT1,LUPRI)
      READ (LUIT1)
C     ... skip "new orbitals", and write energies label:
C
      WRITE (LABENR(3),'(I8)') NCROOT
      WRITE (LUIT1) LABENR
      NCR4 = MAX(4,NCROOT)
      WRITE (LUIT1) (ECI(I),I=1,NCR4)
      WRITE (LUIT1) (RESID(I),I=1,NCR4)
      CALL GETDAT(LABEOD(2),LABEOD(3))
      WRITE (LUIT1) LABEOD
C
      CALL QEXIT('CISAVE')
      RETURN
      END
C  /* Deck cist */
      SUBROUTINE CIST (ICI1,NCSIM,CMO,HDIAG,EMY,WRK,LFREE)
C
C Written 20-Jun-1985 by Hans Joergen Aa. Jensen
C
C Purpose:
C   Obtain initial CI trial vectors for a CI calculation.
C   To write molecular orbitals on LUIT1.
C   The routine is controlled by the ICI1 parameter
C   (similar to the ICI0 parameter for OPTST).
C
C ICI1<0 : select -ICI1 as start configuration
C
C ICI1=1 : Compute startguess of CI-vector from H-diagonal (HDIAG)
C
C ICI1=2 : Start CI-vector(s) are already on LUIT1 (this is checked).
C
C ICI1=4 : Same as ICI1=2, for compatibility with MCSCF CISTRT
C
#include <implicit.h>
C
      DIMENSION CMO(*),HDIAG(*),WRK(*)
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFINP : FLAG(),ISTACI,?
C   INFVAR : NCONF,?
C   INFORB : NCMOT,?
C   INFIND : ?
C   INFDIM : ?
C   INFTAP : LUIT1,?
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
C
C *** Externals:
C
      LOGICAL FNDLAB
C
C *** Local variables:
C
      CHARACTER*8 TABLE(5), LAB123(3), CHRDUM
C
C *** data:
C
      DATA LAB123(1)/'********'/
      DATA TABLE/'********','OLDORB  ','STARTVEC','CIDIAG2 ',
     *           'NEWORB  '/
C
      CALL QENTER('CIST  ')
      CALL GETDAT(LAB123(2),LAB123(3))
      LAB123(3) = '( CIST )'
c     WRITE(LUW4,*) 'CIST called with ICI1 = ',ICI1
      KFREE  = 1
      NC4    = MAX(4,NCONF)
      MCSIM  = 0
      NCSIM0 = NCSIM
      JCI1   = ICI1
      IF (JCI1 .EQ. 4 .OR. JCI1 .EQ. 5) JCI1 = 2
   10 CONTINUE
      IF (JCI1 .NE. 2 .AND. MCSIM .EQ. 0) THEN
C
         CALL NEWIT1
         WRITE (LUIT1) LAB123,TABLE(2)
         CALL WRITT(LUIT1,NCMOT,CMO)
C
C        If RHF or NCONF.eq.1, simply write C0(1) = D1 to LUIT1
C
         IF (FLAG(21) .OR. NCONF.EQ.1) THEN
            WRK(1) = D1
            WRITE (LAB123(2),'(2I4)') 1,1
            WRITE (LUIT1) LAB123,TABLE(3)
            CALL WRITT(LUIT1,4,WRK)
            GO TO 9000
         END IF
      END IF
C
C *** Go to appropriate section to obtain start guess for
C     CI vectors (as specified with the input parameter JCI1)
C
C *** JCI1 < 0 ***
C === start guess = CSF no. -JCI1
C
      IF (JCI1.LT.0) THEN
         IF (-JCI1.GT.NCONF) THEN
            WRITE (LUPRI,5010) -JCI1,NCONF
            CALL QTRACE(LUPRI)
            CALL QUIT('*** ERROR (CIST) illegal JCI1 parameter')
         END IF
         IF (NCSIM.GT.1) THEN
            WRITE (LUPRI,5020)
            CALL QTRACE(LUPRI)
            CALL QUIT('*** ERROR (CIST) illegal JCI1 parameter')
         END IF
         CALL DZERO(WRK,NCONF)
         WRK(-JCI1) = D1
         WRITE (LAB123(2),'(2I4)') 1,-1
         WRITE (LUIT1) LAB123,TABLE(3)
         CALL WRITT(LUIT1,NC4,WRK)
         GO TO 9000
      END IF
 5010 FORMAT(/' (CIST) ERROR, specified start configuration no.',
     *  I10,' is non-existent (last CSF is no.',I10,')')
 5020 FORMAT(/' (CIST) ERROR, negative JCI1 option not allowed ',
     *  'when NCSIM is greater then one.')
C
      IF (JCI1.LT.1 .OR. JCI1.GT.2) THEN
         WRITE (LUPRI,5030) JCI1
         CALL QTRACE(LUPRI)
         CALL QUIT('*** ERROR (CIST) illegal JCI1 parameter')
      END IF
 5030 FORMAT(/' CIST FATAL ERROR, JCI1 =',I4,' is not defined')
C
C
C
C
#if !defined (VAR_OLDCIST)
      NCSIM = NCSIM0 + MCSIM
C     ... (870804-hjaaj) Until now we used IFRST = MCSIM + 1
C         and no change in NCSIM, however, it gave bad convergence
C         in roots MCSIM+1 to NCSIM not to include the lowest
C         MCSIM diagonal elements explicitly.
C     940511-hjaaj MAERKE: Maybe this can be done with .OLSEN ??
#endif
      GO TO (100,200) JCI1
C
C *** JCI1 = 1 ***
C *** Use PHPGET and PHPVEC to get lowest eigenvectors of subspace
C *** or find smallest diag. elements of H to select start config.
C
  100 CONTINUE
      IF (MAXPHP .GE. 2*NCSIM) THEN
         IF (MCSIM .EQ. 0) THEN
            WRITE (LAB123(2),'(2I4)') NCSIM,-1
            LAB123(3) = 'CIST-PHP'
            WRITE (LUIT1) LAB123,TABLE(3)
         ELSE
C           ... some, but not all, requested CI trial vectors available
            REWIND LUIT1
            CALL MOLLAB(TABLE(3),LUIT1,LUPRI)
            DO 161 I = 1,MCSIM
               READ (LUIT1)
  161       CONTINUE
         END IF
         DO 110 ICSIM = 1,NCSIM
            CALL PHPVEC(1,ICSIM,WRK,NCONF,HDIAG)
C           CALL PHPVEC(NVEC,I1VEC,CVECS,NCONF,DIAG)
            CALL WRITT(LUIT1,NC4,WRK)
  110    CONTINUE
         NCSIM = MCSIM + NCSIM
      ELSE
         CALL CIST1(ISTACI,FLAG(58),NCONF,NCSIM,MCSIM,HDIAG,EMY,
     &              WRK,LFREE)
C        CALL CIST1(ISTACI,PLUSCB,NCONF,NCSIM,MCSIM,HDIAG,EMY,WRK,LWRK)
      END IF
      GO TO 9000
C
C *** JCI1 = 2 ***
C === start guess should already be on LUIT1
C     (but we test to be sure!)
C
  200 CONTINUE
         MCSIM = 0
         REWIND LUIT1
         IF ( .NOT.FNDLAB(TABLE(2),LUIT1) ) GO TO 290
            READ (LUIT1,END=290,ERR=290)
         IF ( .NOT.FNDLAB(TABLE(3),LUIT1) ) GO TO 290
         DO 210 I = 1,NCSIM
            READ (LUIT1,END=290,ERR=290) CHRDUM
            IF (CHRDUM .EQ. TABLE(1)) GO TO 290
            MCSIM = I
  210    CONTINUE
         REWIND LUIT1
      GO TO 9000
C
  290 CONTINUE
         REWIND LUIT1
         IF (MCSIM .EQ. 0) THEN
            WRITE (LUW4,'(//A//)')
     *      ' *** CI control: start vectors not found on LUIT1,'//
     *      ' lowest diagonal elements used.'
         ELSE
#if defined (VAR_OLDCIST)
            WRITE (LUW4,'(//A,I4,A,/2A//)')
     *      ' *** CI control: only',MCSIM,
     *      ' start vectors found on LUIT1,',
     *      ' diagonal elements used for remaining',
     *      ' requested trial vectors.'
         END IF
#else
            WRITE (LUW4,'(//A,I4,A,/A/A//)')
     *      ' *** CI control: only',MCSIM,
     *      ' start vectors found on LUIT1,',
     *      ' In addition to these start vectors, lowest diagonal',
     *      ' elements used for requested trial vectors.'
#endif
         END IF
         JCI1 = 1
      GO TO 10
C
C *** End of subroutine CIST
C
C     Write NEWORB at end of file (so e.g. SIRORB won't overwrite
C     STARTVEC).
C
 9000 CONTINUE
      IF (JCI1 .NE. 2) THEN
         CALL NEWORB('CIST    ',CMO,.FALSE.)
      END IF
      CALL QEXIT('CIST  ')
      RETURN
      END
C  /* Deck cist1 */
      SUBROUTINE CIST1(ISTACI,PLUSCB,NCONF,NCSIM,MCSIM,
     *                 HDIAG,EMY,WRK,LWRK)
C
C  30-Oct-1990 hjaaj
C
C  Block of old CIST, purpose: find start vectors corresponding
C  to lowest diagonal elements and write these to LUIT1.
C
#include <implicit.h>
      LOGICAL PLUSCB
      DIMENSION HDIAG(*), WRK(*)
      PARAMETER (MXNLIS = 220)
      PARAMETER (THRSIM = 0.000 01 D0, THRDEG = 0.000 000 001 D0)
C     ... THRSIM is threshold for two values are similar
C         THRDEG is threshold for two values are degenerate
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFTAP : LUIT1
C   INFPRI : LUW*, IPRI*
C
#include <infpri.h>
#include <inftap.h>
#include <priunit.h>
C
C *** Local variables:
C
      DIMENSION NLIST(MXNLIS)
      CHARACTER*8 LBLIT1(4)
C
C *** data:
C
      DATA LBLIT1/'********','        ','        ','STARTVEC'/
C
      CALL QENTER('CIST1 ')
C
      NC4    = MAX(4,NCONF)
C
      IF (PLUSCB) THEN
C        ... take plus combination of degenerate diagonal elements
C            we thus need more ...
         NRS = MIN( NCONF, 2 * NCSIM + 6)
         IF (NRS .GT. MXNLIS .AND. MXNLIS .GT. NCSIM) NRS = MXNLIS
         WRITE (LUPRI,'(//A/A/)')
     *   ' >>> SIRCI.CIST1: plus combination of all degenerate',
     *   '                 configurations is used as start vectors.'
      ELSE
         NRS = MIN(NCSIM + 3,NCONF)
C        ... we add 3 so we can expand if degeneracy.
      END IF
      IF (NRS .GT. MXNLIS) THEN
         WRITE (LUPRI,'(//A,I4,A/I9,A/)')
     1      ' >>> SIRCI.CIST1 ERROR, LOCAL DIMENSION ONLY ALLOWS',
     2      MXNLIS,' DIAGONAL ELEMENTS,', NRS,' WERE REQUESTED.'
         CALL QTRACE(LUPRI)
         CALL QUIT('>>> SIRCI.CIST1 ERROR, MXNLIS PARAMETER TOO SMALL.')
      END IF
      IF (LWRK .LT. NRS+NCONF) CALL ERRWRK('CIST1',(NRS+NCONF),LWRK)
      CALL CIFNMN(NRS,NLIST,HDIAG,NCONF,WRK,LWRK)
C     CALL CIFNMN(NELMNT,IPLACE,VEC,NVEC,WRK,LWRK)
      IF (IPRI4 .GT. 1) WRITE (LUW4,1410)
     *   NRS,(I,NLIST(I),WRK(I)-EMY,WRK(I),I=1,NRS)
      IF (LUPRI .NE. LUW4 .AND. IPRI4 .GT. 1) WRITE (LUPRI,1410)
     *   NRS,(I,NLIST(I),WRK(I)-EMY,WRK(I),I=1,NRS)
 1410 FORMAT(/' (CIST1)',I3,' lowest diagonal elements:',
     *      //' Element no. Config.no.    Active energy',
     *        '      Total energy',
     *      //,(I10,' : ',I10,2F18.10))
C
C
C
#if defined (VAR_OLDCIST)

... Gave bad convergence, when beginning with 20 roots from, e.g.,
    10 converged roots (i.e. MCSIM = 10)
    (870804-hjaaj)

      IF (PLUSCB .AND. MCSIM .GT. 0) THEN
C        ... if use plus combination of degenerate configurations
C            and extension of CI space, count were to start:
         ICSIM = 0
         IFRST = 1
  142       IFRST = IFRST + 1
            IF ( (WRK(IFRST)-WRK(IFRST-1)) .LT. THRSIM ) GO TO 142
            ICSIM = ICSIM + 1
            IF (ICSIM .LT. MCSIM) GO TO 142
         WRITE (LUPRI,'(/A,I5)')
     *   ' New trial vectors selection starts at element no.',IFRST
      ELSE IF (.NOT.PLUSCB) THEN
C        ... else include the degenerate configurations separately,
C            check if last element is degenerate with next element(s)
         IFRST = MCSIM + 1
         LAST  = NCSIM
         DO 150 I = NCSIM+1,NRS
            IF ((WRK(I)-WRK(NCSIM)).LT.THRSIM) LAST = I
  150    CONTINUE
         IF (LAST .GT. NCSIM) THEN
            WRITE (LUPRI,1500) NCSIM,LAST
            NCSIM = LAST
         END IF
 1500    FORMAT(//' (CIST1) number of start vectors is increased from',
     *          I3,' to',I3,' because of degeneracy')
      ELSE
         IFRST = 1
      END IF
#endif
#if !defined (VAR_OLDCIST)
      IFRST = 1
      IF (.NOT.PLUSCB) THEN
C        ... include the degenerate configurations separately,
C            check if last element is degenerate with next element(s)
         LAST  = NCSIM
         DO 150 I = NCSIM+1,NRS
            IF ((WRK(I)-WRK(NCSIM)).LT.THRSIM) LAST = I
  150    CONTINUE
         IF (LAST .GT. NCSIM) THEN
            WRITE (LUPRI,1500) NCSIM,LAST
            NCSIM = LAST
         END IF
 1500    FORMAT(//' (CIST1) number of start vectors is increased from',
     *          I3,' to',I3,' because of degeneracy')
      END IF
#endif
C
C === Write start guess on LUIT1
C
      LAST  = NCSIM
      IF (MCSIM .EQ. 0) THEN
         WRITE(LBLIT1(2),'(2I4)') NCSIM,ISTACI
         LBLIT1(3) = '(CIST1 )'
         WRITE (LUIT1) LBLIT1
      ELSE
C        ... some, but not all, requested CI trial vectors available
         REWIND LUIT1
         CALL MOLLAB(LBLIT1(4),LUIT1,LUPRI)
         DO 161 I = 1,MCSIM
            READ (LUIT1)
  161    CONTINUE
#if defined (VAR_OLDCIST)
C        ... avoid getting part of a degenerate eigenvalue
  162    IF ((WRK(IFRST)-WRK(IFRST-1)).LT.THRSIM) THEN
            IFRST = IFRST + 1
            LAST  = MIN(LAST+1,NRS)
            IF (IFRST .LT. LAST) GO TO 162
         END IF
#endif
      END IF
      CALL DZERO(WRK(NRS+1),NCONF)
      IF (PLUSCB) THEN
C        ... take plus combination of degenerate diagonal elements
  170    CONTINUE
C           ... 170: repeat until (MCSIM .GE. NCSIM .OR. NRS .GE. LAST)
            LAST  = IFRST
  172       IF ( (WRK(LAST+1)-WRK(IFRST)) .LT. THRDEG ) THEN
               LAST = LAST + 1
               IF (LAST .LT. NRS) GO TO 172
            END IF
            FAC = 1 + LAST - IFRST
            FAC = D1 / SQRT(FAC)
            DO 174 I = IFRST,LAST
               WRK(NRS + NLIST(I)) = FAC
  174       CONTINUE
            CALL WRITT(LUIT1,NC4,WRK(NRS+1))
            MCSIM = MCSIM + 1
            DO 176 I = IFRST,LAST
               WRK(NRS + NLIST(I)) = D0
  176       CONTINUE
            IFRST = LAST + 1
         IF (MCSIM .LT. NCSIM .AND. LAST .LT. NRS) GOTO 170
C           ... end repeat
      ELSE
         DO 185 I = IFRST,LAST
            MCSIM = MCSIM + 1
            WRK(NRS + NLIST(I)) = D1
            CALL WRITT(LUIT1,NC4,WRK(NRS+1))
            WRK(NRS + NLIST(I)) = D0
  185    CONTINUE
      END IF
      IF (NCSIM .GT. MCSIM .AND. LAST .LT. NCONF) THEN
C        ... MXNLIS too small, or NRS too small compared to NCSIM.
         WRITE (LUPRI,1850) MCSIM, NCSIM
 1850    FORMAT(
     *   //' (CIST1) INTERNAL ERROR, we only got',I4,' start vectors',
     *    /'                      ',I4,' were requested',
     *    /' >>> action: revise code in SIRCI.CIST1 or modify problem.')
         CALL QUIT('CIST1 internal error, found too few start vectors.')
      END IF
      NCSIM = MCSIM
      CALL QEXIT('CIST1 ')
      RETURN
      END
C  /* Deck cifnmn */
      SUBROUTINE CIFNMN(NELMNT,IPLACE,VEC,NVEC,WRK,LWRK)
C
C 12-Aug-1990 hjaaj, based on CIST
C (Find minimum elemnts)
C
C Purpose:
C   Return in IPLACE addresses on lowest NELMNT elements in VEC.
C
#include <implicit.h>
      DIMENSION VEC(NVEC), IPLACE(NELMNT), WRK(LWRK)
C
      IF (LWRK .LT. NELMNT) THEN
         CALL QUIT('CIFNMN: Insufficient memory (LWRK .lt. NELMNT)')
      END IF
      IF (NELMNT .GT. NVEC) THEN
         CALL QUIT('CIFNMN ERROR: NELMNT .gt. NVEC')
      END IF
C
C     Sort first NELMNT elements of VEC
C
      WRK(1) = VEC(1)
      IPLACE(1) = 1
      DO 120 I = 2,NELMNT
         VECI = VEC(I)
         DO 115 J = 1,(I-1)
         IF (WRK(J) .GT. VECI) THEN
            DO 112 K = I,(J+1),-1
               WRK(K)    = WRK(K-1)
               IPLACE(K) = IPLACE(K-1)
  112       CONTINUE
            IPLACE(J) = I
            WRK(J)    = VECI
         GO TO 120
         END IF
  115    CONTINUE
         IPLACE(I) = I
         WRK(I)    = VECI
  120 CONTINUE
C
C     Find lowest elements by insertion sort
C
      XHGH = WRK(NELMNT)
      DO 140 I = NELMNT+1,NVEC
      IF (VEC(I).GE.XHGH) GO TO 140
         DO 130 J = 1,NELMNT
         IF (VEC(I) .LT. WRK(J)) THEN
            DO 135 K = NELMNT,(J+1),-1
               WRK(K) = WRK(K-1)
               IPLACE(K) = IPLACE(K-1)
  135       CONTINUE
            IPLACE(J) = I
            WRK(J) = VEC(I)
            XHGH   = WRK(NELMNT)
            GO TO 140
         END IF
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
C  /* Deck cirsps */
      SUBROUTINE CIRSPS(ECIREF,EMY,CMO,FCAC,H2AC,HDIAG,INDXCI,WRK,LFREE)
C
C  17-Jun-1987 hjaaj
C  Revised 13-Jun-1991 hjaaj : also calculate and write DV
C
C Purpose:
C  Interface to RESPONSE for CI response calculations:
C  write information to LUSIFC needed for response calculation.
C
C  The following records are written:
C
C    1) POTNUC,EMY,EACTIV,ECIREF
C    2) NISHT,NASHT,...
C    3) CMO
C    4) CREF
C    5) DV
C    6) FCAC
C    7) H2AC
C    8) label "CIDIAG2 "
C    9) HDIAG (diagonal of CI matrix)
C   10) label "SIRFLAGS"
C   11) FLAG(1:NFLAG)
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION CMO(*), FCAC(*), H2AC(*), HDIAG(*), INDXCI(*), WRK(*)
#include <iratdef.h>
#include <dummy.h>
C
C Used from common blocks:
C
C   INFINP : ISTACI,ISPIN,NACTEL,FLAG(1:NFLAG)
C   INFORB : NSYM,NISHT,NASHT,NNASHX,NNASHY,NOCCT,NNOCCX,NCMOT,...
C   INFVAR : NCONF
C   INFTAP : LUIT1, LUSIFC
C   CBGETD : IADH2
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpri.h>
#include <cbgetdis.h>
C
      CHARACTER*8 LAB123(3), TABLE(5)
      DATA LAB123/'********','********','********'/
      DATA TABLE /'CIDIAG2 ','CIRESPON','SIRFLAGS','EODATA  ',
     *            'CREFDETS'/
C
      CALL QENTER('CIRSPS')
      CALL GETDAT(LAB123(2),LAB123(3))
C     ... place date in LAB123(2) and time in LAB123(3)
C
C *** Now create LUSIFC and write ...
C
      CALL GPOPEN(LUSIFC,FNSIFC,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
C     Retrieve CREF
C     Calculate DV
C
      KCREF = 1
      KDV   = KCREF + NCONF
      KWRK1 = KDV   + NNASHX
      LWRK1 = LFREE + 1 - KWRK1
C
      REWIND LUIT1
      CALL MOLLAB('STARTVEC ',LUIT1,LUPRI)
      DO 100 I = 1, (ISTACI-1)
         READ (LUIT1)
  100 CONTINUE
      CALL READT(LUIT1,NCONF,WRK(KCREF))
C
      CALL MAKDV(WRK(KCREF),WRK(KDV),INDXCI,WRK(KWRK1),LWRK1)
C
C     Write LUSIFC
C
      WRITE (LUSIFC) LAB123,TABLE(2)
      ECACTV = ECIREF - EMY
      WRITE (LUSIFC) POTNUC,EMY,ECACTV,ECIREF,ISTACI,ISPIN,NACTEL,LSYM
      WRITE (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,
     *               NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT,
     &               NSYM, MULD2H, NRHF, NFRO,
     &               NISH, NASH, NORB, NBAS,
     &               NELMN1, NELMX1, NELMN3, NELMX3, MCTYPE,
     &               NAS1, NAS2, NAS3
C
      NC4    = MAX(4,NCONF)
      NCMOT4 = MAX(4,NCMOT)
      MMASHX = MAX(4,NNASHX)
C
      CALL WRITT (LUSIFC,NCMOT4,CMO)
      CALL WRITT (LUSIFC,NC4,WRK(KCREF))
      CALL WRITT (LUSIFC,MMASHX,WRK(KDV))
      CALL WRITT (LUSIFC,MMASHX,FCAC)
      IF (IADH2 .GE. 0) THEN
#if defined (VAR_NEWCODE)
... 27-Jun-1987/hjaaj -- noget lignende : MAERKE
         DO 100 IJ = 1,NNASHX
            CALL READDX(LUH2AC,IADH2+IJ,IRAT*NNASHX,H2AC)
            CALL WRITT (LUSIFC,MMASHX,H2AC)
  100    CONTINUE
... men der mangler en label som fortaeller response at H2AC on disk.
#else
         CALL QUIT('SIRCI.CIRSPS: ".DISKH2" not implemented here yet.')
#endif
      ELSE
         CALL WRITT (LUSIFC,(MMASHX*MMASHX),H2AC)
      END IF
      WRITE (LUSIFC) LAB123,TABLE(1)
      CALL WRITT (LUSIFC,NC4,HDIAG)
      WRITE (LUSIFC) LAB123,TABLE(3)
      WRITE (LUSIFC) (FLAG(I),I=1,NFLAG)
C
C     Write CREF in determinants here, if NCDETS .gt. NCONF
C
      IF (NCDETS .GT. NCONF) THEN
         IF (FLAG(27)) THEN
            WRITE (LUPRI,*)
     &      'WR_SIRIFC ERROR, .DETERMINANTS but NCDETS.gt.NCONF:',
     &      NCDETS, NCONF
            CALL QUIT('CIRSPS: NCDETS .gt. NCONF for .DETERMINANTS')
         END IF
         WRITE (LUSIFC) LAB123,TABLE(5)
         KCDET = KDV
         KWRK1 = KCDET + NCDETS
         LWRK1 = LFREE - KWRK1
         CALL GETDETS(LSYM,NCONF,NCDETS,INDXCI,WRK(KCREF),WRK(KCDET),
     &      WRK(KWRK1),LWRK1)
         NC4 = MAX(4,NCDETS)
         CALL WRITT(LUSIFC,NC4,WRK(KCDET))
      END IF
C
C     Write EODATA
C
      WRITE (LUSIFC) LAB123,TABLE(4)
C
C *** end of subroutine CIRSPS
C
      CALL GPCLOSE(LUSIFC,'KEEP')
      CALL QEXIT('CIRSPS')
      RETURN
      END
      SUBROUTINE GETDVREF(DVREF,IST,INDXCI,WRK,LWRK)
C
C     Get DV_ref for CI-DFT and MC-DFT models.
C
C     (c) jkp + hjaaj July 2003
C
C     =====================================================
C     Construct active one-electron density matrix
C     DVREF for CREF vector no. IST
C     =====================================================
C
#include "implicit.h"
C
      DIMENSION DVREF(*), INDXCI(*), WRK(LWRK)
C
C Used from common blocks:
C  infvar: NCONF
C  inftap: LUIT1
C
#include "infvar.h"
#include "inftap.h"
C
C
      IF (IST .LE. 0) THEN
         CALL QUIT('GETDVREF called with IST .le. 0')
      END IF
      KCREF = 1
      KW1   = KCREF + NCONF
      LW1   = LWRK  - KW1
      IF (LW1 .LE. 0) CALL ERRWRK('GETDVREF',-KW1,LWRK)
C
      REWIND LUIT1
      CALL MOLLAB('STARTVEC ',LUIT1,LUPRI)
      DO I = 1, (IST-1)
         READ (LUIT1)
      END DO
      CALL READT(LUIT1,NCONF,WRK(KCREF))
      CALL MAKDV(WRK(KCREF),DVREF,INDXCI,WRK(KW1),LW1)
      RETURN
      END
