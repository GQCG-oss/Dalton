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
C  /* Deck abarsp */
      SUBROUTINE ABARSP(CICLC,HFCLC,TRPCLC,OOTV,IOPSYM,EXCLC,
     *                  EXVAL,NEXVAL,NABATY,NABAOP,LABEL,LUGDVE,
     *                  LUSOVE,LUREVE,THRCLC,MAXCLC,IPRCLC,MXRM,
     *                  MXPHP,WRK,LWRK)
C
C      CICLC  - true for CI calculations
C      HFCLC  - true for RHF - closed shell or one electron in one
C               active orbital
C      TRPCLC - true for triplet perturbation operators
C      OOTV   - true for optimal orbital trial vectors
C      IOPSYM - symmetry of perturbation operators
C      EXCLC  - true for excitation energy calculations,
C               false for linear response equations
C      EXVAL  - calculated excitation energies (output)
C             - frequency for linear response (input)
C      NEXVAL - number of excitation energies/frequencies
C      NABATY - 1 for real operator, -1 for imaginary operator (for each
C               operator)
C      NABAOP - number of right-hand sides
C      LUGDVE - unit number for right-hand sides
C      LUSOVE - unit number for solutions
C      LUREVE - unit number for residuals
C      THRCLC - threshold for convergence
C      MAXCLC - maximum number of iterations
C      IPRCLC - print level
C      MXRM   - maximum size of reduced space
C      MXPHP  - maximum size for explicit subblock of configuration
C               Hessian
C  In common:
C
C      FABRHS - true : transform integrals to MO basis;
C               normally true first time routine is called
C               at a new geometry.
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
#include <maxorb.h>
C
      DIMENSION WRK(LWRK),EXVAL(NEXVAL),NABATY(NABAOP)
C
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      CHARACTER*8 LABEL
C
C Used from common blocks:
C ABAINF : DODRCT
C INFLIN : NCONRF
C INFTRA : USEDRC
C
#include <abainf.h>
#include <exeinf.h>
#include <inflin.h>
#include <infinp.h>
#include <cbisol.h>
#include <inforb.h>
#include <infpri.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infdim.h>
#include <inftra.h>
#include <inftap.h>
#include <rspprp.h>
#include <inflr.h>
#include <infpp.h>
#include <infhyp.h>
#include <priunit.h>
#include <infsop.h>
C
      LOGICAL CICLC,HFCLC,TRPCLC,EXCLC,OOTV,FLAGSV,EX
C
      CALL QENTER('ABARSP')
      TIMRSP = SECOND()
      IPRRSP = IPRCLC - 2
      IF (IPRRSP .GE. 0) THEN
         WRITE (LUPRI,'(//A/A)')
     *   '  THIS IS OUTPUT FROM THE MCSCF AND SOPPA RESPONSE SOLVER',
     *   ' ---------------------------------------------------------'
         IF (EXCLC) THEN
            WRITE(LUPRI,'(2(/A,I6),/A,1P,D13.2)')
     *      ' Symmetry of excitation operator',IOPSYM,
     *      ' Number of excitations:',NEXVAL,
     *      ' Convergence threshold:',THRCLC
         ELSE IF (TRPCLC) THEN
            WRITE(LUPRI,'(/3(/A,I6),/A,1P,D13.2)')
     *      ' Symmetry of triplet property operator', IOPSYM,
     *     ' Number of operators for triplet linear response equations:'
     *      ,NABAOP, ' Number of response frequencies:',NEXVAL,
     *      ' Convergence threshold:',THRCLC
         ELSE
            WRITE(LUPRI,'(/3(/A,I6),/A,1P,D13.2)')
     *      ' Symmetry of singlet property operator', IOPSYM,
     *     ' Number of operators for singlet linear response equations:'
     *      ,NABAOP, ' Number of response frequencies:',NEXVAL,
     *      ' Convergence threshold:',THRCLC
         ENDIF
      END IF
      CALL FLSHFO(LUPRI)
C
      IF (DODRCT .AND. NCONRF.LE.1) FABRHS = .FALSE.
C     ... no integral transformation for HF or DFT direct calculations.
C
      REWIND LUGDVE
      REWIND LUSOVE
      REWIND LUREVE
C
C CONSTRUCTION OF ORBITAL DIAGONAL WITH AVDIA
C
      AVDIA = .TRUE.
C
C     *************************
C     ***** INPUT SECTION *****
C     *************************
C
C DEFINE CALCULATION BY FLAGS AND LOGICALS
C
      IF (TRPCLC) THEN
         TRPLET = .TRUE.
         TRPFLG = .TRUE.
      ELSE
         TRPLET = .FALSE.
         TRPFLG = .FALSE.
      ENDIF
C
      IF (CICLC) THEN
         RSPCI = .TRUE.
      ELSE
         RSPCI = .FALSE.
      ENDIF
C
      IF (CICLC.AND.HFCLC) THEN
         WRITE (LUPRI,'(//,A,/,A,/)')
     *   ' ERROR: THE CALCULATION IS SPECIFIED AS BOTH CI AND HF'
     *   ,' CONSEQUENTLY THE CALCULATION MUST STOP'
         CALL QUIT('ABARSP: INCORRECT SPECIFICATION OF CICLC AND HFCLC')
      ENDIF
C
      KSYMOP = IOPSYM
      OPTORB = OOTV
      MAXRM  = MXRM
      MAXPHP = MXPHP
C
      IF (EXCLC) THEN
         THCPP  = THRCLC
         MAXITP = MAXCLC
         IPRPP  = IPRCLC - 2
         NPPCNV(KSYMOP) = NEXVAL
         NPPSIM(KSYMOP) = NEXVAL
         NPPSTV(KSYMOP) = NEXVAL
      ELSE
         THCLR  = THRCLC
         MAXITL = MAXCLC
         IPRLR  = IPRCLC - 2
         NFREQ  = NEXVAL
         DO 120 I = 1,NEXVAL
            FREQ(I) = EXVAL(I)
  120    CONTINUE
      ENDIF
C
C
      IF (NEXVAL.GT.0) CALL RSPSET
C
      CALL FLSHFO(LUPRI)
C
C     PERFORM INTEGRAL TRANSFORMATION
C
cSPAS:29.11.2002 debugging newtra for SOPPA
cSPAS IF (FABRHS) THEN
      IF (FABRHS .AND. 
     &    .NOT. (SOPPA .AND. NEWTRA .AND. .NOT. CCPPA)) THEN
cKeinSPASMehr:29.11.2002 debugging newtra for SOPPA
         KCMO   = 1
         KWTRA  = KCMO   + NCMOT
         LWTRA  = LWRK   - KWTRA
         IF (LWTRA.LT.0) CALL ERRWRK('ABARSP 1',-(KWTRA-1),LWRK)
         REWIND LUSIFC
         IF (RSPCI) THEN
            CALL MOLLAB('CIRESPON',LUSIFC,LUERR)
         ELSE
            CALL MOLLAB(LBSIFC,LUSIFC,LUERR)
         END IF
         READ (LUSIFC)
         READ (LUSIFC)
         CALL READT (LUSIFC,NCMOT,WRK(KCMO))
         IF ( .NOT. RSPCI ) THEN
            USEDRC = .TRUE.
         ELSE
            USEDRC = .FALSE.
         ENDIF
cs
         IF (SOPPA .OR. NEXVAL.EQ.0) THEN
            JTRLVL = -10
         ELSE IF (NEWTRA) THEN
            JTRLVL = -5 ! level 5 needed for occ-sec distributions in NEWTRA
         ELSE
            JTRLVL = -4
         ENDIF
cs
         FLAGSV = DORSP
         DORSP = .TRUE.
         CALL SIR_INTOPEN
         DORSP = FLAGSV
         CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KWTRA),LWTRA) 
      END IF
      FABRHS = .FALSE.
C
C     ORGANIZE CALCULATION FOR EACH PERTURBATION OPERATOR
C
      IF (NEXVAL.EQ.0) THEN
         KWSYM = 1
         LWSYM = LWRK  - KWSYM
         IF (LWSYM.LT.0) CALL ERRWRK('ABARSP  ',KWSYM,LWRK)
         CALL RSPSYM(WRK(KWSYM),LWSYM)
      END IF
C
C     ALLOCATE WORK SPACE FOR MATRICES THAT WILL BE KEPT DURING THE
C     WHOLE RESPONSE CALCULATION AND READ IN THE MATRICES
C
      KFREE  = 1
      LFREE  = LWRK
      CALL MEMGET('REAL',KINDX,LCINDX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO ,NCMOT ,WRK,KFREE,LFREE)
      LUDV   = NACTT*NACTT
      CALL MEMGET('REAL',KUDV ,LUDV  ,WRK,KFREE,LFREE)
      IF (RSPCI) THEN
         CALL MEMGET('REAL',KPVX ,0     ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFOCK,0     ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFC  ,0     ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFV  ,0     ,WRK,KFREE,LFREE)
      ELSE
         IF (TRPFLG) THEN
C NEED BOTH TRIPLET AND SINGLET TWO ELECTRON DENSITY MATRICES
            LPVX = 2*LPVMAT
         ELSE
C NEED ONLY SINGLET TWO ELECTRON DENSITY MATRIX
            LPVX = LPVMAT
         END IF
         CALL MEMGET('REAL',KPVX ,LPVX  ,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFOCK,N2ORBT,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFC  ,NNORBT,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFV  ,NNORBT,WRK,KFREE,LFREE)
      END IF
      CALL MEMGET('REAL',KFCAC,NNASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2AC,NNASHX*NNASHX,WRK,KFREE,LFREE)
      KTOT  =  KFREE
      KWRK1  = KFREE
      LWRK1  = LFREE
C
      IPRRSP = IPRCLC - 5
C
C     For SOPPA :
C
      IF (SOPPA) THEN
C
C        Initialize XINDX
C
         A2EXIST=.FALSE.
         CALL DZERO(WRK(KINDX),LCINDX)
C
C        Find address array's for SOPPA calculation
C
         CALL SET2SOPPA(WRK(KINDX+KABSAD-1),WRK(KINDX+KABTAD-1),
     *                  WRK(KINDX+KIJSAD-1),WRK(KINDX+KIJTAD-1),
     *                  WRK(KINDX+KIJ1AD-1),WRK(KINDX+KIJ2AD-1),
     *                  WRK(KINDX+KIJ3AD-1),WRK(KINDX+KIADR1-1))
C
      ENDIF
C
C     ... suppress printing of SCF/MCSCF energy for IPRRSP .gt. 0
C         for normal print levels (.le. 5) /960905-hjaaj
      CALL RSPMC(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),WRK(KFC),
     *           WRK(KFV),WRK(KFCAC),WRK(KH2AC),WRK(KINDX),WRK(KWRK1),
     *           LWRK1)
C     CALL RSPMC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,IADR1,WRK,LWRK)
      IPRRSP = IPRCLC - 2
C
      IF (TDHF .NEQV. HFCLC) THEN
         WRITE(LUPRI,'(//A,L10/A,L10)')
     &   ' ABARSP WARNING: "HFCLC" in parameter list is',HFCLC,
     &   '                 but "TDHF" after RSPMC    is',TDHF
         NWARN = NWARN + 1
      END IF
C
      IF ( .NOT.TDHF ) THEN
         CALL GETCIX(WRK(KINDX),IREFSY,IREFSY,WRK(KTOT),LFREE,0)
      END IF
C
      IF ( (.NOT.TDHF) .AND. (.NOT.RSPCI) ) THEN
C
C     ... CALCULATE ONE- AND TWO- BODY DENSITY MATRICES
C         ( THE SYMMETRIC MC TWO BODY DENSITY MATRIX CANNOT BE USED IN
C           RESPONSE CALCULATION )
C
         KCREF  = KWRK1
         KTOT   = KCREF + NCREF
         LFREE  = LWRK  - KTOT
         IF (LFREE.LT.0) CALL ERRWRK('ABARSP 2',-(KTOT-1),LWRK)
         KFREE  = 1
C
         CALL GETREF(WRK(KCREF),NCREF)
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/A)')' ***** ONE BODY DENSITY MATRIX '//
     &                         'from SIRIFC file'
            CALL OUTPUT(WRK(KUDV),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
         ISPIN1 = 0
         ISPIN2 = 0
         CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     *              WRK(KUDV),WRK(KPVX),
     *              ISPIN1,ISPIN2,.FALSE.,.FALSE.,
     *              WRK(KINDX),WRK(KTOT),KFREE,LFREE)
C        CALL RSPDM(ILSYM,IRSYM,NCLDIM,NCRDIM,CL,CR, RHO1,RHO2,
C    *              ISPIN1,ISPIN2,TDM,NORHO2,XNDXCI,WORK,
C    *              KFREE,LFREE)
C
         IF (TRPLET) THEN
            ISPIN1 = 1
            ISPIN2 = 1
            CALL RSPDM(IREFSY,IREFSY,NCREF,NCREF,WRK(KCREF),WRK(KCREF),
     *                 WRK(KUDV),WRK(KPVX+LPVMAT),
     *                 ISPIN1,ISPIN2,.FALSE.,.FALSE.,
     *                 WRK(KINDX),WRK(KTOT),KFREE,LFREE)
         END IF
C
         IF ( IPRRSP.GT.110 ) THEN
            WRITE(LUPRI,'(/A)') ' ** ONE BODY DENSITY MATRIX FROM RSPDM'
            CALL OUTPUT(WRK(KUDV),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
C
C     For a single open shell, one can use the TDHF code, because for 
C     the open shell the one electron density matrix is 1 and the two
C     electron density matrix is 0.
C
      ELSE IF (TDHF .AND. NASHT .EQ. 1) THEN
         WRK(KUDV) = D1
         WRK(KPVX) = D0
         IF (TRPLET) WRK(KPVX+1) = D0
      END IF
C
      CALL FLSHFO(LUPRI)
C
C Open RSPVEC for response vectors, existing vectors may exist for restart
C
      CALL GPINQ('RSPVEC','EXIST',EX)
      IF (.NOT.EX) THEN
         CALL GPOPEN(LURSP,'RSPVEC','NEW',' ','UNFORMATTED',IDUMMY,
     &        .FALSE.)
         WRITE (LURSP) 'EOFLABEL'
      ELSE
         CALL GPOPEN(LURSP,'RSPVEC','OLD',' ','UNFORMATTED',IDUMMY,
     &        .FALSE.)
      END IF
C
C Open files for trial and sigma vectors
C
      CALL GPOPEN(LURSP3,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LURSP4,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPOPEN(LURSP5,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
C
cs
C      IF (NEXVAL.EQ.0) THEN
Ckr+hjaaj-nov 96: are we sure that we may not have a linear
C                 response where accidentally NEXVAL .eq. 0
C                 e.g. because of symmetry ?
cs
cs       SUBROUTINE QRCALC(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
cs     *                  XINDX,WRK,LWRK)
cs      praticamente chiama la QRVEC di rspvec.F
cs      supponiamo che parta a calcolare 
cs
cs         CALL QRCALC(WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFOCK),
cs     *               WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
cs     *               WRK(KINDX),WRK(KWRK1),LWRK1) 
C      ELSE
        IF (IPRRSP.GE.0) THEN
           WRITE(LUPRI,'(//A,I5)')
     *     ' >>>>>>>>>> Symmetry of excitation/property operator',
     *     KSYMOP
           IF (EXCLC) THEN
              WRITE(LUPRI,'(/A,I5)') ' Number of excitations: ',NEXVAL
           ELSE
              WRITE(LUPRI,'(2(/A,I5))')
     *        ' Number of operators for linear response equations: '
     *        ,NABAOP, ' Number of response frequencies: ',NEXVAL
           ENDIF
        END IF
C
C       DEFINE VARIABLES THAT DEPEND ON SYMMETRY
C
        CALL RSPVAR(WRK(KUDV),WRK(KFOCK),WRK(KFC),WRK(KFV),
     *              WRK(KFCAC),WRK(KH2AC),WRK(KINDX),WRK(KWRK1),LWRK1)
C       CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
        IF ( KZVAR.EQ.0) THEN
           NWARN = NWARN + 1
           WRITE(LUPRI,'(/2A)')' ****ABARSP WARNING******',
     *           ' NUMBER OF VARIABLES IN THIS SYMMETRY IS ZERO'
           GO TO 100
        END IF
C
        IF (EXCLC) THEN
C
C       ... FIND EXCITATION ENERGIES AND TRANSITION MOMENTS
C
           IF (IPRRSP .GT. 0) TIMPP = SECOND()
           CALL RSPLEX(LUSOVE,EXVAL,
     *                 WRK(KCMO),WRK(KUDV),WRK(KPVX),
     *                 WRK(KFC),WRK(KFV),WRK(KFCAC),WRK(KH2AC),
     *                 WRK(KINDX),WRK(KWRK1),LWRK1)
           IF (IPRRSP .GT. 0) THEN
              TIMPP = SECOND() - TIMPP
              WRITE (LUPRI,'(//A,F10.2,A,I2)')
     *           'Time used in polarization propagator calculation is',
     *           TIMPP,' CPU seconds for symmetry',KSYMOP
           END IF
C
        ELSE
           IF (IPRRSP .GT. 0) TIMLR = SECOND()
           CALL RSPLLE(LUGDVE,LUSOVE,LUREVE,NABAOP,NABATY,LABEL,
     *                 WRK(KCMO),WRK(KUDV),WRK(KPVX),WRK(KFC),WRK(KFV),
     *                 WRK(KFCAC),WRK(KH2AC),WRK(KINDX),WRK(KFOCK),
     *                 WRK(KWRK1),LWRK1)
           IF (IPRRSP .GT. 0) THEN
              TIMLR = SECOND() - TIMLR
              WRITE (LUPRI,'(//A,F10.2,A,I2)')
     *           ' Time used in linear response calculation is',
     *           TIMLR,' CPU seconds for symmetry',KSYMOP
           END IF
        END IF
C
C      ENDIF
C      ... end if for QRCALC test above
C
      CALL FLSHFO(LUPRI)
C
  100 CONTINUE
C
      CALL GPCLOSE(LURSP3,'KEEP')
      CALL GPCLOSE(LURSP4,'DELETE')
      CALL GPCLOSE(LURSP5,'KEEP')
      CALL GPCLOSE(LURSP,'KEEP')
      IF (LUSOL.GT.0) CALL GPCLOSE(LUSOL,'KEEP')
C
      TIMRSP = SECOND() - TIMRSP
      IF (IPRRSP .GT. 0) WRITE (LUPRI,'(/A,F13.2,A)')
     *   ' >>> Total time used in response solver is',TIMRSP,
     *   ' CPU seconds.'
C
C     *******************************************
C
      CALL QEXIT('ABARSP')
      RETURN
      END
C  /* Deck rsplex */
      SUBROUTINE RSPLEX(LUSOVE,EXVAL,CMO,UDV,PV,FC,FV,FCAC,
     *                  H2AC,XINDX,WRK,LWRK)
C
C  Purpose:
C     CONTROL CALCULATION OF EXCITATION ENERGIES AND WRITE SOLUTION
C     VECTORS ON DISK (LUSOVE)
C
#include <implicit.h>
#include <dummy.h>
      CHARACTER*8 BLANK
      DIMENSION EXVAL(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(LWRK)
C
      PARAMETER ( D100 = 100.0D0, D0 = 0.0D0, BLANK = '        ' )
#include <codata.h>
#include <priunit.h>
C
C
C Used from common blocks:
C  /INFRSP/ : most items (/INFRSP/ gives control information for
C                         the response calculation(s) )
C  /WRKRSP/ :
C
#include <infrsp.h>
#include <inftap.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infpp.h>
#include <inforb.h>
#include <infpri.h>
C
      CALL QENTER('RSPLEX')
C     space allocation for reduced E(2) and reduced S(2)
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KWRK1  = KEIVEC + MAXRM*MAXRM
      LWRK1  = LWRK + 1 - KWRK1
      IF (IPRPP .GT. 2) THEN
         WRITE(LUPRI,*)' IN RSPLEX: MAXRM      ',MAXRM
         WRITE(LUPRI,*)' IN RSPLEX: LWRK ,LWRK1',LWRK,LWRK1
         WRITE(LUPRI,*)' IN RSPLEX: THCPP      ',THCPP
      END IF
C
C     971201-hjaaj: changed test from 3*KZYVAR to 2*KZYVAR;
C     i.e. we let RSPCTL do the testing but we are sure of
C     enough memory for RSPEVE below.
C
      IF (LWRK1 .LT. 2*KZYVAR) THEN
         WRITE (LUPRI,9000) LWRK1,2*KZYVAR
         CALL QTRACE(LUPRI)
         CALL QUIT('ERROR, RSPLEX: INSUFFICIENT WORK SPACE')
      ENDIF
 9000 FORMAT(/' RSPLEX, work space too small for 2 (z,y)-vectors',
     *       /'         had:',I10,', need more than:',I10)
C
      KZRED  = 0
      KZYRED = 0
      THCRSP = THCPP
      IPRRSP = IPRPP
      MAXIT  = MAXITP
C
C     Call RSPCTL to solve propagator eigen problem
C
      CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *            .FALSE.,BLANK,BLANK,DUMMY,DUMMY,WRK(KREDE),WRK(KREDS),
     *            WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *            XINDX,WRK(KWRK1),LWRK1)
C     CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C    *            LINEQ,GD,REDGD,REDE,REDS,
C    *            IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
C
C CALCULATE EIGENVECTORS AND store on LUSOVE for later
C eigenvalues are stored in EXVAL
C
C  1) MAXIMUM NUMBER OF TRIAL VECTORS
C  2) ALLOCATE WORK SPACE
C
      NSIM   = MIN(KEXCNV, (LWRK1-KZYVAR)/KZYVAR)
      KBVECS = KWRK1
      KWRK2  = KBVECS + NSIM*KZYVAR
      LWRK2  = LWRK   - KWRK2
C
      DO 500 ISIM = 1,KEXCNV,NSIM
         NBX = MIN( NSIM,(KEXCNV+1-ISIM) )
         CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),WRK(KBVECS),
     *               WRK(KWRK2),NBX,(ISIM-1))
C        CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
         DO 550 INUM = 1,NBX
            EXVAL(ISIM-1+INUM) = WRK(KEIVAL-1+ISIM-1+INUM)
            CALL WRITT(LUSOVE,KZYVAR,WRK(KBVECS+(INUM-1)*KZYVAR))
            IF (IPRPP .GT. 0) THEN
               WRITE(LUPRI,'(/A,I5,//A,1P,G16.8,A,3(/20X,G16.8,A),/)')
     *             ' STATE NO:',(ISIM-1+INUM),
     *             ' EXCITATION ENERGY :',WRK(KEIVAL-1+ISIM-1+INUM),
     *             ' au',
     *             WRK(KEIVAL-1+ISIM-1+INUM)*XTEV,  ' eV',
     *             WRK(KEIVAL-1+ISIM-1+INUM)*XTKAYS,' cm-1',
     *             WRK(KEIVAL-1+ISIM-1+INUM)*XKJMOL,' kj / mole'
               IF (SOPPA) THEN
                  TXNRM = DNORM2(KZCONF,WRK(KBVECS+(INUM-1)*KZYVAR),1)
                  TYNRM = DNORM2(KZCONF,WRK(KBVECS+(INUM-1)
     &                           *KZYVAR+KZVAR),1)
                  T2P2HN = D100*(TXNRM*TXNRM + TYNRM*TYNRM)
                  TPHNRM = D100 - T2P2HN
                  WRITE(LUPRI,'(2(/A,F8.2,A))')
     &                ' SOPPA  p-h  weight in excitation operator:',
     &                TPHNRM,' %',
     &                ' SOPPA 2p-2h weight in excitation operator:',
     &                T2P2HN,' %'
               END IF
            END IF
            IF (IPRPP .GT. 2) THEN
               WRITE (LUPRI,'(/A,I3)')
     &             ' EIGENVECTOR FOR STATE NO.',(ISIM-1+INUM)
               CALL RSPPRO(WRK(KBVECS+(INUM-1)*KZYVAR+KZCONF),KZVAR,
     *             LUPRI)
C             CALL RSPPRC(WRK(KBVECS+(INUM-1)*KZYVAR),KZCONF,KZVAR,LUPRI)
               CALL RSPANC(WRK(KBVECS+(INUM-1)*KZYVAR),KZCONF,KZVAR,
     *             MULD2H(KSYMOP,IREFSY),XINDX,MULD2H,LUPRI)
            END IF
            CALL WRTRSP(LURSP,KZYVAR,WRK(KBVECS + (INUM-1)*KZYVAR),
     &                  'EXCITLAB',BLANK,WRK(KEIVAL-1+ISIM-1+INUM),D0,
     &                  KSYMOP,0,WRK(KRESID-1+ISIM-1+INUM),D0)
  550    CONTINUE
  500 CONTINUE
C
C *** END OF RSPLEX
C
      CALL QEXIT('RSPLEX')
      RETURN
      END
C  /* Deck rsplle */
      SUBROUTINE RSPLLE(LUGDVE,LUSOVE,LUREVE,NABAOP,NABATY,LABEL,CMO,
     *                  UDV,PV,FC,FV,FCAC,H2AC,XINDX,FOCK,WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
#include <mxcent.h>
C
      CHARACTER*8 BLANK, LABEL
      PARAMETER (BLANK = '        ')
      LOGICAL RESFLG(8)
      DIMENSION NABATY(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),FOCK(*),WRK(LWRK)
C
#include <priunit.h>
#include <gdvec.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <inflr.h>
#include <inftap.h>
#include <infdim.h>
#include <inforb.h>
#include <abainf.h>
#include <absorp.h>
#include <infvar.h>
#include <dftcom.h>
C
      PARAMETER ( DM1 = -1.0D0 )
C
C     Save flag to control that RESPONANT is only called once
C     for each symmetry
C
      SAVE RESFLG
      DATA RESFLG/.FALSE.,.FALSE.,.FALSE.,.FALSE.,
     &            .FALSE.,.FALSE.,.FALSE.,.FALSE./
C
C DETERMINE SECOND ORDER MOLECULAR PROPERTIES
C
      CALL QENTER('RSPLLE')
C
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KREDE,MAXRM*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDS,MAXRM*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDZ,2*MAXRM*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBTYP,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIVAL,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESID,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIVEC,2*MXFREQ*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDGD,MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDZGD,2*MAXRM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KGD,KZYVAR,WRK,KFREE,LFREE)
C
      IF (ABSORP) THEN
C
C     Nonzero damping parameter corresponds to complex response including 
C     absorption in the system
C
         IF (.NOT.RESFLG(KSYMOP)) THEN
            CALL RESONANT(WRK(KREDE),WRK(KREDS),WRK(KIBTYP),
     &           WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),CMO,UDV,PV,
     &           FOCK,FC,FV,FCAC,H2AC,XINDX,DUMMY,WRK(KFREE),LFREE)
            RESFLG(KSYMOP) = .TRUE.
         END IF
C
         CALL MEMGET('REAL',KRES,2*3*3*NFREQ_INTERVAL,WRK,KFREE,LFREE)
         CALL ABSCTL(KOPER,LABEL,WRK(KREDE),WRK(KREDS),
     &        WRK(KREDZ),WRK(KREDGD),WRK(KREDZGD),
     &        WRK(KEIVEC),WRK(KIBTYP),
     &        CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,
     &        WRK(KRES),WRK(KFREE),LFREE)
C
      ELSE
C
C     Infinite lifetime of the excited states,
C     corresponds to the "regular" real response
C
         IF (LFREE.LT.3*KZYVAR) THEN
            WRITE (LUPRI,9100) LFREE,3*KZYVAR
            CALL QTRACE(LUPRI)
            CALL QUIT('RSPLLE, INSUFFICIENT SPACE TO SOLVE '//
     &           'LINEAR EQUATIONS')
         ENDIF
 9100    FORMAT(/' RSPLLE, work space too small for 3 (z,y)-vectors',
     &          /'         had:',I10,', need more than:',I10)
C     WORK SPACE FOR RSPEVE AND RSPRES
C     WORK SPACE EACH TIME RSPRES IS CALLED
         MWRK   = KZYVAR
         IF (.NOT. RSPCI ) MWRK = MWRK + NCREF +  LACIMX
C     Work space for each residual vector
         MRES   = MAX(NORBT*NORBT,2*N2ASHX) + KZYVAR
C     NUMBER OF SIMULTANEOUS VECTORS
         MXSIM = (LFREE-MWRK)/MRES
         MXSIM = MIN(NFREQ,MXSIM)
         IF (MXSIM .LE. 0) THEN
            WRITE(LUPRI,*) 'RSPLLE: insufficient memory for RSPRES'
            WRITE(LUPRI,*) ' had:',LFREE,' need more than:',MWRK+MRES
            CALL QUIT('RSPLLE: insufficient memory for RSPRES')
         END IF
C     
         KWRKE  = KFREE
         KBVECS = KWRKE + KZYVAR
C     
         KZRED  = 0
         KZYRED = 0
         THCRSP = THCLR
         IPRRSP = IPRLR
         MAXIT  = MAXITL
C     
         DO IOP = 1,NABAOP
            IF (IPRLR .GE. 0)
     &           WRITE (LUPRI,'(//A,I3,/A,I3,A,I3/A,(T25,5F10.6))')
     &           ' RSPLLE -- linear response calculation for symmetry',
     &           KSYMOP,
     &           ' RSPLLE -- operator no.',IOP,' out of',NABAOP,
     &           ' RSPLLE -- frequencies :',(FREQ(I),I=1,NFREQ)
            IF (ABASOP .AND. IPRLR .GE. 0)
     &           WRITE (LUPRI,'(/A/A,I8,/A,I8,/A,I8)')
     &           '  SOPPA calculation:'
     &           ,'  p-h variables.            - KZWOPT :',KZWOPT
     &           ,'  2p-2h variables.          - KZCONF :',KZCONF
     &           ,'  Total number of variables - KZVAR  :',KZVAR
C     
C     READ IN GD VECTOR FROM LUGDVE
C     
            CALL READT(LUGDVE,KZVAR,WRK(KGD))
C     
            IF (NABATY(IOP).EQ.-1) THEN
                DFTIMG = .TRUE.
            ELSE
                DFTIMG = .FALSE.
            END IF
            IF (NABATY(IOP).EQ.1) THEN
               CALL DCOPY(KZVAR,WRK(KGD),1,WRK(KGD+KZVAR),1)
               CALL DSCAL(KZVAR,DM1,WRK(KGD+KZVAR),1)
            ELSE IF (NABATY(IOP).EQ.-1) THEN
               CALL DCOPY(KZVAR,WRK(KGD),1,WRK(KGD+KZVAR),1)
            ELSE
               WRITE(LUPRI,'(/2A,I5,A,I5)')
     *              ' NABATY(IOP) NOT AN ALLOWED VALUE'
     *              ,' IOP=',IOP,'  NABATY(IOP)=',NABATY(IOP)
               CALL QUIT(' RSPLLE: INCORRECT NABATY(IOP)')
            END IF
            DO I = 1,NFREQ
               WRK(KEIVAL-1+I) = FREQ(I)
            END DO
            KEXSIM = NFREQ
            KEXCNV = NFREQ
C     
C     Call RSPCTL to solve linear set of response equations
C     
            CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *           .TRUE.,LABEL,BLANK,WRK(KGD),WRK(KREDGD),
     *           WRK(KREDE),WRK(KREDS),WRK(KIBTYP),WRK(KEIVAL),
     *           WRK(KRESID),WRK(KEIVEC),XINDX,WRK(KFREE),LFREE)
            DFTIMG = .FALSE.
C     
            DO IFREQ = 1,NFREQ,MXSIM
               NBX = MIN(MXSIM,(NFREQ+1-IFREQ))
               IBOFF = IFREQ - 1
               CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *              WRK(KBVECS),WRK(KWRKE),NBX,IBOFF)
               DO IVEC = 1,NBX
                  IBV    = (IVEC-1)*KZYVAR + KBVECS
                  CALL WRITT(LUSOVE,KZYVAR,WRK(IBV))
                  IF (IPRLR .GE. 2) THEN
                     WRITE (LUPRI,'(/A,I5)') ' EIGENVECTOR NUMBER ',
     &                    IFREQ-1+IVEC
                     CALL RSPPRO(WRK(IBV+KZCONF),KZVAR,LUPRI)
                     CALL RSPANC(WRK(IBV),KZCONF,KZVAR,
     *                    MULD2H(KSYMOP,IREFSY),XINDX,MULD2H,LUPRI)
                  END IF
                  NATSYM = 1.0D0
                  CALL WRTRSP(LURSP,KZYVAR,WRK(IBV),LABEL,
     &                 BLANK,WRK(KEIVAL-1+IVEC),D0,KSYMOP,0,
     &                 WRK(KRESID-1+IVEC),ANTSYM)
               END DO
               CALL RSPRES(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *              WRK(KGD),WRK(KFREE),XINDX,LFREE,UDV,NBX,IBOFF)
               DO IVEC = 1,NBX
                  IBV    = (IVEC-1)*KZYVAR + KFREE
                  CALL DSCAL(KZYVAR,DM1,WRK(IBV),1)
                  CALL WRITT(LUREVE,KZYVAR,WRK(IBV))
               END DO
            END DO
C
C     End loop over B operators
C
         END DO
      END IF
C
      CALL QEXIT('RSPLLE')
      RETURN
      END
C  /* Deck prpaba */
      SUBROUTINE PRPABA(WORD,GPLON,LWRK)
C
C GET RIGHT HAND SIDE LONDON VECTOR in GPLON(1)
C for RESPONSE module
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
C
      CHARACTER*(*)WORD
      DIMENSION GPLON(LWRK)
C
C Used from common blocks:
C   INFDIM : NVARMA
C   INFLIN : NVARPT
C
#include <mxcent.h>
#include <maxmom.h>
#include <maxorb.h>
#include <maxaqn.h>
C
#include <inftap.h>
#include <nuclei.h>
#include <symmet.h>
C
#include <infdim.h>
#include <inflin.h>
C
      LOGICAL OLDDX
C
      CALL GPOPEN(LUGDI,ABAGDI,'OLD','DIRECT',' ',IRAT*NVARMA,OLDDX)
C
      IF (WORD(1:7).EQ.'XLONMAG') THEN
         IREC = 3*NUCDEP + IPTAX(1,2)
      ELSE IF (WORD(1:7).EQ.'YLONMAG') THEN
         IREC = 3*NUCDEP + IPTAX(2,2)
      ELSE IF (WORD(1:7).EQ.'ZLONMAG') THEN
         IREC = 3*NUCDEP + IPTAX(3,2)
      ELSE
         WRITE(LUPRI,'(2A)') 'Wrong label in PRPABA : ',WORD
         CALL QUIT('Wrong Label in PRPABA')
      ENDIF
C
      CALL READDX(LUGDI,IREC,IRAT*NVARPT,GPLON)
C     imaginary operator, i.e. Y = Z:
      CALL DCOPY(NVARPT,GPLON(1),1,GPLON(1+NVARPT),1)
C
      CALL GPCLOSE(LUGDI,'KEEP')
C
C END OF PRPABA
C
      RETURN
      END
