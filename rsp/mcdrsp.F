!  /* Deck rspctlmcd */
      SUBROUTINE RSPCTLMCD(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &                  LINEQ,LAB1,LAB2,GD,REDGD,REDE,REDS,
     &                  IBTYP,EIVAL,RESID,EIVEC,XINDX,
     &                  Xf,S2Xf,WRK,LWRK)
!
! Sonia 9/12/2009
! Special version of RSPCTL for projecting redundant states
! out  in MCD calculation
! LINEQ=.false. is disabled!!!!!!!
!
! KEXSTV: NUMBER OF START VECTORS
! KEXSIM: DESIRED NUMBER OF SIMULTANEOUS ROOTS IN A
!         MICROITERATION
! JEXSIM: THE NUMBER OF TRIAL VECTORS IN A MACROITERATION
!  ( TRIAL VECTORS ARE ONLY ADDED FOR THE NON CONVERGED VECTORS AND
!    LINEAR DEPENDENCE IS REMOVED)
!
! MAXRIT: MAXIMUM NUMBER OF MICROITERATIONS
!
#include <implicit.h>
#include <dummy.h>
!
      PARAMETER ( MAXVEC = 12, DM1 = -1.0D0 )
      DIMENSION CMO(*),UDV(NASHDI,*),PV(*),FC(*),FV(*)
      DIMENSION FCAC(*), H2AC(*), GD(*),REDGD(*)
      DIMENSION REDE(*),REDS(*),IBTYP(*),EIVAL(*),RESID(*),EIVEC(*)
      DIMENSION XINDX(*),WRK(*)
      !Sonia: the excitation vector X_f and its  S[2]X_f transformation
      DIMENSION Xf(*),S2Xf(*)
!
#include <ibndxdef.h>
!
! Used from common blocks:
!  /WRKRSP/: ??
!  priunit.h: LUERR
!
#include <priunit.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infdim.h>
#include <infpri.h>
#include <inftap.h>
#include <infsop.h>
!
      LOGICAL BINMEM, LINEQ, OPTSAV, LMAXIT
      CHARACTER*8 LAB123(3),RESTAR, LAB1, LAB2
      DATA LAB123/'********','********','RESPONPP'/
      DATA RESTAR/'RESTART '/
!
      CALL QENTER('RSPCTLMCD')
      IF (IPRRSP .GE. 0) THEN
      IF (LINEQ) THEN
         WRITE (LUPRI,'(///2A//A,I3,A,L3/)')
     *   ' <<<  SOLVING PROJECTED SETS OF LINEAR EQUATIONS ',
     *   'FOR LINEAR RESPONSE PROPERTIES at RESONANT FREQ >>>',
     &   ' Operator symmetry =',KSYMOP,'; triplet = ',TRPLET
      ELSE
         CALL QUIT('RSPCTLMCD NOT TO BE USED FOR EXCI energies')
!         WRITE (LUPRI,'(///2A//A,I3,A,L3/)')
!     *   ' <<< EXCITATION ENERGIES',
!     *   ' AND TRANSITION MOMENT CALCULATION (MCTDHF) >>>',
!     &   ' Operator symmetry =',KSYMOP,'; triplet = ',TRPLET
      END IF
      END IF
!
!  DEFINE MAXIMUM NUMBER OF TRIAL VECTORS IN RSPLIN AND RSPNEX
!
!  WORK SPACE REQUIREMENT IN RSPLIN
!
!  SPECIAL REQUIREMENT FOR AN ORBITAL TRIAL VECTOR
!
!  FROM RSPELI + TRIAL VECTORS
!
      KORB   =  KZYWOP + KZYVAR + 3*NORBT*NORBT + 2*NASHT*NORBT
!
!  FROM RSPOLI
!
      KORB   = KORB + 2*NORBT*NASHT + N2ASHX*NNASHX + 2*NORBT*NORBT
      IF (RSPCI) KORB = 0
!
!  SPECIAL REQUIREMENT FOR A CONFIGURATION TRIAL VECTOR
!
!  FROM RSPELI + TRIAL VECTORS
!
      KCSF   =  KZCONF + KZYVAR + NORBT*NORBT + 2*NORBT*NASHT
!
!  FROM RSPOLI
!
      IF (.NOT.RSPCI) KCSF = KCSF + N2ASHX + N2ASHX*N2ASHX
!
!  MAXIMUM SPACE FOR ONE TRIAL VECTOR
!
      KMAXVE = MAX(KORB,KCSF)
!
!  WORK SPACE REQUIREMENT FOR EACH RSPLIN CALL
!
      KUSE   = 2*NORBT*NORBT + 2*LBINTM
!
!  FROM RSPOLI
!
      KUSE   = KUSE + 2*N2ASHX
!                     2 :TO ASSURE SPACE FOR ALL PVCD MATRICES
!
!  WORK SPACE USED EACH TIME A CSF LINEAR TRANSFORMATION
!  OR A DENSITY MATRIX IS CONSTRUCTED
!
      KLIDEN = N2ASHX  + LACIMX
!  REFERENCE VECTOR ONLY NEEDED FOR MCSCF CALCULATION
      IF (.NOT.RSPCI) KLIDEN = KLIDEN + MAX(KZCONF,NCREF)
!
      KTOUSE = LWRK - KUSE - KLIDEN
!     SOPCLI reads 2p-2h diagonal (D(0) matrix)
      IF (SOPPA) KTOUSE = KTOUSE - KZCONF
      MAXSIM = MIN(MAXVEC,(KTOUSE-100)/KMAXVE)
      IF (IPRRSP.GT.37 .OR. MAXSIM.LE.0) THEN
         WRITE(LUPRI,*)
     *     ' KLIDEN,KUSE,KMAXVE,LBINTM,NORBT,NASHT,KZYVAR,LWRK '
         WRITE(LUPRI,*)
     *       KLIDEN,KUSE,KMAXVE,LBINTM,NORBT,NASHT,KZYVAR,LWRK
         WRITE(LUPRI,*)' LACIMX',LACIMX
      END IF
      IF (IPRRSP.GE.10) WRITE(LUPRI,'(/A,I8)')
     *     ' MAXIMUM NUMBER OF SIMULTANEOUS TRIAL VECTORS:',MAXSIM
!
      IF (MAXSIM.LE.0) THEN
         WRITE (LUPRI ,8000) LWRK,(100+KMAXVE-KTOUSE)
         WRITE (LUERR,8000) LWRK,(100+KMAXVE-KTOUSE)
         CALL QTRACE(LUERR)
         CALL QUIT('RSPCTL: INSUFFICIENT MEMORY ALLOCATION')
      ENDIF
 8000 FORMAT(/' RSPCTLmcd, not enough memory for one trial vector',
     *       /'         LWRK =',I10,
     *       /'         need >',I10,' additional work memory space.')
!
! Space for RSPNEX:
!
      MAXNEX = 30
      MAXORP = 100
!
! REQUIREMENT FOR EACH TRIAL VECTOR
!
      KNEX   = KZYVAR + KZYWOP
!
! REQUIREMENT FOR EACH CALL TO RSPNEX
!
      KUSE   = KZYVAR + NORBT*NORBT + 2*N2ASHX
      IF (.NOT.RSPCI) KUSE = KUSE + MAX(NCREF,KZCONF) + LACIMX
      IF (OPTORB) THEN
         KOPT = 2*MAXORP*(MAXORP+1)/2 + 2*MAXORP + 3*KZYWOP*MAXORP
         KOPT = KOPT + KZYWOP
!
!  MAXIMUM FOR ORPLIN AND ORPSVE
!
         KTOT = KORB - KZYCON  -N2ASHX*NNASHX
         KOPT = KOPT + MAX(KZYVAR,KTOT)
      ELSE
         KOPT = 0
      END IF
! REQUIREMENT FOR EACH CALL TO RSPNEX
      IF ( KOPT .GT. KUSE ) KUSE = KOPT
      KLEFT  = LWRK - KUSE
! RSPNEX reads SOPPA 2p-2h diagonal of length KZCONF
      IF (SOPPA) KLEFT = KLEFT - KZCONF
      NSIM3 = MIN(MAXNEX,(KLEFT-100)/KNEX)
      IF(IPRRSP.GE.10 .OR. NSIM3 .LE. 0)  WRITE(LUPRI,'(/A,I8)')
     *' MAXIMUM NUMBER OF SIMULTANEOUS TRIAL VECTORS IN RSPNEX:',NSIM3
      IF (NSIM3.LE.0) THEN
         WRITE (LUPRI ,8000) LWRK,(100+KLEFT+KNEX)
         WRITE (LUERR,8000) LWRK,(100+KLEFT+KNEX)
         CALL QTRACE(LUERR)
         CALL QUIT('RSPCTLMCD: INSUFFICIENT SPACE')
      ENDIF
      KVECS = 1
!
! CONSTRUCT E2 AND S2 EXPLICITLY USING UNIT TRIAL VECTORS
!
      IF (ABCHK ) THEN
         CALL RSPES2(GD,REDGD,REDE,REDS,EIVAL,EIVEC,
     *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               XINDX,WRK,LWRK)
!        CALL RSPES2(GD,REDGD,REDE,REDS,EIVAL,EIVEC,
!    *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
!    *               XINDX,WRK,LWRK)
         WRITE(LUPRI,'(/A)')' .ABCHK E2 AND S2 CONSTUCTED : STOP'
         CALL QUIT('RSPCTLMCD: ABCHK FINISHED')
      END IF
!
!     Initialize micro iteration counter
!
      ITMIC = 0
!
!     SET UP INITIAL TRIAL VECTORS
!
      IF (IPRRSP.GT.10)
     &   WRITE(LUPRI,*)'RSPCTLMCD: LINEQ,RESTPP,ORBSPC',
     &                  LINEQ,RESTPP,ORBSPC
!
      IF (LINEQ) THEN
         WRITE(LUPRI,*)'Entering LRSTMCD, MAXSIM=',MAXSIM
         !we start with 4 trial vectors: Xf, Xf+, GDp, GDp+
         !WRITE(LUPRI,*)'Entering LRSTMCD, EIVAL=',EIVAL
         CALL LRSTMCD(BINMEM,MAXSIM,IBTYP,EIVAL,EIVEC,GD,
     &             REDGD,REDS,REDE,
     &             CMO,PV,FC,FV,
     &             UDV,FCAC,H2AC,XINDX,
     &             LAB1,LAB2,Xf,S2Xf,WRK,LWRK)
         WRITE(LUPRI,*)'Exiting LRSTMCD, BINMEM=',BINMEM
         IF ((KZRED.GT.0) .AND. (KEXSTV.LE.0)) THEN
            ICTL = 2
            KBCVEC = 1
            KBOVEC = 1
            KEVECS = 1
            KWRK2  = 1
            LWRK2  = LWRK
            NSIM   = 0
            GO TO 1000
         END IF
         IF (ORBSPC) THEN
            DO 313 II = 1,KZWOPT
               GD(II) = GD(KZCONF+II)
 313        CONTINUE
!           CALL DCOPY(KZWOPT,GD(1+KZCONF),1,GD,1)
            DO 314 II = 1,KZWOPT
               GD(KZWOPT+II) = GD(KZVAR+KZCONF+II)
 314        CONTINUE
!           CALL DCOPY(KZWOPT,GD(1+KZVAR+KZCONF),1,
!    *                 GD(1+KZWOPT),1)
            CALL OROPTI(EIVAL,IBTYP,GD,CMO,UDV,PV,FC,FV,FCAC,XINDX,
     &                  WRK,LWRK)
            WRITE(LUPRI,'(/A)')
     & ' STOP: CHECK OF ORBITAL BLOCK OF E(2) AND S(2) FINISHED'
            CALL QUIT(' END OF ORBSPC CHECK')
         END IF
      ELSE
         write(lupri,*)'SONIA, RSPCTLMCD: Lineq=false disabled!!!'
         if (.false.) then
         IF (RESTPP) THEN
!           Restart in PP is only possible in first call
!           (in any later call of RSPCTL the restart vectors will be of
!            another symmetry).
            RESTPP = .FALSE.
            CALL  PPRST(IBTYP,REDS,REDE,UDV,XINDX,WRK,LWRK)
            IF (KZRED .LT. KEXSIM) THEN
               WRITE (LUPRI,'(//A,/A,I5,/A,I5)')
     *         ' RSPCTLMCD: ".RESTPP" impossible, too few trial vects.',
     *         '     number of trial vectors on file :',KZRED,
     *         '     number of roots (".ROOTS")     :',KEXSIM
               CALL QUIT('RESTART PP impossible,too few trial vectors.')
            END IF
            ICTL = 2
            KBCVEC = 1
            KBOVEC = 1
            KEVECS = 1
            KWRK2  = 1
            LWRK2  = LWRK
            NSIM   = 0
            GO TO 1000
         ELSE
            CALL PPST(IBTYP,XINDX,FCAC,H2AC,WRK,MAXSIM,BINMEM,LWRK)
         END IF
         end if !disabling of lineq=false
      ENDIF
!
      IF(IPRRSP.GT.10) WRITE(LUPRI,*)'AFTER START ROUTINE BINMEM',BINMEM

      JEXSIM = KEXSTV 
      write(lupri,*)'SONIA RSPCTLMCD JEXSIM, KEXSTV', JEXSIM, KEXSTV
!     ... JEXSIM is redefined in RSPNEX
!
!     --- START OF MICROITERATION LOOP ---
!
 100  CONTINUE
         CALL FLSHFO(LUPRI)
         ITMIC = ITMIC + 1
         IF (IPRRSP.GE.10)
     *   WRITE(LUPRI,'(/A,I5)') ' ** RSPCTLMCD MICROITER NUMBER',ITMIC
!
!        CHECK IF ALL B VECTORS IN MEMORY
!        IF NOT POSITION LURSP3
!
 200     CONTINUE
         NSIM = MIN(MAXSIM,JEXSIM)
         !IF (IPRRSP.GT.15) WRITE(LUPRI,'(A,3I5)')' NSIM,MAXSIM,JEXSIM',
         IF (IPRRSP.GT.10) WRITE(LUPRI,'(A,3I5)')' NSIM,MAXSIM,JEXSIM',
     *      NSIM,MAXSIM,JEXSIM
         NCSIM = 0
         NOSIM = 0
         DO 105 I=1,NSIM
            IF (IBTYP(KZRED+I+KOFFTY).EQ.JBCNDX) THEN
               NCSIM = NCSIM + 1
            ELSE
               NOSIM = NOSIM +1
            ENDIF
 105     CONTINUE
         IF (IPRRSP.GE.15) WRITE(LUPRI,*)' NCSIM,NOSIM',NCSIM,NOSIM
         KBCVEC= KVECS
         KBOVEC= KBCVEC + KZCONF*NCSIM
         KWRK1 = KBOVEC + KZYWOP*NOSIM
         LWRK1 = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPCTLMCD 1',KWRK1-1,LWRK)
         IF (.NOT.BINMEM) THEN
            CALL REWSPL(LURSP3)
#if !defined (VAR_SPLITFILES)
            IF (KOFFTY.EQ.1) READ (LURSP3)
#else
            JRSP3 = 0
            IF (KOFFTY.EQ.1)
     &           CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#endif
            DO 210 I=1,KZRED
#if defined (VAR_SPLITFILES)
 210           CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
 210           READ (LURSP3)
#endif
            ISTBC = KBCVEC
            ISTBO = KBOVEC
            DO 400 ISIM=1,NSIM
               IF (IBTYP(KZRED+ISIM+KOFFTY).EQ.JBCNDX) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBC),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZCONF) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in RSPCTLMCD')
#else
                  CALL READT(LURSP3,KZCONF,WRK(ISTBC))
#endif
                  ISTBC = ISTBC + KZCONF
               ELSE
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP3,LRSP3,WRK(ISTBO),JRSP3,.TRUE.)
                  IF (LRSP3 .NE. KZYWOP) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP3 in RSPCTLMCD')
#else
                  CALL READT(LURSP3,KZYWOP,WRK(ISTBO))
#endif
                  ISTBO = ISTBO + KZYWOP
               ENDIF
 400        CONTINUE
         END IF
C
         IF (IPRRSP.GT.10) THEN
            WRITE(LUPRI,'(/A)')' BEFORE RSPLIN in RSPCTLMCD'
            IF (NOSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NOSIM,' ORBITAL TRIAL VECTORS'
               CALL OUTPUT(WRK(KBOVEC),1,KZYWOP,1,NOSIM,
     *                     KZYWOP,NOSIM,1,LUPRI)
            ENDIF
            IF (NCSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NCSIM,' CONFIG. TRIAL VECTORS'
               CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,
     *                     KZCONF,NCSIM,1,LUPRI)
            ENDIF
         END IF
         CALL FLSHFO(LUPRI)
         CALL RSPELI(NCSIM,NOSIM,WRK(KBCVEC),WRK(KBOVEC),
     *               CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *               XINDX,WRK(KWRK1),LWRK1)
         CALL FLSHFO(LUPRI)

         KEVECS = KWRK1
         KSOP1  = KEVECS + (NCSIM+NOSIM)*KZYVAR
         IF (SOPPA) THEN
            KWRK2 = KSOP1 + KZYWOP
         ELSE
            KWRK2 = KSOP1
         ENDIF
         LWRK2  = LWRK   - KWRK2
         IF (LWRK2.LT.0) CALL ERRWRK('RSPCTLMCD 2',KWRK2-1,LWRK)
!
!        WRITE E(2)*X ON LURSP5
!
         CALL REWSPL(LURSP5)
         JRSP5 = 0
         IF (KOFFTY.EQ.1) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#else
            READ (LURSP5)
#endif
         ENDIF
         DO 215 I=1,KZRED
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#else
            READ (LURSP5)
#endif
 215     CONTINUE
         ISTEC = KEVECS
         ISTEO = ISTEC + NCSIM*KZYVAR
         DO 216 ISIM=1,NSIM
            IF (IBTYP(KZRED+ISIM+KOFFTY).EQ.JBCNDX) THEN
!
!   If SOPPA only write p-h part of transformed vector to disk
!   and construct 2p-2h part when necessary
!
               IF (SOPPA) THEN
                  CALL DCOPY(KZWOPT,WRK(ISTEC+KZCONF),1,
     *                       WRK(KSOP1),1)
                  CALL DCOPY(KZWOPT,WRK(ISTEC+KZVAR+KZCONF),1,
     *                       WRK(KSOP1+KZWOPT),1)
#if defined (VAR_SPLITFILES)
                  CALL WRITST(LURSP5,KZYWOP,WRK(KSOP1),JRSP5)
               ELSE
                  CALL WRITST(LURSP5,KZYVAR,WRK(ISTEC),JRSP5)
#else
                  CALL WRITT(LURSP5,KZYWOP,WRK(KSOP1))
               ELSE
                  CALL WRITT(LURSP5,KZYVAR,WRK(ISTEC))
#endif
               ENDIF
               ISTEC = ISTEC + KZYVAR
            ELSE
#if defined (VAR_SPLITFILES)
               CALL WRITST(LURSP5,KZYVAR,WRK(ISTEO),JRSP5)
#else
               CALL WRITT(LURSP5,KZYVAR,WRK(ISTEO))
#endif
               ISTEO = ISTEO + KZYVAR
            ENDIF
 216     CONTINUE

         KZRED =KZRED+NSIM
         KZYRED=2*KZRED
         write(lupri,*) 'At << prepare for restart >>, KZRED=', KZRED
!
! PREPARE FOR RESTART
!
         IF (IPRRSP.GT.10) THEN
            WRITE(LUPRI,*)' BEFORE WRITING IBTYP ON LURSP1 LINEQ:',LINEQ
         END IF
         KTOT = KZRED + KOFFTY
         IF (IPRRSP.GT.7)
     *      WRITE(LUPRI,'(/A,I5)')' TRIAL VECTORS FOR RESTART:',KTOT
         LURSP1 = -1
         CALL GPOPEN(LURSP1,'RSPRST.E2C','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND (LURSP1)
         WRITE(LURSP1)LAB123,RESTAR
         WRITE(LURSP1)KSYMOP,KTOT,(IBTYP(I),I=1,KTOT)
!        flush any output buffers
         REWIND (LURSP1)
         CALL GPCLOSE(LURSP1,'KEEP')
!
         IF (IPRRSP.GT.10) THEN
            IF (NOSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NOSIM,' ORBITAL TRIAL VECTORS'
               CALL OUTPUT(WRK(KBOVEC),1,KZYWOP,1,NOSIM,
     *                     KZYWOP,NOSIM,1,LUPRI)
            ENDIF
            IF (NCSIM .GT. 0) THEN
               WRITE(LUPRI,'(/I5,A)')NCSIM,' CONFIG. TRIAL VECTORS'
               CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,
     *                     KZCONF,NCSIM,1,LUPRI)
            ENDIF
            WRITE(LUPRI,'(/I5,A)')NSIM, ' E(2) LINEAR TRANSF. VECTORS'
            CALL OUTPUT(WRK(KEVECS),1,KZYVAR,1,NSIM,
     *                     KZYVAR,NSIM,1,LUPRI)
         ENDIF
         CALL FLSHFO(LUPRI)
         JEXSIM =JEXSIM-NSIM
         IF (JEXSIM.GT.0) THEN
!
!           CALL RSPRED(1,..) INCREASE DIMENSION OF REDUCED RSP EQUATION
!
            !iprsav=iprrsp
            !iprrsp = 25
            CALL RSPRED(1,LINEQ,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *                  EIVAL,EIVEC,WRK(KBCVEC),WRK(KBOVEC),
     *                  UDV,WRK(KEVECS),XINDX,WRK(KWRK2),LWRK2)
            !iprrsp = iprsav
!           CALL RSPRED (ICTL,LINEQ,N,IBTYP,GD,REDGD,REDE,REDS,
!    *                 EIVAL,EIVEC,BCVEC,BOVEC,UDV,EVECS,
!    *                 XINDX,WRK,LWRK)
            GO TO 200
!       ^-----------
         END IF

         !Sonia: we are now increasing reduced space dimension _and_ solving
         ICTL = 3
 1000    CONTINUE
!
!        CALL RSPRED(,3,..) INCREASE DIMENSION OF REDUCED RSP EQUATION
!        AND SOLVE FOR EIGENVALUES AND EIGENVECTORS
!
         write(lupri,*)'RSPCTLMCD: at call REDMCD, KZYRED=', KZYRED
         CALL FLSHFO(LUPRI)
         !iprsav=iprrsp
         !iprrsp = 25
         CALL RSPREDMCD(ICTL,LINEQ,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     *               EIVAL,EIVEC,WRK(KBCVEC),WRK(KBOVEC),
     *               UDV,WRK(KEVECS),XINDX,WRK(KWRK2),LWRK2)
         !iprrsp = iprsav

         !write(lupri,*)'RSPCTLMCD: at exit REDMCD, KZYRED=', KZYRED
         CALL FLSHFO(LUPRI)

C        CALCULATE TRANSITION MOMENT
C        CALL RSPMOM(IBTYP,EIVAL,EIVEC,BVECS,SVECS,PRVEC,NSIM)
C
C SET UP REDUCED MATRIX EXPLICLY AND CHECK BLOCK STRUCTURE IN E2 AND S2
C
         IF (ABSYM) THEN
            CALL E2SYM(IBTYP,CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     &                 WRK,LWRK)
         END IF
C        CREATE KEXSIM NEW LINEAR INDEPENDENT TRIAL VECTORS IN NEXT()
C        FROM THE EIGENVECTORS OF THE REDUCED  RSP EQUATION
C
         IF (KEXSIM.LE.MAXSIM) THEN
            BINMEM=.TRUE.
         ELSE
            BINMEM=.FALSE.
         ENDIF
         OPTSAV = OPTORB
         ITCIN = 0
         DO 510 IIB = 1,(KOFFTY+KZRED)
            IF (IBTYP(IIB).EQ.JBCNDX) ITCIN = ITCIN + 1
 510     CONTINUE
         ITCIN = ITCIN - KOFFTY
         IF (.NOT.LINEQ .AND. ITCIN.LE.2) OPTORB = .FALSE.
C        880226, PJ+HJAAJ: big problems with optorb for
C        RSPPP if no previous orbital trial vectors.
         LMAXIT = (ITMIC .GE. MAXIT)
         !
         !SONIA: modified RSPNEX routine for MCD calculations
         !
         write(lupri,*)'RSPCTLMCD: at enter RSPNEXMCD=', KZYRED
         CALL FLSHFO(LUPRI)
         CALL RSPNEXMCD(LINEQ,LMAXIT,NSIM3,IBTYP,EIVAL,RESID,EIVEC,GD,
     &               CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,
     &               WRK,LWRK,Xf,S2Xf)
         write(lupri,*)'RSPCTLMCD: at exit RSPNEXMCD=', KZYRED
         CALL FLSHFO(LUPRI)
         OPTORB = OPTSAV
         IF (KCONV.EQ.-2) THEN
C           (MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED)
            NWARN = NWARN + 1
            WRITE (LUERR,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED.'
            WRITE (LUPRI,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED.'
         ELSE IF (KCONV.LT.0) THEN
C           (LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTOR )
            NWARN = NWARN + 1
            WRITE (LUERR,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     OF LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTORS'
            WRITE (LUPRI,'(/A/A)')
     *          ' *** RSPCTL WARNING-MICROITERATIONS STOPPED BECAUSE',
     *          '     OF LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTORS'
         ELSE IF(KCONV.GT.0)THEN
C           (CONVERGED)
            IF (IPRRSP .GE. 10)
     &      WRITE(LUPRI,'(/A)')'*** RSPCTLMCD MICROITERATS CONVERGED'
            CALL FLSHFO(LUPRI)
         ELSE
C           (NOT CONVERGED)
            IF (ITMIC.GE.MAXIT) THEN
C              (MAX NO OF MICROITERATIONS REACHED)
               NWARN = NWARN + 1
               WRITE(LUERR,'(/A,I4,A)')
     *         ' *** RSPCTLMCD WARNING-MAXIMUM NUMBER OF MICROITERnS,',
     *         ITMIC,', REACHED'
               WRITE(LUPRI,'(/A,I4,A)')
     *         ' *** RSPCTLMCD WARNING-MAXIMUM NUMBER OF MICROITERnS,',
     *         ITMIC,', REACHED'
            ELSE
               GO TO 100
C     ^-----------------
            END IF
         END IF
!
!     --- END OF MICROITERATION LOOP ---
!
!
!     END OF RSPCTLMCD
!
      CALL QEXIT('RSPCTLMCD')
      RETURN
      END
!---
!  /* Deck lrstmcd */
      SUBROUTINE LRSTmcd(BinMem,MAXSIM,IBTYP,EIVAL,EIVEC,GD,
     *                REDGD,REDS,REDE,
     &                CMO,PV,FC,FV,
     *                UDV,FCAC,H2AC,XINDX,LAB1,LAB2,
     &                Xf,S2Xf,
     &                WRK,LWRK)
!
! PURPOSE: CREATE START VECTOR(S) FOR SOLUTION OF A LINEAR
!          SET OF EQUATIONS in the MCD projected case
!          USE Euclidean normalized Xf vector and 
!          projected GRADIENT VECTOR MULTIPLIED WITH INVERSE
!          DIAGONAL HESSIAN MATRIX ELEMENTS
!
!          IF NTYPE.EQ.1 USE OLD TRIAL VECTORS
!
#include <implicit.h>
#include <dummy.h>
      CHARACTER*8 LAB1, LAB2
      DIMENSION IBTYP(*),EIVAL(*),EIVEC(*),GD(*),REDGD(*)
      DIMENSION REDS(*),REDE(*),UDV(*),FCAC(*),H2AC(*),XINDX(*),WRK(*)
      !Sonia
      DIMENSION Xf(*),S2Xf(*)
      DIMENSION CMO(*),PV(*),FC(*),FV(*)
      PARAMETER ( DTEST = 1.0D-4, D0=0.0D0, D1=1.0D0, D1TEST = 1.0D-6)
      PARAMETER ( DM1 = -1.0D0 )
#include <thrldp.h>
#include <ibndxdef.h>
C
C Used from common blocks:
C   INFRSP : ???,SOPPA
C
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <inftap.h>
#include <infpri.h>
#include <infdim.h>
#include <infopt.h>
#include <phpinf.h>
C
      LOGICAL BinMem   !is output in this context
      LOGICAL FNDLAB, FOUND, CONV
C
      KZRED = 0
      IF ( RESTLR ) THEN
         call quit('LRSTMCD: RESTART NOT YET AN OPTION HERE!!!')
      END IF
      write(lupri,*)'LRSTMCD: KZCONF,KZYWOP,KZYVAR,KZVAR',
     &               KZCONF,KZYWOP,KZYVAR,KZVAR

      JLRSTV = KZRED
      NCSIM  = 0
      NOSIM  = 1
      KBCVEC = 1
      KBOVEC = KBCVEC + NCSIM*KZCONF
      KECVEC = KBOVEC + NOSIM*KZYWOP
      KEOVEC = KECVEC + NCSIM*KZYVAR
      KDIAE  = KEOVEC + NOSIM*KZYVAR
      KWRK   = KDIAE
      LWRKRE = LWRK - KWRK

      ISTBC  = KBCVEC
      ISTBO  = KBOVEC
      ISTEC  = KECVEC
      ISTEO  = KEOVEC
      !call dcopy(KZYVAR,Xf,1,wrk(ISTBO),1)
      call dcopy(KZYWOP,Xf,1,wrk(ISTBO),1)
      if (iprrsp.ge.10) then
         write(lupri,*)'LRSTMCD: Xf in input'
         call rsppro(Xf,KZVAR,-LUPRI)
         write(lupri,*)'LRSTMCD: istbo copied'
         call rsppro(wrk(istbo),KZVAR,-LUPRI)
      end if

      !trial vector is on LURPS3, E[2] transformed on LURSP5
      write(lupri,*)'LRSTMCD: normalize Xf and put on file'
      write(lupri,*)'ISTBO, IBTYP,THRLDP,LURSP3,KOFFTY',
     &     ISTBO, IBTYP(1),THRLDP,LURSP3, KOFFTY
      call flshfo(lupri)
      CALL RSPELI(NCSIM,NOSIM,WRK(ISTBC),WRK(ISTBO),
     &            CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &            XINDX,WRK(KWRK),LWRKRE)
      write(lupri,*)'LRSTMCD: E[2]Xf of unnormalized Eigenvec'
      call rsppro(wrk(KWRK),KZVAR,-LUPRI)
C
C CREATE first NEW ORBITAL TRIAL VECTOR BY normalization of 
C Xf excitation vector, dump on LURSP5 file
C
      nbn = 1  !number of new b's
      nbp = 0  !number of previous b's
      !This step is probably not necessary!!!!
      call rspnor(wrk(ISTBO),nbn,nbp,IBTYP,THRLDP,
     &            WRK(KWRK),LURSP3)
      write(lupri,*)'LRSTMCD: exiting RSPNOR', nbn, nbp
      call rsppro(wrk(istbo),KZVAR,-LUPRI)
      call flshfo(lupri)
      !E[2] transform the normalized vector (again, it is redundant!!)
      CALL RSPELI(NCSIM,NOSIM,WRK(ISTBC),WRK(ISTBO),
     &            CMO,UDV,PV,FC,FV,FCAC,H2AC,
     &            XINDX,WRK(KECVEC),LWRKRE)
      write(lupri,*)'SONIA: LRSTMCD, out of E[2]Xf'
      call rsppro(wrk(kecvec),KZVAR,-LUPRI)
      !dump E[2]Xf on file
      CALL REWSPL(LURSP5)
      JRSP5 = 0
      ISTEC = KECVEC
      ISTEO = ISTEC + NCSIM*KZYVAR
      call rsppro(wrk(ISTEO),KZVAR,-LUPRI)
      !DO 216 ISIM=1,NSIM
      ISIM = 1
      IF (IBTYP(KZRED+ISIM+KOFFTY).EQ.JBONDX) THEN
#if defined (VAR_SPLITFILES)
         CALL WRITST(LURSP5,KZYVAR,WRK(ISTEO),JRSP5)
#else
         CALL WRITT(LURSP5,KZYVAR,WRK(ISTEO))
#endif
         !ISTEO = ISTEO + KZYVAR
      ENDIF
      write(lupri,*)'SONIA : what is NSIM?', NSIM
      NSIM = NOSIM + NCSIM
      write(lupri,*)'SONIA : and what is NSIM now?', NSIM
      KZRED=KZRED+NSIM
      KZYRED=2*KZRED
      JLRSTV=JLRSTV-NSIM
!
!  CALL RSPRED(1,..) TO INCREASE DIMENSION OF REDUCED RSP EQUATION
!
      CALL RSPRED(1,.TRUE.,NSIM,IBTYP,GD,REDGD,REDE,REDS,
     &              EIVAL,EIVEC,WRK(KBCVEC),WRK(KBOVEC),
     &              UDV,WRK(KECVEC),XINDX,WRK(KWRK),LWRKRE)
      IF (JLRSTV.LE.0) THEN
         IF (IPRRSP.GE.0) WRITE (LUPRI,'(//A,I5,A/)')
     &             ' ***',KZRED,' FIRST TRIAL VECTOR is X_f'
      ELSE
               write(lupri,*)'SONIA: I HONESTLY DO NOT KNOW'
      END IF

      IF ((.NOT. RESTLR). OR. (KZRED .EQ. 0)) THEN
         IF (KOFFTY.EQ.1) THEN
            write(lupri,*)'Sonia: do we enter here?'
            call flshfo(lupri)
            CALL REWSPL(LURSP3)
            CALL REWSPL(LURSP5)
#if defined (VAR_SPLITFILES)
            JRSP3 = 0
            JRSP5 = 0
#endif
            CALL GETREF(WRK,KZCONF)
#if defined (VAR_SPLITFILES)
            CALL WRITST(LURSP3,KZCONF,WRK,JRSP3)
            CALL WRITST(LURSP5,1,DUMMY,JRSP5)
#else
            CALL WRITT(LURSP3,KZCONF,WRK)
            WRITE (LURSP5) D0, D0, D0, D0
#endif
Chj         ... rec. 1 on LURSP5 should never be read when KOFFTY .eq. 1
            IBTYP(1) = JBCNDX
            IF (IPRRSP.GT.75) THEN
               WRITE(LUPRI,'(/A)')
     *         ' ** LRSTMCD ** REF. VECTOR WRITTEN ON LURSP3 rec. # 1'
            END IF
         ENDIF
      ENDIF

C CREATE second NEW ORBITAL (AND CSF) TRIAL VECTORS BY DIVIDING projected 
C GRADIENTS WITH DIAGONAL HESSIAN ELEMENTS

      KDIAE  = 1
      KDIASO = KDIAE  + LPHPMX
      KTOT   = KDIASO + KZWOPT
      LTOT   = LWRK   - KTOT
      IF (LTOT.LT.0) CALL ERRWRK('LRSTMCD 2',KTOT-1,LWRK)
      IF (KZCONF.GT.0) THEN
         CALL PHPDSK(WRK(KDIAE),LPHPMX,.FALSE.)
      ELSE
         CALL RSPEDG(WRK(KDIAE))
      END IF
      IF (KZWOPT.GT.0) CALL RSPSOD(WRK(KDIASO))

      NLOAD = 0
      NTSIM = 1

      IF (IPRRSP.GT.10) THEN
         WRITE(LUPRI,*)'  GRADIENT VECTOR (Z-part and Y-part)'
         call rsppro(GD,KZVAR,-LUPRI)
      ENDIF
      IF (IPRRSP.GT.10) THEN
         WRITE(LUPRI,*)' DIAGONAL HESSIAN MATRIX WVAL=0.0'
         CALL OUTPUT(WRK(KDIAE),1,KZVAR,1,1,KZVAR,1,1,LUPRI)
      ENDIF

      write(lupri,*)'LRSTMCD at start KEXSIM=', KEXSIM

      DO 150 ISIM = 1,KEXSIM,MAXSIM
         NBX = MIN(MAXSIM,(KEXSIM+1-ISIM))
         write(lupri,*)'LRSTMCD at start nbx=', nbx
         KBCVEC = KTOT
         KBOVEC = KBCVEC + 2*NBX*KZCONF
         KWRK1  = KBOVEC + NBX*KZYWOP
         LWRK1  = LWRK   - KWRK1 - KZYVAR
C
C     We might read in a response vector from file, and do this in
C     WRK(KWRK1).
C
         IF (LWRK1.LT.0) CALL ERRWRK('LRSTMCD 3',KWRK1-1,LWRK)

         ISTBC  = KBCVEC
         ISTBO  = KBOVEC 

         NLOAD = NLOAD + 1
         NCSIM = 0
         NOSIM = 0

         write(lupri,*)'KBOVEC', KBOVEC
         DO 160 IR = 1,NBX
            IROFF = ISIM-1+IR
            WVAL  = EIVAL(IROFF)
!            CALL REARSP(LURSP,KLEN,WRK(KWRK1),LAB1,LAB2,WVAL,D0,
!     &                  KSYMOP,0,THCRSP,FOUND,CONV,ANTSYM)
!            IF (FOUND .AND. .NOT. CONV) THEN
!               write(lupri,*)'SONIA LRSTMCD: I found smthg!'
!               IF (KZCONF .GT. 0) THEN
!                  CALL DCOPY(KZCONF,WRK(KWRK1),1,WRK(ISTBC),1)
!                  NCSIM = NCSIM + 1
!                  ISTBC = ISTBC + KZCONF
!                  IF (WVAL .NE. D0) THEN
!                     CALL DCOPY(KZCONF,WRK(KWRK1+KZVAR),1,WRK(ISTBC),1)
!                     NCSIM = NCSIM + 1
!                     ISTBC = ISTBC + KZCONF
!                  END IF
!               END IF
!               IF (KZWOPT .GT. 0) THEN
!                  CALL DCOPY(KZWOPT,WRK(KWRK1+KZCONF),1,WRK(ISTBO),1)
!                  CALL DCOPY(KZWOPT,WRK(KWRK1+KZVAR+KZCONF),1,
!     &                 WRK(ISTBO+KZWOPT),1)
!                  NOSIM = NOSIM + 1
!                  ISTBO = ISTBO + KZYWOP
!               END IF
!               IF (WVAL .LT. 0) CALL DSWAP(KZWOPT,WRK(KWRK1),1,
!     &                                     WRK(KWRK1 + KZWOPT),1)
!               GO TO 160
!            END IF
!removed config part
            IF (KZWOPT.GE.1) THEN
               write(lupri,*)'SONIA LRSTMCD: prepare preconditioner'
               write(lupri,*)'KZWOPT;', KZWOPT
               DO 105 I=1,KZWOPT
                  IZ = I + KZCONF
                  ZDIA = WRK(KDIAE-1+IZ) - WVAL * WRK(KDIASO-1+I)
                  YDIA = WRK(KDIAE-1+IZ) + WVAL * WRK(KDIASO-1+I)
!                  write(lupri,*)'LRSTMCD: precond', ZDIA,YDIA, ISTBO
                  IF (ABS(ZDIA).LT.DTEST ) THEN
                     WRK(ISTBO-1+I) = GD(KZCONF+I)
     *                           /SIGN(DTEST,ZDIA)
                  ELSE
                     WRK(ISTBO-1+I) = GD(KZCONF+I) / ZDIA
                  ENDIF
!                  write(lupri,*)'WRK(ISTBO-1+I)', WRK(ISTBO-1+I)
                  IF (ABS(YDIA).LT.DTEST ) THEN
                     WRK(ISTBO-1+KZWOPT+I) = GD(KZVAR+KZCONF+I)
     *                          /SIGN(DTEST,YDIA)
                  ELSE
                     WRK(ISTBO-1+KZWOPT+I) = GD(I+KZVAR+KZCONF) / YDIA
                  ENDIF
!                  write(lupri,*)'WRK(ISTBO-1+KZWOPT+I)', 
!     &                           WRK(ISTBO-1+KZWOPT+I)
 105           CONTINUE
               write(lupri,*)'LRSTMCD: NOSIM,ISTBO =',NOSIM,ISTBO
               NOSIM  = NOSIM  + 1
               ISTBO = ISTBO + KZYWOP
               write(lupri,*)'LRSTMCD: NOSIM,ISTBO =',NOSIM,ISTBO
            ENDIF
 160     CONTINUE
         NLSIM = NOSIM + NCSIM
         DO 170 IR = 1 , NLSIM
            IF (IR.LE.NCSIM) IBTYP(KOFFTY+KZRED+NTSIM+IR)=JBCNDX
            IF (IR.GT.NCSIM) IBTYP(KOFFTY+KZRED+NTSIM+IR)=JBONDX
 170     CONTINUE
!         IF ((IPRRSP.GT.100).AND.(NCSIM.GT.0)) THEN
         IF ((IPRRSP.GT.10).AND.(NCSIM.GT.0)) THEN
            WRITE(LUPRI,'(/A,I5)')' NCSIM CSF TRIAL VECTORS',NCSIM
!            CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,
!     &                  KZCONF,NCSIM,1,LUPRI)
         ENDIF
!         IF ((IPRRSP.GT.100).AND.(NOSIM.GT.0)) THEN
         IF ((IPRRSP.GT.10).AND.(NOSIM.GT.0)) THEN
            WRITE(LUPRI,'(/A,I5)')' NOSIM (Z,Y) ORBITAL TRIAL VECTORS',
     &           NOSIM
            !CALL OUTPUT(WRK(KBOVEC),1,KZWOPT,1,2*NOSIM,
            CALL OUTPUT(WRK(KBOVEC),1,KZWOPT,1,2*NOSIM,
     &                  KZWOPT,2*NOSIM,1,LUPRI)
         ENDIF
         IF (IPRRSP.GT.80) THEN
            WRITE(LUPRI,'(/A,I5)')
     *         '    I       IBTYP(I)          I=1,NLSIM NLSIM',NLSIM
            DO 3007 I=1+KOFFTY,NLSIM+KOFFTY
               WRITE(LUPRI,'(I10,I10)')I,IBTYP(I)
 3007       CONTINUE
         END IF
         THRLDV = KZYVAR*THRLDP
!
!        hj-aug2000: we must make sure that configuration and orbital
!        trial vectors are consecutive in memory before CALL RSPORT
!
         ISTBO = KBCVEC+NCSIM*KZCONF 
         IF (NOSIM .GT. 0 .AND. KBOVEC .NE. ISTBO) THEN
            DO I = 1, NOSIM*KZYWOP
               WRK(ISTBO-1+I) = WRK(KBOVEC-1+I)
            END DO
         END IF
         !Orthonormalize 
         nold=1
         CALL RSPORTmcd(WRK(KBCVEC),NLSIM,NOLD,IBTYP,
     *               THRLDV,WRK(KWRK1),Xf,S2Xf,LURSP3)
         IF (IPRRSP.GT.80) THEN
            WRITE(LUPRI,'(/A)')' AFTER RSPORTMCD IN LRSTMCD'
            WRITE(LUPRI,'(/A,I5)')
     *         '    I       IBTYP(I)          I=1,NLSIM NLSIM=',NLSIM
            DO 3008 I=1+KOFFTY,NLSIM+KOFFTY
               WRITE(LUPRI,'(I5,I10)')I,IBTYP(I)
 3008       CONTINUE
         END IF
         NCSIM = 0
         NOSIM = 0
         DO 3010 IX = 1,NLSIM
            IF (IBTYP(KZRED+KOFFTY+NTSIM+IX).EQ.JBCNDX) THEN
               NCSIM = NCSIM + 1
            ELSE
               NOSIM = NOSIM + 1
            ENDIF
 3010    CONTINUE
         IF ((IPRRSP.GT.110).AND.(NCSIM.GT.0)) THEN
            WRITE(LUPRI,*)' NCSIM CSF TRIAL VECTORS',NCSIM
            CALL OUTPUT(WRK(KBCVEC),1,KZCONF,1,NCSIM,KZCONF,NCSIM,1,
     &                  LUPRI)
         ENDIF
         IF ((IPRRSP.GT.110).AND.(NOSIM.GT.0)) THEN
            WRITE(LUPRI,*)' NOSIM ORBITAL TRIAL VECTORS',NOSIM
            CALL OUTPUT(WRK(KBCVEC+NCSIM*KZCONF),1,KZYWOP,
     *                1,NOSIM,KZYWOP,NOSIM,1,LUPRI)
         ENDIF
         NTSIM = NTSIM + NLSIM
 150  CONTINUE
      !SONIA
      !here we decide how many vectors we start with!
      !Info is passed back to rspctl
      !KEXSTV = NTSIM
      KEXSTV = NTSIM - 1
      write(lupri,*)'LRSTMCD: we start with KEXSTV=', KEXSTV
      IF ((NLOAD.EQ.1).AND.(KEXSTV.LE.MAXSIM)) THEN
         BINMEM = .TRUE.
         NTOT = NCSIM*KZCONF + NOSIM*KZYWOP
         IF (KBCVEC .NE. 1) THEN
            DO 156 II = 1,NTOT
               WRK(II) = WRK(KBCVEC-1+II)
 156        CONTINUE
C           CALL DCOPY(NTOT,WRK(KBCVEC),1,WRK(1),1)
         END IF
         write(lupri,*)'LRSTMCD: does it apply for us?'
         write(lupri,*)'B in Mem, Ntot, ', BINMEM, NTOT
         write(lupri,*)'KEXSTV', KEXSTV
      ELSE
         BINMEM = .FALSE.
      ENDIF
C
      IF (KEXSTV .LE. 0) THEN
         IF (IPRRSP.GE.3) WRITE (LUPRI,'(//A,I5)')
     *     ' LRSTMCD, START VECTOR IS NOT LINEAR INDEPEND ',KEXSTV
         IF (KZRED.LE.0) THEN
       CALL QUIT(' LRSTMCD, START VECTORS NOT LINEAR INDEPENDENT.')
         END IF
      END IF
C
C     END OF LRSTMCD.
C
      RETURN
      END
!-------------------------------------------
!  /* Deck rspnor */
      SUBROUTINE RSPNOR(BVECS,NBNEW,NBPREV,IBTYP,THRLDP,OLDVEC,LU3)
!
! Purpose:
! Normalize the given Xf vector. It's based on RSPORT and still contains
! a lot of RSPORT crap (including comments below).
! Sonia, Jan 2010

!  Orthogonalize new b-vectors against all previous b-vectors
!  and among themselves, and renormalize.
!  Each b-vector is in (Z, Y) form, the (Y, Z) vector to be obtained
!  by switching first and second half. Each of the new b-vectors are
!  orthogonalized both against (Z, Y) and (Y, Z) for each previous
!  b-vector.
!  (Orthogonalization is performed twice if round-off is large,
!   if larger than THRRND).
!
! Input:
!  BVECS,  non-orthogonal b-vectors
!  NBNEW,  number of b-vectors in BVECS
!  NBPREV, number of previous b-vectors on LU3
!  IBTYP,  b-vector type (CI or obital)
!  THRLDP, threshold for linear dependence
!
! Output:
!  BVECS,  orthogonal b-vectors (also written to LU3)
!  NBNEW,  number of b-vectors written to LU3
!
! Scratch:
!  OLDVEC(KZYVAR), scratch array for old b-vectors on LU3
!
#include <implicit.h>
#include <dummy.h>
      DIMENSION BVECS(*), IBTYP(*), OLDVEC(*)
C
#include <ibndxdef.h>
C
C    Change value of ZEQLY for problem with linear equation
C    95-02-08
C     PARAMETER (ZEQLY = 1.D-6, T1MIN = 1.D-8, THRRND = 1.D-4)
C
      PARAMETER (OVLMIN = 1.D-4, T1MIN = 1.D-8, THRRND = 1.D-4)
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C  INFRSP: RSPSUP
C  WRKRSP: KZCONF,KZWOPT,KSYMOP
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C     Local arrays:
C
      PARAMETER (MAXVEC = 120)
      DIMENSION JLDP(MAXVEC)
C
      WRITE (LUPRI,*) 'I AM INSIDE RSPNOR: NBNEW ', NBNEW
      IF (NBNEW .GT. MAXVEC) THEN
         WRITE (LUPRI,*) 'ERROR in RSPNOR: NBNEW .gt. MAXVEC'
         WRITE (LUPRI,*) 'NBNEW  =',NBNEW
         WRITE (LUPRI,*) 'MAXVEC =',MAXVEC
         WRITE (LUPRI,*) 'Reduce problem size or recompile RSPORT '//
     *                  'with larger MAXVEC parameter'
         CALL QUIT('RSPNOR error: NBNEW .gt. MAXVEC parameter')
      END IF
!
! *** Average new trial vectors for supersymmetry
!
      IF ( RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         IBVEC = 1
         DO 2600 IRX = 1, NBNEW
            IF (IBTYP (NBPREV+IRX) .EQ. JBCNDX) THEN
!           -- Configuration trial vector; shift the starting point
               IBVEC = IBVEC + KZCONF
            ELSE IF (IBTYP(NBPREV+IRX) .EQ. JBONDX) THEN
               WRITE (LUPRI,*) 'RSPNOR: Supersymmetry\?'
!           -- Orbital trial vector; do the averaging
               CALL RSPAVE(BVECS(IBVEC),KZWOPT,2)
               IBVEC = IBVEC + KZYWOP
            END IF
 2600    CONTINUE
      END IF
!
!     For zero frequency linear response, the result will
!     be Z = Y for real perturbations and Z = -Y for imaginary
!     perturbations; trial vectors will have the same structure.
!
!     However, if Z = Y or Z = -Y then (Z Y) and (Y Z) are linear
!     dependent and we want instead (Z 0) and (0 Z) as trial vectors.
!
      IBVEC = 1
      NCNEW = 0
      NONEW = 0
      DO 1400 IRX = 1,NBNEW
         JLDP(IRX) = 1
         IF (IBTYP(NBPREV+IRX) .EQ. JBCNDX) THEN
            NCNEW = NCNEW + 1
            IBVEC = IBVEC + KZCONF
         ELSE IF (IBTYP(NBPREV+IRX) .EQ. JBONDX) THEN
            NONEW = NONEW + 1
            T1    = DDOT(KZYWOP,BVECS(IBVEC),1,BVECS(IBVEC),1)
            JBVEC = IBVEC + KZWOPT
            OVLPI = D2*DDOT(KZWOPT,BVECS(IBVEC),1,BVECS(JBVEC),1)
!           test if Z + Y is zero
            TZPY  = T1 + OVLPI
            IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A/,(A,1P,D15.8))')
     *         ' ** TEST OF Z = Y OR Z = -Y:',
     *         ' (Z Y) norm squared =',T1,
     *         ' Z + Y norm squared =',TZPY
            IF (TZPY .LE. OVLMIN*T1) THEN
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,I4/A)')
     *            ' Z = -Y in trial vector no.',IRX,
     *            ' Y component removed.'
               CALL DZERO(BVECS(JBVEC),KZWOPT)
            ELSE IF (TZPY .GT. T1) THEN
!              test if Z - Y is zero;
!              only if TZPY = 2 * T1 can Z - Y be zero,
!              so TZPY .gt. T1 is a safe test.
               TZMY = T1 - OVLPI
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(A,1P,D15.8)')
     *            ' Z - Y norm squared =',TZMY
               IF (TZMY .LE. OVLMIN*T1) THEN
                  IF (IPRRSP .GE.10) WRITE (LUPRI,'(/A,I4/A)')
     *            ' Z = Y in trial vector no.',IRX,
     *            ' Y component removed.'
                  CALL DZERO(BVECS(JBVEC),KZWOPT)
               END IF
            END IF
            IBVEC = IBVEC + KZYWOP
         ELSE
            GO TO 9100
         END IF
 1400 CONTINUE
C
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
!
!        Orthogonalize new b-vectors agains previous b-vectors
!        (both (Z, Y) and (Y, Z))
!
         CALL REWSPL(LU3)
#if defined (VAR_SPLITFILES)
         JRSP3 = 0
#endif
         DO 2000 K = 1,NBPREV
            write(lupri,*)'RSPNOR: I should not enter here!!!'
            IBTYPK = IBTYP(K)
            IF (IBTYPK .EQ. JBCNDX) THEN
               IF (NCNEW .EQ. 0) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LU3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
                  READ (LU3)
#endif
                  GO TO 2000
!        v------------------
               END IF
               LBZ = KZCONF
               LBZY= KZCONF
               LBSKIP = 2*KZWOPT
            ELSE IF (IBTYPK .EQ. JBONDX) THEN
               IF (NONEW .EQ. 0) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LU3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
                  READ (LU3)
#endif
                  GO TO 2000
!        v------------------
               END IF
               LBZ = KZWOPT
               LBZY= 2*KZWOPT
               LBSKIP = KZCONF
            ELSE
               GO TO 9100
            END IF
#if defined (VAR_SPLITFILES)
            CALL READST(LU3,LRSP3,OLDVEC,JRSP3,.TRUE.)
            IF (LRSP3 .NE. LBZY) CALL QUIT('Inconsistent record '//
     &           'length on unit LU(RSP)3 in RSPORT')
#else
            CALL READT(LU3,LBZY,OLDVEC)
#endif
            IBVEC = 1
            DO 1900 IRX = 1,NBNEW
               IF (IBTYP(NBPREV+IRX) .NE. IBTYPK) THEN
                  IBVEC = IBVEC + LBSKIP
               ELSE IF (JLDP(IRX).EQ.0) THEN
                  IBVEC = IBVEC + LBZY
               ELSE
!                 Orthogonalize to (Z Y)
                  TT = -DDOT(LBZY,   OLDVEC,1,BVECS(IBVEC),1)
                  CALL DAXPY(LBZY,TT,OLDVEC,1,BVECS(IBVEC),1)
!                 Orthogonalize to (Y Z)
                  IF ( IBTYPK .EQ. JBONDX ) THEN
                     JBVEC = IBVEC + LBZ
                     TT = -DDOT(LBZ,   OLDVEC,       1,BVECS(JBVEC),1)
     *                    -DDOT(LBZ,   OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
                     CALL DAXPY(LBZ,TT,OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
                     CALL DAXPY(LBZ,TT,OLDVEC,       1,BVECS(JBVEC),1)
                  END IF
                  IBVEC = IBVEC + LBZY
               END IF
 1900       CONTINUE
 2000    CONTINUE
!
!        Orthogonalize new vectors against each other
!
         IBVEC  = 1
         DO 2400 IRX = 1,NBNEW
            IBTYPI = IBTYP(NBPREV+IRX)
            IF (IBTYPI .EQ. JBCNDX) THEN
               LBZ = KZCONF
               LBZY= KZCONF
               LBSKIP = 2*KZWOPT
            ELSE
               LBZ = KZWOPT
               LBZY= 2*KZWOPT
               LBSKIP = KZCONF
            END IF
         IF (JLDP(IRX) .NE. 0) THEN
            JBVEC = 1
            DO 2300 JRX = 1,(IRX-1)
            IF (IBTYP(NBPREV+JRX) .NE. IBTYPI) THEN
               JBVEC = JBVEC + LBSKIP
            ELSE IF (JLDP(JRX).EQ.0) THEN
               JBVEC = JBVEC + LBZY
            ELSE
               T1 =  DDOT(LBZY,BVECS(JBVEC),1,BVECS(JBVEC),1)
               TT = -DDOT(LBZY,BVECS(JBVEC),1,BVECS(IBVEC),1)/T1
               write(lupri,*)'RSPNOR: should I enter here\?'
               CALL DAXPY(LBZY,TT,BVECS(JBVEC),1,BVECS(IBVEC),1)
               IF (IBTYPI.EQ.JBONDX) THEN
               write(lupri,*)'RSPNOR: and should I enter here\?'
                  TT=(-DDOT(LBZ,BVECS(JBVEC+LBZ),1,BVECS(IBVEC),1)
     *                -DDOT(LBZ,BVECS(JBVEC),1,BVECS(IBVEC+LBZ),1))/T1
                  CALL DAXPY(LBZ,TT,BVECS(JBVEC+LBZ),1,BVECS(IBVEC),1)
                  CALL DAXPY(LBZ,TT,BVECS(JBVEC),1,BVECS(IBVEC+LBZ),1)
               END IF
               JBVEC = JBVEC + LBZY
            END IF
 2300       CONTINUE
         END IF
!
! NORMALIZE VECTOR NUMBER IRX
!
         IF (JLDP(IRX) .NE. 0) THEN
            TT = DDOT(LBZY,BVECS(IBVEC),1,BVECS(IBVEC),1)
            IF (TT .LE. THRLDP) THEN
               IF (IPRRSP.GT.10) WRITE (LUPRI,3250) IRX,JLDP(IRX),TT
               JLDP(IRX) = 0
            ELSE IF (TT .LT. THRRND) THEN
               IF (ITURN .EQ. 1) THEN
                  IROUND = IROUND + 1
               ELSE
                  IF (IPRRSP.GT.10) WRITE (LUPRI,3255) IRX,JLDP(IRX),TT
                  JLDP(IRX) = 0
               END IF
            END IF
            IF (JLDP(IRX) .NE. 0) THEN
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,2I6,1P,D12.5)')
     &          ' RSPNOR: b-vector, root no., norm**2:',IRX,JLDP(IRX),TT
               IF (TT .LT. T1MIN) THEN
                  TT = D1 / SQRT(TT)
                  CALL DSCAL(LBZY,TT,BVECS(IBVEC),1)
                  TT = DDOT(LBZY,BVECS(IBVEC),1,BVECS(IBVEC),1)
               END IF
               TT = D1 / SQRT(TT)
               CALL DSCAL(LBZY,TT,BVECS(IBVEC),1)
               TT = DDOT(LBZY,BVECS(IBVEC),1,BVECS(IBVEC),1)
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,2I6,1P,D12.5)')
     &          ' RSPNOR 2: b-vec, root no., norm**2:',IRX,JLDP(IRX),TT
               call flshfo(lupri)
            END IF
         END IF
!
! Perform symmetric orthonormalization of (Z Y) and (Y Z) pair
! for vector number IRX
!
!   -1/2       ( C1   C2 )               (  1     OVLPI )
!  S      =   (           ) where S  =  (                )
!              ( C2   C1 )               ( OVLPI     1  )
!
         IF ((JLDP(IRX) .NE. 0).AND.(IBTYPI.EQ.JBONDX)) THEN
            JBVEC = IBVEC + LBZ
            JTURN = 0
 2350       JTURN = JTURN + 1
            OVLPI = D2*DDOT(LBZ,BVECS(IBVEC),1,BVECS(JBVEC),1)
            IF (IPRRSP .GE. 0) THEN
               X1 = DNORM2(LBZY,BVECS(IBVEC),1)
               WRITE (LUPRI,'(A,1P,D14.6,D19.10)')
     *         ' RSPNOR, norm and <ZY/YZ> overlap before S-1/2',X1,OVLPI
            call flshfo(lupri)
            END IF
            X2 = ABS(OVLPI)
            IF (ABS(D1-X2) .LT. OVLMIN) THEN
               IF (IPRRSP.GT.0) WRITE (LUPRI,3260) IRX,JLDP(IRX),OVLPI
               call flshfo(lupri)
               JLDP(IRX) = 0
               GO TO 2390
            END IF
            X1 = D1+OVLPI
            X2 = D1-OVLPI
            IF (ABS(X1 - X2) .GT. THRLDP) THEN
               X1 = DP5 / SQRT(X1)
               X2 = DP5 / SQRT(X2)
               C1 = X1 + X2
               C2 = X1 - X2
               !CALL DCOPY(LBZY,BVECS(IBVEC),1,OLDVEC,1)
               !CALL DSCAL(LBZY,C1,BVECS(IBVEC),1)
               !CALL DAXPY(LBZ,C2,OLDVEC,1,BVECS(JBVEC),1)
               !CALL DAXPY(LBZ,C2,OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
            IF (JTURN .EQ. 1) GO TO 2350
!           ...930712-hjaaj: an output with LBZY only = 480 showed
!              that OVLPI ca. 1.0D-14 after first S-1/2;
!              but ca. 1.0D-16 after second S-1/2
               IF (IPRRSP .GE. 0) THEN
                  X1 = DNORM2(LBZY,BVECS(IBVEC),1)
                  OVLPI= D2*DDOT(LBZ,BVECS(IBVEC),1,BVECS(JBVEC),1)
                  WRITE (LUPRI,'(A,1P,D14.6,D19.10)')
     *            ' RSPNOR, norm and <ZY/YZ> overlap after  S-1/2',
     *            X1,OVLPI
               END IF
            END IF
         END IF
C
C
 2390    IBVEC = IBVEC + LBZY
 2400    CONTINUE
         IF (IROUND.GT.0) THEN
            IROUND = 0
            IF (ITURN.EQ.1) GO TO 1500
C     ^-------------------------------
            CALL QUIT('RSPNOR sw-error; IROUND .gt. 0 and ITURN .gt. 1')
         END IF
 3250 FORMAT(/' (RSPNOR) b-vector no.',I3,' corresponding to ',
     *        'root no.',I3/,14X,'is removed because of linear ',
     *        'dependence; ',1P,D12.5)
 3255 FORMAT(/' (RSPNOR) b-vector no.',I3,' corresponding to root ',
     *        'no.',I3,' is removed',/,14X,'because of small',
     *        'norm**2 in 2. Gram-Schmidt orthonormalization:',1P,D12.5)
 3260 FORMAT(/' (RSPNOR) b-vector no.',I3,' corresponding to ',
     *        'root no.',I3,' is removed',/,14X,'because of lin.dep.,',
     *        '<Z Y|Y Z> overlap =',1P,D14.6)
C
C Save new vector together with old ones on LU3
C
      IBVECX = 1
      IBVEC  = 1
      IBPREV = NBPREV
      DO 3300 IRX = 1,NBNEW
         IBTYPI = IBTYP(NBPREV+IRX)
         IBTYP(NBPREV+IRX) = 0
         IF (IBTYPI .EQ. JBCNDX) THEN
            LBZY= KZCONF
         ELSE
            LBZY= 2*KZWOPT
         END IF
         IF (JLDP(IRX) .NE. 0) THEN
            IBPREV = IBPREV + 1
            IBTYP(IBPREV) = IBTYPI
C
#if defined (VAR_SPLITFILES)
            CALL WRITST(LU3,LBZY,BVECS(IBVEC),JRSP3)
#else
            CALL WRITT(LU3,LBZY,BVECS(IBVEC))
#endif
            IF (IBVECX .NE. IBVEC) THEN
               DO 5000 I = 1,LBZY
                  BVECS(IBVECX-1+I) = BVECS(IBVEC-1+I)
 5000          CONTINUE
C    *         CALL DCOPY(LBZY,BVECS(IBVEC),1,BVECS(IBVECX),1)
            END IF
            IBVECX = IBVECX + LBZY
         END IF
         IBVEC = IBVEC + LBZY
 3300 CONTINUE
!
!     Set NBNEW to actual number of new trial vectors
!
      IF ((IBPREV-NBPREV).NE. NBNEW) THEN
         WRITE(LUPRI,'(A,I5,A,I5,A)')
     &      'RSPNOR:',NBNEW - IPREV + NBPREV,' out of',
     &       NBNEW,' new trial vectors linear dependent'
             call flshfo(lupri)
      END IF
      NBNEW = IBPREV - NBPREV
      RETURN
!
!     ERROR : Illegal IBTYP found:
!
 9100 CONTINUE
         WRITE (LUPRI,'(//A,/,(5X,20I3))')
     *      ' RSPNOR, illegal IBTYP found, dump of IBTYP:',
     *      (IBTYP(I),I=1,(NBPREV+NBNEW))
         CALL QTRACE(LUERR)
         CALL QUIT('RSPNOR ERROR, illegal IBTYP found.')
!
! *** End of subroutine RSPNOR
!
      END
C--
!  /* Deck rsportmcd */
      SUBROUTINE RSPORTmcd(BVECS,NBNEW,NBPREV,
     &                   IBTYP,THRLDP,OLDVEC,Xf,S2Xf,LU3)
!
! Written Nov 2009 by sonia, based on RSPORT by Hans Jorgen Aa. Jensen
! Purpose:
! For the special case of the MCD response equations at resonant frequency
! orthogonalize new b-vectors against all previous b-vectors
! and among themselves, and renormalize.
! Special care is taken NOT to orthonormalize vs btrial_1 and btrial_2
! which correspond to X_f and X_f,pair. The routine is also NOT used
! for btrial_1 and btrial_1paired.
! Each b-vector is in (Z, Y) form, the (Y, Z) vector to be obtained
! by switching first and second half. Each of the new b-vectors are
! orthogonalized both against (Z, Y) and (Y, Z) for each previous
! b-vector, except against btrial_1 and btrial_2.
! (Orthogonalization is performed twice if round-off is large,
! if larger than THRRND).
!
! Input:
!  BVECS,  non-orthogonal b-vectors
!  NBNEW,  number of b-vectors in BVECS
!  NBPREV, number of previous b-vectors on LU3
!  IBTYP,  b-vector type (CI or obital)
!  THRLDP, threshold for linear dependence
!
! Output:
!  BVECS,  orthogonal b-vectors (also written to LU3)
!  NBNEW,  number of b-vectors written to LU3
!
! Scratch:
!  OLDVEC(KZYVAR), scratch array for old b-vectors on LU3
!
#include <implicit.h>
#include <dummy.h>
      DIMENSION BVECS(*), IBTYP(*), OLDVEC(*), Xf(*), S2Xf(*)
C
#include <ibndxdef.h>
C
C    Change value of ZEQLY for problem with linear equation
C    95-02-08
C     PARAMETER (ZEQLY = 1.D-6, T1MIN = 1.D-8, THRRND = 1.D-4)
C
      PARAMETER (OVLMIN = 1.D-4, T1MIN = 1.D-8, THRRND = 1.D-4)
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C  INFRSP: RSPSUP
C  WRKRSP: KZCONF,KZWOPT,KSYMOP
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C     Local arrays:
C
      PARAMETER (MAXVEC = 120)
      DIMENSION JLDP(MAXVEC)
C
      WRITE (LUPRI,*) 'I AM INSIDE RSPORTmcd: NBNEW ', NBNEW
      CALL FLSHFO(lupri)
      IF (NBNEW .GT. MAXVEC) THEN
         WRITE (LUPRI,*) 'ERROR in RSPORTMCD: NBNEW2 .gt. MAXVEC'
         WRITE (LUPRI,*) 'NBNEW  =',NBNEW
         WRITE (LUPRI,*) 'MAXVEC =',MAXVEC
         WRITE (LUPRI,*) 'Reduce problem size or recompile RSPORTMCD '//
     *                  'with larger MAXVEC parameter'
         CALL QUIT('RSPORTMCD error: NBNEW .gt. MAXVEC parameter')
      END IF
!
! *** Average new trial vectors for supersymmetry
!
      IF ( RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         IBVEC = 1
         DO 2600 IRX = 1, NBNEW
            IF (IBTYP (NBPREV+IRX) .EQ. JBCNDX) THEN
!           -- Configuration trial vector; shift the starting point
               IBVEC = IBVEC + KZCONF
            ELSE IF (IBTYP(NBPREV+IRX) .EQ. JBONDX) THEN
!           -- Orbital trial vector; do the averaging
               CALL RSPAVE(BVECS(IBVEC),KZWOPT,2)
               IBVEC = IBVEC + KZYWOP
            END IF
 2600    CONTINUE
      END IF
!
!     For zero frequency linear response, the result will
!     be Z = Y for real perturbations and Z = -Y for imaginary
!     perturbations; trial vectors will have the same structure.
!
!     However, if Z = Y or Z = -Y then (Z Y) and (Y Z) are linear
!     dependent and we want instead (Z 0) and (0 Z) as trial vectors.
!
!check linear dependencies
!
      IBVEC = 1
      NCNEW = 0
      NONEW = 0
      DO 1400 IRX = 1,NBNEW
         JLDP(IRX) = 1
         IF (IBTYP(NBPREV+IRX) .EQ. JBCNDX) THEN
            NCNEW = NCNEW + 1
            IBVEC = IBVEC + KZCONF
         ELSE IF (IBTYP(NBPREV+IRX) .EQ. JBONDX) THEN
            NONEW = NONEW + 1
            T1    = DDOT(KZYWOP,BVECS(IBVEC),1,BVECS(IBVEC),1)
            JBVEC = IBVEC + KZWOPT
            OVLPI = D2*DDOT(KZWOPT,BVECS(IBVEC),1,BVECS(JBVEC),1)
!           test if Z + Y is zero
            TZPY  = T1 + OVLPI
            IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A/,(A,1P,D15.8))')
     *         ' ** TEST OF Z = Y OR Z = -Y:',
     *         ' (Z Y) norm squared =',T1,
     *         ' Z + Y norm squared =',TZPY
            IF (TZPY .LE. OVLMIN*T1) THEN
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,I4/A)')
     *            ' Z = -Y in trial vector no.',IRX,
     *            ' Y component removed.'
               CALL DZERO(BVECS(JBVEC),KZWOPT)
            call flshfo(lupri)
            ELSE IF (TZPY .GT. T1) THEN
!              test if Z - Y is zero;
!              only if TZPY = 2 * T1 can Z - Y be zero,
!              so TZPY .gt. T1 is a safe test.
               TZMY = T1 - OVLPI
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(A,1P,D15.8)')
     *            ' Z - Y norm squared =',TZMY
               IF (TZMY .LE. OVLMIN*T1) THEN
                  IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,I4/A)')
     *            ' Z = Y in trial vector no.',IRX,
     *            ' Y component removed.'
                  CALL DZERO(BVECS(JBVEC),KZWOPT)
            call flshfo(lupri)
               END IF
            END IF
            IBVEC = IBVEC + KZYWOP
         ELSE
            GO TO 9100
         END IF
 1400 CONTINUE
C
      IROUND = 0
      ITURN  = 0
 1500 ITURN  = ITURN + 1
!
!        Orthogonalize new b-vectors agains previous b-vectors
!        (both (Z, Y) and (Y, Z))
!
         CALL REWSPL(LU3)
#if defined (VAR_SPLITFILES)
         JRSP3 = 0
#endif
         DO 2000 K = 1,NBPREV
            IBTYPK = IBTYP(K)
            IF (IBTYPK .EQ. JBCNDX) THEN
               IF (NCNEW .EQ. 0) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LU3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
                  READ (LU3)
#endif
                  GO TO 2000
!        v------------------
               END IF
               LBZ = KZCONF
               LBZY= KZCONF
               LBSKIP = 2*KZWOPT
            ELSE IF (IBTYPK .EQ. JBONDX) THEN
               IF (NONEW .EQ. 0) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LU3,LRSP3,DUMMY,JRSP3,.FALSE.)
#else
                  READ (LU3)
#endif
                  GO TO 2000
!        v------------------
               END IF
               LBZ = KZWOPT
               LBZY= 2*KZWOPT
               LBSKIP = KZCONF
            ELSE
               GO TO 9100
            END IF
#if defined (VAR_SPLITFILES)
            CALL READST(LU3,LRSP3,OLDVEC,JRSP3,.TRUE.)
            IF (LRSP3 .NE. LBZY) CALL QUIT('Inconsistent record '//
     &           'length on unit LU(RSP)3 in RSPORT')
#else
            CALL READT(LU3,LBZY,OLDVEC)
#endif
            IBVEC = 1
            DO 1900 IRX = 1,NBNEW
               IF (IBTYP(NBPREV+IRX) .NE. IBTYPK) THEN
                  IBVEC = IBVEC + LBSKIP
               ELSE IF (JLDP(IRX).EQ.0) THEN
                  IBVEC = IBVEC + LBZY
               ELSE
                  if (k.eq.1) then
                     !double project, k==1 identifies b_1=Xf and b_2=X_-f
                     !<s2b1|b3>
                     tt = DDOT(LBZY,S2Xf,1,BVECS(IBVEC),1)
                     CALL DAXPY(LBZY,-TT,Xf,1,BVECS(IBVEC),1)
                     JBVEC = IBVEC + LBZ
                     tt1= DDOT(LBZ,S2Xf(1+LBZ),1,BVECS(IBVEC),1)
     &                  + DDOT(LBZ,S2Xf(1),1,BVECS(JBVEC),1)
                     CALL DAXPY(LBZ,-TT1,Xf(1+LBZ),1,BVECS(IBVEC),1)
                     CALL DAXPY(LBZ,-TT1,Xf(1),1,BVECS(JBVEC),1)
                  else
!                   Orthogonalize to (Z Y)_old
                    TT = -DDOT(LBZY,   OLDVEC,1,BVECS(IBVEC),1)
                    CALL DAXPY(LBZY,TT,OLDVEC,1,BVECS(IBVEC),1)
!                   Orthogonalize to (Y Z)_old
                    IF ( IBTYPK .EQ. JBONDX ) THEN
                       JBVEC = IBVEC + LBZ
                       TT = -DDOT(LBZ,   OLDVEC,       1,BVECS(JBVEC),1)
     *                      -DDOT(LBZ,   OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
                       CALL DAXPY(LBZ,TT,OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
                       CALL DAXPY(LBZ,TT,OLDVEC,       1,BVECS(JBVEC),1)
                    END IF
                  end if
                  IBVEC = IBVEC + LBZY
               END IF
 1900       CONTINUE
 2000    CONTINUE
!
!        Orthogonalize new vectors against each other
!
         IBVEC  = 1
         DO 2400 IRX = 1,NBNEW
            IBTYPI = IBTYP(NBPREV+IRX)
            IF (IBTYPI .EQ. JBCNDX) THEN
               LBZ = KZCONF
               LBZY= KZCONF
               LBSKIP = 2*KZWOPT
            ELSE
               LBZ = KZWOPT
               LBZY= 2*KZWOPT
               LBSKIP = KZCONF
            END IF
         IF (JLDP(IRX) .NE. 0) THEN
            JBVEC = 1
            DO 2300 JRX = 1,(IRX-1)
            IF (IBTYP(NBPREV+JRX) .NE. IBTYPI) THEN
               JBVEC = JBVEC + LBSKIP
            ELSE IF (JLDP(JRX).EQ.0) THEN
               JBVEC = JBVEC + LBZY
            ELSE
               T1 =  DDOT(LBZY,BVECS(JBVEC),1,BVECS(JBVEC),1)
               TT = -DDOT(LBZY,BVECS(JBVEC),1,BVECS(IBVEC),1)/T1
               CALL DAXPY(LBZY,TT,BVECS(JBVEC),1,BVECS(IBVEC),1)
               IF (IBTYPI.EQ.JBONDX) THEN
                  TT=(-DDOT(LBZ,BVECS(JBVEC+LBZ),1,BVECS(IBVEC),1)
     *                -DDOT(LBZ,BVECS(JBVEC),1,BVECS(IBVEC+LBZ),1))/T1
                  CALL DAXPY(LBZ,TT,BVECS(JBVEC+LBZ),1,BVECS(IBVEC),1)
                  CALL DAXPY(LBZ,TT,BVECS(JBVEC),1,BVECS(IBVEC+LBZ),1)
               END IF
               JBVEC = JBVEC + LBZY
            END IF
 2300       CONTINUE
         END IF
!
! NORMALIZE VECTOR NUMBER IRX
!
         IF (JLDP(IRX) .NE. 0) THEN
            TT = DDOT(LBZY,BVECS(IBVEC),1,BVECS(IBVEC),1)
            IF (TT .LE. THRLDP) THEN
               IF (IPRRSP.GT.2) WRITE (LUPRI,3250) IRX,JLDP(IRX),TT
               JLDP(IRX) = 0
            ELSE IF (TT .LT. THRRND) THEN
               IF (ITURN .EQ. 1) THEN
                  IROUND = IROUND + 1
               ELSE
                  IF (IPRRSP.GT.10) WRITE (LUPRI,3255) IRX,JLDP(IRX),TT
                  JLDP(IRX) = 0
               END IF
            END IF
            IF (JLDP(IRX) .NE. 0) THEN
               IF (IPRRSP .GE. 10) WRITE (LUPRI,'(/A,2I6,1P,D12.5)')
     &    ' RSPORTMCD: b-vector, root no., norm**2:',IRX,JLDP(IRX),TT
               call flshfo(lupri)
               IF (TT .LT. T1MIN) THEN
                  TT = D1 / SQRT(TT)
                  CALL DSCAL(LBZY,TT,BVECS(IBVEC),1)
                  TT = DDOT(LBZY,BVECS(IBVEC),1,BVECS(IBVEC),1)
               END IF
               TT = D1 / SQRT(TT)
               CALL DSCAL(LBZY,TT,BVECS(IBVEC),1)
            END IF
         END IF
!
! Perform symmetric orthonormalization of (Z Y) and (Y Z) pair
! for vector number IRX
!
!   -1/2       ( C1   C2 )               (  1     OVLPI )
!  S      =   (           ) where S  =  (                )
!              ( C2   C1 )               ( OVLPI     1  )
!
         IF ((JLDP(IRX) .NE. 0).AND.(IBTYPI.EQ.JBONDX)) THEN
            JBVEC = IBVEC + LBZ
            JTURN = 0
 2350       JTURN = JTURN + 1
            OVLPI = D2*DDOT(LBZ,BVECS(IBVEC),1,BVECS(JBVEC),1)
            IF (IPRRSP .GE. 10) THEN
               X1 = DNORM2(LBZY,BVECS(IBVEC),1)
               WRITE (LUPRI,'(A,1P,D14.6,D19.10)')
     *    ' RSPORTMCD, norm and <ZY/YZ> overlap before S-1/2',X1,OVLPI
            call flshfo(lupri)
            END IF
            X2 = ABS(OVLPI)
            IF (ABS(D1-X2) .LT. OVLMIN) THEN
               IF (IPRRSP.GT.10) WRITE (LUPRI,3260) IRX,JLDP(IRX),OVLPI
               call flshfo(lupri)
               JLDP(IRX) = 0
               GO TO 2390
            END IF
            X1 = D1+OVLPI
            X2 = D1-OVLPI
            IF (ABS(X1 - X2) .GT. THRLDP) THEN
               X1 = DP5 / SQRT(X1)
               X2 = DP5 / SQRT(X2)
               C1 = X1 + X2
               C2 = X1 - X2
               !note here oldvec is simply scratch area
               CALL DCOPY(LBZY,BVECS(IBVEC),1,OLDVEC,1)
               CALL DSCAL(LBZY,C1,BVECS(IBVEC),1)
               CALL DAXPY(LBZ,C2,OLDVEC,1,BVECS(JBVEC),1)
               CALL DAXPY(LBZ,C2,OLDVEC(1+LBZ),1,BVECS(IBVEC),1)
            IF (JTURN .EQ. 1) GO TO 2350
!           ...930712-hjaaj: an output with LBZY only = 480 showed
!              that OVLPI ca. 1.0D-14 after first S-1/2;
!              but ca. 1.0D-16 after second S-1/2
               !IF (IPRRSP .GE. 22) THEN
               IF (IPRRSP .GE. 10) THEN
                  X1 = DNORM2(LBZY,BVECS(IBVEC),1)
                  OVLPI= D2*DDOT(LBZ,BVECS(IBVEC),1,BVECS(JBVEC),1)
                  WRITE (LUPRI,'(A,1P,D14.6,D19.10)')
     *      ' RSPORTMCD, norm and <ZY/YZ> overlap after  S-1/2',
     *            X1,OVLPI
               END IF
            END IF
         END IF
C
C
 2390    IBVEC = IBVEC + LBZY
 2400    CONTINUE
         IF (IROUND.GT.0) THEN
            IROUND = 0
            IF (ITURN.EQ.1) GO TO 1500
C     ^-------------------------------
            CALL QUIT('RSPORTMCD sw-error; IROUND.gt.0 and ITURN.gt.1')
         END IF
 3250 FORMAT(/' (RSPORTMCD) b-vector no.',I3,' corresponding to ',
     *        'root no.',I3/,14X,'is removed because of linear ',
     *        'dependence; ',1P,D12.5)
 3255 FORMAT(/' (RSPORTMCD) b-vector no.',I3,' corresponding to root ',
     *        'no.',I3,' is removed',/,14X,'because of small',
     *        'norm**2 in 2. Gram-Schmidt orthonormalization:',1P,D12.5)
 3260 FORMAT(/' (RSPORTMCD) b-vector no.',I3,' corresponding to ',
     *        'root no.',I3,' is removed',/,14X,'because of lin.dep.,',
     *        '<Z Y|Y Z> overlap =',1P,D14.6)
C
C Save new vector together with old ones on LU3
C
      IBVECX = 1
      IBVEC  = 1
      IBPREV = NBPREV
      DO 3300 IRX = 1,NBNEW
         IBTYPI = IBTYP(NBPREV+IRX)
         IBTYP(NBPREV+IRX) = 0
         IF (IBTYPI .EQ. JBCNDX) THEN
            LBZY= KZCONF
         ELSE
            LBZY= 2*KZWOPT
         END IF
         IF (JLDP(IRX) .NE. 0) THEN
            IBPREV = IBPREV + 1
            IBTYP(IBPREV) = IBTYPI
C
#if defined (VAR_SPLITFILES)
            CALL WRITST(LU3,LBZY,BVECS(IBVEC),JRSP3)
#else
            CALL WRITT(LU3,LBZY,BVECS(IBVEC))
#endif
            IF (IBVECX .NE. IBVEC) THEN
               DO 5000 I = 1,LBZY
                  BVECS(IBVECX-1+I) = BVECS(IBVEC-1+I)
 5000          CONTINUE
C    *         CALL DCOPY(LBZY,BVECS(IBVEC),1,BVECS(IBVECX),1)
            END IF
            IBVECX = IBVECX + LBZY
         END IF
         IBVEC = IBVEC + LBZY
 3300 CONTINUE
!
!     Set NBNEW to actual number of new trial vectors
!
      IF (IBPREV-NBPREV .NE. NBNEW) THEN
         WRITE(LUPRI,'(A,I5,A,I5,A)')
     &      'RSPORTMCD:',NBNEW - IPREV + NBPREV,' out of',
     &       NBNEW,' new trial vectors linear dependent'
      END IF
      NBNEW = IBPREV - NBPREV
!
      RETURN
!
!     ERROR : Illegal IBTYP found:
!
 9100 CONTINUE
         WRITE (LUPRI,'(//A,/,(5X,20I3))')
     *      ' RSPORTMCD, illegal IBTYP found, dump of IBTYP:',
     *      (IBTYP(I),I=1,(NBPREV+NBNEW))
         CALL QTRACE(LUERR)
         CALL QUIT('RSPORTMCD ERROR, illegal IBTYP found.')
!
! *** End of subroutine RSPORTMCD
!
      END

!  /* Deck rspnexmcd */
      SUBROUTINE RSPNEXmcd(LINEQ,LMAXIT,NSIM,IBTYP,EIVAL,TEST,EIVEC,GD,
     *                  CMO,UDV,PVX,FC,FV,FCAC,H2AC,XINDX,
     *                  BVECS,LWRK,Xf,S2Xf)
!
! PURPOSE: 1)CONSTRUCT RESIDUAL (E(2)-W(I)*S(2))*X(I)
!            AND MODIFIED GRADIENT  (E(2)-W(I)*S(2))*X(I)(csf part)
!            FOR KEXSIM EIGENVECTORS X(I) OF REDUCED RSP PROBLEM
!          2)TEST FOR CONVERGENCE OF KEXCNV EIGENVECTORS,
!            CONVERGENCE CRITERIUM:
!            //(E(2)-W(I)*S(2))*X(I)// .LE. THCRSP * //X(I)//
!          3)GENERATE NEW TRIAL VECTORS
!             A) CSF TRIAL VECTORS, USE GENERALIZED DAVIDSON ALGORITHM
!                CSF TRIAL VECTORS CONTAIN KZCONF ELEMENTS
!             B) ORBITAL TRIAL VECTORS
!                IF OPTORB.EQ.TRUE, USE OPTIMAL ORBITAL ALGORITHM ELSE
!                USE GENERALIZED DAVIDSON ALGORITHM
!
! Parameters: LMAXIT true: max iterations reached, do not
!                          calculate new trial vectors
! PJ NOV-1984
! LAST REV 12 NOV 1984, 22 APR 1988, 28 Oct 93 hjaaj
! SONIA'S SPECIAL VERSION FOR MCD projection, Jan. 2010
!
#include <implicit.h>
#include <mxcent.h>
#include <dummy.h>
      PARAMETER ( MAXVEC = 120 )
      DIMENSION THCORP(MAXVEC),EOVAL(MAXVEC)
      DIMENSION IBTYP(*),EIVAL(*),EIVEC(KZYRED,*)
      DIMENSION GD(*),BVECS(*)
      DIMENSION CMO(*),UDV(*),PVX(*),FC(*),FV(*),FCAC(*)
      DIMENSION XINDX(*),H2AC(*)
      DIMENSION Xf(*),S2Xf(*)
#include <ibndxdef.h>
      PARAMETER ( D0=0.0D0 , D1=1.0D0 , DM1=-1.0D0 )
      PARAMETER ( REDFAC = 0.3D0 , DTEST=1.0D-4  )
C
C Used from common blocks:
C  INFRSP: ?
C  WRKRSP: KZVAR,KZYVAR,KEXSIM,KEXCNV,
C          LURSP5
C  INFOPT: EACTIV
C
#include <priunit.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <inftap.h>
#include <infopt.h>
#include <infdim.h>
#include <abares.h>
#include <thrldp.h>
C
      LOGICAL LINEQ,OLSEN1,LMAXIT
      DIMENSION KEX(MAXVEC), TEST(*)
C
      IF (KEXSIM .GT. MAXVEC) THEN
         WRITE (LUPRI,*) 'ERROR in RSPNEXMCD: KEXSIM .gt. MAXVEC'
         WRITE (LUPRI,*) 'KEXSIM =',KEXSIM
         WRITE (LUPRI,*) 'MAXVEC =',MAXVEC
         WRITE (LUPRI,*) 'Reduce problem size or recompile RSPNEXMCD '//
     *                  'with larger MAXVEC parameter'
         CALL QUIT('RSPNEXMCD error: KEXSIM .gt. MAXVEC parameter')
      END IF
C
C     KEX(*) is used for defining which trial vectors are
C     configuration type (=1) and which are orbital type (=0).
C     920904-hjaaj: KEX(i) = -1 means this root is converged.
C     TEST(*) is used for saving residuals for final print.
C
      NTEST = 0
      IBOFF = 0
      NTCONV = 0
      DO 4000 ISIMC = 1,KEXSIM,NSIM
         NBX  = MIN(NSIM,(KEXSIM+1-ISIMC))
         NFIN = 0
C
C RSPEVE constructs solution vectors from reduced space
C CONSTRUCT -W(I) * S[2]*X(I) IN FIRST NBX ELEMENTS OF
C BVECS. THE PART WHICH IS USED IN THE OPTIMAL ORBITAL
C TRIAL VECTOR ALGORITHM WHICH IS STORED IN THE CONSECUTIVE
C NBX*KZYWOP ELEMENTS
C
C First check if all frequencies zero (W(i)=0, all i) because then
C we may skip the -W(I) * S[2] term completely
C If SOPPA we will construct the 2p-2h diagonal times the 2p-2h trial
C vector explicitly so we need the solution vectors
C
         IF (SOPPA) GOTO 25
         DO 24 I = ISIMC, ISIMC-1+NBX
            IF ( EIVAL(I) .NE. D0 ) GO TO 25
   24    CONTINUE
         LEN = NBX*KZYVAR
         IF (OPTORB) LEN = LEN + NBX*KZYWOP
         CALL DZERO(BVECS,LEN)
         GO TO 200
C
   25    CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS(1+KZYVAR),
     *               BVECS(1),NBX,(ISIMC-1))
         IF (IPRRSP.GT.10) THEN
            WRITE (LUPRI,*) ' RSPNEXMCD: Solution vectors no.',
     *         ISIMC,' to',(ISIMC-1+NBX)
            CALL OUTPUT(BVECS(1+KZYVAR),1,KZYVAR,1,NBX,
     *                  KZYVAR,NBX,-1,LUPRI)
         ENDIF
         KRES   = 1
         KINFIX = KRES + (NBX+1)*KZYVAR
         KBORB  = KINFIX + NBX*KZYWOP
         KDIA   = KBORB + KZYWOP
C
C Make room for the D-matrix in SOPPA
C
         IF (SOPPA) THEN
            KWRK1 = KDIA + KZCONF
         ELSE
            KWRK1 = KDIA
         ENDIF
         LWRK1 = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPNEXMCD 1',KWRK1-1,LWRK)
         IF (OPTORB) CALL DZERO(BVECS(KINFIX),NBX*KZYWOP)
C
C  Read in the D-matrix in SOPPA from LURSP4
C
         IF (SOPPA) THEN
            CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
            JRSP4 = 0
            CALL READST(LURSP4,LRSP4,BVECS(KDIA),JRSP4,.TRUE.)
            IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP4 in RSPNEXMCD')
#else
            CALL READT(LURSP4,KZCONF,BVECS(KDIA))
#endif
         ENDIF
         DO 50 IBX = 1,NBX
            IRESO = (IBX-1)*KZYVAR + 1
            IFIXOF= KINFIX + (IBX-1)*KZYWOP
            WIBX  = -EIVAL(ISIMC-1+IBX)
            CALL DZERO(BVECS(IRESO),KZYVAR)
            IF (WIBX .EQ. D0) GO TO 49
C           ... skip S[2] term if zero frequency W(i) = 0
C           We still need to do the D-matrix term in SOPPA
            IF ( KZCONF .GT. 0 ) THEN
               CALL RSPSLI(1,0,BVECS(IRESO+KZYVAR+KZVAR),DUMMY,
     *                  UDV,BVECS(IRESO),XINDX,BVECS(KWRK1),LWRK1)
               CALL DSWAP(KZVAR,BVECS(IRESO),1,BVECS(IRESO+KZVAR),1)
               CALL DSCAL(KZYVAR,DM1,BVECS(IRESO),1)
               CALL RSPSLI(1,0,BVECS(IRESO+KZYVAR),DUMMY,
     *                     UDV,BVECS(IRESO),XINDX,BVECS(KWRK1),LWRK1)
               IF ( OPTORB ) THEN
                  CALL DAXPY(KZWOPT,WIBX,BVECS(IRESO+KZCONF),1,
     *                       BVECS(IFIXOF),1)
                  CALL DZERO(BVECS(IRESO+KZCONF),KZWOPT)
                  CALL DAXPY(KZWOPT,WIBX,BVECS(IRESO+KZVAR+KZCONF),1,
     *                       BVECS(IFIXOF+KZWOPT),1)
                  CALL DZERO(BVECS(IRESO+KZVAR+KZCONF),KZWOPT)
               END IF
            END IF
            IF ( KZWOPT.GT.0 ) THEN
               CALL DCOPY(KZWOPT,BVECS(IRESO+KZYVAR+KZCONF),1,
     *                    BVECS(KBORB),1)
               CALL DCOPY(KZWOPT,BVECS(IRESO+KZYVAR+KZVAR+KZCONF),1,
     *                    BVECS(KBORB+KZWOPT),1)
               CALL RSPSLI(0,1,DUMMY,BVECS(KBORB),
     *                     UDV,BVECS(IRESO),XINDX,BVECS(KWRK1),LWRK1)
            ENDIF
            CALL DSCAL(KZYVAR,WIBX,BVECS(IRESO),1)
C
 49         IF (SOPPA) THEN
               DO I=0,KZCONF-1
                  BVECS(IRESO+I) = BVECS(IRESO+I) +
     *                             BVECS(KDIA+I) * BVECS(IRESO+KZYVAR+I)
                  BVECS(IRESO+KZVAR+I) = BVECS(IRESO+KZVAR+I) +
     *                       BVECS(KDIA+I) * BVECS(IRESO+KZVAR+KZYVAR+I)
               ENDDO
            ENDIF
 50      CONTINUE
         IF (OPTORB) THEN
            NTT = NBX*KZYWOP
            DO 52 II = 1,NTT
               BVECS(NBX*KZYVAR+II) = BVECS(KINFIX-1+II)
 52         CONTINUE
C           CALL DCOPY(NBX*KZYWOP,BVECS(KINFIX),1,
C    *                            BVECS(1+NBX*KZYVAR),1)
         END IF
C
C ALLOCATE WORK SPACE
C
  200    KBVECS = 1
         KFIXOR = KBVECS + NBX*KZYVAR
         KWRK1  = KFIXOR + NBX*KZYWOP
         LWRK1  = LWRK   - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPNEXMCD 2',KWRK1-1,LWRK)
C
C CONSTRUCT  E[2]*X(I) WHERE X(I) IS THE I'TH
C EIGENVECTOR AND EIGENVALUE RESPECTIVELY OF THE REDUCED
C RSP PROBLEM
C
         CALL REWSPL(LURSP5)
#if !defined (VAR_SPLITFILES)
         IF (KOFFTY.EQ.1) READ(LURSP5)
#else
         JRSP5 = 0
         IF (KOFFTY .EQ. 1)
     &        CALL READST(LURSP5,LRSP5,DUMMY,JRSP5,.FALSE.)
#endif
         DO 600 K = 1,KZRED
C
C   If SOPPA and the b-vector is 2p2h then only read the p-h part
C   from LURSP5 since the 2p-2h is constructed explicitly above
C
            IF (SOPPA) THEN
               IF (IBTYP(K+KOFFTY) .EQ. JBCNDX) THEN
#if defined (VAR_SPLITFILES)
                  CALL READST(LURSP5,LRSP5,BVECS(KWRK1),JRSP5,.TRUE.)
                  IF (LRSP5 .NE. KZYWOP) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP5 in RSPNEXMCD')
               ELSE
                  CALL READST(LURSP5,LRSP5,BVECS(KWRK1),JRSP5,.TRUE.)
                  IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP5 in RSPNEXMCD (2)')
#else
                  CALL READT(LURSP5,KZYWOP,BVECS(KWRK1))
               ELSE
                  CALL READT(LURSP5,KZYVAR,BVECS(KWRK1))
#endif
               ENDIF
            ELSE
#if defined (VAR_SPLITFILES)
               CALL READST(LURSP5,LRSP5,BVECS(KWRK1),JRSP5,.TRUE.)
                  IF (LRSP5 .NE. KZYVAR) CALL QUIT('Inconsistent '//
     &                 'record length on unit LURSP5 in RSPNEXMCD (3)')
#else
               CALL READT(LURSP5,KZYVAR,BVECS(KWRK1))
#endif
            ENDIF
            DO 700 JR = 1,NBX
               JROFF  = (JR-1)*KZYVAR
               JROOTJ = IBOFF+JR
               FAC1   = EIVEC(2*K-1,JROOTJ)
               FAC2   = EIVEC(2*K,JROOTJ)
               IF (SOPPA .AND. IBTYP(K+KOFFTY).EQ.JBCNDX
     *                   .AND. .NOT.OPTORB) THEN
                  CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1),1,
     *                       BVECS(1+JROFF+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1+KZWOPT),1,
     *                       BVECS(1+JROFF+KZVAR+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1+KZWOPT),1,
     *                       BVECS(1+JROFF+KZCONF),1)
                  CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1),1,
     *                       BVECS(1+JROFF+KZVAR+KZCONF),1)
               ELSE IF ( IBTYP(K+KOFFTY) .EQ. JBONDX
     *                              .OR. .NOT.OPTORB) THEN
                  CALL DAXPY(KZYVAR,FAC1,BVECS(KWRK1),1,
     *                       BVECS(1+JROFF),1)
                  CALL DAXPY(KZVAR,FAC2,BVECS(KWRK1+KZVAR),1,
     *                                       BVECS(1+JROFF),1)
                  CALL DAXPY(KZVAR,FAC2,BVECS(KWRK1),1,
     *                                       BVECS(1+KZVAR+JROFF),1)
               ELSE
C
C CSF PART OF LINEAR TRANSFORMED VECTOR
C
                  IF (.NOT. SOPPA) THEN
                     CALL DAXPY(KZCONF,FAC1,BVECS(KWRK1),1,
     *                          BVECS(1+JROFF),1)
                     CALL DAXPY(KZCONF,FAC1,BVECS(KWRK1+KZVAR),1,
     *                                      BVECS(1+KZVAR+JROFF),1)
                     CALL DAXPY(KZCONF,FAC2,BVECS(KWRK1),1,
     *                                      BVECS(1+KZVAR+JROFF),1)
                     CALL DAXPY(KZCONF,FAC2,BVECS(KWRK1+KZVAR),1,
     *                                      BVECS(1+JROFF),1)
                  ENDIF
C
C                 ORBITAL PART OF LINEAR TRANSFORMED VECTOR
C
                  IF (KZWOPT.GT.0) THEN
                     KKAPZY = KFIXOR + (JR-1)*KZYWOP
                     KKAPYZ = KKAPZY + KZWOPT
                     CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1+KZCONF),1,
     *                                      BVECS(KKAPZY),1)
                     CALL DAXPY(KZWOPT,FAC1,BVECS(KWRK1+KZVAR+KZCONF),
     *                                      1,BVECS(KKAPYZ),1)
                     CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1+KZCONF),1,
     *                                      BVECS(KKAPYZ),1)
                     CALL DAXPY(KZWOPT,FAC2,BVECS(KWRK1+KZVAR+KZCONF),
     *                                      1,BVECS(KKAPZY),1)
                  ENDIF
               ENDIF
 700        CONTINUE
 600     CONTINUE
C
         IF (OPTORB) THEN
            KFIXOF = KFIXOR
            DO 1000 JR = 1,NBX
               JRESOF = (JR-1)*KZYVAR + 1
               CALL DAXPY(KZWOPT,D1,BVECS(KFIXOF),1,
     *                              BVECS(JRESOF+KZCONF),1)
               CALL DAXPY(KZWOPT,D1,BVECS(KFIXOF+KZWOPT),1,
     *                              BVECS(JRESOF+KZVAR+KZCONF),1)
               KFIXOF = KFIXOF + KZYWOP
 1000       CONTINUE
         END IF
         IF (LINEQ) THEN
            KFIXOF = KFIXOR
            DO 1050 JR=1,NBX
               JROFF = (JR-1)*KZYVAR
               CALL DAXPY(KZYVAR,DM1,GD(1),1,
     *                                  BVECS(1+JROFF),1)
               IF (OPTORB) THEN
                  CALL DAXPY(KZWOPT,DM1,GD(1+KZCONF),1,
     *                        BVECS(KFIXOF),1)
                  CALL DAXPY(KZWOPT,DM1,GD(1+KZVAR+KZCONF),1,
     *                        BVECS(KFIXOF+KZWOPT),1)
               ENDIF
            KFIXOF = KFIXOF + KZYWOP
 1050       CONTINUE
         ENDIF
         IF (IPRRSP.GT.101) THEN
            WRITE (LUPRI,*) ' RSPNEXMCD: (Z,Y) residual vectors no.',
     *         ISIMC,' to',(ISIMC-1+NBX)
            CALL OUTPUT(BVECS,1,KZYVAR,1,2*NBX,
     *                  KZYVAR,2*NBX,-1,LUPRI)
         ENDIF
!
!        THE RESIDUAL IS NOW CONSTRUCTED
!
!        TEST FOR CONVERGENCE OF THE EIGENVECTORS
!        AND FORM NEW TRIAL VECTORS
!
         CALL PHPINI(LPHPMX,KZCONF,KZWOPT,MAXPHP)
         KDIAE   = KWRK1
         KDIASO  = KDIAE  + LPHPMX
         IF (OLSEN) THEN
            KCVEC = KDIASO + KZWOPT
            KTOT  = KCVEC  + KZCONF
         ELSE
            KCVEC = KDIASO
            KTOT  = KDIASO + KZWOPT
         END IF
         LTOT  = LWRK  - KTOT
         IF (LTOT.LT.0) CALL ERRWRK('RSPNEXMCD 3',KTOT-1,LWRK)
         ECORE  = -EACTIV
         IF (PHPRES) THEN
            IPWAY = 2
         ELSE IF (KZCONF.GT.0) THEN
            CALL PHPDSK(BVECS(KDIAE),LPHPMX,.FALSE.)
         ELSE
            CALL RSPEDG(BVECS(KDIAE))
         END IF
         CALL RSPSOD(BVECS(KDIASO))
         NCSIM = 0
         NOSIM = 0
         NOTCNV = 0
         KCOFF = 0
         IF (IPRRSP.GE.3 .AND. IPRRSP.LE.5) THEN
            IF (LINEQ) THEN
               WRITE(LUPRI,'(/6X,A)')
     &' No.  Residual tot.,    conf., and orb.    Bnorm      Frequency '
            ELSE
               WRITE(LUPRI,'(/6X,A)')
     &'Root  Residual tot.,    conf., and orb.    Bnorm      Eigenvalue'
            END IF
            WRITE(LUPRI,'(6X,A)')
     &'----------------------------------------------------------------'
         END IF
         DO 2000 JR = 1,NBX
            JROFF = (JR-1)*KZYVAR
            JROOTJ = IBOFF+JR
            QONORM = DDOT(KZWOPT,BVECS(1+KZCONF+JROFF),1,
     *                               BVECS(1+KZCONF+JROFF),1)
     *             + DDOT(KZWOPT,BVECS(1+KZVAR+KZCONF+JROFF),1,
     *                               BVECS(1+KZVAR+KZCONF+JROFF),1)
            QCZNOR = DDOT(KZCONF,BVECS(1+JROFF),1,
     *                               BVECS(1+JROFF),1)
            QCYNOR = DDOT(KZCONF,BVECS(1+KZVAR+JROFF),1,
     *                               BVECS(1+KZVAR+JROFF),1)
            QCNORM = QCZNOR + QCYNOR
            QNORM  = SQRT(QONORM+QCNORM)
            QCNORM = SQRT(QCNORM)
            QONORM = SQRT(QONORM)
            QCZNOR = SQRT(QCZNOR)
            QCYNOR = SQRT(QCYNOR)
            BNORM  = DDOT(KZYRED,EIVEC(1,JROOTJ),1,EIVEC(1,JROOTJ),1)
            BNORM  = SQRT(BNORM)
            IF (IPRRSP.GE.6) THEN
               WRITE(LUPRI,'(/A,I3,1P,3(A,D12.5),/T10,2(A,D12.5))')
     *           ' ROOT:',JROOTJ,' RESIDUAL TOT:',QNORM,
     *           ' CONF:',QCNORM,' ORB:',QONORM,
     *           ' BNORM:',BNORM,' EIVAL:',EIVAL(JROOTJ)
               IF (QCNORM .NE. D0) WRITE(LUPRI,'(1P,T10,2(A,D12.5))')
     *           ' RESIDUAL CONF Z-COMP',QCZNOR,'    Y-COMP',QCYNOR
            ELSE IF (IPRRSP.GE.3) THEN
               WRITE(LUPRI,'(I10,1P,D15.5,3D10.2,D15.5)')
     *           JROOTJ,QNORM,QCNORM,QONORM,BNORM,EIVAL(JROOTJ)
            END IF
!HJAAJ      BNORM = MIN(BNORM,D1)
!           ... 880309: ask for absolute convergence if BNORM
!               is greater than 1, otherwise relative convergence.
!               930910-hjaaj: NO, gives problems with number of digits
!               when e.g. BNORM = 1.0D6 for a FC term (remember E[2]
!               is asymmetric of the order of the mcscf gradient).
            QBTEST = QNORM/BNORM
            IF (JROOTJ.LE.MAXVEC) THEN
               TEST(JROOTJ) = QBTEST
            END IF
            IF (QBTEST .LE. THCRSP) THEN
C           ... this root is converged
               KEX(JR) = -1
            ELSE
C           ... this root is not converged
               IF (JROOTJ .LE. KEXCNV) NOTCNV = NOTCNV + 1
               EXCITA  =  EIVAL(JROOTJ)
              IF (.NOT.LMAXIT) THEN
               IF (QCNORM.GT.QONORM) THEN
                  NCSIM = NCSIM + 1
                  KEX(IBOFF+JR) = 1
!
!              USE DAVIDSON ALGORITHM TO FORM NEW TRIAL VECTORS
!
                  IF ( QCZNOR .GT. QCYNOR ) THEN
                     IOFF = 0
                     EXITA =   EXCITA
                  ELSE
                     IOFF = KZVAR
                     EXITA =  -EXCITA
                  END IF
                  IF ((.NOT.LINEQ).AND.OLSEN) THEN
                     OLSEN1 = .TRUE.
                  ELSE
                     OLSEN1 = .FALSE.
                  ENDIF
                  IF (OLSEN1) THEN
                     CALL RSPCVE(IOFF,JROOTJ,IBTYP,EIVAL,EIVEC,
     *                  BVECS(KCVEC),BVECS(KTOT),LTOT)
                  END IF
                  IF (PHPRES) THEN
                     DO 1003 J = 1,KZCONF
                        BVECS(KDIAE-1+J) = BVECS(JROFF+IOFF+J) /
     *                                     BVECS(KDIAE-1+J)
 1003                CONTINUE
                     CALL PHPGET(KSYMST,KZCONF,XINDX,FCAC,
     *                           H2AC,BVECS(KDIAE),
     *                           ECORE,IPWAY,IPRRSP,BVECS(KTOT),LTOT)
                     CALL RSPEDG(BVECS(KDIAE))
                  END IF
                  CALL NEXCI(OLSEN1,EXITA,KZCONF,BVECS(KCOFF+1),
     *            BVECS(KCVEC),BVECS(JROFF+IOFF+1),BVECS(KDIAE),IPRRSP,
     *            BVECS(KTOT),LTOT)
!                 CALL NEXCI(OLSEN,ENER,NCVAR,D,XVEC,RES,
!    &                 DIAG,IPRPHP,WRK,LWRK)
                  KCOFF = KCOFF + KZCONF
               ELSE
                  NOSIM = NOSIM + 1
                  KEX(IBOFF+JR) = 0
                  IF (PHPRES) CALL RSPEDG(BVECS(KDIAE))
                  IF (OPTORB) THEN
!
!        USE OPTIMAL ORBITAL ALGORITHM TO FORM NEW TRIAL VECTORS
!        WILL BE DONE IN LATER SECTION
!
                     IF (QCNORM .GT. THCRSP) THEN
                        THCORP(NOSIM) = REDFAC * QCNORM/BNORM
                     ELSE
!                       ... CI part is converged, attempt
!                           quadratic convergence in orbital part.
                        THCORP(NOSIM) = REDFAC *
     *                     MAX(THCRSP,QBTEST*QBTEST)
                     END IF
                     EOVAL(NOSIM)  = EIVAL(JROOTJ)
                  ELSE
!
!        USE DAVIDSON ALGORITHM TO FORM NEW TRIAL VECTORS
!
                     KSJR = KFIXOR - 1 + (JR-1)*KZYWOP
                     DO 1500 I=1,KZWOPT
                        IZ = I+KZCONF
                        ZDIA = BVECS(KDIAE-1+IZ) -
     &                         EXCITA*BVECS(KDIASO-1+I)
                        YDIA = BVECS(KDIAE-1+IZ) +
     &                         EXCITA*BVECS(KDIASO-1+I)
                        IF (ABS(ZDIA).LE.DTEST) THEN
                           BVECS(KSJR+I)= BVECS(IZ+JROFF) /
     &                     SIGN(DTEST,ZDIA)
                        ELSE
                           BVECS(KSJR+I) = BVECS(IZ+JROFF) / ZDIA
                        ENDIF
                        IF (ABS(YDIA).LE.DTEST) THEN
                           BVECS(KSJR+KZWOPT+I) =
     &                           BVECS(IZ+KZVAR+JROFF)
     &                           /SIGN(DTEST,YDIA)
                        ELSE
                           BVECS(KSJR+KZWOPT+I) =
     &                           BVECS(IZ+KZVAR+JROFF)/ YDIA
                        ENDIF
 1500                CONTINUE
                  ENDIF
               ENDIF
              ENDIF
            ENDIF
 2000    CONTINUE
         NTCONV = NTCONV + NOTCNV
         IF (NTCONV.EQ.0 .OR. NCSIM+NOSIM.EQ.0) GO TO 3999
C        ... NTCONV is number of non-converged response vectors
C
      IF ((IPRRSP.GT.110).AND.(NCSIM.GT.0)) THEN
         WRITE(LUPRI,*)' NCSIM CSF TRIAL VECTORS',NCSIM
         CALL OUTPUT(BVECS,1,KZCONF,1,NCSIM,KZCONF,NCSIM,1,LUPRI)
      ENDIF
      IF ((IPRRSP.GT.110).AND.(NOSIM.GT.0)) THEN
         WRITE(LUPRI,*)' NBX ORBITAL VECTORS WRITTEN OUT',NBX
         WRITE(LUPRI,*)' NOSIM ORBITAL TRIAL VECTORS    ',NOSIM
         WRITE(LUPRI,*)' KZYWOP                         ',KZYWOP
         CALL OUTPUT(BVECS(KFIXOR),1,KZYWOP,
     *             1,NBX,KZYWOP,NBX,1,LUPRI)
      ENDIF
         KCVEC = 1
         KOVEC = KCVEC + NCSIM * KZCONF
         KWRK1 = KOVEC + NOSIM * KZYWOP
         LWRK1 = LWRK  - KWRK1
         IF (LWRK1.LT.0) CALL ERRWRK('RSPNEXMCD 4',KWRK1-1,LWRK)
         KOOFF = KOVEC
         DO 2010 JR = 1,NBX
            IF (KEX(IBOFF+JR) .EQ. 0) THEN
               DO 2011 II = 1,KZYWOP
                  BVECS(KOOFF-1+II) = BVECS(KFIXOR+(JR-1)*KZYWOP-1+II)
 2011          CONTINUE
C              CALL DCOPY(KZYWOP,BVECS(KFIXOR+(JR-1)*KZYWOP),1,
C    *                           BVECS(KOOFF),1)
               KOOFF = KOOFF + KZYWOP
            ENDIF
            IF (JR.LE.NCSIM) THEN
               IBTYP(KOFFTY+KZRED+NTEST+JR) = JBCNDX
            ELSE
               IBTYP(KOFFTY+KZRED+NTEST+JR) = JBONDX
            ENDIF
 2010    CONTINUE
         IF (OPTORB.AND.(NOSIM.GT.0)) THEN
C
C GET DIAGONAL ORBITAL PART OF E(2) AND S(2)
C
            CALL DCOPY(KZWOPT,BVECS(KDIAE+KZCONF),1,BVECS(KTOT),1)
            CALL DCOPY(KZWOPT,BVECS(KDIAE+KZCONF),1,
     *                        BVECS(KTOT+KZWOPT),1)
            CALL DCOPY(KZWOPT,BVECS(KDIASO),1,BVECS(KTOT+KZYWOP),1)
            CALL DCOPY(KZWOPT,BVECS(KDIASO),1,
     *                        BVECS(KTOT+KZYWOP+KZWOPT),1)
            CALL DSCAL(KZWOPT,DM1,BVECS(KTOT+KZYWOP+KZWOPT),1)
            LWRKPA = LTOT - 2*KZYWOP
            CALL ORPPAR(NOSIM,THCORP,EOVAL,IBTYP,BVECS(KOVEC),
     *                  BVECS(KTOT), BVECS(KTOT+KZYWOP),
     *                  CMO,UDV,PVX,FC,FV,FCAC,
     *                  XINDX,BVECS(KTOT+2*KZYWOP),LWRKPA)
C
C           CALL ORPPAR(NOSIM,THCORP,EOVAL,IBTYP,A1,ORBDIE,ORBDIS,
C    *                  CMO,UDV,PVX,FC,FV,FCAC,XINDX,WRK,LWRK)
C
         ENDIF
C
C        ORTHOGONALIZE TRIAL VECTORS AND EXAMINE FOR LINEAR DEPENDENCE
C
         NFIN = NCSIM + NOSIM
         CALL RSPORTmcd(BVECS,NFIN,(KZRED+NTEST+KOFFTY),IBTYP,
     *                THRLDP,BVECS(KTOT),Xf,S2Xf,LURSP3)
C        CALL RSPORT (BVECS,NBX,NBPREV,IBTYP,THRLDP,OLDVEC,LU3)
C
C        NUMBER OF NEW TRIAL VECTORS (NFIN is updated by RSPORT)
C
      IF (IPRRSP.GT.105) THEN
         NCSIM = 0
         NOSIM = 0
         NTSIM = KZRED + KOFFTY + NTEST
         DO 3010 IX = 1,NFIN
            IF (IBTYP(NTSIM+IX).EQ.JBCNDX) THEN
               NCSIM = NCSIM + 1
            ELSE
               NOSIM = NOSIM + 1
            ENDIF
 3010    CONTINUE
      ENDIF
      IF ((IPRRSP.GT.105).AND.(NCSIM.GT.0)) THEN
         WRITE(LUPRI,*)' NCSIM CSF trial vectors after RSPORT',NCSIM
         CALL OUTPUT(BVECS,1,KZCONF,1,NCSIM,KZCONF,NCSIM,1,LUPRI)
      ENDIF
      IF ((IPRRSP.GT.105).AND.(NOSIM.GT.0)) THEN
         WRITE(LUPRI,'(/A,I5/A)')
     +      ' NOSIM orbital trial vectors after RSPORT',NOSIM,
     +      ' ( columns: Z_1 Y_1 Z_2 Y_2 ... Z_NOSIM Y_NOSIM)'
         CALL OUTPUT(BVECS(1+NCSIM*KZCONF),1,KZWOPT,
     *             1,NOSIM*2,KZWOPT,NOSIM*2,1,LUPRI)
      ENDIF
C
 3999 CONTINUE
         IF (IPRRSP.GT.5) THEN
            WRITE(LUPRI,'(/A,I5,A,I5)')
     *      ' NUMBER OF TRIAL VECTORS IN THIS LOAD',NBX,
     *      '   OFFSET FOR THIS LOAD',IBOFF
            WRITE(LUPRI,'(I5,A,I5)') NOTCNV,
     *      ' NON-CONVERGED SOLUTION VECTORS THIS LOAD, TOTAL:',NTCONV
            WRITE(LUPRI,'(I5,A,I5)') NFIN,
     *      ' LINEAR INDEPENDENT TRIAL VECTORS ADDED TO PREVIOUS',NTEST
         END IF
         NTEST  = NTEST + NFIN
         IBOFF  = IBOFF + NBX
         IF (KZYRED + 2*NTEST .GT. MAXRM) GO TO 4010
C        ... test if maximum dimension of reduced space exceeded
 4000 CONTINUE
 4010 CONTINUE
      JEXSIM = NTEST
C
C JEXSIM: THE ACTUAL NUMBER OF NEW TRIAL VECTORS IN THE
C        NEXT MICROITERATION
C
C ---
C
C Output:
C
      IF (KEXSIM.GT.MAXVEC ) WRITE (LUPRI,95) KEXSIM,MAXVEC
 95   FORMAT(/,'  ** RSPNEXMCD ** WARNING, NUMBER OF SIMULTAN TRIAL',
     * ' VECTORS',I5,/,' EXCEEDS THE PRINT MAXIMUM ',I5)
C
      IF (NTCONV.EQ.0) THEN
C
C        EIGENVECTORS CONVERGED
C
         IF (IPRRSP .GE. 0) WRITE (LUPRI,5030) KEXCNV
 5030    FORMAT(/' *** THE REQUESTED',I5,' SOLUTION VECTORS CONVERGED')
         KCONV=1
      ELSE IF (NTEST.EQ.0) THEN
C
C        LINEAR DEPENDENCE BETWEEN TRIAL VECTORS
C
         IF (LMAXIT) THEN
            WRITE(LUPRI,5020)
            KCONV=0
         ELSE
            WRITE(LUPRI,5010)
            KCONV=-1
         END IF
 5010    FORMAT(/' *** MICROITERATIONS STOPPED DUE TO LINEAR',
     *           ' DEPENDENCE BETWEEN NEW TRIAL VECTORS')
 5020    FORMAT(/' *** MICROITERATIONS STOPPED DUE TO MAX',
     *           ' ITERATIONS REACHED.')
         WRITE (LUPRI,5031) KEXCNV
 5031    FORMAT(/' *** REQUESTED',I5,' SOLUTION VECTORS NOT CONVERGED')
      ELSE IF (KZYRED + 2*NTEST .GT. MAXRM) THEN
C
C        MAXIMUM DIMENSION OF REDUCED SPACE EXCEEDED
C
         WRITE (LUPRI,'(/A/A,I5,A/)')
     *   ' *** MICROITERATIONS STOPPED BEFORE CONVERGENCE BECAUSE',
     *   '     MAXIMUM DIMENSION OF REDUCED SPACE',MAXRM,' EXCEEDED.'
         WRITE (LUPRI,5031) KEXCNV
         KCONV = -2
      ELSE
C
C        MICROITERATIONS NOT CONVERGED
C
         KCONV=0
      ENDIF
      IF (LMAXIT .OR. KCONV .LT. 0 .OR. IPRRSP .GE. 11
     &    .OR. (IPRRSP .GE. 2 .AND. KCONV .EQ. 1) ) THEN
         WRITE (LUPRI,96) THCRSP,KZYRED,(I,TEST(I),I=1,KEXCNV)
   96    FORMAT(/' Convergence of RSP solution vectors,',
     &           ' threshold =',1P,D9.2,
     *          /' ------------------------------------',
     &           '---------------------------',
     *          /' (dimension of paired reduced space:',I5,')',
     *          /,(' RSP solution vector no.',I5,'; norm of residual',
     *          1P,D11.2))
      END IF
C
C     Transfer information about the convergence back to ABACUS
C
      IF (KCONV .NE. 1) THEN
         FINRES = TEST(1)
         NOCONV = .TRUE.
      ELSE
         FINRES = TEST(1)
         NOCONV = .FALSE.
      END IF
      RETURN
C
C     END OF RSPNEXMCD.
C
      END
!  /* Deck rspredmcd */
      SUBROUTINE RSPREDMCD(ICTL,LINEQ,N,IBTYP,GD,REDGD,REDE,REDS,
     *                   EIVAL,EIVEC,BCVEC,BOVEC,UDV,EVECS,
     *                   XINDX,WRK,LWRK)
!
! Modified version of RSPRED for MCD calculations, Sonia, Nov 2009
! based on RSPRED of 30-oct-1984 by Jeppe Olsen.
! The reduced space is solved for dimension 2n-1 instead 
! of 2n, b_1=Xf is removed
! Modified 6-june-1986 to handle split orbital and csf
! trial vectors
!
! PURPOSE:
! Solve mcrsp/rpa problem in reduced subspace.
! the structure of E[2] and S[2] is used to obtain
! a reduced problem with pairing of eigensolutions.
! Allthough only KZRED vectors are known, a reduced
! problem with KZRED*2-1 vectors is used.
!
! Input:
!  ICTL, flow control
!         =1, extend the reduced L-matrix
!         =2, find the new eigenvalues and -vectors
!         =3, both
!  REDE,REDS, the old reduced matrices (dimension: KZYRED-2*N)
!  BCVEC,BOVEC, the N new vectors
!  EVECS, E times the N b-vectors
!  S times the N b-vectors are created in RSPSLI and
!  stored in EVECS after EVECS is used
!  N, number of new vectors
!
!
! Output:
!  REDE, the new, extended reduced E[2]-matrix (dimension: KZYRED)
!  REDS, the new, extended reduced S[2]-matrix (dimension: KZYRED)
!  EIVEC, eigenvectors of the new reduced mcrsp/rpa eq.
!  EVALR, real part of eigenvalues of the new reduced  mcrsp/rpa eq.
!
! Scratch:
!  WRK, real scratch array dimension kzyvar
!  bcvec ( in second part) dimension 3*kzyred**2+2*kzyred
#include <implicit.h>
#include <dummy.h>
      DIMENSION IBTYP(*),GD(*),REDGD(*),REDE(*),REDS(*),
     *          UDV(*),BOVEC(*),BCVEC(*),EVECS(KZYVAR,*)
      DIMENSION EIVEC(KZYRED,*),EIVAL(*),XINDX(*),WRK(*)
!
#include <infrsp.h>
#include <ibndxdef.h>
!
      PARAMETER (D0 = 0.0D0 , D1 = 1.0D0 , DRTEST = 1.0D-5)
      PARAMETER (MAXNVE = 200 )
!
! Used from common blocks:
!  /WRKRSP/: KZVAR,KZYRED,KZRED,KZYVAR,?
!  /INFIND/: IROW(*). ?
!  /INFPRI/: ??
!
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inftap.h>
#include <wrkrsp.h>
#include <infind.h>
#include <infpri.h>
#include <absorp.h>
      DIMENSION IVENU(MAXNVE)
!
      LOGICAL LINEQ
!
      DIMENSION ISNDX(3),DET(2)
!
      CALL QENTER('RSPREDMCD')
      IF (KZYRED.GT.MAXRM) THEN
         WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *   ' >>> ERROR IN RSPRED >>>',
     *   ' DIMENSION OF REDUCED SPACE IS  ',KZYRED,
     *   ' WHICH EXCEEDS ALLOWED DIMENSION',MAXRM
         CALL QUIT('RSPRED: TOO LARGE DIMENSION OF REDUCED SPACE')
      END IF
!SPAS:01.05.2009 bug fix for large number of  excitation energies
      IF (KZYRED.GT.LIROW) THEN
         WRITE(LUPRI,'(//A/A,A,I5,A,I5,/A)')
     *   ' >>> ERROR IN RSPRED >>>',
     *   ' DIMENSION OF INDEX VECTOR IROW IS TO SMALL; ',
     *   ' NEEDED ',KZYRED,', GIVEN :',LIROW,
     *   ' CHANGE DIMENSION IN infind.h AND RECOMPILE'
         CALL QUIT('RSPRED: TOO SMALL DIMENSION OF IROW: RECOMPILE')
      END IF
!KeinSPASmehr
!
! Section 1: extend reduced E[2]-AND S[2]-matrices
!            with N new b-vectors
!
      IF (ICTL.LT.1 .OR. ICTL.GT.3) GO TO 2000
      IF (ICTL.EQ.2) GO TO 1000
      IF (N .GT. MAXNVE) THEN
         WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *   ' >>> ERROR IN RSPRED >>>',
     *   ' Number of new trial vectors is ',N,
     *   ' which exceeds the allowed maximum, MAXNVE =',MAXNVE
         CALL QUIT('RSPRED: TOO MANY NEW TRIAL VECTORS')
      END IF
      IF (N.LT.1) GO TO 1000

      IZRED = KZRED - N
!
! Memory check
!
      LNEED = MAX(KZYWOP,KZCONF,NCREF)
!     ... CREF needed in RSPSLI (+ more)
      IF (LNEED .GT. LWRK) CALL ERRWRK('RSPRED',-LNEED,LWRK)
!
! CREATE POINTER VECTOR FOR E(2) LINEAR TRANSFORMED VECTORS
! IN CORE
!
      NCTOT = 0
      NOTOT = 0
      DO 4000 INUM = 1,N
         IF ( IBTYP(KOFFTY+IZRED+INUM) .EQ. JBONDX ) THEN
            NOTOT = NOTOT + 1
         ELSE
            NCTOT = NCTOT + 1
         END IF
 4000 CONTINUE
      ICTOT = 0
      IOTOT = 0
      DO 4010 INUM = 1,N
         IF ( IBTYP(KOFFTY+IZRED+INUM) .EQ. JBONDX ) THEN
            IOTOT = IOTOT + 1
            IVENU(INUM) = NCTOT + IOTOT
         ELSE
            ICTOT = ICTOT + 1
            IVENU(INUM) = ICTOT
         END IF
 4010 CONTINUE
      !sonia
      IF ( IPRRSP .GT. 10 ) THEN
         WRITE(LUPRI,' (/A)')
     *   ' INDEX VECTOR FOR E(2) TRANSFORMED TRIAL VECTORS'
         WRITE(LUPRI,*) ( IVENU(INUM),INUM=1,N )
         WRITE(LUPRI,' (/A)')
     *   ' IBTYP VECTOR FOR E(2) TRANSFORMED TRIAL VECTORS'
         WRITE(LUPRI,*)
     *   ( IBTYP(INUM),INUM=KOFFTY+IZRED+1,KOFFTY+IZRED+N )
      END IF
      IF (LINEQ) THEN
!
! EXTEND REDUCED GRADIENTS
!
         IZYRED= 2*IZRED
         IF (IPRRSP.GT.100) THEN
            WRITE(LUPRI,*)'  GRADIENT VECTOR (Z-part and Y-part)'
            CALL OUTPUT(GD,1,KZVAR,1,2,KZVAR,2,1,LUPRI)
            WRITE(LUPRI,*)KEXSIM,' REDUCED GRADIENT VECTORS DIM IZYRED',
     &         IZYRED
            IF (IZYRED.GT.0)
     *         CALL OUTPUT(REDGD,1,IZYRED,1,KEXSIM,IZYRED,KEXSIM,1,
     &                     LUPRI)
         ENDIF
         DO 14 I=1,KEXSIM
            IF (IZYRED.GT.0)
     *         CALL DCOPY(IZYRED,REDGD((I-1)*IZYRED+1),1,EIVEC(1,I),1)
 14      CONTINUE
         DO 12 J = 1,KEXSIM
            ISTBC = 1
            ISTBO = 1
            DO 10 I=1,N
               IF (IBTYP(IZRED+KOFFTY+I).EQ.JBCNDX) THEN
                  PZY = DDOT(KZCONF,BCVEC(ISTBC),1,GD(1),1)
                  PYZ = DDOT(KZCONF,BCVEC(ISTBC),1,GD(1+KZVAR),1)
                  ISTBC = ISTBC + KZCONF
               ELSE
                  PZY = DDOT(KZWOPT,BOVEC(ISTBO),1,GD(1+KZCONF),1)
     *                 +DDOT(KZWOPT,BOVEC(ISTBO+KZWOPT),1,
     *                       GD(1+KZVAR+KZCONF),1)
                  PYZ = DDOT(KZWOPT,BOVEC(ISTBO),1,
     *                       GD(1+KZVAR+KZCONF),1)
     *             + DDOT(KZWOPT,BOVEC(ISTBO+KZWOPT),1,GD(1+KZCONF),1)
                  ISTBO = ISTBO + KZYWOP
               ENDIF
               EIVEC(IZYRED+2*I-1,J) = PZY
               EIVEC(IZYRED+2*I,J)   = PYZ
 10         CONTINUE
 12      CONTINUE
         NTOT = KEXSIM*KZYRED
         CALL DCOPY(NTOT,EIVEC,1,REDGD,1)
         IF (IPRRSP.GE.15) THEN
            WRITE(LUPRI,'(/I3,A,I5)') KEXSIM,
     &         ' REDUCED GRADIENT VECTORS DIM KZYRED', KZYRED
            CALL OUTPUT(REDGD,1,KZYRED,1,KEXSIM,KZYRED,KEXSIM,1,LUPRI)
         ENDIF
      ENDIF
C
C  EXTEND REDUCED S-AND E-MATRICES
C
C  b-vectors times new S-and E-vectors
C
      CALL REWSPL(LURSP3)
#if !defined (VAR_SPLITFILES)
      IF (KOFFTY.EQ.1) READ (LURSP3)
#else
      JRSP3 = 0
      IF (KOFFTY.EQ.1) CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#endif
      KBOVEC = 1
      KBCVEC = 1
      DO 15 J = 1,KZRED
         IBTYPJ = IBTYP(J+KOFFTY)
         IF (IBTYPJ.EQ.JBONDX) THEN
            JZYVAR = KZYWOP
            JZVAR = KZWOPT
            JZOFF = 1+KZCONF
            JYOFF = KZVAR+JZOFF
         ELSE
            JZYVAR = KZCONF
            JZVAR = KZCONF
            JZOFF = 1
            JYOFF = KZVAR+1
         ENDIF
         IF(J.LE.IZRED) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,WRK,JRSP3,.TRUE.)
            IF (LRSP3 .NE. JZYVAR) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP3 in RSPRED')
#else
            CALL READT(LURSP3,JZYVAR,WRK)
#endif
            KMIN = 1
         ELSE
            IF (IBTYP(J+KOFFTY).EQ.JBONDX) THEN
               CALL DCOPY(JZYVAR,BOVEC(KBOVEC),1,WRK,1)
               KBOVEC = KBOVEC + JZYVAR
            ELSE
               CALL DCOPY(JZYVAR,BCVEC(KBCVEC),1,WRK,1)
               KBCVEC = KBCVEC + JZYVAR
            ENDIF
            KMIN = J-IZRED
         END IF
         DO 13 K=KMIN,N
            KEFF = IVENU(K)
            X1 = DDOT(JZVAR,WRK,         1,EVECS(JZOFF,KEFF),1)
            IF ( IBTYPJ .EQ. JBONDX ) X1 = X1
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JYOFF,KEFF),1)
            X2 = DDOT(JZVAR,WRK,         1,EVECS(JYOFF,KEFF),1)
            IF ( IBTYPJ .EQ. JBONDX ) X2 = X2
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JZOFF,KEFF),1)
            REDE(2*J  +IROW(2*(IZRED+K)))    = X1
            REDE(2*J-1+IROW(2*(IZRED+K)-1) ) = X1
            IF(J.NE.(K+IZRED)) REDE(2*J+IROW(2*(IZRED+K)-1)) = X2
            REDE(2*J-1+IROW(2*(IZRED+K)))    = X2
   13    CONTINUE
   15 CONTINUE
C
C *********************************************************
C
      NCSIM = 0
      NOSIM = 0
      DO 200 I = 1,N
         IF ( IBTYP(KOFFTY+IZRED+I) .EQ. JBONDX ) THEN
            NOSIM = NOSIM + 1
         ELSE
            NCSIM = NCSIM + 1
         END IF
 200  CONTINUE
      NTOT = (NOSIM+NCSIM) * KZYVAR
      CALL DZERO(EVECS,NTOT)
      CALL RSPSLI(NCSIM,NOSIM,BCVEC,BOVEC,
     *                  UDV,EVECS,XINDX,WRK,LWRK)
C     CALL RSPSLI(NCSIM,NOSIM,ZYCVEC,ZYOVEC,
C    *                  UDV,SVEC,XINDX,WRK,LWRK)
      CALL REWSPL(LURSP3)
#if !defined (VAR_SPLITFILES)
      IF (KOFFTY.EQ.1) READ (LURSP3)
#else
      JRSP3 = 0
      IF (KOFFTY.EQ.1) CALL READST(LURSP3,LRSP3,DUMMY,JRSP3,.FALSE.)
#endif
      KBOVEC = 1
      KBCVEC = 1
      DO 25 J = 1,KZRED
         IBTYPJ = IBTYP(J+KOFFTY)
         IF (IBTYPJ.EQ.JBONDX) THEN
            JZYVAR = KZYWOP
            JZVAR = KZWOPT
            JZOFF = 1+KZCONF
            JYOFF = KZVAR+JZOFF
         ELSE
            JZYVAR = KZCONF
            JZVAR = KZCONF
            JZOFF = 1
            JYOFF = KZVAR+1
         ENDIF
         IF(J.LE.IZRED) THEN
#if defined (VAR_SPLITFILES)
            CALL READST(LURSP3,LRSP3,WRK,JRSP3,.TRUE.)
            IF (LRSP3 .NE. JZYVAR) CALL QUIT('Inconsistent record '//
     &           'length on unit LURSP3 in RSPRED (2)')
#else
            CALL READT(LURSP3,JZYVAR,WRK)
#endif
            KMIN = 1
         ELSE
            IF (IBTYP(J+KOFFTY).EQ.JBONDX) THEN
               CALL DCOPY(JZYVAR,BOVEC(KBOVEC),1,WRK,1)
               KBOVEC = KBOVEC + JZYVAR
            ELSE
               CALL DCOPY(JZYVAR,BCVEC(KBCVEC),1,WRK,1)
               KBCVEC = KBCVEC + JZYVAR
            ENDIF
            KMIN = J-IZRED
         END IF
         DO 23 K=KMIN,N
            KEFF = IVENU(K)
            X1 = DDOT(JZVAR,WRK,         1,EVECS(JZOFF,KEFF),1)
            IF (IBTYPJ .EQ. JBONDX) X1 = X1
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JYOFF,KEFF),1)
            X2 = DDOT(JZVAR,WRK,         1,EVECS(JYOFF,KEFF),1)
            IF (IBTYPJ .EQ. JBONDX) X2 = X2
     *         + DDOT(JZVAR,WRK(1+JZVAR),1,EVECS(JZOFF,KEFF),1)
            REDS(2*J-1+IROW(2*(IZRED+K)-1)) =  X1
            REDS(2*J  +IROW(2*(IZRED+K)))   = -X1
            IF(J.NE.(IZRED+K)) REDS(2*J+IROW(2*(IZRED+K)-1)) = X2
            REDS(2*J-1+IROW(2*(IZRED+K)))   = -X2
   23    CONTINUE
   25 CONTINUE
C
C *********************************************************
C Section 2: find eigenvalues and -vectors of reduced L-matrix
C
 1000 CONTINUE
      IF (IPRRSP .GE. 0) THEN
         WRITE (LUPRI,'(/A)') ' RSPRED: REDUCED HESSIAN MATRIX:'
         CALL OUTPAK(REDE,KZYRED,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' RSPRED: REDUCED METRIC MATRIX:'
         CALL OUTPAK(REDS,KZYRED,1,LUPRI)
      END IF
      IF (ICTL .EQ. 1) GO TO 2000
C
      IF (LINEQ) THEN
C
C     Solve reduced linear response problem in subspace
C
         DO 65 IGD = 1,KEXSIM
            AFREQ = EIVAL(IGD)
            IJ=0
            IJ2=0
            !extract the relevant blocks of the LHS
            DO 70 I=1,KZYRED
               DO 75 J=1,I
                  IJ = IJ + 1
                  if ((I.ne.1).and.(j.ne.1)) then
                     IJ2=IJ2+1
                     BCVEC(IJ2) = REDE(IJ) - AFREQ*REDS(IJ)
                  end if
   75          CONTINUE
   70       CONTINUE
            length = KZYRED - 1
            !Sonia: set first element  to zero to make sure b_1=Xf does
            !not contribute when we go from red space to full space
            EIVEC(1,IGD) = 0.0d0
            CALL DCOPY(length,REDGD((IGD-1)*KZYRED+2),1,
     *                  EIVEC(2,IGD),1)
         IF (IPRRSP.GE.10) THEN
            WRITE(LUPRI,*)' AFREQ,KZYRED,length,IGD',AFREQ,KZYRED,
     &                      length,IGD
         ENDIF
         IF (IPRRSP.GT.10) THEN
            WRITE(LUPRI,'(/A)') ' 2n-1 REDUCED E(2)-W*S(2) MATRIX:'
            CALL OUTPAK(BCVEC,length,1,LUPRI)
          WRITE(LUPRI,*)'REDUCED GRADIENT VECTORS DIM KZYRED-1',length
            CALL OUTPUT(EIVEC(2,IGD),1,length,1,1,length,1,1,LUPRI)
         ENDIF
C
C           USE DSPSLI ( WHICH CALLS LINPACK ROUTINES ) TO SOLVE
C           SYSTEM OF LINEAR EQUATIONS
C
            NSIM = 1
            KK1  = 1 + length*(length+1)/2
            CALL DSPSLI(length,NSIM,BCVEC,EIVEC(2,IGD),BCVEC(KK1),INFO,
     *                  DET,ISNDX)
C
C           CALL DSPSLI(N,NSIM,AP,B,KPIVOT,INFO,DET,INERT)
C
            IF (INFO .NE. 0) THEN
               WRITE (LUPRI,'(/A,I5)')
     *              ' *** ERROR IN RSPRED.DSPSLI, INFO =',INFO
               CALL QUIT('*** ERROR IN RSPRED.DSPSLI')
            END IF
C
C
            INEG     = ISNDX(2)
            ISNDX(2) = ISNDX(3)
            ISNDX(3) = INEG
            IF (ISNDX(3).GT.0 .AND. .NOT.ABSORP) 
     *           WRITE(LUPRI,'(/A)')
     *         ' *** INFO, negative eigenvalues in reduced matrix.'
            IF (IPRRSP.GE.10.OR.(.NOT.ABSORP .AND. ISNDX(3).GT.0)) THEN
               !POL=DDOT(KZYRED,EIVEC(1,IGD),1,REDGD((IGD-1)*KZYRED+1),1)
               !Sonia: changed for 2n-1 case
               POL=DDOT(length,EIVEC(2,IGD),1,REDGD((IGD-1)*KZYRED+2),1)
               WRITE(LUPRI,1800) AFREQ,KZRED,POL,ISNDX,DET
            END IF
 1800       FORMAT(/' GP * SOLUTION vector at frequency',F12.5,' a.u.',
     *             /' after',I5,' linear transformations is ',F15.8,
     *             /' INERTIA (POS,ZER,NEG) of reduced matrix is',3I5,
     *             /' Determinant of reduced matrix is',F6.2,
     *              ' * 10 **',F5.1)
C
            IF (IPRRSP .GE. 10) THEN
               WRITE (LUPRI,'(/A)')
     *            ' REDUCED PROPERTY VECTOR AND REDUCED SOLUTION:'
               WRITE (LUPRI,'(1P,/I5,2D15.6/I5,2D15.6)')
     *           (K,REDGD((IGD-1)*KZYRED+K+1), EIVEC(K,IGD), K=1,length)
            END IF
 65      CONTINUE
      ELSE
         !SONIA---------------------------------------------------
         call quit('RSPREDMCD: NOT TO BE USED FOR EIGEN EQUATION')
C        ---------------------------------------------------------
C Solve reduced EIGENVALUE general problem in subspace
C Expand reduced matrices to square matrices
C
         KZYRDQ = KZYRED*KZYRED
         IJ=0
         DO 77 I=1,KZYRED
            DO 78 J=1,I
               IJ=IJ+1
               IJE=(I-1)*KZYRED+J
               JIE=(J-1)*KZYRED+I
               BCVEC(IJE)=REDE(IJ)
               BCVEC(JIE)=REDE(IJ)
               BCVEC(IJE+KZYRDQ)=REDS(IJ)
               BCVEC(JIE+KZYRDQ)=REDS(IJ)
   78       CONTINUE
   77    CONTINUE
C
         KK1 = 1   + 2*KZYRDQ
         KK2 = KK1 + KZYRED
         KK3 = KK2 + KZYRED
         KK4 = KK3 + KZYRDQ
C
C        Use EISPACK routine for real general matrices in
C        generalized eigenvalue problem
C
         CALL RGG(KZYRED,KZYRED,BCVEC(1),BCVEC(1+KZYRDQ),EIVAL,
     *            BCVEC(KK1),BCVEC(KK2),1,EIVEC,IERR)
C        CALL RGG(NM,N,A,B,ALFR,ALFI,BETA,MATZ,Z,IERR)
C
         IF ( IERR.NE.0 ) THEN
            WRITE(LUPRI,'(/A)')
     *      ' EIGENVALUE PROBLEM NOT CONVERGED '
            CALL QUIT(' RSPRED: EIGENVALUE EQUATION NOT CONVERGED ')
         END IF
C
      IF (IPRRSP .GE. 25) THEN
         WRITE (LUPRI,'(/A)') ' REDUCED EIGENVECTORS BEFORE ORDERING:'
         CALL OUTPUT(EIVEC,1,KZYRED,1,KZYRED,KZYRED,KZYRED,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES:'
         CALL OUTPUT(EIVAL,1,KZYRED,1,1,KZYRED,1,1,LUPRI)
      END IF
         IJ = 0
         DO 72 I=1,KZYRED
            DO 72 J=1,I
               IJ=IJ+1
               IJE=(I-1)*KZYRED+J
               JIE=(J-1)*KZYRED+I
               BCVEC(IJE)=REDS(IJ)
               BCVEC(JIE)=REDS(IJ)
   72    CONTINUE
C
         CALL RSPORD(BCVEC(1),EIVEC,EIVAL,BCVEC(KK1),BCVEC(KK2),
     *            BCVEC(KK3),BCVEC(1+KZYRDQ),ISNDX)
C        CALL RSPORD(SRED,EIVEC,ALFAR,ALFI,BETA,BVECS1,BVECS2,ISNDX)
C
         IF (IPRRSP .GE. 10) THEN
            IF (IPRRSP .GE. 25) THEN
               WRITE (LUPRI,'(/A)')
     *               ' REDUCED EIGENVECTORS AFTER ORDERING:'
               CALL OUTPUT(EIVEC,1,KZYRED,1,KZYRED,KZYRED,KZYRED,1,
     &                     LUPRI)
            END IF
            WRITE(LUPRI,'(//2A,3I5//A/A)')
     *         ' NUMBER OF EIGENVALUES WITH POSITIVE,',
     *         ' ZERO, AND NEGATIVE METRIC:',ISNDX,
     *         ' THE EIGENVALUES WITH POSITIVE METRIC:',
     *         ' NUMBER         EIGENVALUE'
            IPOS = ISNDX(1)
            DO 85 I=1,IPOS
               IF (ABS(EIVAL(I)) .LT. DRTEST) THEN
                 WRITE(LUPRI,'(I10,1P,D15.8,A/)') I,EIVAL(I),
     *            ' *** RSPRED  WARNING **** ZERO EIGENVALUE.'
               ELSE
                  WRITE(LUPRI,'(I10,1P,D15.8)') I,EIVAL(I)
               END IF
 85         CONTINUE
         END IF
      ENDIF
C
C *** End of subroutine RSPREDMCD
C
 2000 CONTINUE
      CALL QEXIT('RSPREDMCD')
      RETURN
      END
