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
C  /* Deck broydn */
      SUBROUTINE BROYDN(LUIT7,KIT,NTVAR,VECK,VECG,VECT,VECQ,
     *                  WRK1,WRK2,BROYIT,BCKSTP)
C
C LAST REVISION 23-OCT-1985
C
C PURPOSE:
C  Carry out Broyden inverse update iteration no. k
C
C    H(k) = H(0) - SUM(i=1,k) t(i) * q(i)T
C
C  where
C
C    t(k) = H(k-1) * l(k) - k(k)
C         if k(k) is Broyden step then
C         = H(k-1) * g(k)
C
C    q(k) = H(k-1)T * k(k) / (NORM)
C
C    l(k) = g(k) - g(k-1)
C
C    k(k) = x(k) - x(k-1)
C
C    NORM = k(k)T * H(k) * l(k)
C         if k(k) is Broyden step then
C         = k(k)T * ( t(k) + k(k) )
C
C   and where
C
C   H(0) IS THE INITIAL APPROXIMATION TO THE INVERSE JACOBIAN
C   H(k) IS THE APPROXIMATION TO THE INVERSE JACOBIAN IN
C        ITERATION KIT
C   x(k) IS MC POINT IN ITERATION k
C   g(k) = VECG IS THE GRADIENT AT ITERATION k
C
C   -----
C
C   delta x(k) IS THE BROYDEN STEP TO x(k+1):
C
C    delta x(k) = k(k+1) = -H(k) * g(k)
C         if k(k) is Broyden step then
C         = (- 1 + q(k)T * g(k) ) t(k)
C
C   -----
C
C   If (BCKSTP) the previous step, k(k), has been rejected by the
C   controlling program.  We update the Jacobian information using the
C   rejected step, reset the expansion point (i.e. x(k) := x(k-1)),
C   and evaluate a revised step, k(k+1), from x(k) = x(k-1).
C
C INPUT:
C   KIT   = k, THE ITERATION NUMBER IN BROYDEN UPDATE
C   NTVAR = DIMENSIONS OF VECTOR SPACE
C   VECK  = k(k), THE STEP INCREMENT TO REACH ITERATION POINT
C   VECG  = g(k), THE GRADIENT AT ITERATION POINT KIT
C   VECT  = H(0)*VECG
C   VECQ  = H(0)*VECK
C   LUIT7 HAS H(0) AS FIRST RECORD,
C         THE BROYDEN STEP -H(k-1)*g(k-1) AS SECOND RECORD, AND
C         THE GRADIENT g(k-1) AS THIRD RECORD
C   t(1),q(1),t(2),q(2),..,t(KIT-1),q(KIT-1) ARE STORED AS
C         (NREDH +I, I=1,KIT-1) RECORDS ON LUIT3 AND LUIT5 RESPECTIVELY.
C         (LUIT3 is positioned at t(1) and LUIT5 at q(1) when BROYDN
C          called.)
C
C OUTPUT:
C   k(k+1) , the Broyden step, in VECT
C   t(k)   appended to LUIT3 and q(k) to LUIT5
C   k(k+1) and g(k) on LUIT7
C   VECK   unchanged
C   VECG   if not BCKSTP then unchanged (i.e. g(k)),
C          else reset ot g(k-1).
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION VECK(*),VECG(*),VECT(*),VECQ(*),WRK1(*),WRK2(*)
      LOGICAL BROYIT,BCKSTP
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0 )
C
#include <infpri.h>
#include <inftap.h>
C
      IF (IPRSTAT .GE. 5)
     *   WRITE(LUSTAT,'(/A,I4)')' ******** BROYDN BEGIN, KIT =',KIT
      IF (IPRSTAT .GT. 16) THEN
         WRITE(LUSTAT,'(/A,I5)')' VECK, NTVAR=',NTVAR
         WRITE(LUSTAT,'(1X,5F15.8)')(VECK(I),I=1,NTVAR)
         WRITE(LUSTAT,'(/A)')' VECG'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECG(I),I=1,NTVAR)
         WRITE(LUSTAT,'(/A)')' G(0) g(k)'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECT(I),I=1,NTVAR)
         WRITE(LUSTAT,'(/A)')' G(0)T k(k)'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECQ(I),I=1,NTVAR)
      END IF
C
C     CONSTRUCT VECT = H(k-1) * g(k) AND VECQ = H(k-1)T * k(k)
C
      DO 200 I = 1, (KIT-1)
         CALL READT(LUIT3,NTVAR,WRK1)
         CALL READT(LUIT5,NTVAR,WRK2)
         XT = -DDOT(NTVAR,WRK2,1,VECG,1)
         XQ = -DDOT(NTVAR,WRK1,1,VECK,1)
         CALL DAXPY(NTVAR,XT,WRK1,1,VECT,1)
         CALL DAXPY(NTVAR,XQ,WRK2,1,VECQ,1)
 200  CONTINUE
      IF (IPRSTAT .GT. 15) THEN
         WRITE(LUSTAT,'(/A)')' H(k-1) g(k)'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECT(I),I=1,NTVAR)
         WRITE(LUSTAT,'(/A)')' H(k-1)T k(k)'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECQ(I),I=1,NTVAR)
      END IF
C
C     Finish t(k) in VECT and q(k) in VECQ.
C
      IF (BROYIT) THEN
         XQNOR =  DDOT(NTVAR,VECK,1,VECT,1) + DDOT(NTVAR,VECK,1,VECK,1)
         IF (IPRSTAT .GE. 5) THEN
            WRITE(LUSTAT,'(/A)')' Last step was a Broyden step'
            WRITE(LUSTAT,*)'Variable metric :',XQNOR
         END IF
      ELSE
         REWIND(LUIT7)
         READ  (LUIT7)
         CALL DCOPY(NTVAR,VECT,1,WRK1,1)
         CALL READT(LUIT7,NTVAR,WRK2)
         CALL DAXPY(NTVAR,D1,WRK2,1,VECT,1)
         XQNOR = DDOT(NTVAR,VECT,1,VECK,1)
         IF (IPRSTAT .GE. 5) THEN
            WRITE(LUSTAT,'(/A)')' Last step was NOT the Broyden step'
            WRITE(LUSTAT,*)'Variable metric :',XQNOR
         END IF
         CALL DAXPY(NTVAR,DM1,VECK,1,VECT,1)
      END IF
      XQNOR = D1/XQNOR
      CALL DSCAL(NTVAR,XQNOR,VECQ,1)
C
C     Write t(k) and q(k) vectors on LUIT3 and LUIT5, resp.
C
      CALL WRITT(LUIT3,NTVAR,VECT)
      CALL WRITT(LUIT5,NTVAR,VECQ)
C
C     CONSTRUCT BROYDEN STEP VECTOR IN VECT.
C     If (BCKSTP) reset VECG to g(k-1).
C
      IF (BCKSTP) THEN
         IF (BROYIT) THEN
            REWIND (LUIT7)
            READ   (LUIT7)
            CALL READT(LUIT7,NTVAR,WRK2)
         END IF
         IF (IPRSTAT .GE. 5) THEN
            XQNOR = DDOT(NTVAR,VECQ,1,VECG,1)
            WRITE(LUSTAT,'(/A)')'*** THIS IS A BACK STEP ***'
            WRITE(LUSTAT,*)' < q(k) / g(k) >   =',XQNOR
         END IF
         CALL READT(LUIT7,NTVAR,VECG)
      END IF
      XQNOR = DDOT(NTVAR,VECQ,1,VECG,1)
      IF (BCKSTP) THEN
         CALL DSCAL(NTVAR,XQNOR,VECT,1)
         CALL DAXPY(NTVAR,D1,WRK2,1,VECT,1)
      ELSE IF (BROYIT) THEN
         CALL DSCAL(NTVAR,(XQNOR-D1),VECT,1)
      ELSE
         CALL DSCAL(NTVAR,XQNOR,VECT,1)
         CALL DAXPY(NTVAR,DM1,WRK1,1,VECT,1)
      END IF
      IF (IPRSTAT .GE. 5) THEN
         IF (BCKSTP) THEN
            WRITE(LUSTAT,*)' < q(k) / g(k-1) > =',XQNOR
         ELSE
            WRITE(LUSTAT,*)' < q(k) / g(k) > =',XQNOR
         END IF
      END IF
C
C     WRITE BROYDEN STEP VECTOR AND THE GRADIENT ON LUIT7
C
      REWIND(LUIT7)
      READ  (LUIT7)
      CALL WRITT(LUIT7,NTVAR,VECT)
      CALL WRITT(LUIT7,NTVAR,VECG)
C
      IF (IPRSTAT .GT. 11) THEN
         WRITE(LUSTAT,'(/A)')' Broyden step returned'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECT(I),I=1,NTVAR)
         WRITE(LUSTAT,'(/A)')
     *      ' q vector returned (has been divided by var. metric)'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECQ(I),I=1,NTVAR)
      END IF
C
      RETURN
C
C     End of BROYDN.
C
      END
C  /* Deck gctopc */
      SUBROUTINE GCTOPC(NTVAR,CREF,C,GRD,CSCAL,GSCAL)
C
C PROJECT REFERENCE VECTOR COMPONENTS OUT OF C-VECTOR AND GRADIENT
C AND SCALE TO X REPRESENTATION'
C
C CSCAL = norm of C after scaling to unit CREF component.
C
C GSCAL = gradient scale factor to give correct first order
C         energy prediction in PC representation.
C
#include <implicit.h>
      PARAMETER( D1=1.0D0 , DM1=-1.0D0 )
      DIMENSION CREF(*),C(*),GRD(*)
C
      FAC   = DDOT(NTVAR,CREF,1,C,1)
      CSCAL = D1/FAC
      CALL DSCAL(NTVAR,CSCAL,C,1)
      CALL DAXPY(NTVAR,DM1,CREF,1,C,1)
C
      CALL DSCAL(NTVAR,FAC,GRD,1)
      GSCAL = D1
C
      FAC = -DDOT(NTVAR,CREF,1,GRD,1)
      CALL DAXPY(NTVAR,FAC,CREF,1,GRD,1)
      RETURN
      END
C  /* Deck pctoc */
      SUBROUTINE PCTOC(NTVAR,CREF,C)
C
C TRANSFORM FROM X REPRESENTATION TO C REPRESENTATION
C
#include <implicit.h>
      DIMENSION CREF(*),C(*)
      PARAMETER ( D1=1.0D0 )
C
      CALL DAXPY(NTVAR,D1,CREF,1,C,1)
      XNOR = D1/DNORM2(NTVAR,C,1)
      CALL DSCAL(NTVAR,XNOR,C,1)
      RETURN
      END
C  /* Deck sirupd */
      SUBROUTINE SIRUPD(ICONV,MAXIT,REF1,GRD,CMO1,INDXCI,WRK,LFREE)
C
C LAST REVISION 23-OCT-1985
C
C PURPOSE:
C  Do a Broyden inverse update optimization.
C
#include <implicit.h>
      DIMENSION REF1(*),GRD(*),CMO1(*),INDXCI(*),WRK(*)
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0, DP5=0.5D0 )
      PARAMETER ( XRNEG=1.D-3, XRED=0.8D0, XRBCK=0.5D0 )
      PARAMETER ( XLIM=1.D-7, DP1LIM=1.D-10, XRSTP=10.0D0 )
C     PARAMETER ( XCFROZ=12.D0 )
#include <dummy.h>
C
C used from common blocks:
C    INFINP : ISTATE,POTNUC,FLAG(*),?
C    INFOPT : EMCOLD,EMCSCF,GNORM(3)
C    INFDIM : MAXRL
C    INFPRI : P4FLAG(*),P6FLAG(*)
C
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <infinp.h>
#include <inforb.h>
#include <infopt.h>
#include <infdim.h>
#include <infpri.h>
#include <inftap.h>
#include <inftra.h>
C
      LOGICAL BCKSTP, BROYIT, DOTRA, HESFIX, LSAVE
C
      CALL QENTER('SIRUPD')
      TIMUPD = SECOND()
      WRITE (LUW4,'(///A//)')
     *   '   <<<< MC optimization control switched to update >>>>'
      IF (LUPRI .NE. LUW4) WRITE (LUPRI,'(A//)')
     *   '1  <<<< MC optimization control switched to update >>>>'
C
C
      HESFIX = FLAG(56)
C
C     ALLOCATE WORK SPACE
C
      KREDI =1
      KIBNDX=KREDI+MAXRL*MAXRL
      KREF0 =KIBNDX+MAXRL
      KCMO0 =KREF0+NVAR
      KDV   =KCMO0+NCMOT
      KPV   =KDV  +NNASHX
      KCOR  =KPV  +NNASHX*NNASHX
      LCOR  =LFREE -KCOR
C
C     ALLOCATE WORK SPACE FOR UPDST AND BROYDN
C
      KVECT =KCOR
      KVECQ =KVECT+NVAR
      KWRK1 =KVECQ+NVAR
      KWRK2 =KWRK1+NVAR
      KTOT  =KWRK2+MAX(NVAR,4*MAXRL,MAXRL*(MAXRL+1)/2)
      IF(KTOT.GT.LFREE) CALL ERRWRK('SIRUPD',KTOT,LFREE)
C
      KCREF0=KCOR
      KFC   =KCREF0 +NVAR
      KFV   =KFC    +NNORBT
      KFQ   =KFV    +NNORBT
      KXKAP =KFQ    +NASHT*NORBT
      KCORX =KXKAP  +NWOPT
      LCORX =LFREE  -KCORX
C
C     SET UP INITIAL CONDITIONS IN UPDST
C
      LUIT7 = -1
      CALL GPOPEN(LUIT7,'BROYDEN','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      TIM = SECOND()
      CALL UPDST(LUIT7,WRK(KREDI),WRK(KIBNDX),WRK(KREF0),WRK(KCMO0),
     $           REF1,GRD,WRK(KVECT),WRK(KVECQ),
     $           WRK(KWRK1),WRK(KWRK2),CSCAL,GSCAL)
C     CALL UPDST(LUIT7,REDINV,IBNDX,REF0,CMO,VECK,VECG,
C    *           VECT,VECQ,G0,VECX,CSCAL,GSCAL)
C
      TIM = SECOND() - TIM
      WRITE (LUPRI,'(/A,T40,F12.2)') ' >>>> time for UPDST: ',TIM
C
C     INITIALIZE VARIABLES USED IN BROYDEN ITERATIONS;
C     NO INTEGRAL TRANSFORMATION IF RHF.
C
      KIT    = 1
      IGRTYP = 0
      DP1EC  = D0
      DP1EO  = D0
      XCUPD  = D0
      XOUPD  = D0
      BROYIT = .FALSE.
      BCKSTP = .FALSE.
      DOTRA = .NOT. FLAG(34)
      IF (ISTATE .EQ. 1) THEN
         D1RLOW = 0.D0
         D1RHGH = 1.D5
         GRISE  = 5.D0
      ELSE
         D1RLOW = 0.2D0
         D1RHGH = 1.2D0
         GRISE  = 1.5D0
      END IF
C
C DETERMINE INITIAL MAXIMUM ALLOWED STEP IN UPDATE
C (= norm of step squared because iteration 1 was a second order step
C    provided gradient also has gone down).
C
      XSTP = DNORM2(NVAR,REF1,1) ** 2
      XSTP = MAX(XSTP,GNORM(3))
      XSTP = XRSTP * XSTP
      XMAX = XSTP
      GNOLD= D0
C     ... to avoid compiler messages
C
C START NEXT UPDATE ITERATION
C
 500  KIT = KIT + 1
      TIMIT = SECOND()
      WRITE(LUPRI,'(/A,I5)')' *** MCSCF UPDATE ITERATION NUMBER',KIT
C
C     CALCULATE FIRST ORDER PREDICTED ENERGY CHANGE, NO SCALING
C     TEST IF STEP IS SMALLER THAN XSTP
C
      IF (NCONF .GT. 1) THEN
         DP1EC = GSCAL * DDOT(NCONF,GRD,1,WRK(KVECT),1)
         XCUPD = DNORM2(NCONF,WRK(KVECT),1)
      END IF
      IF (NWOPT .GT. 0) THEN
         DP1EO = DDOT(NWOPT,GRD(1+NCONF),1,WRK(KVECT+NCONF),1)
         XOUPD = DNORM2(NWOPT,WRK(KVECT+NCONF),1)
      END IF
      DP1EMC = DP1EC + DP1EO
      XUPD  = SQRT(XCUPD**2 + XOUPD**2)
      WRITE (LUPRI,'(2(/A,T36,2F18.12))')
     *   ' UPDATE STEP LENGTH, MAX STEP',XUPD,XMAX,
     *   ' - CSF AND ORBITAL STEP LENGTHS',XCUPD,XOUPD
C
C     UPDATE MAXIMUM STEP
C
      IF (XUPD .LT. XMAX) THEN
         XSTP = MAX(XRED*XSTP,XUPD)
         IF (GNORM(3) .GT. D0) XSTP = MIN(XSTP,GNORM(3))
      END IF
C
      XSCAL = D1
      IF ( ISTATE.EQ.1 .AND. GNORM(3).GT.XLIM .AND. DP1EMC.GE.D0 ) THEN
         XSCAL = - MAX(XRNEG*XUPD,MIN(XUPD,THRGRD)) / XUPD
         XUPD  = - XSCAL * XUPD
         WRITE (LUPRI,'(/A,F18.12/A,T36,F18.12)')
     *      ' - PREDICTED ENERGY CHANGE POSITIVE',DP1EMC,
     *      '   STEP REVERSED AND REDUCED TO',XUPD
      END IF
      IF (XUPD .GT. XMAX) THEN
         WRITE (LUPRI,'(A,T36,F18.12)') ' - STEP REDUCED TO',XMAX
         XSCAL  = XMAX/XUPD
         XUPD   = XMAX
      END IF
      IF (XSCAL .NE. D1) THEN
         DP1EMC = XSCAL * DP1EMC
         CALL DSCAL(NVAR,XSCAL,WRK(KVECT),1)
         BROYIT = .FALSE.
      END IF
      IF (.NOT. BCKSTP) XMAX = XSTP
C
C     FIND NEW ORBITAL AND CSF PARAMETERS
C
      IF (IPRSTAT.GT.14) THEN
         WRITE(LUSTAT,'(/A)')' STEP VECTOR IN VECT SIRUPD'
         WRITE(LUSTAT,'(1X,5F15.8)')(WRK(KVECT-1+I),I=1,NVAR)
         WRITE(LUSTAT,'(/A)')' REF0, THE CURRENT MC POINT (BEFORE STEP)'
         WRITE(LUSTAT,'(1X,5F15.8)')(WRK(KREF0-1+I),I=1,NVAR)
      END IF
      DO 400 I=1,NVAR
         REF1(I) = WRK(KVECT-1+I) + WRK(KREF0-1+I)
 400  CONTINUE
C
C     Read in reference vector
C
      REWIND(LUIT3)
      CALL READT(LUIT3,NCONF,WRK(KCREF0))
      CALL PCTOC(NCONF,WRK(KCREF0),REF1)
      CALL DCOPY(NCONF,REF1,1,WRK(KCORX),1)
      CALL DCOPY(NWOPT,REF1(1+NCONF),1,WRK(KXKAP),1)
C
C     Save current expansion point
C     890705: disable SIRCNO call, transformation to natural orbitals,
C     which would make update algorithm meaningless.  Thus INDXCI not
C     needed in SIRSAV.
C
C     CMO1 is updated to new MO coefficients in SIRSAV, based on CMO0.
C     We must therefore restore CMO0 in CMO1 before calling SIRSAV.
C
      CALL DCOPY(NCMOT,WRK(KCMO0),1,CMO1,1)
C
      LSAVE = FLAG(15)
      FLAG(15) = .FALSE.
      CALL SIRSAV('UPDSAVE',CMO1,DUMMY,DUMMY,DUMMY,WRK(KXKAP),
     1            DUMMY,WRK(KCORX),LCORX)
      FLAG(15) = LSAVE
C     CALL SIRSAV(KEYWRD,CMO,IBNDX,REDL,EVEC,XKAP,INDXCI,WRK,LFREE)
C
      IF (DOTRA) CALL TRACTL(ITRLVL,CMO1,WRK(KCORX),LCORX)
C                CALL TRACTL(ITRLVL,CMO,WRK,LFREE)
C
      IF (.NOT. BCKSTP) GNOLD = GNORM(3)
      CALL UPDGRD(IGRTYP,REF1,CMO1,WRK(KDV),WRK(KPV),GRD,EMY,EACTIV,
     *            WRK(KFC),WRK(KFV),WRK(KFQ),
     *            DUMMY,INDXCI,WRK(KCORX),LCORX)
C
      EMCOLD = EMCSCF
      EMCSCF = POTNUC + EMY + EACTIV
      DEMC   = EMCSCF - EMCOLD
      D1RAT  = DEMC / DP1EMC
      WRITE (LUW4,1010) EMCSCF,KIT
      IF (LUPRI.NE.LUW4) WRITE (LUPRI,1010) EMCSCF,KIT
      IF (P4FLAG(2)) WRITE (LUW4,1015) POTNUC,EMY,EACTIV
      IF (P6FLAG(1)) THEN
         IF (LUPRI.NE.LUW4 .OR. .NOT.P4FLAG(2))
     *      WRITE (LUPRI,1015) POTNUC,EMY,EACTIV
      END IF
      WRITE (LUW4,1020) GNORM(3),GNORM(1),GNORM(2)
      IF (LUPRI.NE.LUW4) WRITE (LUPRI,1020) GNORM(3),GNORM(1),GNORM(2)
 1010 FORMAT(//' Total MCSCF energy       :',F25.15,T60,
     *         '(UPD    ',I3,')')
 1015 FORMAT(//' - Nuclear repulsion      :',F25.15,
     *        /' - Inactive energy        :',F25.15,
     *        /' - Active energy          :',F25.15)
 1020 FORMAT(//' Norm of total gradient   :',F22.12,
     *        /' -    of CI gradient      :',F22.12,
     *        /' -    of orbital gradient :',F22.12)
C
      WRITE (LUW4,'(/A,F20.15/A,F10.5/A,F10.5)')
     *' First order predicted energy change                   :',DP1EMC,
     *' Ratio between energy change and first order prediction:',D1RAT,
     *' Value predicted if Jacobian was energy Hessian        :',
     *(D1-DP5*XSCAL),
     *' Norm of wave function in PC representation            :',CSCAL
      IF (XUPD .LT. XLIM .OR. ABS(DP1EMC) .LT. DP1LIM) THEN
         D1RAT = D1
         WRITE (LUW4,'(/A/A)') ' Ratio reset to one because energy'//
     *      ' difference unreliable',' this close to convergence.'
      END IF
C
C CHECK IF BACKSTEP
C
      IF (D1RAT .LE. D1RLOW .OR. D1RAT .GE. D1RHGH .OR.
     *    GNORM(3) .GT. GRISE*GNOLD) THEN
         BCKSTP = .TRUE.
         WRITE (LUW4,'(//A/)') ' *** BACK STEP -- '//
     *      'THIS UPDATE STEP IS REJECTED ***'
         XMAX = XRBCK * XMAX
         EMCSCF = EMCOLD
      ELSE
         BCKSTP = .FALSE.
         XMAX = XSTP
      END IF
C
C PRINT GRADIENT, IF REQUESTED
C
      IF (P4FLAG(15) .AND. NCONF.GT.1) THEN
CHJ: in future version maybe CALL PRWF here (and above for CREF)
         PTEST = MAX(1.D-10,GNORM(1)*THRCGR)
         WRITE (LUW4,3010) KIT, PTEST
         CALL PRVEC(NCONF,GRD,1,-PTEST,200,LUW4)
         CALL PRMGN(NCONF,GRD,1,12,LUW4)
#if defined (VAR_OLDCODE) || defined (VAR_SECSEC)
         DO 301 I = 1,NCONF
            GCI = GRD(I)
            IF (ABS(GCI) .GE. PTEST) WRITE (LUW4,3011) I,GCI
  301    CONTINUE
 3011 FORMAT(I16,F20.10)
#endif
      END IF
      IF (P4FLAG(14)) THEN
         WRITE (LUW4,3020) KIT
         IF (MPRI4.GT.100) THEN
            PRFAC = 0.0D0
         ELSE IF (MPRI4.GT.10) THEN
            PRFAC = 0.1D0
         ELSE
            PRFAC = 0.2D0
         END IF
         CALL PRKAP (NWOPT,GRD(NCONF+1),PRFAC,LUW4)
      END IF
 3010 FORMAT(/' Configuration gradient update it.no.',I3,
     *       /' ---------------------------------------'
     *      //' Cutoff for print:',1P,D10.2,
     *      //' Configuration no.           value'
     *       /' -----------------           -----')
 3020 FORMAT(/' Orbital gradient update it.no.',I3,
     *       /' ---------------------------------')
C
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
C
C TEST FOR CONVERGENCE
C   GNORM(1) = GCI  NORM
C   GNORM(2) = GORB NORM
C   GNORM(3) = GRAD NORM
C
      IF(GNORM(3) .LT. THRGRD) GO TO 8000
      IF(KIT .GE. MAXIT) GO TO 8010
C
C TRANSFORM CSF PART OF GRADIENT AND CSF COEFFICIENTS
C TO PC REPRESENTATION
C
      CALL GCTOPC(NCONF,WRK(KCREF0),REF1,GRD,CSCAL,GSCAL)
C
C PUT ITERATED PARAMETER SET IN REF0 AND STEP INCREMENT IN REF1
C
      IF (BCKSTP) THEN
         CALL DAXPY(NVAR,DM1,WRK(KREF0),1,REF1,1)
      ELSE
         CALL DSWAP(NVAR,WRK(KREF0),1,REF1,1)
         CALL DAXPY(NVAR,DM1,WRK(KREF0),1,REF1,1)
         CALL DSCAL(NVAR,DM1,REF1,1)
      END IF
C
C *** DETERMINE STEP
C
C  A) H0 (inv) * g
C
      CALL UPDH0(LUIT7,WRK(KIBNDX),WRK(KREDI),REF1,GRD,
     *           WRK(KVECT),WRK(KVECQ),WRK(KWRK1),WRK(KWRK2))
C     CALL UPDH0(LUIT7,IBNDX,REDINV,VECK,VECG,VECT,VECQ,WRKX,WRK)
C
C
C B1) FIXED HESSIAN ( i.e.  - H0 (inv) * g ) or
C B2) BROYDEN UPDATE
C
      IF ( HESFIX ) THEN
         CALL DSCAL(NVAR,DM1,WRK(KVECT),1)
      ELSE
         CALL BROYDN(LUIT7,KIT,NVAR,REF1,GRD,WRK(KVECT),WRK(KVECQ),
     $               WRK(KWRK1),WRK(KWRK2),BROYIT,BCKSTP)
C        CALL BROYDN(LUIT7,KIT,NTVAR,VECK,VECG,VECT,VECQ,
C    *               WRK1,WRK2,BROYIT,BCKSTP)
      END IF
      CALL GPCLOSE(LUIT7,'KEEP')
C
C
C IF (XCFROZ*GCINORM .LT. GORBNORM) THEN FREEZE CI
C
CHJ-S-851025; frozen CI made convergence worse
C     IF (XCFROZ*GNORM(1) .LT. GNORM(2) .AND. .NOT.BCKSTP
C    *    .AND. IGRTYP .EQ. 0) THEN
C        WRITE (LUW4,'(/A)') ' *** CI frozen in next update step'
C        CALL DZERO(WRK(KVECT),NCONF)
C        IGRTYP = 1
C        BROYIT = .FALSE.
C     ELSE
         IGRTYP = 0
         BROYIT = .TRUE.
C     END IF
CHJ-E
      TIMIT = SECOND() - TIMIT
      WRITE (LUPRI,'(/A,T40,F12.2)')
     *   ' >>>> time for this update iteration:',TIMIT
      GO TO 500
C
C
 8000 CONTINUE
      WRITE(LUW4,7000) KIT,GNORM(3)
      IF (LUPRI .NE. LUW4) WRITE(LUPRI,7000) KIT,GNORM(3)
 7000 FORMAT(/' MCSCF has converged in',I5,' update iterations ',
     *       /' norm of gradient is = ',1P,D15.6)
      ICONV = 1
      GO TO 9000
C
C
 8010 CONTINUE
      NWARN = NWARN + 1
      WRITE(LUW4,7010)MAXIT,GNORM(3)
      IF (LUPRI .NE. LUW4) WRITE(LUPRI,7010)MAXIT,GNORM(3)
 7010 FORMAT(/' **WARNING**  MCSCF did not converge in',I5,
     *        ' update iterations',
     *       /' norm of gradient is = ',1P,D15.6)
      ICONV = 0
C
C
C
 9000 CONTINUE
      TIMUPD = SECOND() - TIMUPD
      WRITE (LUPRI,'(/A,T40,F12.2)')
     *    ' >>>> total time for update: ',TIMUPD
      CALL QEXIT('SIRUPD')
      RETURN
C
C     End of SIRUPD.
C
      END
C  /* Deck updst */
      SUBROUTINE UPDST(LUIT7,REDINV,IBNDX,REF0,CMO,VECK,VECG,
     *                 VECT,VECQ,G0,VECX,CSCAL,GSCAL)
C
C THIS ROUTINE ASSUME THAT A LOCAL SECOND ORDER ITERATION IS CARRIED OUT
C
C ITERATION POINT ZERO IS THE EXPANSION POINT
C ITERATION POINT ONE  IS THE POINT THAT RESULT FROM CARRYING OUT THE
C                      SECOND ORDER ITERATION
C
C AS APPROXIMATE INVERSE HESSIAN IS USED THE INVERSE REDUCED HESSIAN
C AT IT POINT ZERO TO WHICH IS ADDED THE INVERSE DIAGONAL HESSIAN.
C THE INVERSE DIAGONAL HESSIAN IS PROJECTED FOR REFERENCE COMPONENT
C AND COMPONENTS CONTAINED IN REDUCED HESSIAN.
C
C PURPOSE:
C   1)CONSTRUCT INFORMATION WHICH ARE REQUIRED TO START UPDATE
C   2)PUT T AND Q VECTORS ON LUIT3 AND LUIT5
C   3)RETURN STEP VECTOR IN VECT .GIVES ITERATION POINT TWO WHEN
C     ADDED TO PC PARAMETERS AT ITERATION POINT 2
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION REDINV(*),IBNDX(*),REF0(*),CMO(*)
      DIMENSION VECK(*),VECG(*),VECT(*),VECQ(*),VECX(*),G0(*)
      PARAMETER ( D0=0.0D0 , D1=1.0D0, D2=2.0D0, DM1=-1.0D0)
C
C  INFINP : FLAG(*)
C  INFTAP : LUIT1,?
C
#include <maxorb.h>
#include <infinp.h>
#include <inftap.h>
#include <infvar.h>
#include <inforb.h>
#include <infpri.h>
C
C
      LOGICAL HESFIX
C
C
      HESFIX = FLAG(56)
C
C READ IN MOLECULAR ORBITALS AT EXPANSION POINT (KAPPA=0)
C
      REWIND LUIT1
      CALL MOLLAB('OLDORB  ',LUIT1,lupri)
      CALL READT(LUIT1,NCMOT,CMO)
C
C READ IN KAPPA VECTOR WHICH REPRESENT THE NEW EXPANSION
C
      IF (NWOPT.GT.0)THEN
         CALL MOLLAB('ROTKAPPA',LUIT1,lupri)
         CALL READT(LUIT1,NWOPT,VECK(1+NCONF))
      END IF
C
C READ IN CSF COEFFICIENTS AT EXPANSION POINT (KAPPA=0)
C
C
      REWIND LUIT3
      IF (NCONF.GT.1) THEN
         CALL READT(LUIT3,NCONF,REF0)
      ELSE
         REF0(1) = D1
      END IF
C
C READ IN GRADIENT AT EXPANSION POINT IN VECT
C
      REWIND LUIT5
      IF (NCONF .GT. 1) THEN
         CALL READT(LUIT5,NCONF,VECT)
      ELSE
         VECT(1) = D0
      END IF
      IF (NWOPT .GT. 0) CALL READT(LUIT5,NWOPT,VECT(1+NCONF))
C
C
C TRANSFORM CSF COEFFICIENTS AND GRADIENT AT ITERATION POINT ONE
C TO PC REPRESENTATION
C
      IF (NCONF.GT.1) THEN
         CALL GCTOPC(NCONF,REF0,VECK,VECG,CSCAL,GSCAL)
      ELSE
         REF0(1) = D0
         VECK(1) = D0
         VECG(1) = D0
      END IF
C
C CONSTRUCT INFORMATION WHICH ARE REQUIRED BEFORE UPDATE CAN BE STARTED
C
C
      CALL UPDINI(LUIT7,REDINV,IBNDX,G0,VECX)
C     CALL UPDINI(LUIT7,REDINV,IBNDX,G0,WRK)
C
C CONSTRUCT T AND Q VECTORS AND WRITE ON LUIT3 AND LUIT5
C
C IF (FIXED HESSIAN)  step = - H0(inv) * g(1)
C
C ELSE                t(1) = H0(inv) * ( g(1) - g(0) )
C                     g(1) is saved in REF0, g(0) was saved in VECT.
C
      IF (.NOT.HESFIX) THEN
         CALL DCOPY(NVAR,VECG,1,REF0,1)
         CALL DAXPY(NVAR,DM1,VECT,1,VECG,1)
      END IF
      CALL UPDH0(LUIT7,IBNDX,REDINV,VECK,VECG,VECT,VECQ,G0,VECX)
C     CALL UPDH0(LUIT7,IBNDX,REDINV,VECK,VECG,VECT,VECQ,WRKX,WRK)
C
      IF (HESFIX) THEN
         CALL DSCAL(NVAR,DM1,VECT,1)
         GO TO 725
      END IF
C
C PUT GRADIENT AND PARAMETERS AT NEW POINT IN VECG AND REF0
C
      XQNOR = D1/DDOT(NVAR,VECT,1,VECK,1)
      CALL DSCAL(NVAR,XQNOR,VECQ,1)
      CALL DAXPY(NVAR,DM1,VECK,1,VECT,1)
      CALL WRITT(LUIT3,NVAR,VECT)
      CALL WRITT(LUIT5,NVAR,VECQ)
      IF (IPRSTAT .GT. 11) THEN
         WRITE(LUSTAT,'(/A)')' UPDST: VECT 1 WRITTEN ON LUIT3'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECT(I),I=1,NVAR)
         WRITE(LUSTAT,'(/A)')' UPDST: VECQ 1 WRITTEN ON LUIT5'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECQ(I),I=1,NVAR)
      END IF
C
C FIND BROYDN STEP AND WRITE STEP AS SECOND RECORD ON LUIT7
C
      CALL DCOPY(NVAR,REF0,1,VECG,1)
      XQG = DDOT(NVAR,VECQ,1,VECG,1)
      IF (IPRSTAT.GE.8) WRITE (LUSTAT,'(/A,F15.8)') ' < Q1 / F1 > =',XQG
      CALL DSCAL(NVAR,XQG,VECT,1)
      CALL UPDH0(LUIT7,IBNDX,REDINV,VECK,VECG,REF0,VECQ,G0,VECX)
      CALL DAXPY(NVAR,DM1,REF0,1,VECT,1)
      REWIND(LUIT7)
      READ  (LUIT7)
      CALL WRITT(LUIT7,NVAR,VECT)
      CALL WRITT(LUIT7,NVAR,VECG)
      IF (IPRSTAT .GT. 11) THEN
         WRITE(LUSTAT,'(/A)')' UPDST: STEP VECTOR 1 WRITTEN ON LUIT7'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECT(I),I=1,NVAR)
      END IF
C
C     STEP VECTOR IS IN VECK IN PC REPRESENTATION
C     ( AT EXPANSION POINT KAPPA AND PC COMPONENTS ARE ZERO )
C     COPY COEFFICIENTS AT NEW ITERATION POINT TO REF0
C
 725  CONTINUE
      CALL DCOPY(NVAR,VECK,1,REF0,1)
C
      RETURN
C     End of UPDST.
      END
C  /* Deck updh0 */
      SUBROUTINE UPDH0(LUIT7,IBNDX,REDINV,VECK,VECG,VECT,VECQ,WRKX,WRK)
C
C PURPOSE:
C   MULTIPLY ZEROTH ORDER APPROXIMATION TO THE INVERSE HESSIAN
C   ONTO VECG (RETURNED IN VECT) AND ONTO VECK (RETURNED IN VECQ)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION IBNDX(*),REDINV(NREDL,*),VECK(*),VECG(*),VECT(*),VECQ(*)
      DIMENSION WRKX(*),WRK(*)
C
#include <maxorb.h>
#include <ibndxdef.h>
#include <infpri.h>
#include <infopt.h>
#include <infvar.h>
#include <inftap.h>
C
C     DEFINE WORK SPACE
C
      KYVECG=1
      KYVECK=KYVECG+NREDL
      KRVECG=KYVECK+NREDL
      KRVECK=KRVECG+NREDL
C
C     MULTIPLY VECG AND VECK ON VECTORS ON LUIT3
C
      REWIND(LUIT3)
      DO 1000 I=1,NREDL
         IF (IBNDX(I).EQ.JBONDX) THEN
            NTOT=NWOPT
            IOFF=NCONF+1
         ELSE
            NTOT=NCONF
            IOFF=1
         END IF
         CALL READT(LUIT3,NTOT,WRKX)
         WRK(KYVECG-1+I)=DDOT(NTOT,VECG(IOFF),1,WRKX,1)
         WRK(KYVECK-1+I)=DDOT(NTOT,VECK(IOFF),1,WRKX,1)
 1000 CONTINUE
C
C     MULTIPLY REDINV ON YVECG AND YVECK
C     (we have used YVECK (RVECK) follows YVECG (RVECG) in core)
C
      CALL DGEMM('N','N',NREDL,2,NREDL,1.D0,
     &           REDINV,NREDL,
     &           WRK(KYVECG),NREDL,0.D0,
     &           WRK(KRVECG),NREDL)
      IF (IPRSTAT .GE. 10) THEN
         WRITE(LUSTAT,'(/A)')' (UPDH0) YVECG'
         WRITE(LUSTAT,'(1X,5F15.8)')(WRK(KYVECG-1+I),I=1,NREDL)
         WRITE(LUSTAT,'(/A)')' YVECK'
         WRITE(LUSTAT,'(1X,5F15.8)')(WRK(KYVECK-1+I),I=1,NREDL)
         WRITE(LUSTAT,'(/A)')' RVECG'
         WRITE(LUSTAT,'(1X,5F15.8)')(WRK(KRVECG-1+I),I=1,NREDL)
         WRITE(LUSTAT,'(/A)')' RVECK'
         WRITE(LUSTAT,'(1X,5F15.8)')(WRK(KRVECK-1+I),I=1,NREDL)
      END IF
C
C READ IN DIAGONAL INVERSE HESSIAN
C
      REWIND(LUIT7)
      CALL READT(LUIT7,NVAR,WRKX)
C
C MULTIPLY ZEROTH ORDER INVERSE HESSIAN MULTIPLIED ON VECG (RETURNED
C IN VECT) AND ON VECK (RETURNED IN VECQ)
C
      DO 1030 I=1,NVAR
         VECT(I)=VECG(I)*WRKX(I)
         VECQ(I)=VECK(I)*WRKX(I)
 1030 CONTINUE
      IF (IPRSTAT .GT. 15) THEN
         WRITE(LUSTAT,'(/A)')' T VECTOR AFTER DIAGONAL INVERSE HESSIAN'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECT(I),I=1,NVAR)
         WRITE(LUSTAT,'(/A)')' Q VECTOR AFTER DIAGONAL INVERSE HESSIAN'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECQ(I),I=1,NVAR)
      END IF
      REWIND LUIT3
      REWIND LUIT5
      DO 1040 I=1,NREDL
         IF (IBNDX(I).EQ.JBONDX )THEN
            NTOT=NWOPT
            IOFF=NCONF+1
         ELSE
            NTOT=NCONF
            IOFF=1
         END IF
         CALL READT(LUIT5,NTOT,WRKX)
         CALL DAXPY(NTOT,-WRK(KYVECG-1+I),WRKX,1,VECT(IOFF),1)
         CALL DAXPY(NTOT,-WRK(KYVECK-1+I),WRKX,1,VECQ(IOFF),1)
         FACT=WRK(KRVECG-1+I)-DDOT(NTOT,WRKX,1,VECG(IOFF),1)
         FACQ=WRK(KRVECK-1+I)-DDOT(NTOT,WRKX,1,VECK(IOFF),1)
         CALL READT(LUIT3,NTOT,WRKX)
         CALL DAXPY(NTOT,FACT,WRKX,1,VECT(IOFF),1)
         CALL DAXPY(NTOT,FACQ,WRKX,1,VECQ(IOFF),1)
 1040 CONTINUE
      IF (IPRSTAT .GT. 14) THEN
         WRITE(LUSTAT,'(/A)')' T VECTOR AFTER UPDH0'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECT(I),I=1,NVAR)
         WRITE(LUSTAT,'(/A)')' Q VECTOR AFTER UPDH0'
         WRITE(LUSTAT,'(1X,5F15.8)')(VECQ(I),I=1,NVAR)
      END IF
      RETURN
C     End of UPDH0.
      END
C  /* Deck updini */
      SUBROUTINE UPDINI(LUIT7,REDINV,IBNDX,G0,WRK)
C
C 17-Oct-1985 PJ
C
C PURPOSE:
C  1) STORE THE INVERSE DIAGONAL HESSIAN AS FIRST RECORD ON LUIT7
C  2) STORE THE INVERSE DIAGONAL HESSIAN MULTIPLIED WITH THE TRIAL
C     VECTORS ON LUIT3  ON LUIT5
C  3) ADD TO INVERSE REDUCED HESSIAN THE INVERSE DIAGONAL
C     HESSIAN
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION REDINV(*),IBNDX(*),G0(*),WRK(*)
C
#include <ibndxdef.h>
C
      PARAMETER ( D0=0.0D0 , D1=1.0D0 ,D2=2.0D0 )
      PARAMETER ( DTEST=1.0D-3 )
C
C   INFINP : FLAG(*)
C   INFTAP : LUIT1,?
C
#include <maxorb.h>
#include <infinp.h>
#include <infvar.h>
#include <infopt.h>
#include <inftap.h>
#include <infpri.h>
C
C
      LOGICAL HESIVD
      CHARACTER*8 TABLE(3)
C
C --  data:
C
      DATA TABLE/'xxxxxxxx','RESTART ','LREDUCED'/
C
C
      HESIVD = FLAG(57)
C
C     RECOVER REDUCED (PROJECTED) HESSIAN FROM PREVIOUS MACRO ITERATION
C
      REWIND LUIT1
      CALL MOLLAB(TABLE(2),LUIT1,lupri)
      IF (HESIVD) THEN
         NREDL = 1
         READ (LUIT1)
      ELSE
         READ (LUIT1) DUM,DUM,DUM,DUM,DUM,IDUM,
     *                MREDL,(IBNDX(I),I=1,MREDL)
         IF (MREDL.NE.NREDL) THEN
            WRITE (LUSTAT,'(/1X,A,2I5/1X,A)')
     *      'UPDINI: NREDL on LUIT1 inconsistent with NREDL in common:',
     *      MREDL,NREDL,
     *      '--> INVERSE DIAGONAL USED AS APPROXIMATE INVERSE HESSIAN.'
C           CALL QTRACE(LUSTAT)
C           CALL QUIT('***ERROR*** UPDINI: NREDL on LUIT1 inconsistent')
            NREDL = 1
         END IF
      END IF
      NNREDL=NREDL*(NREDL+1)/2
      READ (LUIT1) DUM,EACT
      IF (IPRSTAT .GE. 5)
     *   WRITE (LUSTAT,'(/A,F25.15)') ' UPDINI: Active energy:',EACT
      IF (HESIVD) GO TO 525
      CALL MOLLAB(TABLE(3),LUIT1,lupri)
      CALL READT(LUIT1,NNREDL,WRK(1))
      IF (IPRSTAT .GE. 8) THEN
         WRITE (LUSTAT,*) 'Reduced L matrix, dimension ',NREDL
         CALL OUTPAK(WRK,NREDL,1,LUSTAT)
      END IF
C
C     PACK REDUCED HESSIAN IN TWO DIMENSIONAL ARRAY
C     (PROJECT REFERENCE VECTOR OUT)
C
      IJ=0
      DO 500 I=1,NREDL
         DO 501 J=1,I
            IJ=IJ+1
            KIJ=(I-1)*NREDL+J
            KJI=(J-1)*NREDL+I
            REDINV(KIJ)=WRK(IJ)
            REDINV(KJI)=WRK(IJ)
 501     CONTINUE
         REDINV((I-1)*NREDL+1)=D0
         REDINV(I)=D0
 500  CONTINUE
      REDINV(1)=D1
      IF (IPRSTAT .GE. 10) THEN
         WRITE(LUSTAT,'(/A)')' REDUCED HESSIAN REF0 PROJECTED OUT'
         CALL OUTPUT(REDINV,1,NREDL,1,NREDL,NREDL,NREDL,1,LUSTAT)
      END IF
C
C     INVERS OF REDUCED HESSIAN
C
      CALL DGEINV(NREDL,REDINV,REDINV,WRK(1),WRK(1+NREDL),INFO)
      IF (INFO .NE. 0) THEN
         WRITE (LUSTAT,'(//A/T6,A,I5/T6,A/)') ' *** ERROR (UPDINI) ***',
     *      ' Reduced Hessian singular, INFO from DGEINV =',INFO,
     *      ' Only diagonal elements used for H0.'
         NREDL = 1
      END IF
 525  CONTINUE
      REDINV(1)=D0
      IF (IPRSTAT .GE. 10) THEN
         WRITE(LUSTAT,'(/A)')
     *      ' REDUCED INVERSE HESSIAN, REF0 PROJECTED OUT'
         CALL OUTPUT(REDINV,1,NREDL,1,NREDL,NREDL,NREDL,1,LUSTAT)
      END IF
C
C     READ IN DIAGONAL L-MATRIX (ACTIVE ENERGY IS NOT SUBTRACTED)
C
      REWIND LUIT2
      IF (NWOPT.GT.0) THEN
         CALL MOLLAB('ORBDIAG ',LUIT2,lupri)
         CALL READT(LUIT2,NWOPT,G0(NCONF+1))
      END IF
      IF (NCONF.GT.1)THEN
         CALL MOLLAB('CIDIAG2 ',LUIT2,lupri)
         CALL READT(LUIT2,NCONF,G0(1))
      ELSE
         G0(1) = D0
      END IF
      IF (IPRSTAT .GT. 12) THEN
         WRITE(LUSTAT,'(/A)') ' Diagonal L matrix'
         WRITE(LUSTAT,'(1X,5F15.8)')(G0(I),I=1,NVAR)
      END IF
C PHPMAERKE : Oct90/hjaaj: implement PHP here ?
C
C     CONSTRUCT INVERSE DIAGONAL HESSIAN AND STORE ON LUIT7
CPHPMAERKE: PHP could be used here 901028/hjaaj
C
      DO 1025 J=1,NCONF
         D=(G0(J)-EACT)*D2
         IF (ABS(D).LT.DTEST) THEN
            G0(J)=D1/SIGN(DTEST,D)
         ELSE
            G0(J)=D1/D
         END IF
 1025 CONTINUE
      DO 1026 J=(1+NCONF),NVAR
         D=G0(J)
         IF (ABS(D).LT.DTEST) THEN
            G0(J)=D1/SIGN(DTEST,D)
         ELSE
            G0(J)=D1/D
         END IF
 1026 CONTINUE
      IF (IPRSTAT .GT. 11) THEN
         WRITE(LUSTAT,'(/A)')'DIAGONAL INVERSE  HESSIAN'
         WRITE(LUSTAT,'(1X,5F15.8)')(G0(I),I=1,NVAR)
      END IF
      REWIND LUIT7
      CALL WRITT(LUIT7,NVAR,G0)
C
C     LUIT5 IS USED TO SAVE INVERSE DIAGONAL HESSIAN MULTIPLIED BY
C     THE TRIAL VECTORS STORED ON LUIT3
C     VECTORS ON LUIT3 DENOTED Y
C     VECTORS ON LUIT5 DENOTED X
C
      REWIND LUIT3
      REWIND LUIT5
      DO 1000 I=1,NREDL
         IF (IBNDX(I).EQ.JBONDX) THEN
            NTOT=NWOPT
            IOFF=NCONF
         ELSE
            NTOT=NCONF
            IOFF=0
         END IF
         CALL READT(LUIT3,NTOT,WRK(1))
         DO 1010 J=1,NTOT
            WRK(J)=WRK(J)*G0(IOFF+J)
 1010    CONTINUE
         CALL WRITT(LUIT5,NTOT,WRK(1))
 1000 CONTINUE
C
C     ADD <X(I)/Y(J)> CONTRIBUTION TO INVERSE REDUCED HESSIAN
C     (i.e. to REDINV(I,J))
C
      REWIND LUIT3
      DO 1030 I=1,NREDL
         IF (IBNDX(I).EQ.JBONDX) THEN
            NTOT=NWOPT
         ELSE
            NTOT=NCONF
         END IF
         CALL READT(LUIT3,NTOT,WRK)
         REWIND LUIT5
         DO 1040 J=1,I
            IF (IBNDX(I).EQ.IBNDX(J)) THEN
               CALL READT(LUIT5,NTOT,G0)
               IJ=(J-1)*NREDL+I
               JI=(I-1)*NREDL+J
               REDINV(JI)=REDINV(JI)+DDOT(NTOT,WRK(1),1,G0,1)
               REDINV(IJ)=REDINV(JI)
            ELSE
               READ(LUIT5)
            END IF
 1040    CONTINUE
 1030 CONTINUE
      IF (IPRSTAT .GE. 10) THEN
         WRITE(LUSTAT,'(/A)')
     *      ' REDUCED INVERSE HESSIAN WITH YX CONTRIBUTION ADDED'
         CALL OUTPUT(REDINV,1,NREDL,1,NREDL,NREDL,NREDL,1,LUSTAT)
      END IF
      RETURN
C     End of UPDINI.
      END
