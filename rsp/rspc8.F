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
#ifdef revlog
!===========================================================================
!921011-pj+hjaaj: EFPOLL added C6IFC test
!920722-Hinne Hettema + HJAaJ
!inserted Hinnes changes in EFPOLL,POLIM (DIFMAX test of convergence,
!changed print levels, new format for LURSP7, fixed CFREQ bug)
!===========================================================================
#endif
C  /* Deck rspc6 */
      SUBROUTINE RSPC6(CMO,UDV,PV,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
#include <implicit.h>
C
      LOGICAL IDONE
#include <iratdef.h>
C
      CHARACTER*8 BLANK
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(LWRK)
C
#include <priunit.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <infc6.h>
#include <inflr.h>
#include <infdim.h>
#include <inforb.h>
C
      PARAMETER ( ZERO = 0.0D0  , DM1 = -1.0D0 , BIG = 1.0D10 )
      PARAMETER ( DLRTST = 1.0D-6, BLANK = '        ' )
C
C DETERMINE POLARIZABILITIES AT IMAGINARY FREQUENCIES.
C
C EXPAND V.(E[2]-OMEGA*S[2])**(-1).Vt = alpha AS CAUCHY-TYPE SERIES
C AND FIND A BASIS OF PSEUDO-STATES BY RECURSIVE SOLUTION OF
C LINEAR EQUATIONS ARISING AT EACH ORDER OF OMEGA.
C FINALLY DIAGONALISE WTHIN REDUCED BASIS AND OBTAIN 16
C POLARISABILITIES AS SUMS OF PDOSD CONTRIBUTIONS.
C SIZE OF PROBLEM IS GOVERNED BY A CUTOFF ON THE NORM OF THE
C PSEUDOSTATE VECTORS S[2]*LAMBDA[N] OR BY A PRESET MAXIMUM
C NUMBER OF MOMENTS IN THE EXPANSION.
C
C NUMBER OF VECTORS OF SYMMETRY KSYMOP
C
      NGDNUM = NGPC6(KSYMOP)
      NGRDMX = MAX(16,NGRID)
C
      KREDE  = 1
      KREDS  = KREDE  + MAXRM*MAXRM
      KIBTYP = KREDS  + MAXRM*MAXRM
      KEIVAL = KIBTYP + MAXRM
      KRESID = KEIVAL + MAXRM
      KEIVEC = KRESID + MAXRM
      KREDGD = KEIVEC + MAXRM*MAXRM
      KGD    = KREDGD + MAXRM
      KSOL   = KGD    + KZYVAR*NGDNUM
      KXMOM  = KSOL   + KZYVAR*NGDNUM
      KOLDAL = KXMOM  + MAXRM*NGDNUM
      KWRK1  = KOLDAL + NGDNUM*(NGRDMX+1)
      LWRK1  = LWRK   - KWRK1
C
      LBVMAX = MAX(KZCONF,KZYWOP)
      IF (LWRK1.LT.2*KZYVAR+LBVMAX) THEN
         WRITE (LUERR,9100) LWRK1,2*KZYVAR+LBVMAX
         CALL QTRACE(LUERR)
         CALL QUIT('RSPC6 ERROR, INSUFFICIENT SPACE TO SOLVE LIN. EQ.S')
      ENDIF
 9100 FORMAT(/' RSPC6, work space too small for 2.5 (z,y)-vectors',
     *       /'        had:',I10,', need more than:',I10)
C
C WORK SPACE FOR RSPEVE
C
      KBVECS = KWRK1
      KWRKE  = KBVECS + KZYVAR
      LWRKE  = LWRK   - KWRKE
      IF (LWRKE.LT.0) CALL ERRWRK('RSPC6',KWRKE-1,LWRK)
C
      KZRED  = 0
      KZYRED = 0
      IPRRSP = IPRC6
      MAXIT  = MAXITC
C
      CALL DZERO(WRK(KOLDAL),NGDNUM*(NGRDMX+1))
C
      IF ( IPRRSP.GT.83) THEN
         WRITE(LUPRI,*)' THCC6,IPRC6,MAXITC',THCC6,IPRC6,MAXITC
         WRITE(LUPRI,*)' NGDNUM',NGDNUM
         WRITE(LUPRI,*)' MAXMOM',MAXMOM
      END IF
C
C     Call RSPCTL to solve linear set of response equations
C
      KGDOFF = KGD
      KSOLOF = KSOL
      DO 600 IOP = 1,NGDNUM
         WRITE(LUPRI,'(//A)') ' Solving RSPC6 for'
         WRITE(LUPRI,'(A,A/A,I4)')
     *      ' OPERATOR TYPE:    ',LBLC6(KSYMOP,IOP),
     *      ' OPERATOR SYMMETRY:',KSYMOP
         CALL GETGPV(LBLC6(KSYMOP,IOP),FC,FV,CMO,UDV,PV,XINDX,ANTSYM,
     *               WRK(KGDOFF),LWRK1)
         KGDOFF = KGDOFF + KZYVAR
 600  CONTINUE
      NTOTAL = KZYVAR*NGDNUM
      DO 610 I = 1,NTOTAL
         WRK(KSOL-1+I) = WRK(KGD-1+I)
 610  CONTINUE
C     CALL DCOPY(KZYVAR*NGDNUM,WRK(KGD),1,WRK(KSOL),1)
      KEXSIM = 1
      KEXCNV = 1
      ITC6   = 0
      C6NEW  = BIG
      C6DIF  = BIG
      THCRSP = THCC6
      IDONE  = .FALSE.
 100     CONTINUE
         WRITE(LUPRI,'(//A,I3/)')
     &      ' **** RSPC6 moment iteration no.',ITC6
         KSOLOF = KSOL
         KGDOFF = KGD
         DO 200 IOP = 1,NGDNUM
               WRITE(LUPRI,'(A,A/A,I4)')
     *          ' OPERATOR TYPE:    ',LBLC6(KSYMOP,IOP),
     *          ' OPERATOR SYMMETRY:',KSYMOP
               WRK(KEIVAL) = ZERO
               IF (IPRRSP.GT.110) THEN
                  WRITE(LUPRI,'(/A,I5)')
     *           ' GRADIENT VECTOR FOR CAUCHY SERIES NUMBER',ITC6
                  CALL OUTPUT(WRK(KSOLOF),1,KZYVAR,1,1,KZYVAR,1,1,LUPRI)
               END IF
               XLRTST = DDOT(KZYVAR,WRK(KSOLOF),1,WRK(KSOLOF),1)
               IF( SQRT(XLRTST).LE.DLRTST) THEN
                  WRITE(LUPRI,'(/A)')
     *            ' *** RIGHT HAND SIDE OF LINEAR EQUATION IS ZERO'
     *            ,' *** THE LINEAR EQUATION IS NOT SOLVED         '
                  GO TO 200
               END IF
               CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
     *            .TRUE.,LBLC6(KSYMOP,IOP),BLANK,WRK(KSOLOF),
     *            WRK(KREDGD),WRK(KREDE),WRK(KREDS),
     *            WRK(KIBTYP),WRK(KEIVAL),WRK(KRESID),WRK(KEIVEC),
     *            XINDX,WRK(KWRK1),LWRK1)
C              CALL RSPCTL(CMO,UDV,PV,FC,FV,FCAC,H2AC,
C                       LINEQ,GD,REDGD,REDE,REDS,
C    *                  IBTYP,EIVAL,EIVEC,XINDX,WRK,LWRK)
C
               RESTLR = .TRUE.
               IF ( IPRRSP.GT.5) THEN
                  WRITE(LUPRI,'(/A,G15.8)')
     *            ' threshold for linear equations THCRSP:',THCRSP
               END IF
               CALL RSPEVE(WRK(KIBTYP),WRK(KEIVAL),WRK(KEIVEC),
     *                  WRK(KBVECS),WRK(KWRKE),1,0)
C              CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
               IF (IPRRSP.GT.100) THEN
                  WRITE(LUPRI,'(/A,I5,/A,A)')
     *           ' SOLUTION VECTOR FOR lamda(k) k =',ITC6,
     *           ' OPERATOR ',LBLC6(KSYMOP,IOP)
                  CALL OUTPUT(WRK(KBVECS),1,KZYVAR,1,1,KZYVAR,1,1,LUPRI)
               END IF
               GDNORM =  DDOT(KZYVAR,WRK(KSOLOF),1,WRK(KBVECS),1)
               IF (IPRRSP.GT.42) THEN
                  WRITE(LUPRI,'(/A,1G15.8)')
     *            ' GDNORM ',GDNORM
               END IF
               IF (IPRRSP.GT.81) THEN
                  CAUCMO =  DDOT(KZYVAR,WRK(KGDOFF),1,WRK(KBVECS),1)
                  WRITE(LUPRI,'(/A,I5,A,G15.8)')
     *           'CAUCHY MOMENT itc6=:',ITC6,' CAUCMO = ',CAUCMO
               END IF
C
C CONSTRUCT  S[2]*X WHERE X IS THE SOLUTION VECTOR
C
               CALL DZERO(WRK(KSOLOF),KZYVAR)
               IF (KZCONF.GT.0) THEN
                  CALL RSPSLI(1,0,WRK(KBVECS+KZVAR),DUMMY,
     *                     UDV,WRK(KSOLOF),XINDX,WRK(KWRKE),LWRKE)
                  CALL DSWAP(KZVAR,WRK(KSOLOF),1,WRK(KSOLOF+KZVAR),1)
                  CALL DSCAL(KZYVAR,DM1,WRK(KSOLOF),1)
                  CALL RSPSLI(1,0,WRK(KBVECS),DUMMY,
     *                        UDV,WRK(KSOLOF),XINDX,WRK(KWRKE),LWRKE)
               END IF
               IF (KZWOPT.GT.0) THEN
                   DO 191 II = 1,KZWOPT
                      WRK(KBVECS+KZVAR-1+II) =
     *                   WRK(KBVECS+KZVAR+KZCONF-1+II)
 191               CONTINUE
C                  CALL DCOPY(KZWOPT,WRK(KBVECS+KZVAR+KZCONF),1,
C    *                    WRK(KBVECS+KZVAR),1)
                 CALL RSPSLI(0,1,WRK(KBVECS+KZCONF),WRK(KBVECS+KZCONF),
     *                    UDV,WRK(KSOLOF),XINDX,WRK(KWRKE),LWRKE)
               END IF
               KGDOFF = KGDOFF + KZYVAR
               KSOLOF = KSOLOF + KZYVAR
 200     CONTINUE
         ITC6 = ITC6 + 1
         IF (ITC6.GT.MAXMOM) THEN
            WRITE(LUPRI,'(A,A,/A,I4,A,I4)')
     *      ' *** WARNING: Stopping ITC6 iterations, but ',
     *      ' polarizabilities not converged.',
     *      ' ITC6 = ', ITC6, '   MAXMOM = ', MAXMOM
            NWARN = NWARN + 1
            GO TO 400
         ENDIF
         CALL POLIM(WRK(KIBTYP),WRK(KREDS),WRK(KREDE),
     *           WRK(KXMOM),WRK(KEIVAL),WRK(KEIVEC),WRK(KGD),
     *           WRK(KOLDAL),UDV,XINDX,WRK(KWRK1),LWRK1,IDONE)
         IF (IDONE) THEN
            WRITE(LUPRI,'(A,/A,I4,/A)') '  ',
     *      ' Number of ITC6 iterations used : ', ITC6,
     *      ' ============================================='
            GOTO 400
         ENDIF
         GO TO 100
C loop 100 finish
 400  CONTINUE
C
      IDONE = .TRUE.
      IF ( MOD(ITC6,2).EQ.0 ) ITC6 = ITC6 - 1
      CALL POLIM(WRK(KIBTYP),WRK(KREDS),WRK(KREDE),
     *           WRK(KXMOM),WRK(KEIVAL),WRK(KEIVEC),WRK(KGD),
     *           WRK(KOLDAL),UDV,XINDX,WRK(KWRK1),LWRK1,IDONE)

      RESTLR = .FALSE.
C
C *** end of RSPC6 --
C
      RETURN
      END
C
C  /* Deck polim */
      SUBROUTINE POLIM(IBTYP,REDS,REDE,XMOM,EIVAL,EIVEC,
     *                 GP,OLDAL,UDV,XINDX,WRK,LWRK,IDONE)
C
C
C     IDONE IS A FLAG. IDONE=.false. THE QUADRATURE IS CARRIED OUT.
C                             .true. THE CAUCHY PROCEDURE HAS
C                             CONVERGED AND THE PSEUDOSTATES ARE USED
C                             TO EVALUATE A SET OF REAL-FREQUENCY
C                             POLARISABILITIES, OSCILLATOR SUMS
C                             AND THE W4 COEFFICIENT.
C
C OUTPUT EIVAL(KZRED) AND XMOM(KZRED,NOP)
C
#include <implicit.h>
C
#include <priunit.h>
#include <inforb.h>
#include <infdim.h>
#include "rspprp.h"
#include <infc6.h>
#include <infpri.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      LOGICAL IDONE
      DIMENSION XMOM(MAXRM,*)
      DIMENSION IBTYP(*),REDS(*),REDE(*),UDV(*),XINDX(*),WRK(LWRK)
      DIMENSION EIVAL(*),EIVEC(*),GP(KZYVAR,*)
      DIMENSION OLDAL(*)
C
      PARAMETER ( ZERO = 0.0D0 ,  D1 = 1.0D0 ,
     *D2 = 2.0D0 , D3 = 3.0D0, D10= 10.0D0, DUMMY = 1.0D20 )
      PARAMETER ( COMPLX = 1.0D7 )
C
      IF (IPRRSP.GE.15)
     *   WRITE (LUPRI,'(A//A)') '1',' *** Output from POLIM ***'
      CALL PPRST(IBTYP,REDS,REDE,UDV,XINDX,WRK,LWRK)
      CALL RSPRED(2,.FALSE.,0,IBTYP,DUMMY,DUMMY,REDE,REDS,
     *            EIVAL,EIVEC,WRK,WRK,UDV,WRK,
     *            XINDX,WRK,LWRK)
C
C CALCULATE EIGENVECTORS AND TRANSITION MOMENTS
C
C ALLOCATE WORK SPACE
C
C MAXIMUM NUMBER OF TRIAL VECTORS
C
      NSIMMA =  (LWRK-KZYVAR)/KZYVAR
      IF (NSIMMA.GT.KZRED) THEN
         NSIM = KZRED
      ELSE
         NSIM = NSIMMA
      END IF
      IF (NSIM .LE. 0) THEN
         WRITE (LUPRI,*) 'ERROR in POLIM, NSIM = 0'
         WRITE (LUPRI,*) 'More memory required.'
         CALL QUIT('Insufficient memory in POLIM')
      END IF
      KBVECS = 1
      KWRK2  = KBVECS + NSIM*KZYVAR
      LWRK2  = LWRK   - KWRK2
      IF (LWRK2.LT.0) CALL ERRWRK('POLIM',KWRK2-1,LWRK)
C
      DO 500 ISIM = 1,KZRED,NSIM
         NBX = MIN( NSIM,(KZRED+1-ISIM) )
         CALL RSPEVE(IBTYP,EIVAL,EIVEC,WRK(KBVECS),
     *               WRK(KWRK2),NBX,(ISIM-1))
C        CALL RSPEVE(IBTYP,EIVAL,EIVEC,BVECS,WRK,NBX,IBOFF)
         DO 550 INUM = 1,NBX
            DO 560 IOP = 1,NGPC6(KSYMOP)
               XMOM(ISIM-1+INUM,IOP) = DDOT(KZYVAR,GP(1,IOP),1,
     *                     WRK(KBVECS+(INUM-1)*KZYVAR),1)
 560        CONTINUE
 550     CONTINUE
 500  CONTINUE
      IF (IPRRSP.GE.30) THEN
         DO 600 IOP = 1,NGPC6(KSYMOP)
            WRITE(LUPRI,'(/A/2A/A,I4)')
     *       ' Test output of pseudo spectrum for',
     *       ' OPERATOR TYPE:    ',LBLC6(KSYMOP,IOP),
     *       ' OPERATOR SYMMETRY:',KSYMOP
            WRITE(LUPRI,'(/2A)')
     *       ' STATE NO:  *TRANSITION MOMENT:     ENERGY(au)'
            DO 700 IST = 1,KZRED
              WRITE(LUPRI,'(1X,I10,2F15.5)')IST,XMOM(IST,IOP),EIVAL(IST)
 700        CONTINUE
 600     CONTINUE
      END IF
C
C     Print of pseudo-spectrum
C
      IF (IDONE) THEN
         LIMPRI =  5
      ELSE
         LIMPRI = 20
      END IF
      IF (IPRRSP .GE. LIMPRI) THEN
         DO 800 IOP = 1,NGPC6(KSYMOP)
            IF (LBLC6(KSYMOP,IOP)(2:4) .EQ. 'DIP') THEN
               CALL C6PRPS(IDONE,LBLC6(KSYMOP,IOP),KSYMOP,
     &                     KZRED,XMOM,EIVAL)
            END IF
  800    CONTINUE
      END IF
C
      NVEL = 0
      NCAR = 0
      NSPH = 0
      DO 900 IOP = 1,NGPC6(KSYMOP)
         IF(LBLC6(KSYMOP,IOP)(2:5).EQ.'DIPL' ) NCAR = NCAR + 1
         IF(LBLC6(KSYMOP,IOP)(3:6).EQ.'QUAD' ) NCAR = NCAR + 1
         IF(LBLC6(KSYMOP,IOP)(4:6).EQ.'MOM'  ) NCAR = NCAR + 1
         IF(LBLC6(KSYMOP,IOP)(2:5).EQ.'DIPV' ) NVEL = NVEL + 1
         IF(LBLC6(KSYMOP,IOP)(1:2).EQ.'SM'   ) NSPH = NSPH + 1
 900  CONTINUE
      IF (IPRRSP.GT.15) THEN
         WRITE(LUPRI,'(/A,I3)')
     *   ' NUMBER OF CARTESIAN Cn LENGTH OPERATORS   : NCAR=',NCAR
         WRITE(LUPRI,'(/A,I3)')
     *   ' NUMBER OF SPHERICAL Cn LENGTH OPERATORS   : NSPH=',NSPH
         WRITE(LUPRI,'(/A,I3)')
     *   ' NUMBER OF C6 VELOCITY OPERATORS : NVEL=',NVEL
      END IF
      NGRDMX = MAX(16,NGRID)
      NMOMC  = 20
C --  NMOMC is the number of Cauchy moments computed
C
      KGRID  = 1
      KALPHA = KGRID  + NGRDMX+1
      KCAUCH = KALPHA + (NGRDMX+1) * NGPC6(KSYMOP) * NGPC6(KSYMOP)
      KWRK1  = KCAUCH + (NMOMC+1) *  NGPC6(KSYMOP) * NGPC6(KSYMOP)
      LWRK1  = LWRK   - KWRK1
      IF (NCAR.EQ.NGPC6(KSYMOP).OR.(NSPH.EQ.NGPC6(KSYMOP))) THEN
         CALL EFPOLL(IDONE,WRK(KGRID),WRK(KALPHA),WRK(KCAUCH),NGRDMX,
     *               NMOMC,NGPC6(KSYMOP),EIVAL,XMOM,OLDAL)
      ELSE IF (NVEL.EQ.NGPC6(KSYMOP)) THEN
         CALL EFPOLV(IDONE,WRK(KGRID),WRK(KALPHA),WRK(KCAUCH),NGRDMX,
     *               NMOMC,NGPC6(KSYMOP),EIVAL,XMOM)
      ELSE
         WRITE(LUPRI,*)' BOTH DIPLEN AND DIPVEL OPERATORS SPECIFIED'
         WRITE(LUPRI,*)' ILLEGAL OPTION. NCAR=',NCAR,' NVEL=',NVEL,
     *                  'NSPH= ',NSPH
         CALL QUIT(' POLIM: Illegal option, operator types are mixed')
      END IF
      RETURN
      END
C  /* Deck efpoll */
      SUBROUTINE EFPOLL(IDONE,GRIDSQ,ALPHA,CAUCHY,NGRDMX,NMOMC,NSYMOP,
     *                  EIVAL,XMOM,OLDAL)
c------------------------------------------------------------------------
c This subroutine calculates the frequency dependent polarisabilities on
c a precomputed grid and optionally on real frequencies.
c The default is the Gauss Chebychev grid, the alternative is the Gauss
c Legendre grid. Details can be found in:
c default:      W.Rijks and P.E.S.Wormer JCP Vol.88, 5704 (1988).
c G-Legendre:   Amos et al, J.Phys.Chem. Vol.89, 2186 (1985)
c
c We predefine 10 Cauchy moments. On change also change the memory layout
c in the C8 module.  Now NMOMC input parameter (900718/hjaaj).
c
c Note:   loops to NGRID go from 0(!)  to NGRID, because in ALPHA(I,J,0)
c we keep the static polarisability. We do not need this in the computa-
c tion of the Cn coefficients.
c
c When this subroutine is called, we assume that we only have operators
c of one type present, so that the interface file to the coupling program
c is consistent. This is ensured by the IF THEN ELSE statement around the
c call.
c
c We require too much space for the polarisability, but it facilitates
c addressing and allows at the same time for properties which are not
c symmetric in the indices.
c------------------------------------------------------------------------
c
#include <implicit.h>
c
#include <pi.h>
#include <dummy.h>
c
#include <priunit.h>
#include <infrsp.h>
#include "rspprp.h"
#include <infc6.h>
#include <wrkrsp.h>
C
#include <inforb.h>
#include <infpri.h>
c
      DIMENSION EIVAL(*),XMOM(MAXRM,*)
      DIMENSION GRIDSQ(0:NGRDMX),ALPHA(1:NSYMOP,1:NSYMOP,0:NGRDMX)
      DIMENSION CAUCHY(1:NSYMOP,1:NSYMOP,1:NMOMC)
      DIMENSION OLDAL(NSYMOP,0:NGRDMX)
      DIMENSION OMEGA(0:16)
      CHARACTER*8 IC,JC
      LOGICAL IDONE,END
      LOGICAL FIRST, LC6IFC
      DATA FIRST /.TRUE./
      SAVE FIRST
      DATA ITYPE /3/
      DATA THRESHW /1.D-08/
c
      DATA OMEGA/0.0D0,
     *           0.0000011354D0,0.0000324954D0,0.0002074939D0,
     *           0.0007766098D0,0.0022314002D0,0.0055272185D0,
     *           0.0125684274D0,0.0273216537D0,0.0585616091D0,
     *           0.1273031181D0,0.2894765239D0,0.7170385627D0,
     *           2.0602366551D0,7.7110713698D0,49.2377677794D0,
     *           1409.1898779477D0/
c -- parameters
      PARAMETER ( ZERO = 0.0D0 , D1 = 1.0D0, D2 = 2.0D0 )
c
c compute squared gridpoints
c
      IF(GSLEGN) THEN
c -- use Gauss Legendre grid from old C6 module
         NGRID=16
         DO 5 I=0,NGRID
            GRIDSQ(I) = OMEGA(I)
   5     CONTINUE
      ELSE
c-- compute gridpoints according to Gauss Chebychev scheme
         GRIDSQ(0) = ZERO
         N= NGRID
         NN=2*N
         FAC=PI/(2*NN)
         DO 15 I=1,N
            X=FAC*(2*I-1)
            GRIDSQ(N-I+1)=(D1/TAN(X))**2
  15     CONTINUE
      ENDIF
c
c write header output file when first call
c
      LURSP7 = -1
      CALL GPOPEN(LURSP7,'RESPONSE.C8','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      IF ( FIRST .AND. C6IFC ) THEN
         WRITE(LURSP7,'(A)')    '* SIRIUS OUTPUT FILE'
         WRITE(LURSP7,'(A,I8)') '* NOCC',NOCCT
         WRITE(LURSP7,'(A,I8)') '* NACT',NASHT
         WRITE(LURSP7,'(A,I8)') '* NVIR',NSSHT
         WRITE(LURSP7,'(A)')    '* GROUP  NONE'
         WRITE(LURSP7,'(A,I8)') '* NOMEGA',NGRID
         DO I=0,NGRID+1,5
            JMAX = MIN(I+4,NGRID)
            WRITE(LURSP7,'(A,5(1X,F11.6))')
     &                   '*',(SQRT(GRIDSQ(J)),J=I,JMAX)
         ENDDO
         FIRST = .FALSE.
      ENDIF
c
c initialise alpha and cauchy moments to zero
c
      CALL DZERO(ALPHA,NSYMOP*NSYMOP*(NGRID+1))
      CALL DZERO(CAUCHY,NSYMOP*NSYMOP*NMOMC)
c
      END   = IDONE
      IDONE = .TRUE.
c
c Compute the polarisability and cauchy moments from the effective spectrum.
c We only need the lower triangle, hence we have j<= i in the second loop
c
c We use that the LBLC6 has the following format: "SM02-02" for
c the l=2,m=-2 multipole matrix to extract the l,m
c
      IF (IPRRSP.GE.1) WRITE(LUPRI,'(/1X,A,/1X,A)')
     &   'Static polarisabilities (length)',
     &   '------------------------------------------------'
      DO 200 I=1,NSYMOP
         IC=LBLC6(KSYMOP,I)
         IF (C6IFC .AND. IC(1:2) .EQ. 'SM') THEN
            READ(IC,'(2X,I2,I3)') L,M
            LC6IFC = .TRUE.
         ELSE
            LC6IFC = .FALSE.
         END IF
         DO 210 J=1,I
            JC=LBLC6(KSYMOP,J)
            DO 100 IEFF=1,KZRED
               OSCSTR = D2 * (XMOM(IEFF,I)*XMOM(IEFF,J)*EIVAL(IEFF))
               DO 110 NG=0,NGRID
                  ALPHA(I,J,NG) = ALPHA(I,J,NG)
     &                   + OSCSTR/(EIVAL(IEFF)**2 + GRIDSQ(NG))
 110           CONTINUE
               DO 111 KK=1,NMOMC
                  CAUCHY(I,J,KK) = CAUCHY(I,J,KK)
     &                   + OSCSTR/(EIVAL(IEFF)**(2*KK) )
 111           CONTINUE
 100        CONTINUE
            IF (END .AND. LC6IFC) THEN
               IF (JC(1:2) .EQ. 'SM' .AND.
     &             ALPHA(I,J,0) .GT. THRESHW ) THEN
                  READ(JC,'(2X,I2,I3)') LP,MP
                  DO 112 K=0, NGRID
                     WRITE(LURSP7,'(4(1X,I2),2X,I3,2X,1P,D15.8,3X,I3)')
     &                 L,LP,M,MP,K,ALPHA(I,J,K),ITYPE
  112             CONTINUE
               ENDIF
            ENDIF
            IF (IPRRSP.GE.1) THEN
               IF (I.EQ.J) WRITE(LUPRI,'(A)') '  '
               WRITE(LUPRI,'(1X,A8,1X,A8,F20.10)') IC,JC,ALPHA(I,J,0)
            END IF
            IF(END.AND.IPRRSP.GE.5) THEN
               WRITE(LUPRI,*) '        GRIDSQ               ALPHA'
               DO 115 II=0,NGRID
                 WRITE(LUPRI,'(1X,F14.7,F20.7)')GRIDSQ(II),ALPHA(I,J,II)
 115           CONTINUE
            ENDIF
c
c  compare old and new polarizabilities and decide on convergence
c
            IF (I.EQ.J .AND. .NOT.END)  THEN
               DIFFMX = ZERO
               DO 211 II=0,NGRID
                  DIFF = ABS(OLDAL(I,II)-ALPHA(I,I,II))
                  DIFFMX = MAX(DIFF,DIFFMX)
                  IDONE=(IDONE.AND.(DIFF.LT.THCC6))
                  OLDAL(I,II) = ALPHA(I,I,II)
 211           CONTINUE
               IF (C6IFC) WRITE(LURSP7,'(A,F14.7,3X,A)')
     *               '* DIFFERENCE FOUND: ', DIFFMX, LBLC6(KSYMOP,I)
               IF (IPRRSP .GT. 1) THEN
                  WRITE(LUPRI,'(A,F14.7,3X,A/)')
     *               '   Max difference in grid polarizabilities: ',
     *               DIFFMX, LBLC6(KSYMOP,I)
               ENDIF
            ENDIF
c
 210     CONTINUE
 200  CONTINUE
      CALL GPCLOSE(LURSP7,'KEEP')
c
c Now do something very much the same to print the Cauchy moments
c
      IF ((END.AND.IPRRSP.GT.1) .OR. IPRRSP .GT. 10) THEN
       WRITE(LUPRI,'(//1X,A/,1X,A)') 'Cauchy moments (length)',
     &   '------------------------------------------------'
       DO 300 I=1,NSYMOP
         IC=LBLC6(KSYMOP,I)
         DO 310 J=1,I
            JC=LBLC6(KSYMOP,J)
            WRITE(LUPRI,'(//1X,A8,1X,A8)') IC,JC
            WRITE(LUPRI,'(/1X,A4,5X,A25)') '  k ','Moment = S(-k-2)'
            DO 320 KK=1,NMOMC
               WRITE(LUPRI,'(1X,I4,5X,1P,E30.20)') 2*KK-2,CAUCHY(I,J,KK)
 320        CONTINUE
 310     CONTINUE
 300   CONTINUE
      END IF
c
c Now compute the polarisabilities at real frequencies
c
      IF((NCFREQ.GT.0).AND. (END.OR.IPRRSP.GT.10)) THEN
         WRITE(LUPRI,'(/1X,A/1X,A/1X,A,F10.4,A)')
     &      'Polarisabilities at real frequencies (length)',
     &      '---------------------------------------------------',
     &      '(lowest pole in effective spectrum:',EIVAL(1),'au )'
         DO 400 I=1,NSYMOP
            IC=LBLC6(KSYMOP,I)
            DO 410 J=1,I
               JC=LBLC6(KSYMOP,J)
               WRITE(LUPRI,'(/1X,A8,2X,A8)') IC,JC
               DO 420 NG=1,NCFREQ
                  POLRL = ZERO
                  DO 430 IEFF=1,KZRED
                     POLRL = POLRL
     &                   + D2 * (XMOM(IEFF,I)*XMOM(IEFF,J)*EIVAL(IEFF))
     &                   /(EIVAL(IEFF)**2 - CFREQ(NG)**2)
 430              CONTINUE
                  WRITE(LUPRI,'(1X,2F20.10)') CFREQ(NG),POLRL
 420           CONTINUE
 410        CONTINUE
 400     CONTINUE
      END IF
c
      RETURN
      END
C  /* Deck efpolv */
      SUBROUTINE EFPOLV(IDONE,GRIDSQ,ALPHA,CAUCHY,NGRDMX,NMOMC,NSYMOP,
     *                  EIVAL,XMOM)
c-----------------------------------------------------------------------
c This subroutine calculates the frequency dependent polarisabilities in
c the velocity representation. For details on grid etc. see subroutine
c EFPOLL
c The Cauchy moments begin with S(0) here, instead of S(-2) as in EFPOLL
c-----------------------------------------------------------------------
c
#include <implicit.h>
c
#include <pi.h>
c
#include <priunit.h>
#include <infrsp.h>
#include "rspprp.h"
#include <infc6.h>
#include <wrkrsp.h>
#include <infpri.h>
c
      DIMENSION EIVAL(*),XMOM(MAXRM,*)
      DIMENSION GRIDSQ(0:NGRDMX),ALPHA(1:NSYMOP,1:NSYMOP,0:NGRDMX)
      DIMENSION CAUCHY(1:NSYMOP,1:NSYMOP,0:NMOMC)
      DIMENSION OMEGA(0:16)
      CHARACTER*4 IC,JC
      LOGICAL IDONE
c
      DATA OMEGA/0.0D0,
     *           0.0000011354D0,0.0000324954D0,0.0002074939D0,
     *           0.0007766098D0,0.0022314002D0,0.0055272185D0,
     *           0.0125684274D0,0.0273216537D0,0.0585616091D0,
     *           0.1273031181D0,0.2894765239D0,0.7170385627D0,
     *           2.0602366551D0,7.7110713698D0,49.2377677794D0,
     *           1409.1898779477D0/
c -- parameters
      PARAMETER ( ZERO = 0.0D0 , D1 = 1.0D0, D2 = 2.0D0 )
c
c compute squared gridpoints
c
      IF(GSLEGN) THEN
c -- use Gauss Legendre grid from old C6 module
         NGRID=16
         DO 5 I=0,NGRID
            GRIDSQ(I) = OMEGA(I)
   5     CONTINUE
      ELSE
c-- compute gridpoints according to Gauss Chebychev scheme
         GRIDSQ(0) = ZERO
         N= NGRID
         NN=2*N
         FAC=PI/(2*NN)
         DO 15 I=1,N
            X=FAC*(2*I-1)
            GRIDSQ(N-I+1)=(D1/TAN(X))**2
  15     CONTINUE
      ENDIF
c
c initialise alpha and cauchy moments to zero
c
      CALL DZERO(ALPHA ,NSYMOP*NSYMOP*(NGRID+1))
      CALL DZERO(CAUCHY,NSYMOP*NSYMOP*(NMOMC+1))
c
c Compute the polarisability and cauchy moments from the effective spectrum.
c We only need the lower triangle, hence we have j<= i in the second loop
c
      IF (IPRRSP .GE. 1) WRITE(LUPRI,'(//1X,A/,1X,A)')
     &   'S(0) and static polarisabilities (velocity)',
     &   '------------------------------------------------'
      IC = '    '
      JC = '    '
      DO 200 I=1,NSYMOP
         IF(LBLC6(KSYMOP,I)(2:5).EQ.'DIPV' ) IC=LBLC6(KSYMOP,I)(1:1)
         DO 210 J=1,I
            IF(LBLC6(KSYMOP,J)(2:5).EQ.'DIPV' ) JC=LBLC6(KSYMOP,J)(1:1)
            DO 100 IEFF=1,KZRED
               OSCSTR = D2 * (XMOM(IEFF,I)*XMOM(IEFF,J)/EIVAL(IEFF))
               DO 110 NG=0,NGRID
                  ALPHA(I,J,NG) = ALPHA(I,J,NG)
     &                   + OSCSTR/(EIVAL(IEFF)**2 + GRIDSQ(NG))
 110           CONTINUE
               CAUCHY(I,J,0) = CAUCHY(I,J,0) + OSCSTR
               DO 111 KK=1,NMOMC
                  CAUCHY(I,J,KK) = CAUCHY(I,J,KK)
     &                   + OSCSTR/(EIVAL(IEFF)**(2*KK) )
 111           CONTINUE
 100        CONTINUE
            IF (IPRRSP.GE.1) THEN
               WRITE(LUPRI,'(1X,A4,1X,A4,F20.10)') IC,JC,CAUCHY(I,J,0)
               WRITE(LUPRI,'(1X,A4,1X,A4,F20.10)') IC,JC,ALPHA(I,J,0)
            END IF
            IF(IDONE .AND. IPRRSP.GE.10) THEN
               WRITE(LUPRI,*) '        GRIDSQ               ALPHA'
               DO 115 II=0,NGRID
                 WRITE(LUPRI,'(1X,F14.7,F20.7)')GRIDSQ(II),ALPHA(I,J,II)
 115           CONTINUE
            ENDIF
 210     CONTINUE
 200  CONTINUE
c
c Now do something very much the same to print the Cauchy moments
c
      IF ((IDONE.AND.IPRRSP.GT.1) .OR. IPRRSP .GT. 10) THEN
       WRITE(LUPRI,'(//1X,A/,1X,A)') 'Cauchy moments (velocity)',
     &   '------------------------------------------------'
       IC = '    '
       JC = '    '
       DO 300 I=1,NSYMOP
         IF(LBLC6(KSYMOP,I)(2:5).EQ.'DIPV' ) IC=LBLC6(KSYMOP,I)(1:1)
         DO 310 J=1,I
            IF(LBLC6(KSYMOP,J)(2:5).EQ.'DIPV' ) JC=LBLC6(KSYMOP,J)(1:1)
            WRITE(LUPRI,'(//1X,A4,1X,A4)') IC,JC
            WRITE(LUPRI,'(/1X,A4,5X,A25)') '  k ','Moment = S(-k)'
            DO 320 KK=0,NMOMC
               WRITE(LUPRI,'(1X,I4,5X,1P,E30.20)') 2*KK,CAUCHY(I,J,KK)
 320        CONTINUE
 310     CONTINUE
 300   CONTINUE
      END IF
c
c Now compute the polarisabilities at real frequencies
c
      IF((NCFREQ.GT.0).AND. (IDONE.OR.IPRRSP.GT.10)) THEN
         WRITE(LUPRI,'(/1X,A/,1X,A/1X,A,F10.4,A)')
     &      'Polarisabilities at real frequencies (velocity)',
     &      '---------------------------------------------------',
     &      '(lowest pole in effective spectrum:',EIVAL(1),'au )'
         DO 400 I=1,NSYMOP
            IF(LBLC6(KSYMOP,I)(2:5).EQ.'DIPV' ) IC=LBLC6(KSYMOP,I)(1:1)
            DO 410 J=1,I
               IF(LBLC6(KSYMOP,J)(2:5).EQ.'DIPV')JC=LBLC6(KSYMOP,J)(1:1)
               WRITE(LUPRI,'(/1X,A4,3X,A4)') IC,JC
               DO 420 NG=1,NCFREQ
                  POLRL = ZERO
                  DO 430 IEFF=1,KZRED
                     POLRL = POLRL
     &                   + D2 * (XMOM(IEFF,I)*XMOM(IEFF,J))
     &                   /(EIVAL(IEFF)*(EIVAL(IEFF)**2 - CFREQ(NG)**2))
 430              CONTINUE
                  WRITE(LUPRI,'(1X,2F20.10)') CFREQ(NG),POLRL
 420           CONTINUE
 410        CONTINUE
 400     CONTINUE
      END IF
c
      RETURN
      END
C  /* Deck c6prps */
      SUBROUTINE C6PRPS(IDONE,OPTYPE,KSYMOP,KZRED,XMOM,EIVAL)
C
C 25-Jun-1990 PJ+HJAAJ
C
C     Print pseudo-spectrum
C
#include <implicit.h>
#include <codata.h>
      LOGICAL     IDONE
      CHARACTER*8 OPTYPE
      DIMENSION   XMOM(KZRED),EIVAL(KZRED)
      PARAMETER ( NSUM = 15 )
      DIMENSION   SUMOSC(NSUM)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
#include <priunit.h>
#include <infpri.h>
C
      IF (IDONE) THEN
         WRITE(LUPRI,'(//A)') ' Final pseudo-spectrum for'
      ELSE
         WRITE(LUPRI,'(//A)') ' Intermediate pseudo-spectrum for'
      END IF
      WRITE(LUPRI,'(/A,A/A,I4)')
     *   ' operator type:    ',OPTYPE,
     *   ' operator symmetry:',KSYMOP
      IF (OPTYPE(2:4) .NE. 'DIP') THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(//3A/A)')
     &      ' C6PRPS WARNING, oscillator strength is not',
     &      ' defined for ',OPTYPE,
     &      ' C6PRPS cannot continue.'
         RETURN
      END IF
      CALL DZERO(SUMOSC,NSUM)
      WRITE(LUPRI,'(/2A)')
     *   ' State no:   Transition moment:',
     *   '    3*Oscillator strength:   Energy (eV):'
C:State no:   Transition moment:    3*Oscillator strength:   Energy (eV):
C:    5        0.123456789012          0.123456789012      0.123456789012
C(I5,1P,G26.12,G24.12,G20.12)
      DO 500 INUM = 1,KZRED
         XTEMP=XMOM(INUM)
         ETEMP=EIVAL(INUM)
         IF (OPTYPE(2:5) .EQ. 'DIPL') THEN
            OSCSTR=D2*XTEMP*XTEMP*ETEMP
         ELSE IF (OPTYPE(2:5) .EQ. 'DIPV') THEN
            OSCSTR=D2*XTEMP*XTEMP/ETEMP
         ELSE
            OSCSTR=D0
         END IF
         WRITE(LUPRI,'(I5,1P,G26.12,G24.12,G20.12)')
     *      INUM,XTEMP,OSCSTR, ETEMP*XTEV
         DO 400 ISUM = 1,NSUM
            SUMOSC(ISUM) = SUMOSC(ISUM)
     *                   + D2*XTEMP*XTEMP*ETEMP**(-ISUM)
  400    CONTINUE
 500  CONTINUE
      IF (OPTYPE(2:5) .EQ. 'DIPL') THEN
         WRITE (LUPRI,'(/A)') ' Sum rules (length)'
         WRITE (LUPRI,'(/A)') '   k           S(k)'
         WRITE (LUPRI,'(I5,1P,G20.12)')
     *      ( (-ISUM-1), SUMOSC(ISUM), ISUM = 1,NSUM)
      ELSE IF (OPTYPE(2:5) .EQ. 'DIPV') THEN
         WRITE (LUPRI,'(/A)') ' Sum rules (velocity)'
         WRITE (LUPRI,'(/A)') '   k           S(k)'
         WRITE (LUPRI,'(I5,1P,G20.12)')
         WRITE (LUPRI,'(I5,1P,G20.12)')
     *      ( (-ISUM+1), SUMOSC(ISUM), ISUM = 1,NSUM)
      END IF
      RETURN
      END
C  /* Deck paddy */
#if defined (VAR_PADDY)
C920722: Old C6 code by Patrick Fowler
      SUBROUTINE PADDY(NCAUCH,NI,CAUCHY,COMI,WRK )
#include <implicit.h>
#include <priunit.h>
#include <infpri.h>
#include <infrsp.h>
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 )
C
      DIMENSION PLVAL(16)
      DIMENSION CAUCHY(*),COMI(*),WRK(*)
C
C DETERMINE DIMENSION OF LINEAR SET OF EQUATIONS
C
      IF (MOD(NCAUCH,2).GT.0) THEN
         K =  0
         N = (NCAUCH-1)/2
         M = N
      ELSE
         K = -1
         N = NCAUCH/2
         M = N - 1
      END IF
C
      KCAUNE = 1
      KCVEC  = KCAUNE + NCAUCH
      KCMAT  = KCVEC  + N
      KBVEC  = KCMAT  + N*N
      KAVEC  = KBVEC  + N +1
      KWRK1  = KAVEC  + M + 1
C
      DO 50 I = 1,NCAUCH
         WRK(KCAUNE-1+I) = CAUCHY(I)
 50   CONTINUE
      DO 100 I = 1,N
         WRK(KCVEC-1+I) = -WRK(KCAUNE+N+K+I)
         DO 200 J = 1,N
            WRK(KCMAT-1+(J-1)*N+I) = WRK(KCAUNE+N+K+I-J)
 200     CONTINUE
 100  CONTINUE
C
      CALL INVMAT(WRK(KCMAT),WRK(KWRK1),N,N)
C
      CALL DZERO(WRK(KBVEC),N+1)
      CALL DGEMM('N','N',N,1,N,1.D0,
     &           WRK(KCMAT),N,
     &           WRK(KCVEC),N,1.D0,
     &           WRK(KBVEC+1),N)
      WRK(KBVEC) = D1
      DO 300 I = 1,N+K+1
         WRK(KAVEC-1+I) = D0
         DO 400 J = 1,I
            WRK(KAVEC-1+I) = WRK(KAVEC-1+I) + WRK(KBVEC-1+J)
     *                      *WRK(KCAUNE+I-J)
 400     CONTINUE
 300  CONTINUE
      WRITE(LUPRI,'(/A,I2,A,I2,A,I3,A/)')
     *' [',N,',',M,'] Pade approximant for ',NCAUCH,' moments'
      IF (IPRRSP.GT.5) THEN
         WRITE(LUPRI,'(/A)')
     *   ' M+1 coefficients in P numerator'
         CALL OUTPUT(WRK(KAVEC),1,M+1,1,1,M+1,1,1,LUPRI)
         WRITE(LUPRI,'(/A)')
     *   ' N+1 coefficients in Q denominator'
         CALL OUTPUT(WRK(KBVEC),1,N+1,1,1,N+1,1,1,LUPRI)
      END IF
      DO 500 I = 1,NI
         CALL POLVAL(M+1,WRK(KAVEC),COMI(I),VALP)
         CALL POLVAL(N+1,WRK(KBVEC),COMI(I),VALQ)
         PLVAL(I) = VALP/VALQ
         WRITE(LUPRI,'(I5,A,1P,G16.8,A,G16.8)')
     *  I, ' alpha at imaginary frequency',SQRT(COMI(I)),' = '
     *  ,PLVAL(I)
 500  CONTINUE
      CALL C6INT(PLVAL,PLVAL,C6TEMP)
      WRITE(LUPRI,'(/A,I2,A,I2,A,G16.8) ')
     *' C6 coefficient for [',N,',',M,'] Pade approximant = ',
     *C6TEMP
      RETURN
      END
      SUBROUTINE POLVAL(NUM,POLCOF,FREQ,VALUE)
#include <implicit.h>
      DIMENSION POLCOF(NUM)
      VALUE = 0.0D0
      DO 100 I = 1,NUM
         VALUE = VALUE + POLCOF(I)*((-FREQ)**(I-1))
 100  CONTINUE
      RETURN
      END
      SUBROUTINE C6INT(A1,A2,SUM)
#include <implicit.h>
C
#include <pi.h>
C
      DIMENSION RLOW(8),WLOW(8)
      DIMENSION W(16),A1(16),A2(16)
      DATA RLOW/
     &0.97891421016235D+00,0.89222197421380D+00,
     &0.74931737854740D+00,0.57063582016217D+00,0.38177105339712D+00,
     &0.20977936861551D+00,0.79300559811486D-01,0.90273770256471D-02/
      DATA WLOW/
     &0.27152459411755D-01,0.62253523938648D-01,
     &0.95158511682493D-01,0.12462897125553D+00,0.14959598881658D+00,
     &0.16915651939500D+00,0.18260341504492D+00,0.18945061045507D+00/
C
C     THIS PROGRAM CALCULATES THE C6 COEFFICIENT FOR THE INTERACTION
C     OF SYSTEM 1 WITH SYSTEM 2, GIVEN THE POLARISABILITY OF EACH
C     SYSTEM AT A SET OF 16 PRE-DEFINED IMAGINARY FREQUENCIES:
C
C     C6 = (3/PI) INTEGRAL(0,INFINITY) ALPHA1(IW) ALPHA2(IW) DW
C
      OMEGA0=0.2D0
C
C     FOR EACH PARTNER READ THE POLARISABILITY AT THE 16 FREQUENCIES.
C     THE FREQUENCIES REQUIRED ARE LISTED AS THE NEGATIVES OF THEIR
C     SQUARES BELOW.
C
C 0.0000011354
C 0.0000324954
C 0.0002074939
C 0.0007766098
C 0.0022314002
C 0.0055272185
C 0.0125684274
C 0.0273216537
C 0.0585616091
C 0.1273031181
C 0.2894765239
C 0.7170385627
C 2.0602366551
C 7.7110713698
C 49.2377677794
C 1409.1898779475
C
C     FOR A DESCRIPTION OF THE QUADRATURE SCHEME, SEE AMOS ET AL.
C     J. PHYS. CHEM. 1985 VOL. 89 PAGE 2186.
C
      DO 1 I=1,8
         WW=WLOW(I)*OMEGA0/PI
         RR=SQRT(RLOW(I))
         I2=16+1-I
         W(I)=WW/(1.0D0+RR)**2
         W(I2)=WW/(1.0D0-RR)**2
1     CONTINUE
      SUM=0.0D0
      DO 2 I=1,16
2     SUM=SUM+A1(I)*A2(I)*W(I)
      SUM=6.0D0*SUM
      RETURN
      END
#endif
C  /* Deck invmat */
#if defined (VAR_PADDY)
C920722: Old C6 code by Patrick Fowler
C        These two routines are used to invert a matrix in PADDY
      SUBROUTINE INVMAT(A,B,MATDIM,NDIM)
C FIND INVERSE OF MATRIX A
C INPUT
C        A : MATRIX TO BE INVERTED
C        B : SCRATCH ARRAY
C        MATDIM : PHYSICAL DIMENSION OF MATRICES
C        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
C
C OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
C WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION A(MATDIM,MATDIM),B(MATDIM,MATDIM)
C
      DETERM=0.0D0
      EPSIL=0.0D0
      ITEST=0
      CALL BNDINV(A,B,NDIM,DETERM,EPSIL,ITEST,MATDIM)
C
      IF( ITEST .NE. 0 ) THEN
        WRITE (LUPRI,'(A,I3)') ' INVERSION PROBLEM NUMBER..',ITEST
      END IF
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
        WRITE(LUPRI,*) ' INVERTED MATRIX '
        CALL WRTMAT(A,NDIM,NDIM,MATDIM,MATDIM)
      END IF
C
      RETURN
      END
        SUBROUTINE BNDINV(A,EL,N,DETERM,EPSIL,ITEST,NSIZE)
C
C       DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
C       FROM "DLYTAP".
C
C*      DOUBLE PRECISION E,F
C*      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
#include <implicit.h>
        DIMENSION A(NSIZE,1),EL(NSIZE,1)
        IF(N.LT.2)GO TO 140
        ISL2=0
        K000FX=2
        IF(ISL2.EQ.0)INDSNL=2
        IF(ISL2.EQ.1)INDSNL=1
C       CALL SLITET(2,INDSNL)
C       CALL OVERFL(K000FX)
C       CALL DVCHK(K000FX)
C
C       SET EL = IDENTITY MATRIX
        DO 30 I=1,N
        DO 10 J=1,N
 10     EL(I,J)=0.0D0
 30     EL(I,I)=1.0D0
C
C       TRIANGULARIZE A, FORM EL
C
        N1=N-1
        M=2
        DO 50 J=1,N1
        DO 45 I=M,N
        IF(A(I,J).EQ.0.0D0)GO TO 45
        D=SQRT(A(J,J)*A(J,J)+A(I,J)*A(I,J))
        C=A(J,J)/D
        S=A(I,J)/D
 38     DO 39 K=J,N
        D=C*A(J,K)+S*A(I,K)
        A(I,K)=C*A(I,K)-S*A(J,K)
        A(J,K)=D
 39     CONTINUE
        DO 40 K=1,N
        D=C*EL(J,K)+S*EL(I,K)
        EL(I,K)=C*EL(I,K)-S*EL(J,K)
        EL(J,K)=D
 40     CONTINUE
 45     CONTINUE
 50     M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,51),K000FX
C
C       CALCULATE THE DETERMINANT
 51     DETERP=A(1,1)
        DO 52 I=2,N
 52     DETERP=DETERP*A(I,I)
        DETERM=DETERP
C       CALL OVERFL(K000FX)
C       GO TO (140,520,520),K000FX
C
C       IS MATRIX SINGULAR
 520    F=A(1,1)
        E=A(1,1)
        DO 58 I=2,N
        IF(ABS(F).LT.ABS(A(I,I)))F=A(I,I)
        IF(ABS(E).GT.ABS(A(I,I)))E=A(I,I)
 58     CONTINUE
        EPSILP=EPSIL
        IF(EPSILP.LE.0)EPSILP=1.0E-8
        RAT=E/F
        IF(ABS(RAT).LT.EPSILP)GO TO 130
C
C       INVERT TRIANGULAR MATRIX
        J=N
        DO 100 J1=1,N
C       CALL SLITE(2)
        I=J
        ISL2=1
        DO 90 I1=1,J
C       CALL SLITET(2,K000FX)
        IF(ISL2.EQ.0)K000FX=2
        IF(ISL2.EQ.1)K000FX=1
        IF(ISL2.EQ.1)ISL2=0
        GO TO (70,75),K000FX
 70     A(I,J)=1.0D0/A(I,I)
        GO TO 90
 75     KS=I+1
        D=0.0D0
        DO 80 K=KS,J
 80     D=D+A(I,K)*A(K,J)
        A(I,J)=-D/A(I,I)
 90     I=I-1
 100    J=J-1
C       CALL OVERFL(K000FX)
C       GO TO (140,103,103),K000FX

C103    CALL DVCHK(K000FX)
C       GO TO (140,105),K000FX
C
C       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
 105    M=1
        DO 120 I=1,N
        DO 118 J=1,N
        D=0.0D0
        DO 107 K=M,N
 107    D=D+A(I,K)*EL(K,J)
        EL(I,J)=D
 118    CONTINUE
 120    M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,123,123),K000FX
C
C       RECOPY EL TO A
 123    DO 124 I=1,N
        DO 124 J=1,N
 124    A(I,J)=EL(I,J)
        ITEST=0
C126    IF(INDSNL.EQ.1)CALL SLITE(2)
 126    IF(INDSNL.EQ.1)ISL2=1
        RETURN
C
 130    ITEST=1
        GO TO 126
 140    ITEST=-1
        GO TO 126
        END
#endif
