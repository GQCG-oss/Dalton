C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
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
#include <single.h>

#undef SRDFT_DEBUG

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRDFT1JT(ND_SIM,EXCMAT,EDFTY,DMAT,DOLND,CORX,CORY,CORZ,
     &                  WEIGHT,NBUF,GSO,NCNT,DMAGAOSR,DMAGAOT,DOGGA,
     &                  DOERG,DOGRD,DOATR,WORK,LWORK,IPRINT)
C*****************************************************************************
C     Calculate DFT energy and potential contributions for
C     SR-DFT hybrid methods
C
C     Based on the routine SRDFT1 of J. K. Pedersen
C
C     Created : 15-02-05, J. Toulouse
C     Modified: 17-08-09, J. Toulouse, add spin density
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pi.h>
#include <dummy.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpar.h>
#include <nuclei.h>
#include <dftinf.h>
#include <dftcom.h>
#include <gnrinf.h>
#include <chrnos.h>

      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)

      LOGICAL   DOLND,DOGGA,DOERG,DOGRD,DOATR,FROMVX,DO_SPINDNS

      DIMENSION EXCMAT(NBAST,NBAST,*),DMAT(N2BASX,*),CORX(NBUF),
     &          CORY(NBUF),CORZ(NBUF),WEIGHT(NBUF),GSO(NBAST*NTYPSO),
     &          NCNT(NBAST),DMAGAOSR(NBAST),DMAGAOT(NBAST),WORK(LWORK),
     &          RHG(3),COR(3)

      CHARACTER*16 QUADNAME
#include <sdpre.h>

#ifdef SRDFT_DEBUG
!     Debug
      WRITE(LUPRI,*) 'SRDFT1JT: entering, ND_SIM = ',ND_SIM
      WRITE(LUPRI,*) 'SRDFT1JT: DOATR,DOGRD,DOGGA  =',DOATR,DOGRD,DOGGA
      WRITE(LUPRI,*) 'SRDFT1JT: ERFEXP =',ERFEXP
      WRITE(LUPRI,*) 'SRDFT1JT: LUPRI  =',LUPRI, IPRINT
      WRITE(LUPRI,*) 'SRDFT1JT: HFXFAC =',HFXFAC
      WRITE(LUPRI,*) 'SRDFT1JT: DOSRX_PBERI=',DOSRX_PBERI, DOERG
C
#endif
C     Initializations -------------------------------------------------------

      DO_SPINDNS = ISPIN .NE. 1  !JT Aug 09

!JT      WRITE(LUPRI,*) 'DO_SPINDNS=',DO_SPINDNS
!JT      WRITE(LUPRI,*) 'DMAT(1,1)=',DMAT(1,1)
!JT      WRITE(LUPRI,*) 'DMAT(1,3)=',DMAT(1,3)

      IF (ND_SIM .NE. 1) THEN
         WRITE(LUPRI,*) 'SRDFT1: illegal ND_SIM ',ND_SIM
         CALL QUIT('SRDFT1: illegal ND_SIM value')
      END IF

      IF (DOGRD) THEN
         CALL QUIT('SRDFTJT: DOGRD not implemented yet')
      END IF

!     (Integrated) exchange and correlation energies
      ExInt  = D0
      EcInt  = D0
      IF (DOERG) EDFTY  = D0

C     DFT grid (NTOTAL is the total number of grid points)
      IF (.NOT.GRDONE) THEN
         CALL DFTGRD(WORK,LWORK,NTOTAL,1,.FALSE.)
         CALL CONDFT
         GRDONE = .TRUE.
      END IF

C     Electron number
      ELCTRN = D0

C     Integrated spin density
      SPINDEN = D0

C     Thresholds for screening
C     THRELCTRN makes sure that error in no. of electrons
C     from grid integration less than DFTHR
      THRELCTOT = 1.0D-7
      THRELCTRN = D0
      THRESRHO  = 1.0D-12
      THRDRHO   = D0
      DFTHR0    = D0
      DFTHRI    = D0
      DFTHRL    = D0

C     Skipped grid points and nb of electrons
      NSKIP1 = 0     
      ELCTRSKIP = D0 

C     For quadrature integration
C
C     Make quadname: Can take 9999 procs
C
      QUADNAME = 'DALTON.QUAD.'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
C
      LUQUAD = -1
      CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)
C
      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
      kpoint = npoint/10

      IF (NPOINT.GT.0) THEN
C      jpnt = 0  !debug
      NPNTS = NPNTS + NPOINT
      CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)

C   Loop over grid points -----------------------------------------------
      DO 300 IPNT = 1, NPOINT

C     Print grid information
      IF (IPRINT .GT. 100) THEN
        WRITE (LUPRI,'(2X,I6,4F12.6)') 
     &        IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
      END IF

C     Quadrature weight
      WGHT  = WEIGHT(IPNT)
      WGHTX = (1.0D0-HFXFAC)*WGHT
      WGHTC = WGHT

C     Get AOs
      THRINT = DFTHRI/WGHT
      COR(1) = CORX(IPNT)
      COR(2) = CORY(IPNT)
      COR(3) = CORZ(IPNT)
      CALL GETSOS(GSO,NCNT,COR,WORK,LWORK,NBAST,DOLND,DOGGA,
     &            THRINT,IPRINT)

C     Density
      CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHO,RHO13,DMAGAOSR,
     &                THRINT,IPRINT)
      CALL DCOPY(NBAST,DMAGAOSR,1,DMAGAOT,1)

C     Contribution to electron number (N_alpha + N_beta)
      DELCTRN = WGHT*RHO

C     Spin density
      IF (DO_SPINDNS) THEN
C        DMAT(1,3) is the spin density matrix
         CALL GETRHO_OLD(DMAT(1,3),GSO(KSO0),RHOS,RHOS13,DMAGAOSR,
     &                   THRINT,IPRINT)
         
C        Contribution to spin density (N_alpha - N_beta)
         SPINDEN = SPINDEN + WGHT*RHOS
      ELSE
         RHOS = 0.0D0
      END IF

C     Screening
!      print *,'DELCTRN=',DELCTRN
!      print *,'THRELCTRN=',THRELCTRN
      IF (ABS(DELCTRN).LE.THRELCTRN .OR. RHO.LE.THRESRHO) THEN
         NSKIP1 = NSKIP1 + 1
         ELCTRSKIP = ELCTRSKIP + DELCTRN
         THRELCTRN = (THRELCTOT-ELCTRSKIP)/100
       ELSE

C     Gradient of density
      IF (DOGGA) THEN
        CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAOSR,1,D0,RHG,1)
        RHOGRD = SQRT(RHG(1)**2 + RHG(2)**2 + RHG(3)**2)
      ENDIF

C     Hessian of density (for second-order optimization)
      IF (DOATR) THEN
       CALL DFTRHH(DMAT,DMAGAOSR,GSO(KSO0),GSO(KSO1),
     >            GSO(KSO2),RHG,RHOLAP,RHOGHG)
      ENDIF

C     Number of electrons
      ELCTRN = ELCTRN + DELCTRN


      IF (DOERG) THEN

C     Exchange and correlation energies and derivatives
      NORDER = 1  ! first-order derivatives
      CALL SRDFTEXCJT(RHO,RHOS,RHOGRD,CHIVAL,NORDER,ERFEXP,
     >           Ex,d1Exdrhot,d1Exdrhos,d1Exdgrd,d1Exdgrda2,
     >           d2Exdrhotdrhot,d2Exdrhotdrhos,d2Exdrhosdrhos,
     >           d2Exdgrddgrd,d2Exdrhodgrd,
     >           Ec,d1Ecdrhot,d1Ecdrhos,d1Ecdgrd,d1Ecdgrda2,
     >           d2Ecdrhotdrhot,d2Ecdrhotdrhos,d2Ecdrhosdrhos,
     >           d2Ecdgrddgrd,d2Ecdrhodgrd)


C      Summation with quadrature weigths
       ExInt     = ExInt + WGHTX*Ex
       EcInt     = EcInt + WGHTC*Ec

C      Potentials with quadrature weigths
       d1Exdrhot  = WGHTX*d1Exdrhot
       d1Exdrhos  = WGHTX*d1Exdrhos
       d1Ecdrhot  = WGHTC*d1Ecdrhot
       d1Ecdrhos  = WGHTC*d1Ecdrhos
       d1Exdgrda2  = WGHTX*d1Exdgrda2
       d1Ecdgrda2  = WGHTC*d1Ecdgrda2

C      Sum of exchange and correlation
       d1Excdrhot = d1Exdrhot + d1Ecdrhot
       d1Excdrhos = d1Exdrhos + d1Ecdrhos
       d1Excdgrda2 = d1Exdgrda2 + d1Ecdgrda2

C      Exchange-correlation contribution to Kohn-Sham matrix (charge and spin components)
       FROMVX = .FALSE.
       CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHG,d1Excdrhot,
     &             d1Excdgrda2,DOGGA,FROMVX,DFTHRL)
       IF (DO_SPINDNS) THEN
         CALL DFTKSM(EXCMAT(1,1,2),GSO(KSO0),GSO(KSO1),RHG,d1Excdrhos,
     &             d1Excdgrda2,DOGGA,FROMVX,DFTHRL)
       ENDIF

C
       ENDIF  !DOERG

C      Hessian transformation for second-order terms
       IF (DOATR) THEN
C        matrix 1 is DTAO
C        matrix 2 is DXTAO
         JDXAO = 2
         CALL SRDFTLTRJT(JWOPSY,DMAT(1,JDXAO),EXCMAT,
     &                 WGHTX,WGHTC, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                 ERFEXP,RHO,RHOS,
     &                 RHOGRD,RHOLAP,RHOGHG,RHG,DMAGAOT,DOGGA)
       ENDIF

       ENDIF  !ABS(DELCTRN).LE.THRELCTRN
C
  300  CONTINUE
C End of loop over points ------------------------------------------------
C
       GO TO 200
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
      END IF  !NPOINT.GT.0

C Termination --------------------------------------------------------------

      CALL GPCLOSE(LUQUAD,'KEEP')

      IF (DOERG) THEN
         EDFTY = ExInt + EcInt
         FAC = DP5
      ELSE IF (DOATR) THEN
         FAC = D1
      ELSE
         CALL QUIT('both DOERG and DOATR are false')
      END IF

      IF (DO_SPINDNS) THEN
        NXCMAT = 2
      ELSE
        NXCMAT = 1
      ENDIF

      DO K = 1, NXCMAT
        DO I = 1, NBAST
          DO J = 1, I - 1
            AVERAG = FAC*(EXCMAT(I,J,K) + EXCMAT(J,I,K))
            EXCMAT(I,J,K) = AVERAG
            EXCMAT(J,I,K) = AVERAG
          END DO
        END DO
      END DO

C     Print section
!     IF (IPRINT .GE. 3) THEN
      IF (DOERG) THEN
         ELCTOT = ELCTRN + ELCTRSKIP
         WRITE (LUPRI,'(4(/1X,A,F18.10))')
     &'No. of non-skipped electrons from num. int.  :',ELCTRN,
     &'No. of     skipped electrons   in num. int.  :',ELCTRSKIP,
     &'Sum of non-skipped and skipped electrons     :',ELCTOT,
     &'Integrated spin density (N_alpha-N_beta)     :',SPINDEN

         WRITE (LUPRI,'(/1X,A,F18.10)')
     &'Exchange energy                                     :',ExInt
         WRITE (LUPRI,'(1X,A,F18.10)')
     &'Correlation energy                                  :',EcInt
         WRITE (LUPRI,'(1X,A,F18.10)')
     &'DFT exchange-correlation energy                     :',EDFTY

         IF (IPRINT.GE.5) THEN
            WRITE(LUPRI,'(//A//)') ' <<< Kohn-Sham matrix >>>'
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF

      END IF
!     END IF  ! IPRINT .GE. 3

      RETURN 
      END

C****************************************************************************
      SUBROUTINE SRDFTEXCJT(rhot,rhos,RHOGRD,CHIVAL,NORDER,ERFEXP,
     >           Ex,d1Exdrhot,d1Exdrhos,d1Exdgrd,d1Exdgrda2,
     >           d2Exdrhotdrhot,d2Exdrhotdrhos,d2Exdrhosdrhos,
     >           d2Exdgrddgrd,d2Exdrhodgrd,
     >           Ec,d1Ecdrhot,d1Ecdrhos,d1Ecdgrd,d1Ecdgrda2,
     >           d2Ecdrhotdrhot,d2Ecdrhotdrhos,d2Ecdrhosdrhos,
     >           d2Ecdgrddgrd,d2Ecdrhodgrd)
C*****************************************************************************
C    Driver for short-range DFT exchange and correlation
C    energies and derivatives
C
C    Based on the routine SRDFTEXC of J. K. Pedersen
C
C    Created : 15-02-05, J. Toulouse
C    Modified: 17-08-09, J. Toulouse, add spin density
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <inforb.h>
#include <infinp.h>
#include <dftinf.h>
#include <dftcom.h>

      double precision rhot,rhos,RHOGRD,CHIVAL
      integer NORDER
      logical ERFEXP, DO_SPINDNS
      double precision Ex,d1Exdrhot,d1Exdrhos,d1Exdgrd,d1Exdgrda2
      double precision d2Exdrhotdrhot,d2Exdrhotdrhos,d2Exdrhosdrhos
      double precision d2Exdgrddgrd,d2Exdrhodgrd
      double precision Ec,d1Ecdrhot,d1Ecdrhos,d1Ecdgrd,d1Ecdgrda2
      double precision d2Ecdrhotdrhot,d2Ecdrhotdrhos,d2Ecdrhosdrhos
      double precision d2Ecdgrddgrd,d2Ecdrhodgrd

C     functionals
      external ESRX_PBEHSEERF,  ESRX_PBETCSERF, ESRX_PBERIERF
      external ESRX_PBEGWSERF
      external ESRC_PBETCSJERF, ESRC_PBETCSERF, ESRC_PBERIERF
      external ESRC_PBEWIERF,   ESRC_PBEGWSERF
C     spin functionals
      external ESRX_LDAERFSPIN, ESRC_LDAERFSPIN

      DO_SPINDNS = ISPIN .NE. 1  !JT Aug 09

c  Exchange  ------------------------------------------------------------

!     Initialisation
      Ex     = 0.d0
      d1Exdrhot  = 0.d0
      d1Exdrhos  = 0.d0
      d1Exdgrd  = 0.d0
      d1Exdgrda2  = 0.d0
      d2Exdrhotdrhot = 0.d0
      d2Exdrhotdrhos = 0.d0
      d2Exdrhosdrhos = 0.d0
      d2Exdgrddgrd = 0.d0
      d2Exdrhodgrd = 0.d0

!     short-range exchange PBE of HSE
      IF(DOSRX_PBEHSE) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRXPBEHSE not implemented for erfgau')
       ELSE
        call ESRX_PBEHSEERF(rhot,RHOGRD,CHIVAL,Ex)
        call DESR(ESRX_PBEHSEERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
       ENDIF
      ENDIF

!     short-range exchange PBE of TCS
      IF(DOSRX_PBETCS) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRXPBETCS not implemented for erfgau')
       ELSE
        call ESRX_PBETCSERF(rhot,RHOGRD,CHIVAL,Ex)
        call DESR(ESRX_PBETCSERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
       ENDIF
      ENDIF

!     short-range exchange PBE of GWS ! JT 11-08-09
      IF(DOSRX_PBEGWS) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRXPBEGWS not implemented for erfgau')
       ELSE
        call ESRX_PBEGWSERF(rhot,RHOGRD,CHIVAL,Ex)
        call DESR(ESRX_PBEGWSERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
       ENDIF
      ENDIF

C  Manu 01-02-2006 begin
!     short-range exchange RI
      IF(DOSRX_PBERI) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRXRI not implemented for erfgau')
       ELSE
        call ESRX_PBERIERF(rhot,RHOGRD,CHIVAL,Ex)
        call DESR(ESRX_PBERIERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
       ENDIF
      ENDIF
C  Manu 01-02-2006 end

!     short-range spin-dependent exchange LDA ! JT 11-08-09
      IF(DOSRX_LDAS) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRXLDAS not implemented for erfgau')
       ELSE
        call ESRX_LDAERFSPIN(rhot,rhos,CHIVAL,DO_SPINDNS,Ex)
        call DESRSPIN(ESRX_LDAERFSPIN,rhot,rhos,CHIVAL,DO_SPINDNS,
     >             NORDER,
     >             d1Exdrhot,d1Exdrhos,
     >             d2Exdrhotdrhot,d2Exdrhotdrhos,d2Exdrhosdrhos)
       ENDIF
      ENDIF
 
C  Correlation  ----------------------------------------------------------

C     Initialisations
      Ec           = 0.d0
      d1Ecdrhot    = 0.d0
      d1Ecdrhos    = 0.d0
      d1Ecdgrd     = 0.d0
      d1Ecdgrda2   = 0.d0
      d2Ecdrhotdrhot = 0.d0
      d2Ecdrhotdrhos = 0.d0
      d2Ecdrhosdrhos = 0.d0
      d2Ecdgrddgrd = 0.d0
      d2Ecdrhodgrd = 0.d0

!     short-range correlation PBE of TCS
      IF(DOSRC_PBETCS) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRCPBETCS not implemented for erfgau')
       ELSE
!        call VSRCPBEERF(ec,d1ecdrho,d1ecdgrd,d1ecdgrda2,d2ecdrhodrho,
!     >                  d2ecdgrddgrd,d2ecdrhodgrd,rho,rhogrd,chival,
!     >                  norder)
        call ESRC_PBETCSERF(rhot,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBETCSERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
       ENDIF
      ENDIF

!     short-range correlation PBE of TCSJ
      IF(DOSRC_PBETCSJ) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRCPBETCSJ not implemented for erfgau')
       ELSE
        call ESRC_PBETCSJERF(rhot,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBETCSJERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
       ENDIF
      ENDIF

!     short-range correlation PBE of GWS
      IF(DOSRC_PBEGWS) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRCPBEGWS not implemented for erfgau')
       ELSE
        call ESRC_PBEGWSERF(rhot,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBEGWSERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
       ENDIF
      ENDIF

!     rational interpolation
      IF(DOSRC_PBERI) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRCRI not implemented for erfgau')
       ELSE
        call ESRC_PBERIERF(rhot,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBERIERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
       ENDIF
      ENDIF

!     weighted interpolation
      IF(DOSRC_PBEWI) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRCWI not implemented for erfgau')
       ELSE
        call ESRC_PBEWIERF(rhot,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBEWIERF,rhot,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
       ENDIF
      ENDIF

!     short-range spin-dependent correlation LDA ! JT 18-08-09
      IF(DOSRC_LDAS) THEN
       IF(ERFEXP) THEN
        CALL QUIT('SRCLDAS not implemented for erfgau')
       ELSE
        call ESRC_LDAERFSPIN(rhot,rhos,CHIVAL,DO_SPINDNS,Ec)
        call DESRSPIN(ESRC_LDAERFSPIN,rhot,rhos,CHIVAL,DO_SPINDNS,
     >             NORDER,
     >             d1Ecdrhot,d1Ecdrhos,
     >             d2Ecdrhotdrhot,d2Ecdrhotdrhos,d2Ecdrhosdrhos)
       ENDIF
      ENDIF

      RETURN
      END

C**************************************************************************
      SUBROUTINE SRDFTLTRJT(KSYMOP,DTRMAT,EXCMAT,WGHTX,WGHTC,GAO,GAO1,
     &                    CHIVAL,ERFEXP,RHO,RHOS,RHOGRD,RHOLAP,RHOGHG,
     &                    RHG,DTGAO,DOGGA)
C**************************************************************************
C     Generate the 1-index transformed DFT Hessian
C     needed in second-order optimization for short-range DFT hybrids.
C
C     Based on the routine SRDFTLT of J. K. Pedersen
C
C     Created : 15-02-05, J. Toulouse
C     Modified: 17-08-09, J. Toulouse, add spin density
C***************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pi.h>
#include <maxorb.h>
#include <inforb.h>
#include <infinp.h>
#include <nuclei.h>
#include <dftinf.h>
#include <dftcom.h>

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0,
     >           D13=1.0D0/3.0D0)
      LOGICAL DOGGA,ERFEXP

      DIMENSION DTRMAT(NBAST,NBAST),GAO(NBAST),GAO1(NBAST,3),
     &          EXCMAT(NBAST,NBAST),DTGAO(NBAST),RHG(3),
     &          DTGTAO(NBAST),B3(3)
C

      CALL DGEMV('N',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
      B0 = DDOT(NBAST,DTGAO,1,GAO,1)

      IF (DOGGA) THEN
C        B3 = GAO1'*DTGAO
         CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D0,B3,1)
C        DTGAO= DTRMAT'*GAO
         CALL DGEMV('T',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
C        B3 = B3 + GAO1'*DTGAO
         CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D1,B3,1)
         BMAX = MAX(DABS(B0),DABS(B3(1)),DABS(B3(2)),DABS(B3(3)))
      ELSE
!        check this case
         BMAX = DABS(B0)
         CALL DZERO(B3,3)
      END IF
C
      IF (BMAX.GT.DFTHRL) THEN 
C

C Exchange and correlation derivatives ------------------------------------

      NORDER = 2  ! first-order and second-order derivatives
      CALL SRDFTEXCJT(RHO,RHOS,RHOGRD,CHIVAL,NORDER,ERFEXP,
     >           Ex,d1Exdrhot,d1Exdrhos,d1Exdgrd,d1Exdgrda2,
     >           d2Exdrhotdrhot,d2Exdrhotdrhos,d2Exdrhosdrhos,
     >           d2Exdgrddgrd,d2Exdrhodgrd,
     >           Ec,d1Ecdrhot,d1Ecdrhos,d1Ecdgrd,d1Ecdgrda2,
     >           d2Ecdrhotdrhot,d2Ecdrhotdrhos,d2Ecdrhosdrhos,
     >           d2Ecdgrddgrd,d2Ecdrhodgrd)


C Total contribution -------------------------------------------------------

         ZNV = D1/RHOGRD
         FZ0 = ZNV*(WGHTX*d1Exdgrd + WGHTC*d1Ecdgrd)
         FRR = WGHTX*d2Exdrhotdrhot + WGHTC*d2Ecdrhotdrhot
         FRZ = WGHTX*d2Exdrhodgrd + WGHTC*d2Ecdrhodgrd
         FZZ = WGHTX*d2Exdgrddgrd + WGHTC*d2Ecdgrddgrd
C
         RX = ZNV*RHG(1)
         RY = ZNV*RHG(2)
         RZ = ZNV*RHG(3)

C Linear transformation -----------------------------------------------------

         IF (DOHFSRDFT) THEN
            IF (DOGGA) THEN
               BR = B3(1)*RX + B3(2)*RY + B3(3)*RZ
               FAC0 = FRR*B0 + FRZ*BR
               FACR = FRZ*B0 + FZZ*BR
            ELSE
               FAC0 = FRR*B0
            END IF
            IF (NSYM.EQ.1) THEN
               DO I = 1, NBAST
                  G0 = GAO(I)
                  IF (DOGGA) THEN
                     GX = GAO1(I,1)
                     GY = GAO1(I,2)
                     GZ = GAO1(I,3)
                  END IF
                  DO J = 1, I 
                     EXCMAT(J,I) = EXCMAT(J,I)+FAC0*G0*GAO(J)
                     IF (DOGGA) THEN
                        AX = GX*GAO(J) + G0*GAO1(J,1)
                        AY = GY*GAO(J) + G0*GAO1(J,2)
                        AZ = GZ*GAO(J) + G0*GAO1(J,3)
                        AR = AX*RX + AY*RY + AZ*RZ
                        AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) - AR*BR
                        EXCMAT(J,I) = EXCMAT(J,I)+FACR*AR+FZ0*AB
                     END IF
                  END DO
               END DO
            ELSE
               DO ISYM = 1, NSYM
                  ISTR = IBAS(ISYM) + 1
                  IEND = IBAS(ISYM) + NBAS(ISYM)
                  JSYM = MULD2H(ISYM,KSYMOP)
                  IF (ISYM.GE.JSYM) THEN
                     JSTR = IBAS(JSYM) + 1
                     JEND = IBAS(JSYM) + NBAS(JSYM)
                     DO I = ISTR, IEND
                        G0 = GAO(I)
                        IF (DOGGA) THEN
                           GX = GAO1(I,1)
                           GY = GAO1(I,2)
                           GZ = GAO1(I,3)
                        END IF
                        DO J = JSTR, MIN(I,JEND) 
                           EXCMAT(J,I) = EXCMAT(J,I) + FAC0*G0*GAO(J)
                           IF (DOGGA) THEN
                              AX = GX*GAO(J) + G0*GAO1(J,1)
                              AY = GY*GAO(J) + G0*GAO1(J,2)
                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
                              AR = AX*RX + AY*RY + AZ*RZ
                              AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) -
     &                             AR*BR
                              EXCMAT(J,I) = EXCMAT(J,I)
     &                                     + FACR*AR + FZ0*AB
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END IF
         END IF
      END IF
c
      RETURN
      END
