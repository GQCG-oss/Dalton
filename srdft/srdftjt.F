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
     &                  WEIGHT,NBUF,GSO,NCNT,DMAGAOSR,DOGGA,
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

      LOGICAL   DOLND,DOGGA,DOERG,DOGRD,DOATR,FROMVX

      DIMENSION EXCMAT(NBAST,NBAST,*),DMAT(N2BASX,*),
     &          CORX(NBUF),CORY(NBUF),CORZ(NBUF),WEIGHT(NBUF),
     &          GSO(NBAST*NTYPSO),NCNT(NBAST),
     &          DMAGAOSR(NBAST,2),
     &          COR(3),WORK(LWORK)
C
      DIMENSION RHO(4), RHG(3,4)
C     ... RHO(1) = RHOT  = rho_a + rho_b
C     ... RHO(2) = RHOS  = rho_a - rho_b
C     ... RHO(3) = RHO2A = 2 rho_a
C     ... RHO(4) = RHO2B = 2 rho_b
C     ... RHG(1:3,1) = RHGT(1:3)  = grad ( rho_a + rho_b ) = grad ( rho_t )
C     ... RHG(1:3,2) = RHGS(1:3)  = grad ( rho_a - rho_b ) = grad ( rho_s )
C     ... RHG(1:3,3) = RHG2A(1:3) = 2 grad ( rho_a )
C     ... RHG(1:3,4) = RHG2B(1:3) = 2 grad ( rho_b )

      CHARACTER*16 QUADNAME
#include <sdpre.h>

#define SRDFT_DEBUG
#ifdef SRDFT_DEBUG
!     Debug
      WRITE(LUPRI,*) 'SRDFT1JT: entering, ND_SIM = ',ND_SIM
      WRITE(LUPRI,*) 'SRDFT1JT: DOATR,DOGRD,DOGGA  =',DOATR,DOGRD,DOGGA
      WRITE(LUPRI,*) 'SRDFT1JT: ERFEXP =',ERFEXP
      WRITE(LUPRI,*) 'SRDFT1JT: LUPRI  =',LUPRI, IPRINT
      WRITE(LUPRI,*) 'SRDFT1JT: HFXFAC =',HFXFAC
      WRITE(LUPRI,*) 'SRDFT1JT: DOSRX_PBERI=',DOSRX_PBERI, DOERG
      WRITE(LUPRI,*) 'DFT_SPINDNS=',DFT_SPINDNS
      WRITE(LUPRI,*) 'D_TOT_AO matrix:'
      call output(DMAT(1,1),1,nbast,1,nbast,nbast,nbast,-1,lupri)
!     IF (DFT_SPINDNS) THEN
      WRITE(LUPRI,*) 'D_SPIN_AO matrix:'
      call output(DMAT(1,3),1,nbast,1,nbast,nbast,nbast,-1,lupri)
!     END IF
C
#endif
C     Initializations -------------------------------------------------------


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
      CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHO(1),RHOT13,DMAGAOSR(1,1),
     &                THRINT,IPRINT)

C     Contribution to electron number (N_alpha + N_beta)
      DELCTRN = WGHT*RHO(1)

C     Screening
!      print *,'DELCTRN=',DELCTRN
!      print *,'THRELCTRN=',THRELCTRN
      IF (ABS(DELCTRN).LE.THRELCTRN .OR. RHO(1).LE.THRESRHO) THEN
         NSKIP1 = NSKIP1 + 1
         ELCTRSKIP = ELCTRSKIP + DELCTRN
         THRELCTRN = (THRELCTOT-ELCTRSKIP)/100
         GO TO 300
C        ... skip this point, go to next grid point
      END IF

C     Spin density
      IF (DFT_SPINDNS) THEN
C        DMAT(1,3) is the spin density matrix
         
         CALL GETRHO_OLD(DMAT(1,3),GSO(KSO0),RHO(2),RHOS13,
     &                   DMAGAOSR(1,2),THRINT,IPRINT)
         
C        Contribution to spin density (N_alpha - N_beta)
         SPINDEN = SPINDEN + WGHT*RHO(2)
      ELSE
         RHO(2) = 0.0D0
      END IF

C     Gradient of density
      IF (DOGGA) THEN
C
        CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAOSR(1,1),
     &     1,D0,RHG(1,1),1)
C
        IF (DFT_SPINDNS) THEN
           CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAOSR,
     &        1,D0,RHG(1,2),1)
        ELSE
           RHG(1,2) = 0.0D0 ! RHGS is zero
           RHG(2,2) = 0.0D0
           RHG(3,2) = 0.0D0
        ENDIF
        DO I = 1,3
           RHG(I,3) = RHG(I,1) + RHG(I,2) ! RHG2A
           RHG(I,4) = RHG(I,1) - RHG(I,2) ! RHG2B
        END DO
      ENDIF
C
CFRAN I DO THE HESSIAN LATER!!
C     Hessian of density (for second-order optimization)
      IF (DOATR) THEN
       CALL DFTRHH(DMAT,DMAGAOSR,GSO(KSO0),GSO(KSO1),
     >            GSO(KSO2),RHG,RHOLAP,RHOGHG)
      ENDIF

C     Number of electrons
      ELCTRN = ELCTRN + DELCTRN


      IF (DOERG) THEN

C
C     Exchange and correlation energies and derivatives
      NORDER = 1  ! first-order derivatives
C
      CALL SRDFTEXCJT(RHO,RHG,CHIVAL,NORDER,ERFEXP,
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
       IF (DFT_SPINDNS) THEN
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
         rhogrdt = sqrt( rhg(1,1)**2 + rhg(2,1)**2 + rhg(3,1)**2 )
         CALL SRDFTLTRJT(JWOPSY,DMAT(1,JDXAO),EXCMAT,
     &                 WGHTX,WGHTC, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                 ERFEXP,RHO(1),RHO(2),
     &                 RHOGRDT,RHOLAP,RHOGHG,RHG,DMAGAOSR(1,1),DOGGA)
       ENDIF

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

      IF (DFT_SPINDNS) THEN
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
      IF (IPRINT .GE. 2) THEN
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
      END IF 

      IF (IPRINT.GE.5) THEN
            WRITE(LUPRI,'(//A/)') ' <<< Kohn-Sham matrix >>>'
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
            IF (DFT_SPINDNS) THEN
               WRITE(LUPRI,'(//A/)')
     &         ' <<< Kohn-Sham matrix from spin-density>>>'
               CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,NBAST,NBAST,NBAST,
     &         -1,LUPRI)
            END IF
      END IF

      END IF  ! IPRINT .GE. 2

      RETURN 
      END

C****************************************************************************
      SUBROUTINE SRDFTEXCJT(RHO,RHG, CHIVAL,NORDER,ERFEXP,
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

      double precision rho(4),rhg(3,4),CHIVAL ! input
      double precision rhogrd(4) ! local variable for norm of RHG(1:3,*)
      integer NORDER
      logical ERFEXP
      double precision Ex,d1Exdrhot,d1Exdrhos,d1Exdgrd,d1Exdgrda2
      double precision d2Exdrhotdrhot,d2Exdrhotdrhos,d2Exdrhosdrhos
      double precision d2Exdgrddgrd,d2Exdrhodgrd
      double precision Ec,d1Ecdrhot,d1Ecdrhos,d1Ecdgrd,d1Ecdgrda2
      double precision d2Ecdrhotdrhot,d2Ecdrhotdrhos,d2Ecdrhosdrhos
      double precision d2Ecdgrddgrd,d2Ecdrhodgrd
C The new terms
      double precision Exa,Exb d1Exdrhoa,d1Exdrhob,d1Exdgrda
      double precision d1Exdgrdb,d1Exdgrda2a,d1Exdgrda2b
      double precision d2Exdrhoadrhoa,d2Exdrhobdrhob,d2Exdrhoadgrda
      double precision d2Exdrhobdgrdb
      save icall
      data icall/0/
C     functionals
      external ESRX_PBEHSEERF,  ESRX_PBETCSERF, ESRX_PBERIERF
      external ESRX_PBEGWSERF
      external ESRC_PBETCSJERF, ESRC_PBETCSERF, ESRC_PBERIERF
      external ESRC_PBEWIERF,   ESRC_PBEGWSERF
C     spin functionals
      external ESRX_LDAERFSPIN, ESRC_LDAERFSPIN

c  Exchange  ------------------------------------------------------------

      icall = icall + 1
      IF (ERFEXP) THEN
         write (lupri,*)
     &   'WARNING srdftexcjt: erfexp true in call no. ',icall
         erfexp = .false.
!        CALL QUIT('erfgau not implemented in SRDFTEXCJT')
      END IF

      IF (DFT_SPINDNS) THEN
         do i = 1,4
            rhogrd(i) = sqrt(rhg(1,i)**2 + rhg(2,i)**2 + rhg(3,i)**2)
         end do
      ELSE
         rhogrd(1) = sqrt(rhg(1,1)**2 + rhg(2,1)**2 + rhg(3,1)**2)
      END IF

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
CFRAN Extra terms initialized!
      Exa = 0.d0
      Exb = 0.d0
      d1Exdrhoa  = 0.d0
      d1Exdrhob  = 0.d0
      d1Exdgrda = 0.d0
      d1Exdgrdb = 0.d0
      d1Exdgrda2a  = 0.d0
      d1Exdgrda2b  = 0.d0
      d2Exdrhoadrhoa = 0.d0
      d2Exdrhobdrhob = 0.d0
      d2Exdrhoadgrda = 0.d0
      d2Exdrhobdgrdb = 0.d0
CFRAN END

!     short-range exchange PBE of HSE
      IF(DOSRX_PBEHSE) THEN
        IF(.not. DFT_SPINDNS) THEN     
           call ESRX_PBEHSEERF(rho,RHOGRD,CHIVAL,Ex)
           call DESR(ESRX_PBEHSEERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
        ELSE
! It should be divided by half but since we have to multiply 
! later by 2 ...   
           call ESRX_PBEHSEERF(rho(3),RHOGRD(3),CHIVAL,Exa)
!          call ESRX_PBEHSEERF(rhoa2,RHOGRDA2,CHIVAL,Exa)
           call ESRX_PBEHSEERF(rho(4),RHOGRD(4),CHIVAL,Exa)
!          call ESRX_PBEHSEERF(rhob2,RHOGRDB2,CHIVAL,Exb)
           call DESR(ESRX_PBEHSEERF,rho(3),RHOGRD(3),CHIVAL,NORDER,
     >             d1Exdrhoa,d1Exdgrda,d1Exdgrda2a,
     >             d2Exdrhoadrhoa,d2Exdgrdadgrda,d2Exdrhoadgrda)
           call DESR(ESRX_PBEHSEERF,rho(4),RHOGRD(4),CHIVAL,NORDER,
     >             d1Exdrhob,d1Exdgrdb,d1Exdgrda2b,
     >             d2Exdrhorbdrhob,d2Exdgrdbdgrdb,d2Exdrhobdgrdb)
C
CFRAN Get TOTAL first and second derivatives TODO
C
           d1Exdrhot      = 0.5d0*(d1Exdrhoa+ d1Exdrhob)
           d1Exdgrd       = 0.5d0*(d1Exdgrda+ d1Exdgrdb)
           d1Exdgrda2     = 0.5d0*(d1Exdgrda2a    + d1Exdgrda2b)
           d2Exdrhotdrhot = 0.5d0*(d2Exdrhoadrhoa + d2Exdrhorbdrhob)
           d2Exdgrddgrd   = 0.5d0*(d2Exdgrdadgrda + d2Exdgrdbdgrdb)
           d2Exdrhodgrd   = 0.5d0*(d2Exdrhoadgrda + d2Exdrhobdgrdb)
C
CFRAN Get SPIN dependent first and second derivatives (???) TODO
C
           d1Exdrhos      = 0.5d0*(d1Exdrhoa      - d1Exdrhob)
           d2Exdrhosdrhos = 0.5d0*(d2Exdrhoadrhoa - d2Exdrhorbdrhob)
C
CFRAN The exchange energy...
C
           Exa = 0.5d0*Exa
           Exb = 0.5d0*Exb
           Ex = Exa + Exb 
        ENDIF
      ENDIF !  (DOSRX_PBEHSE)

!     short-range exchange PBE of TCS
      IF(DOSRX_PBETCS) THEN
        call ESRX_PBETCSERF(rho,RHOGRD,CHIVAL,Ex)
        call DESR(ESRX_PBETCSERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
      ENDIF

!     short-range exchange PBE of GWS ! JT 11-08-09
      IF(DOSRX_PBEGWS) THEN
        call ESRX_PBEGWSERF(rho,RHOGRD,CHIVAL,Ex)
        call DESR(ESRX_PBEGWSERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
      ENDIF

C  Manu 01-02-2006 begin
!     short-range exchange RI
      IF(DOSRX_PBERI) THEN
        call ESRX_PBERIERF(rho,RHOGRD,CHIVAL,Ex)
        call DESR(ESRX_PBERIERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Exdrhot,d1Exdgrd,d1Exdgrda2,
     >             d2Exdrhotdrhot,d2Exdgrddgrd,d2Exdrhodgrd)
      ENDIF
C  Manu 01-02-2006 end

!     short-range spin-dependent exchange LDA ! JT 11-08-09
      IF(DOSRX_LDAS) THEN
        call ESRX_LDAERFSPIN(rho(1),rho(2),CHIVAL,DFT_SPINDNS,Ex)
        call DESRSPIN(ESRX_LDAERFSPIN,rho(1),rho(2),
     >             CHIVAL,DFT_SPINDNS,NORDER,
     >             d1Exdrhot,d1Exdrhos,
     >             d2Exdrhotdrhot,d2Exdrhotdrhos,d2Exdrhosdrhos)
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
!        call VSRCPBEERF(ec,d1ecdrho,d1ecdgrd,d1ecdgrda2,d2ecdrhodrho,
!     >                  d2ecdgrddgrd,d2ecdrhodgrd,rho,rhogrd,chival,
!     >                  norder)
        call ESRC_PBETCSERF(rho,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBETCSERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
      ENDIF

!     short-range correlation PBE of TCSJ
      IF(DOSRC_PBETCSJ) THEN
        call ESRC_PBETCSJERF(rho,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBETCSJERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
      ENDIF

!     short-range correlation PBE of GWS
      IF(DOSRC_PBEGWS) THEN
        call ESRC_PBEGWSERF(rho,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBEGWSERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
      ENDIF

!     rational interpolation
      IF(DOSRC_PBERI) THEN
        call ESRC_PBERIERF(rho,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBERIERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
      ENDIF

!     weighted interpolation
      IF(DOSRC_PBEWI) THEN
        call ESRC_PBEWIERF(rho,RHOGRD,CHIVAL,Ec)
        call DESR(ESRC_PBEWIERF,rho,RHOGRD,CHIVAL,NORDER,
     >             d1Ecdrhot,d1Ecdgrd,d1Ecdgrda2,
     >             d2Ecdrhotdrhot,d2Ecdgrddgrd,d2Ecdrhodgrd)
      ENDIF

!     short-range spin-dependent correlation LDA ! JT 18-08-09
      IF(DOSRC_LDAS) THEN
        call ESRC_LDAERFSPIN(rho(1),rho(2),CHIVAL,DFT_SPINDNS,Ec)
        call DESRSPIN(ESRC_LDAERFSPIN,rho(1),rho(2),
     >             CHIVAL,DFT_SPINDNS,NORDER,
     >             d1Ecdrhot,d1Ecdrhos,
     >             d2Ecdrhotdrhot,d2Ecdrhotdrhos,d2Ecdrhosdrhos)
      ENDIF

      RETURN
      END

C**************************************************************************
      SUBROUTINE SRDFTLTRJT(KSYMOP,DTRMAT,EXCMAT,WGHTX,WGHTC,GAO,GAO1,
     &                    CHIVAL,ERFEXP,RHOT,RHOS,RHOGRD,RHOLAP,RHOGHG,
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
CFRAN RHOSGRD ADDED!
      CALL SRDFTEXCJT(RHO,RHG,CHIVAL,NORDER,ERFEXP,
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
