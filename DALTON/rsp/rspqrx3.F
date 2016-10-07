!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
#define DFTQRC dft_qr_respons
c#define DFTQRC dftqrcf
C
#ifdef REVLOG
C=========================================================================
CRevision 1.2  2000/05/24 19:06:44  hjj
Cnew getref calls with appropriate NDREF instead of NCREF
C(fixing error for triplet with CSF)
C941223-hjaaj: new QFOCK call
C931015-hjaaj: GTZYMT, TRZYMT: stop if NSIM .ne. 1
C===========================================================================
#endif
C  /* Deck t3drv */
      SUBROUTINE T3DRV(NSIM,ISYMA,ISYMB,ISYMC,VECB,VECC,ATEST,VECA,
     *                  OMEGAB,OMEGAC,XINDX,UDV,PV,MJWOP,WRK,LWRK,
     *                  CMO,FC,FV)
C
C     Purpose:
C     Compute E[3] times two vectors
C     Compute S[3] times two vectors, multiply with omega, and add to result
C
C     OUTPUT is put in WRK:
C     First  NSIM * KZYVA elements are result vector
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0 , DM1 = -1.0D0 , DEM10 = 1.0D-10 )
C
#include "priunit.h"
#include "thrzer.h"
#include "infvar.h"
#include "qrinf.h"
#include "rspprp.h"
#include "infhyp.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
#include "infspi.h"
#include "gnrinf.h"
C
      LOGICAL ATEST, TEST
      DATA TEST /.FALSE./
C
      DIMENSION WRK(*)
      DIMENSION XINDX(*)
      DIMENSION UDV(*)
      DIMENSION PV(*), MJWOP(2,MAXWOP,8)
      DIMENSION VECB(*), VECC(*), VECA(*)
      DIMENSION CMO(*),FC(*)
C
      CALL QENTER('T3DRV')
C
C     Compute the length of storage needed
C     KZYVA is the length of the resulting vector
C
      KZYVA  = MZYVAR(ISYMA)
      KZYVB  = MZYVAR(ISYMB)
      KZYVC  = MZYVAR(ISYMC)
C
      IBEQC = 0
      IF ( HYPCAL .AND. ISYMB .EQ. ISYMC .AND. ISPINB.EQ.ISPINC .AND.
     *     ABS(OMEGAB) .EQ. ABS(OMEGAC) .AND. .NOT. RSPCI ) THEN
         IF (OMEGAB .EQ. OMEGAC) THEN
            BMINC = D0
            DO 110 I = 1,KZYVB
               BMINC = BMINC + ABS(VECB(I) - VECC(I))
  110       CONTINUE
            IF (ATEST) WRITE(LUPRI,*) '1norm(B - C) = BMINC =',BMINC
            IF (BMINC .LE. THRZER) IBEQC = 1
         ELSE
            BMINC = D0
            KZVB = MZVAR(ISYMB)
            DO 120 I = 1,KZVB
               BMINC = MAX(BMINC,ABS(VECB(I) + VECC(KZVB+I)
     *                             + VECC(I) + VECB(KZVB+I)) )
  120       CONTINUE
            IF (ATEST) WRITE(LUPRI,*) 'max(B + Ct) = BMINC =',BMINC
            IF (BMINC .LE. DEM10) IBEQC = -1
C        MAERKE: note that this will not catch dipole velocity type
C        where VECB(I) = VECC(KZVB+I) /910905-pj+hjaaj
         END IF
      END IF
C
C Assumed order of spin operators in triplet quadratic response
C
C
      IF( IPRRSP .GT. 50) THEN
         WRITE(LUPRI, '(//A)')
     *   ' Characteristics in T3DRV routine'
         WRITE(LUPRI, '( //A,2I10 )')
     *   ' Symmetry of vectors b and c  : ', ISYMB,ISYMC
         WRITE(LUPRI, '( A,2I10 )')
     *   ' Spin of vectors b and c      : ', ISPINB,ISPINC
         WRITE(LUPRI, '( A,2I10 )')
     *   ' Length of orbital part       : ',
     *     MZWOPT(ISYMB),MZWOPT(ISYMC)
         WRITE(LUPRI, '( A,2I10 )')
     *   ' Length of configuration part : ', MZCONF(ISYMB),MZCONF(ISYMC)
         WRITE(LUPRI, '( A,2I10 )')
     *   ' Length of vectors            : ', KZYVB,KZYVC
         WRITE(LUPRI, '( A,2F10.7 )')
     *   ' Frequencies of vectors       : ', OMEGAB,OMEGAC
         WRITE(LUPRI, '( A,I10 )')
     *   ' vector b equal to vector c ? : ', IBEQC
      END IF
C
      KE3    = 1
      KS3    = KE3  + KZYVA
      KWRK   = KS3  + KZYVA
      LWRKE  = LWRK - KS3
      LWRKS  = LWRK - KWRK
      IF (LWRKS.LT.0) CALL ERRWRK('T3DRV',KWRK-1,LWRK)
C
C     (Keep WRK as the true workspace and split off the computed arrays)
C
      IF( (NSIM .GT. 1) ) THEN
         WRITE(LUERR,'(//A,/A,I5,/A)')
     *   ' --> Error : Routine T3DRV not implemented  with NSIM > 1',
     *   '     NSIM = ', NSIM,
     *   '     Halting execution.'
          CALL QTRACE(LUERR)
          CALL QUIT(' T3DRV ERROR: NSIM GREATER THAN 1 IN T3DRV')
      END IF
C
C     Initialise the result vector
C
      CALL DZERO(WRK(KE3),KZYVA)
C
C     Return if we have a CI response calculation
C
      IF ( RSPCI ) RETURN
C
      IF ( TEST ) THEN
         CALL DCOPY(MZVAR(ISYMB),VECB,1,VECB(1+MZVAR(ISYMB)),1)
         CALL DSCAL(MZVAR(ISYMB),DM1,VECB(1+MZVAR(ISYMB)),1)
         CALL DCOPY(MZVAR(ISYMC),VECC,1,VECC(1+MZVAR(ISYMC)),1)
         CALL DSCAL(MZVAR(ISYMC),DM1,VECC(1+MZVAR(ISYMC)),1)
      END IF
C
C     Compute ( E3(jkl) + E3(jlk) ) N(l) N(k)
C
      CALL E3INIT(VECB,VECC,VECA,ATEST,IBEQC,
     *            WRK(KE3),XINDX,UDV,PV,WRK(KS3),LWRKE,
     *            KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,
     *            ISPINA,ISPINB,ISPINC,CMO,FC,FV,MJWOP)
C
C     Compute S3(jlk) N(l) N(k)
C
      IF (OMEGAB .NE. D0 ) THEN
         CALL S3INIT(KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,
     *               ISPINA,ISPINB,ISPINC,
     *               VECB,VECC,VECA,ATEST,
     *               WRK(KS3),XINDX,UDV,MJWOP,WRK(KWRK),LWRKS)
C
C        Multiply with omega-b
C
         CALL DAXPY(KZYVA,OMEGAB,WRK(KS3),1,WRK(KE3),1)
      END IF
C
C     Compute S3(jkl) N(l) N(k)
C
      IF (OMEGAC .NE. D0 ) THEN
         CALL S3INIT(KZYVA,KZYVC,KZYVB,ISYMA,ISYMC,ISYMB,
     *               ISPINA,ISPINC,ISPINB,
     *               VECC,VECB,VECA,ATEST,
     *               WRK(KS3),XINDX,UDV,MJWOP,WRK(KWRK),LWRKS)
C
C        Multiply with omega-c
C
         CALL DAXPY(KZYVA,OMEGAC,WRK(KS3),1,WRK(KE3),1)
      END IF
C
C     Result is now in WRK(KE3)
C
      IF (IPRRSP .GT. 100 ) THEN
         KZVA = MZVAR(ISYMA)
         WRITE(LUPRI,'(//A/A)')
     &     ' Final result in T3DRV    (Col 1 = Z,  Col 2 = Y)',
     &     ' ======================'
         CALL OUTPUT(WRK(KE3),1,KZVA,1,2,KZVA,2,1,LUPRI)
      END IF
C
      CALL QEXIT('T3DRV')
      RETURN
      END
C  /* Deck dft3drv */
      SUBROUTINE DFT3DRV(VECB, VECC,FI,FA,ZYMB,ZYMC,
     &                   KZYVA,KZYVB,KZYVC,
     &                   ISYMA,ISYMB,ISYMC,
     &                   ISPINA,ISPINB,ISPINC,
     &                   CMO,MJWOP,WRK,LFREE,ADDFOCK)
#include "implicit.h"
#include "infvar.h"
#include "inforb.h"
#include "maxorb.h"
#include "infinp.h"
#include "dftcom.h"
      DIMENSION VECB(*),VECC(*),ZYMB(*),ZYMC(*),WRK(*),CMO(*)
      DIMENSION FI(NORBT,NORBT),FA(NORBT,NORBT),MJWOP(2,MAXWOP,8)
      LOGICAL ADDFOCK
c     unpack response vectors.
      NSIM = 1
      CALL GTZYMT(NSIM,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      CALL GTZYMT(NSIM,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
c     compute dft contribution
C ZR testing dft_qr_ab
      IF (NASHT.GT.0) THEN
         CALL dft_qr_ab(FI,FA,CMO,ZYMB,ISYMB,ISPINB,ZYMC,ISYMC,ISPINC,
     &                  ADDFOCK,WRK,LFREE,IPRDFT)
      ELSE
         CALL DFTQRC(FI,CMO,ZYMB,ISYMB,ISPINB,ZYMC,ISYMC,ISPINC,
     &               ADDFOCK,WRK,LFREE,IPRDFT)
      END IF
c     pack stuff back.
      END
C  /* Deck e3init */
      SUBROUTINE E3INIT(VECB,VECC,VECA,ATEST,IBEQC,
     *                  ETRS,XINDX,UDV,PV,WRK,LWRK,KZYVA,KZYVB,KZYVC,
     *                  ISYMA,ISYMB,ISYMC,ISPINA,ISPINB,ISPINC,CMO,FC,
     *                  FV,MJWOP)

      use qmcmm_response, only: qmcmm_qr
      use pelib_interface, only: use_pelib, pelib_ifc_qro,
     &                           pelib_ifc_gspol
C
C
C
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( D1 = 1.0D0 )
C
C infinp.h : HSROHF, DODFT, ?
C
#include "infdim.h"
#include "maxorb.h"
#include "maxash.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "rspprp.h"
#include "infcr.h"
#include "inftpa.h"
#include "inftmo.h"
#include "pcmlog.h"
#include "infpp.h"
#include "infsmo.h"
#include "infhyp.h"
#include "esg.h"
#include "dftcom.h"
#include "absorp.h"
#include "abslrs.h"
#include "abscrs.h"
#include "gnrinf.h"
#include "infpar.h"
C
      LOGICAL   ATEST, TRPLETSAVE, DFTSAV, ADDFOCK
      DIMENSION VECB(*),VECC(*),VECA(*)
      DIMENSION ETRS(KZYVA)
      DIMENSION XINDX(*)
      DIMENSION UDV(*), MJWOP(2,MAXWOP,8)
      DIMENSION PV(*)
      DIMENSION WRK(*)
      DIMENSION CMO(*),FC(*), FV(*)
C
C     Initialise variables and layout some workspace
C
      CALL QENTER('E3INIT')
      IF (DOHFSRDFT .OR. DOMCSRDFT) THEN
         CALL QUIT('srDFT not implemented in E3INIT')
      END IF
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KZYMB,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZYMC,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC,NORBT*NORBT,WRK,KFREE,LFREE)
      IF (NASHT.GT.0) THEN
         CALL MEMGET('REAL',KFO,NORBT*NORBT,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KFO,0,WRK,KFREE,LFREE)
      END IF
C
C Zero Q matrix for ORBEX (and FO for closed shell)
C
      CALL MEMGET('REAL',KDUM,N2ORBX,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KDUM),N2ORBX)
      IF (NASHT.EQ.0) KFO=KDUM
C
C Fock matrices
C
      IF (TDHF.AND.(NASHT.EQ.0.OR.HSROHF)) THEN ! excluding single open shell RHF with .OPEN SHELL
         IF (CRCAL.OR.TOMOM.OR.TPAMP.OR.ESG) THEN
            CALL C3FOCK(IBEQC,WRK(KFC),VECB,VECC,
     &                 WRK(KZYMB),WRK(KZYMC),
     &                 KZYVB,KZYVC,
     &                 ISYMA,ISYMB,ISYMC,
     &                 CMO,FC,MJWOP,WRK(KFREE),LFREE)
         ELSE IF (HYPCAL.OR.SOMOM.OR.EXMOM.OR.ABS_BETA
     &        .OR.ABS_MCD.OR.ABSLRS_MCD.OR.ABSLRS_MCHD 
     &        .OR.ABSLRS_NSCD.OR.ABSLRS_GAMMA.OR.ABSLRS_IDRI) THEN

            DFTSAV=DODFT
            DODFT=.FALSE.
            CALL Q3FOCK( VECB,VECC,
     &         KZYVB,KZYVC, ISYMB,ISYMC, ISPINB,ISPINC,
     &         CMO, MJWOP, FC, FV,
     &         WRK(KFC),WRK(KFO),
     &         WRK(KFREE),LFREE )
            DODFT=DFTSAV
         ELSE
            CALL QUIT('ERROR in E3INIT for HF or DFT')
         END IF
C
         TRPLETSAVE=TRPLET
         TRPLET=MOD(ISPINA+ISPINB+ISPINC,2).NE.0
         CALL ORBEX(
     &         ISYMA,1,KZYVA,WRK(KFC),WRK(KFO),WRK(KDUM),WRK(KDUM),
     &         ETRS,D1,UDV,MJWOP
     &         )  
         TRPLET=TRPLETSAVE
C
C Special consideration, ingoing FC does contain the DFT contribution
C which is transformed in Q3FOCK as [k,k,FC] so the corresponding
C contribution must be skipped in DFT3DRV
C
         ADDFOCK=.FALSE.
C
      ELSE
         ADDFOCK=.NOT.DIRFCK
         CALL MEMGET('REAL',KDEN1,NASHT*NASHT,WRK,KFREE,LFREE)
         IF (TRPFLG) THEN
            CALL MEMGET('REAL',KDEN2,2*N2ASHX*N2ASHX,WRK,KFREE,LFREE)
         ELSE
            CALL MEMGET('REAL',KDEN2,N2ASHX*N2ASHX,WRK,KFREE,LFREE)
         END IF
         CALL E3DRV(VECB, VECC, VECA, ATEST, IBEQC,
     *              ETRS,XINDX,WRK(KFC),WRK(KZYMB),WRK(KZYMC),
     *              WRK(KDEN1),WRK(KDEN2),UDV,PV,WRK(KFREE),LFREE,
     *              KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,FC,MJWOP)
      END IF
C
C DFT contribution
C
      IF(DODFT) THEN
         CALL DZERO(WRK(KFC),N2ORBX)
         IF (NASHT.GT.0) CALL DZERO(WRK(KFO),N2ORBX)
         CALL DFT3DRV(VECB,VECC,WRK(KFC),WRK(KFO),WRK(KZYMB),WRK(KZYMC),
     *                   KZYVA,KZYVB,KZYVC,
     *                   ISYMA,ISYMB,ISYMC,
     *                   ISPINA,ISPINB,ISPINC,
     *                   CMO,MJWOP,WRK(KFREE),LFREE,ADDFOCK)
         IF (E3TEST) THEN
            E3 = -DDOT(KZYVA,ETRS,1,VECA,1)
         END IF
         TRPLETSAVE=TRPLET
         TRPLET=MOD(ISPINA+ISPINB+ISPINC,2).NE.0
         CALL ORBEX(ISYMA,1,KZYVA,WRK(KFC),WRK(KFO),WRK(KDUM),
     &        WRK(KDUM),ETRS,D1,UDV,MJWOP)
         TRPLET=TRPLETSAVE
         IF (E3TEST) THEN
            E3TOT = -DDOT(KZYVA,ETRS,1,VECA,1)
            E3DFT = E3TOT - E3
            WRITE(LUPRI,'(/A,2F20.8)')
     *         '   COU CONTRIBUTION TO HYPVAL',E3,E3
            WRITE(LUPRI,'(/A,2F20.8)')
     *         '   DFT CONTRIBUTION TO HYPVAL',E3DFT,E3TOT
         END IF
      END IF

      IF (USE_PELIB()) THEN
         IF (.NOT. TDHF .AND. .NOT. TRPLET) THEN 
!           write(lupri,*) '*** WARNING PE-MCSCF QR ***  '
!           write(lupri,*) 'Highly experimental version'
!           write(lupri,*) 'Results *WILL* be wrong'
            CALL QUIT('ERROR: PE-MCSCF quadratic response not'//
     &                ' implemented using the PE library.')
!            CALL PE_RSPMCQR(VECB, VECC, ETRS, XINDX, 
!     &                      WRK(KZYMB), WRK(KZYMC), WRK(KDEN1), UDV, 
!     &                      WRK(KFREE), LFREE, KZYVA, KZYVB, KZYVC,
!     &                      ISYMA, ISYMB, ISYMC, CMO, MJWOP)
         ELSE IF (NASHT > 0) THEN
            CALL QUIT('ERROR: quadratic response not implemented for'//
     &                ' open-shell systems using the PE library.')
         ELSE IF (TRPLET) THEN
            CALL QUIT('ERROR: triplets not implemented for quadratic'//
     &                ' response using the PE library.')
         ENDIF
         IF (.NOT. PELIB_IFC_GSPOL()) THEN
          CALL PELIB_IFC_QRO(VECB, VECC, ETRS, XINDX, WRK(KZYMB),
     &                       WRK(KZYMC), UDV, WRK(KFREE), LFREE,
     &                       KZYVA, KZYVB, KZYVC, ISYMA, ISYMB,
     &                       ISYMC, CMO, MJWOP)
         END IF
      END IF

      IF (QM3) THEN
         IF (TDHF) THEN
#if defined(VAR_MPI)
            IF (NODTOT.GE.1) THEN
               CALL QM3QRO_P(VECB,VECC,
     *                 ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,MJWOP)
               CALL QM3QRO_P(VECC,VECB,
     *                 ETRS,XINDX,WRK(KZYMC),WRK(KZYMB),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVC,KZYVB,ISYMA,ISYMC,ISYMB,CMO,MJWOP)
            ELSE
#endif
               CALL QM3QRO(VECB,VECC,
     *                 ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,MJWOP)
               CALL QM3QRO(VECC,VECB,
     *                 ETRS,XINDX,WRK(KZYMC),WRK(KZYMB),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVC,KZYVB,ISYMA,ISYMC,ISYMB,CMO,MJWOP)
#if defined(VAR_MPI)
            ENDIF
#endif
         END IF
      END IF

      IF (QMMM) THEN
         IF (TDHF) THEN
            CALL QMMMQRO(VECB,VECC,
     *                   ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),
     *                   UDV,WRK(KFREE),LFREE,
     *                   KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,MJWOP,
     *                   ISPINA,ISPINB,ISPINC)
            CALL QMMMQRO(VECC,VECB,
     *                   ETRS,XINDX,WRK(KZYMC),WRK(KZYMB),
     *                   UDV,WRK(KFREE),LFREE,
     *                   KZYVA,KZYVC,KZYVB,ISYMA,ISYMC,ISYMB,CMO,MJWOP,
     *                   ISPINA,ISPINC,ISPINB)
         END IF
      END IF
C
      IF (QMNPMM) THEN
         IF (TDHF) THEN
            CALL qmcmm_qr(VECB,VECC,ETRS,XINDX,WRK(KZYMB),
     *                        WRK(KZYMC),UDV,WRK(KFREE),LFREE, KZYVA,
     *                        KZYVB,KZYVC, ISYMA,ISYMB,ISYMC,CMO,MJWOP,
     *                        ISPINA,ISPINB,ISPINC)
               CALL qmcmm_qr(VECC,VECB,ETRS,XINDX,WRK(KZYMC),
     *                        WRK(KZYMB), UDV,WRK(KFREE),LFREE,KZYVA,
     *                        KZYVC,KZYVB,ISYMA,ISYMC,ISYMB,CMO,MJWOP,
     *                        ISPINA,ISPINC,ISPINB) 
            END IF   
       END IF
C
      IF (FLAG(16)) THEN
         CALL MEMGET('INTE',KSYRLM,(2*LSOLMX+1),WRK,KFREE,LFREE)
         IF (TDHF) THEN
            CALL C3SOL(VECA, VECB, VECC,
     *                 ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,MJWOP,
     *                 WRK(KSYRLM))
            CALL C3SOL(VECA, VECC, VECB,
     *                 ETRS,XINDX,WRK(KZYMC),WRK(KZYMB),
     *                 UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVC,KZYVB,ISYMA,ISYMC,ISYMB,CMO,MJWOP,
     *                 WRK(KSYRLM))
         ELSE
            CALL E3SOL(VECA, VECB, VECC,
     *                 ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),
     *                 WRK(KDEN1),UDV,WRK(KFREE),LFREE,
     *                 KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,MJWOP,
     *                 WRK(KSYRLM))
         END IF
C
C
Clf 10/06/2002
C quadratic PCM response (for the time beeing without memory efficent SCF)
C
      ELSE IF (PCM) THEN
         IF (.NOT.NEWQR) THEN
            CALL E3IEF(VECA, VECB, VECC,
     *           ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),
     *           WRK(KDEN1),UDV,WRK(KFREE),LFREE,
     *           KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,MJWOP)
         ELSE
            CALL E3IEF2(VECA, VECB, VECC,
     *           ETRS,XINDX,WRK(KZYMB),WRK(KZYMC),
     *           WRK(KDEN1),UDV,WRK(KFREE),LFREE,
     *           KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,CMO,MJWOP)
         END IF
      END IF
C
      CALL QEXIT('E3INIT')
      RETURN
      END
C  /* Deck e3drv */
      SUBROUTINE E3DRV(VECB, VECC, VECA, ATEST, IBEQC,
     *                 ETRS,XINDX,FI,ZYMB,ZYMC,DEN1,DEN2,
     *                 UDV,PV,WRK,LWRK,KZYVA,KZYVB,KZYVC,
     *                 ISYMA,ISYMB,ISYMC,CMO,FC,MJWOP)
C
C      Purpose:
C      Outer driver routine for E[3] times two vectors. This subroutine
C      calls the setup of the integral transformation, constructs a den-
C      sity matrix if necessary  and calls RSP2GR to compute the gradient.
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, DH = 0.5D0 , D1 = 1.0D0 , D2 = 2.0D0 )
      PARAMETER ( DM1=-1.0D0 )
C
#include "priunit.h"
#include "infdim.h"
#include "inforb.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infpri.h"
#include "infvar.h"
#include "qrinf.h"
#include "infspi.h"
#include "infden.h"
C
      DIMENSION ETRS(KZYVA), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION PV(N2ASHX*N2ASHX,*)
      DIMENSION DEN1(NASHDI,NASHDI), DEN2(N2ASHX*N2ASHX,*)
      DIMENSION FI(NORBT,NORBT)
      DIMENSION ZYMB(*), ZYMC(*)
      DIMENSION WRK(*)
      DIMENSION CMO(*),FC(*)
      DIMENSION VECB(KZYVB), VECC(KZYVC), VECA(KZYVA)
C
      LOGICAL   ATEST, LCON, LORB, LREFST
      LOGICAL   TDM, NORHO2, LREF
      LOGICAL   LORBB, LORBC, LCONB, LCONC
C
C     Initialise variables and layout some workspace
C
      CALL QENTER('E3DRV')
      IGRSYM = ISYMA
      IGRSPI = ISPINA
      TDM    = .TRUE.
      NORHO2 = .FALSE.
      IDAX0  = 0
C
C     for E3TEST
C
      E3TMZC = D0
      E3TMZO = D0
      E3TMYC = D0
      E3TMYO = D0
      E3TERM = D0
      KZCA   = MZCONF(ISYMA)
      KZOA   = MZWOPT(ISYMA)
      KAZC   = 1
      KAZO   = KAZC + KZCA
      KAYC   = KAZO + KZOA
      KAYO   = KAYC + KZCA
C
      KH1     = 1
      KH1X    = KH1     + N2ORBX
      KWOIT   = KH1X    + N2ORBX
      LWOIT   = LWRK    - KWOIT
      IF (LWOIT.LT.0) CALL ERRWRK('E3DRV',KWOIT-1,LWRK)
      KFROIT  = 1
      KCREF   = KWOIT
      KFREE   = KCREF   + MZCONF(1)
      LFREE   = LWRK    - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('E3DRV',KFREE-1,LWRK)
C
      IPROIT  = MAX(IPRRSP/5,2)
C
      CALL RSPH1(WRK(KH1),CMO,WRK(KWOIT),LWOIT)
      IF (IPRRSP. GT. 200 ) THEN
         WRITE(LUPRI,'(//A)') ' One electron matrix in E3DRV'
         WRITE(LUPRI,'(A)')   ' ============================'
         CALL OUTPUT(WRK(KH1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C Gradient flags
C

      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( IGRSYM .EQ. 1 ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      IF ( ISYMB .EQ. 1 ) THEN
         LCONB = ( MZCONF(ISYMB) .GT. 1 )
      ELSE
         LCONB = ( MZCONF(ISYMB) .GT. 0 )
      END IF
      IF ( ISYMC .EQ. 1 ) THEN
         LCONC = ( MZCONF(ISYMC) .GT. 1 )
      ELSE
         LCONC = ( MZCONF(ISYMC) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LORBB  = ( MZWOPT(ISYMB) .GT. 0 )
      LORBC  = ( MZWOPT(ISYMC) .GT. 0 )
C
C     Case 1:
C     =======
C     IF ( MZWOPT(ISYMB) .EQ. 0 ) GO TO 3000
      IF (.NOT. LORBB) GO TO 3000
C     Transform the integrals with kappa(1)
C
      NSIM = 1
      CALL GTZYMT(NSIM,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      JTRLVL= 2
      IF (.NOT.DIROIT) THEN
         CALL RSPOIT(IDAX0,IDAX1,JTRLVL,ISYMB,ZYMB,
     *            WRK(KWOIT),KFROIT,LWOIT,IPROIT)
      END IF
      CALL DZERO(WRK(KH1X),N2ORBX)
      CALL OITH1(ISYMB,ZYMB,WRK(KH1),WRK(KH1X),1)
C     CALL OITH1(LSYOIT,OITMAT,H1,H1X,IH1SYM)
      CALL DCOPY(N2ORBX,WRK(KH1X),1,FI,1)
C
C     IF (MZCONF(ISYMC) .EQ. 0) GO TO 2000
      IF (.NOT.LCONC) GO TO 2000
C
C     Construct the density matrix <02L|..|0> + <0|..|02R>
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMC)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMC)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMC)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      NSIM = 1
      ISPIN1 = MULSP(IGRSPI,ISPINB)
      ISPIN2 = 0
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VECC,OVLAP,DEN1,DEN2,ISPIN1,ISPIN2,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
      IF (IGRSPI.EQ.1) THEN
         ISPIN1 = ISPINB
         ISPIN2 = IGRSPI
         CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VECC,OVLAP,DEN1,DEN2(1,2),ISPIN1,ISPIN2,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
      END IF
      CALL IPSET(MULSP(IGRSPI,ISPINB),0,ISPINB,IGRSPI)
C
C     Make the gradient
C
      ISYMV  = ISYMC
      INTSYM = ISYMB
      LREFST = .FALSE.
      NZYVEC = MZYVAR(ISYMC)
      NZCVEC = MZCONF(ISYMC)
      IKLVL = 1
      IF (DIROIT) THEN
         ISPIN1 = ISPINB
         ISPIN2 = 0
         IDAX = IDAX0
         KINTSY = 1
      ELSE
         IDAX = IDAX1
         KINTSY = INTSYM
      END IF
      CALL RSP2CR(IGRSYM,IGRSPI,ETRS,VECC,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *            XINDX,OVLAP,DEN1,DEN2,FI,WRK(KWOIT),LWOIT,KZYVA,
     *            LCON,LORB,LREFST,IDAX,KINTSY,ISPIN1,ISPIN2,ISPIN3,
     *            IKLVL,ZYMB,ISYMB,DUM,IDUM,DUM,IDUM,CMO,FC,MJWOP)
C
      IF (ATEST) THEN
         E3NWZC = -DDOT(KZCA,VECA(KAZC),1,ETRS(KAZC),1)
         E3NWZO = -DDOT(KZOA,VECA(KAZO),1,ETRS(KAZO),1)
         E3NWYC = -DDOT(KZCA,VECA(KAYC),1,ETRS(KAYC),1)
         E3NWYO = -DDOT(KZOA,VECA(KAYO),1,ETRS(KAYO),1)
         E3NEW  = E3NWZC + E3NWZO + E3NWYC + E3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 1    :',E3NEW-E3TERM
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 1 ZC :',E3NWZC-E3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 1 ZO :',E3NWZO-E3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 1 YC :',E3NWYC-E3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 1 YO :',E3NWYO-E3TMYO
         E3TMZC = E3NWZC
         E3TMZO = E3NWZO
         E3TMYC = E3NWYC
         E3TMYO = E3NWYO
         E3TERM = E3NEW
      END IF
      IF (IPRRSP .GT. 100) THEN
         WRITE(LUPRI,'(A)') ' Case 1 for E3 term'
         CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
C
C     Case 2:
C     =======
C2000 IF ((MZWOPT(ISYMC) .EQ. 0).OR.(MZWOPT(ISYMB) .EQ. 0)) THEN
 2000 IF (.NOT.LORBC .OR. .NOT.LORBB) THEN
         IF (.NOT.DIROIT) CALL OITCLO(IDAX1,'DELETE')
         GO TO 3000
      END IF
C
C     Transform the integrals with 2k,1k
C
C     Put the factor of one half present in this term into the
C     ZY matrix used for transforming the integrals
C
      NSIM = 1
      CALL GTZYMT(NSIM,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
C
C This half would dissappear were we to use the commutator
C [1k,[2k,H]] = [2k,[1k,H]] + [[2k,1k],H]
C
      CALL DSCAL(NORBT*NORBT,DH,ZYMC,1)
      JTRLVL= 1
      IF (.NOT.DIROIT) THEN
         CALL RSPOIT(IDAX1,IDAX21,JTRLVL,ISYMC,ZYMC,
     *            WRK(KWOIT),KFROIT,LWOIT,IPROIT)
         CALL OITCLO(IDAX1,'DELETE')
      END IF
      CALL DZERO(FI,N2ORBX)
      CALL OITH1(ISYMC,ZYMC,WRK(KH1X),FI,ISYMB)
C     CALL OITH1(LSYOIT,OITMAT,H1,H1X,IH1SYM)
C
C     We have the density matrices over the reference state already
C
      OVLAP = D1
      ISYMDN = 1
C
C     Make the gradient
C
      ISYMV  = IREFSY
      INTSYM = MULD2H(ISYMB,ISYMC)
      LREFST = .TRUE.
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      VECDUM = D0
C
C     Vector may be dummy if LREFST = .TRUE.
C
C (~~| )e(SbSc,+) + (~|~)e(Sb,Sc)
C
      IKLVL=2
      IF (DIROIT) THEN
         IDAX = IDAX0
         KINTSY = 1
      ELSE
         IDAX = IDAX21
         KINTSY = INTSYM
      END IF
      CALL IPSET(0,0,1,1)
      CALL RSP2CR(IGRSYM,IGRSPI,ETRS,VECDUM,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *              XINDX,OVLAP,UDV,PV,FI,WRK(KWOIT),LWOIT,KZYVA,
     *              LCON,LORB,LREFST,IDAX,KINTSY,ISPINB,ISPINC,ISPIN3,
     *              IKLVL,ZYMB,ISYMB,ZYMC,ISYMC,DUM,IDUM,CMO,FC,MJWOP)
C
      IF (.NOT.DIROIT) CALL OITCLO(IDAX21,'DELETE')
C
      IF (ATEST) THEN
         E3NWZC = -DDOT(KZCA,VECA(KAZC),1,ETRS(KAZC),1)
         E3NWZO = -DDOT(KZOA,VECA(KAZO),1,ETRS(KAZO),1)
         E3NWYC = -DDOT(KZCA,VECA(KAYC),1,ETRS(KAYC),1)
         E3NWYO = -DDOT(KZOA,VECA(KAYO),1,ETRS(KAYO),1)
         E3NEW  = E3NWZC + E3NWZO + E3NWYC + E3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 2    :',E3NEW-E3TERM
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 2 ZC :',E3NWZC-E3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 2 ZO :',E3NWZO-E3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 2 YC :',E3NWYC-E3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 2 YO :',E3NWYO-E3TMYO
         E3TMZC = E3NWZC
         E3TMZO = E3NWZO
         E3TMYC = E3NWYC
         E3TMYO = E3NWYO
         E3TERM = E3NEW
      END IF
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 2 for E3 term'
         CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
C
C     Case 3:
C     =======
C3000  IF ( MZWOPT(ISYMC) .EQ. 0 ) GO TO 5000
3000  IF ( .NOT. LORBC ) GO TO 5000
      IF ( IBEQC .EQ. 1 .AND. .NOT.E3TEST) THEN
         CALL DSCAL(KZYVA,D2,ETRS,1)
         IF (ATEST .OR. IPRRSP .GT. 100) THEN
           WRITE(LUPRI,*)
     &           ' E3TEST Case 3 = Case 1 because vecb .eq. vecc'
           WRITE(LUPRI,*)
     &          ' E3TEST Case 4 = Case 2 because vecb .eq. vecc'
           E3TMZC = D2*E3TMZC
           E3TMZO = D2*E3TMZO
           E3TMYC = D2*E3TMYC
           E3TMYO = D2*E3TMYO
           E3TERM = D2*E3TERM
         END IF
         GO TO 5000
      ELSE IF ( IBEQC .EQ. -1 .AND. .NOT.E3TEST) THEN
         KZVA = MZVAR(ISYMA)
         CALL DAXPY(KZVA,DM1,ETRS(KZVA+1),1,ETRS,1)
         CALL DCOPY(KZVA,ETRS,1,ETRS(KZVA+1),1)
         CALL DSCAL(KZVA,DM1,ETRS(KZVA+1),1)
         IF (ATEST .OR. IPRRSP .GT. 100) THEN
           WRITE(LUPRI,*)
     *     ' E3TEST Case 3 = (Case 1)T because vecb .eq. -(vecc)T'
           WRITE(LUPRI,*)
     *     ' E3TEST Case 4 = (Case 2)T because vecb .eq. -(vecc)T'
           E3TMZC = -DDOT(KZCA,VECA(KAZC),1,ETRS(KAZC),1)
           E3TMZO = -DDOT(KZOA,VECA(KAZO),1,ETRS(KAZO),1)
           E3TMYC = -DDOT(KZCA,VECA(KAYC),1,ETRS(KAYC),1)
           E3TMYO = -DDOT(KZOA,VECA(KAYO),1,ETRS(KAYO),1)
           E3TERM = E3TMZC + E3TMZO + E3TMYC + E3TMYO
         END IF
         GO TO 5000
      END IF
C
C     Transform the integrals with kappa(2)
C
      NSIM = 1
      CALL GTZYMT(NSIM,VECC,KZYVC,ISYMC,ZYMC,MJWOP)
      JTRLVL= 2
      IF (.NOT.DIROIT) THEN
         CALL RSPOIT(IDAX0,IDAX2,JTRLVL,ISYMC,ZYMC,
     *            WRK(KWOIT),KFROIT,LWOIT,IPROIT)
      END IF
      CALL DZERO(WRK(KH1X),N2ORBX)
      CALL OITH1(ISYMC,ZYMC,WRK(KH1),WRK(KH1X),1)
C     CALL OITH1(LSYOIT,OITMAT,H1,H1X,IH1SYM)
      CALL DCOPY(N2ORBX,WRK(KH1X),1,FI,1)
C
C     IF (MZCONF(ISYMB) .EQ. 0) GO TO 4000
      IF (.NOT.LCONB) GO TO 4000
C
C     Construct the density matrix <01L|..|0> + <0|..|01R>
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(ISYMB,IREFSY)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMB)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMB)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      NSIM = 1
      ISPIN1 = MULSP(IGRSPI,ISPINC)
      ISPIN2 = 0
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VECB,OVLAP,DEN1,DEN2,ISPIN1,ISPIN2,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
      IF (IGRSPI.EQ.1) THEN
         ISPIN1 = ISPINC
         ISPIN2 = IGRSPI
         CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VECB,OVLAP,DEN1,DEN2(1,2),ISPIN1,ISPIN2,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
      END IF
      CALL IPSET(MULSP(IGRSPI,ISPINC),0,ISPINC,IGRSPI)
C
C     Make the gradient
C
      ISYMV  = ISYMB
      INTSYM = ISYMC
      LREFST = .FALSE.
      NZYVEC = MZYVAR(ISYMB)
      NZCVEC = MZCONF(ISYMB)
      IKLVL = 1
      IF (DIROIT) THEN
         ISPIN1 = ISPINC
         ISPIN2 = 0
         IDAX = IDAX0
         KINTSY = 1
      ELSE
         IDAX = IDAX2
         KINTSY = INTSYM
      END IF
      CALL RSP2CR(IGRSYM,IGRSPI,ETRS,VECB,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *            XINDX,OVLAP,DEN1,DEN2,FI,WRK(KWOIT),LWOIT,KZYVA,
     *            LCON,LORB,LREFST,IDAX,KINTSY,ISPIN1,ISPIN2,ISPIN3,
     *            IKLVL,ZYMC,ISYMC,DUM,IDUM,DUM,IDUM,CMO,FC,MJWOP)
C
      IF (ATEST) THEN
         E3NWZC = -DDOT(KZCA,VECA(KAZC),1,ETRS(KAZC),1)
         E3NWZO = -DDOT(KZOA,VECA(KAZO),1,ETRS(KAZO),1)
         E3NWYC = -DDOT(KZCA,VECA(KAYC),1,ETRS(KAYC),1)
         E3NWYO = -DDOT(KZOA,VECA(KAYO),1,ETRS(KAYO),1)
         E3NEW  = E3NWZC + E3NWZO + E3NWYC + E3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 3    :',E3NEW-E3TERM
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 3 ZC :',E3NWZC-E3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 3 ZO :',E3NWZO-E3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 3 YC :',E3NWYC-E3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 3 YO :',E3NWYO-E3TMYO
         E3TMZC = E3NWZC
         E3TMZO = E3NWZO
         E3TMYC = E3NWYC
         E3TMYO = E3NWYO
         E3TERM = E3NEW
      END IF
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 3 for E3 term'
         CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
C
C      Case 4:
C      =======
C4000 IF ((MZWOPT(ISYMB) .EQ. 0).OR.(MZWOPT(ISYMC) .EQ. 0 )) THEN
 4000 IF (.NOT.LORBB .OR. .NOT.LORBC) THEN
         IF (.NOT.DIROIT) CALL OITCLO(IDAX2,'DELETE')
         GO TO 5000
      END IF
C
C
C     Put the factor of one half present in this term into the
C     ZY matrix used for transforming the integrals
C
      NSIM = 1
      CALL GTZYMT(NSIM,VECB,KZYVB,ISYMB,ZYMB,MJWOP)
      CALL DSCAL(NORBT*NORBT,DH,ZYMB,1)
      JTRLVL= 1
      IF (.NOT.DIROIT) THEN
         CALL RSPOIT(IDAX2,IDAX12,JTRLVL,ISYMB,ZYMB,
     *            WRK(KWOIT),KFROIT,LWOIT,IPROIT)
         CALL OITCLO(IDAX2,'DELETE')
      END IF
      CALL DZERO(FI,N2ORBX)
      CALL OITH1(ISYMB,ZYMB,WRK(KH1X),FI,ISYMC)
C     CALL OITH1(LSYOIT,OITMAT,H1,H1X,IH1SYM)
C
C     We have the density matrices over the reference state already
C
      OVLAP  = D1
      ISYMDN = 1
C
C     Make the gradient
C
      ISYMV  = IREFSY
      INTSYM = MULD2H(ISYMB,ISYMC)
      LREFST = .TRUE.
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      VECDUM = D0
C
C (~~| )e(ScSb,+) + (~|~)e(Sc,Sb)
C
      IKLVL=2
      IF (DIROIT) THEN
         IDAX = IDAX0
         KINTSY = 1
      ELSE
         IDAX = IDAX12
         KINTSY = INTSYM
      END IF
      CALL IPSET(0,0,1,1)
      CALL RSP2CR(IGRSYM,IGRSPI,ETRS,VECDUM,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *            XINDX,OVLAP,UDV,PV,FI,WRK(KWOIT),LWOIT,KZYVA,
     *            LCON,LORB,LREFST,IDAX,KINTSY,ISPINC,ISPINB,ISPIN3,
     *            IKLVL,ZYMC,ISYMC,ZYMB,ISYMB,DUM,IDUM,CMO,FC,MJWOP)
C
      IF (.NOT.DIROIT) CALL OITCLO(IDAX12,'DELETE')
C
      IF (ATEST) THEN
         E3NWZC = -DDOT(KZCA,VECA(KAZC),1,ETRS(KAZC),1)
         E3NWZO = -DDOT(KZOA,VECA(KAZO),1,ETRS(KAZO),1)
         E3NWYC = -DDOT(KZCA,VECA(KAYC),1,ETRS(KAYC),1)
         E3NWYO = -DDOT(KZOA,VECA(KAYO),1,ETRS(KAYO),1)
         E3NEW  = E3NWZC + E3NWZO + E3NWYC + E3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 4    :',E3NEW-E3TERM
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 4 ZC :',E3NWZC-E3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 4 ZO :',E3NWZO-E3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 4 YC :',E3NWYC-E3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 4 YO :',E3NWYO-E3TMYO
         E3TMZC = E3NWZC
         E3TMZO = E3NWZO
         E3TMYC = E3NWYC
         E3TMYO = E3NWYO
         E3TERM = E3NEW
      END IF
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 4 for E3 term'
         CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
C
C     Case 5:
C     =======
C5000  IF (MZCONF(ISYMB) .EQ. 0 .OR. MZCONF(ISYMC) .EQ. 0) GO TO 6000
5000  IF (.NOT.LCONB .OR. .NOT.LCONC) GO TO 6000
C
      CALL DCOPY(N2ORBX,WRK(KH1),1,FI,1)
C
C     Construct <01L|..|02R> + <02L|..|01R>
C
      ILSYM  = MULD2H(IREFSY,ISYMB)
      IRSYM  = MULD2H(IREFSY,ISYMC)
      NCL    = MZCONF(ISYMB)
      NCR    = MZCONF(ISYMC)
      KZVARL = MZYVAR(ISYMB)
      KZVARR = MZYVAR(ISYMC)
      LREF   = .FALSE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      NSIM = 1
      ISPIN1 = IGRSPI
      ISPIN2 = 0
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         VECB,VECC,OVLAP,DEN1,DEN2,ISPIN1,ISPIN2,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
      IF (IGRSPI.EQ.1) THEN
         ISPIN1 = 0
         ISPIN2 = IGRSPI
         CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         VECB,VECC,OVLAP,DEN1,DEN2(1,2),ISPIN1,ISPIN2,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
      END IF
      CALL IPSET(IGRSPI,0,0,IGRSPI)
C
C     Make the gradient
C
      KFREE  = 1
      LFREE  = LWRK  - KFREE
C
      LCON   = .FALSE.
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREFST = .FALSE.
      IDUM   = 1
      VECDUM = D0
C     Calls for configuration part of lt are dummies:
      INTSYM = 1
      IKLVL = 0
      ISPIN1 = 0
      ISPIN2 = 0
      CALL RSP2CR(IGRSYM,IGRSPI,ETRS,VECDUM,IDUM,IDUM,IDUM,ISYMDN,
     *            XINDX,OVLAP,DEN1,DEN2,FI,WRK(KFREE),LFREE,KZYVA,
     *            LCON,LORB,LREFST,IDAX0,INTSYM,ISPIN1,ISPIN2,ISPIN3,
     *            IKLVL,DUM,IDUM,DUM,IDUM,DUM,IDUM,CMO,FC,MJWOP)
C
      IF (ATEST) THEN
         E3NWZC = -DDOT(KZCA,VECA(KAZC),1,ETRS(KAZC),1)
         E3NWZO = -DDOT(KZOA,VECA(KAZO),1,ETRS(KAZO),1)
         E3NWYC = -DDOT(KZCA,VECA(KAYC),1,ETRS(KAYC),1)
         E3NWYO = -DDOT(KZOA,VECA(KAYO),1,ETRS(KAYO),1)
         E3NEW  = E3NWZC + E3NWZO + E3NWYC + E3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 5    :',E3NEW-E3TERM
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 5 ZC :',E3NWZC-E3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 5 ZO :',E3NWZO-E3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 5 YC :',E3NWYC-E3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Case 5 YO :',E3NWYO-E3TMYO
         E3TMZC = E3NWZC
         E3TMZO = E3NWZO
         E3TMYC = E3NWYC
         E3TMYO = E3NWYO
         E3TERM = E3NEW
      END IF
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 5 for E3 term added'
         CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
C
C     End of transformation
C
 6000 CONTINUE
      IF (ATEST) THEN
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Final result    :',E3TERM
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Final result ZC :',E3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Final result ZO :',E3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Final result YC :',E3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' E3TEST Final result YO :',E3TMYO
      END IF
C
      CALL QEXIT('E3DRV')
      RETURN
      END
C  /* Deck rsp2cr */
      SUBROUTINE RSP2CR(IGRSYM,IGRSPI,
     *                  ETRS,VEC,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *                  XINDX,OVLAP,DEN1,DEN2,FI,WRK,LWRK,KZYVA,
     *                  LCON,LORB,LREFST,IDAX,INTSYM,
     *                  ISPIN1,ISPIN2,ISPIN3,
     *                  IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,ZYMAT3,ISYM3,
     *                  CMO,FC,MJWOP)
C
C     Layout the core for RSP2GR
C
#include "implicit.h"
C
#include "infdim.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "qrinf.h"
#include "infpri.h"
#include "infspi.h"
C
      DIMENSION WRK(*)
      DIMENSION ETRS(KZYVA)
      DIMENSION XINDX(*)
      DIMENSION DEN1(*), DEN2(N2ASHX*N2ASHX,*)
      DIMENSION FI(*)
      DIMENSION VEC(NZYVEC),CMO(*),FC(*)
      DIMENSION ZYMAT1(*), ZYMAT2(*), ZYMAT3(*)
C
      LOGICAL LCON, LORB, LREFST, HSORUN
C
      CALL QENTER('RSP2CR')
C
C     Layout of WRK:
C     We keep the H2AX  array at the beginning, in order to
C     release extra workspace in the gradient routine after processing
C     the orbital part of the linear transformation.
C
      KH2AX  = 1
      IF (LCON) THEN
         KFA    = KH2AX  + N2ASHX * N2ASHX
         IF (DIROIT) THEN
            IF (IKLVL.EQ.2) KFA    = KH2AX  + N2ASHX * N2ASHX * 2
            IF (IKLVL.EQ.3) KFA    = KH2AX  + N2ASHX * N2ASHX * 4
         END IF
      ELSE
         KFA    = KH2AX
      END IF
      KQA    = KFA    + NORBT  * NORBT
      KQB    = KQA    + NORBT  * NASHDI
      KPVD   = KQB    + NORBT  * NASHDI
      KH2    = KPVD   + N2ASHX * N2ASHX
      KH2X   = KH2    + N2ORBX
      KWRKO  = KH2X   + N2ORBX
C
C     Get free workspace for orbital part of linear transformation
C
      LWRKO  = LWRK   - KWRKO
      IF (LWRKO.LT.0) CALL ERRWRK('RSP2CR 1',KWRKO-1,LWRK)
C
C     Get space that can be used in configuration part
C     We need the arrays H2XAC and FI, so keep them intact
C
      KWRKC  = KFA
      LWRKC  = LWRK   - KWRKC
      IF (LWRKC.LT.0) CALL ERRWRK('RSP2CR 2',KWRKC-1,LWRK)
      HSORUN = .FALSE.
      CALL RSP2GR (IGRSYM,IGRSPI,ETRS,VEC,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *             FI,WRK(KFA),WRK(KQA),WRK(KQB),WRK(KH2),WRK(KH2X),
     *             OVLAP,DEN1,DEN2,WRK(KPVD),WRK(KH2AX),XINDX,
     *             WRK(KWRKO),LWRKO,WRK(KWRKC),LWRKC,KZYVA,
     *             LCON,LORB,LREFST,IDAX,INTSYM,ISPIN1,ISPIN2,ISPIN3,
     *             IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,ZYMAT3,ISYM3,
     *             HSORUN,CMO,FC,MJWOP)
C
      CALL QEXIT('RSP2CR')
      RETURN
      END
C  /* Deck rsp2gr */
      SUBROUTINE RSP2GR(IGRSYM,IGRSPI,ETRS,VEC,NZYVEC,NZCVEC,ISYMV,
     *                  ISYMDN,
     *                  FI,FA,QA,QB,H2,H2X,
     *                  OVLAP,DEN1,DEN2,PVDEN,H2AX,XINDX,
     *                  WRKO,LWRKO,WRKC,LWRKC,KZYVA,
     *                  LCON,LORB,LREFST,IDAX,INTSYM,
     *                  ISPIN1,ISPIN2,ISPIN3,
     *                  IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,ZYMAT3,ISYM3,
     *                  HSORUN,CMO,FC,MJWOP)
C
C     This routine computes the gradient, given a density matrix,
C     the two vectors and a set of (transformed) integrals. In the deriva-
C     vation of the implemented formulas it is assumed that the integrals
C     have particle permutation symmetry.
C
C     The result is returned added to the values present in ETRS upon
C     calling.
C
C     Memory layout:
C     1.  WRKO is the workarray for the orbital part of the lt.
C     2.  WRKC is the workarray for the CI part of the lt. Since we do
C         not need the active fock matrix and the Q matrices here,
C         some space can be released. Organisation should be
C         done by the calling routine. (See E3NBNC for example)
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, DH = 0.50D0 )
C
C INFINP : DIRFCK
C
#include "dftcom.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infvar.h"
#include "qrinf.h"
#include "inforb.h"
#include "infind.h"
#include "infdim.h"
#include "infpri.h"
#include "infspi.h"
#include "infinp.h"
C
      DIMENSION ETRS(KZYVA), MJWOP(2,MAXWOP,8)
      DIMENSION FI(NORBT,NORBT),FA(NORBT,NORBT)
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION H2 (NORBT,NORBT)
      DIMENSION H2X (NORBT,NORBT)
      DIMENSION H2AX(NASHDI,NASHDI,NASHDI,NASHDI,*)
      DIMENSION DEN1(NASHDI,NASHDI), DEN2(NASHDI,NASHDI,NASHDI,NASHDI,*)
      DIMENSION PVDEN(NASHDI,NASHDI)
      DIMENSION VEC(NZYVEC),CMO(*),FC(*)
      DIMENSION XINDX(*)
      DIMENSION WRKO(*), WRKC(*)
C
      LOGICAL LCON, LORB, LREFST, HSORUN, DFTSAV, ADDFOCK
C
      CALL QENTER('RSP2GR')
C
      IF (DOHFSRDFT .OR. DOMCSRDFT) THEN
         CALL QUIT('srDFT not implemented in RSP2GR')
      END IF
C
      CALL DZERO(FA,NORBT*NORBT)
      IF (NASHT .GT. 0) THEN
         CALL DZERO(QA,NASHT*NORBT)
         CALL DZERO(QB,NASHT*NORBT)
         IF (LCON) THEN
            LH2AX = N2ASHX*N2ASHX
            IF (DIROIT) THEN
               IF (IKLVL.EQ.2) LH2AX = LH2AX*2
               IF (IKLVL.EQ.3) LH2AX = LH2AX*4
            END IF
            CALL DZERO(H2AX,LH2AX)
         END IF
      END IF
C
C
      ONEFAC  = D0
      IF (ISPIN1.EQ.ISPIN2) THEN
         DO 105 IW = 1, NISHT
            IX = ISX(IW)
            ONEFAC = ONEFAC + FI(IX,IX)
105      CONTINUE
      END IF
C
      IF( IPRRSP .GT. 200 ) THEN
         WRITE(LUPRI,'(/A,F15.8)') ' One electron factor : ',ONEFAC
         WRITE(LUPRI,'(/A)') ' One electron matrix in RSP2GR'
         WRITE(LUPRI,'(A)')  ' ============================='
         CALL OUTPUT(FI,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C     Construct FI, FA, QA, QB, H2AX matrices
C
      KFREE = 1
      LFREE = LWRKO
      IF (DIROIT) THEN
         IF (HSORUN) THEN
            CALL HSOFXD (FI,FA,QA,QB,H2AX,
     *            IDAX,INTSYM, ISYMDN,DEN1,DEN2,
     *            PVDEN,H2,H2X,WRKO,KFREE,LFREE,
     *            LCON,LORB,IPRRSP,IGRSPI,ISPIN1,ISPIN2,
     *            IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2)
         ELSE IF (DIRFCK .AND. NASHT.EQ.0) THEN
            DFTSAV = DFTADD
            DFTADD = .FALSE.
c
c     note that QFOCK returns complete FOCK/Kohn-Sham matrix, including
c     DFT contribution, as opposed to RSPFXD/RSPFX which return electronic
c     part only (w/o DFT contribution). If this is the case (i.e. .NOT.DIRFCK) 
c     we add the DFT contribution to KS matrix later in DFTQRC.
c
            IF (NASHT.GT.0) THEN
               !This point in code is reached for the special case
               !single open shell direct HF
               WRITE(LUPRI,*)'Input combination not implemented'
               WRITE(LUPRI,*)
     &        'Remove .DIRECT or reformulate problem as DFT GGAKey HF=1'
               CALL QUIT('QFOCK not implemented for open shells')
            END IF
            CALL QFOCK(1,ISYM1,ISYM2,IGRSPI,ISPIN1,ISPIN2,
     &                 ZYMAT1,ZYMAT2,FC,CMO,FI,WRKO,KFREE,LFREE)
            DFTADD = DFTSAV
         ELSE
            CALL RSPFXD (FI,FA,QA,QB,H2AX,
     *            IDAX,INTSYM, ISYMDN,DEN1,DEN2,
     *            PVDEN,H2,H2X,WRKO,KFREE,LFREE,
     *            LCON,LORB,IPRRSP,IGRSPI,ISPIN1,ISPIN2,ISPIN3,
     *            IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,ZYMAT3,ISYM3)
         END IF
      ELSE
         CALL RSPFX (FI,FA,QA,QB,H2AX,
     *            IDAX,INTSYM, ISYMDN,DEN1,DEN2,
     *            PVDEN,H2X,WRKO,KFREE,LFREE,
     *            LCON,LORB,IPRRSP)
      END IF
C
      IF( IPRRSP .GT. 200 ) THEN
         WRITE(LUPRI,'(/A)') ' COMPLETED MATRICES IN  RSP2GR'
         WRITE(LUPRI,'(A)')  ' ============================='
         WRITE(LUPRI,'(///A)') ' Inactive fock matrix'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(///A)') ' Active fock matrix'
         CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(///A)') ' Qa matrix'
         CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(///A)') ' Qb matrix'
         CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
         CALL PRIAC2(H2AX,NASHT,LUPRI)
      END IF
C
C     Now distribute the computed matrices over the orbital part
C     of the gradient:
C
      IF ( LORB ) THEN
        TRPLET = IGRSPI.NE.MULSP(ISPIN1,ISPIN2)
        CALL ORBEX(IGRSYM,ISYMDN,KZYVA,FI,FA,QA,QB,ETRS,OVLAP,DEN1,
     &             MJWOP)
      END IF
C
C     and over the configuration part:
C
      IF ( LCON ) THEN
         ISGN1 = (-1)**ISPIN1
         ISGN2 = (-1)**ISPIN2
         ISYMJ  = MULD2H( IGRSYM, IREFSY )
         NZCSTJ = MZCONF(IGRSYM)
         IF  ( LREFST ) THEN
            LFREE  = LWRKC - MZCONF(1)
            IF (LFREE.LT.0) CALL ERRWRK('RSP2GR',MZCONF(1),LWRKC)
            CALL GETREF(WRKC,MZCONF(1))
            CALL CONEX(KZYVA,IGRSYM,FI,ONEFAC,H2AX,WRKC,MZCONF(1),
     *                 MZCONF(1),IREFSY,NZCSTJ,ISYMJ,LREFST,ETRS,
     *                 XINDX,ISPIN1,ISPIN2,WRKC(NCREF+1),LFREE)
         ELSE
            CALL CONEX(KZYVA,IGRSYM,FI,ONEFAC,H2AX,VEC,
     *                 NZYVEC,NZCVEC,
     *                 ISYMV, NZCSTJ,ISYMJ,LREFST,ETRS,XINDX,
     *                 ISPIN1,ISPIN2,WRKC,LWRKC)
         END IF
      END IF
C
C     and that's it.
C
      CALL QEXIT('RSP2GR')
      RETURN
      END
C  /* Deck rspfx */
      SUBROUTINE RSPFX (FI,FA,QA,QB,H2AX,
     *                  IDAX,INTSYM, ISYMDN,DEN1,DEN2,
     *                  PVDEN,H2X,WRK,KFREE,LFREE,
     *                  LCON,LORB,IPRFX)
C
C
C 15-Nov-1991 hjaaj (pulled out of RSP2GR)
C
C     This routine computes the FX matrices, that is, the Fock
C     type matrices FI, FA, QA, and QB with transformed ("X")
C     integrals.  In addition, If LCON true then the H2AX matrix
C     with transformed integrals is extracted for the CI routines.
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, DH = 0.50D0 )
C
C Used from common blocks:
C  ?
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "infdim.h"
#include "infpri.h"
C
C
      DIMENSION FI(NORBT,NORBT),FA(NORBT,NORBT)
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION H2X (NORBT,NORBT)
      DIMENSION H2AX(NASHDI,NASHDI,NASHDI,NASHDI)
      DIMENSION DEN1(NASHDI,NASHDI), DEN2(NASHDI,NASHDI,NASHDI,NASHDI)
      DIMENSION PVDEN(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION NEEDED(-4:6)
C
      LOGICAL LCON, LORB
C
      CALL QENTER('RSPFX')
C
C     Put up the structure for reading the (transformed) integrals:
C     IEND < 0 flags the completeness of the distributions
C
      NEEDED(-4:6) = 0
      NEEDED( 1:3) = 1
      IFCKSY = MULD2H(INTSYM,ISYMDN)
      IDIST = 0
1000  CALL NXTH2X(IP,IQ,H2X,IDAX,NEEDED,WRK,KFREE,LFREE,IDIST)
      IF ( IDIST .LT. 0) GO TO 2000
C
      IPTYP = IOBTYP(IP)
      IQTYP = IOBTYP(IQ)
C
      IF (IPRFX .GT. 200) THEN
         WRITE(LUPRI,'(//A,/A,/A,I5,/A,I5,/A)')
     *   '     Getting from NXTH2X: ',
     *   '     ==================== ',
     *   '     IP    = ',IP,
     *   '     IQ    = ',IQ,
     *   '     ===================='
      END IF
C
C     Determine type of distribution
C
C     1 : inactive   - inactive
C     2 : active     - inactive
C     3 : inactive   - active
C     4 : active     - active
C
      IPQTYP = 0
      IF (( IPTYP .EQ. JTINAC) .AND. ( IQTYP .EQ. JTINAC)) IPQTYP = 1
      IF (( IPTYP .EQ. JTACT ) .AND. ( IQTYP .EQ. JTINAC)) IPQTYP = 2
      IF (( IPTYP .EQ. JTINAC) .AND. ( IQTYP .EQ. JTACT )) IPQTYP = 3
      IF (( IPTYP .EQ. JTACT ) .AND. ( IQTYP .EQ. JTACT )) IPQTYP = 4
C
      IF (IPQTYP .EQ. 0) THEN
         WRITE(LUPRI,'(//A,/A,I5,A,I5,/A)')
     *   ' --> Warning : Unidentifyable distribution found ; ',
     *   '     IPTYP = ', IPTYP,' IQTYP = ',IQTYP,
     *   '     Continuing execution by skipping this load'
         GO TO 1000
      END IF
C
C     Generally, an integral is labeled by (pq|rs). We know that p and q
C     are occupied. When going into the specific routines, we substitute
C     p,q,r,s for the special cases. We get symmetry characteristics:
C
      IPSYM = ISMO( IP )
      IQSYM = ISMO( IQ )
      IPQSYM = MULD2H(IPSYM,IQSYM)
C
      IF ( IPRFX .GT. 200) THEN
         WRITE(LUPRI,'(//A)') 'Characteristics of distribution:'
         WRITE(LUPRI,'(A)')   '================================'
         WRITE(LUPRI,'(A,I5)') 'IPTYP   = ', IPTYP
         WRITE(LUPRI,'(A,I5)') 'IQTYP   = ', IQTYP
         WRITE(LUPRI,'(A,I5)') 'IPQTYP  = ', IPQTYP
         WRITE(LUPRI,'(A,I5)') 'IPSYM   = ', IPSYM
         WRITE(LUPRI,'(A,I5)') 'IQSYM   = ', IQSYM
         WRITE(LUPRI,'(A,I5)') 'IPQSYM  = ', IPQSYM
         WRITE(LUPRI,'(A)')   '================================'
C
         WRITE(LUPRI,'(//A/A/)')
     *        ' Two-electron integrals from this distribution',
     *        ' ============================================='
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IF (IPQTYP .LE. 3) THEN
         CALL FIXADD(INTSYM,IP,IQ,IPQTYP,IPQSYM,IPSYM,IQSYM,
     *        FI,H2X )
      END IF
      IF ((IPQTYP .NE. 1) .AND. LORB) THEN
         CALL FAXADD(IFCKSY,INTSYM,IP,IQ,IPQTYP,IPQSYM,IPSYM,
     *        IQSYM,FA,H2X,DEN1 )
      END IF
      IF ( (IPQTYP .EQ. 4) ) THEN
         ISWIP = ISW( IP ) - NISHT
         ISWIQ = ISW( IQ ) - NISHT
         IF (LORB) THEN
            CALL QXADD(IFCKSY,INTSYM,ISWIP,ISWIQ,IP,IQ,IPQTYP,
     *           IPQSYM,IPSYM,IQSYM,
     *           QA,QB,H2X,DEN2,PVDEN)
         END IF
C     
C     Store (xy|uv) in H2ACX for use in configuration part of lt.
C     
         IF ( LCON ) THEN
            IF ( IPRFX .GT. 200) THEN
               WRITE(LUPRI,'(//A,/A,2(/A,I5),/A)' )
     *              'Storing in H2AX :',
     *              '=================',
     *              'ISWIP = ',ISWIP,
     *              'ISWIQ = ',ISWIQ,
     *              '================='
            END IF
            IUVSYM = MULD2H(IPQSYM, INTSYM)
            CALL DSH2AX(H2AX(1,1,ISWIP,ISWIQ),H2X,IUVSYM)
         END IF
C     
         IF( IPRFX .GT. 200 ) THEN
            WRITE(LUPRI,'(/A)') ' PARTIAL MATRICES IN  RSP2GR'
            WRITE(LUPRI,'(A)')  ' ==========================='
            WRITE(LUPRI,'(//A)') ' Inactive fock matrix'
            CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Active fock matrix'
            CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Qa matrix'
            CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Qb matrix'
            CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         END IF
C
      END IF
C
C     We have now completed processing this load, so get the next
C
      GO TO 1000
C
2000  CONTINUE
C
C     Account for the fact that we calculated 2 times FA
      CALL DSCAL(NORBT*NORBT,DH,FA,1)
C
      CALL QEXIT('RSPFX')
      RETURN
      END
C  /* Deck fixadd */
      SUBROUTINE FIXADD(INTSYM,IP,IQ,IPQTYP,IPQSYM,IPSYM,IQSYM,
     *                 FI,H2X )
C
C     Add contribution presently in H2X  to the inactive Fock matrix.
C     Generic indices are i1,i2, for the rest we adhere to convention.
C     Expression:
C
C     iF(i1,i2) = iF(i1,i2) + sum(j) [ 2 (jj|i1 i2) - (i1 j|j i2) ]
C
C     iF(i1,i2) = iF(i1,i2) + sum(j) [ 2 (jj|i1 i2) - (j i2|i1 j) ]
C
C     Procedure:
C     ----------
C     1) Check if for the inactive - inactive distributions we have
C        i = j. If so, add the direct contribution to the inactive Fock
C        matrix.
C
C     2) For occupied - inactive distributions we have to subtract
C        (oj|js) from iF(o,i2) for all s in the distribution
C
C     3) For secondary - inactive distributions we have to subtract
C        (jo|aj) from iF(i1,o)
C
C     Needed inactive Fock matrices and corresponding contributions:
C     (First column exchange = case 2, second = case 3)
C----------------------------------------------------------------------
C     Label  direct     exchange    Label  direct     exchange
C----------------------------------------------------------------------
C     ia     (jj|ia)    (ij|ja)     ai     (jj|ai)    (aj|ji) = (ji|aj)
C     ti     (jj|ti)    (tj|ji)     at     (jj|at)    (aj|jt) = (jt|aj)
C     it     (jj|it)    (ij|jt)
C     ta     (jj|ta)    (tj|ja)
C     tu     (jj|tu)    (tj|ju)
C----------------------------------------------------------------------
C
#include "implicit.h"
C
      PARAMETER ( D2= 2.0D0, DM1 = -1.0D0 )
C
#include "priunit.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infhyp.h"
#include "infdim.h"
#include "inforb.h"
#include "infpri.h"
C
      DIMENSION FI(NORBT,NORBT),H2X (NORBT,NORBT)
C
      IF (IPQTYP .EQ. 4) THEN
C     (We have no contributing distributions, but print a warning)
         WRITE(LUPRI,'(/A)')
     *   ' --> Warning: invalid call of FIXADD, IPQTYP = 4'
         RETURN
      END IF
      IF (( IP .EQ. IQ) .AND. (IPQTYP .EQ. 1)) THEN
C     Process case 1)  inactive - inactive distributions coulomb part
C
C        iF(i1,i2) = iF(i1,i2) + 2.0d0 (jj|i1 i2)
C
         CALL DAXPY(NORBT*NORBT,D2,H2X,1,FI,1)
      END IF
C
C     Now process all exchange contributions
C
      IF (( IPQTYP .EQ. 1) .OR. (IPQTYP .EQ. 2) )  THEN
C     Process case 2)    : occupied - inactive distribution
C
C        iF(o,i2) = iF(o,i2) -  (oj|j i2)
C
         I2SYM   = MULD2H(IPSYM,INTSYM)
         NORBI2  = NORB(I2SYM)
         IOFFI2  = IORB(I2SYM) + 1
         IF( NORBI2 .GT. 0) THEN
            CALL DAXPY(NORBI2,DM1,H2X(IQ,IOFFI2),NORBT,
     *                 FI(IP,IOFFI2),NORBT)
         END IF
      END IF
      IF (( IPQTYP .EQ. 1) .OR. (IPQTYP .EQ. 3) )  THEN
C     Process case 3)   :  inactive - occupied distribution
C
C        iF(a,o) = iF(a,o) -  (jo|aj)
C
         IASYM  = MULD2H(IQSYM,INTSYM)
         NSSHA = NSSH(IASYM)
         IOFFA = IORB(IASYM) + 1 + NOCC(IASYM)
         IF(NSSHA .GT. 0) THEN
            CALL DAXPY(NSSHA,DM1,H2X (IOFFA,IP),1,
     *                 FI(IOFFA,IQ),1)
         END IF
      END IF
      RETURN
      END
C  /* Deck faxadd */
      SUBROUTINE FAXADD(IFCKSY,INTSYM,IP,IQ,IPQTYP,
     *                  IPQSYM,IPSYM,IQSYM,FA,H2X,DEN1)
C
C     Add contribution presently in H2X  to the active Fock matrix.
C     Note: twice the active Fock matrix is computed.
C
C     aF(i1,i2) = aF(i1,i2) + [ (xy|i1 i2) - (i1 y|x i2) ] D(x,y)
C
C     aF(i1,i2) = aF(i1,i2) + [ (xy|i1 i2) - (x i2|i1 y) ] D(x,y)
C
C     Procedure:
C     ----------
C     1) Check for the active - active distributions and add the
C        direct contribution to the inactive Fock matrix.
C
C     2) For occupied - active distributions we have to subtract
C        (oy|xs) D(x,y) from aF(o,i2) for all s in the distribution
C
C     3) For active - occupied distributions we have to subtract
C        (xo|ay) D(x,y) from aF(i1,o)
C
C     Needed active Fock matrices and corresponding contributions:
C     (First column exchange = case 2, second = case 3)
C----------------------------------------------------------------------
C     Label  direct     exchange    Label  direct     exchange
C----------------------------------------------------------------------
C     ia     (xy|ia)    (iy|xa)     ai     (xy|ai)    (ay|xi) = (xi|ay)
C     ti     (xy|ti)    (ty|xi)
C     it     (xy|it)    (iy|xt)
C----------------------------------------------------------------------
C
#include "implicit.h"
C
      PARAMETER ( D2= 2.0D0 )
C
#include "wrkrsp.h"
#include "rspprp.h"
#include "infhyp.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infdim.h"
#include "infind.h"
#include "inforb.h"
#include "infpri.h"
C
      DIMENSION FA(NORBT,NORBT),H2X (NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
C
      IF (IPQTYP .EQ. 1) THEN
C     (We have no contributing distributions, but print a warning)
         WRITE(LUPRI,'(/A)')
     *   ' --> Warning: invalid call of FAXADD, IPQTYP = 1'
         RETURN
      END IF
      IF ( IPQTYP .EQ. 4 ) THEN
C     Process case 1)  active - active distributions coulomb part
C
C        aF(i1,i2) = aF(i1,i2) + (xy|i1 i2) D(x,y)
C
         IXDEN = ISW(IP) - NISHT
         IYDEN = ISW(IQ) - NISHT
         FACT = D2 * DEN1(IXDEN,IYDEN)
         CALL DAXPY(NORBT*NORBT,FACT,H2X,1,FA,1)
      END IF
C
C     Now process all exchange contributions
C
      IF (( IPQTYP .EQ. 3) .OR. (IPQTYP .EQ. 4) )  THEN
C     Process case 2)    : occupied - active distribution
C
C        aF(o,i2) = aF(o,i2) -  sum(x) D(x,y) (oy|x i2)
C
         I2SYM   = MULD2H(IPSYM,IFCKSY)
         IXSYM   = MULD2H( MULD2H( IPQSYM, I2SYM), INTSYM)
         NASHX   = NASH(IXSYM)
         NORBI2  = NORB(I2SYM)
         IF ( (NORBI2.GT.0) .AND. (NASHX.GT.0) ) THEN
            IOFFI2  = IORB(I2SYM) + 1
            IOFFX   = IORB(IXSYM) + 1 + NISH(IXSYM)
            IDENX   = ISW( IOFFX ) - NISHT
            IDENY   = ISW( IQ    ) - NISHT
            CALL DGEMM('T','N',1,NORBI2,NASHX,-1.D0,
     &                 DEN1(IDENX,IDENY),NASHT,
     &                 H2X(IOFFX,IOFFI2),NORBT,1.D0,
     &                 FA(IP,IOFFI2),NORBT)
         END IF
      END IF
      IF ( ( IPQTYP .EQ. 2) .OR. ( IPQTYP .EQ. 4) ) THEN
C     Process case 3)   :  active - occupied distribution
C
C        aF(a,o) = aF(a,o) -  sum(y) (xo|ay) D(x,y)
C
         IASYM   = MULD2H(IQSYM,IFCKSY)
         IYSYM   = MULD2H( MULD2H( IPQSYM, IASYM), INTSYM)
         NASHY   = NASH(IYSYM)
         NSSHA   = NSSH(IASYM)
         IF ( (NSSHA.GT.0) .AND. (NASHY.GT.0) ) THEN
            IOFFA   = IORB(IASYM) + 1 + NOCC(IASYM)
            IOFFY   = IORB(IYSYM) + 1 + NISH(IYSYM)
            IDENX   = ISW( IP    ) - NISHT
            IDENY   = ISW( IOFFY ) - NISHT
            CALL DGEMM('N','T',NSSHA,1,NASHY,-1.D0,
     &                 H2X(IOFFA,IOFFY),NORBT,
     &                 DEN1(IDENX,IDENY),NASHT,1.D0,
     &                 FA(IOFFA,IQ),NORBT)
         END IF
      END IF
      RETURN
      END
C  /* Deck qxadd */
      SUBROUTINE QXADD(IFCKSY,INTSYM,ISWIP,ISWIQ,IP,IQ,IPQTYP,
     *                IPQSYM,IPSYM,IQSYM,QA,QB,
     *                H2X,DEN2,PVDEN)
C
C     Process contributions from the Q matrices. Formulas:
C
C     aQ(i1,y) = aQ(i1,y) + sum(w) (xv|w i1) d(xv;wy)
C
C     bQ(i1,y) = bQ(i1,y) + sum(w) (xv|i1 w) d(xv;yw)
C
C     The sum over x and v is taken by reading in all active contribu-
C     tions.
C
#include "implicit.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "inforb.h"
#include "infdim.h"
#include "infind.h"
#include "wrkrsp.h"
#include "infpri.h"
C
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION DEN2(N2ASHX,N2ASHX), PVDEN(NASHDI,NASHDI)
      DIMENSION H2X(NORBT,NORBT)
C
      IF (IPQTYP .NE. 4) THEN
C     (We have no contributing distributions, but print a warning)
         WRITE(LUPRI,'(/A)')
     *   ' --> Warning: invalid call of QXADD, IPQTYP is not 4'
         RETURN
      END IF
C
C
C     First get d(xv;**)
C
      CALL PVXDIS(ISWIP,ISWIQ,PVDEN,DEN2,4)
C
      IXSYM  = IPSYM
      IVSYM  = IQSYM
      DO 200 IWSYM = 1, NSYM
         II1SYM = MULD2H( MULD2H( IPQSYM, IWSYM) , INTSYM)
         IYSYM  = MULD2H(II1SYM,IFCKSY)
         NORBI1 = NORB(II1SYM)
         NASHW  = NASH(IWSYM)
         NASHY  = NASH(IYSYM)
         IF((NORBI1.GT.0) .AND. (NASHW .GT. 0) .AND.
     *                          (NASHY .GT. 0) ) THEN
            IOFFI1 = IORB(II1SYM) + 1
            IOFFW  = IORB(IWSYM) + NISH(IWSYM) + 1
            IOFFY  = IORB(IYSYM) + NISH(IYSYM) + 1
            IDENW  = ISW( IOFFW ) - NISHT
            IDENY  = ISW( IOFFY ) - NISHT
C
C     aQ(i1,y) = aQ(i1,y) + sum(w) (xv|w i1) d(xv;wy)
C
            CALL DGEMM('T','N',NORBI1,NASHY,NASHW,1.D0,
     &                 H2X(IOFFW,IOFFI1),NORBT,
     &                 PVDEN(IDENW,IDENY),NASHT,1.D0,
     &                 QA(IOFFI1,IDENY),NORBT)
C
C     bQ(i1,y) = bQ(i1,y) + sum(w) (xv|i1 w) d(xv;yw)
C
            CALL DGEMM('N','T',NORBI1,NASHY,NASHW,1.D0,
     &                 H2X(IOFFI1,IOFFW),NORBT,
     &                 PVDEN(IDENY,IDENW),NASHT,1.D0,
     &                 QB(IOFFI1,IDENY),NORBT)
C
         END IF
200   CONTINUE
      RETURN
      END
C  /* Deck orbex */
      SUBROUTINE ORBEX(IGRSYM,ISYMDN,KZYVA,
     *                 FI,FA,QA,QB,ETRS,OVLAP,DEN1,MJWOP)
C
C     Pupose:
C     Get the gradient vector from the active ,inactive, Qa and Qb ma-
C     trices. We have the expression:
C
C     F(l,k) = <0| [ E(k,l), K ] |0>
C
C    We have the expressions:
C
C    1) F(a,i) = 2 * OVLAP * iF(a,i) + 2 aF(a,i)
C
C    2) F(i,a) =  - 2 * OVLAP * iF(i,a) - 2 aF(i,a)
C
C    3) F(t,i) = 2 * OVLAP * iF(t,i) + 2 aF(t,i) - sum(x) iF(x,i) D(t,i)
C                - aQ(i,t)
C
C    4) F(i,t) = sum(x) iF(t,x) + bQ(i,t) - 2 * OVLAP * iF(i,t)
C                - 2 aF(i,t)
C
C    5) F(t,a) = -sum(x) iF(x,a) D(x,t) - aQ(a,t)
C
C    6) F(a,t) = sum(x) iF(a,x) D(t,x) + bQ(a,t)
C
C    7) F(u,t) = sum(x) iF(u,x) D(t,x) + bQ(u,t)
C              - sum(x) iF(x,t) D(x,u) - aQ(t,u)
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, DM1 = -1.0D0 )
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "infdim.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "qrinf.h"
#include "infspi.h"
C
      DIMENSION ETRS(KZYVA)
      DIMENSION FI(NORBT,NORBT),FA(NORBT,NORBT)
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION DEN1(NASHDI,NASHDI), MJWOP(2,MAXWOP,8)
C
C
         IF( IPRRSP .GT. 200 ) THEN
            WRITE(LUPRI,'(/A)') ' Matrices in ORBEX'
            WRITE(LUPRI,'(A)')  ' ================='
            WRITE(LUPRI,'(/A)') ' Density matrix in ORBEX'
            CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,-1,LUPRI)
            WRITE(LUPRI,'(/A)') ' Inactive fock matrix in ORBEX'
            CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
            WRITE(LUPRI,'(/A)') ' Active fock matrix in ORBEX'
            CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
            WRITE(LUPRI,'(/A)') ' Qa matrix in ORBEX'
            CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,-1,LUPRI)
            WRITE(LUPRI,'(/A)') ' Qb matrix in ORBEX'
            CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,-1,LUPRI)
         END IF
C
C
C     Compute offset for the F(k,l) elements
C
      NZCONF = MZCONF(IGRSYM)
      NYCONF = MZCONF(IGRSYM) + MZVAR(IGRSYM)
C
C     distribute fock matrices in the gradient:
C
         KSYM1 = 0
         DO 1300 IG = 1,MZWOPT(IGRSYM)
            K     = MJWOP(1,IG,IGRSYM)
            L     = MJWOP(2,IG,IGRSYM)
            KSYM  = ISMO(K)
            LSYM  = ISMO(L)
            IF( KSYM.NE.KSYM1 ) THEN
               KSYM1 = KSYM
               NORBK = NORB(KSYM)
               IORBK = IORB(KSYM)
               NASHK = NASH(KSYM)
               NISHK = NISH(KSYM)
               IASHK = IASH(KSYM)
               KXSYM = MULD2H(KSYM,ISYMDN)
               NORBKX= NORB(KXSYM)
               IORBKX= IORB(KXSYM)
               NASHKX= NASH(KXSYM)
               NISHKX= NISH(KXSYM)
               IASHKX= IASH(KXSYM)
               IORBL = IORB(LSYM)
               NASHL = NASH(LSYM)
               NISHL = NISH(LSYM)
               NORBL = NORB(LSYM)
               IASHL = IASH(LSYM)
               LXSYM = MULD2H(LSYM,ISYMDN)
               NORBLX= NORB(LXSYM)
               IORBLX= IORB(LXSYM)
               NASHLX= NASH(LXSYM)
               NISHLX= NISH(LXSYM)
               IASHLX= IASH(LXSYM)
            END IF
            NK    = K - IORBK
            NL    = L - IORBL
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF ( ITYPK.EQ.JTINAC )THEN
C
C           Do F(a,i) and F(i,a) and a part of F(t,i) and F(i,t)
C
                  IF (.NOT.TRPLET) THEN
                     ETRS(NZCONF+IG) = ETRS(NZCONF+IG)
     *                    + D2 * OVLAP * FI(L,K)
                     ETRS(NYCONF+IG) = ETRS(NYCONF+IG)
     *                    - D2 * OVLAP * FI(K,L)
                  END IF
                  IF ( NASHT . GT . 0 ) THEN
                     ETRS(NZCONF+IG) = ETRS(NZCONF+IG) +  D2 * FA(L,K)
                     ETRS(NYCONF+IG) = ETRS(NYCONF+IG) -  D2 * FA(K,L)
                  END IF
                  IF ( ITYPL.EQ.JTACT ) THEN
C
C                 Complete F(t,i) and F(i,t)
C
                     NWL = ISW(L) - NISHT
                     IF (NASHLX .GT. 0) THEN
                       ETRS(NZCONF+IG) = ETRS(NZCONF+IG)
     *                  - DDOT(NASHLX,FI(IORBLX+NISHLX+1,K),1,
     *                                   DEN1(IASHLX+1,NWL),1)
                       DO IX = 1,NASHLX
                         ETRS(NYCONF+IG) = ETRS(NYCONF+IG)
     *                     + FI(K,IORBLX+NISHLX+IX)*DEN1(NWL,IASHLX+IX)
                       END DO
                     END IF
                     ETRS(NZCONF+IG) = ETRS(NZCONF+IG) - QA(K,NWL)
                     ETRS(NYCONF+IG) = ETRS(NYCONF+IG) + QB(K,NWL)
                  END IF
            ELSE
C
C           We now know that K labels an active orbital
C
               IF (ITYPL.EQ.JTACT) THEN
C
C                 Do F(t,u)
C
                  NWL = ISW(L) - NISHT
                  NWK = ISW(K) - NISHT
                  IF (NASHLX .GT. 0) THEN
                    ETRS(NZCONF+IG) = ETRS(NZCONF+IG)
     *               -DDOT(NASHLX,FI(IORBLX+NISHLX+1,K),1,
     *                           DEN1(IASHLX+1,NWL),1)
                    DO IX = 1,NASHLX
                      ETRS(NYCONF+IG) = ETRS(NYCONF+IG)
     *                  + FI(K,IORBLX+NISHLX+IX)*DEN1(NWL,IASHLX+IX)
                    END DO
                  END IF
                  IF (NASHKX .GT. 0) THEN
                    ETRS(NYCONF+IG) = ETRS(NYCONF+IG)
     *                  -DDOT(NASHKX,FI(IORBKX+NISHKX+1,L),1,
     *                               DEN1(IASHKX+1,NWK),1)
                    DO IX = 1,NASHKX
                      ETRS(NZCONF+IG) = ETRS(NZCONF+IG)
     *                  + FI(L,IORBKX+NISHKX+IX)*DEN1(NWK,IASHKX+IX)
                    END DO
                  END IF
                  ETRS(NZCONF+IG) = ETRS(NZCONF+IG)
     *                  +QB(L,NWK) - QA(K,NWL)
                  ETRS(NYCONF+IG) = ETRS(NYCONF+IG)
     *                  +QB(K,NWL) - QA(L,NWK)
               ELSE
C
C                 Do F(a,t) and F(t,a)
C
                  NWK = ISW(K) - NISHT
                  IF (NASHKX .GT. 0) THEN
                    ETRS(NYCONF+IG) = ETRS(NYCONF+IG)
     *                  - DDOT(NASHKX,FI(IORBKX+NISHKX+1,L),1,
     *                              DEN1(IASHKX+1,NWK),1)
                    DO IX = 1,NASHKX
                     ETRS(NZCONF+IG) = ETRS(NZCONF+IG)
     *                  + FI(L,IORBKX+NISHKX+IX)*DEN1(NWK,IASHKX+IX)
                    END DO
                  END IF
                  ETRS(NYCONF+IG) = ETRS(NYCONF+IG) - QA(L,NWK)
                  ETRS(NZCONF+IG) = ETRS(NZCONF+IG) + QB(L,NWK)
               ENDIF
            ENDIF
 1300    CONTINUE
C
      IF (IPRRSP .GT. 200 ) THEN
         WRITE(LUPRI,'(//A/A)')
     &   ' ORBEX: Orbital part of linear transformation',
     &   ' ============================================'
         CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck dsh2ax */
      SUBROUTINE DSH2AX(H2AX,H2X,IUVSYM)
C
C 891013-hjaaj
C Modified for no symmetry and name changed 901028 hh
C
C Add contributions from (xy) distribution to
C H2AC(uv,xy) from H2XY(uv); (uv) distribution has symmetry IUVSYM
C
C Subroutine can also be used for contributions to H2XAC.
C
#include "implicit.h"
C
C Used from common blocks:
C   INFORB : NSYM, MULD2H, NORBT, NASHT, IASH, NASH
C   INFIND : ISX, NSM, IROW
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infpri.h"
C
      DIMENSION H2AX(NASHT,NASHT),H2X(NORBT,NORBT)
C
      DO 200 NVW = 1,NASHT
         IVSYM = NSM(NVW)
         IUSYM = MULD2H(IVSYM, IUVSYM)
         IV    = ISX(NISHT+NVW)
         NUWST = IASH(IUSYM) + 1
         NUWEND= IASH(IUSYM) + NASH(IUSYM)
         DO 100 NUW = NUWST,NUWEND
            IU   = ISX(NISHT+NUW)
            H2AX(NUW,NVW) = H2AX(NUW,NVW) + H2X(IU,IV)
  100    CONTINUE
  200 CONTINUE
C
      IF( IPRRSP .GT. 1000 ) THEN
         WRITE(LUPRI,'(//A)') '  Two electron integrals in H2ACXY'
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)') '  Active distribution in H2AX '
         CALL OUTPUT(H2AX,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck conex */
      SUBROUTINE CONEX(KZYVA,IGRSYM,FI,ONEFAC,H2AX,VEC,NZYVEC,
     *                 NZCVEC,ISYMV,NZCSTJ,ISYMJ,LREFST,ETRS,XINDX,
     *                 ISPIN1,ISPIN2,WRK,LWRK)
C
C     Purpose:
C     Perform the configuration part of the linear transformation.
C     Expressions:
C
C           < j    | K | 0(R) >
C         - < 0(L) | K | j    >
C
C
C     The expression of the operator K between two arbitrary CSF's, in
C     terms of iF where appropriate, is:
C
C     <CSF1 | K | CSF2 > = E(inactive) +
C
C        sum(xy) iF(x,y) D(x,y) + 0.5 sum(xyvw) (xy|vw) d(xy;vw)
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, D2 = 2.0D0 )
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infopt.h"
#include "inforb.h"
#include "infind.h"
#include "infpri.h"
#include "rspprp.h"
#include "infhyp.h"
#include "infvar.h"
#include "qrinf.h"
#include "cbgetdis.h"
#include "infspi.h"
#include "infdis.h"
C
      DIMENSION FI(NORBT,NORBT)
      DIMENSION VEC(NZYVEC)
      DIMENSION WRK(*)
      DIMENSION H2AX(N2ASHX,N2ASHX)
      DIMENSION ETRS(KZYVA)
      DIMENSION XINDX(*)
C
      LOGICAL NOH2, IH8SM, LREFST,ITEST
      DATA ITEST /.FALSE./
C
C     Allocate work space for copying inactive Fock matrix and keeping
C     the CI part of the vector, initialise if necessary
C
      KFIAC  = 1
      KTRS   = KFIAC + N2ASHX
      KFREE  = KTRS  + MZCONF(IGRSYM)
      LFREE  = LWRK - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('CONEX',KFREE-1,LWRK)
C
      IADH2  = -1
C     ... IADH2 .lt. 0: H2AC in core
      NOH2   = .FALSE.
      IH8SM  = .FALSE.
      IF ( LREFST ) THEN
         ISYMCI = IREFSY
      ELSE
         ISYMCI = MULD2H(IREFSY,ISYMV)
      END IF
C
C
      IF( IPRRSP .GT.  200) THEN
         WRITE(LUPRI,'(/A/)') ' Inactive Fock matrix in CONEX'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT, 1,LUPRI)
         WRITE(LUPRI,'(/A/)') ' Vector in CONEX'
         CALL OUTPUT(VEC,1,NZYVEC,1,1,NZYVEC,1, 1,LUPRI)
      END IF
C
C    Do     < j    | K | 0(R) >
C    ==========================
C
C    Extract inactive contribution from FI in EINFI
C    Copy transpose active-active blocks of inactive Fock matrix to WRK
C
         CALL DZERO(WRK(KTRS),MZCONF(IGRSYM))
         EINFI = ONEFAC
         IF (ISPIN1.EQ.ISPIN2) THEN
            DO 105 IW = 1, NISHT
               IX = ISX(IW)
               EINFI = EINFI + FI(IX,IX)
  105       CONTINUE
         END IF
C
         DO 110 IW = 1, NASHT
            IX = ISX(NISHT + IW)
            DO 111 JW = 1, NASHT
               JX = ISX(NISHT + JW)
               IND = (IW-1) * NASHT + JW + KFIAC - 1
               WRK(IND) = FI(IX,JX)
111         CONTINUE
110      CONTINUE
C
         IF( IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(/A,F20.10)')
     *      ' Inactive contribution from FI in CONEX',EINFI
            WRITE(LUPRI,'(/A/)') ' Active parts of FI in CONEX'
            CALL OUTPUT(WRK(KFIAC),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
C
C        Tell CISIGD to use the transposed H2AX array
         IF (DIROIT) THEN
            DISTYP = INFDIS(1)
         ELSE
            DISTYP = 9
         END IF
C
         CALL CISIGD(ISYMCI,ISYMJ,NZCVEC,NZCSTJ,VEC,
     *               WRK(KTRS),WRK(KFIAC),H2AX,
     *               NOH2,IH8SM,XINDX,ISPIN1,ISPIN2,WRK(KFREE),
     *               LFREE)
C
C        Add inactive contribution
C
         IF (ISYMCI .EQ. ISYMJ .AND. ISPIN1.EQ.ISPIN2) THEN
            IF ( NISHT .GT. 0 )
     *         CALL DAXPY(NZCVEC,EINFI,VEC,1,WRK(KTRS),1)
         END IF
C
C
C        If we are testing, control for the CREF component
C
         IF( ITEST .AND. ISPIN1 .EQ. 0 .AND. ISPIN2 .EQ. 0 ) THEN
            IF (LREFST .AND. ISYMCI.EQ.ISYMJ) THEN
               T1 = DDOT(NZCVEC,VEC,1,WRK(KTRS),1)
               IF( IPRRSP .GT.  200) THEN
                  WRITE(LUPRI,*) ' CREF factor in CONEX, T1 =',T1
                  WRITE(LUPRI,'(/A/)')
     *            ' Z result in CONEX before removing CREF component'
                  CALL OUTPUT(WRK(KTRS),1,NZCVEC,1,1,NZCVEC,1,1,LUPRI)
               END IF
               CALL DAXPY(NZCVEC,(-T1),VEC,1,WRK(KTRS),1)
            END IF
         END IF
C
         IF( IPRRSP .GT.  200) THEN
            WRITE(LUPRI,'(/A/)') 'Final Z result in CONEX'
            CALL OUTPUT(WRK(KTRS),1,NZCVEC,1,1,NZCVEC,1,1,LUPRI)
         END IF
C
C        Take care of minus sign in CR if we have anything but the
C        reference state
C
         IF (LREFST) THEN
             FACT = D1
         ELSE
             FACT = DM1
         END IF
         CALL DAXPY(MZCONF(IGRSYM),FACT,WRK(KTRS),1,ETRS,1)
C
         IF (IPRRSP .GT. 200 ) THEN
            WRITE(LUPRI,'(/A,/A,/)')
     *      ' Y component of configuration part of linear',
     *      ' transformation added. Gradient is now:'
            CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
         END IF
C
C    Do   - < 0(L) | K | j    >
C    ==========================
C
         CALL DZERO(WRK(KTRS),MZCONF(IGRSYM))
C        Tell CISIGD to use the H2AX array directly
         IF (DIROIT) THEN
            DISTYP = INFDIS(2)
         ELSE
            DISTYP = 10
         END IF
C
         IF( LREFST ) THEN
            ICOFF = 1
         ELSE
            ICOFF = MZVAR(ISYMV) + 1
         END IF
         IEOFF = MZVAR(IGRSYM) + 1
C
C    Copy  active-active blocks of inactive Fock matrix to WRK
C
         DO 210 IW = 1, NASHT
            IX = ISX(NISHT + IW)
            DO 211 JW = 1, NASHT
               JX = ISX(NISHT + JW)
               IND = (IW-1) * NASHT + JW + KFIAC - 1
               WRK(IND) = FI(JX,IX)
211         CONTINUE
210      CONTINUE
C
         IF( IPRRSP .GT. 200) THEN
            WRITE(LUPRI,'(/A,F20.10)')
     *      ' Inactive contribution from FI in CONEX',EINFI
            WRITE(LUPRI,'(/A/)')
     *                     ' Transposed active parts of FI in CONEX'
            CALL OUTPUT(WRK(KFIAC),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
C
         CALL CISIGD(ISYMCI,ISYMJ,NZCVEC,NZCSTJ,VEC(ICOFF),
     *               WRK(KTRS),WRK(KFIAC),H2AX,
     *               NOH2,IH8SM,XINDX,ISPIN1,ISPIN2,WRK(KFREE),
     *               LFREE)
C
C
C        Add inactive contribution
C
         IF (ISYMCI .EQ. ISYMJ .AND. ISPIN1.EQ.ISPIN2) THEN
            IF ( NISHT .GT. 0 )
     *      CALL DAXPY(NZCVEC,EINFI,VEC(ICOFF),
     *                 1,WRK(KTRS),1)
         END IF
C
         IF( ITEST .AND. ISPIN1 .EQ. 0 .AND. ISPIN2 .EQ. 0 ) THEN
C
            IF (LREFST .AND. ISYMCI.EQ.ISYMJ) THEN
               T1 = DDOT(NZCVEC,VEC(ICOFF),1,WRK(KTRS),1)
               IF( IPRRSP .GT.  200) THEN
                  WRITE(LUPRI,*) ' CREF factor in CONEX, T1 =',T1
                  WRITE(LUPRI,'(/A/)')
     *            ' Y result in CONEX before removing CREF component'
                  CALL OUTPUT(WRK(KTRS),1,NZCVEC,1,1,NZCVEC,1,1,LUPRI)
               END IF
               CALL DAXPY(NZCVEC,(-T1),VEC(ICOFF),1,WRK(KTRS),1)
            END IF
         END IF
C
         IF( IPRRSP .GT.  200) THEN
            WRITE(LUPRI,'(/A/)') ' Final result in CONEX'
            CALL OUTPUT(WRK(KTRS),1,NZCVEC,1,1,NZCVEC,1,1,LUPRI)
         END IF
C
C        Do the minus sign of the whole term
C
         CALL DAXPY(MZCONF(IGRSYM),DM1,WRK(KTRS),1,ETRS(IEOFF),1)
C
         IF (IPRRSP .GT. 200 ) THEN
            WRITE(LUPRI,'(/A,/A,/)')
     *      ' Z component of configuration part of linear',
     *      ' transformation added. Gradient is now:'
            CALL OUTPUT(ETRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
         END IF
C
      RETURN
      END
C  /* Deck s3init */
      SUBROUTINE S3INIT(KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,
     *                  IGRSPI,ISPIN1,ISPIN2,
     *                  VECB,VECC,VECA,ATEST,
     *                  S3TRS,XINDX,UDV,MJWOP,WRK,LWRK)
C
C     Layout the core for the calculation of S3 times two vectors
C
#include "implicit.h"
#include "infdim.h"
#include "inforb.h"
#include "maxash.h"
#include "infvar.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infhyp.h"
#include "qrinf.h"
#include "infpri.h"
C
      LOGICAL   ATEST
      DIMENSION WRK(*)
      DIMENSION S3TRS(KZYVA)
      DIMENSION VECB(KZYVB)
      DIMENSION VECC(KZYVC)
      DIMENSION VECA(KZYVA)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI), MJWOP(2,MAXWOP,8)
C
C     Initialise the gradient to zero.
C
      CALL QENTER('S3INIT')
      CALL DZERO(S3TRS,KZYVA)
C
      KZYMAT  = 1
      KDEN1   = KZYMAT + NORBT * NORBT
      KFREE   = KDEN1  + NASHT * NASHT
      LWRKF   = LWRK - KFREE
      IF (LWRKF.LT.0) CALL ERRWRK('S3INIT',KFREE-1,LWRK)

      CALL S3DRV(KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,
     *           IGRSPI,ISPIN1,ISPIN2,
     *           S3TRS,VECB,VECC,VECA,ATEST,
     *           UDV,WRK(KZYMAT),WRK(KDEN1),XINDX,MJWOP,WRK(KFREE),
     *           LWRKF)
C
      CALL QEXIT('S3INIT')
      RETURN
      END
C  /* Deck s3drv */
      SUBROUTINE S3DRV(KZYVA,KZYVB,KZYVC,ISYMA,ISYMB,ISYMC,
     *                 IGRSPI,ISPIN1,ISPIN2,
     *                 S3TRS,VECB,VECC,VECA,ATEST,
     *                 UDV,ZYMAT,DEN1,XINDX,MJWOP,WRK,LWRK)
C
C     Drive the computation of S[3] times two vectors
C
#include "implicit.h"
#include "priunit.h"
#include "infdim.h"
#include "inforb.h"
#include "maxorb.h"
#include "maxash.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infhyp.h"
#include "infvar.h"
#include "infind.h"
#include "qrinf.h"
#include "infpri.h"
#include "infspi.h"
C
      PARAMETER(DM1=-1.0D0, DMH=-0.5D0, D0=0.0D0, D1=1.0D0, D2=2.0D0 )
C
      DIMENSION WRK(*)
      DIMENSION S3TRS(KZYVA)
      DIMENSION VECB(KZYVB)
      DIMENSION VECC(KZYVC)
      DIMENSION VECA(KZYVA)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION ZYMAT(NORBT,NORBT)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI), MJWOP(2,MAXWOP,8)
C
      LOGICAL ATEST, NORHO2, TDM, LORB, LCON, LREFST, LREF
C
      CALL QENTER('S3DRV')
      CALL DZERO(S3TRS, KZYVA)
C
C     Layout some workspace
C
      KCREF  = 1
      KFREE  = KCREF + MZCONF(1)
      LFREE  = LWRK  - KFREE + 1
      IF (LFREE.LT.0) CALL ERRWRK('S3DRV 1',KFREE-1,LWRK)
C
      NSIM = 1
      IGRSYM = ISYMA
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      S3TMZC = D0
      S3TMZO = D0
      S3TMYC = D0
      S3TMYO = D0
      S3TERM = D0
      KZCA   = MZCONF(ISYMA)
      KZOA   = MZWOPT(ISYMA)
      KAZC   = 1
      KAZO   = KAZC + KZCA
      KAYC   = KAZO + KZOA
      KAYO   = KAYC + KZCA
C
C     Case 1
C     ======
      IF ((MZWOPT(ISYMB) .EQ. 0).OR.(MZWOPT(ISYMC).EQ. 0)) GO TO 2000
C     Transform the kappa's
C
      CALL TRZYMT(NSIM,VECC,VECB,KZYVC,KZYVB,ISYMC,ISYMB,ZYMAT,MJWOP,
     *            WRK,LWRK)
C
C     Make copies of the MC density matrix in DEN1
C
      CALL DCOPY(NASHT*NASHT,UDV,1,DEN1,1)
      OVLAP = D1
      ISYMDN = 1
C
C     Put the factor of minus one half into the operator matrix
C
      CALL DSCAL(NORBT*NORBT,DMH,ZYMAT,1)
C
C     Make the gradient
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      INTSYM = MULD2H(ISYMB,ISYMC)
      INTSPI = MULSP(ISPIN1,ISPIN2)
      IF ( INTSYM .NE. IGRSYM ) THEN
         WRITE(LUPRI,'(/A,I5)') ' INTSYM   = ', INTSYM
         WRITE(LUPRI,'(/A,I5)') ' IGRSYM   = ', IGRSYM
         CALL QUIT(' INTSYM is not equal to IGRSYM in S3DRV')
      ENDIF
      IF ( INTSPI .NE. IGRSPI ) THEN
         WRITE(LUPRI,'(/A,I5)') ' INTSPI   = ', INTSPI
         WRITE(LUPRI,'(/A,I5)') ' IGRSPI   = ', IGRSPI
         CALL QUIT(' INTSPI is not equal to IGRSPI in S3DRV')
      END IF
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREFST = .TRUE.
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      CALL RSP1GR(NSIM,KZYVA,INTSYM,INTSPI,IGRSYM,IGRSPI,ISYMV,S3TRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,ZYMAT,
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREFST)
C
      IF (ATEST) THEN
         S3NWZC = DDOT(KZCA,VECA(KAZC),1,S3TRS(KAZC),1)
         S3NWZO = DDOT(KZOA,VECA(KAZO),1,S3TRS(KAZO),1)
         S3NWYC = DDOT(KZCA,VECA(KAYC),1,S3TRS(KAYC),1)
         S3NWYO = DDOT(KZOA,VECA(KAYO),1,S3TRS(KAYO),1)
         S3NEW  = S3NWZC + S3NWZO + S3NWYC + S3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 1    :',S3NEW-S3TERM
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 1 ZC :',S3NWZC-S3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 1 ZO :',S3NWZO-S3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 1 YC :',S3NWYC-S3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 1 YO :',S3NWYO-S3TMYO
         S3TMZC = S3NWZC
         S3TMZO = S3NWZO
         S3TMYC = S3NWYC
         S3TMYO = S3NWYO
         S3TERM = S3NEW
      END IF
C
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 1 for S3 term'
         CALL OUTPUT(S3TRS,1,KZYVA,1,1,KZYVA,1,1,LUPRI)
      END IF
C
C     Case 2
C     ======
 2000 IF((MZWOPT(ISYMB).EQ.0).OR.(MZCONF(ISYMC).EQ.0)) GO TO 3000
C
C     Construct the density matrix <02L|..|0> + <0|..|02R>
C
      OVLAP = D0
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL DZERO(DEN1,NASHT*NASHT)
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(ISYMC,IREFSY)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMC)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMC)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *            WRK(KCREF),VECC,OVLAP,DEN1,DUMMY,IGRSPI,ISPIN1,TDM,
     *            NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
C     Construct the one electron matrix
C
      CALL GTZYMT(NSIM,VECB,KZYVB,ISYMB,ZYMAT,MJWOP)
C
C     Put the factor of minus one into the operator matrix
C
      CALL DSCAL(NORBT*NORBT,DM1,ZYMAT,1)
C
C     Make the gradient
C
      INTSYM = ISYMB
      ISYMV  = ISYMC
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREFST = .FALSE.
      NZCVEC = MZCONF(ISYMC)
      NZYVEC = MZYVAR(ISYMC)
      CALL RSP1GR(NSIM,KZYVA,INTSYM,ISPIN1,IGRSYM,IGRSPI,ISYMV,S3TRS,
     *            VECC,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,ZYMAT,
     *            XINDX,MJWOP,WRK,LWRK,LORB,LCON,LREFST)
C
C
      IF (ATEST) THEN
         S3NWZC = DDOT(KZCA,VECA(KAZC),1,S3TRS(KAZC),1)
         S3NWZO = DDOT(KZOA,VECA(KAZO),1,S3TRS(KAZO),1)
         S3NWYC = DDOT(KZCA,VECA(KAYC),1,S3TRS(KAYC),1)
         S3NWYO = DDOT(KZOA,VECA(KAYO),1,S3TRS(KAYO),1)
         S3NEW  = S3NWZC + S3NWZO + S3NWYC + S3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 2    :',S3NEW-S3TERM
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 2 ZC :',S3NWZC-S3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 2 ZO :',S3NWZO-S3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 2 YC :',S3NWYC-S3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 2 YO :',S3NWYO-S3TMYO
         S3TMZC = S3NWZC
         S3TMZO = S3NWZO
         S3TMYC = S3NWYC
         S3TMYO = S3NWYO
         S3TERM = S3NEW
      END IF
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 2 for S3 term'
         CALL OUTPUT(S3TRS,1,KZYVA,1,NSIM,KZYVA,NSIM,1,LUPRI)
      END IF
C
C     Case 3
C     ======
 3000 IF((MZCONF(ISYMB).EQ.0).OR.(MZCONF(ISYMC).EQ.0)) GO TO 4000
C
C     Construct <01L|..|02R> + <02L|..|01R>
C
      OVLAP = D0
      CALL DZERO(DEN1,NASHT*NASHT)
C
      ILSYM  = MULD2H(IREFSY,ISYMB)
      IRSYM  = MULD2H(IREFSY,ISYMC)
      NCL    = MZCONF(ISYMB)
      NCR    = MZCONF(ISYMC)
      KZVARL = MZYVAR(ISYMB)
      KZVARR = MZYVAR(ISYMC)
      LREF   = .FALSE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      IF (IGRSYM .EQ. ISYMDN) THEN
         CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *               VECB,VECC,OVLAP,DEN1,DUMMY,IGRSPI,0,TDM,
     *               NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
C        Make the gradient on the fly:
C
         NZCONF = MZCONF(IGRSYM)
         NYCONF = MZCONF(IGRSYM) + MZVAR(IGRSYM)
C
         DO 500 IC = 1, MZWOPT(IGRSYM)
            K = MJWOP(1,IC,IGRSYM)
            L = MJWOP(2,IC,IGRSYM)
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF(ITYPK.EQ.JTACT .AND. ITYPL. EQ. JTACT) THEN
               NWK  = ISW(K) - NISHT
               NWL  = ISW(L) - NISHT
                  S3TRS(NZCONF+IC) = S3TRS(NZCONF+IC)
     *                                         + DEN1(NWK,NWL)
                  S3TRS(NYCONF+IC) = S3TRS(NYCONF+IC)
     *                                         + DEN1(NWL,NWK)
            END IF
500      CONTINUE
      END IF
C
      IF (ATEST) THEN
         S3NWZC = DDOT(KZCA,VECA(KAZC),1,S3TRS(KAZC),1)
         S3NWZO = DDOT(KZOA,VECA(KAZO),1,S3TRS(KAZO),1)
         S3NWYC = DDOT(KZCA,VECA(KAYC),1,S3TRS(KAYC),1)
         S3NWYO = DDOT(KZOA,VECA(KAYO),1,S3TRS(KAYO),1)
         S3NEW  = S3NWZC + S3NWZO + S3NWYC + S3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 3    :',S3NEW-S3TERM
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 3 ZC :',S3NWZC-S3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 3 ZO :',S3NWZO-S3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 3 YC :',S3NWYC-S3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 3 YO :',S3NWYO-S3TMYO
         S3TMZC = S3NWZC
         S3TMZO = S3NWZO
         S3TMYC = S3NWYC
         S3TMYO = S3NWYO
         S3TERM = S3NEW
      END IF
C
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 3 for S3 term'
         CALL OUTPUT(S3TRS,1,KZYVA,1,NSIM,KZYVA,NSIM,1,LUPRI)
      END IF
C
C     Case 4
C     ======
 4000 IF((MZWOPT(ISYMB).EQ.0).OR.(MZCONF(ISYMC).EQ.0)) GO TO 5000
C     Do Sj(2) <0 |-kappa(1)| 0>
C
      IF (ISYMB .EQ. 1) THEN
C
C        Compute the expectation value
C
         KSYMOP = ISYMB
         CALL GTZYMT(NSIM,VECB,KZYVB,ISYMB,ZYMAT,MJWOP)
C
         IPRONE = 75
         TRPLET = ISPIN1.EQ.1
         CALL PRPONE(ZYMAT,UDV,FACT,IPRONE,'<0 |-kappa(1)| 0>')
C
         IF (FACT .NE. D0) THEN
            NZCONF = MZCONF(IGRSYM)
            NZVAR  = MZVAR(IGRSYM)
               CALL DAXPY(NZCONF,-FACT,VECC,1,S3TRS,1)
               CALL DAXPY(NZCONF,-FACT,VECC(NZVAR+1),1,
     *                    S3TRS(NZVAR+1),1)
         END IF
      END IF
C
      IF (ATEST) THEN
         S3NWZC = DDOT(KZCA,VECA(KAZC),1,S3TRS(KAZC),1)
         S3NWZO = DDOT(KZOA,VECA(KAZO),1,S3TRS(KAZO),1)
         S3NWYC = DDOT(KZCA,VECA(KAYC),1,S3TRS(KAYC),1)
         S3NWYO = DDOT(KZOA,VECA(KAYO),1,S3TRS(KAYO),1)
         S3NEW  = S3NWZC + S3NWZO + S3NWYC + S3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 4    :',S3NEW-S3TERM
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 4 ZC :',S3NWZC-S3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 4 ZO :',S3NWZO-S3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 4 YC :',S3NWYC-S3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 4 YO :',S3NWYO-S3TMYO
         S3TMZC = S3NWZC
         S3TMZO = S3NWZO
         S3TMYC = S3NWYC
         S3TMYO = S3NWYO
         S3TERM = S3NEW
      END IF
C
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 4 for S3 term'
         CALL OUTPUT(S3TRS,1,KZYVA,1,NSIM,KZYVA,NSIM,1,LUPRI)
      END IF
C
C     Case 5:
C     =======
 5000 IF((MZCONF(ISYMB).EQ.0).OR.(MZCONF(ISYMC).EQ.0)) GO TO 6000
C
      IF ( (ISYMB .EQ. ISYMC) .AND. (NASHT .GT. 0) ) THEN
            FACT = D0
            ICOFF = MZVAR(ISYMC)
            DO 730 IC = 1, MZCONF(ISYMC)
               FACT = FACT + VECB(IC) * VECC(ICOFF + IC)
               FACT = FACT + VECB(ICOFF + IC) * VECC(IC)
 730        CONTINUE
            IF (IPRRSP .GT. 100) THEN
               WRITE(LUPRI,'(/A,F20.10)') ' Factor in case 5 = ',FACT
            END IF
C
            NZCONF = MZCONF(IGRSYM)
            NYCONF = MZCONF(IGRSYM) + MZVAR(IGRSYM)
C
            IF (IGRSPI.EQ.1) THEN
               KCREF   = 1
               KFREE   = KCREF + MZCONF(1)
               LFREE   = LWRK  - KFREE + 1
               IF (LFREE.LT.0) CALL ERRWRK('S3DRV 2',KFREE-1,LWRK)
               CALL GETREF(WRK(KCREF),MZCONF(1))
               CALL RSPDM(IREFSY,IREFSY,MZCONF(1),MZCONF(1),
     *                 WRK(KCREF),WRK(KCREF), DEN1,RHO2,
     *                 IGRSPI,0,.FALSE.,.TRUE.,XINDX,WRK, KFREE,LFREE)
            ELSE
               CALL DCOPY(N2ASHX,UDV,1,DEN1,1)
            END IF
C
            DO 710 IC = 1, MZWOPT(IGRSYM)
               K = MJWOP(1,IC,IGRSYM)
               L = MJWOP(2,IC,IGRSYM)
               ITYPK = IOBTYP(K)
               ITYPL = IOBTYP(L)
               IF(ITYPK.EQ.JTACT .AND. ITYPL. EQ. JTACT) THEN
                  NWK  = ISW(K) - NISHT
                  NWL  = ISW(L) - NISHT
                  S3TRS(NZCONF+IC) = S3TRS(NZCONF+IC)
     *                               + FACT * DEN1(NWK,NWL)
                  S3TRS(NYCONF+IC) = S3TRS(NYCONF+IC)
     *                               + FACT * DEN1(NWL,NWK)
               END IF
710         CONTINUE
      END IF
C
      IF (ATEST) THEN
         S3NWZC = DDOT(KZCA,VECA(KAZC),1,S3TRS(KAZC),1)
         S3NWZO = DDOT(KZOA,VECA(KAZO),1,S3TRS(KAZO),1)
         S3NWYC = DDOT(KZCA,VECA(KAYC),1,S3TRS(KAYC),1)
         S3NWYO = DDOT(KZOA,VECA(KAYO),1,S3TRS(KAYO),1)
         S3NEW  = S3NWZC + S3NWZO + S3NWYC + S3NWYO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 5    :',S3NEW-S3TERM
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 5 ZC :',S3NWZC-S3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 5 ZO :',S3NWZO-S3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 5 YC :',S3NWYC-S3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Case 5 YO :',S3NWYO-S3TMYO
         S3TMZC = S3NWZC
         S3TMZO = S3NWZO
         S3TMYC = S3NWYC
         S3TMYO = S3NWYO
         S3TERM = S3NEW
      END IF
C
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 5 for S3 term'
         CALL OUTPUT(S3TRS,1,KZYVA,1,NSIM,KZYVA,NSIM,1,LUPRI)
      END IF
C
C
6000  CONTINUE
C
C     End of transformation; print final result if required
C
      IF (ATEST) THEN
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Final result    :',S3TERM
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Final result ZC :',S3TMZC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Final result ZO :',S3TMZO
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Final result YC :',S3TMYC
         WRITE(LUPRI,'(A,F24.8)') ' S3TEST Final result YO :',S3TMYO
      END IF
      CALL QEXIT('S3DRV')
      RETURN
      END
C  /* Deck gtzymt */
      SUBROUTINE GTZYMT(NSIM,VEC,KZYV,ISYMV,ZYMAT,MJWOP)
C
C     This subroutine unpacks the ZY matrix from the vector. It uses
C     the Z and the Y part of the vector.
C
#include "implicit.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "qrinf.h"
#include "infpri.h"
C
      DIMENSION VEC(KZYV)
      DIMENSION ZYMAT(NORBT,NORBT), MJWOP(2,MAXWOP,8)
C
      IF (NSIM .NE. 1) THEN
         CALL QUIT('GTZYMT ERROR: NSIM .ne. 1 not implemented')
      END IF
      IF( IPRRSP .GT. 200 ) THEN
         WRITE(LUPRI,'(/A)') ' GTZYMT called with the vector'
         WRITE(LUPRI,'(A)')  ' ============================='
         CALL OUTPUT(VEC,1,KZYV,1,NSIM,KZYV,NSIM,1,LUPRI)
      END IF
C
      IOB  = MZCONF(ISYMV)
      IOBT = MZVAR(ISYMV) + MZCONF(ISYMV)
      CALL DZERO(ZYMAT,NORBT*NORBT*NSIM)
      DO 100 ISIM = 1, NSIM
         DO 200 IG = 1, MZWOPT(ISYMV)
            I = MJWOP(1,IG,ISYMV)
            J = MJWOP(2,IG,ISYMV)
            ZYMAT(J,I) = VEC(IOB  + IG)
            ZYMAT(I,J) = VEC(IOBT + IG)
200      CONTINUE
C
         IF( IPRRSP .GT. 200 ) THEN
            WRITE(LUPRI,'(/A)') ' Vector unpacked in matrix'
            WRITE(LUPRI,'(A)')  ' ========================='
            CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
         END IF
C
100   CONTINUE
C
      RETURN
      END
C  /* Deck trzymt */
      SUBROUTINE TRZYMT(NSIM,VEC1,VEC2,KZYV1,KZYV2,ISYM1,ISYM2,ZYMAT,
     *                  MJWOP,WRK,LWRK)
C
C     This subroutine unpacks the ZY matrices from the two vectors and
C     does the transformation
C
#include "implicit.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "qrinf.h"
#include "infpri.h"
C
      DIMENSION VEC1(KZYV1), VEC2(KZYV2)
      DIMENSION ZYMAT(NORBT,NORBT), MJWOP(2,MAXWOP,8)
      DIMENSION WRK(*)
C
      IF (NSIM .NE. 1) THEN
         CALL QUIT('TRZYMT ERROR: NSIM .ne. 1 not implemented')
      END IF
      KZY1  = 1
      KZY2  = KZY1  + NSIM * NORBT * NORBT
      KFREE = KZY2  + NSIM * NORBT * NORBT
      LFREE = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('TRZYMT',KFREE-1,LWRK)
C
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYM1,WRK(KZY1),MJWOP )
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYM2,WRK(KZY2),MJWOP )
C
      DO 100 ISIM = 1, NSIM
         IZY1 = (ISIM - 1) * NORBT * NORBT + KZY1
         IZY2 = (ISIM - 1) * NORBT * NORBT + KZY2
C
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,1.D0,
     &              WRK(IZY1),NORBT,
     &              WRK(IZY2),NORBT,0.D0,
     &              ZYMAT,NORBT)
C
         CALL DGEMM('N','N',NORBT,NORBT,NORBT,-1.D0,
     &              WRK(IZY2),NORBT,
     &              WRK(IZY1),NORBT,1.D0,
     &              ZYMAT,NORBT)
C
         IF( IPRRSP .GT. 200 ) THEN
            WRITE(LUPRI,'(/A)') ' Final result in TRZYMT'
            WRITE(LUPRI,'(A)')  ' ======================'
            CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
         END IF
C
100   CONTINUE
C
      RETURN
      END
C  /* Deck x2init */
      SUBROUTINE X2INIT(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *                 GRVEC,VEC,X2TRS,XINDX,UDV,PV,OPLBL,IOPSYM,IOPSPI,
     *                 CMO,MJWOP,WRK,LWRK)
C
C     Layout the core for the calculation of X2 times one vector
C
#include "implicit.h"
#include "infdim.h"
#include "inforb.h"
#include "maxash.h"
#include "priunit.h"
#include "infvar.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infhyp.h"
#include "qrinf.h"
#include "infpri.h"
#include "infspi.h"
C
      CHARACTER*8 OPLBL
C
      DIMENSION WRK(*)
      DIMENSION X2TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION GRVEC(KZYVR), VEC(KZYV)
      DIMENSION XINDX(*)
      DIMENSION CMO(*)
      DIMENSION UDV(NASHDI,NASHDI), PV(*)
C
C     Any variables in this symmetry?
C
      IF (KZYVR .EQ. 0) RETURN
C
      IF( (NSIM .GT. 1) ) THEN
         WRITE(LUERR,'(//A,/A,I5,/A)')
     *   ' --> Error : Routine X2INIT not implemented  with NSIM > 1',
     *   '     NSIM = ', NSIM,
     *   '     Halting execution.'
          CALL QTRACE(LUERR)
          CALL QUIT(' RSPQR ERROR: NSIM GREATER THAN 1 IN X2INIT')
      END IF
C
C     Initialise the gradient to zero.
C
      CALL DZERO(X2TRS,KZYVR)
C
      KOPMAT  = 1
      KDEN1   = KOPMAT + NORBT * NORBT
      KDEN2   = KDEN1  + NASHT * NASHT
      IF (OPLBL(3:8).EQ.'SPNORB' .OR. OPLBL(3:8).EQ.'MNF-SO'
     &   .OR. OPLBL(3:8).EQ.'SPNSCA') THEN
         KFREE = KDEN2 + N2ASHX*N2ASHX*2
      ELSE
         KFREE = KDEN2
      END IF
      LWRKF   = LWRK - KFREE + 1
      IF (LWRKF.LT.0) CALL ERRWRK('X2INIT',KFREE-1,LWRK)
C
      IF (IPRRSP .GT. 50) THEN
         WRITE(LUPRI,'(//A)')  ' Characteristics of X2 gradient'
         WRITE(LUPRI,'(A)')    ' =============================='
         WRITE(LUPRI,'(A,I8)') ' Gradient symmetry :',IGRSYM
         WRITE(LUPRI,'(A,I8)') ' Gradient spin     :',IGRSPI
         WRITE(LUPRI,'(A,I8)') ' Gradient length   :',KZYVR
         WRITE(LUPRI,'(A,I8)') ' Vector symmetry   :',ISYMV
         WRITE(LUPRI,'(A,I8)') ' Vector spin       :',ISPINV
         WRITE(LUPRI,'(A,I8)') ' Vector length     :',KZYV
         WRITE(LUPRI,'(A,I8)') ' Operator symmetry :',IOPSYM
         WRITE(LUPRI,'(A,I8)') ' Operator spin     :',IOPSPI
         WRITE(LUPRI,'(A,A8)') ' Operator label    :',OPLBL
         WRITE(LUPRI,'(A//)')  ' =============================='
      END IF
C
      IF (OPLBL(3:8).EQ.'SPNORB' .OR. OPLBL(3:8).EQ.'MNF-SO'
     &   .OR. OPLBL(3:8).EQ.'SPNSCA') THEN
         CALL X2HSO(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *           OPLBL,IOPSYM,IOPSPI,GRVEC,X2TRS,
     *           VEC,UDV,PV,WRK(KOPMAT),WRK(KDEN1),WRK(KDEN2),
     *           XINDX,CMO,MJWOP,WRK(KFREE),LWRKF)
      ELSE
         CALL X2DRV(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *           OPLBL,IOPSYM,IOPSPI,GRVEC,X2TRS,
     *           VEC,UDV,WRK(KOPMAT),WRK(KDEN1),
     *           XINDX,CMO,MJWOP,WRK(KFREE),LWRKF)
      END IF
C
      RETURN
      END
C  /* Deck x2drv */
      SUBROUTINE X2DRV(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *                 OPLBL,IOPSYM,IOPSPI,
     *                 GRVEC,X2TRS,VEC,UDV,OPMAT,DEN1,XINDX,CMO,
     *                 MJWOP,WRK,LWRK)
C
C     Drive the computation of X[2] times one vector
C
#include "implicit.h"
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0 )
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "infvar.h"
#include "infdim.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infhyp.h"
#include "qrinf.h"
#include "infpri.h"
#include "infspi.h"
C
      CHARACTER*8 OPLBL
C
      DIMENSION WRK(*)
      DIMENSION X2TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION GRVEC(KZYVR), VEC(KZYV)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION XINDX(*)
      DIMENSION CMO(*)
      DIMENSION UDV(NASHDI,NASHDI)
C
      LOGICAL NORHO2, TDM, LORB, LCON, LREFST, LREF
C
C     Layout some workspace
C
      KCREF  = 1
      KZYM   = KCREF + MZCONF(1)
      KX2X   = KZYM  + N2ORBX
      KFREE  = KX2X  + N2ORBX
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('X2DRV 1',KFREE-1,LWRK)
C
      IF (X2TEST) THEN
         X2TMZC = D0
         X2TMZO = D0
         X2TMYC = D0
         X2TMYO = D0
         X2TM = D0
         KZX2O = MZWOPT(IGRSYM)
         KZX2C = MZCONF(IGRSYM)
         KZC = 1
         KZO = KZC + KZX2C
         KYC = KZO + KZX2O
         KYO = KYC + KZX2C
      END IF
C
      KSYMOP = IOPSYM
      TDM    = .TRUE.
      NORHO2 = .TRUE.
C
C     Get the operator matrix
C     Put the minus sign in the whole term into the operator
C
      KSYMP = -1
      CALL PRPGET(OPLBL,CMO,OPMAT,KSYMP,ANTSYM,WRK,LWRK,IPRRSP)
      IF (KSYMP .NE. KSYMOP) CALL QUIT(
     &   'ERROR: Unexpected symmetry of property matrix')
      CALL DSCAL(NORBT*NORBT,DM1,OPMAT,1)
C
C     Case 1
C     ======
      IF (MZWOPT(ISYMV) .EQ. 0 ) GOTO 2000
C
      ISPIN = MULSP(ISPINV,IOPSPI)
      IF (NASHT.EQ.0 .AND. ISPIN.NE.IGRSPI) THEN
         WRITE(LUPRI,*) "FOUND: ", ISPIN, " Expected:", IGRSPI
         WRITE(LUPRI,*) "parameters: ", ispinv, iopspi
         CALL QUIT('X2DRV: SPIN ERROR')
      ENDIF
C
C     Transform the operator
C
      CALL GTZYMT(NSIM,VEC,KZYV,ISYMV,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KX2X),N2ORBX)
      CALL OITH1(ISYMV,WRK(KZYM),OPMAT,WRK(KX2X),IOPSYM)
C
C     Make copies of the MC density matrix in DEN1
C
      IF (ISPIN.EQ.IGRSPI) THEN
         CALL DCOPY(N2ASHX,UDV,1,DEN1,1)
      ELSE
         CALL DCOPY(N2ASHX,UDV(1+N2ASHX,1),1,DEN1,1)
      END IF
      OVLAP = D1
      ISYMDN = 1
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Make the gradient
C
      ISYMR  = IREFSY
      INTSYM = MULD2H(KSYMOP,ISYMV)
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      LREFST = .TRUE.
      CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,X2TRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,WRK(KX2X),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREFST)
C
      IF (X2TEST) THEN
         X2ZO = DDOT(KZX2O,GRVEC(KZO),1,X2TRS(KZO),1)
         X2ZC = DDOT(KZX2C,GRVEC(KZC),1,X2TRS(KZC),1)
         X2YO = DDOT(KZX2O,GRVEC(KYO),1,X2TRS(KYO),1)
         X2YC = DDOT(KZX2C,GRVEC(KYC),1,X2TRS(KYC),1)
         X2TOT = X2ZO+X2ZC+X2YO+X2YC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 ZO :',X2ZO - X2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 ZC :',X2ZC - X2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 YO :',X2YO - X2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 YC :',X2YC - X2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1    :',X2TOT - X2TM
         X2TMZO = X2ZO
         X2TMZC = X2ZC
         X2TMYO = X2YO
         X2TMYC = X2YC
         X2TM = X2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 1 for X2 term'
         CALL OUTPUT(X2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Case 2
C     ======
 2000 IF (MZCONF(ISYMV) .EQ. 0 ) GOTO 3000
C     Construct the density matrix <02L|..|0> + <0|..|02R>
C
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL DZERO(DEN1,NSIM*NASHT*NASHT)
      OVLAP = D0
C
      ISPIN1 = IGRSPI
      ISPIN2 = IOPSPI
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *            WRK(KCREF),VEC,OVLAP,DEN1,DUMMY,ISPIN1,ISPIN2,
     *            TDM,NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
C     Make the gradient
C
      ISYMR  = ISYMV
      INTSYM = KSYMOP
      ISPIN = IOPSPI
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      NZYVEC = MZYVAR(ISYMV)
      NZCVEC = MZCONF(ISYMV)
      LREFST = .FALSE.
      CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,X2TRS,
     *            VEC,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,OPMAT,
     *            XINDX,MJWOP,WRK,LWRK,LORB,LCON,LREFST)
C
      IF (X2TEST) THEN
         X2ZO = DDOT(KZX2O,GRVEC(KZO),1,X2TRS(KZO),1)
         X2ZC = DDOT(KZX2C,GRVEC(KZC),1,X2TRS(KZC),1)
         X2YO = DDOT(KZX2O,GRVEC(KYO),1,X2TRS(KYO),1)
         X2YC = DDOT(KZX2C,GRVEC(KYC),1,X2TRS(KYC),1)
         X2TOT = X2ZO+X2ZC+X2YO+X2YC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 ZO :',X2ZO - X2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 ZC :',X2ZC - X2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 YO :',X2YO - X2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 YC :',X2YC - X2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2    :',X2TOT - X2TM
         X2TMZO = X2ZO
         X2TMZC = X2ZC
         X2TMYO = X2YO
         X2TMYC = X2YC
         X2TM = X2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 2 for X2 term'
         CALL OUTPUT(X2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Case 3
C     ======
C     Do Sj(2) <0 |op| 0>
C
C     Compute the expectation value
C
      IF (KSYMOP .EQ. 1) THEN
C
         ISPIN1 = IOPSPI
         ISPIN2 = 0
         TRPLET = IOPSPI.EQ.1
         IF (RSPCI.OR.TRPLET) THEN
            KCREF   = 1
            KDEN1   = KCREF + MZCONF(1)
            KFREE   = KDEN1 + N2ASHX
            LFREE   = LWRK  - KFREE
            IF (KFREE.LT.0) CALL ERRWRK('X2DRV 2',KFREE-1,LWRK)
            CALL GETREF(WRK(KCREF),MZCONF(1))
            CALL RSPDM(IREFSY,IREFSY,MZCONF(1),MZCONF(1),
     *                 WRK(KCREF),WRK(KCREF),WRK(KDEN1),RHO2,
     *                 ISPIN1,ISPIN2,.FALSE.,.TRUE.,XINDX,WRK,
     *                 KFREE,LFREE)
            IPRONE = 75
            CALL PRPONE(OPMAT,WRK(KDEN1),FACT,IPRONE,' Sj(2) <0|op|0>')
         ELSE
            IPRONE = 75
            CALL PRPONE(OPMAT,UDV,FACT,IPRONE,' Sj(2) <0|op|0>')
         END IF
C
         IF (FACT .NE. D0) THEN
            NZCONF = MZCONF(IGRSYM)
            NZVAR  = MZVAR(IGRSYM)
            CALL DAXPY(NZCONF,FACT,VEC         ,1,X2TRS         ,1)
            CALL DAXPY(NZCONF,FACT,VEC(NZVAR+1),1,X2TRS(NZVAR+1),1)
         END IF
      END IF
C
      IF (X2TEST) THEN
         X2ZO = DDOT(KZX2O,GRVEC(KZO),1,X2TRS(KZO),1)
         X2ZC = DDOT(KZX2C,GRVEC(KZC),1,X2TRS(KZC),1)
         X2YO = DDOT(KZX2O,GRVEC(KYO),1,X2TRS(KYO),1)
         X2YC = DDOT(KZX2C,GRVEC(KYC),1,X2TRS(KYC),1)
         X2TOT = X2ZO+X2ZC+X2YO+X2YC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 3 ZO :',X2ZO - X2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 3 ZC :',X2ZC - X2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 3 YO :',X2YO - X2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 3 YC :',X2YC - X2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 3    :',X2TOT - X2TM
         X2TMZO = X2ZO
         X2TMZC = X2ZC
         X2TMYO = X2YO
         X2TMYC = X2YC
         X2TM = X2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 3 for X2 term'
         CALL OUTPUT(X2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
 3000 CONTINUE
      IF (X2TEST) THEN
         WRITE(LUPRI,'(/A,F24.8)')' X2TEST Final result:  ZO', X2ZO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:  ZC', X2ZC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:  YO', X2YO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:  YC', X2YC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:    ', X2TOT
      END IF
      RETURN
      END
C  /* Deck a2init */
      SUBROUTINE A2INIT(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *                  GRVEC,VEC,A2TRS,XINDX,UDV,PV,
     *                  OPLBL,IOPSYM,IOPSPI,
     *                  CMO,MJWOP,WRK,LWRK)
C
C     Layout the core for the calculation of A[2] times one vector
C
#include "implicit.h"
#include "maxash.h"
#include "priunit.h"
#include "infdim.h"
#include "inforb.h"
#include "infvar.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "qrinf.h"
#include "infpri.h"
C
      CHARACTER*8 OPLBL
C
      DIMENSION WRK(*)
      DIMENSION A2TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION GRVEC(KZYVR), VEC(KZYV)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION CMO(*)
C
C     Any variables in this symmetry?
C
      IF (KZYVR .EQ. 0) RETURN
C
      IF( (NSIM .GT. 1) ) THEN
         WRITE(LUERR,'(//A,/A,I5,/A)')
     *   ' --> Error : Routine A2INIT not implemented  with NSIM > 1',
     *   '     NSIM = ', NSIM,
     *   '     Halting execution.'
          CALL QTRACE(LUERR)
          CALL QUIT(' RSPQR ERROR: NSIM GREATER THAN 1 IN A2INIT')
      END IF
C
C     Initialise the gradient to zero.
C
      CALL DZERO(A2TRS,KZYVR)
C
      KOPMAT  = 1
      KFREE   = KOPMAT + NORBT * NORBT
      LWRKF   = LWRK - KFREE + 1
      IF (LWRKF.LT.0) CALL ERRWRK('A2INIT',KFREE-1,LWRK)
      CALL DZERO(WRK,KFREE-1)
C
      IF (IPRRSP .GT. 50 ) THEN
         WRITE(LUPRI,'(//A)')  ' Characteristics of A2 gradient'
         WRITE(LUPRI,'(A)')    ' =============================='
         WRITE(LUPRI,'(A,I8)') ' Gradient symmetry :',IGRSYM
         WRITE(LUPRI,'(A,I8)') ' Gradient spin     :',IGRSPI
         WRITE(LUPRI,'(A,I8)') ' Gradient length   :',KZYVR
         WRITE(LUPRI,'(A,I8)') ' Vector symmetry   :',ISYMV
         WRITE(LUPRI,'(A,I8)') ' Vector spin       :',ISPINV
         WRITE(LUPRI,'(A,I8)') ' Vector length     :',KZYV
         WRITE(LUPRI,'(A,I8)') ' Operator symmetry :',IOPSYM
         WRITE(LUPRI,'(A,I8)') ' Operator spin     :',IOPSPI
         WRITE(LUPRI,'(A,A8)') ' Operator label    :',OPLBL
         WRITE(LUPRI,'(A//)')  ' =============================='
      END IF
C
      IF (OPLBL(3:8).EQ.'SPNORB' .OR. OPLBL(3:8).EQ.'MNF-SO'
     &   .OR. OPLBL(3:8).EQ.'SPNSCA') THEN
         CALL A2HSO(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *           OPLBL,IOPSYM,IOPSPI,
     *           GRVEC,A2TRS,VEC,UDV,PV,WRK(KOPMAT),XINDX,CMO,
     *           MJWOP,WRK(KFREE),LWRKF)
      ELSE
         CALL A2DRV(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *           OPLBL,IOPSYM,IOPSPI,
     *           GRVEC,A2TRS,VEC,UDV,WRK(KOPMAT),XINDX,CMO,
     *           MJWOP,WRK(KFREE),LWRKF)
C
      END IF
      RETURN
      END
C  /* Deck a2drv */
      SUBROUTINE A2DRV(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *                 OPLBL,IOPSYM,IOPSPI,
     *                 GRVEC,A2TRS,VEC,UDV,OPMAT,XINDX,CMO,
     *                 MJWOP,WRK,LWRK)
C
C     Drive the computation of A[2] times one vector
C
#include "implicit.h"
C
      PARAMETER( D0 = 0.0D0, D1 = 1.0D0, DH = 0.5D0 )
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infdim.h"
#include "inforb.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "infvar.h"
#include "infind.h"
#include "qrinf.h"
#include "infpri.h"
#include "infspi.h"
C
      CHARACTER*8 OPLBL
C
      DIMENSION WRK(*)
      DIMENSION A2TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION GRVEC(KZYVR), VEC(KZYV)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION CMO(*)
C
      LOGICAL NORHO2, TDM, LORB, LCON, LREFST, LREF
C
C     Layout some workspace
C
      KCREF  = 1
      KZYM   = KCREF + MZCONF(1)
      KA2X   = KZYM  + N2ORBX
      KFREE  = KA2X  + N2ORBX
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('A2DRV 1',KFREE-1,LWRK)
C
      TDM    = .TRUE.
      NORHO2 = .TRUE.
C
      IF (A2TEST) THEN
         A2TMZC = D0
         A2TMZO = D0
         A2TMYC = D0
         A2TMYO = D0
         A2TM = D0
         KZA2O = MZWOPT(IGRSYM)
         KZA2C = MZCONF(IGRSYM)
         KZC = 1
         KZO = KZC + KZA2C
         KYC = KZO + KZA2O
         KYO = KYC + KZA2C
      END IF
C     Get the operator matrix
C     910408-hjaaj: transpose the operator matrix to
C                   get right sign for imaginary operators
C
      KSYMOP = IOPSYM
      KSYMP = -1
      CALL PRPGET(OPLBL,CMO,OPMAT,KSYMP,ANTSYM,WRK,LWRK,IPRRSP)
      IF (KSYMP .NE. KSYMOP) CALL QUIT(
     &   'ERROR: Unexpected symmetry of property matrix')
      CALL DGETRN(OPMAT,NORBT,NORBT)
C
C     Case 1:
C     =======
C     / 1/2 <0| [qj+,A(K)] |0> \
C     |   - <0| A(K) |j>       |
C     | 1/2 <0| [qj ,A(K)] |0> |
C     \     <j| A(K) |0>       /
C
      IF ( MZWOPT(ISYMV) .EQ. 0 ) GOTO 2000
C
      ISPIN = MULSP(ISPINV,IOPSPI)
      IF (NASHT.EQ.0 .AND. ISPIN.NE.IGRSPI) THEN
         WRITE(LUPRI,*) "FOUND: ", ISPIN, " Expected:", IGRSPI
         WRITE(LUPRI,*) "parameters: ", ispinv, iopspi
        CALL QUIT('A2DRV: SPIN ERROR')
      END IF
C
C     Transform the operator
C
      INTSYM = MULD2H(ISYMV,IOPSYM)
      CALL GTZYMT(NSIM,VEC,KZYV,ISYMV,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KA2X),N2ORBX)
      CALL OITH1(ISYMV,WRK(KZYM),OPMAT,WRK(KA2X),IOPSYM)
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Make the gradient
C
      OVLAP  = D1
      ISYMDN = 1
      ISYMR  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF(ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB = (MZWOPT(IGRSYM) .GT. 0 )
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      LREFST = .TRUE.
      IF (LCON) THEN
         CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,A2TRS,
     *               WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,
     *               WRK(KA2X),XINDX,MJWOP,WRK(KFREE),LFREE,
     *               .FALSE.,LCON,LREFST)
      END IF
      IF (LORB) THEN
C        multiply A(K) with 1/2 for orbital part
         CALL DSCAL(N2ORBX,DH,WRK(KA2X),1)
         CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,A2TRS,
     *               WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,
     *               WRK(KA2X),XINDX,MJWOP,WRK(KFREE),LFREE,
     *               LORB,.FALSE.,LREFST)
      END IF
C
      IF (A2TEST) THEN
         A2ZO = DDOT(KZA2O,GRVEC(KZO),1,A2TRS(KYO),1)
         A2ZC = DDOT(KZA2C,GRVEC(KZC),1,A2TRS(KYC),1)
         A2YO = DDOT(KZA2O,GRVEC(KYO),1,A2TRS(KZO),1)
         A2YC = DDOT(KZA2C,GRVEC(KYC),1,A2TRS(KZC),1)
         A2TOT = A2ZO+A2ZC+A2YO+A2YC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 ZO :',A2ZO - A2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 ZC :',A2ZC - A2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 YO :',A2YO - A2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 YC :',A2YC - A2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1    :',A2TOT - A2TM
         A2TMZO = A2ZO
         A2TMZC = A2ZC
         A2TMYO = A2YO
         A2TMYC = A2YC
         A2TM   = A2TOT
      END IF
      IF ( IPRRSP .GT. 30 ) THEN
         WRITE(LUPRI,'(A)') ' Case 1 for A2 term'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Case 2:
C     =======
2000  IF (MZCONF(ISYMV) .EQ. 0 ) GO TO 4000
C
C     Multiply the factor of one half into the operator matrix
C     for evaluation of Case 2 and 3.
C
      CALL DSCAL(NSIM*NORBT*NORBT,DH,OPMAT,1)
C
      ISPIN = IOPSPI
C
C     Make the gradient
C
      OVLAP  = D0
      ISYMDN = 1
      INTSYM = IOPSYM
      ISYMR  = ISYMV
      NZYVEC = MZYVAR(ISYMV)
      NZCVEC = MZCONF(ISYMV)
      LORB   = .FALSE.
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF(ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LREFST = .FALSE.
      CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,A2TRS,
     *            VEC,NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,OPMAT,
     *            XINDX,MJWOP,WRK,LWRK,LORB,LCON,LREFST)
C
      IF (A2TEST) THEN
         A2ZO = DDOT(KZA2O,GRVEC(KZO),1,A2TRS(KYO),1)
         A2ZC = DDOT(KZA2C,GRVEC(KZC),1,A2TRS(KYC),1)
         A2YO = DDOT(KZA2O,GRVEC(KYO),1,A2TRS(KZO),1)
         A2YC = DDOT(KZA2C,GRVEC(KYC),1,A2TRS(KZC),1)
         A2TOT = A2ZO+A2ZC+A2YO+A2YC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 ZO :',A2ZO - A2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 ZC :',A2ZC - A2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 YO :',A2YO - A2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 YC :',A2YC - A2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2    :',A2TOT - A2TM
         A2TMZO = A2ZO
         A2TMZC = A2ZC
         A2TMYO = A2YO
         A2TMYC = A2YC
         A2TM   = A2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 2 for A2 term'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Case 3:
C     =======
C     Do Sj(2) <0 |op| 0>
C
C     Compute the expectation value (factor of 0.5 multiplied into OPMAT)
C
      IF ( IOPSYM .EQ. 1 ) THEN
C
         TRPLET = IOPSPI.EQ.1
         ISPIN1 = IOPSPI
         ISPIN2 = 0
         IF (RSPCI.OR.TRPLET) THEN
            KCREF   = 1
            KDEN1   = KCREF + MZCONF(1)
            KFREE   = KDEN1 + N2ASHX
            LFREE   = LWRK  - KFREE
            IF (LFREE.LT.0) CALL ERRWRK('A2DRV 2',KFREE-1,LWRK)
            CALL GETREF(WRK(KCREF),MZCONF(1))
            CALL RSPDM(IREFSY,IREFSY,MZCONF(1),MZCONF(1),
     *                 WRK(KCREF),WRK(KCREF),WRK(KDEN1),RHO2,
     *                 ISPIN1,ISPIN2,.FALSE.,.TRUE.,XINDX,WRK,
     *                 KFREE,LFREE)
            IPRONE = 75
            CALL PRPONE(OPMAT,WRK(KDEN1),FACT,IPRONE,
     *         ' A(2) TERM <0|op|0>')
         ELSE
            IPRONE = 75
            CALL PRPONE(OPMAT,UDV,FACT,IPRONE,
     *         ' A(2) TERM <0|op|0>')
         END IF
C
         IF ( FACT .NE. D0 ) THEN
            NZCONF = MZCONF(IGRSYM)
            NZVAR  = MZVAR(IGRSYM)
            CALL DAXPY(NZCONF,FACT,VEC,1,A2TRS,1)
            CALL DAXPY(NZCONF,FACT,VEC(NZVAR+1),1,A2TRS(NZVAR+1),1)
         END IF
      END IF
C
      IF (A2TEST) THEN
         A2ZO = DDOT(KZA2O,GRVEC(KZO),1,A2TRS(KYO),1)
         A2ZC = DDOT(KZA2C,GRVEC(KZC),1,A2TRS(KYC),1)
         A2YO = DDOT(KZA2O,GRVEC(KYO),1,A2TRS(KZO),1)
         A2YC = DDOT(KZA2C,GRVEC(KYC),1,A2TRS(KZC),1)
         A2TOT = A2ZO+A2ZC+A2YO+A2YC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 3 ZO :',A2ZO - A2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 3 ZC :',A2ZC - A2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 3 YO :',A2YO - A2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 3 YC :',A2YC - A2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 3 Tot:',A2TOT - A2TM
         A2TMZO = A2ZO
         A2TMZC = A2ZC
         A2TMYO = A2YO
         A2TMYC = A2YC
         A2TM   = A2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 3 for A2 term'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
4000  CONTINUE
      IF (A2TEST) THEN
         WRITE(LUPRI,'(/A,F24.8)')' A2TEST Final result:  ZO', A2ZO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:  ZC', A2ZC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:  YO', A2YO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:  YC', A2YC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:    ', A2TOT
      END IF
      IF ( IPRRSP .GT.  150 ) THEN
         WRITE(LUPRI,'(//A)') ' Gradient before swapping in A2DRV'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Swap the final result to conform with the notation for A[2]
C
      CALL DSWAP(MZVAR(IGRSYM), A2TRS,1, A2TRS(MZVAR(IGRSYM)+1),1)
C
      IF ( IPRRSP .GT.  100 ) THEN
         WRITE(LUPRI,'(//A)') ' Gradient after swapping in A2DRV'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck rsp1gr */
      SUBROUTINE RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,
     *                  ISYMV,OTRS,
     *                  VEC,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,OPMAT,
     *                  XINDX,MJWOP,WRK,LWRK,LORB,LCON,LREFST)
C
C     Compute the gradient resulting from multiplying the S matrix
C     with one ore more vectors.
C
#include "implicit.h"
#include "infdim.h"
#include "inforb.h"
#include "priunit.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "rspprp.h"
#include "infhyp.h"
#include "infvar.h"
#include "qrinf.h"
#include "infpri.h"
C
      LOGICAL LORB, LCON, LREFST
C
      DIMENSION OTRS (KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION VEC(NZYVEC)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION XINDX(*)
      DIMENSION WRK(*)
C
      IF ( IPRRSP .GT. 150 ) THEN
         WRITE(LUPRI,'(A)') ' Vector in RSP1GR'
         IF (LREFST) WRITE(LUPRI,'(A)') ' (Reference state)'
         CALL OUTPUT(VEC,1,NZYVEC,1,NSIM,NZYVEC,NSIM,1,LUPRI)
         IF ( LORB ) THEN
            WRITE(LUPRI,'(//A)') ' Density matrix in RSP1GR'
            CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
         WRITE(LUPRI,'(//A)') ' One electron matrix in RSP1GR'
         CALL OUTPUT(OPMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IF ( LORB ) THEN
         TRPLET = IGRSPI.NE.ISPIN
         CALL ORBSX(NSIM,IGRSYM,KZYVR,OTRS,OPMAT,OVLAP,
     *              ISYMDN,DEN1,MJWOP,WRK,LWRK)
      END IF
C
      IF ( LCON ) THEN
         ISYMJ  = MULD2H( IGRSYM, IREFSY )
         NZCSTJ = MZCONF( IGRSYM )
C
         TRPLET = ISPIN.EQ.1
         CALL CONSX(NSIM,KZYVR,IGRSYM,OPMAT,VEC,NZYVEC,
     *              NZCVEC,ISYMV,NZCSTJ,ISYMJ,LREFST,OTRS,XINDX,
     *              WRK,LWRK)
      END IF
C
      RETURN
      END
C  /* Deck orbsx */
      SUBROUTINE ORBSX(NSIM,IGRSYM,KZYVR,OTRS,OPMAT,OVLAP,ISYMDN,DEN1,
     *                 MJWOP,WRK,LWRK)
C
C     Pupose:
C     Get the gradient vector from the zymat one -electron matrix.
C     We have the expression:
C
C     k(l,k) = <0| [ E(k,l), zymat ] |0>
C
C    We have the expressions:
C
C    1) k(a,i) = 2 * OVLAP * zymat(a,i)
C
C    2) k(i,a) = - 2 * OVLAP * zymat(i,a)
C
C    3) k(t,i) = 2 * OVLAP * zymat(t,i)- sum(x) zymat(x,i) D(x,t)
C
C    4) k(i,t) = sum(x) zymat(i,x) D(t,x) - 2 * OVLAP * zymat(i,t)
C
C    5) k(t,a) = - sum(x) zymat(x,a) D(t,x)
C
C    6) k(a,t) = sum(x) zymat(a,x) D(t,x)
C
C    7) k(u,t) = sum(x) zymat(u,x) D(t,x) - zymat(x,t) D(x,u)
C
#include "implicit.h"
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "infdim.h"
#include "infpri.h"
#include "infrsp.h"
#include "wrkrsp.h"
#include "qrinf.h"
C
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, DM1 = -1.0D0 )
C
      DIMENSION OTRS(KZYVR)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION OPMAT(NORBT,NORBT), MJWOP(2,MAXWOP,8)
      DIMENSION WRK(*)
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(//A)') ' One electron matrix in ORBSX'
         CALL OUTPUT(OPMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(//A)')
     *             ' Sigma vector before one electron operator'
         CALL OUTPUT(OTRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Compute offset for the k(k,l) elements
C
      NZCONF = MZCONF(IGRSYM)
      NYCONF = MZCONF(IGRSYM) + MZVAR(IGRSYM)
C
C     distribute one electron matrix in the gradient:
C
         KSYM1 = 0
         DO 1300 IG = 1,MZWOPT(IGRSYM)
            K     = MJWOP(1,IG,IGRSYM)
            L     = MJWOP(2,IG,IGRSYM)
            KSYM  = ISMO(K)
            LSYM  = ISMO(L)
            IF( KSYM.NE.KSYM1 ) THEN
               KSYM1 = KSYM
               IORBK = IORB(KSYM)
               NASHK = NASH(KSYM)
               NISHK = NISH(KSYM)
               IASHK = IASH(KSYM)
               IIORBK= IIORB(KSYM)
               KXSYM = MULD2H(KSYM,ISYMDN)
               NORBKX= NORB(KXSYM)
               IORBKX= IORB(KXSYM)
               NASHKX= NASH(KXSYM)
               NISHKX= NISH(KXSYM)
               IASHKX= IASH(KXSYM)
               IORBL = IORB(LSYM)
               NASHL = NASH(LSYM)
               NISHL = NISH(LSYM)
               IASHL = IASH(LSYM)
               IIORBL= IIORB(LSYM)
               LXSYM = MULD2H(LSYM,ISYMDN)
               NORBLX= NORB(LXSYM)
               IORBLX= IORB(LXSYM)
               NASHLX= NASH(LXSYM)
               NISHLX= NISH(LXSYM)
               IASHLX= IASH(LXSYM)
            END IF
            ITYPK = IOBTYP(K)
            ITYPL = IOBTYP(L)
            IF ( ITYPK.EQ.JTINAC )THEN
C
C           Do k(a,i) and k(i,a)
C
               IF (.NOT.TRPLET) THEN
                  OTRS(NZCONF+IG) = OTRS(NZCONF+IG)
     *                 + D2 * OVLAP * OPMAT(L,K)
                  OTRS(NYCONF+IG) = OTRS(NYCONF+IG)
     *                 - D2 * OVLAP * OPMAT(K,L)
               END IF
               IF ( ITYPL.EQ.JTACT .AND. NASHLX .GT. 0) THEN
C
C              Do k(t,i) and k(i,t)
C
                  NWL = ISW(L) - NISHT
                     OTRS(NZCONF+IG) = OTRS(NZCONF+IG)
     *               - DDOT(NASHLX,OPMAT(IORBLX+NISHLX+1,K),1,
     *                                  DEN1(IASHLX+1,NWL),1)
                  DO 730 IX = 1,NASHLX
                     OTRS(NYCONF+IG) = OTRS(NYCONF+IG)
     *               + OPMAT(K,IORBLX+NISHLX+IX)*
     *                                  DEN1(NWL,IASHLX+IX)
 730              CONTINUE
               END IF
            ELSE
C
C           We now know that K labels an active orbital
C
               IF (ITYPL.EQ.JTACT) THEN
C
C                 Do k(t,u)
C
                  NWL = ISW(L) - NISHT
                  NWK = ISW(K) - NISHT
                  IF (NASHLX .GT. 0) THEN
                     OTRS(NZCONF+IG) = OTRS(NZCONF+IG)
     *               -DDOT(NASHLX,OPMAT(IORBLX+NISHLX+1,K),1,
     *                                     DEN1(IASHLX+1,NWL),1)
                     DO 740 IX = 1,NASHLX
                        OTRS(NYCONF+IG) = OTRS(NYCONF+IG)
     *                    +OPMAT(K,IORBLX+NISHLX+IX)*DEN1(NWL,IASHLX+IX)
 740                 CONTINUE
                  END IF
                  IF (NASHKX .GT. 0) THEN
                     OTRS(NYCONF+IG) = OTRS(NYCONF+IG)
     *                  -DDOT(NASHKX,OPMAT(IORBKX+NISHKX+1,L),1,
     *                                     DEN1(IASHKX+1,NWK),1)
                     DO 750 IX = 1,NASHKX
                        OTRS(NZCONF+IG) = OTRS(NZCONF+IG)
     *                    +OPMAT(L,IORBKX+NISHKX+IX)*DEN1(NWK,IASHKX+IX)
 750                 CONTINUE
                  END IF
               ELSE IF (NASHKX .GT. 0) THEN
C
C                 Do k(a,t) and k(t,a)
C
                  NWK = ISW(K) - NISHT
                  OTRS(NYCONF+IG) = OTRS(NYCONF+IG)
     *                  -DDOT(NASHKX,OPMAT(IORBKX+NISHKX+1,L),1,
     *                                     DEN1(IASHKX+1,NWK),1)
                  DO 760 IX = 1,NASHKX
                     OTRS(NZCONF+IG) = OTRS(NZCONF+IG)
     *                  +OPMAT(L,IORBKX+NISHKX+IX)*DEN1(NWK,IASHKX+IX)
 760              CONTINUE
               ENDIF
            ENDIF
 1300    CONTINUE
C
      IF (IPRRSP.GT.100) THEN
         WRITE(LUPRI,'(//A)') ' Sigma vector from one electron operator'
         CALL OUTPUT(OTRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck consx */
      SUBROUTINE CONSX(NSIM,KZYVR,IGRSYM,OPMAT,VEC,NZYVEC,NZCVEC,
     *                 ISYMV,NZCSTJ,ISYMJ,LREFST,OTRS,XINDX,WRK,LWRK)
C
C     Purpose:
C     Perform the configuration part of the linear transformation.
C     Expressions:
C
C           < j    | zymat | 0(R) >
C         - < 0(L) | zymat | j    >
C
C
C     The expression of the operator zymat between two arbitrary
C     CSF's, in terms of zymat, is:
C
C     <CSF1 | K | CSF2 > = sum(xy)  zymat(x,y) D(x,y)
C
#include "implicit.h"
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, D2 = 2.0D0 )
C
#include "maxorb.h"
#include "maxash.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "infopt.h"
#include "inforb.h"
#include "infind.h"
#include "infpri.h"
#include "rspprp.h"
#include "infhyp.h"
#include "infvar.h"
#include "qrinf.h"
#include "cbgetdis.h"
C
      DIMENSION OPMAT(NORBT,NORBT), VEC(NZYVEC)
      DIMENSION WRK(*), OTRS(KZYVR), XINDX(*)
C
      LOGICAL NOH2,IH8SM,LREFST, ITEST
      DATA ITEST /.FALSE./
C
C     Allocate work space for copying inactive part of one
C     electron matrix
C
      KZYAC  = 1
      KFREE  = KZYAC + N2ASHX
      LFREE  = LWRK - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('CONSX',KFREE-1,LWRK)
C
      NOH2 = .TRUE.
      IH8SM = .TRUE.
      IF (TRPLET) THEN
         ISPIN1 = 1
         ISPIN2 = 0
      ELSE
         ISPIN1 = 0
         ISPIN2 = 0
      END IF
      IF ( LREFST ) THEN
         ISYMCI = IREFSY
      ELSE
         ISYMCI = MULD2H(IREFSY,ISYMV)
      END IF
C
      IF( IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A/)') ' Vector in CONSX'
         CALL OUTPUT(VEC,1,NZYVEC,1,NSIM,NZYVEC,NSIM,1,LUPRI)
      END IF
C
C
      IF( IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A/)') ' One electron matrix in CONSX'
         CALL OUTPUT(OPMAT,1,NORBT,1,NORBT,NORBT,NORBT,
     *   1,LUPRI)
      END IF
C
C    Do <j | A | CR>
C    ===============
C
C    Copy  active-active blocks of one electron matrix to WRK
C    Put the sign of CR into the one electron matrix if .NOT. LREFST
C
      IF (. NOT. LREFST ) THEN
         SIGN = DM1
      ELSE
         SIGN = D1
      END IF
C
      OPINAC = D0
      IF (.NOT.TRPLET) THEN
         DO 105 IW = 1,NISHT
            IX = ISX(IW)
            OPINAC = OPINAC + OPMAT(IX,IX)
  105    CONTINUE
         OPINAC = SIGN * D2 * OPINAC
      END IF
C
C     We need transposed matrix because
C     CISIGD calculates <0|A|j>, and <j|A|0> = <0|At|j>
C
      DO 110 IW = 1, NASHT
         IX = ISX(NISHT + IW)
         DO 111 JW = 1, NASHT
            JX = ISX(NISHT + JW)
            IND = (IW-1) * NASHT + JW + KZYAC - 1
            WRK(IND) = SIGN * OPMAT(IX,JX)
111      CONTINUE
110   CONTINUE
C
      IF( IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A/)')
     *          ' Active part of one electron matrix'
         CALL OUTPUT(WRK(KZYAC),1,NASHT,1,NASHT,NASHT,NASHT,
     *   1,LUPRI)
      END IF
C
      CALL CISIGD(ISYMCI,ISYMJ,NZCVEC,NZCSTJ,VEC,
     *               OTRS,WRK(KZYAC),DUMMY,
     *               NOH2,IH8SM,XINDX,ISPIN1,ISPIN2,WRK(KFREE),
     *               LFREE)
C
      IF (ISYMCI .EQ. ISYMJ .AND..NOT.TRPLET) THEN
         IF ( NISHT .GT. 0 )
     *   CALL DAXPY(NZCVEC,OPINAC,VEC,1,
     *              OTRS,1)
      END IF
C
      IF( ITEST .AND. .NOT.TRPLET ) THEN
         IF ( LREFST .AND. (ISYMCI .EQ .ISYMJ) ) THEN
            WRITE(LUPRI,'(//A,/A)')
     *          ' Sigma vector after Z part of linear transformation',
     *          ' Reference state included'
            CALL OUTPUT(OTRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
            T1 = DDOT(NZCVEC,VEC,1,OTRS,1)
            WRITE(LUPRI,'(//A)') ' Reference state projected out'
            WRITE(LUPRI,'(A,F14.7)') ' Average value T1 =', T1
            CALL DAXPY(NZCVEC,(-T1),VEC,1,OTRS,1)
         END IF
      END IF
C
      IF ( IPRRSP .GT. 200 ) THEN
         WRITE(LUPRI,'(//A)')
     *         ' Sigma vector after Z part of linear transformation'
         CALL OUTPUT(OTRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C    Do -<CL | A | j>
C    ===============
C
C    Copy  active-active blocks of one electron  matrix to WRK
C    Put the minus sign into the operator
C
      IF ( LREFST ) OPINAC = -OPINAC
      DO 210 IW = 1, NASHT
         IX = ISX(NISHT + IW)
         DO 211 JW = 1, NASHT
            JX = ISX(NISHT + JW)
            IND = (IW-1) * NASHT + JW + KZYAC - 1
            WRK(IND) = DM1 * OPMAT(JX,IX)
211      CONTINUE
210   CONTINUE
C
      IF( IPRRSP .GT. 200) THEN
         WRITE(LUPRI,'(//A/)')
     *             ' Transposed active part of one electron matrix'
         CALL OUTPUT(WRK(KZYAC),1,NASHT,1,NASHT,NASHT,NASHT,
     *      1,LUPRI)
      END IF
C
      IF( LREFST ) THEN
         ICOFF = 1
      ELSE
         ICOFF = MZCONF(ISYMV) + MZWOPT(ISYMV) + 1
      ENDIF
      ISOFF = MZCONF(IGRSYM) + MZWOPT(IGRSYM) + 1
C
      CALL CISIGD(ISYMCI,ISYMJ,NZCVEC,NZCSTJ,VEC(ICOFF),
     *            OTRS(ISOFF),WRK(KZYAC),DUMMY,
     *            NOH2,IH8SM,XINDX,ISPIN1,ISPIN2,WRK(KFREE),
     *            LFREE)
C
      IF (ISYMCI .EQ. ISYMJ .AND. .NOT.TRPLET) THEN
         IF ( NISHT .GT. 0 )
     *      CALL DAXPY(NZCVEC,OPINAC,VEC(ICOFF),1,
     *                 OTRS(ISOFF),1)
      END IF
C
      IF( ITEST .AND. .NOT.TRPLET ) THEN
         IF ((IREFSY .EQ .ISYMJ).AND. LREFST ) THEN
            WRITE(LUPRI,'(//A,/A)')
     *          ' Sigma vector after Y part of linear transformation',
     *          ' Reference state included'
            CALL OUTPUT(OTRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
            T1 = DDOT(NZCVEC,VEC(ICOFF),1,
     *                   OTRS(ISOFF),1)
            WRITE(LUPRI,'(//A)') ' Reference state projected out'
            WRITE(LUPRI,'(A,F14.7)') ' Average value T1 =', T1
            CALL DAXPY(NZCVEC,(-T1),VEC(ICOFF),1,
     *                    OTRS(ISOFF),1)
         END IF
      END IF
C
      IF ( IPRRSP .GT. 200 ) THEN
         WRITE(LUPRI,'(//A)')
     *        ' Sigma vector after Y part of linear transformation'
         CALL OUTPUT(OTRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck rsph1 */
      SUBROUTINE RSPH1(H1,CMO,WRK,LWRK)
C
C     Get the one-index transformed one electron integrals and put
C     them in H1
C
#include "implicit.h"
      DIMENSION H1(NORBT,NORBT), CMO(NCMOT), WRK(LWRK)
C
C Used from common blocks:
C  INFORB : NORBT,NBAST,...
C  INFRSP : IPRRSP
C
#include "priunit.h"
#include "inforb.h"
#include "infrsp.h"
C
      LNEED = N2BASX + NNBAST
      IF (LNEED .GT. LWRK) THEN
         CALL QUIT('Insufficient WRK space in RSPH1')
      END IF
      KH1AO = 1
      KWRK1 = KH1AO + N2BASX
      LWRK1 = LWRK + 1 - KWRK1
      IF (LWRK1.LT.0) CALL ERRWRK('RSPH1',KWRK1-1,LWRK)
C     use SIRH1, not RDONEL, in order get field dependent terms added
C     CALL RDONEL('ONEHAMIL',.TRUE.,WRK,NNBAST)
      CALL SIRH1(WRK(KH1AO),WRK(KWRK1),LWRK1)
      CALL PKSYM1(WRK(KWRK1),WRK,NBAS,NSYM,-1)
      CALL DSPTSI(NBAST,WRK(KWRK1),WRK(KH1AO))
C
      CALL DZERO(H1,N2ORBX)
      DO 600 ISYM=1,NSYM
         IF(NORB(ISYM).LE.0)GO TO 600
         CALL UTHV(CMO(ICMO(ISYM)+1),WRK(KH1AO),CMO(ICMO(ISYM)+1),
     &             ISYM,ISYM,NBAS(ISYM),NBAS(ISYM),
     &             H1,WRK(KWRK1))
C        CALL UTHV(U,PRPAO,V,ISYM,JSYM,NBASI,NBASJ,PRPMO,WRK)
 600  CONTINUE
      IF (IPRRSP.GT.1000) THEN
         WRITE(LUPRI,'(/A)')' One-electron Hamiltonian in MO basis'
         CALL OUTPUT(H1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      RETURN
      END
      SUBROUTINE QRGP(WORD,GP,CMO,XINDX,ANTSYM_GP,WORK,LWORK)
#include "implicit.h"
      CHARACTER*8 WORD
      DIMENSION GP(*), CMO(*), XINDX(*), WORK(LWORK)
#include "inforb.h"
#include "infrsp.h"
#include "infvar.h"
#include "rspprp.h"
#include "inflr.h"
#include "wrkrsp.h"
#include "priunit.h"
#include "qrinf.h"
#include "inflin.h"
#include "infrank.h"
      LOGICAL LREFST, LORB, LCON, TRPSAVE
      PARAMETER (D1     = 1.0D0)
C
      CALL QENTER('QRGP')
C
      KWORK = 1
      KFREE = KWORK
      LFREE = LWORK
C
C Operator
C
      CALL MEMGET('REAL',KOP,N2ORBX,WORK,KFREE,LFREE)
      KSYMP=KSYMOP
      IPR=IPRRSP
      CALL PRPGET(WORD,CMO,WORK(KOP),KSYMP,ANTSYM,WORK(KFREE),LFREE,IPR)
      IOPSYM=KSYMOP
      IOPSPI=OPRANK(INDPRP(WORD))
C
C     Whereas ANTSYM from PRPGET refers to matrix symmetry, ANTSYM_GP will
C     from here on refer to the antisymmetry of the response property vector,
C     thus we change sign on it.
C
      ANTSYM_GP = -ANTSYM
C
C Density 
C
      IGRSYM=KSYMOP
      IF (TRPLET) THEN
         IGRSPI=1
      ELSE
         IGRSPI=0
      END IF
      CALL MEMGET('REAL',KD,N2ASHX,WORK,KFREE,LFREE)
      IF (TDHF) THEN
         CALL DUNIT(WORK(KD),NASHT)
         CALL MEMGET('REAL',KCREF,0,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KCREF,MZCONF(1),WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KP,1,WORK,KFREE,LFREE)
         CALL GETREF(WORK(KCREF),MZCONF(1))
         CALL RSPDM(IREFSY,IREFSY,MZCONF(1),MZCONF(1),
     &      WORK(KCREF),WORK(KCREF),WORK(KD),WORK(KP),
     &      IGRSPI,IOPSPI,.FALSE.,.TRUE.,XINDX,WORK,KFREE,LFREE)
      END IF
C
C Gradient
C
      LGP    = KZYVAR
      IGPSYM = IGRSYM
      IDSYM  = 1
      LREFST = .TRUE.
      LORB   = KZWOPT.GT.0
      IF (IGRSYM.EQ.1) THEN
         LCON   = KZCONF.GT.1
      ELSE
         LCON   = KZCONF.GT.0
      END IF
      CALL MEMGET('INTE',KMJWOP,2*MAXWOP*8,WORK,KFREE,LFREE)
      CALL SETZY(WORK(KMJWOP))
      CALL DZERO(GP,LGP)
      TRPSAVE=TRPLET
      CALL RSP1GR(1,LGP,IOPSYM,IOPSPI,IGRSYM,IGRSPI,
     &   IGPSYM,GP,
     &   WORK(KCREF),MZCONF(1),MZCONF(1),D1,IDSYM,WORK(KD),WORK(KOP),
     &   XINDX,WORK(KMJWOP),WORK(KFREE),LFREE,LORB,LCON,.TRUE.)
      TRPLET=TRPSAVE
#ifdef DEBUG
      CALL HEADER('QRGP:Operator '//WORD,0)
      CALL OUTPUT(WORK(KOP),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      CALL HEADER('QRGP:One electron density',0)
      CALL OUTPUT(WORK(KD),1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      CALL HEADER('QRGP:Property vector',0)
      CALL OUTPUT(GP,1,LGP/2,1,2,LGP/2,2,1,LUPRI)
#endif
      CALL MEMREL('QRGP',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('QRGP')
      END
