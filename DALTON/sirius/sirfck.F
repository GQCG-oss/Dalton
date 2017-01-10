#define SIRFCK_DEBUG -1
#define SIRFCK_DEBUG_HERMIT -1
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
! -- FILE: DALTON/sirius/sirfck.F --
C
!#define noSIRFCK_MANU_SRHFX 
C
C  /* Deck sirfck */
      SUBROUTINE SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,AO_DIRECT,
     &                  WRK,LWORK)
C
C     March 1997 - tsaue Screening
C     921201-hjaaj; partly based on GETFCK
C     written by Henrik Koch and Trygve Helgaker 22-November-1991.
C
C PURPOSE : Driver routine for the calculation of the two-electron
C           part of the fock matrices from AO integrals.
C           We assume the densities and fock matrices are full
C           squares and without symmtry reduction.
C           The subroutine offers the possibility of calculating
C           the fock matrices directly or by reading the integrals
C           from disk.
C
C output:
C     FMAT   :  NDMAT Fock matrices with 2-el. cont. added
C  input:
C     FMAT   :  NDMAT Fock matrices
C     DMAT   :  NDMAT density matrices
C     NDMAT  :  number of density matrices and Fock matrices
C     ISYMDM :  point group symmetry of DMAT and FMAT
C     IFCTYP = +/-XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C       + sign: alpha + beta matrix (singlet)
C       - sign: alpha - beta matrix (triplet)
C
C     AO_DIRECT = .TRUE.   Direct calculation using hermit.
C                 .FALSE.  Read integrals from AOTWOINT/AOSR2INT
Cholesky       Overridden by the CHOINT flag
C
#include "implicit.h"
      LOGICAL   AO_DIRECT, DOERG, DOATR, LSRHYBR
      LOGICAL   DO_SPINDNS, DFT_SPINDNS_SAVE
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,*)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT)
      DIMENSION ISYMDM_local(NDMAT),IFCTYP_local(NDMAT) ! automatic, temporary arrays
      DIMENSION WRK(*)
      DIMENSION DFCOEFF(NBAST) !density-fitting coefficients
C
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D4 = 4.0D0, 
     &           DMP5 = -0.5D0)
C
C Used from common blocks
C   GNRINF: SRINTS, CHIVAL
C   INFORB: NSYM, NBAS(:), N2BASX
C   INFTAP: LUINTA
C   INFPRI: IPRFCK
C   CBIREA: ATOMDF, NODDYDF
C
#include "maxorb.h"
#include "priunit.h"
#include "gnrinf.h"
#include "infinp.h"
#include "inforb.h"
#include "inftap.h"
#include "infpri.h"
#include "dftcom.h"
#include "dfterg.h"
#include "cbirea.h"
C
Cholesky
#include "ccdeco.h"
#include "sirftim.h"
#include "choscf.h"
Cholesky
C

C
      LOGICAL SFIRST, PREHF, IOPTVP, DO_CHOINT
      SAVE SFIRST
      DATA SFIRST  /.TRUE./
C
C
#include "memint.h"
      CALL QENTER('SIRFCK')
      CALL GETTIM(T0,W0)
C
      iprfck = max(iprfck, SIRFCK_DEBUG)
      iprsrdft = iprfck
      IF (NDMAT.EQ.0) GO TO 9999

      IF (DODFT .OR. ADDSRI) THEN
         DOERG = DFTADD
         DOATR = .NOT.DFTADD
C        ... DFTADD tells if add normal DFT Vxc or second order Vxc
      END IF
C
C     Find out if we have spin-density terms
C     Find out if we shall and can use Cholesky
C
      DO_SPINDNS = .FALSE.
      DO_CHOINT  = CHOINT .and. CHIVAL .lt. 0.0D0 ! no Cholesky for short-range integrals (yet)
      DO I = 1,NDMAT
         IF (MOD(IFCTYP(I),10) .EQ. 0) IFCTYP(I) = 0 ! "X" does not matter if no contribution
         IF (IFCTYP(I) .LT. 0) DO_SPINDNS = .TRUE.
         IF (IFCTYP(I) .ne. 13 .and. IFCTYP(I) .ne. 12 .and.
     &       IFCTYP(I) .ne.  3 .and. IFCTYP(I) .ne.  2 .and.
     &       IFCTYP(I) .ne. 22 .and. IFCTYP(I) .ne. 11 .and.
     &       IFCTYP(I) .ne.  1 .and. IFCTYP(I) .ne.  0)
     &       DO_CHOINT = .FALSE.     ! actually all types are OK now ! /hjaaj Mar 2012
      END DO
      IF (DO_SPINDNS) THEN
         DO_CHOINT = .FALSE.
         IF (.NOT. ADDSRI) THEN
            WRITE(LUPRI,'(//A/A,20I5)')
     &      'ERROR: SIRFCK can only do spin density for srDFT',
     &      'ERROR: IFCTYP =',IFCTYP(1:NDMAT)
            CALL QUIT('SIRFCK can only do spin dens. for srDFT')
         END IF
      END IF
      IF (CHOINT .AND. .NOT. DO_CHOINT) THEN
         NINFO = NINFO + 1
         WRITE(LUPRI,'(A,(10I5))')
     &   ' INFO: IFCTYP(1:NDMAT) =',(IFCTYP(I),I=1,NDMAT)
         WRITE(LUPRI,'(/A/A//A)')
     &   ' INFO: Cholesky decomposed two-electron integrals can'//
     &   ' only be applied for IFCTYP equal to 2-3, 11-13 or 22.',
     &   ' You may have to implement the request yourself !',
     &   ' Proceeding with standard AO-integral direct '//
     &     'Fock matrix construction.'
      END IF

#if SIRFCK_DEBUG > 0
      iprsrdft = max(1,iprsrdft)  ! development code, always print
      write(lupri,*)'SIRFCK: entering; dumping qenter call tree'
      call qdump(lupri)
      write(lupri,*)'iprfck,iprsrdft =',iprfck,iprsrdft
      write(lupri,*)'SIRFCK: NDMAT=', NDMAT
      write(lupri,*)'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
      write(lupri,*)'SIRFCK: DO_SPINDNS,DFT_SPINDNS = ',
     & DO_SPINDNS,DFT_SPINDNS
#endif
      iprsrdft = iprfck
C
      IF (DO_SPINDNS) THEN
C        DC, DV, DS or DXC, DXV, DXS
         NSR_DMAT = NDMAT / 3
         NDCandDV = NSR_DMAT * 2
         IF (NSR_DMAT*3 .ne. NDMAT) THEN
C        ... check if programming was OK
            write(lupri,*)'Inconsistency(1), NSR_DMAT*3 ne NDMAT'
            write(lupri,*)'SIRFCK: NDMAT,NSR_DMAT=', NDMAT,NSR_DMAT
            write(lupri,*)'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
            CALL QUIT('Inconsistency(1), NSR_DMAT*3 ne NDMAT')
         END IF
         IF (NSR_DMAT .gt. 1) THEN
C        off-sets below assume DXC1, DXC2, DXV1, DXV2
C        but the correct order is DXC1, DXV1, DXC2, DXV2
           call quit('hjaaj: not correct off-sets in sirfck DO_SPINDNS')
C        REMEMBER ALSO zero the spin-density FMAT, because that is not done
C        in the calculation of the long-range contribution
         END IF
      ELSE IF (NASHT .EQ. 0 .OR. NDMAT .EQ. 1 .OR. .NOT. ADDSRI) THEN
C        ... NDMAT typically .eq. 1 for CI or CI-srDFT,
C        DC or DXC (closed shell SCF) / Dtot
         NSR_DMAT = NDMAT
         NDCandDV = NDMAT
      ELSE
C        DC, DV or DXC, DXV
         NSR_DMAT = NDMAT / 2
         NDCandDV = NDMAT
         IF (NSR_DMAT*2 .ne. NDMAT) THEN
            write(lupri,*)'Inconsistency(2), NSR_DMAT*2 ne NDMAT'
            write(lupri,*)'SIRFCK: NDMAT,NSR_DMAT=', NDMAT,NSR_DMAT
            write(lupri,*)'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
            CALL QUIT('Inconsistency(2), NSR_DMAT*2 ne NDMAT')
         END IF
      END IF
      IF (ADDSRI) THEN
         IF (DOERG .AND. NSR_DMAT .GT. 1) THEN
            CALL QUIT('SIRFCK ERROR: nsr_dmat > 1 for sr gradient call')
         END IF
      END IF

C     ...PREHF: Find out if this is a HF calculation preceding
C        a Multiconfigurational SR-DFT calculation.
      PREHF = DOSCF.AND.(DOMCSRDFT.OR.DOCISRDFT).AND.
     &                  (.NOT.DOHFSRDFT).AND.
     &                  (MCTYPE.EQ.0.OR.MCTYPE.EQ.-1)
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT .AND. .NOT.PREHF
C
#if SIRFCK_DEBUG > 1
CFRAN 17/02/2010
      WRITE(LUPRI,*) 'PREHF   = ', PREHF
      WRITE(LUPRI,*) 'LSRHYBR = ', LSRHYBR
CFRAN 17/02/2010 END
#endif
      HFXFACSAVE = HFXFAC
      CHIVALSAVE = CHIVAL
      IF (DO_CHOINT) THEN
         IF (IPRFCK .GT. 0) THEN
            WRITE(LUPRI,'(/A,I4,A/A,(20I5))')
     &      'SIRFCK: Cholesky construction of',NDMAT,' Fock matrices',
     &      '        IFCTYP:',(IFCTYP(I),I=1,NDMAT)
         END IF
         DTIME  = SECOND()
         ICOUN3 = 0
C
C        *************** Cholesky section *************** 
C
C        ... Outermost loop over density matrices in order
C            to maximize work space for batches of Cholesky
C            vectors. We may want to change this to avoid
C            reading through the Cholesky file(s) NDMAT times.
C
         DO I = 1,NDMAT
C
           JSYMDM = ISYMDM(I)
           IFCTYP_CHO = IFCTYP(I)

           HFXFAC_SAVE = HFXFAC
           IF (MOD(IFCTYP_CHO,10) .EQ. 1) THEN       ! No exchange
              ITSLDA = .TRUE. ! flag in Cholesky for no exchange
              HFXFAC = 0.0D0
           ELSE
              ITSLDA = .FALSE.
           END IF
C
C          hjaaj note 120327:
C            NOTE that CHO_FCKDS1 normally does not use the DMAT in input for the exchange part,
C            but rather the wave function density from CMO
C            (i.e. not DXCAO, DXVAO, DTVAO even if they are given in input!!!)
C
           IF (IFCTYP_CHO .EQ. 13 .OR. IFCTYP_CHO .EQ. 11) THEN
               IOPTVP = .TRUE.                ! Coulomb included
C           ... singlet, (tot.) symmetric density matrices
               CALL CHO_FCKDS1(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         JSYMDM,IOPTVP)
C
           ELSE IF (IFCTYP_CHO .EQ. 12) THEN
               IOPTVP = .FALSE.
C           ... singlet, (tot.) symmetric density matrices (exchange only)
               CALL CHO_FCKDS1(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         JSYMDM,IOPTVP)
C
           ELSE IF (IFCTYP_CHO .EQ.  3 .OR. IFCTYP_CHO .EQ. 1) THEN
               IOPTVP = .TRUE.
C           ...  singlet, non-symmetric (general) density matrices
               ICOUN3 = ICOUN3 + 1
               CALL CHO_FCKDS3(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         JSYMDM,IOPTVP)
C
           ELSE IF ((IFCTYP_CHO .EQ. 2) .OR. (IFCTYP_CHO .EQ.  22)) THEN
               IOPTVP = .FALSE.
C           ...  singlet, non-symmetric density matrices (exchange only)
               ICOUN3 = ICOUN3 + 1
               CALL CHO_FCKDS3(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         JSYMDM,IOPTVP)
C
           ELSE IF (IFCTYP_CHO .NE.  0) THEN
               GO TO 8000
           END IF
C
           HFXFAC = HFXFAC_SAVE

         ENDDO ! I = 1,NDMAT
         DTIME = SECOND() - DTIME
         TSIRF = TSIRF    + DTIME
         IF (IPRFCK .GT. 0 .AND. ICOUN3 .EQ. NDMAT) THEN
            WRITE(LUPRI,'(/A,I3,A,F10.2,A,/)')
     &     ' --- SIRFCK Cholesky: ',NDMAT,
     &     ' non-symmetric Fock matrices calculated in',DTIME,' seconds'
         ENDIF

      ELSE IF (AO_DIRECT) THEN ! not Cholesky, AO-direct?
C
C
C     ****************** AO direct Fock matrices ***************
C
C
         IF (IPRFCK .GT. 0) THEN
            WRITE(LUPRI,'(/A,I4,A/A,(20I5))')
     &      'SIRFCK: AO-direct construction of',NDMAT,' Fock matrices',
     &      '        IFCTYP:',(IFCTYP(I),I=1,NDMAT)
            flush(lupri)
         END IF
         IF (SFIRST) THEN
            IPRHER = 2*IPRFCK
            CALL SETHER(IPRHER,.FALSE.,WRK,LWORK)
C           CALL SETHER(IPRHER,NEWGEO,WORK,LWORK)
            SFIRST = .FALSE.
         END IF
         DO I = 1,NDMAT
            ISYMDM_local(I) = ISYMDM(I) - 1 ! 1:8 -> 0:7 range for Hermit
         END DO
         IPRHER = MAX(IPRFCK/100, SIRFCK_DEBUG_HERMIT/10)
C
C        Special code for short range/long range integrals
C        ... never short-range integrals in first call of HERFCK
C        ... if normal HF/MCSCF then normal 2-el ints (CHIVAL=0)
         SRINTS = .FALSE.
         IFCTYP_local(1:NDMAT) = IFCTYP(1:NDMAT)
         IF (ADDSRI .OR. PREHF) THEN
C        ... full exchange for long range part and for PREHF
            HFXFAC = D1
            DO I = 1,NDCandDV
            IF (IFCTYP_local(I) .NE. 0) THEN ! skip "FTC" with IFCTYP=0
               IFC_X = IFCTYP_local(I) / 10
               IFCTYP_local(I) = 10*IFC_X + 3
C              ... make sure the "do exchange" bit is set in IFCTYP
            END IF
            END DO
         END IF
         IF (PREHF) THEN
C           this is e.g. PREHF, i.e. standard HF
            CHIVAL = -1.0D0
         END IF
#if SIRFCK_DEBUG > 0
         write(lupri,*) 'sirfck: before calling herfck ...'
         write(lupri,*) 'SIRFCK: ADDSRI, SRINTS =', ADDSRI,SRINTS
         write(lupri,*) 'SIRFCK: calling HERFCK with CHIVAL ',CHIVAL
         write(lupri,*) 'SIRFCK: calling HERFCK with HFXFAC ',HFXFAC
         write(lupri,*) 'SIRFCK: ISYMDM_local', ISYMDM_local(1:NDMAT)
         write(lupri,*) 'SIRFCK: IFCTYP_local', IFCTYP_local(1:NDMAT)
#if SIRFCK_DEBUG > 4
         write(lupri,*) 'SIRFCK: first DMAT'
         call output(DMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
         flush(lupri)
#endif
C
       IF (CHIVAL .ne. D0) THEN
C        ... No long-range integrals if CHIVAL .eq. 0

         LU_DMAT = -1
         IF (NSYM.GT.1 .AND. (DODFT.OR.DOHFSRDFT.OR.DOMCSRDFT)) THEN
C           HERFCK destroys DMAT when called with symmetry
C           and DMAT needed below for short-range and DFT contributions
            CALL GPOPEN(LU_DMAT,' ','UNKNOWN',' ','UNFORMATTED',
     &         IDUMMY,.FALSE.)
            REWIND (LU_DMAT)
            CALL WRITT(LU_DMAT,NDCandDV*N2BASX,DMAT)
            REWIND (LU_DMAT)
         END IF

         if (atomdf .AND. CHIVAL .LT. D0) then
C          hjaaj mar 2006: dens_fit not for short range integrals (yet)
#if SIRFCK_DEBUG > 0
           write(lupri,*) 'calling dens_fit'; flush(lupri)
#endif
           call dens_fit(FMAT,DMAT,NDCandDV,ISYMDM_local,IFCTYP_local,
     &                   IPRHER,WRK,LWORK)
         else
#if SIRFCK_DEBUG > 0
           write(lupri,*) 'calling HERFCK, IPRHER=',IPRHER; flush(lupri)
#endif
           CALL HERFCK  (FMAT,DMAT,NDCandDV,ISYMDM_local,IFCTYP_local,
     &                   IPRHER,WRK,LWORK)
         end if

         IF (LU_DMAT .NE. -1) THEN
            CALL READT(LU_DMAT,NDCandDV*N2BASX,DMAT)
            CALL GPCLOSE(LU_DMAT,'DELETE')
         END IF

       END IF ! (CHIVAL .ne. D0)
#if SIRFCK_DEBUG > 4
       do i = 1,NDCandDV
         write(lupri,*) 'SIRFCK: FMAT no.',i,' after long-range'
         call output(FMAT(1,i),1,nbast,1,nbast,nbast,nbast,-1,lupri)
       end do
#endif
         HFXFAC = HFXFACSAVE
         CHIVAL = CHIVALSAVE
C
C      ***********************************************************
C        Add Short-range integrals (HF/MC-SRDFT hybrid model)
C      ***********************************************************
C
#if SIRFCK_DEBUG > 0
         write(lupri,*) 'add direct sr integrals ? ADDSRI = ', ADDSRI 
#endif 
         IF (ADDSRI .AND. (DOHFSRDFT.OR.DOMCSRDFT)) THEN

            SRINTS = .TRUE.

            NFSRTYP = 1
            NFSRAO = NSR_DMAT ! FT
            IFC_X = IFCTYP(NSR_DMAT) / 10 ! IFCTYP(1) might be 0 for "FTC", thus we go for the FXC/FTV value
               ! as representative for all 1:NSR_DMAT matrices
            IFC_X = 10*IFC_X ! DMAT symmetric, antisymmetric or unsymmetric ?
            IF (HFXFAC .EQ. 0.0D0) THEN
               IFCTYP_local(1:NSR_DMAT) = IFC_X + 1 ! only Coulomb
            ELSE
               IFCTYP_local(1:NSR_DMAT) = IFC_X + 3 ! Coulomb + HFXFAC*exchange
            END IF
            IF (DO_SPINDNS) THEN
               IF (HFXFAC .NE. D0) THEN
                  NFSRAO = NFSRAO + NSR_DMAT ! also FS
                  IFCTYP_local(NSR_DMAT+1:NFSRAO) = IFC_X + 2 ! only exchange
                  NFSRTYP = NFSRTYP + 1
               END IF
               JFMAT = 3
            ELSE IF (NASHT .GT. 0 .AND. (.NOT.DOCISRDFT)) THEN
               JFMAT = 2
            ELSE
C           ... total DMAT = DC when NASHT.eq.0 or CISRDFT
               JFMAT = 1
            END IF
            ISYMDM_local(1:NFSRAO) = ISYMDM(1) - 1
            KFSRAO = 1
            KDSRAO = KFSRAO + NFSRAO*N2BASX
            KLAST  = KDSRAO + NFSRAO*N2BASX
            LWRK   = LWORK - KLAST + 1
            IF (KLAST.GT.LWORK)
     &         CALL STOPIT('SIRFCK','FSRAO+DSRAO',KLAST,LWORK)

            IF (NASHT.GT.0 .AND. (.NOT.DOCISRDFT)) THEN
C              ... get total DMAT = DT = DC + DV if NASHT.gt.0
C                  for Hartree SR term and - only if singlet - HFXFAC*(exact exchange SR term)
               CALL DCOPY(NSR_DMAT*N2BASX,DMAT,1,WRK(KDSRAO),1)
               CALL DAXPY(NSR_DMAT*N2BASX,D1,DMAT(1,NSR_DMAT+1),1,   ! DT = DC + DV
     &                    WRK(KDSRAO),1)
               IF ( NFSRTYP .EQ. 2 ) THEN ! DO_SPINDNS .and. HFXFAC .ne. 0
                  KDS_SRAO = KDSRAO + NSR_DMAT*N2BASX
                  CALL DCOPY(NSR_DMAT*N2BASX,DMAT(1,2*NSR_DMAT+1),1, ! DS for exchange
     &                       WRK(KDS_SRAO),1)
               END IF
            ELSE ! NASHT.eq.0 .or. DOCISRDFT
               CALL DCOPY(NSR_DMAT*N2BASX,DMAT,1,WRK(KDSRAO),1)      ! DT = DC
            END IF

#if SIRFCK_DEBUG > 0
            write(lupri,*) 'SIRFCK: ADDSRI, SRINTS =', ADDSRI,SRINTS
            write(lupri,*) 'SIRFCK: calling HERFCK with CHIVAL ',CHIVAL
            write(lupri,*) 'SIRFCK: calling HERFCK with HFXFAC ',HFXFAC
            write(lupri,*) 'JFMAT    =', JFMAT
            write(lupri,*) 'nsr_dmat,nfsrao =', NSR_DMAT,NFSRAO
            write(lupri,*) 'ISYMDM_local',(ISYMDM_local(I),I=1,NFSRAO)
            write(lupri,*) 'IFCTYP_local',(IFCTYP_local(I),I=1,NFSRAO)
#if SIRFCK_DEBUG > 4
         do i = 1,NFSRAO
            write(lupri,*) 'SIRFCK: DMAT_SR no.',i,' out of',NFSRAO
            JDSRAO = KDSRAO + (I-1)*N2BASX
          call output(wrk(jdsrAO),1,nbast,1,nbast,nbast,nbast,-1,lupri)
         end do
#endif
#endif
!           if (atomdf .AND. CHIVAL .LT. D0) then
!               hjaaj mar 2006: dens_fit not for short range integrals (yet)
!               call dens_fit(wrk(KFSRAO),wrk(kdmat),NSR_DMAT,
!    &                      ISYMDM,IFCTYP_local,IPRHER,WRK(KLAST),LWRK)
!           else
              CALL HERFCK(WRK(KFSRAO),WRK(KDSRAO),NFSRAO,
     &             ISYMDM_local,IFCTYP_local,IPRHER,WRK(KLAST),LWRK)
!             NB! WRK(KDSRAO) may be destroyed inside HERFCK, do not use below!
!           end if

            SRINTS = .FALSE.

C
         IF (DOERG) THEN
            IF (NSR_DMAT .GT. 1)
     &         CALL QUIT('SIRFCK SRINT ERROR: DOERG and NSR_DMAT>1')
            EJSR = DP5*DDOT(N2BASX,DMAT,1,WRK(KFSRAO),1) ! Inactive part of EJsr+HFXFAC*EKsr (DMAT is DC)
            IF (NFSRTYP .EQ. 1) THEN
               EKSR = D0
               EKSR_S = D0
            ELSE IF (NFSRTYP .EQ. 2) THEN
               JFSRAO = KFSRAO + N2BASX
               EKSR_S = DP5*DDOT(N2BASX,DMAT(1,3),1,WRK(JFSRAO),1) ! EKsr_spin from 0.5 * DS * HFXFAC*Ksr(DS)
            ELSE
               CALL QUIT('NFSRTYP neither 1 nor 2')
            END IF
            IF (NASHT.GT.0) THEN
              EJKVSR = -DP5*DDOT(N2BASX,DMAT(1,2),1,WRK(KFSRAO),1) ! -0.5 * (DV * Fsr)
            ELSE
              EJKVSR =  D0
            END IF
            EJSR = EJSR - EJKVSR ! Now EJsr or EJsr + HFXFAC*Eksr from DT=DC+DV
         END IF
 
         CALL DAXPY(NSR_DMAT*N2BASX,D1,WRK(KFSRAO),1,FMAT,1) ! adding short-range Jsr+Ksr terms to FC
         IF (NFSRTYP .EQ. 2) THEN
            JFSRAO = KFSRAO + NSR_DMAT*N2BASX
            CALL DAXPY(NSR_DMAT*N2BASX,1.0d0,WRK(JFSRAO),1,FMAT(1,3),1) ! add Ksr_spin to FS
!           write (lupri,'(/A)') 'SIRFCK SPIN DENSITY EXX TERM DIRECT'
         END IF

         END IF   ! IF (ADDSRI .AND. (DOHFSRDFT.OR.DOMCSRDFT)) THEN

#if SIRFCK_DEBUG > 3
         if(addsri)then
           do i = 1,nsr_dmat
              write(lupri,*)
     &        'SIRFCK: direct short-range contribution to FMAT no.', i
           call output(WRK(KFSRAO+n2basx*(i-1)),1,nbast,1,nbast,
     &                 nbast,nbast,-1,lupri)
           end do
           write(lupri,*) 'SIRFCK: FMAT (lr + sr)'
           call output(FMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
         end if
#endif

      ELSE ! not Cholesky and not AO-direct; i.e. conventional from AOTWOINT and AOSRINT files

         IFCTYP_local(1:NDMAT) = IFCTYP(1:NDMAT)
         IF (IPRFCK .GT. 0) THEN
            WRITE(LUPRI,'(/A,I4,A/A,(20I5))')
     &      'SIRFCK: Construction of',NDMAT,
     &      ' Fock matrices from integrals on file',
     &      '        Input IFCTYP:',(IFCTYP(I),I=1,NDMAT)
         END IF
C
C     ****************** "conventional" AO Fock matrices ***************
C                        (from AO integrals on file)
C
         IF (HFXMU .NE. D0) CALL QUIT(
     &       'SIRFCK2: nonzero HFXMU requires .DIRECT calculation.') 
C
C        ... if SRHYBR (DFT-SRDFT hybrid) only include Coulomb terms 
C            (core is treated with BLYP -> hence no exchange here!)
         IF(LSRHYBR) THEN
            IFCTYP_local(1) = 11
C           ... FC,lr = JC,lr = J,lr * DC for SR-hybrid
            IF (NDMAT .EQ. 2) THEN
               IFCTYP_local(2) = 13
            END IF
C           ... FV,lr = G,lr * DV
            IF (NDMAT .GT. 2) CALL QUIT('NDMAT.gt.2 for SRHYBR & lr')
         END IF ! LSRHYBR

         IF (ADDSRI .OR. PREHF) THEN
C        ... full exchange for long range part and for PREHF
            HFXFAC = D1
            DO I = 1,NDCandDV
            IF (IFCTYP_local(I) .GT. 0) THEN ! skip "FTC"
               J = MOD(IFCTYP(I),10)
               IF (IAND(J,2) .eq. 0) IFCTYP(I) = IFCTYP(I) + 2
C              ... make sure the "do exchange" bit is set in IFCTYP
            END IF
            END DO
         END IF
#if SIRFCK_DEBUG > 0
         write(lupri,*)'SIRFCK: calling SIRFCK2 "AOTWOINT" with'
         write(lupri,*)'SIRFCK: NDMAT =',NDMAT
         write(lupri,*)'SIRFCK: IFCTYP=',(IFCTYP(I),I=1,NDMAT)
#endif
         IF (CHIVAL .ne. D0) THEN
C        ... no long-range integrals if chival .eq. 0
            CALL SIRFCK2(LUINTA,'AOTWOINT',
     &                   FMAT,DMAT,NDCandDV,ISYMDM,IFCTYP,WRK,LWORK)
         END IF ! (CHIVAL .ne. D0)

         if (lsrhybr) then
            do i = 1,ndmat
               ifctyp(i) = IFCTYP_local(i)
            end do
         end if
         HFXFAC = HFXFACSAVE
C
C        Add Short-range integrals to FC part of FMAT (SRDFT hybrid models)
C
         IF (ADDSRI) THEN
            IF (NASHT.GT.0 .AND. (.NOT.DOCISRDFT))THEN
C              ... get total DMAT = DC + DV if NASHT.gt.0
C                  for Hartree SR term
               DO I = 1,NSR_DMAT
                  CALL DAXPY(N2BASX,D1,DMAT(1,NSR_DMAT+I),1,DMAT(1,I),1)
C                 ... add DV(I) to DC(I)
               END DO
            END IF
            IF (NSR_DMAT.GT.1) CALL QUIT('NSR_DMAT .GT. 1 for ADDSRI')
C
            IF (HFXFAC.NE.D0) THEN
C           ... for normal/preceding HF and for SRDFT with exact exch.
C               (add SR coul. and exch. to LR to give full 2-el part)
                IFCTYP_local(1) = 13
            ELSE 
C           ... add only Coul
C               (only SR Hartree term for SRDFT hybrids)
                IFCTYP_local(1) = 11
            END IF
            KFSRAO = 1
            KLAST  = KFSRAO + NSR_DMAT*N2BASX
            LWRK   = LWORK  - KLAST + 1
#if SIRFCK_DEBUG > 0
            write(lupri,*)'SIRFCK: calling SIRFCK2 "AOSR2INT" with'
            write(lupri,*)'SIRFCK: NSR_DMAT,NDMAT=',NSR_DMAT,NDMAT
            write(lupri,*)'SIRFCK: IFCTYP_local=',
     &         (IFCTYP_local(I),I=1,NDMAT)
#endif
            CALL DZERO(WRK(KFSRAO),NSR_DMAT*N2BASX)
            CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFSRAO),
     &         DMAT,NSR_DMAT,ISYMDM,IFCTYP_local,WRK(KLAST),LWRK)
            CALL DAXPY(N2BASX,D1,WRK(KFSRAO),1,FMAT,1) ! add F_sr to FC_lr
            EJSR = DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KFSRAO),1) ! EJsr+HFXFAC*EKsr from DC+DV
            IF (NASHT.GT.0) THEN
               EJKVSR = -DP5*
     &            DDOT(N2BASX,DMAT(1,NSR_DMAT+1),1,WRK(KFSRAO),1) ! correction from DV
            ELSE
               EJKVSR = D0
            END IF
!            
#ifdef SIRFCK_MANU_SRHFX
! Manu: write the srHFX energy when the HFEXCH option is
! used in srDFT
!
            iprsrdft_save=iprsrdft
            iprsrdft = 1
#endif
!
            if (iprsrdft .gt. 0 .and.
     &         nsr_dmat.eq.1 .and. HFXFAC .NE. 0.0D0) then
               IFCTYP_local(1) = 12 ! exchange only
               CALL DZERO(WRK(KFSRAO),N2BASX)
               CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFSRAO),
     &            DMAT,NSR_DMAT,ISYMDM,IFCTYP_local,WRK(KLAST),LWRK)
               EKSR = DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KFSRAO),1)
               EJSR = EJSR - EKSR
            else
               EKSR = D0
            end if
           
#ifdef SIRFCK_MANU_SRHFX
            write(lupri,*) 'srHFX energy =', EKSR
            iprsrdft = iprsrdft_save ! restore
#endif
!
! Finished with short-range Hartree and normal HF exchange term,
! now for any spin density exchange term:
!
            IF (HFXFAC .NE. D0 .and. DO_SPINDNS) THEN
               write (lupri,'(/A)') 'SIRFCK2 SPIN DENSITY EXX TERM'
               IF (NSR_DMAT .NE. 1)
     &         call quit('HFXFAC & SPINDNS only for NSR_DMAT.eq.1')
               ifctyp_local(1) = 12 ! only exchange
               isymdm_local(1) = 1
               CALL DZERO(WRK(KFSRAO),N2BASX)
               CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFSRAO),
     &                   DMAT(1,3), 1 ,isymdm_local,ifctyp_local,
     &                   WRK(KLAST),LWRK)
#if SIRFCK_DEBUG > 0
               write (lupri,'(/A)')
     &         'SIRFCK SPIN DENSITY sr-DFT exact exchange Fock matrix'
           call output(wrk(kfsrao),1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
               CALL DAXPY(N2BASX,1.0D0,WRK(KFSRAO),1,FMAT(1,3),1)  ! add HFXFAC*Ksr_spin to FS
               EKSR_S = DP5*DDOT(N2BASX,DMAT(1,3),1,WRK(KFSRAO),1) ! HFXFAC*EKsr_spin from HFXFAC * DS * Ksr(DS)
            ELSE
               EKSR_S = D0
            END IF
!
            IF (NASHT.GT.0) THEN
               DO I = 1,NSR_DMAT
                  CALL DAXPY(N2BASX,-D1,DMAT(1,NSR_DMAT+I),1,
     &                                  DMAT(1,I)         ,1)
C                 ... restore DC(I) = DTOT(I) -  DV(I)
               END DO
            END IF
         END IF
      END IF ! if (do_choint) then ... else if (direct) then ... else ... end if
      HFXFAC = HFXFACSAVE
#if SIRFCK_DEBUG > 0
         write(lupri,*) 'SIRFCK: final first FMAT'
         call output(FMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
      IF (IPRFCK .GT. 0) THEN 
         CALL GETTIM(T1,W1)
         WRITE(LUPRI,'(/A,2F10.3)')
     &   ' SIRFCK: CPU and WALL time used in Fock matrix constructions',
     &   T1-T0,W1-W0
      END IF
C
#ifdef OLD_DFT_CODE
!hjaaj: can probably be deleted, but want to be sure about srdft first. March 2006.
C
C     Add DFT contribution
C
      IF (DODFT) THEN
         IF(DOATR)THEN
            J = NDMAT + 1
C           DTOTAO is put in DMAT(1,NDMAT+1) in TR1FCK.
C           We now need to copy DXAO to DMAT(1,NDMAT+2) for EXCMAT
            CALL DCOPY(N2BASX,DMAT,1,DMAT(1,J+1),1)
            NDMTMP = 2
            KSRAO  = 1
            KLAST  = KSRAO  + N2BASX
            LWRK  = LWORK - KLAST + 1
            CALL DZERO(WRK(KSRAO),N2BASX)
            CALL DFTEXC(DMAT(1,J),NDMTMP,WRK(KSRAO),1,
     &           .FALSE.,DOERG,.FALSE.,.FALSE.,DOATR,
     &           .FALSE.,DUMMY,.FALSE.,WRK(KLAST),LWRK,IPRFCK)
C           CALL HEADER('AO Fock matrix contribution in DFTATR /SIRFCK',-1)
C           CALL OUTPUT(WRK(KSRAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL DAXPY(N2BASX,D1,WRK(KSRAO),1,FMAT,1)
         ELSE
            J = 1
            NDMTMP = 1
            CALL DFTEXC(DMAT(1,J),NDMTMP,FMAT,1,.FALSE.,DOERG,.FALSE.,
     &           .FALSE.,DOATR,.FALSE.,DUMMY,.FALSE.,WRK,LWORK,IPRFCK)
         ENDIF
         IF (NDMAT.GT.1) CALL QUIT('DODFT and NDMAT.gt.1')
      END IF
#endif

C
C     ***************** ADD SR-DFT contribution if MCSCF-DFT or if HF-SRDFT
C
      IF (ADDSRI.AND.(DOHFSRDFT.OR.DOMCSRDFT).AND.(.NOT.PREHF)) THEN
#ifndef MOD_SRDFT
         call quit('srdft not included in this version')
#else
#if SIRFCK_DEBUG > 4
         do i = 1,ndmat+1
            write(lupri,*) 'SIRFCK srdft: DMAT no.',i,' out of',ndmat+1
            call output(DMAT(1,i),1,nbast,1,nbast,nbast,nbast,-1,lupri)
         end do
#endif
         IF (LSRHYBR .AND. DO_SPINDNS) THEN
            call quit('hybrid srDFT only for singlet MC-srDFT')
         END IF
         IF (DO_SPINDNS) THEN
            NSRAO = 2
            NDOFF = 3
C           3 because DC, DV, DS or DXC, DXV, DXS
         ELSE IF (NASHT.GT.0 .AND. (.NOT.DOCISRDFT))THEN
C           ... get total DMAT = DC + DV if NASHT.gt.0
C               for xc sr term
            NSRAO = 1
            NDOFF = 2
C           2 because DC, DV or DXC, DXV
         ELSE
            NSRAO = 1
            NDOFF = 1
C           1 because DC or DXC
         END IF
         IF (LSRHYBR) NSRAO = NSRAO + 1
C        ... one more EXCMAT from SRDFT when LSRHYBR
         KSRAO = 1
         KLAST = KSRAO + NSRAO*N2BASX
         IF (DOERG .AND. NSR_DMAT .GT. 1) THEN
            WRITE(LUPRI,*) 'SRFCK.SRDFT only for one SRMAT when DOERG',
     &         NSR_DMAT
            CALL QUIT('SRFCK.SRDFT only for one SRMAT when DOERG')
         END IF
         IF (LSRHYBR) THEN
            IF (NDMAT .NE. 2) THEN
               write(lupri,*) 'ERROR LSRHYBR: NDMAT.ne.2 but is ',NDMAT
               CALL QUIT('NDMAT .ne. 2 for LSRHYBR')
            END IF
            NDOFF  = 2
         END IF
         LWRK  = LWORK - KLAST + 1

         JDAO   = 0 ! initialize JD*AO to avoid compiler warnings
         JDXAO  = 0
         JDXSAO = 0
         IF (DOATR) THEN
            JDAO = NDMAT + 1
C           ... DTOTAO is in DMAT(1,NDMAT+1) from TR1FCK or RSPOLI
C           For SRHYBR we also need DVAO in DMAT(1,NDMAT+2) from TR1FCK.:
            IF (LSRHYBR) THEN
               JDXAO = JDAO + 2
            ELSE
               JDXAO = JDAO + 1
            END IF
            IF (DO_SPINDNS) THEN
C              DT, DXT, DS, DXS
               JDSAO  = JDXAO + 1
               JDXSAO = JDSAO + 1
            ELSE
C              DT, DXT
            END IF
         END IF   ! IF (DOATR)
C
!
! Manu november 2011: rewrote specific code for MCSCF-srDFT linear
!                     response when we need to compute srxc Kernel
!                     contributions for more than one root. In this case
!                     DMAT looks like
!
!  DXCAO1 DXCAO2 ... DXVAO1 DXVAO2 ... DTCAO1 DTCAO2 ... DTVAO1 DTVAO2 ... DTOT 
!  <---------------> <---------------> <---------------> <--------------->
!      NOSIM              NOSIM             NCSIM             NCSIM         1
!
!  where NOSIM > 1 or NCSIM > 1 
!
! Note: spin polarization and LSRHYBR option are not considered. Should be done
! later on ...
!
! NB! specific code for more than one root removed again April 2012 / hjaaj
!
!
C            loop over NSRMAT output matrices (later do the loop inside SRDFT for
C            efficiency!!!!!) /hjaaj June 2009
C
         DO JSRMAT = 1, NSR_DMAT
          IF (DOATR) THEN
#if SIRFCK_DEBUG > 0
              write(lupri,*) 'DOATR in sirfck: ndoff =', ndoff
              write(lupri,*) 'JSRMAT = ', JSRMAT
              write(lupri,*) 'JDXAO  = ', JDXAO
#if SIRFCK_DEBUG > 4
              WRITE(LUPRI,*) 'DMAT(1,JSRMAT) = DCMAT(JSRMAT)'
              call output(DMAT(1,JSRMAT),1,NBAST,1,NBAST,NBAST,
     &                    NBAST,-1,lupri)
              WRITE(LUPRI,*) 'DMAT(1,JSRMAT+NSRMAT) = DVMAT(JSRMAT)'
              call output(DMAT(1,JSRMAT+NSRMAT),1,NBAST,1,NBAST,NBAST,
     &                NBAST,-1,lupri)
#endif
#endif
C           we need the symmetric part of DXTOTAO in DMAT(1,JDXAO) for SRDFT :
                    IF (NDOFF .EQ. 1) THEN
C              closed shell SCF case
                       CALL DGETSI(NBAST,DMAT(1,JSRMAT),DMAT(1,JDXAO))
                    ELSE
                       CALL DCOPY(N2BASX,DMAT(1,JSRMAT),1,
     &                            DMAT(1,JDXAO),1)
                       CALL DAXPY(N2BASX,D1,DMAT(1,JSRMAT+NSR_DMAT),1,
     &                              DMAT(1,JDXAO),1)
C              ... DXTOTAO = DXCAO + DXVAO
                       CALL DGETSI(NBAST,DMAT(1,JDXAO),DMAT(1,JDXAO))
C              ... extract symmetric component (matrix is assumed
C              symmetric in SRDFT, and it will in general not have any
C              matrix symmetry e.g. when SIRFCK is called from RSPOLI)
                       IF (DO_SPINDNS) THEN
                          CALL DGETSI(NBAST,DMAT(1,JSRMAT+2*NSR_DMAT),
     &                                DMAT(1,JDXSAO))
C              ... extract symmetric component of DXSAO
                       END IF
                       IF (LSRHYBR) THEN
                          CALL DGETSI(NBAST,DMAT(1,JSRMAT+1),
     &                                DMAT(1,JDXAO+1))
C                 ... copy symmetric part of DXVAO
               END IF
            END IF
          ELSE ! (.not. DOATR)
            JDAO = JSRMAT
            IF (NDOFF.NE.1)
     &      CALL DAXPY(N2BASX,D1,DMAT(1,JDAO+1),1,DMAT(1,JDAO),1)
C           ... get total DMAT = DC + DV if NASHT .gt. 0
C           ... spin DMAT: DS is expected in JDAO + 2
          END IF
C         SRDFT(NDMAT,EXCMAT,DMAT,EDFTY,DOERG,DOGRD,DOATR,WORK,LWORK,IPRINT)
#if SIRFCK_DEBUG > 0
          WRITE(LUPRI,*) 'SIRFCK: calling SRDFT, DOERG,DOATR,NSRAO',
     &             DOERG,DOATR,NSRAO  !JT
#if SIRFCK_DEBUG > 1
          WRITE(LUPRI,*) 'DT: DMAT(1,JDAO), JDAO =',JDAO
          call output(DMAT(1,JDAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,lupri)
          IF (DOATR) THEN
             WRITE(LUPRI,*) 'DXT: DMAT(1,JDAO+1), JDAO+1 =',JDAO+1
             call output(DMAT(1,JDAO+1),1,NBAST,1,NBAST,NBAST,
     &                NBAST,-1,lupri)
                END IF
#endif
          WRITE(LUPRI,*) 'SIRFCK: DO_SPINDNS = ', DO_SPINDNS, 
     &                   'BEFORE SRDFT'
                IF (DO_SPINDNS) THEN
                 WRITE(LUPRI,*) 'DS: DMAT(1,JDAO+2), JDAO+2 =',JDAO+2
                 call output(DMAT(1,JDAO+2),1,NBAST,1,NBAST,NBAST,NBAST,
     &                      -1,lupri)
                 IF (DOATR) THEN
                  WRITE(LUPRI,*) 'DXS: DMAT(1,JDAO+3), JDAO+3 =',JDAO+3
                  call output(DMAT(1,JDAO+3),1,NBAST,1,NBAST,NBAST,
     &                        NBAST,-1,lupri)
                 END IF
                END IF
#endif
                CALL DZERO(WRK(KSRAO),NSRAO*N2BASX)
                DFT_SPINDNS_SAVE = DFT_SPINDNS
                DFT_SPINDNS = DO_SPINDNS 
                CALL SRDFT(1,WRK(KSRAO),DMAT(1,JDAO),ESRDFT,
     &                     DOERG,.FALSE.,DOATR,WRK(KLAST),LWRK,IPRSRDFT)
C                    SRDFT(NDMAT,EXCMAT,DMAT,ESRDFT,
C                          DOERG,DOGRD,DOATR,WORK,LWORK,IPRINT)
                DFT_SPINDNS = DFT_SPINDNS_SAVE
C         IF (DOATR)
C    &      ESRLTR = ESRLTR + DDOT(N2BASX,WRK(KSRAO),1,DMAT(1,JDAO),1)
C           ... for MCSRDFT csf sigma-vectors : - 2 ESRLTR cref(i)
C               NOT NEEDED because we project out CREF component in sigma vectors
C
                CALL DAXPY(N2BASX,D1,WRK(KSRAO),1,FMAT(1,JSRMAT),1) ! add VXsr/VTsr to FXCAO
                IF (DO_SPINDNS)
     &             CALL DAXPY(N2BASX,1.0d0,WRK(KSRAO+N2BASX),1,
     &                        FMAT(1,JSRMAT+2*NSR_DMAT),1)              !  add VSXsr to FXSAO

                IF (LSRHYBR) 
     &             CALL DCOPY(N2BASX,WRK(KSRAO+N2BASX),1,
     &                        FMAT(1,JSRMAT+1),1)
C
            END DO  ! JSRMAT = 1, NSR_DMAT
C
C
         IF (DOERG) THEN
            IF (NDOFF.NE.1)
     &         CALL DAXPY(N2BASX,-D1,DMAT(1,2),1,DMAT(1,1),1)
C        ... restore DC = DTOT - DV
!
!        Short range energy contributions:
!           EKSR_S =  0.5*Tr( Ksr(Ds) Ds ) * HFXFAC
!           ESRDFT = Exsr + Ecsr + EKSR_S
!           EJSR   =  0.5*Tr( Jsr(Dtot) Dtot )
!           EKSR   =  0.5*Tr( Ksr(Dtot) Dtot ) * HFXFAC
!           (EJSR and EKSR are included in EMCSCF via modified FCAC and therefore not included in ESRDFTY)
!           IF (iprsrdft .GT. 0) then
               WRITE (LUPRI,'(/A,F20.10)')
     &            '   Ex-sr + Ec-sr                  ',ESRDFT
               IF (EKSR .ne. 0.0D0) THEN
                  WRITE (LUPRI,'(A,F20.10)')
     &            ' + EJsr = sr Coulomb energy       ',EJSR,
     &            ' + EKsr = sr exchange energy      ',EKSR
               ELSE IF (HFXFAC .eq. 0.0D0) THEN
                  WRITE (LUPRI,'(A,F20.10)')
     &            ' + EJsr = sr Coulomb energy       ',EJSR
               ELSE
                  WRITE (LUPRI,'(A,F20.10)')
     &            ' + EJKsr = sr Coulomb + exchange  ',EJSR
               END IF
               IF (DO_SPINDNS .and. HFXFAC .ne. D0)
     &            WRITE (LUPRI,'(A,F20.10)')
     &            ' + EKsr_S from spin density       ',EKSR_S
               WRITE (LUPRI,'(A,F20.10/)')
     &            ' = Total E(srDFT)                 ',
     &            ESRDFT+EJSR+EKSR+EKSR_S
               flush(lupri)
!           end if
            ESRDFT  = ESRDFT + EKSR_S
            ESRDFTY = ESRDFT
            ESRDFT  = ESRDFT + EJSR + EKSR

!           ESRDFTY is ESRDFT + a correction energy needed because
!           the JKsr and Vxcsr matrices are added to FC to get correct MC-srDFt gradient.
!           However, this means FCAC will be modified and this gives a contribution to
!           the energy in CIGRAD and ENRACT, which is not the same as ESRDFT because of
!           the non-linearity and which has to be removed again to get correct energy
!           from the MCSCF code.
!           The EJsr and any EKsr from DC is calculated via the modified FC and saved in EMY.
!           In CIGRAD we get the correct JKsr(DC)*DV contribution but
!           also JKsr(DV) * DV, i.e. double counting, and we need to subtract half of it again
!           (this correction is in EJKVsr)
            ECSR = -DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KSRAO),1)
            IF (NASHT.GT.0) THEN
               EVSR = -DDOT(N2BASX,DMAT(1,2),1,WRK(KSRAO),1)
            ELSE
               EVSR = D0
            END IF
            ESRDFTY  = ESRDFTY + ECSR + EVSR + EJKVSR

C           ECSR   = -0.5*Tr( Vxcsr Dc )
C           EVSR   =     -Tr( Vxcsr Dv )
C           EJKVSR = -0.5*Tr( (Jsr+HFXFAC*Ksr) Dv)
C           note: ECSR + EVSR + EJKVSR corrections are added
C                 to get right contributions from MCSCF code.
C
!           IF (iprsrdft .GT. 0) then
!
               WRITE (LUPRI,'(/A)')
     &         ' Corrections needed for correct CI energy evaluation:'
               WRITE (LUPRI,'(A,F20.10)')
     &         ' - 0.5 Tr(Vxc-sr Dcore)           ',ECSR,
     &         ' -     Tr(Vxc-sr Dval)            ',EVSR,
     &         ' - 0.5 Tr( (Jsr+HFXFAC*Ksr) Dval )',EJKVSR
               flush(lupri)
!           END IF
         END IF ! (DOERG)
         IF (IPRFCK .GT. 0) THEN 
            CALL GETTIM(T2,W2)
            WRITE(LUPRI,'(/A,2F10.3)')
     &   ' SIRFCK: CPU and WALL time used for srDFT terms             ',
     &      T2-T1,W2-W1
         END IF
#endif  /* MOD_SRDFT */
      END IF ! (ADDSRI.AND.(DOHFSRDFT.OR.DOMCSRDFT).AND.(.NOT.PREHF))
C
      IF (NODDYDF) THEN
         IF (IFCTYP(1).NE.11) THEN
            WRITE(LUPRI,*)'IFCTYP(1) = ',IFCTYP(1)
            WRITE(LUPRI,*)'DenFit noddy code works only for Coulomb'
            CALL QUIT('Ooops... Probably trying to calculate exchange '
     &                //'with density-fitting noddy code.')
         END IF
         IF (NDMAT.GT.1) THEN
            WRITE(LUPRI,*)'NDMAT = ',NDMAT
            WRITE(LUPRI,*)'DenFit noddy code works only for NDMAT=1'
            CALL QUIT('Illegal NDMAT value in DenFit noddy code')
         END IF
C
         call dzero(FMAT(1,1),N2BASX)
C
         !debugging tool... Simulates real J_ab (no DenFit)
*        call jab_noddy(nbast,dmat(1,1),fmat(1,1))
C
         !calculate density-fitting coefficeints...
         CALL DENFIT_COEFFS_NODDY(NBAST,NBAST,DMAT(1,1),DFCOEFF)
c
         WRITE (LUPRI,'(//,1X,A,I3)') ' Noddy fitting coeffs'
         CALL OUTPUT(DFCOEFF,1,NBAST,1,1,NBAST,1,-1,LUPRI)
c
         !...and use them to build the (fitted) Fock matrix
         CALL JAB_DENFIT_NODDY(NBAST,NBAST,DFCOEFF,FMAT(1,1))
      END IF
C
#if defined (VAR_OLDCODE)
--- keep print option here ??? 921201-hjaaj
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRFCK.GT.2) THEN
         CALL HEADER('Density and Fock matrices in SIRFCK',-1)
         DO I = 1, NDMAT
            WRITE (LUPRI,'(//A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
            WRITE (LUPRI,'(//A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
         END DO
      END IF
#endif
C
C
#if SIRFCK_DEBUG > 0
      write(lupri,*) 'Leaving SIRFCK'; flush(lupri)
#endif
 9999 CALL QEXIT('SIRFCK')
      RETURN
C
C
 8000 CONTINUE
      WRITE (LUPRI,'(/A/,(8I10))')
     &   ' FATAL ERROR SIRFCK: Unknown IFCTYP; dump of IFCTYP:',
     &   (IFCTYP(I),I=1,NDMAT)
C
      CALL QTRACE(LUPRI)
      CALL QUIT('SIRFCK ERROR: Unknown IFCTYP')
C
      END
C  /* Deck sirfck2 */
      SUBROUTINE SIRFCK2(LU2INT,FN2INT,
     &                   FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,WRK,LWORK)
C
C     Extracted non-Direct code from SIRFCK into SIRFCK2,
C     with file unit and name in parameter list.
C     This means that this can also be used e.g. for
C     the short range Hartree terms in the MC-DFT model.
C                \jkp+hjaaj
C
C output:
C     FMAT   :  NDMAT Fock matrices with 2-el. cont. added
C  input:
C     FMAT   :  NDMAT Fock matrices
C     DMAT   :  NDMAT density matrices
C     NDMAT  :  number of density matrices and Fock matrices
C     ISYMDM :  symmetry of DMAT and FMAT
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
#include "implicit.h"
      CHARACTER*(*) FN2INT
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT), MBAS(8)
      DIMENSION WRK(LWORK)
C
#include "iratdef.h"
#include "dummy.h"
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D4 = 4.0D0, DMP5 = -0.5D0)
C
C Used from common blocks
C   INFORB: NSYM, NBAS(:), N2BASX
C   INFPRI: IPRFCK
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infpri.h"
C
C for HFXFAC:
#include "dftcom.h"


      CALL QENTER('SIRFCK2')

#include "memint.h"
#if SIRFCK_DEBUG > 0
      write(lupri,*)'SIRFCK2: NDMAT=', NDMAT
      write(lupri,*)'SIRFCK2: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
#endif
C
C     Transfer constant factors to DMAT
C     IFCTYP 11,12,13,22; symmetrization or antisymmetrization:
C          factor 2 because factor 0.5 in FCKEN1 and FCKEN2 in order
C          to preserve previous content of FMAT
C     IFCTYP 11,13: factor  2 (see note in FCKD11/FCKD13)
C     IFCTYP 12,22: factor  0.5 from exchange integral
C     IFCTYP 02   : factor -0.5 from exchange integral
      DO I = 1,NDMAT
         IF (HFXFAC .EQ. D0) THEN
C           remove any "do exchange" bits
            J = MOD(IFCTYP(I),10)
            IF (IAND(J,2) .ne. 0) IFCTYP(I) = IFCTYP(I) - 2
         END IF
         IF (IFCTYP(I) .EQ. 11 .OR. IFCTYP(I) .EQ. 13) THEN
            CALL DSCAL(N2BASX,D4,DMAT(1,I),1)
         ELSE IF (IFCTYP(I) .EQ. 2) THEN
            IF (HFXFAC .EQ. D0) THEN
C              This exchange-only matrix is zero
               IFCTYP(I) = 0
            ELSE
               CALL DSCAL(N2BASX,HFXFAC*DMP5,DMAT(1,I),1)
            END IF
         END IF
         IF (MOD(IFCTYP(I),10) .EQ. 0) IFCTYP(I) = 0 ! no contribution, X value does not matter
      END DO
C
C
      IF (LU2INT .GT. 0) CALL GPCLOSE(LU2INT,'KEEP')
      CALL GPOPEN(LU2INT,FN2INT,'OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
      REWIND(LU2INT)
      CALL MOLLAB('BASINFO ',LU2INT,LUPRI)
      READ (LU2INT) MSYM, MBAS, LBUF, NIBUF, NBITS, LENINT4
      IF (IPRFCK .GE. 5) THEN
         WRITE(LUPRI,'(/3A,3I5)')
     &   'SIRFCK2: File "',FN2INT,'" opened. LBUF,NIBUF,NBITS:',
     &   LBUF,NIBUF,NBITS
      END IF
      IF (MSYM .NE. NSYM) THEN
         CALL QUIT('SIRFCK2: NSYM from AOxxxINT .ne. NSYM in common')
      END IF
      IF (NBITS .NE. NIBUF*8) THEN
         CALL QUIT('SIRFCK2: NBITS from AOxxxINT .ne. NIBUF*8')
      END IF
      DO I = 1,MSYM
         IF (NBAS(I) .NE. MBAS(I)) GO TO 110
      END DO
      GO TO 120
  110    WRITE (LUPRI,'(/A/3A,8I5)')
     &      ' SIRFCK2 FATAL ERROR:',
     &      ' NBAS from ',FN2INT,' :',(MBAS(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)')
     &      ' NBAS from /INFORB/ :',(NBAS(I),I=1,NSYM)
         CALL QUIT('SIRFCK2: NBAS(:) from AOxxxINT not consistent')
  120 CONTINUE
C
      IF (LENINT4 .ne. 2*LBUF + NIBUF*LBUF + 1) THEN  ! LENINT4 = length in integer*4 units
         CALL QUIT('Error in LENINT4 from AO 2-el integral file')
      END IF

      KFREE = 1
      LFREE = LWORK
      CALL MEMGET2('INT4','INT',KINT ,LENINT4,WRK,KFREE,LFREE)
      KIINT = KINT + LBUF
      CALL MEMGET2('INTE','IIN4',KIIN4,4*LBUF,WRK,KFREE,LFREE)
C
      CALL MOLLAB('BASTWOEL',LU2INT,LUPRI)
C
      NREC = 0
      NINTS_TOT = 0
    1 CONTINUE
C
         CALL READI4(LU2INT,LENINT4,WRK(KINT))
         NREC = NREC + 1
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,WRK(KIIN4),NINTS)
         IF (NINTS .EQ. 0)  GOTO 1
         IF (NINTS .EQ. -1) GOTO 2
         NINTS_TOT = NINTS_TOT + NINTS
C
         DO 200 I = 1,NDMAT
            IF (IFCTYP(I) .EQ. 13) THEN
C           ... e.g. singlet, symmetric density matrix
             CALL FCKD13(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINTS)
            ELSE IF (IFCTYP(I) .EQ. 22
     &          .OR. IFCTYP(I) .EQ. 23
     &          .OR. IFCTYP(I) .EQ. 12) THEN
C           ... Symmetric or antisymmetric, Only Exchange
C               (type 23: Coulomb term is zero because antisymmetric)
C           ... e.g. singlet, antisymmetric density matrix OR
C           ... or   triplet,     symmetric density matrix OR
C           ... or   triplet, antisymmetric density matrix
             CALL FCKD12(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINTS)
            ELSE IF (IFCTYP(I) .EQ. 11) THEN
C           ... Symmetric, Only Coulomb
             CALL FCKD11(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINTS)
            ELSE IF (IFCTYP(I) .EQ. 3) THEN
C           ... Non-symmetric, Coulomb+Exchange
C           ... e.g. singlet, non-symmetric (general) density matrix
             CALL FCKD03(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINTS)
            ELSE IF (IFCTYP(I) .EQ. 2) THEN
C           ... Non-symmetric, Only Exchange
C           ... e.g. triplet, non-symmetric (general) density matrix
             CALL FCKD02(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINTS)
            ELSE IF (IFCTYP(I) .EQ. 1) THEN
C           ... Non-symmetric, Only Coulomb
             CALL FCKD01(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINTS)
            ELSE IF (MOD(IFCTYP(I),10).EQ. 0 .OR. IFCTYP(I).EQ.21) THEN
C              *0:  nothing
C              21:  antisymmetric, only Coulomb gives a zero matrix
              GOTO 200
            ELSE
              GO TO 8000
            END IF
  200    CONTINUE
C
      GOTO 1
C
    2 CONTINUE
      CALL GPCLOSE(LU2INT,'KEEP')
      IF (IPRFCK .GE. 5) THEN
         WRITE(LUPRI,'(/A,I6,A/I16,A)')
     &   'SIRFCK2: ', NREC,' records with integrals read.',
     &   NINTS_TOT,' integrals processed.'
      END IF
C
C     Finish symmetric and antisymmetric Fock matrices
C     Restore original density matrices
C
      DO 300 I = 1,NDMAT
         IF (IFCTYP(1) .EQ. 11 .OR. IFCTYP(I) .EQ. 12 .OR.
     &       IFCTYP(I) .EQ. 13) THEN
            CALL FCKEN1(FMAT(1,I),NBAST)
         ELSE IF (IFCTYP(I) .EQ. 22 .OR. IFCTYP(I) .EQ. 23) THEN
            CALL FCKEN2(FMAT(1,I),NBAST)
         END IF
C
         IF (IFCTYP(I) .EQ. 11 .OR. IFCTYP(I) .EQ. 13) THEN
            DFAC = D1 / D4
            CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
         ELSE IF (IFCTYP(I) .EQ. 2) THEN
            IF (HFXFAC .NE. D0) THEN
               DFAC = D1 / (DMP5*HFXFAC)
               CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
            END IF
         END IF
C
  300 CONTINUE
C
      CALL QEXIT('SIRFCK2')
      RETURN
C
 8000 CONTINUE
      WRITE (LUPRI,'(/A/,(8I10))')
     &   ' FATAL ERROR SIRFCK2: Unknown IFCTYP; dump of IFCTYP:',
     &   (IFCTYP(I),I=1,NDMAT)
      CALL QTRACE(LUPRI)
      CALL QUIT('SIRFCK2 ERROR: Unknown IFCTYP')
      END
C  /* Deck fckd13 */
      SUBROUTINE FCKD13(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen October 1994: general routine for
C symmetric singlet Fock matrix from symmetric singlet density matrix
C (based on FCKD03)
C
C This subroutine adds (derivative) two-electron integrals to
C Fock matrices. The Fock matrices are assumed
C to be square matrices in full dimension without symmetry reduction
C in size. Remember to zero out the Fock matrices before starting
C to accumulate.
C
C tsaue Feb 1997: indices permuted for exchange
C HJAaJ October 1994: NOTE:
C     Symmetry is implicitly taken into account by the fact that
C     only non-zero integrals are considered.  The symmetry of DMAT
C     is not taken into account.
C DFT modifications T. Helgaker
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DP25 = 0.25D0)
      INTEGER P, Q, R, S
#include "inforb.h"
#include "dftcom.h"
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      EFAC = DP25*HFXFAC
C     ... 0.25 because DMAT is multiplied by extra 2 in SIRFCK
C         in order to get full Coulomb contribution:
C         DINT*(DMAT(R,S) + DMAT(S,R)) = 2*DINT*DMAT(R,S)
C         Compare general code in FCKD03.
      DO INT = 1, LENGTH
C        factor "4" is multiplied on DMAT in SIRFCK
C        DINT = 4*BUF(INT)
         DINT =   BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DINT / 2
         IF (R.EQ.S)              DINT = DINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DINT / 2
         EINT = EFAC*DINT
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(R,S)
         FMAT(P,R) = FMAT(P,R) - EINT*DMAT(Q,S)
         FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(P,S)
         FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(P,R)
         FMAT(P,S) = FMAT(P,S) - EINT*DMAT(Q,R)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(P,Q)
C        ... DMAT(i,j) = DMAT(j,i) has been used above
      END DO
      RETURN
      END
C  /* Deck fckd12 */
      SUBROUTINE FCKD12(FMAT,XMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen and Rika Kobayashi 24-Apr-1992
C     (based on FCKDIR by Henrik Koch and Trygve Helgaker)
C     970304-tsaue: index permuted
C     941011-hjaaj: renamed from LFKDIR to FCKD12
C                   more than halved the work by using
C                   the antisymmetry also in the exchange part
C     DFT modifications T. Helgaker
C
C     This subroutine adds two-electron integrals multiplied
C     by antisymmetric effective density matrices to Fock matrices
C     (e.g. from London reorthonormalization one-index transformations).
C     The Fock matrices are assumed to be square matrices
C     in full dimension without symmetry reduction in size.
C     Remember to zero out the fock matrices before starting
C     to accumulate.
C
C     The matrix is finished by antisymmetrization in FCKES2.
C
C     FMAT(P,R)  = - 0.5 SUM(q,s) (p q | r s) X(s,q)
C                = + 0.5 SUM(q,s) (p q | r s) X(q,s)
C     Note that X(q,s) = -X(s,q)
C
#include "implicit.h"
#include "priunit.h"
      INTEGER P, Q, R, S
#include "inforb.h"
#include "dftcom.h"
      DIMENSION FMAT(NBAST,NBAST), XMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      EFAC = HFXFAC
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = EFAC*BUF(INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
         FMAT(P,R) = FMAT(P,R) - GINT*XMAT(Q,S)
         FMAT(P,S) = FMAT(P,S) - GINT*XMAT(Q,R)
         FMAT(Q,R) = FMAT(Q,R) - GINT*XMAT(P,S)
         FMAT(Q,S) = FMAT(Q,S) - GINT*XMAT(P,R)
CFCKD03 MAERK 941013: bemaerk 0.5 faktor mangler sammenlignet med FCKD03
CFCKD03               and XMAT transposed using XMAT(i,j) = -XMAT(j,i)
CFCKD03  FMAT(P,R) = FMAT(P,R) - EINT*DMAT(S,Q)
CFCKD03  FMAT(P,S) = FMAT(P,S) - EINT*DMAT(R,Q)
CFCKD03  FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(S,P)
CFCKD03  FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd11 */
      SUBROUTINE FCKD11(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen March 2003: add to FMAT for
C IX = 1 (symmetric DMAT), IY = 1 (only Coulomb)
C (based on FCKD13)
C
C The Fock matrices are assumed to be square matrices
C in full dimension without symmetry reduction in size.
C Remember to zero out the fock matrices before starting
C to accumulate.
C
C HJAaJ October 1994: NOTE:
C     Symmetry is implicitly taken into account by the fact that
C     only non-zero integrals are considered.  The symmetry of DMAT
C     is not taken into account.
C
#include "implicit.h"
#include "priunit.h"
      INTEGER P, Q, R, S
#include "inforb.h"
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      DO INT = 1, LENGTH
C        factor "4" is multiplied on DMAT in SIRFCK
C        DINT = 4*BUF(INT)
         DINT =   BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DINT / 2
         IF (R.EQ.S)              DINT = DINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DINT / 2
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(R,S)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(P,Q)
      END DO
      RETURN
      END
C  /* Deck fckd03 */
      SUBROUTINE FCKD03(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C     970303-tsaue : index permutation
C     941011-hjaaj: renamed from FCKDI1 to FCKD03
C     DFT modifications T. Helgaker
C
C     This subroutine adds derivative two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include "implicit.h"
#include "priunit.h"
      INTEGER P, Q, R, S
#include "inforb.h"
#include "dftcom.h"
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)

C
      EFAC = HFXFAC / 2
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
C           coulomb:
         FADD = GINT*(DMAT(R,S)+DMAT(S,R))
         FMAT(P,Q) = FMAT(P,Q) + FADD
         FMAT(Q,P) = FMAT(Q,P) + FADD
         FADD = GINT*(DMAT(P,Q)+DMAT(Q,P))
         FMAT(R,S) = FMAT(R,S) + FADD
         FMAT(S,R) = FMAT(S,R) + FADD
C           exchange:
         GINT = EFAC*GINT
         FMAT(P,R) = FMAT(P,R) - GINT*DMAT(Q,S)
         FMAT(R,P) = FMAT(R,P) - GINT*DMAT(S,Q)
         FMAT(P,S) = FMAT(P,S) - GINT*DMAT(Q,R)
         FMAT(S,P) = FMAT(S,P) - GINT*DMAT(R,Q)
         FMAT(Q,R) = FMAT(Q,R) - GINT*DMAT(P,S)
         FMAT(R,Q) = FMAT(R,Q) - GINT*DMAT(S,P)
         FMAT(Q,S) = FMAT(Q,S) - GINT*DMAT(P,R)
         FMAT(S,Q) = FMAT(S,Q) - GINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd02 */
      SUBROUTINE FCKD02(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen 11-Oct-1994
C     General routine for triplet Fock matrices, based on FCKD03.
C     970304 - tsaue: Index permutation
C
C     This subroutine adds two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include "implicit.h"
#include "priunit.h"
      INTEGER P, Q, R, S
#include "inforb.h"
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)

C
      DO 100 INT = 1, LENGTH
C        factor "-0.5" is multiplied on DMAT in SIRFCK
C        factor HFXFAC is multiplied on DMAT in SIRFCK
C        GINT = - ( HFXFAC * BUF(INT) ) / 2
         GINT = BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
         FMAT(P,R) = FMAT(P,R) + GINT*DMAT(Q,S)
         FMAT(Q,R) = FMAT(Q,R) + GINT*DMAT(P,S)
         FMAT(R,P) = FMAT(R,P) + GINT*DMAT(S,Q)
         FMAT(S,P) = FMAT(S,P) + GINT*DMAT(R,Q)
         FMAT(P,S) = FMAT(P,S) + GINT*DMAT(Q,R)
         FMAT(Q,S) = FMAT(Q,S) + GINT*DMAT(P,R)
         FMAT(R,Q) = FMAT(R,Q) + GINT*DMAT(S,P)
         FMAT(S,Q) = FMAT(S,Q) + GINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd01 */
      SUBROUTINE FCKD01(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     HJAaJ Mar 2003: based on FCKD03
C
C     The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include "implicit.h"
      INTEGER P, Q, R, S
#include "inforb.h"
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
C           coulomb:
         FADD = GINT*(DMAT(R,S)+DMAT(S,R))
         FMAT(P,Q) = FMAT(P,Q) + FADD
         FMAT(Q,P) = FMAT(Q,P) + FADD
         FADD = GINT*(DMAT(P,Q)+DMAT(Q,P))
         FMAT(R,S) = FMAT(R,S) + FADD
         FMAT(S,R) = FMAT(S,R) + FADD
  100 CONTINUE
      RETURN
      END
C  /* Deck fcken1 */
      SUBROUTINE FCKEN1(FMAT,NBAST)
C
C     Hans Joergen Aa. Jensen October 1994.
C
C     This subroutine completes symmetric Fock matrices
C     accumulated by FCKD13 or FCKD02 by symmetrization.
C
C     Note: 0.5 factor is required in order to preserve
C     symmetric content of FMAT from before SIRFCK call.
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION FMAT(NBAST,NBAST)
      INTEGER   P, Q
C
      DO P = 1,NBAST
         DO Q = 1,P
            FPQ = ( FMAT(P,Q) + FMAT(Q,P) ) / 2
            FMAT(P,Q) = FPQ
            FMAT(Q,P) = FPQ
         END DO
      END DO
      RETURN
      END
C  /* Deck fcken2 */
      SUBROUTINE FCKEN2(FMAT,NBAST)
C
C     Hans Joergen Aa. Jensen October 1994.
C
C     This subroutine completes antisymmetric Fock matrices
C     accumulated by FCKD12 by antisymmetrization.
C
C     Note: 0.5 factor is required in order to preserve
C     antisymmetric content of FMAT from before SIRFCK call.
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION FMAT(NBAST,NBAST)
      INTEGER   P, Q
C
      DO P = 1,NBAST
         DO Q = 1,P
            FPQ = ( FMAT(P,Q) - FMAT(Q,P) ) / 2
            FMAT(P,Q) =  FPQ
            FMAT(Q,P) = -FPQ
         END DO
      END DO
      RETURN
      END
C
C  /* Deck cho_fckds1 */
      SUBROUTINE CHO_FCKDS1(FOCK,DEN,WORK,LWORK,ISYDEN,IOPTVP)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Calculate two-electron integral contribution to
C        FOCK matrix using Cholesky decomposition of
C        (1) two-electron integrals, AND
C        if CMO not available on disk,
C        (2) density matrix DEN (decomposed locally).
C
C     Notes:
C        (a) The density matrix *must* be tot. sym. (ISYDEN = 1).
C        (b) FOCK and DEN are assumed stored as full square
C            without symmetry reduction.
C
#include "implicit.h"
      DIMENSION FOCK(*), DEN(*), WORK(LWORK)
#include "infpri.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccisao.h"
#include "ccsdsym.h"
#include "choscf.h"
#include "priunit.h"
#include "ccdeco.h"

      INTEGER NVECD(8), IVECD(8,8)

      PARAMETER (HALF = 0.50D0, ONE = 1.00D0)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_FCKDS1')

      LOGICAL LOCDBG, TINFO
      PARAMETER (LOCDBG = .FALSE., TINFO = .FALSE.)

      LOGICAL IOPTVP
C
C
C     Start timing.
C     -------------

      TIMT = SECOND()

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Check ISYDEN.
C     -------------

      IF (ISYDEN .NE. 1) THEN
         WRITE(LUPRI,'(//,A,A)')
     &   SECNAM,': FATAL ERROR: Density *must* be tot. symmetric !'
         WRITE(LUPRI,'(A,I10,/)')
     &   'Symmetry of DEN on input: ',ISYDEN
         CALL QUIT('SYMMETRY MISMATCH IN '//SECNAM)
      ENDIF

C     Allocation.
C     -----------

      KFOCK = 1
      KDEN  = KFOCK + N2BST(ISYDEN)
      KEND1 = KDEN  + N2BST(ISYDEN)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A)')
     &   'Insufficient memory in ',SECNAM,':',
     &   'Allocation: Symmetry pack'
         WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &   'Need      : ',KEND1-1,
     &   'Available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Symmetry pack DEN and scale.
C     ----------------------------

      TIMPK = SECOND()
      CALL CHO_DENPK(DEN,WORK(KDEN),ISYDEN)
      CALL DSCAL(N2BST(ISYDEN),HALF,WORK(KDEN),1)
      TIMPK = SECOND() - TIMPK

C     Cholesky decompose symmetry blocks of DEN if necessary;
C     else read occupied part of CMO from disk.
C     Set up index array for the CMO 'vectors'.
C     -------------------------------------------------------

      TIMDC = SECOND()

!     IF (NEWSCF) CCMODSK = .FALSE.   ! hjaaj: NEWSCF is never defined in this version of Dalton /120327-hjaaj
      IF (.NOT. CCMODSK) THEN

         CALL DCOPY(N2BST(ISYDEN),WORK(KDEN),1,WORK(KFOCK),1)

         KVECD  = KEND1
         ICOUNT = 0

         DO ISYM = 1,NSYM

            IVECD(ISYM,ISYM) = ICOUNT

            KEND2 = KVECD + IVECD(ISYM,ISYM)
            LWRK2 = LWORK - KEND2 + 1
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A,I1)')
     &         'Insufficient memory in ',SECNAM,':',
     &         'Allocation: Density decomposition, sym. block: ',ISYM
               WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &         'Need      : ',KEND2-1,
     &         'Available : ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            KOFFD = KFOCK + IAODIS(ISYM,ISYM)
            CALL SIR_CHODECO(WORK(KOFFD),NBAS(ISYM),WORK(KEND2),LWRK2,
     &                       THRDC,NVECD(ISYM))

            ICOUNT = ICOUNT + NBAS(ISYM)*NVECD(ISYM)

         ENDDO

         KEND2 = KVECD + ICOUNT
         LWRK2 = LWORK - KEND2

      ELSE

         ICOUNT = 0
         DO ISYM = 1,NSYM
            IVECD(ISYM,ISYM) = ICOUNT
            NVECD(ISYM) = NDVCS(ISYM)
            ICOUNT = ICOUNT + NBAS(ISYM)*NVECD(ISYM)
         ENDDO
         LEN  = ICOUNT

         KVECD = KEND1
         KEND2 = KVECD + LEN
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A)')
     &      'Insufficient memory in ',SECNAM,':',
     &      'Allocation: occ. CMO read'
            WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &      'Need      : ',KEND2-1,
     &      'Available : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL WOPEN2(LUCCMO,FCCMO,64,0)
         IADR = 1
         CALL GETWA2(LUCCMO,FCCMO,WORK(KVECD),IADR,LEN)
         CALL WCLOSE2(LUCCMO,FCCMO,'KEEP')

      ENDIF

      TIMDC = SECOND() - TIMDC

C     Debug: Print number of vectors.
C     Debug: Check decomposition/CMO read.
C     ------------------------------------

      IF (LOCDBG) THEN
         IF (.NOT. CCMODSK) THEN
            WRITE(LUPRI,'(/,2X,A,A,/)')
     &      SECNAM,': Cholesky decomposition of density:'
         ELSE
            WRITE(LUPRI,'(/,2X,A,A,/)')
     &      SECNAM,': occ. CMO read:'
         ENDIF
         WRITE(LUPRI,'(2X,A,/,2X,A)')
     &   'Sym.  No. vectors  Offset',
     &   '-------------------------'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(3X,I1,4X,I6,7X,I6)')
     &      ISYM,NVECD(ISYM),IVECD(ISYM,ISYM)
         ENDDO
         WRITE(LUPRI,'(2X,A,/)')
     &   '-------------------------'
      ENDIF

C     Calculate contributions.
C     ------------------------

      TIMCL = SECOND()
      CALL DZERO(WORK(KFOCK),N2BST(ISYDEN))
      CALL SCF_CHOFCK1(WORK(KDEN),WORK(KVECD),WORK(KFOCK),
     &                 WORK(KEND2),LWRK2,NVECD,IVECD,TINFO,IOPTVP)
      TIMCL = SECOND() - TIMCL

C     Add result to FOCK matrix.
C     --------------------------

      TIMAD = SECOND()
      DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYDEN)
         DO IB = 1,NBAS(ISYMB)
            B = IBAS(ISYMB) + IB
            KOFF1 = KFOCK + IAODIS(ISYMA,ISYMB)
     &            + NBAS(ISYMA)*(IB - 1)
            KOFF2 = NBAST*(B - 1) + IBAS(ISYMA) + 1
            CALL DAXPY(NBAS(ISYMA),ONE,WORK(KOFF1),1,FOCK(KOFF2),1)
         ENDDO
      ENDDO
      TIMAD = SECOND() - TIMAD

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Print timing.
C     -------------

      IF (TINFO .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sym. packing density : ',TIMPK,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for decomposing density  : ',TIMDC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating contr.   : ',TIMCL,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for adding to Fock matr. : ',TIMAD,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck scf_chofck1 */
      SUBROUTINE SCF_CHOFCK1(DEN,CMO,FOCK,WORK,LWORK,NVECD,IVECD,
     &                       TINFO,IOPTVP)
C
C     Written by Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Calculate two-electron integral contributions to the AO Fock
C        matrix based on Cholesky decomposed integrals.
C        Result is added into Fock matrix which is assumed allocated
C        as symmetry packed (!) full square.
C        The density matrix is used for the Coulomb part
C        and the MO coefficient matrix (or the Cholesky vector
C        representing the density matrix) for the exchange part.
C
C     Version 1:
C        The density and the MO coefficients are assumed
C        TOTAL SYMMETRIC, as appropriate for an SCF ground state
C        closed-shell optimization.
C
C     NB: NO SPARSITY IS USED WHATSOEVER !!!!!
C
#include "implicit.h"
      DIMENSION DEN(*), CMO(*), FOCK(*), WORK(LWORK)
      INTEGER IVECD(8,8), NVECD(8)
      LOGICAL TINFO
#include "iratdef.h"
#include "infpri.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccdeco.h"
#include "sirftim.h"
#include "priunit.h"
Casm
#include "dftcom.h"
#include "choscf.h"
Casm

      INTEGER IHOLE(8,8), NACHH(8), IALPHK(8,8), NALPHK(8)

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (ISTORE = 2)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'SCF_CHOFCK1')

      LOGICAL IOPTVP


C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMO = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIMX = ZERO

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(DEN(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': MO Coefficient Matrix on Entry')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IVECD(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            NUMK = NVECD(ISYM)
            CALL OUTPUT(CMO(KOFF),1,NUMB,1,NUMK,NUMB,NUMK,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(FOCK(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Set up local index arrays.
C     --------------------------

      DO ISYCHO = 1,NSYM

         ICOUNT = 0

         DO ISYMK = 1,NSYM

            ISYMA = MULD2H(ISYMK,ISYCHO)

            IALPHK(ISYMA,ISYMK) = ICOUNT

            ICOUNT = ICOUNT + NBAS(ISYMA)*NVECD(ISYMK)

         ENDDO

         NALPHK(ISYCHO) = ICOUNT

      ENDDO

C     Start of Cholesky loop.
C     -----------------------

      DO ISYCHO = 1,NSYM

         IF (N2BST(ISYCHO)  .LE. 0)  GO TO 999
         IF (NUMCHO(ISYCHO) .LE. 0)  GO TO 999
         IF (NALPHK(ISYCHO) .LE. 0)  GO TO 999

C        Allocation: Cholesky.
C        ---------------------

         LREAD = MEMRD(1,ISYCHO,ISTORE)
         LSCR1 = MAX(LREAD,NALPHK(ISYCHO))

         KCHOL = KEND0
         KSCR1 = KCHOL + N2BST(ISYCHO)
         KEND1 = KSCR1 + LSCR1
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - allocation: Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,A,I10,/)')
     &      'Need     : ',KEND1,
     &      'Available: ',LWORK
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NALPHK(ISYCHO)
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: Batch'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Need per vector    : ',MINMEM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available for batch: ',LWRK1
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

         KHOLE = KEND1
         KEND2 = KHOLE + NVEC*NALPHK(ISYCHO)
         LWRK2 = LWORK - KEND2 + 1

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

         IF (LOCDBG) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'BATCH INFO FROM ',SECNAM,':'
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total number of vectors          : ',NUMCHO(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'AO vector dimension              : ',N2BST(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total memory                     : ',LWORK
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Memory used for AO vectors       : ',KEND1-1
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Memory left for (alpha,k) vectors: ',LWRK1
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Min. memory for (alpha,k) vectors: ',MINMEM
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of batches                : ',NBATCH
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of vectors per batch      : ',NVEC
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Memory remainder                 : ',LWRK2
         ENDIF

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

C           Set up index arrays for HOLE.
C           -----------------------------

            ICOUNH = 0
            DO ISYMK = 1,NSYM

               ISYMA = MULD2H(ISYMK,ISYCHO)

               IHOLE(ISYMA,ISYMK) = ICOUNH

               NACHH(ISYMA) = NBAS(ISYMA)*NUMV

               ICOUNH = ICOUNH + NACHH(ISYMA)*NVECD(ISYMK)

            ENDDO

            JCHOL1 = NVEC*(IBATCH - 1) + 1

            DO JVEC = 1,NUMV

               JCHOL = JCHOL1 + JVEC - 1

C              Read vector.
C              ------------

               DTIME = SECOND()
               CALL CHO_READN(WORK(KCHOL),JCHOL,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,ISTORE,WORK(KSCR1),LSCR1)
c              xnorm = ddot(nnbst(isycho),work(kchol),1,work(kchol),1)
c              write(lupri,'(a,i4,a,i2,a,d12.4)') 'Vector',jchol,
c    &               ' of symmetry',isycho,'. Norm: ',xnorm 
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

C              Coulomb part.
C              -------------

               IF ((ISYCHO .EQ. 1) .AND. IOPTVP) THEN
                  DTIME = SECOND()
                  CONST = TWO*DDOT(N2BST(ISYCHO),DEN,1,WORK(KCHOL),1)
                  CALL DAXPY(N2BST(ISYCHO),CONST,WORK(KCHOL),1,
     &                       FOCK,1)
                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME
               ENDIF

Casm
C              No exchange for LDA. Thus, no transformation
C
               IF (ITSLDA .OR. HFXFAC .EQ. 0.0D0) GOTO 788
Casm

C              Hole transformation.
C              --------------------

               DTIME = SECOND()
               DO ISYMD = 1,NSYM

                  ISYMA = MULD2H(ISYMD,ISYCHO)
                  ISYMK = ISYMD

                  NA    = NBAS(ISYMA)
                  NTOTA = MAX(NA,1)
                  ND    = NBAS(ISYMD)
                  NTOTD = MAX(ND,1)
                  NK    = NVECD(ISYMK)

                  KOFFC = KCHOL + IAODIS(ISYMD,ISYMA)
                  KOFFH = IVECD(ISYMD,ISYMK) + 1
                  KOFFS = KSCR1 + IALPHK(ISYMA,ISYMK)

                  CALL DGEMM('T','N',NA,NK,ND,ONE,
     &                       WORK(KOFFC),NTOTD,CMO(KOFFH),NTOTD,
     &                       ZERO,WORK(KOFFS),NTOTA)

               ENDDO
               DTIME = SECOND() - DTIME
               TIMO  = TIMO     + DTIME

C              Copy into HOLE array.
C              ---------------------

               DTIME = SECOND()
               DO ISYMK = 1,NSYM

                  ISYMA = MULD2H(ISYMK,ISYCHO)

                  DO K = 1,NVECD(ISYMK)

                     KOFFS = KSCR1 + IALPHK(ISYMA,ISYMK)
     &                     + NBAS(ISYMA)*(K - 1)
                     KOFFH = KHOLE + IHOLE(ISYMA,ISYMK)
     &                     + NACHH(ISYMA)*(K - 1)
     &                     + NBAS(ISYMA)*(JVEC - 1)

                     CALL DCOPY(NBAS(ISYMA),WORK(KOFFS),1,
     &                          WORK(KOFFH),1)

                  ENDDO

               ENDDO

  788          CONTINUE

               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME

            ENDDO
Casm
C           Do not include exchange for DFT/LDA
C
            IF (ITSLDA .OR. HFXFAC .EQ. 0.0D0) GOTO 789
Casm

C           Exchange part.
C           --------------

            DTIME = SECOND()
            DO ISYMK = 1,NSYM

               ISYMA = MULD2H(ISYMK,ISYCHO)
               ISYMB = ISYMA

               NA    = NBAS(ISYMA)
               NTOTA = MAX(NA,1)
               NB    = NBAS(ISYMB)
               NTOTB = MAX(NB,1)
               NVK   = NUMV*NVECD(ISYMK)

               KOFFF = IAODIS(ISYMA,ISYMB) + 1
               KOFFH = KHOLE + IHOLE(ISYMA,ISYMK)
Casm
C              Very dirty hack to solve weird bug when linking
C              against some ATLAS libraries.
C              See the end of cc_cholesky.F for more info
Casm
               CALL OTHER_DGEMM('N','T',NA,NB,NVK,
     &                    HFXFAC*XMONE,WORK(KOFFH),NTOTA,WORK(KOFFH),
     &                    NTOTB,ONE,FOCK(KOFFF),NTOTA)

            ENDDO
            DTIME = SECOND() - DTIME
            TIMX  = TIMX     + DTIME
 
  789       CONTINUE

         ENDDO

  999    CONTINUE

      ENDDO

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(DEN(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': MO Coefficient Matrix on Exit')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IVECD(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            NUMK = NVECD(ISYM)
            CALL OUTPUT(CMO(KOFF),1,NUMB,1,NUMK,NUMB,NUMK,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(FOCK(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Print timing.
C     -------------

      CSIRF = CSIRF + TIMC
      XSIRF = XSIRF + TIMX
      RSIRF = RSIRF + TIMR
      OSIRF = OSIRF + TIMO
      SSIRF = SSIRF + TIMS

      IF (TINFO .OR. LOCDBG) THEN
c      IF (.true.) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O and squaring     : ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Coulomb  part        : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for partial MO transform.: ',TIMO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sorting (copy)       : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for exchange part        : ',TIMX,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cho_denpk */
      SUBROUTINE CHO_DENPK(DEN,DENPK,ISYDEN)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Symmetry pack DEN.
C
#include "implicit.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccsdsym.h"

      DIMENSION DEN(*), DENPK(*)

      DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYDEN)
         DO IB = 1,NBAS(ISYMB)
            B = IBAS(ISYMB) + IB
            DO IA = 1,NBAS(ISYMA)
               A = IBAS(ISYMA) + IA
               KABPK = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(IB - 1) + IA
               KABFL = NBAST*(B - 1) + A
               DENPK(KABPK) = DEN(KABFL)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck sir_chodeco */
      SUBROUTINE SIR_CHODECO(XMAT,NDIM,WORK,LWORK,THRS,NVECS)
C
C     Henrik Koch   10-may-2002
C
#include "implicit.h"
#include "infpri.h"
#include "priunit.h"
C
      DIMENSION XMAT(NDIM,NDIM),WORK(NDIM,*)
C
      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'SIR_CHODECO')
      PARAMETER (TOL = 5.00D-13)
C
c     DO I = 1,NDIM
c        WRITE(LUPRI,*) 'Diagonal :',I,XMAT(I,I)
c     ENDDO
C
      IF (LWORK .LT. NDIM*NDIM) CALL QUIT('Not enough memory')
C
      DO I = 1,NDIM
C
         XMAX = XMAT(1,1)
         IMAX = 1
         DO J = 2,NDIM
            IF (ABS(XMAT(J,J)) .GT. XMAX) THEN
               XMAX = XMAT(J,J)
               IMAX = J
            ENDIF
         ENDDO
C 
         IF (XMAX .LT. 0.0D0) THEN
            WRITE(LUPRI,'(A,A,A)')
     &      '*** ',SECNAM,': NOTICE ***'
            WRITE(LUPRI,'(A,1P,D30.15)')
     &      '    Negative diagonal encountered: ',XMAX
            IF (ABS(XMAX) .GT. TOL) THEN
               WRITE(LUPRI,'(A,/)')
     &         '    - unable to continue!'
               CALL QUIT('ERROR IN '//SECNAM)
            ELSE
               WRITE(LUPRI,'(A)')
     &         '    - program continues nevertheless....'
            ENDIF
         ENDIF
C
         IF (ABS(XMAX) .LT. THRS) THEN
            NVECS = I-1
            RETURN
         ENDIF
C
c        WRITE(LUPRI,*)  'Largest diagonal  : ',NDIM,I,IMAX,XMAX
C
         DO J = 1,NDIM
            WORK(J,I) = XMAT(J,IMAX)/SQRT(XMAX)
         ENDDO
C
         DO J = 1,NDIM
            DO K = 1,NDIM
               XMAT(K,J) = XMAT(K,J) - WORK(K,I)*WORK(J,I)
            ENDDO
         ENDDO
C
         DO J = 1,NDIM
            XMAT(J,IMAX) = 0.0D0
            XMAT(IMAX,J) = 0.0D0
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck cho_fckds3 */
      SUBROUTINE CHO_FCKDS3(FOCK,DEN,WORK,LWORK,ISYDEN,IOPTVP)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Calculate two-electron integral contribution to
C        FOCK matrix using Cholesky decomposition of
C        the two-electron integrals.
C
C     Notes:
C        (a) FOCK and DEN are assumed stored as full square
C            without symmetry reduction.
C        (b) This is a pilot version: density driven but
C            with use of neither sparsity of the Cholesky vectors
C            nor screening.
C
#include "implicit.h"
      DIMENSION FOCK(*), DEN(*), WORK(LWORK)
#include "infpri.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccisao.h"
#include "ccsdsym.h"
#include "priunit.h"

      PARAMETER (HALF = 0.50D0, ONE = 1.00D0)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_FCKDS3')

      LOGICAL LOCDBG, TINFO, IOPTVP
      PARAMETER (LOCDBG = .FALSE., TINFO = .FALSE.)

C     Start timing.
C     -------------

      TIMT = SECOND()

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Allocation.
C     -----------

      KFOCK = 1
      KDEN  = KFOCK + N2BST(ISYDEN)
      KEND1 = KDEN  + N2BST(ISYDEN)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A)')
     &   'Insufficient memory in ',SECNAM,':',
     &   'Allocation: Symmetry pack'
         WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &   'Need      : ',KEND1-1,
     &   'Available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Symmetry pack DEN and scale.
C     ----------------------------

      TIMPK = SECOND()
      CALL CHO_DENPK(DEN,WORK(KDEN),ISYDEN)
      CALL DSCAL(N2BST(ISYDEN),HALF,WORK(KDEN),1)
      TIMPK = SECOND() - TIMPK

C     Calculate contributions.
C     ------------------------

      TIMCL = SECOND()
      CALL DZERO(WORK(KFOCK),N2BST(ISYDEN))
      CALL SCF_CHOFCK3(WORK(KDEN),WORK(KFOCK),WORK(KEND1),LWRK1,
     &                 ISYDEN,TINFO,IOPTVP)
      TIMCL = SECOND() - TIMCL

C     Add result to FOCK matrix.
C     --------------------------

      TIMAD = SECOND()
      DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYDEN)
         DO IB = 1,NBAS(ISYMB)
            B = IBAS(ISYMB) + IB
            KOFF1 = KFOCK + IAODIS(ISYMA,ISYMB)
     &            + NBAS(ISYMA)*(IB - 1)
            KOFF2 = NBAST*(B - 1) + IBAS(ISYMA) + 1
            CALL DAXPY(NBAS(ISYMA),ONE,WORK(KOFF1),1,FOCK(KOFF2),1)
         ENDDO
      ENDDO
      TIMAD = SECOND() - TIMAD

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Print timing.
C     -------------

      IF (TINFO .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sym. packing density : ',TIMPK,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating contr.   : ',TIMCL,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for adding to Fock matrix: ',TIMAD,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck scf_chofck3 */
      SUBROUTINE SCF_CHOFCK3(DEN,FOCK,WORK,LWORK,ISYDEN,TINFO,IOPTVP)
C
C     Written by Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Calculate two-electron integral contributions to the AO Fock
C        matrix based on Cholesky decomposed integrals.
C        Result is added into Fock matrix which is assumed allocated
C        as symmetry packed (!) full square.
C        The density matrix is used for both Coulomb and exchange part.
C
C     Formula used:
C        F(ab) = F(ab) +   2    * Sum(gd) D(dg) (ab|dg)
C                      - HFXFAC * Sum(gd) D(dg) (ad|bg)
C
C     NB: NO SPARSITY IS USED WHATSOEVER !!!!!
C
#include "implicit.h"
      DIMENSION DEN(*), FOCK(*), WORK(LWORK)
      LOGICAL TINFO
#include "iratdef.h"
#include "infpri.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccdeco.h"
#include "priunit.h"
#include "dftcom.h"

      INTEGER IOFFC(8)

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'SCF_CHOFCK3')

      LOGICAL IOPTVP

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIMX = ZERO

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(DEN(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(FOCK(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Start Cholesky loop.
C     --------------------

      DO ISYCHO = 1,NSYM

         IF (N2BST(ISYCHO)  .LE. 0) GO TO 999
         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 999

         ISYMX = MULD2H(ISYCHO,ISYDEN)

C        Allocation: Cholesky.
C        ---------------------

         LREAD = 2*NNBST(ISYCHO) + (NNBST(ISYCHO) - 1)/2 + 2

         KCSCR = KEND0
         KREAD = KCSCR + N2BST(ISYCHO)
         KEND1 = KREAD + LREAD
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - allocation: Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,A,I10,/)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MAXAG = -1
         DO ISYMG = 1,NSYM
            ISYMA = MULD2H(ISYMG,ISYMX)
            MAXAG = MAX(MAXAG,NBAS(ISYMA)*NBAS(ISYMG))
         ENDDO

         MINMEM = N2BST(ISYCHO) + MAXAG
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: Batch'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Need per vector    : ',MINMEM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available for batch: ',LWRK1
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

         IF (LOCDBG) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'BATCH INFO FROM ',SECNAM,':'
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total number of vectors          : ',NUMCHO(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'AO Cholesky vector dimension     : ',N2BST(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'AO intermediate vector dimension : ',MAXAG
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total memory                     : ',LWORK
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Min. memory needed for batch     : ',MINMEM
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of batches                : ',NBATCH
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of vectors per batch      : ',NVEC
         ENDIF

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JCHOL1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            KCHOL = KEND1
            KXMAT = KCHOL + N2BST(ISYCHO)*NUMV

C           Set up index array IOFFC.
C           -------------------------

            ICOUNT = 0
            DO ISYMG = 1,NSYM
               ISYMB = MULD2H(ISYMG,ISYCHO)
               IOFFC(ISYMG) = ICOUNT
               ICOUNT = ICOUNT + NBAS(ISYMB)*NUMV*NBAS(ISYMG)
            ENDDO

            DO JVEC = 1,NUMV

               JCHOL = JCHOL1 + JVEC - 1

C              Read vector.
C              ------------

               DTIME = SECOND()
               IOPT  = 2
               CALL CHO_READN(WORK(KCSCR),JCHOL,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPT,WORK(KREAD),LREAD)
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

C              Coulomb part.
C              -------------

               IF (IOPTVP .AND. (ISYCHO .EQ. ISYDEN)) THEN
                  DTIME = SECOND()
                  CONST = TWO*DDOT(N2BST(ISYDEN),DEN,1,WORK(KCSCR),1)
                  CALL DAXPY(N2BST(ISYDEN),CONST,WORK(KCSCR),1,
     &                       FOCK,1)
                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME
               ENDIF

C              Prepare for exchange part.
C              --------------------------

             IF (HFXFAC .ne. 0.0d0) THEN

C              Reorder: L(bg,J) -> L(b,J,g).
C              -----------------------------

               DTIME = SECOND()
               DO ISYMG = 1,NSYM

                  ISYMB = MULD2H(ISYMG,ISYCHO)

                  IF (NBAS(ISYMB) .GT. 0) THEN

                     DO G = 1,NBAS(ISYMG)

                        KOFF1 = KCSCR + IAODIS(ISYMB,ISYMG)
     &                        + NBAS(ISYMB)*(G - 1)
                        KOFF2 = KCHOL + IOFFC(ISYMG)
     &                        + NBAS(ISYMB)*NUMV*(G - 1)
     &                        + NBAS(ISYMB)*(JVEC - 1)

                        CALL DCOPY(NBAS(ISYMB),WORK(KOFF1),1,
     &                                         WORK(KOFF2),1)

                     ENDDO

                  ENDIF

               ENDDO
               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME

             END IF   ! (HFXFAC .ne. 0.0d0)

            ENDDO   ! JVEC = 1, NUMV

           IF (HFXFAC .NE. 0.0D0) THEN

            DTIME = SECOND()
            DO ISYMG = 1,NSYM

               ISYMB = MULD2H(ISYMG,ISYCHO)
               ISYMD = MULD2H(ISYMG,ISYDEN)
               ISYMA = MULD2H(ISYMB,ISYDEN)

C              Calculate intermediate:
C              X(aJ,g) = Sum(d) L(aJ,d) * D(g,d).
C              ----------------------------------

               NTOAJ = MAX(NBAS(ISYMA)*NUMV,1)
               NTOTD = MAX(NBAS(ISYMD),1)

               KOFF1 = KCHOL + IOFFC(ISYMD)
               KOFF2 = IAODIS(ISYMD,ISYMG) + 1

               CALL DGEMM('N','N',
     &                    NBAS(ISYMA)*NUMV,NBAS(ISYMG),NBAS(ISYMD),
     &                    ONE,WORK(KOFF1),NTOAJ,DEN(KOFF2),NTOTD,
     &                    ZERO,WORK(KXMAT),NTOAJ)

C              Calculate exchange contribution:
C              F(a,b) = F(a,b) - Sum(g) X(a,Jg) * L(b,Jg).
C              -------------------------------------------

               NTOTA = MAX(NBAS(ISYMA),1)
               NTOTB = MAX(NBAS(ISYMB),1)

               KOFF3 = KCHOL + IOFFC(ISYMG)
               KOFF4 = IAODIS(ISYMA,ISYMB) + 1

               CALL DGEMM('N','T',
     &                    NBAS(ISYMA),NBAS(ISYMB),NUMV*NBAS(ISYMG),
     &                    HFXFAC*XMONE,WORK(KXMAT),NTOTA,WORK(KOFF3),
     &                    NTOTB,ONE,FOCK(KOFF4),NTOTA)

            ENDDO
            DTIME = SECOND() - DTIME
            TIMX  = TIMX     + DTIME
           END IF   ! (HFXFAC .ne. 0.0d0)

         ENDDO


  999    CONTINUE

      ENDDO

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(DEN(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(FOCK(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Print timing.
C     -------------

      IF (TINFO .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O and squaring     : ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Coulomb  part        : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sorting (copy)       : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for exchange part        : ',TIMX,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
! -- end of sirius/sirfck.F --
