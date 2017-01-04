C
C  /* Deck so_excit1 */
      SUBROUTINE SO_EXCIT1(TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                     BSRLON,SECMAT,EXENG,
CClark:11/01/2016
     &                     BETHE,STOPP,
CClark:end
     &                     WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     SPAS :4/12-1997 : call to ccsd_init1 included. This should allow
C                       to run the RPA part without going through the
C                       CC program.
C     SPAS & PFP: 5/11-2013: triplet excitation energies included
C     RF: Jan. 2016   : Rewritten to loop over models instead of
C                       repeating code for each.
C
C     PURPOSE: Main driver routine for the calculation of excitation
C              energies and transition moments with the atomic integral
C              direct SOPPA program.
C
      use so_info, only: sop_num_models,
     &                   sop_models, sop_mod_fullname,
     &                   so_num_active_models,
     &                   so_get_active_models,
     &                   so_has_doubles,
     &                   so_singles_second,
     &                   so_needs_densai,
     &                   sop_dens_label,
     &                   sop_model_rpa, sop_model_rpad,
     &                   sop_conv_thresh
#ifdef VAR_MPI
      use so_parutils, only: soppa_initialize_slaves,
     &                       soppa_release_slaves
#endif
C
#include "implicit.h"
C
#include "priunit.h"
C
#include "cbiexc.h"
CPi for NT2AMX, NT1AM(), NT2AM() (t2 amplitude vectors?)
#include "ccsdsym.h"
#include "infdim.h"
#include "inforb.h"
C for irat
#include "iratdef.h"
#include "maxaqn.h"
C for mxshel, mxprim
#include "maxorb.h"
C for mxcont
#include "aovec.h"
#include "mxcent.h"
#include "nuclei.h"
#include "soppinf.h"
#include "symmet.h"
#include "wrkrsp.h"

C
      PARAMETER (HALF = 0.5D0,ESUDIP = 64604.885D0,ESUECD = 471.44360D0)
C
CSPAS:5/11-13: HRPA included
CRF  Defer last dimension size...
CRF  actual size will be set to NUM_ACTIVE later
      DIMENSION TRLEN(3,NSYM,MXNEXI,*),
     &          TRVEL(3,NSYM,MXNEXI,*)
      DIMENSION TQLEN(3,3,NSYM,MXNEXI,*),
     &          TQVEL(3,3,NSYM,MXNEXI,*)
      DIMENSION TTLEN(10,NSYM,MXNEXI,*)
      DIMENSION TRLON(3,NSYM,MXNEXI,*),
     &          TRMAG(3,NSYM,MXNEXI,*)
      DIMENSION BSRLON(3,NSYM,MXNEXI,*),
     &          EXENG(NSYM,MXNEXI,*)
CClark:11/01/2016
      DIMENSION BETHE(3,LQ,*),STOPP(3,LVEL,2,*)
CClark:end
CKeinSPASmehr
      DIMENSION SECMAT(3,MXNEXI,NSYM)
      DIMENSION WORK(LWORK)
C
      CHARACTER*5  MODEL
      CHARACTER*11 FULLNAME
C

      integer :: num_active
      logical :: active_models(sop_num_models)
      logical :: doubles, need_t2, rpa_only
#ifdef VAR_MPI
!
! This variable ensures that common blocks are only sent to the slaves
! once.
      LOGICAL update_common_blocks, get_mxcall
      LOGICAL so_get_herdir, herdir, so_get_direct

C      herdir = so_get_herdir()
      update_common_blocks = .true.
      get_mxcall = .false.
#endif
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_EXCIT1')
C
C
C     Get the active models and their number
      CALL SO_GET_ACTIVE_MODELS (ACTIVE_MODELS)
      NUM_ACTIVE = COUNT ( ACTIVE_MODELS )
      RPA_ONLY = ACTIVE_MODELS(SOP_MODEL_RPA) .AND.
     &           (NUM_ACTIVE.EQ.1)
C
C
C     Set the correct convergence threshold -- THREXC from
C     cbiexc.h
      sop_conv_thresh = THREXC
C
C-------------------------------------------
C     Start timing the AO-SOPPA calculation.
C-------------------------------------------
C
      DTIME   = SECOND()
      TIMTOT  = DTIME
      CALL GETTIM (DUMMY,WTIME)
      TIMWTO  = WTIME
C
CPFP 2009-Jun
C
C---------------------------------------------------------
C     Initialize TRIPLET keyword used inside SOPPA module.
C---------------------------------------------------------
C
      IF (EXCTRP) THEN
        TRIPLET = .TRUE.
      ELSE
        TRIPLET = .FALSE.
      END IF
C
Cend-PFP
C
C
C-------------------------------------------
C     Initializes memory statistics routine.
C-------------------------------------------
C
      CALL SO_MEMMAX ('START',LWORK)
C
C---------------------------------------------------------------
C     If only RPA initialize some pointers and sort AO-integrals
C     into distributions.
C---------------------------------------------------------------
C
      IF ( RPA_ONLY ) THEN
         CALL RP_INIT(WORK,LWORK)
#ifdef VAR_MPI
C     In pure aorpa calculations, eridi has not been initialized
         herdir = so_get_herdir()
         get_mxcall = so_get_direct()
#endif
      ENDIF

#ifdef VAR_MPI
C
C In AORPA we cannot be sure that the number of integral distributions
C is known in advance. This is needed in order to allocate memory for
C load-distribution later on
      if (get_mxcall) then
         if (herdir) then
            call herdi1(work(1), lwork, iprint)
         else
            kccfbt  = 1
            kindxbt = kccfbt + mxprim*mxcont
            kend3   = kindxbt + 8*mxshel*mxcont/irat
            lwork3 = lwork - kend3
            call eridi1(kodcl1, kodcl2, kodbc1, kodbc2,
     &                  krdbc1, krdbc2, kodpp1, kodpp2,
     &                  krdpp1, krdpp2, kfree,  lfree,
     &                  kend3, work(kccfbt), work(kindxbt),
     &                  work(kend3), lwork3, iprint)
         endif
         get_mxcall = .false.
      endif
#endif
C
C----------------------------------
C     Set print level for AO-SOPPA.
C----------------------------------
C
      IPRSOP = IPREXC
C
C---------------------------------------------------
C     Initialize a few pointerarrays and dimensions.
C---------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_INIT
      DTIME     = SECOND()  - DTIME
      SOTIME(2) = SOTIME(2) + DTIME
C
C---------------------------------
C     1. allocation of work space.
C---------------------------------
C
      LFOCKD  = NORBT
C
      KFOCKD  = 1
      KEND1   = KFOCKD  + LFOCKD
      LWORK1  = LWORK   - KEND1
C
CRF if anything but AORPA
      IF ( .NOT.RPA_ONLY ) THEN
C
CSPAS:5/11-2013: Patricio had the following statement before the IF statement
         LT2AM   = NT2AMX
CPFP
C        For triplet calculations LT2AM has to be twice as large,
C        because there are two version of the correlation coefficients
C
         IF (TRIPLET) LT2AM = NT2AMX*2
Cend-PFP
         LDENSIJ = NIJDEN(1)
         LDENSAB = NABDEN(1)
         LDENSAI = NAIDEN(1)
C
         KT2AM   = KFOCKD  + LFOCKD
         KDENSIJ = KT2AM   + LT2AM
         KDENSAB = KDENSIJ + LDENSIJ
         KDENSAI = KDENSAB + LDENSAB
         KEND1   = KDENSAI + LDENSAI
         LWORK1  = LWORK   - KEND1
      ELSE
         LT2AM   = 0
         LDENSAI = 0
         LDENSAB = 0
         LDENSIJ = 0
C        Force a crash of this is used
         KT2AM   = LWORK + 1
         KDENSIJ = LWORK + 1
         KDENSAB = LWORK + 1
         KDENSAI = LWORK + 1
C
      END IF
C
      CALL SO_MEMMAX ('SO_EXCIT1.1',LWORK1)
      IF (LWORK1 .LT .0) CALL STOPIT('SO_EXCIT1.1',' ',KEND1,LWORK)
C
C------------------------------------------------
C     Get MO-energies (the fock matrix diagonal).
C------------------------------------------------
C
      DTIME     = SECOND()
      CALL SO_MOENERGY(WORK(KFOCKD),WORK(KEND1),LWORK1)
      DTIME     = SECOND()  - DTIME
      SOTIME(3) = SOTIME(3) + DTIME
C
C------------------------------------------------------
C     Construct property-integrals and write to LUPROP.
C------------------------------------------------------
C
      CALL SO_PRPINT('EXCITA',NLBTOT,WORK(KEND1),LWORK1)
C
C---------------------------------------------------------------------
C     Initialize arrays of transition moments and excitation energies.
C---------------------------------------------------------------------
C
CSPAS:5/11-13: HRPA included
      CALL DZERO(TRLEN, 3*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(TRVEL, 3*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(TQLEN, 3*3*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(TQVEL, 3*3*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(TTLEN, 10*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(TRLON, 3*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(TRMAG, 3*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(BSRLON,3*NSYM*MXNEXI*NUM_ACTIVE)
      CALL DZERO(SECMAT,3*NSYM*MXNEXI)
      CALL DZERO(EXENG, NSYM*MXNEXI*NUM_ACTIVE)
CKeinSPASmehr
CClark:09/03/2016
      CALL DZERO(BETHE,3*LQ*NUM_ACTIVE)
      CALL DZERO(STOPP,3*LVEL*2*NUM_ACTIVE)
CClark:end
C
C========================================================
C     Determine excitation energies, response vectors and
C     linear transition moments.
C========================================================
C
C----------------------------------------------------------
C     Adjust number of trialvectors for each excitation and
C     write information to output.
C----------------------------------------------------------
C
      IF ( NSAVMX .LT. 2 ) THEN
C
         NSAVMX = 2
C
         WRITE(LUPRI,'(1X,A,/,A,I2,A)')
     &   'NOTICE: Maximum number of trial vectors for each'//
     &   ' excitation',' is raised to',NSAVMX, ' as this is'//
     &   ' allowed.'
C
      END IF
C
      WRITE(LUPRI,'(/,1X,A,I2,/)')
     &'Maximum number of trial vectors for each'//
     &' excitation is ',NSAVMX
C
C-------------------------------------------------------------
C     Loop over symmetry of the excitation vector.
C     NOTE: The excited state has the symmtry of the reference
C     state times the symmetry of the excitation vector.
C-------------------------------------------------------------
C
      DO 100 ISYM = 1, NSYM
C
         CALL RSPSET
C
         IF (NEXCIT(ISYM) .GT. 0) THEN
C
C---------------------------------------
C           2. allocation of work space.
C---------------------------------------
C
            KEXVAL = KEND1
            KEND2  = KEXVAL  + NEXCIT(ISYM)
            LWORK2 = LWORK   - KEND2
C
            CALL SO_MEMMAX ('SO_EXCIT1.2',LWORK2)
            IF (LWORK2 .LT.0) CALL STOPIT('SO_EXCIT1.2',' ',KEND2,LWORK)
C
            ISYMTR = ISYM
            LEXVAL = NEXCIT(ISYM)
C
C=========================================
C           Loop over the possible methods
C=========================================
C
            IOUT = 0
            DO IMODEL = 1, SOP_NUM_MODELS
               IF (.NOT. ACTIVE_MODELS(IMODEL) ) CYCLE

C              Increment output pointer
               IOUT = IOUT + 1
C----------------------------------------
C              Look up info on this model
C----------------------------------------
               MODEL = SOP_MODELS(IMODEL)
               FULLNAME = SOP_MOD_FULLNAME(IMODEL)
               DOUBLES = SO_HAS_DOUBLES(IMODEL)
               NEED_T2 = DOUBLES.OR.SO_SINGLES_SECOND(MODEL)
C
               IF (IPRSOP .GE. 1) THEN
                  WRITE(LUPRI,9000)
                  WRITE(LUPRI,'(1X,A)') TRIM(FULLNAME)//':'
                  WRITE(LUPRI,'(A,I8,/,A,I8)')
     &            ' Perturbation symmetry     :',ISYM,
     &            ' p-h variables             :',2*NT1AM(ISYM)
                  IF (DOUBLES) THEN
CPi 01-04-16
C                     NVARPT  = 2*(NT1AM(ISYM)+NT2AM(ISYM))
                     NVARPT  = 2*(NT1AM(ISYM)+N2P2HOP(ISYM))
                     WRITE(LUPRI,'(A,I8,/,A,I8)')
C     &               ' 2p-2h variables           :',2*NT2AM(ISYM),
     &               ' 2p-2h variables           :',2*N2P2HOP(ISYM),
Cend-Pi
     &               ' Total number of variables :',NVARPT
                  ENDIF
                  WRITE(LUPRI,9001)
               END IF
C
C              Adjust number of excitations of this symmetry down to
C              the highest possible number
C
               IF (DOUBLES .AND.(.NOT.IMODEL.EQ.SOP_MODEL_RPAD) ) THEN
C                  NEXCI = MIN(NEXCIT(ISYM),NT1AM(ISYM)+NT2AM(ISYM))
                  NEXCI = MIN(NEXCIT(ISYM),NT1AM(ISYM)+N2P2HOP(ISYM))
               ELSE
                  NEXCI = MIN(NEXCIT(ISYM),NT1AM(ISYM))
               ENDIF
C
               IF ( NEXCI .EQ. 0 ) THEN
                  WRITE (LUPRI,'(1X,A)')
     &               'There are no possible excitations of this '//
     &               'symmetry!'
                  CYCLE
               ELSEIF ( NEXCI .LT. NEXCIT(ISYM) ) THEN
                  WRITE(LUPRI,'(1X,A,I5,/,1X,A)')
     &            'NOTICE: The number of excitations are reduced to '
     &            ,NEXCI,'as this is the full dimension of the '//
     &            'excitation space'
               END IF
C
C-----------------------------------------------------
C              Get T2 amplitudes and density matrices.
C-----------------------------------------------------
C
               IF (NEED_T2) THEN
                  CALL SO_GETT2(SOP_DENS_LABEL(IMODEL),
     &                          WORK(KT2AM),LT2AM,WORK(KDENSAI),
     &                          LDENSAI,WORK(KDENSIJ),LDENSIJ,
     &                          WORK(KDENSAB),LDENSAB,
     &                          WORK(KEND2),LWORK2)
               ENDIF
C
C---------------------------------------------------------
C              Ready the slaves for parallel calculations.
C---------------------------------------------------------
C
#ifdef VAR_MPI
               call soppa_initialize_slaves ( update_common_blocks,
     &                                       WORK(KT2AM), LT2AM, model)
               update_common_blocks = .false.
#endif
C
C-----------------------------------------
C              Initialize DENSAI if needed
C-----------------------------------------
C
               IF ( SO_NEEDS_DENSAI(MODEL))
     &            CALL DZERO(WORK(KDENSAI),LDENSAI)
C
C------------------------------------------------------------------
C              Determine excitation energies and excitation
C              vectors. The excitation vectors are written to file.
C------------------------------------------------------------------
C
               IF ( IMODEL .NE. SOP_MODEL_RPAD) THEN
                  CALL SO_RSPLEX(MODEL,ISYMTR,NEXCI,WORK(KEXVAL),LEXVAL,
     &                           WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                           LDENSAB,WORK(KDENSAI),LDENSAI,
     &                           WORK(KT2AM),
     &                           LT2AM,WORK(KFOCKD),LFOCKD,
     &                           WORK(KEND2),LWORK2)
               ELSE ! RPA(D) still has its own driver...
                  CALL DC_RSPLEX(ISYMTR,NEXCI,WORK(KEXVAL),LEXVAL,
     &                           WORK(KDENSIJ),LDENSIJ,WORK(KDENSAB),
     &                           LDENSAB,WORK(KDENSAI),LDENSAI,
     &                           WORK(KT2AM),
     &                           LT2AM,WORK(KFOCKD),LFOCKD,
     &                           WORK(KEND2),LWORK2)
               ENDIF
C
C-----------------------------------------------------
C              Copy exciation energies to EXENG.
C-----------------------------------------------------
C
               DO IEXVAL = 1,NEXCI
                  IF (IEXVAL .LE. MXNEXI)
     &               EXENG(ISYM,IEXVAL,IOUT) = WORK(KEXVAL-1+IEXVAL)
                  IF (IEXVAL .GT. MXNEXI)
     &               WRITE(LUPRI,*)
     &                    'WARNING: IEXVAL greater than MXNEXI.'
               END DO
C
C-------------------------------------------------
C              Determine transition moments.
C-------------------------------------------------
C
               CALL SO_TRMOM(MODEL,ISYMTR,NLBTOT,NEXCI,
     &                       WORK(KT2AM),LT2AM,WORK(KDENSIJ),LDENSIJ,
     &                       WORK(KDENSAB),LDENSAB,
     &                       WORK(KDENSAI),LDENSAI,
     &                       TRLEN(1,1,1,IOUT),TRVEL(1,1,1,IOUT),
     &                       TQLEN(1,1,1,1,IOUT),TQVEL(1,1,1,1,IOUT),
     &                       TTLEN(1,1,1,IOUT),
     &                       TRLON(1,1,1,IOUT),TRMAG(1,1,1,IOUT),
     &                       BSRLON(1,1,1,IOUT),EXENG(1,1,IOUT),
CClark:11/01/2016
     &                       BETHE(1,1,IOUT),STOPP(1,1,1,IOUT),
CClark:end
     &                       WORK(KEND2),LWORK2)

C
C------------------------------------------------------------------
C              Determine charge radius of excited states (RPA only)
C------------------------------------------------------------------
C
               IF ( IMODEL .EQ. SOP_MODEL_RPA)
     &            CALL RP_CHARGE(ISYMTR,NEXCI,SECMAT,WORK(KEND2),LWORK2)
#ifdef VAR_MPI
C-------------------------------------------------------
C              Release slaves to the global node-driver.
C-------------------------------------------------------
               call soppa_release_slaves()
#endif

            END DO
         END IF
C
  100 CONTINUE
C
C---------------------------------
C     3. allocation of work space.
C---------------------------------
C
      LPARRA = LSOTIM
C
      KPARRA = KEND1
      KEND3  = KPARRA + LPARRA
      LWORK3 = LWORK  - KEND3
C
      CALL SO_MEMMAX ('SO_EXCIT1.3     ',LWORK3)
      IF (LWORK3 .LT.0) CALL STOPIT('SO_EXCIT1.3',' ',KEND3,LWORK)
C
C---------------------------------------------------
C     Print memory statistics for SOPPA subroutines.
C---------------------------------------------------
C
      CALL SO_MEMMAX('STATISTICS      ',0)
C
C-----------------------------------------
C     Print timings for SOPPA subroutines.
C-----------------------------------------
C
      CALL SO_TIME(TIMTOT,TIMWTO,WORK(KPARRA),LPARRA)
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_EXCIT1')
C
      RETURN
C
 9000 FORMAT(/' -----------------------------------')
 9001 FORMAT(' -----------------------------------')
      END

#ifdef VAR_MPI
      LOGICAL FUNCTION SO_GET_HERDIR()
C        WORK-AROUND TO THE ISSUE THAT TOO MANY COMMON BLOCKS HAS
C        THE VARIABLE "SKIP"
C        RETURNS THE VALUE OF HERDIR from ccsdinp.h

#include "ccsdinp.h"
      SO_GET_HERDIR = HERDIR
      RETURN
      END
C
      LOGICAL FUNCTION SO_GET_DIRECT()
C        WORK-AROUND TO THE ISSUE THAT TOO MANY COMMON BLOCKS HAS
C        THE VARIABLE "SKIP"
C        RETURNS THE VALUE OF DIRECT from ccsdinp.h

#include "ccsdinp.h"
      SO_GET_DIRECT = DIRECT
      RETURN
      END
#endif

