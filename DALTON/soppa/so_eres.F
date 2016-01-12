C  /* Deck so_eres */
      SUBROUTINE SO_ERES(MODEL,  NOLDTR, NNEWTR,  DENSIJ,  LDENSIJ, 
     &                   DENSAB, LDENSAB, T2MP,    LT2MP,   FOCKD,
     &                   LFOCKD, DENSAI,  LDENSAI, NIT,     ISYMTR,
#ifdef VAR_MPI
     &                   AssignedIndices, maxnumjobs,     
#endif
     &                   WORK,   LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Frederik Beyer & Stephan P. A. Sauer: 27.08.2013: call to ERI corrected
C
C     PURPOSE: Driver routine for making a linear transformation of
C              a trialvector with the SOPPA hessian matricx E[2]. 
C              The trial vector consists of four parts TR1E, TR1D, 
C              TR2E, and TR2D. E refers to excitations and D to 
C              de-excitations. 1 refer to the one-particle part and
C              2 to the two-particle part. The linear transformed 
C              trialvector is refered to as the resultvector and is
C              kept in four corresponding arrays. For the linear
C              transformation with E[2] the result vector is in RES1E,
C              RES1D, RES2E, and RES2D.
C              The linear transformation is driven over atomic orbitals,
C              and E[2] is not constructed explicitly.
C
      use so_info, only: so_singles_first, so_singles_second, 
     &                   so_has_doubles

#ifdef VAR_MPI
      use so_parutils, only: soppa_comm_active, soppa_num_active,
     &                       soppa_nint 
#endif

#include "implicit.h"

#ifdef VAR_MPI
#include "mpif.h"
#endif

#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), T2MP(LT2MP)
      DIMENSION FOCKD(LFOCKD)   
      DIMENSION DENSAI(LDENSAI) !intent(out) 
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*5
      integer :: inewtr, nnewtr, isymd1, ntosym
      integer :: thisinewtr
#ifdef VAR_MPI
      integer :: maxnumjobs, nloopidx
      integer :: AssignedIndices(maxnumjobs)
#endif
C
#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "soppinf.h"


      logical :: singles_first, singles_second, doubles  
#ifdef VAR_MPI
#include "iprtyp.h"
#include "infpar.h"
C#include "parsoppa.h"
C      integer, save :: numprocs
       integer :: numprocs
       logical :: loadbal_dyn
       double precision ::  timeini, timefin
       numprocs = nodtot + 1
C
C  When to do dynamic load_balancing..?
C  The sooner the better, but is the first iteration close to the
C  average?
       loadbal_dyn = .false.
       if (nit .eq. 1) loadbal_dyn = .true.
#endif

C
C------------------
C     Add to trace.
C------------------
C
#ifdef VAR_MPI
      if (mynum .eq. 0 ) then
#endif      
      CALL QENTER('SO_ERES')
#ifdef VAR_MPI
C     Slaves need to zero DENSAI if it is (re)calculated.      
      elseif ( (MODEL .EQ.'AOSOP'.OR.MODEL .EQ. 'DCRPA')
     &                   .AND. NIT .EQ. 1) THEN
         CALL DZERO(DENSAI,LDENSAI)
      endif
#endif
C
C-------------------------------------------------------------
C     Determine which terms to incude
C-------------------------------------------------------------      
C
C      singles_first = model .ne. 'DCRPA'
      singles_first = so_singles_first(model)
C      singles_second = model .ne. 'AORPA'
      singles_second = so_singles_second(model)
C      doubles = (model .eq. 'DCRPA') .or. (model .eq. 'AOSOP') .or.
C     &          (model .eq. 'AOSOC')
      doubles = so_has_doubles(model)
C
C------------------------------------------------------------------
C     Determine the symmetry of the result vector from the symmetry
C     of the trial vector ISYMTR, and the opperator symmtry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      LCMO   = NLAMDT
C
      KCMO    = 1
      KEND1   = KCMO  + LCMO
      LWORK1  = LWORK - KEND1
C
      CALL SO_MEMMAX ('SO_ERES.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_ERES.1',' ',KEND1,LWORK)
C
C-------------------------------------------------------
C     Get the matrix which contains the MO coefficients.
C-------------------------------------------------------
C
#ifdef VAR_MPI
C Only master reads...
      IF (MYNUM .EQ. 0) THEN
#endif
         DTIME      = SECOND()
         CALL SO_GETMO(WORK(KCMO),LCMO,WORK(KEND1),LWORK1)
         DTIME      = SECOND()   - DTIME
         SOTIME(1)  = SOTIME(1) + DTIME
#ifdef VAR_MPI
      ENDIF
C Should probably use non-blocking collectives, where implemented
!      IF (MPI_VERSION.GE.3) THEN
!         MPI_IBCAST(WORK(KCMO),LCMO, my_MPI_INTEGER, 0,
!     &              MPI_COMM_WORLD, IERR)
!      ELSE
      CALL MPI_BCAST(WORK(KCMO), LCMO, MPI_REAL8, 0,
     &               SOPPA_COMM_ACTIVE, IERR)

!      ENDIF
#endif
C
C---------------------------------
C     Work space allocation no. 2.
C---------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRES1E  = NT1AM(ISYMTR)
      LRES1D  = NT1AM(ISYMTR)
      LFOCK   = N2BST(ISYRES)
      LDENS   = N2BST(ISYMTR)
      LBTR1E  = NT1AO(ISYMTR)
      LBTR1D  = NT1AO(ISYMTR)
      LBTJ1E  = NMATAV(ISYMTR)
      LBTJ1D  = NMATAV(ISYMTR)
                 
      IF(DOUBLES) THEN
         LTR2E   = NT2AM(ISYMTR)
         LTR2D   = NT2AM(ISYMTR)
         LRES2E  = NT2AM(ISYMTR)
         LRES2D  = NT2AM(ISYMTR)
      ELSE
         LTR2E   = 0
         LTR2D   = 0
         LRES2E  = 0
         LRES2D  = 0
      ENDIF 

      IF (singles_second) THEN
         LSIGAI1 = NT1AO(ISYRES)
         LSIGAI2 = NT1AO(ISYRES)
         LSIGDA1 = NMATAV(ISYRES)
         LSIGDA2 = NMATAV(ISYRES)
         LAIJ    = NRHFT*NRHFT
         LAAB    = NVIRT*NVIRT
      ELSE
         LSIGAI1 = 0
         LSIGAI2 = 0
         LSIGDA1 = 0
         LSIGDA2 = 0
         LAIJ    = 0
         LAAB    = 0
      ENDIF
C
      KTR1E   = KEND1
      KTR1D   = KTR1E   + LTR1E
      KTR2E   = KTR1D   + LTR1D
      KTR2D   = KTR2E   + LTR2E

      KRES1E  = KTR2D   + LTR2D
      KRES1D  = KRES1E  + LRES1E
      KRES2E  = KRES1D  + LRES1D
      KRES2D  = KRES2E  + LRES2E

      KFOCK   = KRES2D  + LRES2D
      KDENS   = KFOCK   + LFOCK
      KBTR1E  = KDENS   + LDENS
      KBTR1D  = KBTR1E  + LBTR1E
      KBTJ1E  = KBTR1D  + LBTR1D
      KBTJ1D  = KBTJ1E  + LBTJ1E

      KSIGAI1 = KBTJ1D  + LBTJ1D
      KSIGAI2 = KSIGAI1 + LSIGAI1
      KSIGDA1 = KSIGAI2 + LSIGAI2
      KSIGDA2 = KSIGDA1 + LSIGDA1

      KAIJ    = KSIGDA2 + LSIGDA2
      KAAB    = KAIJ    + LAIJ
      KEND2   = KAAB    + LAAB
#ifdef VAR_MPI
C     MPI -- Allocate timings array
      if ( loadbal_dyn ) then
         KTIMING = KEND2
         KEND2   = KTIMING + SOPPA_NINT
         CALL DZERO(WORK(KTIMING), SOPPA_NINT)
      endif 
#endif
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('SO_ERES.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_ERES.2',' ',KEND2,LWORK)
C
C----------------------------
C     Initialize AIJ and AAB.
C----------------------------
C
      CALL DZERO(WORK(KAIJ),LAIJ)
      CALL DZERO(WORK(KAAB),LAAB)
#ifdef VAR_MPI
C MPI -- ONLY MASTER DOES THE READING
      IF ( MYNUM .EQ. 0 ) THEN
#endif
C
C----------------------------------------------
C     Open files with trial and result vectors.
C----------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
C
      IF (DOUBLES) THEN
         CALL SO_OPEN(LUTR2E,FNTR2E,LTR2E)
         CALL SO_OPEN(LUTR2D,FNTR2D,LTR2D)
         CALL SO_OPEN(LURS2E,FNRS2E,LRES2E)
         CALL SO_OPEN(LURS2D,FNRS2D,LRES2D)
      ENDIF
#ifdef VAR_MPI
      ENDIF
#endif
C

      IF ( IPRSOP. GE. 7 
#ifdef VAR_MPI
     &      .AND. (MYNUM .EQ.0) 
#endif
     &                     ) THEN ! Only printing related.
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
         DO 50 INEWTR = 1,NNEWTR
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
            INEW = NOLDTR + INEWTR
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
            CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,INEW)
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
               CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,INEW)
            ENDIF
C
            WRITE(LUPRI,'(/,I3,A)') INEWTR,'. new trial vector'
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
            IF (DOUBLES) THEN
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &            (I,WORK(KTR2E+I-1),WORK(KTR2D+I-1),I=1,LTR2E)
            ENDIF
   50    CONTINUE
      END IF
C
C================================================
C     Loop over number of excitations considered.
C================================================
C
      DO 100 INEWTR = 1,NNEWTR
C
C-------------------------------------------------
C        Determine pointer to INEWTR trial vector.
C-------------------------------------------------
C
         INEW = NOLDTR + INEWTR
C
C--------------------------------------------------------------
C        Initialize RES1E, RES1D, RES2E, RES2D, SIGAI1, SIGAI2,
C                   SIGDA1, SIGDA2 and FOCK
C--------------------------------------------------------------
C
         CALL DZERO(WORK(KRES1E),LRES1E)
         CALL DZERO(WORK(KRES1D),LRES1D)
        
         IF (DOUBLES) CALL DZERO(WORK(KRES2E),LRES2E)
         IF (DOUBLES) CALL DZERO(WORK(KRES2D),LRES2D)
C
         IF (SINGLES_SECOND) THEN
            CALL DZERO(WORK(KSIGAI1),LSIGAI1)
            CALL DZERO(WORK(KSIGAI2),LSIGAI2)
            CALL DZERO(WORK(KSIGDA1),LSIGDA1)
            CALL DZERO(WORK(KSIGDA2),LSIGDA2)
         ENDIF
C
         CALL DZERO(WORK(KFOCK),LFOCK)
C
C--------------------------
C        Read trial vector.
C--------------------------
C
#ifdef VAR_MPI
         IF ( MYNUM .EQ. 0 ) THEN
#endif
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
            CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KTR2E),LTR2E,LUTR2E,FNTR2E,INEW)
               CALL SO_READ(WORK(KTR2D),LTR2D,LUTR2D,FNTR2D,INEW)
            ENDIF
#ifdef VAR_MPI
         ENDIF
C Rememember that the trial-vectors are contigous in memory         
         LTRTOT = LTR1E + LTR1D + LTR2E + LTR2D
         CALL MPI_BCAST(WORK(KTR1E), LTRTOT, MPI_REAL8, 
     &                  0, SOPPA_COMM_ACTIVE, IERR)
C Note for future adjustment: 
C It may be worthwile using non-blocking communication here
C since the next two calls only use the singles vectors.
C
C Also, depending on the communication overhead, it may be better
C to transform on host only and then send it to the slaves..?
#endif         
C
C---------------------------------------------------
C        Calculate RPA-density matrices in AO basis.
C---------------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_AODENS(WORK(KDENS),LDENS,WORK(KCMO),LCMO,
     &                  WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,ISYMTR,
     &                  WORK(KEND2),LWORK2)
         DTIME     = SECOND()  - DTIME
         SOTIME(6) = SOTIME(6) + DTIME
C
C--------------------------------------------
C        Backtransformation of trial vectors.
C--------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_BCKTR(WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                 WORK(KBTR1E),LBTR1E,WORK(KBTR1D),LBTR1D,
     &                 WORK(KBTJ1E),LBTJ1E,WORK(KBTJ1D),LBTJ1D,
     &                 WORK(KCMO),LCMO,ISYMTR)
         DTIME     = SECOND()  - DTIME
         SOTIME(7) = SOTIME(7) + DTIME

C
C=======================================================
C        Start the loop over distributions of integrals.
C=======================================================
C
         IF (DIRECT) THEN
            NTOSYM = 1
            DTIME  = SECOND()
            IF (HERDIR) THEN
               CALL HERDI1(WORK(KEND2),LWORK2,IPRINT)
               KINDXB = KEND2
            ELSEIF(INEWTR.EQ.1) THEN
C Can we get away with doing this once?
C Because we CAN'T move KEND2 forwards in each loop cycle,
C since that'll be leaking memory               
               KCCFB1 = KEND2
               KINDXB = KCCFB1 + MXPRIM*MXCONT
               KEND2  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
               LWORK2 = LWORK  - KEND2

               CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                     KODPP1,KODPP2,KRDPP1,KRDPP2,KFREE,LFREE,
     &                     KEND2,WORK(KCCFB1),WORK(KINDXB),WORK(KEND2),
     &                     LWORK2,IPRINT)

               KEND2  = KFREE
               LWORK2 = LFREE
            ENDIF
            DTIME     = SECOND()  - DTIME
            SOTIME(8) = SOTIME(8) + DTIME
         ELSE
            NTOSYM = NSYM
         ENDIF
C
         ICDEL1  = 0
C
#ifdef VAR_MPI
C-----------------------------------------------------------------
C        In MPI calculations, we need to distribute the indices.
C        On first pass, we set this using the pre-sorted approach.
C-----------------------------------------------------------------
C         IF ( (nit .ne. 1 .or. inewtr .ne. 1) 
         IF ( .not. (nit .eq. 1 .and. inewtr .eq. 1) 
     &                        .or.  numprocs .eq. 1 ) THEN
C           After first pass, load-balancing has been done
            CONTINUE
         ELSEIF( mynum .eq. 0 )THEN
C           Master does the balancing
            call presortloadbal_parsoppa(AssignedIndices,   maxnumjobs,
     &                     work(kindxb), work(kend2), lwork2)
         ELSE 
C           Slaves receives the scatter, send arguments should be ignored
            call mpi_scatter( AssignedIndices, maxnumjobs, mpi_integer,
     &                        AssignedIndices, maxnumjobs, mpi_integer,
     &                        0, soppa_comm_active, ierr )
         ENDIF
#endif
         DO 210 ISYMD1 = 1,NTOSYM
C
            IF (DIRECT) THEN
               IF (HERDIR) THEN
                  NTOT = MAXSHL
               ELSE
                  NTOT = MXCALL
               ENDIF
            ELSE
               NTOT = NBAS(ISYMD1)
            ENDIF
#ifdef VAR_MPI
            IF(numprocs .gt. 1) THEN
               NLOOPIDX = maxnumjobs
            ELSE
               NLOOPIDX = NTOT
               loadbal_dyn = .false.
            ENDIF
#endif
C

C-------------------------------------------------
C           Main loop over integral-distributions.
C-------------------------------------------------   

#ifdef VAR_MPI
C------------------------------------------------------------------
C           For_parallel calculations, we have som stuff to set up.
C------------------------------------------------------------------
            DO 220 ILLLDUMMY = 1, nloopidx
               if ( numprocs .gt. 1 ) then 
                  ILLL = assignedIndices(illldummy)
C                 A zero indicates that we have no more work, exit loop
                  IF (ILLL .eq. 0) exit
                  if ( loadbal_dyn ) timeini = mpi_wtime()
               else
                  ILLL = ILLLDUMMY
               endif
#else
            DO 220 ILLL = 1,NTOT
#endif
C               print *, ILLL
C------------------------------------------------
C              If direct calculate the integrals.
C------------------------------------------------
               IF (DIRECT) THEN
C
                  DTIME  = SECOND()
                  IF (HERDIR) THEN
                    CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
     &                          IPRINT)
                  ELSE
C
                     CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                           WORK(KODCL1),WORK(KODCL2),
     &                           WORK(KODBC1),WORK(KODBC2),
     &                           WORK(KRDBC1),WORK(KRDBC2),
     &                           WORK(KODPP1),WORK(KODPP2),
     &                           WORK(KRDPP1),WORK(KRDPP2),
     &                           WORK(KCCFB1),WORK(KINDXB), 
     &                           WORK(KEND2),LWORK2,IPRINT)
C
                  ENDIF
                  DTIME     = SECOND()  - DTIME
                  SOTIME(9) = SOTIME(9) + DTIME
C
                  LRECNR  = ( (NBUFX(0) -1) / IRAT ) + 1
                  KRECNR  = KEND2
                  KEND2B   = KRECNR + LRECNR
                  LWORK2B  = LWORK  - KEND2B
C
                  CALL SO_MEMMAX ('SO_ERES.2B',LWORK2B)
                  IF (LWORK2 .LT. 0) 
     &                CALL STOPIT('SO_ERES.2B',' ',KEND2B,LWORK)
C
               ELSE
                  NUMDIS = 1
                  KEND2B = KEND2
               ENDIF
C
C-------------------------------------------------------------------
C   Loop over number of distributions in disk.
C   In the case of ERI there are more than one distribution and IDEL2
C   loops over them and the actual index of the delta orbital IDEL is 
C   then obtained from the array INDEXA. In the case of a not direct 
C   calculation there is only one distribution on the disk, which 
C   implies that IDEL2 is always 1 and that IDEL is systematically 
C   incremented by one each time.
C--------------------------------------------------------------------
C
               DO 230 IDEL2 = 1,NUMDIS
C

                  IF (DIRECT) THEN
                     IDEL  = INDEXA(IDEL2)
                     ISYMD = ISAO(IDEL) !keeps track of current symmetry
                  ELSE
                     IDEL  = IBAS(ISYMD1) + ILLL
                     ISYMD = ISYMD1
                  ENDIF
C
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
C
                  IT2DEL(IDEL) = ICDEL1
                  ICDEL1       = ICDEL1 + NT2BCD(ISYDIS)
C
C---------------------------------------------
C                 Work space allocation no. 3.
C---------------------------------------------
C
                  LXINT  = NDISAO(ISYDIS)
C
                  KXINT   = KEND2B
                  KEND3   = KXINT + LXINT
                  LWORK3  = LWORK - KEND3
C
                  CALL SO_MEMMAX ('SO_ERES.3',LWORK3)
                  IF (LWORK3 .LT. 0) 
     &                CALL STOPIT('SO_ERES.3',' ',KEND3,LWORK)
C
C--------------------------------------------
C                 Read in batch of integrals.
C--------------------------------------------
C
                  DTIME      = SECOND()
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                        WORK(KRECNR),DIRECT)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(10) = SOTIME(10) + DTIME
C                  
C
C---------------------------------------------
C                    Work space allocation no. 4.
C---------------------------------------------
C
                  ISAIJ = MULD2H(ISYMD,1)
C
                  IF ((NVIR(ISYMD) .GT. 0) .AND.
     &                (DOUBLES.OR.SINGLES_SECOND) ) THEN
                     LT2M1 = NT2BCD(ISAIJ)
                  ELSE
                     LT2M1 = 0
                  END IF
C
                  KT2M1   = KEND3
                  KEND4   = KT2M1  + LT2M1
                  LWORK4  = LWORK  - KEND4
C
                  IF (DOUBLES.OR.SINGLES_SECOND) THEN
                     CALL SO_MEMMAX ('SO_ERES.4',LWORK4)
                     IF (LWORK4 .LT. 0) 
     &                   CALL STOPIT('SO_ERES.4',' ',KEND4,LWORK)
C                   
C
C------------------------------------------------------------
C                    Construct the partially back-transformed T2
C                    MP-amplitudes.
C------------------------------------------------------------
C
                     DTIME      = SECOND()
                     CALL SO_T2M1(WORK(KT2M1),LT2M1,T2MP,LT2MP,
     &                         WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                         WORK(KEND4),LWORK4)
                     DTIME      = SECOND()   - DTIME
                     SOTIME(12) = SOTIME(12) + DTIME
                  ENDIF
C
C---------------------------------------------
C                 Work space allocation no. 5.
C---------------------------------------------
C
                  LDSRHF = NDSRHF(ISYMD)
C
                  KDSRHF  = KEND4
                  KEND5   = KDSRHF + LDSRHF
                  LWORK5  = LWORK  - KEND5
C
                  CALL SO_MEMMAX ('SO_ERES.5',LWORK5)
                  IF (LWORK5 .LT. 0) 
     &                CALL STOPIT('SO_ERES.5',' ',KEND5,LWORK)
C
C----------------------------------------------------------------
C                 Transform one index in the integral batch to an
C                 occupied index.
C----------------------------------------------------------------
C
                  DTIME      = SECOND()
                  ISYMLP = 1
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     &                        ISYMLP,WORK(KEND5),LWORK5,ISYDIS)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(13) = SOTIME(13) + DTIME
C
C-------------------------------------------------------------------
C                 Calculate part of the second order density matrix.
C-------------------------------------------------------------------
C
                  DTIME      = SECOND()
                  IF ( (MODEL.EQ.'AOSOP'.OR.MODEL.EQ.'DCRPA') .AND. 
     &                           (NIT.EQ.1) .AND. (  INEWTR.EQ.1) ) THEN
                     CALL SO_DENSAI1(DENSAI,LDENSAI,WORK(KDSRHF),LDSRHF,
     &                               WORK(KCMO),LCMO,WORK(KT2M1),LT2M1,
     &                               ISYMD,ISYDIS,WORK(KEND5),
     &                               LWORK5)

                  END IF
                  DTIME      = SECOND()   - DTIME
                  SOTIME(41) = SOTIME(41) + DTIME
C
C----------------------------------------------
C                 Calculate the AO-Fock matrix.
C----------------------------------------------
C
                  DTIME      = SECOND()
                  IF (TRIPLET) THEN
C
                     CALL CC_AOFOCK3(WORK(KXINT),WORK(KDENS),
     &                               WORK(KFOCK),WORK(KEND5),
     &                               LWORK5,IDEL,ISYMD,ISYMTR)
                  ELSE
                     CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),WORK(KFOCK),
     &                              WORK(KEND5),LWORK5,
     &                              IDEL,ISYMD,.FALSE.,
     &                              'crashifTrue',ISYMTR)
                  ENDIF
                  DTIME      = SECOND()   - DTIME
                  SOTIME(11) = SOTIME(11) + DTIME
C
                  IF (singles_second) THEN
C
C----------------------------------------------------------------------
C                 Calculate part of the result vectors RES1E and RES1D,
C                 specifically the first and the second term in eqs.
C                 (34,35). Also calculate Aij and Aab in eqs. (43,44).
C----------------------------------------------------------------------
C


                     DTIME      = SECOND()
                     CALL SO_RES_A(WORK(KRES1E),LRES1E,
     &                             WORK(KRES1D),LRES1D,
     &                             WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                             WORK(KDSRHF),LDSRHF,WORK(KCMO),LCMO,
     &                             WORK(KT2M1),LT2M1,WORK(KAIJ),LAIJ,
     &                             WORK(KAAB),LAAB,INEWTR,ISYMD,ISYDIS,
     &                             ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                     DTIME      = SECOND()   - DTIME
                     SOTIME(14) = SOTIME(14) + DTIME
                  ENDIF
C
                  IF (DOUBLES) THEN
C
C-------------------------------------------------------------------
C                 Calculate the part of the result vectors RES1E and
C                 RES1D which originate from the C matrices. See 
C                 eqs. (72) and (73).
C-------------------------------------------------------------------
C
                     DTIME      = SECOND()
                     CALL SO_RES_TCB(WORK(KRES1E),LRES1E,WORK(KRES1D),
     &                               LRES1D,WORK(KTR2E),LTR2E,
     &                               WORK(KTR2D),LTR2D,
     &                               WORK(KDSRHF),LDSRHF,
     &                               WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                               ISYMTR,WORK(KEND5),LWORK5)
                     DTIME      = SECOND()   - DTIME
                     SOTIME(29) = SOTIME(29) + DTIME

C
C----------------------------------------------------------------------
C                 Construct C-contribution to 2p2h result vectors RES2E
C                 and RES2D.
C----------------------------------------------------------------------
C
                     DTIME      = SECOND()
                     CALL SO_RES_CB(WORK(KRES2E),LRES2E,WORK(KRES2D),
     &                             LRES2D,
     &                             WORK(KDSRHF),LDSRHF,
     &                             WORK(KBTR1E),LBTR1E,
     &                             WORK(KBTR1D),LBTR1D,
     &                             WORK(KBTJ1E),LBTJ1E,
     &                             WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,
     &                             IDEL,ISYMD,ISYDIS,ISYMTR,WORK(KEND5),
     &                             LWORK5)
                     DTIME      = SECOND()   - DTIME
                     SOTIME(15) = SOTIME(15) + DTIME
C                  
                  ENDIF
C
                  IF (singles_second) THEN
C
C--------------------------------------------------------------------
C                 Construct SIGMAI1(ALFA,I) and SIGMAI2(ALFA,I) which
C                 are used in SO_RES_B.
C--------------------------------------------------------------------
C
                     DTIME      = SECOND()
                     CALL SO_SIGAI(WORK(KSIGAI1),LSIGAI1,WORK(KSIGAI2),
     &                             LSIGAI2,
     &                             WORK(KT2M1),LT2M1,WORK(KXINT),LXINT,
     &                             WORK(KBTR1E),LBTR1E,
     &                             WORK(KBTR1D),LBTR1D,
     &                             WORK(KBTJ1E),LBTJ1E,
     &                             WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,
     &                             ISYMD,ISYDIS,
     &                             ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                     DTIME      = SECOND()   - DTIME
                     SOTIME(16) = SOTIME(16) + DTIME
C
C--------------------------------------------------------------------
C                    Construct SIGDA1(DELTA,A) and SIGDA2(DELTA,A) which
C                 are used    SO_RES_C.
C--------------------------------------------------------------------
C
                     DTIME      = SECOND()
                     CALL SO_SIGDA(WORK(KSIGDA1),LSIGDA1,
     &                             WORK(KSIGDA2),LSIGDA2,
     &                             T2MP,LT2MP,WORK(KDSRHF),LDSRHF,
     &                             WORK(KBTR1E),LBTR1E,
     &                             WORK(KBTR1D),LBTR1D,
     &                             WORK(KBTJ1E),LBTJ1E,
     &                             WORK(KBTJ1D),LBTJ1D,
     &                             WORK(KCMO),LCMO,IDEL,ISYMD,ISYDIS,
     &                             ISYRES,ISYMTR,WORK(KEND5),LWORK5)
                     DTIME      = SECOND()   - DTIME
                     SOTIME(17) = SOTIME(17) + DTIME

                  ENDIF
C
  230          CONTINUE  ! End of IDEL2 loop
C
#ifdef VAR_MPI
               IF (loadbal_dyn) then
                   timefin = mpi_wtime()
                   itimeilll = ktiming + illl - 1
                   work(itimeilll) = work(itimeilll) + timefin - timeini
               ENDIF
#endif              
C
  220       CONTINUE !End of ILLL loop
C
  210    CONTINUE !end of ISYMD1 loop
C
C  888    continue !Continuation point in case of parallel run.
C====================================================
C        End of loop over distributions of integrals.
C====================================================
C
#ifdef VAR_MPI
C-------------------------------------------------
C        Communicate the result vectors to master.
C-------------------------------------------------
C  Note for further development: 
C  In the following, non-blocking reductions could be used,
C  since the one-particle result-vector is first needed in the
C  second call, sigai and sigda in the third and fourth call, 
C  and the two-particle vector is just written to file.
C  The alternative approach, is to do all these calculations
C  one the slaves, and only reduce in the end, this would save
C  communication of sigai and sigda.
C
C  Use the fact that memory are together
         lrestot = lres1e + lres1d + lres2e + lres2d
         lsigtot = lsigai1 + lsigai2 + lsigda1 + lsigda2
         latot   = laij + laab
!         write(*,*) lrestot, lsigtot, latot
C  Master use inplace operations       
         if (mynum .eq. 0 ) then
C
C           Fock-matrix
            call mpi_reduce( mpi_in_place, work(kfock), lfock, 
     &                       MPI_REAL8, MPI_SUM, 0, 
     &                       soppa_comm_active, ierr)
C
C           Result-vectors
            call mpi_reduce( mpi_in_place, work(kres1e), lrestot,
     &                       MPI_REAL8, MPI_SUM, 0, 
     &                       soppa_comm_active, ierr)
C
            if (singles_second) then
C
C           Sigai / Sigda
               call mpi_reduce( mpi_in_place, work(ksigai1), lsigtot,
     &                          MPI_REAL8, MPI_SUM, 0,
     &                          soppa_comm_active, ierr)
C
C           Aij / Aab -- only calculated in first pass  
               if ( inewtr .eq. 1 ) then
C               write(lupri,*) 'Aij'            
                  call mpi_reduce ( mpi_in_place, work(kaij), latot,
     &                              MPI_REAL8, MPI_SUM, 0,
     &                              soppa_comm_active, ierr)
               endif
            endif
         else
C  Slaves pass the same buffer as the recieve-buffer...
C
C           Fock-matrix
            call mpi_reduce( work(kfock), work(kfock), lfock, 
     &                       MPI_REAL8, MPI_SUM, 0, 
     &                       soppa_comm_active, ierr)
C
C           Result-vectors            
            call mpi_reduce( work(kres1e), work(kres1e), lrestot,
     &                       MPI_REAL8, MPI_SUM, 0, 
     &                       soppa_comm_active, ierr)
C
            if (singles_second) then
C              Sigai / Sigda
               call mpi_reduce( work(ksigai1), work(ksigai1), lsigtot,
     &                          MPI_REAL8, MPI_SUM, 0,
     &                          soppa_comm_active, ierr)
C
C              Aij / Aab            
               if ( inewtr .eq. 1 ) then
                  call mpi_reduce ( work(kaij), work(kaij), latot,
     &                              MPI_REAL8, MPI_SUM, 0,
     &                              soppa_comm_active, ierr)
               endif
            endif
C
C  After the reductions, the slaves are done; cycle loop
            goto 100
         endif
#endif
C
C---------------------------------------------
C        Transform AO Fock matrix to MO basis.
C---------------------------------------------
C
         DTIME      = SECOND()
         CALL CC_FCKMO(WORK(KFOCK),WORK(KCMO),WORK(KCMO),
     &                    WORK(KEND2),LWORK2,ISYRES,1,1)
         DTIME      = SECOND()   - DTIME
         SOTIME(24) = SOTIME(24) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA two-particle parts to the result 
C        vectors.
C------------------------------------------------------------------
C

         DTIME      = SECOND()
         CALL SO_TWOFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                   WORK(KFOCK),LFOCK,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(25) = SOTIME(25) + DTIME
C
         if (singles_second) then
C
C-------------------------------------------------------------------------
C           Add contribution from sigma1(alfa,i) to RES1E and from 
C           sigma2(alfa,i) to RES1D. I.e. the third terms in eqs. (34) and 
C           (35).
C-------------------------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_B(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                    WORK(KSIGAI1),LSIGAI1,WORK(KSIGAI2),LSIGAI2,
     &                    WORK(KCMO),LCMO,ISYRES)
            DTIME      = SECOND()   - DTIME
            SOTIME(18) = SOTIME(18) + DTIME
C
C-----------------------------------------------------------------------
C           Add contribution from sigda1(delta,a) to RES1E and from
C           sigda2(delta,a) to RES1D. I.e. the fourth terms in eqs. (34)
C           and (35).
C-----------------------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_C(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                    WORK(KSIGDA1),LSIGDA1,WORK(KSIGDA2),LSIGDA2,
     &                    WORK(KCMO),LCMO,ISYRES)
            DTIME      = SECOND()   - DTIME
            SOTIME(19) = SOTIME(19) + DTIME

C-----------------------------------------------------------------
C           Calculate and add the symmetry correcting term to A in
C           eq. (44).
C-----------------------------------------------------------------
C
            DTIME      = SECOND()

            CALL SO_RES_SYM(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                      WORK(KAIJ),LAIJ,WORK(KAAB),LAAB,WORK(KTR1E),
     &                      LTR1E,WORK(KTR1D),LTR1D,ISYRES)
            DTIME      = SECOND()   - DTIME
            SOTIME(20) = SOTIME(20) + DTIME
C
C---------------------------------------------------------
C           Calculate and add the Fock-term to A in eq. (40).
C---------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_FCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                      WORK(KTR1E),LTR1E,WORK(KTR1D),
     &                      LTR1D,FOCKD,LFOCKD,DENSIJ,LDENSIJ,DENSAB,
     &                      LDENSAB,ISYRES,ISYMTR)

            DTIME      = SECOND()   - DTIME
            SOTIME(21) = SOTIME(21) + DTIME
C
         endif
C
C------------------------------------------------------------------
C        Calculate and add the RPA one-particle parts to the result
C        vectors.
C------------------------------------------------------------------
C
         DTIME      = SECOND()

         CALL SO_ONEFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,FOCKD,
     &                   LFOCKD,WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                   ISYRES,ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(26) = SOTIME(26) + DTIME
C
         IF (DOUBLES) THEN
C
C-----------------------------------------------------------------
C        Construct D-contribution to 2p2h result vectors RES2E and
C        RES2D.
C-----------------------------------------------------------------
C
            DTIME      = SECOND()
            CALL SO_RES_CD(WORK(KRES2E),LRES2E,WORK(KRES2D),LRES2D,
     &                     WORK(KTR2E),LTR2E,WORK(KTR2D),LTR2D,
     &                     FOCKD,LFOCKD,ISYRES,WORK(KEND2),LWORK2)
         
            DTIME      = SECOND()   - DTIME
            SOTIME(30) = SOTIME(30) + DTIME
         ENDIF
#ifdef VAR_MPI
C Slaves are done
         IF (MYNUM .NE. 0) GOTO 100
#endif
C
C----------------------------------------
C        Write new result vectors to file.
C----------------------------------------
C
         CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
         CALL SO_WRITE(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,INEW)
         IF (DOUBLES) THEN
            CALL SO_WRITE(WORK(KRES2E),LRES2E,LURS2E,FNRS2E,INEW)
            CALL SO_WRITE(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,INEW)
         ENDIF
C
  100 CONTINUE 
C
C==================================
C     End of loop over excitations.
C==================================
C
#ifdef VAR_MPI
C-----------------------------------------------------------
C     Communicate the second order density matrix if needed.
C-----------------------------------------------------------
      IF ( (MODEL .EQ. 'AOSOP') .AND. 
     &     (NIT .EQ. 1)              ) THEN
         IF (MYNUM.EQ.0) THEN
            CALL MPI_REDUCE( MPI_IN_PLACE, DENSAI, LDENSAI, MPI_REAL8,
     &                       MPI_SUM, 0, SOPPA_COMM_ACTIVE, IERR)
         ELSE 
            CALL MPI_REDUCE( DENSAI, DENSAI, LDENSAI, MPI_REAL8,
     &                       MPI_SUM, 0, SOPPA_COMM_ACTIVE, IERR)
         ENDIF
      ENDIF
C--------------------------------------
C     Communicate the timings if needed
C--------------------------------------
      if ( loadbal_dyn ) then
         if (mynum .eq.0) then
C  Master
C --------           
C Recieve the timings            
            call mpi_reduce( mpi_in_place, work(ktiming), soppa_nint,
     &                       MPI_REAL8, MPI_SUM, 0, SOPPA_COMM_ACTIVE,
     &                       IERR)
C
C Redo the loadbalancing based on the timings and distribute them
            ksorted = kend2
            ktmp    = (soppa_num_active*maxnumjobs)
            knasjob = ksorted + (ktmp+1)/irat
            kswork  = knasjob + (soppa_num_active+1)/irat
            kendf   = kswork + soppa_num_active
C
            call dynloadbal_parsoppa( AssignedIndices, maxnumjobs,
     &                                work(ktiming), soppa_nint,
     &                                work(ksorted), work(knasjob),
     &                                work(kswork) )  
            !                    add empty work-space...
            !                    currently it can start at kend2
         else
C   Slave
C  ---------          
C  Send the timings            
            call mpi_reduce( work(ktiming), work(ktiming), soppa_nint,
     &                       MPI_REAL8, MPI_SUM, 0, SOPPA_COMM_ACTIVE,
     &                       IERR)
C  Recieve the new indices
            call mpi_scatter( AssignedIndices, maxnumjobs, mpi_integer,
     &                        AssignedIndices, maxnumjobs, mpi_integer,
     &                        0, soppa_comm_active, ierr)
         endif
      endif
C Slaves are done
      IF (MYNUM .NE. 0) RETURN
#endif
C----------------------------------------------------------------
C     Calculate the last part of the second order density matrix.
C----------------------------------------------------------------
C
      DTIME      = SECOND()
      IF ( (MODEL .EQ. 'AOSOP')   .AND.
     &     (NIT .EQ. 1)         )  THEN
         CALL SO_DENSAI2(DENSAI,LDENSAI,FOCKD,LFOCKD)
      END IF
C
      DTIME      = SECOND()   - DTIME
      SOTIME(41) = SOTIME(41) + DTIME
C
      IF ( IPRSOP. GE. 7 ) THEN
C------------------------------------------
C        Write new resultvectors to output.
C------------------------------------------
         DO 400 INEWTR = 1,NNEWTR
            INEW = NOLDTR + INEWTR
            WRITE(LUPRI,'(/,I3,A)') INEWTR,
     &                '. new E[2] linear transformed trial vector'
            CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
            CALL SO_READ(WORK(KRES2E),LRES2E,LURS2E,FNRS1D,INEW)
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KRES1E+I-1),WORK(KRES1D+I-1),I=1,LRES1E)
            IF (DOUBLES) THEN
               CALL SO_READ(WORK(KRES1D),LRES1D,LURS1D,FNRS2E,INEW)
               CALL SO_READ(WORK(KRES2D),LRES2D,LURS2D,FNRS2D,INEW)
               WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &             (I,WORK(KRES2E+I-1),WORK(KRES2D+I-1),I=1,LRES2E)
            ENDIF
  400    CONTINUE
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C   
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
C
      IF (DOUBLES) THEN
         CALL SO_CLOSE(LUTR2E,FNTR2E,'KEEP')
         CALL SO_CLOSE(LUTR2D,FNTR2D,'KEEP')
         CALL SO_CLOSE(LURS2E,FNRS2E,'KEEP')
         CALL SO_CLOSE(LURS2D,FNRS2D,'KEEP')
      ENDIF
C
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_ERES')
C
      RETURN
      END SUBROUTINE

C
#ifdef VAR_MPI
C
C
      subroutine dynloadbal_parsoppa( localIndices, maxnumjobs,
     &                                timings, ltimings,
     &                                sortedindices, numassignjobs,
     &                                sumofwork)
C    Dynamic Load Balancing for the parallel SOPPA calculations (and parallel RPA). 
C    The routine assumed that the MPI processes that does the load balancing is the master. If any other routine enters, there will be issues with the updated ILLL indices in the AssignedIndices array.  Right now, any slave that enters is immediately evicted from the routine.
C
C The routine takes an array of timings as input. The timings are creatd on the fly by every parallel process and the index corresponds to the ILLL index in the parallel calculation. The time is given in integers from calling system_clock. The timings are not given in any human-readable form. This choice was made so that one can reuse the subroutine getallocsize, which relies on an integer input array. 
C The work is resorted once the actual time associated with an ILLL index is known and the indices are rebalanced once more and sent to all slaves.
C The improvement over the presorted balancing scheme is expected to be very small, but gets better the larger the basis set.
C 
C F.Beyer Oct. 2014.

      use so_parutils, only : soppa_comm_active, soppa_num_active,
     &                        soppa_nint

#include "implicit.h"
#include "mpif.h"
!#include"parsoppa.h" !fetches ntot
#include "maxorb.h"
#include "distcl.h"
#include "priunit.h"

      ! Dummy parameters
      double precision, dimension(ltimings), intent(inout) :: timings
      integer, intent(inout) :: localIndices(maxnumjobs)
      integer  :: sortedindices(maxnumjobs,soppa_num_active)
      integer  :: numassignjobs(soppa_num_active)
      dimension :: sumofwork(soppa_num_active)
      integer  :: ltimings, maxnumjobs

      ! Bookkeeping
      integer :: maxrows, maxcols, numprocs, numrecipients
      integer :: getnumjobs, targetID, myid
      integer :: colindex, rowindex, assignILLL, col,i
      integer :: ntot
  
      ntot = soppa_nint
      ! In case the amount of work is smaller than
      ! the number of MPI processes...
      maxcols = soppa_num_active

      ! Explicitly set the maximum number of jobs a single MPI process can be allocated.
      maxrows = maxnumjobs

      sortedindices(:,:) = 0
      numassignjobs(:) = 0
      sumofwork(:) = 0.0D0


      DO i=1, ntot
              ! Find the largest chunk of available work
              assignILLL = maxloc(timings,dim=1, mask=timings.ge.0.0D0)

              ! Find the laziest slave
              colindex = minloc(sumofwork, dim=1)

              ! Write the ILLL index to the correct row in the sortedmatrix
              numassignjobs(colindex) = numassignjobs(colindex) + 1
              rowindex = numassignjobs(colindex)
              sortedindices(rowindex, colindex) = assignILLL
              
              ! Update the slave's expected work/walltime and
              ! remove the work from the timings array
              sumofwork(colindex) = sumofwork(colindex) 
     &                             +timings(assignILLL)
              timings(assignILLL) = -1.0D0

      ENDDO

C     PRINT HOW WORK IS DISTRIBUTED
      write (LUPRI,'(a)') 'AOSOPPA Work distribution'
      write (LUPRI,'(a)') 'NODE       Expected time '
      do i = 1, maxcols
           write (LUPRI, '(i5,f20.5)') i, sumofwork(i) 
      enddo 
      ! Send the info to every slave that does computation (some slaves might be stalled in the polling barrier in case there are too many MPI processes compared to the number of tasks).
      call mpi_scatter( sortedindices, maxnumjobs, mpi_integer, 
     &                  localIndices, maxnumjobs, mpi_integer,
     &                  0, soppa_comm_active, ierr) 

      ! Update the master's own array of ILLL indices
      ! should be done by the scatter

      return
      end subroutine






      subroutine presortloadbal_parsoppa(localIndices, maxnumjobs,
     &           indexb, 
     &           work, lwork) 
C     This subroutine load balances parallel SOPPA/RPA 
C     calculation.
C
C     The routine makes a best guess of loadbalancing by giving 
C     every MPI process an equal number of distributions to handle.
C     Testing shows this is a better first guess that giving every
C     MPI process an equal number of ILLL indices to work with 
C     (since there is a different number of distributions 
C     associated with every ILLL index).
C
C     F.Beyer Oct. 2014.

      use so_parutils, only: soppa_comm_active, soppa_nint,
     &                       soppa_num_active

C#include "implicit.h"
#include "priunit.h"
#include "mpif.h"
C#include "parsoppa.h" ! to get the value for NTOT and KINDXB
#include "maxorb.h"
#include "iratdef.h"
C fetch HERDIR
#include "ccsdinp.h"
C     Dummy parameters
      integer, intent(in) :: lwork, maxnumjobs
      double precision :: work(lwork) !intent in
      integer :: indexb(*)
      integer, dimension(maxnumjobs), intent(out) :: localIndices

C     Pre-sorting load balancing variables
      integer :: getnumjobs, getindices, myid, numprocs, ierr            
C      integer :: presortarray, finalsorted
      integer :: col, ntot
      integer :: kpresortarray, kend, kfinalsorted, ktmp

      ntot = soppa_nint
      
      ! numprocs = soppa_num_active !! not consistent with other use 
!      call mpi_comm_size(soppa_comm_active, numprocs, ierr)  
      numprocs = soppa_num_active


      IF (.NOT. HERDIR) THEN
C     ERI code

C     FIND THE AMOUNT OF WORK ASSOCIATED WITH EVERY AO INDEX
         kpresortarray = 1
         kend = kpresortarray + 2* ntot / irat      
         call presortaodist(ntot, indexb, work(kpresortarray) )


C     CREATE THE FINAL MATRIX OF PRE-SORTED AO INDICES
         kfinalsorted = kend 
         ktmp = kfinalsorted + (maxnumjobs * numprocs + 1) /irat
         kend = ktmp + 2*numprocs / irat
         if( kend .gt. lwork ) then
            call quit('Insufficient memory in presorted loadbalancer!')
         endif
         call partitionAOindices(ntot, maxnumjobs, numprocs, 
     &              work(kpresortarray), work(kfinalsorted),
     &              work(ktmp) )
      else
C     HERDIR code
         kfinalsorted = 1
         kend = kfinalsorted + ( maxnumjobs * numprocs + 
     &             mod(maxnumjobs*numprocs,irat) )/irat
         if( kend .gt. lwork ) then
            call quit('Insufficient memory in presorted loadbalancer!')
         endif
         call herdir_presort( work(kfinalsorted) ) 
      endif


C     TRANSFER AO INDICES TO SLAVES
C     Use scatter 
      call mpi_scatter( work(kfinalsorted), maxnumjobs, mpi_integer,
     &                  localIndices, maxnumjobs, mpi_integer,
     &                  0, soppa_comm_active, ierr )

      return
      end subroutine 


      subroutine herdir_presort ( sorted, maxnumjobs )
C     This is a subroutine for doing the initial distribution of the    
C     AO-indices in a parallel SOPPA calculation.
C     This subroutine works with Hermite (HERDIR) integral generator,
C     as the code by F. Beyer only supports using ERI
C
      use so_parutils, only: soppa_num_active, soppa_nint
      implicit none
      integer, intent(out) :: sorted( maxnumjobs, soppa_num_active )
      integer, intent(in)  :: maxnumjobs
      integer              :: intnum, inext, icol, inum

      sorted = 0
      inext = 0
      icol = 1
C     For now just do it the stupid way ( round robin )

C     initial implementation, is this what chokes ifort?
C      do intnum = 1, soppa_nint
C         inext = inext + 1
C         if ( inext .gt. soppa_num_active ) then
C             inext = 1
C             icol = icol + 1
C         endif
C         sorted ( icol, inext ) = intnum
C      enddo

C     Distribute first an even number      
      inum = soppa_nint/soppa_num_active
      do icol = 1, soppa_num_active
         do inext = 1, inum
            sorted( inext, icol ) = inext + inum*(icol-1)
         enddo
      enddo
C     Distribute the remainder
      do icol = 1, mod(soppa_nint, soppa_num_active)
         sorted ( inum+1, soppa_num_active-icol+1) =
     &            inum*soppa_num_active + icol       
      enddo
C     Debug print
      do inext = 1, soppa_num_active 
         print '(10i3)', sorted(:,inext)
      enddo
      return
      end subroutine 


      
      subroutine pollingbarrier(pollinginterval)
C     This subroutine is implemented in case there is ever a need to 
C     remove certain processes from a calculation.
C
C     The processes will repeatedly poll until they receive a non-blocking 
C     send from the master with a logical value equal to .true.
C
C     F.Beyer Oct. 2014
      
#include "implicit.h"
#include "mpif.h"

      integer, dimension(MPI_STATUS_SIZE) :: mpistatus
      integer, intent(in)                 :: pollinginterval !in milliseconds
      integer                             :: myid, request
      logical                             :: flag, exitbarrier
      volatile                            :: exitbarrier

      call mpi_comm_rank(mpi_comm_world, myid, ierr)
      print *, 'I am in the polling barrier ', myid
           
      ! Initiate the polling variables and ask for the non-blocking update
      exitbarrier = .false. 
      call mpi_irecv(exitbarrier, 1, MPI_LOGICAL, 0, myid,
     &               mpi_comm_world, request, ierr) 

      ! Polling barrier, the process will cycle repeatedly until released.
130   continue 
#ifdef VAR_IFORT
      call sleepqq(pollinginterval)  
#endif
      call mpi_test(request, flag, mpistatus, ierr)
      if (.not. flag) goto 130
      ! Warning, seems that we'll go into an infinite loop,
      ! if exitbarrier is ever sent as .false.
      if (.not.exitbarrier) then
         goto 130 ! Cycle to the top of the barrier.
      else
         return
      endif

      return
      end subroutine

C /* deck presortaodist*/
      subroutine presortaodist(Nindex, indxbt, outlist)
C A subroutine associated with the atomic integral parallel RPA/SOPPA calculations.
C Pre-calculate IDEL2 indexes before starting a parallel calculation.
C In other words, this subroutine calculates how many distributions are associated with the calculation.
C
C This routine assembles a matrix that counts the number of AOs
C associated with an ILLL distribution index. This array is used by getallocsize
C and partitionAOindices to pre-sort the integrals that need to be done 
C
C The first row in outlist is the number of distributions
C The second row in outlist is the associated ILLL index
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "eridst.h"

      integer,                       intent(in) :: Nindex
      integer, dimension(*),         intent(in) :: indxbt
      integer, dimension(2, Nindex), intent(out):: outlist
      integer :: i

      do i=1, Nindex
         call getdst(i, 0, 0)
         call pickao(0)
         call eridsi(indxbt, 0)
         outlist(1, i) = ndistr
         outlist(2, i) = i 
      enddo

      return

      end subroutine


C     /* deck getallocsize */
      SUBROUTINE getallocsize(ntot, originalsort, maxnumjobs)
C A subroutine associated with the atomic integral parallel RPA/SOPPA calculations.
C 
C This subroutine is used to get the amount of storage that needs
C to be allocated for a parallel SO_ERES run.
C 
C The subroutine calculates which process will be assigned the
C most single jobs (not the largest total amount of work) for a
C parallel RPA/SOPPA calculation. Its output is used to allocate 
C the right amount of storage by the master when it starts pre-sorting
C the integrals for the parallel calculation of the E matrix.
C

#include "implicit.h"
#include "mpif.h"
      integer,                     intent(in)  :: ntot
      integer, dimension(2, ntot), intent(in)  :: originalsort
      integer,                     intent(out) :: maxnumjobs

      integer, dimension(:,:), allocatable  :: copysort
      integer, dimension(:,:), allocatable  :: sumofwork
      integer, dimension(2) :: temploc, tempwork, tempout
      integer :: allocstatus, deallocstatus, numprocs



      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      allocate( copysort(2, ntot), sumofwork(2, numprocs)
     &         ,stat=allocstatus)
      if(.not.(allocstatus.eq.0) ) then
         call quit('Allocation error in GETALLOCSIZE')
      endif

      call izero(sumofwork, (2*numprocs) )
      copysort = originalsort

      DO i=1, ntot
         ! Find location of largest chunk of work and the work itself
         temploc = maxloc(copysort, DIM=2, mask=copysort.gt.0) 
         addwork = copysort(1, temploc(1) )
         copysort( 1,temploc(1) ) = 0

         ! Find laziest slave and simulate the workload on the slave
         tempwork = minloc(sumofwork, DIM=2) 
         sumofwork(1, tempwork(1)) = sumofwork(1, tempwork(1)) + addwork! adding total number of distributions
         sumofwork(2, tempwork(1)) = sumofwork(2, tempwork(1)) + 1  !adding total number of assigned indexes
      ENDDO

      tempout = maxloc(sumofwork, dim=2)
      maxnumjobs = sumofwork( 2,tempout(2) )

      deallocate(copysort, sumofwork, stat=deallocstatus)
      if(.not.(deallocstatus.eq.0) ) then
         call quit('Deallocation error in GETALLOCSIZE')
      endif

      return

      END SUBROUTINE


C     /* deck partitionAOindices */
      SUBROUTINE  partitionAOindices(ntot, rows, cols, presortedarray,
     &                               sorted, sumofwork)
C A subroutine associated with the atomic integral parallel RPA/SOPPA calculations.
C
C The output from this routine is the array 'sorted'.
C The sorted matrix contains AO integral indexes for two-electron integrals.
C Every column contains a list of ILLL indexes that are to be assigned as work 
C to a single process. The total amount of work per column is estimated based
C on the number of distributions that are related to the ILLL indexes in the column. 
C
C This pre-sorting of ILLL indexes for a parallel calculation approximates an
C even distribution of total work for all processes when performing two-electron 
C integrals in parallel for AOSOPPA and AORPA. 
C
      !use mpi
#include "implicit.h"
#include "mpif.h"

      integer,                        intent(in)    :: ntot, rows, cols
      integer, dimension(2, ntot),    intent(inout) :: presortedarray
      integer, dimension(rows, cols), intent(out)   :: sorted

      integer, dimension(2, cols), intent(out) :: sumofwork
      integer, dimension(2) :: tempavail, templazy
      integer :: numprocs, availloc, targetrow
      integer :: allocstatus, deallocstatus

      numprocs = cols
 
      call izero(sorted, (rows*cols) )
      call izero(sumofwork, (2*numprocs) )

      DO i=1, ntot
C        FIND LARGEST CHUNK OF AVAILABLE WORK AND ITS INDEX
         tempavail = maxloc(presortedarray, DIM=2, 
     &                      mask=presortedarray.gt.0) 
         availloc = tempavail(1)
         numdists = presortedarray(1, availloc) ! amount of available work
         aoindex  = presortedarray(2, availloc) ! the index to be passed to process

         presortedarray(1, availloc) = 0
         presortedarray(2, availloc) = 0

C        FIND THE LAZIEST PROCESS, GIVE IT WORK & INCREMENT THE ROW COUNTER
         templazy = minloc(sumofwork, DIM=2) 
         lazyloc = templazy(1) ! This is equal to: MYID+1 
         sumofwork( 1,lazyloc ) = sumofwork( 1,lazyloc ) + numdists 
         sumofwork( 2,lazyloc ) = sumofwork( 2,lazyloc ) + 1 
         targetrow = sumofwork( 2,lazyloc )

C        ADD AO-INDEX TO FIRST AVAILABLE ROW IN THE LAZY SLAVE'S COLUMN
         sorted( targetrow,lazyloc ) = aoindex
      ENDDO
      
      return

      END SUBROUTINE


#endif
!VAR_MPI at the beginning of dynloadbal_parsoppa
