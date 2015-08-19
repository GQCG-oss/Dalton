C
C  /* Deck rp_eres */
      SUBROUTINE RP_ERES(NOLDTR, NNEWTR, FOCKD, LFOCKD, ISYMTR, 
     &                   WORK,   LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     PFP & SPAS, November 2013: triplet excitation energies
C
C     PURPOSE: Driver routine for making a linear transformation of
C              a trialvector with the SOPPA hessian matricx E[2]. 
C              The trial vector consists of two parts TR1E and TR1D. 
C              E refers to excitations and D to 
C              de-excitations. 1 refer to the one-particle part and
C              2 to the two-particle part. The linear transformed 
C              trialvector is refered to as the resultvector and is
C              kept in four corresponding arrays. For the linear
C              transformation with E[2] the result vector is in RES1E,
C              RES1D.
C              The linear transformation is driven over atomic orbitals,
C              and E[2] is not constructed explicitly.
C
C >>> RP_ERES - Data dictionary: <<<
C NOLDTR: Number of Old trial vectors
C NNEWTR: Numver of new trial vectors
C FOCKD: The Fock Matrix
C LFOCKD: Number of entries in the Fock Matrix
C ISYMTR:
C WORK: temporary work-vector 
C LWORK: Total needed length of the work-vector.
#include "implicit.h"
#ifdef VAR_MPI
         !use mpi
#include "mpif.h"
#endif
#include "maxorb.h" 
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "eritap.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "soppinf.h"

#ifdef VAR_MPI
#include "iprtyp.h"
#include "infpar.h"
#include "parsoppa.h"
      integer :: numprocs
      call mpi_comm_size(mpi_comm_world, numprocs, ierr) 
      copyisymtr = isymtr
      forceupdate = .true. !Implicitly saved. Do NOT move the initialization to the main block of code.
#endif

C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_ERES')
C
C------------------------------------------------------------------
C     Determine the symmetry of the result vector from the symmetry
C     of the trial vector ISYMTR, and the operator symmetry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      LCMO    = NLAMDT
C
      KCMO    = 1
      KEND1   = KCMO   + LCMO  
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_ERES.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_ERES.1',' ',KEND1,LWORK)
C
C-------------------------------------------------------
C     Get the matrix which contains the MO coefficients.
C-------------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_GETMO(WORK(KCMO),LCMO,WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(1)  = SOTIME(1) + DTIME
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
C
      KTR1E   = KEND1
      KTR1D   = KTR1E   + LTR1E
      KRES1E  = KTR1D   + LTR1D
      KRES1D  = KRES1E  + LRES1E
      KFOCK   = KRES1D  + LRES1D
      KDENS   = KFOCK   + LFOCK
      KBTR1E  = KDENS   + LDENS
      KBTR1D  = KBTR1E  + LBTR1E
      KBTJ1E  = KBTR1D  + LBTR1D
      KBTJ1D  = KBTJ1E  + LBTJ1E
      KEND2   = KBTJ1D  + LBTJ1D
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('RP_ERES.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('RP_ERES.2',' ',KEND2,LWORK)
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
C
      IF ( IPRSOP. GE. 7 ) THEN !This in merely printing related.
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
         DO 50 INEWTR = 1,NNEWTR
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
            INEW = NOLDTR + INEWTR
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
            CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
            WRITE(LUPRI,'(/,I3,A)') 
     &            INEWTR,'. new trial vector in RP_ERES'
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
   50    CONTINUE
C
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
C----------------------------------------
C        Initialize RES1E, RES1D and FOCK
C----------------------------------------
C
         CALL DZERO(WORK(KRES1E),LRES1E)
         CALL DZERO(WORK(KRES1D),LRES1D)
         CALL DZERO(WORK(KFOCK),LFOCK)
C
C--------------------------
C        Read trial vector.
C--------------------------
C
         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
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
         CALL SO_BCKTR(WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,WORK(KBTR1E),
     &                 LBTR1E,WORK(KBTR1D),LBTR1D,WORK(KBTJ1E),LBTJ1E,
     &                 WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,ISYMTR)
         DTIME     = SECOND()  - DTIME
         SOTIME(7) = SOTIME(7) + DTIME
C
C=======================================================
C        Start the loop over distributions of integrals.
C=======================================================
C
         IF (DIRECT) THEN
            NTOSYM = 1
            DTIME     = SECOND()
            IF (HERDIR) THEN
               CALL HERDI1(WORK(KEND2),LWRK2,IPRINT)
            ELSE
               KCCFB1 = KEND2
               KINDXB = KCCFB1 + MXPRIM*MXCONT
               KEND2  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
               LWORK2  = LWORK  - KEND2

               CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                     KODPP1,KODPP2,KRDPP1,KRDPP2,KFREE,LFREE,
     &                     KEND2,WORK(KCCFB1),WORK(KINDXB),WORK(KEND2),
     &                     LWORK2,IPRINT)

               KEND2  = KFREE
               LWORK2 = LFREE
               DTIME     = SECOND()  - DTIME
               SOTIME(8) = SOTIME(8) + DTIME
            ENDIF
         ELSE
            NTOSYM = NSYM
         ENDIF
C
         KENDSV  = KEND2
         LWORKSV = LWORK2
C
         ICDEL1 = 0
         DO 210 ISYMD1 = 1,NTOSYM
C
            IF (DIRECT) THEN  
               NTOT = MXCALL
            ELSE
               NTOT = NBAS(ISYMD1)
            ENDIF
C
#ifdef VAR_MPI
      ! Start a parallel calculation of the fock matrix in case more than 1 process is available.
      ! The result is kept in the master process' work array at work(kfock)
      if (numprocs.gt.1) then
          call mpixbcast(85, 1, 'INTEGE', 0)
          call mpixbcast(0, 1, 'INTEGE', 0)
          CALL PAR_RP_ERES(WORK, LWORK, 0)
          goto 888
      endif
#endif
            DO 220 ILLL = 1,NTOT
C
C------------------------------------------------
C              If direct calculate the integrals.
C------------------------------------------------
C
C
C
               IF (DIRECT) THEN
C
                  KEND2  = KENDSV
                  LWORK2 = LWORKSV
C
                  DTIME     = SECOND()
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
                     DTIME     = SECOND()  - DTIME
                     SOTIME(9) = SOTIME(9) + DTIME
                  ENDIF
C
                  LRECNR = ((NBUFX(0) - 1)/IRAT) +1
C
                  KRECNR  = KEND2
                  KEND2   = KRECNR + LRECNR
                  LWORK2  = LWORK  - KEND2
C
                  CALL SO_MEMMAX ('RP_ERES.2B',LWORK2)
                  IF (LWORK2 .LT. 0) 
     &               CALL STOPIT('RP_ERES.2B',' ',KEND2,LWORK)
C
               ELSE
                  NUMDIS = 1
               ENDIF
C
C------------------------------------------------------------------------
C  Loop over number of distributions in disk.
C  In the case of ERI there are more than one distribution and 
C  IDEL2 loops over them and the actual index of the delta orbital IDEL is 
C  then obtain from the array INDEXA.
C  In the case of a not direct calculation there is only one distribution
C  on the disk, which implies that IDEL2 is always 1 and that IDEL is 
C  systematically incremented by one each time.
C------------------------------------------------------------------------
C
               DO 230 IDEL2 = 1,NUMDIS
C
                  IF (DIRECT) THEN
                     IDEL  = INDEXA(IDEL2)
                     ISYMD = ISAO(IDEL)
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
                  KXINT   = KEND2
                  KEND3   = KXINT  + LXINT
                  LWORK3  = LWORK  - KEND3
C
                  CALL SO_MEMMAX ('RP_ERES.3',LWORK3)
                  IF (LWORK3 .LT. 0) 
     &               CALL STOPIT('RP_ERES.3',' ',KEND3,LWORK)
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
C---------------------------------------------
C                 Work space allocation no. 4.
C---------------------------------------------
C
                  LDSRHF = NDSRHF(ISYMD)
C
                  KDSRHF  = KEND3
                  KEND4   = KDSRHF + LDSRHF
                  LWORK4  = LWORK  - KEND4
C
                  CALL SO_MEMMAX ('RP_ERES.4',LWORK4)
                  IF (LWORK4 .LT. 0) 
     &               CALL STOPIT('RP_ERES.4',' ',KEND4,LWORK)
C
C
C----------------------------------------------
C                 Calculate the AO-Fock matrix.
C----------------------------------------------
C
                  DTIME      = SECOND()
                  IF (TRIPLET) THEN
C
                     CALL CC_AOFOCK3(WORK(KXINT),WORK(KDENS),
     &                               WORK(KFOCK),WORK(KEND4),
     &                               LWORK4,IDEL,ISYMD,ISYMTR)
C
                  ELSE
C
                     CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),
     &                              WORK(KFOCK),WORK(KEND4),
     &                              LWORK4,IDEL,ISYMD,.FALSE.,
     &                              'CrashifNecessary',ISYMTR)
C Jan 30. 2014. F.Beyer:
C Ive altered the penultimate variable to 'CrashifNecessary' to make sure CC_aofock crashes
C If the boolean is ever true for any reason.
C If you use .TRUE. instead of .FALSE. in the call to CC_AOFOCK then make sure you 
C correct the variable that follows. The last variable that was in place before I changed it to 'CrashifNecessary' was uninitialized
C and could potentially cause a lot of havoc if the calculation didn't fail every time.
                  END IF
                  DTIME      = SECOND()   - DTIME
                  SOTIME(11) = SOTIME(11) + DTIME
C
  230          CONTINUE ! End of IDEL2 loop.
  220       CONTINUE ! End of ILLL loop.
  210    CONTINUE

  888    continue !Continuation point in case the calculation is done in parallel.

C
C====================================================
C        End of loop over distributions of integrals.
C====================================================
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
C----------------------------------------
C        Write new resultvectors to file.
C----------------------------------------
C
         CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
         CALL SO_WRITE(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,INEW)
C
  100 CONTINUE ! End of loop over number of excitations
C
C==================================
C     End of loop over excitations.
C==================================
C
      IF ( IPRSOP. GE. 7 ) THEN
C
C------------------------------------------
C        Write new resultvectors to output.
C------------------------------------------
C
         DO 400 INEWTR = 1,NNEWTR
C
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
C
            INEW = NOLDTR + INEWTR
C
            WRITE(LUPRI,'(/,I3,A)') INEWTR,
     &         '. new E[2] linear transformed'//
     &         ' trial vector'
C
            CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
            CALL SO_READ(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,INEW)
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KRES1E+I-1),WORK(KRES1D+I-1),I=1,LRES1E)
C
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
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_ERES')
C
      RETURN
      END SUBROUTINE










#ifdef VAR_MPI
C  /* Deck PAR_RP_ERES */
      SUBROUTINE PAR_RP_ERES(WORK, LWORK, MYID)
C PAR_RP_ERES is a Parallel subroutine for RP_ERES.F
C It calculates the AO-FOCK matrix in parallel by distributing the 
C two-electron integrals over several processes and assembling the 
C partial Fock matrices from every process before returning.
C
C WORK is the array where master and every slave assembles the fock matrix 
C incrementally. This is the 'output' of this subroutine.
C
C The load balancing is done by examining the number of distributions associated 
C with every AO index in the two-electron integrals, before the integrals are 
C performed. The total load is calculated and the work is spread as evenly as 
C possible among the master and all slaves. The AO indices with the biggest 
C associated load of distributions is assigned first; then the smaller indices 
C are assigned.
C
C A large number of common blocks need to be transferred from the master to 
C the slaves by using the subroutine getbytespan and mpi_bcast. 
C
C                       DATA DICTIONARY:
C WORK intent(inout) temporary storage for several calculations. The master
C                    returns from par_rp_eres with the correct result stored
C                    in the work vector and processes the result in rp_eres.
C LWORK intent(in)   Length of work vector.
C MYID intent(in)    Process ID in MPI_COMM_WORLD.

         !use mpi
      use dyn_iadrpk

#ifdef VAR_IFORT
      use IFPORT, ONLY: SLEEPQQ
#endif
! Parameter-only include files
#include "implicit.h"

#include "mpif.h"

#include "priunit.h" 
#include "maxorb.h" 
#include "maxash.h"
#include "mxcent.h"
#include "iratdef.h"
#include "iprtyp.h" 
#include "maxaqn.h" 
#include "chrnos.h" 
#include "ibtpar.h"
!#include "ibtfun.h"

! These include files depend on previous include files
#include "infpar.h" 
#include "eritap.h" 
#include "eribuf.h" 

#include "inftap.h"
#include "ccorb.h"
#include "infind.h" 
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h" 
#include "soppinf.h"
#include "parsoppa.h"
#include "aobtch.h"
#include "odclss.h"
#include "ccom.h" 
#include "ericom.h" 
#include "eridst.h" 
#include "erithr.h"   
#include "erimem.h"  
#include "odbtch.h"
#include "nuclei.h"
#include "symmet.h"
#include "r12int.h" 
#include "hertop.h"

C Common blocks from  PICKAO
#include "cbirea.h"
#include "erisel.h"
#include "symsq.h"
#include "gnrinf.h"
#include "ccpack.h"
#include "ccinftap.h"

C     Intent of Dummy arguments
      integer, intent(in) :: myid, lwork
      dimension work(lwork)  !intent(inout)

C     Variables for updating the slaves' /common/ blocks.
      integer  :: bytesize, allocstatus, deallocstatus
      integer  :: isymtr

C     Pre-sorting load balancing variables
      logical, save :: forceloadbalance        =.true.
      logical, save :: forcedynamicloadbalance =.false.
      logical, save :: slavesupdatetasks       =.false.
      integer, save :: maxnumjobs
      logical       :: assignnow, workreceived
      integer       :: getnumjobs, getindices
      integer, allocatable, dimension(:), save :: AssignedIndices

C     Timing variables used for load monitoring
      integer :: color
      integer :: counti, countf, count_rate
      real    :: timespent
      real, dimension(:), allocatable, save :: timings

C     MPI related and load balancing variables
      integer :: numprocs, request, allocsize
      integer, dimension(ntot)            :: temparray
      integer, dimension(MPI_STATUS_SIZE) :: mpistatus

      integer, allocatable, dimension(:,:) :: presortarray !only used by master
      integer, allocatable, dimension(:,:) :: finalsorted !only used by master
      real*8,  allocatable, dimension(:)   :: tempfock
      integer :: workindex
      
C     Miscellaneous
      DIMENSION INDEXA(MXCORB)

      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      call mpi_bcast(forceupdate,1,mpi_logical,0,mpi_comm_world, ierr)

C**********************************************************
C     Update all necessary common blocks
C     This large sequence of bcasts will only execute
C     once per call to RP_ERES
C
C     The <name> of the common block can be read from the <name>LAST
C     variable in the getbytespan call.
C**********************************************************
      commonblocktransfers: if (forceupdate) then

      call getbytespan(lbuf, eribufLAST, bytesize)
      call mpi_bcast(lbuf, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(luaorc, eritapLAST, bytesize)
      call mpi_bcast(luaorc, bytesize, mpi_byte,0, mpi_comm_world, ierr)

      call getbytespan(nsym, ccorbLAST, bytesize)
      call mpi_bcast(nsym, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(irow, infindLAST, bytesize)
      call mpi_bcast(irow, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(centsh, blocksLAST, bytesize)
      call mpi_bcast(centsh, bytesize, mpi_byte,0, mpi_comm_world, ierr)

      call getbytespan(skip, ccsdgninpLAST, bytesize)
      call mpi_bcast(skip, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(ccs, ccmodelsLAST, bytesize)
      call mpi_bcast(ccs, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(etmp, etmpLAST, bytesize)
      call mpi_bcast(etmp, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(nckijmax, ccsdmaxLAST, bytesize)
      call mpi_bcast(nckijmax, bytesize, mpi_byte,0,mpi_comm_world,ierr)

      call getbytespan(nt1amx, ccsdsymLAST, bytesize)
      call mpi_bcast(nt1amx, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(it2del, ccsdioLAST, bytesize)
      call mpi_bcast(it2del, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(mxcall, distclLAST, bytesize)
      call mpi_bcast(mxcall, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(runeri, cbieriLAST, bytesize)
      call mpi_bcast(runeri, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(sotime, soppinfLAST, bytesize)
      call mpi_bcast(sotime, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nexci2, soppexcLAST, bytesize)
      call mpi_bcast(nexci2, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(soorwc, rwinfLAST, bytesize)
      call mpi_bcast(soorwc, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(ntot, parsoppaLAST, bytesize)
      call mpi_bcast(ntot, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(expbt, aobtchLAST, bytesize)
      call mpi_bcast(expbt, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(nitcl, odclssLAST, bytesize)
      call mpi_bcast(nitcl, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(thrs, ccomLAST, bytesize)
      call mpi_bcast(thrs, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(gtotyp, ccomcLAST, bytesize)
      call mpi_bcast(gtotyp, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(scrmab, ericomLAST, bytesize)
      call mpi_bcast(scrmab, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(ndistr, eridstLAST, bytesize)
      call mpi_bcast(ndistr, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(thrsh, erithrLAST, bytesize)
      call mpi_bcast(thrsh, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(memadd, erimemLAST, bytesize)
      call mpi_bcast(memadd, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nodbch, odbtchLAST, bytesize)
      call mpi_bcast(nodbch, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(charge, nucleiLAST, bytesize)
      call mpi_bcast(charge, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(namn, nuclecLAST, bytesize)
      call mpi_bcast(namn, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(maxrep, symmtiLAST, bytesize)
      call mpi_bcast(maxrep, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(fmult, symmtrLAST, bytesize)
      call mpi_bcast(fmult, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(gamac, comr12LAST, bytesize)
      call mpi_bcast(gamac, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(mbas1, cmmmulLAST, bytesize)
      call mpi_bcast(mbas1, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(jtop, hertopLAST, bytesize)
      call mpi_bcast(jtop, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

      call getbytespan(zcmval, cbireaLAST, bytesize)
      call mpi_bcast(zcmval, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(mulnam, cbirea_cLAST, bytesize)
      call mpi_bcast(mulnam, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(nmulbs, cmmbasLAST, bytesize)
      call mpi_bcast(nmulbs, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(i2bst, symsqLAST, bytesize)
      call mpi_bcast(i2bst, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call getbytespan(thrpckint, ccpackLAST, bytesize)
      call mpi_bcast(thrpckint,bytesize,mpi_byte,0,mpi_comm_world, ierr)

      call getbytespan(luiajb, cc_tapLAST, bytesize)
      call mpi_bcast(luiajb,bytesize,mpi_byte,0,mpi_comm_world, ierr)

      call getbytespan(gradml, gnrinfLAST, bytesize) 
      call mpi_bcast(gradml, bytesize, mpi_byte, 0, mpi_comm_world,ierr)

      call stupid_isao_bcast_routine()

      ! Common blocks with explicitly set sizes, no need to calculate bytespan
      call mpi_bcast(basdir,612, mpi_character, 0, mpi_comm_world, ierr)
      call mpi_bcast(fnvajkl, 10, mpi_character, 0, mpi_comm_world,ierr)
      call mpi_bcast(vclthr, 24, mpi_byte, 0, mpi_comm_world, ierr)

      ! The slaves need to crate the iadrpk array by 
      ! calling the module functions rather than with a bcast
      call mpi_bcast(iadrpk_dim, 1, mpi_integer, 0,mpi_comm_world,ierr) 
      if (.not. allocated(iadrpk) ) then
         call get_iadrpk(lupri,nsym,muld2h,nbas,
     &                   nbast,i2bst,iaodis,iaodpk)
      endif

      call mpi_bcast(work, kendsv, mpi_real8, 0, mpi_comm_world, ierr)

      forceupdate = .false. ! This avoids redundant updating of common blocks. RP_ERES sets this to true when necessary.

      endif commonblocktransfers


C     This data needs to be updated every time PAR_RP_ERES is called 
      call mpi_bcast(work(kdens),ldens,mpi_real8,0,mpi_comm_world, ierr)
      call mpi_bcast(isympar, 1, mpi_integer, 0, mpi_comm_world, ierr)
      isymtr = copyisymtr 
      call dzero(work(kfock), lfock)
      

C**************************************************    
C     Start load balancing
C**************************************************    

C   In case there are too many processes compared
C   to the number of ILLL indices, the surplus processes are
C   redirected. They will return from the polling barrier when
C   the master sends a non-blocking logical equal to .true.
      if (myid.ge.ntot) color = 1
      if (myid.lt.ntot) color = 0 !these are the ones performing timings, other processes stall in the polling barrier.
      call mpi_comm_split(mpi_comm_world, color, myid, 
     &                    soppa_comm_active,ierr)
      if (myid.ge.ntot) then
         call pollingbarrier(1)
         goto  800 !skip to assembly of matrices. 
      endif


C     Start the dynamic load balancing
      if(forcedynamicloadbalance) then

          forcedynamicloadbalance = .false. ! This is the last time work is redistributed in an AOSOPPA computation.
          slavesupdatetasks       = .true.  ! Tels the slaves they have to receive work from the master.
          if (myid.eq.0) then
             call dynloadbal_parsoppa(timings, allocsize, temparray)
             AssignedIndices = temparray(1:allocsize)
          endif
      endif

C     Start the presorting load balancing.
      if (forceloadbalance) then

         forceloadbalance        = .false. ! Don't use the pre-sorting load balancing more than once, next time use dynamic instead.
         forcedynamicloadbalance = .true.  ! Monitor the tasks and re-balance once the actual walltimes are known.
         slavesupdatetasks       = .true.

         !performs balancing and sends tasks to slaves.
         if (myid.eq.0) then 
            call presortloadbal_parsoppa(work,lwork,allocsize,temparray)
            maxnumjobs = allocsize
            if (.not. allocated(AssignedIndices) ) then
               allocate(AssignedIndices(maxnumjobs), stat=allocstatus)
            endif
            AssignedIndices = temparray(1:allocsize)
         endif 

         if (.not. allocated(timings) ) then
            allocate(timings(ntot) )
            timings = 0.0
         endif

      endif !forceloadbalance.


C************************************************** 
C   Slaves receive the indices they have to work with
C**************************************************
         if (slavesupdatetasks) then
         slavesupdatetasks = .false.
         IF (myid.ne.0) then

C           Get the number of tasks and allocate memory for the indices
C            call mpi_irecv(maxnumjobs, 1, mpi_integer, 0, myid,
C     &                  mpi_comm_world, getnumjobs, ierr)
C001         continue
C            call mpi_test(getnumjobs, assignnow, mpistatus, ierr)
C            if (assignnow)  then
            maxnumjobs = ntot - min(nodtot+1,ntot) + 1
               if (.not. allocated(AssignedIndices) ) then
                  allocate(AssignedIndices(maxnumjobs),stat=allocstatus)
                  if(.not.(allocstatus.eq.0) ) then
                     call quit('Allocation error; AssignedIndices ')
                  endif
               endif
C            else 
C               goto 001
C            endif

C           Get the pre-sorted ILLL indices.
            call mpi_scatter(AssignedIndices, maxnumjobs, mpi_integer, 
     &                       AssignedIndices, maxnumjobs, mpi_integer,
     &                       0, soppa_comm_active, ierr)
C002         continue 
C            call mpi_test(getindices, workreceived, mpistatus, ierr)
C            if (workreceived)  then
C               continue 
C            else 
C               goto 002
C            endif

         ENDIF 
         endif


C *******************************************
C            Primary parallel loop. 
C *******************************************
      doilll: DO 300 workindex=1, maxnumjobs
         ILLL = AssignedIndices(workindex) 
         if (ILLL.eq.0) exit

         if (forcedynamicloadbalance) then
            call system_clock(counti, count_rate)
         endif

                  if (.not.direct) then
                     call quit('RPA in parallel must use .DIRECT. ',
     &                         'integral transformations')
                  endif

                  KEND2  = KENDSV
                  LWORK2 = LWORKSV
                  IF (HERDIR) THEN 
                    CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
     &                          IPRINT)
                  ELSE


                     CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                           WORK(KODCL1),WORK(KODCL2),
     &                           WORK(KODBC1),WORK(KODBC2),
     &                           WORK(KRDBC1),WORK(KRDBC2),
     &                           WORK(KODPP1),WORK(KODPP2),
     &                           WORK(KRDPP1),WORK(KRDPP2),
     &                           WORK(KCCFB1),WORK(KINDXB), 
     &                           WORK(KEND2),LWORK2,IPRINT)
                  ENDIF
                  LRECNR = ((NBUFX(0) - 1)/IRAT) + 1 
                  KRECNR  = KEND2
                  KEND2   = KRECNR + LRECNR
                  LWORK2  = LWORK  - KEND2

                  CALL SO_MEMMAX ('RP_ERES.2B',LWORK2)
                  IF (LWORK2 .LT. 0) 
     &               CALL STOPIT('RP_ERES.2B',' ',KEND2,LWORK)

C------------------------------------------------------------------------
C              Loop over number of distributions in disk.
C------------------------------------------------------------------------
               idel2loop: DO 310 IDEL2 = 1,NUMDIS

                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
                  IT2DEL(IDEL) = ICDEL1
                  ICDEL1       = ICDEL1 + NT2BCD(ISYDIS)

C---------------------------------------------
C                 Work space allocation no. 3.
C---------------------------------------------
                  LXINT  = NDISAO(ISYDIS)
                  KXINT   = KEND2
                  KEND3   = KXINT  + LXINT
                  LWORK3  = LWORK  - KEND3

                  CALL SO_MEMMAX ('RP_ERES.3',LWORK3)
                  IF (LWORK3 .LT. 0) 
     &               CALL STOPIT('RP_ERES.3',' ',KEND3,LWORK)

C--------------------------------------------
C                 Read in batch of integrals.
C--------------------------------------------
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                        WORK(KRECNR),DIRECT)
C---------------------------------------------
C                 Work space allocation no. 4.
C---------------------------------------------
                  LDSRHF = NDSRHF(ISYMD)
                  KDSRHF  = KEND3
                  KEND4   = KDSRHF + LDSRHF
                  LWORK4  = LWORK  - KEND4

                  CALL SO_MEMMAX ('RP_ERES.4',LWORK4)
                  IF (LWORK4 .LT. 0) 
     &               CALL STOPIT('RP_ERES.4',' ',KEND4,LWORK)
C----------------------------------------------
C                 Calculate the AO-Fock matrix.
C----------------------------------------------
                  IF (TRIPLET) THEN
                     CALL CC_AOFOCK3(WORK(KXINT),WORK(KDENS),
     &                               WORK(KFOCK),WORK(KEND4),
     &                               LWORK4,IDEL,ISYMD,ISYMTR)
                  ELSE
                     CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),
     &                              WORK(KFOCK),WORK(KEND4),
     &                              LWORK4,IDEL,ISYMD,.FALSE.,
     &                              'Dummy',ISYMTR)

                  END IF
  310    enddo idel2loop

         if (forcedynamicloadbalance) then
            call system_clock(countf)
            timespent =   ( countf-counti )/REAL(count_rate)
            timings(ILLL) = timespent
         endif


  300 enddo doilll


      if (forcedynamicloadbalance ) then
              !assemble the array with timings 
              if (myid.eq.0) then
                 call mpi_reduce(mpi_in_place, timings, ntot, MPI_REAL,
     &                           MPI_SUM, 0, soppa_comm_active, ierr)
              else
                 call mpi_reduce(timings, temparray, ntot, MPI_REAL,
     &                           MPI_SUM, 0, soppa_comm_active, ierr)
              endif
      endif


C**************************************************
C    End of parallel calculation
C**************************************************


C     The master sends a release signal to all 
C     processes participating in the spinning barrier.
      if (numprocs.gt.ntot) then
         if (myid.eq.0) then
            do i=ntot, (numprocs-1)
               call mpi_isend(.true., 1, MPI_LOGICAL, i,
     &         i, mpi_comm_world, request, ierr)
            enddo
         endif
      endif
800   continue  ! redirection point for slaves from the spinning barrier 


C     Assemble the Fock Matrix in the master's work array

      allocate(tempfock(lfock))
      tempfock = 0.00000

      if  (myid.eq.0) then
         call mpi_reduce(MPI_IN_PLACE, work(kfock), lfock, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      else
         call mpi_reduce(work(kfock), tempfock, lfock, MPI_REAL8,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      endif

      call mpi_barrier(mpi_comm_world, ierr)


      deallocate(tempfock)

      !if (allocated(AssignedIndices)) then
      !deallocate( AssignedIndices, stat=deallocstatus)
      !if(.not.(allocstatus.eq.0) ) then
         !call quit('Deallocation error in PAR_RP_ERES')
      !endif
      !endif



      if (allocated(presortarray)) then
      deallocate(presortarray, stat=deallocstatus)
      if(.not.(deallocstatus.eq.0) ) then
         call quit('Presortarray was not deallocated in PAR_RP_ERES')
      endif
      endif


      if (allocated(finalsorted)) then
      deallocate(finalsorted, stat=deallocstatus)
      if(.not.(deallocstatus.eq.0) ) then
         call quit('Finalsorted was not deallocated in PAR_RP_ERES')
      endif
      endif
         
      RETURN
      END SUBROUTINE

#endif
!!! end of VAR_MPI 





#ifdef VAR_MPI
C /* deck getbytespan */
      subroutine getbytespan(firstvar, lastvar, bytespan)
C Frederik Beyer, March 2014.
C
C This subroutine calculates the memory span in bytes between two variables.
C
C This is used for easy updating of common block in parallel calculations.
C The former approaches relied on counting the number of occurrences
C of variables of a particular type and then transferring a common in
C a series of mpi_bcasts; one for every datatype in the common block.
C This approach utilizes the contiguous storage of variables in a common
C block to calculate the span in bytes from the first variable
C in the common block and up to, but not including, the last variable.
C
C Point the firstvar to the first variable in the common block  
C and point lastvar to the last variable in the common block.
C
C Lastvar should have a name like <commonblockname>LAST
C ie. CCSDGNINPLAST (see include/ccsdinp.h for an example)
C
C These <name>LAST variables are only there to facilitate easy common block
C transfers. They are never explicitly needed in a calculation. 
C They should always be of type int.
C
C Getbytespan calculates the total amount of bytes needed for an MPI_BCAST 
C to transfer the whole block in one go (with datatype = mpi_byte), 
C including the first but excluding the <name>last variable.
C
C Example of use:
C to update the common block /eribuf/:
C
C      call getbytespan(lbuf, eribufLAST, bytesize)
C      call mpi_bcast(lbuf, bytesize, mpi_byte, 0, mpi_comm_world, ierr)

        
         !use mpi
#include "implicit.h"

#include "priunit.h"
#include "ccsdinp.h"
#include "mpif.h"

      integer, intent(out) :: bytespan
      intent(in)           :: firstvar, lastvar
CRF This should be mpi_address_kind
      integer(mpi_address_kind) :: firstmem, lastmem
      integer              :: totaltransfer = 0
      integer, parameter   :: approxeager = 4096 ! Approximate upper limit for eager protocol transfers. Implementation dependent.
      integer              :: numeagersends=0, numrendezsends=0 



      call mpi_get_address(firstvar, firstmem, ierr) 
      call mpi_get_address(lastvar, lastmem, ierr) 
      bytespan = lastmem - firstmem 
      if (bytespan.lt.1) then
          call quit('subroutine getbytespan calculated a non-sensical ',
     &    'size for common block transfer.')
      endif


      if (debug) then
         totaltransfer = bytespan + totaltransfer
         if (bytespan.lt.approxeager) then
            numeagersends   = numeagersends + 1
         else
            numrendezsends  = numrendezsends + 1
         endif

         write(lupri, '(a, i8)') 
     &   "Running amount of Bytes transferred via getbytespan: "
     &   ,totaltransfer
         write(lupri, '(a, i8)') "Estimated transfers using the ",
     &               "eager protocol: ", numeagersends
         write(lupri, '(a, i8)') "Estimated transfers using the ",
     &                "rendezvous protocol: ", numrendezsends
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



C /* deck getbigjob*/
      subroutine getbigjob(list, cols, nextindex)
C
C NOTE !!!! As of April 2014 this routine is not used by par_rp_eres as indicated. 
C A pre-sorting algorithm has been implemented in its stead. 
C
C
C
C
C F.Beyer Mar. 2014.
C Selection routine for dynamic load balancing for the par_rp_eres subroutine.
C Before doing two-electron integrals, the amount of work associated with every 
C AO index (ILLL) is estimated and a list of work and corresponding indexes is created 
C by the subroutine presortaodist.
C
C Fetch the largest available job from the supplied list.
C Once the job is found, set the matrix entry that describes the amount of work to 0, 
C update the nextindex variable to the corresponding ILLL index and return.
C
C First row in the list is the number of distributions associated with the integral.
C Second row in the list is the ILLL index that corresponds to the amount of work.

#include "implicit.h"
      integer, intent(in) :: cols
      integer, dimension(2,cols), intent(inout) :: list
      integer, intent(out)  :: nextindex
      integer, dimension(2) :: temparray
      integer :: highdex
     
      temparray = maxloc(list, dim=2, mask=list.gt.0)
      highdex = temparray(1)
C TODO Insert error checking in case no value is found!
      nextindex = list(2, highdex) 
      list(1, highdex ) = 0 

      return
     
      end subroutine


C /* deck getsmalljob*/
      subroutine getsmalljob(list, cols, nextindex)
C
C NOTE !!!! As of April 2014 this routine is not used by par_rp_eres as indicated. 
C A pre-sorting algorithm has been implemented in its stead. 
C
C
C
C F.Beyer Mar. 2014.
C Selection routine for dynamic load balancing for the par_rp_eres subroutine.
C Before doing two-electron integrals, the amount of work associated with every 
C distribution index (ILLL) is estimated and a list of work and corresponding indexes is created 
C by the subroutine presortaodist.
C
C Fetch the smallest available job from the supplied list.
C Once the job is found, set the matrix entry that describes the amount of work to 0, 
C update the nextindex variable to the corresponding ILLL index and return.
C
C First row in the list is the number of AOs associated with the distribution.
C Second row in the list is the ILLL index that corresponds to the amount of work.

#include "implicit.h"
      integer, intent(in) :: cols
      integer, dimension(2,cols), intent(inout) :: list
      integer, intent(out)  :: nextindex
      integer, dimension(2) :: temparray
      integer :: lowdex
      
      temparray = minloc(list, dim=2, mask=list.gt.0)
      lowdex = temparray(1)
C TODO Insert error checking in case no value is found!
      nextindex = list(2, lowdex) 
      list(1, lowdex ) = 0 

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
     &                               sorted)
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

      integer, dimension(:,:), allocatable :: sumofwork
      integer, dimension(2) :: tempavail, templazy
      integer :: numprocs, availloc, targetrow
      integer :: allocstatus, deallocstatus

      call mpi_comm_size(mpi_comm_world, numprocs, ierr)
      allocate (sumofwork(2, numprocs) , stat=allocstatus)
      if(.not.(allocstatus.eq.0) ) then
         call quit('Allocation error in partitionAOindices')
      endif
 
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
      
      deallocate(sumofwork, stat=deallocstatus)
      if(.not.(deallocstatus.eq.0) ) then
         call quit('Deallocation error in partitionAOindices')
      endif

      return

      END SUBROUTINE



#endif
!!! End of VAR_MPI preprocessing flag from subroutine getbytespan
