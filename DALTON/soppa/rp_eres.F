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
C RP_ERES - Data dictionary:
C NOLDTR: Number of Old trial vectors
C NNEWTR: Numver of new trial vectors
C FOCKD: The Fock Matrix
C LFOCKD: Number of entries in the Fock Matrix
C ISYMTR:
C WORK: This is the name of the temporary work-vector (notice that the first element is 1 in the memory allocation)
C LWORK: Total needed length of the work-vector.

#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "eribuf.h"
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
#include "parrpa.h" !This initializes most parameters for ERIDI1 and connects with the parallel calculation.
c         integer :: numprocs, ierr
c      call mpi_comm_size(mpi_comm_world, numprocs, ierr) !This call to comm_size makes loop 220 crash with a write statement!!!!! What the actual fuck???
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
C TODO BEYER
C        Check if AODENS contains stuff that need to be broadcast to all slaves.
C        SO_AODENS makes use of DGEMM. Insert a call to utilize the slaves for this calculation
C TODO ... BEYER
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
C TODO BEYER
C        Check if BCKTR contains stuff that need to be broadcast to all slaves.
C TODO ... BEYER
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
C Experimental beyer todo
            call sleepslaves(300)
C end experiment beyer todo 
            DO 220 ILLL = 1,NTOT
C
C------------------------------------------------
C              If direct calculate the integrals.
C------------------------------------------------
C
C TODO delete/move/remove/rewrite these lines
#ifdef VAR_MPI
C qwer
C asdfasdf
      if (.false.) then
          !ISYMPAR = ISYMTR
          !call mpixbcast(85, 1, 'INTEGE', 0)
          !call mpixbcast(0, 1, 'INTEGE', 0)
          !CALL RPA_E_LIN_TRANSFORM(WORK, LWORK, 0, LUPRI, IPRINT)
          !goto 888

      endif
#endif
      
C ... F.Beyer
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
C
  220       CONTINUE ! End of ILLL loop.
C
C Experimental beyer todo
         call wakeslaves()
C end experiment beyer todo 
  210    CONTINUE


  888    continue
c BEYER: hele KFOCK matricen skal vaere samlet her! 
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
      END












#ifdef VAR_MPI
C  /* Deck RPA_E_LIN_TRANSFORM */
      SUBROUTINE RPA_E_LIN_TRANSFORM(WORK, LWORK, MYID)
C AS OF DEC 2013 THIS IS NOT A FUNCTIONING ROUTINE!!!!
C Parallel subroutine for RP_ERES.F
C Calculates the AO-FOCK matrix in parallel by distributing the two-electron integrals over several processes and assembling the partial Fock matrices from every process before returning.
C
C WORK is the array where master and every slave assembles the fock matrix incrementally
C fock is the final assembled matrix 
C Every slave calculates it own batch of intgrals depending on the last index in the two-electron integral.
C Every iteration is added as a direct sum for the final matrix and all Fock matrices are sent back to the 
C master for assembly after completing all iterations.
C 
C The CCRDAO routine is called by the master and the result is sent to every slave.
C TODO: if there is an excess of processes compared to the number of tasks, you need to set all elements in the fock matric for the surplus tasks to 0 so they can participate in the call mpi_reduce functio later on.
C TODO:  There are calls to reduce later on and this collective operation will stall because not all processes participate in the parallel calculation.
C TODO: implement a dynamic load balancing with a pool of jobs and maybe a mutex solution to handle job allocation
C TODO: Implement a synchronization routine that uses iterative non-blocking sends to avoid the barrier/broadcast call when synching arrays and common block variables.
C
C AS OF FEB 2013 THIS IS NOT A FUNCTIONING ROUTINE!!!!
C
C IMPORTANT! This routine is called from rp_eres. rp_eres takes the variable ISYMTR as an input parameter. To move ISYMTR to the common block /parrpa/ I have copied it to isympar (which is in the common block) right before calling PAR_RP_ERES. That is why ISYMTR is set equal to ISYMPAR at the beginning of this routine. If you do not update ISYMTR before calling PAR_RP_ERES the routine will not calculate two-electron integrals for the correct symmetries.





      use mpi 
#ifdef VAR_IFORT
!      subroutine sleepqq is in IFPORT module. Sleeping time is in milliseconds
      use IFPORT, ONLY: SLEEPQQ
#endif
c aaincludes xxincludes
C common blocks from RP_ERES
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "eribuf.h"
#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "soppinf.h"



#include "maxaqn.h" ! NOTE: maxaqn has to be included before the remaining include files because some of them are dependent on the initialization order.
#include "aobtch.h"
#include "odclss.h"
#include "ccom.h"
#include "ericom.h"
#include "eridst.h"
#include "r12int.h"
#include "hertop.h"
#include "nuclei.h"
#include "infpar.h"
C Common blocks from  PICKAO
#include "cbirea.h"
c#include "shells.h" ! this conflicts with ccorb.h
#include "erisel.h"


#include "iprtyp.h"
#include "parrpa.h"

C     Intent of Dummy arguments
      integer, intent(in) :: myid, lwork
      dimension work(lwork) !intent(inout) !! This uses DALTON's implicit type rather than hardcoding it to real(kind=8)
      logical :: reset, release, flag, exitspinningbarrier
      integer :: numprocs, sharedwork, rest, firstindex, lastindex
      integer :: isymtr, request 
      integer, dimension(MPI_STATUS_SIZE) :: mpistatus
      DIMENSION INDEXA(MXCORB)
      DIMENSION FOCKD(LFOCK)
      logical :: debugging = .false.
      isymtr = isympar !Inherit the symmetry from a common block.

C      CALL QENTER('RPA_E_LIN_TRANSFORM')
      call mpi_comm_size(mpi_comm_world, numprocs, ierr)


C     First broadcast updates every variable in the /common/ block PARRPA. 
C     These are all the pointers to the work array. The second bcast 
C     updates the contents of the work array. Right now the unneeded trial vectors
C     and results vectors are broadcast as well. 
C     Be aware that the mpi_real8 data type is due to the value of REAL in
C     implicit.h. 
      call mpi_bcast(NTOT, lparrp, mpi_integer, 0, mpi_comm_world, ierr)
      call mpi_bcast(work(1), kendsv, mpi_real8, 0,mpi_comm_world, ierr)


C **********************************************************
C     The load balancing depends on the master participating in the
C     parallel calculation! 
C
C     The number of indexes to be handled is distributed evenly over 
C     available processes. The tasks that are leftover once all processes
C     have received their share of jobs is distributed among the slaves
C     only. The master handles the minimum number of jobs.
C
C     !If the number of processes is larger than the 
C     number of tasks, the remaining processes will 
C     sleep to release CPU time for threaded MKL applications.
C     ERIDI2 makes use of DGEMM which will run in parallel if there are
C     available processes for MKL threading.
C
C     Once a scheme has been chosen, it redirects with a GOTO to label 666
C     This was done to it's easier to alter any of the balancing schemes later on.
C
C     TODO The excess processes need a fock matrix with only 0's so that
C     they can participate in the reduce call without adding anything.
C
C     TODO consider dynamic load balancing instead of static load balancing.
C **********************************************************
      sharedwork = ntot/numprocs 
      rest = mod(ntot, numprocs)
      if (sharedwork.eq.1) goto 200 !load balancing for low task/process ratio
      if (sharedwork.lt.1) goto 210 !load balancing when there are excess processes 


C     Normal static load balancing
      staticbalance: if (myid.eq.0) then
         firstindex = 1
         lastindex = sharedwork
      elseif ((myid.le.rest).and.(myid.ne.0)) then 
         firstindex = myid*sharedwork + myid
         lastindex = firstindex + sharedwork  
      elseif ((myid.gt.rest).and.(myid.lt.ntot)) then
         firstindex = myid*sharedwork + 1 + rest
         lastindex = firstindex + sharedwork - 1
      endif staticbalance
      goto 666 !end of normal balance


200   continue
C     In case there is only one task to be solved per process (very special case), this loadbalance is used. 
      lowbalance: if (myid.eq.0) then
         firstindex = 1
         lastindex = firstindex 
      elseif (myid.le.rest) then
         firstindex = myid*sharedwork + myid
         lastindex = firstindex + 1 
      elseif (myid.gt.rest) then
         firstindex = myid*sharedwork + rest + 1
         lastindex =  firstindex
      endif lowbalance
      goto 666 !end of sparse balance


210   continue
C     If there are more processes than tasks, this load balance is used.
      excessbalance: if (myid.lt.rest) then      
         firstindex = myid + 1 
         lastindex = myid + 1
      endif excessbalance








C   TODO is it a better idea to create a communicator for the spinning slaves that includes the master?
C   is it faster to do it that way rather than with a irecv?

C   This spinning barrier will only work if INTEL is used.
C   The sleepqq command is not available otherwise

C   Spinning barrier in case there are too many processes compared
C   to number of tasks.  The processes will sleep, enabling the 
C   remaining processes to use them for threaded MKL applications.
      exitspinningbarrier = .false. 
      spinningbarrier: if (myid.ge.ntot) then
         if (debugging) then
         endif
C        TODO set all Fock matrix elements to zero here before continuing
         call mpi_irecv(exitspinningbarrier, 1, MPI_LOGICAL, 0, myid,
     &   mpi_comm_world, request, ierr) 
230      continue 
#ifdef VAR_IFORT
         call sleepqq(100)  !sleep 100 milliseconds to release CPU time for the MKL library for working slaves.
#endif
         call mpi_test(request, flag, mpistatus, ierr)
         if (.not.exitspinningbarrier) then
            goto 230
         else
             goto 800 !Skip to the assembly of the Fock matrix 
         endif

      endif spinningbarrier



666   continue !load balancing schemes redirect here.




      if (debugging) then
         write(*, '(a, i3, a, i2, a, i2, a, i3)' )
     &   "Slave with ID:", myid, " has first index: "
     &   ,firstindex, " and last index: ", lastindex," out of: ",ntot,
     &   " total jobs."
         write(lupri, '(a, i3, a, i2, a, i2, a, i3)' )
     &   "Slave with ID:", myid, " has first index: "
     &   ,firstindex, " and last index: ", lastindex," out of: ",ntot,
     &   " total jobs."
      endif



C *******************************************
C            Primary parallel loop. 
C *******************************************
            doilll: DO 300 ILLL = firstindex, lastindex
                  write(*,'(a, i2, a, i2)')
     &     "MYID: ", myid, " Entered ILLL: ", ILLL
                  if (.not.direct) then
                     call quit('RPA in parallel must use .DIRECT.',
     & ' integral transformations')
                  endif

                  KEND2  = KENDSV
                  LWORK2 = LWORKSV

                  IF (HERDIR) THEN 
                    CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
     &                          IPRINT)
                  ELSE
C
C TODO ERIDI2 dies when it enters pickao
                     write(*,*) "about to call eridi2"
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
C
                  write(*,*) "done calling eridi2"
                  LRECNR = ((NBUFX(0) - 1)/IRAT) + 1 
C
                  KRECNR  = KEND2
                  KEND2   = KRECNR + LRECNR
                  LWORK2  = LWORK  - KEND2
C
                  CALL SO_MEMMAX ('RP_ERES.2B',LWORK2)
                  IF (LWORK2 .LT. 0) 
     &               CALL STOPIT('RP_ERES.2B',' ',KEND2,LWORK)
C
C
C------------------------------------------------------------------------
C              Loop over number of distributions in disk.
C------------------------------------------------------------------------
C
               idel2loop: DO 310 IDEL2 = 1,NUMDIS
C
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
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
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                        WORK(KRECNR),DIRECT)
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
                  IF (TRIPLET) THEN
                     CALL CC_AOFOCK3(WORK(KXINT),WORK(KDENS),
     &                               WORK(KFOCK),WORK(KEND4),
     &                               LWORK4,IDEL,ISYMD,ISYMTR)
                  ELSE

                     CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),
     &                              WORK(KFOCK),WORK(KEND4),
     &                              LWORK4,IDEL,ISYMD,.FALSE.,
     &                              'CrashIfNecessary',ISYMTR)
C Jan 30. 2014. F.Beyer:
C Ive altered the penultimate variable to 'CrashifNecessary' to make sure CC_aofock crashes
C If the boolean is ever true for any reason.
C If you use .TRUE. instead of .FALSE. in the call to CC_AOFOCK then make sure you 
C correct the variable that follows. The last variable that was in place before I changed it to 'CrashifNecessary' was uninitialized
C and could potentially cause a lot of havoc if the calculation didn't fail every time.
                  END IF
C
  310          enddo idel2loop
  300       enddo doilll

C


C     The master sends a termination signal to all surplus 
C     processes participating in the spinning barrier.
C     all IDs .le. ntot
      if (numprocs.ge.ntot) then
         if (myid.eq.0) then
            exitspinningbarrier = .true.
            do i=ntot, (numprocs-1)
               call mpi_isend(exitspinningbarrier, 1, MPI_LOGICAL, i,
     &         i, mpi_comm_world, request, ierr)
            enddo
         endif
      endif

800   continue !Slaves exiting from the spinning barrier are redirected here


C TODO make sure that the processes participating in the spinning barrier will contribute with 0 to this reduce call
C     !!! IMPORTANT !!!
C     This reduce call will only reduce for an intracommunicator call.
C     If you have split processes into different communicators, you will
C     need to reduce by other means than MPI_IN_PLACE in the master's call.
c      goto 900
      if  (myid.eq.0) then
         call mpi_reduce(MPI_IN_PLACE, fockd(1), lfock 
     & , MPI_REAL, MPI_SUM,0, MPI_COMM_WORLD, ierr)
      else
         call mpi_reduce(work(kfock), fockd(1), lfock
     &      , MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      endif
900   continue

c      write(*,*) "from parallel kfock:"
c      do i=0, 20
c         write(*,*) work(kfock+i)
c      enddo
c      write(*,*) "end parallel kfock:"

990   continue

c      deallocate(sync) 
      RETURN
      END SUBROUTINE

#endif
