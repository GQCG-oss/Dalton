#define IPRPAR_DEBUG -1
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
C
C FILE: dalpar.F
C
C (until Jan. 2017 these subroutines were in eri/eri2par.F)
C
C  /* Deck dalton_nodedriver */
#if defined (VAR_MPI)
      SUBROUTINE DALTON_NODEDRIVER(MWORK,WRKDLM)
C
C Slaves are sent to this nodedriver from abacus/dalton.F where they
C idle until the keyword IPRTYP is set with an MPI_BCAST (this is
C mediated via the subroutine MPIXBCAST from /gp/mpimacro.F)
C Once the slaves have received the IPRTYP.ne.0 broadcast, they need a
C broadcast with the print level.
C IPRTYP is set in iprtyp.h. This file associates every parallel
C computation in Dalton with an specific integer and a character string.
C This iprtyp.h file has to be imported in every subroutine that wants to
C communicate with this nodedriver.
C
C If you add a parallel routine to Dalton it needs to be added to
C iprtyp.h
C
C  This idling works because MPIXBCAST in a collective and blocking
C communication. Every slave will be blocked by the broadcast at label
C 100 until the master participates in the broadcast. As soon as an
C IPRTYP broadcast is made, every slave executes the computation
C associated with the broadcast once, exits the if-statement and is
C blocked/idles once again by the MPI_COMM_WORLD broadcast at label 100.
C
C If IPRTYP is set to 0 the nodes exit the loop.
C
C As of Aug 2013, you have to broadcast the IPRTYP from the master node
C using MPIXBCAST or MPI_BCAST, to communicate with the nodedriver.
C This mean you can change IPRTYP locally if you need it for whatever
C reason. To start a parallel calculation the master HAS to broadcast
C the IPRTYP and the needed parameters for the calculation.
C
C IMPORTANT: The communication is always via the communicator
C MPI_COMM_WORLD, so if another communicator is made, only those nodes
C included in the new communicator will receive the broadcast IPRTYP.
C This is because the nodedriver uses the MPIXBCAST subroutine to
C communicate and as of Nov 2014 this is hardcoded to MPI_COMM_WORLD.
C
#ifdef VAR_IFORT
!     subroutine sleep is in IFPORT module
      use IFPORT, ONLY: SLEEPQQ
#endif
#ifdef MOD_CGTO_DIFF_ERI
      use eri_basis_loops
#endif
#ifdef ENABLE_XCINT
      use dalton_xcint_interface, only: xcint_integrate_worker
#endif
      use so_parutils, only: soppa_nodedriver
      use pelib_interface, only: pelib_ifc_slave
#if defined (HAS_PCMSOLVER)
      use pcm_parallel, only: j1ints_pcm
#endif

#include "implicit.h"
#include "priunit.h"
#include "mtags.h"
C infpar.h: MASTER,MYNUM, PARHER
#include "maxorb.h"
#include "infpar.h"
C iprtyp.h: define parallel calculation types
#include "iprtyp.h"
#include "mpif.h"

C gnrinf.h : PARCAL, DIRCAL
#include "gnrinf.h"

      real(8), allocatable  :: work(:)
      integer(4), parameter :: TEN_4 = 10

! ======================================================

      CALL QENTER('DALTON_NODEDRIVER')
      CALL GETTIM(CSTR,WSTR)

      PARCAL = .TRUE.
      DIRCAL = .TRUE.
      PARHER = .TRUE.

      N_TASKS = 0
C
C     Get calculation type (IPRTYP)
  100 CONTINUE
#ifdef VAR_IFORT
!     Nov 2013 hjaaj: sleep 10 millisecond before next check for a task
!     - this saves a lot of cpu time on other cores in a sequential section
!     - may make it possible to go into turbo mode
!     - may make it possible to exploit e.g. parallel DGEMM in MKL in serial sections
C      CALL SLEEP(1)
c
!     - only intel ifort has a sleepqq counting in milliseconds (that I know of),
!       the more common sleep(n) sleeps in at least one second, which is too long.
!       (One could use an interface to the C routine usleep.)
      CALL SLEEPQQ(TEN_4)
#endif
#ifdef VAR_GFORTRAN
!      CALL SLEEP(1)
#endif
#ifdef COPIED_FROM_NET_CHECK_IT_OUT
! milliseconds wait
use iso_c_binding
implicit none
interface
   subroutine usleep(useconds) bind(C)
!  integer(c_int) function usleep(useconds) bind(C)
    use iso_c_binding
    implicit none
    integer(c_int32_t), value :: useconds
!  end function
   end subroutine
end interface
  print *, 'before'
  call usleep(500000_c_int32_t) ! sleep for 1/2 s
  print *, 'after'
end
#endif
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
C
      IF (IPRTYP.GT.0) THEN
         N_TASKS = N_TASKS + 1

         CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)
         IF (IPRINT .GT. 1)
     &      write (lupri,'(/A,I5,A,I8,A,I5,A,I5/A)')
     &      'Node',MYNUM,' received task',N_TASKS,' IPRTYP =',IPRTYP,
     &      ' Print level is',IPRINT,
     &      '(See DALTON/include/iprtyp.h for IPRTYP codes)'
            flush(lupri)

!        allocate memory here if we do not call hermit (more to follow...).
!        idea: move allocation down to main driver routines of modules as we might want to use
!              dynamic memory allocation down there (for MPI performance/MPI standard compliance reasons)
!              which is severely hampered by a "static" WORK for everything.
!              sknecht - feb 2013.
         if(iprtyp /= HER_WORK .and. .not.(allocated(work)))then
!          allocate memory using f90 utilities
           allocate(work(0:mwork+1),stat=i)
!          Set memory traps
           work(0)       = wrkdlm
           work(1+mwork) = wrkdlm

           if(i /= 0)then
             write (lupri,*) mynum,
     &        ': ALLOCATE command failed to allocate'//
     &        ' the requested memory ',mwork,' on slave. Error code:',i
             write (lupri,*)
     &        'Reduce the memory demands and be welcome back'
             call quit('Failed to allocate memory on slave')
           end if
         end if

C  new calculation type needs to be added in include/iprtyp.h
         IF      (IPRTYP.EQ.HER_WORK) THEN
            CALL HER_NODSTR(MWORK,wrkdlm,IPRINT)
         ELSE IF (IPRTYP.EQ.ERI_WORK) THEN
            CALL ERI_NODSTR(WORK(1),MWORK,IPRINT)
C        ELSE IF (IPRTYP.EQ.GEM_WORK) THEN
C           CALL GEMNOD(WORK(1),MWORK,MASTER,MYNUM,MYNUM,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_KSM_WORK) THEN
            CALL DFT_KSMSLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_C_WORK) THEN
            CALL DFT_CSLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_GRAD_WORK) THEN
            CALL DFT_GRAD_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_EXPGRAD_WORK) THEN
            CALL DFT_EXPG_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_BRHS_WORK) THEN
            CALL DFT_BRHS_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_HESSTAT_WORK) THEN
            CALL DFT_HESSTAT_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFTHED_WORK) THEN
            CALL DFTHED_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_SUSCEP_WORK) THEN
             CALL DFT_SUSCEP_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.NPETES2_WORK .OR.
     &            IPRTYP.EQ.PCMBSOL_WORK) THEN
            CALL J1INTS(WORK(1),MWORK,IPRTYP,IPRINT)
         ELSE IF (IPRTYP.EQ.J1XP_WORK) THEN
            CALL J1XS(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.J2XP_WORK) THEN
            CALL J2XS(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_MOMO_WORK) THEN
            CALL DFT_MOMO_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3_NSP_WORK) THEN
            CALL QM3_NSPS(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3LNO_WORK) THEN
            CALL QM3LNO_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3QRO_WORK) THEN
            CALL QM3QRO_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3FIRST_1_WORK) THEN
            CALL QM3FIRST_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3FIRST_2_WORK) THEN
            CALL QM3FIRST_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.PARQMMM__WORK) THEN
            CALL PARQMMM_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MM_FIELD_1_WORK) THEN
            CALL MM_FIELD_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MM_FIELD_2_WORK) THEN
            CALL MM_FIELD_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MM_POLAR_CONTR_WORK) THEN
            CALL MM_POLAR_CONTR_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MMITER_INDDIP_WORK) THEN
            CALL MMITER_INDDIP_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMLNO_1_WORK) THEN
            CALL QMMMLNO_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMLNO_2_WORK) THEN
            CALL QMMMLNO_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMQRO_1_WORK) THEN
            CALL QMMMQRO_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMQRO_2_WORK) THEN
            CALL QMMMQRO_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMFIRST_1_WORK) THEN
            CALL QMMMFIRST_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMFIRST_2_WORK) THEN
            CALL QMMMFIRST_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMFIRST_3_WORK) THEN
            CALL QMMMFIRST_S3(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMFIRST_4_WORK) THEN
            CALL QMMMFIRST_S4(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_1_WORK) THEN
            CALL QMMMB2_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_2_WORK) THEN
            CALL QMMMB2_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_3_WORK) THEN
            CALL QMMMB2_S3(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_4_WORK) THEN
            CALL QMMMB2_S4(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.LUCITA_WORK) THEN
            call lucita_coworker_main(work(1),mwork)
         ELSE IF (IPRTYP.EQ.MEP_WORK) THEN
            CALL MEP_SLAVE(WORK(1),MWORK)
         ELSE IF (IPRTYP == POLARIZABLE_EMBEDDING) THEN
            CALL PELIB_IFC_SLAVE(IPRINT)
#if defined(BUILD_GEN1INT)
C added by Bin Gao, May 12, 2012
         ELSE IF (IPRTYP.EQ.GEN1INT_INIT) THEN
            CALL gen1int_worker_init()
         ELSE IF (IPRTYP.EQ.GEN1INT_GET_INT) THEN
            CALL gen1int_worker_get_int(MWORK,WORK(1),LUPRI,IPRINT)
         ELSE IF (IPRTYP.EQ.GEN1INT_GET_EXPT) THEN
            CALL gen1int_worker_get_expt(MWORK,WORK(1),LUPRI,IPRINT)
         ELSE IF (IPRTYP.EQ.GEN1INT_GET_CUBE) THEN
            CALL gen1int_worker_get_cube(MWORK,WORK(1),LUPRI,IPRINT)
#endif
         ELSE IF (IPRTYP.EQ.PARA_SO_ERES) THEN
            CALL SOPPA_NODEDRIVER(work(1), mwork, iprint )
#ifdef MOD_CGTO_DIFF_ERI
         ! added by Bin Gao, Oct. 27, 2012
         ! cgto-diff-eri code to calculate arbitray order geometric derivatives
         ! of two-electron integrals
         ELSE IF (IPRTYP.EQ.CGTO_DIFF_ERI_INIT) THEN
            CALL unopt_geodiff_loop_worker(LUPRI,IPRINT)
#endif

#ifdef ENABLE_XCINT
         else if (iprtyp == XCINT_MPI_WAKEUP_SIGNAL) then
            call xcint_integrate_worker(MPI_COMM_WORLD)
#endif
         ELSE IF (IPRTYP.EQ.QMCMM_WORK) THEN
#ifdef ENABLE_VPOTDAMP
            CALL vpotdamped_slave()
#else
            call quit('VPOTDAMP not compiled in this version')
#endif

#if defined(ENABLE_EFS)
         else if (IPRTYP.eq.efs_add_atomtype_work) then
            call efs_add_atomtype_slave()
         else if (IPRTYP.eq.efs_init_work) then
            call efs_init_slave()
         else if (IPRTYP.eq.efs_generate_basis_work) then
            call efs_generate_basis_product_slave()
         else if (IPRTYP.eq.efs_init_2efock_work) then
            call efs_init_2efock_slave()
         else if (IPRTYP.eq.efs_fock_update_work) then
            call efs_fock_update_slave()
#endif
#if defined (HAS_PCMSOLVER)
         else if (iprtyp.eq.pcmsolver_npetes2_work) then
            call j1ints_pcm(work(1),mwork,iprtyp,iprint)
         else if (iprtyp.eq.pcmsolver_pcmbsol_work) then
            call j1ints_pcm(work(1),mwork,iprtyp,iprint)
#endif

         else if (IPRTYP .eq. CALL_FCKTRA_DISTRIBUTED) THEN
            call FCKTRA_DISTRIBUTED_NODE(WORK(1),MWORK,IPRINT)
         ELSE
            WRITE(LUPRI,'(/A,I5,A,I10)')
     &         '(DALTON_NODEDRIVER) ERROR on slave',mynum,
     &         ': Unknown calculation type request',IPRTYP
            WRITE(*,'(/A,I5,A,I10)')
     &         '(DALTON_NODEDRIVER) ERROR on slave',mynum,
     &         ': Unknown calculation type request',IPRTYP
            CALL QUIT('DALTON_NODEDRIVER: Unknown calculation type')
         END IF
C
         GO TO 100
      END IF ! IPRTYP .gt. 0

!     release memory
      if(allocated(work)) deallocate(work)
C
      IF (MYNUM .LE. 20) THEN
         CALL GETTIM(CEND,WEND)
         CTOT = CEND - CSTR
         WTOT = WEND - WSTR
C
         WRITE (LUPRI,'(/A,I10/)') 'Number of tasks performed:',N_TASKS
         CALL TIMTXT(' Total CPU  time used:',CTOT,LUPRI)
         CALL TIMTXT(' Total wall time used:',WTOT,LUPRI)
C
C        Stamp date and time and hostname to output
C
         CALL TSTAMP(' ',LUPRI)
      END IF
C
      CALL QEXIT('DALTON_NODEDRIVER')
      RETURN
      END
#endif
C  /* Deck sendname */
      SUBROUTINE SENDNAME(NASTER,NYTID)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
C
#include "infpar.h"
C
      MASTER = NASTER
      MYTID  = NYTID
C
C     Get my name
C
      CALL GETNAME(MYNAME)
C
C     Return name to master
C
#if defined (VAR_MPI)
      CALL MPIXSEND(MYNAME,20,'STRING',MASTER,MTAG8)
#endif
#if defined (VAR_PVM)
C
C     NCODE is the first variable PVM-slaves must receive
C
      CALL PVMXRECV(MASTER,MTAG0)
      CALL PVMXUNPACK(NCODE,1,'INTEGER')
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(MYTID,1,'INTEGER')
      CALL PVMXPACK(MYNAME,20,'STRING')
      CALL PVMXSEND(NASTER,MTAG8)
#endif
C
      RETURN
      END
C  /* Deck recvnames */
      SUBROUTINE RECVNAMES
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
      CHARACTER NAME*20
#include "infpar.h"
C
#if defined (VAR_PVM)
C
C     NCODE must be the first variable PVM-slaves must receive
C
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(NCODE,1,'INTEGER')
      CALL PVMXMCAST(NODTOT,NODEID(1),MTAG0)
#endif
C     Get my name
C
      CALL GETNAME(MYNAME)
      NODNAM(0) = MYNAME
C
C     Return name to master
C
      DO I = 1,NODTOT
#if defined (VAR_MPI)
         NWHO = -1
         CALL MPIXRECV(NAME,20 ,'STRING',NWHO,MTAG8)
         NODNAM(NWHO) = NAME
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG8)
         CALL PVMXUNPACK(NTID,1,'INTEGER')
         CALL PVMXUNPACK(NAME,20,'STRING')
         DO J=1, NODTOT
            IF (NODEID(J) .EQ. NTID) NWHO = J
         END DO
         NODNAM(NWHO) = NAME
#endif
      END DO
C
      IF (IPRPAR.GT.0) THEN
         WRITE (LUPRI,'(3X,A,A)') ' Slaves are : ',NODNAM(1)
         WRITE (LUPRI,'(17X,A)') (NODNAM(I),I=2,NODTOT)
      END IF
C
      RETURN
      END
C  /* Deck sdtim */
      SUBROUTINE SDTIM(CPU1,WALL1,IPRINT)
C
C Used from common blocks:
C
C MXORB  : MXSHEL (for INFPAR)
C INFPAR : MASTER, MYNUM, NCODE
C
#include "implicit.h"
#include "maxorb.h"
#include "mtags.h"
      REAL*8 TMGARR(2)
#include "infpar.h"
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from SDTIM','*',103)
C
      CALL GETTIM(CPU2,WALL2)
C
      TMGARR(1) = CPU2  - CPU1
      TMGARR(2) = WALL2 - WALL1
C
#if defined (VAR_MPI)
      CALL MPIXSEND(TMGARR,2 ,'DOUBLE',MASTER,MTAG8)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(MYNUM ,1,'INTEGER')
      CALL PVMXPACK(TMGARR,2,'DOUBLE')
      CALL PVMXSEND(MASTER,MTAG8)
#endif
C
      RETURN
      END
C  /* Deck rvtim */
      SUBROUTINE RVTIM(NSTAT,CPUS,WALLS,NODNAM,NODTOT,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "mtags.h"
      CHARACTER*20 NAME, NODNAM(NODTOT)
      DIMENSION NSTAT(NODTOT),CPUS(NODTOT),WALLS(NODTOT)
      REAL*8 TMGARR(2)
C
      SAVE TSOFAR
      DATA TSOFAR/0.0D0/
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from RVTIM','*',103)
C

      DO 100 I = 1,NODTOT
#if defined (VAR_MPI)
         NWHO = -1
         CALL MPIXRECV(TMGARR,2,'DOUBLE',NWHO,MTAG8)
#elif defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG8)
         CALL PVMXUNPACK(NWHO  ,1,'INTEGER')
         CALL PVMXUNPACK(TMGARR,2,'DOUBLE')
#else
         NWHO = I
#endif
         CPUS(NWHO)   = TMGARR(1)
         WALLS(NWHO)  = TMGARR(2)
         TSOFAR = TSOFAR + CPUS(NWHO)
C
  100 CONTINUE
C
      IF (IPRINT .LE. 0) GOTO 9000
      IF (IPRINT .GT. 1) THEN
         CALL AROUND('Overall statistics for distribution of batches')
C
         WRITE(LUPRI,'(9X,5(A,4X))')   'Nodename',
     &                                 'Batches',
     &                                 'CPU time',
     &                                 'Wall time',
     &                                 'Efficiency'
         WRITE(LUPRI,'(9X,5(A,4X)/)')  '--------',
     &                                 '-------',
     &                                 '--------',
     &                                 '---------',
     &                                 '----------'
C
         DO 200 NWHO = 1, NODTOT
            IF (WALLS(NWHO) .NE. 0.0D0) THEN
               EFFI = (CPUS(NWHO)/WALLS(NWHO))*100
            ELSE
               EFFI = 0.0D0
            END IF
            WRITE(LUPRI,'(9X,A10,I7,6X,F8.2,5X,F8.2,6X,F6.2)')
     &           NODNAM(NWHO),NSTAT(NWHO),CPUS(NWHO),WALLS(NWHO),EFFI
 200     CONTINUE
      END IF
C
      IMINS  = INT(NINT(TSOFAR)/60)
      IHOURS = INT(IMINS/60)
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(TSOFAR) - 3600*IHOURS - 60*IMINS
C
      WRITE(LUPRI,'(/A,I5.2,A,I2.2,A,I2.2)')
     &     ' Total CPU  time used in NODES  so far   :',
     &     IHOURS,':',IMINS,':',ISECS
      CALL FLSHFO(LUPRI)
C
 9000 RETURN
      END
C  /* Deck pardrv */
      SUBROUTINE PARDRV(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,WORK,NSTAT,
     &                  HESSEE,LWORK,ITYPE,MAXDIF,IATOM,NODV,NOPV,
     &                  NOCONT,TKTIME,RETUR,IPRINT,IPRTYP,
     &                  ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,GENCTR)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,GENCTR
      DIMENSION FMAT(*),DMAT(*),ISYMDM(*),IFCTYP(*),WORK(LWORK),
     &          NSTAT(NODTOT), HESSEE(*)
CTROND
      DIMENSION GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)
CTROND
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C INFPAR : MASTER, NODTOT, IPRPAR,  NCODE, DBG
C iprtyp.h : HER_WORK, ERI_WORK
C GNRINF : BASDIR, LENBAS
C
#include "infpar.h"
#include "iprtyp.h"
#include "gnrinf.h"
      DATA TOTWAL /0./
      SAVE TOTWAL
C
      CALL QENTER('PARDRV')
C
!     IPRPAR = 11
      IPRINT = MAX(IPRINT,IPRPAR)
C
      IF (IPRPAR .GT. 0) CALL GETTIM(CPU1,WALL1)
C
      IF (INFPAR_DEBUG) THEN
         IPRNOD = IPRINT
      ELSE
         IPRNOD = 0
      END IF
      IPRNOD = MAX(IPRNOD,IPRPAR_DEBUG)
      IF (IPRINT .GT. 3) THEN
         CALL TITLER('Output from PARDRV','*',103)
         WRITE(LUPRI,'(/A,2I5)')
     &      'Parallel types IPRTYP and ITYPE are',IPRTYP,ITYPE,
     &      'Print level on nodes is',IPRNOD
      END IF
C
C     Tell nodes whether this is an ERI or HERMIT run. Send IPRNOD (PVM: NCODE)
C
#if defined (VAR_MPI)
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRNOD,1,'INTEGER',MASTER)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(IPRTYP,1,'INTEGER')
      CALL PVMXPACK(IPRNOD,1,'INTEGER')
      CALL PVMXMCAST(NODTOT,NODEID(1),MTAG1)
#endif
      IF (IPRTYP.EQ.HER_WORK) THEN
         if (iprint .gt. 3)
     &      write (lupri,*) 'Calling parallel HERMIT routines.'
         CALL HER_PARDRV(FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                   WORK,LWORK,ITYPE,MAXDIF,IATOM,NODV,
     &                   NOPV,NOCONT,TKTIME,RETUR,IPRINT,
     &                   ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,
     &                   DINTSKP,RELCAL,GENCTR)
      ELSE IF (IPRTYP.EQ.ERI_WORK) THEN
         if (iprint .gt. 3)
     &      write (lupri,*) 'Calling parallel ERI routines.'
         CALL ERI_PARDRV(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                   WORK,LWORK,IPRINT)
      ELSE
         write (lupri,*) 'PARDRV: Invalid IPRTYP = ',IPRTYP
         CALL QUIT('PARDRV: Invalid IPRTYP here')
      END IF
C
      KCPUS  = 1
      KWALLS = KCPUS  + NODTOT
      KLAST  = KWALLS + NODTOT
      IF (KLAST .GT. LWORK)
     &   CALL STOPIT('PARDRV','RVTIM',KLAST,LWORK)
C
      CALL RVTIM(NSTAT,WORK(KCPUS),WORK(KWALLS),NODNAM(1),NODTOT,IPRINT)
C
      IF (IPRPAR .GT. 0) THEN
         CALL GETTIM(CPU2,WALL2)
         CPU    = CPU2 - CPU1
         WALL   = WALL2 - WALL1
         TOTWAL = TOTWAL + WALL
C
         IMINS  = INT(NINT(TOTWAL)/60)
         IHOURS = INT(IMINS/60)
         IMINS  = IMINS - 60*IHOURS
         ISECS  = NINT(TOTWAL) - 3600*IHOURS - 60*IMINS
C
         WRITE(LUPRI,'(A,I5.2,A,I2.2,A,I2.2)')
     &     ' Total wall time used in PARDRV so far   :',
     &      IHOURS,':',IMINS,':',ISECS
         IF (IPRPAR .GT. 2) THEN
            WRITE(LUPRI,'(A,F11.2,/A,F11.2)')
     &        ' CPU  time used in PARDRV last iteration :',
     &          CPU,
     &        ' Wall time used in PARDRV last iteration :',
     &          WALL
         END IF
      END IF
C
      CALL QEXIT('PARDRV')
      RETURN
      END
! --- end of main/dalpar.F ---
